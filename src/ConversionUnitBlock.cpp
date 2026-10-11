/*--------------------------------------------------------------------------*/
/*---------------------- File ConversionUnitBlock.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the ConversionUnitBlock class, which derives from
 * UnitBlock [see UnitBlock.h] and implements a unit with one commitment
 * and several coupled generators, each possibly on its own node, whose
 * operating region is given in H-form by the data [see
 * ConversionUnitBlock.h], and of its Modification and Solution.
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>
#include <cmath>
#include <map>
#include <set>

#include "ConversionUnitBlock.h"

#include "LinearFunction.h"

#include "DQuadFunction.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

using Index = Block::Index;
using Subset = Block::Subset;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register ConversionUnitBlock to the Block factory
SMSpp_insert_in_factory_cpp_1( ConversionUnitBlock );

// register ConversionUnitBlockSolution to the Solution factory
SMSpp_insert_in_factory_cpp_0( ConversionUnitBlockSolution );

/*--------------------------------------------------------------------------*/
/*-------------------------------- FUNCTIONS -------------------------------*/
/*--------------------------------------------------------------------------*/

static const double INF = Inf< double >();

/// true if the matrix has no element (the datum is absent)

static bool absent( const boost::multi_array< double , 2 > & a )
{
 return( a.num_elements() == 0 );
 }

/*--------------------------------------------------------------------------*/
/// copies b into a, whatever the shapes

static void copy_into( boost::multi_array< double , 2 > & a ,
                       const boost::multi_array< double , 2 > & b )
{
 a.resize( boost::extents[ b.shape()[ 0 ] ][ b.shape()[ 1 ] ] );
 a = b;
 }

/*--------------------------------------------------------------------------*/
/* Sorts the indices of a Subset given unordered, together with the values
 * that go with them, as in ThermalUnitBlock: the values are copied, in the
 * new order, in sorted, and values is set to its beginning. */

static void sort_by_index( Subset & subset ,
                           std::vector< double >::const_iterator & values ,
                           std::vector< double > & sorted )
{
 std::vector< std::pair< Index , double > > iv( subset.size() );
 for( Index i = 0 ; i < subset.size() ; ++i )
  iv[ i ] = std::make_pair( subset[ i ] , *( values + i ) );
 std::stable_sort( iv.begin() , iv.end() ,
                   []( const auto & a , const auto & b ) {
                    return( a.first < b.first ); } );
 sorted.resize( iv.size() );
 for( Index i = 0 ; i < iv.size() ; ++i ) {
  subset[ i ] = iv[ i ].first;
  sorted[ i ] = iv[ i ].second;
  }
 values = sorted.cbegin();
 }

/*--------------------------------------------------------------------------*/
/* The minimum up (down) time m of a unit with initial up/down time init over
 * a horizon of T instants, at least 1 and at most T + max( 1 , k ), k the
 * number of instants the unit has been on (off) before the horizon, as in
 * ThermalUnitBlock. */

static Index clamp_min_time( Index m , Index T , int init , bool up )
{
 const Index k = up ? ( init > 0 ? Index( init ) : 0 )
                    : ( init <= 0 ? Index( - init ) : 0 );
 return( std::min( std::max( m , Index( 1 ) ) ,
                   T + std::max( k , Index( 1 ) ) ) );
 }

/*--------------------------------------------------------------------------*/
// the tolerance and the kind of violation that a Configuration carries

static void extract_tolerance( Configuration * fsbc , BlockConfig * bcfg ,
                               double & tol , bool & rel_viol )
{
 auto extract = [ & tol , & rel_viol ]( Configuration * c ) -> bool {
  if( auto tc = dynamic_cast< SimpleConfiguration< double > * >( c ) ) {
   tol = tc->f_value;
   return( true );
   }
  if( auto tc = dynamic_cast< SimpleConfiguration<
                                std::pair< double , int > > * >( c ) ) {
   tol = tc->f_value.first;
   rel_viol = tc->f_value.second;
   return( true );
   }
  return( false );
  };

 if( ( ! extract( fsbc ) ) && bcfg )
  extract( bcfg->f_is_feasible_Configuration );
 }

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS OF ConversionUnitBlock ----------------------*/
/*--------------------------------------------------------------------------*/

ConversionUnitBlock::~ConversionUnitBlock()
{
 Constraint::clear( Logical_Const );
 Constraint::clear( MinUp_Const );
 Constraint::clear( MinDown_Const );
 Constraint::clear( MinPower_Const );
 Constraint::clear( MaxPower_Const );
 Constraint::clear( Operating_Const );
 Constraint::clear( RampUp_Const );
 Constraint::clear( RampDown_Const );
 Constraint::clear( PrimaryRho_Const );
 Constraint::clear( SecondaryRho_Const );
 Constraint::clear( ReserveBound_Const );
 Constraint::clear( ReserveOperating_Const );
 Constraint::clear( ReserveRamp_Const );
 Constraint::clear( Commitment_Bound );
 ShutDownZero_Bound.clear();
 Constraint::clear( ActivePower_Bound );
 Constraint::clear( StartUp_Binary_Bound );
 Constraint::clear( ShutDown_Binary_Bound );

 objective.clear();
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::deserialize( const netCDF::NcGroup & group )
{
 static const std::string who = "ConversionUnitBlock::deserialize";

 UnitBlock::deserialize( group );
 const Index T = f_time_horizon;

 // the dimensions - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! deserialize_dim( group , "NumberGenerators" , f_G ) )
  f_G = 1;
 if( f_G == 0 )
  throw( std::invalid_argument( who + ": NumberGenerators is 0" ) );

 if( ! deserialize_dim( group , "NumberOperatingRows" , f_M ) )
  f_M = 0;

 f_commitment_generator = 0;
 ::deserialize( group , f_commitment_generator , "CommitmentGenerator" );
 if( f_commitment_generator >= f_G )
  throw( std::invalid_argument( who + ": CommitmentGenerator is " +
                                std::to_string( f_commitment_generator ) +
                                " with " + std::to_string( f_G ) +
                                " generators" ) );

 // a datum indexed over the instants and over cols columns
 auto read = [ & ]( const std::string & name , Index cols , MAdbl & data ,
                    bool optional ) {
  return( ::deserialize( group , name , { T , cols } , data , optional ,
                         true , v_change_intervals ) );
  };

 // a datum indexed over the instants only
 auto read_t = [ & ]( const std::string & name ,
                      std::vector< double > & data ) {
  ::deserialize( group , name , T , data , true , true ,
                 v_change_intervals );
  };

 // the bounds- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 read( "MaxPower" , f_G , v_MaxPower , false );
 if( ! read( "MinPower" , f_G , v_MinPower , true ) )
  v_MinPower.resize( boost::extents[ T ][ f_G ] );

 // the operating rows - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( f_M ) {
  ::deserialize( group , "OperatingMatrix" , { f_M , f_G } ,
                 v_OperatingMatrix , false , false );
  read( "OperatingLHS" , f_M , v_OperatingLHS , true );
  read( "OperatingRHS" , f_M , v_OperatingRHS , true );
  }
 else {
  v_OperatingMatrix.resize( boost::extents[ 0 ][ 0 ] );
  v_OperatingLHS.resize( boost::extents[ 0 ][ 0 ] );
  v_OperatingRHS.resize( boost::extents[ 0 ][ 0 ] );
  }

 // the ramps- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 v_RampMatrix.resize( boost::extents[ 0 ][ 0 ] );
 if( ! group.getVar( "RampMatrix" ).isNull() ) {
  if( ! deserialize_dim( group , "NumberRampRows" , f_K ) )
   throw( std::invalid_argument( who + ": RampMatrix without "
                                 "NumberRampRows" ) );
  ::deserialize( group , "RampMatrix" , { f_K , f_G } , v_RampMatrix ,
                 false , false );
  }
 else
  f_K = f_G;

 if( f_K ) {
  read( "DeltaRampUp" , f_K , v_DeltaRampUp , true );
  read( "DeltaRampDown" , f_K , v_DeltaRampDown , true );
  }
 else {
  v_DeltaRampUp.resize( boost::extents[ 0 ][ 0 ] );
  v_DeltaRampDown.resize( boost::extents[ 0 ][ 0 ] );
  }

 // the start-up and shut-down limits- - - - - - - - - - - - - - - - - - - -

 read( "StartUpLimit" , f_G , v_StartUpLimit , true );
 read( "StartUpLowerLimit" , f_G , v_StartUpLowerLimit , true );
 read( "ShutDownLimit" , f_G , v_ShutDownLimit , true );
 read( "ShutDownLowerLimit" , f_G , v_ShutDownLowerLimit , true );

 // the reserves - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 auto any_nonzero = []( const MAdbl & a ) {
  return( std::any_of( a.data() , a.data() + a.num_elements() ,
                       []( double v ) { return( v != 0 ); } ) );
  };

 read( "PrimaryRho" , f_G , v_PrimaryRho , true );
 if( ! any_nonzero( v_PrimaryRho ) )
  v_PrimaryRho.resize( boost::extents[ 0 ][ 0 ] );
 read( "SecondaryRho" , f_G , v_SecondaryRho , true );
 if( ! any_nonzero( v_SecondaryRho ) )
  v_SecondaryRho.resize( boost::extents[ 0 ][ 0 ] );

 ::deserialize( group , "ReserveDirection" , { f_G , f_G } ,
                v_ReserveDirection , true , false );

 read( "PrimarySpinningReserveCost" , f_G , v_PrimaryCost , true );
 read( "SecondarySpinningReserveCost" , f_G , v_SecondaryCost , true );

 // the costs- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 read( "LinearTerm" , f_G , v_LinearTerm , true );
 read( "QuadTerm" , f_G , v_QuadTerm , true );
 read_t( "ConstTerm" , v_ConstTerm );
 read_t( "StartUpCost" , v_StartUpCost );
 read_t( "ShutDownCost" , v_ShutDownCost );

 // the terms of the linking rows of the UCBlock - - - - - - - - - - - - - -

 read_t( "FixedConsumption" , v_FixedConsumption );
 read_t( "InertiaCommitment" , v_InertiaCommitment );

 MAdbl ip;
 v_InertiaPower.assign( f_G , {} );
 if( read( "InertiaPower" , f_G , ip , true ) )
  for( Index g = 0 ; g < f_G ; ++g ) {
   bool nz = false;
   for( Index t = 0 ; t < T ; ++t )
    nz |= ( ip[ t ][ g ] != 0 );
   if( nz ) {
    v_InertiaPower[ g ].resize( T );
    for( Index t = 0 ; t < T ; ++t )
     v_InertiaPower[ g ][ t ] = ip[ t ][ g ];
    }
   }

 // the initial state and the minimum times- - - - - - - - - - - - - - - - -

 if( ! ::deserialize( group , "InitialPower" , f_G , v_InitialPower , true ,
                      true ) )
  v_InitialPower.assign( f_G , 0 );
 else
  if( v_InitialPower.size() == 1 )
   v_InitialPower.resize( f_G , v_InitialPower.front() );

 f_MinUpTime = f_MinDownTime = 1;
 if( ::deserialize( group , f_MinUpTime , "MinUpTime" ) )
  f_MinUpTime = std::max( f_MinUpTime , Index( 1 ) );
 if( ::deserialize( group , f_MinDownTime , "MinDownTime" ) )
  f_MinDownTime = std::max( f_MinDownTime , Index( 1 ) );

 if( ! ::deserialize( group , f_InitUpDownTime , "InitUpDownTime" ) ) {
  const bool zero = std::all_of( v_InitialPower.begin() ,
                                 v_InitialPower.end() ,
                                 []( double p ) { return( p == 0 ); } );
  f_InitUpDownTime = zero ? - int( f_MinDownTime ) : int( f_MinUpTime );
  }

 f_MinUpTime = clamp_min_time( f_MinUpTime , T , f_InitUpDownTime , true );
 f_MinDownTime = clamp_min_time( f_MinDownTime , T , f_InitUpDownTime ,
                                 false );

 f_scale = 1;
 ::deserialize( group , f_scale , "Scale" );

 compute_reserve_signs( who );
 check_data( who );

 }  // end( ConversionUnitBlock::deserialize )

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG

std::vector< std::string > ConversionUnitBlock::expected_dims( void ) const
{
 static const std::vector< std::string > ed =
  { "NumberGenerators" , "NumberOperatingRows" , "NumberRampRows" };

 auto ret = UnitBlock::expected_dims();
 ret.insert( ret.end() , ed.begin() , ed.end() );
 return( ret );
 }

/*--------------------------------------------------------------------------*/

std::vector< std::string > ConversionUnitBlock::expected_vars( void ) const
{
 static const std::vector< std::string > ev =
  { "CommitmentGenerator" , "MinPower" , "MaxPower" , "OperatingMatrix" ,
    "OperatingLHS" , "OperatingRHS" , "RampMatrix" , "DeltaRampUp" ,
    "DeltaRampDown" , "StartUpLimit" , "StartUpLowerLimit" ,
    "ShutDownLimit" , "ShutDownLowerLimit" , "PrimaryRho" , "SecondaryRho" ,
    "ReserveDirection" , "PrimarySpinningReserveCost" ,
    "SecondarySpinningReserveCost" , "LinearTerm" , "QuadTerm" ,
    "ConstTerm" , "StartUpCost" , "ShutDownCost" , "FixedConsumption" ,
    "InertiaCommitment" , "InertiaPower" , "InitialPower" ,
    "InitUpDownTime" , "MinUpTime" , "MinDownTime" , "Scale" };

 auto ret = UnitBlock::expected_vars();
 ret.insert( ret.end() , ev.begin() , ev.end() );
 return( ret );
 }

#endif

/*--------------------------------------------------------------------------*/
/*------------------------ CHECKING THE DATA -------------------------------*/
/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::compute_reserve_signs( const std::string & who )
{
 const Index T = f_time_horizon;
 v_sign.assign( f_G , 0 );
 f_has_primary = f_has_secondary = false;
 for( Index g = 0 ; g < f_G ; ++g ) {
  bool pr = false , sc = false;
  for( Index t = 0 ; t < T ; ++t ) {
   pr |= ( ! absent( v_PrimaryRho ) ) && ( v_PrimaryRho[ t ][ g ] != 0 );
   sc |= ( ! absent( v_SecondaryRho ) ) &&
         ( v_SecondaryRho[ t ][ g ] != 0 );
   }
  f_has_primary |= pr;
  f_has_secondary |= sc;
  if( ! ( pr || sc ) )
   continue;
  bool nonneg = true , nonpos = true;
  for( Index t = 0 ; t < T ; ++t ) {
   nonneg &= ( v_MinPower[ t ][ g ] >= 0 );
   nonpos &= ( v_MaxPower[ t ][ g ] <= 0 );
   }
  if( nonneg )
   v_sign[ g ] = 1;
  else
   if( nonpos )
    v_sign[ g ] = -1;
   else
    throw( std::invalid_argument( who + ": generator " +
                                  std::to_string( g ) + " offers a reserve "
                                  "and its power may have either sign" ) );
  }
 }

/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::check_data( const std::string & who ) const
{
 const Index T = f_time_horizon;
 auto fail = [ & ]( const std::string & what ) {
  throw( std::invalid_argument( who + ": " + what ) );
  };
 auto tg = [ & ]( Index t , Index g ) {
  return( " at instant " + std::to_string( t ) + " of generator " +
          std::to_string( g ) );
  };

 for( Index t = 0 ; t < T ; ++t )
  for( Index g = 0 ; g < f_G ; ++g ) {
   const double mn = v_MinPower[ t ][ g ];
   const double mx = v_MaxPower[ t ][ g ];
   if( ( ! std::isfinite( mn ) ) || ( ! std::isfinite( mx ) ) )
    fail( "a bound is not finite" + tg( t , g ) );
   if( mn > mx )
    fail( "MinPower > MaxPower" + tg( t , g ) );
   if( ( ! absent( v_StartUpLimit ) ) && ( ! absent( v_StartUpLowerLimit ) )
       && ( v_StartUpLowerLimit[ t ][ g ] > v_StartUpLimit[ t ][ g ] ) )
    fail( "StartUpLowerLimit > StartUpLimit" + tg( t , g ) );
   if( ( ! absent( v_ShutDownLimit ) ) &&
       ( ! absent( v_ShutDownLowerLimit ) ) &&
       ( v_ShutDownLowerLimit[ t ][ g ] > v_ShutDownLimit[ t ][ g ] ) )
    fail( "ShutDownLowerLimit > ShutDownLimit" + tg( t , g ) );
   if( ( ! absent( v_QuadTerm ) ) && ( v_QuadTerm[ t ][ g ] < 0 ) )
    fail( "QuadTerm < 0" + tg( t , g ) );
   for( const auto * rho : { & v_PrimaryRho , & v_SecondaryRho } )
    if( ( ! absent( * rho ) ) && ( ( * rho )[ t ][ g ] < 0 ) )
     fail( "a reserve fraction is negative" + tg( t , g ) );
   // a reserve generator keeps the sign it had when the reserve rows were
   // written [see compute_reserve_signs()]
   if( ( v_sign.size() == f_G ) && v_sign[ g ] ) {
    if( ( v_sign[ g ] > 0 ) && ( mn < 0 ) )
     fail( "a nonnegative reserve generator gets MinPower < 0" +
           tg( t , g ) );
    if( ( v_sign[ g ] < 0 ) && ( mx > 0 ) )
     fail( "a nonpositive reserve generator gets MaxPower > 0" +
           tg( t , g ) );
    }
   }

 for( Index m = 0 ; m < f_M ; ++m ) {
  bool finite = false;
  for( Index t = 0 ; t < T ; ++t ) {
   const double lo = get_operating_lhs( t , m );
   const double hi = get_operating_rhs( t , m );
   if( lo > hi )
    fail( "OperatingLHS > OperatingRHS at instant " + std::to_string( t ) +
          " of row " + std::to_string( m ) );
   if( ( lo == INF ) || ( hi == - INF ) )
    fail( "an infinite side of the wrong sign in row " +
          std::to_string( m ) );
   finite |= ( lo > - INF ) || ( hi < INF );
   }
  if( ! finite )
   fail( "operating row " + std::to_string( m ) + " has no finite side" );
  }

 for( const auto * d : { & v_DeltaRampUp , & v_DeltaRampDown } )
  if( ! absent( * d ) )
   for( Index i = 0 ; i < d->num_elements() ; ++i )
    if( ! ( d->data()[ i ] >= 0 ) )
     fail( "a ramp is negative or not a number" );

 if( ! absent( v_ReserveDirection ) )
  for( Index g = 0 ; g < f_G ; ++g )
   if( v_sign[ g ] && ( v_ReserveDirection[ g ][ g ] != 1 ) )
    fail( "ReserveDirection[ " + std::to_string( g ) + " ][ " +
          std::to_string( g ) + " ] is not 1 for a reserve generator" );

 for( auto fc : v_FixedConsumption )
  if( fc < 0 )
   fail( "FixedConsumption < 0" );

 if( v_InitialPower.size() != f_G )
  fail( "InitialPower has not NumberGenerators entries" );

 // a unit on before the horizon is in the region of the instant 0
 if( f_InitUpDownTime > 0 && T ) {
  auto eps = []( double v ) { return( 1e-9 * std::max( 1.0 ,
                                                       std::abs( v ) ) ); };
  for( Index g = 0 ; g < f_G ; ++g ) {
   const double p = v_InitialPower[ g ];
   if( ( p < v_MinPower[ 0 ][ g ] - eps( p ) ) ||
       ( p > v_MaxPower[ 0 ][ g ] + eps( p ) ) )
    fail( "InitialPower of generator " + std::to_string( g ) + " is "
          "outside its bounds at the instant 0 of a unit on before the "
          "horizon" );
   }
  for( Index m = 0 ; m < f_M ; ++m ) {
   double ap = 0;
   for( Index g = 0 ; g < f_G ; ++g )
    ap += v_OperatingMatrix[ m ][ g ] * v_InitialPower[ g ];
   if( ( ap < get_operating_lhs( 0 , m ) - eps( ap ) ) ||
       ( ap > get_operating_rhs( 0 , m ) + eps( ap ) ) )
    fail( "InitialPower violates the operating row " + std::to_string( m ) +
          " at the instant 0 of a unit on before the horizon" );
   }
  }
 }

/*--------------------------------------------------------------------------*/
/*------------------------ READING THE DATA --------------------------------*/
/*--------------------------------------------------------------------------*/

double ConversionUnitBlock::get_operating_lhs( Index t , Index m ) const
{
 return( absent( v_OperatingLHS ) ? - INF : v_OperatingLHS[ t ][ m ] );
 }

/*--------------------------------------------------------------------------*/

double ConversionUnitBlock::get_operating_rhs( Index t , Index m ) const
{
 return( absent( v_OperatingRHS ) ? INF : v_OperatingRHS[ t ][ m ] );
 }

/*--------------------------------------------------------------------------*/

double ConversionUnitBlock::get_ramp_coefficient( Index k , Index g ) const
{
 if( absent( v_RampMatrix ) )
  return( k == g ? 1 : 0 );
 return( v_RampMatrix[ k ][ g ] );
 }

/*--------------------------------------------------------------------------*/

Index ConversionUnitBlock::first_free_instant( void ) const
{
 const int init = f_InitUpDownTime;
 Index t0;
 if( init > 0 )
  t0 = ( Index( init ) >= f_MinUpTime ? 0 : f_MinUpTime - Index( init ) );
 else
  t0 = ( Index( - init ) >= f_MinDownTime ? 0 :
         f_MinDownTime - Index( - init ) );
 return( std::min( t0 , f_time_horizon ) );
 }

/*--------------------------------------------------------------------------*/

bool ConversionUnitBlock::has_upper_limit( Index ) const
{
 return( ( ! absent( v_StartUpLimit ) ) || ( ! absent( v_ShutDownLimit ) ) );
 }

/*--------------------------------------------------------------------------*/

bool ConversionUnitBlock::has_lower_limit( Index ) const
{
 return( ( ! absent( v_StartUpLowerLimit ) ) ||
         ( ! absent( v_ShutDownLowerLimit ) ) );
 }

/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::limit_box( bool su , Index t , Index g ,
                                     double & lo , double & hi ) const
{
 // the power the limits refer to: that of t at a start-up at t, that of
 // t - 1 at a shut-down at t (that of 0 if t = 0)
 const Index tt = ( su || ( t == 0 ) ) ? t : t - 1;
 lo = v_MinPower[ tt ][ g ];
 hi = v_MaxPower[ tt ][ g ];
 const auto & up = su ? v_StartUpLimit : v_ShutDownLimit;
 const auto & dn = su ? v_StartUpLowerLimit : v_ShutDownLowerLimit;
 if( ! absent( up ) )
  hi = std::min( hi , up[ t ][ g ] );
 if( ! absent( dn ) )
  lo = std::max( lo , dn[ t ][ g ] );
 }

/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::combination_range( bool su , Index t , Index k ,
                                             double & lo ,
                                             double & hi ) const
{
 lo = hi = 0;
 for( Index g = 0 ; g < f_G ; ++g ) {
  const double c = get_ramp_coefficient( k , g );
  if( c == 0 )
   continue;
  double l , h;
  limit_box( su , t , g , l , h );
  lo += std::min( c * l , c * h );
  hi += std::max( c * l , c * h );
  }
 }

/*--------------------------------------------------------------------------*/
/*---------------------- GENERATING THE ABSTRACT PART ----------------------*/
/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::generate_abstract_variables( Configuration * stvv )
{
 if( variables_generated() )
  return;

 UnitBlock::generate_abstract_variables( stvv );
 const Index T = f_time_horizon;

 v_commitment.resize( T );
 for( auto & var : v_commitment )
  var.set_type( ColVariable::kBinary );
 add_static_variable( v_commitment , "u_conversion" );

 v_start_up.resize( T );
 for( auto & var : v_start_up )
  var.set_type( ColVariable::kBinary );
 add_static_variable( v_start_up , "v_conversion" );

 v_shut_down.resize( T );
 for( auto & var : v_shut_down )
  var.set_type( ColVariable::kBinary );
 add_static_variable( v_shut_down , "w_conversion" );

 v_active_power.resize( f_G );
 for( auto & pg : v_active_power ) {
  pg.resize( T );
  for( auto & var : pg )
   var.set_type( ColVariable::kContinuous );
  }
 add_static_variable( v_active_power , "p_conversion" );

 // the reserve Variable of the generators that offer it, if the UCBlock
 // asks for that reserve
 auto reserve = [ & ]( const MAdbl & rho ,
                       std::vector< std::vector< ColVariable > > & r ,
                       std::string && name ) {
  r.assign( f_G , {} );
  if( absent( rho ) )
   return;
  for( Index g = 0 ; g < f_G ; ++g ) {
   bool nz = false;
   for( Index t = 0 ; t < T ; ++t )
    nz |= ( rho[ t ][ g ] != 0 );
   if( nz ) {
    r[ g ].resize( T );
    for( auto & var : r[ g ] )
     var.set_type( ColVariable::kNonNegative );
    }
   }
  add_static_variable( r , std::move( name ) );
  };

 v_primary_reserve.clear();
 v_secondary_reserve.clear();
 if( reserve_vars & 1u )
  reserve( v_PrimaryRho , v_primary_reserve , "pr_conversion" );
 if( reserve_vars & 2u )
  reserve( v_SecondaryRho , v_secondary_reserve , "sc_conversion" );

 set_variables_generated();
 }

/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::generate_abstract_constraints( Configuration * stcc )
{
 if( constraints_generated() )
  return;

 bool zo = false;
 if( ( ! stcc ) && f_BlockConfig )
  stcc = f_BlockConfig->f_static_constraints_Configuration;
 if( auto sci = dynamic_cast< SimpleConfiguration< int > * >( stcc ) )
  zo = sci->f_value;

 // the rows are written in two passes, the first counting them, since a
 // group cannot grow once it is filled (a Constraint is not copied)
 f_row_count.clear();
 f_counting = true;
 try {
  build_rows();
  }
 catch( ... ) {
  f_counting = false;
  throw;
  }
 f_counting = false;
 for( auto & [ rows , n ] : f_row_count )
  rows->resize( n );
 f_row_count.clear();  // now the position of the next row of each group
 build_rows();

 if( zo ) {
  StartUp_Binary_Bound.resize( f_time_horizon );
  ShutDown_Binary_Bound.resize( f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   StartUp_Binary_Bound[ t ].set_variable( & v_start_up[ t ] );
   ShutDown_Binary_Bound[ t ].set_variable( & v_shut_down[ t ] );
   }
  add_static_constraint( StartUp_Binary_Bound ,
                         "StartUp_Binary_Bound_Conversion" );
  add_static_constraint( ShutDown_Binary_Bound ,
                         "ShutDown_Binary_Bound_Conversion" );
  }

 set_constraints_generated();
 }

/*--------------------------------------------------------------------------*/
/* What update_rows() collects while build_rows() compares the rows: the
 * changes of the rows and of the boxes, and the number of rows written in
 * each group. */

struct ConversionUnitBlock::RowCmp {
 struct RowChange {
  FRowConstraint * row;
  Subset idx;
  std::vector< double > coef;
  double lhs;
  double rhs;
  };
 struct BoxChange {
  BoxConstraint * box;
  double lhs;
  double rhs;
  };

 std::vector< RowChange > rows;
 std::vector< BoxChange > boxes;
 std::map< std::vector< FRowConstraint > * , Index > pushed;
 };

/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::push_row( std::vector< FRowConstraint > & rows ,
                                    LinearFunction::v_coeff_pair && vars ,
                                    double lhs , double rhs )
{
 if( f_counting ) {
  if( std::find( f_row_groups.begin() , f_row_groups.end() , & rows ) ==
      f_row_groups.end() )
   f_row_groups.push_back( & rows );
  ++f_row_count[ & rows ];
  vars.clear();
  return;
  }

 if( generating_rows() ) {
  auto & row = rows[ f_row_count[ & rows ]++ ];
  row.set_lhs( lhs );
  row.set_rhs( rhs );
  row.set_function( new LinearFunction( std::move( vars ) ) );
  return;
  }

 const Index i = f_row_cmp->pushed[ & rows ]++;
 if( i >= rows.size() )
  throw( std::logic_error( "ConversionUnitBlock::update_rows: a group of " +
                           std::to_string( rows.size() ) + " rows would "
                           "have more of them" ) );

 auto & row = rows[ i ];
 const auto & old = static_cast< const LinearFunction * >(
                                         row.get_function() )->get_v_var();
 if( old.size() != vars.size() )
  throw( std::logic_error( "ConversionUnitBlock::update_rows: a row would "
                           "have another number of Variable" ) );

 RowCmp::RowChange c{ & row , {} , {} , lhs , rhs };
 for( Index k = 0 ; k < vars.size() ; ++k ) {
  if( old[ k ].first != vars[ k ].first )
   throw( std::logic_error( "ConversionUnitBlock::update_rows: a row would "
                            "have other Variable" ) );
  if( old[ k ].second != vars[ k ].second ) {
   c.idx.push_back( k );
   c.coef.push_back( vars[ k ].second );
   }
  }

 if( ( ! c.idx.empty() ) || ( row.get_lhs() != lhs ) ||
     ( row.get_rhs() != rhs ) )
  f_row_cmp->rows.push_back( std::move( c ) );

 vars.clear();
 }

/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::put_box( BoxConstraint & box , ColVariable * var ,
                                   double lhs , double rhs )
{
 if( f_counting )
  return;
 if( generating_rows() ) {
  box.set_lhs( lhs );
  box.set_rhs( rhs );
  box.set_variable( var );
  }
 else
  if( ( box.get_lhs() != lhs ) || ( box.get_rhs() != rhs ) )
   f_row_cmp->boxes.push_back( { & box , lhs , rhs } );
 }

/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::add_rows( std::vector< FRowConstraint > & rows ,
                                    std::string && name )
{
 if( generating_rows() && ( ! f_counting ) && ( ! rows.empty() ) )
  add_static_constraint( rows , std::move( name ) );
 }

/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::build_rows( void )
{
 const Index T = f_time_horizon;
 const Index t0 = first_free_instant();
 const bool on0 = ( f_InitUpDownTime > 0 );
 auto p_init = [ & ]( Index g ) { return( on0 ? v_InitialPower[ g ] : 0 ); };

 LinearFunction::v_coeff_pair vars;
 auto add = [ & vars ]( ColVariable & v , double c ) {
  vars.push_back( std::make_pair( & v , c ) );
  };
 auto pos = []( double v ) { return( std::max( v , 0.0 ) ); };

 // the reserve generators and the direction matrix- - - - - - - - - - - - -
 std::vector< Index > R;
 for( Index g = 0 ; g < f_G ; ++g )
  if( get_primary_spinning_reserve( g ) ||
      get_secondary_spinning_reserve( g ) )
   R.push_back( g );
 auto D = [ & ]( Index h , Index g ) -> double {
  if( absent( v_ReserveDirection ) )
   return( h == g ? 1 : 0 );
  return( v_ReserveDirection[ h ][ g ] );
  };
 // adds coef times r[ t ][ g ] = coef ( pr[ t ][ g ] + sc[ t ][ g ] )
 auto add_r = [ & ]( Index t , Index g , double coef ) {
  if( auto pr = get_primary_spinning_reserve( g ) )
   add( pr[ t ] , coef );
  if( auto sc = get_secondary_spinning_reserve( g ) )
   add( sc[ t ] , coef );
  };
 // whether s * w[ g ] has a positive entry over the reserve generators
 auto pos_entry = [ & ]( const std::vector< double > & w , double s ) {
  return( std::any_of( R.begin() , R.end() ,
                       [ & ]( Index g ) { return( s * w[ g ] > 0 ); } ) );
  };

 // (1)-(3), the commitment- - - - - - - - - - - - - - - - - - - - - - - - -

 for( Index t = 0 ; t < T ; ++t ) {
  add( v_commitment[ t ] , 1 );
  if( t )
   add( v_commitment[ t - 1 ] , -1 );
  add( v_start_up[ t ] , -1 );
  add( v_shut_down[ t ] , 1 );
  const double rhs = ( t == 0 ) && on0 ? 1 : 0;
  push_row( Logical_Const , std::move( vars ) , rhs , rhs );
  }
 add_rows( Logical_Const , "Logical_Const_Conversion" );

 for( Index t = 0 ; t < T ; ++t ) {
  for( Index s = ( t + 1 > f_MinUpTime ? t + 1 - f_MinUpTime : 0 ) ;
       s <= t ; ++s )
   add( v_start_up[ s ] , 1 );
  add( v_commitment[ t ] , -1 );
  push_row( MinUp_Const , std::move( vars ) , - INF , 0 );
  }
 add_rows( MinUp_Const , "MinUp_Const_Conversion" );

 for( Index t = 0 ; t < T ; ++t ) {
  for( Index s = ( t + 1 > f_MinDownTime ? t + 1 - f_MinDownTime : 0 ) ;
       s <= t ; ++s )
   add( v_shut_down[ s ] , 1 );
  add( v_commitment[ t ] , 1 );
  push_row( MinDown_Const , std::move( vars ) , - INF , 1 );
  }
 add_rows( MinDown_Const , "MinDown_Const_Conversion" );

 // the bounds of u and w[ 0 ]- - - - - - - - - - - - - - - - - - - - - - -

 if( generating_rows() )
  Commitment_Bound.resize( T );
 for( Index t = 0 ; t < T ; ++t )
  put_box( Commitment_Bound[ t ] , & v_commitment[ t ] ,
           ( on0 && ( t < t0 ) ) ? 1 : 0 ,
           ( ( ! on0 ) && ( t < t0 ) ) ? 0 : 1 );
 if( generating_rows() && ( ! f_counting ) && T )
  add_static_constraint( Commitment_Bound , "Commitment_Bound_Conversion" );

 if( T ) {
  bool sd0 = true;  // whether the unit may shut down at 0
  if( on0 )
   for( Index g = 0 ; g < f_G ; ++g ) {
    double lo , hi;
    limit_box( false , 0 , g , lo , hi );
    const double p = p_init( g );
    if( ( p < lo - 1e-9 * std::max( 1.0 , std::abs( lo ) ) ) ||
        ( p > hi + 1e-9 * std::max( 1.0 , std::abs( hi ) ) ) )
     sd0 = false;
    }
  put_box( ShutDownZero_Bound , & v_shut_down[ 0 ] , 0 , sd0 ? 1 : 0 );
  if( generating_rows() && ( ! f_counting ) )
   add_static_constraint( ShutDownZero_Bound ,
                          "ShutDownZero_Bound_Conversion" );
  }

 // (4)-(7), the bounds of the powers - - - - - - - - - - - - - - - - - - - -

 auto hat = [ & ]( bool upper , bool su , Index t , Index g ) {
  // U^su, U^sd (upper) or L^su, L^sd at t, the shut-down one being of the
  // shut-down at t + 1 (the bound itself at T - 1)
  if( ( ! su ) && ( t + 1 >= T ) )
   return( upper ? v_MaxPower[ t ][ g ] : v_MinPower[ t ][ g ] );
  double lo , hi;
  limit_box( su , su ? t : t + 1 , g , lo , hi );
  return( upper ? hi : lo );
  };

 for( Index t = 0 ; t < T ; ++t )
  for( Index g = 0 ; g < f_G ; ++g ) {
   const double mn = v_MinPower[ t ][ g ];
   const double mx = v_MaxPower[ t ][ g ];
   const bool last = ( t + 1 >= T );
   auto & p = v_active_power[ g ][ t ];

   if( ! has_lower_limit( g ) ) {
    add( p , 1 );
    add( v_commitment[ t ] , - mn );
    push_row( MinPower_Const , std::move( vars ) , 0 , INF );
    }
   else {
    const double ls = hat( false , true , t , g );
    const double ld = hat( false , false , t , g );
    if( f_MinUpTime >= 2 ) {
     add( p , 1 );
     add( v_commitment[ t ] , - mn );
     add( v_start_up[ t ] , - ( ls - mn ) );
     if( ! last )
      add( v_shut_down[ t + 1 ] , - ( ld - mn ) );
     push_row( MinPower_Const , std::move( vars ) , 0 , INF );
     }
    else {
     if( ! last ) {
      add( p , 1 );
      add( v_commitment[ t ] , - mn );
      add( v_shut_down[ t + 1 ] , - ( ld - mn ) );
      add( v_start_up[ t ] , - pos( ls - ld ) );
      push_row( MinPower_Const , std::move( vars ) , 0 , INF );
      }
     add( p , 1 );
     add( v_commitment[ t ] , - mn );
     if( ! last )
      add( v_shut_down[ t + 1 ] , - pos( ld - ls ) );
     add( v_start_up[ t ] , - ( ls - mn ) );
     push_row( MinPower_Const , std::move( vars ) , 0 , INF );
     }
    }

   if( ! has_upper_limit( g ) ) {
    add( p , 1 );
    add( v_commitment[ t ] , - mx );
    push_row( MaxPower_Const , std::move( vars ) , - INF , 0 );
    }
   else {
    const double us = hat( true , true , t , g );
    const double ud = hat( true , false , t , g );
    if( f_MinUpTime >= 2 ) {
     add( p , 1 );
     add( v_commitment[ t ] , - mx );
     add( v_start_up[ t ] , mx - us );
     if( ! last )
      add( v_shut_down[ t + 1 ] , mx - ud );
     push_row( MaxPower_Const , std::move( vars ) , - INF , 0 );
     }
    else {
     if( ! last ) {
      add( p , 1 );
      add( v_commitment[ t ] , - mx );
      add( v_shut_down[ t + 1 ] , mx - ud );
      add( v_start_up[ t ] , pos( ud - us ) );
      push_row( MaxPower_Const , std::move( vars ) , - INF , 0 );
      }
     add( p , 1 );
     add( v_commitment[ t ] , - mx );
     if( ! last )
      add( v_shut_down[ t + 1 ] , pos( us - ud ) );
     add( v_start_up[ t ] , mx - us );
     push_row( MaxPower_Const , std::move( vars ) , - INF , 0 );
     }
    }
   }
 add_rows( MinPower_Const , "MinPower_Const_Conversion" );
 add_rows( MaxPower_Const , "MaxPower_Const_Conversion" );

 if( generating_rows() ) {
  ActivePower_Bound.resize( f_G );
  for( auto & b : ActivePower_Bound )
   b.resize( T );
  }
 for( Index g = 0 ; g < f_G ; ++g )
  for( Index t = 0 ; t < T ; ++t )
   put_box( ActivePower_Bound[ g ][ t ] , & v_active_power[ g ][ t ] ,
            std::min( 0.0 , double( v_MinPower[ t ][ g ] ) ) ,
            std::max( 0.0 , double( v_MaxPower[ t ][ g ] ) ) );
 if( generating_rows() && ( ! f_counting ) && T )
  add_static_constraint( ActivePower_Bound ,
                         "ActivePower_Bound_Conversion" );

 // (8), the operating rows - - - - - - - - - - - - - - - - - - - - - - - - -

 auto add_ap = [ & ]( Index t , Index m ) {
  for( Index g = 0 ; g < f_G ; ++g )
   if( v_OperatingMatrix[ m ][ g ] != 0 )
    add( v_active_power[ g ][ t ] , v_OperatingMatrix[ m ][ g ] );
  };

 for( Index t = 0 ; t < T ; ++t )
  for( Index m = 0 ; m < f_M ; ++m ) {
   const double lo = get_operating_lhs( t , m );
   const double hi = get_operating_rhs( t , m );
   if( lo == hi ) {
    add_ap( t , m );
    add( v_commitment[ t ] , - hi );
    push_row( Operating_Const , std::move( vars ) , 0 , 0 );
    continue;
    }
   if( hi < INF ) {
    add_ap( t , m );
    add( v_commitment[ t ] , - hi );
    push_row( Operating_Const , std::move( vars ) , - INF , 0 );
    }
   if( lo > - INF ) {
    add_ap( t , m );
    add( v_commitment[ t ] , - lo );
    push_row( Operating_Const , std::move( vars ) , 0 , INF );
    }
   }
 add_rows( Operating_Const , "Operating_Const_Conversion" );

 // (9), the ramps- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // the terms in p of y[ t ] - y[ t - 1 ] times sgn, plus the reserve terms
 // rs * sum_h c[ k ][ h ] D[ h ][ g ] r[ t ][ g ] if rs != 0; returns the
 // constant of y[ -1 ] times sgn at t = 0
 auto add_dy = [ & ]( Index t , Index k , double sgn , double rs ) {
  double cst = 0;
  for( Index g = 0 ; g < f_G ; ++g ) {
   const double c = get_ramp_coefficient( k , g );
   if( c == 0 )
    continue;
   add( v_active_power[ g ][ t ] , sgn * c );
   if( t )
    add( v_active_power[ g ][ t - 1 ] , - sgn * c );
   else
    cst += sgn * c * p_init( g );
   }
  if( rs != 0 )
   for( auto g : R ) {
    double cd = 0;
    for( Index h = 0 ; h < f_G ; ++h )
     cd += get_ramp_coefficient( k , h ) * D( h , g );
    if( cd != 0 )
     add_r( t , g , rs * cd );
    }
  return( cst );
  };

 // the up row (9) of k at t, with the reserve terms times rs
 auto ramp_up = [ & ]( std::vector< FRowConstraint > & rows , Index t ,
                       Index k , double rs ) {
  double sl , sh , dl , dh;
  combination_range( true , t , k , sl , sh );
  combination_range( false , t , k , dl , dh );
  const double dlt = v_DeltaRampUp[ t ][ k ];
  const double cst = add_dy( t , k , 1 , rs );
  add( v_commitment[ t ] , - dlt );
  add( v_start_up[ t ] , dlt - sh );
  add( v_shut_down[ t ] , dl );
  push_row( rows , std::move( vars ) , - INF , cst );
  };

 // the down row (9) of k at t, with the reserve terms times rs
 auto ramp_dn = [ & ]( std::vector< FRowConstraint > & rows , Index t ,
                       Index k , double rs ) {
  double sl , sh , dl , dh;
  combination_range( true , t , k , sl , sh );
  combination_range( false , t , k , dl , dh );
  const double dlt = v_DeltaRampDown[ t ][ k ];
  const double cst = add_dy( t , k , -1 , rs );
  add( v_commitment[ t ] , - dlt );
  add( v_start_up[ t ] , dlt + sl );
  add( v_shut_down[ t ] , - dh );
  push_row( rows , std::move( vars ) , - INF , cst );
  };

 if( ! absent( v_DeltaRampUp ) ) {
  for( Index t = 0 ; t < T ; ++t )
   for( Index k = 0 ; k < f_K ; ++k )
    ramp_up( RampUp_Const , t , k , 0 );
  add_rows( RampUp_Const , "RampUp_Const_Conversion" );
  }

 if( ! absent( v_DeltaRampDown ) ) {
  for( Index t = 0 ; t < T ; ++t )
   for( Index k = 0 ; k < f_K ; ++k )
    ramp_dn( RampDown_Const , t , k , 0 );
  add_rows( RampDown_Const , "RampDown_Const_Conversion" );
  }

 if( R.empty() )
  return;

 // (10), the reserve fractions - - - - - - - - - - - - - - - - - - - - - - -

 auto rho_rows = [ & ]( const MAdbl & rho , std::vector< ColVariable > * r ,
                        Index g , std::vector< FRowConstraint > & rows ) {
  for( Index t = 0 ; t < T ; ++t ) {
   add( ( * r )[ t ] , 1 );
   add( v_active_power[ g ][ t ] , - rho[ t ][ g ] * v_sign[ g ] );
   push_row( rows , std::move( vars ) , - INF , 0 );
   }
  };

 for( auto g : R ) {
  if( ( g < v_primary_reserve.size() ) && ! v_primary_reserve[ g ].empty() )
   rho_rows( v_PrimaryRho , & v_primary_reserve[ g ] , g , PrimaryRho_Const );
  if( ( g < v_secondary_reserve.size() ) &&
      ! v_secondary_reserve[ g ].empty() )
   rho_rows( v_SecondaryRho , & v_secondary_reserve[ g ] , g ,
             SecondaryRho_Const );
  }
 add_rows( PrimaryRho_Const , "PrimaryRho_Const_Conversion" );
 add_rows( SecondaryRho_Const , "SecondaryRho_Const_Conversion" );

 // (11), the bounds of the deployed powers - - - - - - - - - - - - - - - - -

 // the terms of q^s[ t ][ h ]
 auto add_q = [ & ]( Index t , Index h , double s ) {
  add( v_active_power[ h ][ t ] , 1 );
  for( auto g : R )
   if( D( h , g ) != 0 )
    add_r( t , g , s * D( h , g ) );
  };

 for( Index t = 0 ; t < T ; ++t )
  for( double s : { 1.0 , -1.0 } )
   for( Index h = 0 ; h < f_G ; ++h ) {
    std::vector< double > dh( f_G , 0 );
    for( auto g : R )
     dh[ g ] = D( h , g );
    const double mn = v_MinPower[ t ][ h ];
    const double mx = v_MaxPower[ t ][ h ];
    const bool last = ( t + 1 >= T );

    if( pos_entry( dh , s ) ) {  // the upper side
     if( ! has_upper_limit( h ) ) {
      add_q( t , h , s );
      add( v_commitment[ t ] , - mx );
      push_row( ReserveBound_Const , std::move( vars ) , - INF , 0 );
      }
     else {
      add_q( t , h , s );
      add( v_commitment[ t ] , - mx );
      add( v_start_up[ t ] , mx - hat( true , true , t , h ) );
      push_row( ReserveBound_Const , std::move( vars ) , - INF , 0 );
      if( ! last ) {
       add_q( t , h , s );
       add( v_commitment[ t ] , - mx );
       add( v_shut_down[ t + 1 ] , mx - hat( true , false , t , h ) );
       push_row( ReserveBound_Const , std::move( vars ) , - INF , 0 );
       }
      }
     }

    if( pos_entry( dh , - s ) ) {  // the lower side
     if( ! has_lower_limit( h ) ) {
      add_q( t , h , s );
      add( v_commitment[ t ] , - mn );
      push_row( ReserveBound_Const , std::move( vars ) , 0 , INF );
      }
     else {
      add_q( t , h , s );
      add( v_commitment[ t ] , - mn );
      add( v_start_up[ t ] , - ( hat( false , true , t , h ) - mn ) );
      push_row( ReserveBound_Const , std::move( vars ) , 0 , INF );
      if( ! last ) {
       add_q( t , h , s );
       add( v_commitment[ t ] , - mn );
       add( v_shut_down[ t + 1 ] , - ( hat( false , false , t , h ) - mn ) );
       push_row( ReserveBound_Const , std::move( vars ) , 0 , INF );
       }
      }
     }
    }
 add_rows( ReserveBound_Const , "Reserve_Bound_Const_Conversion" );

 // (12), the operating rows of the deployed powers - - - - - - - - - - - - -

 for( Index t = 0 ; t < T ; ++t )
  for( Index m = 0 ; m < f_M ; ++m ) {
   std::vector< double > ad( f_G , 0 );
   for( auto g : R )
    for( Index h = 0 ; h < f_G ; ++h )
     ad[ g ] += v_OperatingMatrix[ m ][ h ] * D( h , g );
   const double lo = get_operating_lhs( t , m );
   const double hi = get_operating_rhs( t , m );

   if( lo == hi ) {
    if( pos_entry( ad , 1 ) || pos_entry( ad , -1 ) ) {
     for( auto g : R )
      if( ad[ g ] != 0 )
       add_r( t , g , ad[ g ] );
     push_row( ReserveOperating_Const , std::move( vars ) , 0 , 0 );
     }
    continue;
    }

   for( double s : { 1.0 , -1.0 } ) {
    auto row = [ & ]( double b , double lhs , double rhs ) {
     add_ap( t , m );
     for( auto g : R )
      if( ad[ g ] != 0 )
       add_r( t , g , s * ad[ g ] );
     add( v_commitment[ t ] , - b );
     push_row( ReserveOperating_Const , std::move( vars ) , lhs , rhs );
     };
    if( ( hi < INF ) && pos_entry( ad , s ) )
     row( hi , - INF , 0 );
    if( ( lo > - INF ) && pos_entry( ad , - s ) )
     row( lo , 0 , INF );
    }
   }
 add_rows( ReserveOperating_Const , "Reserve_Operating_Const_Conversion" );

 // (13), the ramps of the deployed powers- - - - - - - - - - - - - - - - - -

 for( Index t = 0 ; t < T ; ++t )
  for( Index k = 0 ; k < f_K ; ++k ) {
   std::vector< double > cd( f_G , 0 );
   for( auto g : R )
    for( Index h = 0 ; h < f_G ; ++h )
     cd[ g ] += get_ramp_coefficient( k , h ) * D( h , g );
   for( double s : { 1.0 , -1.0 } ) {
    if( ( ! absent( v_DeltaRampUp ) ) && pos_entry( cd , s ) )
     ramp_up( ReserveRamp_Const , t , k , s );
    if( ( ! absent( v_DeltaRampDown ) ) && pos_entry( cd , - s ) )
     ramp_dn( ReserveRamp_Const , t , k , - s );
    }
   }
 add_rows( ReserveRamp_Const , "Reserve_Ramp_Const_Conversion" );

 }  // end( ConversionUnitBlock::build_rows )

/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::update_rows( ModParam issueAMod )
{
 if( ! constraints_generated() )
  return;

 // write the rows anew, comparing them with those there are: nothing is
 // changed until all of them are written, so that a throw changes nothing
 RowCmp cmp;
 f_row_cmp = & cmp;
 try {
  build_rows();
  for( auto * g : f_row_groups )
   if( cmp.pushed[ g ] != g->size() )
    throw( std::logic_error( "ConversionUnitBlock::update_rows: a group of "
                             + std::to_string( g->size() ) + " rows would "
                             "have " + std::to_string( cmp.pushed[ g ] ) +
                             " of them" ) );
  }
 catch( ... ) {
  f_row_cmp = nullptr;
  throw;
  }
 f_row_cmp = nullptr;

 if( cmp.rows.empty() && cmp.boxes.empty() )
  return;  // nothing changes

 // all the Modification in a single GroupModification
 auto nAM = un_ModBlock( make_par( par2mod( issueAMod ) ,
                                   open_channel( par2chnl( issueAMod ) ) ) );

 for( auto & c : cmp.rows ) {
  if( ! c.idx.empty() )
   static_cast< LinearFunction * >( c.row->get_function()
                                    )->modify_coefficients(
                     std::move( c.coef ) , std::move( c.idx ) , true , nAM );
  const bool lc = ( c.row->get_lhs() != c.lhs );
  const bool rc = ( c.row->get_rhs() != c.rhs );
  if( lc && rc && ( c.lhs == c.rhs ) )
   c.row->set_both( c.lhs , nAM );
  else
   if( lc && ( c.lhs > c.row->get_rhs() ) ) {  // the rhs first
    c.row->set_rhs( c.rhs , nAM );
    c.row->set_lhs( c.lhs , nAM );
    }
   else {
    if( lc )
     c.row->set_lhs( c.lhs , nAM );
    if( rc )
     c.row->set_rhs( c.rhs , nAM );
    }
  }

 for( auto & b : cmp.boxes ) {
  if( b.lhs > b.box->get_rhs() ) {  // the rhs first
   b.box->set_rhs( b.rhs , nAM );
   b.box->set_lhs( b.lhs , nAM );
   }
  else {
   if( b.box->get_lhs() != b.lhs )
    b.box->set_lhs( b.lhs , nAM );
   if( b.box->get_rhs() != b.rhs )
    b.box->set_rhs( b.rhs , nAM );
   }
  }

 close_channel( par2chnl( nAM ) );

 }  // end( ConversionUnitBlock::update_rows )

/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::generate_objective( Configuration * objc )
{
 if( objective_generated() )
  return;

 const Index T = f_time_horizon;
 DQuadFunction::v_coeff_triple vars;

 for( Index t = 0 ; t < T ; ++t )
  vars.push_back( std::make_tuple( & v_start_up[ t ] ,
                                   f_scale * get_start_up_cost( t ) , 0 ) );
 for( Index t = 0 ; t < T ; ++t )
  vars.push_back( std::make_tuple( & v_shut_down[ t ] ,
                                   f_scale * get_shut_down_cost( t ) , 0 ) );
 for( Index t = 0 ; t < T ; ++t )
  vars.push_back( std::make_tuple( & v_commitment[ t ] ,
                                   f_scale * get_const_term( t ) , 0 ) );
 for( Index g = 0 ; g < f_G ; ++g )
  for( Index t = 0 ; t < T ; ++t )
   vars.push_back( std::make_tuple( & v_active_power[ g ][ t ] ,
                                    f_scale * get_linear_term( t , g ) ,
                                    f_scale * get_quad_term( t , g ) ) );
 auto reserve = [ & ]( std::vector< std::vector< ColVariable > > & r ,
                       const MAdbl & cost ) {
  for( Index g = 0 ; g < r.size() ; ++g )
   for( Index t = 0 ; t < r[ g ].size() ; ++t )
    vars.push_back( std::make_tuple( & r[ g ][ t ] , absent( cost ) ? 0 :
                                     f_scale * cost[ t ][ g ] , 0 ) );
  };
 reserve( v_primary_reserve , v_PrimaryCost );
 reserve( v_secondary_reserve , v_SecondaryCost );

 objective.set_function( new DQuadFunction( std::move( vars ) ) );
 objective.set_sense( Objective::eMin );
 set_objective( & objective , eNoMod );

 set_objective_generated();
 }

/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::update_objective( const Subset & instants ,
                                            ModParam issueAMod )
{
 if( ! objective_generated() )
  return;

 const Index T = f_time_horizon;
 Subset nms;
 std::vector< double > lin , quad;
 auto put = [ & ]( Index i , double l , double q ) {
  nms.push_back( i );
  lin.push_back( l );
  quad.push_back( q );
  };

 for( auto t : instants ) {
  put( t , f_scale * get_start_up_cost( t ) , 0 );
  put( T + t , f_scale * get_shut_down_cost( t ) , 0 );
  put( 2 * T + t , f_scale * get_const_term( t ) , 0 );
  for( Index g = 0 ; g < f_G ; ++g )
   put( ( 3 + g ) * T + t , f_scale * get_linear_term( t , g ) ,
        f_scale * get_quad_term( t , g ) );
  }
 Index off = ( 3 + f_G ) * T;
 for( const auto * rv : { & v_primary_reserve , & v_secondary_reserve } ) {
  const auto & cost = ( rv == & v_primary_reserve ) ? v_PrimaryCost
                                                    : v_SecondaryCost;
  for( Index g = 0 ; g < rv->size() ; ++g ) {
   if( ( * rv )[ g ].empty() )
    continue;
   for( auto t : instants )
    put( off + t , absent( cost ) ? 0 : f_scale * cost[ t ][ g ] , 0 );
   off += T;
   }
  }

 static_cast< DQuadFunction * >( objective.get_function() )->modify_terms(
                quad.cbegin() , lin.cbegin() , std::move( nms ) , false ,
                un_ModBlock( issueAMod ) );
 }

/*--------------------------------------------------------------------------*/
/*-------------------- CHECKING THE ConversionUnitBlock --------------------*/
/*--------------------------------------------------------------------------*/

bool ConversionUnitBlock::is_feasible( bool useabstract ,
                                       Configuration * fsbc )
{
 double tol = DefaultFeasTol;
 bool rel_viol = true;
 extract_tolerance( fsbc , f_BlockConfig , tol , rel_viol );

 return( ColVariable::is_feasible( v_commitment , tol )
  && ColVariable::is_feasible( v_start_up , tol )
  && ColVariable::is_feasible( v_shut_down , tol )
  && ColVariable::is_feasible( v_active_power , tol )
  && ColVariable::is_feasible( v_primary_reserve , tol )
  && ColVariable::is_feasible( v_secondary_reserve , tol )
  && RowConstraint::is_feasible( Logical_Const , tol , rel_viol )
  && RowConstraint::is_feasible( MinUp_Const , tol , rel_viol )
  && RowConstraint::is_feasible( MinDown_Const , tol , rel_viol )
  && RowConstraint::is_feasible( MinPower_Const , tol , rel_viol )
  && RowConstraint::is_feasible( MaxPower_Const , tol , rel_viol )
  && RowConstraint::is_feasible( Operating_Const , tol , rel_viol )
  && RowConstraint::is_feasible( RampUp_Const , tol , rel_viol )
  && RowConstraint::is_feasible( RampDown_Const , tol , rel_viol )
  && RowConstraint::is_feasible( PrimaryRho_Const , tol , rel_viol )
  && RowConstraint::is_feasible( SecondaryRho_Const , tol , rel_viol )
  && RowConstraint::is_feasible( ReserveBound_Const , tol , rel_viol )
  && RowConstraint::is_feasible( ReserveOperating_Const , tol , rel_viol )
  && RowConstraint::is_feasible( ReserveRamp_Const , tol , rel_viol )
  && RowConstraint::is_feasible( Commitment_Bound , tol , rel_viol )
  && ( ( ! constraints_generated() ) || ( ! f_time_horizon ) ||
       RowConstraint::is_feasible( ShutDownZero_Bound , tol , rel_viol ) )
  && RowConstraint::is_feasible( ActivePower_Bound , tol , rel_viol ) );
 }

/*--------------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/

Solution * ConversionUnitBlock::get_Solution( Configuration * csolc ,
                                              bool emptys )
{
 auto sol = UnitBlock::get_Solution( csolc , true );
 if( ! emptys )
  sol->read( this );
 return( sol );
 }

/*--------------------------------------------------------------------------*/

UnitBlockSolution * ConversionUnitBlock::new_Solution( void ) const
{
 return( new ConversionUnitBlockSolution() );
 }

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR SAVING THE UNIT ------------------------*/
/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::serialize( netCDF::NcGroup & group ) const
{
 UnitBlock::serialize( group );

 auto TH = group.getDim( "TimeHorizon" );
 auto NG = group.addDim( "NumberGenerators" , f_G );
 netCDF::NcDim NM , NK;
 if( f_M )
  NM = group.addDim( "NumberOperatingRows" , f_M );
 if( ! absent( v_RampMatrix ) )
  NK = group.addDim( "NumberRampRows" , f_K );
 else
  NK = NG;

 auto put = [ & ]( const std::string & name , const MAdbl & data ,
                   const netCDF::NcDim & d0 , const netCDF::NcDim & d1 ) {
  if( absent( data ) )
   return;
  group.addVar( name , netCDF::NcDouble() , { d0 , d1 } ).putVar(
                                                             data.data() );
  };
 auto put_t = [ & ]( const std::string & name ,
                     const std::vector< double > & data ) {
  if( data.empty() )
   return;
  group.addVar( name , netCDF::NcDouble() , TH ).putVar( data.data() );
  };

 if( f_commitment_generator )
  ::serialize( group , "CommitmentGenerator" , netCDF::NcUint() ,
               f_commitment_generator );
 put( "MinPower" , v_MinPower , TH , NG );
 put( "MaxPower" , v_MaxPower , TH , NG );
 if( f_M ) {
  put( "OperatingMatrix" , v_OperatingMatrix , NM , NG );
  put( "OperatingLHS" , v_OperatingLHS , TH , NM );
  put( "OperatingRHS" , v_OperatingRHS , TH , NM );
  }
 put( "RampMatrix" , v_RampMatrix , NK , NG );
 put( "DeltaRampUp" , v_DeltaRampUp , TH , NK );
 put( "DeltaRampDown" , v_DeltaRampDown , TH , NK );
 put( "StartUpLimit" , v_StartUpLimit , TH , NG );
 put( "StartUpLowerLimit" , v_StartUpLowerLimit , TH , NG );
 put( "ShutDownLimit" , v_ShutDownLimit , TH , NG );
 put( "ShutDownLowerLimit" , v_ShutDownLowerLimit , TH , NG );
 put( "PrimaryRho" , v_PrimaryRho , TH , NG );
 put( "SecondaryRho" , v_SecondaryRho , TH , NG );
 put( "ReserveDirection" , v_ReserveDirection , NG , NG );
 put( "PrimarySpinningReserveCost" , v_PrimaryCost , TH , NG );
 put( "SecondarySpinningReserveCost" , v_SecondaryCost , TH , NG );
 put( "LinearTerm" , v_LinearTerm , TH , NG );
 put( "QuadTerm" , v_QuadTerm , TH , NG );
 put_t( "ConstTerm" , v_ConstTerm );
 put_t( "StartUpCost" , v_StartUpCost );
 put_t( "ShutDownCost" , v_ShutDownCost );
 put_t( "FixedConsumption" , v_FixedConsumption );
 put_t( "InertiaCommitment" , v_InertiaCommitment );

 if( std::any_of( v_InertiaPower.begin() , v_InertiaPower.end() ,
                  []( const auto & r ) { return( ! r.empty() ); } ) ) {
  MAdbl ip( boost::extents[ f_time_horizon ][ f_G ] );
  for( Index g = 0 ; g < f_G ; ++g )
   for( Index t = 0 ; t < v_InertiaPower[ g ].size() ; ++t )
    ip[ t ][ g ] = v_InertiaPower[ g ][ t ];
  put( "InertiaPower" , ip , TH , NG );
  }

 group.addVar( "InitialPower" , netCDF::NcDouble() , NG ).putVar(
                                                    v_InitialPower.data() );
 ::serialize( group , "InitUpDownTime" , netCDF::NcInt() ,
              f_InitUpDownTime );
 ::serialize( group , "MinUpTime" , netCDF::NcUint() , f_MinUpTime );
 ::serialize( group , "MinDownTime" , netCDF::NcUint() , f_MinDownTime );
 if( f_scale != 1 )
  ::serialize( group , "Scale" , netCDF::NcDouble() , f_scale );

 }  // end( ConversionUnitBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::issue_mod( int type , Index index ,
                                     Subset && subset , ModParam issuePMod )
{
 if( issue_pmod( issuePMod ) )
  Block::add_Modification( std::make_shared< ConversionUnitBlockMod >(
                            this , type , index , std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );
 }

/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::guts_of_set( MAdbl & data , double dflt , int type ,
                                       Index col , MF_dbl_it values ,
                                       Subset && subset , bool ordered ,
                                       ModParam issuePMod ,
                                       ModParam issueAMod , bool rows ,
                                       const std::string & who )
{
 if( subset.empty() )
  return;

 std::vector< double > sorted;
 if( ! ordered )
  sort_by_index( subset , values , sorted );

 if( subset.back() >= f_time_horizon )
  throw( std::invalid_argument( who + ": invalid instant in the subset" ) );

 Index ncol = 0;
 switch( type ) {
  case( ConversionUnitBlockMod::eSetOpRHS ):
  case( ConversionUnitBlockMod::eSetOpLHS ): ncol = f_M; break;
  case( ConversionUnitBlockMod::eSetRampUp ):
  case( ConversionUnitBlockMod::eSetRampDown ): ncol = f_K; break;
  default: ncol = f_G;
  }
 if( col >= ncol )
  throw( std::invalid_argument( who + ": invalid index " +
                                std::to_string( col ) ) );

 if( absent( data ) ) {
  if( std::all_of( values , values + subset.size() ,
                   [ dflt ]( double v ) { return( v == dflt ); } ) )
   return;
  if( ( type == ConversionUnitBlockMod::eSetRampUp ) ||
      ( type == ConversionUnitBlockMod::eSetRampDown ) )
   throw( std::invalid_argument( who + ": the unit has no such ramp" ) );
  }
 else {
  bool same = true;
  auto it = values;
  for( auto t : subset )
   same &= ( data[ t ][ col ] == *( it++ ) );
  if( same )
   return;
  }

 if( ( ( type == ConversionUnitBlockMod::eSetPrCost ) &&
       ( v_primary_reserve.size() <= col ||
         v_primary_reserve[ col ].empty() ) && objective_generated() ) ||
     ( ( type == ConversionUnitBlockMod::eSetScCost ) &&
       ( v_secondary_reserve.size() <= col ||
         v_secondary_reserve[ col ].empty() ) && objective_generated() ) )
  throw( std::invalid_argument( who + ": generator " +
                                std::to_string( col ) + " has no reserve "
                                "Variable" ) );

 if( not_dry_run( issuePMod ) ) {
  MAdbl old;
  copy_into( old , data );
  if( absent( data ) ) {
   data.resize( boost::extents[ f_time_horizon ][ ncol ] );
   std::fill( data.data() , data.data() + data.num_elements() , dflt );
   }
  auto it = values;
  for( auto t : subset )
   data[ t ][ col ] = *( it++ );

  try {
   check_data( who );
   if( rows && not_dry_run( issueAMod ) )
    update_rows( issueAMod );
   }
  catch( ... ) {
   copy_into( data , old );
   throw;
   }

  if( ( ! rows ) && not_dry_run( issueAMod ) )
   update_objective( subset , issueAMod );
  }

 issue_mod( type , col , std::move( subset ) , issuePMod );
 }

/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::guts_of_set( std::vector< double > & data ,
                                       int type , MF_dbl_it values ,
                                       Subset && subset , bool ordered ,
                                       ModParam issuePMod ,
                                       ModParam issueAMod ,
                                       const std::string & who )
{
 if( subset.empty() )
  return;

 std::vector< double > sorted;
 if( ! ordered )
  sort_by_index( subset , values , sorted );

 if( subset.back() >= f_time_horizon )
  throw( std::invalid_argument( who + ": invalid instant in the subset" ) );

 if( data.empty() ) {
  if( std::all_of( values , values + subset.size() ,
                   []( double v ) { return( v == 0 ); } ) )
   return;
  }
 else {
  bool same = true;
  auto it = values;
  for( auto t : subset )
   same &= ( data[ t ] == *( it++ ) );
  if( same )
   return;
  }

 if( not_dry_run( issuePMod ) ) {
  if( data.empty() )
   data.assign( f_time_horizon , 0 );
  auto it = values;
  for( auto t : subset )
   data[ t ] = *( it++ );
  if( not_dry_run( issueAMod ) )
   update_objective( subset , issueAMod );
  }

 issue_mod( type , 0 , std::move( subset ) , issuePMod );
 }

/*--------------------------------------------------------------------------*/
// the subset of a range of instants, cut at the time horizon

static Subset range_subset( Block::Range rng , Index T )
{
 rng.second = std::min( rng.second , T );
 Subset s;
 for( Index t = rng.first ; t < rng.second ; ++t )
  s.push_back( t );
 return( s );
 }

/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::set_linear_term( MF_dbl_it values , Index g ,
                                           Subset && subset ,
                                           const bool ordered ,
                                           ModParam issuePMod ,
                                           ModParam issueAMod )
{
 guts_of_set( v_LinearTerm , 0 , ConversionUnitBlockMod::eSetLinT , g ,
              values , std::move( subset ) , ordered , issuePMod ,
              issueAMod , false , "ConversionUnitBlock::set_linear_term" );
 }

void ConversionUnitBlock::set_linear_term( MF_dbl_it values , Index g ,
                                           Range rng , ModParam issuePMod ,
                                           ModParam issueAMod )
{
 set_linear_term( values , g , range_subset( rng , f_time_horizon ) , true ,
                  issuePMod , issueAMod );
 }

/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::set_quad_term( MF_dbl_it values , Index g ,
                                         Subset && subset ,
                                         const bool ordered ,
                                         ModParam issuePMod ,
                                         ModParam issueAMod )
{
 guts_of_set( v_QuadTerm , 0 , ConversionUnitBlockMod::eSetQuadT , g ,
              values , std::move( subset ) , ordered , issuePMod ,
              issueAMod , false , "ConversionUnitBlock::set_quad_term" );
 }

void ConversionUnitBlock::set_quad_term( MF_dbl_it values , Index g ,
                                         Range rng , ModParam issuePMod ,
                                         ModParam issueAMod )
{
 set_quad_term( values , g , range_subset( rng , f_time_horizon ) , true ,
                issuePMod , issueAMod );
 }

/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::set_primary_spinning_reserve_cost(
                         MF_dbl_it values , Index g , Range rng ,
                         ModParam issuePMod , ModParam issueAMod )
{
 guts_of_set( v_PrimaryCost , 0 , ConversionUnitBlockMod::eSetPrCost , g ,
              values , range_subset( rng , f_time_horizon ) , true ,
              issuePMod , issueAMod , false ,
              "ConversionUnitBlock::set_primary_spinning_reserve_cost" );
 }

/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::set_secondary_spinning_reserve_cost(
                         MF_dbl_it values , Index g , Range rng ,
                         ModParam issuePMod , ModParam issueAMod )
{
 guts_of_set( v_SecondaryCost , 0 , ConversionUnitBlockMod::eSetScCost , g ,
              values , range_subset( rng , f_time_horizon ) , true ,
              issuePMod , issueAMod , false ,
              "ConversionUnitBlock::set_secondary_spinning_reserve_cost" );
 }

/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::set_const_term( MF_dbl_it values ,
                                          Subset && subset ,
                                          const bool ordered ,
                                          ModParam issuePMod ,
                                          ModParam issueAMod )
{
 guts_of_set( v_ConstTerm , ConversionUnitBlockMod::eSetConstT , values ,
              std::move( subset ) , ordered , issuePMod , issueAMod ,
              "ConversionUnitBlock::set_const_term" );
 }

void ConversionUnitBlock::set_const_term( MF_dbl_it values , Range rng ,
                                          ModParam issuePMod ,
                                          ModParam issueAMod )
{
 set_const_term( values , range_subset( rng , f_time_horizon ) , true ,
                 issuePMod , issueAMod );
 }

/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::set_startup_costs( MF_dbl_it values ,
                                             Subset && subset ,
                                             const bool ordered ,
                                             ModParam issuePMod ,
                                             ModParam issueAMod )
{
 guts_of_set( v_StartUpCost , ConversionUnitBlockMod::eSetSUC , values ,
              std::move( subset ) , ordered , issuePMod , issueAMod ,
              "ConversionUnitBlock::set_startup_costs" );
 }

void ConversionUnitBlock::set_startup_costs( MF_dbl_it values , Range rng ,
                                             ModParam issuePMod ,
                                             ModParam issueAMod )
{
 set_startup_costs( values , range_subset( rng , f_time_horizon ) , true ,
                    issuePMod , issueAMod );
 }

/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::set_shutdown_costs( MF_dbl_it values ,
                                              Subset && subset ,
                                              const bool ordered ,
                                              ModParam issuePMod ,
                                              ModParam issueAMod )
{
 guts_of_set( v_ShutDownCost , ConversionUnitBlockMod::eSetSDC , values ,
              std::move( subset ) , ordered , issuePMod , issueAMod ,
              "ConversionUnitBlock::set_shutdown_costs" );
 }

void ConversionUnitBlock::set_shutdown_costs( MF_dbl_it values , Range rng ,
                                              ModParam issuePMod ,
                                              ModParam issueAMod )
{
 set_shutdown_costs( values , range_subset( rng , f_time_horizon ) , true ,
                     issuePMod , issueAMod );
 }

/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::set_maximum_power( MF_dbl_it values , Index g ,
                                             Subset && subset ,
                                             const bool ordered ,
                                             ModParam issuePMod ,
                                             ModParam issueAMod )
{
 guts_of_set( v_MaxPower , 0 , ConversionUnitBlockMod::eSetMaxP , g ,
              values , std::move( subset ) , ordered , issuePMod ,
              issueAMod , true , "ConversionUnitBlock::set_maximum_power" );
 }

void ConversionUnitBlock::set_maximum_power( MF_dbl_it values , Index g ,
                                             Range rng , ModParam issuePMod ,
                                             ModParam issueAMod )
{
 set_maximum_power( values , g , range_subset( rng , f_time_horizon ) ,
                    true , issuePMod , issueAMod );
 }

/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::set_minimum_power( MF_dbl_it values , Index g ,
                                             Subset && subset ,
                                             const bool ordered ,
                                             ModParam issuePMod ,
                                             ModParam issueAMod )
{
 guts_of_set( v_MinPower , 0 , ConversionUnitBlockMod::eSetMinP , g ,
              values , std::move( subset ) , ordered , issuePMod ,
              issueAMod , true , "ConversionUnitBlock::set_minimum_power" );
 }

void ConversionUnitBlock::set_minimum_power( MF_dbl_it values , Index g ,
                                             Range rng , ModParam issuePMod ,
                                             ModParam issueAMod )
{
 set_minimum_power( values , g , range_subset( rng , f_time_horizon ) ,
                    true , issuePMod , issueAMod );
 }

/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::set_operating_rhs( MF_dbl_it values , Index m ,
                                             Range rng , ModParam issuePMod ,
                                             ModParam issueAMod )
{
 guts_of_set( v_OperatingRHS , INF , ConversionUnitBlockMod::eSetOpRHS , m ,
              values , range_subset( rng , f_time_horizon ) , true ,
              issuePMod , issueAMod , true ,
              "ConversionUnitBlock::set_operating_rhs" );
 }

/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::set_operating_lhs( MF_dbl_it values , Index m ,
                                             Range rng , ModParam issuePMod ,
                                             ModParam issueAMod )
{
 guts_of_set( v_OperatingLHS , - INF , ConversionUnitBlockMod::eSetOpLHS ,
              m , values , range_subset( rng , f_time_horizon ) , true ,
              issuePMod , issueAMod , true ,
              "ConversionUnitBlock::set_operating_lhs" );
 }

/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::set_delta_ramp_up( MF_dbl_it values , Index k ,
                                             Range rng , ModParam issuePMod ,
                                             ModParam issueAMod )
{
 guts_of_set( v_DeltaRampUp , 0 , ConversionUnitBlockMod::eSetRampUp , k ,
              values , range_subset( rng , f_time_horizon ) , true ,
              issuePMod , issueAMod , true ,
              "ConversionUnitBlock::set_delta_ramp_up" );
 }

/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::set_delta_ramp_down( MF_dbl_it values , Index k ,
                                               Range rng ,
                                               ModParam issuePMod ,
                                               ModParam issueAMod )
{
 guts_of_set( v_DeltaRampDown , 0 , ConversionUnitBlockMod::eSetRampDown ,
              k , values , range_subset( rng , f_time_horizon ) , true ,
              issuePMod , issueAMod , true ,
              "ConversionUnitBlock::set_delta_ramp_down" );
 }

/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::set_initial_power( MF_dbl_it values , Range rng ,
                                             ModParam issuePMod ,
                                             ModParam issueAMod )
{
 static const std::string who = "ConversionUnitBlock::set_initial_power";
 rng.second = std::min( rng.second , f_G );
 if( rng.second <= rng.first )
  return;
 if( std::equal( values , values + ( rng.second - rng.first ) ,
                 v_InitialPower.begin() + rng.first ) )
  return;

 if( not_dry_run( issuePMod ) ) {
  const auto old = v_InitialPower;
  std::copy( values , values + ( rng.second - rng.first ) ,
             v_InitialPower.begin() + rng.first );
  try {
   check_data( who );
   if( not_dry_run( issueAMod ) )
    update_rows( issueAMod );
   }
  catch( ... ) {
   v_InitialPower = old;
   throw;
   }
  }

 Subset gens;
 for( Index g = rng.first ; g < rng.second ; ++g )
  gens.push_back( g );
 issue_mod( ConversionUnitBlockMod::eSetInitP , 0 , std::move( gens ) ,
            issuePMod );
 }

/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::set_init_updown_time( MF_int_it values ,
                                                Range rng ,
                                                ModParam issuePMod ,
                                                ModParam issueAMod )
{
 static const std::string who = "ConversionUnitBlock::set_init_updown_time";
 if( rng.second <= rng.first )
  return;
 if( *values == f_InitUpDownTime )
  return;

 if( not_dry_run( issuePMod ) ) {
  const auto old = f_InitUpDownTime;
  f_InitUpDownTime = *values;
  try {
   check_data( who );
   if( not_dry_run( issueAMod ) )
    update_rows( issueAMod );
   }
  catch( ... ) {
   f_InitUpDownTime = old;
   throw;
   }
  }

 issue_mod( ConversionUnitBlockMod::eSetInitUD , 0 , Subset() ,
            issuePMod );
 }

/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::scale( MF_dbl_it values , Subset && subset ,
                                 const bool ordered ,
                                 c_ModParam issuePMod ,
                                 c_ModParam issueAMod )
{
 if( subset.empty() || ( f_scale == *values ) )
  return;

 if( not_dry_run( issuePMod ) ) {
  f_scale = *values;
  if( not_dry_run( issueAMod ) && objective_generated() )
   update_objective( range_subset( Range( 0 , f_time_horizon ) ,
                                   f_time_horizon ) , issueAMod );
  }

 if( issue_pmod( issuePMod ) )
  Block::add_Modification( std::make_shared< UnitBlockMod >(
                            this , UnitBlockMod::eScale ) ,
                           Observer::par2chnl( issuePMod ) );
 else
  if( auto fb = get_f_Block() )
   fb->add_Modification( std::make_shared< UnitBlockMod >(
                          this , UnitBlockMod::eScale ) ,
                         Observer::par2chnl( issuePMod ) );
 }

/*--------------------------------------------------------------------------*/
/*------------- THE SETTERS OVER THE ENTRIES [ g ][ t ] --------------------*/
/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::flat_set( MAdbl & data , double dflt , int type ,
                                    Index ncol , MF_dbl_it values ,
                                    Subset && subset , bool ordered ,
                                    ModParam issuePMod , ModParam issueAMod ,
                                    bool rows , const std::string & who )
{
 const Index T = f_time_horizon;
 // the instants and the values of each column, in the order given
 std::map< Index , std::pair< Subset , std::vector< double > > > cols;
 for( auto i : subset ) {
  if( ( ! T ) || ( i >= ncol * T ) )
   throw( std::invalid_argument( who + ": invalid index " +
                                 std::to_string( i ) + " in the subset" ) );
  auto & c = cols[ i / T ];
  c.first.push_back( i % T );
  c.second.push_back( *( values++ ) );
  }
 for( auto & [ col , c ] : cols )
  guts_of_set( data , dflt , type , col , c.second.cbegin() ,
               std::move( c.first ) , ordered , issuePMod , issueAMod ,
               rows , who );
 }

/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::set_linear_term(
 MF_dbl_it values , Subset && subset , const bool ordered ,
 ModParam issuePMod , ModParam issueAMod )
{
 flat_set( v_LinearTerm , 0 , ConversionUnitBlockMod::eSetLinT , f_G ,
           values , std::move( subset ) , ordered , issuePMod , issueAMod ,
           false , "ConversionUnitBlock::set_linear_term" );
 }

void ConversionUnitBlock::set_linear_term(
 MF_dbl_it values , Range rng , ModParam issuePMod , ModParam issueAMod )
{
 set_linear_term( values , range_subset( rng , f_G * f_time_horizon ) ,
  true , issuePMod , issueAMod );
 }
/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::set_quad_term(
 MF_dbl_it values , Subset && subset , const bool ordered ,
 ModParam issuePMod , ModParam issueAMod )
{
 flat_set( v_QuadTerm , 0 , ConversionUnitBlockMod::eSetQuadT , f_G ,
           values , std::move( subset ) , ordered , issuePMod , issueAMod ,
           false , "ConversionUnitBlock::set_quad_term" );
 }

void ConversionUnitBlock::set_quad_term(
 MF_dbl_it values , Range rng , ModParam issuePMod , ModParam issueAMod )
{
 set_quad_term( values , range_subset( rng , f_G * f_time_horizon ) ,
  true , issuePMod , issueAMod );
 }
/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::set_primary_spinning_reserve_cost(
 MF_dbl_it values , Subset && subset , const bool ordered ,
 ModParam issuePMod , ModParam issueAMod )
{
 flat_set( v_PrimaryCost , 0 , ConversionUnitBlockMod::eSetPrCost , f_G ,
           values , std::move( subset ) , ordered , issuePMod , issueAMod ,
           false , "ConversionUnitBlock::set_primary_spinning_reserve_cost" );
 }

void ConversionUnitBlock::set_primary_spinning_reserve_cost(
 MF_dbl_it values , Range rng , ModParam issuePMod , ModParam issueAMod )
{
 set_primary_spinning_reserve_cost( values ,
  range_subset( rng , f_G * f_time_horizon ) , true , issuePMod ,
  issueAMod );
 }
/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::set_secondary_spinning_reserve_cost(
 MF_dbl_it values , Subset && subset , const bool ordered ,
 ModParam issuePMod , ModParam issueAMod )
{
 flat_set( v_SecondaryCost , 0 , ConversionUnitBlockMod::eSetScCost , f_G ,
           values , std::move( subset ) , ordered , issuePMod , issueAMod ,
           false ,
           "ConversionUnitBlock::set_secondary_spinning_reserve_cost" );
 }

void ConversionUnitBlock::set_secondary_spinning_reserve_cost(
 MF_dbl_it values , Range rng , ModParam issuePMod , ModParam issueAMod )
{
 set_secondary_spinning_reserve_cost( values ,
  range_subset( rng , f_G * f_time_horizon ) , true , issuePMod ,
  issueAMod );
 }
/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::set_maximum_power(
 MF_dbl_it values , Subset && subset , const bool ordered ,
 ModParam issuePMod , ModParam issueAMod )
{
 flat_set( v_MaxPower , 0 , ConversionUnitBlockMod::eSetMaxP , f_G ,
           values , std::move( subset ) , ordered , issuePMod , issueAMod ,
           true , "ConversionUnitBlock::set_maximum_power" );
 }

void ConversionUnitBlock::set_maximum_power(
 MF_dbl_it values , Range rng , ModParam issuePMod , ModParam issueAMod )
{
 set_maximum_power( values , range_subset( rng , f_G * f_time_horizon ) ,
  true , issuePMod , issueAMod );
 }
/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::set_minimum_power(
 MF_dbl_it values , Subset && subset , const bool ordered ,
 ModParam issuePMod , ModParam issueAMod )
{
 flat_set( v_MinPower , 0 , ConversionUnitBlockMod::eSetMinP , f_G ,
           values , std::move( subset ) , ordered , issuePMod , issueAMod ,
           true , "ConversionUnitBlock::set_minimum_power" );
 }

void ConversionUnitBlock::set_minimum_power(
 MF_dbl_it values , Range rng , ModParam issuePMod , ModParam issueAMod )
{
 set_minimum_power( values , range_subset( rng , f_G * f_time_horizon ) ,
  true , issuePMod , issueAMod );
 }
/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::set_operating_rhs(
 MF_dbl_it values , Subset && subset , const bool ordered ,
 ModParam issuePMod , ModParam issueAMod )
{
 flat_set( v_OperatingRHS , INF , ConversionUnitBlockMod::eSetOpRHS , f_M ,
           values , std::move( subset ) , ordered , issuePMod , issueAMod ,
           true , "ConversionUnitBlock::set_operating_rhs" );
 }

void ConversionUnitBlock::set_operating_rhs(
 MF_dbl_it values , Range rng , ModParam issuePMod , ModParam issueAMod )
{
 set_operating_rhs( values , range_subset( rng , f_M * f_time_horizon ) ,
  true , issuePMod , issueAMod );
 }
/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::set_operating_lhs(
 MF_dbl_it values , Subset && subset , const bool ordered ,
 ModParam issuePMod , ModParam issueAMod )
{
 flat_set( v_OperatingLHS , - INF , ConversionUnitBlockMod::eSetOpLHS , f_M ,
           values , std::move( subset ) , ordered , issuePMod , issueAMod ,
           true , "ConversionUnitBlock::set_operating_lhs" );
 }

void ConversionUnitBlock::set_operating_lhs(
 MF_dbl_it values , Range rng , ModParam issuePMod , ModParam issueAMod )
{
 set_operating_lhs( values , range_subset( rng , f_M * f_time_horizon ) ,
  true , issuePMod , issueAMod );
 }
/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::set_delta_ramp_up(
 MF_dbl_it values , Subset && subset , const bool ordered ,
 ModParam issuePMod , ModParam issueAMod )
{
 flat_set( v_DeltaRampUp , 0 , ConversionUnitBlockMod::eSetRampUp , f_K ,
           values , std::move( subset ) , ordered , issuePMod , issueAMod ,
           true , "ConversionUnitBlock::set_delta_ramp_up" );
 }

void ConversionUnitBlock::set_delta_ramp_up(
 MF_dbl_it values , Range rng , ModParam issuePMod , ModParam issueAMod )
{
 set_delta_ramp_up( values , range_subset( rng , f_K * f_time_horizon ) ,
  true , issuePMod , issueAMod );
 }
/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::set_delta_ramp_down(
 MF_dbl_it values , Subset && subset , const bool ordered ,
 ModParam issuePMod , ModParam issueAMod )
{
 flat_set( v_DeltaRampDown , 0 , ConversionUnitBlockMod::eSetRampDown , f_K ,
           values , std::move( subset ) , ordered , issuePMod , issueAMod ,
           true , "ConversionUnitBlock::set_delta_ramp_down" );
 }

void ConversionUnitBlock::set_delta_ramp_down(
 MF_dbl_it values , Range rng , ModParam issuePMod , ModParam issueAMod )
{
 set_delta_ramp_down( values , range_subset( rng , f_K * f_time_horizon ) ,
  true , issuePMod , issueAMod );
 }
/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::set_initial_power( MF_dbl_it values ,
                                             Subset && subset ,
                                             const bool ordered ,
                                             ModParam issuePMod ,
                                             ModParam issueAMod )
{
 // one change per generator, so that a refused one leaves the others
 for( auto g : subset ) {
  if( g >= f_G )
   throw( std::invalid_argument( "ConversionUnitBlock::set_initial_power: "
                                 "invalid generator " +
                                 std::to_string( g ) ) );
  set_initial_power( values++ , Range( g , g + 1 ) , issuePMod ,
                     issueAMod );
  }
 }

/*--------------------------------------------------------------------------*/

void ConversionUnitBlock::set_init_updown_time( MF_int_it values ,
                                                Subset && subset ,
                                                const bool ordered ,
                                                ModParam issuePMod ,
                                                ModParam issueAMod )
{
 if( ! subset.empty() )
  set_init_updown_time( values , Range( 0 , 1 ) , issuePMod , issueAMod );
 }

/*--------------------------------------------------------------------------*/
/*---------------- METHODS OF ConversionUnitBlockSolution ------------------*/
/*--------------------------------------------------------------------------*/

void ConversionUnitBlockSolution::deserialize( const netCDF::NcGroup & group )
{
 UnitBlockSolution::deserialize( group );
 ::deserialize( group , "ConversionStartUp" , f_time_horizon , v_start_up ,
                true , false );
 ::deserialize( group , "ConversionShutDown" , f_time_horizon , v_shut_down ,
                true , false );
 }

/*--------------------------------------------------------------------------*/

void ConversionUnitBlockSolution::read( const Block * block )
{
 auto CUB = dynamic_cast< const ConversionUnitBlock * >( block );
 if( ! CUB )
  throw( std::invalid_argument( "ConversionUnitBlockSolution::read: block "
                                "is not a ConversionUnitBlock" ) );

 UnitBlockSolution::read( CUB );

 auto cub = const_cast< ConversionUnitBlock * >( CUB );
 const Index T = CUB->get_time_horizon();
 if( ( ! get_commitment().empty() ) && cub->get_start_up() ) {
  v_start_up.resize( T );
  v_shut_down.resize( T );
  for( Index t = 0 ; t < T ; ++t ) {
   v_start_up[ t ] = cub->get_start_up()[ t ].get_value();
   v_shut_down[ t ] = cub->get_shut_down()[ t ].get_value();
   }
  }
 else {
  v_start_up.clear();
  v_shut_down.clear();
  }
 }

/*--------------------------------------------------------------------------*/

void ConversionUnitBlockSolution::write( Block * block )
{
 auto CUB = dynamic_cast< ConversionUnitBlock * >( block );
 if( ! CUB )
  throw( std::invalid_argument( "ConversionUnitBlockSolution::write: block "
                                "is not a ConversionUnitBlock" ) );

 UnitBlockSolution::write( CUB );

 const Index T = CUB->get_time_horizon();
 if( v_start_up.size() == T )
  if( auto su = CUB->get_start_up() )
   for( Index t = 0 ; t < T ; ++t )
    su[ t ].set_value( v_start_up[ t ] );
 if( v_shut_down.size() == T )
  if( auto sd = CUB->get_shut_down() )
   for( Index t = 0 ; t < T ; ++t )
    sd[ t ].set_value( v_shut_down[ t ] );
 }

/*--------------------------------------------------------------------------*/

void ConversionUnitBlockSolution::serialize( netCDF::NcGroup & group ) const
{
 UnitBlockSolution::serialize( group );
 auto TH = group.getDim( "TimeHorizon" );
 if( TH.isNull() )
  TH = group.addDim( "TimeHorizon" , f_time_horizon );
 if( ! v_start_up.empty() )
  group.addVar( "ConversionStartUp" , netCDF::NcDouble() , TH ).putVar(
                                                        v_start_up.data() );
 if( ! v_shut_down.empty() )
  group.addVar( "ConversionShutDown" , netCDF::NcDouble() , TH ).putVar(
                                                       v_shut_down.data() );
 }

/*--------------------------------------------------------------------------*/

ConversionUnitBlockSolution * ConversionUnitBlockSolution::scale(
                                                       double factor ) const
{
 auto sol = clone();
 if( factor == 1 )
  return( sol );
 guts_of_scale( sol , factor );
 for( auto & v : sol->v_start_up )
  v *= factor;
 for( auto & v : sol->v_shut_down )
  v *= factor;
 return( sol );
 }

/*--------------------------------------------------------------------------*/

void ConversionUnitBlockSolution::sum( const Solution * solution ,
                                       double multiplier )
{
 UnitBlockSolution::sum( solution , multiplier );

 auto CS = dynamic_cast< const ConversionUnitBlockSolution * >( solution );
 if( ! CS )
  throw( std::invalid_argument( "ConversionUnitBlockSolution::sum: "
                                "solution not a "
                                "ConversionUnitBlockSolution" ) );
 if( ( v_start_up.size() != CS->v_start_up.size() ) ||
     ( v_shut_down.size() != CS->v_shut_down.size() ) )
  throw( std::invalid_argument( "ConversionUnitBlockSolution::sum: "
                                "inconsistent start-up or shut-down "
                                "indicators" ) );
 for( Index i = 0 ; i < v_start_up.size() ; ++i )
  v_start_up[ i ] += CS->v_start_up[ i ] * multiplier;
 for( Index i = 0 ; i < v_shut_down.size() ; ++i )
  v_shut_down[ i ] += CS->v_shut_down[ i ] * multiplier;
 }

/*--------------------------------------------------------------------------*/

ConversionUnitBlockSolution * ConversionUnitBlockSolution::clone(
                                                         bool empty ) const
{
 auto sol = new ConversionUnitBlockSolution();
 if( ! empty ) {
  guts_of_clone( sol );
  sol->v_start_up = v_start_up;
  sol->v_shut_down = v_shut_down;
  }
 return( sol );
 }

/*--------------------------------------------------------------------------*/
/*------------------- End File ConversionUnitBlock.cpp ---------------------*/
/*--------------------------------------------------------------------------*/
