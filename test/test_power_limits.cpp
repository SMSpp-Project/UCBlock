/*--------------------------------------------------------------------------*/
/*------------------------ File test_power_limits.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The change of the data of a ThermalUnitBlock (or NuclearUnitBlock) that
 * enter the power limits, i.e., MaxPower and Availability, and of the
 * initial state, i.e., InitialPower and InitUpDownTime, after the abstract
 * representation is generated, compared with the same unit read afresh from
 * the changed data.
 *
 * For every unit of a small set of instances (with and without ramps,
 * with the reserves, FixToMaximum, the design variable, given or default
 * StartUpLimit and ShutDownLimit, the initial state fixing the first
 * instants, and the initial power at MinPower[ 0 ], which decides whether
 * the shut-down interval at 0 exists), for every formulation (3bin, T, pt,
 * DP, SU, SD, SUSD, with and without the Perspective Cuts, and for a
 * NuclearUnitBlock the T and 3bin formulations with the bits TightRules,
 * TightRamp and TightCuts), and for every change of a list (one instant
 * and a range, through the Subset and the Range versions of the setters,
 * up and down, availability 0 and 1, a sequence of three changes, values
 * that have to be refused), the test
 *
 * - reads the unit A from the data, generates the formulation, attaches the
 *   Solver of a BlockSolverConfig read from file (a :MILPSolver first, which
 *   writes the model it solves in the file of its strOutputFile, and the
 *   dynamic programming Solver of the unit), and solves once, so that every
 *   Solver has loaded the model and the change reaches it as Modification;
 *
 * - applies the change to A and solves again;
 *
 * - reads the unit B from the changed data (the data up to the change that
 *   A refused, if it refused one), generates the same formulation and
 *   solves it;
 *
 * - compares the two models written by the :MILPSolver in canonical form
 *   [see lp_compare.h], the differences being reported per group of rows,
 *   the optimal values of the :MILPSolver, the value of every dynamic
 *   programming Solver with that of the :MILPSolver (on A and on B), and
 *   the data of A and B that the change concerns (MaxPower, Availability,
 *   the operational bounds, StartUpLimit, ShutDownLimit, InitialPower,
 *   InitUpDownTime).
 *
 * A step of a change of MaxPower or Availability that gives valid data (B
 * can be read and generated) has to be applied, and A has to equal B; a
 * step that gives invalid data has to be refused, leaving A as the steps
 * before it left it; a change of InitialPower may be refused, and then it
 * has to leave A as it was, and so may a change of InitUpDownTime that
 * changes whether the unit is on before the horizon or the first instant
 * at which it may switch (which decide the Variable), while one that
 * changes neither has to be applied. Without a :MILPSolver of
 * the BlockSolverConfig in the factory the test is skipped.
 *
 * The configuration files, read from the directory given to
 * test_power_limits(), are the BlockConfig of each formulation,
 * TUBCfg-<formulation>.txt for both kinds of unit and
 * NUBCfg-<formulation>.txt for those with the bits of a NuclearUnitBlock,
 * and the BlockSolverConfig of the Solver, LPCmpBSCfg.txt (a
 * ThermalUnitBlock), LPCmpBSCfg-nuc.txt (a NuclearUnitBlock) and
 * LPCmpMILPCfg.txt (the ComputeConfig of the :MILPSolver, which both
 * include).
 *
 * The whole matrix of units, formulations and changes is the full test;
 * the reduced one, which the unit test runs by default, has all the units
 * and all the changes in the formulations 3bin, T and SUSD (with and
 * without the Perspective Cuts) and DP, and in every formulation of a
 * NuclearUnitBlock [see test_power_limits()]. */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cctype>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <netcdf>

#if defined( __unix__ ) || defined( __APPLE__ )
#include <sys/wait.h>
#include <unistd.h>
#endif

#include "lp_compare.h"

#include "BlockSolverConfig.h"
#include "NuclearUnitBlock.h"
#include "ThermalUnitBlock.h"

/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

namespace
{

using Index = Block::Index;
using Range = Block::Range;
using Subset = Block::Subset;

/*--------------------------------------------------------------------------*/
/*------------------------------ GLOBALS -----------------------------------*/
/*--------------------------------------------------------------------------*/

int n_fail = 0;            // failed cases
bool verbose = false;      // one line per case, also when it passes
std::string only;          // only the cases whose "unit | form | change"
                           // contains it, if not empty
std::string out_dir;       // where the models of the failed cases go
int n_kept = 0;            // models kept so far
int n_case = 0;            // cases so far, which name the models
const int max_kept = 40;   // at most this many failed cases keep the models
std::string keep_dir;      // if not empty, where the models of every case
                           // go, named after the case (unit, formulation,
                           // change), so that two builds can be compared

/*--------------------------------------------------------------------------*/
/*------------------------------ THE DATA ----------------------------------*/
/*--------------------------------------------------------------------------*/
/* The data of a unit as written in its netCDF group: the vectors over the
 * instants ("NumberIntervals" = T), the scalars, and for a NuclearUnitBlock
 * the two PowerBands. A datum that is not in the maps is not written. */

struct UnitData {
 std::string name;                                   // for the report
 bool nuclear = false;
 bool reserves = false;                              // set_reserve_vars( 3 )
 Index T = 6;
 std::map< std::string , std::vector< double > > vec;
 std::map< std::string , double > dbl;               // NcDouble scalars
 std::map< std::string , int > num;                  // NcInt scalars
 std::vector< double > bands;                        // PowerBands

 double at( const std::string & n , Index t , double dflt = 0 ) const {
  const auto it = vec.find( n );
  return( it == vec.end() ? dflt : it->second[ t ] );
  }

 double quad( void ) const {
  const auto it = vec.find( "QuadTerm" );
  if( it == vec.end() )
   return( 0 );
  double s = 0;
  for( auto v : it->second )
   s += std::abs( v );
  return( s );
  }
 };

/*--------------------------------------------------------------------------*/
/// the unit read from the data d (a dataset in memory, nothing on disk)

ThermalUnitBlock * new_unit( const UnitData & d )
{
 int ncid;
 if( nc_create( "power_limits_test.nc4" ,
                NC_DISKLESS | NC_NETCDF4 | NC_CLOBBER , & ncid ) != NC_NOERR )
  throw( std::runtime_error( "new_unit: cannot create a netCDF dataset" ) );

 Block * b = nullptr;
 try {
  auto g = netCDF::NcGroup( ncid ).addGroup( "unit" );
  g.putAtt( "type" , d.nuclear ? "NuclearUnitBlock" : "ThermalUnitBlock" );
  g.addDim( "TimeHorizon" , d.T );
  auto NI = g.addDim( "NumberIntervals" , d.T );
  for( const auto & v : d.vec )
   g.addVar( v.first , netCDF::NcDouble() , NI ).putVar( v.second.data() );
  for( const auto & v : d.dbl )
   g.addVar( v.first , netCDF::NcDouble() ).putVar( & v.second );
  for( const auto & v : d.num )
   g.addVar( v.first , netCDF::NcInt() ).putVar( & v.second );
  if( ! d.bands.empty() ) {
   auto nb = g.addDim( "NumberPowerBands" , d.bands.size() );
   g.addVar( "PowerBands" , netCDF::NcDouble() , nb ).putVar(
                                                          d.bands.data() );
   }
  b = Block::new_Block( g );
  }
 catch( ... ) {
  nc_close( ncid );
  throw;
  }
 nc_close( ncid );

 auto u = dynamic_cast< ThermalUnitBlock * >( b );
 if( ! u ) {
  delete b;
  throw( std::logic_error( "new_unit: no ThermalUnitBlock built" ) );
  }
 return( u );
 }

/*--------------------------------------------------------------------------*/
/// the BlockConfig of a formulation, read from the file cfg of the
/// directory of the configuration files [see read_cfg()], is applied to the
/// unit, whose abstract representation is then generated

std::string cfg_dir;  // where the configuration files are

Configuration * read_cfg( const std::string & cfg )
{
 auto prefix = Configuration::get_filename_prefix();
 Configuration::set_filename_prefix( std::string( cfg_dir ) );
 Configuration * c = nullptr;
 try {
  c = Configuration::deserialize( cfg );
  }
 catch( ... ) {
  Configuration::set_filename_prefix( std::move( prefix ) );
  throw;
  }
 Configuration::set_filename_prefix( std::move( prefix ) );
 if( ! c )
  throw( std::runtime_error( "test_power_limits: cannot read " + cfg_dir +
                             cfg ) );
 return( c );
 }

void generate( ThermalUnitBlock * u , const std::string & cfg ,
               bool reserves )
{
 if( reserves )
  u->set_reserve_vars( 3 );
 auto c = read_cfg( cfg );
 auto bc = dynamic_cast< BlockConfig * >( c );
 if( ! bc ) {
  delete c;
  throw( std::runtime_error( "test_power_limits: " + cfg +
                             " is not a BlockConfig" ) );
  }
 bc->apply( u );
 delete bc;
 u->generate_abstract_variables();
 u->generate_abstract_constraints();
 u->generate_objective();
 }

/*--------------------------------------------------------------------------*/
/*------------------------------ THE CHANGES -------------------------------*/
/*--------------------------------------------------------------------------*/

enum Kind { eMaxP , eAv , eInitP , eInitUD };

/// one call of a setter: through the Subset version if sub, through the
/// Range version [ first , first + values.size() ) otherwise
struct Step {
 Kind kind;
 bool sub;
 Subset idx;                  // the instants (Subset version)
 Index first = 0;             // the first instant (Range version)
 std::vector< double > val;
 };

struct Change {
 std::string name;
 std::vector< Step > steps;   // empty: the identity
 };

Step sub_step( Kind k , Subset idx , std::vector< double > v )
{
 return( Step{ k , true , std::move( idx ) , 0 , std::move( v ) } );
 }

Step rng_step( Kind k , Index first , std::vector< double > v )
{
 return( Step{ k , false , {} , first , std::move( v ) } );
 }

/// the instants of a Step
std::vector< Index > instants( const Step & s )
{
 if( s.sub )
  return( s.idx );
 std::vector< Index > r;
 for( Index i = 0 ; i < s.val.size() ; ++i )
  r.push_back( s.first + i );
 return( r );
 }

/// the data after the Step
void apply_data( UnitData & d , const Step & s )
{
 const auto ts = instants( s );
 switch( s.kind ) {
  case( eMaxP ):
   for( Index i = 0 ; i < ts.size() ; ++i )
    d.vec[ "MaxPower" ][ ts[ i ] ] = s.val[ i ];
   break;
  case( eAv ):
   if( ! d.vec.count( "Availability" ) )
    d.vec[ "Availability" ].assign( d.T , 1 );
   for( Index i = 0 ; i < ts.size() ; ++i )
    d.vec[ "Availability" ][ ts[ i ] ] = s.val[ i ];
   break;
  case( eInitP ):
   for( Index i = 0 ; i < ts.size() ; ++i )
    if( ts[ i ] == 0 )
     d.dbl[ "InitialPower" ] = s.val[ i ];
   break;
  case( eInitUD ):
   for( Index i = 0 ; i < ts.size() ; ++i )
    if( ts[ i ] == 0 )
     d.num[ "InitUpDownTime" ] = int( s.val[ i ] );
   break;
  }
 }

/// the Step on the unit
void apply_block( ThermalUnitBlock * u , const Step & s )
{
 std::vector< int > iv( s.val.begin() , s.val.end() );
 const Range r( s.first , s.first + s.val.size() );
 switch( s.kind ) {
  case( eMaxP ):
   if( s.sub )
    u->set_maximum_power( s.val.cbegin() , Subset( s.idx ) );
   else
    u->set_maximum_power( s.val.cbegin() , r );
   break;
  case( eAv ):
   if( s.sub )
    u->set_availability( s.val.cbegin() , Subset( s.idx ) );
   else
    u->set_availability( s.val.cbegin() , r );
   break;
  case( eInitP ):
   if( s.sub )
    u->set_initial_power( s.val.cbegin() , Subset( s.idx ) );
   else
    u->set_initial_power( s.val.cbegin() , r );
   break;
  case( eInitUD ):
   if( s.sub )
    u->set_init_updown_time( iv.cbegin() , Subset( s.idx ) );
   else
    u->set_init_updown_time( iv.cbegin() , r );
   break;
  }
 }

/*--------------------------------------------------------------------------*/
/// the changes for the unit d: the identity, MaxPower up and down at one
/// instant and on a range, below the given StartUpLimit and below
/// MinPower, Availability to 0 at 0 and at an inner instant, to 0.6, on a
/// range mixing 0 and 1, back from 0 to 1 in a sequence of three changes,
/// InitialPower up and to MinPower[ 0 ], InitUpDownTime across 0, longer
/// by 4 (the first free instant moves unless it is 0 already), to 1 (-1),
/// and longer by 4 followed by InitialPower up

std::vector< Change > changes_for( const UnitData & d )
{
 const Index T = d.T;
 const double M = d.at( "MaxPower" , 0 );
 const double m = d.at( "MinPower" , 0 );
 auto cst3 = []( double v ) { return( std::vector< double >( 3 , v ) ); };
 std::vector< Change > c;

 c.push_back( { "identity" , {} } );
 c.push_back( { "maxP up t2 (sub)" ,
                { sub_step( eMaxP , { 2 } , { 1.5 * M } ) } } );
 c.push_back( { "maxP down t2 (sub)" ,
                { sub_step( eMaxP , { 2 } ,
                            { std::max( m , 0.6 * M ) } ) } } );
 c.push_back( { "maxP up [1,4) (rng)" ,
                { rng_step( eMaxP , 1 , cst3( 1.3 * M ) ) } } );
 c.push_back( { "maxP down [T-3,T) (rng)" ,
                { rng_step( eMaxP , T - 3 ,
                            cst3( std::max( m , 0.7 * M ) ) ) } } );
 c.push_back( { "maxP up/down t1,t4 (sub)" ,
                { sub_step( eMaxP , { 4 , 1 } , { std::max( m , 0.6 * M ) ,
                                                  1.4 * M } ) } } );
 if( d.vec.count( "StartUpLimit" ) && ( d.at( "StartUpLimit" , 1 ) > m ) )
  c.push_back( { "maxP t1 below SU (sub)" ,
                 { sub_step( eMaxP , { 1 } ,
                             { ( m + d.at( "StartUpLimit" , 1 ) ) / 2 } )
                 } } );
 c.push_back( { "maxP t2 below minP (sub)" ,
                { sub_step( eMaxP , { 2 } , { m - 1 } ) } } );
 c.push_back( { "av 0 t0 (sub)" , { sub_step( eAv , { 0 } , { 0 } ) } } );
 c.push_back( { "av 0 t3 (rng)" , { rng_step( eAv , 3 , { 0 } ) } } );
 c.push_back( { "av 0.6 t2 (sub)" , { sub_step( eAv , { 2 } , { 0.6 } ) } } );
 c.push_back( { "av mix [1,5) (rng)" ,
                { rng_step( eAv , 1 , { 0 , 1 , 0.5 , 0 } ) } } );
 c.push_back( { "av 0 t4, maxP t1, av 1 t4 (seq)" ,
                { sub_step( eAv , { 4 } , { 0 } ) ,
                  rng_step( eMaxP , 1 , { 1.2 * M } ) ,
                  rng_step( eAv , 4 , { 1 } ) } } );
 if( d.num.count( "InitUpDownTime" ) &&
     ( d.num.at( "InitUpDownTime" ) > 0 ) ) {
  const double ip = d.dbl.count( "InitialPower" ) ?
                    d.dbl.at( "InitialPower" ) : m;
  c.push_back( { "initP +10 (sub)" ,
                 { sub_step( eInitP , { 0 } , { ip + 10 } ) } } );
  c.push_back( { "initP = minP (rng)" ,
                 { rng_step( eInitP , 0 , { m } ) } } );
  const int iud = d.num.at( "InitUpDownTime" );
  c.push_back( { "initUD -> -2 (sub)" ,
                 { sub_step( eInitUD , { 0 } , { -2 } ) } } );
  c.push_back( { "initUD +4 (rng)" ,
                 { rng_step( eInitUD , 0 , { double( iud + 4 ) } ) } } );
  if( iud != 1 )
   c.push_back( { "initUD -> 1 (sub)" ,
                  { sub_step( eInitUD , { 0 } , { 1 } ) } } );
  c.push_back( { "initUD +4, initP +10 (seq)" ,
                 { rng_step( eInitUD , 0 , { double( iud + 4 ) } ) ,
                   sub_step( eInitP , { 0 } , { ip + 10 } ) } } );
  }
 else {
  const int iud = d.num.count( "InitUpDownTime" ) ?
                  d.num.at( "InitUpDownTime" ) : 0;
  c.push_back( { "initUD -> 3 (sub)" ,
                 { sub_step( eInitUD , { 0 } , { 3 } ) } } );
  c.push_back( { "initUD -4 (rng)" ,
                 { rng_step( eInitUD , 0 , { double( iud - 4 ) } ) } } );
  if( iud != -1 )
   c.push_back( { "initUD -> -1 (sub)" ,
                  { sub_step( eInitUD , { 0 } , { -1 } ) } } );
  }
 return( c );
 }

/// whether the unit of the data d is on before the horizon, and the first
/// instant at which it may switch [see ThermalUnitBlock::first_free_instant()]
/// (the minimum times of the units here are below the bound that
/// deserialize() puts on them)
std::pair< bool , Index > init_state( const UnitData & d )
{
 const auto get = [ & ]( const std::string & n , int dflt ) {
  return( d.num.count( n ) ? d.num.at( n ) : dflt );
  };
 const int iud = get( "InitUpDownTime" , 0 );
 const int up = std::max( get( "MinUpTime" , 1 ) , 1 );
 const int dn = std::max( get( "MinDownTime" , 1 ) , 1 );
 int t0 = 0;
 if( iud > 0 )
  t0 = ( iud >= up ? 0 : up - iud );
 else
  t0 = ( - iud >= dn ? 0 : dn + iud );
 return( std::make_pair( iud > 0 , std::min( Index( t0 ) , d.T ) ) );
 }

/// true if the change of the data d has to be applied when the changed data
/// are valid: always, but a change of InitialPower, or of InitUpDownTime
/// that changes the initial state [see init_state()]
bool must_apply( const UnitData & d , const Change & c )
{
 UnitData e = d;
 for( const auto & s : c.steps ) {
  if( s.kind == eInitP )
   return( false );
  const auto before = init_state( e );
  apply_data( e , s );
  if( ( s.kind == eInitUD ) && ( init_state( e ) != before ) )
   return( false );
  }
 return( true );
 }

/*--------------------------------------------------------------------------*/
/*------------------------------ THE UNITS ---------------------------------*/
/*--------------------------------------------------------------------------*/

std::vector< double > cst( Index T , double v )
{
 return( std::vector< double >( T , v ) );
 }

/// prices that make the unit worth running at some instants only
std::vector< double > prices( Index T )
{
 static const double p[] = { -5 , 3 , -8 , -2 , 4 , -6 , -1 , -7 , 2 , -4 };
 std::vector< double > v( T );
 for( Index t = 0 ; t < T ; ++t )
  v[ t ] = p[ t % 10 ];
 return( v );
 }

UnitData base( const std::string & name , Index T , double minP ,
               double maxP , int initUD , double initP , int minUp ,
               int minDown )
{
 UnitData d;
 d.name = name;
 d.T = T;
 d.vec[ "MinPower" ] = cst( T , minP );
 d.vec[ "MaxPower" ] = cst( T , maxP );
 d.vec[ "LinearTerm" ] = prices( T );
 d.vec[ "ConstTerm" ] = cst( T , 10 );
 d.vec[ "StartUpCost" ] = cst( T , 50 );
 d.num[ "InitUpDownTime" ] = initUD;
 if( initUD > 0 )
  d.dbl[ "InitialPower" ] = initP;
 d.num[ "MinUpTime" ] = minUp;
 d.num[ "MinDownTime" ] = minDown;
 return( d );
 }

void ramps( UnitData & d , double ru , double rd )
{
 d.vec[ "DeltaRampUp" ] = cst( d.T , ru );
 d.vec[ "DeltaRampDown" ] = cst( d.T , rd );
 }

void limits( UnitData & d , double su , double sd )
{
 d.vec[ "StartUpLimit" ] = cst( d.T , su );
 d.vec[ "ShutDownLimit" ] = cst( d.T , sd );
 }

std::vector< UnitData > thermal_units( void )
{
 std::vector< UnitData > u;

 // no ramps, default limits, on before
 u.push_back( base( "plain" , 6 , 20 , 100 , 3 , 50 , 2 , 2 ) );

 // ramps with MinUpTime 4, given limits: in the T formulation the rows 5
 // and 6 of the maximum power appear and disappear with MaxPower
 {
  auto d = base( "ramps-T56" , 8 , 20 , 100 , -3 , 0 , 4 , 2 );
  ramps( d , 15 , 15 );
  limits( d , 40 , 40 );
  u.push_back( d );
  }

 // primary and secondary reserve with ramps, on before
 {
  auto d = base( "reserves" , 8 , 30 , 100 , 2 , 60 , 2 , 2 );
  ramps( d , 30 , 30 );
  limits( d , 50 , 50 );
  d.vec[ "PrimaryRho" ] = cst( d.T , 0.1 );
  d.vec[ "SecondaryRho" ] = cst( d.T , 0.15 );
  d.vec[ "PrimarySpinningReserveCost" ] = cst( d.T , -0.5 );
  d.vec[ "SecondarySpinningReserveCost" ] = cst( d.T , -0.3 );
  d.reserves = true;
  u.push_back( d );
  }

 // FixToMaximum, no ramps, prices of both signs: the unit is on at its
 // operational maximum at every instant where that is positive, whatever
 // the price, in the rows p_t >= opmax(t) of the abstract representation
 // and in the DP Solvers alike
 {
  auto d = base( "fixmax" , 6 , 10 , 80 , 2 , 80 , 1 , 1 );
  d.num[ "FixToMaximum" ] = 1;
  u.push_back( d );
  }

 // design variable (off before), with ramps and given limits
 {
  auto d = base( "design" , 6 , 20 , 90 , -1 , 0 , 2 , 1 );
  ramps( d , 40 , 40 );
  limits( d , 50 , 50 );
  d.dbl[ "InvestmentCost" ] = 100;
  u.push_back( d );
  }

 // ramps with default limits and a quadratic cost
 {
  auto d = base( "dflt-limits-ramps" , 8 , 25 , 100 , -2 , 0 , 3 , 2 );
  ramps( d , 25 , 25 );
  d.vec[ "QuadTerm" ] = cst( d.T , 0.01 );
  u.push_back( d );
  }

 // on before at MinPower[ 0 ] with the default ShutDownLimit: the
 // shut-down interval at 0, i.e., the arc (0,0), exists
 u.push_back( base( "arc00" , 6 , 30 , 100 , 5 , 30 , 2 , 2 ) );
 {
  auto d = base( "arc00-ramps" , 6 , 30 , 100 , 5 , 30 , 2 , 2 );
  ramps( d , 20 , 20 );
  u.push_back( d );
  }

 // small ramps, so that the numbers of ramp steps of SUSD move
 {
  auto d = base( "susd-steps" , 8 , 20 , 100 , 3 , 40 , 2 , 2 );
  ramps( d , 10 , 10 );
  limits( d , 30 , 30 );
  d.vec[ "QuadTerm" ] = cst( d.T , 0.02 );
  u.push_back( d );
  }

 // the initial state fixes the first two instants, on and off
 {
  auto d = base( "init_t-on" , 6 , 20 , 100 , 1 , 50 , 3 , 2 );
  ramps( d , 30 , 30 );
  limits( d , 40 , 40 );
  u.push_back( d );
  }
 {
  auto d = base( "init_t-off" , 6 , 20 , 100 , -1 , 0 , 2 , 3 );
  ramps( d , 30 , 30 );
  limits( d , 40 , 40 );
  u.push_back( d );
  }

 return( u );
 }

std::vector< UnitData > nuclear_units( void )
{
 std::vector< UnitData > u;

 // banded output, deep decreases, modulations of two instants
 {
  auto d = base( "nuc-bands" , 8 , 50 , 100 , 10 , 80 , 2 , 2 );
  ramps( d , 20 , 20 );
  limits( d , 60 , 60 );
  d.nuclear = true;
  d.num[ "ModulationTime" ] = 2;
  d.num[ "InitModulation" ] = 2;
  d.vec[ "ModulationDeltaRampUp" ] = cst( d.T , 5 );
  d.vec[ "ModulationDeltaRampDown" ] = cst( d.T , 5 );
  d.num[ "MaxModulationLength" ] = 2;
  d.bands = { 65 , 85 };
  d.vec[ "DeepDecreaseThreshold" ] = cst( d.T , 70 );
  d.vec[ "DeepDecreaseGradient" ] = cst( d.T , 15 );
  u.push_back( d );
  }

 // the plain operating rules, off before
 {
  auto d = base( "nuc-plain" , 8 , 40 , 100 , -2 , 0 , 2 , 2 );
  ramps( d , 30 , 30 );
  limits( d , 50 , 50 );
  d.nuclear = true;
  d.num[ "ModulationTime" ] = 2;
  d.num[ "InitModulation" ] = 1;
  d.vec[ "ModulationDeltaRampUp" ] = cst( d.T , 5 );
  d.vec[ "ModulationDeltaRampDown" ] = cst( d.T , 5 );
  u.push_back( d );
  }

 return( u );
 }

/*--------------------------------------------------------------------------*/
/*---------------------------- FORMULATIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

/* A formulation is the BlockConfig file that selects it: TUBCfg-<name>.txt
 * for the seven formulations of a ThermalUnitBlock, with and without the
 * Perspective Cuts (-PC), which serve a NuclearUnitBlock as well, and
 * NUBCfg-<name>.txt for those that also have the bits of a
 * NuclearUnitBlock (TR TightRules, TM TightRamp, TC TightCuts). */

struct Form {
 std::string name;       // for the report, also the suffix of the file
 std::string cfg;        // the BlockConfig file
 bool pc;                // with the Perspective Cuts
 bool reduced;           // in the reduced matrix
 };

std::vector< Form > thermal_forms( void )
{
 std::vector< Form > f;
 for( bool pc : { false , true } )
  for( std::string n : { "3bin" , "T" , "pt" , "DP" , "SU" , "SD" ,
                         "SUSD" } ) {
   const auto name = n + ( pc ? "-PC" : "" );
   f.push_back( { name , "TUBCfg-" + name + ".txt" , pc ,
                  ( n == "3bin" ) || ( n == "T" ) || ( n == "SUSD" ) ||
                  ( ( n == "DP" ) && ( ! pc ) ) } );
   }
 return( f );
 }

std::vector< Form > nuclear_forms( void )
{
 std::vector< Form > f = { { "T" , "TUBCfg-T.txt" , false , true } ,
                           { "T-PC" , "TUBCfg-T-PC.txt" , true , true } ,
                           { "3bin" , "TUBCfg-3bin.txt" , false , true } };
 for( std::string n : { "T-TR" , "T-TM" , "T-TR-TC" , "T-TR-TM-TC" ,
                        "3bin-TR" , "3bin-TR-TM" } )
  f.push_back( { n , "NUBCfg-" + n + ".txt" , false , true } );
 return( f );
 }

/*--------------------------------------------------------------------------*/
/*------------------------------- SOLVING ----------------------------------*/
/*--------------------------------------------------------------------------*/
/* The Solver come from the BlockSolverConfig LPCmpBSCfg[-nuc].txt: the
 * first is the :MILPSolver, which writes the model in the file of its
 * strOutputFile, the others are the dynamic programming Solver. */

BlockSolverConfig * read_bsc( bool nuclear )
{
 auto c = read_cfg( nuclear ? "LPCmpBSCfg-nuc.txt" : "LPCmpBSCfg.txt" );
 auto bsc = dynamic_cast< BlockSolverConfig * >( c );
 if( ! bsc ) {
  delete c;
  throw( std::runtime_error( "test_power_limits: cannot read the "
                             "BlockSolverConfig in " + cfg_dir ) );
  }
 return( bsc );
 }

/// attaches the Solver of the BlockSolverConfig, which is returned
/// clear()-ed, ready to detach them
std::vector< std::string > solver_names;  // of the last BlockSolverConfig

BlockSolverConfig * attach( ThermalUnitBlock * u , bool nuclear )
{
 auto bsc = read_bsc( nuclear );
 solver_names = bsc->get_SolverNames();
 bsc->apply( u );
 bsc->clear();
 return( bsc );
 }

void detach( ThermalUnitBlock * u , BlockSolverConfig * bsc )
{
 bsc->apply( u );
 delete bsc;
 }

struct Values {
 double milp = std::nan( "" );        // INF if infeasible, NaN if none
 std::vector< double > dp;            // the same for each DP Solver
 std::vector< std::string > dp_name;
 bool dp_sched_ok = true;             // DP schedules satisfy the rows
 std::string dp_sched_who;
 };

double value_of( Solver * s , std::string & err )
{
 try {
  const auto st = s->compute();
  if( st == Solver::kOK )
   return( s->get_ub() );
  if( st == Solver::kInfeasible )
   return( lp_compare::INF );
  err = "status " + std::to_string( st );
  }
 catch( std::exception & e ) {
  err = e.what();
  }
 return( std::nan( "" ) );
 }

/// solves with every Solver; the model of the :MILPSolver goes in lp_file
/// (removed if empty); with check_sched, the schedule of every DP Solver is
/// written in the Variable and checked against the rows of the unit
Values solve_all( ThermalUnitBlock * u , const std::string & lp_file ,
                  bool check_sched , std::string & err )
{
 Values v;
 const auto & sl = u->get_registered_solvers();
 auto it = sl.begin();
 Solver * milp = * it;
 v.milp = value_of( milp , err );
 const std::string written = milp->get_str_par(
                                  milp->str_par_str2idx( "strOutputFile" ) );
 // CPXMILPSolver writes nothing when it finds a row with lhs > rhs
 if( ! std::filesystem::exists( written ) )
  err += " no model written;";
 else
  if( lp_file.empty() )
   std::remove( written.c_str() );
  else
   std::filesystem::rename( written , lp_file );

 for( ++it ; it != sl.end() ; ++it ) {
  Solver * dp = * it;
  const auto i = v.dp.size() + 1;
  const std::string name = i < solver_names.size() ?
                           solver_names[ i ] : "DP" + std::to_string( i );
  std::string e;
  v.dp.push_back( value_of( dp , e ) );
  v.dp_name.push_back( name );
  if( ! e.empty() )
   err += " DP " + name + ": " + e;
  if( check_sched && ( v.dp.back() < lp_compare::INF ) &&
      dp->has_var_solution() ) {
   dp->get_var_solution();
   SimpleConfiguration< double > tol( 1e-7 );
   if( ! u->is_feasible( true , & tol ) ) {
    v.dp_sched_ok = false;
    v.dp_sched_who += " " + name;
    }
   }
  }
 return( v );
 }

bool same_value( double a , double b )
{
 if( std::isnan( a ) || std::isnan( b ) )
  return( false );
 return( lp_compare::same( a , b , 1e-6 ) );
 }

/*--------------------------------------------------------------------------*/
/// the data of A and B that the changes concern, compared: an empty string
/// if they are equal

std::string data_delta( ThermalUnitBlock * a , ThermalUnitBlock * b )
{
 std::ostringstream s;
 auto cmp = [ & ]( const std::string & n , double x , double y ) {
  if( ! lp_compare::same( x , y , 1e-9 ) )
   s << " " << n << " A " << x << " B " << y << ";";
  };
 const Index T = a->get_time_horizon();
 for( Index t = 0 ; t < T ; ++t ) {
  const auto ts = "[" + std::to_string( t ) + "]";
  cmp( "MaxPower" + ts , a->get_max_power( t ) , b->get_max_power( t ) );
  cmp( "Availability" + ts , a->get_availability( t ) ,
       b->get_availability( t ) );
  cmp( "opmin" + ts , a->get_operational_min_power( t ) ,
       b->get_operational_min_power( t ) );
  cmp( "opmax" + ts , a->get_operational_max_power( t ) ,
       b->get_operational_max_power( t ) );
  if( ( a->get_start_up_limit().size() > t ) &&
      ( b->get_start_up_limit().size() > t ) )
   cmp( "StartUpLimit" + ts , a->get_start_up_limit()[ t ] ,
        b->get_start_up_limit()[ t ] );
  if( ( a->get_shut_down_limit().size() > t ) &&
      ( b->get_shut_down_limit().size() > t ) )
   cmp( "ShutDownLimit" + ts , a->get_shut_down_limit()[ t ] ,
        b->get_shut_down_limit()[ t ] );
  }
 cmp( "InitialPower" , a->get_initial_power() , b->get_initial_power() );
 cmp( "InitUpDownTime" , a->get_init_up_down_time() ,
      b->get_init_up_down_time() );
 return( s.str() );
 }

/*--------------------------------------------------------------------------*/
/*------------------------------- ONE CASE ---------------------------------*/
/*--------------------------------------------------------------------------*/

struct Outcome {
 std::string kind;    // EXACT, DIFF, REFUSED, REFUSED-VALID, REFUSED-DIRTY,
                      // ACCEPTED-INVALID, NOGEN, ERROR
 bool fail = false;
 std::string groups;  // the groups of rows (columns) that differ
 std::string detail;
 };

/// a unit read from d and generated, or nullptr (and why) if it cannot be
ThermalUnitBlock * make( const UnitData & d , const Form & wf ,
                         std::string & why )
{
 ThermalUnitBlock * u = nullptr;
 try {
  u = new_unit( d );
  generate( u , wf.cfg , d.reserves );
  return( u );
  }
 catch( std::exception & e ) {
  why = e.what();
  }
 delete u;
 return( nullptr );
 }

std::string tmp_lp( const std::string & tag )
{
 return( out_dir + "/" + tag + ".lp" );
 }

Outcome run_case( const UnitData & d , const Form & wf , const Change & c ,
                  const std::string & tag )
{
 Outcome o;
 const bool exact_milp = ! ( wf.pc && ( d.quad() > 0 ) );
 std::string why;

 // A: read, generate, solve once, change, solve again - - - - - - - - - - -
 auto a = make( d , wf , why );
 if( ! a ) {
  o.kind = "NOGEN";
  o.detail = why;
  return( o );
  }
 auto bsa = attach( a , d.nuclear );
 std::string err;
 solve_all( a , "" , false , err );

 UnitData dref = d;  // the data up to the refused step
 UnitData drfs = d;  // the data the refused step would give
 std::string refused;
 for( const auto & s : c.steps ) {
  try {
   apply_block( a , s );
   }
  catch( std::exception & e ) {
   refused = e.what();
   apply_data( drfs , s );
   break;
   }
  apply_data( dref , s );
  apply_data( drfs , s );
  }
 const auto lpa = tmp_lp( tag + "_A" );
 // the schedule of every DP Solver is written in the Variable, in every
 // formulation [see ThermalUnitBlock::set_solution()], and checked against
 // the rows of the unit
 const bool sched = true;
 auto va = solve_all( a , lpa , sched , err );

 // validity of the refused step: can the data it gives be read and
 // generated? (a step of a sequence may give invalid data, e.g., the
 // availability 0 with a StartUpLimit given, which a later step would
 // make valid again)
 bool valid = true;
 std::string invalid_why;
 if( ! refused.empty() ) {
  auto f = make( drfs , wf , invalid_why );
  valid = ( f != nullptr );
  delete f;
  }

 // B: the unit read afresh from the data A should have - - - - - - - - - - -
 auto b = make( dref , wf , why );
 if( ! b ) {
  // the data that A took cannot be read: A accepted invalid data
  o.kind = "ACCEPTED-INVALID";
  o.fail = true;
  o.detail = why;
  detach( a , bsa );
  delete a;
  return( o );
  }
 auto bsb = attach( b , d.nuclear );
 const auto lpb = tmp_lp( tag + "_B" );
 auto vb = solve_all( b , lpb , sched , err );

 // compare - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 std::ostringstream det;
 lp_compare::Diff diff;
 bool lp_same = false;
 if( ! ( std::filesystem::exists( lpa ) && std::filesystem::exists( lpb ) ) )
  det << " the :MILPSolver wrote no model of "
      << ( std::filesystem::exists( lpa ) ? "B" : "A" ) << ";";
 else
  try {
   diff = lp_compare::compare( lp_compare::read_lp( lpa ) ,
                               lp_compare::read_lp( lpb ) );
   lp_same = diff.empty();
   }
  catch( std::exception & e ) {
   o.kind = "ERROR";
   o.fail = true;
   det << " " << e.what();
   }
 const auto dd = data_delta( a , b );

 const bool milp_same = same_value( va.milp , vb.milp );
 bool dp_ok = true;
 if( exact_milp )
  for( std::size_t i = 0 ; i < va.dp.size() ; ++i )
   if( ( ! same_value( va.dp[ i ] , va.milp ) ) ||
       ( ( i < vb.dp.size() ) && ! same_value( vb.dp[ i ] , vb.milp ) ) ) {
    dp_ok = false;
    det << " DP" << i << " A " << va.dp[ i ] << " B "
        << ( i < vb.dp.size() ? vb.dp[ i ] : std::nan( "" ) )
        << " vs MILP A " << va.milp << " B " << vb.milp << ";";
    }

 if( o.kind.empty() ) {
  if( refused.empty() )
   o.kind = ( lp_same && dd.empty() ) ? "EXACT" : "DIFF";
  else
   if( ! valid )
    o.kind = ( lp_same && dd.empty() ) ? "REFUSED" : "REFUSED-DIRTY";
   else
    if( must_apply( d , c ) )
     o.kind = "REFUSED-VALID";
    else
     o.kind = ( lp_same && dd.empty() ) ? "REFUSED" : "REFUSED-DIRTY";
  }
 if( ( o.kind != "EXACT" ) && ( o.kind != "REFUSED" ) )
  o.fail = true;
 if( ( ! milp_same ) || ( ! dp_ok ) || ( ! va.dp_sched_ok ) ||
     ( ! vb.dp_sched_ok ) )
  o.fail = true;

 o.groups = lp_same || ! diff.empty() ? diff.groups() : "(no model)";
 if( ! dd.empty() )
  o.groups += std::string( o.groups.empty() ? "" : "," ) + "data";
 if( ! milp_same )
  det << " MILP A " << va.milp << " B " << vb.milp << ";";
 if( ! va.dp_sched_ok )
  det << " DP schedule violates the rows of A:" << va.dp_sched_who << ";";
 if( ! vb.dp_sched_ok )
  det << " DP schedule violates the rows of B:" << vb.dp_sched_who << ";";
 if( ! refused.empty() )
  det << " refused: " << refused << ";";
 if( ! valid )
  det << " (invalid data: " << invalid_why << ")";
 if( ! dd.empty() )
  det << "\n  data:" << dd;
 if( ! err.empty() )
  det << "\n  solver: " << err;
 if( ! lp_same )
  det << "\n" << diff.report();
 o.detail = det.str();

 detach( a , bsa );
 detach( b , bsb );
 delete a;
 delete b;

 return( o );
 }

/*--------------------------------------------------------------------------*/
/// the case run in a child process where there is one, so that a crash of
/// the code under test is reported as such and does not stop the sweep

std::string pack( const Outcome & o )
{
 const char sep = '\x1f';
 return( o.kind + sep + ( o.fail ? "1" : "0" ) + sep + o.groups + sep +
         o.detail );
 }

Outcome unpack( const std::string & s )
{
 Outcome o;
 std::vector< std::string > f( 1 );
 for( char c : s )
  if( ( c == '\x1f' ) && ( f.size() < 4 ) )
   f.emplace_back();
  else
   f.back() += c;
 f.resize( 4 );
 o.kind = f[ 0 ];
 o.fail = ( f[ 1 ] == "1" );
 o.groups = f[ 2 ];
 o.detail = f[ 3 ];
 return( o );
 }

Outcome run_case_safe( const UnitData & d , const Form & wf ,
                       const Change & c ,
                       const std::string & tag )
{
 auto guarded = [ & ]() {
  try {
   return( run_case( d , wf , c , tag ) );
   }
  catch( std::exception & e ) {
   Outcome o;
   o.kind = "ERROR";
   o.fail = true;
   o.detail = e.what();
   return( o );
   }
  };

#if defined( __unix__ ) || defined( __APPLE__ )
 int fd[ 2 ];
 std::cout.flush();
 std::cerr.flush();
 if( pipe( fd ) == 0 ) {
  const pid_t pid = fork();
  if( pid == 0 ) {               // the child
   close( fd[ 0 ] );
   const auto msg = pack( guarded() );
   std::size_t done = 0;
   while( done < msg.size() ) {
    const auto w = write( fd[ 1 ] , msg.data() + done , msg.size() - done );
    if( w <= 0 )
     break;
    done += std::size_t( w );
    }
   close( fd[ 1 ] );
   std::cout.flush();
   std::cerr.flush();
   _exit( 0 );
   }
  if( pid > 0 ) {                // the parent
   close( fd[ 1 ] );
   std::string msg;
   char buf[ 4096 ];
   ssize_t r;
   while( ( r = read( fd[ 0 ] , buf , sizeof( buf ) ) ) > 0 )
    msg.append( buf , std::size_t( r ) );
   close( fd[ 0 ] );
   int status = 0;
   waitpid( pid , & status , 0 );
   if( WIFSIGNALED( status ) || msg.empty() ) {
    Outcome o;
    o.kind = "CRASH";
    o.fail = true;
    if( WIFSIGNALED( status ) )
     o.detail = " signal " + std::to_string( WTERMSIG( status ) );
    return( o );
    }
   return( unpack( msg ) );
   }
  close( fd[ 0 ] );
  close( fd[ 1 ] );
  }
#endif
 return( guarded() );
 }

/*--------------------------------------------------------------------------*/
/*------------------------------- THE SWEEP --------------------------------*/
/*--------------------------------------------------------------------------*/

/// counts of the outcomes, by formulation x change
struct Summary {
 std::map< std::string , std::map< std::string , int > > by_form;
 std::map< std::string , std::map< std::string , int > > groups_by_form;
 std::map< std::string , int > total;
 };

/// the name of a case as a file name
std::string file_name( const std::string & label )
{
 std::string f;
 for( char c : label )
  f += std::isalnum( static_cast< unsigned char >( c ) ) || ( c == '-' ) ?
       c : '_';
 return( f );
 }

void sweep( const std::vector< UnitData > & units ,
            const std::vector< Form > & forms , bool only_identity ,
            bool full , Summary & sum )
{
 for( const auto & d : units )
  for( const auto & wf : forms ) {
   if( ( ! full ) && ( ! wf.reduced ) )
    continue;
   const auto fn = wf.name;
   for( const auto & c : changes_for( d ) ) {
    if( only_identity && ! c.steps.empty() )
     continue;
    const auto tag = "case" + std::to_string( n_case++ );
    const auto label = d.name + " | " + fn + " | " + c.name;
    if( ( ! only.empty() ) && ( label.find( only ) == std::string::npos ) )
     continue;
    const auto o = run_case_safe( d , wf , c , tag );
    // the models of every case go in keep_dir, if given, those of the
    // first failed cases are kept otherwise
    if( ! keep_dir.empty() )
     for( const auto & s : { "_A" , "_B" } ) {
      const auto from = tmp_lp( tag + s );
      if( std::filesystem::exists( from ) )
       std::filesystem::rename( from , keep_dir + "/" + file_name( label ) +
                                       s + ".lp" );
      }
    else
     if( ( o.fail || ! only.empty() ) && ( n_kept < max_kept ) )
      ++n_kept;
     else
      for( const auto & s : { "_A" , "_B" } )
       std::remove( tmp_lp( tag + s ).c_str() );
    ++sum.total[ o.kind ];
    ++sum.by_form[ fn ][ o.kind ];
    if( ! o.groups.empty() ) {
     std::stringstream gs( o.groups );
     std::string g;
     while( std::getline( gs , g , ',' ) )
      ++sum.groups_by_form[ fn ][ g ];
     }
    if( o.fail )
     ++n_fail;
    if( o.fail || verbose )
     std::cout << ( o.fail ? "FAIL " : "ok   " ) << tag << " " << label
               << " | " << o.kind
               << ( o.groups.empty() ? "" : " [" + o.groups + "]" )
               << ( o.fail ? " " + o.detail : "" ) << std::endl;
    }
   }
 }

void print_summary( const Summary & sum )
{
 std::cout << "power limits: outcomes by formulation" << std::endl;
 for( const auto & f : sum.by_form ) {
  std::cout << "  " << f.first << ":";
  for( const auto & k : f.second )
   std::cout << " " << k.first << " " << k.second;
  std::cout << std::endl;
  const auto it = sum.groups_by_form.find( f.first );
  if( it != sum.groups_by_form.end() ) {
   std::cout << "    differing groups (cases):";
   for( const auto & g : it->second )
    std::cout << " " << g.first << " " << g.second;
   std::cout << std::endl;
   }
  }
 std::cout << "  total:";
 for( const auto & k : sum.total )
  std::cout << " " << k.first << " " << k.second;
 std::cout << std::endl;
 }

}  // end( anonymous namespace )

/*--------------------------------------------------------------------------*/
/*---------------------------- ENTRY POINT ---------------------------------*/
/*--------------------------------------------------------------------------*/
/// runs the test with the configuration files in dir and the options in
/// args: "--full" the whole matrix (the reduced one otherwise, see the
/// description of the file), "-i" only the identity, i.e., each unit
/// compared with itself read again, "-v" one line per case, "-k s" only
/// the cases whose "unit | formulation | change" contains s (e.g.,
/// "arc00 | pt | av"), "-d dir" the models of every case kept in dir;
/// returns the number of failed cases, 0 if the test is skipped

int test_power_limits( const std::string & dir ,
                       const std::vector< std::string > & args )
{
 cfg_dir = dir;
 if( ( ! cfg_dir.empty() ) && ( cfg_dir.back() != '/' ) )
  cfg_dir += '/';
 bool only_identity = false , full = false;
 verbose = false;
 only.clear();
 keep_dir.clear();
 for( std::size_t i = 0 ; i < args.size() ; ++i )
  if( args[ i ] == "--full" )
   full = true;
  else
   if( args[ i ] == "-i" )
    only_identity = true;
   else
    if( args[ i ] == "-v" )
     verbose = true;
    else
     if( ( args[ i ] == "-k" ) && ( i + 1 < args.size() ) )
      only = args[ ++i ];
     else
      if( ( args[ i ] == "-d" ) && ( i + 1 < args.size() ) )
       keep_dir = args[ ++i ];
      else
       throw( std::invalid_argument( "test_power_limits: unknown option " +
                                     args[ i ] ) );
 n_fail = 0;
 n_kept = 0;
 n_case = 0;

 // skip if a Solver of the configuration is not in the factory
 for( bool nuc : { false , true } ) {
  auto bsc = read_bsc( nuc );
  for( const auto & n : bsc->get_SolverNames() )
   if( ! Solver::has_Solver( n ) ) {
    std::cout << "test_power_limits: " << n << " not in the factory, "
              << "skipped" << std::endl;
    delete bsc;
    return( 0 );
    }
  delete bsc;
  }

 out_dir = "power_limits_lp";
 std::filesystem::create_directories( out_dir );
 if( ! keep_dir.empty() )
  std::filesystem::create_directories( keep_dir );

 Summary sum;
 sweep( thermal_units() , thermal_forms() , only_identity , full , sum );
 sweep( nuclear_units() , nuclear_forms() , only_identity , full , sum );
 std::cout << "power limits (" << ( full ? "full" : "reduced" )
           << " matrix)" << std::endl;
 print_summary( sum );

 if( n_kept == 0 )
  std::filesystem::remove_all( out_dir );
 else
  std::cout << "the models of " << n_kept << " failed cases are in "
            << out_dir << std::endl;
 return( n_fail );
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- STANDALONE MAIN -------------------------------*/
/*--------------------------------------------------------------------------*/
/* With POWER_LIMITS_MAIN defined the file is a program of its own:
 *   test_power_limits <dir of the configuration files> [options]
 * with the options of test_power_limits(). */

#ifdef POWER_LIMITS_MAIN

int main( int argc , char ** argv )
{
 if( argc < 2 ) {
  std::cerr << "usage: " << argv[ 0 ] << " <dir of the configuration files> "
            << "[--full] [-i] [-v] [-k s] [-d dir]" << std::endl;
  return( 1 );
  }
 const int f = test_power_limits( argv[ 1 ] ,
                                  std::vector< std::string >( argv + 2 ,
                                                              argv + argc ) );
 std::cout << ( f ? std::to_string( f ) + " cases FAILED" :
                    std::string( "all cases passed" ) ) << std::endl;
 return( f ? 1 : 0 );
 }

#endif

/*--------------------------------------------------------------------------*/
/*--------------------- End File test_power_limits.cpp ---------------------*/
/*--------------------------------------------------------------------------*/
