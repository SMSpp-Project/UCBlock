/*--------------------------------------------------------------------------*/
/*--------------------- File ThermalUnitBlock.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the ThermalUnitBlock class.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Ali Ghezelsoflu \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Rafael Durbano Lobato \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Tiziano Bacci \n
 *         Istituto di Analisi di Sistemi e Informatica "Antonio Ruberti" \n
 *         Consiglio Nazionale delle Ricerche \n
 *
 * \copyright &copy; by Antonio Frangioni, Ali Ghezelsoflu,
 *                      Rafael Durbano Lobato, Donato Meoli, Tiziano Bacci
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <list>
#include <map>
#include <numeric>
#include <set>
#include <tuple>

#include "LinearFunction.h"

#include "DQuadFunction.h"

#include "ThermalUnitBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

template< typename T >
static bool identical( std::vector< T > & vec , const Block::Subset sbst ,
                       typename std::vector< T >::const_iterator it ) {
 // returns true if the sub-vector of vec[] corresponding to the indices
 // in sbst is identical to the vector starting at it
 for( auto t : sbst )
  if( vec[ t ] != *( it++ ) )
   return( false );
 return( true );
}

/*--------------------------------------------------------------------------*/

template< typename T >
static void assign( std::vector< T > & vec , const Block::Subset sbst ,
                    typename std::vector< T >::const_iterator it ) {
 // assign to the sub-vector of vec[] corresponding to the indices in sbst
 // the values found in vector starting at it
 for( auto t : sbst )
  vec[ t ] = *( it++ );
}

/*--------------------------------------------------------------------------*/
/* Sorts the indices of a Subset that the caller gives unordered, together
 * with the values that go with them: the values are copied, in the new
 * order, in sorted, and values is set to its beginning. The sort is
 * stable, so that an index given more than once takes its last value as it
 * would unsorted. */

static void sort_by_index( Block::Subset & subset ,
                           std::vector< double >::const_iterator & values ,
                           std::vector< double > & sorted )
{
 std::vector< std::pair< Block::Index , double > > iv( subset.size() );
 for( Block::Index i = 0 ; i < subset.size() ; ++i )
  iv[ i ] = std::make_pair( subset[ i ] , *( values + i ) );
 std::stable_sort( iv.begin() , iv.end() ,
                   []( const auto & a , const auto & b ) {
                    return( a.first < b.first ); } );
 sorted.resize( iv.size() );
 for( Block::Index i = 0 ; i < iv.size() ; ++i ) {
  subset[ i ] = iv[ i ].first;
  sorted[ i ] = iv[ i ].second;
  }
 values = sorted.cbegin();
 }

/*--------------------------------------------------------------------------*/

Block::Subset subset_add( const Block::Subset & sbst , Block::Index dlt ) {
 Block::Subset ret = sbst;
 for( auto & t : ret )
  t += dlt;
 return( ret );
}

/*--------------------------------------------------------------------------*/

Block::Subset subset_sbtrct( const Block::Subset & sbst , Block::Index dlt ) {
 Block::Subset ret = sbst;
 for( auto & t : ret )
  t -= dlt;
 return( ret );
}

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register ThermalUnitBlock to the Block factory
SMSpp_insert_in_factory_cpp_1( ThermalUnitBlock );

// register ThermalUnitBlockSolution to the Solution factory
SMSpp_insert_in_factory_cpp_0( ThermalUnitBlockSolution );

/*--------------------------------------------------------------------------*/
/*-------------------------------- FUNCTIONS -------------------------------*/
/*--------------------------------------------------------------------------*/
/* The minimum up (down) time m of a unit with initial up/down time init over
 * a horizon of T instants, taken at least 1 and at most T + max( 1 , k ),
 * where k is the number of instants the unit has been on (off) before the
 * horizon, i.e. init if up and the unit is on, - init if down and the unit
 * is off, 0 otherwise. The unit has then to stay as it is for m - k more
 * instants, which is at most T when m is taken to be the bound, i.e. the
 * unit never switches within the horizon, and anything larger means the
 * same; a smaller bound would free the unit before the minimum time. */

static Block::Index clamp_min_time( Block::Index m , Block::Index T ,
                                    int init , bool up )
{
 const Block::Index k = up ? ( init > 0 ? Block::Index( init ) : 0 )
                           : ( init <= 0 ? Block::Index( - init ) : 0 );
 return( std::min( std::max( m , Block::Index( 1 ) ) ,
                   T + std::max( k , Block::Index( 1 ) ) ) );
 }

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS OF ThermalUnitBlock ----------------------*/
/*--------------------------------------------------------------------------*/

ThermalUnitBlock::~ThermalUnitBlock()
{
 Constraint::clear( CommitmentDesign_Const );
 Constraint::clear( StartUp_ShutDown_Variables_Const );
 Constraint::clear( StartUp_Const );
 Constraint::clear( ShutDown_Const );
 Constraint::clear( RampUp_Const );
 Constraint::clear( RampDown_Const );
 Constraint::clear( MinPower_Const );
 Constraint::clear( MaxPower_Const );
 Constraint::clear( PrimaryRho_Const );
 Constraint::clear( SecondaryRho_Const );
 Constraint::clear( Reserve_Const );

 Constraint::clear( Eq_ActivePower_Const );
 Constraint::clear( Eq_Commitment_Const );
 Constraint::clear( Eq_StartUp_Const );
 Constraint::clear( Eq_ShutDown_Const );
 Constraint::clear( Network_Const );

 Constraint::clear( Init_PC_Const );
 Constraint::clear( Eq_PC_Const );

 Constraint::clear( PC_cuts );
 Constraint::clear( ShutDownZero_Const );
 Constraint::clear( MaxPower5_Const.rows );
 Constraint::clear( MaxPower6_Const.rows );
 Constraint::clear( RampUpSUSD_Const.rows );
 Constraint::clear( RampDownSUSD_Const.rows );

 Constraint::clear( Commitment_bound_Const );
 Constraint::clear( StartUp_Binary_bound_Const );
 Constraint::clear( ShutDown_Binary_bound_Const );

 Constraint::clear( Commitment_fixed_to_One_Const );

 objective.clear();
}

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::deserialize( const netCDF::NcGroup & group )
{
 // Deserialize data from the base class
 UnitBlock::deserialize( group );

 // Mandatory variables

 ::deserialize( group , "MaxPower" , f_time_horizon , v_MaxPower ,
                false , true , v_change_intervals );

 // Optional variables

 // a continuous design is not representable here, as it would multiply the
 // binary commitment variable in the maximum power constraint
 for( const auto & design : { "MaxCapacityDesign" , "MinCapacityDesign" } )
  if( ! group.getVar( design ).isNull() )
   throw( std::logic_error(
    std::string( "ThermalUnitBlock::deserialize: " ) + design +
    " is not supported. A continuous design variable would multiply the "
    "binary commitment variable in the maximum power constraint, making the "
    "model bilinear. The available options are: \"InvestmentCost\", for the "
    "binary decision of building the unit; \"Scale\" = N, for a fleet of N "
    "identical modules committed together; or N replicated ThermalUnitBlock, "
    "for an integer install with independent commitments." ) );

 ::deserialize( group , f_InvestmentCost , "InvestmentCost" );

 ::deserialize( group , f_Capacity , "Capacity" );

 ::deserialize( group , f_scale , "Scale" );

 if( ::deserialize( group , f_MinUpTime , "MinUpTime" ) )
  f_MinUpTime = std::max( f_MinUpTime , static_cast< Index >( 1 ) );

 if( ::deserialize( group , f_MinDownTime , "MinDownTime" ) )
  f_MinDownTime = std::max( f_MinDownTime , static_cast< Index >( 1 ) );

 ::deserialize( group , f_InitialPower , "InitialPower" );

 if( ! ::deserialize( group , f_InitUpDownTime , "InitUpDownTime" ) ) {
  if( f_InitialPower == 0 )
   f_InitUpDownTime = -f_MinDownTime;
  else
   f_InitUpDownTime = f_MinUpTime;
 }

 // the bound of the minimum times depends on the initial state
 f_MinUpTime = clamp_min_time( f_MinUpTime , f_time_horizon ,
                               f_InitUpDownTime , true );
 f_MinDownTime = clamp_min_time( f_MinDownTime , f_time_horizon ,
                                 f_InitUpDownTime , false );

 if( ! ::deserialize( group , "MinPower" , f_time_horizon , v_MinPower ,
                      true , true , v_change_intervals ) )
  v_MinPower.resize( f_time_horizon );

 if( ! ::deserialize( group , "Availability" , f_time_horizon ,
		      v_Availability , true , true , v_change_intervals ) )
  v_Availability.resize( f_time_horizon , 1 );

 if( ! ::deserialize( group , "LinearTerm" , f_time_horizon , v_LinearTerm ,
                      true , true , v_change_intervals ) )
  v_LinearTerm.resize( f_time_horizon );

 // the spinning-reserve linear costs (the dual of the reserve demand pushed
 // onto the unit, the counterpart of LinearTerm); optional -- left empty when
 // absent, i.e. when the unit prices no reserve
 ::deserialize( group , "PrimarySpinningReserveCost" , f_time_horizon ,
                v_PrimarySpinningReserveCost , true , true , v_change_intervals );
 ::deserialize( group , "SecondarySpinningReserveCost" , f_time_horizon ,
                v_SecondarySpinningReserveCost , true , true ,
                v_change_intervals );

 // the reactive-power linear cost (the dual of the reactive nodal balance
 // pushed onto the unit, the counterpart of LinearTerm); optional -- left
 // empty when absent, i.e. when the unit prices no reactive power
 ::deserialize( group , "ReactiveLinearTerm" , f_time_horizon ,
                v_ReactiveLinearTerm , true , true , v_change_intervals );

 if( ! ::deserialize( group , "QuadTerm" , f_time_horizon , v_QuadTerm ,
                      true , true , v_change_intervals ) )
  v_QuadTerm.resize( f_time_horizon );

 if( ! ::deserialize( group , "ConstTerm" , f_time_horizon , v_ConstTerm ,
                      true , true , v_change_intervals ) )
  v_ConstTerm.resize( f_time_horizon );

 if( ! ::deserialize( group , "StartUpCost" , f_time_horizon , v_StartUpCost ,
                      true , true , v_change_intervals ) )
  v_StartUpCost.resize( f_time_horizon );

 // a unit that pays nothing to shut down keeps the vector empty, and then
 // the Objective has no term for the shut-down variables
 if( ::deserialize( group , "ShutDownCost" , f_time_horizon , v_ShutDownCost ,
                    true , true , v_change_intervals ) &&
     std::all_of( v_ShutDownCost.begin() , v_ShutDownCost.end() ,
                  []( double cst ) { return( cst == 0 ); } ) )
  v_ShutDownCost.clear();

 ::deserialize( group , "DeltaRampUp" , f_time_horizon , v_DeltaRampUp ,
                true , true , v_change_intervals );

 ::deserialize( group , "DeltaRampDown" , f_time_horizon , v_DeltaRampDown ,
                true , true , v_change_intervals );

 ::deserialize( group , "FixedConsumption" , f_time_horizon ,
                v_FixedConsumption , true , true , v_change_intervals );

 ::deserialize( group , "InertiaCommitment" , f_time_horizon ,
                v_InertiaCommitment , true , true , v_change_intervals );

 if( ! ( f_ignore_netcdf_vars & 1 ) ) {
  ::deserialize( group , "PrimaryRho" , f_time_horizon , v_PrimaryRho ,
		 true , true , v_change_intervals );
  if( std::all_of( v_PrimaryRho.begin() , v_PrimaryRho.end() ,
		   []( double i ) { return( i == 0 ); } ) )
   v_PrimaryRho.clear();

  ::deserialize( group , "SecondaryRho" , f_time_horizon , v_SecondaryRho ,
                 true , true , v_change_intervals );
  if( std::all_of( v_SecondaryRho.begin() , v_SecondaryRho.end() ,
		   []( double i ) { return( i == 0 ); } ) )
   v_SecondaryRho.clear();
  }

 if( ::deserialize( group , f_fixToMax , "FixToMaximum" ) )
  f_fixToMax = std::max( f_fixToMax , 0 );

 // variables for AC elements
 if( ::deserialize( group , "MaxReactivePower" , f_time_horizon ,
		    v_MaxReactivePower , true , true , v_change_intervals ) )
  if( std::all_of( v_MaxReactivePower.begin() , v_MaxReactivePower.end() ,
		   []( double i ) { return( i == 0 ); } ) )
   v_MaxReactivePower.clear();

 if( ::deserialize( group , "MinReactivePower" , f_time_horizon ,
		    v_MinReactivePower , true , true , v_change_intervals ) )
  if( std::all_of( v_MinReactivePower.begin() , v_MinReactivePower.end() ,
		   []( double i ) { return( i == 0 ); } ) )
   v_MinReactivePower.clear();

 // commitment-gated reactive coefficients (optional): the u[t] terms of the
 // state-dependent bound Qmin_off + Qmin_on u <= q <= Qmax_off + Qmax_on u;
 // absent or all-zero means the plain box above (backward-compatible)
 if( ::deserialize( group , "MaxReactivePowerOn" , f_time_horizon ,
		    v_MaxReactivePowerOn , true , true , v_change_intervals ) )
  if( std::all_of( v_MaxReactivePowerOn.begin() , v_MaxReactivePowerOn.end() ,
		   []( double i ) { return( i == 0 ); } ) )
   v_MaxReactivePowerOn.clear();

 if( ::deserialize( group , "MinReactivePowerOn" , f_time_horizon ,
		    v_MinReactivePowerOn , true , true , v_change_intervals ) )
  if( std::all_of( v_MinReactivePowerOn.begin() , v_MinReactivePowerOn.end() ,
		   []( double i ) { return( i == 0 ); } ) )
   v_MinReactivePowerOn.clear();

 // variables for the reference schedule
 ::deserialize( group , "ReferenceSchedule" , f_time_horizon ,
		v_RefSchedule , true , true , v_change_intervals );

 // a unit on at instant -1 produces there at least its minimum power
 // MinPower[ 0 ], whatever its availability at 0: a smaller InitialPower is
 // raised to it with a warning, before anything is derived from it, so that
 // all the formulations and the DP solvers see the same initial power, which
 // is also the one check_data_consistency() checks; an InitialPower above
 // MaxPower[ 0 ] is only warned about [see "InitialPower" in the .h]
 if( f_InitUpDownTime > 0 ) {
  if( f_InitialPower < v_MinPower.front() ) {
   std::cerr << "ThermalUnitBlock::deserialize: warning: initially-on unit "
                "(InitUpDownTime = " << f_InitUpDownTime << ") has InitialPower "
             << f_InitialPower << " < MinPower[ 0 ] = " << v_MinPower.front()
             << "; InitialPower is raised to MinPower[ 0 ]" << std::endl;
   f_InitialPower = v_MinPower.front();
   }
  else
   if( f_InitialPower > v_MaxPower.front() )
    std::cerr << "ThermalUnitBlock::deserialize: warning: initially-on unit "
                 "(InitUpDownTime = " << f_InitUpDownTime << ") has "
                 "InitialPower " << f_InitialPower << " > MaxPower[ 0 ] = "
              << v_MaxPower.front() << std::endl;
  }

 // start-up, shut-down limits: when not given, the operational minimum
 // power, which they follow when the availability changes [see
 // set_default_limits()]
 f_default_start_up_limit = ! ::deserialize( group , "StartUpLimit" ,
                                             f_time_horizon , v_StartUpLimit ,
                                             true , true ,
                                             v_change_intervals );
 if( f_default_start_up_limit )
  v_StartUpLimit.resize( f_time_horizon );

 f_default_shut_down_limit = ! ::deserialize( group , "ShutDownLimit" ,
                                              f_time_horizon ,
                                              v_ShutDownLimit , true , true ,
                                              v_change_intervals );
 if( f_default_shut_down_limit )
  v_ShutDownLimit.resize( f_time_horizon );

 set_default_limits();

 // the numbers of ramp steps of the SUSD formulation are computed by the
 // ThermalUnitBlock [see compute_ramp_steps()]: a datum has T entries, while
 // the formulation reads T + 1 of them
 if( ::deserialize( group , "MaxRampUpSteps" , f_time_horizon ,
                    v_MaxRampSteps , true , true , v_change_intervals ) )
  throw( std::invalid_argument( "ThermalUnitBlock::deserialize: "
                                "MaxRampUpSteps is computed by the "
                                "ThermalUnitBlock and "
                                "cannot be given" ) );

 if( ::deserialize( group , "MaxRampDownSteps" , f_time_horizon ,
                    v_MaxRampDownSteps , true , true , v_change_intervals ) )
  throw( std::invalid_argument( "ThermalUnitBlock::deserialize: "
                                "MaxRampDownSteps is computed by the "
                                "ThermalUnitBlock and "
                                "cannot be given" ) );

 // with InitialPower as it is after the clamp above
 compute_ramp_steps();

 check_data_consistency();

}  // end( ThermalUnitBlock::deserialize )

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG

/*
std::vector< std::string > ThermalUnitBlock::expected_dims( void )
 const {
 static const std::vector< std::string > ed = { };

 auto ret = UnitBlock::expected_dims();
 ret.insert( ret.end() , ed.begin() , ed.end() );

 return( ret );
 }

----------------------------------------------------------------------------*/

std::vector< std::string > ThermalUnitBlock::expected_vars( void )
 const {
 static const std::vector< std::string > ev =
  { "InvestmentCost" , "Capacity" , "Scale" , "MinPower" , "MaxPower" ,
    "DeltaRampUp" , "DeltaRampDown" , "PrimaryRho" , "SecondaryRho" ,
    "PrimarySpinningReserveCost" , "SecondarySpinningReserveCost" ,
    "ReactiveLinearTerm" ,
    "LinearTerm" , "QuadTerm" , "ConstTerm" , "StartUpCost" ,
    "ShutDownCost" ,
    "FixedConsumption" , "InertiaCommitment" , "InitialPower" , "MinUpTime" ,
    "MinDownTime" , "InitUpDownTime" , "Availability" , "StartUpLimit" ,
    "ShutDownLimit" , "MaxRampUpSteps" , "MaxRampDownSteps" ,
    "InitialReactivePower", "MaxReactivePower" , "MinReactivePower" ,
    "MaxReactivePowerOn" , "MinReactivePowerOn" ,
    "ReferenceSchedule" , "FixToMaximum"
    };

 auto ret = UnitBlock::expected_vars();
 ret.insert( ret.end() , ev.begin() , ev.end() );

 return( ret );
 }

#endif

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::check_data_consistency( void ) const
{
 // InvestmentCost- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ( f_InvestmentCost != 0 ) && ( f_InitUpDownTime >= 0 ) )
  throw( std::logic_error( "ThermalUnitBlock::check_data_consistency: the "
                           "presence of the investment cost of the thermal "
                           "allows the model to switch into the strategic "
                           "scenario mode, but the presence of a positive "
                           "initial up/down time, typical of the operative "
                           "scenario, is incompatible." ) );

 // Minimum and maximum power - - - - - - - - - - - - - - - - - - - - - - - -
 assert( v_MinPower.size() == f_time_horizon );
 assert( v_MaxPower.size() == f_time_horizon );

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {
  if( v_MinPower[ t ] > v_MaxPower[ t ] )
   throw( std::logic_error( "ThermalUnitBlock::check_data_consistency: "
                            "minimum power at time " + std::to_string( t ) +
                            " is " + std::to_string( v_MinPower[ t ] ) +
                            ", which is greater than the maximum power, which "
                            "is " + std::to_string( v_MaxPower[ t ] ) ) );

  if( v_MinPower[ t ] < 0 )
   throw( std::logic_error( "ThermalUnitBlock::check_data_consistency: "
                            "minimum power at time " + std::to_string( t ) +
			                         " is " + std::to_string( v_MinPower[ t ] ) +
                            ", but it must be nonnegative" ) );
  }

 // Availability- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_Availability.empty() ) {
  assert( v_Availability.size() == f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   if( ( v_Availability[ t ] < 0 ) || ( v_Availability[ t ] > 1 ) )
    throw( std::logic_error( "ThermalUnitBlock::check_data_consistency: "
                             "availability at time " + std::to_string( t ) +
			                          " is " + std::to_string( v_Availability[ t ] ) +
                             ", but it must be between 0 and 1" ) );
  }

 // Delta ramp-up - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_DeltaRampUp.empty() ) {
  assert( v_DeltaRampUp.size() == f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   if( v_DeltaRampUp[ t ] < 0 )
    throw( std::logic_error( "ThermalUnitBlock::check_data_consistency: "
                             "delta ramp-up at time " +
                             std::to_string( t ) + " is " +
                             std::to_string( v_DeltaRampUp[ t ] ) +
                             ", but it must be nonnegative" ) );
  }

 // Delta ramp-down - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_DeltaRampDown.empty() ) {
  assert( v_DeltaRampDown.size() == f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   if( v_DeltaRampDown[ t ] < 0 )
    throw( std::logic_error( "ThermalUnitBlock::check_data_consistency: "
                             "delta ramp-down at time " +
                             std::to_string( t ) + " is " +
                             std::to_string( v_DeltaRampDown[ t ] ) +
                             ", but it must be nonnegative" ) );
  }

 // QuadTerm - - - - - - - - - - - - - - - -- - - - - - - - - - - - - - - - -
 if( ! v_QuadTerm.empty() ) {
  assert( v_QuadTerm.size() == f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   if( v_QuadTerm[ t ] < 0 )
    throw( std::logic_error( "ThermalUnitBlock::check_data_consistency: "
                             "quadratic term at time " +
                             std::to_string( t ) + " is " +
                             std::to_string( v_QuadTerm[ t ] ) +
                             ", but it must be nonnegative" ) );
  }

 // FixedConsumption- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // UCBlock subtracts it from the injection of an off unit: a negative value
 // would be a production when off
 for( Index t = 0 ; t < v_FixedConsumption.size() ; ++t )
  if( v_FixedConsumption[ t ] < 0 )
   throw( std::logic_error( "ThermalUnitBlock::check_data_consistency: "
                            "fixed consumption at time " +
                            std::to_string( t ) + " is " +
                            std::to_string( v_FixedConsumption[ t ] ) +
                            ", but it must be nonnegative" ) );

 // InitialPower- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( f_InitialPower < 0 )
  throw( std::logic_error( "ThermalUnitBlock::check_data_consistency: "
                           "initial power is " +
                           std::to_string( f_InitialPower ) +
                           ", but it must be nonnegative" ) );

 // if the unit is initially committed, the initial power must lie within
 // the operational band of the unit: formulations and solvers (e.g. the
 // DP ones) rely on this documented contract, and violating it silently
 // produces spurious must-run behaviour rather than a clean error
 if( ( f_InitUpDownTime > 0 ) && ( f_InitialPower < v_MinPower.front() ) )
  throw( std::logic_error( "ThermalUnitBlock::check_data_consistency: the "
                           "unit is initially committed but the initial "
                           "power is " + std::to_string( f_InitialPower ) +
                           ", which is below the minimum power, which is " +
                           std::to_string( v_MinPower.front() ) ) );

 // StartUpLimit- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 for( Index t = 0 ; t < f_time_horizon ; ++t )
  if( ( v_StartUpLimit[ t ] < get_operational_min_power( t ) ) ||
      ( v_StartUpLimit[ t ] > get_operational_max_power( t ) ) )
   throw( std::logic_error( "ThermalUnitBlock::check_data_consistency: "
                            "start-up limit at time " +
                            std::to_string( t ) + " is " +
                            std::to_string( v_StartUpLimit[ t ] ) +
                            ", but it must be between " +
                            std::to_string( get_operational_min_power( t ) )
                            + " and " +
                            std::to_string( get_operational_max_power( t ) ) ) );

 // ShutDownLimit - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 for( Index t = 0 ; t < f_time_horizon ; ++t )
  if( ( v_ShutDownLimit[ t ] < get_operational_min_power( t ) ) ||
      ( v_ShutDownLimit[ t ] > get_operational_max_power( t ) ) )
   throw( std::logic_error( "ThermalUnitBlock::check_data_consistency: "
                            "shut-down limit at time " +
                            std::to_string( t ) + " is " +
                            std::to_string( v_ShutDownLimit[ t ] ) +
                            ", but it must be between " +
                            std::to_string( get_operational_min_power( t ) )
                            + " and " +
                            std::to_string( get_operational_max_power( t ) )
                            ) );

 }  // end( ThermalUnitBlock::check_data_consistency )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::generate_abstract_variables( Configuration * stvv )
{
 if( variables_generated() )  // variables have already been generated
  return;                     // nothing to do

 UnitBlock::generate_abstract_variables( stvv );

 // read Configuration to set the formulation - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
 Index wf = 1;  // T formulation
 if( ( ! stvv ) && f_BlockConfig )
  stvv = f_BlockConfig->f_static_variables_Configuration;
 if( auto sci = dynamic_cast< SimpleConfiguration< int > * >( stvv ) )
  wf = sci->f_value;

 init_t = first_free_instant( f_InitUpDownTime );

 // Design Binary Variable- - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( f_InvestmentCost != 0 ) {
  design.set_type( ColVariable::kBinary );
  add_static_variable( design , "x_thermal" );
  }
 else
  design.set_value( std::numeric_limits< double >::quiet_NaN() );

 // Commitment Variables- - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 v_commitment.resize( f_time_horizon );
 for( auto & var : v_commitment )
  var.set_type( ColVariable::kBinary );
 add_static_variable( v_commitment , "u_thermal" );

 // Active Power Variables- - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 v_active_power.resize( f_time_horizon );
 for( auto & var : v_active_power )
  var.set_type( ColVariable::kNonNegative );
 add_static_variable( v_active_power , "p_thermal" );

 // Reactive Power Variables, if any- - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( f_reactive_power ) {
  v_reactive_power.resize( f_time_horizon );
  for( auto & var : v_reactive_power )  // absorbed if negative, the bounds
   var.set_type( ColVariable::kContinuous );  // are rows
  add_static_variable( v_reactive_power , "q_thermal" );
  }

 // Start-Up and Shut-Down Binary Variables - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 auto startup_shutdown_size = f_time_horizon - init_t;

 if( startup_shutdown_size > 0 ) {
  auto vartype = ( wf & ZWCont ) ? ColVariable::kPosUnitary
                                 : ColVariable::kBinary;

  v_start_up.resize( startup_shutdown_size );
  for( auto & var : v_start_up )
   var.set_type( vartype );
  add_static_variable( v_start_up , "v_thermal" );

  v_shut_down.resize( startup_shutdown_size );
  for( auto & var : v_shut_down )
   var.set_type( vartype );
  add_static_variable( v_shut_down , "w_thermal" );
  }

 // Primary Spinning Reserve Variables- - - - - - - - - - - - - - - - - - - -
 if( reserve_vars & 1u )  // if UCBlock has primary demand variables
  if( ! v_PrimaryRho.empty() ) {  // if unit produces any primary reserve
   v_primary_spinning_reserve.resize( f_time_horizon );
   for( auto & var : v_primary_spinning_reserve )
    var.set_type( ColVariable::kNonNegative );
   add_static_variable( v_primary_spinning_reserve , "pr_thermal" );
   }

 // Secondary Spinning Reserve Variables- - - - - - - - - - - - - - - - - - -
 if( reserve_vars & 2u )  // if UCBlock has secondary demand variables
  if( ! v_SecondaryRho.empty() ) {  // if unit produces any secondary reserve
   v_secondary_spinning_reserve.resize( f_time_horizon );
   for( auto & var : v_secondary_spinning_reserve )
    var.set_type( ColVariable::kNonNegative );
   add_static_variable( v_secondary_spinning_reserve , "sc_thermal" );
   }

 // possibly fix the commitment variables to 0 or 1 - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( f_InitUpDownTime > 0 ) {  // InitUpDownTime > 0- - - - - - - - - - - - -
  for( Index t = 0 ; t < init_t ; ++t ) {
   if( ! v_commitment.empty() ) {
    v_commitment[ t ].set_value( 1 );
    v_commitment[ t ].is_fixed( true , eNoMod );
    }
   }

  for( Index t = init_t ;
       t < std::min( init_t + f_MinDownTime , f_time_horizon ) ; ++t ) {
   v_start_up[ t - init_t ].set_value( 0 );
   v_start_up[ t - init_t ].is_fixed( true , eNoMod );
   }
  }
 else {  // InitUpDownTime <= 0 - - - - - - - - - - - - - - - - - - - - - - -
  for( Index t = 0 ; t < init_t ; ++t ) {
   if( ! v_active_power.empty() ) {
    v_active_power[ t ].set_value( 0 );
    v_active_power[ t ].is_fixed( true , eNoMod );
    }

   if( ! v_commitment.empty() ) {
    v_commitment[ t ].set_value( 0 );
    v_commitment[ t ].is_fixed( true , eNoMod );
    }

   if( ! v_primary_spinning_reserve.empty() ) {
    v_primary_spinning_reserve[ t ].set_value( 0 );
    v_primary_spinning_reserve[ t ].is_fixed( true , eNoMod );
    }

   if( ! v_secondary_spinning_reserve.empty() ) {
    v_secondary_spinning_reserve[ t ].set_value( 0 );
    v_secondary_spinning_reserve[ t ].is_fixed( true , eNoMod );
    }
   }

  for( Index t = init_t ;
       t < std::min( init_t + f_MinUpTime , f_time_horizon ) ; ++t ) {
   v_shut_down[ t - init_t ].set_value( 0 );
   v_shut_down[ t - init_t ].is_fixed( true , eNoMod );
   }
  }

 bool f_cuts = wf & PCuts;

 // Prospective Cuts Variables- - - - - - - - - - - - - - - - - - - - - - - -

 if( f_cuts ) {
  AR |= PCuts;

  prevpbar.resize( f_time_horizon , 0 );
  v_cut.resize( f_time_horizon );
  for( auto & var : v_cut )
   var.set_type( ColVariable::kNonNegative );
  add_static_variable( v_cut , "z_thermal" );
  }

 // Active Power & Prospective Cuts Variables for DP, SU and SD formulations-
 // (with auxiliary structures initialization)- - - - - - - - - - - - - - - -

 switch( wf & FormMsk ) {
  case( tbinForm ):  // 3bin formulation- - - - - - - - - - - - - - - - - - -
   // AR |= tbinForm;  // does nothing
   break;

  case( TForm ):  // T formulation- - - - - - - - - - - - - - - - - - - - - -
   AR |= TForm;
   break;

  case( ptForm ):  // pt formulation- - - - - - - - - - - - - - - - - - - - -
   // fall through
  case( DPForm ):  // DP formulation- - - - - - - - - - - - - - - - - - - - -
   // fall through
  case( SUForm ):  // SU formulation- - - - - - - - - - - - - - - - - - - - -
   // fall through
  case( SDForm ):  // SD formulation- - - - - - - - - - - - - - - - - - - - -

  case( SUSDForm ):  // SUSD formulation- - - - - - - - - - - - - - - - - - -

   if( f_InitUpDownTime > 0 ) {  // if initial committed, OFF_0

    // intervals ( 0 , k ) represent the pre-horizon ON run terminating (i.e.,
    // the unit shutting down) at k; k == init_t == 0 is the degenerate ( 0 , 0
    // ) interval covering no in-horizon period, i.e., the shut-down at
    // t == 0 (possible only when init_t == 0, the minimum up time being
    // already met). Since ( 0 , 0 ) covers no period, no ramp or shut-down-
    // limit row gates it: the shut-down at 0 is allowed only if
    // InitialPower <= ShutDownLimit[ 0 ], whatever the ramps, as in every
    // formulation and in the DP solvers, and the interval is fixed to 0 by
    // a row otherwise [see "ShutDownZero_Const_Thermal" in build_rows()],
    // so that the Variable do not depend on InitialPower nor on the limits
    for( Index k = init_t ; k <= f_time_horizon + 1 ; ++k )
     v_Y_plus.push_back( std::make_pair( 0 , k ) );

    for( Index k = init_t ; k <= f_time_horizon ; ++k ) {
     if( k <= f_time_horizon - f_MinDownTime - 1 )
      for( Index h = ( k + f_MinDownTime + 1 ) ;  // OFF_h
           h <= f_time_horizon ; ++h )
       v_Y_minus.push_back( std::make_pair( k , h ) );
     v_Y_minus.push_back( std::make_pair( k , f_time_horizon + 1 ) );
     v_nodes_minus.push_back( k );
    }

    for( Index h = ( init_t + f_MinDownTime + 1 ) ;  // OFF_h
         h <= f_time_horizon ; ++h ) {
     if( h <= f_time_horizon - f_MinUpTime + 1 )
      for( Index k = ( h + f_MinUpTime - 1 ) ;  // ON_k
           k <= f_time_horizon ; ++k )
       v_Y_plus.push_back( std::make_pair( h , k ) );
     v_Y_plus.push_back( std::make_pair( h , f_time_horizon + 1 ) );
     v_nodes_plus.push_back( h );
    }

   } else {

    for( Index h = init_t + 1 ; h <= f_time_horizon + 1 ; ++h )  // ON_k
     v_Y_minus.push_back( std::make_pair( 0 , h ) );

    for( Index h = init_t + 1 ; h <= f_time_horizon ; ++h ) {  // OFF_h
     if( h <= f_time_horizon - f_MinUpTime + 1 )
      for( Index k = ( h + f_MinUpTime - 1 ) ;  // ON_k
           k <= f_time_horizon ; ++k )
       v_Y_plus.push_back( std::make_pair( h , k ) );
     v_Y_plus.push_back( std::make_pair( h , f_time_horizon + 1 ) );
     v_nodes_plus.push_back( h );
    }

    for( Index k = init_t + f_MinUpTime ;
         k <= f_time_horizon ; ++k ) {  // ON_k
     if( k <= f_time_horizon - f_MinDownTime - 1 )
      for( Index h = ( k + f_MinDownTime + 1 ) ;  // OFF_h
           h <= f_time_horizon ; ++h )
       v_Y_minus.push_back( std::make_pair( k , h ) );
     v_Y_minus.push_back( std::make_pair( k , f_time_horizon + 1 ) );
     v_nodes_minus.push_back( k );
    }
   }

   v_commitment_plus.resize( v_Y_plus.size() );
   for( auto & var : v_commitment_plus )
    var.set_type( ColVariable::kBinary );
   add_static_variable( v_commitment_plus , "y_plus_thermal" );

   v_commitment_minus.resize( v_Y_minus.size() );
   for( auto & var : v_commitment_minus )
    var.set_type( ColVariable::kBinary );
   add_static_variable( v_commitment_minus , "y_minus_thermal" );

   for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
    for( Index t = 0 ; t < f_time_horizon ; ++t )
     if( ( v_Y_plus[ i ].first <= t + 1 ) &&
	 ( t + 1 <= v_Y_plus[ i ].second ) )
      v_P_h_k.push_back(
	      std::make_pair( t , std::make_pair( v_Y_plus[ i ].first ,
						  v_Y_plus[ i ].second ) ) );

   if( ( wf & FormMsk ) == ptForm )  // pt formulation- - - - - - - - - - - -
    AR |= ptForm;

   if( ( wf & FormMsk ) == DPForm ) {  // DP formulation- - - - - - - - - - -
    AR |= DPForm;

    v_active_power_h_k.resize( v_P_h_k.size() );
    for( auto & var : v_active_power_h_k )
     var.set_type( ColVariable::kNonNegative );
    add_static_variable( v_active_power_h_k , "p_h_k_thermal" );

    if( f_cuts ) {

     for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
      for( Index t = 0 ; t < f_time_horizon ; ++t )
       if( ( v_Y_plus[ i ].first <= t + 1 ) &&
           ( t + 1 <= v_Y_plus[ i ].second ) )
        v_Z_h_k.push_back(
         std::make_pair( t , std::make_pair( v_Y_plus[ i ].first ,
                                             v_Y_plus[ i ].second ) ) );

     v_cut_h_k.resize( v_Z_h_k.size() );
     for( auto & var : v_cut_h_k )
      var.set_type( ColVariable::kNonNegative );
     add_static_variable( v_cut_h_k , "z_h_k_thermal" );
    }
   }

   if( ( wf & FormMsk ) == SUForm || ( wf & FormMsk ) == SUSDForm )
   {  // SU or SUSD formulation - - - - - - - - - - - - - - - - - - - - - - -

     if( ( wf & FormMsk ) == SUForm )
       AR |= SUForm;
     if( ( wf & FormMsk ) == SUSDForm )
       AR |= SUSDForm;

     for( Index i = 0 ; i < v_Y_plus.size() ; ++i ) {
     bool check_var = true;
     if( i > 0 )
      for( Index j = 0 ; j < v_P_h.size() ; ++j )
       if( v_P_h[ j ].second == v_Y_plus[ i ].first )
        check_var = false;
     if( check_var )
      for( Index t = 0 ; t < f_time_horizon ; ++t )
       if( v_Y_plus[ i ].first <= t + 1 ) {
        v_P_h.push_back( std::make_pair( t , v_Y_plus[ i ].first ) );

        if( f_cuts )
         v_Z_h.push_back( std::make_pair( t , v_Y_plus[ i ].first ) );
       }
    }

    v_active_power_h.resize( v_P_h.size() );
    for( auto & var : v_active_power_h )
     var.set_type( ColVariable::kNonNegative );
    add_static_variable( v_active_power_h , "p_h_thermal" );

    if( f_cuts ) {
     v_cut_h.resize( v_Z_h.size() );
     for( auto & var : v_cut_h )
      var.set_type( ColVariable::kNonNegative );
     add_static_variable( v_cut_h , "z_h_thermal" );

     // add variables for SUSD formulation maximizing perspective cuts
     // variables of SU and SD formulations - - - - - - - - - - - - - - - - -
     if( ( wf & FormMsk ) == SUSDForm ) {
      v_cut_teta.resize( f_time_horizon );
      for( auto & var : v_cut_teta )
       var.set_type( ColVariable::kNonNegative );
      add_static_variable( v_cut_teta , "teta_thermal" );
     }
    }
   }

   if( ( wf & FormMsk ) == SDForm || ( wf & FormMsk ) == SUSDForm )
   {  // SD or SUSD formulation - - - - - - - - - - - - - - - - - - - - - - -

    if( ( wf & FormMsk ) == SDForm )
     AR |= SDForm;
    if( ( wf & FormMsk ) == SUSDForm )
     AR |= SUSDForm;
    for( Index i = 0 ; i < v_Y_plus.size() ; ++i ) {
     bool check_var = true;
     if( i > 0 )
      for( Index j = 0 ; j < v_P_k.size() ; ++j )
       if( v_P_k[ j ].second == v_Y_plus[ i ].second )
        check_var = false;
     if( check_var )
      for( Index t = 0 ; t < f_time_horizon ; ++t )
       if( v_Y_plus[ i ].second >= t + 1 ) {
        v_P_k.push_back( std::make_pair( t , v_Y_plus[ i ].second ) );

        if( f_cuts )
         v_Z_k.push_back( std::make_pair( t , v_Y_plus[ i ].second ) );
       }
    }

    v_active_power_k.resize( v_P_k.size() );
    for( auto & var : v_active_power_k )
     var.set_type( ColVariable::kNonNegative );
    add_static_variable( v_active_power_k , "p_k_thermal" );

    if( f_cuts ) {

     v_cut_k.resize( v_Z_k.size() );
     for( auto & var : v_cut_k )
      var.set_type( ColVariable::kNonNegative );
     add_static_variable( v_cut_k , "z_k_thermal" );
    }
   }

   break;

  default:
   throw( std::invalid_argument(
    "ThermalUnitBlock::generate_abstract_variables: invalid formulation" ) );

  }  // end( switch )

 // the variables for the reference schedule, if there- - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! v_RefSchedule.empty() ) {
  v_abs_ref_schedule.resize( f_time_horizon );
  for( auto & var : v_abs_ref_schedule )
   var.set_type( ColVariable::kNonNegative );
  add_static_variable( v_abs_ref_schedule , "v_abs_refschd" );
  }

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 set_variables_generated();

 }  // end( ThermalUnitBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::generate_abstract_constraints( Configuration * stcc )
{
 if( constraints_generated() )  // constraints have already been generated
  return;                       // nothing to do

 bool generate_ZOConstraints = false;
 if( ( ! stcc ) && f_BlockConfig )
  stcc = f_BlockConfig->f_static_constraints_Configuration;
 if( auto sci = dynamic_cast< SimpleConfiguration< int > * >( stcc ) )
  generate_ZOConstraints = sci->f_value;

 build_rows( generate_ZOConstraints );

 set_constraints_generated();

 }  // end( ThermalUnitBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::build_rows( bool generate_ZOConstraints )
{
 // when comparing, the ZOConstraints are those generated, which nothing
 // changes
 if( ! generating_rows() )
  generate_ZOConstraints = false;

 LinearFunction::v_coeff_pair vars;

 // commitment design binary variable constraints - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( f_InvestmentCost != 0 ) {
  size_rows( CommitmentDesign_Const , f_time_horizon );

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   vars.push_back( std::make_pair( & v_commitment[ t ] , 1 ) );
   vars.push_back( std::make_pair( & design , -1 ) );

   put_row( CommitmentDesign_Const , t , std::move( vars ) ,
            -Inf< double >() , 0 );
   }

  add_rows( CommitmentDesign_Const ,
                         "CommitmentDesign_Const_Thermal" );
  }

 // the cap on the last on-power p_t of an on-interval that closes at t: the
 // ShutDownLimit of the shut-down instant t + 1, none (the operational
 // maximum power) if t is the last period, as in the 3bin and T formulations
 // and in the DP solvers
 auto sd_cap = [ this ]( Index t ) -> double {
  const double pmax = get_operational_max_power( t );
  return( t + 1 < f_time_horizon ?
          std::min( double( v_ShutDownLimit[ t + 1 ] ) , pmax ) : pmax );
  };

 // compute the \psi constants for DP-related formulations- - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 std::vector< double > v_psi;

 if( ( ( AR & FormMsk ) == ptForm ) || ( ( AR & FormMsk ) == SUForm ) ||
     ( ( AR & FormMsk ) == SDForm ) || ( ( AR & FormMsk ) == SUSDForm ) ) {
  v_psi.resize( v_P_h_k.size() );

  for( Index j = 0 ; j < v_P_h_k.size() ; ++j ) {
   const Index t = v_P_h_k[ j ].first;
   const Index h = v_P_h_k[ j ].second.first;
   const Index k = v_P_h_k[ j ].second.second;
   // the interval ( h , k ) is on in the periods h - 1 , ... , k - 1: it
   // starts up at h - 1 (before the horizon if h == 0) and shuts down at k
   // (beyond the horizon if k >= T); every instant of it, the first and
   // the last included, is capped by the maximum power, by the start-up
   // limit plus the ramps up from the start-up (the initial power plus the
   // ramps from -1 for h == 0) and by the shut-down limit plus the ramps
   // down to the shut-down, the limits alone if the ramps are not given
   v_psi[ j ] = get_operational_max_power( t );

   if( k < f_time_horizon ) {
    if( t + 1 == k )
     v_psi[ j ] = std::min( v_psi[ j ] , double( v_ShutDownLimit[ k ] ) );
    else
     if( ! v_DeltaRampDown.empty() )
      v_psi[ j ] = std::min( v_psi[ j ] , v_ShutDownLimit[ k ] +
                                          ramp_down_sum( t + 1 , k - 1 ) );
    }

   if( h > 0 ) {
    if( t + 1 == h )
     v_psi[ j ] = std::min( v_psi[ j ] , double( v_StartUpLimit[ t ] ) );
    else
     if( ! v_DeltaRampUp.empty() )
      v_psi[ j ] = std::min( v_psi[ j ] , v_StartUpLimit[ h - 1 ] +
                                          ramp_up_sum( h , t ) );
    }
   else
    if( ( f_InitUpDownTime > 0 ) && ( ! v_DeltaRampUp.empty() ) )
     v_psi[ j ] = std::min( v_psi[ j ] ,
                            f_InitialPower + ramp_up_sum( 0 , t ) );
   // a cap below the minimum power means that the run cannot be on at t,
   // which the maximum and minimum power rows then say together
   }
  }

 // the positions in v_P_h_k (hence in v_psi) of the instant t of each run
 // ( h , k ), in increasing order, so that the rows below find them without
 // a scan of v_P_h_k
 std::map< std::tuple< Index , Index , Index > , std::vector< Index > >
  psi_pos;
 for( Index j = 0 ; j < v_psi.size() ; ++j )
  psi_pos[ std::make_tuple( v_P_h_k[ j ].first ,
                            v_P_h_k[ j ].second.first ,
                            v_P_h_k[ j ].second.second ) ].push_back( j );
 const std::vector< Index > no_pos;
 auto psi_index = [ & ]( Index t , Index i ) -> const std::vector< Index > & {
  const auto it = psi_pos.find( std::make_tuple( t , v_Y_plus[ i ].first ,
                                                 v_Y_plus[ i ].second ) );
  return( it == psi_pos.end() ? no_pos : it->second );
  };

 // add the "fixed to maximum generation" constraints - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( f_fixToMax > 0 ) {
  size_rows( fixed_to_max_Power_Const , f_time_horizon );

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   // P_t >= Pmax(t)
   vars.push_back( std::make_pair( & v_active_power[ t ] , 1.0 ) );
   put_row( fixed_to_max_Power_Const , t , std::move( vars ) ,
            get_operational_max_power( t ) , Inf< double >() );
   }
  add_rows( fixed_to_max_Power_Const, "FixedGeneration" );
  }

 switch( AR & FormMsk ) {
  case( tbinForm ):  // 3bin formulation- - - - - - - - - - - - - - - - - - -
   // fall through
  case( TForm ): {  // T formulation- - - - - - - - - - - - - - - - - - - - -

   // Initializing start-up and shut-down variables connection constraints- -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   auto startup_shutdown_const_size = f_time_horizon - init_t;

   if( startup_shutdown_const_size > 0 ) {

    size_rows( StartUp_ShutDown_Variables_Const ,
               startup_shutdown_const_size );

    for( Index t = init_t , cnstr_idx = 0 ; t < f_time_horizon ;
         ++t , ++cnstr_idx ) {

     vars.push_back( std::make_pair( &v_commitment[ t ] , 1.0 ) );
     vars.push_back( std::make_pair( &v_start_up[ t - init_t ] , -1.0 ) );
     vars.push_back( std::make_pair( &v_shut_down[ t - init_t ] , 1.0 ) );

     if( t > init_t )
      vars.push_back( std::make_pair( &v_commitment[ t - 1 ] , -1.0 ) );

     const double both = ( ( t == init_t ) && ( f_InitUpDownTime > 0 ) ) ?
                         1.0 : 0.0;
     put_row( StartUp_ShutDown_Variables_Const , cnstr_idx ,
              std::move( vars ) , both , both );
    }

    add_rows( StartUp_ShutDown_Variables_Const ,
                           "StartUp_ShutDown_Variables_Const_Thermal" );
   }

   // Initializing turn on constraints (start-up constraints) - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   // the rows start at init_t + MinUpTime - 1, the first instant whose
   // window of MinUpTime instants lies in [ init_t , T - 1 ]; when the
   // horizon ends before it, the only row is the one at T - 1 with the
   // window cut to [ init_t , T - 1 ], since otherwise nothing forbids a
   // run shorter than MinUpTime that ends within the horizon (when a full
   // row exists, it implies the cut ones for binary variables)

   if( f_time_horizon > init_t ) {

    const Index su_first = std::min( init_t + f_MinUpTime - 1 ,
                                     f_time_horizon - 1 );
    size_rows( StartUp_Const , f_time_horizon - su_first );

    for( Index t = su_first , cnstr_idx = 0 ;
        t < f_time_horizon ; ++t, ++cnstr_idx ) {

     for( Index s = ( t + 1 >= init_t + f_MinUpTime ?
                      t + 1 - init_t - f_MinUpTime : 0 ) ;
          s <= t - init_t ; ++s )
      vars.push_back( std::make_pair( &v_start_up[ s ] , -1.0 ) );

     vars.push_back( std::make_pair( &v_commitment[ t ] , 1.0 ) );

     put_row( StartUp_Const , cnstr_idx , std::move( vars ) ,
              0.0 , Inf< double >() );
    }

    add_rows( StartUp_Const ,
                           "StartUp_Commitment_Const_Thermal" );
   }

   // Initializing turn off constraints (shut-down constraints) - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   // as for the start-up rows, the row at T - 1 with the window cut to
   // [ init_t , T - 1 ] when no full window fits in the horizon

   if( f_time_horizon > init_t ) {

    const Index sd_first = std::min( init_t + f_MinDownTime - 1 ,
                                     f_time_horizon - 1 );
    size_rows( ShutDown_Const , f_time_horizon - sd_first );

    for( Index t = sd_first , cnstr_idx = 0 ;
	 t < f_time_horizon ; ++t, ++cnstr_idx ) {
     for( Index s = ( t + 1 >= init_t + f_MinDownTime ?
                      t + 1 - init_t - f_MinDownTime : 0 ) ;
          s <= t - init_t ; ++s )
      vars.push_back( std::make_pair( &v_shut_down[ s ] , 1.0 ) );

     vars.push_back( std::make_pair( &v_commitment[ t ] , 1.0 ) );

     put_row( ShutDown_Const , cnstr_idx , std::move( vars ) ,
              -Inf< double >() , 1.0 );
    }

    add_rows( ShutDown_Const ,
                           "ShutDown_Commitment_Const_Thermal" );
   }

   break;
  }

  case( ptForm ):  // pt formulation- - - - - - - - - - - - - - - - - - - - -
   // fall through
  case( DPForm ):  // DP formulation- - - - - - - - - - - - - - - - - - - - -
   // fall through
  case( SUForm ):  // SU formulation- - - - - - - - - - - - - - - - - - - - -
   // fall through
  case( SDForm ):  // SD formulation- - - - - - - - - - - - - - - - - - - - -
   // fall through
  case ( SUSDForm ):  {  // SUSD formulation- - - - - - - - - - - - - - - - -

   // Constraints connecting power variables of 3bin with those of DP, SU and
   // SD formulations - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if( ( AR & FormMsk ) == ptForm ) {  // pt formulation- - - - - - - - - - -
    ;  // does nothing
   } else {

    if( ( AR & FormMsk ) == SUSDForm )
     size_rows( Eq_ActivePower_Const , 2 * f_time_horizon );
    else
     size_rows( Eq_ActivePower_Const , f_time_horizon );

    if( ( AR & FormMsk ) == DPForm ) {  // DP formulation - - - - - - - - - -

     for( Index t = 0 ; t < f_time_horizon ; ++t ) {

      vars.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );

      for( Index j = 0 ; j < v_P_h_k.size() ; ++j )
       if( v_P_h_k[ j ].first == t )
        vars.push_back( std::make_pair( &v_active_power_h_k[ j ] , -1.0 ) );

      put_row( Eq_ActivePower_Const , t , std::move( vars ) , 0.0 , 0.0 );
     }
    }

    if( ( AR & FormMsk ) == SUForm || ( AR & FormMsk ) == SUSDForm )
    {  // SU or SUSD formulation- - - - - - - - - - - - - - - - - - - - - - -

     for( Index t = 0 ; t < f_time_horizon ; ++t ) {

      vars.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );

      for( Index j = 0 ; j < v_P_h.size() ; ++j )
       if( v_P_h[ j ].first == t )
        vars.push_back( std::make_pair( &v_active_power_h[ j ] , -1.0 ) );

      put_row( Eq_ActivePower_Const , t , std::move( vars ) , 0.0 , 0.0 );
     }
    }

    if( ( AR & FormMsk ) == SDForm || ( AR & FormMsk ) == SUSDForm )
    {  // SD or SUSD formulation- - - - - - - - - - - - - - - - - - - - - - -

     for( Index t = 0 ; t < f_time_horizon ; ++t ) {

      vars.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );

      for( Index j = 0 ; j < v_P_k.size() ; ++j )
       if( v_P_k[ j ].first == t )
        vars.push_back( std::make_pair( &v_active_power_k[ j ] , -1.0 ) );

      if( ( AR & FormMsk ) == SUSDForm ) {
       put_row( Eq_ActivePower_Const , f_time_horizon + t , std::move( vars ) ,
                0.0 , 0.0 );
      } else {
       put_row( Eq_ActivePower_Const , t , std::move( vars ) , 0.0 , 0.0 );
      }
     }
    }

    add_rows( Eq_ActivePower_Const ,
                           "Eq_ActivePower_Const_Thermal" );
   }

   // Constraints connecting commitment variables of 3bin with those of pt,
   // DP, SU and SD formulations- - - - - - - - - - - - - - - - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   size_rows( Eq_Commitment_Const , f_time_horizon );

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {

    vars.push_back( std::make_pair( &v_commitment[ t ] , 1.0 ) );

    for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
     if( ( v_Y_plus[ i ].first <= t + 1 ) && ( t + 1 <= v_Y_plus[ i ].second ) )
      vars.push_back( std::make_pair( &v_commitment_plus[ i ] , -1.0 ) );

    put_row( Eq_Commitment_Const , t , std::move( vars ) , 0.0 , 0.0 );
   }

   add_rows( Eq_Commitment_Const ,
                          "Eq_Commitment_Const_Thermal" );

   // Constraints connecting start-up variables of 3bin with those of pt,
   // DP, SU and SD formulations- - - - - - - - - - - - - - - - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   size_rows( Eq_StartUp_Const , f_time_horizon - init_t );

   for( Index t = init_t , cnstr_idx = 0 ; t < f_time_horizon ;
       ++t, ++cnstr_idx ) {

    vars.push_back( std::make_pair( &v_start_up[ t - init_t ] , 1.0 ) );

    for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
     if( ( v_Y_plus[ i ].first == t + 1 ) && ( t + 1 <= v_Y_plus[ i ].second ) )
      vars.push_back( std::make_pair( &v_commitment_plus[ i ] , -1.0 ) );

    put_row( Eq_StartUp_Const , cnstr_idx , std::move( vars ) , 0.0 , 0.0 );
   }

   add_rows( Eq_StartUp_Const , "Eq_StartUp_Const_Thermal" );

   // Constraints connecting shut-down variables of 3bin with those of pt,
   // DP, SU and SD formulations- - - - - - - - - - - - - - - - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   size_rows( Eq_ShutDown_Const , f_time_horizon - init_t );

   for( Index t = init_t , cnstr_idx = 0 ; t < f_time_horizon ;
       ++t, ++cnstr_idx ) {

    vars.push_back( std::make_pair( &v_shut_down[ t - init_t ] , 1.0 ) );

    for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
     if( ( v_Y_plus[ i ].first <= t ) && ( t == v_Y_plus[ i ].second ) )
      vars.push_back( std::make_pair( &v_commitment_plus[ i ] , -1.0 ) );

    put_row( Eq_ShutDown_Const , cnstr_idx , std::move( vars ) , 0.0 , 0.0 );
   }

   add_rows( Eq_ShutDown_Const ,
                          "Eq_ShutDown_Const_Thermal" );

   // Network Constraints - - - - - - - - - - - - - - - - - - - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   size_rows( Network_Const , v_nodes_plus.size() + v_nodes_minus.size() + 2 );

   auto cnstr_idx = 0;

  for( Index t = 0 ; t < 1 ; ++t, ++cnstr_idx ) {

    if( f_InitUpDownTime > 0 )
     for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
      if( v_Y_plus[ i ].first == t )
       vars.push_back( std::make_pair( &v_commitment_plus[ i ] , -1.0 ) );

    if( f_InitUpDownTime <= 0 )
     for( Index i = 0 ; i < v_Y_minus.size() ; ++i )
      if( v_Y_minus[ i ].first == t )
       vars.push_back( std::make_pair( &v_commitment_minus[ i ] , -1.0 ) );

    put_row( Network_Const , cnstr_idx , std::move( vars ) , -1.0 , -1.0 );
   }

  for( Index t = 0 ; t < v_nodes_plus.size() ; ++t, ++cnstr_idx ) {

    for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
     if( v_Y_plus[ i ].first == v_nodes_plus[ t ] )
      vars.push_back( std::make_pair( &v_commitment_plus[ i ] , -1.0 ) );

    for( Index i = 0 ; i < v_Y_minus.size() ; ++i )
     if( v_Y_minus[ i ].second == v_nodes_plus[ t ] )
      vars.push_back( std::make_pair( &v_commitment_minus[ i ] , 1.0 ) );

    put_row( Network_Const , cnstr_idx , std::move( vars ) , 0.0 , 0.0 );
   }

   for( Index t = 0 ; t < v_nodes_minus.size() ; ++t , ++cnstr_idx ) {

    for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
     if( v_Y_plus[ i ].second == v_nodes_minus[ t ] )
      vars.push_back( std::make_pair( &v_commitment_plus[ i ] , 1.0 ) );

    for( Index i = 0 ; i < v_Y_minus.size() ; ++i )
     if( v_Y_minus[ i ].first == v_nodes_minus[ t ] )
      vars.push_back( std::make_pair( &v_commitment_minus[ i ] , -1.0 ) );

    put_row( Network_Const , cnstr_idx , std::move( vars ) , 0.0 , 0.0 );
   }

   for( Index t = f_time_horizon + 1 ;
        t < f_time_horizon + 2 ; ++t , ++cnstr_idx ) {

    for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
     if( v_Y_plus[ i ].second == t )
      vars.push_back( std::make_pair( &v_commitment_plus[ i ] , 1.0 ) );

    for( Index i = 0 ; i < v_Y_minus.size() ; ++i )
     if( v_Y_minus[ i ].second == t )
      vars.push_back( std::make_pair( &v_commitment_minus[ i ] , 1.0 ) );

    put_row( Network_Const , cnstr_idx , std::move( vars ) , 1.0 , 1.0 );
   }

   add_rows( Network_Const , "Network_Const_Thermal" );

   break;
  }

  default:

   exit( 1 );

 }  // end( switch )

 // the shut-down at 0- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // the output of the instant before a shut-down is at most the shut-down
 // limit of the shut-down instant: a unit on before the horizon, which may
 // shut down at 0 (init_t == 0), can therefore do it only if InitialPower
 // <= ShutDownLimit[ 0 ], whatever the ramps. In the 3bin and T formulations
 // the ramp-down row at 0 says so if DeltaRampDown is given, and the row
 // w_0 <= 0 or 1 is there otherwise; in the other formulations the row says
 // the same of the interval ( 0 , 0 ) [see generate_abstract_variables()].
 // Its right-hand side follows InitialPower and ShutDownLimit[ 0 ]
 if( ( f_InitUpDownTime > 0 ) && ( init_t == 0 ) && ( f_time_horizon > 0 ) ) {
  ColVariable * w0 = nullptr;
  if( ( ( AR & FormMsk ) == tbinForm ) || ( ( AR & FormMsk ) == TForm ) ) {
   if( v_DeltaRampDown.empty() )
    w0 = & v_shut_down[ 0 ];
   }
  else
   for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
    if( ( v_Y_plus[ i ].first == 0 ) && ( v_Y_plus[ i ].second == 0 ) )
     w0 = & v_commitment_plus[ i ];

  if( w0 ) {
   size_rows( ShutDownZero_Const , 1 );
   vars.push_back( std::make_pair( w0 , 1.0 ) );
   put_row( ShutDownZero_Const , 0 , std::move( vars ) , -Inf< double >() ,
            f_InitialPower > v_ShutDownLimit[ 0 ] ? 0.0 : 1.0 );
   add_rows( ShutDownZero_Const , "ShutDownZero_Const_Thermal" );
   }
  }

 // Initializing ramp-up constraints- - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( ! v_DeltaRampUp.empty() ) && ( f_InitUpDownTime > 0 ) &&
     ( f_InitialPower + v_DeltaRampUp[ 0 ] < get_operational_min_power( 0 ) ) )
  throw( std::logic_error(
   "ThermalUnitBlock::generate_abstract_constraints: when InitUpDownTime > 0 "
   "it must be InitialPower + DeltaRampUp[ 0 ] >= the operational minimum "
   "power at 0" ) );

 if( ( ! v_DeltaRampDown.empty() ) && ( f_InitUpDownTime > 0 ) &&
     ( f_InitialPower - v_DeltaRampDown[ 0 ] >
       get_operational_max_power( 0 ) ) )
  throw( std::logic_error(
   "ThermalUnitBlock::generate_abstract_constraints: when InitUpDownTime > 0 "
   "it must be InitialPower - DeltaRampDown[ 0 ] <= the operational maximum "
   "power at 0" ) );

 if( ! v_DeltaRampUp.empty() ) {

  if( ( AR & FormMsk ) == tbinForm ) {  // 3bin formulation - - - - - - - - -

   auto ramp_up_cnstrs_size = f_time_horizon;
   size_rows( RampUp_Const , ramp_up_cnstrs_size );
   auto cnstr_idx = 0;

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {

    vars.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );

    if( t > 0 ) {
     vars.push_back( std::make_pair( &v_active_power[ t - 1 ] , -1.0 ) );
     vars.push_back( std::make_pair( &v_commitment[ t - 1 ] ,
                                     -v_DeltaRampUp[ t ] ) );
    }

    if( t >= init_t )
     vars.push_back( std::make_pair( &v_start_up[ t - init_t ] ,
                                     -v_StartUpLimit[ t ] ) );

    put_row( RampUp_Const , cnstr_idx++ , std::move( vars ) ,
             -Inf< double >() ,
             ( ( t == 0 ) && ( f_InitUpDownTime > 0 ) ) ?
             f_InitialPower + v_DeltaRampUp[ t ] : 0.0 );
   }

  }
  if( ( AR & FormMsk ) == TForm ) {  // T formulation - - - - - - - - - - - -

   size_rows( RampUp_Const , f_time_horizon );

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {

    vars.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );

    // on at t - 1 and t the row is p_t - p_{t-1} <= DeltaRampUp[ t ], at a
    // shut-down at t it is p_{t-1} >= MinPower[ t - 1 ]: the shut-down term
    // makes up for the difference of the minimum powers, and it is there
    // (with coefficient 0) also when they are equal, so that the Variable
    // of the row do not depend on the data [see update_rows()]
    if( t > 0 ) {
     vars.push_back( std::make_pair( &v_active_power[ t - 1 ] , -1.0 ) );
     vars.push_back( std::make_pair( &v_commitment[ t - 1 ] ,
                                     get_operational_min_power( t ) ) );
     if( t >= init_t )
      vars.push_back( std::make_pair( &v_shut_down[ t - init_t ] ,
                                      get_operational_min_power( t - 1 ) -
                                      get_operational_min_power( t ) ) );
    }

    if( t >= init_t )
     vars.push_back( std::make_pair( &v_start_up[ t - init_t ] ,
                                     -( v_StartUpLimit[ t ] -
                                        get_operational_min_power( t ) -
                                        v_DeltaRampUp[ t ] ) ) );

    vars.push_back( std::make_pair( &v_commitment[ t ] ,
                                    -( v_DeltaRampUp[ t ] +
                                       get_operational_min_power( t ) ) ) );

    put_row( RampUp_Const , t , std::move( vars ) , -Inf< double >() ,
             ( ( t == 0 ) && ( f_InitUpDownTime > 0 ) ) ?
             f_InitialPower - get_operational_min_power( t ) : 0.0 );
   }

  }
  if( ( AR & FormMsk ) == ptForm ) {  // pt formulation - - - - - - - - - - -

   if( f_InitUpDownTime > 0 )
    size_rows( RampUp_Const , f_time_horizon );
   if( f_InitUpDownTime <= 0 )
    size_rows( RampUp_Const , f_time_horizon - 1 );

   auto cnstr_idx = 0;

   for( Index t = 0 ; t < f_time_horizon ; ++t )
    if( ( t > 0 ) || ( ( t == 0 ) && ( f_InitUpDownTime > 0 ) ) ) {

     vars.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );

     if( t > 0 ) {
      vars.push_back( std::make_pair( &v_active_power[ t - 1 ] , -1.0 ) );
      for( Index i = 0 ; i < v_Y_plus.size() ; ++i ) {
       if( ( v_Y_plus[ i ].first <= t ) && ( t < v_Y_plus[ i ].second ) )
        vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                        -v_DeltaRampUp[ t ] ) );
       if( ( v_Y_plus[ i ].first <= t ) && ( t == v_Y_plus[ i ].second ) )
        vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                        get_operational_min_power( t - 1 ) ) );
       if( ( v_Y_plus[ i ].first == t + 1 ) &&
           ( t + 1 <= v_Y_plus[ i ].second ) )
        vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                        -v_StartUpLimit[ t ] ) );
      }
     }

     if( ( t == 0 ) && ( f_InitUpDownTime > 0 ) )
      for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
       if( ( v_Y_plus[ i ].first == 0 ) && ( t + 1 <= v_Y_plus[ i ].second ) )
        vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                        -f_InitialPower - v_DeltaRampUp[ t ] ) );

     put_row( RampUp_Const , cnstr_idx , std::move( vars ) ,
              -Inf< double >() , 0.0 );

     cnstr_idx++;
    }

  }
  if( ( AR & FormMsk ) == DPForm ) {  // DP formulation - - - - - - - - - - -

   auto ramp_up_cnstrs_size = 0;

   for( Index j = 0 ; j < v_P_h_k.size() ; ++j ) {
    auto t = v_P_h_k[ j ].first;
    if( ( ( t == 0 ) && ( f_InitUpDownTime > 0 ) ) ||
        ( ( t > 0 ) && ( v_P_h_k[ j ].second.first + 1 <= t + 1 ) &&
          ( v_P_h_k[ j ].second.second >= t + 1 ) ) )
     ramp_up_cnstrs_size++;
   }

   if( ramp_up_cnstrs_size > 0 ) {

    size_rows( RampUp_Const , ramp_up_cnstrs_size );

    auto cnstr_idx = 0;

    for( Index j = 0 ; j < v_P_h_k.size() ; ++j ) {

     auto t = v_P_h_k[ j ].first;
     if( ( ( t == 0 ) && ( f_InitUpDownTime > 0 ) ) ||
         ( ( t > 0 ) && ( v_P_h_k[ j ].second.first + 1 <= t + 1 ) &&
           ( v_P_h_k[ j ].second.second >= t + 1 ) ) ) {

      vars.push_back( std::make_pair( &v_active_power_h_k[ j ] , 1.0 ) );

      if( t > 0 )
       for( Index s = 0 ; s < v_P_h_k.size() ; ++s )
        if( ( v_P_h_k[ j ].second.first == v_P_h_k[ s ].second.first ) &&
            ( v_P_h_k[ j ].second.second == v_P_h_k[ s ].second.second ) )
         if( v_P_h_k[ s ].first == t - 1 )
          vars.push_back( std::make_pair( &v_active_power_h_k[ s ] , -1.0 ) );

      for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
       if( ( v_P_h_k[ j ].second.first == v_Y_plus[ i ].first ) &&
           ( v_P_h_k[ j ].second.second == v_Y_plus[ i ].second ) ) {
        if( t == 0 )
         vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                         -v_DeltaRampUp[ t ] - f_InitialPower ) );
        else if( t > 0 )
         vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                         -v_DeltaRampUp[ t ] ) );
       }

      put_row( RampUp_Const , cnstr_idx , std::move( vars ) ,
               -Inf< double >() , 0.0 );

      cnstr_idx++;
     }
    }
   }

  }
  if( ( AR & FormMsk ) == SUForm ) {  // SU formulation- - - - - - - - - - - -

   auto ramp_up_cnstrs_size = 0;

   for( Index t = 0 ; t < f_time_horizon ; ++t )
    for( Index j = 0 ; j < v_P_h.size() ; ++j )
     if( v_P_h[ j ].first == t ) {
      if( ( f_InitUpDownTime > 0 ) && ( t == 0 ) )
       ramp_up_cnstrs_size++;
      if( ( t > 0 ) && ( t + 1 > v_P_h[ j ].second ) )
       ramp_up_cnstrs_size++;
     }

   size_rows( RampUp_Const , ramp_up_cnstrs_size );

   auto cnstr_idx = 0;

   for( Index t = 0 ; t < f_time_horizon ; ++t )
    for( Index j = 0 ; j < v_P_h.size() ; ++j )
     if( v_P_h[ j ].first == t ) {
      if( ( f_InitUpDownTime > 0 ) && ( t == 0 ) ) {

       vars.push_back( std::make_pair( &v_active_power_h[ j ] , 1.0 ) );

       for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
        if( v_P_h[ j ].second == v_Y_plus[ i ].first )
         if( t + 1 <= v_Y_plus[ i ].second )
          vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                          -v_DeltaRampUp[ t ] - f_InitialPower ) );

       put_row( RampUp_Const , cnstr_idx , std::move( vars ) ,
                -Inf< double >() , 0.0 );

       cnstr_idx++;
      }
      if( ( t > 0 ) && ( t + 1 > v_P_h[ j ].second ) ) {

       vars.push_back( std::make_pair( &v_active_power_h[ j ] , 1.0 ) );

       for( Index s = 0 ; s < v_P_h.size() ; ++s )
        if( ( v_P_h[ j ].first - 1 == v_P_h[ s ].first ) &&
            ( v_P_h[ j ].second == v_P_h[ s ].second ) )
         vars.push_back( std::make_pair( &v_active_power_h[ s ] , -1.0 ) );
       for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
        if( v_P_h[ j ].second == v_Y_plus[ i ].first ) {
         if( t + 1 <= v_Y_plus[ i ].second )
          vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                          -v_DeltaRampUp[ t ] ) );
         if( t == v_Y_plus[ i ].second )
          vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                          get_operational_min_power( t - 1 ) ) );
        }

       put_row( RampUp_Const , cnstr_idx , std::move( vars ) ,
                -Inf< double >() , 0.0 );

       cnstr_idx++;
      }
     }
  }

  if( ( AR & FormMsk ) == SUSDForm ) {  // SUSD formulation - - - - - - - - -

   // the rows of the ramps of k steps, k being up to the number of steps
   // v_MaxRampSteps[] that the bounds leave, whose number depends on the
   // data: they are dynamic, keyed by the indices of the loops, and
   // update_rows() adds and removes them

   if( f_InitUpDownTime > 0 )
    for( Index k = 0 ; k < v_MaxRampSteps[ 0 ] ; ++k )
     for( Index j = 0 ; j < v_P_h.size() ; ++j )
      if( ( v_P_h[ j ].first == k ) && ( v_P_h[ j ].second == 0 ) ) {

       vars.push_back( std::make_pair( &v_active_power_h[ j ] , 1.0 ) );

       for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
        if( v_P_h[ j ].second == v_Y_plus[ i ].first )
         if( k + 1 <= v_Y_plus[ i ].second )
          vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                          -ramp_up_sum( 0 , k ) - f_InitialPower ) );

       put_dyn_row( RampUpSUSD_Const , { 0 , int( k ) , int( j ) } ,
                    std::move( vars ) , -Inf< double >() , 0.0 );
      }

   for( Index t = 0 ; t < f_time_horizon ; ++t )
    for( Index j = 0 ; j < v_P_h.size() ; ++j )
     if( v_P_h[ j ].first == t )
      for( Index s = 0 ; s < v_P_h.size() ; ++s )
       if( v_P_h[ j ].second == v_P_h[ s ].second )
        for( Index k = 1 ; k <= v_MaxRampSteps[ t + 1 ] ; ++k )
         if( v_P_h[ s ].first == t + k ) {

          vars.push_back( std::make_pair( &v_active_power_h[ s ] , 1.0 ) );
          vars.push_back( std::make_pair( &v_active_power_h[ j ] , -1.0 ) );

          for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
           if( v_P_h[ j ].second == v_Y_plus[ i ].first ) {
            if( t + k + 1 <= v_Y_plus[ i ].second )
             vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                             -ramp_up_sum( t + 1 , t + k ) ) );
            if( ( t + 1 <= v_Y_plus[ i ].second ) &&
                ( t + k + 1 > v_Y_plus[ i ].second ) )
             vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                             get_operational_min_power( t ) ) );
           }

          put_dyn_row( RampUpSUSD_Const ,
                       { 1 , int( t ) , int( j ) , int( s ) , int( k ) } ,
                       std::move( vars ) , -Inf< double >() , 0.0 );
         }

   add_dyn_rows( RampUpSUSD_Const , "RampUp_SUSD_Const_Thermal" );
  }
  if( ( ( AR & FormMsk ) == SDForm ) ||  // SD and SUSD formulations - - - -
      ( ( AR & FormMsk ) == SUSDForm ) ) {

   auto ramp_up_cnstrs_size = 0;

   for( Index t = 0 ; t < f_time_horizon ; ++t )
    for( Index j = 0 ; j < v_P_k.size() ; ++j )
     if( v_P_k[ j ].first == t ) {
      if( ( f_InitUpDownTime > 0 ) && ( t == 0 ) )
       ramp_up_cnstrs_size++;
      if( ( t > 0 ) && ( t + 1 <= v_P_k[ j ].second ) )
       ramp_up_cnstrs_size++;
     }

   size_rows( RampUp_Const , ramp_up_cnstrs_size );

   auto cnstr_idx = 0;

   for( Index t = 0 ; t < f_time_horizon ; ++t )
    for( Index j = 0 ; j < v_P_k.size() ; ++j )
     if( v_P_k[ j ].first == t ) {
      if( ( f_InitUpDownTime > 0 ) && ( t == 0 ) ) {

       vars.push_back( std::make_pair( &v_active_power_k[ j ] , 1.0 ) );

       for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
        if( v_P_k[ j ].second == v_Y_plus[ i ].second )
         if( v_Y_plus[ i ].first <= t )
          vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                          -v_DeltaRampUp[ t ] - f_InitialPower ) );

       put_row( RampUp_Const , cnstr_idx , std::move( vars ) ,
                -Inf< double >() , 0.0 );

       cnstr_idx++;
      }
      if( ( t > 0 ) && ( t + 1 <= v_P_k[ j ].second ) ) {

       vars.push_back( std::make_pair( &v_active_power_k[ j ] , 1.0 ) );

       for( Index s = 0 ; s < v_P_k.size() ; ++s )
        if( ( v_P_k[ j ].first - 1 == v_P_k[ s ].first ) &&
            ( v_P_k[ j ].second == v_P_k[ s ].second ) )
         vars.push_back( std::make_pair( &v_active_power_k[ s ] , -1.0 ) );
       for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
        if( v_P_k[ j ].second == v_Y_plus[ i ].second ) {
         if( v_Y_plus[ i ].first <= t )
          vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                          -v_DeltaRampUp[ t ] ) );
         if( t + 1 == v_Y_plus[ i ].first )
          vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                          -v_StartUpLimit[ t ] ) );
        }

       put_row( RampUp_Const , cnstr_idx , std::move( vars ) ,
                -Inf< double >() , 0.0 );

       cnstr_idx++;
      }
     }
  }

  add_rows( RampUp_Const , "RampUp_Const_Thermal" );
  }

 // initializing ramp-down constraints- - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! v_DeltaRampDown.empty() ) {
  if( ( AR & FormMsk ) == tbinForm ) {  // 3bin formulation - - - - - - - - -

   size_rows( RampDown_Const , f_time_horizon );

   auto cnstr_idx = 0;

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {

    vars.push_back( std::make_pair( &v_active_power[ t ] , -1.0 ) );
    vars.push_back( std::make_pair( &v_commitment[ t ] ,
                                    -v_DeltaRampDown[ t ] ) );

    if( t > 0 )
     vars.push_back( std::make_pair( &v_active_power[ t - 1 ] , 1.0 ) );

    if( t >= init_t )
     vars.push_back( std::make_pair( &v_shut_down[ t - init_t ] ,
                                     -v_ShutDownLimit[ t ] ) );

    put_row( RampDown_Const , cnstr_idx++ , std::move( vars ) ,
             -Inf< double >() ,
             ( ( t == 0 ) && ( f_InitUpDownTime > 0 ) ) ?
             -f_InitialPower : 0.0 );
   }

  }
  if( ( AR & FormMsk ) == TForm ) {  // T formulation - - - - - - - - - - - -

   size_rows( RampDown_Const , f_time_horizon );

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {

    vars.push_back( std::make_pair( &v_active_power[ t ] , -1.0 ) );

    if( t > 0 ) {
     vars.push_back( std::make_pair( &v_active_power[ t - 1 ] , 1.0 ) );
     vars.push_back( std::make_pair( &v_commitment[ t - 1 ] ,
                                     -( v_DeltaRampDown[ t ] +
                                        get_operational_min_power( t ) ) ) );
    }

    if( t >= init_t )
     vars.push_back( std::make_pair( &v_shut_down[ t - init_t ] ,
                                     -( v_ShutDownLimit[ t ] -
                                        get_operational_min_power( t ) -
                                        v_DeltaRampDown[ t ] ) ) );

    vars.push_back( std::make_pair( &v_commitment[ t ] ,
                                    get_operational_min_power( t ) ) );

    put_row( RampDown_Const , t , std::move( vars ) , -Inf< double >() ,
             ( ( t == 0 ) && ( f_InitUpDownTime > 0 ) ) ?
             -( f_InitialPower - v_DeltaRampDown[ t ] -
                get_operational_min_power( t ) ) : 0.0 );
   }

  }
  if( ( AR & FormMsk ) == ptForm ) {  // pt formulation - - - - - - - - - - -

   if( f_InitUpDownTime > 0 )
    size_rows( RampDown_Const , f_time_horizon );
   if( f_InitUpDownTime <= 0 )
    size_rows( RampDown_Const , f_time_horizon - 1 );

   auto cnstr_idx = 0;

   for( Index t = 0 ; t < f_time_horizon ; ++t )
    if( ( t > 0 ) || ( ( t == 0 ) && ( f_InitUpDownTime > 0 ) ) ) {

     vars.push_back( std::make_pair( &v_active_power[ t ] , -1.0 ) );

     if( t > 0 ) {
      vars.push_back( std::make_pair( &v_active_power[ t - 1 ] , 1.0 ) );
      for( Index i = 0 ; i < v_Y_plus.size() ; ++i ) {
       if( ( v_Y_plus[ i ].first <= t ) && ( t < v_Y_plus[ i ].second ) )
        vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                        -v_DeltaRampDown[ t ] ) );
       if( ( v_Y_plus[ i ].first <= t ) && ( t == v_Y_plus[ i ].second ) )
        vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                        -v_ShutDownLimit[ t ] ) );
       if( ( v_Y_plus[ i ].first == t + 1 ) &&
           ( t + 1 <= v_Y_plus[ i ].second ) )
        vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                        get_operational_min_power( t ) ) );
      }
     }

     if( ( t == 0 ) && ( f_InitUpDownTime > 0 ) )
      for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
       if( ( v_Y_plus[ i ].first == 0 ) && ( t + 1 <= v_Y_plus[ i ].second ) )
        vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                        f_InitialPower - v_DeltaRampDown[ t ] ) );

     put_row( RampDown_Const , cnstr_idx , std::move( vars ) ,
              -Inf< double >() , 0.0 );

     cnstr_idx++;
    }

  }
  if( ( AR & FormMsk ) == DPForm ) {  // DP formulation - - - - - - - - - - -

   auto ramp_up_cnstrs_size = 0;

   for( Index j = 0 ; j < v_P_h_k.size() ; ++j ) {
    auto t = v_P_h_k[ j ].first;
    if( ( ( t == 0 ) && ( f_InitUpDownTime > 0 ) ) ||
        ( ( t > 0 ) && ( v_P_h_k[ j ].second.first + 1 <= t + 1 ) &&
          ( v_P_h_k[ j ].second.second >= t + 1 ) ) )
     ramp_up_cnstrs_size++;
   }

   if( ramp_up_cnstrs_size > 0 ) {

    size_rows( RampDown_Const , ramp_up_cnstrs_size );

    auto cnstr_idx = 0;

    for( Index j = 0 ; j < v_P_h_k.size() ; ++j ) {

     auto t = v_P_h_k[ j ].first;
     if( ( ( t == 0 ) && ( f_InitUpDownTime > 0 ) ) ||
         ( ( t > 0 ) && ( v_P_h_k[ j ].second.first + 1 <= t + 1 ) &&
           ( v_P_h_k[ j ].second.second >= t + 1 ) ) ) {

      vars.push_back( std::make_pair( &v_active_power_h_k[ j ] , -1.0 ) );

      if( t > 0 )
       for( Index s = 0 ; s < v_P_h_k.size() ; ++s )
        if( ( v_P_h_k[ j ].second.first == v_P_h_k[ s ].second.first ) &&
            ( v_P_h_k[ j ].second.second == v_P_h_k[ s ].second.second ) )
         if( v_P_h_k[ s ].first == t - 1 )
          vars.push_back( std::make_pair( &v_active_power_h_k[ s ] , 1.0 ) );

      for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
       if( ( v_P_h_k[ j ].second.first == v_Y_plus[ i ].first ) &&
           ( v_P_h_k[ j ].second.second == v_Y_plus[ i ].second ) ) {
        if( t == 0 )
         vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                         -v_DeltaRampDown[ t ] + f_InitialPower ) );
        else if( t > 0 )
         vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                         -v_DeltaRampDown[ t ] ) );
       }

      put_row( RampDown_Const , cnstr_idx , std::move( vars ) ,
               -Inf< double >() , 0.0 );

      cnstr_idx++;
     }
    }
   }

  }
  if( ( ( AR & FormMsk ) == SUForm ) ||  // SU and SUSD formulations - - - -
      ( ( AR & FormMsk ) == SUSDForm ) ) {

   auto ramp_down_cnstrs_size = 0;

   for( Index t = 0 ; t < f_time_horizon ; ++t )
    for( Index j = 0 ; j < v_P_h.size() ; ++j )
     if( v_P_h[ j ].first == t ) {
      if( ( f_InitUpDownTime > 0 ) && ( t == 0 ) )
       ramp_down_cnstrs_size++;
      if( ( t > 0 ) && ( t + 1 > v_P_h[ j ].second ) )
       ramp_down_cnstrs_size++;
     }

   size_rows( RampDown_Const , ramp_down_cnstrs_size );

   auto cnstr_idx = 0;

   for( Index t = 0 ; t < f_time_horizon ; ++t )
    for( Index j = 0 ; j < v_P_h.size() ; ++j )
     if( v_P_h[ j ].first == t ) {
      if( ( f_InitUpDownTime > 0 ) && ( t == 0 ) ) {

       vars.push_back( std::make_pair( &v_active_power_h[ j ] , -1.0 ) );

       for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
        if( v_P_h[ j ].second == v_Y_plus[ i ].first )
         if( t + 1 <= v_Y_plus[ i ].second )
          vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                          -v_DeltaRampDown[ t ] + f_InitialPower ) );

       put_row( RampDown_Const , cnstr_idx , std::move( vars ) ,
                -Inf< double >() , 0.0 );

       cnstr_idx++;
      }

      if( ( t > 0 ) && ( t + 1 > v_P_h[ j ].second ) ) {

       vars.push_back( std::make_pair( &v_active_power_h[ j ] , -1.0 ) );

       for( Index s = 0 ; s < v_P_h.size() ; ++s )
        if( ( v_P_h[ j ].first - 1 == v_P_h[ s ].first ) &&
            ( v_P_h[ j ].second == v_P_h[ s ].second ) )
         vars.push_back( std::make_pair( &v_active_power_h[ s ] , 1.0 ) );
       for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
        if( v_P_h[ j ].second == v_Y_plus[ i ].first ) {
         if( t + 1 <= v_Y_plus[ i ].second )
          vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                          -v_DeltaRampDown[ t ] ) );
         if( t == v_Y_plus[ i ].second )
          vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                          -v_ShutDownLimit[ t ] ) );
        }

       put_row( RampDown_Const , cnstr_idx , std::move( vars ) ,
                -Inf< double >() , 0.0 );

       cnstr_idx++;
      }
     }

  }
  if( ( AR & FormMsk ) == SUSDForm ) {  // SUSD formulation - - - - - - - - -

   // dynamic, as the ramp-up rows [see above]

   if( f_InitUpDownTime > 0 )
    for( Index k = 0 ; k < v_MaxRampDownSteps[ 0 ] ; ++k )
     for( Index j = 0 ; j < v_P_k.size() ; ++j )
      if( v_P_k[ j ].first == k ) {

       vars.push_back( std::make_pair( &v_active_power_k[ j ] , -1.0 ) );

       for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
   if( v_P_k[ j ].second == v_Y_plus[ i ].second )
         if( 0 == v_Y_plus[ i ].first )
          vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                          -ramp_down_sum( 0 , k ) + f_InitialPower ) );

       put_dyn_row( RampDownSUSD_Const , { 0 , int( k ) , int( j ) } ,
                    std::move( vars ) , -Inf< double >() , 0.0 );
      }

   for( Index t = 0 ; t < f_time_horizon ; ++t )
    for( Index j = 0 ; j < v_P_k.size() ; ++j )
     if( v_P_k[ j ].first == t )
      for( Index s = 0 ; s < v_P_k.size() ; ++s )
       if( v_P_k[ j ].second == v_P_k[ s ].second )
        for( Index k = 1 ; k <= v_MaxRampDownSteps[ t + 1 ] ; ++k )
         if( v_P_k[ s ].first == t + k ) {

          vars.push_back( std::make_pair( &v_active_power_k[ s ] , -1.0 ) );
          vars.push_back( std::make_pair( &v_active_power_k[ j ] , 1.0 ) );

          for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
           if( v_P_k[ j ].second == v_Y_plus[ i ].second ) {
            if( t + 1 >= v_Y_plus[ i ].first )
             vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                             -ramp_down_sum( t + 1 , t + k ) ) );
            if( ( t + 1 < v_Y_plus[ i ].first ) &&
                ( t + k + 1 >= v_Y_plus[ i ].first ) )
             vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                    get_operational_min_power( t + k ) ) );
           }

          put_dyn_row( RampDownSUSD_Const ,
                       { 1 , int( t ) , int( j ) , int( s ) , int( k ) } ,
                       std::move( vars ) , -Inf< double >() , 0.0 );
         }

   add_dyn_rows( RampDownSUSD_Const , "RampDown_SUSD_Const_Thermal" );
  }
  if( ( AR & FormMsk ) == SDForm ) {  // SD formulation - - - - - - - - - - -
   auto ramp_down_cnstrs_size = 0;

   for( Index t = 0 ; t < f_time_horizon ; ++t )
    for( Index j = 0 ; j < v_P_k.size() ; ++j )
     if( v_P_k[ j ].first == t ) {
      if( ( f_InitUpDownTime > 0 ) && ( t == 0 ) )
       ramp_down_cnstrs_size++;
      if( ( t > 0 ) && ( t + 1 <= v_P_k[ j ].second ) )
       ramp_down_cnstrs_size++;
     }

   size_rows( RampDown_Const , ramp_down_cnstrs_size );

   auto cnstr_idx = 0;

   for( Index t = 0 ; t < f_time_horizon ; ++t )
    for( Index j = 0 ; j < v_P_k.size() ; ++j )
     if( v_P_k[ j ].first == t ) {
      if( ( f_InitUpDownTime > 0 ) && ( t == 0 ) ) {

       vars.push_back( std::make_pair( &v_active_power_k[ j ] , -1.0 ) );
       for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
        if( v_P_k[ j ].second == v_Y_plus[ i ].second )
         if( v_Y_plus[ i ].first <= t )
          vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                          -v_DeltaRampDown[ t ] + f_InitialPower ) );

       put_row( RampDown_Const , cnstr_idx , std::move( vars ) ,
                -Inf< double >() , 0.0 );

       cnstr_idx++;
      }
      if( ( t > 0 ) && ( t + 1 <= v_P_k[ j ].second ) ) {

       vars.push_back( std::make_pair( &v_active_power_k[ j ] , -1.0 ) );

       for( Index s = 0 ; s < v_P_k.size() ; ++s )
        if( ( v_P_k[ j ].first - 1 == v_P_k[ s ].first ) &&
            ( v_P_k[ j ].second == v_P_k[ s ].second ) )
         vars.push_back( std::make_pair( &v_active_power_k[ s ] , 1.0 ) );
       for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
        if( v_P_k[ j ].second == v_Y_plus[ i ].second ) {
         if( v_Y_plus[ i ].first <= t )
          vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                          -v_DeltaRampDown[ t ] ) );
         if( t + 1 == v_Y_plus[ i ].first )
          vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                          get_operational_min_power( t ) ) );
        }

       put_row( RampDown_Const , cnstr_idx , std::move( vars ) ,
                -Inf< double >() , 0.0 );

       cnstr_idx++;
      }
     }
  }

  add_rows( RampDown_Const , "RampDown_Const_Thermal" );
 }

 // Initializing minimum power constraints- - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( ( AR & FormMsk ) == tbinForm ) ||  // 3bin formulation - - - - - - - -
     ( ( AR & FormMsk ) == TForm ) ) {  // T formulation- - - - - - - - - - -

  size_rows( MinPower_Const , f_time_horizon );

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {

   vars.push_back( std::make_pair( &v_commitment[ t ] ,
                                   -get_operational_min_power( t ) ) );
   vars.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );



   put_row( MinPower_Const , t , std::move( vars ) , 0.0 , Inf< double >() );
  }

 }
 if( ( AR & FormMsk ) == ptForm ) {  // pt formulation - - - - - - - -

  size_rows( MinPower_Const , f_time_horizon );

  for( Index t = 0 , constraint_index = 0 ; t < f_time_horizon ;
       ++t , ++constraint_index ) {

   vars.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );



   for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
    if( ( v_Y_plus[ i ].first <= t + 1 ) &&
        ( t + 1 <= v_Y_plus[ i ].second ) )
     vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                     -get_operational_min_power( t ) ) );

   put_row( MinPower_Const , constraint_index , std::move( vars ) ,
            0.0 , Inf< double >() );
  }

 }
 if( ( AR & FormMsk ) == DPForm ) {  // DP formulation - - - - - - - -

  size_rows( MinPower_Const , v_P_h_k.size() );

  for( Index j = 0 ; j < v_P_h_k.size() ; ++j ) {

   auto t = v_P_h_k[ j ].first;
   for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
    if( ( v_P_h_k[ j ].second.first == v_Y_plus[ i ].first ) &&
        ( v_P_h_k[ j ].second.second == v_Y_plus[ i ].second ) )
     vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                     -get_operational_min_power( t ) ) );

   vars.push_back( std::make_pair( &v_active_power_h_k[ j ] , 1.0 ) );

   put_row( MinPower_Const , j , std::move( vars ) , 0.0 , Inf< double >() );
  }
 }

 if( ( AR & FormMsk ) == SUForm || ( AR & FormMsk ) == SUSDForm )
 {  // SU and SUSD formulation- - - - - - - - - - - - - - - - - - - - - - - -

  if( ( AR & FormMsk ) == SUSDForm )
   size_rows( MinPower_Const , v_P_h.size() + v_P_k.size() );
  else
   size_rows( MinPower_Const , v_P_h.size() );

  for( Index j = 0 ; j < v_P_h.size() ; ++j ) {

   auto t = v_P_h[ j ].first;
   for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
    if( ( v_P_h[ j ].second == v_Y_plus[ i ].first ) &&
        ( t + 1 <= v_Y_plus[ i ].second ) )
     vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                     -get_operational_min_power( t ) ) );

   vars.push_back( std::make_pair( &v_active_power_h[ j ] , 1.0 ) );

   put_row( MinPower_Const , j , std::move( vars ) , 0.0 , Inf< double >() );
  }

 }
 if( ( AR & FormMsk ) == SDForm || ( AR & FormMsk ) == SUSDForm )
 {  // SD and SUSD formulation- - - - - - - - - - - - - - - - - - - - - - - -

  if( ( AR & FormMsk ) == SDForm )
   size_rows( MinPower_Const , v_P_k.size() );

  for( Index j = 0 ; j < v_P_k.size() ; ++j ) {

   auto t = v_P_k[ j ].first;
   for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
    if( ( v_P_k[ j ].second == v_Y_plus[ i ].second ) &&
        ( v_Y_plus[ i ].first <= t + 1 ) )
     vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                     -get_operational_min_power( t ) ) );

   vars.push_back( std::make_pair( &v_active_power_k[ j ] , 1.0 ) );

   if( ( AR & FormMsk ) == SUSDForm ) {
    put_row( MinPower_Const , v_P_h.size() + j , std::move( vars ) ,
             0.0 , Inf< double >() );
   } else {
    put_row( MinPower_Const , j , std::move( vars ) , 0.0 , Inf< double >() );
   }
  }
 }

 add_rows( MinPower_Const , "MinPower_Const_Thermal" );

 // initializing maximum power constraints- - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( AR & FormMsk ) == tbinForm ) {  // 3bin formulation- - - - - - - - - -

  // one row per instant, and a second one at each init_t <= t < T - 1 when
  // the minimum up time is 1; the start-up at t is capped by StartUpLimit[ t ]
  // and the last power before a shut-down at t + 1 by sd_cap( t ), i.e.,
  // ShutDownLimit[ t + 1 ] if not above the maximum power at t, at t == 0 as
  // at any other instant (with T == 1 there is no shut-down)
  Index n_max_rows = f_time_horizon;
  if( f_MinUpTime == 1 )
   for( Index t = init_t ; t + 1 < f_time_horizon ; ++t )
    ++n_max_rows;
  size_rows( MaxPower_Const , n_max_rows );

  for( Index t = 0 , cnstr_idx = 0 ; t < f_time_horizon ; ++t , ++cnstr_idx ) {

   if( t >= init_t ) {
    if( t + 1 < f_time_horizon ) {
     vars.push_back( std::make_pair( &v_shut_down[ t + 1 - init_t ] ,
                                     sd_cap( t ) -
                                     get_operational_max_power( t ) ) );
     vars.push_back( std::make_pair( &v_start_up[ t - init_t ] ,
                                     f_MinUpTime == 1 ?
                                     -std::max( 0.0 , sd_cap( t ) - std::min(
                                      double( v_StartUpLimit[ t ] ) ,
                                      get_operational_max_power( t ) ) ) :
                                     v_StartUpLimit[ t ] -
                                     get_operational_max_power( t ) ) );
     }
    else
     vars.push_back( std::make_pair( &v_start_up[ t - init_t ] ,
                                     v_StartUpLimit[ t ] -
                                     get_operational_max_power( t ) ) );
    }
   else
    if( ( t + 1 == init_t ) && ( init_t < f_time_horizon ) )
     // the last instant of the run the state before the horizon imposes,
     // capped by sd_cap( t ) if the unit shuts down at init_t
     vars.push_back( std::make_pair( &v_shut_down[ 0 ] , sd_cap( t ) -
                                     get_operational_max_power( t ) ) );

   vars.push_back( std::make_pair( &v_commitment[ t ] ,
                                   get_operational_max_power( t ) ) );
   vars.push_back( std::make_pair( &v_active_power[ t ] , -1.0 ) );

   put_row( MaxPower_Const , cnstr_idx , std::move( vars ) ,
            0.0 , Inf< double >() );

   if( ( t >= init_t ) && ( f_MinUpTime == 1 ) &&
       ( t + 1 < f_time_horizon ) ) {
    vars.push_back( std::make_pair( &v_shut_down[ t + 1 - init_t ] ,
                                    -std::max( 0.0 , std::min(
                                     double( v_StartUpLimit[ t ] ) ,
                                     get_operational_max_power( t ) ) -
                                     sd_cap( t ) ) ) );
    vars.push_back( std::make_pair( &v_start_up[ t - init_t ] ,
                                    v_StartUpLimit[ t ] -
                                    get_operational_max_power( t ) ) );
    vars.push_back( std::make_pair( &v_commitment[ t ] ,
                                    get_operational_max_power( t ) ) );
    vars.push_back( std::make_pair( &v_active_power[ t ] , -1.0 ) );

    ++cnstr_idx;
    put_row( MaxPower_Const , cnstr_idx , std::move( vars ) ,
             0.0 , Inf< double >() );
    }
   }

 }
 if( ( AR & FormMsk ) == TForm ) {  // T formulation - - - - - - - - -

  Index max_power_cnstrs_size = 0;

  std::vector< int > v_T_RU;
  std::vector< int > v_K_SU;
  std::vector< int > v_K_SD;

  if( ( ! v_DeltaRampUp.empty() ) || ( ! v_DeltaRampDown.empty() ) ) {
   std::vector< int > v_T_RD;

   if( ( ! v_DeltaRampUp.empty() ) ) {
    v_T_RU.resize( f_time_horizon );
    v_K_SU.resize( f_time_horizon );
   }
   if( ( ! v_DeltaRampDown.empty() ) ) {
    v_T_RD.resize( f_time_horizon );
    v_K_SD.resize( f_time_horizon );
   }

   // the number of whole ramps in a distance, as an int; a ramp 0 (or one
   // so small that the count exceeds the horizon) never covers it, which is
   // the count T + MinUpTime, larger than every bound it is compared with
   const int no_cover = static_cast< int >( f_time_horizon + f_MinUpTime );
   auto ramps_in = [ no_cover ]( double dist , double ramp ) {
    if( ( ramp <= 0 ) || ( dist >= no_cover * ramp ) )
     return( no_cover );
    return( static_cast< int >( std::floor( dist / ramp ) ) );
    };

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {

    // TRU is the length of the start-up ramping trajectory (climb from the
    // start-up level to the cap), so it uses the StartUpLimit and DeltaRampUp;
    // TRD is the shut-down trajectory (descent to the shut-down level), using
    // the ShutDownLimit and DeltaRampDown
    if( ( ! v_DeltaRampUp.empty() ) )
     v_T_RU[ t ] = ramps_in( get_operational_max_power( t ) -
                             v_StartUpLimit[ t ] , v_DeltaRampUp[ t ] );
    if( ( ! v_DeltaRampDown.empty() ) ) {
     v_T_RD[ t ] = ramps_in( get_operational_max_power( t ) -
                             v_ShutDownLimit[ t ] , v_DeltaRampDown[ t ] );
     // the window of shut-down / start-up trajectory steps that can jointly
     // constrain a period is bounded by the minimum up time (the shortest an
     // on-interval can be), not by the initial up/down time -- the strengthening
     // must not depend on the unit's state before the horizon.
     v_K_SD[ t ] = std::min( static_cast< int >( f_MinUpTime ) - 1 ,
                             v_T_RD[ t ] );
     v_K_SD[ t ] = std::min( static_cast< int >( f_time_horizon - t ) - 2 ,
                             v_K_SD[ t ] );
    }
    if( ( ! v_DeltaRampUp.empty() ) && ( ! v_DeltaRampDown.empty() ) ) {
     v_K_SU[ t ] = std::min( static_cast< int >( f_MinUpTime ) - 2 -
                             std::max( 0 , v_K_SD[ t ] ) , v_T_RU[ t ] );
     v_K_SU[ t ] = std::min( static_cast< int >( t ) -
                             static_cast< int >( init_t ) , v_K_SU[ t ] );
    }
    if( ( ! v_DeltaRampUp.empty() ) && ( v_DeltaRampDown.empty() ) )
     v_K_SU[ t ] = std::min( static_cast< int >( t ) -
                             static_cast< int >( init_t ) , v_T_RU[ t ] );
   }

   // size bound constraints 4; the bound constraints 5 and 6, whose
   // number depends on the data, are dynamic [see below]
   if( ( ! v_DeltaRampUp.empty() ) )
    max_power_cnstrs_size += f_time_horizon - init_t;
  }

  // size bound constraints 0
  max_power_cnstrs_size += f_time_horizon;

  // size bound constraints 1
  if( f_MinUpTime > 1 )
   max_power_cnstrs_size += f_time_horizon - init_t;
  // size bound constraints 2 + 3
  else
   max_power_cnstrs_size += 2 * f_time_horizon - 2 * init_t;


  size_rows( MaxPower_Const , max_power_cnstrs_size );

  auto cnstr_idx = 0;

  // Bound constraints 0- - - - - - - - - - - - - - - - - - - - - - - - - - -

  for( Index t = 0 ; t < f_time_horizon ; ++t , ++cnstr_idx ) {

   vars.push_back( std::make_pair( &v_commitment[ t ] ,
                                   get_operational_max_power( t ) ) );
   vars.push_back( std::make_pair( &v_active_power[ t ] , -1.0 ) );

   // the last instant of the run the state before the horizon imposes,
   // capped by sd_cap( t ) if the unit shuts down at init_t
   if( ( t + 1 == init_t ) && ( init_t < f_time_horizon ) )
    vars.push_back( std::make_pair( &v_shut_down[ 0 ] , sd_cap( t ) -
                                    get_operational_max_power( t ) ) );

   put_row( MaxPower_Const , cnstr_idx , std::move( vars ) ,
            0.0 , Inf< double >() );
  }

  // Bound constraints 1- - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( f_MinUpTime > 1 )
   for( Index t = init_t ; t < f_time_horizon ; ++t , ++cnstr_idx ) {

    vars.push_back( std::make_pair( &v_commitment[ t ] ,
                                    get_operational_max_power( t ) ) );
    vars.push_back( std::make_pair( &v_active_power[ t ] , -1.0 ) );



    if( t >= init_t ) {
     vars.push_back( std::make_pair( &v_start_up[ t - init_t ] ,
                                     -( get_operational_max_power( t ) -
                                        v_StartUpLimit[ t ] ) ) );
     if( t < ( f_time_horizon - 1 ) )
      vars.push_back( std::make_pair( &v_shut_down[ t + 1 - init_t ] ,
                                      -( get_operational_max_power( t ) -
                                         sd_cap( t ) ) ) );
    }

    put_row( MaxPower_Const , cnstr_idx , std::move( vars ) ,
             0.0 , Inf< double >() );
   }

  if( f_MinUpTime == 1 ) {

   // Bound constraints 2 - - - - - - - - - - - - - - - - - - - - - - - - - -

   for( Index t = init_t ; t < f_time_horizon ; ++t , ++cnstr_idx ) {

    vars.push_back( std::make_pair( &v_commitment[ t ] ,
                                    get_operational_max_power( t ) ) );
    vars.push_back( std::make_pair( &v_active_power[ t ] , -1.0 ) );



    if( t >= init_t ) {
     vars.push_back( std::make_pair( &v_start_up[ t - init_t ] ,
                                     -( get_operational_max_power( t ) -
                                        v_StartUpLimit[ t ] ) ) );
     // a single on period t (v_t = w_{t+1} = 1) is capped by the smaller of
     // the start-up and shut-down caps (the term is there, with coefficient
     // 0, also when there is nothing to cap)
     if( t < ( f_time_horizon - 1 ) ) {
      const double gap = std::min( double( v_StartUpLimit[ t ] ) ,
                                   get_operational_max_power( t ) ) -
                         sd_cap( t );
      vars.push_back( std::make_pair( &v_shut_down[ t + 1 - init_t ] ,
                                      -std::max( gap , 0.0 ) ) );
      }
    }

    put_row( MaxPower_Const , cnstr_idx , std::move( vars ) ,
             0.0 , Inf< double >() );
   }

   // Bound constraints 3 - - - - - - - - - - - - - - - - - - - - - - - - - -

   for( Index t = init_t ; t < f_time_horizon ; ++t , ++cnstr_idx ) {

    vars.push_back( std::make_pair( &v_commitment[ t ] ,
                                    get_operational_max_power( t ) ) );
    vars.push_back( std::make_pair( &v_active_power[ t ] , -1.0 ) );



    if( t >= init_t ) {
     if( t < ( f_time_horizon - 1 ) ) {
      vars.push_back( std::make_pair( &v_shut_down[ t + 1 - init_t ] ,
                                      -( get_operational_max_power( t ) -
                                         sd_cap( t ) ) ) );
      const double gap = sd_cap( t ) -
                         std::min( double( v_StartUpLimit[ t ] ) ,
                                   get_operational_max_power( t ) );
      vars.push_back( std::make_pair( &v_start_up[ t - init_t ] ,
                                      -std::max( gap , 0.0 ) ) );
      }
    }

    put_row( MaxPower_Const , cnstr_idx , std::move( vars ) ,
             0.0 , Inf< double >() );
   }
  }

  if( ( ! v_DeltaRampUp.empty() ) || ( ! v_DeltaRampDown.empty() ) ) {

   // Bound constraints 4 - - - - - - - - - - - - - - - - - - - - - - - - - -
  if( ( ! v_DeltaRampUp.empty() ) )
   for( Index t = init_t ; t < f_time_horizon ; ++t , ++cnstr_idx ) {

    vars.push_back( std::make_pair( &v_commitment[ t ] ,
                                    get_operational_max_power( t ) ) );
    vars.push_back( std::make_pair( &v_active_power[ t ] , -1.0 ) );



    // the start-ups of the last MinUpTime - 1 instants (s = 0, ...,
    // MinUpTime - 2) are there, those that cannot bind with coefficient 0,
    // so that the Variable of the row do not depend on the data [see
    // update_rows()]; a unit started at t - s reaches at most the start-up
    // limit plus s ramps at t
    const int min_RU = static_cast< int >( f_MinUpTime ) - 2;

    if( t >= init_t ) {
     if( t < ( f_time_horizon - 1 ) )
       vars.push_back( std::make_pair( &v_shut_down[ t + 1 - init_t ] ,
                                       -( get_operational_max_power( t ) -
                                          sd_cap( t ) ) ) );

     for( int s = 0 ; s <= min_RU ; ++s )
      if( t - init_t >= Index( s ) )
       vars.push_back( std::make_pair(
        &v_start_up[ t - s - init_t ] ,
        -std::max( 0.0 , get_operational_max_power( t ) -
                         v_StartUpLimit[ t - s ] -
                         ramp_up_sum( t - s + 1 , t ) ) ) );
    }

    put_row( MaxPower_Const , cnstr_idx , std::move( vars ) ,
             0.0 , Inf< double >() );
   }

   // Bound constraints 5 (dynamic) - - - - - - - - - - - - - - - - - - - - -
 if( ( ! v_DeltaRampUp.empty() ) )
   for( Index t = init_t ; t < f_time_horizon ; ++t )
    if( ( f_MinUpTime - 2 ) < v_T_RU[ t ] ) {

     vars.push_back( std::make_pair( &v_commitment[ t ] ,
                                     get_operational_max_power( t ) ) );
     vars.push_back( std::make_pair( &v_active_power[ t ] , -1.0 ) );



     int min_RU = std::min( static_cast< int >( f_MinUpTime ) - 1 , v_T_RU[ t ] );

     if( t >= init_t ) {
      for( int s = 0 ; s <= min_RU ; ++s )
       if( t - init_t >= Index( s ) )
        vars.push_back( std::make_pair(
         &v_start_up[ t - s - init_t ] ,
         -std::max( 0.0 , get_operational_max_power( t ) -
                          v_StartUpLimit[ t - s ] -
                          ramp_up_sum( t - s + 1 , t ) ) ) );
     }

     put_dyn_row( MaxPower5_Const , { int( t ) } , std::move( vars ) ,
                  0.0 , Inf< double >() );
    }

   // Bound constraints 6 (dynamic) - - - - - - - - - - - - - - - - - - - - -
  if( ( ! v_DeltaRampUp.empty() ) && ( ! v_DeltaRampDown.empty() ) )
   for( Index t = init_t ; t < f_time_horizon ; ++t )
    if( v_K_SD[ t ] > 0 ) {

     vars.push_back( std::make_pair( &v_commitment[ t ] ,
                                     get_operational_max_power( t ) ) );
     vars.push_back( std::make_pair( &v_active_power[ t ] , -1.0 ) );



     if( t >= init_t ) {

      for( int s = 0 ; s <= v_K_SD[ t ] ; ++s )
       if( ( t + 1 + s ) < f_time_horizon )
        vars.push_back( std::make_pair(
         &v_shut_down[ t + 1 + s - init_t ] ,
         -std::max( 0.0 , get_operational_max_power( t ) -
                          v_ShutDownLimit[ t + 1 + s ] -
                          ramp_down_sum( t + 1 , t + s ) ) ) );

      for( int s = 0 ; s <= v_K_SU[ t ] ; ++s )
       if( t - init_t >= Index( s ) )
        vars.push_back( std::make_pair(
         &v_start_up[ t - s - init_t ] ,
         -std::max( 0.0 , get_operational_max_power( t ) -
                          v_StartUpLimit[ t - s ] -
                          ramp_up_sum( t - s + 1 , t ) ) ) );
     }

     put_dyn_row( MaxPower6_Const , { int( t ) } , std::move( vars ) ,
                  0.0 , Inf< double >() );
    }
  }

  // the bound constraints 5 and 6 exist at the instants where they are not
  // dominated by the bound constraints 4, which depends on the data: they
  // are dynamic, and update_rows() adds and removes them
  if( ! v_DeltaRampUp.empty() )
   add_dyn_rows( MaxPower5_Const , "MaxPower5_Const_Thermal" );
  if( ( ! v_DeltaRampUp.empty() ) && ( ! v_DeltaRampDown.empty() ) )
   add_dyn_rows( MaxPower6_Const , "MaxPower6_Const_Thermal" );

 } else if( ( AR & FormMsk ) == DPForm ) {  // DP formulation - - - - - - - -

  size_rows( MaxPower_Const , v_P_h_k.size() );

  for( Index j = 0 ; j < v_P_h_k.size() ; ++j ) {

   auto t = v_P_h_k[ j ].first;
   for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
    if( ( v_P_h_k[ j ].second.first == v_Y_plus[ i ].first ) &&
        ( v_P_h_k[ j ].second.second == v_Y_plus[ i ].second ) ) {
     // the first instant of the run is capped by the start-up limit, the
     // last by the shut-down one, a run of one instant by both
     if( v_Y_plus[ i ].first == t + 1 )
      vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                      v_Y_plus[ i ].second == t + 1 ?
                                      std::min( double( v_StartUpLimit[ t ] ) ,
                                                sd_cap( t ) ) :
                                      double( v_StartUpLimit[ t ] ) ) );
     else if( v_Y_plus[ i ].second == t + 1 )
      vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                      sd_cap( t ) ) );
     else
      vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                      get_operational_max_power( t ) ) );
    }

   vars.push_back( std::make_pair( &v_active_power_h_k[ j ] , -1.0 ) );

   put_row( MaxPower_Const , j , std::move( vars ) , 0.0 , Inf< double >() );
  }

 } else {  // pt, SU or SD formulations - - - - - - - - - - - - - - - - - - -


  if( ( AR & FormMsk ) == ptForm ) {  // pt formulation - - - - - - - - - - -

   size_rows( MaxPower_Const , f_time_horizon );

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {

    vars.push_back( std::make_pair( &v_active_power[ t ] , -1.0 ) );



    // every run on at t, with the cap \psi of the instant t of the run
    for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
     for( Index j : psi_index( t , i ) )
      vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                      v_psi[ j ] ) );

    put_row( MaxPower_Const , t , std::move( vars ) , 0.0 , Inf< double >() );
   }
  }

  if( ( AR & FormMsk ) == SUForm || ( AR & FormMsk ) == SUSDForm )
  {  // SU and SUSD formulation - - - - - - - - - - - - - - - - - - - - - - -

   if( ( AR & FormMsk ) == SUSDForm )
    size_rows( MaxPower_Const , v_P_h.size() + v_P_k.size() );
   else
    size_rows( MaxPower_Const , v_P_h.size() );

   for( Index j = 0 ; j < v_P_h.size() ; ++j ) {

    auto t = v_P_h[ j ].first;
    // every run started at h - 1 and on at t, with the cap \psi of the
    // instant t of the run
    for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
     if( ( v_P_h[ j ].second == v_Y_plus[ i ].first ) &&
         ( t + 1 <= v_Y_plus[ i ].second ) )
      for( Index s : psi_index( t , i ) )
       vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                       v_psi[ s ] ) );

    vars.push_back( std::make_pair( &v_active_power_h[ j ] , -1.0 ) );

    put_row( MaxPower_Const , j , std::move( vars ) , 0.0 , Inf< double >() );
   }
  }

  if( ( AR & FormMsk ) == SDForm || ( AR & FormMsk ) == SUSDForm )
  {  // SD and SUSD formulation - - - - - - - - - - - - - - - - - - - - - - -

   if( ( AR & FormMsk ) == SDForm )
    size_rows( MaxPower_Const , v_P_k.size() );

   for( Index j = 0 ; j < v_P_k.size() ; ++j ) {

    auto t = v_P_k[ j ].first;
    // every run shut down at k and on at t, with the cap \psi of the
    // instant t of the run
    for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
     if( ( v_P_k[ j ].second == v_Y_plus[ i ].second ) &&
         ( v_Y_plus[ i ].first <= t + 1 ) )
      for( Index s : psi_index( t , i ) )
       vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                       v_psi[ s ] ) );

    vars.push_back( std::make_pair( &v_active_power_k[ j ] , -1.0 ) );

    if( ( AR & FormMsk ) == SUSDForm ) {
     put_row( MaxPower_Const , v_P_h.size() + j , std::move( vars ) ,
              0.0 , Inf< double >() );
    } else {
     put_row( MaxPower_Const , j , std::move( vars ) , 0.0 , Inf< double >() );
    }
   }
  }
 }

 add_rows( MaxPower_Const , "MaxPower_Const_Thermal" );

 if( reserve_vars & 1u )  // if UCBlock has primary demand variables
  if( ! v_PrimaryRho.empty() ) {  // if unit produces any primary reserve
   // initializing primary rho fraction constraints - - - - - - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   size_rows( PrimaryRho_Const , f_time_horizon );

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {
    vars.push_back( std::make_pair( & v_active_power[ t ] ,
				    v_PrimaryRho[ t ] ) );

    vars.push_back( std::make_pair( & v_primary_spinning_reserve[ t ] ,
				    -1.0 ) );

    put_row( PrimaryRho_Const , t , std::move( vars ) ,
             0.0 , Inf< double >() );
    }

   add_rows( PrimaryRho_Const ,
                          "PrimaryRho_Const_Thermal" );
   }

 if( reserve_vars & 2u )  // if UCBlock has secondary demand variables
  if( ! v_SecondaryRho.empty() ) {  // if unit produces any secondary reserve
   // initializing secondary rho fraction constraints - - - - - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   size_rows( SecondaryRho_Const , f_time_horizon );

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {
    vars.push_back( std::make_pair( & v_active_power[ t ] ,
				    v_SecondaryRho[ t ] ) );

    vars.push_back( std::make_pair( & v_secondary_spinning_reserve[ t ] ,
                                    -1.0 ) );

    put_row( SecondaryRho_Const , t , std::move( vars ) ,
             0.0 , Inf< double >() );
   }

   add_rows( SecondaryRho_Const ,
			  "SecondaryRho_Const_Thermal" );
   }

 // spinning-reserve band constraints (deliverability model)- - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // The reserve r_t = pr_t + sr_t must fit the boundary capacity band around
 // p_t AND be deliverable within one ramp step from the previous realised
 // output (the "residual ramp"). These rows are uniform across all
 // formulations -- they use only the commitment u, the start-up v and the
 // shut-down w variables, present in every formulation -- so the reserve model
 // is formulation-independent and coincides with the one the DP solvers
 // implement. See the "Spinning reserve constraints" part of this method's
 // documentation. The fraction caps pr <= rho^p p, sr <= rho^s p above stay
 // separate; here r_t denotes pr_t + sr_t and every row has reserve coeff +1.
 {
  const bool has_pr = ( reserve_vars & 1u ) && ( ! v_PrimaryRho.empty() );
  const bool has_sr = ( reserve_vars & 2u ) && ( ! v_SecondaryRho.empty() );
  if( has_pr || has_sr ) {
   const bool has_ru = ! v_DeltaRampUp.empty();
   const bool has_rd = ! v_DeltaRampDown.empty();
   // the transitions t - 1 -> t, including -1 -> 0 if the unit is on before
   // the horizon (p_{-1} = InitialPower, u_{-1} = 1)
   const Index t_res = ( f_InitUpDownTime > 0 ) ? 0 : 1;
   const Index ntrans = ( f_time_horizon > t_res ) ?
                        f_time_horizon - t_res : 0;
   size_rows( Reserve_Const , 3 * f_time_horizon +
                         ( has_ru ? ntrans : 0 ) + ( has_rd ? ntrans : 0 ) );

   // start-up (v) / shut-down (w) indicator at t, or nullptr outside the free
   // window [ init_t , time_horizon ) where those variables live
   auto v_su = [ & ]( Index t ) -> ColVariable * {
    if( v_start_up.empty() || ( t < init_t ) ) return( nullptr );
    const Index j = t - init_t;
    return( j < v_start_up.size() ? &v_start_up[ j ] : nullptr );
    };
   auto v_sd = [ & ]( Index t ) -> ColVariable * {
    if( v_shut_down.empty() || ( t < init_t ) ) return( nullptr );
    const Index j = t - init_t;
    return( j < v_shut_down.size() ? &v_shut_down[ j ] : nullptr );
    };
   // append the reserve variable(s) with unit coefficient to the current row
   auto add_res = [ & ]( Index t ) {
    if( has_pr )
     vars.push_back( std::make_pair( &v_primary_spinning_reserve[ t ] , 1.0 ) );
    if( has_sr )
     vars.push_back( std::make_pair( &v_secondary_spinning_reserve[ t ] ,
                                     1.0 ) );
    };

   Index ci = 0;

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {
    const double pmin = get_operational_min_power( t );
    const double pmax = get_operational_max_power( t );

    // (1) foot-room: p_t - r_t >= pmin u_t, i.e. pmin u_t - p_t + r_t <= 0
    vars.push_back( std::make_pair( &v_commitment[ t ] , pmin ) );
    vars.push_back( std::make_pair( &v_active_power[ t ] , -1.0 ) );
    add_res( t );
    put_row( Reserve_Const , ci++ , std::move( vars ) ,
             -Inf< double >() , 0.0 );

    // (2) head-room, start-up side:
    //     p_t + r_t <= pmax u_t + ( SUlim_t - pmax ) v_t
    vars.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );
    vars.push_back( std::make_pair( &v_commitment[ t ] , -pmax ) );
    if( auto * v = v_su( t ) )
     vars.push_back( std::make_pair( v , -( v_StartUpLimit[ t ] - pmax ) ) );
    add_res( t );
    put_row( Reserve_Const , ci++ , std::move( vars ) ,
             -Inf< double >() , 0.0 );

    // (3) head-room, shut-down side:
    //     p_t + r_t <= pmax u_t + ( SDlim_{t+1} - pmax ) w_{t+1}
    // (w_{t+1} = 1 marks the unit off at t+1, i.e. the shut-down event is at
    // t+1, so the cap is the shut-down limit indexed at t+1 -- the same index
    // the T formulation and the DP solvers use, not SDlim_t)
    vars.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );
    vars.push_back( std::make_pair( &v_commitment[ t ] , -pmax ) );
    if( auto * w = v_sd( t + 1 ) )
     vars.push_back( std::make_pair( w , -( v_ShutDownLimit[ t + 1 ] - pmax ) ) );
    add_res( t );
    put_row( Reserve_Const , ci++ , std::move( vars ) ,
             -Inf< double >() , 0.0 );
    }

   // (4) ramp-up deliverability (interior residual; relaxed at a start-up,
   //     where the head-room row already caps at SUlim):
   //     p_t + r_t - p_{t-1} + ( pmax_t - DRU_t ) u_{t-1} <= pmax_t
   //     and at t = 0, if the unit is on before the horizon,
   //     p_0 + r_0 <= InitialPower + DRU_0
   if( has_ru )
    for( Index t = t_res ; t < f_time_horizon ; ++t ) {
     const double pmax = get_operational_max_power( t );
     vars.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );
     if( t ) {
      vars.push_back( std::make_pair( &v_active_power[ t - 1 ] , -1.0 ) );
      vars.push_back( std::make_pair( &v_commitment[ t - 1 ] ,
                                      pmax - v_DeltaRampUp[ t ] ) );
      }
     add_res( t );
     put_row( Reserve_Const , ci++ , std::move( vars ) ,
              -Inf< double >() ,
              t ? pmax : f_InitialPower + v_DeltaRampUp[ 0 ] );
     }

   // (5) ramp-down deliverability:
   //     p_{t-1} - p_t + r_t + ( pmax_{t-1} - DRD_t ) u_t <= pmax_{t-1}
   //     and at t = 0, if the unit is on before the horizon, the same with
   //     p_{-1} = InitialPower in place of pmax_{-1} as well, so that a
   //     shut-down at 0 leaves the row slack:
   //     - p_0 + r_0 + ( InitialPower - DRD_0 ) u_0 <= 0
   //     [see update_initial_power_in_cnstrs()]
   if( has_rd )
    for( Index t = t_res ; t < f_time_horizon ; ++t ) {
     const double pmaxm = t ? get_operational_max_power( t - 1 )
                            : f_InitialPower;
     if( t )
      vars.push_back( std::make_pair( &v_active_power[ t - 1 ] , 1.0 ) );
     vars.push_back( std::make_pair( &v_active_power[ t ] , -1.0 ) );
     vars.push_back( std::make_pair( &v_commitment[ t ] ,
                                     pmaxm - v_DeltaRampDown[ t ] ) );
     add_res( t );
     put_row( Reserve_Const , ci++ , std::move( vars ) ,
              -Inf< double >() , t ? pmaxm : 0.0 );
     }

   add_rows( Reserve_Const , "Reserve_Const_Thermal" );
   }
  }

 // ZOConstraints - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( generate_ZOConstraints ) {
  // the commitment bound constraints
  Commitment_bound_Const.resize( f_time_horizon );

  for( Index t = 0 ; t < f_time_horizon ; ++t )
   Commitment_bound_Const[ t ].set_variable( & v_commitment[ t ] );

  add_rows( Commitment_bound_Const ,
			 "Commitment_bound_Thermal" );

  auto startup_shutdown_size = f_time_horizon - init_t;

  // the startup binary bound constraints
  StartUp_Binary_bound_Const.resize( startup_shutdown_size );

  for( Index t = 0 ; t < startup_shutdown_size ; ++t )
   StartUp_Binary_bound_Const[ t ].set_variable( & v_start_up[ t ] );

  add_rows( StartUp_Binary_bound_Const ,
                         "StartUp_binary_bound_Thermal" );

  // the shut-down binary bound constraints
  ShutDown_Binary_bound_Const.resize( startup_shutdown_size );

  for( Index t = 0 ; t < startup_shutdown_size ; ++t )
   ShutDown_Binary_bound_Const[ t ].set_variable( & v_shut_down[ t ] );

  add_rows( ShutDown_Binary_bound_Const ,
                         "ShoutDown_binary_bound_Thermal" );
  }

 // BoxConstraint - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( generating_rows() && ( init_t > 0 ) && ( f_InitUpDownTime > 0 ) &&
     ( f_InitUpDownTime < f_MinUpTime ) ) {

  // the commitment fixed to one BoxConstraints
  Commitment_fixed_to_One_Const.resize( init_t );

  for( Index t = 0 ; t < init_t ; ++t ) {
   Commitment_fixed_to_One_Const[ t ].set_both( 1 );
   Commitment_fixed_to_One_Const[ t ].set_variable( &v_commitment[ t ] );
  }

  add_rows( Commitment_fixed_to_One_Const ,
                         "Commitment_fixed_to_one_Thermal" );
 }

 if( AR & PCuts ) {

  // Initial perspective cuts constraints - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  auto cnstr_idx = 0;

  if( ( ( AR & FormMsk ) == tbinForm ) ||  // 3bin formulation- - - - - - - -
      ( ( AR & FormMsk ) == TForm ) ) {  // T formulation - - - - - - - - - -

   size_rows( Init_PC_Const , 2 * f_time_horizon );

   for( Index t = 0 ; t < f_time_horizon ; ++t )
    for( Index k = 0 ; k <= 1 ; ++k ) {

     auto value = ( k == 0 ? get_operational_min_power( t )
                           : get_operational_max_power( t ) );

     vars.push_back( std::make_pair( &v_active_power[ t ] , 2 * value ) );
     vars.push_back( std::make_pair( &v_cut[ t ] , -1.0 ) );
     vars.push_back( std::make_pair( &v_commitment[ t ] ,
                                     -std::pow( value , 2 ) ) );

     put_row( Init_PC_Const , cnstr_idx , std::move( vars ) ,
              -Inf< double >() , 0.0 );

     cnstr_idx++;
    }

  } else if( ( AR & FormMsk ) == ptForm ) {  // pt formulation- - - - - - - -

   size_rows( Init_PC_Const , 2 * f_time_horizon );

   for( Index t = 0 ; t < f_time_horizon ; ++t )
    for( Index k = 0 ; k <= 1 ; ++k ) {

     auto value = ( k == 0 ? get_operational_min_power( t )
                           : get_operational_max_power( t ) );

     vars.push_back( std::make_pair( &v_active_power[ t ] , 2 * value ) );
     vars.push_back( std::make_pair( &v_cut[ t ] , -1.0 ) );

     for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
      if( ( v_Y_plus[ i ].first <= t + 1 ) &&
          ( t + 1 <= v_Y_plus[ i ].second ) )
       vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                       -std::pow( value , 2 ) ) );

     put_row( Init_PC_Const , cnstr_idx , std::move( vars ) ,
              -Inf< double >() , 0.0 );

     cnstr_idx++;
    }

  } else if( ( AR & FormMsk ) == DPForm ) {  // DP formulation- - - - - - - -

   size_rows( Init_PC_Const , 2 * v_P_h_k.size() );

   for( Index j = 0 ; j < v_P_h_k.size() ; ++j )
    for( Index k = 0 ; k <= 1 ; ++k ) {

     auto t = v_P_h_k[ j ].first;
     auto value = ( k == 0 ? get_operational_min_power( t )
                           : get_operational_max_power( t ) );

     vars.push_back( std::make_pair( &v_active_power_h_k[ j ] , 2 * value ) );

     for( Index s = 0 ; s < v_P_h_k.size() ; ++s )
      if( ( v_P_h_k[ j ].second.first == v_Z_h_k[ s ].second.first ) &&
          ( v_P_h_k[ j ].second.second == v_Z_h_k[ s ].second.second ) )
       if( v_Z_h_k[ s ].first == t )
        vars.push_back( std::make_pair( &v_cut_h_k[ s ] , -1.0 ) );

     for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
      if( ( v_Y_plus[ i ].first == v_P_h_k[ j ].second.first ) &&
          ( v_Y_plus[ i ].second == v_P_h_k[ j ].second.second ) )
       vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                       -std::pow( value , 2 ) ) );

     put_row( Init_PC_Const , cnstr_idx , std::move( vars ) ,
              -Inf< double >() , 0.0 );

     cnstr_idx++;
    }
  }

  if( ( AR & FormMsk ) == SUForm || ( AR & FormMsk ) == SUSDForm )
  {  // SU and SUSD formulation - - - - - - - - - - - - - - - - - - - - - - -

   if( ( AR & FormMsk ) == SUSDForm )
    size_rows( Init_PC_Const , 2 * v_P_h.size() + 2 * v_P_k.size() );
   else
    size_rows( Init_PC_Const , 2 * v_P_h.size() );

   for( Index j = 0 ; j < v_P_h.size() ; ++j )
    for( Index k = 0 ; k <= 1 ; ++k ) {

     auto t = v_P_h[ j ].first;
     auto value = ( k == 0 ? get_operational_min_power( t )
                           : get_operational_max_power( t ) );

     vars.push_back( std::make_pair( &v_active_power_h[ j ] , 2 * value ) );

     for( Index s = 0 ; s < v_P_h.size() ; ++s )
      if( ( v_P_h[ j ].second == v_Z_h[ s ].second ) &&
          ( v_P_h[ j ].first == v_Z_h[ s ].first ) )
       vars.push_back( std::make_pair( &v_cut_h[ s ] , -1.0 ) );

     for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
      if( v_P_h[ j ].second == v_Y_plus[ i ].first )
       if( t + 1 <= v_Y_plus[ i ].second )
        vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                        -std::pow( value , 2 ) ) );

     put_row( Init_PC_Const , cnstr_idx , std::move( vars ) ,
              -Inf< double >() , 0.0 );

     cnstr_idx++;
    }
  }

  if( ( AR & FormMsk ) == SDForm || ( AR & FormMsk ) == SUSDForm )
  {  // SD and SUSD formulation - - - - - - - - - - - - - - - - - - - - - - -

   if( ( AR & FormMsk ) == SDForm )
    size_rows( Init_PC_Const , 2 * v_P_k.size() );

   for( Index j = 0 ; j < v_P_k.size() ; ++j )
    for( Index k = 0 ; k <= 1 ; ++k ) {

     auto t = v_P_k[ j ].first;
     auto value = ( k == 0 ? get_operational_min_power( t )
                           : get_operational_max_power( t ) );

     vars.push_back( std::make_pair( &v_active_power_k[ j ] , 2 * value ) );

     for( Index s = 0 ; s < v_P_k.size() ; ++s )
      if( ( v_P_k[ j ].second == v_Z_k[ s ].second ) &&
          ( v_P_k[ j ].first == v_Z_k[ s ].first ) )
       vars.push_back( std::make_pair( &v_cut_k[ s ] , -1.0 ) );

     for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
      if( v_P_k[ j ].second == v_Y_plus[ i ].second )
       if( v_Y_plus[ i ].first <= t + 1 )
        vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                        -std::pow( value , 2 ) ) );

     put_row( Init_PC_Const , cnstr_idx , std::move( vars ) ,
              -Inf< double >() , 0.0 );

     cnstr_idx++;
    }
  }

  add_rows( Init_PC_Const , "Init_PC_Const_Thermal" );

  // Constraints connecting variables of the SUSD formulations with the
  // maximum of the perspective function of the SU and the SD formulations- -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( AR & FormMsk ) == SUSDForm ) {

   size_rows( Max_SUSD_PC_Const , 2 * f_time_horizon );

   if( ( AR & FormMsk ) == SUSDForm )
    for( Index t = 0 ; t < f_time_horizon ; ++t ) {

     vars.push_back( std::make_pair( &v_cut_teta[ t ] , -1.0 ) );

     for( Index j = 0 ; j < v_Z_h.size() ; ++j )
      if( v_Z_h[ j ].first == t )
       vars.push_back( std::make_pair( &v_cut_h[ j ] , 1.0 ) );

     put_row( Max_SUSD_PC_Const , t , std::move( vars ) ,
              -Inf< double >() , 0.0 );
    }

   if( ( AR & FormMsk ) == SUSDForm )
    for( Index t = 0 ; t < f_time_horizon ; ++t ) {

     vars.push_back( std::make_pair( &v_cut_teta[ t ] , -1.0 ) );

     for( Index j = 0 ; j < v_Z_k.size() ; ++j )
      if( v_Z_k[ j ].first == t )
       vars.push_back( std::make_pair( &v_cut_k[ j ] , 1.0 ) );

     put_row( Max_SUSD_PC_Const , f_time_horizon + t , std::move( vars ) ,
              -Inf< double >() , 0.0 );
    }
   add_rows( Max_SUSD_PC_Const , "Max_SUSD_PC_Const_Thermal" );
  }

  // Constraints connecting perspective cuts variables of 3bin with those of
  // DP, SU and SD formulations- - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( ( AR & FormMsk ) == tbinForm ) ||  // 3bin formulation- - - - - - - -
      ( ( AR & FormMsk ) == TForm ) ||  // T formulation- - - - - - - - - - -
      ( ( AR & FormMsk ) == ptForm ) ) {  // pt formulation - - - - - - - - -
   ;  // does nothing
  } else {  // DP, SU, SD and SUSD formulations - - - - - - - - - - - - - - -
   if( ( AR & FormMsk ) == SUSDForm )
    size_rows( Eq_PC_Const , 2 * f_time_horizon );
   else
    size_rows( Eq_PC_Const , f_time_horizon );

   if( ( AR & FormMsk ) == DPForm ) {  // DP formulation- - - - - - - - - - -

    for( Index t = 0 ; t < f_time_horizon ; ++t ) {

     vars.push_back( std::make_pair( &v_cut[ t ] , 1.0 ) );

     for( Index j = 0 ; j < v_Z_h_k.size() ; ++j )
      if( v_Z_h_k[ j ].first == t )
       vars.push_back( std::make_pair( &v_cut_h_k[ j ] , -1.0 ) );

     put_row( Eq_PC_Const , t , std::move( vars ) , 0.0 , 0.0 );
    }
   }

   if( ( AR & FormMsk ) == SUForm || ( AR & FormMsk ) == SUSDForm )
   {  // SU and SUSD formulation- - - - - - - - - - - - - - - - - - - - - - -

    for( Index t = 0 ; t < f_time_horizon ; ++t ) {

     vars.push_back( std::make_pair( &v_cut[ t ] , 1.0 ) );

     for( Index j = 0 ; j < v_Z_h.size() ; ++j )
      if( v_Z_h[ j ].first == t )
       vars.push_back( std::make_pair( &v_cut_h[ j ] , -1.0 ) );

     put_row( Eq_PC_Const , t , std::move( vars ) , 0.0 , 0.0 );
    }
   }

   if( ( AR & FormMsk ) == SDForm || ( AR & FormMsk ) == SUSDForm )
   {  // SU and SUSD formulation- - - - - - - - - - - - - - - - - - - - - - -

    for( Index t = 0 ; t < f_time_horizon ; ++t ) {

     vars.push_back( std::make_pair( &v_cut[ t ] , 1.0 ) );

     for( Index j = 0 ; j < v_Z_k.size() ; ++j )
      if( v_Z_k[ j ].first == t )
       vars.push_back( std::make_pair( &v_cut_k[ j ] , -1.0 ) );

     if( ( AR & FormMsk ) == SUSDForm ) {
      put_row( Eq_PC_Const , f_time_horizon + t , std::move( vars ) ,
               0.0 , 0.0 );
     } else {
      put_row( Eq_PC_Const , t , std::move( vars ) , 0.0 , 0.0 );
     }
    }
   }
  }

  add_rows( Eq_PC_Const , "Eq_PC_Const_Thermal" );
  }

 if ( ! v_RefSchedule.empty() ) {
  size_rows( Reference_Schedule_Const , 2 * f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   // | P - Pref | <= v_abs_ref_schedule
   vars.push_back( std::make_pair( & v_active_power[ t ] , 1.0 ) );
   vars.push_back( std::make_pair( & v_abs_ref_schedule[ t ] , -1.0 ) );
   put_row( Reference_Schedule_Const , t , std::move( vars ) ,
            -Inf< double >() , v_RefSchedule[ t ] );
   vars.push_back( std::make_pair( & v_active_power[ t ] , -1.0 ) );
   vars.push_back( std::make_pair( & v_abs_ref_schedule[ t ] , -1.0 ) );
   put_row( Reference_Schedule_Const , f_time_horizon + t ,
            std::move( vars ) , -Inf< double >() , - v_RefSchedule[ t ] );
   }

  add_rows( Reference_Schedule_Const ,
			 "Norm1_Reference_Schedule" );
  }

 // reactive power bounds constraints (if any) - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // the bound is Qmin_off + Qmin_on u <= q <= Qmax_off + Qmax_on u; when the
 // commitment coefficients Qmin/max_on are all zero (the default) it is the
 // plain box on q and a BoxConstraint suffices, otherwise q is coupled to the
 // commitment u and two FRowConstraints are needed. If the unit has a design
 // variable x, an unbuilt unit has no reactive power either: the off bounds
 // Qmin/max_off then multiply x (when finite), which also needs the rows
 const bool reactive_design = ( f_InvestmentCost != 0 );
 const bool reactive_gated = ( ! v_MinReactivePowerOn.empty() ) ||
                             ( ! v_MaxReactivePowerOn.empty() ) ||
                             reactive_design;

 // an absent bound is 0, as the getters and the DP solvers take it
 if( f_reactive_power ) {

  if( ! reactive_gated ) {
   // a box, whose data nothing changes after the generation
   if( generating_rows() ) {
    ReactivePower_Bound_Const.resize( f_time_horizon );

    for( Index t = 0 ; t < f_time_horizon ; ++t ) {
     ReactivePower_Bound_Const[ t ].set_rhs( get_max_reactive_power( t ) );
     ReactivePower_Bound_Const[ t ].set_lhs( get_min_reactive_power( t ) );
     ReactivePower_Bound_Const[ t ].set_variable( & v_reactive_power[ t ] );
     }

    add_rows( ReactivePower_Bound_Const , "ReactivePowerBound_thermal" );
    }
   }
  else {
   size_rows( ReactivePowerMax_Const , f_time_horizon );
   size_rows( ReactivePowerMin_Const , f_time_horizon );

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {
    // q[t] - Qmax_on[t] u[t] - Qmax_off[t] x <= 0, or <= Qmax_off[t]
    const double qmax = get_max_reactive_power( t );
    const bool dmax = reactive_design && ( qmax < Inf< double >() );
    LinearFunction::v_coeff_pair vmax;
    vmax.push_back( std::make_pair( & v_reactive_power[ t ] , 1.0 ) );
    vmax.push_back( std::make_pair( & v_commitment[ t ] ,
                                    - get_max_reactive_power_on( t ) ) );
    if( dmax && ( qmax != 0 ) )
     vmax.push_back( std::make_pair( & design , - qmax ) );
    put_row( ReactivePowerMax_Const , t , std::move( vmax ) ,
             - Inf< double >() , dmax ? 0.0 : qmax );

    // q[t] - Qmin_on[t] u[t] - Qmin_off[t] x >= 0, or >= Qmin_off[t]
    const double qmin = get_min_reactive_power( t );
    const bool dmin = reactive_design && ( qmin > - Inf< double >() );
    LinearFunction::v_coeff_pair vmin;
    vmin.push_back( std::make_pair( & v_reactive_power[ t ] , 1.0 ) );
    vmin.push_back( std::make_pair( & v_commitment[ t ] ,
                                    - get_min_reactive_power_on( t ) ) );
    if( dmin && ( qmin != 0 ) )
     vmin.push_back( std::make_pair( & design , - qmin ) );
    put_row( ReactivePowerMin_Const , t , std::move( vmin ) ,
             dmin ? 0.0 : qmin , Inf< double >() );
    }

   add_rows( ReactivePowerMax_Const ,
                          "ReactivePowerMax_thermal" );
   add_rows( ReactivePowerMin_Const ,
                          "ReactivePowerMin_thermal" );
   }
  }

 // the box of the active power, 0 <= p[ t ] <= the operational maximum
 // power: the rows that tie the power to the commitment imply it, but a
 // Solver that only reads the boxes (say, a BoxSolver bounding the
 // Objective) needs it to see the power bounded
 if( generating_rows() )
  ActivePower_Bound_Const.resize( f_time_horizon );

 for( Index t = 0 ; t < f_time_horizon ; ++t )
  put_box( ActivePower_Bound_Const[ t ] , & v_active_power[ t ] , 0 ,
           get_operational_max_power( t ) );

 add_rows( ActivePower_Bound_Const , "ActivePowerBound_thermal" );

 }  // end( ThermalUnitBlock::build_rows )

/*--------------------------------------------------------------------------*/
/* What update_rows() collects while build_rows() compares the rows: the
 * changes of the static rows and of the boxes, the number of rows written
 * by push_row() in each group, and the rows written in each dynamic
 * group. */

struct ThermalUnitBlock::RowCmp {
 struct RowChange {                   // a static row that changes
  FRowConstraint * row;
  Subset idx;                         // the positions of the coefficients
  std::vector< double > coef;         // ... and their new values
  double lhs;
  double rhs;
  };
 struct BoxChange {                   // a box that changes
  BoxConstraint * box;
  double lhs;
  double rhs;
  };
 struct DynRow {                      // a row written in a dynamic group
  std::vector< int > key;
  LinearFunction::v_coeff_pair vars;
  double lhs;
  double rhs;
  };

 std::vector< RowChange > rows;
 std::vector< BoxChange > boxes;
 std::map< std::vector< FRowConstraint > * , Index > pushed;
 std::set< std::vector< FRowConstraint > * > written;  // groups written
 std::map< DynRows * , std::vector< DynRow > > dyn;
 };

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::size_rows( std::vector< FRowConstraint > & rows ,
                                  Index n )
{
 if( generating_rows() ) {
  register_row_group( rows );
  rows.resize( n );
  }
 else {
  f_row_cmp->written.insert( & rows );
  if( rows.size() != n )
   throw( std::logic_error( "ThermalUnitBlock::update_rows: a group of " +
                            std::to_string( rows.size() ) + " rows would "
                            "have " + std::to_string( n ) + " of them" ) );
  }
 }

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::register_row_group(
                                     std::vector< FRowConstraint > & rows )
{
 if( std::find( f_row_groups.begin() , f_row_groups.end() , & rows ) ==
     f_row_groups.end() )
  f_row_groups.push_back( & rows );
 }

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::put_row( std::vector< FRowConstraint > & rows ,
                                Index i ,
                                LinearFunction::v_coeff_pair && vars ,
                                double lhs , double rhs )
{
 if( generating_rows() ) {
  rows[ i ].set_lhs( lhs );
  rows[ i ].set_rhs( rhs );
  rows[ i ].set_function( new LinearFunction( std::move( vars ) ) );
  return;
  }

 if( i >= rows.size() )
  throw( std::logic_error( "ThermalUnitBlock::update_rows: a group of " +
                           std::to_string( rows.size() ) + " rows would "
                           "have more of them" ) );

 auto & row = rows[ i ];
 const auto & old = static_cast< const LinearFunction * >(
                                         row.get_function() )->get_v_var();
 if( old.size() != vars.size() )
  throw( std::logic_error( "ThermalUnitBlock::update_rows: a row would have "
                           "another number of Variable" ) );

 RowCmp::RowChange c{ & row , {} , {} , lhs , rhs };
 for( Index k = 0 ; k < vars.size() ; ++k ) {
  if( old[ k ].first != vars[ k ].first )
   throw( std::logic_error( "ThermalUnitBlock::update_rows: a row would "
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

void ThermalUnitBlock::push_row( std::vector< FRowConstraint > & rows ,
                                 LinearFunction::v_coeff_pair && vars ,
                                 double lhs , double rhs )
{
 if( generating_rows() ) {
  register_row_group( rows );
  rows.emplace_back();
  put_row( rows , rows.size() - 1 , std::move( vars ) , lhs , rhs );
  }
 else {
  f_row_cmp->written.insert( & rows );
  put_row( rows , f_row_cmp->pushed[ & rows ]++ , std::move( vars ) , lhs ,
           rhs );
  }
 }

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::put_dyn_row( DynRows & rows ,
                                    std::vector< int > && key ,
                                    LinearFunction::v_coeff_pair && vars ,
                                    double lhs , double rhs )
{
 if( generating_rows() ) {
  rows.rows.emplace_back();
  rows.rows.back().set_lhs( lhs );
  rows.rows.back().set_rhs( rhs );
  rows.rows.back().set_function( new LinearFunction( std::move( vars ) ) );
  rows.keys.push_back( std::move( key ) );
  }
 else {
  if( std::find( f_dyn_groups.begin() , f_dyn_groups.end() , & rows ) ==
      f_dyn_groups.end() )
   throw( std::logic_error( "ThermalUnitBlock::update_rows: a dynamic group "
                            "that is not generated would have rows" ) );
  f_row_cmp->dyn[ & rows ].push_back( { std::move( key ) , std::move( vars ) ,
                                        lhs , rhs } );
  }
 }

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::put_box( BoxConstraint & box , ColVariable * var ,
                                double lhs , double rhs )
{
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

void ThermalUnitBlock::update_rows( ModParam issueAMod )
{
 if( ! constraints_generated() )
  return;

 // write the rows anew, comparing them with those there are: nothing is
 // changed until all of them are written, so that a throw changes nothing
 RowCmp cmp;
 f_row_cmp = & cmp;
 try {
  build_rows( false );
  for( const auto & p : cmp.pushed )
   if( p.second != p.first->size() )
    throw( std::logic_error( "ThermalUnitBlock::update_rows: a group of " +
                             std::to_string( p.first->size() ) + " rows "
                             "would have " + std::to_string( p.second ) +
                             " of them" ) );
  // a group of rows none of which is written any longer
  for( auto * g : f_row_groups )
   if( ( ! g->empty() ) && ( ! cmp.written.count( g ) ) )
    throw( std::logic_error( "ThermalUnitBlock::update_rows: a group of " +
                             std::to_string( g->size() ) + " rows would "
                             "have none of them" ) );
  }
 catch( ... ) {
  f_row_cmp = nullptr;
  throw;
  }
 f_row_cmp = nullptr;

 // the dynamic rows: those of a key that is in both and has the same
 // Variable are compared, the others are removed or added
 using RowIt = std::list< FRowConstraint >::iterator;
 using KeyIt = std::list< std::vector< int > >::iterator;
 struct DynChange {
  DynRows * rows;
  std::vector< RowIt > rmv;                 // the rows to be removed
  std::vector< KeyIt > rmv_key;             // ... and their keys
  std::vector< RowCmp::DynRow * > add;      // the rows to be added
  };
 std::vector< DynChange > dyn;

 for( auto & d : cmp.dyn ) {
  DynChange dc{ d.first , {} , {} , {} };
  std::map< std::vector< int > , std::pair< RowIt , KeyIt > > there;
  auto kit = d.first->keys.begin();
  for( auto rit = d.first->rows.begin() ; rit != d.first->rows.end() ;
       ++rit , ++kit )
   there.emplace( * kit , std::make_pair( rit , kit ) );

  for( auto & r : d.second ) {
   auto it = there.find( r.key );
   bool same = false;
   if( it != there.end() ) {
    const auto & old = static_cast< const LinearFunction * >(
                          it->second.first->get_function() )->get_v_var();
    same = ( old.size() == r.vars.size() );
    for( Index k = 0 ; same && ( k < old.size() ) ; ++k )
     same = ( old[ k ].first == r.vars[ k ].first );
    if( same ) {
     RowCmp::RowChange c{ & * it->second.first , {} , {} , r.lhs , r.rhs };
     for( Index k = 0 ; k < old.size() ; ++k )
      if( old[ k ].second != r.vars[ k ].second ) {
       c.idx.push_back( k );
       c.coef.push_back( r.vars[ k ].second );
       }
     if( ( ! c.idx.empty() ) || ( it->second.first->get_lhs() != r.lhs ) ||
         ( it->second.first->get_rhs() != r.rhs ) )
      cmp.rows.push_back( std::move( c ) );
     }
    else {
     dc.rmv.push_back( it->second.first );
     dc.rmv_key.push_back( it->second.second );
     }
    there.erase( it );
    }
   if( ! same )
    dc.add.push_back( & r );
   }

  for( auto & t : there ) {   // the rows of the keys no longer there
   dc.rmv.push_back( t.second.first );
   dc.rmv_key.push_back( t.second.second );
   }

  if( ( ! dc.rmv.empty() ) || ( ! dc.add.empty() ) )
   dyn.push_back( std::move( dc ) );
  }

 // a dynamic group none of whose rows is written any longer
 for( auto * d : f_dyn_groups )
  if( ( ! d->rows.empty() ) && ( ! cmp.dyn.count( d ) ) ) {
   DynChange dc{ d , {} , {} , {} };
   auto kit = d->keys.begin();
   for( auto rit = d->rows.begin() ; rit != d->rows.end() ; ++rit , ++kit ) {
    dc.rmv.push_back( rit );
    dc.rmv_key.push_back( kit );
    }
   dyn.push_back( std::move( dc ) );
   }

 if( cmp.rows.empty() && cmp.boxes.empty() && dyn.empty() )
  return;  // nothing changes

 // now change: all the Modification in a single GroupModification
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

 for( auto & d : dyn ) {
  if( ! d.rmv.empty() ) {
   // remove_dynamic_constraints() wants the rows in the order of the list
   std::map< const FRowConstraint * , Index > pos;
   Index i = 0;
   for( const auto & r : d.rows->rows )
    pos[ & r ] = i++;
   std::vector< Index > ord( d.rmv.size() );
   std::iota( ord.begin() , ord.end() , 0 );
   std::sort( ord.begin() , ord.end() , [ & ]( Index a , Index b ) {
    return( pos[ & * d.rmv[ a ] ] < pos[ & * d.rmv[ b ] ] ); } );
   std::vector< RowIt > rmv( d.rmv.size() );
   for( i = 0 ; i < ord.size() ; ++i )
    rmv[ i ] = d.rmv[ ord[ i ] ];
   for( auto kit : d.rmv_key )
    d.rows->keys.erase( kit );
   remove_dynamic_constraints( d.rows->rows , rmv , nAM );
   }
  if( ! d.add.empty() ) {
   std::list< FRowConstraint > nr;
   for( auto * r : d.add ) {
    nr.emplace_back();
    nr.back().set_lhs( r->lhs );
    nr.back().set_rhs( r->rhs );
    nr.back().set_function( new LinearFunction( std::move( r->vars ) ) );
    d.rows->keys.push_back( std::move( r->key ) );
    }
   add_dynamic_constraints( d.rows->rows , nr , nAM );
   }
  }

 close_channel( par2chnl( nAM ) );  // at the end close the channel

 }  // end( ThermalUnitBlock::update_rows )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::generate_dynamic_constraints( Configuration * dycc )
{
 if( AR & PCuts ) {
  double tol = 1e-7;  // threshold parameter for P/C separation
  double eps = 1e-6;  // tolerance value to consider a binary variable

  bool check_loop = false;
  // the violation of the cut at the point, p^2 / u - z = 2 p^2 / u - psum
  // with psum = z + p^2 / u, is compared with tol times max{ 1 , psum }
  double part1 , psum;

  auto extract_parameters = [ & tol , & eps ]( Configuration * c )
   -> bool {
   if( auto tc = dynamic_cast< SimpleConfiguration< double > * >( c ) ) {
    tol = tc->f_value;
    return( true );
   }
   if( auto tc = dynamic_cast<
    SimpleConfiguration< std::pair< double , double > > * >( c ) ) {
    tol = tc->f_value.first;
    eps = tc->f_value.second;
    return( true );
   }
   return( false );
  };

  if( ( ! extract_parameters( dycc ) ) && f_BlockConfig )
   // if the given Configuration is not valid, try the one from the BlockConfig
   extract_parameters( f_BlockConfig->f_dynamic_constraints_Configuration );

  LinearFunction::v_coeff_pair vars;

  // tbin and T formulations- - - - - - - - - - - - - - - - - - - - - - - - -
  if( ( AR & FormMsk ) == tbinForm || ( AR & FormMsk ) == TForm ) {
   bool check_cut = false;
   double pbar = 0;  // read only by the check_loop variant

   for( Index t = 0 ; t < f_time_horizon ; ++t )

    if( v_commitment[ t ].get_value() > eps ) {

     check_cut = false;

     if( check_loop ) {
      pbar = v_active_power[ t ].get_value() / v_commitment[ t ].get_value();
      double value = ( v_active_power[ t ].get_value() *
       v_active_power[ t ].get_value() ) / v_commitment[ t ].get_value();
      if( v_cut[ t ].get_value() < value - tol )
       if( prevpbar[ t ] == 0 || ( prevpbar[ t ] != 0 && std::abs(
         ( prevpbar[ t ] - pbar ) / prevpbar[ t ] ) > eps ) )
        check_cut = true;
     } else {
      psum = v_cut[ t ].get_value() + std::pow(
             v_active_power[ t ].get_value() /
               v_commitment[ t ].get_value() , 2 ) * v_commitment[ t ].get_value();
      part1 = std::max( 1.0 , psum );

      if( 2.0 * ( v_active_power[ t ].get_value() / v_commitment[ t ].get_value() ) *
          v_active_power[ t ].get_value() - psum >= tol * part1 )
       check_cut = true;
     }

     if( check_cut ) {
      prevpbar[ t ] = pbar;

      std::list< FRowConstraint > cut( 1 );

      vars.push_back( std::make_pair( &v_active_power[ t ] ,
                                      2 * ( v_active_power[ t ].get_value() /
                                            v_commitment[ t ].get_value() ) ) );
      vars.push_back( std::make_pair( &v_cut[ t ] , -1.0 ) );
      vars.push_back(
       std::make_pair( &v_commitment[ t ] ,
                       -( std::pow( v_active_power[ t ].get_value() , 2 ) /
                          std::pow( v_commitment[ t ].get_value() , 2 ) ) ) );

      cut.front().set_lhs( -Inf< double >() );
      cut.front().set_rhs( 0.0 );
      cut.front().set_function(
       new LinearFunction( std::move( vars ) ) , eNoMod );

      add_dynamic_constraints( PC_cuts , cut , eNoBlck );
     }
    }
  }

  // pt formulation - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if( ( AR & FormMsk ) == ptForm ) {
   for( Index t = 0 ; t < f_time_horizon ; ++t ) {
    double sumy = 0.0;
    for( int j = 0 ; j < v_Y_plus.size() ; j++ )
     if( v_Y_plus[ j ].first <= t + 1 && t + 1 <= v_Y_plus[ j ].second ) {
      sumy += v_commitment_plus[ j ].get_value();
     }

    if( sumy > eps ) {
     psum = v_cut[ t ].get_value() +
       std::pow( v_active_power[ t ].get_value() / sumy , 2 ) * sumy;
     part1 = std::max( 1.0 , psum );

     if( 2.0 * ( v_active_power[ t ].get_value() / sumy ) *
                 v_active_power[ t ].get_value() - psum >= tol * part1 ) {
      std::list< FRowConstraint > cut( 1 );

      vars.push_back( std::make_pair( &v_active_power[ t ] ,
                                      2 * ( v_active_power[ t ].get_value() / sumy ) ) );
      vars.push_back( std::make_pair( &v_cut[ t ] , -1.0 ) );

      for( int j = 0 ; j < v_Y_plus.size() ; j++ )
       if( v_Y_plus[ j ].first <= t + 1 && t + 1 <= v_Y_plus[ j ].second ) {
        vars.push_back(
         std::make_pair( &v_commitment_plus[ j ] ,
                         -( std::pow( v_active_power[ t ].get_value() , 2 ) /
                            std::pow( sumy , 2 ) ) ) );
       }

      cut.front().set_lhs( -Inf< double >() );
      cut.front().set_rhs( 0.0 );
      cut.front().set_function(
       new LinearFunction( std::move( vars ) ) , eNoMod );

      add_dynamic_constraints( PC_cuts , cut , eNoBlck );
     }
    }
   }
  }

  // DP formulation - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if( ( AR & FormMsk ) == DPForm ) {
   for( Index i = 0 ; i < v_P_h_k.size() ; ++i ) {
    for( int j = 0 ; j < v_Y_plus.size() ; j++ )
     if( ( v_Y_plus[ j ].first == v_P_h_k[ i ].second.first ) &&
         ( v_P_h_k[ i ].second.second == v_Y_plus[ j ].second ) )
      if( v_commitment_plus[ j ].get_value() > eps )
       for( Index s = 0 ; s < v_P_h_k.size() ; ++s )
        if( ( v_P_h_k[ i ].first == v_Z_h_k[ s ].first ) &&
            ( v_P_h_k[ i ].second.first == v_Z_h_k[ s ].second.first ) &&
            ( v_P_h_k[ i ].second.second == v_Z_h_k[ s ].second.second ) ) {
         psum = v_cut_h_k[ s ].get_value() + std::pow(
            v_active_power_h_k[ i ].get_value() /
            v_commitment_plus[ j ].get_value() , 2 ) *
           v_commitment_plus[ j ].get_value();
         part1 = std::max( 1.0 , psum );

         if( 2.0 * ( v_active_power_h_k[ i ].get_value() /
             v_commitment_plus[ j ].get_value() ) *
             v_active_power_h_k[ i ].get_value() - psum >= tol * part1 ) {
          std::list< FRowConstraint > cut( 1 );

          vars.push_back(
           std::make_pair( &v_active_power_h_k[ i ] ,
            2 * ( v_active_power_h_k[ i ].get_value() /
                  v_commitment_plus[ j ].get_value() ) ) );
          vars.push_back( std::make_pair( &v_cut_h_k[ s ] , -1.0 ) );

          vars.push_back(
           std::make_pair(
            &v_commitment_plus[ j ] ,
            -( std::pow( v_active_power_h_k[ i ].get_value() , 2 ) /
               std::pow( v_commitment_plus[ j ].get_value() , 2 ) ) ) );

          cut.front().set_lhs( -Inf< double >() );
          cut.front().set_rhs( 0.0 );
          cut.front().set_function(
           new LinearFunction( std::move( vars ) ) , eNoMod );

          add_dynamic_constraints( PC_cuts , cut , eNoBlck );
         }
        }
   }
  }

  // SU, SD and SUSD formulations - - - - - - - - - - - - - - - - - - - - - -
  if( ( AR & FormMsk ) == SUForm || ( AR & FormMsk ) == SDForm ||
      ( AR & FormMsk ) == SUSDForm ) {
   if( ( AR & FormMsk ) == SUForm || ( AR & FormMsk ) == SUSDForm )
    for( Index i = 0 ; i < v_P_h.size() ; ++i ) {
     double sumy = 0.0;
     for( int j = 0 ; j < v_Y_plus.size() ; j++ )
      if( ( v_Y_plus[ j ].first == v_P_h[ i ].second ) &&
          ( v_P_h[ i ].first + 1 <= v_Y_plus[ j ].second ) )
       sumy += v_commitment_plus[ j ].get_value();

     if( sumy > eps )
      for( Index s = 0 ; s < v_P_h.size() ; ++s )
       if( ( v_P_h[ i ].first == v_Z_h[ s ].first ) &&
           ( v_P_h[ i ].second == v_Z_h[ s ].second ) )
        if( v_Z_h[ s ].first == v_P_h[ i ].first ) {
         psum = v_cut_h[ s ].get_value() +
           std::pow( v_active_power_h[ i ].get_value() / sumy , 2 ) * sumy;
         part1 = std::max( 1.0 , psum );

         if( 2.0 * ( v_active_power_h[ i ].get_value() / sumy ) *
             v_active_power_h[ i ].get_value() - psum >= tol * part1 ) {
          std::list< FRowConstraint > cut( 1 );

          vars.push_back(
           std::make_pair( &v_active_power_h[ i ] ,
                           2 * ( v_active_power_h[ i ].get_value() / sumy ) ) );
          vars.push_back( std::make_pair( &v_cut_h[ s ] , -1.0 ) );

          for( int j = 0 ; j < v_Y_plus.size() ; j++ )
           if( ( v_Y_plus[ j ].first == v_P_h[ i ].second ) &&
               ( v_P_h[ i ].first + 1 <= v_Y_plus[ j ].second ) ) {
            vars.push_back(
             std::make_pair( &v_commitment_plus[ j ] ,
                             -( std::pow( v_active_power_h[ i ].get_value() , 2 ) /
                                std::pow( sumy , 2 ) ) ) );
           }

          cut.front().set_lhs( -Inf< double >() );
          cut.front().set_rhs( 0.0 );
          cut.front().set_function(
           new LinearFunction( std::move( vars ) ) , eNoMod );

          add_dynamic_constraints( PC_cuts , cut , eNoBlck );
         }
        }
    }

   // SD and SUSD formulations- - - - - - - - - - - - - - - - - - - - - - - -
   if( ( AR & FormMsk ) == SDForm || ( AR & FormMsk ) == SUSDForm )
    for( Index i = 0 ; i < v_P_k.size() ; ++i ) {
     double sumy = 0.0;
     for( int j = 0 ; j < v_Y_plus.size() ; j++ )
      if( ( v_Y_plus[ j ].second == v_P_k[ i ].second ) &&
          ( v_P_k[ i ].first + 1 >= v_Y_plus[ j ].first ) )
       sumy += v_commitment_plus[ j ].get_value();

     if( sumy > eps )
      for( Index s = 0 ; s < v_P_k.size() ; ++s )
       if( ( v_P_k[ i ].first == v_Z_k[ s ].first ) &&
           ( v_P_k[ i ].second == v_Z_k[ s ].second ) )
        if( v_Z_k[ s ].first == v_P_k[ i ].first ) {
         psum = v_cut_k[ s ].get_value() +
           std::pow( v_active_power_k[ i ].get_value() / sumy , 2 ) * sumy;
         part1 = std::max( 1.0 , psum );

         if( 2.0 * ( v_active_power_k[ i ].get_value() / sumy ) *
          v_active_power_k[ i ].get_value() - psum >= tol * part1 ) {

          std::list< FRowConstraint > cut( 1 );

          vars.push_back(
           std::make_pair( &v_active_power_k[ i ] ,
                           2 * ( v_active_power_k[ i ].get_value() / sumy ) ) );
          vars.push_back( std::make_pair( &v_cut_k[ s ] , -1.0 ) );

          for( int j = 0 ; j < v_Y_plus.size() ; j++ )
           if( ( v_Y_plus[ j ].second == v_P_k[ i ].second ) &&
               ( v_P_k[ i ].first + 1 >= v_Y_plus[ j ].first ) ) {
            vars.push_back(
             std::make_pair( &v_commitment_plus[ j ] ,
                             -( std::pow( v_active_power_k[ i ].get_value() , 2 ) /
                                std::pow( sumy , 2 ) ) ) );
           }

          cut.front().set_lhs( -Inf< double >() );
          cut.front().set_rhs( 0.0 );
          cut.front().set_function(
           new LinearFunction( std::move( vars ) ) , eNoMod );

          add_dynamic_constraints( PC_cuts , cut , eNoBlck );
         }
        }
    }
  }

  add_dynamic_constraint( PC_cuts , "PC_cuts_Thermal" );
  }
 }  // end( ThermalUnitBlock::generate_dynamic_constraints )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::generate_objective( Configuration * objc )
{
 if( objective_generated() )  // Objective has already been generated
  return;                     // nothing to do

 // initialize Objective
 //
 // the order of the variables in the Objective Function is:
 //
 // - first f_time_horizon - init_t start-up variables
 //
 // - then, if the unit pays to shut down, f_time_horizon - init_t shut-down
 //   variables
 //
 // - then f_time_horizon active power variables (which may have the
 //   nonzero quadratic cost coefficient, while the others do not)
 //
 // - then f_time_horizon commitment variables
 //
 // - then possibly f_time_horizon primary reserve variables
 //
 // - then possibly f_time_horizon secondary reserve variables
 //
 // - then possibly the perspective-cut and reactive power variables
 //
 // - and FINALLY, if present, the single design (investment) variable
 //
 // this arrangement is exploited in add_Modification to easily map
 // indices in the coefficients of the Objective Function back into
 // indices of the original variables (and figure out the kind of variable).
 // The design variable is kept last so that the time-indexed sections above
 // start at index 0 with no leading offset

 if( v_commitment.size() != f_time_horizon )
  throw( std::logic_error( "ThermalUnitBlock::generate_objective: "
			   "v_commitment must have time horizon" ) );

 if( v_active_power.size() != f_time_horizon )
  throw( std::logic_error( "ThermalUnitBlock::generate_objective: "
			   "v_active_power must have time horizon" ) );

 if( v_start_up.size() != f_time_horizon - init_t )
  throw( std::logic_error( "ThermalUnitBlock::generate_objective: "
			   "v_start_up must have size horizon - init_t" ) );

 DQuadFunction::v_coeff_triple vars;

 // NOTE: the design (investment) variable, if present, is added LAST (after
 // the reactive power variables); see the comment there. This keeps the
 // index mapping done by handle_objective_change() for all the other
 // (time-indexed) variable sections free of any leading offset.

 // add the start-up variables- - - - - - - - - - - - - - - - - - - - - - - -
 // add start-up variables for tbin and T formulations
 //if( ( AR & FormMsk ) == tbinForm || ( AR & FormMsk ) == TForm )
 for( Index t = init_t ; t < f_time_horizon ; ++t )
  vars.push_back( std::make_tuple( & v_start_up[ t - init_t ] ,
                                   f_scale * v_StartUpCost[ t ] , 0 ) );

 // add the shut-down variables, if the unit pays anything to shut down - - -
 f_shut_down_in_obj = ! v_ShutDownCost.empty();
 if( f_shut_down_in_obj )
  for( Index t = init_t ; t < f_time_horizon ; ++t )
   vars.push_back( std::make_tuple( & v_shut_down[ t - init_t ] ,
                                    f_scale * v_ShutDownCost[ t ] , 0 ) );
 // add start-up variables for pt, DP, SU, SD and SUSD formulations
 /*
 if( ( AR & FormMsk ) == ptForm || ( AR & FormMsk ) == DPForm ||
     ( AR & FormMsk ) == SUForm || ( AR & FormMsk ) == SDForm ||
     ( AR & FormMsk ) == SUSDForm )
   for( Index j = 0 ; j < v_Y_minus.size() ; ++j )
   if( v_Y_minus[ j ].second != f_time_horizon + 1 )
    vars.push_back( std::make_tuple( &v_commitment_minus[ j ] ,
                                     f_scale *
                                     v_StartUpCost[ v_Y_minus[ j ].second - 1 ] ,
                                     0.0 ) );
 */

 // add the active power variables- - - - - - - - - - - - - - - - - - - - - -
 for( Index t = 0 ; t < f_time_horizon ; ++t )
  vars.push_back( std::make_tuple( & v_active_power[ t ] ,
                                   f_scale * v_LinearTerm[ t ] ,
                                   AR & PCuts ? 0 : f_scale * v_QuadTerm[ t ] ) );

 // add the commitment variables- - - - - - - - - - - - - - - - - - - - - - -
 for( Index t = 0 ; t < f_time_horizon ; ++t )
  vars.push_back( std::make_tuple( &v_commitment[ t ] ,
                                   f_scale * v_ConstTerm[ t ] , 0 ) );

 // the schedule is that of one unit, from which each of the copies deviates
 // on its own, so that the deviation costs the scale factor times what one
 // copy pays [see UnitBlock::scale()]
 if ( ! v_RefSchedule.empty() )
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   vars.push_back( std::make_tuple( & v_abs_ref_schedule[ t ] , f_scale ,
                                    0 ) );

 // the reserve terms: a reserve enters the Objective if its Variable exist
 // and either its cost is nonzero at some instant or the Configuration of
 // the Objective asks for it (bit 0 the primary, bit 1 the secondary), \p
 // objc or else the one of the BlockConfig; an absent cost is 0
 int objr = 0;
 auto extract_objr = [ & objr ]( Configuration * c ) {
  if( auto sc = dynamic_cast< SimpleConfiguration< int > * >( c ) ) {
   objr = sc->f_value;
   return( true );
   }
  return( false );
  };
 if( ( ! extract_objr( objc ) ) && f_BlockConfig )
  extract_objr( f_BlockConfig->f_objective_Configuration );

 auto nonzero = []( const std::vector< double > & cost ) {
  return( std::any_of( cost.begin() , cost.end() ,
                       []( double c ) { return( c != 0 ); } ) );
  };

 f_primary_in_obj = ( reserve_vars & 1u ) &&
                    ( ! v_primary_spinning_reserve.empty() ) &&
                    ( ( objr & 1 ) ||
                      nonzero( v_PrimarySpinningReserveCost ) );

 f_secondary_in_obj = ( reserve_vars & 2u ) &&
                      ( ! v_secondary_spinning_reserve.empty() ) &&
                      ( ( objr & 2 ) ||
                        nonzero( v_SecondarySpinningReserveCost ) );

 if( f_primary_in_obj ) {
  // add the primary spinning reserve variables - - - - - - - - - - - - - - -
  if( v_primary_spinning_reserve.size() != f_time_horizon )
   throw( std::logic_error( "ThermalUnitBlock::generate_objective: v_primary_"
                            "spinning_reserve must have size equal to the "
                            "time horizon." ) );

  for( Index t = 0 ; t < f_time_horizon ; ++t )
   vars.push_back( std::make_tuple( &v_primary_spinning_reserve[ t ] ,
                                    v_PrimarySpinningReserveCost.empty() ?
                                    0.0 : f_scale *
                                    v_PrimarySpinningReserveCost[ t ] ,
                                    0.0 ) );
  }

 if( f_secondary_in_obj ) {
  // add the secondary spinning reserve variables - - - - - - - - - - - - - -
  if( v_secondary_spinning_reserve.size() != f_time_horizon )
   throw( std::logic_error( "ThermalUnitBlock::generate_objective: v_secondary"
                            "_spinning_reserve must have size equal to the "
                            "time horizon." ) );

  for( Index t = 0 ; t < f_time_horizon ; ++t )
   vars.push_back( std::make_tuple( &v_secondary_spinning_reserve[ t ] ,
                                    v_SecondarySpinningReserveCost.empty() ?
                                    0.0 : f_scale *
                                    v_SecondarySpinningReserveCost[ t ] ,
                                    0.0 ) );
  }

 if( AR & PCuts ) {
  // add the perspective cuts variables - - - - - - - - - - - - - - - - - - -
  // tbin, T and pt formulations- - - - - - - - - - - - - - - - - - - - - - -
  if( ( AR & FormMsk ) == tbinForm || ( AR & FormMsk ) == TForm ||
   ( AR & FormMsk ) == ptForm )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    vars.push_back( std::make_tuple( &v_cut[ t ] ,
                                     f_scale * v_QuadTerm[ t ] , 0.0 ) );
  // DP formulation - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if( ( AR & FormMsk ) == DPForm )
   for( Index i = 0 ; i < v_Z_h_k.size() ; ++i )
    vars.push_back( std::make_tuple( &v_cut_h_k[ i ] ,
                                     f_scale *
                                     v_QuadTerm[ v_Z_h_k[ i ].first ] , 0.0 ) );
  // SU formulation - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if( ( AR & FormMsk ) == SUForm )
   for( Index i = 0 ; i < v_Z_h.size() ; ++i )
    vars.push_back( std::make_tuple( &v_cut_h[ i ] ,
                                     f_scale * v_QuadTerm[ v_Z_h[ i ].first ] ,
                                     0.0 ) );
  // SD formulation - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if( ( AR & FormMsk ) == SDForm )
   for( Index i = 0 ; i < v_Z_k.size() ; ++i )
    vars.push_back( std::make_tuple( &v_cut_k[ i ] ,
                                     f_scale * v_QuadTerm[ v_Z_k[ i ].first ] ,
                                     0.0 ) );
  // SUSD formulation - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if( ( AR & FormMsk ) == SUSDForm )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    vars.push_back( std::make_tuple( &v_cut_teta[ t ] ,
                                     f_scale * v_QuadTerm[ t ] , 0.0 ) );
  }

 // add the reactive power variables LAST - - - - - - - - - - - - - - - - - -
 // q[t] carries no cost in the unit objective, so its coefficient is 0 unless
 // a dualizing Solver sets v_ReactiveLinearTerm (see set_reactive_linear_term):
 // keeping them last lets handle_objective_change() recognise a reactive
 // coefficient change as the trailing [ num_active_var - th , num_active_var )
 // index block, independent of which optional sections precede it.
 if( f_reactive_power )
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   vars.push_back( std::make_tuple(
    & v_reactive_power[ t ] ,
    f_scale * ( v_ReactiveLinearTerm.empty() ? 0.0 : v_ReactiveLinearTerm[ t ] ) ,
    0.0 ) );

 // add the design (investment) variable LAST - - - - - - - - - - - - - - - -
 // x carries the (per-module, scaled) investment cost; a dualizing Solver may
 // change this coefficient (e.g. the non-anticipativity multiplier in a nested
 // Lagrangian). Keeping it as the trailing, single-variable block lets
 // handle_objective_change() recognise a design-coefficient change as the last
 // index num_active_var - 1, peel it off, and route it to
 // update_objective_investment() (which issues a eSetInvCost
 // ThermalUnitBlockMod for the DP Solvers), independent of all the optional
 // sections that may precede it.
 if( f_InvestmentCost != 0 ) {
  vars.push_back( std::make_tuple( & design ,
                                   f_scale * f_InvestmentCost , 0 ) );
  // seed the change-detection cache with the initial design cost, so that the
  // first (no-op) rewrite of the coefficient vector by a dualizing Solver does
  // not issue a useless eSetInvCost (the DP Solvers read the initial value
  // directly at setup via get_design_cost()); unscaled, to match get_design_cost()
  f_last_design_cost = f_InvestmentCost;
  }

 objective.set_function( new DQuadFunction( std::move( vars ) ) );
 objective.set_sense( Objective::eMin );

 // set Block objective
 set_objective( &objective , eNoMod );

 set_objective_generated();

}  // end( ThermalUnitBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*---------------- METHODS FOR CHECKING THE ThermalUnitBlock ---------------*/
/*--------------------------------------------------------------------------*/

// the tolerance and the kind of violation that a Configuration carries: the
// one it is given, or the one of the BlockConfig if that carries none

static void extract_tolerance( Configuration * fsbc , BlockConfig * bcfg ,
			       double & tol , bool & rel_viol )
{
 auto extract_parameters = [ & tol , & rel_viol ]( Configuration * c )
  -> bool {
  if( auto tc = dynamic_cast< SimpleConfiguration< double > * >( c ) ) {
   tol = tc->f_value;
   return( true );
  }
  if( auto tc = dynamic_cast< SimpleConfiguration< std::pair< double , int > > * >( c ) ) {
   tol = tc->f_value.first;
   rel_viol = tc->f_value.second;
   return( true );
  }
  return( false );
 };

 if( ( ! extract_parameters( fsbc ) ) && bcfg )
  extract_parameters( bcfg->f_is_feasible_Configuration );
 }

/*--------------------------------------------------------------------------*/

bool ThermalUnitBlock::is_feasible( bool useabstract , Configuration * fsbc )
{
 double tol = DefaultFeasTol;
 bool rel_viol = true;
 extract_tolerance( fsbc , f_BlockConfig , tol , rel_viol );

 return(
  UnitBlock::is_feasible( useabstract )
  // Variables
  && ColVariable::is_feasible( v_start_up , tol )
  && ColVariable::is_feasible( v_shut_down , tol )
  && ColVariable::is_feasible( v_primary_spinning_reserve , tol )
  && ColVariable::is_feasible( v_secondary_spinning_reserve , tol )
  && ColVariable::is_feasible( v_commitment , tol )
  && ColVariable::is_feasible( v_commitment_plus , tol )
  && ColVariable::is_feasible( v_commitment_minus , tol )
  && ColVariable::is_feasible( v_active_power , tol )
  && ColVariable::is_feasible( v_reactive_power , tol )
  && ColVariable::is_feasible( v_active_power_h_k , tol )
  && ColVariable::is_feasible( v_active_power_h , tol )
  && ColVariable::is_feasible( v_active_power_k , tol )
  && ColVariable::is_feasible( v_cut , tol )
  && ColVariable::is_feasible( v_cut_h_k , tol )
  && ColVariable::is_feasible( v_cut_h , tol )
  && ColVariable::is_feasible( v_cut_k , tol )
  && ColVariable::is_feasible( v_cut_teta , tol )
  && ColVariable::is_feasible( v_abs_ref_schedule , tol )
  && design.is_feasible( tol )
  // Constraints: notice that the ZOConstraints are not checked, since the
  // corresponding check is made on the ColVariable
  && RowConstraint::is_feasible( CommitmentDesign_Const , tol , rel_viol )
  && RowConstraint::is_feasible( StartUp_ShutDown_Variables_Const , tol , rel_viol )
  && RowConstraint::is_feasible( StartUp_Const , tol , rel_viol )
  && RowConstraint::is_feasible( ShutDown_Const , tol , rel_viol )
  && RowConstraint::is_feasible( RampUp_Const , tol , rel_viol )
  && RowConstraint::is_feasible( RampDown_Const , tol , rel_viol )
  && RowConstraint::is_feasible( MinPower_Const , tol , rel_viol )
  && RowConstraint::is_feasible( MaxPower_Const , tol , rel_viol )
  && RowConstraint::is_feasible( PrimaryRho_Const , tol , rel_viol )
  && RowConstraint::is_feasible( SecondaryRho_Const , tol , rel_viol )
  && RowConstraint::is_feasible( Reserve_Const , tol , rel_viol )
  && RowConstraint::is_feasible( Eq_ActivePower_Const , tol , rel_viol )
  && RowConstraint::is_feasible( Eq_Commitment_Const , tol , rel_viol )
  && RowConstraint::is_feasible( Eq_StartUp_Const , tol , rel_viol )
  && RowConstraint::is_feasible( Eq_ShutDown_Const , tol , rel_viol )
  && RowConstraint::is_feasible( Network_Const , tol , rel_viol )
  && RowConstraint::is_feasible( Init_PC_Const , tol , rel_viol )
  && RowConstraint::is_feasible( Eq_PC_Const , tol , rel_viol )
  && RowConstraint::is_feasible( Max_SUSD_PC_Const , tol , rel_viol )
  && RowConstraint::is_feasible( fixed_to_max_Power_Const , tol , rel_viol )
  && RowConstraint::is_feasible( PC_cuts , tol , rel_viol )
  && RowConstraint::is_feasible( ShutDownZero_Const , tol , rel_viol )
  && RowConstraint::is_feasible( MaxPower5_Const.rows , tol , rel_viol )
  && RowConstraint::is_feasible( MaxPower6_Const.rows , tol , rel_viol )
  && RowConstraint::is_feasible( RampUpSUSD_Const.rows , tol , rel_viol )
  && RowConstraint::is_feasible( RampDownSUSD_Const.rows , tol , rel_viol )
  && RowConstraint::is_feasible( Commitment_fixed_to_One_Const , tol ,
				 rel_viol )
  && RowConstraint::is_feasible( Reference_Schedule_Const , tol , rel_viol )
  && RowConstraint::is_feasible( ReactivePower_Bound_Const , tol
				 , rel_viol )
  && RowConstraint::is_feasible( ReactivePowerMax_Const , tol , rel_viol )
  && RowConstraint::is_feasible( ReactivePowerMin_Const , tol , rel_viol )
  && RowConstraint::is_feasible( ActivePower_Bound_Const , tol , rel_viol ) );

}  // end( ThermalUnitBlock::is_feasible )

/*--------------------------------------------------------------------------*/

bool ThermalUnitBlock::is_sol_feasible_physical( void ) const
{
 /* The check of is_sol_feasible() reads the schedule of the unit and tests it
  * against the data of the unit, the Variable of the schedule that are fixed
  * included: a unit that carries something the schedule does not answer for,
  * i.e. a dimensioning variable, the reactive power, a reference schedule, a
  * scale of its own or a fixed Variable of a formulation that the schedule
  * only implies, is left to the check of the base class, which goes through
  * the Variable. */
 if( ( f_scale != 1 ) || ( f_InvestmentCost != 0 ) || ( f_Capacity != 0 ) ||
     ( ! v_RefSchedule.empty() ) ||
     ( ! v_MinReactivePower.empty() ) || ( ! v_MaxReactivePower.empty() ) ||
     ( ! v_MinReactivePowerOn.empty() ) || ( ! v_MaxReactivePowerOn.empty() ) )
  return( false );

 const auto fixed = []( const std::vector< ColVariable > & v ) {
  return( std::any_of( v.begin() , v.end() ,
		       []( const ColVariable & x ) { return( x.is_fixed() ); } ) );
  };
 for( const auto * v : { & v_commitment_plus , & v_commitment_minus ,
			 & v_reactive_power , & v_active_power_h_k ,
			 & v_active_power_h , & v_active_power_k , & v_cut ,
			 & v_cut_h_k , & v_cut_h , & v_cut_k , & v_cut_teta ,
			 & v_abs_ref_schedule } )
  if( fixed( *v ) )
   return( false );

 return( true );

 }  // end( ThermalUnitBlock::is_sol_feasible_physical )

/*--------------------------------------------------------------------------*/

bool ThermalUnitBlock::is_sol_feasible( Solution * sol , Configuration * fsbc )
{
 if( ! is_sol_feasible_physical() )
  // the schedule does not answer for this unit
  return( Block::is_sol_feasible( sol , fsbc ) );

 auto usol = dynamic_cast< UnitBlockSolution * >( sol );
 if( ! usol )
  throw( std::invalid_argument( "ThermalUnitBlock::is_sol_feasible: the "
				"Solution is not a UnitBlockSolution" ) );

 // the feasible region of a unit is bounded, hence it has no rays
 if( usol->is_direction() )
  return( false );

 const auto & PP = usol->get_active_power();
 const auto & UU = usol->get_commitment();
 if( ( PP.shape()[ 0 ] < 1 ) || ( PP.shape()[ 1 ] < f_time_horizon ) ||
     ( UU.shape()[ 0 ] < 1 ) || ( UU.shape()[ 1 ] < f_time_horizon ) )
  return( false );  // it holds no schedule of this unit

 /* The reserves are part of the schedule when the unit produces them: a
  * Solution that does not carry those the unit has cannot answer for the
  * constraints they are in. */
 const auto & R1 = usol->get_primary_spinning_reserve();
 const auto & R2 = usol->get_secondary_spinning_reserve();
 const bool has_r1 = ( R1.shape()[ 0 ] >= 1 ) &&
                     ( R1.shape()[ 1 ] >= f_time_horizon );
 const bool has_r2 = ( R2.shape()[ 0 ] >= 1 ) &&
                     ( R2.shape()[ 1 ] >= f_time_horizon );
 if( ( ( ! v_primary_spinning_reserve.empty() ) && ( ! has_r1 ) ) ||
     ( ( ! v_secondary_spinning_reserve.empty() ) && ( ! has_r2 ) ) )
  return( false );

 double tol = DefaultFeasTol;
 bool rel_viol = true;
 extract_tolerance( fsbc , f_BlockConfig , tol , rel_viol );

 /* lhs <= rhs, up to the tolerance and, in any case, up to a few ulp of the
  * numbers at hand: the data of the unit and the values of the schedule are
  * floating point, hence a constraint that is tight comes out violated by the
  * rounding of one difference, and the check that goes through the abstract
  * representation sees it violated or not depending on the order in which the
  * Function of the row sums its terms. */
 auto le = [ tol , rel_viol ]( double lhs , double rhs ) {
  const double big = std::max( { 1.0 , std::abs( lhs ) , std::abs( rhs ) } );
  return( lhs - rhs <= std::max( rel_viol ? tol * big : tol , 1e-12 * big ) );
  };

 /* The constraints below are those of the unit for an integral commitment,
  * which is what every formulation of it encodes; a commitment that is not
  * integral is one the formulations do not agree on, and what cannot be told
  * is not declared feasible. */
 std::vector< bool > on( f_time_horizon );
 for( Index t = 0 ; t < f_time_horizon ; ++t ) {
  const double u = UU[ 0 ][ t ];
  if( ( u < - tol ) || ( u > 1 + tol ) )
   return( false );
  if( ( u > tol ) && ( u < 1 - tol ) )
   return( false );
  on[ t ] = u > 0.5;
  }

 /* A Variable that is fixed only holds the value it is fixed to, which is
  * how the state the unit comes from is written when it leaves it no choice
  * for the first time instants; a schedule that says otherwise is none of
  * this unit. The Variable are only read, never written. */
 for( Index t = 0 ; t < f_time_horizon ; ++t ) {
  if( ( t < v_commitment.size() ) && v_commitment[ t ].is_fixed() &&
      ( ! le( std::abs( UU[ 0 ][ t ] - v_commitment[ t ].get_value() ) ,
	      0 ) ) )
   return( false );
  if( ( t < v_active_power.size() ) && v_active_power[ t ].is_fixed() &&
      ( ! le( std::abs( PP[ 0 ][ t ] - v_active_power[ t ].get_value() ) ,
	      0 ) ) )
   return( false );
  if( has_r1 && ( t < v_primary_spinning_reserve.size() ) &&
      v_primary_spinning_reserve[ t ].is_fixed() &&
      ( ! le( std::abs( R1[ 0 ][ t ] -
			v_primary_spinning_reserve[ t ].get_value() ) , 0 ) ) )
   return( false );
  if( has_r2 && ( t < v_secondary_spinning_reserve.size() ) &&
      v_secondary_spinning_reserve[ t ].is_fixed() &&
      ( ! le( std::abs( R2[ 0 ][ t ] -
			v_secondary_spinning_reserve[ t ].get_value() ) , 0 ) ) )
   return( false );
  }

 /* The same for the start-up and shut-down indicators, which the Solution
  * holds if it has saved them and the commitment implies otherwise [see
  * derive_start_up()]. */
 const auto fixed = []( const ColVariable & v ) { return( v.is_fixed() ); };
 if( std::any_of( v_start_up.begin() , v_start_up.end() , fixed ) ||
     std::any_of( v_shut_down.begin() , v_shut_down.end() , fixed ) ) {
  std::vector< double > su;
  std::vector< double > sd;
  const auto tsol = dynamic_cast< const ThermalUnitBlockSolution * >( sol );
  if( tsol && ( tsol->get_start_up().size() == v_start_up.size() ) &&
      ( tsol->get_shut_down().size() == v_start_up.size() ) ) {
   su = tsol->get_start_up();
   sd = tsol->get_shut_down();
   }
  else {
   std::vector< double > u( f_time_horizon );
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    u[ t ] = UU[ 0 ][ t ];
   derive_start_up( u , su , sd );
   }

  for( Index j = 0 ; j < v_start_up.size() ; ++j )
   if( v_start_up[ j ].is_fixed() &&
       ( ! le( std::abs( su[ j ] - v_start_up[ j ].get_value() ) , 0 ) ) )
    return( false );
  for( Index j = 0 ; ( j < v_shut_down.size() ) && ( j < sd.size() ) ; ++j )
   if( v_shut_down[ j ].is_fixed() &&
       ( ! le( std::abs( sd[ j ] - v_shut_down[ j ].get_value() ) , 0 ) ) )
    return( false );
  }

 // the power and the reserves against the operational bounds of the unit
 for( Index t = 0 ; t < f_time_horizon ; ++t ) {
  const double p = PP[ 0 ][ t ];
  const double r1 = has_r1 ? R1[ 0 ][ t ] : 0;
  const double r2 = has_r2 ? R2[ 0 ][ t ] : 0;

  if( ( ! le( 0 , r1 ) ) || ( ! le( 0 , r2 ) ) )
   return( false );

  if( ! on[ t ] ) {  // the unit is off: it produces nothing
   if( ( ! le( p , 0 ) ) || ( ! le( 0 , p ) ) ||
       ( ! le( r1 + r2 , 0 ) ) )
    return( false );
   continue;
   }

  // the reserves take room from the power, on both sides
  if( ( ! le( get_operational_min_power( t ) , p - r1 - r2 ) ) ||
      ( ! le( p + r1 + r2 , get_operational_max_power( t ) ) ) )
   return( false );

  // and each of them is a fraction of the power produced
  if( has_r1 && ( ! v_PrimaryRho.empty() ) &&
      ( ! le( r1 , v_PrimaryRho[ t ] * p ) ) )
   return( false );
  if( has_r2 && ( ! v_SecondaryRho.empty() ) &&
      ( ! le( r2 , v_SecondaryRho[ t ] * p ) ) )
   return( false );
  }

 // the ramps, with the limits of the start-up and of the shut-down and the
 // state the unit comes from
 for( Index t = 0 ; t < f_time_horizon ; ++t ) {
  const double p = PP[ 0 ][ t ];
  const bool prev_on = t > 0 ? on[ t - 1 ] : ( f_InitUpDownTime > 0 );
  const double prev_p = t > 0 ? PP[ 0 ][ t - 1 ]
			      : ( f_InitUpDownTime > 0 ? f_InitialPower : 0 );

  if( on[ t ] ) {
   if( prev_on ) {
    if( ( ! v_DeltaRampUp.empty() ) &&
	( ! le( p - prev_p , v_DeltaRampUp[ t ] ) ) )
     return( false );
    if( ( ! v_DeltaRampDown.empty() ) &&
	( ! le( prev_p - p , v_DeltaRampDown[ t ] ) ) )
     return( false );
    }
   else                                   // the unit starts up at t
    if( ! le( p , v_StartUpLimit[ t ] ) )
     return( false );
   }
  else
   if( prev_on )                          // the unit shuts down at t
    if( ! le( prev_p , v_ShutDownLimit[ t ] ) )
     return( false );
  }

 // the minimum up and down times, the state the unit comes from included
 if( ( f_MinUpTime > 1 ) || ( f_MinDownTime > 1 ) ) {
  const int init = f_InitUpDownTime;

  if( init > 0 ) {  // the unit is on, and has been for init time instants
   for( Index t = 0 ; ( t < f_time_horizon ) &&
	  ( int( t ) + init < int( f_MinUpTime ) ) ; ++t )
    if( ! on[ t ] )
     return( false );
   }
  else              // it is off, and has been for -init time instants
   for( Index t = 0 ; ( t < f_time_horizon ) &&
	  ( int( t ) - init < int( f_MinDownTime ) ) ; ++t )
    if( on[ t ] )
     return( false );

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   const bool prev_on = t > 0 ? on[ t - 1 ] : ( f_InitUpDownTime > 0 );

   if( on[ t ] && ( ! prev_on ) )         // the unit starts up at t
    for( Index s = t + 1 ;
	 ( s < f_time_horizon ) && ( s < t + f_MinUpTime ) ; ++s )
     if( ! on[ s ] )
      return( false );

   if( ( ! on[ t ] ) && prev_on )         // the unit shuts down at t
    for( Index s = t + 1 ;
	 ( s < f_time_horizon ) && ( s < t + f_MinDownTime ) ; ++s )
     if( on[ s ] )
      return( false );
   }
  }

 return( true );

 }  // end( ThermalUnitBlock::is_sol_feasible )

/*--------------------------------------------------------------------------*/
/*-------- METHODS FOR LOADING, PRINTING & SAVING THE ThermalUnitBlock -----*/
/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::serialize( netCDF::NcGroup & group ) const
{
 UnitBlock::serialize( group );

 // Serialize scalar variables

 if( f_InvestmentCost != 0 )
  ::serialize( group , "InvestmentCost" , netCDF::NcDouble() ,
               f_InvestmentCost );

 if( f_Capacity != 0 )
  ::serialize( group , "Capacity" , netCDF::NcDouble() , f_Capacity );

 if( f_scale != 1 )
  ::serialize( group , "Scale" , netCDF::NcDouble() , f_scale );

 ::serialize( group , "InitialPower" , netCDF::NcDouble() , f_InitialPower );
 ::serialize( group , "MinUpTime" , netCDF::NcUint() , f_MinUpTime );
 ::serialize( group , "MinDownTime" , netCDF::NcUint() , f_MinDownTime );
 ::serialize( group , "InitUpDownTime" , netCDF::NcInt() , f_InitUpDownTime );

 // Serialize one-dimensional variables

 auto TimeHorizon = group.getDim( "TimeHorizon" );
 auto NumberIntervals = group.getDim( "NumberIntervals" );

 /* This lambda identifies the appropriate dimension for the given variable
  * (whose name is "var_name") and serializes the variable. The variable may
  * have any of the following dimensions: TimeHorizon, NumberIntervals,
  * 1. "allow_scalar_var" indicates whether the variable can be serialized as
  * a scalar variable (in which case the variable must have dimension 1). */
 auto serialize = [ &group , &TimeHorizon , &NumberIntervals ]
  ( const std::string & var_name , const std::vector< double > & data ,
    const netCDF::NcType & ncType = netCDF::NcDouble() ,
    bool allow_scalar_var = true ) {
  if( data.empty() )
   return;
  netCDF::NcDim dimension;
  if( data.size() == TimeHorizon.getSize() )
   dimension = TimeHorizon;
  else if( data.size() == NumberIntervals.getSize() )
   dimension = NumberIntervals;
  else if( data.size() != 1 )
   throw( std::logic_error( "ThermalUnitBlock::serialize: invalid dimension "
                            "for variable " + var_name + ": " +
                            std::to_string( data.size() ) +
                            ". Its dimension must be one of the following: "
                            "TimeHorizon, NumberIntervals, 1." ) );

  ::serialize( group , var_name , ncType , dimension , data ,
               allow_scalar_var );
 };

 serialize( "MinPower" , v_MinPower );
 serialize( "MaxPower" , v_MaxPower );
 serialize( "MaxReactivePower" , v_MaxReactivePower );
 serialize( "MinReactivePower" , v_MinReactivePower );
 serialize( "MaxReactivePowerOn" , v_MaxReactivePowerOn );
 serialize( "MinReactivePowerOn" , v_MinReactivePowerOn );
 serialize( "Availability" , v_Availability );
 serialize( "DeltaRampUp" , v_DeltaRampUp );
 serialize( "DeltaRampDown" , v_DeltaRampDown );
 serialize( "PrimaryRho" , v_PrimaryRho );
 serialize( "SecondaryRho" , v_SecondaryRho );
 serialize( "QuadTerm" , v_QuadTerm );
 serialize( "LinearTerm" , v_LinearTerm );
 // the spinning-reserve linear costs are, like LinearTerm, the Lagrangian dual
 // of a relaxed system constraint (here the reserve demand) pushed onto the
 // unit; serialise them too so a dumped unit carries its reserve prices (empty,
 // hence skipped, when the unit offers no reserve or none is priced)
 serialize( "PrimarySpinningReserveCost" , v_PrimarySpinningReserveCost );
 serialize( "SecondarySpinningReserveCost" , v_SecondarySpinningReserveCost );
 serialize( "ReactiveLinearTerm" , v_ReactiveLinearTerm );
 serialize( "ConstTerm" , v_ConstTerm );
 serialize( "StartUpCost" , v_StartUpCost );
 serialize( "ShutDownCost" , v_ShutDownCost );
 serialize( "FixedConsumption" , v_FixedConsumption );
 serialize( "InertiaCommitment" , v_InertiaCommitment );
 serialize( "StartUpLimit" , v_StartUpLimit );
 serialize( "ShutDownLimit" , v_ShutDownLimit );

}  // end( ThermalUnitBlock::serialize )

/*--------------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/

Solution * ThermalUnitBlock::get_Solution( Configuration * csolc ,
             bool emptys )
{
 Index wsol = 15;
 if( ( ! csolc ) && f_BlockConfig )
  csolc = f_BlockConfig->f_solution_Configuration;

 if( auto config = dynamic_cast< SimpleConfiguration< int > * >( csolc ) )
  wsol = config->f_value;

 // call the method of the base class
 auto * sol = dynamic_cast< ThermalUnitBlockSolution * >(
                     UnitBlock::get_Solution( csolc , emptys ) );
 assert( sol );

 if( ! emptys )
  sol->read( this );

 return( sol );
 }

/*--------------------------------------------------------------------------*/

UnitBlockSolution * ThermalUnitBlock::new_Solution( void ) const {
 return( new ThermalUnitBlockSolution() );
 }

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::add_Modification( sp_Mod mod , ChnlName chnl )
{
 if( mod->concerns_Block() ) {
  mod->concerns_Block( false );
  guts_of_add_Modification( mod.get() , chnl );
 }

 Block::add_Modification( mod , chnl );

}  // end( ThermalUnitBlock::add_Modification )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::check_power_limits( const std::string & who ,
                                           Index t , double max_power ,
                                           double availability ) const
{
 const auto at = " at time " + std::to_string( t );

 if( ( availability < 0 ) || ( availability > 1 ) )
  throw( std::logic_error( who + ": availability (" +
                           std::to_string( availability ) + ")" + at +
                           " is not between 0 and 1" ) );

 if( v_MinPower[ t ] > max_power )
  throw( std::logic_error( who + ": maximum power (" +
                           std::to_string( max_power ) + ")" + at +
                           " is below the minimum power (" +
                           std::to_string( v_MinPower[ t ] ) + ")" ) );

 const double opmin = compute_operational_min_power( v_MinPower[ t ] ,
                                                     availability );
 const double opmax = compute_operational_max_power( max_power ,
                                                     availability );
 if( opmin > opmax )
  throw( std::logic_error( who + ": availability (" +
                           std::to_string( availability ) + ")" + at +
                           " is not consistent" ) );

 // a given limit has to stay within the operational bounds, a default one
 // follows them [see set_default_limits()]
 auto check = [ & ]( bool dflt , double lim , const char * name ) {
  if( ( ! dflt ) && ( ( lim < opmin ) || ( lim > opmax ) ) )
   throw( std::logic_error( who + ": " + name + at + " is " +
                            std::to_string( lim ) + ", but the operational "
                            "bounds would be " + std::to_string( opmin ) +
                            " and " + std::to_string( opmax ) ) );
  };
 check( f_default_start_up_limit , v_StartUpLimit[ t ] , "start-up limit" );
 check( f_default_shut_down_limit , v_ShutDownLimit[ t ] ,
        "shut-down limit" );

 }  // end( ThermalUnitBlock::check_power_limits )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_default_limits( void )
{
 if( f_default_start_up_limit )
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   v_StartUpLimit[ t ] = get_operational_min_power( t );

 if( f_default_shut_down_limit )
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   v_ShutDownLimit[ t ] = get_operational_min_power( t );

 }  // end( ThermalUnitBlock::set_default_limits )

/*--------------------------------------------------------------------------*/

bool ThermalUnitBlock::guts_of_set_power_limits( const std::string & who ,
                                                 bool availability ,
                                                 const Subset & ts ,
                                                 MF_dbl_it values ,
                                                 ModParam issuePMod ,
                                                 ModParam issueAMod )
{
 // if nothing changes, return
 auto val = values;
 bool identical = true;
 for( auto t : ts )
  if( ( availability ? get_availability( t ) : v_MaxPower[ t ] ) !=
      *( val++ ) ) {
   identical = false;
   break;
   }
 if( identical )
  return( false );

 // the new data are checked before anything changes; an instant that
 // appears more than once takes its last value, as the assignment does
 std::map< Index , double > last;
 val = values;
 for( auto t : ts )
  last[ t ] = *( val++ );
 for( const auto & tv : last )
  check_power_limits( who , tv.first ,
                      availability ? v_MaxPower[ tv.first ] : tv.second ,
                      availability ? tv.second :
                                     get_availability( tv.first ) );

 if( ! not_dry_run( issuePMod ) )
  return( true );

 // change the physical representation: the data, the default limits and
 // the numbers of ramp steps, which depend on the operational bounds
 const auto old_max_power = v_MaxPower;
 const auto old_availability = v_Availability;
 const auto old_start_up_limit = v_StartUpLimit;
 const auto old_shut_down_limit = v_ShutDownLimit;

 auto & v = availability ? v_Availability : v_MaxPower;
 if( v.empty() )
  v.assign( f_time_horizon , availability ? 1.0 : 0.0 );
 val = values;
 for( auto t : ts )
  v[ t ] = *( val++ );
 set_default_limits();
 compute_ramp_steps();

 // change the abstract representation: every row of every formulation is
 // written anew, and those that differ change [see update_rows()]; if it
 // cannot be, the physical representation is restored, so that the two do
 // not disagree
 if( not_dry_run( issueAMod ) && constraints_generated() )
  try {
   update_rows( un_ModBlock( issueAMod ) );
   }
  catch( ... ) {
   v_MaxPower = old_max_power;
   v_Availability = old_availability;
   v_StartUpLimit = old_start_up_limit;
   v_ShutDownLimit = old_shut_down_limit;
   compute_ramp_steps();
   throw;
   }

 return( true );

 }  // end( ThermalUnitBlock::guts_of_set_power_limits )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_availability( MF_dbl_it values ,
                                         Subset && subset ,
                                         const bool ordered ,
                                         ModParam issuePMod ,
                                         ModParam issueAMod )
{
 if( subset.empty() )
  return;

 std::vector< double > sorted;  // the values in the order of the subset
 if( ! ordered )
  sort_by_index( subset , values , sorted );

 if( *std::max_element( subset.begin() , subset.end() ) >= f_time_horizon )
  throw( std::invalid_argument(
   "ThermalUnitBlock::set_availability: invalid index in subset." ) );

 if( ! guts_of_set_power_limits( "ThermalUnitBlock::set_availability" ,
                                 true , subset , values , issuePMod ,
                                 issueAMod ) )
  return;  // nothing changes

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockSbstMod >(
                            this , ThermalUnitBlockMod::eSetAv ,
                            std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_availability( subset ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_availability( MF_dbl_it values ,
                                         Range rng ,
                                         ModParam issuePMod ,
                                         ModParam issueAMod )
{
 rng.second = std::min( rng.second , f_time_horizon );
 if( rng.second <= rng.first )
  return;

 Subset ts( rng.second - rng.first );
 std::iota( ts.begin() , ts.end() , rng.first );
 if( ! guts_of_set_power_limits( "ThermalUnitBlock::set_availability" ,
                                 true , ts , values , issuePMod ,
                                 issueAMod ) )
  return;  // nothing changes

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockRngdMod >(
                            this , ThermalUnitBlockMod::eSetAv , rng ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_availability( range ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_maximum_power( MF_dbl_it values ,
                                          Subset && subset ,
                                          const bool ordered ,
                                          ModParam issuePMod ,
                                          ModParam issueAMod )
{
 if( subset.empty() )
  return;

 std::vector< double > sorted;  // the values in the order of the subset
 if( ! ordered )
  sort_by_index( subset , values , sorted );

 if( *std::max_element( subset.begin() , subset.end() ) >= f_time_horizon )
  throw( std::invalid_argument(
   "ThermalUnitBlock::set_maximum_power: invalid index in subset." ) );

 if( ! guts_of_set_power_limits( "ThermalUnitBlock::set_maximum_power" ,
                                 false , subset , values , issuePMod ,
                                 issueAMod ) )
  return;  // nothing changes

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockSbstMod >(
                            this , ThermalUnitBlockMod::eSetMaxP ,
                            std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_maximum_power( subset ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_maximum_power( MF_dbl_it values ,
                                          Range rng ,
                                          ModParam issuePMod ,
                                          ModParam issueAMod )
{
 rng.second = std::min( rng.second , f_time_horizon );
 if( rng.second <= rng.first )
  return;

 Subset ts( rng.second - rng.first );
 std::iota( ts.begin() , ts.end() , rng.first );
 if( ! guts_of_set_power_limits( "ThermalUnitBlock::set_maximum_power" ,
                                 false , ts , values , issuePMod ,
                                 issueAMod ) )
  return;  // nothing changes

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockRngdMod >(
                            this , ThermalUnitBlockMod::eSetMaxP , rng ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_maximum_power( range ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::update_initial_power_in_cnstrs( ModParam issueAMod )
{
 // InitialPower is p_{-1}, which no Variable depends on (the interval
 // ( 0 , 0 ) of the shut-down at 0 is there whenever the unit may shut down
 // at 0, see generate_abstract_variables()): the rows that contain it are
 // written anew [see update_rows()], those of a derived class included
 update_rows( issueAMod );

}  // end( ThermalUnitBlock::update_initial_power_in_cnstrs )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::compute_ramp_steps( void )
{
 // the number of steps within which the ramps do not cover the distance
 // between the bounds, from InitialPower ( [ 0 ], see
 // compute_initial_ramp_steps() ) or from the power at t - 1 ( [ t ] ):
 // upwards, from the minimum power at t - 1 to the largest maximum power
 // at t - 1 or later, downwards from the maximum power at t - 1 to the
 // smallest minimum power at t - 1 or later, so that a row beyond the count
 // is implied by the bounds also when they vary in time (an outage at
 // t - 1 does not remove the rows that start there); the ramp of the step
 // to t is DeltaRampUp[ t ] ( DeltaRampDown[ t ] )
 v_MaxRampSteps.assign( f_time_horizon + 1 , 0 );
 v_MaxRampDownSteps.assign( f_time_horizon + 1 , 0 );
 f_derived_ramp_up_steps = f_derived_ramp_down_steps = true;

 for( Index t = 1 ; t <= f_time_horizon ; ++t ) {
  double maxp = 0;
  double minp = get_operational_min_power( t - 1 );
  for( Index s = t - 1 ; s < f_time_horizon ; ++s ) {
   maxp = std::max( maxp , get_operational_max_power( s ) );
   minp = std::min( minp , get_operational_min_power( s ) );
   }

  const double range_up = maxp - get_operational_min_power( t - 1 );
  int steps = 0;
  double reach = 0;
  while( ( steps < static_cast< int >( f_time_horizon - t ) ) &&
         ( ( reach += get_delta_ramp_up( t + steps ) ) <= range_up ) )
   ++steps;
  v_MaxRampSteps[ t ] = steps;

  const double range_dn = get_operational_max_power( t - 1 ) - minp;
  steps = 0;
  reach = 0;
  while( ( steps < static_cast< int >( f_time_horizon - t ) ) &&
         ( ( reach += get_delta_ramp_down( t + steps ) ) <= range_dn ) )
   ++steps;
  v_MaxRampDownSteps[ t ] = steps;
  }

 compute_initial_ramp_steps();

}  // end( ThermalUnitBlock::compute_ramp_steps )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::compute_initial_ramp_steps( void )
{
 // if the unit is on before the horizon, the number of the first instants
 // j that contain every j at which the ramps from InitialPower do not reach
 // the maximum (minimum) power at j, i.e., one more than the last such j (0
 // if none), so that the rows of the instants after it are implied by the
 // bounds also when they vary in time; only where the data do not give them
 if( f_derived_ramp_up_steps && ( ! v_MaxRampSteps.empty() ) ) {
  int steps = -1;
  if( f_InitUpDownTime > 0 ) {
   steps = 0;
   double reach = f_InitialPower;
   for( Index j = 0 ; j < f_time_horizon ; ++j )
    if( ( reach += get_delta_ramp_up( j ) ) <=
        get_operational_max_power( j ) )
     steps = j + 1;
   }
  v_MaxRampSteps[ 0 ] = steps;
  }

 if( f_derived_ramp_down_steps && ( ! v_MaxRampDownSteps.empty() ) ) {
  int steps = -1;
  if( f_InitUpDownTime > 0 ) {
   steps = 0;
   double reach = f_InitialPower;
   for( Index j = 0 ; j < f_time_horizon ; ++j )
    if( ( reach -= get_delta_ramp_down( j ) ) >=
        get_operational_min_power( j ) )
     steps = j + 1;
   }
  v_MaxRampDownSteps[ 0 ] = steps;
  }

}  // end( ThermalUnitBlock::compute_initial_ramp_steps )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::check_initial_power( double value ) const
{
 // the same refusals as check_data_consistency(), made before any change:
 // a unit on before the horizon below its minimum power would be off at 0
 // for some formulations and solvers and not for others
 if( value < 0 )
  throw( std::logic_error( "ThermalUnitBlock::set_initial_power: the "
                           "initial power " + std::to_string( value ) +
                           " is negative" ) );

 if( ( f_InitUpDownTime > 0 ) && ( value < v_MinPower.front() ) )
  throw( std::logic_error( "ThermalUnitBlock::set_initial_power: the unit "
                           "is on before the horizon, and the initial "
                           "power " + std::to_string( value ) + " is below "
                           "the minimum power " +
                           std::to_string( v_MinPower.front() ) ) );

 }  // end( ThermalUnitBlock::check_initial_power )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_initial_power( MF_dbl_it values ,
                                          Subset && subset ,
                                          const bool ordered ,
                                          ModParam issuePMod ,
                                          ModParam issueAMod )
{
 if( subset.empty() )
  return;

 // Find the last index 0
 auto index_it = std::find( subset.rbegin() , subset.rend() , 0 );

 if( index_it == subset.rend() )
  return;  // 0 is not in subset; return

 std::advance( values , std::distance( index_it , subset.rend() ) - 1 );

 if( f_InitialPower == *values )
  return;  // nothing changes; return

 check_initial_power( *values );

 const double old_power = f_InitialPower;
 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation
  f_InitialPower = *values;
  compute_initial_ramp_steps();
  }

 if( not_dry_run( issueAMod ) && variables_generated() )
  // Change the abstract representation; if it cannot be, the physical one
  // is restored, so that the two do not disagree
  try {
   update_initial_power_in_cnstrs( un_ModBlock( issueAMod ) );
   }
  catch( ... ) {
   f_InitialPower = old_power;
   compute_initial_ramp_steps();
   throw;
   }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockMod >(
                            this , ThermalUnitBlockMod::eSetInitP ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_initial_power( subset ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_initial_power( MF_dbl_it values ,
                                          Range rng ,
                                          ModParam issuePMod ,
                                          ModParam issueAMod )
{
 rng.second = std::min( rng.second , static_cast< Index >( 1 ) );
 if( ! ( ( rng.first <= 0 ) && ( 0 < rng.second ) ) )
  return;  // 0 does not belong to the range; return

 std::advance( values , -rng.first );

 if( f_InitialPower == *values )
  return;  // nothing changes; return

 check_initial_power( *values );

 const double old_power = f_InitialPower;
 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation
  f_InitialPower = *values;
  compute_initial_ramp_steps();
  }

 if( not_dry_run( issueAMod ) && variables_generated() )
  // Change the abstract representation; if it cannot be, the physical one
  // is restored, so that the two do not disagree
  try {
   update_initial_power_in_cnstrs( un_ModBlock( issueAMod ) );
   }
  catch( ... ) {
   f_InitialPower = old_power;
   compute_initial_ramp_steps();
   throw;
   }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockMod >(
                            this , ThermalUnitBlockMod::eSetInitP ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_initial_power( range ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_startup_costs( MF_dbl_it values ,
                                          Subset && subset ,
                                          const bool ordered ,
                                          ModParam issuePMod ,
                                          ModParam issueAMod )
{
 if( subset.empty() )
  return;

 if( v_StartUpCost.empty() ) {
  if( std::all_of( values ,
                   values + subset.size() ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  v_StartUpCost.assign( f_time_horizon , 0 );
 }

 std::vector< double > sorted;  // the values in the order of the subset
 if( ! ordered )
  sort_by_index( subset , values , sorted );

 if( subset.back() >= v_StartUpCost.size() )
  throw( std::invalid_argument(
   "ThermalUnitBlock::set_startup_costs: invalid index in subset." ) );

 if( identical( v_StartUpCost , subset , values ) )  // if nothing changes
  return;                                                 // return

 // the order of the variables in the Objective Function is:
 //
 // - first f_time_horizon - init_t start-up variables
 //
 // - then the rest ...
 //
 // this means that the start_up variable t is in position t - init_t
 // hence, those in the range [ 0 , init_t ) do not exist and their cost
 // cannot be changed
 if( subset.front() < init_t )
  throw( std::invalid_argument(
   "ThermalUnitBlock::set_startup_costs: invalid starting index in subset." ) );

 if( not_dry_run( issuePMod ) )
  // Change the physical representation
  assign( v_StartUpCost , subset , values );

 if( not_dry_run( issueAMod ) && objective_generated() ) {
  // Change the abstract representation
  Subset tmps = subset_sbtrct( subset , init_t );
  DQuadFunction::Vec_FunctionValue tmpv( values , values + subset.size() );
  for( auto & c : tmpv )  // the Objective carries f_scale times the cost
   c *= f_scale;
  static_cast< DQuadFunction * >( objective.get_function()
  )->modify_linear_coefficients( std::move( tmpv ) , std::move( tmps ) ,
                                 true , un_ModBlock( issueAMod ) );
 }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockSbstMod >(
                            this , ThermalUnitBlockMod::eSetSUC ,
                            std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_startup_costs( subset ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_startup_costs( MF_dbl_it values ,
                                          Range rng ,
                                          ModParam issuePMod ,
                                          ModParam issueAMod )
{
 rng.second = std::min( rng.second , f_time_horizon );
 if( rng.second <= rng.first )
  return;

 c_Index sz = rng.second - rng.first;
 if( v_StartUpCost.empty() ) {
  if( std::all_of( values ,
                   values + sz ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  v_StartUpCost.assign( f_time_horizon , 0 );
 }

 // If nothing changes, return
 if( std::equal( values , values + sz , v_StartUpCost.begin() + rng.first ) )
  return;

 // the order of the variables in the Objective Function is:
 //
 // - first f_time_horizon - init_t start-up variables
 //
 // - then the rest ...
 //
 // this means that the start_up variable t is in position t - init_t
 // hence, those in the range [ 0 , init_t ) do not exist and their cost
 // cannot be changed
 if( rng.first < init_t )
  throw( std::invalid_argument( "ThermalUnitBlock::set_startup_costs: invalid"
                                " starting index in range." ) );

 if( not_dry_run( issuePMod ) )
  // Change the physical representation
  std::copy( values ,
             values + sz ,
             v_StartUpCost.begin() + rng.first );

 if( not_dry_run( issueAMod ) && objective_generated() ) {
  // Change the abstract representation
  DQuadFunction::Vec_FunctionValue tmpv( values , values + sz );
  for( auto & c : tmpv )  // the Objective carries f_scale times the cost
   c *= f_scale;
  static_cast< DQuadFunction * >( objective.get_function()
  )->modify_linear_coefficients( std::move( tmpv ) ,
                                 Range( rng.first - init_t ,
                                        rng.second - init_t ) ,
                                 un_ModBlock( issueAMod ) );
 }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockRngdMod >(
                            this , ThermalUnitBlockMod::eSetSUC , rng ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_startup_costs( range ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_shutdown_costs( MF_dbl_it values ,
                                           Subset && subset ,
                                           const bool ordered ,
                                           ModParam issuePMod ,
                                           ModParam issueAMod )
{
 if( subset.empty() )
  return;

 // a unit that was built with no shut-down cost has no shut-down variable
 // in its Objective [see generate_objective()], hence its costs can only
 // stay zero (checked before anything is changed)
 if( objective_generated() && ( ! f_shut_down_in_obj ) &&
     std::any_of( values , values + subset.size() ,
                  []( double cst ) { return( cst != 0 ); } ) )
  throw( std::logic_error( "ThermalUnitBlock::set_shutdown_costs: the "
                           "Objective has no shut-down variables, the unit "
                           "having been generated with no shut-down cost" ) );

 if( v_ShutDownCost.empty() ) {
  if( std::all_of( values ,
                   values + subset.size() ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  v_ShutDownCost.assign( f_time_horizon , 0 );
  }

 std::vector< double > sorted;  // the values in the order of the subset
 if( ! ordered )
  sort_by_index( subset , values , sorted );

 if( subset.back() >= v_ShutDownCost.size() )
  throw( std::invalid_argument(
   "ThermalUnitBlock::set_shutdown_costs: invalid index in subset." ) );

 if( identical( v_ShutDownCost , subset , values ) )  // if nothing changes
  return;                                            // return

 // the shut-down variable t is in position ( f_time_horizon - init_t ) +
 // ( t - init_t ): those in the range [ 0 , init_t ) do not exist and their
 // cost cannot be changed
 if( subset.front() < init_t )
  throw( std::invalid_argument(
   "ThermalUnitBlock::set_shutdown_costs: invalid starting index in "
   "subset." ) );

 if( not_dry_run( issuePMod ) )
  // Change the physical representation
  assign( v_ShutDownCost , subset , values );

 if( not_dry_run( issueAMod ) && objective_generated() &&
     f_shut_down_in_obj ) {
  // Change the abstract representation
  Subset tmps = subset_sbtrct( subset , 2 * init_t );
  tmps = subset_add( tmps , f_time_horizon );
  DQuadFunction::Vec_FunctionValue tmpv( values , values + subset.size() );
  for( auto & c : tmpv )  // the Objective carries f_scale times the cost
   c *= f_scale;
  static_cast< DQuadFunction * >( objective.get_function()
  )->modify_linear_coefficients( std::move( tmpv ) , std::move( tmps ) ,
                                 true , un_ModBlock( issueAMod ) );
  }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockSbstMod >(
                            this , ThermalUnitBlockMod::eSetSDC ,
                            std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_shutdown_costs( subset ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_shutdown_costs( MF_dbl_it values ,
                                           Range rng ,
                                           ModParam issuePMod ,
                                           ModParam issueAMod )
{
 rng.second = std::min( rng.second , f_time_horizon );
 if( rng.second <= rng.first )
  return;

 c_Index sz = rng.second - rng.first;

 // a unit that was built with no shut-down cost has no shut-down variable
 // in its Objective [see generate_objective()], hence its costs can only
 // stay zero (checked before anything is changed)
 if( objective_generated() && ( ! f_shut_down_in_obj ) &&
     std::any_of( values , values + sz ,
                  []( double cst ) { return( cst != 0 ); } ) )
  throw( std::logic_error( "ThermalUnitBlock::set_shutdown_costs: the "
                           "Objective has no shut-down variables, the unit "
                           "having been generated with no shut-down cost" ) );

 if( v_ShutDownCost.empty() ) {
  if( std::all_of( values ,
                   values + sz ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  v_ShutDownCost.assign( f_time_horizon , 0 );
  }

 // If nothing changes, return
 if( std::equal( values , values + sz , v_ShutDownCost.begin() + rng.first ) )
  return;


 // see the comment in the Subset version
 if( rng.first < init_t )
  throw( std::invalid_argument( "ThermalUnitBlock::set_shutdown_costs: "
                                "invalid starting index in range." ) );

 if( not_dry_run( issuePMod ) )
  // Change the physical representation
  std::copy( values , values + sz , v_ShutDownCost.begin() + rng.first );

 if( not_dry_run( issueAMod ) && objective_generated() &&
     f_shut_down_in_obj ) {
  // Change the abstract representation
  c_Index shift = f_time_horizon - 2 * init_t;
  DQuadFunction::Vec_FunctionValue tmpv( values , values + sz );
  for( auto & c : tmpv )  // the Objective carries f_scale times the cost
   c *= f_scale;
  static_cast< DQuadFunction * >( objective.get_function()
  )->modify_linear_coefficients( std::move( tmpv ) ,
                                 Range( rng.first + shift ,
                                        rng.second + shift ) ,
                                 un_ModBlock( issueAMod ) );
  }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockRngdMod >(
                            this , ThermalUnitBlockMod::eSetSDC , rng ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_shutdown_costs( range ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_const_term( MF_dbl_it values ,
                                       Subset && subset ,
                                       const bool ordered ,
                                       ModParam issuePMod ,
                                       ModParam issueAMod )
{
 if( subset.empty() )
  return;

 if( v_ConstTerm.empty() ) {
  if( std::all_of( values ,
                   values + subset.size() ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  v_ConstTerm.assign( f_time_horizon , 0 );
 }

 std::vector< double > sorted;  // the values in the order of the subset
 if( ! ordered )
  sort_by_index( subset , values , sorted );

 if( subset.back() >= v_ConstTerm.size() )
  throw( std::invalid_argument(
   "ThermalUnitBlock::set_const_term: invalid index in subset." ) );

 if( identical( v_ConstTerm , subset , values ) )  // if nothing changes
  return;                                               // return

 if( not_dry_run( issuePMod ) )
  // Change the physical representation
  assign( v_ConstTerm , subset , values );

 if( not_dry_run( issueAMod ) && objective_generated() ) {
  // Change the abstract representation
  // the order of the variables in the Objective Function is:
  //
  // - first f_time_horizon - init_t start-up variables
  //
  // - then f_time_horizon active power variables
  //
  // - then f_time_horizon commitment variables
  //
  // - then possibly the rest
  //
  // hence, the commitment variables, whose coefficient is the fixed
  // cost, start from position 2 * f_time_horizon - init_t
  const Index dpos = 2 * f_time_horizon - init_t + shut_down_offset();

  Subset tmps = subset_add( subset , dpos );
  DQuadFunction::Vec_FunctionValue tmpv( values , values + subset.size() );
  for( auto & c : tmpv )  // the Objective carries f_scale times the cost
   c *= f_scale;
  static_cast< DQuadFunction * >( objective.get_function()
  )->modify_linear_coefficients( std::move( tmpv ) , std::move( tmps ) ,
                                 true , un_ModBlock( issueAMod ) );
 }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockSbstMod >(
                            this , ThermalUnitBlockMod::eSetConstT ,
                            std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_const_term( subset ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_const_term( MF_dbl_it values ,
                                       Range rng ,
                                       ModParam issuePMod ,
                                       ModParam issueAMod )
{
 rng.second = std::min( rng.second , f_time_horizon );
 if( rng.second <= rng.first )
  return;

 c_Index sz = rng.second - rng.first;
 if( v_ConstTerm.empty() ) {
  if( std::all_of( values ,
                   values + sz ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  v_ConstTerm.assign( f_time_horizon , 0 );
 }

 // If nothing changes, return
 if( std::equal( values , values + sz , v_ConstTerm.begin() + rng.first ) )
  return;

 if( not_dry_run( issuePMod ) )
  // Change the physical representation
  std::copy( values ,
             values + sz ,
             v_ConstTerm.begin() + rng.first );

 if( not_dry_run( issueAMod ) && objective_generated() ) {
  // Change the abstract representation
  // the order of the variables in the Objective Function is:
  //
  // - first f_time_horizon - init_t start-up variables
  //
  // - then f_time_horizon active power variables
  //
  // - then f_time_horizon commitment variables
  //
  // - then possibly the rest
  //
  // hence, the commitment variables, whose coefficient is the fixed
  // cost, start from position 2 * f_time_horizon - init_t
  const Index dpos = 2 * f_time_horizon - init_t + shut_down_offset();

  DQuadFunction::Vec_FunctionValue tmpv( values , values + sz );
  for( auto & c : tmpv )  // the Objective carries f_scale times the cost
   c *= f_scale;
  static_cast< DQuadFunction * >( objective.get_function()
  )->modify_linear_coefficients( std::move( tmpv ) ,
                                 Range( rng.first + dpos ,
                                        rng.second + dpos ) ,
                                 un_ModBlock( issueAMod ) );
 }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockRngdMod >(
                            this , ThermalUnitBlockMod::eSetConstT , rng ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_const_term( range ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_linear_term( MF_dbl_it values ,
                                        Subset && subset ,
                                        const bool ordered ,
                                        ModParam issuePMod ,
                                        ModParam issueAMod )
{
 if( subset.empty() )
  return;

 if( v_LinearTerm.empty() ) {
  if( std::all_of( values ,
                   values + subset.size() ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  v_LinearTerm.assign( f_time_horizon , 0 );
 }

 std::vector< double > sorted;  // the values in the order of the subset
 if( ! ordered )
  sort_by_index( subset , values , sorted );

 if( subset.back() >= v_LinearTerm.size() )
  throw( std::invalid_argument(
   "ThermalUnitBlock::set_linear_term: invalid index in subset." ) );

 if( identical( v_LinearTerm , subset , values ) )  // if nothing changes
  return;                                                // return

 if( not_dry_run( issuePMod ) )
  // Change the physical representation
  assign( v_LinearTerm , subset , values );

 if( not_dry_run( issueAMod ) && objective_generated() ) {
  // Change the abstract representation
  // the order of the variables in the Objective Function is:
  //
  // - first f_time_horizon - init_t start-up variables
  //
  // - then f_time_horizon active power variables
  //
  // - then possibly the rest
  //
  // hence, the active power  variables, whose coefficient is the linear
  // term of the cost, start from position f_time_horizon - init_t
  const Index dpos = f_time_horizon - init_t + shut_down_offset();

  Subset tmps = subset_add( subset , dpos );
  DQuadFunction::Vec_FunctionValue tmpv( values , values + subset.size() );
  for( auto & c : tmpv )  // the Objective carries f_scale times the cost
   c *= f_scale;
  static_cast< DQuadFunction * >( objective.get_function()
  )->modify_linear_coefficients( std::move( tmpv ) , std::move( tmps ) ,
                                 true , un_ModBlock( issueAMod ) );
 }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockSbstMod >(
                            this , ThermalUnitBlockMod::eSetLinT ,
                            std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_linear_term( subset ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_linear_term( MF_dbl_it values ,
                                        Range rng ,
                                        ModParam issuePMod ,
                                        ModParam issueAMod )
{
 rng.second = std::min( rng.second , f_time_horizon );
 if( rng.second <= rng.first )
  return;

 c_Index sz = rng.second - rng.first;
 if( v_LinearTerm.empty() ) {
  if( std::all_of( values ,
                   values + sz ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  v_LinearTerm.assign( f_time_horizon , 0 );
 }

 // If nothing changes, return
 if( std::equal( values , values + sz , v_LinearTerm.begin() + rng.first ) )
  return;

 if( not_dry_run( issuePMod ) )
  // Change the physical representation
  std::copy( values , values + sz , v_LinearTerm.begin() + rng.first );

 if( not_dry_run( issueAMod ) && objective_generated() ) {
  // Change the abstract representation
  // the order of the variables in the Objective Function is:
  //
  // - first f_time_horizon - init_t start-up variables
  //
  // - then f_time_horizon active power variables
  //
  // - then possibly the rest
  //
  // hence, the active power variables, whose coefficient is the linear
  // term of the cost, start from position f_time_horizon - init_t
  const Index dpos = f_time_horizon - init_t + shut_down_offset();

  DQuadFunction::Vec_FunctionValue tmpv( values , values + sz );
  for( auto & c : tmpv )  // the Objective carries f_scale times the cost
   c *= f_scale;
  static_cast< DQuadFunction * >( objective.get_function()
  )->modify_linear_coefficients( std::move( tmpv ) ,
                                 Range( rng.first + dpos ,
                                        rng.second + dpos ) ,
                                 un_ModBlock( issueAMod ) );
 }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockRngdMod >(
                            this , ThermalUnitBlockMod::eSetLinT , rng ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_linear_term( range ) )

/*--------------------------------------------------------------------------*/
// the reactive power variables q[t] are the LAST section of the Objective
// (see generate_objective): they are added only when f_reactive_power, after
// every other section, so they occupy the active-variable indices
// [ num_active_var - f_time_horizon , num_active_var ). The DP solvers price
// q[t] over [Qmin,Qmax] using this coefficient.

void ThermalUnitBlock::set_reactive_linear_term( MF_dbl_it values ,
                                                 Subset && subset ,
                                                 const bool ordered ,
                                                 ModParam issuePMod ,
                                                 ModParam issueAMod )
{
 if( subset.empty() )
  return;

 if( v_ReactiveLinearTerm.empty() ) {
  if( std::all_of( values ,
                   values + subset.size() ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  v_ReactiveLinearTerm.assign( f_time_horizon , 0 );
 }

 std::vector< double > sorted;  // the values in the order of the subset
 if( ! ordered )
  sort_by_index( subset , values , sorted );

 if( subset.back() >= v_ReactiveLinearTerm.size() )
  throw( std::invalid_argument(
   "ThermalUnitBlock::set_reactive_linear_term: invalid index in subset." ) );

 if( identical( v_ReactiveLinearTerm , subset , values ) )  // nothing changes
  return;

 if( not_dry_run( issuePMod ) )
  assign( v_ReactiveLinearTerm , subset , values );

 if( not_dry_run( issueAMod ) && objective_generated() && f_reactive_power ) {
  auto * qf = static_cast< DQuadFunction * >( objective.get_function() );
  const Index dpos = reactive_objective_start( qf );
  Subset tmps = subset_add( subset , dpos );
  DQuadFunction::Vec_FunctionValue tmpv( values , values + subset.size() );
  for( auto & c : tmpv )  // the Objective carries f_scale times the cost
   c *= f_scale;
  qf->modify_linear_coefficients( std::move( tmpv ) , std::move( tmps ) ,
                                  true , un_ModBlock( issueAMod ) );
 }

 if( issue_pmod( issuePMod ) )
  Block::add_Modification( std::make_shared< ThermalUnitBlockSbstMod >(
                            this , ThermalUnitBlockMod::eSetReactiveLinT ,
                            std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_reactive_linear_term( subset ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_reactive_linear_term( MF_dbl_it values ,
                                                 Range rng ,
                                                 ModParam issuePMod ,
                                                 ModParam issueAMod )
{
 rng.second = std::min( rng.second , f_time_horizon );
 if( rng.second <= rng.first )
  return;

 c_Index sz = rng.second - rng.first;
 if( v_ReactiveLinearTerm.empty() ) {
  if( std::all_of( values ,
                   values + sz ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  v_ReactiveLinearTerm.assign( f_time_horizon , 0 );
 }

 if( std::equal( values , values + sz ,
                 v_ReactiveLinearTerm.begin() + rng.first ) )
  return;

 if( not_dry_run( issuePMod ) )
  std::copy( values , values + sz , v_ReactiveLinearTerm.begin() + rng.first );

 if( not_dry_run( issueAMod ) && objective_generated() && f_reactive_power ) {
  auto * qf = static_cast< DQuadFunction * >( objective.get_function() );
  const Index dpos = reactive_objective_start( qf );
  DQuadFunction::Vec_FunctionValue tmpv( values , values + sz );
  for( auto & c : tmpv )  // the Objective carries f_scale times the cost
   c *= f_scale;
  qf->modify_linear_coefficients( std::move( tmpv ) ,
                                  Range( rng.first + dpos , rng.second + dpos ) ,
                                  un_ModBlock( issueAMod ) );
 }

 if( issue_pmod( issuePMod ) )
  Block::add_Modification( std::make_shared< ThermalUnitBlockRngdMod >(
                            this , ThermalUnitBlockMod::eSetReactiveLinT ,
                            rng ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_reactive_linear_term( range ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_quad_term( MF_dbl_it values ,
                                      Subset && subset ,
                                      const bool ordered ,
                                      ModParam issuePMod ,
                                      ModParam issueAMod )
{
 if( subset.empty() )
  return;

 if( v_QuadTerm.empty() ) {
  if( std::all_of( values ,
                   values + subset.size() ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  v_QuadTerm.assign( f_time_horizon , 0 );
 }

 std::vector< double > sorted;  // the values in the order of the subset
 if( ! ordered )
  sort_by_index( subset , values , sorted );

 if( subset.back() >= v_QuadTerm.size() )
  throw( std::invalid_argument(
   "ThermalUnitBlock::set_quad_term: invalid index in subset." ) );

 if( identical( v_QuadTerm , subset , values ) )  // if nothing changes
  return;                                              // return

 if( not_dry_run( issuePMod ) )
  // Change the physical representation
  assign( v_QuadTerm , subset , values );

 if( not_dry_run( issueAMod ) && objective_generated() ) {
  // the Objective carries f_scale times the cost
  DQuadFunction::Vec_FunctionValue svalues( values , values + subset.size() );
  for( auto & c : svalues )
   c *= f_scale;

  // Change the abstract representation
  // the order of the variables in the Objective Function is:
  //
  // - first f_time_horizon - init_t start-up variables
  //
  // - then f_time_horizon active power variables (which may have the
  //   nonzero quadratic cost coefficient, while the others do not)
  //
  // - then possibly the rest
  //
  // hence, if no perspective cuts are used, then the active power variables,
  // whose quadratic coefficient is the quadratic term of the cost, start
  // from position f_time_horizon - init_t, else the quadratic coefficient
  // becomes the *linear* coefficient of the perspective-cut variables,
  // whose section starts at cut_section_start(); the section is indexed by
  // time for the tbin / T / pt / SUSD formulations, and by the arcs of the
  // disaggregated graph (v_Z_h_k / v_Z_h / v_Z_k, whose .first is the time
  // instant) for the DP / SU / SD formulations

  if( ! ( AR & PCuts ) ) {

   Subset tmps = subset_add( subset ,
                             f_time_horizon - init_t + shut_down_offset() );
   DQuadFunction::Vec_FunctionValue tmplv( subset.size() , 0 );
   if( ! v_LinearTerm.empty() ) {
    auto tmplvit = tmplv.begin();
    for( auto t : subset )
     *( tmplvit++ ) = f_scale * v_LinearTerm[ t ];
   }

   static_cast< DQuadFunction * >( objective.get_function()
   )->modify_terms( svalues.begin() , tmplv.begin() , std::move( tmps ) ,
                    true , un_ModBlock( issueAMod ) );

  } else {

   const Index dpos = cut_section_start();
   const auto form = AR & FormMsk;

   if( ( form == tbinForm ) || ( form == TForm ) ||
       ( form == ptForm ) || ( form == SUSDForm ) ) {
    // time-indexed cut variables
    Subset tmps = subset_add( subset , dpos );
    DQuadFunction::Vec_FunctionValue tmplv( svalues );
    static_cast< DQuadFunction * >( objective.get_function()
    )->modify_linear_coefficients( std::move( tmplv ) , std::move( tmps ) ,
                                   true , un_ModBlock( issueAMod ) );
   } else {
    // arc-indexed cut variables: change every arc whose time instant
    // belongs to the (sorted) subset
    auto arc_change = [ & ]( const auto & Z ) {
     Subset nms;
     DQuadFunction::Vec_FunctionValue tmplv;
     for( Index i = 0 ; i < Z.size() ; ++i ) {
      auto it = std::lower_bound( subset.begin() , subset.end() ,
                                  Z[ i ].first );
      if( ( it != subset.end() ) && ( *it == Z[ i ].first ) ) {
       nms.push_back( dpos + i );
       tmplv.push_back( *( svalues.begin() +
                           std::distance( subset.begin() , it ) ) );
       }
      }
     if( ! nms.empty() )
      static_cast< DQuadFunction * >( objective.get_function()
      )->modify_linear_coefficients( std::move( tmplv ) , std::move( nms ) ,
                                     true , un_ModBlock( issueAMod ) );
     };

    if( form == DPForm )
     arc_change( v_Z_h_k );
    else
     if( form == SUForm )
      arc_change( v_Z_h );
     else
      arc_change( v_Z_k );
   }
  }
 }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockSbstMod >(
                            this , ThermalUnitBlockMod::eSetQuadT ,
                            std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_quad_term( subset ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_quad_term( MF_dbl_it values ,
                                      Range rng ,
                                      ModParam issuePMod ,
                                      ModParam issueAMod )
{
 rng.second = std::min( rng.second , f_time_horizon );
 if( rng.second <= rng.first )
  return;

 c_Index sz = rng.second - rng.first;
 if( v_QuadTerm.empty() ) {
  if( std::all_of( values ,
                   values + sz ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  v_QuadTerm.assign( f_time_horizon , 0 );
 }

 // If nothing changes, return
 if( std::equal( values , values + sz , v_QuadTerm.begin() + rng.first ) )
  return;

 if( not_dry_run( issuePMod ) )
  // Change the physical representation
  std::copy( values , values + sz , v_QuadTerm.begin() + rng.first );

 if( not_dry_run( issueAMod ) && objective_generated() ) {
  // the Objective carries f_scale times the cost
  DQuadFunction::Vec_FunctionValue svalues( values , values + sz );
  for( auto & c : svalues )
   c *= f_scale;

  // Change the abstract representation
  // the order of the variables in the Objective Function is:
  //
  // - first f_time_horizon - init_t start-up variables
  //
  // - then f_time_horizon active power variables (which may have the
  //   nonzero quadratic cost coefficient, while the others do not)
  //
  // - then possibly the rest
  //
  // hence, if no perspective cuts are used, then the active power variables,
  // whose quadratic coefficient is the quadratic term of the cost, start
  // from position f_time_horizon - init_t, else the quadratic coefficient
  // becomes the *linear* coefficient of the perspective-cut variables,
  // whose section starts at cut_section_start(); the section is indexed by
  // time for the tbin / T / pt / SUSD formulations, and by the arcs of the
  // disaggregated graph (v_Z_h_k / v_Z_h / v_Z_k, whose .first is the time
  // instant) for the DP / SU / SD formulations

  if( ! ( AR & PCuts ) ) {

   const Index dpos = f_time_horizon - init_t + shut_down_offset();
   DQuadFunction::Vec_FunctionValue tmplv( sz , 0 );
   if( ! v_LinearTerm.empty() )
    std::transform( v_LinearTerm.begin() + rng.first ,
                    v_LinearTerm.begin() + rng.second , tmplv.begin() ,
                    [ this ]( double c ) { return( f_scale * c ); } );

   static_cast< DQuadFunction * >( objective.get_function()
   )->modify_terms( svalues.begin() , tmplv.begin() ,
                    Range( rng.first + dpos , rng.second + dpos ) ,
                    un_ModBlock( issueAMod ) );

  } else {

   const Index dpos = cut_section_start();
   const auto form = AR & FormMsk;

   if( ( form == tbinForm ) || ( form == TForm ) ||
       ( form == ptForm ) || ( form == SUSDForm ) ) {
    // time-indexed cut variables
    DQuadFunction::Vec_FunctionValue tmplv( svalues );
    static_cast< DQuadFunction * >( objective.get_function()
    )->modify_linear_coefficients( std::move( tmplv ) ,
                                   Range( rng.first + dpos ,
                                          rng.second + dpos ) ,
                                   un_ModBlock( issueAMod ) );
   } else {
    // arc-indexed cut variables: change every arc whose time instant
    // falls in the range
    auto arc_change = [ & ]( const auto & Z ) {
     Subset nms;
     DQuadFunction::Vec_FunctionValue tmplv;
     for( Index i = 0 ; i < Z.size() ; ++i )
      if( ( Z[ i ].first >= rng.first ) && ( Z[ i ].first < rng.second ) ) {
       nms.push_back( dpos + i );
       tmplv.push_back( *( svalues.begin() + ( Z[ i ].first - rng.first ) ) );
       }
     if( ! nms.empty() )
      static_cast< DQuadFunction * >( objective.get_function()
      )->modify_linear_coefficients( std::move( tmplv ) , std::move( nms ) ,
                                     true , un_ModBlock( issueAMod ) );
     };

    if( form == DPForm )
     arc_change( v_Z_h_k );
    else
     if( form == SUForm )
      arc_change( v_Z_h );
     else
      arc_change( v_Z_k );
   }
  }
 }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockRngdMod >(
                            this , ThermalUnitBlockMod::eSetQuadT , rng ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_quad_term( range ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_primary_spinning_reserve_cost( MF_dbl_it values ,
                                                          Subset && subset ,
                                                          const bool ordered ,
                                                          ModParam issuePMod ,
                                                          ModParam issueAMod )
{
 // after the Variable are generated, a unit without primary reserve
 // Variable ignores the cost; before, the cost is kept, whether or not the
 // Variable will be there
 if( variables_generated() && v_primary_spinning_reserve.empty() )
  return;

 if( subset.empty() )
  return;

 // a unit whose Objective has been generated without this reserve term
 // [see generate_objective()] can only have zero costs
 if( objective_generated() && ( ! f_primary_in_obj ) &&
     std::any_of( values , values + subset.size() ,
                  []( double cst ) { return( cst != 0 ); } ) )
  throw( std::logic_error( "ThermalUnitBlock::set_primary_spinning_"
                           "reserve_cost: the Objective has no "
                           "primary reserve term" ) );

 if( v_PrimarySpinningReserveCost.empty() ) {
  // The primary spinning reserve costs are currently all zero.
  if( std::all_of( values ,
                   values + subset.size() ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;  // The given values are zero: nothing to do

  v_PrimarySpinningReserveCost.assign( f_time_horizon , 0 );
  }

 std::vector< double > sorted;  // the values in the order of the subset
 if( ! ordered )
  sort_by_index( subset , values , sorted );

 if( subset.back() >= v_PrimarySpinningReserveCost.size() )
  throw( std::invalid_argument(
   "ThermalUnitBlock::set_primary_spinning_reserve_cost: "
   "invalid index in subset." ) );

 if( identical( v_PrimarySpinningReserveCost , subset , values ) )
  return;

 if( not_dry_run( issuePMod ) )
  // Change the physical representation
  assign( v_PrimarySpinningReserveCost , subset , values );

 if( not_dry_run( issueAMod ) && objective_generated() &&
     f_primary_in_obj ) {
  // Change the abstract representation
  const Index dpos = reserve_section_start( false );

  Subset tmps = subset_add( subset , dpos );
  DQuadFunction::Vec_FunctionValue tmpv( values , values + subset.size() );
  for( auto & c : tmpv )  // the Objective carries f_scale times the cost
   c *= f_scale;
  static_cast< DQuadFunction * >( objective.get_function()
  )->modify_linear_coefficients( std::move( tmpv ) , std::move( tmps ) ,
                                 true , un_ModBlock( issueAMod ) );
 }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockSbstMod >(
                            this , ThermalUnitBlockMod::eSetPrSpResCost ,
                            std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_primary_spinning_reserve_cost( subset ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_primary_spinning_reserve_cost( MF_dbl_it values ,
                                                          Range rng ,
                                                          ModParam issuePMod ,
                                                          ModParam issueAMod )
{
 // after the Variable are generated, a unit without primary reserve
 // Variable ignores the cost; before, the cost is kept, whether or not the
 // Variable will be there
 if( variables_generated() && v_primary_spinning_reserve.empty() )
  return;

 rng.second = std::min( rng.second , f_time_horizon );
 if( rng.second <= rng.first )
  return;  // Empty range. Return.

 c_Index sz = rng.second - rng.first;

 // a unit whose Objective has been generated without this reserve term
 // [see generate_objective()] can only have zero costs
 if( objective_generated() && ( ! f_primary_in_obj ) &&
     std::any_of( values , values + sz ,
                  []( double cst ) { return( cst != 0 ); } ) )
  throw( std::logic_error( "ThermalUnitBlock::set_primary_spinning_"
                           "reserve_cost: the Objective has no "
                           "primary reserve term" ) );

 if( v_PrimarySpinningReserveCost.empty() ) {
  // The primary spinning reserve costs are currently all zero.
  if( std::all_of( values ,
                   values + sz ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;  // The given values are zero. So, there is nothing to be changed.

  v_PrimarySpinningReserveCost.assign( f_time_horizon , 0 );
 }

 // If nothing changes, return
 if( std::equal( values ,
                 values + sz ,
                 v_PrimarySpinningReserveCost.begin() + rng.first ) )
  return;

 if( not_dry_run( issuePMod ) )
  // Change the physical representation
  std::copy( values ,
             values + sz ,
             v_PrimarySpinningReserveCost.begin() + rng.first );

 if( not_dry_run( issueAMod ) && objective_generated() &&
     f_primary_in_obj ) {
  // Change the abstract representation
  const Index dpos = reserve_section_start( false );

  DQuadFunction::Vec_FunctionValue tmpv( values , values + sz );
  for( auto & c : tmpv )  // the Objective carries f_scale times the cost
   c *= f_scale;
  static_cast< DQuadFunction * >( objective.get_function()
  )->modify_linear_coefficients( std::move( tmpv ) ,
                                 Range( rng.first + dpos ,
                                        rng.second + dpos ) ,
                                 un_ModBlock( issueAMod ) );
 }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockRngdMod >(
                            this , ThermalUnitBlockMod::eSetPrSpResCost , rng ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_primary_spinning_reserve_cost( range ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_secondary_spinning_reserve_cost( MF_dbl_it values ,
                                                            Subset && subset ,
                                                            const bool ordered ,
                                                            ModParam issuePMod ,
                                                            ModParam issueAMod )
{
 // after the Variable are generated, a unit without secondary reserve
 // Variable ignores the cost; before, the cost is kept, whether or not the
 // Variable will be there
 if( variables_generated() && v_secondary_spinning_reserve.empty() )
  return;

 if( subset.empty() )
  return;

 // a unit whose Objective has been generated without this reserve term
 // [see generate_objective()] can only have zero costs
 if( objective_generated() && ( ! f_secondary_in_obj ) &&
     std::any_of( values , values + subset.size() ,
                  []( double cst ) { return( cst != 0 ); } ) )
  throw( std::logic_error( "ThermalUnitBlock::set_secondary_spinning_"
                           "reserve_cost: the Objective has no "
                           "secondary reserve term" ) );

 if( v_SecondarySpinningReserveCost.empty() ) {
  // The secondary spinning reserve costs are currently all zero.
  if( std::all_of( values ,
                   values + subset.size() ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;  // The given values are zero: nothing to do

  v_SecondarySpinningReserveCost.assign( f_time_horizon , 0 );
 }

 std::vector< double > sorted;  // the values in the order of the subset
 if( ! ordered )
  sort_by_index( subset , values , sorted );

 if( subset.back() >= v_SecondarySpinningReserveCost.size() )
  throw( std::invalid_argument( "ThermalUnitBlock::set_secondary_spinning_"
                                "reserve_cost: invalid index in subset." ) );

 if( identical( v_SecondarySpinningReserveCost , subset , values ) )
  return;

 if( not_dry_run( issuePMod ) )
  // Change the physical representation
  assign( v_SecondarySpinningReserveCost , subset , values );

 if( not_dry_run( issueAMod ) && objective_generated() &&
     f_secondary_in_obj ) {
  // Change the abstract representation
  const Index dpos = reserve_section_start( true );

  Subset tmps = subset_add( subset , dpos );
  DQuadFunction::Vec_FunctionValue tmpv( values , values + subset.size() );
  for( auto & c : tmpv )  // the Objective carries f_scale times the cost
   c *= f_scale;
  static_cast< DQuadFunction * >( objective.get_function()
  )->modify_linear_coefficients( std::move( tmpv ) , std::move( tmps ) ,
                                 true , un_ModBlock( issueAMod ) );
 }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockSbstMod >(
                            this , ThermalUnitBlockMod::eSetSecSpResCost ,
                            std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_secondary_spinning_reserve_cost( subset ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_secondary_spinning_reserve_cost( MF_dbl_it values ,
                                                            Range rng ,
                                                            ModParam issuePMod ,
                                                            ModParam issueAMod )
{
 // after the Variable are generated, a unit without secondary reserve
 // Variable ignores the cost; before, the cost is kept, whether or not the
 // Variable will be there
 if( variables_generated() && v_secondary_spinning_reserve.empty() )
  return;

 rng.second = std::min( rng.second , f_time_horizon );
 if( rng.second <= rng.first )
  return;  // Empty range. Return.

 c_Index sz = rng.second - rng.first;

 // a unit whose Objective has been generated without this reserve term
 // [see generate_objective()] can only have zero costs
 if( objective_generated() && ( ! f_secondary_in_obj ) &&
     std::any_of( values , values + sz ,
                  []( double cst ) { return( cst != 0 ); } ) )
  throw( std::logic_error( "ThermalUnitBlock::set_secondary_spinning_"
                           "reserve_cost: the Objective has no "
                           "secondary reserve term" ) );

 if( v_SecondarySpinningReserveCost.empty() ) {
  // The secondary spinning reserve costs are currently all zero.
  if( std::all_of( values ,
                   values + sz ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;  // The given values are zero. So, there is nothing to be changed.

  v_SecondarySpinningReserveCost.assign( f_time_horizon , 0 );
 }

 // If nothing changes, return
 if( std::equal( values ,
                 values + sz ,
                 v_SecondarySpinningReserveCost.begin() + rng.first ) )
  return;

 if( not_dry_run( issuePMod ) )
  // Change the physical representation
  std::copy( values ,
             values + sz ,
             v_SecondarySpinningReserveCost.begin() + rng.first );

 if( not_dry_run( issueAMod ) && objective_generated() &&
     f_secondary_in_obj ) {
  // Change the abstract representation
  const Index dpos = reserve_section_start( true );

  DQuadFunction::Vec_FunctionValue tmpv( values , values + sz );
  for( auto & c : tmpv )  // the Objective carries f_scale times the cost
   c *= f_scale;
  static_cast< DQuadFunction * >( objective.get_function()
  )->modify_linear_coefficients( std::move( tmpv ) ,
                                 Range( rng.first + dpos ,
                                        rng.second + dpos ) ,
                                 un_ModBlock( issueAMod ) );
 }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockRngdMod >(
                            this , ThermalUnitBlockMod::eSetSecSpResCost , rng ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_secondary_spinning_reserve_cost( range ) )

/*--------------------------------------------------------------------------*/

ThermalUnitBlock::Index ThermalUnitBlock::first_free_instant(
                                                         int init ) const
{
 Index t0;
 if( init > 0 )
  t0 = ( Index( init ) >= f_MinUpTime ? 0 : f_MinUpTime - Index( init ) );
 else
  t0 = ( Index( - init ) >= f_MinDownTime ? 0 :
         f_MinDownTime - Index( - init ) );

 // beyond the horizon the unit never switches within it: the commitment is
 // fixed everywhere and there is no start-up or shut-down Variable
 return( std::min( t0 , f_time_horizon ) );

 }  // end( ThermalUnitBlock::first_free_instant )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::guts_of_set_init_updown_time( int value ,
                                                     ModParam issuePMod ,
                                                     ModParam issueAMod )
{
 if( f_InitUpDownTime == value )
  return;  // nothing changes; return

 // the initial state decides the first free instant t0, hence which
 // Variable exist and which are fixed, and, through t0 and whether the
 // unit is on, every row: once the Variable are generated, it may change
 // only if neither changes, nor the bound that it puts on the minimum
 // times [see clamp_min_time()], of which the one at the bound of the old
 // state may have been cut by it and hence may differ in the data; the
 // rows are then written anew [see update_rows()], which finds them equal
 if( not_dry_run( issueAMod ) && variables_generated() ) {
  auto same_bound = [ & ]( Index m , bool up ) {
   const auto b_old = clamp_min_time( std::numeric_limits< Index >::max() ,
                                      f_time_horizon , f_InitUpDownTime ,
                                      up );
   const auto b_new = clamp_min_time( std::numeric_limits< Index >::max() ,
                                      f_time_horizon , value , up );
   return( ( b_old == b_new ) || ( ( m < b_old ) && ( m < b_new ) ) );
   };
  if( ( ( value > 0 ) != ( f_InitUpDownTime > 0 ) ) ||
      ( first_free_instant( value ) != init_t ) ||
      ( ! same_bound( f_MinUpTime , true ) ) ||
      ( ! same_bound( f_MinDownTime , false ) ) )
   throw( std::logic_error( "ThermalUnitBlock::set_init_updown_time: "
                            "InitUpDownTime " + std::to_string( value ) +
                            " changes the Variable of the unit (whether it "
                            "is on before the horizon, or the first instant "
                            "at which it may switch), which cannot change "
                            "once they are generated" ) );
  }

 const int old_value = f_InitUpDownTime;
 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation
  f_InitUpDownTime = value;
  compute_initial_ramp_steps();
  }

 if( not_dry_run( issueAMod ) && variables_generated() )
  // Change the abstract representation; if it cannot be, the physical one
  // is restored, so that the two do not disagree
  try {
   update_rows( un_ModBlock( issueAMod ) );
   }
  catch( ... ) {
   f_InitUpDownTime = old_value;
   compute_initial_ramp_steps();
   throw;
   }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockMod >(
                            this , ThermalUnitBlockMod::eSetInitUD ) ,
                           Observer::par2chnl( issuePMod ) );

 }  // end( ThermalUnitBlock::guts_of_set_init_updown_time )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_init_updown_time( MF_int_it values ,
                                             Subset && subset ,
                                             const bool ordered ,
                                             ModParam issuePMod ,
                                             ModParam issueAMod )
{
 if( subset.empty() )
  return;

 // Find the last index 0
 auto index_it = std::find( subset.rbegin() , subset.rend() , 0 );

 if( index_it == subset.rend() )
  return;  // 0 is not in subset; return

 std::advance( values , std::distance( index_it , subset.rend() ) - 1 );

 guts_of_set_init_updown_time( *values , issuePMod , issueAMod );

}  // end( ThermalUnitBlock::set_init_updown_time( subset ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_init_updown_time( MF_int_it values ,
                                             Range rng ,
                                             ModParam issuePMod ,
                                             ModParam issueAMod )
{
 rng.second = std::min( rng.second , static_cast< decltype( rng.second ) >( 1 ) );
 if( ! ( ( rng.first <= 0 ) && ( 0 < rng.second ) ) )
  return;  // 0 does not belong to the range; return

 std::advance( values , -rng.first );

 guts_of_set_init_updown_time( *values , issuePMod , issueAMod );

}  // end( ThermalUnitBlock::set_init_updown_time( range ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_min_up_down_time( Index min_up_time ,
                                             Index min_down_time ,
                                             ModParam issuePMod )
{
 if( variables_generated() )
  throw( std::logic_error( "ThermalUnitBlock::set_min_up_down_time: the "
                           "Variable have been generated already, and the "
                           "minimum times decide how many there are" ) );

 f_MinUpTime = clamp_min_time( min_up_time , f_time_horizon ,
                               f_InitUpDownTime , true );
 f_MinDownTime = clamp_min_time( min_down_time , f_time_horizon ,
                                 f_InitUpDownTime , false );

 if( issue_pmod( issuePMod ) )
  Block::add_Modification( std::make_shared< ThermalUnitBlockMod >(
                            this , ThermalUnitBlockMod::eSetInitUD ) ,
                           Observer::par2chnl( issuePMod ) );

 }  // end( ThermalUnitBlock::set_min_up_down_time )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::scale( MF_dbl_it values ,
                              Subset && subset ,
                              const bool ordered ,
                              c_ModParam issuePMod ,
                              c_ModParam issueAMod )
{
 if( subset.empty() )
  return;  // Since the given Subset is empty, no operation is performed

 if( f_scale == *values )
  return;  // The scale factor does not change: nothing to do

 if( not_dry_run( issuePMod ) ) {
  f_scale = *values;  // Update the scale factor

  if( not_dry_run( issueAMod ) ) {
   // Update the abstract representation
   if( objective_generated() )
    // Update the Objective
    update_objective( Range( 0 , Inf< Index >() ) , issuePMod , issueAMod );
  }
 }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< UnitBlockMod >(
                            this , UnitBlockMod::eScale ) ,
                           Observer::par2chnl( issuePMod ) );
 else if( auto f_Block = get_f_Block() )
  f_Block->add_Modification( std::make_shared< UnitBlockMod >(
                              this , UnitBlockMod::eScale ) ,
                             Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::scale )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::update_objective_start_up( const Subset & subset ,
                                                  c_ModParam issueAMod ) const
{
 if( ! objective_generated() )
  return;  // the Objective has not been generated: nothing to be done

 auto function = dynamic_cast< DQuadFunction * >( objective.get_function() );

 if( ! function )
  return;

 for( auto t : subset ) {
  if( t < init_t )
   continue;

  auto var_index = function->is_active( &v_start_up[ t - init_t ] );
  assert( var_index < function->get_num_active_var() );
  function->modify_linear_coefficient( var_index ,
                                       f_scale * v_StartUpCost[ t ] ,
                                       un_ModBlock( issueAMod ) );
 }
}  // end( ThermalUnitBlock::update_objective_start_up )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::update_objective_active_power( const Subset & subset ,
                                                      c_ModParam issueAMod )
{
 if( ! objective_generated() )
  return;  // the Objective has not been generated: nothing to be done

 auto function = dynamic_cast< DQuadFunction * >( objective.get_function() );

 if( ! function )
  return;

 // one term of the Objective per instant, hence one abstract Modification
 // each: they all go into a single GroupModification, so that a Solver able
 // to write a whole set of them in one operation does that instead of one
 // call per instant [see MILPSolver::process_group_modification()]
 auto nAM = un_ModBlock( make_par( par2mod( issueAMod ) ,
                                   open_channel( par2chnl( issueAMod ) ) ) );

 for( auto t : subset ) {
  auto var_index = function->is_active( &v_active_power[ t ] );
  assert( var_index < function->get_num_active_var() );
  function->modify_term( var_index ,
                         f_scale * v_LinearTerm[ t ] ,
                         AR & PCuts ? 0.0 : f_scale * v_QuadTerm[ t ] ,
                         nAM );
 }

 close_channel( par2chnl( nAM ) );
}  // end( ThermalUnitBlock::update_objective_active_power )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::update_objective_commitment( const Subset & subset ,
                                                    c_ModParam issueAMod ) const
{
 if( ! objective_generated() )
  return;  // the Objective has not been generated: nothing to be done

 auto function = dynamic_cast< DQuadFunction * >( objective.get_function() );

 if( ! function )
  return;

 for( auto t : subset ) {
  auto var_index = function->is_active( &v_commitment[ t ] );
  assert( var_index < function->get_num_active_var() );
  function->modify_linear_coefficient( var_index ,
                                       f_scale * v_ConstTerm[ t ] ,
                                       un_ModBlock( issueAMod ) );
 }
}  // end( ThermalUnitBlock::update_objective_commitment )

/*--------------------------------------------------------------------------*/

double ThermalUnitBlock::get_design_cost( void ) const
{
 if( f_InvestmentCost == 0 )
  return( 0 );  // no design variable: no design cost

 if( objective_generated() )
  if( auto function = dynamic_cast< DQuadFunction * >(
                                                  objective.get_function() ) ) {
   auto var_index = function->is_active( & design );
   if( var_index < function->get_num_active_var() )
    // the current (possibly dualized) coefficient of the design variable;
    // the Objective stores f_scale * cost, so divide it out to return the
    // cost in the same unscaled units as get_investment_cost()
    return( function->get_linear_coefficient( var_index ) / f_scale );
   }

 // the Objective is not available: fall back to the "original" cost
 return( f_InvestmentCost );

}  // end( ThermalUnitBlock::get_design_cost )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::update_objective_investment( ModParam issuePMod ,
                                                    ModParam issueAMod )
{
 if( ! objective_generated() )
  return;  // the Objective has not been generated: nothing to be done

 if( f_InvestmentCost == 0 )
  return;  // no design term in the Objective: nothing to be done

 if( not_dry_run( issueAMod ) ) {
  // change the abstract representation, i.e., the coefficient of the design
  // variable in the Objective
  auto function = dynamic_cast< DQuadFunction * >( objective.get_function() );
  if( ! function )
   return;
  auto var_index = function->is_active( & design );
  assert( var_index < function->get_num_active_var() );
  function->modify_linear_coefficient( var_index ,
                                       f_scale * f_InvestmentCost ,
                                       un_ModBlock( issueAMod ) );
  }

 // only issue the physical Modification if the design cost actually changed
 // since the last one was issued: a dualizing Solver rewrites the whole
 // coefficient vector (design included) at every iteration, but the design
 // coefficient itself changes only when the dual multiplier on it does. Issuing
 // a eSetInvCost every time would force the dualizing Bundle to invalidate this
 // component's linearizations at each iteration (preventing convergence) and
 // flood the (DP) Solvers with useless global-pool rechecks.
 const double cur = get_design_cost();
 if( ! ( cur == f_last_design_cost ) ) {  // ( cur != cached ), NaN-safe
  f_last_design_cost = cur;

  if( issue_pmod( issuePMod ) )
   // Issue a Physical Modification so that Solvers that consume the structural
   // data (the DP Solvers) refresh their copy of the design cost. Note that the
   // design cost lives only in the abstract Objective (there is no separate
   // physical field), so the Solvers re-read it via get_design_cost()
   Block::add_Modification( std::make_shared< ThermalUnitBlockMod >(
                             this , ThermalUnitBlockMod::eSetInvCost ) ,
                            Observer::par2chnl( issuePMod ) );
  }

}  // end( ThermalUnitBlock::update_objective_investment )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::update_objective( const Subset & subset ,
                                         ModParam issuePMod ,
                                         c_ModParam issueAMod ) {
 update_objective_investment( issuePMod , issueAMod );
 update_objective_start_up( subset , issueAMod );
 update_objective_active_power( subset , issueAMod );
 update_objective_commitment( subset , issueAMod );
 update_objective_other_terms( subset , issueAMod );
 update_objective_tail( subset , issueAMod );
}  // end( ThermalUnitBlock::update_objective( subset ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::update_objective_other_terms( const Subset & subset ,
                                                     c_ModParam issueAMod )
{
 if( ! objective_generated() )
  return;  // the Objective has not been generated: nothing to be done

 auto function = dynamic_cast< DQuadFunction * >( objective.get_function() );

 if( ! function )
  return;

 // the shut-down, reserve, perspective-cut and reactive power terms, each
 // f_scale times the cost of one copy as in generate_objective(): they all
 // go into a single Modification
 Subset nms;
 DQuadFunction::Vec_FunctionValue coeff;
 auto add = [ & ]( const ColVariable & var , double cost ) {
  const auto idx = function->is_active( & var );
  assert( idx < function->get_num_active_var() );
  nms.push_back( idx );
  coeff.push_back( f_scale * cost );
  };

 if( f_shut_down_in_obj )
  for( auto t : subset )
   if( t >= init_t )
    add( v_shut_down[ t - init_t ] , v_ShutDownCost[ t ] );

 if( f_primary_in_obj && ( ! v_PrimarySpinningReserveCost.empty() ) )
  for( auto t : subset )
   add( v_primary_spinning_reserve[ t ] , v_PrimarySpinningReserveCost[ t ] );

 if( f_secondary_in_obj && ( ! v_SecondarySpinningReserveCost.empty() ) )
  for( auto t : subset )
   add( v_secondary_spinning_reserve[ t ] ,
        v_SecondarySpinningReserveCost[ t ] );

 if( AR & PCuts ) {
  const auto form = AR & FormMsk;
  if( ( form == tbinForm ) || ( form == TForm ) || ( form == ptForm ) )
   for( auto t : subset )
    add( v_cut[ t ] , v_QuadTerm[ t ] );
  else
   if( form == SUSDForm )
    for( auto t : subset )
     add( v_cut_teta[ t ] , v_QuadTerm[ t ] );
   else {
    // arc-indexed cut variables, each carrying the quadratic term of the
    // time instant of its arc
    auto arcs = [ & ]( const auto & Z , const std::vector< ColVariable > & cut ) {
     for( Index i = 0 ; i < Z.size() ; ++i )
      if( std::binary_search( subset.begin() , subset.end() , Z[ i ].first ) )
       add( cut[ i ] , v_QuadTerm[ Z[ i ].first ] );
     };
    if( form == DPForm )
     arcs( v_Z_h_k , v_cut_h_k );
    else
     if( form == SUForm )
      arcs( v_Z_h , v_cut_h );
     else
      if( form == SDForm )
       arcs( v_Z_k , v_cut_k );
    }
  }

 if( ! v_RefSchedule.empty() )
  for( auto t : subset )
   add( v_abs_ref_schedule[ t ] , 1 );

 if( f_reactive_power && ( ! v_ReactiveLinearTerm.empty() ) )
  for( auto t : subset )
   add( v_reactive_power[ t ] , v_ReactiveLinearTerm[ t ] );

 if( ! nms.empty() )
  function->modify_linear_coefficients( std::move( coeff ) , std::move( nms ) ,
                                        false , un_ModBlock( issueAMod ) );

}  // end( ThermalUnitBlock::update_objective_other_terms )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::update_objective( Range rng ,
                                         ModParam issuePMod ,
                                         c_ModParam issueAMod ) {
 rng.second = std::min( rng.second , f_time_horizon );
 if( rng.second <= rng.first )
  return;

 Subset subset( rng.second - rng.first );
 std::iota( subset.begin() , subset.end() , rng.first );

 update_objective( subset , issuePMod , issueAMod );
}  // end( ThermalUnitBlock::update_objective( range ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::guts_of_add_Modification( p_Mod mod , ChnlName chnl )
{
 // process abstract Modification - - - - - - - - - - - - - - - - - - - - - -
 /* This requires to patiently sift through the possible Modification types
  * to find what this Modification exactly is and appropriately mirror the
  * changes to the "abstract representation" to the "physical one".
  *
  * Note that since ThermalUnitBlock is a "leaf" Block (has no sub-Block),
  * this method does not have to deal with GroupModification since these
  * are produced by Block::add_Modification(), but this method is called
  * *before* that one is.
  *
  * As an important consequence,
  *
  *   THE STATE OF THE DATA STRUCTURE IN ThermalUnitBlock WHEN THIS METHOD
  *   IS EXECUTED IS PRECISELY THE ONE IN WHICH THE Modification WAS
  *   ISSUED: NO COMPLICATED OPERATIONS (Variable AND/OR Constraint BEING
  *   ADDED/REMOVED ...) CAN HAVE BEEN PERFORMED IN THE MEANTIME
  *
  * This assumption drastically simplifies some logic here. */

 // VariableMod - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( const auto tmod = dynamic_cast< VariableMod * >( mod ) ) {
  // the only supported change is fixing/unfixing a Variable: the abstract
  // representation is already up to date (the fixed status lives in the
  // Variable itself, and it has no physical counterpart in the
  // ThermalUnitBlock data), so all that remains to be done is telling the
  // Solvers that consume the physical representation about the change; it
  // is then their business to deal with it (or to complain if they can't)
  if( Variable::is_fixed( tmod->old_state() ) ==
      Variable::is_fixed( tmod->new_state() ) )
   throw( std::logic_error( "ThermalUnitBlock::add_Modification: VariableMod "
			    "changing anything but the fixed status is not "
			    "supported" ) );

  if( anyone_there() )
   Block::add_Modification( std::make_shared< ThermalUnitBlockMod >(
			     this , ThermalUnitBlockMod::eFixVars ) , chnl );
  return;
 }

 // BlockMod - Generic modification - - - - - - - - - - - - - - - - - - - - -
 if( const auto tmod = dynamic_cast< BlockMod * >( mod ) ) {
  // changing the Objective is not supported, but the Modification is issued
  // when it is first set, in which case it must be ignored
  // THE Modification SHOULD NOT BE ISSUED WHEN IT IS CREATED!
  //if( ! objective_generated() )
  // return;

  throw( std::logic_error( "ThermalUnitBlock - BlockMod not supported." ) );

  // TODO: BlockMod - obj changed
  return;
  }

 // FunctionMod - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( const auto tmod = dynamic_cast< FunctionMod * >( mod ) ) {
  auto f = tmod->function();
  if( f == static_cast< FRealObjective * >( get_objective()
					    )->get_function() ) {
   handle_objective_change( tmod , chnl );
   return;
   }

  std::ostringstream em;
  em << *mod;
  throw( std::invalid_argument( "ThermalUnitBlock::add_Modification: "
				"unsupported " + em.str() ) );
  return;
  }

 // any other Modification is not supported - - - - - - - - - - - - - - - - -

 std::ostringstream em;
 em << *mod;
 throw( std::invalid_argument( "ThermalUnitBlock::add_Modification: "
			       "unsupported " + em.str() ) );

 }  // end( ThermalUnitBlock::guts_of_add_Modification )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::derive_start_up( const std::vector< double > & u ,
                                       std::vector< double > & su ,
                                       std::vector< double > & sd ) const
{
 const auto n = v_start_up.size();
 su.assign( n , 0.0 );
 sd.assign( n , 0.0 );
 if( ( ! n ) || ( u.size() < f_time_horizon ) )
  return;

 // the first instant is decided by the state the unit was in before the
 // horizon, the others by the profile itself
 if( init_t == 0 ) {
  su[ 0 ] = ( ( f_InitUpDownTime <= 0 ) && ( u[ 0 ] > 0.5 ) ) ? 1.0 : 0.0;
  sd[ 0 ] = ( ( f_InitUpDownTime > 0 ) && ( u[ 0 ] <= 0.5 ) ) ? 1.0 : 0.0;
  }

 for( Index t = std::max( init_t , Index( 1 ) ) ; t < f_time_horizon ; ++t ) {
  su[ t - init_t ] = ( ( u[ t ] > 0.5 ) && ( u[ t - 1 ] <= 0.5 ) ) ? 1.0 : 0.0;
  sd[ t - init_t ] = ( ( u[ t ] <= 0.5 ) && ( u[ t - 1 ] > 0.5 ) ) ? 1.0 : 0.0;
  }
 }  // end( ThermalUnitBlock::derive_start_up )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_solution( void )
{
 // canonical part: the caller has already set v_active_power[t] and
 // v_commitment[t] for all t. We read them back to derive the
 // formulation-specific auxiliaries below.
 auto Pi = get_const_active_power( 0 );
 auto Ci = get_const_commitment( 0 );
 if( ( ! Pi ) || ( ! Ci ) )
  return;             // no canonical variables: nothing to derive from

 // start_up[ t ] = 1 iff commitment goes off->on at t. start_up is
 // indexed from init_t onwards (size = time_horizon - init_t); the
 // boundary case t == init_t == 0 is handled via the pre-horizon
 // state in f_InitUpDownTime: if the unit was off before t = 0
 // (f_InitUpDownTime <= 0) and is on at t = 0, that counts as a
 // start-up at t = 0; otherwise start_up[ 0 ] = 0
 if( auto sup_it = get_start_up() ) {
  if( init_t == 0 )
   sup_it[ 0 ].set_value(
    ( ( f_InitUpDownTime <= 0 ) && ( Ci[ 0 ].get_value() > 0.5 ) )
    ? 1.0 : 0.0 );
  for( Index t = std::max( init_t , Index( 1 ) ) ;
       t < f_time_horizon ; ++t )
   sup_it[ t - init_t ].set_value(
    ( ( Ci[ t ].get_value() > 0.5 ) && ( Ci[ t - 1 ].get_value() <= 0.5 ) )
    ? 1.0 : 0.0 );
  }

 // shut_down[ t ] = 1 iff commitment goes on->off at t (symmetric)
 if( auto sdn_it = get_shut_down() ) {
  if( init_t == 0 )
   sdn_it[ 0 ].set_value(
    ( ( f_InitUpDownTime > 0 ) && ( Ci[ 0 ].get_value() <= 0.5 ) )
    ? 1.0 : 0.0 );
  for( Index t = std::max( init_t , Index( 1 ) ) ;
       t < f_time_horizon ; ++t )
   sdn_it[ t - init_t ].set_value(
    ( ( Ci[ t ].get_value() <= 0.5 ) && ( Ci[ t - 1 ].get_value() > 0.5 ) )
    ? 1.0 : 0.0 );
  }

 // the deviation from the reference schedule, if any
 for( Index t = 0 ; t < v_abs_ref_schedule.size() ; ++t )
  v_abs_ref_schedule[ t ].set_value( std::abs( Pi[ t ].get_value() -
                                               v_RefSchedule[ t ] ) );

 // the variables z_t of the perspective cuts, in every formulation that has
 // them: the perspective value p^2 / u, which is p^2 for u = 1 and 0 for
 // u = 0, satisfies every cut and gives a_t z_t = a_t p_t^2, the cost of the
 // quadratic term, also for a fractional (convex combination) commitment
 if( auto cut_it = get_cut() )
  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   const double pt = Pi[ t ].get_value();
   const double ut = Ci[ t ].get_value();
   cut_it[ t ].set_value( ( ut > 1e-6 ) ? pt * pt / ut : 0.0 );
   }

 // the path formulations (pt, DP, SU, SD, SUSD): the arcs of the path of
 // the schedule, and the power of each run on its own Variable
 if( ! v_commitment_plus.empty() )
  set_path_solution();

 }  // end( ThermalUnitBlock::set_solution )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_path_solution( void )
{
 auto Pi = get_const_active_power( 0 );
 auto Ci = get_const_commitment( 0 );
 const Index T = f_time_horizon;

 // the path is defined only for an integer commitment: a fractional one
 // (e.g., a convex combination of schedules) is a combination of paths
 // that it does not identify, and the arc variables are left as they are
 for( Index t = 0 ; t < T ; ++t ) {
  const double ut = Ci[ t ].get_value();
  if( std::min( std::abs( ut ) , std::abs( 1 - ut ) ) > 1e-6 )
   return;
  }

 // the runs of the schedule, with the indices of the arcs [see
 // generate_abstract_variables()]: an on-arc ( h , k ) is on at the
 // instants h - 1 , ... , k - 1 and shuts down at k, the run in progress
 // before the horizon having h = 0, and an off-arc ( k , h ) is off at the
 // instants k , ... , h - 2; a run that the end of the horizon cuts ends
 // at T + 1
 std::vector< std::pair< Index , Index > > on_arcs , off_arcs;
 bool on = ( f_InitUpDownTime > 0 );
 Index from = 0;
 for( Index t = 0 ; t < T ; ++t ) {
  const bool ont = ( Ci[ t ].get_value() > 0.5 );
  if( ont == on )
   continue;
  if( on ) {               // a shut-down at t
   on_arcs.emplace_back( from , t );
   from = t;
   }
  else {                   // a start-up at t
   off_arcs.emplace_back( from , t + 1 );
   from = t + 1;
   }
  on = ont;
  }
 if( on )
  on_arcs.emplace_back( from , T + 1 );
 else
  off_arcs.emplace_back( from , T + 1 );

 auto in = []( const std::vector< std::pair< Index , Index > > & arcs ,
               const std::pair< Index , Index > & a ) {
  return( std::find( arcs.begin() , arcs.end() , a ) != arcs.end() );
  };

 for( Index i = 0 ; i < v_commitment_plus.size() ; ++i )
  v_commitment_plus[ i ].set_value( in( on_arcs , v_Y_plus[ i ] ) ? 1 : 0 );
 for( Index i = 0 ; i < v_commitment_minus.size() ; ++i )
  v_commitment_minus[ i ].set_value( in( off_arcs , v_Y_minus[ i ] ) ? 1
                                                                      : 0 );

 // the power of the run that contains t, on the Variable of that run, and
 // its square on the Variable of the perspective cuts
 auto pw = [ & ]( Index t ) { return( Pi[ t ].get_value() ); };
 auto sq = [ & ]( Index t ) { return( pw( t ) * pw( t ) ); };
 auto run_hk = [ & ]( Index t , Index h , Index k ) {
  return( in( on_arcs , std::make_pair( h , k ) ) &&
          ( h <= t + 1 ) && ( t + 1 <= k ) );
  };
 auto run_h = [ & ]( Index t , Index h ) {  // started at h - 1, on at t
  for( const auto & a : on_arcs )
   if( ( a.first == h ) && ( h <= t + 1 ) && ( t + 1 <= a.second ) )
    return( true );
  return( false );
  };
 auto run_k = [ & ]( Index t , Index k ) {  // shut down at k, on at t
  for( const auto & a : on_arcs )
   if( ( a.second == k ) && ( a.first <= t + 1 ) && ( t + 1 <= k ) )
    return( true );
  return( false );
  };

 for( Index j = 0 ; j < v_active_power_h_k.size() ; ++j ) {
  const auto & e = v_P_h_k[ j ];
  v_active_power_h_k[ j ].set_value(
   run_hk( e.first , e.second.first , e.second.second ) ? pw( e.first ) : 0 );
  }
 for( Index j = 0 ; j < v_cut_h_k.size() ; ++j ) {
  const auto & e = v_Z_h_k[ j ];
  v_cut_h_k[ j ].set_value(
   run_hk( e.first , e.second.first , e.second.second ) ? sq( e.first ) : 0 );
  }
 for( Index j = 0 ; j < v_active_power_h.size() ; ++j )
  v_active_power_h[ j ].set_value(
   run_h( v_P_h[ j ].first , v_P_h[ j ].second ) ? pw( v_P_h[ j ].first )
                                                  : 0 );
 for( Index j = 0 ; j < v_cut_h.size() ; ++j )
  v_cut_h[ j ].set_value(
   run_h( v_Z_h[ j ].first , v_Z_h[ j ].second ) ? sq( v_Z_h[ j ].first )
                                                  : 0 );
 for( Index j = 0 ; j < v_active_power_k.size() ; ++j )
  v_active_power_k[ j ].set_value(
   run_k( v_P_k[ j ].first , v_P_k[ j ].second ) ? pw( v_P_k[ j ].first )
                                                  : 0 );
 for( Index j = 0 ; j < v_cut_k.size() ; ++j )
  v_cut_k[ j ].set_value(
   run_k( v_Z_k[ j ].first , v_Z_k[ j ].second ) ? sq( v_Z_k[ j ].first )
                                                  : 0 );
 for( Index t = 0 ; t < v_cut_teta.size() ; ++t )
  v_cut_teta[ t ].set_value( ( Ci[ t ].get_value() > 0.5 ) ? sq( t ) : 0 );

 }  // end( ThermalUnitBlock::set_path_solution )


/*--------------------------------------------------------------------------*/

Block::Index ThermalUnitBlock::cut_section_start( void ) const
{
 // the start-up, active power and commitment variables always come first
 // (see the layout in generate_objective()); the shut-down, the
 // schedule-deviation and the reserve sections are present only when
 // generate_objective() has put them there
 return( reserve_section_start( true ) +
         ( f_secondary_in_obj ? f_time_horizon : 0 ) );

}  // end( ThermalUnitBlock::cut_section_start )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::handle_objective_change( FunctionMod * mod ,
                                                ChnlName chnl )
{
 const auto * qf = static_cast< const DQuadFunction * >( mod->function() );
 auto par = make_par( eNoBlck , chnl );
 // the method is implemented by calling the physical change methods
 // set_startup_costs(), set_linear_term() etc.; these need be called with
 // eNoBlck for PMod (the physical representation need be changed, but the
 // corresponding Modification has to be ignored by the ThermalUnitBlock
 // since it has been self-inflicted), and with eDryRun for AMod (the
 // abstract representation has been changed already)
 Index th = f_time_horizon;

 // the Objective carries f_scale times the cost of one copy of the unit
 // [see UnitBlock::scale()], the physical representation the cost of one
 // copy, hence the coefficients are divided back before being stored; with
 // a zero scale the Objective tells nothing of the cost of a copy, which is
 // then kept as it is
 if( f_scale == 0 )
  return;
 auto coef = [ qf , this ]( Index i ) {
  return( qf->get_linear_coefficient( i ) / f_scale );
  };
 auto qcoef = [ qf , this ]( Index i ) {
  return( qf->get_quadratic_coefficient( i ) / f_scale );
  };

 // C05FunctionModLinRngd / DQuadFunctionModRngd - - - - - - - - - - - - - - -
 // split the Modification in up to 5 physical Modification by calling the
 // appropriate set_*() methods (ranged version) for those among startup,
 // power, commitment, primary/secondary reserve variables whose coefficient
 // change. This heavily relies on the fact that variables of the same type
 // are consecutive (and ordered in the obvious way) when set as coefficients
 // in the Objective. C05FunctionModLinRngd only reports changes to linear
 // coefficients, while DQuadFunctionModRngd also reports changes to quadratic
 // coefficients; since only the active power variables carry a non-zero
 // quadratic coefficient, set_quad_term() is called exclusively for that
 // section, and only when the Modification is a DQuadFunctionModRngd (the
 // with_quad flag).

 const Range * rng = nullptr;
 bool with_quad = false;
 if( const auto tmod = dynamic_cast< C05FunctionModLinRngd * >( mod ) )
  rng = & tmod->range();
 else
  if( const auto tmod = dynamic_cast< DQuadFunctionModRngd * >( mod ) ) {
   rng = & tmod->range();
   with_quad = true;
   }

 if( rng ) {
  Index l = rng->first;
  Index r = rng->second;

  if( r > qf->get_num_active_var() )
   throw( std::invalid_argument(
		  "ThermalUnitBlock::add_Modification: invalid Range [" +
		  std::to_string( l ) + ", " + std::to_string( r ) + ")" ) );

  // the Variable appended by a derived class [see objective_tail()] are the
  // very last ones: peel off a change to them, and route it to the derived
  // class; afterwards work with the count of the ThermalUnitBlock ones
  const Index ntail = objective_tail();
  const Index ntub = qf->get_num_active_var() - ntail;
  if( ntail && ( r > ntub ) ) {
   objective_tail_change( qf , std::max( l , ntub ) - ntub , r - ntub );
   r = ntub;
   if( r <= l )
    return;  // only the appended coefficients changed
   }

  // the design (investment) variable, if present, is the trailing single-var
  // block at index num_active_var - 1: peel off a change to it and route it to
  // update_objective_investment(), which (re)issues a eSetInvCost
  // ThermalUnitBlockMod for the DP Solvers (the abstract Objective has been
  // changed already, hence eDryRun). Afterwards work with the design-free count
  // nav, which makes the time-indexed sections below start at index 0.
  const Index has_design = ( f_InvestmentCost != 0 ) ? Index( 1 ) : Index( 0 );
  const Index nav = ntub - has_design;
  if( has_design && ( r > nav ) ) {
   update_objective_investment( par , eDryRun );
   r = nav;
   if( r <= l )
    return;  // only the design coefficient changed
   }

  // reactive power variables are the trailing section: peel off any index in
  // [ q_start , nav ) and route it to set_reactive_linear_term(),
  // then let the code below handle the remaining (non-reactive) prefix only.
  if( f_reactive_power ) {
   const Index q_start = nav - th;
   if( r > q_start ) {
    const Index rl = std::max( l , q_start );
    std::vector< double > qv( r - rl );
    auto qvit = qv.begin();
    for( Index i = rl ; i < r ; )
     *( qvit++ ) = coef( i++ );
    set_reactive_linear_term( qv.begin() ,
                              Range( rl - q_start , r - q_start ) ,
                              par , eDryRun );
    if( l >= q_start )
     return;  // the whole range was reactive
    r = q_start;
    }
   }

  std::vector< double > nv( r - l + 1 );
  std::vector< double > nvq;
  if( with_quad )
   nvq.resize( r - l + 1 );
  Index gl = 0;
  Index gr = th - init_t;

  if( l < gr ) {  // startup variables
   Index r2 = std::min( r , gr );
   auto nvit = nv.begin();
   for( Index i = l ; i < r2 ; )
    *( nvit++ ) = coef( i++ );
   // set_startup_costs( Range ) expects rng in time-space [ init_t ,
   // f_time_horizon ); active-var indices [ l , r2 ) correspond to time
   // indices [ l + init_t , r2 + init_t )
   set_startup_costs( nv.begin() ,
                      Range( l + init_t , r2 + init_t ) ,
                      par , eDryRun );
   l = r2;
   if( l == r )
    return;
   }

  if( f_shut_down_in_obj ) {  // shut-down variables
   gl = gr;
   gr += th - init_t;

   if( l < gr ) {
    Index r2 = std::min( r , gr );
    auto nvit = nv.begin();
    for( Index i = l ; i < r2 ; )
     *( nvit++ ) = coef( i++ );
    // as for the start-up ones, the active-var indices [ gl , gr ) map to
    // the time indices [ init_t , f_time_horizon )
    set_shutdown_costs( nv.begin() ,
                        Range( l - gl + init_t , r2 - gl + init_t ) ,
                        par , eDryRun );
    l = r2;
    if( l == r )
     return;
    }
   }

  gl = gr;
  gr += th;

  if( l < gr ) {  // active power variables
   Index r2 = std::min( r , gr );
   auto nvit = nv.begin();
   for( Index i = l ; i < r2 ; )
    *( nvit++ ) = coef( i++ );
   set_linear_term( nv.begin() , Range( l - gl , r2 - gl ) , par , eDryRun );
   // with perspective cuts the active power variables carry no quadratic
   // coefficient (the quadratic term lives on the cut variables, see below)
   if( with_quad && ! ( AR & PCuts ) ) {
    auto nvqit = nvq.begin();
    for( Index i = l ; i < r2 ; )
     *( nvqit++ ) = qcoef( i++ );
    set_quad_term( nvq.begin() , Range( l - gl , r2 - gl ) , par , eDryRun );
    }

   l = r2;
   if( l == r )
    return;
   }

  gl = gr;
  gr += th;

  if( l < gr ) {  // commitment variables
   Index r2 = std::min( r , gr );
   auto nvit = nv.begin();
   for( Index i = l ; i < r2 ; )
    *( nvit++ ) = coef( i++ );
   set_const_term( nv.begin() , Range( l - gl , r2 - gl ) , par , eDryRun );

   l = r2;
   if( l == r )
    return;
   }

  if( ! v_RefSchedule.empty() ) {  // schedule-deviation variables
   gl = gr;
   gr += th;

   if( l < gr ) {
    // their coefficient is the scale factor: a Modification whose range
    // crosses this section is fine as long as it leaves it there, which is
    // what a restore of the original costs does
    Index r2 = std::min( r , gr );
    for( Index i = l ; i < r2 ; ++i )
     if( ( qf->get_linear_coefficient( i ) != f_scale ) ||
	 ( with_quad && ( qf->get_quadratic_coefficient( i ) != 0 ) ) )
      throw( std::invalid_argument( "ThermalUnitBlock::add_Modification: the "
       "coefficients of the schedule-deviation variables cannot change" ) );

    l = r2;
    if( l == r )
     return;
    }
   }

  if( f_primary_in_obj ) {
   gl = gr;
   gr += th;

   if( l < gr ) {  // primary spinning reserve variables
    Index r2 = std::min( r , gr );
    auto nvit = nv.begin();
    for( Index i = l ; i < r2 ; )
     *( nvit++ ) = coef( i++ );
    set_primary_spinning_reserve_cost( nv.begin() ,
                                       Range( l - gl , r2 - gl ) ,
                                       par , eDryRun );
    l = r2;
    if( l == r )
     return;
    }
   }

  if( f_secondary_in_obj ) {
   gl = gr;
   gr += th;

   if( l < gr ) {  // secondary spinning reserve variables
    Index r2 = std::min( r , gr );
    auto nvit = nv.begin();
    for( Index i = l ; i < r2 ; )
     *( nvit++ ) = coef( i++ );
    set_secondary_spinning_reserve_cost( nv.begin() ,
                                         Range( l - gl , r2 - gl ) ,
                                         par , eDryRun );
    l = r2;
    if( l == r )
     return;
    }
   }

  if( AR & PCuts ) {  // perspective-cut variables
   // their linear coefficient is the quadratic term of the cost (see
   // generate_objective()); map the change back to set_quad_term()
   gl = gr;
   const auto form = AR & FormMsk;

   if( ( form == tbinForm ) || ( form == TForm ) ||
       ( form == ptForm ) || ( form == SUSDForm ) ) {
    // time-indexed cut variables
    gr += th;
    if( l < gr ) {
     Index r2 = std::min( r , gr );
     auto nvit = nv.begin();
     for( Index i = l ; i < r2 ; )
      *( nvit++ ) = coef( i++ );
     set_quad_term( nv.begin() , Range( l - gl , r2 - gl ) , par , eDryRun );
     if( r2 == r )
      return;
     }
    }
   else {
    // arc-indexed cut variables: map each arc back to its time instant
    auto arc_decode = [ & ]( const auto & Z ) {
     gr += Z.size();
     if( l >= gr )
      return( false );
     Index r2 = std::min( r , gr );
     std::map< Index , double > tv;
     for( Index i = l ; i < r2 ; ++i )
      tv[ Z[ i - gl ].first ] = coef( i );
     Subset nms( tv.size() );
     std::vector< double > tvv( tv.size() );
     auto nmsit = nms.begin();
     auto tvvit = tvv.begin();
     for( const auto & p : tv ) {
      *( nmsit++ ) = p.first;
      *( tvvit++ ) = p.second;
      }
     set_quad_term( tvv.begin() , std::move( nms ) , true , par , eDryRun );
     return( r2 == r );
     };

    bool done;
    if( form == DPForm )
     done = arc_decode( v_Z_h_k );
    else
     if( form == SUForm )
      done = arc_decode( v_Z_h );
     else
      done = arc_decode( v_Z_k );
    if( done )
     return;
    }
   }

  throw( std::invalid_argument( "ThermalUnitBlock::add_Modification: invalid "
				"variable in Modification" ) );
  return;

  }  // end( Rngd )

 // C05FunctionModLinSbst / DQuadFunctionModSbst - - - - - - - - - - - - - - -
 // split the Modification in up to 5 physical Modification by calling the
 // appropriate set_*() methods (subset version); see the comment on the
 // ranged case above: the same considerations apply, with the quadratic
 // coefficient being handled only for the active power section and only
 // when the Modification is a DQuadFunctionModSbst.

 const Subset * sbs = nullptr;
 with_quad = false;
 if( const auto tmod = dynamic_cast< C05FunctionModLinSbst * >( mod ) )
  sbs = & tmod->subset();
 else
  if( const auto tmod = dynamic_cast< DQuadFunctionModSbst * >( mod ) ) {
   sbs = & tmod->subset();
   with_quad = true;
   }

 Subset tail_reduced;  // storage when the appended indices are peeled off
 Subset des_reduced;  // storage when the trailing design index is peeled off
 Subset reduced;  // storage when the reactive tail has to be peeled off
 if( sbs ) {
  if( sbs->back() > qf->get_num_active_var() )
   throw( std::invalid_argument( "ThermalUnitBlock::add_Modification: "
				 "invalid Subset" ) );

  // the Variable appended by a derived class (see the ranged case)
  const Index ntail = objective_tail();
  const Index ntub = qf->get_num_active_var() - ntail;
  if( ntail && ( ! sbs->empty() ) && ( sbs->back() >= ntub ) ) {
   auto it = std::lower_bound( sbs->begin() , sbs->end() , ntub );
   objective_tail_change( qf , *it - ntub , sbs->back() + 1 - ntub );
   tail_reduced.assign( sbs->begin() , it );
   sbs = & tail_reduced;
   if( sbs->empty() )
    return;  // only the appended coefficients changed
   }

  // the design (investment) variable, if present, is the trailing single-var
  // block at index num_active_var - 1 (i.e., nav): peel a change to it off the
  // (sorted) Subset and route it to update_objective_investment() (see the
  // ranged case). Afterwards work with the design-free count nav.
  const Index has_design = ( f_InvestmentCost != 0 ) ? Index( 1 ) : Index( 0 );
  const Index nav = ntub - has_design;
  if( has_design && ( ! sbs->empty() ) && ( sbs->back() >= nav ) ) {
   update_objective_investment( par , eDryRun );
   des_reduced.assign( sbs->begin() , sbs->end() - 1 );
   sbs = & des_reduced;
   if( sbs->empty() )
    return;  // only the design coefficient changed
   }

  // peel off the trailing reactive section (see the ranged case): entries in
  // [ q_start , nav ) go to set_reactive_linear_term(); the
  // remaining prefix is handled by the section walk below unchanged.
  if( f_reactive_power ) {
   const Index q_start = nav - th;
   auto qit = std::lower_bound( sbs->begin() , sbs->end() , q_start );
   if( qit != sbs->end() ) {
    Subset qms( std::distance( qit , sbs->end() ) );
    std::vector< double > qv( qms.size() );
    auto qvit = qv.begin();
    auto qmsit = qms.begin();
    for( auto it = qit ; it != sbs->end() ; ++it ) {
     *( qvit++ ) = coef( *it );
     *( qmsit++ ) = *it - q_start;
     }
    set_reactive_linear_term( qv.begin() , std::move( qms ) , true ,
                              par , eDryRun );
    if( qit == sbs->begin() )
     return;  // the whole subset was reactive
    reduced.assign( sbs->begin() , qit );
    sbs = & reduced;
    }
   }

  std::vector< double > nv( sbs->size() );
  std::vector< double > nvq;
  if( with_quad )
   nvq.resize( sbs->size() );
  auto l = sbs->begin();
  Index gl = 0;
  Index gr = th - init_t;

  if( *l < gr ) {  // startup variables
   auto r = l;
   for( ++r ; ( r != sbs->end() ) && ( *r < gr ) ; )
    ++r;
   Subset nms( std::distance( l , r ) );
   auto nvit = nv.begin();
   auto nmsit = nms.begin();
   // set_startup_costs( Subset ) expects subset entries in time-space
   // [ init_t , f_time_horizon ); active-var indices [ 0 , th - init_t )
   // map to time indices by adding init_t
   while( l != r ) {
    *( nvit++ ) = coef( *l );
    *( nmsit++ ) = *( l++ ) + init_t;
    }
   set_startup_costs( nv.begin() , std::move( nms ) , true , par , eDryRun );
   if( r == sbs->end() )
    return;
   }

  if( f_shut_down_in_obj ) {  // shut-down variables
   gl = gr;
   gr += th - init_t;

   if( *l < gr ) {
    auto r = l;
    for( ++r ; ( r != sbs->end() ) && ( *r < gr ) ; )
     ++r;
    Subset nms( std::distance( l , r ) );
    auto nvit = nv.begin();
    auto nmsit = nms.begin();
    while( l != r ) {
     *( nvit++ ) = coef( *l );
     *( nmsit++ ) = *( l++ ) - gl + init_t;
     }
    set_shutdown_costs( nv.begin() , std::move( nms ) , true , par ,
                        eDryRun );
    if( r == sbs->end() )
     return;
    }
   }

  gl = gr;
  gr += th;

  if( *l < gr ) {  // active power variables
   auto r = l;
   for( ++r ; ( r != sbs->end() ) && ( *r < gr ) ; )
    ++r;
   Subset nms( std::distance( l , r ) );
   auto nvit = nv.begin();
   auto nmsit = nms.begin();
   auto nvqit = nvq.begin();
   while( l != r ) {
    *( nvit++ ) = coef( *l );
    if( with_quad )
     *( nvqit++ ) = qcoef( *l );
    *( nmsit++ ) = *( l++ ) - gl;
    }
   // with perspective cuts the active power variables carry no quadratic
   // coefficient (the quadratic term lives on the cut variables, see below)
   if( with_quad && ! ( AR & PCuts ) ) {
    Subset nmsq = nms;  // copy, as it is moved
    set_quad_term( nvq.begin() , std::move( nmsq ) , true , par , eDryRun );
    }
   set_linear_term( nv.begin() , std::move( nms ) , true , par , eDryRun );
   if( r == sbs->end() )
    return;
   }

  gl = gr;
  gr += th;

  if( *l < gr ) {  // commitment variables
   auto r = l;
   for( ++r ; ( r != sbs->end() ) && ( *r < gr ) ; )
    ++r;
   Subset nms( std::distance( l , r ) );
   auto nvit = nv.begin();
   auto nmsit = nms.begin();
   while( l != r ) {
    *( nvit++ ) = coef( *l );
    *( nmsit++ ) = *( l++ ) - gl;
    }
   set_const_term( nv.begin() , std::move( nms ) , true , par , eDryRun );
   if( r == sbs->end() )
    return;
   }

  if( ! v_RefSchedule.empty() ) {  // schedule-deviation variables
   gl = gr;
   gr += th;

   if( *l < gr ) {
    // their coefficient is the scale factor [see the ranged case]
    auto r = l;
    for( ++r ; ( r != sbs->end() ) && ( *r < gr ) ; )
     ++r;
    while( l != r ) {
     if( ( qf->get_linear_coefficient( *l ) != f_scale ) ||
	 ( with_quad && ( qf->get_quadratic_coefficient( *l ) != 0 ) ) )
      throw( std::invalid_argument( "ThermalUnitBlock::add_Modification: the "
       "coefficients of the schedule-deviation variables cannot change" ) );
     ++l;
     }
    if( r == sbs->end() )
     return;
    }
   }

  if( f_primary_in_obj ) {
   gl = gr;
   gr += th;

   if( *l < gr ) {  // primary spinning reserve variables
    auto r = l;
    for( ++r ; ( r != sbs->end() ) && ( *r < gr ) ; )
     ++r;
    Subset nms( std::distance( l , r ) );
    auto nvit = nv.begin();
    auto nmsit = nms.begin();
    while( l != r ) {
     *( nvit++ ) = coef( *l );
     *( nmsit++ ) = *( l++ ) - gl;
     }
    set_primary_spinning_reserve_cost( nv.begin() , std::move( nms ) ,
                                       true , par , eDryRun );
    if( r == sbs->end() )
     return;
    }
   }

  if( f_secondary_in_obj ) {
   gl = gr;
   gr += th;

   if( *l < gr ) {  // secondary spinning reserve variables
    auto r = l;
    for( ++r ; ( r != sbs->end() ) && ( *r < gr ) ; )
     ++r;
    Subset nms( std::distance( l , r ) );
    auto nvit = nv.begin();
    auto nmsit = nms.begin();
    while( l != r ) {
     *( nvit++ ) = coef( *l );
     *( nmsit++ ) = *( l++ ) - gl;
     }
    set_secondary_spinning_reserve_cost( nv.begin() , std::move( nms ) ,
                                         true , par , eDryRun );
    if( r == sbs->end() )
     return;
    }
   }

  if( AR & PCuts ) {  // perspective-cut variables
   // their linear coefficient is the quadratic term of the cost (see
   // generate_objective()); map the change back to set_quad_term()
   gl = gr;
   const auto form = AR & FormMsk;

   if( ( form == tbinForm ) || ( form == TForm ) ||
       ( form == ptForm ) || ( form == SUSDForm ) ) {
    // time-indexed cut variables
    gr += th;
    if( *l < gr ) {
     auto r = l;
     for( ++r ; ( r != sbs->end() ) && ( *r < gr ) ; )
      ++r;
     Subset nms( std::distance( l , r ) );
     auto nvit = nv.begin();
     auto nmsit = nms.begin();
     while( l != r ) {
      *( nvit++ ) = coef( *l );
      *( nmsit++ ) = *( l++ ) - gl;
      }
     set_quad_term( nv.begin() , std::move( nms ) , true , par , eDryRun );
     }
    }
   else {
    // arc-indexed cut variables: map each arc back to its time instant
    auto arc_decode = [ & ]( const auto & Z ) {
     gr += Z.size();
     if( *l >= gr )
      return;
     std::map< Index , double > tv;
     while( ( l != sbs->end() ) && ( *l < gr ) ) {
      tv[ Z[ *l - gl ].first ] = coef( *l );
      ++l;
      }
     Subset nms( tv.size() );
     std::vector< double > tvv( tv.size() );
     auto nmsit = nms.begin();
     auto tvvit = tvv.begin();
     for( const auto & p : tv ) {
      *( nmsit++ ) = p.first;
      *( tvvit++ ) = p.second;
      }
     set_quad_term( tvv.begin() , std::move( nms ) , true , par , eDryRun );
     };

    if( form == DPForm )
     arc_decode( v_Z_h_k );
    else
     if( form == SUForm )
      arc_decode( v_Z_h );
     else
      arc_decode( v_Z_k );
    }
   }

  if( l != sbs->end() )
   throw( std::invalid_argument( "ThermalUnitBlock::add_Modification: invalid "
				 "variable in Modification" ) );
  return;

  }  // end( Sbst )

 throw( std::invalid_argument( "ThermalUnitBlock:: unsupported FunctionMod "
			       "from Objective" ) );

 }  // end( ThermalUnitBlock::handle_objective_change )

/*--------------------------------------------------------------------------*/
/*------------------ METHODS OF ThermalUnitBlockSolution -------------------*/
/*--------------------------------------------------------------------------*/

void ThermalUnitBlockSolution::deserialize( const netCDF::NcGroup & group )
{
 // call the method of the base class
 UnitBlockSolution::deserialize( group );

 if( f_number_generators != 1 )
  throw( std::logic_error( "ThermalUnitBlockSolution::deserialize: "
         "thermals have only one generator" ) );

 // deserialize the design - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! ::deserialize< double >( group , f_design , "ThermalDesign" ) )
  f_design = dNaN;

 // deserialize the start-up and shut-down indicators, if there - - - - - - -
 // a Solution written before they were saved simply has none, and the
 // indicators are derived from the commitment as they used to be
 v_start_up.clear();
 v_shut_down.clear();
 ::deserialize( group , v_start_up , "ThermalStartUp" , "NumberStartUp" );
 ::deserialize( group , v_shut_down , "ThermalShutDown" , "NumberStartUp" );

 }  // end( ThermalUnitBlockSolution::deserialize )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlockSolution::read( const Block * block )
{
 auto TUB = dynamic_cast< const ThermalUnitBlock * >( block );
 if( ! TUB )
  throw( std::invalid_argument( "ThermalUnitBlockSolution::read: block "
        "is not a ThermalUnitBlock" ) );

 UnitBlockSolution::read( TUB );  // call the method of the base class

 // read the design- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 f_design = TUB->get_const_design().get_value();

 // read the start-up and shut-down indicators - - - - - - - - - - - - - - -
 // they are saved rather than left to be derived from the commitment when
 // the Solution is written back, because deriving them only works for one
 // schedule [see get_start_up() in ThermalUnitBlockSolution]
 const auto nsu = TUB->get_number_start_up();
 if( nsu ) {
  auto tub = const_cast< ThermalUnitBlock * >( TUB );
  v_start_up.resize( nsu );
  auto su = tub->get_start_up();
  for( Index i = 0 ; i < nsu ; ++i )
   v_start_up[ i ] = su[ i ].get_value();
  v_shut_down.resize( nsu );
  if( auto sd = tub->get_shut_down() )
   for( Index i = 0 ; i < nsu ; ++i )
    v_shut_down[ i ] = sd[ i ].get_value();
  else
   v_shut_down.clear();
  }
 else {
  v_start_up.clear();
  v_shut_down.clear();
  }

 }  // end( ThermalUnitBlockSolution::read )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlockSolution::write( Block * block )
{
 UnitBlockSolution::write( block );  // call the method of the base class

 auto TUB = dynamic_cast< ThermalUnitBlock * >( block );
 if( ! TUB )
  throw( std::invalid_argument( "ThermalUnitBlockSolution::write: block "
        "is not a ThermalUnitBlock" ) );

 // write the design - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 TUB->get_design().set_value( f_design );

 // (p, u) have just been restored into the Block by UnitBlockSolution::
 // write(); delegate all the formulation-specific bookkeeping (start_up /
 // shut_down, perspective-cut auxiliaries, ...) to the Block itself. This
 // shares the implementation with the inner DP Solvers and ensures that
 // LagBFunction::get_linearization_constant(), which reads variable
 // values to recompute f(x*) via its saved CostMatrix, sees a state
 // consistent with the saved (p, u).
 TUB->set_solution();

 // ... unless the indicators were saved, in which case they are the ones to
 // use: in a convex combination of Solution they have been averaged like
 // everything else, while deriving them from the averaged commitment gives
 // the start-ups of the average, which are fewer and therefore cheaper than
 // the average of the start-ups
 if( ! v_start_up.empty() ) {
  const auto nsu = std::min( Index( v_start_up.size() ) ,
                             TUB->get_number_start_up() );
  if( auto su = TUB->get_start_up() )
   for( Index i = 0 ; i < nsu ; ++i )
    su[ i ].set_value( v_start_up[ i ] );
  if( ! v_shut_down.empty() )
   if( auto sd = TUB->get_shut_down() )
    for( Index i = 0 ; i < std::min( Index( v_shut_down.size() ) , nsu ) ;
         ++i )
     sd[ i ].set_value( v_shut_down[ i ] );
  }

 }  // end( ThermalUnitBlockSolution::write )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlockSolution::serialize( netCDF::NcGroup & group ) const
{
 UnitBlockSolution::serialize( group );  // call the method of the base class

 // serialize the design- - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! std::isnan( f_design ) )
  ::serialize< double >( group , "ThermalDesign" , netCDF::NcDouble() ,
       f_design );

 // serialize the start-up and shut-down indicators, if any - - - - - - - - -
 // they are as many as the instants in which the unit can change state,
 // which is its own dimension [see ThermalUnitBlock::get_number_start_up()]
 if( ! v_start_up.empty() ) {
  auto nsu = group.addDim( "NumberStartUp" , v_start_up.size() );
  group.addVar( "ThermalStartUp" , netCDF::NcDouble() , nsu ).putVar(
                    { 0 } , { v_start_up.size() } , v_start_up.data() );
  if( ! v_shut_down.empty() )
   group.addVar( "ThermalShutDown" , netCDF::NcDouble() , nsu ).putVar(
                    { 0 } , { v_shut_down.size() } , v_shut_down.data() );
  }

 }  // end( ThermalUnitBlockSolution::serialize )

/*--------------------------------------------------------------------------*/

ThermalUnitBlockSolution * ThermalUnitBlockSolution::scale( double factor )
 const
{
 auto sol = clone();

 if( factor == 1 )
  return( sol );

 guts_of_scale( sol , factor );

 if( ! std::isnan( f_design ) )
  sol->f_design *= factor;

 for( auto & v : sol->v_start_up )
  v *= factor;
 for( auto & v : sol->v_shut_down )
  v *= factor;

 return( sol );

 }  // end( ThermalUnitBlockSolution::scale )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlockSolution::sum( const Solution * solution ,
           double multiplier )
{
 // call the method of the base class
 UnitBlockSolution::sum( solution , multiplier );

 auto TUBS = dynamic_cast< const ThermalUnitBlockSolution * >( solution );
 if( ! TUBS )
  throw( std::invalid_argument( "ThermalUnitBlockSolution::sum: "
        "solution not a ThermalUnitBlockSolution" ) );

 if( ! std::isnan( f_design ) )
  f_design += TUBS->f_design * multiplier;

 if( v_start_up.size() != TUBS->v_start_up.size() )
  throw( std::invalid_argument( "ThermalUnitBlockSolution::sum: "
        "inconsistent start-up indicators, " +
        std::to_string( v_start_up.size() ) + " against " +
        std::to_string( TUBS->v_start_up.size() ) ) );
 for( Index i = 0 ; i < Index( v_start_up.size() ) ; ++i )
  v_start_up[ i ] += TUBS->v_start_up[ i ] * multiplier;

 if( v_shut_down.size() != TUBS->v_shut_down.size() )
  throw( std::invalid_argument( "ThermalUnitBlockSolution::sum: "
        "inconsistent shut-down indicators, " +
        std::to_string( v_shut_down.size() ) + " against " +
        std::to_string( TUBS->v_shut_down.size() ) ) );
 for( Index i = 0 ; i < Index( v_shut_down.size() ) ; ++i )
  v_shut_down[ i ] += TUBS->v_shut_down[ i ] * multiplier;

 }  // end( ThermalUnitBlockSolution::sum )

/*--------------------------------------------------------------------------*/

ThermalUnitBlockSolution * ThermalUnitBlockSolution::clone( bool empty ) const
{
 auto * sol = new ThermalUnitBlockSolution();

 if( ! empty ) {
  guts_of_clone( sol );
  sol->f_design = f_design;
  sol->v_start_up = v_start_up;
  sol->v_shut_down = v_shut_down;
  }

 return( sol );

 }  // end( ThermalUnitBlockSolution::clone )

/*--------------------------------------------------------------------------*/
/*------------------- End File ThermalUnitBlock.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
