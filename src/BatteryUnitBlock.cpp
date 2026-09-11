/*--------------------------------------------------------------------------*/
/*----------------------- File BatteryUnitBlock.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the BatteryUnitBlock class.
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
 * \copyright &copy; by Antonio Frangioni, Ali Ghezelsoflu,
 *                      Rafael Durbano Lobato, Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <iostream>

#include <random>

#include <map>

#include "BatteryUnitBlock.h"

#include "LinearFunction.h"

#include "FRealObjective.h"

#include "UnitBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register BatteryUnitBlock to the Block factory
SMSpp_insert_in_factory_cpp_1( BatteryUnitBlock );

// register BatteryUnitBlockSolution to the Solution factory
SMSpp_insert_in_factory_cpp_0( BatteryUnitBlockSolution );

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS OF BatteryUnitBlock ----------------------*/
/*--------------------------------------------------------------------------*/

BatteryUnitBlock::~BatteryUnitBlock()
{
 Constraint::clear( active_power_bounds_Const );
 Constraint::clear( active_power_bounds_design_Const );
 Constraint::clear( intake_outtake_upper_bounds_design_Const );
 Constraint::clear( storage_level_bounds_design_Const );
 Constraint::clear( intake_outtake_binary_Const );
 Constraint::clear( power_intake_outtake_Const );
 Constraint::clear( ramp_up_Const );
 Constraint::clear( ramp_down_Const );
 Constraint::clear( demand_Const );
 Constraint::clear( Reference_Schedule_Const );
 Constraint::clear( ReactivePower_Bound_Const );
 Constraint::clear( intake_outtake_bounds_Const );
 Constraint::clear( primary_upper_bound_Const );
 Constraint::clear( secondary_upper_bound_Const );
 Constraint::clear( storage_level_bounds_Const );

 Constraint::clear( battery_binary_bound_Const );

 batt_design_bound_Const.clear();
 conv_design_bound_Const.clear();

 objective.clear();
}

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::deserialize( const netCDF::NcGroup & group )
{
 // Deserialize data from the base class
 UnitBlock::deserialize( group );

 // Mandatory variables

 ::deserialize( group , "MinStorage" , f_time_horizon , v_MinStorage ,
                false , true , v_change_intervals );

 ::deserialize( group , "MaxStorage" , f_time_horizon , v_MaxStorage ,
                false , true , v_change_intervals );

 ::deserialize( group , "MaxPower" , f_time_horizon , v_MaxPower ,
                false , true , v_change_intervals );

 // Optional variables

 if( ! ::deserialize( group , "MinPower" , f_time_horizon , v_MinPower ,
                      true , true , v_change_intervals ) ) {
  v_MinPower.resize( v_MaxPower.size() );
  std::copy( v_MaxPower.begin() , v_MaxPower.end() , v_MinPower.begin() );
  std::transform( v_MinPower.cbegin() , v_MinPower.cend() , v_MinPower.begin() ,
                  []( double p ) { return( -p ); } );
 }

 if( ! ::deserialize( group , "ConverterMaxPower" , f_time_horizon ,
                      v_ConvMaxPower , true , true , v_change_intervals ) ) {
  v_ConvMaxPower.resize( v_MaxPower.size() );
  std::copy( v_MaxPower.begin() , v_MaxPower.end() , v_ConvMaxPower.begin() );
 }

 ::deserialize( group , f_MaxCRateCharge , "MaxCRateCharge" );
 ::deserialize( group , f_MaxCRateDischarge , "MaxCRateDischarge" );

 ::deserialize( group , f_InitialStorage , "InitialStorage" );

 ::deserialize( group , f_InitialPower , "InitialPower" );

 ::deserialize( group , f_kappa , "Kappa" );

 ::deserialize( group , "MaxPrimaryPower" , f_time_horizon ,
                v_MaxPrimaryPower , true , true , v_change_intervals );
 ::deserialize( group , "MaxSecondaryPower" , f_time_horizon ,
                v_MaxSecondaryPower , true , true , v_change_intervals );

 ::deserialize( group , "DeltaRampUp" , f_time_horizon , v_DeltaRampUp ,
                true , true , v_change_intervals );
 ::deserialize( group , "DeltaRampDown" , f_time_horizon , v_DeltaRampDown ,
                true , true , v_change_intervals );

 ::deserialize( group , "Demand" , f_time_horizon , v_Demand ,
                true , true , v_change_intervals );

 ::deserialize( group , "StoringBatteryRho" , f_time_horizon ,
                v_StoringBatteryRho , true , true , v_change_intervals );
 ::deserialize( group , "ExtractingBatteryRho" , f_time_horizon ,
                v_ExtractingBatteryRho , true , true , v_change_intervals );
 ::deserialize( group , "StandingBatteryRho" , f_time_horizon ,
                v_StandingBatteryRho , true , true , v_change_intervals );

 if( ! ::deserialize( group , "Cost" , f_time_horizon , v_Cost ,
                      true , true , v_change_intervals ) )
  v_Cost.resize( f_time_horizon );

 if( ::deserialize( group , f_BattInvestmentCost , "BatteryInvestmentCost" ) ) {

  ::deserialize( group , f_BattMinCapacityDesign , "BatteryMinCapacityDesign" );

  ::deserialize( group , f_BattMaxCapacityDesign , "BatteryMaxCapacityDesign" );
 }

 if( ::deserialize( group , f_ConvInvestmentCost , "ConverterInvestmentCost" ) ) {

  ::deserialize( group , f_ConvMinCapacityDesign , "ConverterMinCapacityDesign" );

  ::deserialize( group , f_ConvMaxCapacityDesign , "ConverterMaxCapacityDesign" );
 }

 ::deserialize( group , f_scale , "Scale" );

 // variables for AC elements
 if( ::deserialize( group , "MaxReactivePower" , f_time_horizon ,
		    v_MaxReactivePower , true , true , v_change_intervals )
     )
  if( std::all_of( v_MaxReactivePower.begin() , v_MaxReactivePower.end() ,
		   []( double i ) { return( i == 0 ); } ) )
   v_MaxReactivePower.clear();

 if( ::deserialize( group , "MinReactivePower" , f_time_horizon ,
		    v_MinReactivePower , true , true , v_change_intervals )
     )
  if( std::all_of( v_MinReactivePower.begin() , v_MinReactivePower.end() ,
		   []( double i ) { return( i == 0 ); } ) )
   v_MinReactivePower.clear();

 check_data_consistency();  // final checks

 }  // end( BatteryUnitBlock::deserialize )

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG

std::vector< std::string > BatteryUnitBlock::expected_dims( void ) const {
 auto ret = UnitBlock::expected_dims();
 ret.push_back( "NumberIntervals" );

 return( ret );
 }

/*--------------------------------------------------------------------------*/

std::vector< std::string > BatteryUnitBlock::expected_vars( void ) const {
 static const std::vector< std::string > ev =
 { "MinStorage" , "MaxStorage" , "InitialStorage" , "MinPower" , "MaxPower" ,
   "InitialPower" , "ConverterMaxPower" , "MaxPrimaryPower" ,
   "MaxSecondaryPower" , "DeltaRampUp" , "DeltaRampDown" ,
   "StoringBatteryRho" , "ExtractingBatteryRho" , "StandingBatteryRho",
   "Cost" , "Demand" ,
   "Kappa" , "MaxCRateCharge" , "MaxCRateDischarge" , "BatteryMaxCapacity" ,
   "ConverterMaxCapacity" , "BatteryInvestmentCost" ,
   "ConverterInvestmentCost" , "BatteryMinCapacityDesign" ,
   "ConverterMinCapacityDesign" , "BatteryMaxCapacityDesign" ,
   "ConverterMaxCapacityDesign" , "MinReactivePower", "MaxReactivePower",
   "ReferenceSchedule" , "Scale" };

 auto ret = UnitBlock::expected_vars();
 ret.insert( ret.end() , ev.begin() , ev.end() );

 return( ret );
 }

#endif

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::check_data_consistency( void ) const
{
 // InvestmentCost

 if( ( ( f_BattInvestmentCost != 0 ) || ( f_ConvInvestmentCost != 0 ) ) &&
     ( f_InitialStorage > 0 ) )
  throw( std::logic_error( "BatteryUnitBlock::check_data_consistency: the "
                           "presence of the investment cost of the battery "
                           "allows the model to switch into the strategic "
                           "scenario mode, but the presence of also a positive "
                           "initial storage, typical of the operative "
                           "scenario, is incompatible." ) );

 // Min/Max battery capacity design

 if( f_BattMinCapacityDesign < 0 )
  throw( std::logic_error( "BatteryMinCapacityDesign must be >= 0." ) );

 if( f_BattMaxCapacityDesign > 0 ) {

  if( f_BattMinCapacityDesign > f_BattMaxCapacityDesign )
   throw( std::logic_error(
    "BatteryMinCapacityDesign must be <= BatteryMaxCapacityDesign." ) );

  if( ( std::abs( f_BattMaxCapacityDesign ) == 1 ) && ( f_BattMinCapacityDesign > 1 ) )
   throw( std::logic_error(
    "BatteryMinCapacityDesign must be <= 1 when |BatteryMaxCapacityDesign| = 1." ) );
 }

 if( f_BattMaxCapacityDesign < 0 ) {
  if( f_BattMinCapacityDesign > 1 )
   throw( std::logic_error( "BatteryMinCapacityDesign must be <= 1 when "
                            "BatteryMaxCapacityDesign < 0 (binary design)." ) );
 }

 // Scale and design granularity are two equivalent ways to model multiple
 // identical modules; combining them with both > 1 over-counts the fleet.
 if( ( f_scale != 1 ) && ( std::abs( f_BattMaxCapacityDesign ) > 1 ) )
  throw( std::logic_error( "BatteryUnitBlock::check_data_consistency: "
                           "cannot combine Scale != 1 with "
                           "|BatteryMaxCapacityDesign| > 1; they represent "
                           "the same multi-module fleet and would double-count." ) );

 // Min/Max converter capacity design

 if( f_ConvMinCapacityDesign < 0 )
  throw( std::logic_error( "ConverterMinCapacityDesign must be >= 0." ) );

 if( f_ConvMaxCapacityDesign > 0 ) {

  if( f_ConvMinCapacityDesign > f_ConvMaxCapacityDesign )
   throw( std::logic_error(
   "ConverterMinCapacityDesign must be <= ConverterMaxCapacityDesign." ) );

  if( ( std::abs( f_ConvMaxCapacityDesign ) == 1 ) && ( f_ConvMinCapacityDesign > 1 ) )
   throw( std::logic_error( "ConverterMinCapacityDesign must be <= 1 when "
                            "|ConverterMaxCapacityDesign| = 1." ) );
 }

 if( f_ConvMaxCapacityDesign < 0 ) {
  if( f_ConvMinCapacityDesign > 1 )
   throw( std::logic_error( "ConverterMinCapacityDesign must be <= 1 when "
                            "ConverterMaxCapacityDesign < 0 (binary design)." ) );
 }

 // Same Scale-vs-design exclusion as the battery (see note above).
 if( ( f_scale != 1 ) && ( std::abs( f_ConvMaxCapacityDesign ) > 1 ) )
  throw( std::logic_error( "BatteryUnitBlock::check_data_consistency: "
                           "cannot combine Scale != 1 with "
                           "|ConverterMaxCapacityDesign| > 1; they represent "
                           "the same multi-module fleet and would double-count." ) );

 // Minimum and maximum power

 assert( v_MinPower.size() == f_time_horizon );
 assert( v_MaxPower.size() == f_time_horizon );

 for( Index t = 0 ; t < f_time_horizon ; ++t )
  if( v_MinPower[ t ] > v_MaxPower[ t ] )
   throw( std::logic_error( "BatteryUnitBlock::check_data_consistency: minimum "
                            "power for time " + std::to_string( t ) + " is " +
                            std::to_string( v_MinPower[ t ] ) + ", which is "
                            "greater than the maximum power, which is " +
                            std::to_string( v_MaxPower[ t ] ) + "." ) );

 // Minimum and maximum storage levels

 assert( v_MinStorage.size() == f_time_horizon );
 assert( v_MaxStorage.size() == f_time_horizon );

 for( Index t = 0 ; t < f_time_horizon ; ++t )
  if( ( v_MinStorage[ t ] > v_MaxStorage[ t ] ) )
   throw( std::logic_error( "BatteryUnitBlock::check_data_consistency: maximum "
                            "and minimum storage levels must be such that "
                            "maximum_storage >= minimum_storage." ) );

 // Inefficiency of storing, extracting and standing energy

 if( ! v_StoringBatteryRho.empty() ) {
  assert( v_StoringBatteryRho.size() == f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   if( v_StoringBatteryRho[ t ] < 0. )
    throw( std::logic_error( "BatteryUnitBlock::check_data_consistency: invalid "
                             "efficiency of storing energy for time step " +
                             std::to_string( t ) + ": " +
                             std::to_string( v_StoringBatteryRho[ t ] ) +
                             ". It must not be lower than 0." ) );
 }

 if( ! v_ExtractingBatteryRho.empty() ) {
  assert( v_ExtractingBatteryRho.size() == f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   if( v_ExtractingBatteryRho[ t ] < 0. )
    throw( std::logic_error( "BatteryUnitBlock::check_data_consistency: invalid "
                             "inefficiency of extracting energy for time step " +
                             std::to_string( t ) + ": " +
                             std::to_string( v_ExtractingBatteryRho[ t ] ) +
                             ". It must not be lower than 0." ) );
 }

 // Roundtrip efficiency: the energy retrieved over a charge-discharge cycle
 // must not exceed the energy put in, i.e. StoringRho / ExtractingRho <= 1.
 // It is checked as StoringRho <= ExtractingRho to avoid a division (and to
 // correctly reject the degenerate ExtractingRho == 0 with StoringRho > 0). A
 // missing vector defaults to 1, consistently with how these coefficients are
 // used when building the storage balance constraints. Note that, when a time
 // step represents several elementary periods, an individual rho may exceed 1:
 // only the roundtrip ratio is constrained.

 if( ( ! v_StoringBatteryRho.empty() ) || ( ! v_ExtractingBatteryRho.empty() ) )
  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   const double storing = v_StoringBatteryRho.empty()
                          ? 1. : v_StoringBatteryRho[ t ];
   const double extracting = v_ExtractingBatteryRho.empty()
                             ? 1. : v_ExtractingBatteryRho[ t ];
   if( storing > extracting )
    throw( std::logic_error( "BatteryUnitBlock::check_data_consistency: invalid "
                             "roundtrip efficiency of battery for time step " +
                             std::to_string( t ) + ": StoringBatteryRho (" +
                             std::to_string( storing ) +
                             ") must not be greater than ExtractingBatteryRho ("
                             + std::to_string( extracting ) + ")." ) );
  }

 if( ! v_StandingBatteryRho.empty() ) {
  assert( v_StandingBatteryRho.size() == f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   if( v_StandingBatteryRho[ t ] > 1. ) 
    throw( std::logic_error( "BatteryUnitBlock::check_data_consistency: invalid"
                             " inefficiency of standing energy for time "
                             "step " + std::to_string( t ) + ": " +
                             std::to_string( v_StandingBatteryRho[ t ] ) +
                             ". It must not be greater than 1." ) );
   else if( v_StandingBatteryRho[ t ] < 0. )
    throw( std::logic_error( "BatteryUnitBlock::check_data_consistency: invalid"
                             " inefficiency of standing energy for time "
                             "step " + std::to_string( t ) + ": " +
                             std::to_string( v_StandingBatteryRho[ t ] ) +
                             ". It must not be lower than 0." ) );
 }

 // Delta ramp-up

 if( ! v_DeltaRampUp.empty() ) {
  assert( v_DeltaRampUp.size() == f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   if( v_DeltaRampUp[ t ] < 0 )
    throw( std::invalid_argument( "BatteryUnitBlock::check_data_consistency: "
                                  "wrong DeltaRampUp for time step " +
                                  std::to_string( t ) + ": " +
                                  std::to_string( v_DeltaRampUp[ t ] ) ) );
 }

 // Delta ramp-down

 if( ! v_DeltaRampDown.empty() ) {
  assert( v_DeltaRampDown.size() == f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   if( v_DeltaRampDown[ t ] < 0 )
    throw( std::invalid_argument( "BatteryUnitBlock::check_data_consistency: "
                                  "wrong DeltaRampDown for time step " +
                                  std::to_string( t ) + ": " +
                                  std::to_string( v_DeltaRampDown[ t ] ) ) );
 }

 // Maximum active power that can be used as primary reserve

 if( ! v_MaxPrimaryPower.empty() ) {
  assert( v_MaxPrimaryPower.size() == f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   if( v_MaxPrimaryPower[ t ] < 0 )
    throw( std::invalid_argument( "BatteryUnitBlock::check_data_consistency: "
                                  "the maximum power that can be used as "
                                  "primary reserve for time " +
                                  std::to_string( t ) + " is " +
                                  std::to_string( v_MaxPrimaryPower[ t ] ) +
                                  ", but it must be nonnegative." ) );
 }

 // Maximum active power that can be used as secondary reserve

 if( ! v_MaxSecondaryPower.empty() ) {
  assert( v_MaxSecondaryPower.size() == f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   if( v_MaxSecondaryPower[ t ] < 0 )
    throw( std::invalid_argument(
     "BatteryUnitBlock::check_data_consistency: the maximum power that "
     "can be used as secondary reserve for time " + std::to_string( t ) +
     " is " + std::to_string( v_MaxSecondaryPower[ t ] ) +
     ", but it must be nonnegative." ) );
 }

 // Demand

 if( ! v_Demand.empty() ) {
  assert( v_Demand.size() == f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   if( v_Demand[ t ] < 0 )
    throw( std::invalid_argument( "BatteryUnitBlock::check_data_consistency: "
                                  "demand for time " + std::to_string( t ) +
                                  " is " + std::to_string( v_Demand[ t ] ) +
                                  ", but it must be nonnegative." ) );
 }

 // Kappa

 if( f_kappa < 0 )
  throw( std::logic_error( "BatteryUnitBlock::check_data_consistency: "
                           "kappa must be nonnegative, but it is " +
                           std::to_string( f_kappa ) + "." ) );

}  // end( BatteryUnitBlock::check_data_consistency )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::generate_abstract_variables( Configuration * stvv )
{
 if( variables_generated() )  // variables have already been generated
  return;                     // nothing to do

 UnitBlock::generate_abstract_variables( stvv );

 // check if negative prices may occur- - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 bool negative_prices = false;
 if( ( ! stvv ) && f_BlockConfig )
  stvv = f_BlockConfig->f_static_variables_Configuration;
 if( auto sci = dynamic_cast< SimpleConfiguration< int > * >( stvv ) )
  negative_prices = sci->f_value;

 // Binary variables must be generated if negative prices may occur and if
 // there is some t with a roundtrip loss, i.e. StoringBatteryRho[ t ] <
 // ExtractingBatteryRho[ t ] (a missing vector defaults to 1); otherwise the
 // battery could charge and discharge at the same time to exploit the
 // efficiencies.
 //
 // The test is the direct comparison storing < extracting rather than the
 // ratio storing / extracting < 1: the two are equivalent (for extracting > 0),
 // but the direct form needs no division and gracefully handles the degenerate
 // extracting == 0. It is also scale-invariant w.r.t. the temporal resolution:
 // when a step aggregates several elementary periods both coefficients carry
 // the same period-duration factor (storing_rho = eta_c * tau, extracting_rho =
 // tau / eta_d), so an individual rho may exceed 1 while the ratio that drives
 // the loss, storing / extracting = eta_c * eta_d, stays independent of tau and
 // the inequality keeps detecting the loss correctly.

 bool generate_binary_variables = false;

 if( negative_prices && ( ( ! v_StoringBatteryRho.empty() ) ||
                          ( ! v_ExtractingBatteryRho.empty() ) ) ) {
  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   const double storing = v_StoringBatteryRho.empty()
                          ? 1. : v_StoringBatteryRho[ t ];
   const double extracting = v_ExtractingBatteryRho.empty()
                             ? 1. : v_ExtractingBatteryRho[ t ];
   if( storing < extracting ) {
    generate_binary_variables = true;
    break;
   }
  }
 }

 static constexpr double dNaN = std::numeric_limits< double >::quiet_NaN();

 // battery Design Variable - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( f_BattInvestmentCost != 0 ) {
  if( f_BattMaxCapacityDesign < 0 ) {
   // integer design in {0, ..., |BatteryMaxCapacityDesign|}; reduces to
   // binary {0,1} when |BatteryMaxCapacityDesign| == 1.
   if( std::abs( f_BattMaxCapacityDesign ) == 1.0 )
    batt_design.set_type( ColVariable::kBinary );
   else
    batt_design.set_type( ColVariable::kInteger );
  }
  else
   batt_design.set_type( ColVariable::kNonNegative );
  add_static_variable( batt_design , "x_battery" );
  }
 else
  batt_design.set_value( dNaN );

 // converter Design Variable - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( f_ConvInvestmentCost != 0 ) {
  if( f_ConvMaxCapacityDesign < 0 ) {
   // integer design in {0, ..., |ConverterMaxCapacityDesign|}; reduces to
   // binary {0,1} when |ConverterMaxCapacityDesign| == 1.
   if( std::abs( f_ConvMaxCapacityDesign ) == 1.0 )
    conv_design.set_type( ColVariable::kBinary );
   else
    conv_design.set_type( ColVariable::kInteger );
  }
  else
   conv_design.set_type( ColVariable::kNonNegative );
  add_static_variable( conv_design , "x_converter" );
  }
 else
  conv_design.set_value( dNaN );

 // storage level, intake, outtake- - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 v_storage_level.resize( f_time_horizon );
 for( auto & var : v_storage_level )
  var.set_type( ColVariable::kContinuous );
 add_static_variable( v_storage_level , "sl_battery" );

 if( needs_intake_outtake() ) {
  v_intake_level.resize( f_time_horizon );
  for( auto & var : v_intake_level )
   var.set_type( ColVariable::kNonNegative );
  add_static_variable( v_intake_level , "il_battery" );

  v_outtake_level.resize( f_time_horizon );
  for( auto & var : v_outtake_level )
   var.set_type( ColVariable::kNonNegative );
  add_static_variable( v_outtake_level , "ol_battery" );
  }

 // binary variables to avoid simultaneous charge and discharge - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( generate_binary_variables ) {
  v_battery_binary.resize( f_time_horizon );
  for( auto & var : v_battery_binary )
   var.set_type( ColVariable::kBinary );
  add_static_variable( v_battery_binary , "b_battery" );
  }

 // Active Power Variable - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 v_active_power.resize( f_time_horizon );
 for( auto & var : v_active_power )
  var.set_type( ColVariable::kContinuous );
 add_static_variable( v_active_power , "p_battery" );

 // Reactive Power Variable, if any - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( f_reactive_power ) {
  v_reactive_power.resize( f_time_horizon );
  for( auto & var : v_reactive_power )
   var.set_type( ColVariable::kNonNegative );
  add_static_variable( v_reactive_power , "q_battery" );
  }

 // Primary Spinning Reserve Variable - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( reserve_vars & 1u )  // if UCBlock has primary demand variables
  // if unit produces any primary reserve
  if( ! v_MaxPrimaryPower.empty() ) {
   v_primary_spinning_reserve.resize( f_time_horizon );
   for( auto & var : v_primary_spinning_reserve )
    var.set_type( ColVariable::kNonNegative );
   add_static_variable( v_primary_spinning_reserve , "pr_battery" );
   }

 // Secondary Spinning Reserve Variable - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( reserve_vars & 2u )  // if UCBlock has secondary demand variables
  // if unit produces any secondary reserve
  if( ! v_MaxSecondaryPower.empty() ) {
   v_secondary_spinning_reserve.resize( f_time_horizon );
   for( auto & var : v_secondary_spinning_reserve )
    var.set_type( ColVariable::kNonNegative );
   add_static_variable( v_secondary_spinning_reserve , "sc_battery" );
   }

 // variables for reference schedule, if there- - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! v_RefSchedule.empty() ) {
  v_abs_ref_schedule.resize( f_time_horizon );
  for( auto & var : v_abs_ref_schedule )
   var.set_type( ColVariable::kNonNegative );
  add_static_variable( v_abs_ref_schedule , "v_absb_refschd" );
  }

 set_variables_generated();

 }  // end( BatteryUnitBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::generate_abstract_constraints( Configuration * stcc )
{
 if( constraints_generated() )  // constraints have already been generated
  return;                       // nothing to do

 bool generate_ZOConstraints = false;
 if( ( ! stcc ) && f_BlockConfig )
  stcc = f_BlockConfig->f_static_constraints_Configuration;
 if( auto sci = dynamic_cast< SimpleConfiguration< int > * >( stcc ) )
  generate_ZOConstraints = sci->f_value;

 LinearFunction::v_coeff_pair vars;

 if( f_BattInvestmentCost == 0 ) {
  // no battery design variable - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // Intake outtake bounds constraints- - - - - - - - - - - - - - - - - - - -

  if( ! v_intake_level.empty() ) {

   intake_outtake_bounds_Const.resize(
    boost::multi_array< BoxConstraint , 2 >::extent_gen()
    [ 2 ][ f_time_horizon ] );  // 2 dims, i.e., intake and outtake

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {

    // set the maximum dispatch of converter not to exceed the C-rate of the
    // battery in charge
    intake_outtake_bounds_Const[ 0 ][ t ].set_rhs(
                              - f_kappa * f_MaxCRateCharge * v_MinPower[ t ] );
    intake_outtake_bounds_Const[ 0 ][ t ].set_variable( &v_intake_level[ t ] );

    // set the maximum dispatch of converter not to exceed the C-rate of the
    // battery in discharge
    intake_outtake_bounds_Const[ 1 ][ t ].set_rhs(
     f_kappa * f_MaxCRateDischarge * v_MaxPower[ t ] );
    intake_outtake_bounds_Const[ 1 ][ t ].set_variable( &v_outtake_level[ t ] );
    }
   }
  else {

   // the pair is the negative and the positive part of the active power, so
   // the two one-sided fences are the two sides of a single bound on it

   intake_outtake_bounds_Const.resize(
    boost::multi_array< BoxConstraint , 2 >::extent_gen()
    [ 1 ][ f_time_horizon ] );

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {
    intake_outtake_bounds_Const[ 0 ][ t ].set_lhs(
     f_kappa * f_MaxCRateCharge * v_MinPower[ t ] );
    intake_outtake_bounds_Const[ 0 ][ t ].set_rhs(
     f_kappa * f_MaxCRateDischarge * v_MaxPower[ t ] );
    intake_outtake_bounds_Const[ 0 ][ t ].set_variable( &v_active_power[ t ] );
    }
   }

  add_static_constraint( intake_outtake_bounds_Const ,
			 "IntakeOuttake_Battery" );

  // Active power bounds constraints- - - - - - - - - - - - - - - - - - - - -

  active_power_bounds_Const.resize(
   boost::multi_array< FRowConstraint , 2 >::extent_gen()
   [ 2 ][ f_time_horizon ] );  // 2 dims, i.e., the lower and upper bounds

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   vars.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );

   if( reserve_vars & 1u )  // if UCBlock has primary demand variables
    if( ! v_MaxPrimaryPower.empty() )
     // if this unit produces any primary reserve
     vars.push_back( std::make_pair( &v_primary_spinning_reserve[ t ] ,
                                     -1.0 ) );

   if( reserve_vars & 2u )  // if UCBlock has secondary demand variables
    if( ! v_MaxSecondaryPower.empty() )
     // if unit produces any secondary reserve
     vars.push_back( std::make_pair( &v_secondary_spinning_reserve[ t ] ,
                                     -1.0 ) );

   active_power_bounds_Const[ 0 ][ t ].set_lhs( f_kappa * v_MinPower[ t ] );
   active_power_bounds_Const[ 0 ][ t ].set_rhs( Inf< double >() );
   active_power_bounds_Const[ 0 ][ t ].set_function(
    new LinearFunction( std::move( vars ) ) );

   vars.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );

   if( reserve_vars & 1u )  // if UCBlock has primary demand variables
    if( ! v_MaxPrimaryPower.empty() )
     // if this unit produces any primary reserve
     vars.push_back( std::make_pair( &v_primary_spinning_reserve[ t ] ,
                                     1.0 ) );

   if( reserve_vars & 2u )  // if UCBlock has secondary demand variable
    if( ! v_MaxSecondaryPower.empty() )
     // if this unit produces any secondary reserve
     vars.push_back( std::make_pair( &v_secondary_spinning_reserve[ t ] ,
                                     1.0 ) );

   active_power_bounds_Const[ 1 ][ t ].set_lhs( -Inf< double >() );
   active_power_bounds_Const[ 1 ][ t ].set_rhs( f_kappa * v_MaxPower[ t ] );
   active_power_bounds_Const[ 1 ][ t ].set_function(
    new LinearFunction( std::move( vars ) ) );
   }

  add_static_constraint( active_power_bounds_Const , "ActivePower_Battery" );
  }
 else {
  // battery design variable- - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // Intake outtake upper bounds design constraints - - - - - - - - - - - - -

  const bool split = ! v_intake_level.empty();

  // the pair is fenced by three one-sided rows; folded onto the active power
  // the two C-rate rows keep their shape, while the converter row, which
  // reads il + ol <= C, becomes the two-sided -C <= p <= C: one more row
  // when the converter is designed, a plain bound when it is not
  const Index conv_rows = ( f_ConvInvestmentCost != 0 ) ? 2 : 0;

  intake_outtake_upper_bounds_design_Const.resize(
   boost::multi_array< FRowConstraint , 2 >::extent_gen()
   [ split ? 3 : ( 2 + conv_rows ) ][ f_time_horizon ] );

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   // Upper bound of the intake level design constraints:
   //
   //     v_intake_level <= ( -k C_ch v_MinPower ) x_b
   // => v_intake_level + ( k C_ch v_MinPower ) x_b <= 0
   //
   // folded onto the active power, v_intake_level = max( 0 , -p ):
   //
   //  => -p + ( k C_ch v_MinPower ) x_b <= 0

   // set the maximum dispatch of converter not to exceed the C-rate of the
   // battery in charge
   if( split )
    vars.push_back( std::make_pair( & v_intake_level[ t ] , 1.0 ) );
   else
    vars.push_back( std::make_pair( & v_active_power[ t ] , -1.0 ) );

   vars.push_back( std::make_pair( & batt_design ,
			   f_kappa * f_MaxCRateCharge * v_MinPower[ t ] ) );

   intake_outtake_upper_bounds_design_Const[ 0 ][ t ].set_lhs( -Inf< double >() );
   intake_outtake_upper_bounds_design_Const[ 0 ][ t ].set_rhs( 0.0 );
   intake_outtake_upper_bounds_design_Const[ 0 ][ t ].set_function(
                                    new LinearFunction( std::move( vars ) ) );

   // Upper bound of the outtake level design constraints:
   //
   //      v_outtake_level <= (k C_dis v_MaxPower ) x_b
   // => v_outtake_level - ( k C_dis v_MaxPower ) x_b <= 0
   //
   // folded onto the active power, v_outtake_level = max( 0 , p ):
   //
   //  => p - ( k C_dis v_MaxPower ) x_b <= 0

   // set the maximum dispatch of converter not to exceed the C-rate of the
   // battery in discharge
   if( split )
    vars.push_back( std::make_pair( & v_outtake_level[ t ] , 1.0 ) );
   else
    vars.push_back( std::make_pair( & v_active_power[ t ] , 1.0 ) );

   vars.push_back( std::make_pair( & batt_design ,
                        -f_kappa * f_MaxCRateDischarge * v_MaxPower[ t ] ) );

   intake_outtake_upper_bounds_design_Const[ 1 ][ t ].set_lhs( -Inf< double >() );
   intake_outtake_upper_bounds_design_Const[ 1 ][ t ].set_rhs( 0.0 );
   intake_outtake_upper_bounds_design_Const[ 1 ][ t ].set_function(
                                    new LinearFunction( std::move( vars ) ) );

   // Upper bound of the intake + outtake level design constraints:
   //
   //      v_intake_level + v_outtake_level <= ( k v_ConvMaxPower ) x_c
   // => v_intake_level + v_outtake_level - ( k v_ConvMaxPower ) x_c <= 0
   //
   // This constraint is generated only if the converter has a design variable
   // with a positive maximum power. Otherwise, a numerical bound or a
   // non-binding constraint is used to avoid disabling the battery.

   const double conv_power = v_ConvMaxPower.empty() ? 0.
                                                    : v_ConvMaxPower[ t ];

   if( split ) {
    if( ( f_ConvInvestmentCost != 0 ) && ( conv_power > 0 ) ) {
     // Case 1: converter is a design variable -> standard constraint with conv_design
     vars.push_back( std::make_pair( & v_intake_level[ t ] , 1 ) );
     vars.push_back( std::make_pair( & v_outtake_level[ t ] , 1 ) );
     vars.push_back( std::make_pair( & conv_design ,
                                     -f_kappa * conv_power ) );

     intake_outtake_upper_bounds_design_Const[ 2 ][ t ].set_lhs(
                                                           -Inf< double >() );
     intake_outtake_upper_bounds_design_Const[ 2 ][ t ].set_rhs( 0 );
     intake_outtake_upper_bounds_design_Const[ 2 ][ t ].set_function(
      new LinearFunction( std::move( vars ) ) );
     }
    else
     if( ( f_ConvInvestmentCost == 0 ) && ( conv_power > 0 ) ) {
      // Case 2: no converter design variable -> use numerical bound
      // v_intake_level + v_outtake_level <= k * v_ConvMaxPower[t]
      vars.push_back( std::make_pair( &v_intake_level[ t ] , 1 ) );
      vars.push_back( std::make_pair( &v_outtake_level[ t ] , 1 ) );

      intake_outtake_upper_bounds_design_Const[ 2 ][ t ].set_lhs(
                                                          -Inf< double >() );
      intake_outtake_upper_bounds_design_Const[ 2 ][ t ].set_rhs(
                                                 f_kappa * conv_power );
      intake_outtake_upper_bounds_design_Const[ 2 ][ t ].set_function(
                                   new LinearFunction( std::move( vars ) ) );
      }
     else {
      // Case 3: zero converter power -> deactivate constraint
      // Avoid generating "il + ol <= 0" which would switch off the battery
      intake_outtake_upper_bounds_design_Const[ 2 ][ t ].set_lhs(
                                                         -Inf< double >() );
      intake_outtake_upper_bounds_design_Const[ 2 ][ t ].set_rhs(
                                                          Inf< double >() );
      intake_outtake_upper_bounds_design_Const[ 2 ][ t ].set_function(
                                  new LinearFunction( std::move( vars ) ) );
      }
    }
   else
    if( f_ConvInvestmentCost != 0 ) {
     // folded onto the active power the converter fence is |p| <= C x_c,
     // i.e., the two rows p - ( k v_ConvMaxPower ) x_c <= 0 and
     // -p - ( k v_ConvMaxPower ) x_c <= 0; a zero converter power leaves
     // them non-binding, as it leaves the single row above
     const bool binding = ( conv_power > 0 );

     vars.push_back( std::make_pair( & v_active_power[ t ] , 1.0 ) );
     if( binding )
      vars.push_back( std::make_pair( & conv_design ,
                                      -f_kappa * conv_power ) );

     intake_outtake_upper_bounds_design_Const[ 2 ][ t ].set_lhs(
                                                          -Inf< double >() );
     intake_outtake_upper_bounds_design_Const[ 2 ][ t ].set_rhs(
                                     binding ? 0. : Inf< double >() );
     intake_outtake_upper_bounds_design_Const[ 2 ][ t ].set_function(
                                   new LinearFunction( std::move( vars ) ) );

     vars.push_back( std::make_pair( & v_active_power[ t ] , -1.0 ) );
     if( binding )
      vars.push_back( std::make_pair( & conv_design ,
                                      -f_kappa * conv_power ) );

     intake_outtake_upper_bounds_design_Const[ 3 ][ t ].set_lhs(
                                                          -Inf< double >() );
     intake_outtake_upper_bounds_design_Const[ 3 ][ t ].set_rhs(
                                     binding ? 0. : Inf< double >() );
     intake_outtake_upper_bounds_design_Const[ 3 ][ t ].set_function(
                                   new LinearFunction( std::move( vars ) ) );
     }
   }

  add_static_constraint( intake_outtake_upper_bounds_design_Const ,
                         "IntakeOuttake_Design_Battery" );

  // with no converter design variable the folded converter fence is the
  // plain bound -k v_ConvMaxPower <= p <= k v_ConvMaxPower

  if( ( ! split ) && ( f_ConvInvestmentCost == 0 ) &&
      std::any_of( v_ConvMaxPower.begin() , v_ConvMaxPower.end() ,
                   []( double power ) { return( power > 0 ); } ) ) {

   intake_outtake_bounds_Const.resize(
    boost::multi_array< BoxConstraint , 2 >::extent_gen()
    [ 1 ][ f_time_horizon ] );

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {
    if( v_ConvMaxPower[ t ] > 0 ) {
     intake_outtake_bounds_Const[ 0 ][ t ].set_lhs(
                                       -f_kappa * v_ConvMaxPower[ t ] );
     intake_outtake_bounds_Const[ 0 ][ t ].set_rhs(
                                        f_kappa * v_ConvMaxPower[ t ] );
     }
    else {
     intake_outtake_bounds_Const[ 0 ][ t ].set_lhs( -Inf< double >() );
     intake_outtake_bounds_Const[ 0 ][ t ].set_rhs( Inf< double >() );
     }

    intake_outtake_bounds_Const[ 0 ][ t ].set_variable( &v_active_power[ t ] );
    }

   add_static_constraint( intake_outtake_bounds_Const ,
                          "IntakeOuttake_Battery" );
   }

  // Active power bounds design constraints - - - - - - - - - - - - - - - - -

  active_power_bounds_design_Const.resize(
   boost::multi_array< FRowConstraint , 2 >::extent_gen()
   [ 2 ][ f_time_horizon ] );  // 2 dims, i.e., the lower and upper bounds

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   // Lower bound of the active power design constraints:
   //
   //      v_minimum_power x_b <= v_active_power
   // => 0 <= v_active_power - v_minimum_power x_b

   vars.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );

   if( reserve_vars & 1u )  // if UCBlock has primary demand variables
    if( ! v_MaxPrimaryPower.empty() )
     // if this unit produces any primary reserve
      vars.push_back( std::make_pair( &v_primary_spinning_reserve[ t ] ,
                                      -1.0 ) );

   if( reserve_vars & 2u )  // if UCBlock has secondary demand variables
    if( ! v_MaxSecondaryPower.empty() )
     // if unit produces any secondary reserve
      vars.push_back( std::make_pair( &v_secondary_spinning_reserve[ t ] ,
                                      -1.0 ) );

   vars.push_back( std::make_pair( &batt_design ,
                                   -f_kappa * v_MinPower[ t ] ) );

   active_power_bounds_design_Const[ 0 ][ t ].set_lhs( 0.0 );
   active_power_bounds_design_Const[ 0 ][ t ].set_rhs( Inf< double >() );
   active_power_bounds_design_Const[ 0 ][ t ].set_function(
    new LinearFunction( std::move( vars ) ) );

   // Upper bound of the active power design constraints:
   //
   //      v_active_power <= v_maximum_power x_b
   // => v_active_power - v_maximum_power x_b <= 0

   vars.push_back( std::make_pair( & v_active_power[ t ] , 1.0 ) );

   if( reserve_vars & 1u )  // if UCBlock has primary demand variables
    if( ! v_MaxPrimaryPower.empty() )
     // if this unit produces any primary reserve
      vars.push_back( std::make_pair( & v_primary_spinning_reserve[ t ] ,
                                      1.0 ) );

   if( reserve_vars & 2u )  // if UCBlock has secondary demand variable
    if( ! v_MaxSecondaryPower.empty() )
     // if this unit produces any secondary reserve
      vars.push_back( std::make_pair( & v_secondary_spinning_reserve[ t ] ,
                                      1.0 ) );

   vars.push_back( std::make_pair( &batt_design ,
                                   -f_kappa * v_MaxPower[ t ] ) );

   active_power_bounds_design_Const[ 1 ][ t ].set_lhs( -Inf< double >() );
   active_power_bounds_design_Const[ 1 ][ t ].set_rhs( 0.0 );
   active_power_bounds_design_Const[ 1 ][ t ].set_function(
                                 new LinearFunction( std::move( vars ) ) );
   }

  add_static_constraint( active_power_bounds_design_Const ,
                         "ActivePower_Design_Battery" );

  // Battery design bounds. f_BattMaxCapacityDesign < 0 selects integer
  // design with bound |f_BattMaxCapacityDesign|; >= 0 selects continuous
  // design with bound f_BattMaxCapacityDesign.
  const double lb_b = std::max( 0.0 , f_BattMinCapacityDesign );
  const bool is_integer_b = ( f_BattMaxCapacityDesign < 0.0 );
  const double ub_b = std::abs( f_BattMaxCapacityDesign );

  if( ( lb_b == 1.0 ) && ( ub_b == 1.0 ) )
   batt_design.is_unitary( true , eNoMod );
  else {
   batt_design_bound_Const.set_lhs( lb_b );
   batt_design_bound_Const.set_rhs( ub_b );
   batt_design_bound_Const.set_variable( &batt_design );

   add_static_constraint( batt_design_bound_Const ,
			  "BattDesignBound_Battery" );

   if( is_integer_b )
    batt_design.is_integer( true , eNoMod );
   }
  }  // end( battery design variables ) - - - - - - - - - - - - - - - - - - -
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( f_ConvInvestmentCost != 0 ) {
  // converter design variables - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // Converter design bounds. f_ConvMaxCapacityDesign < 0 selects integer
  // design with bound |f_ConvMaxCapacityDesign|; >= 0 selects continuous
  // design with bound f_ConvMaxCapacityDesign.
  const double lb_c = std::max( 0.0 , f_ConvMinCapacityDesign );
  const bool is_integer_c = ( f_ConvMaxCapacityDesign < 0.0 );
  const double ub_c = std::abs( f_ConvMaxCapacityDesign );

  if( ( lb_c == 1.0 ) && ( ub_c == 1.0 ) )
   conv_design.is_unitary( true , eNoMod );
  else {
   conv_design_bound_Const.set_lhs( lb_c );
   conv_design_bound_Const.set_rhs( ub_c );
   conv_design_bound_Const.set_variable( &conv_design );

   add_static_constraint( conv_design_bound_Const ,
			  "ConvDesignBound_Battery" );

   if( is_integer_c )
    conv_design.is_integer( true , eNoMod );
   }
  }

 // power_intake_outtake_Const- - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //
 // the row that defines the active power out of the pair is only there when
 // the pair is

 if( ! v_intake_level.empty() ) {

  power_intake_outtake_Const.resize( f_time_horizon );

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {

   vars.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );
   vars.push_back( std::make_pair( &v_intake_level[ t ] , 1.0 ) );
   vars.push_back( std::make_pair( &v_outtake_level[ t ] , -1.0 ) );

   power_intake_outtake_Const[ t ].set_both( 0.0 );
   power_intake_outtake_Const[ t ].set_function(
                                   new LinearFunction( std::move( vars ) ) );
   }

  add_static_constraint( power_intake_outtake_Const ,
                         "PowerIntakeOuttake_Battery" );
  }

 // Ramp-up constraints - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! v_DeltaRampUp.empty() ) {
  ramp_up_Const.resize( f_time_horizon );

  vars.push_back( std::make_pair( &v_active_power[ 0 ] , 1.0 ) );

  ramp_up_Const[ 0 ].set_lhs( -Inf< double >() );
  ramp_up_Const[ 0 ].set_rhs( v_DeltaRampUp[ 0 ] + f_InitialPower );
  ramp_up_Const[ 0 ].set_function( new LinearFunction( std::move( vars ) ) );

  for( Index t = 1 ; t < f_time_horizon ; ++t ) {
   vars.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );
   vars.push_back( std::make_pair( &v_active_power[ t - 1 ] , -1.0 ) );

   ramp_up_Const[ t ].set_lhs( -Inf< double >() );
   ramp_up_Const[ t ].set_rhs( v_DeltaRampUp[ t ] );
   ramp_up_Const[ t ].set_function( new LinearFunction( std::move( vars ) ) );
   }

  add_static_constraint( ramp_up_Const , "RampUp_Battery" );
  }

 // Ramp-down constraints - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! v_DeltaRampDown.empty() ) {
  ramp_down_Const.resize( f_time_horizon );

  vars.push_back( std::make_pair( &v_active_power[ 0 ] , 1.0 ) );

  ramp_down_Const[ 0 ].set_lhs( -v_DeltaRampDown[ 0 ] + f_InitialPower );
  ramp_down_Const[ 0 ].set_rhs( Inf< double >() );
  ramp_down_Const[ 0 ].set_function( new LinearFunction( std::move( vars ) ) );

  for( Index t = 1 ; t < f_time_horizon ; ++t ) {

   vars.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );
   vars.push_back( std::make_pair( &v_active_power[ t - 1 ] , -1.0 ) );

   ramp_down_Const[ t ].set_lhs( -v_DeltaRampDown[ t ] );
   ramp_down_Const[ t ].set_rhs( Inf< double >() );
   ramp_down_Const[ t ].set_function( new LinearFunction(
						      std::move( vars ) ) );
   }

  add_static_constraint( ramp_down_Const , "RampDown_Battery" );
  }

 // demand constraints- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 demand_Const.resize( f_time_horizon );

 vars.push_back( std::make_pair( &v_storage_level[ 0 ] , 1.0 ) );

 // -rho_st il + rho_ex ol is the energy the pair puts in the level; with
 // the pair folded onto the active power the two coefficients agree and
 // that is rho_ex ( ol - il ) = rho_ex p

 double intake_coeff = 1.0;
 if( ! v_StoringBatteryRho.empty() )
  intake_coeff = -v_StoringBatteryRho[ 0 ];

 double outtake_coeff = 1.0;
 if( ! v_ExtractingBatteryRho.empty() )
  outtake_coeff = v_ExtractingBatteryRho[ 0 ];

 if( ! v_intake_level.empty() ) {
  vars.push_back( std::make_pair( &v_intake_level[ 0 ] , intake_coeff ) );
  vars.push_back( std::make_pair( &v_outtake_level[ 0 ] , outtake_coeff ) );
  }
 else
  vars.push_back( std::make_pair( &v_active_power[ 0 ] , outtake_coeff ) );

 double standing_coeff = 1.0;
  if( ! v_StandingBatteryRho.empty() )
    standing_coeff = v_StandingBatteryRho[ 0 ];

 if( f_InitialStorage < 0 )  // cyclical notation
  vars.push_back( std::make_pair( &v_storage_level[ f_time_horizon - 1 ] ,
                                  -standing_coeff ) );

 if( ! v_Demand.empty() )
  demand_Const[ 0 ].set_both(
        ( f_InitialStorage < 0 ? 0.0 : f_InitialStorage ) - v_Demand[ 0 ] );
 else
  demand_Const[ 0 ].set_both(
        ( f_InitialStorage < 0 ? 0.0 : f_InitialStorage ) );

 demand_Const[ 0 ].set_function( new LinearFunction( std::move( vars ) ) );

 for( Index t = 1 ; t < f_time_horizon ; ++t ) {  
  vars.push_back( std::make_pair( &v_storage_level[ t ] , 1.0 ) );
  
  standing_coeff = 1.0;
  if( ! v_StandingBatteryRho.empty() )
    standing_coeff = v_StandingBatteryRho[ t ];
  vars.push_back( std::make_pair( &v_storage_level[ t - 1 ] , -standing_coeff ) );

  intake_coeff = 1.0;
  if( ! v_StoringBatteryRho.empty() )
   intake_coeff = -v_StoringBatteryRho[ t ];

  outtake_coeff = 1.0;
  if( ! v_ExtractingBatteryRho.empty() )
   outtake_coeff = v_ExtractingBatteryRho[ t ];

  if( ! v_intake_level.empty() ) {
   vars.push_back( std::make_pair( &v_intake_level[ t ] , intake_coeff ) );
   vars.push_back( std::make_pair( &v_outtake_level[ t ] , outtake_coeff ) );
   }
  else
   vars.push_back( std::make_pair( &v_active_power[ t ] , outtake_coeff ) );

  if( ! v_Demand.empty() )
   demand_Const[ t ].set_both( -v_Demand[ t ] );
  else
   demand_Const[ t ].set_both( 0.0 );

  demand_Const[ t ].set_function( new LinearFunction( std::move( vars ) ) );
  }

 add_static_constraint( demand_Const , "Demand_Battery" );

 if( f_BattInvestmentCost == 0 ) {
  // no battery design variables- - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // Storage level bound constraints- - - - - - - - - - - - - - - - - - - - -

  storage_level_bounds_Const.resize( f_time_horizon );

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   storage_level_bounds_Const[ t ].set_lhs( f_kappa * v_MinStorage[ t ] );
   storage_level_bounds_Const[ t ].set_rhs( f_kappa * v_MaxStorage[ t ] );
   storage_level_bounds_Const[ t ].set_variable( &v_storage_level[ t ] );
   }

  add_static_constraint( storage_level_bounds_Const ,
			 "StorageLevel_Battery" );
  }
 else {
  // battery design variables - - - - - - - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  // Storage level bound design constraints - - - - - - - - - - - - - - - - -

  storage_level_bounds_design_Const.resize(
   boost::multi_array< FRowConstraint , 2 >::extent_gen()
   [ 2 ][ f_time_horizon ] );  // 2 dims, i.e., the lower and upper bounds

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   // Lower bound of the storage level design constraints:
   //
   //      v_MinStorage x_b <= v_storage_level
   // => 0 <= v_storage_level - v_MinStorage x_b

   vars.push_back( std::make_pair( & v_storage_level[ t ] , 1.0 ) );
   vars.push_back( std::make_pair( & batt_design ,
                                   -f_kappa * v_MinStorage[ t ] ) );

   storage_level_bounds_design_Const[ 0 ][ t ].set_lhs( 0.0 );
   storage_level_bounds_design_Const[ 0 ][ t ].set_rhs( Inf< double >() );
   storage_level_bounds_design_Const[ 0 ][ t ].set_function(
                                  new LinearFunction( std::move( vars ) ) );

   // Upper bound of the storage level design constraints:
   //
   //      v_storage_level <= v_MaxStorage x_b
   // => v_storage_level - v_MaxStorage x_b <= 0
   vars.push_back( std::make_pair( & v_storage_level[ t ] , 1.0 ) );
   vars.push_back( std::make_pair( & batt_design ,
                                   -f_kappa * v_MaxStorage[ t ] ) );

   storage_level_bounds_design_Const[ 1 ][ t ].set_lhs( -Inf< double >() );
   storage_level_bounds_design_Const[ 1 ][ t ].set_rhs( 0.0 );
   storage_level_bounds_design_Const[ 1 ][ t ].set_function(
                                  new LinearFunction( std::move( vars ) ) );
   }

  add_static_constraint( storage_level_bounds_design_Const ,
                         "StorageLevel_Design_Battery" );

  }  // end( battery design variables ) - - - - - - - - - - - - - - - - - - -
     // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // intake_outtake_binary constraints - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! v_battery_binary.empty() ) {
  intake_outtake_binary_Const.resize(
   boost::multi_array< FRowConstraint , 2 >::extent_gen()
   [ 2 ][ f_time_horizon ] );  // 2 dims, i.e., the intake and outtake bounds

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   //      v_intake_level <= - v_MinPower b
   // => v_intake_level + v_MinPower b <= 0
   // (v_MinPower < 0 by battery convention, so -v_MinPower b is the
   //  positive intake bound when b == 1)
   vars.push_back( std::make_pair( & v_intake_level[ t ] , 1.0 ) );
   vars.push_back( std::make_pair( & v_battery_binary[ t ] ,
                                   f_kappa * v_MinPower[ t ] ) );

   intake_outtake_binary_Const[ 0 ][ t ].set_lhs( -Inf< double >() );
   intake_outtake_binary_Const[ 0 ][ t ].set_rhs( 0.0 );
   intake_outtake_binary_Const[ 0 ][ t ].set_function(
                                    new LinearFunction( std::move( vars ) ) );

   //      v_outtake_level <= v_MaxPower ( 1 - b )
   // => v_outtake_level + v_MaxPower b <= v_MaxPower
   vars.push_back( std::make_pair( & v_outtake_level[ t ] , 1.0 ) );
   vars.push_back( std::make_pair( & v_battery_binary[ t ] ,
                                   f_kappa * v_MaxPower[ t ] ) );

   intake_outtake_binary_Const[ 1 ][ t ].set_lhs( -Inf< double >() );
   intake_outtake_binary_Const[ 1 ][ t ].set_rhs( f_kappa * v_MaxPower[ t ] );
   intake_outtake_binary_Const[ 1 ][ t ].set_function(
                                    new LinearFunction( std::move( vars ) ) );
   }

  add_static_constraint( intake_outtake_binary_Const ,
                         "Intake_Outtake_Binary_Battery" );
  }

 // primary_upper_bound constraints - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( reserve_vars & 1u )  // if UCBlock has primary demand variables
  if( ! v_MaxPrimaryPower.empty() ) {
   // if this unit produces any primary reserve

   primary_upper_bound_Const.resize( f_time_horizon );

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {
    primary_upper_bound_Const[ t ].set_rhs(
					 f_kappa * v_MaxPrimaryPower[ t ] );
    primary_upper_bound_Const[ t ].set_variable(
					& v_primary_spinning_reserve[ t ] );
    }

   add_static_constraint( primary_upper_bound_Const ,
                          "Primary_UpperBound_Battery" );
   }

 // secondary_upper_bound constraints - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( reserve_vars & 2u )  // if UCBlock has secondary demand variables
  if( ! v_MaxSecondaryPower.empty() ) {
   // if this unit produces any secondary reserve

   secondary_upper_bound_Const.resize( f_time_horizon );

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {
    secondary_upper_bound_Const[ t ].set_rhs(
					f_kappa * v_MaxSecondaryPower[ t ] );
    secondary_upper_bound_Const[ t ].set_variable(
				       & v_secondary_spinning_reserve[ t ] );
    }

   add_static_constraint( secondary_upper_bound_Const ,
                          "Secondary_UpperBound_Battery" );
   }

 // ZOConstraint- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // these are useless, but some approaches do require the bounds to be
 // "physically" there, e.g., to compute dual variables
 
 if( ! v_battery_binary.empty() )
  if( generate_ZOConstraints ) {  // the battery binary bound constraints
   battery_binary_bound_Const.resize( f_time_horizon );

   for( Index t = 0 ; t < f_time_horizon ; ++t )
    battery_binary_bound_Const[ t ].set_variable( &v_battery_binary[ t ] );

   add_static_constraint( battery_binary_bound_Const , "Binary_Battery" );
   }

 // reference schedule constraints- - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! v_RefSchedule.empty() ) {
  Reference_Schedule_Const.resize( 2 * f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   // | P - Pref | <= v_abs_ref_schedule
   auto lfunc_1 = new LinearFunction();
   lfunc_1->add_variable( & v_active_power[ t ], 1.0 );
   lfunc_1->add_variable( & v_abs_ref_schedule[ t ], -1.0 );
   Reference_Schedule_Const[ t ].set_lhs( -Inf< double >() );
   Reference_Schedule_Const[ t ].set_rhs( v_RefSchedule[ t ] );
   Reference_Schedule_Const[ t ].set_function( lfunc_1 );
   //
   auto lfunc_2 = new LinearFunction();
   lfunc_2->add_variable( & v_active_power[ t ], -1.0 );
   lfunc_2->add_variable( & v_abs_ref_schedule[ t ], -1.0 );
   Reference_Schedule_Const[ f_time_horizon + t ].set_lhs( -Inf< double >() );
   Reference_Schedule_Const[ f_time_horizon + t ].set_rhs( -v_RefSchedule[ t ] );
   Reference_Schedule_Const[ f_time_horizon + t ].set_function( lfunc_2 );
   }

  add_static_constraint( Reference_Schedule_Const ,
			 "Norm1B_Reference_Schedule" );
  }

 // reactive power bounds constraints (if any)- - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( f_reactive_power && 
     ( ( ! v_MinReactivePower.empty() ) || ( ! v_MaxReactivePower.empty() ) )
     ) {  // reactive-related stuff- - - - - - - - - - - - - - - - - - - - -
  // Reactive power bounds constraints
  if( ReactivePower_Bound_Const.size() != f_time_horizon )
   ReactivePower_Bound_Const.resize( f_time_horizon );

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   ReactivePower_Bound_Const[ t ].set_rhs( get_max_reactive_power( t ) );
   ReactivePower_Bound_Const[ t ].set_lhs( get_min_reactive_power( t ) );
   ReactivePower_Bound_Const[ t ].set_variable( & v_reactive_power[ t ] );
   }

  add_static_constraint( ReactivePower_Bound_Const , "ReactivePowerBound" );

  /*!! Link between active and reactive power
  Reactive_2_Active_Const.resize( f_time_horizon );

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   // Q(t) - P(t) <= 0
   auto lfunc = new LinearFunction();
   lfunc->add_variable( &v_active_power[ t ], -1.0 );
   lfunc->add_variable( &v_reactive_power[ t ], 1.0 );
    
   Reactive_2_Active_Const[ t ].set_lhs( -Inf< double >() );
   Reactive_2_Active_Const[ t ].set_rhs( 0.0 );
   Reactive_2_Active_Const[ t ].set_function( lfunc );
   }

  add_static_constraint( Reactive_2_Active_Const, "QandPbattery" );
  !!*/
  }

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 set_constraints_generated();

 }  // end( BatteryUnitBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::generate_objective( Configuration *objc )
{
 if( objective_generated() )  // Objective has already been generated
  return;                     // nothing to do

 auto lf = new LinearFunction();

 // the two operating costs are opposite, hence they fold onto the active
 // power as -Cost p whenever the pair is not there

 if( v_RefSchedule.empty() ) {
  if( ! v_intake_level.empty() )
   for( Index t = 0 ; t < f_time_horizon ; ++t ) {
    lf->add_variable( &v_intake_level[ t ] , f_scale * v_Cost[ t ] , eNoMod );
    lf->add_variable( &v_outtake_level[ t ] , -f_scale * v_Cost[ t ] ,
                      eNoMod );
    }
  else
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    lf->add_variable( &v_active_power[ t ] , -f_scale * v_Cost[ t ] , eNoMod );
  }


 if( f_BattInvestmentCost != 0 )
  lf->add_variable( &batt_design , f_scale * f_BattInvestmentCost );

 if( f_ConvInvestmentCost != 0 )
  lf->add_variable( &conv_design , f_scale * f_ConvInvestmentCost );

 if( ! v_RefSchedule.empty() ) {
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   lf->add_variable( &v_abs_ref_schedule[ t ] , 1.0 );
 }

 objective.set_function( lf );
 objective.set_sense( Objective::eMin );

 // Set Block objective
 this->set_objective( &objective , eNoMod );

 set_objective_generated();

 }  // end( BatteryUnitBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*---------------- METHODS FOR CHECKING THE BatteryUnitBlock ---------------*/
/*--------------------------------------------------------------------------*/

bool BatteryUnitBlock::is_feasible( bool useabstract , Configuration * fsbc )
{
 // Retrieve the tolerance and the type of violation.
 double tol = 0;
 bool rel_viol = true;

 // Try to extract, from "c", the parameters that determine feasibility.
 // If it succeeds, it sets the values of the parameters and returns
 // true. Otherwise, it returns false.
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

 if( ( ! extract_parameters( fsbc ) ) && f_BlockConfig )
  // if the given Configuration is not valid, try the one from the BlockConfig
  extract_parameters( f_BlockConfig->f_is_feasible_Configuration );

 return(
  UnitBlock::is_feasible( useabstract )
  // Variables
  && ColVariable::is_feasible( v_storage_level , tol )
  && ColVariable::is_feasible( v_intake_level , tol )
  && ColVariable::is_feasible( v_outtake_level , tol )
  && ColVariable::is_feasible( v_battery_binary , tol )
  && ColVariable::is_feasible( v_active_power , tol )
  && ColVariable::is_feasible( v_reactive_power , tol )
  && ColVariable::is_feasible( v_primary_spinning_reserve , tol )
  && ColVariable::is_feasible( v_secondary_spinning_reserve , tol )
  // Constraints: notice that the ZOConstraint are not checked, since the
  // corresponding check is made on the ColVariable
  && RowConstraint::is_feasible( batt_design_bound_Const , tol , rel_viol )
  && RowConstraint::is_feasible( conv_design_bound_Const , tol , rel_viol )
  && RowConstraint::is_feasible( active_power_bounds_Const , tol , rel_viol )
  && RowConstraint::is_feasible( intake_outtake_upper_bounds_design_Const , tol , rel_viol )
  && RowConstraint::is_feasible( active_power_bounds_design_Const , tol , rel_viol )
  && RowConstraint::is_feasible( storage_level_bounds_design_Const , tol , rel_viol )
  && RowConstraint::is_feasible( intake_outtake_binary_Const , tol , rel_viol )
  && RowConstraint::is_feasible( power_intake_outtake_Const , tol , rel_viol )
  && RowConstraint::is_feasible( ramp_up_Const , tol , rel_viol )
  && RowConstraint::is_feasible( ramp_down_Const , tol , rel_viol )
  && RowConstraint::is_feasible( demand_Const , tol , rel_viol )
  && RowConstraint::is_feasible( storage_level_bounds_Const , tol , rel_viol )
  && RowConstraint::is_feasible( intake_outtake_bounds_Const , tol , rel_viol )
  && RowConstraint::is_feasible( primary_upper_bound_Const , tol , rel_viol )
  && RowConstraint::is_feasible( secondary_upper_bound_Const , tol , rel_viol )
  && RowConstraint::is_feasible( Reference_Schedule_Const , tol , rel_viol )
  && RowConstraint::is_feasible( ReactivePower_Bound_Const , tol , rel_viol ) );

} // end( BatteryUnitBlock::is_feasible )

/*--------------------------------------------------------------------------*/
/*------- METHODS FOR LOADING, PRINTING & SAVING THE BatteryUnitBlock ------*/
/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::serialize( netCDF::NcGroup & group ) const {

 UnitBlock::serialize( group );

 // Serialize scalar variables

 ::serialize( group , "InitialPower" , netCDF::NcDouble() , f_InitialPower );
 ::serialize( group , "InitialStorage" , netCDF::NcDouble() , f_InitialStorage );
 ::serialize( group , "Kappa" , netCDF::NcDouble() , f_kappa );

 if( f_BattInvestmentCost != 0 ) {
  ::serialize( group , "BatteryInvestmentCost" , netCDF::NcDouble() ,
               f_BattInvestmentCost );
  if( f_BattMinCapacityDesign != 0 )
   ::serialize( group , "BatteryMinCapacityDesign" , netCDF::NcDouble() ,
                f_BattMinCapacityDesign );
  if( f_BattMaxCapacityDesign != 1 )
   ::serialize( group , "BatteryMaxCapacityDesign" , netCDF::NcDouble() ,
                f_BattMaxCapacityDesign );
 }

 if( f_ConvInvestmentCost != 0 ) {
  ::serialize( group , "ConverterInvestmentCost" , netCDF::NcDouble() ,
               f_ConvInvestmentCost );
  if( f_ConvMinCapacityDesign != 0 )
   ::serialize( group , "ConverterMinCapacityDesign" , netCDF::NcDouble() ,
                f_ConvMinCapacityDesign );
  if( f_ConvMaxCapacityDesign != 1 )
   ::serialize( group , "ConverterMaxCapacityDesign" , netCDF::NcDouble() ,
                f_ConvMaxCapacityDesign );
 }

 if( f_scale != 1 )
  ::serialize( group , "Scale" , netCDF::NcDouble() , f_scale );

 if( f_MaxCRateCharge != 1 )
  ::serialize( group , "MaxCRateCharge" , netCDF::NcDouble() ,
               f_MaxCRateCharge );

 if( f_MaxCRateDischarge != 1 )
  ::serialize( group , "MaxCRateDischarge" , netCDF::NcDouble() ,
               f_MaxCRateDischarge );

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
   throw( std::logic_error(
    "BatteryUnitBlock::serialize: invalid dimension for variable " +
    var_name + ": " + std::to_string( data.size() ) +
    ". Its dimension must be one of the following: TimeHorizon, "
    "NumberIntervals, 1." ) );

  ::serialize( group , var_name , ncType , dimension , data ,
               allow_scalar_var );
 };

 serialize( "MinStorage" , v_MinStorage );
 serialize( "MaxStorage" , v_MaxStorage );
 serialize( "MinPower" , v_MinPower );
 serialize( "MaxPower" , v_MaxPower );
 serialize( "ConverterMaxPower" , v_ConvMaxPower );
 serialize( "MaxPrimaryPower" , v_MaxPrimaryPower );
 serialize( "MaxSecondaryPower" , v_MaxSecondaryPower );
 serialize( "DeltaRampUp" , v_DeltaRampUp );
 serialize( "DeltaRampDown" , v_DeltaRampDown );
 serialize( "StoringBatteryRho" , v_StoringBatteryRho );
 serialize( "ExtractingBatteryRho" , v_ExtractingBatteryRho );
 serialize( "StandingBatteryRho" , v_StandingBatteryRho );
 serialize( "Cost" , v_Cost );
 serialize( "Demand" , v_Demand );

}  // end( BatteryUnitBlock::serialize )

/*--------------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/

Solution * BatteryUnitBlock::get_Solution( Configuration * csolc ,
					   bool emptys )
{
 Index wsol = 63;
 if( ( ! csolc ) && f_BlockConfig )
  csolc = f_BlockConfig->f_solution_Configuration;

 if( auto config = dynamic_cast< SimpleConfiguration< int > * >( csolc ) )
  wsol = config->f_value;

 // call the method of the base class
 auto * sol = dynamic_cast< BatteryUnitBlockSolution * >(
		                  UnitBlock::get_Solution( csolc , emptys ) );
 assert( sol );

 if( wsol & 16 )
  sol->v_storage.resize( get_time_horizon() );

 if( wsol & 32 )
  sol->v_intake.resize( get_time_horizon() );

 if( ! emptys )
  sol->read( this );

 return( sol );
 }

/*--------------------------------------------------------------------------*/

UnitBlockSolution * BatteryUnitBlock::new_Solution( void ) const {
 return( new BatteryUnitBlockSolution() );
 }

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::update_initial_storage_in_cnstrs( c_ModParam issueAMod )
{
 if( demand_Const.empty() )
  return;

 if( ! v_Demand.empty() )
  demand_Const[ 0 ].set_both(
   ( f_InitialStorage < 0 ? 0.0 : f_InitialStorage ) - v_Demand[ 0 ] ,
   issueAMod );
 else
  demand_Const[ 0 ].set_both(
   ( f_InitialStorage < 0 ? 0.0 : f_InitialStorage ) ,
   issueAMod );

}  // end( BatteryUnitBlock::update_initial_storage_in_cnstrs )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::set_initial_storage( MF_dbl_it values ,
                                            Subset && subset ,
                                            const bool ordered ,
                                            c_ModParam issuePMod ,
                                            c_ModParam issueAMod )
{
 if( subset.empty() )
  return;

 // Find the last index 0
 auto index_it = std::find( subset.rbegin() , subset.rend() , 0 );

 if( index_it == subset.rend() )
  return;  // 0 is not in subset; return

 std::advance( values , std::distance( index_it , subset.rend() ) - 1 );

 if( f_InitialStorage == *values )
  return;

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation
  f_InitialStorage = *values;

  if( not_dry_run( issueAMod ) && constraints_generated() )
   // Change the abstract representation
   update_initial_storage_in_cnstrs( issueAMod );
 }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< BatteryUnitBlockMod >(
                            this , BatteryUnitBlockMod::eSetInitS ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( BatteryUnitBlock::set_initial_storage( subset ) )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::set_initial_storage( MF_dbl_it values ,
                                            Range rng ,
                                            c_ModParam issuePMod ,
                                            c_ModParam issueAMod )
{
 rng.second = std::min( rng.second , static_cast< decltype( rng.second ) >( 1 ) );
 if( ! ( ( rng.first <= 0 ) && ( 0 < rng.second ) ) )
  return;  // 0 does not belong to the range; return

 std::advance( values , -rng.first );

 if( f_InitialStorage == *values )
  return;  // nothing changes; return

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation
  f_InitialStorage = *values;

  if( not_dry_run( issueAMod ) && constraints_generated() )
   // Change the abstract representation
   update_initial_storage_in_cnstrs( issueAMod );
 }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< BatteryUnitBlockMod >(
                            this , BatteryUnitBlockMod::eSetInitS ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( BatteryUnitBlock::set_initial_storage( range ) )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::update_initial_power_in_cnstrs( c_ModParam issueAMod )
{
 if( ! ( ramp_up_Const.empty() || v_DeltaRampUp.empty() ) )
  ramp_up_Const[ 0 ].set_rhs( v_DeltaRampUp[ 0 ] + f_InitialPower ,
                              issueAMod );

 if( ! ( ramp_down_Const.empty() || v_DeltaRampDown.empty() ) )
  ramp_down_Const[ 0 ].set_lhs( -v_DeltaRampDown[ 0 ] + f_InitialPower ,
                                issueAMod );

}  // end( BatteryUnitBlock::update_initial_power_in_cnstrs )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::set_initial_power( MF_dbl_it values ,
                                          Subset && subset ,
                                          const bool ordered ,
                                          c_ModParam issuePMod ,
                                          c_ModParam issueAMod )
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

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation
  f_InitialPower = *values;

  if( not_dry_run( issueAMod ) && constraints_generated() )
   // Change the abstract representation
   update_initial_power_in_cnstrs( issueAMod );
 }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< BatteryUnitBlockMod >(
                            this , BatteryUnitBlockMod::eSetInitP ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( BatteryUnitBlock::set_initial_power( subset ) )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::set_initial_power( MF_dbl_it values ,
                                          Range rng ,
                                          c_ModParam issuePMod ,
                                          c_ModParam issueAMod )
{
 rng.second = std::min( rng.second , static_cast< decltype( rng.second ) >( 1 ) );
 if( ! ( ( rng.first <= 0 ) && ( 0 < rng.second ) ) )
  return;  // 0 does not belong to the range; return

 std::advance( values , -rng.first );

 if( f_InitialPower == *values )
  return;  // nothing changes; return

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation
  f_InitialPower = *values;

  if( not_dry_run( issueAMod ) && constraints_generated() )
   // Change the abstract representation
   update_initial_power_in_cnstrs( issueAMod );
 }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< BatteryUnitBlockMod >(
                            this , BatteryUnitBlockMod::eSetInitP ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( BatteryUnitBlock::set_initial_power( range ) )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::set_cost( MF_dbl_it values ,
                                 Subset && subset ,
                                 const bool ordered ,
                                 ModParam issuePMod ,
                                 ModParam issueAMod )
{
 if( subset.empty() )
  return;  // Since the given Subset is empty, no operation is performed

 if( v_Cost.empty() ) {
  if( std::all_of( values ,
                   values + subset.size() ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  v_Cost.assign( f_time_horizon , 0.0 );
 }

 Subset actual_subset = subset;
 if( ! ordered )
  std::sort( actual_subset.begin() , actual_subset.end() );

 if( actual_subset.back() >= v_Cost.size() )
  throw( std::invalid_argument(
   "BatteryUnitBlock::set_cost: invalid index in subset." ) );

 auto values_it = values;
 bool identical = true;
 for( auto t : actual_subset ) {
  if( v_Cost[ t ] != *( values_it++ ) ) {
   identical = false;
   break;
  }
 }

 if( identical )
  return;  // nothing changes; return

 if( not_dry_run( issuePMod ) ) {
  values_it = values;
  for( auto t : actual_subset )
   v_Cost[ t ] = *( values_it++ );

  if( not_dry_run( issueAMod ) && objective_generated() )
   // Update the abstract representation
   update_objective( issueAMod );
 }

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification
  if( ! ordered )
   std::sort( subset.begin() , subset.end() );

  Block::add_Modification( std::make_shared< BatteryUnitBlockSbstMod >(
                            this , BatteryUnitBlockMod::eSetCost ,
                            std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );
 }
}  // end( BatteryUnitBlock::set_cost( subset ) )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::set_cost( MF_dbl_it values ,
                                 Range rng ,
                                 ModParam issuePMod ,
                                 ModParam issueAMod )
{
 rng.second = std::min( rng.second , f_time_horizon );
 if( rng.second <= rng.first )
  return;  // An empty Range was given: no operation is performed.

 c_Index sz = rng.second - rng.first;

 if( v_Cost.empty() ) {
  if( std::all_of( values ,
                   values + sz ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  v_Cost.assign( f_time_horizon , 0.0 );
 }

 if( std::equal( values , values + sz , v_Cost.begin() + rng.first ) )
  return;  // nothing changes; return

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation
  std::copy( values ,
             values + sz ,
             v_Cost.begin() + rng.first );

  if( not_dry_run( issueAMod ) && objective_generated() )
   // Update the abstract representation
   update_objective( issueAMod );
 }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
   Block::add_Modification( std::make_shared< BatteryUnitBlockRngdMod >(
                             this , BatteryUnitBlockMod::eSetCost , rng ) ,
                            Observer::par2chnl( issuePMod ) );

}  // end( BatteryUnitBlock::set_cost( range ) )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::scale( MF_dbl_it values ,
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
    update_objective( issueAMod );
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

}  // end( BatteryUnitBlock::scale )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::update_kappa_in_cnstrs( ModParam issueAMod )
{
 // with no intake/outtake pair the fences are stated on the active power:
 // the C-rate pair is a single bound when the battery is not designed, and
 // the converter fence is a bound when the converter is not
 const bool split = ! v_intake_level.empty();

 if( ! intake_outtake_bounds_Const.empty() ) {

  if( split )
   for( Index t = 0 ; t < f_time_horizon ; ++t ) {
    intake_outtake_bounds_Const[ 0 ][ t ].set_rhs(
     -f_kappa * f_MaxCRateCharge * v_MinPower[ t ] , issueAMod );
    intake_outtake_bounds_Const[ 1 ][ t ].set_rhs(
     f_kappa * f_MaxCRateDischarge * v_MaxPower[ t ] , issueAMod );
   }
  else
   if( f_BattInvestmentCost == 0 )
    for( Index t = 0 ; t < f_time_horizon ; ++t ) {
     intake_outtake_bounds_Const[ 0 ][ t ].set_lhs(
      f_kappa * f_MaxCRateCharge * v_MinPower[ t ] , issueAMod );
     intake_outtake_bounds_Const[ 0 ][ t ].set_rhs(
      f_kappa * f_MaxCRateDischarge * v_MaxPower[ t ] , issueAMod );
    }
   else
    for( Index t = 0 ; t < f_time_horizon ; ++t )
     if( v_ConvMaxPower[ t ] > 0 ) {
      intake_outtake_bounds_Const[ 0 ][ t ].set_lhs(
       -f_kappa * v_ConvMaxPower[ t ] , issueAMod );
      intake_outtake_bounds_Const[ 0 ][ t ].set_rhs(
       f_kappa * v_ConvMaxPower[ t ] , issueAMod );
     }
  }

 if( ! intake_outtake_upper_bounds_design_Const.empty() )

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {

   auto f0 = static_cast< LinearFunction * >(
    intake_outtake_upper_bounds_design_Const[ 0 ][ t ].get_function() );

   const auto batt_design_idx0 = f0->is_active( &batt_design );

   if( batt_design_idx0 == Inf< Index >() )
    throw( std::logic_error( "BatteryUnitBlock::update_kappa_in_cnstrs: "
                             "expected Variable not found in "
                             "intake_outtake_upper_bounds_design_Const." ) );

   f0->modify_coefficient( batt_design_idx0 ,
                           f_kappa * f_MaxCRateCharge * v_MinPower[ t ] ,
                           issueAMod );

   auto f1 = static_cast< LinearFunction * >(
    intake_outtake_upper_bounds_design_Const[ 1 ][ t ].get_function() );

   const auto batt_design_idx1 = f1->is_active( &batt_design );

   if( batt_design_idx1 == Inf< Index >() )
    throw( std::logic_error( "BatteryUnitBlock::update_kappa_in_cnstrs: "
                             "expected Variable not found in "
                             "intake_outtake_upper_bounds_design_Const." ) );

   f1->modify_coefficient( batt_design_idx1 ,
                           -f_kappa * f_MaxCRateDischarge * v_MaxPower[ t ] ,
                           issueAMod );

   if( v_ConvMaxPower.empty() || ( v_ConvMaxPower[ t ] <= 0 ) )
    continue;  // the converter fence is not binding

   if( f_ConvInvestmentCost != 0 ) {
    // one row on the pair, two on the active power, all of them holding
    // the same coefficient of the converter design variable
    const Index rows = split ? 1 : 2;

    for( Index row = 2 ; row < 2 + rows ; ++row ) {
     auto f2 = static_cast< LinearFunction * >(
      intake_outtake_upper_bounds_design_Const[ row ][ t ].get_function() );

     const auto conv_design_idx2 = f2->is_active( &conv_design );

     if( conv_design_idx2 == Inf< Index >() )
      throw( std::logic_error( "BatteryUnitBlock::update_kappa_in_cnstrs: "
                               "expected Variable not found in "
                               "intake_outtake_upper_bounds_design_Const." ) );

     f2->modify_coefficient( conv_design_idx2 ,
                             -f_kappa * v_ConvMaxPower[ t ] ,
                             issueAMod );
     }
    }
   else
    if( split )
     intake_outtake_upper_bounds_design_Const[ 2 ][ t ].set_rhs(
      f_kappa * v_ConvMaxPower[ t ] , issueAMod );
   }

 if( ! active_power_bounds_Const.empty() )

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   active_power_bounds_Const[ 0 ][ t ].set_lhs(
    f_kappa * v_MinPower[ t ] , issueAMod );
   active_power_bounds_Const[ 1 ][ t ].set_rhs(
    f_kappa * v_MaxPower[ t ] , issueAMod );
  }

 else if( ! active_power_bounds_design_Const.empty() )

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   auto f0 = static_cast< LinearFunction * >(
    active_power_bounds_design_Const[ 0 ][ t ].get_function() );

   const auto batt_design_idx0 = f0->is_active( &batt_design );

   if( batt_design_idx0 == Inf< Index >() )
    throw( std::logic_error( "BatteryUnitBlock::update_kappa_in_cnstrs: "
                             "expected Variable not found in "
                             "active_power_bounds_design_Const." ) );

   f0->modify_coefficient( batt_design_idx0 ,
                           -f_kappa * v_MinPower[ t ] ,
                           issueAMod );

   auto f1 = static_cast< LinearFunction * >(
    active_power_bounds_design_Const[ 1 ][ t ].get_function() );

   const auto batt_design_idx1 = f1->is_active( &batt_design );

   if( batt_design_idx1 == Inf< Index >() )
    throw( std::logic_error( "BatteryUnitBlock::update_kappa_in_cnstrs: "
                             "expected Variable not found in "
                             "active_power_bounds_design_Const." ) );

   f1->modify_coefficient( batt_design_idx1 ,
                           -f_kappa * v_MaxPower[ t ] ,
                           issueAMod );
  }

 if( ! storage_level_bounds_Const.empty() )

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   storage_level_bounds_Const[ t ].set_lhs(
    f_kappa * v_MinStorage[ t ] , issueAMod );
   storage_level_bounds_Const[ t ].set_rhs(
    f_kappa * v_MaxStorage[ t ] , issueAMod );
  }

 else if( ! storage_level_bounds_design_Const.empty() )

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {

   auto f0 = static_cast< LinearFunction * >(
    storage_level_bounds_design_Const[ 0 ][ t ].get_function() );

   const auto batt_design_idx0 = f0->is_active( &batt_design );

   if( batt_design_idx0 == Inf< Index >() )
    throw( std::logic_error( "BatteryUnitBlock::update_kappa_in_cnstrs: "
                             "expected Variable not found in "
                             "storage_level_bounds_design_Const." ) );

   f0->modify_coefficient( batt_design_idx0 ,
                           -f_kappa * v_MinStorage[ t ] ,
                           issueAMod );

   auto f1 = static_cast< LinearFunction * >(
    storage_level_bounds_design_Const[ 1 ][ t ].get_function() );

   const auto batt_design_idx1 = f1->is_active( &batt_design );

   if( batt_design_idx1 == Inf< Index >() )
    throw( std::logic_error( "BatteryUnitBlock::update_kappa_in_cnstrs: "
                             "expected Variable not found in "
                             "storage_level_bounds_design_Const." ) );

   f1->modify_coefficient( batt_design_idx1 ,
                           -f_kappa * v_MaxStorage[ t ] ,
                           issueAMod );
  }

 if( ( ! intake_outtake_binary_Const.empty() ) &&
     ( ! intake_outtake_binary_Const[ 0 ].empty() ) )

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {

   auto f = static_cast< LinearFunction * >(
    intake_outtake_binary_Const[ 0 ][ t ].get_function() );

   const auto index = f->is_active( &v_battery_binary[ t ] );

   if( index == Inf< Index >() )
    throw( std::logic_error( "BatteryUnitBlock::set_kappa: expected Variable"
                             "not found in intake_binary_Const." ) );

   // intake bound is - v_MinPower b (v_MinPower is the charging-side limit,
   // negative by convention): coefficient on v_battery_binary is + kappa
   // v_MinPower, see the construction at the top of the file
   f->modify_coefficient( index , f_kappa * v_MinPower[ t ] , issueAMod );
  }

 if( ( ! intake_outtake_binary_Const.empty() ) &&
     ( ! intake_outtake_binary_Const[ 1 ].empty() ) )

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {

   auto f = static_cast< LinearFunction * >
    ( intake_outtake_binary_Const[ 1 ][ t ].get_function() );

   const auto index = f->is_active( &v_battery_binary[ t ] );

   if( index == Inf< Index >() )
    throw( std::logic_error( "BatteryUnitBlock::set_kappa: expected Variable"
                             "not found in outtake_binary_Const." ) );

   // outtake bound is v_MaxPower (1 - b): coefficient on v_battery_binary
   // is + kappa v_MaxPower and the RHS is kappa v_MaxPower, see the
   // construction at the top of the file
   f->modify_coefficient( index , f_kappa * v_MaxPower[ t ] , issueAMod );

   intake_outtake_binary_Const[ 1 ][ t ].set_rhs(
    f_kappa * v_MaxPower[ t ] , issueAMod );
  }

 if( ! primary_upper_bound_Const.empty() )
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   primary_upper_bound_Const[ t ].set_rhs( f_kappa * v_MaxPrimaryPower[ t ] ,
                                           issueAMod );

 if( ! secondary_upper_bound_Const.empty() )
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   secondary_upper_bound_Const[ t ].set_rhs( f_kappa * v_MaxSecondaryPower[ t ] ,
                                             issueAMod );

 }  // end( BatteryUnitBlock::update_kappa_in_cnstrs )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::set_kappa( MF_dbl_it values ,
                                  Subset && subset ,
                                  const bool ordered ,
                                  ModParam issuePMod ,
                                  ModParam issueAMod )
{
 if( subset.empty() )
  return;  // Since the given Subset is empty, no operation is performed

 if( f_kappa == *values )
  return;  // The kappa constant does not change: nothing to do

 if( not_dry_run( issuePMod ) ) {
  f_kappa = *values;  // Update the kappa constant

  if( not_dry_run( issueAMod ) )
   // Update the abstract representation
   if( constraints_generated() )
    // Update the constraints
    update_kappa_in_cnstrs( issueAMod );
 }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< BatteryUnitBlockMod >(
                            this , BatteryUnitBlockMod::eSetKappa ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( BatteryUnitBlock::set_kappa( subset ) )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::set_kappa( MF_dbl_it values ,
                                  Range rng ,
                                  ModParam issuePMod ,
                                  ModParam issueAMod )
{
 if( rng.first >= rng.second )
  return;  // An empty Range was given: no operation is performed.

 Subset subset( 1 , 0 );

 set_kappa( values , std::move( subset ) , true , issuePMod , issueAMod );

}  // end( BatteryUnitBlock::set_kappa( range ) )

/*--------------------------------------------------------------------------*/

double BatteryUnitBlock::get_kappa_linearization( void ) const {

 /* The kappa constant of a BatteryUnitBlock appears in the following
  * Constraints for each time instant t:
  *
  * - Minimum and maximum power output constraint (lambda):
  *
  *   kappa * P^{min}_{t} <= p^{ac}_{t} - p^{pr}_{t} - p^{sc}_{t}  [lambda_min]
  *
  *   p^{ac}_{t} + p^{pr}_{t} + p^{sc}_{t} <= kappa * P^{max}_{t}  [lambda_max]
  *
  * - Intake and outtake level bounds (alpha), in the form the unit has: one
  *   one-sided bound per variable when intake and outtake are kept apart,
  *
  *   p^{+}_{t} <= - kappa * C^{c} * P^{min}_{t}                   [alpha_in]
  *
  *   p^{-}_{t} <= kappa * C^{d} * P^{max}_{t}                     [alpha_out]
  *
  *   or, when they are not, the single two-sided bound they collapse into,
  *
  *   kappa * C^{c} * P^{min}_{t} <= p^{ac}_{t} <= kappa * C^{d} * P^{max}_{t}
  *
  *   which for a battery carrying an investment is instead the converter
  *   fence, where the converter has a power at all,
  *
  *   - kappa * P^{conv}_{t} <= p^{ac}_{t} <= kappa * P^{conv}_{t}
  *
  *   whose multiplier belongs to the side its sign points at [alpha]
  *
  *   p^{+}_{t} <= kappa * u^{+}_t * P^{max}_{t}                   [alpha_max_u]
  *
  *   p^{-}_{t} <= - kappa * (1 - u^{+}_t) * P^{min}_{t}           [alpha_min_u]
  *
  * - Storage level bounds (beta):
  *
  *   kappa * V^{min}_t <= v_t                                     [beta_min]
  *
  *   v_t <= kappa * V^{max}_t                                     [beta_max]
  *
  * - Primary and secondary reserves bounds (gamma):
  *
  *   p^{pr}_t <= kappa P^{pr max}_t                               [gamma_pr]
  *
  *   p^{sc}_t <= kappa P^{sc max}_t                               [gamma_sc]
  *
  * The name between [] represents the dual variable associated with each
  * constraint. The linearization coefficient is
  *
  *   P^{min} ' (lambda_min + C^{c} alpha_in + (1 - u^+) * alpha_min_u) -
  *   P^{max} ' (lambda_max + C^{d} alpha_out + u^+ * alpha_max_u) +
  *   V^{min} ' beta_min - V^{max} ' beta_max -
  *   P^{pr max} ' gamma_pr - P^{sc max} ' gamma_sc
  *
  * with the C-rates there because they multiply kappa in the bounds.
  */

 double linearization = 0;

 // Minimum and maximum power output constraint

 const auto & min_power_constraints = get_min_power_constraints();

 const auto & max_power_constraints = get_max_power_constraints();

 // Intake and outtake level bounds

 const auto & intake_bound_constraints = get_max_intake_bounds();

 const auto & outtake_bound_constraints = get_max_outtake_bounds();

 const auto & max_intake_binary_constraints =
  get_max_intake_binary_constraints();

 const auto & max_outtake_binary_constraints =
  get_max_outtake_binary_constraints();

 const auto & u = get_intake_outtake_binary_variables();

 // Storage level bounds

 const auto & storage_level_bound_constraints = get_storage_level_bounds();

 // Primary and secondary reserves bounds

 const auto & primary_reserve_bounds = get_primary_reserve_bounds();

 const auto & secondary_reserve_bounds = get_secondary_reserve_bounds();

 /* The dual value of a constraint that has both finite lower and upper bounds
  * is associated with either the lower bound or the upper bound
  * constraint. This will help determine to which bound the dual is associated
  * with. */
 const auto obj_sign =
  ( get_objective_sense() == Objective::eMin ) ? - 1 : 1;

 const auto time_horizon = get_time_horizon();

 for( Index t = 0 ; t < time_horizon ; ++t ) {

  const auto min_power = get_min_power( t );
  const auto max_power = get_max_power( t );
  const auto min_storage = get_min_storage()[ t ];
  const auto max_storage = get_max_storage()[ t ];

  // Minimum and maximum power output constraint

  const auto lambda_min = std::abs( min_power_constraints[ t ].get_dual() );
  const auto lambda_max = std::abs( max_power_constraints[ t ].get_dual() );

  linearization += min_power * lambda_min - max_power * lambda_max;

  /* Intake and outtake level bounds, in the very form the unit states them
   * [see update_kappa_in_cnstrs()]: two one-sided bounds, one per variable,
   * when intake and outtake are kept apart; one two-sided bound on the
   * active power when they are not, whose fences are the C-rate pair for a
   * battery that carries no investment and the converter fence for one that
   * does. The coefficient is the multiplier times the derivative of the
   * fence the Constraint actually holds, which is why the C-rates and the
   * converter power belong here: they multiply kappa in it. */

  if( intake_bound_constraints ) {
   if( outtake_bound_constraints ) {
    linearization += intake_bound_constraints[ t ].get_dual() *
                     f_MaxCRateCharge * min_power;

    linearization += - outtake_bound_constraints[ t ].get_dual() *
                     f_MaxCRateDischarge * max_power;
    }
   else {
    const bool converter = ( f_BattInvestmentCost != 0 ) &&
                           ( t < v_ConvMaxPower.size() ) &&
                           ( v_ConvMaxPower[ t ] > 0 );

    const auto lower = converter ? - v_ConvMaxPower[ t ]
                                 : f_MaxCRateCharge * min_power;
    const auto upper = converter ? v_ConvMaxPower[ t ]
                                 : f_MaxCRateDischarge * max_power;

    const auto alpha = intake_bound_constraints[ t ].get_dual();
    linearization += - alpha * ( ( obj_sign * alpha > 0 ) ? lower : upper );
    }
   }

  if( max_intake_binary_constraints ) {
   const auto alpha_max_u =
    std::abs( max_intake_binary_constraints[ t ].get_dual() );
   linearization += - alpha_max_u * u[ t ].get_value() * max_power;
  }

  if( max_outtake_binary_constraints ) {
   const auto alpha_min_u =
    std::abs( max_outtake_binary_constraints[ t ].get_dual() );
   linearization += ( 1.0 - u[ t ].get_value() ) * alpha_min_u * min_power;
  }

  // Storage level bounds

  const auto dual = storage_level_bound_constraints[ t ].get_dual();
  auto bound = max_storage;
  if( obj_sign * dual > 0 ) {
   // The bound is associated with the lower bound constraint
   bound = min_storage;
  }

  // The bound is associated with the upper bound constraint.
  linearization += - dual * bound;

  // Primary and secondary reserves bounds

  if( ! primary_reserve_bounds.empty() ) {
   const auto gamma_pr = std::abs( primary_reserve_bounds[ t ].get_dual() );
   linearization += - get_max_primary_power()[ t ] * gamma_pr;
  }

  if( ! secondary_reserve_bounds.empty() ) {
   const auto gamma_sc = std::abs( secondary_reserve_bounds[ t ].get_dual() );
   linearization += - get_max_secondary_power()[ t ] * gamma_sc;
  }

 }

 return( linearization );
 }  // end( BatteryUnitBlock::get_kappa_linearization )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::update_objective( c_ModParam issueAMod ) const {

 if( ! objective_generated() )
  return;  // the Objective has not been generated: nothing to be done

 auto function = static_cast< LinearFunction * >( objective.get_function() );

 // intake / outtake operating costs are emitted only when v_RefSchedule is
 // empty (see generate_objective); refresh them only in that case
 if( v_RefSchedule.empty() ) {
  const bool split = ! v_intake_level.empty();
  const Index terms = ( split ? 2 : 1 ) * f_time_horizon;

  LinearFunction::Vec_FunctionValue coefficients;
  coefficients.reserve( terms );

  for( Index t = 0 ; t < f_time_horizon ; ++t )
   if( split ) {
    coefficients.push_back( f_scale * v_Cost[ t ] );
    coefficients.push_back( -f_scale * v_Cost[ t ] );
    }
   else
    coefficients.push_back( -f_scale * v_Cost[ t ] );

  function->modify_coefficients( std::move( coefficients ) ,
                                 Range( 0 , terms ) ,
                                 issueAMod );
 }

 // refresh the scale-aware design coefficients (battery / converter)
 if( f_BattInvestmentCost != 0 ) {
  const auto idx = function->is_active( & batt_design );
  assert( idx < function->get_num_active_var() );
  function->modify_coefficient( idx ,
                                f_scale * f_BattInvestmentCost ,
                                issueAMod );
 }

 if( f_ConvInvestmentCost != 0 ) {
  const auto idx = function->is_active( & conv_design );
  assert( idx < function->get_num_active_var() );
  function->modify_coefficient( idx ,
                                f_scale * f_ConvInvestmentCost ,
                                issueAMod );
 }

}  // end( BatteryUnitBlock::update_objective )

/*--------------------------------------------------------------------------*/
/*------------------ METHODS OF BatteryUnitBlockSolution -------------------*/
/*--------------------------------------------------------------------------*/

void BatteryUnitBlockSolution::deserialize( const netCDF::NcGroup & group )
{
 // call the method of the base class
 UnitBlockSolution::deserialize( group );

 if( f_number_generators != 1 )
  throw( std::logic_error( "BatteryUnitBlockSolution::deserialize: "
			   "batteries have only one generator" ) );

 // deserialize the storage levels - - - - - - - - - - - - - - - - - - - - -
 ::deserialize< double >( group , "StorageLevel" , v_storage , false );

 // deserialize the intakes- - - - - - - - - - - - - - - - - - - - - - - - -
 ::deserialize< double >( group , "InOutTake" , v_intake , false );

 // deserialize the battery design - - - - - - - - - - - - - - - - - - - - -

 if( ! ::deserialize< double >( group , f_b_design , "BatteryDesign" ) )
  f_b_design = dNaN;

 // deserialize the converter design - - - - - - - - - - - - - - - - - - - -

 if( ! ::deserialize< double >( group , f_c_design , "ConverterDesign" ) )
  f_c_design = dNaN;

 }  // end( BatteryUnitBlockSolution::deserialize )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlockSolution::read( const Block * block )
{
 auto BUB = dynamic_cast< const BatteryUnitBlock * >( block );
 if( ! BUB )
  throw( std::invalid_argument(
       "BatteryUnitBlockSolution::read: block is not a BatteryUnitBlock" ) );

 UnitBlockSolution::read( BUB );  // call the method of the base class

 if( ! v_storage.empty() ) {
  // read the storage levels - - - - - - - - - - - - - - - - - - - - - - - -
  auto SLit = BUB->get_const_storage_level().begin();
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   v_storage[ t ] = ( SLit++ )->get_value();
  }

 if( ! v_intake.empty() ) {
  // read the intakes- - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // the intake is the net one, hence it is minus the active power when the
  // pair is not there
  if( ! BUB->get_const_intake_level().empty() ) {
   auto Iit = BUB->get_const_intake_level().begin();
   auto Oit = BUB->get_const_outtake_level().begin();
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    v_intake[ t ] = ( Iit++ )->get_value() - ( Oit++ )->get_value();
   }
  else {
   auto Pit = BUB->get_const_active_power( 0 );
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    v_intake[ t ] = - ( Pit++ )->get_value();
   }
  }

 // read the battery design- - - - - - - - - - - - - - - - - - - - - - - - -
 f_b_design = BUB->get_const_batt_design().get_value();

 // read the converter design- - - - - - - - - - - - - - - - - - - - - - - -
 f_c_design = BUB->get_const_conv_design().get_value();

 }  // end( BatteryUnitBlockSolution::read )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlockSolution::write( Block * block )
{
 UnitBlockSolution::write( block );  // call the method of the base class

 auto BUB = dynamic_cast< BatteryUnitBlock * >( block );
 if( ! BUB )
  throw( std::invalid_argument(
       "BatteryUnitBlockSolution::read: block is not a BatteryUnitBlock" ) );

 if( ! v_storage.empty() ) {
  // write the storage levels- - - - - - - - - - - - - - - - - - - - - - - -
  auto SLit = BUB->get_storage_level().begin();
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   ( SLit++ )->set_value( v_storage[ t ] );
  }

 if( ! v_intake.empty() ) {
  // write the intakes - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if( ! BUB->get_intake_level().empty() ) {
   auto Iit = BUB->get_intake_level().begin();
   auto Oit = BUB->get_outtake_level().begin();
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    if( v_intake[ t ] >= 0 ) {
     ( Iit++ )->set_value( v_intake[ t ] );
     ( Oit++ )->set_value( 0 );
     }
    else {
     ( Iit++ )->set_value( 0 );
     ( Oit++ )->set_value( -v_intake[ t ] );
     }
   }
  else {
   // with no pair the net intake is minus the active power
   auto Pit = BUB->get_active_power( 0 );
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    ( Pit++ )->set_value( -v_intake[ t ] );
   }
  }

 // write the battery design - - - - - - - - - - - - - - - - - - - - - - - -
 BUB->get_batt_design().set_value( f_b_design );

 // write the converter design - - - - - - - - - - - - - - - - - - - - - - -
 BUB->get_conv_design().set_value( f_c_design );

 }  // end( BatteryUnitBlockSolution::write )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlockSolution::serialize( netCDF::NcGroup & group ) const
{
 UnitBlockSolution::serialize( group );  // call the method of the base class

 // recover the just serialized time horizon
 netCDF::NcDim th = group.getDim( "TimeHorizon" );

 // serialize the storage levels- - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_storage.empty() )
  ::serialize< double >( group , "StorageLevel" , netCDF::NcDouble() , th ,
			 v_storage );

 // serialize the intake- - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_intake.empty() )
  ::serialize< double >( group , "InOutTake" , netCDF::NcDouble() , th ,
			 v_intake );

 // serialize the battery design- - - - - - - - - - - - - - - - - - - - - - -
 if( ! std::isnan( f_b_design ) )
  ::serialize< double >( group , "BatteryDesign" , netCDF::NcDouble() ,
			 f_b_design );

 // serialize the converter design- - - - - - - - - - - - - - - - - - - - - -
 if( ! std::isnan( f_c_design ) )
  ::serialize< double >( group , "ConverterDesign" , netCDF::NcDouble() ,
			 f_c_design );

 }  // end( BatteryUnitBlockSolution::serialize )

/*--------------------------------------------------------------------------*/

BatteryUnitBlockSolution * BatteryUnitBlockSolution::scale( double factor )
 const
{
 auto sol = clone();

 if( factor == 1 )
  return( sol );

 guts_of_scale( sol , factor );

 if( ! v_storage.empty() )
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   sol->v_storage[ t ] *= factor;

 if( ! v_intake.empty() )
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   sol->v_intake[ t ] *= factor;

 if( ! std::isnan( f_b_design ) )
  sol->f_b_design *= factor;

 if( ! std::isnan( f_c_design ) )
  sol->f_c_design *= factor;

 return( sol );

 }  // end( BatteryUnitBlockSolution::scale )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlockSolution::sum( const Solution * solution ,
				    double multiplier )
{
 // call the method of the base class
 UnitBlockSolution::sum( solution , multiplier );

 auto BUBS = dynamic_cast< const BatteryUnitBlockSolution * >( solution );
 if( ! BUBS )
  throw( std::invalid_argument( "BatteryUnitBlockSolution::sum: solution "
				"not a BatteryUnitBlockSolution" ) );

 if( ! v_storage.empty() )
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   v_storage[ t ] += BUBS->v_storage[ t ] * multiplier;

 if( ! v_intake.empty() )
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   v_intake[ t ] += BUBS->v_intake[ t ] * multiplier;

 if( ! std::isnan( f_b_design ) )
  f_b_design += BUBS->f_b_design * multiplier;

 if( ! std::isnan( f_c_design ) )
  f_c_design += BUBS->f_c_design * multiplier;

 }  // end( BatteryUnitBlockSolution::sum )

/*--------------------------------------------------------------------------*/

BatteryUnitBlockSolution * BatteryUnitBlockSolution::clone( bool empty )
 const
{
 auto * sol = new BatteryUnitBlockSolution();

 if( ! empty ) {
  guts_of_clone( sol );
  sol->v_storage = v_storage;
  sol->v_intake = v_intake;
  sol->f_b_design = f_b_design;
  sol->f_c_design = f_c_design;
  }

 return( sol );

 }  // end( BatteryUnitBlockSolution::clone )

/*--------------------------------------------------------------------------*/
/*----------------- End File BatteryUnitBlock.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
