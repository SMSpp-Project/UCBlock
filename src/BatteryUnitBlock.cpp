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

 Constraint::clear( intake_outtake_bounds_Const );
 Constraint::clear( primary_upper_bound_Const );
 Constraint::clear( secondary_upper_bound_Const );
 Constraint::clear( storage_level_bounds_Const );

 Constraint::clear( battery_binary_bound_Const );

 objective.clear();
}

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::deserialize( const netCDF::NcGroup & group )
{

#ifndef NDEBUG
 std::vector< std::string > expected_dims = { "TimeHorizon" ,
                                              "NumberIntervals" };
 check_dimensions( group , expected_dims , std::cerr );

 std::vector< std::string > expected_vars = { "MinStorage" , "MaxStorage" ,
                                              "InitialStorage" ,
                                              "MinPower" , "MaxPower" ,
                                              "InitialPower" ,
                                              "ConverterMaxPower" ,
                                              "MaxPrimaryPower" ,
                                              "MaxSecondaryPower" ,
                                              "DeltaRampUp" , "DeltaRampDown" ,
                                              "StoringBatteryRho" ,
                                              "ExtractingBatteryRho" ,
                                              "Cost" , "Demand" , "Kappa" ,
                                              "MaxCRateCharge" ,
                                              "MaxCRateDischarge" ,
                                              "BatteryMaxCapacity" ,
                                              "ConverterMaxCapacity" ,
                                              "MaxIntakePower" ,
                                              "MaxOuttakePower" ,
                                              "BatteryInvestmentCost" ,
                                              "ConverterInvestmentCost" };
 check_variables( group , expected_vars , std::cerr );
#endif

 // Deserialize data from the base class
 UnitBlock::deserialize( group );

 // Mandatory variables

 ::deserialize( group , "MinStorage" , v_MinStorage , false );
 ::deserialize( group , "MaxStorage" , v_MaxStorage , false );

 ::deserialize( group , "MaxPower" , v_MaxPower , false );

 // Optional variables

 if( ! ::deserialize( group , "MinPower" , v_MinPower ) ) {
  v_MinPower.resize( v_MaxPower.size() );
  std::copy( v_MaxPower.begin() , v_MaxPower.end() , v_MinPower.begin() );
  std::transform( v_MinPower.cbegin() , v_MinPower.cend() , v_MinPower.begin() ,
                  []( double p ) { return -p; } );
 }

 if( ! ::deserialize( group , "ConverterMaxPower" , v_ConvMaxPower ) ) {
  v_ConvMaxPower.resize( v_MaxPower.size() );
  std::copy( v_MaxPower.begin() , v_MaxPower.end() , v_ConvMaxPower.begin() );
 }

 ::deserialize( group , f_MaxCRateCharge , "MaxCRateCharge" );
 ::deserialize( group , f_MaxCRateDischarge , "MaxCRateDischarge" );

 ::deserialize( group , f_InitialStorage , "InitialStorage" );

 ::deserialize( group , f_InitialPower , "InitialPower" );

 ::deserialize( group , f_kappa , "Kappa" );

 ::deserialize( group , "MaxPrimaryPower" , v_MaxPrimaryPower );
 ::deserialize( group , "MaxSecondaryPower" , v_MaxSecondaryPower );

 ::deserialize( group , "DeltaRampUp" , v_DeltaRampUp );
 ::deserialize( group , "DeltaRampDown" , v_DeltaRampDown );

 ::deserialize( group , "Demand" , v_Demand );

 ::deserialize( group , "StoringBatteryRho" , v_StoringBatteryRho );
 ::deserialize( group , "ExtractingBatteryRho" , v_ExtractingBatteryRho );

 if( ! ::deserialize( group , "Cost" , v_Cost ) )
  v_Cost.resize( 1 );

 ::deserialize( group , f_BattInvestmentCost , "BatteryInvestmentCost" );
 ::deserialize( group , f_ConvInvestmentCost , "ConverterInvestmentCost" );

 ::deserialize( group , f_BattMaxCapacity , "BatteryMaxCapacity" );
 ::deserialize( group , f_ConvMaxCapacity , "ConverterMaxCapacity" );

 // Decompress vectors

 decompress_vector( v_MinPower );
 decompress_vector( v_MaxPower );
 decompress_vector( v_ConvMaxPower );
 decompress_vector( v_MinStorage );
 decompress_vector( v_MaxStorage );
 decompress_vector( v_MaxPrimaryPower );
 decompress_vector( v_MaxSecondaryPower );
 decompress_vector( v_DeltaRampUp );
 decompress_vector( v_DeltaRampDown );
 decompress_vector( v_StoringBatteryRho );
 decompress_vector( v_ExtractingBatteryRho );
 decompress_vector( v_Demand );
 decompress_vector( v_Cost );

 check_data_consistency();

}  // end( BatteryUnitBlock::deserialize )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::check_data_consistency( void ) const {

 // InvestmentCost

 if( ( ( f_BattInvestmentCost != 0 ) || ( f_ConvInvestmentCost != 0 ) ) &&
     ( f_InitialStorage >= 0 ) )
  throw( std::logic_error( "BatteryUnitBlock::check_data_consistency: the "
                           "presence of the investment cost of the battery "
                           "allows the model to switch into the strategic "
                           "scenario mode, but the presence of also a positive "
                           "initial storage, typical of the operative "
                           "scenario, is incompatible." ) );

 // Minimum and maximum power

 assert( v_MinPower.size() == f_time_horizon );
 assert( v_MaxPower.size() == f_time_horizon );

 for( Index t = 0 ; t < f_time_horizon ; ++t )
  if( v_MinPower[ t ] > v_MaxPower[ t ] )
   throw( std::logic_error( "BatteryUnitBlock::check_data_consistency: minimum "
                            "power for time " + std::to_string( t ) + " is " +
                            std::to_string( v_MinPower[ t ] ) + ", which "
                            "greater than the maximum power, which is " +
                            std::to_string( v_MaxPower[ t ] ) + "." ) );

 // Minimum and maximum storage levels

 assert( v_MinStorage.size() == f_time_horizon );
 assert( v_MaxStorage.size() == f_time_horizon );

 for( Index t = 0 ; t < f_time_horizon ; ++t )
  if( ( v_MinStorage[ t ] > v_MaxStorage[ t ] ) ||
      ( v_MinStorage[ t ] < 0 ) )
   throw( std::logic_error( "BatteryUnitBlock::check_data_consistency: maximum "
                            "and minimum storage levels must be such that "
                            "maximum_storage >= minimum_storage >= 0." ) );

 // Inefficiency of storing and extracting energy

 if( ! v_StoringBatteryRho.empty() ) {
  assert( v_StoringBatteryRho.size() == f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   if( v_StoringBatteryRho[ t ] > 1 )
    throw( std::logic_error( "BatteryUnitBlock::check_data_consistency: invalid"
                             " inefficiency of storing energy for time step " +
                             std::to_string( t ) + ": " +
                             std::to_string( v_StoringBatteryRho[ t ] ) +
                             ". It must not be greater than 1." ) );
 }

 if( ! v_ExtractingBatteryRho.empty() ) {
  assert( v_ExtractingBatteryRho.size() == f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   if( v_ExtractingBatteryRho[ t ] < 1 )
    throw( std::logic_error( "BatteryUnitBlock::check_data_consistency: invalid"
                             " inefficiency of extracting energy for time "
                             "step " + std::to_string( t ) + ": " +
                             std::to_string( v_ExtractingBatteryRho[ t ] ) +
                             ". It must not be less than 1." ) );
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
                                  ", but is must be nonnegative." ) );
 }

 // Kappa

 if( f_kappa < 0 )
  throw( std::logic_error( "BatteryUnitBlock::check_data_consistency: "
                           "kappa must be nonnegative, but it is" +
                           std::to_string( f_kappa ) + "." ) );

}  // end( BatteryUnitBlock::check_data_consistency )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::generate_abstract_variables( Configuration * stvv )
{
 if( variables_generated() )  // variables have already been generated
  return;                     // nothing to do

 UnitBlock::generate_abstract_variables( stvv );

 // Check if negative prices may occur
 bool negative_prices = false;
 if( ( ! stvv ) && f_BlockConfig )
  stvv = f_BlockConfig->f_static_variables_Configuration;
 if( auto sci = dynamic_cast< SimpleConfiguration< int > * >( stvv ) )
  negative_prices = sci->f_value;

 // Binary variables must be generated if negative prices may occur and if
 // there is some t such that
 // StoringBatteryRho[ t ] < 1 < ExtractingBatterRho[ t ]

 bool generate_binary_variables = false;

 if( negative_prices && ( ! v_StoringBatteryRho.empty() ) &&
     ( ! v_ExtractingBatteryRho.empty() ) ) {
  assert( v_StoringBatteryRho.size() == f_time_horizon );
  assert( v_ExtractingBatteryRho.size() == f_time_horizon );
  for( Index t = 0 ; t < v_StoringBatteryRho.size() ; ++t )
   if( ( v_StoringBatteryRho[ t ] < 1 ) &&
       ( v_ExtractingBatteryRho[ t ] > 1 ) ) {
    generate_binary_variables = true;
    break;
   }
 }

 // Add the static variables

 v_storage_level.resize( f_time_horizon );
 for( auto & var : v_storage_level )
  var.set_type( ColVariable::kNonNegative );
 add_static_variable( v_storage_level , "sl_battery" );

 v_intake_level.resize( f_time_horizon );
 for( auto & var : v_intake_level )
  var.set_type( ColVariable::kNonNegative );
 add_static_variable( v_intake_level , "il_battery" );

 v_outtake_level.resize( f_time_horizon );
 for( auto & var : v_outtake_level )
  var.set_type( ColVariable::kNonNegative );
 add_static_variable( v_outtake_level , "ol_battery" );

 if( generate_binary_variables ) {
  v_battery_binary.resize( f_time_horizon );
  for( auto & var : v_battery_binary )
   var.set_type( ColVariable::kBinary );
  add_static_variable( v_battery_binary , "b_battery" );
 }

 // Battery Design Variable
 if( f_BattInvestmentCost != 0 ) {
  batt_design.set_type( ColVariable::kPosUnitary );
  add_static_variable( batt_design , "x_battery" );
 }

 // Converter Design Variable
 if( f_ConvInvestmentCost != 0 ) {
  conv_design.set_type( ColVariable::kPosUnitary );
  add_static_variable( conv_design , "x_converter" );
 }

 // Active Power Variable
 v_active_power.resize( f_time_horizon );
 for( auto & var : v_active_power )
  var.set_type( ColVariable::kContinuous );
 add_static_variable( v_active_power , "p_battery" );

 // Primary Spinning Reserve Variable
 if( reserve_vars & 1u )  // if UCBlock has primary demand variables
  // if unit produces any primary reserve
  if( ! v_MaxPrimaryPower.empty() ) {
   v_primary_spinning_reserve.resize( f_time_horizon );
   for( auto & var : v_primary_spinning_reserve )
    var.set_type( ColVariable::kNonNegative );
   add_static_variable( v_primary_spinning_reserve , "pr_battery" );
  }

 // Secondary Spinning Reserve Variable
 if( reserve_vars & 2u )  // if UCBlock has secondary demand variables
  // if unit produces any secondary reserve
  if( ! v_MaxSecondaryPower.empty() ) {
   v_secondary_spinning_reserve.resize( f_time_horizon );
   for( auto & var : v_secondary_spinning_reserve )
    var.set_type( ColVariable::kNonNegative );
   add_static_variable( v_secondary_spinning_reserve , "sc_battery" );
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

  // Intake outtake bounds constraints

  intake_outtake_bounds_Const.resize(
   boost::multi_array< BoxConstraint , 2 >::extent_gen()
   [ 2 ][ f_time_horizon ] );  // 2 dims, i.e., the intake outtake upper bounds

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {

   // set the maximum dispatch of converter not to exceed the C-rate of the
   // battery in discharge
   intake_outtake_bounds_Const[ 0 ][ t ].set_rhs(
    -f_kappa * f_MaxCRateDischarge * v_MinPower[ t ] );
   intake_outtake_bounds_Const[ 0 ][ t ].set_variable( &v_intake_level[ t ] );

   // set the maximum dispatch of converter not to exceed the C-rate of the
   // battery in charge
   intake_outtake_bounds_Const[ 1 ][ t ].set_rhs(
    f_kappa * f_MaxCRateCharge * v_MaxPower[ t ] );
   intake_outtake_bounds_Const[ 1 ][ t ].set_variable( &v_outtake_level[ t ] );
  }

  add_static_constraint( intake_outtake_bounds_Const , "IntakeOuttake_Battery" );

  // Active power bounds constraints

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

 } else {

  // Intake outtake upper bounds design constraints

  intake_outtake_upper_bounds_design_Const.resize(
   boost::multi_array< FRowConstraint , 2 >::extent_gen()
    // 3 dims, i.e., intake, outtake and intake + outtake upper bounds
   [ 3 ][ f_time_horizon ] );

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {

   // Upper bound of the intake level design constraints:
   //
   //      v_intake_level <= - v_MinPower x   x \in {0,1}, for all t
   // => v_intake_level + v_MinPower x <= 0   x \in {0,1}, for all t

   // set the maximum dispatch of converter not to exceed the C-rate of the
   // battery in discharge
   vars.push_back( std::make_pair( &v_intake_level[ t ] , 1.0 ) );
   vars.push_back( std::make_pair( &batt_design ,
                                   f_kappa * f_MaxCRateDischarge *
                                   v_MinPower[ t ] ) );

   intake_outtake_upper_bounds_design_Const[ 0 ][ t ].set_lhs( -Inf< double >() );
   intake_outtake_upper_bounds_design_Const[ 0 ][ t ].set_rhs( 0.0 );
   intake_outtake_upper_bounds_design_Const[ 0 ][ t ].set_function(
    new LinearFunction( std::move( vars ) ) );

   // Upper bound of the outtake level design constraints:
   //
   //      v_outtake_level <= v_MaxPower x     x \in {0,1}, for all t
   // => v_outtake_level - v_MaxPower x <= 0   x \in {0,1}, for all t

   // set the maximum dispatch of converter not to exceed the C-rate of the
   // battery in charge
   vars.push_back( std::make_pair( &v_outtake_level[ t ] , 1.0 ) );
   vars.push_back( std::make_pair( &batt_design ,
                                   -f_kappa * f_MaxCRateCharge *
                                   v_MaxPower[ t ] ) );

   intake_outtake_upper_bounds_design_Const[ 1 ][ t ].set_lhs( -Inf< double >() );
   intake_outtake_upper_bounds_design_Const[ 1 ][ t ].set_rhs( 0.0 );
   intake_outtake_upper_bounds_design_Const[ 1 ][ t ].set_function(
    new LinearFunction( std::move( vars ) ) );

   // Upper bound of the intake + outtake level design constraints:
   //
   //      v_intake_level + v_outtake_level <= v_ConvMaxPower x
   // => v_intake_level + v_outtake_level - v_ConvMaxPower x <= 0
   //                 x \in {0,1}, for all t

   vars.push_back( std::make_pair( &v_intake_level[ t ] , 1.0 ) );
   vars.push_back( std::make_pair( &v_outtake_level[ t ] , 1.0 ) );
   vars.push_back( std::make_pair( &conv_design ,
                                   -f_kappa * v_ConvMaxPower[ t ] ) );

   intake_outtake_upper_bounds_design_Const[ 2 ][ t ].set_lhs( -Inf< double >() );
   intake_outtake_upper_bounds_design_Const[ 2 ][ t ].set_rhs( 0.0 );
   intake_outtake_upper_bounds_design_Const[ 2 ][ t ].set_function(
    new LinearFunction( std::move( vars ) ) );
  }

  add_static_constraint( intake_outtake_upper_bounds_design_Const ,
                         "IntakeOuttake_Design_Battery" );

  // Active power bounds design constraints

  active_power_bounds_design_Const.resize(
   boost::multi_array< FRowConstraint , 2 >::extent_gen()
   [ 2 ][ f_time_horizon ] );  // 2 dims, i.e., the lower and upper bounds

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {

   // Lower bound of the active power design constraints:
   //
   //      v_minimum_power x <= v_active_power     x \in {0,1}, for all t
   // => 0 <= v_active_power - v_minimum_power x   x \in {0,1}, for all t

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
   //      v_active_power <= v_maximum_power x     x \in {0,1}, for all t
   // => v_active_power - v_maximum_power x <= 0   x \in {0,1}, for all t

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

   vars.push_back( std::make_pair( &batt_design ,
                                   -f_kappa * v_MaxPower[ t ] ) );

   active_power_bounds_design_Const[ 1 ][ t ].set_lhs( -Inf< double >() );
   active_power_bounds_design_Const[ 1 ][ t ].set_rhs( 0.0 );
   active_power_bounds_design_Const[ 1 ][ t ].set_function(
    new LinearFunction( std::move( vars ) ) );
  }

  add_static_constraint( active_power_bounds_design_Const ,
                         "ActivePower_Design_Battery" );
 }

 // Initializing power_intake_outtake_Const

 power_intake_outtake_Const.resize( f_time_horizon );

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {

  vars.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );
  vars.push_back( std::make_pair( &v_intake_level[ t ] , -1.0 ) );
  vars.push_back( std::make_pair( &v_outtake_level[ t ] , 1.0 ) );

  power_intake_outtake_Const[ t ].set_both( 0.0 );
  power_intake_outtake_Const[ t ].set_function(
   new LinearFunction( std::move( vars ) ) );
 }

 add_static_constraint( power_intake_outtake_Const ,
                        "PowerIntakeOuttake_Battery" );

 // Ramp-up constraints

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
 }

 add_static_constraint( ramp_up_Const , "RampUp_Battery" );

 // Ramp-down constraints

 if( ! v_DeltaRampDown.empty() ) {

  ramp_down_Const.resize( f_time_horizon );

  vars.push_back( std::make_pair( &v_active_power[ 0 ] , 1.0 ) );

  ramp_down_Const[ 0 ].set_lhs( -v_DeltaRampDown[ 0 ] + f_InitialPower );
  ramp_down_Const[ 0 ].set_rhs( Inf< double >() );
  ramp_down_Const[ 0 ].set_function( new LinearFunction( std::move( vars ) ) );

  for( Index t = 1 ; t < f_time_horizon ; ++t ) {

   vars.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );
   vars.push_back( std::make_pair( &v_active_power[ t - 1 ] , -1.0 ) );

   ramp_down_Const[ t ].set_lhs( -v_DeltaRampDown[ 0 ] );
   ramp_down_Const[ t ].set_rhs( Inf< double >() );
   ramp_down_Const[ t ].set_function( new LinearFunction( std::move( vars ) ) );
  }
 }

 add_static_constraint( ramp_down_Const , "RampDown_Battery" );

 // Initializing demand_Const

 demand_Const.resize( f_time_horizon );

 vars.push_back( std::make_pair( &v_storage_level[ 0 ] , 1.0 ) );

 if( f_InitialStorage < 0 )  // cyclical notation
  vars.push_back( std::make_pair( &v_storage_level[ f_time_horizon - 1 ] ,
                                  -1.0 ) );

 double outtake_coeff = -1;
 if( ! v_StoringBatteryRho.empty() )
  outtake_coeff = -v_StoringBatteryRho[ 0 ];
 vars.push_back( std::make_pair( &v_outtake_level[ 0 ] , outtake_coeff ) );

 double intake_coeff = 1;
 if( ! v_ExtractingBatteryRho.empty() )
  intake_coeff = v_ExtractingBatteryRho[ 0 ];
 vars.push_back( std::make_pair( &v_intake_level[ 0 ] , intake_coeff ) );

 if( ! v_Demand.empty() )
  demand_Const[ 0 ].set_both(
   ( f_InitialStorage < 0 ? 0.0 : f_InitialStorage ) - v_Demand[ 0 ] );
 else
  demand_Const[ 0 ].set_both(
   ( f_InitialStorage < 0 ? 0.0 : f_InitialStorage ) );

 demand_Const[ 0 ].set_function( new LinearFunction( std::move( vars ) ) );

 for( Index t = 1 ; t < f_time_horizon ; ++t ) {

  double intake_coeff = 1;
  if( ! v_ExtractingBatteryRho.empty() )
   intake_coeff = v_ExtractingBatteryRho[ t ];
  vars.push_back( std::make_pair( &v_intake_level[ t ] , intake_coeff ) );

  double outtake_coeff = -1;
  if( ! v_StoringBatteryRho.empty() )
   outtake_coeff = -v_StoringBatteryRho[ t ];
  vars.push_back( std::make_pair( &v_outtake_level[ t ] , outtake_coeff ) );

  vars.push_back( std::make_pair( &v_storage_level[ t ] , 1.0 ) );
  vars.push_back( std::make_pair( &v_storage_level[ t - 1 ] , -1.0 ) );

  if( ! v_Demand.empty() )
   demand_Const[ t ].set_both( -v_Demand[ t ] );
  else
   demand_Const[ t ].set_both( 0.0 );

  demand_Const[ t ].set_function( new LinearFunction( std::move( vars ) ) );
 }

 add_static_constraint( demand_Const , "Demand_Battery" );

 if( f_BattInvestmentCost == 0 ) {

  // Storage level bound constraints

  storage_level_bounds_Const.resize( f_time_horizon );

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   storage_level_bounds_Const[ t ].set_lhs( f_kappa * v_MinStorage[ t ] );
   storage_level_bounds_Const[ t ].set_rhs( f_kappa * v_MaxStorage[ t ] );
   storage_level_bounds_Const[ t ].set_variable( &v_storage_level[ t ] );
  }

  add_static_constraint( storage_level_bounds_Const , "StorageLevel_Battery" );

 } else {

  // Storage level bound design constraints

  storage_level_bounds_design_Const.resize(
   boost::multi_array< FRowConstraint , 2 >::extent_gen()
   [ 2 ][ f_time_horizon ] );  // 2 dims, i.e., the lower and upper bounds

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {

   // Lower bound of the storage level design constraints:
   //
   //      v_MinStorage x <= v_storage_level     x \in {0,1}, for all t
   // => 0 <= v_storage_level - v_MinStorage x   x \in {0,1}, for all t

   vars.push_back( std::make_pair( &v_storage_level[ t ] , 1.0 ) );
   vars.push_back( std::make_pair( &batt_design ,
                                   -f_kappa * v_MinStorage[ t ] ) );

   storage_level_bounds_design_Const[ 0 ][ t ].set_lhs( 0.0 );
   storage_level_bounds_design_Const[ 0 ][ t ].set_rhs( Inf< double >() );
   storage_level_bounds_design_Const[ 0 ][ t ].set_function(
    new LinearFunction( std::move( vars ) ) );

   // Upper bound of the storage level design constraints:
   //
   //      v_storage_level <= v_MaxPower x     x \in {0,1}, for all t
   // => v_storage_level - v_MaxPower x <= 0   x \in {0,1}, for all t

   vars.push_back( std::make_pair( &v_storage_level[ t ] , 1.0 ) );
   vars.push_back( std::make_pair( &batt_design ,
                                   -f_kappa * v_MaxStorage[ t ] ) );

   storage_level_bounds_design_Const[ 1 ][ t ].set_lhs( -Inf< double >() );
   storage_level_bounds_design_Const[ 1 ][ t ].set_rhs( 0.0 );
   storage_level_bounds_design_Const[ 1 ][ t ].set_function(
    new LinearFunction( std::move( vars ) ) );
  }

  add_static_constraint( storage_level_bounds_design_Const ,
                         "StorageLevel_Design_Battery" );
 }

 // Initializing intake_outtake_binary_Const

 if( ! v_battery_binary.empty() ) {

  intake_outtake_binary_Const.resize(
   boost::multi_array< FRowConstraint , 2 >::extent_gen()
   [ 2 ][ f_time_horizon ] );  // 2 dims, i.e., the intake and outtake bounds

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {

   // v_intake_level <= v_MaxPower b   b \in {0,1}, for all t

   vars.push_back( std::make_pair( &v_intake_level[ t ] , 1.0 ) );
   vars.push_back( std::make_pair( &v_battery_binary[ t ] ,
                                   -f_kappa * v_MaxPower[ t ] ) );

   intake_outtake_binary_Const[ 0 ][ t ].set_lhs( -Inf< double >() );
   intake_outtake_binary_Const[ 0 ][ t ].set_rhs( 0.0 );
   intake_outtake_binary_Const[ 0 ][ t ].set_function(
    new LinearFunction( std::move( vars ) ) );

   //      v_outtake_level <= - v_MinPower ( 1 - b )
   // => v_outtake_level - v_MinPower b <= - v_MinPower
   //             b \in [0,1], for all t

   vars.push_back( std::make_pair( &v_outtake_level[ t ] , 1.0 ) );
   vars.push_back( std::make_pair( &v_battery_binary[ t ] ,
                                   -f_kappa * v_MinPower[ t ] ) );

   intake_outtake_binary_Const[ 1 ][ t ].set_lhs( -Inf< double >() );
   intake_outtake_binary_Const[ 1 ][ t ].set_rhs( -f_kappa * v_MinPower[ t ] );
   intake_outtake_binary_Const[ 1 ][ t ].set_function(
    new LinearFunction( std::move( vars ) ) );
  }

  add_static_constraint( intake_outtake_binary_Const ,
                         "Intake_Outtake_Binary_Battery" );
 }

 // Initializing primary_upper_bound_Const

 if( reserve_vars & 1u )  // if UCBlock has primary demand variables
  if( ! v_MaxPrimaryPower.empty() ) {
   // if this unit produces any primary reserve

   primary_upper_bound_Const.resize( f_time_horizon );

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {

    primary_upper_bound_Const[ t ].set_rhs(
     f_kappa * v_MaxPrimaryPower[ t ] );
    primary_upper_bound_Const[ t ].set_variable(
     &v_primary_spinning_reserve[ t ] );
   }

   add_static_constraint( primary_upper_bound_Const ,
                          "Primary_UpperBound_Battery" );
  }

 // Initializing secondary_upper_bound_Const

 if( reserve_vars & 2u )  // if UCBlock has secondary demand variables
  if( ! v_MaxSecondaryPower.empty() ) {
   // if this unit produces any secondary reserve

   secondary_upper_bound_Const.resize( f_time_horizon );

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {

    secondary_upper_bound_Const[ t ].set_rhs(
     f_kappa * v_MaxSecondaryPower[ t ] );
    secondary_upper_bound_Const[ t ].set_variable(
     &v_secondary_spinning_reserve[ t ] );
   }

   add_static_constraint( secondary_upper_bound_Const ,
                          "Secondary_UpperBound_Battery" );
  }

/*------------------------------ ZOConstraint ------------------------------*/

 if( ! v_battery_binary.empty() )

  if( generate_ZOConstraints ) {

   // the battery binary bound constraints
   battery_binary_bound_Const.resize( f_time_horizon );

   for( Index t = 0 ; t < f_time_horizon ; ++t )
    battery_binary_bound_Const[ t ].set_variable( &v_battery_binary[ t ] );

   add_static_constraint( battery_binary_bound_Const , "Binary_Battery" );
  }

 set_constraints_generated();

}  // end( BatteryUnitBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::generate_objective( Configuration *objc )
{
 if( objective_generated() )  // Objective has already been generated
  return;                     // nothing to do

 auto lf = new LinearFunction();

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {
  lf->add_variable( &v_intake_level[ t ] , f_scale * v_Cost[ t ] , eDryRun );
  lf->add_variable( &v_outtake_level[ t ] , f_scale * v_Cost[ t ] , eDryRun );
 }

 if( f_BattInvestmentCost != 0 )
  lf->add_variable( &batt_design , f_BattInvestmentCost );

 if( f_ConvInvestmentCost != 0 )
  lf->add_variable( &conv_design , f_ConvInvestmentCost );

 objective.set_function( lf );
 objective.set_sense( Objective::eMin );

 // Set Block objective
 this->set_objective( &objective );

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
  && ColVariable::is_feasible( v_primary_spinning_reserve , tol )
  && ColVariable::is_feasible( v_secondary_spinning_reserve , tol )
  // Constraints: notice that the ZOConstraint are not checked, since the
  // corresponding check is made on the ColVariable
  && RowConstraint::is_feasible( active_power_bounds_Const , tol , rel_viol )
  && RowConstraint::is_feasible( intake_outtake_upper_bounds_design_Const , tol , rel_viol )
  && RowConstraint::is_feasible( storage_level_bounds_design_Const , tol , rel_viol )
  && RowConstraint::is_feasible( intake_outtake_binary_Const , tol , rel_viol )
  && RowConstraint::is_feasible( power_intake_outtake_Const , tol , rel_viol )
  && RowConstraint::is_feasible( ramp_up_Const , tol , rel_viol )
  && RowConstraint::is_feasible( ramp_down_Const , tol , rel_viol )
  && RowConstraint::is_feasible( demand_Const , tol , rel_viol )
  && RowConstraint::is_feasible( storage_level_bounds_Const , tol , rel_viol )
  && RowConstraint::is_feasible( intake_outtake_bounds_Const , tol , rel_viol )
  && RowConstraint::is_feasible( primary_upper_bound_Const , tol , rel_viol )
  && RowConstraint::is_feasible( secondary_upper_bound_Const , tol , rel_viol ) );

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

 if( f_BattInvestmentCost != 0 )
  ::serialize( group , "BatteryInvestmentCost" , netCDF::NcDouble() ,
               f_BattInvestmentCost );

 if( f_BattMaxCapacity != 0 )
  ::serialize( group , "BatteryMaxCapacity" , netCDF::NcDouble() ,
               f_BattMaxCapacity );

 if( f_ConvInvestmentCost != 0 )
  ::serialize( group , "ConverterInvestmentCost" , netCDF::NcDouble() ,
               f_ConvInvestmentCost );

 if( f_ConvMaxCapacity != 0 )
  ::serialize( group , "ConverterMaxCapacity" , netCDF::NcDouble() ,
               f_ConvMaxCapacity );

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
 serialize( "MaxPrimaryPower" , v_MaxPrimaryPower );
 serialize( "MaxSecondaryPower" , v_MaxSecondaryPower );
 serialize( "DeltaRampUp" , v_DeltaRampUp );
 serialize( "DeltaRampDown" , v_DeltaRampDown );
 serialize( "StoringBatteryRho" , v_StoringBatteryRho );
 serialize( "ExtractingBatteryRho" , v_ExtractingBatteryRho );
 serialize( "Cost" , v_Cost );
 serialize( "Demand" , v_Demand );

}  // end( BatteryUnitBlock::serialize )

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
 rng.second = std::min( rng.second , decltype( rng.second )( 1 ) );
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
 rng.second = std::min( rng.second , decltype( rng.second )( 1 ) );
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
 if( ! intake_outtake_bounds_Const.empty() )

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   intake_outtake_bounds_Const[ 0 ][ t ].set_rhs(
    -f_kappa * f_MaxCRateDischarge * v_MinPower[ t ] , issueAMod );
   intake_outtake_bounds_Const[ 1 ][ t ].set_rhs(
    f_kappa * f_MaxCRateCharge * v_MaxPower[ t ] , issueAMod );
  }

 else if( ! intake_outtake_upper_bounds_design_Const.empty() )

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {

   auto f0 = static_cast< LinearFunction * >(
    intake_outtake_upper_bounds_design_Const[ 0 ][ t ].get_function() );

   const auto batt_design_idx0 = f0->is_active( &batt_design );

   if( batt_design_idx0 == Inf< Index >() )
    throw( std::logic_error( "BatteryUnitBlock::update_kappa_in_cnstrs: "
                             "expected Variable not found in "
                             "intake_outtake_upper_bounds_design_Const." ) );

   f0->modify_coefficient( batt_design_idx0 ,
                           f_kappa * f_MaxCRateDischarge * v_MinPower[ t ] ,
                           issueAMod );

   auto f1 = static_cast< LinearFunction * >(
    intake_outtake_upper_bounds_design_Const[ 1 ][ t ].get_function() );

   const auto batt_design_idx1 = f1->is_active( &batt_design );

   if( batt_design_idx1 == Inf< Index >() )
    throw( std::logic_error( "BatteryUnitBlock::update_kappa_in_cnstrs: "
                             "expected Variable not found in "
                             "intake_outtake_upper_bounds_design_Const." ) );

   f1->modify_coefficient( batt_design_idx1 ,
                           -f_kappa * f_MaxCRateCharge * v_MaxPower[ t ] ,
                           issueAMod );

   auto f2 = static_cast< LinearFunction * >(
    intake_outtake_upper_bounds_design_Const[ 2 ][ t ].get_function() );

   const auto conv_design_idx2 = f2->is_active( &conv_design );

   if( conv_design_idx2 == Inf< Index >() )
    throw( std::logic_error( "BatteryUnitBlock::update_kappa_in_cnstrs: "
                              "expected Variable not found in "
                              "intake_outtake_upper_bounds_design_Const." ) );

   f2->modify_coefficient( conv_design_idx2 ,
                           -f_kappa * v_ConvMaxPower[ t ] ,
                           issueAMod );
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

   f->modify_coefficient( index , -f_kappa * v_MaxPower[ t ] , issueAMod );
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

   f->modify_coefficient( index , -f_kappa * v_MinPower[ t ] , issueAMod );

   intake_outtake_binary_Const[ 1 ][ t ].set_rhs(
    -f_kappa * v_MinPower[ t ] , issueAMod );
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

void BatteryUnitBlock::update_objective( c_ModParam issueAMod ) {

 if( ! objective_generated() )
  return;  // the Objective has not been generated: nothing to be done

 auto function = static_cast< LinearFunction * >( objective.get_function() );

 LinearFunction::Vec_FunctionValue coefficients;
 coefficients.reserve( 2 * f_time_horizon );

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {
  coefficients.push_back( f_scale * v_Cost[ t ] );
  coefficients.push_back( f_scale * v_Cost[ t ] );
 }

 function->modify_coefficients( std::move( coefficients ) ,
                                Range( 0 , Inf< Index >() ) , issueAMod );

}  // end( BatteryUnitBlock::update_objective )

/*--------------------------------------------------------------------------*/
/*----------------- End File BatteryUnitBlock.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
