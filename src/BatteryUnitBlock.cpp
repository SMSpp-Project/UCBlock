/*--------------------------------------------------------------------------*/
/*----------------------- File BatteryUnitBlock.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the BatteryStorageUnitBlock class.
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
 * \copyright &copy; by Antonio Frangioni, Ali Ghezelsoflu,
 *                      Rafael Durbano Lobato
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

BatteryUnitBlock::~BatteryUnitBlock() {

 Constraint::clear( active_power_bounds_Const );
 Constraint::clear( active_power_bounds_design_Const );
 Constraint::clear( ramp_up_Const );
 Constraint::clear( ramp_down_Const );
 Constraint::clear( power_intake_outtake_Const );
 Constraint::clear( storage_intake_outtake_Const );
 Constraint::clear( intake_binary_Const );
 Constraint::clear( outtake_binary_Const );
 Constraint::clear( demand_Const );

 Constraint::clear( storage_level_bounds_Const );
 Constraint::clear( intake_upper_bound_Const );
 Constraint::clear( outtake_upper_bound_Const );
 Constraint::clear( primary_upper_bound_Const );
 Constraint::clear( secondary_upper_bound_Const );

 Constraint::clear( battery_binary_bound_Const );

 objective.clear();
}

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::deserialize( const netCDF::NcGroup & group ) {

#ifndef NDEBUG
 std::vector< std::string > expected_dims =
  { "TimeHorizon" , "NumberIntervals" };
 check_dimensions( group, expected_dims, std::cerr );

 std::vector< std::string > expected_vars =
  { "MinStorage" , "MaxStorage" , "MinPower" , "MaxPower" , "InitialPower" ,
    "MaxPrimaryPower" , "MaxSecondaryPower" , "DeltaRampUp" , "DeltaRampDown" ,
    "StoringBatteryRho" , "ExtractingBatteryRho" , "InitialStorage" , "Cost" ,
    "Demand" , "Kappa" , "BatteryInvestmentCost" , "ConverterInvestmentCost" };
 check_variables( group, expected_vars, std::cerr );
#endif

 // Deserialize data from the base class
 UnitBlock::deserialize( group );

 // Mandatory variables

 ::deserialize( group , "MinStorage" , v_minimum_storage , false );
 ::deserialize( group , "MaxStorage" , v_maximum_storage , false );
 ::deserialize( group , "MaxPower" , v_maximum_power , false );

 // Optional variables

 if( ! ::deserialize( group , "MinPower" , v_minimum_power ) )
  v_minimum_power.assign( f_time_horizon , 0 );

 ::deserialize( group , f_initial_storage , "InitialStorage" );

 ::deserialize( group , f_initial_power , "InitialPower" );

 ::deserialize( group , f_kappa , "Kappa" );

 ::deserialize( group , "MaxPrimaryPower" , v_maximum_primary_rho );
 ::deserialize( group , "MaxSecondaryPower" , v_maximum_secondary_rho );
 ::deserialize( group , "DeltaRampUp" , v_delta_ramp_up );
 ::deserialize( group , "DeltaRampDown" , v_delta_ramp_down );
 ::deserialize( group , "Demand" , v_demand );
 ::deserialize( group , "StoringBatteryRho" , v_storing_battery_rho );
 ::deserialize( group , "ExtractingBatteryRho" , v_extracting_battery_rho );

 if( ! ::deserialize( group , "Cost" , v_cost ) )
  v_cost.resize( 1 , 0 );

 ::deserialize( group , f_batt_investment_cost , "BatteryInvestmentCost" );
 ::deserialize( group , f_conv_investment_cost , "ConverterInvestmentCost" );

 // Decompress vectors

 decompress_vector( v_minimum_power );
 decompress_vector( v_maximum_power );
 decompress_vector( v_minimum_storage );
 decompress_vector( v_maximum_storage );
 decompress_vector( v_maximum_primary_rho );
 decompress_vector( v_maximum_secondary_rho );
 decompress_vector( v_delta_ramp_up );
 decompress_vector( v_delta_ramp_down );
 decompress_vector( v_storing_battery_rho );
 decompress_vector( v_extracting_battery_rho );
 decompress_vector( v_demand );
 decompress_vector( v_cost );

 check_data_consistency();

}  // end( BatteryUnitBlock::deserialize )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::check_data_consistency( void ) const {

 // Minimum and maximum power

 assert( v_minimum_power.size() == f_time_horizon );
 assert( v_maximum_power.size() == f_time_horizon );

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {
  if( v_minimum_power[ t ] > v_maximum_power[ t ] ) {
   throw( std::logic_error( "BatteryUnitBlock::check_data_consistency: minimum "
                            "power for time " + std::to_string( t ) + " is " +
                            std::to_string( v_minimum_power[ t ] ) + ", which "
                            "greater than the maximum power, which is " +
                            std::to_string( v_maximum_power[ t ] ) + "." ) );
  }
 }

 // Minimum and maximum storage levels

 assert( v_minimum_storage.size() == f_time_horizon );
 assert( v_maximum_storage.size() == f_time_horizon );

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {
  if( ( v_minimum_storage[ t ] > v_maximum_storage[ t ] ) ||
      ( v_minimum_storage[ t ] < 0 ) ) {
   throw( std::logic_error( "BatteryUnitBlock::check_data_consistency: maximum "
                            "and minimum storage levels must be such that "
                            "maximum_storage >= minimum_storage >= 0." ) );
  }
 }

 // Inefficiency of storing and extracting energy

 if( ! v_storing_battery_rho.empty() ) {
  assert( v_storing_battery_rho.size() == f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   if( v_storing_battery_rho[ t ] > 1 ) {
    throw( std::logic_error( "BatteryUnitBlock::check_data_consistency: invalid"
                             " inefficiency of storing energy for time step " +
                             std::to_string( t ) + ": " +
                             std::to_string( v_storing_battery_rho[ t ] ) +
                             ". It must not be greater than 1." ) );
   }
  }
 }

// if( ! v_extracting_battery_rho.empty() ) {
//  assert( v_extracting_battery_rho.size() == f_time_horizon );
//  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
//   if( v_extracting_battery_rho[ t ] < 1 ) {
//    throw( std::logic_error( "BatteryUnitBlock::check_data_consistency: invalid"
//                             " inefficiency of extracting energy for time "
//                             "step " + std::to_string( t ) + ": " +
//                             std::to_string( v_extracting_battery_rho[ t ] ) +
//                             ". It must not be less than 1." ) );
//   }
//  }
// }

 if( ( ! v_storing_battery_rho.empty() ) &&
     ( ! v_extracting_battery_rho.empty() ) ) {
  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   if( v_extracting_battery_rho[ t ] < v_storing_battery_rho[ t ] ) {
    throw( std::logic_error( "BatteryUnitBlock::check_data_consistency: the "
                             "inefficiency of storing energy must not be greater "
                             "than the inefficiency of extracting energy." ) );
   }
  }
 }

 // Delta ramp-up

 if( ! v_delta_ramp_up.empty() ) {
  assert( v_delta_ramp_up.size() == f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   if( v_delta_ramp_up[ t ] < 0 )
    throw( std::invalid_argument( "BatteryUnitBlock::check_data_consistency: "
                                  "wrong DeltaRampUp for time step " +
                                  std::to_string( t ) + ": " +
                                  std::to_string( v_delta_ramp_up[ t ] ) ) );
 }

 // Delta ramp-down

 if( ! v_delta_ramp_down.empty() ) {
  assert( v_delta_ramp_down.size() == f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   if( v_delta_ramp_down[ t ] < 0 )
    throw( std::invalid_argument( "BatteryUnitBlock::check_data_consistency: "
                                  "wrong DeltaRampDown for time step " +
                                  std::to_string( t ) + ": " +
                                  std::to_string( v_delta_ramp_down[ t ] ) ) );
 }

 // Maximum active power that can be used as primary reserve

 if( ! v_maximum_primary_rho.empty() ) {
  assert( v_maximum_primary_rho.size() == f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   if( v_maximum_primary_rho[ t ] < 0 )
    throw( std::invalid_argument( "BatteryUnitBlock::check_data_consistency: "
                                  "the maximum power that can be used as "
                                  "primary reserve for time " +
                                  std::to_string( t ) + " is " +
                                  std::to_string( v_maximum_primary_rho[ t ] ) +
                                  ", but it must be nonnegative." ) );
 }

 // Maximum active power that can be used as secondary reserve

 if( ! v_maximum_secondary_rho.empty() ) {
  assert( v_maximum_secondary_rho.size() == f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   if( v_maximum_secondary_rho[ t ] < 0 )
    throw( std::invalid_argument
           ( "BatteryUnitBlock::check_data_consistency: the maximum power that "
             "can be used as secondary reserve for time " +
             std::to_string( t ) + " is " +
             std::to_string( v_maximum_secondary_rho[ t ] ) +
             ", but it must be nonnegative." ) );
 }

 // Demand

 if( ! v_demand.empty() ) {
  assert( v_demand.size() == f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   if( v_demand[ t ] < 0 )
    throw( std::invalid_argument( "BatteryUnitBlock::check_data_consistency: "
                                  "demand for time " + std::to_string( t ) +
                                  " is " + std::to_string( v_demand[ t ] ) +
                                  ", but is must be nonnegative." ) );
 }

 // Initial storage

 if( f_initial_storage < 0 ) {
  throw( std::invalid_argument( "BatteryUnitBlock::check_data_consistency: "
                                "initial storage is " +
                                std::to_string( f_initial_storage ) +
                                ", but it must be nonnegative." ) );
 }

 // Kappa

 if( f_kappa < 0 ) {
  throw( std::logic_error( "BatteryUnitBlock::check_data_consistency: "
                           "kappa must be nonnegative, but it is" +
                           std::to_string( f_kappa ) + "." ) );
 }

}  // end( BatteryUnitBlock::check_data_consistency )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::generate_abstract_variables( Configuration * stvv ) {

 auto battery_type = get_battery_type();

 if( variables_generated() )
  return; // variables have already been generated

 UnitBlock::generate_abstract_variables( stvv );

 v_storage_level.resize( f_time_horizon );
 for( auto & var : v_storage_level )
  var.set_type( ColVariable::kNonNegative );
 add_static_variable( v_storage_level , "SL_battery" );

 v_intake_level.resize( f_time_horizon );
 for( auto & var : v_intake_level )
  var.set_type( ColVariable::kNonNegative );
 add_static_variable( v_intake_level , "IL_battery" );

 v_outtake_level.resize( f_time_horizon );
 for( auto & var : v_outtake_level )
  var.set_type( ColVariable::kNonNegative );
 add_static_variable( v_outtake_level , "OL_battery" );

 int relax_binary = 0;
 auto config = dynamic_cast<SimpleConfiguration< int > *>( stvv );
 if( ( ! config ) && f_BlockConfig &&
     f_BlockConfig->f_static_variables_Configuration )
  config = dynamic_cast< SimpleConfiguration< int > * >
  ( f_BlockConfig->f_static_variables_Configuration );
 if( config )
  relax_binary = config->f_value;

 v_battery_binary.resize( f_time_horizon );
 for( auto & var : v_battery_binary ) {
  if( relax_binary )
   var.set_type( ColVariable::kPosUnitary );
  else
   var.set_type( ColVariable::kBinary );
 }
 if( battery_type == Binary_Variables_Constraints )
  add_static_variable( v_battery_binary , "BB_battery" );

 // Battery Design Variable
 if( f_batt_investment_cost != 0 ) {
  if( relax_binary )
   v_batt_design.set_type( ColVariable::kPosUnitary );
  else
   v_batt_design.set_type( ColVariable::kBinary );
  add_static_variable( v_batt_design , "D_battery" );
 }

 // Converter Design Variable
 if( f_conv_investment_cost != 0 ) {
  if( relax_binary )
   v_conv_design.set_type( ColVariable::kPosUnitary );
  else
   v_conv_design.set_type( ColVariable::kBinary );
  add_static_variable( v_conv_design , "D_converter" );
 }

 // Active Power Variable
 v_active_power.resize( f_time_horizon );
 for( auto & var : v_active_power )
  var.set_type( ColVariable::kContinuous );
 add_static_variable( v_active_power , "p_battery" );

 // Primary Spinning Reserve Variable
 if( reserve_vars & 1u ) { // if UCBlock has primary demand variables
  // if unit produces any primary reserve
  if( ! v_maximum_primary_rho.empty() ) {
   v_primary_spinning_reserve.resize( f_time_horizon );
   for( auto & var : v_primary_spinning_reserve )
    var.set_type( ColVariable::kNonNegative );
   add_static_variable( v_primary_spinning_reserve , "pr_battery" );
  }
 }

 // Secondary Spinning Reserve Variable
 if( reserve_vars & 2u ) { // if UCBlock has secondary demand variables
  // if unit produces any secondary reserve
  if( ! v_maximum_secondary_rho.empty() ) {
   v_secondary_spinning_reserve.resize( f_time_horizon );
   for( auto & var : v_secondary_spinning_reserve )
    var.set_type( ColVariable::kNonNegative );
   add_static_variable( v_secondary_spinning_reserve , "sc_battery" );
  }
 }

 set_variables_generated();
}  // end( BatteryUnitBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::generate_abstract_constraints( Configuration * stcc ) {

 const auto battery_type = get_battery_type();

 if( constraints_generated() )
  return; // constraints have already been generated

 if( ! variables_generated() )
  throw( std::logic_error( "BatteryUnitBlock::generate_abstract_constraints: "
                           "variables need be generated for constraints "
                           "to be." ) );

 int generate_ZOConstraint = 0;
 auto config = dynamic_cast<SimpleConfiguration< int > *>( stcc );
 if( ( ! config ) && f_BlockConfig &&
     f_BlockConfig->f_static_constraints_Configuration )
  config = dynamic_cast< SimpleConfiguration< int > * >
  ( f_BlockConfig->f_static_constraints_Configuration );
 if( config )
  generate_ZOConstraint = config->f_value;

 // Active power bound constraints

 active_power_bounds_Const.resize(
  boost::multi_array< FRowConstraint , 2 >::extent_gen()
  [ f_time_horizon ][ 2 ] ); // 2 dims, i.e., the lower and upper bounds

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {

  LinearFunction::v_coeff_pair lower_vars;

  lower_vars.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );

  if( reserve_vars & 1u ) { // if UCBlock has primary demand variables
   if( ! v_maximum_primary_rho.empty() ) {
    // if this unit produces any primary reserve
    lower_vars.push_back( std::make_pair( &v_primary_spinning_reserve[ t ] ,
                                          -1.0 ) );
   }
  }

  if( reserve_vars & 2u ) { // if UCBlock has secondary demand variables
   if( ! v_maximum_secondary_rho.empty() ) {
    // if unit produces any secondary reserve
    lower_vars.push_back( std::make_pair( &v_secondary_spinning_reserve[ t ] ,
                                          -1.0 ) );
   }
  }

  active_power_bounds_Const[ t ][ 0 ].set_lhs( f_kappa * v_minimum_power[ t ] );
  active_power_bounds_Const[ t ][ 0 ].set_rhs( Inf< double >() );
  active_power_bounds_Const[ t ][ 0 ].set_function(
   new LinearFunction( std::move( lower_vars ) ) );

  LinearFunction::v_coeff_pair upper_vars;

  upper_vars.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );

  if( reserve_vars & 1u ) // if UCBlock has primary demand variables
   if( ! v_maximum_primary_rho.empty() )
    // if this unit produces any primary reserve
    upper_vars.push_back( std::make_pair( &v_primary_spinning_reserve[ t ] ,
                                          1.0 ) );

  if( reserve_vars & 2u ) // if UCBlock has secondary demand variable
   if( ! v_maximum_secondary_rho.empty() )
    // if this unit produces any secondary reserve
    upper_vars.push_back( std::make_pair( &v_secondary_spinning_reserve[ t ] ,
                                          1.0 ) );

  active_power_bounds_Const[ t ][ 1 ].set_lhs( -Inf< double >() );
  active_power_bounds_Const[ t ][ 1 ].set_rhs( f_kappa * v_maximum_power[ t ] );
  active_power_bounds_Const[ t ][ 1 ].set_function(
   new LinearFunction( std::move( upper_vars ) ) );
 }

 add_static_constraint( active_power_bounds_Const ,
                        "ActivePower_Bounds_Battery" );

 // Active power bound design constraints

 if( f_batt_investment_cost != 0 ) {

  active_power_bounds_design_Const.resize(
   boost::multi_array< FRowConstraint , 2 >::extent_gen()
   [ f_time_horizon ][ 2 ] ); // 2 dims, i.e., the lower and upper bounds

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {

   // Lower bound of the active power design constraints:
   //
   //      v_minimum_power z <= v_active_power     z \in {0,1}, for all t
   // => 0 <= v_active_power - v_minimum_power z   z \in {0,1}, for all t

   LinearFunction::v_coeff_pair lower_vars;

   lower_vars.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );
   lower_vars.push_back( std::make_pair( &v_batt_design ,
                                         -f_kappa * v_minimum_power[ t ] ) );

   active_power_bounds_design_Const[ t ][ 0 ].set_lhs( 0.0 );
   active_power_bounds_design_Const[ t ][ 0 ].set_rhs( Inf< double >() );
   active_power_bounds_design_Const[ t ][ 0 ].set_function(
    new LinearFunction( std::move( lower_vars ) ) );

   // Upper bound of the active power design constraints:
   //
   //      v_active_power <= v_maximum_power z     z \in {0,1}, for all t
   // => v_active_power - v_maximum_power z <= 0   z \in {0,1}, for all t

   LinearFunction::v_coeff_pair upper_vars;

   upper_vars.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );
   upper_vars.push_back( std::make_pair( &v_batt_design ,
                                         -f_kappa * v_maximum_power[ t ] ) );

   active_power_bounds_design_Const[ t ][ 1 ].set_lhs( -Inf< double >() );
   active_power_bounds_design_Const[ t ][ 1 ].set_rhs( 0.0 );
   active_power_bounds_design_Const[ t ][ 1 ].set_function(
    new LinearFunction( std::move( upper_vars ) ) );
  }

  add_static_constraint( active_power_bounds_design_Const ,
                         "ActivePower_Bounds_Design_Battery" );
 }

 // Initializing ramp-up constraints

 if( ! v_delta_ramp_up.empty() ) {

  ramp_up_Const.resize( f_time_horizon );

  LinearFunction::v_coeff_pair vars_1;

  vars_1.push_back( std::make_pair( &v_active_power[ 0 ] , 1.0 ) );

  ramp_up_Const[ 0 ].set_lhs( -Inf< double >() );
  ramp_up_Const[ 0 ].set_rhs( v_delta_ramp_up[ 0 ] + f_initial_power );
  ramp_up_Const[ 0 ].set_function( new LinearFunction( std::move( vars_1 ) ) );

  for( Index t = 1 ; t < f_time_horizon ; ++t ) {

   LinearFunction::v_coeff_pair vars_2;

   vars_2.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );
   vars_2.push_back( std::make_pair( &v_active_power[ t - 1 ] , -1.0 ) );

   ramp_up_Const[ t ].set_lhs( -Inf< double >() );
   ramp_up_Const[ t ].set_rhs( v_delta_ramp_up[ t ] );
   ramp_up_Const[ t ].set_function( new LinearFunction( std::move( vars_2 ) ) );
  }
 }

 add_static_constraint( ramp_up_Const , "RampUp_Battery" );

 // Initializing ramp-down constraints

 if( ! v_delta_ramp_down.empty() ) {

  ramp_down_Const.resize( f_time_horizon );

  LinearFunction::v_coeff_pair vars_1;

  vars_1.push_back( std::make_pair( &v_active_power[ 0 ] , 1.0 ) );

  ramp_down_Const[ 0 ].set_lhs( -v_delta_ramp_down[ 0 ] +
                                      f_initial_power );
  ramp_down_Const[ 0 ].set_rhs( Inf< double >() );
  ramp_down_Const[ 0 ].set_function(
   new LinearFunction( std::move( vars_1 ) ) );

  for( Index t = 1 ; t < f_time_horizon ; ++t ) {

   LinearFunction::v_coeff_pair vars_2;

   vars_2.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );
   vars_2.push_back( std::make_pair( &v_active_power[ t - 1 ] , -1.0 ) );

   ramp_down_Const[ t ].set_lhs( -v_delta_ramp_down[ 0 ] );
   ramp_down_Const[ t ].set_rhs( Inf< double >() );
   ramp_down_Const[ t ].set_function(
    new LinearFunction( std::move( vars_2 ) ) );
  }
 }

 add_static_constraint( ramp_down_Const ,
                        "RampDown_Battery" );

 // Initializing power_intake_outtake_Const

 power_intake_outtake_Const.resize( f_time_horizon );

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {

  LinearFunction::v_coeff_pair vars;

  vars.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );
  vars.push_back( std::make_pair( &v_intake_level[ t ] , -1.0 ) );
  vars.push_back( std::make_pair( &v_outtake_level[ t ] , 1.0 ) );

  power_intake_outtake_Const[ t ].set_both( 0.0 );
  power_intake_outtake_Const[ t ].set_function(
   new LinearFunction( std::move( vars ) ) );
 }

 add_static_constraint( power_intake_outtake_Const ,
                        "Power_Intake_Outtake_Battery" );

 // Initializing intake_upper_bound_Const

 intake_upper_bound_Const.resize( f_time_horizon );

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {
  intake_upper_bound_Const[ t ].set_lhs( 0.0 );
  intake_upper_bound_Const[ t ].set_rhs( f_kappa * v_maximum_power[ t ] );
  intake_upper_bound_Const[ t ].set_variable( &v_intake_level[ t ] );
 }

 add_static_constraint( intake_upper_bound_Const ,
                        "Intake_UpperBound_Battery" );

 // Initializing outtake_upper_bound_Const

 outtake_upper_bound_Const.resize( f_time_horizon );

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {
  outtake_upper_bound_Const[ t ].set_lhs( 0.0 );
  outtake_upper_bound_Const[ t ].set_rhs( f_kappa * v_maximum_power[ t ] );
  outtake_upper_bound_Const[ t ].set_variable( &v_outtake_level[ t ] );
 }

 add_static_constraint( outtake_upper_bound_Const ,
                        "Outtake_UpperBound_Battery" );

 // Initializing demand_Const

 demand_Const.resize( f_time_horizon );

 LinearFunction::v_coeff_pair vars_1;

 vars_1.push_back( std::make_pair( &v_storage_level[ 0 ] , 1.0 ) );
 vars_1.push_back( std::make_pair( &v_storage_level[ f_time_horizon - 1 ] ,
                                   -1.0 ) );

 double outtake_coeff = -1;
 if( ! v_storing_battery_rho.empty() )
  outtake_coeff = -v_storing_battery_rho[ 0 ];
 vars_1.push_back( std::make_pair( &v_outtake_level[ 0 ] , outtake_coeff ) );

 double intake_coeff = 1;
 if( ! v_extracting_battery_rho.empty() )
  intake_coeff = v_extracting_battery_rho[ 0 ];
 vars_1.push_back( std::make_pair( &v_intake_level[ 0 ] , intake_coeff ) );

 double rhs = f_initial_storage;
 if( ! v_demand.empty() )
  rhs -= v_demand[ 0 ];
 demand_Const[ 0 ].set_both( rhs );

 demand_Const[ 0 ].set_function( new LinearFunction( std::move( vars_1 ) ) );

 for( Index t = 1 ; t < f_time_horizon ; ++t ) {

  LinearFunction::v_coeff_pair vars_2;

  double intake_coeff = 1;
  if( ! v_extracting_battery_rho.empty() )
   intake_coeff = v_extracting_battery_rho[ t ];
  vars_2.push_back( std::make_pair( &v_intake_level[ t ] , intake_coeff ) );

  double outtake_coeff = -1;
  if( ! v_storing_battery_rho.empty() )
   outtake_coeff = -v_storing_battery_rho[ t ];
  vars_2.push_back( std::make_pair( &v_outtake_level[ t ] , outtake_coeff ) );

  vars_2.push_back( std::make_pair( &v_storage_level[ t ] , 1.0 ) );
  vars_2.push_back( std::make_pair( &v_storage_level[ t - 1 ] , -1.0 ) );

  double rhs = 0;
  if( ! v_demand.empty() )
   rhs = -v_demand[ t ];
  demand_Const[ t ].set_both( rhs );

  demand_Const[ t ].set_function( new LinearFunction( std::move( vars_2 ) ) );
 }

 add_static_constraint( demand_Const , "demand_Battery" );

 // Initializing storage_level_bounds_Const

 storage_level_bounds_Const.resize( f_time_horizon );

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {

  storage_level_bounds_Const[ t ].set_lhs( f_kappa * v_minimum_storage[ t ] );
  storage_level_bounds_Const[ t ].set_rhs( f_kappa * v_maximum_storage[ t ] );
  storage_level_bounds_Const[ t ].set_variable( &v_storage_level[ t ] );
 }

 add_static_constraint( storage_level_bounds_Const ,
                        "StorageLevel_Bounds_Battery" );

 // Initializing intake_binary_Const

 if( battery_type == Binary_Variables_Constraints ) {

  intake_binary_Const.resize( f_time_horizon );

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {

   LinearFunction::v_coeff_pair vars;

   vars.push_back( std::make_pair( &v_intake_level[ t ] , 1.0 ) );
   vars.push_back( std::make_pair( &v_battery_binary[ t ] ,
                                   -f_kappa * v_maximum_power[ t ] ) );

   intake_binary_Const[ t ].set_lhs( -Inf<double>() );
   intake_binary_Const[ t ].set_rhs( 0.0 );
   intake_binary_Const[ t ].set_function(
    new LinearFunction( std::move( vars ) ) );
  }

  add_static_constraint( intake_binary_Const ,
                         "Intake_Binary_Battery" );

  // Initializing outtake_binary_Const

  outtake_binary_Const.resize( f_time_horizon );

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {

   LinearFunction::v_coeff_pair vars;

   vars.push_back( std::make_pair( &v_outtake_level[ t ] , 1.0 ) );
   vars.push_back( std::make_pair( &v_battery_binary[ t ] ,
                                   -f_kappa * v_minimum_power[ t ] ) );

   outtake_binary_Const[ t ].set_lhs( -Inf<double>() );
   outtake_binary_Const[ t ].set_rhs( -f_kappa * v_minimum_power[ t ] );
   outtake_binary_Const[ t ].set_function(
    new LinearFunction( std::move( vars ) ) );
  }

  add_static_constraint( outtake_binary_Const ,
                         "Outtake_Binary_Battery" );
 }

 // Initializing primary_upper_bound_Const

 if( reserve_vars & 1u ) { // if UCBlock has primary demand variables
  if( ! v_maximum_primary_rho.empty() ) {
   // if this unit produces any primary reserve

   primary_upper_bound_Const.resize( f_time_horizon );

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {

    primary_upper_bound_Const[ t ].set_lhs( 0.0 );
    primary_upper_bound_Const[ t ].set_rhs(
     f_kappa * v_maximum_primary_rho[ t ] );
    primary_upper_bound_Const[ t ].set_variable(
     &v_primary_spinning_reserve[ t ] );
   }

   add_static_constraint( primary_upper_bound_Const ,
                          "Primary_UpperBound_Battery" );
  }
 }

 // Initializing secondary_upper_bound_Const

 if( reserve_vars & 2u ) { // if UCBlock has secondary demand variables
  if( ! v_maximum_secondary_rho.empty() ) {
   // if this unit produces any secondary reserve

   secondary_upper_bound_Const.resize( f_time_horizon );

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {

    secondary_upper_bound_Const[ t ].set_lhs( 0.0 );
    secondary_upper_bound_Const[ t ].set_rhs(
     f_kappa * v_maximum_secondary_rho[ t ] );
    secondary_upper_bound_Const[ t ].set_variable(
     &v_secondary_spinning_reserve[ t ] );
   }

   add_static_constraint( secondary_upper_bound_Const ,
                          "Secondary_UpperBound_Battery" );
  }
 }

/*------------------------------ ZOConstraint ------------------------------*/

 if( battery_type == Binary_Variables_Constraints ) {

  if( generate_ZOConstraint ) {

   // the battery binary bound constraints
   battery_binary_bound_Const.resize( f_time_horizon );

   for( Index t = 0 ; t < f_time_horizon ; ++t )
    battery_binary_bound_Const[ t ].set_variable( &v_battery_binary[ t ] );

   add_static_constraint( battery_binary_bound_Const ,
                          "BB_bound_battery" );
  }
 }

 set_constraints_generated();

}  // end( BatteryUnitBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::generate_objective( Configuration *objc ) {

 if( objective_generated() )
  return; // Objective has already been generated

 if( ! variables_generated() )
  throw( std::logic_error( "BatteryUnitBlock::generate_objective: variables "
                           "need be generated for constraints to be." ) );

 if( get_objective() != nullptr )  // an objective is there already
  return;                         // cowardly (and silently) return

 auto linear_function = new LinearFunction();

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {
  linear_function->add_variable( &v_intake_level[ t ] ,
                                 f_scale * v_cost[ t ] , eDryRun );
  linear_function->add_variable( &v_outtake_level[ t ] ,
                                 f_scale * v_cost[ t ] , eDryRun );
 }

 if( f_batt_investment_cost != 0 )
  linear_function->add_variable( &v_batt_design , f_batt_investment_cost );

 if( f_conv_investment_cost != 0 )
  linear_function->add_variable( &v_conv_design , f_conv_investment_cost );

 objective.set_function( linear_function );
 objective.set_sense( Objective::eMin );

 // Set Block objective
 this->set_objective( &objective );

 set_objective_generated();

}  // end( BatteryUnitBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*------- METHODS FOR LOADING, PRINTING & SAVING THE BatteryUnitBlock ------*/
/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::serialize( netCDF::NcGroup & group ) const {

 UnitBlock::serialize( group );

 // Serialize scalar variables.

 ::serialize( group , "InitialPower" , netCDF::NcDouble() , f_initial_power );
 ::serialize( group , "InitialStorage" , netCDF::NcDouble() ,
              f_initial_storage );
 ::serialize( group , "Kappa" , netCDF::NcDouble() , f_kappa );

 // Serialize one-dimensional variables.

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
  else if( data.size() != 1 ) {
   throw( std::logic_error
          ( "BatteryUnitBlock::serialize: invalid dimension for variable " +
            var_name + ": " + std::to_string( data.size() ) + ". Its dimension "
            "must be one of the following: TimeHorizon, NumberIntervals, 1.") );
  }

  ::serialize( group , var_name , ncType , dimension , data ,
               allow_scalar_var );
 };

 serialize( "MinStorage" , v_minimum_storage );
 serialize( "MaxStorage" , v_maximum_storage );
 serialize( "MinPower" , v_minimum_power );
 serialize( "MaxPower" , v_maximum_power );
 serialize( "MaxPrimaryRho" , v_maximum_primary_rho );
 serialize( "MaxSecondaryRho" , v_maximum_secondary_rho );
 serialize( "DeltaRampUp" , v_delta_ramp_up );
 serialize( "DeltaRampDown" , v_delta_ramp_down );
 serialize( "StoringBatteryRho" , v_storing_battery_rho );
 serialize( "ExtractingBatteryRho" , v_extracting_battery_rho );
 serialize( "Cost" , v_cost );
 serialize( "Demand" , v_demand );
}  // end( BatteryUnitBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::update_initial_storage_in_constraints
( c_ModParam issueAMod ) {

 if( demand_Const.empty() )
  return;

 if( ! v_demand.empty() )
  demand_Const[ 0 ].set_both( f_initial_storage - v_demand[ 0 ] ,
                                    issueAMod );
 else
  demand_Const[ 0 ].set_both( f_initial_storage , issueAMod );
}  // end( BatteryUnitBlock::update_initial_storage_in_constraints )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::set_initial_storage
( std::vector< double >::const_iterator it , Block::Subset && subset ,
  const bool ordered , c_ModParam issuePMod , c_ModParam issueAMod ) {

 if( subset.empty() )
  return;

 // Find the last index 0
 auto index_it = std::find( subset.rbegin() , subset.rend() , 0 );

 if( index_it == subset.rend() )
  return; // 0 is not in subset; return

 std::advance( it , std::distance( index_it , subset.rend() ) - 1 );

 if( f_initial_storage == *it )
  return;

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation
  f_initial_storage = *it;

  if( not_dry_run( issueAMod ) && constraints_generated() ) {
   // Change the abstract representation
   update_initial_storage_in_constraints( issueAMod );
  }
 }

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< BatteryUnitBlockMod >
                           ( this , BatteryUnitBlockMod::eSetInitS ) ,
                           Observer::par2chnl( issuePMod ) );
 }
}  // end( BatteryUnitBlock::set_initial_storage )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::set_initial_storage
( std::vector< double >::const_iterator it , Block::Range rng ,
  c_ModParam issuePMod , c_ModParam issueAMod ) {

 rng.second = std::min( rng.second , decltype( rng.second )( 1 ) );
 if( ! ( rng.first <= 0 && 0 < rng.second ) )
  return; // 0 does not belong to the range; return

 std::advance( it , - rng.first );

 if( f_initial_storage == *it )
  return; // nothing changes; return

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation
  f_initial_storage = *it;

  if( not_dry_run( issueAMod ) && constraints_generated() ) {
   // Change the abstract representation
   update_initial_storage_in_constraints( issueAMod );
  }
 }

 if( issue_pmod( issuePMod ) ) {
  Block::add_Modification( std::make_shared< BatteryUnitBlockMod >
                           ( this , BatteryUnitBlockMod::eSetInitS ) ,
                           Observer::par2chnl( issuePMod ) );
 }
}  // end( BatteryUnitBlock::set_initial_storage )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::update_initial_power_in_constraints
( c_ModParam issueAMod ) {
 if( ! ( ramp_up_Const.empty() || v_delta_ramp_up.empty() ) )
  ramp_up_Const[ 0 ].set_rhs( v_delta_ramp_up[ 0 ] + f_initial_power ,
                                    issueAMod );

 if( ! ( ramp_down_Const.empty() || v_delta_ramp_down.empty() ) )
  ramp_down_Const[ 0 ].set_lhs( -v_delta_ramp_down[ 0 ] +
                                      f_initial_power , issueAMod );
}  // end( BatteryUnitBlock::update_initial_power_in_constraints )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::set_initial_power
( std::vector< double >::const_iterator it , Block::Subset && subset ,
  const bool ordered , c_ModParam issuePMod , c_ModParam issueAMod ) {

 if( subset.empty() )
  return;

 // Find the last index 0
 auto index_it = std::find( subset.rbegin() , subset.rend() , 0 );

 if( index_it == subset.rend() )
  return; // 0 is not in subset; return

 std::advance( it , std::distance( index_it , subset.rend() ) - 1 );

 if( f_initial_power == *it )
  return; // nothing changes; return

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation
  f_initial_power = *it;

  if( not_dry_run( issueAMod ) && constraints_generated() ) {
   // Change the abstract representation
   update_initial_power_in_constraints( issueAMod );
  }
 }

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< BatteryUnitBlockMod >
                           ( this , BatteryUnitBlockMod::eSetInitP ) ,
                           Observer::par2chnl( issuePMod ) );
 }

}  // end( BatteryUnitBlock::set_initial_power )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::set_initial_power
( std::vector< double >::const_iterator it , Block::Range rng ,
  c_ModParam issuePMod , c_ModParam issueAMod ) {

 rng.second = std::min( rng.second , decltype( rng.second )( 1 ) );
 if( ! ( rng.first <= 0 && 0 < rng.second ) )
  return; // 0 does not belong to the range; return

 std::advance( it , - rng.first );

 if( f_initial_power == *it )
  return; // nothing changes; return

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation
  f_initial_power = *it;

  if( not_dry_run( issueAMod ) && constraints_generated() ) {
   // Change the abstract representation
   update_initial_power_in_constraints( issueAMod );
  }
 }

 if( issue_pmod( issuePMod ) ) {
  Block::add_Modification( std::make_shared< BatteryUnitBlockMod >
                           ( this , BatteryUnitBlockMod::eSetInitP ) ,
                           Observer::par2chnl( issuePMod ) );
 }
}  // end( BatteryUnitBlock::set_initial_power )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::scale
( std::vector< double >::const_iterator values , Subset && subset ,
  const bool ordered , c_ModParam issuePMod , c_ModParam issueAMod ) {

 if( subset.empty() )
  return; // Since the given Subset is empty, no operation is performed

 if( f_scale == *values )
  return; // The scale factor does not change: nothing to do

 if( not_dry_run( issuePMod ) ) {
  f_scale = *values; // Update the scale factor

  if( not_dry_run( issueAMod ) ) {
   // Update the abstract representation
   if( objective_generated() )
    // Update the Objective
    update_objective( issueAMod );
  }
 }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< UnitBlockMod >
                           ( this , UnitBlockMod::eScale ) ,
                           Observer::par2chnl( issuePMod ) );
 else if( auto f_Block = get_f_Block() )
  f_Block->add_Modification( std::make_shared< UnitBlockMod >
                             ( this , UnitBlockMod::eScale ) ,
                             Observer::par2chnl( issuePMod ) );

}  // end( BatteryUnitBlock::scale )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::update_kappa_in_constraints( ModParam issueAMod ) {

 if( ! active_power_bounds_Const.empty() ) {
  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   active_power_bounds_Const[ t ][ 0 ].set_lhs(
    f_kappa * v_minimum_power[ t ] , issueAMod );
   active_power_bounds_Const[ t ][ 1 ].set_rhs(
    f_kappa * v_maximum_power[ t ] , issueAMod );
  }
 }

 if( ! intake_upper_bound_Const.empty() ) {
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   intake_upper_bound_Const[ t ].set_rhs(
    f_kappa * v_maximum_power[ t ] , issueAMod );
 }

 if( ! outtake_upper_bound_Const.empty() ) {
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   outtake_upper_bound_Const[ t ].set_rhs(
    f_kappa * v_maximum_power[ t ] , issueAMod );
 }

 if( ! storage_level_bounds_Const.empty() ) {
  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   storage_level_bounds_Const[ t ].set_lhs(
    f_kappa * v_minimum_storage[ t ] , issueAMod );
   storage_level_bounds_Const[ t ].set_rhs(
    f_kappa * v_maximum_storage[ t ] , issueAMod );
  }
 }

 if( ! intake_binary_Const.empty() ) {
  for( Index t = 0 ; t < f_time_horizon ; ++t ) {

   auto f = static_cast< LinearFunction * >
    ( intake_binary_Const[ t ].get_function() );

   const auto index = f->is_active( & v_battery_binary[ t ] );

   if( index == Inf< Index >() )
    throw( std::logic_error( "BatteryUnitBlock::set_kappa: expected Variable"
                             "not found in intake_binary_Const." ) );

   f->modify_coefficient( index , -f_kappa * v_maximum_power[ t ] ,
                          issueAMod );
  }
 }

 if( ! outtake_binary_Const.empty() ) {

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {

   auto f = static_cast< LinearFunction * >
    ( outtake_binary_Const[ t ].get_function() );

   const auto index = f->is_active( & v_battery_binary[ t ] );

   if( index == Inf< Index >() )
    throw( std::logic_error( "BatteryUnitBlock::set_kappa: expected Variable"
                             "not found in outtake_binary_Const." ) );

   f->modify_coefficient( index , -f_kappa * v_minimum_power[ t ] ,
                          issueAMod );

   outtake_binary_Const[ t ].set_rhs( - f_kappa * v_minimum_power[ t ] ,
                                            issueAMod );
  }
 }

 if( ! primary_upper_bound_Const.empty() ) {
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   primary_upper_bound_Const[ t ].set_rhs
    ( f_kappa * v_maximum_primary_rho[ t ] , issueAMod );
 }

 if( ! secondary_upper_bound_Const.empty() ) {
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   secondary_upper_bound_Const[ t ].set_rhs
    ( f_kappa * v_maximum_secondary_rho[ t ] , issueAMod );
 }
}

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::set_kappa
( std::vector< double >::const_iterator values , Subset && subset ,
  const bool ordered , ModParam issuePMod , ModParam issueAMod ) {

 if( subset.empty() )
  return; // Since the given Subset is empty, no operation is performed

 if( f_kappa == *values )
  return; // The kappa constant does not change: nothing to do

 if( not_dry_run( issuePMod ) ) {
  f_kappa = *values; // Update the kappa constant

  if( not_dry_run( issueAMod ) ) {
   // Update the abstract representation
   if( constraints_generated() ) {
    // Update the constraints
    update_kappa_in_constraints( issueAMod );
   }  // end( constraints_generated )
  }  // end( if( not_dry_run( issueAMod ) )
 }  // end( if( not_dry_run( issuePMod ) )

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< BatteryUnitBlockMod >
                           ( this , BatteryUnitBlockMod::eSetKappa ) ,
                           Observer::par2chnl( issuePMod ) );
 }
}  // end( BatteryUnitBlock::set_kappa )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::set_kappa
( std::vector< double >::const_iterator values , Range rng ,
  ModParam issuePMod , ModParam issueAMod ) {

 if( rng.first >= rng.second )
  return; // An empty Range was given: no operation is performed.

 Subset subset( 1 , 0 );

 set_kappa( values , std::move( subset ) , true , issuePMod , issueAMod );

}  // end( BatteryUnitBlock::set_kappa )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::update_objective( c_ModParam issueAMod ) {

 if( ! objective_generated() )
  return; // the Objective has not been generated: nothing to be done

 auto function = static_cast< LinearFunction * >( objective.get_function() );

 LinearFunction::Vec_FunctionValue coefficients;
 coefficients.reserve( 2 * f_time_horizon );

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {
  coefficients.push_back( f_scale * v_cost[ t ] );
  coefficients.push_back( f_scale * v_cost[ t ] );
 }

 function->modify_coefficients( std::move( coefficients ) ,
                                Range( 0 , Inf< Index >() ) , issueAMod );

}  // end( BatteryUnitBlock::update_objective )

/*--------------------------------------------------------------------------*/
/*----------------- End File BatteryUnitBlock.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
