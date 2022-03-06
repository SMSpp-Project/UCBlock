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

#include "FRealObjective.h"
#include "LinearFunction.h"
#include "BatteryUnitBlock.h"
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
 auto clear_constraints =
  []( std::vector< FRowConstraint > & constraints ) {
   for( auto & constraint : constraints )
    constraint.clear();
  };

 clear_constraints( active_power_upper_bound_Constraints );
 clear_constraints( active_power_lower_bound_Constraints );
 clear_constraints( ramp_up_Constraints );
 clear_constraints( ramp_down_Constraints );
 clear_constraints( power_intake_outtake_Constraints );
 clear_constraints( storage_intake_outtake_Constraints );
 clear_constraints( intake_binary_Constraints );
 clear_constraints( outtake_binary_Constraints );
 clear_constraints( demand_Constraints );

 auto clear_boxconstraints =
  []( std::vector< BoxConstraint > & constraints ) {
   for( auto & constraint : constraints )
    constraint.clear();
  };

 clear_boxconstraints( intake_upper_bound_Constraints );
 clear_boxconstraints( storage_level_bounds_Constraints );
 clear_boxconstraints( primary_upper_bound_Constraints );
 clear_boxconstraints( secondary_upper_bound_Constraints );

 auto clear_ZOConstraints =
  []( std::vector< ZOConstraint > & constraints ) {
   for( auto & constraint : constraints )
    constraint.clear();
  };

 clear_ZOConstraints( battery_binary_bound_Constraints );

 objective.clear();
}

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::deserialize( const netCDF::NcGroup & group ) {

#ifndef NDEBUG
 std::vector< std::string > expected_dims =
  { "TimeHorizon" , "NumberIntervals" };
 check_dimensions( group , expected_dims , std::cerr );

 std::vector< std::string > expected_vars =
  { "MinStorage" , "MaxStorage" , "MinPower" , "MaxPower" , "InitialPower" ,
    "MaxPrimaryPower" , "MaxSecondaryPower" , "DeltaRampUp" , "DeltaRampDown" ,
    "StoringBatteryRho" , "ExtractingBatteryRho" , "InitialStorage" , "Cost" ,
    "Demand" };
 check_variables( group , expected_vars , std::cerr );
#endif

 // Mandatory variables

 ::deserialize( group , "MinStorage" , v_minimum_storage , false );
 ::deserialize( group , "MaxStorage" , v_maximum_storage , false );
 ::deserialize( group , "MinPower" , v_minimum_power , false );
 ::deserialize( group , "MaxPower" , v_maximum_power , false );
 ::deserialize( group , f_initial_storage , "InitialStorage" , false );

 // Optional variables

 if( !::deserialize( group , f_initial_power , "InitialPower" ) )
  f_initial_power = 0;

 ::deserialize( group , "MaxPrimaryPower" , v_maximum_primary_rho );
 ::deserialize( group , "MaxSecondaryPower" , v_maximum_secondary_rho );
 ::deserialize( group , "DeltaRampUp" , v_delta_ramp_up );
 ::deserialize( group , "DeltaRampDown" , v_delta_ramp_down );
 ::deserialize( group , "Cost" , v_cost );
 ::deserialize( group , "Demand" , v_demand );
 ::deserialize( group , "StoringBatteryRho" , v_storing_battery_rho );
 ::deserialize( group , "ExtractingBatteryRho" , v_extracting_battery_rho );

 // Deserialize data from the base class

 UnitBlock::deserialize( group );

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

 check_data_consistency();

} // end( BatteryUnitBlock::deserialize )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::check_data_consistency() const {

 // Minimum and maximum power

 assert( v_minimum_power.size() == f_time_horizon );
 assert( v_maximum_power.size() == f_time_horizon );

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {
  if( v_minimum_power[ t ] > v_maximum_power[ t ] ) {
   throw ( std::logic_error(
    "BatteryUnitBlock::check_data_consistency: minimum "
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
   throw ( std::logic_error(
    "BatteryUnitBlock::check_data_consistency: maximum "
    "and minimum storage levels must be such that "
    "maximum_storage >= minimum_storage >= 0." ) );
  }
 }

 // Inefficiency of storing and extracting energy

 if( !v_storing_battery_rho.empty() ) {
  assert( v_storing_battery_rho.size() == f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   if( v_storing_battery_rho[ t ] > 1 ) {
    throw ( std::logic_error(
     "BatteryUnitBlock::check_data_consistency: invalid"
     " inefficiency of storing energy for time step " +
     std::to_string( t ) + ": " +
     std::to_string( v_storing_battery_rho[ t ] ) +
     ". It must not be greater than 1." ) );
   }
  }
 }

 if( !v_extracting_battery_rho.empty() ) {
  assert( v_extracting_battery_rho.size() == f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   if( v_extracting_battery_rho[ t ] < 1 ) {
    throw ( std::logic_error(
     "BatteryUnitBlock::check_data_consistency: invalid"
     " inefficiency of extracting energy for time "
     "step " + std::to_string( t ) + ": " +
     std::to_string( v_extracting_battery_rho[ t ] ) +
     ". It must not be less than 1." ) );
   }
  }
 }

 if( ( !v_storing_battery_rho.empty() ) &&
     ( !v_extracting_battery_rho.empty() ) ) {
  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   if( v_extracting_battery_rho[ t ] < v_storing_battery_rho[ t ] ) {
    throw ( std::logic_error(
     "BatteryUnitBlock::check_data_consistency: the ine"
     "fficiency of storing energy must not be greater "
     "than the inneficiency of extracting energy." ) );
   }
  }
 }

 // Delta ramp-up

 if( !v_delta_ramp_up.empty() ) {
  assert( v_delta_ramp_up.size() == f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   if( v_delta_ramp_up[ t ] < 0 )
    throw ( std::invalid_argument( "BatteryUnitBlock::check_data_consistency: "
                                   "wrong DeltaRampUp for time step " +
                                   std::to_string( t ) + ": " +
                                   std::to_string( v_delta_ramp_up[ t ] ) ) );
 }

 // Delta ramp-down

 if( !v_delta_ramp_down.empty() ) {
  assert( v_delta_ramp_down.size() == f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   if( v_delta_ramp_down[ t ] < 0 )
    throw ( std::invalid_argument( "BatteryUnitBlock::check_data_consistency: "
                                   "wrong DeltaRampDown for time step " +
                                   std::to_string( t ) + ": " +
                                   std::to_string( v_delta_ramp_down[ t ] ) ) );
 }

 // Maximum active power that can be used as primary reserve

 if( !v_maximum_primary_rho.empty() ) {
  assert( v_maximum_primary_rho.size() == f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   if( v_maximum_primary_rho[ t ] < 0 )
    throw ( std::invalid_argument( "BatteryUnitBlock::check_data_consistency: "
                                   "the maximum power that can be used as "
                                   "primary reserve for time " +
                                   std::to_string( t ) + " is " +
                                   std::to_string(
                                    v_maximum_primary_rho[ t ] ) +
                                   ", but it must be nonnegative." ) );
 }

 // Maximum active power that can be used as secondary reserve

 if( !v_maximum_secondary_rho.empty() ) {
  assert( v_maximum_secondary_rho.size() == f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   if( v_maximum_secondary_rho[ t ] < 0 )
    throw ( std::invalid_argument
     ( "BatteryUnitBlock::check_data_consistency: the maximum power that "
       "can be used as secondary reserve for time " +
       std::to_string( t ) + " is " +
       std::to_string( v_maximum_secondary_rho[ t ] ) +
       ", but it must be nonnegative." ) );
 }

 // Demand

 if( !v_demand.empty() ) {
  assert( v_demand.size() == f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   if( v_demand[ t ] < 0 )
    throw ( std::invalid_argument( "BatteryUnitBlock::check_data_consistency: "
                                   "demand for time " + std::to_string( t ) +
                                   " is " + std::to_string( v_demand[ t ] ) +
                                   ", but is must be nonnegative." ) );
 }

 // Initial storage

 if( f_initial_storage < 0 ) {
  throw ( std::invalid_argument( "BatteryUnitBlock::check_data_consistency: "
                                 "initial storage is " +
                                 std::to_string( f_initial_storage ) +
                                 ", but it must be nonnegative." ) );
 }

} // end( BatteryUnitBlock::check_data_consistency )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::generate_abstract_variables( Configuration * stvv ) {
 auto battery_type = get_battery_type();

 if( variables_generated() )
  return; // variables have already been generated

 UnitBlock::generate_abstract_variables( stvv );

 int relax_binary = 0;
 auto config = dynamic_cast<SimpleConfiguration< int > *>( stvv );
 if( ( !config ) && f_BlockConfig &&
     f_BlockConfig->f_static_variables_Configuration )
  config = dynamic_cast< SimpleConfiguration< int > * >
  ( f_BlockConfig->f_static_variables_Configuration );
 if( config )
  relax_binary = config->f_value;

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


 v_battery_binary.resize( f_time_horizon );
 for( auto & var : v_battery_binary ) {
  if( relax_binary )
   var.set_type( ColVariable::kPosUnitary );
  else
   var.set_type( ColVariable::kBinary );
 }
 if( battery_type == Binary_Variables_Constraints ) {

  add_static_variable( v_battery_binary , "BB_battery" );
 }
 // Active Power Variable

 v_active_power.resize( f_time_horizon );
 for( auto & var : v_active_power )
  var.set_type( ColVariable::kContinuous );
 add_static_variable( v_active_power , "p_battery" );

 // Primary Spinning Reserve Variable
 if( reserve_vars & 1u ) { // if UCBlock has primary demand variables
  if( !v_maximum_primary_rho.empty() ) { // if unit produces any primary reserve
   v_primary_spinning_reserve.resize( f_time_horizon );
   for( auto & var : v_primary_spinning_reserve ) {
    var.set_type( ColVariable::kNonNegative );
   }
   add_static_variable( v_primary_spinning_reserve , "pr_battery" );
  }
 }

 // Secondary Spinning Reserve Variable
 if( reserve_vars & 2u ) { // if UCBlock has secondary demand variables
  if( !v_maximum_secondary_rho.empty() ) { // if unit produces any secondary reserve
   v_secondary_spinning_reserve.resize( f_time_horizon );
   for( auto & var : v_secondary_spinning_reserve ) {
    var.set_type( ColVariable::kNonNegative );
   }
   add_static_variable( v_secondary_spinning_reserve , "sc_battery" );
  }
 }

 set_variables_generated();
} // end( BatteryUnitBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::generate_abstract_constraints( Configuration * stcc ) {
 auto battery_type = get_battery_type();

 if( constraints_generated() )
  return; // constraints have already been generated

 if( !variables_generated() )
  throw ( std::logic_error( "BatteryUnitBlock::generate_abstract_constraints: "
                            "variables need be generated for constraints "
                            "to be." ) );

 int generate_ZOConstraint = 0;
 auto config = dynamic_cast<SimpleConfiguration< int > *>( stcc );
 if( ( !config ) && f_BlockConfig &&
     f_BlockConfig->f_static_constraints_Configuration )
  config = dynamic_cast< SimpleConfiguration< int > * >
  ( f_BlockConfig->f_static_constraints_Configuration );
 if( config )
  generate_ZOConstraint = config->f_value;

 // Initializing minimum power constraints
 active_power_lower_bound_Constraints.resize( f_time_horizon );

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {

  auto linear_function = new LinearFunction();

  linear_function->add_variable( &v_active_power[ t ] , 1.0 );
  if( reserve_vars & 1u ) { // if UCBlock has primary demand variables
   if( !v_maximum_primary_rho.empty() ) { // if unit produces any primary reserve
    linear_function->add_variable( &v_primary_spinning_reserve[ t ] , -1.0 );
   }
  }
  if( reserve_vars & 2u ) { // if UCBlock has secondary demand variables
   if( !v_maximum_secondary_rho.empty() ) { // if unit produces any secondary reserve
    linear_function->add_variable( &v_secondary_spinning_reserve[ t ] , -1.0 );
   }
  }
  active_power_lower_bound_Constraints[ t ].set_lhs( v_minimum_power[ t ] );
  active_power_lower_bound_Constraints[ t ].set_rhs( Inf< double >() );
  active_power_lower_bound_Constraints[ t ].set_function( linear_function );

 }
 add_static_constraint( active_power_lower_bound_Constraints ,
                        "ActivePower_LowerBound_Constraints_Battery" );

 // Initializing maximum power constraints

 active_power_upper_bound_Constraints.resize( f_time_horizon );

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {

  auto linear_function = new LinearFunction();

  linear_function->add_variable( &v_active_power[ t ] , 1.0 );
  if( reserve_vars & 1u ) { // if UCBlock has primary demand variables
   if( !v_maximum_primary_rho.empty() ) { // if unit produces any primary reserve
    linear_function->add_variable( &v_primary_spinning_reserve[ t ] , 1.0 );
   }
  }
  if( reserve_vars & 2u ) { // if UCBlock has secondary demand variable
   if( !v_maximum_secondary_rho.empty() ) { // if unit produces any secondary reserve
    linear_function->add_variable( &v_secondary_spinning_reserve[ t ] , 1.0 );
   }
  }
  active_power_upper_bound_Constraints[ t ].set_lhs( -Inf< double >() );
  active_power_upper_bound_Constraints[ t ].set_rhs( v_maximum_power[ t ] );
  active_power_upper_bound_Constraints[ t ].set_function( linear_function );
 }

 add_static_constraint( active_power_upper_bound_Constraints ,
                        "ActivePower_UpperBound_Constraints_Battery" );

/*--------------------------------------------------------------------------*/
 // Initializing ramp-up constraints
 if( !v_delta_ramp_up.empty() ) {

  ramp_up_Constraints.resize( f_time_horizon );

  auto linear_function = new LinearFunction();

  linear_function->add_variable( &v_active_power[ 0 ] , 1.0 );

  ramp_up_Constraints[ 0 ].set_lhs( -Inf< double >() );
  ramp_up_Constraints[ 0 ].set_rhs( v_delta_ramp_up[ 0 ] + f_initial_power );
  ramp_up_Constraints[ 0 ].set_function( linear_function );

  for( Index t = 1 , constraint_index = 1 ;
       t < f_time_horizon ; ++t , ++constraint_index ) {

   auto lf = new LinearFunction();

   lf->add_variable( &v_active_power[ t ] , 1.0 );
   lf->add_variable( &v_active_power[ t - 1 ] , -1.0 );

   ramp_up_Constraints[ constraint_index ].set_lhs( -Inf< double >() );
   ramp_up_Constraints[ constraint_index ].set_rhs( v_delta_ramp_up[ t ] );
   ramp_up_Constraints[ constraint_index ].set_function( lf );
  }
 }

 add_static_constraint( ramp_up_Constraints , "RampUp_Constraints_Battery" );

 // Initializing ramp-down constraints
 if( !v_delta_ramp_down.empty() ) {

  ramp_down_Constraints.resize( f_time_horizon );

  auto linear_function = new LinearFunction();

  linear_function->add_variable( &v_active_power[ 0 ] , 1.0 );

  ramp_down_Constraints[ 0 ].set_lhs(
   -v_delta_ramp_down[ 0 ] + f_initial_power );
  ramp_down_Constraints[ 0 ].set_rhs( Inf< double >() );
  ramp_down_Constraints[ 0 ].set_function( linear_function );

  for( Index t = 1 , constraint_index = 1 ;
       t < f_time_horizon ; ++t , ++constraint_index ) {

   auto lf = new LinearFunction();

   lf->add_variable( &v_active_power[ t ] , 1.0 );
   lf->add_variable( &v_active_power[ t - 1 ] , -1.0 );

   ramp_down_Constraints[ constraint_index ].set_lhs( -v_delta_ramp_down[ 0 ] );
   ramp_down_Constraints[ constraint_index ].set_rhs( Inf< double >() );
   ramp_down_Constraints[ constraint_index ].set_function( lf );
  }

 }
 add_static_constraint( ramp_down_Constraints ,
                        "RampDown_Constraints_Battery" );

/*--------------------------------------------------------------------------*/
// Initializing power_intake_outtake_Constraints

 power_intake_outtake_Constraints.resize( f_time_horizon );

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {

  auto linear_function = new LinearFunction();

  linear_function->add_variable( &v_active_power[ t ] , 1.0 );
  linear_function->add_variable( &v_intake_level[ t ] , -1.0 );
  linear_function->add_variable( &v_outtake_level[ t ] , 1.0 );

  power_intake_outtake_Constraints[ t ].set_both( 0.0 );
  power_intake_outtake_Constraints[ t ].set_function( linear_function );

 }

 add_static_constraint( power_intake_outtake_Constraints ,
                        "Power_Intake_Outtake_Constraints_Battery" );


// Initializing intake_upper_bound_Constraints


 intake_upper_bound_Constraints.resize( f_time_horizon );

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {

  intake_upper_bound_Constraints[ t ].set_lhs( 0.0 );
  intake_upper_bound_Constraints[ t ].set_rhs( v_maximum_power[ t ] );
  intake_upper_bound_Constraints[ t ].set_variable( &v_intake_level[ t ] );

 }
 add_static_constraint( intake_upper_bound_Constraints ,
                        "Intake_UpperBound_Constraints_Battery" );

 /*--------------------------------------------------------------------------*/
// Initializing demand_Constraints

 {
  demand_Constraints.resize( f_time_horizon );


  auto linear_fun = new LinearFunction();

  linear_fun->add_variable( &v_storage_level[ 0 ] , 1.0 );
  if( !v_storing_battery_rho.empty() ) {
   linear_fun->add_variable( &v_outtake_level[ 0 ] ,
                             -v_storing_battery_rho[ 0 ] );
  } else {
   linear_fun->add_variable( &v_outtake_level[ 0 ] , -1 );

  }
  if( !v_extracting_battery_rho.empty() ) {

   linear_fun->add_variable( &v_intake_level[ 0 ] ,
                             v_extracting_battery_rho[ 0 ] );
  } else {
   linear_fun->add_variable( &v_intake_level[ 0 ] , 1 );

  }
  if( !v_demand.empty() ) {
   demand_Constraints[ 0 ].set_both( ( f_initial_storage - v_demand[ 0 ] ) );
  } else {
   demand_Constraints[ 0 ].set_both( ( f_initial_storage ) );

  }
  demand_Constraints[ 0 ].set_function( linear_fun );

  for( Index t = 1 , constraint_index = 1 ;
       t < f_time_horizon ; ++t , ++constraint_index ) {

   auto linear_function = new LinearFunction();
   if( !v_extracting_battery_rho.empty() ) {
    linear_function->add_variable( &v_intake_level[ t ] ,
                                   v_extracting_battery_rho[ t ] );
   } else {
    linear_function->add_variable( &v_intake_level[ t ] , 1 );

   }
   if( !v_storing_battery_rho.empty() ) {
    linear_function->add_variable( &v_outtake_level[ t ] ,
                                   -v_storing_battery_rho[ t ] );
   } else {
    linear_function->add_variable( &v_outtake_level[ t ] , -1 );

   }
   linear_function->add_variable( &v_storage_level[ t ] , 1.0 );
   linear_function->add_variable( &v_storage_level[ t - 1 ] , -1.0 );

   if( !v_demand.empty() ) {
    demand_Constraints[ constraint_index ].set_both( -v_demand[ t ] );
   } else {
    demand_Constraints[ constraint_index ].set_both( 0.0 );

   }
   demand_Constraints[ constraint_index ].set_function( linear_function );
  }
 }

 add_static_constraint( demand_Constraints , "demand_Constraints_Battery" );

 /*--------------------------------------------------------------------------*/
// Initializing storage_level_bounds_Constraints

 storage_level_bounds_Constraints.resize( f_time_horizon );

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {

  storage_level_bounds_Constraints[ t ].set_lhs( v_minimum_storage[ t ] );
  storage_level_bounds_Constraints[ t ].set_rhs( v_maximum_storage[ t ] );
  storage_level_bounds_Constraints[ t ].set_variable( &v_storage_level[ t ] );

 }
 add_static_constraint( storage_level_bounds_Constraints ,
                        "StorageLevel_Bounds_Constraints_Battery" );

/*--------------------------------------------------------------------------*/

 // Initializing intake_binary_Constraints

 if( battery_type == Binary_Variables_Constraints ) {

  intake_binary_Constraints.resize( f_time_horizon );

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_intake_level[ t ] , 1.0 );
   linear_function->add_variable( &v_battery_binary[ t ] ,
                                  -v_maximum_power[ t ] );


   intake_binary_Constraints[ t ].set_lhs( -Inf< double >() );
   intake_binary_Constraints[ t ].set_rhs( 0.0 );
   intake_binary_Constraints[ t ].set_function( linear_function );

  }

  add_static_constraint( intake_binary_Constraints ,
                         "Intake_Binary_Constraints_Battery" );

/*--------------------------------------------------------------------------*/

  // Initializing outtake_binary_Constraints
  outtake_binary_Constraints.resize( f_time_horizon );

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_outtake_level[ t ] , 1.0 );
   linear_function->add_variable( &v_battery_binary[ t ] ,
                                  -v_minimum_power[ t ] );


   outtake_binary_Constraints[ t ].set_lhs( -Inf< double >() );
   outtake_binary_Constraints[ t ].set_rhs( -v_minimum_power[ t ] );
   outtake_binary_Constraints[ t ].set_function( linear_function );

  }
  add_static_constraint( outtake_binary_Constraints ,
                         "Outtake_Binary_Constraints_Battery" );
 }
/*--------------------------------------------------------------------------*/

 // Initializing primary_upper_bound_Constraints
 if( reserve_vars & 1u ) { // if UCBlock has primary demand variables
  if( !v_maximum_primary_rho.empty() ) { // if unit produces any primary reserve
   primary_upper_bound_Constraints.resize( f_time_horizon );

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {

    primary_upper_bound_Constraints[ t ].set_lhs( 0.0 );
    primary_upper_bound_Constraints[ t ].set_rhs( v_maximum_primary_rho[ t ] );
    primary_upper_bound_Constraints[ t ].set_variable(
     &v_primary_spinning_reserve[ t ] );

   }
   add_static_constraint( primary_upper_bound_Constraints ,
                          "Primary_UpperBound_Constraints_Battery" );
  }
 }
 // Initializing secondary_upper_bound_Constraints
 if( reserve_vars & 2u ) { // if UCBlock has secondary demand variables
  if( !v_maximum_secondary_rho.empty() ) { // if unit produces any secondary reserve
   secondary_upper_bound_Constraints.resize( f_time_horizon );

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {

    secondary_upper_bound_Constraints[ t ].set_lhs( 0.0 );
    secondary_upper_bound_Constraints[ t ].set_rhs(
     v_maximum_secondary_rho[ t ] );
    secondary_upper_bound_Constraints[ t ].set_variable(
     &v_secondary_spinning_reserve[ t ] );

   }
   add_static_constraint( secondary_upper_bound_Constraints ,
                          "Secondary_UpperBound_Constraints_Battery" );
  }
 }
/*-------------------------------ZOConstraint-------------------------------*/
 if( battery_type == Binary_Variables_Constraints ) {

  if( generate_ZOConstraint ) {
   // the battery binary bound constraints
   battery_binary_bound_Constraints.resize( f_time_horizon );
   for( Index t = 0 ; t < f_time_horizon ; ++t ) {
    battery_binary_bound_Constraints[ t ].set_variable(
     &v_battery_binary[ t ] );
   }
   add_static_constraint( battery_binary_bound_Constraints ,
                          "BB_bound_battery" );
  }
 }
 set_constraints_generated();
} // end( BatteryUnitBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::generate_objective( Configuration * objc ) {
 if( objective_generated() )
  return; // Objective has already been generated

 if( !variables_generated() )
  throw ( std::logic_error( "BatteryUnitBlock::generate_objective: variables "
                            "need be generated for constraints to be." ) );

 if( get_objective() != nullptr )  // an objective is there already
  return;                         // cowardly (and silently) return

 // Construct the costs vector

 std::vector< double > cost = v_cost;
 if( cost.size() == 1 )
  cost.resize( f_time_horizon , cost[ 0 ] );
 else if( cost.size() < f_time_horizon ) {
  assert( cost.size() == v_change_intervals.size() );
  cost.resize( f_time_horizon , 0 );
  Index t = 0;
  for( Index k = 0 ; k < v_change_intervals.size() ; ++k ) {
   auto upper_endpoint = v_change_intervals[ k ];
   if( k == v_change_intervals.size() - 1 )
    // The upper endpoint of the last interval must be time_horizon - 1. Since
    // it may not be provided in v_change_intervals (the value for the last
    // element of v_change_intervals is not required), we manually set it
    // here.
    upper_endpoint = f_time_horizon - 1;
   for( ; t <= upper_endpoint ; ++t ) {
    cost[ t ] = v_cost[ k ];
   }
  }
 }

 // initialize objective function - - - - - - - - - - - - - - - - - - - - - -

 auto linear_function = new LinearFunction();

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {
  linear_function->add_variable( &v_intake_level[ t ] , cost[ t ] , 0.0 );
  linear_function->add_variable( &v_outtake_level[ t ] , cost[ t ] , 0.0 );
 }

 objective.set_function( linear_function );
 objective.set_sense( Objective::eMin );

 // Set Block objective
 this->set_objective( &objective );

 set_objective_generated();

}  // end( BatteryUnitBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*--- METHODS FOR LOADING, PRINTING & SAVING THE BatteryUnitBlock ---*/
/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::serialize( netCDF::NcGroup & group ) const {

 UnitBlock::serialize( group );

 // Serialize scalar variables.

 ::serialize( group , "InitialPower" , netCDF::NcDouble() , f_initial_power );
 ::serialize( group , "InitialStorage" , netCDF::NcDouble() ,
              f_initial_storage );

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
   throw ( std::logic_error
    ( "BatteryUnitBlock::serialize: invalid dimension for variable " +
      var_name + ": " + std::to_string( data.size() ) + ". Its dimension " +
      "must be one of the following: TimeHorizon, NumberIntervals, 1." ) );
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

 if( demand_Constraints.empty() )
  return;

 if( !v_demand.empty() )
  demand_Constraints[ 0 ].set_both( f_initial_storage - v_demand[ 0 ] ,
                                    issueAMod );
 else
  demand_Constraints[ 0 ].set_both( f_initial_storage , issueAMod );
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
 if( !( rng.first <= 0 && 0 < rng.second ) )
  return; // 0 does not belong to the range; return

 std::advance( it , -rng.first );

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
 if( !( ramp_up_Constraints.empty() || v_delta_ramp_up.empty() ) )
  ramp_up_Constraints[ 0 ].set_rhs( v_delta_ramp_up[ 0 ] + f_initial_power ,
                                    issueAMod );

 if( !( ramp_down_Constraints.empty() || v_delta_ramp_down.empty() ) )
  ramp_down_Constraints[ 0 ].set_lhs( -v_delta_ramp_down[ 0 ] +
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
 if( !( rng.first <= 0 && 0 < rng.second ) )
  return; // 0 does not belong to the range; return

 std::advance( it , -rng.first );

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
/*----------------- End File BatteryUnitBlock.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
