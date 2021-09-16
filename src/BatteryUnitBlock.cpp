/*--------------------------------------------------------------------------*/
/*----------------- File BatteryUnitBlock.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the BatteryStorageUnitBlock class.
 *
 * \version 0.11
 *
 * \date 19 - 08 - 2021
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 *
 * \author Ali Ghezelsoflu \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 *
 * Copyright &copy by Antonio Frangioni, Ali Ghezelsoflu
 */

/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/

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
 std::vector< std::string > expected_dims = { "TimeHorizon",
                                              "NumberIntervals" };
 check_dimensions( group, expected_dims, std::cerr );
 std::vector< std::string > expected_vars = { "MinStorage",
                                              "MaxStorage",
                                              "MinPower",
                                              "MaxPower",
                                              "InitialPower",
                                              "MaxPrimaryPower",
                                              "MaxSecondaryPower",
                                              "DeltaRampUp",
                                              "DeltaRampDown",
                                              "StoringBatteryRho",
                                              "ExtractingBatteryRho",
                                              "InitialStorage",
                                              "Cost",
                                              "Demand" };
 check_variables( group, expected_vars, std::cerr );
#endif


 UnitBlock::deserialize_time_horizon( group );
 UnitBlock::deserialize_change_intervals( group );

 ::deserialize( group, "MinStorage",f_time_horizon, v_minimum_storage, true, true );
 ::deserialize( group, "MaxStorage", f_time_horizon, v_maximum_storage, true, true );
 ::deserialize( group, "MinPower", f_time_horizon, v_minimum_power, true,true);
 ::deserialize( group, "MaxPower", f_time_horizon, v_maximum_power, true,true);
 ::deserialize( group, "InitialPower", &f_initial_power );
 ::deserialize( group, "MaxPrimaryPower",f_time_horizon, v_maximum_primary_rho, true,true);
 ::deserialize( group, "MaxSecondaryPower",f_time_horizon, v_maximum_secondary_rho, true,true );
 ::deserialize( group, "DeltaRampUp",f_time_horizon, v_delta_ramp_up, true,true );
 ::deserialize( group, "DeltaRampDown",f_time_horizon, v_delta_ramp_down, true,true);
 ::deserialize( group, "StoringBatteryRho",f_time_horizon, v_storing_battery_rho, true,true );
 ::deserialize( group, "ExtractingBatteryRho",f_time_horizon, v_extracting_battery_rho, true,true );
 ::deserialize( group, "Cost", f_time_horizon, v_cost, true,true);
 ::deserialize( group, "Demand", f_time_horizon, v_demand, true, false );

 ::deserialize( group, "InitialStorage", &f_initial_storage );

 decompress_vector(v_minimum_power);
 decompress_vector(v_maximum_power);
 decompress_vector(v_minimum_storage);
 decompress_vector(v_maximum_storage);
 decompress_vector(v_maximum_primary_rho);
 decompress_vector(v_maximum_secondary_rho);
 decompress_vector(v_delta_ramp_up);
 decompress_vector(v_delta_ramp_down);
 decompress_vector(v_storing_battery_rho);
 decompress_vector(v_extracting_battery_rho);
 decompress_vector(v_demand);

 UnitBlock::deserialize( group );

}// end( BatteryUnitBlock::deserialize )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::generate_abstract_variables( Configuration *stvv ) {
 auto battery_type = get_battery_type();

 if( variables_generated() )
  return; // variables have already been generated

 UnitBlock::generate_abstract_variables( stvv );

 int relax_binary = 0;
 auto config = dynamic_cast<SimpleConfiguration<int> *>( stvv );
 if( ( ! config ) && f_BlockConfig &&
     f_BlockConfig->f_static_variables_Configuration )
  config = dynamic_cast< SimpleConfiguration< int > * >
   ( f_BlockConfig->f_static_variables_Configuration );
 if( config )
  relax_binary = config->f_value;

 v_storage_level.resize( f_time_horizon );
 for( auto & var : v_storage_level )
  var.set_type( ColVariable::kNonNegative );
 add_static_variable( v_storage_level, "SL_battery" );


 v_intake_level.resize( f_time_horizon );
 for( auto & var : v_intake_level )
  var.set_type( ColVariable::kNonNegative );
 add_static_variable( v_intake_level, "IL_battery" );


 v_outtake_level.resize( f_time_horizon );
 for( auto & var : v_outtake_level )
  var.set_type( ColVariable::kNonNegative );
 add_static_variable( v_outtake_level, "OL_battery" );


 v_battery_binary.resize( f_time_horizon );
 for( auto & var : v_battery_binary ) {
  if( relax_binary )
   var.set_type( ColVariable::kPosUnitary );
  else
   var.set_type( ColVariable::kBinary );
 }
 if ( battery_type == Binary_Variables_Constraints ) {

  add_static_variable( v_battery_binary, "BB_battery" );
 }
  // Active Power Variable

 v_active_power.resize( f_time_horizon );
 for( auto & var : v_active_power )
  var.set_type( ColVariable::kContinuous );
 add_static_variable( v_active_power, "p_battery" );

  // Primary Spinning Reserve Variable
 if( reserve_vars & 1u ) { // if UCBlock has primary demand variables
  if (!v_maximum_primary_rho.empty() ) { // if unit produces any primary reserve
   v_primary_spinning_reserve.resize( f_time_horizon );
   for( auto & var : v_primary_spinning_reserve ) {
    var.set_type( ColVariable::kNonNegative );
   }
   add_static_variable( v_primary_spinning_reserve, "pr_battery" );
  }
 }

 // Secondary Spinning Reserve Variable
 if( reserve_vars & 2u ) { // if UCBlock has secondary demand variables
  if (!v_maximum_secondary_rho.empty() ) { // if unit produces any secondary reserve
   v_secondary_spinning_reserve.resize( f_time_horizon );
   for( auto & var : v_secondary_spinning_reserve ) {
    var.set_type( ColVariable::kNonNegative );
   }
   add_static_variable( v_secondary_spinning_reserve, "sc_battery" );
  }
 }

 set_variables_generated();
} // end( BatteryUnitBlock::generate_abstract_variables )
/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::generate_abstract_constraints ( Configuration * stcc ) {
 auto battery_type = get_battery_type();

 if( constraints_generated() )
  return; // constraints have already been generated

 int generate_ZOConstraint = 0;
 auto config = dynamic_cast<SimpleConfiguration<int> *>( stcc );
 if( ( ! config ) && f_BlockConfig &&
     f_BlockConfig->f_static_constraints_Configuration )
  config = dynamic_cast< SimpleConfiguration< int > * >
   ( f_BlockConfig->f_static_constraints_Configuration );
 if( config )
  generate_ZOConstraint = config->f_value;

 // Initial data check
   for( Index t = 0; t < f_time_horizon; ++t ) {
    if( v_minimum_power[t] >= v_maximum_power[t] ) {
     throw ( std::logic_error
             ( "BatteryUnitBlock::maximum and minimum power output constraints: "
               "it must be that v_maximum_power > v_minimum_power." ));
    }
   }

 for( Index t = 0; t < f_time_horizon; ++t ) {
  if( v_minimum_storage[t] > v_maximum_storage[t] ||
          v_minimum_storage[t] < 0  || v_maximum_storage[t] < 0 ) {
   throw ( std::logic_error
           ( "BatteryUnitBlock::maximum and minimum storage output constraints: "
             "it must be that v_maximum_storage >= v_minimum_storage >= 0." ));
  }
 }

 if (!v_storing_battery_rho.empty() & !v_extracting_battery_rho.empty()) {
  for( Index t = 0; t < f_time_horizon; ++t ) {
   if( v_extracting_battery_rho[t] < 1 || v_storing_battery_rho[t] > 1 ||
           v_extracting_battery_rho[t]  < v_storing_battery_rho[t] ) {
    throw ( std::logic_error
            ( "BatteryUnitBlock::storing_battery_rho(SBR) and "
              "extracting_battery_rho(EBR): it must be that "
              "EBR[ t ] >= 1 [>= SBR[ t ]] " ));
   }
  }
 }

 if (v_maximum_primary_rho.empty() & !v_maximum_secondary_rho.empty()) {
    throw ( std::logic_error
            ( "BatteryUnitBlock:: if MaxPrimaryPower is not defined then "
              "neither should MaxSecondaryPower." ));
 }
 // Initializing minimum power constraints
  active_power_lower_bound_Constraints.resize( f_time_horizon );

  for( Index t = 0; t < f_time_horizon; ++t ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_active_power[t], 1.0 );
   if( reserve_vars & 1u ) { // if UCBlock has primary demand variables
    if( !v_maximum_primary_rho.empty()) { // if unit produces any primary reserve
     linear_function->add_variable( &v_primary_spinning_reserve[t], -1.0 );
    }
   }
   if( reserve_vars & 2u ) { // if UCBlock has secondary demand variables
    if( !v_maximum_secondary_rho.empty()) { // if unit produces any secondary reserve
     linear_function->add_variable( &v_secondary_spinning_reserve[t], -1.0 );
    }
   }
    active_power_lower_bound_Constraints[t].set_lhs( v_minimum_power[t] );
    active_power_lower_bound_Constraints[t].set_rhs( Inf< double >());
    active_power_lower_bound_Constraints[t].set_function( linear_function );

  }
  add_static_constraint( active_power_lower_bound_Constraints, "ActivePower_LowerBound_Constraints_Battery" );

  // Initializing maximum power constraints

   active_power_upper_bound_Constraints.resize( f_time_horizon );

   for( Index t = 0; t < f_time_horizon; ++t ) {

    auto linear_function = new LinearFunction();

    linear_function->add_variable( &v_active_power[t], 1.0 );
    if( reserve_vars & 1u ) { // if UCBlock has primary demand variables
     if( !v_maximum_primary_rho.empty()) { // if unit produces any primary reserve
      linear_function->add_variable( &v_primary_spinning_reserve[t], 1.0 );
     }
    }
    if( reserve_vars & 2u ) { // if UCBlock has secondary demand variable
     if( !v_maximum_secondary_rho.empty()) { // if unit produces any secondary reserve
      linear_function->add_variable( &v_secondary_spinning_reserve[t], 1.0 );
     }
    }
    active_power_upper_bound_Constraints[t].set_lhs( -Inf< double >());
    active_power_upper_bound_Constraints[t].set_rhs( v_maximum_power[t] );
    active_power_upper_bound_Constraints[t].set_function( linear_function );
   }

  add_static_constraint( active_power_upper_bound_Constraints, "ActivePower_UpperBound_Constraints_Battery" );

/*--------------------------------------------------------------------------*/
 // Initializing ramp-up constraints
 if( !v_delta_ramp_up.empty()) {

   ramp_up_Constraints.resize( f_time_horizon );

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_active_power[0], 1.0 );

   ramp_up_Constraints[0].set_lhs( -Inf< double >());
   ramp_up_Constraints[0].set_rhs( v_delta_ramp_up[0] + f_initial_power );
   ramp_up_Constraints[0].set_function( linear_function );

   for( Index t = 1, constraint_index = 1; t < f_time_horizon; ++t, ++constraint_index ) {

    auto lf = new LinearFunction();

    lf->add_variable( &v_active_power[t], 1.0 );
    lf->add_variable( &v_active_power[t - 1], -1.0 );

    ramp_up_Constraints[constraint_index].set_lhs( -Inf< double >() );
    ramp_up_Constraints[constraint_index].set_rhs( v_delta_ramp_up[t] );
    ramp_up_Constraints[constraint_index].set_function( lf );
   }
 }

 add_static_constraint( ramp_up_Constraints, "RampUp_Constraints_Battery" );

 // Initializing ramp-down constraints
 if( !v_delta_ramp_down.empty() ) {

  ramp_down_Constraints.resize( f_time_horizon );

  auto linear_function = new LinearFunction();

  linear_function->add_variable( &v_active_power[ 0 ], 1.0 );

   ramp_down_Constraints[ 0 ].set_lhs ( -v_delta_ramp_down[ 0 ] + f_initial_power );
   ramp_down_Constraints[ 0 ].set_rhs( Inf< double >() );
   ramp_down_Constraints[ 0 ].set_function( linear_function );

  for( Index t = 1, constraint_index = 1; t < f_time_horizon; ++t, ++constraint_index ) {

   auto lf = new LinearFunction();

   lf->add_variable( &v_active_power[t], 1.0 );
   lf->add_variable( &v_active_power[t - 1], -1.0 );

   ramp_down_Constraints[constraint_index].set_lhs( -v_delta_ramp_down[ 0 ] );
   ramp_down_Constraints[constraint_index].set_rhs( Inf< double >() );
   ramp_down_Constraints[constraint_index].set_function( lf );
  }

 }
 add_static_constraint( ramp_down_Constraints, "RampDown_Constraints_Battery" );

/*--------------------------------------------------------------------------*/
// Initializing power_intake_outtake_Constraints

  power_intake_outtake_Constraints.resize( f_time_horizon );

  for( Index t = 0; t < f_time_horizon; ++t ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_active_power[t], 1.0 );
   linear_function->add_variable( &v_intake_level[t], -1.0 );
   linear_function->add_variable( &v_outtake_level[t], 1.0 );

   power_intake_outtake_Constraints[t].set_both( 0.0);
   power_intake_outtake_Constraints[t].set_function( linear_function );

  }

 add_static_constraint( power_intake_outtake_Constraints, "Power_Intake_Outtake_Constraints_Battery" );


// Initializing intake_upper_bound_Constraints


  intake_upper_bound_Constraints.resize( f_time_horizon );

  for( Index t = 0; t < f_time_horizon; ++t ) {

   intake_upper_bound_Constraints[t].set_lhs( 0.0);
   intake_upper_bound_Constraints[t].set_rhs( v_maximum_power[ t ]);
   intake_upper_bound_Constraints[t].set_variable(&v_intake_level[t]);

  }
 add_static_constraint( intake_upper_bound_Constraints, "Intake_UpperBound_Constraints_Battery" );

 /*--------------------------------------------------------------------------*/
// Initializing demand_Constraints

 {
  demand_Constraints.resize( f_time_horizon );


  auto linear_fun = new LinearFunction();

  linear_fun->add_variable( &v_storage_level[0], 1.0 );
  if (!v_storing_battery_rho.empty()) {
   linear_fun->add_variable( &v_outtake_level[0], -v_storing_battery_rho[0] );
  } else {
   linear_fun->add_variable( &v_outtake_level[0], -1 );

  }
  if (!v_extracting_battery_rho.empty()) {

   linear_fun->add_variable( &v_intake_level[0], v_extracting_battery_rho[0] );
  } else{
   linear_fun->add_variable( &v_intake_level[0], 1 );

  }
  if (!v_demand.empty()) {
   demand_Constraints[0].set_both(( f_initial_storage - v_demand[0] ));
  } else {
   demand_Constraints[0].set_both(( f_initial_storage ));

  }
  demand_Constraints[0].set_function( linear_fun );

  for( Index t = 1, constraint_index = 1; t < f_time_horizon; ++t, ++constraint_index ) {

   auto linear_function = new LinearFunction();
   if (!v_extracting_battery_rho.empty()) {
    linear_function->add_variable( &v_intake_level[t], v_extracting_battery_rho[t] );
   } else{
    linear_function->add_variable( &v_intake_level[t], 1 );

   }
   if (!v_storing_battery_rho.empty()) {
    linear_function->add_variable( &v_outtake_level[t], -v_storing_battery_rho[t] );
   } else{
    linear_function->add_variable( &v_outtake_level[t], -1 );

   }
   linear_function->add_variable( &v_storage_level[t], 1.0 );
   linear_function->add_variable( &v_storage_level[t-1], -1.0 );

   if (!v_demand.empty()) {
    demand_Constraints[constraint_index].set_both( -v_demand[t] );
   } else {
    demand_Constraints[constraint_index].set_both( 0.0);

   }
   demand_Constraints[constraint_index].set_function( linear_function );
  }
 }

 add_static_constraint( demand_Constraints, "demand_Constraints_Battery" );

 /*--------------------------------------------------------------------------*/
// Initializing storage_level_bounds_Constraints

  storage_level_bounds_Constraints.resize( f_time_horizon );

  for( Index t = 0; t < f_time_horizon; ++t ) {

   storage_level_bounds_Constraints[t].set_lhs( v_minimum_storage[ t ]);
   storage_level_bounds_Constraints[t].set_rhs( v_maximum_storage[ t ]);
   storage_level_bounds_Constraints[t].set_variable(&v_storage_level[t]);

  }
 add_static_constraint( storage_level_bounds_Constraints, "StorageLevel_Bounds_Constraints_Battery" );

/*--------------------------------------------------------------------------*/

 // Initializing intake_binary_Constraints

 if ( battery_type == Binary_Variables_Constraints ) {

  intake_binary_Constraints.resize( f_time_horizon );

  for( Index t = 0; t < f_time_horizon; ++t ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_intake_level[t], 1.0 );
   linear_function->add_variable( &v_battery_binary[t],  -v_maximum_power[ t ] );


   intake_binary_Constraints[t].set_lhs( -Inf<double>());
   intake_binary_Constraints[t].set_rhs( 0.0 );
   intake_binary_Constraints[t].set_function( linear_function );

  }

 add_static_constraint( intake_binary_Constraints, "Intake_Binary_Constraints_Battery" );

/*--------------------------------------------------------------------------*/

 // Initializing outtake_binary_Constraints
  outtake_binary_Constraints.resize( f_time_horizon );

  for( Index t = 0; t < f_time_horizon; ++t ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_outtake_level[ t ], 1.0 );
   linear_function->add_variable( &v_battery_binary[ t ],  - v_minimum_power[ t ] );


   outtake_binary_Constraints[t].set_lhs( -Inf<double>());
   outtake_binary_Constraints[t].set_rhs( -v_minimum_power[ t ] );
   outtake_binary_Constraints[t].set_function( linear_function );

  }
  add_static_constraint( outtake_binary_Constraints, "Outtake_Binary_Constraints_Battery" );
 }
/*--------------------------------------------------------------------------*/

  // Initializing primary_upper_bound_Constraints
 if( reserve_vars & 1u ) { // if UCBlock has primary demand variables
  if( !v_maximum_primary_rho.empty()) { // if unit produces any primary reserve
   primary_upper_bound_Constraints.resize( f_time_horizon );

   for( Index t = 0; t < f_time_horizon; ++t ) {

    primary_upper_bound_Constraints[t].set_lhs( 0.0 );
    primary_upper_bound_Constraints[t].set_rhs( v_maximum_primary_rho[t] );
    primary_upper_bound_Constraints[t].set_variable( &v_primary_spinning_reserve[t] );

   }
   add_static_constraint( primary_upper_bound_Constraints, "Primary_UpperBound_Constraints_Battery" );
  }
 }
  // Initializing secondary_upper_bound_Constraints
 if( reserve_vars & 2u ) { // if UCBlock has secondary demand variables
  if( !v_maximum_secondary_rho.empty()) { // if unit produces any secondary reserve
   secondary_upper_bound_Constraints.resize( f_time_horizon );

   for( Index t = 0; t < f_time_horizon; ++t ) {

    secondary_upper_bound_Constraints[t].set_lhs( 0.0 );
    secondary_upper_bound_Constraints[t].set_rhs( v_maximum_secondary_rho[t] );
    secondary_upper_bound_Constraints[t].set_variable( &v_secondary_spinning_reserve[t] );

   }
   add_static_constraint( secondary_upper_bound_Constraints,
                          "Secondary_UpperBound_Constraints_Battery" );
  }
 }
/*-------------------------------ZOConstraint-------------------------------*/
 if ( battery_type == Binary_Variables_Constraints ) {

  if( generate_ZOConstraint ) {
   // the battery binary bound constraints
   battery_binary_bound_Constraints.resize( f_time_horizon );
   for( Index t = 0; t < f_time_horizon; ++t ) {
    battery_binary_bound_Constraints[t].set_variable( &v_battery_binary[t] );
   }
   add_static_constraint( battery_binary_bound_Constraints,
                          "BB_bound_battery" );
  }
 }
 set_constraints_generated();
} // end( BatteryUnitBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::generate_objective( Configuration *objc )
{
 if( objective_generated() )
  return; // Objective has already been generated

 // Initial condition of each vector

 std::vector<double> cost = v_cost;
 if (cost.size() == 1) {
  cost.resize(f_time_horizon, cost[0]);

 } else if (cost.size() < f_time_horizon) {
  cost.resize(f_time_horizon);
  int j = 0;
  for (unsigned long i = 0; i < v_change_intervals.size(); ++i) {
   Index sup;
   if (i == v_change_intervals.size() - 1 ) {
    sup = f_time_horizon;
   } else {
    sup = v_change_intervals[i];
   }
   for (; j < sup; ++j) {
    cost[j] = v_cost[i];
   }
  }
 }
 // initialize objective function - - - - - - - - - - - - - - - - - - - - - -

 if( get_objective() != nullptr )  // an objective is there already
  return;                         // cowardly (and silently) return


 if( v_intake_level.size() != f_time_horizon  ) {
  throw ( std::logic_error
          ( "BatteryUnitBlock::generate_objective: v_intake_level must have "
            "size equal to the time horizon." ));
 }

 if(  v_outtake_level.size() != f_time_horizon) {
  throw ( std::logic_error
          ( "BatteryUnitBlock::generate_objective: v_outtake_level must have "
            "size equal to the time horizon." ));
 }
  auto linear_function = new LinearFunction();

  for( Index t = 0; t < f_time_horizon; ++t ) {
   linear_function->add_variable( &v_intake_level[ t ],
                                 cost[ t  ],
                                 0.0 );

   linear_function->add_variable( &v_outtake_level[ t ],
                                 cost[ t  ],
                                 0.0 );
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

 auto NumberIntervals = group.getDim( "NumberIntervals" );
 auto TimeHorizon = group.getDim( "TimeHorizon" );

 ::serialize( group, "InitialPower", netCDF::NcDouble(), f_initial_power );
 ::serialize( group, "InitialStorage", netCDF::NcDouble(), f_initial_storage );

 ::serialize( group, "MinStorage", netCDF::NcDouble(),
              NumberIntervals, v_minimum_storage, true );

 ::serialize( group, "MaxStorage", netCDF::NcDouble(),
              NumberIntervals, v_maximum_storage, true );

 ::serialize( group, "MinPower", netCDF::NcDouble(),
              NumberIntervals, v_minimum_power, true );

 ::serialize( group, "MaxPower", netCDF::NcDouble(),
              NumberIntervals, v_maximum_power, true );

 ::serialize( group, "MaxPrimaryRho", netCDF::NcDouble(),
              NumberIntervals, v_maximum_primary_rho, true );

 ::serialize( group, "MaxSecondaryRho", netCDF::NcDouble(),
              NumberIntervals, v_maximum_secondary_rho, true );

 ::serialize( group, "DeltaRampUp", netCDF::NcDouble(),
              NumberIntervals, v_delta_ramp_up, true );

 ::serialize( group, "DeltaRampDown", netCDF::NcDouble(),
              NumberIntervals, v_delta_ramp_down, true );

 ::serialize( group, "StoringBatteryRho", netCDF::NcDouble(),
              NumberIntervals, v_storing_battery_rho, true );

 ::serialize( group, "ExtractingBatteryRho", netCDF::NcDouble(),
              NumberIntervals, v_extracting_battery_rho, true );

 ::serialize( group, "Cost", netCDF::NcDouble(),
              NumberIntervals, v_cost, true );

 ::serialize( group, "Demand", netCDF::NcDouble(),
              TimeHorizon, v_demand, true );
}  // end( BatteryUnitBlock::serialize )

/*--------------------------------------------------------------------------*/

template< typename T >
void BatteryUnitBlock::decompress_vector( std::vector< T > & v ) {

 if ( v.empty() )
  return;

 if( v.size() == 1 ) {
  // The given vector has a single element. Thus, for each time instant, the
  // value is equal to that single given element.
  v.resize( f_time_horizon , v[ 0 ] );
 }
 else if( v.size() < f_time_horizon ) {
  // Since the number of elements is greater than 1 and less than the time
  // horizon, it must be equal to the number of change intervals.
  if( v.size() != v_change_intervals.size() ) {
   throw ( std::logic_error
           ( "BatteryUnitBlock::decompress_vector: invalid number of elements"
             " (" + std::to_string( v.size() ) + ") for some variable. It "
             "should be equal to the number of change intervals (" +
             std::to_string( v_change_intervals.size() ) + ")" ) );
  }

  // For each time instant t, the value associated with time t is equal to
  // given_vector[ k ], where k is such that t belongs to the closed interval
  // [i_{k-1} + 1, i_k] and i_k is the k-th element of v_change_intervals
  // (starting from k = 0) and i_{-1} = -1 by definition. We resize the vector
  // so that its size becomes f_time_horizon and copy the given data.

  std::vector< T > given_vector = v;
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
    v[ t ] = given_vector[ k ];
   }
  }
 }
}

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::update_initial_storage_in_constraints
( c_ModParam issueAMod ) {

 if( demand_Constraints.empty() )
  return;

 if( ! v_demand.empty() )
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
 if( ! ( ramp_up_Constraints.empty() || v_delta_ramp_up.empty() ) )
  ramp_up_Constraints[ 0 ].set_rhs( v_delta_ramp_up[ 0 ] + f_initial_power ,
                                    issueAMod );

 if( ! ( ramp_down_Constraints.empty() || v_delta_ramp_down.empty() ) )
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
/*----------------- End File BatteryUnitBlock.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
