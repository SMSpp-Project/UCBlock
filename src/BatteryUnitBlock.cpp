/*--------------------------------------------------------------------------*/
/*----------------- File BatteryUnitBlock.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the BatteryStorageUnitBlock class.
 *
 * \version 0.11
 *
 * \date 18 - 07 - 2019
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
#include "BatteryUnitBlock.h"
#include "LinearFunction.h"
#include <map>
#include "FRowConstraint.h"
#include "DQuadFunction.h"
#include "FRealObjective.h"
#include "UnitBlock.h"
#include "UCBlock.h"

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
/*------------------- METHODS OF BatteryUnitBlock -------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
void BatteryUnitBlock::deserialize( netCDF::NcGroup & group ) {

 UnitBlock::deserialize( group );

 ::deserialize( group, "MinStorage", f_number_intervals, v_minimum_storage, true, true );

 ::deserialize( group, "MaxStorage", f_number_intervals, v_maximum_storage, true, true );

 ::deserialize( group, "MinPower", f_number_intervals, v_minimum_power, true, true );

 ::deserialize( group, "MaxPower", f_number_intervals, v_maximum_power, true, true );

 ::deserialize( group, "InitialPower", &f_initial_power );

 ::deserialize( group, "MaxPrimaryPower", f_number_intervals, v_maximum_primary_rho, true, true );

 ::deserialize( group, "MaxSecondaryPower", f_number_intervals, v_maximum_secondary_rho, true, true );

 ::deserialize( group, "DeltaRampUp", f_number_intervals, v_delta_ramp_up, true, true );

 ::deserialize( group, "DeltaRampDown", f_number_intervals, v_delta_ramp_down, true, true );

 ::deserialize( group, "StoringBatteryRho", f_number_intervals, v_storing_battery_rho, true, true );

 ::deserialize( group, "ExtractingBatterRho", f_number_intervals, v_extracting_battery_rho, true, true );

 ::deserialize( group, "InitialStorage", &f_initial_storage );

 ::deserialize( group, "Cost", f_number_intervals, v_cost, true, true );

 ::deserialize( group, "Demand", f_time_horizon, v_demand, true, false );

}// end( BatteryUnitBlock::deserialize )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::generate_abstract_variables  ( Configuration *stvv ) {
 UnitBlock::generate_abstract_variables( stvv );

 if( f_time_horizon == 0 ) {
  // there are no variables to be generated
  return;
 }

 if( f_time_horizon > 0 ) {

  if( v_storage_level.size() != f_time_horizon ) {
   assert( v_storage_level.empty()); // this should only happen once
   v_storage_level.resize( f_time_horizon );
   int n = 0;
   for( auto & i : v_storage_level ) {
    i.set_type( ColVariable::kNonNegative );
    //TODO look ThermalUnitBlock variables
    add_static_variable( i, "StorageLevel " + std::to_string( n++ ));
   }
  }

  if( v_intake_level.size() != f_time_horizon ) {
   assert( v_intake_level.empty()); // this should only happen once
   v_intake_level.resize( f_time_horizon );
   int n = 0;
   for( auto & i : v_intake_level ) {
    i.set_type( ColVariable::kNonNegative );
    //TODO look ThermalUnitBlock variables
    add_static_variable( i, "IntakeLevel " + std::to_string( n++ ));
   }
  }

  if( v_outtake_level.size() != f_time_horizon ) {
   assert( v_outtake_level.empty()); // this should only happen once
   v_outtake_level.resize( f_time_horizon );
   int n = 0;
   for( auto & i : v_outtake_level ) {
    i.set_type( ColVariable::kNonNegative );
    //TODO look ThermalUnitBlock variables
    add_static_variable( i, "OuttakeLevel " + std::to_string( n++ ));
   }
  }

  if( v_battery_binary.size() != f_time_horizon ) {
   assert( v_battery_binary.empty()); // this should only happen once
   v_battery_binary.resize( f_time_horizon );
   int n = 0;
   for( auto & i : v_battery_binary ) {
    i.set_type( ColVariable::kBinary );
    //TODO look ThermalUnitBlock variables
    add_static_variable( i, "BatteryBinary " + std::to_string( n++ ));
   }
  }
 } // end( BatteryUnitBlock::generate_abstract_variables )
}
/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::generate_abstract_constraints ( Configuration *stcc )
{
 // initial condition of each vector
 std::vector<double> min_power = v_minimum_power;
 if (min_power.size() == 1) {
  min_power.resize(f_time_horizon, min_power[0]);

 } else if (min_power.size() < f_time_horizon) {
  min_power.resize(f_time_horizon);
  int j = 0;
  for (unsigned long i = 0; i < v_change_intervals.size(); ++i) {
   Index sup;
   if (i == v_change_intervals.size() - 1 ) {
    sup = f_time_horizon;
   } else {
    sup = v_change_intervals[i];
   }
   for (; j < sup; ++j) {
    min_power[j] = v_minimum_power[i];
   }
  }
 }

 std::vector<double> max_power = v_maximum_power;
 if (max_power.size() == 1) {
  max_power.resize(f_time_horizon, max_power[0]);

 } else if (max_power.size() < f_time_horizon) {
  max_power.resize(f_time_horizon);
  int j = 0;
  for (unsigned long i = 0; i < v_change_intervals.size(); ++i) {
   Index sup;
   if (i == v_change_intervals.size() - 1 ) {
    sup = f_time_horizon;
   } else {
    sup = v_change_intervals[i];
   }
   for (; j < sup; ++j) {
    max_power[j] = v_maximum_power[i];
   }
  }
 }

 std::vector<double> min_storage = v_minimum_storage;
 if (min_storage.size() == 1) {
  min_storage.resize(f_time_horizon, min_storage[0]);

 } else if (min_storage.size() < f_time_horizon) {
  min_storage.resize(f_time_horizon);
  int j = 0;
  for (unsigned long i = 0; i < v_change_intervals.size(); ++i) {
   Index sup;
   if (i == v_change_intervals.size() - 1 ) {
    sup = f_time_horizon;
   } else {
    sup = v_change_intervals[i];
   }
   for (; j < sup; ++j) {
    min_storage[j] = v_minimum_storage[i];
   }
  }
 }

 std::vector<double> max_storage = v_maximum_storage;
 if (max_storage.size() == 1) {
  max_storage.resize(f_time_horizon, max_storage[0]);

 } else if (max_storage.size() < f_time_horizon) {
  max_storage.resize(f_time_horizon);
  int j = 0;
  for (unsigned long i = 0; i < v_change_intervals.size(); ++i) {
   Index sup;
   if (i == v_change_intervals.size() - 1 ) {
    sup = f_time_horizon;
   } else {
    sup = v_change_intervals[i];
   }
   for (; j < sup; ++j) {
    max_storage[j] = v_maximum_storage[i];
   }
  }
 }

 std::vector<double> max_primary_rho = v_maximum_primary_rho;
 if (max_primary_rho.size() == 1) {
  max_primary_rho.resize(f_time_horizon, max_primary_rho[0]);

 } else if (max_primary_rho.size() < f_time_horizon) {
  max_primary_rho.resize(f_time_horizon);
  int j = 0;
  for (unsigned long i = 0; i < v_change_intervals.size(); ++i) {
   Index sup;
   if (i == v_change_intervals.size() - 1 ) {
    sup = f_time_horizon;
   } else {
    sup = v_change_intervals[i];
   }
   for (; j < sup; ++j) {
    max_primary_rho[j] = v_maximum_primary_rho[i];
   }
  }
 }

 std::vector<double> max_secondary_rho = v_maximum_secondary_rho;
 if (max_secondary_rho.size() == 1) {
  max_secondary_rho.resize(f_time_horizon, max_secondary_rho[0]);

 } else if (max_secondary_rho.size() < f_time_horizon) {
  max_secondary_rho.resize(f_time_horizon);
  int j = 0;
  for (unsigned long i = 0; i < v_change_intervals.size(); ++i) {
   Index sup;
   if (i == v_change_intervals.size() - 1 ) {
    sup = f_time_horizon;
   } else {
    sup = v_change_intervals[i];
   }
   for (; j < sup; ++j) {
    max_secondary_rho[j] = v_maximum_secondary_rho[i];
   }
  }
 }

 std::vector<double> delta_ramp_up = this->v_delta_ramp_up;
 if( delta_ramp_up.size() == 1 ) {
  delta_ramp_up.resize( f_time_horizon, delta_ramp_up[0] );
 }else if (delta_ramp_up.size() < f_time_horizon) {
  delta_ramp_up.resize(f_time_horizon);
  int j = 0;
  for (unsigned long i = 0; i < v_change_intervals.size(); ++i) {
   Index sup;
   if (i == v_change_intervals.size() - 1 ) {
    sup = f_time_horizon;
   } else {
    sup = v_change_intervals[i];
   }
   for (; j < sup; ++j) {
    delta_ramp_up[j] = v_delta_ramp_up[i];
   }
  }
 }

 std::vector<double> delta_ramp_down = this->v_delta_ramp_down;
 if( delta_ramp_down.size() == 1 ) {
  delta_ramp_down.resize( f_time_horizon, delta_ramp_down[0] );
 }else if (delta_ramp_down.size() < f_time_horizon) {
  delta_ramp_down.resize(f_time_horizon);
  int j = 0;
  for (unsigned long i = 0; i < v_change_intervals.size(); ++i) {
   Index sup;
   if (i == v_change_intervals.size() - 1 ) {
    sup = f_time_horizon;
   } else {
    sup = v_change_intervals[i];
   }
   for (; j < sup; ++j) {
    delta_ramp_down[j] = v_delta_ramp_down[i];
   }
  }
 }


 std::vector<double> storing_battery_rho = v_storing_battery_rho;
 if (storing_battery_rho.size() == 1) {
  storing_battery_rho.resize(f_time_horizon, storing_battery_rho[0]);

 } else if (storing_battery_rho.size() < f_time_horizon) {
  storing_battery_rho.resize(f_time_horizon);
  int j = 0;
  for (unsigned long i = 0; i < v_change_intervals.size(); ++i) {
   Index sup;
   if (i == v_change_intervals.size() - 1 ) {
    sup = f_time_horizon;
   } else {
    sup = v_change_intervals[i];
   }
   for (; j < sup; ++j) {
    storing_battery_rho[j] = v_storing_battery_rho[i];
   }
  }
 }

 std::vector<double> extracting_battery_rho = v_extracting_battery_rho;
 if (extracting_battery_rho.size() == 1) {
  extracting_battery_rho.resize(f_time_horizon, extracting_battery_rho[0]);

 } else if (extracting_battery_rho.size() < f_time_horizon) {
  extracting_battery_rho.resize(f_time_horizon);
  int j = 0;
  for (unsigned long i = 0; i < v_change_intervals.size(); ++i) {
   Index sup;
   if (i == v_change_intervals.size() - 1 ) {
    sup = f_time_horizon;
   } else {
    sup = v_change_intervals[i];
   }
   for (; j < sup; ++j) {
    extracting_battery_rho[j] = v_extracting_battery_rho[i];
   }
  }
 }


 //TODO Demand
/*--------------------------------------------------------------------------*/

 // Initializing minimum power constraints
 if( !v_maximum_primary_rho.empty() & !v_maximum_secondary_rho.empty()) {
  active_power_lower_bound_Constraints.resize( f_time_horizon );

  for( Index t = 0; t < f_time_horizon; ++t ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_active_power[t][0], 1.0 );
   linear_function->add_variable( &v_primary_spinning_reserve[t][0], -1.0 );
   linear_function->add_variable( &v_secondary_spinning_reserve[t][0], -1.0 );

   active_power_lower_bound_Constraints[t].set_lhs( v_minimum_power[ t ]);
   active_power_lower_bound_Constraints[t].set_rhs( Inf< double >() );
   active_power_lower_bound_Constraints[t].set_function( linear_function );
  }
 } else {

  active_power_lower_bound_Constraints.resize( f_time_horizon );

  for( Index t = 0; t < f_time_horizon; ++t ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_active_power[t][0], 1.0 );

   active_power_lower_bound_Constraints[t].set_lhs( v_minimum_power[ t ] );
   active_power_lower_bound_Constraints[t].set_rhs( Inf< double >() );
   active_power_lower_bound_Constraints[t].set_function( linear_function );
  }
 }
 add_static_constraint( active_power_lower_bound_Constraints, "Active Power Lower Bound Constraints" );

 // Initializing maximum power constraints
 if( !v_maximum_primary_rho.empty() & !v_maximum_secondary_rho.empty()) {
  active_power_upper_bound_Constraints.resize( f_time_horizon );

  for( Index t = 0; t < f_time_horizon; ++t ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_active_power[t][0], 1.0 );
   linear_function->add_variable( &v_primary_spinning_reserve[ t ][ 0 ], 1.0 );
   linear_function->add_variable( &v_secondary_spinning_reserve[ t ][ 0 ], 1.0 );

   active_power_upper_bound_Constraints[t].set_lhs( 0.0 );
   active_power_upper_bound_Constraints[t].set_rhs( v_maximum_power[ t ] );
   active_power_upper_bound_Constraints[t].set_function( linear_function );
  }
 } else {

  active_power_upper_bound_Constraints.resize( f_time_horizon );

  for( Index t = 0; t < f_time_horizon; ++t ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_active_power[t][0], 1.0 );

   active_power_upper_bound_Constraints[t].set_lhs( 0.0 );
   active_power_upper_bound_Constraints[t].set_rhs( v_maximum_power[ t ] );
   active_power_upper_bound_Constraints[t].set_function( linear_function );
  }
 }
 add_static_constraint( active_power_upper_bound_Constraints, "Active Power Upper Bound Constraints" );

/*--------------------------------------------------------------------------*/
 // Initializing ramp-up constraints
 if( !v_delta_ramp_up.empty() & !v_delta_ramp_down.empty()) {
  {
   ramp_up_Constraints.resize( f_time_horizon );

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_active_power[0][0], 1.0 );

   ramp_up_Constraints[0].set_lhs( 0.0 );
   ramp_up_Constraints[0].set_rhs( delta_ramp_up[0] + f_initial_power );
   ramp_up_Constraints[0].set_function( linear_function );

   for( Index t = 1, constraint_index = 1; t < f_time_horizon; ++t, ++constraint_index ) {

    auto lf = new LinearFunction();

    lf->add_variable( &v_active_power[t][0], 1.0 );
    lf->add_variable( &v_active_power[t - 1][0], -1.0 );

    ramp_up_Constraints[constraint_index].set_lhs( 0.0 );
    ramp_up_Constraints[constraint_index].set_rhs( delta_ramp_up[t] );
    ramp_up_Constraints[constraint_index].set_function( lf );
   }
  }
  add_static_constraint( ramp_up_Constraints, "Ramp Up Constraints" );

  // Initializing ramp-down constraints

  {
  ramp_down_Constraints.resize( f_time_horizon );

  auto linear_function = new LinearFunction();

  linear_function->add_variable( &v_active_power[ 0 ][ 0 ], 1.0 );

   ramp_down_Constraints[ 0 ].set_lhs ( -delta_ramp_down[ 0 ] + f_initial_power );
   ramp_down_Constraints[ 0 ].set_rhs( Inf< double >() );
   ramp_down_Constraints[ 0 ].set_function( linear_function );

  for( Index t = 1, constraint_index = 1; t < f_time_horizon; ++t, ++constraint_index ) {

   auto lf = new LinearFunction();

   lf->add_variable( &v_active_power[t][0], 1.0 );
   lf->add_variable( &v_active_power[t - 1][0], -1.0 );

   ramp_down_Constraints[constraint_index].set_lhs( -delta_ramp_down[ 0 ] );
   ramp_down_Constraints[constraint_index].set_rhs( Inf< double >() );
   ramp_down_Constraints[constraint_index].set_function( lf );
  }
  }
  add_static_constraint( ramp_down_Constraints, "Ramp Down Constraints" );

 }
/*--------------------------------------------------------------------------*/
// Initializing power_intake_outtake_Constraints
 {
  power_intake_outtake_Constraints.resize( f_time_horizon );

  for( Index t = 1; t < f_time_horizon; ++t ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_active_power[t][0], 1.0 );
   linear_function->add_variable( &v_intake_level[t], -1.0 );
   linear_function->add_variable( &v_outtake_level[t], 1.0 );

   power_intake_outtake_Constraints[t].set_both( 0.0);
   power_intake_outtake_Constraints[t].set_function( linear_function );

  }
 }

 add_static_constraint( power_intake_outtake_Constraints, "Power_Intake_Outtake Constraints" );


// Initializing intake_upper_bound_Constraints

 {
  intake_upper_bound_Constraints.resize( f_time_horizon );

  for( Index t = 1; t < f_time_horizon; ++t ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_intake_level[t], 1.0 );

   intake_upper_bound_Constraints[t].set_lhs( 0.0);
   intake_upper_bound_Constraints[t].set_rhs( v_maximum_power[ t ]);
   intake_upper_bound_Constraints[t].set_function( linear_function );

  }
 }
 add_static_constraint( intake_upper_bound_Constraints, "Intake UpperBound Constraints" );

 /*--------------------------------------------------------------------------*/
//TODO DEMAND CONSTRAINT

 /*--------------------------------------------------------------------------*/
// Initializing storage_level_bounds_Constraints

 {
  storage_level_bounds_Constraints.resize( f_time_horizon );

  for( Index t = 1; t < f_time_horizon; ++t ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_storage_level[t], 1.0 );

   storage_level_bounds_Constraints[t].set_lhs( v_minimum_storage[ t ]);
   storage_level_bounds_Constraints[t].set_rhs( v_maximum_storage[ t ]);
   storage_level_bounds_Constraints[t].set_function( linear_function );

  }
 }
 add_static_constraint( storage_level_bounds_Constraints, "Storage Level Bounds Constraints" );

/*--------------------------------------------------------------------------*/

 // Initializing intake_binary_Constraints

 {
  intake_binary_Constraints.resize( f_time_horizon );

  for( Index t = 1; t < f_time_horizon; ++t ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_intake_level[t], 1.0 );
   linear_function->add_variable( &v_battery_binary[t], - v_maximum_power[ t ] );


   intake_binary_Constraints[t].set_both( 0.0);
   intake_binary_Constraints[t].set_function( linear_function );

  }
 }
 add_static_constraint( intake_binary_Constraints, "Intake Binary Constraints" );

/*--------------------------------------------------------------------------*/

 // Initializing outtake_binary_Constraints

 {
  outtake_binary_Constraints.resize( f_time_horizon );

  for( Index t = 1; t < f_time_horizon; ++t ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_outtake_level[ t ], 1.0 );
   linear_function->add_variable( &v_battery_binary[ t ],  v_minimum_power[ t ] );


   outtake_binary_Constraints[t].set_lhs( 0.0);
   outtake_binary_Constraints[t].set_rhs( -v_minimum_power[ t ] );
   outtake_binary_Constraints[t].set_function( linear_function );

  }
 }
 add_static_constraint( outtake_binary_Constraints, "Outtake Binary Constraints" );
/*--------------------------------------------------------------------------*/

 if( !v_maximum_primary_rho.empty() & !v_maximum_secondary_rho.empty()) {

  // Initializing primary_upper_bound_Constraints

  primary_upper_bound_Constraints.resize( f_time_horizon );

  for( Index t = 1; t < f_time_horizon; ++t ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_primary_spinning_reserve[ t ][ 0 ], 1.0 );


   primary_upper_bound_Constraints[t].set_lhs( 0.0);
   primary_upper_bound_Constraints[t].set_rhs( v_maximum_primary_rho[ t ] );
   primary_upper_bound_Constraints[t].set_function( linear_function );

  }
  add_static_constraint( primary_upper_bound_Constraints, "Primary Upper Bound Constraints" );

  // Initializing secondary_upper_bound_Constraints

  secondary_upper_bound_Constraints.resize( f_time_horizon );

  for( Index t = 1; t < f_time_horizon; ++t ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_secondary_spinning_reserve[ t ][ 0 ], 1.0 );


   secondary_upper_bound_Constraints[t].set_lhs( 0.0);
   secondary_upper_bound_Constraints[t].set_rhs( v_maximum_secondary_rho[ t ] );
   secondary_upper_bound_Constraints[t].set_function( linear_function );

  }
  add_static_constraint( secondary_upper_bound_Constraints, "Secondary Upper Bound Constraints" );
 }

 } // end( BatteryUnitBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::generate_objective( Configuration *objc )
{
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

 if( get_objective() == nullptr )  // an objective is there already
  return;                         // cowardly (and silently) return


 if( v_intake_level.size() != f_time_horizon  || v_outtake_level.size() != f_time_horizon) {
  throw ( std::logic_error
          ( "BatteryUnitBlock::generate_objective: v_intake_level and  must v_outtake_level have "
            "size equal to the time horizon." ));
 }
  auto dquad_function = new DQuadFunction();

  for( Index t = 0; t < f_time_horizon; ++t ) {
   dquad_function->add_variable( &v_intake_level[ t ],
                                 v_cost[ t  ],
                                 0.0 );

   dquad_function->add_variable( &v_outtake_level[ t ],
                                 v_cost[ t  ],
                                 0.0 );
  }

 objective.set_function( dquad_function );
 objective.set_sense( Objective::eMin );

 // Set Block objective
 this->set_objective( &objective );

}  // end( BatteryUnitBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*--- METHODS FOR LOADING, PRINTING & SAVING THE BatteryUnitBlock ---*/
/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::serialize( netCDF::NcGroup & group ) const {

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

 ::serialize( group, "ExtractingBatterRho", netCDF::NcDouble(),
              NumberIntervals, v_extracting_battery_rho, true );

 ::serialize( group, "Cost", netCDF::NcDouble(),
              NumberIntervals, v_cost, true );

 ::serialize( group, "Demand", netCDF::NcDouble(),
              TimeHorizon, v_cost, false );
}  // end( BatteryUnitBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------- End File BatteryUnitBlock.cpp -----------------------*/
/*--------------------------------------------------------------------------*/