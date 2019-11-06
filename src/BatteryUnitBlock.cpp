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
    i.set_type( ColVariable::kBinary );
    //TODO look ThermalUnitBlock variables
    add_static_variable( i, "IntakeLevel " + std::to_string( n++ ));
   }
  }

  if( v_outtake_level.size() != f_time_horizon ) {
   assert( v_outtake_level.empty()); // this should only happen once
   v_outtake_level.resize( f_time_horizon );
   int n = 0;
   for( auto & i : v_outtake_level ) {
    i.set_type( ColVariable::kBinary );
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


 //TODO Cost
/*--------------------------------------------------------------------------*/





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
 if( get_objective() == nullptr )  // an objective is there already
  return;                         // cowardly (and silently) return

// initialize objective function - - - - - - - - - - - - - - - - - - - - - -
 //TODO
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