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
#include <map>

#include "DQuadFunction.h"
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
/*------------------- METHODS OF BatteryUnitBlock -------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
void BatteryUnitBlock::deserialize( netCDF::NcGroup & group ) {

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
 ::deserialize( group, "StoringBatteryrho",f_time_horizon, v_storing_battery_rho, true,true );
 ::deserialize( group, "ExtractingBatteryrho",f_time_horizon, v_extracting_battery_rho, true,true );
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

void BatteryUnitBlock::generate_abstract_variables  ( Configuration *stvv ) {

 UnitBlock::generate_abstract_variables( stvv );

 auto var_size = f_time_horizon ;

 if( var_size > 0 ) {

  if( v_storage_level.size() != var_size ) {
   assert( v_storage_level.empty()); // this should only happen once
   v_storage_level.resize( var_size );
   int n = 0;
   for( auto & i : v_storage_level ) {
    i.set_type( ColVariable::kNonNegative );
   }
   add_static_variable( v_storage_level, "SL" );
  }

  if( v_intake_level.size() != var_size ) {
   assert( v_intake_level.empty()); // this should only happen once
   v_intake_level.resize( var_size );
   int n = 0;
   for( auto & i : v_intake_level ) {
    i.set_type( ColVariable::kNonNegative );
   }
   add_static_variable( v_intake_level, "IL" );
  }

  if( v_outtake_level.size() != var_size ) {
   assert( v_outtake_level.empty()); // this should only happen once
   v_outtake_level.resize( var_size );
   int n = 0;
   for( auto & i : v_outtake_level ) {
    i.set_type( ColVariable::kNonNegative );
   }
   add_static_variable( v_outtake_level, "OL" );
  }

  if( v_battery_binary.size() != var_size ) {
   assert( v_battery_binary.empty()); // this should only happen once
   v_battery_binary.resize( var_size );
   int n = 0;
   for( auto & i : v_battery_binary ) {
    i.set_type( ColVariable::kBinary );
   }
   add_static_variable( v_battery_binary, "BB" );

  }

  // Active Power Variable

  if( v_active_power.size() != f_time_horizon ) {
   assert( v_active_power.empty() ); // this should only happen once
   v_active_power.resize( f_time_horizon );
   int n = 0;
   for( auto & i : v_active_power ) {
    i.set_type( ColVariable::kContinuous );
   }
   add_static_variable( v_active_power, "p" );
  }

  // Primary Spinning Reserve Variable

  if( v_primary_spinning_reserve.size() != f_time_horizon ) {
   assert( v_primary_spinning_reserve.empty() ); // this should only happen once
   v_primary_spinning_reserve.resize( f_time_horizon );
   int n = 0;
   for( auto & i : v_primary_spinning_reserve ) {
    i.set_type( ColVariable::kNonNegative );
   }
   add_static_variable( v_primary_spinning_reserve, "pr"  );

  }

  // Secondary Spinning Reserve Variable

  if( v_secondary_spinning_reserve.size() != f_time_horizon ) {
   assert( v_secondary_spinning_reserve.empty() ); // this should only happen once
   v_secondary_spinning_reserve.resize( f_time_horizon );
   int n = 0;
   for( auto & i : v_secondary_spinning_reserve ) {
    i.set_type( ColVariable::kNonNegative );
   }
   add_static_variable( v_secondary_spinning_reserve, "sr"  );

  }

 }
} // end( BatteryUnitBlock::generate_abstract_variables )
/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::generate_abstract_constraints ( Configuration * stcc ) {



 // Initial data check
   for( Index t = 0; t < f_time_horizon; ++t ) {
    if( v_minimum_power[t] >= v_maximum_power[t] ) {
     throw ( std::logic_error
             ( "BatteryUnitBlock::maximum and minimum power output constraints: "
               "it must be that v_maximum_power > v_minimum_power." ));
    }
   }


 for( Index t = 0; t < f_time_horizon; ++t ) {
  if( v_minimum_storage[t] >= v_maximum_storage[t] ||
          v_minimum_storage[t] < 0  || v_maximum_storage[t] < 0 ) {
   throw ( std::logic_error
           ( "BatteryUnitBlock::maximum and minimum storage output constraints: "
             "it must be that v_maximum_storage > v_minimum_storage >= 0." ));
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
   if( !v_maximum_primary_rho.empty()) {
    linear_function->add_variable( &v_primary_spinning_reserve[t], -1.0 );
    if( !v_maximum_secondary_rho.empty()) {
     linear_function->add_variable( &v_secondary_spinning_reserve[t], -1.0 );
    }
    active_power_lower_bound_Constraints[t].set_lhs( v_minimum_power[t] );
    active_power_lower_bound_Constraints[t].set_rhs( Inf< double >());
    active_power_lower_bound_Constraints[t].set_function( linear_function );
   }
  }
  add_static_constraint( active_power_lower_bound_Constraints, "Active Power Lower Bound Constraints" );

  // Initializing maximum power constraints

   active_power_upper_bound_Constraints.resize( f_time_horizon );

   for( Index t = 0; t < f_time_horizon; ++t ) {

    auto linear_function = new LinearFunction();

    linear_function->add_variable( &v_active_power[t], 1.0 );
    if( !v_maximum_primary_rho.empty()) {
     linear_function->add_variable( &v_primary_spinning_reserve[t], 1.0 );
    }
    if (!v_maximum_secondary_rho.empty()) {
     linear_function->add_variable( &v_secondary_spinning_reserve[t], 1.0 );
    }
    active_power_upper_bound_Constraints[t].set_lhs( -Inf< double >());
    active_power_upper_bound_Constraints[t].set_rhs( v_maximum_power[t] );
    active_power_upper_bound_Constraints[t].set_function( linear_function );
   }

  add_static_constraint( active_power_upper_bound_Constraints, "Active Power Upper Bound Constraints" );

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

 add_static_constraint( ramp_up_Constraints, "Ramp Up Constraints" );

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
 add_static_constraint( ramp_down_Constraints, "Ramp Down Constraints" );

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

 add_static_constraint( power_intake_outtake_Constraints, "Power_Intake_Outtake Constraints" );


// Initializing intake_upper_bound_Constraints


  intake_upper_bound_Constraints.resize( f_time_horizon );

  for( Index t = 0; t < f_time_horizon; ++t ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_intake_level[t], 1.0 );

   intake_upper_bound_Constraints[t].set_lhs( 0.0);
   intake_upper_bound_Constraints[t].set_rhs( v_maximum_power[ t ]);
   intake_upper_bound_Constraints[t].set_function( linear_function );

  }
 add_static_constraint( intake_upper_bound_Constraints, "Intake UpperBound Constraints" );

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

 add_static_constraint( demand_Constraints, "demand Constraints" );

 /*--------------------------------------------------------------------------*/
// Initializing storage_level_bounds_Constraints

  storage_level_bounds_Constraints.resize( f_time_horizon );

  for( Index t = 0; t < f_time_horizon; ++t ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_storage_level[t], 1.0 );

   storage_level_bounds_Constraints[t].set_lhs( v_minimum_storage[ t ]);
   storage_level_bounds_Constraints[t].set_rhs( v_maximum_storage[ t ]);
   storage_level_bounds_Constraints[t].set_function( linear_function );

  }
 add_static_constraint( storage_level_bounds_Constraints, "Storage Level Bounds Constraints" );

/*--------------------------------------------------------------------------*/

 // Initializing intake_binary_Constraints

  if (!v_storing_battery_rho.empty() & !v_extracting_battery_rho.empty()){

  intake_binary_Constraints.resize( f_time_horizon );

  for( Index t = 0; t < f_time_horizon; ++t ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_intake_level[t], 1.0 );
   linear_function->add_variable( &v_battery_binary[t],  -v_maximum_power[ t ] );


   intake_binary_Constraints[t].set_lhs( -Inf<double>());
   intake_binary_Constraints[t].set_rhs( 0.0 );
   intake_binary_Constraints[t].set_function( linear_function );

  }
 }
 add_static_constraint( intake_binary_Constraints, "Intake Binary Constraints" );

/*--------------------------------------------------------------------------*/

 // Initializing outtake_binary_Constraints

 if (!v_storing_battery_rho.empty() & !v_extracting_battery_rho.empty()){
  outtake_binary_Constraints.resize( f_time_horizon );

  for( Index t = 0; t < f_time_horizon; ++t ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_outtake_level[ t ], 1.0 );
   linear_function->add_variable( &v_battery_binary[ t ],  - v_minimum_power[ t ] );


   outtake_binary_Constraints[t].set_lhs( -Inf<double>());
   outtake_binary_Constraints[t].set_rhs( -v_minimum_power[ t ] );
   outtake_binary_Constraints[t].set_function( linear_function );

  }
 }
 add_static_constraint( outtake_binary_Constraints, "Outtake Binary Constraints" );
/*--------------------------------------------------------------------------*/

  // Initializing primary_upper_bound_Constraints
 if ( !v_maximum_primary_rho.empty() ) {
  primary_upper_bound_Constraints.resize( f_time_horizon );

  for( Index t = 0; t < f_time_horizon; ++t ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_primary_spinning_reserve[t], 1.0 );


   primary_upper_bound_Constraints[t].set_lhs( 0.0 );
   primary_upper_bound_Constraints[t].set_rhs( v_maximum_primary_rho[t] );
   primary_upper_bound_Constraints[t].set_function( linear_function );

  }
  add_static_constraint( primary_upper_bound_Constraints, "Primary Upper Bound Constraints" );
 }

  // Initializing secondary_upper_bound_Constraints
 if ( !v_secondary_spinning_reserve.empty() ) {

 secondary_upper_bound_Constraints.resize( f_time_horizon );

  for( Index t = 0; t < f_time_horizon; ++t ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_secondary_spinning_reserve[ t ], 1.0 );


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

 if( get_objective() != nullptr )  // an objective is there already
  return;                         // cowardly (and silently) return


 if( v_intake_level.size() != f_time_horizon  ) {
  throw ( std::logic_error
          ( "BatteryUnitBlock::generate_objective: v_intake_level and  must v_intake_level have "
            "size equal to the time horizon." ));
 }

 if(  v_outtake_level.size() != f_time_horizon) {
  throw ( std::logic_error
          ( "BatteryUnitBlock::generate_objective: v_intake_level and  must v_outtake_level have "
            "size equal to the time horizon." ));
 }
  auto dquad_function = new DQuadFunction();

  for( Index t = 0; t < f_time_horizon; ++t ) {
   dquad_function->add_variable( &v_intake_level[ t ],
                                 cost[ t  ],
                                 0.0 );

   dquad_function->add_variable( &v_outtake_level[ t ],
                                 cost[ t  ],
                                 0.0 );
  }

 objective.set_function( dquad_function );
 objective.set_sense( Objective::eMin );

 // Set Block objective
 this->set_objective( &objective );
 AR |= HasObj;

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

template< typename T >
void BatteryUnitBlock::decompress_vector( std::vector< T > & v ) {
 if (v.empty()) {
  return;
 }

 if( v.size() == 1 ) {
  v.resize( f_time_horizon, v[ 0 ] );
 } else if( v.size() < f_time_horizon ) {
  std::vector< T > temp = v;
  v.resize( f_time_horizon );
  int j = 0;
  for( unsigned long i = 0; i < v_change_intervals.size(); ++i ) {
   Index sup;
   if( i == v_change_intervals.size() - 1 ) {
    sup = f_time_horizon;
   } else {
    sup = v_change_intervals[ i ];
   }
   for( ; j < sup; ++j ) {
    v[ j ] = temp[ i ];
   }
  }
 }
}

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/
void
BatteryUnitBlock::set_initial_storage( std::vector< double >::const_iterator it,
                                       Block::Subset && subset,
                                       const bool ordered,
                                       c_ModParam issuePMod,
                                       c_ModParam issueAMod ) {
 // TODO PUT STUFF HERE
}

void
BatteryUnitBlock::set_initial_storage( std::vector< double >::const_iterator it,
                                       Block::Range rng,
                                       c_ModParam issuePMod,
                                       c_ModParam issueAMod ) {
 // TODO PUT STUFF HERE
}

/*--------------------------------------------------------------------------*/
/*------------- End File BatteryUnitBlock.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
