/*--------------------------------------------------------------------------*/
/*--------------------- File ThermalUnitBlock.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the ThermalUnitBlock class.
 *
 * \version 0.11
 *
 * \date 19 - 06 - 2019
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Ali Ghezelsoflu \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Rafael Durbano Lobato \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Antonio Frangioni, Ali Ghezelsoflu, and Rafael
 * Durbano Lobato
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
#include "ThermalUnitBlock.h"
#include "UnitBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;


/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register ThermalUnitBlock to the Block factory

SMSpp_insert_in_factory_cpp_1( ThermalUnitBlock );

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS OF ThermalUnitBlock ----------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::deserialize( netCDF::NcGroup & group ) {

 UnitBlock::deserialize( group );

 ::deserialize( group, "MinPower", f_number_intervals, v_MinPower, true, true );

 ::deserialize( group, "MaxPower", f_number_intervals, v_MaxPower, true, true );

 ::deserialize( group, "DeltaRampUp", f_number_intervals, v_DeltaRampUp, true, true );

 ::deserialize( group, "DeltaRampDown", f_number_intervals, v_DeltaRampDown, true, true );

 ::deserialize( group, "PrimaryRho", f_number_intervals, v_PrimaryRho, true, true );

 ::deserialize( group, "SecondaryRho", f_number_intervals, v_SecondaryRho, true, true );

 ::deserialize( group, "LinearTerm", f_number_intervals, v_LinearTerm, true, true );

 ::deserialize( group, "QuadTerm", f_number_intervals, v_QuadTerm, true, true );

 ::deserialize( group, "ConstTerm", f_number_intervals, v_ConstTerm, true, true );

 ::deserialize( group, "StartUpCost", f_number_intervals, v_StartUpCost, true, true );

 ::deserialize( group, "FixedConsumption", v_fixed_consumption, true, true );
 ::deserialize( group, "InertiaCommitment", v_inertia_commitment, true, true );

 if (v_fixed_consumption.empty()) {
  v_fixed_consumption.resize( boost::extents[ f_time_horizon ][ 1 ] );
 }
 if (v_inertia_commitment.empty()) {
  v_inertia_commitment.resize( boost::extents[ f_time_horizon ][ 1 ] );
 }

 ::deserialize( group, "InitialPower", &f_initial_power );
 ::deserialize( group, "MinUpTime", &f_MinUpTime );
 ::deserialize( group, "MinDownTime", &f_MinDownTime );
 ::deserialize( group, "InitUpDownTime", &f_InitUpDownTime );

}  // end( ThermalUnitBlock::deserialize )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::generate_abstract_variables( Configuration * stvv ) {

 UnitBlock::generate_abstract_variables( stvv );

 if( f_InitUpDownTime > 0 ) {
  init_t = ( f_InitUpDownTime >= f_MinUpTime ? 0 :
             f_MinUpTime - f_InitUpDownTime );
 } else {
  init_t = ( -f_InitUpDownTime >= f_MinDownTime ? 0 :
             f_MinDownTime + f_InitUpDownTime );
 }
/*--------------------------------------------------------------------------*/
if ( f_time_horizon > 0 ){

 // Commitment Variable

 if( v_commitment.size() != f_time_horizon ) {
  assert( v_commitment.empty() ); // this should only happen once
  v_commitment.resize( f_time_horizon );
  int n = 0;
  for( auto & i : v_commitment ) {
   i.set_type( ColVariable::kBinary );
   add_static_variable( i, "u_" + std::to_string( n++ ) );
  }
 }

 // Active Power Variable

 if( v_active_power.size() != f_time_horizon ) {
  assert( v_active_power.empty() ); // this should only happen once
  v_active_power.resize( f_time_horizon );
  int n = 0;
  for( auto & i : v_active_power ) {
   i.set_type( ColVariable::kNonNegative );
   add_static_variable( i, "p_" + std::to_string( n++ ) );
  }
 }

 // Primary Spinning Reserve Variable

 if( v_primary_spinning_reserve.size() != f_time_horizon ) {
  assert( v_primary_spinning_reserve.empty() ); // this should only happen once
  v_primary_spinning_reserve.resize( f_time_horizon );
  int n = 0;
  for( auto & i : v_primary_spinning_reserve ) {
   i.set_type( ColVariable::kNonNegative );
   add_static_variable( i, "pr_" + std::to_string( n++ ) );
  }
 }

 // Secondary Spinning Reserve Variable

 if( v_secondary_spinning_reserve.size() != f_time_horizon ) {
  assert( v_secondary_spinning_reserve.empty() ); // this should only happen once
  v_secondary_spinning_reserve.resize( f_time_horizon );
  int n = 0;
  for( auto & i : v_secondary_spinning_reserve ) {
   i.set_type( ColVariable::kNonNegative );
   add_static_variable( i, "sr_" + std::to_string( n++ ) );
  }
 }

}

/*--------------------------------------------------------------------------*/
 // START UP AND SHUT DOWN BINARY VARIABLES
 auto startup_shutdown_size = f_time_horizon - init_t;

 if( startup_shutdown_size > 0 ) {

  if( v_start_up.size() != startup_shutdown_size ) {
   assert( v_start_up.empty() ); // this should only happen once
   v_start_up.resize( startup_shutdown_size );
   int n = 0;
   for( auto & i : v_start_up ) {
    i.set_type( ColVariable::kBinary );
    add_static_variable( i, "v_" + std::to_string( n++ ) );
   }
   // add_static_variable( v_start_up, "v" );
  }

  if( v_shut_down.size() != startup_shutdown_size ) {
   assert( v_shut_down.empty() ); // this should only happen once
   v_shut_down.resize( startup_shutdown_size );
   int n = 0;
   for( auto & i : v_shut_down ) {
    i.set_type( ColVariable::kBinary );
    add_static_variable( i, "w_" + std::to_string( n++ ) );
   }
   // add_static_variable( v_shut_down, "w" );
  }
 }

/*--------------------------------------------------------------------------*/
 // POSSIBLY FIXING THE COMMITMENT VARIABLES TO 0 OR 1

 if( init_t > 0 ) {

  double commitment_variable_value = -1.0;

  if( f_InitUpDownTime < 0 && -f_InitUpDownTime < f_MinDownTime ) {

   commitment_variable_value = 0.0;

   // This unit must remain off from time 0 to init_t - 1. Therefore,
   // it should not produce any power.

   if( !v_active_power.empty()) {
    for( Index t = 0; t < init_t; ++t ) {
     v_active_power[t].set_value( 0.0 );
     v_active_power[t].is_fixed( true );
    }
   }

   if( !v_primary_spinning_reserve.empty()) {
    for( Index t = 0; t < init_t; ++t ) {
     v_primary_spinning_reserve[t].set_value( 0.0 );
     v_primary_spinning_reserve[t].is_fixed( true );
    }
   }

   if( !v_secondary_spinning_reserve.empty()) {
    for( Index t = 0; t < init_t; ++t ) {
     v_secondary_spinning_reserve[t].set_value( 0.0 );
     v_secondary_spinning_reserve[t].is_fixed( true );
    }
   }
  } else if( f_InitUpDownTime > 0 && f_InitUpDownTime < f_MinUpTime ) {
   // This unit must remain on from time 0 to init_t - 1.
   commitment_variable_value = 1.0;

  }

  if( !v_commitment.empty()) {
   if( commitment_variable_value >= 0.0 ) {
    // Fixing the commitment variables v_commitment to 0 or 1 from
    // time 0 to init_t - 1, depending on whether this unit must
    // remain off or on for the first init_t time steps.
    for( Index t = 0; t < init_t; ++t ) {
     v_commitment[t].set_value( commitment_variable_value );
     v_commitment[t].is_fixed( true );
    }
   }
  }

  if (f_InitUpDownTime > 0){
   for ( Index t = 0; t < f_MinDownTime ; ++t) {
    v_start_up[t].set_value( 0.0 );
    v_start_up[t].is_fixed( true );
   }
  } else if (f_InitUpDownTime < 0){
   for ( Index t = 0; t < f_MinUpTime  ; ++t) {
    v_shut_down[t].set_value( 0.0 );
    v_shut_down[t].is_fixed( true );
   }
  }
 }

 if (init_t == 0){

  if (f_InitUpDownTime > 0){

   for ( Index t = 0; t < f_MinDownTime ; ++t) {
    v_start_up[t].set_value( 0.0 );
    v_start_up[t].is_fixed( true );

   }
  } else if (f_InitUpDownTime < 0){

   for ( Index t = 0; t < f_MinUpTime; ++t) {
    v_shut_down[t].set_value( 0.0 );
    v_shut_down[t].is_fixed( true );
   }
  }
 }
} // end( ThermalUnitBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::generate_abstract_constraints( Configuration * stcc ) {

 // initial condition of each vector
 std::vector<double> min_power = v_MinPower;
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
    min_power[j] = v_MinPower[i];
   }
  }
 }

 std::vector<double> max_power = this->v_MaxPower;
 if( max_power.size() == 1 ) {
  max_power.resize( f_time_horizon, max_power[0] );
 }else if (max_power.size() < f_time_horizon) {
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
    max_power[j] = v_MaxPower[i];
   }
  }
 }


 std::vector<double> delta_ramp_up = this->v_DeltaRampUp;
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
    delta_ramp_up[j] = v_DeltaRampUp[i];
   }
  }
 }

 std::vector<double> delta_ramp_down = this->v_DeltaRampDown;
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
    delta_ramp_down[j] = v_DeltaRampDown[i];
   }
  }
 }

 std::vector<double> primary_rho = this->v_PrimaryRho;
 if( primary_rho.size() == 1 ) {
  primary_rho.resize( f_time_horizon, primary_rho[0] );
 }else if (primary_rho.size() < f_time_horizon) {
  primary_rho.resize(f_time_horizon);
  int j = 0;
  for (unsigned long i = 0; i < v_change_intervals.size(); ++i) {
   Index sup;
   if (i == v_change_intervals.size() - 1 ) {
    sup = f_time_horizon;
   } else {
    sup = v_change_intervals[i];
   }
   for (; j < sup; ++j) {
    primary_rho[j] = v_PrimaryRho[i];
   }
  }
 }

 std::vector<double> secondary_rho = this->v_SecondaryRho;
 if( secondary_rho.size() == 1 ) {
  secondary_rho.resize( f_time_horizon, secondary_rho[0] );
 }else if (secondary_rho.size() < f_time_horizon) {
  secondary_rho.resize(f_time_horizon);
  int j = 0;
  for (unsigned long i = 0; i < v_change_intervals.size(); ++i) {
   Index sup;
   if (i == v_change_intervals.size() - 1 ) {
    sup = f_time_horizon;
   } else {
    sup = v_change_intervals[i];
   }
   for (; j < sup; ++j) {
    secondary_rho[j] = v_SecondaryRho[i];
   }
  }
 }
/*--------------------------------------------------------------------------*/

 // Initializing start up and shut down variables connection constraints
 if( init_t == 0 ) {

  if( f_time_horizon > 0 ) {

   StartUp_ShutDown_Variables_Constraints.resize( f_time_horizon );

   // Initial condition

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_commitment[0], 1.0 );
   linear_function->add_variable( &v_start_up[0], -1.0 );
   linear_function->add_variable( &v_shut_down[0], 1.0 );

   auto initial_commitment = ( f_InitUpDownTime > 0 ? 1.0 : 0.0 );

   StartUp_ShutDown_Variables_Constraints[0].set_both( initial_commitment );
   StartUp_ShutDown_Variables_Constraints[0].set_function( linear_function );


   for( Index t = 1, constraint_index = 1; t < f_time_horizon;
        ++t, ++constraint_index ) {

    auto lf = new LinearFunction();

    lf->add_variable( &v_commitment[t], 1.0 );
    lf->add_variable( &v_start_up[t], -1.0 );
    lf->add_variable( &v_shut_down[t], 1.0 );
    lf->add_variable( &v_commitment[t - 1], -1.0 );
    StartUp_ShutDown_Variables_Constraints[constraint_index].set_both( 0.0 );

    StartUp_ShutDown_Variables_Constraints[constraint_index].
            set_function( lf );
   }
  }
 }

 if( init_t > 0 ) {

  auto startup_shutdown_const_size = static_cast<int>(f_time_horizon - init_t );

  if( startup_shutdown_const_size > 0 ) {

   StartUp_ShutDown_Variables_Constraints.resize( startup_shutdown_const_size );

   // Initial condition

   auto l_function = new LinearFunction();

   l_function->add_variable( &v_commitment[init_t], 1.0 );
   l_function->add_variable( &v_start_up[0], -1.0 );
   l_function->add_variable( &v_shut_down[0], 1.0 );

   if( f_InitUpDownTime < 0 && -f_InitUpDownTime < f_MinDownTime ) {

    StartUp_ShutDown_Variables_Constraints[0].set_both( 0.0 );
    StartUp_ShutDown_Variables_Constraints[0].set_function( l_function );

   } else if( f_InitUpDownTime > 0 && f_InitUpDownTime < f_MinUpTime ) {

    StartUp_ShutDown_Variables_Constraints[0].set_both( 1.0 );
    StartUp_ShutDown_Variables_Constraints[0].set_function( l_function );
   }

   for( Index t = init_t + 1 , constraint_index = 1; t < f_time_horizon;
        ++t, ++constraint_index ) {

    auto linear_function = new LinearFunction();

    linear_function->add_variable( &v_commitment[t], 1.0 );
    linear_function->add_variable( &v_start_up[t - init_t], -1.0 );
    linear_function->add_variable( &v_shut_down[t - init_t], 1.0 );
    linear_function->add_variable( &v_commitment[t - 1], -1.0 );
    StartUp_ShutDown_Variables_Constraints[constraint_index].set_both( 0.0 );

    StartUp_ShutDown_Variables_Constraints[constraint_index].
            set_function( linear_function );
   }
  }
 }
 add_static_constraint( StartUp_ShutDown_Variables_Constraints,
                        "startup_shutdown_vars_c" );


 // Initializing turn on constraints (start up constraints)

 auto startup_const_size = static_cast<int>(f_time_horizon - init_t - f_MinUpTime);
 if( startup_const_size > 0 ) {

  StartUp_Constraints.resize( startup_const_size );

  for( Index t = init_t + f_MinUpTime, constraint_index = 0;
       t < f_time_horizon; ++t, ++constraint_index ) {

   auto linear_function = new LinearFunction();

   for( Index s = t - f_MinUpTime; s < t; ++s ) {
    linear_function->add_variable( &v_start_up[s - init_t + 1], -1.0 );
   }

   linear_function->add_variable( &v_commitment[t], 1.0 );
   StartUp_Constraints[constraint_index].set_lhs( 0.0 );
   StartUp_Constraints[constraint_index].set_rhs( Inf< double >());
   StartUp_Constraints[constraint_index].set_function( linear_function );
  }

  add_static_constraint( StartUp_Constraints, "startup_c" );
 }

 // Initializing turn off constraints (shut down constraints)

 auto shutdown_const_size = static_cast<int>(f_time_horizon - init_t - f_MinDownTime);
 if( shutdown_const_size > 0 ) {

  ShutDown_Constraints.resize( shutdown_const_size );

  for( Index t = init_t + f_MinDownTime, constraint_index = 0;
       t < f_time_horizon; ++t, ++constraint_index ) {

   auto linear_function = new LinearFunction();

   for( Index s = t - f_MinDownTime; s < t; ++s ) {
    linear_function->add_variable( &v_shut_down[s - init_t + 1], 1.0 );
   }

   linear_function->add_variable( &v_commitment[t], 1.0 );
   ShutDown_Constraints[constraint_index].set_lhs( 0.0 );
   ShutDown_Constraints[constraint_index].set_rhs( 1.0 );
   ShutDown_Constraints[constraint_index].set_function( linear_function );
  }

  add_static_constraint( ShutDown_Constraints, "shutdown_c" );
 }

/*--------------------------------------------------------------------------*/
 // Initializing ramp-up constraints with 3-Binary Variables

 if( !v_DeltaRampUp.empty() & !v_DeltaRampDown.empty()) {

  RampUp_Constraints.resize( f_time_horizon );

  // Initial condition


  if( init_t == 0 ) {

   // Initial condition

   auto linear_function = new LinearFunction();
   linear_function->add_variable( &v_active_power[0], 1.0 );
   linear_function->add_variable( &v_start_up[0], -min_power[0] );
   auto initial_commitment = ( f_InitUpDownTime > 0 ? 1.0 : 0.0 );

   RampUp_Constraints[0].set_lhs
           ( 0.0 );
   RampUp_Constraints[0].set_rhs( (delta_ramp_up[0] * initial_commitment) + f_initial_power );
   RampUp_Constraints[0].set_function( linear_function );

   for( Index t = 1, constraint_index = 1; t < f_time_horizon; ++t, ++constraint_index ) {

    auto lf = new LinearFunction();

    lf->add_variable( &v_active_power[t], -1.0 );
    lf->add_variable( &v_active_power[t - 1], 1.0 );
    lf->add_variable( &v_start_up[t], min_power[t] );

    lf->add_variable( &v_commitment[t - 1], ( delta_ramp_up[t] ));

    RampUp_Constraints[constraint_index].set_lhs( 0.0 );
    RampUp_Constraints[constraint_index].set_rhs( Inf< double >());
    RampUp_Constraints[constraint_index].set_function( lf );
   }

  }

  if( init_t > 0 ) {

   auto linear_function = new LinearFunction();
   linear_function->add_variable( &v_active_power[0], 1.0 );

   auto initial_commitment = ( f_InitUpDownTime > 0 ? 1.0 : 0.0 );

   RampUp_Constraints[0].set_lhs
           ( 0.0 );
   RampUp_Constraints[0].set_rhs( (delta_ramp_up[0] * initial_commitment) + f_initial_power );
   RampUp_Constraints[0].set_function( linear_function );

   for( Index t = 1, constraint_index = 1; t < init_t;
        ++t, ++constraint_index ) {

    auto lf = new LinearFunction();

    lf->add_variable( &v_active_power[t], -1.0 );
    lf->add_variable( &v_active_power[t - 1], 1.0 );
    lf->add_variable
            ( &v_commitment[t - 1], ( delta_ramp_up[t] ));

    RampUp_Constraints[constraint_index].set_lhs( 0.00 );
    RampUp_Constraints[constraint_index].set_rhs( Inf< double >());
    RampUp_Constraints[constraint_index].set_function( lf );
   }


   if( f_InitUpDownTime < 0 && -f_InitUpDownTime < f_MinDownTime ) {

    auto LFunction = new LinearFunction();
    LFunction->add_variable( &v_active_power[init_t], 1.0 );
    LFunction->add_variable( &v_start_up[0], -min_power[init_t] );

    RampUp_Constraints[init_t].set_lhs( -Inf< double >() );
    RampUp_Constraints[init_t].set_rhs( 0.0 );
    RampUp_Constraints[init_t].set_function( LFunction );

   } else if( f_InitUpDownTime > 0 && f_InitUpDownTime < f_MinUpTime ) {

    auto LFunction = new LinearFunction();
    LFunction->add_variable( &v_active_power[init_t], 1.0 );
    LFunction->add_variable( &v_active_power[init_t - 1], -1.0 );
    LFunction->add_variable( &v_start_up[0], -min_power[init_t] );

    RampUp_Constraints[init_t].set_lhs( -Inf< double >() );
    RampUp_Constraints[init_t].set_rhs( delta_ramp_up[init_t] );
    RampUp_Constraints[init_t].set_function( LFunction );
   }

   for( Index t = init_t + 1, constraint_index = init_t + 1; t < f_time_horizon; ++t, ++constraint_index ) {

    auto lf = new LinearFunction();

    lf->add_variable( &v_active_power[t], -1.0 );
    lf->add_variable( &v_active_power[t - 1], 1.0 );
    lf->add_variable( &v_start_up[t - init_t], min_power[t] );

    lf->add_variable
            ( &v_commitment[t - 1], ( delta_ramp_up[t] ));

    RampUp_Constraints[constraint_index].set_lhs( 0.0 );
    RampUp_Constraints[constraint_index].set_rhs( Inf< double >());
    RampUp_Constraints[constraint_index].set_function( lf );
   }

  }
  add_static_constraint( RampUp_Constraints, "rampup_c" );


  // Initializing ramp down constraints


  RampDown_Constraints.resize( f_time_horizon );

  // Initial condition
  if( init_t == 0 ) {
   // Remaining constraints
   auto linear_function = new LinearFunction();
   linear_function->add_variable( &v_active_power[0], 1.0 );
   linear_function->add_variable
           ( &v_commitment[0], delta_ramp_down[0] );
   linear_function->add_variable
           ( &v_shut_down[0], min_power[0] );

   RampDown_Constraints[0].set_lhs
           ( f_initial_power );
   RampDown_Constraints[0].set_rhs( Inf< double >());
   RampDown_Constraints[0].set_function( linear_function );
   for( Index t = 1, constraint_index = 1; t < f_time_horizon; ++t, ++constraint_index ) {

    auto lf = new LinearFunction();

    lf->add_variable( &v_active_power[t], 1.0 );
    lf->add_variable( &v_active_power[t - 1], -1.0 );
    lf->add_variable
            ( &v_shut_down[t], min_power[t] );
    lf->add_variable
            ( &v_commitment[t], ( delta_ramp_down[t] ));

    RampDown_Constraints[constraint_index].set_lhs( 0.0 );
    RampDown_Constraints[constraint_index].set_rhs( Inf< double >());
    RampDown_Constraints[constraint_index].set_function( lf );
   }
  }

  if( init_t > 0 ) {
   // Initial condition
   auto linear_function = new LinearFunction();
   linear_function->add_variable( &v_active_power[0], 1.0 );
   linear_function->add_variable
           ( &v_commitment[0], delta_ramp_down[0] );

   RampDown_Constraints[0].set_lhs
           ( f_initial_power );
   RampDown_Constraints[0].set_rhs( Inf< double >());
   RampDown_Constraints[0].set_function( linear_function );

   for( Index t = 1, constraint_index = 1; t < init_t;
        ++t, ++constraint_index ) {
    auto lf = new LinearFunction();

    lf->add_variable( &v_active_power[t - 1], -1.0 );
    lf->add_variable( &v_active_power[t], 1.0 );
    lf->add_variable
            ( &v_commitment[t], delta_ramp_down[t] );
    RampDown_Constraints[constraint_index].set_lhs( 0.0 );
    RampDown_Constraints[constraint_index].set_rhs( Inf< double >());
    RampDown_Constraints[constraint_index].set_function( lf );
   }
   // Remaining constraints

   for( Index t = init_t, constraint_index = init_t; t < f_time_horizon; ++t, ++constraint_index ) {

    auto lf = new LinearFunction();

    lf->add_variable( &v_active_power[t], 1.0 );
    lf->add_variable( &v_active_power[t - 1], -1.0 );
    lf->add_variable
            ( &v_shut_down[t - init_t], min_power[t] );
    lf->add_variable
            ( &v_commitment[t], ( delta_ramp_down[t] ));

    RampDown_Constraints[constraint_index].set_lhs( 0.0 );
    RampDown_Constraints[constraint_index].set_rhs( Inf< double >());
    RampDown_Constraints[constraint_index].set_function( lf );
   }
  }

  add_static_constraint( RampDown_Constraints, "rampdown_c" );

 }
/*--------------------------------------------------------------------------*/

 // Initializing minimum power constraints
 if( !v_PrimaryRho.empty() & !v_SecondaryRho.empty()) {
  MinPower_Constraints.resize( f_time_horizon );

  for( Index t = 0; t < f_time_horizon; ++t ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_active_power[t], 1.0 );
   linear_function->add_variable( &v_primary_spinning_reserve[t], -1.0 );
   linear_function->add_variable( &v_secondary_spinning_reserve[t], -1.0 );
   linear_function->add_variable( &v_commitment[t], -min_power[t] );

   MinPower_Constraints[t].set_rhs( Inf< double >());
   MinPower_Constraints[t].set_lhs( 0.0 );
   MinPower_Constraints[t].set_function( linear_function );
  }
 } else {

  MinPower_Constraints.resize( f_time_horizon );

  for( Index t = 0; t < f_time_horizon; ++t ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_active_power[t], 1.0 );
   linear_function->add_variable( &v_commitment[t], -min_power[t] );

   MinPower_Constraints[t].set_rhs( Inf< double >());
   MinPower_Constraints[t].set_lhs( 0.0 );
   MinPower_Constraints[t].set_function( linear_function );
  }
 }
 add_static_constraint( MinPower_Constraints, "minpower_c" );

 // Initializing maximum power constraints
 if( !v_PrimaryRho.empty() & !v_SecondaryRho.empty()) {
  MaxPower_Constraints.resize( f_time_horizon );

  for( Index t = 0; t < f_time_horizon; ++t ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_active_power[t], -1.0 );
   linear_function->add_variable( &v_primary_spinning_reserve[ t ], -1.0 );
   linear_function->add_variable( &v_secondary_spinning_reserve[ t ], -1.0 );
   linear_function->add_variable( &v_commitment[t], max_power[t] );

   MaxPower_Constraints[t].set_lhs( 0.0 );
   MaxPower_Constraints[t].set_rhs( Inf< double >());
   MaxPower_Constraints[t].set_function( linear_function );
  }
 } else {

  MaxPower_Constraints.resize( f_time_horizon );

  for( Index t = 0; t < f_time_horizon; ++t ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_active_power[t], -1.0 );
   linear_function->add_variable( &v_commitment[t], max_power[t] );

   MaxPower_Constraints[t].set_lhs( 0.0 );
   MaxPower_Constraints[t].set_rhs( Inf< double >());
   MaxPower_Constraints[t].set_function( linear_function );
  }
 }
 add_static_constraint( MaxPower_Constraints, "maxpower_c" );


 // Initializing primary rho fraction constraints

 if( !v_PrimaryRho.empty() ) {
  PrimaryRho_Constraints.resize( f_time_horizon );

  for( Index t = 0; t < f_time_horizon; ++t ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_active_power[ t ], primary_rho[ t ] );
   linear_function->add_variable( &v_primary_spinning_reserve[ t ], -1.0 );

   PrimaryRho_Constraints[ t ].set_lhs( 0.0 );
   PrimaryRho_Constraints[ t ].set_rhs( Inf< double >() );
   PrimaryRho_Constraints[ t ].set_function( linear_function );
  }

  add_static_constraint( PrimaryRho_Constraints, "primaryrho_c" );
 }
 // Initializing secondary rho fraction constraints

 if( !v_SecondaryRho.empty() ) {
  SecondaryRho_Constraints.resize( f_time_horizon );

  for( Index t = 0; t < f_time_horizon; ++t ) {

   auto linear_function = new LinearFunction();
   linear_function->add_variable( &v_active_power[ t ],
                                  secondary_rho[ t ] );
   linear_function->add_variable( &v_secondary_spinning_reserve[ t ],
                                  -1.0 );

   SecondaryRho_Constraints[ t ].set_lhs( 0.0 );
   SecondaryRho_Constraints[ t ].set_rhs( Inf< double >() );
   SecondaryRho_Constraints[ t ].set_function( linear_function );
  }

  add_static_constraint( SecondaryRho_Constraints, "secondaryrho_c" );
 }

} // end( ThermalUnitBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::generate_objective( Configuration * objc ) {

// Initial condition of each vector
 std::vector<double> start_up_cost = this->v_StartUpCost;
 if( start_up_cost.size() == 1 ) {
  start_up_cost.resize( f_time_horizon - init_t, start_up_cost[ 0 ] );
 }else if (start_up_cost.size() < f_time_horizon - init_t) {
  start_up_cost.resize(f_time_horizon - init_t);
  int j = 0;
  for (unsigned long i = 0; i < v_change_intervals.size(); ++i) {
   Index sup;
   if (i == v_change_intervals.size() - 1 ) {
    sup = f_time_horizon;
   } else {
    sup = v_change_intervals[i];
   }
   for (; j < sup; ++j) {
    start_up_cost[j] = v_StartUpCost[i];
   }
  }
 }


 std::vector<double> linear_term = this->v_LinearTerm;
 if( linear_term.size() == 1 ) {
  linear_term.resize( f_time_horizon, linear_term[ 0 ] );
 }else if (linear_term.size() < f_time_horizon) {
  linear_term.resize(f_time_horizon);
  int j = 0;
  for (unsigned long i = 0; i < v_change_intervals.size(); ++i) {
   Index sup;
   if (i == v_change_intervals.size() - 1 ) {
    sup = f_time_horizon;
   } else {
    sup = v_change_intervals[i];
   }
   for (; j < sup; ++j) {
    linear_term[j] = v_LinearTerm[i];
   }
  }
 }

 std::vector<double> quad_term = this->v_QuadTerm;
 if( quad_term.size() == 1 ) {
  quad_term.resize( f_time_horizon, quad_term[ 0 ] );
 }else if (quad_term.size() < f_time_horizon) {
  quad_term.resize(f_time_horizon);
  int j = 0;
  for (unsigned long i = 0; i < v_change_intervals.size(); ++i) {
   Index sup;
   if (i == v_change_intervals.size() - 1 ) {
    sup = f_time_horizon;
   } else {
    sup = v_change_intervals[i];
   }
   for (; j < sup; ++j) {
    quad_term[j] = v_QuadTerm[i];
   }
  }
 }

 std::vector<double> const_term = this->v_ConstTerm;
 if( const_term.size() == 1 ) {
  const_term.resize( f_time_horizon, const_term[ 0 ] );
 }else if (const_term.size() < f_time_horizon) {
  const_term.resize(f_time_horizon);
  int j = 0;
  for (unsigned long i = 0; i < v_change_intervals.size(); ++i) {
   Index sup;
   if (i == v_change_intervals.size() - 1 ) {
    sup = f_time_horizon;
   } else {
    sup = v_change_intervals[i];
   }
   for (; j < sup; ++j) {
    const_term[j] = v_ConstTerm[i];
   }
  }
 }
/*--------------------------------------------------------------------------*/
 if( get_objective() != nullptr )  // an objective is there already
  return;                        // cowardly (and silently) return

 // Initialize objective function

 if( v_commitment.size() != f_time_horizon ) {
  throw ( std::logic_error
          ( "ThermalUnitBlock::generate_objective: v_commitment must have "
            "size equal to the time horizon." ) );
 }
 if( v_active_power.size() != f_time_horizon ) {
  throw ( std::logic_error
          ( "ThermalUnitBlock::generate_objective: v_active_power must have "
            "size equal to the time horizon." ) );
 }

 if( v_start_up.size() != f_time_horizon - init_t ) {
  throw ( std::logic_error
          ( "ThermalUnitBlock::generate_objective: v_start_up must have "
            "size equal to the time horizon - init_t." ) );
 }


 auto dquad_function = new DQuadFunction();

 for( Index t = init_t; t < f_time_horizon; ++t ) {
  dquad_function->add_variable( &v_start_up[ t - init_t ],
                                start_up_cost[ t - init_t ],
                                0.0 );
 }

 for( Index t = 0; t < f_time_horizon; ++t ) {
  dquad_function->add_variable( &v_active_power[ t ],
                                linear_term[ t ],
                                quad_term[ t ] );
  dquad_function->add_variable( &v_commitment[ t ],
                                const_term[ t ],
                                0.0 );
 }

 objective.set_function( dquad_function );
 objective.set_sense( Objective::eMin );

 // Set Block objective
 this->set_objective( &objective );

}  // end( ThermalUnitBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*-------- METHODS FOR LOADING, PRINTING & SAVING THE ThermalUnitBlock -----*/
/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::serialize( netCDF::NcGroup & group ) const {

 UnitBlock::serialize( group );

 ::serialize( group, "InitialPower", netCDF::NcDouble(), f_initial_power );
 ::serialize( group, "MinUpTime", netCDF::NcUint64(), f_MinUpTime );
 ::serialize( group, "MinDownTime", netCDF::NcUint64(), f_MinDownTime );
 ::serialize( group, "InitUpDownTime", netCDF::NcInt64(), f_InitUpDownTime );

 auto NumberIntervals = group.getDim( "NumberIntervals" );

 ::serialize( group, "MinPower", netCDF::NcDouble(),
              NumberIntervals, v_MinPower, true );

 ::serialize( group, "MaxPower", netCDF::NcDouble(),
              NumberIntervals, v_MaxPower, true );

 ::serialize( group, "DeltaRampUp", netCDF::NcDouble(),
              NumberIntervals, v_DeltaRampUp, true );

 ::serialize( group, "DeltaRampDown", netCDF::NcDouble(),
              NumberIntervals, v_DeltaRampDown, true );

 if( !v_PrimaryRho.empty() ) {
  ::serialize( group, "PrimaryRho", netCDF::NcDouble(),
               NumberIntervals, v_PrimaryRho, true );
 }

 if( !v_SecondaryRho.empty() ) {
  ::serialize( group, "SecondaryRho", netCDF::NcDouble(),
               NumberIntervals, v_SecondaryRho, true );
 }

 ::serialize( group, "QuadTerm", netCDF::NcDouble(),
              NumberIntervals, v_QuadTerm, true );

 ::serialize( group, "LinearTerm", netCDF::NcDouble(),
              NumberIntervals, v_LinearTerm, true );

 ::serialize( group, "ConstTerm", netCDF::NcDouble(),
              NumberIntervals, v_ConstTerm, true );
 ::serialize( group, "StartUpCost", netCDF::NcDouble(),
              NumberIntervals, v_StartUpCost, true );

 ::serialize( group, "FixedConsumption", netCDF::NcDouble(),
              {NumberIntervals}, v_fixed_consumption, true );

 ::serialize( group, "InertiaCommitment", netCDF::NcDouble(),
              {NumberIntervals}, v_inertia_commitment, true );
}  // end( ThermalUnitBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------------- End File ThermalUnitBlock.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
