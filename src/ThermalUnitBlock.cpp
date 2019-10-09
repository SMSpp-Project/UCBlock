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
#include "FRowConstraint.h"
#include "LinearFunction.h"
#include "ThermalUnitBlock.h"
#include "UCBlock.h"
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


 ::deserialize( group, "InitialPower", &f_initial_power );
 ::deserialize( group, "MinUpTime", &f_MinUpTime );
 ::deserialize( group, "MinDownTime", &f_MinDownTime );
 ::deserialize( group, "InitUpDownTime", &f_InitUpDownTime );

}  // end( ThermalUnitBlock::deserialize )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::generate_abstract_variables( Configuration * stvv ) {

 UnitBlock::generate_abstract_variables( stvv );

 // auto var = get_static_variable<boost::multi_array<ColVariable, 2>*>(0);


 if( f_InitUpDownTime > 0 ) {
  init_t = ( f_InitUpDownTime >= f_MinUpTime ? 0 :
             f_MinUpTime - f_InitUpDownTime );
 } else {
  init_t = ( -f_InitUpDownTime >= f_MinDownTime ? 0 :
             f_MinDownTime + f_InitUpDownTime );
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
    add_static_variable( i, "Startup " + std::to_string( n++ ) );
   }
  }

  if( v_shut_down.size() != startup_shutdown_size ) {
   assert( v_shut_down.empty() ); // this should only happen once
   v_shut_down.resize( startup_shutdown_size );
   int n = 0;
   for( auto & i : v_shut_down ) {
    i.set_type( ColVariable::kBinary );
    add_static_variable( i, "Shutdown " + std::to_string( n++ ) );
   }
  }
 }

/*--------------------------------------------------------------------------*/
 // POSSIBLY FIXING THE COMMITMENT VARIABLES TO 0 OR 1

 double commitment_variable_value = -1.0;

 if( f_InitUpDownTime < 0 && -f_InitUpDownTime < f_MinDownTime ) {

  commitment_variable_value = 0.0;

  // This unit must remain off from time 0 to init_t - 1. Therefore,
  // it should not produce any power.

  for( Index t = 0; t < init_t; ++t ) {
   v_active_power[ t ][ 0 ].set_value( 0.0 );
   v_active_power[ t ][ 0 ].is_fixed( true );
  }

  if( !v_primary_spinning_reserve.empty() ) {
   for( Index t = 0; t < init_t; ++t ) {
    v_primary_spinning_reserve[ t ][ 0 ].set_value( 0.0 );
    v_primary_spinning_reserve[ t ][ 0 ].is_fixed( true );
   }
  }

  if( !v_secondary_spinning_reserve.empty() ) {
   for( Index t = 0; t < init_t; ++t ) {
    v_secondary_spinning_reserve[ t ][ 0 ].set_value( 0.0 );
    v_secondary_spinning_reserve[ t ][ 0 ].is_fixed( true );
   }
  }
 } else if( f_InitUpDownTime > 0 && f_InitUpDownTime < f_MinUpTime ) {
  // This unit must remain on from time 0 to init_t - 1.
  commitment_variable_value = 1.0;
 }

 if( !v_commitment.empty() ) {
  if( commitment_variable_value >= 0.0 ) {
   // Fixing the commitment variables v_commitment to 0 or 1 from
   // time 0 to init_t - 1, depending on whether this unit must
   // remain off or on for the first init_t time steps.
   for( Index t = 0; t < init_t; ++t ) {
    v_commitment[ t ][ 0 ].set_value( commitment_variable_value );
    v_commitment[ t ][ 0 ].is_fixed( true );
   }
  }
 }

} // end( ThermalUnitBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::generate_abstract_constraints( Configuration * stcc ) {

 // FIXME: This should be implemented better
 std::vector v_MinPower = this->v_MinPower;
 if( v_MinPower.size() == 1 ) {
  v_MinPower.resize( f_time_horizon, v_MinPower[ 0 ] );
 }
 std::vector v_MaxPower = this->v_MaxPower;
 if( v_MaxPower.size() == 1 ) {
  v_MaxPower.resize( f_time_horizon, v_MaxPower[ 0 ] );
 }
 std::vector v_DeltaRampUp = this->v_DeltaRampUp;
 if( v_DeltaRampUp.size() == 1 ) {
  v_DeltaRampUp.resize( f_time_horizon, v_DeltaRampUp[ 0 ] );
 }
 std::vector v_DeltaRampDown = this->v_DeltaRampDown;
 if( v_DeltaRampDown.size() == 1 ) {
  v_DeltaRampDown.resize( f_time_horizon, v_DeltaRampDown[ 0 ] );
 }

 // MINIMUM UP AND DOWN TIME CONSTRAINTS

 // Initializing start up and shut down variables connection constraints

 if( f_time_horizon - init_t > 0 ) {

  StartUp_ShutDown_Variables_Constraints.resize( f_time_horizon - init_t );

  for( Index t = init_t, constraint_index = 0; t < f_time_horizon;
       ++t, ++constraint_index ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_commitment[ t ][ 0 ], 1.0 );
   linear_function->add_variable( &v_start_up[ t - init_t ], -1.0 );
   linear_function->add_variable( &v_shut_down[ t - init_t ], 1.0 );

   if( t > 0 ) {
    linear_function->add_variable( &v_commitment[ t - 1 ][ 0 ], -1.0 );
    StartUp_ShutDown_Variables_Constraints[ constraint_index ].set_both( 0.0 );
   } else {
    StartUp_ShutDown_Variables_Constraints[ constraint_index ].
     set_both( f_InitUpDownTime > 0 ? 1.0 : 0.0 );
   }

   StartUp_ShutDown_Variables_Constraints[ constraint_index ].
    set_function( linear_function );
  }
  add_static_constraint( StartUp_ShutDown_Variables_Constraints,
                         "Startup/Shutdown Variables Constraints" );
 }

 // Initializing turn on constraints (start up constraints)

 if( f_time_horizon - init_t - f_MinUpTime > 0 ) {

  StartUp_Constraints.resize( f_time_horizon - init_t - f_MinUpTime );

  for( Index t = init_t + f_MinUpTime, constraint_index = 0;
       t < f_time_horizon; ++t, ++constraint_index ) {

   auto linear_function = new LinearFunction();

   for( Index s = t - f_MinUpTime; s <= t; ++s ) {
    linear_function->add_variable( &v_start_up[ s - init_t ], 1.0 );
   }

   linear_function->add_variable( &v_commitment[ t ][ 0 ], -1.0 );
   StartUp_Constraints[ constraint_index ].set_lhs( -Inf< double >() );
   StartUp_Constraints[ constraint_index ].set_rhs( 0.0 );
   StartUp_Constraints[ constraint_index ].set_function( linear_function );
  }

  add_static_constraint( StartUp_Constraints, "Startup Constraints" );
 }

 // Initializing turn off constraints (shut down constraints)
 if( f_time_horizon - init_t - f_MinDownTime > 0 ) {

  ShutDown_Constraints.resize( f_time_horizon - init_t - f_MinDownTime );

  for( Index t = init_t + f_MinDownTime, constraint_index = 0;
       t < f_time_horizon; ++t, ++constraint_index ) {

   auto linear_function = new LinearFunction();

   for( Index s = t - f_MinDownTime; s <= t; ++s ) {
    linear_function->add_variable( &v_shut_down[ s - init_t ], 1.0 );
   }

   linear_function->add_variable( &v_commitment[ t ][ 0 ], 1.0 );
   ShutDown_Constraints[ constraint_index ].set_lhs( -Inf< double >() );
   ShutDown_Constraints[ constraint_index ].set_rhs( 1.0 );
   ShutDown_Constraints[ constraint_index ].set_function( linear_function );
  }

  add_static_constraint( ShutDown_Constraints, "Shutdown constraints" );
 }

/*--------------------------------------------------------------------------*/
 // RAMP UP AND RAMP DOWN CONSTRAINTS

 // Initializing ramp-up constraints

 {
  RampUp_Constraints.resize( f_time_horizon );

  // Initial condition
  for( Index t = 0; t < init_t; ++t ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_active_power[ t + 1 ][ 0 ], -1.0 );
   linear_function->add_variable( &v_active_power[ t ][ 0 ], 1.0 );

   RampUp_Constraints[ t ].set_lhs( -Inf< double >() );
   RampUp_Constraints[ t ].set_rhs( -v_DeltaRampUp[ t ] );
   RampUp_Constraints[ t ].set_function( linear_function );
  }
  // Remaining constraints

  for( Index t = init_t; t < f_time_horizon - 1; ++t ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_active_power[ t + 1 ][ 0 ], -1.0 );
   linear_function->add_variable( &v_active_power[ t ][ 0 ], 1.0 );
   linear_function->add_variable( &v_start_up[ t - init_t ], -v_DeltaRampUp[ t - init_t ] );

   linear_function->add_variable
    ( &v_commitment[ t + 1 ][ 0 ], ( v_MinPower[ t ] + v_DeltaRampUp[ t ] ) );

   linear_function->add_variable( &v_commitment[ t ][ 0 ], -v_MinPower[ t ] );

   RampUp_Constraints[ t ].set_lhs( 0.0 );
   RampUp_Constraints[ t ].set_rhs( Inf< double >() );
   RampUp_Constraints[ t ].set_function( linear_function );
  }

  // THIS IS NOT GOOD -----
  auto linear_function = new LinearFunction();
  RampUp_Constraints[ f_time_horizon - 1 ].set_lhs( -Inf< double >() );
  RampUp_Constraints[ f_time_horizon - 1 ].set_rhs( Inf< double >() );
  RampUp_Constraints[ f_time_horizon - 1 ].set_function( linear_function );
  // --------------------

  add_static_constraint( RampUp_Constraints, "Ramp Up Constraints" );
 }

 // Initializing ramp down constraints

 {
  RampDown_Constraints.resize( f_time_horizon );

  // Initial condition
  for( Index t = 0; t < init_t; ++t ) {
   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_active_power[ t ][ 0 ], 1.0 );
   linear_function->add_variable( &v_active_power[ t + 1 ][ 0 ], -1.0 );

   RampDown_Constraints[ t ].set_lhs( 0.0 );
   RampDown_Constraints[ t ].set_rhs( v_DeltaRampDown[ t ] );
   RampDown_Constraints[ t ].set_function( linear_function );
  }
  // Remaining constraints

  for( Index t = init_t; t < f_time_horizon - 1; ++t ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_active_power[ t + 1 ][ 0 ], 1.0 );
   linear_function->add_variable( &v_active_power[ t ][ 0 ], -1.0 );
   linear_function->add_variable
    ( &v_shut_down[ t - init_t ], -v_DeltaRampDown[ t - init_t ] );
   linear_function->add_variable
    ( &v_commitment[ t ][ 0 ], ( v_MinPower[ t ] + v_DeltaRampDown[ t ] ) );
   linear_function->add_variable( &v_commitment[ t + 1 ][ 0 ], -v_MinPower[ t ] );

   RampDown_Constraints[ t ].set_lhs( 0.0 );
   RampDown_Constraints[ t ].set_rhs( Inf< double >() );
   RampDown_Constraints[ t ].set_function( linear_function );
  }

  // THIS IS NOT GOOD -----
  auto linear_function = new LinearFunction();
  RampDown_Constraints[ f_time_horizon - 1 ].set_lhs( -Inf< double >() );
  RampDown_Constraints[ f_time_horizon - 1 ].set_rhs( Inf< double >() );
  RampDown_Constraints[ f_time_horizon - 1 ].set_function( linear_function );
  // --------------------

  add_static_constraint( RampDown_Constraints, "Ramp Down Constraints" );
 }

/*--------------------------------------------------------------------------*/

 // POWER OUTPUT CONSTRAINTS

 // Initial condition for power output startup and shutdown constraints

 if( f_MinUpTime >= 2 ) {

  if( f_time_horizon - init_t > 0 ) {

   Power_StartUp_ShutDown_Variables_Constraints.resize
    ( f_time_horizon - init_t - 2 );

   // Initializing power output startup and shutdown constraints

   for( Index t = init_t + 1, constraint_index = 0; t < f_time_horizon - 1;
        ++t, ++constraint_index ) {

    auto linear_function = new LinearFunction();

    linear_function->add_variable( &v_active_power[ t - 1 ][ 0 ], -1.0 );
    linear_function->add_variable( &v_commitment[ t - 1 ][ 0 ], v_MaxPower[ t - 1 ] );
    linear_function->add_variable( &v_start_up[ t - init_t - 1 ],
                                   -( v_MaxPower[ t - 1 ] - v_MinPower[ t - 1 ] ) );
    linear_function->add_variable( &v_shut_down[ t - init_t ],
                                   ( v_MaxPower[ t - 1 ] - v_MinPower[ t - 1 ] ) );

    Power_StartUp_ShutDown_Variables_Constraints[ constraint_index ].
     set_lhs( 0.0 );
    Power_StartUp_ShutDown_Variables_Constraints[ constraint_index ].
     set_rhs( Inf< double >() );

    Power_StartUp_ShutDown_Variables_Constraints[ constraint_index ].
     set_function( linear_function );
   }

   add_static_constraint( Power_StartUp_ShutDown_Variables_Constraints,
                          "Power Output Startup and Shutdown Constraints" );
  }
 } // end power output startup and shutdown constraints


 if( f_MinUpTime == 1 ) {

  // Initializing power output startup constraints

  if( f_time_horizon - init_t > 0 ) {

   Power_StartUp_Variable_Constraints.resize
    ( f_time_horizon - init_t - 2 );

   for( Index t = init_t + 2, constraint_index = 0; t < f_time_horizon;
        ++t, ++constraint_index ) {

    auto linear_function = new LinearFunction();

    linear_function->add_variable( &v_active_power[ t ][ 0 ], -1.0 );
    linear_function->add_variable( &v_commitment[ t ][ 0 ],
                                   v_MaxPower[ t ] );
    linear_function->add_variable( &v_start_up[ t ],
                                   -( v_MaxPower[ t ]
                                      - v_MinPower[ t ] ) );

    Power_StartUp_Variable_Constraints[ constraint_index ].
     set_lhs( 0.0 );
    Power_StartUp_Variable_Constraints[ constraint_index ].
     set_rhs( Inf< double >() );

    Power_StartUp_Variable_Constraints[ constraint_index ].
     set_function( linear_function );
   }

   add_static_constraint( Power_StartUp_Variable_Constraints,
                          "Power Output Startup Constraints" );

  } // end power output startup  constraints


  // power output shut down constraints

  if( f_time_horizon - init_t > 0 ) {

   Power_ShutDown_Variable_Constraints.resize
    ( f_time_horizon - init_t - 2 );


   // Initializing power output shutdown constraints

   for( Index t = init_t + 2, constraint_index = 0; t < f_time_horizon;
        ++t, ++constraint_index ) {

    auto linear_function = new LinearFunction();

    linear_function->add_variable( &v_active_power[ t ][ 0 ], -1.0 );
    linear_function->add_variable( &v_commitment[ t ][ 0 ],
                                   v_MaxPower[ t ] );
    linear_function->add_variable( &v_shut_down[ t + 1 ],
                                   -( v_MaxPower[ t + 1 ]
                                      - v_MinPower[ t + 1 ] ) );

    Power_ShutDown_Variable_Constraints[ constraint_index ].
     set_lhs( 0.0 );
    Power_ShutDown_Variable_Constraints[ constraint_index ].
     set_rhs( Inf< double >() );

    Power_ShutDown_Variable_Constraints[ constraint_index ].
     set_function( linear_function );
   }

   add_static_constraint( Power_ShutDown_Variable_Constraints,
                          "Power Output Shutdown Constraints" );

  } // end power output shut down constraints

 }

/*--------------------------------------------------------------------------*/

 // Initializing minimum power constraints

 MinPower_Constraints.resize( f_time_horizon );

 for( Index t = 0; t < f_time_horizon; ++t ) {

  auto linear_function = new LinearFunction();

  linear_function->add_variable( &v_active_power[ t ][ 0 ], 1.0 );
  linear_function->add_variable( &v_primary_spinning_reserve[ t ][ 0 ], -1.0 );
  linear_function->add_variable( &v_secondary_spinning_reserve[ t ][ 0 ], -1.0 );
  linear_function->add_variable( &v_commitment[ t ][ 0 ], -v_MinPower[ t ] );

  MinPower_Constraints[ t ].set_rhs( Inf< double >() );
  MinPower_Constraints[ t ].set_lhs( 0.0 );
  MinPower_Constraints[ t ].set_function( linear_function );
 }

 add_static_constraint( MinPower_Constraints, "Min Power Constraints" );

 // Initializing maximum power constraints

 MaxPower_Constraints.resize( f_time_horizon );

 for( Index t = 0; t < f_time_horizon; ++t ) {

  auto linear_function = new LinearFunction();

  linear_function->add_variable( &v_active_power[ t ][ 0 ], 1.0 );
  linear_function->add_variable( &v_primary_spinning_reserve[ t ][ 0 ], 1.0 );
  linear_function->add_variable( &v_secondary_spinning_reserve[ t ][ 0 ], 1.0 );
  linear_function->add_variable( &v_commitment[ t ][ 0 ], -v_MaxPower[ t ] );

  MaxPower_Constraints[ t ].set_lhs( -Inf< double >() );
  MaxPower_Constraints[ t ].set_rhs( 0.0 );
  MaxPower_Constraints[ t ].set_function( linear_function );
 }

 add_static_constraint( MaxPower_Constraints, "Max Power Constraints" );

 // Initializing primary rho fraction constraints

 if( !v_PrimaryRho.empty() ) {
  PrimaryRho_Constraints.resize( f_time_horizon );

  for( Index t = 0; t < f_time_horizon; ++t ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_active_power[ t ][ 0 ], v_PrimaryRho[ t ] );
   linear_function->add_variable( &v_primary_spinning_reserve[ t ][ 0 ], -1.0 );

   PrimaryRho_Constraints[ t ].set_lhs( 0.0 );
   PrimaryRho_Constraints[ t ].set_rhs( Inf< double >() );
   PrimaryRho_Constraints[ t ].set_function( linear_function );
  }

  add_static_constraint( PrimaryRho_Constraints, "Primary Rho Constraints" );
 }
 // Initializing secondary rho fraction constraints

 if( !v_SecondaryRho.empty() ) {
  SecondaryRho_Constraints.resize( f_time_horizon );

  for( Index t = 0; t < f_time_horizon; ++t ) {

   auto linear_function = new LinearFunction();
   linear_function->add_variable( &v_active_power[ t ][ 0 ],
                                  v_SecondaryRho[ t ] );
   linear_function->add_variable( &v_secondary_spinning_reserve[ t ][ 0 ],
                                  -1.0 );

   SecondaryRho_Constraints[ t ].set_lhs( 0.0 );
   SecondaryRho_Constraints[ t ].set_rhs( Inf< double >() );
   SecondaryRho_Constraints[ t ].set_function( linear_function );
  }

  add_static_constraint( SecondaryRho_Constraints, "Secondary Rho Constraints" );
 }

/*--------------------------------------------------------------------------*/
 // TIME DEPENDENT START UP COSTS CONSTRAINTS

 //TODO if any exist

} // end( ThermalUnitBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::generate_objective( Configuration * objc ) {

 if( get_objective() != nullptr )  // an objective is there already
  return;                        // cowardly (and silently) return

 // Initialize objective function

 if( v_commitment.size() != f_time_horizon ) {
  throw ( std::logic_error
   ( "ThermalUnitBlock::generate_objective: v_commitment must have "
     "size equal to the time horizon." ) );
 }

 // FIXME: This should be implemented better
 std::vector v_StartUpCost = this->v_StartUpCost;
 if( v_StartUpCost.size() == 1 ) {
  v_StartUpCost.resize( f_time_horizon - init_t, v_StartUpCost[ 0 ] );
 }
 std::vector v_LinearTerm = this->v_LinearTerm;
 if( v_LinearTerm.size() == 1 ) {
  v_LinearTerm.resize( f_time_horizon, v_LinearTerm[ 0 ] );
 }
 std::vector v_QuadTerm = this->v_QuadTerm;
 if( v_QuadTerm.size() == 1 ) {
  v_QuadTerm.resize( f_time_horizon, v_QuadTerm[ 0 ] );
 }
 std::vector v_ConstTerm = this->v_ConstTerm;
 if( v_ConstTerm.size() == 1 ) {
  v_ConstTerm.resize( f_time_horizon, v_ConstTerm[ 0 ] );
 }

 auto dquad_function = new DQuadFunction();

 for( Index t = init_t; t < f_time_horizon; ++t ) {
  dquad_function->add_variable( &v_start_up[ t - init_t ],
                                v_StartUpCost[ t - init_t ],
                                0.0 );
 }

 for( Index t = 0; t < f_time_horizon; ++t ) {
  dquad_function->add_variable( &v_active_power[ t ][ 0 ],
                                v_LinearTerm[ t ],
                                v_QuadTerm[ t ] );
  dquad_function->add_variable( &v_commitment[ t ][ 0 ],
                                v_ConstTerm[ t ],
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

}  // end( ThermalUnitBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------------- End File ThermalUnitBlock.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
