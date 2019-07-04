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

//void ThermalUnitBlock::load( )
//{ }

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::deserialize( netCDF::NcGroup & group ) {

  UnitBlock::deserialize( group );

  ::deserialize( group, "MinPower",      f_number_intervals, v_MinPower );
  ::deserialize( group, "MaxPower",      f_number_intervals, v_MaxPower );
  ::deserialize( group, "DeltaRampUp",   f_number_intervals, v_DeltaRampUp );
  ::deserialize( group, "DeltaRampDown", f_number_intervals, v_DeltaRampDown );
  ::deserialize( group, "PrimaryRho",    f_number_intervals, v_PrimaryRho );
  ::deserialize( group, "SecondaryRho",  f_number_intervals, v_SecondaryRho );
  ::deserialize( group, "LinearTerm",    f_number_intervals, v_LinearTerm );
  ::deserialize( group, "QuadTerm",      f_number_intervals, v_QuadTerm );
  ::deserialize( group, "ConstTerm",     f_number_intervals, v_ConstTerm );
  ::deserialize( group, "StartUpCost",   f_number_intervals,  v_StartUpCost );

  ::deserialize( group, "InitialPower",         & f_initial_power );
  ::deserialize( group, "InitialMinPower",      & f_initial_min_power );
  ::deserialize( group, "InitialDeltaRampUp",   & f_initial_delta_ramp_up );
  ::deserialize( group, "InitialDeltaRampDown", & f_initial_delta_ramp_down );
  ::deserialize( group, "MinUpTime",            & f_MinUpTime );
  ::deserialize( group, "MinDownTime",          & f_MinDownTime );
  ::deserialize( group, "InitUpDownTime",       & f_InitUpDownTime );

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
  // START UP AND SHUT DOWN BINARY VARIABLES
auto startup_shutdown_size = f_time_horizon - init_t;

 if( startup_shutdown_size > 0 ) {

  if( v_start_up.size() != startup_shutdown_size &&
      v_shut_down.size() != startup_shutdown_size ) {

   assert( v_start_up.empty() ); // this should only happen once
   v_start_up.resize( startup_shutdown_size );

   assert( v_shut_down.empty() ); // this should only happen once
   v_shut_down.resize( startup_shutdown_size );

   for( Index i = 0; i < startup_shutdown_size; ++i ) {
    v_start_up[ i ].set_type( ColVariable::kBinary );
    v_shut_down[ i ].set_type( ColVariable::kBinary );
   }

   add_static_variable ( v_start_up );
   add_static_variable ( v_shut_down );
  }
 }

/*--------------------------------------------------------------------------*/
 // POSSIBLY FIXING THE COMMITMENT VARIABLES TO 0 OR 1

 double commitment_variable_value = - 1.0;

 if( f_InitUpDownTime < 0 && - f_InitUpDownTime < f_MinDownTime ) {

  commitment_variable_value = 0.0;

  // This unit must remain off from time 0 to init_t - 1. Therefore,
  // it should not produce any power.

  for( Index t = 0; t < init_t; ++t ) {
   v_active_power[ t ].set_value( 0.0 );
   v_active_power[ t ].is_fixed( true );
  }

  if( !v_primary_spinning_reserve.empty() ) {
   for( Index t = 0; t < init_t; ++t ) {
    v_primary_spinning_reserve[ t ].set_value( 0.0 );
    v_primary_spinning_reserve[ t ].is_fixed( true );
   }
  }

  if( !v_secondary_spinning_reserve.empty() ) {
   for( Index t = 0; t < init_t; ++t ) {
    v_secondary_spinning_reserve[ t ].set_value( 0.0 );
    v_secondary_spinning_reserve[ t ].is_fixed( true );
   }
  }
 }

 else if( f_InitUpDownTime > 0 && f_InitUpDownTime < f_MinUpTime ) {
  // This unit must remain on from time 0 to init_t - 1.
  commitment_variable_value = 1.0;
 }

 if( commitment_variable_value >= 0.0 ) {
  // Fixing the commitment variables v_commitment to 0 or 1 from
  // time 0 to init_t - 1, depending on whether this unit must
  // remain off or on for the first init_t time steps.
  for( Index t = 0; t < init_t ; ++t ) {
   v_commitment[ t ].set_value( commitment_variable_value );
   v_commitment[ t ].is_fixed( true );
  }
 }

} // end( ThermalUnitBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::generate_abstract_constraints( Configuration *stcc ) {
 // FIXME: multiple linear_function declarations shadow the local variable

 // MINIMUM UP AND DOWN TIME CONSTRAINTS

 // Initializing start up and shut down variables connection constraints

 if( f_time_horizon - init_t > 0 ) {

  StartUp_ShutDown_Variables_Constraints.resize( f_time_horizon - init_t );

  for( Index t = init_t, constraint_index = 0; t < f_time_horizon;
       ++t, ++constraint_index ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( & v_commitment[ t ],  1.0 );
   linear_function->add_variable( & start_up( t ),     -1.0 );
   linear_function->add_variable( & shut_down( t ),     1.0 );

   if( t > 0 ) [[likely]] {
    linear_function->add_variable( & v_commitment[ t - 1 ], -1.0 );
    StartUp_ShutDown_Variables_Constraints[ constraint_index ].
            set_both( 0.0 );
   }
   else {
    StartUp_ShutDown_Variables_Constraints[ constraint_index ].
            set_both( ( f_InitUpDownTime > 0 ? 1.0 : 0.0 ) );
   }

   StartUp_ShutDown_Variables_Constraints[ constraint_index ].
           set_function( linear_function );
  }

  add_static_constraint( StartUp_ShutDown_Variables_Constraints );
 }

 // Initializing turn on constraints (start up constraints)

 if( f_time_horizon - init_t - f_MinUpTime > 0 ) {

  StartUp_Constraints.resize( f_time_horizon - init_t - f_MinUpTime );

  for( Index t = init_t + f_MinUpTime, constraint_index = 0;
       t < f_time_horizon; ++t, ++constraint_index ) {

   auto linear_function = new LinearFunction();

   for( Index s = t - f_MinUpTime; s <= t; ++s ) {
    linear_function->add_variable( & start_up( s ), 1.0 );
   }

   linear_function->add_variable( & v_commitment[ t ], -1.0 );

   StartUp_Constraints[ constraint_index ].set_lhs( - 1.0 );
   StartUp_Constraints[ constraint_index ].set_rhs(   0.0 );
   StartUp_Constraints[ constraint_index ].set_function( linear_function );
  }

  add_static_constraint( StartUp_Constraints );
 }

 // Initializing turn off constraints (shut down constraints)
 if( f_time_horizon - init_t - f_MinDownTime > 0 ) {

  ShutDown_Constraints.resize( f_time_horizon - init_t - f_MinDownTime );

  for( Index t = init_t + f_MinDownTime, constraint_index = 0;
       t < f_time_horizon; ++t, ++constraint_index ) {

   auto linear_function = new LinearFunction();

   for( Index s = t - f_MinDownTime; s <= t; ++s ) {
    linear_function->add_variable( & shut_down( s ), 1.0 );
   }

   linear_function->add_variable( & v_commitment[ t ], 1.0 );
   ShutDown_Constraints[ constraint_index ].set_lhs( 0.0 );
   ShutDown_Constraints[ constraint_index ].set_rhs( 1.0 );
   ShutDown_Constraints[ constraint_index ].set_function( linear_function );
  }

  add_static_constraint( ShutDown_Constraints );
 }

/*--------------------------------------------------------------------------*/
 // RAMP UP AND RAMP DOWN CONSTRAINTS

 // Initializing ramp-up constraints

 {
  RampUp_Constraints.resize( f_time_horizon );

  // Initial condition

  auto linear_function = new LinearFunction();

  linear_function->add_variable( & v_active_power[ 0 ], -1.0 );
  linear_function->add_variable( & start_up( 0 ), - f_initial_delta_ramp_up );
  linear_function->add_variable
          ( &v_commitment[ 0 ], ( f_initial_min_power + f_initial_delta_ramp_up ) );

  auto initial_commitment = ( f_InitUpDownTime > 0 ? 1.0 : 0.0 );

  RampUp_Constraints[ 0 ].set_lhs
          ( f_initial_min_power * initial_commitment - f_initial_power );
  RampUp_Constraints[ 0 ].set_rhs( Inf<double>() );
  RampUp_Constraints[ 0 ].set_function( linear_function );

  // Remaining constraints

  for( Index t = 0, constraint_index = 1; t < f_time_horizon - 1;
       ++t, ++constraint_index ) {

   auto linear_function = new LinearFunction();
   RampUp_Constraints[ constraint_index ].set_function( linear_function );

   linear_function->add_variable( & v_active_power[ t + 1 ], -1.0 );
   linear_function->add_variable( & v_active_power[ t ],     +1.0 );
   linear_function->add_variable
           ( & start_up( t + 1 ), - v_DeltaRampUp[ t ] );

   linear_function->add_variable
           ( & v_commitment[ t + 1 ], ( v_MinPower[ t ] + v_DeltaRampUp[ t ] ) );

   linear_function->add_variable( & v_commitment[ t ], - v_MinPower[ t ] );

   RampUp_Constraints[ constraint_index ].set_lhs( 0.0 );
   RampUp_Constraints[ constraint_index ].set_rhs( Inf<double>() );
   RampUp_Constraints[ constraint_index ].set_function( linear_function );
  }

  add_static_constraint( RampUp_Constraints );
 }

 // Initializing ramp down constraints

 {
  RampDown_Constraints.resize( f_time_horizon );

  // Initial condition

  auto linear_function = new LinearFunction();

  linear_function->add_variable( & v_active_power[ 0 ], 1.0 );
  linear_function->add_variable
          ( & shut_down( 0 ), - f_initial_delta_ramp_down );
  linear_function->add_variable( & v_commitment[ 0 ], - f_initial_min_power );

  auto initial_commitment = ( f_InitUpDownTime > 0 ? 1.0 : 0.0 );

  RampDown_Constraints[ 0 ].set_lhs
          ( f_initial_power - initial_commitment *
                              ( f_initial_min_power + f_initial_delta_ramp_down ) );
  RampDown_Constraints[ 0 ].set_rhs( Inf<double>() );
  RampDown_Constraints[ 0 ].set_function( linear_function );

  // Remaining constraints

  for( Index t = 0, constraint_index = 1; t < f_time_horizon - 1;
       ++t, ++constraint_index ) {

   auto linear_function = new LinearFunction();
   RampDown_Constraints[ constraint_index ].set_function( linear_function );

   linear_function->add_variable( & v_active_power[ t + 1], 1.0 );
   linear_function->add_variable( & v_active_power[ t ],   -1.0 );
   linear_function->add_variable
           ( & shut_down( t + 1 ), - v_DeltaRampDown[ t ] );
   linear_function->add_variable
           ( & v_commitment[ t ], ( v_MinPower[ t ] + v_DeltaRampDown[ t ] ) );
   linear_function->add_variable( &v_commitment[ t + 1 ], - v_MinPower[ t ]);

   RampDown_Constraints[ constraint_index ].set_lhs( 0.0 );
   RampDown_Constraints[ constraint_index ].set_rhs( Inf<double>() );
   RampDown_Constraints[ constraint_index ].set_function( linear_function );
  }

  add_static_constraint( RampDown_Constraints );
 }

/*--------------------------------------------------------------------------*/

 // POWER OUTPUT CONSTRAINTS

 // Initial condition for power output startup and shutdown constraints

 if ( f_MinUpTime >= 2 ) {

  if( f_time_horizon - init_t > 0 ) {

   Power_StartUp_ShutDown_Variables_Constraints.resize
           ( f_time_horizon - init_t );

   auto linear_function = new LinearFunction();

   // Constraints at last time step

   linear_function->add_variable( & v_active_power[ f_time_horizon ], - 1.0);
   linear_function->add_variable( & v_commitment[ f_time_horizon ],
                                  v_MaxPower[ f_time_horizon ] );
   linear_function->add_variable( & start_up( f_time_horizon ),
                                  -( v_MaxPower[ f_time_horizon ]
                                     - v_MinPower[ f_time_horizon ] ) );

   Power_StartUp_ShutDown_Variables_Constraints[ f_time_horizon ].
           set_lhs( 0.0 );
   Power_StartUp_ShutDown_Variables_Constraints[ f_time_horizon ].
           set_rhs( Inf<double>());
   Power_StartUp_ShutDown_Variables_Constraints[ f_time_horizon ].
           set_function( linear_function );

   // Initializing power output startup and shutdown constraints

   for( Index t = init_t , constraint_index = 0; t < f_time_horizon - 1;
        ++t, ++constraint_index ) {

    auto linear_function = new LinearFunction();

    linear_function->add_variable( & v_active_power[ t ],  - 1.0 );
    linear_function->add_variable( & v_commitment[ t ],  v_MaxPower[ t ] );
    linear_function->add_variable( & start_up( t ),
                                   -( v_MaxPower[ t ] - v_MinPower[ t ] ) );
    linear_function->add_variable( & shut_down( t + 1 ),
                                   ( v_MaxPower[ t + 1 ] - v_MinPower[ t + 1] ) );

    Power_StartUp_ShutDown_Variables_Constraints[ constraint_index ].
            set_lhs( 0.0 );
    Power_StartUp_ShutDown_Variables_Constraints[ constraint_index ].
            set_rhs( Inf<double>());

    Power_StartUp_ShutDown_Variables_Constraints[ constraint_index ].
            set_function( linear_function );
   }

   add_static_constraint( Power_StartUp_ShutDown_Variables_Constraints );
  }
 } // end power output startup and shutdown constraints


 if ( f_MinUpTime == 1 ) {

  // Initializing power output startup constraints

  if( f_time_horizon - init_t > 0 ) {

   Power_StartUp_Variable_Constraints.resize
           ( f_time_horizon - init_t );

   for( Index t = init_t , constraint_index = 0; t < f_time_horizon ;
        ++t, ++constraint_index ) {

    auto linear_function = new LinearFunction();

    linear_function->add_variable( & v_active_power[ t ],  - 1.0 );
    linear_function->add_variable( & v_commitment[ t ],
                                   v_MaxPower[ t ] );
    linear_function->add_variable( & start_up( t ),
                                   -( v_MaxPower[ t ]
                                      - v_MinPower[ t ] ) );

    Power_StartUp_Variable_Constraints[ constraint_index ].
            set_lhs( 0.0 );
    Power_StartUp_Variable_Constraints[ constraint_index ].
            set_rhs( Inf<double>());

    Power_StartUp_Variable_Constraints[ constraint_index ].
            set_function( linear_function );
   }

   add_static_constraint( Power_StartUp_Variable_Constraints );

  } // end power output startup  constraints


  // power output shut down constraints

  if( f_time_horizon - init_t > 0 ) {

   Power_ShutDown_Variable_Constraints.resize
           ( f_time_horizon - init_t );

   // Initial condition

   auto linear_function = new LinearFunction();

   // Constraints at last time step

   linear_function->add_variable( & v_active_power[ f_time_horizon ],
                                  - 1.0);
   linear_function->add_variable( & v_commitment[ f_time_horizon ],
                                  v_MaxPower[ f_time_horizon ] );
   linear_function->add_variable( & start_up( f_time_horizon ),
                                  -( v_MaxPower[ f_time_horizon ]
                                     - v_MinPower[ f_time_horizon ] ) );

   Power_ShutDown_Variable_Constraints[ f_time_horizon ].
           set_lhs( 0.0 );
   Power_ShutDown_Variable_Constraints[ f_time_horizon ].
           set_rhs( Inf<double>());
   Power_ShutDown_Variable_Constraints[ f_time_horizon ].

           set_function( linear_function );

   // Initializing power output shutdown constraints

   for( Index t = init_t, constraint_index = 0; t < f_time_horizon - 1;
        ++t, ++constraint_index ) {

    auto linear_function = new LinearFunction();

    linear_function->add_variable( & v_active_power[ t ],  - 1.0 );
    linear_function->add_variable( & v_commitment[ t ],
                                   v_MaxPower[ t ] );
    linear_function->add_variable( & shut_down( t + 1 ),
                                   -( v_MaxPower[ t + 1 ]
                                      - v_MinPower[ t + 1 ] ) );

    Power_ShutDown_Variable_Constraints[ constraint_index ].
            set_lhs( 0.0 );
    Power_ShutDown_Variable_Constraints[ constraint_index ].
            set_rhs( Inf<double>());

    Power_ShutDown_Variable_Constraints[ constraint_index ].
            set_function( linear_function );
   }

   add_static_constraint( Power_ShutDown_Variable_Constraints );

  } // end power output shut down constraints

 }

/*--------------------------------------------------------------------------*/

 // Initializing minimum power constraints

 MinPower_Constraints.resize( f_time_horizon );

 for( Index t = 0; t < f_time_horizon; ++t ) {

  auto linear_function = new LinearFunction();

  linear_function->add_variable( & v_active_power[ t ],                1.0 );
  linear_function->add_variable( & v_primary_spinning_reserve[ t ] ,  -1.0 );
  linear_function->add_variable( & v_secondary_spinning_reserve[ t ], -1.0 );
  linear_function->add_variable( & v_commitment[ t ],    - v_MinPower[ t ] );

  MinPower_Constraints[ t ].set_rhs( Inf<double>() );
  MinPower_Constraints[ t ].set_lhs( 0.0 );
  MinPower_Constraints[ t ].set_function( linear_function );
 }

 add_static_constraint( MinPower_Constraints );

 // Initializing maximum power constraints

 MaxPower_Constraints.resize( f_time_horizon );

 for( Index t = 0; t < f_time_horizon; ++t ) {

  auto linear_function = new LinearFunction();

  linear_function->add_variable( & v_active_power[ t ],               1.0 );
  linear_function->add_variable( & v_primary_spinning_reserve[ t ],   1.0 );
  linear_function->add_variable( & v_secondary_spinning_reserve[ t ], 1.0 );
  linear_function->add_variable( & v_commitment[ t ],   - v_MaxPower[ t ] );

  MaxPower_Constraints[ t ].set_lhs( -Inf<double>() );
  MaxPower_Constraints[ t ].set_rhs( 0.0 );
  MaxPower_Constraints[ t ].set_function( linear_function );

 }

 add_static_constraint( MaxPower_Constraints );

 // Initializing primary rho fraction constraints

 PrimaryRho_Constraints.resize( f_time_horizon );

 for( Index t = 0; t < f_time_horizon; ++t ) {

  auto linear_function = new LinearFunction();

  linear_function->add_variable( & v_active_power[ t ], v_PrimaryRho[ t ] );
  linear_function->add_variable( & v_primary_spinning_reserve[ t ], - 1.0 );

  PrimaryRho_Constraints[ t ].set_lhs( 0.0 );
  PrimaryRho_Constraints[ t ].set_rhs( Inf<double>() );
  PrimaryRho_Constraints[ t ].set_function( linear_function );
 }

 add_static_constraint( PrimaryRho_Constraints );

 // Initializing secondary rho fraction constraints

 SecondaryRho_Constraints.resize( f_time_horizon );

 for( Index t = 0; t < f_time_horizon; ++t ) {

  auto linear_function = new LinearFunction();
  linear_function->add_variable( & v_active_power[ t ],
                                 v_SecondaryRho[ t ] );
  linear_function->add_variable( & v_secondary_spinning_reserve[ t ],
                                 - 1.0 );

  SecondaryRho_Constraints[ t ].set_lhs( 0.0 );
  SecondaryRho_Constraints[ t ].set_rhs( Inf<double>() );
  SecondaryRho_Constraints[ t ].set_function( linear_function );
 }

 add_static_constraint( SecondaryRho_Constraints );

/*--------------------------------------------------------------------------*/
 // TIME DEPENDENT START UP COSTS CONSTRAINTS

 //TODO if any exist

} // end( ThermalUnitBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::generate_objective( Configuration *objc ) {

 if( ! get_objective().empty() )  // an objective is there already
  return;                        // cowardly (and silently) return

 // Initialize objective function

 if( v_commitment.size() != f_time_horizon ) {
  throw( std::logic_error
          ( "ThermalUnitBlock::generate_objective: v_commitment must have "
            "size equal to the time horizon." ) );
 }

 auto dquad_function = new DQuadFunction();

 for( Index t = init_t; t < f_time_horizon; ++t ) {
  dquad_function->add_variable( & start_up( t ), v_StartUpCost[ t ], 0.0 );
 }

 for( Index t = 0; t < f_time_horizon; ++t ) {
  dquad_function->add_variable
          ( & v_active_power[ t ], v_LinearTerm[ t ], v_QuadTerm[ t ] );
  dquad_function->add_variable( & v_commitment[ t ], v_ConstTerm[ t ], 0.0 );
 }

 objective.set_function( dquad_function );
 objective.set_sense( Objective::eMin );
 objective.set_Block( this );

}  // end( ThermalUnitBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*-------- METHODS FOR LOADING, PRINTING & SAVING THE ThermalUnitBlock -----*/
/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::serialize( netCDF::NcGroup & group ) const {

 UnitBlock::serialize( group );

 group.putAtt( "type" , "ThermalUnitBlock" );
 group.addDim( "TimeHorizon" , f_time_horizon );

 ::serialize( group, "InitialPower",   netCDF::NcDouble(), f_initial_power );
 ::serialize( group, "MinUpTime",      netCDF::NcUint64(), f_MinUpTime );
 ::serialize( group, "MinDownTime",    netCDF::NcUint64(), f_MinDownTime );
 ::serialize( group, "InitUpDownTime", netCDF::NcUint64(), f_InitUpDownTime );

 // ::serialize( group, "ShutDownCapability", netCDF::NcDouble(),
 //              f_shut_down_capability );
 // ::serialize( group, "StartUpCapability", netCDF::NcDouble(),
 //             f_start_up_capability );

 ::serialize( group, "InitialMinPower",
              netCDF::NcDouble(), f_initial_min_power );
 ::serialize( group, "InitialDeltaRampUp",
              netCDF::NcDouble(), f_initial_delta_ramp_up );
 ::serialize( group, "InitialDeltaRampDown",
              netCDF::NcDouble(), f_initial_delta_ramp_down );

 auto NumberIntervals = group.getDim( "NumberIntervals" );

 ::serialize( group, "MinPower", netCDF::NcDouble(),
              NumberIntervals, v_MinPower );

 ::serialize( group, "MaxPower", netCDF::NcDouble(),
              NumberIntervals, v_MaxPower );

 ::serialize( group, "DeltaRampUp", netCDF::NcDouble(),
              NumberIntervals, v_DeltaRampUp );

 ::serialize( group, "DeltaRampDown", netCDF::NcDouble(),
              NumberIntervals, v_DeltaRampDown );

 ::serialize( group, "PrimaryRho", netCDF::NcDouble(),
              NumberIntervals, v_PrimaryRho );

 ::serialize( group, "SecondaryRho", netCDF::NcDouble(),
              NumberIntervals, v_SecondaryRho );

 ::serialize( group, "QuadTerm", netCDF::NcDouble(),
              NumberIntervals, v_QuadTerm );

 ::serialize( group, "LinearTerm", netCDF::NcDouble(),
              NumberIntervals, v_LinearTerm );

 ::serialize( group, "ConstTerm", netCDF::NcDouble(),
              NumberIntervals, v_ConstTerm );
 ::serialize( group, "StartUpCost",    netCDF::NcDouble(),
              NumberIntervals, v_StartUpCost );

}  // end( ThermalUnitBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------------- End File ThermalUnitBlock.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
