/*--------------------------------------------------------------------------*/
/*--------------------- File ThermalUnitBlock.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the ThermalUnitBlock class.
 *
 * \version 0.11
 *
 * \date 03 - 10 - 2020
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

ThermalUnitBlock::~ThermalUnitBlock() {
 auto clear_constraints =
  []( std::vector< FRowConstraint > & constraints ) {
   for( auto & constraint : constraints )
    constraint.clear();
  };

 clear_constraints( Power_StartUp_ShutDown_Variables_Constraints );
 clear_constraints( Power_StartUp_Variable_Constraints );
 clear_constraints( Power_ShutDown_Variable_Constraints );
 clear_constraints( StartUp_ShutDown_Variables_Constraints );
 clear_constraints( StartUp_Constraints );
 clear_constraints( ShutDown_Constraints );
 clear_constraints( RampUp_Constraints );
 clear_constraints( RampDown_Constraints );
 clear_constraints( PrimaryRho_Constraints );
 clear_constraints( SecondaryRho_Constraints );
 clear_constraints( MinPower_Constraints );
 clear_constraints( MaxPower_Constraints );
}

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::deserialize( netCDF::NcGroup & group ) {

#ifndef NDEBUG
 std::cerr << "[DEBUG] ThermalUnitBlock::deserialize() - Checking Dims"
           << std::endl;
 std::vector< std::string > expected_dims = { "TimeHorizon",
                                              "NumberIntervals" };
 check_dimensions( group, expected_dims, std::cerr );

 std::cerr << "[DEBUG] ThermalUnitBlock::deserialize() - Checking Vars"
           << std::endl;
 std::vector< std::string > expected_vars = { "MinPower",
                                              "MaxPower",
                                              "DeltaRampUp",
                                              "DeltaRampDown",
                                              "PrimaryRho",
                                              "SecondaryRho",
                                              "LinearTerm",
                                              "QuadTerm",
                                              "ConstTerm",
                                              "StartUpCost",
                                              "FixedConsumption",
                                              "InertiaCommitment",
                                              "InitialPower",
                                              "MinUpTime",
                                              "MinDownTime",
                                              "InitUpDownTime" };
 check_variables( group, expected_vars, std::cerr );
#endif

 UnitBlock::deserialize_time_horizon( group );
 UnitBlock::deserialize_change_intervals( group );

 ::deserialize( group, "MinPower", f_time_horizon, v_MinPower, true, true );
 ::deserialize( group, "MaxPower", f_time_horizon, v_MaxPower, true, true );
 ::deserialize( group, "DeltaRampUp", f_time_horizon, v_DeltaRampUp, true, true );
 ::deserialize( group, "DeltaRampDown", f_time_horizon, v_DeltaRampDown, true, true );
 ::deserialize( group, "PrimaryRho", f_time_horizon, v_PrimaryRho, true, true );
 ::deserialize( group, "SecondaryRho", f_time_horizon, v_SecondaryRho, true, true );
 ::deserialize( group, "LinearTerm", f_time_horizon, v_LinearTerm, true, true );
 ::deserialize( group, "QuadTerm", f_time_horizon, v_QuadTerm, true, true );
 ::deserialize( group, "ConstTerm", f_time_horizon, v_ConstTerm, true, true );
 ::deserialize( group, "StartUpCost", f_time_horizon, v_StartUpCost, true, true );
 ::deserialize( group, "FixedConsumption", f_time_horizon, v_fixed_consumption, true, true );
 ::deserialize( group, "InertiaCommitment", f_time_horizon, v_inertia_commitment, true, true );

 ::deserialize( group, "InitialPower", &f_initial_power );
 ::deserialize( group, "MinUpTime", &f_MinUpTime );
 ::deserialize( group, "MinDownTime", &f_MinDownTime );
 ::deserialize( group, "InitUpDownTime", &f_InitUpDownTime );

 if( ! ::deserialize( group, "Availability", f_time_horizon,
                      v_Availability, true, true ) ) {
  v_Availability.resize( get_time_horizon() , 1.0 );
 }

 decompress_vector( v_MinPower );
 decompress_vector( v_MaxPower );
 decompress_vector( v_Availability );
 if( ! v_DeltaRampUp.empty() ) {
  decompress_vector( v_DeltaRampUp );
 }
 if( ! v_DeltaRampDown.empty() ) {
  decompress_vector( v_DeltaRampDown );
 }
 decompress_vector( v_PrimaryRho );
 decompress_vector( v_SecondaryRho );
 decompress_vector( v_LinearTerm );
 decompress_vector( v_QuadTerm );
 decompress_vector( v_ConstTerm );
 decompress_vector( v_StartUpCost );
 decompress_vector( v_fixed_consumption );
 decompress_vector( v_inertia_commitment );

 UnitBlock::deserialize( group );
}  // end( ThermalUnitBlock::deserialize )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::generate_abstract_variables( Configuration * stvv ) {

 if( variables_generated() )
  return; // variables have already been generated

 UnitBlock::generate_abstract_variables( stvv );

 if( f_InitUpDownTime > 0 ) {
  init_t = ( f_InitUpDownTime >= f_MinUpTime ? 0 :
             f_MinUpTime - f_InitUpDownTime );
 } else {
  init_t = ( -f_InitUpDownTime >= f_MinDownTime ? 0 :
             f_MinDownTime + f_InitUpDownTime );
 }

/*--------------------------------------------------------------------------*/

 // Commitment Variable

 v_commitment.resize( f_time_horizon );
 for( auto & var : v_commitment )
  var.set_type( ColVariable::kBinary );
 add_static_variable( v_commitment, "u" );

 // Active Power Variable

 v_active_power.resize( f_time_horizon );
 for( auto & var : v_active_power )
  var.set_type( ColVariable::kNonNegative );
 add_static_variable( v_active_power, "p" );

 // Primary Spinning Reserve Variable

 v_primary_spinning_reserve.resize( f_time_horizon );
 for( auto & var : v_primary_spinning_reserve )
  var.set_type( ColVariable::kNonNegative );
 add_static_variable( v_primary_spinning_reserve, "pr" );

 // Secondary Spinning Reserve Variable

 v_secondary_spinning_reserve.resize( f_time_horizon );
 for( auto & var : v_secondary_spinning_reserve )
  var.set_type( ColVariable::kNonNegative );
 add_static_variable( v_secondary_spinning_reserve, "sr" );

/*--------------------------------------------------------------------------*/
 // START UP AND SHUT DOWN BINARY VARIABLES
 auto startup_shutdown_size = f_time_horizon - init_t;

 if( startup_shutdown_size > 0 ) {

  v_start_up.resize( startup_shutdown_size );
  for( auto & var : v_start_up )
   var.set_type( ColVariable::kBinary );
  add_static_variable( v_start_up, "v" );

  v_shut_down.resize( startup_shutdown_size );
  for( auto & var : v_shut_down )
   var.set_type( ColVariable::kBinary );
  add_static_variable( v_shut_down, "w" );

 }

/*--------------------------------------------------------------------------*/
 // POSSIBLY FIXING THE COMMITMENT VARIABLES TO 0 OR 1

 if( init_t > 0 ) {

  double commitment_variable_value = -1.0;

  if( f_InitUpDownTime <= 0 && -f_InitUpDownTime < f_MinDownTime ) {

   commitment_variable_value = 0.0;

   // This unit must remain off from time 0 to init_t - 1. Therefore,
   // it should not produce any power.

   if( ! v_active_power.empty() ) {
    for( Index t = 0 ; t < init_t ; ++t ) {
     v_active_power[ t ].set_value( 0.0 );
     v_active_power[ t ].is_fixed( true );
    }
   }

   if( ! v_primary_spinning_reserve.empty() ) {
    for( Index t = 0 ; t < init_t ; ++t ) {
     v_primary_spinning_reserve[ t ].set_value( 0.0 );
     v_primary_spinning_reserve[ t ].is_fixed( true );
    }
   }

   if( ! v_secondary_spinning_reserve.empty() ) {
    for( Index t = 0 ; t < init_t ; ++t ) {
     v_secondary_spinning_reserve[ t ].set_value( 0.0 );
     v_secondary_spinning_reserve[ t ].is_fixed( true );
    }
   }
  } else if( f_InitUpDownTime > 0 && f_InitUpDownTime < f_MinUpTime ) {
   // This unit must remain on from time 0 to init_t - 1.
   commitment_variable_value = 1.0;

  }

  if( ! v_commitment.empty() ) {
   if( commitment_variable_value >= 0.0 ) {
    // Fixing the commitment variables v_commitment to 0 or 1 from
    // time 0 to init_t - 1, depending on whether this unit must
    // remain off or on for the first init_t time steps.
    for( Index t = 0 ; t < init_t ; ++t ) {
     v_commitment[ t ].set_value( commitment_variable_value );
     v_commitment[ t ].is_fixed( true );
    }
   }
  }

  if( f_InitUpDownTime > 0 ) {

   if( f_MinDownTime < startup_shutdown_size ) {
    for( Index t = 0 ; t < f_MinDownTime ; ++t ) {
     v_start_up[ t ].set_value( 0.0 );
     v_start_up[ t ].is_fixed( true );
    }
   }
   else if ( f_MinDownTime >= startup_shutdown_size ) {
    for( Index t = 0 ; t < startup_shutdown_size ; ++t ) {
     v_start_up[ t ].set_value( 0.0 );
     v_start_up[ t ].is_fixed( true );
    }
   }
  }
  else if( f_InitUpDownTime <= 0 ) {
   if( f_MinUpTime < startup_shutdown_size ) {
    for( Index t = 0 ; t < f_MinUpTime ; ++t ) {
     v_shut_down[ t ].set_value( 0.0 );
     v_shut_down[ t ].is_fixed( true );
    }
   }
   else if( f_MinUpTime >= startup_shutdown_size ) {
    for( Index t = 0 ; t < startup_shutdown_size ; ++t ) {
     v_shut_down[ t ].set_value( 0.0 );
     v_shut_down[ t ].is_fixed( true );
    }
   }
  }
 }

 else if ( init_t == 0 ) {

  if( f_InitUpDownTime > 0 ) {

   if( f_MinDownTime < f_time_horizon ) {
    for( Index t = 0 ; t < f_MinDownTime ; ++t ) {
     v_start_up[ t ].set_value( 0.0 );
     v_start_up[ t ].is_fixed( true );
    }
   }
   else if( f_MinDownTime >= f_time_horizon ) {
    for( Index t = 0 ; t < f_time_horizon ; ++t ) {
     v_start_up[ t ].set_value( 0.0 );
     v_start_up[ t ].is_fixed( true );
    }
   }
  }
  if( f_InitUpDownTime <= 0 ) {
   if( f_MinUpTime < f_time_horizon ) {
    for( Index t = 0 ; t < f_MinUpTime ; ++t ) {
     v_shut_down[ t ].set_value( 0.0 );
     v_shut_down[ t ].is_fixed( true );
    }
   }
   else if( f_MinUpTime >= f_time_horizon ) {
    for( Index t = 0 ; t < f_time_horizon ; ++t ) {
     v_shut_down[ t ].set_value( 0.0 );
     v_shut_down[ t ].is_fixed( true );
    }
   }
  }
 }

 set_variables_generated();
} // end( ThermalUnitBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::generate_abstract_constraints( Configuration * stcc ) {

 if( constraints_generated() )
  return; // constraints have already been generated

 // Initializing start up and shut down variables connection constraints
 if( init_t == 0 ) {

  if( f_time_horizon > 0 ) {

   StartUp_ShutDown_Variables_Constraints.resize( f_time_horizon );

   // Initial condition

   auto linear_function = new LinearFunction();

   linear_function->add_variable( & v_commitment[ 0 ] ,  1.0 );
   linear_function->add_variable( & v_start_up[ 0 ]   , -1.0 );
   linear_function->add_variable( & v_shut_down[ 0 ]  ,  1.0 );

   auto initial_commitment = ( f_InitUpDownTime > 0 ? 1.0 : 0.0 );

   StartUp_ShutDown_Variables_Constraints[ 0 ].set_both( initial_commitment );
   StartUp_ShutDown_Variables_Constraints[ 0 ].set_function( linear_function );

   for( Index t = 1 ; t < f_time_horizon ; ++t ) {

    auto lf = new LinearFunction();

    lf->add_variable( & v_commitment[ t ]     ,  1.0 );
    lf->add_variable( & v_start_up[ t ]       , -1.0 );
    lf->add_variable( & v_shut_down[ t ]      ,  1.0 );
    lf->add_variable( & v_commitment[ t - 1 ] , -1.0 );
    StartUp_ShutDown_Variables_Constraints[ t ].set_both( 0.0 );
    StartUp_ShutDown_Variables_Constraints[ t ].set_function( lf );
   }
  }
 }

 if( init_t > 0 ) {

  auto startup_shutdown_const_size = static_cast<int>( f_time_horizon - init_t );

  if( startup_shutdown_const_size > 0 ) {

   StartUp_ShutDown_Variables_Constraints.resize( startup_shutdown_const_size );

   // Initial condition

   auto l_function = new LinearFunction();

   l_function->add_variable( & v_commitment[ init_t ] ,  1.0 );
   l_function->add_variable( & v_start_up[ 0 ]        , -1.0 );
   l_function->add_variable( & v_shut_down[ 0 ]       ,  1.0 );

   if( f_InitUpDownTime < 0 && -f_InitUpDownTime < f_MinDownTime ) {
    StartUp_ShutDown_Variables_Constraints[ 0 ].set_both( 0.0 );
    StartUp_ShutDown_Variables_Constraints[ 0 ].set_function( l_function );
   }
   else if( f_InitUpDownTime > 0 && f_InitUpDownTime < f_MinUpTime ) {
    StartUp_ShutDown_Variables_Constraints[ 0 ].set_both( 1.0 );
    StartUp_ShutDown_Variables_Constraints[ 0 ].set_function( l_function );
   }

   for( Index t = init_t + 1 , constraint_index = 1 ; t < f_time_horizon ;
        ++t , ++constraint_index ) {

    auto linear_function = new LinearFunction();

    linear_function->add_variable( & v_commitment[ t ]         ,  1.0 );
    linear_function->add_variable( & v_start_up[ t - init_t ]  , -1.0 );
    linear_function->add_variable( & v_shut_down[ t - init_t ] ,  1.0 );
    linear_function->add_variable( & v_commitment[ t - 1 ]     , -1.0 );
    StartUp_ShutDown_Variables_Constraints[ constraint_index ].set_both( 0.0 );
    StartUp_ShutDown_Variables_Constraints[ constraint_index ].
     set_function( linear_function );
   }
  }
 }
 add_static_constraint( StartUp_ShutDown_Variables_Constraints,
                        "StartUp_ShutDown_Commitment_Constraints" );

 // Initializing turn on constraints (start up constraints)

 auto startup_const_size = static_cast<int>
  ( f_time_horizon - init_t - f_MinUpTime );

 if( startup_const_size > 0 ) {

  StartUp_Constraints.resize( startup_const_size );

  for( Index t = init_t + f_MinUpTime , constraint_index = 0;
       t < f_time_horizon ; ++t , ++constraint_index ) {

   auto linear_function = new LinearFunction();

   for( Index s = t - f_MinUpTime ; s < t ; ++s ) {
    linear_function->add_variable( & v_start_up[ s - init_t + 1 ] , -1.0 );
   }

   linear_function->add_variable( & v_commitment[ t ] , 1.0 );
   StartUp_Constraints[ constraint_index ].set_lhs( 0.0 );
   StartUp_Constraints[ constraint_index ].set_rhs( Inf< double >());
   StartUp_Constraints[ constraint_index ].set_function( linear_function );
  }

  add_static_constraint( StartUp_Constraints , "StartUp_Commitment_Constraints" );
 }

 // Initializing turn off constraints (shut down constraints)

 auto shutdown_const_size = static_cast<int>
  ( f_time_horizon - init_t - f_MinDownTime );

 if( shutdown_const_size > 0 ) {

  ShutDown_Constraints.resize( shutdown_const_size );

  for( Index t = init_t + f_MinDownTime , constraint_index = 0 ;
       t < f_time_horizon ; ++t , ++constraint_index ) {

   auto linear_function = new LinearFunction();

   for( Index s = t - f_MinDownTime ; s < t ; ++s ) {
    linear_function->add_variable( & v_shut_down[ s - init_t + 1 ] , 1.0 );
   }

   linear_function->add_variable( & v_commitment[ t ] , 1.0 );
   ShutDown_Constraints[ constraint_index ].set_lhs( 0.0 );
   ShutDown_Constraints[ constraint_index ].set_rhs( 1.0 );
   ShutDown_Constraints[ constraint_index ].set_function( linear_function );
  }

  add_static_constraint( ShutDown_Constraints , "ShutDown_Commitment_Constraints" );
 }

/*--------------------------------------------------------------------------*/
 // Initializing ramp-up constraints with 3-Binary Variables

 if( ! v_DeltaRampUp.empty() && ! v_DeltaRampDown.empty() ) {

  if( f_InitUpDownTime > 0 ) {
   for( Index t = 0; t < f_time_horizon; ++t ) {
    if( f_initial_power + v_DeltaRampUp[0] < get_operational_min_power( 0 ) ||
        f_initial_power - v_DeltaRampDown[0] > get_operational_max_power( 0 )) {
     throw ( std::logic_error
             ( "ThermalUnitBlock::Ramp Constraints: when f_InitUpDownTime > 0,"
               "it must be that"
               "f_initial_power + v_DeltaRampUp[ 0 ] >= get_operational_min_power( 0 )"
               "f_initial_power - v_DeltaRampDown[ 0 ] <= get_operational_max_power( 0 ) " ));
    }
   }
  }
 }

 if( ! v_DeltaRampUp.empty() ) {

  RampUp_Constraints.resize( f_time_horizon );

  // Initial condition

  if( init_t == 0 ) {

   // Initial condition

   auto linear_function = new LinearFunction();
   linear_function->add_variable( &v_active_power[0], 1.0 );
   linear_function->add_variable( &v_start_up[0],
                                  -get_operational_min_power( 0 ));
   auto initial_commitment = ( f_InitUpDownTime > 0 ? 1.0 : 0.0 );

   RampUp_Constraints[0].set_lhs( 0.0 );

   if( f_InitUpDownTime > 0 ) {
    RampUp_Constraints[0].set_rhs
            (( v_DeltaRampUp[0] * initial_commitment ) + f_initial_power );
   } else if( f_InitUpDownTime <= 0 ) {
    RampUp_Constraints[0].set_rhs
            (( v_DeltaRampUp[0] * initial_commitment ));
   }
   RampUp_Constraints[0].set_function( linear_function );

   for( Index t = 1; t < f_time_horizon; ++t ) {

    auto lf = new LinearFunction();

    lf->add_variable( &v_active_power[t], -1.0 );
    lf->add_variable( &v_active_power[t - 1], 1.0 );
    lf->add_variable( &v_start_up[t], get_operational_min_power( t ));
    lf->add_variable( &v_commitment[t - 1], v_DeltaRampUp[t] );

    RampUp_Constraints[t].set_lhs( 0.0 );
    RampUp_Constraints[t].set_rhs( Inf< double >());
    RampUp_Constraints[t].set_function( lf );
   }
  } else if( init_t > 0 ) {

   auto linear_function = new LinearFunction();
   linear_function->add_variable( &v_active_power[0], 1.0 );

   auto initial_commitment = ( f_InitUpDownTime > 0 ? 1.0 : 0.0 );

   RampUp_Constraints[0].set_lhs( 0.0 );

   if( f_InitUpDownTime > 0 ) {
    RampUp_Constraints[0].set_rhs
            (( v_DeltaRampUp[0] * initial_commitment ) + f_initial_power );
   } else if( f_InitUpDownTime <= 0 ) {
    RampUp_Constraints[0].set_rhs( v_DeltaRampUp[0] * initial_commitment );
   }
   RampUp_Constraints[0].set_function( linear_function );

   for( Index t = 1; t < init_t; ++t ) {

    auto lf = new LinearFunction();

    lf->add_variable( &v_active_power[t], -1.0 );
    lf->add_variable( &v_active_power[t - 1], 1.0 );
    lf->add_variable( &v_commitment[t - 1], v_DeltaRampUp[t] );

    RampUp_Constraints[t].set_lhs( 0.0 );
    RampUp_Constraints[t].set_rhs( Inf< double >());
    RampUp_Constraints[t].set_function( lf );
   }

   if( f_InitUpDownTime < 0 && -f_InitUpDownTime < f_MinDownTime ) {

    auto LFunction = new LinearFunction();
    LFunction->add_variable( &v_active_power[init_t], 1.0 );
    LFunction->add_variable( &v_start_up[0],
                             -get_operational_min_power( init_t ));

    RampUp_Constraints[init_t].set_lhs( -Inf< double >());
    RampUp_Constraints[init_t].set_rhs( 0.0 );
    RampUp_Constraints[init_t].set_function( LFunction );

   } else if( f_InitUpDownTime > 0 && f_InitUpDownTime < f_MinUpTime ) {

    auto LFunction = new LinearFunction();
    LFunction->add_variable( &v_active_power[init_t], 1.0 );
    LFunction->add_variable( &v_active_power[init_t - 1], -1.0 );
    LFunction->add_variable( &v_start_up[0],
                             -get_operational_min_power( init_t ));

    RampUp_Constraints[init_t].set_lhs( -Inf< double >());
    RampUp_Constraints[init_t].set_rhs( v_DeltaRampUp[init_t] );
    RampUp_Constraints[init_t].set_function( LFunction );
   }

   for( Index t = init_t + 1; t < f_time_horizon; ++t ) {

    auto lf = new LinearFunction();

    lf->add_variable( &v_active_power[t], -1.0 );
    lf->add_variable( &v_active_power[t - 1], 1.0 );
    lf->add_variable( &v_start_up[t - init_t],
                      get_operational_min_power( t ));
    lf->add_variable( &v_commitment[t - 1], v_DeltaRampUp[t] );

    RampUp_Constraints[t].set_lhs( 0.0 );
    RampUp_Constraints[t].set_rhs( Inf< double >());
    RampUp_Constraints[t].set_function( lf );
   }
  }

  add_static_constraint( RampUp_Constraints, "RampUp_Constraints_Thermal" );
 }
  // Initializing ramp down constraints
 if( ! v_DeltaRampDown.empty() ) {

  RampDown_Constraints.resize( f_time_horizon );

  // Initial condition
  if( init_t == 0 ) {
   auto linear_function = new LinearFunction();
   linear_function->add_variable( & v_active_power[ 0 ], 1.0 );
   linear_function->add_variable( & v_commitment[ 0 ] , v_DeltaRampDown[ 0 ] );
   linear_function->add_variable( & v_shut_down[ 0 ] ,
                                  get_operational_min_power( 0 ) );
   if( f_InitUpDownTime > 0 ) {
    RampDown_Constraints[ 0 ].set_lhs( f_initial_power );
   }
   else if( f_InitUpDownTime <= 0 ){
    RampDown_Constraints[ 0 ].set_lhs( 0.0 );
   }
   RampDown_Constraints[ 0 ].set_rhs( Inf< double >() );
   RampDown_Constraints[ 0 ].set_function( linear_function );

   // Remaining constraints
   for( Index t = 1 ; t < f_time_horizon ; ++t ) {

    auto lf = new LinearFunction();

    lf->add_variable( & v_active_power[ t ]     , 1.0 );
    lf->add_variable( & v_active_power[ t - 1 ] , -1.0 );
    lf->add_variable( & v_shut_down[ t ] , get_operational_min_power( t ) );
    lf->add_variable( & v_commitment[ t ] , v_DeltaRampDown[ t ] );

    RampDown_Constraints[ t ].set_lhs( 0.0 );
    RampDown_Constraints[ t ].set_rhs( Inf< double >() );
    RampDown_Constraints[ t ].set_function( lf );
   }
  }

  if( init_t > 0 ) {
   // Initial condition
   auto linear_function = new LinearFunction();
   linear_function->add_variable( & v_active_power[ 0 ] , 1.0 );
   linear_function->add_variable( & v_commitment[ 0 ] , v_DeltaRampDown[ 0 ] );
   if ( f_InitUpDownTime > 0 ) {
    RampDown_Constraints[ 0 ].set_lhs( f_initial_power );
   }
   else if ( f_InitUpDownTime <= 0 ) {
    RampDown_Constraints[ 0 ].set_lhs( 0.0 );
   }
   RampDown_Constraints[ 0 ].set_rhs( Inf< double >() );
   RampDown_Constraints[ 0 ].set_function( linear_function );

   for( Index t = 1 ; t < init_t ; ++t ) {
    auto lf = new LinearFunction();

    lf->add_variable( & v_active_power[ t - 1 ] , -1.0 );
    lf->add_variable( & v_active_power[ t ]     ,  1.0 );
    lf->add_variable( & v_commitment[ t ] , v_DeltaRampDown[ t ] );
    RampDown_Constraints[ t ].set_lhs( 0.0 );
    RampDown_Constraints[ t ].set_rhs( Inf< double >() );
    RampDown_Constraints[ t ].set_function( lf );
   }

   // Remaining constraints

   for( Index t = init_t ; t < f_time_horizon ; ++t ) {

    auto lf = new LinearFunction();

    lf->add_variable( & v_active_power[ t ]       ,  1.0 );
    lf->add_variable( & v_active_power[ t - 1 ]   , -1.0 );
    lf->add_variable( & v_shut_down[ t - init_t ] ,
                      get_operational_min_power( t ) );
    lf->add_variable( & v_commitment[ t ]         , v_DeltaRampDown[ t ] );

    RampDown_Constraints[ t ].set_lhs( 0.0 );
    RampDown_Constraints[ t ].set_rhs( Inf< double >() );
    RampDown_Constraints[ t ].set_function( lf );
   }
  }

  add_static_constraint( RampDown_Constraints, "RampDown_Constraints_Thermal" );

 }

/*--------------------------------------------------------------------------*/

 // Initializing minimum power constraints

 MinPower_Constraints.resize( f_time_horizon );

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {

  auto linear_function = new LinearFunction();

  linear_function->add_variable( & v_active_power[ t ] , 1.0 );
  linear_function->add_variable( & v_primary_spinning_reserve[ t ] , -1.0 );
  linear_function->add_variable( & v_secondary_spinning_reserve[ t ] , -1.0 );
  linear_function->add_variable( & v_commitment[ t ] ,
                                 - get_operational_min_power( t ) );

  MinPower_Constraints[t].set_rhs( Inf< double >() );
  MinPower_Constraints[t].set_lhs( 0.0 );
  MinPower_Constraints[t].set_function( linear_function );
 }
 add_static_constraint( MinPower_Constraints, "MinPower_Constraints_Thermal" );

 // Initializing maximum power constraints
 MaxPower_Constraints.resize( f_time_horizon );

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {

  auto linear_function = new LinearFunction();

  linear_function->add_variable( & v_active_power[ t ] , -1.0 );
  linear_function->add_variable( & v_primary_spinning_reserve[ t ] , -1.0 );
  linear_function->add_variable( & v_secondary_spinning_reserve[ t ] , -1.0 );

  linear_function->add_variable( & v_commitment[ t ] ,
                                 get_operational_max_power( t ) );

  MaxPower_Constraints[ t ].set_lhs( 0.0 );
  MaxPower_Constraints[ t ].set_rhs( Inf< double >() );
  MaxPower_Constraints[ t ].set_function( linear_function );
 }
 add_static_constraint( MaxPower_Constraints, "MaxPower_Constraints_Thermal" );


 // Initializing primary rho fraction constraints

  PrimaryRho_Constraints.resize( f_time_horizon );

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {

   auto linear_function = new LinearFunction();

   if( !v_PrimaryRho.empty() ) {

    linear_function->add_variable( &v_active_power[t], v_PrimaryRho[t] );
   } else {
    linear_function->add_variable( &v_active_power[t], 0.0);

   }
   linear_function->add_variable( & v_primary_spinning_reserve[ t ] , -1.0 );

   PrimaryRho_Constraints[ t ].set_lhs( 0.0 );
   PrimaryRho_Constraints[ t ].set_rhs( Inf< double >() );
   PrimaryRho_Constraints[ t ].set_function( linear_function );
  }

  add_static_constraint( PrimaryRho_Constraints, "PrimaryRho_Constraints_Thermal" );

 // Initializing secondary rho fraction constraints

  SecondaryRho_Constraints.resize( f_time_horizon );

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {

   auto linear_function = new LinearFunction();
   if( !v_SecondaryRho.empty()  ) {
    linear_function->add_variable( &v_active_power[t], v_SecondaryRho[t] );
   }else{
    linear_function->add_variable( &v_active_power[t], 0.0);
   }
   linear_function->add_variable( & v_secondary_spinning_reserve[ t ], -1.0 );

   SecondaryRho_Constraints[ t ].set_lhs( 0.0 );
   SecondaryRho_Constraints[ t ].set_rhs( Inf< double >() );
   SecondaryRho_Constraints[ t ].set_function( linear_function );
  }

  add_static_constraint( SecondaryRho_Constraints, "SecondaryRho_Constraints_Thermal" );


 set_constraints_generated();
} // end( ThermalUnitBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::generate_objective( Configuration * objc ) {

 if( objective_generated() )
  return; // Objective has already been generated

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
  dquad_function->add_variable( &v_start_up[ t - init_t ] ,
                                start_up_cost[ t - init_t ] ,
                                0.0 );
 }

 for( Index t = 0; t < f_time_horizon; ++t ) {
  dquad_function->add_variable
   ( & v_active_power[ t ] , linear_term[ t ] , quad_term[ t ] );

  dquad_function->add_variable
   ( & v_commitment[ t ] , const_term[ t ] , 0.0 );
 }

 objective.set_function( dquad_function );
 objective.set_sense( Objective::eMin );

 // Set Block objective
 this->set_objective( &objective );

 set_objective_generated();
}  // end( ThermalUnitBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*-------- METHODS FOR LOADING, PRINTING & SAVING THE ThermalUnitBlock -----*/
/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::serialize( netCDF::NcGroup & group ) const {

 UnitBlock::serialize( group );

 ::serialize( group, "InitialPower", netCDF::NcDouble(), f_initial_power );
 ::serialize( group, "MinUpTime", netCDF::NcUint(), f_MinUpTime );
 ::serialize( group, "MinDownTime", netCDF::NcUint(), f_MinDownTime );
 ::serialize( group, "InitUpDownTime", netCDF::NcInt(), f_InitUpDownTime );

 auto NumberIntervals = group.getDim( "NumberIntervals" );

 ::serialize( group, "MinPower", netCDF::NcDouble(),
              NumberIntervals, v_MinPower, true );

 ::serialize( group, "MaxPower", netCDF::NcDouble(),
              NumberIntervals, v_MaxPower, true );

 ::serialize( group, "Availability", netCDF::NcDouble(),
              NumberIntervals, v_Availability, true );

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
              { NumberIntervals }, v_fixed_consumption, true );

 ::serialize( group, "InertiaCommitment", netCDF::NcDouble(),
              { NumberIntervals }, v_inertia_commitment, true );
}

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::update_availability_dependents
( Index t , c_ModParam issueAMod ) {

 if( ! constraints_generated() )
  return;

 // MaxPower_Constraints
 {
  auto f = static_cast<LinearFunction *>( MaxPower_Constraints[ t ].
                                          get_function() );
  auto var_index = f->is_active( & v_commitment[ t ] );
  assert( var_index < f->get_num_active_var() );
  f->modify_coefficient( var_index , get_operational_max_power( t ) ,
                         issueAMod );
 }

 // MinPower_Constraints
 {
  auto f = static_cast<LinearFunction *>( MinPower_Constraints[ t ].
                                          get_function() );
  auto var_index = f->is_active( & v_commitment[ t ] );
  assert( var_index < f->get_num_active_var() );
  f->modify_coefficient( var_index , - get_operational_min_power( t ) ,
                         issueAMod );
 }

 // RampUp_Constraints

 if( init_t == 0 ) {

  double coefficient = get_operational_min_power( t );
  if( t == 0 )
   coefficient *= -1.0;

  auto f = static_cast<LinearFunction *>( RampUp_Constraints[ t ].
                                          get_function() );
  auto var_index = f->is_active( & v_start_up[ t ] );
  assert( var_index < f->get_num_active_var() );
  f->modify_coefficient( var_index , coefficient , issueAMod );
 }
 else if( init_t > 0 ) {

  auto depends_on_min_power = ( t > init_t );
  depends_on_min_power |= ( t == init_t ) &&
   ( ( f_InitUpDownTime < 0 && -f_InitUpDownTime < f_MinDownTime ) ||
     ( f_InitUpDownTime > 0 &&  f_InitUpDownTime < f_MinUpTime ) );

  if( depends_on_min_power ) {

   auto coefficient = get_operational_min_power( t );
   if( t == init_t )
    coefficient *= -1.0;

   auto f = static_cast<LinearFunction *>( RampUp_Constraints[ t ].
                                           get_function() );
   auto var_index = f->is_active( & v_start_up[ t - init_t ] );
   assert( var_index < f->get_num_active_var() );
   f->modify_coefficient( var_index , coefficient , issueAMod );
  }
 }

 // RampDown_Constraints
 if( ( init_t == 0 && t == 0 ) || ( init_t > 0 && t >= init_t ) ) {

  auto f = static_cast<LinearFunction *>( RampDown_Constraints[ t ].
                                          get_function() );
  auto var_index = f->is_active( & v_shut_down[ t - init_t ] );
  assert( var_index < f->get_num_active_var() );
  auto coefficient = get_operational_min_power( t );
  f->modify_coefficient( var_index , coefficient , issueAMod );
 }
}  // end( ThermalUnitBlock::update_availability_dependents )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_availability
( std::vector< double >::const_iterator values, Block::Subset && subset,
  const bool ordered, c_ModParam issuePMod, c_ModParam issueAMod ) {

 if( subset.empty() ) {
  return;
 }

 if( v_Availability.empty() ) {
  if( std::all_of( values,
                   values + subset.size(),
                   []( double cst ) {
                    return ( cst == 1.0 );
                   } ) ) {
   return;
  }

  v_Availability.assign( get_time_horizon() , 1.0 );
 }

 // If nothing changes, return
 bool identical = true;
 auto availability = values;
 for( auto t : subset ) {

  if( t >= v_Availability.size() ) {
   throw( std::invalid_argument
          ( "ThermalUnitBlock::set_availability: invalid index in subset: "
            + std::to_string( t ) ) );
  }

  // Check change

  if( v_Availability[ t ] != *availability ) {
   identical = false;
  }

  // Check consistency

  if( ! availability_is_consistent( t , *availability ) )
   throw( std::logic_error
          ( "ThermalUnitBlock::set_availability: availability (" +
            std::to_string( *availability ) + ") at time " +
            std::to_string( t ) + " is not consistent." ) );

  std::advance( availability , 1 );
 }

 if( identical )
  return; // nothing changes

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  availability = values;
  for( auto t : subset ) {
   v_Availability[ t ] = *( availability++ );
  }

  if( not_dry_run( issueAMod ) && constraints_generated() ) {
   // Change the abstract representation
   for( auto t : subset )
    update_availability_dependents( t , issueAMod );
  }
 }

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification
  if( ! ordered ) {
   std::sort( subset.begin(), subset.end() );
  }
  Block::add_Modification( std::make_shared< ThermalUnitBlockSbstMod >( this ,
                                               ThermalUnitBlockMod::eSetAv ,
                                               std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );
 }
}  // end( ThermalUnitBlock::set_availability )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_availability
( std::vector< double >::const_iterator values, Block::Range rng,
  c_ModParam issuePMod, c_ModParam issueAMod ) {

 rng.second = std::min( rng.second, f_time_horizon );
 if( rng.second <= rng.first ) {
  return;
 }

 if( v_Availability.empty() ) {
  if( std::all_of( values,
                   values + ( rng.second - rng.first ),
                   []( double cst ) {
                    return ( cst == 1.0 );
                   } ) ) {
   return;
  }

  v_Availability.assign( get_time_horizon() , 1.0 );
 }

 if( rng.first >= v_MaxPower.size() ) {
  throw( std::invalid_argument
         ( "ThermalUnitBlock::set_availability: invalid first endpoint of "
           "range: " + std::to_string( rng.first ) ) );
 }

 // If nothing changes, return
 if( std::equal( values,
                 values + ( rng.second - rng.first ),
                 v_Availability.begin() + rng.first ) ) {
  return;
 }

 // Check consistency

 auto availability = values;
 for( Index t = rng.first ; t < rng.second ; ++t ) {
  if( ! availability_is_consistent( t , *availability ) )
   throw( std::logic_error
          ( "ThermalUnitBlock::set_availability: availability (" +
            std::to_string( *availability ) + ") at time " +
            std::to_string( t ) + " is not consistent." ) );
  std::advance( availability , 1 );
 }

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  std::copy( values,
             values + ( rng.second - rng.first ),
             v_Availability.begin() + rng.first );

  if( not_dry_run( issueAMod ) && constraints_generated() ) {
   // Change the abstract representation
   for( Index t = rng.first ; t < rng.second ; ++t )
    update_availability_dependents( t , issueAMod );
  }
 }

 if( issue_pmod( issuePMod ) ) {
  Block::add_Modification( std::make_shared< ThermalUnitBlockRngdMod >( this ,
                                               ThermalUnitBlockMod::eSetAv ,
                                               rng ) ,
                           Observer::par2chnl( issuePMod ) );
 }
}  // end( ThermalUnitBlock::set_availability )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_maximum_power
( std::vector< double >::const_iterator values, Block::Subset && subset,
  const bool ordered, c_ModParam issuePMod, c_ModParam issueAMod ) {

 if( subset.empty() ) {
  return;
 }

 if( v_MaxPower.empty() ) {
  if( std::all_of( values,
                   values + subset.size(),
                   []( double cst ) {
                    return ( cst == 0 );
                   } ) ) {
   return;
  }

  v_MaxPower.assign( get_time_horizon() , 0 );
 }

 // If nothing changes, return
 bool identical = true;
 auto availability = values;
 for( auto t : subset ) {
  if( t >= v_MaxPower.size() ) {
   throw( std::invalid_argument
          ( "ThermalUnitBlock::set_maximum_power: invalid index in subset: "
            + std::to_string( t ) ) );
  }
  if( v_MaxPower[ t ] != *( availability++ ) ) {
   identical = false;
   break;
  }
 }
 if( identical )
  return;

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  availability = values;
  for( auto t : subset ) {
   v_MaxPower[ t ] = *( availability++ );
  }

  if( not_dry_run( issueAMod ) && constraints_generated() ) {
   // Change the abstract representation

   for( auto t : subset ) {
    auto f = dynamic_cast<LinearFunction *>
     ( MaxPower_Constraints[ t ].get_function() );
    auto var_index = f->is_active( & v_commitment[ t ] );
    assert( var_index < f->get_num_active_var() );
    f->modify_coefficient( var_index , get_operational_max_power( t ) ,
                           issueAMod );
   }
  }
 }

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification
  if( !ordered ) {
   std::sort( subset.begin(), subset.end() );
  }
  Block::add_Modification( std::make_shared< ThermalUnitBlockSbstMod >( this ,
                                               ThermalUnitBlockMod::eSetMaxP ,
                                               std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );
 }
}  // end( ThermalUnitBlock::set_maximum_power )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_maximum_power
( std::vector< double >::const_iterator values, Block::Range rng,
  c_ModParam issuePMod, c_ModParam issueAMod ) {

 rng.second = std::min( rng.second, f_time_horizon );
 if( rng.second <= rng.first ) {
  return;
 }

 if( v_MaxPower.empty() ) {
  if( std::all_of( values,
                   values + ( rng.second - rng.first ),
                   []( double cst ) {
                    return ( cst == 0 );
                   } ) ) {
   return;
  }

  v_MaxPower.assign( get_time_horizon() , 0 );
 }

 if( rng.first >= v_MaxPower.size() ) {
  throw( std::invalid_argument
         ( "ThermalUnitBlock::set_maximum_power: invalid first endpoint of "
           "range: " + std::to_string( rng.first ) ) );
 }

 // If nothing changes, return
 if( std::equal( values,
                 values + ( rng.second - rng.first ),
                 v_MaxPower.begin() + rng.first ) ) {
  return;
 }

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  std::copy( values,
             values + ( rng.second - rng.first ),
             v_MaxPower.begin() + rng.first );

  if( not_dry_run( issueAMod ) && constraints_generated() ) {
   // Change the abstract representation

   for( Index t = rng.first; t < rng.second; ++t ) {
    auto f = dynamic_cast<LinearFunction *>
     ( MaxPower_Constraints[ t ].get_function() );
    auto var_index = f->is_active( & v_commitment[ t ] );
    assert( var_index < f->get_num_active_var() );
    f->modify_coefficient( var_index , get_operational_max_power( t ) ,
                           issueAMod );
   }
  }
 }

 if( issue_pmod( issuePMod ) ) {
  Block::add_Modification( std::make_shared< ThermalUnitBlockRngdMod >( this ,
                                               ThermalUnitBlockMod::eSetMaxP ,
                                               rng ) ,
                           Observer::par2chnl( issuePMod ) );
 }
}  // end( ThermalUnitBlock::set_maximum_power )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_initial_power(
 std::vector< double >::const_iterator values,
 Block::Subset && subset,
 const bool ordered,
 c_ModParam issuePMod,
 c_ModParam issueAMod ) {

 if( subset.size() != 1 ) {
  return;
 }

 if( f_initial_power == *values ) {
  return;
 }

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation
  f_initial_power = *values;

  if( not_dry_run( issueAMod ) && constraints_generated() ) {
   auto initial_commitment = f_InitUpDownTime > 0 ? 1.0 : 0.0;
   RampUp_Constraints[ 0 ].set_rhs(
    v_DeltaRampUp[ 0 ] * initial_commitment + f_initial_power,
    issueAMod );
   RampDown_Constraints[ 0 ].set_lhs( f_initial_power, issueAMod );
  }
 }

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification
  Block::add_Modification(
   std::make_shared< ThermalUnitBlockSbstMod >( this,
                                                ThermalUnitBlockMod::eSetInitP,
                                                std::move( subset ) ),
   Observer::par2chnl( issuePMod ) );
 }
}  // end( ThermalUnitBlock::set_initial_power )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_initial_power(
 std::vector< double >::const_iterator values,
 Block::Range rng,
 c_ModParam issuePMod,
 c_ModParam issueAMod ) {

 rng.second = std::min( rng.second, f_time_horizon );
 if( rng.second != rng.first ) {
  return;
 }

 if( f_initial_power == *values ) {
  return;
 }

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation
  f_initial_power = *values;

  if( not_dry_run( issueAMod ) && constraints_generated() ) {
   auto initial_commitment = f_InitUpDownTime > 0 ? 1.0 : 0.0;
   RampUp_Constraints[ 0 ].set_rhs(
    v_DeltaRampUp[ 0 ] * initial_commitment + f_initial_power,
    issueAMod );
   RampDown_Constraints[ 0 ].set_lhs( f_initial_power, issueAMod );
  }
 }

 if( issue_pmod( issuePMod ) ) {
  Block::add_Modification(
   std::make_shared< ThermalUnitBlockRngdMod >( this,
                                                ThermalUnitBlockMod::eSetInitP,
                                                rng ),
   Observer::par2chnl( issuePMod ) );
 }
}  // end( ThermalUnitBlock::set_initial_power )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_init_updown_time(
 std::vector< int >::const_iterator values,
 Block::Subset && subset,
 const bool ordered,
 c_ModParam issuePMod,
 c_ModParam issueAMod ) {

 if( subset.size() != 1 ) {
  return;
 }

 if( f_InitUpDownTime == *values ) {
  return;
 }

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation
  f_InitUpDownTime = *values;

  if( not_dry_run( issueAMod ) && variables_generated() ) {
   // TODO Nuclear option
  }
 }

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification
  Block::add_Modification(
   std::make_shared< ThermalUnitBlockSbstMod >( this,
                                                ThermalUnitBlockMod::eSetInitUD,
                                                std::move( subset ) ),
   Observer::par2chnl( issuePMod ) );
 }
}  // end( ThermalUnitBlock::set_init_updown_time )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_init_updown_time(
 std::vector< int >::const_iterator values,
 Block::Range rng,
 c_ModParam issuePMod,
 c_ModParam issueAMod ) {

 rng.second = std::min( rng.second, f_time_horizon );
 if( rng.second != rng.first ) {
  return;
 }

 if( f_InitUpDownTime == *values ) {
  return;
 }

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation
  f_InitUpDownTime = *values;

  if( not_dry_run( issueAMod ) && variables_generated() ) {
   // TODO Nuclear option
  }
 }

 if( issue_pmod( issuePMod ) ) {
  Block::add_Modification(
   std::make_shared< ThermalUnitBlockRngdMod >( this,
                                                ThermalUnitBlockMod::eSetInitP,
                                                rng ),
   Observer::par2chnl( issuePMod ) );
 }
}  // end( ThermalUnitBlock::set_init_updown_time )

/*--------------------------------------------------------------------------*/

template< typename T >
void ThermalUnitBlock::decompress_vector( std::vector< T > & v ) {
 if( v.size() == 1 ) {
  v.resize( f_time_horizon, v[ 0 ] );
 } else if( v.size() < f_time_horizon ) {
  std::vector< T > temp = v;
  v.resize( f_time_horizon );
  Index j = 0;
  for( decltype( v_change_intervals )::size_type i = 0;
       i < v_change_intervals.size(); ++i ) {
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
} // end( ThermalUnitBlock::decompress_vector )

/*--------------------------------------------------------------------------*/
/*------------------- End File ThermalUnitBlock.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
