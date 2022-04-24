/*--------------------------------------------------------------------------*/
/*--------------------- File ThermalUnitBlock.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the ThermalUnitBlock class.
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
 * Copyright &copy by Antonio Frangioni, Ali Ghezelsoflu,
 *                    Rafael Durbano Lobato
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "LinearFunction.h"
#include "ThermalUnitBlock.h"

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

 auto clear_ZOconstraints =
  []( std::vector< ZOConstraint > & constraints ) {
   for( auto & constraint : constraints )
    constraint.clear();
  };
 clear_ZOconstraints( Commitment_bound_Constraints );
 clear_ZOconstraints( StartUp_Binary_bound_Constraints );
 clear_ZOconstraints( ShoutDown_Binary_bound_Constraints );

 auto clear_Boxconstraints =
         []( std::vector< BoxConstraint > & constraints ) {
          for( auto & constraint : constraints )
           constraint.clear();
         };
 clear_Boxconstraints(Commitment_fixed_to_One_Constraints);


 objective.clear();

}

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::deserialize( const netCDF::NcGroup & group )
{
#ifndef NDEBUG
 std::vector< std::string > expected_dims =
  { "TimeHorizon" , "NumberIntervals" };

 check_dimensions( group, expected_dims, std::cerr );

 std::vector< std::string > expected_vars = { "MinPower" , "MaxPower" ,
  "DeltaRampUp" , "DeltaRampDown" , "PrimaryRho" , "SecondaryRho" ,
  "LinearTerm" , "QuadTerm" , "ConstTerm" , "StartUpCost" , "FixedConsumption" ,
  "InertiaCommitment" , "InitialPower" , "MinUpTime" , "MinDownTime" ,
  "InitUpDownTime" , "Availability" };

 check_variables( group , expected_vars , std::cerr );
#endif

 UnitBlock::deserialize( group );

 // Mandatory variables

 ::deserialize( group , "MinPower" , v_MinPower , false );
 ::deserialize( group , "MaxPower" , v_MaxPower , false );

 // Optional variables

 if( ! ::deserialize( group , f_MinUpTime , "MinUpTime" ) )
  f_MinUpTime = 0;

 if( ! ::deserialize( group , f_MinDownTime , "MinDownTime" ) )
  f_MinDownTime = 0;

 if( ! ::deserialize( group, f_initial_power  , "InitialPower" ) )
  f_initial_power = 0;

 if( ! ::deserialize( group , f_InitUpDownTime , "InitUpDownTime" ) ) {
  if( f_initial_power == 0 )
   f_InitUpDownTime = - f_MinDownTime;
  else
   f_InitUpDownTime = f_MinUpTime;
  }

 if( ! ::deserialize( group , "Availability" , v_Availability ) )
  v_Availability.resize( get_time_horizon() , 1.0 );

 ::deserialize( group , "DeltaRampUp" , v_DeltaRampUp );
 ::deserialize( group , "DeltaRampDown" , v_DeltaRampDown );
 ::deserialize( group , "PrimaryRho" , v_PrimaryRho );
 ::deserialize( group , "SecondaryRho" , v_SecondaryRho );
 ::deserialize( group , "LinearTerm" , v_LinearTerm );
 ::deserialize( group , "QuadTerm" , v_QuadTerm );
 ::deserialize( group , "ConstTerm" , v_ConstTerm );
 ::deserialize( group , "StartUpCost" , v_StartUpCost );
 ::deserialize( group , "FixedConsumption" , v_fixed_consumption );
 ::deserialize( group , "InertiaCommitment" , v_inertia_commitment );

 // Decompress vectors
 decompress_vector( v_MinPower );
 decompress_vector( v_MaxPower );
 decompress_vector( v_Availability );
 decompress_vector( v_DeltaRampUp );
 decompress_vector( v_DeltaRampDown );
 decompress_vector( v_PrimaryRho );
 decompress_vector( v_SecondaryRho );
 decompress_vector( v_LinearTerm );
 decompress_vector( v_QuadTerm );
 decompress_vector( v_ConstTerm );
 decompress_vector( v_StartUpCost );
 decompress_vector( v_fixed_consumption );
 decompress_vector( v_inertia_commitment );

 check_data_consistency();

 }  // end( ThermalUnitBlock::deserialize )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::check_data_consistency() const {

 // Minimum and maximum power

 assert( v_MinPower.size() == f_time_horizon );
 assert( v_MaxPower.size() == f_time_horizon );

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {
  if( v_MinPower[ t ] > v_MaxPower[ t ] )
   throw( std::logic_error( "ThermalUnitBlock::check_data_consistency: "
                            "minimum power at time " + std::to_string( t ) +
                            " is " + std::to_string( v_MinPower[ t ] ) +
                            ", which is greater than the maximum power, which "
                            "is " + std::to_string( v_MaxPower[ t ] ) + "." ) );

  if( v_MinPower[ t ] < 0 )
   throw( std::logic_error( "ThermalUnitBlock::check_data_consistency: "
                            "minimum power for time step "
                            + std::to_string( t ) + " is " +
                            std::to_string( v_MinPower[ t ] ) +
                            ", but it must be nonnegative." ) );
 }

 // Availability

 if( ! v_Availability.empty() ) {
  assert( v_Availability.size() == f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   if( ( v_Availability[ t ] < 0 ) || ( v_Availability[ t ] > 1 ) )
    throw( std::logic_error( "ThermalUnitBlock::check_data_consistency: "
                             "availability for time step " +
                             std::to_string( t ) + " is " +
                             std::to_string( v_Availability[ t ] ) +
                             ", but it must be between 0 and 1." ) );
 }

 // Delta ramp-up

 if( ! v_DeltaRampUp.empty() ) {
  assert( v_DeltaRampUp.size() == f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   if( v_DeltaRampUp[ t ] < 0 )
    throw( std::logic_error( "ThermalUnitBlock::check_data_consistency: "
                             "delta ram pup for time step " +
                             std::to_string( t ) + " is " +
                             std::to_string( v_DeltaRampUp[ t ] ) +
                             ", but it must be nonnegative." ) );
 }

 // Delta ramp-down

 if( ! v_DeltaRampDown.empty() ) {
  assert( v_DeltaRampDown.size() == f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   if( v_DeltaRampDown[ t ] < 0 )
    throw( std::logic_error( "ThermalUnitBlock::check_data_consistency: "
                             "delta ramp down for time step " +
                             std::to_string( t ) + " is " +
                             std::to_string( v_DeltaRampDown[ t ] ) +
                             ", but it must be nonnegative" ) );
 }

 // Quadratic term of the objective function

 if( ! v_QuadTerm.empty() ) {
  assert( v_QuadTerm.size() == f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   if( v_QuadTerm[ t ] < 0 )
    throw( std::logic_error( "ThermalUnitBlock::check_data_consistency: "
                             "quadratic term for time " +
                             std::to_string( t ) + " is " +
                             std::to_string( v_QuadTerm[ t ] ) +
                             ", but it must be nonnegative." ) );
 }

 // MinUpTime

 if( f_MinUpTime < 0 )
  throw( std::logic_error( "ThermalUnitBlock::check_data_consistency: "
                           "minimum up time is "
                           + std::to_string( f_MinUpTime ) +
                           ", but it must be nonnegative." ) );

 // MinDownTime

 if( f_MinDownTime < 0 )
  throw( std::logic_error( "ThermalUnitBlock::check_data_consistency: "
                           "minimum down time is " +
                           std::to_string( f_MinDownTime ) +
                           ", but it must be nonnegative." ) );

 // InitialPower

 if( f_initial_power < 0 )
  throw( std::logic_error( "ThermalUnitBlock::check_data_consistency: "
                           "initial power is " +
                           std::to_string( f_initial_power ) +
                           ", but it must be nonnegative." ) );
}

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

 int relax_binary = 0;
 auto config = dynamic_cast<SimpleConfiguration<int> *>( stvv );
 if( ( ! config ) && f_BlockConfig &&
     f_BlockConfig->f_static_variables_Configuration )
  config = dynamic_cast< SimpleConfiguration< int > * >
   ( f_BlockConfig->f_static_variables_Configuration );
 if( config )
  relax_binary = config->f_value;

/*--------------------------------------------------------------------------*/

 // Commitment Variable

 v_commitment.resize( f_time_horizon );
 for( auto & var : v_commitment ) {
  if( relax_binary )
   var.set_type( ColVariable::kPosUnitary );
  else
   var.set_type( ColVariable::kBinary );
 }
 add_static_variable( v_commitment, "u_thermal" );

 // Active Power Variable

 v_active_power.resize( f_time_horizon );
 for( auto & var : v_active_power )
  var.set_type( ColVariable::kNonNegative );
 add_static_variable( v_active_power, "p_thermal" );

 // Primary Spinning Reserve Variable
 if( reserve_vars & 1u ) { // if UCBlock has primary demand variables
  if( !v_PrimaryRho.empty() ) { // if unit produces any primary reserve
   v_primary_spinning_reserve.resize( f_time_horizon );
   for( auto & var : v_primary_spinning_reserve ) {
    var.set_type( ColVariable::kNonNegative );
   }
   add_static_variable( v_primary_spinning_reserve, "pr_thermal" );
  }
 }

 // Secondary Spinning Reserve Variable
 if( reserve_vars & 2u ) { // if UCBlock has secondary demand variables
  if( !v_SecondaryRho.empty() ) { // if unit produces any secondary reserve
   v_secondary_spinning_reserve.resize( f_time_horizon );
   for( auto & var : v_secondary_spinning_reserve ) {
    var.set_type( ColVariable::kNonNegative );
   }
   add_static_variable( v_secondary_spinning_reserve, "sc_thermal" );
  }
 }
/*--------------------------------------------------------------------------*/
 // START UP AND SHUT DOWN BINARY VARIABLES
 auto startup_shutdown_size = f_time_horizon - init_t;

 if( startup_shutdown_size > 0 ) {

  v_start_up.resize( startup_shutdown_size );
  for( auto & var : v_start_up ) {
   if( relax_binary )
    var.set_type( ColVariable::kPosUnitary );
   else
    var.set_type( ColVariable::kBinary );
  }
  add_static_variable( v_start_up, "v" );

  v_shut_down.resize( startup_shutdown_size );
  for( auto & var : v_shut_down ) {
   if( relax_binary )
    var.set_type( ColVariable::kPosUnitary );
   else
    var.set_type( ColVariable::kBinary );
  }
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
     v_active_power[ t ].is_fixed( true , eNoMod );
    }
   }

   if( ! v_primary_spinning_reserve.empty() ) {
    for( Index t = 0 ; t < init_t ; ++t ) {
     v_primary_spinning_reserve[ t ].set_value( 0.0 );
     v_primary_spinning_reserve[ t ].is_fixed( true , eNoMod );
    }
   }

   if( ! v_secondary_spinning_reserve.empty() ) {
    for( Index t = 0 ; t < init_t ; ++t ) {
     v_secondary_spinning_reserve[ t ].set_value( 0.0 );
     v_secondary_spinning_reserve[ t ].is_fixed( true , eNoMod );
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
     v_commitment[ t ].is_fixed( true , eNoMod );
    }
   }
  }

  if( f_InitUpDownTime > 0 ) {

   if( f_MinDownTime < startup_shutdown_size ) {
    for( Index t = 0 ; t < f_MinDownTime ; ++t ) {
     v_start_up[ t ].set_value( 0.0 );
     v_start_up[ t ].is_fixed( true , eNoMod );
    }
   }
   else {
    for( Index t = 0 ; t < startup_shutdown_size ; ++t ) {
     v_start_up[ t ].set_value( 0.0 );
     v_start_up[ t ].is_fixed( true , eNoMod );
    }
   }
  }
  else {
   if( f_MinUpTime < startup_shutdown_size ) {
    for( Index t = 0 ; t < f_MinUpTime ; ++t ) {
     v_shut_down[ t ].set_value( 0.0 );
     v_shut_down[ t ].is_fixed( true , eNoMod );
    }
   }
   else {
    for( Index t = 0 ; t < startup_shutdown_size ; ++t ) {
     v_shut_down[ t ].set_value( 0.0 );
     v_shut_down[ t ].is_fixed( true , eNoMod );
    }
   }
  }
 }

 else if ( init_t == 0 ) {

  if( f_InitUpDownTime > 0 ) {

   if( f_MinDownTime < f_time_horizon ) {
    for( Index t = 0 ; t < f_MinDownTime ; ++t ) {
     v_start_up[ t ].set_value( 0.0 );
     v_start_up[ t ].is_fixed( true , eNoMod );
    }
   }
   else {
    for( Index t = 0 ; t < f_time_horizon ; ++t ) {
     v_start_up[ t ].set_value( 0.0 );
     v_start_up[ t ].is_fixed( true , eNoMod );
    }
   }
  }
  if( f_InitUpDownTime <= 0 ) {
   if( f_MinUpTime < f_time_horizon ) {
    for( Index t = 0 ; t < f_MinUpTime ; ++t ) {
     v_shut_down[ t ].set_value( 0.0 );
     v_shut_down[ t ].is_fixed( true , eNoMod );
    }
   }
   else {
    for( Index t = 0 ; t < f_time_horizon ; ++t ) {
     v_shut_down[ t ].set_value( 0.0 );
     v_shut_down[ t ].is_fixed( true , eNoMod );
    }
   }
  }
 }

 set_variables_generated();
} // end( ThermalUnitBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::generate_abstract_constraints( Configuration * stcc ) {

 if( constraints_generated())
  return; // constraints have already been generated

 int generate_ZOConstraint = 0;
 auto config = dynamic_cast<SimpleConfiguration< int > *>( stcc );
 if(( !config ) && f_BlockConfig &&
    f_BlockConfig->f_static_constraints_Configuration )
  config = dynamic_cast< SimpleConfiguration< int > * >
  ( f_BlockConfig->f_static_constraints_Configuration );
 if( config )
  generate_ZOConstraint = config->f_value;

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

   for( Index t = 1; t < f_time_horizon; ++t ) {

    auto lf = new LinearFunction();

    lf->add_variable( &v_commitment[t], 1.0 );
    lf->add_variable( &v_start_up[t], -1.0 );
    lf->add_variable( &v_shut_down[t], 1.0 );
    lf->add_variable( &v_commitment[t - 1], -1.0 );
    StartUp_ShutDown_Variables_Constraints[t].set_both( 0.0 );
    StartUp_ShutDown_Variables_Constraints[t].set_function( lf );
   }
  }
 }

 if( init_t > 0 ) {

  auto startup_shutdown_const_size = static_cast<int>( f_time_horizon - init_t );

  if( startup_shutdown_const_size > 0 ) {

   StartUp_ShutDown_Variables_Constraints.resize( startup_shutdown_const_size );

   // Initial condition

   auto l_function = new LinearFunction();

   l_function->add_variable( &v_commitment[init_t], 1.0 );
   l_function->add_variable( &v_start_up[0], -1.0 );
   l_function->add_variable( &v_shut_down[0], 1.0 );

   if( f_InitUpDownTime <= 0 ) { // -f_InitUpDownTime < f_MinDownTime
    StartUp_ShutDown_Variables_Constraints[0].set_both( 0.0 );
    StartUp_ShutDown_Variables_Constraints[0].set_function( l_function );
   } else { // f_InitUpDownTime > 0 && f_InitUpDownTime < f_MinUpTime
    StartUp_ShutDown_Variables_Constraints[0].set_both( 1.0 );
    StartUp_ShutDown_Variables_Constraints[0].set_function( l_function );
   }

   for( Index t = init_t + 1, constraint_index = 1; t < f_time_horizon;
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
                        "StartUp_ShutDown_Variables_Constraints" );

 // Initializing turn on constraints (start up constraints)

 auto startup_const_size = static_cast<int>
 ( f_time_horizon - init_t - f_MinUpTime );

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

  add_static_constraint( StartUp_Constraints, "StartUp_Commitment_Constraints" );
 }

 // Initializing turn off constraints (shut down constraints)

 auto shutdown_const_size = static_cast<int>
 ( f_time_horizon - init_t - f_MinDownTime );

 if( shutdown_const_size > 0 ) {

  ShutDown_Constraints.resize( shutdown_const_size );

  for( Index t = init_t + f_MinDownTime, constraint_index = 0;
       t < f_time_horizon; ++t, ++constraint_index ) {

   auto linear_function = new LinearFunction();

   for( Index s = t - f_MinDownTime; s < t; ++s ) {
    linear_function->add_variable( &v_shut_down[s - init_t + 1], 1.0 );
   }

   linear_function->add_variable( &v_commitment[t], 1.0 );
   ShutDown_Constraints[constraint_index].set_lhs( -Inf< double >());
   ShutDown_Constraints[constraint_index].set_rhs( 1.0 );
   ShutDown_Constraints[constraint_index].set_function( linear_function );
  }

  add_static_constraint( ShutDown_Constraints, "ShutDown_Commitment_Constraints" );
 }

/*--------------------------------------------------------------------------*/
 // Initializing ramp-up constraints with 3-Binary Variables

 if( !v_DeltaRampUp.empty() && !v_DeltaRampDown.empty()) {

  if( f_InitUpDownTime > 0 ) {
   for( Index t = 0; t < f_time_horizon; ++t ) {
    if( f_initial_power + v_DeltaRampUp[0] < get_operational_min_power( 0 ) ||
        f_initial_power - v_DeltaRampDown[0] > get_operational_max_power( 0 ) ) {
     throw( std::logic_error
            ( "ThermalUnitBlock::Ramp Constraints: when f_InitUpDownTime > 0,"
              " it must be that f_initial_power + v_DeltaRampUp[ 0 ] >= "
              "get_operational_min_power( 0 ) and f_initial_power - "
              "v_DeltaRampDown[ 0 ] <= get_operational_max_power( 0 ) " ) );
    }
   }
  }
 }

 if( !v_DeltaRampUp.empty()) {

  RampUp_Constraints.resize( f_time_horizon );

  // Initial condition

  if( init_t == 0 ) {

   // Initial condition

   auto linear_function = new LinearFunction();
   linear_function->add_variable( &v_active_power[0], 1.0 );
   linear_function->add_variable( &v_start_up[0],
                                  -get_operational_min_power( 0 ));
   auto initial_commitment = ( f_InitUpDownTime > 0 ? 1.0 : 0.0 );

   RampUp_Constraints[0].set_lhs( -Inf< double >());

   if( f_InitUpDownTime > 0 ) {
    RampUp_Constraints[0].set_rhs
            (( v_DeltaRampUp[0] * initial_commitment ) + f_initial_power );
   } else {
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

   RampUp_Constraints[0].set_lhs( -Inf< double >());

   if( f_InitUpDownTime > 0 ) {
    RampUp_Constraints[0].set_rhs
            (( v_DeltaRampUp[0] * initial_commitment ) + f_initial_power );
   } else {
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

   if( f_InitUpDownTime <= 0 ) { // -f_InitUpDownTime < f_MinDownTime

    auto LFunction = new LinearFunction();
    LFunction->add_variable( &v_active_power[init_t], 1.0 );
    LFunction->add_variable( &v_start_up[0],
                             -get_operational_min_power( init_t ));

    RampUp_Constraints[init_t].set_lhs( -Inf< double >());
    RampUp_Constraints[init_t].set_rhs( 0.0 );
    RampUp_Constraints[init_t].set_function( LFunction );

   } else { // f_InitUpDownTime > 0 && f_InitUpDownTime < f_MinUpTime

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
 if( !v_DeltaRampDown.empty()) {

  RampDown_Constraints.resize( f_time_horizon );

  // Initial condition
  if( init_t == 0 ) {
   auto linear_function = new LinearFunction();
   linear_function->add_variable( &v_active_power[0], 1.0 );
   linear_function->add_variable( &v_commitment[0], v_DeltaRampDown[0] );
   linear_function->add_variable( &v_shut_down[0],
                                  get_operational_min_power( 0 ));
   if( f_InitUpDownTime > 0 ) {
    RampDown_Constraints[0].set_lhs( f_initial_power );
   } else {
    RampDown_Constraints[0].set_lhs( 0.0 );
   }
   RampDown_Constraints[0].set_rhs( Inf< double >());
   RampDown_Constraints[0].set_function( linear_function );

   // Remaining constraints
   for( Index t = 1; t < f_time_horizon; ++t ) {

    auto lf = new LinearFunction();

    lf->add_variable( &v_active_power[t], 1.0 );
    lf->add_variable( &v_active_power[t - 1], -1.0 );
    lf->add_variable( &v_shut_down[t], get_operational_min_power( t ));
    lf->add_variable( &v_commitment[t], v_DeltaRampDown[t] );

    RampDown_Constraints[t].set_lhs( 0.0 );
    RampDown_Constraints[t].set_rhs( Inf< double >());
    RampDown_Constraints[t].set_function( lf );
   }
  }

  if( init_t > 0 ) {
   // Initial condition
   auto linear_function = new LinearFunction();
   linear_function->add_variable( &v_active_power[0], 1.0 );
   linear_function->add_variable( &v_commitment[0], v_DeltaRampDown[0] );
   if( f_InitUpDownTime > 0 ) {
    RampDown_Constraints[0].set_lhs( f_initial_power );
   } else {
    RampDown_Constraints[0].set_lhs( 0.0 );
   }
   RampDown_Constraints[0].set_rhs( Inf< double >());
   RampDown_Constraints[0].set_function( linear_function );

   for( Index t = 1; t < init_t; ++t ) {
    auto lf = new LinearFunction();

    lf->add_variable( &v_active_power[t - 1], -1.0 );
    lf->add_variable( &v_active_power[t], 1.0 );
    lf->add_variable( &v_commitment[t], v_DeltaRampDown[t] );
    RampDown_Constraints[t].set_lhs( 0.0 );
    RampDown_Constraints[t].set_rhs( Inf< double >());
    RampDown_Constraints[t].set_function( lf );
   }

   // Remaining constraints

   for( Index t = init_t; t < f_time_horizon; ++t ) {

    auto lf = new LinearFunction();

    lf->add_variable( &v_active_power[t], 1.0 );
    lf->add_variable( &v_active_power[t - 1], -1.0 );
    lf->add_variable( &v_shut_down[t - init_t],
                      get_operational_min_power( t ));
    lf->add_variable( &v_commitment[t], v_DeltaRampDown[t] );

    RampDown_Constraints[t].set_lhs( 0.0 );
    RampDown_Constraints[t].set_rhs( Inf< double >());
    RampDown_Constraints[t].set_function( lf );
   }
  }

  add_static_constraint( RampDown_Constraints, "RampDown_Constraints_Thermal" );

 }

/*--------------------------------------------------------------------------*/

 // Initializing minimum power constraints

 MinPower_Constraints.resize( f_time_horizon );

 for( Index t = 0; t < f_time_horizon; ++t ) {

  auto linear_function = new LinearFunction();

  linear_function->add_variable( &v_active_power[t], 1.0 );
  if( reserve_vars & 1u ) { // if UCBlock has primary demand variables
   if( !v_PrimaryRho.empty()) { // if unit produces any primary reserve
    linear_function->add_variable( &v_primary_spinning_reserve[t], -1.0 );
   }
  }
  if( reserve_vars & 2u ) { // if UCBlock has secondary demand variables
   if( !v_SecondaryRho.empty()) { // if unit produces any secondary reserve
    linear_function->add_variable( &v_secondary_spinning_reserve[t], -1.0 );
   }
  }
  linear_function->add_variable( &v_commitment[t],
                                 -get_operational_min_power( t ));

  MinPower_Constraints[t].set_rhs( Inf< double >());
  MinPower_Constraints[t].set_lhs( 0.0 );
  MinPower_Constraints[t].set_function( linear_function );
 }
 add_static_constraint( MinPower_Constraints, "MinPower_Constraints_Thermal" );

 // Initializing maximum power constraints
 MaxPower_Constraints.resize( f_time_horizon );

 for( Index t = 0; t < f_time_horizon; ++t ) {

  auto linear_function = new LinearFunction();

  linear_function->add_variable( &v_active_power[t], -1.0 );
  if( reserve_vars & 1u ) { // if UCBlock has primary demand variables
   if( !v_PrimaryRho.empty()) { // if unit produces any primary reserve
    linear_function->add_variable( &v_primary_spinning_reserve[t], -1.0 );
   }
  }
  if( reserve_vars & 2u ) { // if UCBlock has secondary demand variables
   if( !v_SecondaryRho.empty()) { // if unit produces any secondary reserve
    linear_function->add_variable( &v_secondary_spinning_reserve[t], -1.0 );
   }
  }
  linear_function->add_variable( &v_commitment[t],
                                 get_operational_max_power( t ));

  MaxPower_Constraints[t].set_lhs( 0.0 );
  MaxPower_Constraints[t].set_rhs( Inf< double >());
  MaxPower_Constraints[t].set_function( linear_function );
 }
 add_static_constraint( MaxPower_Constraints, "MaxPower_Constraints_Thermal" );

 if( reserve_vars & 1u ) { // if UCBlock has primary demand variables
  if( !v_PrimaryRho.empty()) { // if unit produces any primary reserve

   // Initializing primary rho fraction constraints

   PrimaryRho_Constraints.resize( f_time_horizon );

   for( Index t = 0; t < f_time_horizon; ++t ) {

    auto linear_function = new LinearFunction();

    if( !v_PrimaryRho.empty()) {

     linear_function->add_variable( &v_active_power[t], v_PrimaryRho[t] );
    } else {
     linear_function->add_variable( &v_active_power[t], 0.0 );

    }
    linear_function->add_variable( &v_primary_spinning_reserve[t], -1.0 );

    PrimaryRho_Constraints[t].set_lhs( 0.0 );
    PrimaryRho_Constraints[t].set_rhs( Inf< double >());
    PrimaryRho_Constraints[t].set_function( linear_function );
   }

   add_static_constraint( PrimaryRho_Constraints, "PrimaryRho_Constraints_Thermal" );
  }
 }

 if( reserve_vars & 2u ) { // if UCBlock has secondary demand variables
  if( !v_SecondaryRho.empty()) { // if unit produces any secondary reserve
  // Initializing secondary rho fraction constraints

  SecondaryRho_Constraints.resize( f_time_horizon );

  for( Index t = 0; t < f_time_horizon; ++t ) {

   auto linear_function = new LinearFunction();
   if( !v_SecondaryRho.empty()) {
    linear_function->add_variable( &v_active_power[t], v_SecondaryRho[t] );
   } else {
    linear_function->add_variable( &v_active_power[t], 0.0 );
   }
   linear_function->add_variable( &v_secondary_spinning_reserve[t], -1.0 );

   SecondaryRho_Constraints[t].set_lhs( 0.0 );
   SecondaryRho_Constraints[t].set_rhs( Inf< double >());
   SecondaryRho_Constraints[t].set_function( linear_function );
  }

  add_static_constraint( SecondaryRho_Constraints, "SecondaryRho_Constraints_Thermal" );
 }
}
/*-------------------------------ZOConstraint-------------------------------*/

  if( generate_ZOConstraint ) {

   // the commitment bound constraints
   Commitment_bound_Constraints.resize( f_time_horizon );
   for( Index t = 0 ; t < f_time_horizon ; ++t ) {
    Commitment_bound_Constraints[ t ].set_variable(&v_commitment[t]);
   }
   add_static_constraint( Commitment_bound_Constraints,
                          "Commitment_bound_Thermal" );

   // the startup binary bound constraints
   auto startup_shutdown_size = f_time_horizon - init_t;

   StartUp_Binary_bound_Constraints.resize( startup_shutdown_size );
   for( Index t = 0 ; t < startup_shutdown_size ; ++t ) {
    StartUp_Binary_bound_Constraints[ t ].set_variable(&v_start_up[t]);
   }
   add_static_constraint( StartUp_Binary_bound_Constraints,
                          "StartUp_binary_bound_Thermal" );
   // the shut down binary bound constraints

   ShoutDown_Binary_bound_Constraints.resize( startup_shutdown_size );
   for( Index t = 0 ; t < startup_shutdown_size ; ++t ) {
    ShoutDown_Binary_bound_Constraints[ t ].set_variable(&v_shut_down[t]);
   }
   add_static_constraint( ShoutDown_Binary_bound_Constraints,
                          "ShoutDown_binary_bound_Thermal" );
  }
/*-------------------------------BoxConstraint-------------------------------*/
 if( init_t > 0 && f_InitUpDownTime > 0 && f_InitUpDownTime < f_MinUpTime ) {

  // the commitment fixed to one BoxConstraints
  Commitment_fixed_to_One_Constraints.resize( f_time_horizon );
  for( Index t = 0 ; t < init_t ; ++t ) {
   Commitment_fixed_to_One_Constraints[ t ].set_lhs( 1);
   Commitment_fixed_to_One_Constraints[ t ].set_rhs( 1);
   Commitment_fixed_to_One_Constraints[ t ].set_variable(&v_commitment[t]);
  }
  add_static_constraint( Commitment_fixed_to_One_Constraints,
                         "Commitment_fixed_to_one_Thermal" );
 }


 set_constraints_generated();

} // end( ThermalUnitBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

/// verifies whether the current solution is feasible for the given constraints
/** This function checks whether the relative violation of each RowConstraint
 * in the given group of RowConstraint is not greater than the provided
 * tolerance.
 *
 * @return This function returns true if and only if the relative violation of
 *         each RowConstraint in the given group is not greater than the given
 *         tolerance. */

template<class C>
static std::enable_if_t< std::is_base_of_v< RowConstraint , C > , bool >
is_feasible( std::vector< C > & constraints , double tolerance ) {
 for( auto & constraint : constraints ) {
  if( constraint.is_relaxed() )
   continue;
  constraint.compute();
  if( constraint.rel_viol() > tolerance )
   return false;
 }
 return true;
}

/*--------------------------------------------------------------------------*/

/// verifies whether the given ColVariable are feasible
/** This function returns true if and only if each given ColVariable is
 * feasible with respect to the given tolerance (see
 * ColVariable::is_feasible()).
 *
 * @return This function returns true if and only if each of the given
 *         ColVariable is feasible considering the given tolerance. */

template<class V>
static std::enable_if_t< std::is_base_of_v< ColVariable , V > , bool >
is_feasible( const std::vector< V > & variables , double tolerance ) {
 for( const auto & variable : variables ) {
  if( ! variable.is_feasible( tolerance ) )
   return false;
 }
 return true;
}

/*--------------------------------------------------------------------------*/

bool ThermalUnitBlock::is_feasible( bool useabstract , Configuration * fsbc ) {

 // Retrieve the tolerance.

 auto config = dynamic_cast< SimpleConfiguration< double > * >( fsbc );

 if( ( ! config ) && f_BlockConfig )
  config = dynamic_cast< SimpleConfiguration< double > * >
   ( f_BlockConfig->f_is_feasible_Configuration );

 // If a tolerance has not been provided, use the default tolerance.
 const auto tolerance = config ? config->f_value : 1.0e-8;

 // Notice that the ZOConstraint are not checked, since the corresponding
 // check is made on the ColVariable.

 return
  UnitBlock::is_feasible( useabstract )
  // Constraints
  && ::is_feasible( Power_StartUp_ShutDown_Variables_Constraints , tolerance )
  && ::is_feasible( Power_StartUp_Variable_Constraints , tolerance )
  && ::is_feasible( Power_ShutDown_Variable_Constraints , tolerance )
  && ::is_feasible( StartUp_ShutDown_Variables_Constraints , tolerance )
  && ::is_feasible( StartUp_Constraints , tolerance )
  && ::is_feasible( ShutDown_Constraints , tolerance )
  && ::is_feasible( RampUp_Constraints , tolerance )
  && ::is_feasible( RampDown_Constraints , tolerance )
  && ::is_feasible( PrimaryRho_Constraints , tolerance )
  && ::is_feasible( SecondaryRho_Constraints , tolerance )
  && ::is_feasible( MinPower_Constraints , tolerance )
  && ::is_feasible( MaxPower_Constraints , tolerance )
  && ::is_feasible( Commitment_fixed_to_One_Constraints , tolerance )
  // Variables
  && ::is_feasible( v_start_up , tolerance )
  && ::is_feasible( v_shut_down , tolerance )
  && ::is_feasible( v_commitment , tolerance )
  && ::is_feasible( v_active_power , tolerance )
  && ::is_feasible( v_primary_spinning_reserve , tolerance )
  && ::is_feasible( v_secondary_spinning_reserve , tolerance );

} // end( ThermalUnitBlock::is_feasible )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::generate_objective( Configuration * objc )
{
 if( objective_generated() )
  return;  // Objective has already been generated

 // initialize Objective
 //
 // the order of the variables in the DQuadFunction is:
 //
 // - first f_time_horizon - init_t start-up variables
 //
 // - then f_time_horizon active power variables (which may have the
 //   nonzero quadratic cost coefficient, while the others do not)
 //
 // - then f_time_horizon commitment variables
 //
 // - then possibly f_time_horizon primary reserve variables
 //
 // - then possibly f_time_horizon secondary reserve variables
 //
 // this arrangement is exploited in add_Modification to easily map
 // indices in the coefficients of the DQuadFunction back into indices
 // of the original variables (and figure out the kind of variable)

 if( v_commitment.size() != f_time_horizon )
  throw( std::logic_error(
            "ThermalUnitBlock::generate_objective: v_commitment must have "
            "size equal to the time horizon." ) );

 if( v_active_power.size() != f_time_horizon )
  throw( std::logic_error(
            "ThermalUnitBlock::generate_objective: v_active_power must have "
            "size equal to the time horizon." ) );

 if( v_start_up.size() != f_time_horizon - init_t )
  throw( std::logic_error(
            "ThermalUnitBlock::generate_objective: v_start_up must have "
            "size equal to the time horizon - init_t." ) );

 auto dquad_function = new DQuadFunction();

 for( Index t = init_t ; t < f_time_horizon ; ++t )
  dquad_function->add_variable( & v_start_up[ t - init_t ] ,
                                get_start_up_cost( t ) , 0.0 );

 for( Index t = 0 ; t < f_time_horizon ; ++t )
  dquad_function->add_variable( & v_active_power[ t ] ,
                                get_linear_term( t ) , get_quad_term( t ) );

 for( Index t = 0 ; t < f_time_horizon ; ++t )
  dquad_function->add_variable( & v_commitment[ t ] ,
                                get_const_term( t ) , 0.0 );

 // possibly add the primary and secondary spinning reserve variables

 bool add_primary_reserve = false;
 bool add_secondary_reserve = false;
 auto config = dynamic_cast<SimpleConfiguration<int> *>( objc );
 if( ( ! config ) && f_BlockConfig )
  config = dynamic_cast< SimpleConfiguration< int > * >(
                                 f_BlockConfig->f_objective_Configuration );
 if( config ) {
  add_primary_reserve = config->f_value & 1u;
  add_secondary_reserve = config->f_value & 2u;
  }

 if( ( ! v_primary_spinning_reserve.empty() ) && add_primary_reserve ) {
  // Add the primary spinning reserve variables

  if( v_primary_spinning_reserve.size() != f_time_horizon )
   throw( std::logic_error( "ThermalUnitBlock::generate_objective: v_primary_"
                            "spinning_reserve must have size equal to the "
                            "time horizon." ) );

  if( v_primary_spinning_reserve_cost.empty() )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    dquad_function->add_variable( & v_primary_spinning_reserve[ t ] , 0 , 0 );
  else
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    dquad_function->add_variable( & v_primary_spinning_reserve[ t ] ,
                                  v_primary_spinning_reserve_cost[ t ] , 0 );
  }

 if( ( ! v_secondary_spinning_reserve.empty() ) && add_secondary_reserve ) {
  // Add the secondary spinning reserve variables

  if( v_secondary_spinning_reserve.size() != f_time_horizon )
   throw( std::logic_error( "ThermalUnitBlock::generate_objective: v_secondary"
                            "_spinning_reserve must have size equal to the "
                            "time horizon." ) );

  if( v_secondary_spinning_reserve_cost.empty() )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    dquad_function->add_variable( & v_secondary_spinning_reserve[ t ] , 0 , 0 );
  else
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    dquad_function->add_variable( & v_secondary_spinning_reserve[ t ] ,
                                  v_secondary_spinning_reserve_cost[ t ] , 0 );
  }

 if( f_scale != 1.0 )
  // Update the Objective to take into account the scale factor
  update_objective( Range( 0 , Inf<Index>() ) , eNoMod );

 objective.set_function( dquad_function );
 objective.set_sense( Objective::eMin );

 // set Block objective
 this->set_objective( & objective );

 set_objective_generated();

 }  // end( ThermalUnitBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*-------- METHODS FOR LOADING, PRINTING & SAVING THE ThermalUnitBlock -----*/
/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::serialize( netCDF::NcGroup & group ) const {

 UnitBlock::serialize( group );

 // Serialize scalar variables.

 ::serialize( group , "InitialPower" , netCDF::NcDouble() , f_initial_power );
 ::serialize( group , "MinUpTime" , netCDF::NcUint() , f_MinUpTime );
 ::serialize( group , "MinDownTime" , netCDF::NcUint() , f_MinDownTime );
 ::serialize( group , "InitUpDownTime" , netCDF::NcInt() , f_InitUpDownTime );

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
          ( "ThermalUnitBlock::serialize: invalid dimension for variable " +
            var_name + ": " + std::to_string( data.size() ) + ". Its dimension "
            "must be one of the following: TimeHorizon, NumberIntervals, 1.") );
  }

  ::serialize( group , var_name , ncType , dimension , data ,
               allow_scalar_var );
 };

 serialize( "MinPower" , v_MinPower );
 serialize( "MaxPower" , v_MaxPower );
 serialize( "Availability" , v_Availability );
 serialize( "DeltaRampUp" , v_DeltaRampUp );
 serialize( "DeltaRampDown" , v_DeltaRampDown );
 serialize( "PrimaryRho" , v_PrimaryRho );
 serialize( "SecondaryRho" , v_SecondaryRho );
 serialize( "QuadTerm" , v_QuadTerm );
 serialize( "LinearTerm" , v_LinearTerm );
 serialize( "ConstTerm" , v_ConstTerm );
 serialize( "StartUpCost" , v_StartUpCost );
 serialize( "FixedConsumption" , v_fixed_consumption );
 serialize( "InertiaCommitment" , v_inertia_commitment );
}

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::add_Modification( sp_Mod mod, ChnlName chnl )
{
 if( mod->concerns_Block() ) {
  mod->concerns_Block( false );
  guts_of_add_Modification( mod.get() , chnl );
  }

 Block::add_Modification( mod, chnl );
 }

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::update_availability_dependents( Index t ,
                                                       c_ModParam issueAMod )
{
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

void ThermalUnitBlock::update_initial_power_in_constraints
( c_ModParam issueAMod ) {

 if( ! ( RampUp_Constraints.empty() || v_DeltaRampUp.empty() ) ) {
  if( f_InitUpDownTime > 0 )
   RampUp_Constraints[ 0 ].set_rhs( v_DeltaRampUp[ 0 ] + f_initial_power ,
                                    issueAMod );
 }

 if( ! RampDown_Constraints.empty() )
  if( f_InitUpDownTime > 0 )
   RampDown_Constraints[ 0 ].set_lhs( f_initial_power , issueAMod );

}  // end( ThermalUnitBlock::update_initial_power_in_constraints )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_initial_power
( std::vector< double >::const_iterator values , Block::Subset && subset ,
  const bool ordered , c_ModParam issuePMod , c_ModParam issueAMod ) {

 if( subset.empty() )
  return;

 // Find the last index 0
 auto index_it = std::find( subset.rbegin() , subset.rend() , 0 );

 if( index_it == subset.rend() )
  return; // 0 is not in subset; return

 std::advance( values , std::distance( index_it , subset.rend() ) - 1 );

 if( f_initial_power == *values )
  return; // nothing changes; return

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation
  f_initial_power = *values;

  if( not_dry_run( issueAMod ) && constraints_generated() ) {
   // Change the abstract representation
   update_initial_power_in_constraints( issueAMod );
  }
 }

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockMod >
                           ( this , ThermalUnitBlockMod::eSetInitP ) ,
                           Observer::par2chnl( issuePMod ) );
 }
}  // end( ThermalUnitBlock::set_initial_power )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_initial_power
( std::vector< double >::const_iterator values , Block::Range rng ,
  c_ModParam issuePMod , c_ModParam issueAMod ) {

 rng.second = std::min( rng.second , decltype( rng.second )( 1 ) );
 if( ! ( rng.first <= 0 && 0 < rng.second ) )
  return; // 0 does not belong to the range; return

 std::advance( values , - rng.first );

 if( f_initial_power == *values )
  return; // nothing changes; return

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation
  f_initial_power = *values;

  if( not_dry_run( issueAMod ) && constraints_generated() ) {
   // Change the abstract representation
   update_initial_power_in_constraints( issueAMod );
  }
 }

 if( issue_pmod( issuePMod ) ) {
  Block::add_Modification( std::make_shared< ThermalUnitBlockMod >
                           ( this , ThermalUnitBlockMod::eSetInitP ) ,
                           Observer::par2chnl( issuePMod ) );
 }
}  // end( ThermalUnitBlock::set_initial_power )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_startup_costs(
 std::vector< double >::const_iterator values,
 Subset && subset,
 const bool ordered,
 c_ModParam issuePMod,
 c_ModParam issueAMod ) {

 if( subset.empty() ) {
  return;
 }

 if( v_StartUpCost.empty() ) {
  if( std::all_of( values, values + subset.size(),
                   []( double cst ) {
                    return ( cst == 0 );
                   } ) ) {
   return;
  }

  v_StartUpCost.assign( get_time_horizon() , 0 );
 }

 // If nothing changes, return
 bool identical = true;
 auto startup_cost = values;

 for( auto t : subset ) {
  if( t >= v_StartUpCost.size() ) {
   throw std::invalid_argument
    ( "ThermalUnitBlock::set_startup_costs: invalid index in subset: "
      + std::to_string( t ) );
  }

  if( v_StartUpCost[ t ] != *( startup_cost++ ) ) {
   identical = false;
   break;
  }
 }

 if( identical ) {
  return;
 }

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  startup_cost = values;
  for( auto t : subset ) {
   v_StartUpCost[ t ] = *( startup_cost++ );
  }

  if( not_dry_run( issueAMod ) && objective_generated() ) {
   // Change the abstract representation

   auto qf = dynamic_cast<DQuadFunction *>(objective.get_function());

   for( auto t : subset ) {

    auto var_index = qf->is_active( &v_start_up[ t ] );
    assert( var_index < qf->get_num_active_var() );
    qf->modify_linear_coefficient( var_index,
                                   get_start_up_cost( t ),
                                   issueAMod );
   }
  }
 }

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification
  if( !ordered ) {
   std::sort( subset.begin(), subset.end() );
  }

  Block::add_Modification(
   std::make_shared< ThermalUnitBlockSbstMod >(
    this,
    ThermalUnitBlockMod::eSetSUC,
    std::move( subset ) ),
   Observer::par2chnl( issuePMod ) );
 }
}

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_startup_costs(
 std::vector< double >::const_iterator values,
 Range rng,
 c_ModParam issuePMod,
 c_ModParam issueAMod ) {

 rng.second = std::min( rng.second, f_time_horizon );
 if( rng.second <= rng.first ) {
  return;
 }

 if( v_StartUpCost.empty() ) {
  if( std::all_of( values,
                   values + ( rng.second - rng.first ),
                   []( double cst ) {
                    return ( cst == 0 );
                   } ) ) {
   return;
  }

  v_StartUpCost.assign( get_time_horizon() , 0 );
 }

 if( rng.first >= v_StartUpCost.size() ) {
  throw std::invalid_argument
   ( "ThermalUnitBlock::set_startup_costs: invalid first endpoint of "
     "range: " + std::to_string( rng.first ) );
 }

 // If nothing changes, return
 if( std::equal( values,
                 values + ( rng.second - rng.first ),
                 v_StartUpCost.begin() + rng.first ) ) {
  return;
 }

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  std::copy( values,
             values + ( rng.second - rng.first ),
             v_StartUpCost.begin() + rng.first );

  if( not_dry_run( issueAMod ) && objective_generated() ) {
   // Change the abstract representation

   auto qf = dynamic_cast<DQuadFunction *>(objective.get_function());

   for( Index t = rng.first; t < rng.second; ++t ) {

    auto var_index = qf->is_active( &v_start_up[ t ] );
    assert( var_index < qf->get_num_active_var() );
    qf->modify_linear_coefficient( var_index,
                                   get_start_up_cost( t ),
                                   issueAMod );
   }
  }
 }

 if( issue_pmod( issuePMod ) ) {
  Block::add_Modification(
   std::make_shared< ThermalUnitBlockRngdMod >(
    this,
    ThermalUnitBlockMod::eSetSUC,
    rng ),
   Observer::par2chnl( issuePMod ) );
 }
}

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_const_term(
 std::vector< double >::const_iterator values,
 Subset && subset,
 const bool ordered,
 c_ModParam issuePMod,
 c_ModParam issueAMod ) {

 if( subset.empty() ) {
  return;
 }

 if( v_ConstTerm.empty() ) {
  if( std::all_of( values, values + subset.size(),
                   []( double cst ) {
                    return ( cst == 0 );
                   } ) ) {
   return;
  }

  v_ConstTerm.assign( get_time_horizon() , 0 );
 }

 // If nothing changes, return
 bool identical = true;
 auto const_term = values;

 for( auto t : subset ) {
  if( t >= v_ConstTerm.size() ) {
   throw std::invalid_argument
    ( "ThermalUnitBlock::set_const_term: invalid index in subset: "
      + std::to_string( t ) );
  }

  if( v_ConstTerm[ t ] != *( const_term++ ) ) {
   identical = false;
   break;
  }
 }

 if( identical ) {
  return;
 }

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  const_term = values;
  for( auto t : subset ) {
   v_ConstTerm[ t ] = *( const_term++ );
  }

  if( not_dry_run( issueAMod ) && objective_generated() ) {
   // Change the abstract representation

   auto qf = dynamic_cast<DQuadFunction *>(objective.get_function());

   for( auto t : subset ) {

    auto var_index = qf->is_active( &v_commitment[ t ] );
    assert( var_index < qf->get_num_active_var() );
    qf->modify_linear_coefficient( var_index,
                                   get_const_term()[ t ],
                                   issueAMod );
   }
  }
 }

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification
  if( !ordered ) {
   std::sort( subset.begin(), subset.end() );
  }

  Block::add_Modification(
   std::make_shared< ThermalUnitBlockSbstMod >(
    this,
    ThermalUnitBlockMod::eSetConstT,
    std::move( subset ) ),
   Observer::par2chnl( issuePMod ) );
 }
}

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_const_term(
 std::vector< double >::const_iterator values,
 Range rng,
 c_ModParam issuePMod,
 c_ModParam issueAMod ) {

 rng.second = std::min( rng.second, f_time_horizon );
 if( rng.second <= rng.first ) {
  return;
 }

 if( v_ConstTerm.empty() ) {
  if( std::all_of( values,
                   values + ( rng.second - rng.first ),
                   []( double cst ) {
                    return ( cst == 0 );
                   } ) ) {
   return;
  }

  v_ConstTerm.assign( get_time_horizon() , 0 );
 }

 if( rng.first >= v_ConstTerm.size() ) {
  throw std::invalid_argument
   ( "ThermalUnitBlock::set_const_term: invalid first endpoint of "
     "range: " + std::to_string( rng.first ) );
 }

 // If nothing changes, return
 if( std::equal( values,
                 values + ( rng.second - rng.first ),
                 v_ConstTerm.begin() + rng.first ) ) {
  return;
 }

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  std::copy( values,
             values + ( rng.second - rng.first ),
             v_ConstTerm.begin() + rng.first );

  if( not_dry_run( issueAMod ) && objective_generated() ) {
   // Change the abstract representation

   auto qf = dynamic_cast<DQuadFunction *>(objective.get_function());

   for( Index t = rng.first; t < rng.second; ++t ) {

    auto var_index = qf->is_active( &v_commitment[ t ] );
    assert( var_index < qf->get_num_active_var() );
    qf->modify_linear_coefficient( var_index,
                                   get_const_term()[ t ],
                                   issueAMod );
   }
  }
 }

 if( issue_pmod( issuePMod ) ) {
  Block::add_Modification(
   std::make_shared< ThermalUnitBlockRngdMod >(
    this,
    ThermalUnitBlockMod::eSetConstT,
    rng ),
   Observer::par2chnl( issuePMod ) );
 }
}

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_linear_term(
 std::vector< double >::const_iterator values,
 Subset && subset,
 const bool ordered,
 c_ModParam issuePMod,
 c_ModParam issueAMod ) {

 if( subset.empty() ) {
  return;
 }

 if( v_LinearTerm.empty() ) {
  if( std::all_of( values, values + subset.size(),
                   []( double cst ) {
                    return ( cst == 0 );
                   } ) ) {
   return;
  }

  v_LinearTerm.assign( get_time_horizon() , 0 );
 }

 // If nothing changes, return
 bool identical = true;
 auto linear_term = values;

 for( auto t : subset ) {
  if( t >= v_LinearTerm.size() ) {
   throw std::invalid_argument
    ( "ThermalUnitBlock::set_linear_term: invalid index in subset: "
      + std::to_string( t ) );
  }

  if( v_LinearTerm[ t ] != *( linear_term++ ) ) {
   identical = false;
   break;
  }
 }

 if( identical ) {
  return;
 }

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  linear_term = values;
  for( auto t : subset ) {
   v_LinearTerm[ t ] = *( linear_term++ );
  }

  if( not_dry_run( issueAMod ) && objective_generated() ) {
   // Change the abstract representation

   auto qf = dynamic_cast<DQuadFunction *>(objective.get_function());

   for( auto t : subset ) {

    auto var_index = qf->is_active( &v_active_power[ t ] );
    assert( var_index < qf->get_num_active_var() );
    qf->modify_linear_coefficient( var_index,
                                   get_linear_term()[ t ],
                                   issueAMod );
   }
  }
 }

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification
  if( !ordered ) {
   std::sort( subset.begin(), subset.end() );
  }

  Block::add_Modification(
   std::make_shared< ThermalUnitBlockSbstMod >(
    this,
    ThermalUnitBlockMod::eSetLinT,
    std::move( subset ) ),
   Observer::par2chnl( issuePMod ) );
 }
}

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_linear_term(
 std::vector< double >::const_iterator values,
 Range rng,
 c_ModParam issuePMod,
 c_ModParam issueAMod ) {

 rng.second = std::min( rng.second, f_time_horizon );
 if( rng.second <= rng.first ) {
  return;
 }

 if( v_LinearTerm.empty() ) {
  if( std::all_of( values,
                   values + ( rng.second - rng.first ),
                   []( double cst ) {
                    return ( cst == 0 );
                   } ) ) {
   return;
  }

  v_LinearTerm.assign( get_time_horizon() , 0 );
 }

 if( rng.first >= v_LinearTerm.size() ) {
  throw std::invalid_argument
   ( "ThermalUnitBlock::set_const_term: invalid first endpoint of "
     "range: " + std::to_string( rng.first ) );
 }

 // If nothing changes, return
 if( std::equal( values,
                 values + ( rng.second - rng.first ),
                 v_LinearTerm.begin() + rng.first ) ) {
  return;
 }

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  std::copy( values,
             values + ( rng.second - rng.first ),
             v_LinearTerm.begin() + rng.first );

  if( not_dry_run( issueAMod ) && objective_generated() ) {
   // Change the abstract representation

   auto qf = dynamic_cast<DQuadFunction *>(objective.get_function());

   for( Index t = rng.first; t < rng.second; ++t ) {

    auto var_index = qf->is_active( &v_active_power[ t ] );
    assert( var_index < qf->get_num_active_var() );
    qf->modify_linear_coefficient( var_index,
                                   get_linear_term()[ t ],
                                   issueAMod );
   }
  }
 }

 if( issue_pmod( issuePMod ) ) {
  Block::add_Modification(
   std::make_shared< ThermalUnitBlockRngdMod >(
    this,
    ThermalUnitBlockMod::eSetLinT,
    rng ),
   Observer::par2chnl( issuePMod ) );
 }
}

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_quad_term(
 std::vector< double >::const_iterator values,
 Subset && subset,
 const bool ordered,
 c_ModParam issuePMod,
 c_ModParam issueAMod ) {

 if( subset.empty() ) {
  return;
 }

 if( v_QuadTerm.empty() ) {
  if( std::all_of( values, values + subset.size(),
                   []( double cst ) {
                    return ( cst == 0 );
                   } ) ) {
   return;
  }

  v_QuadTerm.assign( get_time_horizon() , 0 );
 }

 // If nothing changes, return
 bool identical = true;
 auto quad_term = values;

 for( auto t : subset ) {
  if( t >= v_QuadTerm.size() ) {
   throw std::invalid_argument
    ( "ThermalUnitBlock::set_quad_term: invalid index in subset: "
      + std::to_string( t ) );
  }

  if( v_QuadTerm[ t ] != *( quad_term++ ) ) {
   identical = false;
   break;
  }
 }

 if( identical ) {
  return;
 }

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  quad_term = values;
  for( auto t : subset ) {
   v_QuadTerm[ t ] = *( quad_term++ );
  }

  if( not_dry_run( issueAMod ) && objective_generated() ) {
   // Change the abstract representation

   auto qf = dynamic_cast<DQuadFunction *>(objective.get_function());

   for( auto t : subset ) {

    auto var_index = qf->is_active( &v_active_power[ t ] );
    assert( var_index < qf->get_num_active_var() );
    qf->modify_term( var_index,
                     get_linear_term()[ t ],
                     get_quad_term()[ t ],
                     issueAMod );
   }
  }
 }

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification
  if( !ordered ) {
   std::sort( subset.begin(), subset.end() );
  }

  Block::add_Modification(
   std::make_shared< ThermalUnitBlockSbstMod >(
    this,
    ThermalUnitBlockMod::eSetQuadT,
    std::move( subset ) ),
   Observer::par2chnl( issuePMod ) );
 }
}

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_quad_term(
 std::vector< double >::const_iterator values,
 Range rng,
 c_ModParam issuePMod,
 c_ModParam issueAMod ) {


 rng.second = std::min( rng.second, f_time_horizon );
 if( rng.second <= rng.first ) {
  return;
 }

 if( v_QuadTerm.empty() ) {
  if( std::all_of( values,
                   values + ( rng.second - rng.first ),
                   []( double cst ) {
                    return ( cst == 0 );
                   } ) ) {
   return;
  }

  v_QuadTerm.assign( get_time_horizon() , 0 );
 }

 if( rng.first >= v_QuadTerm.size() ) {
  throw std::invalid_argument
   ( "ThermalUnitBlock::set_quad_term: invalid first endpoint of "
     "range: " + std::to_string( rng.first ) );
 }

 // If nothing changes, return
 if( std::equal( values,
                 values + ( rng.second - rng.first ),
                 v_QuadTerm.begin() + rng.first ) ) {
  return;
 }

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  std::copy( values,
             values + ( rng.second - rng.first ),
             v_QuadTerm.begin() + rng.first );

  if( not_dry_run( issueAMod ) && objective_generated() ) {
   // Change the abstract representation

   auto qf = dynamic_cast<DQuadFunction *>(objective.get_function());

   for( Index t = rng.first; t < rng.second; ++t ) {

    auto var_index = qf->is_active( &v_active_power[ t ] );
    assert( var_index < qf->get_num_active_var() );
    qf->modify_term( var_index,
                     get_linear_term()[ t ],
                     get_quad_term()[ t ],
                     issueAMod );
   }
  }
 }

 if( issue_pmod( issuePMod ) ) {
  Block::add_Modification(
   std::make_shared< ThermalUnitBlockRngdMod >(
    this,
    ThermalUnitBlockMod::eSetQuadT,
    rng ),
   Observer::par2chnl( issuePMod ) );
 }
}

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_primary_spinning_reserve_cost
( std::vector< double >::const_iterator values , Subset && subset ,
  const bool ordered , c_ModParam issuePMod , c_ModParam issueAMod ) {

 if( subset.empty() )
  return;

 if( v_primary_spinning_reserve_cost.empty() ) {
  // The primary spinning reserve costs are currently all zero.
  if( std::all_of( values , values + subset.size() ,
                   []( double cst ) { return ( cst == 0 ); } ) ) {
   return; // The given values are zero. Nothing to do.
  }

  v_primary_spinning_reserve_cost.assign( get_time_horizon() , 0 );
 }

 // If nothing changes, return
 bool identical = true;
 auto cost_coefficient = values;

 for( auto t : subset ) {
  if( t >= v_primary_spinning_reserve_cost.size() ) {
   throw std::invalid_argument( "ThermalUnitBlock::set_primary_spinning_"
                                "reserve_cost: invalid index in subset: "
                                + std::to_string( t ) );
  }

  if( v_primary_spinning_reserve_cost[ t ] != *( cost_coefficient++ ) ) {
   identical = false;
   break;
  }
 }

 if( identical )
  return; // The given coefficients are equal to the ones already
          // here. So, there is nothing to be changed.

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  cost_coefficient = values;
  for( auto t : subset ) {
   v_primary_spinning_reserve_cost[ t ] = *( cost_coefficient++ );
  }

  if( not_dry_run( issueAMod ) && objective_generated() ) {
   // Change the abstract representation

   auto qf = dynamic_cast<DQuadFunction *>( objective.get_function() );

   for( auto t : subset ) {
    auto var_index = qf->is_active( &v_primary_spinning_reserve[ t ] );
    assert( var_index < qf->get_num_active_var() );
    qf->modify_linear_coefficient
     ( var_index , v_primary_spinning_reserve_cost[ t ] , issueAMod );
   }
  }
 }

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification
  if( !ordered ) {
   std::sort( subset.begin(), subset.end() );
  }

  Block::add_Modification( std::make_shared< ThermalUnitBlockSbstMod >
                           ( this , ThermalUnitBlockMod::eSetPrSpResCost ,
                             std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );
 }
} // end( ThermalUnitBlock::set_primary_spinning_reserve_cost )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_primary_spinning_reserve_cost
( std::vector< double >::const_iterator values , Range rng ,
  c_ModParam issuePMod , c_ModParam issueAMod ) {

 rng.second = std::min( rng.second, f_time_horizon );
 if( rng.second <= rng.first ) {
  return; // Empty range. Return.
 }

 if( v_primary_spinning_reserve_cost.empty() ) {
  // The primary spinning reserve costs are currently all zero.
  if( std::all_of( values , values + ( rng.second - rng.first ) ,
                   []( double cst ) { return ( cst == 0 ); } ) ) {
   return; // The given values are zero. So, there is nothing to be changed.
  }

  v_primary_spinning_reserve_cost.assign( get_time_horizon() , 0 );
 }

 if( rng.first >= v_primary_spinning_reserve_cost.size() ) {
  throw std::invalid_argument( "ThermalUnitBlock::set_primary_spinning_reserve"
                               "_cost: invalid first endpoint of range: " +
                               std::to_string( rng.first ) );
 }

 // If nothing changes, return
 if( std::equal( values , values + ( rng.second - rng.first ) ,
                 v_primary_spinning_reserve_cost.begin() + rng.first ) ) {
  return;
 }

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  std::copy( values , values + ( rng.second - rng.first ) ,
             v_primary_spinning_reserve_cost.begin() + rng.first );

  if( not_dry_run( issueAMod ) && objective_generated() ) {
   // Change the abstract representation

   auto qf = dynamic_cast<DQuadFunction *>( objective.get_function() );

   for( Index t = rng.first ; t < rng.second ; ++t ) {

    auto var_index = qf->is_active( &v_commitment[ t ] );
    assert( var_index < qf->get_num_active_var() );
    qf->modify_linear_coefficient
     ( var_index , v_primary_spinning_reserve_cost[ t ] , issueAMod );
   }
  }
 }

 if( issue_pmod( issuePMod ) ) {
  Block::add_Modification( std::make_shared< ThermalUnitBlockRngdMod >
                           ( this, ThermalUnitBlockMod::eSetPrSpResCost , rng ),
                           Observer::par2chnl( issuePMod ) );
 }
} // end( ThermalUnitBlock::set_primary_spinning_reserve_cost )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_secondary_spinning_reserve_cost
( std::vector< double >::const_iterator values , Subset && subset ,
  const bool ordered , c_ModParam issuePMod , c_ModParam issueAMod ) {

 if( subset.empty() )
  return;

 if( v_secondary_spinning_reserve_cost.empty() ) {
  // The secondary spinning reserve costs are currently all zero.
  if( std::all_of( values , values + subset.size() ,
                   []( double cst ) { return ( cst == 0 ); } ) ) {
   return; // The given values are zero. Nothing to do.
  }

  v_secondary_spinning_reserve_cost.assign( get_time_horizon() , 0 );
 }

 // If nothing changes, return
 bool identical = true;
 auto cost_coefficient = values;

 for( auto t : subset ) {
  if( t >= v_secondary_spinning_reserve_cost.size() ) {
   throw std::invalid_argument( "ThermalUnitBlock::set_secondary_spinning_"
                                "reserve_cost: invalid index in subset: "
                                + std::to_string( t ) );
  }

  if( v_secondary_spinning_reserve_cost[ t ] != *( cost_coefficient++ ) ) {
   identical = false;
   break;
  }
 }

 if( identical )
  return; // The given coefficients are equal to the ones already
          // here. So, there is nothing to be changed.

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  cost_coefficient = values;
  for( auto t : subset ) {
   v_secondary_spinning_reserve_cost[ t ] = *( cost_coefficient++ );
  }

  if( not_dry_run( issueAMod ) && objective_generated() ) {
   // Change the abstract representation

   auto qf = dynamic_cast<DQuadFunction *>( objective.get_function() );

   for( auto t : subset ) {
    auto var_index = qf->is_active( &v_secondary_spinning_reserve[ t ] );
    assert( var_index < qf->get_num_active_var() );
    qf->modify_linear_coefficient
     ( var_index , v_secondary_spinning_reserve_cost[ t ] , issueAMod );
   }
  }
 }

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification
  if( !ordered ) {
   std::sort( subset.begin(), subset.end() );
  }

  Block::add_Modification( std::make_shared< ThermalUnitBlockSbstMod >
                           ( this , ThermalUnitBlockMod::eSetSecSpResCost ,
                             std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );
 }
} // end( ThermalUnitBlock::set_secondary_spinning_reserve_cost )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_secondary_spinning_reserve_cost
( std::vector< double >::const_iterator values , Range rng ,
  c_ModParam issuePMod , c_ModParam issueAMod ) {

 rng.second = std::min( rng.second, f_time_horizon );
 if( rng.second <= rng.first ) {
  return; // Empty range. Return.
 }

 if( v_secondary_spinning_reserve_cost.empty() ) {
  // The secondary spinning reserve costs are currently all zero.
  if( std::all_of( values , values + ( rng.second - rng.first ) ,
                   []( double cst ) { return ( cst == 0 ); } ) ) {
   return; // The given values are zero. So, there is nothing to be changed.
  }

  v_secondary_spinning_reserve_cost.assign( get_time_horizon() , 0 );
 }

 if( rng.first >= v_secondary_spinning_reserve_cost.size() ) {
  throw std::invalid_argument( "ThermalUnitBlock::set_secondary_spinning_"
                               "reserve_cost: invalid first endpoint of "
                               "range: " + std::to_string( rng.first ) );
 }

 // If nothing changes, return
 if( std::equal( values , values + ( rng.second - rng.first ) ,
                 v_secondary_spinning_reserve_cost.begin() + rng.first ) ) {
  return;
 }

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  std::copy( values , values + ( rng.second - rng.first ) ,
             v_secondary_spinning_reserve_cost.begin() + rng.first );

  if( not_dry_run( issueAMod ) && objective_generated() ) {
   // Change the abstract representation

   auto qf = dynamic_cast<DQuadFunction *>( objective.get_function() );

   for( Index t = rng.first ; t < rng.second ; ++t ) {

    auto var_index = qf->is_active( &v_commitment[ t ] );
    assert( var_index < qf->get_num_active_var() );
    qf->modify_linear_coefficient
     ( var_index , v_secondary_spinning_reserve_cost[ t ] , issueAMod );
   }
  }
 }

 if( issue_pmod( issuePMod ) ) {
  Block::add_Modification( std::make_shared< ThermalUnitBlockRngdMod >
                           ( this , ThermalUnitBlockMod::eSetSecSpResCost ,
                             rng ) , Observer::par2chnl( issuePMod ) );
 }
} // end( ThermalUnitBlock::set_secondary_spinning_reserve_cost )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_init_updown_time
( std::vector< int >::const_iterator values , Block::Subset && subset ,
  const bool ordered , c_ModParam issuePMod , c_ModParam issueAMod ) {

 if( subset.empty() )
  return;

 // Find the last index 0
 auto index_it = std::find( subset.rbegin() , subset.rend() , 0 );

 if( index_it == subset.rend() )
  return; // 0 is not in subset; return

 std::advance( values , std::distance( index_it , subset.rend() ) - 1 );

 if( f_InitUpDownTime == *values )
  return; // nothing changes; return

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation
  f_InitUpDownTime = *values;

  if( not_dry_run( issueAMod ) && variables_generated() ) {
   // TODO
   throw( std::logic_error( "ThermalUnitBlock::set_init_updown_time: it is "
                            "currently not possible to update the abstract "
                            "representation." ) );
  }
 }

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockMod >
                           ( this , ThermalUnitBlockMod::eSetInitUD ) ,
                           Observer::par2chnl( issuePMod ) );
 }
}  // end( ThermalUnitBlock::set_init_updown_time )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_init_updown_time(
              std::vector< int >::const_iterator values , Block::Range rng ,
              c_ModParam issuePMod , c_ModParam issueAMod )
{
 rng.second = std::min( rng.second , decltype( rng.second )( 1 ) );
 if( ! ( rng.first <= 0 && 0 < rng.second ) )
  return; // 0 does not belong to the range; return

 std::advance( values , - rng.first );

 if( f_InitUpDownTime == *values )
  return; // nothing changes; return

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation
  f_InitUpDownTime = *values;

  if( not_dry_run( issueAMod ) && variables_generated() ) {
   // TODO
   throw( std::logic_error( "ThermalUnitBlock::set_init_updown_time: it is "
                            "currently not possible to update the abstract "
                            "representation." ) );
  }
 }

 if( issue_pmod( issuePMod ) ) {
  Block::add_Modification( std::make_shared< ThermalUnitBlockMod >
                           ( this , ThermalUnitBlockMod::eSetInitUD ),
                           Observer::par2chnl( issuePMod ) );
 }
}  // end( ThermalUnitBlock::set_init_updown_time )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::scale
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
    update_objective( Range( 0 , Inf<Index>() ) , issueAMod );
  }
 }

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< UnitBlockMod >
                           ( this , UnitBlockMod::eScale ) ,
                           Observer::par2chnl( issuePMod ) );
 }
}  // end( ThermalUnitBlock::scale )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::update_objective_start_up( const Subset & subset ,
                                                  c_ModParam issueAMod ) {

 if( ! objective_generated() )
  return; // the Objective has not been generated: nothing to be done

 auto function = dynamic_cast<DQuadFunction *>( objective.get_function() );

 if( ! function )
  return;

 for( auto t : subset ) {
  if( t < init_t )
   continue;
  auto var_index = function->is_active( &v_start_up[ t ] );
  assert( var_index < function->get_num_active_var() );
  function->modify_linear_coefficient
   ( var_index , f_scale * get_start_up_cost( t ) , issueAMod );
 }
}  // end( ThermalUnitBlock::update_objective_start_up )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::update_objective_active_power( const Subset & subset ,
                                                      c_ModParam issueAMod ) {

 if( ! objective_generated() )
  return; // the Objective has not been generated: nothing to be done

 auto function = dynamic_cast<DQuadFunction *>( objective.get_function() );

 if( ! function )
  return;

 for( auto t : subset ) {
  auto var_index = function->is_active( &v_active_power[ t ] );
  assert( var_index < function->get_num_active_var() );
  function->modify_term( var_index , f_scale *  get_linear_term( t ),
                         f_scale *  get_quad_term( t ), issueAMod );
 }
}  // end( ThermalUnitBlock::update_objective_active_power )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::update_objective_commitment( const Subset & subset ,
                                                    c_ModParam issueAMod ) {

 if( ! objective_generated() )
  return; // the Objective has not been generated: nothing to be done

 auto function = dynamic_cast<DQuadFunction *>( objective.get_function() );

 if( ! function )
  return;

 for( auto t : subset ) {
  auto var_index = function->is_active( &v_commitment[ t ] );
  assert( var_index < function->get_num_active_var() );
  function->modify_linear_coefficient
   ( var_index , f_scale * get_const_term( t ) , issueAMod );
 }
}  // end( ThermalUnitBlock::update_objective_commitment )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::update_objective( const Subset & subset ,
                                         c_ModParam issueAMod ) {
 update_objective_start_up( subset , issueAMod );
 update_objective_active_power( subset , issueAMod );
 update_objective_commitment( subset , issueAMod );
}  // end( ThermalUnitBlock::update_objective )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::update_objective( Range rng , c_ModParam issueAMod ) {

 rng.second = std::min( rng.second , f_time_horizon );
 if( rng.second <= rng.first )
  return;

 Subset subset( rng.second - rng.first );
 std::iota( subset.begin() , subset.end() , rng.first );

 update_objective( subset , issueAMod );
}  // end( ThermalUnitBlock::update_objective )

/*--------------------------------------------------------------------------*/

template< typename T >
void ThermalUnitBlock::decompress_vector( std::vector< T > & v )
{
 if( v.empty() )
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
           ( "ThermalUnitBlock::decompress_vector: invalid number of elements"
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
  v.resize( f_time_horizon );
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
} // end( ThermalUnitBlock::decompress_vector )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::guts_of_add_Modification( p_Mod mod , ChnlName chnl )
{
 // process abstract Modification - - - - - - - - - - - - - - - - - - - - - -
 /* This requires to patiently sift through the possible Modification types
  * to find what this Modification exactly is and appropriately mirror the
  * changes to the "abstract representation" to the "physical one".
  *
  * Note that since ThermalUnitBlock is a "leaf" Block (has no sub-Block),
  * this method does not have to deal with GroupModification since these
  * are produced by Block::add_Modification(), but this method is called
  * *before* that one is.
  *
  * As an important consequence,
  *
  *   THE STATE OF THE DATA STRUCTURE IN ThermalUnitBlock WHEN THIS METHOD
  *   IS EXECUTED IS PRECISELY THE ONE IN WHICH THE Modification WAS
  *   ISSUED: NO COMPLICATED OPERATIONS (Variable AND/OR Constraint BEING
  *   ADDED/REMOVED ...) CAN HAVE BEEN PERFORMED IN THE MEANTIME
  *
  * This assumption drastically simplifies some logic here.*/

 // VariableMod - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( const auto tmod = dynamic_cast< VariableMod * >( mod ) ) {
  // changing the Variable is not supported, but the Modification is issued
  // when they are first generated, in which case it must be ignored
  if( ! variables_generated() )
   return;

  throw( std::logic_error( "ThermalUnitBlock - VariableMod not supported"
                           ) );
  /*
  auto v = dynamic_cast< ColVariable * const >( tmod->variable() );

  if( v->is_fixed() ) {
   // TODO: Do something to the physical representation

   } else {
   // TODO: Do something to the physical representation
   }
  */
  return;
  }

 // BlockMod - Generic modification - - - - - - - - - - - - - - - - - - - - -
 if( const auto tmod = dynamic_cast< BlockMod * >( mod ) ) {
  // changing the Objective is not supported, but the Modification is issued
  // when it is first set, in which case it must be ignored
  if( ! objective_generated() )
   return;

  throw( std::logic_error( "ThermalUnitBlock - BlockMod not supported" ) );

  // TODO: BlockMod - obj changed
  return;
  }

 // FunctionMod - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( const auto tmod = dynamic_cast< FunctionMod * >( mod ) ) {
  auto f = tmod->function();
  if( f == static_cast< FRealObjective * >( get_objective() )->get_function()
      ) {
   handle_objective_change( tmod , chnl );
   return;
   }

  std::ostringstream em;
  em << *mod;
  throw( std::invalid_argument( "ThermalUnitBlock: unsupported " + em.str()
                                ) );
  return;
  }

 // any other Modification is not supported - - - - - - - - - - - - - - - - -

 std::ostringstream em;
 em << *mod;
 throw( std::invalid_argument( "ThermalUnitBlock: unsupported " + em.str()
                               ) );

 }  // end( ThermalUnitBlock::guts_of_add_Modification )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::handle_objective_change( FunctionMod * mod ,
                                                ChnlName chnl )
{
 const auto * qf = static_cast< const DQuadFunction * >( mod->function() );
 auto par = make_par( eNoBlck , chnl );
 Index th = get_time_horizon();

 // C05FunctionModLinRngd - - - - - - - - - - - - - - - - - - - - - - - - - -
 // split the C05FunctionModLinRngd in up to 5 physical Modification by
 // calling the appropriate set_*() methods (ranged version) for those among
 // startup, power, commitment, primary/secondary reserve variables whose
 // coefficient change. this heavily relies on the fact that variables of
 // the same type are consecutive (and ordered in the obvious way) when
 // set as coefficients in the Objective

 if( const auto tmod = dynamic_cast< C05FunctionModLinRngd * >( mod ) ) {

  Index l = tmod->range().first;
  Index r = tmod->range().second;

  if( tmod->range().second > qf->get_num_active_var() )
   throw( std::invalid_argument(
         "ThermalUnitBlock: invalid Range [" + std::to_string( l ) + ", "
         + std::to_string( r ) + ") in C05FunctionModLinRngd" ) );

  std::vector< double > nv( r - l + 1 );
  Index gl = 0;
  Index gr = th - init_t;

  if( l < gr ) {  // startup variables
   Index r2 = std::min( r , gr );
   auto nvit = nv.begin();
   for( Index i = l ; i < r2 ; )
    *(nvit++) = qf->get_linear_coefficient( i++ );
   set_startup_costs( nv.begin() , Range( l , r2 ) , par , eDryRun );
   l = r2;
   if( l == r )
    return;
   }

  gl = gr;
  gr = 2 * th - init_t;

  if( l < gr ) {  // active power variables
   Index r2 = std::min( r , gr );
   auto nvit = nv.begin();
   for( Index i = l ; i < r2 ; )
    *(nvit++) = qf->get_linear_coefficient( i++ );
   set_linear_term( nv.begin() , Range( l - gl , r2 - gl ) , par , eDryRun );
   l = r2;
   if( l == r )
    return;
   }

  gl = gr;
  gr = 3 * th - init_t;

  if( l < gr ) {  // commitment variables
   Index r2 = std::min( r , gr );
   auto nvit = nv.begin();
   for( Index i = l ; i < r2 ; )
    *(nvit++) = qf->get_linear_coefficient( i++ );
   set_const_term( nv.begin() , Range( l - gl , r2 - gl ) , par , eDryRun );
   l = r2;
   if( l == r )
    return;
   }

  gl = gr;
  gr = 4 * th - init_t;

  if( l < gr ) {  // primary spinning reserve variables
   Index r2 = std::min( r , gr );
   auto nvit = nv.begin();
   for( Index i = l ; i < r2 ; )
    *(nvit++) = qf->get_linear_coefficient( i++ );
   set_primary_spinning_reserve_cost( nv.begin() , Range( l - gl , r2 - gl ) ,
                                      par , eDryRun );
   l = r2;
   if( l == r )
    return;
   }

  gl = gr;
  gr = 5 * th - init_t;

  if( l < gr ) {  // secondary spinning reserve variables
   Index r2 = std::min( r , gr );
   auto nvit = nv.begin();
   for( Index i = l ; i < r2 ; )
    *(nvit++) = qf->get_linear_coefficient( i++ );
   set_secondary_spinning_reserve_cost( nv.begin() , Range( l - gl , r2 - gl ) ,
                                        par , eDryRun );
   if( r2 == r )
    return;
   }

  throw( std::invalid_argument(
           "ThermalUnitBlock: invalid variable in C05FunctionModLinRngd" ) );
  return;

  }  // end( C05FunctionModLinRngd )

 // C05FunctionModLinSbst - - - - - - - - - - - - - - - - - - - - - - - - - -
 // split the C05FunctionModLinSbst in up to 5 physical Modification by
 // calling the appropriate set_*() methods (subset version) for those among
 // startup, power, commitment, primary/secondary reserve variables whose
 // coefficient change. this heavily relies on the fact that variables of
 // the same type are consecutive (and ordered in the obvious way) when
 // set as coefficients in the Objective

 if( const auto tmod = dynamic_cast< C05FunctionModLinSbst * >( mod ) ) {

  if( tmod->subset().back() > qf->get_num_active_var() )
   throw( std::invalid_argument(
              "ThermalUnitBlock: invalid Subset in C05FunctionModLinSbst" ) );

  std::vector< double > nv( tmod->subset().size() );
  auto l = tmod->subset().begin();
  Index gl = 0;
  Index gr = th - init_t;

  if( *l < gr ) {  // startup variables
   auto r = l;
   for( ++r ; ( r != tmod->subset().end() ) && ( *r < gr ) ; )
    ++r;
   Subset nms( std::distance( l , r ) );
   auto nvit = nv.begin();
   auto nmsit = nms.begin();
   while( l != r ) {
    *(nvit++) = qf->get_linear_coefficient( *l );
    *(nmsit++) = *(l++);
    }
   set_startup_costs( nv.begin() , std::move( nms ) , true , par , eDryRun );
   if( r == tmod->subset().end() )
    return;
   }

  gl = gr;
  gr = 2 * th - init_t;

  if( *l < gr ) {  // active power variables
   auto r = l;
   for( ++r ; ( r != tmod->subset().end() ) && ( *r < gr ) ; )
    ++r;
   Subset nms( std::distance( l , r ) );
   auto nvit = nv.begin();
   auto nmsit = nms.begin();
   while( l != r ) {
    *(nvit++) = qf->get_linear_coefficient( *l );
    *(nmsit++) = *(l++) - gl;
    }
   set_linear_term( nv.begin() , std::move( nms ) , true , par , eDryRun );
   if( r == tmod->subset().end() )
    return;
   }

  gl = gr;
  gr = 3 * th - init_t;

  if( *l < gr ) {  // commitment variables
   auto r = l;
   for( ++r ; ( r != tmod->subset().end() ) && ( *r < gr ) ; )
    ++r;
   Subset nms( std::distance( l , r ) );
   auto nvit = nv.begin();
   auto nmsit = nms.begin();
   while( l != r ) {
    *(nvit++) = qf->get_linear_coefficient( *l );
    *(nmsit++) = *(l++) - gl;
    }
   set_const_term( nv.begin() , std::move( nms ) , true , par , eDryRun );
   if( r == tmod->subset().end() )
    return;
   }

  gl = gr;
  gr = 4 * th - init_t;

  if( *l < gr ) {  // primary spinning reserve variables
   auto r = l;
   for( ++r ; ( r != tmod->subset().end() ) && ( *r < gr ) ; )
    ++r;
   Subset nms( std::distance( l , r ) );
   auto nvit = nv.begin();
   auto nmsit = nms.begin();
   while( l != r ) {
    *(nvit++) = qf->get_linear_coefficient( *l );
    *(nmsit++) = *(l++) - gl;
    }
   set_primary_spinning_reserve_cost( nv.begin() , std::move( nms ) ,
                                      true , par , eDryRun );
   if( r == tmod->subset().end() )
    return;
   }

  gl = gr;
  gr = 5 * th - init_t;

  if( *l < gr ) {  // secondary spinning reserve variables
   auto r = l;
   for( ++r ; ( r != tmod->subset().end() ) && ( *r < gr ) ; )
    ++r;
   Subset nms( std::distance( l , r ) );
   auto nvit = nv.begin();
   auto nmsit = nms.begin();
   while( l != r ) {
    *(nvit++) = qf->get_linear_coefficient( *l );
    *(nmsit++) = *(l++) - gl;
    }
   set_secondary_spinning_reserve_cost( nv.begin() , std::move( nms ) ,
                                        true , par , eDryRun );
   }

  if( l != tmod->subset().end() )
   throw( std::invalid_argument(
           "ThermalUnitBlock: invalid variable in C05FunctionModLinSbst" ) );
  return;

  }  // end( C05FunctionModLinSbst )


 throw( std::invalid_argument(
             "ThermalUnitBlock:: unsupported FunctionMod from Objective" ) );

 }  // end( ThermalUnitBlock::handle_objective_change )

/*--------------------------------------------------------------------------*/
/*------------------- End File ThermalUnitBlock.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
