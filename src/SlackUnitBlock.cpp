/*--------------------------------------------------------------------------*/
/*------------------------ File SlackUnitBlock.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the SlackUnitBlock class, which derives from UnitBlock
 * [see UnitBlock.h] and implements a "slack" unit; a (typically, fictitious)
 * unit capable of producing (typically, a large amount of) active power
 * and/or primary/secondary reserve and/or inertia at any time period
 * completely indeopendently from each other and from all other time periods,
 * albeit at a (typically, huge) cost.
 *
 * \version 0.11
 *
 * \date 23 - 09 - 2021
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
 * Copyright &copy by Antonio Frangioni, Ali Ghezelsoflu
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "SlackUnitBlock.h"
#include "LinearFunction.h"
#include "FRealObjective.h"
/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register SlackUnitBlock to the Block factory


SMSpp_insert_in_factory_cpp_1( SlackUnitBlock );

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS OF SlackUnitBlock ------------------------*/
/*--------------------------------------------------------------------------*/
SlackUnitBlock::~SlackUnitBlock() {
 auto clear_LB0Constraints =
         []( std::vector< LB0Constraint > & constraints ) {
          for( auto & constraint : constraints )
           constraint.clear();
         };
 clear_LB0Constraints( Secondary_Spinning_Reserve_Bound_Constraints );
 clear_LB0Constraints( Primary_Spinning_Reserve_Bound_Constraints );
 clear_LB0Constraints( ActivePower_Bound_Constraints );

 auto clear_ZOConstraints =
         []( std::vector< ZOConstraint > & constraints ) {
          for( auto & constraint : constraints )
           constraint.clear();
         };
 clear_ZOConstraints( Inertia_Bound_Constraints );

 objective.clear();

}

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
void SlackUnitBlock::deserialize( const netCDF::NcGroup & group ) {



#ifndef NDEBUG
 std::vector< std::string > expected_dims = { "TimeHorizon",
                                              "NumberIntervals" };
 check_dimensions( group, expected_dims, std::cerr );
 std::vector< std::string > expected_vars = { "MaxPower",
                                              "MaxPrimaryPower",
                                              "MaxSecondaryPower",
                                              "ActivePowerCost",
                                              "PrimaryCost",
                                              "SecondaryCost",
                                              "InertiaCost",
                                              "MaxInertia"};
 check_variables( group, expected_vars, std::cerr );
#endif


 UnitBlock::deserialize_time_horizon( group );
 UnitBlock::deserialize_change_intervals( group );

 ::deserialize( group, "MaxPower",f_time_horizon,v_MaxPower, true, true);
 ::deserialize( group, "MaxPrimaryPower",f_time_horizon,v_MaxPrimaryPower, true,true);
 ::deserialize( group, "MaxSecondaryPower",f_time_horizon,v_MaxSecondaryPower, true,true);
 ::deserialize( group, "ActivePowerCost",f_time_horizon,v_active_power_cost, true,true );
 ::deserialize( group, "PrimaryCost",f_time_horizon,v_primary_cost, true,true);
 ::deserialize( group, "SecondaryCost",f_time_horizon,v_secondary_cost, true,true );
 ::deserialize( group, "InertiaCost",f_time_horizon,v_inertia_cost, true,true);
 ::deserialize( group, "MaxInertia",f_time_horizon, v_MaxInertia, true,true );

 decompress_vector(v_MaxPower);
 decompress_vector(v_MaxPrimaryPower);
 decompress_vector(v_MaxSecondaryPower);
 decompress_vector(v_active_power_cost);
 decompress_vector(v_primary_cost);
 decompress_vector(v_secondary_cost);
 decompress_vector(v_inertia_cost);
 decompress_vector(v_MaxInertia);

 UnitBlock::deserialize( group );
}// end( SlackUnitBlock::deserialize )

/*--------------------------------------------------------------------------*/

void SlackUnitBlock::generate_abstract_variables
        ( Configuration *stvv )
{

 if( variables_generated() )
  return; // variables have already been generated

/*--------------------------------------------------------------------------*/
  // Commitment Variable
 if( reserve_vars & 4u ) { // if UCBlock has inertia demand variables
  if( !v_MaxInertia.empty()) { // if unit produces any inertia reserve
   v_commitment.resize( f_time_horizon );
   for( auto & i : v_commitment ) {
    i.set_type( ColVariable::kPosUnitary );
   }
   add_static_variable( v_commitment, "u_inertia_slack" );
  }
 }
  // Active Power Variable
  v_active_power.resize( f_time_horizon );
  for( auto & var : v_active_power )
   var.set_type( ColVariable::kNonNegative );
  add_static_variable( v_active_power, "p_slack" );

  // Primary Spinning Reserve Variable
 if( reserve_vars & 1u ) { // if UCBlock has primary demand variables
  if( !v_MaxPrimaryPower.empty()) { // if unit produces any primary reserve
   v_primary_spinning_reserve.resize( f_time_horizon );
   for( auto & var : v_primary_spinning_reserve ) {
    var.set_type( ColVariable::kNonNegative );
   }
   add_static_variable( v_primary_spinning_reserve, "pr_slack" );
  }
 }
  // Secondary Spinning Reserve Variable
 if( reserve_vars & 2u ) { // if UCBlock has secondary demand variables
  if( !v_MaxSecondaryPower.empty()) { // if unit produces any secondary reserve
   v_secondary_spinning_reserve.resize( f_time_horizon );
   for( auto & var : v_secondary_spinning_reserve ) {
    var.set_type( ColVariable::kNonNegative );
   }
   add_static_variable( v_secondary_spinning_reserve, "sr_slack" );
  }
 }
 set_variables_generated();

} // end( SlackUnitBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void SlackUnitBlock::generate_abstract_constraints
        ( Configuration *stcc ) {

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

 // Initializing active power bounds constraints
 if( ActivePower_Bound_Constraints.size() != f_time_horizon ) {
  // this should only happen once
  assert( ActivePower_Bound_Constraints.empty());

  ActivePower_Bound_Constraints.resize( f_time_horizon );
 }

 for( Index t = 0; t < f_time_horizon; ++t ) {

  if( !v_MaxPower.empty()) {
   ActivePower_Bound_Constraints[t].set_rhs( v_MaxPower[t] );
  } else {
   ActivePower_Bound_Constraints[t].set_rhs( 0.0 );
  }
  ActivePower_Bound_Constraints[t].set_variable( &v_active_power[t] );
 }

 add_static_constraint( ActivePower_Bound_Constraints, "ActivePowerBound_Slack" );
/*--------------------------------------------------------------------------*/

 // Initializing primary spinning reserve bounds constraints
 if( reserve_vars & 1u ) {
  if(!v_MaxPrimaryPower.empty()) {

   if( Primary_Spinning_Reserve_Bound_Constraints.size() != f_time_horizon ) {
   // this should only happen once
   assert( Primary_Spinning_Reserve_Bound_Constraints.empty());

   Primary_Spinning_Reserve_Bound_Constraints.resize( f_time_horizon );
  }

  for( Index t = 0; t < f_time_horizon; ++t ) {
   Primary_Spinning_Reserve_Bound_Constraints[t].set_rhs( v_MaxPrimaryPower[t] );
   Primary_Spinning_Reserve_Bound_Constraints[t].set_variable( &v_primary_spinning_reserve[t] );
  }

  add_static_constraint( Primary_Spinning_Reserve_Bound_Constraints, "PrimarySpinningReserveBound_Slack" );
 }
}
/*--------------------------------------------------------------------------*/

 // Initializing secondary spinning reserve bounds constraints
 if( reserve_vars & 2u ) {
  if( !v_MaxSecondaryPower.empty()) {
   if( Secondary_Spinning_Reserve_Bound_Constraints.size() != f_time_horizon ) {
    // this should only happen once
    assert( Secondary_Spinning_Reserve_Bound_Constraints.empty());

    Secondary_Spinning_Reserve_Bound_Constraints.resize( f_time_horizon );
   }

   for( Index t = 0; t < f_time_horizon; ++t ) {
    Secondary_Spinning_Reserve_Bound_Constraints[t].set_rhs( v_MaxSecondaryPower[t] );
    Secondary_Spinning_Reserve_Bound_Constraints[t].set_variable( &v_secondary_spinning_reserve[t] );
   }

   add_static_constraint( Secondary_Spinning_Reserve_Bound_Constraints, "SecondarySpinningReserveBound_Slack" );
  }
 }
 /*-------------------------------ZOConstraint-------------------------------*/

 if( generate_ZOConstraint ) {

  // the commitment bound constraints
  if( reserve_vars & 4u ) {
   Inertia_Bound_Constraints.resize( f_time_horizon );
   for( Index t = 0; t < f_time_horizon; ++t ) {
    Inertia_Bound_Constraints[t].set_variable( &v_commitment[t] );
   }
   add_static_constraint( Inertia_Bound_Constraints, "Inertia_bound_Slack" );
  }
 }

 set_constraints_generated();
} // end( SlackUnitBlock::generate_abstract_constraints )


/*--------------------------------------------------------------------------*/

void SlackUnitBlock::generate_objective( Configuration *objc )
{
 if( get_objective() != nullptr )  // an objective is there already
  return;                        // cowardly (and silently) return

 // Initialize objective function
 if( reserve_vars & 4u ) {
  if( v_commitment.size() != f_time_horizon ) {
   throw ( std::logic_error
           ( "SlackUnitBlock::generate_objective: v_commitment must have "
             "size equal to the time horizon." ));
  }
 }
 if( v_active_power.size() != f_time_horizon ) {
  throw ( std::logic_error
          ( "SlackUnitBlock::generate_objective: v_active_power must have "
            "size equal to the time horizon." ) );
 }
 if( reserve_vars & 1u ) {
  if( !v_MaxPrimaryPower.empty()) {
   if( v_primary_spinning_reserve.size() != f_time_horizon ) {
    throw ( std::logic_error
            ( "SlackUnitBlock::generate_objective: v_primary_spinning_reserve "
              "must have size equal to the time horizon." ));
   }
  }
 }
 if( reserve_vars & 2u ) {
  if( !v_MaxSecondaryPower.empty()) {
   if( v_secondary_spinning_reserve.size() != f_time_horizon ) {
    throw ( std::logic_error
            ( "SlackUnitBlock::generate_objective: v_secondary_spinning_reserve"
              "must have size equal to the time horizon." ));
   }
  }
 }
 auto linear_function = new LinearFunction();

 for( Index t = 0; t < f_time_horizon; ++t ) {

  if( !v_active_power_cost.empty()) {
   linear_function->add_variable( &v_active_power[t],
                                  v_active_power_cost[t],
                                  0.0 );
  } else {
   linear_function->add_variable( &v_active_power[t],
                                  0.0,
                                  0.0 );
  }
  if( reserve_vars & 1u ) {
   if( !v_MaxPrimaryPower.empty()) {
    if( !v_primary_cost.empty()) {
     linear_function->add_variable( &v_primary_spinning_reserve[t],
                                    v_primary_cost[t],
                                    0.0 );
    } else {
     linear_function->add_variable( &v_primary_spinning_reserve[t],
                                    0.0,
                                    0.0 );
    }
   }
  }
  if( reserve_vars & 2u ) {
   if( !v_MaxSecondaryPower.empty()) {

    if( !v_secondary_cost.empty()) {
     linear_function->add_variable( &v_secondary_spinning_reserve[t],
                                    v_secondary_cost[t],
                                    0.0 );
    } else {
     linear_function->add_variable( &v_secondary_spinning_reserve[t],
                                    0.0,
                                    0.0 );
    }
   }
  }
  if( reserve_vars & 4u ) {
   if( !v_inertia_cost.empty() && !v_MaxInertia.empty()) {
    linear_function->add_variable( &v_commitment[t],
                                   v_inertia_cost[t] * v_MaxInertia[t],
                                   0.0 );
   } else {
    linear_function->add_variable( &v_commitment[t],
                                   0.0,
                                   0.0 );
   }
  }
 }
 objective.set_function( linear_function );
 objective.set_sense( Objective::eMin );

 // Set Block objective
 this->set_objective( &objective );

 set_objective_generated();
}  // end( SlackUnitBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*------- METHODS FOR LOADING, PRINTING & SAVING THE SlackUnitBlock --------*/
/*--------------------------------------------------------------------------*/

void SlackUnitBlock::serialize( netCDF::NcGroup & group ) const {

 UnitBlock::serialize( group );

 auto NumberIntervals = group.getDim( "NumberIntervals" );

 ::serialize( group, "MaxPower", netCDF::NcDouble(),
              v_MaxPower );

 ::serialize( group, "MaxPrimaryPower", netCDF::NcDouble(),
               v_MaxPrimaryPower);

 ::serialize( group, "MaxSecondaryPower", netCDF::NcDouble(),
              NumberIntervals, v_MaxSecondaryPower, true );

 ::serialize( group, "ActivePowerCost", netCDF::NcDouble(),
              NumberIntervals, v_active_power_cost, true );

 ::serialize( group, "PrimaryCost", netCDF::NcDouble(),
              NumberIntervals, v_primary_cost, true );

 ::serialize( group, "SecondaryCost", netCDF::NcDouble(),
              NumberIntervals, v_secondary_cost, true );

 ::serialize( group, "InertiaCost", netCDF::NcDouble(),
              NumberIntervals, v_inertia_cost, true );

 ::serialize( group, "MaxInertia", netCDF::NcDouble(),
              { NumberIntervals }, v_MaxInertia, true );
}  // end( SlackUnitBlock::serialize )

/*--------------------------------------------------------------------------*/

template< typename T >
void SlackUnitBlock::decompress_vector( std::vector< T > & v ) {
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
}

/*--------------------------------------------------------------------------*/
/*----------------------- End File SlackUnitBlock.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
