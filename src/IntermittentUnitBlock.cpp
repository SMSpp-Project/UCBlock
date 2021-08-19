/*--------------------------------------------------------------------------*/
/*------------- File IntermittentUnitBlock.cpp -------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the IntermittentUnitBlock class.
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
 * \copyright &copy; by Antonio Frangioni, Ali Ghezelsoflu, and Rafael Durbano
 * Lobato
 */

/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <iostream>
#include <random>
#include "IntermittentUnitBlock.h"
#include "LinearFunction.h"
#include "FRealObjective.h"
#include <map>
#include "UnitBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;


/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register IntermittentUnitBlock to the Block factory


SMSpp_insert_in_factory_cpp_1( IntermittentUnitBlock );

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS OF IntermittentUnitBlock -------------------*/
/*--------------------------------------------------------------------------*/

IntermittentUnitBlock::~IntermittentUnitBlock() {
 for( auto & constraint : MinPower_Constraints )
  constraint.clear();
 for( auto & constraint : MaxPower_Constraints )
  constraint.clear();
 for( auto & constraint : active_power_bounds_Constraints )
  constraint.clear();

 objective.clear();
}

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::deserialize( const netCDF::NcGroup & group ) {


#ifndef NDEBUG
 std::vector< std::string > expected_dims = { "TimeHorizon",
                                              "NumberIntervals" };
 check_dimensions( group, expected_dims, std::cerr );
 std::vector< std::string > expected_vars = { "MinPower",
                                              "MaxPower",
                                              "InertiaPower",
                                              "Gamma",
                                              "Kappa" };
 check_variables( group, expected_vars, std::cerr );
#endif



 UnitBlock::deserialize_time_horizon( group );
 UnitBlock::deserialize_change_intervals( group );

  ::deserialize( group, "MinPower", f_time_horizon, v_minimum_power, true, true );

  ::deserialize( group, "MaxPower",f_time_horizon,v_maximum_power, true, true );

 ::deserialize( group, "MaxPower",v_maximum_power, true );

 ::deserialize( group, "InertiaPower", v_inertia_power, true );

 ::deserialize( group, "Gamma", &f_gamma, false );

if (! ::deserialize( group, "Kappa", &f_kappa, true )) {
 f_kappa = 1;
}

 decompress_vector( v_minimum_power );
 decompress_vector( v_maximum_power );

 UnitBlock::deserialize( group );
}// end( IntermittentUnitBlock::deserialize )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::generate_abstract_variables( Configuration *stvv ) {

 if( variables_generated() )
  return; // variables have already been generated

 UnitBlock::generate_abstract_variables( stvv );

  // Active Power Variable

 v_active_power.resize( f_time_horizon );
 for( auto & var : v_active_power )
  var.set_type( ColVariable::kNonNegative );
 add_static_variable( v_active_power, "p_intermittent" );

  // Primary Spinning Reserve Variable
 if( reserve_vars & 1u ) { // if UCBlock has primary demand variables
  if ( f_gamma != 0 ) { // if unit produces any reserve
   v_primary_spinning_reserve.resize( f_time_horizon );
   for( auto & var : v_primary_spinning_reserve ) {
    var.set_type( ColVariable::kNonNegative );
   }
   add_static_variable( v_primary_spinning_reserve, "pr_intermittent" );
  }
 }

 // Secondary Spinning Reserve Variable
 if( reserve_vars & 2u ) { // if UCBlock has secondary demand variables
  if ( f_gamma != 0 ) { // if unit produces any reserve
   v_secondary_spinning_reserve.resize( f_time_horizon );
   for( auto & var : v_secondary_spinning_reserve ) {
    var.set_type( ColVariable::kNonNegative );
   }
   add_static_variable( v_secondary_spinning_reserve, "sr_intermittent" );
  }
 }


 set_variables_generated();
} // end( IntermittentUnitBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::generate_abstract_constraints
( Configuration *stcc ) {

 if( constraints_generated() )
  return; // constraints have already been generated

 std::vector<double> max_power = v_maximum_power;
 std::vector<double> min_power = v_minimum_power;

/*--------------------------------------------------------------------------*/
 // Initializing maximum power constraints

 if ( f_gamma != 0 ) { // if unit produces any reserve

  if( MaxPower_Constraints.size() != f_time_horizon ) {
   // this should only happen once
   assert( MaxPower_Constraints.empty());

   MaxPower_Constraints.resize( f_time_horizon );
  }

  for( Index t = 0; t < f_time_horizon; ++t ) {

   auto linear_function = new LinearFunction();
   linear_function->add_variable( &v_active_power[t], f_gamma );
   if( reserve_vars & 1u ) {
    linear_function->add_variable( &v_primary_spinning_reserve[t], 1.0 );
   }
   if( reserve_vars & 2u ) {
    linear_function->add_variable( &v_secondary_spinning_reserve[t], 1.0 );
   }
   MaxPower_Constraints[t].set_lhs( -Inf< double >());
   MaxPower_Constraints[t].set_rhs(( f_gamma * f_kappa * ( max_power[t] )));
   MaxPower_Constraints[t].set_function( linear_function );
  }

  add_static_constraint( MaxPower_Constraints, "MaxPower_Intermittent" );

 }
  // Initializing minimum power constraints

  if( MinPower_Constraints.size() != f_time_horizon ) {
   // this should only happen once
   assert( MinPower_Constraints.empty());

   MinPower_Constraints.resize( f_time_horizon );
  }

  for( Index t = 0; t < f_time_horizon; ++t ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_active_power[t], 1.0 );
   if( reserve_vars & 1u ) {
    if ( f_gamma != 0 ) { // if unit produces any reserve
     linear_function->add_variable( &v_primary_spinning_reserve[t], -1.0 );
    }
   }
   if( reserve_vars & 2u ) {
    if ( f_gamma != 0 ) { // if unit produces any reserve
     linear_function->add_variable( &v_secondary_spinning_reserve[t], -1.0 );
    }
   }

   MinPower_Constraints[t].set_lhs( f_kappa * min_power[t] );
   MinPower_Constraints[t].set_rhs( Inf< double >());
   MinPower_Constraints[t].set_function( linear_function );
  }

  add_static_constraint( MinPower_Constraints, "MinPower_Intermittent" );

 // Initializing active power bounds constraints
 if( active_power_bounds_Constraints.size() != f_time_horizon ) {
  // this should only happen once
  assert( active_power_bounds_Constraints.empty());

  active_power_bounds_Constraints.resize( f_time_horizon );
 }

 for( Index t = 0; t < f_time_horizon; ++t ) {

  active_power_bounds_Constraints[t].set_lhs( f_kappa * min_power[ t ]);
  active_power_bounds_Constraints[t].set_rhs( f_kappa * max_power[ t ] );
  active_power_bounds_Constraints[t].set_variable( &v_active_power[t] );
 }

 add_static_constraint( active_power_bounds_Constraints, "ActivePowerBound_Intermittent" );

 set_constraints_generated();
} // end( IntermittentUnitBlock::generate_abstract_constraints )

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

bool IntermittentUnitBlock::is_feasible( bool useabstract ,
                                         Configuration * fsbc ) {

 // Retrieve the tolerance.

 auto config = dynamic_cast< SimpleConfiguration< double > * >( fsbc );

 if( ( ! config ) && f_BlockConfig )
  config = dynamic_cast< SimpleConfiguration< double > * >
   ( f_BlockConfig->f_is_feasible_Configuration );

 // If a tolerance has not been provided, use the default tolerance.
 const auto tolerance = config ? config->f_value : 1.0e-8;

 return
  UnitBlock::is_feasible( useabstract )
  // Constraints
  && ::is_feasible( MinPower_Constraints , tolerance )
  && ::is_feasible( MaxPower_Constraints , tolerance )
  && ::is_feasible( active_power_bounds_Constraints , tolerance )
  // Variables
  && ::is_feasible( v_active_power , tolerance )
  && ::is_feasible( v_primary_spinning_reserve , tolerance )
  && ::is_feasible( v_secondary_spinning_reserve , tolerance );

} // end( IntermittentUnitBlock::is_feasible )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::generate_objective( Configuration * objc ) {

 if( objective_generated() )
  return; // Objective has already been generated

 if( get_objective() != nullptr )  // an objective is there already
  return;                         // cowardly (and silently) return

 auto linear_function = new LinearFunction();
 objective.set_function( linear_function );
 // Set Block objective
 this->set_objective( &objective );

 set_objective_generated();

}  // end( IntermittentUnitBlock::generate_objective )
/*--------------------------------------------------------------------------*/
/*--- METHODS FOR LOADING, PRINTING & SAVING THE IntermittentUnitBlock -----*/
/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::serialize( netCDF::NcGroup & group ) const {

 UnitBlock::serialize( group );

 auto NumberIntervals = group.getDim( "NumberIntervals" );

 ::serialize( group, "Gamma", netCDF::NcDouble(), f_gamma );

 ::serialize( group, "Kappa", netCDF::NcDouble(), f_kappa );

 ::serialize( group, "MinPower", netCDF::NcDouble(),
              NumberIntervals, v_minimum_power, true );

 ::serialize( group, "MaxPower", netCDF::NcDouble(),
              NumberIntervals, v_maximum_power, true );

 ::serialize( group, "InertiaPower", netCDF::NcDouble(),
              { NumberIntervals }, v_inertia_power, true );
}  // end( IntermittentUnitBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::update_max_power_in_constraints
( const Block::Subset & time , c_ModParam issueAMod ) {

 if( ! MaxPower_Constraints.empty() ) {
  for( auto t : time ) {
   MaxPower_Constraints[ t ].set_rhs
    ( f_kappa * f_gamma * v_maximum_power[ t ] , issueAMod );
   // FIXME: use a GroupModification
  }
 }
 if( ! active_power_bounds_Constraints.empty() ) {
  for( auto t : time ) {
   active_power_bounds_Constraints[ t ].set_rhs
    ( f_kappa * v_maximum_power[ t ] , issueAMod );
   // FIXME: use a GroupModification
  }
 }
}

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::update_max_power_in_constraints
( const Block::Range & time , c_ModParam issueAMod ) {

 if( ! MaxPower_Constraints.empty() ) {
  for( auto t = time.first ; t < time.second ; ++t ) {
   MaxPower_Constraints[ t ].set_rhs
    ( f_kappa * f_gamma * v_maximum_power[ t ] , issueAMod );
   // FIXME: use a GroupModification
  }
 }
 if( ! active_power_bounds_Constraints.empty() ) {
  for( auto t = time.first ; t < time.second ; ++t ) {
   active_power_bounds_Constraints[ t ].set_rhs
    ( f_kappa * v_maximum_power[ t ] , issueAMod );
   // FIXME: use a GroupModification
  }
 }
}

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::set_maximum_power
( std::vector< double >::const_iterator values , Block::Subset && subset ,
  const bool ordered , c_ModParam issuePMod , c_ModParam issueAMod ) {

 if( subset.empty() ) {
  return;
 }

 if( v_maximum_power.empty() ) {
  if( std::all_of( values , values + subset.size() ,
                   []( double cst ) { return ( cst == 0 ); } ) ) {
   return;
  }

  Index max_index = * std::max_element( std::begin( subset ) ,
                                        std::end( subset ) );
  v_maximum_power.assign( max_index , 0 );
 }

 // If nothing changes, return
 bool identical = true;
 for( auto t : subset ) {
  if( t >= v_maximum_power.size() ) {
   throw ( std::invalid_argument( "IntermittentUnitBlock::set_maximum_power:"
                                  " invalid value in subset" ) );
  }
  auto max_power = *( values++ );
  if( v_maximum_power[ t ] != max_power ) {
   identical = false;
   if( not_dry_run( issuePMod ) )
    // Change the physical representation
    v_maximum_power[ t ] = max_power;
  }
 }
 if( identical )
  return;  // nothing changes; return

 if( not_dry_run( issuePMod ) && not_dry_run( issueAMod ) &&
     constraints_generated() ) {
  // Change the abstract representation
  update_max_power_in_constraints( subset , issueAMod );
 }

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification
  if( !ordered ) {
   std::sort( subset.begin(), subset.end() );
  }

  Block::add_Modification( std::make_shared< IntermittentUnitBlockSbstMod >
                           ( this, IntermittentUnitBlockMod::eSetMaxP ,
                             std::move( subset ) ),
                           Observer::par2chnl( issuePMod ) );
 }
}

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::set_maximum_power
( std::vector< double >::const_iterator values , Block::Range rng ,
  c_ModParam issuePMod , c_ModParam issueAMod ) {

 rng.second = std::min( rng.second , f_time_horizon );
 if( rng.second <= rng.first ) {
  return;
 }

 if( v_maximum_power.empty() ) {
  if( std::all_of( values , values + ( rng.second - rng.first ) ,
                   []( double cst ) { return ( cst == 0 ); } ) ) {
   return;
  }

  auto max_index = rng.second;
  v_maximum_power.assign( max_index , 0 );
 }

 // If nothing changes, return
 if( std::equal( values , values + ( rng.second - rng.first ) ,
                 v_maximum_power.begin() + rng.first ) ) {
  return;
 }

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  std::copy( values , values + ( rng.second - rng.first ) ,
             v_maximum_power.begin() + rng.first );

  if( not_dry_run( issueAMod ) && constraints_generated() ) {
   // Change the abstract representation
   update_max_power_in_constraints( rng , issueAMod );
  }
 }

 if( issue_pmod( issuePMod ) ) {
  Block::add_Modification( std::make_shared< IntermittentUnitBlockRngdMod >
                           ( this , IntermittentUnitBlockMod::eSetMaxP , rng ) ,
                           Observer::par2chnl( issuePMod ) );
 }
}

/*--------------------------------------------------------------------------*/

template< typename T >
void IntermittentUnitBlock::decompress_vector( std::vector< T > & v ) {
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
/*------------------- End File IntermittentUnitBlock.cpp -------------------*/
/*--------------------------------------------------------------------------*/
