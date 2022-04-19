/*--------------------------------------------------------------------------*/
/*-------------------- File IntermittentUnitBlock.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the IntermittentUnitBlock class.
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
 * \copyright &copy; by Antonio Frangioni, Ali Ghezelsoflu,
 *                   Rafael Durbano Lobato
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <iostream>
#include <random>
#include <map>

#include "IntermittentUnitBlock.h"
#include "LinearFunction.h"
#include "FRealObjective.h"
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
 clear_constraints( MinPower_Constraints );
 clear_constraints( MaxPower_Constraints );
 clear_constraints( active_power_bounds_Constraints );

 objective.clear();
}

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::deserialize( const netCDF::NcGroup & group ) {

#ifndef NDEBUG
 static std::vector< std::string > expected_dims = { "TimeHorizon" ,
                                                     "NumberIntervals" };
 check_dimensions( group , expected_dims , std::cerr );

 static std::vector< std::string > expected_vars = { "MinPower" ,
                                                     "MaxPower" ,
                                                     "InertiaPower" ,
                                                     "Gamma" ,
                                                     "Kappa" };
 check_variables( group , expected_vars , std::cerr );
#endif

 // Deserialize data that is needed for deserializing the variables
 UnitBlock::deserialize_time_horizon( group );

 // Mandatory variables
 ::deserialize( group , "MinPower" , v_minimum_power , false );
 ::deserialize( group , "MaxPower" , v_maximum_power , false );
 ::deserialize( group , f_gamma , "Gamma" , false );

 // Optional variables
 if( !::deserialize( group , "InertiaPower" , v_inertia_power ) )
  v_inertia_power.assign( f_time_horizon , 0 );

 if( !::deserialize( group , f_kappa , "Kappa" ) )
  f_kappa = 1;

 // Deserialize data from the base class
 UnitBlock::deserialize( group );

 // Decompress vectors
 decompress_vector( v_minimum_power );
 decompress_vector( v_maximum_power );
 decompress_vector( v_inertia_power );

 check_data_consistency();

}  // end( IntermittentUnitBlock::deserialize )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::check_data_consistency( void ) const {
 // Minimum and maximum power

 assert( v_minimum_power.size() == f_time_horizon );
 assert( v_maximum_power.size() == f_time_horizon );

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {
  if( v_minimum_power[ t ] > v_maximum_power[ t ] ) {
   throw ( std::logic_error( "IntermittentUnitBlock::check_data_consistency: "
                             "minimum power at time " + std::to_string( t ) +
                             " is " + std::to_string( v_minimum_power[ t ] ) +
                             ", which is greater than the maximum power, which "
                             "is " + std::to_string( v_maximum_power[ t ] ) +
                             "." ) );
  }

  if( v_minimum_power[ t ] < 0 ) {
   throw ( std::logic_error( "IntermittentUnitBlock::check_data_consistency: "
                             "minimum power at time " + std::to_string( t ) +
                             " is " + std::to_string( v_minimum_power[ t ] ) +
                             ", which is negative." ) );
  }
 }

 // Gamma

 if( ( f_gamma < 0 ) || ( f_gamma > 1 ) ) {
  throw ( std::logic_error( "IntermittentUnitBlock::check_data_consistency: "
                            "gamma must be between 0 and 1, but it is " +
                            std::to_string( f_gamma ) + "." ) );
 }

 // Kappa

 if( f_kappa < 0 ) {
  throw ( std::logic_error( "IntermittentUnitBlock::check_data_consistency: "
                            "kappa must be nonnegative, but it is" +
                            std::to_string( f_kappa ) + "." ) );
 }

 if( !v_inertia_power.empty() ) {
  assert( v_inertia_power.size() == f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   if( v_inertia_power[ t ] < 0 ) {
    throw ( std::logic_error( "IntermittentUnitBlock::check_data_consistency: "
                              "inertia power for time " + std::to_string( t ) +
                              " must be nonnegative, but it is" +
                              std::to_string( v_inertia_power[ t ] ) + "." ) );
   }
  }
 }
}  // end( IntermittentUnitBlock::check_data_consistency )

/*--------------------------------------------------------------------------*/

void
IntermittentUnitBlock::generate_abstract_variables( Configuration * stvv ) {
 if( variables_generated() )
  return; // variables have already been generated

 UnitBlock::generate_abstract_variables( stvv );

 // Active Power Variable

 v_active_power.resize( f_time_horizon );
 for( auto & var : v_active_power )
  var.set_type( ColVariable::kNonNegative );
 add_static_variable( v_active_power , "p_intermittent" );

 // Primary Spinning Reserve Variable

 if( reserve_vars & 1u ) { // if UCBlock has primary demand variables
  if( f_gamma != 0 ) { // if unit produces any reserve
   v_primary_spinning_reserve.resize( f_time_horizon );
   for( auto & var : v_primary_spinning_reserve ) {
    var.set_type( ColVariable::kNonNegative );
   }
   add_static_variable( v_primary_spinning_reserve , "pr_intermittent" );
  }
 }

 // Secondary Spinning Reserve Variable
 if( reserve_vars & 2u ) { // if UCBlock has secondary demand variables
  if( f_gamma != 0 ) { // if unit produces any reserve
   v_secondary_spinning_reserve.resize( f_time_horizon );
   for( auto & var : v_secondary_spinning_reserve ) {
    var.set_type( ColVariable::kNonNegative );
   }
   add_static_variable( v_secondary_spinning_reserve , "sr_intermittent" );
  }
 }

 set_variables_generated();

}  // end( IntermittentUnitBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::generate_abstract_constraints(
 Configuration * stcc ) {
 if( constraints_generated() )
  return; // constraints have already been generated

 std::vector< double > max_power = v_maximum_power;
 std::vector< double > min_power = v_minimum_power;

/*--------------------------------------------------------------------------*/
 // Initializing maximum power constraints

 if( f_gamma != 0 ) { // if unit produces any reserve

  if( MaxPower_Constraints.size() != f_time_horizon ) {
   // this should only happen once
   assert( MaxPower_Constraints.empty() );

   MaxPower_Constraints.resize( f_time_horizon );
  }

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {

   auto linear_function = new LinearFunction();
   linear_function->add_variable( &v_active_power[ t ] , f_gamma );
   if( reserve_vars & 1u ) {
    linear_function->add_variable( &v_primary_spinning_reserve[ t ] , 1.0 );
   }
   if( reserve_vars & 2u ) {
    linear_function->add_variable( &v_secondary_spinning_reserve[ t ] , 1.0 );
   }
   MaxPower_Constraints[ t ].set_lhs( -Inf< double >() );
   MaxPower_Constraints[ t ].set_rhs(
    ( f_gamma * f_kappa * ( max_power[ t ] ) ) );
   MaxPower_Constraints[ t ].set_function( linear_function );
  }

  add_static_constraint( MaxPower_Constraints , "MaxPower_Intermittent" );

 }
 // Initializing minimum power constraints

 if( MinPower_Constraints.size() != f_time_horizon ) {
  // this should only happen once
  assert( MinPower_Constraints.empty() );

  MinPower_Constraints.resize( f_time_horizon );
 }

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {

  auto linear_function = new LinearFunction();

  linear_function->add_variable( &v_active_power[ t ] , 1.0 );
  if( reserve_vars & 1u ) {
   if( f_gamma != 0 ) { // if unit produces any reserve
    linear_function->add_variable( &v_primary_spinning_reserve[ t ] , -1.0 );
   }
  }
  if( reserve_vars & 2u ) {
   if( f_gamma != 0 ) { // if unit produces any reserve
    linear_function->add_variable( &v_secondary_spinning_reserve[ t ] , -1.0 );
   }
  }

  MinPower_Constraints[ t ].set_lhs( f_kappa * min_power[ t ] );
  MinPower_Constraints[ t ].set_rhs( Inf< double >() );
  MinPower_Constraints[ t ].set_function( linear_function );
 }

 add_static_constraint( MinPower_Constraints , "MinPower_Intermittent" );

 // Initializing active power bounds constraints
 if( active_power_bounds_Constraints.size() != f_time_horizon ) {
  // this should only happen once
  assert( active_power_bounds_Constraints.empty() );

  active_power_bounds_Constraints.resize( f_time_horizon );
 }

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {

  active_power_bounds_Constraints[ t ].set_lhs( f_kappa * min_power[ t ] );
  active_power_bounds_Constraints[ t ].set_rhs( f_kappa * max_power[ t ] );
  active_power_bounds_Constraints[ t ].set_variable( &v_active_power[ t ] );
 }

 add_static_constraint( active_power_bounds_Constraints ,
                        "ActivePowerBound_Intermittent" );

 set_constraints_generated();

}  // end( IntermittentUnitBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

/// verifies whether the current solution is feasible for the given constraints
/** This function checks whether the relative violation of each RowConstraint
 * in the given group of RowConstraint is not greater than the provided
 * tolerance.
 *
 * @return This function returns true if and only if the relative violation of
 *         each RowConstraint in the given group is not greater than the given
 *         tolerance. */

template< class C >
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

template< class V >
static std::enable_if_t< std::is_base_of_v< ColVariable , V > , bool >
is_feasible( const std::vector< V > & variables , double tolerance ) {
 return std::all_of( variables.begin() , variables.end() ,
                     [ tolerance ]( const auto & variable ) {
                      return variable.is_feasible( tolerance );
                     } );
}

/*--------------------------------------------------------------------------*/

bool IntermittentUnitBlock::is_feasible( bool useabstract ,
                                         Configuration * fsbc ) {
 // Retrieve the tolerance.

 auto config = dynamic_cast< SimpleConfiguration< double > * >( fsbc );

 if( ( !config ) && f_BlockConfig )
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

}  // end( IntermittentUnitBlock::is_feasible )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::generate_objective( Configuration * objc ) {
 if( objective_generated() )
  return; // Objective has already been generated

 if( get_objective() != nullptr )  // an objective is there already
  return;                          // cowardly (and silently) return

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

 // Serialize scalar variables.
 ::serialize( group , "Gamma" , netCDF::NcDouble() , f_gamma );
 ::serialize( group , "Kappa" , netCDF::NcDouble() , f_kappa );

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
   throw ( std::logic_error
    ( "IntermittentUnitBlock::serialize: invalid dimension for variable " +
      var_name + ": " + std::to_string( data.size() ) + ". Its dimension " +
      "must be one of the following: TimeHorizon, NumberIntervals, 1." ) );
  }

  ::serialize( group , var_name , ncType , dimension , data ,
               allow_scalar_var );
 };

 serialize( "MinPower" , v_minimum_power );
 serialize( "MaxPower" , v_maximum_power );
 serialize( "InertiaPower" , v_inertia_power );

}  // end( IntermittentUnitBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::update_max_power_in_constraints(
 const Subset & time ,
 ModParam issueAMod ) {
 if( !MaxPower_Constraints.empty() ) {
  for( auto t : time ) {
   MaxPower_Constraints[ t ].set_rhs
    ( f_kappa * f_gamma * v_maximum_power[ t ] , issueAMod );
   // FIXME: use a GroupModification
  }
 }

 if( !active_power_bounds_Constraints.empty() ) {
  for( auto t : time ) {
   active_power_bounds_Constraints[ t ].set_rhs
    ( f_kappa * v_maximum_power[ t ] , issueAMod );
   // FIXME: use a GroupModification
  }
 }
}

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::update_max_power_in_constraints(
 const Range & time , ModParam issueAMod ) {
 if( !MaxPower_Constraints.empty() ) {
  for( auto t = time.first ; t < time.second ; ++t ) {
   MaxPower_Constraints[ t ].set_rhs
    ( f_kappa * f_gamma * v_maximum_power[ t ] , issueAMod );
   // FIXME: use a GroupModification
  }
 }

 if( !active_power_bounds_Constraints.empty() ) {
  for( auto t = time.first ; t < time.second ; ++t ) {
   active_power_bounds_Constraints[ t ].set_rhs
    ( f_kappa * v_maximum_power[ t ] , issueAMod );
   // FIXME: use a GroupModification
  }
 }
}

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::set_maximum_power( MF_dbl_it values ,
                                               Subset && subset ,
                                               bool ordered ,
                                               ModParam issuePMod ,
                                               ModParam issueAMod ) {
 if( subset.empty() )
  return;

 if( v_maximum_power.empty() ) {
  if( std::all_of( values , values + subset.size() ,
                   []( double cst ) { return( cst == 0 ); } ) ) {
   return;
  }

  Index max_index = *std::max_element( std::begin( subset ) ,
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
   std::sort( subset.begin() , subset.end() );
  }

  Block::add_Modification( std::make_shared< IntermittentUnitBlockSbstMod >
                            ( this , IntermittentUnitBlockMod::eSetMaxP ,
                              std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );
 }
}

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::set_maximum_power( MF_dbl_it values , Range rng ,
                                               ModParam issuePMod ,
                                               ModParam issueAMod ) {
 rng.second = std::min( rng.second , f_time_horizon );
 if( rng.second <= rng.first )
  return;

 if( v_maximum_power.empty() ) {
  if( std::all_of( values , values + ( rng.second - rng.first ) ,
                   []( double cst ) { return( cst == 0 ); } ) ) {
   return;
  }

  auto max_index = rng.second;
  v_maximum_power.assign( max_index , 0 );
 }

 // If nothing changes, return
 if( std::equal( values , values + ( rng.second - rng.first ) ,
                 v_maximum_power.begin() + rng.first ) )
  return;

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  std::copy( values , values + ( rng.second - rng.first ) ,
             v_maximum_power.begin() + rng.first );

  if( not_dry_run( issueAMod ) && constraints_generated() ) {
   // Change the abstract representation
   update_max_power_in_constraints( rng , issueAMod );
  }
 }

 if( issue_pmod( issuePMod ) )
  Block::add_Modification( std::make_shared< IntermittentUnitBlockRngdMod >(
                            this , IntermittentUnitBlockMod::eSetMaxP , rng ) ,
                           Observer::par2chnl( issuePMod ) );
}

/*--------------------------------------------------------------------------*/
/*------------------- End File IntermittentUnitBlock.cpp -------------------*/
/*--------------------------------------------------------------------------*/
