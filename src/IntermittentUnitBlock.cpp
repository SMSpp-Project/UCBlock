/*--------------------------------------------------------------------------*/
/*------------- File IntermittentUnitBlock.cpp -------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the IntermittentUnitBlock class.
 *
 * \version 0.11
 *
 * \date 15 - 02 - 2021
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

 v_primary_spinning_reserve.resize( f_time_horizon );
 for( auto & var : v_primary_spinning_reserve )
  var.set_type( ColVariable::kNonNegative );
 if ( f_gamma != 0 ) {
  add_static_variable( v_primary_spinning_reserve, "pr_intermittent" );
 }

 // Secondary Spinning Reserve Variable

 v_secondary_spinning_reserve.resize( f_time_horizon );
 for( auto & var : v_secondary_spinning_reserve )
  var.set_type( ColVariable::kNonNegative );
 if ( f_gamma != 0 ) {
  add_static_variable( v_secondary_spinning_reserve, "sr_intermittent" );
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

 if ( f_gamma != 0 ) {

  if( MaxPower_Constraints.size() != f_time_horizon ) {
   // this should only happen once
   assert( MaxPower_Constraints.empty());

   MaxPower_Constraints.resize( f_time_horizon );
  }

  for( Index t = 0; t < f_time_horizon; ++t ) {

   auto linear_function = new LinearFunction();
   linear_function->add_variable( &v_active_power[t], f_gamma );
   linear_function->add_variable( &v_primary_spinning_reserve[t], 1.0 );
   linear_function->add_variable( &v_secondary_spinning_reserve[t], 1.0 );

   MaxPower_Constraints[t].set_lhs( -Inf< double >());
   MaxPower_Constraints[t].set_rhs(( f_gamma * f_kappa * ( max_power[t] )));
   MaxPower_Constraints[t].set_function( linear_function );
  }

  add_static_constraint( MaxPower_Constraints, "MaxPower_Intermittent" );


  // Initializing minimum power constraints

  if( MinPower_Constraints.size() != f_time_horizon ) {
   // this should only happen once
   assert( MinPower_Constraints.empty());

   MinPower_Constraints.resize( f_time_horizon );
  }

  for( Index t = 0; t < f_time_horizon; ++t ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_active_power[t], 1.0 );
   linear_function->add_variable( &v_primary_spinning_reserve[t], -1.0 );
   linear_function->add_variable( &v_secondary_spinning_reserve[t], -1.0 );

   MinPower_Constraints[t].set_lhs( f_kappa * min_power[t] );
   MinPower_Constraints[t].set_rhs( Inf< double >());
   MinPower_Constraints[t].set_function( linear_function );
  }

  add_static_constraint( MinPower_Constraints, "MinPower_Intermittent" );
 }
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
/*------------------- End File IntermittentUnitBlock.cpp -------------------*/
/*--------------------------------------------------------------------------*/
