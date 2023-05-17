/*--------------------------------------------------------------------------*/
/*-------------------- File BusNetworkBlock.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the BusNetworkBlock class.
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
 * \author Kostas Tavlaridis-Gyparakis \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Ali Ghezelsoflu,
 *                      Rafael Durbano Lobato,
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <map>
#include "NetworkBlock.h"
#include "BusNetworkBlock.h"
#include "LinearFunction.h"
#include "FRealObjective.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register NetworkBlock to the Block factory

SMSpp_insert_in_factory_cpp_1( BusNetworkBlock );

/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

void BusNetworkBlock::generate_abstract_variables( Configuration * stvv )
{
 if( variables_generated() )
  return; // variables have already been generated

 // In BusNetworkBlock, number_nodes = 1
 auto active_demand = get_active_demand()[ 0 ];
 v_node_injection[ 0 ].set_value( active_demand );
 v_node_injection[ 0 ].is_fixed( true );
 v_node_injection[ 0 ].set_type( ColVariable::kContinuous );
 add_static_variable( v_node_injection[ 0 ], "S");

 set_variables_generated();
 }

/*--------------------------------------------------------------------------*/

void BusNetworkBlock::generate_abstract_constraints ( Configuration * stcc )
{
 if( constraints_generated())
  return; // constraints have already been generated

 // the node injection bound constraints
 NodeInjection_bound_Constraints.resize( 1 );
 auto active_demand = get_active_demand()[ 0 ];
 NodeInjection_bound_Constraints[ 0 ].set_lhs( active_demand );
 NodeInjection_bound_Constraints[ 0 ].set_rhs( active_demand );
 NodeInjection_bound_Constraints[ 0 ].set_variable(&v_node_injection[0]);

 add_static_constraint( NodeInjection_bound_Constraints,
                        "NodeInjection_bound_BusNetwork" );
 }

/*--------------------------------------------------------------------------*/

void BusNetworkBlock::generate_objective( Configuration * objc )
{
 if( objective_generated() )
  return; // Objective has already been generated

 if( get_objective() != nullptr )  // an objective is there already
  return;                         // cowardly (and silently) return

 auto linear_function = new LinearFunction();

 objective.set_function( linear_function );

 // Set Block objective
 this->set_objective( &objective );

 set_objective_generated();

 }  // end( BusNetworkBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

void BusNetworkBlock::set_active_demand( MF_dbl_it values ,
					 Subset && subset , bool ordered ,
					 ModParam issuePMod ,
					 ModParam issueAMod )
{
 // For BusNetworkBlock, this method degenerates a bit

 if( subset.empty() )
  return;

 if (subset.size() > 1 )
  throw ( std::invalid_argument( "subset is too big" ) );

 if( v_active_demand.empty() ) {
  if( *values == 0 )
   return;

  v_active_demand.assign( 1, 0 );
  }

 if( subset[0] >= 1 )
  throw ( std::invalid_argument( "invalid value in subset" ) );

 // If nothing changes, return
 if( v_active_demand[ 0 ] == *values )
  return;

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  v_active_demand[ 0 ] = *values ;

  if( not_dry_run( issueAMod ) && variables_generated() ) {
   // Change the abstract representation

   v_node_injection[ 0 ].set_value( v_active_demand[ 0 ] );
   }
  }

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification

  Block::add_Modification(
   std::make_shared< NetworkBlockSbstMod >( this ,
                                            NetworkBlockMod::eSetActD,
                                            std::move( subset ) ) ,
   Observer::par2chnl( issuePMod ) );
  }
 }

/*--------------------------------------------------------------------------*/

void BusNetworkBlock::set_active_demand( MF_dbl_it values , Range rng ,
					 ModParam issuePMod ,
					 ModParam issueAMod )
{
 // For BusNetworkBlock, this method degenerates a bit

 if( rng != Range( 0, 1 ) )
  throw ( std::invalid_argument( "invalid value in subset" ) );

 if( v_active_demand.empty() ) {
  if( *values == 0 )
   return;

  v_active_demand.assign( 1, 0 );
  }

 // If nothing changes, return
 if( v_active_demand[ 0 ] == *values )
  return;

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  v_active_demand[ 0 ] = *values ;

  if( not_dry_run( issueAMod ) && variables_generated() )
   // Change the abstract representation
   v_node_injection[ 0 ].set_value( v_active_demand[ 0 ] );
  }

 if( issue_pmod( issuePMod ) ) {
  Block::add_Modification(
   std::make_shared< NetworkBlockRngdMod >( this ,
                                            NetworkBlockMod::eSetActD ,
                                            rng ) ,
   Observer::par2chnl( issuePMod ) );
  }
 }

/*--------------------------------------------------------------------------*/
/*--------------------- End File BusNetworkBlock.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
