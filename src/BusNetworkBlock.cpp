/*--------------------------------------------------------------------------*/
/*-------------------- File BusNetworkBlock.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the BusNetworkBlock class.
 *
 * \version 0.11
 *
 * \date 01 - 07 - 2019
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
 * \author Kostas Tavlaridis-Gyparakis \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy by Antonio Frangioni, Ali Ghezelsoflu, Rafael
 * Durbano Lobato, and Kostas Tavlaridis-Gyparakis
 */

/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <map>
#include "NetworkBlock.h"
#include "BusNetworkBlock.h"

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
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

void BusNetworkBlock::generate_abstract_variables( Configuration * stvv ) {

 if( AR & HasVar )
  return; // variables have already been generated

 // In BusNetworkBlock, number_nodes = 1
 auto active_demand = get_active_demand()[ 0 ];
 v_node_injection[ 0 ].set_value( active_demand );
 v_node_injection[ 0 ].is_fixed( true );
 v_node_injection[ 0 ].set_type( ColVariable::kContinuous );
 add_static_variable( v_node_injection[ 0 ], "S");

 AR |= HasVar;
}

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

void BusNetworkBlock::set_active_demand(
 std::vector< double >::const_iterator values,
 Block::Subset && subset,
 const bool ordered,
 c_ModParam issuePMod,
 c_ModParam issueAMod ) {

 // For BusNetworkBlock, this method degenerates a bit

 if( subset.empty() ) {
  return;
 }

 if (subset.size() > 1 ) {
  throw ( std::invalid_argument( "subset is too big" ) );
 }

 if( v_active_demand.empty() ) {
  if( *values == 0 ) {
   return;
  }
  v_active_demand.assign( 1, 0 );
 }

 if( subset[0] >= 1 ) {
  throw ( std::invalid_argument( "invalid value in subset" ) );
 }

 // If nothing changes, return
 if( v_active_demand[ 0 ] == *values ) {
  return;
 }

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  v_active_demand[ 0 ] = *values ;

  if( not_dry_run( issueAMod ) && AR & HasVar ) {
   // Change the abstract representation

   v_node_injection[ 0 ].set_value( v_active_demand[ 0 ] );
  }
 }

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification

  Block::add_Modification(
   std::make_shared< NetworkBlockSbstMod >( this,
                                            NetworkBlockMod::eSetActD,
                                            std::move( subset ) ),
   Observer::par2chnl( issuePMod ) );
 }
}

void BusNetworkBlock::set_active_demand(
 std::vector< double >::const_iterator values,
 Block::Range rng,
 c_ModParam issuePMod,
 c_ModParam issueAMod ) {

 // For BusNetworkBlock, this method degenerates a bit

 if( rng != Range( 0, 1 ) ) {
  throw ( std::invalid_argument( "invalid value in subset" ) );
 }

 if( v_active_demand.empty() ) {
  if( *values == 0 ) {
   return;
  }
  v_active_demand.assign( 1, 0 );
 }

 // If nothing changes, return
 if( v_active_demand[ 0 ] == *values ) {
  return;
 }

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  v_active_demand[ 0 ] = *values ;

  if( not_dry_run( issueAMod ) && AR & HasVar ) {
   // Change the abstract representation

   v_node_injection[ 0 ].set_value( v_active_demand[ 0 ] );
  }
 }

 if( issue_pmod( issuePMod ) ) {
  Block::add_Modification(
   std::make_shared< NetworkBlockRngdMod >( this,
                                            NetworkBlockMod::eSetActD,
                                            rng ),
   Observer::par2chnl( issuePMod ) );
 }
}

/*--------------------------------------------------------------------------*/
/*--------------------- End File BusNetworkBlock.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
