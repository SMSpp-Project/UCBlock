/*--------------------------------------------------------------------------*/
/*-------------------------- File BusNetworkBlock.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class BusNetworkBlock, which derives from NetworkBlock
 * [see NetworkBlock.h] in order to define a "bus" transmission network in
 * the Unit Commitment problem. There is actually precious little that this
 * class has to do that is not done already by the base NetworkBlock class.
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
 * \copyright &copy; by Antonio Frangioni, Ali Ghezelsoflu, Rafael
 * Durbano Lobato, and Kostas Tavlaridis-Gyparakis
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __BusNetworkBlock
#define __BusNetworkBlock
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "NetworkBlock.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it {

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS BusNetworkBlock ---------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// A NetworkBlock with only one node, i.e., a "bus" transmission network
/** The BusNetworkBlock class, which derives from NetworkBlock [see
 * NetworkBlock.h] implements the Block concept [see Block.h] in order to
 * define a "bus" transmission network in the Unit Commitment problem for a
 * given instant in the time horizon. There is actually precious little that
 * this class has to do that is not done already by the base NetworkBlock
 * class.*/

class BusNetworkBlock : public NetworkBlock {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

 /// constructor of BusNetworkBlock, taking possibly a pointer to its father

 explicit BusNetworkBlock( Block * f_block = nullptr ) :
  NetworkBlock( f_block ) {
  v_node_injection.resize(1);
 }

/*--------------------------------------------------------------------------*/
 /// destructor of BusNetworkBlock, (understandably) does nothing

 ~BusNetworkBlock() override = default;

/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

/// loads the BusNetworkBlock instance from memory
/** Like load( std::istream & ), if there is any Solver attached to this
 *  BusNetworkBlock then a NBModification (the "nuclear option") is issued.
 */

 void load( std::istream & input ) override {
  throw ( std::logic_error( "BusNetworkBlock::load() not implemented yet" ) );
  }

/*--------------------------------------------------------------------------*/
/// generates the static variables of BusNetworkBlock
/** The base BusNetworkBlock class has just the node injection variables.
 * Since a "bus" network has just one node, and therefore a single value D for
 * the demand and a single injection variable s, which can hardly be called a
 * variable since the only possible way to satisfy the constraints is by
 * having s = D which in fact makes the variable a constant.*/

 void generate_abstract_variables( Configuration * stvv ) override;

/**@} ----------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

 void set_active_demand( std::vector< double >::const_iterator values,
                         Subset && subset = { 0 },
                         bool ordered = false,
                         c_ModParam issuePMod = eNoBlck,
                         c_ModParam issueAMod = eNoBlck ) final;

 void set_active_demand( std::vector< double >::const_iterator values,
                         Range rng = Range( 0, 1 ),
                         c_ModParam issuePMod = eNoBlck,
                         c_ModParam issueAMod = eNoBlck ) final;

 static void static_initialization() {
  register_method< BusNetworkBlock >( "BusNetworkBlock::set_active_demand",
                                      &BusNetworkBlock::set_active_demand,
                                      MS_dbl_sbst::args() );

  register_method< BusNetworkBlock >( "BusNetworkBlock::set_active_demand",
                                      &BusNetworkBlock::set_active_demand,
                                      MS_dbl_rngd::args() );
 }

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------data--------------------------------------*/

/*-----------------------------variables------------------------------------*/

/*----------------------------constraints-----------------------------------*/


/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
 private:

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE FIELDS -------------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/

};  // end( class( BusNetworkBlock ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* BusNetworkBlock.h included */

/*--------------------------------------------------------------------------*/
/*----------------------- End File BusNetworkBlock.h -----------------------*/
/*--------------------------------------------------------------------------*/
