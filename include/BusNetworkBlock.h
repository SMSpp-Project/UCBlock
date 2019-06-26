/*--------------------------------------------------------------------------*/
/*-------------------------- File BusNetworkBlock.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class BusNetworkBlock, which derives from NetworkBlock
 * [see NetworkBlock.h] in order to define a "bus" transmission network in
 * the Unit Commitment problem.
 *
 * \version 0.11
 *
 * \date 24 - 06 - 2019
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
 * Copyright &copy; by Antonio Frangioni, Ali Ghezelsoflu, Rafael
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

/// a NetworkBlock with only one node, i.e., a "bus" transmission network
/** The BusNetworkBlock class, which derives from NetworkBlock [see
 * NetworkBlock.h] implements the Block concept [see Block.h] in order to
 * define a "bus" transmission network in the Unit Commitment problem for a
 * given instant in the time horizon. A "bus" network has just one node,
 * and therefore a single value D for the demand and a single injection
 * variable s, which can hardly be called a variable since the only possible
 * way to satisfy the constraints is by having
 *
 *    s = D
 *
 * which in fact makes the variable a constant. */

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

 BusNetworkBlock( Block * f_block = nullptr ): NetworkBlock( f_block ) { }

/*--------------------------------------------------------------------------*/

/// destructor of BusNetworkBlock, (understandably) doea nothing

 virtual ~BusNetworkBlock() {}

/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

/// loads the BusNetworkBlock instance from memory
/** Loads the BusNetworkBlock instance from memory.
 * Like load( std::istream & ), if there is any Solver attached to this
 * BusNetworkBlock then a NBModification (the "nuclear option") is issued.
 * */

 virtual void load( std::istream &input ) override {
  throw( std::logic_error( "BusNetworkBlock::load() not implemented yet" ) );
  }

/*--------------------------------------------------------------------------*/

/// extends Block::deserialize( netCDF::NcGroup )
/** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
 * the BusNetworkBlock. Besides the mandatory "type" attribute of any :Block,
 * the group should contain the following:
 *
 * TODO: no, ActiveDemand is loaded by NetworkBlock already. I think we don't
 *       need to do anything here bacause the base NetworkBlock should have
 *       done it all (see the e-mail about get_NetworkData())
 *
 * - the variable "ActiveDemand", of type double and indexed over the
 *   dimension "NumberNodes" when NumberNodes == 1; the entry of the variable
 *   is assumed to contain the active power demand at node 1 in the network;
 *   the variable is optional and has to be defined when NumberNodes == 1; if
 *   NumberNodes > 1 there is no Bus-Network then this variable need not be
 *   defined, since is not loaded.
 */

 virtual void deserialize( netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
/// generate the static variables of BusNetworkBlock
/** Method that generates the static variables of this BusNetworkBlock. The
 * base BusNetworkBlock class has just the node injection variables. Since,
 * there exists just one node in this class, the variable node injection is
 * fixed to the active demand value for that node. */

 virtual void generate_abstract_variables( Configuration *stvv = nullptr )
    override;

/**@} ----------------------------------------------------------------------*/
/*--------------- METHODS FOR MODIFYING THE BusNetworkBlock ----------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the BusNetworkBlock
 *  @{ */

 // TODO: no
 /// method to set the NetworkData object
 /** Method to set the NetworkData object. This method does nothing in
  * BusNetworkBlock because by definition the network is made by only one
  * node. */
 
 void set_NetworkData( UCBlock::NetworkData * nd = nullptr ) override { }

/**@} ----------------------------------------------------------------------*/
/*-------------------- METHODS FOR SAVING THE BusNetworkBlock --------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the BusNetworkBlock
 *  @{ */

/// extends Block::serialize( netCDF::NcGroup )
/** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
 * BusNetworkBlock. See BusNetworkBlock::deserialize( netCDF::NcGroup ) for
 * details of the format of the created netCDF group.
 */

 virtual void serialize( netCDF::NcGroup & group ) const override;

/**@} ----------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

  protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 // TODO: no, NetworkBlock already has v_active_demand, why should you have
 //       this??
 /// the demand of node
 double f_active_demand;

  };   // end( class( BusNetworkBlock ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  /* end namespace SMSpp_di_unipi_it */

/*--------------------------------------------------------------------------*/

#endif /* BusNetworkBlock.h included */

/*--------------------------------------------------------------------------*/
/*----------------------- End File BusNetworkBlock.h -----------------------*/
/*--------------------------------------------------------------------------*/
