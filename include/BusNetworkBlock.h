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

#include "Block.h"
#include "NetworkBlock.h"
#include "ColVariable.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it {

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS BusNetworkBlock ---------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/

/// implementation of the Block concept for the BusNetworkBlock
/** The BusNetworkBlock class implements the Block concept [see Block.h] for
 *  a "reasonably standard" bus-network of a Unit Commitment Problem.
 *  A BusNetworkBlock defined when the NetworkBlock has just one node. The
 *  demand satisfaction should be constructed as follow:
 *
 * \f[
 *  S_{t} = D^{ac}_{t}) \quad t \in \mathcal{T}
 * \f]
 *   Where \f$ S_{t} \f$  and \f$ D^{ac}_{t} \f$ are the node injection
 *   variable of each node and the active demand in the network.
 */

  class BusNetworkBlock : public NetworkBlock {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

  public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
 *
 * BusNetworkBlock defines a main public type:
 *
 * - Index, the type of indices;
 @{ */

typedef std::size_t Index;

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

/// constructor of BusNetworkBlock, taking possibly a pointer to its father

 BusNetworkBlock( Block * f_block = nullptr ): NetworkBlock( f_block ) { }

/*--------------------------------------------------------------------------*/

/// destructor of BusNetworkBlock

 virtual ~BusNetworkBlock() {}

/*@} -----------------------------------------------------------------------*/
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
 * fixed by the active demand value for that node. */

 virtual void generate_abstract_variables( Configuration *stvv = nullptr )
    override;

/*@} -----------------------------------------------------------------------*/
/*--------------- METHODS FOR MODIFYING THE BusNetworkBlock ----------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the BusNetworkBlock
 *  @{ */

 /// method to set the NetworkData object
 /** Method to set the NetworkData object. This method does nothing in
  * BusNetworkBlock because by definition the network is made by only one
  * node. */
 
 void set_NetworkData( UCBlock::NetworkData * network_data = nullptr )
   override { }

/*@} -----------------------------------------------------------------------*/
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

/*@} -----------------------------------------------------------------------*/
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


 /// the demand of node
 double f_active_demand;

  };   // end( class( BusNetworkBlock ) )

}  /* namespace SMSpp_di_unipi_it */

#endif /* BusNetworkBlock.h included */

/*--------------------------------------------------------------------------*/
/*----------------------- End File BusNetworkBlock.h -----------------------*/
/*--------------------------------------------------------------------------*/
