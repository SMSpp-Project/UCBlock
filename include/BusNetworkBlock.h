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
 * given instant in the time horizon. */

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

/// destructor of BusNetworkBlock, (understandably) does nothing

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
 * When the NetworkData object is not provided by NetworkBlock (basically,
 * "NumberNodes" is not provided or it is == 1) then the transmission
 * network is taken to have only one node (a bus) and variable "ActiveDemand"
 * is long 1 which is provided by NetworkBlock *get_active_demand(). This part
 * does not do anything since the base NetworkBlock should have done it all
 */

 virtual void deserialize( netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
/// generate the static variables of BusNetworkBlock
/** Method that generates the static variables of this BusNetworkBlock. The
 * base BusNetworkBlock class has just the node injection variables. Since, a
 * "bus" network has just one node, and therefore a single value D for the
 * demand and a single injection variable s, which can hardly be called a
 * variable since the only possible way to satisfy the constraints is by
 * having s = D which in fact makes the variable a constant.*/

 virtual void generate_abstract_variables( Configuration *stvv = nullptr )
    override;

/**@} ----------------------------------------------------------------------*/
/*--------------- METHODS FOR MODIFYING THE BusNetworkBlock ----------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the BusNetworkBlock
 *  @{ */

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

  };   // end( class( BusNetworkBlock ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  /* end namespace SMSpp_di_unipi_it */

/*--------------------------------------------------------------------------*/

#endif /* BusNetworkBlock.h included */

/*--------------------------------------------------------------------------*/
/*----------------------- End File BusNetworkBlock.h -----------------------*/
/*--------------------------------------------------------------------------*/
