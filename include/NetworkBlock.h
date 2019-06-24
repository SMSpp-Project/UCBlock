/*--------------------------------------------------------------------------*/
/*--------------------------- File NetworkBlock.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 *
 * Header file for the class NetworkBlock, which derives from the Block, in
 * order to define a vary basic interface for any possible network that can be
 * attached to a UCBlock. It has very basic information that can characterize
 * almost any different kind of network(such as BusNetworkBlock,
 * DCNetworkBlock, and ACNetworkBlock , ..), which includes a set of node
 * injection variables.
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
 * Copyright &copy by Antonio Frangioni, Ali Ghezelsoflu, Rafael
 * Durbano Lobato, and Kostas Tavlaridis-Gyparakis
 */

/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __NetworkBlock
#define __NetworkBlock
/* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Block.h"
#include "UCBlock.h"
#include "ColVariable.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it {

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS NetworkBlock ----------------------------*/
/*--------------------------------------------------------------------------*/
/** The class NetworkBlock, which derives from the Block, defines a base class
 * for any possible network that can be attached to a UCBlock. It has very
 * basic information that can characterize almost any different kind of
 * network(such as BusNetworkBlock, DCNetworkBlock, and ACNetworkBlock , ..),
 * which includes a set of node injection variables. This class has thus been
 * constructed having the following elements:
 *
 * - A virtual public method that is used to initialize and read the data of
 *   any possible derived NetworkBlock class.
 *
 * - A vector of doubles used to store the values of the demand at each node
 *   of the network, which have size equal to the number of nodes or is empty,
 *   in which case the corresponding variables simply do not exist.
 *
 * - A vector of ColVariable objects, that are used to store the information
 *   regarding the node injection variable for each node. */

class NetworkBlock : public Block {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:
/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
 *
 * NetworkBlock defines a main public type:
 *
 * - Index, the type of indices;
 @{ */

typedef std::size_t Index;                 ///< index of parameters

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */
/// constructor, takes the father and the network-class data
/** Constructor of NetworkBlock, taking possibly a pointer of its father
 * Block and the network class data. */

 NetworkBlock( Block * father = nullptr ) : Block( father ) { }

/*--------------------------------------------------------------------------*/
 /// destructor of NetworkBlock

 virtual ~NetworkBlock() {}

/*@} -----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

/// extends Block::deserialize( netCDF::NcGroup )
/** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
 * the NetworkBlock. Besides the mandatory "type" attribute of any :Block, the
 * group should contain the following:
 *
 * - the variable "ActiveDemand", of type double and indexed over the
 *   dimension "NumberNodes"; the i-th entry of the variable is assumed to
 *   contain the active power demand at node i in the network;
 *
 * Furthermore, there may be the dimensions and variables necessary to a
 * NetworkData object, that describe the transmission network. See
 * NetworkData::deserialize() for details. All that is optional, because the
 * NetworkData object can alternatively be passed to the NetworkBlock via a
 * call to set_NetworkData(). Note that if set_NetworkData() is called, but
 * the representation of a NetworkData object is found in the NcGroup, then
 * the NetworkData passed by set_NetworkData() is ignored, and a new
 * NetworkData object is read from the NcGroup and used instead. */

 virtual void deserialize( netCDF::NcGroup & group ) override;
/*--------------------------------------------------------------------------*/
/// generate the static variables of NetworkBlock
/** Method that generates the static variables of this NetworkBlock. The
 * base NetworkBlock class has just the node injection variables. */

 virtual void generate_abstract_variables( Configuration *stvv = nullptr )
   override;
/*--------------------------------------------------------------------------*/

 virtual void load( std::istream &input ) override {
   throw( std::logic_error( "NetworkBlock::load() not implemented yet" ) );
 }

/*@} -----------------------------------------------------------------------*/
/*--------------- METHODS FOR MODIFYING THE NetworkBlock -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the NetworkBlock
 *  @{ */

 /// method to set the NetworkData object
 /** Method to set the NetworkData object.
  * This method can be called *before* that deserialize() is called to provide
  * the NetworkBlock with the data corresponding to the transmission network
  * description. This allows the information not to be duplicated in the
  * netCDF group that describes the NetworkBlock, since usually (but not
  * necessarily) a NetworkBlock is deserialized inside a UCBlock, and all
  * networks have the same data, that can therefore be read once and for all
  * by the father UCBlock.
  *
  * If this method is *not* called, which means that no NetworkData has been
  * provided, then when deserialize() is called the information has to be
  * available by other means, i.e.:
  *
  * (i)  If there is no data for the NetworkData in netCDF input, then the
  *      NetworkBlock must have a father, which must be a UCBlock: the
  *      network data is then taken to be that of the father. If the
  *      NetworkBlock does not have a father (or it is not a UCBlock), then
  *      exception is thrown.
  *
  * (ii) If all the data for the NetworkData is present in the netCDF input
  *      of NetworkBlock, the data provided there is used with no check that
  *      the NetworkBlock has a father at all, or the father is a UCBlock.
  *
  * If this method *is* called, which has to happen before that deserialize()
  * is called, then if the data for the NetworkData is present in netCDF
  * input, then it is used by the NetworkBlock, disregarding the NetworkData
  * object that was passed with this method. If the data for the NetworkData
  * is not present in netCDF input, it must have been passed from outside
  * with this method.
  *
  * If this method is called *after* that deserialize() is called, this is
  * taken to mean that the NetworkBlock is being "reset", and that immediately
  * after deserialize() will be called again. The same rules as above are to
  * be followed for that subsequent call to deserialize().
  *
  * Note that passing a new NetworkData causes all references to any previous
  * NetworkData to be lost. If the NetworkData was an "externally provided"
  * one this is no problem, but it means that it is responsibility of who
  * set it in the first place to delete it. If the NetworkData was created
  * by the NetworkBlock, it is the NetworkBlock's responsibility to delete
  * it during this call. This is not done in the base NetworkBlock class
  * because it has no data structures to hold the NetworkData pointer (in
  * fact, this method is pure virtual), so it is demanded to derived classes.
  */

 virtual void set_NetworkData
 ( UCBlock::NetworkData * network_data = nullptr ) = 0;

/*@} -----------------------------------------------------------------------*/
/*----------- METHODS FOR READING THE DATA OF THE NetworkBlock -------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the NetworkBlock
    @{ */

 /// returns the vector of node injection variables
 const std::vector<ColVariable> & get_node_injection( void ) const {
   return v_node_injection;
 }

/*@} -----------------------------------------------------------------------*/
/*--------------------- METHODS FOR SAVING THE NetworkBlock ----------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the NetworkBlock
 *  @{ */

/// extends Block::serialize( netCDF::NcGroup )
/** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
 * NetworkBlock. See NetworkBlock::deserialize( netCDF::NcGroup ) for
 * details of the format of the created netCDF group.
 */

 virtual void serialize( netCDF::NcGroup & group ) const override;

/*@} -----------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

protected:

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/
 /// the NetworkData object
 UCBlock::NetworkData * f_NetworkData;

 /// vector to store the demand of each node of the network
 std::vector< double > v_active_demand;

 /// power injection at each node
 std::vector< ColVariable > v_node_injection;

 };   // end( class( NetworkBlock ) )

/*--------------------------------------------------------------------------*/

}  /* namespace SMSpp_di_unipi_it */

#endif /* NetworkBlock.h included */

/*--------------------------------------------------------------------------*/
/*------------------------ End File NetworkBlock.h -------------------------*/
/*--------------------------------------------------------------------------*/
