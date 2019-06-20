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
 * \date 20 - 06 - 2019
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
#include "ColVariable.h"
#include "FRowConstraint.h"

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
 * - A virtual public method that is used to initialize and read the
 *   data of any possible derived NetworkBlock class.
 *
 * - The number of nodes and lines of the network.
 *
 * - A vector of doubles used to store the values of the demand at
 *   each node of the network, which have size equal to the number of nodes
 *   or is empty, in which case the corresponding variables simply do
 *   not exist (for instance, the network may not have reserve).
 *
 * - A vector of doubles used to store the values of the susceptance
 *   of each line of the network, which has size equal to the number of lines
 *   or is empty, in which case the corresponding variables simply do
 *   not exist (for instance, the network may not have reserve).
 *
 * - Two vectors of doubles to store the minimum and maximum power
 *   flow in each line of the network, , which have size equal to the number
 *   of nodes or is empty, in which case the corresponding variables simply do
 *   not exist (for instance, the network may not have reserve).
 *
 * - A vector of ColVariable objects, that are used to store the power
 *   injection to each node. */

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

 NetworkBlock( Block * father_block = nullptr ): Block( father_block ) {}

/*--------------------------------------------------------------------------*/

 /// destructor of NetworkBlock

 virtual ~NetworkBlock();

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
 * - the dimension "NumberNodes" containing the number of nodes in the
 *   problem; this dimension is optional, if it is not provided then it is
 *   taken to be == 1;
 *
 * - the dimension "NumberLines" containing the number of arcs in the problem;
 *   if NumberNodes == 1 then this dimension need not to be present since it
 *   is not loaded;
 *
 * - the variable "StartLine", of type int and indexed over the dimension
 *   "NumberNodes"; the i-th entry of the variable is the starting point of
 *   the line (however, lines are not oriented)
 *
 * - the variable "EndLine", of type int and indexed over the dimension
 *   "NumberNodes"; the i-th entry of the variable is the ending point of the
 *   line (however, lines are not oriented)
 *
 * - the variable "ActiveDemand", of type double and indexed over the
 *   dimension "NumberNodes"; the i-th entry of the variable is assumed to
 *   contain the active power demand at node i in the network;
 *
 * - the variable "MinPowerFlow", of type double and indexed over the
 *   dimension "NumberLines"; the i-th entry of the variable is assumed to
 *   contain the minimum power flow at line i; if NumberNodes == 1 then this
 *   variable need not to be present since it is not loaded;
 *
 * - the variable "MaxPowerFlow", of type double and indexed over the
 *   dimension "NumberLines"; the i-th entry of the variable is assumed to
 *   contain the maximum power flow at line i; if NumberNodes == 1 then this
 *   variable need not to be present since it is not loaded;
 *
 * - the variable "Susceptance", of type double and indexed over the
 *   "NumberLines"; the i-th entry of this variable is assumed to contain the
 *   susceptance of line i; if NumberNodes == 1 then this variable need not to
 *   be present since it is not loaded.
 */

 virtual void deserialize( netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
/// generate the static variables of NetworkBlock
 /** Method that generates the static variables of this NetworkBlock. The base
  * NetworkBlock class has just the node injection variables;
  */

 virtual void generate_abstract_variables( Configuration *stvv = nullptr )
   override;

/*--------------------------------------------------------------------------*/


 virtual void load( std::istream &input ) override {
      throw( std::logic_error( "NetworkBlock::load() not implemented yet" ) );
 };

/*@} -----------------------------------------------------------------------*/
/*--------------- METHODS FOR MODIFYING THE NetworkBlock -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the NetworkBlock
 *  @{ */

 /// set the network data method
 /** Set the network data method.
  * This method can be called *before* that deserialize() is called to provide
  * the NetworkBlock with all data needed. This allows the information not to
  * be duplicated in the netCDF group that describes the network, since
  * usually (bit not necessarily) a NetworkBlock is deserialized inside a
  * UCBlock, and all networks have the same data, that can therefore be read
  * once and for all by the father UCBlock.
  *
  * If this method is *not* called, which means that the set_NetworkData() has
  * no object, then when deserialize() is called the information has to be
  * available by other means, i.e.:
  *
  * (i)  If there is no data for the NetworkBlock in netCDF input, then the
  *
  *      NetworkBlock must have a father, which must be a UCBlock: the
  *      network data is then taken to be that of the father. If the
  *      NetworkBlock does not have a father (or it is not a UCBlock), then
  *      exception is thrown.
  *
  * (ii) If all the data is presented in the netCDF input of NetworkBlock,
  *      the data provided there is used with no check that the NetworkBlock
  *      has a father at all, or the father is a UCBlock.
  *
  * If this method *is* called, which has to happen before that deserialize()
  * is called, then if the network data is present in netCDF input, then
  * without checking any thing with father, uses it. If the network data is
  * not presented in netCDF input, it must have been passed from father which
  * is UCBlock.
  *
  * If this method is called *after* that deserialize() is called, this is
  * taken to mean that the NetworkBlock is being "reset", and that immediately
  * after deserialize() will be called again. The same rules as above are to
  * be followed for that subsequent call to deserialize(). */

 void set_NetworkData( );

/*@} -----------------------------------------------------------------------*/
/*----------- METHODS FOR READING THE DATA OF THE NetworkBlock -------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the NetworkBlock
    @{ */

 /// returns the vector of node injection variables
 const std::vector<ColVariable> & get_node_injection() const {
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
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 /// number of nodes of the network
 Index f_number_nodes;

 /// number of lines of the network
 Index f_number_lines;

 /// set starting lines
  std::vector< int > v_start_line;

 /// set ending lines
 std::vector< int > v_end_line;

 /// vector to store the demand of each node of the network
 std::vector< double > v_active_demand;

 /// vector to store the susceptance of each line of the network
 std::vector< double > v_susceptance;

 /// vector to store the minimum power flow at each line
 std::vector< double > v_min_power_flow;

 /// vector to store the maximum power flow at each line
 std::vector< double > v_max_power_flow;

 /// power injection at each node
 std::vector< ColVariable > v_node_injection;

 /// flow limit constraints
// std::vector<FRowConstraint> v_flow_limit_constraints;

};   // end( class( NetworkBlock ) )

}  /* namespace SMSpp_di_unipi_it */

#endif /* NetworkBlock.h included */

/*--------------------------------------------------------------------------*/
/*------------------------ End File NetworkBlock.h -------------------------*/
/*--------------------------------------------------------------------------*/
