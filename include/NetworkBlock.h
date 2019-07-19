/*--------------------------------------------------------------------------*/
/*--------------------------- File NetworkBlock.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 *
 * Header file for the class NetworkBlock, which derives from the Block, in
 * order to define the basic interface for the constraints/optimization
 * problems which describe the behaviour of the transmission network in a
 * specific time instant in the Unit Commitment (UC) problem, as represented
 * in UCBlock.
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
 * Copyright &copy; by Antonio Frangioni, Ali Ghezelsoflu, Rafael
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

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it {

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS NetworkBlock ----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// Block that describes the transmission network in the UC problem
/** The class NetworkBlock, which derives from the Block, defines the basic
 * interface for the constraints/optimization problems which describe the
 * behaviour of the transmission network in a specific time instant in the
 * Unit Commitment (UC) problem, as represented in UCBlock.
 *
 * The base class handles only basic information: it allows to read/set the
 * topology (and capacity/susceptances) of the network, and the active power
 * demand at the different nodes in the given time instant. This information
 * is actually bunched together in a small "passive" NetworkData object (no
 * methods, just a data repository) that can be either de-serialized or
 * passed ready-made (typically, by the UCBlock). Details of the kind of
 * network that is implemented ("bus", DC equations, AC equations, OPF, ...)
 * are entirely demanded to derived objects. The interface between a
 * NetworkBlock and the rest of the UC is just the vector of node injection
 * variables, which will have to satisfy the technical constraints of the
 * transmission network. */

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
 *
 * - NetworkData, a small auxiliary class to bunch together the basic data
 *   (topology and electrical characteristics) of the transmission network.
 *  @{ */

 typedef std::size_t Index;                 ///< index of parameters

/*--------------------------------------------------------------------------*/
/*-------------------- CLASS NetworkBlock::NetworkData ---------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// auxiliary class holding basic data about the transmission network
/** The NetworkData class is a nested sub-class which only serves to have a
 * quick way to load all the basic data (topology and electrical
 * caracteristics) that describe the transmission network. The rationale is
 * that while often the network does not change during the (short) time
 * horizon of UC, it makes sense to allow for this to happen. This means that
 * individual NetworkBlock objects may in principle have different
 * NetworkData, but most often they can share the same. By bunching all the
 * information together we make it easy for this sharing to happen. */

 class NetworkData {

/*--------------------------------------------------------------------------*/
/*----------------- PUBLIC PART OF THE NetworkData CLASS -------------------*/
/*--------------------------------------------------------------------------*/

  public:

/**@} ----------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

  /// constructor of NetworkData, does nothing
  NetworkData();

  /// destructor of NetworkData: it is virtual, and empty
  ~NetworkData() = default;

/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

  /// deserialize a NetworkData out of a netCDF::NcGroup
  /** Deserialize a NetworkData out of a netCDF::NcGroup, which should contain
   * the following:
   *
   * - The dimension "NumberNodes" containing the number of nodes in the
   *   problem; this dimension is optional, if it is not provided then it is
   *   taken to be == 1.
   *
   * If NumberNodes == 1 (equivalently, it is not provided), the network is a
   * "bus" formed of only one node, and therefore all the subsequent
   * information need not to be present since it is not loaded. If
   * NumberNodes > 1, then all the subsequent information is mandatory:
   *
   * - The dimension "NumberLines" containing the number of lines in the
   *   transmission network.
   *
   * - The variable "StartLine", of type int and indexed over the dimension
   *   "NumberNodes"; the i-th entry of the variable is the starting point of
   *   the line (a number in 0, ..., NumberNodes - 1). Note that lines are not
   *   oriented, but the flow of energy is; that is, a positive flow along
   *   line i means that energy is being taken away from StartLine[ i ] and
   *   delivered to EndLine[ i ] (see next), a negative flow means vice-versa.
   *   Note that node names here go from 0 to NNodes.getSize() - 1;
   *
   * - The variable "EndLine", of type int and indexed over the dimension
   *   "NumberNodes"; the i-th entry of the variable is the ending point of the
   *   line (a number in 0, ..., NumberNodes - 1; lines are not oriented, but
   *   see above). StartLine[ i ] == EndLine[ i ] (a self-loop) is not allowed,
   *   but multiple lines between the same pair of nodes are. Note that node
   *   names here go from 0 to NNodes.getSize() - 1;
   *
   * - The variable "MinPowerFlow", of type double and indexed over the
   *   dimension "NumberLines"; the i-th entry of the variable is assumed to
   *   contain the minimum power flow on line i (note that this is typically
   *   a negative number as lines are bi-directional, see above).
   *
   * - The variable "MaxPowerFlow", of type double and indexed over the
   *   dimension "NumberLines"; the i-th entry of the variable is assumed to
   *   contain the maximum power flow at line i (a non-negative number).
   *
   * - The variable "Susceptance", of type double and indexed over the dimension
   *   "NumberLines"; the i-th entry of this variable is assumed to contain the
   *   susceptance of line i. Note that this is strictly a positive value. */
  
  virtual void deserialize( netCDF::NcGroup & group );

/**@} ----------------------------------------------------------------------*/
/*------------- METHODS FOR READING THE DATA OF THE NetworkData ------------*/
/*--------------------------------------------------------------------------*/
  /// returns the number of nodes of the network
  /** Method for returning the number of nodes of the network. When it is
   * equal to one, it means that the transmission network is bus, and
   * therefore all the rest of the data is meaningless. */

  Index get_number_nodes() const { return( f_number_nodes ); }

/*--------------------------------------------------------------------------*/
  /// returns the number of lines of the network
  /** Method for returning the number of lines of the network. When
   * get_number_nodes() == 1 (the network is a bus), get_number_lines() == 0
   * (no self-loops are allowed, hence there is no line to be made with a
   * single node). */

  Index get_number_lines() const { return( f_number_lines ); }
 
/*--------------------------------------------------------------------------*/
  /// returns the vector of start nodes
  /** Method for returning the vector of starting point of each line. This
   *  vector may have empty size (bus network) or the size of number of nodes,
   *  then there are two possible cases:
   *
   *  - if f_number_nodes == 1, this vector has empty size which means there
   *    is no line at network (bus network), and this vector is not needed to
   *    be defined.
   *
   *  - if f_number_nodes > 1, this vector have size of f_number_nodes and each
   *    element of the vectors gives starting point of each line in the network.
   */

  const std::vector< Index > & get_start_line() const { return v_start_line; }

/*--------------------------------------------------------------------------*/
  /// returns vector of end nodes
  /** Method for returning the vector of ending point of each line. This
   *  vector may have empty size (bus network) or the size of number of nodes,
   *  then there are two possible cases:
   *
   *  - if f_number_nodes == 1, this vector has empty size which means there
   *    is no line at network (bus network), and this vector is not needed to
   *    be defined.
   *
   *  - if f_number_nodes > 1, this vector have size of f_number_nodes and each
   *    element of the vectors gives ending point of each line in the network.
   */

  const std::vector< Index > & get_end_line() const { return( v_end_line ); }

/*--------------------------------------------------------------------------*/
  /// returns vector of the minimum power flow
  /** Method for returning the vector of minimum power flow of each line. This
   *  vector may have empty size (bus network) or the size of number of nodes,
   *  then there are two possible cases:
   *
   *  - if f_number_lines == 0, this vector has empty size which means there is
   *    no line at network (bus network).
   *
   *  - if f_number_lines >= 1, this vector have size of f_number_lines and each
   *    element of the vectors gives minimum power flow of each line in the
   *    network. */

  const std::vector< double > & get_min_power_flow() const {
   return v_min_power_flow;
   }

/*--------------------------------------------------------------------------*/
  /// returns vector of the maximum power flow
  /** Method for returning the vector of maximum power flow of each line. This
   *  vector may have empty size (bus network) or the size of number of nodes,
   *  then there are two possible cases:
   *
   *  - if f_number_lines == 0, this vector has empty size which means there is
   *    no line at network (bus network).
   *
   *  - if f_number_lines >= 1, this vector have size of f_number_lines and
   *    each element of the vectors gives maximum power flow of each line in
   *    the network. */

  const std::vector< double > & get_max_power_flow() const {
   return v_max_power_flow;
   }

/*--------------------------------------------------------------------------*/
 /// returns vector of the susceptances
  /** Method for returning the vector of susceptances for each line. This
   *  vector may have empty size (bus network) or the size of number of nodes,
   *  then there are two possible cases:
   *
   *  - if f_number_lines == 0, this vector has empty size which means there is
   *    no line at network (bus network).
   *
   *  - if f_number_lines >= 1, this vector has size of f_number_lines and each
   *    element of the vectors gives the Susceptance value for each line in the
   *    network. */
  
  const std::vector< double > & get_susceptance() const {
   return v_susceptance;
   }

/**@} ----------------------------------------------------------------------*/
/*--------------------- METHODS FOR SAVING THE NetworkData -----------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the NetworkData
 *  @{ */

/// Serialize a NetworkData out of a netCDF::NcGroup
/** Serialize a NetworkData out of a netCDF::NcGroup to the specific format of
 * a NetworkData. See NetworkBlock::deserialize( netCDF::NcGroup ) for details
 * of the format of the created netCDF group. */

  virtual void serialize( netCDF::NcGroup & group ) const;

/**@} ----------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

  protected:

/*--------------------------------------------------------------------------*/
/*----------------- PROTECTED FIELDS OF THE NetworkData --------------------*/
/*--------------------------------------------------------------------------*/

  Index f_number_nodes;    ///< Number of nodes of the network

  Index f_number_lines;    ///< Number of lines of the network

  std::vector< Index > v_start_line;  ///< Vector of starting nodes

  std::vector< Index > v_end_line;    ///< Vector of ending nodes

  /// Vector to store the susceptance of each line of the network
  std::vector< double > v_susceptance;

  /// Vector to store the minimum power flow at each line
  std::vector< double > v_min_power_flow;

  /// Vector to store the maximum power flow at each line
  std::vector< double > v_max_power_flow;

/*--------------------------------------------------------------------------*/

 };   // end( class( NetworkData ) )

/**@} ----------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

/// Constructor, takes the father
/** Constructor of NetworkBlock, taking possibly a pointer of its father
 * Block. */

 explicit NetworkBlock( Block * father = nullptr ) : Block( father ) {}

/*--------------------------------------------------------------------------*/
/// Destructor of NetworkBlock

 ~NetworkBlock() override = default;

/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

/// Extends Block::deserialize( netCDF::NcGroup )
/** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
 * the NetworkBlock. Besides the mandatory "type" attribute of any :Block, the
 * group should contain the following:
 *
 * - Optionally, the dimensions and variables necessary to a NetworkData
 *   object, that describe the transmission network; see
 *   NetworkData::deserialize() for details. All that is optional, because the
 *   NetworkData object can alternatively be passed to the NetworkBlock via a
 *   call to set_NetworkData(). Note that if set_NetworkData() is called, but
 *   the representation of a NetworkData object is found in the NcGroup, then
 *   the NetworkData passed by set_NetworkData() is ignored, and a new
 *   NetworkData object is read from the NcGroup and used instead.
 *
 * - The "ActiveDemand", of type double, and of size "number of nodes". If the
 *   NetworkData object description is present in the NcGroup this is the
 *   dimension "NumberNodes", but this may not happen: if "NumberNodes" is
 *   not there then the NetworkData object must have been passed by
 *   set_NetworkData(), and the number of nodes can be read via
 *   NetworkData::get_number_nodes(). */

 void deserialize( netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
/// generate the static variables of NetworkBlock
/** Method that generates the static variables of this NetworkBlock. The
 * base NetworkBlock class has just the node injection variables, which are
 * mandatory as that's how the NetworkBlock is linked to the rest of the UC
 * model. The size of this variable is the number of nodes, which can be
 * read 
 *
 * - if NetworkData object is not provided (basically, "NumberNodes" is not
 *   provided or it is == 1) then the transmission network is taken to have
 *   only one node (a bus) and there is only one node injection variable.
 *
 * - if the NetworkData object is present (either in the NcGroup or because it
 *   has been passed and NumberNodes > 1) then this variable has size
 *   "NumberNodes", which can be read via NetworkData::get_number_nodes(). */

 void generate_abstract_variables( Configuration * stvv ) override {}

/*--------------------------------------------------------------------------*/
/**
 * @brief It loads a NetworkBlock from a input standard stream.
 * @warning This method is not implemented yet.
 * @param input an input stream
 */

 void load( std::istream & input ) override {
  throw ( std::logic_error( "NetworkBlock::load() not implemented yet" ) );
  }

/**@} ----------------------------------------------------------------------*/
/*--------------- METHODS FOR MODIFYING THE NetworkBlock -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the NetworkBlock
 *  @{ */

 /// method to set the NetworkData object
 /** This method can be called *before* that deserialize() is called to provide
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
  *
  * The default implementation of this method is empty, which is OK for a
  * NetworkBlock which only handles the "bus" case. */

 virtual void set_NetworkData( NetworkData * nd ) {}

/**@} ----------------------------------------------------------------------*/
/*----------- METHODS FOR READING THE DATA OF THE NetworkBlock -------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the NetworkBlock
    @{ */

 /// returns the number of nodes
 /** Returns the number of nodes in the transmission network. This should
  * just be equivalent to get_NetworkData()-> get_number_nodes(), but the
  * base NetworkBlock class does not handle it, and therefore it assumes the
  * network is a bus and returns 1. */

 virtual Index get_number_nodes( void ) const { return( 1 ); }

/*--------------------------------------------------------------------------*/
 /// returns the NetworkData object
 /** The method of the base class always returns nullptr, because the base
  * class does not handle the NetworkData object. This is OK for derived
  * classes that only handle the "bus" case. */

 virtual NetworkData * get_NetworkData() const {
  return nullptr;
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of active demands
 /** Method for returning the active demand for the given node, which is
  * assumed to have size get_number_nodes(). */

 const std::vector< double > & get_active_demand() const {
  return v_active_demand;
  }

/**@} ----------------------------------------------------------------------*/
/*----------- METHODS FOR READING THE Variable OF THE NetworkBlock ---------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the Variable of the NetworkBlock
 *
 * These methods allow to read the just the one set of Variable (which is node
 * injection ones for each node) that NetworkBlock in necessarily has.
 * @{ */

 /// returns the vector of node injection variables
 /** Method for returning vector of node injection variables, which is
  * assumed to have size get_number_nodes(). */

 const std::vector< ColVariable > & get_node_injection() const {
  return v_node_injection;
  }

/**@} ----------------------------------------------------------------------*/
/*--------------------- METHODS FOR SAVING THE NetworkBlock ----------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the NetworkBlock
 *  @{ */

/// Extends Block::serialize( netCDF::NcGroup )
/** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
 * NetworkBlock. See NetworkBlock::deserialize( netCDF::NcGroup ) for
 * details of the format of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/**@} ----------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 /// vector to store the demand of each node of the network
 std::vector< double > v_active_demand;

 /// power injection at each node
 std::vector< ColVariable > v_node_injection;

/*--------------------------------------------------------------------------*/

 };   // end( class( NetworkBlock ) )

/*--------------------------------------------------------------------------*/

}  /* namespace SMSpp_di_unipi_it */

/*--------------------------------------------------------------------------*/

#endif /* NetworkBlock.h included */

/*--------------------------------------------------------------------------*/
/*------------------------ End File NetworkBlock.h -------------------------*/
/*--------------------------------------------------------------------------*/
