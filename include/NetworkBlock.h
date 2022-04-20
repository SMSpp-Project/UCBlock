/*--------------------------------------------------------------------------*/
/*--------------------------- File NetworkBlock.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 *
 * Header file for the class NetworkBlock, which derives from the Block, in
 * order to define the basic interface for the constraints/optimization
 * problems which describe the behaviour of the network in a specific time
 * instant in the Unit Commitment (UC) problem, as represented in UCBlock.
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
 * Copyright &copy; by Antonio Frangioni, Ali Ghezelsoflu,
 *                  Rafael Durbano Lobato
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

namespace SMSpp_di_unipi_it
{
/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS NetworkBlock ----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// Block that describes the network in the UC problem
/** The class NetworkBlock, which derives from the Block, defines the basic
 * interface for the constraints/optimization problems which describe the
 * behaviour of the network in a specific time instant in the Unit Commitment
 * (UC) problem, as represented in UCBlock.
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
 * network. */

class NetworkBlock : public Block
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
 *
 * NetworkBlock defines the NetworkData public type, a small auxiliary class
 * to bunch the basic topological data of the network.
 * @{ */

/*--------------------------------------------------------------------------*/
/*-------------------- CLASS NetworkBlock::NetworkData ---------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/

 /// Auxiliary class holding basic data about the network
 /** The NetworkData class is a nested sub-class which only serves to have a
  * quick way to load all the basic data (topology and electrical
  * characteristics) that describe the network. The rationale is that while
  * often the network does not change during the (short) time horizon of UC,
  * it makes sense to allow for this to happen. This means that individual
  * NetworkBlock objects may in principle have different NetworkData, but
  * most often they can share the same. By bunching all the information
  * together we make it easy for this sharing to happen. */

 class NetworkData
 {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

  public:

/**@} ----------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 * @{ */

  /// constructor of NetworkData, does nothing
  NetworkData();

/*--------------------------------------------------------------------------*/

  /// construct a :NetworkData of specific type using the Block factory
  /** Use the NetworkData factory to construct a :NetworkData object of type
   * specified by classname (a std::string with the name of the class inside).
   * If there is no class with the given name, exception is thrown.
   *
   * Note that the method is static because the factory is static, hence it is
   * to be called as
   *
   *  NetworkData * myNetworkData = NetworkData::new_NetworkData( some_class );
   *
   * i.e., without any reference to any specific NetworkData (and, therefore,
   * it can be used to construct the very first NetworkData if needed).
   *
   * Note that the :NetworkData returned my this method is "empty": it
   * contains no (instance) data, and therefore it has to be explicitly
   * initialized with
   * any of the corresponding methods (operator>>, serialize(), anything that
   * the specific :NetworkData class provides) before it can be used.
   *
   * For this to work, each :NetworkData has to:
   *
   * - add the line
   *
   *       SMSpp_insert_in_factory_h;
   *
   *   to its definition (typically, in the private part in its .h file);
   *
   * - add the line
   *
   *       SMSpp_insert_in_factory_cpp_1( name_of_the_class );
   *
   *   to exactly *one* .cpp file, typically that :NetworkData .cpp file. If
   *   the name of the class contains any parentheses, then one must enclose
   *   the name of the class in parentheses and instead add the line
   *
   *       SMSpp_insert_in_factory_cpp_1( ( name_of_the_class ) );
   *
   * Any whitespaces that the given \p classname may contain is ignored. So,
   * for example, to create an instance of the class MyNetworkData<int> one
   * could pass "MyNetworkData<int>" or "MyNetworkData< int >"
   * (even " M y B l o c k < int > " would work).
   *
   * @param classname The name of the :NetworkData class that must be
   *                  constructed. */

  static NetworkData * new_NetworkData( const std::string & classname ) {
   const std::string classname_( SMSpp_classname_normalise(
    std::string( classname ) ) );
   const auto it = NetworkData::f_factory().find( classname_ );
   if( it == NetworkData::f_factory().end() )
    throw( std::invalid_argument( classname +
                                   " not present in NetworkData factory" ) );
   return( ( it->second )( nullptr ) );
  }

  /// destructor of NetworkData: it is virtual, and empty
  virtual ~NetworkData() = default;

/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 * @{ */

  /// deserialize a NetworkData out of a netCDF::NcGroup
  /** Deserialize a NetworkData out of a netCDF::NcGroup, which should contain
   * the following:
   *
   * - The dimension "NumberNodes" containing the number of nodes in the
   *   problem; this dimension is optional, if it is not provided then it is
   *   taken to be == 1.
   *
   * If NumberNodes == 1 (equivalently, it is not provided), the network is a
   * "bus" formed of only one node, and therefore all the subsequent information
   * need not to be present since it is not loaded. If NumberNodes > 1, then all
   * the subsequent information is mandatory:
   *
   * - The dimension "NumberLines" containing the number of lines in the
   *   network.
   *
   * - The variable "StartLine", of type netCDF::NcUint and indexed over the
   *   dimension "NumberLines"; the l-th entry of the variable is the starting
   *   point of the line (a number in 0, ..., NumberLines - 1). Note that lines
   *   are not oriented, but the flow of energy is; that is, a positive flow
   *   along line l means that energy is being taken away from StartLine[ l ]
   *   and delivered to EndLine[ l ] (see next), a negative flow means
   *   vice-versa. Note that node names here go from 0 to NNodes.getSize() - 1;
   *
   * - The variable "EndLine", of type netCDF::NcUint and indexed over the
   *   dimension "NumberLines"; the l-th entry of the variable is the ending
   *   point of the line (a number in 0, ..., NumberLines - 1; lines are not
   *   oriented, but see above). StartLine[ l ] == EndLine[ l ] (a self-loop) is
   *   not allowed, but multiple lines between the same pair of nodes are. Note
   *   that node names here go from 0 to NNodes.getSize() - 1;
   *
   * - The variable "MinPowerFlow", of type netCDF::NcDouble and indexed over
   *   the dimension "NumberLines". This is meant to represent the vector MxP[ l
   *   ] that, for each line l, contains the minimum power flow at line l (note
   *   that this is typically a negative number as lines are bi-directional, see
   *   above).
   *
   * - The variable "MaxPowerFlow", of type netCDF::NcDouble and indexed over
   *   the dimension "NumberLines". This is meant to represent the vector MxP[ l
   *   ] that, for each line l, contains the maximum power flow at line l (a
   *   non-negative number).
   *
   * - The variable "Susceptance", of type netCDF::NcDouble and indexed over the
   *   dimension "NumberLines". This is meant to represent the vector S[ l ]
   *   that, for each line i contains the susceptance of the network for the
   *   corresponding line i. Note that this variable is optional, for each line
   *   l if it is provided then it is assumed that S[ l ] != 0, otherwise it is
   *   assumed that S[ l ] == 0. In fact, when S[ l ] != 0 this corresponds to a
   *   model with AC lines, and when for each line l, it's not defined or S[ l ]
   *   == 0, then it corresponds to a single connected grid composed of HVDC
   *   lines only which is also known as the Net Transfer Capacity (NTC)
   *   model.
   *
   * - The variable "NetworkCost", of type netCDF::NcDouble and indexed over the
   *   dimension "NumberLines". This is meant to represent the vector NC[ l ]
   *   that, for each line l, contains the monetary cost to exchanges between
   *   nodes or each network. This variable is optional; if it is not provided
   *   then it's taken to be zero. */

  virtual void deserialize( const netCDF::NcGroup & group );

/**@} ----------------------------------------------------------------------*/
/*------------- METHODS FOR READING THE DATA OF THE NetworkData ------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the NetworkData
 * @{ */

  /// returns the number of nodes of the network
  /** Method for returning the number of nodes of the network. When it is equal
   * to one, it means that the network is bus, and therefore all the rest of
   * the data is meaningless. */

  Index get_number_nodes( void ) const { return( f_number_nodes ); }

/*--------------------------------------------------------------------------*/

  /// returns the number of intervals of the network
  /** Method for returning the number of intervals the network refers to. */

  Index get_number_intervals( void ) const { return( f_number_intervals ); }

/*--------------------------------------------------------------------------*/

  /// returns the number of lines of the network
  /** Method for returning the number of lines of the network. When
   * get_number_nodes() == 1 (the network is a bus), get_number_lines() == 0
   * (no self-loops are allowed, hence there is no line to be made with a single
   * node). */

  Index get_number_lines( void ) const { return( f_number_lines ); }

/*--------------------------------------------------------------------------*/

  /// returns the vector of start lines
  /** Method for returning the vector of starting point of each line. This
   *  vector may have empty size (bus network) or the size of number of lines,
   *  then there are two possible cases:
   *
   *  - if f_number_nodes == 1, this vector has empty size which means there is
   *    no line at network (bus network), and this vector is not needed to be
   *    defined.
   *
   *  - if f_number_nodes > 1, this vector have size of f_number_lines and each
   *    element of the vectors gives starting point of each line in the network.
  */

  const std::vector< Index > & get_start_line( void ) const {
   return v_start_line;
  }

/*--------------------------------------------------------------------------*/

  /// returns vector of end lines
  /** Method for returning the vector of ending point of each line. This vector
   * may have empty size (bus network) or the size of number of lines, then
   * there are two possible cases:
   *
   *  - if f_number_nodes == 1, this vector has empty size which means there is
   *    no line at network (bus network), and this vector is not needed to be
   *    defined.
   *
   *  - if f_number_nodes > 1, this vector have size of f_number_lines and each
   *    element of the vectors gives ending point of each line in the network.
   */

  const std::vector< Index > & get_end_line( void ) const {
   return( v_end_line );
  }

/**@} ----------------------------------------------------------------------*/
/*--------------------- METHODS FOR SAVING THE NetworkData -----------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the NetworkData
 * @{ */

  /// Serialize a NetworkData out of a netCDF::NcGroup
  /** Serialize a NetworkData out of a netCDF::NcGroup to the specific format of
   * a NetworkData. See NetworkBlock::deserialize( netCDF::NcGroup ) for details
   * of the format of the created netCDF group. */

  virtual void serialize( netCDF::NcGroup & group ) const;

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

  protected:

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED TYPES ------------------------------*/
/*--------------------------------------------------------------------------*/

  typedef boost::function< NetworkData *( NetworkData * ) > NetworkDataFactory;
  // type of the factory of NetworkData

  typedef std::map< std::string , NetworkDataFactory > NetworkDataFactoryMap;
  // Type of the map between strings and the factory of NetworkData

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED METHODS OF THE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

  /** @name Protected methods for handling static fields
   *
   * These methods allow derived classes to partake into static initialization
   * procedures performed once and for all at the start of the program. These
   * are typically related with factories.
   * @{ */

  /// method encapsulating the NetworkData factory
  /** This method returns the NetworkData factory, which is a static object.
   * The rationale for using a method is that this is the "Construct On
   * First Use Idiom" that solves the "static initialization order problem".
   */

  static NetworkDataFactoryMap & f_factory( void );

/*--------------------------------------------------------------------------*/

  /// empty placeholder for class-specific static initialization
  /** The method static_initialization() is an empty placeholder which is made
   * available to derived classes that need to perform some class-specific
   * static initialization besides these of any :NetworkBlock::NetworkData
   * class, i.e., the management of the factory. This method is invoked by the
   * SMSpp_insert_in_factory_cpp_* macros [see SMSTypedefs.h] during the
   * standard initialization procedures. If a derived class needs to perform
   * any static initialization it just have to do this into its version of
   * this method; if not it just has nothing to do, as the (empty) method of
   * the base class will be called.
   *
   * This mechanism has a potential drawback in that a redefined
   * static_initialization() may be called multiple times. Assume that a
   * derived class X redefines the method to perform something, and that a
   * further class Y is derived from X that has to do nothing, and that
   * therefore will not define Y::static_initialization(): them, within the
   * SMSpp_insert_in_factory_cpp_* of Y, X::static_initialization() will be
   * called again.
   *
   * If this is undesirable, X will have to explicitly instruct derived classes
   * to redefine their (empty) static_initialization(). Alternatively,
   * X::static_initialization() may contain mechanisms to ensure that it will
   * actually do things only the very first time it is called. One standard
   * trick is to do everything within the initialisation of a static local
   * variable of X::static_initialization(): this is guaranteed by the
   * compiler to happen only once, regardless of how many times the function
   * is called. Alternatively, an explicit static boolean could be used (this
   * may just be the same as what the compiler does during the initialization
   * of static variables without telling you). */

  static void static_initialization( void ) {}

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

  /// Number of nodes of the network
  Index f_number_nodes;

  /// Number of lines of the network
  Index f_number_lines;

  /// Number of intervals
  Index f_number_intervals;

  /// Vector of starting lines
  std::vector< Index > v_start_line;

  /// Vector of ending lines
  std::vector< Index > v_end_line;

  private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

  // Definition of State::private_name() (pure virtual)

  virtual const std::string & private_name( void ) const = 0;

 };  // end( class( NetworkData ) )

/**@} ----------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 * @{ */

 /// Constructor, takes the father
 /** Constructor of NetworkBlock, taking possibly a pointer of its father
  * Block. */

 explicit NetworkBlock( Block * father = nullptr ) : Block( father ) {}

/*--------------------------------------------------------------------------*/

 /// Destructor of NetworkBlock

 virtual ~NetworkBlock() override = default;

/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 * @{ */

 /// extends Block::deserialize( netCDF::NcGroup )
 /** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
  * the NetworkBlock. Besides the mandatory "type" attribute of any :Block, the
  * group should contain the following:
  *
  * - Optionally, the dimensions and variables necessary to a NetworkData
  *   object, that describe the network; see NetworkData::deserialize() for
  *   details. All that is optional, because the NetworkData object can
  *   alternatively be passed to the NetworkBlock via a call to
  *   set_NetworkData(). Note that if set_NetworkData() is called, but
  *   the representation of a NetworkData object is found in the NcGroup, then
  *   the NetworkData passed by set_NetworkData() is ignored, and a new
  *   NetworkData object is read from the NcGroup and used instead.
  *
  * - The "ActiveDemand", of type double, and of size "number of nodes". If the
  *   NetworkData object description is present in the NcGroup this is the
  *   dimension "NumberNodes", but the NetworkData object is optional and it
  *   may not be there. Thus, if "NumberNodes" is not there and "ActiveDemand"
  *   is, then the NetworkData object must have been passed by set_NetworkData(),
  *   and the number of nodes can be read via NetworkData::get_number_nodes().
  *   However, "ActiveDemand" itself is optional. If it is not found in the
  *   NcGroup, then it *must* be passed (either before or after the call to
  *   deserialize()) by calling set_active_demand(). Since both groups of data
  *   are optional, the NcGroup  can actually be empty which implies that all
  *   the data will be (or have been) passed by the in-memory interface. In
  *   this case, it would clearly be preferable to *entirely avoid the NcGroup
  *   to be there*, and in fact UCBlock has provisions for the NcGroup
  *   describing the NetworkBlock to be optional [see the comments to
  *   UCBlock::deserialize()]. */

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/

 /// generate the static variables of NetworkBlock
 /** Method that generates the static variables of this NetworkBlock. The
  * base NetworkBlock class has just the node injection variables, which are
  * mandatory as that's how the NetworkBlock is linked to the rest of the UC
  * model. The size of this variable is the number of nodes, which can be
  * read
  *
  * - if NetworkData object is not provided (basically, "NumberNodes" is not
  *   provided or it is == 1) then the network is taken to have only one node
  *   (a bus) and there is only one node injection variable.
  *
  * - if the NetworkData object is present (either in the NcGroup or because it
  *   has been passed and NumberNodes > 1) then this variable has size
  *   "NumberNodes", which can be read via NetworkData::get_number_nodes(). */

 void generate_abstract_variables( Configuration * stvv = nullptr )
 override {}

/*--------------------------------------------------------------------------*/

 /**
  * Loads a NetworkBlock from a input standard stream.
  * @warning This method is not implemented yet.
  * @param input an input stream
  * @param frmt the verbosity level
  */

 void load( std::istream & input , char frmt = 0 ) override {
  throw( std::logic_error( "NetworkBlock::load() not implemented yet" ) );
 }

/**@} ----------------------------------------------------------------------*/
/*--------------- METHODS FOR MODIFYING THE NetworkBlock -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the NetworkBlock
 * @{ */

 /// method to set the NetworkData object
 /** This method can be called *before* that deserialize() is called to provide
  * the NetworkBlock with the data corresponding to the network description.
  * This allows the information not to be duplicated in the netCDF group that
  * describes the NetworkBlock, since usually (but not necessarily) a
  * NetworkBlock is deserialized inside a UCBlock, and all networks have the
  * same data, that can therefore be read once and for all by the father
  * UCBlock.
  *
  * If this method is *not* called, which means that no NetworkData has been
  * provided, then when deserialize() is called the information has to be
  * available by other means, i.e.:
  *
  * (i)  If there is no data for the NetworkData in netCDF input, then the
  *      NetworkBlock must have a father, which must be a UCBlock: the network
  *      data is then taken to be that of the father. If the NetworkBlock does
  *      not have a father (or it is not a UCBlock), then exception is thrown.
  *
  * (ii) If all the data for the NetworkData is present in the netCDF input of
  *      NetworkBlock, the data provided there is used with no check that the
  *      NetworkBlock has a father at all, or the father is a UCBlock.
  *
  * If this method *is* called, which has to happen before that deserialize()
  * is called, then if the data for the NetworkData is present in netCDF input,
  * then it is used by the NetworkBlock, disregarding the NetworkData object
  * that was passed with this method. If the data for the NetworkData is not
  * present in netCDF input, it must have been passed from outside with this
  * method.
  *
  * If this method is called *after* that deserialize() is called, this is
  * taken to mean that the NetworkBlock is being "reset", and that immediately
  * after deserialize() will be called again. The same rules as above are to be
  * followed for that subsequent call to deserialize().
  *
  * Note that passing a new NetworkData causes all references to any previous
  * NetworkData to be lost. If the NetworkData was an "externally provided" one
  * this is no problem, but it means that it is responsibility of who set it in
  * the first place to delete it. If the NetworkData was created by the
  * NetworkBlock, it is the NetworkBlock's responsibility to delete it during
  * this call. This is not done in the base NetworkBlock class because it has
  * no data structures to hold the NetworkData pointer (in fact, this method is
  * pure virtual), so it is demanded to derived classes.
  *
  * The default implementation of this method is empty, which is OK for a
  * NetworkBlock which only handles the "bus" case. */

 virtual void set_NetworkData( NetworkData * nd ) {}

/*--------------------------------------------------------------------------*/

 /// method to set the ActiveDemand
 /** This method can be called either before or after that deserialize() is
  * called to provide the NetworkBlock with the ActiveDemand data. This allows
  * all Active Power Demand data corresponding to some UC problem to be
  * "grouped" together (typically, in UCBlock) rather than "spread" among the
  * different NetworkBlock, which may be convenient for some user.
  *
  * If this method is called *before* deserialize(), the data is just copied.
  * However, when deserialize() is called, if ActiveDemand data is present in
  * the NcGroup then this data is used, replacing (and therefore ignoring) the
  * data set by this method.
  *
  * Similarly, if this method is called *after* deserialize(), but some the
  * ActiveDemand was already present in the NcGroup, then that data is kept and
  * the call to this method does nothing.
  *
  * When this method is called, if it is empty it is written into, otherwise
  * nothing happens. In deserialize(), if the data is there in the NcGroup then
  * it is written in v_active_demand (which therefore is no longer empty),
  * otherwise it is left empty so that it can be set by this method. */

 virtual void set_ActiveDemand( const double * v ) = 0;

/**@} ----------------------------------------------------------------------*/
/*----------- METHODS FOR READING THE DATA OF THE NetworkBlock -------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the NetworkBlock
 * @{ */

 /// returns the number of nodes of the network
 /** Returns the number of nodes in the network. This should just be
  * equivalent to get_NetworkData()->get_number_nodes(), but the base
  * NetworkBlock class does not handle it, and therefore it assumes the network
  * is a bus, i.e., get_number_nodes() == 1, and returns 1. */

 virtual Index get_number_nodes( void ) const { return( 1 ); }

/*--------------------------------------------------------------------------*/

 /// returns the number of lines of the network
 /** Returns the number of lines of the network. This should just
  * be equivalent to get_NetworkData()->get_number_lines(), but the base
  * NetworkBlock class does not handle it, and therefore it assumes the network
  * is a bus, i.e., get_number_nodes() == 1, and returns 0 (no self-loops are
  * allowed, hence there is no line to be made with a single node). */

 virtual Index get_number_lines( void ) const { return( 0 ); }

/*--------------------------------------------------------------------------*/

 /// returns the number of intervals
 /** Returns the number of intervals spanned by this NetworkBlock. This should
  * just be equivalent to get_NetworkData()->get_number_intervals(), but the base
  * NetworkBlock class does not handle it, and therefore it assumes the network
  * handle just one time horizon and returns 1. */

 virtual Index get_number_intervals( void ) const { return( 1 ); }

/*--------------------------------------------------------------------------*/

 /// returns the NetworkData object
 /** The method of the base class always returns nullptr, because the base
  * class does not handle the NetworkData object. This is OK for derived
  * classes that only handle the "bus" case. */

 virtual NetworkData * get_NetworkData( void ) const {
  return nullptr;
 }

/*--------------------------------------------------------------------------*/

 /// returns the matrix of active demands
 /** Method for returning the active demand for the given interval, which is
  * assumed to have size get_number_intervals() per get_number_nodes().
  * There are two possible cases:
  *
  * - if the matrix only has one row (i.e., the first dimension has size 1),
  *   then the active demand for each user u is D[ 0 , u ] for all intervals t,
  *   which means that the second dimension has size get_number_nodes().
  *   This will be the default case;
  *
  * - otherwise, the matrix has size get_number_intervals() per
  *   get_number_nodes(), then the D[ i , u ] represents the active demand
  *   for the problem at time t for each user u, e.g., ECNetwork case;
  *
  * @param i The interval wrt the vector of demands for each user is
  *          returned. */

 virtual const double * get_active_demand( Index i = 0 ) const {
  return nullptr;
 }

/*--------------------------------------------------------------------------*/

 /// returns the constant term

 const double & get_const_term( void ) const {
  return f_const_term;
 }

/**@} ----------------------------------------------------------------------*/
/*----------- METHODS FOR READING THE Variable OF THE NetworkBlock ---------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the Variable of the NetworkBlock
 * @{ */

 /// returns the matrix of node injection variables
 /** Method for returning the node injection variables for the given interval,
  * which is assumed to have size get_number_intervals() per get_number_nodes().
  * There are two possible cases:
  *
  * - if the matrix only has one row (i.e., the first dimension has size 1),
  *   then the node injection for each user u is I[ 0 , u ] for all intervals t,
  *   which means that the second dimension has size get_number_nodes().
  *   This will be the default case;
  *
  * - otherwise, the matrix has size get_number_intervals() per
  *   get_number_nodes(), then the I[ i , u ] represents the node injection
  *   for the problem at time t for each user u, e.g., ECNetwork case;
  *
  * @param i The interval wrt the vector of node injections for each user is
  *          returned. */

 virtual ColVariable * get_node_injection( Index i = 0 ) {
  return nullptr;
 }

/**@} ----------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Solution
 * @{ */

 /// returns a Solution representing the current solution of this NetworkBlock
 /** This method must construct and return a (pointer to a) Solution object
  * representing the current "solution state" of this NetworkBlock. The base
  * NetworkBlock class defaults to ColVariableSolution, RowConstraintSolution,
  * and ColRowSolution, but :NetworkBlock may make different choices.
  *
  * The parameter for deciding which kind of Solution must be returned is a
  * single int value. If this value is
  *
  * - 1, then a RowConstraintSolution is returned;
  *
  * - 2, then a ColRowSolution is returned;
  *
  * - any other value, then a ColVariable Solution is returned.
  *
  * This value is to be found as:
  *
  * - if solc is not nullptr and it is a SimpleConfiguration< int >, then it
  *   is solc->f_value;
  *
  * - otherwise, if f_BlockConfig is not nullptr,
  *   f_BlockConfig->f_solution_Configuration is not nullptr and it is a
  *   SimpleConfiguration< int >, then it is
  *   f_BlockConfig->f_solution_Configuration->f_value;
  *
  * - otherwise, it is 0. */

 Solution * get_Solution( Configuration * solc = nullptr ,
                          bool emptys = true ) override;

/**@} ----------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

 virtual void set_active_demand( std::vector< double >::const_iterator values ,
                                 Subset && subset , bool ordered ,
                                 ModParam issuePMod , ModParam issueAMod ) = 0;

 virtual void set_active_demand( std::vector< double >::const_iterator values ,
                                 Range rng , ModParam issuePMod ,
                                 ModParam issueAMod ) = 0;

/*--------------------------------------------------------------------------*/
/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED METHODS OF THE CLASS ---------------------*/
/*--------------------------------------------------------------------------*/

 /// states that the Variable of the NetworkBlock have been generated
 void set_variables_generated( void ) { AR |= HasVar; }

 /// states that the Constraint of the NetworkBlock have been generated
 void set_constraints_generated( void ) { AR |= HasCst; }

 /// states that the Objective of the NetworkBlock has been generated
 void set_objective_generated( void ) { AR |= HasObj; }

 /// indicates whether the Variable of the NetworkBlock have been generated
 bool variables_generated( void ) const { return( AR & HasVar ); }

 /// indicates whether the Constraint of the NetworkBlock have been generated
 bool constraints_generated( void ) const { return( AR & HasCst ); }

 /// indicates whether the Objective of the NetworkBlock has been generated
 bool objective_generated( void ) const { return( AR & HasObj ); }

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

/*---------------------------------- data ----------------------------------*/

 /// the constant term
 double f_const_term;

/*-------------------------------- variables -------------------------------*/



/*------------------------------- constraints ------------------------------*/



/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE FIELDS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 ///< bit-wise coded: what abstract is there
 unsigned char AR{};

 static constexpr unsigned char HasVar = 1;
 ///< first bit of AR == 1 if the Variables have been constructed
 static constexpr unsigned char HasCst = 2;
 ///< second bit of AR == 1 if the Constraints have been constructed
 static constexpr unsigned char HasObj = 4;
 ///< third bit of AR == 1 if the Objective has been constructed

};  // end( class( NetworkBlock ) )

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS NetworkBlockMod --------------------------*/
/*--------------------------------------------------------------------------*/

/// Derived class from Modification for modifications to a NetworkBlock
class NetworkBlockMod : public Modification
{

 public:

 /// Public enum for the types of NetworkBlockMod
 enum NetB_mod_type
 {
  eSetActD = 0  ///< Set max power values
 };

 /// Constructor, takes the NetworkBlock and the type
 NetworkBlockMod( NetworkBlock * const fblock ,
                  const int type )
  : f_Block( fblock ) , f_type( type ) {}

 ///< Destructor, does nothing
 virtual ~NetworkBlockMod( void ) override = default;

 /// returns the Block to which the Modification refers
 Block * get_Block( void ) const override { return( f_Block ); }

 /// Accessor to the type of modification
 int type( void ) { return( f_type ); }

 protected:

 /// prints the NetworkBlockMod
 void print( std::ostream & output ) const override {
  output << "NetworkBlockMod[" << this << "]: ";
  switch( f_type ) {
   default:
    output << "Set active demand values ";
  }
 }

 NetworkBlock * f_Block{};
 ///< pointer to the Block to which the Modification refers

 int f_type; ///< type of modification

};  // end( class( NetworkBlockMod ) )

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS NetworkBlockRngdMod ------------------------*/
/*--------------------------------------------------------------------------*/
/// derived from NetworkBlockMod for "ranged" modifications
class NetworkBlockRngdMod : public NetworkBlockMod
{

 public:

 /// constructor: takes the NetworkBlock, the type, and the range
 NetworkBlockRngdMod( NetworkBlock * const fblock ,
                      const int type ,
                      Block::Range rng )
  : NetworkBlockMod( fblock , type ) , f_rng( rng ) {}

 /// destructor, does nothing
 virtual ~NetworkBlockRngdMod() override = default;

 /// accessor to the range
 Block::c_Range & rng( void ) { return( f_rng ); }

 protected:

 /// prints the NetworkBlockRngdMod
 void print( std::ostream & output ) const override {
  NetworkBlockMod::print( output );
  output << "[ " << f_rng.first << ", " << f_rng.second << " )" << std::endl;
 }

 Block::Range f_rng; ///< the range
};  // end( class( NetworkBlockRngdMod ) )

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS NetworkBlockSbstMod -----------------------*/
/*--------------------------------------------------------------------------*/

/// derived from NetworkBlockMod for "subset" modifications
class NetworkBlockSbstMod : public NetworkBlockMod
{

 public:

 /// constructor: takes the NetworkBlock, the type, and the subset
 NetworkBlockSbstMod( NetworkBlock * const fblock ,
                      const int type ,
                      Block::Subset && nms )
  : NetworkBlockMod( fblock , type ) , f_nms( std::move( nms ) ) {}

 /// destructor, does nothing
 virtual ~NetworkBlockSbstMod() override = default;

 /// accessor to the subset
 Block::c_Subset & nms( void ) { return( f_nms ); }

 protected:

 /// prints the NetworkBlockSbstMod
 void print( std::ostream & output ) const override {
  NetworkBlockMod::print( output );
  output << "(# " << f_nms.size() << ")" << std::endl;
 }

 Block::Subset f_nms;  ///< the subset

};  // end( class( NetworkBlockSbstMod ) )

/*--------------------------------------------------------------------------*/

}  /* namespace SMSpp_di_unipi_it */

/*--------------------------------------------------------------------------*/

#endif /* NetworkBlock.h included */

/*--------------------------------------------------------------------------*/
/*------------------------ End File NetworkBlock.h -------------------------*/
/*--------------------------------------------------------------------------*/
