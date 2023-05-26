/*--------------------------------------------------------------------------*/
/*--------------------------- File DCNetworkBlock.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 *
 * Header file for class DCNetworkBlock, which derives from NetworkBlock and
 * defines the standard linear constraints corresponding to the "DC model"
 * of the transmission network in the Unit Commitment problem.
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
 * \copyright &copy; by Antonio Frangioni, Ali Ghezelsoflu,
 *                      Rafael Durbano Lobato
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __DCNetworkBlock
 #define __DCNetworkBlock
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Block.h"

#include "FRowConstraint.h"

#include "OneVarConstraint.h"

#include "NetworkBlock.h"

#include "FRealObjective.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)

namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS DCNetworkBlock ---------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// a transmission NetworkBlock, i.e., a "DC" transmission network
/** The DCNetworkBlock class derives from NetworkBlock, and defines the
 * standard linear constraints corresponding to the "DC model" of the
 * transmission network in the Unit Commitment problem. Generally, there exist
 * three different kinds of DCNetworkBlock:
 *
 * - DCNetworkBlock with just HVDC lines; where the susceptance for all lines
 *   is equal to zero. It's also known as the Net Transfer Capacity (NTC)
 *   model.
 *
 * - DCNetworkBlock with just AC lines; where the susceptance for all lines
 *   is a non-zero value.
 *
 * - DCNetworkBlock of an hybrid AC/HVDC grid (both AC and HVDC lines). This
 *   is a combination of first and second cases, where for some lines (not all
 *   of them) may have zero susceptance. */

class DCNetworkBlock : public NetworkBlock
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
 * NetworkBlock defines two main public types:
 *
 * - line_type, an enum defining the types of lines present in the network.
 *
 * - DCNetworkData, a small auxiliary class to bunch the basic electrical data
 *   of the transmission network.
 * @{ */

 /// public enum for defining the types of lines of the network
 enum line_type
 {
  kNone = 0 ,  ///< no line
  kAC ,        ///< AC lines
  kHVDC ,      ///< HVDC lines
  kAC_HVDC     ///< AC and HVDC lines
 };

/*--------------------------------------------------------------------------*/
/*-------------------- CLASS NetworkBlock::NetworkData ---------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
 /// auxiliary class holding basic data about the transmission network
 /** The DCNetworkData class is a nested sub-class which only serves to have a
  * quick way to load all the basic data (topology and electrical
  * characteristics) that describe the transmission network. The rationale is
  * that while often the network does not change during the (short) time
  * horizon of UC, it makes sense to allow for this to happen. This means that
  * individual NetworkBlock objects may in principle have different
  * DCNetworkData, but most often they can share the same. By bunching all the
  * information together we make it easy for this sharing to happen. */

 class DCNetworkData : public NetworkBlock::NetworkData
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

  /// constructor of DCNetworkData, does nothing
  DCNetworkData( void ) : NetworkBlock::NetworkData() {}

  /// copy constructor of DCNetworkData, does nothing
  DCNetworkData( NetworkData * dc_network_data ) {}

  /// destructor of DCNetworkData: it is virtual, and empty
  virtual ~DCNetworkData() = default;

/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 * @{ */

  /// deserialize a DCNetworkData out of a netCDF::NcGroup
  /** Deserialize a DCNetworkData out of a netCDF::NcGroup, which should contain
   * the following:
   *
   * - The dimension "NumberNodes" containing the number of nodes in the
   *   problem; this dimension is optional, if it is not provided then it is
   *   taken to be equal to 1.
   *
   * If NumberNodes == 1 (equivalently, it is not provided), the network is a
   * "bus" formed of only one node, and therefore all the subsequent information
   * need not to be present since it is not loaded. If NumberNodes > 1, then all
   * the subsequent information is mandatory:
   *
   * - The dimension "NumberLines" containing the number of lines in the
   *   transmission network.
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
   *   the dimension "NumberLines". This is meant to represent the vector
   *   MxP[ l ] that, for each line l, contains the minimum power flow at
   *   line l (note that this is typically a negative number as lines are
   *   bi-directional, see above).
   *
   * - The variable "MaxPowerFlow", of type netCDF::NcDouble and indexed over
   *   the dimension "NumberLines". This is meant to represent the vector
   *   MxP[ l ] that, for each line l, contains the maximum power flow at
   *   line l (a non-negative number).
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
   *   then it's taken to be zero.
   *
   * - The variable "NodeName", of type netCDF::NcString() and indexed over the
   *   dimension "NumberNodes". Its i-th entry, namely NodeName[ i ], contains
   *   the name of the i-th node. This variable is optional.
   *
   * - The variable "LineName", of type netCDF::NcString() and indexed over the
   *   dimension "NumberLines". Its i-th entry, namely LineName[ i ], contains
   *   the name of the i-th transmission line. This variable is optional. */

  virtual void deserialize( const netCDF::NcGroup & group ) override;

/**@} ----------------------------------------------------------------------*/
/*------------ METHODS FOR READING THE DATA OF THE DCNetworkData -----------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the DCNetworkData
 * @{ */

  /// returns the number of lines of the network
  /** Method for returning the number of lines of the network. When
   * get_number_nodes() == 1 (the network is a bus), get_number_lines() == 0
   * (no self-loops are allowed, hence there is no line to be made with a single
   * node). */

  Index get_number_lines( void ) const { return( f_number_lines ); }

/*--------------------------------------------------------------------------*/

  /// returns the vector of start lines
  /** Method for returning the vector of starting point of each line. This
   * vector may have empty size (bus network) or the size of number of lines,
   * then there are two possible cases:
   *
   * - if f_number_nodes == 1, this vector has empty size which means there is
   *   no line at network (bus network), and this vector is not needed to be
   *   defined.
   *
   * - if f_number_nodes > 1, this vector have size of f_number_lines and each
   *   element of the vectors gives starting point of each line in the network.
   */

  const std::vector< Index > & get_start_line( void ) const {
   return( v_start_line );
  }

/*--------------------------------------------------------------------------*/

  /// returns vector of end lines
  /** Method for returning the vector of ending point of each line. This vector
   * may have empty size (bus network) or the size of number of lines, then
   * there are two possible cases:
   *
   * - if f_number_nodes == 1, this vector has empty size which means there is
   *   no line at network (bus network), and this vector is not needed to be
   *   defined.
   *
   * - if f_number_nodes > 1, this vector have size of f_number_lines and each
   *   element of the vectors gives ending point of each line in the network.
   */

  const std::vector< Index > & get_end_line( void ) const {
   return( v_end_line );
  }

/*--------------------------------------------------------------------------*/

  /// returns vector of the minimum power flow
  /** Method for returning the vector of minimum power flow of each line. This
   * vector may have empty size (bus network) or the size of number of nodes,
   * then there are two possible cases:
   *
   * - if f_number_lines == 0, this vector has empty size which means there is
   *   no line at network (bus network).
   *
   * - if f_number_lines >= 1, this vector have size of f_number_lines and each
   *   element of the vector gives minimum power flow of each line in the
   *   network. */

  const std::vector< double > & get_min_power_flow( void ) const {
   return( v_min_power_flow );
  }

/*--------------------------------------------------------------------------*/
  /// returns minimum power flow of the given \p line
  /** This method returns the minimum power flow of the given \p line.
   *
   * @return the minimum power flow of the given \p line. */

  double get_min_power_flow( Index line ) const {
   assert( line < get_number_lines() );
   if( v_min_power_flow.empty() )
    return( 0 );
   return( v_min_power_flow[ line ] );
  }

/*--------------------------------------------------------------------------*/
  /// returns vector of the maximum power flow
  /** Method for returning the vector of maximum power flow of each line. This
   * vector may have empty size (bus network) or the size of number of nodes,
   * then there are two possible cases:
   *
   * - if f_number_lines == 0, this vector has empty size which means there is
   *   no line at network (bus network).
   *
   * - if f_number_lines >= 1, this vector have size of f_number_lines and each
   *   element of the vector gives maximum power flow of each line in the
   *   network. */

  const std::vector< double > & get_max_power_flow( void ) const {
   return( v_max_power_flow );
  }

/*--------------------------------------------------------------------------*/
  /// returns maximum power flow of the given \p line
  /** This method returns the maximum power flow of the given \p line.
   *
   * @return the maximum power flow of the given \p line. */

  double get_max_power_flow( Index line ) const {
   assert( line < get_number_lines() );
   if( v_max_power_flow.empty() )
    return( 0 );
   return( v_max_power_flow[ line ] );
  }

/*--------------------------------------------------------------------------*/
  /// returns vector of the susceptances
  /** Method for returning the vector of susceptances for each line. This vector
   * may have empty size (bus network) or the size of number of nodes, then
   * there are two possible cases:
   *
   * - if f_number_lines == 0, this vector has empty size which means there is
   *   no line at network (bus network).
   *
   * - if f_number_lines >= 1, this vector has size of f_number_lines and each
   *   element of the vectors gives the Susceptance value for each line in the
   *   network. */

  const std::vector< double > & get_susceptance( void ) const {
   return( v_susceptance );
  }

/*--------------------------------------------------------------------------*/
  /// returns vector of the network cost
  /** Method for returning the vector of network cost for each line. This vector
   * may have empty size (bus network) or the size of number of lines, then
   * there are two possible cases:
   *
   * - if f_number_lines == 0, this vector has empty size which means there is
   *   no line at network (bus network).
   *
   * - if f_number_lines >= 1, this vector has size of f_number_lines and each
   *   element of the vectors gives the network cost value for each line in the
   *   network. */

  const std::vector< double > & get_network_cost( void ) const {
   return( v_network_cost );
  }

/*--------------------------------------------------------------------------*/
  /// returns the types of lines in the network
  /** This method returns the types of lines present in the network. */

  line_type get_lines_type( void ) const {
   if( get_number_lines() == 0 )
    return( kNone );
   if( std::all_of( v_susceptance.cbegin() , v_susceptance.cend() ,
                    []( double s ) { return( s == 0.0 ); } ) )
    return( kHVDC );
   if( std::all_of( v_susceptance.cbegin() , v_susceptance.cend() ,
                    []( double s ) { return( s != 0.0 ); } ) )
    return( kAC );
   return( kAC_HVDC );
  }

/*--------------------------------------------------------------------------*/

  /// returns the vector containing the name of the nodes
  const std::vector< std::string > & get_node_names( void ) const {
   return( v_node_names );
  }

/*--------------------------------------------------------------------------*/

  /// returns the vector containing the name of the lines
  const std::vector< std::string > & get_line_names( void ) const {
   return( v_line_names );
  }

/**@} ----------------------------------------------------------------------*/
/*-------------------- METHODS FOR SAVING THE DCNetworkData ----------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the DCNetworkData
 * @{ */

  /// serialize a DCNetworkData out of a netCDF::NcGroup
  /** Serialize a DCNetworkData out of a netCDF::NcGroup to the specific
   * format of a DCNetworkData. See NetworkBlock::deserialize( netCDF::NcGroup
   * ) for details of the format of the created netCDF group. */

  virtual void serialize( netCDF::NcGroup & group ) const override;

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

  protected:

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED METHODS OF THE CLASS ---------------------*/
/*--------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

  /// number of lines of the network
  Index f_number_lines{};

  /// vector of starting lines
  std::vector< Index > v_start_line;

  /// vector of ending lines
  std::vector< Index > v_end_line;

  /// vector to store the susceptance of each line of the network
  std::vector< double > v_susceptance;

  /// vector to store the minimum power flow at each line
  std::vector< double > v_min_power_flow;

  /// vector to store the maximum power flow at each line
  std::vector< double > v_max_power_flow;

  /// vector to store the network cost at each line
  std::vector< double > v_network_cost;

  std::vector< std::string > v_node_names;  ///< Node names

  std::vector< std::string > v_line_names;  ///< Line names

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

  private:

/*--------------------------------------------------------------------------*/
/*-------------------- PRIVATE FIELDS OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/



  SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/
/*---------------------- PRIVATE METHODS OF THE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/



 };  // end( class( DCNetworkData ) )

/**@} ----------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 * @{ */

 /// constructor of DCNetworkBlock
 /** Constructor of DCNetworkBlock, taking possibly a pointer of its
  * father Block. */

 explicit DCNetworkBlock( Block * f_block = nullptr )
  : NetworkBlock( f_block ) , f_NetworkData( nullptr ) {}

/*--------------------------------------------------------------------------*/
 /// destructor of DCNetworkBlock

 virtual ~DCNetworkBlock() override;

/*--------------------------------------------------------------------------*/
 /// generate the abstract variables of the DCNetworkBlock
 /** The size of node injection variable is the  number of intervals spanned
  * this DCNetworkBlock, i.e., 1, by the number of nodes, which can be read:
  *
  * - if NetworkData object is not provided (basically, "NumberNodes" is not
  *   provided or it is == 1) then the network is taken to have only one node
  *   (a bus) and there is only one node injection variable.
  *
  * - if the NetworkData object is present (either in the NcGroup or because it
  *   has been passed and NumberNodes > 1) then this variable has size
  *   "NumberNodes", which can be read via NetworkData::get_number_nodes().
  *
  * Depending on the susceptance for each line of the network, the
  * DCNetworkBlock class may have a power flow variable or not. In other word,
  * if the susceptance is equal to zero (or not defined), the corresponding line
  * is a HVDC line and it must have the power flow variable. It means, each
  * HVDC line correspond to a power flow variable, then for the Net Transfer
  * Capacity (NTC) model all lines must have a power flow variable. If the
  * susceptance value is a non-zero value, the corresponding line is called AC
  * and there is no needed to define the power flow variable for that line.
  * Therefore, in the case of pure AC line there is no needed to define power
  * flow variables. Consequently, for the mixed case AC-HVDC, the power flow
  * variable must define just for HVDC lines. Similarly, depending on the
  * NetworkCost for each line of the network, the DCNetworkBlock class may have
  * an auxiliary variable or not. In other word, if the NetworkCost is equal to
  * zero (or not defined), the auxiliary variable and corresponding constraints
  * will not be defined. */

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate abstract constraints of DCNetworkBlock
 /** Three different kinds of DCNetworkBlock constraints are defined as below.
  * The topology of the transmission network is defined by a set of nodes
  * \f$ N \f$ and a set of lines \f$ L \f$. Moreover, it's assumed that
  * \f$ P^{mn}_l \f$ and \f$ P^{mx}_l \f$ are minimum and maximum power flows
  * at each line \f$ l \in L \f$ and \f$ D^{ac}_n \f$ is active power demand
  * at node \f$ n \in N \f$ in the network respectively. The node injection
  * variable of each node \f$ n \in N \f$ and the power flows variable and an
  * auxiliary variable (which is not be defined if there is no network cost), of
  * each line \f$ l \in L \f$ are defined as \f$ S_n \f$, \f$ F_l \f$ and
  * \f$ V_l \f$ respectively.
  *
  * - DCNetworkBlock with just HVDC lines or the Net Transfer Capacity (NTC)
  *   model.
  *   In this special case the susceptance value for each line is equal to
  *   zero. In fact, this corresponds to a model with a single connected grid
  *   composed of HVDC lines only. In this case, the flow limit equations
  *   define as:
  *
  *   \f[
  *    \kappa_l P^{mn}_l \leq F_l  \leq \kappa_l P^{mx}_l
  *                                           \quad l \in \mathcal{L} \quad (1)
  *   \f]
  *
  *   where in each line \f$ l \f$, \f$ n \f$ and \f$ n' \f$ are supposed to
  *   be the start and the end point of that respectively. Besides, the
  *   following link between power flows and injected power at each node of the
  *   grid:
  *
  *   \f[
  *     \sum_{l=(n,n') } F_l - \sum_{l=(n',n)} F_l = S_n - D_n
  *                                           \quad n \in \mathcal{N} \quad (2)
  *   \f]
  *
  *   Moreover, when NetworkCost for each line is not equal to zero, DCNetwork
  *   will have an objective function which is equal to multiplying
  *   NetworkCost by the absolute value of the power flows variable. To relax
  *   the absolute value, an auxiliary variable and constraints as below are
  *   needed:
  *
  *   \f[
  *    F_l \leq V_l                           \quad l \in \mathcal{L} \quad (3)
  *   \f]
  *
  *   \f[
  *    -V_l \leq F_l                          \quad l \in \mathcal{L} \quad (4)
  *   \f]
  *
  * - DCNetworkBlock with just AC lines model.
  *   By considering a \f$ |L| \times |N| \f$ matrix \f$ B \f$ which
  *   constitutes the so-called Power Transfer Distribution Factor matrix
  *   (PTDF-matrix) which represents the linear relationship between power
  *   injections at each node of the grid and active power flows through the
  *   transmission lines.
  *
  *   The flow limit equations can be written as follow:
  *
  *   \f[
  *    P^{mn}_l\leq \sum_{ n \in N} B_{(l , n)}
  *    (S_n - D^{ac}_n) \leq  P^{mx}_l
  *                                           \quad l \in \mathcal{L} \quad (5)
  *   \f]
  *
  * - DCNetworkBlock of an hybrid AC/HVDC grid (both AC and HVDC lines).
  *   This is the case of an hybrid grid constituted of both AC and HVDC (High
  *   Voltage Direct Current) lines. The DC lines are characterized by the
  *   fact that the flow passing through those lines is fully controllable.
  *   However, this flow still has an impact on the flows passing through
  *   connected AC lines. To use the matrix formalism, we first introduce some
  *   additional notations:
  *
  *   - Lines of the grid are indexed by \f$ l = 1,··· ,|L^{ac}| \f$ for AC
  *     lines, while indexes \f$ l = |L^{ac}| + 1,··· ,|L^{ac}|+|L^{dc}| \f$
  *     refer to DC lines.
  *
  *   - For any \f$ l \in \{1,··· ,|L^{ac}| \}\f$ and
  *     \f$ k \in \{1,··· ,|L^{dc}| \} \f$, and put \f$ \ell(l) \f$ the pair of
  *     nodes related by the AC line indexed by \f$ l \f$ and
  *     \f$ \ell(k+|L^{ac}|) \f$ denotes the pair of nodes related by the DC
  *     line indexed by \f$ k+|L^{dc}| \f$.
  *
  *   - \f$ A^{dc} \f$ denotes the \f$ |L^{dc}| \times |N| \f$ incidence
  *     matrix induced by DC lines of the grid and the
  *     \f$ |L| \times (|L^{dc}| + |N|) \f$ matrix of A, where obtained by
  *     concatenation of bloc matrices as follows:
  *
  *     \f[
  *      A = \left[
  *      \begin{array}{cc}
  *      B & -B(A^{dc})^T \\
  *      0_{|L^{dc}| \times |N|} & I_{|L^{dc}| \times |L^{dc}|}
  *      \end{array}\right]                                           \quad (6)
  *     \f]
  *
  *     where \f$ 0_{|L^{dc}| \times |N|} \f$  denotes the
  *     \f$ |L^{dc}| \times |N| \f$ zero matrix and
  *     \f$ I_{|L^{dc}| \times |L^{dc}|}\f$ the \f$|L^{dc}| \times |L^{dc}|\f$
  *     identity matrix. Therefore, the flow limit equations can be transformed
  *     into:
  *
  *     \f[
  *      P^{mn} \leq A \left[
  *      \begin{array}{c}
  *      a \\
  *      b
  *      \end{array}\right]
  *      \leq  P^{mx}                                                 \quad (7)
  *     \f]
  *
  *     where the vector \f$ a = (a_n)_{n = 1, ... , |N| }\f$ and
  *     \f$ b = (b_m)_{m = 1, ... , |L^{dc}| }\f$ are such that for any
  *     \f$ n \in \{ 1, ... , |N|\}\f$ and \f$ m \in \{ 1, ... , |L^{dc}|\}\f$
  *     which \f$ a_n = \sum_{ i \in I_n} p^{ac}_i - D^{ac}_n \f$ and
  *     \f$ b_m = p_{m + |L^{ac}|} = p^{dc}_{\ell(m + |L^{ac}|)}\f$. */

 void generate_abstract_constraints( Configuration * stcc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the objective of the DCNetworkBlock
 /** Method that generates the objective of the DCNetworkBlock.
  *
  * - Objective function: the objective function of the DCNetworkBlock
  *   is given as below:
  *
  *   \f[
  *    \min ( \sum_{ l \in \mathcal{L} } ( NC_l V_l ) )
  *   \f]
  *
  *   where \f$ NC_l \f$, is a network cost and \f$ V_l \f$ is the auxiliary
  *   variable. */

 void generate_objective( Configuration * objc = nullptr ) override;

/**@} ----------------------------------------------------------------------*/
/*---------------- Methods for checking the DCNetworkBlock -----------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for checking solution information in the DCNetworkBlock
 * @{ */

 /// returns true if the current solution is (approximately) feasible
 /** This function returns true if and only if the solution encoded in the
  * current value of the Variable of this DCNetworkBlock is approximately
  * feasible within the given tolerance. That is, a solution is considered
  * feasible if and only if
  *
  * -# each ColVariable is feasible; and
  *
  * -# the violation of each Constraint of this DCNetworkBlock is not
  *    greater than the tolerance.
  *
  * Every Constraint of this DCNetworkBlock is a RowConstraint and its
  * violation is given by either the relative (see RowConstraint::rel_viol())
  * or the absolute violation (see RowConstraint::abs_viol()), depending on
  * the Configuration that is provided.
  *
  * The tolerance and the type of violation can be provided by either \p fsbc
  * or #f_BlockConfig->f_is_feasible_Configuration and they are determined as
  * follows:
  *
  * - If \p fsbc is not a nullptr and it is a pointer to a
  *   SimpleConfiguration< double >, then the tolerance is the value present
  *   in that SimpleConfiguration and the relative violation is considered.
  *
  * - If \p fsbc is not nullptr and it is a
  *   SimpleConfiguration< std::pair< double , int > >, then the tolerance is
  *   fsbc->f_value.first and the type of violation is determined by
  *   fsbc->f_value.second (any nonzero number for relative violation and
  *   zero for absolute violation);
  *
  * - Otherwise, if both #f_BlockConfig and
  *   f_BlockConfig->f_is_feasible_Configuration are not nullptr and the
  *   latter is a pointer to either a SimpleConfiguration< double > or to a
  *   SimpleConfiguration< std::pair< double , int > >, then the values of the
  *   parameters are obtained analogously as above;
  *
  * - Otherwise, by default, the tolerance is 0 and the relative violation
  *   is considered.
  *
  * This function currently considers only the abstract representation to
  * determine if the solution is feasible. So, the parameter \p useabstract is
  * currently ignored. If no abstract Variable has been generated, then this
  * function returns true. Moreover, if no abstract Constraint has been
  * generated, the solution is considered to be feasible with respect to the
  * set of Variable only. Notice also that, before checking if the solution
  * satisfies a Constraint, the Constraint is computed (Constraint::compute()).
  *
  * @param useabstract This parameter is currently ignored.
  *
  * @param fsbc The pointer to a Configuration that specifies the tolerance
  *             and the type of violation that must be considered. */

 bool is_feasible( bool useabstract = false ,
                   Configuration * fsbc = nullptr ) override;

/**@} ----------------------------------------------------------------------*/
/*---------- METHODS FOR READING THE DATA OF THE DCNetworkBlock ------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the NetworkBlock
 * @{ */

 /// returns the number of nodes
 /** Returns the number of nodes in the transmission network. If
  * get_NetworkData() returns nullptr, this is equivalent to
  * get_NetworkData()->get_number_nodes(). Otherwise, it assumes the network
  * is a bus and returns 1.
  *
  * @return the number of nodes in the network. */

 Index get_number_nodes( void ) const override {
  if( ! f_NetworkData )
   return( 1 );
  return( f_NetworkData->get_number_nodes() );
 }

/*--------------------------------------------------------------------------*/
 /// returns the number of lines of the network
 /** This function returns the number of lines in the transmission network. If
  * get_NetworkData() returns nullptr, this is equivalent to
  * get_NetworkData()->get_number_lines(). Otherwise, it returns zero.
  *
  * @return the number of lines in the network. */

 Index get_number_lines( void ) const {
  if( ! f_NetworkData )
   return( 0 );
  return( f_NetworkData->get_number_lines() );
 }

/*--------------------------------------------------------------------------*/
 /// returns the kappa constant associated with the given \p line
 /** This function returns the kappa constant associated with the given \p
  * line. This is the constant that multiplies the minimum and maximum flow in
  * the flow limit constraint associated with the given \p line.
  *
  * @param line The index of a line (between 0 and get_number_lines() - 1).
  *
  * @return The kappa constant associated with the given \p line. */

 double get_kappa( Index line ) const {
  if( v_kappa.empty() )
   return( 1 );
  assert( line < v_kappa.size() );
  return( v_kappa[ line ] );
 }

/*--------------------------------------------------------------------------*/
 /// returns the minimum power flow on the given \p line
 /** This function returns the minimum power flow on the given \p line. If
  * this DCNetworkBlock has no NetworkData, thus function returns
  * 0. Otherwise, it returns the minimum power flow specified by the
  * NetworkData object.
  *
  * @param line The index of a line.
  *
  * @return The minimum power flow on the given \p line. */

 double get_min_power_flow( Index line ) const {
  if( ! f_NetworkData )
   return( 0 );
  return( f_NetworkData->get_min_power_flow( line ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the maximum power flow on the given \p line
 /** This function returns the maximum power flow on the given \p line. If
  * this DCNetworkBlock has no NetworkData, thus function returns
  * 0. Otherwise, it returns the maximum power flow specified by the
  * NetworkData object.
  *
  * @param line The index of a line.
  *
  * @return The maximum power flow on the given \p line. */

 double get_max_power_flow( Index line ) const {
  if( ! f_NetworkData )
   return( 0 );
  return( f_NetworkData->get_max_power_flow( line ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns a pointer to the DCNetworkData
 /** Return a pointer to the DCNetworkData. */

 NetworkData * get_NetworkData( void ) const override {
  return( f_NetworkData );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of active demands
 /** Returns the active demand for the given interval, which is assumed to
  * have size by get_number_nodes().
  *
  * @param interval The interval wrt the vector of demands for each user is
  *                 returned. */

 const double * get_active_demand( Index interval = 0 ) const override {
  if( v_ActiveDemand.empty() )
   return( nullptr );
  return( &( v_ActiveDemand.front() ) );
 }

/**@} ----------------------------------------------------------------------*/
/*---------- METHODS FOR READING THE Variable OF THE DCNetworkBlock --------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the Variable of the DCNetworkBlock
 * @{ */

 /// returns the vector of power flow variables
 /** The returned std::vector< ColVariable >, say F, contains the power flow
  * variables and is indexed over the dimension "NumberLines". There are two
  * possible cases:
  *
  * - if F is empty(), then this variable is not defined;
  *
  * - otherwise, F must have f_number_lines rows and F[ l ] is the power flow
  *   variable for line l. */

 const std::vector< ColVariable > & get_power_flow( void ) const {
  return( v_power_flow );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of auxiliary variables
 /** The returned std::vector< ColVariable >, say V, contains the auxiliary
  * variables and is indexed over the dimension "NumberLines". There are two
  * possible cases:
  *
  * - if V is empty(), then this variable is not defined;
  *
  * - otherwise, V must have f_number_lines rows and V[ l ] is the auxiliary
  *   variable for line l. */

 const std::vector< ColVariable > & get_auxiliary_variable( void ) const {
  return( v_auxiliary_variable );
 }

/**@} ----------------------------------------------------------------------*/
/*--------- METHODS FOR READING THE Constraint OF THE DCNetworkBlock -------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the Constraint of the DCNetworkBlock
 * @{ */

 /// returns the vector of power flow limit constraints
 /** This function returns a const reference to the vector of power flow limit
  * constraints. The i-th element of this vector is a FRowConstraint for the
  * i-th line of the network. */

 const std::vector< FRowConstraint > &
 get_power_flow_limit_constraints( void ) const {
  if( ! f_NetworkData )
   throw( std::logic_error( "DCNetworkBlock::get_power_flow_limit_constraints:"
                            " DCNetworkData has not been set." ) );

  switch( f_NetworkData->get_lines_type() ) {
   case( kAC ):
    return( v_AC_power_flow_limit_const );
   case( kAC_HVDC ):
   default:
    return( v_AC_HVDC_power_flow_limit_const );
  }
 }

/*--------------------------------------------------------------------------*/

 /// returns the vector of power flow limit HVDC bounds

 const std::vector< BoxConstraint > &
 get_power_flow_limit_HVDC_bounds( void ) const {
  if( ! f_NetworkData )
   throw( std::logic_error( "DCNetworkBlock::get_power_flow_limit_HVDC_bounds:"
                            " DCNetworkData has not been set." ) );
  return( v_HVDC_power_flow_limit_const );
 }

/**@} ----------------------------------------------------------------------*/
/*--------------- METHODS FOR MODIFYING THE DCNetworkBlock -----------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the DCNetworkBlock
 * @{ */

 void set_NetworkData( NetworkBlock::NetworkData * nd = nullptr ) override {
  // if there was a previous DCNetworkData, and it was local, delete it
  if( f_NetworkData && f_local_NetworkData )
   delete( f_NetworkData );

  f_NetworkData = static_cast< DCNetworkData * >( nd );
  f_local_NetworkData = false;
 }

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
  * it is written in v_ActiveDemand (which therefore is no longer empty),
  * otherwise it is left empty so that it can be set by this method. */

 void set_ActiveDemand(
  const std::vector< std::vector< double > > & v ) override {
  if( v_ActiveDemand.empty() )
   v_ActiveDemand = v[ 0 ];
 }

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 * @{ */

 /// deserialize a DCNetworkBlock out of a netCDF::NcGroup
 /** Deserialize a DCNetworkBlock out of a netCDF::NcGroup, which should
  * contain all the data necessary to describe a NetworkBlock (see
  * NetworkBlock::deserialize()) and possibly the following variable:
  *
  * - The variable "ActiveDemand", of type netCDF::NcDouble and indexed over
  *   the dimension "NumberNodes".
  *   If the NetworkData object description is present in the NcGroup this is
  *   the dimension "NumberNodes", but the NetworkData object is optional and it
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
  *   UCBlock::deserialize()];
  *
  * - The variable "Kappa", of type netCDF::NcDouble and either being a
  *   scalar or indexed over the number of lines. If this variable is a
  *   scalar, let say k, then it is assumed that Kappa[ l ] = k for each line
  *   l in {0, ..., get_number_lines() - 1}. For each line l in {0, ...,
  *   get_number_lines() - 1}, Kappa[ l ] is the constant that multiplies the
  *   minimum and maximum flow in the flow limit constraints. This variable is
  *   optional. If it is not provided, it is assumed that Kappa[ l ] == 1 for
  *   each line l in {0, ..., get_number_lines() - 1};
  *
  * - The variable "ConstantTerm", of type netCDF::NcDouble and containing the
  *   constant term. */

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
 /// loads the DCNetworkBlock instance from memory
 /** Like load( std::istream & ), if there is any Solver attached to this
  * DCNetworkBlock then a NBModification (the "nuclear option") is issued. */

 void load( std::istream & input , char frmt = 0 ) override {
  throw( std::logic_error( "DCNetworkBlock::load() not implemented yet" ) );
 }

/**@} ----------------------------------------------------------------------*/
/*-------------------- METHODS FOR SAVING THE DCNetworkBlock ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the DCNetworkBlock
 * @{ */

 /// extends Block::serialize( netCDF::NcGroup )
 /** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
  * NetworkBlock. See NetworkBlock::deserialize( netCDF::NcGroup ) for
  * details of the format of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/** @} ---------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

 /// set the kappa constants for the lines specified by \p subset
 /** This function sets the kappa constant of each line in the given \p
  * subset. The kappa constant of each line whose index is specified by the
  * i-th element in \p subset is given by the i-th element of the vector
  * pointed by \p values, i.e., it is given by the value pointed by (values +
  * i). The parameter \p ordered indicates whether the \p subset is ordered.
  *
  * @param values An iterator to a vector containing the kappa constants.
  *
  * @param subset The indices of the lines whose kappa constants are being
  *        modified.
  *
  * @param ordered It indicates whether \p subset is ordered.
  *
  * @param issuePMod It controls how physical Modification are issued.
  *
  * @param issueAMod It controls how abstract Modification are issued. */

 void set_kappa( MF_dbl_it values ,
                 Subset && subset ,
                 const bool ordered = false ,
                 c_ModParam issuePMod = eNoBlck ,
                 c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// set the kappa constants for the lines specified by \p rng
 /** This function sets the kappa constants for the lines in the given Range
  * \p rng. For each i in the given Range (up to the number of lines minus 1),
  * the kappa constant for line i is given by the element of the vector
  * pointed by \p values whose index is (i - rng.first), i.e., it is given by
  * the value pointed by (values + i - rng.first).
  *
  * @param values An iterator to a vector containing the kappa constants.
  *
  * @param rng A Range containing the indices of the lines whose kappa
  *        constants are being modified.
  *
  * @param issuePMod It controls how physical Modification are issued.
  *
  * @param issueAMod It controls how abstract Modification are issued. */

 void set_kappa( MF_dbl_it values ,
                 Range rng = Range( 0 , Inf< Index >() ) ,
                 c_ModParam issuePMod = eNoBlck ,
                 c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// set the active demand at the nodes specified by \p subset
 /** This function sets the active demand at each node in the given \p
  * subset. The active demand at the node whose index is specified by the i-th
  * element in \p subset is given by the i-th element of the vector pointed by
  * \p values, i.e., it is given by the value pointed by (values + i). The
  * parameter \p ordered indicates whether the \p subset is ordered.
  *
  * @param values An iterator to a vector containing the active demand.
  *
  * @param subset The indices of the nodes at which the active demand is being
  *        modified.
  *
  * @param ordered It indicates whether \p subset is ordered.
  *
  * @param issuePMod It controls how physical Modification are issued.
  *
  * @param issueAMod It controls how abstract Modification are issued. */

 void set_active_demand( MF_dbl_it values ,
                         Subset && subset ,
                         const bool ordered = false ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck ) override final;

/*--------------------------------------------------------------------------*/
 /// set the active demand at the nodes specified by \p rng
 /** This function sets the active demand at each node in the given Range \p
  * rng. For each i in the given Range (up to the number of nodes minus 1),
  * the active demand at node i is given by the element of the vector pointed
  * by \p values whose index is (i - rng.first), i.e., it is given by the
  * value pointed by (values + i - rng.first).
  *
  * @param values An iterator to a vector containing the active demand.
  *
  * @param rng A Range containing the indices of the nodes at which the active
  *        demand is being modified.
  *
  * @param issuePMod It controls how physical Modification are issued.
  *
  * @param issueAMod It controls how abstract Modification are issued. */

 void set_active_demand( MF_dbl_it values ,
                         Range rng = Range( 0 , Inf< Index >() ) ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck ) override final;

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED METHODS OF THE CLASS ---------------------*/
/*--------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

/*---------------------------------- data ----------------------------------*/

 /// the DCNetworkData object
 DCNetworkData * f_NetworkData;

 /// vector to store the demand of each node of the network
 std::vector< double > v_ActiveDemand;

 /// the kappa constant for each line
 std::vector< double > v_kappa;

/*-------------------------------- variables -------------------------------*/

 /// the power flow variables
 std::vector< ColVariable > v_power_flow;

 /// the auxiliary network cost variable
 std::vector< ColVariable > v_auxiliary_variable;

/*------------------------------- constraints ------------------------------*/

 /// AC power flow limit constraints
 std::vector< FRowConstraint > v_AC_power_flow_limit_const;

 /// AC/HVDC power flow limit constraints
 std::vector< FRowConstraint > v_AC_HVDC_power_flow_limit_const;

 /// HVDC power flow and node injection constraints
 std::vector< FRowConstraint > v_power_flow_injection_const;

 /// HVDC power flow auxiliary variable constraints
 boost::multi_array< FRowConstraint , 2 > v_power_flow_relax_abs;


 /// HVDC power flow limit constraints
 std::vector< BoxConstraint > v_HVDC_power_flow_limit_const;

 /// the node injection bound constraints
 std::vector< BoxConstraint > node_injection_bounds_const;


 /// the objective function
 FRealObjective objective;

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------- PRIVATE FIELDS OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/



 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/
/*---------------------- PRIVATE METHODS OF THE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

 static void static_initialization( void ) {

  /* Warning: Not all C++ compilers enjoy the template wizardry behind the
   * three-args version of register_method<> with the compact MS_*_*::args(),
   *
   * register_method< DCNetworkBlock >( "DCNetworkBlock::set_active_demand",
   *                                    &DCNetworkBlock::set_active_demand,
   *                                    MS_dbl_sbst::args() );
   *
   * so we just use the slightly less compact one with the explicit argument
   * and be done with it. */

  register_method< DCNetworkBlock , MF_dbl_it , Subset && , bool >(
   "DCNetworkBlock::set_active_demand" , &DCNetworkBlock::set_active_demand );

  register_method< DCNetworkBlock , MF_dbl_it , Range >(
   "DCNetworkBlock::set_active_demand" , &DCNetworkBlock::set_active_demand );
 }

};  // end( class( DCNetworkBlock ) )

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS DCNetworkBlockMod -------------------------*/
/*--------------------------------------------------------------------------*/

/// derived class from NetworkBlockMod for modifications to a DCNetworkBlock
class DCNetworkBlockMod : public NetworkBlockMod
{

 public:

 /// public enum for the types of DCNetworkBlockMod
 enum DCNetB_mod_type
 {
  eSetKappa = eNetBModLastParam ,  ///< set the kappa constants
  eDCNetBModLastParam  ///< first allowed parameter value for derived classes
  /**< Convenience value to easily allow derived classes to extend the set of
   * types of DCNetworkBlockMod. */
 };

 /// constructor, takes the DCNetworkBlock and the type
 DCNetworkBlockMod( DCNetworkBlock * const fblock , const int type )
  : NetworkBlockMod( fblock , type ) {}

 /// destructor, does nothing
 virtual ~DCNetworkBlockMod() override = default;

 /// returns the Block to which the Modification refers
 Block * get_Block( void ) const override { return( f_Block ); }

 protected:

 /// prints the DCNetworkBlockMod
 void print( std::ostream & output ) const override {
  output << "DCNetworkBlockMod[" << this << "]: ";
  switch( f_type ) {
   default:
    output << "Set active demand values ";
  }
 }

 DCNetworkBlock * f_Block{};
 ///< pointer to the Block to which the Modification refers

};  // end( class( DCNetworkBlockMod ) )

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS DCNetworkBlockRngdMod -----------------------*/
/*--------------------------------------------------------------------------*/

/// derived from DCNetworkBlockMod for "ranged" modifications
class DCNetworkBlockRngdMod : public DCNetworkBlockMod
{

 public:

 /// constructor: takes the DCNetworkBlock, the type, and the range
 DCNetworkBlockRngdMod( DCNetworkBlock * const fblock , const int type ,
                        Block::Range rng )
  : DCNetworkBlockMod( fblock , type ) , f_rng( rng ) {}

 /// destructor, does nothing
 virtual ~DCNetworkBlockRngdMod() override = default;

 /// accessor to the range
 Block::c_Range & rng( void ) { return( f_rng ); }

 protected:

 /// prints the DCNetworkBlockRngdMod
 void print( std::ostream & output ) const override {
  DCNetworkBlockMod::print( output );
  output << "[ " << f_rng.first << ", " << f_rng.second << " )" << std::endl;
 }

 Block::Range f_rng;  ///< the range

};  // end( class( DCNetworkBlockRngdMod ) )

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS DCNetworkBlockSbstMod ----------------------*/
/*--------------------------------------------------------------------------*/

/// derived from DCNetworkBlockMod for "subset" modifications
class DCNetworkBlockSbstMod : public DCNetworkBlockMod
{

 public:

 /// constructor: takes the DCNetworkBlock, the type, and the subset
 DCNetworkBlockSbstMod( DCNetworkBlock * const fblock , const int type ,
                        Block::Subset && nms )
  : DCNetworkBlockMod( fblock , type ) , f_nms( std::move( nms ) ) {}

 /// destructor, does nothing
 virtual ~DCNetworkBlockSbstMod() override = default;

 /// accessor to the subset
 Block::c_Subset & nms( void ) { return( f_nms ); }

 protected:

 /// prints the DCNetworkBlockSbstMod
 void print( std::ostream & output ) const override {
  DCNetworkBlockMod::print( output );
  output << "(# " << f_nms.size() << ")" << std::endl;
 }

 Block::Subset f_nms;  ///< the subset

};  // end( class( DCNetworkBlockSbstMod ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* DCNetworkBlock.h included */

/*--------------------------------------------------------------------------*/
/*-------------------- End File DCNetworkBlock.h ---------------------------*/
/*--------------------------------------------------------------------------*/
