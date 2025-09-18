/*--------------------------------------------------------------------------*/
/*--------------------------- File DCNetworkBlock.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 *
 * Header file for class DCNetworkBlock, which derives from NetworkBlock and
 * defines the standard linear constraints corresponding to the "DC model"
 * of the transmission network in the Unit Commitment problem.
 *
 * \author Wim van Ackooij \n
 *         EDF R&D OSIRIS \n
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Ali Ghezelsoflu \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Quentin Jacquet \n
 *         EDF R&D OSIRIS \n
 *
 * \author Rafael Durbano Lobato \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Rafael Durbano Lobato
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

#include "LinearFunction.h"

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
 *   of them) may have zero susceptance.
 */

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
  * information together we make it easy for this sharing to happen.
  */

class DCNetworkData : public NetworkData
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
 DCNetworkData( void ) : f_number_lines( 0 ) , f_reference_node( 0 ) ,
  f_lines_type( -1 ) , f_number_branches( 0 ) {}

 /// copy constructor of DCNetworkData, does nothing
 explicit DCNetworkData( const NetworkData * ) : f_number_lines( 0 ) ,
  f_reference_node( 0 ) , f_lines_type( -1 ) , f_number_branches( 0 ) {}

 /// destructor of DCNetworkData: it is virtual, and empty
 virtual ~DCNetworkData() override = default;

/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 * @{ */

 /// deserialize a DCNetworkData out of a netCDF::NcGroup
 /** Deserialize a DCNetworkData out of a netCDF::NcGroup, which should
  * contain the following:
  *
  * - The dimension "NumberNodes" containing the number of nodes in the
  *   problem; this dimension is optional, if it is not provided then it is
  *   taken to be equal to 1.
  *
  * If NumberNodes == 1 (equivalently, it is not provided), the network is
  * a "bus" formed of only one node, and therefore all the subsequent
  * information need not to be present since it is not loaded. If
  * NumberNodes > 1, then all the subsequent information is considered:
  *
  * - The dimension "NumberLines" containing the number of lines in the
  *   transmission network. Each line can be either a "regular" line / link
  *   / arc (one head node / bus, one tail node / bus) or a hyperarc (still
  *   one head node / bus, but possible multiple tail nodes / buses); see
  *   the (optional) dimension NumberBranches right next. The dimension is
  *   mandatory.
  *
  * - The dimension "NumberBranches" that is used to describe hyperarcs.
  *   The dimension is optional, if it is not defined then it is assumed that
  *   "NumberBranches" == "NumberLines", i.e., all lines are "regular".
  *   Otherwise, "NumberBranches" >= "NumberLines" (in fact, >) must hold
  *   since one single hyperarc is described by its multiple "branches", as
  *   detailed in "HyperArcID".
  *
  * - The variable "StartLine", of type netCDF::NcUint and indexed over the
  *   dimension "NumberBranches" (if it is defined, otherwise "NumberLines");
  *   the l-th entry of the variable is the starting point of the line (a
  *   number in 0, ..., NumberNodes - 1). Note that lines are not oriented,
  *   but the flow of energy is; that is, a positive flow along line l means
  *   that energy is being taken away from StartLine[ l ] and delivered to
  *   EndLine[ l ] (see next), a negative flow means vice versa. The variable
  *   is mandatory.
  *
  * - The variable "EndLine", of type netCDF::NcUint and indexed over the
  *   dimension "NumberBranches" (if it is defined, otherwise "NumberLines");
  *   the l-th entry of the variable is the ending point of the line (a number
  *   in 0, ..., NumberNodes - 1; lines are not oriented, but see above).
  *   StartLine[ l ] == EndLine[ l ] (a self-loop) is not allowed, but multiple
  *   lines between the same pair of nodes are. The variable is mandatory.
  *
  * - The variable "HyperArcID", of type netCDF::NcUint and indexed over the
  *   dimension "NumberBranches". The variable is mandatory if "NumberBranches"
  *   exists, and therefore "NumberBranches" > "NumberLines", and ignored
  *   otherwise. The variable is used to specify which entries of "StartLine"
  *   and "EndLine" are different "branches" that correspond to the same
  *   hyperarc. The entries of the variable are a number in 0, ..., 
  *   NumberLines - 1: HyperArcID[ i ] == l means that StartLine[ i ] and
  *   EndLine[ i ] describe one of the "branches" of the (hyper)line(arc) l.
  *   If a (hyper)line(arc) l has only one branch, i.e., HyperArcID[ i ] == l
  *   happens precisely for one index i, then l is a "regular" line. Note
  *   that, FOR EACH l = 0, ..., NumberLines - 1, THERE MUST BE AT LEAST ONE
  *   INDEX i such that HyperArcID[ i ] == l. If HyperArcID[ i ] == l happens
  *   for more than one index i, then l is a hyperarc (line). It is required
  *   that StartLine[ i ] == StartLine[ j ] and EndLine[ i ] == EndLine[ j ]
  *   for all pairs ( i , j ) such that HyperArcID[ i ] == HyperArcID[ j ],
  *   i.e., ALL "branches" MUST HAVE THE SAME "tail" and different heads.
  *
  * - The variable "MaxPowerFlow", of type netCDF::NcDouble and indexed over
  *   the dimension "NumberLines". This is meant to represent the vector
  *   MxP[ l ] that, for each line l, contains the maximum power flow at
  *   line l (a non-negative number). Note that if line l is a hyperarc (see
  *   "HyperArcID") the capacity is still one number representing the
  *   maximum amount of flow leaving the tail bus, although then some flow
  *   (not necessarily the same amount, see "Efficiency") can reach more than
  *   one head bus. The variable is optional, if not provided it is assumed
  *   that MxP[ l ] == 0 for all line l.
  *
  * - The variable "MinPowerFlow", of type netCDF::NcDouble and indexed over
  *   the dimension "NumberLines". This is meant to represent the vector
  *   MnP[ l ] that, for each line l, contains the minimum power flow at
  *   line l (note that this is typically a negative number as lines are
  *   bi-directional, see above). Note that if line l is a hyperarc (see
  *   "HyperArcID") the capacity is still one number representing the
  *   minimum amount of flow leaving the tail bus, although then some flow
  *   (not necessarily the same amount, see "Efficiency") can reach more than
  *   one head bus. The variable is optional, if not provided it is assumed
  *   that MnP[ l ] == 0 for all line l.
  *
  * - The variable "LineSusceptance", of type netCDF::NcDouble and indexed
  *   over the dimension "NumberLines". This is meant to represent the
  *   vector S[ l ] that, for each line l contains the susceptance of the
  *   network for the corresponding line l. Note that this variable is
  *   optional, for each line l if it is provided then it is assumed that
  *   S[ l ] != 0, otherwise it is assumed that S[ l ] == 0. In fact, when
  *   S[ l ] != 0 this corresponds to a model with AC lines, and when for
  *   each line l, it's not defined or S[ l ] == 0, then it corresponds to
  *   a single connected grid composed of HVDC lines only which is also
  *   known as the Net Transfer Capacity (NTC) model. Also, note that
  *   ALL HYPERARCS MUST HAVE 0 SUSCEPTANCE.
  *
  * - The dimension "ReferenceNode", that specifies which bus gets 0
  *   potential in Kirchhoff's equations. This changes the form of the PTDF
  *   matrix computed for the lines that have a nonzero LineSusceptance;
  *   although the problem should be mathematically equivalent whatever this
  *   choice is, numerically it may make a difference. The choice is
  *   obviously irrelevant for a pure HVDC network (when all LineSusceptance
  *   are 0), and in fact the dimension is optional: if not specified, the
  *   reference bus (if at all significant) is chosen as 0.
  *
  * - The variable "NetworkCost", of type netCDF::NcDouble and indexed over
  *   the dimension "NumberLines". This is meant to represent the vector
  *   NC[ l ] that, for each line l, contains the monetary cost to send one
  *   one unit of flow from StartLine[ l ] to EndLine[ l ]. Note that, if l
  *   is a hyperarc (see "HyperArcID"), the cost is still one number
  *   representing the unitary cost of one unit of flow leaving the tail bus,
  *   although then some flow (not necessarily the same amount, see
  *   "Efficiency") can reach more than one head bus.
  *
  * - The variable "Efficiency", of type netCDF::NcDouble indexed over the
  *   dimension "NumberBranches" (if it is defined, otherwise "NumberLines");
  *   Efficiency[ l ] represents the efficiency of branch l. This means that
  *   if X is the amount of flow leaving StartLine[ l ], then
  *   X * Efficiency[ l ] is the amount of flow reaching EndLine[ l ]. Note
  *   that, if l is a hyperarc (see "HyperArcID"), each branch can have a
  *   different Efficiency: say, an hyperarc with branches 1 --> 2 with
  *   Efficiency 0.5 and 1 --> 3 with Efficiency 0.5 means that one unit of
  *   flow leaves 1 and half of it reaches 2 while the other half reaches 3.
  *   There is no requirements that the efficiencies of the different branches
  *   of the same hyperarc sum to 1: in fact, this variable is optional, if it
  *   is not specified than Efficiency[ l ] == 1 for all branches / lines.
  *
  * - The variable "LineName", of type netCDF::NcString() and indexed over
  *   the dimension "NumberLines". Its i-th entry, namely LineName[ i ],
  *   contains the name of the i-th transmission line. This variable is
  *   optional.
  */
 
 virtual void deserialize( const netCDF::NcGroup & group ) override;

/** @} ---------------------------------------------------------------------*/
/*------------ METHODS FOR READING THE DATA OF THE DCNetworkData -----------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the DCNetworkData
 * @{ */

 /// returns the number of lines of the network
 /** Method for returning the number of lines of the network. When
  * get_number_nodes() == 1 (the network is a bus), get_number_lines() == 0
  * (no self-loops are allowed, hence there is no line to be made with a
  * single node).
  */

 Index get_number_lines( void ) const { return( f_number_lines ); }

/*--------------------------------------------------------------------------*/
 /// returns the reference node of the network
 /** Method for returning the reference node of the network. */

 Index get_reference_node( void ) const { return( f_reference_node ); }

 /*--------------------------------------------------------------------------*/
 /// returns true if the network is a hypergraph
 /** Method for returning true if the network is a hypergraph, i.e., if it has
  * at least one line with multiple head buses. When is_hypergraph() == false
  * the network is a "regular graph" and therefore get_end_line() has to be
  * used, while if is_hypergraph() == true the network is a hypergraph and
  * therefore get_end_lines() has to be used.
  */

 bool is_hypergraph( void ) const {
  return( f_number_branches > f_number_lines );
  }

/*--------------------------------------------------------------------------*/
 /// returns true if line \p is an hyperarc (more than one head bus)

 bool is_hyperarc( Index line ) const {
  if( is_hypergraph() )
   return( v_end_lines[ line ].size() > 1 );
  else
   return( false );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of start buses for all lines
 /** Method for returning the vector of starting point of each line. This
  * vector may have empty size (bus network) or the size of number of lines,
  * then there are two possible cases:
  *
  * - if get_number_nodes() == 1, this vector has empty size which means there
  *   is no line at network (bus network), and this vector is not needed;
  *
  * - if get_number_nodes() > 1, this vector have size of f_number_lines and
  *    get_start_line()[ l ] gives starting (tail) bus of line l.
  */

 const std::vector< Index > & get_start_line( void ) const {
  return( v_start_line );
  }

/*--------------------------------------------------------------------------*/
 /// returns the start bus of line \p line

 Index get_start_line( Index line ) const { return( v_start_line[ line ] ); }

/*--------------------------------------------------------------------------*/
 /// returns the vector of end buses for all lines
 /** Method for returning the vector of ending point of each line. This
  * vector may have empty size (bus network) or the size of number of lines,
  * then there are three possible cases:
  *
  * - if get_number_nodes() == 1, this vector has empty size which means there
  *   is no line at network (bus network), and this vector is not needed to
  *   be defined.
  *
  * - if get_number_nodes() > 1 and get_number_hyperarcs() == 0, then the
  *   network is a "regular graph", this vector have size of f_number_lines,
  *   and  get_end_line()[ l ] gives ending (head) bus of line l.
  *
  * - if get_number_nodes() > 1 and get_number_hyperarcs() > 0, then the
  *   network is a hypergraph and this vector is again empty since
  *   get_end_lines() must be used to get the set of end buses of the lines.
  */

 const std::vector< Index > & get_end_line( void ) const {
  return( v_end_line );
  }

/*--------------------------------------------------------------------------*/
 /// returns the end (first, in the hypergraph case) bus of line \p line 

 Index get_end_line( Index line ) const {
  if( is_hypergraph() )
   return( v_end_lines[ line ].front() );   
  else
   return( v_end_line[ line ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of (sets of) end lines
 /** Method for returning the vector of sets of ending point of each line.
  * This vector is empty if get_number_hyperarcs() == 0, i.e., the network is
  * "regular graph" (which is true in particular if get_number_nodes() == 1),
  * otherwise  get_end_lines()[ l ] is a (const) std::vector< Index >
  * containing the end buses / nodes of line l. Line l is a "regular arc" if
  * get_end_lines()[ l ].size() == 1, and an hyperarc if
  * get_end_lines()[ l ].size() > 1 (it cannot obviously be 0). The number of
  * lines l such that get_end_lines()[ l ].size() > 1 is equal to
  * get_number_hyperarcs(). Each std::vector< Index > is ordered in increasing
  * sense and without repeated elements.
  */

 const std::vector< std::vector< Index > > & get_end_lines( void ) const {
  return( v_end_lines );
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
  * - if f_number_lines >= 1, this vector have size of f_number_lines and
  *   each element of the vector gives minimum power flow of each line in
  *   the network.
  */

 const std::vector< double > & get_min_power_flow( void ) const {
  return( v_min_power_flow );
  }

/*--------------------------------------------------------------------------*/
 /// returns minimum power flow of the given \p line
 /** This method returns the minimum power flow of the given \p line.
  *
  * @return the minimum power flow of the given \p line.
  */

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
  * - if f_number_lines == 0, this vector has empty size which means there
  *   is no line at network (bus network).
  *
  * - if f_number_lines >= 1, this vector have size of f_number_lines and
  *   each element of the vector gives maximum power flow of each line in
  *   the network.
  */

 const std::vector< double > & get_max_power_flow( void ) const {
  return( v_max_power_flow );
  }

/*--------------------------------------------------------------------------*/
 /// returns maximum power flow of the given \p line
 /** This method returns the maximum power flow of the given \p line.
  *
  * @return the maximum power flow of the given \p line.
  */

 double get_max_power_flow( Index line ) const {
  assert( line < get_number_lines() );
  if( v_max_power_flow.empty() )
   return( 0 );
  return( v_max_power_flow[ line ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns vector of the susceptances
 /** Method for returning the vector of susceptances for each line. This
  * vector may have empty size (bus network) or the size of number of nodes,
  * then there are two possible cases:
  *
  * - if f_number_lines == 0, this vector has empty size which means there is
  *   no line at network (bus network).
  *
  * - if f_number_lines >= 1, this vector has size of f_number_lines and each
  *   element of the vectors gives the Susceptance value for each line in the
  *   network.
  */

 const std::vector< double > & get_line_susceptance( void ) const {
  return( v_line_susceptance );
  }

/*--------------------------------------------------------------------------*/
 /// returns vector of the network cost
 /** Method for returning the vector of network cost for each line. This
  * vector may have empty size (bus network) or the size of number of lines,
  * then there are two possible cases:
  *
  * - if f_number_lines == 0, this vector has empty size which means there is
  *   no line at network (bus network).
  *
  * - if f_number_lines >= 1, this vector has size of f_number_lines and each
  *   element of the vectors gives the network cost value for each line in the
  *   network.
  */

 const std::vector< double > & get_network_cost( void ) const {
  return( v_network_cost );
  }

/*--------------------------------------------------------------------------*/
 /// returns the types of lines in the network
 /** This method returns the types of lines present in the network. */

 line_type get_lines_type( void ) {
  if( f_lines_type < 0 ) {
   if( get_number_lines() == 0 )
    f_lines_type = kNone;
   else
    if( std::all_of( v_line_susceptance.cbegin() , v_line_susceptance.cend() ,
		     []( double s ) { return( s == 0.0 ); } ) )
     f_lines_type = kHVDC;
    else
     if( std::all_of( v_line_susceptance.cbegin() , v_line_susceptance.cend() ,
		      []( double s ) { return( s != 0.0 ); } ) )
      f_lines_type = kAC;
     else
      f_lines_type = kAC_HVDC;
   }
  return( line_type( f_lines_type ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the efficiency of \p line (1 if not specified or not a HVDC line)

 double get_line_efficiency( Index line ) const {
  assert( line < get_number_lines() );
  if( is_hypergraph() )
   throw( std::logic_error(
		      "get_line_efficiency() called but hypergraph" ) );
  if( v_efficiency.empty() || get_line_susceptance().empty() ||
      ( get_line_susceptance()[ line ] != 0.0 ) )
   return( 1.0 );
  return( v_efficiency[ line ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the set of efficiencies for all heads of hyperline \p line

 const std::vector< double > & get_line_efficiencies( Index line ) const {
  if( ! is_hypergraph() )
   throw( std::logic_error(
		      "get_line_efficiencies() called but no hypergraph" ) );

  return( v_h_efficiency[ line ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector containing the name of the lines

 const std::vector< std::string > & get_line_names( void ) const {
  return( v_line_names );
  }

/** @} ---------------------------------------------------------------------*/
/*-------------------- METHODS FOR SAVING THE DCNetworkData ----------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the DCNetworkData
 * @{ */

 /// serialize a DCNetworkData out of a netCDF::NcGroup
 /** Serialize a DCNetworkData out of a netCDF::NcGroup to the specific
  * format of a DCNetworkData. See NetworkBlock::deserialize( netCDF::NcGroup
  * ) for details of the format of the created netCDF group.
  */

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

 Index f_number_lines;      ///< number of lines of the network

 Index f_reference_node;    ///< reference node (used in the PTDF matrix)

 int f_lines_type;          ///< the type of the network
 
 Index f_number_branches;  ///< the number of branches of all hyperarcs

 Subset v_start_line;      ///< vector of starting lines

 Subset v_end_line;        ///< vector of ending lines

 std::vector< Subset > v_end_lines;
 ///< vector of (vector of) sets of ending lines for hyperarcs
 
 /// vector to store the susceptance of each line of the network
 std::vector< double > v_line_susceptance;

 /// vector to store the minimum power flow at each line
 std::vector< double > v_min_power_flow;

 /// vector to store the maximum power flow at each line
 std::vector< double > v_max_power_flow;

 /// vector to store the network cost at each line
 std::vector< double > v_network_cost;

 /** vector to store the network efficiency of each line in the graph case,
  * effective only for HVDC lines and ignored otherwise */
 std::vector< double > v_efficiency;

 /** vector to store the network efficiency of each (hyper)line in the
  * hypergraph case, effective only for HVDC lines and ignored otherwise */
 std::vector< std::vector< double > > v_h_efficiency;

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
 : NetworkBlock( f_block ) , f_NetworkData( nullptr ) ,
   f_InvestmentCost( 0 ), f_MinCapacityDesign( 0 ), f_MaxCapacityDesign( 1 ) {}

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
  * - if the NetworkData object is present (either in the NcGroup or because
  *   it has been passed and NumberNodes > 1) then this variable has size
  *   "NumberNodes", which can be read via NetworkData::get_number_nodes().
  *
  * Depending on the susceptance for each line of the network, the
  * DCNetworkBlock class may have a power flow variable or not. In other
  * word, if the susceptance is equal to zero (or not defined), the
  * corresponding line is a HVDC line and it must have the power flow
  * variable. It means, each HVDC line correspond to a power flow variable,
  * then for the Net Transfer Capacity (NTC) model all lines must have a
  * power flow variable. If the susceptance value is a non-zero value, the
  * corresponding line is called AC and there is no needed to define the
  * power flow variable for that line. Therefore, in the case of pure AC line
  * there is no needed to define power flow variables. Consequently, for the
  * mixed case AC-HVDC, the power flow variable must define just for HVDC
  * lines. Similarly, depending on the NetworkCost for each line of the
  * network, the DCNetworkBlock class may have an auxiliary variable or not.
  * In other word, if the NetworkCost is equal to zero (or not defined), the
  * auxiliary variable and corresponding constraints will not be defined.
  *
  * In the design scenario of the UC problem (i.e., when an investment cost
  * is provided), an additional design variable \f$ x \f$ is created. Its
  * type depends on \f$ \mathrm{MaxCapacityDesign} \f$: if
  * \f$ \mathrm{MaxCapacityDesign} < 0 \f$ then \f$ x \f$ is binary; otherwise
  * \f$ x \f$ is nonnegative continuous and bounded by
  * \f$ 0 \le x \le \mathrm{MaxCapacityDesign} \f$.
  *
  * In addition, when \( \mathrm{MaxCapacityDesign} \ge 0 \) a lower bound
  * \( \mathrm{MinCapacityDesign} \) may be provided, yielding
  * \( \mathrm{MinCapacityDesign} \le x \le \mathrm{MaxCapacityDesign} \).
  * When \( \mathrm{MaxCapacityDesign} < 0 \) (binary design), \( x \in \{0,1\} \);
  * if \( \mathrm{MinCapacityDesign} > 0 \), then \( x \) is effectively forced to 1.
  */

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

/*--------------------------------------------------------------------------*/
 ///generate PTDF matrix of DCNetworkBlock
 /** The function constructs an Eigen-type PTDF (Power Transfer Distribution 
  * Factor) matrix of the AC part of the Network.
  * Denoting by \f$\mathfrak{S}_l\f$ the susceptance of the line \f$l\f$, the 
  * PTDF matrix \f$ B \f$ is computed as 
  * \f$ B = (\hat{B}I_{n_0})(I_{n_0}^T\bar{B}I_{n_0})^{-1} \f$, where
  * \f$ n_0 \f$ is the reference node of the network. The matrix
  * \f$ \hat{B} \f$ is defined as
  * \f$ (\hat{B})_{(l=(n,n'),n)} = \mathfrak{S}_l \f$, 
  * \f$ (\hat{B})_{(l=(n',n),n)} = -\mathfrak{S}_l \f$ and \f$ 0 \f$
  * otherwise.
  *    
  * The matrix \f$ \bar{B} \f$ is defined as 
  * 
  * \f[
  *   (\bar{B})_{(n,n')} = \begin{cases}
  *      \mathfrak{S}_{l=(n,n')} \textnormal{ if } (n,n')\in L\	\
  *      \mathfrak{S}_{l=(n,n')} \textnormal{ if } (n',n)\in L\	\
  *      \sum_{l=(n,\cdot)\in L} \mathfrak{S}_{l} 
  *            + \sum_{l=(\cdot,n)\in L} \mathfrak{S}_{l} 
  *      \textnormal{ if } n = n'               \quad n,n' \in N
  *  \f]
  */

 Eigen::MatrixXd get_PTDF( const std::vector< Index > & AC_lines ) const;

 Eigen::MatrixXd get_PTDF() const {
  std::vector< Index > all_lines( get_number_lines() );
  std::iota( all_lines.begin() , all_lines.end() , 0 );
  return( get_PTDF( all_lines ) );
  }

/*--------------------------------------------------------------------------*/

 int get_reducedIdx( int idx ) const;

/*--------------------------------------------------------------------------*/
 /// generate the abstract constraints of the DCNetworkBlock
 /** Method that generates the abstract constraints of the DCNetworkBlock.
  * These are:
  *
  * - **flow bounds for HVDC lines** (NTC model). Each HVDC line
  *   \f$ l \in \mathcal{L} \f$ has controllable flow variable \f$ F_l \f$
  *   with bounds scaled by \f$ \kappa_l \f$:
  *
  *   \f[
  *     \kappa_l P^{mn}_l \;\le\; F_l \;\le\; \kappa_l P^{mx}_l
  *         \quad l \in \mathcal{L} \quad (1a)
  *   \f]
  *
  *   In the design scenario these become:
  *
  *   \f[
  *     x \, ( \kappa_l P^{mn}_l ) \;\le\; F_l \;\le\;
  *     x \, ( \kappa_l P^{mx}_l ) \quad l \in \mathcal{L} \quad (1b)
  *   \f]
  *
  *   where \f$ x \f$ is the design variable.
  *
  * - **nodal balance equations** for active power injections:
  *
  *   \f[
  *     \sum_{l=(n,\cdot)} F_l \;-\; \sum_{l=(\cdot,n)} F_l
  *       \;=\; S_n - D^{ac}_n
  *       \quad n \in \mathcal{N} \quad (2)
  *   \f]
  *
  *   with \f$ S_n \f$ the injection variable and \f$ D^{ac}_n \f$ the demand.
  *
  * - **absolute-value linearization** of HVDC flows when network
  *   costs are active:
  *
  *   \f[
  *     F_l \le V_l \quad l \in \mathcal{L} \quad (3a)
  *   \f]
  *
  *   \f[
  *     -V_l \le F_l \quad l \in \mathcal{L} \quad (3b)
  *   \f]
  *
  *   where \f$ V_l \f$ is the auxiliary variable.
  *
  * - **flow bounds for AC lines.** With PTDF matrix
  *   \f$ B \in \mathbb{R}^{|\mathcal{L}|\times|\mathcal{N}|} \f$:
  *
  *   \f[
  *     P^{mn}_l \;\le\; \sum_{n \in \mathcal{N}} B_{(l,n)} (S_n - D^{ac}_n)
  *      \;\le\; P^{mx}_l \quad l \in \mathcal{L} \quad (4a)
  *   \f]
  *
  *   In the design scenario:
  *
  *   \f[
  *     x \, P^{mn}_l \;\le\; \sum_{n \in \mathcal{N}} B_{(l,n)} (S_n - D^{ac}_n)
  *      \;\le\; x \, P^{mx}_l \quad l \in \mathcal{L} \quad (4b)
  *   \f]
  *
  * - **hybrid AC/HVDC network.** Let \f$ \mathcal{L}^{ac} \f$ and
  *   \f$ \mathcal{L}^{dc} \f$ denote AC and HVDC lines, and define
  *   incidence \f$ A^{dc} \f$ and PTDF \f$ B \f$. Then
  *
  *   \f[
  *     A =
  *     \begin{bmatrix}
  *       B & - B (A^{dc})^\top \\
  *       0 & I
  *     \end{bmatrix} ,
  *   \f]
  *
  *   and flows satisfy
  *
  *   \f[
  *     P^{mn} \;\le\; A \begin{bmatrix} a \\ b \end{bmatrix}
  *      \;\le\; P^{mx} \quad (5a)
  *   \f]
  *
  *   with injections \f$ a_n = S_n - D^{ac}_n \f$ and HVDC flows
  *   \f$ b_m = F_{m+|\mathcal{L}^{ac}|} \f$. In design mode:
  *
  *   \f[
  *     x \, P^{mn} \;\le\; A \begin{bmatrix} a \\ b \end{bmatrix}
  *      \;\le\; x \, P^{mx} \quad (5b)
  *   \f]
  *
  * - **overall balance constraint** ensuring that total injection equals
  *   total demand:
  *
  *   \f[
  *     \sum_{n \in \mathcal{N}} S_n \;=\;
  *     \sum_{n \in \mathcal{N}} D^{ac}_n \quad (6)
  *   \f]
  *
  * - **design bounds** on \f$ x \f$:
  *
  *   \[
  *     x \in
  *     \begin{cases}
  *       \{0,1\} & \text{if } \mathrm{MaxCapacityDesign} < 0 \\
  *       [\,\mathrm{MinCapacityDesign},\,\mathrm{MaxCapacityDesign}\,]
  *         & \text{if } \mathrm{MaxCapacityDesign} \ge 0
  *     \end{cases}
  *   \]
  *
  *   with the convention that if \f$ \mathrm{MinCapacityDesign} > 0 \f$
  *   and binary design, then \f$ x = 1 \f$.
  */

 void generate_abstract_constraints( Configuration * stcc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the objective of the DCNetworkBlock
 /** Method that generates the objective of the DCNetworkBlock. The
  * objective can include:
  *
  * - a linear term on the auxiliary variables associated with network
  *   costs, if the vector "NetworkCost" is provided:
  *   \f[
  *     \min \ \sum_{l \in \mathcal{L}} NC_l \cdot V_l
  *   \f]
  *   where \f$ NC_l \f$ is the unit network cost of line \f$l\f$ and
  *   \f$ V_l \f$ is the corresponding auxiliary variable
  *   (coefficients are also scaled by the Block scale factor, if any);
  *
  * - an investment term in design mode, if "InvestmentCost" \f$ \ne 0 \f$:
  *   \f$ + \ I \cdot x \f$, where \f$I\f$ is the investment cost and
  *   \f$x\f$ is the design variable.
  *
  * Hence, in the design scenario the full objective is:
  * \f[
  *   \min \ \sum_{l \in \mathcal{L}} NC_l \cdot V_l \;+\; I \cdot x \; .
  * \f]
  * If "NetworkCost" is not provided, \f$ NC_l = 0 \f$ and only the
  * investment term remains in design mode.
  */
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
  * satisfies a Constraint, the Constraint is computed
  * (Constraint::compute()).
  *
  * @param useabstract This parameter is currently ignored.
  *
  * @param fsbc The pointer to a Configuration that specifies the tolerance
  *             and the type of violation that must be considered.
  */

 bool is_feasible( bool useabstract = false ,
                   Configuration * fsbc = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
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
  * @return the number of nodes in the network.
  */

 Index get_number_nodes( void ) const override {
  if( ! f_NetworkData )
   return( 1 );
  return( f_NetworkData->get_number_nodes() );
  }

/*--------------------------------------------------------------------------*/
 /// returns the AC lines
 /** This function returns the AC lines in the transmission network.
  *
  * @return the AC lines in the network.
  */

 std::vector< Index > get_AC_lines( void ) const {
  std::vector< Index > AC_lines;
  const auto & susceptance = f_NetworkData->get_line_susceptance();
  if( AC_lines.empty() )
   return( AC_lines );
  for( Index line_id = 0 ; line_id < f_NetworkData->get_number_lines() ;
       ++line_id ) {
   if( susceptance[ line_id ] > 0. )
    AC_lines.push_back( line_id );
   }
  return( AC_lines );
  }

/*--------------------------------------------------------------------------*/
 /// returns the DC lines
 /** This function returns the DC lines in the transmission network.
  *
  * @return the DC lines in the network.
  */

 std::vector< Index > get_DC_lines( void ) const {
  std::vector< Index > DC_lines;
  const auto & susceptance = f_NetworkData->get_line_susceptance();
  for( Index line_id = 0 ; line_id < f_NetworkData->get_number_lines() ;
       ++line_id ) {
   if( ( susceptance.empty() ) || ( susceptance[ line_id ] == 0. ) )
    DC_lines.push_back( line_id );
   }
  return( DC_lines );
  }

/*--------------------------------------------------------------------------*/
 /// returns the number of lines of the network
 /** This function returns the number of lines in the transmission network.
  * If get_NetworkData() returns nullptr, this is equivalent to
  * get_NetworkData()->get_number_lines(). Otherwise, it returns zero.
  *
  * @return the number of lines in the network.
  */

 Index get_number_lines( void ) const {
  if( ! f_NetworkData )
   return( 0 );
  return( f_NetworkData->get_number_lines() );
  }

/*--------------------------------------------------------------------------*/
 /// returns the types of lines in the network
 /** This method returns the types of lines present in the network. */

 line_type get_lines_type( void ) const {
  if( ! f_NetworkData )
   return( kNone );
  return( f_NetworkData->get_lines_type() );
  }

/*--------------------------------------------------------------------------*/
 /// returns the kappa constant associated with the given \p line
 /** This function returns the kappa constant associated with the given \p
  * line. This is the constant that multiplies the minimum and maximum flow in
  * the flow limit constraint associated with the given \p line.
  *
  * @param line The index of a line (between 0 and get_number_lines() - 1).
  *
  * @return The kappa constant associated with the given \p line.
  */

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
  * @return The minimum power flow on the given \p line.
  */

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
  * @return The maximum power flow on the given \p line.
  */

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
  *                 returned.
  */

 const double * get_active_demand( Index interval = 0 ) const override {
  if( v_ActiveDemand.empty() )
   return( nullptr );
  return( &( v_ActiveDemand.front() ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the investment cost
 double get_investment_cost( void ) const { return( f_InvestmentCost ); }

/** @} ---------------------------------------------------------------------*/
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
  *   variable for line l.
  */

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
  *   variable for line l.
  */

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
  * i-th line of the network.
  */

 const std::vector< FRowConstraint > &
 get_power_flow_limit_constraints( void ) const {
  if( ! f_NetworkData )
   throw( std::logic_error(
			 "DCNetworkBlock::get_power_flow_limit_constraints:"
			 " DCNetworkData has not been set" ) );

  switch( f_NetworkData->get_lines_type() ) {
   case( kAC ):      return( v_AC_power_flow_limit_const );
   case( kAC_HVDC ):
   default:          return( v_AC_HVDC_power_flow_limit_const );
   }
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of power flow limit HVDC bounds

 const std::vector< BoxConstraint > &
 get_power_flow_limit_HVDC_bounds( void ) const {
  if( ! f_NetworkData )
   throw( std::logic_error(
			"DCNetworkBlock::get_power_flow_limit_HVDC_bounds:"
			" DCNetworkData has not been set" ) );
  return( v_HVDC_power_flow_limit_const );
  }

/*--------------------------------------------------------------------------*/
 /// returns the dual prices of power flow limits, however the network is

 void get_dual_prices( std::vector< double > & dp ) const {
  auto nl = get_number_lines();
  if( ! nl ) {
   dp.clear();
   return;
   }

  dp.resize( nl );
  auto lt = f_NetworkData->get_lines_type();
  switch( lt ) {
   case( kHVDC ):
    for( Index l = 0 ; l < nl ; ++l )
     dp[ l ] = v_HVDC_power_flow_limit_const[ l ].get_dual();
    break;
   case( kAC ):
    for( Index l = 0 ; l < nl ; ++l )
     dp[ l ] = v_AC_power_flow_limit_const[ l ].get_dual();
    break;
   case( kAC_HVDC ):
    for( Index l = 0 ; l < nl ; ++l )
     dp[ l ] = v_AC_HVDC_power_flow_limit_const[ l ].get_dual();
    break;
   default:
    throw( std::logic_error( "unknown or unhandled get_lines_type()" ) );
   }
  }
 
/** @} ---------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Solution
 * @{ */

 /// returns a Solution representing the current solution of this NetworkBlock
 /** This method must construct and return a (pointer to a) Solution object
  * representing the current "solution state" of this NetworkBlock. This is
  * a DCNetworkBlockSolution extending NetworkBlockSolution with the specific
  * extra solution information of DCNetworkBlock.
  *
  * The parameter for deciding which kind of Solution must be returned is a
  * single int value, coded bitwise:
  *
  * - bit 0 (& 1) is "taken" by the base :NetworkBlock[Solution]
  *
  * - bit 1 (& 2) means "store the flow values"
  *
  * - bit 2 (& 4) means "store the dual prices"
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
  * - otherwise, it is 7 (save everything).
  */

 Solution * get_Solution( Configuration * solc = nullptr ,
                          bool emptys = true ) override;

/*--------------------------------------------------------------------------*/
 /// return the "appropriate" [DC]NetworkBlockSolution
 
 NetworkBlockSolution * new_Solution( void ) const override;

/** @} ---------------------------------------------------------------------*/
/*--------------- METHODS FOR MODIFYING THE DCNetworkBlock -----------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the DCNetworkBlock
 * @{ */

 void set_NetworkData( NetworkData * nd = nullptr ) override {
  // if there was a previous DCNetworkData, and it was local, delete it
  if( f_NetworkData && f_local_NetworkData )
   delete( f_NetworkData );

  f_NetworkData = dynamic_cast< DCNetworkData * >( nd );
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
  * the NcGroup then this data is used, replacing (and therefore ignoring)
  * the data set by this method.
  *
  * Similarly, if this method is called *after* deserialize(), but some the
  * ActiveDemand was already present in the NcGroup, then that data is kept
  * and the call to this method does nothing.
  *
  * When this method is called, if it is empty it is written into, otherwise
  * nothing happens. In deserialize(), if the data is there in the NcGroup 
  *then it is written in v_ActiveDemand (which therefore is no longer empty),
  * otherwise it is left empty so that it can be set by this method.
  */

 void set_ActiveDemand( const boost::multi_array< double , 2 > & v )
  override {
  if( v_ActiveDemand.empty() )
   v_ActiveDemand.assign( v[ 0 ].begin() , v[ 0 ].end() );
  }

/*--------------------------------------------------------------------------*/
 /// sets the power flows

 void set_power_flow( const std::vector< double > & pf ) {
  for( Index l = 0 ; l < v_power_flow.size() ; ++l )
   v_power_flow[ l ].set_value( pf[ l ] );
  }
  
/*--------------------------------------------------------------------------*/
 /// sets the dual prices of power flow limits, however the network is

 void set_dual_prices( const std::vector< double > & dp ) {
  auto nl = get_number_lines();
  if( ! nl )
   return;

  auto lt = f_NetworkData->get_lines_type();
  switch( lt ) {
   case( kHVDC ):
    for( Index l = 0 ; l < nl ; ++l )
     v_HVDC_power_flow_limit_const[ l ].set_dual( dp[ l ] );
    break;
   case( kAC ):
    for( Index l = 0 ; l < nl ; ++l )
     v_AC_power_flow_limit_const[ l ].set_dual( dp[ l ] );
    break;
   case( kAC_HVDC ):
    for( Index l = 0 ; l < nl ; ++l )
     v_AC_HVDC_power_flow_limit_const[ l ].set_dual( dp[ l ] );
    break;
   default:
    throw( std::logic_error( "unknown or unhandled get_lines_type()" ) );
   }
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
  * - The scalar variable "InvestmentCost", of type netCDF::NcDouble and not
  *   indexed over any dimension. When provided and different from 0, the
  *   model enters the design scenario and a design variable \f$ x \f$ is
  *   generated.
  *
  * - The scalar variable "MinCapacityDesign", of type netCDF::NcDouble and
  *   not indexed over any dimension. This sets the lower bound of the design
  *   variable \( x \) in design mode (i.e., when InvestmentCost != 0). If not
  *   provided, the default is 0. Its meaning depends on "MaxCapacityDesign":
  *   - if \( \mathrm{MaxCapacityDesign} < 0 \) (binary design), then
  *     \( x \in \{0,1\} \) and \( \mathrm{MinCapacityDesign} > 0 \) implies
  *     \( x = 1 \);
  *   - otherwise (continuous design), \( x \) is nonnegative continuous with
  *     \( \mathrm{MinCapacityDesign} \le x \le \mathrm{MaxCapacityDesign} \).
  *
  * - The scalar variable "MaxCapacityDesign", of type netCDF::NcDouble and
  *   not indexed over any dimension. This limits the design variable \( x \):
  *   - if \( \mathrm{MaxCapacityDesign} < 0 \) then \( x \in \{0,1\} \) (binary);
  *   - if \( \mathrm{MaxCapacityDesign} = 1 \) then \( x \in [0,1] \) when
  *     \( \mathrm{MinCapacityDesign} = 0 \); otherwise
  *     \( x \in [\,\mathrm{MinCapacityDesign},\,1] \);
  *   - if \( \mathrm{MaxCapacityDesign} > 0 \) then \( x \) is nonnegative
  *     continuous with \( \mathrm{MinCapacityDesign} \le x \le \mathrm{MaxCapacityDesign} \).
  *   If not provided, the default is 1.
  *
  * - The variable "ActiveDemand", of type netCDF::NcDouble and indexed over
  *   the dimension "NumberNodes". If the NetworkData object description is
  *   present in the NcGroup this is the dimension "NumberNodes", but the
  *   NetworkData object is optional and it may not be there. Thus, if
  *   "NumberNodes" is not there and "ActiveDemand" is, then the NetworkData
  *   object must have been passed by set_NetworkData(), and the number of
  *   nodes can be read via NetworkData::get_number_nodes(). However,
  *   "ActiveDemand" itself is optional. If it is not found in the NcGroup,
  *   then it *must* be passed (either before or after the call to
  *   deserialize()) by calling set_active_demand(). Since both groups of data
  *   are optional, the NcGroup  can actually be empty which implies that all
  *   the data will be (or have been) passed by the in-memory interface. In
  *   this case, it would clearly be preferable to *entirely avoid the
  *   NcGroup to be there*, and in fact UCBlock has provisions for the
  *   NcGroup describing the NetworkBlock to be optional [see the comments to
  *   UCBlock::deserialize()].
  *
  * - The variable "Kappa", of type netCDF::NcDouble and either being a
  *   scalar or indexed over the number of lines. If this variable is a
  *   scalar, let say k, then it is assumed that Kappa[ l ] = k for each line
  *   l in {0, ..., get_number_lines() - 1}. For each line l in {0, ...,
  *   get_number_lines() - 1}, Kappa[ l ] is the constant that multiplies the
  *   minimum and maximum flow in the flow limit constraints. This variable is
  *   optional. If it is not provided, it is assumed that Kappa[ l ] == 1 for
  *   each line l in {0, ..., get_number_lines() - 1}.
  *
  * - The variable "ConstantTerm", of type netCDF::NcDouble and containing the
  *   constant term.
  */

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
 /// loads the DCNetworkBlock instance from memory
 /** Like load( std::istream & ), if there is any Solver attached to this
  * DCNetworkBlock then a NBModification (the "nuclear option") is issued.
  */

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
  * details of the format of the created netCDF group.
  */

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
  *               modified.
  *
  * @param ordered It indicates whether \p subset is ordered.
  *
  * @param issuePMod It controls how physical Modification are issued.
  *
  * @param issueAMod It controls how abstract Modification are issued.
  */

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
  *            constants are being modified.
  *
  * @param issuePMod It controls how physical Modification are issued.
  *
  * @param issueAMod It controls how abstract Modification are issued.
  */

 void set_kappa( MF_dbl_it values ,
                 Range rng = Range( 0 , Inf< Index >() ) ,
                 c_ModParam issuePMod = eNoBlck ,
                 c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/

 /// change the abstract representation of the power flow limit constraints
 /** This function changes the abstract representation of the power flow limit 
  * constraints for indices in \p modified_lines.
  *
  * @param modified_lines A vector of the indices of constraints that 
  *                       have to be modified.
  *
  * @param issueAMod It controls how abstract Modification are issued.
  */

 void change_power_flow_limit_constraints( 
                const std::vector< Index > & modified_lines,
                c_ModParam issueAMod);

/*--------------------------------------------------------------------------*/

 /// change the abstract representation of the constraints on the auxiliary variables
 /** This function changes the abstract representation of the constraints on the
  * auxiliary variables for indices in \p modified_lines.
  *
  * @param modified_lines A vector of the indices of constraints that 
  *                       have to be modified.
  *
  * @param issueAMod It controls how abstract Modification are issued.
  */

 void change_relax_abs_constraints(
        const std::vector< Index > & modified_lines , c_ModParam issueAMod );

/*--------------------------------------------------------------------------*/
 /// change the abstract representation of the injection constraints
 /** This function changes the abstract representation of the power flow
  * injection constraints for indices in \p modified_nodes.
  *
  * @param modified_nodes A vector of the indices of nodes that 
  *                       have to be modified.
  *
  * @param issueAMod It controls how abstract Modification are issued.
  */

 void change_DC_power_flow_injection_constraints(
                const std::vector< Index > & modified_nodes ,
                c_ModParam issueAMod );

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
  *               modified.
  *
  * @param ordered It indicates whether \p subset is ordered.
  *
  * @param issuePMod It controls how physical Modification are issued.
  *
  * @param issueAMod It controls how abstract Modification are issued.
  */

 void set_active_demand( MF_dbl_it values , Subset && subset ,
                         bool ordered = false ,
                         c_ModParam issuePMod = eNoBlck ,
                         c_ModParam issueAMod = eNoBlck ) override;

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
  *            demand is being modified.
  *
  * @param issuePMod It controls how physical Modification are issued.
  *
  * @param issueAMod It controls how abstract Modification are issued.
  */

 void set_active_demand( MF_dbl_it values ,
                         Range rng = Range( 0 , Inf< Index >() ) ,
                         c_ModParam issuePMod = eNoBlck ,
                         c_ModParam issueAMod = eNoBlck ) override;

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

 DCNetworkData * f_NetworkData;  ///< the DCNetworkData object

 /// the investment cost
 double f_InvestmentCost;

 /// the minimum capacity design allowed
 double f_MinCapacityDesign;

 /// the maximum capacity design allowed
 double f_MaxCapacityDesign;

 /// vector to store the demand of each node of the network
 std::vector< double > v_ActiveDemand;

 /// the kappa constant for each line
 std::vector< double > v_kappa;

/*-------------------------------- variables -------------------------------*/

 /// the power flow variables
 std::vector< ColVariable > v_power_flow;

 /// the auxiliary network cost variable
 std::vector< ColVariable > v_auxiliary_variable;

 /// the design variable
 ColVariable design;

/*------------------------------- constraints ------------------------------*/

 /// AC power flow limit constraints
 std::vector< FRowConstraint > v_AC_power_flow_limit_const;

 /// AC power flow bounds design constraints
 boost::multi_array< FRowConstraint , 2 > v_AC_power_flow_bounds_design_const;


 /// AC/HVDC power flow limit constraints
 std::vector< FRowConstraint > v_AC_HVDC_power_flow_limit_const;

 /// AC/HVDC power flow bounds design constraints
 boost::multi_array< FRowConstraint , 2 > v_AC_HVDC_power_flow_bounds_design_const;


 /// HVDC power flow limit constraints
 std::vector< BoxConstraint > v_HVDC_power_flow_limit_const;

 /// HVDC power flow bounds design constraints
 boost::multi_array< FRowConstraint , 2 > v_HVDC_power_flow_bounds_design_const;


 /// HVDC power flow and node injection constraints
 std::vector< FRowConstraint > v_power_flow_injection_const;

 /// HVDC power flow auxiliary variable constraints
 boost::multi_array< FRowConstraint , 2 > v_power_flow_relax_abs;

 /// injection equals to demand
 FRowConstraint overall_balanced_const;

 /// the design bound constraint
 BoxConstraint design_bound_const;

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

 /// verify whether the data in this DCNetworkBlock is consistent
 /** This function checks whether the data in this DCNetworkBlock is
  * consistent. The data is consistent if all the following conditions are met.
  *
  * - Design bounds consistency:
  *   - \( \mathrm{MinCapacityDesign} \ge 0 \);
  *   - if \( \mathrm{MaxCapacityDesign} > 0 \), then
  *     \( \mathrm{MinCapacityDesign} \le \mathrm{MaxCapacityDesign} \);
  *   - if \( |\mathrm{MaxCapacityDesign}| = 1 \), then
  *     \( \mathrm{MinCapacityDesign} \le 1 \);
  *   - if \( \mathrm{MaxCapacityDesign} < 0 \) (binary), then
  *     \( \mathrm{MinCapacityDesign} \le 1 \)
  *     (note: \( \mathrm{MinCapacityDesign} > 0 \Rightarrow x = 1 \)).
  *
  * If any of the above conditions are not met, an exception is thrown.
  */
 void check_data_consistency( void ) const;

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
 DCNetworkBlockRngdMod( DCNetworkBlock * const fblock , int type ,
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
 DCNetworkBlockSbstMod( DCNetworkBlock * const fblock , int type ,
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
/*-------------------- CLASS DCNetworkBlockSolution ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// a [NetworkBlock]Solution of a DCNetworkBlock
/** The DCNetworkBlockSolution class derives from NetworkBlockSolution and
 * adds the "standard" information stored in there (the node injection
 * variables) the other information that is typical of the DCNetworkBlock,
 * i.e.,
 *
 * - the flow variables on each link
 *
 * - [if available] the dual prices of the link capacity constraints; since
 *   these are typically interpreted as costs and the sign depends on
 *   whether the "upper" or "lower" capacity is active, but the orientation
 *   of links is arbitrary, the absolute value of the reduced cost of the
 *   corresponding constraints is returned
 *
 * Note that one DCNetworkBlock covers one time instant, so these variables
 * do not need to be indexed over time instants (unlike those of the base
 * NetworkBlockSolution).
 */

class DCNetworkBlockSolution : public NetworkBlockSolution
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*------------------------------- FRIENDS ----------------------------------*/

 friend DCNetworkBlock;  ///< make DCNetworkBlock friend

/*---------- CONSTRUCTING AND DESTRUCTING DCNetworkBlockSolution -----------*/

 explicit DCNetworkBlockSolution( void ) : NetworkBlockSolution() ,
  f_number_lines( 0 ) { }  ///< constructor, it has nothing to do

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// deserialize a DCNetworkBlockSolution from a netCDF::NcGroup

 void deserialize( const netCDF::NcGroup & group ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// deserialize a DCNetworkBlockSolution from a "global" netCDF::NcGroup

 void deserialize( const netCDF::NcGroup & group , size_t idx ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 ~DCNetworkBlockSolution() = default;
 ///< destructor: it is virtual, and empty

/*------ METHODS DESCRIBING THE BEHAVIOR OF A DCNetworkBlockSolution ------*/

 void read( const Block * block ) override final;

 void write( Block * block ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// serialize a DCNetworkBlockSolution into a netCDF::NcGroup
 /** Serialize a DCNetworkBlockSolution into a netCDF::NcGroup. The format is
  * the one of NetworkBlockSolution, cf. the comments in
  * NetworkBlockSolution::serialize( netCDF::NcGroup & ), except that
  *
  *     "NumberNetworks" IS NOT REALLY NEEDED, BECAUSE "TotalNumberInstants"
  *     AND "EndInstant" ARE NOR REQUIRED SINCE DCNetworkBlock ALWAYS HAS
  *     DCNetworkBlock::get_number_intervals() == 1, AND ALL THE
  *     NetworkBlock IN \p group ARE SUPPOSED TO BE DCNetworkBlock
  *
  * In addition, \p group must contain:
  *
  * - The dimension "NumberLines" containing the number of lines in the
  *   transmission network. It is mandatory. Note that
  *
  *       ALL THE DCNetworkBlock MUST HAVE THE SAME NUMBER OF LINES
  *
  * - The variable "FlowValue", of type netCDF::NcDouble and indexed over
  *   the dimension "NumberLines"; FlowValue[ l ] is the optimal value of
  *   the power flow on line l. The variable is optional.
  *
  * - The variable "DualCost", of type netCDF::NcDouble and indexed over
  *   the dimension "NumberLines"; DualCost[ l ] is the absolute value of
  *   the dual variable of the constraint representing the capacity of
  *   line l. The variable is optional.
  */

 void serialize( netCDF::NcGroup & group ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// serialize a DCNetworkBlockSolution into a "global" netCDF::NcGroup
 /** "nonstandard" version of serialize() that loads a DCNetworkBlockSolution
  * from a "global" netCDF::NcGroup, i.e., one where the solution information
  * of multiple DCNetworkBlock are stored together (to avoid performance
  * issues due to the fact that netCDF is not structured to work with a large
  * number of sub-NcGroup in a file). The format is the  "nonstandard" one of
  * NetworkBlockSolution, cf. the comments in
  * NetworkBlockSolution::serialize( netCDF::NcGroup & , size_t ), except
  * that
  *
  *     "NumberNetworks" IS NOT REALLY NEEDED, BECAUSE "TotalNumberInstants"
  *     AND "EndInstant" ARE NOT REQUIRED SINCE DCNetworkBlock ALWAYS HAS
  *     DCNetworkBlock::get_number_intervals() == 1, AND ALL THE
  *     NetworkBlock IN \p group ARE SUPPOSED TO BE DCNetworkBlock
  *
  * In addition, \p group must contain:
  *
  * - The dimension "NumberLines" containing the number of lines in the
  *   transmission network. It is mandatory. Note that there is only one
  *   copy of the dimension, and as a consequence
  *
  *       ALL THE DCNetworkBlock MUST HAVE THE SAME NUMBER OF LINES
  *
  *   (which is of course necessary since they all take their data from
  *   the same variables where "NumberLines" is one of the dimensions)
  *
  * - The variable "FlowValue", of type netCDF::NcDouble and indexed both
  *   over the dimension "NumberNetworks" (which is the same as
  *   "TotalNumberInstants", that does not exist) and the dimension
  *   "NumberLines"; FlowValue[ idx ][ l ] is the optimal value of
  *   the power flow on line l for this DCNetworkBlock. The variable is
  *   optional.
  *
  * - The variable "DualCost", of type netCDF::NcDouble and indexed both
  *   over the dimension "NumberNetworks" (which is the same as
  *   "TotalNumberInstants", that does not exist) and the dimension
  *   "NumberLines"; DualCost[ idx ][ l ] is the absolute value of the
  *   dual variable of the constraint representing the capacity of line l
  *   for this DCNetworkBlock. The variable is optional.
  *
  * Note that the variable are constructed when \p idx == 0 according to
  * the fact that the corresponding DCNetworkBlockSolution has or not been
  * Configure-d to hold them, which means that
  *
  *       ALL THE DCNetworkBlockSolution MUST HAVE BEEN Configure-d IN THE
  *       SAME WAY
  *
  * (although, technically, if some of the DCNetworkBlockSolution that
  * appears when \p idx > 0 is Configure-d with less information than that
  * when idx == 0 the code will not break, but there will be uninitialised
  * values in the netCDF).
  */

 void serialize( netCDF::NcGroup & group , size_t idx ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 DCNetworkBlockSolution * scale( double factor ) const override;

 void sum( const Solution * solution , double multiplier ) override;

 DCNetworkBlockSolution * clone( bool empty = false ) const override;

/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/

 void print( std::ostream &output ) const override {
  output << "DCNetworkBlockSolution [" << this << "]: " << std::endl;
  }

/*---------------------- PRIVATE PART OF THE CLASS -------------------------*/

 private:

/*---------------------------- PRIVATE FIELDS ------------------------------*/

 Index f_number_lines;          ///< the number of lines

 std::vector< double > v_flow;  ///< v_flow[ l ] = flow variable on line l

 std::vector< double > v_cost;  /**< v_cost[ l ] = absolute value of the
				 *                 reduced cost of the
				 * capacity constraint of line l */

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class( DCNetworkBlockSolution ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 }  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* __DCNetworkBlock */

/*--------------------------------------------------------------------------*/
/*-------------------- End File DCNetworkBlock.h ---------------------------*/
/*--------------------------------------------------------------------------*/
