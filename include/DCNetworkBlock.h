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

#include <Eigen/Sparse>

#include <utility>

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)

namespace SMSpp_di_unipi_it
{
 using SpMat = Eigen::SparseMatrix< double >;

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS DCNetworkBlock ---------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// a "DC" transmission NetworkBlock
/** The DCNetworkBlock class derives from NetworkBlock, and defines the
 * linear constraints corresponding to the "DC model" of the transmission
 * network in the Unit Commitment problem. Generally, there exist
 * three different kinds of DCNetworkBlock:
 *
 * - DCNetworkBlock with just HVDC lines, i.e., where the susceptance for
 *   all lines is equal to zero. It's also known as the Net Transfer
 *   Capacity (NTC) model.
 *
 * - DCNetworkBlock with just DC lines, i.e., where the susceptance for all
 *   lines is a non-zero value.
 *
 * - DCNetworkBlock of a hybrid DC-HVDC grid (both DC and HVDC lines).
 *   This is a combination of first and second cases, where some lines but
 *   not all of them have zero susceptance.
 *
 * These can be implemented with at least three different formulations:
 *
 * - The PTDF formulation using the Power Transfer Distribution Factor
 *   matrix for the DC lines and standard flow conservation constraints
 *   for the HVDC lines.
 *
 * - The CYCLE formulation ... TODO: DESCRIBE
 *
 * - The KIRCHHOFF formulation, which directly encodes Kirchhoff's laws
 *   using both power flow variables F_l and voltage angle variables
 *   theta_n. For each DC line l (non-zero susceptance):
 *     F_l = B_l * ( theta_{from(l)} - theta_{to(l)} )   (KVL)
 *   For each node n:
 *     sum_{l:out(n)} F_l - sum_{l:in(n)} eta_l F_l = S_n - D_n  (KCL)
 *   with a reference node angle fixed to zero. For HVDC lines (zero
 *   susceptance), no angle relationship is imposed: the flow is only
 *   constrained by capacity limits and node balance.
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
 * - DCNetworkData, a small auxiliary class to bunch the basic electrical data
 *   of the transmission network.
 * @{ */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 enum formulation_type
 {
  PTDF = 0 ,
  CYCLE ,
  KIRCHHOFF
  };

/*--------------------------------------------------------------------------*/
/*------------------- CLASS DCNetworkBlock::DCNetworkData ------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
 /// auxiliary class holding basic data about the (DC) transmission network
 /** The DCNetworkData class is a nested sub-class which only serves to have
  * a quick way to load all the basic data (topology and electrical
  * characteristics) that describe the transmission network. The rationale
  * is that, while often the network does not change during the (short) time
  * horizon of UC, it makes sense to allow for this to happen. This means
  * that individual NetworkBlock objects may in principle have different
  * DCNetworkData, but most often they can share the same. By bunching all
  * the information together, we make it easy for this sharing to happen. */

class DCNetworkData : public NetworkData
{
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

 public:

/** @} ---------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/** @name Constructor and Destructor
 * @{ */

 /// constructor of DCNetworkData, does nothing
 DCNetworkData( void ) : f_number_lines( 0 ) , f_number_HVDC_lines( 0 ) ,
  f_reference_node( 0 ) , DCDF_was_computed( false ) ,
  f_number_branches( 0 ) , cycle_basis_was_computed( false ) {}

 /// destructor of DCNetworkData: it is virtual, and empty
 ~DCNetworkData() override = default;

/** @} --------------------- OTHER INITIALIZATIONS -------------------------*/
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
  * information need not be present since it is not loaded. If
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
  *   StartLine[ l ] == EndLine[ l ] (a self-loop) is not allowed, but
  *   multiple lines between the same pair of nodes are. The variable is
  *   mandatory.
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
  *   unit of flow from StartLine[ l ] to EndLine[ l ]. Note that, if l
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
  *   There is no requirement that the efficiencies of the different branches
  *   of the same hyperarc sum to 1: in fact, this variable is optional, if it
  *   is not specified then Efficiency[ l ] == 1 for all branches / lines.
  *
  * - The variable "LineName", of type netCDF::NcString() and indexed over
  *   the dimension "NumberLines". Its i-th entry, namely LineName[ i ],
  *   contains the name of the i-th transmission line. This variable is
  *   optional. */

 virtual void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG
 /// extends NetworkData::expected_dims()

 std::vector< std::string > expected_dims( void ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// extends NetworkData::expected_vars()

 std::vector< std::string > expected_vars( void ) const override;

#endif

/** @} ------- METHODS FOR READING THE DATA OF THE DCNetworkData -----------*/
/** @name Reading the data of the DCNetworkData
 * @{ */

 /// returns the number of lines of the network
 /** Method for returning the total number of lines of the network. When
  * get_number_nodes() == 1 (the network is a bus), get_number_lines() == 0
  * (no self-loops are allowed, hence there is no line to be made with a
  * single node). */

 Index get_number_lines( void ) const { return( f_number_lines ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the number of HVDC lines of the network
 /** Method for returning the number of HVDC lines of the network, i.e.,
  * those with 0 susceptance. Clearly, get_number_HVDC_lines() <=
  * get_number_lines(): when the two are equal this is a HVDC grid, when
  * get_number_HVDC_lines() == 0 this is a DC grid, in between this is a
  * DC-HVDC grid. */

 Index get_number_HVDC_lines( void ) const { return( f_number_HVDC_lines ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns true if this is a pure HVDC grid

 bool is_HVDC( void ) { return( v_line_susceptance.empty() ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns true if this is a pure DC grid

 bool is_DC( void ) { return( f_number_HVDC_lines == 0 ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns true if this is a mixed DC - HVDC grid

 bool is_DC_HVDC( void ) {
  return( ( ! v_line_susceptance.empty() ) && ( f_number_HVDC_lines > 0 ) );
  }
 
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
  * therefore get_end_lines() has to be used. */

 bool is_hypergraph( void ) const {
  return( f_number_branches > f_number_lines );
  }

/*--------------------------------------------------------------------------*/
 /// returns true if \p line is an hyperarc (more than one head bus)

 bool is_hyperarc( Index line ) const {
  if( is_hypergraph() )
   return( v_end_lines[ line ].size() > 1 );
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
  *    get_start_line()[ l ] gives starting (tail) bus of line l. */

 const std::vector< Index > & get_start_line( void ) const {
  return( v_start_line );
  }

/*--------------------------------------------------------------------------*/
 /// returns the start bus of line \p line

 Index get_start_line( Index line ) const {
  return( v_start_line[ line ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of end buses for all lines
 /** Method for returning the vector of ending point of each line. This
  * vector may have empty size (bus network) or the size of number of lines,
  * then there are three possible cases:
  *
  * - if get_number_nodes() == 1, this vector has empty size which means there
  *   is no line at network (bus network), and this vector does not need to
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
  * sense and without repeated elements. */

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
  *   the network. */

 const std::vector< double > & get_min_power_flow( void ) const {
  return( v_min_power_flow );
  }

/*--------------------------------------------------------------------------*/
 /// returns minimum power flow of the given \p line
 /** This method returns the minimum power flow of the given \p line.
  *
  * @return the minimum power flow of the given \p line. */

 double get_min_power_flow( Index line ) const {
  assert( line < f_number_lines );
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
  *   the network. */

 const std::vector< double > & get_max_power_flow( void ) const {
  return( v_max_power_flow );
  }

/*--------------------------------------------------------------------------*/
 /// returns maximum power flow of the given \p line
 /** This method returns the maximum power flow of the given \p line.
  *
  * @return the maximum power flow of the given \p line. */

 double get_max_power_flow( Index line ) const {
  assert( line < f_number_lines );
  if( v_max_power_flow.empty() )
   return( 0 );
  return( v_max_power_flow[ line ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the DC lines
 /** This function returns the DC lines in the transmission network, i.e.,
  * those with nonzero susceptance.
  * @return the AC lines in the network. */

 const Subset & get_DC_lines( void ) {
  if( v_DC_lines.empty() && ( f_number_lines > f_number_HVDC_lines ) ) {
    v_DC_lines.reserve( f_number_lines - f_number_HVDC_lines );
    for( Index line_id = 0 ; line_id < f_number_lines ; ++line_id )
     if( v_line_susceptance[ line_id ] != 0 )
      v_DC_lines.push_back( line_id );
   }

  return( v_DC_lines );
  }

/*--------------------------------------------------------------------------*/
 /// returns the HVDC lines
 /** This function returns the HVDC lines in the transmission network, i.e.,
  * those with zero susceptance.
  * @return the HVDC lines in the network. */

 const Subset & get_HVDC_lines( void ) {
  if( v_HVDC_lines.empty() && ( f_number_HVDC_lines > 0 ) ) {
   if( ! v_line_susceptance.empty() ) {
    v_HVDC_lines.reserve( f_number_HVDC_lines );
    for( Index line_id = 0 ; line_id < f_number_lines ; ++line_id )
     if( v_line_susceptance[ line_id ] == 0 )
      v_HVDC_lines.push_back( line_id );
    }
   else {
    v_HVDC_lines.resize( f_number_lines );
    std::iota( v_HVDC_lines.begin() , v_HVDC_lines.end() , 0 );
    }
   }
  return( v_HVDC_lines );
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
  *   network. */

 const std::vector< double > & get_line_susceptance( void ) const {
  return( v_line_susceptance );
  }

/*--------------------------------------------------------------------------*/

 int get_reducedIdx( int idx ) const;

/*--------------------------------------------------------------------------*/
 // the inverse of get_reducedIdx

 int get_originalIdx( int idx ) const;

/*--------------------------------------------------------------------------*/

 void compute_DCDF( c_Subset & HVDC_lines , const SpMat & PTDF_matrix );

/*--------------------------------------------------------------------------*/

 const SpMat & get_DCDF( void ) const { return( DCDF ); }

/*--------------------------------------------------------------------------*/

 bool was_DCDF_computed( void ) const { return( DCDF_was_computed ); }

/*--------------------------------------------------------------------------*/

 SpMat get_PTDF( c_Subset & DC_lines , double tikhonov_coeff = 1e-4 );

/*--------------------------------------------------------------------------*/

 SpMat get_PTDF( void ) {
  Subset all_lines( f_number_lines );
  std::iota( all_lines.begin() , all_lines.end() , 0 );
  return( get_PTDF( all_lines ) );
  }

/*--------------------------------------------------------------------------*/

 std::pair< SpMat , SpMat > get_stored_B2( void ) {
  return( std::make_pair( stored_B2 , stored_B2_inv ) );
  }

/*--------------------------------------------------------------------------*/

 void set_stored_B2( const SpMat & B2 , const SpMat & B2_inv ) {
  stored_B2 = B2;
  stored_B2_inv = B2_inv;
  }

/*--------------------------------------------------------------------------*/

 bool was_cycle_basis_computed( void ) const {
  return( cycle_basis_was_computed );
  }

/*--------------------------------------------------------------------------*/

 void set_cycle_basis_computed( void ) { cycle_basis_was_computed = true; }

/*--------------------------------------------------------------------------*/
 /// compute the decomposition of the graph into cycles and spanning tree
 /** Methods for computing a spanning tree of the network and a cycle basis
  * The functions get_cycle_basis() and get_spanning_tree() return results
  * in terms of node ids, not line ids. To access the line ids in the
  * spanning tree (resp. in the cycles), use get_lines_in_spanning_tree()
  * (resp. get_lines_in_cycles()). Note that \p root is int and not Index,
  * as it can be negative. */

 void compute_cycle_basis( int root = -1 , bool only_DC_lines = true );

/*--------------------------------------------------------------------------*/

 const std::vector< std::vector< Index > > & get_cycle_basis( void ) {
  if( ! cycle_basis_was_computed )
   this->compute_cycle_basis();
  return( v_cycle_basis );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const std::map< Index , std::set< Index > > & get_spanning_tree( void ) {
  if( ! cycle_basis_was_computed )
   this->compute_cycle_basis();
  return( m_spanning_tree );
  }

 const std::vector< Index > & get_spanning_tree_root( void ) {
  return( m_span_root );
  }

/*--------------------------------------------------------------------------*/
/* Return a map where the keys are the line ids involved in the spanning
 * tree and the value is 1 if the directed line is in the tree and -1 if 
 * the reverse directed line is in the tree.*/

 std::map< Index , int > get_lines_in_spanning_tree( void ) {
  if( ! cycle_basis_was_computed )
   this->compute_cycle_basis();

  const auto number_lines = get_number_lines();
  if( number_lines <= 0 )
   throw( std::logic_error( "DCNetworkData::get_lines_in_spanning_tree: "
			    "number of lines of DCNetworkBlock is not set" )
	  );

  const auto & start_line = get_start_line();
  const auto & end_line = get_end_line();

  std::map< Index , int > lines_in_spanning_tree;
  for( Index line_id = 0 ; line_id < number_lines ; ++line_id ) {
   Index i = start_line[ line_id ];
   Index j = end_line[ line_id ];
   if( this->m_spanning_tree[ i ].contains( j ) )  // line in spanning tree
    lines_in_spanning_tree[ line_id ] = 1;
   else
    if( this->m_spanning_tree[ j ].contains( i ) )  // reverse line in spanning tree
      lines_in_spanning_tree[ line_id ] = -1;
   }
  return( lines_in_spanning_tree );
  }

/*--------------------------------------------------------------------------*/
/* Return a vector of map where the keys are the line ids involved in the
 * cycle and the value is 1 if the directed line is in the cycle and -1 if
 * the reverse directed line is in the cycle. */

 std::vector< std::map< Index , int > > get_lines_in_cycles( void ) {
  if( ! cycle_basis_was_computed )
   this->compute_cycle_basis();

  const auto number_nodes = get_number_nodes();
  const auto number_lines = get_number_lines();
  if( number_lines <= 0 )
   throw( std::logic_error( "DCNetworkData::get_lines_in_spanning_tree: "
			    "number of lines of DCNetworkBlock is not set" )
	  );

  const auto & start_line = get_start_line();
  const auto & end_line = get_end_line();

  std::vector< std::map< Index , int > > lines_in_cycles =
   std::vector< std::map< Index , int > >( this->v_cycle_basis.size() );
  int idx_cycle = 0;
  for (auto & cycle : this->v_cycle_basis ) {
    for (Index line_id = 0 ; line_id < number_lines ; ++line_id ) {
      Index i = start_line[ line_id ];
      Index j = end_line[ line_id ];
      auto it_i = std::find( cycle.begin() , cycle.end() , i );
      int pos_i = std::distance( cycle.begin() , it_i );
      if ( it_i != cycle.end() ) {
        // the line or reverse line may be in the cycle
        if ( ( pos_i < cycle.size() - 1 ) && ( cycle[ pos_i + 1 ] == j ) ) {
         lines_in_cycles[ idx_cycle ][ line_id ] = 1; // true line
        }
        if( ( pos_i == cycle.size() - 1 ) && ( cycle[ 0 ] == j ) ) {
         lines_in_cycles[ idx_cycle ][ line_id ] = 1; // true line
        }
        if( ( pos_i > 0 ) && ( cycle[ pos_i - 1 ] == j ) ) {
         lines_in_cycles[ idx_cycle ][ line_id ] = -1; // reverse line
        }
        if( ( pos_i == 0 ) && ( cycle[ cycle.size() - 1 ] == j ) ) {
         lines_in_cycles[ idx_cycle ][ line_id ] = -1; // reverse line
        }
      }
    }
    ++ idx_cycle;
  }

  assert( lines_in_cycles.size() == number_lines - number_nodes + 1 );
  // take advantage of the theory to ensure the size of the cycle basis

  return( lines_in_cycles );
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
  *   network. */

 const std::vector< double > & get_network_cost( void ) const {
  return( v_network_cost );
  }

/*--------------------------------------------------------------------------*/
 /// returns the efficiency of \p line (1 if not specified or not a HVDC line)

 double get_line_efficiency( Index line ) const {
  assert( line < f_number_lines );
  if( is_hypergraph() )
   throw( std::logic_error( "get_line_efficiency() called but hypergraph" ) );
  if( v_efficiency.empty() ||
      ( ( ! v_line_susceptance.empty() ) &&
	( v_line_susceptance[ line ] != 0.0 ) ) )
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

/** @} --------------- METHODS FOR SAVING THE DCNetworkData ----------------*/
/** @name Methods for loading, printing & saving the DCNetworkData
 * @{ */

 /// serialize a DCNetworkData out of a netCDF::NcGroup
 /** Serialize a DCNetworkData out of a netCDF::NcGroup to the specific
  * format of a DCNetworkData. See
  * DCNetworkBlock::deserialize( netCDF::NcGroup ) for details of the format
  * of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/

 protected:

/*--------------------- PROTECTED METHODS OF THE CLASS ---------------------*/

/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/

 Index f_number_lines;           ///< number of lines of the network
 Index f_number_HVDC_lines;      ///< number of HVDC lines of the network

 /// reference node (used in the PTDF matrix)
 Index f_reference_node;

 /// A boolean to avoid forming A^dc multiple times
 bool DCDF_was_computed;

 /// A boolean to avoid recomputing the cycle basis algorithm
 bool cycle_basis_was_computed;

 Index f_number_branches;     ///< the number of branches of all hyperarcs

 Subset v_start_line;         ///< vector of starting lines

 Subset v_end_line;           ///< vector of ending lines

 /// vector of (vector of) sets of ending lines for hyperarcs
 std::vector< Subset > v_end_lines;

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

 Subset v_DC_lines;      ///< the indices of DC lines (susceptance > 0)
 Subset v_HVDC_lines;    ///< the indices of HVDC lines (susceptance == 0)

 /// vector to store the cycle basis
 std::vector< Subset > v_cycle_basis;

 /// vector to store the spanning tree
 std::map< Index , std::set< Index > > m_spanning_tree;
 /// The spanning tree may have a different root than the reference_node, as a result we store it here
 std::vector< Index > m_span_root; // The graph can in fact have disconnected components because we have HVDC lines that connect the various parts or for some other reason
 // then the spanning tree is more like a "spanning forest" which can be uncovered with the knowledge of multiple roots.

 /// to not recompute each time the PTDF
 SpMat stored_B2;
 SpMat stored_B2_inv;

 /** A SparseMatrix resulting from the product of the PTDF and (A^dc)^T,
  * where the latter is the incidence matrix of the pure DC lines */
 SpMat DCDF;

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

/*-------------------- PRIVATE FIELDS OF THE CLASS -------------------------*/

 SMSpp_insert_in_factory_h;

/*---------------------- PRIVATE METHODS OF THE CLASS ----------------------*/

 };  // end( class( DCNetworkData ) )

/** @} ---------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 * @{ */

 /// constructor of DCNetworkBlock
 /** Constructor of DCNetworkBlock, taking possibly a pointer of its
  * father Block. */

 explicit DCNetworkBlock( Block * f_block = nullptr )
  : NetworkBlock( f_block ) , f_NetworkData( nullptr ) , ftype( PTDF ) ,
    v_design( nullptr ) , f_C_v_scal( 1 ) , f_tikhonov_coeff( 1e-4 ) ,
    f_ptdf_round( 1e-16 ) {}

/*--------------------------------------------------------------------------*/
 /// destructor of DCNetworkBlock

 virtual ~DCNetworkBlock() override;

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 * @{ */

 /// deserialize a DCNetworkBlock out of a netCDF::NcGroup
 /** Deserialize a DCNetworkBlock out of a netCDF::NcGroup, which should
  * contain all the data necessary to describe a NetworkBlock (see
  * NetworkBlock::deserialize()) and possibly the following variables:
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
  *   scalar, let's say k, then it is assumed that Kappa[ l ] = k for each line
  *   l in {0, ..., get_number_lines() - 1}. For each line l in {0, ...,
  *   get_number_lines() - 1}, Kappa[ l ] is the constant that multiplies the
  *   minimum and maximum flow in the flow limit constraints. This variable is
  *   optional. If it is not provided, it is assumed that Kappa[ l ] == 1 for
  *   each line l in {0, ..., get_number_lines() - 1}. */

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG
 // extends NetworkBlock::expected_dims()
 /* not necessary since DCNetworkBlock does not have any new dims save those
  * of the DCNetworkData that are automatically taken into account.

 std::vector< std::string > expected_dims( void ) const override;
 */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// extends NetworkBlock::expected_vars()

 std::vector< std::string > expected_vars( void ) const override;

#endif

/*--------------------------------------------------------------------------*/
 /// loads the DCNetworkBlock instance from a stream
 /** Like load( std::istream & ), if there is any Solver attached to this
  * DCNetworkBlock then a NBModification (the "nuclear option") is issued. */

 void load( std::istream & input , char frmt = 0 ) override {
  throw( std::logic_error( "DCNetworkBlock::load() not implemented yet" ) );
  }

/*--------------------------------------------------------------------------*/
 /// generate the abstract variables of the DCNetworkBlock
 /** The size of node injection variable is the number of intervals spanned
  * this DCNetworkBlock, i.e., 1, by the number of nodes, which can be read
  * via NetworkData::get_number_nodes().
  *
  * Depending on the susceptance for each line of the network, the
  * DCNetworkBlock class may have a power flow variable or not. In other
  * words, if the susceptance is equal to zero (or not defined), the
  * corresponding line is a HVDC line, and it must have the power flow
  * variable. It means, each HVDC line corresponds to a power flow variable,
  * then for the Net Transfer Capacity (NTC) model all lines must have a
  * power flow variable. If the susceptance value is a non-zero value, the
  * corresponding line is called DC and there is no need to define the
  * power flow variable for that line. Therefore, in the case of pure DC lines
  * there is no need to define power flow variables. Consequently, for the
  * mixed case DC-HVDC, the power flow variable must be defined only for HVDC
  * lines. Similarly, depending on the NetworkCost for each line of the
  * network, the DCNetworkBlock class may have an auxiliary variable or not.
  * In other words, if the NetworkCost is equal to zero (or not defined), the
  * auxiliary variable and corresponding constraints will not be defined.
  *
  * TODO: IF THE ABOVE DESCRIPTION ONLY APPLIES TO SOME OF THE FORMULATIONS,
  *       MOVE / REFACTOR AS APPROPRIATE 
  *
  * DCNetworkBlock supports three possible different formulations:
  *
  * - The PTDF formulation, which uses the Power Transfer Distribution
  *   Factor matrix to express DC line flows as a linear combination of
  *   nodal injections. Only power flow variables \f$ F_l \f$ are needed
  *   (no angle variables). For HVDC lines, standard flow conservation
  *   constraints are used.
  *
  * - The CYCLE formulation, which uses a spanning tree and a fundamental
  *   cycle basis to express power flows in terms of tree-transfer
  *   coefficients and cycle flow variables \f$ h_c \f$. Requires both
  *   \f$ F_l \f$ and \f$ h_c \f$ variables (no angle variables).
  *
  * - The KIRCHHOFF formulation, which introduces voltage angle variables
  *   \f$ \theta_n \f$ for each node and encodes:
  *   \f[
  *     F_l \;=\; \mathfrak{S}_l \bigl(\theta_{\mathrm{from}(l)}
  *                                   - \theta_{\mathrm{to}(l)}\bigr)
  *     \qquad \forall l \in \mathcal{L}^{DC}
  *     \qquad (10)
  *   \f]
  *   \f[
  *     - S_n
  *     \;+\;
  *     \sum_{l=(n,\cdot)} F_l
  *     \;-\;
  *     \sum_{l=(\cdot,n)} \eta_l\, F_l
  *     \;=\;
  *     - D^{ac}_n
  *     \qquad \forall n \in \mathcal{N}
  *     \qquad (11)
  *   \f]
  *   \f[
  *     \theta_{\mathrm{ref}} = 0
  *     \qquad (12)
  *   \f]
  *   Capacity limits (1) or (1a)\--(1b) apply as in the other formulations
  *
  * The different possible formulations are represented by a the int value
  * "wf" that is obtained as follows:
  *
  * - if either \p stvv is not nullptr and it is a SimpleConfiguration< int >,
  *   or f_BlockConfig is not nullptr,
  *   f_BlockConfig->f_static_variables_Configuration is not nullptr,
  *   and it is a SimpleConfiguration< int >, then wf is the f_value of the
  *   SimpleConfiguration< int >
  *
  * - otherwise, wf is 0
  *
  * The chosen formulation is CYCLE if wf == 1,KIRCHHOFF if wf == 2,
  * is PTDF in all other cases (default). */

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

/*--------------------------------------------------------------------------*/

 void generate_PTDF_variables( void );

/*--------------------------------------------------------------------------*/

 void generate_CYCLE_variables( void );

/*--------------------------------------------------------------------------*/

 void generate_KIRCHHOFF_variables( void );

/*--------------------------------------------------------------------------*/
 /// generate abstract constraints of DCNetworkBlock
 /** This method generates the linear constraints of the DC network according
  * to the internal formulation type #ftype, which can be **PTDF**, **CYCLE**,
  * or **KIRCHHOFF**. The topology of the transmission network is defined by
  * a set of
  * nodes \f$ \mathcal{N} \f$ and a set of lines \f$ \mathcal{L} \f$. For each
  * line \f$ l \in \mathcal{L} \f$, let \f$ P^{mn}_l \f$ and \f$ P^{mx}_l \f$
  * denote the minimum and maximum admissible power flows, and
  * \f$ \kappa_l \f$ a line-specific scaling factor. For each node
  * \f$ n \in \mathcal{N} \f$, \f$ D^{ac}_n \f$ is the active power demand.
  *
  * The Block defines:
  *
  * - Node injection variables \f$ S_n \f$ for \f$ n \in \mathcal{N} \f$;
  * - Line flow variables \f$ F_l \f$ for \f$ l \in \mathcal{L} \f$;
  * - (optional) Auxiliary variables \f$ V_l \f$ (if a per-line NetworkCost is
  *   defined, used for linearizing \f$ |F_l| \f$);
  *
  * - optionally, design variables \f$ x_l \f$ for a subset of lines can be
  *   externally set (see set_design_variables()).
  *
  * \b Capacity \b limits.
  * Each line \f$ l \f$ is constrained either by a static box or by a
  * design-modulated form, depending on the existence of a design variable:
  *
  * - Without design variable:
  *   \f[
  *     \kappa_l P^{mn}_l \;\le\; F_l \;\le\; \kappa_l P^{mx}_l
  *     \qquad (1)
  *   \f]
  *
  * - With design variable \f$ x_l \f$:
  *   \f[
  *     F_l - \kappa_l P^{mn}_l\, x_l \;\ge\; 0
  *     \qquad (1a)
  *   \f]
  *   \f[
  *     F_l - \kappa_l P^{mx}_l\, x_l \;\le\; 0
  *     \qquad (1b)
  *   \f]
  *
  * \b NetworkCost \b term.
  * When NetworkCost is defined for one or more lines, the absolute value of
  * the power flow \f$ |F_l| \f$ is linearized using an auxiliary variable
  * \f$ V_l \f$ as:
  * \f[
  *   0 \;\le\; V_l - F_l \qquad (2)
  *   \qquad
  *   0 \;\le\; V_l + F_l \qquad (3)
  * \f]
  * and the objective contributes
  * \f$ \sum_{l \in \mathcal{L}} \mathrm{NetworkCost}_l\, V_l \f$.
  *
  * \b HVDC-only \b (NTC) \b formulation.
  * For networks composed exclusively of HVDC lines, all susceptances are zero,
  * and the flows are fully controllable. The nodal power balance reads:
  * \f[
  *   \sum_{l=(n,\cdot)} F_l
  *   \;-\;
  *   \sum_{l=(\cdot,n)} \eta_l\, F_l
  *   \;=\;
  *   S_n - D^{ac}_n
  *   \qquad \forall n \in \mathcal{N} \qquad (4)
  * \f]
  * where \f$ \eta_l \f$ denotes the efficiency of line \f$ l \f$ (possibly
  * different for each branch in hypergraph topologies). Capacity limits
  * follow (1) or (1a)–(1b).
  *
  * \b Hybrid \b DC/HVDC \b (PTDF) \b formulation.
  * DC flows are expressed as a linear combination of nodal injections via the
  * PTDF matrix \f$ B \f$, with coupling to DC flows through a DCDF matrix:
  * \f[
  *   F_l
  *   \;=\;
  *   \sum_{n \in \mathcal{N}} B_{l n}\, \bigl(S_n - D^{ac}_n\bigr)
  *   \;+\;
  *   \sum_{k \in \mathcal{L}^{dc}} \mathrm{DCDF}_{l k}\, F_k
  *   \quad \text{(up to a small numerical slack)}
  *   \qquad (5)
  * \f]
  * Nodal balances are imposed only for nodes impacted by HVDC lines,
  * allowing a small tolerance \f$ \varepsilon \f$:
  * \f[
  *   -S_n
  *   \;+\;
  *   \sum_{l=(n,\cdot)} F_l
  *   \;-\;
  *   \sum_{l=(\cdot,n)} \eta_l\, F_l
  *   \;\in\;
  *   [-D^{ac}_n - \varepsilon,\; -D^{ac}_n + \varepsilon]
  *   \qquad (6)
  * \f]
  * HVDC capacities are treated as in (1) or (1a)–(1b).
  *
  * \b Cycle-flow \b (CYCLE) \b formulation.
  * A cycle basis and a spanning tree are computed. Let \f$ h_c \f$ be cycle
  * flows and \f$ C_{l c} \f$ the cycle incidence coefficients. For each line
  * \f$ l \f$:
  * \f[
  *   F_l
  *   \;=\;
  *   \sum_i T_{l i}\, p_i
  *   \;+\;
  *   \sum_c C_{l c}\, h_c
  *   \qquad (7)
  * \f]
  * where \f$ T_{l i} \f$ is the tree-transfer matrix derived from the
  * spanning tree rooted at the reference node. Kirchhoff’s cycle laws are:
  * \f[
  *   \sum_l \frac{C_{l c}}{\mathfrak{S}_l}\, F_l
  *   \;=\; 0
  *   \qquad \forall c
  *   \qquad (8)
  * \f]
  * and the global power balance is enforced as:
  * \f[
  *   \sum_i p_i
  *   \;=\;
  *   \sum_{n \in \mathcal{N}} (S_n - D^{ac}_n)
  *   \;=\;
  *   0
  *   \qquad (9)
  * \f]
  * Capacity limits (1) or (1a)–(1b) apply to each line depending on whether a
  * design variable \f$ x_l \f$ exists.
  *
  * \b Kirchhoff \b (KIRCHHOFF) \b formulation.
  * Voltage angle variables \f$ \theta_n \f$ are introduced for each node.
  * For each DC line \f$ l \f$ (non-zero susceptance \f$ \mathfrak{S}_l \f$),
  * the flow-angle relationship (Kirchhoff's Voltage Law) is imposed:
  * \f[
  *   F_l
  *   \;=\;
  *   \mathfrak{S}_l \bigl(\theta_{\mathrm{from}(l)}
  *                       - \theta_{\mathrm{to}(l)}\bigr)
  *   \qquad \forall l \in \mathcal{L}^{DC}
  *   \qquad (10)
  * \f]
  * For each node \f$ n \f$, a power balance constraint (Kirchhoff's Current
  * Law) is enforced over \e all lines (both DC and HVDC):
  * \f[
  *   -S_n
  *   \;+\;
  *   \sum_{l=(n,\cdot)} F_l
  *   \;-\;
  *   \sum_{l=(\cdot,n)} \eta_l\, F_l
  *   \;=\;
  *   -D^{ac}_n
  *   \qquad \forall n \in \mathcal{N}
  *   \qquad (11)
  * \f]
  * A reference node angle is fixed to zero:
  * \f[
  *   \theta_{\mathrm{ref}} = 0
  *   \qquad (12)
  * \f]
  * HVDC lines (zero susceptance) have no angle relationship and are only
  * constrained by flow limits and node balance. Hypergraph HVDC lines are
  * supported in the node balance. Capacity limits (1) or (1a)--(1b)
  * apply as in the other formulations.
  *
  * Flow balance constraints may have to be scaled for numerical stability
  * reasons.
  *
  *   TODO: PUT THE SCALING FACTOR IN THE RIGTH EQUATIONS OR AT LEAST TELL
  *   WHICH ONES THEY ARE
  *
  * This is why the following scaling constants are defined:
  *
  * - C_v_scal, with default value of 1 (no scaling);
  *
  * - tikhonov_coeff, with default value of 1e-4, which is used to regularise
  *   (obviously, in the Tikhonov sense) the computation of the inverse in
  *  the PTDF matrix.
  *
  * Setting these to non-default values is possible with the Configuration
  * parameter, that is either \p stcc or, if f_BlockConfig is not nullptr,
  * f_BlockConfig->f_static_constraints_Configuration. If the result is not
  * nullptr, then is is a SimpleConfiguration< ... > which can contain up
  * to two numbers, i.e.,
  *
  * - a SimpleConfiguration< double > for setting C_v_scal alone;
  *
  * - a SimpleConfiguration< std::pair< double , double > > for setting
  *   C_v_scal and tikhonov_coeff. */

 void generate_abstract_constraints( Configuration * stcc = nullptr )
  override;

/*--------------------------------------------------------------------------*/

 void generate_PTDF_constraints( Configuration * stcc = nullptr );

/*--------------------------------------------------------------------------*/

 void generate_CYCLE_constraints( Configuration * stcc = nullptr );

/*--------------------------------------------------------------------------*/

 void generate_KIRCHHOFF_constraints( Configuration * stcc = nullptr );

/*--------------------------------------------------------------------------*/
// Generate the nodal balance equations needed to have HVDC lines.
// the additional boolean can be used to simply overload the model with unnecessary constraints
// 
 void generate_HVDC_nodal_constraints( bool full_formulation = false );

/*--------------------------------------------------------------------------*/
 /// generate the NetworkCost auxiliary constraints
 /** Generates the auxiliary constraints for the linearisation of |F_l|
  * (absolute-value relaxation) when the "NetworkCost" vector is provided:
  *   \f[
  *     V_l \ge  F_l, \quad V_l \ge -F_l \qquad \forall\, l \in \mathcal{L}
  *   \f]
  * Does nothing if NetworkCost is empty. This method is intended to be
  * called by generate_KIRCHHOFF_constraints() and overriding classes. */

 void generate_network_cost_constraints( void );

/*--------------------------------------------------------------------------*/
 /// generate the reference-node angle constraint
 /** Fixes the voltage angle of the reference node to zero:
  *   \f[
  *     \theta_{\mathrm{ref}} = 0
  *   \f]
  * Skipped for pure HVDC networks (no angle variables).
  * This method is intended to be called by generate_KIRCHHOFF_constraints()
  * and overriding classes. */

 void generate_reference_angle_constraint( void );

/*--------------------------------------------------------------------------*/
 /// generate the KCL node-balance constraints
 /** Generates Kirchhoff's Current Law at every node:
  *   \f[
  *     -S_n + \sum_{l:\,\mathrm{start}(l)=n} F_l
  *          - \sum_{l:\,\mathrm{end}(l)=n} \eta_l\, F_l = -D_n
  *     \qquad \forall\, n
  *   \f]
  * handling DC lines, HVDC lines and hypergraph topologies.
  * This method is intended to be called by generate_KIRCHHOFF_constraints()
  * and overriding classes. */

 void generate_node_balance_constraints( void );

/*--------------------------------------------------------------------------*/

 void generate_bound_constraints( void );

/*--------------------------------------------------------------------------*/
 /// a bogus function to round nasty coefficients in the DCOPF equations

 static double round_to( double value , double precision = 1.0 ) {
  return( std::round( value / precision ) * precision );
  }

/*--------------------------------------------------------------------------*/
 /// generate the objective of the DCNetworkBlock
 /** Method that generates the objective of the DCNetworkBlock. The objective
  * can include a linear term on the auxiliary variables associated with
  * network costs, if the vector "NetworkCost" is provided:
  *   \f[
  *     \min \ \sum_{l \in \mathcal{L}} NC_l \cdot V_l
  *   \f]
  *
  *   where \f$ NC_l \f$ is the unit network cost of line \f$l\f$ and
  *   \f$ V_l \f$ is the corresponding auxiliary variable
  *   (coefficients are also scaled by the Block scale factor, if any). */

 void generate_objective( Configuration * objc = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
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
  * @return the number of nodes in the network. */

 Index get_number_nodes( void ) const override {
  if( ! f_NetworkData )
   return( 1 );
  return( f_NetworkData->get_number_nodes() );
  }

/*--------------------------------------------------------------------------*/
 /// returns the number of lines of the network
 /** This function returns the number of lines in the transmission network.
  * If get_NetworkData() returns nullptr, this is equivalent to
  * get_NetworkData()->get_number_lines(). Otherwise, it returns zero.
  *
  * @return the number of lines in the network. */

 Index get_number_lines( void ) const {
  if( ! f_NetworkData )
   return( 0 );
  return( f_NetworkData->get_number_lines() );
  }

/*--------------------------------------------------------------------------*/
 /// returns the susceptance of each lines of the network

 double get_line_susceptance( Index l ) const {
  if( ! ( f_NetworkData ) || f_NetworkData->get_line_susceptance().empty() )
   return( 0 );
  return( f_NetworkData->get_line_susceptance()[ l ] );
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
  * this DCNetworkBlock has no NetworkData, this function returns
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
  return( f_NetworkData ? f_NetworkData->get_max_power_flow( line ) : 0 );
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
  *   variable for line l. */

 const std::vector< ColVariable > & get_power_flow( void ) const {
  return( v_power_flow );
  }

/*--------------------------------------------------------------------------*/

 const std::vector< ColVariable > & get_cycle_flow( void ) const {
  return( v_cycle_flow );
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

/*--------------------------------------------------------------------------*/
 /// returns true if there is any design variable associated with some line

 bool has_design( void ) const {
  return( v_design && ( ! v_design->empty() ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns true if all lines have associated design variable

 bool all_design( void ) const {
  return( has_design() && ( v_design->size() == get_number_lines() ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the design variable associated with the given \p line
 /** This function returns a pointer to the design variable corresponding to
  * the physical transmission line indexed by \p line, or nullptr if there
  * is no such design variable.
  *
  * @param line The index of the line, between 0 and get_number_lines() - 1.
  *
  * @return A pointer to the design variable corresponding to \p line,
  *         or nullptr if that line has no design variable.
  *
  * \note The returned pointer refers to a variable owned externally by the
  *       corresponding DesignNetworkBlock; it must **not** be deleted or
  *       modified outside the intended modeling interface. */

 ColVariable * get_design( Index line ) const {
  if( ( ! v_design ) || v_design->empty() )
   return( nullptr );

  if( ! v_dense_design.empty() )
   return( v_dense_design[ line ] );

  if( line >= v_design->size() )
   return( nullptr );

  return( & (*v_design)[ line ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the (const) design variable associated with the given \p line
 /** Const-qualified overload of get_design(). Returns a const pointer to the
  * design variable corresponding to the given physical line index \p line.
  * See get_design() for details. */

 const ColVariable * get_const_design( Index line ) const {
  return( get_design( line ) );
  }

/** @} ---------------------------------------------------------------------*/
/*--------- METHODS FOR READING THE Constraint OF THE DCNetworkBlock -------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the Constraint of the DCNetworkBlock
 * @{ */

 /// returns the vector of power flow limit constraints
 /** This function returns a const reference to the vector of power flow limit
  * constraints. The i-th element of this vector is a BoxConstraint for the
  * i-th line of the network. */

 const std::vector< BoxConstraint > &
  get_power_flow_limit_constraints( void ) const {
  return( v_power_flow_limit_const );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of power flow limit HVDC bounds

 const std::vector< BoxConstraint > &
  get_power_flow_limit_HVDC_bounds( void ) const {
  return( v_power_flow_limit_const );
  }

/*--------------------------------------------------------------------------*/
 /// return the vector of power losses on lines

 virtual std::vector< double > get_line_losses( void ) const {
  return( std::vector( get_number_lines() , 0. ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the dual prices of power flow limits

 void get_dual_prices( std::vector< double > & dp ) const {
  auto nl = get_number_lines();
  if( ! nl ) {
   dp.clear();
   return;
   }

  dp.resize( nl );
  for( Index l = 0 ; l < nl ; ++l ) {
   if( has_design() && get_design( l ) )  {  // design on this line
    dp[ l ] = v_power_flow_limit_design_const[ 1 ][ l ].get_dual() -
              v_power_flow_limit_design_const[ 0 ][ l ].get_dual();
    continue;
    }

   dp[ l ] = v_power_flow_limit_const.empty() ? 0 :
             v_power_flow_limit_const[ l ].get_dual();
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
 /** The method is actually implemented since DCNetworkBlock is a concrete
  * class. Note that DCNetworkBlock always covers one interval only,
  * hence we expect v[] to contain just get_number_nodes() elements. */

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

  for( Index l = 0 ; l < nl ; ++l )
   v_power_flow_limit_const[ l ].set_dual( dp[ l ] );
  }

/*--------------------------------------------------------------------------*/
 /// set the (shared) design variables
 /** Sets the (shared) design variables, that are used to dimension all the
  * lines in the network. These are shared since they are typically decided
  * once and then used throughout all the (short-term) time horizon.
  *
  * Note that DCNetworkBlock retains the pointers to the two vectors \p DV
  * and \p Which, which therefore must not be changed by the caller for all
  * the lifetime of the object; in turn, DCNetworkBlock cannot change them.
  * (since they are const).
  *
  * If \p DV == nullptr or *DV.empty() then no line has design variables.
  * This is the default if this method is never called.
  *
  * If \p WDV == nullptr or WDV->empty(), it is assumed that DV[ i ]
  * refers to line i for all i = 0 , ... , DV.size() - 1. Otherwise,
  * DV[ i ] refers to line WDV[ i ]. If nonempty, \p WDV is supposed to
  * contain numbers in 0, ..., get_number_lines() - 1, be ordered in
  * increasing sense and without repeated elements.
  *
  * Must be called before DCNetworkBlock::generate_abstract_constraints(). */

 void set_design_variables( std::vector< ColVariable > * DV = nullptr ,
			    c_Subset * WDV = nullptr ) {
  if( constraints_generated() )
   throw( std::logic_error( "DCNetworkBlock::set_design_variables: called "
			    "when constraints are already generated" ) );
  v_design = DV;
  if( WDV && ( ! WDV->empty() ) ) {
   #ifndef NDEBUG
    for( Index i = 0 ; i < WDV->size() - 1 ; ++i )
     if( (*WDV)[ i ] >= (*WDV)[ i + 1 ] )
      throw( std::invalid_argument( "DCNetworkBlock::set_design_variables: "
				    "WDV not ordered" ) );
   #endif
   
   if( WDV->back() >= get_number_lines() )
    throw( std::invalid_argument( "DCNetworkBlock::set_design_variables: "
				  "invalid line number in WDV" ) );

   v_dense_design.resize( get_number_lines() , nullptr );
   for( Index i = 0 , j = 0 ; i < v_dense_design.size() ; ++i )
    if( (*WDV)[ j ] == i )
     v_dense_design[ i ] = & (*v_design)[ j++ ];
   }
  else
   v_dense_design.clear();
  }

/** @} ---------------------------------------------------------------------*/
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
  * @param issueAMod It controls how abstract Modification are issued.
  */

 void set_kappa( MF_dbl_it values , Subset && subset , bool ordered = false ,
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
  * @param issueAMod It controls how abstract Modification are issued.
  */

 void set_kappa( MF_dbl_it values , Range rng = Range( 0 , Inf< Index >() ) ,
                 c_ModParam issuePMod = eNoBlck ,
                 c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// change the abstract representation of the power flow limit constraints
 /** This function changes the abstract representation of the power flow limit
  * constraints for indices in \p modified_lines.
  *
  * @param modified_lines A vector of the indices of constraints that
  *                         have to be modified.
  *
  * @param issueAMod It controls how abstract Modification are issued. */

 void change_power_flow_limit_constraints( c_Subset & modified_lines,
					   c_ModParam issueAMod );

/*--------------------------------------------------------------------------*/
 /// change the abstract repr. of the constraints on the auxiliary variables
 /** This function changes the abstract representation of the constraints on
  * the auxiliary variables for indices in \p modified_lines.
  *
  * @param modified_lines A vector of the indices of constraints that
  *                         have to be modified.
  *
  * @param issueAMod It controls how abstract Modification are issued. */

 void change_relax_abs_constraints( c_Subset & modified_lines ,
				    c_ModParam issueAMod );

/*--------------------------------------------------------------------------*/
 /// change the abstract representation of the injection constraints
 /** This function changes the abstract representation of the power flow
  * injection constraints for indices in \p modified_nodes.
  *
  * @param modified_nodes A vector of the indices of nodes that
  *                         have to be modified.
  *
  * @param issueAMod It controls how abstract Modification are issued. */

 void change_DC_power_flow_injection_constraints( c_Subset & modified_nodes ,
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
  *        modified.
  *
  * @param ordered It indicates whether \p subset is ordered.
  *
  * @param issuePMod It controls how physical Modification are issued.
  *
  * @param issueAMod It controls how abstract Modification are issued. */

 void set_active_demand( MF_dbl_it values , Subset && subset ,
                         bool ordered = false ,
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
  * @param issueAMod It controls how abstract Modification are issued.
  */

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

 DCNetworkData * get_new_NetworkData( void ) const override {
  return( new DCNetworkData() );
  }

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

/*---------------------------------- data ----------------------------------*/

 /// the DCNetworkData object
 DCNetworkData * f_NetworkData;

 /// vector to store the demand of each node of the network
 std::vector< double > v_ActiveDemand;

 std::vector< double > v_kappa;   ///< the kappa constant for each line

 formulation_type ftype;          ///< choice of model

 double f_C_v_scal;               ///< scaling factor for flow bounds

 double f_tikhonov_coeff;         ///< regularization for PTDF computation

 double f_ptdf_round ;           ///< a coefficient to round some of the possibly nasty numerical values in the PTDF matrices

/*-------------------------------- variables -------------------------------*/

 /// the power flow variables
 std::vector< ColVariable > v_power_flow;

 /// the power flow variables on cycle basis
 std::vector< ColVariable > v_cycle_flow;

 /// the voltage angle variables (Kirchhoff formulation)
 std::vector< ColVariable > v_voltage_angle;

 /// the auxiliary network cost variable
 std::vector< ColVariable > v_auxiliary_variable;

 /// the design variables for lines
 std::vector< ColVariable > * v_design;

 /// the "densified" version of v_design (if needed)
 std::vector< ColVariable * > v_dense_design;

/*------------------------------- constraints ------------------------------*/

 /// HVDC power flow and node injection constraints
 std::vector< FRowConstraint > v_power_flow_injection_const;

 /// Mixed DC - HVDC node injection constraints
 std::vector< FRowConstraint > v_DC_HVDC_power_flow_const;

 /// HVDC power flow auxiliary variable constraints
 boost::multi_array< FRowConstraint , 2 > v_power_flow_relax_abs;

 /// Definition of power flow
 std::vector< FRowConstraint > v_power_flow_def;

 /// Power flow limit constraints
 std::vector< BoxConstraint > v_power_flow_limit_const;

 /// Power flow limit design constraints
 boost::multi_array< FRowConstraint , 2 > v_power_flow_limit_design_const;

 /// injection equals to demand
 FRowConstraint overall_balanced_const;

 /// definition of the flow
 std::vector< FRowConstraint > v_CYCLE_def_flow_const;

 /// definition of the flow on cycles
 std::vector< FRowConstraint > v_CYCLE_def_cycle_const;

 /// HVDC constraints for cycle formulation
 std::vector< FRowConstraint > v_CYCLE_def_HVDC_const;

 /// flow-angle definition constraints (Kirchhoff formulation)
 /// F_l - B_l * ( theta_from - theta_to ) = 0 for each DC line
 std::vector< FRowConstraint > v_KIRCHHOFF_power_flow_def;

 /// node power balance constraints (Kirchhoff formulation)
 /// for each node n: -S_n + sum_outgoing F_l - sum_incoming eta_l F_l = -D_n
 std::vector< FRowConstraint > v_KIRCHHOFF_node_balance_const;

 /// reference node angle constraint (Kirchhoff formulation)
 BoxConstraint v_reference_angle_const;

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

 static void static_initialization( void )
 {
  register_method< DCNetworkBlock , MF_dbl_it , Subset && , bool >(
   "DCNetworkBlock::set_active_demand" ,
   & DCNetworkBlock::set_active_demand );

  register_method< DCNetworkBlock , MF_dbl_it , Range >(
   "DCNetworkBlock::set_active_demand" ,
   & DCNetworkBlock::set_active_demand );
  }

/*--------------------------------------------------------------------------*/

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
 DCNetworkBlockMod( DCNetworkBlock * const fblock , int type )
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
                        const Block::Range & rng )
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
 * NetworkBlockSolution). */

class DCNetworkBlockSolution : public NetworkBlockSolution
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*------------------------------- FRIENDS ----------------------------------*/

 friend DCNetworkBlock;  ///< make DCNetworkBlock friend

/*---------- CONSTRUCTING AND DESTRUCTING DCNetworkBlockSolution -----------*/

 /// constructor, does nothing
 explicit DCNetworkBlockSolution( void ) : NetworkBlockSolution() ,
  f_number_lines( 0 ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// deserialize a DCNetworkBlockSolution from a netCDF::NcGroup

 void deserialize( const netCDF::NcGroup & group ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// deserialize a DCNetworkBlockSolution from a "global" netCDF::NcGroup

 void deserialize( const netCDF::NcGroup & group , size_t idx ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 ~DCNetworkBlockSolution() override = default;
 ///< destructor: it is virtual, and empty

/*------ METHODS DESCRIBING THE BEHAVIOR OF A DCNetworkBlockSolution ------*/

 void read( const Block * block ) override;

 void write( Block * block ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// serialize a DCNetworkBlockSolution into a netCDF::NcGroup
 /** Serialize a DCNetworkBlockSolution into a netCDF::NcGroup. The format is
  * the one of NetworkBlockSolution, cf. the comments in
  * NetworkBlockSolution::serialize( netCDF::NcGroup & ), except that
  *
  *     "NumberNetworks" IS NOT REALLY NEEDED, BECAUSE "TotalNumberInstants"
  *     AND "EndInstant" ARE NOT REQUIRED SINCE DCNetworkBlock ALWAYS HAS
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
  *   line l. The variable is optional. */

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
  * Note that the variables are constructed when \p idx == 0 according to
  * the fact that the corresponding DCNetworkBlockSolution has or not been
  * Configure-d to hold them, which means that
  *
  *       ALL THE DCNetworkBlockSolution MUST HAVE BEEN Configure-d IN THE
  *       SAME WAY
  *
  * (although, technically, if some of the DCNetworkBlockSolution that
  * appears when \p idx > 0 is Configure-d with less information than that
  * when idx == 0 the code will not break, but there will be uninitialised
  * values in the netCDF). */

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

/*-------------------------- PROTECTED FIELDS ------------------------------*/

 Index f_number_lines;          ///< the number of lines

 std::vector< double > v_flow;  ///< v_flow[ l ] = flow variable on line l

 std::vector< double > v_cost;  /**< v_cost[ l ] = absolute value of the
                                 *   reduced cost of the capacity constraint
                                 *   of line l */

/*---------------------- PRIVATE PART OF THE CLASS -------------------------*/

 private:

/*---------------------------- PRIVATE FIELDS ------------------------------*/

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
