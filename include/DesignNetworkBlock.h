/*--------------------------------------------------------------------------*/
/*---------------------- File DesignNetworkBlock.h -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 *
 * Header file for class DesignNetworkBlock, which derives from Block and
 * defines the standard representation of **per-line design data** and
 * **per-line design variables** for transmission network models in the Unit
 * Commitment problem. The rationale is to centralize once-and-for-all the
 * vectors:
 *
 * - InvestmentCost
 * - MinCapacityDesign
 * - MaxCapacityDesign
 *
 * indexed over the dimension "NumberDesignLines". The Block creates the design
 * variables \f$ x_l \f$ and their bound constraints, contributes the
 * investment term to the objective, and exposes / shares the variables with
 * child NetworkBlock objects (e.g., DCNetworkBlock, ACNetworkBlock) so that
 * all children operate on the very same \f$ x_l \f$ without additional
 * coupling constraints.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __DesignNetworkBlock
 #define __DesignNetworkBlock
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Block.h"

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
/*---------------------- CLASS DesignNetworkBlock --------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// a "design" Block holding investment data and variables shared by children
/** The DesignNetworkBlock class derives from Block, and defines the
 * standard representation of **per-line design** for transmission networks.
 * It reads the design data (InvestmentCost, MinCapacityDesign,
 * MaxCapacityDesign) and creates one design variable \f$ x_l \f$ per line
 * \f$ l \in \mathcal{L} \f$. Bounds of \f$ x_l \f$ are set according to
 * Min/Max; if \f$ \mathrm{MaxCapacityDesign}[l] < 0 \f$ the variable is
 * binary, otherwise it is continuous with
 * \f$ \mathrm{MinCapacityDesign}[l] \le x_l \le \mathrm{MaxCapacityDesign}[l]
 * \f$. The Block also contributes the term
 * \f$ \sum_{l \in \mathcal{L}} I_l x_l \f$ to the objective.
 *
 * Child NetworkBlock objects (e.g., DCNetworkBlock instances) can be made to
 * **share** these variables so that the design appears exactly once in the
 * model while multiple network formulations can reference it. */

class DesignNetworkBlock : public NetworkBlock
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

 /// constructor of DesignNetworkBlock, taking a pointer of its father

 explicit DesignNetworkBlock( Block * f_block = nullptr )
 : NetworkBlock( f_block ), f_NetworkData( nullptr ),
   f_number_subnetworks( 0 ) {}

 /// destructor of DesignNetworkBlock

 virtual ~DesignNetworkBlock() override;

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 * @{ */

 /// deserialize a DesignNetworkBlock out of a netCDF::NcGroup
 /** Deserialize a DesignNetworkBlock out of a netCDF::NcGroup, which should
  * contain the following:
  *
  * - The dimension "NumberDesignLines" containing the number of lines that
  *   have an associated design variable. The dimension may be zero, in which
  *   case no design variable is created unless the in-memory interface
  *   provides the vectors.
  *
  * - The (optional) integer variable "DesignLines" (`NcInt`), indexed over
  *   "NumberDesignLines", listing the **indices of lines** (with respect to
  *   the underlying NetworkData) that are under design. If "DesignLines" is
  *   **absent**, it is assumed that the designed lines are exactly
  *   \f$ \{ 0 , 1 , \ldots , \mathrm{NumberDesignLines}-1 \} \f$ (in this
  *   order).
  *
  * - The variable "InvestmentCost" (`NcDouble`, either scalar or indexed over
  *   "NumberDesignLines"): per-line investment costs \f$ I_l \f$. If provided
  *   as a scalar, the value is replicated over all designed lines. Missing
  *   entries default to 0.
  *
  * - The variables "MinCapacityDesign" and "MaxCapacityDesign" (`NcDouble`)
  *   provided **indexed over "NumberDesignLines"** (one value per designed
  *   line, in the same order as "DesignLines" or the implicit order above),
  *   or as scalars (replicated). Defaults: MinCapacityDesign = 0,
  *   MaxCapacityDesign = 1. The sign of the (per-line) MaxCapacityDesign
  *   determines the nature of \f$ x_l \f$:
  *
  *     = if \f$ \mathrm{MaxCapacityDesign}[ l ] < 0 \f$ then
  *       \f$ x_l \in \{ 0 , 1 \} \f$ (binary);
  *
  *     = otherwise \f$ x_l \f$ is continuous with
  *       \f$ \mathrm{MinCapacityDesign}[ l ] \le x_l \le
  *           \mathrm{MaxCapacityDesign}[ l ] \f$.
  *
  * No network-topological information is handled here; only design data. */

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
 /// serialize a DesignNetworkBlock into a netCDF::NcGroup
 /** Serialize a DesignNetworkBlock into a netCDF::NcGroup to the specific
  * format of a design data provider. See
  * NetworkBlock::deserialize( netCDF::NcGroup ) for details of the format of
  * the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/** @} ---------------------------------------------------------------------*/
/*----------------------- METHODS FOR GENERATION ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for generating Variable / Constraint / Objective
 * @{ */

 /// generate the abstract variables of the DesignNetworkBlock
 /** One variable \f$ x_l \f$ is created for each line \f$ l \f$. The type and
  * bounds follow Min/Max as specified in deserialize(). */

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the abstract constraints of the DesignNetworkBlock
 /** Bound constraints on \f$ x_l \f$ are produced according to Min/Max. If
  * \f$ \mathrm{MaxCapacityDesign}[l] < 0 \f$, the variable is binary and
  * the bound box is set to \f$[0,1]\f$ possibly tightened by
  * \f$\mathrm{MinCapacityDesign}[l]\f$. */

 void generate_abstract_constraints( Configuration * stcc = nullptr )
  override;

/*--------------------------------------------------------------------------*/
 /// generate the objective of the DesignNetworkBlock
 /** The objective contains the investment term
  * \f[
  *   \sum_{ l \in \mathcal{L} } I_l \cdot x_l \; .
  * \f]
  * where \f$ I_l \f$ is the investment cost of line \f$l\f$ and
  * \f$ x_l \f$ is the (per-line) design variable. */

 void generate_objective( Configuration * objc = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
/*-------- METHODS FOR READING THE DATA OF THE DesignNetworkBlock ----------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the NetworkBlock
 * @{ */

 ColVariable * get_node_injection( Index interval = 0 ) override {
  if( v_Block.empty() )
   return( nullptr );

  if( interval < v_Block.size() )
   return( static_cast< NetworkBlock * >(
    v_Block[ interval ] )->get_node_injection() );

  return( nullptr );
  }

/*--------------------------------------------------------------------------*/

 const ColVariable * get_const_node_injection( Index interval = 0 )
 const override {
  if( v_Block.empty() )
   return( nullptr );

  if( interval < v_Block.size() )
   return( static_cast< NetworkBlock * >(
    v_Block[ interval ] )->get_const_node_injection() );

  return( nullptr );
  }

/*--------------------------------------------------------------------------*/

 Index get_number_nodes( void ) const override {
  if( v_Block.empty() )
   return( 0 );
  return( static_cast< NetworkBlock * >(
   v_Block.front() )->get_number_nodes() );
  }

/*--------------------------------------------------------------------------*/

 Index get_number_intervals( void ) const override {
  Index ni = 0;
  for( auto * bi : v_Block )
   ni += static_cast< NetworkBlock * >( bi )->get_number_intervals();

  return( ni );
  }

/*--------------------------------------------------------------------------*/

 NetworkData * get_NetworkData( void ) const override {
  return( f_NetworkData );
  }

/*--------------------------------------------------------------------------*/

 const double * get_active_demand( Index interval ) const override {
  if( v_Block.empty() )
   return( nullptr );

  if( interval < v_Block.size() )
   return( static_cast< NetworkBlock * >(
    v_Block[ interval ] )->get_active_demand( 0 ) );

  return( nullptr );
  }

/*--------------------------------------------------------------------------*/
 /// returns the number of lines for which design is defined

 Index get_number_design_lines( void ) const {
  return( v_design_lines.size() );
  }

/*--------------------------------------------------------------------------*/
 /// returns the investment cost associated with the given \p line

 double get_investment_cost( Index line ) const {
  if( v_InvestmentCost.empty() )
   return( 0 );
  if( const auto idx = get_design_index( line ); idx < Inf< Index >() )
   return( v_InvestmentCost[ idx ] );
  else
   return( 0 );
  }

/*--------------------------------------------------------------------------*/
 /// returns the minimum capacity design associated with the given \p line

 double get_min_capacity_design( Index line ) const {
  if( v_MinCapacityDesign.empty() )
   return( 0 );
  if( const auto idx = get_design_index( line ); idx < Inf< Index >() )
   return( v_MinCapacityDesign[ idx ] );
  else
   return( 0 );
  }

/*--------------------------------------------------------------------------*/
 /// returns the maximum capacity design associated with the given \p line

 double get_max_capacity_design( Index line ) const {
  if( v_MaxCapacityDesign.empty() )
   return( 1 );
  if( const auto idx = get_design_index( line ); idx < Inf< Index >() )
   return( v_MaxCapacityDesign[ idx ] );
  else
   return( 1 );
  }

/** @} ---------------------------------------------------------------------*/
/*---------- METHODS FOR READING THE Variable OF THE Block -----------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the Variable of the DesignNetworkBlock
 * @{ */

/*--------------------------------------------------------------------------*/
 /// returns the index of the line in the vector of design variables
 /** If \p line has an associated design variable then returns the index of
  * that in the vectors of defining them, cf. get_design_variables() and
  * get_design_lines(); that is, get_design_lines( get_design_index( line ) )
  * == line. Otherwise, returns Inf< Index >(). */

 Index get_design_index( Index line ) const {
  if( v_design_lines.empty() )
   return( line );

  const auto it = std::lower_bound( v_design_lines.begin() ,
                                    v_design_lines.end() , line );
  if( ( it == v_design_lines.end() ) || ( *it != line ) )
   return( Inf< Index >() );
  return( static_cast< Index >( std::distance( v_design_lines.begin() , it ) ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the design variable for the given design-line index
 /** Returns a pointer the design variable \f$ x_{\mathrm{line}} \f$ of the
  * given design-line index \p line, or nullptr if the line has no design
  * variable. */

 const ColVariable * get_design( Index line ) const {
  if( v_design_lines.empty() )
   return( line >= v_design.size() ? nullptr : &v_design[ line ] );

  const auto it = std::lower_bound( v_design_lines.begin() ,
                                    v_design_lines.end() , line );

  if( ( it == v_design_lines.end() ) || ( *it != line ) )
   return( nullptr );

  return( &v_design[ static_cast< Index >(
   std::distance( v_design_lines.begin() , it ) ) ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the (const) design variable for the given design-line index
 /** Const-qualified counterpart of get_design(), see it for details. */

 const ColVariable * get_const_design( Index line ) const {
  return( get_design( line ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the (const) vector of design variables

 const std::vector< ColVariable > & get_const_design( void ) const {
  return( v_design );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of design variables

 std::vector< ColVariable > & get_design( void ) {
  return( v_design );
  }

/*--------------------------------------------------------------------------*/
 /// returns the (ordered) list of line indices that have a design variable
 /** Returns the vector of indices of lines that have a design variable,
  * ordered in increasing sense. If it is non-empty, then
  * get_design_variables()[ i ] is the design variable of line
  * get_design_lines()[ i ]. If it is empty, then all lines 0, ...,
  * get_design_variables().size - 1 have a design variable, and
  * get_design_variables()[ i ] is the design variable of line i. */

 c_Subset & get_design_lines( void ) const { return( v_design_lines ); }

/** @} ---------------------------------------------------------------------*/
/*------------ METHODS FOR MODIFYING THE DesignNetworkBlock ----------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the NetworkBlock
 * @{ */

 void set_ActiveDemand( const boost::multi_array< double , 2 > & apd )
  override {
  if( v_Block.empty() )
   return;

  const Index total_intervals = static_cast< Index >( apd.shape()[ 0 ] );
  const Index number_nodes = static_cast< Index >( apd.shape()[ 1 ] );

  Index offset = 0;

  for( auto * nb_ptr : v_Block ) {
   auto * nb = static_cast< NetworkBlock * >( nb_ptr );
   const Index ni = nb->get_number_intervals();

#ifndef NDEBUG
   if( offset + ni > total_intervals )
    throw std::logic_error(
      "DesignNetworkBlock::set_ActiveDemand: "
      "inconsistent number of intervals between UCBlock and subnetworks"
    );
#endif

   boost::multi_array< double , 2 > sub( boost::extents[ ni ][ number_nodes ] );

   for( Index i = 0 ; i < ni ; ++i , ++offset ) {
    auto src_row = apd[ boost::indices[ offset ]
                      [ boost::multi_array_types::index_range( 0 , number_nodes ) ] ];
    std::copy( src_row.begin() , src_row.end() , sub[ i ].begin() );
   }

   nb->set_ActiveDemand( sub );
  }

#ifndef NDEBUG
  if( offset != total_intervals )
   throw std::logic_error(
     "DesignNetworkBlock::set_ActiveDemand: "
     "unused intervals in ActiveDemand matrix"
   );
#endif
  }

/*--------------------------------------------------------------------------*/

 void set_active_demand( MF_dbl_it , Subset && , bool ,
                         ModParam , ModParam ) override final;

/*--------------------------------------------------------------------------*/

 void set_active_demand( MF_dbl_it , Range ,
                         ModParam , ModParam ) override final;

/*--------------------------------------------------------------------------*/

 void set_NetworkData( NetworkData * nd ) override {
  f_NetworkData = nd;
  for( auto * nb : v_Block )
   static_cast< NetworkBlock * >( nb )->set_NetworkData( nd );
  }

/*--------------------------------------------------------------------------*/

 void set_constant_term( const double const_term ) override {
  for( auto * nb : v_Block )
   static_cast< NetworkBlock * >( nb )->set_constant_term( const_term );
  }

/*--------------------------------------------------------------------------*/

 void set_min_node_injection( const double min_injection ,
                              Index node ,
                              Index interval = 0 ) override {
  if( v_Block.empty() )
   return;

  if( interval < v_Block.size() )
   static_cast< NetworkBlock * >(
    v_Block[ interval ] )->set_min_node_injection( min_injection , node );
  }

/*--------------------------------------------------------------------------*/

 void set_max_node_injection( const double max_injection ,
                              Index node ,
                              Index interval = 0 ) override {
  if( v_Block.empty() )
   return;

  if( interval < v_Block.size() )
   static_cast< NetworkBlock * >(
    v_Block[ interval ] )->set_min_node_injection( max_injection , node );
  }

/** @} ---------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Solution
 * @{ */

 /// returns a Solution representing the current solution of this NetworkBlock
 /** This method must construct and return a (pointer to a) Solution object
  * representing the current "solution state" of this NetworkBlock.�, i.e., a
  * DesignNetworkBlockSolution.
  *
  * The parameter for deciding which kind of Solution must be returned is a
  * single int value, coded bitwise:
  *
  * - bit 0 (& 1) means "store the node injection"
  * - bit 1 (& 2) means "store the design variables"
  * - bit 2 (& 4) means "store the sub-Networks": note that if this is 0
  *               then bit 0 is ignored since the node injections are saved
  *               in the sub-NetworkBlockSolution
  * - bit 3 (& 8) means "store the sub-Networks in compressed format"
  * - all subsequent bits, if nonzero, are used to configure the
  *   sub-NetworkBlockSolution of the DesignNetworkBlockSolution. This
  *   requires an int value, that is generated as follows: the first bit
  *   is copied over from [bit 0], i.e., if the  node injections have to be
  *   saved then they are saved in the sub-NetworkBlockSolution; then, all
  *   bits 4 - ... from the int are copied as the bits 1 - ... of the values
  *   passed to the sub-objects.
  *
  * This value is to be found as:
  *
  * - if solc is not nullptr and it is a SimpleConfiguration< int >, then it
  *   is solc->f_value;
  *
  * - otherwise, if f_BlockConfig is not nullptr,
  *   f_BlockConfig->f_solution_Configuration is not nullptr, and it is a
  *   SimpleConfiguration< int >, then it is
  *   f_BlockConfig->f_solution_Configuration->f_value;
  *
  * - otherwise, it is 15 + 24 = 31 (save everything in compressed format
  *   and set bit 1 and 2 of the sub-Networks, which corresponds to the
  *   default value in DCNetworkBlock). */

 Solution * get_Solution( Configuration * solc = nullptr ,
                          bool emptys = true ) override;

/*--------------------------------------------------------------------------*/
 /// return the "appropriate" [Design]NetworkBlockSolution

 NetworkBlockSolution * new_Solution( void ) const override;

/** @} ---------------------------------------------------------------------*/
/*-------------- Methods for checking the DesignNetworkBlock ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for checking solution information in the DesignNetworkBlock
 * @{ */

 /// returns true if the current solution is (approximately) feasible
 /** This function returns true if and only if the solution encoded in the
  * current value of the Variable of this DesignNetworkBlock is approximately
  * feasible within the given tolerance. That is, a solution is considered
  * feasible if and only if
  *
  * -# each ColVariable is feasible; and
  *
  * -# the violation of each Constraint of this DesignNetworkBlock is not
  *    greater than the tolerance.
  *
  * Every Constraint of this DesignNetworkBlock is a RowConstraint and its
  * violation is given by either the relative (see RowConstraint::rel_viol())
  * or the absolute violation (see RowConstraint::abs_viol()), depending on
  * the Configuration that is provided.
  *
  * The tolerance and the type of violation can be provided by either \p fsbc
  * or #f_BlockConfig->f_is_feasible_Configuration, and they are determined as
  * follows:
  *
  * - If \p fsbc is not a nullptr, and it is a pointer to a
  *   SimpleConfiguration< double >, then the tolerance is the value present
  *   in that SimpleConfiguration and the relative violation is considered.
  *
  * - If \p fsbc is not nullptr, and it is a
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
  *             and the type of violation that must be considered. */

 bool is_feasible( bool useabstract = false ,
                   Configuration * fsbc = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
/*-------------------------- READING THE DATA ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the DesignNetworkBlock
 * @{ */

/** @} ---------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED METHODS OF THE CLASS ---------------------*/
/*--------------------------------------------------------------------------*/

 /// verify whether the design data in this Block is consistent
 /** This function checks whether the design data in this DesignNetworkBlock
  * is consistent. The data is consistent if, **for each line** \f$ l \f$,
  * all the following conditions are met:
  *
  * - \f$ \mathrm{MinCapacityDesign}[l] \ge 0 \f$;
  * - if \f$ \mathrm{MaxCapacityDesign}[l] > 0 \f$, then
  *   \f$ \mathrm{MinCapacityDesign}[l] \le \mathrm{MaxCapacityDesign}[l] \f$;
  * - if \f$ |\mathrm{MaxCapacityDesign}[l]| = 1 \f$, then
  *   \f$ \mathrm{MinCapacityDesign}[l] \le 1 \f$;
  * - if \f$ \mathrm{MaxCapacityDesign}[l] < 0 \f$ (binary), then
  *   \f$ \mathrm{MinCapacityDesign}[l] \le 1 \f$
  *   (note: \f$ \mathrm{MinCapacityDesign}[l] > 0 \Rightarrow x_l = 1 \f$).
  *
  * If any of the above conditions are not met for some line \f$ l \f$,
  * an exception is thrown. */

 void check_data_consistency( void ) const;

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

/*---------------------------------- data ----------------------------------*/

 NetworkData * f_NetworkData;            ///< the NetworkData object

 Index f_number_subnetworks;             ///< number of subnetworks

 /// the investment cost for each line
 std::vector< double > v_InvestmentCost;

 /// the minimum capacity design allowed for each line
 std::vector< double > v_MinCapacityDesign;

 /// the maximum capacity design allowed for each line
 std::vector< double > v_MaxCapacityDesign;

/*-------------------------------- variables -------------------------------*/

 /// the design variable for each design line
 std::vector< ColVariable > v_design;

 /// the list of line indices that have an associated design variable
 Subset v_design_lines;

/*------------------------------- constraints ------------------------------*/

 /// the design bound constraint for each line
 std::vector< BoxConstraint > v_design_bound_const;

/*-------------------------------- objective -------------------------------*/

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

 /// deserialize the Network Blocks

 void deserialize_network_blocks( const netCDF::NcGroup & group );

/*--------------------------------------------------------------------------*/

 };  // end( class( DesignNetworkBlock ) )

/*--------------------------------------------------------------------------*/
/*------------------ CLASS DesignNetworkBlockSolution ----------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// a [NetworkBlock]Solution of a DesignNetworkBlock
/** The DesignNetworkBlockSolution class derives from NetworkBlockSolution
 * and adds the "standard" information stored in there (the node injection
 * variables) the other information that is typical of the
 * DesignNetworkBlock, i.e.,
 *
 * - the design variables on (a subset of) the link(s)
 *
 * - the [DC]NetworkBlockSolution information corrseponding to the inner
 *   [DC]NetworkBlock, in basically the same format as that of UCBlock,
 *   i.e., in two possible versions:
  *
  *   = the "standard" one, i.e., sub-groups "NetworkBlock_0",
  *     "NetworkBlock_1", ..., "NetworkBlock_T" with "NetworkBlock_t"
  *     containing the NetworkBlockSolution corresponding to the network
  *     constraints at some specific time instant t;
  *
  *   = all the data corresponding to all the time instants that the
  *     DesignNetworkBlock covers "compressed" in variables in the given
  *     group, see NetworkBlock::serialize( group & , int ) for a
  *     description of the format. */

class DesignNetworkBlockSolution : public NetworkBlockSolution
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*------------------------------- FRIENDS ----------------------------------*/

 friend DesignNetworkBlock;  ///< make DesignNetworkBlock friend

/*-------- CONSTRUCTING AND DESTRUCTING DesignNetworkBlockSolution ---------*/

 /// constructor, does nothing
 explicit DesignNetworkBlockSolution( void )
 : f_design_lines( 0 ), f_number_subnetworks( 0 ), f_compressed( false ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// deserialize a DesignNetworkBlockSolution from a netCDF::NcGroup

 void deserialize( const netCDF::NcGroup & group ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// deserialize a DesignNetworkBlockSolution from a "global" netCDF::NcGroup

 void deserialize( const netCDF::NcGroup & group , size_t idx ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// destructor: it is virtual

 ~DesignNetworkBlockSolution() override  {
   for( auto nbs : v_network_Solution )
    delete nbs;
   }

/*----------- READING THE DATA OF THE DesignNetworkBlockSolution ----------*/

/*---- METHODS DESCRIBING THE BEHAVIOR OF A DesignNetworkBlockSolution ----*/

 void read( const Block * block ) override final;

 void write( Block * block ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// serialize a DesignNetworkBlockSolution into a netCDF::NcGroup
 /** Serialize a DesignNetworkBlockSolution into a netCDF::NcGroup. The
  * format is the one of NetworkBlockSolution, cf. the comments in
  * NetworkBlockSolution::serialize( netCDF::NcGroup & ), plus
  *
  * - The dimension "NumberDesignLines" containing the number of lines in
  *   the transmission network that have design variables. It is opitonal,
  *   but if it's not there then "DesignValue" must not be there.
  *
  * - The variable "DesignValue", of type netCDF::NcDouble and indexed over
  *   the dimension "NumberDesignLines"; DesignValue[ l ] is the optimal
  *   value of the design variable l. The variable is optional, but if
  *   "NumberDesignLines" is there then it must also be there.
  *
  * - All the information relative to the sub-NetworkBlock, either in
  *   "standard" format, i.e.,
  *
  *   = The dimension "NumberSubNetwork" containing the number of
  *     sub-NetworkBlock (and, therefore, their NetworkBlockSolution) in the
  *     DesignNetworkBlock (and therefore DesignNetworkBlockSolution). It is
  *     optional, but if it's not there then the sub-groups (see below)
  *     cannot be there.
  *
  *   = Sub-groups "SubNetworkBlock_0", "SubNetworkBlock_1", ...,
  *     "SubNetworkBlock_T" with T = NumberSubNetworks - 1, with
  *     "SubNetworkBlock_i" containing each the NetworkBlockSolution
  *     corresponding to that network constraints for some specific subset
  *     of time instants.
  *
  *   or in "nonstandard" format, cf. the comments to
  *   NetworkBlockSolution::serialize( netCDF::NcGroup & , * size_t ), where
  *   the sub-NetworkBlock are all serialize()-d in \p group. */

 void serialize( netCDF::NcGroup & group ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// serialize a DesignNetworkBlockSolution into a "global" netCDF::NcGroup
 /** fake "nonstandard" version of serialize() that loads a
  * DesignNetworkBlockSolution from a "global" netCDF::NcGroup, i.e., one
  * where supposedly the solution information of multiple DesignNetworkBlock
  * are stored together. However, this is "fake" in the sense that the
  * format is the "nonstandard" one of NetworkBlockSolution, cf. the
  * comments to NetworkBlockSolution::serialize( netCDF::NcGroup & ,
  * size_t ), plus:
  *
  * - A single new sub-group "DesignNetworkBlock_<idx>" containing all the
  *   information in the same format as serialize( netCDF::NcGroup & ).
  *
  * That is, in the DesignNetworkBlockSolution case the "nonstandard" format
  * is group-based basically ad the "standard" one. Note, however, that
  * inside the group the sub-NetworkBlockSolution can be stored in the
  * "truly nonstandard" (compressed) form -- but this is true even for the
  * "standard" case. */

 void serialize( netCDF::NcGroup & group , size_t idx ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 DesignNetworkBlockSolution * scale( double factor ) const override;

 void sum( const Solution * solution , double multiplier ) override;

 DesignNetworkBlockSolution * clone( bool empty = false ) const override;

/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/

 void print( std::ostream &output ) const override {
  output << "DesignNetworkBlockSolution [" << this << "]: " << std::endl;
  }

/*---------------------- PRIVATE PART OF THE CLASS -------------------------*/

 private:

/*---------------------------- PRIVATE FIELDS ------------------------------*/

 Index f_design_lines;          ///< the number of lines with design variables

 Index f_number_subnetworks;    ///< the number of sub-NetworkBlockSolution

 bool f_compressed;
 ///< true if sub-NetworkBlockSolution are stored in compressed format

 std::vector< double > v_design;
 ///< v_design[ l ] = design variable on constructable line l

 std::vector< NetworkBlockSolution * > v_network_Solution;
 ///< the Solution for each sub-NetworkBlock

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class( DesignNetworkBlockSolution ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* __DesignNetworkBlock */

/*--------------------------------------------------------------------------*/
/*-------------------- End File DesignNetworkBlock.h -----------------------*/
/*--------------------------------------------------------------------------*/
