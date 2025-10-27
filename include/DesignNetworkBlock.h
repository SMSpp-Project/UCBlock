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
 * indexed over the dimension "NumberLines". The Block creates the design
 * variables \f$ x_l \f$ and their bound constraints, contributes the
 * investment term to the objective, and exposes / shares the variables with
 * child NetworkBlock objects (e.g., DCNetworkBlock, ACNetworkBlock) so that
 * all children operate on the very same \f$ x_l \f$ without additional
 * coupling constraints.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
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

#include <unordered_map>

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

 /// constructor of DesignNetworkBlock, taking possibly a pointer of its father
 explicit DesignNetworkBlock( Block * f_block = nullptr )
  : NetworkBlock( f_block ) , f_number_lines( 0 ) {}

 /// destructor of DesignNetworkBlock

 virtual ~DesignNetworkBlock() override;

/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 * @{ */

 /// deserialize a DesignNetworkBlock out of a netCDF::NcGroup
 /** Deserialize a DesignNetworkBlock out of a netCDF::NcGroup, which should
  * contain the following:
  *
  * - The dimension "NumberLines" containing the number of lines for which
  *   design information is provided. The dimension is optional; if it is not
  *   present then it is assumed to be 0, in which case no design variable is
  *   created unless the in-memory interface provides the vectors.
  *
  * - The variable "InvestmentCost" (`NcDouble`, either scalar or indexed over
  *   "NumberLines"): per-line investment costs \f$ I_l \f$. If provided as a
  *   scalar, the value is replicated over all lines. Missing entries default
  *   to 0.
  *
  * - **(Subset-based design selection)** The dimension "NumDesignLines"
  *   containing the number of lines that have a design variable. If present,
  *   the (optional) integer variable "DesignLines", indexed over
  *   "NumDesignLines", lists the **indices of lines** in \f$[0,\ldots,
  *   \mathrm{NumberLines}-1]\f$ that are under design.
  *   If "DesignLines" is **absent**, it is assumed that the designed lines
  *   are exactly \f$ \{ 0 , 1 , \ldots , \mathrm{NumDesignLines}-1 \} \f$
  *   (in this order).
  *
  * - The variables "MinCapacityDesign" and "MaxCapacityDesign" (`NcDouble`)
  *   may be provided **indexed over "NumDesignLines"** (one value per
  *   designed line, in the same order as "DesignLines" or the implicit order
  *   above), or as scalars (replicated). As a backward-compatible option,
  *   they can still be provided indexed over "NumberLines". Defaults:
  *   MinCapacityDesign = 0, MaxCapacityDesign = 1. The sign of the (per-line)
  *   MaxCapacityDesign determines the nature of \f$ x_l \f$:
  *     - if \f$ \mathrm{MaxCapacityDesign}[l] < 0 \f$ then
  *       \f$ x_l \in \{0,1\} \f$ (binary);
  *     - otherwise \f$ x_l \f$ is continuous with
  *       \f$ \mathrm{MinCapacityDesign}[l] \le x_l \le
  *          \mathrm{MaxCapacityDesign}[l] \f$.
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

/**@} ----------------------------------------------------------------------*/
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
 /** Bound constraints on \f$ x_l \f$ are produced according to Min/Max.
  * If \f$ \mathrm{MaxCapacityDesign}[l] < 0 \f$, the variable is binary and
  * the bound box is set to \f$[0,1]\f$ possibly tightened by
  * \f$\mathrm{MinCapacityDesign}[l]\f$. */

 void generate_abstract_constraints( Configuration * stcc = nullptr ) override;

/*--------------------------------------------------------------------------*/

 /// generate the objective of the DesignNetworkBlock
 /** The objective contains the investment term
  * \f[
  *   \sum_{ l \in \mathcal{L} } I_l \cdot x_l \; .
  * \f]
  * where \f$ I_l \f$ is the investment cost of line \f$l\f$ and
  * \f$ x_l \f$ is the (per-line) design variable. */

 void generate_objective( Configuration * objc = nullptr ) override;

 Index get_number_nodes( void ) const override { return( 0 ); }

 NetworkData * get_NetworkData( void ) const override { return( nullptr ); }

 const double * get_active_demand( Index = 0 ) const override { return( nullptr ); }

 void set_ActiveDemand( const boost::multi_array< double , 2 > & ) override {}

 void set_active_demand( MF_dbl_it , Subset && , bool ,
                         ModParam , ModParam ) override {}

 void set_active_demand( MF_dbl_it , Range ,
                         ModParam , ModParam ) override {}

 void set_NetworkData( NetworkData * ) override {}

/**@} ----------------------------------------------------------------------*/
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
  *             and the type of violation that must be considered. */

 bool is_feasible( bool useabstract = false ,
                   Configuration * fsbc = nullptr ) override;

/**@} ----------------------------------------------------------------------*/
/*-------------------------- READING THE DATA ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the DesignNetworkBlock
 * @{ */

 /// returns the number of lines for which design is defined
 Index get_number_lines( void ) const { return( f_number_lines ); }

/*--------------------------------------------------------------------------*/

 /// returns the investment cost associated with the given \p line
 double get_investment_cost( Index line ) const {
  if( v_InvestmentCost.empty() ) return( 0 );
  return( v_InvestmentCost[ line ] );
 }

/*--------------------------------------------------------------------------*/

 /// returns the minimum capacity design associated with the given \p line
 double get_min_capacity_design( Index line ) const {
  if( v_MinCapacityDesign.empty() ) return( 0 );
  return( v_MinCapacityDesign[ line ] );
 }

/*--------------------------------------------------------------------------*/

 /// returns the maximum capacity design associated with the given \p line
 double get_max_capacity_design( Index line ) const {
  if( v_MaxCapacityDesign.empty() ) return( 1 );
  return( v_MaxCapacityDesign[ line ] );
 }

/**@} ---------------------------------------------------------------------*/
/*---------- METHODS FOR READING THE Variable OF THE Block -----------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the Variable of the DesignNetworkBlock
 * @{ */

 /// returns the design variable for the given line
 ColVariable & get_design( Index line ) {
  return( v_design[ line2pos.at( line ) ] );
 }

/*--------------------------------------------------------------------------*/

 /// returns the const design variable for the given line
 const ColVariable & get_const_design( Index line ) const {
  return( v_design[ line2pos.at( line ) ] );
 }

/*--------------------------------------------------------------------------*/

 /// returns the vector of design variables
 const std::vector< ColVariable > & get_design_variables( void ) const {
  return( v_design );
 }

/*--------------------------------------------------------------------------*/

 /// returns true iff the given line has a design variable (belongs to the subset)
 bool is_designed_line( Index l ) const {
  return( line2pos.contains( l ) );
 }

/*--------------------------------------------------------------------------*/

 /// returns the (ordered) list of line indices that have a design variable
 const std::vector< Index > & get_design_lines( void ) const {
  return( v_design_lines );
 }

/**@} ----------------------------------------------------------------------*/
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

 /// number of lines for which design is defined
 Index f_number_lines;

 /// the investment cost for each line
 std::vector< double > v_InvestmentCost;

 /// the minimum capacity design allowed for each line
 std::vector< double > v_MinCapacityDesign;

 /// the maximum capacity design allowed for each line
 std::vector< double > v_MaxCapacityDesign;

/*-------------------------------- variables -------------------------------*/

 /// the design variable for each line (only for the lines listed in v_design_lines)
 std::vector< ColVariable > v_design;

 /// the list of line indices that have an associated design variable
 std::vector< Index > v_design_lines;

 /// map: line index -> position within v_design / v_design_lines
 std::unordered_map< Index , Index > line2pos;

/*------------------------------- constraints ------------------------------*/

 /// the design bound constraint for each line
 std::vector< BoxConstraint > v_design_bound_const;

/*-------------------------------- objective -------------------------------*/

 /// the objective function (investment term)
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

 /// static initialization for method registration (if needed later)
 static void static_initialization( void ) {}

/*--------------------------------------------------------------------------*/

};  // end( class( DesignNetworkBlock ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* __DesignNetworkBlock */

/*--------------------------------------------------------------------------*/
/*-------------------- End File DesignNetworkBlock.h -----------------------*/
/*--------------------------------------------------------------------------*/
