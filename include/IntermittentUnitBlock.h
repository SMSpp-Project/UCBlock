/*--------------------------------------------------------------------------*/
/*----------------------- File IntermittentUnitBlock.h ---------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class IntermittentUnitBlock, which derives from
 * UnitBlock [see UnitBlock.h], in order to define a Unit representing
 * Intermittent Generation in the Unit Commitment Problem.
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
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Ali Ghezelsoflu,
 *                      Rafael Durbano Lobato, Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __IntermittentUnitBlock
 #define __IntermittentUnitBlock
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "ColVariable.h"

#include "FRowConstraint.h"

#include "OneVarConstraint.h"

#include "UnitBlock.h"

#include "FRealObjective.h"

/*--------------------------------------------------------------------------*/
/*------------------------------ NAMESPACE ---------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS IntermittentUnitBlock -----------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// implementation of the Block concept for the Intermittent Generation unit
/** The IntermittentUnitBlock class implements the Block concept
 * [see Block.h] for units representing generation (be it centralized or
 * distributed) by intermittent (= unreliable) sources in the unit
 * commitment problem, such as wind farms, solar parks and run-of-the-river
 * hydroelectricity. Each unit is supposed to be connected to a specific
 * node of the clustered network (which means that the "distributed" case
 * refers to "distributed in a small region", where of course "small"
 * depends on the granularity of the network description). The model relies
 * mainly on historical data of local generation of wind and solar at each
 * node of the grid; these data are used to develop normalized generation
 * profiles associated with wind and solar generators. Intermittent
 * generators are supposed to be able to contribute to primary and
 * secondary reserves. Contribution to the system inertia more specifically
 * concerns run-of-the-river generators. The potential contribution of
 * solar or wind generation to inertia is still the subject of active
 * research. Reserve requirements are specified in order to be
 * symmetrically available to increase or decrease power injected into the
 * grid. Then the technical and physical constraints are mainly divided in
 * three different categories:
 *
 * - the active power bounds;
 *
 * - the maximum power output constraints according to primary and secondary
 *   spinning reserves;
 *
 * - the minimum power output constraints according to primary and secondary
 *   spinning reserves.
 */

class IntermittentUnitBlock : public UnitBlock
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 * @{ */

 /// constructor, takes the father block
 /** Constructor of IntermittentUnitBlock, taking possibly a pointer of its
  * father Block.
  */
 explicit IntermittentUnitBlock( Block * f_block = nullptr )
  : UnitBlock( f_block ), f_InvestmentCost( 0 ), f_MinCapacityDesign( 0 ),
    f_MaxCapacityDesign( 1 ), f_gamma( 0 ), f_kappa( 1 ), f_scale( 1 ),
    f_max_power_epsilon( 0 ) {}

/*--------------------------------------------------------------------------*/
 /// destructor of IntermittentUnitBlock
 virtual ~IntermittentUnitBlock() override;

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 * @{ */

 /// extends Block::deserialize( netCDF::NcGroup )
 /** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
  * the IntermittentUnitBlock. Besides the mandatory "type" attribute of any
  * :Block, the group must contain all the data required by the base
  * UnitBlock, as described in the comments to UnitBlock::deserialize(
  * netCDF::NcGroup ). In particular, we refer to that description for the
  * crucial dimensions "TimeHorizon", "NumberIntervals" and
  * "ChangeIntervals". The netCDF::NcGroup must then also contain:
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
  *
  *   - if \( \mathrm{MaxCapacityDesign} < 0 \) (binary design), then
  *     \( x \in \{ 0 , 1 \} \) and \( \mathrm{MinCapacityDesign} > 0 \)
  *     implies \( x = 1 \);
  *
  *   - otherwise (continuous design), \( x \) is nonnegative continuous with
  *     \( \mathrm{MinCapacityDesign} \le x \le \mathrm{MaxCapacityDesign} \).
  *
  * - The scalar variable "MaxCapacityDesign", of type netCDF::NcDouble and
  *   not indexed over any dimension. This limits the design variable \( x \):
  *
  *   - if \( \mathrm{MaxCapacityDesign} < 0 \) then \( x \in \{ 0 , 1\} \)
  *     (binary);
  *
  *   - if \( \mathrm{MaxCapacityDesign} = 1 \) then \( x \in [ 0 , 1 ] \)
  *     when \( \mathrm{MinCapacityDesign} = 0 \); otherwise
  *     \( x \in [\,\mathrm{MinCapacityDesign},\,1] \);
  *
  *   - if \( \mathrm{MaxCapacityDesign} > 0 \) then \( x \) is nonnegative
  *     continuous with \( \mathrm{MinCapacityDesign} \le x \le
  *     \mathrm{MaxCapacityDesign} \).
  *
  *   If not provided, the default is 1.
  *
  * - The scalar variable "MaxCapacity", of type netCDF::NcDouble and not
  *   indexed over any dimension. This is the maximum installable capacity
  *   chosen by the user (used for consistency checks and/or reporting).
  *
  * - The variable "MinPower", of type netCDF::NcDouble and either of size 1
  *   or indexed over the dimension "NumberIntervals" (if "NumberIntervals"
  *   is not provided, then this variable can also be indexed over
  *   "TimeHorizon"). This is meant to represent the vector
  *   \f$ \mathrm{MinP} [ t ] \f$ that, for each time instant \f$ t \f$,
  *   contains the minimum potential production value of the unit for the
  *   corresponding time step. If "MinPower" has length 1 then
  *   \f$ \mathrm{MinP}[ t ] \f$ contains the same value for all \f$ t \f$.
  *   Otherwise, \f$ \mathrm{MinPower}[ i ] \f$  is the fixed value of
  *   \f$ \mathrm{MinP}[ t ] \f$ for all \f$ t \f$ in the interval
  *   \f$ [ \mathrm{ChangeIntervals}[ i - 1 ] ,
  *         \mathrm{ChangeIntervals}[ i ] ] \f$, with the assumption that
  *   \f$ \mathrm{ChangeIntervals}[ - 1 ] = 0 \f$. If
  *   \f$ \mathrm{NumberIntervals} \le 1 \f$ or
  *   \f$ \mathrm{NumberIntervals} \ge \mathrm{TimeHorizon} \f$, then the
  *   mapping clearly does not require "ChangeIntervals", which in fact is
  *   not loaded. Note that it must be \f$ \mathrm{MinP}[ t ] \ge 0 \f$ for
  *   all \f$ t \f$.
  *
  * - The variable "MaxPower", of type netCDF::NcDouble and either of size 1
  *   or indexed over the dimension "NumberIntervals" (if "NumberIntervals"
  *   is not provided, then this variable can also be indexed over
  *   "TimeHorizon"). This is meant to represent the vector \f$ \mathrm{MaxP}
  *   [ t ] \f$ that, for each time instant \f$ t \f$, contains the maximum
  *   potential production value of the unit for the corresponding time step.
  *   If "MaxPower" has length 1 then \f$ \mathrm{MaxP}[ t ] \f$ contains the
  *   same value for all \f$ t \f$. Otherwise, \f$ \mathrm{MaxPower}[ i ] \f$
  *   is the fixed value of \f$ \mathrm{MaxP}[ t ] \f$ for all \f$ t \f$ in
  *   the interval \f$ [ \mathrm{ChangeIntervals}[ i - 1 ] ,
  *   \mathrm{ChangeIntervals}[ i ] ] \f$, with the assumption that
  *   \f$ \mathrm{ChangeIntervals}[ - 1 ] = 0 \f$. If
  *   \f$ \mathrm{NumberIntervals} \le 1 \f$ or
  *   \f$ \mathrm{NumberIntervals} \ge \mathrm{TimeHorizon} \f$, then the
  *   mapping clearly does not require "ChangeIntervals", which in fact is
  *   not loaded. Note that it must be
  *   \f$ \mathrm{MaxP}[ t ] \ge \mathrm{MinP}[ t ] \ge 0 \f$ for all
  *   \f$ t \f$. Yet, \f$ \mathrm{MaxP}[ t ] = \mathrm{MinP}[ t ] \f$ is
  *   possible: it means that (at time instant \f$ t \f$) the unit cannot be
  *   curtailed and cannot provide any reserve.
  *
  * - The variable "InertiaPower", of type netCDF::NcDouble and either of
  *   size 1 or indexed over the dimension "NumberIntervals" (if
  *   "NumberIntervals" is not provided, then this variable can also be
  *   indexed over "TimeHorizon"). This is meant to represent the vector
  *   \f$ \mathrm{IP}[ t ] \f$ which, for each time instant \f$ t \f$,
  *   contains the contribution that the unit can give to the inertia
  *   constraint which depends on the active power that it is currently
  *   generating (basically, the constant to be multiplied by the active
  *   power variable) at time \f$ t \f$ for this unit. The variable is
  *   optional; if it is not defined, \f$ \mathrm{IP}[ t ] = 0 \f$ for each
  *   time instant \f$ t \f$. If it has size 1 then the entry
  *   \f$ \mathrm{IP}[ 0 ] \f$ is assumed to contain the inertia power value
  *   for this unit and all time instants \f$ t \f$. Otherwise,
  *   \f$ \mathrm{InertiaPower}[ i ] \f$ is the fixed value of
  *   \f$ \mathrm{IP}[ t ] \f$ for all \f$ t \f$ in the interval
  *   \f$ [ \mathrm{ChangeIntervals}[ i - 1 ] ,
  *         \mathrm{ChangeIntervals}[ i ] ] \f$, with the assumption that
  *   \f$ \mathrm{ChangeIntervals}[ - 1 ] = 0 \f$. If
  *   \f$ \mathrm{NumberIntervals} \le 1 \f$ or
  *   \f$ \mathrm{NumberIntervals} \ge \mathrm{TimeHorizon} \f$ then the
  *   mapping clearly does not require "ChangeIntervals", which in fact is
  *   not loaded.
  *
  * - The variable "ActivePowerCost", of type netCDF::NcDouble and either of
  *   size 1 or indexed over the dimension "NumberIntervals" (if
  *   "NumberIntervals" is not provided, then this variable can also be
  *   indexed over "TimeHorizon"). This is meant to represent the vector
  *   \f$ B[ t ] \f$ that, for each time instant \f$ t \f$, contains the
  *   coefficient of the linear active-power cost function of the unit for
  *   the corresponding time step. This variable is optional; if it is not
  *   provided then it is assumed that \f$ B[ t ] = 0 \f$, i.e., the cost of
  *   the unit has no linear dependence on the produced power (say, only the
  *   quadratic one). If "ActivePowerCost" has length 1 then \f$ B[ t ] \f$
  *   contains the same value for all \f$ t \f$. Otherwise,
  *   \f$ \mathrm{ActivePowerCost}[ i ] \f$ is the fixed value of
  *   \f$ B[ t ] \f$ for all \f$ t \f$ in the interval
  *   \f$ [ \mathrm{ChangeIntervals}[ i - 1 ] ,
  *         \mathrm{ChangeIntervals}[ i ] ] \f$, with the assumption that
  *   \f$ \mathrm{ChangeIntervals}[ - 1 ] = 0 \f$. If
  *   \f$ \mathrm{NumberIntervals} \le 1 \f$ or
  *   \f$ \mathrm{NumberIntervals} \ge \mathrm{TimeHorizon} \f$ then the
  *   mapping clearly does not require "ChangeIntervals", which in fact is
  *   not loaded.
  *
  * - The scalar variable "Gamma", of type netCDF::NcDouble and not indexed
  *   over any dimension. This variable is used to take into account an
  *   uncertainty on the maximal potential production. Note that it must be
  *   \f$ 0 \le \Gamma \le 1 \f$; when \f$ \Gamma = 0 \f$, the unit does not
  *   provide any reserve.
  *
  * - The scalar variable "Kappa", of type netCDF::NcDouble and not indexed
  *   over any dimension. This variable multiplies the minimum and maximum
  *   power at each time instant \f$ t \f$. This variable is optional; if it
  *   is not provided it is taken to be \f$ \kappa = 1 \f$. */

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG
 /* extends UnitBlock::expected_dims()
  * not necessary, no new dimensions

 std::vector< std::string > expected_dims( void ) const override;
 */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// extends UnitBlock::expected_vars()

 std::vector< std::string > expected_vars( void ) const override;

#endif

/*--------------------------------------------------------------------------*/
 /// generate the abstract variables of the IntermittentUnitBlock
 /** The IntermittentUnitBlock class has three different variables which are:
  *
  * - the primary spinning reserve variables;
  *
  * - the secondary spinning reserve variables;
  *
  * - the active power variables.
  *
  * All of those variables are optional except the active power variables,
  * in the sense that the model may just not have them; whenever a group of
  * the above variables is created, its size will be the time horizon.
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
  * When \( \mathrm{MaxCapacityDesign} < 0 \) (binary design),
  * \( x \in \{ 0 , 1 \} \); if \( \mathrm{MinCapacityDesign} > 0 \), then
  * \( x \) is effectively forced to 1. */

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the static constraints of the IntermittentUnitBlock
 /** Method that generates the static constraints of the IntermittentUnitBlock.
  * These are:
  *
  * - maximum and minimum power output constraints according to primary and
  *   secondary spinning reserves, given by (1)–(2). Each of them is stored
  *   as a std::vector< FRowConstraint >, with dimension get_time_horizon().
  *   The entry \f$ t \f$, for
  *   \f$ t \in \mathcal{T} = \{ 0 , \ldots ,
  *   \mathrm{get\_time\_horizon}() - 1 \} \f$, refers to time \f$ t \f$.
  *   These constraints ensure the maximum (or minimum) amount of power that
  *   the unit can produce (or reduce) at time \f$ t \f$.
  *
  *   \f[
  *     p^{pr}_t + p^{sc}_t \leq \gamma \big( \kappa P^{mx}_t - p^{ac}_t \big)
  *                                         \quad t \in \mathcal{T} \quad (1)
  *   \f]
  *
  *   \f[
  *     p^{pr}_t + p^{sc}_t \leq p^{ac}_t - \big( \kappa P^{mn}_t \big)
  *                                         \quad t \in \mathcal{T} \quad (2)
  *   \f]
  *
  *   where \f$ P^{mx}_t \f$ and \f$ P^{mn}_t \f$ are the maximum and
  *   minimum power output parameters for each time \f$ t \f$ in
  *   \f$ \mathcal{T} \f$, respectively.
  *
  * - the active power bounds:
  *
  *   \f[
  *     p^{ac}_t \in [ \kappa P^{mn}_t , \kappa P^{mx}_t ]
  *                                         \quad t \in \mathcal{T} \quad (3a)
  *   \f]
  *
  *   which, in the design scenario of the UC problem, become:
  *
  *   \f[
  *     x \, ( \kappa P^{mn}_t ) \leq p^{ac}_t \leq
  *     x \, ( \kappa P^{mx}_t ) \quad t \in \mathcal{T} \quad (3b)
  *   \f]
  *
  *   with the design variable \( x \) constrained as follows:
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
  *   In the binary case, if \( \mathrm{MinCapacityDesign} > 0 \) then \( x = 1 \).
  */

 void generate_abstract_constraints( Configuration * stcc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the objective of the IntermittentUnitBlock
 /** Method that generates the objective of the IntermittentUnitBlock. The
  * objective can include:
  *
  * - a linear term on active power, if the vector "ActivePowerCost" is
  *   provided:
  *   \f[
  *     \min \ \sum_t B[t] \cdot p^{ac}_t
  *   \f]
  *   (coefficients are also scaled by the Block scale factor, if any);
  *
  * - an investment term in design mode, if "InvestmentCost" \f$ \ne 0 \f$:
  *   \f$ + \ I \cdot x \f$.
  *
  * Hence, in the design scenario the full objective is:
  * \f[
  *   \min \ \sum_t B[t] \cdot p^{ac}_t \;+\; I \cdot x \; .
  * \f]
  * If "ActivePowerCost" is not provided, \f$ B[t] = 0 \f$ and only the
  * investment term remains in design mode.
  */

 void generate_objective( Configuration * objc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// setting the BlockConfig
 /** This method sets the BlockConfig of this IntermittentUnitBlock. Besides
  * the Configuration for the is_feasible() function, the
  * IntermittentUnitBlock also considers the extra Configuration of the
  * BlockConfig. If the extra Configuration is a non-null pointer to a
  * SimpleConfiguration< double >, then the value, let us call it
  * \f$ \varepsilon \f$, stored in that Configuration will replace any zero
  * value that may appear as maximum power at any time instant.
  *
  * For instance, if the maximum power provided during deserialization
  * (see IntermittentUnitBlock::deserialize( netCDF::NcGroup )) is zero for
  * some time instant \f$ t \f$, then it will become \f$ \varepsilon \f$ for
  * that time instant. Moreover, if any zero value is provided to
  * set_maximum_power() for some time instant \f$ t \f$, then the maximum
  * power for time instant \f$ t \f$ will become \f$ \varepsilon \f$.
  *
  * When \f$ \varepsilon > 0 \f$, this can be used to prevent the maximum
  * power from being zero. Notice, however, that the actual maximum power may
  * become zero even if \f$ \varepsilon > 0 \f$ if the kappa constant is zero
  * (see set_kappa()).
  *
  * The reason behind this is that some Solver may not be able to handle
  * modifications in the maximum power if it is initially zero and becomes
  * nonzero after a modification. By setting \f$ \varepsilon > 0 \f$, this
  * issue is avoided.
  *
  * Please see the comments to Block::set_BlockConfig() for more details
  * about the BlockConfig.
  */

 void set_BlockConfig( BlockConfig * newBC = nullptr ,
                       bool deleteold = true ) override;

/**@} ----------------------------------------------------------------------*/
/*------------- Methods for checking the IntermittentUnitBlock -------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for checking solution information in the
 *        IntermittentUnitBlock
 * @{ */

 /// returns true if the current solution is (approximately) feasible
 /** This function returns true if and only if the solution encoded in the
  * current value of the Variable of this IntermittentUnitBlock is
  * approximately feasible within the given tolerance. That is, a solution
  * is considered feasible if and only if
  *
  * -# each ColVariable is feasible; and
  *
  * -# the violation of each Constraint of this IntermittentUnitBlock is not
  *    greater than the tolerance.
  *
  * Every Constraint of this IntermittentUnitBlock is a RowConstraint and its
  * violation is given by either the relative (see RowConstraint::rel_viol())
  * or the absolute violation (see RowConstraint::abs_viol()), depending on
  * the Configuration that is provided.
  *
  * The tolerance and the type of violation can be provided by either \p fsbc
  * or #f_BlockConfig->f_is_feasible_Configuration, and they are determined as
  * follows:
  *
  * - If \p fsbc is not nullptr, and it is a pointer to a
  *   SimpleConfiguration< double >, then the tolerance is the value present
  *   in that SimpleConfiguration and the relative violation is considered.
  *
  * - If \p fsbc is not nullptr, and it is a pointer to a
  *   SimpleConfiguration< std::pair< double , int > >, then the tolerance is
  *   fsbc->f_value.first and the type of violation is determined by
  *   fsbc->f_value.second (any nonzero number for relative violation and
  *   zero for absolute violation);
  *
  * - Otherwise, if both #f_BlockConfig and
  *   f_BlockConfig->f_is_feasible_Configuration are not nullptr and the
  *   latter is a pointer to either a SimpleConfiguration< double > or to a
  *   SimpleConfiguration< std::pair< double , int > >, then the values of
  *   the parameters are obtained as above;
  *
  * - Otherwise, by default, the tolerance is 0 and the relative violation is
  *   considered.
  *
  * This function currently considers only the abstract representation to
  * determine if the solution is feasible. So, the parameter \p useabstract
  * is currently ignored. If no abstract Variable has been generated, then
  * this function returns true. Moreover, if no abstract Constraint has been
  * generated, the solution is considered to be feasible with respect to the
  * set of Variables only. Notice also that, before checking if the solution
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

/**@} ----------------------------------------------------------------------*/
/*------- METHODS FOR READING THE DATA OF THE IntermittentUnitBlock --------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the IntermittentUnitBlock
 *
 * These methods allow reading data that must be common to (in principle)
 * all kinds of Intermittent Generation units
 * @{ */

 /// returns the gamma value
 double get_gamma( void ) const { return( f_gamma ); }

 /// returns the kappa value
 double get_kappa( void ) const { return( f_kappa ); }

 /// returns the investment cost
 double get_investment_cost( void ) const { return( f_InvestmentCost ); }

 /// returns the maximum installable capacity by the user
 double get_max_capacity_design( void ) const { return( f_MaxCapacityDesign ); }

/*--------------------------------------------------------------------------*/
 /// returns the minimum power of \p generator at time \p t

 double get_min_power( Index t , Index generator = 0 ) const override {
  return( ( v_MinPower.size() > t ) ? v_MinPower[ t ] : 0 );
  }

/*--------------------------------------------------------------------------*/
 /// returns the maximum power of \p generator at time \p t

 double get_max_power( Index t , Index generator = 0 ) const override {
  return( ( v_MaxPower.size() > t ) ? v_MaxPower[ t ] : 0 );
  }

/*--------------------------------------------------------------------------*/
 /// returns the minimum reactive power of \p generator at time \t

 double get_min_reactive_power( Index t , Index generator = 0 )
  const override {
  return( ( v_MinReactivePower.size() > t ) ? v_MinReactivePower[ t ] : 0 );
  }

/*--------------------------------------------------------------------------*/
 /// returns the maximum reactive power of \p generator at time \t

 double get_max_reactive_power( Index t , Index generator = 0 )
  const override {
  return( ( v_MaxReactivePower.size() > t ) ? v_MaxReactivePower[ t ] : 0 );
  }

/*--------------------------------------------------------------------------*/
 /// returns the matrix of inertia power
 /** The returned value \f$ U = \mathrm{get\_inertia\_power()} \f$ contains
  * the contribution to inertia (basically, the constants to be multiplied
  * by the active power variables returned by get_active_power()) of all the
  * generators at all time instants. There are four possible cases:
  *
  * - if the matrix is empty, then the inertia power is always 0 and this
  *   function returns nullptr;
  * - if the matrix only has one row (i.e., the first dimension has size 1),
  *   then the inertia power for each generator \f$ g \f$ is \f$ U[0,g] \f$
  *   for all \f$ t \f$ which means that the second dimension has size
  *   get_number_generators();
  * - if the matrix only has one column with size get_time_horizon() (i.e.,
  *   the second dimension has size 1), then
  *   \f$ \mathrm{InertiaPower}[t,0] \f$ gives the inertia power for the
  *   problem at time \f$ t \f$. Since in this unit there is only one
  *   electrical generator, this case should happen by assumption;
  * - otherwise, the matrix has size get_time_horizon() per
  *   get_number_generators(), then \f$ \mathrm{InertiaPower}[t,g] \f$
  *   represents the inertia power for the problem at time \f$ t \f$ for
  *   each electrical generator \f$ g \f$. */

 const double * get_inertia_power( Index generator ) const override {
  if( v_InertiaPower.empty() )
   return( nullptr );
  return( &( v_InertiaPower.front() ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of active power cost
 /** The returned vector contains the active power cost at time \f$ t \f$.
  * There are three possible cases:
  *
  * - if the vector is empty, then the linear power cost of the unit is 0;
  * - if the vector has only one element, then that element is the active
  *   power cost for all the time horizon;
  * - otherwise, the vector must have size get_time_horizon() and each
  *   element of the vector represents the active power cost at time 
  *   \f$ t \f$. */

 const std::vector< double > & get_active_power_cost( void ) const {
  return( v_ActivePowerCost );
  }

/*--------------------------------------------------------------------------*/
 /// returns the coefficient of the active power cost function
 /** This function returns the coefficient of the linear active-power cost
  * term that represents the cost of the power produced by the unit at the
  * given time instant.
  *
  * @param t A time instant between 0 and get_time_horizon() - 1.
  *
  * @return The coefficient of the active power cost of the linear function
  *         that represents the cost of the power produced by the unit at
  *         the given time instant. */

 double get_active_power_cost( Index t ) const {
  if( v_ActivePowerCost.empty() )
   return( 0 );
  if( v_ActivePowerCost.size() == 1 )
   return( v_ActivePowerCost.front() );
  assert( v_ActivePowerCost.size() == f_time_horizon );
  if( t >= f_time_horizon )
   throw( std::logic_error(
    "IntermittentUnitBlock::get_active_power_cost: Invalid time index: " +
    std::to_string( t ) ) );
  return( v_ActivePowerCost[ t ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the scale factor

 double get_scale( void ) const override { return( f_scale ); }

/** @} ---------------------------------------------------------------------*/
/*------ METHODS FOR READING THE Variable OF THE IntermittentUnitBlock -----*/
/*--------------------------------------------------------------------------*/
/** @name Reading the Variable of the IntermittentUnitBlock
 *
 * These methods allow to read each group of Variable that any
 * IntermittentUnitBlock in principle has (although some may not):
 *
 * - active and reactive power variables;
 *
 * - primary_spinning_reserve variables;
 *
 * - secondary_spinning_reserve variables.
 * @{ */

 /// returns the vector of active_power variables

 ColVariable * get_active_power( Index generator ) override {
  if( v_active_power.empty() )
   return( nullptr );
  return( &( v_active_power.front() ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of reactive power variables

 ColVariable * get_reactive_power( Index generator ) override {
  if( v_reactive_power.empty() )
   return( nullptr );
  return( &( v_reactive_power.front() ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of primary_spinning_reserve variables

 ColVariable * get_primary_spinning_reserve( Index generator ) override {
  if( v_primary_spinning_reserve.empty() )
   return( nullptr );
  return( &( v_primary_spinning_reserve.front() ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of secondary_spinning_reserve variables

 ColVariable * get_secondary_spinning_reserve( Index generator ) override {
  if( v_secondary_spinning_reserve.empty() )
   return( nullptr );
  return( &( v_secondary_spinning_reserve.front() ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the design variable

 ColVariable & get_design( void ) { return( design ); }

/*--------------------------------------------------------------------------*/
 /// returns the const design variable

 const ColVariable & get_const_design( void ) const { return( design ); }

/*--------------------------------------------------------------------------*/
 /// returns the minimum total power constraints

 const std::vector< FRowConstraint > &
 get_min_power_constraints( void ) const { return( min_power_Const ); }

/*--------------------------------------------------------------------------*/
 /// returns the maximum total power constraints
 
 const std::vector< FRowConstraint > &
 get_max_power_constraints( void ) const { return( max_power_Const ); }

/*--------------------------------------------------------------------------*/
 /// returns the bound constraints on the active power
 
 const std::vector< BoxConstraint > &
 get_active_power_bound_constraints( void ) const {
  return( active_power_bounds_Const );
  }

/** @} ---------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Solution
 * @{ */

 /// returns a Solution storing for this IntermittentUnitBlock
 /** This method must construct and return a (pointer to a) Solution object
  * representing the current "solution state" of this IntermittentUnitBlock.
  * This is a IntermittentUnitBlockSolution extending UnitBlockSolution with
  * the specific extra solution information of BatteryUnitBlock.
  *
  * The parameter for deciding which kind of Solution must be returned is a
  * single int value, coded bitwise:
  *
  * - the first four bits (bit 0 to bit 3) are "taken" by the base
  *   UnitBlock[Solution]
  *
  * This value is to be found as:
  *
  * - if solc is not nullptr, and it is a SimpleConfiguration< int >, then it
  *   is solc->f_value;
  *
  * - otherwise, if f_BlockConfig is not nullptr,
  *   f_BlockConfig->f_solution_Configuration is not nullptr, and it is a
  *   SimpleConfiguration< int >, then it is
  *   f_BlockConfig->f_solution_Configuration->f_value;
  *
  * - otherwise, it is 15 (save everything). */

 Solution * get_Solution( Configuration * solc = nullptr ,
                          bool emptys = true ) override;

/*--------------------------------------------------------------------------*/
 /// return the "appropriate" [Intermittent]UnitBlockSolution

 UnitBlockSolution * new_Solution( void ) const override;

/** @} ---------------------------------------------------------------------*/
/*-------------- METHODS FOR SAVING THE IntermittentUnitBlock---------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for printing & saving the IntermittentUnitBlock
 * @{ */

/// extends Block::serialize( netCDF::NcGroup )
/** Extends Block::serialize( netCDF::NcGroup ) to the specific format of an
 * IntermittentUnitBlock. See IntermittentUnitBlock::deserialize(
 * netCDF::NcGroup ) for details of the format of the created netCDF group.
 */
 void serialize( netCDF::NcGroup & group ) const override;

/** @} ---------------------------------------------------------------------*/
/*----------- METHODS FOR INITIALIZING THE IntermittentUnitBlock -----------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the data of the IntermittentUnitBlock
 * @{ */

 void load( std::istream & input , char frmt = 0 ) override {
  throw( std::logic_error(
   "IntermittentUnitBlock::load() not implemented yet" ) );
 }

/** @} ---------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

 void set_maximum_power( MF_dbl_it values ,
                         Subset && subset ,
                         const bool ordered = false ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck );

 void set_maximum_power( MF_dbl_it values ,
                         Range rng = Range( 0 , Inf< Index >() ) ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// sets the scale factor
 /** This method sets the scale factor.
  *
  * @param values An iterator to a vector containing the scale factor.
  *
  * @param subset If non-empty, the scale factor is set to the value pointed
  *               by \p values. If empty, no operation is performed.
  *
  * @param ordered This parameter is ignored.
  * @param issuePMod Controls how physical Modifications are issued.
  * @param issueAMod Controls how abstract Modifications are issued. */

 void scale( MF_dbl_it values ,
             Subset && subset , bool ordered = false ,
             c_ModParam issuePMod = eNoBlck ,
             c_ModParam issueAMod = eNoBlck ) override;

/*--------------------------------------------------------------------------*/
 /// set the kappa constant
 /** This function sets the kappa constant, which multiplies the minimum and
  * maximum power in the constraints of this IntermittentUnitBlock.
  *
  * @param values  Iterator to a vector containing the kappa constants.
  * @param subset  If non-empty, the kappa constant is set to the value
  *                pointed by \p values. If empty, no operation is performed.
  * @param ordered It indicates whether \p subset is ordered.
  * @param issuePMod Controls how physical Modifications are issued.
  * @param issueAMod Controls how abstract Modifications are issued. */

 void set_kappa( MF_dbl_it values ,
                 Subset && subset , bool ordered = false ,
                 c_ModParam issuePMod = eNoBlck ,
                 c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// set the kappa constant
 /** This function sets the kappa constant, which multiplies the minimum and
  * maximum power in the constraints of this IntermittentUnitBlock.
  *
  * @param values Iterator to a vector containing the kappa constants.
  * @param rng    If non-empty, the kappa constant is set to the value
  *               pointed by \p values. If empty, no operation is performed.
  * @param issuePMod Controls how physical Modifications are issued.
  * @param issueAMod Controls how abstract Modifications are issued.
  */
 void set_kappa( MF_dbl_it values ,
                 Range rng = Range( 0 , Inf< Index >() ) ,
                 c_ModParam issuePMod = eNoBlck ,
                 c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// set the kappa constant
 /** This function sets the kappa constant, which multiplies the minimum and
  * maximum power in the constraints of this IntermittentUnitBlock.
  *
  * @param value     The value of the kappa constant.
  * @param issuePMod Controls how physical Modifications are issued.
  * @param issueAMod Controls how abstract Modifications are issued. */

 void set_kappa( double value , c_ModParam issuePMod = eNoBlck ,
                                c_ModParam issueAMod = eNoBlck ) {
  std::vector< double > vector = { value };
  set_kappa( vector.cbegin() , Range( 0 , Inf< Index >() ) ,
             issuePMod , issueAMod );
  }

/*--------------------------------------------------------------------------*/
 // For the Range version, use the default implementation defined in UnitBlock

 using UnitBlock::scale;

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

 /// the vector of MinPower
 std::vector< double > v_MinPower;

 /// the vector of MaxPower
 std::vector< double > v_MaxPower;

 /// the matrix of inertia power of generators
 std::vector< double > v_InertiaPower;

 /// the vector of ActivePowerCost
 std::vector< double > v_ActivePowerCost;

 /// the vector of MinReactivePower
 std::vector< double > v_MinReactivePower;

 /// the vector of MaxReactivePower
 std::vector< double > v_MaxReactivePower;

 /// the investment cost
 double f_InvestmentCost;

 /// the minimum capacity design allowed
 double f_MinCapacityDesign;

 /// the maximum capacity design allowed
 double f_MaxCapacityDesign;

 /// the gamma value
 double f_gamma;

 /// the kappa value
 double f_kappa;

 /// the scale factor
 double f_scale;

 /// value used to replace any zero value in the maximum power
 double f_max_power_epsilon;

/*-------------------------------- variables -------------------------------*/

 /// the active power variables
 std::vector< ColVariable > v_active_power;

 /// the reactive power variables
 std::vector< ColVariable > v_reactive_power;

 /// the primary spinning reserve variables
 std::vector< ColVariable > v_primary_spinning_reserve;

 /// the secondary spinning reserve variables
 std::vector< ColVariable > v_secondary_spinning_reserve;

 /// the design variable
 ColVariable design;

/*------------------------------- constraints ------------------------------*/

 /// the active power lower bound constraints
 std::vector< FRowConstraint > min_power_Const;

 /// the active power upper bound constraints
 std::vector< FRowConstraint > max_power_Const;

 /// the active power bounds design constraints
 boost::multi_array< FRowConstraint , 2 > active_power_bounds_design_Const;

 /// the active power bounds constraints
 std::vector< BoxConstraint > active_power_bounds_Const;

 /// the reactive power bound constraints
 std::vector< BoxConstraint > ReactivePower_Bound_Const;

 /*!! Q <= P
 std::vector< FRowConstraint > Reactive_2_Active_Const;
 !!*/

 /// the design bound constraint
 BoxConstraint design_bound_Const;

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

 /// updates the constraints for the current maximum power
 /** This function updates the right-hand side of the "maximum power" and
  * the "active power bounds" constraints associated with the time instants
  * given in \p time.
  */
 void update_max_power_in_cnstrs( const Block::Subset & time ,
                                  c_ModParam issueAMod );

/*--------------------------------------------------------------------------*/
 /// updates the constraints for the current maximum power
 /** This function updates the right-hand side of the "maximum power" and
  * the "active power bounds" constraints associated with the time instants
  * given in \p time.
  */
 void update_max_power_in_cnstrs( const Block::Range & time ,
                                  c_ModParam issueAMod );

/*--------------------------------------------------------------------------*/
 /// verify whether the data in this IntermittentUnitBlock is consistent
 /** This function checks whether the data in this IntermittentUnitBlock is
  * consistent. The data is consistent if all the following conditions are met.
  *
  * - The maximum power is greater than or equal to the minimum power.
  *
  * - The minimum power is nonnegative.
  *
  * - \f$ 0 \le \Gamma \le 1 \f$.
  *
  * - \f$ \kappa \ge 0 \f$.
  *
  * - The inertia power is nonnegative.
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
   * three-args version of register_method<> with the compact
   * MS_*_*::args(),
   *
   * register_method< IntermittentUnitBlock >
   *                ( "IntermittentUnitBlock::set_maximum_power",
   *                  &IntermittentUnitBlock::set_maximum_power,
   *                  MS_dbl_sbst::args() );
   *
   * so we just use the slightly less compact one with the explicit argument
   * and be done with it.
   */

  register_method< IntermittentUnitBlock , MF_dbl_it , Subset && , bool >(
   "IntermittentUnitBlock::set_maximum_power" ,
   &IntermittentUnitBlock::set_maximum_power );

  register_method< IntermittentUnitBlock , MF_dbl_it , Range >(
   "IntermittentUnitBlock::set_maximum_power" ,
   &IntermittentUnitBlock::set_maximum_power );

  register_method< IntermittentUnitBlock , MF_dbl_it , Subset && , bool >(
   "IntermittentUnitBlock::scale" ,
   &IntermittentUnitBlock::scale );

  register_method< IntermittentUnitBlock , MF_dbl_it , Range >(
   "IntermittentUnitBlock::scale" ,
   &IntermittentUnitBlock::scale );

  register_method< IntermittentUnitBlock , MF_dbl_it , Subset && , bool >(
   "IntermittentUnitBlock::set_kappa" ,
   &IntermittentUnitBlock::set_kappa );

  register_method< IntermittentUnitBlock , MF_dbl_it , Range >(
   "IntermittentUnitBlock::set_kappa" ,
   &IntermittentUnitBlock::set_kappa );
 }

};  // end( class( IntermittentUnitBlock ) )

/*--------------------------------------------------------------------------*/
/*--------------------- CLASS IntermittentUnitBlockMod ---------------------*/
/*--------------------------------------------------------------------------*/

/// derived class from Modification for changes to an IntermittentUnitBlock

class IntermittentUnitBlockMod : public UnitBlockMod
{
 public:

 /// public enum for the types of IntermittentUnitBlockMod
 enum IUB_mod_type
 {
  eSetMaxP = eUBModLastParam , ///< set max power values
  eSetKappa ,                  ///< set the kappa constant
  eIUBModLastParam         ///< first allowed parameter for derived classes
  /**< Convenience value to easily allow derived classes to extend the set
   * of types of IntermittentUnitBlockMod. */
  };

 /// constructor, takes the IntermittentUnitBlock and the type
 IntermittentUnitBlockMod( IntermittentUnitBlock * const fblock ,
                           const int type )
  : UnitBlockMod( fblock , type ) {}

 /// destructor, does nothing
 virtual ~IntermittentUnitBlockMod() override = default;

 /// returns the Block to which the Modification refers
 Block * get_Block( void ) const override { return( f_Block ); }

 protected:

 /// prints the IntermittentUnitBlockMod
 void print( std::ostream & output ) const override {
  output << "IntermittentUnitBlockMod[" << this << "]: ";
  switch( f_type ) {
   default:
    output << "Set max power values ";
  }
 }

 IntermittentUnitBlock * f_Block{};
 ///< pointer to the Block to which the Modification refers

};  // end( class( IntermittentUnitBlockMod ) )

/*--------------------------------------------------------------------------*/
/*------------------- CLASS IntermittentUnitBlockRngdMod -------------------*/
/*--------------------------------------------------------------------------*/

/// derived from IntermittentUnitBlockMod for "ranged" modifications
class IntermittentUnitBlockRngdMod : public IntermittentUnitBlockMod
{

 public:

 /// constructor: takes the IntermittentUnitBlock, the type, and the range
 IntermittentUnitBlockRngdMod( IntermittentUnitBlock * const fblock ,
                               const int type , const Block::Range & rng )
  : IntermittentUnitBlockMod( fblock , type ) , f_rng( rng ) {}

 /// destructor, does nothing
 virtual ~IntermittentUnitBlockRngdMod() override = default;

 /// accessor to the range
 Block::c_Range & rng( void ) { return( f_rng ); }

 protected:

 /// prints the IntermittentUnitBlockRngdMod
 void print( std::ostream & output ) const override {
  IntermittentUnitBlockMod::print( output );
  output << "[ " << f_rng.first << ", " << f_rng.second << " )"
         << std::endl;
 }

 Block::Range f_rng;  ///< the range

};  // end( class( IntermittentUnitBlockRngdMod ) )

/*--------------------------------------------------------------------------*/
/*------------------- CLASS IntermittentUnitBlockSbstMod -------------------*/
/*--------------------------------------------------------------------------*/

/// derived from IntermittentUnitBlockMod for "subset" modifications
class IntermittentUnitBlockSbstMod : public IntermittentUnitBlockMod
{

 public:

 /// constructor: takes the IntermittentUnitBlock, the type, and the subset
 IntermittentUnitBlockSbstMod( IntermittentUnitBlock * const fblock ,
                               const int type , Block::Subset && nms )
  : IntermittentUnitBlockMod( fblock , type ) , f_nms( std::move( nms ) ) {}

 /// destructor, does nothing
 virtual ~IntermittentUnitBlockSbstMod() override = default;

 /// accessor to the subset
 Block::c_Subset & nms( void ) { return( f_nms ); }

 protected:

 /// prints the IntermittentUnitBlockSbstMod
 void print( std::ostream & output ) const override {
  IntermittentUnitBlockMod::print( output );
  output << "(# " << f_nms.size() << ")" << std::endl;
 }

 Block::Subset f_nms;  ///< the subset

 };  // end( class( IntermittentUnitBlockSbstMod ) )

/*--------------------------------------------------------------------------*/
/*------------------ CLASS IntermittentUnitBlockSolution -------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// a [UnitBlock]Solution of a IntermittentUnitBlock
/** The IntermittentUnitBlockSolution class derives from UnitBlockSolution and
 * adds to the "standard" information stored in there (active power, possibly
 * commitment and primary/secondary reserve) the other information that is
 * typical of the IntermittentUnitBlock, i.e.,
 *
 * - if defined, the value of the Intermittent Design Variable */

class IntermittentUnitBlockSolution : public UnitBlockSolution
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*----------------------------- CONSTANTS ----------------------------------*/

 static constexpr double dNaN = std::numeric_limits< double >::quiet_NaN();
 ///< convenience constexpr for "NaN", *not* to be used with ==

/*------------------------------- FRIENDS ----------------------------------*/

 friend IntermittentUnitBlock;  ///< make IntermittentUnitBlock friend

/*------- CONSTRUCTING AND DESTRUCTING IntermittentUnitBlockSolution -------*/

 /// constructor, it has nothing to do
 explicit IntermittentUnitBlockSolution( void ) : f_design( dNaN ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void deserialize( const netCDF::NcGroup & group ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 ~IntermittentUnitBlockSolution() override = default;
 ///< destructor: it is virtual, and empty

/*--- METHODS DESCRIBING THE BEHAVIOR OF A IntermittentUnitBlockSolution --*/

 void read( const Block * block ) override final;

 void write( Block * block ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// serialize a IntermittentUnitBlockSolution into a netCDF::NcGroup
 /** Serialize a IntermittentUnitBlockSolution into a netCDF::NcGroup.
  * The format is the one of UnitBlockSolution
  * [cf. UnitBlockSolution::serialize()], plus:
  *
  * - The scalar variable "IntermittentDesign", of type netCDF::NcDouble,
  *   that represent the value of the dimensioning variable; the variable
  *   is optional in that the intermittent unit may not have any
  *   dimensioning variable. */

 void serialize( netCDF::NcGroup & group ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 IntermittentUnitBlockSolution * scale( double factor ) const override;

 void sum( const Solution * solution , double multiplier ) override;

 IntermittentUnitBlockSolution * clone( bool empty = false ) const override;

/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/

 void print( std::ostream & output ) const override {
  output << "IntermittentUnitBlockSolution [" << this << "]: " << std::endl;
  }

/*---------------------- PRIVATE PART OF THE CLASS -------------------------*/

 private:

/*---------------------------- PRIVATE FIELDS ------------------------------*/

 double f_design;    ///< the value of the dimensioning variable

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class( IntermittentUnitBlockSolution ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* __IntermittentUnitBlock */

/*--------------------------------------------------------------------------*/
/*------------------ End File IntermittentUnitBlock.h ----------------------*/
/*--------------------------------------------------------------------------*/
