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
/// implementation of the Block concept for an intermittent generation unit
/** IntermittentUnitBlock implements the Block concept [see Block.h] for a
 * unit whose production depends on a primary source that cannot be stored and
 * is only forecast, such as a wind farm, a solar park or a run-of-the-river
 * plant. Such a unit is connected to one node of the network, which may be a
 * node of the transmission network (centralized generation) or a node of a
 * distribution network (distributed generation, e.g., photovoltaic panels on
 * roofs). It may also be a node of an aggregated network that represents a
 * zone, in which case the unit represents all the capacity of its kind
 * installed in the zone, and the equations are the same in all cases. It has
 * one generator, and its place in the model of the complete system (the node
 * balance, the reserve and inertia rows, the scale factor and the investment)
 * is described in \ref ucblock_model.
 *
 * \par Capacity and production
 * We index the instants by \f$ t \in \mathcal{T} = \{ 0 , \ldots , T - 1
 * \} \f$, and each datum is given per instant. We intend the datum
 * \f$ P^{mx}_t \f$ ("MaxPower") as \f$ P^{mx}_t = \kappa^{mx} L_t \f$, the
 * product of the capacity \f$ \kappa^{mx} \f$ that can be installed and of
 * the load factor \f$ L_t \in [ 0 , 1 ] \f$ at instant \f$ t \f$, i.e., the
 * production of a unit of capacity, which is typically obtained from
 * historical series of the source at the node. Accordingly, the datum
 * \f$ \kappa \geq 0 \f$ ("Kappa", 1 by default) is the installed fraction of
 * \f$ \kappa^{mx} \f$, and the unit can produce at most \f$ \kappa P^{mx}_t =
 * \kappa \kappa^{mx} L_t \f$; usually \f$ \kappa \in [ 0 , 1 ] \f$, but
 * larger values are accepted, and \f$ \kappa \f$ multiplies \f$ P^{mn}_t \f$
 * ("MinPower") as well. The active power \f$ p^{ac}_t \f$ ("p_intermittent")
 * lies in \f$ [ \kappa P^{mn}_t , \kappa P^{mx}_t ] \f$ (see (3) below), and
 * the difference \f$ \kappa P^{mx}_t - p^{ac}_t \f$ is the curtailed
 * production; with \f$ P^{mn}_t = P^{mx}_t \f$ the unit can be neither
 * curtailed nor used for reserve at \f$ t \f$. Usually, the datum
 * \f$ P^{mn}_t \f$ is nonnegative; a negative one is accepted, though, and it
 * makes the unit a sink that absorbs a variable power
 * (\f$ P^{mn}_t < 0 = P^{mx}_t \f$), in which case the active power is a free
 * variable whose sign comes from (3) only. The energy produced over the
 * horizon may be bounded by \f$ G^{mn} \f$ ("MinGeneration") and
 * \f$ G^{mx} \f$ ("MaxGeneration") in (4), and these are absolute energies,
 * multiplied neither by \f$ \kappa \f$ nor by the design variable.
 *
 * \par Reserves
 * The primary and secondary reserves \f$ p^{pr}_t , p^{sc}_t \geq 0 \f$
 * ("pr_intermittent", "sr_intermittent") are symmetric, i.e., each must be
 * available both upwards and downwards. Downwards the margin
 * \f$ p^{ac}_t - \kappa P^{mn}_t \f$ is fully available (see (2)). Upwards,
 * since \f$ \kappa P^{mx}_t \f$ is only a forecast, only the fraction
 * \f$ \gamma \in [ 0 , 1 ] \f$ ("Gamma") of the margin
 * \f$ \kappa P^{mx}_t - p^{ac}_t \f$ is (see (1)), and \f$ 1 - \gamma \f$ is
 * a safety margin against the error of the forecast. With \f$ \gamma = 0 \f$,
 * the default, the unit gives no reserve and the reserve variables do not
 * exist; otherwise each of them exists if the enclosing UCBlock has primary
 * (respectively, secondary) zones.
 *
 * \par Inertia and cost
 * The unit contributes \f$ h^p_t p^{ac}_t \f$ to the inertia rows of UCBlock,
 * with \f$ h^p_t \geq 0 \f$ ("InertiaPower", 0 by default). It is thus
 * proportional to the dispatched power, and a curtailed unit contributes less
 * (this is mostly meant for run-of-the-river plants). As for the Objective,
 * it is linear, with a price \f$ b_t \f$ ("ActivePowerCost") per unit of
 * energy produced (see (5)).
 *
 * \par Investment
 * One can choose the capacity in two ways, which are not meant to be
 * combined. One is \f$ \kappa \f$, changed from outside the Block (see
 * set_kappa() and get_kappa_linearization()); the other is a design variable
 * \f$ x \f$ ("x_intermittent"), which exists if \f$ c^{inv} \f$
 * ("InvestmentCost") is nonzero, multiplies \f$ \kappa P^{mn}_t \f$ and
 * \f$ \kappa P^{mx}_t \f$ in (1)-(3), which then take the forms (1d)-(3d),
 * and costs \f$ \sigma c^{inv} \f$ per unit. If \f$ \bar{X} \geq 0 \f$ the
 * design variable lies in \f$ [ \max\{ 0 , \underline{X} \} , \bar{X} ] \f$,
 * with \f$ \underline{X} \f$ and \f$ \bar{X} \f$ the data "MinCapacityDesign"
 * and "MaxCapacityDesign". If \f$ \bar{X} < 0 \f$ it is integer in \f$ \{
 * \max\{ 0 , \underline{X} \} , \ldots , | \bar{X} | \} \f$ (say, a number of
 * modules), and binary if \f$ \bar{X} = -1 \f$; its bounds are always a row,
 * and therefore with both bounds equal to 1 the unit is built. Finally, the
 * scale factor \f$ \sigma \f$ ("Scale") represents \f$ \sigma \f$ identical
 * copies of the unit that share its Variable (see UnitBlock::scale()).
 *
 * \par Rows and Objective
 * With the terms in \f$ p^{pr}_t \f$ and \f$ p^{sc}_t \f$ present only for
 * the reserves that exist, the rows are, for all \f$ t \in \mathcal{T} \f$,
 * \f{align*}{
 *   & p^{pr}_t + p^{sc}_t \leq \gamma ( \kappa P^{mx}_t - p^{ac}_t )
 *     \tag{1} \\
 *   & p^{pr}_t + p^{sc}_t \leq p^{ac}_t - \kappa P^{mn}_t \tag{2} \\
 *   & \kappa P^{mn}_t \leq p^{ac}_t \leq \kappa P^{mx}_t \tag{3} \\
 *   & G^{mn} \leq \sum_{ t \in \mathcal{T} } p^{ac}_t \leq G^{mx} \tag{4}
 * \f}
 * where (1) and (2) exist only if some reserve does, and (4), which is one
 * row, only if \f$ G^{mn} > - \infty \f$ or \f$ G^{mx} < \infty \f$. With
 * the design variable the rows (1)-(3) are
 * \f{align*}{
 *   & p^{pr}_t + p^{sc}_t \leq \gamma ( \kappa P^{mx}_t x - p^{ac}_t )
 *     \tag{1d} \\
 *   & p^{pr}_t + p^{sc}_t \leq p^{ac}_t - \kappa P^{mn}_t x \tag{2d} \\
 *   & \kappa P^{mn}_t x \leq p^{ac}_t \leq \kappa P^{mx}_t x \tag{3d}
 * \f}
 * and the Objective is
 * \f[
 *   \sigma \Bigl( \sum_{ t \in \mathcal{T} } b_t p^{ac}_t + c^{inv} x
 *   \Bigr) , \tag{5}
 * \f]
 * where the term in \f$ x \f$ is present only if the design variable exists.
 *
 * \par Features not modeled
 * We do not model a contribution to the inertia proportional to the potential
 * production \f$ \kappa P^{mx}_t \f$, independent of the dispatch, since it
 * would be a constant that plays no role in the dispatch of the unit. The
 * forecast \f$ P^{mx}_t \f$ is a datum, whose uncertainty is represented
 * within one scenario only by \f$ \gamma \f$, while several scenarios are
 * handled outside the Block, e.g., by a stochastic Block that changes
 * "MaxPower" with set_maximum_power(). Also, the active and the reactive
 * power are not linked. Finally, a biomass plant is dispatchable, and it is
 * therefore represented by a ThermalUnitBlock (with zero emission factors if
 * its emissions are not counted) rather than by an IntermittentUnitBlock. */

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
  * "ChangeIntervals". The symbols are those of the class description.
  * Several data below are time-indexed: such a variable, of type
  * netCDF::NcDouble, is either a scalar, or indexed over the dimension
  * "NumberIntervals" (one value per interval, see "ChangeIntervals"), or
  * indexed over "TimeHorizon", and it is expanded into one value per instant
  * as described in UnitBlock::deserialize(). The netCDF::NcGroup must then
  * also contain:
  *
  * - the time-indexed variable "MaxPower", the maximum power
  *   \f$ P^{mx}_t = \kappa^{mx} L_t \f$ in (1) and (3); if the
  *   SimpleConfiguration< double > \f$ \varepsilon > 0 \f$ of the
  *   BlockConfig is set [see set_BlockConfig()], every zero in it is
  *   replaced by \f$ \varepsilon \f$;
  *
  * and it may contain:
  *
  * - the time-indexed variable "MinPower", the minimum power
  *   \f$ P^{mn}_t \leq P^{mx}_t \f$ in (2) and (3), 0 if absent; it is
  *   normally nonnegative, while a negative value makes the unit absorb
  *   power (see the class description);
  *
  * - the time-indexed variable "InertiaPower", the nonnegative coefficient
  *   \f$ h^p_t \f$ of the active power in the inertia rows of UCBlock, 0 if
  *   absent;
  *
  * - the time-indexed variable "ActivePowerCost", the price \f$ b_t \f$ of
  *   (5), 0 if absent;
  *
  * - the scalar variables "MinGeneration" and "MaxGeneration", the bounds
  *   \f$ G^{mn} \leq G^{mx} \f$ of (4) on the energy produced over the
  *   horizon, \f$ - \infty \f$ and \f$ + \infty \f$ if absent; they are
  *   absolute energies, not multiplied by \f$ \kappa \f$ nor by \f$ x \f$;
  *
  * - the scalar variable "Gamma", the fraction \f$ \gamma \in [ 0 , 1 ] \f$
  *   of the upward margin offered as reserve in (1); 0 if absent, in which
  *   case the unit gives no reserve;
  *
  * - the scalar variable "Kappa", the nonnegative \f$ \kappa \f$; 1 if
  *   absent;
  *
  * - the scalar variable "Scale", the scale factor \f$ \sigma \f$ [see
  *   UnitBlock::scale()], 1 if absent; a value different from 1 is refused
  *   together with \f$ | \bar{X} | > 1 \f$, since the two represent the same
  *   replication of the unit (by deserialize() only, scale() not repeating
  *   this check);
  *
  * - the scalar variable "InvestmentCost", the investment cost
  *   \f$ c^{inv} \f$ of (5), 0 if absent; a nonzero value creates the design
  *   variable \f$ x \f$;
  *
  * - only if "InvestmentCost" is given, the scalar variables
  *   "MinCapacityDesign" and "MaxCapacityDesign", the bounds
  *   \f$ \underline{X} \geq 0 \f$ (0 if absent) and \f$ \bar{X} \f$ (1 if
  *   absent) that give the domain of \f$ x \f$ described in the class
  *   description; it must be \f$ \underline{X} \leq \bar{X} \f$ if
  *   \f$ \bar{X} > 0 \f$, \f$ \underline{X} \leq 1 \f$ if
  *   \f$ | \bar{X} | = 1 \f$ and \f$ \underline{X} \leq | \bar{X} | \f$ if
  *   \f$ \bar{X} < 0 \f$;
  *
  * - the time-indexed variables "MinReactivePower" and "MaxReactivePower",
  *   the bounds of the reactive power, used only if the enclosing UCBlock
  *   asks for it [see UnitBlock::set_reactive_power()]; a vector of zeros is
  *   treated as absent.
  *
  * The scalar variable "MaxCapacity" may be present, and is ignored. */

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
 /** Generates the static Variable of the IntermittentUnitBlock, with the
  * symbols of the class description; each group is a
  * std::vector< ColVariable > of size get_time_horizon(), or empty if it is
  * not generated:
  *
  * - the active power \f$ p^{ac}_t \f$ ("p_intermittent"), continuous and
  *   free (its sign comes from (3));
  *
  * - the primary and secondary reserves \f$ p^{pr}_t \f$
  *   ("pr_intermittent") and \f$ p^{sc}_t \f$ ("sr_intermittent"),
  *   nonnegative, each only if \f$ \gamma \neq 0 \f$ and the enclosing
  *   UCBlock asks for that reserve [see set_reserve_vars()];
  *
  * - the reactive power ("q_intermittent"), nonnegative, only if the
  *   enclosing UCBlock asks for it [see set_reactive_power()];
  *
  * and the single design variable \f$ x \f$ ("x_intermittent"), only if
  * "InvestmentCost" is nonzero: nonnegative if \f$ \bar{X} \geq 0 \f$,
  * integer if \f$ \bar{X} < 0 \f$ (binary if \f$ \bar{X} = -1 \f$), its
  * bounds being set by generate_abstract_constraints(). The parameter
  * \p stvv is only passed to UnitBlock::generate_abstract_variables(). */

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the static constraints of the IntermittentUnitBlock
 /** Generates the static Constraint of the IntermittentUnitBlock, i.e., the
  * rows (1)-(4) of the class description, or (1d)-(3d) and (4) with the
  * design variable, under the conditions said there:
  *
  * - (1) is "MaxPower_Intermittent" and (2) is "MinPower_Intermittent", two
  *   vectors of FRowConstraint [see get_max_power_constraints() and
  *   get_min_power_constraints()], written as \f$ \gamma p^{ac}_t +
  *   p^{pr}_t + p^{sc}_t \leq \gamma \kappa P^{mx}_t \f$ and
  *   \f$ \kappa P^{mn}_t \leq p^{ac}_t - p^{pr}_t - p^{sc}_t \f$; with the
  *   design variable, (1d) and (2d) are the same groups, with
  *   \f$ - \gamma \kappa P^{mx}_t \f$ and \f$ - \kappa P^{mn}_t \f$ as
  *   coefficients of \f$ x \f$ (the latter possibly 0) and 0 as right-hand
  *   side;
  *
  * - (3) is "ActivePower_Intermittent", a vector of BoxConstraint [see
  *   get_active_power_bound_constraints()]; with the design variable (3d)
  *   is "ActivePower_Design_Intermittent", a multi_array of FRowConstraint
  *   with the upper rows \f$ p^{ac}_t - \kappa P^{mx}_t x \leq 0 \f$ in its
  *   last slice and, only if \f$ P^{mn}_t \neq 0 \f$ at some instant, the
  *   lower rows \f$ p^{ac}_t - \kappa P^{mn}_t x \geq 0 \f$ in the slice 0,
  *   while if \f$ P^{mn}_t = 0 \f$ at all instants the lower side is the
  *   bound \f$ p^{ac}_t \geq 0 \f$ in "ActivePower_Intermittent";
  *
  * - (4) is "MaxMinGeneration_intermittent", a single FRowConstraint;
  *
  * - the bounds of \f$ x \f$ are "DesignBound_Intermittent", a
  *   BoxConstraint, present whenever \f$ x \f$ is;
  *
  * - the bounds of the reactive power are "ReactivePowerBound_intermittent",
  *   a vector of BoxConstraint, if the reactive power exists and one of its
  *   bounds is given.
  *
  * The parameter \p stcc is not used. */

 void generate_abstract_constraints( Configuration * stcc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the objective of the IntermittentUnitBlock
 /** Generates the Objective (5) of the class description, to be minimized:
  * a LinearFunction with coefficient \f$ \sigma b_t \f$ on
  * \f$ p^{ac}_t \f$ and \f$ \sigma c^{inv} \f$ on \f$ x \f$, if the latter
  * exists. The parameter \p objc is not used. */

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
  * or f_BlockConfig->f_is_feasible_Configuration, and they are determined as
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
  * - Otherwise, if both f_BlockConfig and
  *   f_BlockConfig->f_is_feasible_Configuration are not nullptr and the
  *   latter is a pointer to either a SimpleConfiguration< double > or to a
  *   SimpleConfiguration< std::pair< double , int > >, then the values of
  *   the parameters are obtained as above;
  *
  * - Otherwise, by default, the tolerance is 0 and the relative violation is
  *   considered.
  *
  * This function considers only the abstract representation to
  * determine if the solution is feasible. So, the parameter \p useabstract
  * is ignored. If no abstract Variable has been generated, then
  * this function returns true. Moreover, if no abstract Constraint has been
  * generated, the solution is considered to be feasible with respect to the
  * set of Variables only. Notice also that, before checking if the solution
  * satisfies a Constraint, the Constraint is computed
  * (Constraint::compute()).
  *
  * @param useabstract This parameter is ignored.
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

 /// returns the fraction gamma of the upward margin offered as reserve
 double get_gamma( void ) const { return( f_gamma ); }

 /// returns the capacity multiplier kappa
 /** Returns \f$ \kappa \f$ ("Kappa"), which multiplies \f$ P^{mn}_t \f$ and
  * \f$ P^{mx}_t \f$ in (1)-(3) of the class description, and not the energy
  * bounds of (4). */

 double get_kappa( void ) const override { return( f_kappa ); }

 /// returns the linearization coefficient of the kappa-parametrized objective
 /** Since \f$ \kappa \f$ only appears in the right-hand sides of (1)-(3) of
  * the class description when no design variable exists, the optimal value
  * of a convex minimization problem containing this unit (e.g., its
  * continuous relaxation) is a convex function of \f$ \kappa \f$, and a
  * subgradient of it is
  * \f[
  *   \sum_{ t \in \mathcal{T} } \Bigl( P^{mn}_t \bigl( \mu^{(3),mn}_t +
  *   \mu^{(2)}_t \bigr) - P^{mx}_t \bigl( \mu^{(3),mx}_t +
  *   \gamma \mu^{(1)}_t \bigr) \Bigr) ,
  * \f]
  * where \f$ \mu^{(3),mn}_t , \mu^{(3),mx}_t \geq 0 \f$ are the absolute
  * values of the dual value of (3) on its lower and upper side (the side
  * being told by the sign of the dual value) and \f$ \mu^{(2)}_t ,
  * \mu^{(1)}_t \geq 0 \f$ those of (2) and (1), if they exist. Writing
  * \f$ P^{mx}_t = \kappa^{mx} L_t \f$, the term in \f$ \mu^{(3),mx}_t \f$
  * is \f$ - \kappa^{mx} L_t \mu^{(3),mx}_t \f$: a larger capacity can only
  * decrease the cost. This method returns that value for the dual solution
  * held by the Constraint. If the unit has a design variable the data
  * multiply it instead, and the method throws std::logic_error. */

 double get_kappa_linearization( void ) const override;

 /// returns the investment cost
 double get_investment_cost( void ) const { return( f_InvestmentCost ); }

 /// returns the upper bound of the design variable ("MaxCapacityDesign")
 double get_max_capacity_design( void ) const {
  return( f_MaxCapacityDesign );
  }

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
 /// returns the minimum reactive power of \p generator at time \p t

 double get_min_reactive_power( Index t , Index generator = 0 )
  const override {
  return( ( v_MinReactivePower.size() > t ) ? v_MinReactivePower[ t ] : 0 );
  }

/*--------------------------------------------------------------------------*/
 /// returns the maximum reactive power of \p generator at time \p t

 double get_max_reactive_power( Index t , Index generator = 0 )
  const override {
  return( ( v_MaxReactivePower.size() > t ) ? v_MaxReactivePower[ t ] : 0 );
  }

/*--------------------------------------------------------------------------*/
 /// returns the coefficients of the active power in the inertia rows
 /** Returns a pointer to the get_time_horizon() coefficients
  * \f$ h^p_t \f$ ("InertiaPower") of \f$ p^{ac}_t \f$ in the inertia rows of
  * UCBlock, the unit having one generator; the vector is filled with zeros
  * if the datum is absent, so that the pointer is never nullptr after
  * deserialize(). */

 const double * get_inertia_power( Index generator ) const override {
  if( v_InertiaPower.empty() )
   return( nullptr );
  return( &( v_InertiaPower.front() ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of the prices of the energy produced
 /** Returns the vector of the prices \f$ b_t \f$ ("ActivePowerCost") of (5)
  * of the class description, of size get_time_horizon() (zeros if the datum
  * is absent). */

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

 double get_design_ub( void ) const override {
  // f_MaxCapacityDesign < 0: integer design in {0, ..., |MaxCapacityDesign|};
  // f_MaxCapacityDesign >= 0: continuous design in [0, MaxCapacityDesign].
  return( std::abs( f_MaxCapacityDesign ) );
 }

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

 /// the unit provides the reserve only if it has any to give
 /** The enclosing UCBlock asking for the reserve is not enough, the unit has
  * to have some to give, i.e., a nonzero fraction of its active power: this
  * is the very condition with which the Variable are generated, said in
  * terms of the data alone. */

 bool has_primary_reserve( void ) const override {
  return( ( reserve_vars & 1u ) && ( f_gamma != 0 ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 bool has_secondary_reserve( void ) const override {
  return( ( reserve_vars & 2u ) && ( f_gamma != 0 ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

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
 /// returns the reserve rows (2) of the class description

 const std::vector< FRowConstraint > &
 get_min_power_constraints( void ) const { return( min_power_Const ); }

/*--------------------------------------------------------------------------*/
 /// returns the reserve rows (1) of the class description
 
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
  * the specific extra solution information of IntermittentUnitBlock.
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
/** @name Methods for changing the data of the IntermittentUnitBlock
 *  @{ */

 /// set the maximum power values
 /** This function sets the maximum powers \f$ P^{mx}_t \f$ ("MaxPower") of
  * this IntermittentUnitBlock (e.g., a new forecast of the load factor), and
  * changes (1) and (3), or (1d) and (3d), of the class description
  * accordingly; a zero is replaced by the \f$ \varepsilon \f$ of
  * set_BlockConfig(), if positive.
  *
  * @param values  Iterator to a vector containing the maximum power values.
  * @param subset  If non-empty, the maximum power values corresponding to the
  *                indices in \p subset are set to the values pointed by
  *                \p values. If empty, no operation is performed.
  * @param ordered It indicates whether \p subset is ordered.
  * @param issuePMod Controls how physical Modifications are issued.
  * @param issueAMod Controls how abstract Modifications are issued.
  */

 void set_maximum_power( MF_dbl_it values ,
                         Subset && subset ,
                         const bool ordered = false ,
                         c_ModParam issuePMod = eNoBlck ,
                         c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// set the maximum power values
 /** This function sets the maximum powers \f$ P^{mx}_t \f$ ("MaxPower") of
  * this IntermittentUnitBlock (e.g., a new forecast of the load factor), and
  * changes (1) and (3), or (1d) and (3d), of the class description
  * accordingly; a zero is replaced by the \f$ \varepsilon \f$ of
  * set_BlockConfig(), if positive.
  *
  * @param values Iterator to a vector containing the maximum power values.
  * @param rng    If non-empty, the maximum power values corresponding to the
  *               indices in \p rng are set to the values pointed by
  *               \p values. If empty, no operation is performed.
  * @param issuePMod Controls how physical Modifications are issued.
  * @param issueAMod Controls how abstract Modifications are issued.
  */

 void set_maximum_power( MF_dbl_it values ,
                         Range rng = Range( 0 , Inf< Index >() ) ,
                         c_ModParam issuePMod = eNoBlck ,
                         c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// set the active power cost values
 /** This function sets the active power cost values of this
  * IntermittentUnitBlock.
  *
  * @param values  Iterator to a vector containing the active power cost
  *                values.
  * @param subset  If non-empty, the active power cost values corresponding
  *                to the indices in \p subset are set to the values pointed
  *                by \p values. If empty, no operation is performed.
  * @param ordered It indicates whether \p subset is ordered.
  * @param issuePMod Controls how physical Modifications are issued.
  * @param issueAMod Controls how abstract Modifications are issued.
  */
 void set_active_power_cost( MF_dbl_it values ,
                             Subset && subset ,
                             bool ordered = false ,
                             c_ModParam issuePMod = eNoBlck ,
                             c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// set the active power cost values
 /** This function sets the active power cost values of this
  * IntermittentUnitBlock.
  *
  * @param values Iterator to a vector containing the active power cost
  *               values.
  * @param rng    If non-empty, the active power cost values corresponding
  *               to the indices in \p rng are set to the values pointed by
  *               \p values. If empty, no operation is performed.
  * @param issuePMod Controls how physical Modifications are issued.
  * @param issueAMod Controls how abstract Modifications are issued.
  */
 void set_active_power_cost( MF_dbl_it values ,
                             Range rng = Range( 0 , Inf< Index >() ) ,
                             c_ModParam issuePMod = eNoBlck ,
                             c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// set the same active power cost at every instant
 /** This function sets the active power cost of this IntermittentUnitBlock
  * to \p value at every instant of the time horizon.
  *
  * @param value     The value of the active power cost.
  * @param issuePMod Controls how physical Modifications are issued.
  * @param issueAMod Controls how abstract Modifications are issued.
  */
 void set_active_power_cost( double value ,
                             c_ModParam issuePMod = eNoBlck ,
                             c_ModParam issueAMod = eNoBlck ) {
  std::vector< double > vector( get_time_horizon() , value );
  set_active_power_cost( vector.cbegin() ,
                         Range( 0 , Inf< Index >() ) ,
                         issuePMod , issueAMod );
 }

/*--------------------------------------------------------------------------*/
 /// sets the scale factor
 /** This method sets the scale factor of this IntermittentUnitBlock.
  *
  * @param values An iterator to a vector containing the scale factor.
  * @param subset If non-empty, the scale factor is set to the value pointed
  *               by \p values. If empty, no operation is performed.
  * @param ordered This parameter is ignored.
  * @param issuePMod Controls how physical Modifications are issued.
  * @param issueAMod Controls how abstract Modifications are issued. */

 void scale( MF_dbl_it values ,
             Subset && subset , bool ordered = false ,
             c_ModParam issuePMod = eNoBlck ,
             c_ModParam issueAMod = eNoBlck ) override;

/*--------------------------------------------------------------------------*/
 /// set the kappa constant
 /** This function sets \f$ \kappa \f$, which multiplies \f$ P^{mn}_t \f$
  * and \f$ P^{mx}_t \f$ in the rows (1)-(3) of the class description, or,
  * with the design variable, in its coefficients in (1d)-(3d).
  *
  * @param values  Iterator to a vector containing the kappa constants.
  * @param subset  If non-empty, the kappa constant corresponding to the
  *                indices in \p subset is set to the values pointed by
  *                \p values. If empty, no operation is performed.
  * @param ordered It indicates whether \p subset is ordered.
  * @param issuePMod Controls how physical Modifications are issued.
  * @param issueAMod Controls how abstract Modifications are issued. */

 void set_kappa( MF_dbl_it values ,
                 Subset && subset , bool ordered = false ,
                 c_ModParam issuePMod = eNoBlck ,
                 c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// set the kappa constant
 /** This function sets \f$ \kappa \f$, which multiplies \f$ P^{mn}_t \f$
  * and \f$ P^{mx}_t \f$ in the rows (1)-(3) of the class description, or,
  * with the design variable, in its coefficients in (1d)-(3d).
  *
  * @param values Iterator to a vector containing the kappa constants.
  * @param rng    If non-empty, the kappa constant corresponding to the
  *               indices in \p rng is set to the values pointed by
  *               \p values. If empty, no operation is performed.
  * @param issuePMod Controls how physical Modifications are issued.
  * @param issueAMod Controls how abstract Modifications are issued.
  */
 void set_kappa( MF_dbl_it values ,
                 Range rng = Range( 0 , Inf< Index >() ) ,
                 c_ModParam issuePMod = eNoBlck ,
                 c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// set the kappa constant
 /** This function sets \f$ \kappa \f$, which multiplies \f$ P^{mn}_t \f$
  * and \f$ P^{mx}_t \f$ in the rows (1)-(3) of the class description, or,
  * with the design variable, in its coefficients in (1d)-(3d).
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

/** @} ---------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED METHODS OF THE CLASS ---------------------*/
/*--------------------------------------------------------------------------*/

 /// refresh the f_scale-aware coefficients in the Objective
 /** Re-applies \f$ \sigma c \f$ to every coefficient in the Objective that
  * is meant to scale with the current scale factor \f$ \sigma \f$ (see
  * UnitBlock::scale()): the active-power costs and the investment cost on
  * the design variable. Called by scale() after f_scale changes so the
  * abstract representation stays in sync.
  *
  * @param issueAMod controls how abstract Modification are issued. */

 void update_objective( c_ModParam issueAMod );

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

 /// the maximum generation allowed
 double f_MaxGeneration;

 /// the minimum generation allowed
 double f_MinGeneration;

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

 /// the reserve rows (2), lower fence of the active power
 std::vector< FRowConstraint > min_power_Const;

 /// the reserve rows (1), upper fence of the active power
 std::vector< FRowConstraint > max_power_Const;

 /// the active power bounds design constraints
 /** The lower-bound rows are only there when the unit has a nonzero minimum
  * power, see generate_abstract_constraints(): with no minimum power the
  * design Variable has a zero coefficient in them and what is left is the
  * sign of the active power, which is a bound. Hence the first index is not
  * the side of the fence: it is design_max_row() for the upper-bound rows,
  * and 0 for the lower-bound ones when they are there. */
 boost::multi_array< FRowConstraint , 2 > active_power_bounds_design_Const;

 /// index of the upper-bound slice of active_power_bounds_design_Const
 [[nodiscard]] Index design_max_row( void ) const {
  return( active_power_bounds_design_Const.empty() ? 0 :
          Index( active_power_bounds_design_Const.shape()[ 0 ] ) - 1 );
  }

 /// true if active_power_bounds_design_Const carries the lower-bound rows
 [[nodiscard]] bool has_design_min_rows( void ) const {
  return( design_max_row() > 0 );
  }

 /// the active power bounds constraints
 std::vector< BoxConstraint > active_power_bounds_Const;

 /// the reactive power bound constraints
 std::vector< BoxConstraint > ReactivePower_Bound_Const;

 /// the active power bounds constraints
 FRowConstraint MaxMinGeneration_Const;

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
 /// updates the reserve row (1), if max_side, or (2) of instant t
 /** Writes the current \f$ \gamma \kappa P^{mx}_t \f$, or
  * \f$ \kappa P^{mn}_t \f$, in the row: as a side of it with no design
  * Variable, as the coefficient (with the opposite sign) of the design
  * Variable otherwise [see generate_abstract_constraints()]. */

 void update_reserve_row( Index t , bool max_side , c_ModParam issueAMod );

/*--------------------------------------------------------------------------*/
 /// verify whether the data in this IntermittentUnitBlock is consistent
 /** Checks the data of this IntermittentUnitBlock, with the symbols of the
  * class description, and throws std::logic_error if one of the following
  * conditions does not hold:
  *
  * - \f$ P^{mn}_t \leq P^{mx}_t \f$ for all \f$ t \f$ (\f$ P^{mn}_t \f$ of
  *   any sign);
  *
  * - \f$ G^{mn} \leq G^{mx} \f$;
  *
  * - \f$ 0 \leq \gamma \leq 1 \f$ and \f$ \kappa \geq 0 \f$ (values of
  *   \f$ \kappa \f$ above 1 are accepted);
  *
  * - \f$ h^p_t \geq 0 \f$ for all \f$ t \f$;
  *
  * - \f$ \underline{X} \geq 0 \f$, \f$ \underline{X} \leq \bar{X} \f$ if
  *   \f$ \bar{X} > 0 \f$, \f$ \underline{X} \leq 1 \f$ if
  *   \f$ | \bar{X} | = 1 \f$, \f$ \underline{X} \leq | \bar{X} | \f$ if
  *   \f$ \bar{X} < 0 \f$, and \f$ | \bar{X} | \leq 1 \f$ if
  *   \f$ \sigma \neq 1 \f$. */

 void check_data_consistency( void ) const;

/*--------------------------------------------------------------------------*/

 static void static_initialization( void )
 {
  register_method< IntermittentUnitBlock , MF_dbl_it , Subset && , bool >(
   "IntermittentUnitBlock::set_maximum_power" ,
   & IntermittentUnitBlock::set_maximum_power );

  register_method< IntermittentUnitBlock , MF_dbl_it , Range >(
   "IntermittentUnitBlock::set_maximum_power" ,
   & IntermittentUnitBlock::set_maximum_power );

  register_method< IntermittentUnitBlock , MF_dbl_it , Subset && , bool >(
   "IntermittentUnitBlock::set_active_power_cost" ,
   & IntermittentUnitBlock::set_active_power_cost );

  register_method< IntermittentUnitBlock , MF_dbl_it , Range >(
   "IntermittentUnitBlock::set_active_power_cost" ,
   & IntermittentUnitBlock::set_active_power_cost );

  register_method< IntermittentUnitBlock , MF_dbl_it , Subset && , bool >(
   "IntermittentUnitBlock::scale" ,
   & IntermittentUnitBlock::scale );

  register_method< IntermittentUnitBlock , MF_dbl_it , Range >(
   "IntermittentUnitBlock::scale" ,
   & IntermittentUnitBlock::scale );

  register_method< IntermittentUnitBlock , MF_dbl_it , Subset && , bool >(
   "IntermittentUnitBlock::set_kappa" ,
   & IntermittentUnitBlock::set_kappa );

  register_method< IntermittentUnitBlock , MF_dbl_it , Range >(
   "IntermittentUnitBlock::set_kappa" ,
   & IntermittentUnitBlock::set_kappa );

  // The same two methods are registered again under names that say which of
  // the two ways of sizing a Block [see Design and scaling of this Block in
  // Block.h] they implement, so that a consumer can pick the operation
  // without knowing the class:
  //
  // - replicate: this UnitBlock stands for k identical copies of itself. Its
  //   Constraint and Variable keep describing one copy and whoever holds a
  //   coupling multiplies by get_scale(); only its Objective is rewritten;
  //
  // - resize: this UnitBlock stands for one unit of k times the size it was
  //   given. Its own rows are rewritten and no consumer multiplies anything.
  //
  // The old names are kept so that instances naming them keep working: a
  // registered name travels inside the instances, and dropping one breaks
  // them at load time with no diagnostic.

  register_method< IntermittentUnitBlock , MF_dbl_it , Subset && , bool >(
   "IntermittentUnitBlock::replicate" ,
   & IntermittentUnitBlock::scale );

  register_method< IntermittentUnitBlock , MF_dbl_it , Range >(
   "IntermittentUnitBlock::replicate" ,
   & IntermittentUnitBlock::scale );

  register_method< IntermittentUnitBlock , MF_dbl_it , Subset && , bool >(
   "IntermittentUnitBlock::resize" ,
   & IntermittentUnitBlock::set_kappa );

  register_method< IntermittentUnitBlock , MF_dbl_it , Range >(
   "IntermittentUnitBlock::resize" ,
   & IntermittentUnitBlock::set_kappa );

  // ... and the getter reading the sensitivity back, registered under a name
  // that says which of the two it answers for. It is the same virtual that
  // UnitBlock declares; what the factory adds is that a consumer can reach it
  // without holding a UnitBlock *, and that a class not answering for it
  // simply does not register the name.

  using qry_sbst = QueryType< MF_dbl_msp , c_Subset & , bool >;
  using qry_rngd = QueryType< MF_dbl_msp , Range >;

  register_method< qry_sbst >(
   "IntermittentUnitBlock::get_resize_linearization" , new qry_sbst(
    []( const Block * blck , MF_dbl_msp msp , c_Subset & idxs , bool ) {
     // one size parameter, so the only index that exists here is 0
     if( idxs.size() > msp.size() )
      throw( std::invalid_argument(
       "IntermittentUnitBlock::get_resize_linearization: the span is shorter "
       "than the subset it is asked to answer for" ) );
     for( Index i = 0 ; i < idxs.size() ; ++i ) {
      if( idxs[ i ] )
       throw( std::invalid_argument(
        "IntermittentUnitBlock::get_resize_linearization: a unit has one size "
        "parameter, and there is no index " + std::to_string( idxs[ i ] ) ) );
      msp[ i ] = static_cast< const IntermittentUnitBlock * >( blck )->
                 get_kappa_linearization();
      }
     } ) );

  register_method< qry_rngd >(
   "IntermittentUnitBlock::get_resize_linearization" , new qry_rngd(
    []( const Block * blck , MF_dbl_msp msp , Range rng ) {
     // one size parameter: the family has one element, of index 0
     rng.second = std::min( rng.second , Index( 1 ) );
     if( rng.first >= rng.second )
      return;
     if( msp.empty() )
      throw( std::invalid_argument(
       "IntermittentUnitBlock::get_resize_linearization: the span is shorter "
       "than the range it is asked to answer for" ) );
     msp[ 0 ] = static_cast< const IntermittentUnitBlock * >( blck )->
                get_kappa_linearization();
     } ) );
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
  eSetActPCost ,               ///< set active power cost values
  eIUBModLastParam             ///< first allowed parameter for derived classes
  /**< Convenience value to easily allow derived classes to extend the set
   * of types of IntermittentUnitBlockMod. */
  };

 /// constructor, takes the IntermittentUnitBlock and the type
 IntermittentUnitBlockMod( IntermittentUnitBlock * const fblock ,
                           const int type )
  : UnitBlockMod( fblock , type ) , f_Block( fblock ) {}

 /// destructor, does nothing
 virtual ~IntermittentUnitBlockMod() override = default;

 /// returns the Block to which the Modification refers
 Block * get_Block( void ) const override { return( f_Block ); }

 protected:

 /// prints the IntermittentUnitBlockMod
 void print( std::ostream & output ) const override {
  output << "IntermittentUnitBlockMod[" << this << "]: ";
  switch( f_type ) {
   case( eSetMaxP ):
    output << "Set max power values ";
    break;
   case( eSetKappa ):
    output << "Set the kappa constant ";
    break;
   case( eSetActPCost ):
    output << "Set active power cost values ";
    break;
   default:;
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
