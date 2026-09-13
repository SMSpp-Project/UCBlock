/*--------------------------------------------------------------------------*/
/*------------------------- File NuclearUnitBlock.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class NuclearUnitBlock, which derives from
 * ThermalUnitBlock [see ThermalUnitBlock.h] and therefore from UnitBlock
 * [see UnitBlock.h], in order to define a "nuclear unit", i.e., "a thermal
 * unit subject to modulation constraints".
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Antonio Frangioni, Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __NuclearUnitBlock
 #define __NuclearUnitBlock
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "ThermalUnitBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------------ NAMESPACE ---------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{
/*--------------------------------------------------------------------------*/
/*----------------------- CLASS NuclearUnitBlock ---------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// implementation of the Block concept for a nuclear (thermal) unit
/** The NuclearUnitBlock class derives from ThermalUnitBlock (and therefore
 * from UnitBlock) and implements the concept of a "nuclear unit", i.e.,
 * a thermal unit (as implemented in ThermalUnitBlock) subject to extra
 * "modulation constraints" that limit the number of time instants in which
 * the energy production can (significantly) change.
 *
 * The class makes some assumptions about the behaviour and implementation of
 * the base class:
 *
 * - whatever formulation is implemented in ThermalUnitBlock, it comprises
 *   the start-up and shut-down variables (v_t and w_t in the standard
 *   notation), where
 *
 *    v_t = v_start_up[ t - init_t ]  ,  w_t = v_shut_down[ t - init_t ]
 *
 *   with init_t a field of the class that indicates the first time instant
 *   in which the commitment variables are free to be chosen (while in all
 *   instants 0 <= t < init_t, if any, they are fixed to either 0 or 1
 *   depending on the state of the unit prior to the beginning of time
 *   (instant 0) due to the minimum up- or down-time constraints
 *
 * - variables_generated() and constraints_generated() can be used to assess
 *   whether or not (static) Variable and Constraint have already been
 *   generated
 *
 * - standard ramp-up and ramp-down deltas are defined, i.e.,
 *   v_DeltaRampUp and v_DeltaRampDown are nonempty
 *
 * The operating rules of the unit are those of a nuclear unit that follows
 * the load: the output of an on unit that is not
 * starting up is either *stable*, in which case it moves by at most the
 * (small) stability ramps ModulationDeltaRampUp / ModulationDeltaRampDown,
 * or it performs a *modulation step*, upwards or downwards. A *modulation*
 * is a maximal sequence of consecutive modulation steps in the same
 * direction, lasting at most MaxModulationLength instants: all its steps
 * but the last move the output by exactly the full ramp (DeltaRampUp or
 * DeltaRampDown), the last one by at most as much (in the same direction).
 * After the end of a modulation no other one can start for ModulationTime
 * - 1 instants, and no modulation can start in the StabilityAfterStartUp
 * instants that follow a start-up; with PowerBands the output is split into
 * three bands, a stable instant stays in its own and a modulation moves to
 * an adjacent one; at most ModulationsPerDay modulations can start in each day
 * of DayLength instants, at most StartUpsPerDay start-ups can happen in
 * each day, and at most DeepDecreasesPerDay "deep decreases" (a decrease by
 * at least DeepDecreaseGradient to an output not larger than
 * DeepDecreaseThreshold) can happen in each day. With the defaults
 * (MaxModulationLength = 1 and no daily limits) every modulation is a
 * single instant in which the full ramp is allowed, two modulations are at
 * least ModulationTime instants apart, and the model is the original one of
 * the class. See generate_abstract_constraints() for the formulation.
 *
 * Current significant limitations of the class are:
 *
 * - the duration of the modulation period and the initial state of
 *   modulation cannot be changed at all, nor can the data of the operating
 *   rules other than the stability ramps
 *
 * - the unit is assumed not to be in the middle of a modulation at the
 *   beginning of the horizon, and to have performed no modulation, start-up
 *   and deep decrease yet in the first day
 *
 * - the costs of the downward modulation steps and of the deep decreases
 *   (DownModulationCost and DeepDecreaseCost) cannot be changed, neither
 *   from the physical nor from the abstract representation
 *
 * - the modulation ramp-up and ramp-down deltas can be changed, but only
 *   from the physical representation: doing that from the abstract
 *   representation is not allowed (it requires changing coefficients in
 *   the corresponding constraints, which is not supported already in
 *   ThermalUnitBlock, so one would get an exception)
 */

class NuclearUnitBlock : public ThermalUnitBlock
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*----------------------------- CONSTANTS ----------------------------------*/

 /// mask for the 6th bit of the formulation code, == 1 if the tight rows
 /** The int Configuration that selects the formulation of the
  * ThermalUnitBlock [see ThermalUnitBlock::generate_abstract_variables()]
  * also selects, with this bit, the formulation of the operating rules: the
  * default (bit off) is the one described in generate_operating_rules(),
  * while with the bit on the stability, the longest modulation and the deep
  * decreases are described by the tight rows discussed there, which have a
  * stronger continuous relaxation at the price of the (continuous) start
  * Variable and of \f$ O( n \tau^M ) \f$ rows. Both describe the same set of
  * schedules, hence the specialised Solver is indifferent to the choice. */

 static constexpr unsigned char TightRules = 32;

 /// mask for the 7th bit of the formulation code, == 1 if the tight big-M
 /** With this bit the rows of the full ramp of a modulation use one
  * coefficient per case rather than one for all of them: see
  * generate_operating_rules(). It is independent from TightRules, so that
  * the two families can be used, and measured, separately. */

 static constexpr unsigned char TightRamp = 64;

 /// mask for the 8th bit of the formulation code, == 1 if the tight rows of
 /// TightRules are separated rather than written
 /** The tight rows of TightRules are \f$ O( n \tau^M ) \f$, which at a long
  * horizon is a model several times larger than the default one; since they
  * are valid inequalities, and not part of the description of the schedules,
  * they can be left out and added only where they are violated, exactly as
  * the Perspective Cuts of the ThermalUnitBlock are [see
  * generate_dynamic_constraints()]. With this bit on (and TightRules on) the
  * three families written on the starts are separated instead of being
  * written; the rows of the deep decreases, which are \f$ O( n ) \f$, stay
  * where they are. */

 static constexpr unsigned char TightCuts = 128;

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

 /// constructor, takes the father and the time horizon
 /** Constructor of NuclearUnitBlock, possibly taking a pointer of its
  * father Block (presumably, but not necessarily, a UCBlock). */

 explicit NuclearUnitBlock( Block * f_block = nullptr ) :
  ThermalUnitBlock( f_block ) , f_initial_modulation( 0 ) ,
  f_modulation_interval( 0 ) {}

/*--------------------------------------------------------------------------*/
 /// destructor of NuclearUnitBlock

 virtual ~NuclearUnitBlock() override;

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

/// Extends Block::deserialize( netCDF::NcGroup )
/** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
 * the NuclearUnitBlock, i.e., that of ThermalUnitBlock plus the extra
 * information needed to handle modulation constraints:
 *
 * - The scalar variable "ModulationTime" of type netCDF::NcUint for the
 *   modulation interval, i.e., the number of consecutive time instants
 *   within which at most one modulation is allowed. The value must be >= 2,
 *   as otherwise the modulation constraint is useless and one can use an
 *   original ThermalUnitBlock. Note that, like with "MinUpTime" and
 *   "MinDownTime", there is no way to change the value of "ModulationTime"
 *   once the NuclearUnitBlock is loaded, so useless once, useless always.
 *   This variable is optional. If it is not provided, then the value is
 *   assumed to be 2.
 *
 * - The scalar variable "InitModulation" of type netCDF::NcUint for the
 *   initial modulation state, i.e., the number of instants before the first
 *   (0) when the unit last modulated. The value must be >= 1, as a value of
 *   0 would mean that the unit modulated at the first interval, which is
 *   nonsensical since it may be forced to be off, and even if it is on is
 *   should be free to decide whether or not modulating then: "the history
 *   is written only for instants in the past, i.e, before 0". This variable
 *   is optional. If it is not provided, then  the value is assumed to be
 *   equal to that of "ModulationTime", which means that lhe unit last
 *   modulated "a long time ago" (possibly never) and it is allowed to start
 *   modulating immediately (provided it is on, see "InitUpDownTime" in
 *   the original ThermalUnitBlock).
 *
 * - The variable "ModulationDeltaRampUp" of type double and either of size 1
 *   or indexed over the dimension "NumberIntervals"; if "NumberIntervals" is
 *   not provided, then this variable can also be indexed over "TimeHorizon".
 *   This is meant to represent the vector MDP[ t ] that, for each time
 *   instant t, contains the modulation ramp-up value of the unit for the
 *   corresponding time step, i.e., the maximum possible increase of active
 *   power production w.r.t. the power that had been produced in time instant
 *   t - 1 *unless a modulation is performed* (note that start-ups are not
 *   modulations, hence for a modulation to be performed the unit necessarily
 *   had to be producing power at t - 1). This variable is optional; if it is
 *   not provided then it is assumed that MDP[ t ] == 0, i.e., the unit cannot
 *   increase its power output unless a modulation is performed. Note that the
 *   value of the variable must always be such that 0 <= MDP[ t ] <= DP[ t ],
 *   where DP[ t ] is the vector of original ramp-up values fot the unit (see
 *   "DeltaRampUp" in the original ThermalUnitBlock). If
 *   "ModulationDeltaRampUp" has length 1, then MDP[ t ] contains the same
 *   value for all t. Otherwise, ModulationDeltaRampUp[ i ] is the fixed
 *   value of MDP[ t ] for all t in the interval [ ChangeIntervals[ i - 1 ] ,
 *   ChangeIntervals[ i ] ], with the assumption that ChangeIntervals[ - 1 ]
 *    = 0. If NumberIntervals <= 1 or NumberIntervals >= TimeHorizon, then
 *   the mapping clearly does not require "ChangeIntervals", which in fact
 *   is not loaded.
 *
 * - The variable "ModulationDeltaRampDown" of type double and either of size
 *   1 or indexed over the dimension "NumberIntervals"; if "NumberIntervals"
 *   is not provided, then this variable can also be indexed over
 *   "TimeHorizon". This is meant to represent the vector MDM[ t ] that, for
 *   each time instant t, contains the modulation ramp-down value of the unit
 *   for the corresponding time step, i.e., the maximum possible decrease of
 *   active power production w.r.t. the power that had been produced in time
 *   instant t - 1 *unless a modulation is performed* (note that start-ups
 *   are not modulations, hence for a modulation to be performed the unit
 *   necessarily had to be producing power at t - 1). This variable is
 *   optional; if it is not provided then it is assumed that MDM[ t ] == 0,
 *   i.e., the unit cannot decrease its power output unless a modulation is
 *   performed. Note that the value of the variable must always be such that
 *   0 <= MDM[ t ] <= DM[ t ], where DM[ t ] is the vector of original
 *   ramp-down values fot the unit (see "DeltaRampDown" in the original
 *   ThermalUnitBlock). If "ModulationDeltaRampDown" has length 1, then
 *   MDM[ t ] contains the same value for all t. Otherwise,
 *   ModulationDeltaRampDown[ i ] is the fixed value of MDM[ t ] for all t in
 *   the interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with
 *   the assumption that ChangeIntervals[ - 1 ]  = 0. If NumberIntervals <= 1
 *   or NumberIntervals >= TimeHorizon, then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded.
 *
 * - The scalar variable "MaxModulationLength" of type netCDF::NcUint for the
 *   maximum number \f$ L^M \geq 1 \f$ of consecutive instants a modulation
 *   may last. This variable is optional; if it is not provided, then
 *   \f$ L^M = 1 \f$, i.e., every modulation is a single instant. For
 *   \f$ L^M > 1 \f$, "ModulationTime" counts from the end of a modulation:
 *   if a modulation ends at \f$ t \f$, then no modulation step can happen
 *   at \f$ t + 1 , \ldots , t + \tau^M - 1 \f$, and "InitModulation" is the
 *   number of instants before 0 at which the last modulation ended.
 *
 * - The variable "PowerBands" of type netCDF::NcDouble and of size 2, the
 *   two breakpoints \f$ \underline{p} < b_1 < b_2 < \overline{p} \f$ that
 *   split the output of the unit into a low band \f$ [ \underline{p} ,
 *   b_1 ] \f$, an intermediate one \f$ [ b_1 , b_2 ] \f$ and a high one
 *   \f$ [ b_2 , \overline{p} ] \f$. This variable is optional; if it is not
 *   provided the output is not banded, which is the default. With the bands,
 *   a stable instant keeps the output in the band it is in and a modulation
 *   moves it to an adjacent band, where it has to land at its last step,
 *   while the instants in between are free [see
 *   generate_abstract_constraints()].
 *
 * - The scalar variable "StabilityAfterStartUp" of type netCDF::NcUint for
 *   the number \f$ A \geq 0 \f$ of instants, the start-up one comprised, in
 *   which a unit that has just started up cannot begin a modulation. This
 *   variable is optional; if it is not provided, then \f$ A = 0 \f$, i.e.,
 *   a modulation may start at the instant that follows a start-up (at the
 *   start-up one it cannot, a start-up not being a modulation).
 *
 * - The scalar variables "ModulationsPerDay", "StartUpsPerDay" and
 *   "DeepDecreasesPerDay" of type netCDF::NcUint for the maximum number of
 *   modulations (counted at their first step), of start-ups and of deep
 *   decreases in each day. Each of them is optional; if it is not provided,
 *   the corresponding number is not limited.
 *
 * - The scalar variable "DayLength" of type netCDF::NcUint for the number
 *   of instants of a day; the days are the consecutive intervals of
 *   DayLength instants starting at 0. This variable is optional; if it is
 *   not provided (or it is 0), then the whole horizon is a single day.
 *
 * - The variable "DownModulationCost" of type double, either of size 1 or
 *   indexed over "NumberIntervals" (or "TimeHorizon" if "NumberIntervals" is
 *   not provided), with the same mapping as "ModulationDeltaRampUp", for the
 *   cost \f$ c^-_t \geq 0 \f$ of each downward modulation step. This
 *   variable is optional; if it is not provided, the cost is 0.
 *
 * - The variables "DeepDecreaseThreshold", "DeepDecreaseGradient" and
 *   "DeepDecreaseCost", with the same shape, for the threshold
 *   \f$ \tilde{p}_t \f$, the gradient \f$ \tilde{\Delta}_t > 0 \f$ and the
 *   cost \f$ c^d_t \geq 0 \f$ of a deep decrease, which happens at an on
 *   instant \f$ t \f$ with \f$ u_{t-1} = 1 \f$ when \f$ p_{t-1} - p_t
 *   \geq \tilde{\Delta}_t \f$ and \f$ p_t \leq \tilde{p}_t \f$. The deep
 *   decreases exist if and only if both "DeepDecreaseThreshold" and
 *   "DeepDecreaseGradient" are provided; "DeepDecreaseCost" is optional, the
 *   cost being 0 if it is not provided.
 *
 * - The variables "ModulationDeltaRampUp" and "ModulationDeltaRampDown" are
 *   optional as described above, being 0 if they are not provided.
 *
 * Important: nuclear units *must* have ramp-up and ramp-down constraints,
 * hence DeltaRampUp and DeltaRampDown must be provided in the \p group
 * (although they are in principle optional for ThermalUnitBlock). */

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG
 /* extends ThermalUnitBlock::expected_dims()
  * not necessary, no new dimensions

 std::vector< std::string > expected_dims( void ) const override;
 */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// extends ThermalUnitBlock::expected_vars()

 std::vector< std::string > expected_vars( void ) const override;

#endif

/*--------------------------------------------------------------------------*/
/// generate the abstract variables of the NuclearUnitBlock
/** Besides those of ThermalUnitBlock (depending on the exact chosen
 * formulation, see the Configuration parameter in the ThermalUnitBlock
 * class), NuclearUnitBlock has the std::vector of size TimeHorizon
 * containing the binary variables m[ t ] indicating whether or not the unit
 * performs a modulation step at time instant t. These variables are
 * mandatory. The following ones, all of size TimeHorizon, only exist when
 * the operating rules need them:
 *
 * - the binary variables d[ t ] indicating whether the modulation step at t
 *   is downwards (m[ t ] - d[ t ] being the upward one), which exist if the
 *   direction of the modulation matters, i.e., if MaxModulationLength > 1,
 *   or DownModulationCost is nonzero, or the deep decreases exist
 *   [see has_modulation_direction()];
 *
 * - the variables s[ t ] in [ 0 , 1 ] that are at least 1 when a modulation
 *   starts at t, which exist if MaxModulationLength > 1 and the number of
 *   modulations per day is limited;
 *
 * - the binary variables \f$ \delta_t \f$, \f$ \delta'_t \f$ and
 *   \f$ \delta''_t \f$ indicating a deep decrease at t, a decrease by more
 *   than the gradient and an output below the threshold, which exist if the
 *   deep decreases exist.
 *
 * The Configuration parameter has no use here and it is just passed up to
 * the method of the base ThermalUnitBlock class that is called first thing
 * here inside. */

 void generate_abstract_variables( Configuration *stvv = nullptr ) override;

/*--------------------------------------------------------------------------*/
/// Generate the static constraint of the NuclearUnitBlock
/** Besides those of ThermalUnitBlock (depending on the exact chosen
 * formulation, see the Configuration parameter in the ThermalUnitBlock
 * class), this method generates the abstract constraints of the
 * NuclearUnitBlock.
 *
 * The method constructs five groups of constraints, oll of which are
 * mandatory and instantiated for all time instants t:
 *
 * - the modulation ramp-up constraint
 *   \f[
 *     p_t - p_{t-1} \leq \Delta^M_{t+} u_{t-1} +
 *     ( \Delta_{t+} - \Delta^M_{t+} ) m_t + \bar{l}_t v_t
 *   \f]
 *
 *   with \f$ \Delta^M_{t+} \f$ the modulation ramp-up value at time t,
 *   \f$ \Delta_{t+} \f$ the original ramp-up value at time t, and
 *   \f$ \bar{l}_t \f$ the start-up limit at time t
 *
 * - the modulation ramp-down constraint
 *   \f[
 *     p_{t-1} - p_t \leq \Delta^M_{t-} u_t +
 *     ( \Delta_{t-} - \Delta^M_{t-} ) m_t + \bar{u}_t w_t
 *   \f]
 *
 *   with \f$ \Delta^M_{t-} \f$ the modulation ramp-down value at time t,
 *   \f$ \Delta_{t-} \f$ the original ramp-down value at time t, and
 *   \f$ \bar{u}_t \f$ the shut-down limit at time t
 *
 * - the logical constraints \f$ m_t \leq u_t \f$ (modulations cannot
 *   happen when the unit is down)
 *
 * - the logical constraints \f$ m_t \leq ( 1 - v_t ) \f$ (start-ups are
 *   not modulations, so the unit cannot modulate while starting up=
 *
 * - the modulation constraint proper
 *   \f[
 *     \sum_{h = \max\{ 0 , t - \tau^M + 1 \}}^t m_h \leq 1
 *   \f]
 *
 *   where \f$ \tau^M \f$ is the modulation interval.
 *
 * Note that another group of constraints is needed to fix to 0 certain
 * modulation variables depending on the initial modulation value; this is
 * done by changing the bounds (and fixing them) and therefore it does not
 * result in a separate group of constraints.
 *
 * When the direction of the modulation matters the two ramp constraints
 * become
 *   \f[
 *     p_t - p_{t-1} \leq \Delta^M_{t+} ( u_{t-1} - m_t ) +
 *     \Delta_{t+} ( m_t - d_t ) + \bar{l}_t v_t
 *     \quad , \quad
 *     p_{t-1} - p_t \leq \Delta^M_{t-} ( u_t - m_t ) +
 *     \Delta_{t-} d_t + \bar{u}_t w_t
 *   \f]
 * (i.e., an upward step cannot decrease the output and a downward one
 * cannot increase it), together with \f$ d_t \leq m_t \f$. When
 * \f$ L^M > 1 \f$ the modulation window is replaced by the constraints
 *   \f[
 *   \begin{array}{ll}
 *     d_{t+1} - d_t \leq 1 - m_t \, , \quad d_t - d_{t+1} \leq 1 - m_{t+1}
 *     & \mbox{(same direction along a modulation)} \\
 *     p_t - p_{t-1} \geq \Delta_{t+} m_{t+1} - M_t ( 1 - m_t + d_t ) &
 *     \mbox{(full ramp up but at the last step)} \\
 *     p_{t-1} - p_t \geq \Delta_{t-} m_{t+1} - M'_t ( 1 - d_t ) &
 *     \mbox{(full ramp down but at the last step)} \\
 *     \sum_{h = 2}^{K_t} m_{t+h} + ( K_t - 1 ) ( m_t - m_{t+1} ) \leq
 *     K_t - 1 & \mbox{(stability)} \\
 *     \sum_{h = t}^{t + L^M} m_h \leq L^M & \mbox{(maximum length)}
 *   \end{array}
 *   \f]
 * with \f$ M_t = \Delta_{t+} + \max\{ \Delta_{t-} , \bar{u}_t \} \f$ and
 * \f$ M'_t = \Delta_{t-} + \max\{ \Delta_{t+} , \bar{l}_t \} \f$ (the
 * modulation ramp constraints bound the decrease of the output by
 * \f$ \max\{ \Delta_{t-} , \bar{u}_t \} \f$ and its increase by
 * \f$ \max\{ \Delta_{t+} , \bar{l}_t \} \f$) and
 * \f$ K_t = \min\{ \tau^M - 1 , T - 1 - t \} \f$ (a modulation ending at t
 * forbids the steps up to \f$ t + \tau^M - 1 \f$). The daily limits, for
 * each day \f$ \mathcal{D} \f$, are
 *   \f[
 *     \sum_{t \in \mathcal{D}} m_t \leq C \mbox{ if } L^M = 1 \, , \quad
 *     s_t \geq m_t - m_{t-1} \, , \, \sum_{t \in \mathcal{D}} s_t \leq C
 *     \mbox{ otherwise} \, , \quad
 *     \sum_{t \in \mathcal{D}} v_t \leq V \, , \quad
 *     \sum_{t \in \mathcal{D}} \delta_t \leq A
 *   \f]
 * and the deep decreases are defined by
 *   \f[
 *     p_t \geq \tilde{p}_t ( 1 - \delta''_t ) \, , \quad
 *     p_{t-1} - p_t \leq \tilde{\Delta}_t + ( \max\{ \Delta_{t-} ,
 *     \bar{u}_t \} - \tilde{\Delta}_t )^+ \delta'_t \, , \quad
 *     \delta_t \geq \delta'_t + \delta''_t + u_t - 2
 *   \f]
 * for the instants t with an on predecessor (\f$ t \geq 1 \f$, or t = 0 if
 * the unit is on at the beginning); a decrease exactly equal to the
 * gradient, or to an output exactly equal to the threshold, is therefore
 * not a deep decrease, which is the closure the dynamic programming Solver
 * uses too.
 *
 * The Configuration parameter has no use here and it is just passed up to
 * the method of the base ThermalUnitBlock class that is called first thing
 * here inside. */

 void generate_abstract_constraints( Configuration * stcc = nullptr ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// generate the Objective of the NuclearUnitBlock
/** The Objective is that of ThermalUnitBlock, plus the costs
 * \f$ c^-_t d_t \f$ of the downward modulation steps and
 * \f$ c^d_t \delta_t \f$ of the deep decreases, when they are nonzero; the
 * corresponding Variable are appended after all those of the
 * ThermalUnitBlock [see ThermalUnitBlock::objective_tail()], and their
 * coefficients cannot be changed via the abstract representation. */

 void generate_objective( Configuration * objc = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
/*--------------- Methods for checking the NuclearUnitBlock ----------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for checking solution information in the NuclearUnitBlock
 *  @{ */

 /// returns true if the current solution is (approximately) feasible
 /** This function returns true if and only if the solution encoded in the
  * current value of the Variable of this NuclearUnitBlock is approximately
  * feasible within the given tolerance. That is, a solution is considered
  * feasible if and only if
  *
  *   -# is is (approximately) feasible for the original ThermalUnitBlock;
  *
  *   -# each new ColVariable of NuclearUnitBlock is feasible;
  *
  *   -# the violation of each new Constraint of this NuclearUnitBlock is not
  *      greater than the tolerance.
  *
  * Every Constraint of this NuclearUnitBlock is a RowConstraint and its
  * violation is given by either the relative (see RowConstraint::rel_viol())
  * or the absolute violation (see RowConstraint::abs_viol()), depending on
  * the Configuration that is provided.
  *
  * The tolerance and the type of violation can be provided by either \p fsbc
  * or #f_BlockConfig->f_is_feasible_Configuration and they are determined as
  * follows:
  *
  *   - If \p fsbc is not a nullptr and it is a pointer to a
  *     SimpleConfiguration< double >, then the tolerance is the value present
  *     in that SimpleConfiguration and the relative violation is considered.
  *
  *   - If \p fsbc is not nullptr and it is a
  *     SimpleConfiguration< std::pair< double , int > >, then the tolerance is
  *     fsbc->f_value.first and the type of violation is determined by
  *     fsbc->f_value.second (any nonzero number for relative violation and
  *     zero for absolute violation);
  *
  *   - Otherwise, if both #f_BlockConfig and
  *     f_BlockConfig->f_is_feasible_Configuration are not nullptr and the
  *     latter is a pointer to either a SimpleConfiguration< double > or to a
  *     SimpleConfiguration< std::pair< double , int > >, then the values of the
  *     parameters are obtained analogously as above;
  *
  *   - Otherwise, by default, the tolerance is 0 and the relative violation
  *     is considered.
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
  *        and the type of violation that must be considered. */

 bool is_feasible( bool useabstract = false ,
                   Configuration * fsbc = nullptr ) override;


/** @} ---------------------------------------------------------------------*/
/*--------- METHODS FOR READING THE DATA OF THE NuclearUnitBlock -----------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the NuclearUnitBlock
 * @{ */

 /// returns the initial modulation value
 double get_initial_modulation( void ) const {
  return( f_initial_modulation );
  }

 /// returns the minimum allowed modulation interval
 Index get_modulation_interval( void ) const {
  return( f_modulation_interval );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of modulation delta ramp-up
 /** The returned vector contains the modulation delta ramp-up at each time.
  * The size of the vector is always get_time_horizon(). */

 const std::vector< double > & get_modulation_ramp_up( void ) const {
  return( v_modulation_ramp_up );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of modulation delta ramp-down
 /** The returned vector contains the modulation delta ramp-down at each
  * time. The size of the vector is always get_time_horizon(). */

 const std::vector< double > & get_modulation_ramp_down( void ) const {
  return( v_modulation_ramp_down );
  }

/*--------------------------------------------------------------------------*/
 /// returns the two breakpoints of the bands of the output, if any
 /** Returns the two breakpoints that split the output of the unit into
  * three bands, an empty vector if the output is not banded. */

 const std::vector< double > & get_power_bands( void ) const {
  return( v_power_bands );
  }

/*--------------------------------------------------------------------------*/
 /// true if the output of the unit is split into bands
 bool has_power_bands( void ) const { return( ! v_power_bands.empty() ); }

/*--------------------------------------------------------------------------*/
 /// returns the instants of stability that follow a start-up ( >= 0 )
 /** Returns the number of instants, the start-up one comprised, in which a
  * unit that has just started up cannot begin a modulation; 0 (the default)
  * leaves only the start-up instant, in which a modulation is impossible
  * anyway [see generate_operating_rules()]. */

 Index get_stability_after_start_up( void ) const {
  return( f_stability_after_start );
  }

/*--------------------------------------------------------------------------*/
 /// returns the maximum number of instants of a modulation ( >= 1 )
 Index get_max_modulation_length( void ) const {
  return( f_max_modulation_length );
  }

 /// returns the maximum number of modulations per day, -1 if unlimited
 int get_modulations_per_day( void ) const {
  return( f_modulations_per_day );
  }

 /// returns the maximum number of start-ups per day, -1 if unlimited
 int get_start_ups_per_day( void ) const { return( f_start_ups_per_day ); }

 /// returns the maximum number of deep decreases per day, -1 if unlimited
 int get_deep_decreases_per_day( void ) const {
  return( f_deep_decreases_per_day );
  }

 /// returns the number of instants of a day, 0 if the horizon is one day
 Index get_day_length( void ) const { return( f_day_length ); }

 /// returns the day of the time instant t
 Index get_day( Index t ) const {
  return( f_day_length ? t / f_day_length : 0 );
  }

/*--------------------------------------------------------------------------*/
 /// returns the cost of the downward modulation steps
 /** The returned vector is either empty, meaning that the cost is 0, or of
  * size get_time_horizon(). */

 const std::vector< double > & get_down_modulation_cost( void ) const {
  return( v_down_modulation_cost );
  }

 /// returns the deep-decrease threshold
 /** The returned vector is either empty, meaning that there are no deep
  * decreases, or of size get_time_horizon(), and so is that of
  * get_deep_decrease_gradient(). */

 const std::vector< double > & get_deep_decrease_threshold( void ) const {
  return( v_deep_threshold );
  }

 /// returns the deep-decrease gradient [see get_deep_decrease_threshold()]
 const std::vector< double > & get_deep_decrease_gradient( void ) const {
  return( v_deep_gradient );
  }

 /// returns the cost of a deep decrease, empty if it is 0
 const std::vector< double > & get_deep_decrease_cost( void ) const {
  return( v_deep_cost );
  }

 /// true if the unit has deep decreases
 bool has_deep_decrease( void ) const { return( ! v_deep_threshold.empty() ); }

 /// true if the direction of a modulation step matters
 /** That is, if the modulations may last more than one instant, or the
  * downward steps have a cost, or the unit has deep decreases; the
  * variables d[ t ] exist exactly in this case. */

 bool has_modulation_direction( void ) const {
  return( ( f_max_modulation_length > 1 ) ||
          ( ! v_down_modulation_cost.empty() ) || has_deep_decrease() );
  }

/** @} ---------------------------------------------------------------------*/
/*--------- METHODS FOR READING THE Variable OF THE NuclearUnitBlock -------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the Variable of the NuclearUnitBlock
 * @{ */

 /// returns the vector of modulation variables

 ColVariable * get_modulation( void ) {
  if( v_modulation.empty() )
   return( nullptr );
  return( &( v_modulation.front() ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// like get_modulation(), but returns a const * so that it can be const

 const ColVariable * get_const_modulation( void ) const {
  return( const_cast< NuclearUnitBlock * >( this )->get_modulation() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the vector of downward modulation variables, if any
 ColVariable * get_modulation_down( void ) {
  if( v_modulation_down.empty() )
   return( nullptr );
  return( &( v_modulation_down.front() ) );
  }

 /// like get_modulation_down(), but const
 const ColVariable * get_const_modulation_down( void ) const {
  return( const_cast< NuclearUnitBlock * >( this )->get_modulation_down() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the vector of deep-decrease variables, if any
 ColVariable * get_deep_decrease( void ) {
  if( v_deep.empty() )
   return( nullptr );
  return( &( v_deep.front() ) );
  }

 /// like get_deep_decrease(), but const
 const ColVariable * get_const_deep_decrease( void ) const {
  return( const_cast< NuclearUnitBlock * >( this )->get_deep_decrease() );
  }

 /// returns the vector of the Variable "decrease by more than the gradient"
 ColVariable * get_deep_drop( void ) {
  return( v_deep_drop.empty() ? nullptr : &( v_deep_drop.front() ) );
  }

 /// like get_deep_drop(), but const
 const ColVariable * get_const_deep_drop( void ) const {
  return( const_cast< NuclearUnitBlock * >( this )->get_deep_drop() );
  }

 /// returns the vector of the Variable "output below the threshold"
 ColVariable * get_deep_low( void ) {
  return( v_deep_low.empty() ? nullptr : &( v_deep_low.front() ) );
  }

 /// like get_deep_low(), but const
 const ColVariable * get_const_deep_low( void ) const {
  return( const_cast< NuclearUnitBlock * >( this )->get_deep_low() );
  }

 /// returns the vector of the modulation start Variable, if any
 ColVariable * get_modulation_start( void ) {
  return( v_modulation_start.empty() ? nullptr
                                     : &( v_modulation_start.front() ) );
  }

/** @} ---------------------------------------------------------------------*/
/*------------------ METHODS FOR SAVING THE NuclearUnitBlock ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the NuclearUnitBlock
 *  @{ */

/// Extends Block::serialize( netCDF::NcGroup )
/** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
 * NuclearUnitBlock (extending that of ThermalUnitBlock). See
 * NuclearUnitBlock::deserialize( netCDF::NcGroup ) (and therefore
 * ThermalUnitBlock::deserialize( netCDF::NcGroup )) for
 * details of the format of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/** @} ---------------------------------------------------------------------*/
/*-------------- METHODS FOR INITIALIZING THE NuclearUnitBlock -------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the data of the NuclearUnitBlock
 *  @{ */

 void load( std::istream & input , char frmt = 0 ) override {
  throw( std::logic_error( "NuclearUnitBlock::load() not implemented yet" ) );
  }

/** @} ---------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for changing the data of the NuclearUnitBlock
 *  @{ */

 /* Method for handling Modification.
  *
  * This method has to intercept any "abstract Modification" that
  * modifies the "abstract representation" of the NuclearUnitBlock, and
  * "translate" them into both changes of the actual data structures and
  * corresponding "physical Modification". These Modification are those
  * for which Modification::concerns_Block() is true.
  *
  * Not needed: NuclearUnitBlock does not allow changes via (its part of)
  * the "abstract representation", any Modification produced in that way
  * will pass through ThermalUnitBlock::add_Modification(), not be
  * recognised as valid, and throw exception. */

 // void add_Modification( sp_Mod mod , ChnlName chnl = 0 ) override;

/*--------------------------------------------------------------------------*/
 /// update a subset of the modulation ramp up values of the unit
 /** This method updates the modulation ramp up values of the unit. The
  * \p subset parameter contains a list of time instants and \p values
  * contains the modulation ramp up values of the unit at those time
  * instants. The modulation ramp up values of the unit at time subset[ i ]
  * is given by std::next( values , i ) for each i in
  * { 0 , ..., subset.size() - 1 }.
  *
  * Recall that the modulation ramp up values must be non-negative and not
  * larger than the original ramp up values, otherwise an exception is
  * thrown.  */

 void set_modulation_ramp_up( MF_dbl_it values ,
			      Subset && subset , bool ordered = false ,
			      ModParam issuePMod = eNoBlck ,
			      ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// update a range of the modulation ramp up values of the unit
 /** This method updates the modulation ramp up values of the unit. The
  * \p rng parameter contains a range of time instants and \p values contains
  * the modulation ramp up values of the unit at those time instants. The
  * modulation ramp up values of the unit at time rng.first + i is given by
  * std::next( values , i ) for each i in {0, ...,
  * ( std::min( rng.second, get_time_horizon() ) - rng.first - 1 )}.
  *
  * Recall that the modulation ramp up values must be non-negative and not
  * larger than the original ramp up values, otherwise an exception is
  * thrown.  */

 void set_modulation_ramp_up( MF_dbl_it values , Range rng = INFRange ,
			      ModParam issuePMod = eNoBlck ,
			      ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// update a subset of the modulation ramp down values of the unit
 /** This method updates the modulation ramp down values of the unit. The
  * \p subset parameter contains a list of time instants and \p values
  * contains the modulation ramp down values of the unit at those time
  * instants. The modulation ramp down values of the unit at time subset[ i ]
  * is given by std::next( values , i ) for each i in
  * { 0 , ..., subset.size() - 1 }.
  *
  * Recall that the modulation ramp down values must be non-negative and not
  * larger than the original ramp down values, otherwise an exception is
  * thrown.  */

 void set_modulation_ramp_down( MF_dbl_it values ,
				Subset && subset , bool ordered = false ,
				ModParam issuePMod = eNoBlck ,
				ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// update a range of the modulation ramp down values of the unit
 /** This method updates the modulation ramp down values of the unit. The
  * \p rng parameter contains a range of time instants and \p values contains
  * the modulation ramp down values of the unit at those time instants. The
  * modulation ramp down values of the unit at time rng.first + i is given by
  * std::next( values , i ) for each i in {0, ...,
  * ( std::min( rng.second, get_time_horizon() ) - rng.first - 1 )}.
  *
  * Recall that the modulation ramp down values must be non-negative and not
  * larger than the original ramp down values, otherwise an exception is
  * thrown.  */

 void set_modulation_ramp_down( MF_dbl_it values , Range rng = INFRange ,
				ModParam issuePMod = eNoBlck ,
				ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// derive the auxiliary Variable of the operating rules
 /** Extends ThermalUnitBlock::set_solution(): once the active power, the
  * commitment and the modulation indicators are set, the start of a
  * modulation and the three indicators of the deep decreases follow from
  * them, and this method computes them, so that a Solution that carries
  * \f$ ( p , u , m , m^- ) \f$ restores a consistent state of the whole
  * abstract representation [see NuclearUnitBlockSolution::write()]. */

 void set_solution( void ) override;

/** @} ---------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*------------------- PROTECTED METHODS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 // void guts_of_add_Modification( p_Mod mod , ChnlName chnl );
 // not required, changes via the "abstract representation" are not supported

 /// verify whether the data in this NuclearUnitBlock is consistent
 /** This function checks whether the data in this NuclearUnitBlock is
  * consistent. The data is consistent if all the following conditions are
  * met.
  *
  * - The delta ramp-up and delta ramp-down vectors are not empty: they are
  *   inherited from ThermalUnitBlock and are mandatory in the nuclear
  *   specialization.
  *
  * - The modulation interval (ModulationTime) is at least 2.
  *
  * - The initial modulation (InitModulation) is at least 1.
  *
  * - For each time step t, the modulation ramp-up and ramp-down are
  *   nonnegative and do not exceed the corresponding ThermalUnitBlock
  *   ramp-up:
  *   \f$ 0 \leq v\_modulation\_ramp\_up[t] \leq v\_DeltaRampUp[t] \f$ and
  *   \f$ 0 \leq v\_modulation\_ramp\_down[t] \leq v\_DeltaRampUp[t] \f$.
  *
  * If any of the above conditions are not met, an exception is thrown. */

 void check_data_consistency( void ) const;

 void update_initial_power_in_cnstrs( c_ModParam issueAMod = eNoBlck )
  override;

 /// generate the constraints of the operating rules
 /** Generates the constraints of the operating rules that the original
  * model does not have (see generate_abstract_constraints()): the direction
  * of the modulation, the modulations lasting more than one instant, the
  * daily limits and the deep decreases, each group only when the data asks
  * for it.
  *
  * Three of these groups have a tight form, which the bit TightRules of the
  * formulation code adds to (in the case of the deep decreases, substitutes
  * for) the default rows; the two describe the same schedules, and differ in
  * their continuous relaxation:
  *
  * - the stability, i.e., that no modulation starts in the \f$ \tau^M - 1
  *   \f$ instants that follow the end of one, which the rows
  *   \f$ \sum_{h = 2}^{K} m_{t+h} + ( K - 1 )( m_t - m_{t+1} ) \leq K - 1
  *   \f$ state on the steps, the end of the modulation being the difference
  *   \f$ m_t - m_{t+1} \f$ and therefore easily fractional. The tight form
  *   adds two families on the starts: the clique \f$ \sum_{h = t}^{t +
  *   \tau^M - 1} s_h \leq 1 \f$, two starts being at least \f$ \tau^M \f$
  *   instants apart, and \f$ m_t - m_{t+1} + \sum_{h = t+1}^{t + \tau^M - 1}
  *   s_h \leq 1 \f$, which is the rule itself with the window aggregated in
  *   one row. Neither implies the rule alone: after a modulation of \f$ k \f$
  *   instants the next start comes \f$ k + \tau^M - 1 \f$ instants after the
  *   previous one, not \f$ \tau^M \f$;
  *
  * - the longest modulation \f$ L^M \f$: besides \f$ \sum_{h = t}^{t + L^M}
  *   m_h \leq L^M \f$, the tight form adds that each step belongs to a
  *   modulation started in the last \f$ L^M \f$ instants, \f$ m_t \leq
  *   \sum_{h = (t - L^M + 1)^+}^{t} s_h \f$;
  *
  * The bit TightRamp does the same for the rows of the full ramp, whose
  * single big-M \f$ M_t = \Delta^+_t + \max\{ \Delta^-_t , \bar{u}_t \} \f$
  * covers every instant that is not a non-final step of a modulation: with
  * the bit on, a stable instant only pays \f$ \Delta^+_t + \Delta^{M-}_t \f$
  * (nothing more than the stability ramp, hence \f$ \Delta^+_t \f$ in the
  * rules where a stable output does not change), a downward step pays
  * \f$ \Delta^+_t + \Delta^-_t \f$, and the shut-down limit is paid on the
  * shut-down Variable alone, where the output really does fall to zero.
  *
  * - the deep decreases: the tight form uses the minimum power in the row
  *   of the threshold, \f$ p_t + ( \tilde{p}_t - \underline{p}_t )^+
  *   dd''_t - \tilde{p}_t u_t \geq 0 \f$ rather than \f$ p_t +
  *   \tilde{p}_t dd''_t \geq \tilde{p}_t \f$, and adds that a deep decrease
  *   is a downward modulation step, \f$ dd_t \leq d_t \f$, which holds
  *   whenever the deep gradient exceeds the stability ramp. */

 void generate_operating_rules( void );

 /// separate the tight rows of the operating rules
 /** Adds the tight rows of TightRules that the current values of the
  * Variable violate, as the ThermalUnitBlock does with the Perspective Cuts
  * (whose separation this method calls first). Only the three families
  * written on the start Variable are separated, each row at most once, and
  * only if the bits TightRules and TightCuts are both on; with TightCuts off
  * they are written by generate_operating_rules() and this method has
  * nothing of its own to do. */

 void generate_dynamic_constraints( Configuration * dycc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// the number of Variable appended to the Objective [see generate_objective]
 Index objective_tail( void ) const override { return( v_obj_tail.size() ); }

 /// accepts a change of the appended coefficients only if it does not
 /// change them
 void objective_tail_change( const DQuadFunction * qf , Index first ,
                             Index last ) override;

 /// the big-M of the full ramp up at t: D+_t + max{ D-_t , SD_t }
 double full_ramp_up_M( Index t ) const {
  return( v_DeltaRampUp[ t ] +
          std::max( v_DeltaRampDown[ t ] , v_ShutDownLimit[ t ] ) );
  }

 /// the big-M of the full ramp down at t: D-_t + max{ D+_t , SU_t }
 double full_ramp_down_M( Index t ) const {
  return( v_DeltaRampDown[ t ] +
          std::max( v_DeltaRampUp[ t ] , v_StartUpLimit[ t ] ) );
  }

 /// the constant of the row of the full ramp up at t, whichever its form
 /** The constant that the row of the full ramp up at t carries, which is
  * full_ramp_up_M() with the default form and the (smaller) one of the
  * stable instant, \f$ \Delta^+_t + \Delta^{M-}_t \f$, with TightRamp; it
  * is what the RHS of the row at t = 0 has to be updated with when the
  * initial power changes. */

 double full_ramp_up_const( Index t ) const {
  return( f_tight_ramp ? v_DeltaRampUp[ t ] + v_modulation_ramp_down[ t ]
                       : full_ramp_up_M( t ) );
  }

 /// the constant of the row of the full ramp down at t, whichever its form
 double full_ramp_down_const( Index t ) const {
  return( f_tight_ramp ? v_DeltaRampDown[ t ] + v_modulation_ramp_up[ t ]
                       : full_ramp_down_M( t ) );
  }

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

//-------------------------------- data -------------------------------------

 /// the vector of modulation delta ramp up
 std::vector< double > v_modulation_ramp_up;

 /// the vector of modulation delta ramp down
 std::vector< double > v_modulation_ramp_down;

 /// the initial modulation value
 int f_initial_modulation;

 /// the minimum allowed modulation interval
 int f_modulation_interval;

 /// the maximum number of instants of a modulation
 Index f_max_modulation_length{ 1 };

 /// the instants of stability that follow a start-up
 Index f_stability_after_start{};

 /// the two breakpoints of the bands of the output, empty if not banded
 std::vector< double > v_power_bands;

 /// the maximum number of modulations per day, -1 if unlimited
 int f_modulations_per_day{ -1 };

 /// the maximum number of start-ups per day, -1 if unlimited
 int f_start_ups_per_day{ -1 };

 /// the maximum number of deep decreases per day, -1 if unlimited
 int f_deep_decreases_per_day{ -1 };

 /// the number of instants of a day, 0 if the horizon is a single day
 Index f_day_length{ 0 };

 /// the cost of the downward modulation steps (empty if 0)
 std::vector< double > v_down_modulation_cost;

 /// the deep-decrease threshold (empty if there are no deep decreases)
 std::vector< double > v_deep_threshold;

 /// the deep-decrease gradient (empty if there are no deep decreases)
 std::vector< double > v_deep_gradient;

 /// the cost of a deep decrease (empty if 0)
 std::vector< double > v_deep_cost;

 /// the coefficients of the Variable appended to the Objective
 std::vector< double > v_obj_tail;

 /// true if the operating rules use their tight rows [see TightRules]
 bool f_tight_rules = false;

 /// true if the full ramp uses its tight big-M [see TightRamp]
 bool f_tight_ramp = false;

 /// true if the tight rows are separated rather than written [see TightCuts]
 bool f_tight_cuts = false;

 /// which of the tight rows have been separated already, three per instant
 std::vector< char > v_cut_done;

//----------------------------- Variable ------------------------------------

 /// the modulation (binary) variables
 std::vector< ColVariable > v_modulation;

 /// the downward modulation (binary) variables, if the direction matters
 std::vector< ColVariable > v_modulation_down;

 /// the modulation start ( [ 0 , 1 ] ) variables, if needed
 std::vector< ColVariable > v_modulation_start;

 /// the band of the output, three binaries per instant, if banded
 std::vector< ColVariable > v_band;

 /// the end of a modulation, one ( [ 0 , 1 ] ) variable per instant, if
 /// the bands need it
 std::vector< ColVariable > v_modulation_end;

 /// the deep-decrease (binary) variables, and the two auxiliary ones
 std::vector< ColVariable > v_deep;
 std::vector< ColVariable > v_deep_drop;
 std::vector< ColVariable > v_deep_low;

//---------------------------- Constraint -----------------------------------

 /// the Modulation RampUp time constraints
 std::vector< FRowConstraint > Modulation_RampUp_Constraints;

 /// the Modulation RampDown time constraints
 std::vector< FRowConstraint > Modulation_RampDown_Constraints;

 /// the NoDownModulation constraints
 std::vector< FRowConstraint > NoDownModulation;

 /// the NoStartUpModulation constraints
 std::vector< FRowConstraint > NoStartUpModulation;

 /// the Modulation constraints constraints proper
 std::vector< FRowConstraint > ModulationConst;

 /// d_t <= m_t, if the direction matters
 std::vector< FRowConstraint > DownModulationLink;

 /// the same direction along a modulation, if L^M > 1
 std::vector< FRowConstraint > ModulationSameDirection;

 /// the full ramp up / down but at the last step of a modulation, if L^M > 1
 std::vector< FRowConstraint > Modulation_FullRampUp;
 std::vector< FRowConstraint > Modulation_FullRampDown;

 /// the stability after a modulation, if L^M > 1
 std::vector< FRowConstraint > ModulationStability;

 /// the stability after a start-up, if there is any
 std::vector< FRowConstraint > StartUpStability;

 /// the maximum length of a modulation, if L^M > 1
 std::vector< FRowConstraint > ModulationMaxLength;

 /// s_t >= m_t - m_{t-1}, if needed
 std::vector< FRowConstraint > ModulationStartLink;

 /// the daily limits on the modulations, the start-ups and the deep
 /// decreases, one per day
 std::vector< FRowConstraint > ModulationsPerDayConst;
 std::vector< FRowConstraint > StartUpsPerDayConst;
 std::vector< FRowConstraint > DeepDecreasesPerDayConst;

 /// the definition of the deep decreases
 std::vector< FRowConstraint > DeepLowConst;
 std::vector< FRowConstraint > DeepDropConst;
 std::vector< FRowConstraint > DeepLinkConst;

 /// a deep decrease is a downward modulation step, with the tight rows
 std::vector< FRowConstraint > DeepDownLink;

 /// the bands of the output: one band per on instant, the output in it,
 /// the band unchanged but at the end of a modulation, and the end of a
 /// modulation moving it to an adjacent one
 std::vector< FRowConstraint > BandChoice;
 std::vector< FRowConstraint > BandPower;
 std::vector< FRowConstraint > BandKeep;
 std::vector< FRowConstraint > BandMove;
 std::vector< FRowConstraint > ModulationEndLink;

 /// the tight rows that have been separated [see TightCuts]
 std::list< FRowConstraint > Nuclear_cuts;

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

private:

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE FIELDS -------------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/
 /// register the methods in the methods factory

 static void static_initialization( void )
 {
  register_method< NuclearUnitBlock , MF_dbl_it , Subset && , bool >(
   "NuclearUnitBlock::set_modulation_ramp_up" ,
   & NuclearUnitBlock::set_modulation_ramp_up );

  register_method< NuclearUnitBlock , MF_dbl_it , Range >(
   "NuclearUnitBlock::set_modulation_ramp_up" ,
   & NuclearUnitBlock::set_modulation_ramp_up );

  register_method< NuclearUnitBlock , MF_dbl_it , Subset && , bool >(
   "NuclearUnitBlock::set_modulation_ramp_down" ,
   & NuclearUnitBlock::set_modulation_ramp_down );

  register_method< NuclearUnitBlock , MF_dbl_it , Range >(
   "NuclearUnitBlock::set_modulation_ramp_down" ,
   & NuclearUnitBlock::set_modulation_ramp_down );
 }

/*--------------------------------------------------------------------------*/
/*--------------------- Methods for handling Solution ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Solution
 *  @{ */

 /// extends ThermalUnitBlock::get_Solution() with the modulation
 /** Extends ThermalUnitBlock::get_Solution() to also save the modulation
  * indicators, which go with the commitment (bit 1 of the Configuration,
  * i.e., wsol & 2) since they are the same kind of information. */

 Solution * get_Solution( Configuration * solc = nullptr ,
                          bool emptys = true ) override;

/*--------------------------------------------------------------------------*/
 /// return the "appropriate" NuclearUnitBlockSolution

 UnitBlockSolution * new_Solution( void ) const override;

/** @} ---------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

};  // end( class( NuclearUnitBlock ) )

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS NuclearUnitBlockMod ------------------------*/
/*--------------------------------------------------------------------------*/

/// Derived class from Modification for modifications to a NuclearUnitBlock
class NuclearUnitBlockMod : public ThermalUnitBlockMod {

 public:

 /// Public enum for the types of ThermalUnitBlockMod
 enum NUB_mod_type {
  eSetModDP = eTUBModLastParam , ///< set modulation delta ramp up
  eSetModDM ,                     ///< set modulation delta ramp down
  eNUBModLastParam  ///< first allowed parameter value for derived classes
  /**< Convenience value to easily allow derived classes to extend the set of
   * types of NuclearUnitBlockMod. */
  };

 /// Constructor, takes the NuclearUnitBlock and the type
 NuclearUnitBlockMod( NuclearUnitBlock * fblock , int type )
  : ThermalUnitBlockMod( fblock , type ) {}

 ///< Destructor, does nothing
 virtual ~NuclearUnitBlockMod() override = default;

 /// prints the NuclearUnitBlockMod
 void print( std::ostream & output ) const override {
  output << "NuclearUnitBlockMod[" << this << "]: ";
  switch( f_type ) {
   case( eSetMaxP ):
    output << "set max power values";
    break;
   case( eSetInitP ):
    output << "set initial power values";
    break;
   case( eSetInitUD ):
    output << "Set initial up/down times";
    break;
   case( eSetAv ):
    output << "Set availability";
    break;
   case( eSetSUC ):
    output << "Set startup costs";
    break;
   case( eSetLinT ):
    output << "Set linear term";
    break;
   case( eSetQuadT ):
    output << "Set quad term";
    break;
   case( eSetConstT ):
    output << "Set constant term";
    break;
   case( eSetModDP ):
    output << "set modulation delta ramp up";
    break;
   case( eSetModDM ):
    output << "set modulation delta ramp down";
    break;
   default:;
   }
  }
 }; // end( class( NuclearUnitBlockMod ) )

/*--------------------------------------------------------------------------*/
/*--------------------- CLASS NuclearUnitBlockRngdMod ----------------------*/
/*--------------------------------------------------------------------------*/
/// derived from NuclearUnitBlockMod for "ranged" changes in a nuclear

class NuclearUnitBlockRngdMod : public NuclearUnitBlockMod {

 public:

 /// constructor: takes the NuclearUnitBlock, the type, and the range
 NuclearUnitBlockRngdMod( NuclearUnitBlock * fblock , int type ,
                          Block::Range rng )
  : NuclearUnitBlockMod( fblock , type ) , f_rng( rng ) {}

 /// destructor, does nothing
 virtual ~NuclearUnitBlockRngdMod() override = default;

 /// accessor to the range
 Block::c_Range & rng( void ) { return( f_rng ); }

 protected:

 /// prints the NuclearUnitBlockRngdMod
 void print( std::ostream & output ) const override {
  NuclearUnitBlockMod::print( output );
  output << "[ " << f_rng.first << ", " << f_rng.second << " )" << std::endl;
  }

 Block::Range f_rng;  ///< the range

 };  // end( class( ThermalUnitBlockRngdMod ) )

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS NuclearUnitBlockSbstMod ---------------------*/
/*--------------------------------------------------------------------------*/
/// derived from NuclearUnitBlockMod for "subset" changes in a nuclear

class NuclearUnitBlockSbstMod : public NuclearUnitBlockMod {

 public:

 /// constructor: takes the NuclearUnitBlock, the type, and the subset
 NuclearUnitBlockSbstMod( NuclearUnitBlock * fblock , int type ,
                          Block::Subset && nms )
  : NuclearUnitBlockMod( fblock , type ) , f_nms( std::move( nms ) ) {}

 /// destructor, does nothing
 virtual ~NuclearUnitBlockSbstMod() override = default;

 /// accessor to the subset
 Block::c_Subset & nms( void ) { return( f_nms ); }

 protected:

 /// prints the NuclearUnitBlockSbstMod
 void print( std::ostream &output ) const override {
  NuclearUnitBlockMod::print( output );
  output << "(# " << f_nms.size() << ")" << std::endl;
  }

 Block::Subset f_nms;  ///< the subset

 };  // end( class( NuclearUnitBlockSbstMod ) )

/*--------------------------------------------------------------------------*/
/*------------------ CLASS NuclearUnitBlockSolution ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// a solution of a NuclearUnitBlock
/** The NuclearUnitBlockSolution class derives from ThermalUnitBlockSolution
 * and adds the only piece of solution information that a nuclear unit has
 * and a thermal one does not, i.e., the modulation indicators m_t. */

class NuclearUnitBlockSolution : public ThermalUnitBlockSolution
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*------------------------------- FRIENDS ----------------------------------*/

 friend NuclearUnitBlock;  ///< make NuclearUnitBlock friend

/*--------- CONSTRUCTING AND DESTRUCTING NuclearUnitBlockSolution ----------*/

 /// constructor, it has nothing to do
 explicit NuclearUnitBlockSolution( void ) : ThermalUnitBlockSolution() {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void deserialize( const netCDF::NcGroup & group ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 ~NuclearUnitBlockSolution() override = default;
 ///< destructor: it is virtual, and empty

/*----- METHODS DESCRIBING THE BEHAVIOR OF A NuclearUnitBlockSolution -----*/

 void read( const Block * block ) override final;

 void write( Block * block ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// serialize a NuclearUnitBlockSolution into a netCDF::NcGroup
 /** Serialize a NuclearUnitBlockSolution into a netCDF::NcGroup. The format
  * is the one of ThermalUnitBlockSolution [cf.
  * ThermalUnitBlockSolution::serialize()], plus:
  *
  * - The variable "Modulation", of type netCDF::NcDouble and indexed over
  *   "TimeHorizon", holding the modulation indicators; the variable is
  *   optional, in that the modulation may not be saved. */

 void serialize( netCDF::NcGroup & group ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 NuclearUnitBlockSolution * scale( double factor ) const override;

 void sum( const Solution * solution , double multiplier ) override;

 NuclearUnitBlockSolution * clone( bool empty = false ) const override;

/*----------- METHODS FOR READING AND WRITING THE SOLUTION -----------------*/

 /// returns the deep-decrease indicators saved in this Solution
 /** Returns the deep-decrease indicators saved in this Solution, an empty
  * vector if they are not saved [see NuclearUnitBlock::get_Solution()]. */

 [[nodiscard]] const std::vector< double > & get_deep_decrease( void ) const {
  return( v_deep );
  }

/*--------------------------------------------------------------------------*/
 /// returns the modulation indicators saved in this Solution
 /** Returns the modulation indicators saved in this Solution, an empty
  * vector if they are not saved. */

 [[nodiscard]] const std::vector< double > & get_modulation( void ) const {
  return( v_modulation );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// sets the modulation indicators saved in this Solution
 /** Sets the modulation indicators saved in this Solution, which is what a
  * Solver filling the Solution out of its own data structures uses [see
  * set_active_power() and the like in UnitBlockSolution]. */

 void set_modulation( std::vector< double > && m ) {
  v_modulation = std::move( m );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// returns the downward modulation indicators saved in this Solution
 /** Returns the downward modulation indicators saved in this Solution, an
  * empty vector if they are not saved (in particular, if the direction of
  * the modulation does not matter for the NuclearUnitBlock). */

 [[nodiscard]] const std::vector< double > & get_modulation_down( void )
  const { return( v_modulation_down ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// sets the downward modulation indicators saved in this Solution
 void set_modulation_down( std::vector< double > && d ) {
  v_modulation_down = std::move( d );
  }

/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/

 void print( std::ostream & output ) const override {
  output << "NuclearUnitBlockSolution [" << this << "]: " << std::endl;
  }

/*---------------------- PRIVATE PART OF THE CLASS -------------------------*/

 private:

/*---------------------------- PRIVATE FIELDS ------------------------------*/

 std::vector< double > v_modulation;  ///< the modulation indicators

 /// the downward modulation indicators, if the direction matters
 std::vector< double > v_modulation_down;

 /// the deep-decrease indicators, and the two auxiliary ones
 std::vector< double > v_deep;
 std::vector< double > v_deep_drop;
 std::vector< double > v_deep_low;

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class( NuclearUnitBlockSolution ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* __NuclearUnitBlock */

/*--------------------------------------------------------------------------*/
/*------------------- End File NuclearUnitBlock.h --------------------------*/
/*--------------------------------------------------------------------------*/
