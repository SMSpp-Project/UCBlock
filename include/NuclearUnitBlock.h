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
 * Copyright &copy; by Antonio Frangioni, Donato Meoli
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
 * from UnitBlock) and implements a "nuclear unit", i.e., a thermal unit (as
 * described in ThermalUnitBlock, whose Variable, Constraint and Objective it
 * keeps unchanged, together with the references given there for each of its
 * formulations) subject to further operating rules that limit the instants
 * in which its output can change, and how. Its composition with the other
 * units into the unit commitment problem is described in \ref ucblock_model.
 * We use the notation of ThermalUnitBlock: the instants are \f$ t \in
 * \mathcal{T} = \{ 0 , \ldots , T - 1 \} \f$, and \f$ p^{ac}_t \f$,
 * \f$ u_t \f$, \f$ v_t \f$ and \f$ w_t \f$ are the active power, the
 * commitment, the start-up and the shut-down at \f$ t \f$, with
 * \f$ v_t = w_t = 0 \f$ for \f$ t < t_0 \f$, where \f$ t_0 \f$ (init_t) is
 * the first instant at which the commitment is free. Also, \f$ u_{-1} = 1 \f$
 * if the unit is on before the horizon (InitUpDownTime \f$ > 0 \f$) and
 * \f$ u_{-1} = 0 \f$ otherwise, while \f$ p_{-1} \f$ is InitialPower if
 * \f$ u_{-1} = 1 \f$ and 0 otherwise. Finally, the ramps \f$ \Delta^+_t \f$
 * and \f$ \Delta^-_t \f$ (DeltaRampUp[ t ], DeltaRampDown[ t ]) bound the
 * change of the output from \f$ t - 1 \f$ to \f$ t \f$. Similarly,
 * \f$ P^{su}_t \f$ (StartUpLimit[ t ]) is the largest output at a start-up at
 * \f$ t \f$ and \f$ P^{sd}_t \f$ (ShutDownLimit[ t ]) the largest output at
 * \f$ t - 1 \f$ of a unit that shuts down at \f$ t \f$.
 *
 * The class relies on the following properties of the base class:
 *
 * - whatever formulation ThermalUnitBlock uses, it has the start-up and
 *   shut-down variables (\f$ v_t \f$ is v_start_up[ t - init_t ] and
 *   \f$ w_t \f$ is v_shut_down[ t - init_t ] for \f$ t \geq t_0 \f$). At the
 *   instants \f$ 0 \leq t < t_0 \f$ (if any) the commitment is instead fixed
 *   to the state of the unit before the horizon, as the minimum up and down
 *   times require;
 *
 * - variables_generated() and constraints_generated() return whether the
 *   (static) Variable and Constraint have been generated already;
 *
 * - the ramps are given, i.e., v_DeltaRampUp and v_DeltaRampDown are
 *   nonempty, which deserialize() checks.
 *
 * We state the rules on the instants in which the output moves. At an instant
 * \f$ t \f$ that is not a start-up, the output of an on unit is either stable
 * or it performs a modulation step. A stable output moves by at most the
 * stability ramps, i.e., \f$ - \Delta^{M-}_t \leq p^{ac}_t - p^{ac}_{t-1}
 * \leq \Delta^{M+}_t \f$ (ModulationDeltaRampDown, ModulationDeltaRampUp,
 * both 0 by default), while a modulation step (\f$ m_t = 1 \f$) is upwards or
 * downwards (\f$ d_t = 1 \f$ for a downward one). A modulation is a maximal
 * sequence of consecutive steps in the same direction, lasting at most
 * \f$ L^M \f$ (MaxModulationLength) instants; all its steps but the last move
 * the output by exactly the full ramp \f$ \Delta^+_t \f$ or
 * \f$ \Delta^-_t \f$, and the last one by at most as much in the same
 * direction. After the end of a modulation no step can happen in the
 * \f$ \tau^M - 1 \f$ instants that follow (\f$ \tau^M \f$ is ModulationTime),
 * and no modulation can start in the \f$ \tau^v \f$ (StabilityAfterStartUp)
 * instants that begin with a start-up. With PowerBands the range of the
 * output is split into three bands. A stable instant keeps the output in its
 * band, while a modulation crosses exactly one boundary of its band in its
 * direction, i.e., the steps that precede its last one keep the output in the
 * band of origin and the last one lands in the adjacent band. Hence, the
 * length of a modulation is fixed by the distance of the output from the
 * breakpoint it crosses (cf. "PowerBands" in deserialize()), and \f$ L^M \f$
 * is only a further cap. Also, no modulation starts upwards from the highest
 * band nor downwards from the lowest one, not even when the end of the
 * horizon cuts it. In each day of \f$ T^{day} \f$ (DayLength) instants at
 * most \f$ N^M \f$ (ModulationsPerDay) modulations can start, at most
 * \f$ N^{su} \f$ (StartUpsPerDay) start-ups and at most \f$ N^{dd} \f$
 * (DeepDecreasesPerDay) deep decreases can happen; a deep decrease is a
 * decrease of the output by more than \f$ \tilde\Delta_t \f$
 * (DeepDecreaseGradient) to an output below \f$ \tilde p_t \f$
 * (DeepDecreaseThreshold). A downward step costs \f$ c^-_t \f$
 * (DownModulationCost) and a deep decrease \f$ c^{dd}_t \f$
 * (DeepDecreaseCost). With the defaults (\f$ L^M = 1 \f$, no bands, no
 * stability after a start-up, no daily limit and no cost) every modulation is
 * a single instant in which the full ramp is allowed, and two of them are at
 * least \f$ \tau^M \f$ instants apart. These rows are written in
 * generate_abstract_constraints() and the costs in generate_objective(),
 * while NuclearUnitExtDPSolver solves the same model by dynamic programming,
 * with the exceptions stated in NuclearUnitExtDPSolver.h (a
 * ReferenceSchedule, a deep decrease fixed to 1 and a fixed auxiliary
 * Variable of the operating rules, which it refuses). Every row of the
 * operating rules uses the operational bounds \f$ P^{mn}_t \f$ and
 * \f$ P^{mx}_t \f$ of ThermalUnitBlock (i.e., those with the availability
 * \f$ \chi_t \f$), as the dynamic programming Solver does.
 *
 * A description by states, in which a modulation is a state lasting
 * \f$ \ell \f$ instants whose first instant leaves the output unchanged,
 * corresponds to the steps above as follows: such a state is a stable instant
 * followed by \f$ \ell - 1 \f$ modulation steps. Hence, a stability of
 * \f$ B \f$ instants after the end of a modulation state corresponds to
 * \f$ \tau^M = B + 2 \f$. Furthermore, a cost paid at every instant of a
 * downward modulation state corresponds to the cost \f$ c^-_t \f$ of the
 * downward steps plus that of the stable instant that precedes the first of
 * them, which the unit does not charge (a data conversion may add it to the
 * cost of the first step). Similarly, a modulation is counted in the day of
 * its first step rather than in that of the stable instant before it. The
 * longest modulation \f$ L^M \f$ is a datum, and any value not smaller than
 * \f$ \lfloor ( \max_t P^{mx}_t - \min_t P^{mn}_t ) / \min_t \min\{
 * \Delta^+_t , \Delta^-_t \} \rfloor + 1 \f$ never binds: the \f$ L^M - 1 \f$
 * full ramps that precede the last step of a longer modulation would move the
 * output by more than the range it has over the horizon. We remark that the
 * range at a single instant is not enough when the bounds change in time,
 * since an output that follows a decreasing maximum power may decrease by
 * many full ramps. Note that the same rules can also describe a thermal unit
 * that is not nuclear (e.g., with DayLength = 0, a limit on the start-ups
 * over the whole horizon), which is then a NuclearUnitBlock as well.
 *
 * A unit may shut down during the stability that follows a modulation or a
 * start-up, since the stability forbids only the modulation steps. Also, a
 * modulation ends at the latest at the instant before a shut-down, since
 * \f$ m_t \leq u_t \f$; with the bands its last step has then to land in the
 * adjacent band, and therefore the unit cannot shut down in the middle of a
 * modulation. The spinning reserves of the unit are those of
 * ThermalUnitBlock, i.e., continuous quantities bounded by a fraction of the
 * output and by the ramp left after the move of the output. In particular,
 * the start of a reserve does not allow the output of a stable instant to
 * move by more than the stability ramps, and no reserve can be held at a step
 * that moves the output by the full ramp.
 *
 * We do not represent the following features:
 *
 * - reserves offered in fixed amounts in discrete reserve states (no
 *   reserve, primary, primary and secondary, the maximum primary reserve of
 *   a nuclear unit), with a secondary reserve that requires the primary one,
 *   and the minimum number of instants a reserve has to last once started;
 *
 * - the changes of the output forced by the start of a reserve, by which a
 *   unit too close to its minimum (maximum) output to hold the reserve moves
 *   upwards (downwards) by at most the full ramp without this being a
 *   modulation;
 *
 * - the rules on the reserves during a modulation (a reserve equal to the
 *   smaller of its values at the start and at the end of the modulation, or a
 *   secondary reserve held only at the first instant of it, depending on the
 *   class of the unit); the only rule is then the one of the ramp left after
 *   the move of the output stated above;
 *
 * - a start-up cost that depends on how long the unit has been off, and the
 *   other features that ThermalUnitBlock does not represent (synchronous
 *   condensers, polyhedral costs; see ThermalUnitBlock);
 *
 * - a modulation, a start-up or a deep decrease in progress at the beginning
 *   of the horizon: the unit is assumed not to be in the middle of a
 *   modulation at instant 0 and to have performed no modulation, start-up and
 *   deep decrease yet in the first day. Its only history is thus
 *   InitModulation, InitUpDownTime and InitialPower;
 *
 * - a stability after a start-up, or after a modulation, that also keeps the
 *   unit on, since the rows only forbid the modulation steps (the unit may
 *   shut down during either; see above).
 *
 * Finally, the class has the following limitations:
 *
 * - ModulationTime, InitModulation and the data of the operating rules
 *   other than the stability ramps and the two costs cannot be changed once
 *   the unit is loaded;
 *
 * - the costs of the downward modulation steps and of the deep decreases
 *   (DownModulationCost and DeepDecreaseCost) can be changed, from the
 *   physical representation [see set_down_modulation_costs()] or from the
 *   abstract one [see generate_objective()], only if the unit is loaded
 *   with a nonzero cost of that kind, a cost absent or all 0 giving no
 *   term in the Objective;
 *
 * - the stability ramps can be changed only from the physical
 *   representation [see set_modulation_ramp_up()], since a change of the
 *   coefficients of the rows via the abstract representation is not
 *   supported (an exception is thrown). */


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
  * default (bit off) is the one described in
  * generate_abstract_constraints(), while with the bit on the stability,
  * the longest modulation, the stability after a start-up and the deep
  * decreases are described by the tight rows discussed there, i.e., the
  * rows (27) and (29) are replaced by their tight forms and the rows (16),
  * (17), (19) and (32) are added, (15) and (18) being kept; they have a
  * stronger continuous relaxation at the price of the (continuous) start
  * Variable and of \f$ O( n \tau^M ) \f$ rows. Both describe the same set of
  * schedules, hence the specialised Solver is indifferent to the choice. */

 static constexpr unsigned char TightRules = 32;

 /// mask for the 7th bit of the formulation code, == 1 if the tight big-M
 /** With this bit the rows of the full ramp of a modulation use one
  * coefficient per case rather than one for all of them: see (12), (13) in
  * generate_abstract_constraints(). It is independent from TightRules, so
  * that
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
  * three families written on the starts, (16), (17) and (19) of
  * generate_abstract_constraints(), are separated instead of being written;
  * the rows of the deep decreases, which are \f$ O( n ) \f$, stay where
  * they are. */

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
 * the NuclearUnitBlock, i.e., that of ThermalUnitBlock plus the data of the
 * operating rules, all of which are optional; the notation is that of the
 * class comment.
 *
 * - The scalar variable "ModulationTime" of type netCDF::NcUint, the
 *   number \f$ \tau^M \geq 2 \f$ that rules the instants between two
 *   modulations: if a modulation ends at \f$ t \f$, then no modulation
 *   step can happen at \f$ t + 1 , \ldots , t + \tau^M - 1 \f$, which for
 *   single-instant modulations (\f$ L^M = 1 \f$) means that at most one
 *   step happens in any \f$ \tau^M \f$ consecutive instants. A value
 *   smaller than 2 is rejected, since the rule would then be void and the
 *   unit is better described by a ThermalUnitBlock. If it is not provided,
 *   then \f$ \tau^M = 2 \f$. In a description by states with a stability of
 *   \f$ B \f$ instants after the end of a modulation state, the value is
 *   \f$ \tau^M = B + 2 \f$ (see the class comment).
 *
 * - The scalar variable "InitModulation" of type netCDF::NcUint, the
 *   number \f$ \tau^M_0 \geq 1 \f$ of instants before instant 0 at which the
 *   last modulation of the unit ended (at which the unit last modulated, if
 *   \f$ L^M = 1 \f$): \f$ \tau^M_0 = 1 \f$ means that it ended at instant
 *   -1, which is the latest possible, the history being written only for
 *   the instants before 0. Then \f$ m_t = 0 \f$ for \f$ 0 \leq t <
 *   \tau^M - \tau^M_0 \f$. If it is not provided, then \f$ \tau^M_0 = \tau^M
 *   \f$, i.e., the unit last modulated "long ago" (possibly never) and can
 *   modulate at instant 0, provided it is on.
 *
 * - The variable "ModulationDeltaRampUp" of type double, either of size 1
 *   or indexed over "NumberIntervals" (or over "TimeHorizon" if
 *   "NumberIntervals" is not provided), with the mapping of the data of
 *   ThermalUnitBlock: if it has size 1 the value is the same at every
 *   instant, otherwise ModulationDeltaRampUp[ i ] is the value at all the
 *   instants of the interval [ ChangeIntervals[ i - 1 ] ,
 *   ChangeIntervals[ i ] ], with ChangeIntervals[ - 1 ] = 0 (if
 *   NumberIntervals <= 1 or NumberIntervals >= TimeHorizon, then
 *   "ChangeIntervals" is not needed, and it is not read). It gives the
 *   stability ramp-up \f$ \Delta^{M+}_t \f$, i.e., the largest increase of
 *   the output from \f$ t - 1 \f$ to \f$ t \f$ of an on unit that does not
 *   perform a modulation step at \f$ t \f$ (a start-up is not a
 *   modulation, hence a modulation step needs the unit on at
 *   \f$ t - 1 \f$). It must be \f$ 0 \leq \Delta^{M+}_t \leq \Delta^+_t \f$;
 *   if it is not provided, then \f$ \Delta^{M+}_t = 0 \f$, i.e., the output
 *   can increase only by a modulation step.
 *
 * - The variable "ModulationDeltaRampDown", with the same shape and
 *   mapping, for the stability ramp-down \f$ \Delta^{M-}_t \f$, the largest
 *   decrease of the output from \f$ t - 1 \f$ to \f$ t \f$ without a
 *   modulation step; it must be \f$ 0 \leq \Delta^{M-}_t \leq \Delta^-_t
 *   \f$, and if it is not provided, then \f$ \Delta^{M-}_t = 0 \f$.
 *
 * - The scalar variable "MaxModulationLength" of type netCDF::NcUint, the
 *   largest number \f$ L^M \geq 1 \f$ of consecutive instants a modulation
 *   may last. If it is not provided, then \f$ L^M = 1 \f$, i.e., every
 *   modulation is a single instant.
 *
 * - The variable "PowerBands" of type netCDF::NcDouble and of size 2 (over
 *   the dimension "NumberPowerBands"), the two breakpoints \f$ P^b_1 < P^b_2
 *   \f$ (checked to be increasing) that split the output of the unit into
 *   the low band \f$ [ P^{mn}_t , P^b_1 ] \f$, the intermediate band
 *   \f$ [ P^b_1 , P^b_2 ] \f$ and the high band \f$ [ P^b_2 , P^{mx}_t ]
 *   \f$, with the operational bounds of ThermalUnitBlock (at an instant
 *   with \f$ \chi_t = 0 \f$, where the unit may be on at the output 0, the
 *   low band is \f$ \{ 0 \} \f$ and the other two are empty). If it is not
 *   provided, the output is not banded, which is the default. With the
 *   bands, a stable instant keeps the output in the band it is in, and a
 *   modulation crosses exactly one boundary of its band, in its own
 *   direction: the steps that precede the last one move the output by the
 *   full ramp and keep it in the band of origin, and the last step lands in
 *   the adjacent band. Hence, with a constant full ramp \f$ \Delta \f$, a
 *   modulation of \f$ k \f$ steps that starts from the output \f$ p \f$
 *   and crosses the breakpoint \f$ P^b \f$ has, with \f$ x = | p - P^b | /
 *   \Delta \f$, \f$ ( k - 1 ) \Delta \leq | p - P^b | \f$ (the full ramps
 *   before the last step stay in the band of origin) and \f$ k \Delta \geq
 *   | p - P^b | \f$ (the last step, of at most \f$ \Delta \f$, reaches the
 *   adjacent band), i.e., \f$ x \leq k \leq x + 1 \f$: it lasts
 *   \f$ \lceil x \rceil \f$ steps if \f$ x \f$ is not an integer, and
 *   \f$ x \f$ or \f$ x + 1 \f$ steps if it is, the last step of the longer
 *   one moving the output by 0 from \f$ P^b \f$, \f$ L^M \f$ being only a
 *   further cap. In particular, an output equal to a breakpoint
 *   (\f$ x = 0 \f$), in the band of which it is the boundary, can pass to
 *   the adjacent band by a single step that moves it by at most
 *   \f$ \Delta \f$, possibly by 0: such a step, which only changes the band,
 *   is a modulation, i.e., it pays \f$ c^-_t \f$ if downward, counts in
 *   \f$ N^M \f$ and starts the stability, in the rows as in the dynamic
 *   programming Solver. For example, with the bands \f$ [ 200 , 268 ] \f$,
 *   \f$ [ 268 , 382 ] \f$, \f$ [ 382 , 450 ] \f$ and \f$ \Delta = 30 \f$, a
 *   downward modulation from 300 is \f$ 300 \to 270 \f$ followed by a last
 *   step into \f$ [ 240 , 268 ] \f$. At the last instant of the horizon a
 *   modulation either ends there and lands, or is cut by the horizon, in
 *   which case it moves the output by the full ramp and keeps it in the
 *   band of origin, with at most \f$ L^M - 1 \f$ steps within the horizon
 *   [see generate_abstract_constraints()]. No row requires that it could
 *   land within \f$ L^M \f$ steps after the horizon, e.g., a downward
 *   modulation that starts five full ramps above the breakpoint with
 *   \f$ L^M = 3 \f$ is admitted if the horizon cuts it after one or two
 *   steps, while within the horizon it could not start. Only the
 *   modulations outwards from the extreme bands, which could never land,
 *   are excluded also there. The dynamic programming Solver admits the
 *   same modulations.
 *   When the reserves of the unit are given as fixed amounts
 *   \f$ R^{pr} \f$ (primary) and \f$ R^{sc} \f$ (secondary), the
 *   natural breakpoints are \f$ P^b_1 = \hat P^{mn} + R^{pr} + R^{sc} \f$
 *   and \f$ P^b_2 = \hat P^{mx} - R^{pr} - R^{sc} \f$, with the nominal
 *   bounds \f$ \hat P^{mn} \f$ and \f$ \hat P^{mx} \f$ (MinPower and
 *   MaxPower), i.e., \f$ P^b_1 \f$ is the lowest output below which the
 *   whole reserve fits (above the minimum output) and \f$ P^b_2 \f$ the
 *   highest one above which it fits (below the maximum output).
 *
 * - The scalar variable "StabilityAfterStartUp" of type netCDF::NcUint, the
 *   number \f$ \tau^v \geq 0 \f$ of instants, the start-up one comprised,
 *   in which a unit that has just started up cannot perform a modulation
 *   step. If it is not provided, then \f$ \tau^v = 0 \f$; the values 0 and
 *   1 both leave only the start-up instant, in which a step is impossible
 *   anyway, a start-up not being a modulation.
 *
 * - The scalar variables "ModulationsPerDay", "StartUpsPerDay" and
 *   "DeepDecreasesPerDay" of type netCDF::NcUint, the largest numbers
 *   \f$ N^M \f$, \f$ N^{su} \f$ and \f$ N^{dd} \f$ of modulations (counted
 *   at their first step), of start-ups and of deep decreases in each day.
 *   Each of them is optional; if it is not provided, the corresponding
 *   number is not limited.
 *
 * - The scalar variable "DayLength" of type netCDF::NcUint, the number
 *   \f$ T^{day} \f$ of instants of a day; the days are the disjoint
 *   intervals \f$ \{ k T^{day} , \ldots , \min\{ ( k + 1 ) T^{day} , T \} - 1
 *   \} \f$, \f$ k = 0 , 1 , \ldots \f$, the last one possibly shorter. If it
 *   is not provided (or it is 0), then the whole horizon is a single day.
 *
 * - The variable "DownModulationCost" of type double, with the shape and
 *   the mapping of "ModulationDeltaRampUp", for the cost \f$ c^-_t \geq 0
 *   \f$ of a downward modulation step at \f$ t \f$. If it is not provided
 *   (or it is all 0), then there is no such cost.
 *
 * - The variables "DeepDecreaseThreshold", "DeepDecreaseGradient" and
 *   "DeepDecreaseCost", with the same shape and mapping, for the threshold
 *   \f$ \tilde p_t \f$, the gradient \f$ \tilde\Delta_t > 0 \f$ and the
 *   cost \f$ c^{dd}_t \geq 0 \f$ of a deep decrease, which happens at an
 *   instant \f$ t \f$ at which the unit is on, with \f$ u_{t-1} = 1 \f$,
 *   when \f$ p^{ac}_{t-1} - p^{ac}_t > \tilde\Delta_t \f$ and \f$ p^{ac}_t <
 *   \tilde p_t \f$ (strict inequalities: a decrease exactly equal to the
 *   gradient, or an output exactly equal to the threshold, is not a deep
 *   decrease). Consecutive instants that satisfy the conditions are each a
 *   deep decrease. The deep decreases exist if and only if both
 *   "DeepDecreaseThreshold" and "DeepDecreaseGradient" are provided;
 *   "DeepDecreaseCost" is optional, the cost being 0 if it is not provided.
 *
 * The ramps "DeltaRampUp" and "DeltaRampDown" of ThermalUnitBlock, which
 * are optional there, are mandatory for a NuclearUnitBlock, and an
 * exception is thrown if they are not in the \p group. */

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
/** Besides those of ThermalUnitBlock (which depend on the formulation, see
 * ThermalUnitBlock::generate_abstract_variables()), NuclearUnitBlock has
 * the binary variables \f$ m_t \f$, \f$ t \in \mathcal{T} \f$ (the group
 * "m_thermal"), which are 1 if the unit performs a modulation step at
 * \f$ t \f$. They are fixed to 0 at the instants \f$ 0 \leq t <
 * \min\{ \tau^M - \tau^M_0 , T \} \f$, in which the end of the last
 * modulation before the horizon forbids a step, and, if the unit is off
 * before the horizon, at the instants \f$ 0 \leq t < t_0 \f$, in which the
 * unit is held off. The following ones, all indexed over
 * \f$ \mathcal{T} \f$, exist only when the operating rules need them:
 *
 * - the binary variables \f$ d_t \f$ ("m_down_nuclear"), which are 1 if the
 *   modulation step at \f$ t \f$ is downwards (\f$ m_t - d_t \f$ being 1 for
 *   an upward one), and which exist if the direction of the modulation
 *   matters, i.e., if \f$ L^M > 1 \f$, or DownModulationCost is nonzero, or
 *   the deep decreases exist, or the output is banded [see
 *   has_modulation_direction()]; they are fixed to 0 wherever \f$ m_t \f$
 *   is;
 *
 * - the variables \f$ s^M_t \in [ 0 , 1 ] \f$ ("m_start_nuclear"), which
 *   are at least 1 when a modulation starts at \f$ t \f$, and which exist if
 *   \f$ L^M > 1 \f$ and either the modulations per day are limited or the
 *   bit TightRules is on;
 *
 * - the binary variables \f$ b^k_t \f$, \f$ k = 1 , 2 , 3 \f$
 *   ("band_nuclear", stored as \f$ 3 T \f$ variables, the one of band
 *   \f$ k \f$ at \f$ t \f$ being in position \f$ ( k - 1 ) T + t \f$, the
 *   code numbering the bands 0, 1, 2), which are 1 if the output at
 *   \f$ t \f$ is in the low, intermediate or high band, and the variables
 *   \f$ e_t \in [ 0 , 1 ] \f$ ("m_end_nuclear"), which are 1 at the last
 *   step of a modulation; both exist if the output is banded [see
 *   get_power_bands()];
 *
 * - the binary variables \f$ \delta_t \f$, \f$ \delta'_t \f$ and
 *   \f$ \delta''_t \f$ ("deep_nuclear", "deep_drop_nuclear",
 *   "deep_low_nuclear"), which indicate a deep decrease at \f$ t \f$, a
 *   decrease of the output by more than the gradient and an output below
 *   the threshold, and which exist if the deep decreases exist; if the unit
 *   is off before the horizon the three are fixed to 0 at \f$ t = 0 \f$,
 *   where the unit has no on predecessor to decrease from.
 *
 * The Configuration parameter is passed to the method of ThermalUnitBlock,
 * which is called first; the bits TightRules, TightRamp and TightCuts of the
 * same int Configuration select the form of the operating rules [see
 * generate_abstract_constraints()]. */

 void generate_abstract_variables( Configuration *stvv = nullptr ) override;

/*--------------------------------------------------------------------------*/
/// generate the static constraints of the NuclearUnitBlock
/** Besides those of ThermalUnitBlock (which depend on the formulation, see
 * ThermalUnitBlock::generate_abstract_constraints(), and are kept as they
 * are), this method generates the rows of the operating rules, in the
 * notation of the class comment; the name of the group of each row is given
 * in parentheses. The rows that are always there come first. For
 * \f$ t \in \mathcal{T} \f$ the two modulation ramps (groups
 * "Modulation_RampUp_Constraints_Nuclear" and
 * "Modulation_RampDown_Constraints_Nuclear") are
 * \f{align*}{
 *   p^{ac}_t - p^{ac}_{t-1} & \leq \Delta^{M+}_t u_{t-1}
 *     + ( \Delta^+_t - \Delta^{M+}_t ) m_t + P^{su}_t v_t , \tag{1} \\
 *   p^{ac}_{t-1} - p^{ac}_t & \leq \Delta^{M-}_t u_t
 *     + ( \Delta^-_t - \Delta^{M-}_t ) m_t + P^{sd}_t w_t , \tag{2}
 * \f}
 * which are replaced, when the direction of the modulation matters (i.e.,
 * when the variables \f$ d_t \f$ exist), by
 * \f{align*}{
 *   p^{ac}_t - p^{ac}_{t-1} & \leq \Delta^{M+}_t ( u_{t-1} - m_t )
 *     + \Delta^+_t ( m_t - d_t ) + P^{su}_t v_t , \tag{3} \\
 *   p^{ac}_{t-1} - p^{ac}_t & \leq \Delta^{M-}_t ( u_t - m_t )
 *     + \Delta^-_t d_t + P^{sd}_t w_t , \tag{4}
 * \f}
 * so that an upward step cannot decrease the output and a downward one
 * cannot increase it. At \f$ t = 0 \f$ the constants \f$ p_{-1} \f$
 * and \f$ u_{-1} \f$ are moved to the right-hand side, and the terms in
 * \f$ v_t \f$ and \f$ w_t \f$ are absent for \f$ t < t_0 \f$. For an on
 * unit that does not start up at \f$ t \f$ these rows give the window
 * \f$ [ - \Delta^{M-}_t , \Delta^{M+}_t ] \f$ of a stable instant and the
 * windows \f$ [ 0 , \Delta^+_t ] \f$ and \f$ [ - \Delta^-_t , 0 ] \f$ of an
 * upward and of a downward step (\f$ [ - \Delta^-_t , \Delta^+_t ] \f$ with
 * (1), (2)), while at a start-up they bound \f$ p^{ac}_t \f$ by
 * \f$ P^{su}_t \f$ and at a shut-down \f$ p^{ac}_{t-1} \f$ by
 * \f$ P^{sd}_t \f$, together with the rows of ThermalUnitBlock. In
 * particular a shut-down at 0 requires \f$ p_{-1} \leq P^{sd}_0 \f$ by (4)
 * (or (2)), the ramps being mandatory here, so that the formulations of
 * ThermalUnitBlock, which differ on this point only without DeltaRampDown,
 * and the dynamic programming Solver agree on it. Then, for
 * \f$ t = t_0 , \ldots , T - 1 \f$,
 * \f[
 *   m_t \leq u_t \tag{5}
 * \f]
 * ("NoDownModulation_Nuclear"), the unit modulating only when it is on, and
 * \f[
 *   m_t + v_t \leq 1 \tag{6}
 * \f]
 * ("NoStartUpModulation_Nuclear"), a start-up not being a modulation step;
 * for \f$ t < t_0 \f$ either \f$ u_t = 1 \f$, and (5) is void, or
 * \f$ u_t = m_t = 0 \f$, and \f$ v_t \f$ does not exist. If \f$ L^M = 1 \f$
 * the stability is the window
 * \f[
 *   \sum_{ h = \max\{ 0 , t - \tau^M + 1 \} }^{ t } m_h \leq 1 \tag{7}
 * \f]
 * ("ModulationConst_Nuclear"), written for \f$ t > \max\{ \tau^M -
 * \tau^M_0 , 0 \} \f$ (\f$ t > \max\{ \tau^M - \tau^M_0 , t_0 \} \f$ if the
 * unit is off before the horizon), the rows of the earlier instants being
 * implied by the variables fixed to 0 [see generate_abstract_variables()].
 *
 * The other rows exist only when the data ask for them. When the direction
 * matters, for \f$ t \in \mathcal{T} \f$,
 * \f[
 *   d_t \leq m_t \tag{8}
 * \f]
 * ("DownModulationLink_Nuclear"). When \f$ L^M > 1 \f$, the window (7) is
 * replaced by the following rows. For \f$ t = 0 , \ldots , T - 2 \f$ the
 * direction stays the same along a modulation,
 * \f[
 *   d_{t+1} - d_t \leq 1 - m_t , \qquad
 *   d_t - d_{t+1} \leq 1 - m_{t+1} \tag{9}
 * \f]
 * ("ModulationSameDirection_Nuclear"); for \f$ t \in \mathcal{T} \f$ every
 * step but the last of a modulation moves the output by the full ramp
 * ("Modulation_FullRampUp_Nuclear", "Modulation_FullRampDown_Nuclear"),
 * \f{align*}{
 *   p^{ac}_t - p^{ac}_{t-1} & \geq \Delta^+_t \hat m_{t+1}
 *     - M_t ( 1 - m_t + d_t ) , \tag{10} \\
 *   p^{ac}_{t-1} - p^{ac}_t & \geq \Delta^-_t \hat m_{t+1}
 *     - M'_t ( 1 - d_t ) , \tag{11}
 * \f}
 * where \f$ M_t = \Delta^+_t + \max\{ \Delta^-_t , P^{sd}_t \} \f$ and
 * \f$ M'_t = \Delta^-_t + \max\{ \Delta^+_t , P^{su}_t \} \f$ (the rows
 * (3), (4) bound the decrease of the output by \f$ \max\{ \Delta^-_t ,
 * P^{sd}_t \} \f$ and its increase by \f$ \max\{ \Delta^+_t , P^{su}_t \}
 * \f$), \f$ \hat m_{t+1} = m_{t+1} \f$ for \f$ t < T - 1 \f$, and at the
 * last instant \f$ \hat m_T = m_{T-1} - e_{T-1} \f$ if the output is
 * banded and \f$ \hat m_T = 0 \f$ otherwise (see (22) below). With the bit
 * TightRamp the two rows use one coefficient per case instead of a single
 * big-M,
 * \f{align*}{
 *   p^{ac}_t - p^{ac}_{t-1} & \geq \Delta^+_t \hat m_{t+1}
 *     - M^{1}_t ( 1 - m_t ) - M^{2}_t d_t - P^{sd}_t w_t , \tag{12} \\
 *   p^{ac}_{t-1} - p^{ac}_t & \geq \Delta^-_t \hat m_{t+1}
 *     - M^{3}_t ( 1 - d_t ) - M^{4}_t ( m_t - d_t ) - P^{su}_t v_t , \tag{13}
 * \f}
 * with \f$ M^{1}_t = \Delta^+_t + \Delta^{M-}_t \f$,
 * \f$ M^{2}_t = \Delta^+_t + \Delta^-_t \f$,
 * \f$ M^{3}_t = \Delta^-_t + \Delta^{M+}_t \f$ and
 * \f$ M^{4}_t = ( \Delta^+_t - \Delta^{M+}_t )^+ \f$ (a stable instant moves
 * the output down by at most the stability ramp, a downward step by the full
 * ramp, and only a shut-down brings it to 0 from as high as \f$ P^{sd}_t \f$);
 * the terms in \f$ w_t \f$ and \f$ v_t \f$ are there for \f$ t \geq t_0 \f$
 * only. Rows (12) and (13) describe the same integer points as (10), (11),
 * as one checks case by case (a stable instant, a continuing and a last step
 * in each direction, a start-up and a shut-down): in each case the row is the
 * full ramp where (10), (11) are, and otherwise its right-hand side is not
 * above the least change of the output that (3), (4) allow, and hence the
 * row is void. If the
 * variables \f$ s^M_t \f$ exist, they are linked to the starts by
 * \f[
 *   s^M_t \geq m_t - m_{t-1} \tag{14}
 * \f]
 * ("ModulationStartLink_Nuclear") for \f$ t \in \mathcal{T} \f$, with
 * \f$ m_{-1} = 0 \f$, the unit not being in the middle of a modulation at
 * the beginning of the horizon. The stability after a modulation is
 * \f[
 *   \sum_{ h = 2 }^{ K_t } m_{t+h} + ( K_t - 1 ) ( m_t - m_{t+1} )
 *   \leq K_t - 1 , \qquad K_t = \min\{ \tau^M - 1 , T - 1 - t \} ,
 *   \tag{15}
 * \f]
 * ("ModulationStability_Nuclear"), written for the \f$ t \f$ with
 * \f$ t + 2 < T \f$ and \f$ K_t \geq 2 \f$ (the other ones being void): a
 * modulation ending at \f$ t \f$ forbids the steps at \f$ t + 2 , \ldots ,
 * t + K_t \f$, the step at \f$ t + 1 \f$ being excluded by the end itself.
 * With the bit TightRules (and \f$ s^M_t \f$) two families on the starts are
 * added, which have a stronger continuous relaxation:
 * \f{align*}{
 *   \sum_{ h = t }^{ t + K - 1 } s^M_h & \leq 1 , \qquad
 *     K = \min\{ \tau^M , T - t \} \geq 2 , \tag{16} \\
 *   m_t - m_{t+1} + \sum_{ h = t + 1 }^{ t + K' } s^M_h & \leq 1 , \qquad
 *     K' = \min\{ \tau^M - 1 , T - 1 - t \} \geq 1 , \tag{17}
 * \f}
 * for \f$ t = 0 , \ldots , T - 2 \f$ ("ModulationStartsApart_Nuclear",
 * "ModulationEndStarts_Nuclear"), i.e., two starts are at least
 * \f$ \tau^M \f$ instants apart, and a modulation ending at \f$ t \f$
 * forbids the starts up to \f$ t + \tau^M - 1 \f$. The rows (16) do not
 * imply (15), since after a modulation of \f$ k \f$ instants the next
 * start comes at least \f$ k + \tau^M - 1 \f$ instants after the previous
 * one, while (17) with (14) does at the integer points (a step at
 * \f$ t + h \f$, \f$ 2 \leq h \leq K_t \f$, after an end at \f$ t \f$
 * needs a start between \f$ t + 2 \f$ and \f$ t + h \f$); (15) is kept
 * all the same, the rows (16), (17) only strengthening the relaxation. The
 * longest modulation is
 * \f[
 *   \sum_{ h = t }^{ t + L^M } m_h \leq L^M \tag{18}
 * \f]
 * ("ModulationMaxLength_Nuclear") for \f$ t = 0 , \ldots , T - 1 - L^M \f$,
 * plus, with the bit TightRules, the condition that each step belongs to a
 * modulation started in the last \f$ L^M \f$ instants,
 * \f[
 *   m_t \leq \sum_{ h = \max\{ 0 , t - L^M + 1 \} }^{ t } s^M_h \tag{19}
 * \f]
 * ("ModulationStepStarted_Nuclear") for \f$ t \in \mathcal{T} \f$. With the
 * bit TightCuts (and TightRules) the rows (16), (17) and (19) are not
 * written, but separated by generate_dynamic_constraints() when they are
 * violated, each at most once.
 *
 * If the output is banded, with \f$ P^b_1 < P^b_2 \f$ the two PowerBands
 * and \f$ P^{mn}_t \f$, \f$ P^{mx}_t \f$ the operational bounds of
 * ThermalUnitBlock, for \f$ t \in \mathcal{T} \f$ every on instant has one
 * band ("BandChoice_Nuclear"), and the output is in it ("BandPower_Nuclear"):
 * \f{align*}{
 *   & b^1_t + b^2_t + b^3_t = u_t , \tag{20} \\
 *   & P^{mn}_t b^1_t + P^b_1 b^2_t + P^b_2 b^3_t \leq p^{ac}_t
 *     \leq P^b_1 b^1_t + P^b_2 b^2_t + P^{mx}_t b^3_t ; \tag{21}
 * \f}
 * since the band changes only at the last step of a modulation (see (23)
 * below), the steps that precede it keep the output in the band of origin,
 * which fixes the length of the modulation as described for "PowerBands"
 * in deserialize().
 * The variables \f$ e_t \f$ are defined by ("ModulationEnd_Nuclear")
 * \f[
 *   e_t \leq m_t \; ( t \in \mathcal{T} ) , \quad
 *   e_t \leq 1 - m_{t+1} , \;\; e_t \geq m_t - m_{t+1} \; ( t < T - 1 ) ,
 *   \quad \sum_{ h = \max\{ 0 , T - L^M \} }^{ T - 1 } m_h - e_{T-1} \leq
 *   L^M - 1 . \tag{22}
 * \f]
 * Here \f$ e_t \f$ is the last step of a modulation, i.e.,
 * \f$ e_t = m_t ( 1 - m_{t+1} ) \f$ inside the horizon, while at the
 * last instant a modulation either ends (\f$ e_{T-1} = 1 \f$) and lands in
 * its band, or is cut by the horizon (\f$ e_{T-1} = 0 \f$), in which case it
 * keeps the output in the band of origin, its step at \f$ T - 1 \f$ is a
 * full ramp by
 * (10), (11) with \f$ \hat m_T = m_{T-1} - e_{T-1} \f$, and it may have only
 * \f$ L^M - 1 \f$ steps within the horizon, since it still owes its last
 * one; for \f$ L^M = 1 \f$ the last row says that a single-step modulation
 * always ends. For \f$ t = 1 , \ldots , T - 1 \f$ the band does not change
 * unless a modulation ends or the unit restarts ("BandKeep_Nuclear"), it
 * does change when a modulation ends, and only to an adjacent band
 * ("BandMove_Nuclear"), in the direction of the modulation, and no
 * modulation starts upwards from the highest band nor downwards from the
 * lowest one:
 * \f{align*}{
 *   & b^k_t - b^k_{t-1} \leq e_t + 1 - u_{t-1} , \quad
 *     b^k_{t-1} - b^k_t \leq e_t + 1 - u_t ,
 *     \qquad k = 1 , 2 , 3 , \tag{23} \\
 *   & b^k_t + b^k_{t-1} + e_t \leq 2 \;\; ( k = 1 , 2 , 3 ) , \quad
 *     b^1_t + b^3_{t-1} \leq 1 , \quad b^3_t + b^1_{t-1} \leq 1 ,
 *     \tag{24} \\
 *   & b^k_t + b^{k+1}_{t-1} \leq 1 + d_t , \quad
 *     b^{k+1}_t + b^k_{t-1} \leq 2 - d_t , \qquad k = 1 , 2 , \tag{25} \\
 *   & m_t - d_t + b^3_{t-1} \leq 1 , \quad d_t + b^1_{t-1} \leq 1 .
 *     \tag{26}
 * \f}
 * The rows (25) matter where the output lands on a breakpoint that both
 * bands contain. Within the horizon the rows (26) are implied by the others
 * at the integer points, since a modulation that leaves the extreme band
 * outwards could never land; they exclude the modulations that the horizon
 * cuts, which need not land, and which would otherwise move outwards from
 * the extreme band within it. At \f$ t = 0 \f$, if the unit is on before the
 * horizon, the rows (23)-(26) are written with \f$ u_{-1} = 1 \f$ and the
 * constant \f$ b^k_{-1} \f$ equal to 1 for the band \f$ k_0 \f$ of
 * \f$ p_{-1} \f$ and to 0 for the other two, \f$ k_0 \f$ being the
 * lowest band that contains \f$ p_{-1} \f$ (\f$ k_0 = 1 \f$ if
 * \f$ p_{-1} \leq P^b_1 \f$, \f$ k_0 = 2 \f$ if \f$ P^b_1 <
 * p_{-1} \leq P^b_2 \f$, \f$ k_0 = 3 \f$ otherwise), the rows that the
 * constants make void being dropped; if it is off, there is no row at
 * \f$ t = 0 \f$, the band of a start-up being free. Since \f$ k_0 \f$ is
 * the lowest band, an initial power equal to \f$ P^b_1 \f$ is in the low
 * band: at \f$ t = 0 \f$ the unit can neither start a downward modulation
 * nor raise its output above \f$ P^b_1 \f$ by a stable step, which a unit
 * at the same output labeled in the intermediate band at a later instant
 * could (and likewise at \f$ P^b_2 \f$); the dynamic programming Solver
 * labels the initial power in the same way.
 *
 * If \f$ \tau^v \geq 2 \f$, a unit that starts up at \f$ t \f$ cannot
 * perform a step before \f$ t + \tau^v \f$ ("StartUpStability_Nuclear"):
 * for \f$ t = t_0 , \ldots , T - 1 \f$ and \f$ h_t = \min\{ t + \tau^v , T
 * \} \f$,
 * \f[
 *   \sum_{ h = t }^{ h_t - 1 } m_h \leq ( h_t - t ) ( 1 - v_t ) , \tag{27}
 * \f]
 * replaced, with the bit TightRules, by the \f$ h_t - t - 1 \f$ rows
 * \f$ m_h + v_t \leq 1 \f$, \f$ h = t + 1 , \ldots , h_t - 1 \f$ (the term
 * \f$ h = t \f$ of (27) is there for uniformity, (6) already excluding a
 * step at the start-up instant).
 *
 * For each day \f$ \mathcal{D} \f$ (see "DayLength" in deserialize()), the
 * daily limits ("ModulationsPerDay_Nuclear", "StartUpsPerDay_Nuclear",
 * "DeepDecreasesPerDay_Nuclear") are, when the corresponding number is
 * given,
 * \f[
 *   \sum_{ t \in \mathcal{D} } m_t \leq N^M \; \text{ if } L^M = 1 , \quad
 *   \sum_{ t \in \mathcal{D} } s^M_t \leq N^M \; \text{ if } L^M > 1 ,
 *   \quad \sum_{ t \in \mathcal{D} , \, t \geq t_0 } v_t \leq N^{su} ,
 *   \quad \sum_{ t \in \mathcal{D} } \delta_t \leq N^{dd} , \tag{28}
 * \f]
 * the row of the start-ups being absent in a day that ends before
 * \f$ t_0 \f$.
 *
 * If the deep decreases exist, for \f$ t = t^\delta , \ldots , T - 1 \f$,
 * where \f$ t^\delta = 0 \f$ if the unit is on before the horizon and
 * \f$ t^\delta = 1 \f$ otherwise (the instants with an on predecessor),
 * \f{align*}{
 *   p^{ac}_t & \geq \tilde p_t ( 1 - \delta''_t ) , \tag{29} \\
 *   p^{ac}_{t-1} - p^{ac}_t & \leq \tilde\Delta_t + G_t \delta'_t ,
 *     \qquad G_t = ( \max\{ \Delta^-_t , P^{sd}_t \} - \tilde\Delta_t )^+ ,
 *     \tag{30} \\
 *   \delta_t & \geq \delta'_t + \delta''_t + u_t - 2 \tag{31}
 * \f}
 * ("DeepLow_Nuclear", "DeepDrop_Nuclear", "DeepLink_Nuclear"), the decrease
 * of the output being at most \f$ \max\{ \Delta^-_t , P^{sd}_t \} \f$ by
 * (4). Hence \f$ \delta''_t = 1 \f$ whenever \f$ p^{ac}_t < \tilde p_t \f$,
 * \f$ \delta'_t = 1 \f$ whenever \f$ p^{ac}_{t-1} - p^{ac}_t >
 * \tilde\Delta_t \f$, and \f$ \delta_t = 1 \f$ whenever both hold and the
 * unit is on at \f$ t \f$ (a shut-down is not a deep decrease), while at
 * the boundary the rows leave the indicators free, so that a decrease
 * exactly equal to the gradient, or to an output exactly equal to the
 * threshold, is not a deep decrease (\f$ \delta_t \f$, which costs and is
 * counted, is 0 there in an optimal solution, as it is in the dynamic
 * programming Solver). With the bit TightRules (29) is replaced by
 * \f$ p^{ac}_t \geq \tilde p_t u_t - ( \tilde p_t - P^{mn}_t )^+
 * \delta''_t \f$, with the operational \f$ P^{mn}_t \f$, which an on unit
 * satisfies with \f$ \delta''_t = 1 \f$ as the default (29) does (also at
 * an instant with \f$ \chi_t = 0 \f$) and an off unit with any
 * \f$ \delta''_t \f$, and the rows
 * \f[
 *   \delta_t \leq d_t \tag{32}
 * \f]
 * ("DeepDownLink_Nuclear") are added at the instants \f$ t \geq t^\delta \f$
 * with \f$ \tilde\Delta_t > \Delta^{M-}_t \f$, where a deep decrease is
 * larger than what a stable instant allows, hence it is a downward step.
 *
 * The rows are written by build_rows(), after those of ThermalUnitBlock,
 * and the Configuration parameter is that of
 * ThermalUnitBlock::generate_abstract_constraints(). */

 void generate_abstract_constraints( Configuration * stcc = nullptr )
  override {
  ThermalUnitBlock::generate_abstract_constraints( stcc );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// generate the Objective of the NuclearUnitBlock
/** The Objective is that of ThermalUnitBlock [see
 * ThermalUnitBlock::generate_objective()], plus the costs of the downward
 * modulation steps and of the deep decreases, scaled as all the others by
 * the scale factor \f$ \sigma \f$ of the unit:
 * \f[
 *   \sigma \sum_{ t \in \mathcal{T} } ( c^-_t d_t + c^{dd}_t \delta_t ) ,
 *   \tag{33}
 * \f]
 * each of the two sums being there only if the corresponding cost is
 * nonzero (and the variables exist). The corresponding Variable are
 * appended after all those of the ThermalUnitBlock [see
 * ThermalUnitBlock::objective_tail()], in this order. Their coefficients
 * are changed by set_down_modulation_costs() and set_deep_decrease_costs(),
 * or via the abstract representation, in which case the new linear
 * coefficients, divided by \f$ \sigma \f$, become the costs of the
 * physical representation [see objective_tail_change()]; a quadratic term
 * on these Variable is refused with std::invalid_argument. */

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
  *   -# it is (approximately) feasible for the ThermalUnitBlock part;
  *
  *   -# each ColVariable that NuclearUnitBlock adds is feasible;
  *
  *   -# the violation of each Constraint that NuclearUnitBlock adds is not
  *      greater than the tolerance.
  *
  * Every Constraint of this NuclearUnitBlock is a RowConstraint and its
  * violation is given by either the relative (see RowConstraint::rel_viol())
  * or the absolute violation (see RowConstraint::abs_viol()), depending on
  * the Configuration that is provided.
  *
  * The tolerance and the type of violation can be provided by either \p fsbc
  * or f_BlockConfig->f_is_feasible_Configuration and they are determined as
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
  *   - Otherwise, if both f_BlockConfig and
  *     f_BlockConfig->f_is_feasible_Configuration are not nullptr and the
  *     latter is a pointer to either a SimpleConfiguration< double > or to a
  *     SimpleConfiguration< std::pair< double , int > >, then the values of
  *     the parameters are obtained analogously as above;
  *
  *   - Otherwise, by default, the tolerance is Block::DefaultFeasTol and the
  *     relative violation is considered.
  *
  * This function considers only the abstract representation to
  * determine if the solution is feasible. So, the parameter \p useabstract is
  * ignored. If no abstract Variable has been generated, then this
  * function returns true. Moreover, if no abstract Constraint has been
  * generated, the solution is considered to be feasible with respect to the
  * set of Variable only. Notice also that, before checking if the solution
  * satisfies a Constraint, the Constraint is computed
  * (Constraint::compute()).
  *
  * @param useabstract This parameter is ignored.
  *
  * @param fsbc The pointer to a Configuration that specifies the tolerance
  *        and the type of violation that must be considered. */

 bool is_feasible( bool useabstract = false ,
                   Configuration * fsbc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// false: the schedule does not answer for the modulation
 /** Returns false, so that ThermalUnitBlock::is_sol_feasible() leaves the
  * check to the base class, which goes through the Variable and hence through
  * is_feasible(): the schedule of a thermal unit says nothing of the
  * modulation of a nuclear one, nor of the constraints it is in. */

 [[nodiscard]] bool is_sol_feasible_physical( void ) const override {
  return( false );
  }


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
  * anyway [see (27) in generate_abstract_constraints()]. */

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
  * downward steps have a cost, or the unit has deep decreases, or the output
  * is banded, a modulation moving the band in its own direction; the
  * variables d[ t ] exist exactly in this case. */

 bool has_modulation_direction( void ) const {
  return( ( f_max_modulation_length > 1 ) ||
          ( ! v_down_modulation_cost.empty() ) || has_deep_decrease() ||
          has_power_bands() );
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

 /// like get_modulation_start(), but const
 const ColVariable * get_const_modulation_start( void ) const {
  return( const_cast< NuclearUnitBlock * >( this )->get_modulation_start() );
  }

 /// returns the vector of the modulation end Variable, if any
 ColVariable * get_modulation_end( void ) {
  return( v_modulation_end.empty() ? nullptr
                                   : &( v_modulation_end.front() ) );
  }

 /// like get_modulation_end(), but const
 const ColVariable * get_const_modulation_end( void ) const {
  return( const_cast< NuclearUnitBlock * >( this )->get_modulation_end() );
  }

 /// returns the band Variable, band k at t in position k T + t, if any
 ColVariable * get_band( void ) {
  return( v_band.empty() ? nullptr : &( v_band.front() ) );
  }

 /// like get_band(), but const
 const ColVariable * get_const_band( void ) const {
  return( const_cast< NuclearUnitBlock * >( this )->get_band() );
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
  * Not needed: the changes of the coefficients of the Objective reach
  * objective_tail_change() through ThermalUnitBlock::add_Modification(). */

 // void add_Modification( sp_Mod mod , ChnlName chnl = 0 ) override;

/*--------------------------------------------------------------------------*/
 /// changes the modulation ramp up at a subset of the instants
 /** Sets the modulation ramp up \f$ \Delta^{M+}_t \f$ at the instant
  * subset[ i ] to std::next( values , i ), for i = 0, ..., subset.size() - 1;
  * if \p ordered is false the instants are sorted together with their values,
  * an instant given more than once taking its last value. Each new value has
  * to be between 0 and DeltaRampUp at its instant, otherwise std::logic_error
  * is thrown and nothing changes. If the Constraint are generated, every row
  * is written anew from the changed data and those that differ change [see
  * ThermalUnitBlock::update_rows()], i.e., (1) or (3) and, with TightRamp,
  * (13); if this cannot be done in place, the data are restored and
  * std::logic_error is thrown. A dry run of \p issuePMod changes nothing, one
  * of \p issueAMod the data only. */

 void set_modulation_ramp_up( MF_dbl_it values ,
			      Subset && subset , bool ordered = false ,
			      ModParam issuePMod = eNoBlck ,
			      ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// changes the modulation ramp up at the instants of a Range
 /** Sets the modulation ramp up \f$ \Delta^{M+}_t \f$ at the instant
  * rng.first + i to std::next( values , i ), for i = 0, ...,
  * std::min( rng.second , get_time_horizon() ) - rng.first - 1. Each new value
  * has to be between 0 and DeltaRampUp at its instant, otherwise
  * std::logic_error is thrown and nothing changes. If the Constraint are
  * generated, every row is written anew from the changed data and those that
  * differ change [see ThermalUnitBlock::update_rows()], i.e., (1) or (3)
  * and, with TightRamp, (13); if this cannot be done in place, the data are
  * restored and std::logic_error is thrown. A dry run of \p issuePMod changes
  * nothing, one of \p issueAMod the data only. */

 void set_modulation_ramp_up( MF_dbl_it values , Range rng = INFRange ,
			      ModParam issuePMod = eNoBlck ,
			      ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// changes the modulation ramp down at a subset of the instants
 /** Sets the modulation ramp down \f$ \Delta^{M-}_t \f$ at the instant
  * subset[ i ] to std::next( values , i ), for i = 0, ..., subset.size() - 1;
  * if \p ordered is false the instants are sorted together with their values,
  * an instant given more than once taking its last value. Each new value has
  * to be between 0 and DeltaRampDown at its instant, otherwise
  * std::logic_error is thrown and nothing changes. If the Constraint are
  * generated, every row is written anew from the changed data and those that
  * differ change [see ThermalUnitBlock::update_rows()], i.e., (2) or (4)
  * and, with TightRamp, (12); if this cannot be done in place, the data are
  * restored and std::logic_error is thrown. With TightRules the ramp also
  * decides at which instants the rows (32) exist (those where the
  * deep-decrease gradient is larger than it): a change that adds or removes
  * some of them is refused in the same way. A dry run of \p issuePMod changes
  * nothing, one of \p issueAMod the data only. */

 void set_modulation_ramp_down( MF_dbl_it values ,
				Subset && subset , bool ordered = false ,
				ModParam issuePMod = eNoBlck ,
				ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// changes the modulation ramp down at the instants of a Range
 /** Sets the modulation ramp down \f$ \Delta^{M-}_t \f$ at the instant
  * rng.first + i to std::next( values , i ), for i = 0, ...,
  * std::min( rng.second , get_time_horizon() ) - rng.first - 1. Each new value
  * has to be between 0 and DeltaRampDown at its instant, otherwise
  * std::logic_error is thrown and nothing changes. If the Constraint are
  * generated, every row is written anew from the changed data and those that
  * differ change [see ThermalUnitBlock::update_rows()], i.e., (2) or (4)
  * and, with TightRamp, (12); if this cannot be done in place, the data are
  * restored and std::logic_error is thrown. With TightRules the ramp also
  * decides at which instants the rows (32) exist (those where the
  * deep-decrease gradient is larger than it): a change that adds or removes
  * some of them is refused in the same way. A dry run of \p issuePMod changes
  * nothing, one of \p issueAMod the data only. */

 void set_modulation_ramp_down( MF_dbl_it values , Range rng = INFRange ,
				ModParam issuePMod = eNoBlck ,
				ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// changes the cost of a downward modulation step
 /** Changes the cost the unit pays for a downward modulation step at each
  * instant of \p rng. Only a unit loaded with a nonzero DownModulationCost
  * has this cost, and the Objective the corresponding term [see
  * generate_objective()]: for a unit that pays nothing (DownModulationCost
  * absent or all 0) the call does nothing. */

 void set_down_modulation_costs( MF_dbl_it values ,
                                 Subset && subset ,
                                 const bool ordered = false ,
                                 ModParam issuePMod = eNoBlck ,
                                 ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/

 void set_down_modulation_costs( MF_dbl_it values , Range rng = INFRange ,
                                 ModParam issuePMod = eNoBlck ,
                                 ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// changes the cost of a deep decrease
 /** Changes the cost the unit pays for a deep decrease at each instant of
  * \p rng, with the same proviso of set_down_modulation_costs(). */

 void set_deep_decrease_costs( MF_dbl_it values ,
                               Subset && subset ,
                               const bool ordered = false ,
                               ModParam issuePMod = eNoBlck ,
                               ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/

 void set_deep_decrease_costs( MF_dbl_it values , Range rng = INFRange ,
                               ModParam issuePMod = eNoBlck ,
                               ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// derive the auxiliary Variable of the operating rules
 /** Extends ThermalUnitBlock::set_solution(): once the active power, the
  * commitment, the modulation steps \f$ m_t \f$ and their directions
  * \f$ d_t \f$ are set, computes the Variable that follow from them: the
  * starts \f$ s^M_t \f$ of the modulations, the three indicators of the
  * deep decreases and, with the bands, the ends \f$ e_t \f$ and the bands
  * \f$ b^k_t \f$, so that a Solution that carries \f$ ( p^{ac} , u , m ,
  * d ) \f$ restores a consistent state of the whole abstract
  * representation [see NuclearUnitBlockSolution::write()]. Inside the
  * horizon \f$ e_t = m_t ( 1 - m_{t+1} ) \f$, while \f$ e_{T-1} \f$ is 1
  * if the step at \f$ T - 1 \f$ is not a full ramp (or the last
  * \f$ L^M \f$ instants are all steps) and is otherwise decided together
  * with the bands; these are found by a sweep over the instants that keeps,
  * for each band, a sequence of bands that contain the outputs, change only
  * where a modulation ends, to the adjacent band in its direction, and
  * start no modulation outwards from an extreme band, i.e., the rows
  * (20)-(26). If no such sequence exists, which happens only for a schedule
  * that is not feasible, each on instant gets the lowest band that contains
  * its output, so that the rows the schedule breaks are the ones that say
  * why. */

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
  * consistent, i.e., that
  *
  * - the ramps \f$ \Delta^+_t \f$ and \f$ \Delta^-_t \f$ (DeltaRampUp and
  *   DeltaRampDown), which are optional in ThermalUnitBlock, are given;
  *
  * - ModulationTime is at least 2 and InitModulation at least 1;
  *
  * - for each \f$ t \f$ the stability ramps do not exceed the ramps of
  *   ThermalUnitBlock, i.e., \f$ 0 \leq \Delta^{M+}_t \leq \Delta^+_t \f$
  *   and \f$ 0 \leq \Delta^{M-}_t \leq \Delta^-_t \f$;
  *
  * - MaxModulationLength is at least 1;
  *
  * - DownModulationCost and DeepDecreaseCost, if given, are nonnegative,
  *   and DeepDecreaseGradient, if given, is positive.
  *
  * If any of the above conditions is not met, an exception is thrown. The
  * two PowerBands, if given, are checked to be increasing (and to be two)
  * by deserialize(). */

 void check_data_consistency( void ) const;

 /// updates the rows that contain the initial power [see
 /// ThermalUnitBlock::set_initial_power()]
 /** Called by ThermalUnitBlock::set_initial_power() once the abstract
  * representation exists, this method calls the one of ThermalUnitBlock,
  * which writes all the rows anew [see build_rows()], and so updates,
  * besides the rows of ThermalUnitBlock, the rows of the operating rules
  * whose constants contain \f$ p_{-1} \f$, i.e., the rows (1)-(4),
  * (10)-(13) and (30) of generate_abstract_constraints() at \f$ t = 0 \f$.
  * If the output is banded and the unit is on before the horizon, the rows
  * (23)-(26) at \f$ t = 0 \f$ depend on the band \f$ k_0 \f$ of the
  * initial power, which decides which rows exist: a new initial power in
  * the same band needs no change of them, while one in another band cannot
  * be handled, and std::logic_error is thrown before anything is changed,
  * after which set_initial_power() restores the previous InitialPower, so
  * that the physical representation is unchanged. The changes of MaxPower
  * and Availability reach the rows of the operating rules that contain the
  * operational bounds and the limits (the modulation ramps (1)-(4), the
  * full ramps (10)-(13), the band powers (21) and the deep decreases
  * (29)-(30)) through build_rows() in the same way [see
  * ThermalUnitBlock::set_maximum_power()]. */

 void update_initial_power_in_cnstrs( c_ModParam issueAMod = eNoBlck )
  override;

 /// writes the rows of the ThermalUnitBlock and those of the nuclear unit
 /** Calls ThermalUnitBlock::build_rows() and then writes the rows (1)-(32)
  * of generate_abstract_constraints() with the same methods, so that they
  * are generated by generate_abstract_constraints() and compared by
  * ThermalUnitBlock::update_rows() alike [see
  * ThermalUnitBlock::build_rows()]. */

 void build_rows( bool generate_ZOConstraints ) override;

 /// generate the constraints of the operating rules
 /** Generates the rows (8)-(32) of generate_abstract_constraints(), i.e.,
  * those that the model with single-instant modulations does not have: the
  * direction of a modulation, the modulations lasting more than one instant,
  * the bands, the stability after a start-up, the daily limits and the deep
  * decreases, each group only when the data ask for it, and in the form
  * selected by the bits TightRules, TightRamp and TightCuts. Both forms of
  * each group describe the same schedules and differ in their continuous
  * relaxation only. */

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

 /// folds a change of the appended coefficients into the costs they come
 /// from [see set_down_modulation_costs()]
 void objective_tail_change( const DQuadFunction * qf , Index first ,
                             Index last ) override;

 /// makes the appended coefficients follow the scale factor
 void update_objective_tail( const Subset & subset ,
                             c_ModParam issueAMod ) override;

/*--------------------------------------------------------------------------*/
 /// changes one of the two costs of the operating rules
 /** The common guts of set_down_modulation_costs() and of
  * set_deep_decrease_costs(): \p cost is the vector of the costs, \p pos the
  * position at which its coefficients start in the tail of the Objective
  * (v_obj_tail.size() if the unit has no such term), and \p type the type of
  * the Modification to issue. */

 void guts_of_set_rule_costs( MF_dbl_it values , Range rng ,
                              std::vector< double > & cost , Index pos ,
                              int type , ModParam issuePMod ,
                              ModParam issueAMod );

/*--------------------------------------------------------------------------*/
 /// changes one of the two costs at the instants of a Subset
 /** The Subset version of guts_of_set_rule_costs(), with the same
  * parameters. */

 void guts_of_set_rule_costs( MF_dbl_it values , Subset && subset ,
                              bool ordered , std::vector< double > & cost ,
                              Index pos , int type , ModParam issuePMod ,
                              ModParam issueAMod );

 /// the power of the unit at instant -1: InitialPower if the unit is on
 /// before the horizon, 0 if it is off (InitialPower is then ignored)
 double power_before( void ) const {
  return( ( f_InitUpDownTime > 0 ) ? f_InitialPower : 0 );
  }

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
  * stable instant, \f$ \Delta^+_t + \Delta^{M-}_t \f$, with TightRamp
  * [see (10) and (12) in generate_abstract_constraints()]. */

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

 /// with the tight rules, at most one start in any tau^M instants
 std::vector< FRowConstraint > ModulationStartsApart;

 /// with the tight rules, no start in the tau^M - 1 instants after the end
 /// of a modulation
 std::vector< FRowConstraint > ModulationEndStarts;

 /// the stability after a start-up, if there is any:
 /// StartUpStability[ t ] has the rows of a start-up at t
 std::vector< std::vector< FRowConstraint > > StartUpStability;

 /// the maximum length of a modulation, if L^M > 1
 std::vector< FRowConstraint > ModulationMaxLength;

 /// with the tight rules, each step belongs to a modulation started in the
 /// last L^M instants
 std::vector< FRowConstraint > ModulationStepStarted;

 /// s^M_t >= m_t - m_{t-1}, if needed
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
 /// modulation moving it to the adjacent one
 std::vector< FRowConstraint > BandChoice;
 std::vector< FRowConstraint > BandPower;
 /// (the last three have the rows of instant t in their entry [ t ])
 std::vector< std::vector< FRowConstraint > > BandKeep;
 std::vector< std::vector< FRowConstraint > > BandMove;
 std::vector< std::vector< FRowConstraint > > ModulationEndLink;

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

  register_method< NuclearUnitBlock , MF_dbl_it , Subset && , bool >(
   "NuclearUnitBlock::set_down_modulation_costs" ,
   & NuclearUnitBlock::set_down_modulation_costs );

  register_method< NuclearUnitBlock , MF_dbl_it , Range >(
   "NuclearUnitBlock::set_down_modulation_costs" ,
   & NuclearUnitBlock::set_down_modulation_costs );

  register_method< NuclearUnitBlock , MF_dbl_it , Subset && , bool >(
   "NuclearUnitBlock::set_deep_decrease_costs" ,
   & NuclearUnitBlock::set_deep_decrease_costs );

  register_method< NuclearUnitBlock , MF_dbl_it , Range >(
   "NuclearUnitBlock::set_deep_decrease_costs" ,
   & NuclearUnitBlock::set_deep_decrease_costs );

  // k identical copies of this unit, as for ThermalUnitBlock, whose scale()
  // this is [see ThermalUnitBlock::static_initialization()]: without its own
  // name, a nuclear unit would be sized as a thermal one by whoever falls
  // back on the class, and refused by whoever asks for the name

  register_method< NuclearUnitBlock , MF_dbl_it , Subset && , bool >(
   "NuclearUnitBlock::replicate" , & NuclearUnitBlock::scale );

  register_method< NuclearUnitBlock , MF_dbl_it , Range >(
   "NuclearUnitBlock::replicate" , & NuclearUnitBlock::scale );
 }

/*--------------------------------------------------------------------------*/
/*--------------------- Methods for handling Solution ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Solution
 *  @{ */

 /// extends ThermalUnitBlock::get_Solution() with the modulation
 /** Extends ThermalUnitBlock::get_Solution() to also save the modulation
  * indicators and, where they exist, the downward ones and the three
  * deep-decrease indicators, which go with the commitment (bit 1 of the
  * Configuration, i.e., wsol & 2) since they are the same kind of
  * information. */

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
  eSetModCost ,                   ///< set the cost of a modulation step
  eSetDeepCost ,                  ///< set the cost of a deep decrease
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
   case( eSetModCost ):
    output << "set the cost of a modulation step";
    break;
   case( eSetDeepCost ):
    output << "set the cost of a deep decrease";
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
 * and adds the solution information that a nuclear unit has and a thermal
 * one does not, i.e., the modulation indicators \f$ m_t \f$, the
 * downward ones \f$ d_t \f$ (if the direction matters) and the
 * three deep-decrease indicators (if there are deep decreases); the other
 * variables of the operating rules are derived from them [see
 * NuclearUnitBlock::set_solution()]. */

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
  *   "TimeHorizon", holding the modulation indicators \f$ m_t \f$;
  *
  * - the variable "ModulationDown", with the same type and dimension,
  *   holding the downward indicators \f$ d_t \f$;
  *
  * - the variables "DeepDecrease", "DeepDrop" and "DeepLow", with the same
  *   type and dimension, holding \f$ \delta_t \f$, \f$ \delta'_t \f$ and
  *   \f$ \delta''_t \f$.
  *
  * Each of them is optional, and it is written only if it is saved [see
  * NuclearUnitBlock::get_Solution()]. */

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
