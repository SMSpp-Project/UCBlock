/*--------------------------------------------------------------------------*/
/*------------------------- File BatteryUnitBlock.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class BatteryUnitBlock, which derives from UnitBlock
 * [see UnitBlock.h], in order to define a "reasonably standard" Battery
 * storage, E-mobility, Centralized demand response, Distributed load
 * management, Distributed storage, and Power-to-gas units in a single class
 * in the Unit Commitment problem.
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

#ifndef __BatteryUnitBlock
 #define __BatteryUnitBlock
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <cmath>

#include <limits>

#include "ColVariable.h"

#include "FRowConstraint.h"

#include "OneVarConstraint.h"

#include "FRealObjective.h"

#include "UnitBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------------ NAMESPACE ---------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS BatteryUnitBlock -------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// implementation of the Block concept for a storage unit
/** BatteryUnitBlock implements the Block concept [see Block.h] for a unit
 * that draws electrical energy from the network, stores it (possibly after
 * converting it into another energy vector) and delivers it back later,
 * possibly with losses and possibly after part of it has been withdrawn for
 * another use. Physical batteries are of this kind, be they large units at a
 * node of the transmission network (battery storage) or small ones at a node
 * of a distribution network (distributed storage, the batteries of a fleet of
 * electric vehicles). So are the units that store an intermediate energy
 * vector (power-to-gas) and the logical mechanisms that shift a consumption
 * in time, i.e., load shifting and load curtailment, which are called
 * centralized demand response at a node of the transmission network and
 * distributed load management at a node of a distribution network. All of
 * them obey the equations below and differ only in their data, as detailed in
 * the special cases at the end of this description. The unit has one
 * generator, and its place in the model of the complete system (the node
 * balance, the reserve and inertia rows, the scale factor and the investment)
 * is described in \ref ucblock_model.
 *
 * \par Time and sign
 * We index the instants by \f$ t \in \mathcal{T} = \{ 0 , \ldots , T - 1
 * \} \f$, and each datum is given per instant: no length of the time step
 * appears in the rows, and the data given per hour are converted as described
 * below for the storage balance. The active power \f$ p^{ac}_t \f$
 * ("p_battery") is positive when the unit injects power into the network and
 * negative when it absorbs power from it. When the two directions have to be
 * told apart (see below) it is split as \f$ p^{ac}_t = p^+_t - p^-_t \f$ with
 * \f$ p^+_t , p^-_t \geq 0 \f$, where \f$ p^+_t \f$, the outtake
 * ("ol_battery"), is the power the unit delivers while discharging and
 * \f$ p^-_t \f$, the intake ("il_battery"), is the power it draws while
 * charging.
 *
 * \par Storage balance
 * We denote by \f$ v^{ba}_t \f$ ("sl_battery") the energy stored at the end
 * of instant \f$ t \f$, which obeys (6) below, i.e., \f$ v^{ba}_t =
 * \rho^{st}_t v^{ba}_{t-1} + \rho^{in}_t p^-_t - \rho^{out}_t p^+_t -
 * d^{ba}_t \f$. Here \f$ \rho^{in}_t \f$ ("StoringBatteryRho") is the energy
 * that enters the storage per unit of power drawn during instant \f$ t \f$,
 * \f$ \rho^{out}_t \f$ ("ExtractingBatteryRho") the energy that leaves it per
 * unit of power delivered during instant \f$ t \f$,
 * \f$ \rho^{st}_t \in [ 0 , 1 ] \f$ ("StandingBatteryRho") the fraction of
 * the level that is kept from one instant to the next, and
 * \f$ d^{ba}_t \geq 0 \f$ ("Demand") an energy withdrawn during instant
 * \f$ t \f$ for another use (say, the driving of the vehicles). Absent data
 * give \f$ \rho^{in}_t = \rho^{out}_t = \rho^{st}_t = 1 \f$ and
 * \f$ d^{ba}_t = 0 \f$. For charge and discharge efficiencies
 * \f$ \eta^c , \eta^d \in ( 0 , 1 ] \f$, a standing loss \f$ s \f$ per hour,
 * powers in MW, levels in MWh and an instant of \f$ \Delta t \f$ hours, one
 * sets \f$ \rho^{in}_t = \eta^c \Delta t \f$,
 * \f$ \rho^{out}_t = \Delta t / \eta^d \f$ and
 * \f$ \rho^{st}_t = ( 1 - s )^{\Delta t} \f$. This is why each of the first
 * two coefficients may exceed 1, and only the round trip condition
 * \f$ \rho^{in}_t \leq \rho^{out}_t \f$ is required (i.e.,
 * \f$ \eta^c \eta^d \leq 1 \f$: a charge followed by a discharge does not
 * create energy). Likewise, the ramp limits and "Demand" are given per
 * instant, while "MinPower", "MaxPower" and "ConverterMaxPower" are powers
 * and the storage bounds are energies.
 *
 * \par Initial and cyclic conditions
 * In the row of instant \f$ t = 0 \f$ we distinguish two cases by the sign of
 * the datum \f$ V^0 \f$ ("InitialStorage"). If \f$ V^0 \geq 0 \f$, it is the
 * level at the start of the horizon and enters the right-hand side as it is,
 * i.e., \f$ v^{ba}_0 = V^0 + \rho^{in}_0 p^-_0 - \rho^{out}_0 p^+_0 -
 * d^{ba}_0 \f$: the standing loss \f$ \rho^{st}_0 \f$ is not applied to
 * \f$ V^0 \f$, and a level computed at the end of a previous horizon (see
 * set_initial_storage()) is passed on unchanged. This is a convention, by
 * which the loss of instant 0 is not counted on the initial level, while the
 * cyclic form below does apply \f$ \rho^{st}_0 \f$ to \f$ v^{ba}_{T-1} \f$;
 * the two forms are therefore not the same row with different data. If
 * \f$ V^0 < 0 \f$ the level is cyclic: \f$ v^{ba}_{-1} \f$ is replaced by
 * \f$ v^{ba}_{T-1} \f$, i.e., \f$ v^{ba}_0 = \rho^{st}_0 v^{ba}_{T-1} +
 * \rho^{in}_0 p^-_0 - \rho^{out}_0 p^+_0 - d^{ba}_0 \f$, and the value of
 * \f$ V^0 \f$ plays no other role (with \f$ T = 1 \f$ the row has the single
 * term \f$ ( 1 - \rho^{st}_0 ) v^{ba}_0 \f$). Which of the two forms the row
 * has is fixed when the Constraints are generated, and set_initial_storage()
 * refuses a change of sign afterwards. In a formulation where the level is
 * indexed by the start of the instant (the power of instant \f$ t - 1 \f$
 * entering the level of instant \f$ t \f$) the index of the level is shifted
 * by one. Then \f$ V^0 \f$ is the level at the start of the horizon and a
 * bound on it is a condition on the datum, while the bounds (7) apply from
 * the end of the first instant to the end of the last one, and the final
 * level is bounded as well.
 *
 * \par Charging and discharging at the same instant
 * When \f$ \rho^{in}_t = \rho^{out}_t =: \rho_t \f$ for all \f$ t \f$ (in
 * particular, when neither is given) no row distinguishes the two directions;
 * hence, the split is not generated, the balance reads \f$ v^{ba}_t =
 * \rho^{st}_t v^{ba}_{t-1} - \rho_t p^{ac}_t - d^{ba}_t \f$, and the rows on
 * \f$ p^+_t \f$ and \f$ p^-_t \f$ become rows on \f$ p^{ac}_t \f$ (see (5)
 * and (10)). Otherwise the split is generated. A solution with
 * \f$ p^+_t > 0 \f$ and \f$ p^-_t > 0 \f$ can then be replaced by the one
 * that keeps only the net value of the two, which has the same
 * \f$ p^{ac}_t \f$, the same cost and, since
 * \f$ \rho^{in}_t \leq \rho^{out}_t \f$, a level at least as high at
 * \f$ t \f$ and at each later instant. The replacement is feasible unless the
 * higher level violates the upper bound (7) at some instant, i.e., unless the
 * unit has to absorb energy without storing it; with
 * \f$ \rho^{in}_t = \rho^{out}_t \f$ the level is the same and the
 * replacement is always feasible. Hence, charging and discharging at the same
 * instant only pays when absorbing energy without storing it has a value,
 * i.e., when the prices may be negative (and the upper bound of the level can
 * then be binding with energy still to absorb), and when
 * \f$ \rho^{in}_t < \rho^{out}_t \f$ for some \f$ t \f$. In that case only
 * (see generate_abstract_variables()) the binary variables \f$ u^{ch}_t \f$
 * ("b_battery", 1 when charging) and the rows (8) are generated. Note that
 * \f$ \rho^{in}_t < \rho^{out}_t \f$ is the condition
 * \f$ \eta^c \eta^d < 1 \f$ on the efficiencies, which does not depend on
 * \f$ \Delta t \f$.
 *
 * \par Power, C-rates and converter
 * The data \f$ P^{mn}_t \f$ ("MinPower", usually nonpositive, and
 * \f$ - P^{mx}_t \f$ by default) and \f$ P^{mx}_t \f$ ("MaxPower") bound the
 * power the unit draws and delivers together with the reserves in (1) and
 * (2). Moreover, the C-rates \f$ C^{ch} \f$ ("MaxCRateCharge") and
 * \f$ C^{dis} \f$ ("MaxCRateDischarge"), 1 by default, bound the intake and
 * the outtake separately in (5), as fractions of the same data. Finally, the
 * converter, whose power \f$ P^{cv}_t \f$ is "ConverterMaxPower", bounds in
 * (10) the sum \f$ p^-_t + p^+_t \f$ of the power through it in the two
 * directions, or \f$ | p^{ac}_t | \f$ when the split is not generated (the
 * two coincide unless the unit charges and discharges at the same instant); a
 * nonpositive \f$ P^{cv}_t \f$ means no converter limit at \f$ t \f$.
 *
 * \par Price of the energy
 * In the Objective, the datum \f$ b_t \f$ ("Cost") is a price on the net
 * injection: the Objective contains \f$ - \sigma b_t p^{ac}_t = \sigma b_t (
 * p^-_t - p^+_t ) \f$ (see (12)), i.e., \f$ b_t \f$ is paid for the energy
 * drawn and earned for the energy delivered at instant \f$ t \f$, where
 * \f$ \sigma \f$ is the scale factor. A cost \f$ c_t \geq 0 \f$ per unit of
 * energy delivered by a unit that never draws (e.g., the cost of curtailing a
 * load) is therefore given as \f$ b_t = - c_t \f$. If a reference schedule
 * \f$ \hat{p}_t \f$ ("ReferenceSchedule") is given, the term in \f$ b_t \f$
 * is dropped, and the Objective contains instead the deviation \f$ \sigma
 * \sum_t | p^{ac}_t - \hat{p}_t | \f$, linearized by (11).
 *
 * \par Reserves and inertia
 * The primary and secondary reserves \f$ p^{pr}_t , p^{sc}_t \geq 0 \f$
 * ("pr_battery", "sc_battery") exist only if the enclosing UCBlock has
 * primary (respectively, secondary) zones and "MaxPrimaryPower"
 * (respectively, "MaxSecondaryPower") is given. They are symmetric, i.e.,
 * each is counted both as an upward margin in (1) and as a downward one in
 * (2), and they are bounded by \f$ \kappa P^{pr}_t \f$ and
 * \f$ \kappa P^{sc}_t \f$ in (9). Also, the unit contributes nothing to the
 * inertia rows of UCBlock (get_inertia_commitment() and get_inertia_power()
 * return nullptr) and has no fixed consumption.
 *
 * \par Capacity: kappa, scale and design
 * The datum \f$ \kappa \geq 0 \f$ ("Kappa", 1 by default) multiplies
 * \f$ P^{mn}_t \f$, \f$ P^{mx}_t \f$, \f$ P^{pr}_t \f$, \f$ P^{sc}_t \f$,
 * \f$ V^{mn}_t \f$ and \f$ V^{mx}_t \f$ wherever they appear, i.e., it
 * resizes the storage together with its power, while the datum
 * \f$ \kappa^{cv} \geq 0 \f$ ("ConverterKappa", equal to \f$ \kappa \f$ if
 * absent) multiplies \f$ P^{cv}_t \f$ in (10), i.e., it resizes the
 * converter. Thus, if only "Kappa" is given the converter follows the
 * storage, and set_kappa() resizes the two together, as one asset; once
 * \f$ \kappa^{cv} \f$ is given (in the data or by set_converter_kappa()),
 * the storage and the converter are two assets of the same unit, whose
 * sizes are chosen apart (e.g., a large storage with a small converter,
 * which a single multiplier cannot represent unless one of the two costs
 * nothing). The ramp limits, \f$ V^0 \f$, \f$ d^{ba}_t \f$ and
 * "InitialPower" are instead absolute quantities. Hence, with
 * \f$ \kappa < 1 \f$ an initial level above \f$ \kappa V^{mx}_0 \f$
 * forces a discharge at \f$ t = 0 \f$, and the instance is infeasible if
 * not even the largest discharge that (1), (5) and (10) allow at
 * \f$ t = 0 \f$ brings the level within \f$ \kappa V^{mx}_0 \f$. The
 * scale factor \f$ \sigma \f$ ("Scale") represents \f$ \sigma \f$
 * identical copies of the unit that share its Variable (see
 * UnitBlock::scale()). If \f$ c^{inv}_b \f$
 * ("BatteryInvestmentCost") is nonzero, a design variable \f$ x_b \f$
 * ("x_battery") multiplies the data of the storage and of the power in (1),
 * (2), (5) and (7), which take the forms (1d), (2d), (5d) and (7d) below.
 * Likewise, if \f$ c^{inv}_c \f$ ("ConverterInvestmentCost") is nonzero a
 * design variable \f$ x_c \f$ ("x_converter") multiplies \f$ P^{cv}_t \f$ in
 * (10), and the Objective contains \f$ \sigma ( c^{inv}_b x_b + c^{inv}_c x_c
 * ) \f$ for the design variables that exist. A design variable lies in \f$ [
 * \max\{ 0 , \underline{X} \} , \bar{X} ] \f$ if \f$ \bar{X} \geq 0 \f$,
 * where \f$ \underline{X} \f$ and \f$ \bar{X} \f$ are the "MinCapacityDesign"
 * and "MaxCapacityDesign" data of the battery or of the converter. If
 * \f$ \bar{X} < 0 \f$ it is integer in \f$ \{ \max\{ 0 , \underline{X} \} ,
 * \ldots , | \bar{X} | \} \f$, say a number of modules, and binary when
 * \f$ \bar{X} = -1 \f$ (an asset that is built or not); its bounds are always
 * a row, and therefore with both bounds equal to 1 the asset is built. With a
 * battery design variable the binary rows (8) cannot multiply
 * \f$ u^{ch}_t \f$ by \f$ x_b \f$, which would be a bilinear term, and they
 * take the form (8d) below, with the largest design \f$ | \bar X | \f$ of the
 * battery in place of \f$ x_b \f$. Hence, they only make the intake and the
 * outtake exclusive, while the size of either is bounded by (1d), (2d) and
 * (5d) (at \f$ u^{ch}_t = 1 \f$, (8d) gives \f$ p^+_t = 0 \f$ and then (2d)
 * gives \f$ p^-_t \leq - \kappa P^{mn}_t x_b \f$, and symmetrically at
 * \f$ u^{ch}_t = 0 \f$). The reserves are bounded by the rows (9d), which are
 * linear in \f$ x_b \f$, while the bounds (9) are kept with
 * \f$ | \bar X | \f$ in place of \f$ x_b \f$. The converter scales with the
 * battery when it has no separate size. Indeed, with a battery design
 * variable, no converter design variable and no "ConverterMaxPower", each
 * module has a converter of power \f$ \kappa^{cv} P^{cv}_t \f$ with
 * \f$ P^{cv}_t = P^{mx}_t \f$, and (10) becomes (10d) below with
 * \f$ x_b \f$; hence, \f$ n \f$ modules have \f$ n \f$ times the converter
 * power of one (and \f$ \kappa^{cv} = \kappa \f$ unless "ConverterKappa"
 * is given). If instead "ConverterMaxPower"
 * is given and the converter has no design variable, the converter is an
 * asset of the given fixed size, and (10) keeps \f$ \kappa^{cv} P^{cv}_t
 * \f$ without \f$ x_b \f$. The two investment mechanisms, i.e.,
 * \f$ \kappa \f$ and \f$ \kappa^{cv} \f$ changed from outside the Block
 * (see set_kappa() and set_converter_kappa()) and the design variables, are
 * not meant to be combined, and get_kappa_linearization() and
 * get_converter_kappa_linearization() throw if a design variable exists.
 *
 * \par Rows and Objective
 * With these data and variables, the rows are, for all
 * \f$ t \in \mathcal{T} \f$,
 * \f{align*}{
 *   & p^{ac}_t + p^{pr}_t + p^{sc}_t \leq \kappa P^{mx}_t \tag{1} \\
 *   & \kappa P^{mn}_t \leq p^{ac}_t - p^{pr}_t - p^{sc}_t \tag{2} \\
 *   & p^{ac}_t - p^{ac}_{t-1} \leq \Delta^+_t , \qquad
 *     p^{ac}_{t-1} - p^{ac}_t \leq \Delta^-_t \tag{3} \\
 *   & p^{ac}_t = p^+_t - p^-_t \tag{4} \\
 *   & p^-_t \leq - \kappa C^{ch} P^{mn}_t , \qquad
 *     p^+_t \leq \kappa C^{dis} P^{mx}_t \tag{5} \\
 *   & v^{ba}_t = \rho^{st}_t v^{ba}_{t-1} + \rho^{in}_t p^-_t -
 *     \rho^{out}_t p^+_t - d^{ba}_t \tag{6} \\
 *   & \kappa V^{mn}_t \leq v^{ba}_t \leq \kappa V^{mx}_t \tag{7} \\
 *   & p^-_t \leq - \kappa P^{mn}_t u^{ch}_t , \qquad
 *     p^+_t \leq \kappa P^{mx}_t ( 1 - u^{ch}_t ) ,
 *     \qquad u^{ch}_t \in \{ 0 , 1 \} \tag{8} \\
 *   & p^{pr}_t \leq \kappa P^{pr}_t , \qquad
 *     p^{sc}_t \leq \kappa P^{sc}_t \tag{9} \\
 *   & p^-_t + p^+_t \leq \kappa^{cv} P^{cv}_t \tag{10} \\
 *   & - \delta^{rs}_t \leq p^{ac}_t - \hat{p}_t \leq \delta^{rs}_t \tag{11}
 * \f}
 * where the terms in \f$ p^{pr}_t \f$, \f$ p^{sc}_t \f$, \f$ p^+_t \f$,
 * \f$ p^-_t \f$ and \f$ u^{ch}_t \f$ are present only when the corresponding
 * variables are. In (3), \f$ p^{ac}_{-1} \f$ is the datum "InitialPower" (0
 * by default), and each of the two families of rows exists only if
 * \f$ \Delta^+_t \f$ ("DeltaRampUp"), respectively \f$ \Delta^-_t \f$
 * ("DeltaRampDown"), is given. The rows (4), (5) as written and (8) exist
 * only when the split is generated; without it, (5) reads \f$ \kappa C^{ch}
 * P^{mn}_t \leq p^{ac}_t \leq \kappa C^{dis} P^{mx}_t \f$ and (6) has
 * \f$ - \rho^{out}_t p^{ac}_t \f$ in place of the two terms in
 * \f$ p^\pm_t \f$. At \f$ t = 0 \f$ the row (6) has the initial or cyclic
 * form described above, and (9) exists only for the reserves that exist. The
 * row (10), or \f$ - \kappa^{cv} P^{cv}_t \leq p^{ac}_t \leq \kappa^{cv}
 * P^{cv}_t \f$ without the split, exists only if "ConverterMaxPower" or
 * "ConverterKappa" is given (the latter also by set_converter_kappa()
 * before the rows are generated) or a design variable exists (with
 * \f$ P^{cv}_t = P^{mx}_t \f$ by default), and only at the instants where
 * \f$ P^{cv}_t > 0 \f$. Finally, (11), with
 * \f$ \delta^{rs}_t \geq 0 \f$ ("v_absb_refschd"), exists only if
 * "ReferenceSchedule" is given. Note that the levels \f$ v^{ba}_t \f$ and the
 * active power are free variables, whose sign is given by (7) and (1)-(2)
 * only. With a design variable the rows (1), (2), (5) and (7) become
 * \f{align*}{
 *   & p^{ac}_t + p^{pr}_t + p^{sc}_t \leq \kappa P^{mx}_t x_b \tag{1d} \\
 *   & \kappa P^{mn}_t x_b \leq p^{ac}_t - p^{pr}_t - p^{sc}_t \tag{2d} \\
 *   & p^-_t \leq - \kappa C^{ch} P^{mn}_t x_b , \qquad
 *     p^+_t \leq \kappa C^{dis} P^{mx}_t x_b \tag{5d} \\
 *   & \kappa V^{mn}_t x_b \leq v^{ba}_t \leq \kappa V^{mx}_t x_b \tag{7d} \\
 *   & p^-_t \leq - \kappa | \bar X | P^{mn}_t u^{ch}_t , \qquad
 *     p^+_t \leq \kappa | \bar X | P^{mx}_t ( 1 - u^{ch}_t ) \tag{8d} \\
 *   & p^{pr}_t \leq \kappa P^{pr}_t x_b , \qquad
 *     p^{sc}_t \leq \kappa P^{sc}_t x_b \tag{9d}
 * \f}
 * (with \f$ - p^{ac}_t \f$ and \f$ p^{ac}_t \f$ in place of \f$ p^-_t \f$
 * and \f$ p^+_t \f$ in (5d) without the split), and (10) becomes
 * \f[
 *   p^-_t + p^+_t \leq \kappa^{cv} P^{cv}_t x \tag{10d}
 * \f]
 * (\f$ - \kappa^{cv} P^{cv}_t x \leq p^{ac}_t \leq \kappa^{cv} P^{cv}_t x
 * \f$ without the split), where \f$ x = x_c \f$ if the converter has a
 * design variable and \f$ x = x_b \f$ if it has none and
 * "ConverterMaxPower" is not given. An
 * infinite bound of the level stays infinite whatever \f$ \kappa \f$ is (also
 * 0), and the side of (7d) that it would give has no term in \f$ x_b \f$, and
 * the row is free on that side. The Objective is
 * \f[
 *   \sum_{ t \in \mathcal{T} } \sigma \bigl( b_t ( p^-_t - p^+_t ) \bigr) +
 *   \sigma ( c^{inv}_b x_b + c^{inv}_c x_c ) \tag{12}
 * \f]
 * with \f$ - b_t p^{ac}_t \f$ in place of \f$ b_t ( p^-_t - p^+_t ) \f$
 * without the split, \f$ \delta^{rs}_t \f$ in place of it if
 * "ReferenceSchedule" is given, and the investment terms only for the design
 * variables that exist.
 *
 * \par Special cases
 * This class represents the following kinds of units, each obtained by the
 * data indicated.
 *
 * - A battery storage, or a distributed storage when the node is one of a
 *   distribution network, uses all the data above as the physical device
 *   dictates. It gives primary and secondary reserve only if
 *   "MaxPrimaryPower" and "MaxSecondaryPower" are given, and its storage
 *   cycle is short w.r.t. the horizon, hence no value is attached to its
 *   final level (see the features not modeled below).
 *
 * - The batteries of a fleet of electric vehicles (e-mobility), aggregated
 *   into one unit, have \f$ V^{mx}_t \f$ the energy the batteries of the
 *   fleet can hold, \f$ - P^{mn}_t \f$ the power at which the vehicles
 *   connected at \f$ t \f$ can be charged and \f$ P^{mx}_t \f$ the power they
 *   can give back to the network (vehicle-to-grid). If there is none,
 *   \f$ P^{mx}_t = 0 \f$, and "MinPower" has to be given, since its default
 *   is \f$ - P^{mx}_t \f$. The energy the vehicles use for driving during
 *   instant \f$ t \f$ is \f$ d^{ba}_t \f$, a minimum level \f$ V^{mn}_t \f$
 *   can represent the energy needed by the vehicles that leave at \f$ t \f$,
 *   and the efficiencies are those of the chargers. A need to charge a given
 *   energy over a given period, with no other use of the batteries, is
 *   instead a load shifting (see below).
 *
 * - A power-to-gas unit converts electrical energy into a gas that it stores
 *   and either converts back or delivers to other uses. Its level
 *   \f$ v^{ba}_t \f$ is the energy content of the gas stored,
 *   \f$ \rho^{in}_t \f$ is the efficiency of the conversion times
 *   \f$ \Delta t \f$, \f$ P^{mx}_t \f$ is the power of the conversion back (0
 *   if there is none, in which case "MinPower" is given) and \f$ d^{ba}_t \f$
 *   is the gas withdrawn for the other uses.
 *
 * - A load shifting flexibility (centralized demand response at a node of the
 *   transmission network, distributed load management at a node of a
 *   distribution network) is a flexible consumption \f$ d^{ls}_t \geq 0 \f$.
 *   Its energy over each sub-period \f$ \mathcal{T}_q \f$,
 *   \f$ q = 1 , \ldots , n_q \f$, of a partition of \f$ \mathcal{T} \f$ into
 *   sub-periods of consecutive instants is a given \f$ E_q \f$, i.e.,
 *   \f$ \sum_{ t \in \mathcal{T}_q } d^{ls}_t = E_q \f$, and its deviation
 *   from a reference profile \f$ D^{ls}_t \geq 0 \f$ is bounded as
 *   \f$ P^{mn}_t \leq D^{ls}_t - d^{ls}_t \leq P^{mx}_t \f$, where the
 *   profile satisfies \f$ \sum_{ t \in \mathcal{T}_q } D^{ls}_t = E_q \f$ and
 *   is included in the demand of the node. Taking
 *   \f$ p^{ac}_t = D^{ls}_t - d^{ls}_t \f$, it is a BatteryUnitBlock with no
 *   "Cost", no reserve and no efficiency data (the powers are then energies
 *   per instant; with powers in MW, both efficiencies are equal to
 *   \f$ \Delta t \f$), and therefore (6) reads
 *   \f$ v^{ba}_t = v^{ba}_{t-1} - p^{ac}_t \f$. Its other data are "MinPower"
 *   \f$ = P^{mn}_t \f$, "MaxPower" \f$ = P^{mx}_t \f$, "InitialStorage"
 *   \f$ = V^0 := \max_q E_q \f$, "Kappa" \f$ = 1 \f$ (since \f$ \kappa \f$
 *   scales the storage bounds but not \f$ V^0 \f$), and "MinStorage"
 *   \f$ = \f$ "MaxStorage" \f$ = V^0 \f$ at the last instant \f$ t^+_q \f$ of
 *   each sub-period. In fact, (6) gives \f$ v^{ba}_{t^+_q} =
 *   v^{ba}_{t^+_{q-1}} - \sum_{ t \in \mathcal{T}_q } p^{ac}_t \f$, with
 *   \f$ v^{ba}_{t^+_0} := V^0 \f$, and fixing the two levels to \f$ V^0 \f$
 *   is \f$ \sum_{ t \in \mathcal{T}_q } p^{ac}_t = 0 \f$, i.e., the energy
 *   condition. The level is fixed at the last instant of a sub-period only:
 *   since it is taken at the end of the instant, fixing it at the first
 *   instant \f$ t^-_q \f$ as well would force \f$ p^{ac}_{t^-_q} = 0 \f$. At
 *   the other instants of \f$ \mathcal{T}_q \f$ the bounds \f$ V^0 - E_q \f$
 *   and \f$ V^0 + E_q \f$ are valid and cut no feasible profile, since
 *   \f$ v^{ba}_t = V^0 - \sum_{ s \in \mathcal{T}_q , s \leq t } D^{ls}_s +
 *   \sum_{ s \in \mathcal{T}_q , s \leq t } d^{ls}_s \f$ and both partial
 *   sums lie in \f$ [ 0 , E_q ] \f$; the choice of \f$ V^0 \f$ keeps the
 *   lower bounds nonnegative. The condition \f$ d^{ls}_t \geq 0 \f$ follows
 *   from \f$ P^{mx}_t \leq D^{ls}_t \f$, and the ramp limits are optional.
 *
 * - A load curtailment flexibility is a contract that allows the consumption
 *   of the node to be reduced by a total energy \f$ V^0 \f$ over its
 *   duration, by at most \f$ P^{mx}_t \f$ at each instant. It is a
 *   BatteryUnitBlock that only injects, with "InitialStorage" \f$ = V^0 \f$,
 *   "MinStorage" \f$ = 0 \f$, "MaxStorage" \f$ = V^0 \f$, "MinPower"
 *   \f$ = 0 \f$, "MaxPower" \f$ = P^{mx}_t \f$, no efficiency data and, if
 *   the curtailment has a cost \f$ c_t \f$ per unit of energy, "Cost"
 *   \f$ = - c_t \f$. "MinPower" has to be 0: its default \f$ - P^{mx}_t \f$
 *   would let the unit increase the consumption, and a positive value would
 *   force a curtailment at each instant (besides, a positive "MinPower" is
 *   refused whenever the split is generated). No condition is imposed on the
 *   final level, which has no value within the horizon. When the contract
 *   spans several stages of a multistage problem, the level at the end of a
 *   stage is the "InitialStorage" of the next one (see
 *   set_initial_storage()).
 *
 * \par Features not modeled
 * We do not represent a cost of wear charged in both directions,
 * \f$ b_t ( p^+_t + p^-_t ) \f$, since "Cost" is a price on the net
 * injection, nor the contribution of a battery to the inertia, which would be
 * proportional to its signed active power (and then negative while it
 * charges). No value is attached to the energy left at the end of the
 * horizon; therefore, neither a future value of the final level of a battery
 * nor a load curtailment treated as a seasonal storage (valued, jointly with
 * the hydro reservoirs, by a cost-to-go function of its final level) is
 * represented, and the level can only be passed to the next stage with
 * set_initial_storage(). Also, the efficiencies depend neither on the level
 * nor on the power, and the active and the reactive power are not linked. */

class BatteryUnitBlock : public UnitBlock
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
 /** Constructor of BatteryUnitBlock, taking possibly a pointer of its father
  * Block. */

 explicit BatteryUnitBlock( Block * f_block = nullptr ) :
  UnitBlock( f_block ), f_BattInvestmentCost( 0 ), f_ConvInvestmentCost( 0 ),
  f_BattMinCapacityDesign( 0 ), f_BattMaxCapacityDesign( 1 ),
  f_ConvMinCapacityDesign( 0 ), f_ConvMaxCapacityDesign( 1 ),
  f_InitialStorage( 0 ), f_InitialPower( 0 ), f_MaxCRateCharge( 1 ),
  f_MaxCRateDischarge( 1 ), f_kappa( 1 ),
  f_conv_kappa( std::numeric_limits< double >::quiet_NaN() ), f_scale( 1 ),
  f_conv_max_power_given( false ) {}

/*--------------------------------------------------------------------------*/
 /// destructor of BatteryUnitBlock
 virtual ~BatteryUnitBlock() override;

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 * @{ */

 /// extends Block::deserialize( netCDF::NcGroup )
 /** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
  * the BatteryUnitBlock. Besides the mandatory "type" attribute of any
  * :Block, the group must contain all the data required by the base
  * UnitBlock, as described in the comments to UnitBlock::deserialize(
  * netCDF::NcGroup ). In particular, we refer to that description for the
  * crucial dimensions "TimeHorizon", "NumberIntervals" and
  * "ChangeIntervals". The symbols are those of the class description.
  * Several data below are time-indexed: such a variable, of type
  * netCDF::NcDouble, is either a scalar, or indexed over the dimension
  * "NumberIntervals" (one value per interval, see "ChangeIntervals"), or
  * indexed over "TimeHorizon", and it is expanded into one value per instant
  * as described in UnitBlock::deserialize(), so that the vector returned by
  * the corresponding method is either empty (an optional datum that is
  * absent, unless a default is said) or of size get_time_horizon(). The
  * netCDF::NcGroup must then also contain:
  *
  * - the time-indexed variable "MinStorage", the minimum level
  *   \f$ V^{mn}_t \f$ in (7); it must be \f$ V^{mn}_t \leq V^{mx}_t \f$,
  *   equality being allowed (and used by a load shifting flexibility), and
  *   \f$ V^{mn}_t \f$ may be negative, or \f$ - \infty \f$ for a level that
  *   is unbounded below, since the level is a free variable;
  *
  * - the time-indexed variable "MaxStorage", the maximum level
  *   \f$ V^{mx}_t \f$ in (7);
  *
  * - the time-indexed variable "MaxPower", the maximum power
  *   \f$ P^{mx}_t \f$ that the unit delivers, in (1), (5) and (8);
  *
  * and it may contain:
  *
  * - the time-indexed variable "MinPower", the minimum active power
  *   \f$ P^{mn}_t \leq P^{mx}_t \f$ in (2), (5) and (8), i.e., minus the
  *   maximum power that the unit draws; if absent, \f$ P^{mn}_t = -
  *   P^{mx}_t \f$. When the split is generated (see the class description)
  *   it must be \f$ P^{mn}_t \leq 0 \f$, since (5) bounds the intake by
  *   \f$ - \kappa C^{ch} P^{mn}_t \f$;
  *
  * - the scalar variables "MaxCRateCharge" and "MaxCRateDischarge", the
  *   C-rates \f$ C^{ch} \f$ and \f$ C^{dis} \f$ in (5); 1 if absent;
  *
  * - the time-indexed variable "ConverterMaxPower", the power
  *   \f$ P^{cv}_t \f$ of the converter in (10): with no design variable the
  *   rows (10) exist only if it or "ConverterKappa" is given, while with one
  *   of them they always exist, with \f$ P^{cv}_t = P^{mx}_t \f$ if it is
  *   absent, in which case
  *   a battery with a design variable and a converter without one has a
  *   converter per module, i.e., (10d) with \f$ x_b \f$; a nonpositive
  *   \f$ P^{cv}_t \f$ means no converter limit at \f$ t \f$. The default is
  *   not written by serialize();
  *
  * - the scalar variable "InitialStorage", the datum \f$ V^0 \f$ (0 if
  *   absent): if \f$ V^0 \geq 0 \f$ it is the level at the start of the
  *   horizon, which enters (6) at \f$ t = 0 \f$ without the standing loss,
  *   while \f$ V^0 < 0 \f$ makes the level cyclic (see the class
  *   description); it is not multiplied by \f$ \kappa \f$, and a positive
  *   value is refused if an investment cost is nonzero (by deserialize()
  *   only, set_initial_storage() not repeating this check);
  *
  * - the scalar variable "InitialPower", the active power
  *   \f$ p^{ac}_{-1} \f$ before the horizon (0 if absent), used only by the
  *   ramp rows (3) at \f$ t = 0 \f$;
  *
  * - the time-indexed variables "DeltaRampUp" and "DeltaRampDown", the
  *   nonnegative ramp limits \f$ \Delta^+_t \f$ and \f$ \Delta^-_t \f$ of
  *   (3) for the step from \f$ t - 1 \f$ to \f$ t \f$, i.e., changes of
  *   the active power per instant; if one of them is absent the
  *   corresponding rows do not exist. They are not multiplied by
  *   \f$ \kappa \f$;
  *
  * - the time-indexed variables "StoringBatteryRho" and
  *   "ExtractingBatteryRho", the nonnegative coefficients
  *   \f$ \rho^{in}_t \f$ and \f$ \rho^{out}_t \f$ of (6), 1 if absent; it
  *   must be \f$ \rho^{in}_t \leq \rho^{out}_t \f$, while each of them may
  *   exceed 1 (see the class description), and if they coincide at all
  *   instants the active power is not split;
  *
  * - the time-indexed variable "StandingBatteryRho", the coefficient
  *   \f$ \rho^{st}_t \in [ 0 , 1 ] \f$ of (6), 1 if absent;
  *
  * - the time-indexed variable "Demand", the nonnegative energy
  *   \f$ d^{ba}_t \f$ withdrawn from the storage during instant \f$ t \f$
  *   for another use, in (6); 0 if absent, and not multiplied by
  *   \f$ \kappa \f$;
  *
  * - the time-indexed variable "Cost", the price \f$ b_t \f$ on the net
  *   injection in (12), paid for the energy drawn and earned for the energy
  *   delivered; 0 if absent;
  *
  * - the time-indexed variable "ReferenceSchedule", the active power
  *   \f$ \hat{p}_t \f$ the unit is asked to follow: if given, (11) is
  *   generated and the Objective (12) contains
  *   \f$ \sigma \sum_t \delta^{rs}_t \f$ instead of the term in "Cost", and
  *   set_kappa() refuses to resize the unit;
  *
  * - the time-indexed variables "MaxPrimaryPower" and "MaxSecondaryPower",
  *   the nonnegative bounds \f$ P^{pr}_t \f$ and \f$ P^{sc}_t \f$ of (9);
  *   the primary (respectively, secondary) reserve variables exist only if
  *   the datum is given and the enclosing UCBlock asks for that reserve
  *   [see UnitBlock::set_reserve_vars()];
  *
  * - the scalar variable "Kappa", the nonnegative \f$ \kappa \f$; 1 if
  *   absent;
  *
  * - the scalar variable "ConverterKappa", the nonnegative
  *   \f$ \kappa^{cv} \f$ of the converter in (10) and (10d); if absent,
  *   \f$ \kappa^{cv} = \kappa \f$ whatever \f$ \kappa \f$ is set to, i.e.,
  *   one multiplier resizes the whole unit, while if given the rows (10)
  *   exist also with no "ConverterMaxPower" and no design variable (with
  *   \f$ P^{cv}_t = P^{mx}_t \f$), and serialize() writes it;
  *
  * - the scalar variable "Scale", the scale factor \f$ \sigma \f$ [see
  *   UnitBlock::scale()]; 1 if absent. A value different from 1 is refused
  *   together with \f$ | \bar{X} | > 1 \f$ for the battery or the
  *   converter, since the two represent the same replication of the unit
  *   (by deserialize() only, scale() not repeating this check);
  *
  * - the scalar variables "BatteryInvestmentCost" and
  *   "ConverterInvestmentCost", the investment costs \f$ c^{inv}_b \f$ and
  *   \f$ c^{inv}_c \f$ of (12); 0 if absent, and a nonzero value creates the
  *   corresponding design variable;
  *
  * - only if the corresponding investment cost is given, the scalar
  *   variables "BatteryMinCapacityDesign" and "BatteryMaxCapacityDesign",
  *   respectively "ConverterMinCapacityDesign" and
  *   "ConverterMaxCapacityDesign", the bounds \f$ \underline{X} \geq 0 \f$
  *   (0 if absent) and \f$ \bar{X} \f$ (1 if absent) that give the domain of
  *   the design variable described in the class description; it must be
  *   \f$ \underline{X} \leq | \bar{X} | \f$;
  *
  * - the time-indexed variables "MinReactivePower" and "MaxReactivePower",
  *   the bounds of the reactive power, used only if the enclosing UCBlock
  *   asks for it [see UnitBlock::set_reactive_power()]; a vector of zeros is
  *   treated as absent.
  *
  * The scalar variables "BatteryMaxCapacity" and "ConverterMaxCapacity" may
  * be present, and are ignored. */

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG
 // extends UnitBlock::expected_dims()

 std::vector< std::string > expected_dims( void ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// extends UnitBlock::expected_vars()

 std::vector< std::string > expected_vars( void ) const override;

#endif

/*--------------------------------------------------------------------------*/
 /// generate the abstract variables of the BatteryUnitBlock
 /** Generates the static Variable of the BatteryUnitBlock, with the symbols
  * of the class description; each group is a std::vector< ColVariable > of
  * size get_time_horizon(), or empty if it is not generated:
  *
  * - the active power \f$ p^{ac}_t \f$ ("p_battery"), continuous and free;
  *
  * - the level \f$ v^{ba}_t \f$ ("sl_battery"), continuous and free (its
  *   sign comes from (7) only);
  *
  * - the intake \f$ p^-_t \f$ ("il_battery") and the outtake \f$ p^+_t \f$
  *   ("ol_battery"), nonnegative, only if "StoringBatteryRho" and
  *   "ExtractingBatteryRho" differ at some instant [see
  *   needs_intake_outtake()];
  *
  * - the binary variables \f$ u^{ch}_t \f$ ("b_battery"), 1 when charging,
  *   only if negative prices may occur (see below) and
  *   \f$ \rho^{in}_t < \rho^{out}_t \f$ at some instant, which implies that
  *   the intake and the outtake exist;
  *
  * - the primary and secondary reserves \f$ p^{pr}_t \f$ ("pr_battery") and
  *   \f$ p^{sc}_t \f$ ("sc_battery"), nonnegative, each only if the
  *   enclosing UCBlock asks for that reserve [see set_reserve_vars()] and
  *   the corresponding "MaxPrimaryPower" or "MaxSecondaryPower" is given;
  *
  * - the deviations \f$ \delta^{rs}_t \f$ ("v_absb_refschd") from the
  *   reference schedule, nonnegative, only if "ReferenceSchedule" is given;
  *
  * - the reactive power ("q_battery"), nonnegative, only if the enclosing
  *   UCBlock asks for it [see set_reactive_power()];
  *
  * and the single design variables \f$ x_b \f$ ("x_battery") and
  * \f$ x_c \f$ ("x_converter"), only if "BatteryInvestmentCost",
  * respectively "ConverterInvestmentCost", is nonzero: each is nonnegative
  * if the corresponding \f$ \bar{X} \geq 0 \f$, integer if \f$ \bar{X} < 0
  * \f$ (binary if \f$ \bar{X} = -1 \f$), its bounds being set by
  * generate_abstract_constraints().
  *
  * Whether negative prices may occur is told by \p stvv or, if it is
  * nullptr, by the f_static_variables_Configuration of the BlockConfig, if
  * any: a SimpleConfiguration< int > with a nonzero value means that they
  * may. By default they may not, and the binary variables are not
  * generated. */

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the static constraints of the BatteryUnitBlock
 /** Generates the static Constraint of the BatteryUnitBlock, i.e., the rows
  * (1)-(11) of the class description, with the forms (1d), (2d), (5d) and
  * (7d) when the battery has a design variable. Each group is a
  * std::vector or a boost::multi_array of size get_time_horizon() along the
  * time index, and it is generated only under the condition the class
  * description gives for the rows:
  *
  * - (1) and (2) are "ActivePower_Battery", a 2 x T multi_array of
  *   FRowConstraint whose first index is 0 for (2) and 1 for (1); with a
  *   battery design variable they are "ActivePower_Design_Battery", with
  *   the same layout, and "ActivePower_Battery" is empty;
  *
  * - (3) are "RampUp_Battery" and "RampDown_Battery", two vectors of
  *   FRowConstraint;
  *
  * - (4) is "PowerIntakeOuttake_Battery", a vector of FRowConstraint;
  *
  * - (5) is "IntakeOuttake_Battery", a 2 x T multi_array of BoxConstraint
  *   whose first index is 0 for the intake and 1 for the outtake, or a
  *   1 x T one on the active power when the split is not generated; with a
  *   battery design variable (5d) is in the rows 0 (intake) and 1 (outtake)
  *   of "IntakeOuttake_Design_Battery", a multi_array of FRowConstraint;
  *
  * - (6) is "Demand_Battery", a vector of FRowConstraint [see
  *   get_storage_balance_constraints()];
  *
  * - (7) is "StorageLevel_Battery", a vector of BoxConstraint, or, with a
  *   battery design variable, (7d) is "StorageLevel_Design_Battery", a
  *   2 x T multi_array of FRowConstraint whose first index is 0 for the
  *   lower and 1 for the upper bound;
  *
  * - (8) is "Intake_Outtake_Binary_Battery", a 2 x T multi_array of
  *   FRowConstraint whose first index is 0 for the intake and 1 for the
  *   outtake;
  *
  * - (9) are "Primary_UpperBound_Battery" and "Secondary_UpperBound_Battery",
  *   two vectors of LB0Constraint;
  *
  * - (10) is, with no battery design variable, "Converter_Battery", a
  *   multi_array of FRowConstraint with one row per instant on the intake
  *   plus the outtake, or two (on \f$ p^{ac}_t \f$ and on
  *   \f$ - p^{ac}_t \f$) when the split is not generated [see
  *   get_converter_bounds()]; with a battery design variable it is the row 2
  *   of "IntakeOuttake_Design_Battery", or its rows 2 and 3 on
  *   \f$ \pm p^{ac}_t \f$ when the split is not generated and (10d) holds
  *   (a converter design variable, or none and no "ConverterMaxPower"), or,
  *   when the split is not generated and the converter power is the
  *   constant "ConverterMaxPower", the bound on \f$ p^{ac}_t \f$ in
  *   "IntakeOuttake_Battery" (present if \f$ P^{cv}_t > 0 \f$ at some
  *   instant); at the instants where \f$ P^{cv}_t \leq 0 \f$ these rows
  *   are present and non-binding;
  *
  * - (11) is "Norm1B_Reference_Schedule", a vector of 2 T FRowConstraint,
  *   the upper sides first;
  *
  * - the bounds of the design variables are "BattDesignBound_Battery" and
  *   "ConvDesignBound_Battery", two BoxConstraint, present whenever the
  *   variable is;
  *
  * - (9d) are "Primary_Design_Battery" and "Secondary_Design_Battery", two
  *   vectors of T FRowConstraint, present when the battery design variable
  *   and the corresponding reserve exist (see
  *   get_primary_reserve_design_rows()); (8d) are the rows (8), with the
  *   coefficients given there;
  *
  * - the bounds of the reactive power are "ReactivePowerBound", a vector of
  *   BoxConstraint, if the reactive power exists and one of its bounds is
  *   given;
  *
  * - "Binary_Battery", a vector of ZOConstraint on the binary variables,
  *   only if \p stcc (or, if it is nullptr, the
  *   f_static_constraints_Configuration of the BlockConfig) is a
  *   SimpleConfiguration< int > with a nonzero value; these rows repeat the
  *   bounds of the variables, which some approaches need in order to have a
  *   dual value for them. */

 void generate_abstract_constraints( Configuration * stcc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the objective of the BatteryUnitBlock
 /** Generates the Objective (12) of the class description, to be minimized:
  * a LinearFunction with \f$ \sigma b_t \f$ on the intake and
  * \f$ - \sigma b_t \f$ on the outtake, or \f$ - \sigma b_t \f$ on the
  * active power when the split is not generated, unless "ReferenceSchedule"
  * is given, in which case these terms are absent and the deviations
  * \f$ \delta^{rs}_t \f$ have coefficient \f$ \sigma \f$; plus
  * \f$ \sigma c^{inv}_b \f$ on \f$ x_b \f$ and \f$ \sigma c^{inv}_c \f$ on
  * \f$ x_c \f$ if they exist. "Cost" is thus a price on the net injection,
  * not a cost of wear charged in both directions. */

 void generate_objective( Configuration * objc = nullptr ) override;

/**@} ----------------------------------------------------------------------*/
/*---------------- Methods for checking the BatteryUnitBlock ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for checking solution information in the BatteryUnitBlock
 * @{ */

 /// returns true if the current solution is (approximately) feasible
 /** This function returns true if and only if the solution encoded in the
  * current value of the Variable of this BatteryUnitBlock is approximately
  * feasible within the given tolerance. That is, a solution is considered
  * feasible if and only if
  *
  * -# each ColVariable is feasible; and
  *
  * -# the violation of each Constraint of this BatteryUnitBlock is not
  *    greater than the tolerance.
  *
  * Every Constraint of this BatteryUnitBlock is a RowConstraint and its
  * violation is given by either the relative (see RowConstraint::rel_viol())
  * or the absolute violation (see RowConstraint::abs_viol()), depending on
  * the Configuration that is provided.
  *
  * The tolerance and the type of violation can be provided by either \p fsbc
  * or f_BlockConfig->f_is_feasible_Configuration and they are determined as
  * follows:
  *
  * - If \p fsbc is not nullptr and it is a pointer to a
  *   SimpleConfiguration< double >, then the tolerance is the value present
  *   in that SimpleConfiguration and the relative violation is considered.
  *
  * - If \p fsbc is not nullptr and it is a pointer to a
  *   SimpleConfiguration< std::pair< double , int > >, then the tolerance is
  *   fsbc->f_value.first and the type of violation is determined by
  *   fsbc->f_value.second (any nonzero number for relative violation and
  *   zero for absolute violation);
  *
  * - Otherwise, if both f_BlockConfig and
  *   f_BlockConfig->f_is_feasible_Configuration are not nullptr and the
  *   latter is a pointer to either a SimpleConfiguration< double > or to a
  *   SimpleConfiguration< std::pair< double , int > >, then the values of the
  *   parameters are obtained analogously as above;
  *
  * - Otherwise, by default, the tolerance is 0 and the relative violation is
  *   considered.
  *
  * This function considers only the abstract representation to
  * determine if the solution is feasible. So, the parameter \p useabstract is
  * ignored. If no abstract Variable has been generated, then this
  * function returns true. Moreover, if no abstract Constraint has been
  * generated, the solution is considered to be feasible with respect to the
  * set of Variables only. Notice also that, before checking if the solution
  * satisfies a Constraint, the Constraint is computed (Constraint::compute()).
  *
  * @param useabstract This parameter is ignored.
  *
  * @param fsbc The pointer to a Configuration that specifies the tolerance
  *             and the type of violation that must be considered. */

 bool is_feasible( bool useabstract = false ,
                   Configuration * fsbc = nullptr ) override;


/** @} ---------------------------------------------------------------------*/
/*--------- METHODS FOR READING THE DATA OF THE BatteryUnitBlock -----------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the BatteryUnitBlock
 *
 * These methods allow to read data that must be common to (in principle) all
 * kinds of battery storage units
 * @{ */

 /// returns the initial storage value
 double get_initial_storage( void ) const { return( f_InitialStorage ); }

 /// returns the initial power value
 double get_initial_power( void ) const { return( f_InitialPower ); }

 /// returns the maximum C-rate of the battery in charge
 double get_max_C_rate_charge( void ) const { return( f_MaxCRateCharge ); }

 /// returns the maximum C-rate of the battery in discharge
 double get_max_C_rate_discharge( void ) const {
  return( f_MaxCRateDischarge );
  }

 /// returns the battery investment cost
 double get_batt_investment_cost( void ) const {
  return( f_BattInvestmentCost );
  }

 /// returns the converter investment cost
 double get_conv_investment_cost( void ) const {
  return( f_ConvInvestmentCost );
  }

 /// returns the maximum battery installable capacity by the user
 double get_batt_max_capacity_design( void ) const {
  return( f_BattMaxCapacityDesign );
  }

 /// returns the maximum converter installable capacity by the user
 double get_conv_max_capacity_design( void ) const {
  return( f_ConvMaxCapacityDesign );
  }

/*--------------------------------------------------------------------------*/
 /// returns the minimum levels ("MinStorage")
 /** Returns the vector of the minimum levels \f$ V^{mn}_t \f$ of (7), of
  * size get_time_horizon(); the methods below that return the vector of a
  * time-indexed datum likewise return it either empty (the datum is absent)
  * or with one value per instant [see deserialize()]. */

 const std::vector< double > & get_min_storage( void ) const {
  return( v_MinStorage );
  }

/*--------------------------------------------------------------------------*/
 /// returns the maximum levels ("MaxStorage")

 const std::vector< double > & get_max_storage( void ) const {
  return( v_MaxStorage );
  }

/*--------------------------------------------------------------------------*/
 /// returns the minimum active power ("MinPower") at instant t

 double get_min_power( Index t , Index generator = 0 ) const override {
  return( v_MinPower[ t ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the maximum active power ("MaxPower") at instant t

 double get_max_power( Index t , Index generator = 0 ) const override {
  return( v_MaxPower[ t ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the powers of the converter ("ConverterMaxPower")
 /** Returns the vector of the powers \f$ P^{cv}_t \f$ of (10), which is
  * filled with "MaxPower" when "ConverterMaxPower" is absent (the rows (10)
  * of a battery with no design variable being then absent). */

 const std::vector< double > & get_converter_max_power( void ) const {
  return( v_ConvMaxPower );
  }

/*--------------------------------------------------------------------------*/
 /// returns the minimum reactive power of the \p generator at time \p t

 double get_min_reactive_power( Index t , Index generator = 0 )
  const override {
  return( ( v_MinReactivePower.size() > t ) ? v_MinReactivePower[ t ] : 0 );
  }

/*--------------------------------------------------------------------------*/
 /// returns the maximum reactive power of the \p generator at time \p t

 double get_max_reactive_power( Index t , Index generator = 0 )
  const override {
  return( ( v_MaxReactivePower.size() > t ) ? v_MaxReactivePower[ t ] : 0 );
  }

/*--------------------------------------------------------------------------*/
 /// returns the bounds of the primary reserve ("MaxPrimaryPower")
 /** Returns the vector of the bounds \f$ P^{pr}_t \f$ of (9), empty if the
  * unit gives no primary reserve. */

 const std::vector< double > & get_max_primary_power( void ) const {
  return( v_MaxPrimaryPower );
  }

/*--------------------------------------------------------------------------*/
 /// returns the bounds of the secondary reserve ("MaxSecondaryPower")
 /** Returns the vector of the bounds \f$ P^{sc}_t \f$ of (9), empty if the
  * unit gives no secondary reserve. */

 const std::vector< double > & get_max_secondary_power( void ) const {
  return( v_MaxSecondaryPower );
  }

/*--------------------------------------------------------------------------*/
 /// returns the reference schedule of the unit, if it has one
 /** Returns the vector of the active power \f$ \hat{p}_t \f$ the unit is
  * asked to follow ("ReferenceSchedule"), empty if it has none; the
  * deviation from it, weighed with the scale factor, replaces the term in
  * "Cost" in the Objective [see generate_objective()]. */

 const std::vector< double > & get_reference_schedule( void ) const {
  return( v_RefSchedule );
  }

/*--------------------------------------------------------------------------*/
 /// returns the ramp-up limits ("DeltaRampUp")
 /** Returns the vector of the limits \f$ \Delta^+_t \f$ of (3), empty if
  * there is no ramp-up row. */

 const std::vector< double > & get_delta_ramp_up( void ) const {
  return( v_DeltaRampUp );
  }

/*--------------------------------------------------------------------------*/
 /// returns the ramp-down limits ("DeltaRampDown")
 /** Returns the vector of the limits \f$ \Delta^-_t \f$ of (3), empty if
  * there is no ramp-down row. */

 const std::vector< double > & get_delta_ramp_down( void ) const {
  return( v_DeltaRampDown );
  }

/*--------------------------------------------------------------------------*/
 /// returns the coefficients of the intake ("StoringBatteryRho")
 /** Returns the vector of the coefficients \f$ \rho^{in}_t \f$ of (6),
  * empty if the datum is absent (all of them being then 1). */

 const std::vector< double > & get_storing_battery_rho( void ) const {
  return( v_StoringBatteryRho );
  }

/*--------------------------------------------------------------------------*/
 /// returns the coefficients of the outtake ("ExtractingBatteryRho")
 /** Returns the vector of the coefficients \f$ \rho^{out}_t \f$ of (6),
  * empty if the datum is absent (all of them being then 1). */

 const std::vector< double > & get_extracting_battery_rho( void ) const {
  return( v_ExtractingBatteryRho );
  }

/*--------------------------------------------------------------------------*/
 /// returns the coefficients of the level ("StandingBatteryRho")
 /** Returns the vector of the coefficients \f$ \rho^{st}_t \f$ of (6),
  * empty if the datum is absent (all of them being then 1). */

 const std::vector< double > & get_standing_battery_rho( void ) const {
  return( v_StandingBatteryRho );
  }

/*--------------------------------------------------------------------------*/
 /// returns the prices on the net injection ("Cost")
 /** Returns the vector of the prices \f$ b_t \f$ of (12), of size
  * get_time_horizon() (zeros if "Cost" is absent): \f$ b_t \f$ is paid for
  * the energy drawn and earned for the energy delivered at instant
  * \f$ t \f$. */

 const std::vector< double > & get_cost( void ) const {
  return( v_Cost );
  }

/*--------------------------------------------------------------------------*/
 /// returns the energies withdrawn for another use ("Demand")
 /** Returns the vector of the energies \f$ d^{ba}_t \f$ of (6), empty if
  * the datum is absent (all of them being then 0). */

 const std::vector< double > & get_demand( void ) const {
  return( v_Demand );
  }

/*--------------------------------------------------------------------------*/
 /// returns the scale factor of this BatteryUnitBlock

 double get_scale( void ) const override { return( f_scale ); }

 double get_design_ub( void ) const override {
  // the active power is bounded by (1d)-(2d), which depend on x_b only (a
  // converter design bounds the power further, never beyond them): the
  // bound is that of x_b, and 1 without a battery design
  return( batt_design_mult() );
 }

/**@} ----------------------------------------------------------------------*/
/*-------- METHODS FOR READING THE Variable OF THE BatteryUnitBlock --------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the Variable of the BatteryUnitBlock
 *
 * These methods give the Variable of the BatteryUnitBlock [see
 * generate_abstract_variables()], each group a std::vector< ColVariable >
 * of size get_time_horizon() or empty if it is not generated, and the
 * groups of Constraint used from outside the Block.
 * @{ */

 /// returns the capacity multiplier kappa
 /** Returns \f$ \kappa \f$ ("Kappa"), which multiplies \f$ P^{mn}_t \f$,
  * \f$ P^{mx}_t \f$, \f$ P^{pr}_t \f$, \f$ P^{sc}_t \f$, \f$ V^{mn}_t \f$
  * and \f$ V^{mx}_t \f$ in the rows of the class description, also
  * \f$ P^{cv}_t \f$ if the converter follows it [see
  * converter_follows_kappa()], and neither the ramp limits nor \f$ V^0 \f$
  * nor \f$ d^{ba}_t \f$. */

 double get_kappa( void ) const override { return( f_kappa ); }

/*--------------------------------------------------------------------------*/
 /// returns the capacity multiplier of the converter
 /** Returns \f$ \kappa^{cv} \f$, which multiplies \f$ P^{cv}_t \f$ in (10)
  * and (10d): "ConverterKappa" if it has been given (in the data or by
  * set_converter_kappa()), and \f$ \kappa \f$ otherwise. */

 double get_converter_kappa( void ) const {
  return( converter_follows_kappa() ? f_kappa : f_conv_kappa );
  }

/*--------------------------------------------------------------------------*/
 /// tells whether the converter is resized by kappa
 /** Returns true if \f$ \kappa^{cv} \f$ has not been given, neither in the
  * data ("ConverterKappa") nor by set_converter_kappa(), in which case
  * \f$ \kappa^{cv} = \kappa \f$ and set_kappa() resizes storage and
  * converter together. */

 bool converter_follows_kappa( void ) const {
  return( std::isnan( f_conv_kappa ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the linearization coefficient of the kappa-parametrized objective
 /** When no design variable exists, \f$ \kappa \f$ appears in right-hand
  * sides and, in the binary rows (8), as a coefficient of
  * \f$ u^{ch}_t \f$; in the continuous relaxation, the change of variable
  * \f$ w_t = \kappa u^{ch}_t \f$ makes all the rows linear in
  * \f$ \kappa \f$ and the other variables jointly (see InvestmentFunction),
  * so that the optimal value of a convex minimization problem containing
  * this unit (e.g., its continuous relaxation) is a convex function of
  * \f$ \kappa \f$, and a subgradient of it is
  * \f[
  *   \sum_{ t \in \mathcal{T} } \Bigl( P^{mn}_t \bigl( \mu^{(2)}_t +
  *   C^{ch} \mu^{(5),-}_t + u^{ch}_t \mu^{(8),-}_t \bigr) - P^{mx}_t
  *   \bigl( \mu^{(1)}_t + C^{dis} \mu^{(5),+}_t + ( 1 - u^{ch}_t )
  *   \mu^{(8),+}_t \bigr) + V^{mn}_t \mu^{(7),mn}_t - V^{mx}_t
  *   \mu^{(7),mx}_t - P^{pr}_t \mu^{(9),pr}_t - P^{sc}_t \mu^{(9),sc}_t
  *   \Bigr) + \beta \, g^{cv} ,
  * \f]
  * with \f$ g^{cv} \f$ the linearization of the converter [see
  * get_converter_kappa_linearization()] and \f$ \beta = 1 \f$ if the
  * converter follows \f$ \kappa \f$ [see converter_follows_kappa()], since
  * then \f$ \kappa^{cv} = \kappa \f$ and the method gives the derivative
  * with respect to the one multiplier of the whole unit, and
  * \f$ \beta = 0 \f$ otherwise, since then the method gives the derivative
  * with respect to the storage alone,
  * where every \f$ \mu \geq 0 \f$ is the absolute value of the dual value
  * of the row of the class description in its superscript, on the side
  * (intake or charge, outtake or discharge, lower or upper bound) that the
  * further superscript says; the side of a two-sided row the dual value
  * belongs to is told by its sign, and a term exists only if its row
  * does. This method returns that value
  * for the dual solution held by the Constraint, with two caveats.
  * First, the terms of the binary rows (8) are evaluated at the current
  * value of \f$ u^{ch}_t \f$, so that the value is a subgradient only when
  * these variables are fixed or relaxed. Second, where \f$ \kappa = 0 \f$
  * the two sides of a pair of rows meet, both dual values may be nonzero
  * on what is one equality, and of the two possible readings the one that
  * keeps the linearization below the value function (the smaller one in a
  * minimization) is taken. If the battery or the converter has a design
  * variable the data multiply it instead, and the method throws
  * std::logic_error. */

 double get_kappa_linearization( void ) const override;

/*--------------------------------------------------------------------------*/
 /// returns the linearization coefficient with respect to the converter
 /** When no design variable exists, \f$ \kappa^{cv} \f$ appears only in
  * the right-hand sides of the rows (10), so that the optimal value of a
  * convex minimization problem containing this unit is a convex function
  * of \f$ \kappa^{cv} \f$ (the other data fixed), and a subgradient of it
  * is
  * \f[
  *   g^{cv} = - \sum_{ t \in \mathcal{T} } P^{cv}_t \mu^{(10)}_t ,
  * \f]
  * where \f$ \mu^{(10)}_t \geq 0 \f$ is the sum of the absolute values of
  * the dual values of the rows (10) of instant \f$ t \f$ (one row with the
  * split, two without it), 0 where \f$ P^{cv}_t \leq 0 \f$ or if the rows
  * do not exist. Where \f$ \kappa^{cv} = 0 \f$ the two folded rows meet
  * on what is one equality, and their net dual value is taken, as in
  * get_kappa_linearization(). The value is the derivative with respect to
  * \f$ \kappa^{cv} \f$ also when the converter follows \f$ \kappa \f$, in
  * which case it is the part of get_kappa_linearization() due to the
  * converter. If the battery or the converter has a design variable the
  * method throws std::logic_error. */

 double get_converter_kappa_linearization( void ) const;

/*--------------------------------------------------------------------------*/
 /// returns the vector of storage level variables
 /** This method returns a vector V containing the storage level
  * variables. There are two possible cases:
  *
  * - if V is empty(), then these variables are not defined;
  *
  * - otherwise, V must have size get_time_horizon() and V[ t ] is the storage
  *   level variable for time step t. */

 std::vector< ColVariable > & get_storage_level( void ) {
  return( v_storage_level );
  }

/*--------------------------------------------------------------------------*/
 /// returns the const vector of storage level variables
 /** This method returns a const vector V containing the storage level
  * variables. There are two possible cases:
  *
  * - if V is empty(), then these variables are not defined;
  *
  * - otherwise, V must have size get_time_horizon() and V[ t ] is the storage
  *   level variable for time step t. */

 const std::vector< ColVariable > & get_const_storage_level( void ) const {
  return( v_storage_level );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of intake level variables
 /** This method returns a vector V containing the intake level variables.
  * There are two possible cases:
  *
  * - if V is empty(), then these variables are not defined;
  *
  * - otherwise, V must have size get_time_horizon() and V[ t ] is the intake
  *   level variable for time step t. */

 std::vector< ColVariable > & get_intake_level( void ) {
  return( v_intake_level );
  }

/*--------------------------------------------------------------------------*/
 /// returns the const vector of intake level variables
 /** This method returns a const vector V containing the intake level
  * variables. There are two possible cases:
  *
  * - if V is empty(), then these variables are not defined;
  *
  * - otherwise, V must have size get_time_horizon() and V[ t ] is the intake
  *   level variable for time step t. */

 const std::vector< ColVariable > & get_const_intake_level( void ) const {
  return( v_intake_level );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of outtake level variables
 /** This method returns a vector V containing the outtake level variables.
  * There are two possible cases:
  *
  * - if V is empty(), then these variables are not defined;
  *
  * - otherwise, V must have size get_time_horizon() and V[ t ] is the outtake
  *   level variable for time step t. */

 std::vector< ColVariable > & get_outtake_level( void ) {
  return( v_outtake_level );
  }

/*--------------------------------------------------------------------------*/
 /// returns the const vector of outtake level variables
 /** This method returns a const vector V containing the outtake level
  * variables. There are two possible cases:
  *
  * - if V is empty(), then these variables are not defined;
  *
  * - otherwise, V must have size get_time_horizon() and V[ t ] is the outtake
  *   level variable for time step t. */

 const std::vector< ColVariable > & get_const_outtake_level( void ) const {
  return( v_outtake_level );
  }

/*--------------------------------------------------------------------------*/
 /// a battery has one storage, its charge

 Index get_number_storages( void ) const override { return( 1 ); }

/*--------------------------------------------------------------------------*/
 /// returns the array of storage level variables [see get_storage_level()]

 ColVariable * get_storage_level( Index storage ) override {
  if( ( storage > 0 ) || v_storage_level.empty() )
   return( nullptr );
  return( v_storage_level.data() );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of active power variables

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
 /// returns the vector of primary spinning reserve variables

 /// the unit provides the reserve only if it has any to give
 /** The enclosing UCBlock asking for the reserve is not enough, the unit has
  * to have some to give: this is the very condition with which the Variable
  * are generated, said in terms of the data alone. */

 bool has_primary_reserve( void ) const override {
  return( ( reserve_vars & 1u ) && ( ! v_MaxPrimaryPower.empty() ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 bool has_secondary_reserve( void ) const override {
  return( ( reserve_vars & 2u ) && ( ! v_MaxSecondaryPower.empty() ) );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 ColVariable * get_primary_spinning_reserve( Index generator ) override {
  if( v_primary_spinning_reserve.empty() )
   return( nullptr );
  return( &( v_primary_spinning_reserve.front() ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of secondary spinning reserve variables

 ColVariable * get_secondary_spinning_reserve( Index generator ) override {
  if( v_secondary_spinning_reserve.empty() )
   return( nullptr );
  return( &( v_secondary_spinning_reserve.front() ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the battery design variable

 ColVariable & get_batt_design( void ) { return( batt_design ); }

/*--------------------------------------------------------------------------*/
 /// returns the const battery design variable

 const ColVariable & get_const_batt_design( void ) const {
  return( batt_design );
  }

/*--------------------------------------------------------------------------*/
 /// returns the converter design variable

 ColVariable & get_conv_design( void ) { return( conv_design ); }

/*--------------------------------------------------------------------------*/
 /// returns the const converter design variable

 const ColVariable & get_const_conv_design( void ) const {
  return( conv_design );
  }

/*--------------------------------------------------------------------------*/
 /// returns the intake/outtake binary variables

 const std::vector< ColVariable > &
 get_intake_outtake_binary_variables( void ) const {
  return( v_battery_binary );
  }

/*--------------------------------------------------------------------------*/
 /// returns the minimum power output constraints
 const FRowConstraint * get_min_power_constraints( void ) const {
  if( active_power_bounds_Const.empty() ||
      active_power_bounds_Const[ 0 ].empty() )
   return( nullptr );
  return( &( active_power_bounds_Const[ 0 ][ 0 ] ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the minimum power output constraint associated with time t
 const FRowConstraint * get_min_power_constraint( Index t ) const {
  if( active_power_bounds_Const.empty() ||
      active_power_bounds_Const[ 0 ].empty() )
   return( nullptr );
  return( &( active_power_bounds_Const[ 0 ][ t ] ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the maximum power output constraints
 const FRowConstraint * get_max_power_constraints( void ) const {
  if( active_power_bounds_Const.empty() ||
      active_power_bounds_Const[ 1 ].empty() )
   return( nullptr );
  return( &( active_power_bounds_Const[ 1 ][ 0 ] ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the maximum power output constraint associated with time t
 const FRowConstraint * get_max_power_constraint( Index t ) const {
  if( active_power_bounds_Const.empty() ||
      active_power_bounds_Const[ 1 ].empty() )
   return( nullptr );
  return( &( active_power_bounds_Const[ 1 ][ t ] ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the intake upper bound constraints with binary variables
 const FRowConstraint * get_max_intake_binary_constraints( void ) const {
  if( intake_outtake_binary_Const.empty() ||
      intake_outtake_binary_Const[ 0 ].empty() )
   return( nullptr );
  return( &( intake_outtake_binary_Const[ 0 ][ 0 ] ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the intake upper bound constraint with binary variables for time t
 const FRowConstraint * get_max_intake_binary_constraint( Index t ) const {
  if( intake_outtake_binary_Const.empty() ||
      intake_outtake_binary_Const[ 0 ].empty() )
   return( nullptr );
  return( &( intake_outtake_binary_Const[ 0 ][ t ] ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the outtake upper bound constraints with binary variables
 const FRowConstraint * get_max_outtake_binary_constraints( void ) const {
  if( intake_outtake_binary_Const.empty() ||
      intake_outtake_binary_Const[ 1 ].empty() )
   return( nullptr );
  return( &( intake_outtake_binary_Const[ 1 ][ 0 ] ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the outtake upper bound constraint + binary variables for time t
 const FRowConstraint * get_max_outtake_binary_constraint( Index t ) const {
  if( intake_outtake_binary_Const.empty() ||
      intake_outtake_binary_Const[ 1 ].empty() )
   return( nullptr );
  return( &( intake_outtake_binary_Const[ 1 ][ t ] ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the storage level bound constraints
 const std::vector< BoxConstraint > & get_storage_level_bounds( void ) const {
  return( storage_level_bounds_Const );
  }

/*--------------------------------------------------------------------------*/
 /// returns the intake upper bound constraints
 /** The first row of the intake/outtake bound Constraints, one per time
  * instant. When the unit keeps intake and outtake apart this bounds the
  * intake level; when it does not there is a single row, one two-sided bound
  * on the active power, and this is it: which of the two it is is told by
  * get_max_outtake_bounds() being nullptr. */

 const BoxConstraint * get_max_intake_bounds( void ) const {
  if( intake_outtake_bounds_Const.num_elements() == 0 )
   return( nullptr );
  return( &( intake_outtake_bounds_Const[ 0 ][ 0 ] ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the intake upper bound constraint associated with time t
 const BoxConstraint * get_max_intake_bound( Index t ) const {
  if( intake_outtake_bounds_Const.num_elements() == 0 )
   return( nullptr );
  return( &( intake_outtake_bounds_Const[ 0 ][ t ] ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the outtake upper bound constraints
 /** The second row of the intake/outtake bound Constraints, one per time
  * instant, or nullptr when the unit keeps a single one [see
  * get_max_intake_bounds()]. */

 const BoxConstraint * get_max_outtake_bounds( void ) const {
  if( intake_outtake_bounds_Const.shape()[ 0 ] < 2 )
   return( nullptr );
  return( &( intake_outtake_bounds_Const[ 1 ][ 0 ] ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the outtake upper bound constraint associated with time t
 const BoxConstraint * get_max_outtake_bound( Index t ) const {
  if( intake_outtake_bounds_Const.shape()[ 0 ] < 2 )
   return( nullptr );
  return( &( intake_outtake_bounds_Const[ 1 ][ t ] ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the primary reserve bound constraints
 const std::vector< LB0Constraint > & get_primary_reserve_bounds( void )
  const {
  return( primary_upper_bound_Const );
  }

/*--------------------------------------------------------------------------*/
 /// returns the secondary reserve bound constraints
 const std::vector< LB0Constraint > & get_secondary_reserve_bounds( void )
  const {
  return( secondary_upper_bound_Const );
  }

/*--------------------------------------------------------------------------*/
 /// returns the primary reserve rows of the design
 /** The rows \f$ p^{pr}_t \leq \kappa P^{pr}_t x_b \f$, one per instant,
  * present if the battery design variable and the primary reserve exist,
  * empty otherwise. */

 const std::vector< FRowConstraint > & get_primary_reserve_design_rows( void )
  const {
  return( primary_upper_bound_design_Const );
  }

/*--------------------------------------------------------------------------*/
 /// returns the secondary reserve rows of the design
 /** The rows \f$ p^{sc}_t \leq \kappa P^{sc}_t x_b \f$, as
  * get_primary_reserve_design_rows(). */

 const std::vector< FRowConstraint > &
  get_secondary_reserve_design_rows( void ) const {
  return( secondary_upper_bound_design_Const );
  }

/*--------------------------------------------------------------------------*/
 /// returns the storage balance constraints (6)
 /** One row per instant, empty if the Constraints are not generated [see
  * generate_abstract_constraints()]. */

 const std::vector< FRowConstraint > & get_storage_balance_constraints( void )
  const {
  return( demand_Const );
  }

/*--------------------------------------------------------------------------*/
 /// returns the converter bound constraints of a battery with no design
 /** The rows (10) of a battery with no design Variable, empty if there are
  * none: the first index is 0 for the row on intake plus outtake, and 0 and
  * 1 for the two rows on the active power when the pair is folded onto it
  * [see generate_abstract_constraints()]. */

 const boost::multi_array< FRowConstraint , 2 > &
 get_converter_bounds( void ) const {
  return( converter_bounds_Const );
  }

/*--------------------------------------------------------------------------*/
 /// returns the intake and outtake rows of a battery with a design variable
 /** The rows "IntakeOuttake_Design_Battery" of a battery with a design
  * Variable, empty if there is none: the first index is 0 and 1 for the
  * two rows (5d), and 2 (and 3, on \f$ \pm p^{ac}_t \f$, when the pair is
  * folded onto the active power) for the rows (10d) or for the bound (10)
  * with the constant converter power [see generate_abstract_constraints()].
  */

 const boost::multi_array< FRowConstraint , 2 > &
 get_intake_outtake_design_rows( void ) const {
  return( intake_outtake_upper_bounds_design_Const );
  }

/** @} ---------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Solution
 * @{ */

 /// returns a Solution storing the current solution of this BatteryUnitBlock
 /** This method must construct and return a (pointer to a) Solution object
  * representing the current "solution state" of this BatteryUnitBlock. This
  * is a BatteryUnitBlockSolution extending UnitBlockSolution with the
  * specific extra solution information of BatteryUnitBlock.
  *
  * The parameter for deciding which kind of Solution must be returned is a
  * single int value, coded bitwise:
  *
  * - the first four bits (bit 0 to bit 3) are "taken" by the base
  *   UnitBlock[Solution]
  *
  * - bit 4 (& 16) means "store the storage levels"
  *
  * - bit 5 (& 32) means "store the intakes and outtakes"
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
  * - otherwise, it is 63 (save everything). */

 Solution * get_Solution( Configuration * solc = nullptr ,
                          bool emptys = true ) override;

/*--------------------------------------------------------------------------*/
 /// return the "appropriate" [Battery]UnitBlockSolution

 UnitBlockSolution * new_Solution( void ) const override;

/** @} ---------------------------------------------------------------------*/
/*---------------- METHODS FOR SAVING THE BatteryUnitBlock------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the BatteryUnitBlock
 * @{ */

 /// extends Block::serialize( netCDF::NcGroup )
 /** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
  * BatteryUnitBlock. See BatteryUnitBlock::deserialize( netCDF::NcGroup ) for
  * details of the format of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/** @} ---------------------------------------------------------------------*/
/*--------------- METHODS FOR INITIALIZING THE BatteryUnitBlock ------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the data of the BatteryUnitBlock
 * @{ */

 void load( std::istream & input , char frmt = 0 ) override {
  throw( std::logic_error( "BatteryUnitBlock::load not implemented yet" ) );
  }

/** @} ---------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for changing the data of the BatteryUnitBlock
 *  @{ */

 /// sets the initial storage ("InitialStorage")
 /** Sets "InitialStorage" to the value associated with the index 0 in
  * \p subset (the last one, if 0 appears more than once; nothing is done if
  * it does not appear), and changes the right-hand side of (6) at
  * \f$ t = 0 \f$ accordingly. This is how a level is passed from one stage
  * to the next of a multistage problem. The sign of the datum decides
  * whether that row is cyclic (see the class description): after the
  * Constraints are generated a value of the other sign than the current
  * one cannot be represented, and std::logic_error is thrown. */

 void set_initial_storage( MF_dbl_it it ,
                           Subset && subset ,
                           const bool ordered = false ,
                           c_ModParam issuePMod = eNoBlck ,
                           c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// sets the initial storage ("InitialStorage")
 /** As the Subset version; the value is the element of the vector pointed
  * by \p it in the position opposite to \p rng.first (which is at most
  * 0 when \p rng contains 0), and nothing is done if \p rng does not
  * contain 0. */

 void set_initial_storage( MF_dbl_it it ,
                           Range rng = Range( 0 , Inf< Index >() ) ,
                           c_ModParam issuePMod = eNoBlck ,
                           c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// sets the initial power
 /** If the given \p subset contains the 0 index, this function sets the
  * initial power. If the given \p subset does not contain the index 0, this
  * function does nothing. Since \p subset can have multiple zeros, only the
  * last one is considered, which means that the value for the initial power
  * will be that in the vector pointed by \p it associated with this last
  * zero. */

 void set_initial_power( MF_dbl_it it ,
                         Subset && subset ,
                         const bool ordered = false ,
                         c_ModParam issuePMod = eNoBlck ,
                         c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// sets the initial power
 /** If the given Range \p rng contains 0, this function sets the initial
  * power. In this case, if the first element of \p rng is 0, the initial
  * power will be set to the value pointed by the given iterator. In general,
  * the initial power will be the one found at position -rng.first in the
  * vector pointed by \p it if this Range contains the 0 index. If the given
  * Range \p rng does not contain the 0 index, this function does nothing. */

 void set_initial_power( MF_dbl_it it ,
                         Range rng = Range( 0 , Inf< Index >() ) ,
                         c_ModParam issuePMod = eNoBlck ,
                         c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// set the cost values
 /** This function sets the prices \f$ b_t \f$ ("Cost") of this
  * BatteryUnitBlock, and the Objective with them [see generate_objective()].
  *
  * @param values  Iterator to a vector containing the cost values.
  * @param subset  If non-empty, the cost value of the instant
  *                \p subset[ i ] is set to the i-th value pointed by
  *                \p values, in the order \p subset gives, which need not be
  *                sorted. If empty, no operation is performed.
  * @param ordered It indicates whether \p subset is ordered.
  * @param issuePMod Controls how physical Modifications are issued.
  * @param issueAMod Controls how abstract Modifications are issued.
  */
 void set_cost( MF_dbl_it values ,
                Subset && subset ,
                bool ordered = false ,
                c_ModParam issuePMod = eNoBlck ,
                c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// set the cost values
 /** This function sets the cost values of this BatteryUnitBlock.
  *
  * @param values Iterator to a vector containing the cost values.
  * @param rng    If non-empty, the cost values corresponding to the indices
  *               in \p rng are set to the values pointed by \p values. If
  *               empty, no operation is performed.
  * @param issuePMod Controls how physical Modifications are issued.
  * @param issueAMod Controls how abstract Modifications are issued.
  */
 void set_cost( MF_dbl_it values ,
                Range rng = Range( 0 , Inf< Index >() ) ,
                c_ModParam issuePMod = eNoBlck ,
                c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// set the same cost at every instant
 /** This function sets the cost of this BatteryUnitBlock to \p value at
  * every instant of the time horizon.
  *
  * @param value     The value of the cost.
  * @param issuePMod Controls how physical Modifications are issued.
  * @param issueAMod Controls how abstract Modifications are issued.
  */
 void set_cost( double value ,
                c_ModParam issuePMod = eNoBlck ,
                c_ModParam issueAMod = eNoBlck ) {
  std::vector< double > vector( get_time_horizon() , value );
  set_cost( vector.cbegin() ,
            Range( 0 , Inf< Index >() ) ,
            issuePMod , issueAMod );
 }

/*--------------------------------------------------------------------------*/
 /// sets the scale factor of this BatteryUnitBlock
 /** This method sets the scale factor of this BatteryUnitBlock.
  *
  * @param values An iterator to a vector containing the scale factor.
  *
  * @param subset If non-empty, the scale factor is set to the value pointed
  *               by \p values. If empty, no operation is performed.
  *
  * @param ordered This parameter is ignored.
  *
  * @param issuePMod Controls how physical Modifications are issued.
  *
  * @param issueAMod Controls how abstract Modifications are issued. */

 void scale( MF_dbl_it values ,
             Subset && subset ,
             const bool ordered = false ,
             c_ModParam issuePMod = eNoBlck ,
             c_ModParam issueAMod = eNoBlck ) override;

/*--------------------------------------------------------------------------*/
 /// set the kappa constant
 /** This function sets \f$ \kappa \f$, which multiplies \f$ P^{mn}_t \f$,
  * \f$ P^{mx}_t \f$, \f$ P^{pr}_t \f$, \f$ P^{sc}_t \f$, \f$ V^{mn}_t \f$
  * and \f$ V^{mx}_t \f$ in the Constraints of this BatteryUnitBlock, and
  * also \f$ P^{cv}_t \f$ if the converter follows it [see get_kappa()]. A
  * unit that follows a reference schedule refuses to be resized
  * (std::logic_error is thrown), since whether the schedule of a resized
  * unit is the same one in absolute terms or one resized with it is not
  * decided by the data.
  *
  * @param values An iterator to a vector containing the kappa constants.
  *
  * @param subset If non-empty, the kappa constant is set to the value pointed
  *        by \p values. If empty, no operation is performed.
  *
  * @param ordered It indicates whether \p subset is ordered.
  *
  * @param issuePMod Controls how physical Modifications are issued.
  *
  * @param issueAMod Controls how abstract Modifications are issued. */

 void set_kappa( MF_dbl_it values ,
                 Subset && subset ,
                 const bool ordered = false ,
                 c_ModParam issuePMod = eNoBlck ,
                 c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// set the kappa constant
 /** This function sets \f$ \kappa \f$, which multiplies \f$ P^{mn}_t \f$,
  * \f$ P^{mx}_t \f$, \f$ P^{pr}_t \f$, \f$ P^{sc}_t \f$, \f$ V^{mn}_t \f$
  * and \f$ V^{mx}_t \f$ in the Constraints of this BatteryUnitBlock, and
  * also \f$ P^{cv}_t \f$ if the converter follows it [see get_kappa()]. A
  * unit that follows a reference schedule refuses to be resized
  * (std::logic_error is thrown), since whether the schedule of a resized
  * unit is the same one in absolute terms or one resized with it is not
  * decided by the data.
  *
  * @param values An iterator to a vector containing the kappa constants.
  *
  * @param rng If non-empty, the kappa constant is set to the value pointed by
  *        \p values. If empty, no operation is performed.
  *
  * @param issuePMod Controls how physical Modifications are issued.
  *
  * @param issueAMod Controls how abstract Modifications are issued. */

 void set_kappa( MF_dbl_it values ,
                 Range rng = Range( 0 , Inf< Index >() ) ,
                 c_ModParam issuePMod = eNoBlck ,
                 c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// set the kappa constant
 /** This function sets \f$ \kappa \f$, see the other versions.
  *
  * @param value The value of the kappa constant.
  *
  * @param issuePMod Controls how physical Modifications are issued.
  *
  * @param issueAMod Controls how abstract Modifications are issued. */

 void set_kappa( double value , c_ModParam issuePMod = eNoBlck ,
                 c_ModParam issueAMod = eNoBlck ) {
  std::vector< double > vector = { value };
  set_kappa( vector.cbegin() , Range( 0 , Inf< Index >() ) ,
             issuePMod , issueAMod );
  }

/*--------------------------------------------------------------------------*/
 /// set the capacity multiplier of the converter
 /** This function sets \f$ \kappa^{cv} \f$, which multiplies
  * \f$ P^{cv}_t \f$ in the rows (10) or (10d) [see get_converter_kappa()];
  * from then on the converter no longer follows \f$ \kappa \f$, i.e.,
  * set_kappa() resizes the storage alone, and the storage and the converter
  * are two assets of the unit. Called before the Constraints are generated,
  * it makes the rows (10) exist also with no "ConverterMaxPower" and no
  * design variable (with \f$ P^{cv}_t = P^{mx}_t \f$); called after, it
  * throws std::logic_error if the rows (10) do not exist, since a static
  * Constraint cannot be added then. As set_kappa(), it throws
  * std::logic_error for a unit that follows a reference schedule, and
  * std::invalid_argument for a negative value.
  *
  * @param values An iterator to a vector containing the multiplier.
  *
  * @param subset If non-empty, the multiplier is set to the value pointed
  *        by \p values. If empty, no operation is performed.
  *
  * @param ordered It indicates whether \p subset is ordered.
  *
  * @param issuePMod Controls how physical Modifications are issued.
  *
  * @param issueAMod Controls how abstract Modifications are issued. */

 void set_converter_kappa( MF_dbl_it values ,
                           Subset && subset ,
                           const bool ordered = false ,
                           c_ModParam issuePMod = eNoBlck ,
                           c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// set the capacity multiplier of the converter
 /** This function sets \f$ \kappa^{cv} \f$, see the other versions.
  *
  * @param values An iterator to a vector containing the multiplier.
  *
  * @param rng If non-empty, the multiplier is set to the value pointed by
  *        \p values. If empty, no operation is performed.
  *
  * @param issuePMod Controls how physical Modifications are issued.
  *
  * @param issueAMod Controls how abstract Modifications are issued. */

 void set_converter_kappa( MF_dbl_it values ,
                           Range rng = Range( 0 , Inf< Index >() ) ,
                           c_ModParam issuePMod = eNoBlck ,
                           c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// set the capacity multiplier of the converter
 /** This function sets \f$ \kappa^{cv} \f$, see the other versions.
  *
  * @param value The value of the multiplier.
  *
  * @param issuePMod Controls how physical Modifications are issued.
  *
  * @param issueAMod Controls how abstract Modifications are issued. */

 void set_converter_kappa( double value , c_ModParam issuePMod = eNoBlck ,
                           c_ModParam issueAMod = eNoBlck ) {
  std::vector< double > vector = { value };
  set_converter_kappa( vector.cbegin() , Range( 0 , Inf< Index >() ) ,
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

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

/*---------------------------------- data ----------------------------------*/

 /// the vector of minimum storage
 std::vector< double > v_MinStorage;

 /// the vector of maximum storage
 std::vector< double > v_MaxStorage;

 /// the vector of MinPower
 std::vector< double > v_MinPower;

 /// the vector of MaxPower
 std::vector< double > v_MaxPower;

 /// the vector of MinReactivePower
 std::vector< double > v_MinReactivePower;

 /// the vector of MaxReactivePower
 std::vector< double > v_MaxReactivePower;

 /// the vector of ConverterMaxPower
 std::vector< double > v_ConvMaxPower;

 /// the vector of MaxPrimaryPower
 std::vector< double > v_MaxPrimaryPower;

 /// the vector of MaxSecondaryPower
 std::vector< double > v_MaxSecondaryPower;

 /// the vector of RampUp
 std::vector< double > v_DeltaRampUp;

 /// the vector of RampDown
 std::vector< double > v_DeltaRampDown;

 /// the vector of StoringBatteryRho
 std::vector< double > v_StoringBatteryRho;

 /// the vector of ExtractingBatteryRho
 std::vector< double > v_ExtractingBatteryRho;

 /// the vector of StandingBatteryRho
 std::vector< double > v_StandingBatteryRho;

 /// the vector of Cost
 std::vector< double > v_Cost;

 /// the vector of demand
 std::vector< double > v_Demand;

 /// the battery investment cost
 double f_BattInvestmentCost;

 /// the converter investment cost
 double f_ConvInvestmentCost;

 /// the minimum battery capacity design (lower bound of x_b), default 0
 double f_BattMinCapacityDesign;

 /// the maximum battery capacity design, default 1
 /** If >= 0, x_b is continuous in [ BatteryMinCapacityDesign ,
  * BatteryMaxCapacityDesign ]; if < 0, x_b is integer in
  * { BatteryMinCapacityDesign , ... , | BatteryMaxCapacityDesign | }. */
 double f_BattMaxCapacityDesign;

 /// the minimum converter capacity design (lower bound of x_c), default 0
 double f_ConvMinCapacityDesign;

 /// the maximum converter capacity design, default 1
 /** As f_BattMaxCapacityDesign, for x_c. */
 double f_ConvMaxCapacityDesign;

 /// the InitialStorage value
 double f_InitialStorage;

 /// the InitialPower value
 double f_InitialPower;

 /// the MaxCRateCharge value
 double f_MaxCRateCharge;

 /// the MaxCRateDischarge value
 double f_MaxCRateDischarge;

 /// the kappa value
 double f_kappa;

 /// the kappa of the converter, NaN if it follows f_kappa
 double f_conv_kappa;

 /// the scale factor
 double f_scale;

 /// true if "ConverterMaxPower" was in the data
 bool f_conv_max_power_given;

 /// the reference Schedule to deviate minimally from if there
 std::vector< double > v_RefSchedule;

/*-------------------------------- variables -------------------------------*/

 /// the vector of storage level variables
 std::vector< ColVariable > v_storage_level;

 /// the vector of intake level variables
 std::vector< ColVariable > v_intake_level;

 /// the vector of outtake level variables
 std::vector< ColVariable > v_outtake_level;

 /// the vector of binary variables
 std::vector< ColVariable > v_battery_binary;

 /// the active power variables
 std::vector< ColVariable > v_active_power;

 /// the reactive power variables
 std::vector< ColVariable > v_reactive_power;

 /// the primary spinning reserve variables
 std::vector< ColVariable > v_primary_spinning_reserve;

 /// the secondary spinning reserve variables
 std::vector< ColVariable > v_secondary_spinning_reserve;

 /// the battery design variable
 ColVariable batt_design;

 /// the converter design variable
 ColVariable conv_design;

 /// the variables for deviation to reference schedule
 std::vector< ColVariable > v_abs_ref_schedule;

/*------------------------------- constraints ------------------------------*/

/// the reference schedule constraints
 std::vector< FRowConstraint > Reference_Schedule_Const;

 /// the active power bounds constraints
 boost::multi_array< FRowConstraint , 2 > active_power_bounds_Const;

 /// the battery design bound constraint
 BoxConstraint batt_design_bound_Const;

 /// the converter design bound constraint
 BoxConstraint conv_design_bound_Const;

 /// the active power bounds design constraints
 boost::multi_array< FRowConstraint , 2 > active_power_bounds_design_Const;

 /// the intake outtake upper bounds design constraints
 boost::multi_array< FRowConstraint , 2 >
                                   intake_outtake_upper_bounds_design_Const;

 /// the converter bound constraints of a battery with no design Variable
 boost::multi_array< FRowConstraint , 2 > converter_bounds_Const;

 /// the storage level bounds design constraints
 boost::multi_array< FRowConstraint , 2 > storage_level_bounds_design_Const;

 /// the intake and outtake binary variable relation constraints
 boost::multi_array< FRowConstraint , 2 > intake_outtake_binary_Const;

 /// the active power, intake and outtake relation constraints
 std::vector< FRowConstraint > power_intake_outtake_Const;

 /// the ramp up constraints
 std::vector< FRowConstraint > ramp_up_Const;

 /// the ramp down constraints
 std::vector< FRowConstraint > ramp_down_Const;

 /// the demand constraints
 std::vector< FRowConstraint > demand_Const;

 /// the storage level bound constraints
 std::vector< BoxConstraint > storage_level_bounds_Const;

 /// the intake and outtake bounds constraints
 /** Two-sided: folded onto the signed active power, the pair of one-sided
  * fences on the intake and the outtake becomes a single bound whose LHS is
  * the (negative) charging side. */
 boost::multi_array< BoxConstraint , 2 > intake_outtake_bounds_Const;

 /// primary upper bound constraints
 std::vector< LB0Constraint > primary_upper_bound_Const;

 /// secondary upper bound constraints
 std::vector< LB0Constraint > secondary_upper_bound_Const;
 
 /// primary reserve rows of the design, p^pr_t <= kappa P^pr_t x_b
 std::vector< FRowConstraint > primary_upper_bound_design_Const;

 /// secondary reserve rows of the design, p^sc_t <= kappa P^sc_t x_b
 std::vector< FRowConstraint > secondary_upper_bound_design_Const;

 /// the vector of binary bound constraints
 std::vector< ZOConstraint > battery_binary_bound_Const;

 /// the reactive power bound constraints
 std::vector< BoxConstraint > ReactivePower_Bound_Const;

 /*!! Q <= P
 std::vector< FRowConstraint > Reactive_2_Active_Const;
 */

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

 /// the factor of the binary rows (8) and of the reserve bounds (9)
 /** \f$ | \bar{X} | \f$ ("BatteryMaxCapacityDesign" in absolute value) if
  * the battery design variable exists, 1 otherwise. */

 double batt_design_mult( void ) const {
  return( f_BattInvestmentCost != 0 ? std::abs( f_BattMaxCapacityDesign )
	                            : 1.0 );
  }

/*--------------------------------------------------------------------------*/
 /// the design Variable that multiplies the converter power, if any
 /** The converter design \f$ x_c \f$ if it exists; otherwise, for a
  * battery with a design \f$ x_b \f$ and no "ConverterMaxPower", the
  * battery design itself, each module having its converter of the power of
  * the module; nullptr if the converter power is a constant. */

 ColVariable * conv_scale_var( void ) {
  if( f_ConvInvestmentCost != 0 )
   return( & conv_design );
  if( ( f_BattInvestmentCost != 0 ) && ( ! f_conv_max_power_given ) )
   return( & batt_design );
  return( nullptr );
  }

/*--------------------------------------------------------------------------*/
 /// updates the constraints for the current initial storage
 /** This function updates both sides of the demand constraint at time 0
  * (which is the constraint that depends on the initial storage). */

 void update_initial_storage_in_cnstrs( c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// throws if a new InitialStorage would change the form of the first row
 /** The sign of "InitialStorage" selects whether the storage balance of
  * instant 0 is closed on the level of the last instant (cyclic case) or
  * has a constant right-hand side; after the Constraints are generated only
  * the latter can follow a change, hence a change of sign throws. */

 void check_initial_storage_sign( double value ) const;

/*--------------------------------------------------------------------------*/
 /// updates the constraints for the current initial power
 /** This function updates the right-hand side of the ramp-up constraints and
  * the left-hand side of the ramp-down constraints at time 0 (which are the
  * constraints that depend on the initial power). */

 void update_initial_power_in_cnstrs( c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// updates the constraints for the current kappa
 /** This function updates the constraints to take into account the current
  * value of the kappa constant, those of the converter only if it follows
  * kappa [see update_conv_kappa_rows()]. */

 void update_kappa_in_cnstrs( ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// updates the rows of the converter for the current converter kappa
 /** This function writes get_converter_kappa() times \f$ P^{cv}_t \f$ into
  * the rows (10) or (10d), as a right-hand side or as the coefficient of the
  * design variable, with the abstract Modification parameter \p nAM, whose
  * channel (if any) is opened and closed by the caller. */

 void update_conv_kappa_rows( ModParam nAM );

/*--------------------------------------------------------------------------*/
 /// tells whether the rows (10) or (10d) of the converter are generated
 /** The converter rows are generated if "ConverterMaxPower" or
  * "ConverterKappa" is given or a design variable exists [see
  * generate_abstract_constraints()]. */

 bool has_converter_rows( void ) const {
  return( f_conv_max_power_given || ( ! converter_follows_kappa() ) ||
          ( f_BattInvestmentCost != 0 ) || ( f_ConvInvestmentCost != 0 ) );
  }

/*--------------------------------------------------------------------------*/
 /// updates the coefficients of the Objective
 /** This method updates the coefficients of the Objective.
  *
  * @param issueAMod Controls how abstract Modifications are issued. */

 void update_objective( c_ModParam issueAMod ) const;

/*--------------------------------------------------------------------------*/
 /// tells whether the active power has to be split in intake and outtake
 /** The intake and the outtake level are the negative and the positive part
  * of the active power, and only the level balance can tell them apart: it
  * weights them with StoringBatteryRho and ExtractingBatteryRho, and the
  * operating cost already prices them symmetrically. When the two weights
  * agree at all times the balance reads the same on the active power alone,
  * ( il , ol ) = ( max( 0 , -p ) , max( 0 , p ) ) satisfies every row that
  * fences the pair, and the split is unnecessary.
  *
  * @return true if the intake and outtake level Variable are needed. */

 bool needs_intake_outtake( void ) const {
  if( v_StoringBatteryRho.empty() && v_ExtractingBatteryRho.empty() )
   return( false );

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   const double storing = v_StoringBatteryRho.empty()
                          ? 1. : v_StoringBatteryRho[ t ];
   const double extracting = v_ExtractingBatteryRho.empty()
                             ? 1. : v_ExtractingBatteryRho[ t ];
   if( storing != extracting )
    return( true );
   }

  return( false );
  }

/*--------------------------------------------------------------------------*/
 /// verify whether the data in this BatteryUnitBlock is consistent
 /** Checks the data of this BatteryUnitBlock, with the symbols of the class
  * description, and throws std::logic_error (or std::invalid_argument) if
  * one of the following conditions does not hold:
  *
  * - \f$ P^{mn}_t \leq P^{mx}_t \f$ for all \f$ t \f$, and
  *   \f$ P^{mn}_t \leq 0 \f$ if the split is generated
  *   [see needs_intake_outtake()];
  *
  * - \f$ V^{mn}_t \leq V^{mx}_t \f$ for all \f$ t \f$ (equality allowed,
  *   and \f$ V^{mn}_t \f$ of any sign);
  *
  * - \f$ \rho^{in}_t \geq 0 \f$, \f$ \rho^{out}_t \geq 0 \f$,
  *   \f$ \rho^{in}_t \leq \rho^{out}_t \f$ (each may exceed 1) and
  *   \f$ \rho^{st}_t \in [ 0 , 1 ] \f$ for all \f$ t \f$, the absent ones
  *   being 1;
  *
  * - \f$ \Delta^+_t \f$, \f$ \Delta^-_t \f$, \f$ P^{pr}_t \f$,
  *   \f$ P^{sc}_t \f$ and \f$ d^{ba}_t \f$ nonnegative, and
  *   \f$ \kappa \geq 0 \f$;
  *
  * - \f$ V^0 \leq 0 \f$ if an investment cost is nonzero;
  *
  * - for the battery and for the converter, \f$ \underline{X} \geq 0 \f$,
  *   \f$ \underline{X} \leq | \bar{X} | \f$, and
  *   \f$ | \bar{X} | \leq 1 \f$ if \f$ \sigma \neq 1 \f$. */

 void check_data_consistency( void ) const;

/*--------------------------------------------------------------------------*/

 static void static_initialization( void )
 {
  register_method< BatteryUnitBlock , MF_dbl_it , Subset && , bool >(
   "BatteryUnitBlock::set_initial_storage" ,
   & BatteryUnitBlock::set_initial_storage );

  register_method< BatteryUnitBlock , MF_dbl_it , Range >(
   "BatteryUnitBlock::set_initial_storage" ,
   & BatteryUnitBlock::set_initial_storage );

  register_method< BatteryUnitBlock , MF_dbl_it , Subset && , bool >(
   "BatteryUnitBlock::set_initial_power" ,
   & BatteryUnitBlock::set_initial_power );

  register_method< BatteryUnitBlock , MF_dbl_it , Range >(
   "BatteryUnitBlock::set_initial_power" ,
   & BatteryUnitBlock::set_initial_power );

  register_method< BatteryUnitBlock , MF_dbl_it , Subset && , bool >(
   "BatteryUnitBlock::set_cost" ,
   & BatteryUnitBlock::set_cost );

  register_method< BatteryUnitBlock , MF_dbl_it , Range >(
   "BatteryUnitBlock::set_cost" ,
   & BatteryUnitBlock::set_cost );

  register_method< BatteryUnitBlock , MF_dbl_it , Subset && , bool >(
   "BatteryUnitBlock::set_kappa" ,
   & BatteryUnitBlock::set_kappa );

  register_method< BatteryUnitBlock , MF_dbl_it , Range >(
   "BatteryUnitBlock::set_kappa" ,
   & BatteryUnitBlock::set_kappa );

  register_method< BatteryUnitBlock , MF_dbl_it , Subset && , bool >(
   "BatteryUnitBlock::set_converter_kappa" ,
   & BatteryUnitBlock::set_converter_kappa );

  register_method< BatteryUnitBlock , MF_dbl_it , Range >(
   "BatteryUnitBlock::set_converter_kappa" ,
   & BatteryUnitBlock::set_converter_kappa );
  }

};  // end( class( BatteryUnitBlock ) )

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS BatteryUnitBlockMod ------------------------*/
/*--------------------------------------------------------------------------*/

/// derived class from Modification for modifications to a BatteryUnitBlock
class BatteryUnitBlockMod : public UnitBlockMod
{

 public:

 /// public enum for the types of BatteryUnitBlockMod
 enum BUB_mod_type
 {
  eSetInitS = eUBModLastParam , ///< set initial storage values
  eSetInitP ,                   ///< set initial power values
  eSetCost ,                    ///< set cost values
  eSetKappa ,                   ///< set the kappa constant
  eSetConvKappa ,               ///< set the kappa of the converter
  eBUBModLastParam  ///< first allowed parameter value for derived classes
  /**< Convenience value to easily allow derived classes to extend the set of
   * types of BatteryUnitBlockMod. */
  };

 /// constructor, takes the BatteryUnitBlock and the type
 BatteryUnitBlockMod( BatteryUnitBlock * const fblock , const int type )
  : UnitBlockMod( fblock , type ) , f_Block( fblock ) {}

 /// destructor, does nothing
 virtual ~BatteryUnitBlockMod() override = default;

 /// returns the Block to which the Modification refers
 Block * get_Block( void ) const override { return( f_Block ); }

 protected:

 /// prints the BatteryUnitBlockMod
 void print( std::ostream & output ) const override {
  output << "BatteryUnitBlockMod[" << this << "]: ";
  switch( f_type ) {
   case( eSetInitS ):
    output << "set initial storage values ";
    break;
   case( eSetInitP ):
    output << "set initial power values ";
    break;
   case( eSetCost ):
    output << "set cost values ";
    break;
   case( eSetKappa ):
    output << "set kappa ";
    break;
   case( eSetConvKappa ):
    output << "set converter kappa ";
    break;
   default:;
  }
 }

 BatteryUnitBlock * f_Block{};
 ///< pointer to the Block to which the Modification refers

};  // end( class( BatteryUnitBlockMod ) )

/*--------------------------------------------------------------------------*/
/*--------------------- CLASS BatteryUnitBlockRngdMod ----------------------*/
/*--------------------------------------------------------------------------*/

/// derived from BatteryUnitBlockMod for "ranged" modifications
class BatteryUnitBlockRngdMod : public BatteryUnitBlockMod
{

 public:

 /// constructor: takes the BatteryUnitBlock, the type, and the range
 BatteryUnitBlockRngdMod( BatteryUnitBlock * const fblock ,
                          const int type ,
                          const Block::Range & rng )
  : BatteryUnitBlockMod( fblock , type ) , f_rng( rng ) {}

 /// destructor, does nothing
 virtual ~BatteryUnitBlockRngdMod() override = default;

 /// accessor to the range
 Block::c_Range & rng( void ) const { return( f_rng ); }

 protected:

 /// prints the BatteryUnitBlockRngdMod
 void print( std::ostream & output ) const override {
  BatteryUnitBlockMod::print( output );
  output << "[ " << f_rng.first << ", " << f_rng.second << " )" << std::endl;
  }

 Block::Range f_rng;  ///< the range

};  // end( class( BatteryUnitBlockRngdMod ) )

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS BatteryUnitBlockSbstMod ---------------------*/
/*--------------------------------------------------------------------------*/

/// derived from BatteryUnitBlockMod for "subset" modifications
class BatteryUnitBlockSbstMod : public BatteryUnitBlockMod
{

 public:

 /// constructor: takes the BatteryUnitBlock, the type, and the subset
 BatteryUnitBlockSbstMod( BatteryUnitBlock * const fblock ,
                          const int type ,
                          Block::Subset && nms )
  : BatteryUnitBlockMod( fblock , type ) , f_nms( std::move( nms ) ) {}

 /// destructor, does nothing
 virtual ~BatteryUnitBlockSbstMod() override = default;

 /// accessor to the subset
 Block::c_Subset & nms( void ) const { return( f_nms ); }

 protected:

 /// prints the BatteryUnitBlockSbstMod
 void print( std::ostream & output ) const override {
  BatteryUnitBlockMod::print( output );
  output << "(# " << f_nms.size() << ")" << std::endl;
  }

 Block::Subset f_nms;  ///< the subset

};  // end( class( BatteryUnitBlockSbstMod ) )

/*--------------------------------------------------------------------------*/
/*-------------------- CLASS BatteryUnitBlockSolution ----------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// a [UnitBlock]Solution of a BatteryUnitBlock
/** The BatteryUnitBlockSolution class derives from UnitBlockSolution and
 * adds to the "standard" information stored in there (active power, possibly
 * commitment and primary/secondary reserve) the other information that is
 * typical of the BatteryUnitBlock, i.e.,
 *
 * - [possibly] the storage level of the battery at each time instant
 *
 * - [possibly] the intake/outtake in the battery at each time instant;
 *   since the battery is supposed to never be charged and discharged at
 *   the same time instant, the value is positive if the battery is being
 *   charged (intake) and negative if it is being discharged (outtake)
 *
 * - if defined, the value of the Battery Design Variable
 *
 * - if defined, the value of the Converter Design Variable */

class BatteryUnitBlockSolution : public UnitBlockSolution
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*----------------------------- CONSTANTS ----------------------------------*/

 static constexpr double dNaN = std::numeric_limits< double >::quiet_NaN();
 ///< convenience constexpr for "NaN", *not* to be used with ==

/*------------------------------- FRIENDS ----------------------------------*/

 friend BatteryUnitBlock;  ///< make BatteryUnitBlock friend

/*--------- CONSTRUCTING AND DESTRUCTING BatteryUnitBlockSolution ----------*/

 /// constructor, it has nothing to do
 explicit BatteryUnitBlockSolution( void ) :
  f_b_design( dNaN ) , f_c_design( dNaN ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void deserialize( const netCDF::NcGroup & group ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 ~BatteryUnitBlockSolution() override = default;
 ///< destructor: it is virtual, and empty

/*----- METHODS DESCRIBING THE BEHAVIOR OF A BatteryUnitBlockSolution -----*/

 void read( const Block * block ) override final;

 void write( Block * block ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// serialize a BatteryUnitBlockSolution into a netCDF::NcGroup
 /** Serialize a BatteryUnitBlockSolution into a netCDF::NcGroup. The format
  * is the one of UnitBlockSolution [cf. UnitBlockSolution::serialize()],
  * plus:
  *
  * - The variable "StorageLevel", of type netCDF::NcDouble and indexed over
  *   the dimension "TimeHorizon"; StorageLevel[ t ] is the optimal value of
  *   the storage level of the battery at time \f$ t \f$. The variable is
  *   optional.
  *
  * - The variable "InOutTake", of type netCDF::NcDouble and indexed over
  *   the dimension "TimeHorizon"; InOutTake[ t ] is the amount of energy
  *   being charged in the battery at time \f$ t \f$ (negative if it is
  *   discharged). The variable is optional.
  *
  * Note that, unlike those of the base class, these variables do not need
  * to be indexed over the dimension "NumberGenerators" since
  * BatteryUnitBlock always has exactly one generator.
  *
  * - The scalar variable "BatteryDesign", of type netCDF::NcDouble, that
  *   represent the value of the dimensioning variable of the battery;
  *   the variable is optional in that the battery may not have any
  *   dimensioning variable.
  *
  * - The scalar variable "ConverterDesign", of type netCDF::NcDouble, that
  *   represent the value of the dimensioning variable of the converter;
  *   the variable is optional in that the converter may not have any
  *   dimensioning variable. */

 void serialize( netCDF::NcGroup & group ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 BatteryUnitBlockSolution * scale( double factor ) const override;

 void sum( const Solution * solution , double multiplier ) override;

 BatteryUnitBlockSolution * clone( bool empty = false ) const override;

/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/

 void print( std::ostream & output ) const override {
  output << "BatteryUnitBlockSolution [" << this << "]: " << std::endl;
  }

/*---------------------- PRIVATE PART OF THE CLASS -------------------------*/

 private:

/*---------------------------- PRIVATE FIELDS ------------------------------*/

 std::vector< double > v_storage;
 ///< v_storage[ t ] = value of stored energy at time t

 std::vector< double > v_intake;  ///< v_intake[ t ] = intake at time t

 double f_b_design;    ///< the value of the battery dimensioning variable

 double f_c_design;    ///< the value of the converter dimensioning variable

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

};  // end( class( BatteryUnitBlockSolution ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* __BatteryUnitBlock */

/*--------------------------------------------------------------------------*/
/*---------------------- End File BatteryUnitBlock.h -----------------------*/
/*--------------------------------------------------------------------------*/
