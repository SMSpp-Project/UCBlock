/*--------------------------------------------------------------------------*/
/*------------------------ File NuclearUnitExtDPSolver.h ----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the NuclearUnitExtDPSolver class (NUEDPS), a Solver for the
 * NuclearUnitBlock that solves the single-Unit Commitment (1UC) problem of a
 * nuclear (thermal + modulation) unit by Dynamic Programming.
 *
 * It is built directly on top of ThermalUnitExtDPSolver (the "Ext" hybrid DP,
 * multi-layer ON / single-layer OFF), reusing its entire piecewise-quadratic
 * machinery (PQFun, sliding_min, add_quadratic, add_pwq, clamp_domain, the
 * domination check) and the spinning-reserve handling (reserve_alloc /
 * build_reserve_discount). The only structural addition is the *modulation
 * lockout counter* needed to handle the extra constraints of a nuclear unit.
 *
 * A nuclear unit (see NuclearUnitBlock) is a thermal unit on which, during a
 * continuous on-run, the power may move between two consecutive on instants by
 * at most the (small) *modulation* ramp Delta^M, unless a *modulation* is
 * performed at that instant (the binary m_t), in which case the full thermal
 * ramp Delta is allowed; at most one modulation may occur in any window of
 * ModulationTime (tau^M) consecutive instants. Start-ups and shut-downs are
 * never modulations and keep their usual bound_on / bound_down caps.
 *
 * The window constraint couples the m_t across time, so (unlike the reserves)
 * it cannot be folded into a per-period effective cost. It is instead enforced
 * by a single bounded counter
 *
 *      ell in { 0, 1, ..., tau^M - 1 }
 *
 * recording how many more instants the unit is locked out from modulating; ell
 * is added to the DP state, giving on-side value functions F^{tau,ell}_t and
 * an ell-indexed "ready-to-restart" off-side scalar. See the design document
 * (reserves-DP, "Extending the DP to nuclear units") for the full derivation.
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

#ifndef __NuclearUnitExtDPSolver
 #define __NuclearUnitExtDPSolver
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "ThermalUnitExtDPSolver.h"

#include "NuclearUnitBlock.h"

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS NuclearUnitExtDPSolver ------------------------*/
/*--------------------------------------------------------------------------*/
/// DP solver for the 1UC of a nuclear unit (thermal + modulation)

class NuclearUnitExtDPSolver : public ThermalUnitExtDPSolver
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

 NuclearUnitExtDPSolver( void ) : ThermalUnitExtDPSolver() {}

 ~NuclearUnitExtDPSolver() override = default;

/*--------------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

 /// sets the (NuclearUnitBlock) that the Solver has to solve
 void set_Block( Block * block ) override;

 /// writes the current solution (incl. the modulation m_t) into the Block
 void get_var_solution( Configuration * solc ) override;

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*----------------------- PROTECTED METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

 /// read the base parameters and, on top, the modulation data
 void load_parameters( void ) override;

 /// intercept the nuclear-specific Modifications (else delegate to base)
 bool guts_of_process_modifications( const p_Mod mod ) override;

 /// run the forward DP with the extra modulation-lockout dimension
 void run_DP( void ) override;

 /// backtrack the DP, recovering (p, u) and the modulation m_t
 void build_solution( void ) override;

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 // ---- modulation data loaded from the NuclearUnitBlock ---------------- //

 Index f_mod_interval{ 2 };   ///< modulation interval tau^M ( >= 2 )
 int   f_init_modulation{ 2 };///< InitModulation kappa^0 ( >= 1 )
 std::vector< double > mod_ramp_up;   ///< modulation Delta^M+ per instant
 std::vector< double > mod_ramp_down; ///< modulation Delta^M- per instant

 // ---- ON-side extra state, parallel to f_F[t] / f_tau[t] / f_on[t] ---- //

 /// f_lock[t][i] : lockout ell of slot i at time t, i.e. the number of
 /// instants (counting from t+1) for which a modulation is still forbidden
 /// *after* the decision taken at t. ell == tau^M - 1 means a modulation was
 /// performed at t; ell == 0 means modulation is free again from t+1.
 std::vector< std::vector< Index > > f_lock;

 /// back-pointer for backtracking: for a tau > 1 slot it is the lockout ell'
 /// of the predecessor on-slot (t-1, tau-1, ell'); for a tau == 1 (restart)
 /// slot it is the off-side "ready" lockout ell' at t-1 from which the
 /// restart was taken. Lets build_solution find the unique predecessor.
 std::vector< std::vector< Index > > f_back_lock;

 /// back-pointer index into the (final, post-prune) slot list at t-1 of the
 /// predecessor on-slot, for tau > 1 slots; unused (BAD) for tau == 1 / seed.
 std::vector< std::vector< std::size_t > > f_back_idx;

 static constexpr std::size_t BAD = std::size_t( -1 );

 // ---- OFF-side state with the lockout dimension ----------------------- //
 // Only the "ready-to-restart" scalar needs the ell index, because it is the
 // only off-side quantity that feeds a restart F^1. c_off_any (base member,
 // 1D) never feeds a restart and so stays ell-agnostic (min over ell).

 /// co_ready[t][ell] : min cost of a schedule off at t, off for at least
 /// min_down_time consecutive instants (legal to restart at t+1), with
 /// lockout ell into t+1.
 std::vector< std::vector< double > > co_ready;

 /// ready_pred[t][ell] : origin time h of the shutdown that produced
 /// co_ready[t][ell] through the long arc, or -1 (initial off trail / +INF).
 std::vector< std::vector< int > > ready_pred;

 /// ready_pred_lock[t][ell] : lockout ell_h *at the shutdown* (into h+1) that
 /// produced co_ready[t][ell]; needed to locate the originating on-slot.
 std::vector< std::vector< Index > > ready_pred_lock;

 /// fresh-shutdown lockout backing c_off_any[t] (base 1D member): the ell_h
 /// of the v_shutdown that contributed; needed for an off-termination walk.
 std::vector< Index > f_any_lock;

 /// vs[h][ell] : cost of reaching the long-shutdown arc at the end of time h
 /// with lockout ell into h+1 (i.e. min over surviving (tau >= mut) on-slots
 /// at h with lockout ell of F^{tau,ell}_h on [P, SD_{h+1}]). vs_tau / vs_p /
 /// vs_idx record the achieving (tau, p, slot index) for backtracking.
 std::vector< std::vector< double > > vs;
 std::vector< std::vector< Index  > > vs_tau;
 std::vector< std::vector< double > > vs_p;
 std::vector< std::vector< std::size_t > > vs_idx;

 // ---- output: recovered modulation profile --------------------------- //

 /// M[t] : the optimal modulation indicator of the recovered schedule
 std::vector< bool > M;

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 private:

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class( NuclearUnitExtDPSolver ) )

};  // end( namespace SMSpp_di_unipi_it )

#endif  /* NuclearUnitExtDPSolver.h included */

/*--------------------------------------------------------------------------*/
/*------------------- End File NuclearUnitExtDPSolver.h -----------------------*/
/*--------------------------------------------------------------------------*/
