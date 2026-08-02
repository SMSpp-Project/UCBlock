/*--------------------------------------------------------------------------*/
/*--------------------- File ThermalUnitExtDPSolver.cpp --------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of ThermalUnitExtDPSolver: a hybrid DP solver for the
 * single Unit Commitment problem in which the on-side of the graph keeps
 * the multi-layer \f$ \tau \f$-counter of Wuijts et al. (2021) while the
 * off-side is collapsed to a single layer and the min-down-time
 * constraint is enforced by a "long" arc that jumps \f$ M_{down} \f$
 * instants ahead.
 *
 * The ON layer stores convex piecewise quadratic functions
 * \f$ F^\tau_t \f$: the sliding minimum and domain clamping follow the
 * standard Frangioni-Gentile three-case analysis. The OFF layer stores
 * two scalars per time step (c_off_ready[t], c_off_any[t]) plus the
 * auxiliary sequence v_shutdown[t] (cost of a schedule on through t plus
 * the shutdown cost, ready to feed the long shutdown arc).
 *
 * Data loading, Modification handling and Solution writing follow the
 * same pattern as ThermalUnitDPSolver, but are implemented independently
 * (no code shared at build time).
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
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

// If TUEDPS_PROFILE > 0, compute() accumulates per-phase wall time
// (run_DP / build_solution) into static counters and prints a cumulative
// summary to cerr at each call. Temporary instrumentation: keep at 0.
// Defined BEFORE the header include so the profiling counters at the end
// of ThermalUnitDPSolverBase.h (under #if TUEDPS_PROFILE) get a single
// shared definition across the two translation units.
#define TUEDPS_PROFILE 0

#include "ThermalUnitExtDPSolver.h"

#include "ThermalUnitBlock.h"

#include <algorithm>
#include <cmath>
#include <functional>

#if TUEDPS_PROFILE
 #include <chrono>
 #include <iostream>
#endif

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_0( ThermalUnitExtDPSolver );

/*--------------------------------------------------------------------------*/
/*--------------------------- Solver INTERFACE -----------------------------*/
/*--------------------------------------------------------------------------*/

// Attach this Solver to a ThermalUnitBlock. The concrete type of the Block
// is checked with typeid() (as opposed to dynamic_cast<>) to intentionally
// *exclude* classes derived from ThermalUnitBlock: such derivatives might
// add constraints or variables (e.g., primary/secondary reserves) whose
// semantics are not understood by this DP. When a matching Block is
// attached, all input parameters are loaded eagerly so that run_DP() can
// use them directly at the next compute().

void ThermalUnitExtDPSolver::set_Block( Block * block )
{
 if( block == f_Block )
  return;

 Solver::set_Block( block );

 if( block ) {
  if( typeid( ThermalUnitBlock ) != typeid( *f_Block ) )
   throw( std::runtime_error(
    "ThermalUnitExtDPSolver::set_Block: ThermalUnitBlock required." ) );
  load_parameters();
  }
 }

/*--------------------------------------------------------------------------*/

// Main entry point of the Solver interface. The computation is organised
// in two stages (forward DP, backtracking) so that repeated calls with no
// intervening Modifications short-circuit on the cached state. Before
// touching the state we drain the Modification queue: any Modification
// that invalidates the cached DP resets stage back to start, forcing the
// relevant stages to re-run. Returns kOK when a feasible solution was
// found, kInfeasible when even the "unit off for the whole horizon"
// schedule is infeasible (should not happen under standard inputs).

int ThermalUnitExtDPSolver::compute( bool changedvars )
{
 lock();

 process_modifications();

#if TUEDPS_PROFILE
 static double t_dp = 0 , t_bs = 0;
 static unsigned long n_dp = 0 , n_bs = 0 , n_call = 0;
 using clk = std::chrono::steady_clock;
 auto tic = clk::now();
 if( stage < dp_OK ) {
  static const long REP = std::getenv( "TUEDPS_REPEAT" )
                          ? std::atol( std::getenv( "TUEDPS_REPEAT" ) ) : 1;
  for( long r = 1 ; r < REP ; ++r ) {  // extra reps for profiling
   run_DP();
   stage = start;
   }
  run_DP();
  t_dp += std::chrono::duration< double >( clk::now() - tic ).count() / REP;
  ++n_dp;
  tic = clk::now();
  }
 if( stage < sol_OK ) {
  build_solution();
  t_bs += std::chrono::duration< double >( clk::now() - tic ).count(); ++n_bs;
  }
 ++n_call;
 std::cerr << "TUEDPS_PROF n=" << time_horizon << " calls=" << n_call
           << " run_DP=" << t_dp << "s/" << n_dp
           << " build_solution=" << t_bs << "s/" << n_bs
           << " | smc=" << g_smc << " gmin=" << g_gmin << " node=" << g_node
           << " Favg=" << ( g_smc ? double( g_Fsum ) / g_smc : 0 )
           << " Fmax=" << g_Fmax
           << " gmin/smc=" << ( g_smc ? double( g_gmin ) / g_smc : 0 )
           << " node/smc=" << ( g_smc ? double( g_node ) / g_smc : 0 )
           << " pcs/smc=" << ( g_smc ? double( g_pieces ) / g_smc : 0 )
           << " | BITE smc_bite%="
           << ( g_smc ? 100.0 * g_smc_bite / g_smc : 0 )
           << " gmin_bite%=" << ( g_gmin ? 100.0 * g_gmin_bite / g_gmin : 0 )
           << " fast/node%=" << ( g_node ? 100.0 * g_fast / g_node : 0 )
           << " | PARAM calls=" << g_param_calls << " bad=" << g_param_bad
           << " maxdev(vsTrue)=" << g_param_maxdev
           << " oracle_maxdev(vsTrue)=" << g_param_ormaxdev
           << " gmin_maxdev(vsTrue)=" << g_gmin_maxdev
           << " pcs/call="
           << ( g_param_calls ? double( g_param_pcs ) / g_param_calls : 0 )
           << std::endl;
#else
 if( stage < dp_OK )
  run_DP();

 if( stage < sol_OK )
  build_solution();
#endif

 unlock();

 assert( stage == sol_OK );
 return( f_best_cost >= TUEDPINF ? kInfeasible : kOK );
 }

/*--------------------------------------------------------------------------*/

// Write the optimal schedule stored in P[] and U[] into the variables of
// the ThermalUnitBlock. The DP only knows the canonical (active power,
// commitment) representation; the Block is responsible for filling in the
// formulation-specific auxiliaries (start_up / shut_down indicators,
// perspective-cut variables when PCuts is on, ...) consistently from
// (P, U) -- see ThermalUnitBlock::set_solution() for the details.

void ThermalUnitExtDPSolver::get_var_solution( Configuration * solc )
{
 bool owned = f_Block->is_owned_by( f_id );
 if( ( ! owned ) && ( ! f_Block->lock( f_id ) ) )
  throw( std::runtime_error(
   "ThermalUnitExtDPSolver::get_var_solution: unable to lock the Block." ) );

 auto b = static_cast< ThermalUnitBlock * >( f_Block );

 // design (investment): write the design variable and, when the unit is not
 // built, force the whole schedule to zero overriding the DP's (P, U) which
 // were computed assuming the unit exists. "built" is true when there is no
 // design at all.
 const bool built = ( ! has_design ) || design_on;
 if( has_design )
  b->get_design().set_value( design_on ? 1 : 0 );

 // canonical part: active power and commitment from the DP's (P, U)
 if( auto pow_it = b->get_active_power( 0 ) )
  for( Index i = 0 ; i < time_horizon ; ++i )
   ( pow_it++ )->set_value( built ? P[ i ] : 0 );

 if( auto com_it = b->get_commitment( 0 ) )
  for( Index i = 0 ; i < time_horizon ; ++i )
   ( com_it++ )->set_value( ( built && U[ i ] ) ? 1 : 0 );

 // spinning reserve variables (if present): the optimal pr/sr provision
 // given the recovered power profile P[]. The reserve band must be the
 // SAME the DP priced, else the recovered solution would disagree with
 // the reported value. It is the residual-ramp band: the capacity band
 // around P[i], with the cap bound_on at a start-up, bound_down[i+1] at
 // a shut-down (off at i+1), and max_power at an interior period /
 // on-to-end, intersected, at an interior (both i-1 and i on)
 // transition, with the ramp room left over after the scheduled move,
 // min( DRU_{i-1} - (P[i]-P[i-1]) , DRD_{i-1} + (P[i]-P[i-1]) ).
 // There is no ramp term at a start-up (no in-interval predecessor) nor
 // at i == 0 (matching the DP, which prices the reserve of the first
 // on-instant by capacity alone). A start-up is the first on-instant of
 // an on-interval.
 auto res_band = [ & ]( Index i ) -> double {
  if( ! ( built && U[ i ] ) )
   return( 0 );                                      // off: no reserve
  const bool is_su = ( i == 0 ) ? ( init_up_down_time <= 0 )
                                : ( ! U[ i - 1 ] );
  double cap = max_power[ i ];                        // interior / on-to-end
  if( is_su )
   cap = bound_on[ i ];                               // start-up
  else if( ( i + 1 < time_horizon ) && ( ! U[ i + 1 ] ) )
   cap = bound_down[ i + 1 ];                        // shut-down (off at i+1)
  double H = std::min( P[ i ] - min_power[ i ] , cap - P[ i ] );
  if( ( i >= 1 ) && U[ i - 1 ] ) {                    // interior transition
   const double d = P[ i ] - P[ i - 1 ];
   H = std::min( H , std::min( delta_ramp_up[ i - 1 ] - d ,
                               delta_ramp_down[ i - 1 ] + d ) );
   }
  return( H );
  };

 if( auto pr_it = b->get_primary_spinning_reserve( 0 ) )
  for( Index i = 0 ; i < time_horizon ; ++i ) {
   double pr , sr;
   reserve_alloc_band( i , built ? P[ i ] : 0 , res_band( i ) , pr , sr );
   ( pr_it++ )->set_value( pr );
   }

 if( auto sr_it = b->get_secondary_spinning_reserve( 0 ) )
  for( Index i = 0 ; i < time_horizon ; ++i ) {
   double pr , sr;
   reserve_alloc_band( i , built ? P[ i ] : 0 , res_band( i ) , pr , sr );
   ( sr_it++ )->set_value( sr );
   }

 // reactive power variables (AC instances): q[t] is a box-bounded
 // variable in [ Qmin(t) , Qmax(t) ] that is separable from the DP.
 // Under a dualizing Solver it carries the linear cost
 // reactive_linear_term[t]: the optimal q*[t] is the box endpoint
 // minimising c*q (lower bound if c>0, upper if c<0; any feasible value,
 // here 0 clamped to the box, if c==0). This matches the contribution
 // run_DP() adds to the value. See ReactivePower_Bound_Const. The box is
 // state-dependent: [Qmin_off,Qmax_off] when off, widened by the
 // commitment coefficients [Qmin_on,Qmax_on] when on (U[i]), the same
 // box the per-period contribution in run_DP() prices. Off/not built
 // collapses to the off box (q forced to 0 when not built).
 if( auto q_it = b->get_reactive_power( 0 ) )
  for( Index i = 0 ; i < time_horizon ; ++i ) {
   double q = 0;  // not built: q forced to 0 with all operational variables
   if( built ) {
    const bool on = U[ i ];
    const double qlo = b->get_min_reactive_power( i ) +
                       ( on ? b->get_min_reactive_power_on( i ) : 0.0 );
    const double qhi = b->get_max_reactive_power( i ) +
                       ( on ? b->get_max_reactive_power_on( i ) : 0.0 );
    const double c = reactive_linear_term.empty() ? 0.0
                                                  : reactive_linear_term[ i ];
    if( c > 0 )
     q = qlo;
    else if( c < 0 )
     q = qhi;
    else
     q = std::min( std::max( 0.0 , qlo ) , qhi );
    }
   ( q_it++ )->set_value( q );
   }

 // formulation-specific bookkeeping is delegated to the Block
 b->set_solution();

 if( ! owned )
  f_Block->unlock( f_id );

 }  // end( ThermalUnitExtDPSolver::get_var_solution )

// NOTE: the reserve model (reserve_alloc / build_reserve_discount /
// add_pwq / ...) and the piecewise-quadratic machinery are implemented
// in the base class ThermalUnitDPSolverBase
// (see ThermalUnitDPSolverBase.cpp).

/*--------------------------------------------------------------------------*/
/*---------------------- PARAMETER LOADING ---------------------------------*/
/*--------------------------------------------------------------------------*/

// Copy every parameter of the ThermalUnitBlock that the DP will need into
// local fields of this Solver. Everything is fetched in one batch under a
// read-lock of the Block so as to get a consistent snapshot. Vectors that
// the Block exposes as "possibly empty" (meaning the default value should
// be used) or "possibly of size 1" (meaning the scalar should be
// broadcast across the whole horizon) go through retrieve_term() for the
// broadcast. Fields that ThermalUnitDPSolver does not understand
// (primary/secondary reserves) abort loading rather than silently drop
// constraints.

void ThermalUnitExtDPSolver::load_parameters( void )
{
 load_common_parameters();

 // reset output buffers and pipeline state
 P.assign( time_horizon , 0 );
 U.assign( time_horizon , false );
 stage = start;
 f_solved = false;
 f_best_cost = TUEDPINF;

 }  // end( ThermalUnitExtDPSolver::load_parameters )

/*--------------------------------------------------------------------------*/

void ThermalUnitExtDPSolver::load_fixings( void )
{
 f_has_fixings = f_must_build = f_no_build = false;

 // locking the Block
 bool owned = f_Block->is_owned_by( f_id );
 if( ( ! owned ) && ( ! f_Block->read_lock() ) )
  throw( std::runtime_error(
   "ThermalUnitExtDPSolver::load_fixings: unable to lock the Block." ) );

 auto b = static_cast< ThermalUnitBlock * >( f_Block );

 // commitment fixings, translated into the "first instant >= t fixed
 // OFF / ON" tables that run_DP() uses to kill the incompatible states
 if( auto u = b->get_commitment( 0 ) ) {
  nxt_off.assign( time_horizon + 1 , time_horizon );
  nxt_on.assign( time_horizon + 1 , time_horizon );
  for( Index t = time_horizon ; t-- > 0 ; ) {
   nxt_off[ t ] = nxt_off[ t + 1 ];
   nxt_on[ t ] = nxt_on[ t + 1 ];
   if( u[ t ].is_fixed() ) {
    f_has_fixings = true;
    if( u[ t ].get_value() >= 0.5 ) {
     nxt_on[ t ] = t;
     f_must_build = true;  // an ON instant requires the unit to exist
     }
    else
     nxt_off[ t ] = t;
    }
   }
  }

 // the design variable can be fixed, too
 if( has_design ) {
  const auto & d = b->get_const_design();
  if( d.is_fixed() ) {
   if( d.get_value() >= 0.5 )
    f_must_build = true;
   else
    f_no_build = true;
   }
  }

 // fixings of any other Variable cannot be honored by the DP, save the
 // structural ones that ThermalUnitBlock itself makes when generating the
 // variables to encode the initial conditions (the pre-t_init instants and
 // the start-up / shut-down windows right after them), which the DP
 // enforces anyway
 auto refuse = [ & ]( const ColVariable * v , Index n , const char * name ,
		      auto && structural ) {
  for( Index t = 0 ; v && ( t < n ) ; ++t )
   if( v[ t ].is_fixed() && ( ! structural( t , v[ t ].get_value() ) ) ) {
    if( ! owned )
     f_Block->read_unlock();
    throw( std::logic_error(
     std::string( "ThermalUnitExtDPSolver: fixed " ) + name +
     " Variable not supported (yet)" ) );
    }
  };

 const bool init_on = init_up_down_time > 0;
 auto never = []( Index , double ) { return( false ); };
 auto pre_horizon_zero = [ & ]( Index t , double val ) {
  return( ( ! init_on ) && ( t < t_init ) && ( val == 0 ) );
  };

 const Index nsd = time_horizon > t_init ? time_horizon - t_init : 0;
 refuse( b->get_active_power( 0 ) , time_horizon , "active power" ,
	 pre_horizon_zero );
 refuse( b->get_reactive_power( 0 ) , time_horizon , "reactive power" ,
	 never );
 refuse( b->get_primary_spinning_reserve( 0 ) , time_horizon ,
	 "primary reserve" , pre_horizon_zero );
 refuse( b->get_secondary_spinning_reserve( 0 ) , time_horizon ,
	 "secondary reserve" , pre_horizon_zero );
 refuse( b->get_start_up() , nsd , "start-up" ,
	 [ & ]( Index k , double val ) {
	  return( init_on && ( k < min_down_time ) && ( val == 0 ) ); } );
 refuse( b->get_shut_down() , nsd , "shut-down" ,
	 [ & ]( Index k , double val ) {
	  return( ( ! init_on ) && ( k < min_up_time ) && ( val == 0 ) ); } );

 // unlock the Block
 if( ! owned )
  f_Block->read_unlock();

 }  // end( ThermalUnitExtDPSolver::load_fixings )

/*--------------------------------------------------------------------------*/

// Drain the Modification queue produced by the Block since the last call
// to compute(). Each Modification is dispatched through
// guts_of_process_modifications(), which either patches the cached state
// locally (returning false) or signals that the cached state is beyond
// repair and must be rebuilt from scratch (returning true, triggering a
// full load_parameters()). The queue is processed under f_mod_lock so
// that the list stays consistent if the Block concurrently posts more
// Modifications.

void ThermalUnitExtDPSolver::process_modifications( void )
{
 bool reload = false;

 while( f_mod_lock.test_and_set( std::memory_order_acquire ) )
  ;

 for( auto mod : v_mod )
  if( guts_of_process_modifications( mod.get() ) ) {
   reload = true;
   break;        // a full reload moots any subsequent Modification
   }

 v_mod.clear();

 f_mod_lock.clear( std::memory_order_release );

 if( reload )
  load_parameters();

 }  // end( ThermalUnitExtDPSolver::process_modifications )

/*--------------------------------------------------------------------------*/

// Dispatch a single Modification. Returns true iff the caller must
// discard its cached state and re-run load_parameters(); returns false
// when the Modification either has been absorbed in-place (input vector
// refreshed, stage reset) or is unrelated to our DP.
//
// The recognised Modification flavours are:
//
//  - NBModification      : the Block has been fully rebuilt; everything
//                          we had is stale;
//  - GroupModification   : bundle of sub-Modifications; recurse, collect
//                          the worst outcome (reload wins);
//  - ThermalUnitBlockMod : fine-grained, per-field notifications. Most
//                          of them simply re-read the affected field
//                          and reset the stage. eSetAv (availability
//                          profile) is currently not handled locally,
//                          so we fall back to a full reload.
//
// Any other Modification is ignored (returns false): ThermalUnitBlock
// installs handle_objective_change(), which translates abstract-level
// changes to the Objective function (C05FunctionModLin*,
// DQuadFunctionMod*) into physical ThermalUnitBlockMod events that will
// reach us through the queue on the next call.

bool ThermalUnitExtDPSolver::guts_of_process_modifications( const p_Mod mod )
{
 if( dynamic_cast< NBModification * >( mod ) )
  return( true );

 if( const auto gm = dynamic_cast< GroupModification * >( mod ) ) {
  bool reload = false;
  for( const auto & submod : gm->sub_Modifications() )
   if( guts_of_process_modifications( submod.get() ) )
    reload = true;
  return( reload );
  }

 if( const auto tubm = dynamic_cast< ThermalUnitBlockMod * >( mod ) ) {
  auto b = static_cast< ThermalUnitBlock * >( f_Block );
  switch( tubm->type() ) {
   case( ThermalUnitBlockMod::eSetMaxP ):
    max_power = b->get_max_power();
    stage = start;
    return( false );

   case( ThermalUnitBlockMod::eSetInitP ):
    initial_power = b->get_initial_power();
    stage = start;
    return( false );

   case( ThermalUnitBlockMod::eSetInitUD ):
    init_up_down_time = b->get_init_up_down_time();
    min_up_time = b->get_min_up_time();
    min_down_time = b->get_min_down_time();
    if( init_up_down_time > 0 )
     if( min_up_time > Index( init_up_down_time ) )
      t_init = std::min( time_horizon ,
                         min_up_time - Index( init_up_down_time ) );
     else
      t_init = 0;
    else
     if( Index( - init_up_down_time ) < min_down_time )
      t_init = std::min( time_horizon ,
                         min_down_time - Index( - init_up_down_time ) );
     else
      t_init = 0;
    stage = start;
    return( false );

   case( ThermalUnitBlockMod::eSetAv ):
    return( true );

   case( ThermalUnitBlockMod::eSetSUC ):
    startup_costs = b->get_start_up_cost();
    stage = start;
    return( false );

   case( ThermalUnitBlockMod::eSetLinT ):
    retrieve_term( linear_term , b->get_linear_term() );
    stage = start;
    return( false );

   case( ThermalUnitBlockMod::eSetQuadT ):
    retrieve_term( quad_term , b->get_quad_term() );
    stage = start;
    return( false );

   case( ThermalUnitBlockMod::eSetConstT ):
    retrieve_term( const_term , b->get_const_term() );
    stage = start;
    return( false );

   case( ThermalUnitBlockMod::eSetPrSpResCost ):
    primary_reserve_cost = b->get_primary_spinning_reserve_cost();
    if( primary_reserve_cost.empty() )
     primary_reserve_cost = primary_rho;
    stage = start;
    return( false );

   case( ThermalUnitBlockMod::eSetSecSpResCost ):
    secondary_reserve_cost = b->get_secondary_spinning_reserve_cost();
    if( secondary_reserve_cost.empty() )
     secondary_reserve_cost = secondary_rho;
    stage = start;
    return( false );

   case( ThermalUnitBlockMod::eSetReactiveLinT ):
    retrieve_term( reactive_linear_term , b->get_reactive_linear_term() );
    stage = start;
    return( false );

   case( ThermalUnitBlockMod::eSetInvCost ): {
    // the cost of the design (investment) variable has (possibly) changed,
    // e.g. because a dualizing Solver pushed a multiplier into its Objective
    // coefficient. Refresh our copy from get_design_cost(); only
    // invalidate the cached state when the value actually changed, since
    // this is issued whenever the Objective is touched (the design index
    // is always in range).
    auto nc = b->get_design_cost();
    if( nc != design_cost ) {
     design_cost = nc;
     stage = start;
     }
    return( false );
    }

   case( ThermalUnitBlockMod::eFixVars ):
    // the fixed status of some Variable changed: the DP has to be re-run
    // from scratch, as run_DP() re-reads the fixings (via load_fixings())
    // and kills the incompatible states
    stage = start;
    return( false );
   }
  return( true );
  }

 return( false );

 }  // end( ThermalUnitExtDPSolver::guts_of_process_modifications )

/*--------------------------------------------------------------------------*/

// NOTE: retrieve_term(), the piecewise-quadratic helpers (eval /
// argmin_piece / min_over / shift_by / add_quadratic / clamp_domain /
// add_pwq / is_dominated_by / sliding_min) and the reserve model
// (reserve_alloc / reserve_alloc_band / reserve_reward /
// build_reserve_discount / sliding_min_corr / reserve_corr_argmin) are
// implemented in the base class ThermalUnitDPSolverBase
// (ThermalUnitDPSolverBase.cpp).

/*--------------------------------------------------------------------------*/

double ThermalUnitExtDPSolver::eval_allon_cost(
                                    const std::vector< double > & P ) const
{
 double c = 0;
 const Index N = P.size();
 for( Index t = 0 ; t < N ; ++t ) {
  c += ( quad_term[ t ] * P[ t ] + linear_term[ t ] ) * P[ t ] +
       const_term[ t ];
  const double A = std::min( P[ t ] - min_power[ t ] ,
                             max_power[ t ] - P[ t ] );
  double H = ( A > 0 ) ? A : 0;   // t=0: capacity band (no in-horizon ramp)
  if( t > 0 ) {
   const double d = P[ t ] - P[ t - 1 ];
   double B = std::min( delta_ramp_up[ t - 1 ] - d ,
                        delta_ramp_down[ t - 1 ] + d );
   if( B < 0 )
    B = 0;
   H = ( A > 0 ) ? std::min( A , B ) : 0;
   }
  c += reserve_reward( t , P[ t ] , H );
  }
 return( c );
 }

/*--------------------------------------------------------------------------*/

std::vector< double > ThermalUnitExtDPSolver::ff_at_traj(
                                    const std::vector< double > & P ) const
{
 std::vector< double > tr( P.size() , TUEDPINF );
 for( Index t = 0 ; ( t < P.size() ) && ( t < f_F.size() ) ; ++t ) {
  double best = TUEDPINF;
  for( const auto & F : f_F[ t ] )
   for( const auto & pc : F )
    if( ( P[ t ] >= pc.left - 1e-9 ) && ( P[ t ] <= pc.right + 1e-9 ) ) {
     const double v = eval_piece( pc , P[ t ] );
     if( v < best )
      best = v;
     break;
     }
  tr[ t ] = best;
  }
 return( tr );
 }

/*--------------------------------------------------------------------------*/

void ThermalUnitExtDPSolver::dump_states_at( Index t , double p ) const
{
 std::cerr.precision( 10 );
 std::cerr << "STATES t=" << t << " p=" << p << " (#states="
           << f_F[ t ].size() << "):";
 for( std::size_t i = 0 ; i < f_F[ t ].size() ; ++i ) {
  const PQFun & F = f_F[ t ][ i ];
  double v = TUEDPINF;
  for( const auto & pc : F )
   if( ( p >= pc.left - 1e-9 ) && ( p <= pc.right + 1e-9 ) ) {
    v = eval_piece( pc , p );
    break;
    }
  std::cerr << "\n  tau=" << f_tau[ t ][ i ] << " dom=["
            << F.front().left << "," << F.back().right << "] npc="
            << F.size() << " v(p)="
            << ( v >= TUEDPINF ? std::string( "INF/uncovered" )
                               : std::to_string( v ) );
  }
 std::cerr << "\n";
 }

/*--------------------------------------------------------------------------*/
/*------------------------- CORE DP ----------------------------------------*/
/*--------------------------------------------------------------------------*/

// Forward dynamic programming on the hybrid "multi-layer ON / single-layer
// OFF" state space. At each time t we maintain:
//
//  - a *sparse* list of surviving F^tau_t value functions (parallel
//    vectors f_F[t], f_tau[t], f_on[t] sorted by tau); these are the
//    ON-side states (on^tau, t) not yet proven irrelevant by RRF+;
//  - two OFF-side scalars c_off_ready[t] (ready to restart) and
//    c_off_any[t] (regardless of duration), plus their back-pointers
//    f_ready_pred[t] / f_any_pred[t];
//  - v_shutdown[t] and (v_shutdown_tau[t], v_shutdown_p[t]): the cost
//    and witness of arriving at the "long shutdown arc" at end of t.
//
// At t = 0 the relevant states are seeded from the initial condition of
// the unit (init_up_down_time): either a single F^{init+1}_0 piece if the
// unit was already on, or a restart F^1_0 plus the initial off-trail
// contribution to c_off_ready[0] / c_off_any[0] if it was off long enough.
//
// For t >= 1 the recurrence is:
//  - c_off_any[t]   = min( c_off_any[t-1], v_shutdown[t-1] )
//  - c_off_ready[t] = min( c_off_ready[t-1],
//                          v_shutdown[t - mdt] if t >= mdt,
//                          0 if init_up_down_time <= 0 and
//                              |init_up_down_time| + t + 1 >= mdt )
//  - new F^1_t            built from c_off_ready[t-1] (restart arc);
//  - new F^{tau+1}_t      built from each surviving F^tau_{t-1} via
//                         sliding_min + add_quadratic + clamp_domain;
//  - prune the new list with RRF+ (pairwise domination check);
//  - finally v_shutdown[t] is the min over surviving tau >= mut of
//    F^tau_t(p) restricted to [P, SD_{t+1}].
//
// The RRF+ pruning is the key to speed: without it the number of F^tau
// entries grows ~ n per time step, making the DP O(n^2) overall; with it
// the surviving set stays O(1) in practice (Wuijts reports k <= 5 on
// the benchmark instances), bringing the DP down to O(n).

void ThermalUnitExtDPSolver::run_DP( void )
{
 if( time_horizon == 0 ) {
  f_best_cost = 0;
  stage = dp_OK;
  return;
  }

 // read the current fixed status of the Variable: the commitment fixings
 // kill the incompatible DP states below, any other fixing makes
 // load_fixings() throw
 load_fixings();
 const auto fixed_on = [ & ]( Index t ) {
  return( f_has_fixings && ( nxt_on[ t ] == t ) );
  };
 const auto fixed_off = [ & ]( Index t ) {
  return( f_has_fixings && ( nxt_off[ t ] == t ) );
  };

 // initialise all the per-time-step state vectors to empty / +INF.
 // recycle the previous solve's PQFun buffers into the pool first, then
 // empty the per-step lists keeping their outer-vector capacity
 for( auto & ft : f_F )
  for( auto & f : ft )
   pool_give( f );
 f_F  .resize( time_horizon );
 f_tau.resize( time_horizon );
 f_on .resize( time_horizon );
 for( auto & v : f_F   ) v.clear();
 for( auto & v : f_tau ) v.clear();
 for( auto & v : f_on  ) v.clear();
 c_off_ready.assign( time_horizon , TUEDPINF );
 c_off_any  .assign( time_horizon , TUEDPINF );
 v_shutdown .assign( time_horizon , TUEDPINF );
 v_shutdown_tau.assign( time_horizon , 0 );
 v_shutdown_p  .assign( time_horizon , 0.0 );
 f_ready_pred.assign( time_horizon , -1 );
 f_any_pred  .assign( time_horizon , -1 );

 // mut / mdt: at least 1, since shut-down and re-start cannot happen in
 // the same instant regardless of the values coming from the Block
 const Index n = time_horizon;
 const Index mut = std::max( min_up_time   , Index( 1 ) );
 const Index mdt = std::max( min_down_time , Index( 1 ) );

 // whether the unit (if on at the initial state) is still inside a start-up
 // or shut-down trajectory at t = 0 and so cannot have been shut down for
 // free at the end of t = -1: it must stay on for at least t = 0. This
 // forbids both the off-at-0 seeding *and* the "free pre-horizon shutdown"
 // init-ready term of the off-side recurrence. See the t = 0 seeding below.
 const bool startup_in_progress =
  has_ramp_up && ( initial_power < min_power[ 0 ] - 1e-9 );
 const bool shutdown_in_progress =
  has_ramp_down && ( initial_power > bound_down[ 0 ] + 1e-9 );

 // precompute the per-period reserve discount g_t(p): the effective cost
 // the DP minimises is f_t + g_t, so g_t is added wherever f_t is. Empty
 // (no-op) unless some reserve price is negative, so the energy-only
 // path is unchanged. The interior discount (cap = max_power) is used at
 // every on->on step, and the start-up discount (cap = bound_on) at the
 // off->on step: the reserve band at a start-up period must fit under
 // the start-up cap bound_on, not the full max_power (the boundary
 // correction). The shut-down variant (cap = bound_down) is built on the
 // fly in compute_v_shutdown, where the actual shut-down cap is known.
 std::vector< PQFun > eff_disc( n ) , eff_disc_su( n );
 for( Index t = 0 ; t < n ; ++t ) {
  eff_disc   [ t ] = build_reserve_discount( t , max_power[ t ] );
  eff_disc_su[ t ] = build_reserve_discount( t , bound_on  [ t ] );
  }

 // reactive on-increment: when the reactive box is commitment-gated
 // ([Qmin_off,Qmax_off] widened by [Qmin_on,Qmax_on] while on), being on
 // earns, on top of the off-state reward accumulated into Q_star
 // (computed on the off box, unchanged), an extra
 // Delta_t = r_on_t - r_off_t per on-period. It is a constant, folded
 // into const_term at each on-period site below. Zero (empty vectors)
 // recovers the plain separable reactive constant.
 std::vector< double > reactive_delta( n , 0.0 );
 if( ( ! reactive_linear_term.empty() ) &&
     ( ( ! reactive_min_on.empty() ) || ( ! reactive_max_on.empty() ) ) )
  for( Index t = 0 ; t < n ; ++t ) {
   const double c = reactive_linear_term[ t ];
   const double lo_off = reactive_min[ t ] , hi_off = reactive_max[ t ];
   const double lo_on = lo_off +
    ( reactive_min_on.empty() ? 0.0 : reactive_min_on[ t ] );
   const double hi_on = hi_off +
    ( reactive_max_on.empty() ? 0.0 : reactive_max_on[ t ] );
   reactive_delta[ t ] = std::min( c * lo_on , c * hi_on ) -
                         std::min( c * lo_off , c * hi_off );
   }

 // RRF+ pairwise-domination pruning (Wuijts et al. 2021, Prop. 6.1): a
 // F^tau_t with tau >= mut is "irrelevant" if it is pointwise dominated
 // by some F^tau'_t with tau' >= mut; if so it cannot appear in any
 // optimal schedule, since its downstream F^{tau+k}_{t+k} is dominated
 // too. For tau < mut the function is kept (it is the unique path
 // through that run-length signature).
 auto prune_RRF_plus = [ & ]( std::vector< PQFun > & v_F ,
                              std::vector< Index > & v_tau ,
                              std::vector< OnSlot > & v_on ) {
  const std::size_t sz = v_F.size();
  if( sz <= 1 )
   return;
#if TUEDPS_PROFILE
  // diagnostic: keep ALL states (no RRF+ pruning)
  if( std::getenv( "TUEDPS_NOPRUNE" ) )
   return;
#endif
  std::vector< char > keep( sz , 1 );
  for( std::size_t i = 0 ; i < sz ; ++i ) {
   if( ! keep[ i ] ) continue;
   if( v_tau[ i ] < mut ) continue;  // keep all tau < mut
   for( std::size_t j = 0 ; j < sz ; ++j ) {
    if( j == i ) continue;
    if( ! keep[ j ] ) continue;
    if( v_tau[ j ] < mut ) continue;
    double domeps = 0.0;
#if TUEDPS_PROFILE
    if( std::getenv( "TUEDPS_DOMEPS" ) )
     domeps = std::atof( getenv( "TUEDPS_DOMEPS" ) );
#endif
    if( is_dominated_by( v_F[ i ] , v_F[ j ] , domeps ) ) {
     keep[ i ] = 0;
     break;
     }
    }
   }
  // compact in-place, preserving order
  std::size_t out = 0;
  for( std::size_t i = 0 ; i < sz ; ++i )
   if( keep[ i ] ) {
    if( out != i ) {
     v_F  [ out ] = std::move( v_F  [ i ] );
     v_tau[ out ] = v_tau[ i ];
     v_on [ out ] = v_on [ i ];
     }
    ++out;
    }
  v_F  .resize( out );
  v_tau.resize( out );
  v_on .resize( out );
  };

 // helper to compute v_shutdown[t]: the min cost of a run that CLOSES at
 // t (unit off at t+1). The closing period's reserve must fit under the
 // shut-down cap sd_hi, NOT the interior max_power. The g0 (per-period
 // band) AND the corr (ramp-residual band) both depend on that cap, and
 // corr is already baked into f_F[t] with cap=max_power, so a scalar
 // delta cannot fix it. When the cap actually bites the closing
 // transition t-1 -> t is re-run from the predecessor f_F[t-1] with
 // acap=sd_hi (exact); otherwise the plain readout of f_F[t] holds.
 auto compute_v_shutdown = [ & ]( Index t , double sd_hi ) {
  double best = TUEDPINF;
  Index best_tau = 0;
  double best_p = 0;
  const double Plo = min_power[ t ];
  const bool cap_bites = ( ! eff_disc[ t ].empty() ) &&
                         ( sd_hi < max_power[ t ] - 1e-12 );
#if TUEDPS_PROFILE
  // diagnostic: skip the shut-down re-run
  if( std::getenv( "TUEDPS_NOSDFIX" ) ) {
   for( std::size_t i = 0 ; i < f_F[ t ].size() ; ++i ) {
    if( f_tau[ t ][ i ] < mut ) continue;
    auto [ v , p ] = min_over( f_F[ t ][ i ] , Plo , sd_hi );
    if( v < best ) { best = v; best_tau = f_tau[ t ][ i ]; best_p = p; } }
   v_shutdown[ t ] = best;
   v_shutdown_tau[ t ] = best_tau;
   v_shutdown_p[ t ] = best_p;
   return;
   }
#endif
  const bool correct = cap_bites && ( t >= 1 );
  if( ! correct ) {                      // plain readout (+ g0 delta at t==0)
   // at t==0 there is NO on->on transition, hence no corr: the reserve
   // cap bites only the per-period g0, which the scalar
   // delta = g_sd - g_int fixes exactly.
   PQFun delta; const bool use_delta = cap_bites && ( t == 0 );
   if( use_delta ) {
    delta = build_reserve_discount( t , sd_hi );
    PQFun neg = eff_disc[ t ];
    for( auto & pc : neg ) { pc.beta = -pc.beta; pc.gamma = -pc.gamma; }
    add_pwq( delta , neg );               // delta = g_sd - g_int
    }
   for( std::size_t i = 0 ; i < f_F[ t ].size() ; ++i ) {
    if( f_tau[ t ][ i ] < mut ) continue;
    std::pair< double , double > vp;
    if( use_delta && ( f_tau[ t ][ i ] > 1 ) ) {  // tau==1 carries g_su, skip
     PQFun F = f_F[ t ][ i ]; add_pwq( F , delta );
     vp = min_over( F , Plo , sd_hi );
     }
    else vp = min_over( f_F[ t ][ i ] , Plo , sd_hi );
    if( vp.first < best ) { best = vp.first; best_tau = f_tau[ t ][ i ];
                            best_p = vp.second; }
    }
   }
  else {                               // re-run the closing transition capped
   const double ru_prev = delta_ramp_up  [ t - 1 ];
   const double rd_prev = delta_ramp_down[ t - 1 ];
   PQFun Fsd;   // g0 under sd_hi is folded into corr (acap=sd_hi)
   for( std::size_t i = 0 ; i < f_F[ t - 1 ].size() ; ++i ) {
    const Index tau = f_tau[ t - 1 ][ i ] + 1;   // this closing period is on
    if( tau < mut ) continue;                    // run too short to shut down
    Fsd.clear();
    // acap = sd_hi (shut-down band)
    sliding_min_corr( f_F[ t - 1 ][ i ] , ru_prev , rd_prev , Plo ,
                      sd_hi , t , Fsd , sd_hi );
    if( Fsd.empty() ) continue;
    add_quadratic( Fsd , quad_term[ t ] , linear_term[ t ] ,
                   const_term[ t ] + reactive_delta[ t ] );
    auto [ v , p ] = min_over( Fsd , Plo , sd_hi );
    if( v < best ) { best = v; best_tau = tau; best_p = p; }
    }
   // tau == 1 (mut == 1: start-up and shut-down in the same period) is not
   // produced by the t-1 -> t transition; read it off f_F[t] as before. Its
   // band is under bound_on already; the extra min(bound_on,sd_hi) tightening
   // is a rare (mut==1) refinement left for a follow-up.
   if( mut <= 1 )
    for( std::size_t i = 0 ; i < f_F[ t ].size() ; ++i ) {
     if( f_tau[ t ][ i ] != 1 ) continue;
     auto [ v , p ] = min_over( f_F[ t ][ i ] , Plo , sd_hi );
     if( v < best ) { best = v; best_tau = 1; best_p = p; }
     }
   }
  v_shutdown[ t ] = best;
  v_shutdown_tau[ t ] = best_tau;
  v_shutdown_p  [ t ] = best_p;
  };

 // -------------------------------------------------------------- t = 0 --
 // Seed the DP state from the initial condition of the unit.

 if( init_up_down_time > 0 ) {
  // the unit was on for init_up_down_time instants strictly before t = 0
  // and is still on at t = 0; the current on-run therefore has length
  // tau0 = init_up_down_time + 1 at t = 0. The power p_0 is constrained
  // to respect the ramp bounds around initial_power, giving a restricted
  // domain [lo, hi] for F^tau0_0(p_0) = f_0(p_0) on that domain.
  Index tau0 = Index( init_up_down_time + 1 );
  double lo = std::max( min_power[ 0 ] ,
                        initial_power - delta_ramp_down[ 0 ] );
  double hi = std::min( max_power[ 0 ] ,
                        initial_power + delta_ramp_up[ 0 ] );
  if( ( ! fixed_off( 0 ) ) && ( lo < hi + 1e-12 ) ) {
   PQFun F;
   F.push_back( { quad_term[ 0 ] , linear_term[ 0 ] ,
                  const_term[ 0 ] + reactive_delta[ 0 ] , lo , hi } );
   add_pwq( F , eff_disc[ 0 ] );
   auto [ v , p ] = min_over( F , lo , hi );
   f_F  [ 0 ].push_back( std::move( F ) );
   f_tau[ 0 ].push_back( tau0 );
   f_on [ 0 ].push_back( { v , p } );
   }
  // when init_up_down_time >= min_up_time the min-up-time is already
  // satisfied at t = 0, and the unit could equivalently have been shut
  // down at the end of t = -1 (a pre-horizon decision, hence "free"
  // within the horizon) and be off at t = 0. We must seed the OFF-side
  // scalars accordingly, otherwise the DP misses any schedule that
  // starts off and pays no cost: with high linear costs this is
  // typically the cheapest option, and overlooking it makes f_best_cost
  // too large (the on-only schedules are forced to pay f_t at t = 0).
  //
  // This "free pre-horizon shutdown" is forbidden, however, when the unit
  // could not actually have shut down at the end of t = -1, i.e. when it is
  // still inside a start-up or a shut-down trajectory at t = 0:
  //  - START-UP: initial_power < min_power[0] (the unit is below minimum
  //    power, ramping up) *and* an actual ramp-up limit is in force
  //    (delta_ramp_up[0] < max_power[0]); ThermalUnitBlock then requires
  //    initial_power + delta_ramp_up[0] >= min_power[0] (RampUpConstraints)
  //    and the MILP keeps the unit on until it reaches min power.
  //  - SHUT-DOWN: initial_power > shut_down_limit[0] (the last on-power was
  //    above the shut-down cap), so the unit must stay on and ramp down to
  //    the cap before it can be switched off.
  // In both cases the MILP keeps the unit on for at least t = 0; the off-at-0
  // seeding is then infeasible and must be skipped (the in-horizon shut-down
  // path through v_shutdown handles the trajectory, possibly multi-period).
  if( ( Index( init_up_down_time ) >= min_up_time ) &&
      ( ! startup_in_progress ) && ( ! shutdown_in_progress ) &&
      ( ! fixed_on( 0 ) ) ) {
   c_off_any[ 0 ] = 0;
   f_any_pred[ 0 ] = -1;
   // ready at t = 0 means the off run ending at t = 0 (just the single
   // instant t = 0, since the shutdown happened at end of t = -1) is at
   // least mdt long; this only happens when mdt <= 1
   if( Index( 1 ) >= min_down_time ) {
    c_off_ready[ 0 ] = 0;
    f_ready_pred[ 0 ] = -1;
    }
   }
  }
 else {
  // init_up_down_time <= 0: the unit has been off for |init| instants
  // strictly before t = 0. We can immediately *restart* at t = 0 if the
  // off period satisfies mdt, i.e. iff |init| >= mdt.
  bool can_restart_t0 = ( Index( - init_up_down_time ) >= mdt );
  if( can_restart_t0 && ( ! fixed_off( 0 ) ) ) {
   // F^1_0(p) = f_0(p) + SUC[0] on [P, min(Pbar, SU)]: the restart arc
   // pays the start-up cost and constrains p by the start-up ramp limit
   double lo = min_power[ 0 ];
   double hi = std::min( max_power[ 0 ] , bound_on[ 0 ] );
   if( lo < hi + 1e-12 ) {
    double suc = startup_costs.empty() ? 0.0 : startup_costs[ 0 ];
    PQFun F;
    F.push_back( { quad_term[ 0 ] , linear_term[ 0 ] ,
                   const_term[ 0 ] + suc + reactive_delta[ 0 ] , lo , hi } );
    add_pwq( F , eff_disc_su[ 0 ] );  // start-up: reserve band under bound_on
    auto [ v , p ] = min_over( F , lo , hi );
    f_F  [ 0 ].push_back( std::move( F ) );
    f_tau[ 0 ].push_back( 1 );
    f_on [ 0 ].push_back( { v , p } );
    }
   }
  // "unit off at t = 0, any duration": cost 0 from the initial off state;
  // f_any_pred[0] = -1 signals that no in-horizon shutdown happened
  if( ! fixed_on( 0 ) ) {
   c_off_any[ 0 ] = 0;
   f_any_pred[ 0 ] = -1;
   // "unit off at t = 0 AND ready": the off trail before t = 0 counts
   // |init| instants and t = 0 itself adds one more, so the condition is
   // |init| + 1 >= mdt (equivalently, the unit is ready right at t = 0)
   if( Index( - init_up_down_time ) + 1 >= mdt ) {
    c_off_ready[ 0 ] = 0;
    f_ready_pred[ 0 ] = -1;
    }
   }
  }

 // prune at t = 0 (with a single seeded entry this is a no-op, but for
 // completeness we apply it to every time step)
 prune_RRF_plus( f_F[ 0 ] , f_tau[ 0 ] , f_on[ 0 ] );

 // v_shutdown[0]: only meaningful if there is at least one more instant
 // in the horizon after t = 0 (so that being off at t = 1 is possible)
 if( n > 1 )
  compute_v_shutdown( 0 , bound_down[ 1 ] );

#if TUEDPS_PROFILE
 // CLEAN-REFERENCE on-run (TUEDPS_CLEANREF): propagate a pure all-on run
 // two ways from the SAME clean input each step: (1) the solver's
 // sliding_min_corr (mpQP if TUEDPS_MPQP), (2) a dense-grid BRUTE of the
 // identical formula min_q[F(q)+reserve_reward(min(A,B))]. Reveals the
 // true per-transition drift (breaks the out=Gmin self-consistency).
 if( std::getenv( "TUEDPS_CLEANREF" ) ) {
  double glo = 1e300 , ghi = -1e300;
  for( Index t = 0 ; t < n ; ++t ) { glo = std::min( glo , min_power[ t ] );
                                     ghi = std::max( ghi , max_power[ t ] ); }
  const int N = std::getenv( "TUEDPS_CRN" )
                ? std::atoi( getenv( "TUEDPS_CRN" ) ) : 2000;
  const int NQ = std::getenv( "TUEDPS_CRNQ" )
                 ? std::atoi( getenv( "TUEDPS_CRNQ" ) ) : 1500;
  std::vector< double > cf( N + 1 );
  auto gp = [ & ]( int i ){ return glo + ( ghi - glo ) * i / N; };
  for( int i = 0 ; i <= N ; ++i ) { double p = gp( i );
   cf[ i ] = ( p >= min_power[ 0 ] - 1e-9 && p <= max_power[ 0 ] + 1e-9 )
           ? ( quad_term[ 0 ] * p + linear_term[ 0 ] ) * p + const_term[ 0 ]
           : TUEDPINF; }
  auto cfEval = [ & ]( double q )->double {
   if( q < glo - 1e-9 || q > ghi + 1e-9 ) return TUEDPINF;
   double x = ( q - glo ) / ( ghi - glo ) * N; int i = ( int )x;
   if( i < 0 ) i = 0;
   if( i >= N ) i = N - 1;
   double a = cf[ i ] , b = cf[ i + 1 ];
   if( a>=TUEDPINF || b>=TUEDPINF ) return std::min( a , b );
   return a + ( b - a ) * ( x - i ); };
  auto cfToPQ = [ & ]()->PQFun { PQFun G;
   for( int i = 0 ; i < N ; ++i ) {
    if( cf[ i ]>=TUEDPINF || cf[ i + 1 ]>=TUEDPINF ) continue;
    double pa = gp( i ) , pb = gp( i + 1 ) ,
           sl = ( cf[ i + 1 ] - cf[ i ] ) / ( pb - pa );
    G.push_back( { 0.0 , sl , cf[ i ] - sl * pa , pa , pb } ); }
   return G; };
  const double crstop = std::getenv( "TUEDPS_CRSTOP" )
                        ? std::atof( getenv( "TUEDPS_CRSTOP" ) ) : 1e30;
  for( Index t = 1 ; t < n ; ++t ) {
   PQFun Fprev = cfToPQ(); if( Fprev.empty() ) break;
   const double dru = delta_ramp_up[ t - 1 ] , drd = delta_ramp_down[ t - 1 ];
   PQFun outM;
   sliding_min_corr( Fprev , dru , drd , min_power[ t ] ,
                     max_power[ t ] , t , outM );
   if( ! outM.empty() )
    add_quadratic( outM , quad_term[ t ] , linear_term[ t ] ,
                   const_term[ t ] );
   std::vector< double > nf( N + 1 , TUEDPINF );
   double maxdev = 0 , mx = 0;
   const double fl = Fprev.front().left , fr = Fprev.back().right;
   for( int i = 0 ; i <= N ; ++i ) { double p = gp( i );
    if( p < min_power[ t ] - 1e-9 || p > max_power[ t ] + 1e-9 ) continue;
    const double A = std::min( p - min_power[ t ] , max_power[ t ] - p );
    const double qlo = std::max( p - dru , fl ) ,
                 qhi = std::min( p + drd , fr );
    if( qlo > qhi + 1e-12 ) continue;
    double best = TUEDPINF;
    for( int j = 0 ; j <= NQ ; ++j ) {
     double q = qlo + ( qhi - qlo ) * j / NQ;
     const double fq = cfEval( q ); if( fq >= TUEDPINF ) continue;
     const double d = p - q;
     double B = std::min( dru - d , drd + d ); if( B<0 )B = 0;
     const double H = ( A<=0 ) ? 0 : std::min( A , B );
     const double v = fq + reserve_reward( t , p , H );
     if( v < best ) best = v; }
    if( best < TUEDPINF ) {
     nf[ i ] = best + ( quad_term[ t ] * p + linear_term[ t ] ) * p +
               const_term[ t ];
     const double vm = eval( outM , p );
     if( vm < TUEDPINF ) { const double dd = std::abs( vm - nf[ i ] );
      if( dd > maxdev ) { maxdev = dd; mx = p; } } } }
   std::cerr.precision( 10 );
   std::cerr << "CLEANREF t=" << t << " maxdev(mpQP-brute)=" << maxdev
             << " at p=" << mx << " |Fprev|=" << Fprev.size()
             << " |outM|=" << outM.size() << "\n";
   cf = nf;   // propagate the CLEAN (brute) F
   if( maxdev > crstop ) {
    std::cerr << "CLEANREF STOP: mpQP diverges from brute at t=" << t
              << " by " << maxdev << " at p=" << mx << "\n";
    break; }
   }
  std::exit( 0 );
  }
#endif

 // ----------------------------------------------------- main loop: t >=1 -
 for( Index t = 1 ; t < n ; ++t ) {

  // -- OFF-side updates ------------------------------------------------ //
  // c_off_any[t] : min of (stay off from t-1) and (just shut down at end
  // of t-1). The freshly-shutdown path uses v_shutdown[t-1], while the
  // stay path carries over whatever last shutdown contributed to
  // c_off_any[t-1] (through f_any_pred[t-1]).
  double stay_any = c_off_any[ t - 1 ];
  double fresh_any = v_shutdown[ t - 1 ];
  if( fresh_any < stay_any ) {
   c_off_any [ t ] = fresh_any;
   f_any_pred[ t ] = int( t - 1 );
   }
  else {
   c_off_any [ t ] = stay_any;
   f_any_pred[ t ] = f_any_pred[ t - 1 ];
   }

  // c_off_ready[t] : min of three contributions,
  //  - stay ready from t-1;
  //  - "long shutdown arc": a shutdown decided at the end of time
  //    (t - mdt) reaches c_off_ready[t] through the long arc spanning
  //    mdt instants (the min-down-time window);
  //  - initial-off trail: when init_up_down_time <= 0 and
  //    |init_up_down_time| + t + 1 >= mdt the unit has been off for at
  //    least mdt consecutive instants ending at t even without any
  //    in-horizon shutdown, with zero accumulated cost.
  double stay_ready = c_off_ready[ t - 1 ];
  double fresh_ready = TUEDPINF;
  int fresh_h = -1;
  // the "long shutdown arc" spans the off instants [ t - mdt + 1 , t ]
  // directly, so it must be checked against the fixed-ON instants (which
  // the per-instant kill below cannot intercept)
  if( ( t >= mdt ) &&
      ( ( ! f_has_fixings ) || ( nxt_on[ t - mdt + 1 ] > t ) ) ) {
   fresh_ready = v_shutdown[ t - mdt ];
   fresh_h = int( t - mdt );
   }
  // initial-off trail covers t when the off run by end of t reaches mdt:
  //  - if init_up_down_time <= 0, the off run includes the |init| pre-
  //    horizon instants plus the in-horizon ones up to and including t,
  //    i.e. -init_ud + t + 1 instants;
  //  - if init_up_down_time > 0 but >= mut, the unit could have shut
  //    down at end of t = -1 (a free pre-horizon decision since mut is
  //    already satisfied), so the in-horizon off run is t + 1 instants
  double init_ready = TUEDPINF;
  // the initial-off trail spans the off instants [ 0 , t ] directly, so it
  // must be checked against the fixed-ON instants too
  if( ( ! f_has_fixings ) || ( nxt_on[ 0 ] > t ) ) {
   if( ( init_up_down_time <= 0 ) &&
       ( Index( - init_up_down_time ) + t + 1 >= mdt ) )
    init_ready = 0;
   else if( ( init_up_down_time > 0 ) &&
	    ( Index( init_up_down_time ) >= min_up_time ) &&
	    ( t + 1 >= mdt ) &&
	    ( ! startup_in_progress ) && ( ! shutdown_in_progress ) )
    // the "free pre-horizon shutdown" is unavailable while the unit is still
    // inside its initial start-up / shut-down trajectory (see t = 0 seeding)
    init_ready = 0;
   }
  double best_ready = std::min( { stay_ready , fresh_ready , init_ready } );
  c_off_ready[ t ] = best_ready;
  // record the source so backtracking knows where the "ready" came from;
  // the init-ready case has no in-horizon shutdown and uses -1
  if( best_ready == init_ready && init_ready < TUEDPINF )
   f_ready_pred[ t ] = -1;
  else if( best_ready == fresh_ready && fresh_ready < TUEDPINF )
   f_ready_pred[ t ] = fresh_h;
  else
   f_ready_pred[ t ] = f_ready_pred[ t - 1 ];

  if( fixed_on( t ) ) {   // the unit cannot be off at t: kill both OFF
   c_off_any   [ t ] = TUEDPINF;   // states (the "stay off" paths through
   f_any_pred  [ t ] = -1;         // t die here, the arcs jumping over t
   c_off_ready [ t ] = TUEDPINF;   // have been checked above)
   f_ready_pred[ t ] = -1;
   }

  // -- ON-side updates ------------------------------------------------- //
  // build a fresh sparse list for time t from (a) a possible restart
  // tau=1 and (b) the increments of each surviving entry at t-1
  double Plo = min_power[ t ];
  double Phi = max_power[ t ];
  double ru_prev = delta_ramp_up  [ t - 1 ];
  double rd_prev = delta_ramp_down[ t - 1 ];

  // per-step build buffers: reused member scratch (cleared, capacity kept)
  m_new_F .clear();
  m_new_tau.clear();
  m_new_on.clear();
  m_new_F  .reserve( f_F[ t - 1 ].size() + 1 );
  m_new_tau.reserve( f_F[ t - 1 ].size() + 1 );
  m_new_on .reserve( f_F[ t - 1 ].size() + 1 );

  // (a) tau = 1 entry (restart arc): F^1_t(p) = f_t(p) + SUC[t]
  //     + c_off_ready[t-1], defined for p in [P, min(Pbar, SU)]. Only
  //     built if c_off_ready[t-1] is finite (otherwise the unit cannot
  //     legally restart at t).
  if( ( ! fixed_off( t ) ) && ( c_off_ready[ t - 1 ] < TUEDPINF ) ) {
   double lo = Plo;
   double hi = std::min( Phi , bound_on[ t ] );
   if( lo < hi + 1e-12 ) {
    double suc = startup_costs.empty() ? 0.0 : startup_costs[ t ];
    PQFun F = pool_take();
    F.push_back( { quad_term[ t ] , linear_term[ t ] ,
                   const_term[ t ] + suc + c_off_ready[ t - 1 ] +
                   reactive_delta[ t ] ,
                   lo , hi } );
    add_pwq( F , eff_disc_su[ t ] );  // start-up: reserve band under bound_on
    auto [ v , p ] = min_over( F , lo , hi );
    m_new_F  .push_back( std::move( F ) );
    m_new_tau.push_back( 1 );
    m_new_on .push_back( { v , p } );
    }
   }

  // (b) tau > 1 entries: from each surviving entry of f_F[t-1] apply
  //     sliding_min (ramp window from t-1 to t) to get the value
  //     function restricted to [Plo, Phi], then add f_t(p). The output
  //     tau is the previous tau + 1 (continuing the on-run). No ON state
  //     exists at t if the commitment there is fixed OFF.
  for( std::size_t i = 0 ;
       ( ! fixed_off( t ) ) && ( i < f_F[ t - 1 ].size() ) ; ++i ) {
   PQFun F = pool_take();
   // interior on->on step with the residual-ramp reserve penalty corr folded
   // into the ramp-window minimisation (falls back to the plain sliding_min
   // when no reserve is rewarded at t, keeping energy-only bit-identical)
   sliding_min_corr( f_F[ t - 1 ][ i ] , ru_prev , rd_prev , Plo , Phi ,
                     t , F );
   if( F.empty() ) {           // intersection with [Plo, Phi] empty
    pool_give( F );            // hand the unused buffer back to the pool
    continue;
    }
   add_quadratic( F , quad_term[ t ] , linear_term[ t ] ,
                  const_term[ t ] + reactive_delta[ t ] );
   // g^0 for period t is folded into corr (=g(p,min(A,B))) inside
   // sliding_min_corr, so no separate eff_disc add here (keeps F
   // exactly convex).
#if TUEDPS_PROFILE
   // DIAGNOSTIC (env TUEDPS_SIMP=eps): lossy compressor of the convex
   // F. Greedily fuses a maximal run of adjacent pieces into one
   // quadratic interpolating the run endpoints+midpoint, accepted only
   // if within eps*(1+|f|) at 11 samples and still convex. Measures the
   // accuracy/speed tradeoff of bounding |F|.
   { static const char * senv = std::getenv( "TUEDPS_SIMP" );
     if( senv && ( F.size() > 2 ) ) {
      const double eps = std::atof( senv );
      auto ev = []( const PieceQuad & pc , double x )
       { return( ( pc.alfa * x + pc.beta ) * x + pc.gamma ); };
      auto evF = [ & ]( std::size_t lo , std::size_t hi , double x ) {
       for( std::size_t k = lo ; k <= hi ; ++k )
        if( ( x >= F[ k ].left - 1e-12 ) && ( x <= F[ k ].right + 1e-12 ) )
         return( ev( F[ k ] , x ) );
       return( ev( F[ hi ] , x ) ); };
      PQFun sout; sout.reserve( F.size() );
      std::size_t si = 0;
      while( si < F.size() ) {
       std::size_t sj = si; PieceQuad cur = F[ si ];
       while( sj + 1 < F.size() ) {
        const double L = F[ si ].left , R = F[ sj + 1 ].right ,
                     Mm = 0.5 * ( L + R );
        const double fL = evF( si , sj + 1 , L ) ,
                     fM = evF( si , sj + 1 , Mm ) ,
                     fR = evF( si , sj + 1 , R );
        const double aa = ( ( fR - fM ) / ( R - Mm ) -
                            ( fM - fL ) / ( Mm - L ) ) / ( R - L );
        const double bb = ( fM - fL ) / ( Mm - L ) - aa * ( L + Mm );
        const double cc = fL - aa * L * L - bb * L;
        if( aa < -1e-12 ) break;
        bool ok = true;
        for( int s = 0 ; s <= 10 && ok ; ++s ) {
         const double x = L + ( R - L ) * s / 10.0 ,
                      fx = evF( si , sj + 1 , x );
         if( std::abs( ( aa * x + bb ) * x + cc - fx ) >
             eps * ( 1 + std::abs( fx ) ) )
          ok = false;
         }
        if( ! ok ) break;
        cur = { aa , bb , cc , L , R }; ++sj;
        }
       sout.push_back( cur ); si = sj + 1;
       }
      F.swap( sout );
      } }
#endif
   auto [ v , p ] = min_over( F , Plo , Phi );
   m_new_F  .push_back( std::move( F ) );
   m_new_tau.push_back( f_tau[ t - 1 ][ i ] + 1 );
   m_new_on .push_back( { v , p } );
   }

  // the new_tau vector is already sorted ascending: the restart entry
  // (tau = 1) was pushed first, and each subsequent entry has
  // tau = (old tau) + 1 >= 2 with the old taus appearing in ascending
  // order; since the old list had tau >= 1, we have new tau >= 2 and
  // monotonicity is preserved.

  // RRF+ pruning: remove functions with tau >= mut that are pointwise
  // dominated by some other function with tau >= mut in the new list.
  prune_RRF_plus( m_new_F , m_new_tau , m_new_on );

  // hand the per-step lists over to f_F[t] (emptied at reset); the build
  // buffers are left empty and regrow on the next step
  f_F  [ t ] = std::move( m_new_F );
  f_tau[ t ] = std::move( m_new_tau );
  f_on [ t ] = std::move( m_new_on );

#if TUEDPS_PROFILE
  // GROWTH TRACE (env TUEDPS_TRACE): for the max-|F| state at selected
  // periods, report |F| and how many pieces survive a value-merge at
  // 1e-10/1e-8/1e-6, distinguishing genuine cost-to-go complexity from
  // representation bloat.
  if( std::getenv( "TUEDPS_TRACE" ) ) {
   std::size_t bi = 0 , bm = 0;
   for( std::size_t i = 0 ; i < f_F[ t ].size() ; ++i )
    if( f_F[ t ][ i ].size() > bm ) { bm = f_F[ t ][ i ].size(); bi = i; }
   static std::size_t gmax = 0 , gmaxt = 0;
   if( bm > gmax ) { gmax = bm; gmaxt = t; }
   if( t + 1 == n )
    std::cerr << "TUEDPS_TRACE GLOBAL-MAX |F|=" << gmax << " at t=" << gmaxt
              << std::endl;
   if( ( bm > 100 ) || ( t % 200 == 0 ) || ( t + 1 == n ) ) {
    const PQFun & Fb = f_F[ t ][ bi ];
    auto ev = []( const PieceQuad & pc , double x )
     { return( ( pc.alfa * x + pc.beta ) * x + pc.gamma ); };
    auto merged = [ & ]( double eps ) -> std::size_t {
     if( Fb.size() <= 1 ) return( Fb.size() );
     std::size_t cnt = 1; PieceQuad b = Fb[ 0 ];
     for( std::size_t i = 1 ; i < Fb.size() ; ++i ) {
      const PieceQuad & c = Fb[ i ];
      const double x1 = c.right , x2 = 0.5 * ( c.left + c.right );
      const double e1 = std::abs( ev( b , x1 ) - ev( c , x1 ) ) ,
                   e2 = std::abs( ev( b , x2 ) - ev( c , x2 ) );
      const double sc = std::max( 1.0 , std::abs( ev( c , x2 ) ) );
      if( ( e1 <= eps * sc ) && ( e2 <= eps * sc ) ) b.right = c.right;
      else { ++cnt; b = c; } }
     return( cnt ); };
    std::cerr << "TUEDPS_TRACE t=" << t << " states=" << f_F[ t ].size()
              << " tau=" << f_tau[ t ][ bi ] << " |F|=" << bm
              << " merge1e-10=" << merged( 1e-10 )
              << " merge1e-8=" << merged( 1e-8 )
              << " merge1e-6=" << merged( 1e-6 ) << std::endl;
    }
   }
#endif

  // v_shutdown[t] is meaningful only if there is at least one more
  // instant in the horizon after t, since an in-horizon off period
  // triggered by shutdown at end of t starts at t + 1
  if( t < n - 1 )
   compute_v_shutdown( t , bound_down[ t + 1 ] );

#if TUEDPS_PROFILE
  // per-t objective trajectory (diff mpQP vs param)
  if( std::getenv( "TUEDPS_TTRACE" ) ) {
   double bon = TUEDPINF; std::size_t nf = 0;
   for( const auto & s : f_on[ t ] ) bon = std::min( bon , s.min_val );
   for( const auto & Fv : f_F[ t ] ) nf += Fv.size();
   std::cerr.precision( 12 );
   std::cerr << "TT t=" << t << " bon=" << bon << " coff=" << c_off_any[ t ]
             << " vsd=" << v_shutdown[ t ] << " states=" << f_F[ t ].size()
             << " sumF=" << nf << "\n";
   }
#endif

  }  // end( for( t ) )

 // finalise best cost at the end of the horizon. Two options:
 //  - terminate *on*: the optimal cost is the minimum of F^tau_{n-1}(p)
 //    over all surviving tau and p (no SD constraint on p: the unit is
 //    free to remain on after the horizon);
 //  - terminate *off*: c_off_any[n-1] aggregates both "stayed off since
 //    before the horizon" and "shut down at some point during the
 //    horizon, then stayed off" without further constraints.
 double best_on = TUEDPINF;
 for( const auto & slot : f_on[ n - 1 ] )
  if( slot.min_val < best_on )
   best_on = slot.min_val;
 double best_off = c_off_any[ n - 1 ];

 f_best_cost = std::min( best_on , best_off );

 // reactive power contribution (AC instances): q[t] in [Qmin,Qmax] is
 // separable from the commitment/active-power DP and carries the
 // dualized linear cost reactive_linear_term[t]. The optimal q*[t] is
 // the box endpoint minimising c*q, contributing min(c*Qmin, c*Qmax);
 // being path-independent it is computed once. It is part of the
 // operational cost (the unit holds q only if it exists), hence folded
 // in BEFORE the design decision below; for a non-built unit
 // get_var_solution() reports q*[t] = 0.
 double Q_star = 0;
 if( f_best_cost < TUEDPINF && ! reactive_linear_term.empty() )
  for( Index t = 0 ; t < time_horizon ; ++t ) {
   const double c = reactive_linear_term[ t ];
   if( c != 0 )
    Q_star += std::min( c * reactive_min[ t ] , c * reactive_max[ t ] );
   }

 // design (investment) decision: the DP above solved the operational
 // problem as if the unit exists (design == 1), so f_best_cost + Q_star
 // is the optimal operational cost p*, the active-power schedule plus
 // the reactive contribution, both available only if the unit is built.
 // Building costs design_cost on top; it is worth building iff
 // p* + design_cost <= 0, otherwise the unit is not built and
 // contributes nothing (cost 0, zero schedule, zero reactive). A unit
 // carrying an investment cost is always initially off
 // (ThermalUnitBlock forbids InitUpDownTime >= 0 with an investment
 // cost), so the design is never forced on. Generalises to an integer
 // design by the same threshold argument; the continuous case does not
 // apply (binary commitment decisions inside).
 if( has_design ) {
  // a commitment fixed ON (or the design variable fixed to 1) forces the
  // unit to be built regardless of the economics; if it cannot (no feasible
  // schedule, or the design variable is fixed to 0) the problem is
  // unfeasible, since the all-zero "not built" solution violates the fixings
  if( ( ! f_no_build ) && ( f_best_cost < TUEDPINF ) &&
      ( f_must_build || ( f_best_cost + Q_star + design_cost <= 0 ) ) ) {
   design_on = true;
   f_best_cost += Q_star + design_cost;
   }
  else if( f_must_build ) {
   design_on = false;
   f_best_cost = TUEDPINF;   // unfeasible: must be built, but cannot
   }
  else {
   design_on = false;
   f_best_cost = 0;          // not built: the unit is absent
   }
  }
 else
  f_best_cost += Q_star;    // no design: the reactive term is always incurred

 f_solved = ( f_best_cost < TUEDPINF );

 stage = dp_OK;

 }  // end( ThermalUnitExtDPSolver::run_DP )

/*--------------------------------------------------------------------------*/
/*------------------------ BUILDING THE SOLUTION ---------------------------*/
/*--------------------------------------------------------------------------*/

// Reconstruct the optimal schedule by walking the DP backward. The
// output is written to P[] (power per time step) and U[] (commitment per
// time step). The walk alternates two phases:
//
//  - ON-phase (walk_on): given a state (t, tau, p), set P[t]=p, U[t]=1,
//    then step to the predecessor (t-1, tau-1, p_{t-1}). The predecessor
//    power is given by Wuijts eq. (16): p_{t-1} is the point of
//    [p - Delta+, p + Delta-] closest to p_star_{t-1, tau-1}. Iterate
//    until tau reaches 1 (meaning the on-run just started at this t) or
//    t reaches 0 (meaning we hit the initial on state from before the
//    horizon). Returns the time at which the on-run starts.
//
//  - OFF-phase (handled in the main loop): a just-finished on-run ends
//    at t_start = 1 means the preceding off-period at t_start - 1 was
//    caused by some earlier shutdown h, recoverable from f_ready_pred.
//    If h >= 0, we jump to (h, v_shutdown_tau[h], v_shutdown_p[h]) and
//    continue with walk_on. If h == -1, the off period came from the
//    initial state and backtracking ends.
//
// When the horizon terminates *off* (best_off < best_on), the chain is
// bootstrapped from f_any_pred[n-1] instead of an end-of-horizon on-state.

void ThermalUnitExtDPSolver::build_solution( void )
{
 std::fill( P.begin() , P.end() , 0.0 );
 std::fill( U.begin() , U.end() , false );

 // fast paths: nothing to recover from
 if( ! f_solved ) {
  stage = sol_OK;
  return;
  }

 const Index n = time_horizon;
 if( n == 0 ) {
  stage = sol_OK;
  return;
  }

 // binary-search helper: find the index i with f_tau[t][i] == tau in the
 // sorted sparse list, or report "not found" (returns the size of the
 // list, analogous to std::lower_bound's end sentinel)
 auto find_tau = [ & ]( Index t , Index tau ) -> std::size_t {
  const auto & v = f_tau[ t ];
  auto it = std::lower_bound( v.begin() , v.end() , tau );
  if( it != v.end() && *it == tau )
   return( std::size_t( it - v.begin() ) );
  return( v.size() );
  };

 // choose the optimal terminal state: either some (tau, p) at t = n-1
 // (for an on-termination) or c_off_any[n-1] (for an off-termination)
 double best_on = TUEDPINF;
 Index  best_tau = 0;
 double best_p = 0;
 for( std::size_t i = 0 ; i < f_on[ n - 1 ].size() ; ++i )
  if( f_on[ n - 1 ][ i ].min_val < best_on ) {
   best_on = f_on[ n - 1 ][ i ].min_val;
   best_tau = f_tau[ n - 1 ][ i ];
   best_p = f_on[ n - 1 ][ i ].argmin_p;
   }
 double best_off = c_off_any[ n - 1 ];

 // walk_on: steps one on-run back from (t, tau, p), writing P[] / U[] as
 // it goes, and returns the time at which the run started. The caller
 // then decides whether (and where) to resume with an earlier on-run via
 // f_ready_pred[].
 // acap0 is the shut-down reserve cap of the CLOSING period
 // (bound_down[t+1]) when this on-run ends in a shut-down: the
 // transition into the closing period was priced with that capped band,
 // so its predecessor must be recovered with the same cap. Only the
 // first predecessor step is the closing one; interior steps (acap = -1)
 // use max_power as usual.
 auto walk_on = [ & ]( Index t , Index tau , double p ,
                       double acap0 ) -> Index {
  double acap = acap0;
  while( true ) {
   P[ t ] = p;
   U[ t ] = true;
   if( tau == 1 )
    return( t );            // first instant of this on-run
   if( t == 0 )
    return( 0 );            // hit initial on state (init_up_down_time > 0)
   // locate the predecessor (t-1, tau-1) in the sparse list
   std::size_t idx = find_tau( t - 1 , tau - 1 );
   if( idx == f_tau[ t - 1 ].size() )
    // defensive: by RRF+ invariants this never happens (a surviving
    // F^tau_t was built from a then-surviving F^{tau-1}_{t-1}), but if
    // numerical effects prune the predecessor anyway we stop cleanly
    return( t );
   // recover the predecessor power. With no reserve this is Wuijts
   // eq. (16), the projection of argmin F^{tau-1}_{t-1} onto the ramp
   // window; with reserve rewarded the transition minimised
   // F(q)+corr_t(q,p), so the corr-aware argmin must be used instead
   // (reserve_corr_argmin subsumes both and clamps to F's own domain,
   // which already carries the tau-1==1 start-up cap).
   const double ru = delta_ramp_up  [ t - 1 ];
   const double rd = delta_ramp_down[ t - 1 ];
   const double p_prev =
    reserve_corr_argmin( f_F[ t - 1 ][ idx ] , ru , rd , t , p , acap );
   acap = -1.0;             // only the closing step carries the shut-down cap
   t   -= 1;
   tau -= 1;
   p = p_prev;
   }
  };

 // The path backward alternates between on-runs and off-periods. Each
 // on-run either ends at n-1 (best_on case) or just before an off-period
 // (best_off or intermediate case); each off-period was introduced by a
 // shutdown event at some index h, unless it is the initial off interval
 // from before the horizon. The loop below walks these segments until
 // either time 0 is reached, or the initial off state is reached.

 Index on_t;
 Index on_tau;
 double on_p;
 double on_acap;   // closing period's shut-down cap (-1 if on-to-end)
 bool have_on;

 // bootstrap the alternating on/off walk
 if( best_on <= best_off ) {
  // on-termination: the last on-run ends at (n-1, best_tau, best_p) -- no
  // shut-down closes it (the unit is on at the horizon end), so acap = -1.
  on_t = n - 1;
  on_tau = best_tau;
  on_p = best_p;
  on_acap = -1.0;
  have_on = true;
  }
 else {
  // off-termination: the optimum ends off at n-1. The off period covers
  // [f_any_pred[n-1] + 1, n-1] (inclusive). We have nothing to write
  // on that range because P[] and U[] are already 0/false. If
  // f_any_pred[n-1] >= 0 we continue with the last shutdown event;
  // otherwise the unit has been off throughout the horizon and we're done.
  have_on = false;
  int h = f_any_pred[ n - 1 ];
  if( h >= 0 ) {
   on_t = Index( h );
   on_tau = v_shutdown_tau[ h ];
   on_p = v_shutdown_p[ h ];
   on_acap = ( Index( h ) + 1 < n ) ? bound_down[ h + 1 ] : -1.0;
   have_on = true;
   }
  }

 // alternate on-runs and off-gaps until we reach the beginning of the
 // horizon (t_start == 0) or an off-gap that was inherited from before
 // the horizon (f_ready_pred[t_start - 1] == -1)
 while( have_on ) {
  Index t_start = walk_on( on_t , on_tau , on_p , on_acap );

  if( t_start == 0 ) {
   have_on = false;
   break;
   }
  // the on-run ended at t_start, so t_start - 1 was off. Find the
  // shutdown event that produced the ready state at t_start - 1
  int h = f_ready_pred[ t_start - 1 ];
  if( h < 0 ) {
   // the ready came from the initial off trail (pre-horizon); no
   // earlier on-run exists, so backtracking is complete
   have_on = false;
   break;
   }
  // jump to the shutdown event and continue the walk; this on-run closes with
  // a shut-down at on_t (off at on_t+1), so its closing period is capped.
  on_t = Index( h );
  on_tau = v_shutdown_tau[ h ];
  on_p = v_shutdown_p[ h ];
  on_acap = ( Index( h ) + 1 < n ) ? bound_down[ h + 1 ] : -1.0;
  }

 stage = sol_OK;

 }  // end( ThermalUnitExtDPSolver::build_solution )

/*--------------------------------------------------------------------------*/
/*-------------- End File ThermalUnitExtDPSolver.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
