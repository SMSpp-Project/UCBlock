/*--------------------------------------------------------------------------*/
/*--------------------- File ThermalUnitExtDPSolver.cpp --------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of ThermalUnitExtDPSolver: a hybrid DP solver for the
 * single Unit Commitment problem in which the on-side of the graph keeps
 * the multi-layer tau-counter of Wuijts et al. (2021) while the off-side
 * is collapsed to a single layer and the min-down-time constraint is
 * enforced by a "long" arc that jumps M_down instants ahead.
 *
 * The ON layer stores convex piecewise quadratic functions F^tau_t : the
 * sliding minimum and domain clamping follow standard Frangioni-Gentile
 * three-case analysis. The OFF layer stores two scalars per time step
 * (c_off_ready[t], c_off_any[t]) plus the auxiliary sequence v_shutdown[t]
 * (cost of a schedule on-through-t plus the shutdown cost, ready to feed
 * the long shutdown arc).
 *
 * Data loading, Modification handling and Solution writing are patterned
 * after ThermalUnitDPSolver, but implemented independently (no code shared
 * at build time) per user request.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Claude Opus 4.7 \n
 *         Anthropic
 *
 * \copyright &copy; by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "ThermalUnitExtDPSolver.h"

#include "ThermalUnitBlock.h"

#include <algorithm>
#include <cmath>

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

 if( stage < dp_OK )
  run_DP();

 if( stage < sol_OK )
  build_solution();

 unlock();

 assert( stage == sol_OK );
 return( f_best_cost >= TUEDPINF ? kInfeasible : kOK );
 }

/*--------------------------------------------------------------------------*/

// Write the optimal schedule stored in P[] and U[] into the variables of
// the ThermalUnitBlock. This matches ThermalUnitDPSolver's semantics:
//
//  - active_power and commitment variables are set uniformly across the
//    whole horizon from P[] and U[];
//  - start_up[t] is set iff the unit transitions off->on between t-1
//    and t, accounting for the "boundary" cases before t_init by
//    consulting init_up_down_time instead of U[-1];
//  - shut_down[t] symmetrically encodes on->off transitions.
//
// The Block is locked for write access and unlocked before returning; if
// the caller already owns the Block (same id), we skip the locking dance.

void ThermalUnitExtDPSolver::get_var_solution( Configuration * solc )
{
 bool owned = f_Block->is_owned_by( f_id );
 if( ( ! owned ) && ( ! f_Block->lock( f_id ) ) )
  throw( std::runtime_error(
   "ThermalUnitExtDPSolver::get_var_solution: unable to lock the Block." ) );

 auto b = static_cast< ThermalUnitBlock * >( f_Block );

 // active power: one value per time step, directly from the DP output
 if( auto pow_it = b->get_active_power( 0 ) )
  for( Index i = 0 ; i < time_horizon ; )
   ( pow_it++ )->set_value( P[ i++ ] );

 // commitment: binary, 1 iff the unit is on at that step
 if( auto com_it = b->get_commitment( 0 ) )
  for( Index i = 0 ; i < time_horizon ; )
   ( com_it++ )->set_value( U[ i++ ] ? 1 : 0 );

 // start-up indicator: 1 iff the unit just turned on at step i. The
 // boundary at i == 0 uses init_up_down_time to decide whether the
 // previous step was off
 if( auto sup_it = b->get_start_up() ) {
  if( ! t_init )
   ( sup_it++ )->set_value( ( init_up_down_time <= 0 ) && ( U[ 0 ] ? 1 : 0 ) );
  for( Index i = std::max( t_init , Index( 1 ) ) ; i < time_horizon ; ++i )
   ( sup_it++ )->set_value( ( U[ i ] ) && ( ! U[ i - 1 ] ) ? 1 : 0 );
  }

 // shut-down indicator: symmetric to start-up
 if( auto sdn_it = b->get_shut_down() ) {
  if( ! t_init )
   ( sdn_it++ )->set_value( ( init_up_down_time > 0 ) && ( ! U[ 0 ] ) ? 1 : 0 );
  for( Index i = std::max( t_init , Index( 1 ) ) ; i < time_horizon ; ++i )
   ( sdn_it++ )->set_value( ( ! U[ i ] ) && ( U[ i - 1 ] ? 1 : 0 ) );
  }

 if( ! owned )
  f_Block->unlock( f_id );

 }  // end( ThermalUnitExtDPSolver::get_var_solution )

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
 bool owned = f_Block->is_owned_by( f_id );
 if( ( ! owned ) && ( ! f_Block->read_lock() ) )
  throw( std::runtime_error(
   "ThermalUnitExtDPSolver::load_parameters: unable to lock the Block." ) );

 auto b = static_cast< ThermalUnitBlock * >( f_Block );

 // guard against unsupported features: the DP does not model reserves
 if( ! b->get_primary_rho().empty() )
  throw( std::invalid_argument( "ThermalUnitExtDPSolver::load_parameters: "
                                "primary reserve not supported." ) );

 if( ! b->get_secondary_rho().empty() )
  throw( std::invalid_argument( "ThermalUnitExtDPSolver::load_parameters: "
                                "secondary reserve not supported." ) );

 // scalar parameters of the unit
 time_horizon      = b->get_time_horizon();
 init_up_down_time = b->get_init_up_down_time();
 min_up_time       = b->get_min_up_time();
 min_down_time     = b->get_min_down_time();
 initial_power     = b->get_initial_power();

 // t_init: first instant at which the commitment decision is genuinely
 // free. If the unit was on before the horizon and min_up_time has not
 // been satisfied yet, the first few instants must stay on; symmetrically
 // for init_up_down_time < 0 and min_down_time. This mirrors the same
 // bookkeeping of ThermalUnitDPSolver and is used only to correctly
 // compute the start_up / shut_down indicator variables in get_var_solution.
 if( init_up_down_time > 0 )
  if( min_up_time > Index( init_up_down_time ) )
   t_init = std::min( time_horizon , min_up_time - Index( init_up_down_time ) );
  else
   t_init = 0;
 else
  if( Index( - init_up_down_time ) < min_down_time )
   t_init = std::min( time_horizon ,
                      min_down_time - Index( - init_up_down_time ) );
  else
   t_init = 0;

 // startup costs, power bounds and shutdown/startup ramp bounds
 startup_costs = b->get_start_up_cost();
 min_power     = b->get_min_power();
 max_power     = b->get_max_power();
 bound_on      = b->get_start_up_limit();
 bound_down    = b->get_shut_down_limit();

 // ramp-up/down limits default to max_power (no effective ramping) when
 // the Block does not set them explicitly
 if( b->get_delta_ramp_up().empty() )
  delta_ramp_up = max_power;
 else
  delta_ramp_up = b->get_delta_ramp_up();

 if( b->get_delta_ramp_down().empty() )
  delta_ramp_down = max_power;
 else
  delta_ramp_down = b->get_delta_ramp_down();

 // cost coefficients f_t(p) = quad_term[t] p^2 + linear_term[t] p + const_term[t]
 retrieve_term( quad_term   , b->get_quad_term() );
 retrieve_term( linear_term , b->get_linear_term() );
 retrieve_term( const_term  , b->get_const_term() );

 if( ! owned )
  f_Block->read_unlock();

 // reset output buffers and pipeline state
 P.assign( time_horizon , 0 );
 U.assign( time_horizon , false );
 stage = start;
 f_solved = false;
 f_best_cost = TUEDPINF;

 }  // end( ThermalUnitExtDPSolver::load_parameters )

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
    min_up_time       = b->get_min_up_time();
    min_down_time     = b->get_min_down_time();
    if( init_up_down_time > 0 )
     if( min_up_time > Index( init_up_down_time ) )
      t_init = std::min( time_horizon , min_up_time - Index( init_up_down_time ) );
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
   }
  return( true );
  }

 return( false );

 }  // end( ThermalUnitExtDPSolver::guts_of_process_modifications )

/*--------------------------------------------------------------------------*/

// Expand a cost-like input vector coming from ThermalUnitBlock into a
// full time_horizon-long vector stored locally. ThermalUnitBlock uses
// three conventions for these vectors:
//   - empty  : the feature is absent / zero for the whole horizon;
//   - size 1 : one scalar value that applies to every time step;
//   - size time_horizon : one value per time step.
// Centralising the expansion here keeps the run_DP() loop free of these
// conditionals: it can always index into the local vector with t and
// get the correct per-instant value.

void ThermalUnitExtDPSolver::retrieve_term(
 std::vector< double > & out , const std::vector< double > & in ) const
{
 if( in.empty() ) {
  out.assign( time_horizon , 0 );
  return;
  }

 if( in.size() == 1 ) {
  out.assign( time_horizon , in[ 0 ] );
  return;
  }

 out = in;
 }

/*--------------------------------------------------------------------------*/
/*------------------- PIECEWISE-QUADRATIC HELPERS --------------------------*/
/*--------------------------------------------------------------------------*/

// Evaluate F at a single point p. Empty F denotes +INF; p outside the
// function's support likewise returns +INF. Linear scan is fine in
// practice since RRF+ keeps the number of pieces very small (the paper
// reports m <= 10 in their benchmarks).

double ThermalUnitExtDPSolver::eval( const PQFun & F , double p )
{
 if( F.empty() )
  return( TUEDPINF );
 for( const auto & pc : F )
  if( pc.left - 1e-12 <= p && p <= pc.right + 1e-12 )
   return( eval_piece( pc , p ) );
 return( TUEDPINF );
 }

/*--------------------------------------------------------------------------*/

// Minimiser of a single quadratic piece on its own [left, right]. Three
// sub-cases:
//  - strictly convex (alfa > 0): unconstrained vertex is -beta/(2 alfa),
//    clamped to [left, right];
//  - linear (alfa == 0, beta != 0): min at the endpoint in the direction
//    of -beta;
//  - constant: any point works; pick left for determinism.

double ThermalUnitExtDPSolver::argmin_piece( const PieceQuad & pc )
{
 if( pc.alfa > 1e-16 ) {
  double p = - pc.beta / ( 2.0 * pc.alfa );
  if( p < pc.left ) return( pc.left );
  if( p > pc.right ) return( pc.right );
  return( p );
  }
 if( pc.beta > 0 ) return( pc.left );
 if( pc.beta < 0 ) return( pc.right );
 return( pc.left );
 }

/*--------------------------------------------------------------------------*/

// min of F over [lo, hi] intersected with dom(F). For each piece, clamp
// the piece's own interval to [lo, hi], run argmin_piece() on the
// resulting sub-piece and keep the best (value, argmin) pair seen.
// Empty F or empty [lo, hi] yields (+INF, 0). The result is used both
// to summarise F^tau_t in the OnSlot associated to each surviving
// function and to compute v_shutdown[h] (min over tau >= M_up of F^tau_h
// restricted to the shutdown range [P, SD]).

std::pair< double , double > ThermalUnitExtDPSolver::min_over(
 const PQFun & F , double lo , double hi )
{
 if( F.empty() || lo > hi + 1e-12 )
  return std::make_pair( TUEDPINF , 0.0 );

 double best_v = TUEDPINF;
 double best_p = 0;

 for( const auto & pc : F ) {
  double l = std::max( pc.left  , lo );
  double r = std::min( pc.right , hi );
  if( l > r + 1e-12 )
   continue;
  // clamp pc's own interval to [lo, hi] before taking its argmin
  PieceQuad cpc = { pc.alfa , pc.beta , pc.gamma , l , r };
  double p = argmin_piece( cpc );
  double v = eval_piece( cpc , p );
  if( v < best_v ) {
   best_v = v;
   best_p = p;
   }
  }

 return std::make_pair( best_v , best_p );
 }

/*--------------------------------------------------------------------------*/

// Pointwise add the scalar c to every piece (equivalent to shifting the
// function vertically). Only the constant coefficient gamma is touched;
// the piece intervals and higher-order coefficients are preserved.

void ThermalUnitExtDPSolver::shift_by( PQFun & F , double c )
{
 for( auto & pc : F )
  pc.gamma += c;
 }

/*--------------------------------------------------------------------------*/

// Pointwise add the quadratic (alfa p^2 + beta p + gamma) to every
// piece. This is the "+ f_t(p)" step of the DP recurrence: after
// sliding_min() gives a function representing the cheapest path to
// reach the power p at time t from some predecessor, we add the cost
// of *producing* p at time t. Convexity is preserved because alfa >= 0
// is added to each alfa-coefficient, which is already >= 0.

void ThermalUnitExtDPSolver::add_quadratic( PQFun & F ,
                                            double alfa , double beta ,
                                            double gamma )
{
 for( auto & pc : F ) {
  pc.alfa  += alfa;
  pc.beta  += beta;
  pc.gamma += gamma;
  }
 }

/*--------------------------------------------------------------------------*/

// Restrict F's domain to [lo, hi], dropping pieces that fall entirely
// outside and truncating pieces that only partially intersect. The
// piece coefficients themselves are preserved; only the interval
// endpoints are updated.

void ThermalUnitExtDPSolver::clamp_domain( PQFun & F , double lo , double hi )
{
 PQFun out;
 out.reserve( F.size() );
 for( const auto & pc : F ) {
  double l = std::max( pc.left  , lo );
  double r = std::min( pc.right , hi );
  if( l < r - 1e-15 )
   out.push_back( { pc.alfa , pc.beta , pc.gamma , l , r } );
  }
 F = std::move( out );
 }

/*--------------------------------------------------------------------------*/

// Return true iff F1 is everywhere dominated by F2 on its own domain,
// i.e., F1(p) >= F2(p) - eps for every p in dom(F1). Used by RRF+ to
// prune irrelevant F^tau functions (Wuijts et al. 2021, Prop. 6.1): a
// dominated function cannot contribute to the optimal solution since
// its downstream F^{tau+k}_{t+k} (obtained by sliding_min + f_t + ...)
// is dominated step-by-step by the analogous propagation of the
// dominating function.
//
// The test is exact: it walks the two functions piece-by-piece with a
// two-pointer sweep and, over each common sub-interval, minimises the
// pointwise difference
//
//    D(p) = (alfa_1 - alfa_2) p^2 + (beta_1 - beta_2) p + (gamma_1 - gamma_2)
//
// analytically. If at some sub-interval min D < -eps the domination
// fails; if F1 reaches a point where F2 is not defined, F2 is +INF
// there and F1 is not dominated (early-return false). Coverage gaps in
// F2's pieces are likewise failures.
//
// The empty function is the +INF function, so:
//  - F1 empty  => F1 = +INF <= +INF: vacuously dominated (return true);
//  - F2 empty  => F2 = +INF: nothing can be dominated by it (return false).

bool ThermalUnitExtDPSolver::is_dominated_by( const PQFun & F1 ,
                                              const PQFun & F2 , double eps )
{
 if( F1.empty() )
  return( true );
 if( F2.empty() )
  return( false );

 const double tol = 1e-12;

 // outer loop: iterate over the pieces of F1 in order
 for( std::size_t i = 0 ; i < F1.size() ; ++i ) {
  double l1 = F1[ i ].left;
  double r1 = F1[ i ].right;
  double p  = l1;
  std::size_t j = 0;
  // advance j to the first F2-piece that overlaps [p, r1] at all
  while( j < F2.size() && F2[ j ].right <= p + tol )
   ++j;
  // consume [p, r1] by stepping through the F2-pieces it touches
  while( p < r1 - tol ) {
   if( j >= F2.size() )
    return( false );     // F2 does not cover [p, r1] on the right
   double l2 = F2[ j ].left;
   double r2 = F2[ j ].right;
   if( l2 > p + tol )
    return( false );     // gap in F2 coverage inside F1's support
   double seg_l = p;
   double seg_r = std::min( r1 , r2 );
   // compute coefficients of the difference D on the common segment
   double a = F1[ i ].alfa  - F2[ j ].alfa;
   double b = F1[ i ].beta  - F2[ j ].beta;
   double c = F1[ i ].gamma - F2[ j ].gamma;
   // analytical minimum of a p^2 + b p + c on [seg_l, seg_r]
   double min_D;
   if( a > tol ) {
    // strictly convex: vertex at -b / (2a), clamped to the segment
    double p_opt = - b / ( 2.0 * a );
    if( p_opt < seg_l ) p_opt = seg_l;
    if( p_opt > seg_r ) p_opt = seg_r;
    min_D = a * p_opt * p_opt + b * p_opt + c;
    }
   else {
    // concave (a < 0) or linear/constant (a ~= 0): the minimum is at
    // one of the endpoints; evaluate both and pick the smaller
    double vl = a * seg_l * seg_l + b * seg_l + c;
    double vr = a * seg_r * seg_r + b * seg_r + c;
    min_D = std::min( vl , vr );
    }
   if( min_D < - eps )
    return( false );     // F1 strictly below F2 somewhere in the segment
   p = seg_r;
   // advance to the next F2-piece if the current one ends before r1
   if( r2 < r1 - tol )
    ++j;
   else
    break;
   }
  }
 return( true );
 }

/*--------------------------------------------------------------------------*/

// Ramp-constrained sliding minimum of a convex piecewise quadratic F.
// Returns the function
//
//    G(p_t) = min_{ q in [p_t - ramp_up, p_t + ramp_down] } F(q)
//
// restricted to the domain [lo, hi] of p_t. Here ramp_up is Delta+ (the
// maximum per-step increase in power from t-1 to t) and ramp_down is
// Delta- (the maximum per-step decrease): given a target p_t at time t,
// the admissible predecessor power q at time t-1 lies in
// [p_t - ramp_up, p_t + ramp_down]. This is the core one-step transition
// of the DP, i.e., Wuijts et al. (2021) eq. (16), which is in turn the
// three-case analysis of Frangioni and Gentile (2006).
//
// The closed-form construction: let (p_star, f_star) be the unconstrained
// minimiser/value of F. Then
//
//    G(p_t) = | F(p_t + ramp_down) , p_t <  p_star - ramp_down
//             | f_star              , p_t in [p_star - ramp_down,
//             |                             p_star + ramp_up   ]
//             | F(p_t - ramp_up)   , p_t >  p_star + ramp_up
//
// so the pieces of F split into three groups when producing G's pieces:
//
//  - pieces whose domain lies strictly left of p_star contribute to G by
//    shifting -ramp_down in the p_t axis (with a corresponding coefficient
//    substitution q = p_t + ramp_down);
//  - a single *flat* piece of value f_star covers [p_star - ramp_down,
//    p_star + ramp_up];
//  - pieces strictly right of p_star shift +ramp_up (q = p_t - ramp_up).
//
// The piece containing p_star is split between the left-shifted and the
// right-shifted groups, touching the flat middle on both sides. The
// resulting pieces are then sorted, clamped to [lo, hi] and the degenerate
// ones (empty intervals) are dropped. The output is still a convex
// piecewise quadratic on [lo, hi].

ThermalUnitExtDPSolver::PQFun ThermalUnitExtDPSolver::sliding_min(
 const PQFun & F , double ramp_up , double ramp_down ,
 double lo , double hi )
{
 PQFun G;
 if( F.empty() || lo > hi + 1e-12 )
  return( G );

 // clamp negative ramps defensively (should not occur in normal use)
 if( ramp_up   < 0 ) ramp_up   = 0;
 if( ramp_down < 0 ) ramp_down = 0;

 // locate the unconstrained minimiser p_star and its value f_star by
 // scanning each piece's own argmin (all pieces are convex so the piece-wise
 // minimum on each interval is reached at its local vertex or endpoint)
 double p_star = 0;
 double f_star = TUEDPINF;
 for( const auto & pc : F ) {
  double p = argmin_piece( pc );
  double v = eval_piece( pc , p );
  if( v < f_star ) {
   f_star = v;
   p_star = p;
   }
  }

 const double tol = 1e-12;

 // build the output in an intermediate buffer, sort + clamp at the end
 PQFun raw;
 raw.reserve( F.size() * 2 + 1 );

 // (1) left-shifted pieces: portions of F with q <= p_star.
 // The substitution q = p_t + ramp_down implies that the interval
 // [pc.left, pc.right] of the original variable q maps to
 // [pc.left - ramp_down, pc.right - ramp_down] of the new variable p_t,
 // and the coefficients transform as
 //   F(q) = alfa q^2 + beta q + gamma
 //        = alfa (p_t + rd)^2 + beta (p_t + rd) + gamma
 //        = alfa p_t^2 + (2 alfa rd + beta) p_t + (alfa rd^2 + beta rd + gamma)
 for( const auto & pc : F ) {
  double a = pc.left, b = pc.right;
  if( b <= p_star + tol ) {
   // piece entirely to the left of (or touching) p_star: shift whole piece
   double na = pc.alfa;
   double nb = 2 * pc.alfa * ramp_down + pc.beta;
   double nc = pc.alfa * ramp_down * ramp_down +
               pc.beta * ramp_down + pc.gamma;
   raw.push_back( { na , nb , nc , a - ramp_down , b - ramp_down } );
   }
  else if( a < p_star - tol ) {
   // piece straddles p_star: keep only the left portion [a, p_star]
   double na = pc.alfa;
   double nb = 2 * pc.alfa * ramp_down + pc.beta;
   double nc = pc.alfa * ramp_down * ramp_down +
               pc.beta * ramp_down + pc.gamma;
   raw.push_back( { na , nb , nc , a - ramp_down , p_star - ramp_down } );
   }
  }

 // (2) flat middle piece at value f_star covers
 // [p_star - ramp_down, p_star + ramp_up]. Skipped if both ramps are zero
 // (no-op sliding window): the left/right branches already include p_star
 // itself and no middle region remains.
 if( ramp_up > tol || ramp_down > tol )
  raw.push_back( { 0 , 0 , f_star , p_star - ramp_down , p_star + ramp_up } );

 // (3) right-shifted pieces: portions of F with q >= p_star.
 // Symmetric to (1), with substitution q = p_t - ramp_up:
 //   F(q) = alfa (p_t - ru)^2 + beta (p_t - ru) + gamma
 //        = alfa p_t^2 + (-2 alfa ru + beta) p_t + (alfa ru^2 - beta ru + gamma)
 for( const auto & pc : F ) {
  double a = pc.left, b = pc.right;
  if( a >= p_star - tol ) {
   // piece entirely to the right of (or touching) p_star
   double na = pc.alfa;
   double nb = - 2 * pc.alfa * ramp_up + pc.beta;
   double nc = pc.alfa * ramp_up * ramp_up -
               pc.beta * ramp_up + pc.gamma;
   raw.push_back( { na , nb , nc , a + ramp_up , b + ramp_up } );
   }
  else if( b > p_star + tol ) {
   // straddling piece: keep only the right portion [p_star, b]
   double na = pc.alfa;
   double nb = - 2 * pc.alfa * ramp_up + pc.beta;
   double nc = pc.alfa * ramp_up * ramp_up -
               pc.beta * ramp_up + pc.gamma;
   raw.push_back( { na , nb , nc , p_star + ramp_up , b + ramp_up } );
   }
  }

 // sort the pieces by left endpoint so that the output respects the
 // PQFun invariant (pieces in increasing order of p)
 std::sort( raw.begin() , raw.end() ,
            []( const PieceQuad & a , const PieceQuad & b ) {
             return( a.left < b.left );
             } );

 // clamp to the output domain [lo, hi] and drop degenerate pieces
 G.reserve( raw.size() );
 for( const auto & pc : raw ) {
  double l = std::max( pc.left  , lo );
  double r = std::min( pc.right , hi );
  if( l < r - 1e-15 )
   G.push_back( { pc.alfa , pc.beta , pc.gamma , l , r } );
  }

 return( G );

 }  // end( ThermalUnitExtDPSolver::sliding_min )

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

 // initialise all the per-time-step state vectors to empty / +INF
 f_F  .assign( time_horizon , {} );
 f_tau.assign( time_horizon , {} );
 f_on .assign( time_horizon , {} );
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
  std::vector< char > keep( sz , 1 );
  for( std::size_t i = 0 ; i < sz ; ++i ) {
   if( ! keep[ i ] ) continue;
   if( v_tau[ i ] < mut ) continue;  // keep all tau < mut
   for( std::size_t j = 0 ; j < sz ; ++j ) {
    if( j == i ) continue;
    if( ! keep[ j ] ) continue;
    if( v_tau[ j ] < mut ) continue;
    if( is_dominated_by( v_F[ i ] , v_F[ j ] ) ) {
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

 // helper to compute v_shutdown[t] from the current f_F[t]
 auto compute_v_shutdown = [ & ]( Index t , double sd_hi ) {
  double best = TUEDPINF;
  Index best_tau = 0;
  double best_p = 0;
  double Plo = min_power[ t ];
  for( std::size_t i = 0 ; i < f_F[ t ].size() ; ++i ) {
   if( f_tau[ t ][ i ] < mut ) continue;
   auto [ v , p ] = min_over( f_F[ t ][ i ] , Plo , sd_hi );
   if( v < best ) { best = v; best_tau = f_tau[ t ][ i ]; best_p = p; }
   }
  v_shutdown[ t ]     = best;
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
  if( lo < hi + 1e-12 ) {
   PQFun F;
   F.push_back( { quad_term[ 0 ] , linear_term[ 0 ] , const_term[ 0 ] ,
                  lo , hi } );
   auto [ v , p ] = min_over( F , lo , hi );
   f_F  [ 0 ].push_back( std::move( F ) );
   f_tau[ 0 ].push_back( tau0 );
   f_on [ 0 ].push_back( { v , p } );
   }
  // OFF-side scalars stay +INF: the unit is on, no off path has been
  // produced yet
  }
 else {
  // init_up_down_time <= 0: the unit has been off for |init| instants
  // strictly before t = 0. We can immediately *restart* at t = 0 if the
  // off period satisfies mdt, i.e. iff |init| >= mdt.
  bool can_restart_t0 = ( Index( - init_up_down_time ) >= mdt );
  if( can_restart_t0 ) {
   // F^1_0(p) = f_0(p) + SUC[0] on [P, min(Pbar, SU)]: the restart arc
   // pays the start-up cost and constrains p by the start-up ramp limit
   double lo = min_power[ 0 ];
   double hi = std::min( max_power[ 0 ] , bound_on[ 0 ] );
   if( lo < hi + 1e-12 ) {
    double suc = startup_costs.empty() ? 0.0 : startup_costs[ 0 ];
    PQFun F;
    F.push_back( { quad_term[ 0 ] , linear_term[ 0 ] ,
                   const_term[ 0 ] + suc , lo , hi } );
    auto [ v , p ] = min_over( F , lo , hi );
    f_F  [ 0 ].push_back( std::move( F ) );
    f_tau[ 0 ].push_back( 1 );
    f_on [ 0 ].push_back( { v , p } );
    }
   }
  // "unit off at t = 0, any duration": cost 0 from the initial off state;
  // f_any_pred[0] = -1 signals that no in-horizon shutdown happened
  c_off_any[ 0 ]  = 0;
  f_any_pred[ 0 ] = -1;
  // "unit off at t = 0 AND ready": the off trail before t = 0 counts
  // |init| instants and t = 0 itself adds one more, so the condition is
  // |init| + 1 >= mdt (equivalently, the unit is ready right at t = 0)
  if( Index( - init_up_down_time ) + 1 >= mdt ) {
   c_off_ready[ 0 ]  = 0;
   f_ready_pred[ 0 ] = -1;
   }
  }

 // prune at t = 0 (with a single seeded entry this is a no-op, but for
 // completeness we apply it to every time step)
 prune_RRF_plus( f_F[ 0 ] , f_tau[ 0 ] , f_on[ 0 ] );

 // v_shutdown[0]: only meaningful if there is at least one more instant
 // in the horizon after t = 0 (so that being off at t = 1 is possible)
 if( n > 1 )
  compute_v_shutdown( 0 , bound_down[ 1 ] );

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
  double stay_ready  = c_off_ready[ t - 1 ];
  double fresh_ready = TUEDPINF;
  int fresh_h = -1;
  if( t >= mdt ) {
   fresh_ready = v_shutdown[ t - mdt ];
   fresh_h = int( t - mdt );
   }
  double init_ready = TUEDPINF;
  if( ( init_up_down_time <= 0 ) &&
      ( Index( - init_up_down_time ) + t + 1 >= mdt ) )
   init_ready = 0;
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

  // -- ON-side updates ------------------------------------------------- //
  // build a fresh sparse list for time t from (a) a possible restart
  // tau=1 and (b) the increments of each surviving entry at t-1
  double Plo     = min_power[ t ];
  double Phi     = max_power[ t ];
  double ru_prev = delta_ramp_up  [ t - 1 ];
  double rd_prev = delta_ramp_down[ t - 1 ];

  std::vector< PQFun > new_F;
  std::vector< Index > new_tau;
  std::vector< OnSlot > new_on;
  new_F  .reserve( f_F[ t - 1 ].size() + 1 );
  new_tau.reserve( f_F[ t - 1 ].size() + 1 );
  new_on .reserve( f_F[ t - 1 ].size() + 1 );

  // (a) tau = 1 entry (restart arc): F^1_t(p) = f_t(p) + SUC[t]
  //     + c_off_ready[t-1], defined for p in [P, min(Pbar, SU)]. Only
  //     built if c_off_ready[t-1] is finite (otherwise the unit cannot
  //     legally restart at t).
  if( c_off_ready[ t - 1 ] < TUEDPINF ) {
   double lo = Plo;
   double hi = std::min( Phi , bound_on[ t ] );
   if( lo < hi + 1e-12 ) {
    double suc = startup_costs.empty() ? 0.0 : startup_costs[ t ];
    PQFun F;
    F.push_back( { quad_term[ t ] , linear_term[ t ] ,
                   const_term[ t ] + suc + c_off_ready[ t - 1 ] ,
                   lo , hi } );
    auto [ v , p ] = min_over( F , lo , hi );
    new_F  .push_back( std::move( F ) );
    new_tau.push_back( 1 );
    new_on .push_back( { v , p } );
    }
   }

  // (b) tau > 1 entries: from each surviving entry of f_F[t-1] apply
  //     sliding_min (ramp window from t-1 to t) to get the value
  //     function restricted to [Plo, Phi], then add f_t(p). The output
  //     tau is the previous tau + 1 (continuing the on-run).
  for( std::size_t i = 0 ; i < f_F[ t - 1 ].size() ; ++i ) {
   PQFun F = sliding_min( f_F[ t - 1 ][ i ] , ru_prev , rd_prev , Plo , Phi );
   if( F.empty() ) continue;   // intersection with [Plo, Phi] empty
   add_quadratic( F , quad_term[ t ] , linear_term[ t ] , const_term[ t ] );
   auto [ v , p ] = min_over( F , Plo , Phi );
   new_F  .push_back( std::move( F ) );
   new_tau.push_back( f_tau[ t - 1 ][ i ] + 1 );
   new_on .push_back( { v , p } );
   }

  // the new_tau vector is already sorted ascending: the restart entry
  // (tau = 1) was pushed first, and each subsequent entry has
  // tau = (old tau) + 1 >= 2 with the old taus appearing in ascending
  // order; since the old list had tau >= 1, we have new tau >= 2 and
  // monotonicity is preserved.

  // RRF+ pruning: remove functions with tau >= mut that are pointwise
  // dominated by some other function with tau >= mut in the new list.
  prune_RRF_plus( new_F , new_tau , new_on );

  f_F  [ t ] = std::move( new_F );
  f_tau[ t ] = std::move( new_tau );
  f_on [ t ] = std::move( new_on );

  // v_shutdown[t] is meaningful only if there is at least one more
  // instant in the horizon after t, since an in-horizon off period
  // triggered by shutdown at end of t starts at t + 1
  if( t < n - 1 )
   compute_v_shutdown( t , bound_down[ t + 1 ] );

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
   best_on  = f_on[ n - 1 ][ i ].min_val;
   best_tau = f_tau[ n - 1 ][ i ];
   best_p   = f_on[ n - 1 ][ i ].argmin_p;
   }
 double best_off = c_off_any[ n - 1 ];

 // walk_on: steps one on-run back from (t, tau, p), writing P[] / U[] as
 // it goes, and returns the time at which the run started. The caller
 // then decides whether (and where) to resume with an earlier on-run via
 // f_ready_pred[].
 auto walk_on = [ & ]( Index t , Index tau , double p ) -> Index {
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
   // compute the predecessor power via Wuijts eq. (16) three-case formula
   double p_prev_star = f_on[ t - 1 ][ idx ].argmin_p;
   double ru = delta_ramp_up  [ t - 1 ];
   double rd = delta_ramp_down[ t - 1 ];
   double p_prev;
   if( p > p_prev_star + ru )
    p_prev = p - ru;        // right branch: closest window endpoint to p*
   else
    if( p < p_prev_star - rd )
     p_prev = p + rd;       // left branch: symmetric
    else
     p_prev = p_prev_star;  // middle branch: p* itself is in the window
   // clamp p_prev to the predecessor's own domain; when tau-1 == 1 the
   // domain is additionally constrained by the start-up ramp SU
   double lo_prev = min_power[ t - 1 ];
   double hi_prev = max_power[ t - 1 ];
   if( tau - 1 == 1 )
    hi_prev = std::min( hi_prev , bound_on[ t - 1 ] );
   if( p_prev < lo_prev ) p_prev = lo_prev;
   if( p_prev > hi_prev ) p_prev = hi_prev;
   t   -= 1;
   tau -= 1;
   p    = p_prev;
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
 bool have_on;

 // bootstrap the alternating on/off walk
 if( best_on <= best_off ) {
  // on-termination: the last on-run ends at (n-1, best_tau, best_p)
  on_t   = n - 1;
  on_tau = best_tau;
  on_p   = best_p;
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
   on_t   = Index( h );
   on_tau = v_shutdown_tau[ h ];
   on_p   = v_shutdown_p[ h ];
   have_on = true;
   }
  }

 // alternate on-runs and off-gaps until we reach the beginning of the
 // horizon (t_start == 0) or an off-gap that was inherited from before
 // the horizon (f_ready_pred[t_start - 1] == -1)
 while( have_on ) {
  Index t_start = walk_on( on_t , on_tau , on_p );

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
  // jump to the shutdown event and continue the walk
  on_t   = Index( h );
  on_tau = v_shutdown_tau[ h ];
  on_p   = v_shutdown_p[ h ];
  }

 stage = sol_OK;

 }  // end( ThermalUnitExtDPSolver::build_solution )

/*--------------------------------------------------------------------------*/
/*-------------- End File ThermalUnitExtDPSolver.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
