/*--------------------------------------------------------------------------*/
/*----------------------- File ThermalUnitDPSolver.cpp ---------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the ThermalUnitDPSolver class.
 *
 * \author Claudio Gentile \n
 *         Istituto di Analisi di Sistemi e Informatica "Antonio Ruberti" \n
 *         Consiglio Nazionale delle Ricerche \n
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Niccolo' Iardella \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Claudio Gentile, Antonio Frangioni, Niccolo' Iardella,
 *                      Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------- MACROS -----------------------------------*/
/*--------------------------------------------------------------------------*/

#define COMPUTE_DUALS 0
/* If COMPUTE_DUALS > 0, the solver allocates more memory and store more
 * information about the solution process in such a way as to make it possible
 * to reconstruct the optimal dual solution in the end. However, this is not
 * implemented yet, so that currently the setting makes no sense. */

#define TUDPS_PROFILE 0
/* If TUDPS_PROFILE > 0, compute() accumulates per-phase wall-clock time
 * (build_graph / compute_EDPs / min_path / compute_solutions) into static
 * counters and prints a cumulative summary to cerr at each call. Only the DP
 * phases are timed, so any other solver running in the same process (e.g. the
 * MILP comparison solver in the test harness) does not pollute the figures.
 * Temporary instrumentation: keep at 0 in committed code. */

/* TUDPS_PARALLEL and the default TUDPS_PAR_MIN_N are defined in the header.
 * Below the time horizon f_par_min_n (the run-time intParMinN parameter,
 * defaulting to TUDPS_PAR_MIN_N) compute_EDPs() stays serial, since the
 * thread-dispatch overhead is not amortised on short horizons. */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "ThermalUnitDPSolver.h"

#include "ThermalUnitBlock.h"

#if TUDPS_PROFILE
 #include <chrono>
 #include <iostream>
#endif

#if TUDPS_PARALLEL
 #include <thread>
 #include <ff/parallel_for.hpp>
 // Keep this include: ff/pipeline.hpp forward-declares the isa2a_get*set()
 // helpers as static and uses them in ff_pipeline members, but their definitions
 // live in ff/graph_utils.hpp, which only ff/ff.hpp pulls in. GCC/Clang tolerate
 // the undefined static (the calling members are never emitted here), but MSVC
 // rejects it with C2129. Including the definitions here is harmless everywhere
 // (the header is include-guarded and self-contained), so we keep it
 // unconditional rather than guarding it on the compiler/platform.
 #include <ff/graph_utils.hpp>
#endif

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register ThermalUnitDPSolver to the Block factory

SMSpp_insert_in_factory_cpp_0( ThermalUnitDPSolver );

/*--------------------------------------------------------------------------*/
/*--------------------------- Solver INTERFACE -----------------------------*/
/*--------------------------------------------------------------------------*/

// defined here (not defaulted in the header) so that the destructor of the
// pimpl'd ff::ParallelFor is instantiated where the type is complete

ThermalUnitDPSolver::~ThermalUnitDPSolver() = default;

/*--------------------------------------------------------------------------*/

void ThermalUnitDPSolver::set_Block( Block * block )
{
 if( block == f_Block )
  return;

 Solver::set_Block( block );

 if( block ) {
  // note that we do not use
  //
  //    if( ! dynamic_cast< ThermalUnitBlock * >( f_Block ) )
  //
  // since we want to avoid that the check succeeds for possible derived
  // classes of ThermalUnitBlock which may have other features / constraints
  // that ThermalUnitDPSolver does not know about and therefore cannot
  // properly handle. typeid() works in this case since ThermalUnitBlock
  // is a  polymorphic object, i.e., it has at least one virtual method
  if( typeid( ThermalUnitBlock ) != typeid( *f_Block ) )
   throw( std::runtime_error(
    "ThermalUnitDPSolver::set_Block: ThermalUnitBlock required." ) );
  load_parameters();
  }
 }

/*--------------------------------------------------------------------------*/

int ThermalUnitDPSolver::compute( bool changedvars )
{
 lock();  // lock the mutex

 process_modifications();

#if TUDPS_PROFILE
 static double t_bg = 0 , t_ed = 0 , t_mp = 0 , t_cs = 0;
 static unsigned long n_bg = 0 , n_ed = 0 , n_mp = 0 , n_cs = 0 , n_call = 0;
 using clk = std::chrono::steady_clock;
 auto tic = clk::now();
 #define TUDPS_TIME( acc , cnt ) { auto now = clk::now(); \
   acc += std::chrono::duration< double >( now - tic ).count(); ++cnt; \
   tic = now; }
 switch( stage ) {
  case( start ):    build_graph();       TUDPS_TIME( t_bg , n_bg )
  case( graph_OK ): compute_EDPs();      TUDPS_TIME( t_ed , n_ed )
  case( edps_OK ):  min_path();          TUDPS_TIME( t_mp , n_mp )
  case( path_OK ):  compute_solutions(); TUDPS_TIME( t_cs , n_cs )
  }
 #undef TUDPS_TIME
 ++n_call;
 std::cerr << "TUDPS_PROF n=" << time_horizon << " calls=" << n_call
           << " build_graph=" << t_bg << "s/" << n_bg
           << " compute_EDPs=" << t_ed << "s/" << n_ed
           << " min_path=" << t_mp << "s/" << n_mp
           << " compute_solutions=" << t_cs << "s/" << n_cs << std::endl;
#else
 switch( stage ) {
  case( start ):    build_graph();
  case( graph_OK ): compute_EDPs();
  case( edps_OK ):  min_path();
  case( path_OK ):  compute_solutions();
  }
#endif

 unlock();  // unlock the mutex

 assert( stage == sol_OK );
 return( f_end.lab == TUDPINF ? kInfeasible : kOK );
 }

/*--------------------------------------------------------------------------*/

void ThermalUnitDPSolver::get_var_solution( Configuration * solc )
{
 // lock the Block
 bool owned = f_Block->is_owned_by( f_id );
 if( ( ! owned ) && ( ! f_Block->lock( f_id ) ) )
  throw( std::runtime_error(
   "ThermalUnitDPSolver::get_var_solution: unable to lock the Block." ) );

 auto b = static_cast< ThermalUnitBlock * >( f_Block );

 // design (investment): write the design variable and, when the unit is not
 // built, force the whole schedule to zero overriding the DP's (P, U) which
 // were computed assuming the unit exists. "built" is true when there is no
 // design at all.
 const bool built = ( ! has_design ) || design_on;
 if( has_design )
  b->get_design().set_value( design_on ? 1 : 0 );

 // canonical part: write the schedule the DP produced, (P, U), into the
 // active power and commitment ColVariables of the Block
 if( auto pow_it = b->get_active_power( 0 ) )
  for( Index i = 0 ; i < time_horizon ; ++i )
   ( pow_it++ )->set_value( built ? P[ i ] : 0 );

 if( auto com_it = b->get_commitment( 0 ) )
  for( Index i = 0 ; i < time_horizon ; ++i )
   ( com_it++ )->set_value( ( built && U[ i ] ) ? 1 : 0 );

 // spinning reserve variables (if present): the optimal pr/sr provision given
 // the active power P[i]. The reserve band must use the SAME cap as the value
 // computation, or the recovered solution would disagree with the reported
 // optimal value: bound_on at a start-up period, bound_down[i+1] at a shut-down
 // period (off at i+1) -- the boundary correction -- and max_power at an
 // interior period (or a unit on to the horizon end, a tail with no shut-down).
 auto band_cap = [ & ]( Index i ) -> double {
  if( ! ( built && U[ i ] ) )
   return( max_power[ i ] );                        // off: reserve is 0 anyway
  const bool is_su = ( i == 0 ) ? ( init_up_down_time <= 0 )
                                : ( ! U[ i - 1 ] );
  if( is_su )
   return( bound_on[ i ] );                          // start-up
  if( ( i + 1 < time_horizon ) && ( ! U[ i + 1 ] ) )
   return( bound_down[ i + 1 ] );                    // shut-down (off at i+1)
  return( max_power[ i ] );                           // interior / on-to-end
  };

 if( auto pr_it = b->get_primary_spinning_reserve( 0 ) )
  for( Index i = 0 ; i < time_horizon ; ++i ) {
   double pr , sr;
   reserve_alloc( i , built ? P[ i ] : 0 , pr , sr , band_cap( i ) );
   ( pr_it++ )->set_value( pr );
   }

 if( auto sr_it = b->get_secondary_spinning_reserve( 0 ) )
  for( Index i = 0 ; i < time_horizon ; ++i ) {
   double pr , sr;
   reserve_alloc( i , built ? P[ i ] : 0 , pr , sr , band_cap( i ) );
   ( sr_it++ )->set_value( sr );
   }

 // reactive power variables (AC instances): q[t] in [ Qmin(t) , Qmax(t) ] is
 // separable from the DP; under a dualizing Solver it carries the linear cost
 // reactive_linear_term[t] and the optimal q*[t] is the box endpoint
 // minimising c*q (lower if c>0, upper if c<0, else feasible 0). This matches
 // the contribution run_DP() adds to the value. See the analogous comment in
 // ThermalUnitExtDPSolver::get_var_solution().
 if( auto q_it = b->get_reactive_power( 0 ) )
  for( Index i = 0 ; i < time_horizon ; ++i ) {
   double q = 0;  // not built: q is forced to 0 with all operational variables
   if( built ) {
    const double qlo = b->get_min_reactive_power( i );
    const double qhi = b->get_max_reactive_power( i );
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

 // formulation-specific bookkeeping (start_up / shut_down indicators,
 // perspective-cut auxiliaries, ...) is delegated to the Block, which
 // knows which variables exist in the current formulation and how they
 // relate to (p, u). This way the DP stays formulation-agnostic.
 b->set_solution();

 // unlock the Block
 if( ! owned )
  f_Block->unlock( f_id );

 }  // end( ThermalUnitDPSolver::get_var_solution )

/*--------------------------------------------------------------------------*/

// Per-period spinning-reserve LP: given the unit on at t with power p, find
// the optimal primary (pr) and secondary (sr) reserve provision. Each unit
// of reserve costs its objective coefficient (cp / cs); since reserve only
// appears in upper-bound constraints, it is worth providing only when that
// coefficient is negative (a Lagrangian reward). Then we greedily fill the
// more rewarding (more negative) reserve first, up to its participation cap
// (rho*p) and the shared headroom. Returns the optimal per-period reserve
// cost g_t(p) = cp*pr + cs*sr (0 in the standalone case, where cp = rho_p >= 0
// and cs = rho_s >= 0).

double ThermalUnitDPSolver::reserve_alloc( Index t , double p ,
                                           double & pr , double & sr ,
                                           double cap ) const
{
 pr = sr = 0;

 // symmetric reserve band beta_t(p) = min( p - min_power , cap - p ): the
 // reserve must fit both above (the cap constraint, p + pr + sr <= cap) and
 // below (min_power constraint, p - pr - sr >= min_power). cap is max_power at
 // an interior period and the tighter start-up/shut-down cap at a boundary one.
 const double H = std::min( p - min_power[ t ] , cap - p );
 if( H <= 0 )
  return( 0 );

 const double cap_p = primary_rho.empty()   ? 0 : primary_rho[ t ]   * p;
 const double cap_s = secondary_rho.empty() ? 0 : secondary_rho[ t ] * p;
 const double cp = primary_reserve_cost.empty()   ? 0
                                                  : primary_reserve_cost[ t ];
 const double cs = secondary_reserve_cost.empty() ? 0
                                                  : secondary_reserve_cost[ t ];

 double rem = H;
 if( ( cp < 0 ) && ( ( cs >= 0 ) || ( cp <= cs ) ) ) {
  pr = std::min( cap_p , rem );  rem -= pr;
  if( cs < 0 )
   sr = std::min( cap_s , rem );
  }
 else if( cs < 0 ) {
  sr = std::min( cap_s , rem );  rem -= sr;
  if( cp < 0 )
   pr = std::min( cap_p , rem );
  }

 return( cp * pr + cs * sr );

 }  // end( ThermalUnitDPSolver::reserve_alloc )

/*--------------------------------------------------------------------------*/

bool ThermalUnitDPSolver::reserve_rewarded( void ) const
{
 for( auto c : primary_reserve_cost )
  if( c < 0 )
   return( true );
 for( auto c : secondary_reserve_cost )
  if( c < 0 )
   return( true );
 return( false );

 }  // end( ThermalUnitDPSolver::reserve_rewarded )

/*--------------------------------------------------------------------------*/

// Build g_t(p) (the reserve "discount") as convex piecewise-linear pieces on
// [ min_power[t] , max_power[t] ]. g_t is nonpositive when some reserve price
// is negative; it is built exactly (no approximation) by evaluating
// reserve_alloc() at the analytic breakpoints where the binding constraint of
// the per-period reserve LP changes (band vs participation caps). Mirrors
// ThermalUnitExtDPSolver::build_reserve_discount.

std::vector< ThermalUnitDPSolver::g_piece >
ThermalUnitDPSolver::build_reserve_discount( Index t , double cap ) const
{
 std::vector< g_piece > G;
 const double cp = primary_reserve_cost.empty()   ? 0
                                                  : primary_reserve_cost[ t ];
 const double cs = secondary_reserve_cost.empty() ? 0
                                                  : secondary_reserve_cost[ t ];
 if( ( cp >= 0 ) && ( cs >= 0 ) )
  return( G );  // no negative price: g_t == 0

 const double lo = min_power[ t ];
 const double hi = cap;  // upper power cap U_t (interior max_power, or the
                         // start-up/shut-down cap at a boundary period)
 if( hi <= lo + 1e-12 )
  return( G );

 // participation caps of the active (negative-price) reserves, ordered
 // most-negative-price first to match reserve_alloc's greedy fill
 const double rp = primary_rho.empty()   ? 0 : primary_rho[ t ];
 const double rs = secondary_rho.empty() ? 0 : secondary_rho[ t ];
 double rho1 = 0 , rho2 = 0;
 if( ( cp < 0 ) && ( cs < 0 ) ) {
  if( cp <= cs ) { rho1 = rp; rho2 = rs; }
  else           { rho1 = rs; rho2 = rp; }
  }
 else if( cp < 0 )
  rho1 = rp;
 else
  rho1 = rs;  // cs < 0

 std::vector< double > bp = { lo , hi , 0.5 * ( lo + hi ) };
 auto add_bp = [ & ]( double d ) {
  if( ( d > lo + 1e-12 ) && ( d < hi - 1e-12 ) ) bp.push_back( d );
  };
 if( rho1 < 1 ) add_bp( lo / ( 1 - rho1 ) );
 add_bp( hi / ( 1 + rho1 ) );
 const double rsum = rho1 + rho2;
 if( rho2 > 0 ) {
  if( rsum < 1 ) add_bp( lo / ( 1 - rsum ) );
  add_bp( hi / ( 1 + rsum ) );
  }
 std::sort( bp.begin() , bp.end() );
 bp.erase( std::unique( bp.begin() , bp.end() ,
            []( double a , double b ) { return( b - a <= 1e-12 ); } ) ,
           bp.end() );

 G.reserve( bp.size() );
 for( std::size_t i = 0 ; i + 1 < bp.size() ; ++i ) {
  const double a = bp[ i ] , b = bp[ i + 1 ];
  if( b - a <= 1e-12 )
   continue;
  double pr , sr;
  const double ga = reserve_alloc( t , a , pr , sr , cap );
  const double gb = reserve_alloc( t , b , pr , sr , cap );
  const double slope = ( gb - ga ) / ( b - a );
  G.push_back( { a , b , slope , ga - slope * a } );
  }
 return( G );

 }  // end( ThermalUnitDPSolver::build_reserve_discount )

/*--------------------------------------------------------------------------*/
/*------------------ BUILDING AND SOLVING THE DP PROBLEM -------------------*/
/*--------------------------------------------------------------------------*/

void ThermalUnitDPSolver::build_graph( void )
{
 // read the current fixed status of the Variable: the commitment fixings
 // prune the arcs below, any other fixing makes load_fixings() throw
 load_fixings();

 // first reset any existing graph; note that node labels will be set to 0,
 // which we use as a way to indicate that the node has not been proved
 // reachable from s yet

 // (re)build the ED-solver pool only when the time horizon changes; otherwise
 // the per-instant solvers and their O(n) buffers are kept alive and reused
 // across re-solves, which avoids the (dominant) allocate/free churn
 if( ed_pool_th != time_horizon ) {
  v_on_eds.clear();
  v_on_eds.resize( time_horizon );  // value-initialised to nullptr
  f_start_ed.reset();
  ed_pool_th = time_horizon;
  }

 v_on_nodes.resize( time_horizon );   // no-op unless the horizon changed
 v_off_nodes.resize( time_horizon );

 // reset the nodes in place (their arc storage and the pooled ED solvers
 // survive, so only the cheap bookkeeping is redone)
 reset_node( f_start );
 for( auto & nde : v_on_nodes )
  reset_node( nde );
 for( auto & nde : v_off_nodes )
  reset_node( nde );

 // now rebuild the graph: start from s- - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( init_up_down_time > 0 ) {
  // the unit is already on- - - - - - - - - - - - - - - - - - - - - - - - -

  // s therefore works as an ON-node: take its ED solver from the pool
  if( ! f_start_ed )
   f_start_ed = std::make_unique< DPEDSolver >( 0 , this );
  f_start.DPS = f_start_ed.get();

  // compute kMin, the first time step the unit can be turned OFF due to
  // the need to reach power bound_down[ i ] from the initial power
  // initial_power respecting the ramp-down constraints

  Index kMin = 0;
  for( auto tmp = initial_power ;
       ( kMin < time_horizon ) && ( tmp >= bound_down[ kMin ] + eps ) ; )
   tmp -= delta_ramp_down[ kMin++ ];

  if( kMin < t_init )  // the ramp-down time is less than the time required
   kMin = t_init;      // by the min up-time constraints: use the latter

  if( kMin > time_horizon )  // weird case: the unit must remain on for
   kMin = time_horizon;      // more than the time horizon, i.e., for all
                             // (and only) the time horizon

  // allocate the set of arcs: these are
  //
  //       time_horizon - kMin + 1
  //
  // (note that kMin <= time_horizon, so at least one arc is there)
  // in particular they are ( s , kMin ) (meaning: the unit remains on
  // at 0, 1, 2, ..., kMin - 1 and is off at kMin, and these are kMin
  // instants), ( s , kMin + 1 ), ..., ( s , time_horizon - 1 ),
  // plus there is the final arc ( s , d ).
  //
  // for illustration, consider time_horizon == 6, kMin == 3: the nodes
  // (all OFF ones, so we don't write) are 0, 1, 2, 3, 4, 5, d. the arcs
  // are ( s , 3 ), ( s, 4 ), ( s, 5 ), ( s, d ) These are 6 - 3 + 1 = 4.
  //
  // note the weird case where kMin == 0, i.e., ( s , 0 ) exists, i.e.,
  // the unit is on, but it is immediately turned off: this "oddball"
  // arc corresponds to an empty ED and always has 0 cost

  double fc = 0;             // compute the fixed-cost component of the cost
  Index j = 0;               // this surely comprises the fixed costs
  while( j < kMin )          // between 0 (included) and kMin (excluded)
   fc += const_term[ j++ ];  // since the unit is on in that period

  // the arc ( s , j ) has ON-run [ 0 , j ), which is allowed only if no
  // instant before j is fixed OFF; the ( s , d ) arc (run [ 0 , T )) is
  // allowed only if no instant at all is fixed OFF. The allowed j form a
  // prefix, so the arc list is just truncated
  const Index lim = f_has_fixings ? nxt_off[ 0 ] : time_horizon;
  const Index jend = std::min( time_horizon , lim + 1 );
  const bool darc = ( lim >= time_horizon );

  f_start.v_arcs.resize( ( jend > kMin ? jend - kMin : 0 ) +
			 ( darc ? 1 : 0 ) );
  auto ai = f_start.v_arcs.begin();

  // construct the "normal" arcs up to ( s , jend - 1 )
  for( ; j < jend ; ++j , ++ai ) {
   ai->cost1 = fc;
   ai->cost2 = 0;
   ai->tail = & v_off_nodes[ j ];
   ai->tail->lab = 1;      // mark the tail node as reachable
   fc += const_term[ j ];  // the next fixed cost will comprise that at j
   }

  // now construct the special last arc ( s , d ), if allowed; note that in
  // this case jend == time_horizon, hence the fixed cost from 0 to n - 1
  // (included) has been computed already
  if( darc ) {
   ai->cost1 = fc;
   ai->cost2 = 0;
   ai->tail = & f_end;
   }
  }
 else {  // init_up_down_time <= 0, the unit was off- - - - - - - - - - - -

  // s therefore works as an OFF-node: f_start.DPS must be nullptr
  f_start.DPS = nullptr;

  // allocate the set of arcs: these are
  //
  //       time_horizon - t_init + 1
  //
  // (note that t_init <= time_horizon, so at least one arc is there)
  // where note that t_init == 0 is now possible meaning that
  // init_up_down_time == min_down_time == 0; this implies that the first
  // arc is ( s , 0 ), i.e., "the unit was off at the beginning, but it
  // starts up immediately". Apart from this the structure of the arcs is
  // analogous as in the init_up_down_time > 0 case, except of course they
  // go to the ON nodes

  // the arc ( s , j ) has OFF-run [ 0 , j ), which is allowed only if no
  // instant before j is fixed ON; the ( s , d ) arc (run [ 0 , T )) is
  // allowed only if no instant at all is fixed ON. The allowed j form a
  // prefix, so the arc list is just truncated
  const Index lim = f_has_fixings ? nxt_on[ 0 ] : time_horizon;
  const Index jend = std::min( time_horizon , lim + 1 );
  const bool darc = ( lim >= time_horizon );

  f_start.v_arcs.resize( ( jend > t_init ? jend - t_init : 0 ) +
			 ( darc ? 1 : 0 ) );
  auto ai = f_start.v_arcs.begin();

  // construct the "normal" arcs up to ( s , jend - 1 )
  for( Index j = t_init ; j < jend ; ++j , ++ai ) {
   ai->cost1 = compute_startup_costs( 0 , j );
   ai->cost2 = 0;
   ai->tail = & v_on_nodes[ j ];
   ai->tail->lab = 1;            // mark the tail node as reachable
   }

  // now construct the special last arc ( s , d ), if allowed; note that
  // the fixed cost is 0 because no startup ever happens during the time
  // horizon
  if( darc ) {
   ai->cost1 = 0;
   ai->cost2 = 0;
   ai->tail = & f_end;
   }

  }  // end( else( the unit was off ) )

 // now build the ON and OFF nodes - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // we do this in the order i = 0, 1, ..., n - 1: since the graph is
 // acyclic, if the lab of the node is still 0 when we process it then the
 // node is unreachable from d, and we need not construct any arc

 const Index mut = std::max( min_up_time , Index( 1 ) );
 // min up-time of 0 makes no sense
 const Index mdt = std::max( min_down_time , Index( 1 ) );
 // min down-time of 0 does make sense, but OFF arcs always go forward by
 // at least one time instant, so we pretend that 1 is the minimum value

 for( Index i = 0 ; i < time_horizon ; ++i ) {
  // process ON node ( i , 1 ) - - - - - - - - - - - - - - - - - - - - - - -
  if( v_on_nodes[ i ].lab ) {  // ... but only if it is reachable

   // take the EDSolver of the node from the pool (allocated on first use)
   v_on_nodes[ i ].DPS = get_on_ed( i );

   // allocate the set of arcs, which are:
   //
   // - if i + min_up_time < time_horizon, then
   //
   //       time_horizon - ( i + min_up_time ) + 1
   //
   //   considering that min_up_time >= 1
   //
   //    in particular they are ( i , i + min_up_time ) (meaning: the unit
   //    remains on i, i + 1, ..., i + min_up_time - 1 and is off at
   //    i + min_up_time, and these are min_up_time instants),
   //    ( i , i + min_up_time + 1 ), ..., ( i , time_horizon - 1 ),
   //    plus there is the final arc ( i , d ).
   //
   //    for illustration, consider time_horizon == 6, i = 1, min_up_time = 2
   //    the nodes (all OFF ones, so we don't write) are 0, 1, 2, 3, 4, 5, d.
   //    the arcs are ( 1 , 4 ), ( 1, 5 ), ( 1, d ). These are
   //    6 - ( 1 + 3 ) + 1 = 2.
   //
   // - if, instead, i + min_up_time >= time_horizon, then there only is the
   //   single arc ( i , d ) corresponding to "the unit remains on from i to
   //   the end of the horizon, and it will have to remain on after (but this
   //   is not our concern)
   //
   // the fixed-cost component of the cost of all these arcs surely comprises
   // the fixed costs between i (included) and i + mut (excluded), since the
   // unit is on in that period, save of course if i + mut > time_horizon,
   // in which case it is only the sum up to time_horizon - 1
   double fc = 0;
   Index j = i;
   Index endi = std::min( time_horizon , i + mut );
   while( j < endi )
    fc += const_term[ j++ ];  // since the

   // the arc ( i , j ) has ON-run [ i , j ), which is allowed only if no
   // instant in it is fixed OFF; the ( i , d ) arc (run [ i , T )) is
   // allowed only if no instant >= i is fixed OFF. The allowed j form a
   // prefix, so the arc list is just truncated
   const Index lim = f_has_fixings ? nxt_off[ i ] : time_horizon;
   const Index jend = std::min( time_horizon , lim + 1 );
   const bool darc = ( lim >= time_horizon );

   v_on_nodes[ i ].v_arcs.resize( ( jend > endi ? jend - endi : 0 ) +
				  ( darc ? 1 : 0 ) );
   auto ai = v_on_nodes[ i ].v_arcs.begin();

   // construct the "normal" arcs up to ( i , jend - 1 )
   for( ; j < jend ; ++j , ++ai ) {
    ai->cost1 = fc;
    ai->cost2 = 0;
    ai->tail = & v_off_nodes[ j ];
    ai->tail->lab = 1;      // mark the tail node as reachable
    fc += const_term[ j ];  // the next fixed cost will comprise that at j
    }

   // now construct the special last arc ( i , d ), if allowed; note that in
   // this case jend == time_horizon, hence the fixed cost from i to n - 1
   // (included) has been computed already
   if( darc ) {
    ai->cost1 = fc;
    ai->cost2 = 0;
    ai->tail = & f_end;
    }

   }  // end( if( reached )

  // process OFF node ( i , 0 )- - - - - - - - - - - - - - - - - - - - - - -
  if( v_off_nodes[ i ].lab ) {  // ... but only if it is reachable
   // v_on_nodes[ i ].DPS is and will always remain nullptr here

   // allocate the set of arcs, which are:
   //
   // - if i + min_down_time < time_horizon, then
   //
   //       time_horizon - ( i + min_down_time ) + 1
   //
   //   considering that min_down_time >= 1; note that min_down_time == 0
   //   is in fact possible, but we know that shutting down a unit only to
   //   powering it up again immediately is never a good idea, so we force
   //   down-time periods to be at least of length one. Thus, the structure
   //   of the arcs is analogous as in the ON nodes, except of course they
   //   go to the ON nodes themselves
   //
   // - if, instead, i + min_down_time >= time_horizon, then there only is
   //   the single arc ( i , d ) corresponding to "the unit remains off
   //   from i to the end of the horizon, and it will have to remain off
   //   after (but this is not our concern)

   Index endi = std::min( time_horizon , i + mdt );

   // the arc ( i , j ) has OFF-run [ i , j ), which is allowed only if no
   // instant in it is fixed ON; the ( i , d ) arc (run [ i , T )) is
   // allowed only if no instant >= i is fixed ON. The allowed j form a
   // prefix, so the arc list is just truncated
   const Index lim = f_has_fixings ? nxt_on[ i ] : time_horizon;
   const Index jend = std::min( time_horizon , lim + 1 );
   const bool darc = ( lim >= time_horizon );

   v_off_nodes[ i ].v_arcs.resize( ( jend > endi ? jend - endi : 0 ) +
				   ( darc ? 1 : 0 ) );
   auto ai = v_off_nodes[ i ].v_arcs.begin();

   // construct the "normal" arcs up to ( i , jend - 1 )
   for( Index j = i + mdt ; j < jend ; ++j , ++ai ) {
    ai->cost1 = compute_startup_costs( i , j );
    ai->cost2 = 0;
    ai->tail = & v_on_nodes[ j ];
    ai->tail->lab = 1;            // mark the tail node as reachable
    }

   // now construct the special last arc ( i , d ), if allowed; note that
   // the fixed cost is 0 because no startup ever happens during the time
   // horizon
   if( darc ) {
    ai->cost1 = 0;
    ai->cost2 = 0;
    ai->tail = & f_end;
    }

   }  // end( if( reached )
  }  // end( for( i ) )

 // the graph is now constructed- - - - - - - - - - - - - - - - - - - - - - -
 // nothing needs be done for the destination d

 stage = graph_OK;  // update stage

 }  // end( ThermalUnitDPSolver::build_graph )

/*--------------------------------------------------------------------------*/

void ThermalUnitDPSolver::compute_EDPs( void )
{
 if( stage < graph_OK )
  throw( std::logic_error(
   "ThermalUnitDPSolver::compute_EDPs: graph not ready." ) );

 // precompute the per-period reserve discount g_t (shared by all ON nodes,
 // it depends only on (t)-indexed data): the effective economic-dispatch cost
 // is the quadratic energy cost plus g_t. Empty unless some reserve is
 // rewarded, in which case compute_costs() runs the plain quadratic dispatch.
 // g_disc is the interior discount (cap = max_power) added at every interior
 // period; g_disc_su is the start-up discount (cap = bound_on) added at the
 // first period of an on-interval (the boundary correction). The shut-down
 // variant (cap = bound_down) is built on the fly at the cost readout, where
 // the actual shut-down cap bound_down[k+1] is known.
 if( reserve_rewarded() ) {
  g_disc.resize( time_horizon );
  g_disc_su.resize( time_horizon );
  for( Index t = 0 ; t < time_horizon ; ++t ) {
   g_disc   [ t ] = build_reserve_discount( t , max_power[ t ] );
   g_disc_su[ t ] = build_reserve_discount( t , bound_on  [ t ] );
   }
  }
 else {
  g_disc.clear();
  g_disc_su.clear();
  }

 std::vector< double > cost( time_horizon );

 // update variable costs in the arcs outgoing from s; note that s may have
 // no arcs at all if the fixings forbid them all (the problem is unfeasible)
 if( f_start.DPS && ( ! f_start.v_arcs.empty() ) ) {  // f_start is a ON node
  f_start.DPS->compute_costs( cost );  // solve EDPs, retrieve optimal costs

  // index of first tail node (note: one arc surely exists)
  Index h = h_of_node( f_start.v_arcs.front().tail );

  // the cost of ( s , h ) is found in cost[ h - 1 ]: however, one must
  // be careful of the weird case where ( s , 0 ) is present, i.e.,
  // the unit is on, but it immediately shuts down. this arc has cost 0
  // as "nothing happens there". this is how the arc cost is initialized,
  // and it is never changed, so one just need to skip it
  auto ai = f_start.v_arcs.begin();
  if( h )
   --h;
  else
   ++ai;

  // set the variable costs in the arcs
  for( ; ai != f_start.v_arcs.end() ; ++ai )
   ai->cost2 = cost[ h++ ];
  }

 // update variable costs in the arcs outgoing from every ON( i ): each node
 // has its own ED solver and writes only its own arcs, so the only shared
 // state is the cost scratch, which is made private per worker below; the
 // rest of the solver state is read-only here. solving one ON node is the
 // body of the loop, factored out so the serial and parallel paths share it
 auto solve_on_node = [ this ]( Index i , std::vector< double > & cost ) {
  auto & nd = v_on_nodes[ i ];
  if( nd.v_arcs.empty() )  // unreachable node: nothing to do
   return;

  nd.DPS->compute_costs( cost );  // solve EDPs, retrieve optimal costs

  // the cost of ( i , h ) is found in cost[ h - 1 ]; note that h > i,
  // and therefore h > 0, and therefore h - 1 is well-defined
  Index h = h_of_node( nd.v_arcs.front().tail ) - 1;

  for( auto & a : nd.v_arcs )  // set the variable costs in the arcs
   a.cost2 = cost[ h++ ];
  };

#if TUDPS_PARALLEL
 // resolve the number of workers from intMaxThread (0 = all cores) and go
 // parallel only if it is worth more than one worker and the horizon is large
 // enough for the thread overhead to pay off
 const long nw = f_max_thread > 0 ? long( f_max_thread )
                : std::max< long >( 1 , std::thread::hardware_concurrency() );
 if( ( nw > 1 ) && ( long( time_horizon ) >= long( f_par_min_n ) ) ) {
  // the ParallelFor (hence its worker threads) is created once per solver
  // instance and reused: re-creating it per call would spawn/join threads
  // every time and tank performance. blocking workers (default) sleep idle
  if( ! f_pf )
   f_pf = std::make_unique< ff::ParallelFor >(
           std::max< long >( 1 , std::thread::hardware_concurrency() ) );
  // one cost scratch per worker, reused across calls and iterations
  if( long( f_tcost.size() ) < nw )
   f_tcost.resize( nw );
  for( auto & c : f_tcost )
   if( c.size() < time_horizon )
    c.resize( time_horizon );
  f_pf->parallel_for_thid( 0 , long( time_horizon ) , 1 , 0 ,
   [ & ]( const long i , const int thid ) {
    solve_on_node( Index( i ) , f_tcost[ thid ] );
    } , nw );
  }
 else
#endif
  for( Index i = 0 ; i < time_horizon ; ++i )
   solve_on_node( i , cost );

 stage = edps_OK;  // update stage

 }  // end( ThermalUnitDPSolver::compute_EDPs )

/*--------------------------------------------------------------------------*/

void ThermalUnitDPSolver::min_path( void )
{
 if( stage < edps_OK )
  throw( std::logic_error(
   "ThermalUnitDPSolver::min_path: graph and/or EDPs not ready." ) );

 // reset labels and predecessors for all nodes (except f_start, that will
 // always have lab == 0 and predecessor == nullptr)

 for( auto & nde : v_on_nodes )
  init_node( nde );
 for( auto & nde : v_off_nodes )
  init_node( nde );

 init_node( f_end );

 // now run the shortest path, exploiting the fact that the graph is
 // acyclic and therefore the order s, i = 0, 1, ..., n - 1 for both
 // ON and OFF node is correct

 process_node( f_start );

 for( Index i = 0 ; i < time_horizon ; ++i ) {
  process_node( v_on_nodes[ i ] );
  process_node( v_off_nodes[ i ] );
  }

 // reactive power contribution (AC instances): q[t] in [Qmin,Qmax] is
 // separable from the commitment/active-power shortest path and carries the
 // dualized linear cost reactive_linear_term[t]; the optimal q*[t] is the box
 // endpoint minimising c*q, contributing min(c*Qmin, c*Qmax). Path-independent,
 // so computed once. It is part of the operational cost (the unit holds q only
 // if it exists), hence it is folded in BEFORE the design decision below, and
 // get_var_solution() reports q*[t] = 0 for a non-built unit.
 double Q_star = 0;
 if( ( f_end.lab < TUDPINF ) && ! reactive_linear_term.empty() )
  for( Index t = 0 ; t < time_horizon ; ++t ) {
   const double c = reactive_linear_term[ t ];
   if( c != 0 )
    Q_star += std::min( c * reactive_min[ t ] , c * reactive_max[ t ] );
   }

 // design (investment) decision: the shortest path above solved the
 // operational problem assuming the unit exists (design == 1), so
 // f_end.lab + Q_star is the optimal operational cost p* -- the active-power
 // schedule plus the reactive contribution, both available only if the unit is
 // built. Building costs design_cost on top; it is worth building iff
 // p* + design_cost <= 0, otherwise the unit is not built and contributes
 // nothing (cost 0, zero schedule, zero reactive). A designable unit is always
 // initially off (ThermalUnitBlock forbids InitUpDownTime >= 0 with an
 // investment cost) so it can always stay off the whole horizon and the design
 // is never forced on. Generalises to an integer design by the same threshold
 // argument; the continuous case does not apply (binary commitments inside).
 if( has_design ) {
  // a commitment fixed ON (or the design variable fixed to 1) forces the
  // unit to be built regardless of the economics; if it cannot (no feasible
  // path, or the design variable is fixed to 0) the problem is unfeasible,
  // since the all-zero "not built" solution violates the fixings
  if( ( ! f_no_build ) && ( f_end.lab < TUDPINF ) &&
      ( f_must_build || ( f_end.lab + Q_star + design_cost <= 0 ) ) ) {
   design_on = true;
   f_end.lab += Q_star + design_cost;
   }
  else if( f_must_build ) {
   design_on = false;
   f_end.lab = TUDPINF;    // unfeasible: must be built, but cannot
   f_end.pred = nullptr;
   }
  else {
   design_on = false;
   f_end.lab = 0;          // not built: the unit is absent
   }
  }
 else
  f_end.lab += Q_star;     // no design: the reactive term is always incurred

 stage = path_OK;  // all done: update stage

 }  // end( ThermalUnitDPSolver::min_path )

/*--------------------------------------------------------------------------*/

void ThermalUnitDPSolver::compute_solutions( void )
{
 if( stage < edps_OK )
  throw( std::logic_error(
   "ThermalUnitDPSolver::compute_solutions: graph and/or path not ready." ) );

 std::fill( P.begin() , P.end() , 0 );
 std::fill( U.begin() , U.end() , false );

 Index k = time_horizon;
 auto n = f_end.pred;

 if( ! n ) {        // no path to the destination: the fixings made the
  stage = sol_OK;   // problem unfeasible, there is no solution to compute
  return;
  }

 // compute the solution by visiting the optimal path backward from f_end

 do {
  Index h = h_of_node( n );  // the current arc is ( h , k )
  if( n->DPS && k ) {
   // n is ON( h ), or the source (if h == 0) that works as an ON node
   // the power and commitment variables of this arc are these with index
   // h, ..., k - 1, comprised if n == f_start (this is why h_of_node()
   // returns 0 for it); however, one has to explicitly avoid the special
   // case of the "empty" arc ( s , 0 ) that has no power and commitment
   // variables
   // get optimal values of power variables out of the EDSolver
   n->DPS->compute_power_variables( k - 1 , P );
   for( Index i = h ; i < k ; )  // set all commitment variables to true
    U[ i++ ] = true;
   }
  // else n is OFF( h ), or the source (if h == 0) that works as an OFF
  // node: P[ i ] = U[ i ] = 0 for i = h, ..., k - 1, but these already
  // have those values

  k = h;        // the previous beginning will be the end
  n = n->pred;  // back one arc

  } while( n );  // ... until we hit f_start that has pred == nullptr

 stage = sol_OK;  // all done: update stage

 }  // end( ThermalUnitDPSolver::compute_solutions )

/*--------------------------------------------------------------------------*/
/*---------------------- PRIVATE METHODS OF THE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

void ThermalUnitDPSolver::load_parameters( void )
{
 // locking the Block
 bool owned = f_Block->is_owned_by( f_id );
 if( ( ! owned ) && ( ! f_Block->read_lock() ) )
  throw( std::runtime_error(
   "ThermalUnitDPSolver::load_parameters: unable to lock the Block." ) );

 // casting has been checked in set_Block() already
 auto b = static_cast< ThermalUnitBlock * >( f_Block );

 // scalar values
 time_horizon = b->get_time_horizon();
 init_up_down_time = b->get_init_up_down_time();
 min_up_time = b->get_min_up_time();
 min_down_time = b->get_min_down_time();
 initial_power = b->get_initial_power();

 // t_init: first instant in which a decision can be made, as all the
 //         instants before are "blocked" by the initial conditions
 if( init_up_down_time > 0 )
  if( min_up_time > init_up_down_time )
   t_init = std::min( time_horizon , min_up_time - init_up_down_time );
  else
   t_init = 0;
 else
  if( min_down_time > - init_up_down_time )
   t_init = std::min( time_horizon , min_down_time + init_up_down_time );
  else
   t_init = 0;

 // power vectors
 startup_costs = b->get_start_up_cost();
 min_power = b->get_min_power();
 max_power = b->get_max_power();
 bound_on = b->get_start_up_limit();
 bound_down = b->get_shut_down_limit();

 if( b->get_delta_ramp_up().empty() )
  delta_ramp_up = max_power;
 else
  delta_ramp_up = b->get_delta_ramp_up();

 if( b->get_delta_ramp_down().empty() )
  delta_ramp_down = max_power;
 else
  delta_ramp_down = b->get_delta_ramp_down();

 retrieve_term( quad_term , b->get_quad_term() );
 retrieve_term( linear_term , b->get_linear_term() );
 retrieve_term( const_term , b->get_const_term() );

 // spinning reserve: participation factors (caps in pr<=rho_p*p / sr<=rho_s*p)
 // and objective cost coefficients. The cost defaults to the participation
 // factor but may carry a (possibly negative) Lagrangian price, so it is read
 // from the separate cost getter. All empty when the reserve is absent.
 primary_rho            = b->get_primary_rho();
 secondary_rho          = b->get_secondary_rho();
 primary_reserve_cost   = b->get_primary_spinning_reserve_cost();
 secondary_reserve_cost = b->get_secondary_spinning_reserve_cost();
 if( primary_reserve_cost.empty() )
  primary_reserve_cost = primary_rho;
 if( secondary_reserve_cost.empty() )
  secondary_reserve_cost = secondary_rho;

 // reactive power (AC instances): read the box [Qmin,Qmax] and the dualized
 // linear cost coefficient on q[t]. q[t] is separable from the DP, so run_DP()
 // adds its optimal contribution as a constant. Empty reactive_linear_term ->
 // the unit has no reactive power and the term is skipped.
 if( b->get_reactive_power( 0 ) ) {
  retrieve_term( reactive_linear_term , b->get_reactive_linear_term() );
  reactive_min.resize( time_horizon );
  reactive_max.resize( time_horizon );
  for( Index t = 0 ; t < time_horizon ; ++t ) {
   reactive_min[ t ] = b->get_min_reactive_power( t );
   reactive_max[ t ] = b->get_max_reactive_power( t );
   }
  }
 else
  reactive_linear_term.clear();

 // design (investment): present iff the unit carries a nonzero investment cost.
 // The DP solves the operational problem assuming the unit exists; min_path()
 // then decides whether to build it. See min_path() for the threshold rule.
 // NOTE: the actual building cost is the *current* coefficient of the design
 // variable in the Objective (get_design_cost()), which may differ from the
 // structural investment cost if a dualizing Solver has changed it; this is
 // kept in sync via the eSetInvCost Modification.
 has_design  = ( b->get_investment_cost() != 0 );
 design_cost = b->get_design_cost();
 design_on   = false;

 // unlock the Block
 if( ! owned )
  f_Block->read_unlock();

 v_on_nodes.resize( time_horizon );
 v_off_nodes.resize( time_horizon );
 P.resize( time_horizon );
 U.resize( time_horizon );
 stage = start;

 }  // end( ThermalUnitDPSolver::load_parameters )

/*--------------------------------------------------------------------------*/

void ThermalUnitDPSolver::load_fixings( void )
{
 f_has_fixings = f_must_build = f_no_build = false;

 // locking the Block
 bool owned = f_Block->is_owned_by( f_id );
 if( ( ! owned ) && ( ! f_Block->read_lock() ) )
  throw( std::runtime_error(
   "ThermalUnitDPSolver::load_fixings: unable to lock the Block." ) );

 auto b = static_cast< ThermalUnitBlock * >( f_Block );

 // commitment fixings, translated into the "first instant >= t fixed
 // OFF / ON" tables that build_graph() uses to prune the arcs
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
     std::string( "ThermalUnitDPSolver: fixed " ) + name +
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

 }  // end( ThermalUnitDPSolver::load_fixings )

/*--------------------------------------------------------------------------*/

void ThermalUnitDPSolver::process_modifications( void )
{
 bool reload = false;

 // note: since processing the Modification is fast, we don't bother with
 // being nice to other processes and do it all with v_mod under lock
  // try to acquire lock, spin on failure
 while( f_mod_lock.test_and_set( std::memory_order_acquire ) )
  ;

 // process all the Modifications
 for( auto mod : v_mod )
  if( guts_of_process_modifications( mod.get() ) ) {
   reload = true;  // a reset must be done
   break;          // ignore all the remaining Modifications
   }

 v_mod.clear();  // all Modifications tackled, clear the list

 f_mod_lock.clear( std::memory_order_release );  // release lock

 if( reload )
  load_parameters();

 }  // end( ThermalUnitDPSolver::process_modifications )

/*--------------------------------------------------------------------------*/

bool ThermalUnitDPSolver::guts_of_process_modifications( const p_Mod mod )
{
 // NBModification
 if( dynamic_cast< NBModification * >( mod ) )
  return( true );

 // GroupModification
 if( const auto gm = dynamic_cast< GroupModification * >( mod ) ) {
  bool reload = false;
  for( const auto & submod : gm->sub_Modifications() )
   if( guts_of_process_modifications( submod.get() ) )
    reload = true;

  return( reload );
  }

 // ThermalUnitBlockMod
 if( const auto tubm = dynamic_cast< ThermalUnitBlockMod * >( mod ) ) {
   auto b = static_cast< ThermalUnitBlock * >( f_Block );

   switch( tubm->type() ) {
    case( ThermalUnitBlockMod::eSetMaxP ):
     max_power = b->get_max_power();
     if( stage > graph_OK )
      stage = graph_OK;
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
      if( min_up_time > init_up_down_time )
       t_init = std::min( time_horizon , min_up_time - init_up_down_time );
      else
       t_init = 0;
     else
      if( min_down_time > - init_up_down_time )
       t_init = std::min( time_horizon , min_down_time + init_up_down_time );
      else
       t_init = 0;
     stage = start;
     return( false );

    case( ThermalUnitBlockMod::eSetAv ):
     return( true );  // TODO

    case( ThermalUnitBlockMod::eSetSUC ):
     startup_costs = b->get_start_up_cost();
     if( stage > edps_OK )
      stage = edps_OK;
     return( false );

    case( ThermalUnitBlockMod::eSetLinT ):
     retrieve_term( linear_term , b->get_linear_term() );
     if( stage > graph_OK )
      stage = graph_OK;
     return( false );

    case( ThermalUnitBlockMod::eSetQuadT ):
     retrieve_term( quad_term , b->get_quad_term() );
     if( stage > graph_OK )
      stage = graph_OK;
     return( false );

    case( ThermalUnitBlockMod::eSetConstT ):
     retrieve_term( const_term , b->get_const_term() );
     stage = start;
     return( false );

    case( ThermalUnitBlockMod::eSetPrSpResCost ):
     primary_reserve_cost = b->get_primary_spinning_reserve_cost();
     if( primary_reserve_cost.empty() )
      primary_reserve_cost = primary_rho;
     if( stage > graph_OK )
      stage = graph_OK;  // reserve price feeds the per-period economic dispatch
     return( false );

    case( ThermalUnitBlockMod::eSetSecSpResCost ):
     secondary_reserve_cost = b->get_secondary_spinning_reserve_cost();
     if( secondary_reserve_cost.empty() )
      secondary_reserve_cost = secondary_rho;
     if( stage > graph_OK )
      stage = graph_OK;
     return( false );

    case( ThermalUnitBlockMod::eSetReactiveLinT ):
     // reactive price only changes the separable constant added at the end of
     // run_DP(); a full re-run recomputes it (the graph/dispatch are unaffected)
     retrieve_term( reactive_linear_term , b->get_reactive_linear_term() );
     stage = start;
     return( false );

    case( ThermalUnitBlockMod::eSetInvCost ): {
     // the cost of the design (investment) variable has (possibly) changed,
     // e.g. because a dualizing Solver pushed a multiplier into its Objective
     // coefficient. Refresh our copy from get_design_cost(); only invalidate
     // the cached state when the value actually changed, since this is issued
     // whenever the Objective is touched (the design index is always in range).
     auto nc = b->get_design_cost();
     if( nc != design_cost ) {
      design_cost = nc;
      stage = start;
      }
     return( false );
     }

    case( ThermalUnitBlockMod::eFixVars ):
     // the fixed status of some Variable changed: the graph has to be
     // rebuilt from scratch, as build_graph() re-reads the fixings (via
     // load_fixings()) and prunes the incompatible arcs
     stage = start;
     return( false );

    }  // end( switch )

  return( true );

  }  // end( ThermalUnitBlockMod )

 return( false );  // any other Modification: I assume it's harmless

 }  // end( ThermalUnitDPSolver::guts_of_process_modifications )

/*--------------------------------------------------------------------------*/

void ThermalUnitDPSolver::retrieve_term( std::vector< double > & out ,
                                         const std::vector< double > & in ) const
{
 if( in.empty() ) {
  out.resize( time_horizon );
  std::fill( out.begin() , out.end() , 0 );
  return;
  }

 if( in.size() == 1 ) {
  out.resize( time_horizon );
  std::fill( out.begin() , out.end() , in[ 0 ] );
  return;
  }

 out = in;
 }

/*--------------------------------------------------------------------------*/

ThermalUnitDPSolver::DPEDSolver * ThermalUnitDPSolver::get_on_ed( Index i )
{
 if( ! v_on_eds[ i ] )
  v_on_eds[ i ] = std::make_unique< DPEDSolver >( i , this );
 return( v_on_eds[ i ].get() );
 }

/*--------------------------------------------------------------------------*/
/*----------- METHODS OF ThermalUnitDPSolver::DPEDSolver -------------------*/
/*--------------------------------------------------------------------------*/

ThermalUnitDPSolver::DPEDSolver::DPEDSolver( Index h ,
                                             ThermalUnitDPSolver * s )
 : EDSolver( h , s )
{
 auto & time_horizon = f_solver->time_horizon;

 // when the unit carries spinning-reserve data the economic dispatch adds the
 // piecewise-linear reserve discount g_t to the per-period cost, splitting the
 // value-function pieces at g_t's (at most ~5 internal) breakpoints; doubling
 // each ping-pong half is therefore more than enough. f_rmul == 1 reproduces
 // the plain quadratic dispatch exactly (g_t absent)
 f_rmul = ( ( ! f_solver->primary_rho.empty() ) ||
            ( ! f_solver->secondary_rho.empty() ) ) ? 2 : 1;

 #if( COMPUTE_DUALS )
  Index coeffsize = time_horizon * time_horizon +f_h * f_h -
                    2 * f_h * time_horizon;
 #else
  // Index coeffsize = 4 * ( time_horizon - f_h + 1 );
  // the theory says it should work, but it does not
  Index coeffsize = f_rmul * 4 * time_horizon;
 #endif
 if( coeffsize != coeffs.size() ) {
  coeffs.resize( coeffsize );
  #if( COMPUTE_DUALS )
   m.resize( coeffsize + time_horizon - f_h );
   v.resize( time_horizon );
   pos.resize( time_horizon );
  #else
   m.resize( coeffsize + 2 );
   v.resize( 2 );
   pos.resize( 2 );
  #endif
  unc_p.resize( time_horizon );
  con_p.resize( time_horizon );
  }
 }

/*--------------------------------------------------------------------------*/

void ThermalUnitDPSolver::DPEDSolver::compute_costs(
 std::vector< double > & costs )
{
 // scalar values
 auto time_horizon = f_solver->time_horizon;
 auto init_up_down_time = f_solver->init_up_down_time;
 auto initial_power = f_solver->initial_power;

 // power vectors
 const auto & min_power = f_solver->min_power;
 const auto & max_power = f_solver->max_power;
 const auto & delta_ramp_up = f_solver->delta_ramp_up;
 const auto & delta_ramp_down = f_solver->delta_ramp_down;
 const auto & bound_on = f_solver->bound_on;
 const auto & bound_down = f_solver->bound_down;

 // coefficients of the objective function
 const auto & quad_term = f_solver->quad_term;
 const auto & linear_term = f_solver->linear_term;

 Index k = f_h;

 coeffs[ 0 ].alfa = quad_term[ k ];
 coeffs[ 0 ].beta = linear_term[ k ];
 coeffs[ 0 ].gamma = 0;

 /* Initialize the vector m containing the endpoints of the pieces.
  * At first, it contains the two endpoints of the individual piece.
  * At startup, power can't exceed the bound-on value \barl_k.
  * However, if the unit is on at the beginning of the time horizon with
  * the given initial value initial_power, then the interval is restricted
  * to take it into account. */

 if( ( f_h == 0 ) && ( init_up_down_time > 0 ) ) {
  m[ 0 ] = std::max( min_power[ k ] , initial_power - delta_ramp_down[ k ] );
  m[ 1 ] = std::min( max_power[ k ] , initial_power + delta_ramp_up[ k ] );
  }
 else {
  m[ 0 ] = min_power[ k ];
  m[ 1 ] = std::min( bound_on[ k ] , max_power[ k ] );  // \bar{l}_k
  }

 #if ( COMPUTE_DUALS )
  Index mcnt = 2;
  Index coeffcnt = 1;
  v[ k ] = 0;
  pos[ k ].begm = 0;
  pos[ k ].begt = 0;
 #else
  Index coeffcnt = f_rmul * 2 * ( time_horizon - f_h );
  // next free position in coeffs[]
  Index mcnt = f_rmul * 2 * ( time_horizon - f_h ) + 1;
  // next free position in m[]
  Index nextk = 1;     // next free position in v[], pos[]
  // since there are only two positions, nextk ping-pongs between 1 and 0

  v[ 0 ] = 0;
  // initialize the vector pos containing the initial indices of the pieces.
  pos[ 0 ].begm = 0;
  pos[ 0 ].begt = 0;
 #endif

 // initialize the vector of unconstrained power values, i.e., power values
 // are not constrained by bound_down[ k ]. When reserves are rewarded the
 // single quadratic piece is first augmented with the reserve discount g_k,
 // turning it into a convex piecewise-quadratic; the optimum and cost are then
 // read from the augmented pieces
 if( ! f_solver->g_disc.empty() ) {
  Index mtmp = 2 , ctmp = 1;
  // k == f_h is the first period of the on-interval: a start-up, so its reserve
  // band uses the start-up cap bound_on (g_disc_su). A single-period interval
  // [h,h] is both start-up and shut-down; like the run-length solver it keeps
  // the start-up band (no shut-down swap at the base case).
  unc_p[ k ] = augment_with_g( f_solver->g_disc_su[ k ] , 0 , 0 ,
                               mtmp , ctmp , v[ 0 ] );

  if( ( k < time_horizon - 1 ) && ( unc_p[ k ] > bound_down[ k + 1 ] ) )
   con_p[ k ] = bound_down[ k + 1 ];
  else
   con_p[ k ] = unc_p[ k ];

  int qq = 0;
  while( ( qq < v[ 0 ] ) && ( con_p[ k ] > m[ qq + 1 ] ) )
   ++qq;
  costs[ k ] = coeffs[ qq ].alfa * con_p[ k ] * con_p[ k ] +
               coeffs[ qq ].beta * con_p[ k ] + coeffs[ qq ].gamma;
  }
 else {
  if( std::abs( coeffs[ 0 ].alfa ) <= 1e-16 )
   unc_p[ k ] = ( coeffs[ 0 ].beta <= 0 ? m[ 1 ] : m[ 0 ] );
  else
   unc_p[ k ] = std::min( m[ 1 ] ,
			  std::max( m[ 0 ] ,
				    -coeffs[ 0 ].beta / ( 2 * coeffs[ 0 ].alfa )
				    ) );

  /* Initialize the vector of constrained power values, that will be
   * computed at each iteration.
   * Constrained means that they must be <= bound_down[ k ] */

  if( ( k < time_horizon - 1 ) && ( unc_p[ k ] > bound_down[ k + 1 ] ) )
   con_p[ k ] = bound_down[ k + 1 ];
  else
   con_p[ k ] = unc_p[ k ];

  costs[ k ] = coeffs[ 0 ].alfa * con_p[ k ] * con_p[ k ] +
               coeffs[ 0 ].beta * con_p[ k ];
  }

 // --- shut-down boundary correction helpers (only used when reserves are
 // rewarded) ----------------------------------------------------------------
 // gpart(g,p): slope and intercept of the piecewise-linear g at p (0 outside).
 auto gpart = []( const std::vector< g_piece > & g , double p ,
                  double & sl , double & ic ) {
  sl = ic = 0;
  for( const auto & gp : g )
   if( ( p >= gp.lo - 1e-9 ) && ( p <= gp.hi + 1e-9 ) ) {
    sl = gp.slope; ic = gp.intercept; return;
    }
  };
 // shutdown_min: min over [lo,hi] of z(p) - g_int(p) + g_sd(p), where z is the
 // piecewise-quadratic coeffs[begt+i] on [m[begm+i],m[begm+i+1]] (i < np) with
 // the interior g_int already baked in. Swaps the interior reserve band for the
 // shut-down one g_sd; read-only on the sweep buffers.
 auto shutdown_min = [ & ]( Index begm , Index begt , int np ,
                            double lo , double hi ,
                            const std::vector< g_piece > & g_sd ,
                            const std::vector< g_piece > & g_int ) -> double {
  std::vector< double > bp = { lo , hi };
  for( int i = 0 ; i <= np ; ++i ) {
   const double x = m[ begm + i ];
   if( ( x > lo + 1e-12 ) && ( x < hi - 1e-12 ) ) bp.push_back( x );
   }
  for( const auto & gp : g_sd )
   if( ( gp.hi > lo + 1e-12 ) && ( gp.hi < hi - 1e-12 ) ) bp.push_back( gp.hi );
  for( const auto & gp : g_int )
   if( ( gp.hi > lo + 1e-12 ) && ( gp.hi < hi - 1e-12 ) ) bp.push_back( gp.hi );
  std::sort( bp.begin() , bp.end() );
  bp.erase( std::unique( bp.begin() , bp.end() ,
             []( double a , double b ) { return( b - a <= 1e-12 ); } ) ,
            bp.end() );
  double best = TUDPINF;
  for( std::size_t j = 0 ; j + 1 < bp.size() ; ++j ) {
   const double a = bp[ j ] , b = bp[ j + 1 ] , mid = 0.5 * ( a + b );
   int zi = 0;
   while( ( zi < np - 1 ) && ( mid > m[ begm + zi + 1 ] ) ) ++zi;
   const double a2 = coeffs[ begt + zi ].alfa ,
                b2 = coeffs[ begt + zi ].beta ,
                c2 = coeffs[ begt + zi ].gamma;
   double ssl , sic , isl , iic;
   gpart( g_sd , mid , ssl , sic ); gpart( g_int , mid , isl , iic );
   const double bb = b2 + ssl - isl , cc = c2 + sic - iic;
   const double pstar = ( std::abs( a2 ) <= 1e-16 )
                      ? ( bb <= 0 ? b : a )
                      : std::min( b , std::max( a , -bb / ( 2 * a2 ) ) );
   const double val = a2 * pstar * pstar + bb * pstar + cc;
   if( val < best ) best = val;
   }
  return( best );
  };

 // outermost loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 for( ; ++k < time_horizon ; ) {
  /* Building pieces: \bar{m}_0 is the first endpoint of the first piece of
   * the z_{hk}(\bar{p}) objective function. Such endpoint will be saved in
   * the m vector. */

  #if( COMPUTE_DUALS )
   pos[ k ].begm = mcnt;
   pos[ k ].begt = coeffcnt;

   m[ mcnt ] = std::max( min_power[ k ] ,
			 m[ pos[ k - 1 ].begm ] - delta_ramp_down[ k - 1 ] );
  #else
   pos[ nextk ].begm = mcnt;
   pos[ nextk ].begt = coeffcnt;

   m[ mcnt ] = std::max( min_power[ k ] ,
			 m[ pos[ 1 - nextk ].begm ] - delta_ramp_down[ k - 1 ]
			 );
  #endif

  double p_bar = m[ mcnt ];  // \bar{m}_0
  Index v_bar = 0;           // after the case 3 will contain v[ k ]

  // compute q, the index of the piece where p^*(\bar{p}) belongs.

  double pstar;  // p^*(\bar{p})

  if( p_bar < unc_p[ k - 1 ] )
   pstar = std::min( unc_p[ k - 1 ] , p_bar + delta_ramp_down[ k - 1 ] );
  else
   pstar = std::max( unc_p[ k - 1 ] , p_bar - delta_ramp_up[ k - 1 ] );

  #if( COMPUTE_DUALS )
   Index qm = pos[ k - 1 ].begm;

   while( ( pstar >= m[ qm + 1 ] ) && ( qm < pos[ k ].begm - 2 ) )
    ++qm;

   Index q = qm - pos[ k - 1 ].begm + pos[ k - 1 ].begt;

   // compute the last endpoint of the piece, \bar{u}
   double u_bar = std::min( max_power[ k ] ,
			    m[ mcnt - 1 ] + delta_ramp_up[ k - 1 ] );
  #else
   Index qm = pos[ 1 - nextk ].begm;
   /*!!
   Index poslim = pos[ 1 - nextk ].begm + ( v[ 1 - nextk ] + 1 ) + 1 - 2;

   while( ( pstar >= m[ qm + 1 ] ) && ( qm < poslim ) )
    ++qm;
    !!*/

   if( ( v[ 1 - nextk ] >= 0 ) ||
       ( pos[ 1 - nextk ].begm > - v[ 1 - nextk ] ) ) {
    Index poslim = pos[ 1 - nextk ].begm + v[ 1 - nextk ];

    while( ( pstar >= m[ qm + 1 ] ) && ( qm < poslim ) )
     ++qm;
    }

   Index q = qm - pos[ 1 - nextk ].begm + pos[ 1 - nextk ].begt;

   // compute the last endpoint of the piece, \bar{u}
   double u_bar = std::min( max_power[ k ] ,
			    m[ pos[ 1 - nextk ].begm + v[ 1 - nextk ] + 1 ]
			    + delta_ramp_up[ k - 1 ] );
  #endif
  ++mcnt;

  bool firstTime = true;

  // CASE 1- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  while( unc_p[ k - 1 ] > p_bar + delta_ramp_down[ k - 1 ] + f_solver->eps )
  {
   // set coeffs fields to compute \bar{z}^{\bar{v}}(p)

   coeffs[ coeffcnt ].alfa = quad_term[ k ] + coeffs[ q ].alfa;
   coeffs[ coeffcnt ].beta = linear_term[ k ] + coeffs[ q ].beta +
                             2 * delta_ramp_down[ k - 1 ] * coeffs[ q ].alfa;
   coeffs[ coeffcnt ].gamma = coeffs[ q ].gamma +
    coeffs[ q ].alfa * delta_ramp_down[ k - 1 ] * delta_ramp_down[ k - 1 ] +
    coeffs[ q ].beta * delta_ramp_down[ k - 1 ];

   /* Compute the maximum value for \bar{p} such that:
    * - p^*_k(\bar{p}) stays in the q-th interval;
    * - unc_p stays out of the admissible range;
    * - \bar{p} stays admissible. */

   if( m[ qm + 1 ] - delta_ramp_down[ k - 1 ] <
       unc_p[ k - 1 ] - delta_ramp_down[ k - 1 ] - f_solver->eps ) {
    p_bar = m[ qm + 1 ] - delta_ramp_down[ k - 1 ];
    ++q;
    ++qm;
    }
   else
    p_bar = unc_p[ k - 1 ] - delta_ramp_down[ k - 1 ];

   if( p_bar > u_bar )
    p_bar = u_bar;

   ++v_bar;
   m[ mcnt++ ] = p_bar;

   // compute unc_p, unconstrained optimal value for z_{hk}

   if( firstTime &&
       ( 2 * coeffs[ coeffcnt ].alfa * p_bar + coeffs[ coeffcnt ].beta > 0 ) ) {
    if( std::abs( coeffs[ coeffcnt ].alfa ) <= 1e-16 ) {
     if( coeffs[ coeffcnt ].beta >= 0 )
      unc_p[ k ] = m[ mcnt - 2 ];
     // else do nothing, the function is still decreasing in the next interval
     }
    else
     unc_p[ k ] = std::max( m[ mcnt - 2 ] ,
	        - coeffs[ coeffcnt ].beta / ( 2 * coeffs[ coeffcnt ].alfa ) );
    firstTime = false;
    }

   ++coeffcnt;

   }  // end( while( CASE 1 ) )

  // CASE 2- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if( unc_p[ k - 1 ] >= p_bar - delta_ramp_up[ k - 1 ] ) {

   // set coeffs fields to compute \bar{z}^{\bar{v}}(p)

   coeffs[ coeffcnt ].alfa = quad_term[ k ];
   coeffs[ coeffcnt ].beta = linear_term[ k ];
   coeffs[ coeffcnt ].gamma =
    coeffs[ q ].alfa * unc_p[ k - 1 ] * unc_p[ k - 1 ] +
    coeffs[ q ].beta * unc_p[ k - 1 ] + coeffs[ q ].gamma;

   /* compute the maximum value for \bar{p} such that:
    * - unc_p stays out of the admissible range;
    * - \bar{p} stays admissible. */

   p_bar = std::min( u_bar , unc_p[ k - 1 ] + delta_ramp_up[ k - 1 ] );

   ++v_bar;
   m[ mcnt++ ] = p_bar;

   if( firstTime &&
       ( 2 * coeffs[ coeffcnt ].alfa * p_bar + coeffs[ coeffcnt ].beta > 0 ) ) {
    if( std::abs( coeffs[ coeffcnt ].alfa ) <= 1e-16 ) {
     if( coeffs[ coeffcnt ].beta >= 0 )
      unc_p[ k ] = m[ mcnt - 2 ];
     // else do nothing, the function is still decreasing in the next interval
     }
    else
     unc_p[ k ] = std::max( m[ mcnt - 2 ] ,
		- coeffs[ coeffcnt ].beta / ( 2 * coeffs[ coeffcnt ].alfa ) );
    firstTime = false;
    }

   ++coeffcnt;

   }  // end( if( CASE 2 ) )

  // CASE 3- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  while( p_bar < u_bar ) {
   // set coeffs fields to compute \bar{z}^{\bar{v}}(p)

   coeffs[ coeffcnt ].alfa = quad_term[ k ] + coeffs[ q ].alfa;
   coeffs[ coeffcnt ].beta = linear_term[ k ] + coeffs[ q ].beta -
                              2 * delta_ramp_up[ k - 1 ] * coeffs[ q ].alfa;
   coeffs[ coeffcnt ].gamma = coeffs[ q ].gamma +
    coeffs[ q ].alfa * delta_ramp_up[ k - 1 ] * delta_ramp_up[ k - 1 ] -
    coeffs[ q ].beta * delta_ramp_up[ k - 1 ];

   /* Compute the maximum value for \bar{p} such that:
    * - p^*_k(\bar{p}) stays in the q-th interval;
    * - \bar{p} stays admissible. */

   p_bar = std::min( m[ qm + 1 ] + delta_ramp_up[ k - 1 ] , u_bar );

   ++v_bar;
   m[ mcnt++ ] = p_bar;
   ++q;
   ++qm;

   if( firstTime &&
       ( 2 * coeffs[ coeffcnt ].alfa * p_bar + coeffs[ coeffcnt ].beta > 0 ) ) {
    if( std::abs( coeffs[ coeffcnt ].alfa ) <= 1e-16 ) {
     if( coeffs[ coeffcnt ].beta >= 0 )
      unc_p[ k ] = m[ mcnt - 2 ];
     // else do nothing, the function is still decreasing in the next interval
     }
    else
     unc_p[ k ] = std::max( m[ mcnt - 2 ] ,
		- coeffs[ coeffcnt ].beta / ( 2 * coeffs[ coeffcnt ].alfa ) );
    firstTime = false;
    }

   ++coeffcnt;

   }  // end( while( CASE 3 ) )
   // end of the tree cases - - - - - - - - - - - - - - - - - - - - - - - - -

  #if( COMPUTE_DUALS )
   v[ k ] = v_bar - 1;
  #else
   v[ nextk ] = v_bar - 1;
  #endif

  if( firstTime )  // function is strictly decreasing
   unc_p[ k ] = u_bar;

  // reserves rewarded: augment z_{h,k} with the discount g_k (splitting the
  // pieces just built at g_k's breakpoints) and recompute the unconstrained
  // optimum; the con_p / costs computation below then reads the augmented
  // pieces. Only the (active) COMPUTE_DUALS == 0 ping-pong layout is supported
  #if( ! COMPUTE_DUALS )
   if( ! f_solver->g_disc.empty() )
    unc_p[ k ] = augment_with_g( f_solver->g_disc[ k ] ,
                                 Index( pos[ nextk ].begm ) ,
                                 Index( pos[ nextk ].begt ) ,
                                 mcnt , coeffcnt , v[ nextk ] );
  #endif

  // shut-down boundary correction: the cost of closing the on-interval at k
  // (cost arc to shut down at k) must price the reserve under the shut-down cap
  // bound_down[k+1], not the interior max_power baked into z_{h,k}. Compute it
  // here on the current (pre-flip) pieces; it overrides costs[k] below. No-op at
  // the horizon end (k == n-1: a tail, no shut-down) or when the cap is not
  // tighter than max_power. The continuing on-run keeps the interior g_k.
  bool sd_use = false; double sd_cost = 0;
  #if( ! COMPUTE_DUALS )
  if( ( ! f_solver->g_disc.empty() ) && ( k < time_horizon - 1 ) &&
      ( bound_down[ k + 1 ] < max_power[ k ] - 1e-12 ) ) {
   const Index  bm = pos[ nextk ].begm , bt = pos[ nextk ].begt;
   const int    nptmp = v[ nextk ];
   const double sdlo = m[ bm ] , sdhi = std::min( double( bound_down[ k + 1 ] ) ,
                                                  m[ bm + nptmp ] );
   if( sdhi > sdlo + 1e-12 ) {
    const auto g_sd = f_solver->build_reserve_discount( k , bound_down[ k + 1 ] );
    sd_cost = shutdown_min( bm , bt , nptmp , sdlo , sdhi , g_sd ,
                            f_solver->g_disc[ k ] );
    sd_use = true;
    }
   }
  #endif

  // compute con_p[ k ], constrained optimal value for the entire function.
  // important note: the case where k == time_horizon - 1 is dealt with
  // in a special way, i.e., by not requiring the power to be at the level
  // that would be required for the unit to stop. this means that a less
  // constrained problem is solved, resulting in a smaller value. this is
  // because there is no point in forcing the unit to shut down at the end
  // of the time instant, since what happens after that is irrelevant to
  // the problem we are solving

  if( ( k < time_horizon - 1 ) && ( unc_p[ k ] > bound_down[ k + 1 ] ) )
   con_p[ k ] = bound_down[ k + 1 ];
  else
   con_p[ k ] = unc_p[ k ];

  // compute the cost for the node (h,k) in costs[]

  #if( COMPUTE_DUALS )
    qm = pos[ k ].begm;
  #else
    qm = pos[ nextk ].begm;
  #endif

  //?? while( ( con_p[ k ] > m[ qm + 1 ] ) && ( m[ qm + 1 ] != 0 ) )
  while( con_p[ k ] > m[ qm + 1 ] )
   ++qm;

  #if( COMPUTE_DUALS )
   q = qm - pos[ k ].begm + pos[ k ].begt;
  #else
   q = qm - pos[ nextk ].begm + pos[ nextk ].begt;
   nextk = 1 - nextk;
   coeffcnt = nextk * f_rmul * ( 2 * time_horizon - f_h );
   mcnt     = nextk * ( f_rmul * ( 2 * time_horizon - f_h ) + 1 );
  #endif

  costs[ k ] = sd_use ? sd_cost
             : ( coeffs[ q ].alfa * con_p[ k ] * con_p[ k ] +
                 coeffs[ q ].beta * con_p[ k ] + coeffs[ q ].gamma );

  }  // end( for( k ) )
 }  // end( ThermalUnitDPSolver::compute_costs )

/*--------------------------------------------------------------------------*/

// Add the convex piecewise-linear reserve discount g to the just-built pieces
// of z_{h,k}, which occupy m[ begm .. mcnt ) and coeffs[ begt .. coeffcnt ) at
// the top of the current ping-pong half. The base function has ( vcnt + 1 )
// pieces; this refines the breakpoint set with g's internal breakpoints, adds
// g's (slope, intercept) on each resulting sub-interval, writes the augmented
// pieces back in place (expanding upward; the half was over-allocated for it),
// updates mcnt / coeffcnt / vcnt and returns the unconstrained minimizer of the
// (still convex) augmented function \hat z = z + g.

double ThermalUnitDPSolver::DPEDSolver::augment_with_g(
 const std::vector< g_piece > & g , Index begm , Index begt ,
 Index & mcnt , Index & coeffcnt , int & vcnt )
{
 const Index npc = Index( vcnt ) + 1;        // number of base pieces
 const double lo0 = m[ begm ] , hi0 = m[ begm + npc ];

 // snapshot the base pieces (they get overwritten in place below)
 std::vector< coeff_t > base( coeffs.begin() + begt , coeffs.begin() + begt + npc );
 std::vector< double >  bm( m.begin() + begm , m.begin() + begm + npc + 1 );

 // merged breakpoint set: base endpoints + g's internal breakpoints
 std::vector< double > bk( bm );
 for( const auto & gp : g ) {
  if( ( gp.lo > lo0 + 1e-12 ) && ( gp.lo < hi0 - 1e-12 ) ) bk.push_back( gp.lo );
  if( ( gp.hi > lo0 + 1e-12 ) && ( gp.hi < hi0 - 1e-12 ) ) bk.push_back( gp.hi );
  }
 std::sort( bk.begin() , bk.end() );
 bk.erase( std::unique( bk.begin() , bk.end() ,
            []( double a , double b ){ return( b - a <= 1e-12 ); } ) , bk.end() );

 auto gval = [ & ]( double x , double & slope , double & icpt ) {
  slope = 0; icpt = 0;
  for( const auto & gp : g )
   if( ( x >= gp.lo - 1e-12 ) && ( x <= gp.hi + 1e-12 ) ) {
    slope = gp.slope; icpt = gp.intercept; return;
    }
  };
 auto basei = [ & ]( double x ) -> Index {
  Index j = 0;
  while( ( j + 1 < npc ) && ( x >= bm[ j + 1 ] ) ) ++j;
  return( j );
  };

 const Index nc = Index( bk.size() ) - 1;    // number of augmented pieces
 m[ begm ] = bk[ 0 ];
 for( Index i = 0 ; i < nc ; ++i ) {
  const double mid = 0.5 * ( bk[ i ] + bk[ i + 1 ] );
  const Index bj = basei( mid );
  double gs , gi; gval( mid , gs , gi );
  coeffs[ begt + i ].alfa  = base[ bj ].alfa;
  coeffs[ begt + i ].beta  = base[ bj ].beta + gs;
  coeffs[ begt + i ].gamma = base[ bj ].gamma + gi;
  m[ begm + i + 1 ] = bk[ i + 1 ];
  }
 mcnt = begm + nc + 1;
 coeffcnt = begt + nc;
 vcnt = int( nc ) - 1;

 // unconstrained minimizer of the convex augmented function
 double bestv = Inf< double >() , bestp = bk[ 0 ];
 for( Index i = 0 ; i < nc ; ++i ) {
  const auto & cc = coeffs[ begt + i ];
  const double lo = bk[ i ] , hi = bk[ i + 1 ];
  double ps;
  if( std::abs( cc.alfa ) <= 1e-16 )
   ps = ( cc.beta <= 0 ? hi : lo );
  else
   ps = std::min( hi , std::max( lo , - cc.beta / ( 2 * cc.alfa ) ) );
  const double val = cc.alfa * ps * ps + cc.beta * ps + cc.gamma;
  if( val < bestv ) { bestv = val; bestp = ps; }
  }
 return( bestp );

 }  // end( ThermalUnitDPSolver::DPEDSolver::augment_with_g )

/*--------------------------------------------------------------------------*/

void ThermalUnitDPSolver::DPEDSolver::compute_power_variables( Index k ,
                                                               std::vector< double > & p )
{
 auto & delta_ramp_up = f_solver->delta_ramp_up;
 auto & delta_ramp_down = f_solver->delta_ramp_down;

 p[ k ] = con_p[ k ];

 for( Index t = k ; t-- > f_h ; ) {
  /* Project unconstrained optimal value unc_p[ t ] on the interval:
   * [ p[ t + 1 ] - delta_ramp_up[ t ] , p[ t + 1 ] + delta_ramp_down[ t ] ]
   *
   * If the unconstrained optimal value is on the left of the interval,
   * then the optimal power value is the left endpoint of the function.
   *
   * If the unconstrained optimal value is inside the interval,
   * then the optimal power value is exactly the unconstrained optimal value.
   *
   * If the unconstrained optimal value is on the right of the interval,
   * then the optimal power value is the right endpoint of the function.
   */

  if( unc_p[ t ] < p[ t + 1 ] - delta_ramp_up[ t ] )
   p[ t ] = p[ t + 1 ] - delta_ramp_up[ t ];
  else
   if( unc_p[ t ] <= p[ t + 1 ] + delta_ramp_down[ t ] )
    p[ t ] = unc_p[ t ];
   else
    p[ t ] = p[ t + 1 ] + delta_ramp_down[ t ];
  }
 }  // end( ThermalUnitDPSolver::compute_power_variables )

/*--------------------------------------------------------------------------*/
/*----------------- End File ThermalUnitDPSolver.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
