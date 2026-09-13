/*--------------------------------------------------------------------------*/
/*-------------------- File ThermalUnitDPSolverBase.cpp --------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of ThermalUnitDPSolverBase: the data loading, the convex
 * piecewise-quadratic value-function machinery and the spinning-reserve
 * model shared by ThermalUnitDPSolver and ThermalUnitExtDPSolver. */
/*--------------------------------------------------------------------------*/

// Diagnostic instrumentation shared with ThermalUnitExtDPSolver.cpp; keep 0
// in production. Defined BEFORE the header include so the profiling counters
// at the end of ThermalUnitDPSolverBase.h (under #if TUEDPS_PROFILE) get a
// single shared definition across the two translation units.
#define TUEDPS_PROFILE 0

#include "ThermalUnitDPSolverBase.h"

#include "ThermalUnitBlock.h"

#include <algorithm>
#include <cmath>
#include <functional>

#if TUEDPS_PROFILE
 #include <chrono>
 #include <iostream>
#endif

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*---------------------- STATIC THREAD-LOCAL SCRATCH -----------------------*/
/*--------------------------------------------------------------------------*/
// Per-thread reusable storage for the PQ machinery: static thread_local so
// the parallel per-source sweeps of ThermalUnitDPSolver (parDP) never race
// on it, while the single-threaded run-length solver is unaffected.

thread_local std::vector< ThermalUnitDPSolverBase::PQFun >
 ThermalUnitDPSolverBase::m_pqpool;
thread_local ThermalUnitDPSolverBase::PQFun
 ThermalUnitDPSolverBase::m_raw;

/*--------------------------------------------------------------------------*/
/*------------------------- SHARED DATA LOADING ----------------------------*/
/*--------------------------------------------------------------------------*/

void ThermalUnitDPSolverBase::load_common_parameters( void )
{
 bool owned = f_Block->is_owned_by( f_id );
 if( ( ! owned ) && ( ! f_Block->read_lock() ) )
  throw( std::runtime_error(
  "ThermalUnitDPSolverBase::load_common_parameters: unable to lock the Block."
   ) );

 auto b = static_cast< ThermalUnitBlock * >( f_Block );

 // scalar parameters of the unit
 time_horizon = b->get_time_horizon();
 init_up_down_time = b->get_init_up_down_time();
 min_up_time = b->get_min_up_time();
 min_down_time = b->get_min_down_time();
 initial_power = b->get_initial_power();

 // t_init: first instant at which the commitment decision is genuinely
 // free. If the unit was on before the horizon and min_up_time has not
 // been satisfied yet, the first few instants must stay on; symmetrically
 // for init_up_down_time < 0 and min_down_time. This mirrors the same
 // bookkeeping of ThermalUnitDPSolver and is used only to correctly
 // compute the start_up / shut_down indicator variables in get_var_solution.
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

 // startup costs, power bounds and shutdown/startup ramp bounds
 startup_costs = b->get_start_up_cost();
 min_power = b->get_min_power();
 max_power = b->get_max_power();
 bound_on = b->get_start_up_limit();
 bound_down = b->get_shut_down_limit();

 // ramp-up/down limits default to max_power (no effective ramping) when
 // the Block does not set them explicitly
 has_ramp_up = ! b->get_delta_ramp_up().empty();
 if( ! has_ramp_up )
  delta_ramp_up = max_power;
 else
  delta_ramp_up = b->get_delta_ramp_up();

 has_ramp_down = ! b->get_delta_ramp_down().empty();
 if( ! has_ramp_down )
  delta_ramp_down = max_power;
 else
  delta_ramp_down = b->get_delta_ramp_down();

 // cost coefficients
 // f_t(p) = quad_term[t] p^2 + linear_term[t] p + const_term[t]
 retrieve_term( quad_term   , b->get_quad_term() );
 retrieve_term( linear_term , b->get_linear_term() );
 retrieve_term( const_term  , b->get_const_term() );

 // spinning reserve: participation factors (caps in pr <= rho_p p and
 // sr <= rho_s p) and objective cost coefficients. The cost defaults to the
 // participation factor but may carry a (possibly negative) Lagrangian
 // price, so it is read from the separate cost getter. All empty when the
 // reserve is absent.
 primary_rho = b->get_primary_rho();
 secondary_rho = b->get_secondary_rho();
 primary_reserve_cost = b->get_primary_spinning_reserve_cost();
 secondary_reserve_cost = b->get_secondary_spinning_reserve_cost();
 if( primary_reserve_cost.empty() )
  primary_reserve_cost = primary_rho;
 if( secondary_reserve_cost.empty() )
  secondary_reserve_cost = secondary_rho;

 // reactive power (AC instances): read the box [Qmin,Qmax] and the
 // (dualized) linear cost coefficient on q[t]. q[t] is separable from the
 // DP, so run_DP() adds its optimal contribution as a constant. Empty
 // reactive_linear_term -> the unit has no reactive power and the term is
 // skipped.
 if( b->get_reactive_power( 0 ) ) {
  retrieve_term( reactive_linear_term , b->get_reactive_linear_term() );
  reactive_min.resize( time_horizon );
  reactive_max.resize( time_horizon );
  reactive_min_on.resize( time_horizon );
  reactive_max_on.resize( time_horizon );
  bool any_on = false;
  for( Index t = 0 ; t < time_horizon ; ++t ) {
   reactive_min[ t ] = b->get_min_reactive_power( t );
   reactive_max[ t ] = b->get_max_reactive_power( t );
   reactive_min_on[ t ] = b->get_min_reactive_power_on( t );
   reactive_max_on[ t ] = b->get_max_reactive_power_on( t );
   any_on = any_on || ( reactive_min_on[ t ] != 0.0 ) ||
            ( reactive_max_on[ t ] != 0.0 );
   }
  if( ! any_on ) {  // no commitment gating: keep the plain separable box
   reactive_min_on.clear();
   reactive_max_on.clear();
   }
  }
 else {
  reactive_linear_term.clear();
  reactive_min_on.clear();
  reactive_max_on.clear();
  }

 // design (investment): present iff the unit carries a nonzero investment
 // cost. The DP solves the operational problem assuming the unit exists;
 // run_DP() then decides whether to build it. See run_DP() for the
 // threshold rule.
 // NOTE: the actual building cost is the *current* coefficient of the
 // design variable in the Objective (get_design_cost()), which may differ
 // from the structural investment cost if a dualizing Solver has changed
 // it; this is kept in sync via the eSetInvCost Modification.
 has_design = ( b->get_investment_cost() != 0 );
 design_cost = b->get_design_cost();
 design_on = false;

 if( ! owned )
  f_Block->read_unlock();

}  // end( ThermalUnitDPSolverBase::load_common_parameters )

/*--------------------------------------------------------------------------*/

double ThermalUnitDPSolverBase::reserve_alloc( Index t , double p ,
                                              double & pr , double & sr ,
                                              double cap ) const
{
 // symmetric reserve band beta_t(p) = min( p - min_power , cap - p ): the
 // reserve must fit both above (the cap constraint, p + pr + sr <= cap) and
 // below (min_power constraint, p - pr - sr >= min_power) the production p.
 // cap is max_power at an interior period and the tighter start-up/shut-down
 // cap at a boundary period.
 return( reserve_alloc_band( t , p ,
                             std::min( p - min_power[ t ] , cap - p ) ,
                             pr , sr ) );

 }  // end( ThermalUnitDPSolverBase::reserve_alloc )

/*--------------------------------------------------------------------------*/

Solution * ThermalUnitDPSolverBase::pack_Solution(
                                       const std::vector< double > & p ,
                                       const std::vector< double > & u ,
                                       const std::vector< double > & pr ,
                                       const std::vector< double > & sr ,
                                       const std::vector< double > & q ,
                                       double design ) const
{
 // the shape of the Solution is a decision of the Block, not of the Solver:
 // ask it for its own one [see UnitBlock::new_Solution()], which is how a
 // NuclearUnitBlock gets a NuclearUnitBlockSolution out of this. note that
 // new_Solution() only builds the object, it reads no Variable
 auto sol = static_cast< ThermalUnitBlockSolution * >(
		     static_cast< ThermalUnitBlock * >( f_Block
						        )->new_Solution() );

 // a ThermalUnitBlock is one generator over the time horizon
 sol->set_dimensions( 1 , time_horizon );

 auto sz = boost::multi_array< double , 2 >::extent_gen()[ 1 ][ time_horizon ];

 auto pack = [ & ]( const std::vector< double > & v ) {
  boost::multi_array< double , 2 > a( sz );
  for( Index i = 0 ; i < time_horizon ; ++i )
   a[ 0 ][ i ] = v[ i ];
  return( a );
  };

 sol->set_active_power( pack( p ) );
 sol->set_commitment( pack( u ) );

 // the Objective pays the start-up through its own Variable, so the
 // indicators travel with the Solution rather than being derived from the
 // commitment when it is written back: deriving them works for one schedule
 // and not for a convex combination of several, where the start-ups of the
 // averaged commitment are fewer than the average of the start-ups
 {
  std::vector< double > su , sd;
  static_cast< ThermalUnitBlock * >( f_Block )->derive_start_up( u , su , sd );
  if( ! su.empty() ) {
   sol->set_start_up( std::move( su ) );
   sol->set_shut_down( std::move( sd ) );
   }
  }

 if( ! pr.empty() )
  sol->set_primary_spinning_reserve( pack( pr ) );

 if( ! sr.empty() )
  sol->set_secondary_spinning_reserve( pack( sr ) );

 if( ! q.empty() )
  sol->set_reactive_power( pack( q ) );

 sol->set_design( design );

 return( sol );

 }  // end( ThermalUnitDPSolverBase::pack_Solution )

/*--------------------------------------------------------------------------*/

bool ThermalUnitDPSolverBase::has_reactive_power( void ) const
{
 auto b = static_cast< ThermalUnitBlock * >( f_Block );

 for( Index i = 0 ; i < time_horizon ; ++i )
  if( b->get_min_reactive_power( i ) || b->get_max_reactive_power( i ) ||
      b->get_min_reactive_power_on( i ) || b->get_max_reactive_power_on( i ) )
   return( true );

 return( false );

 }  // end( ThermalUnitDPSolverBase::has_reactive_power )

/*--------------------------------------------------------------------------*/

double ThermalUnitDPSolverBase::reserve_alloc_band( Index t , double p ,
                                                   double H , double & pr ,
                                                   double & sr ) const
{
 pr = sr = 0;
 if( H <= 0 )
  return( 0 );

 const double cap_p = primary_rho.empty()   ? 0 : primary_rho[ t ]   * p;
 const double cap_s = secondary_rho.empty() ? 0 : secondary_rho[ t ] * p;
 const double cp = primary_reserve_cost.empty()   ? 0
                                                  : primary_reserve_cost[ t ];
 const double cs = secondary_reserve_cost.empty()
                   ? 0 : secondary_reserve_cost[ t ];

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

 }  // end( ThermalUnitDPSolverBase::reserve_alloc_band )

/*--------------------------------------------------------------------------*/

ThermalUnitDPSolverBase::PQFun
ThermalUnitDPSolverBase::build_reserve_discount( Index t , double cap ) const
{
 PQFun G;
 const double cp = primary_reserve_cost.empty()   ? 0
                                                  : primary_reserve_cost[ t ];
 const double cs = secondary_reserve_cost.empty()
                   ? 0 : secondary_reserve_cost[ t ];
 if( ( cp >= 0 ) && ( cs >= 0 ) )
  return( G );  // no negative price: g_t == 0

 const double lo = min_power[ t ];
 const double hi = cap;  // upper power cap U_t (interior up, or
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
  G.push_back( { 0.0 , slope , ga - slope * a , a , b } );
  }
 return( G );

 }  // end( ThermalUnitDPSolverBase::build_reserve_discount )

/*--------------------------------------------------------------------------*/

void ThermalUnitDPSolverBase::add_pwq( PQFun & F , const PQFun & G )
{
 if( G.empty() || F.empty() )
  return;

 const double tol = 1e-12;
 PQFun out;
 out.reserve( F.size() + G.size() );

 std::size_t gi = 0;
 for( const auto & fp : F ) {
  double p = fp.left;
  while( ( gi < G.size() ) && ( G[ gi ].right <= p + tol ) )
   ++gi;
  if( fp.right <= fp.left + tol ) {  // degenerate (zero-width) piece: keep
   double gb = 0 , gc = 0;           // it, adding G at the single point p
   if( ( gi < G.size() ) && ( G[ gi ].left <= p + tol ) ) {
    gb = G[ gi ].beta;
    gc = G[ gi ].gamma;
    }
   out.push_back( { fp.alfa , fp.beta + gb , fp.gamma + gc ,
                    fp.left , fp.right } );
   continue;
   }
  while( p < fp.right - tol ) {
   double r = fp.right , gb = 0 , gc = 0;
   if( gi < G.size() ) {
    if( G[ gi ].left <= p + tol ) {       // G covers p
     r = std::min( fp.right , G[ gi ].right );
     gb = G[ gi ].beta;
     gc = G[ gi ].gamma;
     }
    else                                  // gap: G == 0 up to G[gi].left
     r = std::min( fp.right , G[ gi ].left );
    }
   out.push_back( { fp.alfa , fp.beta + gb , fp.gamma + gc , p , r } );
   p = r;
   if( ( gi < G.size() ) && ( G[ gi ].right <= p + tol ) )
    ++gi;
   }
  }

 // coalesce adjacent pieces the overlay left as the SAME quadratic: G is
 // flat or zero over long stretches, so F+G repeats its coefficients across
 // many of the breakpoints just introduced. This is exact (equal
 // coefficients => equal function), and without it the piece count grows
 // every on-step along a run, reaching a few hundred and making every
 // downstream operation (and especially the residual-ramp Gmin, which scans
 // the pieces) proportionally slower.
 if( out.size() > 1 ) {
  std::size_t w = 0;
  for( std::size_t i = 1 ; i < out.size() ; ++i ) {
   PieceQuad & b = out[ w ];
   const PieceQuad & c = out[ i ];
   if( ( std::abs( b.alfa  - c.alfa  ) <=
         1e-11 * ( 1 + std::abs( b.alfa  ) ) ) &&
       ( std::abs( b.beta  - c.beta  ) <=
         1e-11 * ( 1 + std::abs( b.beta  ) ) ) &&
       ( std::abs( b.gamma - c.gamma ) <=
         1e-9  * ( 1 + std::abs( b.gamma ) ) ) )
    b.right = c.right;             // same quadratic: absorb c into b
   else
    out[ ++w ] = c;
   }
  out.resize( w + 1 );
  }
 F = std::move( out );

 }  // end( ThermalUnitDPSolverBase::add_pwq )

/*--------------------------------------------------------------------------*/

void ThermalUnitDPSolverBase::retrieve_term(
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

double ThermalUnitDPSolverBase::eval( const PQFun & F , double p )
{
 if( F.empty() )
  return( TUEDPINF );
 for( const auto & pc : F )
  if( pc.left - 1e-12 <= p && p <= pc.right + 1e-12 )
   return( eval_piece( pc , p ) );
 return( TUEDPINF );
 }

/*--------------------------------------------------------------------------*/

double ThermalUnitDPSolverBase::argmin_piece( const PieceQuad & pc )
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

std::pair< double , double > ThermalUnitDPSolverBase::min_over(
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

void ThermalUnitDPSolverBase::shift_by( PQFun & F , double c )
{
 for( auto & pc : F )
  pc.gamma += c;
 }

/*--------------------------------------------------------------------------*/

void ThermalUnitDPSolverBase::add_quadratic( PQFun & F ,
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

void ThermalUnitDPSolverBase::clamp_domain( PQFun & F ,
                                            double lo , double hi )
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

bool ThermalUnitDPSolverBase::is_dominated_by( const PQFun & F1 ,
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

  // zero-width (single-point) F1 piece, e.g. a restart whose start-up limit
  // equals the minimum power so its domain is [P, P]. The general sweep
  // below skips it (its guard p < r1 - tol is false), which would wrongly
  // report vacuous domination. Handle it explicitly: F1 is dominated here
  // iff F2 is defined at the point and lies at or below F1 there.
  if( r1 - l1 <= tol ) {
   const double p0 = l1;
   bool covered = false;
   for( const auto & q : F2 )
    if( ( q.left <= p0 + tol ) && ( p0 <= q.right + tol ) ) {
     const double d = ( F1[ i ].alfa  - q.alfa  ) * p0 * p0
                    + ( F1[ i ].beta  - q.beta  ) * p0
                    + ( F1[ i ].gamma - q.gamma );
     if( d < - eps )
      return( false );  // F1 strictly below F2 at the point
     covered = true;
     break;
     }
   if( ! covered )
    return( false );    // F2 is +INF at the point: not dominated
   continue;
   }

  double p = l1;
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

void ThermalUnitDPSolverBase::sliding_min(
 const PQFun & F , double ramp_up , double ramp_down ,
 double lo , double hi , PQFun & out )
{
 out.clear();  // keeps capacity
 if( F.empty() || lo > hi + 1e-12 )
  return;

 // a window [ -ramp_down , ramp_up ] of the move p_t - q that does not
 // contain 0, i.e., a move that must be strictly upwards (ramp_down < 0) or
 // strictly downwards (ramp_up < 0), is a one-sided window translated by its
 // endpoint nearest to 0: with d that endpoint,
 //   min{ F( q ) : p_t - q in [ a , b ] } = H( p_t - d ) ,
 //   H( x ) = min{ F( q ) : x - q in [ a - d , b - d ] }
 // where [ a - d , b - d ] contains 0; the degenerate window a = b is the
 // exact shift p_t = q + a (a move at the full ramp rate). An empty window
 // ( ramp_up + ramp_down < 0 ) admits no move at all
 if( ( ramp_up < 0 ) || ( ramp_down < 0 ) ) {
  if( ramp_up + ramp_down < -1e-12 )
   return;
  const double d = ( ramp_down < 0 ) ? - ramp_down : ramp_up;
  PQFun H;
  if( ramp_down < 0 )
   sliding_min( F , ramp_up - d , 0.0 , lo - d , hi - d , H );
  else
   sliding_min( F , 0.0 , ramp_down + d , lo - d , hi - d , H );
  out.reserve( H.size() );
  for( const auto & pc : H )  // H( p_t - d ) as a function of p_t
   out.push_back( { pc.alfa , pc.beta - 2 * pc.alfa * d ,
                    pc.alfa * d * d - pc.beta * d + pc.gamma ,
                    pc.left + d , pc.right + d } );
  return;
  }

 // the null window is the identity: F itself, restricted to [ lo , hi ]
 // (the construction below would split the piece containing the minimiser);
 // a domain that is a single point stays a single point, rather than being
 // taken for the collapse of [ lo , hi ] handled at the end
 if( ( ramp_up == 0 ) && ( ramp_down == 0 ) ) {
  const double dl = std::max( F.front().left , lo );
  const double dr = std::min( F.back().right , hi );
  if( dr < dl - 1e-12 )
   return;                      // no common point: no feasible move
  if( dr - dl <= 1e-12 ) {       // a single common point
   const double v = eval( F , dl );
   if( v < TUEDPINF )
    out.push_back( { 0 , 0 , v , dl , dl } );
   return;
   }
  out.reserve( F.size() );
  for( const auto & pc : F ) {
   const double l = std::max( pc.left , lo );
   const double r = std::min( pc.right , hi );
   if( l < r - 1e-15 )
    out.push_back( { pc.alfa , pc.beta , pc.gamma , l , r } );
   }
  return;
  }

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

 // build the output in an intermediate buffer, sort + clamp at the end;
 // m_raw is a reused member scratch (cleared, capacity retained)
 PQFun & raw = m_raw;
 raw.clear();
 raw.reserve( F.size() * 2 + 1 );

 // (1) left-shifted pieces: portions of F with q <= p_star.
 // The substitution q = p_t + ramp_down implies that the interval
 // [pc.left, pc.right] of the original variable q maps to
 // [pc.left - ramp_down, pc.right - ramp_down] of the new variable p_t,
 // and the coefficients transform as
 //   F(q) = alfa q^2 + beta q + gamma
 //        = alfa (p_t + rd)^2 + beta (p_t + rd) + gamma
 //        = alfa p_t^2 + (2 alfa rd + beta) p_t
 //          + (alfa rd^2 + beta rd + gamma)
 for( const auto & pc : F ) {
  double a = pc.left , b = pc.right;
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
 //        = alfa p_t^2 + (-2 alfa ru + beta) p_t
 //          + (alfa ru^2 - beta ru + gamma)
 for( const auto & pc : F ) {
  double a = pc.left , b = pc.right;
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
 out.reserve( raw.size() );
 for( const auto & pc : raw ) {
  double l = std::max( pc.left  , lo );
  double r = std::min( pc.right , hi );
  if( l < r - 1e-15 )
   out.push_back( { pc.alfa , pc.beta , pc.gamma , l , r } );
  }

 // a domain that collapses to a single point [min_power == max_power at t,
 // say a unit that at t can only sit at one power] leaves no piece of
 // positive width, and an empty output means to the caller that there is no
 // feasible transition at all. The value function there is a point: emit it
 // as a zero-width piece, which the rest of the machinery handles [see
 // add_pwq() and min_over()]. If no piece of the transformed function
 // covers that point the transition really is infeasible, and the output
 // stays empty
 // (if [ lo , hi ] is not a single point the transformed function only
 // touches it at one of its ends, and the value is at that point alone)
 if( out.empty() && ( lo <= hi + 1e-12 ) )
  for( const double x : { lo , hi } ) {
   double v = TUEDPINF;
   for( const auto & pc : raw )
    if( ( pc.left <= x + 1e-12 ) && ( x <= pc.right + 1e-12 ) )
     v = std::min( v , eval_piece( pc , x ) );
   if( v < TUEDPINF ) {
    if( hi - lo <= 1e-12 )
     out.push_back( { 0 , 0 , v , lo , hi } );
    else
     out.push_back( { 0 , 0 , v , x , x } );
    break;
    }
   }

 }  // end( ThermalUnitDPSolverBase::sliding_min )

/*--------------------------------------------------------------------------*/

double ThermalUnitDPSolverBase::reserve_reward( Index t , double p ,
                                               double H ) const
{
 if( H <= 0 )
  return( 0 );
 const double cap_p = primary_rho.empty()   ? 0 : primary_rho[ t ]   * p;
 const double cap_s = secondary_rho.empty() ? 0 : secondary_rho[ t ] * p;
 const double cp = primary_reserve_cost.empty()   ? 0
                                                  : primary_reserve_cost[ t ];
 const double cs = secondary_reserve_cost.empty()
                   ? 0 : secondary_reserve_cost[ t ];
 double pr = 0 , sr = 0 , rem = H;
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

 }  // end( ThermalUnitDPSolverBase::reserve_reward )

/*--------------------------------------------------------------------------*/

void ThermalUnitDPSolverBase::sliding_min_corr(
 const PQFun & F , double ramp_up , double ramp_down ,
 double lo , double hi , Index t , PQFun & out , double acap ,
 double win_up , double win_down )
{
 out.clear();
 // window of the scheduled move q in [ p - wu , p + wd ], the tent of the
 // reserve deliverability being that of ramp_up / ramp_down
 const double wu = std::isnan( win_up ) ? ramp_up : win_up;
 const double wd = std::isnan( win_down ) ? ramp_down : win_down;
 const double cp = primary_reserve_cost.empty()   ? 0
                                                  : primary_reserve_cost[ t ];
 const double cs = secondary_reserve_cost.empty()
                   ? 0 : secondary_reserve_cost[ t ];
 // no reserve rewarded at t (no negative price) => corr == 0 => exact
 // standard sliding minimum (energy-only stays bit-identical to sliding_min)
 if( ( cp >= 0 ) && ( cs >= 0 ) ) {
  sliding_min( F , wu , wd , lo , hi , out );
  return;
  }
#if TUEDPS_PROFILE
 if( std::getenv( "TUEDPS_NOCORR" ) ) {  // diagnostic: capacity-band model
  sliding_min( F , wu , wd , lo , hi , out );
  return;                               // (g0 still added by add_pwq outside)
  }
#endif
 if( F.empty() || ( lo > hi + 1e-12 ) )
  return;
#if TUEDPS_PROFILE
 ++g_smc; g_Fsum += F.size(); if( F.size() > g_Fmax ) g_Fmax = F.size();
#endif

 // pmax is the CAPACITY-band ceiling of the reserve (A=min(p-pmin,pmax-p)):
 // max_power at an interior step, the shut-down cap when the run closes here.
 const double pmin = min_power[ t ];
 const double pmax = ( acap < 0.0 ) ? max_power[ t ] : acap;
 const double domL = F.front().left , domR = F.back().right;

 // ordered participation factors (most-negative-price reserve first), used to
 // place corr's reward-regime kinks; mirrors build_reserve_discount()
 const double rp = primary_rho.empty()   ? 0 : primary_rho[ t ];
 const double rs = secondary_rho.empty() ? 0 : secondary_rho[ t ];
 double rho1 = 0 , rho2 = 0;
 if( ( cp < 0 ) && ( cs < 0 ) ) {
  if( cp <= cs ) { rho1 = rp; rho2 = rs; } else { rho1 = rs; rho2 = rp; }
  }
 else if( cp < 0 ) rho1 = rp; else rho1 = rs;
 const double rhot = rho1 + rho2;
#if TUEDPS_PROFILE
 bool corr_hit = false;   // did corr>0 ever land at the optimum in this call?
#endif

 auto Aband = [ & ]( double p ) {
  return( std::min( p - pmin , pmax - p ) ); };
 // corr_t(q,p) >= 0 : reward under the residual-ramp band minus that under
 // the capacity band; zero where the ramp tent B(p-q) >= the capacity band
 // A(p). FULL band reward g(p,min(A,B)) folded into the sliding minimum:
 // the capacity credit g^0=g(p,A) (added per-period BEFORE) and the ramp
 // penalty g(B)-g(A) are combined here as one term. This makes the
 // transition value EXACTLY convex (no separate +g^0/-g^0 add on mismatched
 // meshes to cancel). The caller does not add eff_disc after
 // sliding_min_corr. Zero in q where A<=0 (no band).
 auto corr = [ & ]( double q , double p ) -> double {
  const double A = Aband( p );
  if( A <= 0 ) return( 0 );
  const double d = p - q;
  double B = std::min( ramp_up - d , ramp_down + d );
  if( B < 0 ) B = 0;
  return( reserve_reward( t , p , std::min( A , B ) ) );
  };

 // F is a sorted, contiguous partition => the piece containing q is found
 // in O(log|F|); a linear scan here would make every Gmin O(|F|^2) (eval
 // called O(|F|) times, each an O(|F|) scan), the single biggest cost in
 // the residual-ramp construction.
 auto findF = [ & ]( double q ) -> int {
  int a = 0 , b = int( F.size() ) - 1;
  while( a < b ) {
   const int m = ( a + b + 1 ) / 2;
   if( F[ m ].left <= q + 1e-12 ) a = m; else b = m - 1;
   }
  return( a );
  };
 auto evalF = [ & ]( double q ) -> double {
  const PieceQuad & pc = F[ findF( q ) ];
  return( ( pc.alfa * q + pc.beta ) * q + pc.gamma );
  };

 // exact ( G(p) , argmin q* ) = min over the ramp window of F(q) + corr(q,p).
 // phi(q) = F(q)+corr(q,p) is convex; the minimiser is a breakpoint or a
 // combined (F-piece x corr-segment) vertex, all enumerated in qs.
 auto Gmin = [ & ]( double p ) -> std::pair< double , double > {
#if TUEDPS_PROFILE
  ++g_gmin;
#endif
  const double qL = std::max( p - wu , domL );
  const double qR = std::min( p + wd , domR );
  if( qL > qR + 1e-12 ) return std::make_pair( TUEDPINF , 0.0 );
  if( qR - qL <= 1e-12 )                 // degenerate window: single point
   return std::make_pair( evalF( qL ) + corr( qL , p ) , qL );
  const double A = Aband( p );
  static thread_local std::vector< double > qs;
  qs.clear();
  qs.push_back( qL ); qs.push_back( qR );
  // only the F-piece boundaries inside the ramp window matter; reach the
  // first via binary search and walk forward, instead of scanning all of F
  // (which can be hundreds of pieces) on every one of the ~280 Gmin calls
  // per corr build
  for( int k = findF( qL ) ; k < int( F.size() ) && F[ k ].left < qR - 1e-12 ;
       ++k )
   if( F[ k ].left > qL + 1e-12 )
    qs.push_back( F[ k ].left );
  auto add_v = [ & ]( double v ) {  // tent B(p-q)=v => q = p-ru+v or p+rd-v
   if( v <= 1e-12 ) return;
   const double q1 = p - ramp_up + v , q2 = p + ramp_down - v;
   if( ( q1 > qL + 1e-12 ) && ( q1 < qR - 1e-12 ) ) qs.push_back( q1 );
   if( ( q2 > qL + 1e-12 ) && ( q2 < qR - 1e-12 ) ) qs.push_back( q2 );
   };
  if( A > 0 ) { add_v( A ); add_v( rho1 * p );
                if( rho2 > 0 ) add_v( rhot * p ); }
  { const double qpk = p - 0.5 * ( ramp_up - ramp_down );  // tent peak
    if( ( qpk > qL + 1e-12 ) && ( qpk < qR - 1e-12 ) ) qs.push_back( qpk ); }
  std::sort( qs.begin() , qs.end() );
  double best = TUEDPINF , bestq = qL;
  // corr = g(p,min(A,B)) (full band reward, capacity+ramp folded). Evaluate
  // F(q)+corr at every distinct candidate ONCE (caching corr for the slope
  // test), instead of re-evaluating each shared interval endpoint twice.
  auto corrq = [ & ]( double q ) -> double {
   if( A <= 0 ) return( 0.0 );
   const double d = p - q;
   double B = std::min( ramp_up - d , ramp_down + d );
   if( B < 0 ) B = 0;
   return( reserve_reward( t , p , std::min( A , B ) ) );
   };
  auto updv = [ & ]( double q , double v ) {
   if( v < best ) { best = v; bestq = q; } };
  static thread_local std::vector< double > cq; cq.resize( qs.size() );
  static thread_local std::vector< int > kf; kf.resize( qs.size() );
  // qs is sorted and every in-window F-boundary is in qs, so no F-boundary
  // lies strictly inside (qs[i],qs[i+1]): walk a piece-cursor instead of
  // binary-searching findF at each candidate. kf[i] is the F-piece active
  // AT and just-right-of qs[i].
  int kc = findF( qs[ 0 ] );
  for( std::size_t i = 0 ; i < qs.size() ; ++i ) {
   while( ( kc + 1 < int( F.size() ) ) &&
          ( F[ kc + 1 ].left <= qs[ i ] + 1e-12 ) )
    ++kc;
   kf[ i ] = kc;
   cq[ i ] = corrq( qs[ i ] );
   if( ( i > 0 ) && ( qs[ i ] - qs[ i - 1 ] <= 1e-12 ) ) continue;  // dup
   const PieceQuad & pc = F[ kc ];
   updv( qs[ i ] ,
         ( pc.alfa * qs[ i ] + pc.beta ) * qs[ i ] + pc.gamma + cq[ i ] );
   }
  for( std::size_t i = 0 ; i + 1 < qs.size() ; ++i ) {
   const double a = qs[ i ] , b = qs[ i + 1 ];
   if( b - a <= 1e-12 ) continue;
   const PieceQuad & pcm = F[ kf[ i ] ];  // no F-boundary inside (a,b)
   const double aF = pcm.alfa , bF = pcm.beta;
   if( aF > 1e-16 ) {
    const double sc = ( cq[ i + 1 ] - cq[ i ] ) / ( b - a );  // corr slope
    const double qv = - ( bF + sc ) / ( 2.0 * aF );
    if( qv > a && qv < b )
     updv( qv , ( pcm.alfa * qv + pcm.beta ) * qv + pcm.gamma + corrq( qv ) );
    }
   }
#if TUEDPS_PROFILE
  if( corr( bestq , p ) > 1e-9 ) { ++g_gmin_bite; corr_hit = true; }
#endif
  return std::make_pair( best , bestq );
  };
 auto Geval = [ & ]( double p ) -> double { return( Gmin( p ).first ); };

#if TUEDPS_PROFILE
 // EXACTNESS ORACLE (diagnostic, opt-in via TUEDPS_PARAMDIFF): the genuinely
 // trustworthy value of G(p) = min_{q in [qL,qR]} F(q)+corr(q,p), computed
 // INDEPENDENTLY of Gmin's breakpoint algebra so it can validate Gmin itself.
 // phi(q)=F(q)+corr(q,p) is convex in q: F is convex, and corr is
 // (convex-decreasing reward) composed with (concave ramp tent min A), hence
 // convex, so a ternary search converges geometrically to the exact
 // minimum (interval width (2/3)^N * span). A coarse grid min is folded in
 // as a belt-and-suspenders guard against any numerical non-convexity.
 // Values good to ~1e-12.
 auto Gtrue = [ & ]( double p ) -> std::pair< double , double > {
  const double qL = std::max( p - wu , domL );
  const double qR = std::min( p + wd , domR );
  if( qL > qR + 1e-12 ) return std::make_pair( TUEDPINF , 0.0 );
  auto phi = [ & ]( double q ) { return( evalF( q ) + corr( q , p ) ); };
  if( qR - qL <= 1e-13 ) return std::make_pair( phi( qL ) , qL );
  // coarse grid bracket (guard) -----------------------------------------
  const int NG = 4000;
  double gbest = TUEDPINF , gbestq = qL;
  for( int i = 0 ; i <= NG ; ++i ) {
   const double q = qL + ( qR - qL ) * i / NG;
   const double v = phi( q );
   if( v < gbest ) { gbest = v; gbestq = q; }
   }
  // ternary refine on the full convex interval (exact to fp) --------------
  double a = qL , b = qR;
  for( int it = 0 ;
       it < 300 && ( b - a ) > 1e-13 * std::max( 1.0 , std::abs( b ) ) ;
       ++it ) {
   const double m1 = a + ( b - a ) / 3 , m2 = b - ( b - a ) / 3;
   if( phi( m1 ) < phi( m2 ) ) b = m2; else a = m1;
   }
  const double qt = 0.5 * ( a + b ) , vt = phi( qt );
  return( vt <= gbest ? std::make_pair( vt , qt )
                      : std::make_pair( gbest , gbestq ) );
  };
#endif

 // reconstruct out(p) over the finite-G domain. KEY STRUCTURE: the
 // minimiser q*(p) is either a window edge (q = p-wu or p+wd)
 // and then G is quadratic in p (F shifted, minus a pw-linear reserve term)
 // or interior and CONSTANT in p on a fixed configuration, and then
 // G(p) = F(q*)+corr(q*,p) is LINEAR in p. So each piece is fitted by exact
 // interpolation of a genuine line (2 points, interior) or parabola
 // (3 points, edge) on a config-constant sub-interval; a verify probe
 // subdivides where the configuration changes, and at the depth cap a
 // straight segment is emitted (never a spurious curvature).
 const double plo = std::max( lo , domL - wd );
 const double phi = std::min( hi , domR + wu );
 if( plo > phi + 1e-12 )
  return;

 // ====================================================================
 // VERIFY-FREE PARAMETRIC ENGINE (6 affine forms). Builds G(p) by sweeping
 // p: ONE classify probe per piece -> exact analytic event set (all linear
 // solves) -> closed-form piece (curvature anchored, b,c interpolated from
 // 2 endpoints). Legacy engine, gated behind TUEDPS_PARAM; the production
 // default is the mpQP sweep below, the seed+verify+subdivide oracle is
 // gated behind TUEDPS_ORACLE.
 // Forms: PIN (q* const, G linear) + MOVING q*=sl*p+ic (G quad, curv
 // aF*sl^2): WUP(1,-ru) WDN(1,+rd) TENT(1,-dpk) BAND(2,ic)
 // KBAND(1-+cum_m, rd|-ru). See TUDPS.tex 5.7 for the taxonomy.
 const double dpk = 0.5 * ( ramp_up - ramp_down );
 const double Bmax = 0.5 * ( ramp_up + ramp_down );
 const double midA = 0.5 * ( pmin + pmax );
 double segcum[ 2 ] , segc[ 2 ]; int K = 0;
 if( ( cp < 0 ) && ( cs < 0 ) ) {
  if( cp <= cs ) { segc[ 0 ] = cp; segcum[ 0 ] = rp;
                   segc[ 1 ] = cs; segcum[ 1 ] = rp + rs; }
  else           { segc[ 0 ] = cs; segcum[ 0 ] = rs;
                   segc[ 1 ] = cp; segcum[ 1 ] = rs + rp; }
  K = 2; }
 else if( cp < 0 ) { segc[ 0 ] = cp; segcum[ 0 ] = rp; K = 1; }
 else if( cs < 0 ) { segc[ 0 ] = cs; segcum[ 0 ] = rs; K = 1; }
 auto nearestFb = [ & ]( double q ) -> double {  // nearest F breakpoint to q
  double best = F[ 0 ].left , bd = std::abs( q - F[ 0 ].left );
  for( const auto & pc : F ) { const double d2 = std::abs( q - pc.right );
    if( d2 < bd ) { bd = d2; best = pc.right; } }
  return( best ); };
 struct Reg { int form; double qbar , sl , ic; int m; };
 auto qstarR = [ & ]( const Reg & r , double p ) -> double {
  return( r.form == 0 ? r.qbar : r.sl * p + r.ic ); };
 auto GformR = [ & ]( const Reg & r , double p ) -> double {
  const double q = qstarR( r , p ); return( evalF( q ) + corr( q , p ) ); };
 auto curvR = [ & ]( const Reg & r , double p ) -> double {
  if( r.form == 0 ) return( 0.0 );
  return( F[ findF( qstarR( r , p ) ) ].alfa * r.sl * r.sl ); };
 // classify the regime governing (p, .] with ONE Gmin probe. WHICH bound
 // binds gives the form; the window-vs-domain-edge ambiguity is resolved by
 // the CONSTRAINT test (is the ramp window inside the F domain?), reliable
 // because the probe offset eps=1e-5|p| exceeds the position tol=1e-6|p|
 // (a smaller offset would fail just past a transition). The interior forms
 // are told apart by directly testing which convex kink q* rides, tent
 // peak / B=A / B=cum_m p, so no motion probe is needed (a stationary q*
 // landing on such a locus is a single-point coincidence, cut immediately
 // by the event mesh).
 auto classify = [ & ]( double p ) -> Reg {
  Reg r{ -1 , 0 , 0 , 0 , 0 };
  const double pe = p + 1e-5 * std::max( 1.0 , std::abs( p ) );
  auto g = Gmin( pe );
  if( g.first >= TUEDPINF ) return( r );
  const double q = g.second , tol = 1e-6 * std::max( 1.0 , std::abs( pe ) );
  const double qLb = std::max( pe - wu , domL );
  const double qRb = std::min( pe + wd , domR );
  if( std::abs( q - qRb ) < tol ) {          // upper bound binds (1 probe)
   if( pe + wd < domR - tol )
    { r.form = 2; r.sl = 1; r.ic = wd; }               // WDN
   else { r.form = 0; r.qbar = domR; }                 // PIN domR
   return( r ); }
  if( std::abs( q - qLb ) < tol ) {          // lower bound binds (1 probe)
   if( pe - wu > domL + tol )
    { r.form = 1; r.sl = 1; r.ic = -wu; }              // WUP
   else { r.form = 0; r.qbar = domL; }                 // PIN domL
   return( r ); }
  // interior minimiser: a stationary q* can coincide with a moving-locus
  // kink at an isolated p, so here (only here) a second probe distinguishes
  // a genuine moving form (q* tracks the kink) from a pinned stationary
  // point (q* const).
  const double e1 = 1e-7 * std::max( 1.0 , std::abs( p ) );
  const double slope = ( q - Gmin( p + e1 ).second ) / ( pe - ( p + e1 ) );
  const bool moving = std::abs( slope ) > 0.5;
  const double d = pe - q; double B = std::min( ramp_up - d , ramp_down + d );
  if( B < 0 ) B = 0; const double A = std::min( pe - pmin , pmax - pe );
  const bool upper = d > dpk;
  if( moving && ( std::abs( q - ( pe - dpk ) ) < tol ) )
   { r.form = 3; r.sl = 1; r.ic = -dpk; }          // TENT
  else if( moving &&
           ( std::abs( B - A ) < 1e-6 * std::max( 1.0 , std::abs( A ) ) ) ) {
   r.form = 4; r.sl = 2; r.m = upper ? 1 : 0;
   r.ic = upper ? ( -ramp_up - pmin ) : ( ramp_down - pmax ); }   // BAND
  else {
   int km = -1;
   for( int m = 0 ; m < K ; ++m )
    if( moving && ( std::abs( B - segcum[ m ] * pe ) <
                    1e-6 * std::max( 1.0 , segcum[ m ] * pe ) ) )
     { km = m; break; }
   if( km >= 0 ) { r.form = 5; r.m = km;   // KBAND: qc rides B=cum_m*p
    if( upper ) { r.sl = 1 + segcum[ km ]; r.ic = -ramp_up; }
    else        { r.sl = 1 - segcum[ km ]; r.ic = ramp_down; } }
   else { r.form = 0;                      // PIN interior stationary
    const double qb = nearestFb( q );
    r.qbar = ( std::abs( q - qb ) < 1e-6 ) ? qb : q; }
   }
  return( r );
  };
 // directional slope of corr in q at (q,p): corr'(q)=segc[seg(B)]*dB/dq,
 // 0 if B<=0 or B>=A. dB/dq=+1 on the upper tent branch (ru-d<rd+d), -1 on
 // the lower. Used by the pin-at-kink event and the verify-based optimality
 // certificate (formOK).
 auto corrDer = [ & ]( double q , double pe ) -> double {
  const double d = pe - q , Bup = ramp_up - d , Blo = ramp_down + d;
  double B = std::min( Bup , Blo ); if( B < 0 ) B = 0;
  const double A = std::min( pe - pmin , pmax - pe );
  if( ( B <= 1e-12 ) || ( B >= A - 1e-12 ) ) return( 0.0 );
  double marg = 0;
  if( ( K > 0 ) && ( B < segcum[ K - 1 ] * pe - 1e-12 ) ) {
   int m = 0; for( ; m < K ; ++m ) if( B <= segcum[ m ] * pe + 1e-9 ) break;
   if( m >= K ) m = K - 1; marg = segc[ m ]; }
  return( marg * ( ( Bup < Blo ) ? 1.0 : -1.0 ) );
  };
 // complete event set (> pc). Every entry a linear/trivial solve;
 // over-generation is safe (identical halves coalesce), MISSING one is not.
 // Includes the form-transition (stat-entry / window-reach / kink-meet)
 // events per form.
 auto next_event = [ & ]( const Reg & r , double pc ) -> double {
  double best = phi;
  auto rel = [ & ]( double pe ) {
   if( ( pe > pc + 1e-9 ) && ( pe < best ) ) best = pe; };
  rel( midA );                              // A(p) kink (all forms)
  for( int m = 0 ; m < K ; ++m ) {          // A = K_m (rising / falling)
   if( std::abs( 1 - segcum[ m ] ) > 1e-14 )
    rel( pmin / ( 1 - segcum[ m ] ) );
   rel( pmax / ( 1 + segcum[ m ] ) ); }
  if( r.form == 0 ) {                       // PIN q*=qbar
   const double qbar = r.qbar;
   rel( qbar + dpk );                       // tent peak
   rel( 0.5 * ( ramp_up + qbar + pmin ) );  // B=A rising-A / upper-tent
   rel( 0.5 * ( pmax - ramp_down + qbar ) );  // B=A falling-A / lower-tent
   for( int m = 0 ; m < K ; ++m ) {         // B = K_m (lower / upper tent)
    if( std::abs( 1 - segcum[ m ] ) > 1e-14 )
     rel( ( qbar - ramp_down ) / ( 1 - segcum[ m ] ) );
    rel( ( ramp_up + qbar ) / ( 1 + segcum[ m ] ) ); }
   rel( qbar + wu ); rel( qbar - wd );              // window-reach
   }                                                // (-> WUP / WDN)
  else if( ( r.form == 1 ) || ( r.form == 2 ) ) {  // WUP / WDN, q = p - sh
   const double sh = ( r.form == 1 ) ? wu : -wd;
   for( const auto & pc2 : F )
    { rel( pc2.left + sh ); rel( pc2.right + sh ); }        // F-cross
   const int i = findF( pc - sh );
   const double a = F[ i ].alfa , b = F[ i ].beta;
   if( a > 1e-16 )                          // stationary point enters piece
    rel( ( r.form == 1 ) ? ( wu - ( b + segc[ 0 ] ) / ( 2 * a ) )
                         : ( ( segc[ 0 ] - b ) / ( 2 * a ) - wd ) );
   }
  else if( r.form == 3 ) {                  // TENT q*=p-dpk
   for( const auto & pc2 : F )
    { rel( pc2.left + dpk ); rel( pc2.right + dpk ); }      // F-cross
   for( int m = 0 ; m < K ; ++m )
    if( segcum[ m ] > 1e-14 ) rel( Bmax / segcum[ m ] );    // Bmax=K_m
   rel( pmin + Bmax ); rel( pmax - Bmax );  // Bmax = A (TENT on/off)
   const int i = findF( pc - dpk );
   const double a = F[ i ].alfa , b = F[ i ].beta;
   int m = 0;
   for( ; m < K ; ++m ) if( Bmax <= segcum[ m ] * pc + 1e-9 ) break;
   if( m >= K ) m = K - 1; const double sc = segc[ m ];
   if( a > 1e-16 ) { rel( dpk - ( b + sc ) / ( 2 * a ) );
                     rel( dpk - ( b - sc ) / ( 2 * a ) ); }
   }
  else if( ( r.form == 4 ) || ( r.form == 5 ) ) {  // BAND / KBAND:
   for( const auto & pc2 : F )                     // qc = sl*p + ic
    { rel( ( pc2.left - r.ic ) / r.sl );
      rel( ( pc2.right - r.ic ) / r.sl ); }                 // F-cross
   rel( ( domL - r.ic ) / r.sl );
   rel( ( domR - r.ic ) / r.sl );           // qc hits domain
   rel( pmin + Bmax ); rel( pmax - Bmax );  // A = Bmax (BAND/KBAND <-> TENT)
   if( r.form == 4 ) {                      // BAND stat-entry (F'(qc)=0 side)
    const int i = findF( r.sl * pc + r.ic );
    const double a = F[ i ].alfa , b = F[ i ].beta;
    const double A = std::min( pc - pmin , pmax - pc );
    int mA = 0;
    for( ; mA < K ; ++mA ) if( A <= segcum[ mA ] * pc + 1e-9 ) break;
    if( mA >= K ) mA = K - 1;
    const double sgn = ( r.m == 1 ) ? 1.0 : -1.0;
    if( a > 1e-16 )
     { rel( ( ( -segc[ mA ] * sgn - b ) / ( 2 * a ) - r.ic ) / r.sl );
       rel( ( ( 0.0 - b ) / ( 2 * a ) - r.ic ) / r.sl ); }
    }
   if( r.form == 5 ) {                      // KBAND stat-entry + tent-meet
    if( segcum[ r.m ] > 1e-14 ) rel( Bmax / segcum[ r.m ] );
    const int i = findF( r.sl * pc + r.ic );
    const double a = F[ i ].alfa , b = F[ i ].beta;
    const double sgn = ( r.sl > 1.0 ) ? 1.0 : -1.0;
    const double mlo = - segc[ r.m ] * sgn;
    const double mhi = - ( ( r.m + 1 < K ) ? segc[ r.m + 1 ] : 0.0 ) * sgn;
    if( a > 1e-16 ) { rel( ( ( mlo - b ) / ( 2 * a ) - r.ic ) / r.sl );
                      rel( ( ( mhi - b ) / ( 2 * a ) - r.ic ) / r.sl ); }
    }
   }
  return( best );
  };
 // nearest event that may CHANGE THE FORM (not merely the F-piece or reward
 // segment, which GformR/curvR absorb): window-reach, stationary entry, the
 // B=A / Bmax=A corr on-off, tent peak, and qc leaving the domain.
 // Everything NOT listed here (F-cross, A-kink, A=K_m, B=K_m, Bmax=K_m)
 // keeps the form, so the sweep advances through it WITHOUT re-classifying.
 // Over-including is safe (it only re-classifies more often); missing a
 // form change is not.
 auto next_form_event = [ & ]( const Reg & r , double pc ) -> double {
  double best = phi;
  auto rel = [ & ]( double pe ) {
   if( ( pe > pc + 1e-9 ) && ( pe < best ) ) best = pe; };
  if( r.form == 0 ) {                               // PIN: q* starts to move
   const double qbar = r.qbar;
   rel( qbar + wu ); rel( qbar - wd );              // window-reach -> WUP/WDN
   rel( qbar + dpk );                               // tent peak -> TENT
   rel( 0.5 * ( ramp_up + qbar + pmin ) );          // B=A -> BAND
   rel( 0.5 * ( pmax - ramp_down + qbar ) );
   for( int m = 0 ; m < K ; ++m ) {                 // B=K_m -> KBAND
    if( std::abs( 1 - segcum[ m ] ) > 1e-14 )
     rel( ( qbar - ramp_down ) / ( 1 - segcum[ m ] ) );
    rel( ( ramp_up + qbar ) / ( 1 + segcum[ m ] ) ); }
   }
  else {                                            // MOVING q*=sl*p+ic
   rel( ( domL - r.ic ) / r.sl );                   // q* hits domain -> PIN
   rel( ( domR - r.ic ) / r.sl );
   if( std::abs( r.sl - 1.0 ) > 1e-12 ) {  // q* hits a window edge -> WUP/WDN
    rel( ( -wu - r.ic ) / ( r.sl - 1.0 ) );
    rel( ( wd - r.ic ) / ( r.sl - 1.0 ) ); }
   rel( pmin + Bmax ); rel( pmax - Bmax );  // Bmax = A (<-> TENT / corr off)
   // stationary point enters (-> PIN): the moving locus q*(p) sweeps
   // through F as p grows, so the entry can land in ANY F-piece (and any
   // reward segment), not just the one at pc. Check every piece x segment;
   // over-generation is safe.
   for( const auto & fp : F ) {
    const double a = fp.alfa , b = fp.beta;
    if( a <= 1e-16 ) continue;
    if( r.form == 1 )                                                  // WUP
     rel( wu - ( b + segc[ 0 ] ) / ( 2 * a ) );
    else if( r.form == 2 )                                             // WDN
     rel( ( segc[ 0 ] - b ) / ( 2 * a ) - wd );
    else if( r.form == 3 )                                             // TENT
     for( int m = 0 ; m < K ; ++m ) { const double sc = segc[ m ];
      rel( dpk - ( b + sc ) / ( 2 * a ) );
      rel( dpk - ( b - sc ) / ( 2 * a ) ); }
    else if( r.form == 4 ) {                                           // BAND
     const double sgn = ( r.m == 1 ) ? 1.0 : -1.0;
     for( int mA = 0 ; mA < K ; ++mA )
      rel( ( ( -segc[ mA ] * sgn - b ) / ( 2 * a ) - r.ic ) / r.sl );
     rel( ( ( 0.0 - b ) / ( 2 * a ) - r.ic ) / r.sl ); }
    else {                                                            // KBAND
     const double sgn = ( r.sl > 1.0 ) ? 1.0 : -1.0;
     const double mlo = - segc[ r.m ] * sgn;
     const double mhi = - ( ( r.m + 1 < K ) ? segc[ r.m + 1 ] : 0.0 ) * sgn;
     rel( ( ( mlo - b ) / ( 2 * a ) - r.ic ) / r.sl );
     rel( ( ( mhi - b ) / ( 2 * a ) - r.ic ) / r.sl ); }
    }
   // PIN-at-F-kink: the moving locus q*(p) can WEDGE at a convex F-kink q_k
   // (the min pins there -> PIN) as it crosses it, a jump at an
   // F-crossing, invisible to the stationary test. Pin iff 0 in the
   // subgradient [F'(q_k^-)+corr'(q_k^-), F'(q_k^+)+corr'(q_k^+)] at the p
   // where q*(p)=q_k. corr' is DIRECTIONAL (corrDer, hoisted above): corr
   // itself kinks at q_k (e.g. WUP: B=0 for q<q_k, B>0 for q>q_k).
   for( std::size_t i = 0 ; i + 1 < F.size() ; ++i ) {
    const double qk = F[ i ].right;
    const double pe = ( qk - r.ic ) / r.sl;
    if( pe <= pc + 1e-9 ) continue;
    const double eps = 1e-6 * std::max( 1.0 , std::abs( qk ) );
    const double gL = 2 * F[ i ].alfa * qk + F[ i ].beta +
                      corrDer( qk - eps , pe );
    const double gR = 2 * F[ i + 1 ].alfa * qk + F[ i + 1 ].beta +
                      corrDer( qk + eps , pe );
    const double sc = std::max( 1.0 ,
                                std::max( std::abs( gL ) , std::abs( gR ) ) );
    if( ( gL <= 1e-6 * sc ) && ( gR >= -1e-6 * sc ) )
     rel( pe );                                     // straddle -> pin
    }
   }
  return( best );
  };
 // EXACT closed-form (a,b,c) per regime (no 2-endpoint interpolation
 // b=(rr-rl)/(pnext-p), whose subtraction cancels on narrow over-generated
 // pieces). Within a regime corr(q*(p),p)=reward(p,B(p))-reward(p,A(p)) is
 // LINEAR in p (B,A linear; reward linear in its argument with the active
 // segment fixed), so G(p)=aF*(sl*p+ic)^2+bF*(sl*p+ic)+cF
 // + (m_corr*p + k_corr) is exactly quadratic.
 // reward decomposition: in active segment m (kappa_{m-1}<=H<kappa_m,
 // kappa_j=segcum[j]*p) reward(p,H)=Rbase[m]*p+segc[m]*H; beyond the last
 // cap it is Rfull*p.
 double Rbase[ 3 ] = { 0 , 0 , 0 } , Rfull = 0;
 { double prev = 0 , acc = 0;
   for( int j = 0 ; j < K ; ++j ) {
    Rbase[ j ] = acc - segc[ j ] * prev;
    acc += segc[ j ] * ( segcum[ j ] - prev );
    prev = segcum[ j ]; }
   Rfull = acc; }
 // linear coeffs (coefP, cst) of reward(p, H(p)=Hs*p+Hi) at p_eval:
 // reward=coefP*p+cst
 auto rewardLin = [ & ]( double Hs , double Hi , double p_eval ,
                         double & coefP , double & cst ) {
  const double H = Hs * p_eval + Hi;
  if( H <= 1e-12 ) { coefP = 0; cst = 0; return; }
  if( ( K > 0 ) && ( H >= segcum[ K - 1 ] * p_eval - 1e-12 ) )
   { coefP = Rfull; cst = 0; return; }     // beyond total cap: marginal 0
  int m = 0;
  for( ; m < K ; ++m ) if( H <= segcum[ m ] * p_eval + 1e-9 ) break;
  if( m >= K ) m = K - 1;
  coefP = Rbase[ m ] + segc[ m ] * Hs; cst = segc[ m ] * Hi;
  };
 // exact (a,b,c) of the piece for a locus q*(p)=sl*p+ic (mpQP critical
 // region), evaluated at interior pm. Curvature aF*sl^2 anchored to the
 // z-piece at q*(pm); b,c include the corr = g(p,min(A,B)) linear coeffs
 // (rewardLin). Callable directly by the mpQP sweep (with sl,ic from the
 // finite-diff locus id) or via closedABC(Reg).
 auto closedABCst = [ & ]( double sl , double ic , double pm ,
                          double & a , double & b , double & c ) {
  const double qm = sl * pm + ic;
  const PieceQuad & pc = F[ findF( qm ) ];
  const double aF = pc.alfa , bF = pc.beta , cF = pc.gamma;
  a = aF * sl * sl;
  double mcorr = 0 , kcorr = 0;
  const double A = std::min( pm - pmin , pmax - pm );
  if( A > 0 ) {                           // corr = g(p, H(p)), H = min(A,B)
   const double d = pm - qm;
   double B = std::min( ramp_up - d , ramp_down + d ); if( B < 0 ) B = 0;
   double Hs , Hi;
   if( B < A - 1e-12 ) {                  // ramp binds: H = B on its branch
    if( d > dpk ) { Hs = sl - 1; Hi = ramp_up   + ic; }   // upper: B = ru - d
    else          { Hs = 1 - sl; Hi = ramp_down - ic; }   // lower: B = rd + d
    }
   else {                                 // ramp slack: H = A on its side
    if( pm < midA ) { Hs = 1; Hi = -pmin; } else { Hs = -1; Hi = pmax; }
    }
   rewardLin( Hs , Hi , pm , mcorr , kcorr );
   }
  b = 2 * aF * sl * ic + bF * sl + mcorr;
  c = aF * ic * ic + bF * ic + cF + kcorr;
  };
 auto closedABC = [ & ]( const Reg & r , double pm ,
                         double & a , double & b , double & c ) {
  closedABCst( ( r.form == 0 ) ? 0.0 : r.sl ,
               ( r.form == 0 ) ? r.qbar : r.ic , pm , a , b , c );
  };
 // VERIFY certificate: q*(pp) is the true min of convex phi=F+corr over
 // [qL,qR] iff it cannot be improved either way (sound by convexity; O(1),
 // no Gmin). This lets form-reuse re-classify ONLY at genuine form changes
 // -- no need for a provably complete next_form_event.
 // phi'(q^-/+) = F'(q^-/+) + corr'(q^-/+) (directional).
 auto formOK = [ & ]( const Reg & rr , double pp ) -> bool {
  if( rr.form < 0 ) return( false );
  const double qL = std::max( pp - wu , domL );
  const double qR = std::min( pp + wd , domR );
  const double q = qstarR( rr , pp );
  if( ( q < qL - 1e-7 ) || ( q > qR + 1e-7 ) )
   return( false );                                 // q* left the window
  const double e = 1e-7 * std::max( 1.0 , std::abs( q ) );
  // q-e and q+e are in the same F-piece, unless at a kink
  const PieceQuad & pL = F[ findF( q - e ) ];
  const PieceQuad & pR = ( q + e <= pL.right + 1e-12 ) ? pL
                                                       : F[ findF( q + e ) ];
  const double dphiL = 2 * pL.alfa * q + pL.beta +
                       corrDer( q - e , pp );                    // phi'(q^-)
  const double dphiR = 2 * pR.alfa * q + pR.beta +
                       corrDer( q + e , pp );                    // phi'(q^+)
  const double tol = 1e-7 * std::max( 1.0 ,
                                      std::abs( dphiL ) + std::abs( dphiR ) );
  const bool leftok = ( q <= qL + 1e-7 ) ||
                      ( dphiL <= tol );             // can't improve left
  const bool rightok = ( q >= qR - 1e-7 ) ||
                       ( dphiR >= -tol );           // can't improve right
  return( leftok && rightok );
  };
 const double pctol = std::getenv( "TUEDPS_CTOL" )
                      ? std::atof( std::getenv( "TUEDPS_CTOL" ) ) : 1e-6;
 // sweep p left->right, emit one closed-form piece per regime, then coalesce
 // adjacent identical pieces (over-generated events split real pieces).
 const bool bcchk = std::getenv( "TUEDPS_BCCHECK" );  // hoisted out of the
 const bool formreuse = std::getenv( "TUEDPS_FORMREUSE" );      // loop
 const bool frdbg = std::getenv( "TUEDPS_FRDBG" );
 auto sweepParam = [ & ]( PQFun & dst , bool reuse ) {
  dst.clear();
  if( std::getenv( "TUEDPS_FINCHK" ) )       // is the INPUT F convex?
   for( std::size_t i = 1 ; i < F.size() ; ++i ) {
    const double bp = F[ i ].left;
    const double sL = 2 * F[ i - 1 ].alfa * bp + F[ i - 1 ].beta ,
                 sR = 2 * F[ i ].alfa * bp + F[ i ].beta;
    const double vL = eval_piece( F[ i - 1 ] , bp ) ,
                 vR = eval_piece( F[ i ] , bp );
    if( std::abs( vL - vR ) > 1e-6 * std::max( 1.0 , std::abs( vR ) ) )
     std::cerr << "FINCHK input-F DISCONT at q=" << bp
               << " jump=" << ( vR - vL ) << "\n";
    if( sR - sL < -1e-6 * std::max( 1.0 , std::abs( sL ) ) )
     std::cerr << "FINCHK input-F slopedrop at q=" << bp << " sL=" << sL
               << " sR=" << sR << " drop=" << ( sL - sR ) << "\n";
    }
  double p = plo; int guard = 0;
  // FORM-REUSE fast path (reuse=true): carry the classified form across
  // pieces and re-classify ONLY when the O(1) verify certificate (formOK)
  // says it is no longer optimal, i.e. at genuine form changes. Sound by
  // convexity, so no dependence on a provably-complete next_form_event.
  // next_event bounds form-constancy, so a form that verifies at a piece
  // midpoint is the true form on the whole piece.
  Reg r = classify( p );
  int prevform = -99;
  while( ( p < phi - 1e-9 ) && ( guard++ < 200000 ) ) {
   if( reuse ) {                     // verify the carried form; refresh if
    double pn = ( r.form < 0 ) ? phi : next_event( r , p );        // stale
    if( pn > phi ) pn = phi;
    const double w = pn - p;                // verify q* optimal (convex phi
    const int nv = std::getenv( "TUEDPS_V1" ) ? 1 : 3;   // at continuous F
    bool ok = ( r.form >= 0 );                           // => 1 point sound
    if( ok && nv == 1 ) ok = formOK( r , p + 0.5 * w );
    else if( ok ) ok = formOK( r , p + 0.25 * w ) && formOK( r , p + 0.5 * w )
                       && formOK( r , p + 0.75 * w );
    if( ! ok ) r = classify( p );
    }
   else
    r = classify( p );               // trusted: classify every piece
   double pnext = ( r.form < 0 ) ? phi : next_event( r , p );
   if( pnext > phi ) pnext = phi;
   if( pnext <= p + 1e-12 ) break;
   if( r.form >= 0 ) {               // emit exact closed-form piece
    const double pm = 0.5 * ( p + pnext );
    double a , b , c;
    closedABC( r , pm , a , b , c );
    if( reuse && frdbg ) {  // enumerate where the reused form is wrong
     const Reg rc = classify( pm );
     if( rc.form != r.form ) {
      const double gr = std::abs( ( a * pm + b ) * pm + c
                                  - GformR( rc , pm ) );
      if( gr > 1.0 ) {      // real value error -> dump formOK internals
       const double q = qstarR( r , pm );
       const double qL = std::max( pm - wu , domL ) ,
                    qR = std::min( pm + wd , domR );
       const double e = 1e-7 * std::max( 1.0 , std::abs( q ) );
       const double dL = 2 * F[ findF( q - e ) ].alfa * q +
                         F[ findF( q - e ) ].beta + corrDer( q - e , pm );
       const double dR = 2 * F[ findF( q + e ) ].alfa * q +
                         F[ findF( q + e ) ].beta + corrDer( q + e , pm );
       const double qt = qstarR( rc , pm );
       const double h = 1e-4;               // finite-diff of ACTUAL phi
       const double phib = evalF( q ) + corr( q , pm );
       const double fdL = ( phib - ( evalF( q - h )
                                     + corr( q - h , pm ) ) ) / h;
       const double fdR = ( ( evalF( q + h ) + corr( q + h , pm ) )
                            - phib ) / h;   // fdL/fdR = phi'(qbar^-/+)
       std::cerr << "PINBUG " << r.form << "->" << rc.form << " qbar=" << q
                 << " qtrue=" << qt << " | analytic dL=" << dL << " dR=" << dR
                 << " | finitediff fdL=" << fdL << " fdR=" << fdR << "\n"; }
      }
     }
    if( bcchk ) {         // dev: verify vs exact GformR
     const double pa = p + 0.25 * ( pnext - p ) ,
                  pb = p + 0.75 * ( pnext - p );
     const double e0 = std::abs( ( a * pa + b ) * pa + c - GformR( r , pa ) );
     const double e1 = std::abs( ( a * pb + b ) * pb + c - GformR( r , pb ) );
     const double sc = std::max( 1.0 , std::abs( GformR( r , pm ) ) );
     if( ( e0 > 1e-7 * sc ) || ( e1 > 1e-7 * sc ) )
      std::cerr << "BCCHK form=" << r.form << " e0=" << e0 << " e1=" << e1
                << " w=" << ( pnext - p )
                << " [" << p << "," << pnext << "]\n";
     }
    if( std::getenv( "TUEDPS_RCHK" ) ) {   // does the piece SPAN a
     auto sig = [ & ]( double x , int & f , int & as , int & bs , int & bi ,
                       int & br , int & asd ) {   // sub-regime change?
      const double q = qstarR( r , x );
      f = findF( q );
      const double A = std::min( x - pmin , pmax - x );
      asd = ( x < midA ) ? 1 : 0;
      const double d = x - q , Bup = ramp_up - d , Blo = ramp_down + d;
      double B = std::min( Bup , Blo ); if( B < 0 ) B = 0;
      br = ( Bup < Blo ) ? 1 : 0;
      bi = ( ( B > 1e-12 ) && ( B < A - 1e-12 ) && ( A > 0 ) ) ? 1 : 0;
      as = 0; for( ; as < K ; ++as ) if( A <= segcum[ as ] * x + 1e-9 ) break;
      bs = -1; if( bi ) { bs = 0;
       for( ; bs < K ; ++bs ) if( B <= segcum[ bs ] * x + 1e-9 ) break; }
      };
     const double e = 1e-6 * std::max( 1.0 , std::abs( pnext - p ) );
     int f1 , a1 , b1 , i1 , r1 , d1 , f2 , a2 , b2 , i2 , r2 , d2;
     sig( p + e , f1 , a1 , b1 , i1 , r1 , d1 );
     sig( pnext - e , f2 , a2 , b2 , i2 , r2 , d2 );
     const bool tentbr = ( r.form == 3 ) && ( f1 == f2 ) && ( a1 == a2 ) &&
                         ( b1 == b2 ) &&( i1 == i2 ) &&
                         ( d1 == d2 );   // harmless TENT branch flicker
     if( ( ! tentbr ) &&
         ( ( f1 != f2 ) || ( a1 != a2 ) || ( b1 != b2 ) || ( i1 != i2 ) ||
           ( r1 != r2 ) || ( d1 != d2 ) ) )
      std::cerr << "RCHK form=" << r.form << " [" << p << "," << pnext << "]"
                << " F:" << f1 << ">" << f2 << " As:" << a1 << ">" << a2
                << " Bs:" << b1 << ">" << b2 << " bite:" << i1 << ">" << i2
                << " br:" << r1 << ">" << r2 << " Asd:" << d1 << ">" << d2
                << " midA=" << midA
                << " qL=" << std::max( p - wu , domL ) << "\n";
     }
    if( std::getenv( "TUEDPS_JOINT" ) && ! dst.empty() ) {  // convexity at
     const PieceQuad & pv = dst.back();                     // the joint
     const double vL = eval_piece( pv , p ) , vR = ( a * p + b ) * p + c;
     const double sL = 2 * pv.alfa * p + pv.beta , sR = 2 * a * p + b;
     const double sc = std::max( 1.0 , std::abs( vR ) );
     if( ( std::abs( vL - vR ) > 1e-6 * sc ) ||
         ( sR - sL < -1e-6 * std::max( 1.0 , std::abs( sL ) ) ) ) {
      const double hh = 1e-5 * std::max( 1.0 ,
                                       std::abs( p ) );  // true G' via exact
      const double tgL = ( GformR( r , p ) -                     // GformR
                           GformR( r , p - hh ) ) / hh;
      const double tgR = ( GformR( r , p + hh ) - GformR( r , p ) ) / hh;
      std::cerr << "JOINT " << prevform << "->" << r.form << " bp=" << p
                << " vjump=" << ( vR - vL ) << " | closedABC sL=" << sL
                << " sR=" << sR
                << " | trueG' tgL=" << tgL << " tgR=" << tgR << "\n"; }
     }
    prevform = r.form;
    dst.push_back( { a , b , c , p , pnext } );
    }
   p = pnext;
   }
  if( dst.size() > 1 ) {                           // coalesce
   PQFun m; m.reserve( dst.size() ); m.push_back( dst.front() );
   for( std::size_t i = 1 ; i < dst.size() ; ++i ) {
    const PieceQuad & b = m.back() , & c = dst[ i ];
    const double x1 = c.right , x2 = 0.5 * ( c.left + c.right );
    const double e1 = std::abs( eval_piece( b , x1 ) - eval_piece( c , x1 ) );
    const double e2 = std::abs( eval_piece( b , x2 ) - eval_piece( c , x2 ) );
    const double sc = std::max( 1.0 , std::abs( eval_piece( c , x2 ) ) );
    if( ( e1<=pctol * sc ) && ( e2<=pctol * sc ) ) m.back().right = c.right;
    else m.push_back( c );
    }
   dst.swap( m );
   }
  if( std::getenv( "TUEDPS_FCHK" ) )    // convexity audit of the built VF
   for( std::size_t i = 0 ; i < dst.size() ; ++i ) {
    if( dst[ i ].alfa < -1e-9 )
     std::cerr << "FCHK negcurv a=" << dst[ i ].alfa
               << " [" << dst[ i ].left << "," << dst[ i ].right << "]\n";
    if( i > 0 ) {
     const double bp = dst[ i ].left;
     const double vlo = eval_piece( dst[ i - 1 ] , bp ) ,
                  vhi = eval_piece( dst[ i ] , bp );
     const double slo = 2 * dst[ i - 1 ].alfa * bp + dst[ i - 1 ].beta;
     const double shi = 2 * dst[ i ].alfa * bp + dst[ i ].beta;
     const double sc = std::max( 1.0 , std::abs( vhi ) );
     if( std::abs( vlo - vhi ) > 1e-6 * sc )
      std::cerr << "FCHK discont at " << bp
                << " jump=" << ( vhi - vlo ) << "\n";
     if( shi - slo < -1e-6 * std::max( 1.0 , std::abs( slo ) ) )
      std::cerr << "FCHK slopedrop at " << bp << " sL=" << slo
                << " sR=" << shi << " drop=" << ( slo - shi ) << "\n";
     }
    }
  };
 // ==================================================================
 // mpQP SWEEP: principled multiparametric-QP construction validated in
 // Python (scratchpad/mpqp_*.py) to match brute-force to ~1e-13 AND be
 // exactly convex for K=1,2 and cap-variants. Locus id by finite difference
 // of q*(p) (Gmin), a COMPLETE event set (crossings of the current locus
 // with all critical-region loci, z-kinks and regime boundaries), and the
 // exact closedABCst piece. Production default; classify+next_event above
 // is the gated legacy path.
 auto sweepMPQP = [ & ]( PQFun & dst ) {
  dst.clear();
  static thread_local std::vector< std::pair< double , double > > Lall;
  Lall.clear();
  Lall.push_back( { 1.0 , -wu } );                               // window
  Lall.push_back( { 1.0 , wd } );
  Lall.push_back( { 0.0 , domL } );                              // domain
  Lall.push_back( { 0.0 , domR } );
  Lall.push_back( { 1.0 , -dpk } );                              // tent peak
  Lall.push_back( { 2.0 , -pmin - ramp_up } );                   // B=A rise-A
  Lall.push_back( { 0.0 , ramp_down + pmin } );
  Lall.push_back( { 0.0 , pmax - ramp_up } );                    // B=A fall-A
  Lall.push_back( { 2.0 , ramp_down - pmax } );
  for( int m = 0 ; m < K ; ++m ) {                               // B=kappa_m
   Lall.push_back( { 1.0 + segcum[ m ] , -ramp_up } );
   Lall.push_back( { 1.0 - segcum[ m ] , ramp_down } ); }
  // NOTE: the z-piece stationaries are NOT global loci, only the CURRENT
  // piece's matter for the next event (q* must cross an F-kink, already an
  // event, to reach another piece). Adding them per-iteration avoids
  // O(|F|) spurious over-generation.
  static thread_local std::vector< double > ZK; ZK.clear();
  for( const auto & pc : F ) ZK.push_back( pc.left );
  ZK.push_back( F.back().right );
  static thread_local std::vector< double > RB; RB.clear();
  RB.push_back( midA ); RB.push_back( pmin + ramp_up );
  RB.push_back( domR - ramp_down );
  if( ( wu != ramp_up ) || ( wd != ramp_down ) ) {  // a narrower window
   RB.push_back( pmin + wu ); RB.push_back( domL + wu );
   RB.push_back( domR - wd ); }
  RB.push_back( pmin + Bmax ); RB.push_back( pmax - Bmax );
  RB.push_back( pmin ); RB.push_back( pmax );
  for( int m = 0 ; m < K ; ++m ) {
   if( std::abs( 1 - segcum[ m ] ) > 1e-12 )
    RB.push_back( pmin / ( 1 - segcum[ m ] ) );
   RB.push_back( pmax / ( 1 + segcum[ m ] ) );
   if( segcum[ m ] > 1e-12 ) RB.push_back( Bmax / segcum[ m ] ); }
  std::sort( RB.begin() , RB.end() );  // p-thresholds: sorted so the sweep
                                       // can cursor them
  // FORM-CARRY fast path: carry the finite-diff locus (s,t) across F-piece
  // boundaries (z-kinks: form unchanged, only the underlying quadratic
  // differs -> closedABCst re-derives the piece with NO Gmin). A pre-emit
  // certificate checks the locus over the WHOLE piece [p,pnext] (3 probes);
  // a stale carry re-acquires (fresh finite diff), and a fresh locus that
  // still fails an unforeseen sub-event is bisected toward p (progress
  // guaranteed). Re-Gmin fires only at genuine form changes, not at every
  // z-kink.
  auto okPiece = [ & ]( const Reg & rr , double pa , double pb ) -> bool {
   // q*(p)=s*p+t is monotone (s>0), so the carried (F-piece x corr-segment)
   // form is active on a CONTIGUOUS p-interval: if the optimality
   // certificate holds at both endpoints it holds throughout, so two probes
   // bracket validity (midpoint redundant).
   const double w = pb - pa;
   return( formOK( rr , pa + 0.02 * w ) && formOK( rr , pa + 0.98 * w ) ); };
  // ROBUST TAIL FILL: cover [p0,phi] with an adaptive-linear Gmin
  // subdivision. Used when the finite-diff form-ID fails (okPiece rejects /
  // bisection collapses): a plain BREAK would TRUNCATE the value function's
  // domain and silently DROP DP states (landing powers above p0), making
  // the DP miss reachable optimal schedules; falling back to
  // exact-at-samples linear pieces keeps the full [plo,phi] domain always
  // emitted. Gmin is exact wherever the ramp window is non-empty.
  auto fillTail = [ & ]( double p0 ) {
   if( p0 >= phi - 1e-9 ) return;
   const double vmid = Gmin( 0.5 * ( p0 + phi ) ).first;
   if( vmid >= TUEDPINF ) return;   // no reachable landing: nothing to add
   const double tol = 1e-7 * std::max( 1.0 , std::abs( vmid ) );
   static thread_local std::vector< std::pair< double , double > > stk;
   stk.clear();
   stk.push_back( { p0 , phi } ); int g2 = 0;
   while( ! stk.empty() && ( g2++ < 200000 ) ) {
    const auto ab = stk.back(); stk.pop_back();
    const double a = ab.first , b = ab.second , m = 0.5 * ( a + b );
    const double va = Gmin( a ).first , vb = Gmin( b ).first ,
                 vm = Gmin( m ).first;
    if( ( va >= TUEDPINF ) || ( vb >= TUEDPINF ) ) continue;
    const double lin = va + ( vb - va ) * ( m - a ) / ( b - a );
    if( ( b - a < 1e-5 ) ||
        ( ( vm < TUEDPINF ) && ( std::abs( vm - lin ) <= tol ) ) ) {
     const double sl = ( vb - va ) / ( b - a );
     dst.push_back( { 0.0 , sl , va - sl * a , a , b } ); }
    else { stk.push_back( { m , b } );
           stk.push_back( { a , m } ); }   // left first (in order)
    }
   };
  const bool doSnap = ! std::getenv( "TUEDPS_NOSNAP" );  // snap locus to
                                           // exact loci (A/B only)
  double p = plo; int guard = 0; bool haveForm = false; double s = 0 , t = 0;
  int zkptr = 0;   // sweep-line cursor into the sorted z-kinks
                   // (q*=s*p+t is monotone, s>0)
  int rbptr = 0;   // sweep-line cursor into the sorted regime boundaries
                   // (p is monotone non-decr)
  static thread_local std::vector< double > LX;
  LX.clear(); int lxptr = 0;               // locus-crossing p's
  while( ( p < phi - 1e-9 ) && ( guard++ < 200000 ) ) {
#if TUEDPS_PROFILE
   ++g_param_calls;                                    // raw sweep steps
#endif
   bool justAcq = false;
   if( ! haveForm ) {          // (re)acquire the locus by robust finite diff
#if TUEDPS_PROFILE
    ++g_param_bad;                         // re-acquisitions (2 Gmin each)
#endif
    const double e1 = 1e-7 * std::max( 1.0 , std::abs( p ) );
    const auto g1 = Gmin( p + e1 );
    if( g1.first >= TUEDPINF ) {
#if TUEDPS_PROFILE
     if( std::getenv( "TUEDPS_BRKLOG" ) && ( p < phi - 1.0 ) )
      std::cerr << "MPBREAK reason=Gmin-INF p=" << p << " phi=" << phi
                << " domL=" << domL << " domR=" << domR << " ru=" << ramp_up
                << " rd=" << ramp_down << "\n";
#endif
     fillTail( p ); break; }
    const double q1 = g1.second , pe = p + e1;
    const double e2 = 2 * e1 , q2 = Gmin( p + e2 ).second;
    s = ( q2 - q1 ) / e1; t = q1 - s * pe; haveForm = true; justAcq = true;
    // SNAP the finite-diff locus to the EXACT analytic locus: the ~1e-6
    // slope noise would give a ~1e-4 intercept error that ACCUMULATES over
    // the horizon. PIN (|s|~0): the pin sits at the exact argmin q1.
    // MOVING: pick the Lall (slope,intercept) whose line best fits both
    // exact argmin probes (pe,q1),(p+e2,q2).
    if( doSnap ) {
     if( std::abs( s ) < 0.5 ) { s = 0.0; t = q1; }
     else {
      double bd = 1e300 , bs = s , bt = t;
      for( const auto & L : Lall ) {
       const double dd = std::abs( L.first * pe + L.second - q1 )
                       + std::abs( L.first * ( p + e2 ) + L.second - q2 );
       if( dd < bd ) { bd = dd; bs = L.first; bt = L.second; } }
      if( bd <= 1e-6 * std::max( 1.0 , std::abs( q1 ) ) ) { s = bs; t = bt; }
      }
     }
    }
   if( justAcq ) {   // (s,t) just (re)set: cache the sorted locus-crossing
    LX.clear();      // p-values, reset cursor
    for( const auto & L : Lall ) if( std::abs( L.first - s ) > 1e-12 ) {
     const double pp = ( L.second - t ) / ( s - L.first );
     if( pp > p + 1e-9 ) LX.push_back( pp ); }
    std::sort( LX.begin() , LX.end() ); lxptr = 0;
    }
   const double qcur = s * p + t;   // argmin at p under the carried locus
   double pnext = phi;
   while( ( lxptr < ( int ) LX.size() ) &&
          ( LX[ lxptr ] <= p + 1e-9 ) ) ++lxptr;              // drop passed
   if( ( lxptr < ( int ) LX.size() ) && ( LX[ lxptr ] < pnext ) )
    pnext = LX[ lxptr ];                                    // nearest > p
   if( std::abs( s ) > 1e-12 ) {
    // z-kink crossings: q*(p)=s*p+t rises monotonically (s>0), so z-kinks
    // are crossed in increasing order, only the NEXT z-kink above qcur
    // can end this piece. Advance a cursor instead of rescanning all of ZK
    // every step (an O(|F|)-per-step scan would dominate the sweep).
    const int nzk = ( int ) ZK.size();
    if( justAcq ) {          // form (hence q*) jumped: relocate the cursor
     int lo = 0 , hi = nzk;
     while( lo < hi ) { const int m = ( lo + hi ) / 2;
      if( ZK[ m ] <= qcur + 1e-9 ) lo = m + 1; else hi = m; }
     zkptr = lo; }
    else while( ( zkptr < nzk ) && ( ZK[ zkptr ] <= qcur + 1e-9 ) ) ++zkptr;
    if( zkptr < nzk ) {
     const double pp = ( ZK[ zkptr ] - t ) / s;
     if( ( pp > p + 1e-9 ) && ( pp < pnext ) ) pnext = pp; }
    const PieceQuad & zp = F[ findF( qcur ) ];  // current z-piece
    if( zp.alfa > 1e-15 ) {                     // stationaries only
     const double base = -zp.beta / ( 2 * zp.alfa );
     double sq = ( base - t ) / s;
     if( ( sq > p + 1e-9 ) && ( sq < pnext ) ) pnext = sq;
     for( int m = 0 ; m < K ; ++m ) {
      sq = ( base - segc[ m ] / ( 2 * zp.alfa ) - t ) / s;
      if( ( sq > p + 1e-9 ) && ( sq < pnext ) ) pnext = sq;
      sq = ( base + segc[ m ] / ( 2 * zp.alfa ) - t ) / s;
      if( ( sq > p + 1e-9 ) && ( sq < pnext ) ) pnext = sq; }
     }
    }
   while( ( rbptr < ( int ) RB.size() ) &&
          ( RB[ rbptr ] <= p + 1e-9 ) )
    ++rbptr;                                               // drop passed RBs
   if( ( rbptr < ( int ) RB.size() ) && ( RB[ rbptr ] < pnext ) )
    pnext = RB[ rbptr ];                                    // nearest RB > p
   if( pnext > phi ) pnext = phi;
   if( pnext <= p + 1e-12 ) {
#if TUEDPS_PROFILE
    if( std::getenv( "TUEDPS_BRKLOG" ) && ( p < phi - 1.0 ) )
     std::cerr << "MPBREAK reason=noadvance p=" << p << " phi=" << phi
               << " s=" << s << " t=" << t << " pnext=" << pnext << "\n";
#endif
    fillTail( p ); break; }
   double pm = 0.5 * ( p + pnext );
   const Reg rr = ( std::abs( s ) < 1e-12 ) ? Reg{ 0 , t , 0 , 0 , 0 }
                                          : Reg{ 1 , 0 , s , t , 0 };
   bool bisected = false;
   if( ! okPiece( rr , p , pnext ) ) {
    if( ! justAcq )
     { haveForm = false; continue; }  // stale carry: re-acquire fresh at p
    int bg = 0;              // fresh locus valid at p but not across
                             // [p,pnext]: an event the scan cannot see
                             // -> bisect toward p
    while( ( ! okPiece( rr , p , pnext ) ) && ( pnext > p + 1e-9 ) &&
           ( bg++ < 80 ) )
     pnext = 0.5 * ( p + pnext );
    if( pnext <= p + 1e-9 ) {
#if TUEDPS_PROFILE
     if( std::getenv( "TUEDPS_BRKLOG" ) && ( p < phi - 1.0 ) )
      std::cerr << "MPBREAK reason=bisect-collapse p=" << p << " phi=" << phi
                << " s=" << s << " t=" << t << " justAcq=" << justAcq << "\n";
#endif
     fillTail( p ); break; }
    pm = 0.5 * ( p + pnext ); bisected = true;
    }
   double a , b , c; closedABCst( s , t , pm , a , b , c );
#if TUEDPS_PROFILE
   if( std::getenv( "TUEDPS_MPTRACE" )
       && ( F.size() ==
            ( std::size_t ) std::atoi( getenv( "TUEDPS_MPTRACE" ) ) )
       && ( p >= 213.0 ) && ( p <= 215.0 ) )
    std::cerr.precision( 10 ) ,
    std::cerr << "MPTR p=" << p << " s=" << s << " t=" << t
     << " pnext=" << pnext << " a=" << a << " qcur=" << ( s * p + t )
     << " justAcq=" << justAcq << " bisect=" << bisected << "\n";
#endif
   dst.push_back( { a , b , c , p , pnext } );
   haveForm = ! bisected;    // carry the locus; a bisected boundary IS a
   p = pnext;                // form change
   }
  // coalesce ONLY truly-identical adjacent pieces (same quadratic).
  // NB: a VALUE-based merge is UNSOUND here, two pieces that agree in
  // value+slope at their shared kink diverge only QUADRATICALLY over a
  // narrow interval, so different-curvature pieces test "value-close" and
  // merge, keeping the wrong alfa; that wrong curvature then propagates
  // through the horizon. Merge on COEFFICIENTS.
  if( dst.size() > 1 ) {
   PQFun m; m.reserve( dst.size() ); m.push_back( dst.front() );
   for( std::size_t i = 1 ; i < dst.size() ; ++i ) {
    const PieceQuad & bb = m.back() , & cc = dst[ i ];
    const bool same = ( std::abs( bb.alfa - cc.alfa ) <=
                        1e-9 * std::max( 1.0 , std::abs( cc.alfa ) ) )
                   && ( std::abs( bb.beta - cc.beta ) <=
                        1e-7 * std::max( 1.0 , std::abs( cc.beta ) ) )
                   && ( std::abs( bb.gamma - cc.gamma ) <=
                        1e-6 * std::max( 1.0 , std::abs( cc.gamma ) ) );
    if( same ) m.back().right = cc.right; else m.push_back( cc );
    }
   dst.swap( m );
   }
  if( std::getenv( "TUEDPS_FCHK" ) )    // convexity audit of the mpQP VF
   for( std::size_t i = 0 ; i < dst.size() ; ++i ) {
    if( dst[ i ].alfa < -1e-9 )
     std::cerr << "FCHK negcurv a=" << dst[ i ].alfa
               << " [" << dst[ i ].left << "," << dst[ i ].right << "]\n";
    if( i > 0 ) {
     const double bp = dst[ i ].left;
     const double vlo = eval_piece( dst[ i - 1 ] , bp ) ,
                  vhi = eval_piece( dst[ i ] , bp );
     const double slo = 2 * dst[ i - 1 ].alfa * bp + dst[ i - 1 ].beta;
     const double shi = 2 * dst[ i ].alfa * bp + dst[ i ].beta;
     const double sc = std::max( 1.0 , std::abs( vhi ) );
     if( std::abs( vlo - vhi ) > 1e-6 * sc )
      std::cerr << "FCHK discont at " << bp
                << " jump=" << ( vhi - vlo ) << "\n";
     if( shi - slo < -1e-6 * std::max( 1.0 , std::abs( slo ) ) )
      std::cerr << "FCHK slopedrop at " << bp << " sL=" << slo
                << " sR=" << shi << " drop=" << ( slo - shi ) << "\n";
     }
    }
  };
 // PRODUCTION DEFAULT: the principled multiparametric-QP engine (sweepMPQP)
 // builds out. It is the only fast path whose STORED value function is
 // exactly convex (FINCHK=0 over the full day/week/month x u x it grid),
 // where the legacy 6-form sweepParam leaves a materially non-convex F
 // (~820 input-F violations/day-solve). It is also 2.9-5.7x FASTER than
 // sweepParam (day/week/month wall) thanks to the form-carry fast path
 // (carry the finite-diff locus across F-piece boundaries, re-Gmin only at
 // true form changes), and matches golden/MILP (1176/1180 within 1e-4,
 // identical to sweepParam; the residual few are golden=NA or
 // DP-beats-a-gap-limited-MILP). On a finite-diff form-ID failure the sweep
 // never truncates the domain (which would drop every landing power above
 // that point, losing reachable high-power DP states): fillTail covers the
 // whole [plo,phi] with robust adaptive-linear Gmin pieces. Legacy
 // sweepParam is kept gated (TUEDPS_PARAM) for A/B; the seed+subdivide
 // oracle (also exactly convex, but far slower) via TUEDPS_ORACLE.
 if( ! std::getenv( "TUEDPS_ORACLE" ) ) {
  const bool useParam = std::getenv( "TUEDPS_PARAM" );
  if( ! useParam ) sweepMPQP( out );    // DEFAULT: convex + fast + correct
  else sweepParam( out , formreuse );   // legacy 6-form (gated, non-convex)
#if TUEDPS_PROFILE
  // DIRECT construction diff (TUEDPS_MPCMP): here `out` is sweepParam
  // (CORRECT, from correct F). Build mpQP from the SAME F and compare ->
  // isolates mpQP's construction error vs the known-good reference. Dumps
  // the first transition that diverges.
  if( std::getenv( "TUEDPS_MPCMP" ) && ! std::getenv( "TUEDPS_MPQP" ) ) {
   PQFun outM; sweepMPQP( outM );
   double mdm = 0 , mx = plo;
   const int NS = 100;              // vs Gmin (fast, exact pointwise min);
   for( int s = 0 ; s <= NS ; ++s ) {   // `out`=Gmin(p) if reconstruction
    const double x = plo + ( phi - plo ) * s / NS;             // correct
    const double gt = Gmin( x ).first , gm = eval( outM , x );
    if( ( gt < TUEDPINF ) && ( gm < TUEDPINF ) &&
        ( std::abs( gm - gt ) > mdm ) )
     { mdm = std::abs( gm - gt ); mx = x; }
    }
   static double gworst = 0;        // GLOBAL worst over the whole horizon
   const double th = std::getenv( "TUEDPS_MPTH" )
                     ? std::atof( getenv( "TUEDPS_MPTH" ) ) : 1e-3;
   if( ( mdm > gworst + 1e-12 ) && ( mdm > th ) ) {  // dump each NEW
    gworst = mdm; std::cerr.precision( 10 );   // global-worst (last = worst)
    std::cerr << "MPCMP maxdev=" << mdm << " at p=" << mx
              << " |F|=" << F.size()
              << " K=" << K << " ru=" << ramp_up << " rd=" << ramp_down
              << " pmin=" << pmin << " pmax=" << pmax << " domL=" << domL
              << " domR=" << domR << " plo=" << plo << " phi=" << phi
              << "\n F:";
    for( const auto & pc : F ) std::cerr << " [" << pc.left << "," << pc.right
     << "]a=" << pc.alfa << ",b=" << pc.beta << ",c=" << pc.gamma;
    std::cerr << "\n mpQP out:";
    for( const auto & pc : outM )
     std::cerr << "\n   [" << pc.left << "," << pc.right
               << "] a=" << pc.alfa << " b=" << pc.beta << " c=" << pc.gamma;
    std::cerr << "\n samples near maxdev (p | mpQP | Gmin=true | q*_true):";
    for( int s = -6 ; s <= 6 ; ++s ) {
     const double x = mx + ( phi - plo ) * s / 400.0;
     if( x < plo || x > phi ) continue;
     auto gt = Gmin( x );
     std::cerr << "\n   p=" << x << "  M=" << eval( outM , x )
               << "  T=" << gt.first << "  q*=" << gt.second;
     }
    std::cerr << "\n";
    if( std::getenv( "TUEDPS_MPEXIT" ) ) std::exit( 0 );
    }
   }
#endif
#if TUEDPS_PROFILE
  // SELF-CONSISTENT construction probe (TUEDPS_MPSELF): compare the
  // just-built `out` against Gtrue using the SAME (mpQP-propagated) F. A
  // large gap here = a transition CONSTRUCTION bug (not accumulation, which
  // would leave each transition gap ~0).
  if( std::getenv( "TUEDPS_MPSELF" ) ) {
   double mdm = 0 , mx = plo , mgt = 0 , mgx = plo;  // mdm=|out-Gmin|,
   for( const auto & pc : out ) {                    // mgt=|Gmin-Gtrue|;
    // check only at out's breakpoints+mid
    const double xs[ 3 ] = { pc.left , 0.5 * ( pc.left + pc.right ) ,
                             pc.right };
    for( double x : xs ) {
     const double gmin = Gmin( x ).first; if( gmin >= TUEDPINF ) continue;
     const double gm = eval( out , x );
     if( ( gm < TUEDPINF ) && ( std::abs( gm - gmin ) > mdm ) )
      { mdm = std::abs( gm - gmin ); mx = x; }
     const double gtr = Gtrue( x ).first;    // INDEPENDENT ternary: is Gmin
     if( ( gtr < TUEDPINF ) &&               // the TRUE min?
         ( std::abs( gmin - gtr ) > mgt ) )
      { mgt = std::abs( gmin - gtr ); mgx = x; }
     } }
   { static double gwt = 0;
     const double tt = std::getenv( "TUEDPS_MPTH" )
                        ? std::atof( getenv( "TUEDPS_MPTH" ) ) : 1e-4;
     if( ( mgt > gwt + 1e-12 ) && ( mgt > tt ) ) {
      gwt = mgt; std::cerr.precision( 12 );
      auto gg = Gmin( mgx ); auto tg = Gtrue( mgx );
      std::cerr << "GMINBAD |Gmin-Gtrue|=" << mgt << " at p=" << mgx
                << " |F|=" << F.size()
                << " Gmin=" << gg.first << "(q*=" << gg.second
                << ") Gtrue=" << tg.first
                << "(q*=" << tg.second << ") plo=" << plo << " phi=" << phi
                << " acap-note pmax=" << pmax << "\n";
      if( std::getenv( "TUEDPS_MPEXIT" ) ) std::exit( 0 ); } }
   { static double gwg = 0;
     const double th2 = std::getenv( "TUEDPS_MPTH" )
                         ? std::atof( getenv( "TUEDPS_MPTH" ) ) : 1e-3;
     if( ( mgt > gwg + 1e-12 ) && ( mgt > th2 ) ) {
      gwg = mgt; std::cerr.precision( 10 );
      std::cerr << "GMINBAD |Gmin-Gtrue|=" << mgt << " at p=" << mgx
                << " |F|=" << F.size()
                << " Gmin=" << Gmin( mgx ).first
                << " Gtrue=" << Gtrue( mgx ).first
                << " q*mpq=... acap-note plo=" << plo
                << " phi=" << phi << "\n"; } }
   static double gworst = 0;        // GLOBAL worst over the whole horizon
   const double th = std::getenv( "TUEDPS_MPTH" )
                      ? std::atof( getenv( "TUEDPS_MPTH" ) ) : 1e-3;
   if( ( mdm > gworst + 1e-12 ) && ( mdm > th ) ) {  // dump each new
    gworst = mdm; std::cerr.precision( 10 );   // global-worst (last = worst)
    std::cerr << "MPSELF maxdev=" << mdm << " at p=" << mx
              << " |F|=" << F.size()
              << " K=" << K << " ru=" << ramp_up << " rd=" << ramp_down
              << " pmin=" << pmin << " pmax=" << pmax << " domL=" << domL
              << " domR=" << domR << " plo=" << plo << " phi=" << phi
              << "\n F:";
    for( const auto & pc : F ) std::cerr << " [" << pc.left << "," << pc.right
     << "]a=" << pc.alfa << ",b=" << pc.beta << ",c=" << pc.gamma;
    std::cerr << "\n out pieces:";
    for( const auto & pc : out )
     std::cerr << "\n   [" << pc.left << "," << pc.right
               << "] a=" << pc.alfa << " b=" << pc.beta << " c=" << pc.gamma;
    std::cerr
     << "\n samples near maxdev (p | out | Gmin=true | q*_true | gap):";
    for( int s = -8 ; s <= 8 ; ++s ) {
     const double x = mx + ( phi - plo ) * s / 300.0;
     if( x < plo || x > phi ) continue;
     auto gt = Gmin( x );
     std::cerr << "\n   p=" << x << "  O=" << eval( out , x )
               << "  T=" << gt.first << "  q*=" << gt.second
               << "  gap=" << ( eval( out , x ) - gt.first );
     }
    std::cerr << "\n";
    if( std::getenv( "TUEDPS_MPEXIT" ) ) std::exit( 0 );
    }
   }
#endif
  if( std::getenv( "TUEDPS_FRDIFF" ) ) {  // A/B: reuse vs classify-every
   PQFun a2 , b2; sweepParam( a2 , false ); sweepParam( b2 , true );
   double md = 0 , pw = 0;
   for( int k = 0 ; k <= 200 ; ++k ) {    // sample G(p) from both, max |diff|
    const double pp = plo + ( phi - plo ) * ( k / 200.0 );
    double va = TUEDPINF , vb = TUEDPINF;
    for( const auto & q : a2 )
     if( pp >= q.left - 1e-9 && pp <= q.right + 1e-9 )
      { va = eval_piece( q , pp ); break; }
    for( const auto & q : b2 )
     if( pp >= q.left - 1e-9 && pp <= q.right + 1e-9 )
      { vb = eval_piece( q , pp ); break; }
    if( ( va < TUEDPINF ) && ( vb < TUEDPINF ) ) {
     const double e = std::abs( va - vb );
     if( e > md ) { md = e; pw = pp; } }
    }
   const double sc = std::max( 1.0 ,
    std::abs( eval_piece( a2.empty() ? out.front() : a2.front() , plo ) ) );
   if( md > 1e-7 * sc )
    std::cerr << "FRDIFF t=" << t << " md=" << md << " at p=" << pw
              << " napcs=" << a2.size() << " nbpcs=" << b2.size() << "\n";
   }
  // as at the end of sliding_min(): a domain collapsed to a single point is
  // a value, not the absence of one
  if( out.empty() && ( plo <= phi + 1e-12 ) ) {
   const double v = Geval( plo );
   if( v < TUEDPINF )
    out.push_back( { 0 , 0 , v , plo , phi } );
   }
#if TUEDPS_PROFILE
  g_pieces += out.size();                 // (g_smc already counted above)
#endif
  return;
  }

 static thread_local std::vector< double > ps;
 ps.clear();
 ps.push_back( plo ); ps.push_back( phi );
 auto addp = [ & ]( double x ) {
  if( ( x > plo + 1e-12 ) && ( x < phi - 1e-12 ) ) ps.push_back( x );
  };
 for( const auto & pc : F ) {
  addp( pc.left  + wu ); addp( pc.left  - wd );
  addp( pc.right + wu ); addp( pc.right - wd );
  }
 addp( pmin ); addp( 0.5 * ( pmin + pmax ) ); addp( pmax );
 if( rho1 < 1 ) addp( pmin / ( 1 - rho1 ) );
 addp( pmax / ( 1 + rho1 ) );
 if( rho2 > 0 ) {
  if( rhot < 1 ) addp( pmin / ( 1 - rhot ) );
  addp( pmax / ( 1 + rhot ) );
  }
 std::sort( ps.begin() , ps.end() );
 ps.erase( std::unique( ps.begin() , ps.end() ,
            []( double a , double b ) { return( b - a <= 1e-12 ); } ) ,
           ps.end() );

 const double ftol = 1e-6;
 // Values of G at the seed points, computed once and shared by adjacent
 // seed intervals; each subdivision then carries its endpoint values down
 // to its children (a child endpoint is either a parent endpoint or the
 // parent midpoint, all already known), so every node costs one Gmin at its
 // midpoint plus two Geval verify probes, not five Gmin.
 static thread_local std::vector< double > vps;
 vps.resize( ps.size() );
 for( std::size_t i = 0 ; i < ps.size() ; ++i ) vps[ i ] = Geval( ps[ i ] );

 // subdivide each seed interval by an EXPLICIT stack (no recursive
 // std::function, whose type-erased call is costly); pushing the right
 // half before the left keeps the emitted pieces in increasing-p order
 struct Iv { double pl , pr , vl , vr; int depth; };
 static thread_local std::vector< Iv > stk;
 for( std::size_t i = 0 ; i + 1 < ps.size() ; ++i ) {
  if( ps[ i + 1 ] - ps[ i ] <= 1e-12 ) continue;
  stk.clear();
  stk.push_back( { ps[ i ] , ps[ i + 1 ] , vps[ i ] , vps[ i + 1 ] , 0 } );
  while( ! stk.empty() ) {
   const Iv iv = stk.back(); stk.pop_back();
#if TUEDPS_PROFILE
   ++g_node;
#endif
   const double pl = iv.pl , pr = iv.pr , vl = iv.vl , vr = iv.vr;
   if( ( vl >= TUEDPINF ) || ( vr >= TUEDPINF ) )  // outside finite domain
    continue;
   const double pm = 0.5 * ( pl + pr );
   auto [ vm , qm ] = Gmin( pm );
   const double scale = std::max( 1.0 , std::abs( vm ) );
   // classify the minimiser at the midpoint. Only a WINDOW-edge minimiser
   // (q = p-ramp_up or p+ramp_down) makes G quadratic, with curvature
   // EXACTLY the alfa of the F-piece at that q -- anchored, not fitted, so
   // no spurious curvature can be manufactured. An interior or
   // F-domain-edge minimiser is constant in p on the piece, so G is linear.
   double qedge = 0; bool windowedge = false;
   if( std::abs( qm - ( pm - wu ) ) <=
       1e-7 * std::max( 1.0 , std::abs( pm ) ) )
    { qedge = pm - wu; windowedge = true; }
   else if( std::abs( qm - ( pm + wd ) ) <=
            1e-7 * std::max( 1.0 , std::abs( pm ) ) )
    { qedge = pm + wd; windowedge = true; }
   double a , b , c;
   if( windowedge ) {             // G(p) = F(p∓Δ) + (linear reserve term)
    const double aF = F[ findF( qedge ) ].alfa;
    a = aF;                       // curvature bounded by F's own curvature
    const double rl = vl - aF * pl * pl , rr = vr - aF * pr * pr;
    b = ( rr - rl ) / ( pr - pl );
    c = rl - b * pl;
    }
   else {                         // interior / domain-edge => G is linear
    a = 0;
    b = ( vr - vl ) / ( pr - pl );
    c = vl - b * pl;
    }
   const double t1 = 0.5 * ( pl + pm ) , t2 = 0.5 * ( pm + pr );
   const double e1 = std::abs( a * t1 * t1 + b * t1 + c - Geval( t1 ) );
   const double e2 = std::abs( a * t2 * t2 + b * t2 + c - Geval( t2 ) );
   if( ( ( e1 <= ftol * scale ) && ( e2 <= ftol * scale ) ) ||
       ( pr - pl <= 1e-9 ) ) {
    out.push_back( { a , b , c , pl , pr } );
    }
   else if( iv.depth >= 14 ) {   // safety: emit the straight segment
    const double bb = ( vr - vl ) / ( pr - pl );        // (bounded)
    out.push_back( { 0.0 , bb , vl - bb * pl , pl , pr } );
    }
   else {                        // push right then left => left first
    stk.push_back( { pm , pr , vm , vr , iv.depth + 1 } );
    stk.push_back( { pl , pm , vl , vm , iv.depth + 1 } );
    }
   }
  }

 // coalesce adjacent pieces that carry (nearly) the same quadratic: a
 // genuinely quadratic region of G, subdivided by the refinement into
 // several identical pieces, collapses back to one, keeping the piece count
 // bounded (essential: these functions propagate through the whole DP)
 if( out.size() > 1 ) {
  PQFun m;
  m.reserve( out.size() );
  m.push_back( out.front() );
  for( std::size_t i = 1 ; i < out.size() ; ++i ) {
   const PieceQuad & b = m.back();
   const PieceQuad & c = out[ i ];
   // value-based test (robust to fp noise in the fitted coefficients): the
   // two pieces are the same quadratic iff b, extended over c's interval,
   // matches c
   const double x1 = c.right , x2 = 0.5 * ( c.left + c.right );
   const double e1 = std::abs( eval_piece( b , x1 ) - eval_piece( c , x1 ) );
   const double e2 = std::abs( eval_piece( b , x2 ) - eval_piece( c , x2 ) );
   const double sc = std::max( 1.0 , std::abs( eval_piece( c , x2 ) ) );
   const double ctol = std::getenv( "TUEDPS_CTOL" )
                       ? std::atof( std::getenv( "TUEDPS_CTOL" ) ) : 1e-6;
   if( ( e1 <= ctol * sc ) && ( e2 <= ctol * sc ) )
    m.back().right = c.right;       // same quadratic: extend b over c
   else
    m.push_back( c );
   }
  out.swap( m );
  }

 // as in sliding_min(): a domain collapsed to a single point is a value,
 // not the absence of one, and must be emitted as a zero-width piece
 if( out.empty() && ( plo <= phi + 1e-12 ) ) {
  const double v = Geval( plo );
  if( v < TUEDPINF )
   out.push_back( { 0 , 0 , v , plo , phi } );
  }

#if TUEDPS_PROFILE
 // ------------------------------------------------------------------------
 // REGIME MAP (env TUEDPS_REGIME): ground-truth breakpoint discovery for
 // the exact parametric. Sweeps p, gets q*(p) from the trustworthy Gtrue,
 // and brackets every point where the REGIME SIGNATURE changes --
 // (minimiser class, F-piece at q*, tent side, active reward segment of B
 // and of A) -- then bisects to the exact breakpoint (~1e-12). Each smooth
 // segment between breakpoints is one analytic G-piece; this list is the
 // reference the event algebra must reproduce with NO verify. Dumps the
 // call with |F|==TUEDPS_REGIMEF (default 1), the TUEDPS_REGIMEN-th such
 // occurrence (default 1).
 if( std::getenv( "TUEDPS_REGIME" ) ) {
  static int r_seen = 0;
  const std::size_t fwant = std::getenv( "TUEDPS_REGIMEF" )
   ? ( std::size_t ) std::atoi( std::getenv( "TUEDPS_REGIMEF" ) ) : 1;
  const int nwant = std::getenv( "TUEDPS_REGIMEN" )
                    ? std::atoi( std::getenv( "TUEDPS_REGIMEN" ) ) : 1;
  if( ( F.size() == fwant ) && ( ++r_seen == nwant ) ) {
   double segcum2[ 2 ] = { 0 , 0 }; int Kr = 0;
   if( ( cp < 0 ) && ( cs < 0 ) ) {
    if( cp <= cs ) { segcum2[ 0 ] = rp; segcum2[ 1 ] = rp + rs; }
    else           { segcum2[ 0 ] = rs; segcum2[ 1 ] = rs + rp; } Kr = 2; }
   else if( cp < 0 ) { segcum2[ 0 ] = rp; Kr = 1; }
   else if( cs < 0 ) { segcum2[ 0 ] = rs; Kr = 1; }
   const double dpk = 0.5 * ( ramp_up - ramp_down );
   // The signature must be STABLE to Gtrue's q*-localisation noise (~1e-9):
   // raw findF(q*) flips between the two pieces adjacent to a kink when q*
   // is pinned there, manufacturing hundreds of false breaks. So classify
   // by regime FORM (how q* depends on p) with tight tolerances on the
   // moving laws, and encode a PINNED q* by its nearest-F-boundary index
   // (constant), not by findF.
   struct Sig { int form , key , side , segB , segA; bool fin; };
   // 0=PIN(fixed q*)  1=WUP q=p-ru  2=WDN q=p+rd  3=TENT q=p-dpk  4=BAND B=A
   // 5=STAT interior stationary (q* constant, not at a boundary)
   auto boundIdx = [ & ]( double q ) -> int {   // nearest F-boundary id
    int bi = 0; double best = std::abs( q - F[ 0 ].left );
    for( std::size_t k = 0 ; k < F.size() ; ++k ) {
     const double dd = std::abs( q - F[ k ].right );
     if( dd < best ) { best = dd; bi = int( k ) + 1; } }
    return bi; };
   auto sig = [ & ]( double p ) -> Sig {
    auto g = Gtrue( p );
    if( g.first >= TUEDPINF ) return Sig{ 0 , 0 , 0 , 0 , 0 , false };
    const double q = g.second , tol = 1e-7;
    const double d = p - q;
    double B = std::min( ramp_up - d , ramp_down + d );
    if( B < 0 ) B = 0; const double A = Aband( p );
    int form , key;
    // does q* sit on a known F breakpoint (pinned)? distance to nearest
    // bound:
    const int bi = boundIdx( q );
    const double qb = ( bi == 0 ) ? F[ 0 ].left : F[ bi - 1 ].right;
    if( std::abs( q - ( p - ramp_up ) ) < tol )
     { form = 1; key = findF( p - ramp_up ); }
    else if( std::abs( q - ( p + ramp_down ) ) < tol )
     { form = 2; key = findF( p + ramp_down ); }
    else if( std::abs( q - ( p - dpk ) ) < tol )
     { form = 3; key = findF( p - dpk ); }
    else if( std::abs( q - qb ) < tol )
     { form = 0; key = bi; }                     // pinned at boundary
    else if( std::abs( B - A ) < tol * std::max( 1.0 , std::abs( A ) ) )
     { form = 4; key = findF( q ); }
    else { form = 5; key = findF( q ); }         // interior stationary
    auto segH = [ & ]( double H ) -> int {
     int m = 0;
     for( ; m < Kr ; ++m ) if( H <= segcum2[ m ] * p + 1e-9 ) break;
     return m; };
    const int side = ( std::abs( d - dpk ) < 1e-6 ) ? 2 : ( d > dpk ? 1 : 0 );
    return Sig{ form , key , side , segH( std::min( A , B ) ) , segH( A ) ,
                true };
    };
   auto same = []( const Sig & a , const Sig & b ) {
    return( a.fin == b.fin && a.form == b.form && a.key == b.key &&
            a.side == b.side && a.segB == b.segB && a.segA == b.segA ); };
   const char * CLS[] = { "PIN" , "WUP" , "WDN" , "TENT" , "BAND" , "STAT" };
   auto row = [ & ]( double p , const char * tag ) {
    auto g = Gtrue( p ); const double q = g.second , d = p - q;
    double B = std::min( ramp_up - d , ramp_down + d ); if( B < 0 ) B = 0;
    Sig s = sig( p );
    std::cerr << "  " << tag << " p=" << p << " G=" << g.first << " q*=" << q
              << " d=" << d << " B=" << B << " A=" << Aband( p )
              << " k1=" << rho1 * p << " kt=" << rhot * p
              << " form=" << ( s.fin ? CLS[ s.form ] : "INF" )
              << " key=" << s.key << " side=" << s.side
              << " segB=" << s.segB << " segA=" << s.segA << "\n"; };
   std::cerr.precision( 12 );
   std::cerr << "REGIME plo=" << plo << " phi=" << phi << " |F|=" << F.size()
             << " K=" << Kr << " ru=" << ramp_up << " rd=" << ramp_down
             << " pmin=" << pmin << " pmax=" << pmax << " rho1=" << rho1
             << " rhot=" << rhot << " dpk=" << dpk
             << " domL=" << domL << " domR=" << domR << "\n";
   std::cerr << " F pieces:"; for( const auto & pc : F )
    std::cerr << " [" << pc.left << "," << pc.right << "|a=" << pc.alfa
              << ",b=" << pc.beta << ",c=" << pc.gamma << "]";
   std::cerr << "\n";
   const int NSW = std::getenv( "TUEDPS_REGIMENSW" )
                   ? std::atoi( getenv( "TUEDPS_REGIMENSW" ) ) : 8000;
   row( plo , "START" );
   double pprev = plo; Sig sprev = sig( plo ); int nbr = 0;
   for( int s = 1 ; s <= NSW ; ++s ) {
    const double p = plo + ( phi - plo ) * s / NSW;
    const Sig sc = sig( p );
    if( ! same( sc , sprev ) ) {
     double a = pprev , b = p;
     for( int it = 0 ;
          it < 90 && b - a > 1e-13 * std::max( 1.0 , std::abs( b ) ) ;
          ++it ) {
      const double m = 0.5 * ( a + b );
      if( same( sig( m ) , sprev ) ) a = m; else b = m;
      }
     row( 0.5 * ( a + b ) , "BREAK" ); ++nbr;
     sprev = sc;
     }
    pprev = p;
    }
   row( phi , "END" );
   std::cerr << "REGIME nbreak=" << nbr
             << " oracle_pcs=" << out.size() << "\n";
   }
  }

 // ------------------------------------------------------------------------
 // PARAMETRIC CONTINUATION (env TUEDPS_PARAM): build the same G by sweeping
 // p and predicting the next regime breakpoint analytically (a piece per
 // true breakpoint, ~1 Gmin/piece), instead of seed+verify+subdivide. The
 // current `out` above is the validated ORACLE; here `outP` is built and
 // diffed, so a mistake in the (delicate) event algebra is caught, never
 // shipped. Generic in the number K of active reward segments (1 or 2
 // reserves): the events loop over the segment list, so single-reserve just
 // has fewer of them.
 if( std::getenv( "TUEDPS_PARAM" ) ) {
  PQFun outP;
  sweepParam( outP , false );          // the shared production engine
  // mpQP-vs-Gtrue construction probe (env TUEDPS_MPDUMP): here F is the
  // CORRECT oracle-propagated input, so any outM-vs-Gtrue gap is a mpQP
  // CONSTRUCTION bug (isolated from horizon accumulation). Dumps the FIRST
  // transition that diverges.
  if( std::getenv( "TUEDPS_MPDUMP" ) ) {
   PQFun outM; sweepMPQP( outM );
   double mdm = 0 , mx = plo;
   const int NS = 2000;        // vs Gmin (fast) on the CORRECT oracle F
   for( int s = 0 ; s <= NS ; ++s ) {
    const double x = plo + ( phi - plo ) * s / NS;
    const double gt = Gmin( x ).first; if( gt >= TUEDPINF ) continue;
    const double gm = eval( outM , x );
    if( ( gm < TUEDPINF ) && ( std::abs( gm - gt ) > mdm ) )
     { mdm = std::abs( gm - gt ); mx = x; }
    }
   static bool mdumped = false;
   const double mth = std::getenv( "TUEDPS_MPTH" )
                       ? std::atof( getenv( "TUEDPS_MPTH" ) ) : 1e-3;
   if( ( mdm > mth ) && ! mdumped ) {
    mdumped = true; std::cerr.precision( 10 );
    std::cerr << "MPDUMP maxdev=" << mdm << " at p=" << mx
              << " |F|=" << F.size()
              << " K=" << K << " ru=" << ramp_up << " rd=" << ramp_down
              << " pmin=" << pmin << " pmax=" << pmax << " domL=" << domL
              << " domR=" << domR << " rho1=" << rho1 << " rhot=" << rhot
              << " plo=" << plo << " phi=" << phi << "\n F:";
    for( const auto & pc : F ) std::cerr << " [" << pc.left << "," << pc.right
     << "]a=" << pc.alfa << ",b=" << pc.beta << ",c=" << pc.gamma;
    std::cerr << "\n outM pieces:";
    for( const auto & pc : outM )
     std::cerr << "\n   [" << pc.left << "," << pc.right
               << "] a=" << pc.alfa << " b=" << pc.beta << " c=" << pc.gamma;
    std::cerr << "\n samples (p | mpQP | Gtrue | q*_true):";
    for( int s = 0 ; s <= 40 ; ++s ) {
     const double x = plo + ( phi - plo ) * s / 40;
     auto gt = Gtrue( x );
     std::cerr << "\n   p=" << x << "  M=" << eval( outM , x )
               << "  T=" << gt.first << "  q*=" << gt.second
               << "  gap=" << ( eval( outM , x ) - gt.first );
     }
    std::cerr << "\n";
    if( std::getenv( "TUEDPS_MPEXIT" ) )
     std::exit( 0 );                   // stop after first dump (fast)
    }
   }
  static thread_local std::vector< double > swtrace;
  swtrace.clear();                                          // (debug slot)
  // grid probe of the (q*, B, A, kappa) config, to derive the exact events
  static bool gridded = false;
  if( ! gridded && ( F.size() == 1 ) && std::getenv( "TUEDPS_PARAMGRID" ) ) {
   gridded = true;
   std::cerr.precision( 8 );
   std::cerr << "PARAMGRID plo=" << plo << " phi=" << phi
             << " ru=" << ramp_up << " rd=" << ramp_down
             << " pmin=" << pmin << " pmax=" << pmax
             << " rho1=" << rho1 << "\n";
   const int NG = 120;
   for( int s = 0 ; s <= NG ; ++s ) {
    const double p = plo + ( phi - plo ) * s / NG;
    auto g = Gmin( p );
    const double q = g.second , d = p - q;
    double B = std::min( ramp_up - d , ramp_down + d ); if( B < 0 ) B = 0;
    const double A = std::min( p - pmin , pmax - p );
    const double k1 = rho1 * p;
    const char * cls = ( std::abs( q - ( p - ramp_up ) ) < 1e-6 ) ? "WUP"
                     : ( std::abs( q - ( p + ramp_down ) ) < 1e-6 ) ? "WDN"
                     : ( std::abs( q - domL ) < 1e-6 ) ? "DL"
                     : ( std::abs( q - domR ) < 1e-6 ) ? "DR"
                     : ( std::abs( B - A ) < 1e-4 ) ? "BAND"
                     : ( B >= A ) ? "PEN0" : "WING";
    std::cerr << "  p=" << p << " q*=" << q << " B=" << B << " A=" << A
              << " k1=" << k1 << " G=" << g.first << " " << cls << "\n";
    }
   }

  // diff against the TRUE G (Gmin) at dense samples -- NOT against the
  // oracle, whose breakpoint LOCATIONS are only accurate to the subdivision
  // resolution (its values are exact to ftol, but a slope change is placed
  // at a bisection midpoint, not the true breakpoint). A verify-free
  // parametric with exact events can be MORE accurate than the oracle, so
  // the oracle is the wrong ruler.
  double maxdev = 0 , maxdev_or = 0 , maxdev_gm = 0;
  if( std::getenv( "TUEDPS_PARAMDIFF" ) ) {  // expensive Gtrue-based diff,
   const int NS = 100;                       // opt-in
   for( int s = 0 ; s <= NS ; ++s ) {
    const double x = plo + ( phi - plo ) * s / NS;
    const double gt = Gtrue( x ).first;   // trustworthy reference (not Gmin)
    if( gt >= TUEDPINF ) continue;
    const double gp = eval( outP , x ) , go = eval( out , x );
    const double gm = Gmin( x ).first;    // Gmin's own error vs the true G
    if( gp < TUEDPINF ) maxdev = std::max( maxdev , std::abs( gp - gt ) );
    if( go < TUEDPINF )
     maxdev_or = std::max( maxdev_or , std::abs( go - gt ) );
    if( gm < TUEDPINF )
     maxdev_gm = std::max( maxdev_gm , std::abs( gm - gt ) );
    // pinpoint WHERE Gmin fails: on each new global-worst gm-vs-Gtrue gap,
    // dump the two minimisers and the local config so the missing candidate
    // is visible
    if( ( gm < TUEDPINF ) && ( std::abs( gm - gt ) > g_gmin_maxdev + 1e-15 )
        && ( std::abs( gm - gt ) > 1e-6 ) &&
        std::getenv( "TUEDPS_GMINDUMP" ) ) {
     auto [ vg , qg ] = Gmin( x );
     auto [ vtt , qt ] = Gtrue( x );
     const double dg = x - qg , dt = x - qt;
     auto Bof = [ & ]( double d ){
      double B = std::min( ramp_up - d , ramp_down + d );
      return( B<0 ? 0 : B ); };
     std::cerr.precision( 12 );
     std::cerr << "GMINDUMP p=" << x << " |F|=" << F.size()
               << " gap=" << ( vg - vtt )
               << " | Gmin q*=" << qg << " v=" << vg
               << " (d=" << dg << " B=" << Bof( dg ) << ")"
               << "  [eF=" << evalF( qg ) << " corr=" << corr( qg , x )
               << " sum=" << ( evalF( qg ) + corr( qg , x ) ) << "]"
               << " | Gtrue q*=" << qt << " v=" << vtt
               << " (d=" << dt << " B=" << Bof( dt ) << ")"
               << "  [eF=" << evalF( qt ) << " corr=" << corr( qt , x )
               << " sum=" << ( evalF( qt ) + corr( qt , x ) ) << "]"
               << " | A=" << Aband( x ) << " k1=" << rho1 * x
               << " kt=" << rhot * x
               << " ru=" << ramp_up << " rd=" << ramp_down
               << " domL=" << domL << " domR=" << domR
               << " findF(qg)=" << findF( qg )
               << " findF(qt)=" << findF( qt ) << "\n";
     }
    }
   if( maxdev_or > g_param_ormaxdev ) g_param_ormaxdev = maxdev_or;
   if( maxdev_gm > g_gmin_maxdev )    g_gmin_maxdev = maxdev_gm;
   }
  g_param_calls++;
  if( maxdev > g_param_maxdev ) g_param_maxdev = maxdev;
  if( maxdev > 1e-4 ) g_param_bad++;
  g_param_pcs += outP.size();
  static bool dumped = false;
  const std::size_t dumpF = std::getenv( "TUEDPS_DUMPF" )
   ? ( std::size_t ) std::atoi( std::getenv( "TUEDPS_DUMPF" ) ) : 1;
  const double dumpTh = std::getenv( "TUEDPS_DUMPTH" )
                        ? std::atof( std::getenv( "TUEDPS_DUMPTH" ) ) : 1e-4;
  if( ( maxdev > dumpTh ) && ! dumped && ( F.size() <= dumpF )
      && std::getenv( "TUEDPS_PARAMDUMP" ) ) {
   dumped = true;
   std::cerr.precision( 10 );
   std::cerr << "PARAMDUMP plo=" << plo << " phi=" << phi
             << " maxdev=" << maxdev
             << " |F|=" << F.size() << " K=" << K
             << " ramp_up=" << ramp_up << " ramp_down=" << ramp_down
             << " pmin=" << pmin << " pmax=" << pmax
             << " rho1=" << rho1 << " rhot=" << rhot << "\n";
   std::cerr << " F: "; for( const auto & pc : F )
    std::cerr << "[" << pc.left << "," << pc.right << "] a=" << pc.alfa
              << " b=" << pc.beta << " c=" << pc.gamma << "  ";
   std::cerr << "\n ORACLE breakpoints:";
   for( const auto & pc : out ) std::cerr << " " << pc.left;
   std::cerr << " " << out.back().right << "\n PARAM breakpoints:";
   for( const auto & pc : outP ) std::cerr << " " << pc.left;
   std::cerr << " " << outP.back().right << "\n outP pieces:";
   for( const auto & pc : outP )
    std::cerr << "\n   [" << pc.left << "," << pc.right << "] a=" << pc.alfa
              << " b=" << pc.beta << " c=" << pc.gamma;
   std::cerr << "\n sweep steps (p form qbar/kappa pnext):";
   for( std::size_t i = 0 ; i + 3 < swtrace.size() ; i += 4 )
    std::cerr << "\n   p=" << swtrace[ i ]
              << " form=" << ( int )swtrace[ i + 1 ]
              << " q/k=" << swtrace[ i + 2 ] << " pnext=" << swtrace[ i + 3 ];
   std::cerr << "\n samples (p | param | Gtrue | q* | B | A):";
   const double slo = std::getenv( "TUEDPS_SLO" )
                       ? std::atof( getenv( "TUEDPS_SLO" ) ) : plo;
   const double shi = std::getenv( "TUEDPS_SHI" )
                       ? std::atof( getenv( "TUEDPS_SHI" ) ) : phi;
   for( int s = 0 ; s <= 20 ; ++s ) {
    const double x = slo + ( shi - slo ) * s / 20;
    auto gt = Gtrue( x );
    const double d = x - gt.second;
    double B = std::min( ramp_up - d , ramp_down + d );
    if( B<0 )B = 0; const double A = std::min( x - pmin , pmax - x );
    std::cerr << "\n   p=" << x << "  P=" << eval( outP , x )
              << "  T=" << gt.first << "  q*=" << gt.second
              << "  B=" << B << "  A=" << A;
    }
   std::cerr << "\n";
   }
  }
#endif

#if TUEDPS_PROFILE
 g_pieces += out.size();
 if( corr_hit ) ++g_smc_bite;
#endif

 }  // end( ThermalUnitDPSolverBase::sliding_min_corr )

/*--------------------------------------------------------------------------*/

double ThermalUnitDPSolverBase::reserve_corr_argmin(
 const PQFun & F , double ramp_up , double ramp_down , Index t , double p ,
 double acap , double win_up , double win_down ) const
{
 if( F.empty() ) return( p );
 // window of the scheduled move, the tent being that of ramp_up / ramp_down
 const double wu = std::isnan( win_up ) ? ramp_up : win_up;
 const double wd = std::isnan( win_down ) ? ramp_down : win_down;
 const double domL = F.front().left , domR = F.back().right;
 const double qL = std::max( p - wu , domL );
 const double qR = std::min( p + wd , domR );
 if( qR - qL <= 1e-12 ) return( std::min( std::max( p - wu , domL ) ,
                                          domR ) );

 const double cp = primary_reserve_cost.empty()   ? 0
                                                  : primary_reserve_cost[ t ];
 const double cs = secondary_reserve_cost.empty()
                   ? 0 : secondary_reserve_cost[ t ];

 // energy argmin (projection of argmin F onto the window) -- also the answer
 // when no reserve is rewarded
 double pstar = 0 , fstar = TUEDPINF;
 for( const auto & pc : F ) {
  const double q = argmin_piece( pc ) , v = eval_piece( pc , q );
  if( v < fstar ) { fstar = v; pstar = q; }
  }
 double qeng = pstar;
 if( p > pstar + wu )        qeng = p - wu;
 else if( p < pstar - wd ) qeng = p + wd;
 if( qeng < qL ) qeng = qL;
 if( qeng > qR ) qeng = qR;
 if( ( cp >= 0 ) && ( cs >= 0 ) )
  return( qeng );

 const double pmin = min_power[ t ];
 const double pmax = ( acap < 0.0 ) ? max_power[ t ] : acap;  // shut-down cap
 const double rp = primary_rho.empty()   ? 0 : primary_rho[ t ];
 const double rs = secondary_rho.empty() ? 0 : secondary_rho[ t ];
 double rho1 = 0 , rho2 = 0;
 if( ( cp < 0 ) && ( cs < 0 ) ) {
  if( cp <= cs ) { rho1 = rp; rho2 = rs; } else { rho1 = rs; rho2 = rp; }
  }
 else if( cp < 0 ) rho1 = rp; else rho1 = rs;
 const double rhot = rho1 + rho2;

 const double A = std::min( p - pmin , pmax - p );
 auto corr = [ & ]( double q ) -> double {
  if( A <= 0 ) return( 0 );
  const double d = p - q;
  double B = std::min( ramp_up - d , ramp_down + d );
  if( B < 0 ) B = 0;
  if( B >= A ) return( 0 );
  double pr , sr;
  return( reserve_alloc_band( t , p , B , pr , sr ) -
          reserve_alloc_band( t , p , A , pr , sr ) );
  };

 static thread_local std::vector< double > qs;
 qs.clear();
 qs.push_back( qL ); qs.push_back( qR );
 for( const auto & pc : F )
  if( ( pc.left > qL + 1e-12 ) && ( pc.left < qR - 1e-12 ) )
   qs.push_back( pc.left );
 auto add_v = [ & ]( double v ) {
  if( v <= 1e-12 ) return;
  const double q1 = p - ramp_up + v , q2 = p + ramp_down - v;
  if( ( q1 > qL + 1e-12 ) && ( q1 < qR - 1e-12 ) ) qs.push_back( q1 );
  if( ( q2 > qL + 1e-12 ) && ( q2 < qR - 1e-12 ) ) qs.push_back( q2 );
  };
 if( A > 0 ) { add_v( A ); add_v( rho1 * p );
               if( rho2 > 0 ) add_v( rhot * p ); }
 { const double qpk = p - 0.5 * ( ramp_up - ramp_down );
   if( ( qpk > qL + 1e-12 ) && ( qpk < qR - 1e-12 ) ) qs.push_back( qpk ); }
 std::sort( qs.begin() , qs.end() );
 double best = TUEDPINF , bestq = qeng;
 auto upd = [ & ]( double q ) {
  const double v = eval( F , q ) + corr( q );
  if( v < best ) { best = v; bestq = q; } };
 for( std::size_t i = 0 ; i + 1 < qs.size() ; ++i ) {
  const double a = qs[ i ] , b = qs[ i + 1 ];
  if( b - a <= 1e-12 ) continue;
  upd( a ); upd( b );
  const double mid = 0.5 * ( a + b );
  double aF = 0 , bF = 0;
  for( const auto & pc : F )
   if( ( mid >= pc.left - 1e-12 ) && ( mid <= pc.right + 1e-12 ) ) {
    aF = pc.alfa; bF = pc.beta; break;
    }
  const double sc = ( corr( b ) - corr( a ) ) / ( b - a );
  if( aF > 1e-16 ) {
   const double qv = - ( bF + sc ) / ( 2.0 * aF );
   if( qv > a && qv < b ) upd( qv );
   }
  }
 return( bestq );

 }  // end( ThermalUnitDPSolverBase::reserve_corr_argmin )

/*--------------------------------------------------------------------------*/
/*--------------- End File ThermalUnitDPSolverBase.cpp ---------------------*/
/*--------------------------------------------------------------------------*/
