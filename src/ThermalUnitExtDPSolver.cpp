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

void ThermalUnitExtDPSolver::get_var_solution( Configuration * solc )
{
 // lock the Block
 bool owned = f_Block->is_owned_by( f_id );
 if( ( ! owned ) && ( ! f_Block->lock( f_id ) ) )
  throw( std::runtime_error(
   "ThermalUnitExtDPSolver::get_var_solution: unable to lock the Block." ) );

 auto b = static_cast< ThermalUnitBlock * >( f_Block );

 if( auto pow_it = b->get_active_power( 0 ) )
  for( Index i = 0 ; i < time_horizon ; )
   ( pow_it++ )->set_value( P[ i++ ] );

 if( auto com_it = b->get_commitment( 0 ) )
  for( Index i = 0 ; i < time_horizon ; )
   ( com_it++ )->set_value( U[ i++ ] ? 1 : 0 );

 if( auto sup_it = b->get_start_up() ) {
  if( ! t_init )
   ( sup_it++ )->set_value( ( init_up_down_time <= 0 ) && ( U[ 0 ] ? 1 : 0 ) );
  for( Index i = std::max( t_init , Index( 1 ) ) ; i < time_horizon ; ++i )
   ( sup_it++ )->set_value( ( U[ i ] ) && ( ! U[ i - 1 ] ) ? 1 : 0 );
  }

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

void ThermalUnitExtDPSolver::load_parameters( void )
{
 bool owned = f_Block->is_owned_by( f_id );
 if( ( ! owned ) && ( ! f_Block->read_lock() ) )
  throw( std::runtime_error(
   "ThermalUnitExtDPSolver::load_parameters: unable to lock the Block." ) );

 auto b = static_cast< ThermalUnitBlock * >( f_Block );

 if( ! b->get_primary_rho().empty() )
  throw( std::invalid_argument( "ThermalUnitExtDPSolver::load_parameters: "
                                "primary reserve not supported." ) );

 if( ! b->get_secondary_rho().empty() )
  throw( std::invalid_argument( "ThermalUnitExtDPSolver::load_parameters: "
                                "secondary reserve not supported." ) );

 time_horizon      = b->get_time_horizon();
 init_up_down_time = b->get_init_up_down_time();
 min_up_time       = b->get_min_up_time();
 min_down_time     = b->get_min_down_time();
 initial_power     = b->get_initial_power();

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

 startup_costs = b->get_start_up_cost();
 min_power     = b->get_min_power();
 max_power     = b->get_max_power();
 bound_on      = b->get_start_up_limit();
 bound_down    = b->get_shut_down_limit();

 if( b->get_delta_ramp_up().empty() )
  delta_ramp_up = max_power;
 else
  delta_ramp_up = b->get_delta_ramp_up();

 if( b->get_delta_ramp_down().empty() )
  delta_ramp_down = max_power;
 else
  delta_ramp_down = b->get_delta_ramp_down();

 retrieve_term( quad_term   , b->get_quad_term() );
 retrieve_term( linear_term , b->get_linear_term() );
 retrieve_term( const_term  , b->get_const_term() );

 if( ! owned )
  f_Block->read_unlock();

 P.assign( time_horizon , 0 );
 U.assign( time_horizon , false );
 stage = start;
 f_solved = false;
 f_best_cost = TUEDPINF;

 }  // end( ThermalUnitExtDPSolver::load_parameters )

/*--------------------------------------------------------------------------*/

void ThermalUnitExtDPSolver::process_modifications( void )
{
 bool reload = false;

 while( f_mod_lock.test_and_set( std::memory_order_acquire ) )
  ;

 for( auto mod : v_mod )
  if( guts_of_process_modifications( mod.get() ) ) {
   reload = true;
   break;
   }

 v_mod.clear();

 f_mod_lock.clear( std::memory_order_release );

 if( reload )
  load_parameters();

 }  // end( ThermalUnitExtDPSolver::process_modifications )

/*--------------------------------------------------------------------------*/

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

double ThermalUnitExtDPSolver::eval( const PQFun & F , double p )
{
 if( F.empty() )
  return( TUEDPINF );
 // binary search could be used; linear scan is fine for small m
 for( const auto & pc : F )
  if( pc.left - 1e-12 <= p && p <= pc.right + 1e-12 )
   return( eval_piece( pc , p ) );
 return( TUEDPINF );
 }

/*--------------------------------------------------------------------------*/

double ThermalUnitExtDPSolver::argmin_piece( const PieceQuad & pc )
{
 if( pc.alfa > 1e-16 ) {
  double p = - pc.beta / ( 2.0 * pc.alfa );
  if( p < pc.left ) return( pc.left );
  if( p > pc.right ) return( pc.right );
  return( p );
  }
 // linear or constant piece
 if( pc.beta > 0 ) return( pc.left );
 if( pc.beta < 0 ) return( pc.right );
 return( pc.left );  // constant: pick left
 }

/*--------------------------------------------------------------------------*/

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

void ThermalUnitExtDPSolver::shift_by( PQFun & F , double c )
{
 for( auto & pc : F )
  pc.gamma += c;
 }

/*--------------------------------------------------------------------------*/

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

bool ThermalUnitExtDPSolver::is_dominated_by( const PQFun & F1 ,
                                              const PQFun & F2 , double eps )
{
 // returns true iff F1(p) >= F2(p) - eps on all of dom( F1 ); if dom( F1 )
 // has a point outside dom( F2 ) (where F2 = +INF), returns false.

 if( F1.empty() )
  return( true );        // empty = +INF, trivially dominated
 if( F2.empty() )
  return( false );

 const double tol = 1e-12;

 // walk through F1's pieces, verifying F2 covers each and F1 - F2 >= -eps
 for( std::size_t i = 0 ; i < F1.size() ; ++i ) {
  double l1 = F1[ i ].left;
  double r1 = F1[ i ].right;
  double p  = l1;
  std::size_t j = 0;
  // advance j to the first F2-piece that overlaps [p, r1]
  while( j < F2.size() && F2[ j ].right <= p + tol )
   ++j;
  // now consume [p, r1] using F2 pieces
  while( p < r1 - tol ) {
   if( j >= F2.size() )
    return( false );     // F2 doesn't cover [p, r1]
   double l2 = F2[ j ].left;
   double r2 = F2[ j ].right;
   if( l2 > p + tol )
    return( false );     // gap in F2 coverage
   double seg_l = p;
   double seg_r = std::min( r1 , r2 );
   // diff D(p) = (alfa1 - alfa2) p^2 + (beta1 - beta2) p + (gamma1 - gamma2)
   double a = F1[ i ].alfa  - F2[ j ].alfa;
   double b = F1[ i ].beta  - F2[ j ].beta;
   double c = F1[ i ].gamma - F2[ j ].gamma;
   double min_D;
   if( a > tol ) {
    double p_opt = - b / ( 2.0 * a );
    if( p_opt < seg_l ) p_opt = seg_l;
    if( p_opt > seg_r ) p_opt = seg_r;
    min_D = a * p_opt * p_opt + b * p_opt + c;
    }
   else {
    // concave (a < 0) or linear (a ~= 0): min at one of the endpoints
    double vl = a * seg_l * seg_l + b * seg_l + c;
    double vr = a * seg_r * seg_r + b * seg_r + c;
    min_D = std::min( vl , vr );
    }
   if( min_D < - eps )
    return( false );
   p = seg_r;
   if( r2 < r1 - tol )
    ++j;
   else
    break;               // F2[j] covers up to or past r1
   }
  }
 return( true );
 }

/*--------------------------------------------------------------------------*/

ThermalUnitExtDPSolver::PQFun ThermalUnitExtDPSolver::sliding_min(
 const PQFun & F , double ramp_up , double ramp_down ,
 double lo , double hi )
{
 // implements: G(p_t) = min_{ q in [p_t - ramp_up, p_t + ramp_down] } F(q)
 // output domain is clamped to [lo, hi]; empty pieces are dropped.
 //
 // three-case analysis (Wuijts et al. 2021, eq. (16)):
 // let p_star be the global minimiser of F over its domain; then
 //   - for p_t > p_star + ramp_up :  G(p_t) = F(p_t - ramp_up)
 //   - for p_t in [p_star - ramp_down, p_star + ramp_up] : G(p_t) = F(p_star)
 //   - for p_t < p_star - ramp_down : G(p_t) = F(p_t + ramp_down)
 // so pieces of F strictly left  of p_star shift by -ramp_down (new piece),
 // pieces strictly right of p_star shift by +ramp_up, and a flat middle
 // piece is inserted. The piece containing p_star is split and contributes
 // to both the left-shifted and the right-shifted parts, around the middle.

 PQFun G;
 if( F.empty() || lo > hi + 1e-12 )
  return( G );

 if( ramp_up < 0 ) ramp_up = 0;
 if( ramp_down < 0 ) ramp_down = 0;

 // find p_star and f_star over F's domain
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

 PQFun raw;
 raw.reserve( F.size() * 2 + 1 );

 // left-shifted pieces: portions of F with q <= p_star
 // new piece: G(p_t) = F(p_t + ramp_down) on [a - ramp_down, b - ramp_down]
 // coefficient substitution p_t = q - ramp_down, so q = p_t + ramp_down
 //   F(q) = alfa q^2 + beta q + gamma
 //        = alfa (p_t + rd)^2 + beta (p_t + rd) + gamma
 //        = alfa p_t^2 + (2 alfa rd + beta) p_t + (alfa rd^2 + beta rd + gamma)
 for( const auto & pc : F ) {
  double a = pc.left, b = pc.right;
  if( b <= p_star + tol ) {
   // entire piece <= p_star
   double na = pc.alfa;
   double nb = 2 * pc.alfa * ramp_down + pc.beta;
   double nc = pc.alfa * ramp_down * ramp_down +
               pc.beta * ramp_down + pc.gamma;
   raw.push_back( { na , nb , nc , a - ramp_down , b - ramp_down } );
   }
  else if( a < p_star - tol ) {
   // p_star is strictly inside [a, b]: left part [a, p_star]
   double na = pc.alfa;
   double nb = 2 * pc.alfa * ramp_down + pc.beta;
   double nc = pc.alfa * ramp_down * ramp_down +
               pc.beta * ramp_down + pc.gamma;
   raw.push_back( { na , nb , nc , a - ramp_down , p_star - ramp_down } );
   }
  }

 // the middle flat piece at value f_star on [p_star - rd, p_star + ru]
 if( ramp_up > tol || ramp_down > tol )
  raw.push_back( { 0 , 0 , f_star , p_star - ramp_down , p_star + ramp_up } );

 // right-shifted pieces: portions of F with q >= p_star
 // new piece: G(p_t) = F(p_t - ramp_up) on [a + ramp_up, b + ramp_up]
 // coefficient substitution p_t = q + ramp_up, so q = p_t - ramp_up
 //   F(q) = alfa (p_t - ru)^2 + beta (p_t - ru) + gamma
 //        = alfa p_t^2 + (-2 alfa ru + beta) p_t + (alfa ru^2 - beta ru + gamma)
 for( const auto & pc : F ) {
  double a = pc.left, b = pc.right;
  if( a >= p_star - tol ) {
   double na = pc.alfa;
   double nb = - 2 * pc.alfa * ramp_up + pc.beta;
   double nc = pc.alfa * ramp_up * ramp_up -
               pc.beta * ramp_up + pc.gamma;
   raw.push_back( { na , nb , nc , a + ramp_up , b + ramp_up } );
   }
  else if( b > p_star + tol ) {
   double na = pc.alfa;
   double nb = - 2 * pc.alfa * ramp_up + pc.beta;
   double nc = pc.alfa * ramp_up * ramp_up -
               pc.beta * ramp_up + pc.gamma;
   raw.push_back( { na , nb , nc , p_star + ramp_up , b + ramp_up } );
   }
  }

 // sort and clamp to [lo, hi]
 std::sort( raw.begin() , raw.end() ,
            []( const PieceQuad & a , const PieceQuad & b ) {
             return( a.left < b.left );
             } );

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

void ThermalUnitExtDPSolver::run_DP( void )
{
 if( time_horizon == 0 ) {
  f_best_cost = 0;
  stage = dp_OK;
  return;
  }

 // sparse per-time-step storage: f_F[ t ], f_tau[ t ], f_on[ t ] are
 // parallel vectors sorted by tau ascending; only reachable, non-pruned
 // F^tau functions are stored.
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

 if( init_up_down_time > 0 ) {
  // unit was on for init_up_down_time instants before t = 0, and is
  // still on at t = 0: p_0 respects ramping from initial_power.
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
  // c_off_* remain +INF
  }
 else {
  // init_up_down_time <= 0: unit is off at t = 0. Restart possible iff
  // it has been off for >= mdt instants strictly before t = 0.
  bool can_restart_t0 = ( Index( - init_up_down_time ) >= mdt );
  if( can_restart_t0 ) {
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
  c_off_any[ 0 ]  = 0;
  f_any_pred[ 0 ] = -1;
  if( Index( - init_up_down_time ) + 1 >= mdt ) {
   c_off_ready[ 0 ]  = 0;
   f_ready_pred[ 0 ] = -1;
   }
  }

 // prune at t = 0 (can't hurt even if single entry)
 prune_RRF_plus( f_F[ 0 ] , f_tau[ 0 ] , f_on[ 0 ] );

 // v_shutdown[ 0 ]
 if( n > 1 )
  compute_v_shutdown( 0 , bound_down[ 1 ] );

 // ----------------------------------------------------- main loop: t >=1 -
 for( Index t = 1 ; t < n ; ++t ) {

  // OFF-side updates
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
  if( best_ready == init_ready && init_ready < TUEDPINF )
   f_ready_pred[ t ] = -1;
  else if( best_ready == fresh_ready && fresh_ready < TUEDPINF )
   f_ready_pred[ t ] = fresh_h;
  else
   f_ready_pred[ t ] = f_ready_pred[ t - 1 ];

  // ON-side: build the new sparse list for time t
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

  // tau = 1 entry: restart from c_off_ready[t-1]
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

  // tau > 1 entries: from each surviving entry of f_F[t-1]
  for( std::size_t i = 0 ; i < f_F[ t - 1 ].size() ; ++i ) {
   PQFun F = sliding_min( f_F[ t - 1 ][ i ] , ru_prev , rd_prev , Plo , Phi );
   if( F.empty() ) continue;
   add_quadratic( F , quad_term[ t ] , linear_term[ t ] , const_term[ t ] );
   auto [ v , p ] = min_over( F , Plo , Phi );
   new_F  .push_back( std::move( F ) );
   new_tau.push_back( f_tau[ t - 1 ][ i ] + 1 );
   new_on .push_back( { v , p } );
   }

  // note: new_tau is already sorted ascending (tau=1 first, then the
  // previous-step taus in order, each incremented by 1; since a tau=1
  // entry at t-1 becomes tau=2 at t, the new restart (tau=1) correctly
  // leads).

  // RRF+ pruning
  prune_RRF_plus( new_F , new_tau , new_on );

  f_F  [ t ] = std::move( new_F );
  f_tau[ t ] = std::move( new_tau );
  f_on [ t ] = std::move( new_on );

  // v_shutdown[t] if there is a future off instant inside the horizon
  if( t < n - 1 )
   compute_v_shutdown( t , bound_down[ t + 1 ] );

  }  // end( for( t ) )

 // finalise best cost: min over { on at n-1 (no SD constraint),
 //                                off at n-1 regardless }
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

void ThermalUnitExtDPSolver::build_solution( void )
{
 std::fill( P.begin() , P.end() , 0.0 );
 std::fill( U.begin() , U.end() , false );

 if( ! f_solved ) {
  stage = sol_OK;
  return;
  }

 const Index n = time_horizon;
 if( n == 0 ) {
  stage = sol_OK;
  return;
  }

 // lookup helper: find index of entry with the given tau in f_tau[t], or
 // return f_tau[t].size() if not present
 auto find_tau = [ & ]( Index t , Index tau ) -> std::size_t {
  const auto & v = f_tau[ t ];
  auto it = std::lower_bound( v.begin() , v.end() , tau );
  if( it != v.end() && *it == tau )
   return( std::size_t( it - v.begin() ) );
  return( v.size() );
  };

 // pick the optimal final state (over the surviving tau values)
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

 // walk_on: same semantics as before, but now using sparse lookup to
 // find the predecessor's (tau-1) entry
 auto walk_on = [ & ]( Index t , Index tau , double p ) -> Index {
  while( true ) {
   P[ t ] = p;
   U[ t ] = true;
   if( tau == 1 )
    return( t );
   if( t == 0 )
    return( 0 );
   std::size_t idx = find_tau( t - 1 , tau - 1 );
   if( idx == f_tau[ t - 1 ].size() )
    // the predecessor was pruned (should not happen under RRF+ invariants
    // since a surviving F^tau_t must have come from a surviving
    // F^{tau-1}_{t-1}); if it does, we bail out to avoid an infinite loop
    return( t );
   double p_prev_star = f_on[ t - 1 ][ idx ].argmin_p;
   double ru = delta_ramp_up  [ t - 1 ];
   double rd = delta_ramp_down[ t - 1 ];
   double p_prev;
   if( p > p_prev_star + ru )
    p_prev = p - ru;
   else
    if( p < p_prev_star - rd )
     p_prev = p + rd;
    else
     p_prev = p_prev_star;
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

 if( best_on <= best_off ) {
  on_t   = n - 1;
  on_tau = best_tau;
  on_p   = best_p;
  have_on = true;
  }
 else {
  // final state is off at n-1; look up the shutdown that produced it
  have_on = false;
  // the off-period covers [ f_any_pred[ n-1 ] + 1 , n-1 ] (inclusive);
  // nothing to write on that range (P, U already 0/false); but we must
  // backtrack further if f_any_pred[ n-1 ] >= 0.
  int h = f_any_pred[ n - 1 ];
  if( h >= 0 ) {
   on_t   = Index( h );
   on_tau = v_shutdown_tau[ h ];
   on_p   = v_shutdown_p[ h ];
   have_on = true;
   }
  }

 while( have_on ) {
  Index t_start = walk_on( on_t , on_tau , on_p );

  if( t_start == 0 ) {
   // hit the beginning of the horizon: done
   have_on = false;
   break;
   }
  // else: on-run starts at t_start >= 1, preceded by off at t_start - 1.
  // Find the last shutdown before t_start that produced a "ready" state
  // at time t_start - 1, i.e., f_ready_pred[ t_start - 1 ].
  int h = f_ready_pred[ t_start - 1 ];
  if( h < 0 ) {
   // the off state at t_start - 1 came from the initial condition
   // (unit off from before the horizon); nothing more to write.
   have_on = false;
   break;
   }
  on_t   = Index( h );
  on_tau = v_shutdown_tau[ h ];
  on_p   = v_shutdown_p[ h ];
  }

 stage = sol_OK;

 }  // end( ThermalUnitExtDPSolver::build_solution )

/*--------------------------------------------------------------------------*/
/*-------------- End File ThermalUnitExtDPSolver.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
