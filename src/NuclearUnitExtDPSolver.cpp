/*--------------------------------------------------------------------------*/
/*-------------------- File NuclearUnitExtDPSolver.cpp ---------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the NuclearUnitExtDPSolver class: the Ext hybrid DP of
 * ThermalUnitExtDPSolver augmented with the modulation lockout counter that
 * enforces the nuclear modulation constraints of NuclearUnitBlock.
 *
 * See the design document (reserves-DP, "Extending the DP to nuclear units")
 * and the class header for the model and the algorithm.
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
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "NuclearUnitExtDPSolver.h"

#include <algorithm>

/*--------------------------------------------------------------------------*/
/*-------------------------- NAMESPACE AND USING ---------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------------- FACTORY ---------------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_cpp_0( NuclearUnitExtDPSolver );

/*--------------------------------------------------------------------------*/
/*------------------------ DERIVED METHODS OF BASE -------------------------*/
/*--------------------------------------------------------------------------*/

void NuclearUnitExtDPSolver::set_Block( Block * block )
{
 if( block == f_Block )
  return;

 Solver::set_Block( block );

 if( block ) {
  if( typeid( NuclearUnitBlock ) != typeid( *f_Block ) )
   throw( std::runtime_error(
    "NuclearUnitExtDPSolver::set_Block: NuclearUnitBlock required." ) );
  load_parameters();
  }
 }

/*--------------------------------------------------------------------------*/

void NuclearUnitExtDPSolver::load_parameters( void )
{
 // load all the base (thermal + reserve + design) parameters first; this
 // also resets the pipeline state (stage = start, P/U cleared)
 ThermalUnitExtDPSolver::load_parameters();

 bool owned = f_Block->is_owned_by( f_id );
 if( ( ! owned ) && ( ! f_Block->read_lock() ) )
  throw( std::runtime_error(
   "NuclearUnitExtDPSolver::load_parameters: unable to lock the Block." ) );

 auto b = static_cast< NuclearUnitBlock * >( f_Block );

 f_mod_interval = b->get_modulation_interval();
 f_init_modulation = int( b->get_initial_modulation() );
 mod_ramp_up = b->get_modulation_ramp_up();
 mod_ramp_down = b->get_modulation_ramp_down();

 if( ! owned )
  f_Block->read_unlock();

 // the modulation profile output buffer
 M.assign( time_horizon , false );

 }  // end( NuclearUnitExtDPSolver::load_parameters )

/*--------------------------------------------------------------------------*/

bool NuclearUnitExtDPSolver::guts_of_process_modifications( const p_Mod mod )
{
 // a change in the modulation ramp data invalidates everything: force a full
 // reload (modulation ramp changes are rare, so a reload is acceptable)
 if( auto tmod = dynamic_cast< NuclearUnitBlockMod * >( mod ) )
  if( ( tmod->type() == NuclearUnitBlockMod::eSetModDP ) ||
      ( tmod->type() == NuclearUnitBlockMod::eSetModDM ) )
   return( true );

 // everything else: defer to the base handling
 return( ThermalUnitExtDPSolver::guts_of_process_modifications( mod ) );

 }  // end( NuclearUnitExtDPSolver::guts_of_process_modifications )

/*--------------------------------------------------------------------------*/
/*-------------------------------- THE DP ----------------------------------*/
/*--------------------------------------------------------------------------*/

void NuclearUnitExtDPSolver::run_DP( void )
{
 const Index n = time_horizon;
 if( n == 0 ) {
  f_best_cost = 0;
  stage = dp_OK;
  return;
  }

 // the lockout-augmented DP does not honor fixed Variable: unlike the base
 // class it does not kill the incompatible states, so it must refuse them
 // rather than silently ignoring them. The only tolerated fixings are the
 // structural ones with which ThermalUnitBlock encodes the initial
 // conditions (the commitments before t_init fixed to the initial state),
 // which this DP enforces natively anyway
 load_fixings();
 bool foreign = f_no_build ||
                ( f_must_build && ( init_up_down_time <= 0 ) );
 if( f_has_fixings ) {
  if( init_up_down_time > 0 )
   foreign = foreign || ( nxt_off[ 0 ] < n ) || ( nxt_on[ t_init ] < n );
  else
   foreign = foreign || ( nxt_on[ 0 ] < n ) || ( nxt_off[ t_init ] < n );
  }
 if( foreign )
  throw( std::logic_error( "NuclearUnitExtDPSolver: fixed Variable not "
                           "supported (yet)" ) );

 const Index mut = std::max( min_up_time   , Index( 1 ) );
 const Index mdt = std::max( min_down_time , Index( 1 ) );
 const Index L = std::max( f_mod_interval , Index( 2 ) );  // tau^M
 const Index lockmax = L - 1;
 // lockout entering t = 0 (max{ tau^M - InitModulation , 0 }); a negative
 // InitModulation is clamped to 0 (the Index cast would otherwise wrap)
 const Index im = ( f_init_modulation > 0 ) ? Index( f_init_modulation )
                                            : Index( 0 );
 const Index l0 = ( im < L ) ? ( L - im ) : Index( 0 );

 // ell after one idle (no-modulation / off) instant
 auto dec = []( Index ell ) -> Index { return( ell ? ell - 1 : 0 ); };
 // ell after k idle instants
 auto dec_by = []( Index ell , Index k ) -> Index {
  return( ell > k ? ell - k : Index( 0 ) ); };

 // -- per-step effective ramp half-widths (Remark "Effective per-step ramp
 // window"): the base ramp at step (t-1 -> t) is indexed by t-1 (by t for
 // the special t = 0 step), the modulation ramp of the same step by t; the
 // min keeps the DP faithful to both. "mod" uses the full thermal ramp.
 auto nomod_up = [ & ]( Index t ) -> double {
  double base = ( t == 0 ) ? delta_ramp_up[ 0 ] : delta_ramp_up[ t - 1 ];
  return( std::min( base , mod_ramp_up[ t ] ) ); };
 auto nomod_dn = [ & ]( Index t ) -> double {
  double base = ( t == 0 ) ? delta_ramp_down[ 0 ] : delta_ramp_down[ t - 1 ];
  return( std::min( base , mod_ramp_down[ t ] ) ); };
 auto mod_up = [ & ]( Index t ) -> double {
  double base = ( t == 0 ) ? delta_ramp_up[ 0 ] : delta_ramp_up[ t - 1 ];
  return( std::min( base , delta_ramp_up[ t ] ) ); };
 auto mod_dn = [ & ]( Index t ) -> double {
  double base = ( t == 0 ) ? delta_ramp_down[ 0 ] : delta_ramp_down[ t - 1 ];
  return( std::min( base , delta_ramp_down[ t ] ) ); };

 // -- reset all the per-time-step state --------------------------------- //
 f_F  .assign( n , std::vector< PQFun >() );
 f_tau.assign( n , std::vector< Index >() );
 f_on .assign( n , std::vector< OnSlot >() );
 f_lock     .assign( n , std::vector< Index >() );
 f_back_lock.assign( n , std::vector< Index >() );
 f_back_idx .assign( n , std::vector< std::size_t >() );

 co_ready       .assign( n , std::vector< double >( L , TUEDPINF ) );
 ready_pred     .assign( n , std::vector< int >( L , -1 ) );
 ready_pred_lock.assign( n , std::vector< Index >( L , 0 ) );
 vs    .assign( n , std::vector< double >( L , TUEDPINF ) );
 vs_tau.assign( n , std::vector< Index >( L , 0 ) );
 vs_p  .assign( n , std::vector< double >( L , 0.0 ) );
 vs_idx.assign( n , std::vector< std::size_t >( L , BAD ) );

 c_off_any  .assign( n , TUEDPINF );
 f_any_pred .assign( n , -1 );
 f_any_lock .assign( n , 0 );

 // per-period reserve discount g_t(p), added wherever f_t is (energy-only
 // path: all empty, no-op). Shared by all ON nodes at the same t.
 // NOTE: interior cap (max_power); the start-up / shut-down boundary
 // correction of the thermal solver is not mirrored in the nuclear DP.
 // With the common bound_on = bound_down = min_power data it is moot.
 std::vector< PQFun > eff_disc( n );
 for( Index t = 0 ; t < n ; ++t )
  eff_disc[ t ] = build_reserve_discount( t , max_power[ t ] );

 const bool startup_in_progress =
  has_ramp_up && ( initial_power < min_power[ 0 ] - 1e-9 );
 const bool shutdown_in_progress =
  has_ramp_down && ( initial_power > bound_down[ 0 ] + 1e-9 );

 // -- generalised RRF+ domination prune of an on-side slot list --------- //
 // slot i is pruned if some slot j has: cost pointwise <= i's, lockout
 // ell_j <= ell_i (at least as much modulation freedom), and either both
 // tau >= mut (min-up no longer distinguishes them) or tau_j == tau_i.
 auto prune = [ & ]( std::vector< PQFun > & vF , std::vector< Index > & vT ,
                     std::vector< Index > & vL , std::vector< OnSlot > & vO ,
                     std::vector< Index > & vBL ,
                     std::vector< std::size_t > & vBI ) {
  const std::size_t sz = vF.size();
  if( sz <= 1 )
   return;
  std::vector< char > keep( sz , 1 );
  for( std::size_t i = 0 ; i < sz ; ++i ) {
   if( ! keep[ i ] ) continue;
   for( std::size_t j = 0 ; j < sz ; ++j ) {
    if( ( j == i ) || ( ! keep[ j ] ) ) continue;
    const bool tau_ok = ( ( vT[ i ] >= mut ) && ( vT[ j ] >= mut ) ) ||
                        ( vT[ i ] == vT[ j ] );
    if( ! tau_ok ) continue;
    if( vL[ j ] > vL[ i ] ) continue;      // j less flexible: cannot dominate
    if( is_dominated_by( vF[ i ] , vF[ j ] ) ) {
     // tie-break to avoid mutual elimination of identical slots: when j is
     // not strictly better, only let the lower-index one survive
     if( is_dominated_by( vF[ j ] , vF[ i ] ) && ( vL[ j ] == vL[ i ] ) &&
         ( vT[ j ] == vT[ i ] ) && ( j > i ) )
      continue;
     keep[ i ] = 0;
     break;
     }
    }
   }
  std::size_t out = 0;
  for( std::size_t i = 0 ; i < sz ; ++i )
   if( keep[ i ] ) {
    if( out != i ) {
     vF [ out ] = std::move( vF[ i ] );
     vT [ out ] = vT [ i ];  vL [ out ] = vL [ i ];  vO [ out ] = vO [ i ];
     vBL[ out ] = vBL[ i ];  vBI[ out ] = vBI[ i ];
     }
    ++out;
    }
  vF .resize( out );  vT .resize( out );  vL .resize( out );
  vO .resize( out );  vBL.resize( out );  vBI.resize( out );
  };

 // -- v_shutdown at time t: per lockout ell, best (tau >= mut, p) slot -- //
 auto compute_vs = [ & ]( Index t , double sd_hi ) {
  double Plo = min_power[ t ];
  for( std::size_t i = 0 ; i < f_F[ t ].size() ; ++i ) {
   if( f_tau[ t ][ i ] < mut ) continue;
   const Index ell = f_lock[ t ][ i ];
   auto [ v , p ] = min_over( f_F[ t ][ i ] , Plo , sd_hi );
   if( v < vs[ t ][ ell ] ) {
    vs    [ t ][ ell ] = v;
    vs_tau[ t ][ ell ] = f_tau[ t ][ i ];
    vs_p  [ t ][ ell ] = p;
    vs_idx[ t ][ ell ] = i;
    }
   }
  };

 // ============================================================= t = 0 ====
 if( init_up_down_time > 0 ) {
  // unit on at t = 0 with run-length tau0 = init + 1
  const Index tau0 = Index( init_up_down_time + 1 );
  // (a) no-modulation step from the initial power
  {
   double lo = std::max( min_power[ 0 ] , initial_power - nomod_dn( 0 ) );
   double hi = std::min( max_power[ 0 ] , initial_power + nomod_up( 0 ) );
   if( lo < hi + 1e-12 ) {
    PQFun F;
    F.push_back( { quad_term[ 0 ] , linear_term[ 0 ] , const_term[ 0 ] ,
                   lo , hi } );
    add_pwq( F , eff_disc[ 0 ] );
    auto [ v , p ] = min_over( F , lo , hi );
    f_F  [ 0 ].push_back( std::move( F ) );
    f_tau[ 0 ].push_back( tau0 );
    f_lock[ 0 ].push_back( dec( l0 ) );
    f_on [ 0 ].push_back( { v , p } );
    f_back_lock[ 0 ].push_back( 0 );
    f_back_idx [ 0 ].push_back( BAD );
    }
   }
  // (b) modulation step (only if free to modulate entering t = 0)
  if( l0 == 0 ) {
   double lo = std::max( min_power[ 0 ] , initial_power - mod_dn( 0 ) );
   double hi = std::min( max_power[ 0 ] , initial_power + mod_up( 0 ) );
   if( lo < hi + 1e-12 ) {
    PQFun F;
    F.push_back( { quad_term[ 0 ] , linear_term[ 0 ] , const_term[ 0 ] ,
                   lo , hi } );
    add_pwq( F , eff_disc[ 0 ] );
    auto [ v , p ] = min_over( F , lo , hi );
    f_F  [ 0 ].push_back( std::move( F ) );
    f_tau[ 0 ].push_back( tau0 );
    f_lock[ 0 ].push_back( lockmax );
    f_on [ 0 ].push_back( { v , p } );
    f_back_lock[ 0 ].push_back( 0 );
    f_back_idx [ 0 ].push_back( BAD );
    }
   }
  // free pre-horizon shutdown (mirrors the base): only if min-up already
  // satisfied and the unit is not still inside a start-up / shut-down ramp
  if( ( Index( init_up_down_time ) >= min_up_time ) &&
      ( ! startup_in_progress ) && ( ! shutdown_in_progress ) ) {
   c_off_any [ 0 ] = 0;
   f_any_pred[ 0 ] = -1;
   if( Index( 1 ) >= min_down_time ) {
    const Index lr = dec( l0 );
    co_ready  [ 0 ][ lr ] = 0;
    ready_pred[ 0 ][ lr ] = -1;
    }
   }
  }
 else {
  // unit off at t = 0 for |init| pre-horizon instants
  const bool can_restart_t0 = ( Index( - init_up_down_time ) >= mdt );
  if( can_restart_t0 ) {
   double lo = min_power[ 0 ];
   double hi = std::min( max_power[ 0 ] , bound_on[ 0 ] );
   if( lo < hi + 1e-12 ) {
    double suc = startup_costs.empty() ? 0.0 : startup_costs[ 0 ];
    PQFun F;
    F.push_back( { quad_term[ 0 ] , linear_term[ 0 ] ,
                   const_term[ 0 ] + suc , lo , hi } );
    add_pwq( F , eff_disc[ 0 ] );
    auto [ v , p ] = min_over( F , lo , hi );
    f_F  [ 0 ].push_back( std::move( F ) );
    f_tau[ 0 ].push_back( 1 );
    f_lock[ 0 ].push_back( dec( l0 ) );  // m_0 = 0 forced at start-up
    f_on [ 0 ].push_back( { v , p } );
    f_back_lock[ 0 ].push_back( l0 );    // off-ready lockout into 0
    f_back_idx [ 0 ].push_back( BAD );
    }
   }
  c_off_any [ 0 ] = 0;
  f_any_pred[ 0 ] = -1;
  if( Index( - init_up_down_time ) + 1 >= mdt ) {
   const Index lr = dec( l0 );
   co_ready  [ 0 ][ lr ] = 0;
   ready_pred[ 0 ][ lr ] = -1;
   }
  }

 prune( f_F[ 0 ] , f_tau[ 0 ] , f_lock[ 0 ] , f_on[ 0 ] ,
        f_back_lock[ 0 ] , f_back_idx[ 0 ] );

 if( n > 1 )
  compute_vs( 0 , bound_down[ 1 ] );

 // ============================================================ t >= 1 ====
 for( Index t = 1 ; t < n ; ++t ) {

  // -- OFF-side: c_off_any (1D, never feeds a restart) ----------------- //
  {
   double stay = c_off_any[ t - 1 ];
   double fresh = TUEDPINF;  Index fresh_ell = 0;
   for( Index e = 0 ; e < L ; ++e )
    if( vs[ t - 1 ][ e ] < fresh ) {
     fresh = vs[ t - 1 ][ e ];  fresh_ell = e;
     }
   if( fresh < stay ) {
    c_off_any [ t ] = fresh;
    f_any_pred[ t ] = int( t - 1 );
    f_any_lock[ t ] = fresh_ell;
    }
   else {
    c_off_any [ t ] = stay;
    f_any_pred[ t ] = f_any_pred[ t - 1 ];
    f_any_lock[ t ] = f_any_lock[ t - 1 ];
    }
   }

  // -- OFF-side: co_ready[t][.] (scatter into targets) ----------------- //
  // helper that relaxes co_ready[t][tgt] with priority init > fresh > stay
  auto relax_ready = [ & ]( Index tgt , double cost , int h , Index ellh ) {
   if( cost < co_ready[ t ][ tgt ] ) {
    co_ready  [ t ][ tgt ] = cost;
    ready_pred[ t ][ tgt ] = h;
    ready_pred_lock[ t ][ tgt ] = ellh;
    }
   };
  // init-off trail (highest priority on ties): reaches "ready" at t with
  // lockout dec_by( l0 , t + 1 ) (l0 idle instants 0..t)
  {
   bool init_ready = false;
   if( ( init_up_down_time <= 0 ) &&
       ( Index( - init_up_down_time ) + t + 1 >= mdt ) )
    init_ready = true;
   else if( ( init_up_down_time > 0 ) &&
            ( Index( init_up_down_time ) >= min_up_time ) &&
            ( t + 1 >= mdt ) &&
            ( ! startup_in_progress ) && ( ! shutdown_in_progress ) )
    init_ready = true;
   if( init_ready )
    relax_ready( dec_by( l0 , t + 1 ) , 0 , -1 , 0 );
   }
  // fresh long shutdown arc: shutdown at end of h = t - mdt with lockout
  // ell_h into h+1, then mdt idle instants -> dec_by( ell_h , mdt )
  if( t >= mdt ) {
   const Index h = t - mdt;
   for( Index e = 0 ; e < L ; ++e )
    if( vs[ h ][ e ] < TUEDPINF )
     relax_ready( dec_by( e , mdt ) , vs[ h ][ e ] , int( h ) , e );
   }
  // stay ready from t-1 (lowest priority on ties)
  for( Index e = 0 ; e < L ; ++e )
   if( co_ready[ t - 1 ][ e ] < TUEDPINF )
    relax_ready( dec( e ) , co_ready[ t - 1 ][ e ] ,
                 ready_pred[ t - 1 ][ e ] , ready_pred_lock[ t - 1 ][ e ] );

  // -- ON-side build --------------------------------------------------- //
  const double Plo = min_power[ t ];
  const double Phi = max_power[ t ];

  // (a) restart arcs (tau = 1): for each target lockout, take the cheapest
  // ready source mapping into it. dec( ell' ) = tgt, so tgt in [0, lockmax-1]
  for( Index tgt = 0 ; tgt < lockmax ; ++tgt ) {
   double best = TUEDPINF;  Index best_ell = 0;
   for( Index e = 0 ; e < L ; ++e )
    if( ( dec( e ) == tgt ) && ( co_ready[ t - 1 ][ e ] < best ) ) {
     best = co_ready[ t - 1 ][ e ];  best_ell = e;
     }
   if( best >= TUEDPINF )
    continue;
   double lo = Plo;
   double hi = std::min( Phi , bound_on[ t ] );
   if( lo >= hi + 1e-12 )
    continue;
   double suc = startup_costs.empty() ? 0.0 : startup_costs[ t ];
   PQFun F;
   F.push_back( { quad_term[ t ] , linear_term[ t ] ,
                  const_term[ t ] + suc + best , lo , hi } );
   add_pwq( F , eff_disc[ t ] );
   auto [ v , p ] = min_over( F , lo , hi );
   f_F  [ t ].push_back( std::move( F ) );
   f_tau[ t ].push_back( 1 );
   f_lock[ t ].push_back( tgt );
   f_on [ t ].push_back( { v , p } );
   f_back_lock[ t ].push_back( best_ell );  // off-ready lockout used
   f_back_idx [ t ].push_back( BAD );
   }

  // (b) continuation arcs (tau > 1) from each surviving slot at t-1
  for( std::size_t i = 0 ; i < f_F[ t - 1 ].size() ; ++i ) {
   const Index taup = f_tau [ t - 1 ][ i ];
   const Index ellp = f_lock[ t - 1 ][ i ];

   // Option A: no modulation (always available)
   {
    PQFun F;
    sliding_min( f_F[ t - 1 ][ i ] , nomod_up( t ) , nomod_dn( t ) ,
                 Plo , Phi , F );
    if( ! F.empty() ) {
     add_quadratic( F , quad_term[ t ] , linear_term[ t ] , const_term[ t ] );
     add_pwq( F , eff_disc[ t ] );
     auto [ v , p ] = min_over( F , Plo , Phi );
     f_F  [ t ].push_back( std::move( F ) );
     f_tau[ t ].push_back( taup + 1 );
     f_lock[ t ].push_back( dec( ellp ) );
     f_on [ t ].push_back( { v , p } );
     f_back_lock[ t ].push_back( ellp );
     f_back_idx [ t ].push_back( i );
     }
    }

   // Option B: modulation (only if free, ell' == 0); full ramp; ell = lockmax
   if( ellp == 0 ) {
    PQFun F;
    sliding_min( f_F[ t - 1 ][ i ] , mod_up( t ) , mod_dn( t ) ,
                 Plo , Phi , F );
    if( ! F.empty() ) {
     add_quadratic( F , quad_term[ t ] , linear_term[ t ] , const_term[ t ] );
     add_pwq( F , eff_disc[ t ] );
     auto [ v , p ] = min_over( F , Plo , Phi );
     f_F  [ t ].push_back( std::move( F ) );
     f_tau[ t ].push_back( taup + 1 );
     f_lock[ t ].push_back( lockmax );
     f_on [ t ].push_back( { v , p } );
     f_back_lock[ t ].push_back( 0 );
     f_back_idx [ t ].push_back( i );
     }
    }
   }

  prune( f_F[ t ] , f_tau[ t ] , f_lock[ t ] , f_on[ t ] ,
         f_back_lock[ t ] , f_back_idx[ t ] );

  if( t < n - 1 )
   compute_vs( t , bound_down[ t + 1 ] );

  }  // end( for( t ) )

 // -- finalise best cost ------------------------------------------------ //
 double best_on = TUEDPINF;
 for( const auto & slot : f_on[ n - 1 ] )
  if( slot.min_val < best_on )
   best_on = slot.min_val;
 double best_off = c_off_any[ n - 1 ];

 f_best_cost = std::min( best_on , best_off );

 // design (investment) decision, identical rule to the base solver
 if( has_design ) {
  if( f_best_cost + design_cost <= 0 ) {
   design_on = true;
   f_best_cost += design_cost;
   }
  else {
   design_on = false;
   f_best_cost = 0;
   }
  }

 f_solved = ( f_best_cost < TUEDPINF );
 stage = dp_OK;

 }  // end( NuclearUnitExtDPSolver::run_DP )

/*--------------------------------------------------------------------------*/
/*--------------------------- BUILDING SOLUTION ----------------------------*/
/*--------------------------------------------------------------------------*/

void NuclearUnitExtDPSolver::build_solution( void )
{
 std::fill( P.begin() , P.end() , 0.0 );
 std::fill( U.begin() , U.end() , false );
 std::fill( M.begin() , M.end() , false );

 if( ! f_solved ) { stage = sol_OK; return; }

 const Index n = time_horizon;
 if( n == 0 ) { stage = sol_OK; return; }

 const Index L = std::max( f_mod_interval , Index( 2 ) );
 const Index lockmax = L - 1;

 auto dec = []( Index ell ) -> Index { return( ell ? ell - 1 : 0 ); };

 auto nomod_up = [ & ]( Index t ) -> double {
  double base = ( t == 0 ) ? delta_ramp_up[ 0 ] : delta_ramp_up[ t - 1 ];
  return( std::min( base , mod_ramp_up[ t ] ) ); };
 auto nomod_dn = [ & ]( Index t ) -> double {
  double base = ( t == 0 ) ? delta_ramp_down[ 0 ] : delta_ramp_down[ t - 1 ];
  return( std::min( base , mod_ramp_down[ t ] ) ); };
 auto mod_up = [ & ]( Index t ) -> double {
  double base = ( t == 0 ) ? delta_ramp_up[ 0 ] : delta_ramp_up[ t - 1 ];
  return( std::min( base , delta_ramp_up[ t ] ) ); };
 auto mod_dn = [ & ]( Index t ) -> double {
  double base = ( t == 0 ) ? delta_ramp_down[ 0 ] : delta_ramp_down[ t - 1 ];
  return( std::min( base , delta_ramp_down[ t ] ) ); };

 // walk one on-run back from on-slot (t, idx) with power p, writing
 // P/U/M; returns the run-start time and, via off_lock, the off-ready
 // lockout used by the restart (NONE if the run hit t = 0)
 const Index NONE = Index( -1 );
 auto walk_on = [ & ]( Index t , std::size_t idx , double p ,
                       Index & off_lock ) -> Index {
  while( true ) {
   const Index tau = f_tau [ t ][ idx ];
   const Index ell = f_lock[ t ][ idx ];
   P[ t ] = p;
   U[ t ] = true;
   M[ t ] = ( ell == lockmax );   // a slot at lockmax is reached only by mod
   if( tau == 1 ) { off_lock = f_back_lock[ t ][ idx ]; return( t ); }
   if( t == 0 )   { off_lock = NONE; return( 0 ); }
   const std::size_t pidx = f_back_idx[ t ][ idx ];
   if( pidx == BAD ) { off_lock = NONE; return( t ); }  // defensive
   const double p_prev_star = f_on[ t - 1 ][ pidx ].argmin_p;
   const bool was_mod = ( ell == lockmax );
   const double ru = was_mod ? mod_up( t ) : nomod_up( t );
   const double rd = was_mod ? mod_dn( t ) : nomod_dn( t );
   double p_prev;
   if( p > p_prev_star + ru )      p_prev = p - ru;
   else if( p < p_prev_star - rd ) p_prev = p + rd;
   else                            p_prev = p_prev_star;
   double lo_prev = min_power[ t - 1 ];
   double hi_prev = max_power[ t - 1 ];
   if( f_tau[ t - 1 ][ pidx ] == 1 )
    hi_prev = std::min( hi_prev , bound_on[ t - 1 ] );
   if( p_prev < lo_prev ) p_prev = lo_prev;
   if( p_prev > hi_prev ) p_prev = hi_prev;
   t = t - 1;  idx = pidx;  p = p_prev;
   }
  };

 // choose the optimal terminal state
 double best_on = TUEDPINF;  std::size_t best_idx = BAD;
 for( std::size_t i = 0 ; i < f_on[ n - 1 ].size() ; ++i )
  if( f_on[ n - 1 ][ i ].min_val < best_on ) {
   best_on = f_on[ n - 1 ][ i ].min_val;  best_idx = i;
   }
 double best_off = c_off_any[ n - 1 ];

 // bootstrap the alternating on/off walk
 bool have_on = false;
 Index on_t = 0;  std::size_t on_idx = BAD;  double on_p = 0;
 if( best_on <= best_off ) {
  on_t = n - 1;  on_idx = best_idx;
  on_p = f_on[ n - 1 ][ best_idx ].argmin_p;
  have_on = ( best_idx != BAD );
  }
 else {
  int h = f_any_pred[ n - 1 ];
  if( h >= 0 ) {
   const Index ellh = f_any_lock[ n - 1 ];
   on_t = Index( h );  on_idx = vs_idx[ h ][ ellh ];
   on_p = vs_p[ h ][ ellh ];
   have_on = ( on_idx != BAD );
   }
  }

 // alternate on-runs and off-gaps
 while( have_on ) {
  Index off_lock = NONE;
  Index t_start = walk_on( on_t , on_idx , on_p , off_lock );

  if( t_start == 0 ) break;
  if( off_lock == NONE ) break;

  // the on-run started at t_start with a restart; the off-ready state used
  // is co_ready[ t_start - 1 ][ off_lock ], whose shutdown origin is
  const int h = ready_pred[ t_start - 1 ][ off_lock ];
  if( h < 0 ) break;   // ready came from the initial off trail (pre-horizon)
  const Index ellh = ready_pred_lock[ t_start - 1 ][ off_lock ];
  on_t = Index( h );
  on_idx = vs_idx[ h ][ ellh ];
  on_p = vs_p  [ h ][ ellh ];
  if( on_idx == BAD ) break;
  }

 stage = sol_OK;

 }  // end( NuclearUnitExtDPSolver::build_solution )

/*--------------------------------------------------------------------------*/
/*----------------------- WRITING THE SOLUTION -----------------------------*/
/*--------------------------------------------------------------------------*/

void NuclearUnitExtDPSolver::get_var_solution( Configuration * solc )
{
 bool owned = f_Block->is_owned_by( f_id );
 if( ( ! owned ) && ( ! f_Block->lock( f_id ) ) )
  throw( std::runtime_error(
   "NuclearUnitExtDPSolver::get_var_solution: unable to lock the Block." ) );

 auto b = static_cast< NuclearUnitBlock * >( f_Block );

 const bool built = ( ! has_design ) || design_on;
 if( has_design )
  b->get_design().set_value( design_on ? 1 : 0 );

 if( auto pow_it = b->get_active_power( 0 ) )
  for( Index i = 0 ; i < time_horizon ; ++i )
   ( pow_it++ )->set_value( built ? P[ i ] : 0 );

 if( auto com_it = b->get_commitment( 0 ) )
  for( Index i = 0 ; i < time_horizon ; ++i )
   ( com_it++ )->set_value( ( built && U[ i ] ) ? 1 : 0 );

 // the modulation indicators recovered by build_solution()
 if( auto mod_it = b->get_modulation() )
  for( Index i = 0 ; i < time_horizon ; ++i )
   ( mod_it++ )->set_value( ( built && M[ i ] ) ? 1 : 0 );

 if( auto pr_it = b->get_primary_spinning_reserve( 0 ) )
  for( Index i = 0 ; i < time_horizon ; ++i ) {
   double pr , sr;
   reserve_alloc( i , built ? P[ i ] : 0 , pr , sr , max_power[ i ] );
   ( pr_it++ )->set_value( pr );
   }

 if( auto sr_it = b->get_secondary_spinning_reserve( 0 ) )
  for( Index i = 0 ; i < time_horizon ; ++i ) {
   double pr , sr;
   reserve_alloc( i , built ? P[ i ] : 0 , pr , sr , max_power[ i ] );
   ( sr_it++ )->set_value( sr );
   }

 b->set_solution();

 if( ! owned )
  f_Block->unlock( f_id );

 }  // end( NuclearUnitExtDPSolver::get_var_solution )

/*--------------------------------------------------------------------------*/
/*------------------ End File NuclearUnitExtDPSolver.cpp -------------------*/
/*--------------------------------------------------------------------------*/
