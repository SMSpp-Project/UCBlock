/*--------------------------------------------------------------------------*/
/*-------------------- File NuclearUnitExtDPSolver.cpp ---------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the NuclearUnitExtDPSolver class: the labels of the
 * states of the DP of ThermalUnitExtDPSolver are the modulation lockout,
 * which enforces the modulation constraints of NuclearUnitBlock. See the
 * header for the model and the algorithm.
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
 // load all the base (thermal + reserve + reactive + design) parameters
 // first; this also resets the pipeline state (stage = start, P/U cleared)
 ThermalUnitExtDPSolver::load_parameters();

 bool owned = f_Block->is_owned_by( f_id );
 if( ( ! owned ) && ( ! f_Block->read_lock() ) )
  throw( std::runtime_error(
   "NuclearUnitExtDPSolver::load_parameters: unable to lock the Block." ) );

 auto b = static_cast< NuclearUnitBlock * >( f_Block );

 // the data of the operating rules [see NuclearRules]
 auto & r = f_rules;
 r.time_horizon = time_horizon;
 r.mod_interval = b->get_modulation_interval();
 r.init_modulation = int( b->get_initial_modulation() );
 r.mod_ramp_up = b->get_modulation_ramp_up();
 r.mod_ramp_down = b->get_modulation_ramp_down();
 f_can_stay = std::none_of( r.mod_ramp_up.begin() , r.mod_ramp_up.end() ,
                            []( double x ) { return( x < 0 ); } ) &&
              std::none_of( r.mod_ramp_down.begin() ,
                            r.mod_ramp_down.end() ,
                            []( double x ) { return( x < 0 ); } );

 r.max_mod_length = std::max( b->get_max_modulation_length() , Index( 1 ) );
 r.stab_start = b->get_stability_after_start_up();
 r.bands = b->get_power_bands();
 r.mod_per_day = b->get_modulations_per_day();
 r.deep_per_day = b->get_deep_decreases_per_day();
 r.starts_per_day = b->get_start_ups_per_day();
 r.day_length = b->get_day_length();
 r.direction = b->has_modulation_direction();
 r.down_cost = b->get_down_modulation_cost();
 r.deep_thr = b->get_deep_decrease_threshold();
 r.deep_grad = b->get_deep_decrease_gradient();
 r.deep_cost = b->get_deep_decrease_cost();
 r.inf = TUEDPINF;

 if( ! owned )
  f_Block->read_unlock();

 // the sizes of the parts of the labels: the (mode, lockout or steps)
 // pairs, and the ranges of the counters that are limited
 r.set_sizes();

 }  // end( NuclearUnitExtDPSolver::load_parameters )

/*--------------------------------------------------------------------------*/

bool NuclearUnitExtDPSolver::guts_of_process_modifications( const p_Mod mod )
{
 // a change in the modulation ramps or in the costs of the operating rules
 // is rare: reload everything
 if( auto tmod = dynamic_cast< NuclearUnitBlockMod * >( mod ) )
  if( ( tmod->type() == NuclearUnitBlockMod::eSetModDP ) ||
      ( tmod->type() == NuclearUnitBlockMod::eSetModDM ) ||
      ( tmod->type() == NuclearUnitBlockMod::eSetModCost ) ||
      ( tmod->type() == NuclearUnitBlockMod::eSetDeepCost ) )
   return( true );

 return( ThermalUnitExtDPSolver::guts_of_process_modifications( mod ) );

 }  // end( NuclearUnitExtDPSolver::guts_of_process_modifications )

/*--------------------------------------------------------------------------*/

void NuclearUnitExtDPSolver::load_fixings( void )
{
 ThermalUnitExtDPSolver::load_fixings();

 bool owned = f_Block->is_owned_by( f_id );
 if( ( ! owned ) && ( ! f_Block->read_lock() ) )
  throw( std::runtime_error(
   "NuclearUnitExtDPSolver::load_fixings: unable to lock the Block." ) );

 // the read lock is released however the method ends, a throw included
 struct ReadUnlock {
  Block * b;
  ~ReadUnlock() { if( b ) b->read_unlock(); }
  } read_unlock{ owned ? nullptr : f_Block };

 auto b = static_cast< NuclearUnitBlock * >( f_Block );

 // the fixed Variable of the rules: -1 where free, the fixed value where
 // fixed. The structural fixings to 0 of the initial conditions (the
 // instants in which the unit is still locked out, and those before the
 // first free commitment of a unit initially off) are read like any other,
 // the label forbidding those moves anyway
 auto scan = [ & ]( const ColVariable * v , std::vector< signed char > & f ) {
  f.clear();
  if( ! v )
   return;
  for( Index t = 0 ; t < time_horizon ; ++t )
   if( v[ t ].is_fixed() ) {
    const double val = v[ t ].get_value();
    if( ( val < -1e-9 ) || ( val > 1 + 1e-9 ) ||
        ( ( val > 1e-9 ) && ( val < 1 - 1e-9 ) ) )
     throw( std::logic_error( "NuclearUnitExtDPSolver::load_fixings: the "
      "fixed value of a Variable of the operating rules is not 0 or 1" ) );
    if( f.empty() )
     f.assign( time_horizon , -1 );
    f[ t ] = ( val > 0.5 ) ? 1 : 0;
    }
  };
 scan( b->get_const_modulation() , f_rules.fix_mod );
 scan( b->get_const_modulation_down() , f_rules.fix_down );
 scan( b->get_const_deep_decrease() , f_rules.fix_deep );

 // the other Variable of the rules the labels do not read, hence their
 // fixings cannot be honored, save the structural ones at 0 of the first
 // instant of a unit initially off [see generate_abstract_variables()]
 auto refuse = [ & ]( const ColVariable * v , Index n , const char * name ) {
  for( Index t = 0 ; v && ( t < n ) ; ++t )
   if( v[ t ].is_fixed() &&
       ( ! ( ( t == 0 ) && ( init_up_down_time <= 0 ) &&
	     ( v[ t ].get_value() == 0 ) ) ) )
    throw( std::logic_error( std::string( "NuclearUnitExtDPSolver::"
     "load_fixings: fixed " ) + name + " Variable not supported (yet)" ) );
  };
 refuse( b->get_const_deep_drop() , time_horizon , "deep drop" );
 refuse( b->get_const_deep_low() , time_horizon , "deep low" );
 refuse( b->get_const_modulation_start() , time_horizon ,
	 "modulation start" );
 refuse( b->get_const_modulation_end() , time_horizon , "modulation end" );
 refuse( b->get_const_band() , 3 * time_horizon , "band" );

 // a deep decrease fixed to 1 is, in the rows of NuclearUnitBlock, the cost
 // and the count of one at that instant, whatever the unit does there (no
 // row bounds it from above, save (32) of the tight rules): it is not a
 // move of the DP, which only has the decreases that are deep, hence it is
 // refused rather than solved as a stricter problem
 if( std::any_of( f_rules.fix_deep.begin() , f_rules.fix_deep.end() ,
                  []( signed char f ) { return( f > 0 ); } ) )
  throw( std::logic_error( "NuclearUnitExtDPSolver::load_fixings: a deep "
                           "decrease Variable fixed to 1 is not supported" ) );

 // a modulation (or a deep decrease) fixed to 1 needs the unit on at that
 // instant, since m_t <= u_t (and d_t <= m_t): it is a fixing ON of the
 // commitment for the tables of the base DP [see start_label() for the
 // start-up, which cannot happen at that instant either]
 bool forced = false;
 f_last_forced = -1;
 for( Index t = 0 ; t < time_horizon ; ++t )
  if( f_rules.on_forced( t ) ) {
   forced = true;
   f_last_forced = int( t );
   }
 if( forced ) {
  if( nxt_on.size() != time_horizon + 1 ) {
   nxt_on.assign( time_horizon + 1 , time_horizon );
   nxt_off.assign( time_horizon + 1 , time_horizon );
   }
  for( Index t = time_horizon ; t-- > 0 ; )
   if( f_rules.on_forced( t ) )
    nxt_on[ t ] = t;
   else
    if( nxt_on[ t ] != t )
     nxt_on[ t ] = nxt_on[ t + 1 ];
  f_has_fixings = f_must_build = true;
  }

 }  // end( NuclearUnitExtDPSolver::load_fixings )

/*--------------------------------------------------------------------------*/

bool NuclearUnitExtDPSolver::reads_group( const std::string & name ) const
{
 static const std::vector< std::string > read = {
  "m_thermal" , "m_down_nuclear" , "m_start_nuclear" , "band_nuclear" ,
  "m_end_nuclear" , "deep_nuclear" , "deep_drop_nuclear" ,
  "deep_low_nuclear" };

 return( ThermalUnitExtDPSolver::reads_group( name ) ||
         ( std::find( read.begin() , read.end() , name ) != read.end() ) );
 }

/*--------------------------------------------------------------------------*/
// The moves that cost something are the downward steps, possibly deep
// decreases.

double NuclearUnitExtDPSolver::min_move_cost( Index t ) const
{
 const double cdn = f_rules.down_cost.empty() ? 0.0 : f_rules.down_cost[ t ];
 const double cdd = f_rules.deep_cost.empty() ? 0.0 : f_rules.deep_cost[ t ];
 return( std::min( { 0.0 , cdn , cdd , cdn + cdd } ) );
 }

/*--------------------------------------------------------------------------*/

// A stable unit can keep its output by the stable move, whose window
// contains 0 and whose range is the band of the label, unless a modulation
// is imposed later by a fixed Variable of the rules.

bool NuclearUnitExtDPSolver::can_hold( Index t , Index lab ) const
{
 if( f_last_forced > int( t ) )
  return( false );
 if( ! f_can_stay )
  return( false );
 return( f_rules.on_label( lab ).mode == 0 );
 }

/*--------------------------------------------------------------------------*/

// The bound of a relaxation of the rules, computed backwards on a grid of
// cells of the output. The relaxation keeps the stable mode with its
// lockout and the number of modulations of the day of the labels, merges the
// two directions and all the steps of a modulation into one mode whose steps
// move within the ramps and that may last any number of steps, and drops the
// bands, the other counters of the day, the
// costs of the downward steps and of the deep decreases (which are not
// negative, or are bounded by min_move_cost()), the minimum up and down
// times and the fixings. An instant on in a cell costs at least the least
// cost of being on in the cell, the reward of the reserves on the whole band
// included; a move from a cell may land in any cell that meets the window of
// the move widened by the cell. Hence H_t( k , m , c ) below is not larger
// than the cost of the instants after t of any schedule that is on at t in
// the cell c with a label whose core is k and with m modulations in the day
// of t + 1, and Z_t( l , m ) than that of any schedule that is off at t with
// lockout l and m modulations entering t + 1.

// The least cost of being on at s with the output in [ a , b ], the reward
// of the reserves on the whole band (the piecewise linear disc, 0 outside
// its domain) included.

double NuclearUnitExtDPSolver::least_on_cost( Index s , const PQFun & disc ,
                                              double qd , double a ,
                                              double b ) const
{
 const PieceQuad f{ quad_term[ s ] , linear_term[ s ] , const_term[ s ] + qd ,
                    a , b };
 auto least = [ & ]( const PieceQuad & q , double l , double r ) {
  const PieceQuad x{ q.alfa , q.beta , q.gamma , l , r };
  return( eval_piece( x , argmin_piece( x ) ) );
  };
 if( disc.empty() )
  return( least( f , a , b ) );
 double v = TUEDPINF;
 if( a < disc.front().left )
  v = least( f , a , std::min( b , disc.front().left ) );
 if( b > disc.back().right )
  v = std::min( v , least( f , std::max( a , disc.back().right ) , b ) );
 for( const auto & pc : disc ) {
  const double l = std::max( a , pc.left ) , r = std::min( b , pc.right );
  if( l > r )
   continue;
  const PieceQuad g{ f.alfa + pc.alfa , f.beta + pc.beta , f.gamma + pc.gamma ,
                     l , r };
  v = std::min( v , least( g , l , r ) );
  }
 return( v );
 }

/*--------------------------------------------------------------------------*/

void NuclearUnitExtDPSolver::prepare_completion_bound(
 const std::vector< PQFun > & disc , const std::vector< double > & qd )
{
 f_cb_H.clear();
 const Index n = time_horizon;
 if( n < 2 )
  return;
 double lo = TUEDPINF , hi = - TUEDPINF;
 for( Index t = 0 ; t < n ; ++t ) {
  lo = std::min( lo , min_power[ t ] );
  hi = std::max( hi , max_power[ t ] );
  }
 if( ! ( hi > lo ) )
  return;
 // the cores: the stable ones with lockout 0 , ... , L - 1, then the
 // modulation (if any); each with M values of the counter
 const Index L = f_rules.lockout_max() + 1;
 const Index K = ( f_rules.max_mod_length > 1 ) ? 1 : 0;
 const Index NC = L + K;
 const Index M = f_rules.nc;
 // about 2^18 values at most, a cell not finer than 1 / 128 of the range
 const Index G = std::clamp( Index( ( std::size_t( 1 ) << 18 ) /
                                    ( std::size_t( n ) * NC * M ) ) ,
                             Index( 16 ) , Index( 128 ) );
 f_cb_p0 = lo;
 f_cb_h = ( hi - lo ) / G;
 f_cb_g = G;
 f_cb_nc = NC * M;
 f_cb_H.assign( std::size_t( n ) * NC * M * G , 0.0 );
 const Index B = f_rules.mod_lockout();
 const Index Bv = f_rules.stab_start ? f_rules.stab_start - 1 : 0;
 auto H = [ & ]( Index t , Index k , Index m ) {
  return( f_cb_H.data() + ( ( std::size_t( t ) * NC + k ) * M + m ) * G );
  };
 std::vector< double > Z( L * M , 0.0 ) , Zn( L * M ) , cst( G ) , A( G ) ,
                       W( G );
 std::vector< Index > dq( G );
 // W[ c ] = min of A over the cells that meet [ lo_c - a , hi_c + b ), the
 // window widened by the cell (and by a little more, which only weakens it)
 auto slide = [ & ]( double a , double b ) {
  const double h = f_cb_h;
  std::size_t head = 0 , tail = 0;
  Index next = 0;
  for( Index c = 0 ; c < G ; ++c ) {
   const double jl = std::floor( c - a / h - 1e-9 );
   const double jh = std::ceil( c + 1 + b / h + 1e-9 ) - 1;
   if( ( jh < 0 ) || ( jl > double( G - 1 ) ) || ( jl > jh ) ) {
    W[ c ] = TUEDPINF;
    continue;
    }
   const Index jlo = jl < 0 ? 0 : Index( jl );
   const Index jhi = jh >= G - 1 ? G - 1 : Index( jh );
   for( ; next <= jhi ; ++next ) {
    while( ( tail > head ) && ( A[ dq[ tail - 1 ] ] >= A[ next ] ) )
     --tail;
    dq[ tail++ ] = next;
    }
   while( ( tail > head ) && ( dq[ head ] < jlo ) )
    ++head;
   W[ c ] = ( tail > head ) ? A[ dq[ head ] ] : TUEDPINF;
   }
  };
 for( Index t = n - 1 ; t-- > 0 ; ) {
  const Index s = t + 1;
  // the least cost of being on at s in each cell; a move into s costs at
  // least min_move_cost( s ), if that is negative
  const double mmc = std::min( 0.0 , min_move_cost( s ) );
  for( Index c = 0 ; c < G ; ++c ) {
   const double a = std::max( f_cb_p0 + c * f_cb_h , min_power[ s ] );
   const double b = std::min( f_cb_p0 + ( c + 1 ) * f_cb_h , max_power[ s ] );
   if( a > b + 1e-12 ) {
    cst[ c ] = TUEDPINF;
    continue;
    }
   cst[ c ] = least_on_cost( s , disc[ s ] , qd[ s ] , a ,
                             std::max( a , b ) ) + mmc;
   }
  const double sdc = shutdown_costs.empty() ? 0.0 : shutdown_costs[ s ];
  const double suc = startup_costs.empty() ? 0.0 : startup_costs[ s ];
  const double ru = delta_ramp_up[ s ] , rd = delta_ramp_down[ s ];
  const double wu = std::min( ru , f_rules.mod_ramp_up[ s ] );
  const double wd = std::min( rd , f_rules.mod_ramp_down[ s ] );
  // the counter of the labels entering s + 1 is reset if a day begins there
  const bool newday = ( s + 1 < n ) &&
                      ( f_rules.day( s + 1 ) != f_rules.day( s ) );
  auto cnt = [ & ]( Index m ) { return( newday ? Index( 0 ) : m ); };
  auto at = [ & ]( Index k , Index m , Index c ) {
   return( cst[ c ] + H( s , k , cnt( m ) )[ c ] );
   };
  for( Index m = 0 ; m < M ; ++m ) {
   // the off side
   for( Index l = 0 ; l < L ; ++l ) {
    const Index l1 = l ? l - 1 : 0;
    double v = Z[ l1 * M + cnt( m ) ];
    const Index k = std::max( l1 , Bv );
    for( Index c = 0 ; c < G ; ++c )
     v = std::min( v , suc + at( k , m , c ) );
    Zn[ l * M + m ] = v;
    }
   // a modulation may start if the counter allows it
   const bool start = ( f_rules.mod_per_day < 0 ) ||
                      ( m < Index( f_rules.mod_per_day ) );
   const Index m1 = ( M > 1 ) ? std::min( m + 1 , M - 1 ) : 0;
   for( Index k = 0 ; k < NC ; ++k ) {
    double * Ht = H( t , k , m );
    if( k < L ) {                                    // stable, lockout k
     const Index k1 = k ? k - 1 : 0;
     for( Index c = 0 ; c < G ; ++c )
      A[ c ] = at( k1 , m , c );
     slide( wd , wu );
     for( Index c = 0 ; c < G ; ++c )
      Ht[ c ] = std::min( W[ c ] , sdc + Z[ k1 * M + cnt( m ) ] );
     if( ( k == 0 ) && start ) {
      for( Index c = 0 ; c < G ; ++c )
       A[ c ] = std::min( at( B , m1 , c ) ,
                          K ? at( L , m1 , c ) : TUEDPINF );
      slide( rd , ru );
      for( Index c = 0 ; c < G ; ++c )
       Ht[ c ] = std::min( Ht[ c ] , W[ c ] );
      }
     }
    else {                                           // k - L + 1 steps
     for( Index c = 0 ; c < G ; ++c )
      A[ c ] = std::min( at( B , m , c ) , at( k , m , c ) );
     slide( rd , ru );
     for( Index c = 0 ; c < G ; ++c )
      Ht[ c ] = W[ c ];
     }
    }
   }
  Z.swap( Zn );
  }
 }

/*--------------------------------------------------------------------------*/

bool NuclearUnitExtDPSolver::completion_bound( Index t , Index lab ,
                                               const PQFun & F ,
                                               const OnSlot & m ,
                                               double lev ) const
{
 if( f_cb_H.empty() || F.empty() )
  return( true );
 const NuclearRules::Label l = f_rules.on_label( lab );
 const Index k = ( l.mode == 0 ) ? l.lk : f_rules.lockout_max() + 1;
 const double * Ht = f_cb_H.data() +
                     ( std::size_t( t ) * f_cb_nc + k * f_rules.nc + l.c ) *
                     f_cb_g;
 auto cell = [ & ]( double p ) {
  const double f = std::floor( ( p - f_cb_p0 ) / f_cb_h );
  return( f < 0 ? Index( 0 ) : ( f >= f_cb_g ? f_cb_g - 1 : Index( f ) ) );
  };
 // the cell of the minimizer is often enough to keep the state
 if( m.min_val + Ht[ cell( m.argmin_p ) ] <= lev )
  return( true );
 for( const auto & pc : F ) {
  const Index c1 = cell( pc.right );
  for( Index c = cell( pc.left ) ; c <= c1 ; ++c ) {
   if( m.min_val + Ht[ c ] <= lev ) {   // else no point of c is below lev
    // the piece on the closed cell, which only makes the bound smaller
    const double a = std::max( pc.left , f_cb_p0 + c * f_cb_h );
    const double b = std::min( pc.right , f_cb_p0 + ( c + 1 ) * f_cb_h );
    const PieceQuad q{ pc.alfa , pc.beta , pc.gamma , std::min( a , b ) ,
                       std::max( a , b ) };
    if( eval_piece( q , argmin_piece( q ) ) + Ht[ c ] <= lev )
     return( true );
    }
   }
  }
 return( false );
 }

/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------- WRITING THE SOLUTION -----------------------------*/
/*--------------------------------------------------------------------------*/

// m_t = 1 exactly at the on instants reached by the move with modulation,
// which is the second one of on_moves()

Solution * NuclearUnitExtDPSolver::get_Solution( Configuration * solc )
{
 // the base packs everything but the modulation, and it does so into a
 // NuclearUnitBlockSolution because the shape is asked to the Block
 auto sol = static_cast< NuclearUnitBlockSolution * >(
                            ThermalUnitExtDPSolver::get_Solution( solc ) );

 const bool built = ( ! has_design ) || design_on;

 std::vector< double > m( time_horizon );
 for( Index i = 0 ; i < time_horizon ; ++i )
  m[ i ] = ( built && U[ i ] && ( U_move[ i ] >= 0 ) && ( U_move[ i ] & 1 ) )
           ? 1 : 0;
 sol->set_modulation( std::move( m ) );

 if( f_rules.direction ) {
  std::vector< double > d( time_horizon );
  for( Index i = 0 ; i < time_horizon ; ++i )
   d[ i ] = ( built && U[ i ] && ( U_move[ i ] >= 0 ) && ( U_move[ i ] & 2 ) )
            ? 1 : 0;
  sol->set_modulation_down( std::move( d ) );
  }

 return( sol );

 }  // end( NuclearUnitExtDPSolver::get_Solution )

/*--------------------------------------------------------------------------*/

void NuclearUnitExtDPSolver::get_var_solution( Configuration * solc )
{
 bool owned = f_Block->is_owned_by( f_id );
 if( ( ! owned ) && ( ! f_Block->lock( f_id ) ) )
  throw( std::runtime_error(
   "NuclearUnitExtDPSolver::get_var_solution: unable to lock the Block." ) );

 // the lock is released however the method ends, a throw included
 struct Unlock {
  Block * b;
  void * id;
  ~Unlock() { if( b ) b->unlock( id ); }
  } unlock_block{ owned ? nullptr : f_Block , f_id };

 auto b = static_cast< NuclearUnitBlock * >( f_Block );
 const bool built = ( ! has_design ) || design_on;

 auto tag = [ & ]( Index i , int bit ) -> double {
  return( ( built && U[ i ] && ( U_move[ i ] >= 0 ) && ( U_move[ i ] & bit ) )
          ? 1 : 0 );
  };

 if( auto mod_it = b->get_modulation() )
  for( Index i = 0 ; i < time_horizon ; ++i )
   ( mod_it++ )->set_value( tag( i , 1 ) );

 if( auto mod_it = b->get_modulation_down() )
  for( Index i = 0 ; i < time_horizon ; ++i )
   ( mod_it++ )->set_value( tag( i , 2 ) );

 // the start of the modulations
 if( auto s_it = b->get_modulation_start() )
  for( Index i = 0 ; i < time_horizon ; ++i )
   ( s_it++ )->set_value( std::max( tag( i , 1 ) -
                                    ( i ? tag( i - 1 , 1 ) : 0.0 ) , 0.0 ) );

 // the deep decrease and its two auxiliaries (a decrease by more than the
 // gradient, an output below the threshold) out of the recovered schedule,
 // at the instants with an on predecessor; the auxiliaries exist in the
 // default formulation of the rules only
 if( auto dd = b->get_deep_decrease() ) {
  auto ddrop = b->get_deep_drop();
  auto dlow = b->get_deep_low();
  for( Index i = 0 ; i < time_horizon ; ++i ) {
   dd[ i ].set_value( tag( i , 4 ) );
   if( ( ! ddrop ) || ( ! dlow ) )
    continue;
   const bool on_pred = built && U[ i ] &&
                        ( i ? bool( U[ i - 1 ] ) : ( init_up_down_time > 0 ) );
   const double prev = i ? P[ i - 1 ] : initial_power;
   const double tol = 1e-7 * std::max( 1.0 , std::abs( P[ i ] ) );
   ddrop[ i ].set_value( ( on_pred &&
                           ( prev - P[ i ] > f_rules.deep_grad[ i ] + tol ) )
                         ? 1 : 0 );
   dlow[ i ].set_value( ( on_pred &&
                          ( P[ i ] < f_rules.deep_thr[ i ] - tol ) )
                        ? 1 : 0 );
   }
  }

 // all the rest, the Block being already locked by this Solver
 ThermalUnitExtDPSolver::get_var_solution( solc );

 }  // end( NuclearUnitExtDPSolver::get_var_solution )

/*--------------------------------------------------------------------------*/
/*------------------ End File NuclearUnitExtDPSolver.cpp -------------------*/
/*--------------------------------------------------------------------------*/
