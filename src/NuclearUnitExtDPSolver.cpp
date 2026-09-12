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

 f_mod_interval = b->get_modulation_interval();
 f_init_modulation = int( b->get_initial_modulation() );
 mod_ramp_up = b->get_modulation_ramp_up();
 mod_ramp_down = b->get_modulation_ramp_down();

 f_max_mod_length = std::max( b->get_max_modulation_length() , Index( 1 ) );
 f_stab_start = b->get_stability_after_start_up();
 f_bands = b->get_power_bands();
 f_mod_per_day = b->get_modulations_per_day();
 f_deep_per_day = b->get_deep_decreases_per_day();
 f_starts_per_day = b->get_start_ups_per_day();
 f_day_length = b->get_day_length();
 f_direction = b->has_modulation_direction();
 f_down_cost = b->get_down_modulation_cost();
 f_deep_thr = b->get_deep_decrease_threshold();
 f_deep_grad = b->get_deep_decrease_gradient();
 f_deep_cost = b->get_deep_decrease_cost();

 if( ! owned )
  f_Block->read_unlock();

 // the sizes of the parts of the labels: the (mode, lockout or steps)
 // pairs, and the ranges of the counters that are limited
 f_ncore = lockout_max() + 1 + 2 * ( f_max_mod_length - 1 );
 f_nc = ( f_mod_per_day >= 0 ) ? Index( f_mod_per_day ) + 1 : 1;
 f_na = ( ( ! f_deep_thr.empty() ) && ( f_deep_per_day >= 0 ) )
        ? Index( f_deep_per_day ) + 1 : 1;
 f_nv = ( f_starts_per_day >= 0 ) ? Index( f_starts_per_day ) + 1 : 1;
 f_nband = f_bands.empty() ? 1 : 3;
 f_ncount = f_nband * f_nc * f_na * f_nv;

 }  // end( NuclearUnitExtDPSolver::load_parameters )

/*--------------------------------------------------------------------------*/

bool NuclearUnitExtDPSolver::guts_of_process_modifications( const p_Mod mod )
{
 // a change in the modulation ramps is rare: reload everything
 if( auto tmod = dynamic_cast< NuclearUnitBlockMod * >( mod ) )
  if( ( tmod->type() == NuclearUnitBlockMod::eSetModDP ) ||
      ( tmod->type() == NuclearUnitBlockMod::eSetModDM ) )
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
 scan( b->get_const_modulation() , f_fix_mod );
 scan( b->get_const_modulation_down() , f_fix_down );
 scan( b->get_const_deep_decrease() , f_fix_deep );

 if( ! owned )
  f_Block->read_unlock();

 }  // end( NuclearUnitExtDPSolver::load_fixings )

/*--------------------------------------------------------------------------*/
/*------------------------- THE LABELS OF THE STATES -----------------------*/
/*--------------------------------------------------------------------------*/

// The lockout entering t = 0: the unit last modulated InitModulation
// instants before 0, hence it is locked out for tau^M - InitModulation more
// instants (none if that is not positive).

NuclearUnitExtDPSolver::Index NuclearUnitExtDPSolver::init_label( void ) const
{
 // the stable state with the initial lockout and no count yet: its on- and
 // off-code coincide. The lockout the unit enters the horizon with comes
 // from the last modulation, hence from tau^M, and not from the largest
 // lockout a label may carry, which the stability after a start-up may
 // have made larger
 const Index L = mod_lockout() + 1;
 const Index im = ( f_init_modulation > 0 ) ? Index( f_init_modulation )
                                            : Index( 0 );
 Label l;
 l.mode = 0;
 l.lk = ( im < L ) ? L - im : Index( 0 );
 l.c = l.a = l.s = 0;
 l.b = band_of( initial_power );
 return( on_code( l ) );
 }

/*--------------------------------------------------------------------------*/

NuclearUnitExtDPSolver::Label NuclearUnitExtDPSolver::on_label( Index lab )
 const
{
 Label l;
 Index core = lab % f_ncore;
 Index cnt = lab / f_ncore;
 const Index B = lockout_max();
 if( core <= B ) {
  l.mode = 0;
  l.lk = core;
  }
 else {
  core -= B;                                // 1 .. 2 ( L^M - 1 )
  l.mode = ( core < f_max_mod_length ) ? 1 : 2;
  l.lk = ( l.mode == 1 ) ? core : core - ( f_max_mod_length - 1 );
  }
 l.b = cnt % f_nband;
 cnt /= f_nband;
 l.c = cnt % f_nc;
 cnt /= f_nc;
 l.a = cnt % f_na;
 l.s = cnt / f_na;
 return( l );
 }

/*--------------------------------------------------------------------------*/

NuclearUnitExtDPSolver::Label NuclearUnitExtDPSolver::off_label( Index e )
 const
{
 Label l;
 const Index L = lockout_max() + 1;
 l.mode = 0;
 l.lk = e % L;
 Index cnt = e / L;
 l.b = cnt % f_nband;
 cnt /= f_nband;
 l.c = cnt % f_nc;
 cnt /= f_nc;
 l.a = cnt % f_na;
 l.s = cnt / f_na;
 return( l );
 }

/*--------------------------------------------------------------------------*/

NuclearUnitExtDPSolver::Index NuclearUnitExtDPSolver::shut_label(
                                               Index t , Index lab ) const
{
 Label l = on_label( lab );
 if( l.mode != 0 )                    // no shut-down during a modulation
  return( NO_LABEL );
 l.b = 0;              // an off unit has no output, hence no band
 return( off_code( l ) );
 }

/*--------------------------------------------------------------------------*/

NuclearUnitExtDPSolver::Index NuclearUnitExtDPSolver::idle_label(
                                        Index t , Index e , Index k ) const
{
 Label l = off_label( e );
 l.lk = ( l.lk > k ) ? l.lk - k : 0;
 if( day( t + k ) != day( t ) )
  l.c = l.a = l.s = 0;
 return( off_code( l ) );
 }

/*--------------------------------------------------------------------------*/

void NuclearUnitExtDPSolver::start_labels( Index t , Index e ,
    std::vector< std::pair< Index , std::pair< double , double > > > & ls )
 const
{
 ls.clear();
 const Index lab = start_label( t , e );
 if( lab == NO_LABEL )
  return;
 if( f_bands.empty() ) {          // one label over the whole range
  ls.push_back( { lab , { - TUEDPINF , TUEDPINF } } );
  return;
  }
 // one label per band, each over the range of its own band: which band a
 // unit restarts in is decided by the power it restarts at
 Label l = on_label( lab );
 for( Index b = 0 ; b < f_nband ; ++b ) {
  l.b = b;
  ls.push_back( { on_code( l ) ,
                  { b ? f_bands[ b - 1 ] : - TUEDPINF ,
                    ( b + 1 < f_nband ) ? f_bands[ b ] : TUEDPINF } } );
  }
 }

/*--------------------------------------------------------------------------*/

NuclearUnitExtDPSolver::Index NuclearUnitExtDPSolver::start_label(
                                                   Index t , Index e ) const
{
 Label l = off_label( e );
 if( ( f_starts_per_day >= 0 ) && ( l.s >= Index( f_starts_per_day ) ) )
  return( NO_LABEL );                 // the start-ups of the day are over
 if( f_nv > 1 )
  ++l.s;
 // no modulation at a start-up, and none for the A - 1 instants that
 // follow it either, which is the lockout the restart is born with
 l.lk = std::max( l.lk ? l.lk - 1 : Index( 0 ) ,
                  f_stab_start ? f_stab_start - 1 : Index( 0 ) );
 if( day( t + 1 ) != day( t ) )
  l.c = l.a = l.s = 0;
 return( on_code( l ) );
 }

/*--------------------------------------------------------------------------*/

bool NuclearUnitExtDPSolver::label_dominates( Index a , Index b ) const
{
 const Label la = on_label( a );
 const Label lb = on_label( b );
 // two labels of different bands are not comparable: the band says where
 // the output is, not how much history the unit carries
 return( ( la.mode == lb.mode ) && ( la.b == lb.b ) && ( la.lk <= lb.lk ) &&
         ( la.c <= lb.c ) && ( la.a <= lb.a ) && ( la.s <= lb.s ) );
 }

/*--------------------------------------------------------------------------*/

// The moves out of an on-state with label lab at t - 1 [see the table in the
// header]. The ramp of the ThermalUnitBlock for the step t-1 -> t is indexed
// by t-1 (by 0 for the step into t = 0), those of the modulation by t, and
// each window is the intersection of the two; a full-ramp step is only
// possible if the ramp of the ThermalUnitBlock allows it. The tag of a move
// has bit 0 for a modulation step, bit 1 for a downward one and bit 2 for a
// deep decrease.

void NuclearUnitExtDPSolver::on_moves( Index t , Index lab ,
                                      std::vector< OnMove > & mv ) const
{
 const Index k = t ? t - 1 : 0;
 const double ru = delta_ramp_up[ k ];
 const double rd = delta_ramp_down[ k ];
 const double fu = delta_ramp_up[ t ];     // the full ramps of the step
 const double fd = delta_ramp_down[ t ];
 const double wu = std::min( ru , fu );
 const double wd = std::min( rd , fd );
 const double cdn = f_down_cost.empty() ? 0.0 : f_down_cost[ t ];
 const bool deep = ! f_deep_thr.empty();
 const Index B = mod_lockout();      // what a modulation leaves behind
 const bool newday = ( t + 1 < time_horizon ) && ( day( t + 1 ) != day( t ) );

 const Label from = on_label( lab );

 // a fixed Variable of the operating rules [see load_fixings()] admits only
 // the moves that agree with it: bit 0 of the tag of a move says that it is
 // a modulation step, bit 1 that it goes downwards and bit 2 that it is a
 // deep decrease
 auto agree = [ & ]( const std::vector< signed char > & f , bool what ) {
  return( f.empty() || ( f[ t ] < 0 ) || ( ( f[ t ] > 0 ) == what ) );
  };
 const bool no_deep = agree( f_fix_deep , false );
 const bool yes_deep = agree( f_fix_deep , true );

 // the range of the output in each band: the bands only exist if the two
 // breakpoints are there, otherwise there is the one range of everything
 auto band_lo = [ & ]( Index b ) {
  return( f_bands.empty() || ( b == 0 ) ? - TUEDPINF : f_bands[ b - 1 ] );
  };
 auto band_hi = [ & ]( Index b ) {
  return( f_bands.empty() || ( b + 1 >= f_nband ) ? TUEDPINF : f_bands[ b ] );
  };

 // append the move landing in label to, whose landing power is restricted
 // to [ lo , hi ], splitting it for the deep decrease if its window reaches
 // a decrease of the deep-decrease gradient
 auto emit = [ & ]( Label to , double w_up , double w_dn , double cost ,
                    int tag , double lo = - TUEDPINF ,
                    double hi = TUEDPINF ) {
  if( newday )
   to.c = to.a = to.s = 0;
  if( ( ! agree( f_fix_mod , tag & 1 ) ) ||
      ( ! agree( f_fix_down , tag & 2 ) ) )
   return;
  if( ( ! deep ) || ( w_dn < f_deep_grad[ t ] - 1e-9 ) ) {
   if( no_deep )
    mv.push_back( { on_code( to ) , w_up , w_dn , cost , lo , hi , tag } );
   return;
   }
  const double dg = f_deep_grad[ t ];
  const double th = f_deep_thr[ t ];
  const double wu_deep = std::min( w_up , - dg );
  // a decrease of at least the gradient to at most the threshold: deep
  if( yes_deep && ( ( f_deep_per_day < 0 ) ||
                    ( from.a < Index( f_deep_per_day ) ) ) ) {
   Label td = to;
   if( ( f_na > 1 ) && ( ! newday ) )
    ++td.a;
   mv.push_back( { on_code( td ) , wu_deep , w_dn ,
                   cost + ( f_deep_cost.empty() ? 0.0 : f_deep_cost[ t ] ) ,
                   lo , std::min( hi , th ) , tag | 4 } );
   }
  if( ! no_deep )      // the deep decrease is imposed: nothing else is left
   return;
  // the same decrease to at least the threshold
  mv.push_back( { on_code( to ) , wu_deep , w_dn , cost ,
                  std::max( lo , th ) , hi , tag } );
  // a decrease of at most the gradient (or an increase)
  if( w_up >= - dg - 1e-9 )
   mv.push_back( { on_code( to ) , w_up , std::min( w_dn , dg ) , cost ,
                   lo , hi , tag } );
  };

 const bool can_count = ( f_mod_per_day < 0 ) ||
                        ( from.c < Index( f_mod_per_day ) );
 Label counted = from;
 if( f_nc > 1 )
  ++counted.c;

 // with the bands a modulation moves to an adjacent one, hence the two
 // directions never share a move, and the landing power of a stable
 // instant and of the last step of a modulation is that of a band
 const bool banded = ! f_bands.empty();
 const bool split = f_direction || banded;

 if( from.mode == 0 ) {                          // stable
  Label st = from;
  st.lk = from.lk ? from.lk - 1 : 0;
  emit( st , std::min( ru , mod_ramp_up[ t ] ) ,
        std::min( rd , mod_ramp_down[ t ] ) , 0.0 , 0 ,
        band_lo( from.b ) , band_hi( from.b ) );
  if( ( from.lk == 0 ) && can_count ) {          // start a modulation
   Label end = counted;
   end.mode = 0;
   end.lk = B;
   const bool can_up = ( ! banded ) || ( from.b + 1 < f_nband );
   const bool can_dn = ( ! banded ) || ( from.b > 0 );
   Label eu = end , ed = end;
   if( banded ) {
    eu.b = from.b + 1;
    ed.b = from.b ? from.b - 1 : 0;
    }
   if( ! split )                                 // the two merged
    emit( end , wu , wd , 0.0 , 1 );
   else {
    if( can_up )                                         // up, ends
     emit( eu , wu , 0.0 , 0.0 , 1 , band_lo( eu.b ) , band_hi( eu.b ) );
    if( can_dn )                                         // down, ends
     emit( ed , 0.0 , wd , cdn , 3 , band_lo( ed.b ) , band_hi( ed.b ) );
    }
   if( f_max_mod_length > 1 ) {                  // up / down, continues
    Label go = counted;
    go.lk = 1;
    if( can_up && ( fu <= ru + 1e-9 ) ) {
     go.mode = 1;
     emit( go , fu , - fu , 0.0 , 1 );
     }
    if( can_dn && ( fd <= rd + 1e-9 ) ) {
     go.mode = 2;
     emit( go , - fd , fd , cdn , 3 );
     }
    }
   }
  return;
  }

 // in the middle of a modulation: continue it or end it, the end landing
 // in the band next to the one the modulation left
 Label end = from;
 end.mode = 0;
 end.lk = B;
 if( banded )
  end.b = ( from.mode == 1 ) ? from.b + 1 : ( from.b ? from.b - 1 : 0 );
 Label go = from;
 ++go.lk;
 if( from.mode == 1 ) {                          // upward
  if( ( go.lk < f_max_mod_length ) && ( fu <= ru + 1e-9 ) )
   emit( go , fu , - fu , 0.0 , 1 );
  emit( end , wu , 0.0 , 0.0 , 1 , band_lo( end.b ) , band_hi( end.b ) );
  }
 else {                                          // downward
  if( ( go.lk < f_max_mod_length ) && ( fd <= rd + 1e-9 ) )
   emit( go , - fd , fd , cdn , 3 );
  emit( end , 0.0 , wd , cdn , 3 , band_lo( end.b ) , band_hi( end.b ) );
  }
 }

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

 if( f_direction ) {
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
 // at the instants with an on predecessor
 if( auto dd = b->get_deep_decrease() ) {
  auto ddrop = b->get_deep_drop();
  auto dlow = b->get_deep_low();
  for( Index i = 0 ; i < time_horizon ; ++i ) {
   const bool on_pred = built && U[ i ] &&
                        ( i ? bool( U[ i - 1 ] ) : ( init_up_down_time > 0 ) );
   const double prev = i ? P[ i - 1 ] : initial_power;
   const double tol = 1e-7 * std::max( 1.0 , std::abs( P[ i ] ) );
   dd[ i ].set_value( tag( i , 4 ) );
   ddrop[ i ].set_value( ( on_pred &&
                           ( prev - P[ i ] > f_deep_grad[ i ] + tol ) )
                         ? 1 : 0 );
   dlow[ i ].set_value( ( on_pred && ( P[ i ] < f_deep_thr[ i ] - tol ) )
                        ? 1 : 0 );
   }
  }

 // all the rest, the Block being already locked by this Solver
 ThermalUnitExtDPSolver::get_var_solution( solc );

 if( ! owned )
  f_Block->unlock( f_id );

 }  // end( NuclearUnitExtDPSolver::get_var_solution )

/*--------------------------------------------------------------------------*/
/*------------------ End File NuclearUnitExtDPSolver.cpp -------------------*/
/*--------------------------------------------------------------------------*/
