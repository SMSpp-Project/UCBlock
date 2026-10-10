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
 for( Index t = 0 ; t < time_horizon ; ++t )
  forced |= f_rules.on_forced( t );
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
