/*--------------------------------------------------------------------------*/
/*------------------- File NuclearUnitBlockForms.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the formulations of the operating rules of
 * NuclearUnitBlock other than the default one, i.e., F0, F2, F3a and F4
 * (the compact ones) and F5, F6 and F7 (the flows on the label graph of
 * NuclearRules); see NuclearUnitBlock::generate_abstract_constraints() for
 * their rows.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Antonio Frangioni, Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "LinearFunction.h"

#include "NuclearUnitBlock.h"

#include <algorithm>

#include <array>

#include <cmath>

#include <map>

#include <tuple>

#include <unordered_map>

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

using coeff_pair = LinearFunction::coeff_pair;

using v_coeff_pair = LinearFunction::v_coeff_pair;

using Index = Block::Index;

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

namespace {

/// a linear expression: a sum of coefficients times Variable, plus a
/// constant
struct Ex {
 v_coeff_pair t;
 double c = 0;

 Ex( void ) = default;

 explicit Ex( ColVariable * v , double a = 1 ) { t.emplace_back( v , a ); }

 Ex & add( ColVariable * v , double a = 1 ) {
  t.emplace_back( v , a );
  return( *this );
  }

 Ex & add( const Ex & e , double a = 1 ) {
  for( const auto & p : e.t )
   t.emplace_back( p.first , a * p.second );
  c += a * e.c;
  return( *this );
  }

 Ex & cst( double a ) {
  c += a;
  return( *this );
  }
 };

Ex operator+( Ex a , const Ex & b ) { return( std::move( a.add( b ) ) ); }

Ex operator-( Ex a , const Ex & b ) {
 return( std::move( a.add( b , -1 ) ) );
 }

Ex operator*( double k , Ex a ) {
 for( auto & p : a.t )
  p.second *= k;
 a.c *= k;
 return( a );
 }

/*--------------------------------------------------------------------------*/
/* The terms of e with the coefficients of the same Variable summed, in the
 * order of their first appearance; a coefficient that sums to 0 is kept, so
 * that the Variable of a row do not depend on the data. */

v_coeff_pair merged( v_coeff_pair && t )
{
 v_coeff_pair r;
 r.reserve( t.size() );
 if( t.size() <= 16 ) {
  for( auto & p : t ) {
   auto it = std::find_if( r.begin() , r.end() ,
                           [ & ]( const coeff_pair & q ) {
                            return( q.first == p.first ); } );
   if( it == r.end() )
    r.push_back( p );
   else
    it->second += p.second;
   }
  return( r );
  }
 std::unordered_map< ColVariable * , Index > pos;
 for( auto & p : t ) {
  auto it = pos.find( p.first );
  if( it == pos.end() ) {
   pos.emplace( p.first , r.size() );
   r.push_back( p );
   }
  else
   r[ it->second ].second += p.second;
  }
 return( r );
 }

/*--------------------------------------------------------------------------*/
/// the tolerance with which two outputs (or moves) are the same

double ptol( double x ) { return( 1e-9 * std::max( 1.0 , std::abs( x ) ) ); }

/// the tolerance with which a schedule satisfies a window or a range

double stol( double x ) { return( 1e-6 * std::max( 1.0 , std::abs( x ) ) ); }

/*--------------------------------------------------------------------------*/
/// a move of NuclearRules::on_moves(), as it fills the vector

struct Move {
 Index lab;
 double wup;
 double wdn;
 double cost;
 double lo;
 double hi;
 int tag;
 };

}  // end( namespace )

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS OF NuclearUnitBlock ----------------------*/
/*--------------------------------------------------------------------------*/

std::vector< std::pair< Index , Index > > NuclearUnitBlock::day_ranges(
                                                                void ) const
{
 std::vector< std::pair< Index , Index > > days;
 for( Index d0 = 0 ; d0 < f_time_horizon ; ) {
  const Index d1 = f_day_length ? std::min( d0 + f_day_length ,
                                            f_time_horizon )
                                : f_time_horizon;
  days.emplace_back( d0 , d1 );
  d0 = d1;
  }
 return( days );
 }

/*--------------------------------------------------------------------------*/

NuclearRules NuclearUnitBlock::rules_of( bool counters ) const
{
 NuclearRules r;
 r.time_horizon = f_time_horizon;
 r.mod_interval = Index( f_modulation_interval );
 r.init_modulation = f_initial_modulation;
 r.mod_ramp_up = v_modulation_ramp_up;
 r.mod_ramp_down = v_modulation_ramp_down;
 r.max_mod_length = std::max( f_max_modulation_length , Index( 1 ) );
 r.stab_start = f_stability_after_start;
 r.bands = v_power_bands;
 r.mod_per_day = counters ? f_modulations_per_day : -1;
 r.deep_per_day = counters ? f_deep_decreases_per_day : -1;
 r.starts_per_day = counters ? f_start_ups_per_day : -1;
 r.day_length = f_day_length;
 r.direction = has_modulation_direction();
 r.down_cost = v_down_modulation_cost;
 r.deep_thr = v_deep_threshold;
 r.deep_grad = v_deep_gradient;
 r.deep_cost = v_deep_cost;
 r.inf = Inf< double >();
 r.set_sizes();
 return( r );
 }

/*--------------------------------------------------------------------------*/

std::vector< FRowConstraint > & NuclearUnitBlock::form_group( void )
{
 if( generating_rows() ) {
  v_form_rows.emplace_back();
  return( v_form_rows.back() );
  }
 if( f_form_group >= v_form_rows.size() )
  throw( std::logic_error( "NuclearUnitBlock::update_rows: the formulation "
                           "of the rules would have more groups of rows" ) );
 return( v_form_rows[ f_form_group++ ] );
 }

/*--------------------------------------------------------------------------*/

void NuclearUnitBlock::add_form_rows( StagedRows && rows ,
                                      std::string && name )
{
 auto & g = form_group();
 if( generating_rows() )
  g.reserve( rows.size() );
 for( auto & r : rows )
  push_row( g , std::move( std::get< 0 >( r ) ) , std::get< 1 >( r ) ,
            std::get< 2 >( r ) );
 if( ! g.empty() )
  add_rows( g , std::move( name ) );
 }

/*--------------------------------------------------------------------------*/
/*--------------------- THE VARIABLE OF THE FORMULATIONS -------------------*/
/*--------------------------------------------------------------------------*/

void NuclearUnitBlock::generate_form_variables( void )
{
 const Index T = f_time_horizon;
 const int form = f_rules_form;
 const bool graph = form >= NuclearRules::F5Form;

 auto make = [ & ]( std::vector< ColVariable > & v , Index n ,
                    ColVariable::var_type type , const char * name ) {
  v.resize( n );
  for( auto & x : v )
   x.set_type( type );
  if( n )
   add_static_variable( v , name );
  };

 if( graph ) {
  // the flows on the label graph, which needs the unit to be committed
  // freely: the design variable would allow an empty path
  if( f_InvestmentCost != 0 )
   throw( std::logic_error( "NuclearUnitBlock::generate_abstract_variables: "
    "the formulations F5, F6 and F7 do not support the design variable" ) );
  build_label_graph( v_lg_nodes , v_lg_arcs , v_lg_q );
  make( v_arc , v_lg_arcs.size() , ColVariable::kBinary , "arc_nuclear" );
  if( form != NuclearRules::F7Form ) {
   Index nl = 0 , nd = 0;
   for( auto & a : v_lg_arcs ) {
    a.pland = ( ( a.kind == 0 ) || ( a.kind == 1 ) ) ? int( nl++ ) : -1;
    a.pdep = ( ( ( a.kind == 0 ) || ( a.kind == 2 ) ) &&
               ( v_lg_nodes[ a.from ].t >= 0 ) ) ? int( nd++ ) : -1;
    }
   make( v_arc_land , nl , ColVariable::kNonNegative ,
         "arc_power_nuclear" );
   make( v_arc_dep , nd , ColVariable::kNonNegative ,
         "arc_departure_nuclear" );
   }
  return;
  }

 // the compact formulations: the starts and the ends, by direction
 if( form != NuclearRules::F4Form ) {
  make( v_mod_start_up , T , ColVariable::kBinary , "m_start_up_nuclear" );
  make( v_mod_start_dn , T , ColVariable::kBinary , "m_start_dn_nuclear" );
  const Index ne = ( form == NuclearRules::F3aForm ) ? 1 : T;
  make( v_mod_end_up , ne , ColVariable::kBinary , "m_end_up_nuclear" );
  make( v_mod_end_dn , ne , ColVariable::kBinary , "m_end_dn_nuclear" );
  }
 else {
  // the modulations of F4: a direction, a first step a (not before the
  // instants in which a step is forbidden by the initial state), k steps,
  // the last at a + k - 1 <= T - 1; one that ends at T - 1 with at most
  // L^M - 1 steps may also be cut by the horizon
  Index a0 = 0;
  while( ( a0 < T ) && v_modulation[ a0 ].is_fixed() &&
         ( v_modulation[ a0 ].get_value() == 0 ) )
   ++a0;
  const Index L = std::max( f_max_modulation_length , Index( 1 ) );
  v_runs.clear();
  for( int dn = 0 ; dn < 2 ; ++dn )
   for( Index a = a0 ; a < T ; ++a )
    for( Index k = 1 ; ( k <= L ) && ( a + k <= T ) ; ++k ) {
     v_runs.push_back( { dn , a , k , false } );
     if( ( a + k == T ) && ( k < L ) )
      v_runs.push_back( { dn , a , k , true } );
     }
  make( v_run , v_runs.size() , ColVariable::kBinary , "run_nuclear" );
  }

 if( form != NuclearRules::F0Form ) {
  // the exact move: the move of a stable instant, the last steps, the
  // output at a start-up and before a shut-down
  make( v_stable_move , T , ColVariable::kContinuous ,
        "stable_move_nuclear" );
  make( v_last_up , T , ColVariable::kNonNegative , "last_step_up_nuclear" );
  make( v_last_dn , T , ColVariable::kNonNegative , "last_step_dn_nuclear" );
  make( v_start_power , T - init_t , ColVariable::kNonNegative ,
        "start_up_power_nuclear" );
  make( v_stop_power , T - init_t , ColVariable::kNonNegative ,
        "shut_down_power_nuclear" );

  // the split of the downward cases for the deep decreases, at the
  // instants with an on predecessor where a downward step may be deep
  v_split_t.clear();
  if( has_deep_decrease() && ( ! f_deep_f1 ) ) {
   const Index t0 = ( f_InitUpDownTime > 0 ) ? 0 : 1;
   for( Index t = t0 ; t < T ; ++t )
    if( v_deep_gradient[ t ] < v_DeltaRampDown[ t ] )
     v_split_t.push_back( t );
   const Index ns = v_split_t.size();
   make( v_full_deep , ns , ColVariable::kBinary , "full_deep_nuclear" );
   make( v_full_nodeep , ns , ColVariable::kBinary ,
         "full_not_deep_nuclear" );
   make( v_last_deep , ns , ColVariable::kBinary , "last_deep_nuclear" );
   make( v_last_high , ns , ColVariable::kBinary , "last_high_nuclear" );
   make( v_last_small , ns , ColVariable::kBinary , "last_small_nuclear" );
   make( v_xi_deep , ns , ColVariable::kNonNegative ,
         "last_deep_step_nuclear" );
   make( v_xi_high , ns , ColVariable::kNonNegative ,
         "last_high_step_nuclear" );
   make( v_xi_small , ns , ColVariable::kNonNegative ,
         "last_small_step_nuclear" );
   }
  }
 }  // end( NuclearUnitBlock::generate_form_variables )

/*--------------------------------------------------------------------------*/
/*----------------------- THE ROWS OF THE FORMULATIONS ---------------------*/
/*--------------------------------------------------------------------------*/

void NuclearUnitBlock::build_form_rows( void )
{
 f_form_group = 0;
 if( f_rules_form >= NuclearRules::F5Form )
  build_graph_rows();
 else
  build_compact_rows();
 }

/*--------------------------------------------------------------------------*/

void NuclearUnitBlock::build_compact_rows( void )
{
 static const std::string fn = "NuclearUnitBlock::build_rows";

 const Index T = f_time_horizon;
 const int iT = int( T );
 const Index L = std::max( f_max_modulation_length , Index( 1 ) );
 const int B = f_modulation_interval - 1;
 const int Bv = int( f_stability_after_start );
 const int it0 = int( init_t );
 const double INF = Inf< double >();
 const int form = f_rules_form;
 const bool f0 = ( form == NuclearRules::F0Form );
 const bool f3a = ( form == NuclearRules::F3aForm );
 const bool f4 = ( form == NuclearRules::F4Form );
 const auto days = day_ranges();

 // the exact move needs the downward stable move not to be a deep decrease
 if( ( ! f0 ) && has_deep_decrease() && ( ! f_deep_f1 ) )
  for( Index t = 0 ; t < T ; ++t )
   if( v_modulation_ramp_down[ t ] > v_deep_gradient[ t ] )
    throw( std::logic_error( fn + ": the formulations F2, F3a and F4 need "
     "ModulationDeltaRampDown not larger than DeepDecreaseGradient, which "
     "is not the case at time " + std::to_string( t ) ) );

 // append the row lhs <= e <= rhs, the constant of e moved to the sides
 auto row = [ & ]( StagedRows & g , Ex e ,
                   double lhs , double rhs ) {
  g.emplace_back( merged( std::move( e.t ) ) ,
                  ( lhs == - INF ) ? lhs : lhs - e.c ,
                  ( rhs == INF ) ? rhs : rhs - e.c );
  };

 // the natural variables, and the constants before the horizon
 auto u = [ & ]( int t ) -> Ex {
  if( t < 0 )
   return( Ex().cst( ( f_InitUpDownTime > 0 ) ? 1 : 0 ) );
  return( Ex( & v_commitment[ t ] ) );
  };
 auto v = [ & ]( int t ) -> Ex {
  return( ( t >= it0 ) && ( t < iT ) ? Ex( & v_start_up[ t - it0 ] )
                                     : Ex() );
  };
 auto w = [ & ]( int t ) -> Ex {
  return( ( t >= it0 ) && ( t < iT ) ? Ex( & v_shut_down[ t - it0 ] )
                                     : Ex() );
  };
 auto ut = [ & ]( int t ) { return( u( t ) - v( t ) ); };
 auto m = [ & ]( int t ) -> Ex {
  return( ( t < 0 ) || ( t >= iT ) ? Ex() : Ex( & v_modulation[ t ] ) );
  };
 auto dn = [ & ]( int t ) -> Ex {
  return( ( t < 0 ) || ( t >= iT ) ? Ex() : Ex( & v_modulation_down[ t ] ) );
  };
 // the steps in direction d, 0 upwards and 1 downwards
 auto md = [ & ]( int d , int t ) {
  return( d ? dn( t ) : m( t ) - dn( t ) );
  };
 auto p = [ & ]( int t ) -> Ex {
  if( t < 0 )
   return( Ex().cst( power_before() ) );
  return( Ex( & v_active_power[ t ] ) );
  };

 // the runs of F4 by their first step, by their end (those not cut), and
 // by the instants they contain
 std::vector< std::vector< Index > > run_start( 2 * T ) , run_end( 2 * T );
 std::vector< std::vector< Index > > run_in( 2 * T );
 if( f4 )
  for( Index r = 0 ; r < v_runs.size() ; ++r ) {
   const auto & z = v_runs[ r ];
   run_start[ z.dn * T + z.a ].push_back( r );
   if( ! z.cut )
    run_end[ z.dn * T + z.a + z.k - 1 ].push_back( r );
   for( Index h = z.a ; h < z.a + z.k ; ++h )
    run_in[ z.dn * T + h ].push_back( r );
   }

 // the starts s^d_t and the ends f^d_h, h = 1 , ... , T (the last step of
 // the modulation at h - 1)
 auto s = [ & ]( int d , int t ) -> Ex {
  if( ( t < 0 ) || ( t >= iT ) )
   return( Ex() );
  if( f4 ) {
   Ex e;
   for( auto r : run_start[ d * T + t ] )
    e.add( & v_run[ r ] );
   return( e );
   }
  return( Ex( d ? & v_mod_start_dn[ t ] : & v_mod_start_up[ t ] ) );
  };
 auto f = [ & ]( int d , int h ) -> Ex {
  if( ( h <= 0 ) || ( h > iT ) )
   return( Ex() );
  if( f4 ) {
   Ex e;
   for( auto r : run_end[ d * T + h - 1 ] )
    e.add( & v_run[ r ] );
   return( e );
   }
  if( f3a && ( h < iT ) )
   return( md( d , h - 1 ) - md( d , h ) + s( d , h ) );
  auto & fv = d ? v_mod_end_dn : v_mod_end_up;
  return( Ex( & fv[ f3a ? 0 : h - 1 ] ) );
  };

 // the steps: the flow of the starts and the ends, both within the steps -
 if( f4 ) {
  StagedRows gm;
  for( int t = 0 ; t < iT ; ++t ) {
   Ex em = m( t ) , ed = dn( t );
   for( int d = 0 ; d < 2 ; ++d )
    for( auto r : run_in[ d * T + t ] ) {
     em.add( & v_run[ r ] , -1 );
     if( d )
      ed.add( & v_run[ r ] , -1 );
     }
   row( gm , std::move( em ) , 0 , 0 );
   row( gm , std::move( ed ) , 0 , 0 );
   }
  add_form_rows( std::move( gm ) , "RunModulation_Nuclear" );

  // a modulation and the stability that follows it occupy the instants
  // a_r , ... , e_r + B
  StagedRows go;
  for( int t = 0 ; t < iT ; ++t ) {
   Ex e;
   for( Index r = 0 ; r < v_runs.size() ; ++r ) {
    const auto & z = v_runs[ r ];
    if( ( int( z.a ) <= t ) && ( t <= int( z.a + z.k ) - 1 + B ) )
     e.add( & v_run[ r ] );
    }
   if( ! e.t.empty() )
    row( go , std::move( e ) , -INF , 1 );
   }
  add_form_rows( std::move( go ) , "RunOccupation_Nuclear" );
  }
 else {
  if( ! f3a ) {
   StagedRows g;
   for( int d = 0 ; d < 2 ; ++d )
    for( int t = 0 ; t < iT ; ++t )
     row( g , md( d , t ) - md( d , t - 1 ) - s( d , t ) + f( d , t ) ,
          0 , 0 );
   add_form_rows( std::move( g ) , "ModulationFlow_Nuclear" );
   }

  // the starts and the ends are steps; in F3a an end is nonnegative
  StagedRows g;
  for( int d = 0 ; d < 2 ; ++d )
   for( int t = 0 ; t < iT ; ++t ) {
    row( g , s( d , t ) - md( d , t ) , -INF , 0 );
    if( f3a ) {
     if( t + 1 < iT )
      row( g , f( d , t + 1 ) , 0 , INF );
     else
      row( g , f( d , t + 1 ) - md( d , t ) , -INF , 0 );
     }
    else
     row( g , f( d , t + 1 ) - md( d , t ) , -INF , 0 );
    }
  add_form_rows( std::move( g ) , "ModulationStartEnd_Nuclear" );

  // the stability, a turn-off row on the ends of the last B instants
  StagedRows gs;
  for( int t = 1 ; t < iT ; ++t ) {
   Ex e = m( t );
   for( int h = std::max( 1 , t - B + 1 ) ; h <= t ; ++h )
    e.add( f( 0 , h ) ).add( f( 1 , h ) );
   row( gs , std::move( e ) , -INF , 1 );
   }
  add_form_rows( std::move( gs ) , "ModulationTurnOff_Nuclear" );

  // the longest modulation: each step belongs to a modulation started in
  // the last L^M instants, and one that the horizon cuts has L^M - 1 steps
  StagedRows gl;
  for( int d = 0 ; d < 2 ; ++d ) {
   for( int t = 0 ; t < iT ; ++t ) {
    Ex e = md( d , t );
    for( int h = std::max( 0 , t - int( L ) + 1 ) ; h <= t ; ++h )
     e.add( s( d , h ) , -1 );
    row( gl , std::move( e ) , -INF , 0 );
    }
   Ex e = md( d , iT - 1 ) - f( d , iT );
   for( int h = std::max( 0 , iT - int( L ) + 1 ) ; h < iT ; ++h )
    e.add( s( d , h ) , -1 );
   row( gl , std::move( e ) , -INF , 0 );
   }
  add_form_rows( std::move( gl ) , "ModulationLength_Nuclear" );
  }

 // the steps only when the unit stays on, upwards or downwards - - - - - -
 {
  StagedRows g;
  for( int t = 0 ; t < iT ; ++t ) {
   row( g , m( t ) - ut( t ) , -INF , 0 );
   row( g , dn( t ) - m( t ) , -INF , 0 );
   }
  add_form_rows( std::move( g ) , "ModulationOn_Nuclear" );
  }

 // the stability after a start-up- - - - - - - - - - - - - - - - - - - - -
 // a turn-on row if the stability is not longer than the minimum up time
 // (the unit is then on throughout it), else the rows of single instants
 if( Bv >= 2 ) {
  const bool agg = ( ! f_stab_single ) &&
                   ( Bv <= int( std::max( f_MinUpTime , Index( 1 ) ) ) );
  StagedRows g;
  for( int t = it0 ; t < iT ; ++t )
   if( agg ) {
    Ex e = m( t ) - u( t );
    for( int h = std::max( it0 , t - Bv + 1 ) ; h <= t ; ++h )
     e.add( v( h ) );
    row( g , std::move( e ) , -INF , 0 );
    }
   else
    for( int h = std::max( it0 , t - Bv + 1 ) ; h < t ; ++h )
     row( g , m( t ) + v( h ) , -INF , 1 );
  add_form_rows( std::move( g ) , "StartUpTurnOn_Nuclear" );
  }

 // the daily limits on the modulations and on the start-ups - - - - - - - -
 if( f_modulations_per_day >= 0 ) {
  StagedRows g;
  for( auto [ d0 , d1 ] : days ) {
   Ex e;
   for( Index t = d0 ; t < d1 ; ++t )
    e.add( s( 0 , t ) ).add( s( 1 , t ) );
   row( g , std::move( e ) , -INF , double( f_modulations_per_day ) );
   }
  add_form_rows( std::move( g ) , "ModulationsPerDay_Nuclear" );
  }

 if( f_start_ups_per_day >= 0 ) {
  StagedRows g;
  for( auto [ d0 , d1 ] : days ) {
   Ex e;
   for( Index t = std::max( d0 , init_t ) ; t < d1 ; ++t )
    e.add( v( t ) );
   if( ! e.t.empty() )
    row( g , std::move( e ) , -INF , double( f_start_ups_per_day ) );
   }
  add_form_rows( std::move( g ) , "StartUpsPerDay_Nuclear" );
  }

 // the move of the output - - - - - - - - - - - - - - - - - - - - - - - - -
 // the position of t among the instants of the deep split, if any
 std::vector< int > split( T , -1 );
 for( Index k = 0 ; k < v_split_t.size() ; ++k )
  split[ v_split_t[ k ] ] = int( k );

 if( f0 ) {
  // the move linear in the indicators of the state, with the maximum
  // output as the coefficient of the start-ups and shut-downs
  StagedRows g;
  for( int t = 0 ; t < iT ; ++t ) {
   const double Du = v_DeltaRampUp[ t ] , Dd = v_DeltaRampDown[ t ];
   const double Pv = get_operational_max_power( t );
   const double Pw = t ? get_operational_max_power( t - 1 )
                       : power_before();
   const Ex dp = p( t ) - p( t - 1 );
   const Ex st = ut( t ) - m( t );
   row( g , dp - Du * ( md( 0 , t ) - f( 0 , t + 1 ) ) +
            Dd * md( 1 , t ) + v_modulation_ramp_down[ t ] * st +
            Pw * w( t ) , 0 , INF );
   row( g , dp - Du * md( 0 , t ) + Dd * ( md( 1 , t ) - f( 1 , t + 1 ) ) -
            v_modulation_ramp_up[ t ] * st - Pv * v( t ) , -INF , 0 );
   }
  add_form_rows( std::move( g ) , "ModulationMoveSource_Nuclear" );
  }
 else {
  StagedRows ge;
  StagedRows gb;
  for( int t = 0 ; t < iT ; ++t ) {
   const double Du = v_DeltaRampUp[ t ] , Dd = v_DeltaRampDown[ t ];
   Ex e = p( t ) - p( t - 1 );
   e.add( & v_stable_move[ t ] , -1 );
   e.add( Du * ( md( 0 , t ) - f( 0 , t + 1 ) ) , -1 );
   e.add( & v_last_up[ t ] , -1 );
   e.add( Dd * ( md( 1 , t ) - f( 1 , t + 1 ) ) );
   e.add( & v_last_dn[ t ] );
   if( t >= it0 ) {
    e.add( & v_start_power[ t - it0 ] , -1 );
    e.add( & v_stop_power[ t - it0 ] );
    }
   row( ge , std::move( e ) , 0 , 0 );

   // the bounds of the move of each case
   const Ex st = ut( t ) - m( t );
   row( gb , Ex( & v_stable_move[ t ] ) -
             v_modulation_ramp_up[ t ] * st , -INF , 0 );
   row( gb , Ex( & v_stable_move[ t ] ) +
             v_modulation_ramp_down[ t ] * st , 0 , INF );
   row( gb , Ex( & v_last_up[ t ] ) - Du * f( 0 , t + 1 ) , -INF , 0 );
   if( split[ t ] < 0 )
    row( gb , Ex( & v_last_dn[ t ] ) - Dd * f( 1 , t + 1 ) , -INF , 0 );
   if( t >= it0 ) {
    const double lp = get_operational_min_power( t );
    const double lp1 = t ? get_operational_min_power( t - 1 )
                         : power_before();
    Ex su( & v_start_power[ t - it0 ] ) , sd( & v_stop_power[ t - it0 ] );
    row( gb , su - v_StartUpLimit[ t ] * v( t ) , -INF , 0 );
    row( gb , su - lp * v( t ) , 0 , INF );
    row( gb , sd - v_ShutDownLimit[ t ] * w( t ) , -INF , 0 );
    row( gb , sd - lp1 * w( t ) , 0 , INF );
    // the largest output, at a start-up the start-up limit (where the
    // deep decreases are split, the rows of the landing say it)
    if( split[ t ] < 0 )
     row( gb , p( t ) - v_StartUpLimit[ t ] * v( t ) -
               get_operational_max_power( t ) * ( u( t ) - v( t ) ) ,
          -INF , 0 );
    }
   }
  add_form_rows( std::move( ge ) , "ModulationMove_Nuclear" );
  add_form_rows( std::move( gb ) , "ModulationCases_Nuclear" );
  }

 // the reach of the output from a start-up and towards a shut-down - - - -
 // a unit that starts up at h with the output at most P^su_h has at t,
 // h <= t < h + tau^+ (on throughout), at most P^su_h plus the largest
 // increase that the rules allow from h to t; one that shuts down at k,
 // t < k <= t + tau^+ (on from t to k - 1), has at t at most P^sd_k plus
 // the largest decrease that the rules allow from t to k - 1; the window
 // holds at most one start-up (shut-down), hence the rows
 // p_t <= P^mx_t u_t - sum_h ( P^mx_t - psi_{h,t} ) v_h and the same in w
 if( f_reach ) {
  const int tp = int( std::max( f_MinUpTime , Index( 1 ) ) );
  // the largest move from the stable state of lockout lk0 after the
  // instant s0 to the instant s1, upwards or downwards, by the stability
  // and full ramps and the modulations of at most L^M steps
  auto reach = [ & ]( int s0 , int s1 , int lk0 , bool upw ) {
   return( rules_reach( s0 , s1 , lk0 , upw ) );
   };
  StagedRows g;
  for( int t = 0 ; t < iT ; ++t ) {
   const double up = get_operational_max_power( t );
   Ex a = p( t ) - up * u( t );
   bool any = false;
   for( int h = std::max( it0 , t - tp + 1 ) ; h <= t ; ++h ) {
    const double psi = std::min( up , v_StartUpLimit[ h ] +
                       reach( h , t , std::max( Bv - 1 , 0 ) , true ) );
    a.add( ( up - psi ) * v( h ) );
    any = true;
    }
   if( any )
    row( g , std::move( a ) , -INF , 0 );
   Ex b = p( t ) - up * u( t );
   any = false;
   for( int k = std::max( it0 , t + 1 ) ; k <= std::min( iT - 1 , t + tp ) ;
        ++k ) {
    const double psi = std::min( up , v_ShutDownLimit[ k ] +
                                      reach( t , k - 1 , 0 , false ) );
    b.add( ( up - psi ) * w( k ) );
    any = true;
    }
   if( any )
    row( g , std::move( b ) , -INF , 0 );
   }
  add_form_rows( std::move( g ) , "ModulationReach_Nuclear" );
  }

 // the deep decreases - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( has_deep_decrease() ) {
  if( f0 || f_deep_f1 )
   build_deep_rows_f1( days );
  else {
   // the downward cases split, where a downward step may be deep, and no
   // deep decrease elsewhere
   StagedRows gs;
   StagedRows gx;
   StagedRows gp;
   for( int t = 0 ; t < iT ; ++t ) {
    const int k = split[ t ];
    if( k < 0 ) {
     row( gs , Ex( & v_deep[ t ] ) , 0 , 0 );
     continue;
     }
    const double Dd = v_DeltaRampDown[ t ];
    const double dg = v_deep_gradient[ t ];
    const double pd = v_deep_threshold[ t ];
    const double lp = get_operational_min_power( t );
    const double up = get_operational_max_power( t );
    ColVariable * phd = & v_full_deep[ k ];
    ColVariable * phn = & v_full_nodeep[ k ];
    ColVariable * gd = & v_last_deep[ k ];
    ColVariable * gh = & v_last_high[ k ];
    ColVariable * gl = & v_last_small[ k ];
    row( gs , Ex( phd ).add( phn ) - md( 1 , t ) + f( 1 , t + 1 ) , 0 , 0 );
    row( gs , Ex( gd ).add( gh ).add( gl ) - f( 1 , t + 1 ) , 0 , 0 );
    row( gs , Ex( & v_deep[ t ] ).add( phd , -1 ).add( gd , -1 ) , 0 , 0 );

    row( gx , Ex( & v_last_dn[ t ] ).add( & v_xi_deep[ k ] , -1 ).add(
              & v_xi_high[ k ] , -1 ).add( & v_xi_small[ k ] , -1 ) , 0 , 0 );
    row( gx , Ex( & v_xi_small[ k ] ).add( gl , - dg ) , -INF , 0 );
    row( gx , Ex( & v_xi_deep[ k ] ).add( gd , - dg ) , 0 , INF );
    row( gx , Ex( & v_xi_deep[ k ] ).add( gd , - Dd ) , -INF , 0 );
    row( gx , Ex( & v_xi_high[ k ] ).add( gh , - dg ) , 0 , INF );
    row( gx , Ex( & v_xi_high[ k ] ).add( gh , - Dd ) , -INF , 0 );

    // the landing: below the threshold if deep, above it if a full step
    // or a large last step that is not deep
    row( gp , p( t ) - lp * u( t ) + ( lp - pd ) * Ex( phn ).add( gh ) ,
         0 , INF );
    row( gp , p( t ) + ( up - pd ) * Ex( phd ).add( gd ) -
              v_StartUpLimit[ t ] * v( t ) - up * ( u( t ) - v( t ) ) ,
         -INF , 0 );
    }
   add_form_rows( std::move( gs ) , "DeepSplit_Nuclear" );
   add_form_rows( std::move( gx ) , "DeepStep_Nuclear" );
   add_form_rows( std::move( gp ) , "DeepLanding_Nuclear" );

   if( f_deep_decreases_per_day >= 0 ) {
    StagedRows g;
    for( auto [ d0 , d1 ] : days ) {
     Ex e;
     for( Index t = d0 ; t < d1 ; ++t )
      e.add( & v_deep[ t ] );
     row( g , std::move( e ) , -INF , double( f_deep_decreases_per_day ) );
     }
    add_form_rows( std::move( g ) , "DeepDecreasesPerDay_Nuclear" );
    }
   }
  }

 // the bands: those of the default formulation, the last step of a
 // modulation being given by the ends- - - - - - - - - - - - - - - - - - -
 if( has_power_bands() ) {
  StagedRows g;
  for( int t = 0 ; t < iT ; ++t )
   row( g , Ex( & v_modulation_end[ t ] ) - f( 0 , t + 1 ) - f( 1 , t + 1 ) ,
        0 , 0 );
  add_form_rows( std::move( g ) , "ModulationEndDef_Nuclear" );
  build_band_rows( false );
  }

 }  // end( NuclearUnitBlock::build_compact_rows )

/*--------------------------------------------------------------------------*/
/*--------------------------- THE LABEL GRAPH ------------------------------*/
/*--------------------------------------------------------------------------*/

void NuclearUnitBlock::build_label_graph( std::vector< LGNode > & nodes ,
                                          std::vector< LGArc > & arcs ,
                                          std::vector< std::vector< double > >
                                                                    & Q ) const
{
 static const std::string fn = "NuclearUnitBlock::build_label_graph";

 const Index T = f_time_horizon;
 const int form = f_rules_form;
 const bool f7 = ( form == NuclearRules::F7Form );
 const auto R = rules_of( form != NuclearRules::F5Form );
 const Index tp = std::max( f_MinUpTime , Index( 1 ) );
 const Index tm = std::max( f_MinDownTime , Index( 1 ) );
 const double INF = Inf< double >();
 // the largest graph built: F7 is a reference for small instances
 const std::size_t max_arcs = f7 ? 2000000 : 6000000;

 nodes.clear();
 arcs.clear();
 Q.clear();

 // the outputs of F7: the anchors of each instant, moved by the endpoints
 // of the windows forwards and backwards along the horizon, within the
 // bounds of each instant
 if( f7 ) {
  auto lo = [ & ]( Index t ) { return( get_operational_min_power( t ) ); };
  auto hi = [ & ]( Index t ) { return( get_operational_max_power( t ) ); };
  auto anchors = [ & ]( Index j ) {
   std::vector< double > a = { lo( j ) , hi( j ) , v_StartUpLimit[ j ] };
   if( j + 1 < T )
    a.push_back( v_ShutDownLimit[ j + 1 ] );
   if( has_deep_decrease() )
    a.push_back( v_deep_threshold[ j ] );
   for( auto b : v_power_bands )
    a.push_back( b );
   return( a );
   };
  auto ends = [ & ]( Index h ) {
   std::vector< double > e = { 0 , v_DeltaRampUp[ h ] ,
                               - v_DeltaRampDown[ h ] ,
                               v_modulation_ramp_up[ h ] ,
                               - v_modulation_ramp_down[ h ] };
   if( has_deep_decrease() )
    e.push_back( - v_deep_gradient[ h ] );
   return( e );
   };
  // sorted, without values closer than the tolerance, within [ l , u ]
  auto clean = [ & ]( std::vector< double > & x , double l , double u ) {
   std::sort( x.begin() , x.end() );
   std::vector< double > y;
   for( auto z : x )
    if( ( z >= l - ptol( l ) ) && ( z <= u + ptol( u ) ) &&
        ( y.empty() || ( z > y.back() + ptol( z ) ) ) )
     y.push_back( std::min( std::max( z , l ) , u ) );
   x.swap( y );
   if( x.size() > 20000 )
    throw( std::logic_error( fn + ": F7 has too many outputs at an "
                             "instant (" + std::to_string( x.size() ) +
                             "), it is meant for small instances" ) );
   };
  std::vector< std::vector< double > > F( T ) , Bk( T );
  for( Index t = 0 ; t < T ; ++t ) {
   F[ t ] = anchors( t );
   const auto e = ends( t );
   std::vector< double > prev;
   if( t )
    prev = F[ t - 1 ];
   else
    if( f_InitUpDownTime > 0 )
     prev = { power_before() };
   for( auto x : prev )
    for( auto d : e )
     F[ t ].push_back( x + d );
   clean( F[ t ] , lo( t ) , hi( t ) );
   }
  for( Index t = T ; t-- > 0 ; ) {
   Bk[ t ] = anchors( t );
   if( t + 1 < T ) {
    const auto e = ends( t + 1 );
    for( auto x : Bk[ t + 1 ] )
     for( auto d : e )
      Bk[ t ].push_back( x - d );
    }
   clean( Bk[ t ] , lo( t ) , hi( t ) );
   }
  Q.resize( T );
  for( Index t = 0 ; t < T ; ++t ) {
   Q[ t ] = F[ t ];
   Q[ t ].insert( Q[ t ].end() , Bk[ t ].begin() , Bk[ t ].end() );
   clean( Q[ t ] , lo( t ) , hi( t ) );
   }
  }

 // the nodes of each instant, by ( on , tau , label , q )
 std::vector< std::map< std::tuple< bool , Index , Index , int > , Index > >
                                                                   layer( T );
 auto node = [ & ]( int t , bool on , Index tau , Index lab , int q ) {
  auto & L = layer[ t ];
  const auto key = std::make_tuple( on , tau , lab , q );
  auto it = L.find( key );
  if( it != L.end() )
   return( it->second );
  nodes.push_back( { t , on , tau , lab , q } );
  L.emplace( key , Index( nodes.size() - 1 ) );
  return( Index( nodes.size() - 1 ) );
  };

 // the source: the state before the horizon
 const double p_1 = power_before();
 const Index l0 = R.init_label( p_1 );
 if( f_InitUpDownTime > 0 )
  nodes.push_back( { -1 , true ,
                     std::min( Index( f_InitUpDownTime ) , tp ) , l0 , -1 } );
 else
  nodes.push_back( { -1 , false ,
                     std::min( Index( - f_InitUpDownTime ) , tm ) ,
                     R.shut_label( 0 , l0 ) , -1 } );

 auto add_arc = [ & ]( Index from , Index to , char kind , Index t , int tag ,
                       bool start , bool ends , Index band , double wup ,
                       double wdn , double lo , double hi ) {
  arcs.push_back( { from , to , kind , t , tag , start , ends , band , wup ,
                    wdn , lo , hi , -1 , -1 } );
  if( arcs.size() > max_arcs )
   throw( std::logic_error( fn + ": the label graph has more than " +
                            std::to_string( max_arcs ) + " arcs" ) );
  };

 // the output of a node, F7 only
 auto qval = [ & ]( Index n ) {
  const auto & x = nodes[ n ];
  return( x.t < 0 ? p_1 : Q[ x.t ][ x.q ] );
  };

 std::vector< Move > mv;
 std::vector< std::pair< Index , std::pair< double , double > > > ls;
 std::vector< Index > prev = { 0 };
 for( Index t = 0 ; t < T ; ++t ) {
  const double lp = get_operational_min_power( t );
  const double up = get_operational_max_power( t );
  const auto & Qt = f7 ? Q[ t ] : std::vector< double >();
  // the outputs of Q_t in [ a , b ], F7 only
  auto qrange = [ & ]( double a , double b ) {
   const auto i = std::lower_bound( Qt.begin() , Qt.end() , a - ptol( a ) );
   const auto j = std::upper_bound( Qt.begin() , Qt.end() , b + ptol( b ) );
   return( std::make_pair( int( i - Qt.begin() ) , int( j - Qt.begin() ) ) );
   };
  for( auto n : prev ) {
   const auto x = nodes[ n ];   // a copy: nodes may grow
   if( x.on ) {
    // the moves
    mv.clear();
    R.on_moves( t , x.lab , v_DeltaRampUp[ t ] , v_DeltaRampDown[ t ] , mv );
    const auto from = R.on_label( x.lab );
    for( const auto & m : mv ) {
     const auto to = R.on_label( m.lab );
     const bool start = ( m.tag & 1 ) && ( from.mode == 0 );
     const bool ends = ( m.tag & 1 ) && ( to.mode == 0 );
     const Index tau = std::min( x.tau + 1 , tp );
     if( ! f7 ) {
      add_arc( n , node( int( t ) , true , tau , m.lab , -1 ) , 0 , t ,
               m.tag , start , ends , to.b , m.wup , m.wdn , m.lo , m.hi );
      continue;
      }
     const double q = qval( n );
     const double a = std::max( { lp , m.lo , q - m.wdn } );
     const double b = std::min( { up , m.hi , q + m.wup } );
     const auto [ i , j ] = qrange( a , b );
     for( int k = i ; k < j ; ++k )
      add_arc( n , node( int( t ) , true , tau , m.lab , k ) , 0 , t ,
               m.tag , start , ends , to.b , m.wup , m.wdn , m.lo , m.hi );
     }
    // the shut-down, past the minimum up time and not during a modulation
    if( x.tau >= tp ) {
     const Index e = R.shut_label( t ? t - 1 : 0 , x.lab );
     if( ( e != NuclearRules::NO_LABEL ) &&
         ( ( ! f7 ) || ( qval( n ) <= v_ShutDownLimit[ t ] +
                                      ptol( v_ShutDownLimit[ t ] ) ) ) )
      add_arc( n , node( int( t ) , false , 1 , R.idle_label( t , e , 1 ) ,
                         -1 ) , 2 , t , 0 , false , false , 0 , 0 , 0 ,
               -INF , INF );
     }
    continue;
    }
   // an idle instant, and a restart past the minimum down time
   add_arc( n , node( int( t ) , false , std::min( x.tau + 1 , tm ) ,
                      R.idle_label( t , x.lab , 1 ) , -1 ) ,
            3 , t , 0 , false , false , 0 , 0 , 0 , -INF , INF );
   if( x.tau >= tm ) {
    R.start_labels( t , x.lab , ls );
    for( const auto & l : ls ) {
     const Index band = R.on_label( l.first ).b;
     if( ! f7 ) {
      add_arc( n , node( int( t ) , true , 1 , l.first , -1 ) , 1 , t , 0 ,
               false , false , band , 0 , 0 , l.second.first ,
               l.second.second );
      continue;
      }
     const double a = std::max( lp , l.second.first );
     const double b = std::min( v_StartUpLimit[ t ] , l.second.second );
     const auto [ i , j ] = qrange( a , b );
     for( int k = i ; k < j ; ++k )
      add_arc( n , node( int( t ) , true , 1 , l.first , k ) , 1 , t , 0 ,
               false , false , band , 0 , 0 , l.second.first ,
               l.second.second );
     }
    }
   }
  prev.clear();
  for( const auto & kv : layer[ t ] )
   prev.push_back( kv.second );
  }

 // the nodes that reach the last instant, the arcs between them
 std::vector< char > alive( nodes.size() , 0 );
 for( Index n = 0 ; n < nodes.size() ; ++n )
  if( nodes[ n ].t == int( T ) - 1 )
   alive[ n ] = 1;
 for( auto a = arcs.size() ; a-- > 0 ; )
  if( alive[ arcs[ a ].to ] )
   alive[ arcs[ a ].from ] = 1;
 // no path at all (F7 may have none, if no output fits): the graph is
 // kept whole, and its flow is infeasible as the problem is
 if( ! alive[ 0 ] )
  return;
 std::vector< Index > id( nodes.size() );
 std::vector< LGNode > nn;
 for( Index n = 0 ; n < nodes.size() ; ++n )
  if( alive[ n ] ) {
   id[ n ] = nn.size();
   nn.push_back( nodes[ n ] );
   }
 std::vector< LGArc > na;
 for( const auto & a : arcs )
  if( alive[ a.from ] && alive[ a.to ] ) {
   na.push_back( a );
   na.back().from = id[ a.from ];
   na.back().to = id[ a.to ];
   }
 nodes.swap( nn );
 arcs.swap( na );

 }  // end( NuclearUnitBlock::build_label_graph )

/*--------------------------------------------------------------------------*/

void NuclearUnitBlock::build_graph_rows( void )
{
 static const std::string fn = "NuclearUnitBlock::build_rows";

 const Index T = f_time_horizon;
 const int it0 = int( init_t );
 const int form = f_rules_form;
 const bool f7 = ( form == NuclearRules::F7Form );
 const double INF = Inf< double >();
 const double p_1 = power_before();

 // the graph of the data as they are: the one the Variable are built on
 // while the rows are generated, a new one, which must have the same
 // structure, while they are compared
 std::vector< LGNode > nodes_new;
 std::vector< LGArc > arcs_new;
 std::vector< std::vector< double > > Q_new;
 if( ! generating_rows() ) {
  build_label_graph( nodes_new , arcs_new , Q_new );
  bool same = ( nodes_new.size() == v_lg_nodes.size() ) &&
              ( arcs_new.size() == v_lg_arcs.size() ) &&
              ( Q_new == v_lg_q );
  for( Index n = 0 ; same && ( n < nodes_new.size() ) ; ++n ) {
   const auto & a = nodes_new[ n ];
   const auto & b = v_lg_nodes[ n ];
   same = ( a.t == b.t ) && ( a.on == b.on ) && ( a.tau == b.tau ) &&
          ( a.lab == b.lab ) && ( a.q == b.q );
   }
  for( Index k = 0 ; same && ( k < arcs_new.size() ) ; ++k ) {
   auto & a = arcs_new[ k ];
   const auto & b = v_lg_arcs[ k ];
   same = ( a.from == b.from ) && ( a.to == b.to ) && ( a.kind == b.kind ) &&
          ( a.tag == b.tag );
   a.pland = b.pland;
   a.pdep = b.pdep;
   }
  if( ! same )
   throw( std::logic_error( fn + ": the label graph of the formulation of "
                            "the rules would change, which cannot be done "
                            "in place" ) );
  }
 const auto & nodes = generating_rows() ? v_lg_nodes : nodes_new;
 const auto & arcs = generating_rows() ? v_lg_arcs : arcs_new;
 const auto & Q = generating_rows() ? v_lg_q : Q_new;

 auto row = [ & ]( StagedRows & g , Ex e ,
                   double lhs , double rhs ) {
  g.emplace_back( merged( std::move( e.t ) ) ,
                  ( lhs == - INF ) ? lhs : lhs - e.c ,
                  ( rhs == INF ) ? rhs : rhs - e.c );
  };

 // the arcs into and out of each node
 std::vector< std::vector< Index > > in( nodes.size() ) , out( nodes.size() );
 for( Index k = 0 ; k < arcs.size() ; ++k ) {
  out[ arcs[ k ].from ].push_back( k );
  in[ arcs[ k ].to ].push_back( k );
  }

 // the flow: one unit out of the source, conserved at every other node but
 // those of the last instant
 {
  StagedRows g;
  Ex e;
  for( auto k : out[ 0 ] )
   e.add( & v_arc[ k ] );
  row( g , std::move( e ) , 1 , 1 );
  for( Index n = 1 ; n < nodes.size() ; ++n ) {
   if( nodes[ n ].t + 1 >= int( T ) )
    continue;
   Ex c;
   for( auto k : in[ n ] )
    c.add( & v_arc[ k ] );
   for( auto k : out[ n ] )
    c.add( & v_arc[ k ] , -1 );
   row( g , std::move( c ) , 0 , 0 );
   }
  add_form_rows( std::move( g ) , "LabelFlow_Nuclear" );
  }

 // the natural variables as flows on the arcs of each instant
 {
  std::vector< Ex > eu( T ) , ev( T ) , ew( T ) , em( T ) , ed( T ) ,
                    edd( T ) , ee( T ) , ep( T );
  std::vector< std::array< Ex , 3 > > eb( T );
  std::vector< Ex > es( T );   // the starts of the modulations (F5)
  for( Index k = 0 ; k < arcs.size() ; ++k ) {
   const auto & a = arcs[ k ];
   ColVariable * y = & v_arc[ k ];
   const Index t = a.t;
   if( nodes[ a.to ].on ) {
    eu[ t ].add( y );
    if( has_power_bands() )
     eb[ t ][ a.band ].add( y );
    if( f7 )
     ep[ t ].add( y , Q[ t ][ nodes[ a.to ].q ] );
    else
     ep[ t ].add( & v_arc_land[ a.pland ] );
    }
   if( a.kind == 1 )
    ev[ t ].add( y );
   if( a.kind == 2 )
    ew[ t ].add( y );
   if( a.kind == 0 ) {
    if( a.tag & 1 )
     em[ t ].add( y );
    if( a.tag & 2 )
     ed[ t ].add( y );
    if( a.tag & 4 )
     edd[ t ].add( y );
    if( a.ends )
     ee[ t ].add( y );
    if( a.start )
     es[ t ].add( y );
    }
   }

  StagedRows g;
  for( Index t = 0 ; t < T ; ++t ) {
   row( g , Ex( & v_commitment[ t ] ) - eu[ t ] , 0 , 0 );
   if( int( t ) >= it0 ) {
    row( g , Ex( & v_start_up[ t - it0 ] ) - ev[ t ] , 0 , 0 );
    row( g , Ex( & v_shut_down[ t - it0 ] ) - ew[ t ] , 0 , 0 );
    }
   row( g , Ex( & v_modulation[ t ] ) - em[ t ] , 0 , 0 );
   if( ! v_modulation_down.empty() )
    row( g , Ex( & v_modulation_down[ t ] ) - ed[ t ] , 0 , 0 );
   if( ! v_deep.empty() )
    row( g , Ex( & v_deep[ t ] ) - edd[ t ] , 0 , 0 );
   if( has_power_bands() ) {
    for( Index b = 0 ; b < 3 ; ++b )
     row( g , Ex( & v_band[ b * T + t ] ) - eb[ t ][ b ] , 0 , 0 );
    row( g , Ex( & v_modulation_end[ t ] ) - ee[ t ] , 0 , 0 );
    }
   }
  add_form_rows( std::move( g ) , "LabelIndicators_Nuclear" );

  StagedRows gp;
  for( Index t = 0 ; t < T ; ++t )
   row( gp , Ex( & v_active_power[ t ] ) - ep[ t ] , 0 , 0 );
  add_form_rows( std::move( gp ) , "LabelPower_Nuclear" );

  // the daily limits, which the labels of F5 do not count
  if( form == NuclearRules::F5Form ) {
   const auto days = day_ranges();
   if( f_modulations_per_day >= 0 ) {
    StagedRows gd;
    for( auto [ d0 , d1 ] : days ) {
     Ex e;
     for( Index t = d0 ; t < d1 ; ++t )
      e.add( es[ t ] );
     row( gd , std::move( e ) , -INF , double( f_modulations_per_day ) );
     }
    add_form_rows( std::move( gd ) , "ModulationsPerDay_Nuclear" );
    }
   if( f_start_ups_per_day >= 0 ) {
    StagedRows gd;
    for( auto [ d0 , d1 ] : days ) {
     Ex e;
     for( Index t = std::max( d0 , init_t ) ; t < d1 ; ++t )
      e.add( & v_start_up[ t - it0 ] );
     if( ! e.t.empty() )
      row( gd , std::move( e ) , -INF , double( f_start_ups_per_day ) );
     }
    add_form_rows( std::move( gd ) , "StartUpsPerDay_Nuclear" );
    }
   if( ( ! v_deep.empty() ) && ( f_deep_decreases_per_day >= 0 ) ) {
    StagedRows gd;
    for( auto [ d0 , d1 ] : days ) {
     Ex e;
     for( Index t = d0 ; t < d1 ; ++t )
      e.add( & v_deep[ t ] );
     row( gd , std::move( e ) , -INF ,
          double( f_deep_decreases_per_day ) );
     }
    add_form_rows( std::move( gd ) , "DeepDecreasesPerDay_Nuclear" );
    }
   }
  }

 if( f7 )
  return;    // the outputs are those of the nodes

 // the power of the arcs: windows and ranges, and the conservation of the
 // output at the on-nodes
 StagedRows ga;
 for( Index k = 0 ; k < arcs.size() ; ++k ) {
  const auto & a = arcs[ k ];
  const Index t = a.t;
  ColVariable * y = & v_arc[ k ];
  const bool src = ( nodes[ a.from ].t < 0 );
  // the departure power: a Variable, or the initial power from the source
  auto q = [ & ]( void ) -> Ex {
   if( src )
    return( Ex( y , p_1 ) );
   return( Ex( & v_arc_dep[ a.pdep ] ) );
   };
  if( a.kind == 0 ) {
   const double lo = std::max( get_operational_min_power( t ) , a.lo );
   const double hi = std::min( get_operational_max_power( t ) , a.hi );
   Ex pl( & v_arc_land[ a.pland ] );
   row( ga , pl - q() - Ex( y , a.wup ) , -INF , 0 );
   row( ga , q() - pl - Ex( y , a.wdn ) , -INF , 0 );
   row( ga , pl - Ex( y , lo ) , 0 , INF );
   row( ga , pl - Ex( y , hi ) , -INF , 0 );
   }
  if( a.kind == 1 ) {
   const double lo = std::max( get_operational_min_power( t ) , a.lo );
   const double hi = std::min( v_StartUpLimit[ t ] , a.hi );
   Ex pl( & v_arc_land[ a.pland ] );
   row( ga , pl - Ex( y , lo ) , 0 , INF );
   row( ga , pl - Ex( y , hi ) , -INF , 0 );
   }
  if( a.kind == 2 )
   row( ga , q() - Ex( y , v_ShutDownLimit[ t ] ) , -INF , 0 );
  if( ( ( a.kind == 0 ) || ( a.kind == 2 ) ) && ( ! src ) ) {
   row( ga , q() - Ex( y , get_operational_min_power( t - 1 ) ) , 0 , INF );
   row( ga , q() - Ex( y , get_operational_max_power( t - 1 ) ) , -INF , 0 );
   }
  }
 add_form_rows( std::move( ga ) , "LabelArcPower_Nuclear" );

 StagedRows gn;
 for( Index n = 1 ; n < nodes.size() ; ++n ) {
  if( ( ! nodes[ n ].on ) || ( nodes[ n ].t + 1 >= int( T ) ) )
   continue;
  Ex e;
  for( auto k : in[ n ] )
   e.add( & v_arc_land[ arcs[ k ].pland ] );
  for( auto k : out[ n ] )
   e.add( & v_arc_dep[ arcs[ k ].pdep ] , -1 );
  row( gn , std::move( e ) , 0 , 0 );
  }
 add_form_rows( std::move( gn ) , "LabelNodePower_Nuclear" );

 }  // end( NuclearUnitBlock::build_graph_rows )

/*--------------------------------------------------------------------------*/
/*------------------------ THE SOLUTION OF A SCHEDULE ----------------------*/
/*--------------------------------------------------------------------------*/

std::vector< Index > NuclearUnitBlock::label_graph_path( void ) const
{
 const Index T = f_time_horizon;
 const auto & nodes = v_lg_nodes;
 const auto & arcs = v_lg_arcs;
 auto val = [ & ]( const std::vector< ColVariable > & v , Index i ) {
  return( v[ i ].get_value() );
  };
 auto bin = [ & ]( const std::vector< ColVariable > & v , Index i ) {
  return( ( ! v.empty() ) && ( v[ i ].get_value() > 0.5 ) );
  };
 auto P = [ & ]( int t ) {
  return( t < 0 ? power_before() : val( v_active_power , t ) );
  };

 // whether the arc agrees with the schedule
 auto agrees = [ & ]( const LGArc & a ) -> bool {
  const Index t = a.t;
  const bool on = bin( v_commitment , t );
  if( nodes[ a.to ].on != on )
   return( false );
  if( a.kind == 0 ) {
   if( ( ( a.tag & 1 ) != 0 ) != bin( v_modulation , t ) )
    return( false );
   if( ( ! v_modulation_down.empty() ) &&
       ( ( ( a.tag & 2 ) != 0 ) != bin( v_modulation_down , t ) ) )
    return( false );
   if( ( ! v_deep.empty() ) &&
       ( ( ( a.tag & 4 ) != 0 ) != bin( v_deep , t ) ) )
    return( false );
   if( ( ! v_modulation_end.empty() ) &&
       ( a.ends != bin( v_modulation_end , t ) ) )
    return( false );
   const double d = P( t ) - P( int( t ) - 1 );
   if( ( d > a.wup + stol( a.wup ) ) || ( - d > a.wdn + stol( a.wdn ) ) )
    return( false );
   if( ( P( t ) < a.lo - stol( a.lo ) ) || ( P( t ) > a.hi + stol( a.hi ) ) )
    return( false );
   }
  if( a.kind == 1 )
   if( ( P( t ) < a.lo - stol( a.lo ) ) || ( P( t ) > a.hi + stol( a.hi ) ) )
    return( false );
  if( a.kind == 2 )
   if( P( int( t ) - 1 ) > v_ShutDownLimit[ t ] +
                           stol( v_ShutDownLimit[ t ] ) )
    return( false );
  if( nodes[ a.to ].on && has_power_bands() &&
      ( ! bin( v_band , a.band * T + t ) ) )
   return( false );
  if( ( nodes[ a.to ].q >= 0 ) &&
      ( std::abs( v_lg_q[ t ][ nodes[ a.to ].q ] - P( t ) ) >
        stol( P( t ) ) ) )
   return( false );
  return( true );
  };

 // the nodes reachable along the schedule, with the arc they are reached by
 std::vector< int > by( nodes.size() , -2 );
 by[ 0 ] = -1;
 std::vector< std::vector< Index > > out( nodes.size() );
 for( Index k = 0 ; k < arcs.size() ; ++k )
  out[ arcs[ k ].from ].push_back( k );
 std::vector< Index > cur = { 0 } , nxt;
 for( Index t = 0 ; t < T ; ++t ) {
  nxt.clear();
  for( auto n : cur )
   for( auto k : out[ n ] )
    if( ( by[ arcs[ k ].to ] == -2 ) && agrees( arcs[ k ] ) ) {
     by[ arcs[ k ].to ] = int( k );
     nxt.push_back( arcs[ k ].to );
     }
  cur.swap( nxt );
  if( cur.empty() )
   return( std::vector< Index >() );
  }

 std::vector< Index > path( T );
 Index n = cur.front();
 for( Index t = T ; t-- > 0 ; ) {
  path[ t ] = Index( by[ n ] );
  n = arcs[ path[ t ] ].from;
  }
 return( path );
 }

/*--------------------------------------------------------------------------*/

void NuclearUnitBlock::set_form_solution( void )
{
 const int form = f_rules_form;
 if( ( form == NuclearRules::F1Form ) || v_modulation.empty() ||
     v_active_power.empty() )
  return;

 const Index T = f_time_horizon;
 const int iT = int( T );
 const int it0 = int( init_t );

 if( form >= NuclearRules::F5Form ) {
  // the path of the schedule, with its output on the arcs
  for( auto & x : v_arc )
   x.set_value( 0 );
  for( auto & x : v_arc_land )
   x.set_value( 0 );
  for( auto & x : v_arc_dep )
   x.set_value( 0 );
  for( auto k : label_graph_path() ) {
   const auto & a = v_lg_arcs[ k ];
   v_arc[ k ].set_value( 1 );
   if( a.pland >= 0 )
    v_arc_land[ a.pland ].set_value( v_active_power[ a.t ].get_value() );
   if( a.pdep >= 0 )
    v_arc_dep[ a.pdep ].set_value( v_active_power[ a.t - 1 ].get_value() );
   }
  return;
  }

 auto P = [ & ]( int t ) {
  return( t < 0 ? power_before() : v_active_power[ t ].get_value() );
  };
 auto on = [ & ]( int t ) {
  return( t < 0 ? ( f_InitUpDownTime > 0 )
                : ( v_commitment[ t ].get_value() > 0.5 ) );
  };

 // the direction of a step, from the move where it does not matter
 if( ! has_modulation_direction() )
  for( int t = 0 ; t < iT ; ++t )
   v_modulation_down[ t ].set_value(
    ( ( v_modulation[ t ].get_value() > 0.5 ) &&
      ( P( t ) < P( t - 1 ) - stol( P( t - 1 ) ) ) ) ? 1 : 0 );

 // the steps by direction, 0 upwards and 1 downwards
 std::vector< std::array< bool , 2 > > st( T );
 for( int t = 0 ; t < iT ; ++t ) {
  const bool m = v_modulation[ t ].get_value() > 0.5;
  const bool d = v_modulation_down[ t ].get_value() > 0.5;
  st[ t ] = { m && ( ! d ) , m && d };
  }
 auto mx = [ & ]( int d , int t ) {
  return( ( t >= 0 ) && ( t < iT ) && st[ t ][ d ] );
  };

 // the ends: within the horizon where the steps stop, at the last instant
 // by the band if any, else if the step is not a full ramp or the
 // modulation has L^M steps (a cut one, which moves by the full ramp, is
 // as feasible as the one that ends)
 std::vector< std::array< bool , 2 > > fe( T );   // f^d_{t+1}
 for( int t = 0 ; t < iT ; ++t )
  for( int d = 0 ; d < 2 ; ++d )
   fe[ t ][ d ] = mx( d , t ) && ( ! mx( d , t + 1 ) );
 for( int d = 0 ; d < 2 ; ++d ) {
  const int t = iT - 1;
  if( ! mx( d , t ) )
   continue;
  bool ended;
  if( ! v_modulation_end.empty() )
   ended = v_modulation_end[ t ].get_value() > 0.5;
  else {
   int k = 0;
   for( int h = t ; ( h >= 0 ) && mx( d , h ) ; --h )
    ++k;
   const double D = d ? v_DeltaRampDown[ t ] : v_DeltaRampUp[ t ];
   const double mv = d ? P( t - 1 ) - P( t ) : P( t ) - P( t - 1 );
   ended = ( k >= int( std::max( f_max_modulation_length , Index( 1 ) ) ) ) ||
           ( mv < D - stol( D ) );
   }
  fe[ t ][ d ] = ended;
  }
 auto fv = [ & ]( int d , int h ) {        // f^d_h, h = 1 , ... , T
  return( ( h >= 1 ) && ( h <= iT ) && fe[ h - 1 ][ d ] );
  };

 if( form == NuclearRules::F4Form ) {
  for( Index r = 0 ; r < v_runs.size() ; ++r ) {
   const auto & z = v_runs[ r ];
   const int d = z.dn;
   const int a = int( z.a ) , e = int( z.a + z.k ) - 1;
   bool is = ( ! mx( d , a - 1 ) );
   for( int h = a ; is && ( h <= e ) ; ++h )
    is = mx( d , h );
   if( is )
    is = z.cut ? ( ( e == iT - 1 ) && ( ! fv( d , iT ) ) )
               : ( fv( d , e + 1 ) );
   v_run[ r ].set_value( is ? 1 : 0 );
   }
  }
 else {
  for( int t = 0 ; t < iT ; ++t ) {
   v_mod_start_up[ t ].set_value( ( mx( 0 , t ) && ! mx( 0 , t - 1 ) ) ? 1
                                                                       : 0 );
   v_mod_start_dn[ t ].set_value( ( mx( 1 , t ) && ! mx( 1 , t - 1 ) ) ? 1
                                                                       : 0 );
   }
  if( form == NuclearRules::F3aForm ) {
   v_mod_end_up[ 0 ].set_value( fv( 0 , iT ) ? 1 : 0 );
   v_mod_end_dn[ 0 ].set_value( fv( 1 , iT ) ? 1 : 0 );
   }
  else
   for( int t = 0 ; t < iT ; ++t ) {
    v_mod_end_up[ t ].set_value( fe[ t ][ 0 ] ? 1 : 0 );
    v_mod_end_dn[ t ].set_value( fe[ t ][ 1 ] ? 1 : 0 );
    }
  }

 if( form == NuclearRules::F0Form )
  return;

 // the case of each instant gives the value of its own variable of the move
 std::vector< int > split( T , -1 );
 for( Index k = 0 ; k < v_split_t.size() ; ++k )
  split[ v_split_t[ k ] ] = int( k );
 for( int t = 0 ; t < iT ; ++t ) {
  const double d = P( t ) - P( t - 1 );
  double sg = 0 , xu = 0 , xd = 0 , su = 0 , sd = 0;
  const bool u1 = on( t - 1 ) , u0 = on( t );
  if( u0 && ( ! u1 ) )
   su = P( t );
  else
   if( u1 && ( ! u0 ) )
    sd = P( t - 1 );
   else
    if( u0 && u1 ) {
     if( ! ( mx( 0 , t ) || mx( 1 , t ) ) )
      sg = d;
     else
      if( fe[ t ][ 0 ] )
       xu = d;
      else
       if( fe[ t ][ 1 ] )
        xd = - d;
     }
  v_stable_move[ t ].set_value( sg );
  v_last_up[ t ].set_value( std::max( xu , 0.0 ) );
  v_last_dn[ t ].set_value( std::max( xd , 0.0 ) );
  if( t >= it0 ) {
   v_start_power[ t - it0 ].set_value( su );
   v_stop_power[ t - it0 ].set_value( sd );
   }

  const int k = split[ t ];
  if( k < 0 )
   continue;
  const bool dd = v_deep[ t ].get_value() > 0.5;
  const bool full = mx( 1 , t ) && ( ! fe[ t ][ 1 ] );
  const bool last = fe[ t ][ 1 ];
  // a last step that is not deep is a large decrease that lands above the
  // threshold or a small one, the two being both possible only when it
  // drops by the gradient exactly to above the threshold
  const bool high = last && ( ! dd ) &&
                    ( - d > v_deep_gradient[ t ] +
                            stol( v_deep_gradient[ t ] ) );
  v_full_deep[ k ].set_value( ( full && dd ) ? 1 : 0 );
  v_full_nodeep[ k ].set_value( ( full && ( ! dd ) ) ? 1 : 0 );
  v_last_deep[ k ].set_value( ( last && dd ) ? 1 : 0 );
  v_last_high[ k ].set_value( high ? 1 : 0 );
  v_last_small[ k ].set_value( ( last && ( ! dd ) && ( ! high ) ) ? 1 : 0 );
  const double xi = std::max( xd , 0.0 );
  v_xi_deep[ k ].set_value( ( last && dd ) ? xi : 0 );
  v_xi_high[ k ].set_value( high ? xi : 0 );
  v_xi_small[ k ].set_value( ( last && ( ! dd ) && ( ! high ) ) ? xi : 0 );
  }

 // the copies of the rules on the runs, if any
 set_copy_solution();

 }  // end( NuclearUnitBlock::set_form_solution )

/*--------------------------------------------------------------------------*/
/*------------------------ THE RULES ON THE RUNS ---------------------------*/
/*--------------------------------------------------------------------------*/

std::vector< double > NuclearUnitBlock::rules_reach_all( int s0 , int lk0 ,
                                                         bool upw ) const
{
 const int iT = int( f_time_horizon );
 const int B = f_modulation_interval - 1;
 const Index L = std::max( f_max_modulation_length , Index( 1 ) );
 const int nl = std::max( B , lk0 ) + 1;
 const double NEG = - Inf< double >();
 auto best = []( double & x , double v ) { x = std::max( x , v ); };

 // st[ l ] stable with lockout l, md[ j ] in a modulation of j steps
 std::vector< double > st( nl , NEG ) , md( L + 1 , NEG );
 st[ lk0 ] = 0;
 std::vector< double > r( 1 , 0.0 );
 r.reserve( iT - s0 );
 for( int s = s0 + 1 ; s < iT ; ++s ) {
  const double D = upw ? v_DeltaRampUp[ s ] : v_DeltaRampDown[ s ];
  const double DM = upw ? v_modulation_ramp_up[ s ]
                        : v_modulation_ramp_down[ s ];
  std::vector< double > ns( nl , NEG ) , nm( L + 1 , NEG );
  for( int l = 0 ; l < nl ; ++l ) {
   if( st[ l ] == NEG )
    continue;
   best( ns[ std::max( l - 1 , 0 ) ] , st[ l ] + DM );
   if( l == 0 ) {
    best( ns[ B ] , st[ l ] + D );
    if( L > 1 )
     best( nm[ 1 ] , st[ l ] + D );
    }
   }
  for( Index j = 1 ; j < L ; ++j ) {
   if( md[ j ] == NEG )
    continue;
   if( j + 1 < L )
    best( nm[ j + 1 ] , md[ j ] + D );
   best( ns[ B ] , md[ j ] + D );
   }
  st.swap( ns );
  md.swap( nm );
  double m = NEG;
  for( auto x : st )
   m = std::max( m , x );
  for( auto x : md )
   m = std::max( m , x );
  r.push_back( m );
  }
 return( r );
 }

/*--------------------------------------------------------------------------*/

double NuclearUnitBlock::rules_reach( int s0 , int s1 , int lk0 ,
                                      bool upw ) const
{
 return( rules_reach_all( s0 , lk0 , upw )[ s1 - s0 ] );
 }

/*--------------------------------------------------------------------------*/

void NuclearUnitBlock::generate_copy_variables( void )
{
 const Index T = f_time_horizon;
 const auto base = get_formulation();

 // the first and the last instant on of the run i: the run ( h , k ) is on
 // at h - 1 , ... , k - 1, the one in progress before the horizon (h = 0)
 // from 0; a run with no instant on (the shut-down at 0) has no copy
 auto first = [ & ]( Index i ) {
  return( v_Y_plus[ i ].first ? v_Y_plus[ i ].first - 1 : 0 );
  };
 auto empty = [ & ]( Index i ) { return( v_Y_plus[ i ].second == 0 ); };
 auto last = [ & ]( Index i ) {
  return( std::min( v_Y_plus[ i ].second , T ) - 1 );
  };

 v_copies.clear();
 auto add = [ & ]( std::vector< Index > && runs , int fam ) {
  if( runs.empty() )
   return;
  Index a = T , b = 0;
  for( auto i : runs ) {
   a = std::min( a , first( i ) );
   b = std::max( b , last( i ) );
   }
  v_copies.push_back( { std::move( runs ) , a , b , fam , 0 } );
  };

 if( base == DPForm ) {           // one copy per run
  for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
   if( ! empty( i ) )
    add( { i } , 0 );
  }
 if( ( base == SUForm ) || ( base == SUSDForm ) ) {  // one per start-up
  std::map< Index , std::vector< Index > > byh;
  for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
   if( ! empty( i ) )
    byh[ v_Y_plus[ i ].first ].push_back( i );
  for( auto & kv : byh )
   add( std::move( kv.second ) , 0 );
  }
 if( ( base == SDForm ) || ( base == SUSDForm ) ) {  // one per shut-down
  std::map< Index , std::vector< Index > > byk;
  for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
   if( ! empty( i ) )
    byk[ v_Y_plus[ i ].second ].push_back( i );
  for( auto & kv : byk )
   add( std::move( kv.second ) , 1 );
  }

 Index n = 0;
 for( auto & c : v_copies ) {
  c.v0 = n;
  n += c.b - c.a + 1;
  }

 auto make = [ & ]( std::vector< ColVariable > & v ,
                    ColVariable::var_type type , const char * name ) {
  v.resize( n );
  for( auto & x : v )
   x.set_type( type );
  add_static_variable( v , name );
  };
 make( v_cp_mod_up , ColVariable::kPosUnitary , "cp_m_up_nuclear" );
 make( v_cp_mod_dn , ColVariable::kPosUnitary , "cp_m_dn_nuclear" );
 make( v_cp_start_up , ColVariable::kPosUnitary , "cp_s_up_nuclear" );
 make( v_cp_start_dn , ColVariable::kPosUnitary , "cp_s_dn_nuclear" );
 make( v_cp_end_up , ColVariable::kPosUnitary , "cp_f_up_nuclear" );
 make( v_cp_end_dn , ColVariable::kPosUnitary , "cp_f_dn_nuclear" );
 make( v_cp_stable , ColVariable::kContinuous , "cp_stable_move_nuclear" );
 make( v_cp_last_up , ColVariable::kNonNegative , "cp_last_up_nuclear" );
 make( v_cp_last_dn , ColVariable::kNonNegative , "cp_last_dn_nuclear" );
 make( v_cp_start_power , ColVariable::kNonNegative ,
       "cp_start_up_power_nuclear" );
 make( v_cp_stop_power , ColVariable::kNonNegative ,
       "cp_shut_down_power_nuclear" );

 // the steps of the copies need no fixing where those of the unit are
 // fixed to 0: they are nonnegative and they sum to them [see (54) in
 // generate_abstract_constraints()]
 }  // end( NuclearUnitBlock::generate_copy_variables )

/*--------------------------------------------------------------------------*/

void NuclearUnitBlock::build_run_rows( void )
{
 if( ! ( f_run_reach || f_run_rules || f_band_persp ) )
  return;

 const Index T = f_time_horizon;
 const int iT = int( T );
 const double INF = Inf< double >();
 const auto base = get_formulation();
 const int B = f_modulation_interval - 1;
 const Index L = std::max( f_max_modulation_length , Index( 1 ) );
 const int Bv = int( f_stability_after_start );
 const int lock0 = std::max( 0 , std::min( f_modulation_interval -
                                           f_initial_modulation , iT ) );

 auto row = [ & ]( StagedRows & g , Ex e , double lhs , double rhs ) {
  g.emplace_back( merged( std::move( e.t ) ) ,
                  ( lhs == - INF ) ? lhs : lhs - e.c ,
                  ( rhs == INF ) ? rhs : rhs - e.c );
  };

 auto first = [ & ]( Index i ) {
  return( v_Y_plus[ i ].first ? v_Y_plus[ i ].first - 1 : 0 );
  };
 auto last = [ & ]( Index i ) {
  return( int( std::min( v_Y_plus[ i ].second , T ) ) - 1 );
  };
 auto y = [ & ]( Index i ) { return( & v_commitment_plus[ i ] ); };

 // the positions of the private outputs of the network formulations: of
 // the instant first( i ) of the run i in v_P_h_k (in which the instants of
 // a run are consecutive), of the instant 0 of the start-up h in v_P_h and
 // of the instant 0 of the shut-down k in v_P_k (the output of instant t
 // being in position t - h + 1 , respectively t, from there)
 std::vector< Index > off_hk( v_Y_plus.size() , 0 );
 {
  Index j = 0;
  for( Index i = 0 ; i < v_Y_plus.size() ; ++i ) {
   off_hk[ i ] = j;
   for( Index t = 0 ; t < T ; ++t )
    if( ( v_Y_plus[ i ].first <= t + 1 ) && ( t + 1 <= v_Y_plus[ i ].second ) )
     ++j;
   }
  }
 std::map< Index , Index > off_h , off_k;
 for( Index j = v_P_h.size() ; j-- > 0 ; )
  off_h[ v_P_h[ j ].second ] = j;
 for( Index j = v_P_k.size() ; j-- > 0 ; )
  off_k[ v_P_k[ j ].second ] = j;
 auto p_hk = [ & ]( Index i , Index t ) {
  return( & v_active_power_h_k[ off_hk[ i ] + t - first( i ) ] );
  };
 auto p_h = [ & ]( Index h , Index t ) {
  return( & v_active_power_h[ off_h[ h ] + t - ( h ? h - 1 : 0 ) ] );
  };
 auto p_k = [ & ]( Index k , Index t ) {
  return( & v_active_power_k[ off_k[ k ] + t ] );
  };

 // the caps of the runs by the reach of the rules - - - - - - - - - - - - -
 // a run started at s = h - 1 with the output at most P^su_s has at t at
 // most P^su_s plus the largest increase that the rules allow from s to t,
 // the first step not before s + max{ Bv , 1 } (the run in progress before
 // the horizon: p_{-1} plus the largest increase from -1, the first step not
 // before the lockout); one that shuts down at k <= T - 1 has at t at most
 // P^sd_k plus the largest decrease from t to k - 1, the first step not
 // before the same instant
 if( f_run_reach ) {
  const int q0 = lock0;                       // first step, run from before
  std::vector< std::vector< double > > Rup( T + 1 );
  std::map< std::pair< int , int > , std::vector< double > > Rdn;
  auto rup = [ & ]( int s , int t ) {        // from the start-up s (or -1)
   auto & r = Rup[ s + 1 ];
   if( r.empty() )
    r = rules_reach_all( s , s < 0 ? q0 : std::max( Bv - 1 , 0 ) , true );
   return( r[ t - s ] );
   };
  auto rdn = [ & ]( int t , int t1 , int lk ) {
   auto & r = Rdn[ { t , lk } ];
   if( r.empty() )
    r = rules_reach_all( t , lk , false );
   return( r[ t1 - t ] );
   };
  std::vector< std::vector< double > > Rmid( T );
  auto rmid = [ & ]( int s , int t ) {       // from an on-state at s
   auto & r = Rmid[ s ];
   if( r.empty() )
    r = rules_reach_all( s , 0 , true );
   return( r[ t - s ] );
   };
  const bool on0 = ( f_InitUpDownTime > 0 );
  auto psi = [ & ]( Index i , int t ) {
   const int h = int( v_Y_plus[ i ].first );
   const int k = int( v_Y_plus[ i ].second );
   double c = get_operational_max_power( t );
   int q;                                    // the first step of the run
   if( h > 0 ) {
    const int s = h - 1;
    q = s + std::max( Bv , 1 );
    c = std::min( c , v_StartUpLimit[ s ] + ( t > s ? rup( s , t ) : 0.0 ) );
    }
   else {
    q = q0;
    if( on0 )
     c = std::min( c , power_before() + rup( -1 , t ) );
    }
   // every maximum power met along the run, from an instant s at or after
   // the start-up to t, plus the largest increase that the rules allow from
   // s to t from any label (that of a stable instant with no lockout), and
   // from t to an instant s before the shut-down, plus the largest decrease
   // from t to s
   const int s0 = h > 0 ? int( h ) - 1 : 0;
   for( int s = s0 ; s < t ; ++s )
    c = std::min( c , get_operational_max_power( s ) + rmid( s , t ) );
   const int s1 = ( k < iT ) ? int( k ) - 1 : iT - 1;
   for( int s = t + 1 ; s <= s1 ; ++s )
    c = std::min( c , get_operational_max_power( s ) +
                      rdn( t , s , std::max( 0 , q - t - 1 ) ) );
   if( k < iT )
    c = std::min( c , v_ShutDownLimit[ k ] +
                      ( t < k - 1 ? rdn( t , k - 1 ,
                                         std::max( 0 , q - t - 1 ) ) : 0.0 ) );
   return( c );
   };

  StagedRows g;
  if( base == ptForm )
   for( int t = 0 ; t < iT ; ++t ) {
    Ex e( & v_active_power[ t ] );
    for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
     if( ( int( first( i ) ) <= t ) && ( t <= last( i ) ) )
      e.add( y( i ) , - psi( i , t ) );
    row( g , std::move( e ) , -INF , 0 );
    }
  if( base == DPForm )
   for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
    for( int t = int( first( i ) ) ; t <= last( i ) ; ++t )
     row( g , Ex( p_hk( i , t ) ).add( y( i ) , - psi( i , t ) ) , -INF , 0 );
  if( ( base == SUForm ) || ( base == SUSDForm ) )
   for( auto & kv : off_h ) {
    const Index h = kv.first;
    for( int t = h ? int( h ) - 1 : 0 ; t < iT ; ++t ) {
     Ex e( p_h( h , t ) );
     for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
      if( ( v_Y_plus[ i ].first == h ) && ( int( first( i ) ) <= t ) &&
          ( t <= last( i ) ) )
       e.add( y( i ) , - psi( i , t ) );
     row( g , std::move( e ) , -INF , 0 );
     }
    }
  if( ( base == SDForm ) || ( base == SUSDForm ) )
   for( auto & kv : off_k ) {
    const Index k = kv.first;
    for( int t = 0 ; t < int( std::min( k , T ) ) ; ++t ) {
     Ex e( p_k( k , t ) );
     for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
      if( ( v_Y_plus[ i ].second == k ) && ( int( first( i ) ) <= t ) &&
          ( t <= last( i ) ) )
       e.add( y( i ) , - psi( i , t ) );
     row( g , std::move( e ) , -INF , 0 );
     }
    }
  add_form_rows( std::move( g ) , "RunReach_Nuclear" );
  }

 // the rules of F2 on each copy [see RunRules] - - - - - - - - - - - - - - -
 if( f_run_rules ) {
  // the expressions of a copy c at t: its runs on at t (U), started at t
  // (V), on at t - 1 and t (Ut), shut down at t (W), all of them (One),
  // in progress before the horizon (I); and its private output at t
  auto sum = [ & ]( const RunCopy & c , auto pred ) {
   Ex e;
   for( auto i : c.runs )
    if( pred( i ) )
     e.add( y( i ) );
   return( e );
   };
  auto U = [ & ]( const RunCopy & c , int t ) {
   return( sum( c , [ & ]( Index i ) {
    return( ( int( first( i ) ) <= t ) && ( t <= last( i ) ) ); } ) );
   };
  auto V = [ & ]( const RunCopy & c , int t ) {
   return( sum( c , [ & ]( Index i ) {
    return( ( v_Y_plus[ i ].first > 0 ) &&
            ( int( v_Y_plus[ i ].first ) - 1 == t ) ); } ) );
   };
  auto W = [ & ]( const RunCopy & c , int t ) {
   return( sum( c , [ & ]( Index i ) {
    return( int( v_Y_plus[ i ].second ) == t ); } ) );
   };
  auto Ut = [ & ]( const RunCopy & c , int t ) { return( U( c , t ) -
                                                         V( c , t ) ); };
  auto One = [ & ]( const RunCopy & c ) {
   return( sum( c , [ & ]( Index ) { return( true ); } ) );
   };
  auto I0 = [ & ]( const RunCopy & c ) {
   return( sum( c , [ & ]( Index i ) {
    return( v_Y_plus[ i ].first == 0 ); } ) );
   };
  auto P = [ & ]( const RunCopy & c , int t ) -> Ex {
   if( t < int( c.a ) )
    return( Ex() );
   if( base == DPForm )
    return( Ex( p_hk( c.runs.front() , t ) ) );
   if( c.fam == 0 )
    return( Ex( p_h( v_Y_plus[ c.runs.front() ].first , t ) ) );
   return( Ex( p_k( v_Y_plus[ c.runs.front() ].second , t ) ) );
   };
  // the Variable of the copy at t, empty outside its instants; the ends at
  // h in position h - 1
  auto var = [ & ]( std::vector< ColVariable > & v , const RunCopy & c ,
                    int t ) -> Ex {
   if( ( t < int( c.a ) ) || ( t > int( c.b ) ) )
    return( Ex() );
   return( Ex( & v[ c.v0 + t - c.a ] ) );
   };
  auto m = [ & ]( int d , const RunCopy & c , int t ) {
   return( var( d ? v_cp_mod_dn : v_cp_mod_up , c , t ) );
   };
  auto s = [ & ]( int d , const RunCopy & c , int t ) {
   return( var( d ? v_cp_start_dn : v_cp_start_up , c , t ) );
   };
  auto f = [ & ]( int d , const RunCopy & c , int h ) {
   if( h <= int( c.a ) )
    return( Ex() );
   return( var( d ? v_cp_end_dn : v_cp_end_up , c , h - 1 ) );
   };

  StagedRows gf , gse , gs , gl , go , ge , gm , gb;
  for( const auto & c : v_copies ) {
   const int a = int( c.a ) , b = int( c.b );
   // the flow of the starts and the ends, both within the steps
   for( int d = 0 ; d < 2 ; ++d )
    for( int t = a ; t <= b ; ++t ) {
     row( gf , m( d , c , t ) - m( d , c , t - 1 ) - s( d , c , t ) +
               f( d , c , t ) , 0 , 0 );
     row( gse , s( d , c , t ) - m( d , c , t ) , -INF , 0 );
     row( gse , f( d , c , t + 1 ) - m( d , c , t ) , -INF , 0 );
     }
   // the stability, a turn-off row on the ends of the last B instants
   const Ex one = One( c );
   for( int t = a ; t <= b ; ++t ) {
    Ex e = m( 0 , c , t ) + m( 1 , c , t ) - one;
    for( int h = std::max( a + 1 , t - B + 1 ) ; h <= t ; ++h )
     e.add( f( 0 , c , h ) ).add( f( 1 , c , h ) );
    row( gs , std::move( e ) , -INF , 0 );
    }
   // the longest modulation, and the one that the horizon cuts
   for( int d = 0 ; d < 2 ; ++d ) {
    for( int t = a ; t <= b ; ++t ) {
     Ex e = m( d , c , t );
     for( int h = std::max( a , t - int( L ) + 1 ) ; h <= t ; ++h )
      e.add( s( d , c , h ) , -1 );
     row( gl , std::move( e ) , -INF , 0 );
     }
    if( b == iT - 1 ) {
     Ex e = m( d , c , b ) - f( d , c , iT );
     for( int h = std::max( a , iT - int( L ) + 1 ) ; h < iT ; ++h )
      e.add( s( d , c , h ) , -1 );
     row( gl , std::move( e ) , -INF , 0 );
     }
    else   // every run on at b shuts down at b + 1: the modulation ends
     row( gl , f( d , c , b + 1 ) - m( d , c , b ) , 0 , 0 );
    }
   // the steps only when the runs stay on, not within the stability after
   // their start-up
   for( int t = a ; t <= b ; ++t ) {
    Ex e = m( 0 , c , t ) + m( 1 , c , t ) - Ut( c , t );
    if( Bv >= 2 )
     for( auto i : c.runs ) {
      const int h = int( v_Y_plus[ i ].first );
      if( ( h > 0 ) && ( h - 1 < t ) && ( t < h - 1 + Bv ) &&
          ( t <= last( i ) ) )
       e.add( y( i ) );
      }
    row( go , std::move( e ) , -INF , 0 );
    }
   // the exact move of the output of the copy, and the bounds of its cases
   for( int t = a ; t <= b ; ++t ) {
    const double Du = v_DeltaRampUp[ t ] , Dd = v_DeltaRampDown[ t ];
    Ex e = P( c , t ) - P( c , t - 1 );
    if( t == 0 )
     e.add( I0( c ) , - power_before() );
    e.add( var( v_cp_stable , c , t ) , -1 );
    e.add( Du * ( m( 0 , c , t ) - f( 0 , c , t + 1 ) ) , -1 );
    e.add( var( v_cp_last_up , c , t ) , -1 );
    e.add( Dd * ( m( 1 , c , t ) - f( 1 , c , t + 1 ) ) );
    e.add( var( v_cp_last_dn , c , t ) );
    e.add( var( v_cp_start_power , c , t ) , -1 );
    e.add( var( v_cp_stop_power , c , t ) );
    row( ge , std::move( e ) , 0 , 0 );

    const Ex st = Ut( c , t ) - m( 0 , c , t ) - m( 1 , c , t );
    row( gb , var( v_cp_stable , c , t ) - v_modulation_ramp_up[ t ] * st ,
         -INF , 0 );
    row( gb , var( v_cp_stable , c , t ) + v_modulation_ramp_down[ t ] * st ,
         0 , INF );
    row( gb , var( v_cp_last_up , c , t ) - Du * f( 0 , c , t + 1 ) ,
         -INF , 0 );
    row( gb , var( v_cp_last_dn , c , t ) - Dd * f( 1 , c , t + 1 ) ,
         -INF , 0 );
    const Ex vc = V( c , t ) , wc = W( c , t );
    const double lp = get_operational_min_power( t );
    row( gb , var( v_cp_start_power , c , t ) - v_StartUpLimit[ t ] * vc ,
         -INF , 0 );
    row( gb , var( v_cp_start_power , c , t ) - lp * vc , 0 , INF );
    if( t == 0 )   // the shut-down at 0, from the initial power
     row( gb , var( v_cp_stop_power , c , t ) - power_before() * wc , 0 , 0 );
    else {
     row( gb , var( v_cp_stop_power , c , t ) - v_ShutDownLimit[ t ] * wc ,
          -INF , 0 );
     row( gb , var( v_cp_stop_power , c , t ) -
               get_operational_min_power( t - 1 ) * wc , 0 , INF );
     }
    }
   }
  add_form_rows( std::move( gf ) , "CopyFlow_Nuclear" );
  add_form_rows( std::move( gse ) , "CopyStartEnd_Nuclear" );
  add_form_rows( std::move( gs ) , "CopyTurnOff_Nuclear" );
  add_form_rows( std::move( gl ) , "CopyLength_Nuclear" );
  add_form_rows( std::move( go ) , "CopyOn_Nuclear" );
  add_form_rows( std::move( ge ) , "CopyMove_Nuclear" );
  add_form_rows( std::move( gb ) , "CopyCases_Nuclear" );

  // the steps, the starts and the ends of the unit are the sums of those
  // of the copies of each family
  for( int fam = 0 ; fam < 2 ; ++fam ) {
   if( std::none_of( v_copies.begin() , v_copies.end() ,
                     [ & ]( const RunCopy & c ) { return( c.fam == fam ); } ) )
    continue;
   StagedRows g;
   for( int t = 0 ; t < iT ; ++t ) {
    Ex eu = Ex( & v_modulation[ t ] ) - Ex( & v_modulation_down[ t ] );
    Ex ed( & v_modulation_down[ t ] );
    Ex su( & v_mod_start_up[ t ] ) , sd( & v_mod_start_dn[ t ] );
    Ex fu( & v_mod_end_up[ t ] ) , fd( & v_mod_end_dn[ t ] );
    for( const auto & c : v_copies )
     if( c.fam == fam ) {
      eu.add( m( 0 , c , t ) , -1 );
      ed.add( m( 1 , c , t ) , -1 );
      su.add( s( 0 , c , t ) , -1 );
      sd.add( s( 1 , c , t ) , -1 );
      fu.add( f( 0 , c , t + 1 ) , -1 );
      fd.add( f( 1 , c , t + 1 ) , -1 );
      }
    for( auto e : { & eu , & ed , & su , & sd , & fu , & fd } )
     row( g , std::move( *e ) , 0 , 0 );
    }
   add_form_rows( std::move( g ) , "CopyLink_Nuclear" );
   }
  }

 // the perspective cuts on the bands [see BandPersp] - - - - - - - - - - - -
 // for P in { P^mn_t , b_1 , b_2 , P^mx_t }, z_t >= 2 P p_t - sum_k c_k b^k_t
 // with c_k = 2 P x_k - x_k^2 and x_k the point of the band k nearest to P,
 // since p^2 >= 2 P p - c_k for every p in the band k
 if( f_band_persp && has_power_bands() ) {
  StagedRows g;
  for( int t = 0 ; t < iT ; ++t ) {
   const double lo[ 3 ] = { get_operational_min_power( t ) ,
                            v_power_bands[ 0 ] , v_power_bands[ 1 ] };
   const double hi[ 3 ] = { v_power_bands[ 0 ] , v_power_bands[ 1 ] ,
                            get_operational_max_power( t ) };
   for( double Pc : { lo[ 0 ] , lo[ 1 ] , lo[ 2 ] , hi[ 2 ] } ) {
    Ex e( & v_cut[ t ] );
    e.add( & v_active_power[ t ] , - 2 * Pc );
    for( Index k = 0 ; k < 3 ; ++k ) {
     const double x = std::clamp( Pc , lo[ k ] , std::max( lo[ k ] ,
                                                            hi[ k ] ) );
     e.add( & v_band[ k * T + t ] , 2 * Pc * x - x * x );
     }
    row( g , std::move( e ) , 0 , INF );
    }
   }
  add_form_rows( std::move( g ) , "BandPerspective_Nuclear" );
  }

 }  // end( NuclearUnitBlock::build_run_rows )

/*--------------------------------------------------------------------------*/

void NuclearUnitBlock::set_copy_solution( void )
{
 if( v_copies.empty() )
  return;
 const Index T = f_time_horizon;
 const int it0 = int( init_t );
 auto first = [ & ]( Index i ) {
  return( v_Y_plus[ i ].first ? v_Y_plus[ i ].first - 1 : 0 );
  };
 auto last = [ & ]( Index i ) {
  return( int( std::min( v_Y_plus[ i ].second , T ) ) - 1 );
  };
 for( auto * v : { & v_cp_mod_up , & v_cp_mod_dn , & v_cp_start_up ,
                   & v_cp_start_dn , & v_cp_end_up , & v_cp_end_dn ,
                   & v_cp_stable , & v_cp_last_up , & v_cp_last_dn ,
                   & v_cp_start_power , & v_cp_stop_power } )
  for( auto & x : *v )
   x.set_value( 0 );

 // the copy of each run of the path takes the values of the unit at the
 // instants of the run (the arcs being 0 or 1, see set_path_solution())
 for( const auto & c : v_copies )
  for( auto i : c.runs ) {
   if( v_commitment_plus[ i ].get_value() < 0.5 )
    continue;
   const int a = int( first( i ) ) , b = last( i );
   for( int t = a ; t <= b ; ++t ) {
    const Index j = c.v0 + t - c.a;
    const double mt = v_modulation[ t ].get_value();
    const double dt = v_modulation_down[ t ].get_value();
    v_cp_mod_up[ j ].set_value( mt - dt );
    v_cp_mod_dn[ j ].set_value( dt );
    v_cp_start_up[ j ].set_value( v_mod_start_up[ t ].get_value() );
    v_cp_start_dn[ j ].set_value( v_mod_start_dn[ t ].get_value() );
    v_cp_end_up[ j ].set_value( v_mod_end_up[ t ].get_value() );
    v_cp_end_dn[ j ].set_value( v_mod_end_dn[ t ].get_value() );
    v_cp_stable[ j ].set_value( v_stable_move[ t ].get_value() );
    v_cp_last_up[ j ].set_value( v_last_up[ t ].get_value() );
    v_cp_last_dn[ j ].set_value( v_last_dn[ t ].get_value() );
    if( ( t == a ) && ( t >= it0 ) && ( v_Y_plus[ i ].first > 0 ) )
     v_cp_start_power[ j ].set_value( v_start_power[ t - it0 ].get_value() );
    }
   // the shut-down of the run, if within the instants of the copy
   const int k = int( v_Y_plus[ i ].second );
   if( ( k < int( T ) ) && ( k <= int( c.b ) ) && ( k >= int( c.a ) ) &&
       ( k >= it0 ) )
    v_cp_stop_power[ c.v0 + k - c.a ].set_value(
                                     v_stop_power[ k - it0 ].get_value() );
   }
 }  // end( NuclearUnitBlock::set_copy_solution )

/*--------------------------------------------------------------------------*/
/*----------------- End File NuclearUnitBlockForms.cpp ---------------------*/
/*--------------------------------------------------------------------------*/
