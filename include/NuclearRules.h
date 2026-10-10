/*--------------------------------------------------------------------------*/
/*------------------------- File NuclearRules.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class NuclearRules, which holds the operating rules
 * of a NuclearUnitBlock [see NuclearUnitBlock.h] in the form that the
 * dynamic programming Solver of the unit [see NuclearUnitExtDPSolver.h]
 * and the formulations of the unit written on its label graph [see
 * NuclearUnitBlock::generate_abstract_constraints()] share: the codes of
 * the formulations of the rules, the data of the rules, the labels of the
 * states and the moves out of them.
 *
 * The rules are those of NuclearUnitBlock, in its notation: after the end
 * of a modulation no step happens for \f$ \tau^M - 1 \f$ instants, i.e.,
 * the stability is \f$ B = \tau^M - 1 \f$ instants (a description by
 * states, in which a modulation state is a stable instant followed by the
 * steps, has \f$ \tau^M = B + 2 \f$, see NuclearUnitBlock); a modulation
 * has at most \f$ L^M \f$ steps; no step happens in the \f$ \tau^v \f$
 * instants that begin with a start-up; the days limit the modulations, the
 * start-ups and the deep decreases; a deep decrease is a decrease by more
 * than \f$ \tilde\Delta_t \f$ to an output below \f$ \tilde p_t \f$; and,
 * with the bands, a modulation crosses exactly one boundary of its band,
 * the steps before the last one keeping the output in the band of origin.
 * The labels and the moves are described in NuclearUnitExtDPSolver.h
 * (table (1) there), which they implement.
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
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __NuclearRules
 #define __NuclearRules
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <limits>

#include <utility>

#include <vector>

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS NuclearRules ----------------------------*/
/*--------------------------------------------------------------------------*/
/// the operating rules of a nuclear unit: formulations, data, labels, moves
/** The NuclearRules class gathers what the dynamic programming Solver of a
 * NuclearUnitBlock and the formulations of its operating rules have in
 * common. Its constants are the codes, in the int Configuration that
 * selects the formulation of the unit, of the formulations of the rules
 * [see NuclearUnitBlock::generate_abstract_variables()]. Its fields are the
 * data of the rules, which whoever uses it fills from the unit, and its
 * methods give the labels of the states and the moves out of them, as in
 * table (1) of NuclearUnitExtDPSolver.h; the ramps of the unit and its
 * initial power are given to the methods that need them, so that a user
 * that keeps them up to date by itself needs no reload of the rules. */

class NuclearRules
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*------------------------------- TYPES ------------------------------------*/

 using Index = unsigned int;  ///< the index type, that of Block::Index

/*----------------------------- CONSTANTS ----------------------------------*/
/** @name The codes of the formulations of the operating rules
 *  The bits of the int Configuration that selects the formulation of a
 *  NuclearUnitBlock which concern its operating rules [see
 *  NuclearUnitBlock::generate_abstract_variables() and
 *  NuclearUnitBlock::generate_abstract_constraints()]; the bits 0 - 4
 *  select the formulation of the thermal part [see ThermalUnitBlock].
 *  @{ */

 /// bit 5: the tight rows of the default formulation
 static constexpr int TightRules = 32;

 /// bit 6: the full ramp of the default formulation with one coefficient
 /// per case
 static constexpr int TightRamp = 64;

 /// bit 7: the tight rows on the start Variable separated rather than
 /// written; with the formulations of the rules other than the default one,
 /// the rows on the stability and on the longest modulation that are not
 /// needed to describe the schedules
 static constexpr int TightCuts = 128;

 /// the mask of bits 8 - 10, the formulation of the rules
 static constexpr int FormMsk = 7 << 8;

 /// the default formulation of the rules ("F1", with the bits above "F1T")
 static constexpr int F1Form = 0;

 /// the formulation by the state indicators of the modulations ("F0")
 static constexpr int F0Form = 1 << 8;

 /// the modulation-commitment formulation ("F2")
 static constexpr int F2Form = 2 << 8;

 /// F2 with the ends of the modulations projected out ("F3a")
 static constexpr int F3aForm = 3 << 8;

 /// the formulation by the runs of the modulations ("F4")
 static constexpr int F4Form = 4 << 8;

 /// the flow on the label graph without the counters of the day ("F5")
 static constexpr int F5Form = 5 << 8;

 /// the flow on the label graph ("F6")
 static constexpr int F6Form = 6 << 8;

 /// the flow on the label graph with the output discretized ("F7")
 static constexpr int F7Form = 7 << 8;

 /// bit 11: the stability after a start-up by the rows of single
 /// instants also where the aggregated row would do (F2, F3a, F4)
 static constexpr int StartUpStabSingle = 2048;

 /// bit 12: the deep decreases by the rows of the default formulation
 /// rather than by the split of the downward cases (F2, F3a, F4)
 static constexpr int DeepByF1Rows = 4096;

 /// bit 13: the rows that bound the output by what the rules let it reach
 /// from a start-up and towards a shut-down (F0, F2, F3a, F4)
 static constexpr int ReachRows = 8192;

 /// the label that forbids a shut-down or a restart
 static constexpr Index NO_LABEL = Index( -1 );

/** @} ---------------------------------------------------------------------*/
/*-------------------------------- TYPES -----------------------------------*/

 /// a label, decoded: the mode \f$ \varpi \f$ (0 stable, 1 up, 2 down),
 /// the lockout or the steps \f$ n^{lk} \f$, the counters of the day
 /// \f$ n^M \f$, \f$ n^{dd} \f$, \f$ n^{su} \f$ and the band \f$ k \f$
 struct Label {
  int mode;
  Index lk;
  Index c , a , s;
  Index b{};   ///< the band of the output, or the one a modulation left
  };

/*------------------------------- METHODS ----------------------------------*/
/** @name The sizes and the encoding of the labels
 *  @{ */

 /// computes the sizes of the parts of the labels from the data
 /** To be called once the data are filled (and again whenever they
  * change): the (mode, lockout or steps) pairs, and the ranges of the
  * counters that are limited. */

 void set_sizes( void ) {
  ncore = lockout_max() + 1 + 2 * ( max_mod_length - 1 );
  nc = ( mod_per_day >= 0 ) ? Index( mod_per_day ) + 1 : 1;
  na = ( ( ! deep_thr.empty() ) && ( deep_per_day >= 0 ) )
       ? Index( deep_per_day ) + 1 : 1;
  nv = ( starts_per_day >= 0 ) ? Index( starts_per_day ) + 1 : 1;
  nband = bands.empty() ? 1 : 3;
  ncount = nband * nc * na * nv;
  }

 /// number of distinct labels of the on-states
 Index on_labels( void ) const { return( ncore * ncount ); }

 /// number of distinct labels of the off-states
 Index off_labels( void ) const { return( ( lockout_max() + 1 ) * ncount ); }

 /// the lockout a modulation leaves behind, \f$ \tau^M - 1 \f$
 Index mod_lockout( void ) const {
  return( std::max( mod_interval , Index( 2 ) ) - 1 );
  }

 /// the largest lockout a label may carry,
 /// \f$ \max\{ \tau^M , \tau^v \} - 1 \f$
 /** The largest value the lockout of a label can take: \f$ \tau^M - 1 \f$
  * after the end of a modulation [see mod_lockout()] and \f$ \tau^v - 1 \f$
  * after a start-up, \f$ \tau^v \f$ being the stability that follows one.
  * It is the range of the lockout in the encoding of the labels, and
  * nothing else. */

 Index lockout_max( void ) const {
  return( std::max( mod_lockout() + 1 , stab_start ) - 1 );
  }

 /// the band the output p belongs to, the lowest one at a breakpoint, 0 if
 /// the output is not banded (the bands are numbered 0, 1, 2 here and
 /// 1, 2, 3 in NuclearUnitBlock)
 Index band_of( double p ) const {
  if( bands.empty() )
   return( 0 );
  return( ( p <= bands[ 0 ] ) ? 0 : ( ( p <= bands[ 1 ] ) ? 1 : 2 ) );
  }

 /// the day of the time instant t
 Index day( Index t ) const { return( day_length ? t / day_length : 0 ); }

 /// encode an on-label
 Index on_code( const Label & l ) const {
  const Index core = ( l.mode == 0 ) ? l.lk :
   lockout_max() + l.lk + ( l.mode == 2 ? max_mod_length - 1 : 0 );
  return( core + ncore * count_code( l ) );
  }

 /// encode an off-label (always stable)
 Index off_code( const Label & l ) const {
  return( l.lk + ( lockout_max() + 1 ) * count_code( l ) );
  }

 /// the counters part of a label
 Index count_code( const Label & l ) const {
  return( l.b + nband * ( l.c + nc * ( l.a + na * l.s ) ) );
  }

 /// decode an on-label
 Label on_label( Index lab ) const {
  Label l;
  Index core = lab % ncore;
  Index cnt = lab / ncore;
  const Index B = lockout_max();
  if( core <= B ) {
   l.mode = 0;
   l.lk = core;
   }
  else {
   core -= B;                               // 1 .. 2 ( L^M - 1 )
   l.mode = ( core < max_mod_length ) ? 1 : 2;
   l.lk = ( l.mode == 1 ) ? core : core - ( max_mod_length - 1 );
   }
  l.b = cnt % nband;
  cnt /= nband;
  l.c = cnt % nc;
  cnt /= nc;
  l.a = cnt % na;
  l.s = cnt / na;
  return( l );
  }

 /// decode an off-label
 Label off_label( Index e ) const {
  Label l;
  const Index L = lockout_max() + 1;
  l.mode = 0;
  l.lk = e % L;
  Index cnt = e / L;
  l.b = cnt % nband;
  cnt /= nband;
  l.c = cnt % nc;
  cnt /= nc;
  l.a = cnt % na;
  l.s = cnt / na;
  return( l );
  }

/** @} ---------------------------------------------------------------------*/
/** @name The labels of the states and the moves out of them
 *  The rules of table (1) of NuclearUnitExtDPSolver.h. The label of an
 *  on-state at \f$ t \f$ is the one after the decision taken at \f$ t \f$,
 *  i.e., the one entering \f$ t + 1 \f$, and so is that of an off-state.
 *  @{ */

 /// the label of the state entering the time instant 0
 /** The stable state with the initial lockout, no count yet and the band
  * of the initial output \p p_init, as an on-code; where it is needed as
  * an off-state it is read through shut_label(), which drops the band and
  * encodes it as an off-code. The lockout the unit enters the horizon with
  * comes from the last modulation, hence from \f$ \tau^M \f$, and not from
  * the largest lockout a label may carry, which the stability after a
  * start-up may have made larger. */

 Index init_label( double p_init ) const {
  const Index L = mod_lockout() + 1;
  const Index im = ( init_modulation > 0 ) ? Index( init_modulation )
                                           : Index( 0 );
  Label l;
  l.mode = 0;
  l.lk = ( im < L ) ? L - im : Index( 0 );
  l.c = l.a = l.s = 0;
  l.b = band_of( p_init );
  return( on_code( l ) );
  }

 /// the label of the off-state of a unit shutting down after being on at
 /// \p t with label \p lab, NO_LABEL during a modulation
 Index shut_label( Index t , Index lab ) const {
  Label l = on_label( lab );
  if( l.mode != 0 )                    // no shut-down during a modulation
   return( NO_LABEL );
  l.b = 0;              // an off unit has no output, hence no band
  return( off_code( l ) );
  }

 /// the label of an off-state after \p k idle instants from \p t on
 Index idle_label( Index t , Index e , Index k ) const {
  Label l = off_label( e );
  l.lk = ( l.lk > k ) ? l.lk - k : 0;
  if( day( t + k ) != day( t ) )
   l.c = l.a = l.s = 0;
  return( off_code( l ) );
  }

 /// the label of the on-state of a unit restarting at \p t out of the
 /// off-state \p e, NO_LABEL if the restart is forbidden
 Index start_label( Index t , Index e ) const {
  if( on_forced( t ) )     // a modulation fixed at t: m_t <= 1 - v_t
   return( NO_LABEL );
  Label l = off_label( e );
  if( ( starts_per_day >= 0 ) && ( l.s >= Index( starts_per_day ) ) )
   return( NO_LABEL );                // the start-ups of the day are over
  if( nv > 1 )
   ++l.s;
  // no modulation at a start-up, and none for the A - 1 instants that
  // follow it either, which is the lockout the restart is born with
  l.lk = std::max( l.lk ? l.lk - 1 : Index( 0 ) ,
                   stab_start ? stab_start - 1 : Index( 0 ) );
  if( day( t + 1 ) != day( t ) )
   l.c = l.a = l.s = 0;
  return( on_code( l ) );
  }

 /// the labels of the on-states of a unit restarting at \p t, each with
 /// the range of its landing power: one per band if the output is banded
 void start_labels( Index t , Index e ,
                    std::vector< std::pair< Index ,
                                 std::pair< double , double > > > & ls )
  const {
  ls.clear();
  const Index lab = start_label( t , e );
  if( lab == NO_LABEL )
   return;
  if( bands.empty() ) {          // one label over the whole range
   ls.push_back( { lab , { - inf , inf } } );
   return;
   }
  // one label per band, each over the range of its own band: which band a
  // unit restarts in is decided by the power it restarts at
  Label l = on_label( lab );
  for( Index b = 0 ; b < nband ; ++b ) {
   l.b = b;
   ls.push_back( { on_code( l ) ,
                   { b ? bands[ b - 1 ] : - inf ,
                     ( b + 1 < nband ) ? bands[ b ] : inf } } );
   }
  }

 /// true if label \p a is at least as good as label \p b for the future
 bool label_dominates( Index a , Index b ) const {
  const Label la = on_label( a );
  const Label lb = on_label( b );
  // two labels of different bands are not comparable: the band says where
  // the output is, not how much history the unit carries
  return( ( la.mode == lb.mode ) && ( la.b == lb.b ) && ( la.lk <= lb.lk ) &&
          ( la.c <= lb.c ) && ( la.a <= lb.a ) && ( la.s <= lb.s ) );
  }

/*--------------------------------------------------------------------------*/
 /// the moves out of an on-state with label \p lab at \p t - 1
 /** Appends to \p mv the moves out of an on-state with label \p lab at
  * \p t - 1 (for \p t == 0, out of the initial state), as in table (1) of
  * NuclearUnitExtDPSolver.h. Each move is appended as the aggregate
  * { landing label , largest increase , largest decrease , cost ,
  * lowest landing power , highest landing power , tag } of the type MV;
  * the tag has bit 0 for a modulation step, bit 1 for a downward one and
  * bit 2 for a deep decrease. \p ru and \p rd are the ramps of the step
  * \f$ t - 1 \to t \f$, the full ramps of a modulation step and the ones
  * the window of a stable instant is intersected with. The fixed Variable
  * of the rules [see fix_mod] admit only the moves that agree with them. */

 template< class MV >
 void on_moves( Index t , Index lab , double ru , double rd ,
                std::vector< MV > & mv ) const {
  const double fu = ru;                    // the full ramps of the step,
  const double fd = rd;                    // which are the window of a
  const double wu = ru;                    // modulation step (wu, wd) and
  const double wd = rd;                    // the move of a step that is
                                           // not the last one (fu, fd)
  const double cdn = down_cost.empty() ? 0.0 : down_cost[ t ];
  const bool deep = ! deep_thr.empty();
  const Index B = mod_lockout();      // what a modulation leaves behind
  const bool newday = ( t + 1 < time_horizon ) && ( day( t + 1 ) != day( t ) );

  const Label from = on_label( lab );

  // a fixed Variable of the operating rules admits only the moves that agree
  // with it: bit 0 of the tag of a move says that it is a modulation step,
  // bit 1 that it goes downwards and bit 2 that it is a deep decrease
  auto agree = [ & ]( const std::vector< signed char > & f , bool what ) {
   return( f.empty() || ( f[ t ] < 0 ) || ( ( f[ t ] > 0 ) == what ) );
   };
  const bool no_deep = agree( fix_deep , false );
  const bool yes_deep = agree( fix_deep , true );

  // the range of the output in each band: the bands only exist if the two
  // breakpoints are there, otherwise there is the one range of everything
  auto band_lo = [ & ]( Index b ) {
   return( bands.empty() || ( b == 0 ) ? - inf : bands[ b - 1 ] );
   };
  auto band_hi = [ & ]( Index b ) {
   return( bands.empty() || ( b + 1 >= nband ) ? inf : bands[ b ] );
   };

  // append the move landing in label to, whose landing power is restricted
  // to [ lo , hi ], splitting it for the deep decrease if its window reaches
  // a decrease of the deep-decrease gradient
  auto emit = [ & ]( Label to , double w_up , double w_dn , double cost ,
                     int tag , double lo , double hi ) {
   if( newday )
    to.c = to.a = to.s = 0;
   if( ( ! agree( fix_mod , tag & 1 ) ) ||
       ( ! agree( fix_down , tag & 2 ) ) )
    return;
   if( ( ! deep ) || ( w_dn < deep_grad[ t ] - 1e-9 ) ) {
    if( no_deep )
     mv.push_back( { on_code( to ) , w_up , w_dn , cost , lo , hi , tag } );
    return;
    }
   const double dg = deep_grad[ t ];
   const double th = deep_thr[ t ];
   const double wu_deep = std::min( w_up , - dg );
   // a decrease of at least the gradient to at most the threshold: deep
   if( yes_deep && ( ( deep_per_day < 0 ) ||
                     ( from.a < Index( deep_per_day ) ) ) ) {
    Label td = to;
    if( ( na > 1 ) && ( ! newday ) )
     ++td.a;
    mv.push_back( { on_code( td ) , wu_deep , w_dn ,
                    cost + ( deep_cost.empty() ? 0.0 : deep_cost[ t ] ) ,
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

  const bool can_count = ( mod_per_day < 0 ) ||
                         ( from.c < Index( mod_per_day ) );
  Label counted = from;
  if( nc > 1 )
   ++counted.c;

  // with the bands a modulation moves to an adjacent one, hence the two
  // directions never share a move, and the landing power of a stable
  // instant and of the last step of a modulation is that of a band
  const bool banded = ! bands.empty();
  const bool split = direction || banded;

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
    const bool can_up = ( ! banded ) || ( from.b + 1 < nband );
    const bool can_dn = ( ! banded ) || ( from.b > 0 );
    Label eu = end , ed = end;
    if( banded ) {
     eu.b = from.b + 1;
     ed.b = from.b ? from.b - 1 : 0;
     }
    if( ! split )                                 // the two merged
     emit( end , wu , wd , 0.0 , 1 , - inf , inf );
    else {
     if( can_up )                                         // up, ends
      emit( eu , wu , 0.0 , 0.0 , 1 , band_lo( eu.b ) , band_hi( eu.b ) );
     if( can_dn )                                         // down, ends
      emit( ed , 0.0 , wd , cdn , 3 , band_lo( ed.b ) , band_hi( ed.b ) );
     }
    if( max_mod_length > 1 ) {                  // up / down, continues
     // a step that is not the last one keeps the output in the band the
     // modulation starts from, which the label carries until it lands
     Label go = counted;
     go.lk = 1;
     if( can_up ) {
      go.mode = 1;
      emit( go , fu , - fu , 0.0 , 1 , band_lo( from.b ) ,
            band_hi( from.b ) );
      }
     if( can_dn ) {
      go.mode = 2;
      emit( go , - fd , fd , cdn , 3 , band_lo( from.b ) ,
            band_hi( from.b ) );
      }
     }
    }
   return;
   }

  // in the middle of a modulation: continue it, the output staying in the
  // band the modulation left (which the label carries), or end it, landing
  // in the band next to that one
  Label end = from;
  end.mode = 0;
  end.lk = B;
  if( banded )
   end.b = ( from.mode == 1 ) ? from.b + 1 : ( from.b ? from.b - 1 : 0 );
  Label go = from;
  ++go.lk;
  if( from.mode == 1 ) {                          // upward
   if( go.lk < max_mod_length )
    emit( go , fu , - fu , 0.0 , 1 , band_lo( from.b ) , band_hi( from.b ) );
   emit( end , wu , 0.0 , 0.0 , 1 , band_lo( end.b ) , band_hi( end.b ) );
   }
  else {                                          // downward
   if( go.lk < max_mod_length )
    emit( go , - fd , fd , cdn , 3 , band_lo( from.b ) , band_hi( from.b ) );
   emit( end , 0.0 , wd , cdn , 3 , band_lo( end.b ) , band_hi( end.b ) );
   }
  }

/*--------------------------------------------------------------------------*/
 /// true if a Variable of the operating rules is fixed to 1 at t, which
 /// needs the unit on at t and not starting up at t
 bool on_forced( Index t ) const {
  auto one = [ t ]( const std::vector< signed char > & f ) {
   return( ( ! f.empty() ) && ( f[ t ] == 1 ) );
   };
  return( one( fix_mod ) || one( fix_down ) || one( fix_deep ) );
  }

/** @} ---------------------------------------------------------------------*/
/*------------------------------- FIELDS -----------------------------------*/
/** @name The data of the rules
 *  Filled by whoever uses the rules from the NuclearUnitBlock [see its
 *  methods of the same name], after which set_sizes() is called.
 *  @{ */

 Index time_horizon{};                 ///< the time horizon
 Index mod_interval{ 2 };              ///< \f$ \tau^M \f$ (ModulationTime)
 int init_modulation{ 2 };             ///< \f$ \tau^M_0 \f$ (InitModulation)
 std::vector< double > mod_ramp_up;    ///< \f$ \Delta^{M+}_t \f$
 std::vector< double > mod_ramp_down;  ///< \f$ \Delta^{M-}_t \f$

 Index max_mod_length{ 1 };            ///< \f$ L^M \f$
 int mod_per_day{ -1 };                ///< \f$ N^M \f$, -1 if unlimited
 int deep_per_day{ -1 };               ///< \f$ N^{dd} \f$, -1 if unlimited
 int starts_per_day{ -1 };             ///< \f$ N^{su} \f$, -1 if unlimited
 Index day_length{ 0 };                ///< \f$ T^{day} \f$, 0 = horizon
 bool direction{ false };              ///< true if the direction matters
 std::vector< double > down_cost;      ///< \f$ c^-_t \f$ (empty = 0)
 std::vector< double > deep_thr;       ///< \f$ \tilde p_t \f$ (empty = none)
 std::vector< double > deep_grad;      ///< \f$ \tilde\Delta_t \f$
 std::vector< double > deep_cost;      ///< \f$ c^{dd}_t \f$ (empty = 0)

 /// \f$ \tau^v \f$, the instants of stability that begin with a start-up
 Index stab_start{};

 /// the two breakpoints that split the output into bands, empty if there
 /// are no bands [see NuclearUnitBlock::get_power_bands()]
 std::vector< double > bands;

 /// the value of an infinite power, that of the user of the rules
 double inf = std::numeric_limits< double >::infinity();

 /// the fixed Variable \f$ m_t \f$, \f$ d_t \f$ and \f$ \delta_t \f$ of the
 /// operating rules: -1 where free, 0 or 1 where fixed, empty if none of
 /// them is fixed
 std::vector< signed char > fix_mod;
 std::vector< signed char > fix_down;  ///< as fix_mod, for \f$ d_t \f$
 std::vector< signed char > fix_deep;  ///< the same, for \f$ \delta_t \f$

/** @} ---------------------------------------------------------------------*/
/** @name The sizes of the labels [see set_sizes()]
 *  @{ */

 Index ncore{ 2 };   ///< number of (mode, lockout or steps) pairs
 Index nc{ 1 };      ///< range of the counter of the modulations
 Index na{ 1 };      ///< range of the counter of the deep decreases
 Index nv{ 1 };      ///< range of the counter of the start-ups
 Index nband{ 1 };   ///< the number of bands: 3 with the breakpoints, 1 not
 Index ncount{ 1 };  ///< nband * nc * na * nv

/** @} ---------------------------------------------------------------------*/

 };  // end( class( NuclearRules ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif  /* NuclearRules.h included */

/*--------------------------------------------------------------------------*/
/*------------------------- End File NuclearRules.h ------------------------*/
/*--------------------------------------------------------------------------*/
