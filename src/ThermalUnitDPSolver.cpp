/*--------------------------------------------------------------------------*/
/*----------------------- File ThermalUnitDPSolver.cpp ---------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the EDPSolver class.
 *
 * \author Claudio Gentile \n
 *         Istituto di Analisi di Sistemi e Informatica "Antonio Ruberti" \n
 *         Consiglio Nazionale delle Ricerche \n
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Niccolo' Iardella \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; Claudio Gentile, Antonio Frangioni, Niccolo' Iardella
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "ThermalUnitDPSolver.h"

#include "ThermalUnitBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

SMSpp_insert_in_factory_cpp_0( ThermalUnitDPSolver );

/*--------------------------------------------------------------------------*/
/*--------------------------- Solver INTERFACE -----------------------------*/
/*--------------------------------------------------------------------------*/

void ThermalUnitDPSolver::set_Block( Block * block )
{
 if( block == f_Block )
  return;

 Solver::set_Block( block );

 if( block ) {
  if( ! dynamic_cast< ThermalUnitBlock * >( f_Block ) )
   throw( std::runtime_error( "ThermalUnitDPSolver requires ThermalUnitBlock"
			      ) );
  load_parameters();
  }
 }

/*--------------------------------------------------------------------------*/

int ThermalUnitDPSolver::compute( bool changedvars )
{
 process_modifications();

 switch( stage ) {
  case( start ):    build_graph();
  case( graph_OK ): compute_EDPs();
  case( edps_OK ):  min_path();
  case( path_OK ):  compute_solutions();
  }

 assert( stage == sol_OK );
 return( f_end.lab == TUDPINF ? kInfeasible : kOK );
 }

/*--------------------------------------------------------------------------*/

void ThermalUnitDPSolver::get_var_solution( Configuration * solc )
{
 // lock the block
 bool owned = f_Block->is_owned_by( f_id );
 if( ( ! owned ) && ( ! f_Block->lock( f_id ) ) )
  throw( std::runtime_error( "Unable to lock the Block" ) );

 auto b = static_cast< ThermalUnitBlock * >( f_Block );

 // generate abstract representation if necessary
 b->generate_abstract_variables( nullptr );
 // b->generate_objective( nullptr );

 // set active power and unit commitment variables
 auto pow_it = b->get_active_power( 0 );
 auto com_it = b->get_commitment( 0 );

 for( Index i = 0 ; i < time_horizon ; ++i , ++pow_it , ++com_it ) {
  pow_it->set_value( P[ i ] );
  com_it->set_value( U[ i ] ? 1 : 0 );
  }

 // unlock the block
 if( ! owned )
  f_Block->unlock( f_id );

 }  // end( ThermalUnitDPSolver::get_var_solution( )

/*--------------------------------------------------------------------------*/
/*------------------ BUILDING AND SOLVING THE DP PROBLEM -------------------*/
/*--------------------------------------------------------------------------*/

void ThermalUnitDPSolver::build_graph( void )
{
 // first reset any existing graph; note that node labels will be set to 0,
 // which we use as a way to indicate that the node has not been proved
 // reachable from s yet

 delete f_start.DPS;
 
 v_on_nodes.clear();  // this deletes all EDSolver
 v_on_nodes.resize( time_horizon );
 v_off_nodes.resize( time_horizon );

 // now rebuild the graph: start from s- - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
 if( init_up_down_time > 0 ) {
  // the unit is already on- - - - - - - - - - - - - - - - - - - - - - - - -

  // s therefore works as an ON-node: construct the EDSolver
  f_start.DPS = new DPEDSolver( 0 , this );

  // compute kMin, the first time step the unit can be turned OFF due to
  // the need to reach power bound_down[ i ] from the initial power
  // initial_power respecting the ramp-down constraints

  Index kMin = 0;
  for( auto tmp = initial_power ;
       ( kMin < time_horizon ) && ( tmp >= bound_down[ kMin ] + eps ) ; )
   tmp -= delta_ramp_down[ kMin++ ];

  if( kMin < init_t )  // the ramp-down time is less than the time required
   kMin = init_t;      // by the min up-time constraints: use the latter

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
  // the unit is on but it is immediately turned off: this "oddball"
  // arc corresponds to an empty ED and always has 0 cost

  double fc = 0;             // compute the fixed-cost component of the cost
  Index j = 0;               // this surely comprises the fixed costs 
  while( j < kMin )          // between 0 (included) and kMin (excluded)
   fc += const_term[ j++ ];  // since the unit is on in that period

  f_start.v_arcs.resize( time_horizon - kMin + 1 );
  auto ai = f_start.v_arcs.begin();

  // construct the "normal" arcs up to ( s , time_horizon - 1 )
  for( ; j < time_horizon ; ++j , ++ai ) {
   ai->cost1 = fc;
   ai->cost2 = 0;
   ai->tail = & v_off_nodes[ j ];
   ai->tail->lab = 1;      // mark the tail node as reachable
   fc += const_term[ j ];  // the next fixed cost will comprise that at j
   }

  // now construct the special last arc ( s , d ); note that the fixed
  // cost from 0 to n - 1 (included) has been computed already
  ai->cost1 = fc;
  ai->cost2 = 0;
  ai->tail = & f_end;
  }
 else {  // init_up_down_time <= 0, the unit was off- - - - - - - - - - - -

  // s therefore works as an OFF-node: f_start.DPS must be nullptr
  f_start.DPS = nullptr;

  // allocate the set of arcs: these are
  //
  //       time_horizon - init_t + 1
  //
  // (note that init_t <= time_horizon, so at least one arc is there)
  // where note that init_t == 0 is now possible meaning that
  // init_up_down_time == min_down_time == 0; this implies that the first
  // arc is ( s , 0 ), i.e., "the unit was off at the beginning but it
  // starts up immediately". Apart from this the structure of the arcs is
  // analogous as in the init_up_down_time > 0 case, except of course they
  // go to the ON nodes

  f_start.v_arcs.resize( time_horizon - init_t + 1 );
  auto ai = f_start.v_arcs.begin();

  // construct the "normal" arcs up to ( i , time_horizon - 1 )
  for( Index j = init_t ; j < time_horizon ; ++j , ++ai ) {
   ai->cost1 = compute_startup_costs( 0 , j );
   ai->cost2 = 0;
   ai->tail = & v_on_nodes[ j ];
   ai->tail->lab = 1;            // mark the tail node as reachable
   }

  // now construct the special last arc ( s , d ); note that the fixed
  // cost is 0 because no startup ever happens during the time horizon
  ai->cost1 = 0;
  ai->cost2 = 0;
  ai->tail = & f_end;

  }  // end( else( the unit was off ) )

 // now build the ON and OFF nodes - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // we do this in the order i = 0, 1, ..., n - 1: since the graph is
 // acyclic, if the lab of the node is still 0 when we process it then the
 // node is unreachable from d and we need not construct any arc

 const Index mut = std::max( min_up_time , Index( 1 ) );
 // min up-time of 0 makes no sense
 const Index mdt = std::max( min_down_time , Index( 1 ) );
 // min down-time of 0 does make sense, but OFF arcs always go forward by
 // at least one time instant, so we pretend that 1 is the minimum value

 for( Index i = 0 ; i < time_horizon ; ++i ) {
  // process ON node ( i , 1 ) - - - - - - - - - - - - - - - - - - - - - - -
  if( v_on_nodes[ i ].lab ) {  // ... but only if it is reachable

   // allocate and initialise the EDSolver of the node
   v_on_nodes[ i ].DPS = new DPEDSolver( i , this );

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

   v_on_nodes[ i ].v_arcs.resize( time_horizon - endi + 1 );
   auto ai = v_on_nodes[ i ].v_arcs.begin();

   // construct the "normal" arcs up to ( i , time_horizon - 1 )
   for( ; j < time_horizon ; ++j , ++ai ) {
    ai->cost1 = fc;
    ai->cost2 = 0;
    ai->tail = & v_off_nodes[ j ];
    ai->tail->lab = 1;      // mark the tail node as reachable
    fc += const_term[ j ];  // the next fixed cost will comprise that at j
    }

   // now construct the special last arc ( i , d ); note that the fixed
   // cost from i to n - 1 (included) has been computed already
   ai->cost1 = fc;
   ai->cost2 = 0;
   ai->tail = & f_end;

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
   //   down-time periods to be at least of lenght one. Thus, the structure
   //   of the arcs is analogous as in the ON nodes, except of course they
   //   go to the ON nodes themselves
   //
   // - if, instead, i + min_down_time >= time_horizon, then there only is
   //   the single arc ( i , d ) corresponding to "the unit remains off
   //   from i to the end of the horizon, and it will have to remain off
   //   after (but this is not our concern)

   Index endi = std::min( time_horizon , i + mdt );
   v_off_nodes[ i ].v_arcs.resize( time_horizon - endi + 1 );
   auto ai = v_off_nodes[ i ].v_arcs.begin();

   // construct the "normal" arcs up to ( i , time_horizon - 1 )
   for( Index j = i + mdt ; j < time_horizon ; ++j , ++ai ) {
    ai->cost1 = compute_startup_costs( i , j );
    ai->cost2 = 0;
    ai->tail = & v_on_nodes[ j ];
    ai->tail->lab = 1;            // mark the tail node as reachable
    }

   // now construct the special last arc ( i , d ); note that the fixed
   // cost is 0 because no startup ever happens during the time horizon
   ai->cost1 = 0;
   ai->cost2 = 0;
   ai->tail = & f_end;

   }  // end( if( reached )
  }  // end( for( i ) )

 // the graph is now constructed- - - - - - - - - - - - - - - - - - - - - - -
 // nothing needs be done for the destination d

 stage = graph_OK;  // update stage

 }  // end( build_graph )

/*--------------------------------------------------------------------------*/

void ThermalUnitDPSolver::compute_EDPs( void )
{
 if( stage < graph_OK )
  throw( std::logic_error( "compute_EDPs(): graph not ready" ) );

 std::vector< double > cost( time_horizon );

 // update variable costs in the arcs outgoing from s
 if( f_start.DPS ) {                    // f_start is a ON node
  f_start.DPS->compute_costs( cost );   // solve EDPs, retrieve optimal costs

  // index of first tail node (note: one arc surely exists)
  Index h = h_of_node( f_start.v_arcs.front().tail );

  // the cost of ( s , h ) is found in cost[ h - 1 ]: however, one must
  // be careful of the weird case where ( s , 0 ) is present, i.e.,
  // the unit is on but it immediately shuts down. this arc has cost 0
  // as "nothing happens there". this is how the arc cost is initialized
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

 // update variable costs in the arcs outgoing from every ON( i )
 for( Index i = 0 ; i < time_horizon ; ++i ) {
  if( v_on_nodes[ i ].v_arcs.empty() )  // unless it is unreachable
   continue;                            // in which case it is skipped

  // solve EDPs, retrieve optimal costs
  v_on_nodes[ i ].DPS->compute_costs( cost );

  // index of first tail node
  Index h = h_of_node( v_on_nodes[ i ].v_arcs.front().tail );

  // the cost of ( i , h ) is found in cost[ h - 1 ]; note that h > i,
  // and therefore h > 0, and therefore h - 1 is well defined
  --h;

  // set the variable costs in the arcs
  for( auto & a : v_on_nodes[ i ].v_arcs )
   a.cost2 = cost[ h++ ];
  }

 stage = edps_OK;  // update stage

 }  // end( compute_EDPs )

/*--------------------------------------------------------------------------*/

void ThermalUnitDPSolver::min_path( void )
{
 if( stage < edps_OK )
  throw( std::logic_error( "min_path: graph and/or EDPs not ready" ) );

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

 stage = path_OK;  // all done: update stage

 }  // end( min_path )

/*--------------------------------------------------------------------------*/

void ThermalUnitDPSolver::compute_solutions( void )
{
 if( stage < edps_OK )
  throw( std::logic_error( "compute_solutions: graph and/or path not ready"
			   ) );

 std::fill( P.begin() , P.end() , 0 );
 std::fill( U.begin() , U.end() , false );

 Index k = time_horizon;
 auto n = f_end.pred;

 if( ! n )
  throw( std::logic_error(
     "compute_solutions: called when has_var_solution() == false" ) );

 // compute the solution by visiting the optimal path backward from f_end

 do {
  Index h = h_of_node( n );   // the current arc is ( h , k )
  if( n->DPS && k ) {
   // n is ON( h ), or the source (if h == 0) that works as an ON node
   // the power and committment variables of this arc are these with index
   // h, ..., k - 1, comprised if n == f_start (this is why h_of_node()
   // returns 0 for it); however, one has to explicitly avoid the special
   // case of the "empty" arc ( s , 0 ) that has no power and committment
   // variables
   // get optimal values of power variables out of the EDSolver
   n->DPS->compute_power_variables( k - 1 , P );
   for( Index i = h ; i < k ; )   // set all commitment variables to true
    U[ i++ ] = true;
   }
  // else n is OFF( h ), or the source (if h == 0) that works as an OFF
  // node: P[ i ] = U[ i ] = 0 for i = h, ..., k - 1, but these already
  // have those values

  k = h;         // the previous beginning will be the end
  n = n->pred;   // back one arc

  } while( n );  // ... until we hit f_start that has pred == nullptr

 stage = sol_OK;  // all done: update stage

 }  // end( compute_solutions )

/*--------------------------------------------------------------------------*/
/*-------------------- PRIVATE METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

void ThermalUnitDPSolver::load_parameters( void )
{
 // locking the Block
 bool owned = f_Block->is_owned_by( f_id );
 if( ( ! owned ) && ( ! f_Block->read_lock() ) )
  throw( std::runtime_error( "Unable to lock the Block" ) );

 // casting has been checked in set_Block() already
 auto b = static_cast< ThermalUnitBlock * >( f_Block );

 // sanity checks
 if( ! b->get_primary_rho().empty() )
  throw( std::invalid_argument(
	     "ThermalUnitDPSolver does not handle primary reserve yet" ) );
 
 if( ! b->get_secondary_rho().empty() )
  throw( std::invalid_argument(
	   "ThermalUnitDPSolver does not handle secondary reserve yet" ) );

 // scalar values
 time_horizon = b->get_time_horizon();
 init_up_down_time = b->get_init_up_down_time();
 min_up_time = b->get_min_up_time();
 min_down_time = b->get_min_down_time();
 initial_power = b->get_initial_power();

 // init_t: first instant in which a decision can be made, as all the
 //         instants before are "blocked" by the initial conditions
 if( init_up_down_time > 0 )
  if( min_up_time > init_up_down_time )
   init_t = std::min( time_horizon , min_up_time - init_up_down_time );
  else
   init_t = 0;
 else
  if( min_down_time > - init_up_down_time )
   init_t = std::min( time_horizon , min_down_time + init_up_down_time );
  else
   init_t = 0;

 // power vectors
 startup_costs = b->get_start_up_cost();
 min_power = b->get_min_power();
 max_power = b->get_max_power();

 if( b->get_delta_ramp_up().empty() )
  delta_ramp_up = max_power;
 else
  delta_ramp_up = b->get_delta_ramp_up();

 if( b->get_delta_ramp_down().empty() )
  delta_ramp_down = max_power;
 else
  delta_ramp_down = b->get_delta_ramp_down();

 retrieve_term( quad_term, b->get_quad_term() );
 retrieve_term( linear_term, b->get_linear_term() );
 retrieve_term( const_term, b->get_const_term() );

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

double ThermalUnitDPSolver::compute_startup_costs( Index h , Index k )
{
 // one day a time-dependent SUC formula may be easily implemented here
 return( startup_costs[ k ] );
 }

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
   // if a reset is done, all the remaining Modifications can be ignored
   v_mod.clear();
   break;
   }

 f_mod_lock.clear( std::memory_order_release );  // release lock

 if( reload )
  load_parameters();

 }  // end( process_modifications )

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
    case ThermalUnitBlockMod::eSetMaxP:
     max_power = b->get_max_power();
     if( stage > graph_OK )
      stage = graph_OK;
     return( false );

    case ThermalUnitBlockMod::eSetInitP:
     initial_power = b->get_initial_power();
     stage = start;
     return( false );

    case ThermalUnitBlockMod::eSetInitUD:
     init_up_down_time = b->get_init_up_down_time();
     min_up_time = b->get_min_up_time();
     min_down_time = b->get_min_down_time();
     if( init_up_down_time > 0 )
      if( min_up_time > init_up_down_time )
       init_t = std::min( time_horizon , min_up_time - init_up_down_time );
      else
       init_t = 0;
     else
      if( min_down_time > - init_up_down_time )
       init_t = std::min( time_horizon , min_down_time + init_up_down_time );
      else
       init_t = 0;
     stage = start;
     return( false );

    case ThermalUnitBlockMod::eSetAv:
     return( true );  // TODO

    case ThermalUnitBlockMod::eSetSUC:
     startup_costs = b->get_start_up_cost();
     if( stage > edps_OK )
      stage = edps_OK;
     return( false );

    case ThermalUnitBlockMod::eSetLinT:
     retrieve_term( linear_term, b->get_linear_term() );
     if( stage > graph_OK )
      stage = graph_OK;
     return( false );

    case ThermalUnitBlockMod::eSetQuadT:
     retrieve_term( quad_term, b->get_quad_term() );
     if( stage > graph_OK )
      stage = graph_OK;
     return( false );

    case ThermalUnitBlockMod::eSetConstT:
     retrieve_term( const_term, b->get_const_term() );
     stage = start;
     return( false );

    }  // end( switch )

  return( true );

  }  // end( ThermalUnitBlockMod )

 return( false );  // any other Modification: I assume it's harmless

 }  // end( ThermalUnitDPSolver::guts_of_process_modifications )

/*--------------------------------------------------------------------------*/

void ThermalUnitDPSolver::retrieve_term( std::vector< double > & out ,
                                         const std::vector< double > & in )
 const
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
/*----------- METHODS OF ThermalUnitDPSolver::DPEDSolver -------------------*/
/*--------------------------------------------------------------------------*/

ThermalUnitDPSolver::DPEDSolver::DPEDSolver( Index h ,
					     ThermalUnitDPSolver * s )
 : EDSolver( h , s )
{
 auto & time_horizon = f_solver->time_horizon;

 Index coeffsize = time_horizon * time_horizon + f_h * f_h -
                   2 * f_h * time_horizon;
 if( coeffsize != coeffs.size() ) {
  coeffs.resize( coeffsize );
  m.resize( coeffsize + time_horizon - f_h );
  v.resize( time_horizon );
  pos.resize( time_horizon );
  unc_p.resize( time_horizon );
  con_p.resize( time_horizon );
  }
 }

/*--------------------------------------------------------------------------*/

void ThermalUnitDPSolver::DPEDSolver::compute_costs(
					      std::vector< double > & costs )
{
 // scalar values
 auto & time_horizon = f_solver->time_horizon;
 auto & init_up_down_time = f_solver->init_up_down_time;
 auto & initial_power = f_solver->initial_power;

 // power vectors
 auto & min_power = f_solver->min_power;
 auto & max_power = f_solver->max_power;
 auto & delta_ramp_up = f_solver->delta_ramp_up;
 auto & delta_ramp_down = f_solver->delta_ramp_down;
 auto & bound_on = f_solver->bound_on;
 auto & bound_down = f_solver->bound_down;

 // coefficients of the objective function
 auto & quad_term = f_solver->quad_term;
 auto & linear_term = f_solver->linear_term;

 Index k = f_h;

 coeffs[ 0 ].alfa = quad_term[ k ];
 coeffs[ 0 ].beta = linear_term[ k ];
 coeffs[ 0 ].gamma = 0;
 Index coeffcnt = 1;  // next free position in coeffs[]
 v[ k ] = 0;          // because for k = h the number of pieces is 1

 /* Initialize the vector m containing the endpoints of the pieces.
  * At first, it contains the two endpoints of the individual piece.
  * At startup, power can't exceed the bound-on value \barl_k.
  * However, if the unit is on at the beginning of the time horizon with
  * the given initial value initial_power, then the interval is restricted
  * to take it into account. */

 if( ( k == 0 ) && ( init_up_down_time > 0 ) ) {
  m[ 0 ] = std::max( min_power[ k ] , initial_power - delta_ramp_down[ k ] );
  m[ 1 ] = std::min( max_power[ k ] , initial_power + delta_ramp_up[ k ] );
  }
 else {
  m[ 0 ] = min_power[ k ];
  m[ 1 ] = std::min( bound_on[ k ] , max_power[ k ] ); // \bar{l}_k;
  }

 Index mcnt = 2; // Next free position in m[]

 // initialize the vector pos containing the initial indices of the pieces.
 pos[ k ].begm = 0;
 pos[ k ].begt = 0;

 // initialize the vector of unconstrained power values, i.e., 
 // power values are not constrained by bound_down[ k ]
 if( std::abs( coeffs[ 0 ].alfa ) <= 1e-16 )
  if( coeffs[ 0 ].beta <= 0 )
   unc_p[ k ] = m[ 1 ];
  else
   unc_p[ k ] = m[ 0 ];
 else {
  // tmp is p^*_{hk}
  double tmp = -coeffs[ 0 ].beta / ( 2 * coeffs[ 0 ].alfa );
  if( tmp < m[ 0 ] )
   unc_p[ k ] = m[ 0 ];
  else
   if( tmp > m[ 1 ] )
    unc_p[ k ] = m[ 1 ];
   else
    unc_p[ k ] = tmp;
  }

 /* Initialize the vector of constrained power values, that will be
  * computed at each iteration.
  * Constrained means that they must be <= bound_down[ k ] */

 if( ( k < time_horizon - 1 ) && ( unc_p[ k ] > bound_down[ k + 1 ] ) )
  con_p[ k ] = bound_down[ k + 1 ];
 else
  con_p[ k ] = unc_p[ k ];

 costs[ k ] = coeffs[ 0 ].alfa * con_p[ k ] * con_p[ k ] +
              coeffs[ 0 ].beta * con_p[ k ];

 // outermost loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 for( ; ++k < time_horizon ; ) {
  /* Building pieces: \bar{m}_0 is the first endpoint of the first piece of
   * the z_{hk}(\bar{p}) objective function. Such endpoint will be saved in
   * the m vector. */

  pos[ k ].begm = mcnt;
  pos[ k ].begt = coeffcnt;

  if( min_power[ k ] > m[ pos[ k - 1 ].begm ] - delta_ramp_down[ k - 1 ] )
   m[ mcnt ] = min_power[ k ];
  else
   m[ mcnt ] = m[ pos[ k - 1 ].begm ] - delta_ramp_down[ k - 1 ];

  double p_bar = m[ mcnt ]; // \bar{m}_0
  int v_bar = 0;            // After the case 3 will contain v[ k ]

  // compute q, the index of the piece where p^*(\bar{p}) belongs.

  double pstar;  // p^*(\bar{p})

  if( p_bar < unc_p[ k - 1 ] ) {
   pstar = p_bar + delta_ramp_down[ k - 1 ];
   if( pstar > unc_p[ k - 1 ] )
    pstar = unc_p[ k - 1 ];
   }
  else {
   pstar = p_bar - delta_ramp_up[ k - 1 ];
   if( pstar < unc_p[ k - 1 ] )
    pstar = unc_p[ k - 1 ];
   }

  int qm = pos[ k - 1 ].begm;
  while( ( pstar >= m[ qm + 1 ] ) && ( qm < pos[ k ].begm - 2 ) )
   ++qm;

  int q = qm - pos[ k - 1 ].begm + pos[ k - 1 ].begt;

  // compute the last endpoint of the piece, \bar{u}
  double u_bar = std::min( max_power[ k ] ,
                           m[ mcnt - 1 ] + delta_ramp_up[ k - 1 ] );
  ++mcnt;

  bool firstTime = true;

  // CASE 1- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  while( unc_p[ k - 1 ] > p_bar + delta_ramp_down[ k - 1 ] + f_solver->eps )
  {
   // set coeffs fields to compute \bar{z}^{\bar{v}}(p)

   coeffs[ coeffcnt ].alfa = quad_term[ k ] + coeffs[ q ].alfa;
   coeffs[ coeffcnt ].beta = linear_term[ k ] +  coeffs[ q ].beta +
                          2 * delta_ramp_down[ k - 1 ] * coeffs[ q ].alfa;
   coeffs[ coeffcnt ].gamma = coeffs[ q ].gamma +
    coeffs[ q ].alfa * delta_ramp_down[ k - 1 ] * delta_ramp_down[ k - 1 ] +
    coeffs[ q ].beta * delta_ramp_down[ k - 1 ];

   /* Compute the maximum value for \bar{p} such that:
    *  - p^*_k(\bar{p}) stays in the q-th interval;
    *  - unc_p stays out of the admissible range;
    *  - \bar{p} stays admissible. */

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
       2 * coeffs[ coeffcnt ].alfa * p_bar + coeffs[ coeffcnt ].beta > 0 ) {
    if( std::abs( coeffs[ coeffcnt ].alfa ) <= 1e-16 ) {
     if( coeffs[ coeffcnt ].beta >= 0 )
      unc_p[ k ] = m[ mcnt - 2 ];
     // else do nothing, the function is still decreasing in the next interval
     }
    else {
     unc_p[ k ] = -coeffs[ coeffcnt ].beta / ( 2 * coeffs[ coeffcnt ].alfa );
     if( unc_p[ k ] < m[ mcnt - 2 ] ) {
      unc_p[ k ] = m[ mcnt - 2 ];
      }
     }
    firstTime = false;
    }

   ++coeffcnt;
   }

  // CASE 2- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if( unc_p[ k - 1 ] >= p_bar - delta_ramp_up[ k - 1 ] ) {

   // set coeffs fields to compute \bar{z}^{\bar{v}}(p)

   coeffs[ coeffcnt ].alfa = quad_term[ k ];
   coeffs[ coeffcnt ].beta = linear_term[ k ];
   coeffs[ coeffcnt ].gamma =
    coeffs[ q ].alfa * unc_p[ k - 1 ] * unc_p[ k - 1 ] +
    coeffs[ q ].beta * unc_p[ k - 1 ] + coeffs[ q ].gamma;

   /* compute the maximum value for \bar{p} such that:
    *  - unc_p stays out of the admissible range;
    *  - \bar{p} stays admissible. */

   p_bar = std::min( u_bar , unc_p[ k - 1 ] + delta_ramp_up[ k - 1 ] );

   ++v_bar;
   m[ mcnt++ ] = p_bar;

   if( firstTime &&
       2 * coeffs[ coeffcnt ].alfa * p_bar + coeffs[ coeffcnt ].beta > 0 ) {
    if( std::abs( coeffs[ coeffcnt ].alfa ) <= 1e-16 ) {
     if( coeffs[ coeffcnt ].beta >= 0 )
      unc_p[ k ] = m[ mcnt - 2 ];
     // else do nothing, the function is still decreasing in the next interval
     }
    else {
     unc_p[ k ] = -coeffs[ coeffcnt ].beta / ( 2 * coeffs[ coeffcnt ].alfa );
     if( unc_p[ k ] < m[ mcnt - 2 ] ) {
      unc_p[ k ] = m[ mcnt - 2 ];
      }
     }
    firstTime = false;
    }

   ++coeffcnt;
   }

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
    *  - p^*_k(\bar{p}) stays in the q-th interval;
    *  - \bar{p} stays admissible. */

   p_bar = std::min( m[ qm + 1 ] + delta_ramp_up[ k - 1 ] , u_bar );

   ++v_bar;
   m[ mcnt++ ] = p_bar;
   ++q;
   ++qm;

   if( firstTime &&
       2 * coeffs[ coeffcnt ].alfa * p_bar + coeffs[ coeffcnt ].beta > 0 ) {
    if( std::abs( coeffs[ coeffcnt ].alfa ) <= 1e-16 ) {
     if( coeffs[ coeffcnt ].beta >= 0 )
      unc_p[ k ] = m[ mcnt - 2 ];
     // else do nothing, the function is still decreasing in the next interval
     }
    else {
     unc_p[ k ] = -coeffs[ coeffcnt ].beta / ( 2 * coeffs[ coeffcnt ].alfa );
     if( unc_p[ k ] < m[ mcnt - 2 ] ) {
      unc_p[ k ] = m[ mcnt - 2 ];
      }
     }
    firstTime = false;
    }

   ++coeffcnt;
   }

  // end of the tree cases - - - - - - - - - - - - - - - - - - - - - - - - -

  v[ k ] = v_bar - 1;

  if( firstTime )  // function is strictly decreasing
   unc_p[ k ] = u_bar;

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

  qm = pos[ k ].begm;
  while( ( con_p[ k ] > m[ qm + 1 ] ) && ( m[ qm + 1 ] != 0 ) )
   ++qm;

  q = qm - pos[ k ].begm + pos[ k ].begt;

  costs[ k ] = coeffs[ q ].alfa * con_p[ k ] * con_p[ k ] +
               coeffs[ q ].beta * con_p[ k ] + coeffs[ q ].gamma;

  }  // end( for( k ) )
 }  // end( compute_costs )

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
 }  // end( compute_power_variables )

/*--------------------------------------------------------------------------*/
/*----------------- End File ThermalUnitDPSolver.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
