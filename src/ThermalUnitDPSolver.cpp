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
   throw( std::runtime_error(
		  "ThermalUnitDPSolver only supports  ThermalUnitBlocks" ) );

  load_parameters();
  }
 }

/*--------------------------------------------------------------------------*/

int ThermalUnitDPSolver::compute( bool changedvars )
{
 process_modifications();

 switch( stage ) {
  case start:    build_graph();
  case graph_OK: compute_EDPs();
  case edps_OK:  min_path();
  case path_OK:  compute_solutions();
  }

 assert( stage == sol_OK );
 return( kOK );
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

 for( int i = 0 ; i < time_horizon ; ++i ) {
  pow_it->set_value( P[ i ] );
  com_it->set_value( U[ i ] );
  pow_it++;
  com_it++;
  }

 /*!! set startup variables -- I'd frankly avoid it
 auto sup_it = b->get_start_up();
 for( int i = 0 ; i < time_horizon - init_t ; ++i ) {
  sup_it->set_value( startup[ i ] );
  sup_it++;
  }
 */

 // unlock the block
 if( ! owned )
  f_Block->unlock( f_id );

 }  // end( ThermalUnitDPSolver::get_var_solution( )

/*--------------------------------------------------------------------------*/
/*------------------ BUILDING AND SOLVING THE DP PROBLEM -------------------*/
/*--------------------------------------------------------------------------*/

void ThermalUnitDPSolver::build_graph( void )
{
 v_nodes.resize( time_horizon * 2 );

 /*
  * We identify three cases:
  *
  * 1) The unit is already ON and needs to stay ON for some time steps due
  *    to ramp and/or min_up_time constraints;
  * 2) The unit is already OFF and needs to stay OFF for some time steps due
  *    to min_down_time constraints;
  * 3  The unit, regardless its initial state, is not subjected to any
  *    constraint and can be freely change its status from the beginning.
  */

 /*
  * Compute kMin, the first time step the unit can be turned OFF
  */

 kMin = 0;

 if( initial_power >= bound_down[ 0 ] + eps ) {
  double tmp = initial_power;
  tmp -= delta_ramp_down[ kMin ];
  kMin++;
  while( tmp >= bound_down[ kMin ] + eps ) {
   tmp -= delta_ramp_down[ kMin ];
   kMin++;
  }
  kMin--;
 }

 if( kMin >= time_horizon ) {
  kMin = time_horizon - 1;
 }

 /*
  * FIRST CASE: The unit is already ON and needs to stay ON for some
  * time steps due to ramp and/or min_up_time constraints.
  */
 if( init_up_down_time > 0 &&
     ( init_up_down_time < min_up_time || kMin > 0 ) ) {

  /*
   * Compute the time steps the unit must stay ON due to min_up_time (k)
   * and ramp constraints (kMin from before)
   */
  int h = hMin = 0;
  int k = min_up_time - init_up_down_time - 1;
  if( kMin < k ) {
   kMin = k;
  } else {
   k = kMin;
  }

  /*
   * We start from the node (0, ON) and we build arc connections:
   * - To all nodes (kMin, OFF), (kMin+1, OFF), ..., (th -1, OFF),
   *   since kMin is the first time step the unit can be turned OFF;
   * - The node (t, ON) that denotes the case in which the unit stays ON
   *   throughout the whole period. FIXME: This doesn't exist
   */

  {
   double c_i = 0;
   for( int t = h; t < k; ++t ) {
    c_i += const_term[ t ]; // TODO: Check if correct
   }

   v_nodes[ ON( 0 ) ].v_arcs.resize( time_horizon - k );
   int i = 0;

   // Arcs to connect node (0, ON) with
   // nodes (k, OFF), (k+1, OFF), ..., (time_horizon - 1, OFF)
   for( ; k < time_horizon; ++k ) {
    v_nodes[ ON( 0 ) ].v_arcs[ i ].h = h;
    v_nodes[ ON( 0 ) ].v_arcs[ i ].k = k;
    c_i += const_term[ k ]; // TODO: Check if correct
    v_nodes[ ON( 0 ) ].v_arcs[ i ].cost1 = c_i;
    v_nodes[ ON( 0 ) ].v_arcs[ i ].cost2 = 0;
    v_nodes[ ON( 0 ) ].v_arcs[ i ].valid = -1; // ON -> OFF
    i++;
   }

   // FIXME: Why the arc to connect node (0, ON) with node (t, ON) wasn't here?
  }

  /*
   * Build the arcs for the OFF nodes
   * (kMin, OFF), (kMin + 1, OFF), ..., (th - 1, OFF),
   * since the unit couldn't be turned OFF before k. For each node,
   * we identify two cases:
   */
  for( h = kMin; h < time_horizon; ++h ) {

   if( h > time_horizon - 1 - min_down_time ) {

    /*
     * 1) For time steps equal/greater than t - min_down_time, the unit cannot
     *    be turned ON without breaking the min_down_time constraint.
     *    So the only connection is to (t, OFF), that is the unit staying OFF.
     *    // FIXME: It's actually connected to (th-1, OFF). Why?
     */
    v_nodes[ OFF( h ) ].v_arcs.resize( 1 );
    v_nodes[ OFF( h ) ].v_arcs[ 0 ].h = h;
    v_nodes[ OFF( h ) ].v_arcs[ 0 ].k = time_horizon - 1;
    v_nodes[ OFF( h ) ].v_arcs[ 0 ].cost1 = 0;
    v_nodes[ OFF( h ) ].v_arcs[ 0 ].cost2 = 0;
    v_nodes[ OFF( h ) ].v_arcs[ 0 ].valid = -2; // OFF -> OFF

   } else {

    /*
     * 2) For time steps smaller than t - min_down_time, each node (j, OFF)
     *    is connected to all the nodes (j + min_down_time, ON), ..., (t, ON)
     *    since they respect the min_down_time constraint.
     *    Also, they are connected to (t, OFF), that is the unit staying OFF.
     *    // FIXME: It's actually connected to (th-1, OFF). Why?
     */
    k = h + min_down_time - 1;
    v_nodes[ OFF( h ) ].v_arcs.resize( time_horizon - k + 1 );

    int i = 0;
    for( ; k < time_horizon; ++k ) {
     v_nodes[ OFF( h ) ].v_arcs[ i ].h = h;
     v_nodes[ OFF( h ) ].v_arcs[ i ].k = k;
     v_nodes[ OFF( h ) ].v_arcs[ i ].cost1 = 0;
     v_nodes[ OFF( h ) ].v_arcs[ i ].cost2 = 0; // Startup cost
     v_nodes[ OFF( h ) ].v_arcs[ i ].valid = 1; // OFF -> ON
     i++;
    }

    v_nodes[ OFF( h ) ].v_arcs[ i ].h = h;
    v_nodes[ OFF( h ) ].v_arcs[ i ].k = time_horizon - 1;
    v_nodes[ OFF( h ) ].v_arcs[ i ].cost1 = 0;
    v_nodes[ OFF( h ) ].v_arcs[ i ].cost2 = 0;
    v_nodes[ OFF( h ) ].v_arcs[ i ].valid = -2; // OFF -> OFF
   }
  }

  /*
   * Build the arcs for the ON nodes
   * (k + mindowntime, ON), (k + mindowntime + 1, ON), ..., (t, ON),
   * since the unit couldn't be turned ON again before k + mindowntime.
   * For each node, we identify two cases:
   */

  for( h = kMin + min_down_time - 1; h < time_horizon; ++h ) {

   if( h > time_horizon - 1 - min_up_time ) {

    /*
     * 1) For time steps equal/greater than t - min_up_time, the unit cannot
     *    be turned OFF without breaking the min_up_time constraint.
     *    So the only connection is to (t, ON), that is the unit staying ON.
     *    // FIXME: It's actually connected to (th-1, ON). Why?
     */
    v_nodes[ ON( h ) ].v_arcs.resize( 1 );
    v_nodes[ ON( h ) ].v_arcs[ 0 ].h = h;
    v_nodes[ ON( h ) ].v_arcs[ 0 ].k = time_horizon - 1;
    v_nodes[ ON( h ) ].v_arcs[ 0 ].cost1 =
     const_term[ h ] * ( time_horizon - h ); // TODO: Check Niccolò
    v_nodes[ ON( h ) ].v_arcs[ 0 ].cost2 = 0;
    v_nodes[ ON( h ) ].v_arcs[ 0 ].valid = 2; // ON -> ON

   } else {

    /*
     * 2) For time steps smaller than t - min_up_time, each node (j, ON)
     *    is connected to all the nodes (j + min_up_time, OFF), ..., (t, OFF)
     *    since they respect the min_up_time constraint.
     *    Also, they are connected to (t, ON), that is the unit staying ON.
     *    // FIXME: It's actually connected to (th-1, ON). Why?
     */
    k = h + min_up_time - 1;
    double c_i = 0;
    for( int t = h; t < k; ++t ) {
     c_i += const_term[ t ]; // TODO: Check Niccolò
    }

    int i = 0;
    v_nodes[ ON( h ) ].v_arcs.resize( time_horizon - k );
    for( ; k < time_horizon; ++k ) {
     v_nodes[ ON( h ) ].v_arcs[ i ].h = h;
     v_nodes[ ON( h ) ].v_arcs[ i ].k = k;
     c_i += const_term[ k ]; // TODO: Check Niccolò
     v_nodes[ ON( h ) ].v_arcs[ i ].cost1 = c_i;
     v_nodes[ ON( h ) ].v_arcs[ i ].cost2 = 0;  // Startup cost
     v_nodes[ ON( h ) ].v_arcs[ i ].valid = -1; // ON -> OFF
     i++;
    }

    // FIXME: Why the arc to (t, ON) is missing?
   }
  }
 }

  /*
   * SECOND CASE: The unit is already ON and needs to stay OFF for some
   * time steps due to ramp and/or min_down_time constraints.
   */
 else if( init_up_down_time < 0 &&
          -init_up_down_time < min_down_time ) {

  /*
   * Compute the time steps the unit must stay OFF due to min_down_time
   */
  int k = min_down_time + init_up_down_time;
  hMin = k;
  kMin = hMin + min_up_time - 1;

  /*
   * We start from the node (0, OFF) and we build arc connections:
   * - To all nodes (hMin, ON), (hMin+1, ON), ..., (th-1, ON),
   *   since hMin is the first time step the unit can be turned ON;
   * - The node (t, OFF) that denotes the case in which the unit stays OFF
   *   throughout the whole period.
   *   // FIXME: It's actually connected to (th-1, OFF). Why?
   */
  {
   v_nodes[ OFF( 0 ) ].v_arcs.resize( time_horizon - k + 1 );

   int i = 0;
   for( ; k < time_horizon; ++k ) {
    v_nodes[ OFF( 0 ) ].v_arcs[ i ].h = 0;
    v_nodes[ OFF( 0 ) ].v_arcs[ i ].k = k;
    v_nodes[ OFF( 0 ) ].v_arcs[ i ].cost1 = 0;
    v_nodes[ OFF( 0 ) ].v_arcs[ i ].cost2 = 0; // Startup cost
    v_nodes[ OFF( 0 ) ].v_arcs[ i ].valid = 1; // OFF -> ON
    i++;
   }

   v_nodes[ OFF( 0 ) ].v_arcs[ i ].h = 0;
   v_nodes[ OFF( 0 ) ].v_arcs[ i ].k = time_horizon - 1;
   v_nodes[ OFF( 0 ) ].v_arcs[ i ].cost1 = 0;
   v_nodes[ OFF( 0 ) ].v_arcs[ i ].cost2 = 0;
   v_nodes[ OFF( 0 ) ].v_arcs[ i ].valid = -2; // OFF -> OFF
  }

  /*
   * Build the arcs for the ON nodes
   * (hMin, ON), (hMin + 1, ON), ..., (th - 1, ON),
   * since the unit couldn't be turned ON before h.
   * For each node, we identify two cases:
   */

  for( int h = hMin; h < time_horizon; ++h ) {

   if( h > time_horizon - 1 - min_up_time ) {

    /*
     * 1) For time steps equal/greater than t - min_up_time, the unit cannot
     *    be turned OFF without breaking the min_up_time constraint.
     *    So the only connection is to (t, ON), that is the unit staying ON.
     *    // FIXME: It's actually connected to (th-1, ON). Why?
     */
    v_nodes[ ON( h ) ].v_arcs.resize( 1 );
    v_nodes[ ON( h ) ].v_arcs[ 0 ].h = h;
    v_nodes[ ON( h ) ].v_arcs[ 0 ].k = time_horizon - 1;
    v_nodes[ ON( h ) ].v_arcs[ 0 ].cost1 =
     const_term[ h ] * ( time_horizon - h ); // TODO: Check Niccolò
    v_nodes[ ON( h ) ].v_arcs[ 0 ].cost2 = 0;
    v_nodes[ ON( h ) ].v_arcs[ 0 ].valid = 2; // ON -> ON

   } else {

    /*
     * 2) For time steps smaller than t - min_up_time, each node (j, ON)
     *    is connected to all the nodes (j + min_up_time, OFF), ..., (t, OFF)
     *    since they respect the min_up_time constraint.
     *    Also, they are connected to (t, ON), that is the unit staying ON.
     *    // FIXME: It's actually connected to (th-1, ON). Why?
     */
    k = h + min_up_time - 1;
    double c_i = 0;
    for( int t = h; t < k; ++t ) {
     c_i += const_term[ t ]; // TODO: Check Niccolò
    }

    int i = 0;
    v_nodes[ ON( h ) ].v_arcs.resize( time_horizon - k );
    for( ; k < time_horizon; ++k ) {
     v_nodes[ ON( h ) ].v_arcs[ i ].h = h;
     v_nodes[ ON( h ) ].v_arcs[ i ].k = k;
     c_i += const_term[ k ]; // TODO: Check Niccolò
     v_nodes[ ON( h ) ].v_arcs[ i ].cost1 = c_i;
     v_nodes[ ON( h ) ].v_arcs[ i ].cost2 = 0;  // Startup cost
     v_nodes[ ON( h ) ].v_arcs[ i ].valid = -1; // ON -> OFF
     i++;
    }

    // FIXME: Why the arc to (t, ON) is missing?
   }
  }

  /*
   * Build the arcs for the OFF nodes
   * (k + minuptime, OFF), (k + minuptime + 1, OFF), ..., (t, OFF),
   * since the unit couldn't be turned OFF again before k + minuptime.
   * For each node, we identify two cases:
   */

  for( int h = hMin + min_up_time - 1; h < time_horizon - 1; ++h ) {

   if( h > time_horizon - 1 - min_down_time ) {

    /*
     * 1) For time steps equal/greater than t - min_down_time, the unit cannot
     *    be turned ON without breaking the min_up_time min_down_time.
     *    So the only connection is to (t, OFF), that is the unit staying OFF.
     *    // FIXME: It's actually connected to (th-1, OFF). Why?
     */
    v_nodes[ OFF( h ) ].v_arcs.resize( 1 );
    v_nodes[ OFF( h ) ].v_arcs[ 0 ].h = h;
    v_nodes[ OFF( h ) ].v_arcs[ 0 ].k = time_horizon - 1;
    v_nodes[ OFF( h ) ].v_arcs[ 0 ].cost1 = 0;
    v_nodes[ OFF( h ) ].v_arcs[ 0 ].cost2 = 0;
    v_nodes[ OFF( h ) ].v_arcs[ 0 ].valid = -2; // OFF -> OFF

   } else {

    /*
     * 2) For time steps smaller than t - min_down_time, each node (j, OFF)
     *    is connected to all the nodes (j + min_down_time, ON), ..., (t, ON)
     *    since they respect the min_down_time constraint.
     *    Also, they are connected to (t, OFF), that is the unit staying OFF.
     *    // FIXME: It's actually connected to (th-1, OFF). Why?
     */
    k = h + min_down_time - 1;
    v_nodes[ OFF( h ) ].v_arcs.resize( time_horizon - k + 1 );

    int i = 0;
    for( ; k < time_horizon; ++k ) {
     v_nodes[ OFF( h ) ].v_arcs[ i ].h = h;
     v_nodes[ OFF( h ) ].v_arcs[ i ].k = k;
     v_nodes[ OFF( h ) ].v_arcs[ i ].cost1 = 0;
     v_nodes[ OFF( h ) ].v_arcs[ i ].cost2 = 0; // Startup cost
     v_nodes[ OFF( h ) ].v_arcs[ i ].valid = 1; // OFF -> ON
     i++;
    }

    v_nodes[ OFF( h ) ].v_arcs[ i ].h = h;
    v_nodes[ OFF( h ) ].v_arcs[ i ].k = time_horizon - 1;
    v_nodes[ OFF( h ) ].v_arcs[ i ].cost1 = 0;
    v_nodes[ OFF( h ) ].v_arcs[ i ].cost2 = 0;
    v_nodes[ OFF( h ) ].v_arcs[ i ].valid = -2; // OFF -> OFF
   }
  }
 }

  /*
   * THIRD CASE: The unit, regardless its initial state, is not subjected to any
   * constraint and can be freely change its status from the beginning.
   * In this case, we build all the arc connections for all the time steps.
   */
 else if( init_up_down_time != 0 ) {

  hMin = 0;

  for( int h = 0; h < time_horizon; ++h ) {

   if( h > 0 && init_up_down_time > 0 && h < min_down_time ) {
    continue;
   }

   /*
    * For the OFF nodes, we have two cases:
    */
   if( h > time_horizon - 1 - min_down_time ) {

    /*
     * 1) For time steps equal/greater than t - min_down_time, the unit cannot
     *    be turned ON without breaking the min_down_time constraint.
     *    So the only connection is to (t, OFF), that is the unit staying OFF.
     *    // FIXME: It's actually connected to (th-1, OFF). Why?
     */
    v_nodes[ OFF( h ) ].v_arcs.resize( 1 );
    v_nodes[ OFF( h ) ].v_arcs[ 0 ].h = h;
    v_nodes[ OFF( h ) ].v_arcs[ 0 ].k = time_horizon - 1;
    v_nodes[ OFF( h ) ].v_arcs[ 0 ].cost1 = 0;
    v_nodes[ OFF( h ) ].v_arcs[ 0 ].cost2 = 0;
    v_nodes[ OFF( h ) ].v_arcs[ 0 ].valid = -2; // OFF -> OFF

   } else {

    /*
     * 2) For time steps smaller than t - min_down_time, each node (j, OFF)
     *    is connected to all the nodes (j + min_down_time, ON), ..., (t, ON)
     *    since they respect the min_down_time constraint.
     *    Also, they are connected to (t, OFF), that is the unit staying OFF.
     *    // FIXME: It's actually connected to (th-1, OFF). Why?
     */
    int k = 0;
    if( h != 0 ) {
     k = h + min_down_time - 1;
    }

    v_nodes[ OFF( h ) ].v_arcs.resize( time_horizon - k + 1 );
    int i = 0;
    for( ; k < time_horizon; k++ ) {
     v_nodes[ OFF( h ) ].v_arcs[ i ].h = h;
     v_nodes[ OFF( h ) ].v_arcs[ i ].k = k;
     v_nodes[ OFF( h ) ].v_arcs[ i ].cost1 = 0;
     v_nodes[ OFF( h ) ].v_arcs[ i ].cost2 = 0; // Startup cost
     v_nodes[ OFF( h ) ].v_arcs[ i ].valid = 1; // OFF -> ON
     i++;
    }

    v_nodes[ OFF( h ) ].v_arcs[ i ].h = h;
    v_nodes[ OFF( h ) ].v_arcs[ i ].k = time_horizon - 1;
    v_nodes[ OFF( h ) ].v_arcs[ i ].cost1 = 0;
    v_nodes[ OFF( h ) ].v_arcs[ i ].cost2 = 0;
    v_nodes[ OFF( h ) ].v_arcs[ i ].valid = -2; // OFF -> OFF
   }

   /*
    * For the ON nodes, we have two cases:
    */
   if( h > time_horizon - 1 - min_up_time ) {

    /*
     * 1) For time steps equal/greater than t - min_up_time, the unit cannot
     *    be turned OFF without breaking the min_up_time constraint.
     *    So the only connection is to (t, ON), that is the unit staying ON.
     *    // FIXME: It's actually connected to (th-1, ON). Why?
     */
    v_nodes[ ON( h ) ].v_arcs.resize( 1 );
    v_nodes[ ON( h ) ].v_arcs[ 0 ].h = h;
    v_nodes[ ON( h ) ].v_arcs[ 0 ].k = time_horizon - 1;
    v_nodes[ ON( h ) ].v_arcs[ 0 ].cost1 =
     const_term[ h ] * ( time_horizon - h ); // TODO: Check Niccolò
    v_nodes[ ON( h ) ].v_arcs[ 0 ].cost2 = 0;
    v_nodes[ ON( h ) ].v_arcs[ 0 ].valid = 2; // ON -> ON

   } else {

    /*
     * 2) For time steps smaller than t - min_up_time, each node (j, ON)
     *    is connected to all the nodes (j + min_up_time, OFF), ..., (t, OFF)
     *    since they respect the min_up_time constraint.
     *    Also, they are connected to (t, ON), that is the unit staying ON.
     *    // FIXME: It's actually connected to (th-1, ON). Why?
     */
    int k = 0;
    if( h != 0 ) {
     k = h + min_up_time - 1;
    }
    double c_i = 0;
    for( int t = h; t < k; ++t ) {
     c_i += const_term[ t ]; // TODO: Check Niccolò
    }

    int i = 0;
    v_nodes[ ON( h ) ].v_arcs.resize( time_horizon - k );
    for( ; k < time_horizon; k++ ) {
     v_nodes[ ON( h ) ].v_arcs[ i ].h = h;
     v_nodes[ ON( h ) ].v_arcs[ i ].k = k;
     c_i += const_term[ k ]; // TODO: Check Niccolò
     v_nodes[ ON( h ) ].v_arcs[ i ].cost1 = c_i;
     v_nodes[ ON( h ) ].v_arcs[ i ].cost2 = 0;  // Startup cost
     v_nodes[ ON( h ) ].v_arcs[ i ].valid = -1; // ON -> OFF
     i++;
    }
   }
  }
 } else {
  assert( 0 );
 }

 // Update stage
 stage = graph_OK;
}

/*--------------------------------------------------------------------------*/

void ThermalUnitDPSolver::compute_EDPs( void )
{
 if( stage < graph_OK )
  throw( std::logic_error( "compute_EDPs(): graph not ready" ) );

 std::vector< double > v_cost( time_horizon );
 v_EDP.resize( time_horizon - hMin );

 for( int i = hMin; i < time_horizon; ++i ) {
  if( v_nodes[ i ].v_arcs.empty() )
   continue;

  v_EDP[ i - hMin ].initialize( i, this );
  v_EDP[ i - hMin ].compute_costs( v_cost );

  for( auto & v_arc : v_nodes[ i ].v_arcs ) {
   if( v_arc.valid == -1 || v_arc.valid == 2 ) {
    v_arc.cost2 = v_cost[ v_arc.k ];
   }
  }
 }

 // Update stage
 stage = edps_OK;
}

/*--------------------------------------------------------------------------*/

void ThermalUnitDPSolver::min_path() {

 if( stage < edps_OK ) {
  throw std::logic_error( "min_path(): graph and/or EDPs not ready" );
 }

 v_route.resize( time_horizon + 1 );

 /*
  * FIRST CASE: Unit is initially ON, start up costs are zero.
  */
 if( init_up_down_time > 0 ) {
  if( init_up_down_time < min_up_time ) {

   /*
    * Unit is already ON for less than min_up_time.
    */
   const int h = 0;

   // Initialize all nodes that can be connected with the source
   for( int k = kMin; k < time_horizon; ++k ) {
    v_route[ k ].h = h;
    v_route[ k ].pred = -1;
    v_route[ k ].lab = v_nodes[ h ].v_arcs[ k - kMin ].cost1 +
                       v_nodes[ h ].v_arcs[ k - kMin ].cost2;
   }

   // The target node can not be directly connected to the source node
   v_route[ time_horizon ].lab = Inf<double>();
  } else {

   /*
    * Unit is already ON for more than min_up_time.
    * We take into account the ramp constraint as well.
    */

   if( initial_power < bound_down[ 0 ] + eps ) {

    /*
     * The ramp constraint is not violated,
     * we create the connection between source and target node.
     */

    v_route[ time_horizon ].h = -1;
    v_route[ time_horizon ].lab = 0;
    v_route[ time_horizon ].pred = -1;

    // Initialize the other connections to inf
    for( int k = 0; k < time_horizon; ++k ) {
     v_route[ k ].lab = Inf<double>();
    }

    /*
     * Initialize the possible connections from the case where
     * the unit is being turned off at the beginning.
     */

    int h = min_down_time;
    for( ; h < time_horizon - min_up_time + 1; ++h ) {

     /*
      * Check all the possible starts that allow the unit
      * to be turned off again in the period.
      */
     const double currentstartupcost = compute_startup_costs( h );

     int i = 0;
     for( int k = h + min_up_time - 1; k < time_horizon; ++k, ++i ) {
      const double label = currentstartupcost +
                           v_nodes[ h ].v_arcs[ i ].cost1 +
                           v_nodes[ h ].v_arcs[ i ].cost2;

      if( label < v_route[ k ].lab ) {
       v_route[ k ].h = h;
       v_route[ k ].lab = label;
       v_route[ k ].pred = -1;
      }
     }
    }

    for( int i = 0; h < time_horizon; ++h, ++i ) {

     /*
      * Check the possible starts for which the unit has to remain
      * on until the end of the period.
      */
     const int k = time_horizon - 1;

     const double label = compute_startup_costs( h ) +
                          v_nodes[ h ].v_arcs[ i ].cost1 +
                          v_nodes[ h ].v_arcs[ i ].cost2;

     if( label < v_route[ k ].lab ) {
      v_route[ k ].h = h;
      v_route[ k ].lab = label;
      v_route[ k ].pred = -1;
     }
    }

    /*
     * Check the case where the unit remains on at the beginning,
     * creating pairs of (0, kMin), (0, kMin + 1), ..., (0, n-1).
     */

    for( int k = kMin; k < time_horizon; ++k ) {
     const double label = v_nodes[ 0 ].v_arcs[ k - kMin ].cost1 +
                          v_nodes[ 0 ].v_arcs[ k - kMin ].cost2;
     if( label < v_route[ k ].lab ) {
      v_route[ k ].h = 0;
      v_route[ k ].lab = label;
      v_route[ k ].pred = -1;
     }
    }
   } else {
    // The unit cannot be turned off immediately due to ramp constraints

    // No direct connection between source and target nodes
    v_route[ time_horizon ].h = -1;
    v_route[ time_horizon ].lab = Inf<double>();
    v_route[ time_horizon ].pred = -1;

    /*
     * Initialize all the possible pairs, starting with the first time
     * the unit can be turned off: (0, kMin), (0, kMin + 1), ..., (0, n-1).
     */

    const int h = 0;
    for( int k = kMin; k < time_horizon; ++k ) {
     v_route[ k ].h = h;
     v_route[ k ].pred = -1;
     v_route[ k ].lab = v_nodes[ 0 ].v_arcs[ k - kMin ].cost1 +
                        v_nodes[ 0 ].v_arcs[ k - kMin ].cost2;
    }
   }
  }

 } else {
  /*
   * SECOND CASE: Unit is initially OFF.
   */

  int idxcs = -init_up_down_time < min_down_time ?
              min_down_time :
              -init_up_down_time;
  //calculating the time-steps that the unit was off
  int h = hMin;

  /*
   * Initialize all the pairs that start from the first time the
   * unit can be turned on and are connected with source node.
   */
  {
   const double currentstartupcost = compute_startup_costs( idxcs );

   for( int k = kMin; k < time_horizon; ++k ) {
    const double label = currentstartupcost +
                         v_nodes[ h ].v_arcs[ k - kMin ].cost1 +
                         v_nodes[ h ].v_arcs[ k - kMin ].cost2;
    v_route[ k ].h = h;
    v_route[ k ].lab = label;
    v_route[ k ].pred = -1;
   }
  }

  ++h;
  ++idxcs;

  /*
   * Initialize all the other feasible pairs.
   */

  for( ; h < time_horizon - min_up_time + 1; ++h, ++idxcs ) {
   const double currentstartupcost = compute_startup_costs( idxcs );

   int i = 0;
   for( int k = h + min_up_time - 1; k < time_horizon - 1; ++k ) {
    const double label = currentstartupcost +
                         v_nodes[ h ].v_arcs[ i ].cost1 +
                         v_nodes[ h ].v_arcs[ i ].cost2;
    i++;
    if( label < v_route[ k ].lab ) {
     v_route[ k ].h = h;
     v_route[ k ].lab = label;
     v_route[ k ].pred = -1;
    }
   }

   // Nodes that terminate after the end of the interval of definition
   {
    const int k = time_horizon - 1;
    const double label = compute_startup_costs( idxcs ) +
                         v_nodes[ h ].v_arcs[ i ].cost1 +
                         v_nodes[ h ].v_arcs[ i ].cost2;

    if( label < v_route[ k ].lab ) {
     v_route[ k ].h = h;
     v_route[ k ].lab = label;
     v_route[ k ].pred = -1;
    }
   }
  }

  // Initialize the pair (s, d)
  v_route[ time_horizon ].lab = 0;
  v_route[ time_horizon ].pred = -1;
 }


 /*
  * All other arcs
  *
  * 1) for (k = Kmin; k <time_horizon - get_min_down_time() - min_up_time; k++)
  *
  * 2) for (; k <time_horizon - get_min_down_time() - 1; k++)
  *
  * 3) for (; k <time_horizon; k++)
  *
  * In case 2) the first of two cycles of r currently present is always skipped.
  * The division of the cycle eliminates so min_up_time-1 test operations.
  * In case 3) the first two cycles could be removed (thus leaving only
  * If the bow ((h, k), d)) thus saving '2 * (get_min_down_time() + 1) test operations.
  */

 for( int k = kMin; k < time_horizon; ++k ) {

  int r = k + min_down_time + 1;
  for( ; r < time_horizon - min_up_time + 1; ++r ) {
   const double costbeforenode = v_route[ k ].lab +
                                 compute_startup_costs( r - 1 - k );

   int i = 0;
   for( int q = r + min_up_time - 1; q < time_horizon; ++q ) {
    const double label = costbeforenode +
                         v_nodes[ r ].v_arcs[ i ].cost1 +
                         v_nodes[ r ].v_arcs[ i ].cost2;
    i++;

    if( v_route[ q ].lab > label ) {
     v_route[ q ].h = r;
     v_route[ q ].lab = label;
     v_route[ q ].pred = k;
    }
   }
  }


  for( ; r < time_horizon; ++r ) {
   const int q = time_horizon - 1;
   const double label = v_route[ k ].lab +
                        compute_startup_costs( r - 1 - k ) +
                        v_nodes[ r ].v_arcs[ 0 ].cost1 +
                        v_nodes[ r ].v_arcs[ 0 ].cost2;
   if( v_route[ q ].lab > label ) {
    v_route[ q ].h = r;
    v_route[ q ].lab = label;
    v_route[ q ].pred = k;
   }
  }

  // Arc ((h,k), d)
  {
   const double label = v_route[ k ].lab;
   if( label < v_route[ time_horizon ].lab ) {
    v_route[ time_horizon ].lab = label;
    v_route[ time_horizon ].pred = k;
   }
  }
 }

 // Update stage
 stage = path_OK;
}

/*--------------------------------------------------------------------------*/

void ThermalUnitDPSolver::compute_solutions( void )
{
 if( stage < edps_OK )
  throw( std::logic_error( "compute_solutions(): graph and/or path not ready"
			   ) );

 std::fill( P.begin() , P.end() , 0 );
 std::fill( U.begin() , U.end() , 0 );
 std::fill( startup.begin() , startup.end() , 0 );

 int k = v_route[ time_horizon ].pred;
 while( k != -1 ) {
  int h = v_route[ k ].h;

  // Compute active power values
  v_EDP[ h - hMin ].compute_power_variables( k, P );

  // Fill startup variable values
  if( h != 0 || init_up_down_time <= 0 ) {
   startup[ h ] = 1;
  }

  // Fill commitment variable values
  for( int t = h; t <= k; ++t ) {
   U[ t ] = 1;
  }

  k = v_route[ k ].pred;
 }

 // Get total cost
 total_cost = v_route[ time_horizon ].lab;

 // Update stage
 stage = sol_OK;
}

/*--------------------------------------------------------------------------*/
/*-------------------- PRIVATE FIELDS OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

void ThermalUnitDPSolver::load_parameters( void )
{
 // Locking the Block
 bool owned = f_Block->is_owned_by( f_id );
 if( ( ! owned ) && ( ! f_Block->read_lock() ) )
  throw( std::runtime_error( "Unable to lock the Block" ) );

 // Casting should have be checked in set_Block() already
 auto b = static_cast< ThermalUnitBlock * >( f_Block );

 // Scalar values
 time_horizon = ( int ) b->get_time_horizon();
 init_up_down_time = b->get_init_up_down_time();
 min_up_time = ( int ) b->get_min_up_time();
 min_down_time = ( int ) b->get_min_down_time();
 initial_power = b->get_initial_power();

 // Init_t (useful for startup variables)
 if( init_up_down_time > 0 ) {
  init_t = init_up_down_time >= min_up_time ?
           0 : min_up_time - init_up_down_time;
 } else {
  init_t = -init_up_down_time >= min_down_time ?
           0 : min_down_time + init_up_down_time;
 }

 // Power vectors
 startup_costs = b->get_start_up_cost();
 min_power = b->get_min_power();
 max_power = b->get_max_power();

 if( b->get_delta_ramp_up().empty() ) {
  delta_ramp_up = max_power;
 } else {
  delta_ramp_up = b->get_delta_ramp_up();
 }

 if( b->get_delta_ramp_down().empty() ) {
  delta_ramp_down = max_power;
 } else {
  delta_ramp_down = b->get_delta_ramp_down();
 }

 retrieve_term( quad_term, b->get_quad_term() );
 retrieve_term( linear_term, b->get_linear_term() );
 retrieve_term( const_term, b->get_const_term() );

 // Unlock the Block
 if( !owned ) {
  f_Block->read_unlock();
 }

 P.resize( time_horizon );
 U.resize( time_horizon );
 startup.resize( time_horizon );
 stage = start;
}

/*--------------------------------------------------------------------------*/

double ThermalUnitDPSolver::compute_startup_costs( int t )
{
 return t > min_down_time ?
        startup_costs[ min_down_time ] :
        startup_costs[ t - min_down_time ];

 /** se l'unita' e' spenta da piu' di fMaxStartupLevel istanti allora
 i costi di start-up sono costanti */

 // if( t > min_down_time ) {
 //  t = min_down_time;
 // }

 /** altrimenti restituisci il costo corrispondente al numero di ore
 di spegnimento, perche' il primo valore di t nella tabella
 coincide con get_min_down_time, che e' proprio il minimo numero di
 istanti per cui l'unita' puo' essere spenta */

 // t -= min_down_time;
 // return ( startup_costs[ t ] );
 }

/*--------------------------------------------------------------------------*/

void ThermalUnitDPSolver::process_modifications( void )
{
 bool reload = false;

 // A function like this is needed to be called
 // recursively with GroupModifications
 // -------------------------------------------
 std::function< void( sp_Mod ) > f;

 f = [ this, &f, &reload ]( const sp_Mod & mod ) {

  // Group modification
  if( const auto gm = std::dynamic_pointer_cast< GroupModification >( mod ) ) {
   for( const auto & submod : gm->sub_Modifications() ) {
    f( submod );
   }
   return;
  }

  // ThermalUnitBlockMod
  if( const auto tubm = std::dynamic_pointer_cast< ThermalUnitBlockMod >( mod ) ) {
   auto b = static_cast< ThermalUnitBlock * >(f_Block);

   switch( tubm->type() ) {
    case ThermalUnitBlockMod::eSetMaxP:
     max_power = b->get_max_power();
     if( stage > graph_OK ) {
      stage = graph_OK;
     }
     break;

    case ThermalUnitBlockMod::eSetInitP:
     initial_power = b->get_initial_power();
     stage = start;
     break;

    case ThermalUnitBlockMod::eSetInitUD:
     init_up_down_time = b->get_init_up_down_time();
     min_up_time = ( int ) b->get_min_up_time();
     min_down_time = ( int ) b->get_min_down_time();
     if( init_up_down_time > 0 ) {
      init_t = init_up_down_time >= min_up_time ?
               0 : min_up_time - init_up_down_time;
     } else {
      init_t = -init_up_down_time >= min_down_time ?
               0 : min_down_time + init_up_down_time;
     }
     stage = start;
     break;

    case ThermalUnitBlockMod::eSetAv:
     // TODO
     reload = true;
     break;

    case ThermalUnitBlockMod::eSetSUC:
     startup_costs = b->get_start_up_cost();
     if( stage > edps_OK ) {
      stage = edps_OK;
     }
     break;

    case ThermalUnitBlockMod::eSetLinT:
     retrieve_term( linear_term, b->get_linear_term() );
     if( stage > graph_OK ) {
      stage = graph_OK;
     }
     break;

    case ThermalUnitBlockMod::eSetQuadT:
     retrieve_term( quad_term, b->get_quad_term() );
     if( stage > graph_OK ) {
      stage = graph_OK;
     }
     break;

    case ThermalUnitBlockMod::eSetConstT:
     retrieve_term( const_term, b->get_const_term() );
     stage = start;
     break;

    default:
     reload = true;
   }

   return;
  }

  // ThermalUnitBlockMod
  if( std::dynamic_pointer_cast< NBModification >( mod ) ) {
   reload = true;
   return;
  }
 };
 // -------------------------------------------

 // Process the Modifications
 for( auto mod = front(); mod; mod = front() ) {
  f( mod );
  pop_front();

  if( std::dynamic_pointer_cast< NBModification >( mod ) ) {
   // An NBModification has just been handled.
   // All the remaining Modifications must be ignored.
   while( front() )
    pop_front();
   break;
  }
 }

 if( reload )
  load_parameters();
 }

/*--------------------------------------------------------------------------*/

void
ThermalUnitDPSolver::retrieve_term( std::vector< double > & out,
                                         const std::vector< double > & in ) const {
 if( in.empty() ) {
  out.resize( time_horizon );
  std::fill( out.begin(), out.end(), 0 );
  return;
 }

 if( in.size() == 1 ) {
  out.resize( time_horizon );
  std::fill( out.begin(), out.end(), in[ 0 ] );
  return;
 }

 assert( in.size() == time_horizon );
 out = in;
}

/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

void EDPSolver::initialize( int k , ThermalUnitDPSolver * s )
{
 solver = s;
 auto & time_horizon = solver->time_horizon;

 h = k;
 kMax = time_horizon;

 int coeffsize = time_horizon * time_horizon + h * h - 2 * h * time_horizon;
 if( coeffsize != coeffs.size() ) {
  coeffs.resize( coeffsize );

  int msize = coeffsize + time_horizon - h;
  m.resize( msize );

  v.resize( time_horizon );
  pos.resize( time_horizon );
  unc_p.resize( time_horizon );
  con_p.resize( time_horizon );
 }
}

/*--------------------------------------------------------------------------*/

void EDPSolver::compute_costs( std::vector< double > & costs )
{
 // Scalar values
 auto & time_horizon = solver->time_horizon;
 auto & init_up_down_time = solver->init_up_down_time;
 auto & initial_power = solver->initial_power;

 // Power vectors
 auto & min_power = solver->min_power;
 auto & max_power = solver->max_power;
 auto & delta_ramp_up = solver->delta_ramp_up;
 auto & delta_ramp_down = solver->delta_ramp_down;
 auto & bound_on = solver->bound_on;
 auto & bound_down = solver->bound_down;

 // Coefficients of the objective function
 auto & quad_term = solver->quad_term;
 auto & linear_term = solver->linear_term;

 int k = h;

 coeffs[ 0 ].alfa = quad_term[ k ];
 coeffs[ 0 ].beta = linear_term[ k ];
 coeffs[ 0 ].gamma = 0;
 int coeffcnt = 1; // Next free position in coeffs[]
 v[ k ] = 0;       // Because for k = h the number of pieces is 1

 /* Initialize the vector m containing the endpoints of the pieces.
  * At first, it contains the two endpoints of the individual piece.
  * At startup, power can't exceed the bound-on value \barl_k.
  * However, if the unit is on at the beginning of the time horizon with
  * the given initial value initial_power, then the interval is restricted
  * to take it into account. */

 if( ( h == 0 ) && ( init_up_down_time > 0 ) ) {
  m[ 0 ] = std::max( min_power[ k ] , initial_power - delta_ramp_down[ k ] );
  m[ 1 ] = std::min( max_power[ k ] , initial_power + delta_ramp_up[ k ] );
  }
 else {
  m[ 0 ] = min_power[ k ];
  m[ 1 ] = std::min( bound_on[ k ] , max_power[ k ] ); // \bar{l}_k;
  }

 int mcnt = 2; // Next free position in m[]

 /*
  * Initialize the vector pos containing the initial indices of the pieces.
  */

 pos[ k ].begm = 0;
 pos[ k ].begt = 0;

 /*
  * Initialize the vector of unconstrained power values.
  * Unconstrained means that power values are not constrained by bound_down[k].
  */
 if( std::abs( coeffs[ 0 ].alfa ) <= 1e-16 )
  if( coeffs[ 0 ].beta <= 0 )
   unc_p[ k ] = m[ 1 ];
  else
   unc_p[ k ] = m[ 0 ];
 else {
  // tmp is p^{*}_{hk}
  double tmp = -coeffs[ 0 ].beta / ( 2 * coeffs[ 0 ].alfa );
  if( tmp < m[ 0 ] ) {
   unc_p[ k ] = m[ 0 ];
  } else if( tmp > m[ 1 ] ) {
   unc_p[ k ] = m[ 1 ];
  } else {
   unc_p[ k ] = tmp;
  }
 }

 /* Initialize the vector of constrained power values, that will be
  * computed at each iteration.
  * Constrained means that they must be <= bound_down[ k ]. */

 if( ( k < time_horizon - 1 ) && ( unc_p[ k ] > bound_down[ k + 1 ] ) )
  con_p[ k ] = bound_down[ k + 1 ];
 else
  con_p[ k ] = unc_p[ k ];

 costs[ k ] = coeffs[ 0 ].alfa * con_p[ k ] * con_p[ k ] +
              coeffs[ 0 ].beta * con_p[ k ];

 // Outermost loop
 for( k = h + 1; k < kMax; ++k ) {

  /*
   * Building pieces: \bar{m}_0 is the first endpoint of the first piece of
   * the z_{hk}(\bar{p}) objective function. Such endpoint will be saved in
   * the m vector.
   */

  pos[ k ].begm = mcnt;
  pos[ k ].begt = coeffcnt;


  if( min_power[ k ] > m[ pos[ k - 1 ].begm ] - delta_ramp_down[ k - 1 ] ) {
   m[ mcnt ] = min_power[ k ];
  } else {
   m[ mcnt ] = m[ pos[ k - 1 ].begm ] - delta_ramp_down[ k - 1 ];
  }

  double p_bar = m[ mcnt ]; // \bar{m}_0
  int v_bar = 0;            // After the case 3 will contain v[k]

  /*
   * Compute q, the index of the piece where p^*(\bar{p}) belongs.
   */

  double pstar; // p^*(\bar{p})

  if( p_bar < unc_p[ k - 1 ] ) {
   pstar = p_bar + delta_ramp_down[ k - 1 ];
   if( pstar > unc_p[ k - 1 ] ) {
    pstar = unc_p[ k - 1 ];
   }
  } else {
   pstar = p_bar - delta_ramp_up[ k - 1 ];
   if( pstar < unc_p[ k - 1 ] ) {
    pstar = unc_p[ k - 1 ];
   }
  }

  int qm = pos[ k - 1 ].begm;
  while( pstar >= m[ qm + 1 ] && qm < pos[ k ].begm - 2 ) {
   ++qm;
  }

  int q = qm - pos[ k - 1 ].begm + pos[ k - 1 ].begt;

  /*
   * Compute the last endpoint of the piece, \bar{u}.
   */

  double u_bar = std::min( max_power[ k ],
                           m[ mcnt - 1 ] + delta_ramp_up[ k - 1 ] );
  // if( max_power[ k ] < m[ mcnt - 1 ] + delta_ramp_up[ k - 1 ] ) {
  //  u_bar = max_power[ k ];
  // } else {
  //  u_bar = m[ mcnt - 1 ] + delta_ramp_up[ k - 1 ];
  // }
  ++mcnt;


  bool firstTime = true;

  // CASE 1
  while( unc_p[ k - 1 ] > p_bar + delta_ramp_down[ k - 1 ] + eps ) {

   /*
    * Set coeffs fields to compute \bar{z}^{\bar{v}}(p).
    */

   coeffs[ coeffcnt ].alfa = quad_term[ k ] + coeffs[ q ].alfa;
   coeffs[ coeffcnt ].beta =
    linear_term[ k ] +
    coeffs[ q ].beta +
    2 * delta_ramp_down[ k - 1 ] * coeffs[ q ].alfa;
   coeffs[ coeffcnt ].gamma =
    coeffs[ q ].gamma +
    coeffs[ q ].alfa * delta_ramp_down[ k - 1 ] * delta_ramp_down[ k - 1 ] +
    coeffs[ q ].beta * delta_ramp_down[ k - 1 ];

   /*
    * Compute the maximum value for \bar{p} such that:
    *  - p^*_k(\bar{p}) stays in the q-th interval;
    *  - unc_p stays out of the admissible range;
    *  - \bar{p} stays admissible.
    */

   if( m[ qm + 1 ] - delta_ramp_down[ k - 1 ] <
       unc_p[ k - 1 ] - delta_ramp_down[ k - 1 ] - eps ) {
    p_bar = m[ qm + 1 ] - delta_ramp_down[ k - 1 ];
    ++q;
    ++qm;
   } else {
    p_bar = unc_p[ k - 1 ] - delta_ramp_down[ k - 1 ];
   }
   if( p_bar > u_bar ) {
    p_bar = u_bar;
   }
   ++v_bar;
   m[ mcnt++ ] = p_bar;

   /*
    * Compute unc_p, unconstrained optimal value for z_{hk}.
    */

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

  // CASE 2
  if( unc_p[ k - 1 ] >= p_bar - delta_ramp_up[ k - 1 ] ) {

   /*
    * Set coeffs fields to compute \bar{z}^{\bar{v}}(p).
    */

   coeffs[ coeffcnt ].alfa = quad_term[ k ];
   coeffs[ coeffcnt ].beta = linear_term[ k ];
   coeffs[ coeffcnt ].gamma =
    coeffs[ q ].alfa * unc_p[ k - 1 ] * unc_p[ k - 1 ] +
    coeffs[ q ].beta * unc_p[ k - 1 ] +
    coeffs[ q ].gamma;

   /*
    * Compute the maximum value for \bar{p} such that:
    *  - unc_p stays out of the admissible range;
    *  - \bar{p} stays admissible.
    */

   if( ( unc_p[ k - 1 ] + delta_ramp_up[ k - 1 ] ) < u_bar ) {
    p_bar = unc_p[ k - 1 ] + delta_ramp_up[ k - 1 ];
   } else {
    p_bar = u_bar;
   }
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


  // CASE 3
  while( p_bar < u_bar ) {

   /*
    * Set coeffs fields to compute \bar{z}^{\bar{v}}(p).
    */

   coeffs[ coeffcnt ].alfa = quad_term[ k ] + coeffs[ q ].alfa;
   coeffs[ coeffcnt ].beta =
    linear_term[ k ] + coeffs[ q ].beta -
    2 * delta_ramp_up[ k - 1 ] * coeffs[ q ].alfa;
   coeffs[ coeffcnt ].gamma =
    coeffs[ q ].gamma +
    coeffs[ q ].alfa * delta_ramp_up[ k - 1 ] * delta_ramp_up[ k - 1 ] -
    coeffs[ q ].beta * delta_ramp_up[ k - 1 ];

   /*
    * Compute the maximum value for \bar{p} such that:
    *  - p^*_k(\bar{p}) stays in the q-th interval;
    *  - \bar{p} stays admissible.
    */

   if( m[ qm + 1 ] + delta_ramp_up[ k - 1 ] < u_bar ) {
    p_bar = m[ qm + 1 ] + delta_ramp_up[ k - 1 ];
   } else {
    p_bar = u_bar;
   }
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

  // End of the tree cases

  v[ k ] = v_bar - 1;


  if( firstTime )  // Function is strictly decreasing
   unc_p[ k ] = u_bar;

  /* Compute con_p[k], constrained optimal value for the entire function.
   */

  if( ( k < time_horizon - 1 ) && ( unc_p[ k ] > bound_down[ k + 1 ] ) )
   con_p[ k ] = bound_down[ k + 1 ];
  else
   con_p[ k ] = unc_p[ k ];
 
  /*
   * Compute the cost for the node (h,k) in costs[].
   */

  qm = pos[ k ].begm;
  while( con_p[ k ] > m[ qm + 1 ] && m[ qm + 1 ] != 0 ) {
   ++qm;
  }
  q = qm - pos[ k ].begm + pos[ k ].begt;

  costs[ k ] =
   coeffs[ q ].alfa * con_p[ k ] * con_p[ k ] +
   coeffs[ q ].beta * con_p[ k ] +
   coeffs[ q ].gamma;

  }  // end( for( k ) )
 }  // end( compute_costs )

/*--------------------------------------------------------------------------*/

void EDPSolver::compute_power_variables( int k , std::vector< double > & p )
{
 auto & delta_ramp_up = solver->delta_ramp_up;
 auto & delta_ramp_down = solver->delta_ramp_down;

 p[ k ] = con_p[ k ];
 for( int t = k - 1 ; t >= h ; --t ) {
  /* Project unconstrained optimal value unc_p[t] on the interval:
   * [ p[t+1] - delta_ramp_up[t], p[t+1] + delta_ramp_down[t] ]
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
 }

/*--------------------------------------------------------------------------*/
/*----------------- End File ThermalUnitDPSolver.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
