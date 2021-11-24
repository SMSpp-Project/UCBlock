/*--------------------------------------------------------------------------*/
/*----------------------- File ThermalUnitDPSolver.h -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the ThermalUnitDPSolver class.
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
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __ThermalUnitDPSolver
 #define __ThermalUnitDPSolver
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

//#include <SMSTypedefs.h>

#include <Solver.h>

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it {

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS ThermalUnitDPSolver ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/

/// Class for solving a Single Unit Commitment problem with a DP approach.
/** The ThermalUnitDPSolver is a Solver for tackling the Single Unit
 * Commitment problem with ramp-up and ramp-down (as well as minimum up-
 * and down-time) constraints and convex quadratic separable objective.
 * The solver uses a Dynamic Programming approach, recasting the problem as
 * a shortest path one on a graph with the following structure
 *
 *       (0,1)  (1,1)  (2,1)  (3,1)  (4,1)  (5,1)  ...  (n-1,1)
 *  (s)                                                           (t)
 *       (0,0)  (1,0)  (2,0)  (3,0)  (4,0)  (5,0)  ...  (n-1,0)
 *
 * n being the length of the time horizon.
 *
 * The fundamental tool for solving this problem is the ability of efficiently
 * solving Economic Dispatch problems that find the min-cost energy production
 * of the unit if it is on for a continuous time interval. In particular we
 * denote by ED( h , k ) the total cost (power + fixed + start-up) of the
 * unit if brought online exactly at the beginning of time h >= 0 and brought
 * offline exactly at the end of time k >= h, i.e., being online for all
 * the time instants h, h + 1, ..., h (note that h = k is possible).
 * Similarly fe denote bu SUC( h , k ) the cost of having the unit offlie
 * from the beginning of h to the end of k; this is typically easy.
 *
 * The problem is therefore reduced to a Shortest Path between (s) and (t)
 * on the acyclic graph constructed as follows:
 *
 * - Each arc from an ON node ( i , 1 ) to an OFF node ( j , 0 ), for
 *   0 <= i < j <= n - 1, means that the unit is started at the beginning of
 *   time i and shut down at the end of time j - 1, so that it is off at
 *   time j. The cost of the arc is therefore ED( i , j - 1 ). However, such
 *   an arc exists only if j - i = number of consecutive periods the unit
 *   remains on is >= min up-time (in particular, i == j, i.e., the unit is
 *   brough online and shut down at the same period, is only possible if
 *   min up-time <= 1, i.e., there is no min up-time requirement).
 *
 * - Each arc from an OFF node ( i , 0 ) to an ON node ( j , 1 ), for
 *   0 <= i < j <= n - 1, means that the unit is shut down at the beginning
 *   of time i and remains down up until the end of time j - 1, then it is
 *   started up at j. Hence, the cost of the arc is the (possibly,
 *   time-variable) start-up costs SUC( i , j - 1 ). However, such an arc
 *   exists only if j - i = number of consecutive periods the unit remains
 *   off is >= min down-time (in particular, i == j, i.e., the unit is
 *   brough offline and online again in the same period, is only possible if
 *   min down-time == 0, i.e., there is no min down-time requirement).

 for the unit after being down
 * for j - i periods. For instance, ( 1 , 0 ) -> ( 5 , 1 ) means that
 * the start-up costs for the unit being off for 4 time steps must be 
 * computed. Note that, symmetrically, in general the arc ( i , 0 ) ->
 * ( j , i ) can only exist if j - i = the number of consecutive periods the
 * unit remains off is >= than the min up-time.

 *
 * - 

computed by solving am
 * an Economic Dispatch Problem for the involved time steps. For instance,
 * ( 1 , 1 ) -> ( 5 , 0 ) is an EDP where the unit has been producing
 * energy for the time-steps 1 - 4. This means that an arc ( i , 1 ) ->
 * ( n , 0 ) means "the unit is producing from i to the end of the horizon";
 * the fact that the unit is actually shut down at n (the initial node of
 * the next horizon) is irrelevant. What is *not* irrelevant, however, is
 * that in general the arc ( i , 1 ) -> ( j , 0 ) can only exist if
 * j - i = number of consecutive periods the unit remains on is >= min
 * up-time. However, *this is not true when i >= n - min up-time*, since
 * of course the unitcan then remain on for the necessary extra time (or
 * more) after the end of the time horizon, but this is of no concern here.
 *
 * Similarly, the arc transition from an OFF node ( i , 0 ) to an ON node
 * ( j , 1 ) means that the unit is shut down at i and remains down up until
 * j - 1, then it is started up at j. Hence, the cost of the arc is the
 * (possibly, time-variable) start-up costs for the unit after being down
 * for j - i periods. For instance, ( 1 , 0 ) -> ( 5 , 1 ) means that
 * the start-up costs for the unit being off for 4 time steps must be 
 * computed. Note that, symmetrically, in general the arc ( i , 0 ) ->
 * ( j , i ) can only exist if j - i = the number of consecutive periods the
 * unit remains off is >= than the min up-time. However, *this is not true
 * when i >= n - min down-time*, because of course the unit will have to
 * remain off for some time after the end of the time horizon, but this is
 * of no concern here. In particular, then *the cost of such an arc must
 * not be the start-up cost but 0*, as the unit will start-up "somewhere in
 * the far future" and that cost is not paid within this time horizon.
 *
 * Hence, the DPSolver needs to compute all the EDPs for the arc
 * transitions that comply with the constraints and the start-up costs,
 * then solve the min-path problem for the graph.
 *
 * ThermalUnitDPSolver first builds the graph, then uses EDPSolver
 * to solve EDPs, then uses a min-path algorithm to solve the commitment
 * problem for the unit.
 */

 class ThermalUnitDPSolver : public Solver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*------------------------------ PUBLIC TYPES ------------------------------*/
/*--------------------------------------------------------------------------*/

  static constexpr TUDPINF = Inf< double >();  ///< the INF value

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and destructor
 *  @{ */

 ThermalUnitDPSolver() : Solver() {};

 ~ThermalUnitDPSolver() override = default;

/** @} ---------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public methods derived from base classes
 *  @{ */

 /// sets the Block that the Solver has to solve
 void set_Block( Block * block ) override;

 /// solves the constructed problem
 int compute( bool changedvars = true ) override;

 /// writes the current solution in the Block
 void get_var_solution( Configuration * solc ) override;

 /// returns a valid lower bound on the optimal objective function value
 OFValue get_lb( void ) override { return( total_cost ); }

 /// returns a valid upper bound on the optimal objective function value
 OFValue get_ub( void ) override { return( total_cost ); }

 /// returns the value of the current solution, if any
 OFValue get_var_value( void ) override { return( total_cost ); }

/** @} ---------------------------------------------------------------------*/
/*------------------- PROTECTED FIELDS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/

 /// builds the graph
 void build_graph( void );

 /// computes the EDPs
 void compute_EDPs( void );

 /// implements the min-path algorithm
 void min_path( void );

 /// computes the variable values and the total cost

 void compute_solutions( void );

/*--------------------------------------------------------------------------*/
/*-------------------- PRIVATE FIELDS OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE TYPES ---------------------------------*/
/*--------------------------------------------------------------------------*/

 /// Stage of the computation
 enum stage {
  start    = 0 ,
  graph_OK = 1 ,
  edps_OK  = 2 ,
  path_OK  = 3 ,
  sol_OK   = 4
  };

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS EDSolver --------------------------------*/
/*--------------------------------------------------------------------------*/
/// base class for the Economic Dispatch Solver
/** EDSolver is a base classe that defines a minimal interface between the
 * ThermalUnitDPSolver and the solvers of the individual Economic Dispatch
 * Problems that give the cost of the arc in the DP. This is geared towards
 * solvers that can cheapily compute all the costs of all the arcs
 * ( h , h ), ( h , h + 1 ), ..., ( h , n ) in one blow, as the DP
 * solver does. However, it being virtual other implementations may be
 * considered. */

class EDPSolver {

/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

 public:

/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/

 EDPSolver( void ) : f_k( 0 ) , f_solver( nullptr ) {}

 virtual ~EDPSolver() = default;

/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/

 /// initializes the ED: h is the starting time, which is fixed

 virtual void initialize( int h , ThermalUnitDPSolver * s ) {
  f_k = h;
  f_solver = s;
  }

/*--------------------------------------------------------------------------*/
 /// compute the cost vector z_h[ h ], z_h[ h + 1 ], ..., z_h[ t - 1 ]
 /** Compute the costs of the arcs ( h , h ), ( h , h + 1 ), ...
  * ( h , t - 1 ), where t is the end of the horizon. These are t - h
  * values, each corresponding to one of the remaining time instants to
  * decide upon (h comprised), and are written in the positions h, h + 1,
  * ..., t - 1 of the vector cost.
  *
  * Note: the last value, that of ( h , t - 1 ), should correspond to the
  * ED in which the unit shuts down right at the end of the time horizon.
  * However, constraining the unit to do that is a bad idea, in that it
  * forces it to enter in a "shutdown trajectory" in the previos time
  * instants, due to the ramp-down constraints, so that the power at t - 1
  * is the right one to stop. This constrains the ED, as opposed to what is
  * the cost of the arc ( h , s ) = "the unit is on from h up until the end
  * of the horizon, and then we'll see" in which the final power of the unit
  * can be whatever (subject to the other constraints). Clearly, choosing
  * the more constrained problem (= larger cost) is never convenient. Hence,
  * we assume ED to directly produce in cost[ t - 1 ] the cost of the less
  * constrained ED corresponding to the arc ( h , s ), so that the cost of
  * both ( h , t - 1 ) and ( h , s ) will be equal. */

 virtual void compute_costs( std::vector< double > & costs ) = 0;

/*--------------------------------------------------------------------------*/
 /// compute optimal power values p_h[ h ], p_h[ h + 1 ], ..., p_h[ k - 1 ]
 /** After compute_costs() have been called once, it is possible to call
  * compute_power_variables( k ) for h <= k <= t - 1 to get the optimal
  * power values corresponding to the arc ( h , k - 1 ); note that this also
  * works for the arc ( h , s ) by passing k = t, as we cheat so that the two
  * correspond to the same ED (see compute_costs()). The optimal values of
  * the power variables are written in the positions h, h + 1, ..., k - 1 of
  * the vector p. The solution depends on k, but this method is typically
  * only called for one particular value of k >= h during the final
  * computation of the optimal solution to the whole 1UC. */
 
 virtual void compute_power_variables( int k , std::vector< double > & p ) = 0;

/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/

 protected:

 /// the ThermalUnitDPSolver using this EDSolver
 ThermalUnitDPSolver * solver;

 int f_k;

 };  // end( class( EDPSolver ) )

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS DPEDSolver ------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// class solving the Economic Dispatch problem via Dynamic Programming
/** DPEDSolver derives from EDPSolver and solves the Economic Dispatch
 * problem by means of a Dynamic Programming approach. */

class DPEDSolver {

/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

 public:

/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/

 DPEDSolver( void ) : EDSolver() = default;

 virtual ~EDPSolver() = default;

/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/

 void initialize( int h , ThermalUnitDPSolver * s ) override;

 void compute_costs( std::vector< double > & costs ) override;

 void compute_power_variables( int k , std::vector< double > & p ) override;

/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/

 protected:

 /// coefficients for a variable of the objective function
 struct coeff_t {
  double alfa;
  double beta;
  double gamma;
  } __attribute__((aligned(32)));

 /// cost coefficients of the objective function
 std::vector< coeff_t > coeffs;

 /// indices for a piece of the (piece-wise) objective function
 struct pos_t {
  int begt;
  int begm;
  } __attribute__((aligned(8)));

 /** For each k = h, ..., n - 1 the vector contains the indices of the pieces
  * of the objective function.  */
 std::vector< pos_t > pos;

 /// unconstrained optimal power values
 std::vector< double > unc_p;

 /// constrained optimal power values
 std::vector< double > con_p;

 int kMax{};

 std::vector< double > m;
 std::vector< int > v;

 }  // end( class( DPEDSolver ) );

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 class node;  // forward declaration of node

/*--------------------------------------------------------------------------*/
 /// an arc

 class arc {
  public:

  arc( void ) : cost1( 0 ) , cost2( 0 ) , tail( nullptr ) , DPS( nullptr ) {}

  ~arc() { delete DPS; }
 
  double cost1;
  double cost2;
  node * tail;
  DPEDSolver * DPS;

  }  // end( class( arc ) )

/*--------------------------------------------------------------------------*/
 /// a node

 class node {
  public:

  node( void ) : lab( 0 ) , pred( nullptr ) {}

  ~node() = default;
  
  double lab;                 /// the label of the node
  node * pred;                /// the predecessor of the node in the path
  std::vector< arc > v_arcs;  /// the Forward Star of the node

  }  // end( class( node ) )

/*--------------------------------------------------------------------------*/
/*---------------------- PRIVATE METHODS OF THE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

 int h_of_node( node * n ) {
  if( n == f_start )  // the source
   return( 0 );
  if( n->DPS )        // an ON-node
   return( std::distance( n , &( v_on_nodes.begin() ) ) );
  // else it must be an OFF-node, this is never called on the destination
  return( std::distance( n , &( v_off_nodes.begin() ) ) );	   
  }

/*--------------------------------------------------------------------------*/

 void load_parameters( void );

/*--------------------------------------------------------------------------*/

 void process_modifications( void );

/*--------------------------------------------------------------------------*/

 void retrieve_term( std::vector< double > & out,
                     const std::vector< double > & in ) const;

/*--------------------------------------------------------------------------*/

 double compute_startup_costs( int t );

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE FIELDS OF THE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

 int time_horizon{};       ///< time horizon
 int init_up_down_time{};  ///< initial up/down time
 int min_up_time{};        ///< minimum up time
 int min_down_time{};      ///< minimum down time
 double initial_power{};   ///< initial power
 int init_t{};

 std::vector< double > startup_costs;
 std::vector< double > delta_ramp_up;
 std::vector< double > delta_ramp_down;
 std::vector< double > min_power;
 std::vector< double > max_power;

 // these two are identical to min_power but we keep them for readability
 std::vector< double > & bound_on = min_power;
 std::vector< double > & bound_down = min_power;

 std::vector< double > quad_term;
 std::vector< double > linear_term;
 std::vector< double > const_term;

 double eps{ 1e-10 };  ///< tolerance

 node f_start;         ///< starting node
 node f_end;           ///< ending node

 std::vector< node > v_on_nodes;   ///< vector of ON nodes
 std::vector< node > v_off_nodes;  ///< vector of OFF nodes
 
 std::vector< double > P;    ///< Power values
 std::vector< int > U;       ///< Commitment values
 std::vector< int > startup; ///< Startup values
 double total_cost{};        ///< Total cost

 int hMin{}; ///< First time step the unit can be turned ON
 int kMin{}; ///< First time step the unit can be turned OFF

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 };  // end( class( ThermalUnitDPSolver ) )

/*--------------------------------------------------------------------------*/

 };  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* ThermalUnitDPSolver.h included */

/*--------------------------------------------------------------------------*/
/*------------------ End File ThermalUnitDPSolver.h ------------------------*/
/*--------------------------------------------------------------------------*/
