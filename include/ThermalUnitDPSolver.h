/*--------------------------------------------------------------------------*/
/*----------------------- File ThermalUnitDPSolver.h -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the ThermalUnitDPSolver class, that solves the
 * ThermalUnitBlock (without primary and secondary reserve variables, so far)
 * using a Dynamic Programming algorithm.
 *
 * \author Claudio Gentile \n
 *         Istituto di Analisi di Sistemi e Informatica "Antonio Ruberti" \n
 *         Consiglio Nazionale delle Ricerche \n
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Niccolo' Iardella \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Claudio Gentile, Antonio Frangioni, Niccolo' Iardella,
 *                      Donato Meoli
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

#include "Solver.h"

#include "ThermalUnitBlock.h"

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)

namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS ThermalUnitDPSolver ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// class for solving a Single Unit Commitment problem with a DP approach
/** The ThermalUnitDPSolver is a Solver for tackling the single-Unit
 * Commitment problem with ramp-up and ramp-down (as well as minimum up-
 * and down-time) constraints and convex quadratic separable objective.
 * The solver uses a Dynamic Programming approach, recasting the problem as
 * a shortest path one on a graph with the following structure
 *
 *       (0,1)  (1,1)  (2,1)  (3,1)  (4,1)  (5,1)  ...  (n-1,1)
 *  (s)                                                           (d)
 *       (0,0)  (1,0)  (2,0)  (3,0)  (4,0)  (5,0)  ...  (n-1,0)
 *
 * n being the length of the time horizon.
 *
 * The fundamental tool for solving this problem is the ability of efficiently
 * solving Economic Dispatch problems that find the min-cost energy production
 * of the unit if it is on for a continuous time interval. In particular we
 * denote by ED( h , k ) the total power cost (but not the fixed and start-up
 * ones, that are computed separately) of the unit if started up exactly at
 * the beginning of time h >= 0 and shut down exactly at the end of time
 * h <= k <= n - 1, i.e., being online for all the time instants h, h + 1,
 * ..., k (note that h = k is possible). Similarly, we denote bu SUC( h , k )
 * the cost of having the unit off from the beginning of h to the end of k
 * and then starting up at k + 1: this is typically easy to compute.
 *
 * The problem is therefore reduced to a Shortest Path between (s) and (d)
 * on the acyclic graph constructed as follows:
 *
 * - Each arc from an ON node ( i , 1 ) to an OFF node ( j , 0 ), for
 *   0 <= i < j <= n - 1, means that the unit is started at the beginning of
 *   time i and shut down at the end of time j - 1, so that it is off at
 *   time j. The cost of the arc is therefore ED( i , j - 1 ). Note that
 *   this problem concerns the power variables p[ h ], p[ h + 1 ], ...
 *   p[ j - 1 ], but also implicitly p[ j ] that will be necessarily
 *   fixed to 0 (as the unit is down at j). However, such an arc exists only
 *   if j - i = number of consecutive periods the unit remains on is
 *   >= min up-time. In particular, j == i + 1, i.e., the unit is started
 *   up and shut down at the end of the same period, is only possible if
 *   min up-time <= 1, i.e., there is no min up-time requirement. Note that
 *   min up-time need necessarily be >= 1 as the unit cannot remain on for
 *   less than one time instant, so min up-time == 0 hardly makes sense.
 *
 * - Each arc from an OFF node ( i , 0 ) to an ON node ( j , 1 ), for
 *   0 <= i < j <= n - 1, means that the unit is shut down at the beginning
 *   of time i and remains down up until the end of time j - 1, then it is
 *   started up at j. Hence, the cost of the arc is the (possibly,
 *   time-variable) start-up costs SUC( i , j - 1 ). This basically fixes
 *   p[ h ] = p[ h + 1 ] = p[ j - 1 ] = 0, but leaves p[ j ] free to be
 *   anything (it will be decided by the outgoing arcs of ( j , 1 )).
 *   Such an arc exists only if j - i = number of consecutive periods the
 *   unit remains off is >= min down-time. In particular, in this case it
 *   would even be possible i == j, i.e., the unit was shut down at the
 *   end of period i - 1 and it immediately re-started at the beginning
 *   of period i, if min down-time == 0, i.e., there is no min down-time
 *   requirement. Note that, unlike for min up-time, in this case the
 *   value 0 in principle makes sense and it is different from the value
 *   1. However, we avoid "vertical" arcs ( i , 0 ) --> ( i , 1 ) since
 *   a solution where the unit is shut down and immediately started up is
 *   never economical w.r.t. one where the unit is never shut down in the
 *   first place, since shutting down entails a "shutdown trajectory" that
 *   brings the unit to the right stopping power that further constrains
 *   the unit (but this is not prohibited if the unit remains on, which
 *   means that remaining on is always at least as cheap).
 *
 * - From each ON node ( i , 1 ) there always is one arc to the destination
 *   d, meaning that the unit remains on in all the time instants between i
 *   and n - 1, and it is *not* shut down at the end of the period. The
 *   cost of this arc is the optimal cost of a "special" ED( i , n - 1 ), 
 *   deciding on all variables  p[ h ], p[ h + 1 ], ..., p[ n - 1 ] and
 *   *not* (implicitly) fixing p[ n - 1 ] = 0 as ED( i , n - 2 ),
 *   corresponding to the arc ( i , 1 ) --> ( n - 1 , 0 ) does. The reason
 *   why ED( i , n - 1 ) is "special" is that, due to the ramp-down
 *   constraints, if the unit has to be down at time k, then it must enter
 *   in a "shutdown trajectory" in the previous time instants, so that the
 *   final power p[ k - 1 ] is the right one to stop. This constrains the ED,
 *   resulting in a higher cost. This means that forcing the shut down at the
 *   end of n - 1 is never economical: it is in principle better to allow the
 *   unit do what it wants (which may comprise autonomously entering in a
 *   shutdown trajectory if this is the optimal thing to do, as this is not
 *   prohibited). This is why the constraints of ED( h , n - 1 ) do not
 *   include the one forcing the power of the unit at the last time instant to
 *   be the shutdown one, unlike for all the other ED( h , k ).
 *
 * - From each OFF node ( i , 0 ) there always is one arc to the destination
 *   d, meaning that the unit remains off in all the time instants between i
 *   and n - 1. Ordinarily with would imply that the unit is started up right
 *   at the beginning of the next horizon of operations, but this is not of
 *   our concern for the current problem. This means that all these arcs have
 *   *zero cost*, as any startup cost will be accounted for in the next
 *   horizon of operations, if any.
 *
 * - From the node (s) there are arcs going to either ( i , 1 ) or ( i, 0 )
 *   nodes depending on the value of init_up_down_time (the amount of time
 *   the unit has been on or of prior to the initial time instant 0), the
 *   min_up_time, min_down_time, initial_power and delta_ramp_dow values,
 *   as applicable. The rules are the following:
 *
 *   = If init_up_down_time > 0, then the unit has been on for
 *     init_up_down_time periods before the initial time instant 0 and it
 *     is at power initial_power at the beginning of time instant 0. Hence,
 *     by the ramp-down constraints, there is a minimum number of time
 *     instants, t_ramp_min, that are necessary to bring the unit to the
 *     power level required to stop (t_ramp_min could be 0 if initial_power
 *     happens to be exactly the right power level). Then, there will be
 *     arcs between s and all nodes ( i , 0 ) for "sufficiently large"
 *     i >= min_node, where:
 *
 *     * if init_up_down_time >= min_up_time, i.e., the unit is already
 *       on since long enough to satisfy the min up-time constraint, then
 *       min_node = t_ramp_min;
 *
 *     * if init_up_down_time < min_up_time, i.e., the unit has to remain
 *       on anyway for at least other ( min_up_time - init_up_down_time )
 *       instants, then
 *       min_node = max( t_ramp_min , min_up_time - init_up_down_time
 *
 *     The cost of each arc ( s , i ) will be ED( 0 , i - 1 ), where these
 *     ED are also "special" in the sense that, for the sake of ramp-up and
 *     ramp-down constraints, the initial_power value is used as reference.
 *     Note however the very special case where init_up_down_time >=
 *     min_up_time and t_ramp_min == 0, i.e., min_node == 0. This
 *     corresponds to the case where the unit is "on but on the brink of
 *     shutting down" at the beginning of the time horizon. This means that
 *     there will be *both* an arc s --> ( 0 , 1 ) saying "the unit is on
 *     and will remain on for a while", *and* an arc s --> ( 0 , 0 ) saying
 *     "the unit is on but I'll shut it down immediately and will remain
 *     off for w while".
 *
 *     Note that "sufficiently large" i includes i == n, i.e., the arc
 *     s --> d corresponding to "the unit was on at the beginning and
 *     remains on for the whole period".
 *
 *   = If init_up_down_time <= 0, then the unit has been off for
 *     init_up_down_time periods before the initial time instant 0; note
 *     that the value 0 is included, meaning "the unit has just been
 *     shut down when the time begins". Then, there will be arcs between
 *     s and all nodes ( i , 1 ) for "sufficiently large"
 *     i >= max( min down-time + init_up_down_time , 0 ). That is, we
 *     wait for the remaining  min down-time - ( - init_up_down_time )
 *     time periods (if any) required by the min down-time constraint
 *     before allowing the unit to be started up again. These arcs will
 *     have cost SUC( 0 , i ) since the unit will be started up at i.
 *     Note the that if min down-time == 0, i.e., there is no min
 *     down-time requirement, this means that i == 0 is always possible,
 *     i.e., the unit is started up immediately at the beginning. This
 *     potentially yields the "double strange" case where
 *     init_up_down_time == 0, i.e., the unit had just been shut down
 *     and it is immediately restarted. While this is not economical, we
 *     cannot (and have no reason to) avoid it, as in the case of
 *     OFF -> ON arcs, because the decision to shut down the unit right
 *     at the end of the previous interval (encoded by init_up_down_time
 *     == 0) is not in our hands as it has been taken before out time has
 *     come. Yet, this is not impossible not logically contradictory, so
 *     there is no problem (and min down-time == 0 is unlikely anyway).
 *
 *     Note that "sufficiently large" i includes i == n, i.e., the arc
 *     s --> d corresponding to "the unit was off at the beginning and
 *     remains off for the whole period". Like all arcs OFF -> d, this
 *     does not really imply that the unit will necessarily be restarted
 *     immediately after, and anyway this is outside the boundaries of
 *     the current problem; hence, this arc has *zero cost*.
 *
 * ThermalUnitDPSolver first builds the graph, then uses one EDSolver for
 * each ON node (comprised s if the unit is on at the beginning, and therefore
 * is it equivalent to a ON node) to solve EDs to compute the arc costs, then 
 * uses a( acyclic) min-path algorithm to solve the problem. */

class ThermalUnitDPSolver : public Solver
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*------------------------------ PUBLIC TYPES ------------------------------*/
/*--------------------------------------------------------------------------*/

 static constexpr auto TUDPINF = Inf< double >();  ///< the INF value

 using Index = Block::Index;

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and destructor
 * @{ */

 ThermalUnitDPSolver( void ) : Solver() {};

 ~ThermalUnitDPSolver() override = default;

/** @} ---------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public methods derived from base classes
 * @{ */

 /// sets the Block that the Solver has to solve
 void set_Block( Block * block ) override;

 /// solves the constructed problem
 int compute( bool changedvars = true ) override;

 /// tells whether a solution is available
 bool has_var_solution( void ) override { return( f_end.pred ); }

 /// writes the current solution in the Block
 void get_var_solution( Configuration * solc ) override;

 /// returns a valid lower bound on the optimal objective function value
 OFValue get_lb( void ) override { return( f_end.lab ); }

 /// returns a valid upper bound on the optimal objective function value
 OFValue get_ub( void ) override { return( f_end.lab ); }

 /// returns the value of the current solution, if any
 OFValue get_var_value( void ) override { return( f_end.lab ); }

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 protected:

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

 /// stage of the computation
 enum stage_value
 {
  start = 0 ,
  graph_OK = 1 ,
  edps_OK = 2 ,
  path_OK = 3 ,
  sol_OK = 4
 };

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS EDSolver --------------------------------*/
/*--------------------------------------------------------------------------*/
 /// base class for the Economic Dispatch Solver
 /** EDSolver is a base class that defines a minimal interface between the
  * ThermalUnitDPSolver and the solvers of the individual Economic Dispatch
  * Problems that give the cost of the arc in the DP. This is geared towards
  * solvers that can cheaply compute all the costs of all the arcs
  * ( h , h ), ( h , h + 1 ), ..., ( h , n ) in one blow, as the DP
  * solver does. However, it being virtual other implementations may be
  * considered. */

 class EDSolver
 {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

  public:

/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/

  EDSolver( Index h , ThermalUnitDPSolver * s )
   : f_h( h ) , f_solver( s ) {}

  virtual ~EDSolver() = default;

/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/

  /// compute the cost vector z_h[ h ], ..., z_h[ n - 1 ], z_h[ d ]
  /** Compute the costs of the arcs ( h , h ), ( h , h + 1 ), ...
   * ( h , n - 1 ), ( h , d ), where t is the end of the horizon and d is the
   * end node. The picture for h == 2 and n == 6 is
   *
   *       (0,1)   (1,1)   (2,1) -------+-------+-------+
   *                             \       \       \       \-> (t)
   *       (0,0)   (1,0)   (2,0)   (3,0)   (4,0)   (5,0)
   *
   * That is, these are t - h [ = 6 - 2 = 4 ] values corresponding to the
   * costs of the arcs in the DP graph of type ( h, h + 1 ), ( h, h + 2 ),
   * ..., ( h, n - 1 ), and finally the special arc ( h, d ) [ ( 2, 3 ),
   * ( 2, 4 ), ( 2, 5 ), ( 2, t ) ]. These are written in the positions h,
   * h + 1, ..., t - 1 [ 2 , 3 , 4 , 5 ] of the vector cost.
   *
   * The cost of each arc ( h , k ) for h < k <= n - 1 is ED( h , k - 1 ),
   * corresponding to the fact that the unit remains on from h to k - 1
   * included, but it is off at k. The cost of the special arc ( h , t )
   * corresponds to a "special" ED( h , n - 1 ) in which the unit remains
   * on from h to the end of the time horizon, comprised the last instant.
   * The difference is that in this last ED we do *not* assume the unit will
   * be shut down at n, as this is outside of the time horizon and whatever
   * happens to the unit then is of no concern here.
   *
   * More specifically, the point is that, due to the ramp-down constraints,
   * if the unit has to be down at time k, then it must enter in a "shutdown
   * trajectory" in the previous time instants, so that the power at k - 1
   * is the right one to stop. This constrains the ED, resulting in a higher
   * cost. This means that forcing the shut down at the end of n - 1 is never
   * economical: it is in principle better to allow the unit do what it wants
   * (which may comprise autonomously entering in a shutdown trajectory if
   * this is the optimal thing to do, as this is not prohibited). This is
   * why the constraints of the "special" ED( h , n - 1 ) do not include the
   * one forcing the power of the unit at the last time instant to be the
   * shutdown one, unlike for all the other ED( h , k ). */

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

  virtual void compute_power_variables( Index k ,
                                        std::vector< double > & p ) = 0;

/*--------------------------------------------------------------------------*/
/*-------------------- PRIVATE FIELDS OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

  protected:

  /// the ThermalUnitDPSolver using this EDSolver
  ThermalUnitDPSolver * f_solver;

  Index f_h;  ///< the initial time instant for this EDSolver

 };  // end( class( EDSolver ) )

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS DPEDSolver ------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
 /// class solving the Economic Dispatch problem via Dynamic Programming
 /** DPEDSolver derives from EDSolver and solves the Economic Dispatch
  * problem by means of a Dynamic Programming approach. */

 class DPEDSolver : public EDSolver
 {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

  public:

/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/

  DPEDSolver( Index h , ThermalUnitDPSolver * s );

  virtual ~DPEDSolver() = default;

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC METHODS OF THE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

  void compute_costs( std::vector< double > & costs ) override;

  void compute_power_variables( Index k , std::vector< double > & p ) override;

/*--------------------------------------------------------------------------*/
/*-------------------- PRIVATE FIELDS OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

  protected:

  /// coefficients for a variable of the objective function
  struct coeff_t {
      double alfa;
      double beta;
      double gamma;
  };

  /// cost coefficients of the objective function
  std::vector< coeff_t > coeffs;

  /// indices for a piece of the (piece-wise) objective function
  struct pos_t {
      int begt;
      int begm;
  };

  /** For each k = h, ..., n - 1 the vector contains the indices of the pieces
   * of the objective function. */
  std::vector< pos_t > pos;

  /// unconstrained optimal power values
  std::vector< double > unc_p;

  /// constrained optimal power values
  std::vector< double > con_p;

  std::vector< double > m;
  std::vector< int > v;

 };  // end( class( DPEDSolver ) );

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 class node;  // forward declaration of node

/*--------------------------------------------------------------------------*/
 /// an arc

 class arc
 {
  public:

  arc( void ) : cost1( 0 ) , cost2( 0 ) , tail( nullptr ) {}

  ~arc() = default;

  double cost1;  ///< the now-power-dependent part of the cost (fixed, SUC)
  double cost2;  ///< the power-dependent part of the cost
  node * tail;   ///< (pointer to) the tail node

 };  // end( class( arc ) )

/*--------------------------------------------------------------------------*/
 /// a node

 class node
 {
  public:

  node( void ) : lab( 0 ) , pred( nullptr ) , DPS( nullptr ) {}

  ~node() { delete( DPS ); }

  double lab;                 ///< the label of the node
  node * pred;                ///< the predecessor of the node in the path
  EDSolver * DPS;             ///< the Economic Dispatch solver of the node
  std::vector< arc > v_arcs;  ///< the Forward Star of the node

 };  // end( class( node ) )

/*--------------------------------------------------------------------------*/
/*---------------------- PRIVATE METHODS OF THE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

 Index h_of_node( node * n ) {
  if( n == &f_start )  // the source should be "-1", but we make it 0
   return( 0 );
  if( n == &f_end )    // the destination
   return( time_horizon );
  if( n->DPS )          // an ON-node
   return( n - v_on_nodes.data() );
  // else it must be an OFF-node, this is never called on the destination
  return( n - v_off_nodes.data() );
 }

/*--------------------------------------------------------------------------*/

 // reset label and predecessor of a node

 static void init_node( node & nde ) {
  nde.lab = TUDPINF;
  nde.pred = nullptr;
 }

/*--------------------------------------------------------------------------*/

 // do the scanning of the forward star of a node

 static void process_node( node & nde ) {
  for( auto a : nde.v_arcs ) {
   const auto nl = nde.lab + a.cost1 + a.cost2;
   if( ( a.tail )->lab > nl ) {
    ( a.tail )->lab = nl;
    ( a.tail )->pred = &nde;
   }
  }
 }

/*--------------------------------------------------------------------------*/

 void load_parameters( void );

/*--------------------------------------------------------------------------*/

 void process_modifications( void );

 // returns true if everything need be reset
 bool guts_of_process_modifications( const p_Mod mod );

/*--------------------------------------------------------------------------*/

 void retrieve_term( std::vector< double > & out ,
                     const std::vector< double > & in ) const;

/*--------------------------------------------------------------------------*/

 double compute_startup_costs( Index h , Index k ) {
  // one day a time-dependent SUC formula may be easily implemented here
  if( startup_costs.empty() )
   return( 0 );
  return( startup_costs[ k ] );
 }

/*--------------------------------------------------------------------------*/
/*-------------------- PRIVATE FIELDS OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 Index time_horizon;       ///< time horizon
 int init_up_down_time;    ///< initial up/down time (it can be < 0)
 Index min_up_time;        ///< minimum up time
 Index min_down_time;      ///< minimum down time
 double initial_power;     ///< initial power
 Index t_init;             ///< the first instant in which commitment is free

 std::vector< double > startup_costs;
 std::vector< double > delta_ramp_up;
 std::vector< double > delta_ramp_down;
 std::vector< double > min_power;
 std::vector< double > max_power;
 std::vector< double > bound_on;
 std::vector< double > bound_down;

 std::vector< double > quad_term;
 std::vector< double > linear_term;
 std::vector< double > const_term;

 double eps{ 1e-10 };              ///< tolerance

 char stage;                       ///< what has been computed

 node f_start;                     ///< starting node
 node f_end;                       ///< ending node

 std::vector< node > v_on_nodes;   ///< vector of ON nodes
 std::vector< node > v_off_nodes;  ///< vector of OFF nodes

 std::vector< double > P;          ///< power values
 std::vector< bool > U;            ///< commitment values

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/
/*---------------------- PRIVATE METHODS OF THE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/



};  // end( class( ThermalUnitDPSolver ) )

};  // end( namespace SMSpp_di_unipi_it )

#endif  /* ThermalUnitDPSolver.h included */

/*--------------------------------------------------------------------------*/
/*------------------ End File ThermalUnitDPSolver.h ------------------------*/
/*--------------------------------------------------------------------------*/
