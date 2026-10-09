/*--------------------------------------------------------------------------*/
/*----------------------- File ThermalUnitDPSolver.h -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the ThermalUnitDPSolver class, that solves the
 * ThermalUnitBlock, comprised its spinning-reserve variables, using a
 * Dynamic Programming algorithm.
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

#include <memory>

#include "Solver.h"

#include "ThermalUnitBlock.h"

#include "ThermalUnitDPSolverBase.h"

#define TUDPS_PARALLEL 1
/* If TUDPS_PARALLEL > 0, the (independent) per-ON-node Economic Dispatch
 * solves in compute_EDPs() can be run in parallel with FastFlow; the actual
 * number of workers is controlled at run time by the intMaxThread parameter
 * (0 = use all available cores, 1 = serial, k = k workers). Set to 0 in build
 * setups where FastFlow is not available (e.g. the plain makefiles). */

#ifndef TUDPS_PAR_MIN_N
 #define TUDPS_PAR_MIN_N 768
#endif
/* Default time horizon below which compute_EDPs() stays serial even with
 * intMaxThread != 1: the thread-dispatch overhead is not amortised on short
 * horizons (benchmarks on 8 cores put the serial/parallel break-even around
 * 600 time steps, so this is set conservatively above it). It is only the
 * default of the run-time intParMinN parameter (see below); set intParMinN in
 * the ComputeConfig to override it per solver instance (0 forces the parallel
 * path on every horizon), or override this compile-time default with
 * -DTUDPS_PAR_MIN_N=<n>. */

#if TUDPS_PARALLEL
namespace ff { class ParallelFor; }  // forward declaration (FastFlow)
#endif

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
/// dynamic programming Solver of a ThermalUnitBlock over its on-intervals
/** The ThermalUnitDPSolver solves the single-unit commitment problem of a
 * ThermalUnitBlock, i.e., the problem of ThermalUnitDPSolverBase (cf. the
 * comment of ThermalUnitDPSolverBase.h for the model, the notation and the
 * reserves). We formulate it as a shortest path on an acyclic graph whose
 * arcs are the on-intervals and the off-intervals of the unit, and we price
 * each on-interval by the economic dispatch of the unit while it is on. The
 * graph has a source \f$ s \f$, a destination \f$ d \f$ and, for each instant
 * \f$ i \in \mathcal{T} \f$, an on node \f$ ( i , 1 ) \f$ (the unit starts up
 * at \f$ i \f$) and an off node \f$ ( i , 0 ) \f$ (the unit shuts down at
 * \f$ i \f$, i.e., it is on at \f$ i - 1 \f$ and off at \f$ i \f$):
 * \verbatim
 *        (0,1)  (1,1)  (2,1)  ...  (T-1,1)
 *   s                                        d
 *        (0,0)  (1,0)  (2,0)  ...  (T-1,0)
 * \endverbatim
 * With \f$ \tau^+ \f$ and \f$ \tau^- \f$ at least 1, the arcs are the
 * following.
 *
 * - From \f$ ( i , 1 ) \f$ to \f$ ( j , 0 ) \f$ for
 *   \f$ i + \tau^+ \leq j \leq T - 1 \f$: the unit is on at
 *   \f$ i , \ldots , j - 1 \f$ and off at \f$ j \f$; the cost is
 *   \f[
 *     \mathrm{ED}( i , j - 1 ) + \sum_{ t = i }^{ j - 1 } c_t + c^{sd}_j .
 *     \tag{1}
 *   \f]
 *
 * - From \f$ ( i , 1 ) \f$ to \f$ d \f$: the unit is on from \f$ i \f$ to
 *   the end of the horizon, and the cost is
 *   \f$ \mathrm{ED}( i , T - 1 ) + \sum_{ t = i }^{ T - 1 } c_t \f$.
 *
 * - From \f$ ( i , 0 ) \f$ to \f$ ( j , 1 ) \f$ for
 *   \f$ i + \tau^- \leq j \leq T - 1 \f$: the unit is off at
 *   \f$ i , \ldots , j - 1 \f$ and starts up at \f$ j \f$, at the cost
 *   \f$ c^{su}_j \f$, which depends on the instant of the start-up but not on
 *   how long the unit has been off.
 *
 * - From \f$ ( i , 0 ) \f$ to \f$ d \f$, at no cost: the unit stays off to
 *   the end of the horizon.
 *
 * - From \f$ s \f$, if the unit is off before the horizon
 *   (\f$ \tau_0 \leq 0 \f$), to \f$ ( j , 1 ) \f$ for
 *   \f$ t_0 \leq j \leq T - 1 \f$, at the cost \f$ c^{su}_j \f$, and to
 *   \f$ d \f$ at no cost. Here \f$ t_0 = \min\{ T , \max\{ 0 , \tau^- +
 *   \tau_0 \} \} \f$ is the first instant at which the unit may start up.
 *
 * - From \f$ s \f$, if the unit is on before the horizon
 *   (\f$ \tau_0 > 0 \f$), to \f$ ( j , 0 ) \f$ for \f$ \max\{ t_0 , t^{sd}_0
 *   \} \leq j \leq T - 1 \f$, with the cost (1) for \f$ i = 0 \f$, and to
 *   \f$ d \f$ with the cost of the arc from an on node to \f$ d \f$ for
 *   \f$ i = 0 \f$. Here \f$ t_0 = \min\{ T , \max\{ 0 , \tau^+ - \tau_0 \}
 *   \} \f$, and \f$ t^{sd}_0 \f$ is the least \f$ k \f$ such that
 *   \f$ p_{-1} - \sum_{ r = 0 }^{ k - 1 } \Delta^-_r \f$ is below
 *   \f$ P^{sd}_k \f$ (plus a tolerance), i.e., the first instant at which the
 *   unit can be off after ramping down from InitialPower. For
 *   \f$ j \geq 1 \f$ this is only a pruning, since the descent is enforced by
 *   the economic dispatch of the arc, which has
 *   \f$ p^{ac}_{j-1} \leq P^{sd}_j \f$ (cf. (2)). For \f$ j = 0 \f$, instead,
 *   it is the only condition, because the arc to \f$ ( 0 , 0 ) \f$ (i.e., the
 *   shut-down of the unit at instant 0 at the cost \f$ c^{sd}_0 \f$) has an
 *   empty dispatch; this arc exists if \f$ \tau_0 \geq \tau^+ \f$ and
 *   \f$ p_{-1} \leq P^{sd}_0 \f$, whether or not DeltaRampDown is given, as
 *   in all the formulations of ThermalUnitBlock (see
 *   ThermalUnitBlock::generate_abstract_constraints()). Without DeltaRampDown
 *   the descent is not computed, and the unit can be off at 1 if it cannot at
 *   0.
 *
 * An arc whose interval contains an instant at which the commitment Variable
 * is fixed to the opposite state does not exist [see load_fixings()], and
 * therefore the fixings are honored. In (1), \f$ \mathrm{ED}( i , k ) \f$ is
 * the optimal value of the economic dispatch of the run
 * \f$ i , \ldots , k \f$,
 * \f[
 *   \mathrm{ED}( i , k ) = \min \Bigl\{ \, \sum_{ t = i }^{ k } \bigl(
 *     a_t ( p^{ac}_t )^2 + b_t p^{ac}_t + g_t( p^{ac}_t , \mathcal{H}_t(
 *     p^{ac}_t , p^{ac}_{t-1} ) ) \bigr) \, : \,
 *     p^{ac}_t \in [ P^{mn}_t , P^{mx}_t ] \;\, ( t = i , \ldots , k ) ,
 *     \;\, -\Delta^-_t \leq p^{ac}_t - p^{ac}_{t-1} \leq \Delta^+_t \;\,
 *     ( t = i + 1 , \ldots , k ) \, \Bigr\} \tag{2}
 * \f]
 * with \f$ g_t \f$ and \f$ \mathcal{H}_t \f$ those of ThermalUnitDPSolverBase
 * (at \f$ t = i \f$ the band has no ramp terms, save for the run of the unit
 * on before the horizon, which has \f$ p^{ac}_{-1} = p_{-1} \f$; see below).
 * In addition, \f$ p^{ac}_i \leq P^{su}_i \f$ if \f$ ( i , 1 ) \f$ is a
 * start-up, \f$ -\Delta^-_0 \leq p^{ac}_0 - p_{-1} \leq \Delta^+_0 \f$ if the
 * run is the one of the unit on before the horizon, and
 * \f$ p^{ac}_k \leq P^{sd}_{k+1} \f$ if the unit is off at
 * \f$ k + 1 \leq T - 1 \f$ (the run that lasts to the end of the horizon has
 * no such row). It is solved for all \f$ k \f$ at once by the forward
 * recursion on the convex piecewise quadratic functions \f$ z_{i,k} \f$ of
 * the power at \f$ k \f$,
 * \f{align*}{
 *   z_{i,i}( p ) &= a_i p^2 + b_i p + g_i( p , \min\{ p - P^{mn}_i ,
 *     K_i - p \} ) \, , \\
 *   z_{i,k}( p ) &= a_k p^2 + b_k p + \min \bigl\{ \, z_{i,k-1}( q ) +
 *     g_k( p , \mathcal{H}_k( p , q ) ) \, : \, p - \Delta^+_k \leq q \leq
 *     p + \Delta^-_k \, \bigr\} \, , \quad k > i \, ,
 * \f}
 * on the domain of the bounds above. For the run of the unit on before the
 * horizon one has instead \f$ z_{0,0}( p ) = a_0 p^2 + b_0 p + g_0( p ,
 * \mathcal{H}_0( p , p_{-1} ) ) \f$, since the reserve of instant 0 is
 * bounded by the ramp left after the move from InitialPower. Then,
 * \f$ \mathrm{ED}( i , k ) = \min_p z_{i,k}( p ) \f$, where for
 * \f$ k < T - 1 \f$ the minimum is over \f$ p \leq P^{sd}_{k+1} \f$ and the
 * cap \f$ K_k \f$ of the reserve band of instant \f$ k \f$ is lowered to
 * \f$ P^{sd}_{k+1} \f$ (the last step of the recursion is redone with it when
 * it matters [see ThermalUnitDPSolverBase]). When no reserve cost is negative
 * \f$ g \equiv 0 \f$, and the recursion is the plain sliding minimum of the
 * quadratic costs. The labels of the shortest path are then computed in the
 * topological order \f$ s , ( 0 , 1 ) , ( 0 , 0 ) , ( 1 , 1 ) , \ldots ,
 * d \f$ by \f$ \Lambda( s ) = 0 \f$ and \f$ \Lambda( v ) = \min \{ \Lambda( w
 * ) + \mbox{cost}( w , v ) \} \f$ over the arcs into \f$ v \f$. Finally, the
 * schedule is recovered from the predecessors and from the minimizers of the
 * recursion, and the reserves are those of (2) of ThermalUnitDPSolverBase at
 * the recovered schedule.
 *
 * The arcs enumerate exactly the commitments that satisfy the minimum up and
 * down times from the initial state and the fixings. In fact, an on-arc that
 * ends within the horizon lasts at least \f$ \tau^+ \f$ instants, and one
 * that reaches \f$ d \f$ may be shorter; these are the commitments that the
 * rows (2) of ThermalUnitBlock allow (a run shorter than \f$ \tau^+ \f$ that
 * ends within the horizon violates the row at its last instant, or at
 * \f$ T - 1 \f$), and the same holds for the off-arcs and the rows (3). Since
 * the economic dispatch of each on-interval is solved exactly over the
 * continuous power (no discretization is involved), the value is the optimal
 * value of the rows of ThermalUnitBlock, with the exceptions stated in
 * ThermalUnitDPSolverBase.h (a ReferenceSchedule and the fixed Variable other
 * than the commitment and the design, which are refused, and the parts of the
 * transitions computed by interpolation). The reactive power, when the unit
 * has it, is separable and is added as a constant, while its
 * commitment-dependent part is added to \f$ c_t \f$. If the unit has an
 * InvestmentCost, the problem above is the one of the unit built, and the
 * unit is built if its value plus the current coefficient of the design
 * variable is not positive (or if a fixing forces it); the value is then
 * multiplied by the scale factor. Note that the economic dispatches of the
 * different on nodes are independent, and hence they are solved in parallel
 * when the horizon is long enough [see intParMinN]. */

class ThermalUnitDPSolver : public ThermalUnitDPSolverBase
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

 ThermalUnitDPSolver( void );  // defined in the .cpp (pimpl'd FastFlow)

 ~ThermalUnitDPSolver() override;  // defined in the .cpp (pimpl'd FastFlow)

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

/*--------------------------------------------------------------------------*/
 /// returns the schedule the DP has found as a ThermalUnitBlockSolution
 /** Returns the schedule the dynamic programming has found as a
  * ThermalUnitBlockSolution [see ThermalUnitBlock.h], filled straight out of
  * the data structures of the Solver rather than by writing it in the
  * Variable of the ThermalUnitBlock and having it read back from there: no
  * abstract representation is therefore required to exist, and the
  * ThermalUnitBlock is not written into at all, hence it is not lock()-ed
  * and any number of Solver attached to it can produce their own Solution at
  * the same time.
  *
  * The active power, the commitment and the dimensioning variable are always
  * saved; the spinning reserves and the reactive power only if the unit has
  * them, which is the same rule get_var_solution() follows by only writing
  * the Variable that exist. */

 [[nodiscard]] Solution * get_Solution( Configuration * solc = nullptr )
  override;

/*--------------------------------------------------------------------------*/
 /// the Solution is filled from the data of the DP, not from the Variable

 [[nodiscard]] bool is_get_Solution_physical( void ) const override {
  return( true );
  }

 /// returns a valid lower bound on the optimal objective function value
 OFValue get_lb( void ) override { return( scaled_value( f_end.lab ) ); }

 /// returns a valid upper bound on the optimal objective function value
 OFValue get_ub( void ) override { return( scaled_value( f_end.lab ) ); }

 /// returns the value of the current solution, if any
 OFValue get_var_value( void ) override {
  return( scaled_value( f_end.lab ) );
  }

/*--------------------------------------------------------------------------*/

 /// extends Solver::int_par_type_S with the ThermalUnitDPSolver parameters
 enum int_par_type_TUDPS {
  intParMinN = intLastAlgPar , ///< min time horizon for the parallel DP path
  /**< The per-ON-node Economic Dispatch sweep in compute_EDPs() is run in
   * parallel (over the FastFlow workers set by intMaxThread) only when the
   * time horizon is at least intParMinN; below it the thread-dispatch
   * overhead is not amortised and the solver stays serial. Defaults to the
   * compile-time TUDPS_PAR_MIN_N; set it to 0 to force the parallel path on
   * every horizon. */
  intLastAlgParTUDPS ///< 1st allowed new int parameter for derived classes
  };

 using Solver::set_par;  // keep the other set_par() overloads visible

 /// honored parameters: intMaxThread (workers) and intParMinN (threshold)
 void set_par( idx_type par , int value ) override {
  if( par == intMaxThread ) { f_max_thread = value; return; }
  if( par == intParMinN ) { f_par_min_n = value; return; }
  Solver::set_par( par , value );
  }

 [[nodiscard]] int get_int_par( idx_type par ) const override {
  if( par == intMaxThread ) return( f_max_thread );
  if( par == intParMinN ) return( f_par_min_n );
  return( Solver::get_int_par( par ) );
  }

 [[nodiscard]] idx_type get_num_int_par( void ) const override {
  return( Solver::get_num_int_par() + intLastAlgParTUDPS - intLastAlgPar );
  }

 [[nodiscard]] int get_dflt_int_par( idx_type par ) const override {
  return( par == intParMinN ? TUDPS_PAR_MIN_N
                            : Solver::get_dflt_int_par( par ) );
  }

 [[nodiscard]] idx_type int_par_str2idx( const std::string & name )
  const override {
  return( name == "intParMinN" ? intParMinN
                               : Solver::int_par_str2idx( name ) );
  }

 [[nodiscard]] const std::string & int_par_idx2str( idx_type idx )
  const override {
  static const std::string name = "intParMinN";
  return( idx == intParMinN ? name : Solver::int_par_idx2str( idx ) );
  }

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

  /// compute the dispatch costs of all the runs that start at h
  /** Writes in \p costs[ k ], for \f$ k = h , \ldots , T - 1 \f$, the
   * value \f$ \mathrm{ED}( h , k ) \f$ of (2) of the class comment, i.e.,
   * the power-dependent cost of the arc from \f$ ( h , 1 ) \f$ to
   * \f$ ( k + 1 , 0 ) \f$ for \f$ k < T - 1 \f$, with
   * \f$ p_k \leq P^{sd}_{k+1} \f$, and of the arc from \f$ ( h , 1 ) \f$
   * to \f$ d \f$ for \f$ k = T - 1 \f$, with no shut-down limit since the
   * unit is not shut down within the horizon. For the source of a unit on
   * before the horizon the run starts at \f$ h = 0 \f$ from InitialPower,
   * rather than from a start-up. */

  virtual void compute_costs( std::vector< double > & costs ) = 0;

/*--------------------------------------------------------------------------*/
  /// compute the optimal powers of the run from h to k
  /** After compute_costs(), writes in \p p[ h ], ..., \p p[ k ] an optimal
   * dispatch of the run \f$ h , \ldots , k \f$, i.e., a minimizer of
   * \f$ \mathrm{ED}( h , k ) \f$, by walking the recursion of the class
   * comment backwards from the minimizer of \f$ z_{h,k} \f$. */

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
 /// the economic dispatches of the runs from one instant, by recursion
 /** DPEDSolver computes \f$ \mathrm{ED}( h , k ) \f$ for all \f$ k \f$ by
  * the forward recursion on \f$ z_{h,k} \f$ of the class comment. */

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

  void compute_power_variables( Index k ,
                                std::vector< double > & p ) override;

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

  /// size multiplier of the buffers: 1 without reserves, 2 with them, the
  /// reserve term splitting the pieces of the value functions
  Index f_rmul{ 1 };

  /// the recursion of the class comment when a reserve cost is negative
  /** Used by compute_costs() when some reserve cost is negative: it carries
  * \f$ z_{h,k} \f$ as a PQFun and computes each step by
  * ThermalUnitDPSolverBase::sliding_min_corr(), which includes the reserve
  * term \f$ g_k \f$ with the ramp terms of its band. It fills the same
  * \p costs[ k ] \f$ = \mathrm{ED}( h , k ) \f$. */
  void compute_costs_reserve( std::vector< double > & costs );

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

  double cost1;  ///< the part of the cost not on the power (fixed, start-up)
  double cost2;  ///< the power-dependent part of the cost
  node * tail;   ///< (pointer to) the tail node

 };  // end( class( arc ) )

/*--------------------------------------------------------------------------*/
 /// a node

 class node
 {
  public:

  node( void ) : lab( 0 ) , pred( nullptr ) , DPS( nullptr ) {}

  ~node() = default;

  double lab;                 ///< the label of the node
  node * pred;                ///< the predecessor of the node in the path
  EDSolver * DPS;             ///< the ED solver of the node (non-owning: the
                              ///< solvers are owned by the solver's pool)
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

 // reset a node for a fresh build_graph(), reusing its arc storage and the
 // pooled ED solvers; lab == 0 means "not yet proved reachable from s"

 static void reset_node( node & nde ) {
  nde.lab = 0;
  nde.pred = nullptr;
  nde.DPS = nullptr;
  nde.v_arcs.clear();  // keeps the allocated capacity
 }

/*--------------------------------------------------------------------------*/

 // return the pooled ED solver for the ON node at instant i, allocating it
 // the first time (and reusing it, buffers included, on later re-solves)

 DPEDSolver * get_on_ed( Index i );

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

 /// read the fixed status of the Variable of the ThermalUnitBlock
 /** Reads which commitment (and design) Variable are fixed, translating the
  * commitment fixings into the nxt_off / nxt_on tables that build_graph()
  * uses to prune the incompatible arcs; throws std::logic_error if any
  * Variable that the DP cannot honor (active power, reserves, start-up,
  * shut-down, reactive, or a Variable of an extended formulation [see
  * fixed_extended_variable()]) is fixed. */
 void load_fixings( void );

/*--------------------------------------------------------------------------*/

 void process_modifications( void );

 // returns true if everything needs to be reset
 bool guts_of_process_modifications( const p_Mod mod );

/*--------------------------------------------------------------------------*/

 /// recovers the schedule the DP has found
 /** Writes the active power \p p, the commitment \p u, the reserves \p pr
  * and \p sr at the recovered schedule, the reactive power \p q, and
  * whether the unit is \p built. */
 void recover_schedule( std::vector< double > & p ,
                        std::vector< double > & u ,
                        std::vector< double > & pr ,
                        std::vector< double > & sr ,
                        std::vector< double > & q ,
                        bool & built ) const;

 /// true if some reserve cost is negative
 /** When true the economic dispatch takes compute_costs_reserve(), when
  * false the plain sweep of the quadratic costs. */
 bool reserve_rewarded( void ) const;

/*--------------------------------------------------------------------------*/

 double compute_startup_costs( Index h , Index k ) {
  // the cost of the off-run [ h , k ), which starts up at k: it depends on
  // the instant of the start-up only
  if( startup_costs.empty() )
   return( 0 );
  return( startup_costs[ k ] );
 }

/*--------------------------------------------------------------------------*/

 double compute_shutdown_costs( Index h , Index k ) {
  // the cost of the ON-run [ h , k ), which shuts down at k; a run that
  // reaches the end of the horizon never shuts down and pays nothing
  if( shutdown_costs.empty() )
   return( 0 );
  return( shutdown_costs[ k ] );
 }

/*--------------------------------------------------------------------------*/
/*-------------------- PRIVATE FIELDS OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 // NOTE: the data loaded from the ThermalUnitBlock (power/ramp/cost bounds,
 // the spinning-reserve factors and prices, the design and commitment-fixing
 // data, the scalar parameters), the PieceQuad / PQFun value-function type
 // and its operations, and the whole reserve model (reserve_alloc /
 // reserve_alloc_band / reserve_reward / build_reserve_discount /
 // reserve_corr_argmin / sliding_min_corr / ...) live in the base class
 // ThermalUnitDPSolverBase and are inherited. Only the solver-specific
 // reactive increment, the on/off graph, the ED-solver pool and the output
 // are declared here.

 // reactive power q[t] in [reactive_min[t] + reactive_min_on[t] u[t],
 // reactive_max[t] + reactive_max_on[t] u[t]], separable from the
 // active-power DP but gated by the commitment u[t]; reactive_linear_term[t]
 // is its dualized linear cost (empty if the unit has no reactive power).
 // The off-box reward is priced as the constant Q_star; when the box is
 // gated (the _on vectors are non-empty) the per-on-period increment
 // reactive_delta[t] = r_on - r_off is added to the fixed cost of every
 // on-period. fill_reactive_delta() (re)builds reactive_delta from the
 // current reactive price and boxes.
 std::vector< double > reactive_delta;

 void fill_reactive_delta( void );

 /// the reserve terms of the instants with no on predecessor, as PQFuns
 /** eff_disc[ t ] and eff_disc_su[ t ] are
  * ThermalUnitDPSolverBase::build_reserve_discount() at \p t with the cap
  * \f$ P^{mx}_t \f$ and \f$ P^{su}_t \f$, respectively, computed once per
  * compute_EDPs(); both are empty when no reserve cost is negative. */
 std::vector< PQFun > eff_disc;
 std::vector< PQFun > eff_disc_su;  ///< the start-up variant of eff_disc

 char stage;                       ///< what has been computed

 node f_start;                     ///< starting node
 node f_end;                       ///< ending node

 std::vector< node > v_on_nodes;   ///< vector of ON nodes
 std::vector< node > v_off_nodes;  ///< vector of OFF nodes

 /// pool of ED solvers, one slot per time instant, reused across re-solves
 /// so that the (dominant) cost of allocating their O(n) buffers is paid
 /// only once per time horizon rather than at every structural re-solve
 std::vector< std::unique_ptr< DPEDSolver > > v_on_eds;
 std::unique_ptr< DPEDSolver > f_start_ed;  ///< ED solver of the source when
                                            ///< it acts as an ON node
 Index ed_pool_th{ 0 };            ///< time horizon the ED pool was built for

 int f_max_thread{ 0 };           ///< intMaxThread: 0 = all cores, 1 = serial
 int f_par_min_n{ TUDPS_PAR_MIN_N }; ///< intParMinN: min horizon for parallel

#if TUDPS_PARALLEL
 /// FastFlow parallel-for engine for compute_EDPs(), one per solver instance
 /// (created on first parallel use); its worker threads are reused across
 /// re-solves, and the per-worker cost scratch lives in f_tcost
 std::unique_ptr< ff::ParallelFor > f_pf;
 std::vector< std::vector< double > > f_tcost;
#endif

 std::vector< double > P;          ///< power values
 std::vector< bool > U;            ///< commitment values

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

};  // end( class( ThermalUnitDPSolver ) )

};  // end( namespace SMSpp_di_unipi_it )

#endif  /* ThermalUnitDPSolver.h included */

/*--------------------------------------------------------------------------*/
/*------------------ End File ThermalUnitDPSolver.h ------------------------*/
/*--------------------------------------------------------------------------*/
