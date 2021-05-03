/*--------------------------------------------------------------------------*/
/*----------------------- File ThermalUnitDPSolver.h -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the ThermalUnitDPSolver class.
 *
 * ThermalUnitDPSolver is the DPSolver for the Thermal Unit with
 * Ramp Constraints.
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Niccolò Iardella \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Kostas Tavlaridis-Gyparakis \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; Antonio Frangioni, Niccolò Iardella, Kostas Tavlaridis-Gyparakis
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __THERMALUNITDPSOLVER_H
#define __THERMALUNITDPSOLVER_H

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <SMSTypedefs.h>
#include <Solver.h>
#include "EDPSolver.h"

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it {

class Block;

class EDPSolver;

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS ThermalUnitDPSolver ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/

/// Class for solving a Single Unit Commitment problem with a DP approach.
/**
 * The ThermalUnitDPSolver is a DPSolver for tackling the
 * 1UC Academic Problem with Ramp Constraints. The solver takes advantage
 * of an improved graph with the following structure:
 *
 * (0,1) (1,1) (2,1) (3,1) (4,1)  (5,1) . . . (n,1)
 * (0,0) (1,0) (2,0) (3,0) (4,0)  (5,0) . . . (n,0)
 *
 * Each arc transition from an ON node (i, 1) to an OFF node (j, 0)
 * represents an Economic Dispatch Problem for the involved time steps.
 * For example, (1,1) -> (5,0) is an EDP (1,4) where the unit has been
 * producing energy for the time-steps 1-4.
 * In a similar manner, the arc transition from an OFF node (i, 0) to
 * an ON node (j, 1) is represented by the start-up costs for the unit.
 * For example, (1,0) -> (5,1) means that start-up costs for the unit
 * being off for 4 time steps must be calculated.
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
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

/**
 * @name Constructor and destructor
 * @{
 */

 ThermalUnitDPSolver() : Solver() {};

 ~ThermalUnitDPSolver() override = default;
 /// @}

/*--------------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

/**
 * @name Public methods derived from base classes
 * @{
 */

 /// It sets the Block that the Solver has to solve
 void set_Block( Block * block ) override;

 /// It solves the constructed problem
 int compute( bool changedvars = true ) override;

 /// Writes the current solution in the Block
 void get_var_solution( Configuration * solc ) override;

 /// Returns a valid lower bound on the optimal objective function value
 OFValue get_lb() override;

 /// Returns a valid upper bound on the optimal objective function value
 OFValue get_ub() override;

 /// Returns the value of the current solution, if any
 OFValue get_var_value() override;
 /// @}

/*--------------------------------------------------------------------------*/
/*------------------- PROTECTED FIELDS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/**
 * @name Methods for building and tackling the DP
 * @{
 */

 /// Builds the graph
 void build_graph();

 /// Computes the EDPs
 void compute_EDPs();

 /// Implements the min-path algorithm
 void min_path();

 /// Computes the variable values and the total cost
 void compute_solutions();
 /// @}

/*--------------------------------------------------------------------------*/
/*-------------------- PRIVATE FIELDS OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/**
 * @name Thermal unit physical parameters
 * @{
 */
 friend class EDPSolver;

 int time_horizon{};      /// Time horizon
 int init_up_down_time{}; /// Initial up/down time
 int min_up_time{};       /// Minimum up time
 int min_down_time{};     /// Minimum down time
 double initial_power{};  /// Initial power
 int init_t{};

 std::vector< double > startup_costs;
 std::vector< double > delta_ramp_up;
 std::vector< double > delta_ramp_down;
 std::vector< double > min_power;
 std::vector< double > max_power;

 /* These two are identical to min_power but we keep them for readability */
 std::vector< double > & bound_on = min_power;
 std::vector< double > & bound_down = min_power;

 std::vector< double > quad_term;
 std::vector< double > linear_term;
 std::vector< double > const_term;

 /// Loads the parameters from the ThermalUnitBlock
 void load_parameters();
 /// @}

 /// Computes the start-up costs
 double compute_startup_costs( int t );

 /// Tolerance
 double eps{ 1e-10 };

 /// An arc
 struct arc {
  int h{};
  int k{};
  double cost1{};
  double cost2{};
  int valid{};
 } __attribute__((aligned(32)));

 /// A node, i.e., a vector of arcs
 struct node {
  std::vector< arc > v_arcs{};
 } __attribute__((aligned(32)));

 /// Vector of the nodes
 /**
  * The nodes are in the following order:
  * [ (0,1), (1,1), (2,1), ... , (t,1),
  *   (0,0), (1,0), (2,0), ... , (t,0) ]
  *
  * The first t nodes (i, 1) are the ON nodes, the others (i, 0) are the
  * OFF nodes. Each node contains a vector that holds all the valid arc
  * connections to other nodes.
  * ON nodes are connected to the OFF nodes of those time steps that can
  * be reached without violating the constraints. For these connections,
  * EDPs will be calculated.
  * Similarly, OFF nodes are connected to the ON nodes of those time steps
  * that can be reached without violating the constraints.
  * For these connections, start-up costs will be calculated.
  */
 std::vector< node > v_nodes;

 /// A piece of a route
 struct route_t {
  int h{};
  double lab{};
  int pred{};
 } __attribute__((aligned(16)));

 /// Vector that stores the route obtained from the min-path algorithm
 std::vector< route_t > v_route;

 /// EDPSolvers for the arcs that need them
 std::vector< EDPSolver > v_EDP;

 std::vector< double > P;    ///< Power values
 std::vector< int > U;       ///< Commitment values
 std::vector< int > startup; ///< Startup values
 double total_cost{};        ///< Total cost

 int hMin{}; ///< First time step the unit can be turned ON
 int kMin{}; ///< First time step the unit can be turned OFF

/*--------------------------------------------------------------------------*/
/*----------------- INTERFACE FOR SUPPORTING MODIFICATIONS ---------------- */
/*--------------------------------------------------------------------------*/

 /// Stage of the computation
 enum {
  start = 0,
  graph_OK = 1,
  edps_OK = 2,
  path_OK = 3,
  sol_OK = 4
 } stage{};

 /// It processes all the pending modifications
 void process_modifications();

 /*--------------------------------------------------------------------------*/

 void retrieve_term( std::vector< double > & out,
                     const std::vector< double > & in ) const;

 SMSpp_insert_in_factory_h;
};

}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif
