/*--------------------------------------------------------------------------*/
/*-------------------- File AcadThemalUnitGraphSolver.h --------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the Solver of Academic ThermalUnit with Ramp Constraints
 *
 *
 *
 * \version 0.10
 *
 * \date 14 - 11 - 2016
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Kostas Tavlaridis-Gyparakis \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Antonio Frangioni, Kostas Tavlaridis-Gyparakis
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __AcadThemalUnitGraphSolver
 #define __AcadThemalUnitGraphSolver
   /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "SMSTypedefs.h"
#include "Solver.h"
#include "EDPSolver.h"

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

 class Block;  ///< forward definition of class Block

/*--------------------------------------------------------------------------*/
/*------------------ CLASS AcadThemalUnitGraphSolver -----------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// derived class for DPSolver the Academic 1UC with ramp Constraints
/** The AcadThemalUnitGraphSolver following a similar procedure with the one
 * described in AcadThemalUnitWRSolver represents a DPSolver for tackling the
 * 1UC Academic Problem with Ramp-Constraints. The difference between the two
 * is the way the graph is constructed in this case the upgraded graph that is
 * of size 2*t (where t is the time-horizon) exceeds the following structure:
 *
 * (0,1) (1,1) (2,1) (3,1) (4,1)  (5,1) . . . (n,1)
 * (0,0) (1,0) (2,0) (3,0) (4,0)  (5,0) . . . (n,0)
 *
 * Where each arc transition from unit being on to unit being of denotes an 
 * EDP for all the corresponding time steps e.g. (1,1) -> (5,0) is EDP (1,4),
 * where the unit has been producing energy for time-steps 1-4. And in a similar
 * manner the arc transition from unit being off to on denotes the the start-up
 * costs of the unit e.g. (1,0) -> (5,1) means we have to calculate the start-up 
 * costs for the unit being off for 4 time steps.
 * And so one has to compute all the EDP for the valid arc transitions based on 
 * the constraints and also all the different start-up costs and then solve a 
 * min-path problem for the graph.
 *
 * The Class is caring out this procedure by constructing the graph vie the use
 * of the method build_graph() and storing the graph in the vector of structs
 * denoted as v_nodes. 
 * Then all the different EDP's for the valid arc connections are calculated
 * and then a min-path algorithm is used (via the method MinPath()) in order
 * to obtain the best production schedule for the unit. 
*/

class AcadThemalUnitGraphSolver : public Solver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

public:
 AcadThemalUnitGraphSolver( void ) : Solver() {} ;

 virtual ~AcadThemalUnitGraphSolver() {}; 

/*@}------------------------------------------------------------------------*/
/*--------------------- DERIVED METHODS OF BASE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

/** @name Public Methods derived of the Base Class
 *  @{ */

    void set_Block( Block *block );
    ///< method for setting the Block and initializing the DP
    int solve (); 
    ///< Method for solving the constructed DP

    bool new_var_solution( void ) {return (false);}

	void get_var_solution( void );
    ///< Method for retrieving writing the solution in Block

/*@}------------------------------------------------------------------------*/
/*---------------------- PUBLIC METHODS OF THE CLASS  ----------------------*/
/*--------------------------------------------------------------------------*/

/** @name Public Methods for constructing and tackling the DP
 *  @{ */

  void build_graph ();
  ///< method for building the DP-Graph
  void MinPath();
  ///< method for implementing the min-path algorithm

  double ComputeStartupCosts( int t);
  ///< method for computing the start-up costs

  
  double get_fTotalCost() {return fTotalCost;}
  ///< method for returing the total cost of the solution


/*@}------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Protected Fields of the AcadThemalUnitWRSolver
 *  @{ */


protected:

  struct arc{
      arc()
          :h(-32000)
          ,k(-32000)
          ,cost1(-32000)
          ,cost2(-32000)
          ,valid(-32000)
      {}
    int h;
    int k;
    double cost1;
    double cost2;
    int valid;
  };


    struct node {
        node()
            :v_arcs()
    {}
        std::vector< arc > v_arcs;
    };

    std::vector< node > v_nodes;
    /**<
      The vector that denotes all the different nodes of the graph, having the 
      following order:
      [(0,1), (1,1), (2,1), ... , (t,1), (0,0), (1,0), (2,0), ... , (t,0)]
      Each of these nodes store a vector of arcs (v_arcs) that hold all the valid arc
      connections for the corresponding node, where in the case of node where the unit
      is active the arcs are connecting this node with the possible time steps where
      it will be shut down and thus an EDP needs to be calculated and for the case of
      nodes where the unit is off, the arcs are connecting this node with the possible
      time steps where the unit will be turned on and thus the corresponding start up
      costs need to be calculated.
     */

  struct t_route {
        t_route()
         :h(-32000)
         ,lab(-32000)
         ,pred(-32000)
    {};
    
        int h;
        double lab;
        int pred;
  };
    
  std::vector< t_route > v_route;
  ///< vector that stores the obtained from the min-path algorithm route


  std::vector< EDPSolver > v_EDP;
  ///< vector for the EDPSolvers of all the corresponding valid arcs

  std::vector< double > startup_costs; 
  std::vector< double > v_cost;
  ///< vectors storing the costs

  std::vector< double > fPowerVector;
  std::vector< int > fStatusVector;
  ///< vectors for storing the commitement and power variable values
    

  int    hmin ;
  ///< variable denoting firt time period the unit can be turned on
  int    kMin ; 
  ///< variable denoting firt time period the unit can be turned off
  double fTotalCost;
  ///< variable denoting the total costs for the unit
  int    fMaxStartupLevel;
  /**< variable denoting the timesteps after which unit being offline
       has stable costs
    */


/*@}*/


};   // end( class AcadThemalUnitGraphSolver )

}; // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* AcadThemalUnitGraphSolver.h included */
