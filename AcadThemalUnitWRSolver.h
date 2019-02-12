/*--------------------------------------------------------------------------*/
/*---------------------- File AcadThemalUnitWRSolver.h ---------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the 1UC DPSolver of Academic ThermalUnit with Ramp Constraints
 *
 *
 *
 * \version 0.10
 *
 * \date 24 - 09 - 2016
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

#ifndef __AcadThemalUnitWRSolver
 #define __AcadThemalUnitWRSolver  
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
/*------------------- CLASS AcadThemalUnitWRSolver -------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// derived class for DPSolver the Academic 1UC with ramp Constraints
/** The AcadThemalUnitWRSolver represents the DPSolver for tackling the 1UC
 * Academic Problem with Ramp-Constraints via the following DP procedure:
 * A graph is constructed that is of t^2 size (where t denotes time-horizon),
 * where (t,k) represents that the unit in time step t can be in state k, which
 * means that the weather is online (k=1), or offline for a series of time step
 * (k=0,-1,-2,...,-t). The arcs that are connecting the different nodes of the
 * graph wether represent a start up cost if the unit goes from offline to on-
 * line, or an EDP if the unit goes from online to offline state. So for exam-
 * ple the arc (3,1)->(7,0) represents the EDP (3,6) where unit has been produ-
 * cing energy for timesteps 3-6. And in a similar fashion the arc (3,-4)->(4,1)
 * represents the start up cost to start the unit in time step 4 after being
 * switched off for 4 time-steps.
 *
 * The Class is caring out this procedure by constructing the graph vie the use
 * of the method build_graph() and storing the graph in the vector of structs
 * denoted as v_nodes. 
 * Then all the different EDP's for the valid arc connections are calculated
 * and then a min-path algorithm is used (via the method MinPath()) in order
 * to obtain the best production schedule for the unit. 
*/

class AcadThemalUnitWRSolver : public Solver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */
public:
 AcadThemalUnitWRSolver( void ) : Solver() {} ;

 virtual ~AcadThemalUnitWRSolver() {}; 
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

  int position(int h, int k);
  ///< method for tracking the position of nodes in graph
  
  double get_fTotalCost() {return fTotalCost;}
  ///< method for returing the total cost of the solution


/*@}------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Protected Fields of the AcadThemalUnitWRSolver
 *  @{ */

protected:


  struct t_node{
        t_node()
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
  
  std::vector< t_node > v_nodes;
  ///< vector that stores all the nodes of the problem each on including:
     
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
};   // end( class AcadThermalUnitSolver )

}; // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* AcadThemalUnitWRSolver.h included */
