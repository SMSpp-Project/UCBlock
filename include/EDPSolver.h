/*--------------------------------------------------------------------------*/
/*----------------------------- File EDPSolver.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the EDPSolver class.
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
 * \author Kostas Tavlaridis-Gyparakis \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; Claudio Gentile, Antonio Frangioni, Niccolo' Iardella,
 *                   Kostas Tavlaridis-Gyparakis
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __EDPSOLVER_H
 #define __EDPSOLVER_H

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <SMSTypedefs.h>

#include "ThermalUnitDPSolver.h"

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it {

class Block;

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS EDPSolver -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// Class containing the Economic Dispatch Problem Solver.
/** EDPSolver is used in order to calculate the Economic Dispatch Problems
 * arising from a ThermalUnitDPSolver attached to the Thermal Unit with
 * Ramp Constraints (implemented by SMS++ ThermalUnitBlock).
 * Once ThermalUnitDPSolver builds the graph with all the nodes and
 * arc connections,then all the EDPs for the valid arcs are precalculated,
 * to be used in the min-path algorithm.
 * The EDP takes into consideration the ramp constraints of the problem,
 * and provides the cost that will be used in said min-path algorithm.
 */

class EDPSolver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/**
 * @name Constructor and Destructor
 * @{
 */

 EDPSolver() = default;

 virtual ~EDPSolver() = default;

/** @} ---------------------------------------------------------------------*/
/*--------------------- PUBLIC METHODS OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
/**
 * @name Public Methods for building and solving the EDP.
 * @{
 */

 /// initializes the EDP

 void initialize( int k , ThermalUnitDPSolver * s );

 /// compute the cost vector z_h[ h ], z_h[ h + 1 ], ..., z_h[ t - 1]

 void compute_costs( std::vector< double > & costs );

 /// compute optimal power values p_h[ h ], p_h[ h + 1 ], ..., p_h[ k - 1 ]

 void compute_power_variables( int k , std::vector< double > & p );

/** @} ---------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 protected:

 /// the ThermalUnitDPSolver using this EDPSolver
 ThermalUnitDPSolver * solver{};

 /// Coefficients for a variable of the objective function
 struct coeff_t {
  double alfa;
  double beta;
  double gamma;
  } __attribute__((aligned(32)));

 /// Cost coefficients of the objective function
 std::vector< coeff_t > coeffs;

 /// Indices for a piece of the (piece-wise) objective function
 struct pos_t {
  int begt;
  int begm;
  } __attribute__((aligned(8)));

 /**
  * For each k = h, ..., n-1 the vector contains the indices of the pieces
  * of the objective function.
  */
 std::vector< pos_t > pos;

 /// unconstrained optimal power values
 std::vector< double > unc_p;

 /// constrained optimal power values
 std::vector< double > con_p;

 /// Tolerance
 double eps{ 1e-10 };

 int h{};
 int kMax{};

 std::vector< double > m;
 std::vector< int > v;

 };

}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif
