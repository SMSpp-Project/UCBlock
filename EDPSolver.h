/*--------------------------------------------------------------------------*/
/*----------------------------- File EDPSolver.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the Economic Dispatch Solver of Academic ThermalUnit with 
 * Ramp Constraints
 *
 *
 * \version 0.10
 *
 * \date 03 - 10 - 2016
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

#ifndef __EDPSolver
 #define __EDPSolver  
   /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "SMSTypedefs.h"
#include "UCBlock.h"
#include "AcadThermalUnitBlock.h"

/*--------------------------------------------------------------------------*/
/*----------------------------- NAMESPACE ----------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

 class Block;  ///< forward definition of class Block

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS EDPSolver -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// class containing the EDPSolver
/**
 * The following class is used in order to calculate the EDPs arising from 
 * the DPSolvers attached to the Academic ThermalUnit with Ramp Constraints.
 * For each DP once the Graph with all the nodes and arc connections is cons-
 * tructed then all the EDPs for the valid arcs are precalculated prior to
 * the min-path algorithm. 
 * The EDP is taking under consideration the ramp constraints of the problem
 * and provides the cost that will be used in the min-path algorithm that 
 * will be followed.
 */


class EDPSolver {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

public:
 EDPSolver()  {} ;

 virtual ~EDPSolver() {}; 

/*@}------------------------------------------------------------------------*/
/*---------------------- PUBLIC METHODS OF THE CLASS  ----------------------*/
/*--------------------------------------------------------------------------*/

/** @name Public Methods for constructing and tackling the DP
 *  @{ */

    /* method for initializing the fields of the EDPSolver. */ 
    void initialize(int k, UCBlock *f_UCBlock);

    /* method for the calculation of the costs z_h[h],z_h[h+1],..., z_h[k-1]. */ 
    void ComputeCosts(std::vector< double > &costVector, AcadThermalUnitBlock *unit_block, UCBlock *f_UCBlock);

    void ComputePowerVariables(int k, vector<double> &powerVector, AcadThermalUnitBlock *unit_block);
    /* method for calculating the potential optimal p_h[h], p_h[h+1],...,
    p_h[k-1]. Setta quindi i campi powerVector delle potenze comprese
    nel nodo (h,k) */

/*@}------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Protected Fields of the AcadThemalUnitWRSolver
 *  @{ */

protected:

    struct tCoeff{
        double alfa;
        double beta;
        double gamma;
    };

     std::vector<tCoeff> Coeff;
    ///< vector storing the cost coefficients of the objective function

    struct tPosizioni{
        int begt;
        int begm;
    };

    std::vector<tPosizioni> Position;
    /**< Structure containing for each k = h, ..., n-1 the initial indices of
         various sections of which is composed of the objective cost function */


    std::vector< double > unconstrPowerOpt;
    /** vector storing the unconstrainted power optimum values  */

    std::vector< double > constrPowerOpt;
    /** vector storing the unconstrainted power optimum values  */

    double eps ;
    ///< tolerance

    int h;
    int kMax;

    std::vector< double > m;
    std::vector<int> v;

};   // end( class EDPSolver )

}; // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* AcadThemalUnitWRSolver.h included */

