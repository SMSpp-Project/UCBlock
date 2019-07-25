/*--------------------------------------------------------------------------*/
/*------------- File CentralizedDemandResponseUnitBlock.cpp ----------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the CentralizedDemandResponseUnitBlock class.
 *
 * \version 0.11
 *
 * \date 25 - 07 - 2019
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 *
 * \author Ali Ghezelsoflu \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 *
 * Copyright &copy by Antonio Frangioni, Ali Ghezelsoflu
 */

/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <iostream>
#include <random>
#include "CentralizedDemandResponseUnitBlock.h"
#include "LinearFunction.h"
#include <map>
#include "FRowConstraint.h"
#include "DQuadFunction.h"
#include "FRealObjective.h"
#include "UnitBlock.h"
#include "UCBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;


/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register CentralizedDemandResponseUnitBlock to the Block factory


SMSpp_insert_in_factory_cpp_1( CentralizedDemandResponseUnitBlock );

/*--------------------------------------------------------------------------*/
/*------------- METHODS OF CentralizedDemandResponseUnitBlock --------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
void CentralizedDemandResponseUnitBlock::deserialize( netCDF::NcGroup & group ) {

 UnitBlock::deserialize( group );

 //TODO

}// end( CentralizedDemandResponseUnitBlock::deserialize )

/*--------------------------------------------------------------------------*/

void CentralizedDemandResponseUnitBlock::generate_abstract_variables
        ( Configuration *stvv )
{
 UnitBlock::generate_abstract_variables( stvv );



} // end( CentralizedDemandResponseUnitBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void CentralizedDemandResponseUnitBlock::generate_abstract_constraints
        ( Configuration *stcc )
{


} // end( CentralizedDemandResponseUnitBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void CentralizedDemandResponseUnitBlock::generate_objective( Configuration *objc )
{
 if( get_objective() == nullptr )  // an objective is there already
  return;                         // cowardly (and silently) return

// initialize objective function - - - - - - - - - - - - - - - - - - - - - -
 //TODO
}  // end( CentralizedDemandResponseUnitBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*-- METHODS FOR LOADING, PRINTING & SAVING THE CentralizedDemandResponse --*/
/*--------------------------------------------------------------------------*/

void CentralizedDemandResponseUnitBlock::serialize( netCDF::NcGroup & group ) const {

}  // end( CentralizedDemandResponseUnitBlock::serialize )

/*--------------------------------------------------------------------------*/
/*--------- End File CentralizedDemandResponseUnitBlock.cpp ----------------*/
/*--------------------------------------------------------------------------*/