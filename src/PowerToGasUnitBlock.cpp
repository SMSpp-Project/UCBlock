/*--------------------------------------------------------------------------*/
/*------------------- File PowerToGasUnitBlock.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the PowerToGasUnitBlock class.
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
#include "PowerToGasUnitBlock.h"
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

// register PowerToGasUnitBlock to the Block factory


SMSpp_insert_in_factory_cpp_1( PowerToGasUnitBlock );

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS OF PowerToGasUnitBlock ---------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
void PowerToGasUnitBlock::deserialize( netCDF::NcGroup & group ) {

 UnitBlock::deserialize( group );

 //TODO

}// end( PowerToGasUnitBlock::deserialize )

/*--------------------------------------------------------------------------*/

void PowerToGasUnitBlock::generate_abstract_variables
        ( Configuration *stvv )
{
 UnitBlock::generate_abstract_variables( stvv );



} // end( PowerToGasUnitBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void PowerToGasUnitBlock::generate_abstract_constraints
        ( Configuration *stcc )
{


} // end( PowerToGasUnitBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void PowerToGasUnitBlock::generate_objective( Configuration *objc )
{
 if( get_objective() == nullptr )  // an objective is there already
  return;                         // cowardly (and silently) return

// initialize objective function - - - - - - - - - - - - - - - - - - - - - -
 //TODO
}  // end( PowerToGasUnitBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*----- METHODS FOR LOADING, PRINTING & SAVING THE PowerToGasUnitBlock -----*/
/*--------------------------------------------------------------------------*/

void PowerToGasUnitBlock::serialize( netCDF::NcGroup & group ) const {

}  // end( PowerToGasUnitBlock::serialize )

/*--------------------------------------------------------------------------*/
/*----------------- End File PowerToGasUnitBlock.cpp -----------------------*/
/*--------------------------------------------------------------------------*/