/*--------------------------------------------------------------------------*/
/*------------- File IntermittentUnitBlock.cpp -------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the IntermittentUnitBlock class.
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
#include "IntermittentUnitBlock.h"
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

// register IntermittentUnitBlock to the Block factory


SMSpp_insert_in_factory_cpp_1( IntermittentUnitBlock );

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS OF IntermittentUnitBlock -------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
void IntermittentUnitBlock::deserialize( netCDF::NcGroup & group ) {

 UnitBlock::deserialize( group );

 //TODO

}// end( IntermittentUnitBlock::deserialize )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::generate_abstract_variables
        ( Configuration *stvv )
{
 UnitBlock::generate_abstract_variables( stvv );



} // end( IntermittentUnitBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::generate_abstract_constraints
        ( Configuration *stcc )
{


} // end( IntermittentUnitBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::generate_objective( Configuration *objc )
{
 if( get_objective() == nullptr )  // an objective is there already
  return;                         // cowardly (and silently) return

// initialize objective function - - - - - - - - - - - - - - - - - - - - - -
 //TODO
}  // end( IntermittentUnitBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*--- METHODS FOR LOADING, PRINTING & SAVING THE IntermittentUnitBlock -----*/
/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::serialize( netCDF::NcGroup & group ) const {

}  // end( IntermittentUnitBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------------- End File IntermittentUnitBlock.cpp -------------------*/
/*--------------------------------------------------------------------------*/