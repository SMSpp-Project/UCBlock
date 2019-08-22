/*--------------------------------------------------------------------------*/
/*----------------- File BatteryUnitBlock.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the BatteryStorageUnitBlock class.
 *
 * \version 0.11
 *
 * \date 18 - 07 - 2019
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
#include "BatteryUnitBlock.h"
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

// register BatteryUnitBlock to the Block factory


SMSpp_insert_in_factory_cpp_1( BatteryUnitBlock );

/*--------------------------------------------------------------------------*/
/*------------------- METHODS OF BatteryUnitBlock -------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
void BatteryUnitBlock::deserialize( netCDF::NcGroup & group ) {

 UnitBlock::deserialize( group );

 //TODO

}// end( BatteryUnitBlock::deserialize )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::generate_abstract_variables
        ( Configuration *stvv )
{
 UnitBlock::generate_abstract_variables( stvv );



} // end( BatteryUnitBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::generate_abstract_constraints
        ( Configuration *stcc )
{


} // end( BatteryUnitBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::generate_objective( Configuration *objc )
{
 if( get_objective() == nullptr )  // an objective is there already
  return;                         // cowardly (and silently) return

// initialize objective function - - - - - - - - - - - - - - - - - - - - - -
 //TODO
}  // end( BatteryUnitBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*--- METHODS FOR LOADING, PRINTING & SAVING THE BatteryUnitBlock ---*/
/*--------------------------------------------------------------------------*/

void BatteryUnitBlock::serialize( netCDF::NcGroup & group ) const {

}  // end( BatteryUnitBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------- End File BatteryUnitBlock.cpp -----------------------*/
/*--------------------------------------------------------------------------*/