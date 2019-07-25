/*--------------------------------------------------------------------------*/
/*----------------- File BatteryStorageUnitBlock.cpp -----------------------*/
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
#include "BatteryStorageUnitBlock.h"
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

// register BatteryStorageUnitBlock to the Block factory


SMSpp_insert_in_factory_cpp_1( BatteryStorageUnitBlock );

/*--------------------------------------------------------------------------*/
/*------------------- METHODS OF BatteryStorageUnitBlock -------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
void BatteryStorageUnitBlock::deserialize( netCDF::NcGroup & group ) {

 UnitBlock::deserialize( group );

 //TODO

}// end( BatteryStorageUnitBlock::deserialize )

/*--------------------------------------------------------------------------*/

void BatteryStorageUnitBlock::generate_abstract_variables
        ( Configuration *stvv )
{
 UnitBlock::generate_abstract_variables( stvv );



} // end( BatteryStorageUnitBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void BatteryStorageUnitBlock::generate_abstract_constraints
        ( Configuration *stcc )
{


} // end( BatteryStorageUnitBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void BatteryStorageUnitBlock::generate_objective( Configuration *objc )
{
 if( get_objective() == nullptr )  // an objective is there already
  return;                         // cowardly (and silently) return

// initialize objective function - - - - - - - - - - - - - - - - - - - - - -
 //TODO
}  // end( BatteryStorageUnitBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*--- METHODS FOR LOADING, PRINTING & SAVING THE BatteryStorageUnitBlock ---*/
/*--------------------------------------------------------------------------*/

void BatteryStorageUnitBlock::serialize( netCDF::NcGroup & group ) const {

}  // end( BatteryStorageUnitBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------- End File BatteryStorageUnitBlock.cpp -----------------------*/
/*--------------------------------------------------------------------------*/