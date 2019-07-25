/*--------------------------------------------------------------------------*/
/*-------------------- File EMobilityUnitBlock.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the EMobilityUnitBlock class.
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
#include "EMobilityUnitBlock.h"
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

// register EMobilityUnitBlock to the Block factory


SMSpp_insert_in_factory_cpp_1( EMobilityUnitBlock );

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS OF EMobilityUnitBlock ----------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
void EMobilityUnitBlock::deserialize( netCDF::NcGroup & group ) {

 UnitBlock::deserialize( group );

 //TODO

}// end( EMobilityUnitBlock::deserialize )

/*--------------------------------------------------------------------------*/

void EMobilityUnitBlock::generate_abstract_variables
        ( Configuration *stvv )
{
 UnitBlock::generate_abstract_variables( stvv );



} // end( EMobilityUnitBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void EMobilityUnitBlock::generate_abstract_constraints
        ( Configuration *stcc )
{


} // end( EMobilityUnitBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void EMobilityUnitBlock::generate_objective( Configuration *objc )
{
 if( get_objective() == nullptr )  // an objective is there already
  return;                         // cowardly (and silently) return

// initialize objective function - - - - - - - - - - - - - - - - - - - - - -
 //TODO
}  // end( EMobilityUnitBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*----- METHODS FOR LOADING, PRINTING & SAVING THE EMobilityUnitBlock ------*/
/*--------------------------------------------------------------------------*/

void EMobilityUnitBlock::serialize( netCDF::NcGroup & group ) const {

}  // end( EMobilityUnitBlock::serialize )

/*--------------------------------------------------------------------------*/
/*---------------- End File EMobilityUnitBlock.cpp -------------------------*/
/*--------------------------------------------------------------------------*/