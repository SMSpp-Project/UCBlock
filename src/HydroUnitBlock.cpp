/*--------------------------------------------------------------------------*/
/*--------------------- File HydroUnitBlock.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the HydroUnitBlock class.
 *
 * \version 0.11
 *
 * \date 11 - 07 - 2019
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
#include "HydroUnitBlock.h"
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

// register HydroUnitBlock to the Block factory


SMSpp_insert_in_factory_cpp_1( HydroUnitBlock );

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS OF HydroUnitBlock -----------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
void HydroUnitBlock::deserialize( netCDF::NcGroup & group ) {

 UnitBlock::deserialize( group );

 //TODO

}// end( HydroUnitBlock::deserialize )

/*--------------------------------------------------------------------------*/

void HydroUnitBlock::generate_abstract_variables( Configuration *stvv )
{
 UnitBlock::generate_abstract_variables( stvv );

 if( f_time_horizon == 0 ) {
  // there are no variables to be generated
  return;
 }

 if( !v_volumetric.empty() ||
     !v_flow_rate.empty() ) {
  // the abstract variables should be generated only once
  return;
 }

  v_volumetric.resize(boost::extents[f_time_horizon][f_number_arcs]);
  v_flow_rate.resize(boost::extents[f_time_horizon][f_number_arcs]);

 for( Index t = 0; t < f_time_horizon; ++t ) {
  for( Index g = 0; g < f_number_arcs; ++g ) {
   v_volumetric[ t ][ g ].set_type( ColVariable::kNonNegative );
   v_flow_rate[ t ][ g ].set_type( ColVariable::kContinuous );

  }
 }
 add_static_variable ( v_volumetric );
 add_static_variable ( v_flow_rate );


} // end( HydroUnitBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void HydroUnitBlock::generate_abstract_constraints( Configuration *stcc )
{


} // end( HydroUnitBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void HydroUnitBlock::generate_objective( Configuration *objc )
{
    if( get_objective() == nullptr )  // an objective is there already
        return;                         // cowardly (and silently) return

// initialize objective function - - - - - - - - - - - - - - - - - - - - - -
    //TODO
}  // end( HydroUnitBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*-------- METHODS FOR LOADING, PRINTING & SAVING THE HydroUnitBlock -------*/
/*--------------------------------------------------------------------------*/

void HydroUnitBlock::serialize( netCDF::NcGroup & group ) const {

}  // end( HydroUnitBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------------- End File HydroUnitBlock.cpp --------------------------*/
/*--------------------------------------------------------------------------*/