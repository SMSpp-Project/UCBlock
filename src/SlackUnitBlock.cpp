/*--------------------------------------------------------------------------*/
/*------------------------ File SlackUnitBlock.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the SlackUnitBlock class.
 *
 * \version 0.11
 *
 * \date 10 - 12 - 2019
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
#include "SlackUnitBlock.h"
#include "LinearFunction.h"
#include <map>
#include "UnitBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;


/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register SlackUnitBlock to the Block factory


SMSpp_insert_in_factory_cpp_1( SlackUnitBlock );

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS OF SlackUnitBlock ------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
void SlackUnitBlock::deserialize( netCDF::NcGroup & group ) {

 UnitBlock::deserialize( group );


}// end( SlackUnitBlock::deserialize )

/*--------------------------------------------------------------------------*/

void SlackUnitBlock::generate_abstract_variables
        ( Configuration *stvv )
{
 UnitBlock::generate_abstract_variables( stvv );
 if( f_time_horizon == 0 ) {
  // there are no variables to be generated
  return;
 }

 if( !v_commitment.empty()  ) {
  // the abstract variables should be generated only once
  return;
 }
 v_commitment.resize(boost::extents[ f_time_horizon ][ get_number_generators() ]);

 for( Index t = 0; t < f_time_horizon; ++t ) {
  for( Index g = 0; g < get_number_generators(); ++g ) {
   auto & Commitment = v_commitment[ t ][ g ];

   Commitment.set_type( ColVariable::kPosUnitary );

   add_static_variable ( Commitment, "u_" +
                                     std::to_string( t ) + "_" +
                                     std::to_string( g )  );
  }
 }
} // end( SlackUnitBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void SlackUnitBlock::generate_abstract_constraints
        ( Configuration *stcc )
{

} // end( SlackUnitBlock::generate_abstract_constraints )


/*--------------------------------------------------------------------------*/

void SlackUnitBlock::generate_objective( Configuration *objc )
{
}  // end( SlackUnitBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*------- METHODS FOR LOADING, PRINTING & SAVING THE SlackUnitBlock --------*/
/*--------------------------------------------------------------------------*/

void SlackUnitBlock::serialize( netCDF::NcGroup & group ) const {

 UnitBlock::serialize( group );

}  // end( SlackUnitBlock::serialize )

/*--------------------------------------------------------------------------*/
/*----------------------- End File SlackUnitBlock.cpp ----------------------*/
/*--------------------------------------------------------------------------*/