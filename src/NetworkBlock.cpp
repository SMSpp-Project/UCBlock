/*--------------------------------------------------------------------------*/
/*--------------------- File NetworkBlock.cpp ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the NetworkBlock class.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Ali Ghezelsoflu \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Rafael Durbano Lobato \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Kostas Tavlaridis-Gyparakis \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy by Antonio Frangioni, Ali Ghezelsoflu,
 *                  Rafael Durbano Lobato
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <map>

#include "LinearFunction.h"
#include "NetworkBlock.h"
#include "RowConstraintSolution.h"
#include "ColRowSolution.h"
#include "ColVariableSolution.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

NetworkBlock::NetworkData::NetworkData( void ) {
 f_number_lines = 0;
 f_number_nodes = 0;
 f_number_intervals = 0;
}

/*--------------------------------------------------------------------------*/

void NetworkBlock::NetworkData::deserialize( const netCDF::NcGroup & group ) {

 // Optional variables

 if( ! ::deserialize_dim( group , "NumberNodes" ,
                         f_number_nodes , true ) )
  f_number_nodes = 1;

 if( ! ::deserialize_dim( group , "NumberIntervals" ,
                         f_number_intervals , true ) )
  f_number_intervals = 1;
}

/*--------------------------------------------------------------------------*/

void NetworkBlock::deserialize( const netCDF::NcGroup & group ) {

 Block::deserialize( group );

 // Optional variables

 if( ! ::deserialize( group , f_const_term , "ConstTerm" ) )
  f_const_term = 0;
}

/*--------------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/

Solution * NetworkBlock::get_Solution( Configuration * csolc , bool emptys ) {
 Index solution_type = 0;
 if( ( ! csolc ) && f_BlockConfig )
  csolc = f_BlockConfig->f_solution_Configuration;
 if( auto config = dynamic_cast< SimpleConfiguration< int > * >( csolc ) )
  solution_type = config->f_value;

 Solution * sol;
 switch( solution_type ) {
  case 1:
   sol = new RowConstraintSolution;
   break;
  case 2:
   sol = new ColRowSolution;
   break;
  default:
   sol = new ColVariableSolution;
 }

 if( ! emptys )
  sol->read( this );

 return( sol );
}

/*--------------------------------------------------------------------------*/
/*---------- METHODS FOR LOADING, PRINTING & SAVING THE NetworkBlock -------*/
/*--------------------------------------------------------------------------*/

void NetworkBlock::NetworkData::serialize( netCDF::NcGroup & group ) const {
 group.addDim( "NumberNodes" , f_number_nodes );
 group.addDim( "NumberIntervals" , f_number_intervals );

 if( f_number_nodes > 1 ) {
  auto NumberLines = group.addDim( "NumberLines" , f_number_lines );

  ::serialize( group , "StartLine" , netCDF::NcUint() , NumberLines ,
               v_start_line );

  ::serialize( group , "EndLine" , netCDF::NcUint() , NumberLines ,
               v_end_line );
 }
}

/*--------------------------------------------------------------------------*/

NetworkBlock::NetworkData::NetworkDataFactoryMap &
NetworkBlock::NetworkData::f_factory( void ) {
 static NetworkDataFactoryMap s_factory;
 return( s_factory );
}

/*--------------------------------------------------------------------------*/
/*--------------------- End File NetworkBlock.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
