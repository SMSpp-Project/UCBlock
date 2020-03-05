/*--------------------------------------------------------------------------*/
/*------------------------- File UnitBlock.cpp -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the UnitBlock class.
 *
 * \version 0.11
 *
 * \date 08 - 07 - 2019
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Ali Ghezelsoflu \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Rafael Durbano Lobato \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Kostas Tavlaridis-Gyparakis \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy by Antonio Frangioni, Ali Ghezelsoflu, Rafael
 * Durbano Lobato, and Kostas Tavlaridis-Gyparakis
 */

/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "UCBlock.h"
#include "UnitBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register UnitBlock to the Block factory

SMSpp_insert_in_factory_cpp_1( UnitBlock );

/*--------------------------------------------------------------------------*/
/*--------------------------- METHODS OF UnitBlock -------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
UnitBlock::UnitBlock( Block * father_block, UnitBlock::Index t )
 : Block( father_block ), f_time_horizon( t ) {
 f_number_intervals = 0;
}

void UnitBlock::deserialize_time_horizon( netCDF::NcGroup & group ) {
 netCDF::NcDim TimeHorizon = group.getDim( "TimeHorizon" );
 if( TimeHorizon.isNull() ) {
  // dimension TimeHorizon is not present in the netCDF input

  if( f_time_horizon == 0 ) {
   auto f_B = dynamic_cast< UCBlock * >( get_f_Block() );
   if( f_B )
    // The father Block is available. Take time horizon from it.
    this->set_time_horizon( f_B->get_time_horizon() );
   else
    throw ( std::invalid_argument(
     "UnitBlock::deserialize: TimeHorizon is not present in the "
     "netCDF input and UnitBlock does not have a father." ) );
  }
 } else {
  // dimension TimeHorizon is present in the netCDF input

  auto th = TimeHorizon.getSize();
  if( f_time_horizon == 0 )
   this->set_time_horizon( th );
  else if( f_time_horizon != th )
   throw ( std::logic_error(
    "UnitBlock::deserialize: TimeHorizon is not present in the "
    "netCDF. The (nonzero) time horizon of UnitBlock is different "
    "from that of its father, but they should be equal." ) );
 }
}

/*--------------------------------------------------------------------------*/

void UnitBlock::deserialize_change_intervals( netCDF::NcGroup & group ) {

 auto NumberIntervals = group.getDim( "NumberIntervals" );
 if( NumberIntervals.isNull() )
  f_number_intervals = 1;
 else {
  f_number_intervals = NumberIntervals.getSize();
  if( ( f_number_intervals < 1 ) || ( f_number_intervals > f_time_horizon ) )
   throw ( std::invalid_argument
    ( "UnitBlock::deserialize: invalid NumberIntervals. "
      "It must be between 1 and TimeHorizon." ) );
 }

 if( ( f_number_intervals > 1 ) && ( f_number_intervals < f_time_horizon ) ) {

  ::deserialize( group, "ChangeIntervals", f_number_intervals,
                 v_change_intervals );

  // Check that all numbers are between 1 and f_time_horizon, that
  // the last number is == f_time_horizon, and that they are ordered
  // in increasing sense

  if( v_change_intervals.back() != f_time_horizon ) {
   throw ( std::invalid_argument
    ( "UnitBlock::deserialize: invalid value in ChangeIntervals: "
      "the last element must be TimeHorizon." ) );
  }

  Index previous_t = 0;

  for( auto t : v_change_intervals ) {
   if( !( t > previous_t && t < f_time_horizon - 1 ) )
    throw ( std::invalid_argument
     ( "UnitBlock::deserialize: invalid value in ChangeIntervals: " +
       std::to_string( t ) + ". All values must be between 1 and "
                             "TimeHorizon and in strictly increasing order." ) );

   previous_t = t;
  }
 }
}

/*--------------------------------------------------------------------------*/

void UnitBlock::deserialize( netCDF::NcGroup & group ) {
 // deserialize_time_horizon( group );
 // deserialize_change_intervals( group );

 Block::deserialize( group );
}

/*--------------------------------------------------------------------------*/
/*------------------ METHODS FOR MODIFYING THE UnitBlock -------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR SAVING THE UnitBlock -------------------*/
/*--------------------------------------------------------------------------*/

void UnitBlock::serialize( netCDF::NcGroup & group ) const {

 Block::serialize( group );

 group.addDim( "TimeHorizon", f_time_horizon );

 auto NumberIntervals = group.addDim( "NumberIntervals", f_number_intervals );

 if( !v_change_intervals.empty() ) {
  ::serialize( group, "ChangeInterval", netCDF::NcUint64(),
               NumberIntervals, v_change_intervals );
 }
}

/*--------------------------------------------------------------------------*/
/*---------------------- End File UnitBlock.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
