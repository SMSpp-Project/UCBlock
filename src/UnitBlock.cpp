/*--------------------------------------------------------------------------*/
/*------------------------- File UnitBlock.cpp -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the UnitBlock class.
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

#include "UCBlock.h"
#include "UnitBlock.h"
#include "RowConstraintSolution.h"
#include "ColRowSolution.h"
#include "ColVariableSolution.h"

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

UnitBlock::UnitBlock( Block * father_block, Index t )
 : Block( father_block ) , f_time_horizon( t ) , f_number_intervals( 0 ) {}

/*--------------------------------------------------------------------------*/

void UnitBlock::deserialize_time_horizon( const netCDF::NcGroup & group )
{
 netCDF::NcDim TimeHorizon = group.getDim( "TimeHorizon" );
 if( TimeHorizon.isNull() ) {
  // dimension TimeHorizon is not present in the netCDF input

  if( f_time_horizon == 0 ) {
   if( auto f_B = dynamic_cast< UCBlock * >( get_f_Block() ) )
    // The father Block is available. Take time horizon from it.
    this->set_time_horizon( f_B->get_time_horizon() );
   else if( auto f_B = dynamic_cast< UnitBlock * >( get_f_Block() ) )
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

void UnitBlock::deserialize_change_intervals( const netCDF::NcGroup & group )
{
 auto NumberIntervals = group.getDim( "NumberIntervals" );
 if( NumberIntervals.isNull() )
  f_number_intervals = 1;
 else {
  f_number_intervals = NumberIntervals.getSize();
  if( ( f_number_intervals < 1 ) || ( f_number_intervals > f_time_horizon ) )
   throw( std::invalid_argument( "UnitBlock::deserialize: NumberIntervals "
				 "not between 1 and TimeHorizon." ) );
  }

 if( ( f_number_intervals > 1 ) && ( f_number_intervals < f_time_horizon ) ) {
  ::deserialize( group, "ChangeIntervals", f_number_intervals,
                 v_change_intervals );

  // Check that the numbers are ordered in increasing sense. Notice that the
  // upper endpoint of the last interval must necessarily be f_time_horizon -
  // 1. Since it is not required to be provided, we set it.
  v_change_intervals.back() = f_time_horizon - 1;

  for( Index k = 0 ; k < v_change_intervals.size() ; ++k ) {
   const auto t = v_change_intervals[ k ];
   if( ! ( ( t < f_time_horizon ) &&
           ( k == 0 || t > v_change_intervals[ k - 1 ] ) ) )
    throw ( std::invalid_argument( "UnitBlock::deserialize: invalid value in ChangeIntervals: " +
              std::to_string( t ) + ". All values must be between 0 and "
              "TimeHorizon - 1 and in strictly increasing order." ) );
   }
  }
 else
  v_change_intervals.clear();
 }

/*--------------------------------------------------------------------------*/

void UnitBlock::deserialize( const netCDF::NcGroup & group )
{
 Block::deserialize( group );
 deserialize_time_horizon( group );
 deserialize_change_intervals( group );
 }

/*--------------------------------------------------------------------------*/
/*------------------ METHODS FOR MODIFYING THE UnitBlock -------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/

Solution * UnitBlock::get_Solution( Configuration * csolc , bool emptys )
{
 Index solution_type = 0;
 if( ( ! csolc ) && f_BlockConfig )
  csolc = f_BlockConfig->f_solution_Configuration;

 if( auto config = dynamic_cast< SimpleConfiguration< int > * >( csolc ) )
  solution_type = config->f_value;

 Solution * sol = nullptr;
 switch( solution_type ) {
  case 1:  sol = new RowConstraintSolution; break;
  case 2:  sol = new ColRowSolution; break;
  default: sol = new ColVariableSolution;
  }

 if( ! emptys )
  sol->read( this );

 return( sol );
 }

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS FOR SAVING THE UnitBlock -------------------*/
/*--------------------------------------------------------------------------*/

void UnitBlock::serialize( netCDF::NcGroup & group ) const
{
 Block::serialize( group );

 group.addDim( "TimeHorizon", f_time_horizon );

 if( ( f_number_intervals > 1 ) && ( f_number_intervals < f_time_horizon ) ) {
  auto NI = group.addDim( "NumberIntervals" , f_number_intervals );

  ::serialize( group , "ChangeInterval", netCDF::NcUint64() , NI ,
	       v_change_intervals );
  }
 }

/*--------------------------------------------------------------------------*/

template< typename T >
void UnitBlock::decompress_vector( std::vector< T > & v ) {
 if ( v.empty() )
  return;

 if( v.size() == 1 ) {
  // The given vector has a single element. Thus, for each time instant, the
  // value is equal to that single given element.
  v.resize( f_time_horizon , v[ 0 ] );
 }
 else if( v.size() < f_time_horizon ) {
  // Since the number of elements is greater than 1 and less than the time
  // horizon, it must be equal to the number of change intervals.
  if( v.size() != v_change_intervals.size() ) {
   throw ( std::logic_error
    ( "UnitBlock::decompress_vector: invalid number of elements"
      " (" + std::to_string( v.size() ) + ") for some variable. It "
                                          "should be equal to the number of change intervals (" +
      std::to_string( v_change_intervals.size() ) + ")" ) );
  }

  // For each time instant t, the value associated with time t is equal to
  // given_vector[ k ], where k is such that t belongs to the closed interval
  // [i_{k-1} + 1, i_k] and i_k is the k-th element of v_change_intervals
  // (starting from k = 0) and i_{-1} = -1 by definition. We resize the vector
  // so that its size becomes f_time_horizon and copy the given data.

  std::vector< T > given_vector = v;
  v.resize( f_time_horizon );
  Index t = 0;
  for( Index k = 0 ; k < v_change_intervals.size() ; ++k ) {
   auto upper_endpoint = v_change_intervals[ k ];
   if( k == v_change_intervals.size() - 1 )
    // The upper endpoint of the last interval must be time_horizon - 1. Since
    // it may not be provided in v_change_intervals (the value for the last
    // element of v_change_intervals is not required), we manually set it
    // here.
    upper_endpoint = f_time_horizon - 1;
   for( ; t <= upper_endpoint ; ++t ) {
    v[ t ] = given_vector[ k ];
   }
  }
 }
}

/*--------------------------------------------------------------------------*/
/*---------------------- End File UnitBlock.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
