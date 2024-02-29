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
 * \copyright &copy; by Antonio Frangioni, Ali Ghezelsoflu,
 *                      Rafael Durbano Lobato
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
    throw( std::invalid_argument(
     classname() + "::deserialize: TimeHorizon is not present in the "
                   "netCDF input and UnitBlock does not have a father." ) );
  }
 } else {
  // dimension TimeHorizon is present in the netCDF input

  auto th = TimeHorizon.getSize();
  if( f_time_horizon == 0 )
   this->set_time_horizon( th );
  else if( f_time_horizon != th )
   throw( std::logic_error(
    classname() + "::deserialize: TimeHorizon is not present in the "
                  "netCDF. The (nonzero) time horizon of UnitBlock is different "
                  "from that of its father, but they should be equal." ) );
 }
}

/*--------------------------------------------------------------------------*/

void UnitBlock::deserialize_change_intervals( const netCDF::NcGroup & group )
{
 if( ! ::deserialize_dim( group , "NumberIntervals" , f_number_intervals ) )
  f_number_intervals = 1;
 else
  if( ( f_number_intervals < 1 ) || ( f_number_intervals > f_time_horizon ) )
   throw( std::invalid_argument(
    classname() + "::deserialize: NumberIntervals not between 1 and "
                  "TimeHorizon." ) );

 if( ( f_number_intervals > 1 ) && ( f_number_intervals < f_time_horizon ) ) {
  ::deserialize( group , "ChangeIntervals" , f_number_intervals ,
                 v_change_intervals );

  // Check that the numbers are ordered in increasing sense. Notice that the
  // upper endpoint of the last interval must necessarily be f_time_horizon -
  // 1. Since it is not required to be provided, we set it.
  v_change_intervals.back() = f_time_horizon - 1;

  for( Index k = 0 ; k < v_change_intervals.size() ; ++k ) {
   const auto t = v_change_intervals[ k ];
   if( ! ( ( t < f_time_horizon ) &&
           ( ( k == 0 ) || ( t > v_change_intervals[ k - 1 ] ) ) ) )
    throw( std::invalid_argument(
     classname() + "::deserialize: invalid value in ChangeIntervals: " +
     std::to_string( t ) + ". All values must be between 0 and " +
     "TimeHorizon - 1 and in strictly increasing order." ) );
  }
 } else
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

void UnitBlock::scale( MF_dbl_it values ,
                       Range rng ,
                       c_ModParam issuePMod ,
                       c_ModParam issueAMod )
{
 if( rng.first >= rng.second )
  return;  // An empty Range was given: no operation is performed.

 Subset subset;

 if( rng.second == Inf< Index >() ) {
  // If we decide to scale the generators individually rather than the whole
  // unit, then, when rng.second is Inf< Index >(), we could interpret it as
  // changing the scale factor of all generators and the vector containing the
  // scale factor would be expected to have size at least equal to the number
  // of generators. In this case, the subset would have size equal to the
  // number of generators. Alternatively, we could have scale_generators() and
  // leave scale() for scaling the whole unit.
  subset.resize( 1 , 0 );
 }
 else {
  subset.resize( rng.second - rng.first );
  std::iota( subset.begin() , subset.end() , rng.first );
 }

 scale( values , std::move( subset ) , true , issuePMod , issueAMod );
}

/*--------------------------------------------------------------------------*/

void UnitBlock::scale( double scale_factor ,
                       c_ModParam issuePMod ,
                       c_ModParam issueAMod )
{
 Subset subset = { 0 };
 std::vector< double > values = { scale_factor };
 scale( values.cbegin() , std::move( subset ) , true , issuePMod , issueAMod );
}

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

 Solution * sol;
 switch( solution_type ) {
  case( 1 ):
   sol = new RowConstraintSolution;
   break;
  case( 2 ):
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
/*--------------------- METHODS FOR SAVING THE UnitBlock -------------------*/
/*--------------------------------------------------------------------------*/

void UnitBlock::serialize( netCDF::NcGroup & group ) const
{
 Block::serialize( group );

 group.addDim( "TimeHorizon" , f_time_horizon );

 if( ( f_number_intervals > 1 ) && ( f_number_intervals < f_time_horizon ) ) {
  auto NI = group.addDim( "NumberIntervals" , f_number_intervals );

  ::serialize( group , "ChangeInterval" , netCDF::NcUint64() , NI ,
               v_change_intervals );
 }
}

/*--------------------------------------------------------------------------*/
/*---------------------- End File UnitBlock.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
