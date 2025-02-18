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
 Index wsol = 15;
 if( ( ! csolc ) && f_BlockConfig )
  csolc = f_BlockConfig->f_solution_Configuration;

 if( auto config = dynamic_cast< SimpleConfiguration< int > * >( csolc ) )
  wsol = config->f_value;

 auto sol = new_Solution();

 using mad2 = boost::multi_array< double , 2 >;

 if( wsol & 1 )
  sol->v_active_power.resize(
       mad2::extent_gen()[ get_number_generators() ][ get_time_horizon() ] );

 // note: we assume that either all generators have commitment, or none has
 if( ( wsol & 2 ) && get_commitment( 0 ) )
  sol->v_commitment.resize(
       mad2::extent_gen()[ get_number_generators() ][ get_time_horizon() ] );

 // note: we assume that either all generators have primary, or none has
 if( ( wsol & 4 ) && get_primary_spinning_reserve( 0 ) )
  sol->v_primary_reserve.resize(
       mad2::extent_gen()[ get_number_generators() ][ get_time_horizon() ] );

 // note: we assume that either all generators have secondary, or none has
 if( ( wsol & 8 ) && get_secondary_spinning_reserve( 0 ) )
  sol->v_secondary_reserve.resize(
       mad2::extent_gen()[ get_number_generators() ][ get_time_horizon() ] );

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

 if( ( f_number_intervals > 1 ) &&
     ( f_number_intervals < f_time_horizon ) ) {
  auto NI = group.addDim( "NumberIntervals" , f_number_intervals );

  ::serialize( group , "ChangeInterval" , netCDF::NcUint64() , NI ,
               v_change_intervals );
  }
 }

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS OF UnitBlockSolution ------------------------*/
/*--------------------------------------------------------------------------*/

void UnitBlockSolution::deserialize( const netCDF::NcGroup & group )
{
 // "TimeHorizon" is mandatory- - - - - - - - - - - - - - - - - - - - - - - -
 ::deserialize_dim( group , "TimeHorizon" , f_time_horizon , false );

 if( ! ::deserialize_dim( group , "NumberGenerators" , f_time_horizon ,
			  true ) )
  f_number_generators = 1;

 // deserialize the Active Power- - - - - - - - - - - - - - - - - - - - - - -
 ::deserialize< double , 2 >( group , "ActivePower" , v_active_power ,
			      false );

 // deserialize the Commitment- - - - - - - - - - - - - - - - - - - - - - - -
 ::deserialize< double , 2 >( group , "Commitment" , v_commitment , true );

 // deserialize the Primary Reserve - - - - - - - - - - - - - - - - - - - - -
 ::deserialize< double , 2 >( group , "PrimaryReserve" , v_primary_reserve ,
			      true );

 // deserialize the Secondary Reserve- - - - - - - - - - - - - - - - - - - -
 ::deserialize< double , 2 >( group , "SecondaryReserve" ,
			      v_secondary_reserve , true );

 }  // end( UnitBlockSolution::deserialize )

/*--------------------------------------------------------------------------*/

void UnitBlockSolution::read( const Block * block )
{
 auto UB = dynamic_cast< const UnitBlock * >( block );
 if( ! UB )
  throw( std::invalid_argument(
		     "UnitBlockSolution::read: block is not a UnitBlock" ) );

 f_time_horizon = UB->get_time_horizon();
 f_number_generators = UB->get_number_generators();

 if( ! v_active_power.empty() )
  // read the active power variables - - - - - - - - - - - - - - - - - - - -
  for( Index i = 0 ; i < f_number_generators ; ++i ) {
   auto APi = UB->get_active_power( i );
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    v_active_power[ i ][ t ] = APi[ t ].get_value();
   }

 if( ! v_commitment.empty() )
  // read the commitment variables - - - - - - - - - - - - - - - - - - - - -
  for( Index i = 0 ; i < f_number_generators ; ++i )
   if( auto Ci = UB->get_commitment( i ) )
    for( Index t = 0 ; t < f_time_horizon ; ++t )
     v_commitment[ i ][ t ] = Ci[ t ].get_value();

 if( ! v_primary_reserve.empty() )
  // read the primary reserve variables- - - - - - - - - - - - - - - - - - -
  for( Index i = 0 ; i < f_number_generators ; ++i )
   if( auto PRi = UB->get_primary_spinning_reserve( i ) )
    for( Index t = 0 ; t < f_time_horizon ; ++t )
     v_primary_reserve[ i ][ t ] = PRi[ t ].get_value();

 if( ! v_secondary_reserve.empty() )
  // read the secondary reserve variables- - - - - - - - - - - - - - - - - -
  for( Index i = 0 ; i < f_number_generators ; ++i )
   if( auto SRi = UB->get_secondary_spinning_reserve( i ) )
    for( Index t = 0 ; t < f_time_horizon ; ++t )
     v_secondary_reserve[ i ][ t ] = SRi[ t ].get_value();

 }  // end( UnitBlockSolution::read )

/*--------------------------------------------------------------------------*/

void UnitBlockSolution::write( Block * block )
{
 auto UB = dynamic_cast< const UnitBlock * >( block );
 if( ! UB )
  throw( std::invalid_argument(
		   "UnitBlockSolution::write: block is not a UnitBlock" ) );

 if( f_time_horizon != UB->get_time_horizon() )
  throw( std::invalid_argument(
		   "UnitBlockSolution::write: inconsistent time horizon" ) );

 if( f_number_generators != UB->get_number_generators() )
  throw( std::invalid_argument(
	      "UnitBlockSolution::write: inconsistent generators number" ) );

 if( ! v_active_power.empty() )
  // write the active power variables- - - - - - - - - - - - - - - - - - - -
  for( Index i = 0 ; i < f_number_generators ; ++i ) {
   auto APi = UB->get_active_power( i );
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    APi[ t ].set_value( v_active_power[ i ][ t ] );
   }

 if( ! v_commitment.empty() )
  // write the commitment variables- - - - - - - - - - - - - - - - - - - - -
  for( Index i = 0 ; i < f_number_generators ; ++i )
   if( auto Ci = UB->get_commitment( i ) )
    for( Index t = 0 ; t < f_time_horizon ; ++t )
     Ci[ t ].set_value( v_commitment[ i ][ t ] );
   else
    throw( std::invalid_argument(
	  "UnitBlockSolution::write: provided non-existent commitment" ) );
    
 if( ! v_primary_reserve.empty() )
  // write the primary reserve variables - - - - - - - - - - - - - - - - - -
  for( Index i = 0 ; i < f_number_generators ; ++i )
   if( auto PRi = UB->get_primary_spinning_reserve( i ) )
    for( Index t = 0 ; t < f_time_horizon ; ++t )
     PRi[ t ].set_value( v_primary_reserve[ i ][ t ] );
   else
    throw( std::invalid_argument(
	  "UnitBlockSolution::write: provided non-existent primary" ) );

 if( ! v_secondary_reserve.empty() )
  // write the secondary reserve variables - - - - - - - - - - - - - - - - -
  for( Index i = 0 ; i < f_number_generators ; ++i )
   if( auto SRi = UB->get_secondary_spinning_reserve( i ) )
    for( Index t = 0 ; t < f_time_horizon ; ++t )
     SRi[ t ].set_value( v_secondary_reserve[ i ][ t ] );
   else
    throw( std::invalid_argument(
	  "UnitBlockSolution::write: provided non-existent secondary" ) );

 }  // end( UnitBlockSolution::write )

/*--------------------------------------------------------------------------*/

void UnitBlockSolution::serialize( const netCDF::NcGroup & group )
{
 // "TimeHorizon" is mandatory- - - - - - - - - - - - - - - - - - - - - - - -
 auto th = group.addDim( "TimeHorizon" , f_time_horizon );

 netCDF::NcDim ng;
 if( f_number_generators > 1 )
  ng = group.addDim( "NumberGenerators" , f_number_generators );

 // serialize the Active Power- - - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_active_power.empty() )
  ::serialize< double , 2 >( group , "ActivePower" , netCDF::NcDouble() ,
			     { ng , th } , v_active_power.data() ,
			     { f_number_generators , f_time_horizon } );

 // serialize the Commitment- - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_commitment.empty() )
  ::serialize< double , 2 >( group , "Commitment" , netCDF::NcDouble() ,
			     { ng , th } , v_commitment.data() ,
			     { f_number_generators , f_time_horizon } );


 // serialize the Primary Reserve - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_primary_reserve.empty() )
  ::serialize< double , 2 >( group , "Commitment" , netCDF::NcDouble() ,
			     { ng , th } , v_primary_reserve.data() ,
			     { f_number_generators , f_time_horizon } );

 // serialize the Secondary Reserve - - - - - - - - - - - - - - - - - - - - -
 if( ! v_secondary_reserve.empty() )
  ::serialize< double , 2 >( group , "Commitment" , netCDF::NcDouble() ,
			     { ng , th } , v_secondary_reserve.data() ,
			     { f_number_generators , f_time_horizon } );

 }  // end( UnitBlockSolution::serialize )

/*--------------------------------------------------------------------------*/

UnitBlockSolution * UnitBlockSolution::scale( double factor ) const
{
 auto sol = clone();

 if( factor == 1 )
  return( sol );

 if( ! v_active_power.empty() )
  for( Index i = 0 ; i < f_number_generators ; ++i )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    sol->v_active_power[ i ][ t ] *= factor;

 if( ! v_commitment.empty() )
  for( Index i = 0 ; i < f_number_generators ; ++i )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    sol->v_commitment[ i ][ t ] *= factor;

 if( ! v_primary_reserve.empty() )
  for( Index i = 0 ; i < f_number_generators ; ++i )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    sol->v_primary_reserve[ i ][ t ] *= factor;

 if( ! v_secondary_reserve.empty() )
  for( Index i = 0 ; i < f_number_generators ; ++i )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    sol->v_secondary_reserve[ i ][ t ] *= factor;

 return( sol );

 }  // end( UnitBlockSolution::scale )

/*--------------------------------------------------------------------------*/

void UnitBlockSolution::sum( const Solution * solution , double multiplier )
{
 auto UBS = dynamic_cast< const UnitBlockSolution * >( solution );
 if( ! UBS )
  throw( std::invalid_argument(
	      "UnitBlockSolution::sum: solution not a UnitBlockSolution" ) );

 if( f_time_horizon != UBS->f_time_horizon )
  throw( std::invalid_argument(
		     "UnitBlockSolution::sum: inconsistent time horizon" ) );

 if( f_number_generators != UBS->f_number_generators )
  throw( std::invalid_argument(
	        "UnitBlockSolution::sum: inconsistent generators number" ) );

 if( ! v_active_power.empty() )
  for( Index i = 0 ; i < f_number_generators ; ++i )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    v_active_power[ i ][ t ] += UBS->v_active_power[ i ][ t ] * multiplier;

 if( ! v_commitment.empty() )
  for( Index i = 0 ; i < f_number_generators ; ++i )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    v_commitment[ i ][ t ] += UBS->v_commitment[ i ][ t ] * multiplier;

 if( ! v_primary_reserve.empty() )
  for( Index i = 0 ; i < f_number_generators ; ++i )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    v_primary_reserve[ i ][ t ] +=
     UBS->v_primary_reserve[ i ][ t ] * multiplier;

 if( ! v_secondary_reserve.empty() )
  for( Index i = 0 ; i < f_number_generators ; ++i )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    v_secondary_reserve[ i ][ t ] +=
     UBS->v_secondary_reserve[ i ][ t ] * multiplier;

 }  // end( UnitBlockSolution::sum )

/*--------------------------------------------------------------------------*/

UnitBlockSolution * UnitBlockSolution::clone( bool empty ) const
{
 auto * sol = new_Solution();

 if( ! empty ) {
  sol->f_time_horizon = f_time_horizon;
  sol->f_number_generators = f_number_generators;

  sol->v_active_power = v_active_power;
  sol->v_commitment = v_commitment;
  sol->v_primary_reserve = v_primary_reserve;
  sol->v_secondary_reserve = v_secondary_reserve;
  }

 return( sol );

 }  // end( UnitBlockSolution::clone )

/*--------------------------------------------------------------------------*/
/*---------------------- End File UnitBlock.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
