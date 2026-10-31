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

// register UnitBlockSolution to the Solution factory
SMSpp_insert_in_factory_cpp_0( UnitBlockSolution );

/*--------------------------------------------------------------------------*/
/*--------------------------- STATIC FUNCTIONS -----------------------------*/
/*--------------------------------------------------------------------------*/

template< class T , std::size_t K >
static void copy_multi_array( boost::multi_array< T , K > & to ,
			      const boost::multi_array< T , K > & from )
{
 std::vector< size_t > extent;
 auto shape = from.shape();
 extent.assign( shape , shape + from.num_dimensions() );
 to.resize( extent );
 to = from;
 }

/*--------------------------------------------------------------------------*/
/*--------------------------- METHODS OF UnitBlock -------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void UnitBlock::issue_unit_mod( sp_Mod mod , ModParam issuePMod )
{
 if( issue_pmod( issuePMod , mod->changes() ) )
  Block::add_Modification( std::move( mod ) ,
                           Observer::par2chnl( issuePMod ) );
 else
  if( auto fb = get_f_Block() ) {
   const auto umod = dynamic_cast< const UnitBlockMod * >( mod.get() );
   if( umod && umod->changes_power_range() )
    fb->add_Modification( std::move( mod ) ,
                          Observer::par2chnl( issuePMod ) );
   }
 }

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
                   "netCDF input and UnitBlock does not have a father" ) );
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
                  "from that of its father, but they should be equal" ) );
  }
 }

/*--------------------------------------------------------------------------*/

void UnitBlock::deserialize_change_intervals( const netCDF::NcGroup & group )
{
 if( ! deserialize_dim( group , "NumberIntervals" , f_number_intervals ) )
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
  }
 else
  v_change_intervals.clear();
 }

/*--------------------------------------------------------------------------*/

void UnitBlock::deserialize_outputs( const netCDF::NcGroup & group )
{
 static const std::string who = "::deserialize: ";

 f_number_outputs = 0;
 f_output_generators = 0;
 v_output_carrier.clear();
 v_output_efficiency.resize( boost::extents[ 0 ][ 0 ][ 0 ] );

 // without "NumberOutputs", the only output kElectricity, with efficiency 1;
 // the dimension is that of the group, not of the UCBlock above it
 Index K = 0;
 if( ! deserialize_local_dim( group , "NumberOutputs" , K ) )
  return;
 if( ! K )
  throw( std::invalid_argument( classname() + who + "NumberOutputs is 0" ) );

 // OutputCarrier, [ outputs ] or [ generators ][ outputs ]
 auto oc = group.getVar( "OutputCarrier" );
 if( oc.isNull() )
  throw( std::invalid_argument( classname() + who + "OutputCarrier "
				"mandatory with NumberOutputs" ) );
 const auto ocd = oc.getDimCount();
 if( ( ocd < 1 ) || ( ocd > 2 ) || ( oc.getDim( ocd - 1 ).getSize() != K ) )
  throw( std::invalid_argument( classname() + who + "OutputCarrier must "
				"be indexed over [ generators , ] "
				"NumberOutputs" ) );
 const Index G = ( ocd == 2 ) ? oc.getDim( 0 ).getSize() : 1;
 if( ! G )
  throw( std::invalid_argument( classname() + who + "OutputCarrier over 0 "
				"generators" ) );
 std::vector< char * > names( G * K );
 oc.getVar( names.data() );
 v_output_carrier.assign( names.begin() , names.end() );
 for( auto n : names )
  free( n );
 f_number_outputs = K;
 f_output_generators = G;

 // OutputEfficiency, [ 1 or time ][ outputs ] or [ 1 or time ][ 1 or
 // generators ][ outputs ]
 auto oe = group.getVar( "OutputEfficiency" );
 if( oe.isNull() )
  throw( std::invalid_argument( classname() + who + "OutputEfficiency "
				"mandatory with NumberOutputs" ) );
 const auto oed = oe.getDimCount();
 if( ( oed < 2 ) || ( oed > 3 ) || ( oe.getDim( oed - 1 ).getSize() != K ) )
  throw( std::invalid_argument( classname() + who + "OutputEfficiency must "
				"be indexed over TimeHorizon , [ generators , ] "
				"NumberOutputs" ) );
 const Index Tn = oe.getDim( 0 ).getSize();
 const Index Gn = ( oed == 3 ) ? oe.getDim( 1 ).getSize() : 1;
 if( ( ( Tn != 1 ) && ( Tn != f_time_horizon ) ) || ( ! Gn ) )
  throw( std::invalid_argument( classname() + who + "OutputEfficiency must "
				"have size [ 1 or TimeHorizon , [ 1 or "
				"generators , ] NumberOutputs ]" ) );
 v_output_efficiency.resize( boost::extents[ Tn ][ Gn ][ K ] );
 oe.getVar( v_output_efficiency.data() );
 }

/*--------------------------------------------------------------------------*/

void UnitBlock::check_outputs( void ) const
{
 const auto G = get_number_generators();
 if( ( ( f_output_generators > 1 ) && ( f_output_generators != G ) ) ||
     ( ( ! v_output_efficiency.empty() ) &&
       ( v_output_efficiency.shape()[ 1 ] > 1 ) &&
       ( v_output_efficiency.shape()[ 1 ] != G ) ) )
  throw( std::invalid_argument( classname() + "::deserialize: the outputs "
				"are given over a number of generators other "
				"than " + std::to_string( G ) ) );
 }

/*--------------------------------------------------------------------------*/

void UnitBlock::deserialize( const netCDF::NcGroup & group )
{
 Block::deserialize( group );

 // "Scale", the number of copies, is both the weights [see replicate()]
 double copies;
 if( ::deserialize( group , copies , "Scale" ) ) {
  if( ( ! std::isfinite( copies ) ) || ( copies < 0 ) )
   throw( std::invalid_argument( classname() + "::deserialize: Scale must "
                                 "be finite and not negative" ) );
  f_objective_weight = f_variable_weight = copies;
  }

 deserialize_time_horizon( group );
 deserialize_change_intervals( group );
 deserialize_outputs( group );
 }

/*--------------------------------------------------------------------------*/

std::vector< std::string > UnitBlock::expected_dims( void ) const {
 auto ret = Block::expected_dims();
 ret.push_back( "TimeHorizon" );
 ret.push_back( "NumberIntervals" );
 ret.push_back( "NumberOutputs" );

 return( ret );
 }

/*--------------------------------------------------------------------------*/

std::vector< std::string > UnitBlock::expected_vars( void ) const {
 auto ret = Block::expected_vars();
 ret.push_back( "ChangeIntervals" );
 ret.push_back( "OutputCarrier" );
 ret.push_back( "OutputEfficiency" );

 return( ret );
 }

/*--------------------------------------------------------------------------*/
/*------------------ METHODS FOR MODIFYING THE UnitBlock -------------------*/
/*--------------------------------------------------------------------------*/

void UnitBlock::generate_abstract_variables( Configuration * stvv )
{
 if( variables_generated() )  // variables have already been generated
  return;                     // nothing to do

 // Reactive Power Variables- - - - - - - - - - - - - - - - - - - - - - - - -

}  // end( UnitBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void UnitBlock::set_weight( double k , c_ModParam issueMod )
{
 const auto copies = get_variable_weight();
 Block::set_weight( k , issueMod );
 if( get_variable_weight() == copies )
  return;

 /* The UCBlock may not be the father: while a LagrangianDualSolver is
  * attached, the father is the LagBFunction holding the unit alone, whose
  * father is the UCBlock, and the UnitBlock may be in a unit of the UCBlock
  * (the HydroUnitBlock of a HydroSystemUnitBlock). What the UCBlock is told
  * of is its sub-Block that contains this one, a Function holding it not
  * being one. */
 const Block * unit = this;
 for( auto b = get_f_Block() ; b ; b = b->get_f_Block() ) {
  if( auto ucb = dynamic_cast< UCBlock * >( b ) ) {
   ucb->weight_changed( unit );
   break;
   }
  if( ! dynamic_cast< const Function * >( b ) )
   unit = b;
  }
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

 auto sz = boost::multi_array< double , 2 >::extent_gen()
                           [ get_number_generators() ][ get_time_horizon() ];
 /* Which parts the Solution is made of is asked to the *data*, not to the
  * Variable: a Solution has to be possible even if the abstract
  * representation has not been generated, in which case the Variable are not
  * there and looking at them would say that the UnitBlock has nothing. */

 if( wsol & 1 ) {
  sol->v_active_power.resize( sz );
  if( has_reactive_power() )
   sol->v_reactive_power.resize( sz );
  }

 if( ( wsol & 2 ) && has_commitment() )
  sol->v_commitment.resize( sz );

 if( ( wsol & 4 ) && has_primary_reserve() )
  sol->v_primary_reserve.resize( sz );

 if( ( wsol & 8 ) && has_secondary_reserve() )
  sol->v_secondary_reserve.resize( sz );

 if( ! emptys )
  sol->read( this );

 return( sol );
 }

/*--------------------------------------------------------------------------*/

UnitBlockSolution * UnitBlock::new_Solution( void ) const {
  return( new UnitBlockSolution() );
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

  ::serialize( group , "ChangeIntervals" , netCDF::NcUint64() , NI ,
               v_change_intervals );
  }

 serialize_outputs( group );
 }

/*--------------------------------------------------------------------------*/

void UnitBlock::serialize_outputs( netCDF::NcGroup & group ) const
{
 if( v_output_carrier.empty() )
  return;  // the only output kElectricity, with efficiency 1

 auto NO = group.addDim( "NumberOutputs" , f_number_outputs );
 const auto shp = v_output_efficiency.shape();

 // the dimension of the generators, if anything is given per generator
 netCDF::NcDim NG;
 if( ( f_output_generators > 1 ) || ( shp[ 1 ] > 1 ) )
  NG = group.addDim( "NumberOutputGenerators" ,
		     std::max< Index >( f_output_generators , shp[ 1 ] ) );

 if( f_output_generators > 1 ) {
  std::vector< const char * > names( v_output_carrier.size() );
  std::transform( v_output_carrier.begin() , v_output_carrier.end() ,
		  names.begin() ,
		  []( const std::string & n ) { return( n.c_str() ); } );
  group.addVar( "OutputCarrier" , netCDF::NcString() , { NG , NO } ).putVar(
							       names.data() );
  }
 else
  ::serialize( group , "OutputCarrier" , NO , v_output_carrier );

 // the first dimension is a singleton if the efficiencies do not depend on
 // time
 auto TH = group.getDim( "TimeHorizon" );
 if( shp[ 1 ] > 1 )
  ::serialize( group , "OutputEfficiency" , netCDF::NcDouble() ,
	       { TH , NG , NO } , v_output_efficiency , false , true );
 else {
  MAdouble eff( ( MAdouble_ext()[ shp[ 0 ] ][ f_number_outputs ] ) );
  std::copy( v_output_efficiency.data() ,
	     v_output_efficiency.data() + v_output_efficiency.num_elements() ,
	     eff.data() );
  ::serialize( group , "OutputEfficiency" , netCDF::NcDouble() ,
	       { TH , NO } , eff , false , true );
  }
 }

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS OF UnitBlockSolution ------------------------*/
/*--------------------------------------------------------------------------*/

void UnitBlockSolution::deserialize( const netCDF::NcGroup & group )
{
 // "TimeHorizon" is mandatory- - - - - - - - - - - - - - - - - - - - - - - -
 deserialize_dim( group , "TimeHorizon" , f_time_horizon , false );

 // with no "NumberGenerators" there is one generator, and the variables
 // are indexed over the instants only; a scalar gives the same value to
 // every generator and instant
 if( ! deserialize_dim( group , "NumberGenerators" , f_number_generators ,
			true ) )
  f_number_generators = 1;

 // deserialize the Active Power - - - - - - - - - - - - - - - - - - - - - -
 ::deserialize< double , 2 >( group , "ActivePower" ,
                              { f_number_generators , f_time_horizon } ,
                              v_active_power , true , true );

 ::deserialize< double , 2 >( group , "ReactivePower" ,
                              { f_number_generators , f_time_horizon } ,
                              v_reactive_power , true , true );

 // deserialize the Commitment - - - - - - - - - - - - - - - - - - - - - - -
 ::deserialize< double , 2 >( group , "Commitment" ,
                              { f_number_generators , f_time_horizon } ,
                              v_commitment , true , true );

 // deserialize the Primary Reserve- - - - - - - - - - - - - - - - - - - - -
 ::deserialize< double , 2 >( group , "PrimaryReserve" ,
                              { f_number_generators , f_time_horizon } ,
                              v_primary_reserve , true , true );

 // deserialize the Secondary Reserve - - - - - - - - - - - - - - - - - - -
 ::deserialize< double , 2 >( group , "SecondaryReserve" ,
                              { f_number_generators , f_time_horizon } ,
                              v_secondary_reserve , true , true );
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
   auto APi = UB->get_const_active_power( i );
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    v_active_power[ i ][ t ] = APi[ t ].get_value();
   }

 if( ! v_reactive_power.empty() )
  // read the active power variables - - - - - - - - - - - - - - - - - - - -
  for( Index i = 0 ; i < f_number_generators ; ++i ) {
   if( auto RPi = UB->get_const_reactive_power( i ) )
    for( Index t = 0 ; t < f_time_horizon ; ++t )
     v_reactive_power[ i ][ t ] = RPi[ t ].get_value();
   else
    for( Index t = 0 ; t < f_time_horizon ; ++t )
     v_reactive_power[ i ][ t ] = 0;
   }

 if( ! v_commitment.empty() )
  // read the commitment variables - - - - - - - - - - - - - - - - - - - - -
  for( Index i = 0 ; i < f_number_generators ; ++i )
   if( auto Ci = UB->get_const_commitment( i ) )
    for( Index t = 0 ; t < f_time_horizon ; ++t )
     v_commitment[ i ][ t ] = Ci[ t ].get_value();

 if( ! v_primary_reserve.empty() )
  // read the primary reserve variables- - - - - - - - - - - - - - - - - - -
  for( Index i = 0 ; i < f_number_generators ; ++i )
   if( auto PRi = UB->get_const_primary_spinning_reserve( i ) )
    for( Index t = 0 ; t < f_time_horizon ; ++t )
     v_primary_reserve[ i ][ t ] = PRi[ t ].get_value();

 if( ! v_secondary_reserve.empty() )
  // read the secondary reserve variables- - - - - - - - - - - - - - - - - -
  for( Index i = 0 ; i < f_number_generators ; ++i )
   if( auto SRi = UB->get_const_secondary_spinning_reserve( i ) )
    for( Index t = 0 ; t < f_time_horizon ; ++t )
     v_secondary_reserve[ i ][ t ] = SRi[ t ].get_value();

 }  // end( UnitBlockSolution::read )

/*--------------------------------------------------------------------------*/

void UnitBlockSolution::write( Block * block )
{
 auto UB = dynamic_cast< UnitBlock * >( block );
 if( ! UB )
  throw( std::invalid_argument(
		   "UnitBlockSolution::write: block is not a UnitBlock" ) );

 if( f_time_horizon != UB->get_time_horizon() )
  throw( std::invalid_argument(
		   "UnitBlockSolution::write: inconsistent time horizon" ) );

 if( f_number_generators != UB->get_number_generators() )
  throw( std::invalid_argument(
	      "UnitBlockSolution::write: inconsistent generators number" ) );

 /* Which parts the Solution is made of is decided by the data [see
  * UnitBlock::get_Solution()], so a Solution can carry a part that some
  * generator, or the whole UnitBlock, has no Variable for: as read() does,
  * those generators are skipped, writing being possible on what is there. */

 if( ! v_active_power.empty() )
  // write the active power variables- - - - - - - - - - - - - - - - - - - -
  for( Index i = 0 ; i < f_number_generators ; ++i )
   if( auto APi = UB->get_active_power( i ) )
    for( Index t = 0 ; t < f_time_horizon ; ++t )
     APi[ t ].set_value( v_active_power[ i ][ t ] );

 if( ! v_reactive_power.empty() )
  // write the reactive power variables- - - - - - - - - - - - - - - - - - -
  for( Index i = 0 ; i < f_number_generators ; ++i )
   if( auto RPi = UB->get_reactive_power( i ) )
    for( Index t = 0 ; t < f_time_horizon ; ++t )
     RPi[ t ].set_value( v_reactive_power[ i ][ t ] );

 if( ! v_commitment.empty() )
  // write the commitment variables- - - - - - - - - - - - - - - - - - - - -
  for( Index i = 0 ; i < f_number_generators ; ++i )
   if( auto Ci = UB->get_commitment( i ) )
    for( Index t = 0 ; t < f_time_horizon ; ++t )
     Ci[ t ].set_value( v_commitment[ i ][ t ] );

 if( ! v_primary_reserve.empty() )
  // write the primary reserve variables - - - - - - - - - - - - - - - - - -
  for( Index i = 0 ; i < f_number_generators ; ++i )
   if( auto PRi = UB->get_primary_spinning_reserve( i ) )
    for( Index t = 0 ; t < f_time_horizon ; ++t )
     PRi[ t ].set_value( v_primary_reserve[ i ][ t ] );

 if( ! v_secondary_reserve.empty() )
  // write the secondary reserve variables - - - - - - - - - - - - - - - - -
  for( Index i = 0 ; i < f_number_generators ; ++i )
   if( auto SRi = UB->get_secondary_spinning_reserve( i ) )
    for( Index t = 0 ; t < f_time_horizon ; ++t )
     SRi[ t ].set_value( v_secondary_reserve[ i ][ t ] );

 }  // end( UnitBlockSolution::write )

/*--------------------------------------------------------------------------*/

void UnitBlockSolution::serialize( netCDF::NcGroup & group ) const
{
 Solution::serialize( group );

 // "TimeHorizon" is mandatory- - - - - - - - - - - - - - - - - - - - - - - -
 auto th = group.addDim( "TimeHorizon" , f_time_horizon );

 netCDF::NcDim ng;
 if( f_number_generators > 1 ) {
  ng = group.addDim( "NumberGenerators" , f_number_generators );

  // serialize the Active Power - - - - - - - - - - - - - - - - - - - - - - -
  ::serialize< double , 2 >( group , "ActivePower" , netCDF::NcDouble() ,
			     { ng , th } , v_active_power );

  // serialize the Reactive Power - - - - - - - - - - - - - - - - - - - - - -
  ::serialize< double , 2 >( group , "ReactivePower" , netCDF::NcDouble() ,
			     { ng , th } , v_reactive_power );

  // serialize the Commitment - - - - - - - - - - - - - - - - - - - - - - - -
  ::serialize< double , 2 >( group , "Commitment" , netCDF::NcDouble() ,
			     { ng , th } , v_commitment );

  // serialize the Primary Reserve- - - - - - - - - - - - - - - - - - - - - -
  ::serialize< double , 2 >( group , "PrimaryReserve" , netCDF::NcDouble() ,
			     { ng , th } , v_primary_reserve );

  // serialize the Secondary Reserve- - - - - - - - - - - - - - - - - - - - -
  ::serialize< double , 2 >( group , "SecondaryReserve" , netCDF::NcDouble() ,
			     { ng , th } , v_secondary_reserve );
  }
 else {
  // serialize the Active Power - - - - - - - - - - - - - - - - - - - - - - -
  if( ! v_active_power.empty() )
   group.addVar( "ActivePower" , netCDF::NcDouble() , th ).putVar(
			{ 0 } , { f_time_horizon } , v_active_power.data() );

  // serialize the Reactive Power - - - - - - - - - - - - - - - - - - - - - -
  if( ! v_reactive_power.empty() )
   group.addVar( "ReactivePower" , netCDF::NcDouble() , th ).putVar(
		      { 0 } , { f_time_horizon } , v_reactive_power.data() );

  // serialize the Commitment - - - - - - - - - - - - - - - - - - - - - - - -
  if( ! v_commitment.empty() )
   group.addVar( "Commitment" , netCDF::NcDouble() , th ).putVar(
			  { 0 } , { f_time_horizon } , v_commitment.data() );

  // serialize the Primary Reserve- - - - - - - - - - - - - - - - - - - - - -
  if( ! v_primary_reserve.empty() )
   group.addVar( "PrimaryReserve" , netCDF::NcDouble() , th ).putVar(
		     { 0 } , { f_time_horizon } , v_primary_reserve.data() );

  // serialize the Secondary Reserve- - - - - - - - - - - - - - - - - - - - -
  if( ! v_secondary_reserve.empty() )
   group.addVar( "SecondaryReserve" , netCDF::NcDouble() , th ).putVar(
		  { 0 } , { f_time_horizon } , v_secondary_reserve.data() );
  }
 }  // end( UnitBlockSolution::serialize )

/*--------------------------------------------------------------------------*/

UnitBlockSolution * UnitBlockSolution::scale( double factor ) const
{
 auto sol = clone();

 if( factor == 1 )
  return( sol );

 guts_of_scale( sol , factor );

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

 if( ! v_active_power.empty() ) {
  if( UBS->v_active_power.empty() )
   throw( std::invalid_argument(
	             "UnitBlockSolution::sum: inconsistent active power" ) );

  for( Index i = 0 ; i < f_number_generators ; ++i )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    v_active_power[ i ][ t ] += UBS->v_active_power[ i ][ t ] * multiplier;
  }

 if( ! v_reactive_power.empty() ) {
  if( UBS->v_reactive_power.empty() )
   throw( std::invalid_argument(
	           "UnitBlockSolution::sum: inconsistent reactive power" ) );

  for( Index i = 0 ; i < f_number_generators ; ++i )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    v_reactive_power[ i ][ t ] +=
                                UBS->v_reactive_power[ i ][ t ] * multiplier;
  }

 if( ! v_commitment.empty() ) {
  if( UBS->v_commitment.empty() )
   throw( std::invalid_argument(
		       "UnitBlockSolution::sum: inconsistent commitment" ) );

  for( Index i = 0 ; i < f_number_generators ; ++i )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    v_commitment[ i ][ t ] += UBS->v_commitment[ i ][ t ] * multiplier;
  }

 if( ! v_primary_reserve.empty() ) {
  if( UBS->v_primary_reserve.empty() )
   throw( std::invalid_argument(
		 "UnitBlockSolution::sum: inconsistent primary reserve" ) );

  for( Index i = 0 ; i < f_number_generators ; ++i )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    v_primary_reserve[ i ][ t ] +=
     UBS->v_primary_reserve[ i ][ t ] * multiplier;
  }

 if( ! v_secondary_reserve.empty() ) {
  if( UBS->v_secondary_reserve.empty() )
   throw( std::invalid_argument(
	       "UnitBlockSolution::sum: inconsistent secondary reserve" ) );

  for( Index i = 0 ; i < f_number_generators ; ++i )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    v_secondary_reserve[ i ][ t ] +=
     UBS->v_secondary_reserve[ i ][ t ] * multiplier;
  }
 }  // end( UnitBlockSolution::sum )

/*--------------------------------------------------------------------------*/

UnitBlockSolution * UnitBlockSolution::clone( bool empty ) const
{
 auto * sol = new UnitBlockSolution();

 if( ! empty )
  guts_of_clone( sol );

 return( sol );

 }  // end( UnitBlockSolution::clone )

/*--------------------------------------------------------------------------*/

void UnitBlockSolution::guts_of_clone( UnitBlockSolution * sol ) const
{
 sol->f_time_horizon = f_time_horizon;
 sol->f_number_generators = f_number_generators;

 copy_multi_array( sol->v_active_power , v_active_power );
 copy_multi_array( sol->v_reactive_power , v_reactive_power );
 copy_multi_array( sol->v_commitment , v_commitment );
 copy_multi_array( sol->v_primary_reserve , v_primary_reserve );
 copy_multi_array( sol->v_secondary_reserve , v_secondary_reserve );

 }  // end( UnitBlockSolution::guts_of_clone )

/*--------------------------------------------------------------------------*/

void UnitBlockSolution::guts_of_scale( UnitBlockSolution * sol ,
				       double factor ) const
{
 if( ! v_active_power.empty() )
  for( Index i = 0 ; i < f_number_generators ; ++i )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    sol->v_active_power[ i ][ t ] *= factor;

 if( ! v_reactive_power.empty() )
  for( Index i = 0 ; i < f_number_generators ; ++i )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    sol->v_reactive_power[ i ][ t ] *= factor;

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

 }  // end( UnitBlockSolution::guts_of_scale )

/*--------------------------------------------------------------------------*/
/*---------------------- End File UnitBlock.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
