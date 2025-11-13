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
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Ali Ghezelsoflu,
 *                      Rafael Durbano Lobato, Kostas Tavlaridis-Gyparakis,
 *                      Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <map>

#include "NetworkBlock.h"

#include "FRowConstraint.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register NetworkData to the NetworkData factory

typedef NetworkBlock::NetworkData NetworkData;

SMSpp_insert_in_factory_cpp_0( NetworkData );

// register NetworkBlockSolution to the Solution factory
SMSpp_insert_in_factory_cpp_0( NetworkBlockSolution );

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
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void NetworkBlock::generate_abstract_variables( Configuration * stvv )
{
 const auto number_nodes = get_number_nodes();
 const auto number_intervals = get_number_intervals();

 if( number_nodes > 1 ) {
  // the node injection variables
  v_node_injection.resize( boost::extents[ number_intervals ][ number_nodes ] );
  for( Index t = 0 ; t < number_intervals ; ++t )
   for( Index node_id = 0 ; node_id < number_nodes ; ++node_id )
    v_node_injection[ t ][ node_id ].set_type( ColVariable::kContinuous );
  add_static_variable( v_node_injection , "s_network" );
 }
}  // end( NetworkBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void NetworkBlock::generate_abstract_constraints( Configuration * stcc )
{
 if( constraints_generated() )  // constraints have already been generated
  return;                       // nothing to do

 const auto number_nodes = get_number_nodes();
 const auto number_intervals = get_number_intervals();

 // node injection bound constraints

 node_injection_bounds_const.resize(
  boost::multi_array< FRowConstraint , 2 >::extent_gen()
  [ number_nodes ][ number_intervals ] );

 for( Index i = 0 ; i < number_intervals ; ++i )

  for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {

   node_injection_bounds_const[ node_id ][ i ].set_lhs(
    v_MinNodeInjection[ i ][ node_id ] );
   node_injection_bounds_const[ node_id ][ i ].set_rhs(
    v_MaxNodeInjection[ i ][ node_id ] );
   node_injection_bounds_const[ node_id ][ i ].set_variable(
    &v_node_injection[ i ][ node_id ] );
  }

 add_static_constraint( node_injection_bounds_const ,
                        "Node_Injection_Bound_Const_Network" );

 }  // end( NetworkBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/

Solution * NetworkBlock::get_Solution( Configuration * csolc , bool emptys )
{
 Index wsol = 1;
 if( ( ! csolc ) && f_BlockConfig )
  csolc = f_BlockConfig->f_solution_Configuration;

 if( auto config = dynamic_cast< SimpleConfiguration< int > * >( csolc ) )
  wsol = config->f_value;

 auto sol = new_Solution();

 using mad2 = boost::multi_array< double , 2 >;

 if( wsol & 1 )
  sol->v_node_injection.resize(
        mad2::extent_gen()[ get_number_intervals() ][ get_number_nodes() ] );

 if( ! emptys )
  sol->read( this );

 return( sol );
 }

/*--------------------------------------------------------------------------*/

NetworkBlockSolution * NetworkBlock::new_Solution( void ) const {
  return( new NetworkBlockSolution() );
  }

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void NetworkBlock::deserialize( const netCDF::NcGroup & group ) {
 ::deserialize( group , f_ConstTerm , "ConstantTerm" );
 }

/*--------------------------------------------------------------------------*/

void NetworkBlock::NetworkData::deserialize( const netCDF::NcGroup & group )
{
 if( ! deserialize_dim( group , "NumberNodes" , f_number_nodes ) )
  f_number_nodes = 1;

 if( f_number_nodes > 1 )
  ::deserialize( group , "NodeName" , f_number_nodes , v_node_names );
 }

/*--------------------------------------------------------------------------*/
/*---------- METHODS FOR LOADING, PRINTING & SAVING THE NetworkBlock -------*/
/*--------------------------------------------------------------------------*/

void NetworkBlock::serialize( netCDF::NcGroup& group ) const {
 if( f_ConstTerm != 0 )
  ::serialize( group , "ConstantTerm" , netCDF::NcDouble() , f_ConstTerm );
 }

/*--------------------------------------------------------------------------*/

void NetworkBlock::NetworkData::serialize( netCDF::NcGroup& group ) const {
 if( f_number_nodes > 1 ) {
  auto NumberNodes = group.getDim( "NumberNodes" );

  if( ! v_node_names.empty() ) {
   auto NodeName = group.addVar( "NodeName" , netCDF::NcString() ,
				 NumberNodes );
   for( Index i = 0 ; i < v_node_names.size() ; ++i )
    NodeName.putVar( { i } , v_node_names[ i ] );
   }
  }
 }

/*--------------------------------------------------------------------------*/

NetworkBlock::NetworkData::NetworkDataFactoryMap &
NetworkBlock::NetworkData::f_factory( void )
{
 static NetworkDataFactoryMap s_factory;
 return( s_factory );
 }

/*--------------------------------------------------------------------------*/
/*------------------- METHODS OF NetworkBlockSolution ----------------------*/
/*--------------------------------------------------------------------------*/

void NetworkBlockSolution::deserialize( const netCDF::NcGroup & group )
{
 // "NumberNodes" is mandatory- - - - - - - - - - - - - - - - - - - - - - - -
 deserialize_dim( group , "NumberNodes" , f_number_nodes , false );

 // "NumberInstants" is optional- - - - - - - - - - - - - - - - - - - - - - -
 if( ! deserialize_dim( group , "NumberInstants" , f_number_intervals , true )
     )
  f_number_intervals = 1;

 // deserialize the Node Injection - - - - - - - - - - - - - - - - - - - - -
 if( ! ::deserialize< double , 2 >( group , "NodeInjection" ,
                                    { f_number_nodes , f_number_intervals } ,
                                    v_node_injection , true ) ) {
  std::vector< boost::multi_array< double , 2 >::index > sizes( 2 , 0 );
  v_node_injection.resize( sizes );
  }
 }  // end( NetworkBlockSolution::deserialize( NcGroup & )

/*--------------------------------------------------------------------------*/

void NetworkBlockSolution::deserialize( const netCDF::NcGroup & group ,
					size_t idx )
{
 // "NumberNodes" is mandatory- - - - - - - - - - - - - - - - - - - - - - - -
 deserialize_dim( group , "NumberNodes" , f_number_nodes , false );

 // "NumberNetworks" is mandatory - - - - - - - - - - - - - - - - - - - - - -
 int nnw;
 deserialize_dim( group , "NumberNetworks" , nnw , false );
 if( idx >= nnw )
  throw( std::invalid_argument(
	                "NetworkBlockSolution::deserialize: invalid idx" ) );

 // "TotalNumberInstants" is optional - - - - - - - - - - - - - - - - - - - -
 int tni;
 size_t start;
 if( deserialize_dim( group , "TotalNumberInstants" , tni , true ) ) {
  // ... but if it is provided, then "EndInstant" is mandatory
  auto EI = group.getVar( "EndInstant" );
  if( EI.isNull() )
   throw( std::invalid_argument(
	          "NetworkBlockSolution::serialize: missing EndInstant" ) );
  start = 0;
  if( idx > 0 ) {
   std::vector< size_t > vidx = { idx - 1 };
   EI.getVar( vidx , & start );
   }

  int ei;
  std::vector< size_t > vidx = { idx };
  EI.getVar( vidx , & ei );
  f_number_intervals = ei - start;
  }
 else {  // each NetworkBlockSolution covers one instant
  f_number_intervals = 1;
  start = idx;
  tni = nnw;
  }

 // deserialize the Node Injection - - - - - - - - - - - - - - - - - - - - -
 using mad2i = boost::multi_array< double , 2 >::index;

 auto ncVar = group.getVar( "NodeInjection" );
 if( ncVar.isNull() ) {
  std::vector< mad2i > sizes = { 0 , 0 };
  v_node_injection.resize( sizes );
  return;
  }

 std::vector< mad2i > sizes = { f_number_intervals , f_number_nodes };
 v_node_injection.resize( sizes );

 std::vector< size_t > strt = { start , 0 };
 std::vector< size_t > cnt = { f_number_intervals , f_number_nodes };
 ncVar.getVar( strt , cnt , v_node_injection.data() );

 }  // end( NetworkBlockSolution::deserialize( NcGroup & , size_t ) )

/*--------------------------------------------------------------------------*/

void NetworkBlockSolution::deserialize( const netCDF::NcGroup & group ,
		  std::vector< NetworkBlockSolution * > & sols , size_t idx )
{
 // get the type of all NetworkBlockSolution
 std::string tmp;
 auto gtype = group.getAtt( "NBSType" );
 if( gtype.isNull() )
  throw( std::invalid_argument( "NetworkBlockSolution::deserialize: NBSType "
				"not present" ) );
 gtype.getValues( tmp );

 // now deserialize them all
 for( auto & soli : sols ) {
  auto result = new_Solution( tmp );
  auto NSi = dynamic_cast< NetworkBlockSolution * >( result );
  if( ! NSi )
   throw( std::invalid_argument( "NetworkBlockSolution::deserialize: "
				 "invalid NetworkBlockSolution " +
				 std::to_string( idx ) ) );
  NSi->deserialize( group , idx++ );
  soli = NSi;
  }
 }  // end( NetworkBlockSolution::deserialize( netCDF::NcGroup & group ,
    //		           std::vector< NetworkBlockSolution * > , size_t ) )

/*--------------------------------------------------------------------------*/

void NetworkBlockSolution::read( const Block * block )
{
 auto NB = dynamic_cast< const NetworkBlock * >( block );
 if( ! NB )
  throw( std::invalid_argument(
	       "NetworkBlockSolution::read: block is not a NetworkBlock" ) );

 f_number_nodes = NB->get_number_nodes();
 f_number_intervals = NB->get_number_intervals();

 if( ! v_node_injection.empty() )
  // read the node injection variables - - - - - - - - - - - - - - - - - - -
  for( Index t = 0 ; t < f_number_intervals ; ++t ) {
   auto NIt = NB->get_const_node_injection( t );
   for( Index i = 0 ; i < f_number_nodes ; ++i )
    v_node_injection[ t ][ i ] = NIt[ i ].get_value();
   }

 }  // end( NetworkBlockSolution::read )

/*--------------------------------------------------------------------------*/

void NetworkBlockSolution::write( Block * block )
{
 auto NB = dynamic_cast< NetworkBlock * >( block );
 if( ! NB )
  throw( std::invalid_argument(
	      "NetworkBlockSolution::write: block is not a NetworkBlock" ) );

 if( f_number_nodes != NB->get_number_nodes() )
  throw( std::invalid_argument(
		  "NetworkBlockSolution::write: inconsistent node number" ) );

 if( f_number_intervals != NB->get_number_intervals() )
  throw( std::invalid_argument(
	      "NetworkBlockSolution::write: inconsistent instants number" ) );

 if( ! v_node_injection.empty() )
  // write the node injection variables- - - - - - - - - - - - - - - - - - -
  for( Index t = 0 ; t < f_number_intervals ; ++t ) {
   auto NIt = NB->get_node_injection( t );
   for( Index i = 0 ; i < f_number_nodes ; ++i )
    NIt[ i ].set_value( v_node_injection[ t ][ i ] );
   }

 }  // end( NetworkBlockSolution::write )

/*--------------------------------------------------------------------------*/

void NetworkBlockSolution::serialize( netCDF::NcGroup & group ) const
{
 Solution::serialize( group );

 // "NumberNodes" is mandatory- - - - - - - - - - - - - - - - - - - - - - - -
 auto nn = group.addDim( "NumberNodes" , f_number_nodes );

 // "NumberInstants" is optional- - - - - - - - - - - - - - - - - - - - - - -
 netCDF::NcDim ni;
 if( f_number_intervals > 1 )
  ni = group.addDim( "NumberInstants" , f_number_intervals );

 // serialize the Node Injection- - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_node_injection.empty() ) {
  if( ni.isNull() ) {
   std::vector< double > tmp_injection( f_number_nodes );
   for( Index i = 0 ; i < f_number_nodes ; ++i )
    tmp_injection[ i ] = v_node_injection[ 0 ][ i ];
   ::serialize< double >( group , "NodeInjection" , netCDF::NcDouble() ,
                          nn , tmp_injection );
   }
  else
   ::serialize< double , 2 >( group , "NodeInjection" , netCDF::NcDouble() ,
                              { ni , nn } , v_node_injection );
  }
 }  // end( NetworkBlockSolution::serialize( NcGroup & ) )

/*--------------------------------------------------------------------------*/

void NetworkBlockSolution::serialize( netCDF::NcGroup & group , size_t idx )
 const
{
 // "NumberNetworks" is mandatory, and it must be there already - - - - - - -
 auto nnw = group.getDim( "NumberNetworks" );
 if( nnw.isNull() )
  throw( std::invalid_argument(
	       "NetworkBlockSolution::serialize: missing NumberNetworks" ) );

 if( idx >= nnw.getSize() )
  throw( std::invalid_argument(
	                  "NetworkBlockSolution::serialize: invalid idx" ) );

 // "NumberNodes" is mandatory, and if it is not there already it must be
 // created right now - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 auto nn = group.getDim( "NumberNodes" );
 if( nn.isNull() ) {
  if( idx == 0 )  // first call, initialize it
   nn = group.addDim( "NumberNodes" , f_number_nodes );
  else
   throw( std::invalid_argument(
	          "NetworkBlockSolution::serialize: missing NumberNodes" ) );
  }

 size_t start = idx;
 // "TotalNumberInstants" is optional, but it must be there already - - - - -
 auto tni = group.getDim( "TotalNumberInstants" );
 if( tni.isNull() )
  tni = nnw;
 else {
  if( tni.getSize() < nnw.getSize() )
   throw( std::invalid_argument( "NetworkBlockSolution::serialize: less "
				 "instants than networks" ) );

  // if TotalNumberInstants > NumberNetworks, EndInstant is mandatory
  auto EI = group.getVar( "EndInstant" );
  if( EI.isNull() ) {  // if it's not there
   if( idx == 0 )      // and idx == 0, create it now
    EI = group.addVar( "EndInstant" , netCDF::NcInt() , { nnw } );
   else
    throw( std::invalid_argument( "NetworkBlockSolution::serialize: "
				  "missing EndInstant" ) );
   }

  start = 0;
  if( idx > 0 ) {
   std::vector< size_t > vidx = { idx - 1 };
   EI.getVar( vidx , & start );
   }

  // now write EndInstant[ idx ]
  std::vector< size_t > vidx = { idx };
  EI.putVar( vidx , static_cast< int >( start + f_number_intervals ) );
  }

 // now serialize the data structures - - - - - - - - - - - - - - - - - - - -

 if( ! v_node_injection.empty() ) {  // node injections have to be serialised
  netCDF::NcVar NI;  // NodeInjection

  if( idx == 0 )  // first call, have to initialize
   NI = group.addVar( "NodeInjection" , netCDF::NcDouble() , { tni , nn } );
  else            // subsequent cqll, read what is supposedly already there
   NI = group.getVar( "NodeInjection" );

  std::vector< size_t > strt = { start , 0 };
  std::vector< size_t > cnt = { f_number_intervals , f_number_nodes };
  NI.putVar( strt , cnt , v_node_injection.data() );
  }
 }  // end( NetworkBlockSolution::serialize( NcGroup & , size_t ) )

/*--------------------------------------------------------------------------*/

void NetworkBlockSolution::serialize( netCDF::NcGroup & group ,
	    const std::vector< NetworkBlockSolution * > & sols , size_t idx )
{
 // ensure NumberNetworks is there
 auto nnw = group.getDim( "NumberNetworks" );
 if( nnw.isNull() )
  nnw = group.addDim( "NumberNetworks" , sols.size() );

 if( idx >= nnw.getSize() )
  throw( std::invalid_argument( "NetworkBlockSolution::serialize: invalid "
				"idx" ) );

 // ensure TotalNumberInstants is there, if needed
 Index ni = 0;
 for( auto nbsi : sols )
  ni += nbsi->get_number_intervals();

 if( ni < nnw.getSize() )
  throw( std::invalid_argument( "NetworkBlockSolution::serialize: less "
				"instants than networks" ) );

 auto tni = group.getDim( "TotalNumberInstants" );
 if( tni.isNull() && ( ni > nnw.getSize() ) )
  group.addDim( "TotalNumberInstants" , ni );

 // ensure type is there
 auto ty = group.getAtt( "NBSType" );
 if( ty.isNull() )
  group.putAtt( "NBSType" , (sols.front())->classname() );

 // finally, serialize all NetworkBlockSolution
 for( auto nbsi : sols )
  nbsi->serialize( group , idx++ );

 }  // end( NetworkBlockSolution::serialize( NcGroup & ,
    //                   std::vector< NetworkBlockSolution * > & , size_t ) )

/*--------------------------------------------------------------------------*/

NetworkBlockSolution * NetworkBlockSolution::scale( double factor ) const
{
 auto sol = clone();

 if( factor == 1 )
  return( sol );

 if( ! v_node_injection.empty() )
  for( Index t = 0 ; t < f_number_intervals ; ++t )
   for( Index i = 0 ; i < f_number_nodes ; ++i )
    sol->v_node_injection[ t ][ i ] *= factor;

 return( sol );

 }  // end( NetworkBlockSolution::scale )

/*--------------------------------------------------------------------------*/

void NetworkBlockSolution::sum( const Solution * solution ,
				double multiplier )
{
 auto NBS = dynamic_cast< const NetworkBlockSolution * >( solution );
 if( ! NBS )
  throw( std::invalid_argument(
      "NetworkBlockSolution::sum: solution is not a NetworkBlockSolution" ) );

 if( f_number_nodes != NBS->f_number_nodes )
  throw( std::invalid_argument(
		    "NetworkBlockSolution::sum: inconsistent node number" ) );

 if( f_number_intervals != NBS->f_number_intervals )
  throw( std::invalid_argument(
	        "NetworkBlockSolution::sum: inconsistent instants number" ) );

 if( ! v_node_injection.empty() )
  for( Index t = 0 ; t < f_number_intervals ; ++t )
   for( Index i = 0 ; i < f_number_nodes ; ++i )
    v_node_injection[ t ][ i ] +=
     NBS->v_node_injection[ t ][ i ] * multiplier;

 }  // end( NetworkBlockSolution::sum )

/*--------------------------------------------------------------------------*/

NetworkBlockSolution * NetworkBlockSolution::clone( bool empty ) const
{
 auto sol = new NetworkBlockSolution();

 if( ! empty )
  guts_of_clone( sol );

 return( sol );

 }  // end( NetworkBlockSolution::clone )

/*--------------------------------------------------------------------------*/

void NetworkBlockSolution::guts_of_clone( NetworkBlockSolution * sol ) const
{
 sol->f_number_nodes = f_number_nodes;
 sol->f_number_intervals = f_number_intervals;

 copy_multi_array( sol->v_node_injection , v_node_injection );

 }  // end( NetworkBlockSolution::clone )

/*--------------------------------------------------------------------------*/
/*--------------------- End File NetworkBlock.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
