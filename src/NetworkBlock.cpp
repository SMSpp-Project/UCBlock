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

// register NetworkData to the NetworkData factory

typedef NetworkBlock::NetworkData NetworkData;

SMSpp_insert_in_factory_cpp_0( NetworkData );

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
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void NetworkBlock::deserialize( const netCDF::NcGroup & group ) {
 SMSpp_di_unipi_it::deserialize( group , f_ConstTerm , "ConstantTerm" );
}

/*--------------------------------------------------------------------------*/

void NetworkBlock::NetworkData::deserialize( const netCDF::NcGroup & group ) {
 if( ! deserialize_dim( group , "NumberNodes" , f_number_nodes ) )
  f_number_nodes = 1;
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
 if( f_number_nodes > 1 )
  group.addDim( "NumberNodes" , f_number_nodes );
}

/*--------------------------------------------------------------------------*/

NetworkBlock::NetworkData::NetworkDataFactoryMap &
NetworkBlock::NetworkData::f_factory( void )
{
 static NetworkDataFactoryMap s_factory;
 return( s_factory );
}

/*--------------------------------------------------------------------------*/
/*--------------------- End File NetworkBlock.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
