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

NetworkBlock::NetworkData::NetworkData() {
 f_number_lines = 0;
 f_number_nodes = 0;
 }

/*--------------------------------------------------------------------------*/

void NetworkBlock::NetworkData::deserialize( const netCDF::NcGroup & group )
{
 #ifndef NDEBUG
  static std::vector< std::string > expected_dims = { "NumberNodes" ,
						      "NumberLines" };
  check_dimensions( group, expected_dims, std::cerr );
  static std::vector< std::string > expected_vars = { "StartLine" ,
   "EndLine" , "MinPowerFlow" , "MaxPowerFlow" , "Susceptance" ,
   "NetworkCost" , "NodeName" , "LineName"
   };
  check_variables( group, expected_vars, std::cerr );
 #endif

 if( ! ::deserialize_dim( group, "NumberNodes", f_number_nodes, true ) )
  f_number_nodes = 1;

 if ( f_number_nodes > 1  ) {  // DCNetworkBlock

  ::deserialize_dim( group , "NumberLines", f_number_lines , false );

  ::deserialize( group , "StartLine", f_number_lines , v_start_line , false ,
		 true );

  ::deserialize( group , "EndLine", f_number_lines , v_end_line, false, true );

  ::deserialize( group , "MinPowerFlow" , f_number_lines , v_min_power_flow ,
		 true , true );

  ::deserialize( group , "MaxPowerFlow" , f_number_lines , v_max_power_flow ,
		 true , true );

  ::deserialize( group , "Susceptance" , f_number_lines , v_susceptance ,
		 true , true );

  ::deserialize( group , "NetworkCost" , f_number_lines , v_network_cost ,
		 true , true );
  }

 const auto get_string_array =
  [ &group ]( const std::string & var_name ,
              std::vector< std::string > & v_string ,
              Index size = Inf< Index >() ) {
  v_string.clear();
  auto netcdf_var = group.getVar( var_name );
  if( ! netcdf_var.isNull() ) {
   if( netcdf_var.getDimCount() != 1 )
    throw( std::logic_error( "NetworkData::deserialize: the dimension of "
                             "variable'" + var_name + "' must be 1." ) );

   if( ( size < Inf< Index >() ) &&
       ( netcdf_var.getDim( 0 ).getSize() != size ) )
    throw( std::logic_error( "NetworkData::deserialize: the size of "
                             "variable '" + var_name + "' should be " +
                             std::to_string( size ) + "." ) );

   const auto var_size = netcdf_var.getDim( 0 ).getSize();
   v_string.reserve( var_size );

   // TODO The following implementation should change when netCDF provides a
   // better C++ interface.

   for( Index i = 0 ; i < var_size ; ++i ) {
    char * fname = nullptr;
    netcdf_var.getVar( { i } , { 1 } , & fname );
    v_string.push_back( fname );
    free( fname );
    }
   }
  };

 get_string_array( "NodeName" , v_node_names , f_number_nodes );
 get_string_array( "LineName" , v_line_names );

 }

/*--------------------------------------------------------------------------*/

void NetworkBlock::deserialize( const netCDF::NcGroup & group )
{
 Block::deserialize( group );

 #ifndef NDEBUG
  static std::vector< std::string > expected_dims = { "NumberNodes" };
  check_dimensions( group , expected_dims , std::cerr );
  static std::vector< std::string > expected_vars = { "ActiveDemand" };
  check_variables( group , expected_vars , std::cerr );
 #endif

 const auto NumberNodes = group.getDim( "NumberNodes" );

 if( ! NumberNodes.isNull() ) {
  // A NetworkData has been provided. So, the size of the given vector of
  // active demand must be equal to the number of nodes.
  ::deserialize( group , "ActiveDemand" , NumberNodes.getSize() ,
                 v_active_demand );
  }
 else {
  // A NetworkData has not been provided. However, the active demand may still
  // have been provided.
  auto ActiveDemand = group.getVar( "ActiveDemand" );

  if( ! ActiveDemand.isNull() ) {
   // The active demand has indeed been provided.

   if( ActiveDemand.getDimCount() != 1 )
    // The active demand must be a one-dimensional array.
    throw( std::invalid_argument( "NetworkBlock::deserialize(): ActiveDemand"
                                  " should have one dimension, but it has " +
                                  std::to_string( ActiveDemand.getDimCount() )
                                  ) );

   // Retrieve the number of nodes from the size of the given netCDF variable.
   const auto number_nodes = ActiveDemand.getDim( 0 ).getSize();

   // Resize the vector of active demand.
   v_active_demand.resize( number_nodes );

   // Retrieve the active demand from the netCDF variable.
   ActiveDemand.getVar( v_active_demand.data() );
   }
  }
 }

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
/*---------- METHODS FOR LOADING, PRINTING & SAVING THE NetworkBlock -------*/
/*--------------------------------------------------------------------------*/

void NetworkBlock::NetworkData::serialize( netCDF::NcGroup & group ) const
{
 auto NumberNodes = group.addDim( "NumberNodes" , f_number_nodes );

 if( f_number_nodes > 1 ) {
  auto NL = group.addDim( "NumberLines" , f_number_lines );

  ::serialize( group , "StartLine" , netCDF::NcUint() , NL , v_start_line );

  ::serialize( group , "EndLine" , netCDF::NcUint() , NL , v_end_line );

  ::serialize( group , "MinPowerFlow" , netCDF::NcDouble() , NL ,
	       v_min_power_flow );

  ::serialize( group , "MaxPowerFlow" , netCDF::NcDouble() , NL ,
	       v_max_power_flow );

  ::serialize( group , "Susceptance" , netCDF::NcDouble() , NL ,
	       v_susceptance );

  ::serialize( group , "NetworkCost" , netCDF::NcDouble() ,  NL ,
	       v_network_cost );

  if( ! v_line_names.empty() ) {
   assert( v_line_names.size() == NL.getSize() );
   auto LineName = group.addVar( "LineName" , netCDF::NcString() , NL );
   for( Index i = 0 ; i < v_line_names.size() ; ++i )
    LineName.putVar( { i } , v_line_names[ i ] );
   }
  }

 if( ! v_node_names.empty() ) {
  assert( v_node_names.size() == NumberNodes.getSize() );
  auto NodeName = group.addVar( "NodeName" , netCDF::NcString() , NumberNodes );
  for( Index i = 0 ; i < v_node_names.size() ; ++i )
   NodeName.putVar( { i } , v_node_names[ i ] );
  }
 }

/*--------------------------------------------------------------------------*/

void NetworkBlock::serialize( netCDF::NcGroup & group ) const
{
 Block::serialize( group );

 if( auto network_data = get_NetworkData() )
  // If a NetworkData is present, serialize it.
  network_data->serialize( group );

 if( ! v_active_demand.empty() ) {
  // This NetworkBlock has active demand, so it is serialized.

  auto NumberNodes = group.getDim( "NumberNodes" );

  if( NumberNodes.isNull() ) {
   /* The dimension "NumberNodes" is not present in the group (which means
    * that a NetworkData is not present). However, the number of nodes can
    * still be obtained from the size of the active demand vector. Notice that
    * the name "NumberNodes" is not used for this new dimension, because it
    * would indicate that a NetworkData is present (which is not the
    * case). Therefore, we create an alternative dimension in order to be able
    * to serialize the active demand. */
   NumberNodes = group.addDim( "__NumberNodes__" , v_active_demand.size() );
   }

  // Finally, serialize the active demand.
  ::serialize( group , "ActiveDemand" , netCDF::NcDouble() ,
               NumberNodes , v_active_demand );
  }
 }  // end( NetworkBlock::serialize )

/*--------------------------------------------------------------------------*/
/*--------------------- End File NetworkBlock.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
