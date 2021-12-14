/*--------------------------------------------------------------------------*/
/*--------------------- File NetworkBlock.cpp ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the NetworkBlock class.
 *
 * \version 0.11
 *
 * \date 14 - 12 - 2021
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
/*--------------------------------------------------------------------------*/

NetworkBlock::NetworkData::NetworkData() {
 f_number_lines = 0;
 f_number_nodes = 0;
}

void NetworkBlock::NetworkData::deserialize( const netCDF::NcGroup & group ) {

#ifndef NDEBUG
 std::vector< std::string > expected_dims = { "NumberNodes" , "NumberLines"};
 check_dimensions( group, expected_dims, std::cerr );
 std::vector< std::string > expected_vars = { "StartLine" , "EndLine" ,
  "MinPowerFlow" , "MaxPowerFlow" , "Susceptance" , "NetworkCost"};
 check_variables( group, expected_vars, std::cerr );
#endif

 if( ! ::deserialize_dim( group, "NumberNodes", f_number_nodes, true ) )
  f_number_nodes = 1;

 if ( f_number_nodes > 1  ) {  // DCNetworkBlock

  ::deserialize_dim( group, "NumberLines", f_number_lines, false );

  ::deserialize( group, "StartLine", f_number_lines, v_start_line, false, true );

  ::deserialize( group, "EndLine", f_number_lines, v_end_line, false, true );

  ::deserialize( group, "MinPowerFlow", f_number_lines, v_min_power_flow, true, true );

  ::deserialize( group, "MaxPowerFlow", f_number_lines, v_max_power_flow, true, true );

  ::deserialize( group, "Susceptance", f_number_lines, v_susceptance, true, true );

  ::deserialize( group, "NetworkCost", f_number_lines, v_network_cost, true, true );

 }

}

/*--------------------------------------------------------------------------*/

void NetworkBlock::deserialize( const netCDF::NcGroup & group ) {

#ifndef NDEBUG
 std::vector< std::string > expected_dims = { "NumberNodes" };
 check_dimensions( group , expected_dims , std::cerr );
 std::vector< std::string > expected_vars = { "ActiveDemand" };
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
                                  std::to_string( ActiveDemand.getDimCount() ) +
                                  "." ) );

   // Retrieve the number of nodes from the size of the given netCDF variable.
   const auto number_nodes = ActiveDemand.getDim( 0 ).getSize();

   // Resize the vector of active demand.
   v_active_demand.resize( number_nodes );

   // Retrieve the active demand from the netCDF variable.
   ActiveDemand.getVar( v_active_demand.data() );
  }
 }

 Block::deserialize( group );
}

/*--------------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/

Solution * NetworkBlock::get_Solution( Configuration * csolc, bool emptys )
{
 auto config = dynamic_cast< SimpleConfiguration< int > * >( csolc );

 if( ( ! config ) && f_BlockConfig )
  config = dynamic_cast<SimpleConfiguration< int > *>(
          f_BlockConfig->f_solution_Configuration );

 auto solution_type = config ? config->f_value : 0;

 Solution * sol = nullptr;
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
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------- METHODS FOR LOADING, PRINTING & SAVING THE NetworkBlock -------*/
/*--------------------------------------------------------------------------*/

void NetworkBlock::NetworkData::serialize( netCDF::NcGroup & group ) const {

 auto NumberNodes = group.addDim( "NumberNodes" , f_number_nodes );

 if( f_number_nodes > 1 ) { // DCNetworkBlock
  auto NumberLines = group.addDim( "NumberLines" , f_number_lines );

  ::serialize( group , "StartLine" , netCDF::NcUint() ,
               NumberLines , v_start_line );

  ::serialize( group , "EndLine" , netCDF::NcUint() ,
               NumberLines , v_end_line );

  ::serialize( group , "MinPowerFlow" , netCDF::NcDouble() ,
               NumberLines , v_min_power_flow );

  ::serialize( group , "MaxPowerFlow" , netCDF::NcDouble() ,
               NumberLines , v_max_power_flow );

  ::serialize( group , "Susceptance" , netCDF::NcDouble() ,
               NumberLines , v_susceptance );

  ::serialize( group , "NetworkCost" , netCDF::NcDouble() ,
               NumberLines , v_network_cost );
 }
}

/*--------------------------------------------------------------------------*/

void NetworkBlock::serialize( netCDF::NcGroup & group ) const {

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
} // end( NetworkBlock::serialize )

/*--------------------------------------------------------------------------*/
/*--------------------- End File NetworkBlock.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
