/*--------------------------------------------------------------------------*/
/*--------------------- File NetworkBlock.cpp ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the NetworkBlock class.
 *
 * \version 0.11
 *
 * \date 23 - 06 - 2019
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

void NetworkBlock::NetworkData::deserialize( netCDF::NcGroup & group ) {


 ::deserialize_dim( group, "NumberNodes", f_number_nodes, true );

 if ( f_number_nodes > 1  ) {  // DCNetworkBlock

  ::deserialize_dim( group, "NumberLines", f_number_lines, false );

  ::deserialize( group, "StartLine", f_number_lines, v_start_line, false, true );

  ::deserialize( group, "EndLine", f_number_lines, v_end_line, false, true );

  ::deserialize( group, "MinPowerFlow", f_number_lines, v_min_power_flow, true, true );

  ::deserialize( group, "MaxPowerFlow", f_number_lines, v_max_power_flow, true, true );

  ::deserialize( group, "Susceptance", f_number_lines, v_susceptance, true, true );
 }

}

/*--------------------------------------------------------------------------*/
void NetworkBlock::deserialize( netCDF::NcGroup & group ) {

 auto network_data = new NetworkBlock::NetworkData();
 network_data->deserialize( group );

 auto dim_number_nodes = group.getDim( "NumberNodes" );

 if( !dim_number_nodes.isNull() )
  ::deserialize( group, "ActiveDemand", dim_number_nodes.getSize(),
                 v_active_demand );

}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*---------- METHODS FOR LOADING, PRINTING & SAVING THE NetworkBlock -------*/
/*--------------------------------------------------------------------------*/

void NetworkBlock::NetworkData::serialize( netCDF::NcGroup & group ) const {

 auto dim_number_nodes = group.addDim( "NumberNodes", f_number_nodes );

 if( f_number_nodes > 1 ) { // DCNetworkBlock
  auto dim_number_lines = group.addDim( "NumberLines", f_number_lines );

  ::serialize( group, "StartLine", netCDF::NcUint64(),
               { dim_number_lines }, v_start_line );

  ::serialize( group, "EndLine", netCDF::NcUint64(),
               { dim_number_lines }, v_end_line );

  ::serialize( group, "MinPowerFlow", netCDF::NcDouble(),
               { dim_number_lines }, v_min_power_flow );

  ::serialize( group, "MaxPowerFlow", netCDF::NcDouble(),
               { dim_number_lines }, v_max_power_flow );

  ::serialize( group, "Susceptance", netCDF::NcDouble(),
               { dim_number_lines }, v_susceptance );
 }

}

/*--------------------------------------------------------------------------*/
void NetworkBlock::serialize( netCDF::NcGroup & group ) const {

 group.putAtt( "type", name() );
 auto dim_number_nodes = group.getDim( "NumberNodes" );

 if( !dim_number_nodes.isNull() )
  ::serialize( group, "ActiveDemand", netCDF::NcDouble(),
               { dim_number_nodes }, v_active_demand );

}

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

void
NetworkBlock::set_active_demand( std::vector< double >::const_iterator values,
                                 Block::Subset && subset,
                                 const bool ordered,
                                 c_ModParam issuePMod,
                                 c_ModParam issueAMod ) {
 if( subset.empty() ) {
  return;
 }

 if( v_active_demand.empty() ) {
  if( std::all_of( values,
                   values + subset.size(),
                   []( double cst ) {
                    return ( cst == 0 );
                   } ) ) {
   return;
  }

  Index max_index = *max_element( std::begin( subset ), std::end( subset ) );
  v_active_demand.assign( max_index, 0 );
 }

 // If nothing changes, return
 bool identical = true;
 for( auto i : subset ) {
  if( i >= v_active_demand.size() ) {
   throw ( std::invalid_argument( "invalid value in subset" ) );
  }
  if( v_active_demand[ i ] != *( values++ ) ) {
   identical = false;
  }
 }
 if( identical ) {
  return;
 }

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  for( auto i : subset ) {
   v_active_demand[ i ] = *( values++ );
  }

  if( not_dry_run( issueAMod ) && AR & HasCst ) {
   // Change the abstract representation

   // TODO: Change the maxpower values where they are used!
  }

  if( issue_pmod( issuePMod ) ) {
   // Issue a Physical Modification
   if( !ordered ) {
    std::sort( subset.begin(), subset.end() );
   }

   Block::add_Modification(
    std::make_shared< NetworkBlockSbstMod >( this,
                                             NetworkBlockMod::eSetActD,
                                             std::move( subset ) ),
    Observer::par2chnl( issuePMod ) );
  }
 }
}

void
NetworkBlock::set_active_demand( std::vector< double >::const_iterator values,
                                 Block::Range rng,
                                 c_ModParam issuePMod,
                                 c_ModParam issueAMod ) {

 rng.second = std::min( rng.second, get_number_nodes() );
 if( rng.second <= rng.first ) {
  return;
 }

 if( v_active_demand.empty() ) {
  if( std::all_of( values,
                   values + ( rng.second - rng.first ),
                   []( double cst ) {
                    return ( cst == 0 );
                   } ) ) {
   return;
  }

  Index max_index = rng.second;
  v_active_demand.assign( max_index, 0 );
 }

 // If nothing changes, return
 if( std::equal( values,
                 values + ( rng.second - rng.first ),
                 v_active_demand.begin() + rng.first ) ) {
  return;
 }

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  std::copy( values,
             values + ( rng.second - rng.first ),
             v_active_demand.begin() + rng.first );

  if( AR & HasCst ) {
   // Change the abstract representation

   // TODO: Change the maxpower values where they are used!
  }

  if( issue_pmod( issuePMod ) ) {
   Block::add_Modification(
    std::make_shared< NetworkBlockRngdMod >( this,
                                             NetworkBlockMod::eSetActD,
                                             rng ),
    Observer::par2chnl( issuePMod ) );
  }
 }
}
// end( NetworkBlock::serialize )

/*--------------------------------------------------------------------------*/
/*--------------------- End File NetworkBlock.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
