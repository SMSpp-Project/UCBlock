/*--------------------------------------------------------------------------*/
/*----------------------- File UCBlock.cpp ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the UCBlock class.
 *
 * \version 0.11
 *
 * \date 10 - 05 - 2019
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
 * Copyright &copy by Antonio Frangioni, Ali Ghezelsoflu, Rafael
 * Durbano Lobato, and Kostas Tavlaridis-Gyparakis
 */

/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <iostream>
#include <vector>
#include "FRowConstraint.h"
#include "LinearFunction.h"
#include "UCBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;


/*@} -----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

//void UCBlock::load( ) { }

/*--------------------------------------------------------------------------*/

template<class T>
void deserialize_dim( const netCDF::NcGroup & group,
                      const std::string & dim_name, T & data ,
                      bool optional = true ) {

  netCDF::NcDim ncDim = group.getDim( dim_name );

  if( ncDim.isNull() ) {
    if( optional )
      return;
    throw( std::invalid_argument( "UCBlock::deserialize: " +
                                  dim_name + " is not present" ) );
  }

  data = ncDim.getSize();
  if( data <= 0 )
    throw( std::invalid_argument( "UCBlock::deserialize: " +
                                  dim_name + " must be positive" ) );
}

/*--------------------------------------------------------------------------*/

template<class T>
void deserialize( const netCDF::NcGroup & group,
                  const std::string & var_name,
                  std::vector<T> & data,
                  const std::vector<size_t> & sizes,
                  bool optional = true ) {

  auto total_size = std::accumulate( begin( sizes ), end( sizes ), 1,
                                     std::multiplies<size_t>() );

  if( total_size == 0 ) {
    data.resize( 0 );
    return;
  }

  auto ncVar = group.getVar( var_name );
  if( ncVar.isNull() ) {
    if( optional ) {
      data.resize( 0 );
      return;
    }
    throw( std::invalid_argument( "UCBlock::deserialize: " +
                                  var_name + " is not present" ) );
  }

  data.resize( total_size );

  std::vector<size_t> start;
  start.assign( sizes.size(), 0 );

  ncVar.getVar( start, sizes, data.data() );
}

/*--------------------------------------------------------------------------*/

void UCBlock::deserialize_sub_blocks( const netCDF::NcGroup & group ) {

  for( auto block : v_Block )
    delete block;

  v_Block.clear();

  deserialize_sub_blocks( group, "UnitBlock_", f_number_units );
  deserialize_sub_blocks( group, "NetworkBlock_", f_time_horizon );
}

/*--------------------------------------------------------------------------*/

void UCBlock::deserialize_sub_blocks( const netCDF::NcGroup & group,
                                      const std::string sub_group_name_prefix,
                                      const int num_sub_blocks ) {

  for( int i = 0; i < num_sub_blocks; ++i ) {

    std::string sub_group_name = sub_group_name_prefix + std::to_string( i );
    auto sub_group = group.getGroup( sub_group_name );

    if( sub_group.isNull() ) {
      throw( std::invalid_argument( "UCBlock::deserialize: " +
                                    sub_group_name + " is not present" ) );
    }

    auto class_name_attribute = sub_group.getAtt( "ClassName" );

    if( class_name_attribute.isNull() ) {
      throw( std::invalid_argument
             ( "UCBlock::deserialize: ClassName attribute "
               "is not present in group " + sub_group_name ) );
    }

    std::string class_name;
    class_name_attribute.getValues( class_name );
    auto sub_block = new_Block( class_name , this );
    sub_block->deserialize( sub_group );
    v_Block.push_back( sub_block );
  }

}

/*--------------------------------------------------------------------------*/

void UCBlock::deserialize( netCDF::NcGroup & group ) {

  ::deserialize_dim( group, "TimeHorizon", f_time_horizon, false );
  ::deserialize_dim( group, "NumberUnits", f_number_units, false );

  // Default values for optional dimensions
  f_number_nodes           = 1;
  f_number_primary_zones   = 0;
  f_number_secondary_zones = 0;
  f_number_inertia_zones   = 0;
  f_number_pollutants      = 0;

  ::deserialize_dim( group, "NumberNodes",          f_number_nodes );
  ::deserialize_dim( group, "NumberPrimaryZones",   f_number_primary_zones );
  ::deserialize_dim( group, "NumberSecondaryZones", f_number_secondary_zones );
  ::deserialize_dim( group, "NumberInertiaZones",   f_number_inertia_zones );
  ::deserialize_dim( group, "NumberPollutants",     f_number_pollutants );

  if( f_number_primary_zones >= 1 ) {
    ::deserialize( group, "PrimaryZones", v_primary_zones,
                   { f_number_nodes } );

    ::deserialize( group, "PrimaryDemand", v_primary_demand,
                   { f_number_primary_zones, f_time_horizon } );
  }

  if( f_number_secondary_zones >= 1 ) {
    ::deserialize( group, "SecondaryZones", v_secondary_zones,
                   { f_number_nodes } );

    ::deserialize( group, "SecondaryDemand", v_secondary_demand,
                   { f_number_secondary_zones, f_time_horizon } );
  }

  if( f_number_inertia_zones >= 1 ) {
    ::deserialize( group, "InertiaZones", v_inertia_zones, { f_number_nodes } );

    ::deserialize( group, "InertiaDemand", v_inertia_demand,
                   { f_number_inertia_zones, f_time_horizon } );
  }

  if( f_number_pollutants >= 1 ) {

    ::deserialize( group, "NumberPollutantZones", v_number_pollutant_zones,
                   { f_number_pollutants } );

    ::deserialize( group, "PollutantZones", v_pollutant_zones,
                   { f_number_pollutants, f_number_nodes } );

    ::deserialize( group, "PollutantDemand", v_pollutant_demand,
                   { f_number_pollutants } );

    ::deserialize( group, "PollutantRho", v_pollutant_rho,
                   { f_time_horizon, f_number_pollutants, f_number_units } );
  }

  if( f_number_nodes > 1 ) {
    ::deserialize( group, "Node", v_node, { f_number_units } );
  }

  deserialize_sub_blocks( group );

}  // end( UCBlock::deserialize )

/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

void UCBlock::generate_abstract_constraints( Configuration *stcc ) {

  if( v_node_injection_constraints.size() != f_time_horizon ) {
    // this should only happen once
    assert( v_node_injection_constraints.size() == 0 );

    v_node_injection_constraints.resize
      ( boost::multi_array<FRowConstraint *, 2>::
        extent_gen()[f_time_horizon][f_number_nodes] );
  }

  // Node injection constraints.

  for( Index t = 0; t < f_time_horizon; ++t ) {

    auto node_injection = v_network_blocks[t]->get_node_injection();

    for( Index node_id = 0; node_id < f_number_nodes; ++node_id ) {

      auto linear_function = new LinearFunction();

      linear_function->add_variable( & node_injection[ node_id ], - 1.0 );

      v_node_injection_constraints[t][node_id]->set_both( 0.0 );
      v_node_injection_constraints[t][node_id]->set_function( linear_function );
    }

    for( Index unit_id = 0; unit_id < f_number_units; ++unit_id ) {
      auto node_id = v_node[ unit_id ];

      auto linear_function = static_cast<LinearFunction *>
        (v_node_injection_constraints[ t ][ node_id ]->get_function());
      linear_function->
        add_variable( get_unit_block( unit_id )->get_power( t ), 1.0 );
    }
  }

  add_static_constraint( v_node_injection_constraints );

  // Primary demand constraints.

  if( v_PrimaryDemand_Const.size() != f_time_horizon ) {
    // this should only happen once
    assert(v_PrimaryDemand_Const.size() == 0);

    v_PrimaryDemand_Const.resize
            (boost::multi_array<FRowConstraint *, 2>::
             extent_gen()[f_time_horizon][f_number_primary_zones]);
  }

//TODO

  // Secondary demand constraints.

  if( v_SecondaryDemand_Const.size() != f_time_horizon ) {
    // this should only happen once
    assert(v_SecondaryDemand_Const.size() == 0);

    v_SecondaryDemand_Const.resize
            (boost::multi_array<FRowConstraint *, 2>::
             extent_gen()[f_time_horizon][f_number_secondary_zones]);
  }

//TODO

  // Inertia demand constraints.

    if( v_InertiaDemand_Const.size() != f_time_horizon ) {
      // this should only happen once
      assert(v_InertiaDemand_Const.size() == 0);

      v_InertiaDemand_Const.resize
              (boost::multi_array<FRowConstraint *, 2>::
               extent_gen()[f_time_horizon][f_number_inertia_zones]);
    }

//TODO

  // Pollutant demand constraints.
    /*
  if( v_PollutantDemand_Const.size() != f_number_pollutants ) {
    // this should only happen once
    assert(v_PollutantDemand_Const.size() == 0);

    v_PollutantDemand_Const.resize
            (boost::multi_array<FRowConstraint *, 2>::
             extent_gen()[f_number_pollutants][f_number_pollutant_zones]);
  }
    */

//TODO

}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------ METHODS FOR LOADING, PRINTING & SAVING THE UCBlock ----------*/
/*--------------------------------------------------------------------------*/

template<class T>
void serialize( netCDF::NcGroup & group, const std::string & var_name,
                const netCDF::NcType & ncType,
                const std::vector<netCDF::NcDim> & ncDim,
                const std::vector<T> & data ) {

  std::vector<size_t> start;
  start.assign( ncDim.size(), 0 );

  std::vector<size_t> sizes;
  sizes.resize( ncDim.size() );
  for( size_t i = 0; i < sizes.size(); ++i )
    sizes[i] = ncDim[i].getSize();

  auto total_size = std::accumulate( begin( sizes ), end( sizes ), 1,
                                     std::multiplies<size_t>() );

  if( total_size == 0 )
    return;

  group.addVar( var_name , ncType , ncDim )
    .putVar( start , sizes , data.data() );
}

/*--------------------------------------------------------------------------*/

void UCBlock::serialize( netCDF::NcGroup & group ) const {

  group.putAtt( "type" , "UCBlock" );

  auto dim_time_horizon = group.addDim( "TimeHorizon", f_time_horizon );
  auto dim_number_units = group.addDim( "NumberUnits", f_number_units );
  auto dim_number_nodes = group.addDim( "NumberNodes", f_number_nodes );

  auto dim_number_primary_zones =
    group.addDim( "NumberPrimaryZones", f_number_primary_zones );
  auto dim_number_secondary_zones =
    group.addDim( "NumberSecondaryZones", f_number_secondary_zones );
  auto dim_number_inertia_zones =
    group.addDim( "NumberInertiaZones", f_number_inertia_zones );
  auto dim_number_pollutants =
    group.addDim( "NumberPollutants", f_number_pollutants );

  if( f_number_primary_zones >= 1 ) {
    ::serialize( group, "PrimaryZones", netCDF::NcUint64(),
                 { dim_number_nodes },
                 v_primary_zones );

    ::serialize( group, "PrimaryDemand", netCDF::NcDouble(),
                 { dim_number_primary_zones, dim_time_horizon },
                 v_primary_demand );
  }

  if( f_number_secondary_zones >= 1 ) {
    ::serialize( group, "SecondaryZones", netCDF::NcUint64(),
                 { dim_number_nodes },
                 v_secondary_zones );

    ::serialize( group, "SecondaryDemand", netCDF::NcDouble(),
                 { dim_number_secondary_zones, dim_time_horizon },
                 v_secondary_demand );
  }

  if( f_number_inertia_zones >= 1 ) {
    ::serialize( group, "InertiaZones", netCDF::NcUint64(),
                 { dim_number_nodes }, v_inertia_zones );

    ::serialize( group, "InertiaDemand", netCDF::NcDouble(),
                 { dim_number_inertia_zones, dim_time_horizon } ,
                 v_inertia_demand );
  }

  if( f_number_pollutants >= 1 ) {

    ::serialize( group, "NumberPollutantZones", netCDF::NcUint64(),
                 { dim_number_pollutants } ,
                 v_number_pollutant_zones );

    ::serialize( group, "PollutantZones", netCDF::NcUint64(),
                 { dim_number_pollutants, dim_number_nodes },
                 v_pollutant_zones );

    ::serialize( group, "PollutantDemand", netCDF::NcDouble(),
                 { dim_number_pollutants }, v_pollutant_demand );

    ::serialize( group, "PollutantRho", netCDF::NcDouble(),
                 { dim_time_horizon, dim_number_pollutants, dim_number_units },
                 v_pollutant_rho );
  }

  if( f_number_nodes > 1 ) {
    ::serialize( group, "Node", netCDF::NcUint64(),
                 { dim_number_units }, v_node );
  }

  // Serialize sub-blocks

  for( Index i = 0; i < f_number_units; ++i ) {
    auto sub_block = get_unit_block( i );
    auto sub_group = group.addGroup( "UnitBlock_" + std::to_string( i ) );
    sub_block->serialize( sub_group );
  }

  for( Index t = 0; t < f_time_horizon; ++t ) {
    auto sub_block = get_network_block( t );
    auto sub_group = group.addGroup( "NetworkBlock_" + std::to_string( t ) );
    sub_block->serialize( sub_group );
  }

}  // end( UCBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------------------ End File UCBlock.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
