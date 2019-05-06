/*--------------------------------------------------------------------------*/
/*----------------------- File UCBlock.cpp ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the UCBlock class.
 *
 * \version 0.11
 *
 * \date 03 - 05 - 2019
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
#include "NetworkBlock.h"
#include "NetworkNode.h"
#include "UCBlock.h"
#include "UnitBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;


/*@} -----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */



//void UCBlock::load( )
//{


//}


/*--------------------------------------------------------------------------*/

template<class T>
void deserialize_dim( const netCDF::NcGroup & group,
                      const std::string & dim_name, T & data ) {

  netCDF::NcDim ncDim = group.getDim( dim_name );

  if( ncDim.isNull() )
    throw( std::invalid_argument( "UCBlock::deserialize: " +
                                  dim_name + " is not present" ) );

  data = ncDim.getSize();
  if( data <= 0 )
    throw( std::invalid_argument( "UCBlock::deserialize: " +
                                  dim_name + " must be positive" ) );
}

/*--------------------------------------------------------------------------*/

template<class T>
void deserialize( const netCDF::NcGroup & group,
                  const std::string & var_name,
                  std::vector<T> & data ,
                  const typename std::vector<T>::size_type & size ) {

  auto ncVar = group.getVar( var_name );
  if( ncVar.isNull() )
    throw( std::invalid_argument( "UCBlock::deserialize: " +
                                  var_name + " is not present" ) );

  data.resize( size );
  ncVar.getVar( { 0 }, { size }, data.data() );
}

/*--------------------------------------------------------------------------*/

void UCBlock::deserialize( netCDF::NcGroup & group ) {

  ::deserialize_dim( group, "TimeHorizon",    f_time_horizon );
  ::deserialize_dim( group, "NumberUnits",    f_number_units );
  ::deserialize_dim( group, "NumberNodes",    f_number_nodes );
  ::deserialize_dim( group, "PrimaryZones",   f_number_primary_zones );
  ::deserialize_dim( group, "SecondaryZones", f_number_secondary_zones );
  ::deserialize_dim( group, "InertiaZones",   f_number_inertia_zones );
  ::deserialize_dim( group, "EmissionZones",  f_number_emission_zones );
  ::deserialize_dim( group, "PollutantSet",   f_number_pollutants );

  ::deserialize( group, "PrimaryDemand",   f_primary_demand,   f_number_primary_zones );
  ::deserialize( group, "SecondaryDemand", f_secondary_demand, f_number_secondary_zones );
  ::deserialize( group, "InertiaDemand",   f_inertia_demand,   f_number_inertia_zones );
  ::deserialize( group, "PollutantDemand", f_pollutant_demand, f_number_pollutants );

  // TODO PollutantRho

}  // end( UCBlock::deserialize )

/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

void UCBlock::generate_abstract_constraints( Configuration *stcc ) {

  auto num_nodes = f_network.get_num_nodes();

  if( v_node_injection_constraints.size() != f_time_horizon ) {
    // this should only happen once
    assert( v_node_injection_constraints.size() == 0 );

    v_node_injection_constraints.resize
      ( boost::multi_array<FRowConstraint *, 2>::
        extent_gen()[f_time_horizon][num_nodes] );
  }

  // Node injection constraints.

  for( int t = 0; t < f_time_horizon; ++t ) {

    auto node_injection = v_network_blocks[t]->get_node_injection();

    for( int node_id = 0; node_id < num_nodes; ++node_id ) {

      auto linear_function = new LinearFunction();

      for( auto unit_block : f_network.get_node( node_id )->get_unit_blocks() )
        linear_function->add_variable( unit_block->get_power( t ), 1.0 );

      linear_function->add_variable( & node_injection[ node_id ], - 1.0 );

      v_node_injection_constraints[t][node_id]->set_both( 0.0 );
      v_node_injection_constraints[t][node_id]->set_function( linear_function );

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

  if( v_PollutantDemand_Const.size() != f_number_pollutants ) {
    // this should only happen once
    assert(v_PollutantDemand_Const.size() == 0);

    v_PollutantDemand_Const.resize
            (boost::multi_array<FRowConstraint *, 2>::
             extent_gen()[f_number_pollutants][f_number_emission_zones]);
  }

//TODO







}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------ METHODS FOR LOADING, PRINTING & SAVING THE UCBlock ----------*/
/*--------------------------------------------------------------------------*/

template<class T>
void serialize( netCDF::NcGroup & group, const std::string & var_name,
                const netCDF::NcType & ncType, const netCDF::NcDim & ncDim,
                const std::vector<T> & data ) {

  group.addVar( var_name , ncType , ncDim )
    .putVar( { 0 } , { data.size() } , data.data() );
}

/*--------------------------------------------------------------------------*/

void UCBlock::serialize( netCDF::NcGroup & group ) const {

  group.putAtt( "type" , "UCBlock" );

  group.addDim( "TimeHorizon",    f_time_horizon );
  group.addDim( "NumberUnits",    f_number_units );
  group.addDim( "NumberNodes",    f_number_nodes );

  ::serialize( group, "PrimaryDemand", netCDF::NcDouble(),
               group.addDim( "PrimaryZones", f_number_primary_zones ),
               f_primary_demand );

  ::serialize( group, "SecondaryDemand", netCDF::NcDouble(),
               group.addDim( "SecondaryZones", f_number_secondary_zones ),
               f_secondary_demand );

  ::serialize( group, "InertiaDemand", netCDF::NcDouble(),
               group.addDim( "InertiaZones", f_number_inertia_zones ),
               f_inertia_demand );

  ::serialize( group, "PollutantDemand", netCDF::NcDouble(),
               group.addDim( "PollutantSet", f_number_pollutants ),
               f_pollutant_demand );

  // TODO PollutantRho
  //group.addDim( "EmissionZones",  f_number_emission_zones );

}  // end( UCBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------------------ End File UCBlock.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
