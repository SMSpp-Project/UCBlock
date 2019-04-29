/*--------------------------------------------------------------------------*/
/*----------------------- File UCBlock.cpp ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the UCBlock class.
 *
 * \version 0.11
 *
 * \date 24 - 04 - 2019
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



/*--------------------------------------------------------------------------*/

void UCBlock::deserialize( netCDF::NcGroup & group , Block * father ) {


// check the data- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  netCDF::NcDim TH = group.getDim( "TimeHorizon" );
  if( TH.isNull() )
    throw( std::invalid_argument( "TimeHorizon not present" ) );

  f_time_horizon = TH.getSize();
  if( f_time_horizon <= 0 )
    throw( std::invalid_argument( "TimeHorizon <= 0" ) );

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  netCDF::NcDim NU = group.getDim( "NumberUnits" );
  if( NU.isNull() )
    throw( std::invalid_argument( "NumberUnits not present" ) );

  f_number_units = NU.getSize();
  if( f_number_units <= 0 )
    throw( std::invalid_argument( "NumberUnits <= 0" ) );
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  netCDF::NcDim NN = group.getDim( "NumberNodes" );
  if( NN.isNull() )
    throw( std::invalid_argument( "NumberNodes not present" ) );

  f_number_nodes = NN.getSize();
  if( f_number_nodes <= 0 )
    throw( std::invalid_argument( "NumberNodes <= 0" ) );

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  netCDF::NcDim PZ = group.getDim( "PrimaryZones" );
  if( PZ.isNull() )
    throw( std::invalid_argument( "PrimaryZones not present" ) );

  f_number_Pr_zones = PZ.getSize();
  if( f_number_Pr_zones <= 0 )
    throw( std::invalid_argument( "NumberNodes <= 0" ) );
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  netCDF::NcDim SZ = group.getDim( "SecondaryZones" );
  if( SZ.isNull() )
    throw( std::invalid_argument( "SecondaryZones not present" ) );

  f_number_Se_zones = SZ.getSize();
  if( f_number_Se_zones <= 0 )
    throw( std::invalid_argument( "SecondaryZones <= 0" ) );
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  netCDF::NcDim IZ = group.getDim( "InertiaZones" );
  if( IZ.isNull() )
    throw( std::invalid_argument( "InertiaZones not present" ) );

  f_number_In_zones = IZ.getSize();
  if( f_number_In_zones <= 0 )
    throw( std::invalid_argument( "InertiaZones <= 0" ) );
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  netCDF::NcDim EZ = group.getDim( "EmissionZones" );
  if( EZ.isNull() )
    throw( std::invalid_argument( "EmissionZones not present" ) );

  f_number_Em_zones = EZ.getSize();
  if( f_number_Em_zones <= 0 )
    throw( std::invalid_argument( "EmissionZones <= 0" ) );
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  netCDF::NcDim NP = group.getDim( "PollutantSet" );
  if( NP.isNull() )
    throw( std::invalid_argument( "PollutantSet not present" ) );

  f_number_pollutants = NP.getSize();
  if( f_number_pollutants <= 0 )
    throw( std::invalid_argument( "PollutantSet <= 0" ) );

// read problem data- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  std::vector < size_t > start = { 0 };
  std::vector < size_t > countime = { (size_t)f_time_horizon };
  std::vector < size_t > coununit = { (size_t)f_number_units };
  std::vector < size_t > counode = { (size_t)f_number_nodes };
  std::vector < size_t > countprim = { (size_t)f_number_Pr_zones };
  std::vector < size_t > countsecond = { (size_t)f_number_Se_zones};
  std::vector < size_t > counit = { (size_t)f_number_In_zones };
  std::vector < size_t > countunem = { (size_t)f_number_Em_zones };
  std::vector < size_t > countpoll = { (size_t)f_number_pollutants };



/*--------------------PrimaryDemand-deserialize-----------------------------*/

  netCDF::NcVar PrimDemand = group.getVar( "PrimaryDemand" );
  if( PrimDemand.isNull() )
    throw( std::invalid_argument( "Primary Demand not present" ) );

  f_prim_demand.resize( f_number_Pr_zones );

/*--------------------SecondaryDemand-deserialize---------------------------*/

/*--------------------InertiaDemand-deserialize-----------------------------*/

/*-------------------PollutantDemand-deserialize----------------------------*/

/*---------------------PollutantRho-deserialize-----------------------------*/






}  // end( UCBlock::deserialize )

/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

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
}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------ METHODS FOR LOADING, PRINTING & SAVING THE UCBlock ----------*/
/*--------------------------------------------------------------------------*/

void UCBlock::serialize( netCDF::NcGroup & group ) const
{
  group.putAtt( "type" , "UCBlock" );


  netCDF::NcDim time_horizon = group.addDim( "TimeHorizon" , f_time_horizon );
  netCDF::NcDim number_units = group.addDim( "NumberUnits" , f_number_units);

}  // end( UCBlock::serialize )
/*--------------------------------------------------------------------------*/
/*------------------------ End File UCBlock.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/