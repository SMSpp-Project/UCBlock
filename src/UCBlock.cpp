/*--------------------------------------------------------------------------*/
/*----------------------- File UCBlock.cpp ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the UCBlock class.
 *
 * \version 0.11
 *
 * \date 25 - 03 - 2020
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
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <iostream>
#include <vector>
#include "FRowConstraint.h"
#include "HeatBlock.h"
#include "LinearFunction.h"
#include "NetworkBlock.h"
#include "UCBlock.h"
#include "UnitBlock.h"

#include "BusNetworkBlock.h"
#include "DCNetworkBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

SMSpp_insert_in_factory_cpp_1( UCBlock );
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

UCBlock::UCBlock( Block * father ) : Block( father ) {
 f_time_horizon = 0;
 f_number_units = 0;
 f_NetworkData = nullptr;
 f_number_heat_blocks = 0;
 f_number_primary_zones = 0;
 f_number_secondary_zones = 0;
 f_number_inertia_zones = 0;
 f_number_pollutants = 0;
}

UnitBlock * UCBlock::get_unit_block( Index i ) const {
 return dynamic_cast<UnitBlock *>( v_Block[ i ] );
}

NetworkBlock * UCBlock::get_network_block( Index t ) const {
 return dynamic_cast<NetworkBlock *>( v_Block[ f_number_units + t ] );
}

/*--------------------------------------------------------------------------*/

void UCBlock::deserialize_sub_blocks( const netCDF::NcGroup & group ) {

 for( auto block : v_Block )
  delete block;

 v_Block.clear();

 deserialize_sub_blocks( group, "UnitBlock_", f_number_units );
 v_network_blocks.resize( f_time_horizon );
 deserialize_network_blocks( group, f_time_horizon );
 deserialize_sub_blocks( group, "HeatBlock_", f_number_heat_blocks );

}

/*--------------------------------------------------------------------------*/

void UCBlock::deserialize_sub_blocks( const netCDF::NcGroup & group,
                                      const std::string & sub_group_name_prefix,
                                      const int num_sub_blocks ) {

 for( int i = 0; i < num_sub_blocks; ++i ) {

  std::string sub_group_name = sub_group_name_prefix + std::to_string( i );
  auto sub_group = group.getGroup( sub_group_name );

  if( sub_group.isNull() ) {
   throw ( std::invalid_argument( "UCBlock::deserialize: " +
                                  sub_group_name + " is not present" ) );
  }

  // std::string class_name;
  // class_name_attribute.getValues( class_name );
  // auto sub_block = new_Block( class_name, this );
  // sub_block->deserialize( sub_group );
  v_Block.push_back( new_Block( sub_group, this ) );
 }
}

/*--------------------------------------------------------------------------*/

void UCBlock::deserialize_network_blocks( const netCDF::NcGroup & group,
                                          int num_sub_blocks ) {

 for( int i = 0; i < num_sub_blocks; ++i ) {

  std::string sub_group_name = "NetworkBlock_" + std::to_string( i );
  auto sub_group = group.getGroup( sub_group_name );

  if( sub_group.isNull() ) {
   return;
  }

  auto class_name_attribute = sub_group.getAtt( "type" );

  if( class_name_attribute.isNull() ) {
   throw ( std::invalid_argument
           ( "UCBlock::deserialize: type attribute "
             "is not present in group " + sub_group_name ) );
  }

  std::string class_name;
  class_name_attribute.getValues( class_name );
  auto sub_block = new_Block( class_name, this );
  sub_block->deserialize( sub_group );
  v_Block.push_back( sub_block );
  v_network_blocks[ i ] = dynamic_cast<NetworkBlock *>(sub_block);
 }
}

/*--------------------------------------------------------------------------*/

void UCBlock::deserialize( netCDF::NcGroup & group ) {

 Index number_nodes = 1;
 if( ! ::deserialize_dim( group , "NumberNodes" , number_nodes , true ) )
  number_nodes = 1;

 if( number_nodes > 1 ) {
  delete f_NetworkData;
  f_NetworkData = new NetworkBlock::NetworkData();
  f_NetworkData->deserialize( group );
 }

 ::deserialize_dim( group, "TimeHorizon", f_time_horizon, false );
 ::deserialize_dim( group, "NumberUnits", f_number_units, false );

 ::deserialize_dim( group, "NumberHeatGenerators",
                    f_number_heat_generators,  true );
 ::deserialize_dim( group, "TotalNumberPollutantZones",
                    f_total_number_pollutant_zones,        true );

 boost::multi_array< double, 2 > v_active_power_demand;
 bool ap_found = ::deserialize( group, "ActivePowerDemand",
                                v_active_power_demand, true, false );
 transpose(v_active_power_demand);

 // Default values for optional dimensions
 f_number_heat_blocks = 0;
 f_number_primary_zones = 0;
 f_number_secondary_zones = 0;
 f_number_inertia_zones = 0;
 f_number_pollutants = 0;

 ::deserialize_dim( group, "NumberHeatBlocks", f_number_heat_blocks, true );
 ::deserialize_dim( group, "NumberPrimaryZones", f_number_primary_zones, true );
 ::deserialize_dim( group, "NumberSecondaryZones", f_number_secondary_zones, true );
 ::deserialize_dim( group, "NumberInertiaZones", f_number_inertia_zones, true );
 ::deserialize_dim( group, "NumberPollutants", f_number_pollutants, true );


 ::deserialize( group, "HeatSet", { f_number_units, f_number_heat_blocks },
                v_heat_set, true , false);

 ::deserialize( group, "PrimaryZones", number_nodes,
                v_primary_zones, true , true );

 ::deserialize( group, "PrimaryDemand", v_primary_demand, true, false );


 ::deserialize( group, "SecondaryZones", number_nodes,
                v_secondary_zones, true , true );

 ::deserialize( group, "SecondaryDemand",
                v_secondary_demand, true , false );

 ::deserialize( group, "InertiaZones", number_nodes,
                v_inertia_zones, true , true );

 ::deserialize( group, "InertiaDemand",
                v_inertia_demand, true , false );

 ::deserialize( group, "NumberPollutantZones", f_number_pollutants,
                v_number_pollutant_zones, true , true );

 ::deserialize( group, "PollutantZones",
                v_pollutant_zones, true , true );

 ::deserialize( group, "PollutantBudget", {f_total_number_pollutant_zones},
                v_pollutant_budget, true , false );

 ::deserialize( group, "PollutantRho",
                v_pollutant_rho, true , true );
 ::deserialize( group, "PollutantHeatRho", v_pollutant_heat_rho, true , true );


 ::deserialize( group, "PowerHeatRho", f_number_units, v_power_heat_rho, true , true );

 ::deserialize( group, "HeatNode", f_number_heat_blocks, v_heat_node, true , true );

 // TODO
 /* Notice that for units into a HeatBlock that also are electrical
  * units, the v_heat_node variable provides another time an
  * information that is already known, i.e., to which node they
  * belong to. Of course *the two information must agree*,
  * otherwise the input file is ill-defined and exception is
  * thrown.
  */


 deserialize_sub_blocks( group );

 if( !ap_found ) {
  for( auto i : v_network_blocks ) {
   if( i == nullptr )
    throw ( std::invalid_argument
            ( "UCBlock::deserialize: ActivePowerDemand is "
              "mandatory if NetworkBlocks are not specified" ) );
  }
  return;
 }

 // ActivePowerDemand was found, use it to populate NetworkBlocks
 for( Index i = 0; i < f_time_horizon; ++i ) {
  NetworkBlock * sub_block;

  if( number_nodes == 1 ) {
   // BusNetworkBlock

   if( v_network_blocks[ i ] ) {
    // A NetworkBlock already exists
    sub_block = v_network_blocks[ i ];
    if( !sub_block->get_NetworkData() && f_NetworkData ) {
     sub_block->set_NetworkData( f_NetworkData );
    }
    if( sub_block->get_active_demand().empty() ) {
     sub_block->set_ActiveDemand( { v_active_power_demand[ 0 ][ i ] } );
    }

   } else {
    // Create a new BusNetworkBlock
    sub_block = new BusNetworkBlock( this );
    v_Block.push_back( sub_block );
    v_network_blocks[ i ] = dynamic_cast<NetworkBlock *>(sub_block);
    if( f_NetworkData ) {
     sub_block->set_NetworkData( f_NetworkData );
    }
    sub_block->set_ActiveDemand( { v_active_power_demand[ 0 ][ i ] } );
   }


  } else { // number_nodes > 1
   // DCNetworkBlock

   typedef boost::multi_array_types::index_range range;
   auto ap_c = v_active_power_demand[ boost::indices[ range( 0, number_nodes ) ][ i ] ];
   std::vector< double > ap_v( number_nodes );
   std::copy( ap_c.begin(), ap_c.end(), ap_v.begin() );

   if( v_network_blocks[ i ] ) {
    // A NetworkBlock already exists
    sub_block = v_network_blocks[ i ];
    if( !sub_block->get_NetworkData() ) {
     // If we are here, we know that f_NetworkData is not null
     sub_block->set_NetworkData( f_NetworkData );
    }
    if( sub_block->get_active_demand().empty() ) {
     sub_block->set_ActiveDemand( ap_v );
    }

   } else {
    // Create a new DCNetworkBlock
    sub_block = new DCNetworkBlock( this );
    v_Block.push_back( sub_block );
    v_network_blocks[ i ] = dynamic_cast<NetworkBlock *>(sub_block);
    sub_block->set_NetworkData( f_NetworkData );
    sub_block->set_ActiveDemand( ap_v );
   }
  }
 }

 if( ! ::deserialize_dim( group, "NumberElectricalGenerators", f_number_elc_generators, true ) ) {
  f_number_elc_generators = 0;
  for( auto sub_block : get_nested_Blocks() ) {
   if( auto unit_block = dynamic_cast< UnitBlock * >( sub_block ) )
    f_number_elc_generators += unit_block->get_number_generators();
  }
 }

 ::deserialize( group, "GeneratorNode", f_number_elc_generators,
                v_generator_node, true , true );

 Block::deserialize( group );

}  // end( UCBlock::deserialize )

/*--------------------------------------------------------------------------*/

void UCBlock::generate_abstract_constraints( Configuration * stcc ) {

 Block::generate_abstract_constraints(stcc);

 unsigned int number_nodes = f_NetworkData ? f_NetworkData->
         get_number_nodes() : 1;

 // Node injection constraints.

 if( v_node_injection_constraints.size() != f_time_horizon ) {
  // this should only happen once
  assert( v_node_injection_constraints.empty() );

  v_node_injection_constraints.resize
          ( boost::multi_array< FRowConstraint, 2 >::
            extent_gen()[ f_time_horizon ][ number_nodes ] );
 }

 if( number_nodes > 0 ) {

  if ( number_nodes == 1) { //BusNetwork no need to GeneratorNode

   for( Index t = 0; t < f_time_horizon; ++t ) {

    auto & node_injection = v_network_blocks[t]->get_node_injection();
    auto linear_function = new LinearFunction();

    linear_function->add_variable( &node_injection[0], -1.0, eNoMod );

    v_node_injection_constraints[t][0].set_both( 0.0 );

    Index generator_id = 0;
    Index unit_id = 0;

    for( auto block : get_nested_Blocks() ) {
     auto unit_block = dynamic_cast<UnitBlock *>(block);
     if( unit_block == nullptr )
      continue;
     unit_block->get_number_generators();

    for( Index g = 0; g < unit_block->get_number_generators(); ++g ) {

       auto fixed_consumption = unit_block->get_fixed_consumption( g );

       auto ap = unit_block->get_active_power( g );
       auto active_power = &ap[t];

       auto c = unit_block->get_commitment( g );
       auto commitment = &c[t];

       linear_function->add_variable( active_power, 1.0, eNoMod );

       if( c != nullptr ) {
        if( fixed_consumption != nullptr ) {
         linear_function->add_variable( commitment, -fixed_consumption[t], eNoMod );
        } else {
         linear_function->add_variable( commitment, 0.0, eNoMod );
        }
       }
       if( fixed_consumption != nullptr ) {
        v_node_injection_constraints[t][0].set_both
                ( v_node_injection_constraints[t][0].get_rhs()
                  - fixed_consumption[t] );
       } else {
        v_node_injection_constraints[t][0].set_both
                ( v_node_injection_constraints[t][0].get_rhs()
                  - 0.0 );
       }
       generator_id++;
      }
     }
    v_node_injection_constraints[ t ][ 0 ].set_function( linear_function );
    //}
   }
  } else {  //DCNetwork needs GeneratorNode

   for( Index t = 0; t < f_time_horizon; ++t ) {

    auto & node_injection = v_network_blocks[t]->get_node_injection();

    for( Index node_id = 0; node_id < number_nodes; ++node_id ) {

     auto linear_function = new LinearFunction();

     linear_function->add_variable( &node_injection[node_id], -1.0, eNoMod );

     v_node_injection_constraints[t][node_id].set_both( 0.0 );

     Index generator_id = 0;
     Index unit_id = 0;

     for( auto block : get_nested_Blocks()) {
      auto unit_block = dynamic_cast<UnitBlock *>(block);
      if( unit_block == nullptr )
       continue;
      unit_block->get_number_generators();

      for( Index g = 0; g < unit_block->get_number_generators(); ++g  ) { //TODO Cheek it again

       Index generator_node = 0;

       if( node_id == v_generator_node[g] ) {

        generator_node = v_generator_node[g];

        auto fixed_consumption = unit_block->get_fixed_consumption( generator_node );

        auto ap = unit_block->get_active_power( generator_node );
        auto active_power = &ap[t];

        auto c = unit_block->get_commitment( generator_node );
        auto commitment = &c[t];

        linear_function->add_variable( active_power, 1.0, eNoMod );
        if( c != nullptr ) {
         if( fixed_consumption != nullptr ) {
          linear_function->add_variable( commitment, -fixed_consumption[t], eNoMod );
         } else {
          linear_function->add_variable( commitment, 0.0, eNoMod );
         }
        }
        if( fixed_consumption != nullptr ) {
         v_node_injection_constraints[t][node_id].set_both
                 ( v_node_injection_constraints[t][node_id].get_rhs()
                   - fixed_consumption[t] );
        } else {
         v_node_injection_constraints[t][node_id].set_both
                 ( v_node_injection_constraints[t][node_id].get_rhs()
                   - 0.0 );
        }
        generator_node++;
       }
      }
     }
     v_node_injection_constraints[t][node_id].set_function( linear_function );
    }
   }
  }
  add_static_constraint( v_node_injection_constraints, "node_injection_c" );
 }

/*--------------------------------------------------------------------------*/

// initial conditions

 unsigned int number_primary_zones = f_number_primary_zones ;
 unsigned int number_secondary_zones = f_number_secondary_zones ;
 unsigned int number_inertia_zones = f_number_inertia_zones ;


 // Primary demand constraints. For BusNetwork

 if( number_primary_zones > 0 ) {

  if( v_PrimaryDemand_Const.size() != f_time_horizon ) {
   // this should only happen once
   assert( v_PrimaryDemand_Const.empty());

   v_PrimaryDemand_Const.resize
           ( boost::multi_array< FRowConstraint, 2 >::
             extent_gen()[f_time_horizon][number_primary_zones] );
  }

  for( Index t = 0; t < f_time_horizon; ++t ) {

   for( Index zone_id = 0; zone_id < number_primary_zones; ++zone_id ) {

    v_PrimaryDemand_Const[t][zone_id].set_lhs
            ( get_primary_demand()[t][zone_id] );
    v_PrimaryDemand_Const[t][zone_id].set_rhs( Inf< double >());
    v_PrimaryDemand_Const[t][zone_id].set_function( new LinearFunction());
   }

   Index generator_id = 0;
   for( auto block : get_nested_Blocks()) {
    auto unit_block = dynamic_cast<UnitBlock *>(block);
    if( unit_block == nullptr )
     continue;
    unit_block->get_number_generators();

    for( Index g = 0; g < unit_block->get_number_generators(); ++g ) {

     auto primary_s_r = unit_block->get_primary_spinning_reserve( g );
     auto primary_spinning_reserve = &primary_s_r[t];

     auto linear_function = static_cast<LinearFunction *>
     ( v_PrimaryDemand_Const[ t ][ 0 ].get_function() );

     linear_function->
             add_variable( primary_spinning_reserve, 1.0 );

    }
    generator_id++;
   }

  }
  add_static_constraint( v_PrimaryDemand_Const );
 }

 // Secondary demand constraints.

 if( number_secondary_zones > 0 ) {

  if( v_SecondaryDemand_Const.size() != f_time_horizon ) {
   // this should only happen once
   assert( v_SecondaryDemand_Const.empty() );

   v_SecondaryDemand_Const.resize
           ( boost::multi_array< FRowConstraint, 2 >::
             extent_gen()[ f_time_horizon ][ number_secondary_zones ] );
  }

  for( Index t = 0; t < f_time_horizon; ++t ) {

   for( Index zone_id = 0; zone_id < number_secondary_zones; ++zone_id ) {
    v_SecondaryDemand_Const[t][zone_id].set_lhs
            ( get_secondary_demand()[t][zone_id] );
    v_SecondaryDemand_Const[t][zone_id].set_rhs( Inf< double >());
    v_SecondaryDemand_Const[t][zone_id].
            set_function( new LinearFunction());
   }
   Index generator_id = 0;
   for( auto block : get_nested_Blocks()) {
    auto unit_block = dynamic_cast<UnitBlock *>(block);
    if( unit_block == nullptr )
     continue;

    unit_block->get_number_generators();

    for( Index g = 0; g < unit_block->get_number_generators(); ++g ) {

     auto secondary_s_r = unit_block->get_secondary_spinning_reserve( g );
     auto secondary_spinning_reserve = &secondary_s_r[t];


     auto linear_function = dynamic_cast<LinearFunction *>
     ( v_SecondaryDemand_Const[t][0].get_function());
     linear_function->add_variable( secondary_spinning_reserve, 1.0 );
     generator_id++;

    }

   }
  }
  add_static_constraint( v_SecondaryDemand_Const );
 }


 // Inertia demand constraints.

 if( number_inertia_zones > 0 ) {

  if( v_InertiaDemand_Const.size() != f_time_horizon ) {
   // this should only happen once
   assert( v_InertiaDemand_Const.empty() );

   v_InertiaDemand_Const.resize
           ( boost::multi_array< FRowConstraint *, 2 >::
             extent_gen()[ f_time_horizon ][ number_inertia_zones ] );
  }

  for( Index t = 0; t < f_time_horizon; ++t ) {

   for( Index zone_id = 0; zone_id < number_inertia_zones; ++zone_id ) {
    v_InertiaDemand_Const[t][zone_id].set_lhs
            ( get_inertia_demand()[zone_id][t] );
    v_InertiaDemand_Const[t][zone_id].set_rhs( Inf< double >());
    v_InertiaDemand_Const[t][zone_id].set_function
            ( new LinearFunction());
   }

   Index generator_id = 0;
   for( auto block : get_nested_Blocks()) {
    auto unit_block = dynamic_cast<UnitBlock *>(block);
    if( unit_block == nullptr )
     continue;
    for( Index g = 0; g < unit_block->get_number_generators(); ++g ) {

     auto node_id = get_generator_node()[g];

     assert( node_id >= 0 && node_id < number_nodes );

     auto zone_id = get_inertia_zone()[node_id];
     if( zone_id >= number_inertia_zones )
      continue; // this unit does not belong to any zone

     auto ap = unit_block->get_active_power( g );
     auto active_power = &ap[t];

     auto c = unit_block->get_commitment( g );
     auto commitment = &c[t];

     auto inertia_commitment = unit_block->get_inertia_commitment(g);
     auto inertia_power = unit_block->get_inertia_power(g);

     auto linear_function = dynamic_cast<LinearFunction *>
     ( v_InertiaDemand_Const[t][zone_id].get_function());

     if (inertia_commitment != nullptr && c != nullptr ) {

      linear_function->add_variable( commitment , inertia_commitment[t] );
     }

     if( ap != nullptr && inertia_power != nullptr ) {
      linear_function->add_variable( active_power , inertia_power[t] );
     }
    }
   }
  }
  add_static_constraint( v_InertiaDemand_Const );
 }
/*
 // Pollutant budget constraints.

 // TODO This constraint should check again

 if( f_total_number_pollutant_zones > 0 ) {

  if( v_PollutantBudget_Const.size() != f_total_number_pollutant_zones ) {
   // this should only happen once
   assert( v_PollutantBudget_Const.empty() );

   v_PollutantBudget_Const.resize( f_total_number_pollutant_zones );
  }

  for( Index pollutant = 0; pollutant < f_number_pollutants; ++pollutant ) {
   for( Index zone = 0; zone < v_number_pollutant_zones[ pollutant ];
        ++zone ) {
    v_PollutantBudget_Const[ pollutant ][ zone ].set_rhs
     ( v_pollutant_budget[v_number_pollutant_zones[pollutant]][ pollutant ] );
    v_PollutantBudget_Const[ pollutant ][ zone ].set_lhs( -Inf< double >() );
    v_PollutantBudget_Const[ pollutant ][ zone ].set_function
     ( new LinearFunction() );
   }
  }

  for( Index pollutant = 0; pollutant < f_number_pollutants; ++pollutant ) {

   for( Index t = 0; t < f_time_horizon; ++t ) {

    // Terms associated with active power
    for( Index unit_id = 0; unit_id < f_number_units; ++unit_id ) {

     auto node_id = get_generator_node()[ unit_id ];
     auto zone_id = get_pollutant_zone()[ pollutant ][ node_id ];

     if( zone_id >= v_number_pollutant_zones[ pollutant ] )
      continue; // this unit does not belong to any zone

     auto rho = get_pollutant_rho()[ t ][ pollutant ][ unit_id ];
     auto active_power = get_unit_block( unit_id )->get_active_power() [ t ];

     auto linear_function = dynamic_cast<LinearFunction *>
     ( v_PollutantBudget_Const[ pollutant ][ zone_id ].get_function());

     //linear_function->add_variable( &active_power, rho );
    }

    // Terms associated with heat-only generation units

    if( f_number_heat_blocks > 0 ) {

     for( Index h = 0; h < f_number_heat_blocks; ++h ) {

      for( std::vector< Index >::size_type i = 0;
           i <= f_number_units; ++i ) {
       //TODO
       //auto unit_id = v_heat_only_units[ i ];
       auto heat_id = v_heat_set[ i ];

       auto zone_id = get_pollutant_zone()[ pollutant ][ h ];
       if( zone_id >= v_number_pollutant_zones[ pollutant ] )
        continue; // this unit does not belong to any zone

       auto heat = get_heat_block() [ h ]->get_heat()[ t ][ i ];
       auto rho = get_pollutant_heat_rho()[ t ][ pollutant ][ h ];

       auto linear_function = dynamic_cast<LinearFunction *>
       ( v_PollutantBudget_Const[ pollutant ][ zone_id ].get_function());

       linear_function->add_variable( &heat, rho );
      }
     }

    }
   }

   add_static_constraint( v_PollutantBudget_Const[ f_number_pollutants ] );
  }
 }

 // Heat constraints.

 if( f_number_heat_blocks > 0 ) {

  Index num_constraints_per_time = 0;

  // List of units that produce electricity and belong to some
  // HeatBlock.
  std::vector< Index > electricity_generators_inside_a_heat_block;

  if( v_power_Heat_Rho_Const.size() != f_time_horizon ) {

   // this should only happen once
   assert( v_power_Heat_Rho_Const.empty() );

   auto is_electricity_producing_inside_a_heat_block =
    [ this ]( Index unit ) {
     for( Index heat_block_id = 0; heat_block_id < f_number_heat_blocks;
          ++heat_block_id ) {
      if( get_heat_set()[ unit ] <
          v_heat_blocks[ heat_block_id ]->get_number_heat_generators() )
       return true;
     }
     return false;
    };

   for( Index unit_id = 0; unit_id < f_number_units; ++unit_id ) {
    if( is_electricity_producing_inside_a_heat_block( unit_id ) ) {
     electricity_generators_inside_a_heat_block
     [ num_constraints_per_time++ ] = unit_id;
    }
   }

   v_power_Heat_Rho_Const.resize
    ( boost::multi_array< FRowConstraint *, 2 >::
      extent_gen()[ f_time_horizon ][ num_constraints_per_time ] );
  }

  for( Index t = 0; t < f_time_horizon; ++t ) {

   for( std::vector< Index >::size_type constraint_id = 0;
        constraint_id < num_constraints_per_time; ++constraint_id ) {

    auto generator_id = electricity_generators_inside_a_heat_block[ constraint_id ];

    for( Index heat_block_id = 0; heat_block_id < f_number_heat_blocks;
         ++heat_block_id ) {

     auto heat_unit_id = get_heat_set()[ generator_id ];

     if( heat_unit_id >=
         v_heat_blocks[ heat_block_id ]->get_number_heat_generators() )
      continue; // unit_id does not belong to heat_block_id

     if( !v_power_Heat_Rho_Const[ t ][ constraint_id ].get_function() ) {
      v_power_Heat_Rho_Const[ t ][ constraint_id ].
       set_lhs( -Inf< double >() );
      v_power_Heat_Rho_Const[ t ][ constraint_id ].set_rhs( 0.0 );
      v_power_Heat_Rho_Const[ t ][ constraint_id ].set_function
       ( new LinearFunction() );
     }

     auto active_power = get_unit_block( generator_id )->get_active_power();
     auto heat = get_heat_block() [ heat_unit_id ]->get_heat()[ t ][generator_id];
     auto power_heat_rho = get_power_heat_rho()[ generator_id ];

     auto linear_function = dynamic_cast<LinearFunction *>
     ( v_power_Heat_Rho_Const[ t ][ constraint_id ].get_function());
     linear_function->add_variable( &heat, 1, 0 );
     linear_function->add_variable( &active_power[ t ][ generator_id ], -power_heat_rho );
    }
   }
  }

  add_static_constraint( v_power_Heat_Rho_Const );
 }
*/
/*--------------------------------------------------------------------------*/

}  // end( UCBlock::global constraints )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------ METHODS FOR LOADING, PRINTING & SAVING THE UCBlock ----------*/
/*--------------------------------------------------------------------------*/

void UCBlock::serialize( netCDF::NcGroup & group ) const {

 Block::serialize( group );

 auto dim_time_horizon = group.addDim( "TimeHorizon", f_time_horizon );
 auto dim_number_units = group.addDim( "NumberUnits", f_number_units );
 auto dim_number_elc_generators = group.addDim( "NumberElectricalGenerators",
                                                f_number_elc_generators);
 auto dim_number_heat_generators = group.addDim( "NumberHeatGenerators",
                                                 f_number_heat_generators );
 auto dim_total_number_pollutant_zone = group.addDim( "TotalNumberPollutantZones",
                                                      f_total_number_pollutant_zones );
 auto dim_number_nodes = group.addDim( "NumberNodes",
                                       f_NetworkData ? f_NetworkData->get_number_nodes() : 1 );


 auto dim_number_heat_blocks =
         group.addDim( "NumberHeatBlocks", f_number_heat_blocks );
 auto dim_number_primary_zones =
         group.addDim( "NumberPrimaryZones", f_number_primary_zones );
 auto dim_number_secondary_zones =
         group.addDim( "NumberSecondaryZones", f_number_secondary_zones );
 auto dim_number_inertia_zones =
         group.addDim( "NumberInertiaZones", f_number_inertia_zones );
 auto dim_number_pollutants =
         group.addDim( "NumberPollutants", f_number_pollutants );

 ::serialize( group, "HeatSet", netCDF::NcUint64(),
              { dim_number_units, dim_number_heat_blocks }, v_heat_set);


 ::serialize( group, "PrimaryZones", netCDF::NcUint64(),
              {dim_number_nodes}, v_primary_zones, true );


 ::serialize( group, "PrimaryDemand", netCDF::NcDouble(),
              { dim_number_primary_zones, dim_time_horizon },
              v_primary_demand, false );


 ::serialize( group, "SecondaryZones", netCDF::NcUint64(),
              { dim_number_nodes }, v_secondary_zones, true );

 ::serialize( group, "SecondaryDemand", netCDF::NcDouble(),
              { dim_number_secondary_zones, dim_time_horizon },
              v_secondary_demand, false );


 ::serialize( group, "InertiaZones", netCDF::NcUint64(),
              { dim_number_nodes }, v_inertia_zones, true );

 ::serialize( group, "InertiaDemand", netCDF::NcDouble(),
              { dim_number_inertia_zones, dim_time_horizon } ,
              v_inertia_demand, false);

 ::serialize( group, "NumberPollutantZones", netCDF::NcUint64(),
              { dim_number_pollutants }, v_number_pollutant_zones, true );

 ::serialize( group, "PollutantZones", netCDF::NcUint64(),
              { dim_number_pollutants, dim_number_nodes },
              v_pollutant_zones, true );

 ::serialize( group, "PollutantBudget", netCDF::NcDouble(),
              { dim_total_number_pollutant_zone}, v_pollutant_budget, false );

 ::serialize( group, "PollutantRho", netCDF::NcDouble(),
              { dim_time_horizon, dim_number_pollutants, dim_number_units },
              v_pollutant_rho );

 ::serialize( group, "PollutantHeatRho", netCDF::NcDouble(),
              { dim_time_horizon, dim_number_pollutants,
                dim_number_heat_blocks },
              v_pollutant_heat_rho );

 ::serialize( group, "PowerHeatRho", netCDF::NcDouble(),
              { dim_number_units }, v_power_heat_rho );

 ::serialize( group, "GeneratorNode", netCDF::NcUint64(),
              { dim_number_elc_generators }, v_generator_node );


 ::serialize( group, "HeatNode", netCDF::NcUint64(),
              { dim_number_heat_blocks }, v_heat_node );


 // Serialize sub-blocks

 if (f_NetworkData) {
  auto sub_group = group.addGroup( "NetworkData" );
  f_NetworkData->serialize(sub_group);
 }

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

 for( Index i = 0; i < f_number_heat_blocks; ++i ) {
  auto sub_block = get_heat_block() [ i ];
  auto sub_group = group.addGroup( "HeatBlock_" + std::to_string( i ) );
  sub_block->serialize( sub_group );
 }
}  // end( UCBlock::serialize )

template< typename T >
void UCBlock::transpose( boost::multi_array< T, 2 > & a ) {
 long rows = a.shape()[ 0 ];
 long cols = a.shape()[ 1 ];
 if( rows > 1 && cols == 1 ) {
  // The vector must be transposed
  boost::array< typename boost::multi_array< T, 2 >::index, 2 > dims = { { 1, rows } };
  a.reshape( dims );
 }
}
/*--------------------------------------------------------------------------*/
/*------------------------ End File UCBlock.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
