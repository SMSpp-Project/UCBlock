/*--------------------------------------------------------------------------*/
/*----------------------- File UCBlock.cpp ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the UCBlock class.
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

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;


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

 deserialize_sub_blocks( group, "UnitBlock", f_number_units );
 deserialize_sub_blocks( group, "NetworkBlock_", f_time_horizon );
 deserialize_sub_blocks( group, "HeatBlock_", f_number_heat_blocks );

}

/*--------------------------------------------------------------------------*/

void UCBlock::deserialize_sub_blocks
 ( const netCDF::NcGroup & group, const std::string & sub_group_name_prefix,
   const int num_sub_blocks ) {

 for( int i = 0; i < num_sub_blocks; ++i ) {

  std::string sub_group_name = sub_group_name_prefix + std::to_string( i );
  auto sub_group = group.getGroup( sub_group_name );

  if( sub_group.isNull() ) {
   throw ( std::invalid_argument( "UCBlock::deserialize: " +
                                  sub_group_name + " is not present" ) );
  }

  auto class_name_attribute = sub_group.getAtt( "ClassName" );

  if( class_name_attribute.isNull() ) {
   throw ( std::invalid_argument
    ( "UCBlock::deserialize: ClassName attribute "
      "is not present in group " + sub_group_name ) );
  }

  std::string class_name;
  class_name_attribute.getValues( class_name );
  auto sub_block = new_Block( class_name, this );
  sub_block->deserialize( sub_group );
  v_Block.push_back( sub_block );
 }
}

/*--------------------------------------------------------------------------*/

void UCBlock::deserialize( netCDF::NcGroup & group ) {

 auto network_data = new NetworkBlock::NetworkData(); //TODO not Implemented well
 network_data->deserialize( group );

 // FIXME delete f_NetworkData; // This is equivalent to code block below
 // --------------------------------------------------------------------------
 if( network_data ) {  // there is a NetworkData object in the group
  // use it, whatever has happened before
  // if there was a previous NetworkData, delete it
  if( f_NetworkData )
   delete f_NetworkData;
 } else if( !network_data ) {
  // if the NetworkData has not been passed from outside
  throw ( std::logic_error( "UCBlock has no NetworkData access" ) );
 }
 // --------------------------------------------------------------------------

 ::deserialize_dim( group, "TimeHorizon", f_time_horizon, false );
 ::deserialize_dim( group, "NumberUnits", f_number_units, false );
 ::deserialize_dim( group, "NumberElectricalGenerators",
                                 f_number_elc_generators, false );
 ::deserialize_dim( group, "NumberHeatGenerators",
                                f_number_heat_generators, false );

 unsigned int number_nodes = f_NetworkData ? f_NetworkData->
  get_number_nodes() : 1;

 ::deserialize_dim( group, "NumberNodes", number_nodes );

 // Default values for optional dimensions
 f_number_heat_blocks = 0;
 f_number_primary_zones = 0;
 f_number_secondary_zones = 0;
 f_number_inertia_zones = 0;
 f_number_pollutants = 0;

 ::deserialize_dim( group, "NumberHeatBlocks", f_number_heat_blocks );
 ::deserialize_dim( group, "NumberPrimaryZones", f_number_primary_zones );
 ::deserialize_dim( group, "NumberSecondaryZones", f_number_secondary_zones );
 ::deserialize_dim( group, "NumberInertiaZones", f_number_inertia_zones );
 ::deserialize_dim( group, "NumberPollutants", f_number_pollutants );

/*
  if( f_number_heat_blocks >= 1 ) {
    ::deserialize( group, "HeatSet", { f_number_units, f_number_heat_blocks },
                   v_heat_set);
  }
*/
 if( f_number_primary_zones >= 1 ) {
  ::deserialize( group, "PrimaryZones", number_nodes,
                 v_primary_zones );
/*
    ::deserialize( group, "PrimaryDemand",
                   { f_number_primary_zones, f_time_horizon },
                   v_primary_demand );*/
 }

 if( f_number_secondary_zones >= 1 ) {
  ::deserialize( group, "SecondaryZones", number_nodes,
                 v_secondary_zones );
/*
    ::deserialize( group, "SecondaryDemand",
                   { f_number_secondary_zones, f_time_horizon },
                   v_secondary_demand );  */
 }

 if( f_number_inertia_zones >= 1 ) {
  ::deserialize( group, "InertiaZones", number_nodes,
                 v_inertia_zones );
/*
    ::deserialize( group, "InertiaDemand",
                   { f_number_inertia_zones, f_time_horizon },
                   v_inertia_demand );  */
 }

 if( f_number_pollutants >= 1 ) {

  ::deserialize( group, "NumberPollutantZones", f_number_pollutants,
                 v_number_pollutant_zones );
/*
    ::deserialize( group, "PollutantZones",
                   { f_number_pollutants, number_nodes },
                   v_pollutant_zones );
*/
  ::deserialize( group, "PollutantBudget", f_number_pollutants,
                 v_pollutant_budget );
/*
    ::deserialize( group, "PollutantRho",
                   { f_time_horizon, f_number_pollutants, f_number_units },
                   v_pollutant_rho );

    if( f_number_heat_blocks >= 1 ) {

      ::deserialize( group, "PollutantHeatRho",
                     { f_time_horizon, f_number_pollutants,
                         f_number_heat_blocks }, v_pollutant_heat_rho );
    }*/
 }

 if( number_nodes > 1 ) {
  ::deserialize( group, "GeneratorNode", f_number_elc_generators,
                 v_generator_node );
 }

 ::deserialize( group, "PowerHeatRho", f_number_units, v_power_heat_rho );

 if( f_number_heat_blocks > 0 && f_number_pollutants > 0 ) {
  ::deserialize( group, "HeatNode", f_number_heat_blocks, v_heat_node );

  // TODO
  /* Notice that for units into a HeatBlock that also are electrical
   * units, the v_heat_node variable provides another time an
   * information that is already known, i.e., to which node they
   * belong to. Of course *the two information must agree*,
   * otherwise the input file is ill-defined and exception is
   * thrown.
   */
 }

 deserialize_sub_blocks( group );

}  // end( UCBlock::deserialize )

/*--------------------------------------------------------------------------*/

void UCBlock::generate_abstract_constraints( Configuration * stcc ) {

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

  for( Index t = 0; t < f_time_horizon; ++t ) {

   auto node_injection = v_network_blocks[ t ]->get_node_injection();

   for( Index node_id = 0; node_id < number_nodes; ++node_id ) {

    auto linear_function = new LinearFunction();

    linear_function->add_variable( &node_injection[ node_id ], -1.0 );

    v_node_injection_constraints[ t ][ node_id ].set_both( 0.0 );
    v_node_injection_constraints[ t ][ node_id ].
     set_function( linear_function );
   }

   for( Index unit_id = 0; unit_id < f_number_units; ++unit_id ) {

    Index generator_id = get_generator_node()[ generator_id ];
/* //TODO Fixe me
    auto fixed_consumption =
     get_unit_block( unit_id )->get_fixed_consumption()[ t ];

    v_node_injection_constraints[ t ][ node_id ].set_both
     ( v_node_injection_constraints[ t ][ node_id ].get_rhs()
       - fixed_consumption );

    auto linear_function = dynamic_cast<LinearFunction *>
    ( v_node_injection_constraints[ t ][ node_id ].get_function());

    auto power = &get_unit_block( unit_id )->get_active_power() [ t ];
    auto commitment = &( get_unit_block( unit_id )->get_commitment() [ t ] ;

    linear_function->add_variable( power, 1.0 );
    linear_function->add_variable( commitment, -fixed_consumption );*/
   }
  }

  add_static_constraint( v_node_injection_constraints );
 }

/*--------------------------------------------------------------------------*/

 // Primary demand constraints.

 if( f_number_primary_zones > 0 ) {

  if( v_PrimaryDemand_Const.size() != f_time_horizon ) {
   // this should only happen once
   assert( v_PrimaryDemand_Const.empty() );

   v_PrimaryDemand_Const.resize
    ( boost::multi_array< FRowConstraint, 2 >::
      extent_gen()[ f_time_horizon ][ f_number_primary_zones ] );
  }

  for( Index t = 0; t < f_time_horizon; ++t ) {

   for( Index zone_id = 0; zone_id < f_number_primary_zones; ++zone_id ) {
    v_PrimaryDemand_Const[ t ][ zone_id ].set_lhs
     ( get_primary_demand()[ zone_id ][ t ] );
    v_PrimaryDemand_Const[ t ][ zone_id ].set_rhs( Inf< double >() );
    v_PrimaryDemand_Const[ t ][ zone_id ].set_function
     ( new LinearFunction() );
   }

   for( Index unit_id = 0; unit_id < f_number_units; ++unit_id ) {

    auto node_id = get_generator_node()[ unit_id ];

    assert( node_id >= 0 && node_id < number_nodes );

    auto zone_id = get_primary_zone()[ node_id ];
    if( zone_id >= f_number_primary_zones )
     continue; // this unit does not belong to any zone

    auto primary_spinning_reserve = get_unit_block( unit_id )
     ->get_primary_spinning_reserve() [ t ];

    auto linear_function = dynamic_cast<LinearFunction *>
    ( v_PrimaryDemand_Const[ t ][ zone_id ].get_function());
  //  linear_function->
  //   add_variable( &primary_spinning_reserve, 1.0 );
   }
  }

  add_static_constraint( v_PrimaryDemand_Const );
 }

/*--------------------------------------------------------------------------*/

 // Secondary demand constraints.

 if( f_number_secondary_zones > 0 ) {

  if( v_SecondaryDemand_Const.size() != f_time_horizon ) {
   // this should only happen once
   assert( v_SecondaryDemand_Const.empty() );

   v_SecondaryDemand_Const.resize
    ( boost::multi_array< FRowConstraint, 2 >::
      extent_gen()[ f_time_horizon ][ f_number_secondary_zones ] );
  }

  for( Index t = 0; t < f_time_horizon; ++t ) {

   for( Index zone_id = 0; zone_id < f_number_secondary_zones; ++zone_id ) {
    v_SecondaryDemand_Const[ t ][ zone_id ].set_lhs
     ( get_secondary_demand()[ zone_id ][ t ] );
    v_SecondaryDemand_Const[ t ][ zone_id ].set_rhs( Inf< double >() );
    v_SecondaryDemand_Const[ t ][ zone_id ].
     set_function( new LinearFunction() );
   }

   for( Index unit_id = 0; unit_id < f_number_units; ++unit_id ) {

    auto node_id = get_generator_node()[ unit_id ];

    assert( node_id >= 0 && node_id < number_nodes );

    auto zone_id = get_secondary_zone()[ node_id ];
    if( zone_id >= f_number_secondary_zones )
     continue; // this unit does not belong to any zone

    auto secondary_spinning_reserve = get_unit_block( unit_id )->
     get_secondary_spinning_reserve() [ t ];

    auto linear_function = dynamic_cast<LinearFunction *>
    ( v_SecondaryDemand_Const[ t ][ zone_id ].get_function());
//    linear_function->add_variable( &secondary_spinning_reserve, 1.0 );
   }
  }

  add_static_constraint( v_SecondaryDemand_Const );
 }

/*--------------------------------------------------------------------------*/

 // Inertia demand constraints.

 if( f_number_inertia_zones > 0 ) {

  if( v_InertiaDemand_Const.size() != f_time_horizon ) {
   // this should only happen once
   assert( v_InertiaDemand_Const.empty() );

   v_InertiaDemand_Const.resize
    ( boost::multi_array< FRowConstraint *, 2 >::
      extent_gen()[ f_time_horizon ][ f_number_inertia_zones ] );
  }

  for( Index t = 0; t < f_time_horizon; ++t ) {

   for( Index zone_id = 0; zone_id < f_number_inertia_zones; ++zone_id ) {
    v_InertiaDemand_Const[ t ][ zone_id ].set_lhs
     ( get_inertia_demand()[ zone_id ][ t ] );
    v_InertiaDemand_Const[ t ][ zone_id ].set_rhs( Inf< double >() );
    v_InertiaDemand_Const[ t ][ zone_id ].set_function
     ( new LinearFunction() );
   }

   for( Index unit_id = 0; unit_id < f_number_units; ++unit_id ) {

    auto node_id = get_generator_node()[ unit_id ];

    assert( node_id >= 0 && node_id < number_nodes );

    auto zone_id = get_inertia_zone()[ node_id ];
    if( zone_id >= f_number_inertia_zones )
     continue; // this unit does not belong to any zone
/* //TODO FIX ME
    auto commitment_variable =
     &get_unit_block( unit_id )->get_commitment() [ t ];

    auto inertia_commitment =
     get_unit_block( unit_id )->get_inertia_commitment()[][ t ];

    auto active_power_variable =
     &get_unit_block( unit_id )->get_active_power() [ t ];

    auto inertia_power =
     get_unit_block( unit_id )->get_inertia_power()[ t ];

    auto linear_function = dynamic_cast<LinearFunction *>
    ( v_InertiaDemand_Const[ t ][ zone_id ].get_function());

    linear_function->
     add_variable( commitment_variable, inertia_commitment );
    linear_function->add_variable( active_power_variable, inertia_power ); */
   }
  }

  add_static_constraint( v_InertiaDemand_Const );
 }
/*--------------------------------------------------------------------------*/

 // Pollutant budget constraints.

 // TODO This constraint should check again

 if( f_number_pollutants > 0 ) {

  if( v_PollutantBudget_Const.size() != f_number_pollutants ) {
   // this should only happen once
   assert( v_PollutantBudget_Const.empty() );

   v_PollutantBudget_Const.resize( f_number_pollutants );
   for( Index pollutant = 0; pollutant < f_number_pollutants; ++pollutant ) {
    v_PollutantBudget_Const.resize( v_number_pollutant_zones[ pollutant ]
    );
   }
  }

  for( Index pollutant = 0; pollutant < f_number_pollutants; ++pollutant ) {
   for( Index zone = 0; zone < v_number_pollutant_zones[ pollutant ];
        ++zone ) {
    v_PollutantBudget_Const[ pollutant ][ zone ].set_rhs
     ( v_pollutant_budget[ pollutant ] );
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

   add_static_constraint( v_PollutantBudget_Const[ pollutant ] );
  }
 }
/*--------------------------------------------------------------------------*/

 // Heat constraints.

 if( f_number_heat_blocks > 0 ) {

  Index num_constraints_per_time = 0;

  // List of units that produce electricity and belong to some
  // HeatBlock.
  std::vector< Index > electricity_units_inside_a_heat_block;

  if( v_power_Heat_Rho_Const.size() != f_time_horizon ) {

   // this should only happen once
   assert( v_power_Heat_Rho_Const.empty() );

   auto is_electricity_producing_inside_a_heat_block =
    [ this ]( Index unit ) {
     for( Index heat_block_id = 0; heat_block_id < f_number_heat_blocks;
          ++heat_block_id ) {
      if( get_heat_set()[ unit ] <
          v_heat_blocks[ heat_block_id ]->get_number_heat_units() )
       return true;
     }
     return false;
    };

   for( Index unit_id = 0; unit_id < f_number_units; ++unit_id ) {
    if( is_electricity_producing_inside_a_heat_block( unit_id ) ) {
     electricity_units_inside_a_heat_block
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

    auto unit_id = electricity_units_inside_a_heat_block[ constraint_id ];

    for( Index heat_block_id = 0; heat_block_id < f_number_heat_blocks;
         ++heat_block_id ) {

     auto heat_unit_id = get_heat_set()[ unit_id ];

     if( heat_unit_id >=
         v_heat_blocks[ heat_block_id ]->get_number_heat_units() )
      continue; // unit_id does not belong to heat_block_id

     if( !v_power_Heat_Rho_Const[ t ][ constraint_id ].get_function() ) {
      v_power_Heat_Rho_Const[ t ][ constraint_id ].
       set_lhs( -Inf< double >() );
      v_power_Heat_Rho_Const[ t ][ constraint_id ].set_rhs( 0.0 );
      v_power_Heat_Rho_Const[ t ][ constraint_id ].set_function
       ( new LinearFunction() );
     }

     auto active_power = get_unit_block( unit_id )->get_active_power() [ t ];
     auto heat = get_heat_block() [ heat_unit_id ]->get_heat()[ t ][unit_id];
     auto power_heat_rho = get_power_heat_rho()[ unit_id ];

     auto linear_function = dynamic_cast<LinearFunction *>
     ( v_power_Heat_Rho_Const[ t ][ constraint_id ].get_function());
     linear_function->add_variable( &heat, 1, 0 );
//     linear_function->add_variable( &active_power, -power_heat_rho );
    }
   }
  }

  add_static_constraint( v_power_Heat_Rho_Const );
 }

/*--------------------------------------------------------------------------*/

}  // end( UCBlock::global constraints )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------ METHODS FOR LOADING, PRINTING & SAVING THE UCBlock ----------*/
/*--------------------------------------------------------------------------*/

void UCBlock::serialize( netCDF::NcGroup & group ) const {

 group.putAtt( "type", "UCBlock" );

 auto dim_time_horizon = group.addDim( "TimeHorizon", f_time_horizon );
 auto dim_number_units = group.addDim( "NumberUnits", f_number_units );
 auto dim_number_elc_generators = group.addDim( "NumberElectricalGenerators",
         f_number_elc_generators);
 auto dim_number_heat_generators = group.addDim( "NumberHeatGenerators",
                                       f_number_heat_generators );
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

/*
  if( f_number_heat_blocks >= 1 ) {
    ::serialize( group, "HeatSet", netCDF::NcUint64(),
                 { dim_number_units, dim_number_heat_blocks }, v_heat_set );
  }*/

 if( f_number_primary_zones >= 1 ) {
  ::serialize( group, "PrimaryZones", netCDF::NcUint64(),
               { dim_number_nodes }, v_primary_zones );
/*
    ::serialize( group, "PrimaryDemand", netCDF::NcDouble(),
                 { dim_number_primary_zones, dim_time_horizon },
                 v_primary_demand ); */
 }

 if( f_number_secondary_zones >= 1 ) {
  ::serialize( group, "SecondaryZones", netCDF::NcUint64(),
               { dim_number_nodes }, v_secondary_zones );
/*
    ::serialize( group, "SecondaryDemand", netCDF::NcDouble(),
                 { dim_number_secondary_zones, dim_time_horizon },
                 v_secondary_demand ); */
 }

 if( f_number_inertia_zones >= 1 ) {
  ::serialize( group, "InertiaZones", netCDF::NcUint64(),
               { dim_number_nodes }, v_inertia_zones );
/*
    ::serialize( group, "InertiaDemand", netCDF::NcDouble(),
                 { dim_number_inertia_zones, dim_time_horizon } ,
                 v_inertia_demand ); */
 }

 if( f_number_pollutants >= 1 ) {

  ::serialize( group, "NumberPollutantZones", netCDF::NcUint64(),
               { dim_number_pollutants }, v_number_pollutant_zones );
/*
     ::serialize( group, "PollutantZones", netCDF::NcUint64(),
                 { dim_number_pollutants, dim_number_nodes },
                 v_pollutant_zones );

    ::serialize( group, "PollutantBudget", netCDF::NcDouble(),
                 { dim_number_pollutants }, v_pollutant_budget );

    ::serialize( group, "PollutantRho", netCDF::NcDouble(),
                 { dim_time_horizon, dim_number_pollutants, dim_number_units },
                 v_pollutant_rho );

    ::serialize( group, "PollutantHeatRho", netCDF::NcDouble(),
                 { dim_time_horizon, dim_number_pollutants,
                     dim_number_heat_blocks },
                 v_pollutant_heat_rho ); */
 }

 ::serialize( group, "PowerHeatRho", netCDF::NcDouble(),
              { dim_number_units }, v_power_heat_rho );

 if( f_NetworkData->get_number_nodes() > 1 ) {
  ::serialize( group, "GeneratorNode", netCDF::NcUint64(),
               { dim_number_elc_generators }, v_generator_node );
 }

 if( f_number_heat_blocks > 0 && f_number_pollutants > 0 ) {
  ::serialize( group, "HeatNode", netCDF::NcUint64(),
               { dim_number_heat_blocks }, v_heat_node );
 }

 // Serialize sub-blocks

 for( Index i = 0; i < f_number_units; ++i ) {
  auto sub_block = get_unit_block( i );
  auto sub_group = group.addGroup( "UnitBlock" + std::to_string( i ) );
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

/*--------------------------------------------------------------------------*/
/*------------------------ End File UCBlock.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
