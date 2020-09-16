/*--------------------------------------------------------------------------*/
/*--------------------- File HydroUnitBlock.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the HydroUnitBlock class.
 *
 * \version 0.11
 *
 * \date 08 - 09 - 2020
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
 * Copyright &copy by Antonio Frangioni, Ali Ghezelsoflu
 */

/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <random>
#include "HydroUnitBlock.h"
#include "LinearFunction.h"
#include <map>
#include "FRowConstraint.h"
#include "UnitBlock.h"


/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;


/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register HydroUnitBlock to the Block factory


SMSpp_insert_in_factory_cpp_1( HydroUnitBlock );

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS OF HydroUnitBlock -----------------------*/
/*--------------------------------------------------------------------------*/

HydroUnitBlock::~HydroUnitBlock() {
 auto clear_constraints =
  []( boost::multi_array< FRowConstraint, 2 > & constraints ) {
   auto constraint = constraints.data();
   auto n = constraints.num_elements();
   for( decltype( n ) i = 0 ; i < n ; ++i , ++constraint )
    constraint->clear();
  };

 clear_constraints( MaxPowerPrimarySecondary_Const );
 clear_constraints( MinPowerPrimarySecondary_Const );
 clear_constraints( ActivePowerPrimary_Const );
 clear_constraints( ActivePowerSecondary_Const );
 clear_constraints( FlowActivePower_Const );
 clear_constraints( ActivePowerBounds_Const );
 clear_constraints( RampUp_Const );
 clear_constraints( RampDown_Const );
 clear_constraints( FlowRateBounds_Const );
 clear_constraints( FinalVolumeReservoir_Const );
 clear_constraints( VolumetricBounds_Const );
}

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void HydroUnitBlock::deserialize( netCDF::NcGroup & group ) {


#ifndef NDEBUG
 std::cerr << "[DEBUG] HydroUnitBlock::deserialize() - Checking Dims"
           << std::endl;
 std::vector< std::string > expected_dims = { "TimeHorizon",
                                              "NumberIntervals",
                                              "NumberReservoirs",
                                              "NumberArcs"  };
 check_dimensions( group, expected_dims, std::cerr );

 std::cerr << "[DEBUG] HydroUnitBlock::deserialize() - Checking Vars"
           << std::endl;
 std::vector< std::string > expected_vars = { "StartArc",
                                              "EndArc",
                                              "MinFlow",
                                              "MaxFlow",
                                              "MinVolumetric",
                                              "MaxVolumetric",
                                              "Inflows",
                                              "MinPower",
                                              "MaxPower",
                                              "DeltaRampUp",
                                              "DeltaRampDown",
                                              "PrimaryRho",
                                              "SecondaryRho",
                                              "NumberPieces",
                                              "LinearTerm",
                                              "ConstantTerm",
                                              "InertiaPower",
                                              "InitialFlowRate",
                                              "InitialVolumetric",
                                              "UphillFlow",
                                              "DownhillFlow"};
 check_variables( group, expected_vars, std::cerr );
#endif


 UnitBlock::deserialize_time_horizon( group );
 UnitBlock::deserialize_change_intervals( group );

 if( ! ::deserialize_dim( group, "NumberReservoirs", f_number_reservoirs, true ) )
  f_number_reservoirs = 1;

 if( ! ::deserialize_dim( group, "NumberArcs", f_number_arcs, true ) )
  f_number_arcs = 1;

 ::deserialize( group, "NumberPieces", f_number_arcs,
                v_number_pieces, true, true );

 if( ! ::deserialize_dim( group, "TotalNumberPieces", f_total_number_pieces, true ) ) {
  f_total_number_pieces = 0;
  for( const auto & n : v_number_pieces ) {
   f_total_number_pieces += n;
  }
 }
 f_total_number_pieces = f_total_number_pieces ?
                         f_total_number_pieces : f_number_arcs;

 ::deserialize( group, "StartArc", f_number_arcs, v_start_arc );
 ::deserialize( group, "EndArc", f_number_arcs, v_end_arc );

 ::deserialize( group, "Inflows", v_inflows, true, false );
 transpose( v_inflows );

 ::deserialize( group, "MinFlow", v_minimum_flow, true, true );
 transpose( v_minimum_flow );

 ::deserialize( group, "MaxFlow", v_maximum_flow, true, true );
 transpose( v_maximum_flow );


 ::deserialize( group, "MinPower", v_minimum_power, true, true );
 transpose( v_minimum_power );

 ::deserialize( group, "MaxPower", v_maximum_power, true, true );
 transpose( v_maximum_power );

 ::deserialize( group, "DeltaRampUp", v_delta_ramp_up, true, true );
 transpose( v_delta_ramp_up );

 ::deserialize( group, "DeltaRampDown", v_delta_ramp_down, true, true );
 transpose( v_delta_ramp_down );

 ::deserialize( group, "PrimaryRho", v_primary_rho, true, true );
 transpose( v_primary_rho );

 ::deserialize( group, "SecondaryRho", v_secondary_rho, true, true );
 transpose( v_secondary_rho );

 ::deserialize( group, "LinearTerm", f_total_number_pieces,
                v_linear_term, true, true );

 ::deserialize( group, "ConstantTerm", f_total_number_pieces,
                v_const_term, true, true );

 ::deserialize( group, "InertiaPower", v_inertia_power, true, true );

 ::deserialize( group, "InitialFlowRate", f_number_arcs,
                v_initial_flow_rate, true, true );

 ::deserialize( group, "InitialVolumetric", f_number_reservoirs,
                v_initial_volumetric, true, true );

 ::deserialize( group, "UphillFlow", f_number_arcs,
                v_uphill_delay, true, true );

 ::deserialize( group, "DownhillFlow", f_number_arcs,
                v_downhill_delay, true, true );

 ::deserialize( group, "MinVolumetric", v_minimum_volumetric, true, true );
 transpose( v_minimum_volumetric );
 ::deserialize( group, "MaxVolumetric", v_maximum_volumetric, true, true );
 transpose( v_maximum_volumetric );

 decompress_array( v_minimum_flow );
 decompress_array( v_maximum_flow );
 decompress_vol( v_minimum_volumetric );
 decompress_vol( v_maximum_volumetric );
 //decompress_vol( v_inflows );

 decompress_array( v_minimum_power );
 decompress_array( v_maximum_power );
 decompress_array( v_delta_ramp_up );
 decompress_array( v_delta_ramp_down );
 decompress_array( v_primary_rho );
 decompress_array( v_secondary_rho );
 decompress_array( v_inertia_power );

 if( v_linear_term.size() == 1 ) {
  v_linear_term.resize( f_number_arcs, v_linear_term[ 0 ] );
 }
 if( v_const_term.size() == 1 ) {
  v_const_term.resize( f_number_arcs, v_const_term[ 0 ] );
 }

 UnitBlock::deserialize( group );
}// end( HydroUnitBlock::deserialize )

/*--------------------------------------------------------------------------*/

void HydroUnitBlock::generate_abstract_variables( Configuration *stvv )
{
 if( AR & HasVar )
  return; // variables have already been generated

 UnitBlock::generate_abstract_variables( stvv );

 if( f_time_horizon == 0 ) {
  // there are no variables to be generated
  return;
 }

 v_volumetric.resize(boost::extents[f_number_reservoirs][ f_time_horizon]);
 for( Index g = 0; g < f_number_reservoirs; ++g )
  for( Index t = 0; t < f_time_horizon; ++t )
   v_volumetric[ g ][ t ].set_type( ColVariable::kNonNegative );
 add_static_variable ( v_volumetric, "vol" );

 v_flow_rate.resize(boost::extents[f_number_arcs][f_time_horizon]);
 v_active_power.resize(boost::extents[f_number_arcs][f_time_horizon]);
 v_primary_spinning_reserve.resize(boost::extents[f_number_arcs][f_time_horizon]);
 v_secondary_spinning_reserve.resize(boost::extents[f_number_arcs][f_time_horizon]);

 for( Index g = 0; g < f_number_arcs; ++g ) {
  for( Index t = 0; t < f_time_horizon; ++t ) {
   v_flow_rate[ g ][ t ].set_type( ColVariable::kContinuous );
   v_active_power[ g ][ t ].set_type( ColVariable::kContinuous );
   v_primary_spinning_reserve[ g ][ t ].set_type( ColVariable::kNonNegative );
   v_secondary_spinning_reserve[ g ][ t ].set_type( ColVariable::kNonNegative );
  }
 }

 add_static_variable ( v_flow_rate, "F" );
 add_static_variable ( v_active_power, "p" );
 add_static_variable ( v_primary_spinning_reserve, "pr" );
 add_static_variable ( v_secondary_spinning_reserve, "sr" );

 AR |= HasVar;
} // end( HydroUnitBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void HydroUnitBlock::generate_abstract_constraints( Configuration *stcc ) {

 if( AR & HasCst )
  return; // constraints have already been generated

 // final volumes fo each reservoir constraints

 assert( FinalVolumeReservoir_Const.empty());
 FinalVolumeReservoir_Const.resize
  ( boost::multi_array< FRowConstraint, 2 >::
    extent_gen()[f_time_horizon][ f_number_reservoirs ] );

 for( Index n = 0; n < f_number_reservoirs; ++n ) {

  auto l_f = new LinearFunction();

  for( Index l = 0; l < f_number_arcs; ++l ) {

   if( ! v_start_arc.empty() && ! v_end_arc.empty() ) {
    if( v_start_arc[l] == n && v_end_arc[l] <= f_number_reservoirs ) {
     if ( ! v_uphill_delay.empty() ) {
      if( v_uphill_delay[l] == 0 ) {
       l_f->add_variable( get_flow_rate( l , 0 ) , 1.0 );
      }
     }
     else {
      l_f->add_variable( get_flow_rate( l , 0 ) , 1.0 );
     }
    }
    if ( ! v_downhill_delay.empty() ) {
     if( v_downhill_delay[l] == 0 && v_end_arc[l] == n ) {
      l_f->add_variable( get_flow_rate( l , 0 ) , -1.0 );
     }
    }
   }
   else {
    l_f->add_variable( get_flow_rate( l , 0 ) , 1.0 );
   }
  }

  auto volumetric0 = get_volume( n , 0 );

  l_f->add_variable( volumetric0, 1.0 );
  FinalVolumeReservoir_Const[0][n].set_both( v_initial_volumetric[n] +
                                             v_inflows[n][0] );
  FinalVolumeReservoir_Const[0][n].set_function( l_f );

  for( Index t = 1, constraint_index = 1; t < f_time_horizon;
    ++t, ++constraint_index  ) {

   auto linear_function = new LinearFunction();

   for( Index l = 0; l < f_number_arcs; ++l ) {
    if( !v_start_arc.empty() && !v_end_arc.empty() ) {

     if( !v_uphill_delay.empty() ) {

      if( t - v_uphill_delay[l] >= 0 &&
          v_start_arc[l] == n &&
          v_end_arc[l] <= f_number_reservoirs ) {
       auto flow_rate = get_flow_rate( l , t - v_uphill_delay[l] );
       linear_function->add_variable( flow_rate , 1.0 );
      }
     }
     else {
      if( v_start_arc[l] == n && v_end_arc[l] <= f_number_reservoirs ) {
       auto flow_rate = get_flow_rate( l , t );
       linear_function->add_variable( flow_rate, 1.0 );
      }
     }
     if ( !v_downhill_delay.empty() ) {
      if( t - v_downhill_delay[l] >= 0 &&
          t - v_downhill_delay[l] <= f_time_horizon &&
          v_end_arc[l] == n ) {

       auto flow_rate = get_flow_rate( l , t - v_downhill_delay[l] );
       linear_function->add_variable( flow_rate , -1.0 );
      }
     }
     else {
      if( v_end_arc[l] == n ) {
       auto flow_rate = get_flow_rate( l , t );
       linear_function->add_variable( flow_rate , -1.0 );
      }
     }
    }
    else {
     auto flow_rate = get_flow_rate( l , t );
     linear_function->add_variable( flow_rate , 1.0 );
    }
   }

   auto volumetric_t = get_volume( n , t );
   auto volumetric_t_1 = get_volume( n , t-1 );

   linear_function->add_variable( volumetric_t, 1.0 );
   linear_function->add_variable( volumetric_t_1, -1.0 );

   FinalVolumeReservoir_Const[constraint_index][n].
    set_both( v_inflows[n][constraint_index] );

   FinalVolumeReservoir_Const[constraint_index][n].
    set_function( linear_function );
  }
 }
 add_static_constraint( FinalVolumeReservoir_Const, "FinalVolumeReservoir" );

 // maximum power output according to primary-secondary reserves constraints

 // Initial data check
 if ( ( ! v_minimum_power.empty() ) && ( ! v_maximum_power.empty() ) ) {
  for( Index arc = 0; arc < f_number_arcs; ++arc ) {
   for( Index t = 0; t < f_time_horizon; ++t ) {
    if( v_minimum_power[t][arc] > v_maximum_power[t][arc] ) {
     throw ( std::logic_error
             ( "HydroUnitBlock::maximum and minimum power output constraints: "
               "it must be that v_maximum_power >= v_minimum_power." ) );
    }
   }
  }
 }

 // Initial data check
 if ( ( ! v_minimum_flow.empty() ) && ( ! v_maximum_flow.empty() ) ) {
  for( Index arc = 0; arc < f_number_arcs; ++arc ) {
   for( Index t = 0; t < f_time_horizon; ++t ) {
    if( v_maximum_flow[t][arc] <= 0 &&
        v_minimum_flow[t][arc] < 0 &&
        v_number_pieces[arc] > 1 ) {
     throw ( std::logic_error
             ( "HydroUnitBlock::Data Error: it must be that for each pump"
               "when v_minimum_power < 0, then v_number_pieces == 1." ) );
    }
   }
  }
 }

 // Initial data check
 if ( ( ! v_minimum_flow.empty() ) && ( ! v_maximum_flow.empty() ) ) {
  for( Index arc = 0; arc < f_number_arcs; ++arc ) {
   for( Index t = 0; t < f_time_horizon; ++t ) {
    if( v_maximum_flow[t][arc] <= 0 &&
        v_minimum_flow[t][arc] < 0 &&
        !v_primary_rho.empty()) {
     throw ( std::logic_error
             ( "HydroUnitBlock::Data Error: it must be that for each pump"
               " v_primary_rho == 0." ));
    }
   }
  }
 }

 // Initial data check
 if ( ( ! v_minimum_flow.empty() ) && ( ! v_maximum_flow.empty() ) ) {
  for( Index arc = 0; arc < f_number_arcs; ++arc ) {
   for( Index t = 0; t < f_time_horizon; ++t ) {
    if( ( v_maximum_flow[t][arc] <= 0 ) &&
        ( v_minimum_flow[t][arc] < 0 ) &&
        ( ! v_secondary_rho.empty() ) ) {
     throw ( std::logic_error
             ( "HydroUnitBlock::Data Error: it must be that for each pump "
               " then v_secondary_rho == 0." ));
    }
   }
  }
 }

 if( ( ! v_primary_rho.empty() ) && ( ! v_secondary_rho.empty() ) ) {

  assert( MaxPowerPrimarySecondary_Const.empty());

  MaxPowerPrimarySecondary_Const.resize
   ( boost::multi_array< FRowConstraint, 2 >::
     extent_gen()[f_time_horizon][f_number_arcs] );

  for( Index arc = 0; arc < f_number_arcs; ++arc ) {
   for( Index t = 0; t < f_time_horizon; ++t ) {

    auto linear_function = new LinearFunction();

    auto active_power = get_active_power( arc , t );
    auto primary_spinning_reserve = get_primary_spinning_reserve( arc , t );
    auto secondary_spinning_reserve = get_secondary_spinning_reserve( arc , t );

    linear_function->add_variable( active_power, 1.0 );
    linear_function->add_variable( primary_spinning_reserve, 1.0 );
    linear_function->add_variable( secondary_spinning_reserve, 1.0 );

    MaxPowerPrimarySecondary_Const[t][arc].set_lhs( 0.0 );

    if( ! v_maximum_power.empty() ) {
     MaxPowerPrimarySecondary_Const[t][arc].set_rhs( v_maximum_power[t][arc] );
    } else {
     MaxPowerPrimarySecondary_Const[t][arc].set_rhs( 0.0 );
    }
    MaxPowerPrimarySecondary_Const[t][arc].set_function( linear_function );
   }
  }
  add_static_constraint( MaxPowerPrimarySecondary_Const ,
                         "MaxPowerPrimarySecondary" );


  // minimum power output according to primary-secondary reserves constraints

  assert( MinPowerPrimarySecondary_Const.empty());

  MinPowerPrimarySecondary_Const.resize
   ( boost::multi_array< FRowConstraint, 2 >::
     extent_gen()[f_time_horizon][f_number_arcs] );

  for( Index arc = 0; arc < f_number_arcs; ++arc ) {
   for( Index t = 0; t < f_time_horizon; ++t ) {

    auto linear_function = new LinearFunction();

    auto active_power = get_active_power( arc , t );
    auto primary_spinning_reserve = get_primary_spinning_reserve( arc , t );
    auto secondary_spinning_reserve = get_secondary_spinning_reserve( arc , t );

    linear_function->add_variable( active_power, 1.0 );
    linear_function->add_variable( primary_spinning_reserve, -1.0 );
    linear_function->add_variable( secondary_spinning_reserve, -1.0 );

    if( !v_minimum_power.empty()) {
     MinPowerPrimarySecondary_Const[t][arc].set_lhs( v_minimum_power[t][arc] );
    } else {
     MinPowerPrimarySecondary_Const[t][arc].set_lhs( 0.0 );
    }
    MinPowerPrimarySecondary_Const[t][arc].set_rhs( Inf< double >());
    MinPowerPrimarySecondary_Const[t][arc].set_function( linear_function );
   }
  }

  add_static_constraint( MinPowerPrimarySecondary_Const ,
                         "MinPowerPrimarySecondary" );
 } else {

  assert( ActivePowerBounds_Const.empty());

  ActivePowerBounds_Const.resize
   ( boost::multi_array< FRowConstraint, 2 >::
     extent_gen()[f_time_horizon][f_number_arcs] );

  for( Index arc = 0; arc < f_number_arcs; ++arc ) {
   for( Index t = 0; t < f_time_horizon; ++t ) {

    auto linear_function = new LinearFunction();
    auto active_power = get_active_power( arc , t );
    linear_function->add_variable( active_power, 1.0 );

    if( !v_maximum_power.empty()) {

     ActivePowerBounds_Const[t][arc].set_lhs( v_minimum_power[t][arc] );
    } else {
     ActivePowerBounds_Const[t][arc].set_lhs( 0.0 );
    }
    if( !v_maximum_power.empty()) {
     ActivePowerBounds_Const[t][arc].set_rhs( v_maximum_power[t][arc] );
    } else {
     ActivePowerBounds_Const[t][arc].set_rhs( 0.0 );

    }
    ActivePowerBounds_Const[t][arc].set_function( linear_function );
   }
  }
  add_static_constraint( ActivePowerBounds_Const, "ActivePowerBounds" );
 }
 // power output relation with to primary reserves constraints

  if( ! v_primary_rho.empty() ) {

   assert( ActivePowerPrimary_Const.empty());

   ActivePowerPrimary_Const.resize
    ( boost::multi_array< FRowConstraint, 2 >::
      extent_gen()[f_time_horizon][f_number_arcs] );

   for( Index arc = 0; arc < f_number_arcs; ++arc ) {
    for( Index t = 0; t < f_time_horizon; ++t ) {

     if ( v_minimum_flow[t][arc] >= 0 && v_maximum_flow[t][arc] > 0) { //Turbines
      auto linear_func = new LinearFunction();
      auto active_power = get_active_power( arc , t );
      auto primary_spinning_reserve = get_primary_spinning_reserve( arc , t );

      linear_func->add_variable( active_power, v_primary_rho[t][arc] );
      linear_func->add_variable( primary_spinning_reserve, -1.0 );

      ActivePowerPrimary_Const[t][arc].set_lhs( 0.0 );
      ActivePowerPrimary_Const[t][arc].set_rhs( Inf< double >());
      ActivePowerPrimary_Const[t][arc].set_function( linear_func );
     }

     if( ( v_maximum_flow[t][arc] <= 0 ) &&
         ( v_minimum_flow[t][arc] < 0 ) ) { //Pumps
      auto linear_f = new LinearFunction();
      auto primary_spinning_reserve = get_primary_spinning_reserve( arc , t );

      linear_f->add_variable( primary_spinning_reserve, 1.0 );
      ActivePowerPrimary_Const[t][arc].set_both( 0.0 );
      ActivePowerPrimary_Const[t][arc].set_function( linear_f );
     }

     if( ( v_maximum_flow[t][arc] == 0 ) &&
         ( v_minimum_flow[t][arc] == 0 ) ) { //Nothing
      auto linear_function = new LinearFunction();
      auto flow_rate = get_flow_rate( arc , t );
      linear_function->add_variable( flow_rate, 1.0 );
      ActivePowerPrimary_Const[t][arc].set_both( 0.0 );
      ActivePowerPrimary_Const[t][arc].set_function( linear_function );
     }
    }
   }
   add_static_constraint( ActivePowerPrimary_Const, "ActivePowerPrimary" );

  }
  else {
   assert( ActivePowerPrimary_Const.empty());

   ActivePowerPrimary_Const.resize
    ( boost::multi_array< FRowConstraint, 2 >::
      extent_gen()[f_time_horizon][f_number_arcs] );

   for( Index arc = 0; arc < f_number_arcs; ++arc ) {
    for( Index t = 0; t < f_time_horizon; ++t ) {
     auto linear_function = new LinearFunction();
     auto primary_spinning_reserve = get_primary_spinning_reserve( arc , t );

     linear_function->add_variable( primary_spinning_reserve, 1.0 );
     ActivePowerPrimary_Const[t][arc].set_both( 0.0 );
     ActivePowerPrimary_Const[t][arc].set_function( linear_function );
    }
   }
   add_static_constraint( ActivePowerPrimary_Const, "ActivePowerPrimary" );
  }

 // power output relation with to secondary reserves constraints
  if( ! v_secondary_rho.empty() ) {

   assert( ActivePowerSecondary_Const.empty());

   ActivePowerSecondary_Const.resize
    ( boost::multi_array< FRowConstraint, 2 >::
      extent_gen()[f_time_horizon][f_number_arcs] );

   for( Index arc = 0; arc < f_number_arcs; ++arc ) {
    for( Index t = 0; t < f_time_horizon; ++t ) {
     if ( ( v_minimum_flow[t][arc] >= 0 ) &&
          ( v_maximum_flow[t][arc] > 0 ) ) { //Turbines

      auto linear_function = new LinearFunction();
      auto active_power = get_active_power( arc , t );
      auto secondary_spinning_reserve =
       get_secondary_spinning_reserve( arc , t );

      linear_function->add_variable( active_power, v_secondary_rho[t][arc] );
      linear_function->add_variable( secondary_spinning_reserve, -1.0 );
      ActivePowerSecondary_Const[t][arc].set_lhs( 0.0 );
      ActivePowerSecondary_Const[t][arc].set_rhs( Inf< double >());
      ActivePowerSecondary_Const[t][arc].set_function( linear_function );
     }

     if( ( v_maximum_flow[t][arc] <= 0 ) &&
         ( v_minimum_flow[t][arc] < 0 ) ) { //Pumps
      auto linear_f = new LinearFunction();
      auto secondary_spinning_reserve =
       get_secondary_spinning_reserve( arc , t );

      linear_f->add_variable( secondary_spinning_reserve, 1.0 );
      ActivePowerSecondary_Const[t][arc].set_both( 0.0 );
      ActivePowerSecondary_Const[t][arc].set_function( linear_f );
     }

     if( ( v_maximum_flow[t][arc] == 0 ) &&
         ( v_minimum_flow[t][arc] == 0 ) ) { //Nothing
      auto l_function = new LinearFunction();
      auto flow_rate = get_flow_rate( arc , t );
      l_function->add_variable( flow_rate, 1.0 );
      ActivePowerSecondary_Const[t][arc].set_both( 0.0 );
      ActivePowerSecondary_Const[t][arc].set_function( l_function );
     }
    }
   }
   add_static_constraint( ActivePowerSecondary_Const, "ActivePowerSecondary" );
  }
  else {
   assert( ActivePowerSecondary_Const.empty());

   ActivePowerSecondary_Const.resize
    ( boost::multi_array< FRowConstraint, 2 >::
      extent_gen()[f_time_horizon][f_number_arcs] );

   for( Index arc = 0; arc < f_number_arcs; ++arc ) {
    for( Index t = 0; t < f_time_horizon; ++t ) {
     auto linear_function = new LinearFunction();
     auto secondary_spinning_reserve =
      get_secondary_spinning_reserve( arc , t );

     linear_function->add_variable( secondary_spinning_reserve, 1.0 );
     ActivePowerSecondary_Const[t][arc].set_both( 0.0 );
     ActivePowerSecondary_Const[t][arc].set_function( linear_function );
    }
   }
   add_static_constraint( ActivePowerSecondary_Const, "ActivePowerSecondary" );
  }

  assert( FlowActivePower_Const.empty());

  FlowActivePower_Const.resize
   ( boost::multi_array< FRowConstraint, 2 >::
     extent_gen()[f_time_horizon][f_total_number_pieces] );

 if( f_number_arcs > 0 ) {

  if( ( ! v_minimum_flow.empty() ) && ( ! v_maximum_flow.empty() ) ) {

   for( Index t = 0; t < f_time_horizon; ++t ) {
    Index piece = 0;
    Index constraint_index = 0;
    Index end = 0;
    for( Index arc = 0; arc < f_number_arcs; ++arc ) {
     if( ! v_number_pieces.empty() ) {
      end += v_number_pieces[arc];
     }

     if( ( v_minimum_flow[t][arc] >= 0 ) &&
         ( v_maximum_flow[t][arc] > 0 ) ) { //Turbines
      for( ; piece < end; ++piece ) {

       auto linear_function_turbine = new LinearFunction();

       auto active_power = get_active_power( arc , t );
       auto flow_rate = get_flow_rate( arc , t );

       linear_function_turbine->add_variable( active_power, 1.0 );

       if( ! v_linear_term.empty() ) {
        linear_function_turbine->add_variable( flow_rate,
                                               - v_linear_term[piece] );
       }
       else {
        linear_function_turbine->add_variable( flow_rate, 0.0 );
       }
       if( ! v_const_term.empty() ) {
        FlowActivePower_Const[t][piece].set_rhs( v_const_term[piece] );
       }
       else {
        FlowActivePower_Const[t][piece].set_rhs( 0.0 );
       }

       FlowActivePower_Const[t][piece].set_lhs( -Inf< double >() );
       FlowActivePower_Const[t][piece].set_function( linear_function_turbine );
       ++constraint_index;
      }
     }

     if( ( v_maximum_flow[t][arc] <= 0 ) &&
         ( v_minimum_flow[t][arc] < 0 ) ) { //Pumps

      auto linear_function_pump = new LinearFunction();
      auto active_power = get_active_power( arc , t );
      auto flow_rate = get_flow_rate( arc , t );

      linear_function_pump->add_variable( active_power, 1.0 );
      linear_function_pump->add_variable( flow_rate,
                                          - v_linear_term[constraint_index] );
      FlowActivePower_Const[t][constraint_index].set_both( 0.0 );
      FlowActivePower_Const[t][constraint_index].
       set_function( linear_function_pump );
      ++constraint_index;
      piece = constraint_index;
     }

     if( ( v_maximum_flow[t][arc] == 0 ) &&
         ( v_minimum_flow[t][arc] == 0 ) ) { //Nothing

      auto linear_function_nothing = new LinearFunction();
      auto flow_rate = get_flow_rate( arc , t );

      linear_function_nothing->add_variable( flow_rate, 1.0 );
      FlowActivePower_Const[t][constraint_index].set_both( 0.0 );
      FlowActivePower_Const[t][constraint_index].
       set_function( linear_function_nothing );
      ++constraint_index;
      piece = constraint_index;
     }
    }
   }
  }

  if ( ( v_minimum_flow.empty() ) && ( ! v_maximum_flow.empty() ) ) {

   for( Index t = 0; t < f_time_horizon; ++t ) {
    Index piece = 0;
    Index constraint_index = 0;
    Index end = 0;
    for( Index arc = 0; arc < f_number_arcs; ++arc ) {
     if( ! v_number_pieces.empty() ) {
      end += v_number_pieces[arc];
     }

     if( v_maximum_flow[t][arc] > 0 ) { //Turbines
      for( ; piece <= end; ++piece ) {

       auto linear_function_turbine = new LinearFunction();
       auto active_power = get_active_power( arc , t );
       auto flow_rate = get_flow_rate( arc , t );

       linear_function_turbine->add_variable( active_power, 1.0 );

       if( !v_linear_term.empty()) {
        linear_function_turbine->add_variable( flow_rate ,
                                              - v_linear_term[piece] );
       }
       else {
        linear_function_turbine->add_variable( flow_rate, 0.0 );
       }
       if( !v_const_term.empty()) {
        FlowActivePower_Const[t][piece].set_rhs( v_const_term[piece] );
       } else {
        FlowActivePower_Const[t][piece].set_rhs( 0.0 );
       }
       FlowActivePower_Const[t][piece].set_lhs( -Inf< double >() );
       FlowActivePower_Const[t][piece].set_function( linear_function_turbine );
       ++constraint_index;
      }
     }

     if( v_maximum_flow[t][arc] == 0 ) { //Nothing

      auto linear_function_nothing = new LinearFunction();
      auto flow_rate = get_flow_rate( arc , t );

      linear_function_nothing->add_variable( flow_rate, 1.0 );
      FlowActivePower_Const[t][constraint_index].set_both( 0.0 );
      FlowActivePower_Const[t][constraint_index].
       set_function( linear_function_nothing );
      ++constraint_index;
      piece = constraint_index;
     }
    }
   }
  }
  add_static_constraint( FlowActivePower_Const, "FlowActivePower" );
 }

 {
  // Initial data check
  if( ( ! v_minimum_flow.empty() ) && ( ! v_maximum_flow.empty() ) ) {
   for( Index arc = 0; arc < f_number_arcs; ++arc ) {
    for( Index t = 0; t < f_time_horizon; ++t ) {
     if( v_minimum_flow[t][arc] > v_maximum_flow[t][arc] ) {
      throw ( std::logic_error
              ( "HydroUnitBlock::flow rate variable bounds: "
                "it must be that v_maximum_flow >= v_minimum_flow." ));
     }
    }
   }
  }
  if( ( ! v_minimum_flow.empty() ) && ( ! v_maximum_flow.empty() ) ) {

   assert( FlowRateBounds_Const.empty());

   FlowRateBounds_Const.resize
    ( boost::multi_array< FRowConstraint, 2 >::
      extent_gen()[f_time_horizon][f_number_arcs] );

   for( Index t = 0; t < f_time_horizon; ++t ) {
    for( Index arc = 0; arc < f_number_arcs; ++arc ) {
     auto linear_function = new LinearFunction();
     auto flow_rate = get_flow_rate( arc , t );
     linear_function->add_variable( flow_rate, 1.0 );
     FlowRateBounds_Const[t][arc].set_lhs( v_minimum_flow[t][arc] );
     FlowRateBounds_Const[t][arc].set_rhs( v_maximum_flow[t][arc] );
     FlowRateBounds_Const[t][arc].set_function( linear_function );
    }
   }

   add_static_constraint( FlowRateBounds_Const, "FlowRateBounds" );
  }
 }

 // ramp-up constraints
 if( ! v_delta_ramp_up.empty() ) {

  assert( RampUp_Const.empty());

  RampUp_Const.resize
   ( boost::multi_array< FRowConstraint, 2 >::
     extent_gen()[f_time_horizon][f_number_arcs] );

  // Initial condition

  for( Index arc = 0; arc < f_number_arcs; ++arc ) {
   auto linear_f = new LinearFunction();
   auto flow_rate = get_flow_rate( arc , 0 );
   linear_f->add_variable( flow_rate , 1.0 );

   RampUp_Const[0][arc].set_lhs( -Inf< double >());
   RampUp_Const[0][arc].set_rhs( v_delta_ramp_up[0][arc] +
                                 get_initial_flow_rate( arc ) );
   RampUp_Const[0][arc].set_function( linear_f );


   for( Index t = 1, constraint_index = 1; t < f_time_horizon;
        ++t, ++constraint_index  ) {

    auto linear_function = new LinearFunction();
    auto flow_rate_t = get_flow_rate( arc , t );
    auto flow_rate_t_1 = get_flow_rate( arc , t - 1 );

    linear_function->add_variable( flow_rate_t, 1.0 );
    linear_function->add_variable( flow_rate_t_1, -1.0 );

    RampUp_Const[constraint_index][arc].set_lhs( -Inf< double >());
    RampUp_Const[constraint_index][arc].set_rhs( v_delta_ramp_up[t][arc] );
    RampUp_Const[constraint_index][arc].set_function( linear_function );
   }
  }

  add_static_constraint( RampUp_Const,"RampUp");
 }


 // ramp-down constraints
 if( ! v_delta_ramp_down.empty() ) {

  assert( RampDown_Const.empty() );

  RampDown_Const.resize
   ( boost::multi_array< FRowConstraint, 2 >::
     extent_gen()[f_time_horizon][f_number_arcs] );

  // Initial condition

  for( Index arc = 0; arc < f_number_arcs; ++arc ) {

   auto linear_f = new LinearFunction();
   auto flow_rate0 = get_flow_rate( arc , 0 );

   linear_f->add_variable( flow_rate0, 1.0 );

   RampDown_Const[0][arc].set_lhs( get_initial_flow_rate( arc ) -
                                   v_delta_ramp_down[0][arc] );
   RampDown_Const[0][arc].set_rhs( Inf< double >());
   RampDown_Const[0][arc].set_function( linear_f );

   for( Index t = 1, constraint_index = 1; t < f_time_horizon;
        ++t, ++constraint_index  ) {

    auto linear_function = new LinearFunction();
    auto flow_rate_t = get_flow_rate( arc , t );
    auto flow_rate_t_1 = get_flow_rate( arc , t - 1 );

    linear_function->add_variable( flow_rate_t_1, 1.0 );
    linear_function->add_variable( flow_rate_t, -1.0 );

    RampDown_Const[constraint_index][arc].set_lhs( -Inf< double >());
    RampDown_Const[constraint_index][arc].set_rhs( v_delta_ramp_down[t][arc] );
    RampDown_Const[constraint_index][arc].set_function( linear_function );

   }
  }

  add_static_constraint( RampDown_Const, "RampDown");

 }

 // volumetric bounds constraints

 // Initial data check
 if ( ( ! v_minimum_volumetric.empty() ) &&
      ( ! v_maximum_volumetric.empty() ) ) {
  for( Index node = 0; node < f_number_reservoirs; ++node ) {
   for( Index t = 0; t < f_time_horizon; ++t ) {
    if( ( v_minimum_volumetric[node][t] > v_maximum_volumetric[node][t] ) ||
        ( v_minimum_volumetric[node][t] < 0 ) ||
        ( v_maximum_volumetric[node][t] < 0 ) ) {
     throw ( std::logic_error
             ( "HydroUnitBlock::Volumetric Bounds Constraint: "
               "it must be that 0 <= MinV[ r , t ] <= MaxV[ r , t ] ." ) );
    }
   }
  }
 }

 if( ( ! v_minimum_volumetric.empty() ) &&
     ( ! v_maximum_volumetric.empty() ) ) {

  assert( VolumetricBounds_Const.empty());

  VolumetricBounds_Const.resize
   ( boost::multi_array< FRowConstraint, 2 >::
     extent_gen()[f_number_reservoirs][f_time_horizon] );

  for( Index node = 0; node < f_number_reservoirs; ++node ) {
   for( Index t = 0; t < f_time_horizon; ++t ) {

    auto linear_function = new LinearFunction();

    auto volumetric = get_volume( node , t );
    linear_function->add_variable( volumetric, 1.0 );

    VolumetricBounds_Const[node][t].set_lhs( v_minimum_volumetric[node][t] );
    VolumetricBounds_Const[node][t].set_rhs( v_maximum_volumetric[node][t] );
    VolumetricBounds_Const[node][t].set_function( linear_function );
   }
  }

  add_static_constraint( VolumetricBounds_Const, "VolumetricBounds" );
 }

 AR |= HasCst;
} // end( HydroUnitBlock::generate_abstract_constraints )


/*--------------------------------------------------------------------------*/
/*-------- METHODS FOR LOADING, PRINTING & SAVING THE HydroUnitBlock -------*/
/*--------------------------------------------------------------------------*/

void HydroUnitBlock::serialize( netCDF::NcGroup & group ) const {

 UnitBlock::serialize( group );

 auto dim_time_horizon = group.getDim( "TimeHorizon" );
 auto NumberIntervals = group.getDim( "NumberIntervals" );


 auto dim_total_number_pieces = group.addDim( "TotalNumberPieces",
                                              f_total_number_pieces );

 auto dim_number_reservoirs = group.addDim( "NumberReservoirs",
                                            f_number_reservoirs
                                            ? f_number_reservoirs : 1 );

 auto dim_number_arcs = group.addDim( "NumberArcs",
                                      f_number_arcs ? f_number_arcs : 1 );


 ::serialize( group, "NumberPieces", netCDF::NcUint(),
              dim_number_arcs, v_number_pieces, true );

 ::serialize( group, "StartLine", netCDF::NcUint(),
              dim_number_reservoirs, v_start_arc, false );

 ::serialize( group, "EndLine", netCDF::NcUint(),
              dim_number_reservoirs, v_end_arc, false );


 if( !v_minimum_flow.empty() ) {

  ::serialize( group, "MinFlow", netCDF::NcDouble(),
               { NumberIntervals, dim_number_arcs },
               v_minimum_flow, true );
 }

 if( !v_maximum_flow.empty() ) {

  ::serialize( group, "MaxFlow", netCDF::NcDouble(),
               { NumberIntervals, dim_number_arcs },
               v_maximum_flow, true );
 }

 if( !v_minimum_volumetric.empty() ) {

  ::serialize( group, "MinVolumetric", netCDF::NcDouble(),
               { dim_number_reservoirs, NumberIntervals },
               v_minimum_volumetric, true );
 }

 if( !v_maximum_volumetric.empty() ) {

  ::serialize( group, "MaxVolumetric", netCDF::NcDouble(),
               { dim_number_reservoirs, NumberIntervals },
               v_maximum_volumetric, true );
 }

 ::serialize( group, "Inflows", netCDF::NcDouble(),
              { dim_number_reservoirs, dim_time_horizon },
              v_inflows, false );

 ::serialize( group, "MinPower", netCDF::NcDouble(),
              { NumberIntervals, dim_number_arcs },
              v_minimum_power, true );

 ::serialize( group, "MaxPower", netCDF::NcDouble(),
              { NumberIntervals, dim_number_arcs },
              v_maximum_power, true );

 ::serialize( group, "DeltaRampUp", netCDF::NcDouble(),
              { NumberIntervals, dim_number_arcs },
              v_delta_ramp_up, true );

 ::serialize( group, "DeltaRampDown", netCDF::NcDouble(),
              { NumberIntervals, dim_number_arcs },
              v_delta_ramp_down, true );

 if( !v_primary_rho.empty() ) {
  ::serialize( group, "PrimaryRho", netCDF::NcDouble(),
               { NumberIntervals, dim_number_arcs },
               v_primary_rho, true );
 }

 if( !v_secondary_rho.empty() ) {

  ::serialize( group, "SecondaryRho", netCDF::NcDouble(),
               { NumberIntervals, dim_number_arcs },
               v_secondary_rho, true );
 }

 ::serialize( group, "LinearTerm", netCDF::NcDouble(),
              dim_total_number_pieces, v_linear_term, false );

 ::serialize( group, "ConstantTerm", netCDF::NcDouble(),
              dim_total_number_pieces, v_const_term, false );


 ::serialize( group, "InertiaPower", netCDF::NcDouble(),
              { NumberIntervals, dim_number_arcs },
              v_inertia_power, true );

 ::serialize( group, "InitialFlowRate", netCDF::NcDouble(),
              dim_number_arcs, v_initial_flow_rate, false );

 ::serialize( group, "InitialVolumetric", netCDF::NcDouble(),
              dim_number_reservoirs, v_initial_volumetric, false );

 ::serialize( group, "UphillFlow", netCDF::NcUint(),
              dim_number_arcs, v_uphill_delay, true );

 ::serialize( group, "DownhillFlow", netCDF::NcUint(),
              dim_number_arcs, v_downhill_delay, true );
}  // end( HydroUnitBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

void
HydroUnitBlock::set_inflow( std::vector< double >::const_iterator values,
                            Block::Subset && subset,
                            const bool ordered,
                            c_ModParam issuePMod,
                            c_ModParam issueAMod ) {
 if( subset.empty() ) {
  return;
 }

 if( v_inflows.empty() ) {
  if( std::all_of( values,
                   values + subset.size(),
                   []( double cst ) {
                    return ( cst == 0 );
                   } ) ) {
   return;
  }

  v_inflows.resize( boost::extents[ f_number_reservoirs ][ f_time_horizon ] );
 }

 // If nothing changes, return
 bool identical = true;
 for( auto i : subset ) {
  if( i >= v_inflows.size() ) {
   throw ( std::invalid_argument( "invalid value in subset" ) );
  }
  if( *( v_inflows.data() + i ) != *( values++ ) ) {
   identical = false;
  }
 }
 if( identical ) {
  return;
 }
 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  for( auto i : subset ) {
   Index t = i % f_time_horizon;
   Index r = i / f_time_horizon;
   v_inflows[ r ][ t ] = *( values++ );
   // *( v_inflows.data() + i ) = *( values++ );
  }

  if( AR & HasCst ) {
   // Change the abstract representation

   for( auto i : subset ) {
    Index t = i % f_time_horizon;
    Index r = i / f_time_horizon;

    if( t == 0 ) {
     FinalVolumeReservoir_Const[ t ][ r ]
      .set_both( v_initial_volumetric[ r ] + v_inflows[ r ][ t ], issueAMod );
    } else {
     FinalVolumeReservoir_Const[ t ][ r ]
      .set_both( v_inflows[ r ][ t ], issueAMod );
    }
   }
  }
 }

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification
  if( !ordered ) {
   std::sort( subset.begin(), subset.end() );
  }

  Block::add_Modification(
   std::make_shared< HydroUnitBlockSbstMod >( this,
                                              HydroUnitBlockMod::eSetInf,
                                              std::move( subset ) ),
   Observer::par2chnl( issuePMod ) );
 }
}

/*--------------------------------------------------------------------------*/

void
HydroUnitBlock::set_inflow( std::vector< double >::const_iterator values,
                            Block::Range rng,
                            c_ModParam issuePMod,
                            c_ModParam issueAMod ) {
 rng.second = std::min( rng.second,
                        get_time_horizon() * get_number_reservoirs() );
 if( rng.second <= rng.first ) {
  return;
 }

 if( v_inflows.empty() ) {
  if( std::all_of( values,
                   values + ( rng.second - rng.first ),
                   []( double cst ) {
                    return ( cst == 0 );
                   } ) ) {
   return;
  }

  v_inflows.resize( boost::extents[ f_number_reservoirs ][ f_time_horizon ] );
 }

 // If nothing changes, return
 if( std::equal( values,
                 values + ( rng.second - rng.first ),
                 v_inflows.data() + rng.first ) ) {
  return;
 }
 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  std::copy( values,
             values + ( rng.second - rng.first ),
             v_inflows.data() + rng.first );

  if( AR & HasCst ) {
   // Change the abstract representation

   for( Index i = rng.first; i < rng.second; ++i ) {
    Index t = i % f_time_horizon;
    Index r = i / f_time_horizon;

    if( t == 0 ) {
     FinalVolumeReservoir_Const[ t ][ r ]
      .set_both( v_initial_volumetric[ r ] + v_inflows[ r ][ t ], issueAMod );
    } else {
     FinalVolumeReservoir_Const[ t ][ r ]
      .set_both( v_inflows[ r ][ t ], issueAMod );
    }
   }
  }
 }

 if( issue_pmod( issuePMod ) ) {
  Block::add_Modification(
   std::make_shared< HydroUnitBlockRngdMod >( this,
                                              HydroUnitBlockMod::eSetInf,
                                              rng ),
   Observer::par2chnl( issuePMod ) );
 }
}

/*--------------------------------------------------------------------------*/

void
HydroUnitBlock::set_inertia_power( std::vector< double >::const_iterator values,
                                   Subset && subset,
                                   const bool ordered,
                                   c_ModParam issuePMod,
                                   c_ModParam issueAMod ) {
 if( subset.empty() ) {
  return;
 }

 if( v_inertia_power.empty() ) {
  if( std::all_of( values,
                   values + subset.size(),
                   []( double cst ) {
                    return ( cst == 0 );
                   } ) ) {
   return;
  }

  v_inertia_power.resize( boost::extents[ f_time_horizon ][ f_number_arcs ] );
 }

 // If nothing changes, return
 bool identical = true;
 for( auto i : subset ) {
  if( i >= v_inertia_power.size() ) {
   throw ( std::invalid_argument( "invalid value in subset" ) );
  }
  if( *( v_inertia_power.data() + i ) != *( values++ ) ) {
   identical = false;
  }
 }
 if( identical ) {
  return;
 }
 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  for( auto i : subset ) {
   Index a = i % f_number_arcs;
   Index t = i / f_number_arcs;
   v_inertia_power[ t ][ a ] = *( values++ );
  }

  if( AR & HasCst ) {
   // Change the abstract representation
   // FIXME: v_inertia_power is not used
  }
 }

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification
  if( !ordered ) {
   std::sort( subset.begin(), subset.end() );
  }

  Block::add_Modification(
   std::make_shared< HydroUnitBlockSbstMod >( this,
                                              HydroUnitBlockMod::eSetInerP,
                                              std::move( subset ) ),
   Observer::par2chnl( issuePMod ) );
 }
}

/*--------------------------------------------------------------------------*/

void
HydroUnitBlock::set_inertia_power( std::vector< double >::const_iterator values,
                                   Block::Range rng,
                                   c_ModParam issuePMod,
                                   c_ModParam issueAMod ) {
 rng.second = std::min( rng.second, f_number_arcs * get_time_horizon() );
 if( rng.second <= rng.first ) {
  return;
 }

 if( v_inertia_power.empty() ) {
  if( std::all_of( values,
                   values + ( rng.second - rng.first ),
                   []( double cst ) {
                    return ( cst == 0 );
                   } ) ) {
   return;
  }

  v_inertia_power.resize( boost::extents[ f_time_horizon ][ f_number_arcs ] );
 }

 // If nothing changes, return
 if( std::equal( values,
                 values + ( rng.second - rng.first ),
                 v_inertia_power.data() + rng.first ) ) {
  return;
 }
 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  std::copy( values,
             values + ( rng.second - rng.first ),
             v_inertia_power.data() + rng.first );

  if( AR & HasCst ) {
   // Change the abstract representation
   // FIXME: v_inertia_power is not used

  }
 }

 if( issue_pmod( issuePMod ) ) {
  Block::add_Modification(
   std::make_shared< HydroUnitBlockRngdMod >( this,
                                              HydroUnitBlockMod::eSetInerP,
                                              rng ),
   Observer::par2chnl( issuePMod ) );
 }
}

/*--------------------------------------------------------------------------*/

void
HydroUnitBlock::set_initial_volumetric(
 std::vector< double >::const_iterator values,
 Block::Subset && subset,
 const bool ordered,
 c_ModParam issuePMod,
 c_ModParam issueAMod ) {

 if( subset.empty() ) {
  return;
 }

 if( v_initial_volumetric.empty() ) {
  if( std::all_of( values,
                   values + subset.size(),
                   []( double cst ) {
                    return ( cst == 0 );
                   } ) ) {
   return;
  }

  Index max_index = *max_element( std::begin( subset ), std::end( subset ) );
  v_initial_volumetric.assign( max_index, 0 );
 }

 // If nothing changes, return
 bool identical = true;
 auto temp_values = values;
 for( auto i : subset ) {
  if( i >= v_initial_volumetric.size() ) {
   throw ( std::invalid_argument( "invalid value in subset" ) );
  }
  if( v_initial_volumetric[ i ] != *( temp_values++ ) ) {
   identical = false;
  }
 }
 if( identical ) {
  return;
 }

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  temp_values = values;
  for( auto i : subset ) {
   v_initial_volumetric[ i ] = *( temp_values++ );
  }

  if( not_dry_run( issueAMod ) && ( AR & HasObj ) ) {
   // Change the abstract representation
   for( auto i : subset ) {
    Index t = i % f_time_horizon;
    Index r = i / f_time_horizon;

    if( t == 0 ) {
     FinalVolumeReservoir_Const[ t ][ r ]
      .set_both( v_initial_volumetric[ r ] + v_inflows[ r ][ t ], issueAMod );
    }
   }
  }
 }

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification
  if( !ordered ) {
   std::sort( subset.begin(), subset.end() );
  }
  Block::add_Modification(
   std::make_shared< HydroUnitBlockSbstMod >( this,
                                              HydroUnitBlockMod::eSetInitV,
                                              std::move( subset ) ),
   Observer::par2chnl( issuePMod ) );
 }
}

/*--------------------------------------------------------------------------*/

void
HydroUnitBlock::set_initial_volumetric(
 std::vector< double >::const_iterator values,
 Block::Range rng,
 c_ModParam issuePMod,
 c_ModParam issueAMod ) {

 rng.second = std::min( rng.second, f_time_horizon );
 if( rng.second <= rng.first ) {
  return;
 }

 if( v_initial_volumetric.empty() ) {
  if( std::all_of( values,
                   values + ( rng.second - rng.first ),
                   []( double cst ) {
                    return ( cst == 0 );
                   } ) ) {
   return;
  }

  Index max_index = rng.second;
  v_initial_volumetric.assign( max_index, 0 );
 }

 // If nothing changes, return
 if( std::equal( values,
                 values + ( rng.second - rng.first ),
                 v_initial_volumetric.begin() + rng.first ) ) {
  return;
 }

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  std::copy( values,
             values + ( rng.second - rng.first ),
             v_initial_volumetric.begin() + rng.first );

  if( not_dry_run( issueAMod ) && ( AR & HasCst ) ) {
   // Change the abstract representation
   for( Index i = rng.first; i < rng.second; ++i ) {
    Index t = i % f_time_horizon;
    Index r = i / f_time_horizon;

    if( t == 0 ) {
     FinalVolumeReservoir_Const[ t ][ r ]
      .set_both( v_initial_volumetric[ r ] + v_inflows[ r ][ t ], issueAMod );
    }
   }
  }
 }

 if( issue_pmod( issuePMod ) ) {
  Block::add_Modification(
   std::make_shared< HydroUnitBlockRngdMod >( this,
                                              HydroUnitBlockMod::eSetInitV,
                                              rng ),
   Observer::par2chnl( issuePMod ) );
 }
}

/*--------------------------------------------------------------------------*/

template< typename T >
void HydroUnitBlock::transpose( boost::multi_array< T, 2 > & a ) {
 long rows = a.shape()[ 0 ];
 long cols = a.shape()[ 1 ];
 if( rows > 1 && cols == 1 ) {
  // The vector must be transposed
  boost::array< typename boost::multi_array< T, 2 >::index, 2 >
   dims = { { 1, rows } };
  a.reshape( dims );
 }
}

/*--------------------------------------------------------------------------*/

void HydroUnitBlock::decompress_array( boost::multi_array< double, 2 > & a ) {

 if (a.empty()) {
  return;
 }
 boost::multi_array< double, 2 > temp = a;
 a.resize( boost::extents[ f_time_horizon ][ f_number_arcs ] );

 if( a.shape()[ 1 ] == 1 ) {
  for( Index t = 0; t < f_time_horizon; ++t ) {
   for( Index g = 0; g < f_number_arcs; ++g ) {
    a[ t ][ g ] = temp[ 0 ][ g ];
   }
  }

 } else if( a.shape()[ 0 ] < f_time_horizon ) {
  for( Index g = 0; g < f_number_arcs; ++g ) {
   int j = 0;
   for( unsigned long i = 0; i < v_change_intervals.size(); ++i ) {
    Index sup;
    if( i == v_change_intervals.size() - 1 ) {
     sup = f_time_horizon;
    } else {
     sup = v_change_intervals[ i ];
    }
    for( ; j < sup; ++j ) {
     a[ j ][ g ] = temp[ i ][ g ];
    }
   }
  }
 }
}

/*--------------------------------------------------------------------------*/

void HydroUnitBlock::decompress_vol( boost::multi_array< double, 2 > & a ) {

 if (a.empty()) {
  return;
 }
 boost::multi_array< double, 2 > temp = a;
 a.resize( boost::extents[ f_number_reservoirs ][ f_time_horizon ] );

 if( a.shape()[ 0 ] == 1 ) {

  for( Index n = 0; n < f_number_reservoirs; ++n ) {
   for( Index t = 0; t < f_time_horizon; ++t ) {
    a[ n ][ t ] = temp[ n ][ 0 ];
   }
  }

 } else if( a.shape()[ 1 ] < f_time_horizon ) {
  for( Index n = 0; n < f_number_reservoirs; ++n ) {
   int j = 0;
   for( unsigned long i = 0; i < v_change_intervals.size(); ++i ) {
    Index sup;
    if( i == v_change_intervals.size() - 1 ) {
     sup = f_time_horizon;
    } else {
     sup = v_change_intervals[ i ];
    }
    for( ; j < sup; ++j ) {
     a[ n ][ j ] = temp[ n ][ i ];
    }
   }
  }
 }
}

/*--------------------------------------------------------------------------*/
/*------------------- End File HydroUnitBlock.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
