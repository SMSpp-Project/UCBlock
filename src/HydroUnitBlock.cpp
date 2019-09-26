/*--------------------------------------------------------------------------*/
/*--------------------- File HydroUnitBlock.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the HydroUnitBlock class.
 *
 * \version 0.11
 *
 * \date 11 - 07 - 2019
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 *
 * \author Ali Ghezelsoflu \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 *
 * Copyright &copy by Antonio Frangioni, Ali Ghezelsoflu
 */

/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <iostream>
#include <random>
#include "HydroUnitBlock.h"
#include "LinearFunction.h"
#include <map>
#include "FRowConstraint.h"
#include "UnitBlock.h"
#include "UCBlock.h"

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
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
void HydroUnitBlock::deserialize( netCDF::NcGroup & group ) {

 UnitBlock::deserialize( group );

 ::deserialize_dim( group, "NumberReservoirs", f_number_reservoirs, true );

 if( f_number_reservoirs > 1 ) {

  ::deserialize_dim( group, "NumberArcs", f_number_arcs, false );

  ::deserialize( group, "StartLine", f_number_reservoirs, v_start_arc);

  ::deserialize( group, "EndLine", f_number_reservoirs, v_end_arc );  //todo
 }
}// end( HydroUnitBlock::deserialize )

/*--------------------------------------------------------------------------*/

void HydroUnitBlock::generate_abstract_variables( Configuration *stvv )
{
 UnitBlock::generate_abstract_variables( stvv );

 if( f_time_horizon == 0 ) {
  // there are no variables to be generated
  return;
 }

 if( !v_volumetric.empty() ||
     !v_flow_rate.empty() ) {
  // the abstract variables should be generated only once
  return;
 }

  v_volumetric.resize(boost::extents[f_time_horizon][f_number_arcs]);
  v_flow_rate.resize(boost::extents[f_time_horizon][f_number_arcs]);

 for( Index t = 0; t < f_time_horizon; ++t ) {
  for( Index g = 0; g < f_number_arcs; ++g ) {
   v_volumetric[ t ][ g ].set_type( ColVariable::kNonNegative );
   v_flow_rate[ t ][ g ].set_type( ColVariable::kContinuous );

  }
 }
 add_static_variable ( v_volumetric );
 add_static_variable ( v_flow_rate );


} // end( HydroUnitBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void HydroUnitBlock::generate_abstract_constraints( Configuration *stcc )
{

 // maximum power output according to primary-secondary reserves constraints

 if( MaxPowerPrimarySecondary_Const.size() != f_time_horizon ) {
  // this should only happen once
  assert( MaxPowerPrimarySecondary_Const.empty());

  MaxPowerPrimarySecondary_Const.resize
          ( boost::multi_array< FRowConstraint, 2 >::
            extent_gen()[f_time_horizon][f_number_arcs] );
 }

 for( Index t = 0; t < f_time_horizon; ++t ) {
  for( Index arc = 0; arc < f_number_arcs; ++arc ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_active_power[t][arc], 1.0 );
   linear_function->add_variable( &v_primary_spinning_reserve[t][arc], 1.0 );
   linear_function->add_variable( &v_secondary_spinning_reserve[t][arc], 1.0 );

   MaxPowerPrimarySecondary_Const[t][arc].set_lhs( -Inf< double >());
   MaxPowerPrimarySecondary_Const[t][arc].set_rhs( v_maximum_power[t][arc] );
   MaxPowerPrimarySecondary_Const[t][arc].set_function( linear_function );
  }
 }

 add_static_constraint( MaxPowerPrimarySecondary_Const );


 // minimum power output according to primary-secondary reserves constraints

 if( MinPowerPrimarySecondary_Const.size() != f_time_horizon ) {
  // this should only happen once
  assert( MinPowerPrimarySecondary_Const.empty());

  MinPowerPrimarySecondary_Const.resize
          ( boost::multi_array< FRowConstraint, 2 >::
            extent_gen()[f_time_horizon][f_number_arcs] );
 }

 for( Index t = 0; t < f_time_horizon; ++t ) {
  for( Index arc = 0; arc < f_number_arcs; ++arc ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_active_power[t][arc], 1.0 );
   linear_function->add_variable( &v_primary_spinning_reserve[t][arc], -1.0 );
   linear_function->add_variable( &v_secondary_spinning_reserve[t][arc], -1.0 );

   MinPowerPrimarySecondary_Const[t][arc].set_lhs( v_minimum_power[t][arc] );
   MinPowerPrimarySecondary_Const[t][arc].set_rhs( Inf< double >());
   MinPowerPrimarySecondary_Const[t][arc].set_function( linear_function );
  }
 }

 add_static_constraint( MinPowerPrimarySecondary_Const );

 // power output relation with to primary reserves constraints

 if( ActivePowerPrimary_Const.size() != f_time_horizon ) {
  // this should only happen once
  assert( ActivePowerPrimary_Const.empty());

  ActivePowerPrimary_Const.resize
          ( boost::multi_array< FRowConstraint, 2 >::
            extent_gen()[f_time_horizon][f_number_arcs] );
 }

 for( Index t = 0; t < f_time_horizon; ++t ) {
  for( Index arc = 0; arc < f_number_arcs; ++arc ) {

   if( v_minimum_volumetric[t][arc] >= 0 && v_maximum_volumetric[t][arc] >= 0 ) {

    auto linear_function = new LinearFunction();

    linear_function->add_variable( &v_active_power[t][arc], v_primary_rho[t][arc] );
    linear_function->add_variable( &v_primary_spinning_reserve[t][arc], -1.0 );
    ActivePowerPrimary_Const[t][arc].set_lhs( 0.0 );
    ActivePowerPrimary_Const[t][arc].set_rhs( Inf< double >());
    ActivePowerPrimary_Const[t][arc].set_function( linear_function );
   }

  }
 }

 add_static_constraint( ActivePowerPrimary_Const );

 // power output relation with to secondary reserves constraints

 if( ActivePowerSecondary_Const.size() != f_time_horizon ) {
  // this should only happen once
  assert( ActivePowerSecondary_Const.empty());

  ActivePowerSecondary_Const.resize
          ( boost::multi_array< FRowConstraint, 2 >::
            extent_gen()[f_time_horizon][f_number_arcs] );
 }

 for( Index t = 0; t < f_time_horizon; ++t ) {
  for( Index arc = 0; arc < f_number_arcs; ++arc ) {

   if( v_minimum_volumetric[t][arc] >= 0 && v_maximum_volumetric[t][arc] >= 0 ) {

    auto linear_function = new LinearFunction();

    linear_function->add_variable( &v_active_power[t][arc], v_secondary_rho[t][arc] );
    linear_function->add_variable( &v_secondary_spinning_reserve[t][arc], -1.0 );
    ActivePowerSecondary_Const[t][arc].set_lhs( 0.0 );
    ActivePowerSecondary_Const[t][arc].set_rhs( Inf< double >());
    ActivePowerSecondary_Const[t][arc].set_function( linear_function );
   }

  }
 }

 add_static_constraint( ActivePowerSecondary_Const );

 // primary reserves constraints for pumps

 if( PrimaryPumps_Const.size() != f_time_horizon ) {
  // this should only happen once
  assert( PrimaryPumps_Const.empty());

  PrimaryPumps_Const.resize
          ( boost::multi_array< FRowConstraint, 2 >::
            extent_gen()[f_time_horizon][f_number_arcs] );
 }

 for( Index t = 0; t < f_time_horizon; ++t ) {
  for( Index arc = 0; arc < f_number_arcs; ++arc ) {

   if( v_minimum_volumetric[t][arc] < 0 && v_maximum_volumetric[t][arc] < 0 ) {

    auto linear_function = new LinearFunction();

    linear_function->add_variable( &v_primary_spinning_reserve[t][arc], 1.0 );
    PrimaryPumps_Const[t][arc].set_lhs( 0.0 );
    PrimaryPumps_Const[t][arc].set_rhs( 0.0 );
    PrimaryPumps_Const[t][arc].set_function( linear_function );
   }

  }
 }

 add_static_constraint( PrimaryPumps_Const );

 // secondary reserves constraints for pumps

 if( SecondaryPumps_Const.size() != f_time_horizon ) {
  // this should only happen once
  assert( SecondaryPumps_Const.empty());

  SecondaryPumps_Const.resize
          ( boost::multi_array< FRowConstraint, 2 >::
            extent_gen()[f_time_horizon][f_number_arcs] );
 }

 for( Index t = 0; t < f_time_horizon; ++t ) {
  for( Index arc = 0; arc < f_number_arcs; ++arc ) {

   if( v_minimum_volumetric[t][arc] < 0 && v_maximum_volumetric[t][arc] < 0 ) {

    auto linear_function = new LinearFunction();

    linear_function->add_variable( &v_secondary_spinning_reserve[t][arc], 1.0 );
    SecondaryPumps_Const[t][arc].set_lhs( 0.0 );
    SecondaryPumps_Const[t][arc].set_rhs( 0.0 );
    SecondaryPumps_Const[t][arc].set_function( linear_function );
   }

  }
 }

 add_static_constraint( SecondaryPumps_Const );

 // flow to active power function constraints for pumps

 if( FlowActivePowerPumps_Const.size() != f_time_horizon ) {
  // this should only happen once
  assert( FlowActivePowerPumps_Const.empty());

  FlowActivePowerPumps_Const.resize
          ( boost::multi_array< FRowConstraint, 2 >::
            extent_gen()[f_time_horizon][f_number_arcs] );
 }

 for( Index t = 0; t < f_time_horizon; ++t ) {
  for( Index arc = 0; arc < f_number_arcs; ++arc ) {

   if( v_minimum_volumetric[t][arc] < 0 && v_maximum_volumetric[t][arc] < 0 ) {

    auto linear_function = new LinearFunction();

    linear_function->add_variable( &v_active_power[t][arc], 1.0 );
    linear_function->add_variable( &v_flow_rate[t][arc], -v_linear_term[arc] );
    FlowActivePowerPumps_Const[t][arc].set_lhs( 0.0 );
    FlowActivePowerPumps_Const[t][arc].set_rhs( 0.0 );
    FlowActivePowerPumps_Const[t][arc].set_function( linear_function );
   }

  }
 }

 add_static_constraint( FlowActivePowerPumps_Const );



 // flow to active power function constraints for turbines

 if( FlowActivePowerTurbines_Const.size() != f_time_horizon ) {
  // this should only happen once
  assert( FlowActivePowerTurbines_Const.empty());

  //TODO NEED TO WRITE







 }

  // flow rate bounds constraints

  if( FlowRateBounds_Const.size() != f_time_horizon ) {
   // this should only happen once
   assert( FlowRateBounds_Const.empty());

   FlowRateBounds_Const.resize
           ( boost::multi_array< FRowConstraint, 2 >::
             extent_gen()[f_time_horizon][f_number_arcs] );
  }

  for( Index t = 0; t < f_time_horizon; ++t ) {
   for( Index arc = 0; arc < f_number_arcs; ++arc ) {

    auto linear_function = new LinearFunction();

    linear_function->add_variable( &v_flow_rate[t][arc], 1.0 );
    FlowRateBounds_Const[t][arc].set_lhs( v_minimum_flow[t][arc] );
    FlowRateBounds_Const[t][arc].set_rhs( v_maximum_flow[t][arc] );
    FlowRateBounds_Const[t][arc].set_function( linear_function );

   }
  }

  add_static_constraint( FlowRateBounds_Const );

 // ram-up constraints


 if( RampUp_Const.size() != f_time_horizon ) {
  // this should only happen once
  assert( RampUp_Const.empty());

  RampUp_Const.resize
          ( boost::multi_array< FRowConstraint, 2 >::
            extent_gen()[f_time_horizon][f_number_arcs] );
 }
 // Initial condition

  for( Index arc = 0; arc < f_number_arcs; ++arc ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_flow_rate[0][arc], 1.0 );

   RampUp_Const[0][arc].set_lhs( -Inf< double >() );
   RampUp_Const[0][arc].set_rhs( v_delta_ramp_up[0][arc] + v_initial_flow_rate[arc] );
   RampUp_Const[0][arc].set_function( linear_function );

  }

   for( Index t = 1; t < f_time_horizon; ++t ) {
    for( Index arc = 0; arc < f_number_arcs; ++arc ) {
    auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_flow_rate[t][arc], 1.0 );
   linear_function->add_variable( &v_flow_rate[t-1][arc], -1.0 );

   RampUp_Const[t][arc].set_lhs( -Inf< double >());
   RampUp_Const[t][arc].set_rhs( v_delta_ramp_up[t][arc] );
   RampUp_Const[t][arc].set_function( linear_function );

  }
 }

 add_static_constraint( RampUp_Const );



 // ram-down constraints

 if( RampDown_Const.size() != f_time_horizon ) {
  // this should only happen once
  assert( RampDown_Const.empty());

  RampDown_Const.resize
          ( boost::multi_array< FRowConstraint, 2 >::
            extent_gen()[f_time_horizon][f_number_arcs] );
 }

 // Initial condition

 for( Index arc = 0; arc < f_number_arcs; ++arc ) {

  auto linear_function = new LinearFunction();

  linear_function->add_variable( &v_flow_rate[0][arc], 1.0 );

  RampDown_Const[0][arc].set_lhs( v_initial_flow_rate[arc] - v_delta_ramp_down[0][arc]);
  RampDown_Const[0][arc].set_rhs( Inf< double >() );
  RampDown_Const[0][arc].set_function( linear_function );

 }

 for( Index t = 1; t < f_time_horizon; ++t ) {
  for( Index arc = 0; arc < f_number_arcs; ++arc ) {
   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_flow_rate[t-1][arc], 1.0 );
   linear_function->add_variable( &v_flow_rate[t][arc], -1.0 );

   RampDown_Const[t][arc].set_lhs( -Inf< double >());
   RampDown_Const[t][arc].set_rhs( v_delta_ramp_down[t][arc] );
   RampDown_Const[t][arc].set_function( linear_function );

  }
 }

 add_static_constraint( RampDown_Const );



 // final volumes fo each reservoir constraints

 if( FinalVolumeReservoir_Const.size() != f_time_horizon ) {
  // this should only happen once
  assert( FinalVolumeReservoir_Const.empty());

  //TODO NEED TO WRITE







 }

 // volumetric bounds constraints

 if( VolumetricBounds_Const.size() != f_time_horizon ) {
  // this should only happen once
  assert( VolumetricBounds_Const.empty());

  VolumetricBounds_Const.resize
          ( boost::multi_array< FRowConstraint, 2 >::
            extent_gen()[f_number_reservoirs][f_time_horizon] );
 }

 for( Index t = 0; t < f_time_horizon; ++t ) {
  for( Index node = 0; node < f_number_reservoirs; ++node ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_volumetric[node][t], 1.0 );
   VolumetricBounds_Const[node][t].set_lhs( v_minimum_volumetric[node][t] );
   VolumetricBounds_Const[node][t].set_rhs( v_maximum_volumetric[node][t] );
   VolumetricBounds_Const[node][t].set_function( linear_function );

  }
 }

 add_static_constraint( VolumetricBounds_Const );

} // end( HydroUnitBlock::generate_abstract_constraints )


/*--------------------------------------------------------------------------*/
/*-------- METHODS FOR LOADING, PRINTING & SAVING THE HydroUnitBlock -------*/
/*--------------------------------------------------------------------------*/

void HydroUnitBlock::serialize( netCDF::NcGroup & group ) const {

 //TODO I should complete this part after completing the interface

}  // end( HydroUnitBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------------- End File HydroUnitBlock.cpp --------------------------*/
/*--------------------------------------------------------------------------*/