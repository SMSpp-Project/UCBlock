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


 if (v_minimum_flow.empty()) {
  v_minimum_flow.resize( boost::extents[ 0 ][ 0 ] );
  ::deserialize( group, "MinFlow", v_minimum_flow, true, true );

 }

 if (v_minimum_volumetric.empty()) {
  v_minimum_volumetric.resize( boost::extents[ 1 ][ 1 ] );
 }
 if (v_maximum_volumetric.empty()) {
  v_maximum_volumetric.resize( boost::extents[ 1 ][ 1 ] );
 }

 if (v_minimum_power.empty()) {
  v_minimum_power.resize( boost::extents[ 0 ][ 0 ] );
 }
 if (v_maximum_power.empty()) {
  v_maximum_power.resize( boost::extents[ 0 ][ 0 ] );
 }
 if (v_delta_ramp_up.empty()) {
  v_delta_ramp_up.resize( boost::extents[ 0 ][ 0 ] );
 }
 if (v_delta_ramp_down.empty()) {
  v_delta_ramp_down.resize( boost::extents[ 0 ][ 0 ] );

  v_inertia_power.resize( boost::extents[ 0 ][ 0 ] );
 }

 ::deserialize_dim( group, "NumberReservoirs", f_number_reservoirs, true );

 ::deserialize_dim( group, "NumberArcs", f_number_arcs, true );

 ::deserialize( group, "NumberPieces", f_number_arcs ? : 1, v_number_pieces, true, true );

 ::deserialize( group, "StartLine", f_number_arcs ? : 1, v_start_arc);

 ::deserialize( group, "EndLine", f_number_arcs ? : 1, v_end_arc );

 if (v_inflows.empty()) {
  v_inflows.resize( boost::extents[ f_number_reservoirs ? : 1 ][ f_time_horizon ] );

 }
 ::deserialize( group, "Inflows", v_inflows, true, false );

 if (v_maximum_flow.empty()) {
  v_maximum_flow.resize( boost::extents[ 1 ][ 1 ] );
  ::deserialize( group, "MaxFlow", v_maximum_flow, true, true );

 }

 ::deserialize( group, "MinPower", v_minimum_power, true, true );

 ::deserialize( group, "MaxPower", v_maximum_power, true, true );

 ::deserialize( group, "DeltaRampUp", v_delta_ramp_up, true, true );

 ::deserialize( group, "DeltaRampDown", v_delta_ramp_down, true, true );


 ::deserialize( group, "PrimaryRho", v_primary_rho, true, true );


 ::deserialize( group, "SecondaryRho", v_secondary_rho, true, true );


 ::deserialize( group, "LinearTerm", f_number_arcs ? : 1, v_linear_term, true, true );


 ::deserialize( group, "ConstantTerm", f_total_number_pieces, v_const_term, false, true );

 ::deserialize( group, "InertiaPower", v_inertia_power, true, true );


 ::deserialize( group, "InitialFlowRate", f_number_arcs ? : 1, v_initial_flow_rate, true, true );


 ::deserialize( group, "InitialVolumetric", f_number_reservoirs ? : 1, v_initial_volumetric, true, true );


 ::deserialize( group, "UphillFlow", f_number_arcs ? : 1, v_uphill_delay, true, true );


 ::deserialize( group, "DownhillFlow", f_number_arcs ? : 1, v_downhill_delay, true, true );



 ::deserialize( group, "MinVolumetric", v_minimum_volumetric, true, true );

 ::deserialize( group, "MaxVolumetric", v_maximum_volumetric, true, true );

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

 if( !v_active_power.empty()  ) {
  // the abstract variables should be generated only once
  return;
 }
 v_active_power.resize(boost::extents[f_time_horizon][f_number_arcs]);

 for( Index t = 0; t < f_time_horizon; ++t ) {
  for( Index g = 0; g < f_number_arcs; ++g ) {
   v_active_power[ g ][ t ].set_type( ColVariable::kContinuous );

  }
 }
 add_static_variable ( v_active_power );

 if( !v_primary_spinning_reserve.empty()  ) {
  // the abstract variables should be generated only once
  return;
 }
 v_primary_spinning_reserve.resize(boost::extents[f_time_horizon][f_number_arcs]);

 for( Index t = 0; t < f_time_horizon; ++t ) {
  for( Index g = 0; g < f_number_arcs; ++g ) {
   v_primary_spinning_reserve[ g ][ t ].set_type( ColVariable::kNonNegative );

  }
 }
 add_static_variable ( v_primary_spinning_reserve );

 if( !v_secondary_spinning_reserve.empty()  ) {
  // the abstract variables should be generated only once
  return;
 }
 v_secondary_spinning_reserve.resize(boost::extents[f_time_horizon][f_number_arcs]);

 for( Index t = 0; t < f_time_horizon; ++t ) {
  for( Index g = 0; g < f_number_arcs; ++g ) {
   v_secondary_spinning_reserve[ g ][ t ].set_type( ColVariable::kNonNegative );

  }
 }
 add_static_variable ( v_secondary_spinning_reserve );
} // end( HydroUnitBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void HydroUnitBlock::generate_abstract_constraints( Configuration *stcc )
{
 // initial condition of matrix MinFlow
 boost::multi_array< double, 2 > MinFlow = v_minimum_flow;

 if ( MinFlow.size() == f_number_arcs ) {

  for ( Index t = 0; t < f_time_horizon; ++t ) {
   for ( Index g = 0; g < f_number_arcs; ++g ) {

   MinFlow.resize( boost::extents[f_time_horizon][f_number_arcs] );

   MinFlow[t][g] = v_minimum_flow[0][g];

  }
  }
 } else if ( MinFlow.size() < ( f_time_horizon * f_number_arcs ) ) {

  MinFlow.resize( boost::extents[f_time_horizon][f_number_arcs] );

  for( Index g = 0; g < f_number_arcs; ++g ) {

   int j = 0;
   for( unsigned long i = 0; i < v_change_intervals.size(); ++i ) {
    Index sup;
    if( i == v_change_intervals.size() - 1 ) {
     sup = f_time_horizon;
    } else {
     sup = v_change_intervals[i];
    }
    for( ; j < sup; ++j ) {

     MinFlow[j][g] = v_minimum_flow[i][g];
    }
   }
  }
 }


 // initial condition of matrix MaxFlow
 boost::multi_array< double, 2 > MaxFlow = v_maximum_flow;

 if ( MaxFlow.size() == f_number_arcs ) {

  for ( Index t = 0; t < f_time_horizon; ++t ) {
   for ( Index g = 0; g < f_number_arcs; ++g ) {

    MaxFlow.resize( boost::extents[f_time_horizon][f_number_arcs] );

    MaxFlow[t][g] = v_maximum_flow[0][g];

   }
  }
 } else if ( MaxFlow.size() < ( f_time_horizon * f_number_arcs ) ) {

  MaxFlow.resize( boost::extents[f_time_horizon][f_number_arcs] );

  for( Index g = 0; g < f_number_arcs; ++g ) {

   int j = 0;
   for( unsigned long i = 0; i < v_change_intervals.size(); ++i ) {
    Index sup;
    if( i == v_change_intervals.size() - 1 ) {
     sup = f_time_horizon;
    } else {
     sup = v_change_intervals[i];
    }
    for( ; j < sup; ++j ) {

     MaxFlow[j][g] = v_maximum_flow[i][g];
    }
   }
  }
 }

 // initial condition of matrix MinVolumetric
 boost::multi_array< double, 2 > MinVolumetric = v_minimum_volumetric;

 if ( MinVolumetric.size() == f_number_reservoirs ) {

  for ( Index t = 0; t < f_time_horizon; ++t ) {
   for ( Index n = 0; n < f_number_reservoirs; ++n ) {

    MinVolumetric.resize( boost::extents [f_number_reservoirs][f_time_horizon] );

    MinVolumetric[n][t] = v_minimum_volumetric[n][0];

   }
  }
 } else if ( MinFlow.size() < ( f_time_horizon * f_number_reservoirs ) ) {

  MinVolumetric.resize( boost::extents [f_number_reservoirs][f_time_horizon] );

  for( Index n = 0; n < f_number_reservoirs; ++n ) {

   int j = 0;
   for( unsigned long i = 0; i < v_change_intervals.size(); ++i ) {
    Index sup;
    if( i == v_change_intervals.size() - 1 ) {
     sup = f_time_horizon;
    } else {
     sup = v_change_intervals[i];
    }
    for( ; j < sup; ++j ) {

     MinVolumetric[n][j] = v_minimum_volumetric[n][i];
    }
   }
  }
 }

 // initial condition of matrix MaxVolumetric
 boost::multi_array< double, 2 > MaxVolumetric = v_maximum_volumetric;

 if ( MaxVolumetric.size() == f_number_reservoirs ) {

  for ( Index t = 0; t < f_time_horizon; ++t ) {
   for ( Index n = 0; n < f_number_reservoirs; ++n ) {

    MaxVolumetric.resize( boost::extents [f_number_reservoirs][f_time_horizon] );

    MaxVolumetric[n][t] = v_maximum_volumetric[n][0];

   }
  }
 } else if ( MaxVolumetric.size() < ( f_time_horizon * f_number_reservoirs ) ) {

  MaxVolumetric.resize( boost::extents [f_number_reservoirs][f_time_horizon] );

  for( Index n = 0; n < f_number_reservoirs; ++n ) {

   int j = 0;
   for( unsigned long i = 0; i < v_change_intervals.size(); ++i ) {
    Index sup;
    if( i == v_change_intervals.size() - 1 ) {
     sup = f_time_horizon;
    } else {
     sup = v_change_intervals[i];
    }
    for( ; j < sup; ++j ) {

     MaxVolumetric[n][j] = v_maximum_volumetric[n][i];
    }
   }
  }
 }

 // initial condition of matrix MaxVolumetric
 boost::multi_array< double, 2 > Inflows = v_inflows;

  for ( Index t = 0; t < f_time_horizon; ++t ) {
   for ( Index n = 0; n < f_number_reservoirs; ++n ) {

    Inflows.resize( boost::extents [f_number_reservoirs][f_time_horizon] );

    Inflows[n][t] = v_maximum_volumetric[n][t];

   }
  }

 // initial condition of matrix MinPower
 boost::multi_array< double, 2 > MinPower = v_minimum_power;

 if ( MinPower.size() == f_number_reservoirs ) {

  for ( Index t = 0; t < f_time_horizon; ++t ) {
   for ( Index n = 0; n < f_number_reservoirs; ++n ) {

    MinPower.resize( boost::extents [f_number_reservoirs][f_time_horizon] );

    MinPower[n][t] = v_minimum_power[n][0];

   }
  }
 } else if ( MinPower.size() < ( f_time_horizon * f_number_reservoirs ) ) {

  MinPower.resize( boost::extents [f_number_reservoirs][f_time_horizon] );

  for( Index n = 0; n < f_number_reservoirs; ++n ) {

   int j = 0;
   for( unsigned long i = 0; i < v_change_intervals.size(); ++i ) {
    Index sup;
    if( i == v_change_intervals.size() - 1 ) {
     sup = f_time_horizon;
    } else {
     sup = v_change_intervals[i];
    }
    for( ; j < sup; ++j ) {

     MinPower[n][j] = v_minimum_power[n][i];
    }
   }
  }
 }  // initial condition of matrix MaxPower
 boost::multi_array< double, 2 > MaxPower = v_maximum_power;

 if ( MaxPower.size() == f_number_reservoirs ) {

  for ( Index t = 0; t < f_time_horizon; ++t ) {
   for ( Index n = 0; n < f_number_reservoirs; ++n ) {

    MaxPower.resize( boost::extents [f_number_reservoirs][f_time_horizon] );

    MaxPower[n][t] = v_maximum_power[n][0];

   }
  }
 } else if ( MaxPower.size() < ( f_time_horizon * f_number_reservoirs ) ) {

  MaxPower.resize( boost::extents [f_number_reservoirs][f_time_horizon] );

  for( Index n = 0; n < f_number_reservoirs; ++n ) {

   int j = 0;
   for( unsigned long i = 0; i < v_change_intervals.size(); ++i ) {
    Index sup;
    if( i == v_change_intervals.size() - 1 ) {
     sup = f_time_horizon;
    } else {
     sup = v_change_intervals[i];
    }
    for( ; j < sup; ++j ) {

     MaxPower[n][j] = v_maximum_power[n][i];
    }
   }
  }
 }

 // initial condition of matrix DeltaRampUp
 boost::multi_array< double, 2 > DeltaRampUp = v_delta_ramp_up;

 if ( DeltaRampUp.size() == f_number_arcs ) {

  for ( Index t = 0; t < f_time_horizon; ++t ) {
   for ( Index g = 0; g < f_number_arcs; ++g ) {

    DeltaRampUp.resize( boost::extents[f_time_horizon][f_number_arcs] );

    DeltaRampUp[t][g] = v_delta_ramp_up[0][g];

   }
  }
 } else if ( DeltaRampUp.size() < ( f_time_horizon * f_number_arcs ) ) {

  DeltaRampUp.resize( boost::extents[f_time_horizon][f_number_arcs] );

  for( Index g = 0; g < f_number_arcs; ++g ) {

   int j = 0;
   for( unsigned long i = 0; i < v_change_intervals.size(); ++i ) {
    Index sup;
    if( i == v_change_intervals.size() - 1 ) {
     sup = f_time_horizon;
    } else {
     sup = v_change_intervals[i];
    }
    for( ; j < sup; ++j ) {

     DeltaRampUp[j][g] = v_delta_ramp_up[i][g];
    }
   }
  }
 }

 // initial condition of matrix DeltaRampDown
 boost::multi_array< double, 2 > DeltaRampDown = v_delta_ramp_down;

 if ( DeltaRampDown.size() == f_number_arcs ) {

  for ( Index t = 0; t < f_time_horizon; ++t ) {
   for ( Index g = 0; g < f_number_arcs; ++g ) {

    DeltaRampDown.resize( boost::extents[f_time_horizon][f_number_arcs] );

    DeltaRampDown[t][g] = v_delta_ramp_down[0][g];

   }
  }
 } else if ( DeltaRampDown.size() < ( f_time_horizon * f_number_arcs ) ) {

  DeltaRampDown.resize( boost::extents[f_time_horizon][f_number_arcs] );

  for( Index g = 0; g < f_number_arcs; ++g ) {

   int j = 0;
   for( unsigned long i = 0; i < v_change_intervals.size(); ++i ) {
    Index sup;
    if( i == v_change_intervals.size() - 1 ) {
     sup = f_time_horizon;
    } else {
     sup = v_change_intervals[i];
    }
    for( ; j < sup; ++j ) {

     DeltaRampDown[j][g] = v_delta_ramp_down[i][g];
    }
   }
  }
 }

 // initial condition of matrix PrimaryRho
 boost::multi_array< double, 2 > PrimaryRho = v_primary_rho;

 if ( PrimaryRho.size() == f_number_arcs ) {

  for ( Index t = 0; t < f_time_horizon; ++t ) {
   for ( Index g = 0; g < f_number_arcs; ++g ) {

    PrimaryRho.resize( boost::extents[f_time_horizon][f_number_arcs] );

    PrimaryRho[t][g] = v_primary_rho[0][g];

   }
  }
 } else if ( PrimaryRho.size() < ( f_time_horizon * f_number_arcs ) ) {

  PrimaryRho.resize( boost::extents[f_time_horizon][f_number_arcs] );

  for( Index g = 0; g < f_number_arcs; ++g ) {

   int j = 0;
   for( unsigned long i = 0; i < v_change_intervals.size(); ++i ) {
    Index sup;
    if( i == v_change_intervals.size() - 1 ) {
     sup = f_time_horizon;
    } else {
     sup = v_change_intervals[i];
    }
    for( ; j < sup; ++j ) {

     PrimaryRho[j][g] = v_primary_rho[i][g];
    }
   }
  }
 }

 // initial condition of matrix SecondaryRho
 boost::multi_array< double, 2 > SecondaryRho = v_secondary_rho;

 if ( SecondaryRho.size() == f_number_arcs ) {

  for ( Index t = 0; t < f_time_horizon; ++t ) {
   for ( Index g = 0; g < f_number_arcs; ++g ) {

    SecondaryRho.resize( boost::extents[f_time_horizon][f_number_arcs] );

    SecondaryRho[t][g] = v_secondary_rho[0][g];

   }
  }
 } else if ( SecondaryRho.size() < ( f_time_horizon * f_number_arcs ) ) {

  SecondaryRho.resize( boost::extents[f_time_horizon][f_number_arcs] );

  for( Index g = 0; g < f_number_arcs; ++g ) {

   int j = 0;
   for( unsigned long i = 0; i < v_change_intervals.size(); ++i ) {
    Index sup;
    if( i == v_change_intervals.size() - 1 ) {
     sup = f_time_horizon;
    } else {
     sup = v_change_intervals[i];
    }
    for( ; j < sup; ++j ) {

     SecondaryRho[j][g] = v_secondary_rho[i][g];
    }
   }
  }
 }

 // initial condition of vector NumberPieces
 std::vector< Index > NumberPieces = v_number_pieces;

 // initial condition of vector LinearTerm
 std::vector<double> LinearTerm = v_linear_term;

 // initial condition of vector ConstantTerm
 std::vector<double> ConstantTerm = v_const_term;

 // initial condition of matrix InertiaPower
 boost::multi_array< double, 2 > InertiaPower = v_inertia_power;

 if ( InertiaPower.size() == f_number_arcs ) {

  for ( Index t = 0; t < f_time_horizon; ++t ) {
   for ( Index g = 0; g < f_number_arcs; ++g ) {

    InertiaPower.resize( boost::extents[f_time_horizon][f_number_arcs] );

    InertiaPower[t][g] = v_inertia_power[0][g];

   }
  }
 } else if ( InertiaPower.size() < ( f_time_horizon * f_number_arcs ) ) {

  InertiaPower.resize( boost::extents[f_time_horizon][f_number_arcs] );

  for( Index g = 0; g < f_number_arcs; ++g ) {

   int j = 0;
   for( unsigned long i = 0; i < v_change_intervals.size(); ++i ) {
    Index sup;
    if( i == v_change_intervals.size() - 1 ) {
     sup = f_time_horizon;
    } else {
     sup = v_change_intervals[i];
    }
    for( ; j < sup; ++j ) {

     InertiaPower[j][g] = v_inertia_power[i][g];
    }
   }
  }
 }

 // initial condition of vector InitialFlowRate
 std::vector< double > InitialFlowRate = v_initial_flow_rate;

 // initial condition of vector InitialVolumetric
 std::vector< double > InitialVolumetric = v_initial_volumetric;

 // initial condition of vector UphillFlow
 std::vector< Index > UphillFlow = v_uphill_delay;

 // initial condition of vector DownhillFlow
 std::vector< Index > DownhillFlow = v_downhill_delay;

/*--------------------------------------------------------------------------*/

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

   linear_function->add_variable( &v_active_power[arc][t], 1.0 );
   linear_function->add_variable( &v_primary_spinning_reserve[arc][t], 1.0 );
   linear_function->add_variable( &v_secondary_spinning_reserve[arc][t], 1.0 );

   MaxPowerPrimarySecondary_Const[t][arc].set_lhs( -Inf< double >());
   MaxPowerPrimarySecondary_Const[t][arc].set_rhs( MaxPower[t][arc] );
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

   linear_function->add_variable( &v_active_power[arc][t], 1.0 );
   linear_function->add_variable( &v_primary_spinning_reserve[arc][t], -1.0 );
   linear_function->add_variable( &v_secondary_spinning_reserve[arc][t], -1.0 );

   MinPowerPrimarySecondary_Const[t][arc].set_lhs( MinPower[t][arc] );
   MinPowerPrimarySecondary_Const[t][arc].set_rhs( Inf< double >());
   MinPowerPrimarySecondary_Const[t][arc].set_function( linear_function );
  }
 }

 add_static_constraint( MinPowerPrimarySecondary_Const );

 // power output relation with to primary reserves constraints
 if( !v_primary_rho.empty() ) {
  if( ActivePowerPrimary_Const.size() != f_time_horizon ) {
   // this should only happen once
   assert( ActivePowerPrimary_Const.empty());

   ActivePowerPrimary_Const.resize
           ( boost::multi_array< FRowConstraint, 2 >::
             extent_gen()[f_time_horizon][f_number_arcs] );
  }

  for( Index t = 0; t < f_time_horizon; ++t ) {
   for( Index arc = 0; arc < f_number_arcs; ++arc ) {

    if( MinFlow[t][arc] >= 0 ) {

     auto linear_function = new LinearFunction();

     linear_function->add_variable( &v_active_power[arc][t], PrimaryRho[t][arc] );
     linear_function->add_variable( &v_primary_spinning_reserve[arc][t], -1.0 );
     ActivePowerPrimary_Const[t][arc].set_lhs( 0.0 );
     ActivePowerPrimary_Const[t][arc].set_rhs( Inf< double >());
     ActivePowerPrimary_Const[t][arc].set_function( linear_function );
    }

   }
  }

  add_static_constraint( ActivePowerPrimary_Const );
 }
 // power output relation with to secondary reserves constraints
 if( !v_secondary_rho.empty() ) {

  if( ActivePowerSecondary_Const.size() != f_time_horizon ) {
   // this should only happen once
   assert( ActivePowerSecondary_Const.empty());

   ActivePowerSecondary_Const.resize
           ( boost::multi_array< FRowConstraint, 2 >::
             extent_gen()[f_time_horizon][f_number_arcs] );
  }

  for( Index t = 0; t < f_time_horizon; ++t ) {
   for( Index arc = 0; arc < f_number_arcs; ++arc ) {

    if( MinFlow[t][arc] >= 0 ) {

     auto linear_function = new LinearFunction();

     linear_function->add_variable( &v_active_power[arc][t], SecondaryRho[t][arc] );
     linear_function->add_variable( &v_secondary_spinning_reserve[arc][t], -1.0 );
     ActivePowerSecondary_Const[t][arc].set_lhs( 0.0 );
     ActivePowerSecondary_Const[t][arc].set_rhs( Inf< double >());
     ActivePowerSecondary_Const[t][arc].set_function( linear_function );
    }

   }
  }

  add_static_constraint( ActivePowerSecondary_Const );
 }
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

   if( MaxFlow[t][arc] <= 0 ) {

    auto linear_function = new LinearFunction();

    linear_function->add_variable( &v_primary_spinning_reserve[arc][t], 1.0 );
    PrimaryPumps_Const[t][arc].set_both( 0.0 );
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

   if( MaxFlow[t][arc] <= 0 ) {

    auto linear_function = new LinearFunction();

    linear_function->add_variable( &v_secondary_spinning_reserve[arc][t], 1.0 );
    SecondaryPumps_Const[t][arc].set_both( 0.0 );
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

   //for pumps
   if( MaxFlow[t][arc] <= 0 ) {

    auto linear_function = new LinearFunction();

    linear_function->add_variable( &v_active_power[arc][t], 1.0 );
    linear_function->add_variable( &v_flow_rate[t][arc], -LinearTerm[arc] );
    FlowActivePowerPumps_Const[t][arc].set_both( 0.0 );
    FlowActivePowerPumps_Const[t][arc].set_function( linear_function );

   }
   add_static_constraint( FlowActivePowerPumps_Const );


   // flow to active power function constraints for turbines //todo

   if( FlowActivePowerTurbines_Const.size() != f_time_horizon ) {
    // this should only happen once
    assert( FlowActivePowerTurbines_Const.empty());

    FlowActivePowerTurbines_Const.resize
            ( boost::multi_array< FRowConstraint, 2 >::
              extent_gen()[f_time_horizon][f_number_arcs] );
   }

   if( MinFlow[t][arc] >= 0 ) {  //for turbines

    auto linear_function = new LinearFunction();

    linear_function->add_variable( &v_active_power[arc][t], 1.0 );
    linear_function->add_variable( &v_flow_rate[t][arc], -LinearTerm[arc] );
    FlowActivePowerTurbines_Const[t][arc].set_both( 0.0 );
    FlowActivePowerTurbines_Const[t][arc].set_function( linear_function );   }
  }
 }
 add_static_constraint( FlowActivePowerTurbines_Const );


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
    FlowRateBounds_Const[t][arc].set_lhs( MinFlow[t][arc] );
    FlowRateBounds_Const[t][arc].set_rhs( MaxFlow[t][arc] );
    FlowRateBounds_Const[t][arc].set_function( linear_function );

   }
  }

  add_static_constraint( FlowRateBounds_Const );



 // ram-up constraints
 if( !v_delta_ramp_up.empty() ) {

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

   RampUp_Const[0][arc].set_lhs( -Inf< double >());
   RampUp_Const[0][arc].set_rhs( DeltaRampUp[0][arc] + InitialFlowRate[arc] );
   RampUp_Const[0][arc].set_function( linear_function );

  }

  for( Index t = 1; t < f_time_horizon; ++t ) {
   for( Index arc = 0; arc < f_number_arcs; ++arc ) {
    auto linear_function = new LinearFunction();

    linear_function->add_variable( &v_flow_rate[t][arc], 1.0 );
    linear_function->add_variable( &v_flow_rate[t - 1][arc], -1.0 );

    RampUp_Const[t][arc].set_lhs( -Inf< double >());
    RampUp_Const[t][arc].set_rhs( DeltaRampUp[t][arc] );
    RampUp_Const[t][arc].set_function( linear_function );

   }
  }

  add_static_constraint( RampUp_Const );
 }


 // ram-down constraints
 if( !v_delta_ramp_down.empty()) {

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

   RampDown_Const[0][arc].set_lhs( InitialFlowRate[arc] - DeltaRampDown[0][arc] );
   RampDown_Const[0][arc].set_rhs( Inf< double >());
   RampDown_Const[0][arc].set_function( linear_function );

  }

  for( Index t = 1; t < f_time_horizon; ++t ) {
   for( Index arc = 0; arc < f_number_arcs; ++arc ) {
    auto linear_function = new LinearFunction();

    linear_function->add_variable( &v_flow_rate[t - 1][arc], 1.0 );
    linear_function->add_variable( &v_flow_rate[t][arc], -1.0 );

    RampDown_Const[t][arc].set_lhs( -Inf< double >());
    RampDown_Const[t][arc].set_rhs( DeltaRampDown[t][arc] );
    RampDown_Const[t][arc].set_function( linear_function );

   }
  }

  add_static_constraint( RampDown_Const );

 }

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
   VolumetricBounds_Const[node][t].set_lhs( MinVolumetric[node][t] );
   VolumetricBounds_Const[node][t].set_rhs( MaxVolumetric[node][t] );
   VolumetricBounds_Const[node][t].set_function( linear_function );

  }
 }

 add_static_constraint( VolumetricBounds_Const );

} // end( HydroUnitBlock::generate_abstract_constraints )


/*--------------------------------------------------------------------------*/
/*-------- METHODS FOR LOADING, PRINTING & SAVING THE HydroUnitBlock -------*/
/*--------------------------------------------------------------------------*/

void HydroUnitBlock::serialize( netCDF::NcGroup & group ) const {

 UnitBlock::serialize( group );

 auto dim_time_horizon = group.addDim( "TimeHorizon", f_time_horizon );
 auto NumberIntervals = group.getDim( "NumberIntervals" );


  auto dim_total_number_pieces = group.addDim( "TotalNumberPieces",
                                               f_total_number_pieces );

 auto dim_number_reservoirs =
         group.addDim( "NumberReservoirs", f_number_reservoirs );

 auto dim_number_arcs =
         group.addDim( "NumberArcs", f_number_arcs );

 if( f_number_reservoirs > 1 ) {

  ::serialize( group, "StartLine", netCDF::NcInt64(),
               dim_number_reservoirs, v_start_arc, false );

  ::serialize( group, "EndLine", netCDF::NcInt64(),
               dim_number_reservoirs, v_end_arc, false );
 }

 if( !v_minimum_flow.empty() ) {

  ::serialize( group, "MinFlow", netCDF::NcDouble(),
               {NumberIntervals, dim_number_arcs},
               v_minimum_flow, true );
 }

 if( !v_maximum_flow.empty() ) {

  ::serialize( group, "MaxFlow", netCDF::NcDouble(),
               {NumberIntervals, dim_number_arcs},
               v_maximum_flow, true );
 }

 if( !v_minimum_volumetric.empty() ) {

  ::serialize( group, "MinVolumetric", netCDF::NcDouble(),
               {dim_number_reservoirs, NumberIntervals},
               v_minimum_volumetric, true );
 }

 if( !v_maximum_volumetric.empty() ) {

  ::serialize( group, "MaxVolumetric", netCDF::NcDouble(),
               {dim_number_reservoirs, NumberIntervals},
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
               {NumberIntervals, dim_number_arcs},
               v_primary_rho, true );
 }

 if( !v_secondary_rho.empty() ) {

  ::serialize( group, "SecondaryRho", netCDF::NcDouble(),
               {NumberIntervals, dim_number_arcs},
               v_secondary_rho, true );
 }


 ::serialize( group, "NumberPieces", netCDF::NcUint64(),
              dim_number_arcs, v_number_pieces, false );

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

 ::serialize( group, "UphillFlow", netCDF::NcInt64(),
              dim_number_arcs, v_uphill_delay, true );

 ::serialize( group, "DownhillFlow", netCDF::NcUint64(),
              dim_number_arcs, v_downhill_delay, true );
}  // end( HydroUnitBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------------- End File HydroUnitBlock.cpp --------------------------*/
/*--------------------------------------------------------------------------*/