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
  ::deserialize( group, "MinFlow", v_minimum_flow, true, true );
 }

 long rows = v_minimum_flow.shape()[0];
 long cols = v_minimum_flow.shape()[1];
 if (rows > 1 && cols == 1) {
  // The vector must be transposed
  boost::array<boost::multi_array<double, 2>::index, 2> dims = {{1, rows}};
  v_minimum_flow.reshape(dims);
 }

 if (v_maximum_flow.empty()) {
  ::deserialize( group, "MaxFlow", v_maximum_flow, true, true );
 }

 rows = v_maximum_flow.shape()[0];
 cols = v_maximum_flow.shape()[1];
 if (rows > 1 && cols == 1) {
  // The vector must be transposed
  boost::array<boost::multi_array<double, 2>::index, 2> dims = {{1, rows}};
  v_maximum_flow.reshape(dims);
 }

 ::deserialize_dim( group, "NumberReservoirs", f_number_reservoirs, true );

 ::deserialize_dim( group, "NumberArcs", f_number_arcs, true );

 ::deserialize( group, "NumberPieces", f_number_arcs ? f_number_arcs : 1, v_number_pieces, true, true );

 for (auto& n : v_number_pieces) {
  f_total_number_pieces += n;
 }

 ::deserialize( group, "StartArc", f_number_arcs ? f_number_arcs : 1, v_start_arc);

 ::deserialize( group, "EndArc", f_number_arcs ? f_number_arcs : 1, v_end_arc );

 ::deserialize( group, "Inflows", v_inflows, true, false );

 rows = v_inflows.shape()[0];
 cols = v_inflows.shape()[1];
 if (rows > 1 && cols == 1) {
  // The vector must be transposed
  boost::array<boost::multi_array<double, 2>::index, 2> dims = {{1, rows}};
  v_inflows.reshape(dims);
 }

 ::deserialize( group, "MinPower", v_minimum_power, true, true );

 rows = v_minimum_power.shape()[0];
 cols = v_minimum_power.shape()[1];
 if (rows > 1 && cols == 1) {
  // The vector must be transposed
  boost::array<boost::multi_array<double, 2>::index, 2> dims = {{1, rows}};
  v_minimum_power.reshape(dims);
 }

 ::deserialize( group, "MaxPower", v_maximum_power, true, true );

 rows = v_maximum_power.shape()[0];
 cols = v_maximum_power.shape()[1];
 if (rows > 1 && cols == 1) {
  // The vector must be transposed
  boost::array<boost::multi_array<double, 2>::index, 2> dims = {{1, rows}};
  v_maximum_power.reshape(dims);
 }

 ::deserialize( group, "DeltaRampUp", v_delta_ramp_up, true, true );

 rows = v_delta_ramp_up.shape()[0];
 cols = v_delta_ramp_up.shape()[1];
 if (rows > 1 && cols == 1) {
  // The vector must be transposed
  boost::array<boost::multi_array<double, 2>::index, 2> dims = {{1, rows}};
  v_delta_ramp_up.reshape(dims);
 }

 ::deserialize( group, "DeltaRampDown", v_delta_ramp_down, true, true );

 rows = v_delta_ramp_down.shape()[0];
 cols = v_delta_ramp_down.shape()[1];
 if (rows > 1 && cols == 1) {
  // The vector must be transposed
  boost::array<boost::multi_array<double, 2>::index, 2> dims = {{1, rows}};
  v_delta_ramp_down.reshape(dims);
 }

 ::deserialize( group, "PrimaryRho", v_primary_rho, true, true );

 rows = v_primary_rho.shape()[0];
 cols = v_primary_rho.shape()[1];
 if (rows > 1 && cols == 1) {
  // The vector must be transposed
  boost::array<boost::multi_array<double, 2>::index, 2> dims = {{1, rows}};
  v_primary_rho.reshape(dims);
 }

 ::deserialize( group, "SecondaryRho", v_secondary_rho, true, true );

 rows = v_secondary_rho.shape()[0];
 cols = v_secondary_rho.shape()[1];
 if (rows > 1 && cols == 1) {
  // The vector must be transposed
  boost::array<boost::multi_array<double, 2>::index, 2> dims = {{1, rows}};
  v_secondary_rho.reshape(dims);
 }

 ::deserialize( group, "LinearTerm", f_total_number_pieces ? f_total_number_pieces : 1, v_linear_term, true, true );

 ::deserialize( group, "ConstantTerm", f_total_number_pieces ? f_total_number_pieces : 1, v_const_term, true, true );

 ::deserialize( group, "InertiaPower", v_inertia_power, true, true );

 ::deserialize( group, "InitialFlowRate", f_number_arcs ? f_number_arcs : 1, v_initial_flow_rate, true, true );

 ::deserialize( group, "InitialVolumetric", f_number_reservoirs ? f_number_reservoirs : 1, v_initial_volumetric, true, true );

 ::deserialize( group, "UphillFlow", f_number_arcs ? f_number_arcs : 1, v_uphill_delay, true, true );

 ::deserialize( group, "DownhillFlow", f_number_arcs ? f_number_arcs : 1, v_downhill_delay, true, true );

 ::deserialize( group, "MinVolumetric", v_minimum_volumetric, true, true );

 rows = v_minimum_volumetric.shape()[0];
 cols = v_minimum_volumetric.shape()[1];
 if (rows > 1 && cols == 1) {
  // The vector must be transposed
  boost::array<boost::multi_array<double, 2>::index, 2> dims = {{1, rows}};
  v_minimum_volumetric.reshape(dims);
 }
 ::deserialize( group, "MaxVolumetric", v_maximum_volumetric, true, true );

 rows = v_maximum_volumetric.shape()[0];
 cols = v_maximum_volumetric.shape()[1];
 if (rows > 1 && cols == 1) {
  // The vector must be transposed
  boost::array<boost::multi_array<double, 2>::index, 2> dims = {{1, rows}};
  v_maximum_volumetric.reshape(dims);
 }
}// end( HydroUnitBlock::deserialize )

/*--------------------------------------------------------------------------*/

void HydroUnitBlock::generate_abstract_variables( Configuration *stvv )
{
 UnitBlock::generate_abstract_variables( stvv );

 unsigned int number_arcs = f_number_arcs ? f_number_arcs : 1;
 unsigned int number_reservoirs = f_number_reservoirs ? f_number_reservoirs : 1;

 if( f_time_horizon == 0 ) {
  // there are no variables to be generated
  return;
 }

 if( !v_volumetric.empty()  ) {
  // the abstract variables should be generated only once
  return;
 }
 v_volumetric.resize(boost::extents[number_reservoirs][ f_time_horizon]);
 for( Index g = 0; g < number_reservoirs; ++g ) {
 for( Index t = 0; t < f_time_horizon; ++t ) {
   auto & volumetric = v_volumetric[ g ][ t ];

   volumetric.set_type( ColVariable::kNonNegative );
  }
 }
 add_static_variable ( v_volumetric, "vol" );

 if( !v_flow_rate.empty() ) {
  // the abstract variables should be generated only once
  return;
 }
 v_flow_rate.resize(boost::extents[number_arcs][f_time_horizon]);
 for( Index g = 0; g < number_arcs; ++g ) {
 for( Index t = 0; t < f_time_horizon; ++t ) {
   auto & flow_rate = v_flow_rate[ g ][ t ];

   flow_rate.set_type( ColVariable::kContinuous );
  }
 }
 add_static_variable ( v_flow_rate, "F" );

 if( !v_active_power.empty()  ) {
  // the abstract variables should be generated only once
  return;
 }

 int n_gen = f_number_arcs == 0 ? 1 : f_number_arcs;
 v_active_power.resize(boost::extents[n_gen][f_time_horizon]);

 for( Index g = 0; g < n_gen; ++g ) {
  for( Index t = 0; t < f_time_horizon; ++t ) {
   v_active_power[ g ][ t ].set_type( ColVariable::kContinuous );

  }
 }
 add_static_variable ( v_active_power, "p" );

  if( !v_primary_spinning_reserve.empty()  ) {
   // the abstract variables should be generated only once
   return;
  }
  v_primary_spinning_reserve.resize(boost::extents[n_gen][f_time_horizon]);

 for( Index g = 0; g < n_gen; ++g ) {
  for( Index t = 0; t < f_time_horizon; ++t ) {
    v_primary_spinning_reserve[ g ][ t ].set_type( ColVariable::kNonNegative );

   }
  }
  add_static_variable ( v_primary_spinning_reserve, "pr" );

  if( !v_secondary_spinning_reserve.empty()  ) {
   // the abstract variables should be generated only once
  return;
  }
 v_secondary_spinning_reserve.resize(boost::extents[n_gen][f_time_horizon]);

 for( Index g = 0; g < n_gen; ++g ) {
  for( Index t = 0; t < f_time_horizon; ++t ) {
 v_secondary_spinning_reserve[ g ][ t ].set_type( ColVariable::kNonNegative );

  }
 }
  add_static_variable ( v_secondary_spinning_reserve, "sr" );
} // end( HydroUnitBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void HydroUnitBlock::generate_abstract_constraints( Configuration *stcc ) {

 unsigned int number_arcs = f_number_arcs ? f_number_arcs : 1;
 unsigned int number_reservoirs = f_number_reservoirs ? f_number_reservoirs : 1;

 // initial condition of matrix MinFlow
 boost::multi_array< double, 2 > MinFlow = v_minimum_flow;

 if( MinFlow.shape()[ 0 ] == 1 ) {
  MinFlow.resize( boost::extents[ f_time_horizon ][ number_arcs ] );
  for( Index t = 0; t < f_time_horizon; ++t ) {
   for( Index g = 0; g < number_arcs; ++g ) {
    MinFlow[ t ][ g ] = v_minimum_flow[ 0 ][ g ];
   }
  }

 } else if( MinFlow.shape()[ 0 ] < f_time_horizon ) {

  MinFlow.resize( boost::extents[ f_time_horizon ][ number_arcs ] );

  for( Index g = 0; g < number_arcs; ++g ) {

   int j = 0;
   for( unsigned long i = 0; i < v_change_intervals.size(); ++i ) {
    Index sup;
    if( i == v_change_intervals.size() - 1 ) {
     sup = f_time_horizon;
    } else {
     sup = v_change_intervals[ i ];
    }
    for( ; j < sup; ++j ) {

     MinFlow[ j ][ g ] = v_minimum_flow[ i ][ g ];
    }
   }
  }
 }


 // initial condition of matrix MaxFlow
 boost::multi_array< double, 2 > MaxFlow = v_maximum_flow;

 if( MaxFlow.shape()[ 0 ] == 1 ) {
  MaxFlow.resize( boost::extents[ f_time_horizon ][ number_arcs ] );
  for( Index t = 0; t < f_time_horizon; ++t ) {
   for( Index g = 0; g < number_arcs; ++g ) {
    MaxFlow[ t ][ g ] = v_maximum_flow[ 0 ][ g ];

   }
  }
 } else if( MaxFlow.shape()[ 0 ] < f_time_horizon ) {

  MaxFlow.resize( boost::extents[ f_time_horizon ][ number_arcs ] );

  for( Index g = 0; g < number_arcs; ++g ) {

   int j = 0;
   for( unsigned long i = 0; i < v_change_intervals.size(); ++i ) {
    Index sup;
    if( i == v_change_intervals.size() - 1 ) {
     sup = f_time_horizon;
    } else {
     sup = v_change_intervals[ i ];
    }
    for( ; j < sup; ++j ) {

     MaxFlow[ j ][ g ] = v_maximum_flow[ i ][ g ];
    }
   }
  }
 }

 // initial condition of matrix MinVolumetric
 boost::multi_array< double, 2 > MinVolumetric = v_minimum_volumetric;

 if( MinVolumetric.shape()[ 1 ] == 1 ) {
  MinVolumetric.resize( boost::extents[ number_reservoirs ][ f_time_horizon ] );

  for( Index n = 0; n < number_reservoirs; ++n ) {
   for( Index t = 0; t < f_time_horizon; ++t ) {

    MinVolumetric[ n ][ t ] = v_minimum_volumetric[ n ][ 0 ];

   }
  }
 } else if( MinVolumetric.shape()[ 1 ] < f_time_horizon ) {

  MinVolumetric.resize( boost::extents[ number_reservoirs ][ f_time_horizon ] );

  for( Index n = 0; n < number_reservoirs; ++n ) {

   int j = 0;
   for( unsigned long i = 0; i < v_change_intervals.size(); ++i ) {
    Index sup;
    if( i == v_change_intervals.size() - 1 ) {
     sup = f_time_horizon;
    } else {
     sup = v_change_intervals[ i ];
    }
    for( ; j < sup; ++j ) {

     MinVolumetric[ n ][ j ] = v_minimum_volumetric[ n ][ i ];
    }
   }
  }
 }

 // initial condition of matrix MaxVolumetric
 boost::multi_array< double, 2 > MaxVolumetric = v_maximum_volumetric;

 if( MaxVolumetric.shape()[ 1 ] == 1 ) {
  MaxVolumetric.resize( boost::extents[number_reservoirs][f_time_horizon] );

   for( Index n = 0; n < number_reservoirs; ++n ) {
    for( Index t = 0; t < f_time_horizon; ++t ) {


    MaxVolumetric[n][t] = v_maximum_volumetric[n][0];

   }
  }
 } else if( MaxVolumetric.shape()[ 1 ] < f_time_horizon ) {

  MaxVolumetric.resize( boost::extents[number_reservoirs][f_time_horizon] );

  for( Index n = 0; n < number_reservoirs; ++n ) {

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

 // initial condition of matrix Inflows
 boost::multi_array< double, 2 > Inflows;

 if( number_reservoirs == 1 ) {
  Inflows.resize( boost::extents[1][f_time_horizon] );
  for( Index t = 0; t < f_time_horizon; ++t ) {
   Inflows[0][t] = v_inflows[0][t];

  }
 } else {
  Inflows.resize( boost::extents[number_reservoirs][f_time_horizon] );
  for( Index g = 0; g < number_reservoirs; ++g ) {
   for( Index t = 0; t < f_time_horizon; ++t ) {
    Inflows[g][t] = v_inflows[g][t];
   }
  }
 }



 // initial condition of matrix MinPower
 boost::multi_array< double, 2 > MinPower = v_minimum_power;

 if( MinPower.shape()[ 0 ] == 1 ) {
  MinPower.resize( boost::extents[ f_time_horizon ][ number_arcs ] );

  for( Index n = 0; n < number_arcs; ++n ) {
   for( Index t = 0; t < f_time_horizon; ++t ) {
    MinPower[ t ][ n ] = v_minimum_power[ 0 ][ n ];
   }
  }

 } else if( MinPower.shape()[ 0 ] < f_time_horizon ) {
  MinPower.resize( boost::extents[f_time_horizon][number_arcs] );

  for( Index n = 0; n < number_arcs; ++n ) {

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

 if( MaxPower.shape()[0] == 1 ) {
  MaxPower.resize( boost::extents[f_time_horizon][number_arcs] );

  for( Index n = 0; n < number_arcs; ++n ) {
  for( Index t = 0; t < f_time_horizon; ++t ) {

    MaxPower[t][n] = v_maximum_power[0][n];

   }
  }
 } else if( MaxPower.shape()[ 0 ] < f_time_horizon ) {

  MaxPower.resize( boost::extents[number_reservoirs][f_time_horizon] );

  for( Index n = 0; n < number_reservoirs; ++n ) {

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

 if( DeltaRampUp.shape()[0] == 1  ) {

  for( Index t = 0; t < f_time_horizon; ++t ) {
   for( Index g = 0; g < number_arcs; ++g ) {

    DeltaRampUp.resize( boost::extents[f_time_horizon][number_arcs] );

    DeltaRampUp[t][g] = v_delta_ramp_up[0][g];

   }
  }
 } else if( DeltaRampUp.shape()[ 0 ] < f_time_horizon ) {

  DeltaRampUp.resize( boost::extents[f_time_horizon][number_arcs] );

  for( Index g = 0; g < number_arcs; ++g ) {

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

 if( DeltaRampDown.size() == number_arcs ) {

  for( Index t = 0; t < f_time_horizon; ++t ) {
   for( Index g = 0; g < number_arcs; ++g ) {

    DeltaRampDown.resize( boost::extents[f_time_horizon][number_arcs] );

    DeltaRampDown[t][g] = v_delta_ramp_down[0][g];

   }
  }
 } else if( DeltaRampDown.size() > number_arcs ) {

  DeltaRampDown.resize( boost::extents[f_time_horizon][number_arcs] );

  for( Index g = 0; g < number_arcs; ++g ) {

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

 if( PrimaryRho.size() == number_arcs ) {

  for( Index t = 0; t < f_time_horizon; ++t ) {
   for( Index g = 0; g < number_arcs; ++g ) {

    PrimaryRho.resize( boost::extents[f_time_horizon][number_arcs] );

    PrimaryRho[t][g] = v_primary_rho[0][g];

   }
  }
 } else if( PrimaryRho.size() > number_arcs ) {

  PrimaryRho.resize( boost::extents[f_time_horizon][number_arcs] );

  for( Index g = 0; g < number_arcs; ++g ) {

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

 if( SecondaryRho.size() == number_arcs ) {

  for( Index t = 0; t < f_time_horizon; ++t ) {
   for( Index g = 0; g < number_arcs; ++g ) {

    SecondaryRho.resize( boost::extents[f_time_horizon][number_arcs] );

    SecondaryRho[t][g] = v_secondary_rho[0][g];

   }
  }
 } else if( SecondaryRho.size() > number_arcs ) {

  SecondaryRho.resize( boost::extents[f_time_horizon][number_arcs] );

  for( Index g = 0; g < number_arcs; ++g ) {

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
 std::vector< double > LinearTerm = v_linear_term;
 if( LinearTerm.size() == 1 ) {
  LinearTerm.resize( number_arcs, LinearTerm[0] );
 }
 // initial condition of vector ConstantTerm
 std::vector< double > ConstantTerm = v_const_term;
 if( ConstantTerm.size() == 1 ) {
  ConstantTerm.resize( number_arcs, ConstantTerm[0] );
 }
 // initial condition of matrix InertiaPower
 boost::multi_array< double, 2 > InertiaPower = v_inertia_power;

 if( InertiaPower.size() == number_arcs ) {

  for( Index t = 0; t < f_time_horizon; ++t ) {
   for( Index g = 0; g < number_arcs; ++g ) {

    InertiaPower.resize( boost::extents[f_time_horizon][number_arcs] );

    InertiaPower[t][g] = v_inertia_power[0][g];

   }
  }
 } else if( InertiaPower.size() > number_arcs ) {

  InertiaPower.resize( boost::extents[f_time_horizon][number_arcs] );

  for( Index g = 0; g < number_arcs; ++g ) {

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

 // initial condition of vector StartArc
 std::vector< Index > StartArc = v_start_arc;

 // initial condition of vector EndArc
 std::vector< Index > EndArc = v_end_arc;

/*--------------------------------------------------------------------------*/

 // maximum power output according to primary-secondary reserves constraints
  if( MaxPowerPrimarySecondary_Const.size() != f_time_horizon ) {
   // this should only happen once
   assert( MaxPowerPrimarySecondary_Const.empty());

   MaxPowerPrimarySecondary_Const.resize
           ( boost::multi_array< FRowConstraint, 2 >::
             extent_gen()[f_time_horizon][number_arcs] );
  }
 for( Index arc = 0; arc < number_arcs; ++arc ) {

  for( Index t = 0; t < f_time_horizon; ++t ) {
    auto linear_function = new LinearFunction();

    linear_function->add_variable( &v_active_power[arc][t], 1.0 );
    if ( !v_primary_rho.empty() ) {
     linear_function->add_variable( &v_primary_spinning_reserve[arc][t], 1.0 );
    }
    if ( !v_secondary_rho.empty() ) {
     linear_function->add_variable( &v_secondary_spinning_reserve[arc][t], 1.0 );
    }
    if (!v_minimum_power.empty() ) {
     MaxPowerPrimarySecondary_Const[t][arc].set_lhs( MinPower[t][arc] );
    } else {
     MaxPowerPrimarySecondary_Const[t][arc].set_lhs( 0.0 );

    }
   if (!v_maximum_power.empty() ) {
    MaxPowerPrimarySecondary_Const[t][arc].set_rhs( MaxPower[t][arc] );
   } else {
    MaxPowerPrimarySecondary_Const[t][arc].set_rhs( LinearTerm[0] * MaxFlow[t][arc] );

   }
    MaxPowerPrimarySecondary_Const[t][arc].set_function( linear_function );
   }
  }
  add_static_constraint( MaxPowerPrimarySecondary_Const, "MaxPowerPrimarySecondary" );


 // minimum power output according to primary-secondary reserves constraints

  if( MinPowerPrimarySecondary_Const.size() != f_time_horizon ) {
   // this should only happen once
   assert( MinPowerPrimarySecondary_Const.empty());

   MinPowerPrimarySecondary_Const.resize
           ( boost::multi_array< FRowConstraint, 2 >::
             extent_gen()[f_time_horizon][number_arcs] );
  }
  for( Index arc = 0; arc < number_arcs; ++arc ) {

   for( Index t = 0; t < f_time_horizon; ++t ) {

    auto linear_function = new LinearFunction();

    linear_function->add_variable( &v_active_power[arc][t], 1.0 );
    if ( !v_primary_rho.empty() ) {
     linear_function->add_variable( &v_primary_spinning_reserve[arc][t], -1.0 );
    }
    if ( !v_secondary_rho.empty() ) {
     linear_function->add_variable( &v_secondary_spinning_reserve[arc][t], -1.0 );
    }
    if (!v_minimum_power.empty() ) {
     MinPowerPrimarySecondary_Const[t][arc].set_lhs( MinPower[t][arc] );
    } else {
     MinPowerPrimarySecondary_Const[t][arc].set_lhs( 0.0 );

    }
    if (!v_maximum_power.empty() ) {
     MinPowerPrimarySecondary_Const[t][arc].set_rhs( MaxPower[t][arc] );
    } else {
     MinPowerPrimarySecondary_Const[t][arc].set_rhs( LinearTerm[0] * MaxFlow[t][arc] );

    }
    MinPowerPrimarySecondary_Const[t][arc].set_function( linear_function );
   }
  }

  add_static_constraint( MinPowerPrimarySecondary_Const, "MinPowerPrimarySecondary");

 // power output relation with to primary reserves constraints
 if( MaxFlow[0][0] > 0 ) {

  if( !v_primary_rho.empty()) {
   if( ActivePowerPrimary_Const.size() != f_time_horizon ) {
    // this should only happen once
    assert( ActivePowerPrimary_Const.empty());

    ActivePowerPrimary_Const.resize
            ( boost::multi_array< FRowConstraint, 2 >::
              extent_gen()[f_time_horizon][number_arcs] );
   }
   for( Index arc = 0; arc < number_arcs; ++arc ) {

    for( Index t = 0; t < f_time_horizon; ++t ) {

     auto linear_function = new LinearFunction();

     linear_function->add_variable( &v_active_power[arc][t], PrimaryRho[t][arc] );

     linear_function->add_variable( &v_primary_spinning_reserve[arc][t], -1.0 );
     ActivePowerPrimary_Const[t][arc].set_lhs( 0.0 );
     ActivePowerPrimary_Const[t][arc].set_rhs( Inf< double >());
     ActivePowerPrimary_Const[t][arc].set_function( linear_function );
    }

   }
   add_static_constraint( ActivePowerPrimary_Const, "ActivePowerPrimary" );

  }
 }
 // power output relation with to secondary reserves constraints
 if( MaxFlow[0][0] > 0 ) {
  if( !v_secondary_rho.empty()) {
   if( ActivePowerSecondary_Const.size() != f_time_horizon ) {
    // this should only happen once
    assert( ActivePowerSecondary_Const.empty());

    ActivePowerSecondary_Const.resize
            ( boost::multi_array< FRowConstraint, 2 >::
              extent_gen()[f_time_horizon][number_arcs] );
   }
   for( Index arc = 0; arc < number_arcs; ++arc ) {
    for( Index t = 0; t < f_time_horizon; ++t ) {
     auto linear_function = new LinearFunction();
     linear_function->add_variable( &v_active_power[arc][t], SecondaryRho[t][arc] );
     linear_function->add_variable( &v_secondary_spinning_reserve[arc][t], -1.0 );
     ActivePowerSecondary_Const[t][arc].set_lhs( 0.0 );
     ActivePowerSecondary_Const[t][arc].set_rhs( Inf< double >());
     ActivePowerSecondary_Const[t][arc].set_function( linear_function );
    }
   }
   add_static_constraint( ActivePowerSecondary_Const, "ActivePowerSecondary" );
  }
 }
 // primary reserves constraints for pumps
 if( MaxFlow[0][0] <= 0 ) {
  if( PrimaryPumps_Const.size() != f_time_horizon ) {
   // this should only happen once
   assert( PrimaryPumps_Const.empty());

   PrimaryPumps_Const.resize
           ( boost::multi_array< FRowConstraint, 2 >::
             extent_gen()[f_time_horizon][number_arcs] );
  }
  for( Index arc = 0; arc < number_arcs; ++arc ) {
   for( Index t = 0; t < f_time_horizon; ++t ) {


     auto linear_function = new LinearFunction();

     linear_function->add_variable( &v_primary_spinning_reserve[arc][t], 1.0 );
     PrimaryPumps_Const[t][arc].set_both( 0.0 );
     PrimaryPumps_Const[t][arc].set_function( linear_function );
    }

  }

  add_static_constraint( PrimaryPumps_Const, "PrimaryPumps");
 }
 // secondary reserves constraints for pumps
 if( MaxFlow[0][0] <= 0 ) {

  if( SecondaryPumps_Const.size() != f_time_horizon ) {
   // this should only happen once
   assert( SecondaryPumps_Const.empty());

   SecondaryPumps_Const.resize
           ( boost::multi_array< FRowConstraint, 2 >::
             extent_gen()[f_time_horizon][number_arcs] );
  }
  for( Index arc = 0; arc < number_arcs; ++arc ) {
   for( Index t = 0; t < f_time_horizon; ++t ) {

     auto linear_function = new LinearFunction();

     linear_function->add_variable( &v_secondary_spinning_reserve[arc][t], 1.0 );
     SecondaryPumps_Const[t][arc].set_both( 0.0 );
     SecondaryPumps_Const[t][arc].set_function( linear_function );
    }

  }

  add_static_constraint( SecondaryPumps_Const, "SecondaryPumps");
 }

 // flow to active power function constraints for pumps
 if( MaxFlow[0][0] <= 0 ) {

  if( FlowActivePowerPumps_Const.size() != f_time_horizon ) {
   // this should only happen once
   assert( FlowActivePowerPumps_Const.empty());

   FlowActivePowerPumps_Const.resize
           ( boost::multi_array< FRowConstraint, 2 >::
             extent_gen()[f_time_horizon][number_arcs] );
  }
  for( Index arc = 0; arc < number_arcs; ++arc ) {

  for( Index t = 0; t < f_time_horizon; ++t ) {

    auto linear_function = new LinearFunction();

    linear_function->add_variable( &v_active_power[arc][t], 1.0 );
    linear_function->add_variable( &v_flow_rate[arc][t], -LinearTerm[arc] );
    FlowActivePowerPumps_Const[t][arc].set_both( 0.0 );
    FlowActivePowerPumps_Const[t][arc].set_function( linear_function );

   }
  }
  add_static_constraint( FlowActivePowerPumps_Const, "FlowActivePowerPumps");
 }

 // flow to active power function constraints for turbines
 if (MaxFlow[0][0] > 0 ) {
  if (number_arcs > 0 ) {
  int TotalNumberPieces = f_total_number_pieces ? f_total_number_pieces : number_arcs;
  if( FlowActivePowerTurbines_Const.size() != f_time_horizon ) {
   // this should only happen once
   assert( FlowActivePowerTurbines_Const.empty());

   if( TotalNumberPieces == number_arcs ) {

    FlowActivePowerTurbines_Const.resize
            ( boost::multi_array< FRowConstraint, 3 >::
              extent_gen()[f_time_horizon][number_arcs][1] );
   } else if( TotalNumberPieces > number_arcs ) {
    FlowActivePowerTurbines_Const.resize
            ( boost::multi_array< FRowConstraint, 3 >::
              extent_gen()[f_time_horizon][1][TotalNumberPieces] );
   }
  }

  for( Index t = 0; t < f_time_horizon; ++t ) {

   Index piece = 0;
   Index end = 0;
   for( Index arc = 0; arc < number_arcs; ++arc ) {
    if( !NumberPieces.empty()) {
     end += NumberPieces[arc];
    }
    if( TotalNumberPieces == number_arcs ) {

     auto linear_function = new LinearFunction();

     linear_function->add_variable( &v_active_power[arc][t], 1.0 );
     if( !v_linear_term.empty()) {
      linear_function
              ->add_variable( &v_flow_rate[arc][t], -LinearTerm[arc] );
     } else {
      linear_function->add_variable( &v_flow_rate[arc][t], 0.0 );

     }
     if( !v_const_term.empty()) {
      FlowActivePowerTurbines_Const[t][arc][0]
              .set_rhs( ConstantTerm[arc] );
     } else {
      FlowActivePowerTurbines_Const[t][arc][0].set_rhs( 0.0 );

     }
     FlowActivePowerTurbines_Const[t][arc][0].set_lhs( -Inf< double >());
     FlowActivePowerTurbines_Const[t][arc][0]
             .set_function( linear_function );

    } else if( TotalNumberPieces > number_arcs ) {

     for( ; piece < end; ++piece ) {
      auto linear_function = new LinearFunction();

      linear_function->add_variable( &v_active_power[arc][t], 1.0 );
      if( !v_linear_term.empty()) {
       linear_function
               ->add_variable( &v_flow_rate[arc][t], -LinearTerm[piece] );
      } else {
       linear_function->add_variable( &v_flow_rate[arc][t], 0.0 );
      }
      if( !v_const_term.empty()) {
       FlowActivePowerTurbines_Const[t][0][piece]
               .set_rhs( ConstantTerm[piece] );
      } else {
       FlowActivePowerTurbines_Const[t][0][piece].set_rhs( 0.0 );
      }
      FlowActivePowerTurbines_Const[t][0][piece]
              .set_lhs( -Inf< double >());
      FlowActivePowerTurbines_Const[t][0][piece]
              .set_function( linear_function );
     }
    }
   }
  }
  add_static_constraint( FlowActivePowerTurbines_Const, "FlowActivePowerTurbines" );
 }
 }
 // flow rate bounds constraints
 {
  if( FlowRateBounds_Const.size() != f_time_horizon ) {
   // this should only happen once
   assert( FlowRateBounds_Const.empty());

   FlowRateBounds_Const.resize
           ( boost::multi_array< FRowConstraint, 2 >::
             extent_gen()[f_time_horizon][number_arcs] );
  }
  for( Index arc = 0; arc < number_arcs; ++arc ) {

   for( Index t = 0; t < f_time_horizon; ++t ) {

    auto linear_function = new LinearFunction();

    linear_function->add_variable( &v_flow_rate[arc][t], 1.0 );
    FlowRateBounds_Const[t][arc].set_lhs( 0.0 );
    FlowRateBounds_Const[t][arc].set_rhs( MaxFlow[t][arc] );
    FlowRateBounds_Const[t][arc].set_function( linear_function );

   }
  }

  add_static_constraint( FlowRateBounds_Const, "FlowRateBounds" );
 }


 // ram-up constraints
 if( !v_delta_ramp_up.empty() ) {

  if( RampUp_Const.size() != f_time_horizon ) {
   // this should only happen once
   assert( RampUp_Const.empty());

   RampUp_Const.resize
           ( boost::multi_array< FRowConstraint, 2 >::
             extent_gen()[f_time_horizon][number_arcs] );
  }
  // Initial condition

  for( Index arc = 0; arc < number_arcs; ++arc ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_flow_rate[arc][0], 1.0 );

   RampUp_Const[0][arc].set_lhs( -Inf< double >());
   RampUp_Const[0][arc].set_rhs( DeltaRampUp[0][arc] + InitialFlowRate[arc] );
   RampUp_Const[0][arc].set_function( linear_function );

  }

  for( Index t = 1; t < f_time_horizon; ++t ) {
   for( Index arc = 0; arc < number_arcs; ++arc ) {
    auto linear_function = new LinearFunction();

    linear_function->add_variable( &v_flow_rate[arc][t], 1.0 );
    linear_function->add_variable( &v_flow_rate[arc][t - 1], -1.0 );

    RampUp_Const[t][arc].set_lhs( -Inf< double >());
    RampUp_Const[t][arc].set_rhs( DeltaRampUp[t][arc] );
    RampUp_Const[t][arc].set_function( linear_function );

   }
  }

  add_static_constraint( RampUp_Const,"RampUp");
 }


 // ram-down constraints
 if( !v_delta_ramp_down.empty()) {

  if( RampDown_Const.size() != f_time_horizon ) {
   // this should only happen once
   assert( RampDown_Const.empty());

   RampDown_Const.resize
           ( boost::multi_array< FRowConstraint, 2 >::
             extent_gen()[f_time_horizon][number_arcs] );
  }

  // Initial condition

  for( Index arc = 0; arc < number_arcs; ++arc ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_flow_rate[arc][0], 1.0 );

   RampDown_Const[0][arc].set_lhs( InitialFlowRate[arc] - DeltaRampDown[0][arc] );
   RampDown_Const[0][arc].set_rhs( Inf< double >());
   RampDown_Const[0][arc].set_function( linear_function );

  }

  for( Index arc = 0; arc < number_arcs; ++arc ) {
  for( Index t = 1; t < f_time_horizon; ++t ) {
    auto linear_function = new LinearFunction();

    linear_function->add_variable( &v_flow_rate[arc][t - 1 ], 1.0 );
    linear_function->add_variable( &v_flow_rate[arc][t], -1.0 );

    RampDown_Const[t][arc].set_lhs( -Inf< double >());
    RampDown_Const[t][arc].set_rhs( DeltaRampDown[t][arc] );
    RampDown_Const[t][arc].set_function( linear_function );

   }
  }

  add_static_constraint( RampDown_Const, "RampDown");

 }

 // final volumes fo each reservoir constraints

 if( FinalVolumeReservoir_Const.size() != f_time_horizon ) {
  // this should only happen once
  assert( FinalVolumeReservoir_Const.empty());
  FinalVolumeReservoir_Const.resize
          ( boost::multi_array< FRowConstraint, 2 >::
            extent_gen()[f_time_horizon][ number_reservoirs ] );
 }

 for( Index n = 0; n < number_reservoirs; ++n ) {

  auto l_f = new LinearFunction();

  l_f->add_variable( &v_volumetric[n][0], 1.0 );

  for( Index l = 0; l < number_arcs; ++l ) {

   if( !v_start_arc.empty() && !v_end_arc.empty()) {

    if( StartArc[l] == n ) {

     if( StartArc[l] < EndArc[l] ) {

      if( !v_uphill_delay.empty()) {

       l_f->add_variable( &v_flow_rate[l][0 - UphillFlow[l]], 1.0 );

      } else {

       l_f->add_variable( &v_flow_rate[l][0], 1.0 );

      }

     } else if( StartArc[l] > EndArc[l] ) {

      if( !v_downhill_delay.empty()) {

       l_f->add_variable( &v_flow_rate[l][0 - DownhillFlow[l]], -1.0 );

      } else {

       l_f->add_variable( &v_flow_rate[l][0], -1.0 );

      }
     }
    }
    } else {

     l_f->add_variable( &v_flow_rate[n][0], 1.0 );

    }

  }
  FinalVolumeReservoir_Const[0][n].set_both( InitialVolumetric[n] +  Inflows[n][0] );
  FinalVolumeReservoir_Const[0][n].set_function( l_f );

  for( Index t = 1, constraint_index = 1; t < f_time_horizon;
       ++t, ++constraint_index  ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_volumetric[n][t], 1.0 );
   linear_function->add_variable( &v_volumetric[n][t-1], -1.0 );

  // linear_function->add_variable( &v_flow_rate[n][t], 1.0 );

   for( Index l = 0; l < number_arcs; ++l ) {

    if( !v_start_arc.empty() && !v_end_arc.empty()) {

     if( StartArc[l] == n ) {

      if( StartArc[l] < EndArc[l] ) {

      if( !v_uphill_delay.empty()) {

       linear_function->add_variable( &v_flow_rate[l][t - UphillFlow[l]], 1.0 );

      } else {

       linear_function->add_variable( &v_flow_rate[l][t], 1.0 );

      }
     } else if( StartArc[l] > EndArc[l] ) {

      if( !v_downhill_delay.empty()) {

       linear_function->add_variable( &v_flow_rate[l][t - DownhillFlow[l]], -1.0 );

      } else {

       linear_function->add_variable( &v_flow_rate[l][t], -1.0 );
      }
      }
     }
    } else {

     linear_function->add_variable( &v_flow_rate[n][t], 1.0 );

    }
   }
   FinalVolumeReservoir_Const[constraint_index][n].set_both(  Inflows[n][constraint_index] );
   FinalVolumeReservoir_Const[constraint_index][n].set_function( linear_function );
  }
  }

 add_static_constraint( FinalVolumeReservoir_Const, "FinalVolumeReservoir");



 // volumetric bounds constraints

 if( VolumetricBounds_Const.size() != f_time_horizon ) {
  // this should only happen once
  assert( VolumetricBounds_Const.empty());

  VolumetricBounds_Const.resize
          ( boost::multi_array< FRowConstraint, 2 >::
            extent_gen()[number_reservoirs][f_time_horizon] );
 }
 for( Index node = 0; node < number_reservoirs; ++node ) {
  for( Index t = 0; t < f_time_horizon; ++t ) {

   auto linear_function = new LinearFunction();
   linear_function->add_variable( &v_volumetric[node][t], 1.0 );

   VolumetricBounds_Const[node][t].set_lhs( MinVolumetric[node][t] );
   VolumetricBounds_Const[node][t].set_rhs( MaxVolumetric[node][t] );
   VolumetricBounds_Const[node][t].set_function( linear_function );

  }
 }

 add_static_constraint( VolumetricBounds_Const, "VolumetricBounds");

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

 auto dim_number_reservoirs = group.addDim( "NumberReservoirs",
         f_number_reservoirs ? f_number_reservoirs : 1 );

 auto dim_number_arcs = group.addDim( "NumberArcs",
         f_number_arcs ? f_number_arcs : 1 );


 ::serialize( group, "NumberPieces", netCDF::NcUint64(),
              dim_number_arcs, v_number_pieces, true );

  ::serialize( group, "StartLine", netCDF::NcInt64(),
               dim_number_reservoirs, v_start_arc, false );

  ::serialize( group, "EndLine", netCDF::NcInt64(),
               dim_number_reservoirs, v_end_arc, false );


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
