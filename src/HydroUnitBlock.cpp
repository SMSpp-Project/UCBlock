/*--------------------------------------------------------------------------*/
/*--------------------- File HydroUnitBlock.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the HydroUnitBlock class.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Ali Ghezelsoflu \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Rafael Durbano Lobato \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Ali Ghezelsoflu,
 *                      Rafael Durbano Lobato
 */

/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "HydroUnitBlock.h"

#include "LinearFunction.h"

#include "FRowConstraint.h"

#include "FRealObjective.h"

#include "OneVarConstraint.h"

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

HydroUnitBlock::~HydroUnitBlock()
{
 Constraint::clear( MaxPowerPrimarySecondary_Const );
 Constraint::clear( MinPowerPrimarySecondary_Const );
 Constraint::clear( ActivePowerPrimary_Const );
 Constraint::clear( ActivePowerSecondary_Const );
 Constraint::clear( FlowActivePower_Const );
 Constraint::clear( ActivePowerBounds_Const );
 Constraint::clear( RampUp_Const );
 Constraint::clear( RampDown_Const );
 Constraint::clear( FinalVolumeReservoir_Const );

 Constraint::clear( FlowRateBounds_Const );
 Constraint::clear( VolumetricBounds_Const );

 objective.clear();
}

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void HydroUnitBlock::deserialize( const netCDF::NcGroup & group )
{

#ifndef NDEBUG
 std::vector< std::string > expected_dims = { "TimeHorizon" ,
                                              "NumberIntervals" ,
                                              "NumberReservoirs" ,
                                              "NumberArcs" ,
                                              "TotalNumberPieces" };
 check_dimensions( group , expected_dims , std::cerr );

 std::vector< std::string > expected_vars = { "StartArc" ,
                                              "EndArc" ,
                                              "MinFlow" ,
                                              "MaxFlow" ,
                                              "MinVolumetric" ,
                                              "MaxVolumetric" ,
                                              "Inflows" ,
                                              "MinPower" ,
                                              "MaxPower" ,
                                              "DeltaRampUp" ,
                                              "DeltaRampDown" ,
                                              "PrimaryRho" ,
                                              "SecondaryRho" ,
                                              "NumberPieces" ,
                                              "LinearTerm" ,
                                              "ConstantTerm" ,
                                              "InertiaPower" ,
                                              "InitialFlowRate" ,
                                              "InitialVolumetric" ,
                                              "UphillFlow" ,
                                              "DownhillFlow" };
 check_variables( group , expected_vars , std::cerr );
#endif

 UnitBlock::deserialize_time_horizon( group );
 UnitBlock::deserialize_change_intervals( group );

 if( ! ::deserialize_dim( group , "NumberReservoirs" , f_NumberReservoirs ) )
  f_NumberReservoirs = 1;

 if( ! ::deserialize_dim( group , "NumberArcs" , f_NumberArcs ) )
  f_NumberArcs = 1;

 ::deserialize( group , "NumberPieces" , f_NumberArcs ,
                v_NumberPieces , true , true );

 if( ! ::deserialize_dim( group , "TotalNumberPieces" , f_TotalNumberPieces ) ) {
  f_TotalNumberPieces = 0;
  for( const auto & n : v_NumberPieces )
   f_TotalNumberPieces += n;
 }
 f_TotalNumberPieces = f_TotalNumberPieces ?
                       f_TotalNumberPieces : f_NumberArcs;

 ::deserialize( group , "StartArc" , f_NumberArcs , v_StartArc );
 ::deserialize( group , "EndArc" , f_NumberArcs , v_EndArc );

 ::deserialize( group , "Inflows" , v_inflows , true , false );

 ::deserialize( group , "MinFlow" , v_MinFlow , true , true );
 transpose( v_MinFlow );

 ::deserialize( group , "MaxFlow" , v_MaxFlow , true , true );
 transpose( v_MaxFlow );

 ::deserialize( group , "MinPower" , v_MinPower , true , true );
 transpose( v_MinPower );

 ::deserialize( group , "MaxPower" , v_MaxPower , true , true );
 transpose( v_MaxPower );

 ::deserialize( group , "DeltaRampUp" , v_DeltaRampUp , true , true );
 transpose( v_DeltaRampUp );

 ::deserialize( group , "DeltaRampDown" , v_DeltaRampDown , true , true );
 transpose( v_DeltaRampDown );

 ::deserialize( group , "PrimaryRho" , v_PrimaryRho , true , true );
 transpose( v_PrimaryRho );

 ::deserialize( group , "SecondaryRho" , v_SecondaryRho , true , true );
 transpose( v_SecondaryRho );

 ::deserialize( group , "LinearTerm" , f_TotalNumberPieces ,
                v_LinearTerm , true , true );

 ::deserialize( group , "ConstantTerm" , f_TotalNumberPieces ,
                v_ConstTerm , true , true );

 ::deserialize( group , "InertiaPower" , v_InertiaPower , true , true );
 transpose( v_InertiaPower );

 ::deserialize( group , "InitialFlowRate" , f_NumberArcs ,
                v_InitialFlowRate , true , true );

 ::deserialize( group , "InitialVolumetric" , f_NumberReservoirs ,
                v_InitialVolumetric , true , true );

 ::deserialize( group , "UphillFlow" , f_NumberArcs ,
                v_UphillDelay , true , true );

 ::deserialize( group , "DownhillFlow" , f_NumberArcs ,
                v_DownhillDelay , true , true );

 ::deserialize( group , "MinVolumetric" , v_MinVolumetric , true , true );
 ::deserialize( group , "MaxVolumetric" , v_MaxVolumetric , true , true );

 decompress_array( v_MinFlow );
 decompress_array( v_MaxFlow );
 decompress_vol( v_MinVolumetric );
 decompress_vol( v_MaxVolumetric );
 decompress_vol( v_inflows );

 decompress_array( v_MinPower );
 decompress_array( v_MaxPower );
 decompress_array( v_DeltaRampUp );
 decompress_array( v_DeltaRampDown );
 decompress_array( v_PrimaryRho );
 decompress_array( v_SecondaryRho );
 decompress_array( v_InertiaPower );

 if( v_LinearTerm.size() == 1 )
  v_LinearTerm.resize( f_NumberArcs , v_LinearTerm[ 0 ] );

 if( v_ConstTerm.size() == 1 )
  v_ConstTerm.resize( f_NumberArcs , v_ConstTerm[ 0 ] );

 UnitBlock::deserialize( group );

}  // end( HydroUnitBlock::deserialize )

/*--------------------------------------------------------------------------*/

void HydroUnitBlock::generate_abstract_variables( Configuration * stvv )
{
 if( variables_generated() )  // variables have already been generated
  return;                     // nothing to do

 UnitBlock::generate_abstract_variables( stvv );

 if( f_time_horizon == 0 )
  // there are no variables to be generated
  return;

 v_volumetric.resize( boost::extents[ f_NumberReservoirs ][ f_time_horizon ] );
 for( Index g = 0 ; g < f_NumberReservoirs ; ++g )
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   v_volumetric[ g ][ t ].set_type( ColVariable::kNonNegative );
 add_static_variable( v_volumetric , "v_hydro" );

 v_flow_rate.resize( boost::extents[ f_NumberArcs ][ f_time_horizon ] );
 v_active_power.resize( boost::extents[ f_NumberArcs ][ f_time_horizon ] );


 for( Index g = 0 ; g < f_NumberArcs ; ++g ) {
  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   v_flow_rate[ g ][ t ].set_type( ColVariable::kContinuous );
   v_active_power[ g ][ t ].set_type( ColVariable::kContinuous );
  }
 }

 add_static_variable( v_flow_rate , "f_hydro" );
 add_static_variable( v_active_power , "p_hydro" );
 if( reserve_vars & 1u ) {  // if UCBlock has primary demand variables
  if( ! v_PrimaryRho.empty() ) {  // if unit produces any primary reserve
   v_primary_spinning_reserve.resize(
    boost::extents[ f_NumberArcs ][ f_time_horizon ] );
   for( Index g = 0 ; g < f_NumberArcs ; ++g ) {
    for( Index t = 0 ; t < f_time_horizon ; ++t ) {
     v_primary_spinning_reserve[ g ][ t ].set_type( ColVariable::kNonNegative );
    }
   }
   add_static_variable( v_primary_spinning_reserve , "pr_hydro" );
  }
 }

 if( reserve_vars & 2u ) {  // if UCBlock has secondary demand variables
  if( ! v_SecondaryRho.empty() ) {  // if unit produces any secondary reserve
   v_secondary_spinning_reserve.resize(
    boost::extents[ f_NumberArcs ][ f_time_horizon ] );
   for( Index g = 0 ; g < f_NumberArcs ; ++g ) {
    for( Index t = 0 ; t < f_time_horizon ; ++t ) {
     v_secondary_spinning_reserve[ g ][ t ].set_type(
      ColVariable::kNonNegative );
    }
   }
   add_static_variable( v_secondary_spinning_reserve , "sr_hydro" );
  }
 }

 set_variables_generated();

}  // end( HydroUnitBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void HydroUnitBlock::generate_abstract_constraints( Configuration * stcc )
{
 if( constraints_generated() )  // constraints have already been generated
  return;                       // nothing to do

 LinearFunction::v_coeff_pair vars;

 // final volume constraints for each reservoir

 assert( FinalVolumeReservoir_Const.empty() );

 FinalVolumeReservoir_Const.resize(
  boost::multi_array< FRowConstraint , 2 >::
  extent_gen()[ f_time_horizon ][ f_NumberReservoirs ] );

 for( Index n = 0 ; n < f_NumberReservoirs ; ++n ) {

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {

   for( Index l = 0 ; l < f_NumberArcs ; ++l ) {

    if( ( ! v_StartArc.empty() ) && ( ! v_EndArc.empty() ) ) {

     const auto uphill_delay = get_uphill_delay( l );
     if( ( t >= uphill_delay ) && ( t - uphill_delay < f_time_horizon ) &&
         ( v_StartArc[ l ] == n ) )
      vars.push_back( std::make_pair( get_flow_rate( l , t - uphill_delay ) ,
                                      1.0 ) );

     const auto downhill_delay = get_downhill_delay( l );
     if( ( t >= downhill_delay ) && ( v_EndArc[ l ] == n ) )
      vars.push_back( std::make_pair( get_flow_rate( l , t - downhill_delay ) ,
                                      -1.0 ) );

    } else {

     vars.push_back( std::make_pair( get_flow_rate( l , t ) , 1.0 ) );
    }
   }

   vars.push_back( std::make_pair( get_volume( n , t ) , 1.0 ) );

   if( t > 0 )
    vars.push_back( std::make_pair( get_volume( n , t - 1 ) , -1.0 ) );

   double initial_volume = 0.0;
   if( t == 0 )
    initial_volume = v_InitialVolumetric[ n ];

   if( ! v_inflows.empty() )
    FinalVolumeReservoir_Const[ t ][ n ].set_both(
     initial_volume + v_inflows[ n ][ t ] );
   else
    FinalVolumeReservoir_Const[ t ][ n ].set_both( initial_volume );

   FinalVolumeReservoir_Const[ t ][ n ].set_function(
    new LinearFunction( std::move( vars ) ) );
  }
 }

 add_static_constraint( FinalVolumeReservoir_Const ,
                        "FinalVolumeReservoir_HydroUnit" );

 // maximum power output according to primary-secondary reserves constraints

 // Initial data check
 if( ( ! v_MinPower.empty() ) && ( ! v_MaxPower.empty() ) )
  for( Index arc = 0 ; arc < f_NumberArcs ; ++arc )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    if( v_MinPower[ t ][ arc ] > v_MaxPower[ t ][ arc ] )
     throw( std::logic_error( "HydroUnitBlock::maximum and minimum power "
                              "output constraints: it must be that v_MaxPower"
                              " >= v_MinPower." ) );

 // Initial data check
 if( ( ! v_MinFlow.empty() ) && ( ! v_MaxFlow.empty() ) )
  for( Index arc = 0 ; arc < f_NumberArcs ; ++arc )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    if( ( v_MaxFlow[ t ][ arc ] <= 0 ) &&
        ( v_MinFlow[ t ][ arc ] < 0 ) &&
        ( v_NumberPieces[ arc ] > 1 ) )
     throw( std::logic_error( "HydroUnitBlock::Data Error: it must be that "
                              "for each pump when v_MinPower < 0, then "
                              "v_NumberPieces == 1." ) );

 // Initial data check

 if( ( ! v_MinFlow.empty() ) && ( ! v_MaxFlow.empty() ) )
  for( Index arc = 0 ; arc < f_NumberArcs ; ++arc )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    if( ! v_PrimaryRho.empty() )
     if( ( v_MaxFlow[ t ][ arc ] <= 0 ) &&
         ( v_MinFlow[ t ][ arc ] < 0 ) &&
         ( v_PrimaryRho[ t ][ arc ] != 0 ) )
      throw( std::logic_error( "HydroUnitBlock::Data Error: it must be that "
                               "for each pump v_PrimaryRho == 0." ) );

 // Initial data check
 if( ( ! v_MinFlow.empty() ) && ( ! v_MaxFlow.empty() ) )
  for( Index arc = 0 ; arc < f_NumberArcs ; ++arc )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    if( ! v_SecondaryRho.empty() )
     if( ( v_MaxFlow[ t ][ arc ] <= 0 ) &&
         ( v_MinFlow[ t ][ arc ] < 0 ) &&
         ( v_SecondaryRho[ t ][ arc ] != 0 ) )
      throw( std::logic_error( "HydroUnitBlock::Data Error: it must be that "
                               "for each pump then v_SecondaryRho == 0." ) );

 assert( MaxPowerPrimarySecondary_Const.empty() );

 MaxPowerPrimarySecondary_Const.resize(
  boost::multi_array< FRowConstraint , 2 >::
  extent_gen()[ f_time_horizon ][ f_NumberArcs ] );

 for( Index arc = 0 ; arc < f_NumberArcs ; ++arc ) {

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {

   vars.push_back( std::make_pair( get_active_power( arc , t ) , 1.0 ) );

   if( reserve_vars & 1u )  // if UCBlock has primary demand variables
    if( ! v_PrimaryRho.empty() )  // if unit produces any primary reserve
     vars.push_back( std::make_pair( get_primary_spinning_reserve( arc , t ) ,
                                     1.0 ) );

   if( reserve_vars & 2u )  // if UCBlock has secondary demand variables
    if( ! v_SecondaryRho.empty() )  // if unit produces any secondary reserve
     vars.push_back( std::make_pair( get_secondary_spinning_reserve( arc , t ) ,
                                     1.0 ) );

   MaxPowerPrimarySecondary_Const[ t ][ arc ].set_lhs( -Inf< double >() );

   if( ! v_MaxPower.empty() )
    MaxPowerPrimarySecondary_Const[ t ][ arc ].set_rhs(
     v_MaxPower[ t ][ arc ] );
   else
    MaxPowerPrimarySecondary_Const[ t ][ arc ].set_rhs( 0.0 );
   MaxPowerPrimarySecondary_Const[ t ][ arc ].set_function(
    new LinearFunction( std::move( vars ) ) );
  }
 }

 add_static_constraint( MaxPowerPrimarySecondary_Const ,
                        "MaxPowerPrimarySecondary_HydroUnit" );

 // minimum power output according to primary-secondary reserves constraints

 assert( MinPowerPrimarySecondary_Const.empty() );

 MinPowerPrimarySecondary_Const.resize(
  boost::multi_array< FRowConstraint , 2 >::
  extent_gen()[ f_time_horizon ][ f_NumberArcs ] );

 for( Index arc = 0 ; arc < f_NumberArcs ; ++arc ) {

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {

   vars.push_back( std::make_pair( get_active_power( arc , t ) , 1.0 ) );

   if( reserve_vars & 1u )  // if UCBlock has primary demand variables
    if( ! v_PrimaryRho.empty() )  // if unit produces any primary reserve
     vars.push_back( std::make_pair( get_primary_spinning_reserve( arc , t ) ,
                                     -1.0 ) );

   if( reserve_vars & 2u )  // if UCBlock has secondary demand variables
    if( ! v_SecondaryRho.empty() )  // if unit produces any secondary reserve
     vars.push_back( std::make_pair( get_secondary_spinning_reserve( arc , t ) ,
                                     -1.0 ) );

   if( ! v_MinPower.empty() )
    MinPowerPrimarySecondary_Const[ t ][ arc ].set_lhs(
     v_MinPower[ t ][ arc ] );
   else
    MinPowerPrimarySecondary_Const[ t ][ arc ].set_lhs( 0.0 );
   MinPowerPrimarySecondary_Const[ t ][ arc ].set_rhs( Inf< double >() );
   MinPowerPrimarySecondary_Const[ t ][ arc ].set_function(
    new LinearFunction( std::move( vars ) ) );
  }
 }

 add_static_constraint( MinPowerPrimarySecondary_Const ,
                        "MinPowerPrimarySecondary_HydroUnit" );

 // power output relation with to primary reserves constraints

 if( reserve_vars & 1u ) {  // if UCBlock has primary demand variables
  if( ! v_PrimaryRho.empty() ) {  // if unit produces any primary reserve

   assert( ActivePowerPrimary_Const.empty() );

   ActivePowerPrimary_Const.resize(
    boost::multi_array< FRowConstraint , 2 >::
    extent_gen()[ f_time_horizon ][ f_NumberArcs ] );

   if( ( ! v_MinFlow.empty() ) && ( ! v_MaxFlow.empty() ) ) {

    for( Index arc = 0 ; arc < f_NumberArcs ; ++arc ) {

     for( Index t = 0 ; t < f_time_horizon ; ++t ) {

      if( ( v_MinFlow[ t ][ arc ] >= 0 ) &&
          ( v_MaxFlow[ t ][ arc ] > 0 ) ) {  // Turbines

       vars.push_back( std::make_pair( get_active_power( arc , t ) ,
                                       v_PrimaryRho[ t ][ arc ] ) );
       vars.push_back( std::make_pair( get_primary_spinning_reserve( arc , t ) ,
                                       -1.0 ) );

       ActivePowerPrimary_Const[ t ][ arc ].set_lhs( 0.0 );
       ActivePowerPrimary_Const[ t ][ arc ].set_rhs( Inf< double >() );
       ActivePowerPrimary_Const[ t ][ arc ].set_function(
        new LinearFunction( std::move( vars ) ) );
      }

      if( ( v_MaxFlow[ t ][ arc ] <= 0 ) &&
          ( v_MinFlow[ t ][ arc ] < 0 ) ) {  // Pumps

       vars.push_back( std::make_pair( get_primary_spinning_reserve( arc , t ) ,
                                       1.0 ) );
       ActivePowerPrimary_Const[ t ][ arc ].set_both( 0.0 );
       ActivePowerPrimary_Const[ t ][ arc ].set_function(
        new LinearFunction( std::move( vars ) ) );
      }

      if( ( v_MaxFlow[ t ][ arc ] == 0 ) &&
          ( v_MinFlow[ t ][ arc ] == 0 ) ) {  // Nothing

       vars.push_back( std::make_pair( get_flow_rate( arc , t ) , 1.0 ) );
       ActivePowerPrimary_Const[ t ][ arc ].set_both( 0.0 );
       ActivePowerPrimary_Const[ t ][ arc ].set_function(
        new LinearFunction( std::move( vars ) ) );
      }
     }
    }
   }

   if( ( v_MinFlow.empty() ) && ( ! v_MaxFlow.empty() ) ) {

    for( Index arc = 0 ; arc < f_NumberArcs ; ++arc ) {

     for( Index t = 0 ; t < f_time_horizon ; ++t ) {

      if( v_MaxFlow[ t ][ arc ] > 0 ) {  // Turbines

       vars.push_back( std::make_pair( get_active_power( arc , t ) ,
                                       v_PrimaryRho[ t ][ arc ] ) );
       vars.push_back( std::make_pair( get_primary_spinning_reserve( arc , t ) ,
                                       -1.0 ) );

       ActivePowerPrimary_Const[ t ][ arc ].set_lhs( 0.0 );
       ActivePowerPrimary_Const[ t ][ arc ].set_rhs( Inf< double >() );
       ActivePowerPrimary_Const[ t ][ arc ].set_function(
        new LinearFunction( std::move( vars ) ) );
      }

      if( v_MaxFlow[ t ][ arc ] == 0 ) {  // Nothing

       vars.push_back( std::make_pair( get_flow_rate( arc , t ) , 1.0 ) );

       ActivePowerPrimary_Const[ t ][ arc ].set_both( 0.0 );
       ActivePowerPrimary_Const[ t ][ arc ].set_function(
        new LinearFunction( std::move( vars ) ) );
      }
     }
    }
   }

   add_static_constraint( ActivePowerPrimary_Const ,
                          "ActivePowerPrimary_HydroUnit" );
  }
 }

 // power output relation with to secondary reserves constraints
 if( reserve_vars & 2u ) {  // if UCBlock has secondary demand variables
  if( ! v_SecondaryRho.empty() ) {  // if unit produces any secondary reserve

   assert( ActivePowerSecondary_Const.empty() );

   ActivePowerSecondary_Const.resize(
    boost::multi_array< FRowConstraint , 2 >::
    extent_gen()[ f_time_horizon ][ f_NumberArcs ] );

   if( ( ! v_MinFlow.empty() ) && ( ! v_MaxFlow.empty() ) ) {

    for( Index arc = 0 ; arc < f_NumberArcs ; ++arc ) {

     for( Index t = 0 ; t < f_time_horizon ; ++t ) {

      if( ( v_MinFlow[ t ][ arc ] >= 0 ) &&
          ( v_MaxFlow[ t ][ arc ] > 0 ) ) {  // Turbines

       vars.push_back( std::make_pair( get_active_power( arc , t ) ,
                                       v_SecondaryRho[ t ][ arc ] ) );
       vars.push_back( std::make_pair( get_secondary_spinning_reserve( arc , t ) ,
                                       -1.0 ));

       ActivePowerSecondary_Const[ t ][ arc ].set_lhs( 0.0 );
       ActivePowerSecondary_Const[ t ][ arc ].set_rhs( Inf< double >() );
       ActivePowerSecondary_Const[ t ][ arc ].set_function(
        new LinearFunction( std::move( vars ) ) );
      }

      if( ( v_MaxFlow[ t ][ arc ] <= 0 ) &&
          ( v_MinFlow[ t ][ arc ] < 0 ) ) {  // Pumps

       vars.push_back( std::make_pair( get_secondary_spinning_reserve( arc , t ) ,
                                       1.0 ) );

       ActivePowerSecondary_Const[ t ][ arc ].set_both( 0.0 );
       ActivePowerSecondary_Const[ t ][ arc ].set_function(
        new LinearFunction( std::move( vars ) ) );
      }

      if( ( v_MaxFlow[ t ][ arc ] == 0 ) &&
          ( v_MinFlow[ t ][ arc ] == 0 ) ) {  // Nothing

       vars.push_back( std::make_pair( get_flow_rate( arc , t ) , 1.0 ) );

       ActivePowerSecondary_Const[ t ][ arc ].set_both( 0.0 );
       ActivePowerSecondary_Const[ t ][ arc ].set_function(
        new LinearFunction( std::move( vars ) ) );
      }
     }
    }
   }

   if( ( v_MinFlow.empty() ) && ( ! v_MaxFlow.empty() ) ) {

    for( Index arc = 0 ; arc < f_NumberArcs ; ++arc ) {

     for( Index t = 0 ; t < f_time_horizon ; ++t ) {

      if( v_MaxFlow[ t ][ arc ] > 0 ) {  // Turbines

       vars.push_back( std::make_pair( get_active_power( arc , t ) ,
                                       v_SecondaryRho[ t ][ arc ] ) );
       vars.push_back( std::make_pair( get_secondary_spinning_reserve( arc , t ) ,
                                       -1.0 ) );

       ActivePowerSecondary_Const[ t ][ arc ].set_lhs( 0.0 );
       ActivePowerSecondary_Const[ t ][ arc ].set_rhs( Inf< double >() );
       ActivePowerSecondary_Const[ t ][ arc ].set_function(
        new LinearFunction( std::move( vars ) ) );
      }

      if( v_MaxFlow[ t ][ arc ] == 0 ) {  // Nothing

       vars.push_back( std::make_pair( get_flow_rate( arc , t ) ,
                                       1.0 ) );

       ActivePowerSecondary_Const[ t ][ arc ].set_both( 0.0 );
       ActivePowerSecondary_Const[ t ][ arc ].set_function(
        new LinearFunction( std::move( vars ) ) );
      }
     }
    }
   }

   add_static_constraint( ActivePowerSecondary_Const ,
                          "ActivePowerSecondary_HydroUnit" );
  }
 }

 assert( FlowActivePower_Const.empty() );

 FlowActivePower_Const.resize(
  boost::multi_array< FRowConstraint , 2 >::
  extent_gen()[ f_time_horizon ][ f_TotalNumberPieces ] );

 if( f_NumberArcs > 0 ) {

  if( ( ! v_MinFlow.empty() ) && ( ! v_MaxFlow.empty() ) ) {

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {

    Index piece = 0;
    Index cnstr_idx = 0;
    Index end = 0;

    for( Index arc = 0 ; arc < f_NumberArcs ; ++arc ) {

     if( ! v_NumberPieces.empty() )
      end += v_NumberPieces[ arc ];

     if( ( v_MinFlow[ t ][ arc ] >= 0 ) &&
         ( v_MaxFlow[ t ][ arc ] > 0 ) ) {  // Turbines

      for( ; piece < end ; ++piece , ++cnstr_idx ) {

       vars.push_back( std::make_pair( get_active_power( arc , t ) , 1.0 ) );

       if( ! v_LinearTerm.empty() )
        vars.push_back( std::make_pair( get_flow_rate( arc , t ) ,
                                        -v_LinearTerm[ piece ] ) );
       else
        vars.push_back( std::make_pair( get_flow_rate( arc , t ) , 0.0 ) );

       if( ! v_ConstTerm.empty() )
        FlowActivePower_Const[ t ][ piece ].set_rhs( v_ConstTerm[ piece ] );
       else
        FlowActivePower_Const[ t ][ piece ].set_rhs( 0.0 );
       FlowActivePower_Const[ t ][ piece ].set_lhs( -Inf< double >() );
       FlowActivePower_Const[ t ][ piece ].set_function(
        new LinearFunction( std::move( vars ) ) );
      }
     }

     if( ( v_MaxFlow[ t ][ arc ] <= 0 ) &&
         ( v_MinFlow[ t ][ arc ] < 0 ) ) {  // Pumps

      vars.push_back( std::make_pair( get_active_power( arc , t ) , 1.0 ) );
      vars.push_back( std::make_pair( get_flow_rate( arc , t ) ,
                                      -v_LinearTerm[ cnstr_idx ] ) );

      FlowActivePower_Const[ t ][ cnstr_idx ].set_both( 0.0 );
      FlowActivePower_Const[ t ][ cnstr_idx ].set_function(
       new LinearFunction( std::move( vars ) ) );

      ++cnstr_idx;
      piece = cnstr_idx;
     }

     if( ( v_MaxFlow[ t ][ arc ] == 0 ) &&
         ( v_MinFlow[ t ][ arc ] == 0 ) ) {  // Nothing

      vars.push_back( std::make_pair( get_flow_rate( arc , t ) , 1.0 ) );

      FlowActivePower_Const[ t ][ cnstr_idx ].set_both( 0.0 );
      FlowActivePower_Const[ t ][ cnstr_idx ].set_function(
       new LinearFunction( std::move( vars ) ) );

      ++cnstr_idx;
      piece = cnstr_idx;
     }
    }
   }
  }

  if( ( v_MinFlow.empty() ) && ( ! v_MaxFlow.empty() ) ) {

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {

    Index piece = 0;
    Index cnstr_idx = 0;
    Index end = 0;

    for( Index arc = 0 ; arc < f_NumberArcs ; ++arc ) {

     if( ! v_NumberPieces.empty() )
      end += v_NumberPieces[ arc ];
     else
      end = 1;

     if( v_MaxFlow[ t ][ arc ] > 0 ) {  // Turbines

      for( ; piece < end ; ++piece , ++cnstr_idx) {

       vars.push_back( std::make_pair( get_active_power( arc , t ) ,
                                       1.0 ) );

       if( ! v_LinearTerm.empty() )
        vars.push_back( std::make_pair( get_flow_rate( arc , t ) ,
                                        -v_LinearTerm[ piece ] ) );
       else
        vars.push_back( std::make_pair( get_flow_rate( arc , t ) , 0.0 ) );

       if( ! v_ConstTerm.empty() )
        FlowActivePower_Const[ t ][ piece ].set_rhs( v_ConstTerm[ piece ] );
       else
        FlowActivePower_Const[ t ][ piece ].set_rhs( 0.0 );
       FlowActivePower_Const[ t ][ piece ].set_lhs( -Inf< double >() );
       FlowActivePower_Const[ t ][ piece ].set_function(
        new LinearFunction( std::move( vars ) ) );
      }
     }

     if( v_MaxFlow[ t ][ arc ] == 0 ) {  // Nothing

      vars.push_back( std::make_pair( get_flow_rate( arc , t ) , 1.0 ) );

      FlowActivePower_Const[ t ][ cnstr_idx ].set_both( 0.0 );
      FlowActivePower_Const[ t ][ cnstr_idx ].set_function(
       new LinearFunction( std::move( vars ) ) );

      ++cnstr_idx;
      piece = cnstr_idx;
     }
    }
   }
  }

  add_static_constraint( FlowActivePower_Const , "FlowActivePower_HydroUnit" );
 }

 {
  // Initial data check
  if( ( ! v_MinFlow.empty() ) && ( ! v_MaxFlow.empty() ) )
   for( Index arc = 0 ; arc < f_NumberArcs ; ++arc )
    for( Index t = 0 ; t < f_time_horizon ; ++t )
     if( v_MinFlow[ t ][ arc ] > v_MaxFlow[ t ][ arc ] )
      throw( std::logic_error( "HydroUnitBlock::flow rate variable bounds: "
                               "it must be that v_MaxFlow >= v_MinFlow." ) );

  assert( FlowRateBounds_Const.empty() );

  FlowRateBounds_Const.resize(
   boost::multi_array< BoxConstraint , 2 >::
   extent_gen()[ f_time_horizon ][ f_NumberArcs ] );

  for( Index t = 0 ; t < f_time_horizon ; ++t )
   for( Index arc = 0 ; arc < f_NumberArcs ; ++arc ) {
    auto flow_rate = get_flow_rate( arc , t );
    if( ! v_MinFlow.empty() )
     FlowRateBounds_Const[ t ][ arc ].set_lhs( v_MinFlow[ t ][ arc ] );
    else
     FlowRateBounds_Const[ t ][ arc ].set_lhs( 0.0 );
    if( ! v_MaxFlow.empty() )
     FlowRateBounds_Const[ t ][ arc ].set_rhs( v_MaxFlow[ t ][ arc ] );
    else
     FlowRateBounds_Const[ t ][ arc ].set_rhs( 0.0 );
    FlowRateBounds_Const[ t ][ arc ].set_variable( flow_rate );
   }

  add_static_constraint( FlowRateBounds_Const , "FlowRateBounds_HydroUnit" );
 }

 // ramp-up constraints
 if( ! v_DeltaRampUp.empty() ) {

  assert( RampUp_Const.empty() );

  RampUp_Const.resize(
   boost::multi_array< FRowConstraint , 2 >::
   extent_gen()[ f_time_horizon ][ f_NumberArcs ] );

  // Initial condition

  for( Index arc = 0 ; arc < f_NumberArcs ; ++arc ) {

   vars.push_back( std::make_pair( get_flow_rate( arc , 0 ) , 1.0 ) );

   RampUp_Const[ 0 ][ arc ].set_lhs( -Inf< double >() );
   RampUp_Const[ 0 ][ arc ].set_rhs( v_DeltaRampUp[ 0 ][ arc ] +
                                     get_initial_flow_rate( arc ) );
   RampUp_Const[ 0 ][ arc ].set_function(
    new LinearFunction( std::move( vars ) ) );

   for( Index t = 1 , cnstr_idx = 1 ; t < f_time_horizon ; ++t , ++cnstr_idx ) {

    vars.push_back( std::make_pair( get_flow_rate( arc , t ) , 1.0 ) );
    vars.push_back( std::make_pair( get_flow_rate( arc , t - 1 ) , -1.0 ) );

    RampUp_Const[ cnstr_idx ][ arc ].set_lhs( -Inf< double >() );
    RampUp_Const[ cnstr_idx ][ arc ].set_rhs( v_DeltaRampUp[ t ][ arc ] );
    RampUp_Const[ cnstr_idx ][ arc ].set_function(
     new LinearFunction( std::move( vars ) ) );
   }
  }

  add_static_constraint( RampUp_Const , "RampUp_HydroUnit" );
 }

 // ramp-down constraints
 if( ! v_DeltaRampDown.empty() ) {

  assert( RampDown_Const.empty() );

  RampDown_Const.resize(
   boost::multi_array< FRowConstraint , 2 >::
   extent_gen()[ f_time_horizon ][ f_NumberArcs ] );

  // Initial condition

  for( Index arc = 0 ; arc < f_NumberArcs ; ++arc ) {

   vars.push_back( std::make_pair( get_flow_rate( arc , 0 ) , 1.0 ) );

   RampDown_Const[ 0 ][ arc ].set_lhs( get_initial_flow_rate( arc ) -
                                       v_DeltaRampDown[ 0 ][ arc ] );
   RampDown_Const[ 0 ][ arc ].set_rhs( Inf< double >() );
   RampDown_Const[ 0 ][ arc ].set_function(
    new LinearFunction( std::move( vars ) ) );

   for( Index t = 1 , cnstr_idx = 1 ; t < f_time_horizon ; ++t , ++cnstr_idx ) {

    vars.push_back( std::make_pair( get_flow_rate( arc , t - 1 ) , 1.0 ) );
    vars.push_back( std::make_pair( get_flow_rate( arc , t ) , -1.0 ) );

    RampDown_Const[ cnstr_idx ][ arc ].set_lhs( -Inf< double >() );
    RampDown_Const[ cnstr_idx ][ arc ].set_rhs(
     v_DeltaRampDown[ t ][ arc ] );
    RampDown_Const[ cnstr_idx ][ arc ].set_function(
     new LinearFunction( std::move( vars ) ) );
   }
  }

  add_static_constraint( RampDown_Const , "RampDown_HydroUnit" );
 }

 // volumetric bounds constraints

 // Initial data check
 if( ( ! v_MinVolumetric.empty() ) && ( ! v_MaxVolumetric.empty() ) )
  for( Index node = 0 ; node < f_NumberReservoirs ; ++node )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    if( ( v_MinVolumetric[ node ][ t ] > v_MaxVolumetric[ node ][ t ] ) ||
        ( v_MinVolumetric[ node ][ t ] < 0 ) ||
        ( v_MaxVolumetric[ node ][ t ] < 0 ) )
     throw( std::logic_error( "HydroUnitBlock::Volumetric Bounds Constraint: "
                              "it must be that 0 <= MinV[ r , t ] <= "
                              "MaxV[ r , t ] ." ) );

 assert( VolumetricBounds_Const.empty() );

 VolumetricBounds_Const.resize(
  boost::multi_array< BoxConstraint , 2 >::
  extent_gen()[ f_NumberReservoirs ][ f_time_horizon ] );

 for( Index node = 0 ; node < f_NumberReservoirs ; ++node )

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {

   if( ! v_MinVolumetric.empty() )
    VolumetricBounds_Const[ node ][ t ].set_lhs(
     v_MinVolumetric[ node ][ t ] );
   else
    VolumetricBounds_Const[ node ][ t ].set_lhs( 0.0 );
   if( ! v_MaxVolumetric.empty() )
    VolumetricBounds_Const[ node ][ t ].set_rhs(
     v_MaxVolumetric[ node ][ t ] );
   else
    VolumetricBounds_Const[ node ][ t ].set_rhs( 0.0 );
   VolumetricBounds_Const[ node ][ t ].set_variable( get_volume( node , t ) );
  }

 add_static_constraint( VolumetricBounds_Const , "VolumetricBounds_HydroUnit" );

 set_constraints_generated();

}  // end( HydroUnitBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void HydroUnitBlock::generate_objective( Configuration * objc )
{
 if( objective_generated() )  // Objective has already been generated
  return;                     // nothing to do

 objective.set_function( new LinearFunction() );

 // Set Block objective
 this->set_objective( &objective );

 set_objective_generated();

}  // end( HydroUnitBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*----------------- METHODS FOR CHECKING THE HydroUnitBlock ----------------*/
/*--------------------------------------------------------------------------*/

bool HydroUnitBlock::is_feasible( bool useabstract , Configuration * fsbc )
{
 // Retrieve the tolerance and the type of violation.
 double tol = 0;
 bool rel_viol = true;

 // Try to extract, from "c", the parameters that determine feasibility.
 // If it succeeds, it sets the values of the parameters and returns
 // true. Otherwise, it returns false.
 auto extract_parameters = [ & tol , & rel_viol ]( Configuration * c )
  -> bool {
  if( auto tc = dynamic_cast< SimpleConfiguration< double > * >( c ) ) {
   tol = tc->f_value;
   return( true );
  }
  if( auto tc = dynamic_cast< SimpleConfiguration< std::pair< double , int > > * >( c ) ) {
   tol = tc->f_value.first;
   rel_viol = tc->f_value.second;
   return( true );
  }
  return( false );
 };

 if( ( ! extract_parameters( fsbc ) ) && f_BlockConfig )
  // if the given Configuration is not valid, try the one from the BlockConfig
  extract_parameters( f_BlockConfig->f_is_feasible_Configuration );

 return(
  UnitBlock::is_feasible( useabstract )
  // Variables: Notice that there is no check for the v_flow_rate and
  // v_active_power variables, since they are continuous and have no bounds
  && ColVariable::is_feasible( v_volumetric , tol )
  && ColVariable::is_feasible( v_primary_spinning_reserve , tol )
  && ColVariable::is_feasible( v_secondary_spinning_reserve , tol )
  // Constraints
  && RowConstraint::is_feasible( MaxPowerPrimarySecondary_Const , tol , rel_viol )
  && RowConstraint::is_feasible( MinPowerPrimarySecondary_Const , tol , rel_viol )
  && RowConstraint::is_feasible( ActivePowerPrimary_Const , tol , rel_viol )
  && RowConstraint::is_feasible( ActivePowerSecondary_Const , tol , rel_viol )
  && RowConstraint::is_feasible( FlowActivePower_Const , tol , rel_viol )
  && RowConstraint::is_feasible( ActivePowerBounds_Const , tol , rel_viol )
  && RowConstraint::is_feasible( RampUp_Const , tol , rel_viol )
  && RowConstraint::is_feasible( RampDown_Const , tol , rel_viol )
  && RowConstraint::is_feasible( FlowRateBounds_Const , tol , rel_viol )
  && RowConstraint::is_feasible( FinalVolumeReservoir_Const , tol , rel_viol )
  && RowConstraint::is_feasible( VolumetricBounds_Const , tol , rel_viol ) );

} // end( HydroUnitBlock::is_feasible )

/*--------------------------------------------------------------------------*/
/*-------- METHODS FOR LOADING, PRINTING & SAVING THE HydroUnitBlock -------*/
/*--------------------------------------------------------------------------*/

void HydroUnitBlock::serialize( netCDF::NcGroup & group ) const
{
 UnitBlock::serialize( group );

 auto TimeHorizon = group.getDim( "TimeHorizon" );

 auto NumberIntervals = group.getDim( "NumberIntervals" );

 auto TotalNumberPieces = group.addDim( "TotalNumberPieces" ,
                                        f_TotalNumberPieces );

 auto NumberReservoirs = group.addDim( "NumberReservoirs" ,
                                       f_NumberReservoirs ? f_NumberReservoirs
                                                          : 1 );

 auto NumberArcs = group.addDim( "NumberArcs" ,
                                 f_NumberArcs ? f_NumberArcs : 1 );

 // Serialize one-dimensional variables

 ::serialize( group , "StartArc" , netCDF::NcUint() ,
              NumberArcs , v_StartArc , false );

 ::serialize( group , "EndArc" , netCDF::NcUint() ,
              NumberArcs , v_EndArc , false );

 ::serialize( group , "NumberPieces" , netCDF::NcUint() ,
              NumberArcs , v_NumberPieces , false );

 ::serialize( group , "LinearTerm" , netCDF::NcDouble() ,
              TotalNumberPieces , v_LinearTerm , false );

 ::serialize( group , "ConstantTerm" , netCDF::NcDouble() ,
              TotalNumberPieces , v_ConstTerm , false );

 ::serialize( group , "InitialFlowRate" , netCDF::NcDouble() ,
              NumberArcs , v_InitialFlowRate , false );

 ::serialize( group , "InitialVolumetric" , netCDF::NcDouble() ,
              NumberReservoirs , v_InitialVolumetric , false );

 ::serialize( group , "UphillFlow" , netCDF::NcInt() ,
              NumberArcs , v_UphillDelay , false );

 ::serialize( group , "DownhillFlow" , netCDF::NcUint() ,
              NumberArcs , v_DownhillDelay , false );

 // Serialize two-dimensional variables

 ::serialize( group , "MinFlow" , netCDF::NcDouble() ,
              { TimeHorizon , NumberArcs } , v_MinFlow );

 ::serialize( group , "MaxFlow" , netCDF::NcDouble() ,
              { TimeHorizon , NumberArcs } , v_MaxFlow );

 ::serialize( group , "MinVolumetric" , netCDF::NcDouble() ,
              { NumberReservoirs , TimeHorizon } , v_MinVolumetric );

 ::serialize( group , "MaxVolumetric" , netCDF::NcDouble() ,
              { NumberReservoirs , TimeHorizon } , v_MaxVolumetric );

 ::serialize( group , "Inflows" , netCDF::NcDouble() ,
              { NumberReservoirs , TimeHorizon } , v_inflows );

 ::serialize( group , "MinPower" , netCDF::NcDouble() ,
              { TimeHorizon , NumberArcs } , v_MinPower );

 ::serialize( group , "MaxPower" , netCDF::NcDouble() ,
              { TimeHorizon , NumberArcs } , v_MaxPower );

 ::serialize( group , "DeltaRampUp" , netCDF::NcDouble() ,
              { TimeHorizon , NumberArcs } , v_DeltaRampUp );

 ::serialize( group , "DeltaRampDown" , netCDF::NcDouble() ,
              { TimeHorizon , NumberArcs } , v_DeltaRampDown );

 ::serialize( group , "PrimaryRho" , netCDF::NcDouble() ,
              { TimeHorizon , NumberArcs } , v_PrimaryRho );

 ::serialize( group , "SecondaryRho" , netCDF::NcDouble() ,
              { TimeHorizon , NumberArcs } , v_SecondaryRho );

 ::serialize( group , "InertiaPower" , netCDF::NcDouble() ,
              { TimeHorizon , NumberArcs } , v_InertiaPower );

}  // end( HydroUnitBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

void HydroUnitBlock::set_inflow( MF_dbl_it values ,
                                 Block::Subset && subset ,
                                 const bool ordered ,
                                 c_ModParam issuePMod ,
                                 c_ModParam issueAMod )
{
 if( subset.empty() )
  return;

 if( v_inflows.empty() ) {
  if( std::all_of( values ,
                   values + subset.size() ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  v_inflows.resize( boost::extents[ f_NumberReservoirs ][ f_time_horizon ] );
 }

 // If nothing changes, return
 bool identical = true;
 for( auto i : subset ) {
  if( i >= v_inflows.size() )
   throw( std::invalid_argument( "HydroUnitBlock::set_inflow: "
                                 "invalid value in subset." ) );

  if( *( v_inflows.data() + i ) != *(values++) )
   identical = false;
 }

 if( identical )
  return;

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  for( auto i : subset ) {
   Index t = i % f_time_horizon;
   Index r = i / f_time_horizon;
   v_inflows[ r ][ t ] = *(values++);
  }

  if( constraints_generated() )
   // Change the abstract representation

   for( auto i : subset ) {
    Index t = i % f_time_horizon;
    Index r = i / f_time_horizon;

    if( t == 0 )
     FinalVolumeReservoir_Const[ t ][ r ].set_both(
      v_InitialVolumetric[ r ] + v_inflows[ r ][ t ] , issueAMod );
    else
     FinalVolumeReservoir_Const[ t ][ r ].set_both(
      v_inflows[ r ][ t ] , issueAMod );
   }
 }

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification
  if( ! ordered )
   std::sort( subset.begin() , subset.end() );

  Block::add_Modification( std::make_shared< HydroUnitBlockSbstMod >(
                            this , HydroUnitBlockMod::eSetInf , std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );
 }
}  // end( HydroUnitBlock::set_inflow( subset ) )

/*--------------------------------------------------------------------------*/

void HydroUnitBlock::set_inflow( MF_dbl_it values ,
                                 Block::Range rng ,
                                 c_ModParam issuePMod ,
                                 c_ModParam issueAMod )
{
 rng.second = std::min( rng.second ,
                        get_time_horizon() * get_number_reservoirs() );
 if( rng.second <= rng.first )
  return;

 if( v_inflows.empty() ) {
  if( std::all_of( values ,
                   values + ( rng.second - rng.first ) ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  v_inflows.resize( boost::extents[ f_NumberReservoirs ][ f_time_horizon ] );
 }

 // If nothing changes, return
 if( std::equal( values ,
                 values + ( rng.second - rng.first ) ,
                 v_inflows.data() + rng.first ) )
  return;

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation
  std::copy( values ,
             values + ( rng.second - rng.first ) ,
             v_inflows.data() + rng.first );

  if( constraints_generated() ) {
   // Change the abstract representation

   for( Index i = rng.first ; i < rng.second ; ++i ) {
    Index t = i % f_time_horizon;
    Index r = i / f_time_horizon;

    if( t == 0 )
     FinalVolumeReservoir_Const[ t ][ r ].set_both(
      v_InitialVolumetric[ r ] + v_inflows[ r ][ t ] , issueAMod );
    else
     FinalVolumeReservoir_Const[ t ][ r ].set_both(
      v_inflows[ r ][ t ] , issueAMod );
   }
  }
 }

 if( issue_pmod( issuePMod ) )
  Block::add_Modification( std::make_shared< HydroUnitBlockRngdMod >(
                            this , HydroUnitBlockMod::eSetInf , rng ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( HydroUnitBlock::set_inflow( range ) )

/*--------------------------------------------------------------------------*/

void HydroUnitBlock::set_inertia_power( MF_dbl_it values ,
                                        Subset && subset ,
                                        const bool ordered ,
                                        c_ModParam issuePMod ,
                                        c_ModParam issueAMod )
{
 if( subset.empty() )
  return;

 if( v_InertiaPower.empty() ) {
  if( std::all_of( values ,
                   values + subset.size() ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  v_InertiaPower.resize( boost::extents[ f_time_horizon ][ f_NumberArcs ] );
 }

 // If nothing changes, return
 bool identical = true;
 for( auto i : subset ) {
  if( i >= v_InertiaPower.size() )
   throw( std::invalid_argument( "HydroUnitBlock::set_inertia_power: "
                                 "invalid value in subset." ) );

  if( *( v_InertiaPower.data() + i ) != *(values++) )
   identical = false;
 }

 if( identical )
  return;

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  for( auto i : subset ) {
   Index a = i % f_NumberArcs;
   Index t = i / f_NumberArcs;
   v_InertiaPower[ t ][ a ] = *(values++);
  }

  if( constraints_generated() ) {
   // Change the abstract representation
   // FIXME: v_InertiaPower is not used
  }
 }

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification
  if( ! ordered )
   std::sort( subset.begin() , subset.end() );

  Block::add_Modification( std::make_shared< HydroUnitBlockSbstMod >(
                            this , HydroUnitBlockMod::eSetInerP , std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );
 }
}  // end( HydroUnitBlock::set_inertia_power( subset ) )

/*--------------------------------------------------------------------------*/

void HydroUnitBlock::set_inertia_power( MF_dbl_it values ,
                                        Block::Range rng ,
                                        c_ModParam issuePMod ,
                                        c_ModParam issueAMod )
{
 rng.second = std::min( rng.second , f_NumberArcs * get_time_horizon() );
 if( rng.second <= rng.first )
  return;

 if( v_InertiaPower.empty() ) {
  if( std::all_of( values ,
                   values + ( rng.second - rng.first ) ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  v_InertiaPower.resize( boost::extents[ f_time_horizon ][ f_NumberArcs ] );
 }

 // If nothing changes, return
 if( std::equal( values ,
                 values + ( rng.second - rng.first ) ,
                 v_InertiaPower.data() + rng.first ) )
  return;

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation
  std::copy( values ,
             values + ( rng.second - rng.first ) ,
             v_InertiaPower.data() + rng.first );

  if( constraints_generated() ) {
   // Change the abstract representation
   // FIXME: v_InertiaPower is not used
  }
 }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< HydroUnitBlockRngdMod >(
                            this , HydroUnitBlockMod::eSetInerP , rng ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( HydroUnitBlock::set_inertia_power( range ) )

/*--------------------------------------------------------------------------*/

void HydroUnitBlock::set_initial_volume( MF_dbl_it values ,
                                         Block::Subset && subset ,
                                         const bool ordered ,
                                         c_ModParam issuePMod ,
                                         c_ModParam issueAMod )
{
 if( subset.empty() )
  return;

 if( v_InitialVolumetric.empty() ) {
  // The initial volumes are currently zero.
  if( std::all_of( values ,
                   values + subset.size() ,
                   []( double cst ) { return( cst == 0 ); } ) )
   // The initial volumes are still zero. There is nothing to be updated.
   return;

  v_InitialVolumetric.resize( get_number_reservoirs() );
 }

 bool identical = true;
 for( auto r : subset ) {
  if( r >= v_InitialVolumetric.size() )
   throw( std::invalid_argument( "HydroUnitBlock::set_initial_volume: invalid "
                                 "index in subset: " + std::to_string( r ) ) );

  const auto volume = *(values++);
  if( v_InitialVolumetric[ r ] != volume ) {
   identical = false;

   if( not_dry_run( issuePMod ) ) {
    // Change the physical representation
    v_InitialVolumetric[ r ] = volume;
   }
  }
 }

 if( identical )
  // Nothing has changed.
  return;

 if( not_dry_run( issuePMod ) &&
     not_dry_run( issueAMod ) &&
     constraints_generated() ) {
  // Change the abstract representation
  for( auto r : subset ) {
   FinalVolumeReservoir_Const[ 0 ][ r ].set_both
    ( v_InitialVolumetric[ r ] + v_inflows[ r ][ 0 ] , issueAMod );
  }
 }

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification
  if( ! ordered )
   std::sort( subset.begin() , subset.end() );

  Block::add_Modification( std::make_shared< HydroUnitBlockSbstMod >(
                            this , HydroUnitBlockMod::eSetInitV , std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );
 }
}  // end( HydroUnitBlock::set_initial_volume( subset ) )

/*--------------------------------------------------------------------------*/

void HydroUnitBlock::set_initial_volume( MF_dbl_it values ,
                                         Block::Range rng ,
                                         c_ModParam issuePMod ,
                                         c_ModParam issueAMod )
{
 rng.second = std::min( rng.second , get_number_reservoirs() );
 if( rng.second <= rng.first )
  return;

 if( v_InitialVolumetric.empty() ) {
  // The initial volumes are currently zero.
  if( std::all_of( values ,
                   values + ( rng.second - rng.first ) ,
                   []( double cst ) { return( cst == 0 ); } ) )
   // The initial volumes are still zero. There is nothing to be updated.
   return;

  v_InitialVolumetric.resize( get_number_reservoirs() );
 }

 // If nothing changes, return
 if( std::equal( values , values + ( rng.second - rng.first ) ,
                 v_InitialVolumetric.begin() + rng.first ) )
  return;

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation
  std::copy( values ,
             values + ( rng.second - rng.first ) ,
             v_InitialVolumetric.begin() + rng.first );

  if( not_dry_run( issueAMod ) && constraints_generated() ) {
   // Change the abstract representation
   for( Index r = rng.first ; r < rng.second ; ++r ) {
    FinalVolumeReservoir_Const[ 0 ][ r ].set_both
     ( v_InitialVolumetric[ r ] + v_inflows[ r ][ 0 ] , issueAMod );
   }
  }
 }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< HydroUnitBlockRngdMod >(
                            this , HydroUnitBlockMod::eSetInitV , rng ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( HydroUnitBlock::set_initial_volume( range ) )

/*--------------------------------------------------------------------------*/

void HydroUnitBlock::update_initial_flow_rate_in_cnstrs( const Block::Subset & arcs ,
                                                         c_ModParam issueAMod )
{
 if( ! constraints_generated() )
  return;

 // ramp-up constraints
 if( ! ( RampUp_Const.empty() || v_DeltaRampUp.empty() ) ) {
  for( auto arc : arcs )
   RampUp_Const[ 0 ][ arc ].set_rhs( get_initial_flow_rate( arc ) +
                                     v_DeltaRampUp[ 0 ][ arc ] , issueAMod );
 }

 // ramp-down constraints
 if( ! ( RampDown_Const.empty() || v_DeltaRampDown.empty() ) ) {
  for( auto arc : arcs )
   RampDown_Const[ 0 ][ arc ].set_lhs( get_initial_flow_rate( arc ) -
                                       v_DeltaRampDown[ 0 ][ arc ] ,
                                       issueAMod );
 }
}  // end( HydroUnitBlock::update_initial_flow_rate_in_cnstrs )

/*--------------------------------------------------------------------------*/

void HydroUnitBlock::update_initial_flow_rate_in_cnstrs(Block::Range arcs ,
                                                        c_ModParam issueAMod )
{
 if( ! constraints_generated() )
  return;

 // ramp-up constraints
 if( ! ( RampUp_Const.empty() || v_DeltaRampUp.empty() ) ) {
  for( Index arc = arcs.first ; arc < arcs.second ; ++arc )
   RampUp_Const[ 0 ][ arc ].set_rhs( get_initial_flow_rate( arc ) +
                                     v_DeltaRampUp[ 0 ][ arc ] , issueAMod );
 }

 // ramp-down constraints
 if( ! ( RampDown_Const.empty() || v_DeltaRampDown.empty() ) ) {
  for( Index arc = arcs.first ; arc < arcs.second ; ++arc )
   RampDown_Const[ 0 ][ arc ].set_lhs
    ( get_initial_flow_rate( arc ) - v_DeltaRampDown[ 0 ][ arc ] ,
      issueAMod );
 }
}  // end( HydroUnitBlock::update_initial_flow_rate_in_cnstrs )

/*--------------------------------------------------------------------------*/

void HydroUnitBlock::set_initial_flow_rate( MF_dbl_it values ,
                                            Block::Subset && subset ,
                                            const bool ordered ,
                                            c_ModParam issuePMod ,
                                            c_ModParam issueAMod )
{
 if( subset.empty() )
  return;

 if( v_InitialFlowRate.empty() ) {
  if( std::all_of( values ,
                   values + subset.size() ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  auto max_index = *std::max_element( std::begin( subset ) ,
                                      std::end( subset ) );
  v_InitialFlowRate.resize( max_index );
 }

 bool identical = true;
 for( auto i : subset ) {
  if( i >= v_InitialFlowRate.size() )
   throw( std::invalid_argument( "HydroUnitBlock::set_initial_flow_rate: "
                                  "invalid value in subset." ) );
  auto flow_rate = *(values++);
  if( v_InitialFlowRate[ i ] != flow_rate ) {
   identical = false;
   if( not_dry_run( issuePMod ) )
    // Change the physical representation
    v_InitialFlowRate[ i ] = flow_rate;
  }
 }
 if( identical )
  return;  // nothing changes; return

 if( not_dry_run( issuePMod ) &&
     not_dry_run( issueAMod ) &&
     constraints_generated() )
  // Change the abstract representation
  update_initial_flow_rate_in_cnstrs( subset , issueAMod );

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification
  if( ! ordered )
   std::sort( subset.begin() , subset.end() );

  Block::add_Modification( std::make_shared< HydroUnitBlockSbstMod >(
                            this , HydroUnitBlockMod::eSetInitF , std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );
 }
}  // end( HydroUnitBlock::set_initial_flow_rate( subset ) )

/*--------------------------------------------------------------------------*/

void HydroUnitBlock::set_initial_flow_rate( MF_dbl_it values ,
                                            Block::Range rng ,
                                            c_ModParam issuePMod ,
                                            c_ModParam issueAMod )
{
 rng.second = std::min( rng.second , get_number_generators() );
 if( rng.second <= rng.first )
  return;

 if( v_InitialFlowRate.empty() ) {
  if( std::all_of( values ,
                   values + ( rng.second - rng.first ) ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  Index max_index = rng.second;
  v_InitialFlowRate.resize( max_index );
 }

  // If nothing changes, return
 else if( std::equal( values , values + ( rng.second - rng.first ) ,
                      v_InitialFlowRate.begin() + rng.first ) )
  return;

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation
  std::copy( values ,
             values + ( rng.second - rng.first ) ,
             v_InitialFlowRate.begin() + rng.first );

  if( not_dry_run( issueAMod ) && constraints_generated() )
   // Change the abstract representation
   update_initial_flow_rate_in_cnstrs( rng , issueAMod );
 }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< HydroUnitBlockRngdMod >(
                            this , HydroUnitBlockMod::eSetInitF , rng ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( HydroUnitBlock::set_initial_flow_rate( range ) )

/*--------------------------------------------------------------------------*/

template< typename T >
void HydroUnitBlock::transpose( boost::multi_array< T , 2 > & a )
{
 long rows = a.shape()[ 0 ];
 long cols = a.shape()[ 1 ];

 if( ( rows > 1 ) && ( cols == 1 ) && ( f_NumberArcs > 1 ) ) {
  // The given array has dimensions ( number of arcs x 1 ). Therefore, the
  // array must be transposed.
  boost::array< typename boost::multi_array< T , 2 >::index , 2 >
   dims = { { 1 , rows } };
  a.reshape( dims );
 }
}  // end( HydroUnitBlock::transpose )

/*--------------------------------------------------------------------------*/

void HydroUnitBlock::decompress_array( boost::multi_array< double , 2 > & array )
{
 // This function receives an array whose dimensions are N x f_NumberArcs,
 // where N can be 1, time horizon, or the number of change intervals. The
 // array can also be empty, in which case nothing is done.

 if( array.empty() )
  return;

 const auto num_rows = array.shape()[ 0 ];

 if( num_rows == 1 ) {
  // For each arc, the data is the same for every time instant. For arc r, the
  // data at time t is equal to given_array[ 0 ][ r ] for each t in {0, ...,
  // time_horizon - 1}. We resize the array so that its dimensions becomes
  // f_time_horizon x f_NumberArcs and copy the given data.
  boost::multi_array< double , 2 > given_array = array;
  array.resize( boost::extents[ f_time_horizon ][ f_NumberArcs ] );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   for( Index r = 0 ; r < f_NumberArcs ; ++r )
    array[ t ][ r ] = given_array[ 0 ][ r ];
 } else if( num_rows < f_time_horizon ) {
  // Since the number of rows is greater than 1 and less than the time
  // horizon, it must be equal to the number of change intervals.
  if( num_rows != v_change_intervals.size() )
   throw( std::logic_error(
    "HydroUnitBlock::decompress_array: invalid number of rows (" +
    std::to_string( num_rows ) + ") for some variable. It should " +
    "be equal to the number of change intervals (" +
    std::to_string( v_change_intervals.size() ) + ")" ) );

  // For time instant t and arc r, the value for arc r at time t is equal to
  // given_array[ k ][ r ], where k is such that t belongs to the closed
  // interval [i_{k-1} + 1, i_k] and i_k is the k-th element of
  // v_change_intervals (starting from k = 0) and i_{-1} = -1 by
  // definition. We resize the array so that its dimensions becomes
  // f_time_horizon x f_NumberArcs and copy the given data.

  boost::multi_array< double , 2 > given_array = array;
  array.resize( boost::extents[ f_time_horizon ][ f_NumberArcs ] );
  for( Index r = 0 ; r < f_NumberArcs ; ++r ) {
   Index t = 0;
   for( Index k = 0 ; k < v_change_intervals.size() ; ++k ) {
    auto upper_endpoint = v_change_intervals[ k ];
    if( k == v_change_intervals.size() - 1 )
     // The upper endpoint of the last interval must be time_horizon -
     // 1. Since it may not be provided in v_change_intervals (the value for
     // the last element of v_change_intervals is not required), we manually
     // set it here.
     upper_endpoint = f_time_horizon - 1;
    for( ; t <= upper_endpoint ; ++t )
     array[ t ][ r ] = given_array[ k ][ r ];
   }
  }
 }
}  // end( HydroUnitBlock::decompress_array )

/*--------------------------------------------------------------------------*/

void HydroUnitBlock::decompress_vol( boost::multi_array< double , 2 > & array )
{
 // This function receives an array whose dimensions are f_NumberReservoirs x
 // N, where N can be 1, time horizon, or the number of change intervals. The
 // array can also be empty, in which case nothing is done.

 if( array.empty() )
  return;

 const auto num_columns = array.shape()[ 1 ];

 if( num_columns == 1 ) {
  // The maximum or minimum volume is the same for every time instant. For
  // reservatory r, the maximum or minimum volume at time t is equal to
  // given_array[ r ][ 0 ] for each t in {0, ..., time_horizon - 1}. We resize
  // the array so that its dimensions becomes f_NumberReservoirs x
  // f_time_horizon and copy the given data.
  boost::multi_array< double , 2 > given_array = array;
  array.resize( boost::extents[ f_NumberReservoirs ][ f_time_horizon ] );
  for( Index r = 0 ; r < f_NumberReservoirs ; ++r )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    array[ r ][ t ] = given_array[ r ][ 0 ];
 } else if( num_columns < f_time_horizon ) {
  // Since the number of columns is greater than 1 and less than the time
  // horizon, it must be equal to the number of change intervals.
  if( num_columns != v_change_intervals.size() )
   throw( std::logic_error(
    "HydroUnitBlock::decompress_vol: invalid number of columns (" +
    std::to_string( num_columns ) + ") for the maximum or minimum " +
    "volume variable. It should be equal to the number of change intervals"
    " (" + std::to_string( v_change_intervals.size() ) + ")" ) );

  // For each reservatory r and time instant t, the maximum or minimum volume
  // of reservatory r at time t is equal to given_array[ r ][ k ], where k is
  // such that t belongs to the closed interval [i_{k-1} + 1, i_k] and i_k is
  // the k-th element of v_change_intervals (starting from k = 0) and i_{-1} =
  // -1 by definition. We resize the array so that its dimensions becomes
  // f_NumberReservoirs x f_time_horizon and copy the given data.

  boost::multi_array< double , 2 > given_array = array;
  array.resize( boost::extents[ f_NumberReservoirs ][ f_time_horizon ] );
  for( Index r = 0 ; r < f_NumberReservoirs ; ++r ) {
   Index t = 0;
   for( Index k = 0 ; k < v_change_intervals.size() ; ++k ) {
    auto upper_endpoint = v_change_intervals[ k ];
    if( k == v_change_intervals.size() - 1 )
     // The upper endpoint of the last interval must be time_horizon -
     // 1. Since it may not be provided in v_change_intervals (the value for
     // the last element of v_change_intervals is not required), we manually
     // set it here.
     upper_endpoint = f_time_horizon - 1;
    for( ; t <= upper_endpoint ; ++t )
     array[ r ][ t ] = given_array[ r ][ k ];
   }
  }
 }
}  // end( HydroUnitBlock::decompress_vol )

/*--------------------------------------------------------------------------*/
/*------------------- End File HydroUnitBlock.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
