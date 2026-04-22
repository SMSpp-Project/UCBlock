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

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register HydroUnitBlock to the Block factory
SMSpp_insert_in_factory_cpp_1( HydroUnitBlock );

// register HydroUnitBlockSolution to the Solution factory
SMSpp_insert_in_factory_cpp_0( HydroUnitBlockSolution );

/*--------------------------------------------------------------------------*/
/*--------------------------- STATIC FUNCTIONS -----------------------------*/
/*--------------------------------------------------------------------------*/

template< class T , std::size_t K >
static inline void copy_multi_array( boost::multi_array< T , K > & to ,
			       const boost::multi_array< T , K > & from )
{
 std::vector< size_t > extent;
 auto shape = from.shape();
 extent.assign( shape , shape + from.num_dimensions() );
 to.resize( extent );
 to = from;
 }

/*--------------------------------------------------------------------------*/
// should have worked with begin() and end(), but it seems boost::multi_array
// dramatically flunks these definitions somehow

template< class T , std::size_t K >
static inline bool all0( boost::multi_array< T , K > & v )
{
 return( std::all_of( v.data() , v.data() + v.num_elements() ,
		      []( T i ) { return( i == 0 ); } ) );
 }

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
 // Deserialize data from the base class
 UnitBlock::deserialize( group );

 if( ! deserialize_dim( group , "NumberReservoirs" , f_NumberReservoirs ) )
  f_NumberReservoirs = 1;

 if( ! deserialize_dim( group , "NumberArcs" , f_NumberArcs ) )
  f_NumberArcs = 1;

 ::deserialize( group , "NumberPieces" , f_NumberArcs ,
                v_NumberPieces , true , true );

 if( ! deserialize_dim( group , "TotalNumberPieces" , f_TotalNumberPieces )
     ) {
  f_TotalNumberPieces = 0;
  for( const auto & n : v_NumberPieces )
   f_TotalNumberPieces += n;
  }

 f_TotalNumberPieces = f_TotalNumberPieces ?
                       f_TotalNumberPieces : f_NumberArcs;

 ::deserialize( group , "StartArc" , f_NumberArcs , v_StartArc );
 ::deserialize( group , "EndArc" , f_NumberArcs , v_EndArc );

 ::deserialize( group , "Inflows" , { f_NumberReservoirs , f_time_horizon } ,
                v_inflows , true , false , v_change_intervals );

 ::deserialize( group , "MinFlow" , { f_time_horizon , f_NumberArcs } ,
                v_MinFlow , true , true , v_change_intervals );

 ::deserialize( group , "MaxFlow" , { f_time_horizon , f_NumberArcs } ,
                v_MaxFlow , true , true , v_change_intervals );

 ::deserialize( group , "MinPower" , { f_time_horizon , f_NumberArcs } ,
                v_MinPower , true , true , v_change_intervals );

 ::deserialize( group , "MaxPower" , { f_time_horizon , f_NumberArcs } ,
                v_MaxPower , true , true , v_change_intervals );

 ::deserialize( group , "DeltaRampUp" , { f_time_horizon , f_NumberArcs } ,
                v_DeltaRampUp , true , true , v_change_intervals );

 ::deserialize( group , "DeltaRampDown" , { f_time_horizon , f_NumberArcs } ,
                v_DeltaRampDown , true , true , v_change_intervals );

 ::deserialize( group , "PrimaryRho" , { f_time_horizon , f_NumberArcs } ,
                v_PrimaryRho , true , true , v_change_intervals );

 ::deserialize( group , "SecondaryRho" , { f_time_horizon , f_NumberArcs } ,
                v_SecondaryRho , true , true , v_change_intervals );

 ::deserialize( group , "LinearTerm" , f_TotalNumberPieces ,
                v_LinearTerm , true , true );

 ::deserialize( group , "ConstantTerm" , f_TotalNumberPieces ,
                v_ConstTerm , true , true );

 ::deserialize( group , "ActivePowerCost" , f_NumberArcs ,
                v_ActivePowerCost , true , true );

 ::deserialize( group , "InertiaPower" , { f_NumberArcs , f_time_horizon } ,
                v_InertiaPower , true , true , v_change_intervals );

 ::deserialize( group , "InitialFlowRate" , f_NumberArcs ,
                v_InitialFlowRate , true , true );

 ::deserialize( group , "InitialVolumetric" , f_NumberReservoirs ,
                v_InitialVolumetric , true , true );

 ::deserialize( group , "UphillFlow" , f_NumberArcs ,
                v_UphillDelay , true , true );

 ::deserialize( group , "DownhillFlow" , f_NumberArcs ,
                v_DownhillDelay , true , true );

 ::deserialize( group , "MinVolumetric" ,
                { f_NumberReservoirs , f_time_horizon } , v_MinVolumetric ,
                true , true , v_change_intervals );

 ::deserialize( group , "MaxVolumetric" ,
                { f_NumberReservoirs , f_time_horizon } , v_MaxVolumetric ,
                true , true , v_change_intervals );

 /// optional AC variables
 if( ::deserialize( group , "MinReactivePower" ,
		    { f_time_horizon , f_NumberArcs } ,
		    v_MinReactivePower , true , true , v_change_intervals ) )
  if( all0( v_MinReactivePower ) )
   v_MinReactivePower.resize( MAdouble_ext()[ 0 ][ 0 ] );

 if( ::deserialize( group , "MaxReactivePower" ,
		    { f_time_horizon , f_NumberArcs } ,
		    v_MaxReactivePower , true , true , v_change_intervals ) )
  if( all0( v_MaxReactivePower ) )
   v_MaxReactivePower.resize( MAdouble_ext()[ 0 ][ 0 ] );
 
 ::deserialize( group , "ReferenceSchedule" , f_time_horizon , v_RefSchedule ,
                true , true , v_change_intervals );

 }  // end( HydroUnitBlock::deserialize )

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG

std::vector< std::string > HydroUnitBlock::expected_dims( void ) const {
 static const std::vector< std::string > ed =
 { "NumberReservoirs" , "NumberArcs" , "TotalNumberPieces" };

 auto ret = UnitBlock::expected_dims();
 ret.insert( ret.end() , ed.begin() , ed.end() );

 return( ret );
 }

/*--------------------------------------------------------------------------*/

std::vector< std::string > HydroUnitBlock::expected_vars( void ) const {
 static const std::vector< std::string > ev =
 { "StartArc" , "EndArc" , "MinFlow" , "MaxFlow" , "MinVolumetric" ,
   "MaxVolumetric" , "Inflows" , "MinPower" , "MaxPower" , "DeltaRampUp" ,
   "DeltaRampDown" , "PrimaryRho" , "SecondaryRho" , "NumberPieces" ,
   "LinearTerm" , "ConstantTerm" , "ActivePowerCost" , "InertiaPower" ,
   "InitialFlowRate" , "InitialVolumetric" , "UphillFlow" , "DownhillFlow" ,
   "MinReactivePower", "MaxReactivePower", "ReferenceSchedule"
   };

 auto ret = UnitBlock::expected_vars();
 ret.insert( ret.end() , ev.begin() , ev.end() );

 return( ret );
 }

#endif

/*--------------------------------------------------------------------------*/

void HydroUnitBlock::generate_abstract_variables( Configuration * stvv )
{
 if( variables_generated() )  // variables have already been generated
  return;                     // nothing to do

 UnitBlock::generate_abstract_variables( stvv );

 if( f_time_horizon == 0 )  // there are no variables to be generated
  return;

 // volumetric Variable - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 v_volumetric.resize( boost::extents[ f_NumberReservoirs ][ f_time_horizon ]
		      );
 for( Index g = 0 ; g < f_NumberReservoirs ; ++g )
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   v_volumetric[ g ][ t ].set_type( ColVariable::kNonNegative );
 add_static_variable( v_volumetric , "v_hydro" );

 // flow and Active Power Variable- - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 v_flow_rate.resize( boost::extents[ f_NumberArcs ][ f_time_horizon ] );
 v_active_power.resize( boost::extents[ f_NumberArcs ][ f_time_horizon ] );

 for( Index g = 0 ; g < f_NumberArcs ; ++g )
  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   v_flow_rate[ g ][ t ].set_type( ColVariable::kContinuous , eNoMod );
   // put sign constraints on the flow variables in accordance to MinFlow
   // and MaxFlow (basically, turbines have >= 0 flows and pumps have <=
   // flows, although note that the funny case in which both happens and
   // the flow is fixed to 0 is not excluded)
   if( get_min_flow( t , g ) >= 0 )
    v_flow_rate[ g ][ t ].is_positive( true , eNoMod );
   if( get_max_flow( t , g ) <= 0 )
    v_flow_rate[ g ][ t ].is_negative( true , eNoMod );

   v_active_power[ g ][ t ].set_type( ColVariable::kContinuous , eNoMod );
   // put sign constraints on the active power variables in accordance to
   // MinPower and MaxPower
   if( get_min_power( t , g ) >= 0 )
    v_active_power[ g ][ t ].is_positive( true , eNoMod );
   if( get_max_power( t , g ) <= 0 )
    v_active_power[ g ][ t ].is_negative( true , eNoMod );
   }

 add_static_variable( v_flow_rate , "f_hydro" );
 add_static_variable( v_active_power , "p_hydro" );

 // Reactive Power Variables, if any- - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( f_reactive_power ) {
  v_reactive_power.resize( boost::extents[ f_NumberArcs ][ f_time_horizon ] );

  // put sign constraints on the reactive power variables in accordance to
  // MinReactivePower and MaxReactivePower
  for( Index g = 0 ; g < f_NumberArcs ; ++g )
   for( Index t = 0 ; t < f_time_horizon ; ++t ) {
    v_reactive_power[ g ][ t ].set_type( ColVariable::kContinuous , eNoMod );
    if( get_min_reactive_power( t , g ) >= 0 )
     v_reactive_power[ g ][ t ].is_positive( true , eNoMod );
    if( get_max_reactive_power( t , g ) <= 0 )
     v_reactive_power[ g ][ t ].is_negative( true , eNoMod ); 
    }

  add_static_variable( v_reactive_power , "q_hydro" );
  }

 // reserve Variable- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( reserve_vars & 1u )          // if UCBlock has primary demand variables
  if( ! v_PrimaryRho.empty() ) {  // if unit produces any primary reserve
   v_primary_spinning_reserve.resize(
                          boost::extents[ f_NumberArcs ][ f_time_horizon ] );
   for( Index g = 0 ; g < f_NumberArcs ; ++g )
    for( Index t = 0 ; t < f_time_horizon ; ++t )
     v_primary_spinning_reserve[ g ][ t ].set_type(
						 ColVariable::kNonNegative );

   add_static_variable( v_primary_spinning_reserve , "pr_hydro" );
   }

 if( reserve_vars & 2u )         // if UCBlock has secondary demand variables
  if( ! v_SecondaryRho.empty() ) {  // if unit produces any secondary reserve
   v_secondary_spinning_reserve.resize(
                         boost::extents[ f_NumberArcs ][ f_time_horizon ] );
   for( Index g = 0 ; g < f_NumberArcs ; ++g )
    for( Index t = 0 ; t < f_time_horizon ; ++t )
     v_secondary_spinning_reserve[ g ][ t ].set_type(
						ColVariable::kNonNegative );

   add_static_variable( v_secondary_spinning_reserve , "sr_hydro" );
   }

 // reference schedule Variable, if there - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! v_RefSchedule.empty() ) {
  v_abs_ref_schedule.resize( f_time_horizon );
  for( auto & var : v_abs_ref_schedule )
   var.set_type( ColVariable::kNonNegative );
  add_static_variable( v_abs_ref_schedule , "v_absh_refschd" );
  }

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 set_variables_generated();

 }  // end( HydroUnitBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void HydroUnitBlock::generate_abstract_constraints( Configuration * stcc )
{
 if( constraints_generated() )  // constraints have already been generated
  return;                       // nothing to do

 LinearFunction::v_coeff_pair vars;

 using maFRC2 = boost::multi_array< FRowConstraint , 2 >;

 // final volume constraints for each reservoir - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 assert( FinalVolumeReservoir_Const.empty() );

 FinalVolumeReservoir_Const.resize(
             maFRC2::extent_gen()[ f_time_horizon ][ f_NumberReservoirs ] );

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
     }
    else
     vars.push_back( std::make_pair( get_flow_rate( l , t ) , 1.0 ) );
    }

   vars.push_back( std::make_pair( get_volume( n , t ) , 1.0 ) );

   if( t > 0 )
    vars.push_back( std::make_pair( get_volume( n , t - 1 ) , -1.0 ) );

   double initial_volume = 0.0;
   if( t == 0 ) {
    if( v_InitialVolumetric[ n ] >= 0. )
     initial_volume = v_InitialVolumetric[ n ];
    else
     vars.push_back( std::make_pair( get_volume( n , f_time_horizon - 1 ) ,
                                     -1.0 ) );
    }

   if( ! v_inflows.empty() )
    FinalVolumeReservoir_Const[ t ][ n ].set_both( initial_volume +
						   v_inflows[ n ][ t ] );
   else
    FinalVolumeReservoir_Const[ t ][ n ].set_both( initial_volume );

   FinalVolumeReservoir_Const[ t ][ n ].set_function(
				  new LinearFunction( std::move( vars ) ) );
   }
  }

 add_static_constraint( FinalVolumeReservoir_Const ,
                        "FinalVolumeReservoir_HydroUnit" );

 // maximum power output according to primary-secondary reserves constraints
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // Initial data check
 if( ( ! v_MinPower.empty() ) && ( ! v_MaxPower.empty() ) )
  for( Index arc = 0 ; arc < f_NumberArcs ; ++arc )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    if( v_MinPower[ t ][ arc ] > v_MaxPower[ t ][ arc ] )
     throw( std::logic_error( "HydroUnitBlock::maximum and minimum power "
                              "output constraints: it must be that v_MaxPower"
                              " >= v_MinPower" ) );

 // Initial data check
 if( ( ! v_MinFlow.empty() ) && ( ! v_MaxFlow.empty() ) )
  for( Index arc = 0 ; arc < f_NumberArcs ; ++arc )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    if( ( v_MaxFlow[ t ][ arc ] <= 0 ) &&
        ( v_MinFlow[ t ][ arc ] < 0 ) &&
        ( v_NumberPieces[ arc ] > 1 ) )
     throw( std::logic_error( "HydroUnitBlock::Data Error: it must be that "
                              "for each pump when v_MinPower < 0, then "
                              "v_NumberPieces == 1" ) );

 // Initial data check
 if( ( ! v_MinFlow.empty() ) && ( ! v_MaxFlow.empty() ) )
  for( Index arc = 0 ; arc < f_NumberArcs ; ++arc )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    if( ! v_PrimaryRho.empty() )
     if( ( v_MaxFlow[ t ][ arc ] <= 0 ) &&
         ( v_MinFlow[ t ][ arc ] < 0 ) &&
         ( v_PrimaryRho[ t ][ arc ] != 0 ) )
      throw( std::logic_error( "HydroUnitBlock::Data Error: it must be that "
                               "for each pump v_PrimaryRho == 0" ) );

 // Initial data check
 if( ( ! v_MinFlow.empty() ) && ( ! v_MaxFlow.empty() ) )
  for( Index arc = 0 ; arc < f_NumberArcs ; ++arc )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    if( ! v_SecondaryRho.empty() )
     if( ( v_MaxFlow[ t ][ arc ] <= 0 ) &&
         ( v_MinFlow[ t ][ arc ] < 0 ) &&
         ( v_SecondaryRho[ t ][ arc ] != 0 ) )
      throw( std::logic_error( "HydroUnitBlock::Data Error: it must be that "
                               "for each pump then v_SecondaryRho == 0" ) );

 assert( MaxPowerPrimarySecondary_Const.empty() );

 MaxPowerPrimarySecondary_Const.resize(
		   maFRC2::extent_gen()[ f_time_horizon ][ f_NumberArcs ] );

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
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 assert( MinPowerPrimarySecondary_Const.empty() );

 MinPowerPrimarySecondary_Const.resize(
		   maFRC2::extent_gen()[ f_time_horizon ][ f_NumberArcs ] );

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

 // power output relation with primary reserves constraints - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( reserve_vars & 1u ) {  // if UCBlock has primary demand variables
  if( ! v_PrimaryRho.empty() ) {  // if unit produces any primary reserve
   assert( ActivePowerPrimary_Const.empty() );
   ActivePowerPrimary_Const.resize(
                    maFRC2::extent_gen()[ f_time_horizon ][ f_NumberArcs ] );

   for( Index arc = 0 ; arc < f_NumberArcs ; ++arc ) {
    for( Index t = 0 ; t < f_time_horizon ; ++t ) {
     auto MinF = get_MinFlow( t , arc );
     auto MaxF = get_MaxFlow( t , arc );

     if( ( MinF >= 0 ) && ( MaxF > 0 ) ) {  // Turbines
      vars.push_back( std::make_pair( get_active_power( arc , t ) ,
                                      v_PrimaryRho[ t ][ arc ] ) );
      vars.push_back( std::make_pair( get_primary_spinning_reserve( arc , t ) ,
                                      -1.0 ) );

      ActivePowerPrimary_Const[ t ][ arc ].set_lhs( 0.0 );
      ActivePowerPrimary_Const[ t ][ arc ].set_rhs( Inf< double >() );
      }
     else
      if( ( MaxF <= 0 ) && ( MinF < 0 ) ) {  // Pumps: pr_hydro == 0
       vars.push_back( std::make_pair( get_primary_spinning_reserve( arc , t ) ,
                                       1.0 ) );
       ActivePowerPrimary_Const[ t ][ arc ].set_both( 0.0 );
       }
      else {  // MinF == MaxF == 0: f_hydro == pr_hydro [== 0]
       vars.push_back( std::make_pair( get_flow_rate( arc , t ) , 1.0 ) );
       vars.push_back( std::make_pair( get_primary_spinning_reserve( arc , t ) ,
                                       -1.0 ) );
       ActivePowerPrimary_Const[ t ][ arc ].set_both( 0.0 );
       }

     ActivePowerPrimary_Const[ t ][ arc ].set_function(
				  new LinearFunction( std::move( vars ) ) );
     }
    }
   }

  add_static_constraint( ActivePowerPrimary_Const ,
                         "ActivePowerPrimary_HydroUnit" );
  }

 // power output relation with secondary reserves constraints - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( reserve_vars & 2u ) {  // if UCBlock has secondary demand variables
  if( ! v_SecondaryRho.empty() ) {  // if unit produces any secondary reserve
   assert( ActivePowerSecondary_Const.empty() );
   ActivePowerSecondary_Const.resize(
		   maFRC2::extent_gen()[ f_time_horizon ][ f_NumberArcs ] );

   for( Index arc = 0 ; arc < f_NumberArcs ; ++arc ) {
    for( Index t = 0 ; t < f_time_horizon ; ++t ) {
     auto MinF = get_MinFlow( t , arc );
     auto MaxF = get_MaxFlow( t , arc );

     if( ( MinF >= 0 ) && ( MaxF > 0 ) ) {  // Turbines
      vars.push_back( std::make_pair( get_active_power( arc , t ) ,
                                      v_SecondaryRho[ t ][ arc ] ) );
      vars.push_back( std::make_pair( get_secondary_spinning_reserve( arc , t ) ,
                                      -1.0 ) );

      ActivePowerSecondary_Const[ t ][ arc ].set_lhs( 0.0 );
      ActivePowerSecondary_Const[ t ][ arc ].set_rhs( Inf< double >() );
      }
     else
      if( ( MaxF <= 0 ) && ( MinF < 0 ) ) {  // Pumps: sr_hydro == 0
       vars.push_back( std::make_pair(
		       get_secondary_spinning_reserve( arc , t ) , 1.0 ) );

       ActivePowerSecondary_Const[ t ][ arc ].set_both( 0.0 );
       }
      else {  // MinF == MaxF == 0: f_hydro == sr_hydro [== 0]
       vars.push_back( std::make_pair( get_flow_rate( arc , t ) , 1.0 ) );
       vars.push_back( std::make_pair( get_secondary_spinning_reserve( arc , t ) ,
                                       -1.0 ) );

       ActivePowerSecondary_Const[ t ][ arc ].set_both( 0.0 );
       }

     ActivePowerSecondary_Const[ t ][ arc ].set_function(
				 new LinearFunction( std::move( vars ) ) );
     }
    }

   add_static_constraint( ActivePowerSecondary_Const ,
                          "ActivePowerSecondary_HydroUnit" );
   }
  }

 // flow active power constraints - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 assert( FlowActivePower_Const.empty() );
 FlowActivePower_Const.resize(
            maFRC2::extent_gen()[ f_time_horizon ][ f_TotalNumberPieces ] );

 if( f_NumberArcs > 0 ) {
  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   Index piece = 0;
   Index cnstr_idx = 0;
   Index end = 0;

   for( Index arc = 0 ; arc < f_NumberArcs ; ++arc ) {
    if( v_NumberPieces.empty() )
     ++end;
    else
     end += v_NumberPieces[ arc ];

    auto MinF = get_MinFlow( t , arc );
    auto MaxF = get_MaxFlow( t , arc );

    if( ( MinF >= 0 ) && ( MaxF > 0 ) ) {  // Turbines
     for( ; piece < end ; ++piece , ++cnstr_idx ) {
      vars.push_back( std::make_pair( get_active_power( arc , t ) , 1.0 ) );
      if( ! v_LinearTerm.empty() )
       vars.push_back( std::make_pair( get_flow_rate( arc , t ) ,
                                       -v_LinearTerm[ piece ] ) );
      else
       vars.push_back( std::make_pair( get_flow_rate( arc , t ) , 1.0 ) );

      if( ! v_ConstTerm.empty() )
       FlowActivePower_Const[ t ][ piece ].set_rhs( v_ConstTerm[ piece ] );
      else
       FlowActivePower_Const[ t ][ piece ].set_rhs( 0.0 );
      FlowActivePower_Const[ t ][ piece ].set_lhs( -Inf< double >() );
      FlowActivePower_Const[ t ][ piece ].set_function(
				 new LinearFunction( std::move( vars ) ) );
      }
     continue;
     }

    vars.push_back( std::make_pair( get_active_power( arc , t ) , 1.0 ) );

    if( ( MaxF <= 0 ) && ( MinF < 0 ) )  // Pumps
     vars.push_back( std::make_pair( get_flow_rate( arc , t ) ,
                                     -v_LinearTerm[ cnstr_idx ] ) );
    else  // MinF == MaxF == 0:  f_hydro == p_hydro [== 0]
     vars.push_back( std::make_pair( get_flow_rate( arc , t ) , -1.0 ) );

    FlowActivePower_Const[ t ][ cnstr_idx ].set_both( 0.0 );
    FlowActivePower_Const[ t ][ cnstr_idx++ ].set_function(
				 new LinearFunction( std::move( vars ) ) );
    piece = cnstr_idx;
    }
   }

  add_static_constraint( FlowActivePower_Const , "FlowActivePower_HydroUnit" );

  }  // end( if( f_NumberArcs > 0 ) )

 // flow rate bounds- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // Initial data check
 if( ( ! v_MinFlow.empty() ) && ( ! v_MaxFlow.empty() ) )
  for( Index arc = 0 ; arc < f_NumberArcs ; ++arc )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    if( v_MinFlow[ t ][ arc ] > v_MaxFlow[ t ][ arc ] )
     throw( std::logic_error( "HydroUnitBlock::flow rate variable bounds: "
			      "it must be that v_MaxFlow >= v_MinFlow." ) );

 assert( FlowRateBounds_Const.empty() );
 FlowRateBounds_Const.resize(
		   maFRC2::extent_gen()[ f_time_horizon ][ f_NumberArcs ] );

 for( Index t = 0 ; t < f_time_horizon ; ++t )
  for( Index arc = 0 ; arc < f_NumberArcs ; ++arc ) {
   auto flow_rate = get_flow_rate( arc , t );
   FlowRateBounds_Const[ t ][ arc ].set_lhs( get_MinFlow( t , arc ) );
   FlowRateBounds_Const[ t ][ arc ].set_rhs( get_MaxFlow( t , arc ) );
   FlowRateBounds_Const[ t ][ arc ].set_variable( flow_rate );
   }

 add_static_constraint( FlowRateBounds_Const , "FlowRateBounds_HydroUnit" );

 // ramp-up constraints - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! v_DeltaRampUp.empty() ) {
  assert( RampUp_Const.empty() );
  RampUp_Const.resize(
		  maFRC2::extent_gen()[ f_time_horizon ][ f_NumberArcs ] );

  // Initial condition
  for( Index arc = 0 ; arc < f_NumberArcs ; ++arc ) {
   vars.push_back( std::make_pair( get_flow_rate( arc , 0 ) , 1.0 ) );

   RampUp_Const[ 0 ][ arc ].set_lhs( -Inf< double >() );
   RampUp_Const[ 0 ][ arc ].set_rhs( v_DeltaRampUp[ 0 ][ arc ] +
                                     get_initial_flow_rate( arc ) );
   RampUp_Const[ 0 ][ arc ].set_function(
    new LinearFunction( std::move( vars ) ) );

   for( Index t = 1 , cnstr_idx = 1 ; t < f_time_horizon ;
	++t , ++cnstr_idx ) {
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

 // ramp-down constraints - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! v_DeltaRampDown.empty() ) {
  assert( RampDown_Const.empty() );
  RampDown_Const.resize(
		   maFRC2::extent_gen()[ f_time_horizon ][ f_NumberArcs ] );

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
    RampDown_Const[ cnstr_idx ][ arc ].set_rhs( v_DeltaRampDown[ t ][ arc ] );
    RampDown_Const[ cnstr_idx ][ arc ].set_function(
				    new LinearFunction( std::move( vars ) ) );
    }
   }

  add_static_constraint( RampDown_Const , "RampDown_HydroUnit" );
  }

 // volumetric bounds constraints - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // Initial data check
 if( ( ! v_MinVolumetric.empty() ) && ( ! v_MaxVolumetric.empty() ) )
  for( Index node = 0 ; node < f_NumberReservoirs ; ++node )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    if( ( v_MinVolumetric[ node ][ t ] > v_MaxVolumetric[ node ][ t ] ) ||
        ( v_MinVolumetric[ node ][ t ] < 0 ) ||
	( v_MaxVolumetric[ node ][ t ] < 0 ) )
     throw( std::logic_error( "HydroUnitBlock::Volumetric Bounds Constraint "
                              "must be 0 <= MinV[ r , t ] <= MaxV[ r , t ]"
			      ) );

 assert( VolumetricBounds_Const.empty() );
 VolumetricBounds_Const.resize(
	      maFRC2::extent_gen()[ f_NumberReservoirs ][ f_time_horizon ] );

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

 add_static_constraint( VolumetricBounds_Const ,
			"VolumetricBounds_HydroUnit" );

 // reference schedule- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! v_RefSchedule.empty() ) {
   Reference_Schedule_Const.resize( 2 * f_time_horizon );
   for( Index t = 0 ; t < f_time_horizon ; ++t ) {
        // | Sum P - Pref | <= v_abs_ref_schedule
        auto lfunc_1 = new LinearFunction();
        for( Index g = 0 ; g < f_NumberArcs ; ++g ) {
          lfunc_1->add_variable( & v_active_power[ g ][ t ], 1.0 );
        }
        lfunc_1->add_variable( & v_abs_ref_schedule[ t ], -1.0 );
        Reference_Schedule_Const[ t ].set_lhs( -Inf< double >() );
        Reference_Schedule_Const[ t ].set_rhs( v_RefSchedule[ t ] );
        Reference_Schedule_Const[ t ].set_function( lfunc_1 );
        //
        auto lfunc_2 = new LinearFunction();
        for( Index g = 0 ; g < f_NumberArcs ; ++g ) {
          lfunc_2->add_variable( & v_active_power[ g ][ t ], -1.0 );
        }
        lfunc_2->add_variable( & v_abs_ref_schedule[ t ], -1.0 );
        Reference_Schedule_Const[ f_time_horizon + t ].set_lhs( -Inf< double >() );
        Reference_Schedule_Const[ f_time_horizon + t ].set_rhs( -v_RefSchedule[ t ] );
        Reference_Schedule_Const[ f_time_horizon + t ].set_function( lfunc_2 );
   }
   add_static_constraint( Reference_Schedule_Const,
			  "Norm1_H_Reference_Schedule" );
 }

 // bounds on the reactive part - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( f_reactive_power && 
     ( ( ! v_MinReactivePower.empty() ) || ( ! v_MaxReactivePower.empty() ) )
     ) {
  if( ReactivePower_Bound_Const.empty() )
   ReactivePower_Bound_Const.resize(
                          boost::extents[ f_NumberArcs ][ f_time_horizon ] );

 for( Index g = 0 ; g < f_NumberArcs ; ++g )
  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   ReactivePower_Bound_Const[ g ][ t ].set_rhs(
					    get_max_reactive_power( t , g ) );
   ReactivePower_Bound_Const[ g ][ t ].set_lhs(
					    get_min_reactive_power( t , g ) );
   ReactivePower_Bound_Const[ g ][ t ].set_variable(
					       & v_reactive_power[ g ][ t ] );
   }

 add_static_constraint( ReactivePower_Bound_Const , "ReactivePowerBound" );

 /*!! Link between active and reactive power
 Reactive_2_Active_Const.resize(
                          boost::extents[ f_NumberArcs ][ f_time_horizon ] );

 for( Index g = 0 ; g < f_NumberArcs ; ++g )
  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   // Q(g,t) - P(g,t) <= 0
   auto lfunc = new LinearFunction();
   lfunc->add_variable( &v_active_power[ g ][ t ], -1.0 );
   lfunc->add_variable( &v_reactive_power[ g ][ t ], 1.0 );

   Reactive_2_Active_Const[ g ][ t ].set_lhs( -Inf< double >() );
   Reactive_2_Active_Const[ g ][ t ].set_rhs( 0.0 );
   Reactive_2_Active_Const[ g ][ t ].set_function( lfunc );
   }

 add_static_constraint( Reactive_2_Active_Const, "QandPhydro" );
 !!*/
 }

 // all done- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 set_constraints_generated();

 }  // end( HydroUnitBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void HydroUnitBlock::generate_objective( Configuration * objc )
{
 if( objective_generated() )  // Objective has already been generated
  return;                     // nothing to do

 // the variables to fill in - only when the reference schedule is there 
  LinearFunction::v_coeff_pair vars;

 if( ! v_ActivePowerCost.empty() )
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   for( Index arc = 0 ; arc < f_NumberArcs ; ++arc )
    vars.push_back( std::make_pair( get_active_power( arc , t ) ,
                                    v_ActivePowerCost[ arc ] ) );

 if( ! v_RefSchedule.empty() ) {
  for( Index t = 0 ; t < f_time_horizon ; ++t )
      vars.push_back( std::make_pair( &v_abs_ref_schedule[ t ] , 1.0 ) );
 }

 objective.set_function( new LinearFunction( std::move( vars ) ) );
 objective.set_sense( Objective::eMin );

 // Set Block objective
 this->set_objective( &objective , eNoMod );

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
  && ColVariable::is_feasible( v_reactive_power , tol )
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
  && RowConstraint::is_feasible( VolumetricBounds_Const , tol , rel_viol )
  && RowConstraint::is_feasible( Reference_Schedule_Const , tol , rel_viol )
  && RowConstraint::is_feasible( ReactivePower_Bound_Const , tol , rel_viol ) );

 } // end( HydroUnitBlock::is_feasible )

/*--------------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/

Solution * HydroUnitBlock::get_Solution( Configuration * csolc ,
					 bool emptys )
{
 Index wsol = 63;
 if( ( ! csolc ) && f_BlockConfig )
  csolc = f_BlockConfig->f_solution_Configuration;

 if( auto config = dynamic_cast< SimpleConfiguration< int > * >( csolc ) )
  wsol = config->f_value;

 // call the method of the base class
 auto * sol = dynamic_cast< HydroUnitBlockSolution * >(
		                  UnitBlock::get_Solution( csolc , emptys ) );
 assert( sol );

 if( wsol & 16 ) {
  auto sz = boost::multi_array< double , 2 >::extent_gen()
                           [ get_number_reservoirs() ][ get_time_horizon() ];
  sol->v_volume.resize( sz );
  }

 if( wsol & 32 ) {
  auto sz = boost::multi_array< double , 2 >::extent_gen()
                           [ get_number_generators() ][ get_time_horizon() ];
  sol->v_flow.resize( sz );
  }

 if( ! emptys )
  sol->read( this );

 return( sol );
 }

/*--------------------------------------------------------------------------*/

UnitBlockSolution * HydroUnitBlock::new_Solution( void ) const {
 return( new HydroUnitBlockSolution() );
 }

/*--------------------------------------------------------------------------*/
/*-------- METHODS FOR LOADING, PRINTING & SAVING THE HydroUnitBlock -------*/
/*--------------------------------------------------------------------------*/

void HydroUnitBlock::serialize( netCDF::NcGroup & group ) const
{
 UnitBlock::serialize( group );

 auto TimeHorizon = group.getDim( "TimeHorizon" );

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

 ::serialize( group , "ActivePowerCost" , netCDF::NcDouble() ,
              NumberArcs , v_ActivePowerCost , false );

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
              { NumberArcs , TimeHorizon } , v_InertiaPower );

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

  if( *( v_inflows.data() + i ) != *( values++ ) )
   identical = false;
 }

 if( identical )
  return;

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  for( auto i : subset ) {
   Index t = i % f_time_horizon;
   Index r = i / f_time_horizon;
   v_inflows[ r ][ t ] = *( values++ );
  }

  if( constraints_generated() )
   // Change the abstract representation

   for( auto i : subset ) {
    Index t = i % f_time_horizon;
    Index r = i / f_time_horizon;

    if( t == 0 ) {
     const auto volume = v_InitialVolumetric[ r ] >= 0. ?
                          v_InitialVolumetric[ r ] : 0.;
     FinalVolumeReservoir_Const[ t ][ r ].set_both(
      volume + v_inflows[ r ][ t ] , issueAMod );
    }
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

    if( t == 0 ) {
     const auto volume = v_InitialVolumetric[ r ] >= 0. ?
                          v_InitialVolumetric[ r ] : 0.;
     FinalVolumeReservoir_Const[ t ][ r ].set_both(
      volume + v_inflows[ r ][ t ] , issueAMod );
    }
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

  v_InertiaPower.resize( boost::extents[ f_NumberArcs ][ f_time_horizon ] );
 }

 // If nothing changes, return
 bool identical = true;
 auto values_it = values;
 for( auto i : subset ) {
  if( i >= f_NumberArcs * f_time_horizon )
   throw( std::invalid_argument( "HydroUnitBlock::set_inertia_power: "
                                 "invalid value in subset." ) );

  Index a = i / f_time_horizon;
  Index t = i % f_time_horizon;

  if( v_InertiaPower[ a ][ t ] != *( values_it++ ) )
   identical = false;
 }

 if( identical )
  return;

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation
  values_it = values;
  for( auto i : subset ) {
   Index a = i / f_time_horizon;
   Index t = i % f_time_horizon;
   v_InertiaPower[ a ][ t ] = *( values_it++ );
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
                            this , HydroUnitBlockMod::eSetInerP ,
                            std::move( subset ) ) ,
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

  v_InertiaPower.resize( boost::extents[ f_NumberArcs ][ f_time_horizon ] );
 }

 // If nothing changes, return
 bool identical = true;
 auto values_it = values;
 for( Index i = rng.first ; i < rng.second ; ++i ) {
  Index a = i / f_time_horizon;
  Index t = i % f_time_horizon;

  if( v_InertiaPower[ a ][ t ] != *( values_it++ ) ) {
   identical = false;
   break;
  }
 }

 if( identical )
  return;

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation
  values_it = values;
  for( Index i = rng.first ; i < rng.second ; ++i ) {
   Index a = i / f_time_horizon;
   Index t = i % f_time_horizon;
   v_InertiaPower[ a ][ t ] = *( values_it++ );
  }

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

  const auto volume = *( values++ );
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
   const auto volume = v_InitialVolumetric[ r ] >= 0. ?
                         v_InitialVolumetric[ r ] : 0.;
   FinalVolumeReservoir_Const[ 0 ][ r ].set_both
    ( volume + v_inflows[ r ][ 0 ] , issueAMod );
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
    const auto volume = v_InitialVolumetric[ r ] >= 0. ?
                          v_InitialVolumetric[ r ] : 0.;
    FinalVolumeReservoir_Const[ 0 ][ r ].set_both
     ( volume + v_inflows[ r ][ 0 ] , issueAMod );
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

void HydroUnitBlock::update_initial_flow_rate_in_cnstrs( Block::Range arcs ,
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
  v_InitialFlowRate.resize( max_index + 1 );
 }

 bool identical = true;
 for( auto i : subset ) {
  if( i >= v_InitialFlowRate.size() )
   throw( std::invalid_argument( "HydroUnitBlock::set_initial_flow_rate: "
                                  "invalid value in subset." ) );
  auto flow_rate = *( values++ );
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

void HydroUnitBlock::set_active_power_cost( MF_dbl_it values ,
                                            Subset && subset ,
                                            bool ordered ,
                                            c_ModParam issuePMod ,
                                            c_ModParam issueAMod )
{
 if( subset.empty() )
  return;

 if( v_ActivePowerCost.empty() ) {
  if( std::all_of( values ,
                   values + subset.size() ,
                   []( double cst ) { return( cst == 0.0 ); } ) )
   return;

  v_ActivePowerCost.assign( f_NumberArcs , 0.0 );
 }

 for( auto arc : subset )
  if( arc >= v_ActivePowerCost.size() )
   throw( std::invalid_argument(
    "HydroUnitBlock::set_active_power_cost: invalid index in subset." ) );

 auto values_it = values;

 bool identical = true;
 for( auto arc : subset ) {
  if( v_ActivePowerCost[ arc ] != *( values_it++ ) ) {
   identical = false;
   break;
  }
 }

 if( identical )
  return;

 if( not_dry_run( issuePMod ) ) {
  values_it = values;
  for( auto arc : subset )
   v_ActivePowerCost[ arc ] = *( values_it++ );

  if( not_dry_run( issueAMod ) && objective_generated() ) {
   auto * lf = static_cast< LinearFunction * >( objective.get_function() );

   for( auto arc : subset ) {
    for( Index t = 0 ; t < f_time_horizon ; ++t ) {
     const auto idx = lf->is_active( &v_active_power[ arc ][ t ] );

     if( idx == Inf< Index >() )
      throw( std::logic_error(
       "HydroUnitBlock::set_active_power_cost: expected Variable not "
       "found in objective." ) );

     lf->modify_coefficient( idx ,
                             v_ActivePowerCost[ arc ] ,
                             issueAMod );
    }
   }
  }
 }

 if( issue_pmod( issuePMod ) ) {
  if( ! ordered )
   std::sort( subset.begin() , subset.end() );

  Block::add_Modification(
   std::make_shared< HydroUnitBlockSbstMod >(
    this , HydroUnitBlockMod::eSetActPCost , std::move( subset ) ) ,
   Observer::par2chnl( issuePMod ) );
 }
}  // end( HydroUnitBlock::set_active_power_cost( subset ) )

/*--------------------------------------------------------------------------*/

void HydroUnitBlock::set_active_power_cost( MF_dbl_it values ,
                                            Range rng ,
                                            c_ModParam issuePMod ,
                                            c_ModParam issueAMod )
{
 rng.second = std::min( rng.second , f_NumberArcs );
 if( rng.second <= rng.first )
  return;

 c_Index sz = rng.second - rng.first;

 if( v_ActivePowerCost.empty() ) {
  if( std::all_of( values ,
                   values + sz ,
                   []( double cst ) { return( cst == 0.0 ); } ) )
   return;

  v_ActivePowerCost.assign( f_NumberArcs , 0.0 );
 }

 if( std::equal( values ,
                 values + sz ,
                 v_ActivePowerCost.begin() + rng.first ) )
  return;

 if( not_dry_run( issuePMod ) ) {
  std::copy( values ,
             values + sz ,
             v_ActivePowerCost.begin() + rng.first );

  if( not_dry_run( issueAMod ) && objective_generated() ) {
   auto * lf = static_cast< LinearFunction * >( objective.get_function() );

   for( Index arc = rng.first ; arc < rng.second ; ++arc ) {
    for( Index t = 0 ; t < f_time_horizon ; ++t ) {
     const auto idx = lf->is_active( &v_active_power[ arc ][ t ] );

     if( idx == Inf< Index >() )
      throw( std::logic_error(
       "HydroUnitBlock::set_active_power_cost: expected Variable not "
       "found in objective." ) );

     lf->modify_coefficient( idx ,
                             v_ActivePowerCost[ arc ] ,
                             issueAMod );
    }
   }
  }
 }

 if( issue_pmod( issuePMod ) )
  Block::add_Modification(
   std::make_shared< HydroUnitBlockRngdMod >(
    this , HydroUnitBlockMod::eSetActPCost , rng ) ,
   Observer::par2chnl( issuePMod ) );

}  // end( HydroUnitBlock::set_active_power_cost( range ) )

/*--------------------------------------------------------------------------*/
/*------------------ METHODS OF HydroUnitBlockSolution ---------------------*/
/*--------------------------------------------------------------------------*/

void HydroUnitBlockSolution::deserialize( const netCDF::NcGroup & group )
{
 // call the method of the base class
 UnitBlockSolution::deserialize( group );

 if( ! deserialize_dim( group , "NumberReservoirs" , f_reservoirs ,
			true ) )
  f_reservoirs = 1;

 using index = boost::multi_array< double , 2 >::index;
 const std::vector< index > empty = { 0 , 0 };
 const std::vector< index > full = { 1 , f_time_horizon };

 if( f_reservoirs == 1 ) {
  auto ncVar = group.getVar( "VolumetricLevel" );
  if( ncVar.isNull() )
   v_volume.resize( empty );
  else {
   v_volume.resize( full );
   ncVar.getVar( { 0 } , { f_time_horizon } , v_volume.data() );
   }
  }
 else
  ::deserialize< double , 2 >( group , "VolumetricLevel" ,
                               { f_reservoirs , f_time_horizon } ,
                               v_volume , true , true );

 if( f_number_generators == 1 ) {
  auto ncVar = group.getVar( "VolumetricFlow" );
  if( ncVar.isNull() )
   v_flow.resize( empty );
  else {
   v_flow.resize( full );
   ncVar.getVar( { 0 } , { f_time_horizon } , v_flow.data() );
   }
  }
 else
  ::deserialize< double , 2 >( group , "VolumetricFlow" ,
                               { f_number_generators , f_time_horizon } ,
                               v_flow , true , true );

 }  // end( HydroUnitBlockSolution::deserialize )

/*--------------------------------------------------------------------------*/

void HydroUnitBlockSolution::read( const Block * block )
{
 auto HUB = dynamic_cast< const HydroUnitBlock * >( block );
 if( ! HUB )
  throw( std::invalid_argument( "HydroUnitBlockSolution::read: block is "
				"not a HydroUnitBlock" ) );

 UnitBlockSolution::read( HUB );  // call the method of the base class

 f_reservoirs = HUB->get_number_reservoirs();

 if( ! v_volume.empty() )
  for( Index i = 0 ; i < f_reservoirs ; ++i ) {
   auto Vi = HUB->get_const_volumetric( i );
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    v_volume[ i ][ t ] = Vi[ t ].get_value();
   }

 if( ! v_flow.empty() )
  for( Index i = 0 ; i < f_number_generators ; ++i )
   if( auto Fi = HUB->get_const_flow_rate( i ) )
    for( Index t = 0 ; t < f_time_horizon ; ++t )
     v_flow[ i ][ t ] = Fi[ t ].get_value();

 }  // end( HydroUnitBlockSolution::read )

/*--------------------------------------------------------------------------*/

void HydroUnitBlockSolution::write( Block * block )
{
 UnitBlockSolution::write( block );  // call the method of the base class

 auto HUB = dynamic_cast< HydroUnitBlock * >( block );
 if( ! HUB )
  throw( std::invalid_argument( "HydroUnitBlockSolution::write: block is "
				"not a HydroUnitBlock" ) );

 if(  f_reservoirs != HUB->get_number_reservoirs() )
  throw( std::invalid_argument( "HydroUnitBlockSolution::write: "
				"inconsistent number of reservoirs" ) );

 if( ! v_volume.empty() )
  for( Index i = 0 ; i < f_reservoirs ; ++i ) {
   auto Vi = HUB->get_volumetric( i );
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    Vi[ t ].set_value( v_volume[ i ][ t ] );
   }

 if( ! v_flow.empty() )
  for( Index i = 0 ; i < f_number_generators ; ++i ) {
   auto Fi = HUB->get_flow_rate( i );
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    Fi[ t ].set_value( v_flow[ i ][ t ] );
   }

 }  // end( HydroUnitBlockSolution::write )

/*--------------------------------------------------------------------------*/

void HydroUnitBlockSolution::serialize( netCDF::NcGroup & group ) const
{
 UnitBlockSolution::serialize( group );  // call the method of the base class

 // recover the just serialized time horizon
 netCDF::NcDim th = group.getDim( "TimeHorizon" );

 if( f_reservoirs > 1 ) {
  auto nr = group.addDim( "NumberReservoirs" , f_reservoirs );
  ::serialize< double , 2 >( group , "VolumetricLevel" , netCDF::NcDouble() ,
			     { nr , th } , v_volume );
  }
 else
  if( ! v_volume.empty() )
   group.addVar( "VolumetricLevel" , netCDF::NcDouble() , th ).putVar(
			      { 0 } , { f_time_horizon } , v_volume.data() );


 if( f_number_generators > 1 ) {
  // recover the just serialized number of generators
  netCDF::NcDim ng = group.getDim( "NumberGenerators" );
  ::serialize< double , 2 >( group , "VolumetricFlow" , netCDF::NcDouble() ,
			     { ng , th } , v_flow );
  }
 else
  if( ! v_flow.empty() )
   group.addVar( "VolumetricFlow" , netCDF::NcDouble() , th ).putVar(
			        { 0 } , { f_time_horizon } , v_flow.data() );

 }  // end( HydroUnitBlockSolution::serialize )

/*--------------------------------------------------------------------------*/

HydroUnitBlockSolution * HydroUnitBlockSolution::scale( double factor ) const
{
 auto sol = clone();

 if( factor == 1 )
  return( sol );

 guts_of_scale( sol , factor );

 if( ! v_volume.empty() )
  for( Index i = 0 ; i < f_reservoirs ; ++i )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    sol->v_volume[ i ][ t ] *= factor;

 if( ! v_flow.empty() )
  for( Index i = 0 ; i < f_number_generators ; ++i )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    sol->v_flow[ i ][ t ] *= factor;

 return( sol );

 }  // end( HydroUnitBlockSolution::scale )

/*--------------------------------------------------------------------------*/

void HydroUnitBlockSolution::sum( const Solution * solution ,
				  double multiplier )
{
 // call the method of the base class
 UnitBlockSolution::sum( solution , multiplier );

 auto HUBS = dynamic_cast< const HydroUnitBlockSolution * >( solution );
 if( ! HUBS )
  throw( std::invalid_argument( "HydroUnitBlockSolution::sum: solution not "
				"a HydroUnitBlockSolution" ) );

 if( f_reservoirs != HUBS->f_reservoirs )
  throw( std::invalid_argument( "HydroUnitBlockSolution::sum: inconsistent "
				"number of reservoirs" ) );

 if( ! v_volume.empty() )
  for( Index i = 0 ; i < f_reservoirs ; ++i )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    v_volume[ i ][ t ] += HUBS->v_volume[ i ][ t ] * multiplier;

 if( ! v_flow.empty() )
  for( Index i = 0 ; i < f_number_generators ; ++i )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    v_flow[ i ][ t ] += HUBS->v_flow[ i ][ t ] * multiplier;

 }  // end( HydroUnitBlockSolution::sum )

/*--------------------------------------------------------------------------*/

HydroUnitBlockSolution * HydroUnitBlockSolution::clone( bool empty ) const
{
 auto * sol = new HydroUnitBlockSolution();

 if( ! empty ) {
  guts_of_clone( sol );
  copy_multi_array( sol->v_volume , v_volume );
  copy_multi_array( sol->v_flow , v_flow );
  }

 return( sol );

 }  // end( HydroUnitBlockSolution::clone )

/*--------------------------------------------------------------------------*/
/*------------------- End File HydroUnitBlock.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
