/*--------------------------------------------------------------------------*/
/*------------------------- File HeatBlock.cpp -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the HeatBlock class.
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

#include <iostream>

#include <map>

#include <random>

#include "FRealObjective.h"

#include "FRowConstraint.h"

#include "HeatBlock.h"

#include "LinearFunction.h"

#include "UCBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register HeatBlock to the Block factory

SMSpp_insert_in_factory_cpp_1( HeatBlock );

/*--------------------------------------------------------------------------*/
/*--------------------------- METHODS OF HeatBlock -------------------------*/
/*--------------------------------------------------------------------------*/

HeatBlock::~HeatBlock()
{
 Constraint::clear( v_EvolutionStoredHeat_Const );
 Constraint::clear( v_HeatDemand_Const );

 Constraint::clear( v_HeatBounds_Const );
 Constraint::clear( v_HeatStorageBounds_Const );

 objective.clear();
}

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void HeatBlock::deserialize_time_horizon( const netCDF::NcGroup & group )
{
 netCDF::NcDim TimeHorizon = group.getDim( "TimeHorizon" );
 if( TimeHorizon.isNull() ) {
  // dimension TimeHorizon is not present in the netCDF input

  if( f_time_horizon == 0 ) {
   auto f_B = dynamic_cast< UCBlock * >( get_f_Block() );
   if( f_B )
    // The father Block is available. Take time horizon from it.
    this->set_time_horizon( f_B->get_time_horizon() );
   else
    throw( std::invalid_argument(
     "HeatBlock::deserialize: TimeHorizon is not present in the "
     "netCDF input and HeatBlock does not have a father." ) );
  }
 } else {
  // dimension TimeHorizon is present in the netCDF input

  auto th = TimeHorizon.getSize();
  if( f_time_horizon == 0 )
   this->set_time_horizon( th );
  else if( f_time_horizon != th )
   throw( std::logic_error(
    "HeatBlock::deserialize: TimeHorizon is not present in the "
    "netCDF. The (nonzero) time horizon of HeatBlock is different "
    "from that of its father, but they should be equal." ) );
 }
}

/*--------------------------------------------------------------------------*/

void HeatBlock::deserialize_change_intervals( const netCDF::NcGroup & group )
{
 auto NumberIntervals = group.getDim( "NumberIntervals" );
 if( NumberIntervals.isNull() )
  f_number_intervals = 0;
 else {
  f_number_intervals = NumberIntervals.getSize();
  if( ( f_number_intervals < 1 ) || ( f_number_intervals > f_time_horizon ) )
   throw( std::invalid_argument(
    "HeatBlock::deserialize: invalid NumberIntervals. "
    "It must be between 1 and TimeHorizon." ) );
 }

 if( ( f_number_intervals > 1 ) && ( f_number_intervals < f_time_horizon ) ) {

  ::deserialize( group , "ChangeIntervals" , f_number_intervals ,
                 v_change_intervals );

  // Check that all numbers are between 1 and f_time_horizon, that
  // the last number is == f_time_horizon, and that they are ordered
  // in increasing sense

  if( v_change_intervals.back() != f_time_horizon )
   throw( std::invalid_argument(
    "HeatBlock::deserialize: invalid value in ChangeIntervals: "
    "the last element must be TimeHorizon." ) );

  Index previous_t = 0;

  for( auto t : v_change_intervals ) {
   if( ! ( t > previous_t && t < f_time_horizon - 1 ) )
    throw( std::invalid_argument(
     "HeatBlock::deserialize: invalid value in ChangeIntervals: " +
     std::to_string( t ) + ". All values must be between 1 and " +
     "TimeHorizon and in strictly increasing order." ) );

   previous_t = t;
  }
 }
}

/*--------------------------------------------------------------------------*/

void HeatBlock::deserialize( const netCDF::NcGroup & group )
{

#ifndef NDEBUG
 static std::vector< std::string > expected_dims = {};
 check_dimensions( group , expected_dims , std::cerr );

 static std::vector< std::string > expected_vars = {};
 check_variables( group , expected_vars , std::cerr );
#endif

 deserialize_time_horizon( group );
 deserialize_change_intervals( group );

 ::deserialize( group , "TotalHeatDemand" , f_time_horizon ,
                v_heat_demand );
/*
    ::deserialize( group , "CostHeatUnit",
                   { f_number_intervals , f_number_heat_units },
                   v_cost_heat_unit );

    ::deserialize( group , "MinHeatProduction",
                   { f_number_intervals , f_number_heat_units },
                   v_min_heat_production );

    ::deserialize( group , "MaxHeatProduction",
                   { f_number_intervals , f_number_heat_units },
                   v_max_heat_production ); */

 ::deserialize( group , "MinHeatStorage" , f_number_intervals ,
                v_min_heat_storage );

 ::deserialize( group , "MaxHeatStorage" , f_number_intervals ,
                v_max_heat_storage );

 ::deserialize( group , f_storing_heat_rho , "StoringHeatRho" );
 ::deserialize( group , f_extracting_heat_rho , "ExtractingHeatRho" );
 ::deserialize( group , f_keeping_heat_rho , "KeepingHeatRho" );
 ::deserialize( group , f_initial_heat_storage , "InitialHeatAvailable" );

 Block::deserialize( group );

}  // end( HeatBlock::deserialize )

/*--------------------------------------------------------------------------*/

void HeatBlock::generate_abstract_variables( Configuration * stvv )
{
 if( ! v_heat.empty() )  // variables have already been generated
  return;                // nothing to do

 if( f_time_horizon == 0 )
  // there are no variables to be generated
  return;

 // Heat variables

 v_heat.resize( boost::extents[ f_time_horizon ][ f_number_heat_units ] );
 for( Index t = 0 ; t < f_time_horizon ; ++t )
  for( Index g = 0 ; t < f_number_heat_units ; ++g )
   v_heat[ t ][ g ].set_type( ColVariable::kNonNegative );
 add_static_variable( v_heat );

 // Heat added, removed, and available variables

 auto variables_and_types =
  { std::make_pair( v_heat_added , ColVariable::kNonNegative ) ,
    std::make_pair( v_heat_removed , ColVariable::kNonNegative ) ,
    std::make_pair( v_heat_available , ColVariable::kNonNegative )
  };

 // informs which variables must be generated
 int variables_to_be_generated = 0;
 if( ( ! stvv ) && f_BlockConfig )
  stvv = f_BlockConfig->f_static_variables_Configuration;
 if( auto sci = dynamic_cast< SimpleConfiguration< int > * >( stvv ) )
  variables_to_be_generated = sci->f_value;

 unsigned int k = 1;
 for( auto pair : variables_and_types ) {
  if( variables_to_be_generated & k ) {
   auto variables = pair.first;
   variables.resize( f_time_horizon );
   for( auto & variable : variables )
    variable.set_type( pair.second );
   add_static_variable( variables );
  }
  k *= 2;
 }
}  // end( HeatBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void HeatBlock::generate_abstract_constraints( Configuration * stcc )
{
 if( ! v_HeatBounds_Const.empty() )  // constraints have already been generated
  return;                            // nothing to do

 LinearFunction::v_coeff_pair vars;

 // Satisfaction Heat Bounds constraints

 v_HeatBounds_Const.resize(
  boost::multi_array< BoxConstraint , 2 >::
  extent_gen()[ f_time_horizon ][ f_number_heat_units ] );

 for( Index t = 0 ; t < f_time_horizon ; ++t )
  for( Index unit_id = 0 ; unit_id < f_number_heat_units ; ++unit_id ) {
   v_HeatBounds_Const[ t ][ unit_id ].set_lhs(
    v_min_heat_production[ t ][ unit_id ] );
   v_HeatBounds_Const[ t ][ unit_id ].set_rhs(
    v_max_heat_production[ t ][ unit_id ] );
   v_HeatBounds_Const[ t ][ unit_id ].set_variable(
    &v_heat[ t ][ unit_id ] );
  }

 add_static_constraint( v_HeatBounds_Const );

 // Satisfaction Heat Storage Bounds constraints

 if( ! v_heat_available.empty() ) {

  v_HeatStorageBounds_Const.resize( f_time_horizon );

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   v_HeatStorageBounds_Const[ t ].set_lhs( v_min_heat_storage[ t ] );
   v_HeatStorageBounds_Const[ t ].set_rhs( v_max_heat_storage[ t ] );
   v_HeatStorageBounds_Const[ t ].set_variable( &v_heat_available[ t ] );
  }

  add_static_constraint( v_HeatStorageBounds_Const );
 }

 // Evolution Stored Heat Constraints

 if( ! ( v_heat_available.empty() || v_heat_added.empty() ||
         v_heat_removed.empty() ) ) {

  v_EvolutionStoredHeat_Const.resize( f_time_horizon );

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {

   if( t == 0 ) {

    vars.push_back( std::make_pair( &v_heat_available[ t ] , 1.0 ) );
    vars.push_back( std::make_pair( &v_heat_added[ t ] ,
                                    -f_storing_heat_rho ) );
    vars.push_back( std::make_pair( &v_heat_removed[ t ] ,
                                    f_extracting_heat_rho ) );

    v_EvolutionStoredHeat_Const[ t ].set_lhs( 0.0 );
    v_EvolutionStoredHeat_Const[ t ].set_rhs(
     f_initial_heat_storage * f_keeping_heat_rho );

   } else {

    // TODO t+1 is not defined for t = time_horizon - 1
    vars.push_back( std::make_pair( &v_heat_available[ t + 1 ] , 1.0 ) );
    vars.push_back( std::make_pair( &v_heat_available[ t ] ,
                                    -f_keeping_heat_rho ) );
    vars.push_back( std::make_pair( &v_heat_added[ t + 1 ] ,
                                    -f_storing_heat_rho ) );
    vars.push_back( std::make_pair( &v_heat_removed[ t + 1 ] ,
                                    f_extracting_heat_rho ) );

    v_EvolutionStoredHeat_Const[ t ].set_both( 0.0 );
   }

   v_EvolutionStoredHeat_Const[ t ].set_function(
    new LinearFunction( std::move( vars ) ) );
  }

  add_static_constraint( v_EvolutionStoredHeat_Const );
 }

 // Satisfaction Heat Demand constraints

 if( ( ! v_heat_added.empty() ) && ( ! v_heat_removed.empty() ) ) {

  v_HeatDemand_Const.resize( f_time_horizon );

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {

   vars.push_back( std::make_pair( &v_heat_added[ t ] , -1.0 ) );
   vars.push_back( std::make_pair( &v_heat_removed[ t ] , 1.0 ) );
   for( Index unit_id = 0 ; unit_id < f_number_heat_units ; ++unit_id )
    vars.push_back( std::make_pair( &v_heat[ t ][ unit_id ] , 1.0 ) );

   v_HeatDemand_Const[ t ].set_function(
    new LinearFunction( std::move( vars ) ) );
   v_HeatDemand_Const[ t ].set_lhs( v_heat_demand[ t ] );
   v_HeatDemand_Const[ t ].set_rhs( Inf< double >() );
  }

  add_static_constraint( v_HeatDemand_Const );
 }
}  // end( HeatBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void HeatBlock::generate_objective( Configuration * objc )
{
 if( get_objective() )  // Objective has already been generated
  return;               // nothing to do

 if( v_heat.size() != f_time_horizon )
  throw( std::logic_error( "HeatBlock::generate_objective: v_heat must have "
                           "size equal to the time horizon." ) );

 LinearFunction::v_coeff_pair vars;

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {
  for( Index unit_id = 0 ; unit_id < f_number_heat_units ; ++unit_id ) {
   auto cost = get_cost_heat_unit()[ t ][ unit_id ];
   vars.push_back( std::make_pair( &v_heat[ t ][ unit_id ] , cost ) );
  }
 }

 objective.set_function( new LinearFunction( std::move( vars ) ) );
 objective.set_sense( Objective::eMin );

 // Set block objective
 this->set_objective( &objective );

}  // end( HeatBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*---------- METHODS FOR LOADING, PRINTING & SAVING THE HeatBlock ----------*/
/*--------------------------------------------------------------------------*/

void HeatBlock::serialize( netCDF::NcGroup & group ) const
{
 Block::serialize( group );

 group.addDim( "TimeHorizon" , f_time_horizon );

 auto NumberIntervals = group.addDim( "NumberIntervals" , f_number_intervals );

 ::serialize( group , "ChangeInterval" , netCDF::NcUint64() ,
              NumberIntervals , v_change_intervals );

 auto dim_number_units = group.addDim( "NumberHeatUnits" ,
                                       f_number_heat_units );

 ::serialize( group , "StoringHeatRho" , netCDF::NcDouble() ,
              f_storing_heat_rho );

 ::serialize( group , "ExtractingHeatRho" , netCDF::NcDouble() ,
              f_extracting_heat_rho );

 ::serialize( group , "StoringHeatRho" , netCDF::NcDouble() ,
              f_keeping_heat_rho );

 ::serialize( group , "InitialHeatAvailable" , netCDF::NcDouble() ,
              f_initial_heat_storage );

 ::serialize( group , "ChangeInterval" , netCDF::NcUint64() ,
              { NumberIntervals } , v_change_intervals );

 ::serialize( group , "TotalHeatDemand" , netCDF::NcDouble() ,
              { NumberIntervals } , v_heat_demand );

 ::serialize( group , "MinHeatStorage" , netCDF::NcDouble() ,
              { NumberIntervals } , v_min_heat_storage );

 ::serialize( group , "MaxHeatStorage" , netCDF::NcDouble() ,
              { NumberIntervals } , v_max_heat_storage );
/*
    ::serialize( group , "MinHeatProduction" , netCDF::NcDouble() ,
                 {NumberIntervals  , dim_number_units} ,
                 v_min_heat_production);

    ::serialize( group , "MaxHeatProduction" , netCDF::NcDouble() ,
                 {NumberIntervals  , dim_number_units} ,
                 v_max_heat_production); */

}  // end( HeatBlock::serialize )

/*--------------------------------------------------------------------------*/
/*----------------------- End File HeatBlock.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
