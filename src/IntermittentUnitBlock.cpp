/*--------------------------------------------------------------------------*/
/*-------------------- File IntermittentUnitBlock.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the IntermittentUnitBlock class.
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
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Ali Ghezelsoflu,
 *                      Rafael Durbano Lobato, Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <iostream>

#include <random>

#include <map>

#include "IntermittentUnitBlock.h"

#include "LinearFunction.h"

#include "FRealObjective.h"

#include "UnitBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register IntermittentUnitBlock to the Block factory
SMSpp_insert_in_factory_cpp_1( IntermittentUnitBlock );

// register IntermittentUnitBlockSolution to the Solution factory
SMSpp_insert_in_factory_cpp_0( IntermittentUnitBlockSolution );

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS OF IntermittentUnitBlock -------------------*/
/*--------------------------------------------------------------------------*/

IntermittentUnitBlock::~IntermittentUnitBlock()
{
 Constraint::clear( min_power_Const );
 Constraint::clear( max_power_Const );
 Constraint::clear( active_power_bounds_design_Const );
 Constraint::clear( active_power_bounds_Const );

 design_bound_Const.clear();

 objective.clear();
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::deserialize( const netCDF::NcGroup & group )
{
 // Deserialize data from the base class
 UnitBlock::deserialize( group );

 // Mandatory variables

 ::deserialize( group , "MaxPower" , f_time_horizon , v_MaxPower ,
                false , true , v_change_intervals );

 // Optional variables

 if( ::deserialize( group , f_InvestmentCost , "InvestmentCost" ) ) {
  ::deserialize( group , f_MinCapacityDesign , "MinCapacityDesign" );
  ::deserialize( group , f_MaxCapacityDesign , "MaxCapacityDesign" );
  }

 if( ! ::deserialize( group , "MinPower" , f_time_horizon , v_MinPower ,
                      true , true , v_change_intervals ) )
  v_MinPower.resize( f_time_horizon );

 if( ! ::deserialize( group , "InertiaPower" , f_time_horizon ,
		      v_InertiaPower , true , true , v_change_intervals ) )
  v_InertiaPower.resize( f_time_horizon );

 if( ! ::deserialize( group , "ActivePowerCost" , f_time_horizon ,
                      v_ActivePowerCost , true , true , v_change_intervals ) )
  v_ActivePowerCost.resize( f_time_horizon );

 if ( ! ::deserialize( group , f_MaxGeneration, "MaxGeneration" ))
  f_MaxGeneration = Inf< double >();

 if ( ! ::deserialize( group , f_MinGeneration , "MinGeneration" ))
  f_MinGeneration = -Inf< double >();

 ::deserialize( group , f_gamma , "Gamma" );

 ::deserialize( group , f_kappa , "Kappa" );

 // variables for AC elements
 if( ::deserialize( group , "MaxReactivePower" , f_time_horizon ,
		    v_MaxReactivePower , true , true , v_change_intervals ) )
  if( std::all_of( v_MaxReactivePower.begin() , v_MaxReactivePower.end() ,
		   []( double i ) { return( i == 0 ); } ) )
   v_MaxReactivePower.clear();

 if( ::deserialize( group , "MinReactivePower" , f_time_horizon ,
		    v_MinReactivePower , true , true , v_change_intervals ) )
  if( std::all_of( v_MinReactivePower.begin() , v_MinReactivePower.end() ,
		   []( double i ) { return( i == 0 ); } ) )
   v_MinReactivePower.clear();

 if( f_max_power_epsilon > 0 )
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   if( v_MaxPower[ t ] == 0.0 )
    v_MaxPower[ t ] = f_max_power_epsilon;

 check_data_consistency();

 }  // end( IntermittentUnitBlock::deserialize )

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG

/*
std::vector< std::string > IntermittentUnitBlock::expected_dims( void )
 const {
 static const std::vector< std::string > ed = { };

 auto ret = UnitBlock::expected_dims();
 ret.insert( ret.end() , ed.begin() , ed.end() );

 return( ret );
 }

----------------------------------------------------------------------------*/

std::vector< std::string > IntermittentUnitBlock::expected_vars( void )
 const {
 static const std::vector< std::string > ev =
 { "InvestmentCost" , "MinCapacityDesign" , "MaxCapacityDesign" ,
   "MaxCapacity" , "MinPower" , "MaxPower" , "InertiaPower" ,
   "ActivePowerCost", "Gamma" , "Kappa", "MinReactivePower",
   "MaxReactivePower", "MaxGeneration", "MinGeneration"
   };

 auto ret = UnitBlock::expected_vars();
 ret.insert( ret.end() , ev.begin() , ev.end() );

 return( ret );
 }

#endif

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::check_data_consistency( void ) const
{
 // Min/Max capacity design

 if( f_MinCapacityDesign < 0 )
  throw( std::logic_error( "IntermittentUnitBlock::check_data_consistency: "
                           "MinCapacityDesign must be nonnegative." ) );

 // Continue case (MaxCapacityDesign > 0): MinCapacityDesign <= MaxCapacityDesign
 if( ( f_MaxCapacityDesign > 0 ) && ( f_MinCapacityDesign > f_MaxCapacityDesign ) )
  throw( std::logic_error( "IntermittentUnitBlock::check_data_consistency: "
                           "MinCapacityDesign > MaxCapacityDesign." ) );

 // Unitary case (|MaxCapacityDesign| == 1): MinCapacityDesign <= 1
 if( ( std::abs( f_MaxCapacityDesign ) == 1 ) && ( f_MinCapacityDesign > 1.0 ) )
  throw( std::logic_error( "IntermittentUnitBlock::check_data_consistency: "
                           "MinCapacityDesign must be <= 1 when |MaxCapacityDesign| == 1." ) );

 // Binary case (max < 0): MinCapacityDesign <= 1
 if( ( f_MaxCapacityDesign < 0 ) && ( f_MinCapacityDesign > 1.0 ) )
  throw( std::logic_error( "IntermittentUnitBlock::check_data_consistency: "
                           "MinCapacityDesign must be <= 1 for binary design." ) );

 // Minimum and maximum power

 assert( v_MinPower.size() == f_time_horizon );
 assert( v_MaxPower.size() == f_time_horizon );

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {
  if( v_MinPower[ t ] > v_MaxPower[ t ] )
   throw( std::logic_error( "IntermittentUnitBlock::check_data_consistency: "
                             "minimum power at time " + std::to_string( t ) +
                             " is " + std::to_string( v_MinPower[ t ] ) +
                             ", which is greater than the maximum power, which "
                             "is " + std::to_string( v_MaxPower[ t ] ) + "." ) );

 }
 
 // Maximum generation no lower than Minimum Generation
 if( ( f_MaxGeneration < f_MinGeneration ) )
  throw( std::logic_error( "IntermittentUnitBlock::check_data_consistency: "
                           "MaxGeneration must be >= MinGeneration." ) );

 // Gamma

 if( ( f_gamma < 0 ) || ( f_gamma > 1 ) )
  throw( std::logic_error( "IntermittentUnitBlock::check_data_consistency: "
                            "gamma must be between 0 and 1, but it is " +
                            std::to_string( f_gamma ) + "." ) );

 // Kappa

 if( f_kappa < 0 )
  throw( std::logic_error( "IntermittentUnitBlock::check_data_consistency: "
                            "kappa must be nonnegative, but it is" +
                            std::to_string( f_kappa ) + "." ) );

 if( ! v_InertiaPower.empty() ) {
  assert( v_InertiaPower.size() == f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   if( v_InertiaPower[ t ] < 0 )
    throw( std::logic_error( "IntermittentUnitBlock::check_data_consistency: "
                              "inertia power for time " + std::to_string( t ) +
                              " must be nonnegative, but it is" +
                              std::to_string( v_InertiaPower[ t ] ) + "." ) );
 }

}  // end( IntermittentUnitBlock::check_data_consistency )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::generate_abstract_variables( Configuration * stvv )
{
 if( variables_generated() )  // variables have already been generated
  return;                     // nothing to do

 UnitBlock::generate_abstract_variables( stvv );

 // Design Variable
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( f_InvestmentCost != 0 ) {
  if( f_MaxCapacityDesign < 0 )
   design.set_type( ColVariable::kBinary );
  else
   design.set_type( ColVariable::kNonNegative );
  add_static_variable( design , "x_intermittent" );
  }
 else
  design.set_value( std::numeric_limits< double >::quiet_NaN() );

 // Active Power Variable - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 v_active_power.resize( f_time_horizon );
 for( auto & var : v_active_power )
  var.set_type( ColVariable::kContinuous );
 add_static_variable( v_active_power , "p_intermittent" );

 // Reactive Power Variable, if any - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( f_reactive_power ) {
  v_reactive_power.resize( f_time_horizon );
  for( auto & var : v_reactive_power )
   var.set_type( ColVariable::kNonNegative );
  add_static_variable( v_reactive_power , "q_intermittent" );
  }

 // Primary Spinning Reserve Variable - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( reserve_vars & 1u )  // if UCBlock has primary demand variables
  if( f_gamma != 0 ) {  // if unit produces any reserve
   v_primary_spinning_reserve.resize( f_time_horizon );
   for( auto & var : v_primary_spinning_reserve )
    var.set_type( ColVariable::kNonNegative );
   add_static_variable( v_primary_spinning_reserve , "pr_intermittent" );
  }

 // Secondary Spinning Reserve Variable - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( reserve_vars & 2u )  // if UCBlock has secondary demand variables
  if( f_gamma != 0 ) {  // if unit produces any reserve
   v_secondary_spinning_reserve.resize( f_time_horizon );
   for( auto & var : v_secondary_spinning_reserve )
    var.set_type( ColVariable::kNonNegative );
   add_static_variable( v_secondary_spinning_reserve , "sr_intermittent" );
  }

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 set_variables_generated();

 }  // end( IntermittentUnitBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::generate_abstract_constraints( Configuration * stcc )
{
 if( constraints_generated() )  // constraints have already been generated
  return;                       // nothing to do

 LinearFunction::v_coeff_pair vars;

 if( f_InvestmentCost == 0 ) {

  // Minimum power constraints

  min_power_Const.resize( f_time_horizon );

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {

   vars.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );

   if( f_gamma != 0 ) {  // if unit produces any reserve
    if( reserve_vars & 1u )  // if UCBlock has primary demand variables
     vars.push_back( std::make_pair( &v_primary_spinning_reserve[ t ] ,
                                     -1.0 ) );
    if( reserve_vars & 2u )  // if UCBlock has secondary demand variables
     vars.push_back( std::make_pair( &v_secondary_spinning_reserve[ t ] ,
                                     -1.0 ) );
   }

   min_power_Const[ t ].set_lhs( f_kappa * v_MinPower[ t ] );
   min_power_Const[ t ].set_rhs( Inf< double >() );
   min_power_Const[ t ].set_function(
    new LinearFunction( std::move( vars ) ) );
  }

  add_static_constraint( min_power_Const , "MinPower_Intermittent" );

  // Maximum power constraints

  if( f_gamma != 0 ) {  // if unit produces any reserve

   max_power_Const.resize( f_time_horizon );

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {

    vars.push_back( std::make_pair( &v_active_power[ t ] , f_gamma ) );

    if( reserve_vars & 1u )  // if UCBlock has primary demand variables
     vars.push_back( std::make_pair( &v_primary_spinning_reserve[ t ] ,
                                     1.0 ) );
    if( reserve_vars & 2u )  // if UCBlock has secondary demand variables
     vars.push_back( std::make_pair( &v_secondary_spinning_reserve[ t ] ,
                                     1.0 ) );

    max_power_Const[ t ].set_lhs( -Inf< double >() );
    max_power_Const[ t ].set_rhs( f_gamma * f_kappa * v_MaxPower[ t ] );
    max_power_Const[ t ].set_function(
     new LinearFunction( std::move( vars ) ) );
   }

   add_static_constraint( max_power_Const , "MaxPower_Intermittent" );
  }

  // Active power bounds constraints

  active_power_bounds_Const.resize( f_time_horizon );

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   active_power_bounds_Const[ t ].set_lhs( f_kappa * v_MinPower[ t ] );
   active_power_bounds_Const[ t ].set_rhs( f_kappa * v_MaxPower[ t ] );
   active_power_bounds_Const[ t ].set_variable( &v_active_power[ t ] );
  }

  add_static_constraint( active_power_bounds_Const ,
                         "ActivePower_Intermittent" );

 } else {

  // Active power bounds design constraints

  active_power_bounds_design_Const.resize(
   boost::multi_array< FRowConstraint , 2 >::extent_gen()
   [ 2 ][ f_time_horizon ] );  // 2 dims, i.e., the lower and upper bounds

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {

   // Lower bound of the active power design constraints:
   //
   //      v_MinPower x <= v_active_power
   // => 0 <= v_active_power - v_MinPower x

   vars.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );
   vars.push_back( std::make_pair( &design , -f_kappa * v_MinPower[ t ] ) );

   active_power_bounds_design_Const[ 0 ][ t ].set_lhs( 0.0 );
   active_power_bounds_design_Const[ 0 ][ t ].set_rhs( Inf< double >() );
   active_power_bounds_design_Const[ 0 ][ t ].set_function(
    new LinearFunction( std::move( vars ) ) );

   // Upper bound of the active power design constraints:
   //
   //      v_active_power <= v_MaxPower x
   // => v_active_power - v_MaxPower x <= 0

   vars.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );
   vars.push_back( std::make_pair( &design , -f_kappa * v_MaxPower[ t ] ) );

   active_power_bounds_design_Const[ 1 ][ t ].set_lhs( -Inf< double >() );
   active_power_bounds_design_Const[ 1 ][ t ].set_rhs( 0.0 );
   active_power_bounds_design_Const[ 1 ][ t ].set_function(
    new LinearFunction( std::move( vars ) ) );
  }

  add_static_constraint( active_power_bounds_design_Const ,
                         "ActivePower_Design_Intermittent" );

  // Design bounds

  const double lb = std::max( 0.0 , f_MinCapacityDesign );
  const bool is_binary = ( f_MaxCapacityDesign < 0.0 );

  const double ub = is_binary
                     ? 1.0 : ( std::abs( f_MaxCapacityDesign ) == 1.0
                      ? 1.0 : std::abs( f_MaxCapacityDesign ) );

  if( ( lb == 1.0 ) && ( ub == 1.0 ) )
   design.is_unitary( true , eNoMod );
  else {
   design_bound_Const.set_lhs( lb );
   design_bound_Const.set_rhs( ub );
   design_bound_Const.set_variable( &design );

   add_static_constraint( design_bound_Const , "DesignBound_Intermittent" );

   if( is_binary )
    design.is_integer( true , eNoMod );
  }
 }

 // reactive power bounds constraints (if any) - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( f_reactive_power && 
     ( ( ! v_MinReactivePower.empty() ) || ( ! v_MaxReactivePower.empty() ) )
     ) {
  if( ReactivePower_Bound_Const.empty() )
   ReactivePower_Bound_Const.resize( f_time_horizon );

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   ReactivePower_Bound_Const[ t ].set_rhs( get_max_reactive_power( t ) );
   ReactivePower_Bound_Const[ t ].set_lhs( get_min_reactive_power( t ) );
   ReactivePower_Bound_Const[ t ].set_variable( & v_reactive_power[ t ] );
   }

  add_static_constraint( ReactivePower_Bound_Const ,
			 "ReactivePowerBound_intermittent" );

  /*!! Link between active and reactive power
  Reactive_2_Active_Const.resize( f_time_horizon );

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   // Q(t) - P(t) <= 0
   auto lfunc = new LinearFunction();
   lfunc->add_variable( &v_active_power[ t ], -1.0 );
   lfunc->add_variable( &v_reactive_power[ t ], 1.0 );

   Reactive_2_Active_Const[ t ].set_lhs( -Inf< double >() );
   Reactive_2_Active_Const[ t ].set_rhs( 0.0 );
   Reactive_2_Active_Const[ t ].set_function( lfunc );
   }

  add_static_constraint( Reactive_2_Active_Const, "QandP_inter" );
  !!*/
  }

 // Maximum and Minimum generation constraints (if any)
 if( ( f_MaxGeneration < Inf< double >() ) || ( f_MinGeneration > -Inf< double >() ) ) {

  LinearFunction::v_coeff_pair vpair( f_time_horizon ,
    std::make_pair( nullptr , 1 ) );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   vpair[ t ].first = &v_active_power[ t ];
  
  MaxMinGeneration_Const.set_lhs( f_MinGeneration );
  MaxMinGeneration_Const.set_rhs( f_MaxGeneration );
  MaxMinGeneration_Const.set_function( new LinearFunction( std::move( vpair ) ) );
  
  add_static_constraint( MaxMinGeneration_Const , "MaxMinGeneration_intermittent" );
 }

 set_constraints_generated();

 }  // end( IntermittentUnitBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::generate_objective( Configuration * objc )
{
 if( objective_generated() )  // Objective has already been generated
  return;                     // nothing to do

 auto lf = new LinearFunction();

 if( f_InvestmentCost != 0 )
  lf->add_variable( &design , f_InvestmentCost );

 if( ! v_ActivePowerCost.empty() )
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   lf->add_variable( &v_active_power[ t ] ,
		     f_scale * v_ActivePowerCost[ t ] );

 objective.set_function( lf );
 objective.set_sense( Objective::eMin );

 // Set Block objective
 this->set_objective( &objective , eNoMod );

 set_objective_generated();

}  // end( IntermittentUnitBlock::generate_objective )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::set_BlockConfig( BlockConfig * newBC ,
                                             bool deleteold )
{
 UnitBlock::set_BlockConfig( newBC , deleteold );

 if( ! f_BlockConfig )
  return;

 if( auto config = dynamic_cast< SimpleConfiguration< double > * >
     ( f_BlockConfig->f_extra_Configuration ) )
  f_max_power_epsilon = config->f_value;

 }  // end( IntermittentUnitBlock::set_BlockConfig )

/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR CHECKING THE IntermittentUnitBlock -------------*/
/*--------------------------------------------------------------------------*/

bool IntermittentUnitBlock::is_feasible( bool useabstract ,
                                         Configuration * fsbc )
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
  // Variables
  && ColVariable::is_feasible( v_active_power , tol )
  && ColVariable::is_feasible( v_primary_spinning_reserve , tol )
  && ColVariable::is_feasible( v_secondary_spinning_reserve , tol )
  // Constraints
  && RowConstraint::is_feasible( min_power_Const , tol , rel_viol )
  && RowConstraint::is_feasible( design_bound_Const , tol , rel_viol )
  && RowConstraint::is_feasible( max_power_Const , tol , rel_viol )
  && RowConstraint::is_feasible( active_power_bounds_design_Const , tol , rel_viol )
  && RowConstraint::is_feasible( active_power_bounds_Const , tol , rel_viol ) );

}  // end( IntermittentUnitBlock::is_feasible )

/*--------------------------------------------------------------------------*/
/*--- METHODS FOR LOADING, PRINTING & SAVING THE IntermittentUnitBlock -----*/
/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::serialize( netCDF::NcGroup & group ) const
{
 UnitBlock::serialize( group );

 // Serialize scalar variables

 if( f_InvestmentCost != 0 ) {
  ::serialize( group , "InvestmentCost" , netCDF::NcDouble() ,
               f_InvestmentCost );
  if( f_MinCapacityDesign != 0 )
   ::serialize( group , "MinCapacityDesign" , netCDF::NcDouble() ,
                f_MinCapacityDesign );
  if( f_MaxCapacityDesign != 1 )
   ::serialize( group , "MaxCapacityDesign" , netCDF::NcDouble() ,
                f_MaxCapacityDesign );
 }

 ::serialize( group , "Gamma" , netCDF::NcDouble() , f_gamma );
 ::serialize( group , "Kappa" , netCDF::NcDouble() , f_kappa );

 ::serialize( group , "MaxGeneration" , netCDF::NcDouble() , f_MaxGeneration );
 ::serialize( group , "MinGeneration" , netCDF::NcDouble() , f_MinGeneration );

 // Serialize one-dimensional variables

 auto TimeHorizon = group.getDim( "TimeHorizon" );
 auto NumberIntervals = group.getDim( "NumberIntervals" );

 /* This lambda identifies the appropriate dimension for the given variable
  * (whose name is "var_name") and serializes the variable. The variable may
  * have any of the following dimensions: TimeHorizon, NumberIntervals,
  * 1. "allow_scalar_var" indicates whether the variable can be serialized as
  * a scalar variable (in which case the variable must have dimension 1). */
 auto serialize = [ &group , &TimeHorizon , &NumberIntervals ](
    const std::string & var_name , const std::vector< double > & data ,
    const netCDF::NcType & ncType = netCDF::NcDouble() ,
    bool allow_scalar_var = true ) {
  if( data.empty() )
   return;
  netCDF::NcDim dimension;
  if( data.size() == TimeHorizon.getSize() )
   dimension = TimeHorizon;
  else if( data.size() == NumberIntervals.getSize() )
   dimension = NumberIntervals;
  else if( data.size() != 1 )
   throw( std::logic_error(
    "IntermittentUnitBlock::serialize: invalid dimension for variable " +
    var_name + ": " + std::to_string( data.size() ) +
    ". Its dimension must be one of the following: TimeHorizon, "
    "NumberIntervals, 1" ) );

  ::serialize( group , var_name , ncType , dimension , data ,
               allow_scalar_var );
 };

 serialize( "MinPower" , v_MinPower );
 serialize( "MaxPower" , v_MaxPower );
 serialize( "InertiaPower" , v_InertiaPower );
 serialize( "ActivePowerCost" , v_ActivePowerCost );

 }  // end( IntermittentUnitBlock::serialize )

/*--------------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/

Solution * IntermittentUnitBlock::get_Solution( Configuration * csolc ,
						bool emptys )
{
 Index wsol = 15;
 if( ( ! csolc ) && f_BlockConfig )
  csolc = f_BlockConfig->f_solution_Configuration;

 if( auto config = dynamic_cast< SimpleConfiguration< int > * >( csolc ) )
  wsol = config->f_value;

 // call the method of the base class
 auto * sol = dynamic_cast< IntermittentUnitBlockSolution * >(
		                  UnitBlock::get_Solution( csolc , emptys ) );
 assert( sol );

 if( ! emptys )
  sol->read( this );

 return( sol );
 }

/*--------------------------------------------------------------------------*/

UnitBlockSolution * IntermittentUnitBlock::new_Solution( void ) const {
 return( new IntermittentUnitBlockSolution() );
 }

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::update_max_power_in_cnstrs( const Subset & time ,
                                                        c_ModParam issueAMod )
{
 if( ! max_power_Const.empty() )
  for( auto t : time )
   max_power_Const[ t ].set_rhs( f_kappa * f_gamma * v_MaxPower[ t ] ,
                                 issueAMod );

 if( ! active_power_bounds_Const.empty() )
  for( auto t : time )
   active_power_bounds_Const[ t ].set_rhs( f_kappa * v_MaxPower[ t ] ,
                                           issueAMod );
 // FIXME: use a GroupModification
 }  // end( IntermittentUnitBlock::update_max_power_in_cnstrs ( subset ) )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::update_max_power_in_cnstrs( const Range & time ,
                                                        c_ModParam issueAMod )
{
 if( ! max_power_Const.empty() )
  for( auto t = time.first ; t < time.second ; ++t )
   max_power_Const[ t ].set_rhs( f_kappa * f_gamma * v_MaxPower[ t ] ,
                                 issueAMod );
 // FIXME: use a GroupModification
 if( ! active_power_bounds_Const.empty() )
  for( auto t = time.first ; t < time.second ; ++t )
   active_power_bounds_Const[ t ].set_rhs( f_kappa * v_MaxPower[ t ] ,
                                           issueAMod );

 }  // end( IntermittentUnitBlock::update_max_power_in_cnstrs ( range ) )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::set_maximum_power( MF_dbl_it values ,
                                               Subset && subset ,
                                               bool ordered ,
                                               c_ModParam issuePMod ,
                                               c_ModParam issueAMod )
{
 if( subset.empty() )
  return;

 if( v_MaxPower.empty() ) {
  if( std::all_of( values ,
                   values + subset.size() ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  v_MaxPower.assign( f_time_horizon , 0.0 );
 }

 for( auto t : subset )
  if( t >= v_MaxPower.size() )
   throw( std::invalid_argument(
    "IntermittentUnitBlock::set_maximum_power: invalid value in subset." ) );

 auto values_it = values;

 bool identical = true;
 for( auto t : subset ) {
  double max_power = *( values_it++ );
  if( ( f_max_power_epsilon > 0 ) && ( max_power == 0.0 ) )
   max_power = f_max_power_epsilon;

  if( v_MaxPower[ t ] != max_power ) {
   identical = false;
   break;
  }
 }

 if( identical )
  return;

 if( not_dry_run( issuePMod ) ) {
  values_it = values;
  for( auto t : subset ) {
   double max_power = *( values_it++ );
   if( ( f_max_power_epsilon > 0 ) && ( max_power == 0.0 ) )
    max_power = f_max_power_epsilon;

   v_MaxPower[ t ] = max_power;
  }

  if( not_dry_run( issueAMod ) && constraints_generated() )
   update_max_power_in_cnstrs( subset , issueAMod );
 }

 if( issue_pmod( issuePMod ) ) {
  if( ! ordered )
   std::sort( subset.begin() , subset.end() );

  Block::add_Modification(
   std::make_shared< IntermittentUnitBlockSbstMod >(
    this , IntermittentUnitBlockMod::eSetMaxP , std::move( subset ) ) ,
   Observer::par2chnl( issuePMod ) );
 }

}  // end( IntermittentUnitBlock::set_maximum_power( subset ) )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::set_maximum_power( MF_dbl_it values ,
                                               Range rng ,
                                               c_ModParam issuePMod ,
                                               c_ModParam issueAMod )
{
 rng.second = std::min( rng.second , f_time_horizon );
 if( rng.second <= rng.first )
  return;

 c_Index sz = rng.second - rng.first;

 if( v_MaxPower.empty() ) {
  if( std::all_of( values ,
                   values + sz ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  v_MaxPower.assign( f_time_horizon , 0.0 );
 }

 bool identical = true;
 auto values_it = values;
 for( Index t = rng.first ; t < rng.second ; ++t ) {
  double max_power = *( values_it++ );
  if( ( f_max_power_epsilon > 0 ) && ( max_power == 0.0 ) )
   max_power = f_max_power_epsilon;

  if( v_MaxPower[ t ] != max_power ) {
   identical = false;
   break;
  }
 }

 if( identical )
  return;

 if( not_dry_run( issuePMod ) ) {
  values_it = values;
  for( Index t = rng.first ; t < rng.second ; ++t ) {
   double max_power = *( values_it++ );
   if( ( f_max_power_epsilon > 0 ) && ( max_power == 0.0 ) )
    max_power = f_max_power_epsilon;

   v_MaxPower[ t ] = max_power;
  }

  if( not_dry_run( issueAMod ) && constraints_generated() )
   update_max_power_in_cnstrs( rng , issueAMod );
 }

 if( issue_pmod( issuePMod ) )
  Block::add_Modification(
   std::make_shared< IntermittentUnitBlockRngdMod >(
    this , IntermittentUnitBlockMod::eSetMaxP , rng ) ,
   Observer::par2chnl( issuePMod ) );

}  // end( IntermittentUnitBlock::set_maximum_power( range ) )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::set_active_power_cost( MF_dbl_it values ,
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
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  v_ActivePowerCost.assign( f_time_horizon , 0.0 );
 }

 for( auto t : subset )
  if( t >= v_ActivePowerCost.size() )
   throw( std::invalid_argument(
    "IntermittentUnitBlock::set_active_power_cost: invalid index in subset." ) );

 auto values_it = values;

 bool identical = true;
 for( auto t : subset ) {
  if( v_ActivePowerCost[ t ] != *( values_it++ ) ) {
   identical = false;
   break;
  }
 }

 if( identical )
  return;

 if( not_dry_run( issuePMod ) ) {
  values_it = values;
  for( auto t : subset )
   v_ActivePowerCost[ t ] = *( values_it++ );

  if( not_dry_run( issueAMod ) && objective_generated() ) {
   auto * lf = static_cast< LinearFunction * >( objective.get_function() );

   for( auto t : subset ) {
    const auto idx = lf->is_active( &v_active_power[ t ] );

    if( idx == Inf< Index >() )
     throw( std::logic_error(
      "IntermittentUnitBlock::set_active_power_cost: expected Variable not "
      "found in objective." ) );

    lf->modify_coefficient( idx ,
                            f_scale * v_ActivePowerCost[ t ] ,
                            issueAMod );
   }
  }
 }

 if( issue_pmod( issuePMod ) ) {
  if( ! ordered )
   std::sort( subset.begin() , subset.end() );

  Block::add_Modification(
   std::make_shared< IntermittentUnitBlockSbstMod >(
    this , IntermittentUnitBlockMod::eSetActPCost , std::move( subset ) ) ,
   Observer::par2chnl( issuePMod ) );
 }
}  // end( IntermittentUnitBlock::set_active_power_cost( subset ) )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::set_active_power_cost( MF_dbl_it values ,
                                                   Range rng ,
                                                   c_ModParam issuePMod ,
                                                   c_ModParam issueAMod )
{
 rng.second = std::min( rng.second , f_time_horizon );
 if( rng.second <= rng.first )
  return;

 c_Index sz = rng.second - rng.first;

 if( v_ActivePowerCost.empty() ) {
  if( std::all_of( values ,
                   values + sz ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  v_ActivePowerCost.assign( f_time_horizon , 0.0 );
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

   for( Index t = rng.first ; t < rng.second ; ++t ) {
    const auto idx = lf->is_active( &v_active_power[ t ] );

    if( idx == Inf< Index >() )
     throw( std::logic_error(
      "IntermittentUnitBlock::set_active_power_cost: expected Variable not "
      "found in objective." ) );

    lf->modify_coefficient( idx ,
                            f_scale * v_ActivePowerCost[ t ] ,
                            issueAMod );
   }
  }
 }

 if( issue_pmod( issuePMod ) )
  Block::add_Modification(
   std::make_shared< IntermittentUnitBlockRngdMod >(
    this , IntermittentUnitBlockMod::eSetActPCost , rng ) ,
   Observer::par2chnl( issuePMod ) );

}  // end( IntermittentUnitBlock::set_active_power_cost( range ) )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::scale( MF_dbl_it values ,
                                   Subset && subset , bool ordered ,
                                   c_ModParam issuePMod ,
                                   c_ModParam issueAMod )
{
 if( subset.empty() )
  return;  // Since the given Subset is empty, no operation is performed

 if( f_scale == *values )  // the scale factor does not change: nothing to do
  return;

 if( not_dry_run( issuePMod ) )
  f_scale = *values;  // Update the scale factor

 if( issue_pmod( issuePMod ) )  // issue a Physical Modification
  Block::add_Modification( std::make_shared< UnitBlockMod >(
                                           this , UnitBlockMod::eScale ) ,
                                           Observer::par2chnl( issuePMod ) );

 }  // end( IntermittentUnitBlock::scale )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::set_kappa( MF_dbl_it values ,
                                       Subset && subset , bool ordered ,
                                       c_ModParam issuePMod ,
                                       c_ModParam issueAMod )
{
 if( subset.empty() )
  return;  // Since the given Subset is empty, no operation is performed

 if( f_kappa == *values )
  return;  // The kappa constant does not change: nothing to do

 if( not_dry_run( issuePMod ) ) {
  f_kappa = *values;  // Update the kappa constant

  if( not_dry_run( issueAMod ) ) {
   // Update the abstract representation
   if( constraints_generated() ) {
    // Update the constraints

    if( ! active_power_bounds_Const.empty() )

     for( Index t = 0 ; t < f_time_horizon ; ++t ) {
      active_power_bounds_Const[ t ].set_lhs(
       f_kappa * v_MinPower[ t ] , issueAMod );
      active_power_bounds_Const[ t ].set_rhs(
       f_kappa * v_MaxPower[ t ] , issueAMod );
     }

    else if( ! active_power_bounds_design_Const.empty() )

     for( Index t = 0 ; t < f_time_horizon ; ++t ) {
      auto f0 = static_cast< LinearFunction * >(
       active_power_bounds_design_Const[ 0 ][ t ].get_function() );

      const auto design_idx0 = f0->is_active( &design );

      if( design_idx0 == Inf< Index >() )
       throw( std::logic_error( "IntermittentUnitBlock::set_kappa: expected "
                                "Variable not found in "
                                "active_power_bounds_design_Const." ) );

      f0->modify_coefficient( design_idx0 ,
                              -f_kappa * v_MinPower[ t ] ,
                              issueAMod );

      auto f1 = static_cast< LinearFunction * >(
       active_power_bounds_design_Const[ 1 ][ t ].get_function() );

      const auto design_idx1 = f1->is_active( &design );

      if( design_idx1 == Inf< Index >() )
       throw( std::logic_error( "IntermittentUnitBlock::set_kappa: expected "
                                "Variable not found in "
                                "active_power_bounds_design_Const." ) );

      f1->modify_coefficient( design_idx1 ,
                              -f_kappa * v_MaxPower[ t ] ,
                              issueAMod );
     }

    if( ! min_power_Const.empty() )
     for( Index t = 0 ; t < f_time_horizon ; ++t )
      min_power_Const[ t ].set_lhs( f_kappa * v_MinPower[ t ] ,
                                    issueAMod );

    if( ! max_power_Const.empty() )
     for( Index t = 0 ; t < f_time_horizon ; ++t )
      max_power_Const[ t ].set_rhs( f_gamma * f_kappa * v_MaxPower[ t ] ,
                                    issueAMod );
   }  // end( constraints_generated )
  }  // end( if( not_dry_run( issueAMod ) )
 }  // end( if( not_dry_run( issuePMod ) )

 if( issue_pmod( issuePMod ) )  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< IntermittentUnitBlockMod >(
                            this , IntermittentUnitBlockMod::eSetKappa ) ,
                           Observer::par2chnl( issuePMod ) );

 }  // end( IntermittentUnitBlock::set_kappa( subset ) )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::set_kappa( MF_dbl_it values , Range rng ,
                                       c_ModParam issuePMod ,
                                       c_ModParam issueAMod )
{
 if( rng.first >= rng.second )
  return;  // An empty Range was given: no operation is performed

 Subset subset( 1 , rng.first );

 set_kappa( values , std::move( subset ) , true , issuePMod , issueAMod );

}  // end( IntermittentUnitBlock::set_kappa( range ) )

/*--------------------------------------------------------------------------*/
/*--------------- METHODS OF IntermittentUnitBlockSolution -----------------*/
/*--------------------------------------------------------------------------*/

void IntermittentUnitBlockSolution::deserialize(
					      const netCDF::NcGroup & group )
{
 // call the method of the base class
 UnitBlockSolution::deserialize( group );

 if( f_number_generators != 1 )
  throw( std::logic_error( "IntermittentUnitBlockSolution::deserialize: "
			   "intermittents have only one generator" ) );

 // deserialize the design - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! ::deserialize< double >( group , f_design , "IntermittentDesign" ) )
  f_design = dNaN;

 }  // end( IntermittentUnitBlockSolution::deserialize )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlockSolution::read( const Block * block )
{
 auto IUB = dynamic_cast< const IntermittentUnitBlock * >( block );
 if( ! IUB )
  throw( std::invalid_argument( "IntermittentUnitBlockSolution::read: block "
				"is not a IntermittentUnitBlock" ) );

 UnitBlockSolution::read( IUB );  // call the method of the base class

 // read the design- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 f_design = IUB->get_const_design().get_value();

 }  // end( IntermittentUnitBlockSolution::read )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlockSolution::write( Block * block )
{
 UnitBlockSolution::write( block );  // call the method of the base class

 auto IUB = dynamic_cast< IntermittentUnitBlock * >( block );
 if( ! IUB )
  throw( std::invalid_argument( "IntermittentUnitBlockSolution::read: block "
				"is not a IntermittentUnitBlock" ) );

 // write the design - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 IUB->get_design().set_value( f_design );

 }  // end( IntermittentUnitBlockSolution::write )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlockSolution::serialize( netCDF::NcGroup & group ) const
{
 UnitBlockSolution::serialize( group );  // call the method of the base class

 // serialize the design- - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! std::isnan( f_design ) )
   ::serialize< double >( group , "IntermittentDesign" , netCDF::NcDouble() ,
			  f_design );

 }  // end( IntermittentUnitBlockSolution::serialize )

/*--------------------------------------------------------------------------*/

IntermittentUnitBlockSolution * IntermittentUnitBlockSolution::scale(
						        double factor ) const
{
 auto sol = clone();

 if( factor == 1 )
  return( sol );

 guts_of_scale( sol , factor );

 if( ! std::isnan( f_design ) )
  sol->f_design *= factor;

 return( sol );

 }  // end( IntermittentUnitBlockSolution::scale )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlockSolution::sum( const Solution * solution ,
					 double multiplier )
{
 // call the method of the base class
 UnitBlockSolution::sum( solution , multiplier );

 auto IUBS = dynamic_cast< const IntermittentUnitBlockSolution * >(
								 solution );
 if( ! IUBS )
  throw( std::invalid_argument( "IntermittentUnitBlockSolution::sum: "
				"solution not a "
				"IntermittentUnitBlockSolution" ) );

 if( ! std::isnan( f_design ) )
  f_design += IUBS->f_design * multiplier;

 }  // end( IntermittentUnitBlockSolution::sum )

/*--------------------------------------------------------------------------*/

IntermittentUnitBlockSolution * IntermittentUnitBlockSolution::clone(
							   bool empty ) const
{
 auto * sol = new IntermittentUnitBlockSolution();

 if( ! empty ) {
  guts_of_clone( sol );
  sol->f_design = f_design;
  }

 return( sol );

 }  // end( IntermittentUnitBlockSolution::clone )

/*--------------------------------------------------------------------------*/
/*------------------- End File IntermittentUnitBlock.cpp -------------------*/
/*--------------------------------------------------------------------------*/
