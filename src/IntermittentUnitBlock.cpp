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

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS OF IntermittentUnitBlock -------------------*/
/*--------------------------------------------------------------------------*/

IntermittentUnitBlock::~IntermittentUnitBlock()
{
 Constraint::clear( min_power_Const );
 Constraint::clear( max_power_Const );
 Constraint::clear( active_power_bounds_design_Const );

 Constraint::clear( active_power_bounds_Const );

 objective.clear();
}

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::deserialize( const netCDF::NcGroup & group )
{

#ifndef NDEBUG
 static std::vector< std::string > expected_dims = { "TimeHorizon" ,
                                                     "NumberIntervals" };
 check_dimensions( group , expected_dims , std::cerr );

 static std::vector< std::string > expected_vars = { "InvestmentCost" ,
                                                     "MaxCapacity" ,
                                                     "MinPower" , "MaxPower" ,
                                                     "InertiaPower" ,
                                                     "Gamma" , "Kappa" };
 check_variables( group , expected_vars , std::cerr );
#endif

 // Deserialize data from the base class
 UnitBlock::deserialize( group );

 // Mandatory variables

 ::deserialize( group , "MaxPower" , v_MaxPower , false );

 // Optional variables

 ::deserialize( group , f_InvestmentCost , "InvestmentCost" );

 ::deserialize( group , f_MaxCapacity , "MaxCapacity" );

 if( ! ::deserialize( group , "MinPower" , v_MinPower ) )
  v_MinPower.resize( f_time_horizon );

 if( ! ::deserialize( group , "InertiaPower" , v_InertiaPower ) )
  v_InertiaPower.resize( f_time_horizon );

 ::deserialize( group , f_gamma , "Gamma" );

 ::deserialize( group , f_kappa , "Kappa" );

 // Decompress vectors

 decompress_vector( v_MinPower );
 decompress_vector( v_MaxPower );
 decompress_vector( v_InertiaPower );

 if( f_max_power_epsilon > 0 )
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   if( v_MaxPower[ t ] == 0.0 )
    v_MaxPower[ t ] = f_max_power_epsilon;

 check_data_consistency();

}  // end( IntermittentUnitBlock::deserialize )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::check_data_consistency( void ) const
{
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

  if( v_MinPower[ t ] < 0 )
   throw( std::logic_error( "IntermittentUnitBlock::check_data_consistency: "
                             "minimum power at time " + std::to_string( t ) +
                             " is " + std::to_string( v_MinPower[ t ] ) +
                             ", which is negative." ) );
 }

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
 if( f_InvestmentCost != 0 ) {
  design.set_type( ColVariable::kPosUnitary );
  add_static_variable( design , "x_intermittent" );
 }

 // Active Power Variable
 v_active_power.resize( f_time_horizon );
 for( auto & var : v_active_power )
  var.set_type( ColVariable::kNonNegative );
 add_static_variable( v_active_power , "p_intermittent" );

 // Primary Spinning Reserve Variable
 if( reserve_vars & 1u )  // if UCBlock has primary demand variables
  if( f_gamma != 0 ) {  // if unit produces any reserve
   v_primary_spinning_reserve.resize( f_time_horizon );
   for( auto & var : v_primary_spinning_reserve )
    var.set_type( ColVariable::kNonNegative );
   add_static_variable( v_primary_spinning_reserve , "pr_intermittent" );
  }

 // Secondary Spinning Reserve Variable
 if( reserve_vars & 2u )  // if UCBlock has secondary demand variables
  if( f_gamma != 0 ) {  // if unit produces any reserve
   v_secondary_spinning_reserve.resize( f_time_horizon );
   for( auto & var : v_secondary_spinning_reserve )
    var.set_type( ColVariable::kNonNegative );
   add_static_variable( v_secondary_spinning_reserve , "sr_intermittent" );
  }

 set_variables_generated();

}  // end( IntermittentUnitBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::generate_abstract_constraints( Configuration * stcc )
{
 if( constraints_generated() )  // constraints have already been generated
  return;                       // nothing to do

 LinearFunction::v_coeff_pair vars;

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
  min_power_Const[ t ].set_function( new LinearFunction( std::move( vars ) ) );
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
   max_power_Const[ t ].set_function( new LinearFunction( std::move( vars ) ) );
  }

  add_static_constraint( max_power_Const , "MaxPower_Intermittent" );
 }

 if( f_InvestmentCost == 0 ) {

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
   //      v_MinPower x <= v_active_power     x \in {0,1}, for all t
   // => 0 <= v_active_power - v_MinPower x   x \in {0,1}, for all t

   vars.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );
   vars.push_back( std::make_pair( &design , -f_kappa * v_MinPower[ t ] ) );

   active_power_bounds_design_Const[ 0 ][ t ].set_lhs( 0.0 );
   active_power_bounds_design_Const[ 0 ][ t ].set_rhs( Inf< double >() );
   active_power_bounds_design_Const[ 0 ][ t ].set_function(
    new LinearFunction( std::move( vars ) ) );

   // Upper bound of the active power design constraints:
   //
   //      v_active_power <= v_MaxPower x     x \in {0,1}, for all t
   // => v_active_power - v_MaxPower x <= 0   x \in {0,1}, for all t

   vars.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );
   vars.push_back( std::make_pair( &design , -f_kappa * v_MaxPower[ t ] ) );

   active_power_bounds_design_Const[ 1 ][ t ].set_lhs( -Inf< double >() );
   active_power_bounds_design_Const[ 1 ][ t ].set_rhs( 0.0 );
   active_power_bounds_design_Const[ 1 ][ t ].set_function(
    new LinearFunction( std::move( vars ) ) );
  }

  add_static_constraint( active_power_bounds_design_Const ,
                         "ActivePower_Design_Intermittent" );
 }

 set_constraints_generated();

}  // end( IntermittentUnitBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::generate_objective( Configuration * objc )
{
 if( objective_generated() )  // Objective has already been generated
  return;                     // nothing to do

 LinearFunction::v_coeff_pair vars;

 if( f_InvestmentCost != 0 )
  vars.push_back( std::make_pair( &design , f_InvestmentCost ) );

 objective.set_function( new LinearFunction( std::move( vars ) ) );
 objective.set_sense( Objective::eMin );

 // Set Block objective
 this->set_objective( &objective );

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

} // end( IntermittentUnitBlock::set_BlockConfig )

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

 if( f_InvestmentCost != 0 )
  ::serialize( group , "InvestmentCost" , netCDF::NcDouble() ,
               f_InvestmentCost );

 if( f_MaxCapacity != 0 )
  ::serialize( group , "MaxCapacity" , netCDF::NcDouble() , f_MaxCapacity );

 ::serialize( group , "Gamma" , netCDF::NcDouble() , f_gamma );
 ::serialize( group , "Kappa" , netCDF::NcDouble() , f_kappa );

 // Serialize one-dimensional variables

 auto TimeHorizon = group.getDim( "TimeHorizon" );
 auto NumberIntervals = group.getDim( "NumberIntervals" );

 /* This lambda identifies the appropriate dimension for the given variable
  * (whose name is "var_name") and serializes the variable. The variable may
  * have any of the following dimensions: TimeHorizon, NumberIntervals,
  * 1. "allow_scalar_var" indicates whether the variable can be serialized as
  * a scalar variable (in which case the variable must have dimension 1). */
 auto serialize = [ &group , &TimeHorizon , &NumberIntervals ]
  ( const std::string & var_name , const std::vector< double > & data ,
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
    "NumberIntervals, 1." ) );

  ::serialize( group , var_name , ncType , dimension , data ,
               allow_scalar_var );
 };

 serialize( "MinPower" , v_MinPower );
 serialize( "MaxPower" , v_MaxPower );
 serialize( "InertiaPower" , v_InertiaPower );

}  // end( IntermittentUnitBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::update_max_power_in_cnstrs( const Subset & time ,
                                                        ModParam issueAMod )
{
 if( ! max_power_Const.empty() )
  for( auto t : time )
   max_power_Const[ t ].set_rhs( f_kappa * f_gamma * v_MaxPower[ t ] ,
                                 issueAMod );
   // FIXME: use a GroupModification

 if( ! active_power_bounds_Const.empty() )
  for( auto t : time )
   active_power_bounds_Const[ t ].set_rhs( f_kappa * v_MaxPower[ t ] ,
                                           issueAMod );
   // FIXME: use a GroupModification
}  // end( IntermittentUnitBlock::update_max_power_in_cnstrs ( subset ) )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::update_max_power_in_cnstrs( const Range & time ,
                                                        ModParam issueAMod )
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
   // FIXME: use a GroupModification
}  // end( IntermittentUnitBlock::update_max_power_in_cnstrs ( range ) )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::set_maximum_power( MF_dbl_it values ,
                                               Subset && subset ,
                                               const bool ordered ,
                                               ModParam issuePMod ,
                                               ModParam issueAMod )
{
 if( subset.empty() )
  return;

 if( v_MaxPower.empty() ) {
  if( std::all_of( values ,
                   values + subset.size() ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  Index max_index = *std::max_element( std::begin( subset ) ,
                                       std::end( subset ) );
  v_MaxPower.resize( max_index );
 }

 // If nothing changes, return
 bool identical = true;
 for( auto t : subset ) {
  if( t >= v_MaxPower.size() )
   throw( std::invalid_argument( "IntermittentUnitBlock::set_maximum_power:"
                                 " invalid value in subset." ) );
  auto max_power = *(values++);
  if( v_MaxPower[ t ] != max_power ) {
   identical = false;
   if( not_dry_run( issuePMod ) )
    // Change the physical representation
    v_MaxPower[ t ] = max_power;

   if( ( f_max_power_epsilon > 0 ) && ( v_MaxPower[ t ] == 0.0 ) )
    v_MaxPower[ t ] = f_max_power_epsilon;
  }
 }
 if( identical )
  return;  // nothing changes; return

 if( not_dry_run( issuePMod ) &&
     not_dry_run( issueAMod ) &&
     constraints_generated() )
  // Change the abstract representation
  update_max_power_in_cnstrs( subset , issueAMod );

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification
  if( ! ordered )
   std::sort( subset.begin() , subset.end() );

  Block::add_Modification( std::make_shared< IntermittentUnitBlockSbstMod >(
                            this , IntermittentUnitBlockMod::eSetMaxP ,
                            std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );
 }
}  // end( IntermittentUnitBlock::set_maximum_power( subset ) )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::set_maximum_power( MF_dbl_it values ,
                                               Range rng ,
                                               ModParam issuePMod ,
                                               ModParam issueAMod )
{
 rng.second = std::min( rng.second , f_time_horizon );
 if( rng.second <= rng.first )
  return;

 if( v_MaxPower.empty() ) {
  if( std::all_of( values ,
                   values + ( rng.second - rng.first ) ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  auto max_index = rng.second;
  v_MaxPower.resize( max_index );
 }

 // If nothing changes, return
 if( std::equal( values , values + ( rng.second - rng.first ) ,
                 v_MaxPower.begin() + rng.first ) )
  return;

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation
  std::copy( values , values + ( rng.second - rng.first ) ,
             v_MaxPower.begin() + rng.first );

  if( f_max_power_epsilon > 0 )
   for( Index t = rng.first ; t < rng.second ; ++t )
    if( v_MaxPower[ t ] == 0.0 )
     v_MaxPower[ t ] = f_max_power_epsilon;

  if( not_dry_run( issueAMod ) && constraints_generated() )
   // Change the abstract representation
   update_max_power_in_cnstrs( rng , issueAMod );
 }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< IntermittentUnitBlockRngdMod >(
                            this , IntermittentUnitBlockMod::eSetMaxP , rng ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( IntermittentUnitBlock::set_maximum_power( range ) )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::scale( MF_dbl_it values ,
                                   Subset && subset ,
                                   const bool ordered ,
                                   c_ModParam issuePMod ,
                                   c_ModParam issueAMod )
{
 if( subset.empty() )
  return;  // Since the given Subset is empty, no operation is performed

 if( f_scale == *values )
  return;  // The scale factor does not change: nothing to do

 if( not_dry_run( issuePMod ) )
  f_scale = *values;  // Update the scale factor

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< UnitBlockMod >(
                            this , UnitBlockMod::eScale ) ,
                           Observer::par2chnl( issuePMod ) );
 else if( auto f_Block = get_f_Block() )
  f_Block->add_Modification( std::make_shared< UnitBlockMod >(
                              this , UnitBlockMod::eScale ) ,
                             Observer::par2chnl( issuePMod ) );

}  // end( IntermittentUnitBlock::scale )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::set_kappa( MF_dbl_it values ,
                                       Subset && subset ,
                                       const bool ordered ,
                                       ModParam issuePMod ,
                                       ModParam issueAMod )
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
       throw( std::logic_error("IntermittentUnitBlock::set_kappa: expected "
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

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< IntermittentUnitBlockMod >(
                            this , IntermittentUnitBlockMod::eSetKappa ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( IntermittentUnitBlock::set_kappa( subset ) )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::set_kappa( MF_dbl_it values ,
                                       Range rng ,
                                       ModParam issuePMod , ModParam issueAMod )
{
 if( rng.first >= rng.second )
  return;  // An empty Range was given: no operation is performed

 Subset subset( 1 , 0 );

 set_kappa( values , std::move( subset ) , true , issuePMod , issueAMod );

}  // end( IntermittentUnitBlock::set_kappa( range ) )

/*--------------------------------------------------------------------------*/
/*------------------- End File IntermittentUnitBlock.cpp -------------------*/
/*--------------------------------------------------------------------------*/
