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
 * \copyright &copy; by Antonio Frangioni, Ali Ghezelsoflu,
 *                   Rafael Durbano Lobato
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

IntermittentUnitBlock::~IntermittentUnitBlock() {

 Constraint::clear( MinPower_Const );
 Constraint::clear( MaxPower_Const );

 Constraint::clear( active_power_bounds_Const );

 Constraint::clear( active_power_bounds_design_Const );

 objective.clear();
}

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::deserialize( const netCDF::NcGroup & group ) {

#ifndef NDEBUG
 static std::vector< std::string > expected_dims = { "TimeHorizon" ,
                                                     "NumberIntervals" };
 check_dimensions( group , expected_dims , std::cerr );

 static std::vector< std::string > expected_vars =
  { "MinPower" , "MaxPower" , "InertiaPower" , "Gamma" , "Kappa" ,
    "InvestmentCost" };

 check_variables( group , expected_vars , std::cerr );
#endif

 // Deserialize data from the base class
 UnitBlock::deserialize( group );

 // Mandatory variables

 ::deserialize( group , "MaxPower" , v_maximum_power , false );

 // Optional variables

 if( ! ::deserialize( group , "MinPower" , v_minimum_power ) )
  v_minimum_power.assign( f_time_horizon , 0 );

 if( ! ::deserialize( group , "InertiaPower" , v_inertia_power ) )
  v_inertia_power.assign( f_time_horizon , 0 );

 ::deserialize( group , f_gamma , "Gamma" );

 ::deserialize( group , f_kappa , "Kappa" );

 ::deserialize( group , f_investment_cost , "InvestmentCost" );

 // Decompress vectors
 decompress_vector( v_minimum_power );
 decompress_vector( v_maximum_power );
 decompress_vector( v_inertia_power );

 check_data_consistency();

}  // end( IntermittentUnitBlock::deserialize )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::check_data_consistency( void ) const {
 // Minimum and maximum power

 assert( v_minimum_power.size() == f_time_horizon );
 assert( v_maximum_power.size() == f_time_horizon );

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {
  if( v_minimum_power[ t ] > v_maximum_power[ t ] ) {
   throw( std::logic_error( "IntermittentUnitBlock::check_data_consistency: "
                             "minimum power at time " + std::to_string( t ) +
                             " is " + std::to_string( v_minimum_power[ t ] ) +
                             ", which is greater than the maximum power, which "
                             "is " + std::to_string( v_maximum_power[ t ] ) +
                             "." ) );
  }

  if( v_minimum_power[ t ] < 0 ) {
   throw( std::logic_error( "IntermittentUnitBlock::check_data_consistency: "
                             "minimum power at time " + std::to_string( t ) +
                             " is " + std::to_string( v_minimum_power[ t ] ) +
                             ", which is negative." ) );
  }
 }

 // Gamma

 if( ( f_gamma < 0 ) || ( f_gamma > 1 ) ) {
  throw( std::logic_error( "IntermittentUnitBlock::check_data_consistency: "
                            "gamma must be between 0 and 1, but it is " +
                            std::to_string( f_gamma ) + "." ) );
 }

 // Kappa

 if( f_kappa < 0 ) {
  throw( std::logic_error( "IntermittentUnitBlock::check_data_consistency: "
                            "kappa must be nonnegative, but it is" +
                            std::to_string( f_kappa ) + "." ) );
 }

 if( ! v_inertia_power.empty() ) {
  assert( v_inertia_power.size() == f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   if( v_inertia_power[ t ] < 0 ) {
    throw( std::logic_error( "IntermittentUnitBlock::check_data_consistency: "
                              "inertia power for time " + std::to_string( t ) +
                              " must be nonnegative, but it is" +
                              std::to_string( v_inertia_power[ t ] ) + "." ) );
   }
  }
 }
}  // end( IntermittentUnitBlock::check_data_consistency )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::generate_abstract_variables(
 Configuration * stvv ) {

 if( variables_generated() )
  return; // variables have already been generated

 UnitBlock::generate_abstract_variables( stvv );

 int relax_binary = 0;
 auto config = dynamic_cast<SimpleConfiguration< int > *>( stvv );
 if( ( !config ) && f_BlockConfig &&
     f_BlockConfig->f_static_variables_Configuration )
  config = dynamic_cast< SimpleConfiguration< int > * >
  ( f_BlockConfig->f_static_variables_Configuration );
 if( config )
  relax_binary = config->f_value;

 // Design Variable
 if( f_investment_cost != 0 ) {
  if( relax_binary )
   v_design.set_type( ColVariable::kPosUnitary );
  else
   v_design.set_type( ColVariable::kBinary );
  add_static_variable( v_design , "D_intermittent" );
 }

 // Active Power Variable
 v_active_power.resize( f_time_horizon );
 for( auto & var : v_active_power )
  var.set_type( ColVariable::kNonNegative );
 add_static_variable( v_active_power , "p_intermittent" );

 // Primary Spinning Reserve Variable
 if( reserve_vars & 1u ) { // if UCBlock has primary demand variables
  if( f_gamma != 0 ) { // if unit produces any reserve
   v_primary_spinning_reserve.resize( f_time_horizon );
   for( auto & var : v_primary_spinning_reserve )
    var.set_type( ColVariable::kNonNegative );
   add_static_variable( v_primary_spinning_reserve , "pr_intermittent" );
  }
 }

 // Secondary Spinning Reserve Variable
 if( reserve_vars & 2u ) { // if UCBlock has secondary demand variables
  if( f_gamma != 0 ) { // if unit produces any reserve
   v_secondary_spinning_reserve.resize( f_time_horizon );
   for( auto & var : v_secondary_spinning_reserve )
    var.set_type( ColVariable::kNonNegative );
   add_static_variable( v_secondary_spinning_reserve , "sr_intermittent" );
  }
 }

 set_variables_generated();

}  // end( IntermittentUnitBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::generate_abstract_constraints(
 Configuration * stcc ) {

 if( constraints_generated() )
  return; // constraints have already been generated

 // Maximum power constraints

 if( f_gamma != 0 ) { // if unit produces any reserve

  MaxPower_Const.resize( f_time_horizon );

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {

   LinearFunction::v_coeff_pair vars;

   vars.push_back( std::make_pair( &v_active_power[ t ] , f_gamma ) );

   if( reserve_vars & 1u )
    vars.push_back( std::make_pair( &v_primary_spinning_reserve[ t ] ,
                                    1.0 ) );
   if( reserve_vars & 2u )
    vars.push_back( std::make_pair( &v_secondary_spinning_reserve[ t ] ,
                                    1.0 ) );

   MaxPower_Const[ t ].set_lhs( -Inf< double >() );
   MaxPower_Const[ t ].set_rhs( f_gamma * f_kappa * v_maximum_power[ t ] );
   MaxPower_Const[ t ].set_function( new LinearFunction( std::move( vars ) ) );
  }

  add_static_constraint( MaxPower_Const , "MaxPower_Intermittent" );
 }

 // Minimum power constraints

 MinPower_Const.resize( f_time_horizon );

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {

  LinearFunction::v_coeff_pair vars;

  vars.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );

  if( f_gamma != 0 ) { // if unit produces any reserve
   if( reserve_vars & 1u )
    vars.push_back( std::make_pair( &v_primary_spinning_reserve[ t ] ,
                                    -1.0 ) );
   if( reserve_vars & 2u )
    vars.push_back( std::make_pair( &v_secondary_spinning_reserve[ t ] ,
                                    -1.0 ) );
  }

  MinPower_Const[ t ].set_lhs( f_kappa * v_minimum_power[ t ] );
  MinPower_Const[ t ].set_rhs( Inf< double >() );
  MinPower_Const[ t ].set_function( new LinearFunction( std::move( vars ) ) );
 }

 add_static_constraint( MinPower_Const , "MinPower_Intermittent" );

 // Active power bound constraints

 active_power_bounds_Const.resize( f_time_horizon );

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {
  active_power_bounds_Const[ t ].set_lhs( f_kappa * v_minimum_power[ t ] );
  active_power_bounds_Const[ t ].set_rhs( f_kappa * v_maximum_power[ t ] );
  active_power_bounds_Const[ t ].set_variable( &v_active_power[ t ] );
 }

 add_static_constraint( active_power_bounds_Const ,
                        "ActivePower_Bounds_Intermittent" );

 // Active power bound design constraints

 if( f_investment_cost != 0 ) {

  active_power_bounds_design_Const.resize(
   boost::multi_array< FRowConstraint , 2 >::extent_gen()
   [ f_time_horizon ][ 2 ] ); // 2 dims, i.e., the lower and upper bounds

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {

   // Lower bound of the active power design constraints:
   //
   //      v_minimum_power z <= v_active_power     z \in {0,1}, for all t
   // => 0 <= v_active_power - v_minimum_power z   z \in {0,1}, for all t

   LinearFunction::v_coeff_pair lower_vars;

   lower_vars.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );
   lower_vars.push_back( std::make_pair( &v_design ,
                                         -f_kappa * v_minimum_power[ t ] ) );

   active_power_bounds_design_Const[ t ][ 0 ].set_lhs( 0.0 );
   active_power_bounds_design_Const[ t ][ 0 ].set_rhs( Inf< double >() );
   active_power_bounds_design_Const[ t ][ 0 ].set_function(
    new LinearFunction( std::move( lower_vars ) ) );

   // Upper bound of the active power design constraints:
   //
   //      v_active_power <= v_maximum_power z     z \in {0,1}, for all t
   // => v_active_power - v_maximum_power z <= 0   z \in {0,1}, for all t

   LinearFunction::v_coeff_pair upper_vars;

   upper_vars.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );
   upper_vars.push_back( std::make_pair( &v_design ,
                                         -f_kappa * v_maximum_power[ t ] ) );

   active_power_bounds_design_Const[ t ][ 1 ].set_lhs( -Inf< double >() );
   active_power_bounds_design_Const[ t ][ 1 ].set_rhs( 0.0 );
   active_power_bounds_design_Const[ t ][ 1 ].set_function(
    new LinearFunction( std::move( upper_vars ) ) );
  }

  add_static_constraint( active_power_bounds_design_Const ,
                         "ActivePower_Bounds_Design_Intermittent" );
 }

 set_constraints_generated();

}  // end( IntermittentUnitBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

bool IntermittentUnitBlock::is_feasible( bool useabstract ,
                                         Configuration * fsbc ) {
 // Retrieve the tolerance.

 auto config = dynamic_cast< SimpleConfiguration< double > * >( fsbc );

 if( ( ! config ) && f_BlockConfig )
  config = dynamic_cast< SimpleConfiguration< double > * >
  ( f_BlockConfig->f_is_feasible_Configuration );

 // If a tolerance has not been provided, use the default tolerance.
 const auto tol = config ? config->f_value : 1.0e-8;

 return(
  UnitBlock::is_feasible( useabstract )
  // Constraints
  && Constraint::is_feasible( MinPower_Const , tol )
  && Constraint::is_feasible( MaxPower_Const , tol )
  && Constraint::is_feasible( active_power_bounds_Const , tol )
  && Constraint::is_feasible( active_power_bounds_design_Const , tol )
  // Variables
  && ColVariable::is_feasible( v_active_power , tol )
  && ColVariable::is_feasible( v_primary_spinning_reserve , tol )
  && ColVariable::is_feasible( v_secondary_spinning_reserve , tol ) );

}  // end( IntermittentUnitBlock::is_feasible )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::generate_objective( Configuration * objc ) {

 if( objective_generated() )
  return; // Objective has already been generated

 if( get_objective() != nullptr )  // an objective is there already
  return;                          // cowardly (and silently) return

 LinearFunction::v_coeff_pair vars;

 if( f_investment_cost != 0 )
  vars.push_back( std::make_pair( &v_design , f_investment_cost ) );

 objective.set_function( new LinearFunction( std::move( vars ) ) );
 objective.set_sense( Objective::eMin );

 // Set Block objective
 this->set_objective( &objective );

 set_objective_generated();

}  // end( IntermittentUnitBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*--- METHODS FOR LOADING, PRINTING & SAVING THE IntermittentUnitBlock -----*/
/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::serialize( netCDF::NcGroup & group ) const {
 UnitBlock::serialize( group );

 // Serialize scalar variables.
 ::serialize( group , "Gamma" , netCDF::NcDouble() , f_gamma );
 ::serialize( group , "Kappa" , netCDF::NcDouble() , f_kappa );

 // Serialize one-dimensional variables.
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
  else if( data.size() != 1 ) {
   throw( std::logic_error
    ( "IntermittentUnitBlock::serialize: invalid dimension for variable " +
      var_name + ": " + std::to_string( data.size() ) + ". Its dimension "
                                                        "must be one of the following: TimeHorizon, NumberIntervals, 1." ) );
  }

  ::serialize( group , var_name , ncType , dimension , data ,
               allow_scalar_var );
 };

 serialize( "MinPower" , v_minimum_power );
 serialize( "MaxPower" , v_maximum_power );
 serialize( "InertiaPower" , v_inertia_power );

}  // end( IntermittentUnitBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::update_max_power_in_constraints(
 const Subset & time ,
 ModParam issueAMod ) {
 if( ! MaxPower_Const.empty() ) {
  for( auto t : time ) {
   MaxPower_Const[ t ].set_rhs
    ( f_kappa * f_gamma * v_maximum_power[ t ] , issueAMod );
   // FIXME: use a GroupModification
  }
 }

 if( ! active_power_bounds_Const.empty() ) {
  for( auto t : time ) {
   active_power_bounds_Const[ t ].set_rhs
    ( f_kappa * v_maximum_power[ t ] , issueAMod );
   // FIXME: use a GroupModification
  }
 }
}

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::update_max_power_in_constraints(
 const Range & time , ModParam issueAMod ) {
 if( ! MaxPower_Const.empty() ) {
  for( auto t = time.first ; t < time.second ; ++t ) {
   MaxPower_Const[ t ].set_rhs
    ( f_kappa * f_gamma * v_maximum_power[ t ] , issueAMod );
   // FIXME: use a GroupModification
  }
 }

 if( ! active_power_bounds_Const.empty() ) {
  for( auto t = time.first ; t < time.second ; ++t ) {
   active_power_bounds_Const[ t ].set_rhs
    ( f_kappa * v_maximum_power[ t ] , issueAMod );
   // FIXME: use a GroupModification
  }
 }
}

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::set_maximum_power( MF_dbl_it values ,
                                               Subset && subset ,
                                               bool ordered ,
                                               ModParam issuePMod ,
                                               ModParam issueAMod ) {
 if( subset.empty() )
  return;

 if( v_maximum_power.empty() ) {
  if( std::all_of( values , values + subset.size() ,
                   []( double cst ) { return( cst == 0 ); } ) ) {
   return;
  }

  Index max_index = *std::max_element( std::begin( subset ) ,
                                       std::end( subset ) );
  v_maximum_power.assign( max_index , 0 );
 }

 // If nothing changes, return
 bool identical = true;
 for( auto t : subset ) {
  if( t >= v_maximum_power.size() ) {
   throw( std::invalid_argument( "IntermittentUnitBlock::set_maximum_power:"
                                  " invalid value in subset" ) );
  }
  auto max_power = *( values++ );
  if( v_maximum_power[ t ] != max_power ) {
   identical = false;
   if( not_dry_run( issuePMod ) )
    // Change the physical representation
    v_maximum_power[ t ] = max_power;
  }
 }
 if( identical )
  return;  // nothing changes; return

 if( not_dry_run( issuePMod ) && not_dry_run( issueAMod ) &&
     constraints_generated() ) {
  // Change the abstract representation
  update_max_power_in_constraints( subset , issueAMod );
 }

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification
  if( ! ordered ) {
   std::sort( subset.begin() , subset.end() );
  }

  Block::add_Modification( std::make_shared< IntermittentUnitBlockSbstMod >
                            ( this , IntermittentUnitBlockMod::eSetMaxP ,
                              std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );
 }
}

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::set_maximum_power( MF_dbl_it values , Range rng ,
                                               ModParam issuePMod ,
                                               ModParam issueAMod ) {
 rng.second = std::min( rng.second , f_time_horizon );
 if( rng.second <= rng.first )
  return;

 if( v_maximum_power.empty() ) {
  if( std::all_of( values , values + ( rng.second - rng.first ) ,
                   []( double cst ) { return( cst == 0 ); } ) ) {
   return;
  }

  auto max_index = rng.second;
  v_maximum_power.assign( max_index , 0 );
 }

 // If nothing changes, return
 if( std::equal( values , values + ( rng.second - rng.first ) ,
                 v_maximum_power.begin() + rng.first ) )
  return;

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  std::copy( values , values + ( rng.second - rng.first ) ,
             v_maximum_power.begin() + rng.first );

  if( not_dry_run( issueAMod ) && constraints_generated() ) {
   // Change the abstract representation
   update_max_power_in_constraints( rng , issueAMod );
  }
 }

 if( issue_pmod( issuePMod ) )
  Block::add_Modification( std::make_shared< IntermittentUnitBlockRngdMod >(
                            this , IntermittentUnitBlockMod::eSetMaxP , rng ) ,
                           Observer::par2chnl( issuePMod ) );
}

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::scale
 ( std::vector< double >::const_iterator values , Subset && subset ,
   const bool ordered , c_ModParam issuePMod , c_ModParam issueAMod ) {

 if( subset.empty() )
  return; // Since the given Subset is empty, no operation is performed

 if( f_scale == *values )
  return; // The scale factor does not change: nothing to do

 if( not_dry_run( issuePMod ) ) {
  f_scale = *values; // Update the scale factor
 }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< UnitBlockMod >
                            ( this , UnitBlockMod::eScale ) ,
                           Observer::par2chnl( issuePMod ) );
 else if( auto f_Block = get_f_Block() )
  f_Block->add_Modification( std::make_shared< UnitBlockMod >
                             ( this , UnitBlockMod::eScale ) ,
                             Observer::par2chnl( issuePMod ) );

}  // end( IntermittentUnitBlock::scale )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::set_kappa
 ( std::vector< double >::const_iterator values , Subset && subset ,
   const bool ordered , ModParam issuePMod , ModParam issueAMod ) {

 if( subset.empty() )
  return; // Since the given Subset is empty, no operation is performed

 if( f_kappa == *values )
  return; // The kappa constant does not change: nothing to do

 if( not_dry_run( issuePMod ) ) {
  f_kappa = *values; // Update the kappa constant

  if( not_dry_run( issueAMod ) ) {
   // Update the abstract representation
   if( constraints_generated() ) {
    // Update the constraints

    if( ! active_power_bounds_Const.empty() ) {
     for( Index t = 0 ; t < f_time_horizon ; ++t ) {
      active_power_bounds_Const[ t ].set_lhs
       ( f_kappa * v_minimum_power[ t ] , issueAMod );

      active_power_bounds_Const[ t ].set_rhs
       ( f_kappa * v_maximum_power[ t ] , issueAMod );
     }
    }

    if( ! MinPower_Const.empty() ) {
     for( Index t = 0 ; t < f_time_horizon ; ++t ) {
      MinPower_Const[ t ].set_lhs
       ( f_kappa * v_minimum_power[ t ] , issueAMod );
     }
    }

    if( ! MaxPower_Const.empty() ) {
     for( Index t = 0 ; t < f_time_horizon ; ++t ) {
      MaxPower_Const[ t ].set_rhs
       ( f_gamma * f_kappa * v_maximum_power[ t ] , issueAMod );
     }
    }
   }  // end( constraints_generated )
  }  // end( if( not_dry_run( issueAMod ) )
 }  // end( if( not_dry_run( issuePMod ) )

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< IntermittentUnitBlockMod >
                            ( this , IntermittentUnitBlockMod::eSetKappa ) ,
                           Observer::par2chnl( issuePMod ) );
 }
}  // end( IntermittentUnitBlock::set_kappa )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::set_kappa
 ( std::vector< double >::const_iterator values , Range rng ,
   ModParam issuePMod , ModParam issueAMod ) {

 if( rng.first >= rng.second )
  return; // An empty Range was given: no operation is performed.

 Subset subset( 1 , 0 );

 set_kappa( values , std::move( subset ) , true , issuePMod , issueAMod );

}  // end( IntermittentUnitBlock::set_kappa )

/*--------------------------------------------------------------------------*/
/*------------------- End File IntermittentUnitBlock.cpp -------------------*/
/*--------------------------------------------------------------------------*/
