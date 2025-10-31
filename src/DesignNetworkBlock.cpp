/*--------------------------------------------------------------------------*/
/*-------------------- File DesignNetworkBlock.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the DesignNetworkBlock class.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <utility>

#include "DesignNetworkBlock.h"

#include "LinearFunction.h"

#include "OneVarConstraint.h"

#include "FRealObjective.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register DesignNetworkBlock to the Block factory
SMSpp_insert_in_factory_cpp_0( DesignNetworkBlock );

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS OF DesignNetworkBlock --------------------*/
/*--------------------------------------------------------------------------*/

DesignNetworkBlock::~DesignNetworkBlock()
{
 Constraint::clear( v_design_bound_const );

 objective.clear();
}

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void DesignNetworkBlock::deserialize_network_blocks( const netCDF::NcGroup & group )
{
 v_network_blocks.clear();
 Index i = 0; Index found = 0;

 for( ; ; ++i ) {
  std::string sub_group_name = "NetworkBlock_" + std::to_string( i );
  auto sub_group = group.getGroup( sub_group_name );
  if( sub_group.isNull() ) break;

  auto res = almost_new_Block( sub_group , this );
  if( auto nbi = dynamic_cast< NetworkBlock * >( res.second ) ) {
   nbi->set_NetworkData( f_NetworkData );
   nbi->deserialize( res.first );
   v_network_blocks.push_back( nbi );
   ++found;
  }
  else {
   delete res.second;
   throw std::invalid_argument(
     "DesignNetworkBlock::deserialize: " + sub_group_name + " not a NetworkBlock" );
  }
 }

 if( ! found )
  v_network_blocks.clear();
}

/*--------------------------------------------------------------------------*/

void DesignNetworkBlock::deserialize( const netCDF::NcGroup & group )
{
#ifndef NDEBUG
 static std::vector< std::string > expected_dims = { "NumberDesignLines" };
 check_dimensions( group , expected_dims , std::cerr );

 static std::vector< std::string > expected_vars = {
  "InvestmentCost" , "MinCapacityDesign" , "MaxCapacityDesign" };
 check_variables( group , expected_vars , std::cerr );
#endif

 deserialize_dim( group , "NumberDesignLines" , f_number_design_lines , false );

 if( ! ::deserialize( group , "InvestmentCost" , f_number_design_lines ,
                      v_InvestmentCost , true , true ) ) {
  v_InvestmentCost.resize( f_number_design_lines );
 }

 if( f_number_design_lines > 0 ) {
  std::vector< Index > tmp_design_lines;
  if( ::deserialize( group , "DesignLines" , f_number_design_lines ,
                     tmp_design_lines , false , true ) ) {
   v_design_lines = std::move( tmp_design_lines );
   }
  else {
   v_design_lines.resize( f_number_design_lines );
   for( Index p = 0 ; p < f_number_design_lines ; ++p )
    v_design_lines[ p ] = p;
   }
  }

 if( ! ::deserialize( group , "MinCapacityDesign" , f_number_design_lines ,
                      v_MinCapacityDesign , true , true ) )
  v_MinCapacityDesign.resize( f_number_design_lines );

 if( ! ::deserialize( group , "MaxCapacityDesign" , f_number_design_lines ,
                      v_MaxCapacityDesign , true , true ) )
  v_MaxCapacityDesign.resize( f_number_design_lines , 1 );

 deserialize_network_blocks( group );

 // finally call the method of the base class
 NetworkBlock::deserialize( group );

 check_data_consistency();

 }  // end( DesignNetworkBlock::deserialize )

/*--------------------------------------------------------------------------*/

void DesignNetworkBlock::check_data_consistency( void ) const
{
 for( Index l = 0 ; l < f_number_design_lines ; ++l ) {
  // Min/Max capacity design
  if( get_min_capacity_design( l ) < 0 )
   throw( std::logic_error( "DesignNetworkBlock::check_data_consistency: "
                            "MinCapacityDesign must be nonnegative." ) );

  // Continuous case (MaxCapacityDesign > 0):
  // MinCapacityDesign <= MaxCapacityDesign
  if( ( get_max_capacity_design( l ) > 0 ) &&
      ( get_min_capacity_design( l ) > get_max_capacity_design( l ) ) )
   throw( std::logic_error( "DesignNetworkBlock::check_data_consistency: "
                            "MinCapacityDesign > MaxCapacityDesign." ) );

  // Unitary case (|MaxCapacityDesign| == 1): MinCapacityDesign <= 1
  if( ( std::abs( get_max_capacity_design( l ) ) == 1 ) &&
      ( get_min_capacity_design( l ) > 1.0 ) )
   throw( std::logic_error( "DesignNetworkBlock::check_data_consistency: "
                            "MinCapacityDesign must be <= 1 when "
			    "|MaxCapacityDesign| == 1" ) );

  // Binary case (max < 0): MinCapacityDesign <= 1
  if( ( get_max_capacity_design( l ) < 0 ) &&
      ( get_min_capacity_design( l ) > 1.0 ) )
   throw( std::logic_error( "DesignNetworkBlock::check_data_consistency: "
                            "MinCapacityDesign must be <= 1 for binary "
			    "design." ) );
  }
 }  // end( DesignNetworkBlock::check_data_consistency )

/*--------------------------------------------------------------------------*/
/*----------------------- GENERATION METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

void DesignNetworkBlock::generate_abstract_variables( Configuration * stvv )
{
 if( variables_generated() )  // variables have already been generated
  return;                     // nothing to do

 NetworkBlock::generate_abstract_variables( stvv );

 // Create design variables only for the selected ("designed") lines
 if( const auto nd = static_cast< Index >( v_design_lines.size() ) ) {
  v_design.resize( nd );
  for( Index p = 0 ; p < nd ; ++p ) {
   if( get_max_capacity_design( p ) < 0 )
    v_design[ p ].set_type( ColVariable::kBinary );
   else
    v_design[ p ].set_type( ColVariable::kNonNegative );
   }
  add_static_variable( v_design , "x_network" );
 }

 set_variables_generated();
}  // end( DesignNetworkBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void DesignNetworkBlock::generate_abstract_constraints( Configuration * stcc )
{
 if( constraints_generated() )  // constraints have already been generated
  return;                       // nothing to do

 if( const auto nd = static_cast< Index >( v_design_lines.size() ) ) {
  v_design_bound_const.resize( nd );

  for( Index p = 0 ; p < nd ; ++p ) {

    const double lb = std::max( 0.0 , get_min_capacity_design( p ) );
    const double maxd = get_max_capacity_design( p );
    const bool is_binary = ( maxd < 0.0 );
    const double ub = is_binary
                       ? 1.0 : ( std::abs( maxd ) == 1.0
                            ? 1.0 : std::abs( maxd ) );

    if( ( lb == 1.0 ) && ( ub == 1.0 ) ) {
     v_design[ p ].is_unitary( true , eNoMod );
    }
    else {
     v_design_bound_const[ p ].set_lhs( lb );
     v_design_bound_const[ p ].set_rhs( ub );
     v_design_bound_const[ p ].set_variable( &v_design[ p ] );
    }

    if( is_binary )
     v_design[ p ].is_integer( true , eNoMod );
  }

  add_static_constraint( v_design_bound_const , "DesignBound_Network" );
 }

 set_constraints_generated();
}  // end( DesignNetworkBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void DesignNetworkBlock::generate_objective( Configuration * objc )
{
 if( objective_generated() )  // Objective has already been generated
  return;                     // nothing to do

 auto lf = new LinearFunction();

 // Investment term only over the selected ("designed") lines
 const auto nd = static_cast< Index >( v_design_lines.size() );
 for( Index p = 0 ; p < nd ; ++p ) {
  if( get_investment_cost( p ) != 0 )
   lf->add_variable( &v_design[ p ] , get_investment_cost( p ) , eNoMod );
 }

 objective.set_function( lf );
 objective.set_sense( Objective::eMin );

 this->set_objective( & objective );  // set Block objective

 set_objective_generated();
}  // end( DesignNetworkBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*--------------- METHODS FOR CHECKING THE DesignNetworkBlock --------------*/
/*--------------------------------------------------------------------------*/

bool DesignNetworkBlock::is_feasible( bool useabstract , Configuration * fsbc )
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
  NetworkBlock::is_feasible( useabstract )
  // Variables
  && ColVariable::is_feasible( v_design , tol )
  // Constraints
  && RowConstraint::is_feasible( v_design_bound_const , tol , rel_viol ) );
} // end( DesignNetworkBlock::is_feasible )

/*--------------------------------------------------------------------------*/
/*------- METHODS FOR LOADING, PRINTING & SAVING THE DesignNetworkBlock ----*/
/*--------------------------------------------------------------------------*/

void DesignNetworkBlock::serialize( netCDF::NcGroup & group ) const
{
 NetworkBlock::serialize( group );

 auto NumberDesignLines = group.addDim( "NumberDesignLines" ,
                                        f_number_design_lines );

 if( ! v_InvestmentCost.empty() ) {
  ::serialize( group , "InvestmentCost" , netCDF::NcDouble() , NumberDesignLines ,
               v_InvestmentCost );
 }

 const Index nd = static_cast< Index >( v_design_lines.size() );
 if( nd > 0 ) {

  bool is_sequential = true;
  for( Index p = 0 ; p < nd ; ++p )
   if( v_design_lines[ p ] != p ) { is_sequential = false; break; }

  if( ! is_sequential ) {
   ::serialize( group , "DesignLines" , netCDF::NcInt() , NumberDesignLines ,
                v_design_lines );
  }

  std::vector< double > min_on_nd( nd , 0.0 ) , max_on_nd( nd , 1.0 );
  for( Index p = 0 ; p < nd ; ++p ) {
   min_on_nd[ p ] = v_MinCapacityDesign.empty() ? 0.0 : v_MinCapacityDesign[ p ];
   max_on_nd[ p ] = v_MaxCapacityDesign.empty() ? 1.0 : v_MaxCapacityDesign[ p ];
  }

  if( std::any_of( min_on_nd.begin() , min_on_nd.end() ,
                   []( double x ){ return( x != 0.0 ); } ) )
   ::serialize( group , "MinCapacityDesign" , netCDF::NcDouble() ,
                NumberDesignLines , min_on_nd );

  if( std::any_of( max_on_nd.begin() , max_on_nd.end() ,
                   []( double x ){ return( std::abs( x ) != 1.0 ); } ) )
   ::serialize( group , "MaxCapacityDesign" , netCDF::NcDouble() ,
                NumberDesignLines , max_on_nd );
 }
}  // end( DesignNetworkBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

void DesignNetworkBlock::set_active_demand( MF_dbl_it values ,
                                            Subset && subset ,
                                            bool ordered ,
                                            c_ModParam issuePMod ,
                                            c_ModParam issueAMod )
{
  if( v_network_blocks.empty() || subset.empty() )
   return;

  std::vector< double > vals( subset.size() );
  std::copy( values , values + subset.size() , vals.begin() );

  const Subset base_subset = subset;
  for( auto * nb : v_network_blocks ) {
   if( ! nb ) continue;
   Subset sb = base_subset;
   nb->set_active_demand( vals.begin() ,
                          std::move( sb ) ,
                          ordered ,
                          issuePMod ,
                          issueAMod );
  }
}  // end( DesignNetworkBlock::set_active_demand( subset ) )

/*--------------------------------------------------------------------------*/

void DesignNetworkBlock::set_active_demand( MF_dbl_it values ,
                                            Range rng ,
                                            c_ModParam issuePMod ,
                                            c_ModParam issueAMod )
{
  if( v_network_blocks.empty() )
    return;

  if( rng.second <= rng.first )
    return;

  const Index len = rng.second - rng.first;

  std::vector<double> vals( len );
  std::copy( values , values + len , vals.begin() );

  for( auto * nb : v_network_blocks ) {
    if( ! nb ) continue;
    nb->set_active_demand( vals.begin() , rng , issuePMod , issueAMod );
  }
}  // end( DesignNetworkBlock::set_active_demand( range ) )

/*--------------------------------------------------------------------------*/
/*-------------------- End File DesignNetworkBlock.cpp ---------------------*/
/*--------------------------------------------------------------------------*/
