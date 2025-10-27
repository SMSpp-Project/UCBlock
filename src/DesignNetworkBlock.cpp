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

void DesignNetworkBlock::deserialize( const netCDF::NcGroup & group )
{
#ifndef NDEBUG
 static std::vector< std::string > expected_dims = { "NumberLines" };
 check_dimensions( group , expected_dims , std::cerr );

 static std::vector< std::string > expected_vars = {
  "InvestmentCost" , "MinCapacityDesign" , "MaxCapacityDesign" };
 check_variables( group , expected_vars , std::cerr );
#endif

 NetworkBlock::deserialize( group );

 /*--- base cardinality of lines ------------------------------------------*/
 deserialize_dim( group , "NumberLines" , f_number_lines , false );

 /*--- investment cost: per line over NumberLines (or scalar) -------------*/
 if( ! ::deserialize( group , "InvestmentCost" , f_number_lines ,
                      v_InvestmentCost , true , true ) ) {
  v_InvestmentCost.resize( f_number_lines , 0 );
 }

 /*--- subset of designed lines: NumDesignLines [+ optional DesignLines] ---*/
 Index nd = 0;
 deserialize_dim( group , "NumDesignLines" , nd , false );  // optional

 v_design_lines.clear();
 line2pos.clear();

 if( nd == 0 ) {
  // fallback: if NumDesignLines is absent/zero, assume all lines are designed
  nd = f_number_lines;
 }

 if( nd > 0 ) {
  // try read DesignLines[ nd ] (optional)
  std::vector< Index > tmp_design_lines;
  if( ::deserialize( group , "DesignLines" , nd , tmp_design_lines , false , true ) ) {
   v_design_lines = std::move( tmp_design_lines );
  }
  else {
   // implicit convention: {0,1,...,nd-1}
   v_design_lines.resize( nd );
   for( Index p = 0 ; p < nd ; ++p ) v_design_lines[ p ] = p;
  }

  // build map line -> position
  for( Index p = 0 ; p < nd ; ++p )
   line2pos[ v_design_lines[ p ] ] = p;
 }

 /*--- Min/Max capacity design --------------------------------------------*/
 // We store per-line arrays sized to NumberLines (getters expect indexing by line).
 // Defaults for non-designed lines: Min=0, Max=1.
 v_MinCapacityDesign.assign( f_number_lines , 0 );
 v_MaxCapacityDesign.assign( f_number_lines , 1 );

 bool have_Min_over_NumberLines =
  ::deserialize( group , "MinCapacityDesign" , f_number_lines ,
                 v_MinCapacityDesign , true , true );

 bool have_Max_over_NumberLines =
  ::deserialize( group , "MaxCapacityDesign" , f_number_lines ,
                 v_MaxCapacityDesign , true , true );

 if( ! have_Min_over_NumberLines ) {
  // try "NumDesignLines"-indexed (or scalar replicated over nd)
  std::vector< double > tmpMin;
  if( ( nd > 0 ) &&
      ::deserialize( group , "MinCapacityDesign" , nd , tmpMin , true , true ) ) {
   for( Index p = 0 ; p < nd ; ++p ) {
    const Index l = v_design_lines[ p ];
    v_MinCapacityDesign[ l ] = tmpMin[ p ];
   }
  }
  // else: keep defaults (0) already set
 }

 if( ! have_Max_over_NumberLines ) {
  // try "NumDesignLines"-indexed (or scalar replicated over nd)
  std::vector< double > tmpMax;
  if( ( nd > 0 ) &&
      ::deserialize( group , "MaxCapacityDesign" , nd , tmpMax , true , true ) ) {
   for( Index p = 0 ; p < nd ; ++p ) {
    const Index l = v_design_lines[ p ];
    v_MaxCapacityDesign[ l ] = tmpMax[ p ];
   }
  }
  // else: keep defaults (1) already set
 }

 check_data_consistency();
}  // end( DesignNetworkBlock::deserialize )

/*--------------------------------------------------------------------------*/

void DesignNetworkBlock::check_data_consistency( void ) const
{
 for( Index l = 0 ; l < f_number_lines ; ++l ) {

  // Min/Max capacity design
  if( get_min_capacity_design( l ) < 0 )
   throw( std::logic_error( "DesignNetworkBlock::check_data_consistency: "
                            "MinCapacityDesign must be nonnegative." ) );

  // Continuous case (MaxCapacityDesign > 0): MinCapacityDesign <= MaxCapacityDesign
  if( ( get_max_capacity_design( l ) > 0 ) &&
      ( get_min_capacity_design( l ) > get_max_capacity_design( l ) ) )
   throw( std::logic_error( "DesignNetworkBlock::check_data_consistency: "
                            "MinCapacityDesign > MaxCapacityDesign." ) );

  // Unitary case (|MaxCapacityDesign| == 1): MinCapacityDesign <= 1
  if( ( std::abs( get_max_capacity_design( l ) ) == 1 ) &&
      ( get_min_capacity_design( l ) > 1.0 ) )
   throw( std::logic_error( "DesignNetworkBlock::check_data_consistency: "
                            "MinCapacityDesign must be <= 1 when |MaxCapacityDesign| == 1." ) );

  // Binary case (max < 0): MinCapacityDesign <= 1
  if( ( get_max_capacity_design( l ) < 0 ) &&
      ( get_min_capacity_design( l ) > 1.0 ) )
   throw( std::logic_error( "DesignNetworkBlock::check_data_consistency: "
                            "MinCapacityDesign must be <= 1 for binary design." ) );
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
   const Index l = v_design_lines[ p ];
   if( get_max_capacity_design( l ) < 0 )
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

 // Bounds only for the selected ("designed") lines
 if( const auto nd = static_cast< Index >( v_design_lines.size() ) ) {
  v_design_bound_const.resize( nd );
  bool any_bound = false;

  for( Index p = 0 ; p < nd ; ++p ) {
    const Index l = v_design_lines[ p ];

    const double maxd = get_max_capacity_design( l );
    const double lb = std::max( 0.0 , get_min_capacity_design( l ) );
    const double ub = ( std::abs( maxd ) == 1 ? 1.0 : std::abs( maxd ) );

    if( ( lb > 0.0 ) || ( std::abs( maxd ) != 1 ) ) {
     v_design_bound_const[ p ].set_lhs( lb );
     v_design_bound_const[ p ].set_rhs( ub );
     v_design_bound_const[ p ].set_variable( &v_design[ p ] );
     any_bound = true;
    }
    else {
     v_design[ p ].is_unitary( true , eNoMod );
    }

    if( maxd < 0 )
     v_design[ p ].is_integer( true , eNoMod );
  }

  if( any_bound )
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
  const Index l = v_design_lines[ p ];
  if( get_investment_cost( l ) != 0 )
   lf->add_variable( &v_design[ p ] , get_investment_cost( l ) , eNoMod );
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
  && ColVariable::is_feasible( v_node_injection , tol )
  && ColVariable::is_feasible( v_design , tol )
  // Constraints
  && RowConstraint::is_feasible( v_design_bound_const , tol , rel_viol )
  && RowConstraint::is_feasible( node_injection_bounds_const , tol , rel_viol ) );
} // end( DCNetworkBlock::is_feasible )

/*--------------------------------------------------------------------------*/
/*------- METHODS FOR LOADING, PRINTING & SAVING THE DesignNetworkBlock ----*/
/*--------------------------------------------------------------------------*/

void DesignNetworkBlock::serialize( netCDF::NcGroup & group ) const
{
 NetworkBlock::serialize( group );

 // If there is no design to write, return
 if( f_number_lines == 0 )
  return;

 // Always write NumberLines
 auto NumberLines = group.addDim( "NumberLines" , f_number_lines );

 if( ! v_InvestmentCost.empty() ) {
  ::serialize( group , "InvestmentCost" , netCDF::NcDouble() , NumberLines ,
               v_InvestmentCost );
 }

 // subset info
 const Index nd = static_cast< Index >( v_design_lines.size() );
 if( nd > 0 ) {
  auto NumDesignLines = group.addDim( "NumDesignLines" , nd );

  // write DesignLines only if not the implicit sequence 0 to nd-1
  bool is_sequential = true;
  for( Index p = 0 ; p < nd ; ++p )
   if( v_design_lines[ p ] != p ) { is_sequential = false; break; }

  if( ! is_sequential ) {
   ::serialize( group , "DesignLines" , netCDF::NcInt() , NumDesignLines ,
                v_design_lines );
  }

  // write Min/Max indexed over NumDesignLines if any relevant info
  std::vector< double > min_on_nd( nd , 0.0 ) , max_on_nd( nd , 1.0 );
  for( Index p = 0 ; p < nd ; ++p ) {
   const Index l = v_design_lines[ p ];
   min_on_nd[ p ] = v_MinCapacityDesign.empty() ? 0.0 : v_MinCapacityDesign[ l ];
   max_on_nd[ p ] = v_MaxCapacityDesign.empty() ? 1.0 : v_MaxCapacityDesign[ l ];
  }

  if( std::any_of( min_on_nd.begin() , min_on_nd.end() ,
                   []( double x ){ return( x != 0.0 ); } ) )
   ::serialize( group , "MinCapacityDesign" , netCDF::NcDouble() , NumDesignLines ,
                min_on_nd );

  if( std::any_of( max_on_nd.begin() , max_on_nd.end() ,
                   []( double x ){ return( std::abs( x ) != 1.0 ); } ) )
   ::serialize( group , "MaxCapacityDesign" , netCDF::NcDouble() , NumDesignLines ,
                max_on_nd );
 }
}  // end( DesignNetworkBlock::serialize )

/*--------------------------------------------------------------------------*/
/*-------------------- End File DesignNetworkBlock.cpp ---------------------*/
/*--------------------------------------------------------------------------*/
