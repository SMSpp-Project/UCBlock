/*--------------------------------------------------------------------------*/
/*--------------------- File OTSNetworkBlock.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the OTSNetworkBlock class.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Claude 4.6 \n
 *
 * \copyright &copy; by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "OTSNetworkBlock.h"

#include "LinearFunction.h"

#include "OneVarConstraint.h"

#include "FRealObjective.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*------------------------- FACTORY REGISTRATION ---------------------------*/
/*--------------------------------------------------------------------------*/

// register OTSNetworkBlock to the Block factory
SMSpp_insert_in_factory_cpp_1( OTSNetworkBlock );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// register DCNetworkData to the NetworkData factory

using OTSNetworkData = OTSNetworkBlock::OTSNetworkData;

SMSpp_insert_in_factory_cpp_0( OTSNetworkData );

/*--------------------------------------------------------------------------*/
/*-------------------- METHODS OF OTSNetworkData ---------------------------*/
/*--------------------------------------------------------------------------*/

void OTSNetworkData::deserialize( const netCDF::NcGroup & group )
{
 // deserialize the parent DCNetworkData first
 DCNetworkData::deserialize( group );

 // read the optional "SwitchingCost" variable:
 // - if absent, v_switching_cost stays empty (all zeros by default)
 // - if scalar (dimension 0), the value is replicated for every line
 // - if vector of size NumberLines, each entry is the per-line cost
 // The ::deserialize() from SMSTypedefs handles scalar replication
 // automatically when allow_scalar_var = true.
 if( ::deserialize( group , "SwitchingCost" , f_number_lines ,
                    v_switching_cost , true , true ) )
  if( std::all_of( v_switching_cost.begin() , v_switching_cost.end() ,
		   []( auto a ) { return( a == 0 ); } ) )
   // if all costs are zero, clear the vector for efficiency: the rest
   // of the code treats an empty vector as "no switching costs"
   v_switching_cost.clear();

 }  // end( OTSNetworkData::deserialize )

/*--------------------------------------------------------------------------*/
/*------------------- METHODS OF OTSNetworkBlock ---------------------------*/
/*--------------------------------------------------------------------------*/

OTSNetworkBlock::~OTSNetworkBlock()
{
 // clear OTS-specific constraints
 Constraint::clear( v_OTS_KVL_upper );
 Constraint::clear( v_OTS_KVL_lower );
 Constraint::clear( v_OTS_flow_upper );
 Constraint::clear( v_OTS_flow_lower );
 Constraint::clear( v_switching_exclusivity );
 Constraint::clear( v_elastic_precedence );
 Constraint::clear( v_design_coupling );
 }

/*--------------------------------------------------------------------------*/

void OTSNetworkBlock::deserialize( const netCDF::NcGroup & group )
{
 // DCNetworkBlock::deserialize will call get_new_NetworkData() which
 // returns an OTSNetworkData, so "SwitchingCost" will be read
 DCNetworkBlock::deserialize( group );
 }

/*--------------------------------------------------------------------------*/

void OTSNetworkBlock::generate_abstract_variables( Configuration * stvv )
{
 if( variables_generated() )
  return;

 // Force Kirchhoff base formulation: OTS always needs angle variables.
 // We call NetworkBlock::generate_abstract_variables (grandparent), then
 // directly generate Kirchhoff variables, bypassing the parent's
 // formulation switch.
 NetworkBlock::generate_abstract_variables( stvv );

 ftype = KIRCHHOFF;
 generate_PTDF_variables();
 generate_KIRCHHOFF_variables();

 // Read OTS configuration from SimpleConfiguration< int >
 Configuration * cfg = stvv;
 if( ( ! cfg ) && f_BlockConfig )
  cfg = f_BlockConfig->f_static_variables_Configuration;

 Index ots_cfg = 0;
 if( auto sci = dynamic_cast< SimpleConfiguration< int > * >( cfg ) )
  ots_cfg = sci->f_value;

 // bits 0-1: formulation type
 f_ots_type = static_cast< ots_formulation_type >( ots_cfg & 0x3 );

 // compute BigM values for all DC lines
 compute_big_M();

 // generate formulation-specific OTS variables
 switch( f_ots_type ) {
  case kOTS_Standard:
   generate_OTS_Standard_variables();
   break;
  case kOTS_Directional:
   generate_OTS_Directional_variables();
   break;
  case kOTS_Elastic:
   generate_OTS_Elastic_variables();
   break;
  case kOTS_ElasticDirectional:
   generate_OTS_ElasticDirectional_variables();
   break;
  default:
   throw( std::logic_error(
    "OTSNetworkBlock::generate_abstract_variables: unknown OTS type" ) );
  }

 set_variables_generated();

 }  // end( OTSNetworkBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void OTSNetworkBlock::generate_abstract_constraints( Configuration * stcc )
{
 if( constraints_generated() )
  return;

 // 1. NetworkCost auxiliary constraints (|F_l| linearization)
 generate_network_cost_constraints();

 // 2-5. BigM KVL, flow bounds, exclusivity, precedence (per formulation)
 switch( f_ots_type ) {
  case kOTS_Standard:
   generate_OTS_Standard_constraints();
   break;
  case kOTS_Directional:
   generate_OTS_Directional_constraints();
   break;
  case kOTS_Elastic:
   generate_OTS_Elastic_constraints();
   break;
  case kOTS_ElasticDirectional:
   generate_OTS_ElasticDirectional_constraints();
   break;
  }

 // 6. OTS flow bounds (coupled with switching variables)
 generate_OTS_flow_bounds();

 // 7. Design coupling constraints (z <= x or z+ + z- <= x)
 generate_design_coupling_constraints();

 // 8. Reference angle theta_ref = 0
 generate_reference_angle_constraint();

 // 9. KCL node balance at every node
 generate_node_balance_constraints();

 // 10. Standard box bounds for HVDC lines (no switching on HVDC)
 generate_bound_constraints();

 set_constraints_generated();

 }  // end( OTSNetworkBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void OTSNetworkBlock::generate_objective( Configuration * objc )
{
 // call the parent to build the base objective (network cost terms,
 // constant term, sense, registration)
 DCNetworkBlock::generate_objective( objc );

 // the parent only adds NetworkCost terms for the PTDF formulation;
 // since OTS always uses KIRCHHOFF, add them here if present
 auto * lf = dynamic_cast< LinearFunction * >( objective.get_function() );
 if( ! lf )
  return;

 if( ! f_NetworkData->get_network_cost().empty() )
  for( Index l = 0 ; l < get_number_lines() ; ++l )
   lf->add_variable( & v_auxiliary_variable[ l ] ,
                     f_NetworkData->get_network_cost()[ l ] , eNoMod );

 // if no switching costs, nothing more to add
 if( ( ! f_NetworkData ) ||
     ( ! static_cast< OTSNetworkData * >( f_NetworkData
					  )->has_switching_cost() ) )
  return;

 // switching cost terms: c_sw_l * (1 - z_l) = c_sw_l - c_sw_l * z_l
 // The constant part (sum of c_sw_l) is added to the existing constant
 // term; the variable part (-c_sw_l * z_l) is added to the LinearFunction.
 auto & DC_lines = f_NetworkData->get_DC_lines();
 double constant_add = 0.0;

 for( Index idx = 0 ; idx < DC_lines.size() ; ++idx ) {
  auto l = DC_lines[ idx ];
  double c_sw = static_cast< OTSNetworkData * >( f_NetworkData
						 )->get_switching_cost( l );
  if( c_sw == 0.0 )
   continue;

  constant_add += c_sw;

  switch( f_ots_type ) {
   case kOTS_Standard:
   case kOTS_Elastic:
    lf->add_variable( & v_switching[ idx ] , -c_sw , eNoMod );
    break;
   case kOTS_Directional:
   case kOTS_ElasticDirectional:
    lf->add_variable( & v_switching_pos[ idx ] , -c_sw , eNoMod );
    lf->add_variable( & v_switching_neg[ idx ] , -c_sw , eNoMod );
    break;
   }
  }

 if( constant_add > 0.0 )
  lf->set_constant_term( lf->get_constant_term() + constant_add );

 }  // end( OTSNetworkBlock::generate_objective )

/*--------------------------------------------------------------------------*/

bool OTSNetworkBlock::is_feasible( bool useabstract , Configuration * fsbc )
{
 // first check parent feasibility (includes all DCNetworkBlock checks)
 if( ! DCNetworkBlock::is_feasible( useabstract , fsbc ) )
  return( false );

 // retrieve tolerance and violation type (same logic as parent)
 double tol = 0;
 bool rel_viol = true;

 auto extract_parameters = [ & tol , & rel_viol ]( Configuration * c )
  -> bool {
  if( auto tc = dynamic_cast< SimpleConfiguration< double > * >( c ) ) {
   tol = tc->f_value;
   return( true );
   }
  if( auto tc = dynamic_cast< SimpleConfiguration<
                                     std::pair< double , int > > * >( c ) ) {
   tol = tc->f_value.first;
   rel_viol = tc->f_value.second;
   return( true );
   }
  return( false );
  };

 if( ( ! extract_parameters( fsbc ) ) && f_BlockConfig )
  extract_parameters( f_BlockConfig->f_is_feasible_Configuration );

 return(
  // OTS variables
  ColVariable::is_feasible( v_switching , tol )
  && ColVariable::is_feasible( v_switching_pos , tol )
  && ColVariable::is_feasible( v_switching_neg , tol )
  && ColVariable::is_feasible( v_elastic , tol )
  // OTS constraints
  && RowConstraint::is_feasible( v_OTS_KVL_upper , tol , rel_viol )
  && RowConstraint::is_feasible( v_OTS_KVL_lower , tol , rel_viol )
  && RowConstraint::is_feasible( v_OTS_flow_upper , tol , rel_viol )
  && RowConstraint::is_feasible( v_OTS_flow_lower , tol , rel_viol )
  && RowConstraint::is_feasible( v_switching_exclusivity , tol , rel_viol )
  && RowConstraint::is_feasible( v_elastic_precedence , tol , rel_viol )
  && RowConstraint::is_feasible( v_design_coupling , tol , rel_viol )
  );

 }  // end( OTSNetworkBlock::is_feasible )

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED METHODS ----------------------------------*/
/*--------------------------------------------------------------------------*/

void OTSNetworkBlock::compute_big_M( void )
{
 // Fattahi-Lavaei-Atamturk (2019):
 //   M_l = |B_l| * sum_k( f_max_k / |B_k| )
 // where the sum runs over all DC lines (nonzero susceptance).

 auto & DC_lines = f_NetworkData->get_DC_lines();
 if( DC_lines.empty() ) {
  v_big_M.clear();
  return;
  }

 // compute the global sum: sum_k( f_max_k / |B_k| )
 double sum_ratio = 0.0;
 for( auto & k : DC_lines ) {
  double B_k = std::abs( get_line_susceptance( k ) );
  if( B_k > 0.0 )
   sum_ratio += f_NetworkData->get_max_power_flow( k ) / B_k;
  }

 // compute M_l for each DC line
 v_big_M.resize( DC_lines.size() );
 for( Index idx = 0 ; idx < DC_lines.size() ; ++idx ) {
  auto l = DC_lines[ idx ];
  double B_l = std::abs( get_line_susceptance( l ) );
  v_big_M[ idx ] = B_l * sum_ratio;
  }

 // for elastic formulations, also compute alpha and beta
 if( ( f_ots_type == kOTS_Elastic ) ||
     ( f_ots_type == kOTS_ElasticDirectional ) ) {
  v_alpha.resize( DC_lines.size() );
  v_beta.resize( DC_lines.size() );
  for( Index idx = 0 ; idx < DC_lines.size() ; ++idx ) {
   auto l = DC_lines[ idx ];
   double f_max = f_NetworkData->get_max_power_flow( l );
   if( v_big_M[ idx ] > 0.0 ) {
    v_alpha[ idx ] = f_max / v_big_M[ idx ];
    v_beta[ idx ] = 1.0 - v_alpha[ idx ];
    }
   else {
    v_alpha[ idx ] = 1.0;
    v_beta[ idx ] = 0.0;
    }
   }
  }

 }  // end( OTSNetworkBlock::compute_big_M )

/*--------------------------------------------------------------------------*/
/*------------------- VARIABLE GENERATION METHODS --------------------------*/
/*--------------------------------------------------------------------------*/

void OTSNetworkBlock::generate_OTS_Standard_variables( void )
{
 auto & DC_lines = f_NetworkData->get_DC_lines();
 auto n_dc = DC_lines.size();
 if( n_dc == 0 )
  return;

 v_switching.resize( n_dc );
 for( auto & var : v_switching )
  var.set_type( ColVariable::kBinary );
 add_static_variable( v_switching , "OTS_switching" );
 }

/*--------------------------------------------------------------------------*/

void OTSNetworkBlock::generate_OTS_Directional_variables( void )
{
 auto & DC_lines = f_NetworkData->get_DC_lines();
 auto n_dc = DC_lines.size();
 if( n_dc == 0 )
  return;

 v_switching_pos.resize( n_dc );
 for( auto & var : v_switching_pos )
  var.set_type( ColVariable::kBinary );
 add_static_variable( v_switching_pos , "OTS_switching_pos" );

 v_switching_neg.resize( n_dc );
 for( auto & var : v_switching_neg )
  var.set_type( ColVariable::kBinary );
 add_static_variable( v_switching_neg , "OTS_switching_neg" );
 }

/*--------------------------------------------------------------------------*/

void OTSNetworkBlock::generate_OTS_Elastic_variables( void )
{
 // binary z_l
 generate_OTS_Standard_variables();

 // continuous z1_l in [0, 1]
 auto & DC_lines = f_NetworkData->get_DC_lines();
 auto n_dc = DC_lines.size();
 if( n_dc == 0 )
  return;

 v_elastic.resize( n_dc );
 for( auto & var : v_elastic ) {
  var.set_type( ColVariable::kContinuous );
  var.set_value( 0.0 );
  }
 add_static_variable( v_elastic , "OTS_elastic" );
 }

/*--------------------------------------------------------------------------*/

void OTSNetworkBlock::generate_OTS_ElasticDirectional_variables( void )
{
 // binary z+_l, z-_l
 generate_OTS_Directional_variables();

 // continuous z1_l in [0, 1]
 auto & DC_lines = f_NetworkData->get_DC_lines();
 auto n_dc = DC_lines.size();
 if( n_dc == 0 )
  return;

 v_elastic.resize( n_dc );
 for( auto & var : v_elastic ) {
  var.set_type( ColVariable::kContinuous );
  var.set_value( 0.0 );
  }
 add_static_variable( v_elastic , "OTS_elastic" );
 }

/*--------------------------------------------------------------------------*/
/*----------------- CONSTRAINT GENERATION METHODS --------------------------*/
/*--------------------------------------------------------------------------*/

void OTSNetworkBlock::generate_OTS_Standard_constraints( void )
{
 // Standard BigM KVL relaxation for each DC line:
 //   F_l - B_l * (theta_from - theta_to) <= M_l * (1 - z_l)
 //   F_l - B_l * (theta_from - theta_to) >= -M_l * (1 - z_l)
 //
 // Rewritten as:
 //   F_l - B_l*theta_from + B_l*theta_to + M_l*z_l <= M_l
 //   F_l - B_l*theta_from + B_l*theta_to - M_l*z_l >= -M_l

 auto & DC_lines = f_NetworkData->get_DC_lines();
 auto n_dc = DC_lines.size();
 if( n_dc == 0 )
  return;

 const auto & start_line = f_NetworkData->get_start_line();
 const auto & end_line = f_NetworkData->get_end_line();
 LinearFunction::v_coeff_pair vars;

 v_OTS_KVL_upper.resize( n_dc );
 v_OTS_KVL_lower.resize( n_dc );

 for( Index idx = 0 ; idx < n_dc ; ++idx ) {
  auto l = DC_lines[ idx ];
  double B_l = get_line_susceptance( l );
  double M_l = v_big_M[ idx ];

  // upper: F_l - B_l*theta_from + B_l*theta_to + M_l*z_l <= M_l
  vars.emplace_back( & v_power_flow[ l ] , 1.0 );
  vars.emplace_back( & v_voltage_angle[ start_line[ l ] ] , -B_l );
  vars.emplace_back( & v_voltage_angle[ end_line[ l ] ] , B_l );
  vars.emplace_back( & v_switching[ idx ] , M_l );
  v_OTS_KVL_upper[ idx ].set_lhs( -Inf< double >() );
  v_OTS_KVL_upper[ idx ].set_rhs( M_l );
  v_OTS_KVL_upper[ idx ].set_function(
                                   new LinearFunction( std::move( vars ) ) );

  // lower: F_l - B_l*theta_from + B_l*theta_to - M_l*z_l >= -M_l
  vars.emplace_back( & v_power_flow[ l ] , 1.0 );
  vars.emplace_back( & v_voltage_angle[ start_line[ l ] ] , -B_l );
  vars.emplace_back( & v_voltage_angle[ end_line[ l ] ] , B_l );
  vars.emplace_back( & v_switching[ idx ] , -M_l );
  v_OTS_KVL_lower[ idx ].set_lhs( -M_l );
  v_OTS_KVL_lower[ idx ].set_rhs( Inf< double >() );
  v_OTS_KVL_lower[ idx ].set_function(
                                   new LinearFunction( std::move( vars ) ) );
  }

 add_static_constraint( v_OTS_KVL_upper , "OTS_KVL_upper" );
 add_static_constraint( v_OTS_KVL_lower , "OTS_KVL_lower" );

 }  // end( OTSNetworkBlock::generate_OTS_Standard_constraints )

/*--------------------------------------------------------------------------*/

void OTSNetworkBlock::generate_OTS_Directional_constraints( void )
{
 // Directional BigM KVL relaxation for each DC line:
 //   F_l - B_l*Dtheta <= M_l*(1 - z+_l - z-_l)
 //   F_l - B_l*Dtheta >= -M_l*(1 - z+_l - z-_l)
 //
 // Rewritten as:
 //   F_l - B_l*theta_from + B_l*theta_to + M_l*z+_l + M_l*z-_l <= M_l
 //   F_l - B_l*theta_from + B_l*theta_to - M_l*z+_l - M_l*z-_l >= -M_l
 //
 // Plus exclusivity: z+_l + z-_l <= 1

 auto & DC_lines = f_NetworkData->get_DC_lines();
 auto n_dc = DC_lines.size();
 if( n_dc == 0 )
  return;

 const auto & start_line = f_NetworkData->get_start_line();
 const auto & end_line = f_NetworkData->get_end_line();
 LinearFunction::v_coeff_pair vars;

 v_OTS_KVL_upper.resize( n_dc );
 v_OTS_KVL_lower.resize( n_dc );
 v_switching_exclusivity.resize( n_dc );

 for( Index idx = 0 ; idx < n_dc ; ++idx ) {
  auto l = DC_lines[ idx ];
  double B_l = get_line_susceptance( l );
  double M_l = v_big_M[ idx ];

  // upper: F_l - B_l*theta_from + B_l*theta_to + M_l*z+ + M_l*z- <= M_l
  vars.emplace_back( & v_power_flow[ l ] , 1.0 );
  vars.emplace_back( & v_voltage_angle[ start_line[ l ] ] , -B_l );
  vars.emplace_back( & v_voltage_angle[ end_line[ l ] ] , B_l );
  vars.emplace_back( & v_switching_pos[ idx ] , M_l );
  vars.emplace_back( & v_switching_neg[ idx ] , M_l );
  v_OTS_KVL_upper[ idx ].set_lhs( -Inf< double >() );
  v_OTS_KVL_upper[ idx ].set_rhs( M_l );
  v_OTS_KVL_upper[ idx ].set_function(
                                   new LinearFunction( std::move( vars ) ) );

  // lower: F_l - B_l*theta_from + B_l*theta_to - M_l*z+ - M_l*z- >= -M_l
  vars.emplace_back( & v_power_flow[ l ] , 1.0 );
  vars.emplace_back( & v_voltage_angle[ start_line[ l ] ] , -B_l );
  vars.emplace_back( & v_voltage_angle[ end_line[ l ] ] , B_l );
  vars.emplace_back( & v_switching_pos[ idx ] , -M_l );
  vars.emplace_back( & v_switching_neg[ idx ] , -M_l );
  v_OTS_KVL_lower[ idx ].set_lhs( -M_l );
  v_OTS_KVL_lower[ idx ].set_rhs( Inf< double >() );
  v_OTS_KVL_lower[ idx ].set_function(
                                   new LinearFunction( std::move( vars ) ) );

  // exclusivity: z+ + z- <= 1
  vars.emplace_back( & v_switching_pos[ idx ] , 1.0 );
  vars.emplace_back( & v_switching_neg[ idx ] , 1.0 );
  v_switching_exclusivity[ idx ].set_lhs( -Inf< double >() );
  v_switching_exclusivity[ idx ].set_rhs( 1.0 );
  v_switching_exclusivity[ idx ].set_function(
                                   new LinearFunction( std::move( vars ) ) );
  }

 add_static_constraint( v_OTS_KVL_upper , "OTS_KVL_upper" );
 add_static_constraint( v_OTS_KVL_lower , "OTS_KVL_lower" );
 add_static_constraint( v_switching_exclusivity , "OTS_exclusivity" );

 }  // end( OTSNetworkBlock::generate_OTS_Directional_constraints )

/*--------------------------------------------------------------------------*/

void OTSNetworkBlock::generate_OTS_Elastic_constraints( void )
{
 // Elastic BigM KVL relaxation for each DC line:
 //   F_l - B_l*Dtheta <= alpha_l*M_l*(1-z_l) + beta_l*M_l*(1-z1_l)
 //   F_l - B_l*Dtheta >= -(alpha_l*M_l*(1-z_l) + beta_l*M_l*(1-z1_l))
 //
 // Rewritten as:
 //   F_l - B_l*th_from + B_l*th_to + aM*z_l + bM*z1_l <= aM + bM = M_l
 //   F_l - B_l*th_from + B_l*th_to - aM*z_l - bM*z1_l >= -(aM + bM) = -M_l
 //
 //   where aM = alpha_l*M_l, bM = beta_l*M_l
 //
 // Plus precedence: z_l <= z1_l  =>  z_l - z1_l <= 0

 auto & DC_lines = f_NetworkData->get_DC_lines();
 auto n_dc = DC_lines.size();
 if( n_dc == 0 )
  return;

 const auto & start_line = f_NetworkData->get_start_line();
 const auto & end_line = f_NetworkData->get_end_line();
 LinearFunction::v_coeff_pair vars;

 v_OTS_KVL_upper.resize( n_dc );
 v_OTS_KVL_lower.resize( n_dc );
 v_elastic_precedence.resize( n_dc );

 for( Index idx = 0 ; idx < n_dc ; ++idx ) {
  auto l = DC_lines[ idx ];
  double B_l = get_line_susceptance( l );
  double M_l = v_big_M[ idx ];
  double aM = v_alpha[ idx ] * M_l;
  double bM = v_beta[ idx ] * M_l;

  // upper: F_l - B_l*th_from + B_l*th_to + aM*z_l + bM*z1_l <= M_l
  vars.emplace_back( & v_power_flow[ l ] , 1.0 );
  vars.emplace_back( & v_voltage_angle[ start_line[ l ] ] , -B_l );
  vars.emplace_back( & v_voltage_angle[ end_line[ l ] ] , B_l );
  vars.emplace_back( & v_switching[ idx ] , aM );
  vars.emplace_back( & v_elastic[ idx ] , bM );
  v_OTS_KVL_upper[ idx ].set_lhs( -Inf< double >() );
  v_OTS_KVL_upper[ idx ].set_rhs( M_l );
  v_OTS_KVL_upper[ idx ].set_function(
                                   new LinearFunction( std::move( vars ) ) );

  // lower: F_l - B_l*th_from + B_l*th_to - aM*z_l - bM*z1_l >= -M_l
  vars.emplace_back( & v_power_flow[ l ] , 1.0 );
  vars.emplace_back( & v_voltage_angle[ start_line[ l ] ] , -B_l );
  vars.emplace_back( & v_voltage_angle[ end_line[ l ] ] , B_l );
  vars.emplace_back( & v_switching[ idx ] , -aM );
  vars.emplace_back( & v_elastic[ idx ] , -bM );
  v_OTS_KVL_lower[ idx ].set_lhs( -M_l );
  v_OTS_KVL_lower[ idx ].set_rhs( Inf< double >() );
  v_OTS_KVL_lower[ idx ].set_function(
                                   new LinearFunction( std::move( vars ) ) );

  // precedence: z_l - z1_l <= 0
  vars.emplace_back( & v_switching[ idx ] , 1.0 );
  vars.emplace_back( & v_elastic[ idx ] , -1.0 );
  v_elastic_precedence[ idx ].set_lhs( -Inf< double >() );
  v_elastic_precedence[ idx ].set_rhs( 0.0 );
  v_elastic_precedence[ idx ].set_function(
                                   new LinearFunction( std::move( vars ) ) );
  }

 add_static_constraint( v_OTS_KVL_upper , "OTS_KVL_upper" );
 add_static_constraint( v_OTS_KVL_lower , "OTS_KVL_lower" );
 add_static_constraint( v_elastic_precedence , "OTS_elastic_precedence" );

 }  // end( OTSNetworkBlock::generate_OTS_Elastic_constraints )

/*--------------------------------------------------------------------------*/

void OTSNetworkBlock::generate_OTS_ElasticDirectional_constraints( void )
{
 // Elastic Directional BigM KVL relaxation for each DC line:
 //   F_l - B_l*Dtheta <= aM*(1 - z+ - z-) + bM*(1 - z1)
 //   F_l - B_l*Dtheta >= -(aM*(1 - z+ - z-) + bM*(1 - z1))
 //
 // Rewritten as:
 //   F_l - B_l*th_from + B_l*th_to + aM*z+ + aM*z- + bM*z1 <= M_l
 //   F_l - B_l*th_from + B_l*th_to - aM*z+ - aM*z- - bM*z1 >= -M_l
 //
 // Plus exclusivity: z+ + z- <= 1
 // Plus precedence: z+ <= z1, z- <= z1

 auto & DC_lines = f_NetworkData->get_DC_lines();
 auto n_dc = DC_lines.size();
 if( n_dc == 0 )
  return;

 const auto & start_line = f_NetworkData->get_start_line();
 const auto & end_line = f_NetworkData->get_end_line();
 LinearFunction::v_coeff_pair vars;

 v_OTS_KVL_upper.resize( n_dc );
 v_OTS_KVL_lower.resize( n_dc );
 v_switching_exclusivity.resize( n_dc );
 // 2 precedence constraints per line: z+ <= z1, z- <= z1
 v_elastic_precedence.resize( 2 * n_dc );

 for( Index idx = 0 ; idx < n_dc ; ++idx ) {
  auto l = DC_lines[ idx ];
  double B_l = get_line_susceptance( l );
  double M_l = v_big_M[ idx ];
  double aM = v_alpha[ idx ] * M_l;
  double bM = v_beta[ idx ] * M_l;

  // upper: F_l - B_l*th_from + B_l*th_to + aM*z+ + aM*z- + bM*z1 <= M_l
  vars.emplace_back( & v_power_flow[ l ] , 1.0 );
  vars.emplace_back( & v_voltage_angle[ start_line[ l ] ] , -B_l );
  vars.emplace_back( & v_voltage_angle[ end_line[ l ] ] , B_l );
  vars.emplace_back( & v_switching_pos[ idx ] , aM );
  vars.emplace_back( & v_switching_neg[ idx ] , aM );
  vars.emplace_back( & v_elastic[ idx ] , bM );
  v_OTS_KVL_upper[ idx ].set_lhs( -Inf< double >() );
  v_OTS_KVL_upper[ idx ].set_rhs( M_l );
  v_OTS_KVL_upper[ idx ].set_function(
                                   new LinearFunction( std::move( vars ) ) );

  // lower: F_l - B_l*th_from + B_l*th_to - aM*z+ - aM*z- - bM*z1 >= -M_l
  vars.emplace_back( & v_power_flow[ l ] , 1.0 );
  vars.emplace_back( & v_voltage_angle[ start_line[ l ] ] , -B_l );
  vars.emplace_back( & v_voltage_angle[ end_line[ l ] ] , B_l );
  vars.emplace_back( & v_switching_pos[ idx ] , -aM );
  vars.emplace_back( & v_switching_neg[ idx ] , -aM );
  vars.emplace_back( & v_elastic[ idx ] , -bM );
  v_OTS_KVL_lower[ idx ].set_lhs( -M_l );
  v_OTS_KVL_lower[ idx ].set_rhs( Inf< double >() );
  v_OTS_KVL_lower[ idx ].set_function(
                                   new LinearFunction( std::move( vars ) ) );

  // exclusivity: z+ + z- <= 1
  vars.emplace_back( & v_switching_pos[ idx ] , 1.0 );
  vars.emplace_back( & v_switching_neg[ idx ] , 1.0 );
  v_switching_exclusivity[ idx ].set_lhs( -Inf< double >() );
  v_switching_exclusivity[ idx ].set_rhs( 1.0 );
  v_switching_exclusivity[ idx ].set_function(
                                   new LinearFunction( std::move( vars ) ) );

  // precedence 1: z+ - z1 <= 0
  vars.emplace_back( & v_switching_pos[ idx ] , 1.0 );
  vars.emplace_back( & v_elastic[ idx ] , -1.0 );
  v_elastic_precedence[ 2 * idx ].set_lhs( -Inf< double >() );
  v_elastic_precedence[ 2 * idx ].set_rhs( 0.0 );
  v_elastic_precedence[ 2 * idx ].set_function(
                                   new LinearFunction( std::move( vars ) ) );

  // precedence 2: z- - z1 <= 0
  vars.emplace_back( & v_switching_neg[ idx ] , 1.0 );
  vars.emplace_back( & v_elastic[ idx ] , -1.0 );
  v_elastic_precedence[ 2 * idx + 1 ].set_lhs( -Inf< double >() );
  v_elastic_precedence[ 2 * idx + 1 ].set_rhs( 0.0 );
  v_elastic_precedence[ 2 * idx + 1 ].set_function(
                                   new LinearFunction( std::move( vars ) ) );
  }

 add_static_constraint( v_OTS_KVL_upper , "OTS_KVL_upper" );
 add_static_constraint( v_OTS_KVL_lower , "OTS_KVL_lower" );
 add_static_constraint( v_switching_exclusivity , "OTS_exclusivity" );
 add_static_constraint( v_elastic_precedence , "OTS_elastic_precedence" );

 }  // end( OTSNetworkBlock::generate_OTS_ElasticDirectional_constraints )

/*--------------------------------------------------------------------------*/

void OTSNetworkBlock::generate_OTS_flow_bounds( void )
{
 // Generate flow bounds coupled with switching variables.
 // For Standard/Elastic:
 //   -f_max_l * z_l <= F_l <= f_max_l * z_l
 //   i.e.  F_l - f_max*z <= 0   and   -F_l - f_max*z <= 0
 //
 // For Directional/ElasticDirectional:
 //   F_l <= f_max_l * z+  and  F_l >= -f_max_l * z-
 //   i.e.  F_l - f_max*z+ <= 0  and  -F_l - f_max*z- <= 0

 auto & DC_lines = f_NetworkData->get_DC_lines();
 auto n_dc = DC_lines.size();
 if( n_dc == 0 )
  return;

 LinearFunction::v_coeff_pair vars;

 v_OTS_flow_upper.resize( n_dc );
 v_OTS_flow_lower.resize( n_dc );

 for( Index idx = 0 ; idx < n_dc ; ++idx ) {
  auto l = DC_lines[ idx ];
  double f_max = f_NetworkData->get_max_power_flow( l );

  switch( f_ots_type ) {
   case kOTS_Standard:
   case kOTS_Elastic: {
    // upper: F_l - f_max * z_l <= 0
    vars.emplace_back( & v_power_flow[ l ] , 1.0 );
    vars.emplace_back( & v_switching[ idx ] , -f_max );
    v_OTS_flow_upper[ idx ].set_lhs( -Inf< double >() );
    v_OTS_flow_upper[ idx ].set_rhs( 0.0 );
    v_OTS_flow_upper[ idx ].set_function(
                                   new LinearFunction( std::move( vars ) ) );

    // lower: -F_l - f_max * z_l <= 0  =>  F_l + f_max*z >= 0
    vars.emplace_back( & v_power_flow[ l ] , 1.0 );
    vars.emplace_back( & v_switching[ idx ] , f_max );
    v_OTS_flow_lower[ idx ].set_lhs( 0.0 );
    v_OTS_flow_lower[ idx ].set_rhs( Inf< double >() );
    v_OTS_flow_lower[ idx ].set_function(
                                   new LinearFunction( std::move( vars ) ) );
    break;
    }
   case kOTS_Directional:
   case kOTS_ElasticDirectional: {
    // upper: F_l - f_max * z+ <= 0
    vars.emplace_back( & v_power_flow[ l ] , 1.0 );
    vars.emplace_back( & v_switching_pos[ idx ] , -f_max );
    v_OTS_flow_upper[ idx ].set_lhs( -Inf< double >() );
    v_OTS_flow_upper[ idx ].set_rhs( 0.0 );
    v_OTS_flow_upper[ idx ].set_function(
                                   new LinearFunction( std::move( vars ) ) );

    // lower: F_l + f_max * z- >= 0
    vars.emplace_back( & v_power_flow[ l ] , 1.0 );
    vars.emplace_back( & v_switching_neg[ idx ] , f_max );
    v_OTS_flow_lower[ idx ].set_lhs( 0.0 );
    v_OTS_flow_lower[ idx ].set_rhs( Inf< double >() );
    v_OTS_flow_lower[ idx ].set_function(
                                   new LinearFunction( std::move( vars ) ) );
    break;
    }
   }
  }

 add_static_constraint( v_OTS_flow_upper , "OTS_flow_upper" );
 add_static_constraint( v_OTS_flow_lower , "OTS_flow_lower" );

 }  // end( OTSNetworkBlock::generate_OTS_flow_bounds )

/*--------------------------------------------------------------------------*/

void OTSNetworkBlock::generate_design_coupling_constraints( void )
{
 // If design variables exist, couple switching with design:
 //   z_l <= x_l  (Standard/Elastic)
 //   z+_l + z-_l <= x_l  (Directional/ElasticDirectional)

 if( ! is_design() )
  return;

 auto & DC_lines = f_NetworkData->get_DC_lines();
 auto n_dc = DC_lines.size();

 // count how many DC lines actually have a design variable
 Index n_coupling = 0;
 for( Index idx = 0 ; idx < n_dc ; ++idx ) {
  auto l = DC_lines[ idx ];
  if( get_design( l ) )
   ++n_coupling;
  }

 if( n_coupling == 0 )
  return;

 v_design_coupling.resize( n_coupling );
 LinearFunction::v_coeff_pair vars;

 Index c_idx = 0;
 for( Index idx = 0 ; idx < n_dc ; ++idx ) {
  auto l = DC_lines[ idx ];
  ColVariable * x_l = get_design( l );
  if( ! x_l )
   continue;

  switch( f_ots_type ) {
   case kOTS_Standard:
   case kOTS_Elastic:
    // z_l - x_l <= 0
    vars.emplace_back( & v_switching[ idx ] , 1.0 );
    vars.emplace_back( x_l , -1.0 );
    break;
   case kOTS_Directional:
   case kOTS_ElasticDirectional:
    // z+_l + z-_l - x_l <= 0
    vars.emplace_back( & v_switching_pos[ idx ] , 1.0 );
    vars.emplace_back( & v_switching_neg[ idx ] , 1.0 );
    vars.emplace_back( x_l , -1.0 );
    break;
   }

  v_design_coupling[ c_idx ].set_lhs( -Inf< double >() );
  v_design_coupling[ c_idx ].set_rhs( 0.0 );
  v_design_coupling[ c_idx ].set_function(
                                   new LinearFunction( std::move( vars ) ) );
  ++c_idx;
  }

 add_static_constraint( v_design_coupling , "OTS_design_coupling" );

 }  // end( OTSNetworkBlock::generate_design_coupling_constraints )

/*--------------------------------------------------------------------------*/
/*--------------------- End File OTSNetworkBlock.cpp -----------------------*/
/*--------------------------------------------------------------------------*/
