/*--------------------------------------------------------------------------*/
/*--------------------- File DCNetworkBlock.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the DCNetworkBlock class.
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

#include <map>

#include "NetworkBlock.h"

#include "DCNetworkBlock.h"

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

// register DCNetworkBlock to the Block factory

SMSpp_insert_in_factory_cpp_1( DCNetworkBlock );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

// register DCNetworkData to the NetworkData factory

typedef DCNetworkBlock::DCNetworkData DCNetworkData;

SMSpp_insert_in_factory_cpp_1( DCNetworkData );

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS OF DCNetworkBlock ------------------------*/
/*--------------------------------------------------------------------------*/

DCNetworkBlock::~DCNetworkBlock()
{
 Constraint::clear( v_AC_power_flow_limit_const );
 Constraint::clear( v_AC_HVDC_power_flow_limit_const );
 Constraint::clear( v_power_flow_injection_const );
 Constraint::clear( v_power_flow_relax_abs );

 Constraint::clear( v_HVDC_power_flow_limit_const );
 Constraint::clear( node_injection_bounds_const );

 objective.clear();

 // Delete the DCNetworkData if it is local.
 if( f_local_NetworkData )
  delete( f_NetworkData );
}

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void DCNetworkData::deserialize( const netCDF::NcGroup & group )
{

#ifndef NDEBUG
 static std::vector< std::string > expected_dims = { "NumberNodes" ,
                                                     "NumberLines" ,
                                                     // if called from UCBlock:
                                                     "TimeHorizon" ,
                                                     "NumberUnits" ,
                                                     "NumberNetworks" ,
                                                     "NumberElectricalGenerators" };
 check_dimensions( group , expected_dims , std::cerr );

 static std::vector< std::string > expected_vars = { "ActiveDemand" ,
                                                     "StartLine" ,
                                                     "EndLine" ,
                                                     "MinPowerFlow" ,
                                                     "MaxPowerFlow" ,
                                                     "Susceptance" ,
                                                     "NetworkCost" ,
                                                     "NodeName" ,
                                                     "LineName" ,
                                                     // if called from UCBlock:
                                                     "ActivePowerDemand" ,
                                                     "GeneratorNode" ,
                                                     "NetworkConstantTerms" ,
                                                     "NetworkBlockClassname" ,
                                                     "NetworkDataClassname" };
 check_variables( group , expected_vars , std::cerr );
#endif

 // Optional variables

 if( ! ::deserialize_dim( group , "NumberNodes" , f_number_nodes ) )
  f_number_nodes = 1;

 if( f_number_nodes > 1 ) {

  ::deserialize_dim( group , "NumberLines" , f_number_lines , false );

  ::deserialize( group , "StartLine" , f_number_lines , v_start_line , false ,
                 true );

  ::deserialize( group , "EndLine" , f_number_lines , v_end_line , false ,
                 true );

  ::deserialize( group , "MinPowerFlow" , f_number_lines , v_min_power_flow ,
                 true , true );

  ::deserialize( group , "MaxPowerFlow" , f_number_lines , v_max_power_flow ,
                 true , true );

  ::deserialize( group , "Susceptance" , f_number_lines , v_susceptance ,
                 true , true );

  ::deserialize( group , "NetworkCost" , f_number_lines , v_network_cost ,
                 true , true );
 }

 const auto get_string_array =
  [ &group ]( const std::string & var_name ,
              std::vector< std::string > & v_string ,
              Index size = Inf< Index >() ) {
   v_string.clear();
   auto netcdf_var = group.getVar( var_name );
   if( ! netcdf_var.isNull() ) {
    if( netcdf_var.getDimCount() != 1 )
     throw( std::logic_error( "DCNetworkData::deserialize: the dimension of "
                              "variable'" + var_name + "' must be 1." ) );

    if( ( size < Inf< Index >() ) &&
        ( netcdf_var.getDim( 0 ).getSize() != size ) )
     throw( std::logic_error( "DCNetworkData::deserialize: the size of "
                              "variable '" + var_name + "' should be " +
                              std::to_string( size ) + "." ) );

    const auto var_size = netcdf_var.getDim( 0 ).getSize();
    v_string.reserve( var_size );

    // TODO The following implementation should change when netCDF provides a
    // better C++ interface.

    for( Index i = 0 ; i < var_size ; ++i ) {
     char * fname = nullptr;
     netcdf_var.getVar( { i } , { 1 } , &fname );
     v_string.push_back( fname );
     free( fname );
    }
   }
  };

 get_string_array( "NodeName" , v_node_names , f_number_nodes );
 get_string_array( "LineName" , v_line_names );

}  // end( DCNetworkData::deserialize )

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::deserialize( const netCDF::NcGroup & group )
{

#ifndef NDEBUG
 static std::vector< std::string > expected_dims = { "NumberNodes" };
 check_dimensions( group , expected_dims , std::cerr );

 static std::vector< std::string > expected_vars = { "ActiveDemand" ,
                                                     "ConstantTerm" };
 check_variables( group , expected_vars , std::cerr );
#endif

 // Optional variables

 Index NumberNodes;
 if( ::deserialize_dim( group , "NumberNodes" , NumberNodes ) ) {
  // Since the dimension "NumberNodes" has been provided, it means that a
  // DCNetworkData has been provided. Thus, the DCNetworkData is deserialized,
  // and it is marked as being local.
  delete( f_NetworkData );
  f_NetworkData = new DCNetworkData();
  f_NetworkData->deserialize( group );
  f_local_NetworkData = true;
  // A DCNetworkData has been provided. So, the size of the given vector of
  // active demand must be equal to the number of nodes.
  ::deserialize( group , "ActiveDemand" , NumberNodes , v_ActiveDemand );
 } else {
  // A DCNetworkData has not been provided. However, the active demand may still
  // have been provided.

  auto ActiveDemand = group.getVar( "ActiveDemand" );

  if( ! ActiveDemand.isNull() ) {
   // The active demand has indeed been provided.

   if( ActiveDemand.getDimCount() != 1 )
    // The active demand must be a one-dimensional array.
    throw( std::invalid_argument(
     "DCNetworkBlock::deserialize(): ActiveDemand should have one dimension, "
     "but it has " + std::to_string( ActiveDemand.getDimCount() ) ) );

   // Retrieve the number of nodes from the size of the given netCDF variable.
   const auto number_nodes = ActiveDemand.getDim( 0 ).getSize();

   // Resize the vector of active demand.
   v_ActiveDemand.resize( number_nodes );

   // Retrieve the active demand from the netCDF variable.
   ActiveDemand.getVar( v_ActiveDemand.data() );
  }
 }

 ::deserialize( group , f_ConstTerm , "ConstantTerm" );

}  // end( DCNetworkBlock::deserialize )

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::generate_abstract_variables( Configuration * stvv )
{
 if( variables_generated() )  // variables have already been generated
  return;                     // nothing to do

 NetworkBlock::generate_abstract_variables( stvv );

 const auto number_lines = get_number_lines();

 if( number_lines > 0 ) {
  // the power flow Variable
  v_power_flow.resize( number_lines );
  for( auto & var : v_power_flow )
   var.set_type( ColVariable::kContinuous );
  add_static_variable( v_power_flow , "p_flow_network" );

  if( ! f_NetworkData->get_network_cost().empty() ) {
   // the auxiliary Variable
   v_auxiliary_variable.resize( number_lines );
   for( auto & var : v_auxiliary_variable )
    var.set_type( ColVariable::kContinuous );
   add_static_variable( v_auxiliary_variable , "aux_network" );
  }
 }

 set_variables_generated();

}  // end( DCNetworkBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::generate_abstract_constraints( Configuration * stcc )
{
 if( constraints_generated() )  // constraints have already been generated
  return;                       // nothing to do

 const auto number_nodes = get_number_nodes();

 if( number_nodes <= 1 )
  return;

 const auto number_lines = get_number_lines();

 if( number_lines <= 0 )
  throw( std::logic_error( "DCNetworkBlock::generate_abstract_constraints: "
                           "number of lines of DCNetworkBlock is not set" ) );

 const auto & start_line = f_NetworkData->get_start_line();
 const auto & end_line = f_NetworkData->get_end_line();
 const auto lines_type = f_NetworkData->get_lines_type();

 LinearFunction::v_coeff_pair vars;

 if( lines_type == kHVDC ) {  // HVDC power flow limit

  // HVDC power flow limit constraints

  v_HVDC_power_flow_limit_const.resize( number_lines );

  for( Index line_id = 0 ; line_id < number_lines ; ++line_id ) {

   const auto kappa = get_kappa( line_id );

   v_HVDC_power_flow_limit_const[ line_id ].set_lhs(
    kappa * get_min_power_flow( line_id ) );
   v_HVDC_power_flow_limit_const[ line_id ].set_rhs(
    kappa * get_max_power_flow( line_id ) );
   v_HVDC_power_flow_limit_const[ line_id ].set_variable(
    &v_power_flow[ line_id ] );
  }

  add_static_constraint( v_HVDC_power_flow_limit_const ,
                         "HVDC_power_flow_limit" );

  // HVDC power flow and node injection constraints

  v_power_flow_injection_const.resize( number_nodes );

  for( Index n = 0 ; n < number_nodes ; ++n ) {

   vars.push_back( std::make_pair( &v_node_injection[ 0 ][ n ] , -1.0 ) );

   for( Index line_id = 0 ; line_id < number_lines ; ++line_id ) {

    if( start_line[ line_id ] == n )
     vars.push_back( std::make_pair( &v_power_flow[ line_id ] , 1.0 ) );
    if( end_line[ line_id ] == n )
     vars.push_back( std::make_pair( &v_power_flow[ line_id ] , -1.0 ) );
   }

   v_power_flow_injection_const[ n ].set_both( -v_ActiveDemand[ n ] );
   v_power_flow_injection_const[ n ].set_function(
    new LinearFunction( std::move( vars ) ) );
  }

  add_static_constraint( v_power_flow_injection_const ,
                         "HVDC_power_flow_injection" );

  if( ! f_NetworkData->get_network_cost().empty() ) {

   // Auxiliary constraints for relaxing the absolute value

   v_power_flow_relax_abs.resize(
    boost::multi_array< FRowConstraint , 2 >::extent_gen()[ 2 ][ number_lines ] );

   for( Index line_id = 0 ; line_id < number_lines ; ++line_id ) {

    // F_l <= V_l

    vars.push_back( std::make_pair( &v_power_flow[ line_id ] , -1.0 ) );
    vars.push_back( std::make_pair( &v_auxiliary_variable[ line_id ] , 1.0 ) );

    v_power_flow_relax_abs[ 0 ][ line_id ].set_lhs( 0.0 );
    v_power_flow_relax_abs[ 0 ][ line_id ].set_rhs( Inf< double >() );
    v_power_flow_relax_abs[ 0 ][ line_id ].set_function(
     new LinearFunction( std::move( vars ) ) );

    // - F_l <= V_l

    vars.push_back( std::make_pair( &v_power_flow[ line_id ] , 1.0 ) );
    vars.push_back( std::make_pair( &v_auxiliary_variable[ line_id ] , 1.0 ) );

    v_power_flow_relax_abs[ 1 ][ line_id ].set_lhs( 0.0 );
    v_power_flow_relax_abs[ 1 ][ line_id ].set_rhs( Inf< double >() );
    v_power_flow_relax_abs[ 1 ][ line_id ].set_function(
     new LinearFunction( std::move( vars ) ) );
   }

   add_static_constraint( v_power_flow_relax_abs , "power_flow_relax_abs" );
  }  // end( cost not empty )
 }  // end( HVDC_Lines constraints )

/*--------------------------------------------------------------------------*/

 else if( lines_type == kAC ) {  // AC power flow limit

  // AC power flow limit constraints

  /* v_AC_power_flow_limit_const.resize( get_number_lines() );

  for( Index line_id = 0 ; line_id < get_number_lines() ; ++line_id ) {

   const auto kappa = get_kappa( line_id );

   double constant_term = 0;

   for( Index node_id = 0 ; node_id < get_number_nodes() ; ++node_id ) {

    double coefficient = 0.0;

    // Power Transfer Distribution Factor matrix (PTDF-matrix)
    vars.push_back( std::make_pair( &v_node_injection[ 0 ][ node_id ] ,
                                    coefficient ) );
    constant_term -= coefficient * v_ActiveDemand[ node_id ];
   }

   v_AC_power_flow_limit_const[ line_id ].set_lhs(
    kappa * get_min_power_flow( line_id ) - constant_term );
   v_AC_power_flow_limit_const[ line_id ].set_rhs(
    kappa * get_max_power_flow( line_id ) - constant_term );
   v_AC_power_flow_limit_const[ line_id ].set_function(
    new LinearFunction( std::move( vars ) ) );
  }

  add_static_constraint( v_AC_power_flow_limit_const ,
                         "AC_power_flow_limit" ); */

 }  // end( AC_Lines constraints )

/*--------------------------------------------------------------------------*/

 else if( lines_type == kAC_HVDC ) {  // AC/HVDC power flow limit

  // AC/HVDC power flow limit constraints

  /* v_AC_HVDC_power_flow_limit_const.resize( get_number_lines() );


  add_static_constraint( v_AC_HVDC_power_flow_limit_const ,
                         "AC_HVDC_power_flow_limit" ); */

 }  // end( AC/HVDC constraints)

 // node injection bound constraints

 node_injection_bounds_const.resize( number_nodes );

 for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {

  node_injection_bounds_const[ node_id ].set_lhs(
   v_MinNodeInjection[ 0 ][ node_id ] );
  node_injection_bounds_const[ node_id ].set_rhs(
   v_MaxNodeInjection[ 0 ][ node_id ] );
  node_injection_bounds_const[ node_id ].set_variable(
   &v_node_injection[ 0 ][ node_id ] );
 }

 add_static_constraint( node_injection_bounds_const ,
                        "Node_Injection_Bound_Const_Network" );

 set_constraints_generated();

}  // end( DCNetworkBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::generate_objective( Configuration * objc )
{
 if( objective_generated() )  // Objective has already been generated
  return;                     // nothing to do

 auto lf = new LinearFunction();

 if( ! f_NetworkData->get_network_cost().empty() )
  for( Index line_id = 0 ; line_id < get_number_lines() ; ++line_id )
   lf->add_variable( &v_auxiliary_variable[ line_id ] ,
                     f_NetworkData->get_network_cost()[ line_id ] ,
                     eDryRun );

 lf->set_constant_term( f_ConstTerm );

 objective.set_function( lf );
 objective.set_sense( Objective::eMin );

 // Set Block objective
 this->set_objective( &objective );

 set_objective_generated();

}  // end( DCNetworkBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*----------------- METHODS FOR CHECKING THE DCNetworkBlock ----------------*/
/*--------------------------------------------------------------------------*/

bool DCNetworkBlock::is_feasible( bool useabstract , Configuration * fsbc )
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
  && ColVariable::is_feasible( v_power_flow , tol )
  && ColVariable::is_feasible( v_auxiliary_variable , tol )
  // Constraints
  && RowConstraint::is_feasible( v_AC_power_flow_limit_const , tol , rel_viol )
  && RowConstraint::is_feasible( v_AC_HVDC_power_flow_limit_const , tol , rel_viol )
  && RowConstraint::is_feasible( v_power_flow_injection_const , tol , rel_viol )
  && RowConstraint::is_feasible( v_power_flow_relax_abs , tol , rel_viol )
  && RowConstraint::is_feasible( v_HVDC_power_flow_limit_const , tol , rel_viol )
  && RowConstraint::is_feasible( node_injection_bounds_const , tol , rel_viol ) );

} // end( DCNetworkBlock::is_feasible )

/*--------------------------------------------------------------------------*/
/*--------- METHODS FOR LOADING, PRINTING & SAVING THE DCNetworkBlock ------*/
/*--------------------------------------------------------------------------*/

void DCNetworkData::serialize( netCDF::NcGroup & group ) const {

 auto NumberNodes = group.addDim( "NumberNodes" , f_number_nodes );

 if( f_number_nodes > 1 ) {
  auto NumberLines = group.addDim( "NumberLines" );

  ::serialize( group , "StartLine" , netCDF::NcUint() , NumberLines ,
               v_start_line );

  ::serialize( group , "EndLine" , netCDF::NcUint() , NumberLines ,
               v_end_line );

  ::serialize( group , "MinPowerFlow" , netCDF::NcDouble() , NumberLines ,
               v_min_power_flow );

  ::serialize( group , "MaxPowerFlow" , netCDF::NcDouble() , NumberLines ,
               v_max_power_flow );

  ::serialize( group , "Susceptance" , netCDF::NcDouble() , NumberLines ,
               v_susceptance );

  ::serialize( group , "NetworkCost" , netCDF::NcDouble() , NumberLines ,
               v_network_cost );

  if( ! v_line_names.empty() ) {
   assert( v_line_names.size() == NumberLines.getSize() );
   auto LineName = group.addVar( "LineName" , netCDF::NcString() , NumberLines );
   for( Index i = 0 ; i < v_line_names.size() ; ++i )
    LineName.putVar( { i } , v_line_names[ i ] );
  }
 }

 if( ! v_node_names.empty() ) {
  assert( v_node_names.size() == NumberNodes.getSize() );
  auto NodeName = group.addVar( "NodeName" , netCDF::NcString() , NumberNodes );
  for( Index i = 0 ; i < v_node_names.size() ; ++i )
   NodeName.putVar( { i } , v_node_names[ i ] );
 }

}  // end( DCNetworkData::serialize )

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::serialize( netCDF::NcGroup & group ) const
{
 NetworkBlock::serialize( group );

 if( auto network_data = get_NetworkData() )
  // If a DCNetworkData is present, serialize it.
  network_data->serialize( group );

 auto NumberNodes = group.getDim( "NumberNodes" );

 if( ! v_ActiveDemand.empty() ) {
  // This DCNetworkBlock has active demand, so it is serialized.

  if( NumberNodes.isNull() )
   /* The dimension "NumberNodes" is not present in the group (which means
    * that a DCNetworkData is not present). However, the number of nodes can
    * still be obtained from the size of the active demand vector. Notice that
    * the name "NumberNodes" is not used for this new dimension, because it
    * would indicate that a DCNetworkData is present (which is not the
    * case). Therefore, we create an alternative dimension in order to be able
    * to serialize the active demand. */
   NumberNodes = group.addDim( "__NumberNodes__" , v_ActiveDemand.size() );

  // Finally, serialize the active demand.
  ::serialize( group , "ActiveDemand" , netCDF::NcDouble() ,
               NumberNodes , v_ActiveDemand );
 }

 if( f_ConstTerm != 0 )
  ::serialize( group , "ConstantTerm" , netCDF::NcDouble() , f_ConstTerm );

}  // end( DCNetworkBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

void DCNetworkBlock::set_active_demand( MF_dbl_it values ,
                                        Block::Subset && subset ,
                                        const bool ordered ,
                                        c_ModParam issuePMod ,
                                        c_ModParam issueAMod )
{
 if( subset.empty() )
  return;

 if( v_ActiveDemand.empty() ) {
  if( std::all_of( values ,
                   values + subset.size() ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  v_ActiveDemand.resize( get_number_nodes() );
 }

 bool identical = true;
 for( auto i : subset ) {
  if( i >= v_ActiveDemand.size() )
   throw( std::invalid_argument( "DCNetworkBlock::set_active_demand: "
                                 "invalid value in subset." ) );

  auto demand = *(values++);
  if( v_ActiveDemand[ i ] != demand ) {
   identical = false;

   if( not_dry_run( issuePMod ) )
    // Change the physical representation
    v_ActiveDemand[ i ] = demand;
  }
 }

 if( identical )
  return;  // nothing changes; return

 if( not_dry_run( issuePMod ) &&
     not_dry_run( issueAMod ) &&
     constraints_generated() ) {
  // Change the abstract representation

  switch( f_NetworkData->get_lines_type() ) {
   case( kHVDC ): {
    for( auto i : subset )
     v_power_flow_injection_const[ i ].set_both( -v_ActiveDemand[ i ] ,
                                                 issueAMod );
    break;
   }
   case( kAC ):
    // TODO
    break;
   case( kAC_HVDC ):
    // TODO
    break;
   default:
    break;
  }
 }

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification
  if( ! ordered )
   std::sort( subset.begin() , subset.end() );

  Block::add_Modification( std::make_shared< NetworkBlockSbstMod >(
                            this , NetworkBlockMod::eSetActD , std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );
 }
}  // end( DCNetworkData::set_active_demand )

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::set_active_demand( MF_dbl_it values ,
                                        Block::Range rng ,
                                        c_ModParam issuePMod ,
                                        c_ModParam issueAMod )
{
 rng.second = std::min( rng.second , get_number_nodes() );
 if( rng.second <= rng.first )
  return;

 if( v_ActiveDemand.empty() ) {
  if( std::all_of( values ,
                   values + ( rng.second - rng.first ) ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  v_ActiveDemand.resize( get_number_nodes() );
 }

 // If nothing changes, return
 if( std::equal( values , values + ( rng.second - rng.first ) ,
                 v_ActiveDemand.begin() + rng.first ) )
  return;

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation
  std::copy( values ,
             values + ( rng.second - rng.first ) ,
             v_ActiveDemand.begin() + rng.first );

  if( not_dry_run( issueAMod ) && constraints_generated() ) {
   // Change the abstract representation
   switch( f_NetworkData->get_lines_type() ) {
    case( kHVDC ): {
     for( Index i = rng.first ; i < rng.second ; ++i )
      v_power_flow_injection_const[ i ].set_both( -v_ActiveDemand[ i ] ,
                                                        issueAMod );
     break;
    }
    case( kAC ):
     // TODO
     break;
    case( kAC_HVDC ):
     // TODO
     break;
    default:
     break;
   }
  }
 }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< NetworkBlockRngdMod >(
                            this , NetworkBlockMod::eSetActD , rng ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( DCNetworkData::set_active_demand )

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::set_kappa( MF_dbl_it values ,
                                Block::Subset && subset ,
                                const bool ordered ,
                                c_ModParam issuePMod ,
                                c_ModParam issueAMod )
{
 if( subset.empty() )
  return;

 if( v_kappa.empty() ) {
  if( std::all_of( values ,
                   values + subset.size() ,
                   []( double cst ) { return( cst == 1 ); } ) )
   return;

  v_kappa.resize( get_number_lines() , 1 );
 }

 bool identical = true;
 for( auto i : subset ) {
  if( i >= v_kappa.size() )
   throw( std::invalid_argument( "DCNetworkBlock::set_kappa: invalid value in"
                                 " subset: " + std::to_string( i ) + "." ) );

  const auto kappa = *(values++);
  if( v_kappa[ i ] != kappa ) {
   identical = false;
   if( not_dry_run( issuePMod ) )
    // Change the physical representation
    v_kappa[ i ] = kappa;
  }
 }

 if( identical )
  return;  // nothing changes; return

 if( not_dry_run( issuePMod ) &&
     not_dry_run( issueAMod ) &&
     constraints_generated() ) {

  // Change the abstract representation

  switch( f_NetworkData->get_lines_type() ) {
   case( kHVDC ): {
    for( auto i : subset ) {
     v_HVDC_power_flow_limit_const[ i ].set_lhs(
      v_kappa[ i ] * get_min_power_flow( i ) , issueAMod );
     v_HVDC_power_flow_limit_const[ i ].set_rhs(
      v_kappa[ i ] * get_max_power_flow( i ) , issueAMod );
    }
    break;
   }
   case( kAC ): {
    for( auto i : subset ) {
     v_AC_power_flow_limit_const[ i ].set_lhs(
      v_kappa[ i ] * get_min_power_flow( i ) , issueAMod );
     v_AC_power_flow_limit_const[ i ].set_rhs(
      v_kappa[ i ] * get_max_power_flow( i ) , issueAMod );
    }
    break;
   }
   case( kAC_HVDC ):
    // TODO
    break;
   default: break;
  }
 }

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification
  if( ! ordered )
   std::sort( subset.begin() , subset.end() );

  Block::add_Modification( std::make_shared< DCNetworkBlockSbstMod >(
                            this , DCNetworkBlockMod::eSetKappa , std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );
 }
}  // end( DCNetworkData::set_kappa( subset ) )

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::set_kappa( MF_dbl_it values ,
                                Block::Range rng ,
                                c_ModParam issuePMod ,
                                c_ModParam issueAMod )
{
 rng.second = std::min( rng.second , get_number_lines() );
 if( rng.second <= rng.first )
  return;

 if( v_kappa.empty() ) {
  if( std::all_of( values ,
                   values + ( rng.second - rng.first ) ,
                   []( double cst ) { return( cst == 1 ); } ) )
   return;

  v_kappa.resize( get_number_lines() , 1 );
 }

 // If nothing changes, return
 if( std::equal( values , values + ( rng.second - rng.first ) ,
                 v_kappa.begin() + rng.first ) )
  return;

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation
  std::copy( values ,
             values + ( rng.second - rng.first ) ,
             v_kappa.begin() + rng.first );

  if( not_dry_run( issueAMod ) && constraints_generated() ) {
   // Change the abstract representation

   switch( f_NetworkData->get_lines_type() ) {
    case( kHVDC ): {
     for( Index i = rng.first ; i < rng.second ; ++i ) {
      v_HVDC_power_flow_limit_const[ i ].set_lhs(
       v_kappa[ i ] * get_min_power_flow( i ) , issueAMod );
      v_HVDC_power_flow_limit_const[ i ].set_rhs(
       v_kappa[ i ] * get_max_power_flow( i ) , issueAMod );
     }
     break;
    }
    case( kAC ): {
     for( Index i = rng.first ; i < rng.second ; ++i ) {
      v_AC_power_flow_limit_const[ i ].set_lhs(
       v_kappa[ i ] * get_min_power_flow( i ) , issueAMod );
      v_AC_power_flow_limit_const[ i ].set_rhs(
       v_kappa[ i ] * get_max_power_flow( i ) , issueAMod );
     }
     break;
    }
    case( kAC_HVDC ):
     // TODO
     break;
    default: break;
   }
  }
 }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< DCNetworkBlockRngdMod >(
                            this , DCNetworkBlockMod::eSetKappa , rng ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( DCNetworkData::set_kappa( range ) )

/*--------------------------------------------------------------------------*/
/*--------------------- End File DCNetworkBlock.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
