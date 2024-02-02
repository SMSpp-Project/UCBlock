#include <map>
#include "NetworkBlock.h"
#include "ACNetworkBlock.h"
#include "LinearFunction.h"
#include "OneVarConstraint.h"
#include "FRealObjective.h"

using namespace SMSpp_di_unipi_it;

SMSpp_insert_in_factory_cpp_1( ACNetworkBlock );

typedef ACNetworkBlock::ACNetworkData ACNetworkData;

SMSpp_insert_in_factory_cpp_1( ACNetworkData );

ACNetworkBlock::~ACNetworkBlock()
{
 Constraint::clear( v_AC_power_flow_limit_const );
 Constraint::clear( v_AC_HVDC_power_flow_limit_const );
 Constraint::clear( v_power_flow_injection_const );
 Constraint::clear( v_power_flow_relax_abs );

 Constraint::clear( v_HVDC_power_flow_limit_const );
 Constraint::clear( node_injection_bounds_const );

 objective.clear();

 // Delete the ACNetworkData if it is local.
 if( f_local_NetworkData )
  delete( f_NetworkData );
}

void ACNetworkData::deserialize( const netCDF::NcGroup & group )
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

  if( ! ::deserialize_dim( group, "ReferenceNode", f_reference_node, true ) )
    f_reference_node = 0;
 }

 const auto get_string_array =
  [ &group ]( const std::string & var_name ,
              std::vector< std::string > & v_string ,
              Index size = Inf< Index >() ) {
   v_string.clear();
   auto netcdf_var = group.getVar( var_name );
   if( ! netcdf_var.isNull() ) {
    if( netcdf_var.getDimCount() != 1 )
     throw( std::logic_error( "ACNetworkData::deserialize: the dimension of "
                              "variable'" + var_name + "' must be 1." ) );

    if( ( size < Inf< Index >() ) &&
        ( netcdf_var.getDim( 0 ).getSize() != size ) )
     throw( std::logic_error( "ACNetworkData::deserialize: the size of "
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

}  // end( ACNetworkData::deserialize )

void ACNetworkBlock::run_tree_decomposition(){
  std::cout << "coucou" << std::endl;
}

/*--------------------------------------------------------------------------*/

void ACNetworkBlock::deserialize( const netCDF::NcGroup & group )
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
  // ACNetworkData has been provided. Thus, the ACNetworkData is deserialized,
  // and it is marked as being local.
  delete( f_NetworkData );
  f_NetworkData = new ACNetworkData();
  f_NetworkData->deserialize( group );
  f_local_NetworkData = true;
  // A ACNetworkData has been provided. So, the size of the given vector of
  // active demand must be equal to the number of nodes.
  ::deserialize( group , "ActiveDemand" , NumberNodes , v_ActiveDemand );
 } else {
  // A ACNetworkData has not been provided. However, the active demand may still
  // have been provided.

  auto ActiveDemand = group.getVar( "ActiveDemand" );

  if( ! ActiveDemand.isNull() ) {
   // The active demand has indeed been provided.

   if( ActiveDemand.getDimCount() != 1 )
    // The active demand must be a one-dimensional array.
    throw( std::invalid_argument(
     "ACNetworkBlock::deserialize(): ActiveDemand should have one dimension, "
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

}  // end( ACNetworkBlock::deserialize )


/*--------------------------------------------------------------------------*/

void ACNetworkBlock::generate_abstract_variables( Configuration * stvv )
{

}  // end( ACNetworkBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void ACNetworkBlock::generate_abstract_constraints( Configuration * stcc )
{

}  // end( ACNetworkBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void ACNetworkBlock::generate_objective( Configuration * objc )
{

}  // end( ACNetworkBlock::generate_objective )

bool ACNetworkBlock::is_feasible( bool useabstract , Configuration * fsbc )
{
  return true;
} // end( ACNetworkBlock::is_feasible )

void ACNetworkData::serialize( netCDF::NcGroup & group ) const {

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

}  // end( ACNetworkData::serialize )

/*--------------------------------------------------------------------------*/

void ACNetworkBlock::serialize( netCDF::NcGroup & group ) const
{
 NetworkBlock::serialize( group );

 if( auto network_data = get_NetworkData() )
  // If a ACNetworkData is present, serialize it.
  network_data->serialize( group );

 auto NumberNodes = group.getDim( "NumberNodes" );

 if( ! v_ActiveDemand.empty() ) {
  // This ACNetworkBlock has active demand, so it is serialized.

  if( NumberNodes.isNull() )
   /* The dimension "NumberNodes" is not present in the group (which means
    * that a ACNetworkData is not present). However, the number of nodes can
    * still be obtained from the size of the active demand vector. Notice that
    * the name "NumberNodes" is not used for this new dimension, because it
    * would indicate that a ACNetworkData is present (which is not the
    * case). Therefore, we create an alternative dimension in order to be able
    * to serialize the active demand. */
   NumberNodes = group.addDim( "__NumberNodes__" , v_ActiveDemand.size() );

  // Finally, serialize the active demand.
  ::serialize( group , "ActiveDemand" , netCDF::NcDouble() ,
               NumberNodes , v_ActiveDemand );
 }

 if( f_ConstTerm != 0 )
  ::serialize( group , "ConstantTerm" , netCDF::NcDouble() , f_ConstTerm );

}  // end( ACNetworkBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

void ACNetworkBlock::set_active_demand( MF_dbl_it values ,
                                        Block::Subset && subset ,
                                        const bool ordered ,
                                        c_ModParam issuePMod ,
                                        c_ModParam issueAMod )
{

}  // end( ACNetworkData::set_active_demand )

/*--------------------------------------------------------------------------*/

void ACNetworkBlock::set_active_demand( MF_dbl_it values ,
                                        Block::Range rng ,
                                        c_ModParam issuePMod ,
                                        c_ModParam issueAMod )
{

}  // end( ACNetworkData::set_active_demand )

/*--------------------------------------------------------------------------*/

void ACNetworkBlock::set_kappa( MF_dbl_it values ,
                                Block::Subset && subset ,
                                const bool ordered ,
                                c_ModParam issuePMod ,
                                c_ModParam issueAMod )
{

}  // end( ACNetworkData::set_kappa( subset ) )

/*--------------------------------------------------------------------------*/

void ACNetworkBlock::set_kappa( MF_dbl_it values ,
                                Block::Range rng ,
                                c_ModParam issuePMod ,
                                c_ModParam issueAMod )
{

}  // end( ACNetworkData::set_kappa( range ) )

/*--------------------------------------------------------------------------*/
void ACNetworkBlock::change_power_flow_limit_constraints
(const std::vector<Index>& modified_lines, c_ModParam issueAMod){

}

/*--------------------------------------------------------------------------*/
void ACNetworkBlock::change_relax_abs_constraints
(const std::vector<Index>& modified_lines, c_ModParam issueAMod)
{

}

/*--------------------------------------------------------------------------*/
void ACNetworkBlock::change_DC_power_flow_injection_constraints
(const std::vector<Index>& modified_nodes, c_ModParam issueAMod)
{

}
