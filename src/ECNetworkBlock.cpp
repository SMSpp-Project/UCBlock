/*--------------------------------------------------------------------------*/
/*------------------------- File ECNetworkBlock.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the ECNetworkBlock class.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "NetworkBlock.h"

#include "ECNetworkBlock.h"

#include "LinearFunction.h"

#include "UCBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register ECNetworkBlock to the Block factory

SMSpp_insert_in_factory_cpp_1( ECNetworkBlock );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

// register ECNetworkData to the NetworkData factory

typedef ECNetworkBlock::ECNetworkData ECNetworkData;

SMSpp_insert_in_factory_cpp_1( ECNetworkData );

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS OF ECNetworkBlock ------------------------*/
/*--------------------------------------------------------------------------*/

ECNetworkBlock::~ECNetworkBlock()
{
 Constraint::clear( power_balance_const );
 Constraint::clear( power_shared_const );
 Constraint::clear( power_flow_limit_const );

 Constraint::clear( node_injection_bounds_const );

 objective.clear();

 // Delete the ECNetworkData if it is local.
 if( f_local_NetworkData )
  delete( f_NetworkData );
}

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void ECNetworkData::deserialize( const netCDF::NcGroup & group )
{

#ifndef NDEBUG
 static std::vector< std::string > expected_dims = { "NumberNodes" ,
                                                     "NumberIntervals" ,
                                                     // if called from UCBlock:
                                                     "TimeHorizon" ,
                                                     "NumberUnits" ,
                                                     "NumberNetworks" ,
                                                     "NumberElectricalGenerators" };
 check_dimensions( group , expected_dims , std::cerr );

 static std::vector< std::string > expected_vars = { "ActiveDemand" ,
                                                     "BuyPrice" ,
                                                     "SellPrice" ,
                                                     "RewardPrice" ,
                                                     "PeakTariff" ,
                                                     "ConstantTerm" ,
                                                     // if called from UCBlock:
                                                     "ActivePowerDemand" ,
                                                     "GeneratorNode" ,
                                                     "StartNetworkIntervals" ,
                                                     "NetworkConstantTerms" ,
                                                     "NetworkBlockClassname" ,
                                                     "NetworkDataClassname" };
 check_variables( group , expected_vars , std::cerr );
#endif

 // Mandatory variables

 ::deserialize_dim( group , "NumberNodes" , f_number_nodes , false );
 if( f_number_nodes == 1 )
  throw( std::invalid_argument( "ECNetworkBlock::deserialize: cannot create "
                                "an Energy Community with just one user" ) );

 // Optional variables

 ::deserialize( group , f_BuyPrice , "BuyPrice" );
 ::deserialize( group , f_SellPrice , "SellPrice" );
 ::deserialize( group , f_RewardPrice , "RewardPrice" );
 ::deserialize( group , f_PeakTariff , "PeakTariff" );

}  // end( ECNetworkData::deserialize )

/*--------------------------------------------------------------------------*/

void ECNetworkBlock::deserialize( const netCDF::NcGroup & group )
{

#ifndef NDEBUG
 static std::vector< std::string > expected_dims = { "NumberNodes" ,
                                                     "NumberIntervals" };
 check_dimensions( group , expected_dims , std::cerr );

 static std::vector< std::string > expected_vars = { "ActiveDemand" ,
                                                     "BuyPrice" ,
                                                     "SellPrice" ,
                                                     "RewardPrice" ,
                                                     "PeakTariff" ,
                                                     "ConstantTerm" };
 check_variables( group , expected_vars , std::cerr );
#endif

 // Optional variables

 Index NumberNodes;
 if( ::deserialize_dim( group , "NumberNodes" , NumberNodes ) &&
     ::deserialize_dim( group , "NumberIntervals" , f_number_intervals ) ) {
  // Since the dimensions "NumberNodes" and "NumberIntervals" has been provided,
  // it means that a ECNetworkData has been provided. Thus, the ECNetworkData
  // is deserialized, and it is marked as being local.
  delete( f_NetworkData );
  f_NetworkData = new ECNetworkData();
  f_NetworkData->deserialize( group );
  f_local_NetworkData = true;
  // An ECNetworkData has been provided. So, the size of the given vector of
  // active demand must be equal to the number of nodes.
  ::deserialize( group , "ActiveDemand" , v_ActiveDemand );
  // always check if the demand is given in the correct shape
  assert( ( v_ActiveDemand.shape()[ 0 ] == f_number_intervals ) &&
          ( v_ActiveDemand.shape()[ 1 ] == NumberNodes ) );
 }

 ::deserialize( group , "BuyPrice" , f_number_intervals , v_BuyPrice ,
                true , true );
 if( v_BuyPrice.size() == 1 )
  v_BuyPrice.resize( f_number_intervals , v_BuyPrice[ 0 ] );

 ::deserialize( group , "SellPrice" , f_number_intervals , v_SellPrice ,
                true , true );
 if( v_SellPrice.size() == 1 )
  v_SellPrice.resize( f_number_intervals , v_SellPrice[ 0 ] );

 if( ! ::deserialize( group , "RewardPrice" , f_number_intervals ,
                     v_RewardPrice , true , true ) )
  v_RewardPrice.resize( f_number_intervals , 0 );
 else if( v_RewardPrice.size() == 1 )
  v_RewardPrice.resize( f_number_intervals , v_RewardPrice[ 0 ] );

 ::deserialize( group , f_PeakTariff , "PeakTariff" );

 ::deserialize( group , f_ConstTerm , "ConstantTerm" );

}  // end( ECNetworkBlock::deserialize )

/*--------------------------------------------------------------------------*/

void ECNetworkBlock::generate_abstract_variables( Configuration * stvv )
{
 if( variables_generated() )  // variables have already been generated
  return;                     // nothing to do

 NetworkBlock::generate_abstract_variables( stvv );

 const auto number_nodes = get_number_nodes();
 const auto number_intervals = get_number_intervals();

 // the public power injection variables
 v_power_injection.resize(
  boost::extents[ number_intervals ][ number_nodes ] );
 for( Index i = 0 ; i < number_intervals ; ++i )
  for( Index node_id = 0 ; node_id < number_nodes ; ++node_id )
   v_power_injection[ i ][ node_id ].set_type( ColVariable::kNonNegative );
 add_static_variable( v_power_injection , "p_inj_network" );

 // the public power absorption variables
 v_power_absorption.resize(
  boost::extents[ number_intervals ][ number_nodes ] );
 for( Index i = 0 ; i < number_intervals ; ++i )
  for( Index node_id = 0 ; node_id < number_nodes ; ++node_id )
   v_power_absorption[ i ][ node_id ].set_type( ColVariable::kNonNegative );
 add_static_variable( v_power_absorption , "p_abs_network" );

 if( is_cooperative() ) {
  // the microgrid power variables
  v_shared_power.resize( number_intervals );
  for( auto & var : v_shared_power )
   var.set_type( ColVariable::kNonNegative );
  add_static_variable( v_shared_power , "p_shared_network" );
 }

 // the peak power variables
 v_peak_power.resize( number_nodes );
 for( auto & var : v_peak_power )
  var.set_type( ColVariable::kNonNegative );
 add_static_variable( v_peak_power , "p_peak_network" );

 set_variables_generated();

}  // end( ECNetworkBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void ECNetworkBlock::generate_abstract_constraints( Configuration * stcc )
{
 if( constraints_generated() )  // constraints have already been generated
  return;                       // nothing to do

 const auto number_nodes = get_number_nodes();
 const auto number_intervals = get_number_intervals();

/*-------------------------- equality constraints --------------------------*/

 LinearFunction::v_coeff_pair vars;

 // set the power balance, i.e.:
 //
 //    P^+ - P^- - node_injection = - active_demand   for all u, t

 power_balance_const.resize(
  boost::multi_array< FRowConstraint , 2 >::extent_gen()
  [ number_nodes ][ number_intervals ] );

 for( Index i = 0 ; i < number_intervals ; ++i )

  for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {

   vars.push_back( std::make_pair( &v_power_injection[ i ][ node_id ] ,
                                   1.0 ) );
   vars.push_back( std::make_pair( &v_power_absorption[ i ][ node_id ] ,
                                   -1.0 ) );
   vars.push_back( std::make_pair( &v_node_injection[ i ][ node_id ] , -1.0 ) );

   power_balance_const[ node_id ][ i ].set_both(
    -v_ActiveDemand[ i ][ node_id ] );
   power_balance_const[ node_id ][ i ].set_function(
    new LinearFunction( std::move( vars ) ) );
  }

 add_static_constraint( power_balance_const , "Power_Balance_Const_Network" );

/*------------------------- inequality constraints -------------------------*/

 LinearFunction::v_coeff_pair vars_inj;
 LinearFunction::v_coeff_pair vars_abs;

 // max shared power constraints within the microgrid market, i.e.:
 //
 //    P^M <= P^+       for all u, t     (1)
 // => P^M - P^+ <= 0   for all u, t
 //
 //    P^M <= P^-       for all u, t     (2)
 // => P^M - P^- <= 0   for all u, t

 if( is_cooperative() ) {

  power_shared_const.resize(
   boost::multi_array< FRowConstraint , 2 >::extent_gen()
   [ number_intervals ][ 2 ] ); // 2 dims, i.e., injection (+) and absorption (-)

  for( Index i = 0 ; i < number_intervals ; ++i ) {

   for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {

    // case (1)
    vars_inj.push_back( std::make_pair( &v_power_injection[ i ][ node_id ] ,
                                        -1.0 ) );

    // case (2)
    vars_abs.push_back( std::make_pair( &v_power_absorption[ i ][ node_id ] ,
                                        -1.0 ) );
   }

   // case (1)
   vars_inj.push_back( std::make_pair( &v_shared_power[ i ] , 1.0 ) );

   power_shared_const[ i ][ 0 ].set_lhs( -Inf< double >() );
   power_shared_const[ i ][ 0 ].set_rhs( 0.0 );
   power_shared_const[ i ][ 0 ].set_function(
    new LinearFunction( std::move( vars_inj ) ) );

   // case (2)
   vars_abs.push_back( std::make_pair( &v_shared_power[ i ] , 1.0 ) );

   power_shared_const[ i ][ 1 ].set_lhs( -Inf< double >() );
   power_shared_const[ i ][ 1 ].set_rhs( 0.0 );
   power_shared_const[ i ][ 1 ].set_function(
    new LinearFunction( std::move( vars_abs ) ) );
  }

  add_static_constraint( power_shared_const ,
                         "Power_Shared_Const_Network" );
 }

 // set that the dispatch cannot go beyond the maximum dispatch of the
 // corresponding peak power period, i.e.:
 //
 //    P^{max} >= P^+ - P^-         for all u, t     (1)
 // => P^+ - P^- - P^{max} <= 0     for all u, t
 //
 //    P^{max} >= - [ P^+ - P^- ]   for all u, t     (2)
 // => - P^+ + P^- - P^{max} <= 0   for all u, t

 power_flow_limit_const.resize(
  boost::multi_array< FRowConstraint , 3 >::extent_gen()
  [ number_nodes ][ number_intervals ][ 2 ] );  // 2 dims, i.e., the sign (+/-)

 for( Index i = 0 ; i < number_intervals ; ++i )

  for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {

   // case (1)
   vars_inj.push_back( std::make_pair( &v_power_injection[ i ][ node_id ] ,
                                       1.0 ) );
   vars_inj.push_back( std::make_pair( &v_power_absorption[ i ][ node_id ] ,
                                       -1.0 ) );
   vars_inj.push_back( std::make_pair( &v_peak_power[ node_id ] , 1.0 ) );

   // case (2)
   vars_abs.push_back( std::make_pair( &v_power_injection[ i ][ node_id ] ,
                                       -1.0 ) );
   vars_abs.push_back( std::make_pair( &v_power_absorption[ i ][ node_id ] ,
                                       1.0 ) );
   vars_abs.push_back( std::make_pair( &v_peak_power[ node_id ] , 1.0 ) );

   // case (1)
   power_flow_limit_const[ node_id ][ i ][ 0 ].set_lhs( 0.0 );
   power_flow_limit_const[ node_id ][ i ][ 0 ].set_rhs( Inf< double >() );
   power_flow_limit_const[ node_id ][ i ][ 0 ].set_function(
    new LinearFunction( std::move( vars_inj ) ) );

   // case (2)
   power_flow_limit_const[ node_id ][ i ][ 1 ].set_lhs( 0.0 );
   power_flow_limit_const[ node_id ][ i ][ 1 ].set_rhs( Inf< double >() );
   power_flow_limit_const[ node_id ][ i ][ 1 ].set_function(
    new LinearFunction( std::move( vars_abs ) ) );
  }

 add_static_constraint( power_flow_limit_const ,
                        "Power_Flow_Limit_Const_Network" );

 // node injection bound constraints

 node_injection_bounds_const.resize(
  boost::multi_array< FRowConstraint , 2 >::extent_gen()
  [ number_nodes ][ number_intervals ] );

 for( Index i = 0 ; i < number_intervals ; ++i )

  for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {

   node_injection_bounds_const[ node_id ][ i ].set_lhs(
    v_MinNodeInjection[ i ][ node_id ] );
   node_injection_bounds_const[ node_id ][ i ].set_rhs(
    v_MaxNodeInjection[ i ][ node_id ] );
   node_injection_bounds_const[ node_id ][ i ].set_variable(
    &v_node_injection[ i ][ node_id ] );
  }

 add_static_constraint( node_injection_bounds_const ,
                        "Node_Injection_Bound_Const_Network" );

 set_constraints_generated();

}  // end( ECNetworkBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void ECNetworkBlock::generate_objective( Configuration * objc )
{
 if( objective_generated() )  // Objective has already been generated
  return;                     // nothing to do

 const auto is_coop = is_cooperative();

 LinearFunction::v_coeff_pair vars;

 for( Index node_id = 0 ; node_id < get_number_nodes() ; ++node_id ) {

  for( Index t = 0 ; t < get_number_intervals() ; ++t ) {

   vars.push_back( std::make_pair( &v_power_absorption[ t ][ node_id ] ,
                                   get_buy_price( t ) ) );
   vars.push_back( std::make_pair( &v_power_injection[ t ][ node_id ] ,
                                   -get_sell_price( t ) ) );

   if( is_coop )
    vars.push_back( std::make_pair( &v_shared_power[ t ] ,
                                    -get_reward_price( t ) ) );
  }

  vars.push_back( std::make_pair( &v_peak_power[ node_id ] ,
                                  get_peak_tariff() ) );
 }

 auto lf = new LinearFunction( std::move( vars ) );

 lf->set_constant_term( f_ConstTerm );

 objective.set_function( lf );
 objective.set_sense( Objective::eMin );

 // Set block objective
 this->set_objective( &objective );

 set_objective_generated();

}  // end( ECNetworkBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*----------------- METHODS FOR CHECKING THE ECNetworkBlock ----------------*/
/*--------------------------------------------------------------------------*/

bool ECNetworkBlock::is_feasible( bool useabstract , Configuration * fsbc )
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
  && ColVariable::is_feasible( v_node_injection )
  && ColVariable::is_feasible( v_power_injection )
  && ColVariable::is_feasible( v_power_absorption )
  && ColVariable::is_feasible( v_shared_power )
  && ColVariable::is_feasible( v_peak_power )
  // Constraints
  && RowConstraint::is_feasible( power_balance_const , tol , rel_viol )
  && RowConstraint::is_feasible( power_shared_const , tol , rel_viol )
  && RowConstraint::is_feasible( power_flow_limit_const , tol , rel_viol )
  && RowConstraint::is_feasible( node_injection_bounds_const , tol , rel_viol ) );

}  // end( ECNetworkBlock::is_feasible )

/*--------------------------------------------------------------------------*/
/*--------- METHODS FOR LOADING, PRINTING & SAVING THE ECNetworkBlock ------*/
/*--------------------------------------------------------------------------*/

void ECNetworkData::serialize( netCDF::NcGroup & group ) const
{
 group.addDim( "NumberNodes" , f_number_nodes );

 ::serialize( group , "BuyPrice" , netCDF::NcDouble() , f_BuyPrice );
 ::serialize( group , "SellPrice" , netCDF::NcDouble() , f_SellPrice );
 ::serialize( group , "PeakTariff" , netCDF::NcDouble() , f_PeakTariff );

 if( f_RewardPrice != 0 )
  ::serialize( group , "RewardPrice" , netCDF::NcDouble() , f_RewardPrice );

}  // end( ECNetworkData::serialize )

/*--------------------------------------------------------------------------*/

void ECNetworkBlock::serialize( netCDF::NcGroup & group ) const
{
 NetworkBlock::serialize( group );

 auto NumberIntervals = group.getDim( "NumberIntervals" );

 ::serialize( group , "BuyPrice" , netCDF::NcDouble() , NumberIntervals ,
              v_BuyPrice );

 ::serialize( group , "SellPrice" , netCDF::NcDouble() , NumberIntervals ,
              v_SellPrice );

 ::serialize( group , "PeakTariff" , netCDF::NcDouble() , f_PeakTariff );

 if( auto network_data = get_NetworkData() )
  // If an ECNetworkData is present, serialize it.
  network_data->serialize( group );

 auto NumberNodes = group.getDim( "NumberNodes" );

 if( ! v_ActiveDemand.empty() ) {
  // This ECNetworkBlock has active demand, so it is serialized.

  if( NumberNodes.isNull() )
   /* The dimension "NumberNodes" is not present in the group (which means
    * that an ECNetworkData is not present). However, the number of nodes can
    * still be obtained from the size of the active demand vector. Notice that
    * the name "NumberNodes" is not used for this new dimension, because it
    * would indicate that an ECNetworkData is present (which is not the
    * case). Therefore, we create an alternative dimension in order to be able
    * to serialize the active demand. */
   NumberNodes = group.addDim( "__NumberNodes__" , v_ActiveDemand.size() );

  auto NumberIntervals = group.getDim( "NumberIntervals" );

  ::serialize( group , "ActiveDemand" , netCDF::NcDouble() ,
               { NumberIntervals , NumberNodes } , v_ActiveDemand );
 }

 if( std::any_of( v_RewardPrice.begin() , v_RewardPrice.end() ,
                  []( double cst ) { return( cst != 0 ); } ) )
  ::serialize( group , "RewardPrice" , netCDF::NcDouble() , NumberIntervals ,
               v_RewardPrice );

 if( f_ConstTerm != 0 )
  ::serialize( group , "ConstantTerm" , netCDF::NcDouble() , f_ConstTerm );

}  // end( ECNetworkBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

void ECNetworkBlock::set_active_demand( MF_dbl_it values ,
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

  v_ActiveDemand.resize(
   boost::extents[ f_number_intervals ][ get_number_nodes() ] );
 }

 bool identical = true;
 for( auto i : subset ) {
  if( i >= v_ActiveDemand.size() )
   throw( std::invalid_argument( "ECNetworkBlock::set_active_demand: "
                                 "invalid value in subset." ) );

  auto demand = *(values++);
  if( *( v_ActiveDemand.data() + i ) != demand ) {
   identical = false;

   if( not_dry_run( issuePMod ) )
    // Change the physical representation
    *( v_ActiveDemand.data() + i ) = demand;
  }
 }

 if( identical )
  return;  // nothing changes; return

 if( not_dry_run( issuePMod ) &&
     not_dry_run( issueAMod ) &&
     constraints_generated() ) {
  // Change the abstract representation

  for( auto i : subset ) {
   Index t = i % get_number_nodes();
   Index n = i / get_number_nodes();

   power_balance_const[ n ][ t ].set_both( -v_ActiveDemand[ t ][ n ] ,
                                           issueAMod );
  }
 }

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification
  if( ! ordered )
   std::sort( subset.begin() , subset.end() );

  Block::add_Modification( std::make_shared< ECNetworkBlockSbstMod >(
                            this , ECNetworkBlockMod::eSetActD ,
                            std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );
 }
}  // end( ECNetworkBlock::set_active_demand( subset ) )

/*--------------------------------------------------------------------------*/

void ECNetworkBlock::set_active_demand( MF_dbl_it values ,
                                        Block::Range rng ,
                                        c_ModParam issuePMod ,
                                        c_ModParam issueAMod )
{
 rng.second = std::min( rng.second , f_number_intervals * get_number_nodes() );
 if( rng.second <= rng.first )
  return;

 if( v_ActiveDemand.empty() ) {
  if( std::all_of( values ,
                   values + ( rng.second - rng.first ) ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  v_ActiveDemand.resize(
   boost::extents[ f_number_intervals ][ get_number_nodes() ] );
 }

 // If nothing changes, return
 if( std::equal( values ,
                 values + ( rng.second - rng.first ) ,
                 v_ActiveDemand.data() + rng.first ) )
  return;

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation
  std::copy( values ,
             values + ( rng.second - rng.first ) ,
             v_ActiveDemand.data() + rng.first );

  if( not_dry_run( issueAMod ) && constraints_generated() ) {
   // Change the abstract representation

   for( Index i = rng.first ; i < rng.second ; ++i ) {
    Index t = i % get_number_nodes();
    Index n = i / get_number_nodes();

    power_balance_const[ n ][ t ].set_both( -v_ActiveDemand[ t ][ n ] ,
                                            issueAMod );
   }
  }
 }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ECNetworkBlockRngdMod >(
                            this , ECNetworkBlockMod::eSetActD , rng ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ECNetworkBlock::set_active_demand( range ) )

/*--------------------------------------------------------------------------*/
/*----------------------- End File ECNetworkBlock.cpp ----------------------*/
/*--------------------------------------------------------------------------*/