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

SMSpp_insert_in_factory_cpp_0( ECNetworkBlock );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

// register ECNetworkData to the NetworkData factory

typedef ECNetworkBlock::ECNetworkData ECNetworkData;

SMSpp_insert_in_factory_cpp_0( ECNetworkData );

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
 NetworkData::deserialize( group );

 if( f_number_nodes == 1 )
  throw( std::invalid_argument( "ECNetworkBlock::deserialize: cannot create "
                                "an Energy Community with just one user" ) );

 // Optional variables

 if( ! deserialize_dim( group , "NumberIntervals" , f_number_intervals ) )
  f_number_intervals = 1;

 // Mandatory variables

 ::deserialize( group , "BuyPrice" , f_number_intervals , v_BuyPrice ,
               false , true );

 ::deserialize( group , "SellPrice" , f_number_intervals , v_SellPrice ,
                false , true );

 ::deserialize( group , f_PeakTariff , "PeakTariff" , false );

 // Optional variables

 if( ! ::deserialize( group , "RewardPrice" , f_number_intervals ,
                      v_RewardPrice , true , true ) )
  v_RewardPrice.resize( f_number_intervals );

 ::deserialize( group , "PenaltyPrice" , f_number_intervals ,
                v_PenaltyPrice , true , true );

}  // end( ECNetworkData::deserialize )

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG

std::vector< std::string > ECNetworkData::expected_dims( void ) const {
 auto ret = NetworkData::expected_dims();
 ret.push_back( "NumberIntervals" );

 return( ret );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

std::vector< std::string > ECNetworkData::expected_vars( void ) const {
 static const std::vector< std::string > ev =
 { "BuyPrice" , "SellPrice" , "RewardPrice" , "PeakTariff" , "PenaltyPrice" };

 auto ret = NetworkData::expected_vars();
 ret.insert( ret.end() , ev.begin() , ev.end() );

 return( ret );
 }

#endif

/*--------------------------------------------------------------------------*/

void ECNetworkBlock::deserialize( const netCDF::NcGroup & group )
{
 // Optional variables
 Index NumberIntervals;
 if( ! deserialize_dim( group , "NumberIntervals" , NumberIntervals ) )
  NumberIntervals = 1;
 
 Index NumberNodes;
 if( deserialize_dim( group , "NumberNodes" , NumberNodes ) ) {
  // since the dimensions "NumberNodes" has been provided, an ECNetworkData
  // is there: deserialize it and mark it as local
  if( f_local_NetworkData ) {  // if the NetworkData has not been passed
   delete( f_NetworkData );    // from UCBlock, then delete it
   f_NetworkData = nullptr;
   }
  auto ECND = get_new_NetworkData();
  ECND->deserialize( group );
  if( f_NetworkData &&
    ( f_NetworkData->get_number_nodes() != ECND->get_number_nodes() ) )
   throw( std::logic_error( "ECNetworkBlock::deserialize: NumberNodes not "
			    "matching between NetworkData" ) );
  f_NetworkData = ECND;
  f_local_NetworkData = true;
  // an ECNetworkData has been provided, so the size of the given vector of
  // active demand must be equal to [ number of nodes x number of intervals ]
  ::deserialize( group , "ActiveDemand" ,
                 { NumberIntervals , NumberNodes } , v_ActiveDemand );
  }

 // note: NetworkBlock::deserialize() is called after dealing with the
 //       ECNetworkData, so that it's there when the list of expected stuff
 //       is constructed and checked
 NetworkBlock::deserialize( group );

 }  // end( ECNetworkBlock::deserialize )

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG

std::vector< std::string > ECNetworkBlock::expected_vars( void ) const {
 auto ret = NetworkBlock::expected_vars();
 ret.push_back( "ActiveDemand" );

 return( ret );
 }

#endif

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

 // the squilibrium variables are generated only if a PenaltyPrice has been
 // provided to the ECNetworkData; otherwise the model is identical to the
 // one without the imbalance term in the objective
 if( ! f_NetworkData->get_penalty_price().empty() ) {
  // the positive squilibrium variables
  v_power_squilibrium_pos.resize(
   boost::extents[ number_intervals ][ number_nodes ] );
  for( Index i = 0 ; i < number_intervals ; ++i )
   for( Index node_id = 0 ; node_id < number_nodes ; ++node_id )
    v_power_squilibrium_pos[ i ][ node_id ].set_type(
     ColVariable::kNonNegative );
  add_static_variable( v_power_squilibrium_pos , "p_sq_pos_network" );

  // the negative squilibrium variables
  v_power_squilibrium_neg.resize(
   boost::extents[ number_intervals ][ number_nodes ] );
  for( Index i = 0 ; i < number_intervals ; ++i )
   for( Index node_id = 0 ; node_id < number_nodes ; ++node_id )
    v_power_squilibrium_neg[ i ][ node_id ].set_type(
     ColVariable::kNonNegative );
  add_static_variable( v_power_squilibrium_neg , "p_sq_neg_network" );
 }

 set_variables_generated();

}  // end( ECNetworkBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void ECNetworkBlock::generate_abstract_constraints( Configuration * stcc )
{
 if( constraints_generated() )  // constraints have already been generated
  return;                       // nothing to do

 NetworkBlock::generate_abstract_constraints( stcc );

 const auto number_nodes = get_number_nodes();
 const auto number_intervals = get_number_intervals();

/*-------------------------- equality constraints --------------------------*/

 LinearFunction::v_coeff_pair vars;

 const bool has_imbalance = ! v_power_squilibrium_pos.empty();

 // set the power balance, i.e.:
 //
 //    P^+ - P^- - node_injection [ + P_sq^+ - P_sq^- ]
 //        = - active_demand                                   for all u, t
 //
 // the two squilibrium terms (in brackets) appear only when the PenaltyPrice
 // has been provided to the ECNetworkData

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

   if( has_imbalance ) {
    vars.push_back( std::make_pair( &v_power_squilibrium_pos[ i ][ node_id ] ,
                                    1.0 ) );
    vars.push_back( std::make_pair( &v_power_squilibrium_neg[ i ][ node_id ] ,
                                    -1.0 ) );
   }

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

 set_constraints_generated();

}  // end( ECNetworkBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void ECNetworkBlock::generate_objective( Configuration * objc )
{
 if( objective_generated() )  // Objective has already been generated
  return;                     // nothing to do

 const auto is_coop = is_cooperative();
 const bool has_imbalance = ! v_power_squilibrium_pos.empty();

 LinearFunction::v_coeff_pair vars;

 for( Index node_id = 0 ; node_id < get_number_nodes() ; ++node_id ) {

  for( Index t = 0 ; t < get_number_intervals() ; ++t ) {

   vars.push_back( std::make_pair( &v_power_absorption[ t ][ node_id ] ,
                                   get_buy_price( t ) ) );
   vars.push_back( std::make_pair( &v_power_injection[ t ][ node_id ] ,
                                   -get_sell_price( t ) ) );

   if( node_id == 0 && is_coop )
    vars.push_back( std::make_pair( &v_shared_power[ t ] ,
                                    -get_reward_price( t ) ) );

   // the squilibrium variables enter the objective with the same
   // `penalty_price[t]` coefficient on both the positive and the negative
   // one; the term is generated only when the squilibrium variables exist
   if( has_imbalance ) {
    const auto pp = get_penalty_price( t );
    vars.push_back( std::make_pair( &v_power_squilibrium_pos[ t ][ node_id ] ,
                                    pp ) );
    vars.push_back( std::make_pair( &v_power_squilibrium_neg[ t ][ node_id ] ,
                                    pp ) );
   }
  }

  vars.push_back( std::make_pair( &v_peak_power[ node_id ] ,
                                  get_peak_tariff() ) );
 }

 auto lf = new LinearFunction( std::move( vars ) );

 lf->set_constant_term( f_ConstTerm );

 objective.set_function( lf );
 objective.set_sense( Objective::eMin );

 // Set block objective
 this->set_objective( &objective , eNoMod );

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

 NetworkData::serialize( group );

 auto NumberIntervals = group.addDim( "NumberIntervals" , f_number_intervals );

 ::serialize( group , "BuyPrice" , netCDF::NcDouble() , NumberIntervals ,
              v_BuyPrice );

 ::serialize( group , "SellPrice" , netCDF::NcDouble() , NumberIntervals ,
              v_SellPrice );

 ::serialize( group , "PeakTariff" , netCDF::NcDouble() , f_PeakTariff );

 if( std::any_of( v_RewardPrice.begin() , v_RewardPrice.end() ,
                  []( double cst ) { return( cst != 0 ); } ) )
  ::serialize( group , "RewardPrice" , netCDF::NcDouble() , NumberIntervals ,
               v_RewardPrice );

 if( std::any_of( v_PenaltyPrice.begin() , v_PenaltyPrice.end() ,
                  []( double cst ) { return( cst != 0 ); } ) )
  ::serialize( group , "PenaltyPrice" , netCDF::NcDouble() , NumberIntervals ,
               v_PenaltyPrice );

}  // end( ECNetworkData::serialize )

/*--------------------------------------------------------------------------*/

void ECNetworkBlock::serialize( netCDF::NcGroup & group ) const
{
 NetworkBlock::serialize( group );

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
   boost::extents[ get_number_intervals() ][ get_number_nodes() ] );
 }

 bool identical = true;
 for( auto i : subset ) {
  if( i >= v_ActiveDemand.size() )
   throw( std::invalid_argument( "ECNetworkBlock::set_active_demand: "
                                 "invalid value in subset." ) );

  auto demand = *( values++ );
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
 rng.second = std::min( rng.second , get_number_intervals() * get_number_nodes() );
 if( rng.second <= rng.first )
  return;

 if( v_ActiveDemand.empty() ) {
  if( std::all_of( values ,
                   values + ( rng.second - rng.first ) ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  v_ActiveDemand.resize(
   boost::extents[ get_number_intervals() ][ get_number_nodes() ] );
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

void ECNetworkBlock::set_buy_price( MF_dbl_it values ,
                                    Block::Subset && subset ,
                                    const bool ordered ,
                                    c_ModParam issuePMod ,
                                    c_ModParam issueAMod )
{
 if( subset.empty() )
  return;

 auto & buy_price = f_NetworkData->get_buy_price();
 const auto T = get_number_intervals();

 if( buy_price.empty() ) {
  if( std::all_of( values , values + subset.size() ,
                   []( double cst ) { return( cst == 0.0 ); } ) )
   return;
  buy_price.assign( T , 0.0 );
 }

 bool identical = true;
 auto values_it = values;
 for( auto i : subset ) {
  if( i >= T )
   throw( std::invalid_argument(
    "ECNetworkBlock::set_buy_price: invalid value in subset: " +
    std::to_string( i ) ) );
  if( buy_price[ i ] != *( values_it++ ) ) {
   identical = false;
   break;
  }
 }

 if( identical )
  return;

 if( not_dry_run( issuePMod ) ) {
  values_it = values;
  for( auto i : subset )
   buy_price[ i ] = *( values_it++ );

  if( not_dry_run( issueAMod ) && objective_generated() ) {
   auto * lf = static_cast< LinearFunction * >( objective.get_function() );
   const auto N = get_number_nodes();
   for( auto i : subset )
    for( Index n = 0 ; n < N ; ++n ) {
     const auto idx = lf->is_active( &v_power_absorption[ i ][ n ] );
     if( idx == Inf< Index >() )
      throw( std::logic_error(
       "ECNetworkBlock::set_buy_price: expected Variable not found in "
       "objective." ) );
     lf->modify_coefficient( idx , buy_price[ i ] , issueAMod );
    }
  }
 }

 if( issue_pmod( issuePMod ) ) {
  if( ! ordered )
   std::sort( subset.begin() , subset.end() );
  Block::add_Modification( std::make_shared< ECNetworkBlockSbstMod >(
                            this , ECNetworkBlockMod::eSetBuyP ,
                            std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );
 }
}  // end( set_buy_price subset )

/*--------------------------------------------------------------------------*/

void ECNetworkBlock::set_buy_price( MF_dbl_it values ,
                                    Block::Range rng ,
                                    c_ModParam issuePMod ,
                                    c_ModParam issueAMod )
{
 const auto T = get_number_intervals();
 rng.second = std::min( rng.second , T );
 if( rng.second <= rng.first )
  return;

 c_Index sz = rng.second - rng.first;
 auto & buy_price = f_NetworkData->get_buy_price();

 if( buy_price.empty() ) {
  if( std::all_of( values , values + sz ,
                   []( double cst ) { return( cst == 0.0 ); } ) )
   return;
  buy_price.assign( T , 0.0 );
 }

 if( std::equal( values , values + sz , buy_price.begin() + rng.first ) )
  return;

 if( not_dry_run( issuePMod ) ) {
  std::copy( values , values + sz , buy_price.begin() + rng.first );

  if( not_dry_run( issueAMod ) && objective_generated() ) {
   auto * lf = static_cast< LinearFunction * >( objective.get_function() );
   const auto N = get_number_nodes();
   for( Index t = rng.first ; t < rng.second ; ++t )
    for( Index n = 0 ; n < N ; ++n ) {
     const auto idx = lf->is_active( &v_power_absorption[ t ][ n ] );
     if( idx == Inf< Index >() )
      throw( std::logic_error(
       "ECNetworkBlock::set_buy_price: expected Variable not found in "
       "objective." ) );
     lf->modify_coefficient( idx , buy_price[ t ] , issueAMod );
    }
  }
 }

 if( issue_pmod( issuePMod ) )
  Block::add_Modification( std::make_shared< ECNetworkBlockRngdMod >(
                            this , ECNetworkBlockMod::eSetBuyP , rng ) ,
                           Observer::par2chnl( issuePMod ) );
}  // end( set_buy_price range )

/*--------------------------------------------------------------------------*/

void ECNetworkBlock::set_sell_price( MF_dbl_it values ,
                                     Block::Subset && subset ,
                                     const bool ordered ,
                                     c_ModParam issuePMod ,
                                     c_ModParam issueAMod )
{
 if( subset.empty() )
  return;

 auto & sell_price = f_NetworkData->get_sell_price();
 const auto T = get_number_intervals();

 if( sell_price.empty() ) {
  if( std::all_of( values , values + subset.size() ,
                   []( double cst ) { return( cst == 0.0 ); } ) )
   return;
  sell_price.assign( T , 0.0 );
 }

 bool identical = true;
 auto values_it = values;
 for( auto i : subset ) {
  if( i >= T )
   throw( std::invalid_argument(
    "ECNetworkBlock::set_sell_price: invalid value in subset: " +
    std::to_string( i ) ) );
  if( sell_price[ i ] != *( values_it++ ) ) {
   identical = false;
   break;
  }
 }

 if( identical )
  return;

 if( not_dry_run( issuePMod ) ) {
  values_it = values;
  for( auto i : subset )
   sell_price[ i ] = *( values_it++ );

  if( not_dry_run( issueAMod ) && objective_generated() ) {
   auto * lf = static_cast< LinearFunction * >( objective.get_function() );
   const auto N = get_number_nodes();
   for( auto i : subset )
    for( Index n = 0 ; n < N ; ++n ) {
     const auto idx = lf->is_active( &v_power_injection[ i ][ n ] );
     if( idx == Inf< Index >() )
      throw( std::logic_error(
       "ECNetworkBlock::set_sell_price: expected Variable not found in "
       "objective." ) );
     lf->modify_coefficient( idx , -sell_price[ i ] , issueAMod );
    }
  }
 }

 if( issue_pmod( issuePMod ) ) {
  if( ! ordered )
   std::sort( subset.begin() , subset.end() );
  Block::add_Modification( std::make_shared< ECNetworkBlockSbstMod >(
                            this , ECNetworkBlockMod::eSetSellP ,
                            std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );
 }
}  // end( set_sell_price subset )

/*--------------------------------------------------------------------------*/

void ECNetworkBlock::set_sell_price( MF_dbl_it values ,
                                     Block::Range rng ,
                                     c_ModParam issuePMod ,
                                     c_ModParam issueAMod )
{
 const auto T = get_number_intervals();
 rng.second = std::min( rng.second , T );
 if( rng.second <= rng.first )
  return;

 c_Index sz = rng.second - rng.first;
 auto & sell_price = f_NetworkData->get_sell_price();

 if( sell_price.empty() ) {
  if( std::all_of( values , values + sz ,
                   []( double cst ) { return( cst == 0.0 ); } ) )
   return;
  sell_price.assign( T , 0.0 );
 }

 if( std::equal( values , values + sz , sell_price.begin() + rng.first ) )
  return;

 if( not_dry_run( issuePMod ) ) {
  std::copy( values , values + sz , sell_price.begin() + rng.first );

  if( not_dry_run( issueAMod ) && objective_generated() ) {
   auto * lf = static_cast< LinearFunction * >( objective.get_function() );
   const auto N = get_number_nodes();
   for( Index t = rng.first ; t < rng.second ; ++t )
    for( Index n = 0 ; n < N ; ++n ) {
     const auto idx = lf->is_active( &v_power_injection[ t ][ n ] );
     if( idx == Inf< Index >() )
      throw( std::logic_error(
       "ECNetworkBlock::set_sell_price: expected Variable not found in "
       "objective." ) );
     lf->modify_coefficient( idx , -sell_price[ t ] , issueAMod );
    }
  }
 }

 if( issue_pmod( issuePMod ) )
  Block::add_Modification( std::make_shared< ECNetworkBlockRngdMod >(
                            this , ECNetworkBlockMod::eSetSellP , rng ) ,
                           Observer::par2chnl( issuePMod ) );
}  // end( set_sell_price range )

/*--------------------------------------------------------------------------*/

void ECNetworkBlock::set_peak_tariff( MF_dbl_it values ,
                                      Block::Subset && subset ,
                                      const bool ordered ,
                                      c_ModParam issuePMod ,
                                      c_ModParam issueAMod )
{
 if( subset.empty() )
  return;
 if( subset.front() != 0 )
  throw( std::invalid_argument( "ECNetworkBlock::set_peak_tariff: PeakTariff "
                                "is a scalar, only index 0 is allowed." ) );

 auto & peak_tariff = f_NetworkData->get_peak_tariff();
 const auto v = *values;
 if( peak_tariff == v )
  return;

 if( not_dry_run( issuePMod ) ) {
  peak_tariff = v;

  if( not_dry_run( issueAMod ) && objective_generated() ) {
   auto * lf = static_cast< LinearFunction * >( objective.get_function() );
   const auto N = get_number_nodes();
   for( Index n = 0 ; n < N ; ++n ) {
    const auto idx = lf->is_active( &v_peak_power[ n ] );
    if( idx == Inf< Index >() )
     throw( std::logic_error(
      "ECNetworkBlock::set_peak_tariff: expected Variable not found in "
      "objective." ) );
    lf->modify_coefficient( idx , peak_tariff , issueAMod );
   }
  }
 }

 if( issue_pmod( issuePMod ) ) {
  if( ! ordered )
   std::sort( subset.begin() , subset.end() );
  Block::add_Modification( std::make_shared< ECNetworkBlockSbstMod >(
                            this , ECNetworkBlockMod::eSetPeakT ,
                            std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );
 }
}  // end( set_peak_tariff subset )

/*--------------------------------------------------------------------------*/

void ECNetworkBlock::set_peak_tariff( MF_dbl_it values ,
                                      Block::Range rng ,
                                      c_ModParam issuePMod ,
                                      c_ModParam issueAMod )
{
 if( rng.first > 0 || rng.second <= 0 )
  return;

 auto & peak_tariff = f_NetworkData->get_peak_tariff();
 const auto v = *values;
 if( peak_tariff == v )
  return;

 if( not_dry_run( issuePMod ) ) {
  peak_tariff = v;

  if( not_dry_run( issueAMod ) && objective_generated() ) {
   auto * lf = static_cast< LinearFunction * >( objective.get_function() );
   const auto N = get_number_nodes();
   for( Index n = 0 ; n < N ; ++n ) {
    const auto idx = lf->is_active( &v_peak_power[ n ] );
    if( idx == Inf< Index >() )
     throw( std::logic_error(
      "ECNetworkBlock::set_peak_tariff: expected Variable not found in "
      "objective." ) );
    lf->modify_coefficient( idx , peak_tariff , issueAMod );
   }
  }
 }

 if( issue_pmod( issuePMod ) )
  Block::add_Modification( std::make_shared< ECNetworkBlockRngdMod >(
                            this , ECNetworkBlockMod::eSetPeakT ,
                            Range( 0 , 1 ) ) ,
                           Observer::par2chnl( issuePMod ) );
}  // end( set_peak_tariff range )

/*--------------------------------------------------------------------------*/

void ECNetworkBlock::set_const_term( MF_dbl_it values ,
                                     Block::Subset && subset ,
                                     const bool ordered ,
                                     c_ModParam issuePMod ,
                                     c_ModParam issueAMod )
{
 if( subset.empty() )
  return;
 if( subset.front() != 0 )
  throw( std::invalid_argument( "ECNetworkBlock::set_const_term: ConstTerm "
                                "is a scalar, only index 0 is allowed." ) );

 const auto v = *values;
 if( get_const_term() == v )
  return;

 if( not_dry_run( issuePMod ) ) {
  set_constant_term( v );  // updates NetworkBlock::f_ConstTerm
  if( not_dry_run( issueAMod ) && objective_generated() ) {
   auto * lf = static_cast< LinearFunction * >( objective.get_function() );
   lf->set_constant_term( v , issueAMod );
  }
 }

 if( issue_pmod( issuePMod ) ) {
  if( ! ordered )
   std::sort( subset.begin() , subset.end() );
  Block::add_Modification( std::make_shared< ECNetworkBlockSbstMod >(
                            this , ECNetworkBlockMod::eSetConstT ,
                            std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );
 }
}  // end( set_const_term subset )

/*--------------------------------------------------------------------------*/

void ECNetworkBlock::set_const_term( MF_dbl_it values ,
                                     Block::Range rng ,
                                     c_ModParam issuePMod ,
                                     c_ModParam issueAMod )
{
 if( rng.first > 0 || rng.second <= 0 )
  return;

 const auto v = *values;
 if( get_const_term() == v )
  return;

 if( not_dry_run( issuePMod ) ) {
  set_constant_term( v );
  if( not_dry_run( issueAMod ) && objective_generated() ) {
   auto * lf = static_cast< LinearFunction * >( objective.get_function() );
   lf->set_constant_term( v , issueAMod );
  }
 }

 if( issue_pmod( issuePMod ) )
  Block::add_Modification( std::make_shared< ECNetworkBlockRngdMod >(
                            this , ECNetworkBlockMod::eSetConstT ,
                            Range( 0 , 1 ) ) ,
                           Observer::par2chnl( issuePMod ) );
}  // end( set_const_term range )

/*--------------------------------------------------------------------------*/

void ECNetworkBlock::set_penalty_price( MF_dbl_it values ,
                                        Block::Subset && subset ,
                                        const bool ordered ,
                                        c_ModParam issuePMod ,
                                        c_ModParam issueAMod )
{
 if( subset.empty() )
  return;

 auto & penalty_price = f_NetworkData->get_penalty_price();
 const auto T = get_number_intervals();

 if( penalty_price.empty() ) {
  if( std::all_of( values , values + subset.size() ,
                   []( double cst ) { return( cst == 0.0 ); } ) )
   return;
  penalty_price.assign( T , 0.0 );
 }

 bool changed = false;
 auto values_it = values;
 for( auto i : subset ) {
  if( i >= T )
   throw( std::invalid_argument(
    "ECNetworkBlock::set_penalty_price: invalid value in subset: " +
    std::to_string( i ) ) );
  const auto v = *( values_it++ );
  if( penalty_price[ i ] != v ) {
   changed = true;
   if( not_dry_run( issuePMod ) )
    penalty_price[ i ] = v;
  }
 }

 if( ! changed )
  return;

 // when the squilibrium variables are part of the model, the coefficients
 // of `v_power_squilibrium_pos/neg[i][n]` in the objective must be aligned
 // with the new `penalty_price[i]` for every node `n`
 if( ! v_power_squilibrium_pos.empty() && not_dry_run( issueAMod ) &&
     not_dry_run( issuePMod ) && objective_generated() ) {
  auto * lf = static_cast< LinearFunction * >( objective.get_function() );
  const auto N = get_number_nodes();
  for( auto i : subset )
   for( Index n = 0 ; n < N ; ++n ) {
    const auto idx_pos = lf->is_active( &v_power_squilibrium_pos[ i ][ n ] );
    const auto idx_neg = lf->is_active( &v_power_squilibrium_neg[ i ][ n ] );
    if( idx_pos == Inf< Index >() || idx_neg == Inf< Index >() )
     throw( std::logic_error(
      "ECNetworkBlock::set_penalty_price: expected Variable not found in "
      "objective." ) );
    lf->modify_coefficient( idx_pos , penalty_price[ i ] , issueAMod );
    lf->modify_coefficient( idx_neg , penalty_price[ i ] , issueAMod );
   }
 }

 if( issue_pmod( issuePMod ) ) {
  if( ! ordered )
   std::sort( subset.begin() , subset.end() );
  Block::add_Modification( std::make_shared< ECNetworkBlockSbstMod >(
                            this , ECNetworkBlockMod::eSetPenaltyP ,
                            std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );
 }
}  // end( set_penalty_price subset )

/*--------------------------------------------------------------------------*/

void ECNetworkBlock::set_penalty_price( MF_dbl_it values ,
                                        Block::Range rng ,
                                        c_ModParam issuePMod ,
                                        c_ModParam issueAMod )
{
 const auto T = get_number_intervals();
 rng.second = std::min( rng.second , T );
 if( rng.second <= rng.first )
  return;

 c_Index sz = rng.second - rng.first;
 auto & penalty_price = f_NetworkData->get_penalty_price();

 if( penalty_price.empty() ) {
  if( std::all_of( values , values + sz ,
                   []( double cst ) { return( cst == 0.0 ); } ) )
   return;
  penalty_price.assign( T , 0.0 );
 }

 if( std::equal( values , values + sz ,
                 penalty_price.begin() + rng.first ) )
  return;

 if( not_dry_run( issuePMod ) ) {
  std::copy( values , values + sz , penalty_price.begin() + rng.first );

  // when the squilibrium variables are part of the model, the coefficients
  // of `v_power_squilibrium_pos/neg[t][n]` in the objective must be aligned
  // with the new `penalty_price[t]` for every node `n`
  if( ! v_power_squilibrium_pos.empty() && not_dry_run( issueAMod ) &&
      objective_generated() ) {
   auto * lf = static_cast< LinearFunction * >( objective.get_function() );
   const auto N = get_number_nodes();
   for( Index t = rng.first ; t < rng.second ; ++t )
    for( Index n = 0 ; n < N ; ++n ) {
     const auto idx_pos = lf->is_active( &v_power_squilibrium_pos[ t ][ n ] );
     const auto idx_neg = lf->is_active( &v_power_squilibrium_neg[ t ][ n ] );
     if( idx_pos == Inf< Index >() || idx_neg == Inf< Index >() )
      throw( std::logic_error(
       "ECNetworkBlock::set_penalty_price: expected Variable not found in "
       "objective." ) );
     lf->modify_coefficient( idx_pos , penalty_price[ t ] , issueAMod );
     lf->modify_coefficient( idx_neg , penalty_price[ t ] , issueAMod );
    }
  }
 }

 if( issue_pmod( issuePMod ) )
  Block::add_Modification( std::make_shared< ECNetworkBlockRngdMod >(
                            this , ECNetworkBlockMod::eSetPenaltyP , rng ) ,
                           Observer::par2chnl( issuePMod ) );
}  // end( set_penalty_price range )

/*--------------------------------------------------------------------------*/
/*----------------------- End File ECNetworkBlock.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
