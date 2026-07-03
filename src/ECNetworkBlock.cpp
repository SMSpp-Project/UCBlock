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

// register ECNetworkBlockSolution to the Solution factory

SMSpp_insert_in_factory_cpp_0( ECNetworkBlockSolution );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

// register ECNetworkData to the NetworkData factory

typedef ECNetworkBlock::ECNetworkData ECNetworkData;

SMSpp_insert_in_factory_cpp_0( ECNetworkData );

/*--------------------------------------------------------------------------*/
/*--------------------------- STATIC FUNCTIONS -----------------------------*/
/*--------------------------------------------------------------------------*/

template< class T , std::size_t K >
static void copy_multi_array( boost::multi_array< T , K > & to ,
                              const boost::multi_array< T , K > & from )
{
 std::vector< size_t > extent;
 auto shape = from.shape();
 extent.assign( shape , shape + from.num_dimensions() );
 to.resize( extent );
 to = from;
 }

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS OF ECNetworkBlock ------------------------*/
/*--------------------------------------------------------------------------*/

ECNetworkBlock::~ECNetworkBlock()
{
 Constraint::clear( power_balance_const );
 Constraint::clear( power_balance_agg_const );
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

 // the aggregate squilibrium variables and the aggregate declared-dispatch
 // (day-ahead bid) variables: only generated if a PenaltyPrice has been
 // provided to the ECNetworkData. The actual aggregate export then deviates
 // from the declared-dispatch by the squilibrium, penalized in the Objective.
 if( ! f_NetworkData->get_penalty_price().empty() ) {
  v_power_squilibrium_pos.resize( number_intervals );
  for( auto & var : v_power_squilibrium_pos )
   var.set_type( ColVariable::kNonNegative );
  add_static_variable( v_power_squilibrium_pos , "p_sq_pos_network" );

  v_power_squilibrium_neg.resize( number_intervals );
  for( auto & var : v_power_squilibrium_neg )
   var.set_type( ColVariable::kNonNegative );
  add_static_variable( v_power_squilibrium_neg , "p_sq_neg_network" );

  v_power_agg_dec_pos.resize( number_intervals );
  for( auto & var : v_power_agg_dec_pos )
   var.set_type( ColVariable::kNonNegative );
  add_static_variable( v_power_agg_dec_pos , "p_agg_dec_pos_network" );

  v_power_agg_dec_neg.resize( number_intervals );
  for( auto & var : v_power_agg_dec_neg )
   var.set_type( ColVariable::kNonNegative );
  add_static_variable( v_power_agg_dec_neg , "p_agg_dec_neg_network" );
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

 // set the per-node power balance, i.e.:
 //
 //    P^+ - P^- - node_injection = - active_demand            for all u, t

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

 // tie the actual aggregate public-market exchange to the declared-dispatch
 // (day-ahead bid) plus the squilibrium, i.e.:
 //
 //    Σ_n ( P^+_n - P^-_n ) - ( P_agg_dec^+ - P_agg_dec^- )
 //        - ( P_sq^+ - P_sq^- ) = 0                                for all t
 //
 // i.e., the actual aggregate net export equals the declared one plus the
 // imbalance; the squilibrium is penalized in the Objective. With a single
 // short-period scenario the declared-dispatch is free to equal the actual
 // export, so the imbalance is zero. The constraint is generated only when
 // the squilibrium / declared-dispatch variables exist.
 if( has_imbalance ) {
  power_balance_agg_const.resize( number_intervals );

  for( Index i = 0 ; i < number_intervals ; ++i ) {
   for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {
    vars.push_back( std::make_pair( &v_power_injection[ i ][ node_id ] ,
                                    1.0 ) );
    vars.push_back( std::make_pair( &v_power_absorption[ i ][ node_id ] ,
                                    -1.0 ) );
   }

   vars.push_back( std::make_pair( &v_power_agg_dec_pos[ i ] , -1.0 ) );
   vars.push_back( std::make_pair( &v_power_agg_dec_neg[ i ] , 1.0 ) );
   vars.push_back( std::make_pair( &v_power_squilibrium_pos[ i ] , -1.0 ) );
   vars.push_back( std::make_pair( &v_power_squilibrium_neg[ i ] , 1.0 ) );

   power_balance_agg_const[ i ].set_both( 0.0 );
   power_balance_agg_const[ i ].set_function(
    new LinearFunction( std::move( vars ) ) );
  }

  add_static_constraint( power_balance_agg_const ,
                         "Power_Balance_Agg_Const_Network" );
 }

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
  }

  vars.push_back( std::make_pair( &v_peak_power[ node_id ] ,
                                  get_peak_tariff() ) );
 }

 // the squilibrium variables enter the objective with the same
 // `penalty_price[t]` coefficient on both the positive and the negative
 // one; the term is generated only when the squilibrium variables exist
 if( has_imbalance )
  for( Index t = 0 ; t < get_number_intervals() ; ++t ) {
   const auto pp = get_penalty_price( t );
   vars.push_back( std::make_pair( &v_power_squilibrium_pos[ t ] , pp ) );
   vars.push_back( std::make_pair( &v_power_squilibrium_neg[ t ] , pp ) );
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
  && ColVariable::is_feasible( v_power_squilibrium_pos )
  && ColVariable::is_feasible( v_power_squilibrium_neg )
  && ColVariable::is_feasible( v_power_agg_dec_pos )
  && ColVariable::is_feasible( v_power_agg_dec_neg )
  // Constraints
  && RowConstraint::is_feasible( power_balance_const , tol , rel_viol )
  && RowConstraint::is_feasible( power_balance_agg_const , tol , rel_viol )
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
   NumberNodes = group.addDim( "__NumberNodes__" ,
                               v_ActiveDemand.shape()[ 1 ] );

  auto NumberIntervals = group.getDim( "NumberIntervals" );
  if( NumberIntervals.isNull() )
   // no ECNetworkData in the group: like for "NumberNodes" above, an
   // alternative dimension is created to serialize the active demand
   NumberIntervals = group.addDim( "__NumberIntervals__" ,
                                   v_ActiveDemand.shape()[ 0 ] );

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

 // the squilibrium variables exist only when "PenaltyPrice" is provided,
 // so guard the coefficient update accordingly
 if( ! v_power_squilibrium_pos.empty() && not_dry_run( issueAMod ) &&
     not_dry_run( issuePMod ) && objective_generated() ) {
  auto * lf = static_cast< LinearFunction * >( objective.get_function() );
  for( auto i : subset ) {
   const auto idx_pos = lf->is_active( &v_power_squilibrium_pos[ i ] );
   const auto idx_neg = lf->is_active( &v_power_squilibrium_neg[ i ] );
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

  // the squilibrium variables exist only when "PenaltyPrice" is provided,
  // so guard the coefficient update accordingly
  if( ! v_power_squilibrium_pos.empty() && not_dry_run( issueAMod ) &&
      objective_generated() ) {
   auto * lf = static_cast< LinearFunction * >( objective.get_function() );
   for( Index t = rng.first ; t < rng.second ; ++t ) {
    const auto idx_pos = lf->is_active( &v_power_squilibrium_pos[ t ] );
    const auto idx_neg = lf->is_active( &v_power_squilibrium_neg[ t ] );
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

Solution * ECNetworkBlock::get_Solution( Configuration * csolc , bool emptys )
{
 Index wsol = 31;  // by default: save everything
 if( ( ! csolc ) && f_BlockConfig )
  csolc = f_BlockConfig->f_solution_Configuration;

 if( auto config = dynamic_cast< SimpleConfiguration< int > * >( csolc ) )
  wsol = config->f_value;

 // call the method of the base class with empty=true; the base will create
 // the right Solution type via new_Solution(), and size v_node_injection if
 // wsol & 1. We size the EC-specific fields ourselves below, then -- if
 // requested -- call read() at the end
 auto * sol = dynamic_cast< ECNetworkBlockSolution * >(
                          NetworkBlock::get_Solution( csolc , true ) );
 assert( sol );

 const auto ni = get_number_intervals();
 const auto nn = get_number_nodes();

 if( wsol & 2 ) {
  // public-market injection / absorption variables
  sol->v_power_injection.resize(
   boost::multi_array< double , 2 >::extent_gen()[ ni ][ nn ] );
  sol->v_power_absorption.resize(
   boost::multi_array< double , 2 >::extent_gen()[ ni ][ nn ] );
  }

 if( wsol & 4 )
  sol->v_shared_power.resize( ni );

 if( wsol & 8 )
  sol->v_peak_power.resize( nn );

 if( wsol & 16 ) {
  // squilibrium and declared-dispatch variables: only present when the
  // underlying ECNetworkBlock generates them (i.e. when a PenaltyPrice has
  // been provided, so that the vectors are non-empty)
  if( ! v_power_squilibrium_pos.empty() )
   sol->v_power_squilibrium_pos.resize( ni );
  if( ! v_power_squilibrium_neg.empty() )
   sol->v_power_squilibrium_neg.resize( ni );
  if( ! v_power_agg_dec_pos.empty() )
   sol->v_power_agg_dec_pos.resize( ni );
  if( ! v_power_agg_dec_neg.empty() )
   sol->v_power_agg_dec_neg.resize( ni );
  }

 if( ! emptys )
  sol->read( this );

 return( sol );

 }  // end( ECNetworkBlock::get_Solution )

/*--------------------------------------------------------------------------*/

NetworkBlockSolution * ECNetworkBlock::new_Solution( void ) const {
 return( new ECNetworkBlockSolution() );
 }

/*--------------------------------------------------------------------------*/
/*------------------- METHODS OF ECNetworkBlockSolution --------------------*/
/*--------------------------------------------------------------------------*/

void ECNetworkBlockSolution::deserialize( const netCDF::NcGroup & group )
{
 // call the method of the base class
 NetworkBlockSolution::deserialize( group );

 const auto ni = f_number_intervals;
 const auto nn = f_number_nodes;

 using mad2i = MAdouble::index;

 // deserialize PowerInjection - - - - - - - - - - - - - - - - - - - - - - -
 if( ! ::deserialize< double , 2 >( group , "PowerInjection" , { ni , nn } ,
                                    v_power_injection , true ) ) {
  std::vector< mad2i > sizes( 2 , 0 );
  v_power_injection.resize( sizes );
  }

 // deserialize PowerAbsorption- - - - - - - - - - - - - - - - - - - - - - -
 if( ! ::deserialize< double , 2 >( group , "PowerAbsorption" , { ni , nn } ,
                                    v_power_absorption , true ) ) {
  std::vector< mad2i > sizes( 2 , 0 );
  v_power_absorption.resize( sizes );
  }

 // deserialize SharedPower- - - - - - - - - - - - - - - - - - - - - - - - -
 ::deserialize< double >( group , "SharedPower" , v_shared_power , true );

 // deserialize PeakPower- - - - - - - - - - - - - - - - - - - - - - - - - -
 ::deserialize< double >( group , "PeakPower" , v_peak_power , true );

 // deserialize PowerSquilibriumPos / Neg- - - - - - - - - - - - - - - - - -
 ::deserialize< double >( group , "PowerSquilibriumPos" ,
                          v_power_squilibrium_pos , true );

 ::deserialize< double >( group , "PowerSquilibriumNeg" ,
                          v_power_squilibrium_neg , true );

 // deserialize PowerAggDecPos / Neg - - - - - - - - - - - - - - - - - - - - -
 ::deserialize< double >( group , "PowerAggDecPos" ,
                          v_power_agg_dec_pos , true );

 ::deserialize< double >( group , "PowerAggDecNeg" ,
                          v_power_agg_dec_neg , true );

 }  // end( ECNetworkBlockSolution::deserialize( NcGroup & ) )

/*--------------------------------------------------------------------------*/

void ECNetworkBlockSolution::deserialize( const netCDF::NcGroup & group ,
                                          size_t idx )
{
 // call the method of the base class: this fills f_number_nodes,
 // f_number_intervals, and v_node_injection (if present)
 NetworkBlockSolution::deserialize( group , idx );

 const auto ni = f_number_intervals;
 const auto nn = f_number_nodes;

 using mad2i = MAdouble::index;

 // recover the offset in the time-indexed variables for this :Solution: it
 // is the same one used by the base class. If "EndInstant" is defined, the
 // start is EndInstant[ idx - 1 ] (or 0 for idx == 0). Otherwise, each
 // :NetworkBlockSolution covers exactly one instant and the start is idx
 size_t start;
 auto EI = group.getVar( "EndInstant" );
 if( EI.isNull() )
  start = idx;
 else {
  if( idx == 0 )
   start = 0;
  else {
   std::vector< size_t > vidx = { idx - 1 };
   int eitmp = 0;
   EI.getVar( vidx , & eitmp );
   start = eitmp;
   }
  }

 // deserialize PowerInjection - - - - - - - - - - - - - - - - - - - - - - -
 auto ncVar = group.getVar( "PowerInjection" );
 if( ncVar.isNull() ) {
  std::vector< mad2i > sizes = { 0 , 0 };
  v_power_injection.resize( sizes );
  }
 else {
  std::vector< mad2i > sizes = { ni , nn };
  v_power_injection.resize( sizes );
  std::vector< size_t > strt = { start , 0 };
  std::vector< size_t > cnt = { ni , nn };
  ncVar.getVar( strt , cnt , v_power_injection.data() );
  }

 // deserialize PowerAbsorption- - - - - - - - - - - - - - - - - - - - - - -
 ncVar = group.getVar( "PowerAbsorption" );
 if( ncVar.isNull() ) {
  std::vector< mad2i > sizes = { 0 , 0 };
  v_power_absorption.resize( sizes );
  }
 else {
  std::vector< mad2i > sizes = { ni , nn };
  v_power_absorption.resize( sizes );
  std::vector< size_t > strt = { start , 0 };
  std::vector< size_t > cnt = { ni , nn };
  ncVar.getVar( strt , cnt , v_power_absorption.data() );
  }

 // deserialize SharedPower- - - - - - - - - - - - - - - - - - - - - - - - -
 ncVar = group.getVar( "SharedPower" );
 if( ncVar.isNull() )
  v_shared_power.clear();
 else {
  v_shared_power.resize( ni );
  std::vector< size_t > strt = { start };
  std::vector< size_t > cnt = { ni };
  ncVar.getVar( strt , cnt , v_shared_power.data() );
  }

 // deserialize PeakPower (one row per network) - - - - - - - - - - - - - - -
 ncVar = group.getVar( "PeakPower" );
 if( ncVar.isNull() )
  v_peak_power.clear();
 else {
  v_peak_power.resize( nn );
  std::vector< size_t > strt = { idx , 0 };
  std::vector< size_t > cnt = { 1 , nn };
  ncVar.getVar( strt , cnt , v_peak_power.data() );
  }

 // deserialize PowerSquilibriumPos- - - - - - - - - - - - - - - - - - - - -
 ncVar = group.getVar( "PowerSquilibriumPos" );
 if( ncVar.isNull() )
  v_power_squilibrium_pos.clear();
 else {
  v_power_squilibrium_pos.resize( ni );
  std::vector< size_t > strt = { start };
  std::vector< size_t > cnt = { ni };
  ncVar.getVar( strt , cnt , v_power_squilibrium_pos.data() );
  }

 // deserialize PowerSquilibriumNeg- - - - - - - - - - - - - - - - - - - - -
 ncVar = group.getVar( "PowerSquilibriumNeg" );
 if( ncVar.isNull() )
  v_power_squilibrium_neg.clear();
 else {
  v_power_squilibrium_neg.resize( ni );
  std::vector< size_t > strt = { start };
  std::vector< size_t > cnt = { ni };
  ncVar.getVar( strt , cnt , v_power_squilibrium_neg.data() );
  }

 // deserialize PowerAggDecPos - - - - - - - - - - - - - - - - - - - - - - - -
 ncVar = group.getVar( "PowerAggDecPos" );
 if( ncVar.isNull() )
  v_power_agg_dec_pos.clear();
 else {
  v_power_agg_dec_pos.resize( ni );
  std::vector< size_t > strt = { start };
  std::vector< size_t > cnt = { ni };
  ncVar.getVar( strt , cnt , v_power_agg_dec_pos.data() );
  }

 // deserialize PowerAggDecNeg - - - - - - - - - - - - - - - - - - - - - - - -
 ncVar = group.getVar( "PowerAggDecNeg" );
 if( ncVar.isNull() )
  v_power_agg_dec_neg.clear();
 else {
  v_power_agg_dec_neg.resize( ni );
  std::vector< size_t > strt = { start };
  std::vector< size_t > cnt = { ni };
  ncVar.getVar( strt , cnt , v_power_agg_dec_neg.data() );
  }

 }  // end( ECNetworkBlockSolution::deserialize( NcGroup & , size_t ) )

/*--------------------------------------------------------------------------*/

void ECNetworkBlockSolution::read( const Block * block )
{
 // call the method of the base class
 NetworkBlockSolution::read( block );

 auto ECNB = dynamic_cast< const ECNetworkBlock * >( block );
 if( ! ECNB )
  throw( std::invalid_argument(
    "ECNetworkBlockSolution::read: block is not an ECNetworkBlock" ) );

 const auto ni = f_number_intervals;
 const auto nn = f_number_nodes;

 auto NCECNB = const_cast< ECNetworkBlock * >( ECNB );

 // read the public injection / absorption variables (if requested)- - - - -
 if( ! v_power_injection.empty() )
  for( Index t = 0 ; t < ni ; ++t ) {
   auto Pi = NCECNB->get_power_injection( t );
   if( ! Pi ) {
    v_power_injection.resize(
     boost::multi_array< double , 2 >::extent_gen()[ 0 ][ 0 ] );
    break;
    }
   for( Index n = 0 ; n < nn ; ++n )
    v_power_injection[ t ][ n ] = Pi[ n ].get_value();
   }

 if( ! v_power_absorption.empty() )
  for( Index t = 0 ; t < ni ; ++t ) {
   auto Pa = NCECNB->get_power_absorption( t );
   if( ! Pa ) {
    v_power_absorption.resize(
     boost::multi_array< double , 2 >::extent_gen()[ 0 ][ 0 ] );
    break;
    }
   for( Index n = 0 ; n < nn ; ++n )
    v_power_absorption[ t ][ n ] = Pa[ n ].get_value();
   }

 // read the shared power (if requested)- - - - - - - - - - - - - - - - - -
 if( ! v_shared_power.empty() ) {
  const auto & SP = ECNB->get_shared_power();
  if( SP.empty() )
   v_shared_power.clear();
  else {
   if( v_shared_power.size() != SP.size() )
    v_shared_power.resize( SP.size() );
   for( Index t = 0 ; t < SP.size() ; ++t )
    v_shared_power[ t ] = SP[ t ].get_value();
   }
  }

 // read the peak power (if requested) - - - - - - - - - - - - - - - - - - -
 if( ! v_peak_power.empty() ) {
  const auto & PP = ECNB->get_peak_power();
  if( PP.empty() )
   v_peak_power.clear();
  else {
   if( v_peak_power.size() != PP.size() )
    v_peak_power.resize( PP.size() );
   for( Index n = 0 ; n < PP.size() ; ++n )
    v_peak_power[ n ] = PP[ n ].get_value();
   }
  }

 // read the squilibrium variables (if requested and present) - - - - - - -
 if( ! v_power_squilibrium_pos.empty() ) {
  const auto & SQp = ECNB->get_power_squilibrium_pos();
  if( SQp.empty() )
   v_power_squilibrium_pos.clear();
  else {
   if( v_power_squilibrium_pos.size() != SQp.size() )
    v_power_squilibrium_pos.resize( SQp.size() );
   for( Index t = 0 ; t < SQp.size() ; ++t )
    v_power_squilibrium_pos[ t ] = SQp[ t ].get_value();
   }
  }

 if( ! v_power_squilibrium_neg.empty() ) {
  const auto & SQn = ECNB->get_power_squilibrium_neg();
  if( SQn.empty() )
   v_power_squilibrium_neg.clear();
  else {
   if( v_power_squilibrium_neg.size() != SQn.size() )
    v_power_squilibrium_neg.resize( SQn.size() );
   for( Index t = 0 ; t < SQn.size() ; ++t )
    v_power_squilibrium_neg[ t ] = SQn[ t ].get_value();
   }
  }

 // read the declared-dispatch variables (if requested and present) - - - - -
 if( ! v_power_agg_dec_pos.empty() ) {
  const auto & ADp = ECNB->get_power_agg_dec_pos();
  if( ADp.empty() )
   v_power_agg_dec_pos.clear();
  else {
   if( v_power_agg_dec_pos.size() != ADp.size() )
    v_power_agg_dec_pos.resize( ADp.size() );
   for( Index t = 0 ; t < ADp.size() ; ++t )
    v_power_agg_dec_pos[ t ] = ADp[ t ].get_value();
   }
  }

 if( ! v_power_agg_dec_neg.empty() ) {
  const auto & ADn = ECNB->get_power_agg_dec_neg();
  if( ADn.empty() )
   v_power_agg_dec_neg.clear();
  else {
   if( v_power_agg_dec_neg.size() != ADn.size() )
    v_power_agg_dec_neg.resize( ADn.size() );
   for( Index t = 0 ; t < ADn.size() ; ++t )
    v_power_agg_dec_neg[ t ] = ADn[ t ].get_value();
   }
  }

 }  // end( ECNetworkBlockSolution::read )

/*--------------------------------------------------------------------------*/

void ECNetworkBlockSolution::write( Block * block )
{
 // call the method of the base class
 NetworkBlockSolution::write( block );

 auto ECNB = dynamic_cast< ECNetworkBlock * >( block );
 if( ! ECNB )
  throw( std::invalid_argument(
    "ECNetworkBlockSolution::write: block is not an ECNetworkBlock" ) );

 const auto ni = f_number_intervals;
 const auto nn = f_number_nodes;

 // write public injection / absorption variables- - - - - - - - - - - - - -
 if( ! v_power_injection.empty() )
  for( Index t = 0 ; t < ni ; ++t ) {
   auto Pi = ECNB->get_power_injection( t );
   if( Pi )
    for( Index n = 0 ; n < nn ; ++n )
     Pi[ n ].set_value( v_power_injection[ t ][ n ] );
   }

 if( ! v_power_absorption.empty() )
  for( Index t = 0 ; t < ni ; ++t ) {
   auto Pa = ECNB->get_power_absorption( t );
   if( Pa )
    for( Index n = 0 ; n < nn ; ++n )
     Pa[ n ].set_value( v_power_absorption[ t ][ n ] );
   }

 // write the shared power- - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_shared_power.empty() ) {
  auto & SP = const_cast< std::vector< ColVariable > & >(
   ECNB->get_shared_power() );
  if( SP.size() != v_shared_power.size() )
   throw( std::invalid_argument(
     "ECNetworkBlockSolution::write: inconsistent shared power size" ) );
  for( Index t = 0 ; t < SP.size() ; ++t )
   SP[ t ].set_value( v_shared_power[ t ] );
  }

 // write the peak power- - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_peak_power.empty() ) {
  auto & PP = const_cast< std::vector< ColVariable > & >(
   ECNB->get_peak_power() );
  if( PP.size() != v_peak_power.size() )
   throw( std::invalid_argument(
     "ECNetworkBlockSolution::write: inconsistent peak power size" ) );
  for( Index n = 0 ; n < PP.size() ; ++n )
   PP[ n ].set_value( v_peak_power[ n ] );
  }

 // write the squilibrium variables - - - - - - - - - - - - - - - - - - - -
 if( ! v_power_squilibrium_pos.empty() ) {
  auto & SQp = const_cast< std::vector< ColVariable > & >(
   ECNB->get_power_squilibrium_pos() );
  if( SQp.size() != v_power_squilibrium_pos.size() )
   throw( std::invalid_argument(
     "ECNetworkBlockSolution::write: inconsistent squilibrium pos size" ) );
  for( Index t = 0 ; t < SQp.size() ; ++t )
   SQp[ t ].set_value( v_power_squilibrium_pos[ t ] );
  }

 if( ! v_power_squilibrium_neg.empty() ) {
  auto & SQn = const_cast< std::vector< ColVariable > & >(
   ECNB->get_power_squilibrium_neg() );
  if( SQn.size() != v_power_squilibrium_neg.size() )
   throw( std::invalid_argument(
     "ECNetworkBlockSolution::write: inconsistent squilibrium neg size" ) );
  for( Index t = 0 ; t < SQn.size() ; ++t )
   SQn[ t ].set_value( v_power_squilibrium_neg[ t ] );
  }

 // write the declared-dispatch variables - - - - - - - - - - - - - - - - - -
 if( ! v_power_agg_dec_pos.empty() ) {
  auto & ADp = const_cast< std::vector< ColVariable > & >(
   ECNB->get_power_agg_dec_pos() );
  if( ADp.size() != v_power_agg_dec_pos.size() )
   throw( std::invalid_argument(
     "ECNetworkBlockSolution::write: inconsistent agg dec pos size" ) );
  for( Index t = 0 ; t < ADp.size() ; ++t )
   ADp[ t ].set_value( v_power_agg_dec_pos[ t ] );
  }

 if( ! v_power_agg_dec_neg.empty() ) {
  auto & ADn = const_cast< std::vector< ColVariable > & >(
   ECNB->get_power_agg_dec_neg() );
  if( ADn.size() != v_power_agg_dec_neg.size() )
   throw( std::invalid_argument(
     "ECNetworkBlockSolution::write: inconsistent agg dec neg size" ) );
  for( Index t = 0 ; t < ADn.size() ; ++t )
   ADn[ t ].set_value( v_power_agg_dec_neg[ t ] );
  }

 }  // end( ECNetworkBlockSolution::write )

/*--------------------------------------------------------------------------*/

void ECNetworkBlockSolution::serialize( netCDF::NcGroup & group ) const
{
 // call the method of the base class: this writes "NumberNodes",
 // "NumberInstants" (if > 1) and "NodeInjection" (if non-empty)
 NetworkBlockSolution::serialize( group );

 // get / create the dimensions we need - - - - - - - - - - - - - - - - - -
 auto nn = group.getDim( "NumberNodes" );  // base class set it
 auto ni = group.getDim( "NumberInstants" );  // present iff > 1

 // serialize PowerInjection - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_power_injection.empty() ) {
  if( ni.isNull() ) {
   // single-instant: store as a 1D variable indexed over nn
   std::vector< double > tmp( f_number_nodes );
   for( Index n = 0 ; n < f_number_nodes ; ++n )
    tmp[ n ] = v_power_injection[ 0 ][ n ];
   ::serialize< double >( group , "PowerInjection" , netCDF::NcDouble() , nn ,
                          tmp );
   }
  else
   ::serialize< double , 2 >( group , "PowerInjection" , netCDF::NcDouble() ,
                              { ni , nn } , v_power_injection );
  }

 // serialize PowerAbsorption- - - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_power_absorption.empty() ) {
  if( ni.isNull() ) {
   std::vector< double > tmp( f_number_nodes );
   for( Index n = 0 ; n < f_number_nodes ; ++n )
    tmp[ n ] = v_power_absorption[ 0 ][ n ];
   ::serialize< double >( group , "PowerAbsorption" , netCDF::NcDouble() , nn ,
                          tmp );
   }
  else
   ::serialize< double , 2 >( group , "PowerAbsorption" , netCDF::NcDouble() ,
                              { ni , nn } , v_power_absorption );
  }

 // serialize SharedPower- - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_shared_power.empty() ) {
  if( ni.isNull() ) {
   // single-instant: store as a scalar
   auto SP = group.addVar( "SharedPower" , netCDF::NcDouble() );
   SP.putVar( v_shared_power.data() );
   }
  else
   ::serialize< double >( group , "SharedPower" , netCDF::NcDouble() , ni ,
                          v_shared_power );
  }

 // serialize PeakPower- - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_peak_power.empty() )
  ::serialize< double >( group , "PeakPower" , netCDF::NcDouble() , nn ,
                         v_peak_power );

 // serialize PowerSquilibriumPos / Neg- - - - - - - - - - - - - - - - - - -
 if( ! v_power_squilibrium_pos.empty() ) {
  if( ni.isNull() ) {
   auto SP = group.addVar( "PowerSquilibriumPos" , netCDF::NcDouble() );
   SP.putVar( v_power_squilibrium_pos.data() );
   }
  else
   ::serialize< double >( group , "PowerSquilibriumPos" , netCDF::NcDouble() ,
                          ni , v_power_squilibrium_pos );
  }

 if( ! v_power_squilibrium_neg.empty() ) {
  if( ni.isNull() ) {
   auto SP = group.addVar( "PowerSquilibriumNeg" , netCDF::NcDouble() );
   SP.putVar( v_power_squilibrium_neg.data() );
   }
  else
   ::serialize< double >( group , "PowerSquilibriumNeg" , netCDF::NcDouble() ,
                          ni , v_power_squilibrium_neg );
  }

 // serialize PowerAggDecPos / Neg - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_power_agg_dec_pos.empty() ) {
  if( ni.isNull() ) {
   auto SP = group.addVar( "PowerAggDecPos" , netCDF::NcDouble() );
   SP.putVar( v_power_agg_dec_pos.data() );
   }
  else
   ::serialize< double >( group , "PowerAggDecPos" , netCDF::NcDouble() ,
                          ni , v_power_agg_dec_pos );
  }

 if( ! v_power_agg_dec_neg.empty() ) {
  if( ni.isNull() ) {
   auto SP = group.addVar( "PowerAggDecNeg" , netCDF::NcDouble() );
   SP.putVar( v_power_agg_dec_neg.data() );
   }
  else
   ::serialize< double >( group , "PowerAggDecNeg" , netCDF::NcDouble() ,
                          ni , v_power_agg_dec_neg );
  }

 }  // end( ECNetworkBlockSolution::serialize( NcGroup & ) )

/*--------------------------------------------------------------------------*/

void ECNetworkBlockSolution::serialize( netCDF::NcGroup & group ,
                                        size_t idx ) const
{
 // call the method of the base class: this fills/checks NumberNetworks,
 // NumberNodes, TotalNumberInstants, EndInstant, NodeInjection (if any)
 NetworkBlockSolution::serialize( group , idx );

 // recover dimensions; NumberNetworks and NumberNodes must exist
 auto nnw = group.getDim( "NumberNetworks" );
 auto nn = group.getDim( "NumberNodes" );

 // for time-indexed EC variables we use TotalNumberInstants when defined,
 // and fall back to NumberNetworks when each block covers exactly one
 // instant (i.e. when the parent did not create TotalNumberInstants /
 // EndInstant)
 auto tni = group.getDim( "TotalNumberInstants" );
 const bool has_tni = ! tni.isNull();
 auto tdim = has_tni ? tni : nnw;

 // recover the offset in time-indexed variables, exactly as the base class
 // does: 0 for idx == 0, else EndInstant[ idx - 1 ] when present, else idx
 size_t start;
 if( has_tni ) {
  auto EI = group.getVar( "EndInstant" );  // base class created it already
  start = 0;
  if( idx > 0 ) {
   std::vector< size_t > vidx = { idx - 1 };
   EI.getVar( vidx , & start );
   }
  }
 else
  start = idx;

 // create/retrieve the netCDF variables - - - - - - - - - - - - - - - - - -
 netCDF::NcVar PI;  // PowerInjection
 netCDF::NcVar PA;  // PowerAbsorption
 netCDF::NcVar SP;  // SharedPower
 netCDF::NcVar PP;  // PeakPower
 netCDF::NcVar SQp; // PowerSquilibriumPos
 netCDF::NcVar SQn; // PowerSquilibriumNeg
 netCDF::NcVar ADp; // PowerAggDecPos
 netCDF::NcVar ADn; // PowerAggDecNeg

 if( idx == 0 ) {  // first call: initialize variables
  if( ! v_power_injection.empty() )
   PI = group.addVar( "PowerInjection" , netCDF::NcDouble() , { tdim , nn } );

  if( ! v_power_absorption.empty() )
   PA = group.addVar( "PowerAbsorption" , netCDF::NcDouble() , { tdim , nn } );

  if( ! v_shared_power.empty() )
   SP = group.addVar( "SharedPower" , netCDF::NcDouble() , { tdim } );

  if( ! v_peak_power.empty() )
   PP = group.addVar( "PeakPower" , netCDF::NcDouble() , { nnw , nn } );

  if( ! v_power_squilibrium_pos.empty() )
   SQp = group.addVar( "PowerSquilibriumPos" , netCDF::NcDouble() ,
                       { tdim } );

  if( ! v_power_squilibrium_neg.empty() )
   SQn = group.addVar( "PowerSquilibriumNeg" , netCDF::NcDouble() ,
                       { tdim } );

  if( ! v_power_agg_dec_pos.empty() )
   ADp = group.addVar( "PowerAggDecPos" , netCDF::NcDouble() , { tdim } );

  if( ! v_power_agg_dec_neg.empty() )
   ADn = group.addVar( "PowerAggDecNeg" , netCDF::NcDouble() , { tdim } );
  }
 else {  // subsequent call: read what is supposedly already there
  if( ! v_power_injection.empty() )
   PI = group.getVar( "PowerInjection" );

  if( ! v_power_absorption.empty() )
   PA = group.getVar( "PowerAbsorption" );

  if( ! v_shared_power.empty() )
   SP = group.getVar( "SharedPower" );

  if( ! v_peak_power.empty() )
   PP = group.getVar( "PeakPower" );

  if( ! v_power_squilibrium_pos.empty() )
   SQp = group.getVar( "PowerSquilibriumPos" );

  if( ! v_power_squilibrium_neg.empty() )
   SQn = group.getVar( "PowerSquilibriumNeg" );

  if( ! v_power_agg_dec_pos.empty() )
   ADp = group.getVar( "PowerAggDecPos" );

  if( ! v_power_agg_dec_neg.empty() )
   ADn = group.getVar( "PowerAggDecNeg" );
  }

 // write the data - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 std::vector< size_t > strt2 = { start , 0 };
 std::vector< size_t > cnt2 = { f_number_intervals , f_number_nodes };
 std::vector< size_t > strt1 = { start };
 std::vector< size_t > cnt1 = { f_number_intervals };
 std::vector< size_t > strtN = { idx , 0 };
 std::vector< size_t > cntN = { 1 , f_number_nodes };

 if( ! PI.isNull() )
  PI.putVar( strt2 , cnt2 , v_power_injection.data() );

 if( ! PA.isNull() )
  PA.putVar( strt2 , cnt2 , v_power_absorption.data() );

 if( ! SP.isNull() )
  SP.putVar( strt1 , cnt1 , v_shared_power.data() );

 if( ! PP.isNull() )
  PP.putVar( strtN , cntN , v_peak_power.data() );

 if( ! SQp.isNull() )
  SQp.putVar( strt1 , cnt1 , v_power_squilibrium_pos.data() );

 if( ! SQn.isNull() )
  SQn.putVar( strt1 , cnt1 , v_power_squilibrium_neg.data() );

 if( ! ADp.isNull() )
  ADp.putVar( strt1 , cnt1 , v_power_agg_dec_pos.data() );

 if( ! ADn.isNull() )
  ADn.putVar( strt1 , cnt1 , v_power_agg_dec_neg.data() );

 }  // end( ECNetworkBlockSolution::serialize( NcGroup & , size_t ) )

/*--------------------------------------------------------------------------*/

ECNetworkBlockSolution * ECNetworkBlockSolution::scale( double factor ) const
{
 // call the method of the base class, which calls clone() and therefore
 // returns an ECNetworkBlockSolution
 auto sol = dynamic_cast< ECNetworkBlockSolution * >(
                                     NetworkBlockSolution::scale( factor ) );
 assert( sol );

 if( factor == 1 )
  return( sol );

 if( ! v_power_injection.empty() )
  for( Index t = 0 ; t < f_number_intervals ; ++t )
   for( Index n = 0 ; n < f_number_nodes ; ++n )
    sol->v_power_injection[ t ][ n ] *= factor;

 if( ! v_power_absorption.empty() )
  for( Index t = 0 ; t < f_number_intervals ; ++t )
   for( Index n = 0 ; n < f_number_nodes ; ++n )
    sol->v_power_absorption[ t ][ n ] *= factor;

 for( auto & x : sol->v_shared_power ) x *= factor;
 for( auto & x : sol->v_peak_power ) x *= factor;
 for( auto & x : sol->v_power_squilibrium_pos ) x *= factor;
 for( auto & x : sol->v_power_squilibrium_neg ) x *= factor;
 for( auto & x : sol->v_power_agg_dec_pos ) x *= factor;
 for( auto & x : sol->v_power_agg_dec_neg ) x *= factor;

 return( sol );

 }  // end( ECNetworkBlockSolution::scale )

/*--------------------------------------------------------------------------*/

void ECNetworkBlockSolution::sum( const Solution * solution ,
                                  double multiplier )
{
 // call the method of the base class
 NetworkBlockSolution::sum( solution , multiplier );

 auto ECNBS = dynamic_cast< const ECNetworkBlockSolution * >( solution );
 if( ! ECNBS )
  throw( std::invalid_argument(
    "ECNetworkBlockSolution::sum: solution not an ECNetworkBlockSolution" ) );

 // sanity checks - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( v_power_injection.shape()[ 0 ] != ECNBS->v_power_injection.shape()[ 0 ]
  || v_power_injection.shape()[ 1 ] != ECNBS->v_power_injection.shape()[ 1 ] )
  throw( std::invalid_argument(
    "ECNetworkBlockSolution::sum: inconsistent power injection size" ) );

 if( v_power_absorption.shape()[ 0 ] != ECNBS->v_power_absorption.shape()[ 0 ]
  || v_power_absorption.shape()[ 1 ] !=
     ECNBS->v_power_absorption.shape()[ 1 ] )
  throw( std::invalid_argument(
    "ECNetworkBlockSolution::sum: inconsistent power absorption size" ) );

 if( v_shared_power.size() != ECNBS->v_shared_power.size() )
  throw( std::invalid_argument(
    "ECNetworkBlockSolution::sum: inconsistent shared power size" ) );

 if( v_peak_power.size() != ECNBS->v_peak_power.size() )
  throw( std::invalid_argument(
    "ECNetworkBlockSolution::sum: inconsistent peak power size" ) );

 if( v_power_squilibrium_pos.size() != ECNBS->v_power_squilibrium_pos.size() )
  throw( std::invalid_argument(
    "ECNetworkBlockSolution::sum: inconsistent squilibrium pos size" ) );

 if( v_power_squilibrium_neg.size() != ECNBS->v_power_squilibrium_neg.size() )
  throw( std::invalid_argument(
    "ECNetworkBlockSolution::sum: inconsistent squilibrium neg size" ) );

 if( v_power_agg_dec_pos.size() != ECNBS->v_power_agg_dec_pos.size() )
  throw( std::invalid_argument(
    "ECNetworkBlockSolution::sum: inconsistent agg dec pos size" ) );

 if( v_power_agg_dec_neg.size() != ECNBS->v_power_agg_dec_neg.size() )
  throw( std::invalid_argument(
    "ECNetworkBlockSolution::sum: inconsistent agg dec neg size" ) );

 // accumulate - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_power_injection.empty() )
  for( Index t = 0 ; t < f_number_intervals ; ++t )
   for( Index n = 0 ; n < f_number_nodes ; ++n )
    v_power_injection[ t ][ n ] +=
     ECNBS->v_power_injection[ t ][ n ] * multiplier;

 if( ! v_power_absorption.empty() )
  for( Index t = 0 ; t < f_number_intervals ; ++t )
   for( Index n = 0 ; n < f_number_nodes ; ++n )
    v_power_absorption[ t ][ n ] +=
     ECNBS->v_power_absorption[ t ][ n ] * multiplier;

 for( std::size_t t = 0 ; t < v_shared_power.size() ; ++t )
  v_shared_power[ t ] += ECNBS->v_shared_power[ t ] * multiplier;

 for( std::size_t n = 0 ; n < v_peak_power.size() ; ++n )
  v_peak_power[ n ] += ECNBS->v_peak_power[ n ] * multiplier;

 for( std::size_t t = 0 ; t < v_power_squilibrium_pos.size() ; ++t )
  v_power_squilibrium_pos[ t ] +=
   ECNBS->v_power_squilibrium_pos[ t ] * multiplier;

 for( std::size_t t = 0 ; t < v_power_squilibrium_neg.size() ; ++t )
  v_power_squilibrium_neg[ t ] +=
   ECNBS->v_power_squilibrium_neg[ t ] * multiplier;

 for( std::size_t t = 0 ; t < v_power_agg_dec_pos.size() ; ++t )
  v_power_agg_dec_pos[ t ] +=
   ECNBS->v_power_agg_dec_pos[ t ] * multiplier;

 for( std::size_t t = 0 ; t < v_power_agg_dec_neg.size() ; ++t )
  v_power_agg_dec_neg[ t ] +=
   ECNBS->v_power_agg_dec_neg[ t ] * multiplier;

 }  // end( ECNetworkBlockSolution::sum )

/*--------------------------------------------------------------------------*/

ECNetworkBlockSolution * ECNetworkBlockSolution::clone( bool empty ) const
{
 auto sol = new ECNetworkBlockSolution();

 if( ! empty ) {
  NetworkBlockSolution::guts_of_clone( sol );

  copy_multi_array( sol->v_power_injection , v_power_injection );
  copy_multi_array( sol->v_power_absorption , v_power_absorption );
  sol->v_shared_power = v_shared_power;
  sol->v_peak_power = v_peak_power;
  sol->v_power_squilibrium_pos = v_power_squilibrium_pos;
  sol->v_power_squilibrium_neg = v_power_squilibrium_neg;
  sol->v_power_agg_dec_pos = v_power_agg_dec_pos;
  sol->v_power_agg_dec_neg = v_power_agg_dec_neg;
  }

 return( sol );

 }  // end( ECNetworkBlockSolution::clone )

/*--------------------------------------------------------------------------*/
/*----------------------- End File ECNetworkBlock.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
