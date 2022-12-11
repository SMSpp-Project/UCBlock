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
 * \copyright &copy by Antonio Frangioni
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
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

ECNetworkBlock::~ECNetworkBlock() {

 Constraint::clear( micro_power_balance_const );
 Constraint::clear( power_balance_const );
 Constraint::clear( power_flow_limit_const );

 Constraint::clear( node_injection_upper_bound_const );

 objective.clear();
}

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void ECNetworkBlock::ECNetworkData::deserialize(
 const netCDF::NcGroup & group ) {

#ifndef NDEBUG
 static std::vector< std::string > expected_dims = { "NumberNodes" ,
                                                     "NumberIntervals" };
 check_dimensions( group , expected_dims , std::cerr );

 static std::vector< std::string > expected_vars = { "ActiveDemand" ,
                                                     "BuyPrice" ,
                                                     "SellPrice" ,
                                                     "RewardPrice" ,
                                                     "PeakTariff" ,
                                                     "ConstTerm" ,
                                                     "MaxNodeInjection" };
 check_variables( group , expected_vars , std::cerr );
#endif

 // Mandatory variables

 ::deserialize_dim( group , "NumberNodes" , f_number_nodes , false );
 if( f_number_nodes == 1 )
  throw( std::invalid_argument( "ECNetworkBlock::deserialize: cannot create "
                                "an Energy Community with just one user" ) );
}  // end( ECNetworkBlock::ECNetworkData::deserialize )

/*--------------------------------------------------------------------------*/

void ECNetworkBlock::deserialize( const netCDF::NcGroup & group ) {

#ifndef NDEBUG
 static std::vector< std::string > expected_dims = { "NumberNodes" ,
                                                     "NumberIntervals" };
 check_dimensions( group , expected_dims , std::cerr );

 static std::vector< std::string > expected_vars = { "ActiveDemand" ,
                                                     "BuyPrice" ,
                                                     "SellPrice" ,
                                                     "RewardPrice" ,
                                                     "PeakTariff" ,
                                                     "ConstTerm" ,
                                                     "MaxNodeInjection" };
 check_variables( group , expected_vars , std::cerr );
#endif

 // Mandatory variables

 ::deserialize_dim( group , "NumberIntervals" , f_number_intervals , false );

 ::deserialize( group , "BuyPrice" , f_number_intervals , v_BuyPrice ,
                false , true );
 if( v_BuyPrice.size() == 1 )
  v_BuyPrice.resize( f_number_intervals , v_BuyPrice[ 0 ] );

 ::deserialize( group , "SellPrice" , f_number_intervals , v_SellPrice ,
                false , true );
 if( v_SellPrice.size() == 1 )
  v_SellPrice.resize( f_number_intervals , v_SellPrice[ 0 ] );

 ::deserialize( group , "RewardPrice" , f_number_intervals , v_RewardPrice ,
                false , true );
 if( v_RewardPrice.size() == 1 )
  v_RewardPrice.resize( f_number_intervals , v_RewardPrice[ 0 ] );

 ::deserialize( group , f_PeakTariff , "PeakTariff" , false );

 // Optional variables

 Index NumberNodes;
 if( ::deserialize_dim( group , "NumberNodes" , NumberNodes , true ) ) {
  ::deserialize( group , "ActiveDemand" , v_ActiveDemand );
  // always check if the demand is given in the correct shape
  assert( ( v_ActiveDemand.shape()[ 0 ] == f_number_intervals ) &&
          ( v_ActiveDemand.shape()[ 1 ] == NumberNodes ) );
 }

 // it is mandatory ONLY IF we use a Solver that optimize each Block at a
 // time to lower bound the node injection; by default it is set in
 // UCBlock::generate_node_injection_constraints() as the sum of all the
 // maximum powers of the UnitBlock of the problem
 ::deserialize( group , "MaxNodeInjection" , v_MaxNodeInjection );

 ::deserialize( group , f_ConstTerm , "ConstTerm" );
}  // end( ECNetworkBlock::deserialize )

/*--------------------------------------------------------------------------*/
/*--------- METHODS FOR LOADING, PRINTING & SAVING THE ECNetworkBlock ------*/
/*--------------------------------------------------------------------------*/

void ECNetworkBlock::ECNetworkData::serialize( netCDF::NcGroup & group ) const {

 NetworkBlock::NetworkData::serialize( group );

}  // end( ECNetworkBlock::ECNetworkData::serialize )

/*--------------------------------------------------------------------------*/

void ECNetworkBlock::serialize( netCDF::NcGroup & group ) const {

 NetworkBlock::serialize( group );

 ::serialize( group , "PeakTariff" , netCDF::NcDouble() , f_PeakTariff );

 auto NumberIntervals = group.getDim( "NumberIntervals" );

 ::serialize( group , "BuyPrice" , netCDF::NcDouble() , NumberIntervals ,
              v_BuyPrice );

 ::serialize( group , "SellPrice" , netCDF::NcDouble() , NumberIntervals ,
              v_SellPrice );

 ::serialize( group , "RewardPrice" , netCDF::NcDouble() , NumberIntervals ,
              v_RewardPrice );

 if( auto network_data = get_NetworkData() )
  // If an ECNetworkData is present, serialize it.
  network_data->serialize( group );

 if( ! v_ActiveDemand.empty() ) {
  // This ECNetworkBlock has active demand, so it is serialized.

  auto NumberNodes = group.getDim( "NumberNodes" );

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
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

void ECNetworkBlock::generate_abstract_variables( Configuration * stvv ) {

 if( variables_generated() )
  return; // variables have already been generated

 auto number_nodes = get_number_nodes();
 auto number_intervals = get_number_intervals();

 // the node injection variables
 v_node_injection.resize( boost::extents[ number_intervals ][ number_nodes ] );
 for( Index t = 0 ; t < number_intervals ; ++t )
  for( Index node_id = 0 ; node_id < number_nodes ; ++node_id )
   v_node_injection[ t ][ node_id ].set_type( ColVariable::kContinuous );
 add_static_variable( v_node_injection , "S" );

 // the power injected variables
 v_micro_power_injection.resize( boost::extents[ number_intervals ][ number_nodes ] );
 for( Index t = 0 ; t < number_intervals ; ++t )
  for( Index node_id = 0 ; node_id < number_nodes ; ++node_id )
   v_micro_power_injection[ t ][ node_id ].set_type( ColVariable::kNonNegative );
 add_static_variable( v_micro_power_injection , "micro_power_injection" );

 // the power absorbed variables
 v_micro_power_absorption.resize( boost::extents[ number_intervals ][ number_nodes ] );
 for( Index t = 0 ; t < number_intervals ; ++t )
  for( Index node_id = 0 ; node_id < number_nodes ; ++node_id )
   v_micro_power_absorption[ t ][ node_id ].set_type( ColVariable::kNonNegative );
 add_static_variable( v_micro_power_absorption , "micro_power_absorption" );

 // the power injected variables
 v_public_power_injection.resize( boost::extents[ number_intervals ][ number_nodes ] );
 for( Index t = 0 ; t < number_intervals ; ++t )
  for( Index node_id = 0 ; node_id < number_nodes ; ++node_id )
   v_public_power_injection[ t ][ node_id ].set_type( ColVariable::kNonNegative );
 add_static_variable( v_public_power_injection , "public_power_injection" );

 // the power absorbed variables
 v_public_power_absorption.resize( boost::extents[ number_intervals ][ number_nodes ] );
 for( Index t = 0 ; t < number_intervals ; ++t )
  for( Index node_id = 0 ; node_id < number_nodes ; ++node_id )
   v_public_power_absorption[ t ][ node_id ].set_type( ColVariable::kNonNegative );
 add_static_variable( v_public_power_absorption , "public_power_absorption" );

 // the peak power variables
 v_peak_power.resize( number_nodes );
 for( auto & var : v_peak_power )
  var.set_type( ColVariable::kNonNegative );
 add_static_variable( v_peak_power , "peak_power" );

 set_variables_generated();
}  // end( ECNetworkBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void ECNetworkBlock::generate_abstract_constraints( Configuration * stcc ) {

 if( constraints_generated() )
  return; // constraints have already been generated

 auto number_nodes = get_number_nodes();
 auto number_intervals = get_number_intervals();

/*------------------------- inequality constraints -------------------------*/

 // set that the hourly dispatch cannot go beyond the maximum dispatch
 // of the corresponding peak power period, i.e.:
 //
 //    P^{max} >= P^{POD,+} - P^{POD,-}       for all u, t       (1)
 // => P^{POD,+} - P^{POD,-} - P^{max} <= 0   for all u, t
 //
 //    P^{max} >= - [ P^{POD,+} - P^{POD,-} ]   for all u, t     (2)
 // => - P^{POD,+} + P^{POD,-} - P^{max} <= 0   for all u, t
 //
 // where:
 //
 //    P^{POD,+} = P^{P,+} + P^{M,+}
 //    P^{POD,-} = P^{P,-} + P^{M,-}

 power_flow_limit_const.resize(
  boost::multi_array< FRowConstraint , 3 >::extent_gen()
  [ number_nodes ][ number_intervals ][ 2 ] ); // 2 dims, i.e., the sign (+/-)

 for( Index t = 0 ; t < number_intervals ; ++t )

  for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {

   LinearFunction::v_coeff_pair vars_p;
   LinearFunction::v_coeff_pair vars_n;

   // P^{max} vars also depends from P^{M+} and P^{M-}, i.e., the injection
   // and absorption from the microgrid, to give an economic benefit to users
   // that do not contribute to the community by sharing energy since they
   // are unable to install assets due to economic or space reasons;
   // and on which, otherwise, all the costs of the peak powers would be borne
   {
    // case (1)
    vars_p.push_back( std::make_pair( &v_micro_power_injection[ t ][ node_id ] ,
                                      1.0 ) );
    vars_p.push_back( std::make_pair( &v_micro_power_absorption[ t ][ node_id ] ,
                                      -1.0 ) );

    // case (2)
    vars_n.push_back( std::make_pair( &v_micro_power_injection[ t ][ node_id ] ,
                                      -1.0 ) );
    vars_n.push_back( std::make_pair( &v_micro_power_absorption[ t ][ node_id ] ,
                                      1.0 ) );
   }

   // case (1)
   vars_p.push_back( std::make_pair( &v_public_power_injection[ t ][ node_id ] ,
                                     1.0 ) );
   vars_p.push_back( std::make_pair( &v_public_power_absorption[ t ][ node_id ] ,
                                     -1.0 ) );
   vars_p.push_back( std::make_pair( &v_peak_power[ node_id ] , 1.0 ) );

   // case (2)
   vars_n.push_back( std::make_pair( &v_public_power_injection[ t ][ node_id ] ,
                                     -1.0 ) );
   vars_n.push_back( std::make_pair( &v_public_power_absorption[ t ][ node_id ] ,
                                     1.0 ) );
   vars_n.push_back( std::make_pair( &v_peak_power[ node_id ] , 1.0 ) );

   // case (1)
   power_flow_limit_const[ node_id ][ t ][ 0 ].set_rhs( Inf< double >() );
   power_flow_limit_const[ node_id ][ t ][ 0 ].set_lhs( 0.0 );
   power_flow_limit_const[ node_id ][ t ][ 0 ].set_function(
    new LinearFunction( std::move( vars_p ) ) );

   // case (2)
   power_flow_limit_const[ node_id ][ t ][ 1 ].set_rhs( Inf< double >() );
   power_flow_limit_const[ node_id ][ t ][ 1 ].set_lhs( 0.0 );
   power_flow_limit_const[ node_id ][ t ][ 1 ].set_function(
    new LinearFunction( std::move( vars_n ) ) );
  }

 add_static_constraint( power_flow_limit_const , "Power_Flow_Limit_Const" );

/*-------------------------- equality constraints --------------------------*/

 // set the power balance within the microgrid market/network, i.e., the
 // flows within the microgrid to have sum equal to zero:
 //
 //    P^{M,+} = P^{M,-}       for all t
 // => P^{M,+} - P^{M,-} = 0   for all t

 if( micro_power_balance_const.size() != number_intervals ) {
  assert( micro_power_balance_const.empty() );
  micro_power_balance_const.resize( number_intervals );
 }

 for( Index t = 0 ; t < number_intervals ; ++t ) {

  LinearFunction::v_coeff_pair vars;

  for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {

   vars.push_back( std::make_pair( &v_micro_power_injection[ t ][ node_id ] ,
                                   1.0 ) );
   vars.push_back( std::make_pair( &v_micro_power_absorption[ t ][ node_id ] ,
                                   -1.0 ) );
  }

  micro_power_balance_const[ t ].set_both( 0.0 );
  micro_power_balance_const[ t ].set_function(
   new LinearFunction( std::move( vars ) ) );
 }

 add_static_constraint( micro_power_balance_const ,
                        "Micro_Power_Balance_Const" );

 // set the power balance, i.e.:
 //
 //    P^{POD,+} - P^{POD,-} - node_injection = - active_demand   for all t, u
 //
 // where:
 //
 //    P^{POD,+} = P^{P,+} + P^{M,+}
 //    P^{POD,-} = P^{P,-} + P^{M,-}

 power_balance_const.resize(
  boost::multi_array< FRowConstraint , 2 >::extent_gen()
  [ number_nodes ][ number_intervals ] );

 for( Index t = 0 ; t < number_intervals ; ++t )

  for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {

   LinearFunction::v_coeff_pair vars;

   vars.push_back( std::make_pair( &v_public_power_injection[ t ][ node_id ] ,
                                   1.0 ) );
   vars.push_back( std::make_pair( &v_micro_power_injection[ t ][ node_id ] ,
                                   1.0 ) );
   vars.push_back( std::make_pair( &v_public_power_absorption[ t ][ node_id ] ,
                                   -1.0 ) );
   vars.push_back( std::make_pair( &v_micro_power_absorption[ t ][ node_id ] ,
                                   -1.0 ) );
   vars.push_back( std::make_pair( &v_node_injection[ t ][ node_id ] , -1.0 ) );
   power_balance_const[ node_id ][ t ].set_both(
    -v_ActiveDemand[ t ][ node_id ] );
   power_balance_const[ node_id ][ t ].set_function(
    new LinearFunction( std::move( vars ) ) );
  }

 add_static_constraint( power_balance_const , "Power_Balance_Const" );

 // node injection upper bound constraints

 node_injection_upper_bound_const.resize(
  boost::multi_array< FRowConstraint , 2 >::extent_gen()
  [ number_nodes ][ number_intervals ] );

 for( Index t = 0 ; t < number_intervals ; ++t )

  for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {

   node_injection_upper_bound_const[ node_id ][ t ].set_lhs( -Inf< double >() );
   node_injection_upper_bound_const[ node_id ][ t ].set_rhs(
    v_MaxNodeInjection[ t ][ node_id ] );
   node_injection_upper_bound_const[ node_id ][ t ].set_variable(
    &v_node_injection[ t ][ node_id ] );
  }

 add_static_constraint( node_injection_upper_bound_const ,
                        "Node_Injection_Upper_Bound_Const" );

 set_constraints_generated();

}  // end( ECNetworkBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

bool ECNetworkBlock::is_feasible( bool useabstract , Configuration * fsbc ) {

 // Retrieve the tolerance.

 auto config = dynamic_cast< SimpleConfiguration< double > * >( fsbc );

 if( ( ! config ) && f_BlockConfig )
  config = dynamic_cast< SimpleConfiguration< double > * >
  ( f_BlockConfig->f_is_feasible_Configuration );

 // If a tolerance has not been provided, use the default tolerance.
 const auto tol = config ? config->f_value : 1.0e-8;

 return( NetworkBlock::is_feasible( useabstract )
         && Constraint::is_feasible( micro_power_balance_const , tol )
         && Constraint::is_feasible( power_balance_const , tol )
         && Constraint::is_feasible( power_flow_limit_const , tol )
         && Constraint::is_feasible( node_injection_upper_bound_const , tol ) );

}  // end( ECNetworkBlock::is_feasible )

/*--------------------------------------------------------------------------*/

void ECNetworkBlock::generate_objective( Configuration * objc ) {

 if( objective_generated() )
  return; // objective has already been generated

 if( get_objective() != nullptr )  // an objective is there already
  return;                          // cowardly (and silently) return

 LinearFunction::v_coeff_pair vars;

 for( Index node_id = 0 ; node_id < get_number_nodes() ; ++node_id ) {

  for( Index t = 0 ; t < get_number_intervals() ; ++t ) {

   // R_{j}^{U,P}, i.e., the net economic balance wrt the public market
   vars.push_back( std::make_pair( &v_public_power_absorption[ t ][ node_id ] ,
                                   v_BuyPrice[ t ] ) );
   vars.push_back( std::make_pair( &v_micro_power_absorption[ t ][ node_id ] ,
                                   // ECR_{j}, i.e., the reward awarded to the community
                                   v_BuyPrice[ t ] - v_RewardPrice[ t ] ) );
   vars.push_back( std::make_pair( &v_public_power_injection[ t ][ node_id ] ,
                                   -v_SellPrice[ t ] ) );
   vars.push_back( std::make_pair( &v_micro_power_injection[ t ][ node_id ] ,
                                   -v_SellPrice[ t ] ) );
  }

  // C_{j}^{U,P}, i.e., the costs due to the peak power
  vars.push_back( std::make_pair( &v_peak_power[ node_id ] , f_PeakTariff ) );
 }

 auto lf = new LinearFunction( std::move( vars ) );

 lf->set_constant_term( f_ConstTerm );

 objective.set_function( lf );
 objective.set_sense( Objective::eMin );

 // set block objective
 this->set_objective( &objective );

 set_objective_generated();

}  // end( ECNetworkBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

void ECNetworkBlock::set_active_demand(
 // TODO this should be a const ptr to double
 std::vector< double >::const_iterator values , Block::Subset && subset ,
 const bool ordered , c_ModParam issuePMod , c_ModParam issueAMod ) {
 // TODO
}

/*--------------------------------------------------------------------------*/

void ECNetworkBlock::set_active_demand(
 // TODO this should be a const ptr to double
 std::vector< double >::const_iterator values ,
 Block::Range rng , c_ModParam issuePMod , c_ModParam issueAMod ) {
 // TODO
}

/*--------------------------------------------------------------------------*/
/*----------------------- End File ECNetworkBlock.cpp ----------------------*/
/*--------------------------------------------------------------------------*/