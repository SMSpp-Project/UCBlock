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

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

ECNetworkBlock::~ECNetworkBlock() {

 clear_constraints( micro_power_balance_constraints );

 clear_constraints( power_balance_constraints );
 clear_constraints( power_flow_limit_constraints );

 objective.clear();

 // Delete the ECNetworkData if it is local.
 if( f_local_NetworkData )
  delete f_NetworkData;
}

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void ECNetworkBlock::ECNetworkData::deserialize(
 const netCDF::NcGroup & group ) {

 NetworkBlock::NetworkData::deserialize( group );

#ifndef NDEBUG
 static std::vector< std::string > expected_dims = { "NumberNodes" ,
                                                     "NumberLines" ,
                                                     "NumberIntervals" };
 check_dimensions( group , expected_dims , std::cerr );

 static std::vector< std::string > expected_vars = { "StartLine" , "EndLine" ,
                                                     "BuyPrice" ,
                                                     "SellPrice" ,
                                                     "MaxTariff" };
 check_variables( group , expected_vars , std::cerr );
#endif

 // Mandatory variables

 ::deserialize_dim( group , "NumberNodes" , f_number_nodes , false );

 ::deserialize( group , "BuyPrice" , get_number_intervals() ,
                v_buy_price , false , true );

 ::deserialize( group , "SellPrice" , get_number_intervals() ,
                v_sell_price , false , true );

 ::deserialize( group , f_max_tariff , "MaxTariff" , false );
}

/*--------------------------------------------------------------------------*/

void ECNetworkBlock::deserialize( const netCDF::NcGroup & group ) {

 NetworkBlock::deserialize( group );

#ifndef NDEBUG
 static std::vector< std::string > expected_dims = { "NumberNodes" };
 check_dimensions( group , expected_dims , std::cerr );

 static std::vector< std::string > expected_vars = { "ActiveDemand" ,
                                                     "ConstTerm" };
 check_variables( group , expected_vars , std::cerr );
#endif

 delete f_NetworkData;
 f_NetworkData = new ECNetworkData();
 f_NetworkData->deserialize( group );
 f_local_NetworkData = true;

 // Optional variables

 ::deserialize( group , "ActiveDemand" , v_active_demand , true );
}

/*--------------------------------------------------------------------------*/
/*--------- METHODS FOR LOADING, PRINTING & SAVING THE DCNetworkBlock ------*/
/*--------------------------------------------------------------------------*/

void ECNetworkBlock::ECNetworkData::serialize( netCDF::NcGroup & group ) const {

 NetworkBlock::NetworkData::serialize( group );

 ::serialize( group , "MaxTariff" , netCDF::NcDouble() , f_max_tariff );

 auto NumberIntervals = group.getDim( "NumberIntervals" );

 ::serialize( group , "BuyPrice" , netCDF::NcDouble() , NumberIntervals ,
              v_buy_price );

 ::serialize( group , "SellPrice" , netCDF::NcDouble() , NumberIntervals ,
              v_sell_price );
}

/*--------------------------------------------------------------------------*/

void ECNetworkBlock::serialize( netCDF::NcGroup & group ) const {

 NetworkBlock::serialize( group );

 if( auto network_data = get_NetworkData() )
  // If an ECNetworkData is present, serialize it.
  network_data->serialize( group );

 if( !v_active_demand.empty() ) {
  // This DCNetworkBlock has active demand, so it is serialized.

  auto NumberNodes = group.getDim( "NumberNodes" );

  if( NumberNodes.isNull() ) {
   /* The dimension "NumberNodes" is not present in the group (which means
    * that an ECNetworkData is not present). However, the number of nodes can
    * still be obtained from the size of the active demand vector. Notice that
    * the name "NumberNodes" is not used for this new dimension, because it
    * would indicate that an ECNetworkData is present (which is not the
    * case). Therefore, we create an alternative dimension in order to be able
    * to serialize the active demand. */
   NumberNodes = group.addDim( "__NumberNodes__" , v_active_demand.size() );
  }

  auto NumberIntervals = group.getDim( "NumberIntervals" );

  // Finally, serialize the active demand.
  ::serialize( group , "ActiveDemand" , netCDF::NcDouble() ,
               { NumberIntervals , NumberNodes } , v_active_demand );
 }
}

/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

void ECNetworkBlock::generate_abstract_variables(
 Configuration * stvv ) {

 if( variables_generated() )
  return; // variables have already been generated

 auto number_nodes = get_number_nodes();

 // the power injected variables
 v_micro_power_injection.resize( number_nodes );
 for( auto & var : v_micro_power_injection )
  var.set_type( ColVariable::kNonNegative );
 add_static_variable( v_micro_power_injection , "micro_power_injection" );

 // the power absorbed variables
 v_micro_power_absorption.resize( number_nodes );
 for( auto & var : v_micro_power_absorption )
  var.set_type( ColVariable::kNonNegative );
 add_static_variable( v_micro_power_absorption , "micro_power_absorption" );

 // the power injected variables
 v_public_power_injection.resize( number_nodes );
 for( auto & var : v_public_power_injection )
  var.set_type( ColVariable::kNonNegative );
 add_static_variable( v_public_power_injection , "public_power_injection" );

 // the power absorbed variables
 v_public_power_absorption.resize( number_nodes );
 for( auto & var : v_public_power_absorption )
  var.set_type( ColVariable::kNonNegative );
 add_static_variable( v_public_power_absorption , "public_power_absorption" );

 // the max power variables
 v_max_power.resize( number_nodes );
 for( auto & var : v_max_power )
  var.set_type( ColVariable::kNonNegative );
 add_static_variable( v_max_power , "max_power" );

 set_variables_generated();
}

/*--------------------------------------------------------------------------*/

void ECNetworkBlock::generate_abstract_constraints(
 Configuration * stcc ) {

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

 power_flow_limit_constraints.resize(
  boost::multi_array< FRowConstraint , 3 >::extent_gen()
  [ number_nodes ][ number_intervals ][ 2 ] ); // 2 dims, i.e., the sign (+/-)

 for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {

  for( Index t = 0 ; t < number_intervals ; ++t ) {

   LinearFunction::v_coeff_pair vars_p;
   LinearFunction::v_coeff_pair vars_n;

   // P^{max} vars also depends from P^{M+} and P^{M-}
   // vars as specified in the paper
   if( true ) { // default_config

    // case (1)
    vars_p.push_back( std::make_pair( &v_micro_power_injection[ node_id ] ,
                                      1.0 ) );
    vars_p.push_back( std::make_pair( &v_micro_power_absorption[ node_id ] ,
                                      -1.0 ) );

    // case (2)
    vars_n.push_back( std::make_pair( &v_micro_power_injection[ node_id ] ,
                                      -1.0 ) );
    vars_n.push_back( std::make_pair( &v_micro_power_absorption[ node_id ] ,
                                      1.0 ) );
   }

   // case (1)
   vars_p.push_back( std::make_pair( &v_public_power_injection[ node_id ] ,
                                     1.0 ) );
   vars_p.push_back( std::make_pair( &v_public_power_absorption[ node_id ] ,
                                     -1.0 ) );
   vars_p.push_back( std::make_pair( &v_max_power[ node_id ] , -1.0 ) );

   // case (2)
   vars_n.push_back( std::make_pair( &v_public_power_injection[ node_id ] ,
                                     -1.0 ) );
   vars_n.push_back( std::make_pair( &v_public_power_absorption[ node_id ] ,
                                     1.0 ) );
   vars_n.push_back( std::make_pair( &v_max_power[ node_id ] , -1.0 ) );

   // case (1)
   power_flow_limit_constraints[ node_id ][ t ][ 0 ].set_rhs( 0.0 );
   power_flow_limit_constraints[ node_id ][ t ][ 0 ].set_lhs(
    -Inf< double >() );
   power_flow_limit_constraints[ node_id ][ t ][ 0 ].set_function(
    new LinearFunction( std::move( vars_p ) ) );

   // case (2)
   power_flow_limit_constraints[ node_id ][ t ][ 1 ].set_rhs( 0.0 );
   power_flow_limit_constraints[ node_id ][ t ][ 1 ].set_lhs(
    -Inf< double >() );
   power_flow_limit_constraints[ node_id ][ t ][ 1 ].set_function(
    new LinearFunction( std::move( vars_n ) ) );
  }
 }

 add_static_constraint( power_flow_limit_constraints ,
                        "power_flow_limit_constraints" );

/*-------------------------- equality constraints --------------------------*/

 // set the power balance within the microgrid market/network, i.e., the
 // flows within the microgrid to have sum equal to zero:
 //
 //    P^{M,+} = P^{M,-}       for all t
 // => P^{M,+} - P^{M,-} = 0   for all t

 if( micro_power_balance_constraints.size() != number_intervals ) {
  assert( micro_power_balance_constraints.empty() );
  micro_power_balance_constraints.resize( number_intervals );
 }

 for( Index t = 0 ; t < number_intervals ; ++t ) {

  LinearFunction::v_coeff_pair vars;

  for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {

   vars.push_back(
    std::make_pair( &v_micro_power_injection[ node_id ] , 1.0 ) );
   vars.push_back(
    std::make_pair( &v_micro_power_absorption[ node_id ] , -1.0 ) );
  }

  micro_power_balance_constraints[ t ].set_both( 0.0 );
  micro_power_balance_constraints[ t ].set_function(
   new LinearFunction( std::move( vars ) ) );
 }

 add_static_constraint( micro_power_balance_constraints ,
                        "micro_power_balance_constraints" );

 // set the power balance, i.e.:
 //
 //    P^{POD,+} - P^{POD,-} - node_injection = - active_demand   for all t, u
 //
 // where:
 //
 //    P^{POD,+} = P^{P,+} + P^{M,+}
 //    P^{POD,-} = P^{P,-} + P^{M,-}

 power_balance_constraints.resize(
  boost::multi_array< FRowConstraint , 2 >::extent_gen()
  [ number_nodes ][ number_intervals ] );

 for( Index t = 0 ; t < number_intervals ; ++t ) {

  for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {

   LinearFunction::v_coeff_pair vars;

   vars.push_back( std::make_pair( &v_public_power_injection[ node_id ] ,
                                   1.0 ) );
   vars.push_back( std::make_pair( &v_micro_power_injection[ node_id ] ,
                                   1.0 ) );
   vars.push_back( std::make_pair( &v_public_power_absorption[ node_id ] ,
                                   -1.0 ) );
   vars.push_back( std::make_pair( &v_micro_power_absorption[ node_id ] ,
                                   -1.0 ) );
   vars.push_back( std::make_pair( &v_node_injection[ t ][ node_id ] , -1.0 ) );
   power_balance_constraints[ node_id ][ t ].set_both(
    -v_active_demand[ t ][ node_id ] );
   power_balance_constraints[ node_id ][ t ].set_function(
    new LinearFunction( std::move( vars ) ) );
  }
 }

 add_static_constraint( power_balance_constraints ,
                        "power_balance_constraints" );
}

void ECNetworkBlock::generate_objective( Configuration * objc ) {

 if( objective_generated() )
  return; // objective has already been generated

 if( get_objective() != nullptr )  // an objective is there already
  return;                          // cowardly (and silently) return

 LinearFunction::v_coeff_pair vars;

 for( Index node_id = 0 ; node_id < get_number_nodes() ; ++node_id ) {

  for( Index t = 0 ; t < get_number_intervals() ; ++t ) {

   // net economic balance wrt the public market
   vars.push_back( std::make_pair( &v_public_power_absorption[ node_id ] ,
                                   f_NetworkData->get_sell_price()[ t ] ) );
   vars.push_back( std::make_pair( &v_micro_power_absorption[ node_id ] ,
                                   f_NetworkData->get_sell_price()[ t ] ) );
   vars.push_back( std::make_pair( &v_public_power_injection[ node_id ] ,
                                   -f_NetworkData->get_buy_price()[ t ] ) );
   vars.push_back( std::make_pair( &v_micro_power_injection[ node_id ] ,
                                   -f_NetworkData->get_buy_price()[ t ] ) );
  }

  // the costs due to the peak power
  vars.push_back( std::make_pair( &v_max_power[ node_id ] ,
                                  f_NetworkData->get_max_tariff() ) );
 }

 auto lf = new LinearFunction( std::move( vars ) );
 lf->set_constant_term( f_const_term );
 objective.set_function( lf );
 objective.set_sense( Objective::eMax );

 // set block objective
 this->set_objective( &objective );
 set_objective_generated();
}

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

void ECNetworkBlock::set_active_demand(
 // TODO this should be a const ptr to double (?)
 std::vector< double >::const_iterator values ,
 Block::Subset && subset ,
 const bool ordered ,
 c_ModParam issuePMod ,
 c_ModParam issueAMod ) {
 // TODO
}

/*--------------------------------------------------------------------------*/

void ECNetworkBlock::set_active_demand(
 // TODO this should be a const ptr to double (?)
 std::vector< double >::const_iterator values ,
 Block::Range rng ,
 c_ModParam issuePMod ,
 c_ModParam issueAMod ) {
 // TODO
}

/*--------------------------------------------------------------------------*/
/*----------------------- End File ECNetworkBlock.cpp ----------------------*/
/*--------------------------------------------------------------------------*/