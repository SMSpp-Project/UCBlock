/*--------------------------------------------------------------------------*/
/*--------------------- File DCNetworkBlock.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the DCNetworkBlock class.
 *
 * \author Wim van Ackooij \n
 *         EDF R&D OSIRIS \n
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Ali Ghezelsoflu \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Quentin Jacquet \n
 *         EDF R&D OSIRIS \n
 *
 * \author Rafael Durbano Lobato \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Rafael Durbano Lobato
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
SMSpp_insert_in_factory_cpp_0( DCNetworkBlock );

// register DCNetworkBlockSolution to the Solution factory
SMSpp_insert_in_factory_cpp_0( DCNetworkBlockSolution );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

// register DCNetworkData to the NetworkData factory

typedef DCNetworkBlock::DCNetworkData DCNetworkData;

SMSpp_insert_in_factory_cpp_0( DCNetworkData );

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS OF DCNetworkData -------------------------*/
/*--------------------------------------------------------------------------*/

void DCNetworkData::deserialize( const netCDF::NcGroup & group )
{
 #ifndef NDEBUG
  static std::vector< std::string > expected_dims = { "NumberNodes" ,
   "NumberLines" ,  // if called from UCBlock:
   "TimeHorizon" , "NumberUnits" , "NumberNetworks" ,
   "NumberElectricalGenerators" };
  check_dimensions( group , expected_dims , std::cerr );

  static std::vector< std::string > expected_vars = { "ActiveDemand" ,
   "StartLine" , "EndLine" , "MinPowerFlow" ,"MaxPowerFlow" ,
   "LineSusceptance" , "NetworkCost" , "NodeName" , "LineName" ,
   "ConstantTerm" , "Efficiency", // if called from UCBlock:
   "ActivePowerDemand" , "GeneratorNode" , "NetworkConstantTerms" ,
   "NetworkBlockClassname" , "NetworkDataClassname",
   // vars for AC Mode MUST NOT BE HERE!
   "ReactivePowerDemand" , "NodeConductance" , "NodeSusceptance" ,
   "NodeVoltageMagnitude" , "NodeVoltageAngle", "NodeMaxVoltage" ,
   "NodeMinVoltage" , "LineResistance" , "LineReactance" , "LineRatio",
   "LineRATEA" , "LineShiftAngle" , "LineMinAngle", "LineMaxAngle" };
  check_variables( group , expected_vars , std::cerr );
 #endif

 NetworkData::deserialize( group );

 // Optional variables

 if( f_number_nodes > 1 ) {

  deserialize_dim( group , "NumberLines" , f_number_lines , false );

  ::deserialize( group , "StartLine" , f_number_lines , v_start_line , false ,
                 true );

  for( Index i = 0 ; i < f_number_lines ; ++i ) {
   if( ( v_start_line[ i ] < 0 ) || ( v_start_line[ i ] >= f_number_nodes ) )
    throw( std::invalid_argument( "DCNetworkData::deserialize: "
                                  "wrong start node number " +
                                  std::to_string( v_start_line[ i ] ) ) );

  ::deserialize( group , "EndLine" , f_number_lines , v_end_line , false ,
                 true );

   if( ( v_end_line[ i ] < 0 ) || ( v_end_line[ i ] >= f_number_nodes ) )
    throw( std::invalid_argument( "DCNetworkData::deserialize: "
                                  "wrong end node number " +
                                  std::to_string( v_end_line[ i ] ) ) );

   if( v_start_line[ i ] == v_end_line[ i ] )
    throw( std::invalid_argument( "DCNetworkData::deserialize: "
                                  "start node == end node for line " +
                                  std::to_string( v_end_line[ i ] ) ) );
   }

  ::deserialize( group , "MinPowerFlow" , f_number_lines , v_min_power_flow ,
                 true , true );

  ::deserialize( group , "MaxPowerFlow" , f_number_lines , v_max_power_flow ,
                 true , true );

  ::deserialize( group , "NetworkCost" , f_number_lines , v_network_cost ,
                 true , true );

  ::deserialize( group , "Efficiency" , f_number_lines , v_efficiency ,
                true , true );

  if( ! deserialize_dim( group, "ReferenceNode", f_reference_node, true ) )
    f_reference_node = 0;
  }

  // AC vars MUST NOT NE HERE
  ::deserialize( group , "LineSusceptance" , f_number_lines ,
		 v_line_susceptance , true , true );

  ::deserialize( group , "LineReactance" , f_number_lines ,
		 v_line_reactance , true , true );

  ::deserialize( group , "LineResistance" , f_number_lines ,
		 v_line_resistance , true , true );

  ::deserialize( group , "LineRatio" , f_number_lines , v_line_ratio ,
                 true , true );

  ::deserialize( group , "LineRATEA" , f_number_lines , v_line_rate_A ,
                 true , true );

  ::deserialize( group , "LineShiftAngle" , f_number_lines , v_line_angle ,
                 true , true );

  ::deserialize( group , "LineMinAngle" , f_number_lines , v_line_min_angle ,
                 true , true );

  ::deserialize( group , "LineMaxAngle" , f_number_lines , v_line_max_angle ,
                 true , true );

  ::deserialize( group , "NodeConductance" , f_number_nodes ,
		 v_node_conductance , true , true );

  ::deserialize( group , "NodeSusceptance" , f_number_nodes ,
		 v_node_susceptance , true , true );

  ::deserialize( group , "NodeMaxVoltage" , f_number_nodes ,
		 v_node_max_voltage , true , true );

  ::deserialize( group , "NodeMinVoltage" , f_number_nodes ,
		 v_node_min_voltage , true , true );

 ::deserialize( group , "NodeName" , f_number_nodes , v_node_names );
 ::deserialize( group , "LineName" , f_number_lines , v_line_names );

 f_lines_type = -1;

 }  // end( DCNetworkData::deserialize )

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

void DCNetworkBlock::deserialize( const netCDF::NcGroup & group )
{

#ifndef NDEBUG
 static std::vector< std::string > expected_dims = { "NumberNodes" ,
                                                     "NumberLines" };
 check_dimensions( group , expected_dims , std::cerr );

 static std::vector< std::string > expected_vars = { "ActiveDemand" ,
                                                     "StartLine" ,
                                                     "EndLine" ,
                                                     "MinPowerFlow" ,
                                                     "MaxPowerFlow" ,
                                                     "LineSusceptance" ,
                                                     "NodeSusceptance" ,
                                                     "NetworkCost" ,
                                                     "NodeName" ,
                                                     "LineName" ,
                                                     "ConstantTerm",
                                                     "Efficiency" };
 check_variables( group , expected_vars , std::cerr );
#endif

 NetworkBlock::deserialize( group );

 // Optional variables

 Index NumberNodes;
 if( deserialize_dim( group , "NumberNodes" , NumberNodes ) ) {
  // Since the dimension "NumberNodes" has been provided, it means that a
  // DCNetworkData has been provided. Thus, the DCNetworkData is deserialized,
  // and it is marked as being local
  if( f_local_NetworkData )
   // if the NetworkData has not been passed from UCBlock, then delete it
   delete( f_NetworkData );
  auto DCND = new DCNetworkData();
  DCND->deserialize( group );
  if( f_NetworkData &&
    ( f_NetworkData->get_number_nodes() != DCND->get_number_nodes() ) )
   throw( std::logic_error(
    "DCNetworkBlock::deserialize: NumberNodes not matching between NetworkData" ) );
  f_NetworkData = DCND;
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
}  // end( DCNetworkBlock::deserialize )

/*--------------------------------------------------------------------------*/

Eigen::MatrixXd DCNetworkBlock::get_PTDF( const std::vector< Index > & AC_lines ) {
  // get data
  const auto & susceptance = f_NetworkData->get_line_susceptance();
  const auto number_nodes = get_number_nodes();
  const auto number_lines = get_number_lines();
  if( number_lines <= 0 ) {
    throw ( std::logic_error( "DCNetworkBlock::generate_abstract_constraints: "
                              "number of lines of DCNetworkBlock is not set" ) );
  }
  const auto & start_line = f_NetworkData->get_start_line();
  const auto & end_line = f_NetworkData->get_end_line();
  
  // construct the matrix using two sub-matrices B_bar and B_hat
  Eigen::MatrixXd B_hat = Eigen::MatrixXd::Zero(number_lines,number_nodes );
  for( auto & line_id : AC_lines ) {
    B_hat(line_id,start_line[line_id]) =  susceptance[line_id];
    B_hat(line_id,end_line[line_id])   = -susceptance[line_id];
  }
  
  Eigen::MatrixXd B_bar = Eigen::MatrixXd::Zero(number_nodes,number_nodes );
  for( Index node_id = 0; node_id < number_nodes; ++node_id ) {
    for( auto & line_id : AC_lines ) {
      if( start_line[ line_id ] == node_id ) {
        B_bar(node_id,end_line[line_id]) = - susceptance[line_id];
        B_bar(node_id,node_id ) += susceptance[line_id];
      }
      if( end_line[ line_id ] == node_id ) {
        B_bar(node_id,start_line[line_id]) = - susceptance[line_id];
        B_bar(node_id,node_id ) += susceptance[line_id];
      }
    }
  }

  Index ref_node = f_NetworkData->get_reference_node();
  Eigen::MatrixXd I_nref = Eigen::MatrixXd::Zero(number_nodes,number_nodes-1);
  for( Index node_id = 0; node_id < ref_node; ++node_id ) {
    I_nref( node_id , node_id ) = 1.;
  }
  for( Index node_id = ref_node; node_id < number_nodes - 1; ++node_id ) {
    I_nref( node_id + 1 , node_id ) = 1.;
  }

  Eigen::MatrixXd B1 = B_hat*I_nref;
  Eigen::MatrixXd B2 = I_nref.transpose()*B_bar*I_nref;
  Eigen::MatrixXd PTDF_matrix  = B1*B2.inverse();

  return( PTDF_matrix );
}

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::generate_abstract_variables( Configuration * stvv )
{
 if( variables_generated() )  // variables have already been generated
  return;                     // nothing to do

 NetworkBlock::generate_abstract_variables( stvv );

 const auto number_lines = get_number_lines();

 if( number_lines > 0 ) {
  // the power flow Variable
  v_power_flow.resize( number_lines ); // TODO: only needed for DC lines
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

int DCNetworkBlock::get_reducedIdx( int idx ) {
 if( idx > f_NetworkData->get_reference_node() )
  return( idx - 1 );
 return( idx );
 }

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::generate_abstract_constraints( Configuration * stcc )
{
 if( constraints_generated() )  // constraints have already been generated
  return;                       // nothing to do

 NetworkBlock::generate_abstract_constraints( stcc );

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

 // Splitting AC and DC part
 std::vector< Index > AC_lines = get_AC_lines();
 Eigen::MatrixXd PTDF_matrix = get_PTDF( AC_lines );
 std::vector< Index > DC_lines = get_DC_lines();

 // ===== auxiliary variables for nonempty cost
 if( ! f_NetworkData->get_network_cost().empty() ) {
   // 0 <= V_l - F_l && 0 <= V_l + F_l
   v_power_flow_relax_abs.resize(
    boost::multi_array< FRowConstraint , 2 >::extent_gen()[ 2 ][ number_lines ] );

   // DC part 
   for( auto & line_id : DC_lines ) {
    auto lfunc_1 = new LinearFunction();
    lfunc_1->add_variable( &v_power_flow[ line_id ] , -1.0 );
    lfunc_1->add_variable( &v_auxiliary_variable[ line_id ] , 1.0 );
    v_power_flow_relax_abs[ 0 ][ line_id ].set_lhs( 0.0 );
    v_power_flow_relax_abs[ 0 ][ line_id ].set_rhs( Inf< double >() );
    v_power_flow_relax_abs[ 0 ][ line_id ].set_function( lfunc_1 );

    auto lfunc_2 = new LinearFunction();
    lfunc_2->add_variable( &v_power_flow[ line_id ] , 1.0 );
    lfunc_2->add_variable( &v_auxiliary_variable[ line_id ] , 1.0 );
    v_power_flow_relax_abs[ 1 ][ line_id ].set_lhs( 0.0 );
    v_power_flow_relax_abs[ 1 ][ line_id ].set_rhs( Inf< double >() );
    v_power_flow_relax_abs[ 1 ][ line_id ].set_function( lfunc_2 );
   }

   // AC part
   for( auto & line_id : AC_lines ) {
    auto lfunc_1 = new LinearFunction();
    auto lfunc_2 = new LinearFunction();
    double constant_term = 0;
    for( Index node_id = 0; node_id < number_nodes; ++node_id ) {
      if( node_id != f_NetworkData->get_reference_node() ) {
       double coefficient = PTDF_matrix( line_id , get_reducedIdx( node_id ) );
       lfunc_1->add_variable( &v_node_injection[ 0 ][ node_id ] , -coefficient );
       lfunc_2->add_variable( &v_node_injection[ 0 ][ node_id ] , coefficient );
       constant_term -= coefficient * v_ActiveDemand[ node_id ];
     }
    } // for each node
    lfunc_1->add_variable( &v_auxiliary_variable[ line_id ] , 1.0);
    v_power_flow_relax_abs[ 0 ][ line_id ].set_lhs( constant_term );
    v_power_flow_relax_abs[ 0 ][ line_id ].set_rhs( Inf< double >() );
    v_power_flow_relax_abs[ 0 ][ line_id ].set_function( lfunc_1 );
    lfunc_2->add_variable( &v_auxiliary_variable[ line_id ] , 1.0);
    v_power_flow_relax_abs[ 1 ][ line_id ].set_lhs( -constant_term );
    v_power_flow_relax_abs[ 1 ][ line_id ].set_rhs( Inf< double >() );
    v_power_flow_relax_abs[ 1 ][ line_id ].set_function( lfunc_2 );
   }
   add_static_constraint( v_power_flow_relax_abs , "power_flow_relax_abs" );
  } // ===== end( cost not empty )


 // ===== constraints on the DC part
 if( lines_type == kHVDC || lines_type == kAC_HVDC ) {

  // Flow limit constraints
  if( lines_type == kHVDC )
   v_HVDC_power_flow_limit_const.resize( number_lines );
  if( lines_type == kAC_HVDC )
   v_AC_HVDC_power_flow_limit_const.resize( number_lines );

  for( auto & line_id : DC_lines ) {
    const auto kappa = get_kappa( line_id );
    if( lines_type == kHVDC ) {
      v_HVDC_power_flow_limit_const[ line_id ].set_lhs( kappa * get_min_power_flow( line_id ) );
      v_HVDC_power_flow_limit_const[ line_id ].set_rhs( kappa * get_max_power_flow( line_id ) );
      v_HVDC_power_flow_limit_const[ line_id ].set_variable( & v_power_flow[ line_id ] );
      add_static_constraint( v_HVDC_power_flow_limit_const, "HVDC_power_flow_limit" );
    }
    else {
     v_AC_HVDC_power_flow_limit_const[ line_id ].set_lhs( kappa * get_min_power_flow( line_id ) );
     v_AC_HVDC_power_flow_limit_const[ line_id ].set_rhs( kappa * get_max_power_flow( line_id ) );
     auto lfunc = new LinearFunction();
     lfunc->add_variable( & v_power_flow[ line_id ], 1. );
     v_AC_HVDC_power_flow_limit_const[ line_id ].set_function( lfunc );
     // power_flow_limit constraints will be completed in AC part
    }
  }

  // Power flow and node injection constraints
  if( lines_type == kHVDC )
   v_power_flow_injection_const.resize(number_nodes );
  //if( lines_type == kAC_HVDC ) v_AC_HVDC_power_flow_constraints.resize( number_nodes );

  for( Index n = 0 ; n < number_nodes ; ++n ) {
    auto lfunc = new LinearFunction();
    lfunc->add_variable( & v_node_injection[ 0 ][ n ] , -1.0 );

    double eta = 1.0;
    for( auto & line_id : DC_lines ) {
      eta = f_NetworkData->get_line_efficiency( line_id ); // efficiency of the HVDC line
      if( start_line[ line_id ] == n )  lfunc->add_variable( & v_power_flow[ line_id ] , 1.0 );
      if( end_line[ line_id ] == n )    lfunc->add_variable( & v_power_flow[ line_id ] , -eta );
    }
    if( lines_type == kHVDC ) {
      v_power_flow_injection_const[ n ].set_both( -v_ActiveDemand[ n ] );
      v_power_flow_injection_const[ n ].set_function( lfunc );
    }
    /*else {
      v_AC_HVDC_power_flow_const[ n ].set_both( -v_ActiveDemand[ n ] );
      v_AC_HVDC_power_flow_const[ n ].set_function( lfunc );
    }*/

  }
  if( lines_type == kHVDC )
    add_static_constraint( v_power_flow_injection_const, "HVDC_power_flow_injection" );
  /*else
    add_static_constraint(v_AC_HVDC_power_flow_constraints, "AC/HVDC_power_flow_injection");*/

 } // ===== end constraints on HVDC part

 // ===== constraints on AC Part
 if( lines_type == kAC || lines_type == kAC_HVDC ) {
  Eigen::MatrixXd linkingMat;
  if( lines_type == kAC_HVDC ) {
    // linking constraints between AC and HVDC
    Eigen::MatrixXd A_DC = Eigen::MatrixXd::Zero(number_nodes-1, number_lines );
    double eta = 1.0;
    for( auto & line_id : DC_lines ) {
      eta = f_NetworkData->get_line_efficiency( line_id ); // efficiency of the HVDC line
      A_DC(get_reducedIdx( start_line[ line_id ] ) , line_id ) = 1.;
      A_DC(get_reducedIdx( end_line[ line_id ] ) , line_id )   = - eta; // QJ_TOCHECK 1 or -1 ?
    }
    linkingMat = - PTDF_matrix * A_DC;
  }

  // Flow limit constraints
  if( lines_type == kAC )
   v_AC_power_flow_limit_const.resize(number_lines );
  for( auto & line_id : AC_lines ) {
    const auto kappa = get_kappa( line_id );
    auto lfunc = new LinearFunction();
    double constant_term = 0;

    for( Index node_id = 0; node_id < number_nodes; ++node_id ) {
      if( node_id != f_NetworkData->get_reference_node() ) {
       double coefficient = PTDF_matrix(line_id, get_reducedIdx( node_id ) ); // Distribution Factor Matrix
       lfunc->add_variable( &v_node_injection[ 0 ][ node_id ], coefficient );
       constant_term -= coefficient * v_ActiveDemand[ node_id ];
     }
    } // for each node

    if( lines_type == kAC_HVDC ) {
      for( auto & dc_line_id: DC_lines ) {
        lfunc->add_variable( & v_power_flow[ dc_line_id ], linkingMat( line_id , dc_line_id ) );
      } // for each dc line

      // Set the constraint (AC-HVDC)
      v_AC_HVDC_power_flow_limit_const[line_id].set_function( lfunc );
      v_AC_HVDC_power_flow_limit_const[line_id].set_lhs( kappa*get_min_power_flow( line_id ) - constant_term );
      v_AC_HVDC_power_flow_limit_const[line_id].set_rhs( kappa*get_max_power_flow( line_id ) - constant_term );
    }
    else {
      // Set the constraint (AC)
      v_AC_power_flow_limit_const[line_id].set_function( lfunc );
      v_AC_power_flow_limit_const[line_id].set_lhs( kappa*get_min_power_flow( line_id ) - constant_term );
      v_AC_power_flow_limit_const[line_id].set_rhs( kappa*get_max_power_flow( line_id ) - constant_term );
    }
  }

  if( lines_type == kAC_HVDC )
    add_static_constraint( v_AC_HVDC_power_flow_limit_const, "AC/HVDC_power_flow_limits" );
  else
    add_static_constraint( v_AC_power_flow_limit_const, "AC_power_low_limits" );

  auto lfunc = new LinearFunction();
  double constant_term = 0.;
  for( Index node_id = 0; node_id < number_nodes; ++node_id ) {
    lfunc->add_variable( &v_node_injection[ 0 ][ node_id ], 1. );
    constant_term += v_ActiveDemand[ node_id ];
  }
  overall_balanced_const.set_function( lfunc );
  overall_balanced_const.set_lhs( constant_term );
  overall_balanced_const.set_rhs( constant_term );
  add_static_constraint( overall_balanced_const, "overall_balanced_const");

 } // ===== end AC and AC-HVDC constraints

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
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/

Solution * DCNetworkBlock::get_Solution( Configuration * csolc ,
					 bool emptys )
{
 Index wsol = 7;
 if( ( ! csolc ) && f_BlockConfig )
  csolc = f_BlockConfig->f_solution_Configuration;

 if( auto config = dynamic_cast< SimpleConfiguration< int > * >( csolc ) )
  wsol = config->f_value;

 // call the method of the base class
 auto * sol = dynamic_cast< DCNetworkBlockSolution * >(
		               NetworkBlock::get_Solution( csolc , emptys ) );
 assert( sol );

 if( wsol & 2 )
  sol->v_flow.resize( get_number_lines() );

 if( wsol & 4 )
  sol->v_cost.resize( get_number_lines() );

 if( ! emptys )
  sol->read( this );

 return( sol );
 }

/*--------------------------------------------------------------------------*/
 
NetworkBlockSolution * DCNetworkBlock::new_Solution( void ) const {
 return( new DCNetworkBlockSolution() );
 }

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
  && RowConstraint::is_feasible( v_AC_power_flow_limit_const , tol ,
				 rel_viol )
  && RowConstraint::is_feasible( v_AC_HVDC_power_flow_limit_const , tol ,
				 rel_viol )
  && RowConstraint::is_feasible( v_power_flow_injection_const , tol ,
				 rel_viol )
  && RowConstraint::is_feasible( v_power_flow_relax_abs , tol , rel_viol )
  && RowConstraint::is_feasible( v_HVDC_power_flow_limit_const , tol ,
				 rel_viol )
  && RowConstraint::is_feasible( node_injection_bounds_const , tol ,
				 rel_viol ) );

} // end( DCNetworkBlock::is_feasible )

/*--------------------------------------------------------------------------*/
/*--------- METHODS FOR LOADING, PRINTING & SAVING THE DCNetworkBlock ------*/
/*--------------------------------------------------------------------------*/

void DCNetworkData::serialize( netCDF::NcGroup & group ) const
{
 NetworkData::serialize( group );

 if( f_number_nodes > 1 ) {
  auto NumberNodes = group.getDim( "NumberNodes" );
  auto NumberLines = group.addDim( "NumberLines" );

  ::serialize( group , "StartLine" , netCDF::NcUint() , NumberLines ,
               v_start_line );

  ::serialize( group , "EndLine" , netCDF::NcUint() , NumberLines ,
               v_end_line );

  ::serialize( group , "MinPowerFlow" , netCDF::NcDouble() , NumberLines ,
               v_min_power_flow );

  ::serialize( group , "MaxPowerFlow" , netCDF::NcDouble() , NumberLines ,
               v_max_power_flow );

  ::serialize( group , "LineSusceptance" , netCDF::NcDouble() , NumberLines ,
               v_line_susceptance );

  ::serialize( group , "NodeSusceptance" , netCDF::NcDouble() , NumberNodes ,
               v_node_susceptance );

  ::serialize( group , "NetworkCost" , netCDF::NcDouble() , NumberLines ,
               v_network_cost );

  ::serialize( group , "Efficiency" , netCDF::NcDouble() , NumberLines ,
              v_efficiency );

  if( ! v_line_names.empty() ) {
   assert( v_line_names.size() == NumberLines.getSize() );
   auto LineName = group.addVar( "LineName" , netCDF::NcString() ,
				 NumberLines );
   for( Index i = 0 ; i < v_line_names.size() ; ++i )
    LineName.putVar( { i } , v_line_names[ i ] );
   }
  }

 if( ! v_node_names.empty() ) {
  auto NumberNodes = group.getDim( "NumberNodes" );
  assert( v_node_names.size() == NumberNodes.getSize() );
  auto NodeName = group.addVar( "NodeName" , netCDF::NcString() ,
				NumberNodes );
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
 }  // end( DCNetworkBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

void DCNetworkBlock::set_active_demand( MF_dbl_it values , Subset && subset ,
					bool ordered ,
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
      if( f_NetworkData->get_lines_type() == kHVDC ) {
          std::vector< Index > modified_nodes(subset.begin(), subset.end());
          change_DC_power_flow_injection_constraints( modified_nodes , issueAMod );
      }
      if( f_NetworkData->get_lines_type() == kAC || f_NetworkData->get_lines_type() == kAC_HVDC ) {
          std::vector< Index > modified_lines( f_NetworkData->get_number_lines());
          std::iota( modified_lines.begin(), modified_lines.end() , 0 );
          change_relax_abs_constraints(modified_lines, issueAMod ); // all lines are modified
          change_power_flow_limit_constraints(modified_lines, issueAMod );
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

void DCNetworkBlock::set_active_demand( MF_dbl_it values , Range rng ,
                                        c_ModParam issuePMod ,
                                        c_ModParam issueAMod )
{
 rng.second = std::min( rng.second , get_number_nodes() );
 if( rng.second <= rng.first )
  return;

 if( v_ActiveDemand.empty() ) {
  if( std::all_of( values , values + ( rng.second - rng.first ) ,
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
  std::copy( values , values + ( rng.second - rng.first ) ,
             v_ActiveDemand.begin() + rng.first );

  if( not_dry_run( issueAMod ) && constraints_generated() ) {
   // Change the physical representation

  std::copy( values , values + ( rng.second - rng.first ) ,
             v_ActiveDemand.begin() + rng.first );

  if( not_dry_run( issueAMod ) && constraints_generated() ) {
       // Change the abstract representation

       if( f_NetworkData->get_lines_type() == kHVDC ) {
          std::vector< Index > modified_nodes( rng.second - rng.first );
          std::iota( modified_nodes.begin() , modified_nodes.end() ,
		     rng.first );
          change_DC_power_flow_injection_constraints( modified_nodes ,
						      issueAMod );
        }
       if( f_NetworkData->get_lines_type() == kAC ||
	   f_NetworkData->get_lines_type() == kAC_HVDC ) {
	std::vector< Index > modified_lines( f_NetworkData->get_number_lines());
	std::iota( modified_lines.begin() , modified_lines.end(),0);
	change_relax_abs_constraints( modified_lines , issueAMod ); // all lines are modified
	change_power_flow_limit_constraints( modified_lines , issueAMod );
        }
    }
  }
 }

 if( issue_pmod( issuePMod ) )  // issue a Physical Modification
  Block::add_Modification( std::make_shared< NetworkBlockRngdMod >( this ,
					 NetworkBlockMod::eSetActD , rng ) ,
                           Observer::par2chnl( issuePMod ) );

 }  // end( DCNetworkData::set_active_demand )

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::set_kappa( MF_dbl_it values , Subset && subset ,
                                bool ordered ,
                                c_ModParam issuePMod , c_ModParam issueAMod )
{
 if( subset.empty() )
  return;

 if( v_kappa.empty() ) {
  if( std::all_of( values , values + subset.size() ,
                   []( double cst ) { return( cst == 1 ); } ) )
   return;

  v_kappa.resize( get_number_lines() , 1 );
  }

 bool identical = true;
 for( auto i : subset ) {
  if( i >= v_kappa.size() )
   throw( std::invalid_argument( "DCNetworkBlock::set_kappa: invalid value in"
                                 " subset: " + std::to_string( i ) ) );

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
  std::vector< Index > modified_lines( subset.begin() , subset.end() );
  change_power_flow_limit_constraints( modified_lines , issueAMod );
  // kappa only appears in limit constraints
  }

 if( issue_pmod( issuePMod ) ) {  // issue a Physical Modification
  if( ! ordered )
   std::sort( subset.begin() , subset.end() );

  Block::add_Modification( std::make_shared< DCNetworkBlockSbstMod >( this ,
		       DCNetworkBlockMod::eSetKappa , std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );
  }
 }  // end( DCNetworkData::set_kappa( subset ) )

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::set_kappa( MF_dbl_it values , Range rng ,
                                c_ModParam issuePMod , c_ModParam issueAMod )
{
 rng.second = std::min( rng.second , get_number_lines() );
 if( rng.second <= rng.first )
  return;

 if( v_kappa.empty() ) {
  if( std::all_of( values , values + ( rng.second - rng.first ) ,
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
  std::copy( values , values + ( rng.second - rng.first ) ,
             v_kappa.begin() + rng.first );

  if( not_dry_run( issueAMod ) && constraints_generated() ) {
   // Change the abstract representation
   std::vector< Index > modified_lines( rng.second - rng.first );
   std::iota( modified_lines.begin() , modified_lines.end() , rng.first );
   change_power_flow_limit_constraints( modified_lines , issueAMod );
   }
  }

 if( issue_pmod( issuePMod ) )  // issue a Physical Modification
  Block::add_Modification( std::make_shared< DCNetworkBlockRngdMod >( this ,
				       DCNetworkBlockMod::eSetKappa , rng ) ,
                           Observer::par2chnl( issuePMod ) );

 }  // end( DCNetworkData::set_kappa( range ) )

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::change_power_flow_limit_constraints(
         const std::vector< Index > & modified_lines , c_ModParam issueAMod )
{
  Eigen::MatrixXd PTDF_matrix;
  switch( f_NetworkData->get_lines_type() ) {
   case( kHVDC ): {
    for( auto i : modified_lines ) {
     v_HVDC_power_flow_limit_const[ i ].set_lhs(
		       v_kappa[ i ] * get_min_power_flow( i ) , issueAMod );
     v_HVDC_power_flow_limit_const[ i ].set_rhs(
		       v_kappa[ i ] * get_max_power_flow( i ) , issueAMod );
    }
    break;
    }

   case( kAC ): {
    PTDF_matrix = get_PTDF();
    for( auto i : modified_lines ) {
     double constant_term = 0;
     for( Index node_id = 0 ; node_id < f_NetworkData->get_number_nodes() ;
	  ++node_id )
      constant_term -= PTDF_matrix( i , get_reducedIdx( node_id ) )
                     * v_ActiveDemand[ node_id ];

     v_AC_power_flow_limit_const[ i ].set_lhs(
         v_kappa[ i ] * get_min_power_flow( i ) - constant_term , issueAMod );

     v_AC_power_flow_limit_const[ i ].set_rhs(
	 v_kappa[ i ] * get_max_power_flow( i ) - constant_term , issueAMod );
     }
    break;
    }

   case( kAC_HVDC ): {
    std::vector< Index > AC_lines = get_AC_lines();
    PTDF_matrix = get_PTDF( AC_lines );
    for( auto i : modified_lines ) {
     double constant_term = 0;
     for( Index node_id = 0 ; node_id < f_NetworkData->get_number_nodes() ;
	  ++node_id )
      constant_term -= PTDF_matrix( i , get_reducedIdx( node_id ) )
                     * v_ActiveDemand[ node_id ];

     v_AC_HVDC_power_flow_limit_const[ i ].set_lhs(
	v_kappa[ i ] * get_min_power_flow( i ) - constant_term , issueAMod );

     v_AC_HVDC_power_flow_limit_const[ i ].set_rhs(
        v_kappa[ i ] * get_max_power_flow( i ) - constant_term , issueAMod );
     }
    }
   default: break;
   }
 }

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::change_relax_abs_constraints(
         const std::vector< Index > & modified_lines , c_ModParam issueAMod )
{
 if( ( f_NetworkData->get_lines_type() != kAC ) &&
     ( f_NetworkData->get_lines_type() != kAC_HVDC ) )
  return;

 std::vector< Index > AC_lines = get_AC_lines();
 Eigen::MatrixXd PTDF_matrix = get_PTDF( AC_lines );
 for( auto & i : modified_lines ) {
  double constant_term = 0;
  for( Index node_id = 0 ; node_id < f_NetworkData->get_number_nodes() ;
       ++node_id )
   constant_term -= PTDF_matrix( i , get_reducedIdx( node_id ) )
                  * v_ActiveDemand[ node_id ];

  v_power_flow_relax_abs[ 0 ][ i ].set_lhs( constant_term );
  v_power_flow_relax_abs[ 1 ][ i ].set_lhs( -constant_term );
  }
 }

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::change_DC_power_flow_injection_constraints(
         const std::vector< Index > & modified_nodes , c_ModParam issueAMod )
{
 if( f_NetworkData->get_lines_type() != kHVDC )
  return;

 for( auto & n : modified_nodes )
  v_power_flow_injection_const[ n ].set_both( -v_ActiveDemand[ n ] );
 }

/*--------------------------------------------------------------------------*/
/*------------------ METHODS OF DCNetworkBlockSolution ---------------------*/
/*--------------------------------------------------------------------------*/

void DCNetworkBlockSolution::deserialize( const netCDF::NcGroup & group )
{
 // call the method of the base class
 NetworkBlockSolution::deserialize( group );

 // "NumberLines" is mandatory- - - - - - - - - - - - - - - - - - - - - - - -
 deserialize_dim( group , "NumberLines" , f_number_lines , false );

 // deserialize the Flow Variables - - - - - - - - - - - - - - - - - - - - -
 ::deserialize< double >( group , "FlowValue" , v_flow , false );

 // deserialize the Dual Prices- - - - - - - - - - - - - - - - - - - - - - -
 ::deserialize< double >( group , "DualCost" , v_cost , false );

 }  // end( DCNetworkBlockSolution::deserialize( NcGroup & ) )

/*--------------------------------------------------------------------------*/

void DCNetworkBlockSolution::deserialize( const netCDF::NcGroup & group ,
					  int start , int count )
{
 // call the method of the base class
 NetworkBlockSolution::deserialize( group , start , count );

 // "NumberLines" is mandatory- - - - - - - - - - - - - - - - - - - - - - - -
 deserialize_dim( group , "NumberLines" , f_number_lines , false );

 std::vector< size_t > strt = { start , 0 };
 std::vector< size_t > cnt = { count , f_number_lines };

 // deserialize the Flow Variables - - - - - - - - - - - - - - - - - - - - -
 auto ncVar = group.getVar( "FlowValue" );
 if( ncVar.isNull() ) {
  std::vector< boost::multi_array< double , 2 >::index > sizes( 2 , 0 );
  v_flow.resize( sizes );
  }
 else
  ncVar.getVar( strt , cnt , v_flow.data() );

 // deserialize the Dual Prices- - - - - - - - - - - - - - - - - - - - - - -
 ncVar = group.getVar( "DualCost" );
 if( ncVar.isNull() ) {
  std::vector< boost::multi_array< double , 2 >::index > sizes( 2 , 0 );
  v_cost.resize( sizes );
  return;
  }

 ncVar.getVar( strt , cnt , v_cost.data() );

 }  // end( DCNetworkBlockSolution::deserialize( NcGroup & , int , int ) )

/*--------------------------------------------------------------------------*/

void DCNetworkBlockSolution::read( const Block * block )
{
 // call the method of the base class
 NetworkBlockSolution::read( block );

 auto DCNB = dynamic_cast< const DCNetworkBlock * >( block );
 if( ! DCNB )
  throw( std::invalid_argument(
	  "DCNetworkBlockSolution::read: block is not a DCNetworkBlock" ) );

 f_number_lines = DCNB->get_number_lines();

 if( ! v_flow.empty() ) {
  // read the flow power variables - - - - - - - - - - - - - - - - - - - - -
  auto Fl = DCNB->get_power_flow();
  for( Index l = 0 ; l < f_number_lines ; ++l )
   v_flow[ l ] = Fl[ l ].get_value();
  }

 if( ! v_cost.empty() )
  // read the dual prices- - - - - - - - - - - - - - - - - - - - - - - - - -
  DCNB->get_dual_prices( v_cost );

 }  // end( DCNetworkBlockSolution::read )

/*--------------------------------------------------------------------------*/

void DCNetworkBlockSolution::write( Block * block )
{
 // call the method of the base class
 NetworkBlockSolution::write( block );

 auto DCNB = dynamic_cast< DCNetworkBlock * >( block );
 if( ! DCNB )
  throw( std::invalid_argument(
	  "DCNetworkBlockSolution::write: block is not a DCNetworkBlock" ) );

 if( f_number_lines != DCNB->get_number_lines() )
  throw( std::invalid_argument(
	      "DCNetworkBlockSolution::write: inconsistent lines number" ) );

 if( ! v_flow.empty() ) {
  // write the flow power variables- - - - - - - - - - - - - - - - - - - - -
  auto Fl = DCNB->get_power_flow();
  for( Index l = 0 ; l < f_number_lines ; ++l )
   Fl[ l ].set_value( v_flow[ l ] );
  }

 if( ! v_cost.empty() )
  // write the dual prices - - - - - - - - - - - - - - - - - - - - - - - - -
  DCNB->set_dual_prices( v_cost );

 }  // end( DCNetworkBlockSolution::write )

/*--------------------------------------------------------------------------*/

void DCNetworkBlockSolution::serialize( netCDF::NcGroup & group ) const
{
 // call the method of the base class
 NetworkBlockSolution::serialize( group );

 // "NumberLines" is mandatory- - - - - - - - - - - - - - - - - - - - - - - -
 auto nl = group.addDim( "NumberLines" , f_number_lines );

 // serialize the Flow Variables- - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_flow.empty() )
  ::serialize< double >( group , "FlowValue" , netCDF::NcDouble() , nl ,
			 v_flow );

 // serialize the Dual Prices - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_cost.empty() )
  ::serialize< double >( group , "DualCost" , netCDF::NcDouble() , nl ,
			 v_cost );

 }  // end( DCNetworkBlockSolution::serialize( NcGroup & ) )

/*--------------------------------------------------------------------------*/

void DCNetworkBlockSolution::serialize( netCDF::NcGroup & group , int idx )
 const
{
 // call the method of the base class
 NetworkBlockSolution::serialize( group , idx );

 netCDF::NcVar FV;  // FlowValue
 netCDF::NcVar DC;  // DualCost

 if( idx < 0 ) {  // first call, have to initialize everything
  // "NumberLines" is mandatory - - - - - - - - - - - - - - - - - - - - - - -
  auto nl = group.addDim( "NumberLines" , f_number_lines );
  auto ni = group.getDim( "NumberInstants" );

  if( ! v_flow.empty() )
   FV = group.addVar( "FlowValue" , netCDF::NcDouble() , { ni , nl } );

  if( ! v_cost.empty() )
   DC = group.addVar( "DualCost" , netCDF::NcDouble() , { ni , nl } );

  idx = 0;
  }
 else {  // subsequent cqll, read what is supposedly already there
  if( ! v_flow.empty() )
   FV = group.getVar( "FlowValue" );
 
  if( ! v_cost.empty() )
   DC = group.getVar( "DualCost" );
  }
 
 std::vector< size_t > strt = { 0 , idx };
 std::vector< size_t > cnt = { f_number_nodes , f_number_instants };

 if( ! FV.isNull() )  // if Flow Variables have to be serialised
  FV.putVar( strt , cnt , v_flow.data() );

 if( ! DC.isNull() )  // if Flow Variables have to be serialised
  DC.putVar( strt , cnt , v_cost.data() );

 }  // end( DCNetworkBlockSolution::serialize( NcGroup & , int ) )

/*--------------------------------------------------------------------------*/

DCNetworkBlockSolution * DCNetworkBlockSolution::scale( double factor ) const
{
 // call the method of the base class
 auto sol = dynamic_cast< DCNetworkBlockSolution * >(
				     NetworkBlockSolution::scale( factor ) );
 assert( sol );

 if( factor == 1 )
  return( sol );

 if( ! v_flow.empty() )
  for( Index l = 0 ; l < f_number_lines ; ++l )
   sol->v_flow[ l ] *= factor;

 if( ! v_cost.empty() )
  for( Index l = 0 ; l < f_number_lines ; ++l )
   sol->v_cost[ l ] *= factor;

 return( sol );

 }  // end( DCNetworkBlockSolution::scale )

/*--------------------------------------------------------------------------*/

void DCNetworkBlockSolution::sum( const Solution * solution ,
				  double multiplier )
{
 // call the method of the base class
 NetworkBlockSolution::sum( solution , multiplier );

 auto DCNBS = dynamic_cast< const DCNetworkBlockSolution * >( solution );
 if( ! DCNBS )
  throw( std::invalid_argument(
    "DCNetworkBlockSolution::sum: solution not a DCNetworkBlockSolution" ) );

 if( f_number_lines != DCNBS->f_number_lines )
  throw( std::invalid_argument(
		"DCNetworkBlockSolution::sum: inconsistent lines number" ) );

 if( ! v_flow.empty() )
  for( Index l = 0 ; l < f_number_lines ; ++l )
   v_flow[ l ] += DCNBS->v_flow[ l ] * multiplier;

 if( ! v_cost.empty() )
  for( Index l = 0 ; l < f_number_lines ; ++l )
   v_cost[ l ] += DCNBS->v_cost[ l ] * multiplier;

 }  // end( DCNetworkBlockSolution::sum )

/*--------------------------------------------------------------------------*/

DCNetworkBlockSolution * DCNetworkBlockSolution::clone( bool empty ) const
{
 auto sol = new DCNetworkBlockSolution();

 if( ! empty ) {
  NetworkBlockSolution::guts_of_clone( sol );

  sol->f_number_lines = f_number_lines;
  sol->v_flow = v_flow;
  sol->v_cost = v_cost;
  }

 return( sol );

 }  // end( DCNetworkBlockSolution::clone )

/*--------------------------------------------------------------------------*/
/*--------------------- End File DCNetworkBlock.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
