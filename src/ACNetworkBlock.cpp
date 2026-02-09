/*--------------------------------------------------------------------------*/
/*--------------------- File ACNetworkBlock.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/

#include <map>

#include <complex>

# include <cmath>

#include "NetworkBlock.h"

#include "DCNetworkBlock.h"

#include "ACNetworkBlock.h"

#include "LinearFunction.h"

#include "OneVarConstraint.h"

#include "FRealObjective.h"

#include "DQuadFunction.h"

#include "QuadFunction.h"

#include <Eigen/Sparse>

#include <algorithm>

#ifndef PI
#define PI 3.14159265358979323846
#endif

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;
using namespace std::complex_literals;

typedef Eigen::SparseMatrix< std::complex< double > > SpCMat;
typedef Eigen::SparseVector< std::complex< double > > SpCVec;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register ACNetworkBlock to the Block factory

SMSpp_insert_in_factory_cpp_0( ACNetworkBlock );

typedef ACNetworkBlock::ACNetworkData ACNetworkData;

SMSpp_insert_in_factory_cpp_0( ACNetworkData );

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS OF DCNetworkData -------------------------*/
/*--------------------------------------------------------------------------*/

void ACNetworkData::deserialize( const netCDF::NcGroup & group ) {

 DCNetworkData::deserialize( group );

#ifndef NDEBUG
 static std::vector< std::string > expected_vars = { "ReactivePowerDemand" ,
                                                     "NodeConductance" ,
                                                     "NodeSusceptance" ,
                                                     "NodeVoltageMagnitude" ,
                                                     "NodeVoltageAngle" ,
                                                     "NodeMaxVoltage" ,
                                                     "NodeMinVoltage" ,
                                                     "LineResistance" ,
                                                     "LineReactance" ,
                                                     "LineRatio" ,
                                                     "LineRATEA" ,
                                                     "LineShiftAngle" ,
                                                     "LineMinAngle" ,
                                                     "LineMaxAngle" };
 check_variables( group , expected_vars , std::cerr );
#endif

 if( f_number_nodes > 1 ) {
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

  f_lines_type = -1;
 }
} // end( ACNetworkData::deserialize )

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS OF ACNetworkBlock ------------------------*/
/*--------------------------------------------------------------------------*/

void ACNetworkBlock::deserialize( const netCDF::NcGroup & group ) {
 DCNetworkBlock::deserialize( group );
 auto ACND = new ACNetworkData();
 ACND->deserialize( group );
 if( f_NetworkData &&
  ( f_NetworkData->get_number_nodes() != ACND->get_number_nodes() ) )
  throw( std::logic_error(
   "ACNetworkBlock::deserialize: NumberNodes not matching between NetworkData" ) );
 set_NetworkData( ACND );
} // end( ACNetworkBlock::deserialize )

/*--------------------------------------------------------------------------*/

void ACNetworkBlock::generate_abstract_variables( Configuration * stvv )
{
 const auto number_nodes = get_number_nodes();
 const auto number_lines = get_number_lines();

 // generate first the same variables as in the DCNetwork
 DCNetworkBlock::generate_abstract_variables( stvv );

 // ----- complex power flow (real and imaginary part for both directions)
 /*
 S denotes the AC power for each line, therefore it is a 2-dimensional vector
 (real and imaginary part).
 */
 v_power_flow.resize( 2 * number_lines );
 v_reactive_power_flow.resize( 2 * number_lines );
 for( Index line_id = 0 ; line_id < 2 * number_lines ; ++line_id ) {
  v_power_flow[ line_id ].set_type( ColVariable::kContinuous );
  v_reactive_power_flow[ line_id ].set_type( ColVariable::kContinuous );
 }
 add_static_variable( v_power_flow , "v_power_flow_real" );
 add_static_variable( v_reactive_power_flow , "v_reactive_power_flow" );

 // QJ: be carefull, we do not define the reverse value 
 // since the sum is symetric and the diff is anti-symetric
 v_sum_product_voltages.resize( number_lines );
 for( Index line_id = 0 ; line_id < number_lines ; ++line_id ) {
  v_sum_product_voltages[ line_id ].set_type( ColVariable::kContinuous );
 }
 v_diff_product_voltages.resize( number_lines );
 for( Index line_id = 0 ; line_id < number_lines ; ++line_id ) {
  v_diff_product_voltages[ line_id ].set_type( ColVariable::kContinuous );
 }
 v_sqrd_voltages.resize( number_nodes );
 for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {
  v_sqrd_voltages[ node_id ].set_type( ColVariable::kContinuous );
 }
 add_static_variable( v_sum_product_voltages , "v_sum_product_voltages" );
 add_static_variable( v_diff_product_voltages , "v_diff_product_voltages" );
 add_static_variable( v_sqrd_voltages , "v_sqrd_voltages" );
}

/*--------------------------------------------------------------------------*/

void ACNetworkBlock::generate_objective( Configuration * objc ) {
 if( objective_generated() ) // Objective has already been generated
  return; // nothing to do

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
} // end( ACNetworkBlock::generate_objective )

/*--------------------------------------------------------------------------*/

void ACNetworkBlock::generate_abstract_constraints( Configuration * stcc ) {
 if( constraints_generated() ) // constraints have already been generated
  return; // nothing to do

 const auto number_nodes = get_number_nodes();

 if( number_nodes <= 1 )
  return;

 const auto number_lines = get_number_lines();

 if( number_lines <= 0 )
  throw( std::logic_error( "ACNetworkBlock::generate_abstract_constraints: "
   "number of lines of DCNetworkBlock is not set" ) );

 const auto & start_line = f_NetworkData->get_start_line();
 const auto & end_line = f_NetworkData->get_end_line();

 std::vector< Index > AC_lines = f_NetworkData->get_AC_lines();
 int nb_ac_lines = AC_lines.size();

 // recover the DC lines
 std::vector< Index > DC_lines = f_NetworkData->get_DC_lines();
 int nb_dc_lines = DC_lines.size();

 double base_mva = f_NetworkData->get_baseMVA();

 // ----- Voltage bounds
 /*
 We ensure that the voltage magnitude is bounded between min_voltage and max_voltage.
 To this aim, W_{n,n} = (V_n).(V_n)^H = |V_n|^2 so we impose the bounds directly on W_{n,n}.
 */
 v_voltage_bounds_const.resize( number_nodes );
 const auto & min_voltage = f_NetworkData->get_node_min_voltage();
 const auto & max_voltage = f_NetworkData->get_node_max_voltage();
 for( Index n = 0 ; n < number_nodes ; ++n ) {
  v_voltage_bounds_const[ n ].set_lhs( pow( min_voltage[ n ] , 2 ) );
  v_voltage_bounds_const[ n ].set_rhs( pow( max_voltage[ n ] , 2 ) );
  v_voltage_bounds_const[ n ].set_variable( &v_sqrd_voltages[ n ] );
 }
 add_static_constraint( v_voltage_bounds_const , "AC_voltage_bounds_limit" );

 // ----- Angle bounds
 /*
 We aim to incorporate Phase Angle Difference (PAD) constraints for each line = (start,end):
   min_angle_{line} <= angle(V_{end}) - angle(V_{start}) <= max_angle_{line}
 This constraint can be taken into account using:
    tan(min_angle_{line}) Real(W_{start,end})<= Imag(W_{start,end}) <= max_angle_{line} Real(W_{start,end})
 */
 int i_line = 0;
 //v_angle_bounds_const.resize(boost::multi_array< FRowConstraint , 2 >::extent_gen()[ 2 ][ number_lines ] );
 v_angle_bounds_const.resize(
  boost::multi_array< FRowConstraint , 2 >::extent_gen()[ 2 ][ nb_ac_lines ] );
 const auto & min_angle = f_NetworkData->get_line_min_angle();
 const auto & max_angle = f_NetworkData->get_line_max_angle();
 //for ( Index line_id = 0 ; line_id < number_lines ; ++line_id ) {
 for( auto & line_id : AC_lines ) {
  double phi_min = PI * min_angle[ line_id ] / 180.;
  double phi_max = PI * max_angle[ line_id ] / 180.;
  Index p = start_line[ line_id ];
  Index n = end_line[ line_id ];
  // --
  auto lfunc_1 = new LinearFunction();
  lfunc_1->add_variable( &v_diff_product_voltages[ line_id ] , 1.0 );
  lfunc_1->add_variable( &v_sum_product_voltages[ line_id ] , -tan( phi_min ) );
  v_angle_bounds_const[ 0 ][ i_line ].set_lhs( 0.0 );
  v_angle_bounds_const[ 0 ][ i_line ].set_rhs( Inf< double >() );
  v_angle_bounds_const[ 0 ][ i_line ].set_function( lfunc_1 );
  // --
  auto lfunc_2 = new LinearFunction();
  lfunc_2->add_variable( &v_diff_product_voltages[ line_id ] , 1.0 );
  lfunc_2->add_variable( &v_sum_product_voltages[ line_id ] , -tan( phi_max ) );
  v_angle_bounds_const[ 1 ][ i_line ].set_lhs( -Inf< double >() );
  v_angle_bounds_const[ 1 ][ i_line ].set_rhs( 0.0 );
  v_angle_bounds_const[ 1 ][ i_line ].set_function( lfunc_2 );

  ++i_line;
 }
 if( i_line > 0 )
  add_static_constraint( v_angle_bounds_const , "AC_angle_bounds_limit" );

 // ----- Active and Reactive Power conservation:
 // Shunt admittance
 SpCVec Ys = SpCVec( number_nodes );
 for( Index n = 0 ; n < number_nodes ; ++n ) {
  double Gs = f_NetworkData->get_node_conductance().at( n ) / base_mva;
  double Bs = f_NetworkData->get_node_susceptance().at( n ) / base_mva;
  Ys.insert( n ) = Gs + 1i * Bs;
 }
 /*
 The power conservation constraint can be written
   Supply - Demand = \sum_{line_id \in L \cup L^R} power_flow[line_id]

 The complex matrix product <M,W>_F is then decomposed into a real part and an imaginary part.
 */
 v_power_flow_injection_const.resize( 2 * number_nodes );

 // real part of the power flow conservation
 for( Index p = 0 ; p < number_nodes ; ++p ) {
  auto lfunc = new LinearFunction();
  lfunc->add_variable( &v_node_injection[ 0 ][ p ] , -1.0 );
  // QJ: why 0 and not the time step ?
  lfunc->add_variable( &v_sqrd_voltages[ p ] ,
                       -Ys.coeff( p ).real() / base_mva );

  for( Index line_id = 0 ; line_id < number_lines ; ++line_id ) {
   Index i = start_line[ line_id ];
   Index j = end_line[ line_id ];

   if( i == p ) lfunc->add_variable( &v_power_flow[ line_id ] , 1. );
   if( j == p )
    lfunc->add_variable( &v_power_flow[ number_lines + line_id ] ,
                         1. );
  }
  v_power_flow_injection_const[ p ].set_both( -v_ActiveDemand[ p ] / base_mva );
  v_power_flow_injection_const[ p ].set_function( lfunc );
 }
 // imaginary part of the power flow conservation
 for( Index p = 0 ; p < number_nodes ; ++p ) {
  auto lfunc = new LinearFunction();
  lfunc->add_variable( &v_reactive_node_injection[ 0 ][ p ] , -1.0 );
  lfunc->add_variable( &v_sqrd_voltages[ p ] ,
                       -Ys.coeff( p ).imag() / base_mva );

  for( Index line_id = 0 ; line_id < number_lines ; ++line_id ) {
   Index i = start_line[ line_id ];
   Index j = end_line[ line_id ];

   if( i == p ) lfunc->add_variable( &v_reactive_power_flow[ line_id ] , 1. );
   if( j == p )
    lfunc->add_variable(
     &v_reactive_power_flow[ number_lines + line_id ] , 1. );
  }
  v_power_flow_injection_const[ number_nodes + p ].set_both(
   -v_ReactiveDemand[ p ] / base_mva );
  v_power_flow_injection_const[ number_nodes + p ].set_function( lfunc );
 }

 add_static_constraint( v_power_flow_injection_const ,
                        "AC_power_flow_injection" );

 // ----- Since the lines have been duplicated, we need to add for DC ones the link between the two versions
 v_flow_dc.resize( nb_dc_lines );
 int i_dc_line = 0;
 for( auto & line_id : DC_lines ) {
  Index p = start_line[ line_id ];
  Index n = end_line[ line_id ];
  auto lfunc = new LinearFunction();
  lfunc->add_variable( &v_power_flow[ line_id ] , 1.0 );
  lfunc->add_variable( &v_power_flow[ number_lines + line_id ] , 1.0 );
  v_flow_dc[ i_dc_line ].set_both( 0.0 );
  v_flow_dc[ i_dc_line ].set_function( lfunc );
  ++i_dc_line;
 }
 add_static_constraint( v_flow_dc , "DC_flow_links" );

 // ----- Definition of complex power flow
 /*
 We define the complex power flow S_{line} for each line = (start,end) as
   S_{start,end} = Yff_{start,end} W_{start,start} + Yft_{start,end}W_{start,end}
   S_{end,start} = Ytt_{start,end} W_{start,start} + Ytf_{start,end}W_{end,end}
 Onece, again, we then split into two constraints (one for real and one for imaginary part)
 */
 
 // shortcut to recover mathematic notation
 auto* f_net = f_NetworkData;

 auto r     = [f_net](int line_id) {return f_net->get_line_resistance().at( line_id );};
 auto x     = [f_net](int line_id) {return f_net->get_line_reactance().at( line_id );};
 auto b     = [f_net](int line_id) {return f_net->get_line_susceptance().at( line_id );};
 auto tau   = [f_net](int line_id) {return f_net->get_line_ratio().at(line_id);};
 auto theta = [f_net](int line_id) {return PI * f_net->get_line_angle().at( line_id ) / 180;};
 
 auto Y   = [r,x] (int l){ return 1.0/(r(l)+1i*x(l)); }; // common base of matrix (angle = 0, ratio = 1)
 auto Ytt = [Y,b] (int l){ return Y(l) + 0.5i*b(l); };
 auto Yff = [Ytt,tau] (int l){ return Ytt(l) / std::pow(tau(l),2.0); };
 auto Yft = [Y,theta,tau] (int l){ return Y(l) / ( tau(l)*std::exp(-1i*theta(l)) ); };
 auto Ytf = [Y,theta,tau] (int l){ return Y(l) / ( tau(l)*std::exp(1i*theta(l)) ); };

 v_voltage_definition_const.resize(boost::multi_array< FRowConstraint , 2 >::extent_gen()[ 2 ][ 2 * nb_ac_lines ] );
 i_line = 0;
 const auto splitted_lines = f_NetworkData->get_direct_and_reverse_AClines();
 for( Index p = 0 ; p < number_nodes ; ++p ) {

  // 1) direct lines
  for (auto& line_id: splitted_lines[p].first){

    // 1.1) real part
    auto lfunc_1 = new LinearFunction();
    lfunc_1->add_variable( &v_sqrd_voltages[ p ] ,
                           Yff(line_id).real() );
    lfunc_1->add_variable( &v_sum_product_voltages[ line_id ] ,
                           Yft(line_id).real() );
    lfunc_1->add_variable( &v_diff_product_voltages[ line_id ] ,
                           Yft(line_id).imag() );
    lfunc_1->add_variable( &v_power_flow[ line_id ] , -1.0 );
    v_voltage_definition_const[ 0 ][ i_line ].set_both( 0.0 );
    v_voltage_definition_const[ 0 ][ i_line ].set_function( lfunc_1 );

    // 1.2) imag part
    auto lfunc_2 = new LinearFunction();
    lfunc_2->add_variable( &v_sqrd_voltages[ p ] ,
                           -Yff(line_id).imag() );
    lfunc_2->add_variable( &v_sum_product_voltages[ line_id ] ,
                           -Yft(line_id).imag() );
    lfunc_2->add_variable( &v_diff_product_voltages[ line_id ] ,
                           Yft(line_id).real() );
    lfunc_2->add_variable( &v_reactive_power_flow[ line_id ] , -1.0 );
    v_voltage_definition_const[ 1 ][ i_line ].set_both( 0.0 );
    v_voltage_definition_const[ 1 ][ i_line ].set_function( lfunc_2 );

    ++i_line;

  }

  // 2) reverse lines
  for (auto& line_id: splitted_lines[p].second){

    // 2.1) real part
    auto lfunc_1 = new LinearFunction();
    lfunc_1->add_variable( &v_sqrd_voltages[ line_id ] ,
                           Ytt(line_id).real() );
    lfunc_1->add_variable( &v_sum_product_voltages[ line_id ] ,
                           Ytf(line_id).real() );
    lfunc_1->add_variable( &v_diff_product_voltages[ line_id ] ,
                           -Ytf(line_id).imag() ); // be carefull, diff is antisymetric
    lfunc_1->add_variable( &v_power_flow[ line_id ] , -1.0 );
    v_voltage_definition_const[ 0 ][ i_line ].set_both( 0.0 );
    v_voltage_definition_const[ 0 ][ i_line ].set_function( lfunc_1 );

    // 2.2) imag part
    auto lfunc_2 = new LinearFunction();
    lfunc_2->add_variable( &v_sqrd_voltages[ line_id ] ,
                           -Ytt(line_id).imag() );
    lfunc_2->add_variable( &v_sum_product_voltages[ line_id ] ,
                           -Ytf(line_id).imag() );
    lfunc_2->add_variable( &v_diff_product_voltages[ line_id ] ,
                           -Ytf(line_id).real() ); // be carefull, diff is antisymetric
    lfunc_2->add_variable( &v_reactive_power_flow[ line_id ] ,
                           -1.0 );
    v_voltage_definition_const[ 1 ][ i_line ].set_both( 0.0 );
    v_voltage_definition_const[ 1 ][ i_line ].set_function( lfunc_2 );  

    ++i_line; 
  }
 } // i_line should be nb_ac_lines * 2

 if( i_line > 0 )
  add_static_constraint( v_voltage_definition_const ,
                         "AC_voltage_defintion_const" );

 // ----- Thermal limit on lines
 /*
 We impose that |S_{line}| <= rateA_{line}, which corresponds to a thermal limitation.
 Then, to take into account this constraint, we use a DQuadfunction:
   Real(S_{line})^2 + Imag(S_{line})^2 <= rateA_{line}^2
 */
 v_thermal_limit.resize( 2 * number_lines );
 const auto & rate_A = f_NetworkData->get_line_rate_A();
 for( Index line_id = 0 ; line_id < number_lines ; ++line_id ) {
  Index p = start_line[ line_id ];
  Index n = end_line[ line_id ];
  auto qfunc_1 = new DQuadFunction();
  qfunc_1->add_variable( &v_power_flow[ line_id ] , 0.0 , 1.0 );
  qfunc_1->add_variable( &v_reactive_power_flow[ line_id ] , 0.0 , 1.0 );
  v_thermal_limit[ line_id ].set_lhs( -Inf< double >() );
  v_thermal_limit[ line_id ].set_rhs( pow( rate_A[ line_id ] / base_mva , 2 ) );
  v_thermal_limit[ line_id ].set_function( qfunc_1 );
  auto qfunc_2 = new DQuadFunction();
  qfunc_2->add_variable( &v_power_flow[ number_lines + line_id ] , 0.0 , 1.0 );
  qfunc_2->add_variable( &v_reactive_power_flow[ number_lines + line_id ] ,
                         0.0 , 1.0 );
  v_thermal_limit[ number_lines + line_id ].set_lhs( -Inf< double >() );
  v_thermal_limit[ number_lines + line_id ].set_rhs(
   pow( rate_A[ line_id ] / base_mva , 2 ) );
  v_thermal_limit[ number_lines + line_id ].set_function( qfunc_2 );
 }
 add_static_constraint( v_thermal_limit , "AC_thermal_limit_const" );

 // -----QJ: up to now, only SOCP relaxation is available but it could be replaced by something else
 generate_SOCP_relaxation();
}

/*--------------------------------------------------------------------------*/

/*
Links between generic variables
   v_sum_product_voltages,
   v_diff_product_voltages,
   v_sqrd_voltages
 using a SOCP relaxation.
*/
void ACNetworkBlock::generate_SOCP_relaxation( void ) {
 const auto number_nodes = get_number_nodes();
 const auto number_lines = get_number_lines();
 const auto & start_line = f_NetworkData->get_start_line();
 const auto & end_line = f_NetworkData->get_end_line();

 std::vector< Index > AC_lines = f_NetworkData->get_AC_lines();
 int nb_ac_lines = AC_lines.size();

 // ----- Voltage relaxation matrix W.
 /* We aim to impose W = V.V^H, where V is the vector of voltage for each bus/node.
 This nonlinear constraint will be replaced by SOCP relaxation.
 As the matrix W is a complex matrix, we define in the optimization model two matrices:
   - W_voltage (0<=i<n,0<=j<n) for the real part
   - W_woltage (n<=i<2n,n<=j<2n) for the imaginary part
 */

 // ----- Rotated SOCP cone for W matrix
 /*
 As we cannot take into account the true constraint W = V.V^H, we replace it by a SOCP relaxation:
   |W_{ab}|^2 <= W_{aa}W_{bb}
 As we are in complex algebra, we need auxiliary variables to write the SOCP constraints (QJ: maybe can be simplified)
 */

 // v_socp_const.resize( number_lines );
 v_socp_const.resize( nb_ac_lines );
 //for (Index line_id = 0 ; line_id < number_lines ; ++line_id ) {
 int i_line = 0;
 for( auto & line_id : AC_lines ) {
  Index p = start_line[ line_id ];
  Index n = end_line[ line_id ];
  auto qfunc = new QuadFunction();
  qfunc->add_variable( &v_sqrd_voltages[ p ] , 0.0 , 0.0 );
  qfunc->add_variable( &v_sqrd_voltages[ n ] , 0.0 , 0.0 );
  qfunc->add_nd_term( &v_sqrd_voltages[ p ] , &v_sqrd_voltages[ n ] , -1.0 );
  qfunc->add_variable( &v_sum_product_voltages[ line_id ] , 0.0 , 1.0 );
  v_socp_const[ i_line ].set_lhs( -Inf< double >() );
  v_socp_const[ i_line ].set_rhs( 0.0 );
  v_socp_const[ i_line ].set_function( qfunc );

  ++i_line;
 }
 if( i_line > 0 )
  add_static_constraint( v_socp_const , "AC_socp_const" );
}

/*--------------------------------------------------------------------------*/

std::vector< std::pair< double , double > >
ACNetworkBlock::recover_feasible_solution( void ) {
 /*
 Since the solution provided from the AC OPF relaxation problem is not necessary feasible,
 we implement a feasibility recovery algorithm.
 */

 std::vector< std::pair< double , double > > v_feasible_sol;
 // Each pair is the real and imaginary part

 // 1) First, get the solution of the relaxation problem
 std::vector< double > relaxed_power_flow;
 std::vector< double > relaxed_reactive_power_flow;

 std::transform( v_power_flow.begin() , v_power_flow.end() ,
                 relaxed_power_flow.begin() ,
                 []( ColVariable v ) { return( v.get_value() ); }
 );
 std::transform( v_reactive_power_flow.begin() , v_reactive_power_flow.end() ,
                 relaxed_reactive_power_flow.begin() ,
                 []( ColVariable v ) { return( v.get_value() ); }
 );

 // 2) Then compute spanning tree
 auto result = f_NetworkData->get_cycle_basis();

 // 3) Do some magic (TODO)

 return( v_feasible_sol );
}

/*--------------------------------------------------------------------------*/
/*--------------------- End File ACNetworkBlock.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
