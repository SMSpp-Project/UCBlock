/*--------------------------------------------------------------------------*/
/*--------------------- File ACNetworkBlock.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/

#include <map>

#include <complex>

#include <cmath>

#include "ACNetworkBlock.h"

#include "LinearFunction.h"

#include "OneVarConstraint.h"

#include "FRealObjective.h"

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

using SpCMat = Eigen::SparseMatrix< std::complex< double > >;
using SpCVec = Eigen::SparseVector< std::complex< double > >;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/
// register ACNetworkBlock to the Block factory

SMSpp_insert_in_factory_cpp_1( ACNetworkBlock );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// register ACNetworkBlock::ACNetworkData to the NetworkData factory

using ACNetworkData = ACNetworkBlock::ACNetworkData ;

SMSpp_insert_in_factory_cpp_0( ACNetworkData );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// register ACNetworkBlockSolution to the Solution factory

SMSpp_insert_in_factory_cpp_0( ACNetworkBlockSolution );

/*--------------------------------------------------------------------------*/
/*--------------------------- STATIC FUNCTIONS -----------------------------*/
/*--------------------------------------------------------------------------*/

// A rounding function
static inline double round_sig( double value , int digits = 16 )
{
 if( value == 0.0 ) return( 0.0 );

 double abs_v = std::fabs( value );
 int exponent = static_cast< int >( std::floor( std::log10( abs_v ) ) );
 double factor = std::pow( 10.0 , digits - 1 - exponent );

 return( std::round( value * factor ) / factor );
 }

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS OF ACNetworkData -------------------------*/
/*--------------------------------------------------------------------------*/

void ACNetworkData::deserialize( const netCDF::NcGroup & group )
{
#ifndef NDEBUG
 // check all expected variables, comprised those of the base class: see
 // DCNetworkData::deserialize() for the rationale
 static std::vector< std::string > expected_vars = { "ActiveDemand" ,
  "StartLine" , "EndLine" , "HyperArcID" , "MinPowerFlow" , "MaxPowerFlow" ,
  "LineSusceptance" , "NetworkCost" , "NodeName" , "LineName" ,
  "ConstantTerm" , "Efficiency" ,
  // ACNetworkData
  "ReactivePowerDemand" , "NodeConductance" , "NodeSusceptance" ,
  "NodeVoltageMagnitude" , "NodeVoltageAngle" , "NodeMaxVoltage" ,
  "NodeMinVoltage" , "LineResistance" , "LineReactance" , "LineRatio" ,
  "LineRATEA" , "LineShiftAngle" , "LineMinAngle" , "LineMaxAngle" ,
  // if called from UCBlock:
  "ActivePowerDemand" , "GeneratorNode" , "NetworkConstantTerms" ,
  "NetworkBlockClassname" , "NetworkDataClassname"
  };

 check_variables( group , expected_vars , std::cerr );
#endif

 DCNetworkData::deserialize( group );

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
 }  // end( ACNetworkData::deserialize )

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS OF ACNetworkBlock ------------------------*/
/*--------------------------------------------------------------------------*/

void ACNetworkBlock::deserialize( const netCDF::NcGroup & group )
{
 DCNetworkBlock::deserialize( group );

 auto ACND = new ACNetworkData();
 ACND->deserialize( group );
 if( f_NetworkData &&
  ( f_NetworkData->get_number_nodes() != ACND->get_number_nodes() ) )
  throw( std::logic_error( "ACNetworkBlock::deserialize: NumberNodes not "
			   "matching between NetworkData" ) );
 set_NetworkData( ACND );

 }  // end( ACNetworkBlock::deserialize )

/*--------------------------------------------------------------------------*/

void ACNetworkBlock::generate_abstract_variables( Configuration * stvv )
{
 if( variables_generated() )  // variables have already been generated
  return;                     // nothing to do

 DCNetworkBlock::generate_abstract_variables( stvv );

 const auto number_nodes = get_number_nodes();
 const auto number_lines = get_number_lines();

 // the node injection variables
 v_reactive_node_injection.resize( number_nodes );
 for( Index node_id = 0 ; node_id < number_nodes ; ++node_id )
  v_reactive_node_injection[ node_id ].set_type( ColVariable::kContinuous );

 add_static_variable( v_reactive_node_injection , "reactive_s_network" );

 // ----- complex power flow (real and imaginary part for both directions)
 /*
 S denotes the AC power for each line, therefore it is a 2-dimensional vector
 (real and imaginary part).
 */
 v_power_flow.resize( 2 * number_lines );
 for( Index line_id = 0 ; line_id < 2 * number_lines ; ++line_id )
  v_power_flow[ line_id ].set_type( ColVariable::kContinuous );
 add_static_variable( v_power_flow , "v_power_flow_real" );

 v_reactive_power_flow.resize( 2 * number_lines );
 for( Index line_id = 0 ; line_id < 2 * number_lines ; ++line_id )
  v_reactive_power_flow[ line_id ].set_type( ColVariable::kContinuous );
 add_static_variable( v_reactive_power_flow , "v_reactive_power_flow" );

 // be careful, we do not define the reverse value
 // since the sum is symmetric and the diff is anti-symmetric
 v_sum_product_voltages.resize( number_lines );
 for( Index line_id = 0 ; line_id < number_lines ; ++line_id )
  v_sum_product_voltages[ line_id ].set_type( ColVariable::kContinuous );
 add_static_variable( v_sum_product_voltages , "v_sum_product_voltages" );

 v_diff_product_voltages.resize( number_lines );
 for( Index line_id = 0 ; line_id < number_lines ; ++line_id )
  v_diff_product_voltages[ line_id ].set_type( ColVariable::kContinuous );
 add_static_variable( v_diff_product_voltages , "v_diff_product_voltages" );

 v_sqrd_voltages.resize( number_nodes );
 for( Index node_id = 0 ; node_id < number_nodes ; ++node_id )
  v_sqrd_voltages[ node_id ].set_type( ColVariable::kContinuous );
 add_static_variable( v_sqrd_voltages , "v_sqrd_voltages" );

 }  // end( ACNetworkBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void ACNetworkBlock::generate_abstract_constraints( Configuration * stcc )
{
 if( constraints_generated() )  // constraints have already been generated
  return;                       // nothing to do

 DCNetworkBlock::generate_abstract_constraints( stcc );

 const auto number_nodes = get_number_nodes();

 if( number_nodes <= 1 )
  return;

 const auto number_lines = get_number_lines();

 if( number_lines <= 0 )
  throw( std::logic_error( "ACNetworkBlock::generate_abstract_constraints: "
                           "number of lines of DCNetworkBlock is not set" )
	 );

 // node injection bound constraints
 reactive_node_injection_bounds_const.resize( number_nodes );
 for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {
  reactive_node_injection_bounds_const[ node_id ].set_lhs(
				   v_MinReactiveNodeInjection[ node_id ] );
  reactive_node_injection_bounds_const[ node_id ].set_rhs(
				   v_MaxReactiveNodeInjection[ node_id ] );
  reactive_node_injection_bounds_const[ node_id ].set_variable(
				  & v_reactive_node_injection[ node_id ] );
  }

 add_static_constraint( reactive_node_injection_bounds_const ,
                        "Reactive_Node_Injection_Bound_Const_Network" );

 // now start processing lines
 const auto & start_line = f_NetworkData->get_start_line();
 const auto & end_line = f_NetworkData->get_end_line();

 std::vector< Index > AC_lines = f_NetworkData->get_AC_lines();
 int nb_ac_lines = AC_lines.size();

 // recover the DC lines
 std::vector< Index > DC_lines = f_NetworkData->get_DC_lines();
 int nb_dc_lines = DC_lines.size();

 double base_mva = f_NetworkData->get_baseMVA();
 // --- scaling the v_power_flow and v_reactive_power_flow to improve numerical stability
 //     effectively we are swapping out v_power_flow for v_power_flow_tilde with
 //                    v_power_flow_tilde = C * v_power_flow
 //     and likewise v_reactive_power_flow
 constexpr double C_v_scal = 1.0; /* e.g. 100.0 */

 // ----- Voltage bounds
 /*
 We ensure that the voltage magnitude is bounded between min_voltage and max_voltage.
 To this aim, W_{n,n} = (V_n).(V_n)^H = |V_n|^2 so we impose the bounds directly on W_{n,n}.
 */
 v_voltage_bounds_const.resize( number_nodes );
 const auto & min_voltage = ND()->get_node_min_voltage();
 const auto & max_voltage = ND()->get_node_max_voltage();
 for( Index n = 0 ; n < number_nodes ; ++n ) {
  v_voltage_bounds_const[ n ].set_lhs( pow( C_v_scal * min_voltage[ n ] ,
					    2 ) );
  v_voltage_bounds_const[ n ].set_rhs( pow( C_v_scal * max_voltage[ n ] ,
					    2 ) );
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
 v_angle_bounds_const.resize( MAFRC_ext()[ 2 ][ nb_ac_lines ] );
 const auto & min_angle = ND()->get_line_min_angle();
 const auto & max_angle = ND()->get_line_max_angle();

 for( auto & line_id : AC_lines ) {
  double phi_min = PI * min_angle[ line_id ] / 180.;
  double phi_max = PI * max_angle[ line_id ] / 180.;
  // --
  auto lfunc_1 = new LinearFunction();
  lfunc_1->add_variable( & v_diff_product_voltages[ line_id ] , 1.0 );
  lfunc_1->add_variable( & v_sum_product_voltages[ line_id ] ,
			 -tan( phi_min ) );
  v_angle_bounds_const[ 0 ][ i_line ].set_lhs( 0.0 );
  v_angle_bounds_const[ 0 ][ i_line ].set_rhs( Inf< double >() );
  v_angle_bounds_const[ 0 ][ i_line ].set_function( lfunc_1 );
  // --
  auto lfunc_2 = new LinearFunction();
  lfunc_2->add_variable( & v_diff_product_voltages[ line_id ] , 1.0 );
  lfunc_2->add_variable( & v_sum_product_voltages[ line_id ] ,
			 -tan( phi_max ) );
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
  double Gs = ND()->get_node_conductance().at( n ) / base_mva;
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
  // why 0 and not the time step ?
  lfunc->add_variable( &v_sqrd_voltages[ p ] ,
                       -Ys.coeff( p ).real() / base_mva );

  for( Index line_id = 0 ; line_id < number_lines ; ++line_id ) {
   Index i = start_line[ line_id ];
   Index j = end_line[ line_id ];

   if( i == p )
    lfunc->add_variable( &v_power_flow[ line_id ] , 1. / C_v_scal );
   if( j == p )
    lfunc->add_variable( &v_power_flow[ number_lines + line_id ] ,
                         1. / C_v_scal );
   }
  v_power_flow_injection_const[ p ].set_both( -v_ActiveDemand[ p ] /
					      base_mva );
  v_power_flow_injection_const[ p ].set_function( lfunc );
  }

 // imaginary part of the power flow conservation
 for( Index p = 0 ; p < number_nodes ; ++p ) {
  auto lfunc = new LinearFunction();
  lfunc->add_variable( & v_reactive_node_injection[ p ] , -1.0 );
  lfunc->add_variable( & v_sqrd_voltages[ p ] ,
                       -Ys.coeff( p ).imag() / base_mva );

  for( Index line_id = 0 ; line_id < number_lines ; ++line_id ) {
   Index i = start_line[ line_id ];
   Index j = end_line[ line_id ];

   if( i == p ) lfunc->add_variable( &v_reactive_power_flow[ line_id ] ,
                                     1. / C_v_scal );
   if( j == p )
    lfunc->add_variable(
     &v_reactive_power_flow[ number_lines + line_id ] , 1. / C_v_scal );
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
  auto lfunc = new LinearFunction();
  lfunc->add_variable( &v_power_flow[ line_id ] , 1.0 / C_v_scal );
  lfunc->add_variable( &v_power_flow[ number_lines + line_id ] ,
                       1.0 / C_v_scal );
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
 Once, again, we then split into two constraints (one for real and one for imaginary part)
 */

 // shortcut to recover mathematical notation
 auto * f_net = static_cast< ACNetworkData * >( f_NetworkData );

 auto r = [ f_net ]( int line_id ) {
  return( f_net->get_line_resistance().at( line_id ) );
 };
 auto x = [ f_net ]( int line_id ) {
  return( f_net->get_line_reactance().at( line_id ) );
 };
 auto b = [ f_net ]( int line_id ) {
  return( f_net->get_line_susceptance().at( line_id ) );
 };
 auto tau = [ f_net ]( int line_id ) {
  return( f_net->get_line_ratio().at( line_id ) );
 };
 auto theta = [ f_net ]( int line_id ) {
  return( PI * f_net->get_line_angle().at( line_id ) / 180 );
 };

 auto Y = [ r , x ]( int l ) { return( 1.0 / ( r( l ) + 1i * x( l ) ) ); };
 // common base of matrix (angle = 0, ratio = 1)
 auto Ytt = [ Y , b ]( int l ) { return( Y( l ) + 0.5i * b( l ) ); };
 auto Yff = [ Ytt , tau ]( int l ) {
  return( Ytt( l ) / std::pow( tau( l ) , 2.0 ) );
 };
 auto Yft = [ Y , theta , tau ]( int l ) {
  return( -1.0 * Y( l ) / ( tau( l ) * std::exp( -1i * theta( l ) ) ) );
 };
 auto Ytf = [ Y , theta , tau ]( int l ) {
  return( -1.0 * Y( l ) / ( tau( l ) * std::exp( 1i * theta( l ) ) ) );
 };

 v_voltage_definition_const.resize( MAFRC_ext()[ 2 ][ 2 * nb_ac_lines ] );
 i_line = 0;
 const auto splitted_lines = ND()->get_direct_and_reverse_AClines();
 for( Index p = 0 ; p < number_nodes ; ++p ) {
  // --- Slack for AC_voltage_definition_const
  constexpr double f_ACvS = 0.0;
  // --- scaling constant for the AC_voltage_definition_const equations (to improve numeric stability)
  constexpr double f_scale = 1.0;
  // 1) direct lines
  for( auto & line_id : splitted_lines[ p ].first ) {
   // 1.0) Verification that line ratio is not zero as this implies NaN coefficients
   if( ! std::isfinite( Yff( line_id ).real() ) )
    throw( std::logic_error( "Non-finite coefficient, possibly line ratio is "
                             "zero for an AC line" ) );

   // 1.1) real part
   auto lfunc_1 = new LinearFunction();
   lfunc_1->add_variable( &v_sqrd_voltages[ p ] ,
                          round_sig( Yff( line_id ).real() * f_scale ) );
   lfunc_1->add_variable( &v_sum_product_voltages[ line_id ] ,
                          round_sig( Yft( line_id ).real() * f_scale ) );
   lfunc_1->add_variable( &v_diff_product_voltages[ line_id ] ,
                          round_sig( Yft( line_id ).imag() * f_scale ) );
   lfunc_1->add_variable( &v_power_flow[ line_id ] ,
                          -( 1.0 / C_v_scal ) * f_scale );
   //v_voltage_definition_const[ 0 ][ i_line ].set_both( 0.0 );
   v_voltage_definition_const[ 0 ][ i_line ].set_lhs( -f_ACvS );
   v_voltage_definition_const[ 0 ][ i_line ].set_rhs( f_ACvS );
   v_voltage_definition_const[ 0 ][ i_line ].set_function( lfunc_1 );

   // 1.2) imag part
   auto lfunc_2 = new LinearFunction();
   lfunc_2->add_variable( &v_sqrd_voltages[ p ] ,
                          round_sig( -Yff( line_id ).imag() * f_scale ) );
   lfunc_2->add_variable( &v_sum_product_voltages[ line_id ] ,
                          round_sig( -Yft( line_id ).imag() * f_scale ) );
   lfunc_2->add_variable( &v_diff_product_voltages[ line_id ] ,
                          round_sig( Yft( line_id ).real() * f_scale ) );
   lfunc_2->add_variable( &v_reactive_power_flow[ line_id ] ,
                          -( 1.0 / C_v_scal ) * f_scale );
   //v_voltage_definition_const[ 1 ][ i_line ].set_both( 0.0 );
   v_voltage_definition_const[ 1 ][ i_line ].set_lhs( -f_ACvS );
   v_voltage_definition_const[ 1 ][ i_line ].set_rhs( f_ACvS );
   v_voltage_definition_const[ 1 ][ i_line ].set_function( lfunc_2 );

   ++i_line;
  }

  // 2) reverse lines
  for( auto & line_id : splitted_lines[ p ].second ) {
   // 2.1) real part
   auto lfunc_1 = new LinearFunction();
   lfunc_1->add_variable( &v_sqrd_voltages[ p ] ,
                          round_sig( Ytt( line_id ).real() * f_scale ) );
   lfunc_1->add_variable( &v_sum_product_voltages[ line_id ] ,
                          round_sig( Ytf( line_id ).real() * f_scale ) );
   lfunc_1->add_variable( &v_diff_product_voltages[ line_id ] ,
                          round_sig( -Ytf( line_id ).imag() * f_scale ) );
   // be careful, diff is antisymmetric
   lfunc_1->add_variable( &v_power_flow[ number_lines + line_id ] ,
                          -( 1.0 / C_v_scal ) * f_scale );
   //v_voltage_definition_const[ 0 ][ i_line ].set_both( 0.0 );
   v_voltage_definition_const[ 0 ][ i_line ].set_lhs( -f_ACvS );
   v_voltage_definition_const[ 0 ][ i_line ].set_rhs( f_ACvS );
   v_voltage_definition_const[ 0 ][ i_line ].set_function( lfunc_1 );

   // 2.2) imag part
   auto lfunc_2 = new LinearFunction();
   lfunc_2->add_variable( &v_sqrd_voltages[ p ] ,
                          round_sig( -Ytt( line_id ).imag() * f_scale ) );
   lfunc_2->add_variable( &v_sum_product_voltages[ line_id ] ,
                          round_sig( -Ytf( line_id ).imag() * f_scale ) );
   lfunc_2->add_variable( &v_diff_product_voltages[ line_id ] ,
                          round_sig( -Ytf( line_id ).real() * f_scale ) );
   // be careful, diff is antisymmetric
   lfunc_2->add_variable( &v_reactive_power_flow[ number_lines + line_id ] ,
                          -( 1.0 / C_v_scal ) * f_scale );
   //v_voltage_definition_const[ 1 ][ i_line ].set_both( 0.0 );
   v_voltage_definition_const[ 1 ][ i_line ].set_lhs( -f_ACvS );
   v_voltage_definition_const[ 1 ][ i_line ].set_rhs( f_ACvS );
   v_voltage_definition_const[ 1 ][ i_line ].set_function( lfunc_2 );

   ++i_line;
  }
 } // i_line should be nb_ac_lines * 2

 if( i_line > 0 )
  add_static_constraint( v_voltage_definition_const ,
                         "AC_voltage_definition_const" );

 // ----- Thermal limit on lines
 /*
 We impose that |S_{line}| <= rateA_{line}, which corresponds to a thermal limitation.
 Then, to take into account this constraint, we use a DQuadFunction:
   Real(S_{line})^2 + Imag(S_{line})^2 <= rateA_{line}^2
 */
 v_thermal_limit.resize( 2 * number_lines );
 const auto & rate_A = ND()->get_line_rate_A();
 for( Index line_id = 0 ; line_id < number_lines ; ++line_id ) {
  auto qfunc_1 = new DQuadFunction();
  qfunc_1->add_variable( &v_power_flow[ line_id ] , 0.0 , 1.0 );
  qfunc_1->add_variable( &v_reactive_power_flow[ line_id ] , 0.0 , 1.0 );
  v_thermal_limit[ line_id ].set_lhs( -Inf< double >() );
  v_thermal_limit[ line_id ].set_rhs(
   pow( C_v_scal * rate_A[ line_id ] / base_mva , 2 ) );
  v_thermal_limit[ line_id ].set_function( qfunc_1 );
  auto qfunc_2 = new DQuadFunction();
  qfunc_2->add_variable( &v_power_flow[ number_lines + line_id ] , 0.0 , 1.0 );
  qfunc_2->add_variable( &v_reactive_power_flow[ number_lines + line_id ] ,
                         0.0 , 1.0 );
  v_thermal_limit[ number_lines + line_id ].set_lhs( -Inf< double >() );
  v_thermal_limit[ number_lines + line_id ].set_rhs(
   pow( C_v_scal * rate_A[ line_id ] / base_mva , 2 ) );
  v_thermal_limit[ number_lines + line_id ].set_function( qfunc_2 );
 }
 add_static_constraint( v_thermal_limit , "AC_thermal_limit_const" );

 // up to now, only SOCP relaxation is available, but it could be replaced by something else
 generate_SOCP_relaxation();
} // end( ACNetworkBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/
/*
Links between generic variables
   v_sum_product_voltages,
   v_diff_product_voltages,
   v_sqrd_voltages
 using a SOCP relaxation.
*/

void ACNetworkBlock::generate_SOCP_relaxation( void )
{
 const auto & start_line = f_NetworkData->get_start_line();
 const auto & end_line = f_NetworkData->get_end_line();

 std::vector< Index > AC_lines = f_NetworkData->get_AC_lines();
 int nb_ac_lines = AC_lines.size();

 // ----- Voltage relaxation matrix W.
 /* We aim to impose W = V.V^H, where V is the vector of voltage for each bus/node.
 This nonlinear constraint will be replaced by SOCP relaxation.
 As the matrix W is a complex matrix, we define in the optimization model two matrices:
   - W_voltage (0<=i<n,0<=j<n) for the real part
   - W_voltage (n<=i<2n,n<=j<2n) for the imaginary part
 */

 // ----- Rotated SOCP cone for W matrix
 /*
 As we cannot take into account the true constraint W = V.V^H, we replace it by a SOCP relaxation:
   |W_{ab}|^2 <= W_{aa}W_{bb}
 As we are in complex algebra, we need auxiliary variables to write the SOCP constraints (maybe can be simplified)
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
} // end( ACNetworkBlock::generate_SOCP_relaxation )

/*--------------------------------------------------------------------------*/

std::vector< std::pair< double , double > >
                            ACNetworkBlock::recover_feasible_solution( void )
{
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
                 []( const ColVariable & v ) { return( v.get_value() ); }
 );
 std::transform( v_reactive_power_flow.begin() , v_reactive_power_flow.end() ,
                 relaxed_reactive_power_flow.begin() ,
                 []( const ColVariable & v ) { return( v.get_value() ); }
 );

 // 1b) Multiply back the obtained solutions by the earlier scale factor since
 //     indeed we have computed v_power_flow_tilde, and we care for
 //     v_power_flow -> the relation is v_power_flow_tilde = C * v_power_flow

 // 2) Then compute spanning tree
 auto result = f_NetworkData->get_cycle_basis();

 // 3) Do some magic (TODO)

 return( v_feasible_sol );

 } // end( ACNetworkBlock::recover_feasible_solution )

/*--------------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/

Solution * ACNetworkBlock::get_Solution( Configuration * csolc , bool emptys )
{
 Index wsol = 7;
 if( ( ! csolc ) && f_BlockConfig )
  csolc = f_BlockConfig->f_solution_Configuration;

 if( auto config = dynamic_cast< SimpleConfiguration< int > * >( csolc ) )
  wsol = config->f_value;

 // call the method of the base class
 auto * sol = dynamic_cast< ACNetworkBlockSolution * >(
			    DCNetworkBlock::get_Solution( csolc , emptys ) );
 assert( sol );

 if( wsol & 1 )
  sol->v_node_injection_reactive.resize( get_number_nodes() );

 if( wsol & 2 ) {
  sol->v_reactive_flow_from.resize( get_number_lines() );
  sol->v_reactive_flow_to.resize( get_number_lines() );
  }

 if( ! emptys )
  sol->read( this );

 return( sol );

 }  // end( ACNetworkBlock::get_Solution )

/*--------------------------------------------------------------------------*/

NetworkBlockSolution * ACNetworkBlock::new_Solution( void ) const {
 return( new ACNetworkBlockSolution() );
 }

/*--------------------------------------------------------------------------*/
/*------------------ METHODS OF ACNetworkBlockSolution ---------------------*/
/*--------------------------------------------------------------------------*/

void ACNetworkBlockSolution::deserialize( const netCDF::NcGroup & group )
{
 // call the method of the base class
 DCNetworkBlockSolution::deserialize( group );

 // deserialize the (reactive) flow injection - - - - - - - - - - - - - - - -
 ::deserialize< double >( group , "NodeInjectionReactive" ,
			  v_node_injection_reactive , false );

 // deserialize the (reactive) Flow Variables - - - - - - - - - - - - - - - -
 ::deserialize< double >( group , "ReactiveFlowFromValue" ,
			  v_reactive_flow_from , false );

 ::deserialize< double >( group , "ReactiveFlowToValue" ,
			  v_reactive_flow_to , false );

 }  // end( ACNetworkBlockSolution::deserialize( NcGroup & ) )

/*--------------------------------------------------------------------------*/

void ACNetworkBlockSolution::deserialize( const netCDF::NcGroup & group ,
                                          size_t idx )
{
 // call the method of the base class
 DCNetworkBlockSolution::deserialize( group , idx );

 std::vector< size_t > strt = { idx , 0 };

 // deserialize the (reactive) flow injection - - - - - - - - - - - - - - - -
 auto ncVar = group.getVar( "NodeInjectionReactive" );
 if( ncVar.isNull() )
  v_node_injection_reactive.clear();
 else {
  std::vector< size_t > cnt = { 1 , f_number_nodes };
  v_node_injection_reactive.resize( f_number_nodes );
  ncVar.getVar( strt , cnt , v_node_injection_reactive.data() );
  }

 // deserialize the (reactive) Flow Variables- - - - - - - - - - - - - - - -
 ncVar = group.getVar( "ReactiveFlowFromValue" );
 if( ncVar.isNull() ) {
  v_reactive_flow_from.clear();
  v_reactive_flow_to.clear();
  }
 else {
  std::vector< size_t > cnt = { 1 , f_number_lines };
  v_reactive_flow_from.resize( f_number_lines );
  v_reactive_flow_to.resize( f_number_lines );
  ncVar.getVar( strt , cnt , v_reactive_flow_from.data() );

  ncVar = group.getVar( "ReactiveFlowToValue" );
  if( ncVar.isNull() )
   throw( std::logic_error( "DCNetworkBlockSolution::deserialize( idx ): "
			    "ReactiveFlowFromValue present but "
			    "ReactiveFlowToValue not" ) );

  ncVar.getVar( strt , cnt , v_reactive_flow_to.data() );
  }
 }  // end( ACNetworkBlockSolution::deserialize( NcGroup & , size_t ) )

/*--------------------------------------------------------------------------*/

void ACNetworkBlockSolution::read( const Block * block )
{
 // call the method of the base class
 DCNetworkBlockSolution::read( block );

 auto ACNB = dynamic_cast< const ACNetworkBlock * >( block );
 if( ! ACNB )
  throw( std::invalid_argument(
           "ACNetworkBlockSolution::read: block is not a ACNetworkBlock" ) );

 if( ! v_node_injection_reactive.empty() ) {
  // read the (reactive) node injection variables- - - - - - - - - - - - - -
  auto RNI = ACNB->get_const_reactive_node_injection();
  for( Index n = 0 ; n < f_number_nodes ; ++n )
   v_node_injection_reactive[ n ] = (*(RNI++)).get_value();
  }

 if( ! v_reactive_flow_from.empty() ) {
  // read the (reactive) flow power variables- - - - - - - - - - - - - - - -
  auto RFli = ACNB->get_const_reactive_power_flow().begin();
  for( Index l = 0 ; l < f_number_lines ; ++l )
   v_reactive_flow_from[ l ] = (*(RFli++)).get_value();
  for( Index l = 0 ; l < f_number_lines ; ++l )
   v_reactive_flow_to[ l ] = (*(RFli++)).get_value();
  }
 }  // end( ACNetworkBlockSolution::read )

/*--------------------------------------------------------------------------*/

void ACNetworkBlockSolution::write( Block * block )
{
 // call the method of the base class
 DCNetworkBlockSolution::write( block );

 auto ACNB = dynamic_cast< ACNetworkBlock * >( block );
 if( ! ACNB )
  throw( std::invalid_argument(
          "ACNetworkBlockSolution::write: block is not a DCNetworkBlock" ) );

 if( ! v_node_injection_reactive.empty() ) {
  // write the (reactive) node injection variables - - - - - - - - - - - - -
  auto RNI = ACNB->get_reactive_node_injection();
  for( Index n = 0 ; n < f_number_nodes ; ++n )
   (*(RNI++)).set_value( v_node_injection_reactive[ n ] );
  }

 if( ! v_reactive_flow_from.empty() ) {
  // write the (reactive) flow power variables - - - - - - - - - - - - - - -
  auto RFli = ACNB->get_reactive_power_flow().begin();
  for( Index l = 0 ; l < f_number_lines ; ++l )
   (*(RFli++)).set_value( v_reactive_flow_from[ l ] );
  for( Index l = 0 ; l < f_number_lines ; ++l )
   (*(RFli++)).set_value( v_reactive_flow_to[ l ] );
  }
 }  // end( ACNetworkBlockSolution::write )

/*--------------------------------------------------------------------------*/

void ACNetworkBlockSolution::serialize( netCDF::NcGroup & group ) const
{
 // call the method of the base class
 DCNetworkBlockSolution::serialize( group );

 // serialize the (reactive) flow injection - - - - - - - - - - - - - - - - -
 if( ! v_node_injection_reactive.empty() ) {
  auto nn = group.getDim( "NumberNodes" );
  ::serialize< double >( group , "NodeInjectionReactive" ,
			 netCDF::NcDouble() , nn ,
			 v_node_injection_reactive );
  }

 // serialize the (reactive) Flow Variables - - - - - - - - - - - - - - - - -
 if( ! v_reactive_flow_from.empty() ) {
  auto nl = group.getDim( "NumberLines" );
  ::serialize< double >( group , "ReactiveFlowFromValue" ,
			 netCDF::NcDouble() , nl , v_reactive_flow_from );
  ::serialize< double >( group , "ReactiveFlowTOValue" ,
			 netCDF::NcDouble() , nl , v_reactive_flow_to );
  }
 }  // end( ACNetworkBlockSolution::serialize( NcGroup & ) )

/*--------------------------------------------------------------------------*/

void ACNetworkBlockSolution::serialize( netCDF::NcGroup & group ,
                                        size_t idx ) const {
 // call the method of the base class
 DCNetworkBlockSolution::serialize( group , idx );

 // now serialize the data structures - - - - - - - - - - - - - - - - - - - -

 netCDF::NcVar RNI ;  // NodeInjectionReactive
 netCDF::NcVar RFFV;  // ReactiveFlowFromValue
 netCDF::NcVar RFTV;  // ReactiveFlowToValue

 if( idx == 0 ) {  // first call, have to initialize everything
  // "NumberNetworks" is mandatory, and it's checked in the base class
  auto nnw = group.getDim( "NumberNetworks" );

  if( ! v_node_injection_reactive.empty() ) {
   auto nn = group.getDim( "NumberNodes" );

   RNI = group.addVar( "NodeInjectionReactive" , netCDF::NcDouble() ,
		       { nnw , nn } );
  }

  if( ! v_reactive_flow_from.empty() ) {
   auto nl = group.getDim( "NumberLines" );

   RFFV = group.addVar( "ReactiveFlowFromValue" , netCDF::NcDouble() ,
			{ nnw , nl } );
   RFTV = group.addVar( "ReactiveFlowToValue" , netCDF::NcDouble() ,
			{ nnw , nl } );
   }
  }
 else {  // subsequent call, read what is supposedly already there
  if( ! v_node_injection_reactive.empty() )
   RNI =  group.getVar( "NodeInjectionReactive" );

  if( ! v_reactive_flow_from.empty() ) {
   RFFV = group.getVar( "ReactiveFlowFromValue" );
   RFTV = group.getVar( "ReactiveFlowToValue" );
   }
  }

 std::vector< size_t > strt = { idx , 0 };

 if( ! RFFV.isNull() ) {  // reactive power injection have to be serialised
  std::vector< size_t > cnt = { 1 , f_number_nodes };
  RNI.putVar( strt , cnt , v_node_injection_reactive.data() );
  }

 if( ! RFFV.isNull() ) {  // reactive power flows have to be serialised
  std::vector< size_t > cnt = { 1 , f_number_lines };
  RFFV.putVar( strt , cnt , v_reactive_flow_from.data() );
  RFTV.putVar( strt , cnt , v_reactive_flow_to.data() );
  }
 }  // end( ACNetworkBlockSolution::serialize( NcGroup & , size_t ) )

/*--------------------------------------------------------------------------*/

ACNetworkBlockSolution * ACNetworkBlockSolution::scale( double factor ) const
{
 // call the method of the base class, which calls that of
 // NetworkBlockSolution which calls clone() and therefore returns a
 // ACNetworkBlockSolution
 auto sol = dynamic_cast< ACNetworkBlockSolution * >(
                                   DCNetworkBlockSolution::scale( factor ) );
 assert( sol );

 if( factor == 1 )
  return( sol );

 if( ! v_node_injection_reactive.empty() )
  for( Index n = 0 ; n < f_number_nodes ; ++n )
   sol->v_node_injection_reactive[ n ] *= factor;

 if( ! v_reactive_flow_from.empty() )
  for( Index l = 0 ; l < f_number_lines ; ++l ) {
   sol->v_reactive_flow_from[ l ] *= factor;
   sol->v_reactive_flow_to[ l ] *= factor;
   }

 return( sol );

 }  // end( ACNetworkBlockSolution::scale )

/*--------------------------------------------------------------------------*/

void ACNetworkBlockSolution::sum( const Solution * solution ,
                                  double multiplier )
{
 // call the method of the base class
 DCNetworkBlockSolution::sum( solution , multiplier );

 auto ACNBS = dynamic_cast< const ACNetworkBlockSolution * >( solution );
 if( ! ACNBS )
  throw( std::invalid_argument(
    "ACNetworkBlockSolution::sum: solution not a ACNetworkBlockSolution" ) );

 if( ! v_node_injection_reactive.empty() )
  for( Index n = 0 ; n < f_number_nodes ; ++n )
   v_node_injection_reactive[ n ] += ACNBS->v_node_injection_reactive[ n ]
                                     * multiplier;

 if( ! v_reactive_flow_from.empty() )
  for( Index l = 0 ; l < f_number_lines ; ++l ) {
   v_reactive_flow_from[ l ] += ACNBS->v_reactive_flow_from[ l ] * multiplier;
   v_reactive_flow_to[ l ] += ACNBS->v_reactive_flow_to[ l ] * multiplier;
   }

 }  // end( ACNetworkBlockSolution::sum )

/*--------------------------------------------------------------------------*/

ACNetworkBlockSolution * ACNetworkBlockSolution::clone( bool empty ) const
{
 auto sol = new ACNetworkBlockSolution();

 if( ! empty ) {
  NetworkBlockSolution::guts_of_clone( sol );

  //!! kludge: just copy-paste the relevant lines of
  //!! DCNetworkBlockSolution::clone(), better solution needed
  sol->f_number_lines = f_number_lines;
  sol->v_flow = v_flow;
  sol->v_cost = v_cost;
  //!! end of kludge

  sol->v_node_injection_reactive = v_node_injection_reactive;
  sol->v_reactive_flow_from = v_reactive_flow_from;
  sol->v_reactive_flow_to = v_reactive_flow_to;
  }

 return( sol );

 }  // end( ACNetworkBlockSolution::clone )

/*--------------------------------------------------------------------------*/
/*--------------------- End File ACNetworkBlock.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
