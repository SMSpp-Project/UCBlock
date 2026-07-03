/*--------------------------------------------------------------------------*/
/*--------------------- File ACNetworkBlock.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the ACNetworkBlock class.
 *
 * \author Quentin Jacquet \n
 *         EDF R&D OSIRIS \n
 *
 * \author Wim van Ackooij \n
 *         EDF R&D OSIRIS \n
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Quentin Jacquet, Wim van Ackooij
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>

#include <cmath>

#include <complex>

#include <iomanip>

#include <limits>

#include <map>

#include <sstream>

#include "ACNetworkBlock.h"

#include "LinearFunction.h"

#include "OneVarConstraint.h"

#include "FRealObjective.h"

#include "QuadFunction.h"

#include <Eigen/Sparse>

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

// register ACNetworkBlockSolution to the Solution factory
SMSpp_insert_in_factory_cpp_0( ACNetworkBlockSolution );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// register ACNetworkBlock::ACNetworkData to the NetworkData factory

using ACNetworkData = ACNetworkBlock::ACNetworkData;

SMSpp_insert_in_factory_cpp_0( ACNetworkData );

/*--------------------------------------------------------------------------*/
/*--------------------------- STATIC FUNCTIONS -----------------------------*/
/*--------------------------------------------------------------------------*/

// rounding function used to clip "ugly" admittance-matrix coefficients

static double round_sig( double value , int digits = 16 )
{
 if( value == 0.0 )
  return( 0.0 );

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
 DCNetworkData::deserialize( group );

 auto gbaseMVA = group.getAtt( "baseMVA" );
 if( gbaseMVA.isNull() )
  f_base_mva = 1.;
 else {
  std::string tmp_base;
  gbaseMVA.getValues( tmp_base );
  try { f_base_mva = std::stod( tmp_base ); }
  catch( ... ) { f_base_mva = 1.; }
  }

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

  ::deserialize( group , "LineChargingSusceptance" , f_number_lines ,
                 v_line_chargingsusceptance , true , true );

  ::deserialize( group , "NodeConductance" , f_number_nodes ,
                 v_node_conductance , true , true );

  ::deserialize( group , "NodeSusceptance" , f_number_nodes ,
                 v_node_susceptance , true , true );

  ::deserialize( group , "NodeMaxVoltage" , f_number_nodes ,
                 v_node_max_voltage , true , true );

  ::deserialize( group , "NodeMinVoltage" , f_number_nodes ,
                 v_node_min_voltage , true , true );

  // pad optional line-attribute vectors with neutral defaults so that
  // downstream code can always index them safely regardless of whether
  // they were present in the netCDF group:
  //   v_line_chargingsusceptance -> 0
  //   v_line_ratio               -> 1
  //   v_line_angle               -> 0

  if( v_line_chargingsusceptance.empty() )
   v_line_chargingsusceptance.assign( f_number_lines , 0.0 );
  if( v_line_ratio.empty() )
   v_line_ratio.assign( f_number_lines , 1.0 );
  if( v_line_angle.empty() )
   v_line_angle.assign( f_number_lines , 0.0 );

  // uncover the min and max reactive flow if there

  ::deserialize( group , "MinReactivePowerFlow" , f_number_lines ,
                 v_min_reac_power_flow , true , true );

  ::deserialize( group , "MaxReactivePowerFlow" , f_number_lines ,
                 v_max_reac_power_flow , true , true );

  // pre-fill v_DC_lines and v_HVDC_lines to override the logic in the
  // corresponding get_*() of DCNetworkData; this is because for AC
  // networks having 0 susceptance is not enough to declare that a line
  // is "DC" (i.e., not "HVDC"), but it must have 0 reactance and 0
  // resistance too

  v_DC_lines.clear();
  v_HVDC_lines.clear();
  if( v_line_susceptance.empty() ) {
   for( Index i = 0 ; i < f_number_lines ; ++i )
    if( ( ! v_line_reactance[ i ] ) && ( ! v_line_resistance[ i ] ) )
     v_HVDC_lines.push_back( i );
    else
     v_DC_lines.push_back( i );
   }
  else {
   for( Index i = 0 ; i < f_number_lines ; ++i )
    if( ( ! v_line_susceptance[ i ] ) && ( ! v_line_reactance[ i ] ) &&
        ( ! v_line_resistance[ i ] ) )
     v_HVDC_lines.push_back( i );
    else
     v_DC_lines.push_back( i );
   }

  f_number_HVDC_lines = v_HVDC_lines.size();
  v_HVDC_lines.shrink_to_fit();
  v_DC_lines.shrink_to_fit();
  }
 }  // end( ACNetworkData::deserialize )

/*--------------------------------------------------------------------------*/

void ACNetworkData::serialize( netCDF::NcGroup & group ) const
{
 // serialize the parent DCNetworkData first
 DCNetworkData::serialize( group );

 // "baseMVA" is (perhaps surprisingly) a string attribute, see deserialize()
 std::ostringstream base_mva;
 base_mva << std::setprecision( std::numeric_limits< double >::max_digits10 )
	  << f_base_mva;
 group.putAtt( "baseMVA" , base_mva.str() );

 if( f_number_nodes == 1 )
  return;

 auto NumberLines = group.getDim( "NumberLines" );
 auto NumberNodes = group.getDim( "NumberNodes" );

 ::serialize( group , "LineReactance" , netCDF::NcDouble() , NumberLines ,
              v_line_reactance );

 ::serialize( group , "LineResistance" , netCDF::NcDouble() , NumberLines ,
              v_line_resistance );

 ::serialize( group , "LineRatio" , netCDF::NcDouble() , NumberLines ,
              v_line_ratio );

 ::serialize( group , "LineRATEA" , netCDF::NcDouble() , NumberLines ,
              v_line_rate_A );

 ::serialize( group , "LineShiftAngle" , netCDF::NcDouble() , NumberLines ,
              v_line_angle );

 ::serialize( group , "LineMinAngle" , netCDF::NcDouble() , NumberLines ,
              v_line_min_angle );

 ::serialize( group , "LineMaxAngle" , netCDF::NcDouble() , NumberLines ,
              v_line_max_angle );

 ::serialize( group , "LineChargingSusceptance" , netCDF::NcDouble() ,
              NumberLines , v_line_chargingsusceptance );

 ::serialize( group , "NodeConductance" , netCDF::NcDouble() , NumberNodes ,
              v_node_conductance );

 ::serialize( group , "NodeSusceptance" , netCDF::NcDouble() , NumberNodes ,
              v_node_susceptance );

 ::serialize( group , "NodeMaxVoltage" , netCDF::NcDouble() , NumberNodes ,
              v_node_max_voltage );

 ::serialize( group , "NodeMinVoltage" , netCDF::NcDouble() , NumberNodes ,
              v_node_min_voltage );

 ::serialize( group , "MinReactivePowerFlow" , netCDF::NcDouble() ,
              NumberLines , v_min_reac_power_flow );

 ::serialize( group , "MaxReactivePowerFlow" , netCDF::NcDouble() ,
              NumberLines , v_max_reac_power_flow );

 }  // end( ACNetworkData::serialize )

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG
/*
std::vector< std::string > ACNetworkData::expected_dims( void ) const {
 return( DCNetworkData::expected_dims() );
 }
*/

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

std::vector< std::string > ACNetworkData::expected_vars( void ) const
{
 static const std::vector< std::string > ev =
 {
  "ReactivePowerDemand" , "NodeConductance" , "NodeSusceptance" ,
  "NodeVoltageMagnitude" , "NodeVoltageAngle" , "NodeMaxVoltage" ,
  "NodeMinVoltage" , "LineResistance" , "LineReactance" ,
  "LineChargingSusceptance" , "LineRatio" ,
  "LineRATEA" , "LineShiftAngle" , "LineMinAngle" , "LineMaxAngle" ,
  "MinReactivePowerFlow" , "MaxReactivePowerFlow"
 };

 auto ret = DCNetworkData::expected_vars();
 ret.insert( ret.end() , ev.begin() , ev.end() );

 return( ret );
 }

#endif

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS OF ACNetworkBlock ------------------------*/
/*--------------------------------------------------------------------------*/

ACNetworkBlock::~ACNetworkBlock()
{
 Constraint::clear( reactive_node_injection_bounds_const );
 Constraint::clear( v_voltage_bounds_const );
 Constraint::clear( v_angle_bounds_const );
 Constraint::clear( v_voltage_definition_const );
 Constraint::clear( v_thermal_limit );
 Constraint::clear( v_flow_dc );
 Constraint::clear( v_basic_bounds_const );
 Constraint::clear( v_socp_const );
 Constraint::clear( v_diag_const_1 );
 Constraint::clear( v_diag_const_2 );
 Constraint::clear( v_def_alpha_1 );
 Constraint::clear( v_def_alpha_2 );
 Constraint::clear( v_def_beta_1 );
 Constraint::clear( v_def_beta_2 );
 Constraint::clear( v_def_z_1 );
 Constraint::clear( v_def_z_2 );
 Constraint::clear( v_def_z_3 );
 Constraint::clear( v_def_z_4 );
 Constraint::clear( v_def_c_1 );
 Constraint::clear( v_def_c_2 );
 Constraint::clear( v_def_c_3 );
 Constraint::clear( v_def_c_4 );
 Constraint::clear( v_def_s_1 );
 Constraint::clear( v_def_s_2 );
 Constraint::clear( v_def_s_3 );
 Constraint::clear( v_def_s_4 );
 }  // end( ACNetworkBlock::~ACNetworkBlock )

/*--------------------------------------------------------------------------*/

void ACNetworkBlock::deserialize( const netCDF::NcGroup & group )
{
 // just have to call the method of the base class first: most of the
 // difference lies in the [AC]NetworkData, which is automatically
 // deserialize()-d there thanks to get_new_NetworkData()
 DCNetworkBlock::deserialize( group );

 // the reactive demand, if there, mirrors the (1D) active demand
 ::deserialize( group , "ReactiveDemand" , get_number_nodes() ,
                v_ReactiveDemand , true , true );

 }  // end( ACNetworkBlock::deserialize )

/*--------------------------------------------------------------------------*/

void ACNetworkBlock::serialize( netCDF::NcGroup & group ) const
{
 DCNetworkBlock::serialize( group );

 if( ! v_ReactiveDemand.empty() ) {
  // like for "ActiveDemand" in the base class, if no [AC]NetworkData is in
  // the group the number of nodes is recorded by an alternative dimension
  auto NumberNodes = group.getDim( "NumberNodes" );
  if( NumberNodes.isNull() )
   NumberNodes = group.getDim( "__NumberNodes__" );
  if( NumberNodes.isNull() )
   NumberNodes = group.addDim( "__NumberNodes__" , v_ReactiveDemand.size() );

  ::serialize( group , "ReactiveDemand" , netCDF::NcDouble() , NumberNodes ,
	       v_ReactiveDemand );
  }
 }  // end( ACNetworkBlock::serialize )

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG

std::vector< std::string > ACNetworkBlock::expected_vars( void ) const {
 auto ret = DCNetworkBlock::expected_vars();
 ret.push_back( "ReactiveDemand" );

 return( ret );
 }

#endif

/*--------------------------------------------------------------------------*/

void ACNetworkBlock::generate_abstract_variables( Configuration * stvv )
{
 if( variables_generated() )  // variables have already been generated
  return;                     // nothing to do

 DCNetworkBlock::generate_abstract_variables( stvv );

 // read the Configuration (if any)- - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( ! stvv ) && f_BlockConfig )
  stvv = f_BlockConfig->f_static_variables_Configuration;

 if( auto SCdd = dynamic_cast< SimpleConfiguration< int > * >( stvv ) )
  b_strongSOCP = ( SCdd->f_value > 0 );
 else if( auto SCvd = dynamic_cast<
                  SimpleConfiguration< std::vector< int > > * >( stvv ) ) {
  if( SCvd->f_value.size() > 0 )
   b_strongSOCP = ( SCvd->f_value[ 0 ] > 0 );
  }

 const auto number_nodes = get_number_nodes();
 const auto number_lines = get_number_lines();

 // the reactive node injection variables - - - - - - - - - - - - - - - - -
 v_reactive_node_injection.resize( number_nodes );
 for( Index node_id = 0 ; node_id < number_nodes ; ++node_id )
  v_reactive_node_injection[ node_id ].set_type( ColVariable::kContinuous );

 add_static_variable( v_reactive_node_injection , "reactive_s_network" );

 // complex power flow (real and imaginary part for both directions) - - - -
 /* S denotes the AC power for each line: it is a 2-dimensional vector
  * (real and imaginary part). */
 v_power_flow.resize( 2 * number_lines );
 for( Index line_id = 0 ; line_id < 2 * number_lines ; ++line_id )
  v_power_flow[ line_id ].set_type( ColVariable::kContinuous );
 add_static_variable( v_power_flow , "v_power_flow_real" );

 v_reactive_power_flow.resize( 2 * number_lines );
 for( Index line_id = 0 ; line_id < 2 * number_lines ; ++line_id )
  v_reactive_power_flow[ line_id ].set_type( ColVariable::kContinuous );
 add_static_variable( v_reactive_power_flow , "v_reactive_power_flow" );

 // we do not define the reverse value because the sum is symmetric and the
 // diff is anti-symmetric
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

 // if so desired, add the variables for the strengthened SOCP relaxation
 if( b_strongSOCP )
  generate_strengthened_variables();
 }  // end( ACNetworkBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void ACNetworkBlock::generate_strengthened_variables( void )
{
 // uncover some size information
 const auto number_nodes = get_number_nodes();
 auto & DC_lines = f_NetworkData->get_DC_lines();
 int nb_dc_lines = DC_lines.size();
 int i_line;

 // generate auxiliary variables - - - - - - - - - - - - - - - - - - - - - -

 v_voltage.resize( number_nodes );
 for( Index node_id = 0 ; node_id < number_nodes ; ++node_id )
  v_voltage[ node_id ].set_type( ColVariable::kContinuous );
 add_static_variable( v_voltage , "v_voltage" );

 v_theta.resize( number_nodes );
 for( Index node_id = 0 ; node_id < number_nodes ; ++node_id )
  v_theta[ node_id ].set_type( ColVariable::kContinuous );
 add_static_variable( v_theta , "v_theta" );

 v_alpha.resize( nb_dc_lines );
 i_line = 0;
 for( auto & line_id : DC_lines ) {
  v_alpha[ i_line ].set_type( ColVariable::kContinuous );
  ++i_line;
  }
 add_static_variable( v_alpha , "v_alpha" );

 v_beta.resize( nb_dc_lines );
 i_line = 0;
 for( auto & line_id : DC_lines ) {
  v_beta[ i_line ].set_type( ColVariable::kContinuous );
  ++i_line;
  }
 add_static_variable( v_beta , "v_beta" );

 v_z.resize( nb_dc_lines );
 i_line = 0;
 for( auto & line_id : DC_lines ) {
  v_z[ i_line ].set_type( ColVariable::kContinuous );
  ++i_line;
  }
 add_static_variable( v_z , "v_z" );
 }  // end( ACNetworkBlock::generate_strengthened_variables )

/*--------------------------------------------------------------------------*/

void ACNetworkBlock::generate_abstract_constraints( Configuration * stcc )
{
 if( constraints_generated() )  // constraints have already been generated
  return;                       // nothing to do

 // read the Configuration (if any)- - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 f_C_v_scal = 1;
 double f_ACvS = 0.0;  // slack for AC_voltage_definition_const
 double f_scale = 1.0; // scaling for AC_voltage_definition_const equations
 int f_digits = 16;

 if( ( ! stcc ) && f_BlockConfig )
  stcc = f_BlockConfig->f_static_constraints_Configuration;

 if( auto SCdd = dynamic_cast< SimpleConfiguration< double > * >( stcc ) )
  f_C_v_scal = SCdd->f_value;
 else if( auto SCpd = dynamic_cast<
              SimpleConfiguration< std::pair< double , double > > * >( stcc )
          ) {
  f_C_v_scal = SCpd->f_value.first;
  f_ACvS = SCpd->f_value.second;
  }
 else if( auto SCvd = dynamic_cast<
              SimpleConfiguration< std::vector< double > > * >( stcc ) ) {
  if( SCvd->f_value.size() > 0 )
   f_C_v_scal = SCvd->f_value[ 0 ];
  if( SCvd->f_value.size() > 1 )
   f_ACvS = SCvd->f_value[ 1 ];
  if( SCvd->f_value.size() > 2 )
   f_scale = SCvd->f_value[ 2 ];
  if( SCvd->f_value.size() > 3 )
   f_digits = std::max( 0 , std::min( 16 , int( SCvd->f_value[ 3 ] ) ) );
  }

 const auto number_nodes = get_number_nodes();

 if( number_nodes <= 1 )
  return;

 const auto number_lines = get_number_lines();

 if( number_lines <= 0 )
  throw( std::logic_error( "ACNetworkBlock::generate_abstract_constraints: "
                           "number of lines of DCNetworkBlock is not set" )
         );

 // do not call DCNetworkBlock::generate_abstract_constraints( stcc ) since
 // the class entirely redefines its constraints, but directly call
 // generate_bound_constraints() to have the bound constraints generated
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 generate_bound_constraints();

 // node injection bound constraints- - - - - - - - - - - - - - - - - - - -
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

 // now start processing lines- - - - - - - - - - - - - - - - - - - - - - -
 const auto & start_line = f_NetworkData->get_start_line();
 const auto & end_line = f_NetworkData->get_end_line();

 auto & DC_lines = f_NetworkData->get_DC_lines();
 int nb_dc_lines = DC_lines.size();

 auto & HVDC_lines = f_NetworkData->get_HVDC_lines();
 int nb_hvdc_lines = f_NetworkData->get_number_HVDC_lines();

 double base_mva = ND()->get_baseMVA();

 /* Scaling of v_power_flow and v_reactive_power_flow to improve numerical
  * stability; effectively we are swapping out v_power_flow for
  * v_power_flow_tilde with v_power_flow_tilde = C * v_power_flow, and
  * likewise for v_reactive_power_flow. */

 // ----- Voltage bounds- - - - - - - - - - - - - - - - - - - - - - - - - -
 /* Voltage magnitude is bounded between min_voltage and max_voltage. To
  * this aim,
  *    W_{n,n} = (V_n).(V_n)^H = |V_n|^2
  * so we impose the bounds directly on W_{n,n}. */
 v_voltage_bounds_const.resize( number_nodes );
 const auto & min_voltage = ND()->get_node_min_voltage();
 const auto & max_voltage = ND()->get_node_max_voltage();
 for( Index n = 0 ; n < number_nodes ; ++n ) {
  v_voltage_bounds_const[ n ].set_lhs( pow( f_C_v_scal * min_voltage[ n ] ,
                                            2 ) );
  v_voltage_bounds_const[ n ].set_rhs( pow( f_C_v_scal * max_voltage[ n ] ,
                                            2 ) );
  v_voltage_bounds_const[ n ].set_variable( & v_sqrd_voltages[ n ] );
  }
 add_static_constraint( v_voltage_bounds_const , "AC_voltage_bounds_limit" );

 // ----- Angle bounds- - - - - - - - - - - - - - - - - - - - - - - - - - -
 /* Phase Angle Difference (PAD) constraints for each line = (start,end):
  *   min_angle_{line} <= angle(V_{end}) - angle(V_{start}) <= max_angle_{line}
  * which can be reformulated using
  *   tan(min_angle_{line}) Real(W_{start,end}) <= Imag(W_{start,end})
  *                                             <= max_angle_{line}
  *                                                Real(W_{start,end}).
  *
  * Furthermore we induce basic bounds on v_diff_product_voltages and
  * v_sum_product_voltages by leveraging the minimum voltages, maximum
  * voltages and angle bounds. */
 int i_line = 0;
 v_angle_bounds_const.resize( MAFRC_ext()[ 2 ][ nb_dc_lines ] );
 v_basic_bounds_const.resize( MAFRC_ext()[ 2 ][ nb_dc_lines ] );

 const auto & min_angle = ND()->get_line_min_angle();
 const auto & max_angle = ND()->get_line_max_angle();

 for( auto & line_id : DC_lines ) {
  // all angles are typically input as degrees, but we need radians

  double phi_min = PI * min_angle[ line_id ] / 180.;
  double phi_max = PI * max_angle[ line_id ] / 180.;
  double delta_phi = phi_max - phi_min;
  // assuming phi_min <= phi_max evidently

  // classic angle-based bounds on c_{n,n'} and s_{n,n'}
  // tan( phi_min ) c_{n,n'} <= s_{n,n'}
  auto lfunc_1 = new LinearFunction();
  lfunc_1->add_variable( & v_diff_product_voltages[ line_id ] , 1.0 );
  lfunc_1->add_variable( & v_sum_product_voltages[ line_id ] ,
                         -tan( phi_min ) );
  v_angle_bounds_const[ 0 ][ i_line ].set_lhs( 0.0 );
  v_angle_bounds_const[ 0 ][ i_line ].set_rhs( Inf< double >() );
  v_angle_bounds_const[ 0 ][ i_line ].set_function( lfunc_1 );

  // second half
  // s_{n,n'} <= tan( phi_max ) c_{n,n'}
  auto lfunc_2 = new LinearFunction();
  lfunc_2->add_variable( & v_diff_product_voltages[ line_id ] , 1.0 );
  lfunc_2->add_variable( & v_sum_product_voltages[ line_id ] ,
                         -tan( phi_max ) );
  v_angle_bounds_const[ 1 ][ i_line ].set_lhs( -Inf< double >() );
  v_angle_bounds_const[ 1 ][ i_line ].set_rhs( 0.0 );
  v_angle_bounds_const[ 1 ][ i_line ].set_function( lfunc_2 );

  // bounds on v_sum_product_voltages- - - - - - - - - - - - - - - - - - -
  // v_sum_product_voltages = c_{n,n'} = v_n v_n' cos( theta_n - theta_n' );
  // from this relation and the allowed angle bounds (directly bounding
  // theta_n - theta_n') we can deduce proper bounds on these variables.
  v_basic_bounds_const[ 0 ][ i_line ].set_lhs(
   std::min( cos( std::abs( phi_min ) ) , cos( std::abs( phi_max ) ) ) *
   min_voltage[ start_line[ line_id ] ] *
   min_voltage[ end_line[ line_id ] ] *
   pow( f_C_v_scal , 2 ) );
  v_basic_bounds_const[ 0 ][ i_line ].set_rhs(
   max_voltage[ start_line[ line_id ] ] *
   max_voltage[ end_line[ line_id ] ] *
   pow( f_C_v_scal , 2 ) );
  v_basic_bounds_const[ 0 ][ i_line ].set_variable(
   & v_sum_product_voltages[ line_id ] );

  // bounds on v_diff_product_voltages - - - - - - - - - - - - - - - - - -
  // v_diff_product_voltages = s_{n,n'} = v_n v_n' sin( theta_n - theta_n' );
  // from this relation and the allowed angle bounds we can deduce proper
  // bounds on these variables.
  double s_sin = sin( delta_phi );
  v_basic_bounds_const[ 1 ][ i_line ].set_lhs(
   -1.0 * s_sin * max_voltage[ start_line[ line_id ] ] *
   max_voltage[ end_line[ line_id ] ] * pow( f_C_v_scal , 2 ) );
  v_basic_bounds_const[ 1 ][ i_line ].set_rhs(
   s_sin * max_voltage[ start_line[ line_id ] ] *
   max_voltage[ end_line[ line_id ] ] * pow( f_C_v_scal , 2 ) );
  v_basic_bounds_const[ 1 ][ i_line ].set_variable(
   & v_diff_product_voltages[ line_id ] );

  ++i_line;
  }

 if( i_line > 0 ) {
  add_static_constraint( v_angle_bounds_const , "AC_angle_bounds_limit" );
  add_static_constraint( v_basic_bounds_const , "AC_elem_bounds" );
  }

 auto * fnet = static_cast< ACNetworkData * >( f_NetworkData );

 // bounds on reactive flow in HVDC lines
 // HVDC lines have direct bounds both on Active and Reactive power (if given)
 if( fnet->has_reactive_bounds() ) {
  v_reactive_flow_bounds.resize( 2 * nb_hvdc_lines );
  int i_hvdc_line = 0;
  for( auto & line_id : HVDC_lines ) {
   v_reactive_flow_bounds[ i_hvdc_line ].set_lhs(
    f_C_v_scal * fnet->get_min_reac_power_flow( line_id ) );
   v_reactive_flow_bounds[ i_hvdc_line ].set_rhs(
    f_C_v_scal * fnet->get_max_reac_power_flow( line_id ) );
   v_reactive_flow_bounds[ i_hvdc_line ].set_variable(
    & v_reactive_power_flow[ line_id ] );
   // bound on the "to" flow
   v_reactive_flow_bounds[ nb_hvdc_lines + i_hvdc_line ].set_lhs(
    f_C_v_scal * fnet->get_min_reac_power_flow( line_id ) );
   v_reactive_flow_bounds[ nb_hvdc_lines + i_hvdc_line ].set_rhs(
    f_C_v_scal * fnet->get_max_reac_power_flow( line_id ) );
   v_reactive_flow_bounds[ nb_hvdc_lines + i_hvdc_line ].set_variable(
    & v_reactive_power_flow[ number_lines + line_id ] );

   ++i_hvdc_line;
   }

  add_static_constraint( v_reactive_flow_bounds , "Reactive_Flow_Bounds" );
  }

 // ----- Active and reactive power conservation- - - - - - - - - - - - - -
 // shunt admittance
 SpCVec Ys = SpCVec( number_nodes );
 for( Index n = 0 ; n < number_nodes ; ++n ) {
  double Gs = ND()->get_node_conductance().at( n ) / base_mva;
  double Bs = ND()->get_node_susceptance().at( n ) / base_mva;
  Ys.insert( n ) = Gs + 1i * Bs;
  }

 /* The power conservation constraint can be written as
  *
  *   Supply - Demand = sum_{line_id in L cup L^R} power_flow[ line_id ]
  *
  * The complex matrix product <M,W>_F is then decomposed into a real part
  * and an imaginary part. */
 v_power_flow_injection_const.resize( 2 * number_nodes );

 // real part of the power flow conservation- - - - - - - - - - - - - - - -
 for( Index p = 0 ; p < number_nodes ; ++p ) {
  auto lfunc = new LinearFunction();
  lfunc->add_variable( & v_node_injection[ 0 ][ p ] , -1.0 * f_C_v_scal );
  // why 0 and not the time step?
  lfunc->add_variable( & v_sqrd_voltages[ p ] ,
                       -Ys.coeff( p ).real() / f_C_v_scal );

  for( Index line_id = 0 ; line_id < number_lines ; ++line_id ) {
   Index i = start_line[ line_id ];
   Index j = end_line[ line_id ];

   if( i == p )
    lfunc->add_variable( & v_power_flow[ line_id ] , 1.0 );
   if( j == p )
    lfunc->add_variable( & v_power_flow[ number_lines + line_id ] , 1.0 );
   }
  v_power_flow_injection_const[ p ].set_both(
   -v_ActiveDemand[ p ] * f_C_v_scal );
  v_power_flow_injection_const[ p ].set_function( lfunc );
  }

 // imaginary part of the power flow conservation - - - - - - - - - - - - -
 for( Index p = 0 ; p < number_nodes ; ++p ) {
  auto lfunc = new LinearFunction();
  lfunc->add_variable( & v_reactive_node_injection[ p ] ,
                       -1.0 * f_C_v_scal );
  lfunc->add_variable( & v_sqrd_voltages[ p ] ,
                       -Ys.coeff( p ).imag() / f_C_v_scal );

  for( Index line_id = 0 ; line_id < number_lines ; ++line_id ) {
   Index i = start_line[ line_id ];
   Index j = end_line[ line_id ];

   if( i == p )
    lfunc->add_variable( & v_reactive_power_flow[ line_id ] , 1.0 );
   if( j == p )
    lfunc->add_variable(
     & v_reactive_power_flow[ number_lines + line_id ] , 1.0 );
   }
  v_power_flow_injection_const[ number_nodes + p ].set_both(
   -v_ReactiveDemand[ p ] * f_C_v_scal );
  v_power_flow_injection_const[ number_nodes + p ].set_function( lfunc );
  }

 add_static_constraint( v_power_flow_injection_const ,
                        "AC_power_flow_injection" );

 // ----- HVDC flow links - - - - - - - - - - - - - - - - - - - - - - - - -
 /* Since the lines have been duplicated, for DC ones we need to add the
  * link between the two versions: the specific nature of HVDC lines is
  * that the flows are opposite of each other. Observe that this is not
  * true of the "imaginary" (reactive) part of the flow, however. */
 v_flow_dc.resize( nb_hvdc_lines );
 int i_hvdc_line = 0;
 for( auto & line_id : HVDC_lines ) {
  auto lfunc = new LinearFunction();
  lfunc->add_variable( & v_power_flow[ line_id ] , 1.0 );
  lfunc->add_variable( & v_power_flow[ number_lines + line_id ] , 1.0 );
  v_flow_dc[ i_hvdc_line ].set_both( 0.0 );
  v_flow_dc[ i_hvdc_line ].set_function( lfunc );
  ++i_hvdc_line;
  }

 add_static_constraint( v_flow_dc , "HVDC_flow_links" );

 // ----- Definition of complex power flow- - - - - - - - - - - - - - - - -
 /* The complex power flow S_{line} for each line = (start,end) is
 *   S_{start,end} = Yff_{start,end} W_{start,start}
 *                 + Yft_{start,end} W_{start,end}
 *   S_{end,start} = Ytt_{start,end} W_{end,end}
 *                 + Ytf_{start,end} W_{end,start}
 * which is then split into two constraints, one for the real and one for
 * the imaginary part.
 *
 * Note that these equations make it easy to perform an elementary check
 * on flow feasibility: from the bounds on W (lower and upper) and the
 * sign of the coefficients we can deduce that the flow on each line
 * should be at least some value and at most some other value. We can
 * then cross-check this against the thermal limits and warn the user.
 *
 * Indeed with I+ = { i : c_i > 0 } and I- = { i : c_i < 0 }, we have
 * v_flow = sum_{I+} c_i w_i + sum_{I-} c_i w_i and therefore
 *   v_flow >= sum_{I+} c_i underline{w}_i + sum_{I-} c_i overline{w}_i
 *   v_flow <= sum_{I+} c_i overline{w}_i + sum_{I-} c_i underline{w}_i */

 /* Shortcuts to recover the mathematical notation:
  *
  *   Yff_{ab} = [ 1 / ( r_{ab} + i x_{ab} ) + i * 0.5 * b_{ab} ] / tau_{ab}^2
  *   Ytt_{ab} =   1 / ( r_{ab} + i x_{ab} ) + i * 0.5 * b_{ab}
  *   Yft_{ab} = -[ 1 / ( r_{ab} + i x_{ab} ) ] / ( tau_{ab} e^{-i nu_{ab}} )
  *   Ytf_{ab} = -[ 1 / ( r_{ab} + i x_{ab} ) ] / ( tau_{ab} e^{+i nu_{ab}} )
  *
  * Since these coefficients (both real and imaginary parts) are often
  * not "nice" numbers, we have incorporated two options:
  *   a) rounding, using the f_digits value (ensuring 25th-decimal junk
  *      does not interfere);
  *   b) add slack in the equations through f_ACvS.
  *
  * Both parameters seem delicate to fine-tune and perhaps ought to be
  * reconsidered at some stage. */
 auto * f_net = static_cast< ACNetworkData * >( f_NetworkData );

 auto r = [ f_net ]( int line_id ) {
  return( f_net->get_line_resistance().at( line_id ) );
  };
 auto x = [ f_net ]( int line_id ) {
  return( f_net->get_line_reactance().at( line_id ) );
  };
 auto b = [ f_net ]( int line_id ) {
  return( f_net->get_line_chargingsusceptance().at( line_id ) );
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

 // for the check
 const auto & v_l_names = ND()->get_line_names();
 const auto & rate_A = ND()->get_line_rate_A();

 v_voltage_definition_const.resize( MAFRC_ext()[ 2 ][ 2 * nb_dc_lines ] );
 i_line = 0;
 const auto splitted_lines = ND()->get_direct_and_reverse_AClines();
 for( Index p = 0 ; p < number_nodes ; ++p ) {
  // 1) direct lines
  for( auto & line_id : splitted_lines[ p ].first ) {
   // 1.0) verify that the line ratio is not zero (would yield NaN coeffs)
   if( ! std::isfinite( Yff( line_id ).real() ) )
    throw( std::logic_error( "Non-finite coefficient, possibly line ratio "
                             "is zero for an AC line" ) );

   // values for the basic bound check
   double v_flow_lower = 0.0;
   double v_flow_upper = 0.0;
   double phi_min = PI * min_angle[ line_id ] / 180.;
   double phi_max = PI * max_angle[ line_id ] / 180.;
   double delta_phi = phi_max - phi_min;
   // assuming phi_min <= phi_max evidently
   double c_cos = std::min( cos( std::abs( phi_min ) ) ,
                            cos( std::abs( phi_max ) ) );
   double c_sin = sin( delta_phi );

   // 1.1) real part
   auto lfunc_1 = new LinearFunction();
   lfunc_1->add_variable( & v_sqrd_voltages[ p ] ,
                          round_sig( Yff( line_id ).real() * f_scale ,
                                     f_digits ) );

   // update bounds
   v_flow_lower += std::max( Yff( line_id ).real() , 0.0 ) *
    std::pow( min_voltage[ p ] , 2.0 ) +
    std::min( Yff( line_id ).real() , 0.0 ) *
    std::pow( max_voltage[ p ] , 2.0 );
   v_flow_upper += std::max( Yff( line_id ).real() , 0.0 ) *
    std::pow( max_voltage[ p ] , 2.0 ) +
    std::min( Yff( line_id ).real() , 0.0 ) *
    std::pow( min_voltage[ p ] , 2.0 );

   lfunc_1->add_variable( & v_sum_product_voltages[ line_id ] ,
                          round_sig( Yft( line_id ).real() * f_scale ,
                                     f_digits ) );

   v_flow_lower += std::max( Yft( line_id ).real() , 0.0 ) * c_cos *
    min_voltage[ p ] * min_voltage[ end_line[ line_id ] ] +
    std::min( Yft( line_id ).real() , 0.0 ) * max_voltage[ p ] *
    max_voltage[ end_line[ line_id ] ];
   v_flow_upper += std::max( Yft( line_id ).real() , 0.0 ) *
    max_voltage[ p ] * max_voltage[ end_line[ line_id ] ] +
    std::min( Yft( line_id ).real() , 0.0 ) * c_cos *
    min_voltage[ p ] * min_voltage[ end_line[ line_id ] ];

   lfunc_1->add_variable( & v_diff_product_voltages[ line_id ] ,
                          round_sig( Yft( line_id ).imag() * f_scale ,
                                     f_digits ) );

   // since the bounds are symmetric, it suffices to compute this one
   // value, moreover the sum greatly simplifies:
   double w_bound = c_sin * max_voltage[ p ] *
                    max_voltage[ end_line[ line_id ] ];
   v_flow_lower += std::abs( Yft( line_id ).imag() ) * -1.0 * w_bound;
   v_flow_upper += std::abs( Yft( line_id ).imag() ) * w_bound;

   // sanity check vs. min/max power flow and thermal limit
   auto max_p = get_max_power_flow( line_id );
   auto min_p = get_min_power_flow( line_id );
   if( ( v_flow_lower > max_p ) || ( v_flow_upper < min_p ) ) {
    std::cout << " The power line with index = " << line_id << " and name "
              << v_l_names[ line_id ]
              << " has induced bounds from the AC equations that are [ "
              << v_flow_lower << ", " << v_flow_upper << "]"
              << " and imposed bounds [ " << min_p << " , " << max_p
              << " ]" << std::endl;
    }

   if( ( v_flow_lower > 0 ) || ( v_flow_upper < 0 ) ) {
    double min_therm = std::min( std::pow( v_flow_lower , 2.0 ) ,
                                 std::pow( v_flow_upper , 2.0 ) );
    if( min_therm > rate_A[ line_id ] )
     std::cout << " The power line with index = " << line_id
               << " and name " << v_l_names[ line_id ]
               << " has induced bounds from the AC equations that yield a "
               << "minimal thermal limit of " << min_therm
               << " but this exceeds the given limit " << rate_A[ line_id ]
               << std::endl;
    }

   lfunc_1->add_variable( & v_power_flow[ line_id ] ,
                          -1.0 * f_C_v_scal * f_scale );
   v_voltage_definition_const[ 0 ][ i_line ].set_lhs( -f_ACvS );
   v_voltage_definition_const[ 0 ][ i_line ].set_rhs( f_ACvS );
   v_voltage_definition_const[ 0 ][ i_line ].set_function( lfunc_1 );

   // 1.2) imaginary part
   auto lfunc_2 = new LinearFunction();
   lfunc_2->add_variable( & v_sqrd_voltages[ p ] ,
                          round_sig( -Yff( line_id ).imag() * f_scale ,
                                     f_digits ) );
   lfunc_2->add_variable( & v_sum_product_voltages[ line_id ] ,
                          round_sig( -Yft( line_id ).imag() * f_scale ,
                                     f_digits ) );
   lfunc_2->add_variable( & v_diff_product_voltages[ line_id ] ,
                          round_sig( Yft( line_id ).real() * f_scale ,
                                     f_digits ) );
   lfunc_2->add_variable( & v_reactive_power_flow[ line_id ] ,
                          -1.0 * f_C_v_scal * f_scale );
   v_voltage_definition_const[ 1 ][ i_line ].set_lhs( -f_ACvS );
   v_voltage_definition_const[ 1 ][ i_line ].set_rhs( f_ACvS );
   v_voltage_definition_const[ 1 ][ i_line ].set_function( lfunc_2 );

   ++i_line;
   }

  // 2) reverse lines
  for( auto & line_id : splitted_lines[ p ].second ) {
   // 2.1) real part
   auto lfunc_1 = new LinearFunction();
   lfunc_1->add_variable( & v_sqrd_voltages[ p ] ,
                          round_sig( Ytt( line_id ).real() * f_scale ,
                                     f_digits ) );
   lfunc_1->add_variable( & v_sum_product_voltages[ line_id ] ,
                          round_sig( Ytf( line_id ).real() * f_scale ,
                                     f_digits ) );
   lfunc_1->add_variable( & v_diff_product_voltages[ line_id ] ,
                          round_sig( -Ytf( line_id ).imag() * f_scale ,
                                     f_digits ) );
   // be careful: diff is antisymmetric
   lfunc_1->add_variable( & v_power_flow[ number_lines + line_id ] ,
                          -1.0 * f_C_v_scal * f_scale );
   v_voltage_definition_const[ 0 ][ i_line ].set_lhs( -f_ACvS );
   v_voltage_definition_const[ 0 ][ i_line ].set_rhs( f_ACvS );
   v_voltage_definition_const[ 0 ][ i_line ].set_function( lfunc_1 );

   // 2.2) imaginary part
   auto lfunc_2 = new LinearFunction();
   lfunc_2->add_variable( & v_sqrd_voltages[ p ] ,
                          round_sig( -Ytt( line_id ).imag() * f_scale ,
                                     f_digits ) );
   lfunc_2->add_variable( & v_sum_product_voltages[ line_id ] ,
                          round_sig( -Ytf( line_id ).imag() * f_scale ,
                                     f_digits ) );
   lfunc_2->add_variable( & v_diff_product_voltages[ line_id ] ,
                          round_sig( -Ytf( line_id ).real() * f_scale ,
                                     f_digits ) );
   // be careful: diff is antisymmetric
   lfunc_2->add_variable( & v_reactive_power_flow[ number_lines + line_id ] ,
                          -1.0 * f_C_v_scal * f_scale );
   v_voltage_definition_const[ 1 ][ i_line ].set_lhs( -f_ACvS );
   v_voltage_definition_const[ 1 ][ i_line ].set_rhs( f_ACvS );
   v_voltage_definition_const[ 1 ][ i_line ].set_function( lfunc_2 );

   ++i_line;
   }
  }  // at this stage i_line should be nb_dc_lines * 2

 if( i_line > 0 )
  add_static_constraint( v_voltage_definition_const ,
                         "AC_voltage_definition_const" );

 // ----- Thermal limit on lines- - - - - - - - - - - - - - - - - - - - - -
 /* We impose |S_{line}| <= rateA_{line}, a thermal limitation. To take
  * this into account, we use a DQuadFunction:
  *   Real( S_{line} )^2 + Imag( S_{line} )^2 <= rateA_{line}^2 */
 v_thermal_limit.resize( 2 * nb_dc_lines );
 i_line = 0;
 for( auto & line_id : DC_lines ) {
  auto qfunc_1 = new DQuadFunction();
  qfunc_1->add_variable( & v_power_flow[ line_id ] , 0.0 , 1.0 );
  qfunc_1->add_variable( & v_reactive_power_flow[ line_id ] , 0.0 , 1.0 );
  v_thermal_limit[ i_line ].set_lhs( -Inf< double >() );
  v_thermal_limit[ i_line ].set_rhs(
   pow( f_C_v_scal * rate_A[ line_id ] / base_mva , 2 ) );
  v_thermal_limit[ i_line ].set_function( qfunc_1 );

  auto qfunc_2 = new DQuadFunction();
  qfunc_2->add_variable( & v_power_flow[ number_lines + line_id ] ,
                         0.0 , 1.0 );
  qfunc_2->add_variable( & v_reactive_power_flow[ number_lines + line_id ] ,
                         0.0 , 1.0 );
  v_thermal_limit[ nb_dc_lines + i_line ].set_lhs( -Inf< double >() );
  v_thermal_limit[ nb_dc_lines + i_line ].set_rhs(
   pow( f_C_v_scal * rate_A[ line_id ] / base_mva , 2 ) );
  v_thermal_limit[ nb_dc_lines + i_line ].set_function( qfunc_2 );

  ++i_line;
  }
 add_static_constraint( v_thermal_limit , "AC_thermal_limit_const" );

 /* The above equations are always valid. However, in order to have a
  * convex model we use the SOCP relaxation of the following non-convex
  * quadratic relations:
  *     c_{n,n'}  =  c_{n',n}
  *     s_{n,n'}  = -s_{n',n}
  *     c_{n,n'}^2 + s_{n,n'}^2 = c_{n,n} c_{n',n'}
  *
  * Stronger relaxations using Semi-Definite programming exist, but are
  * not implemented yet. Should this become so, the call to the following
  * function can be switched upon. */
 generate_SOCP_relaxation();

 /* Moreover the SOCP relaxation can be made much stronger following the
  * work by Coffrin et al. This can be done by adding multiple McCormick
  * inequalities; it is optional and can be triggered from the
  * BlockConfig file. */
 if( b_strongSOCP )
  strengthen_SOCP_relaxation();

 set_constraints_generated();  // signal all done

 }  // end( ACNetworkBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

/* Links between the generic variables
 *   v_sum_product_voltages, v_diff_product_voltages, v_sqrd_voltages
 * using a SOCP relaxation. */

void ACNetworkBlock::generate_SOCP_relaxation( void )
{
 const auto & start_line = f_NetworkData->get_start_line();
 const auto & end_line = f_NetworkData->get_end_line();

 auto & DC_lines = f_NetworkData->get_DC_lines();
 int nb_dc_lines = DC_lines.size();

 /* Voltage relaxation matrix W.
  *
  * We aim to impose W = V.V^H, where V is the vector of voltages for each
  * bus/node. This nonlinear constraint is replaced by a SOCP relaxation.
  * As W is a complex matrix, we define in the optimization model two
  * matrices:
  *   - W_voltage ( 0 <= i < n , 0 <= j < n ) for the real part
  *   - W_voltage ( n <= i < 2n , n <= j < 2n ) for the imaginary part
  *
  * Since we cannot take into account the true constraint W = V.V^H, we
  * replace it by the SOCP relaxation
  *
  *   | W_{ab} |^2 <= W_{aa} W_{bb}
  *
  * Since we are in complex algebra, auxiliary variables are needed to
  * write the SOCP constraints (this might be simplifiable). */

 v_socp_const.resize( nb_dc_lines );
 int i_line = 0;
 for( auto & line_id : DC_lines ) {
  Index p = start_line[ line_id ];
  Index n = end_line[ line_id ];
  auto qfunc = new QuadFunction();
  qfunc->add_variable( & v_sqrd_voltages[ p ] , 0.0 , 0.0 );
  qfunc->add_variable( & v_sqrd_voltages[ n ] , 0.0 , 0.0 );
  qfunc->add_nd_term( & v_sqrd_voltages[ p ] , & v_sqrd_voltages[ n ] ,
                      -1.0 );
  qfunc->add_variable( & v_sum_product_voltages[ line_id ] , 0.0 , 1.0 );
  qfunc->add_variable( & v_diff_product_voltages[ line_id ] , 0.0 , 1.0 );
  v_socp_const[ i_line ].set_lhs( -Inf< double >() );
  v_socp_const[ i_line ].set_rhs( 0.0 );
  v_socp_const[ i_line ].set_function( qfunc );

  ++i_line;
  }

 if( i_line > 0 )
  add_static_constraint( v_socp_const , "AC_socp_const" );

 }  // end( ACNetworkBlock::generate_SOCP_relaxation )

/*--------------------------------------------------------------------------*/

void ACNetworkBlock::strengthen_SOCP_relaxation( void )
{
 auto * f_net = static_cast< ACNetworkData * >( f_NetworkData );
 const auto & start_line = f_net->get_start_line();
 const auto & end_line = f_net->get_end_line();
 const auto number_nodes = get_number_nodes();

 const auto & min_voltage = f_net->get_node_min_voltage();
 const auto & max_voltage = f_net->get_node_max_voltage();

 const auto & v_line_min_angle = f_net->get_line_min_angle();
 const auto & v_line_max_angle = f_net->get_line_max_angle();

 auto & DC_lines = f_net->get_DC_lines();
 int nb_dc_lines = DC_lines.size();
 int i_line;

 // simple bounds on v_voltage- - - - - - - - - - - - - - - - - - - - - - -
 // these are evident from what the variables represent
 v_volt_bounds.resize( number_nodes );
 for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {
  v_volt_bounds[ node_id ].set_lhs( min_voltage[ node_id ] * f_C_v_scal );
  v_volt_bounds[ node_id ].set_rhs( max_voltage[ node_id ] * f_C_v_scal );
  v_volt_bounds[ node_id ].set_variable( & v_voltage[ node_id ] );
  }
 add_static_constraint( v_volt_bounds , "v_volt_bounds" );

 // simple bounds on v_theta- - - - - - - - - - - - - - - - - - - - - - - -
 v_theta_bounds.resize( nb_dc_lines );
 i_line = 0;
 for( auto & line_id : DC_lines ) {
  Index p = start_line[ line_id ];
  Index n = end_line[ line_id ];

  auto lfunc = new LinearFunction();
  lfunc->add_variable( & v_theta[ p ] , 1.0 );
  lfunc->add_variable( & v_theta[ n ] , -1.0 );
  v_theta_bounds[ i_line ].set_lhs( PI * v_line_min_angle[ line_id ] /
                                    180.0 );
  v_theta_bounds[ i_line ].set_rhs( PI * v_line_max_angle[ line_id ] /
                                    180.0 );
  v_theta_bounds[ i_line ].set_function( lfunc );
  ++i_line;
  }
 add_static_constraint( v_theta_bounds , "v_theta_bounds" );

 // bounds on v_alpha - - - - - - - - - - - - - - - - - - - - - - - - - - -
 v_alpha_bounds.resize( nb_dc_lines );
 i_line = 0;
 for( auto & line_id : DC_lines ) {
  v_alpha_bounds[ i_line ].set_lhs( -1.0 );
  v_alpha_bounds[ i_line ].set_rhs( 1.0 );
  v_alpha_bounds[ i_line ].set_variable( & v_alpha[ i_line ] );
  ++i_line;
  }
 add_static_constraint( v_alpha_bounds , "v_alpha_bounds" );

 // bounds on v_beta- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 v_beta_bounds.resize( nb_dc_lines );
 i_line = 0;
 for( auto & line_id : DC_lines ) {
  v_beta_bounds[ i_line ].set_lhs( -1.0 );
  v_beta_bounds[ i_line ].set_rhs( 1.0 );
  v_beta_bounds[ i_line ].set_variable( & v_beta[ i_line ] );
  ++i_line;
  }
 add_static_constraint( v_beta_bounds , "v_beta_bounds" );

 // generate auxiliary constraints- - - - - - - - - - - - - - - - - - - - -

 /* McCormick envelope of the square term V_n^2, which directly yields
  *   c_{n,n} >= v_n^2
  *   c_{n,n} <= ( overline{v}_n + underline{v}_n ) v_n
  *              - overline{v}_n underline{v}_n
  *
  * This is eq. (21a) combined with (T-CONV) in Coffrin (2016). */
 v_diag_const_1.resize( number_nodes );
 for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {
  auto qfunc = new DQuadFunction();
  qfunc->add_variable( & v_sqrd_voltages[ node_id ] , -1.0 , 0.0 );
  qfunc->add_variable( & v_voltage[ node_id ] , 0.0 , 1.0 );
  v_diag_const_1[ node_id ].set_lhs( -Inf< double >() );
  v_diag_const_1[ node_id ].set_rhs( 0.0 );
  v_diag_const_1[ node_id ].set_function( qfunc );
  }
 add_static_constraint( v_diag_const_1 , "v_diag_const_1" );

 v_diag_const_2.resize( number_nodes );
 for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {
  auto lfunc = new LinearFunction();
  lfunc->add_variable( & v_sqrd_voltages[ node_id ] , -1.0 );
  lfunc->add_variable( & v_voltage[ node_id ] ,
                       min_voltage[ node_id ] * f_C_v_scal +
                       max_voltage[ node_id ] * f_C_v_scal );
  v_diag_const_2[ node_id ].set_rhs( Inf< double >() );
  v_diag_const_2[ node_id ].set_lhs( min_voltage[ node_id ] *
                                     max_voltage[ node_id ] *
                                     pow( f_C_v_scal , 2 ) );
  v_diag_const_2[ node_id ].set_function( lfunc );
  }
 add_static_constraint( v_diag_const_2 , "v_diag_const_2" );

 /* Classic McCormick relaxations for a product of variables:
  *   z_{n,n'} >= underline{v}_n v_{n'}  + underline{v}_{n'} v_n
  *               - underline{v}_n underline{v}_{n'}
  *   z_{n,n'} >= overline{v}_n  v_{n'}  + overline{v}_{n'}  v_n
  *               - overline{v}_n overline{v}_{n'}
  *   z_{n,n'} <= underline{v}_n v_{n'}  + overline{v}_{n'}  v_n
  *               - underline{v}_n overline{v}_{n'}
  *   z_{n,n'} <= overline{v}_n  v_{n'}  + underline{v}_{n'} v_n
  *               - overline{v}_n underline{v}_{n'} */
 v_def_z_1.resize( nb_dc_lines );
 i_line = 0;
 for( auto & line_id : DC_lines ) {
  Index p = start_line[ line_id ];
  Index n = end_line[ line_id ];

  auto lfunc = new LinearFunction();
  lfunc->add_variable( & v_z[ i_line ] , 1.0 );
  lfunc->add_variable( & v_voltage[ p ] , -min_voltage[ n ] * f_C_v_scal );
  lfunc->add_variable( & v_voltage[ n ] , -min_voltage[ p ] * f_C_v_scal );
  v_def_z_1[ i_line ].set_function( lfunc );
  v_def_z_1[ i_line ].set_lhs( -min_voltage[ n ] * min_voltage[ p ] *
                               pow( f_C_v_scal , 2 ) );
  v_def_z_1[ i_line ].set_rhs( Inf< double >() );
  ++i_line;
  }
 add_static_constraint( v_def_z_1 , "v_def_z_1" );

 v_def_z_2.resize( nb_dc_lines );
 i_line = 0;
 for( auto & line_id : DC_lines ) {
  Index p = start_line[ line_id ];
  Index n = end_line[ line_id ];

  auto lfunc = new LinearFunction();
  lfunc->add_variable( & v_z[ i_line ] , 1.0 );
  lfunc->add_variable( & v_voltage[ p ] , -max_voltage[ n ] * f_C_v_scal );
  lfunc->add_variable( & v_voltage[ n ] , -max_voltage[ p ] * f_C_v_scal );
  v_def_z_2[ i_line ].set_function( lfunc );
  v_def_z_2[ i_line ].set_lhs( -max_voltage[ n ] * max_voltage[ p ] *
                               pow( f_C_v_scal , 2 ) );
  v_def_z_2[ i_line ].set_rhs( Inf< double >() );
  ++i_line;
  }
 add_static_constraint( v_def_z_2 , "v_def_z_2" );

 v_def_z_3.resize( nb_dc_lines );
 i_line = 0;
 for( auto & line_id : DC_lines ) {
  Index p = start_line[ line_id ];
  Index n = end_line[ line_id ];

  auto lfunc = new LinearFunction();
  lfunc->add_variable( & v_z[ i_line ] , 1.0 );
  lfunc->add_variable( & v_voltage[ p ] , -max_voltage[ n ] * f_C_v_scal );
  lfunc->add_variable( & v_voltage[ n ] , -min_voltage[ p ] * f_C_v_scal );
  v_def_z_3[ i_line ].set_function( lfunc );
  v_def_z_3[ i_line ].set_lhs( -Inf< double >() );
  v_def_z_3[ i_line ].set_rhs( -min_voltage[ p ] * max_voltage[ n ] *
                               pow( f_C_v_scal , 2 ) );
  ++i_line;
  }
 add_static_constraint( v_def_z_3 , "v_def_z_3" );

 v_def_z_4.resize( nb_dc_lines );
 i_line = 0;
 for( auto & line_id : DC_lines ) {
  Index p = start_line[ line_id ];
  Index n = end_line[ line_id ];

  auto lfunc = new LinearFunction();
  lfunc->add_variable( & v_z[ i_line ] , 1.0 );
  lfunc->add_variable( & v_voltage[ n ] , -max_voltage[ p ] * f_C_v_scal );
  lfunc->add_variable( & v_voltage[ p ] , -min_voltage[ n ] * f_C_v_scal );
  v_def_z_4[ i_line ].set_function( lfunc );
  v_def_z_4[ i_line ].set_lhs( -Inf< double >() );
  v_def_z_4[ i_line ].set_rhs( -min_voltage[ n ] * max_voltage[ p ] *
                               pow( f_C_v_scal , 2 ) );
  ++i_line;
  }
 add_static_constraint( v_def_z_4 , "v_def_z_4" );

 /* alpha_{n,n'} <= 1 - ( 1 - cos( theta^Delta_{n,n'} ) )
  *                     / ( theta^Delta_{n,n'} )^2
  *                     ( theta_n - theta_{n'} )^2 */
 v_def_alpha_1.resize( nb_dc_lines );
 i_line = 0;
 for( auto & line_id : DC_lines ) {
  Index p = start_line[ line_id ];
  Index n = end_line[ line_id ];
  double delta_theta = PI * ( ( std::max )( v_line_max_angle[ line_id ] ,
                                            -v_line_min_angle[ line_id ] ) ) /
   180.0;
  double coeff = ( 1 - cos( delta_theta ) ) / pow( delta_theta , 2 );

  auto qfunc = new QuadFunction();
  qfunc->add_variable( & v_alpha[ i_line ] , 1.0 , 0.0 );
  qfunc->add_variable( & v_theta[ p ] , 0.0 , coeff );
  qfunc->add_variable( & v_theta[ n ] , 0.0 , coeff );
  qfunc->add_nd_term( & v_theta[ p ] , & v_theta[ n ] , -2.0 * coeff );

  v_def_alpha_1[ i_line ].set_function( qfunc );
  v_def_alpha_1[ i_line ].set_lhs( -Inf< double >() );
  v_def_alpha_1[ i_line ].set_rhs( 1.0 );
  ++i_line;
  }
 add_static_constraint( v_def_alpha_1 , "v_def_alpha_1" );

 v_def_alpha_2.resize( nb_dc_lines );
 i_line = 0;
 for( auto & line_id : DC_lines ) {
  double delta_theta = PI * ( ( std::max )( v_line_max_angle[ line_id ] ,
                                            -v_line_min_angle[ line_id ] ) ) /
   180.0;

  auto lfunc = new LinearFunction();
  lfunc->add_variable( & v_alpha[ i_line ] , 1.0 );
  v_def_alpha_2[ i_line ].set_function( lfunc );
  v_def_alpha_2[ i_line ].set_lhs( cos( delta_theta ) );
  v_def_alpha_2[ i_line ].set_rhs( Inf< double >() );
  ++i_line;
  }
 add_static_constraint( v_def_alpha_2 , "v_def_alpha_2" );

 /* McCormick relaxations of c_{n,n'}:
  *
  *  c_{n,n'} >= underline{v}_n underline{v}_{n'} alpha_{n,n'}
  *              + cos( theta^Delta_{n,n'} ) z_{n,n'}
  *              - underline{v}_n underline{v}_{n'} cos( theta^Delta_{n,n'} )
  *  c_{n,n'} >= overline{v}_n overline{v}_{n'} alpha_{n,n'}
  *              + z_{n,n'} - overline{v}_n overline{v}_{n'}
  *  c_{n,n'} <= underline{v}_n underline{v}_{n'} alpha_{n,n'}
  *              + z_{n,n'} - underline{v}_n underline{v}_{n'}
  *  c_{n,n'} <= overline{v}_n overline{v}_{n'} alpha_{n,n'}
  *              + cos( theta^Delta_{n,n'} ) z_{n,n'}
  *              - overline{v}_n overline{v}_{n'} cos( theta^Delta_{n,n'} ) */
 v_def_c_1.resize( nb_dc_lines );
 i_line = 0;
 for( auto & line_id : DC_lines ) {
  Index p = start_line[ line_id ];
  Index n = end_line[ line_id ];
  double delta_theta = PI * std::max( v_line_max_angle[ line_id ] ,
                                      -v_line_min_angle[ line_id ] ) / 180.0;
  double coeff_cos = cos( delta_theta );

  auto lfunc = new LinearFunction();
  lfunc->add_variable( & v_sum_product_voltages[ line_id ] , 1.0 );
  // v_sum is indexed over all lines, whereas the aux variables only over
  // the AC lines
  lfunc->add_variable( & v_alpha[ i_line ] ,
                       -min_voltage[ n ] * min_voltage[ p ] *
                       pow( f_C_v_scal , 2 ) );
  lfunc->add_variable( & v_z[ i_line ] , -coeff_cos );
  v_def_c_1[ i_line ].set_function( lfunc );
  v_def_c_1[ i_line ].set_lhs( -coeff_cos * min_voltage[ n ] *
                               min_voltage[ p ] *
                               pow( f_C_v_scal , 2 ) );
  v_def_c_1[ i_line ].set_rhs( Inf< double >() );
  ++i_line;
  }
 add_static_constraint( v_def_c_1 , "v_def_c_1" );

 v_def_c_2.resize( nb_dc_lines );
 i_line = 0;
 for( auto & line_id : DC_lines ) {
  Index p = start_line[ line_id ];
  Index n = end_line[ line_id ];

  auto lfunc = new LinearFunction();
  lfunc->add_variable( & v_sum_product_voltages[ line_id ] , 1.0 );
  lfunc->add_variable( & v_alpha[ i_line ] ,
                       -max_voltage[ n ] * max_voltage[ p ] *
                       pow( f_C_v_scal , 2 ) );
  lfunc->add_variable( & v_z[ i_line ] , -1.0 );
  v_def_c_2[ i_line ].set_function( lfunc );
  v_def_c_2[ i_line ].set_lhs( -max_voltage[ n ] * max_voltage[ p ] *
                               pow( f_C_v_scal , 2 ) );
  v_def_c_2[ i_line ].set_rhs( Inf< double >() );
  ++i_line;
  }
 add_static_constraint( v_def_c_2 , "v_def_c_2" );

 v_def_c_3.resize( nb_dc_lines );
 i_line = 0;
 for( auto & line_id : DC_lines ) {
  Index p = start_line[ line_id ];
  Index n = end_line[ line_id ];
  double delta_theta = PI * std::max( v_line_max_angle[ line_id ] ,
                                      -v_line_min_angle[ line_id ] ) / 180.0;
  double coeff_cos = cos( delta_theta );

  auto lfunc = new LinearFunction();
  lfunc->add_variable( & v_sum_product_voltages[ line_id ] , 1.0 );
  lfunc->add_variable( & v_alpha[ i_line ] ,
                       -max_voltage[ n ] * max_voltage[ p ] *
                       pow( f_C_v_scal , 2 ) );
  lfunc->add_variable( & v_z[ i_line ] , -coeff_cos );
  v_def_c_3[ i_line ].set_function( lfunc );
  v_def_c_3[ i_line ].set_lhs( -Inf< double >() );
  v_def_c_3[ i_line ].set_rhs( -coeff_cos * max_voltage[ n ] *
                               max_voltage[ p ] *
                               pow( f_C_v_scal , 2 ) );
  ++i_line;
  }
 add_static_constraint( v_def_c_3 , "v_def_c_3" );

 v_def_c_4.resize( nb_dc_lines );
 i_line = 0;
 for( auto & line_id : DC_lines ) {
  Index p = start_line[ line_id ];
  Index n = end_line[ line_id ];

  auto lfunc = new LinearFunction();
  lfunc->add_variable( & v_sum_product_voltages[ line_id ] , 1.0 );
  lfunc->add_variable( & v_alpha[ i_line ] ,
                       -min_voltage[ n ] * min_voltage[ p ] *
                       pow( f_C_v_scal , 2 ) );
  lfunc->add_variable( & v_z[ i_line ] , -1.0 );
  v_def_c_4[ i_line ].set_function( lfunc );
  v_def_c_4[ i_line ].set_lhs( -Inf< double >() );
  v_def_c_4[ i_line ].set_rhs( -min_voltage[ n ] * min_voltage[ p ] *
                               pow( f_C_v_scal , 2 ) );
  ++i_line;
  }
 add_static_constraint( v_def_c_4 , "v_def_c_4" );

 /* The beta variables are involved in the convex relaxation of the sine
  * function:
  *
  *  beta_{n,n'} <= cos( theta^Delta_{n,n'} / 2 )
  *                 ( ( theta_n - theta_{n'} ) - theta^Delta_{n,n'} / 2 )
  *                 + sin( theta^Delta_{n,n'} / 2 )
  *  beta_{n,n'} >= cos( theta^Delta_{n,n'} / 2 )
  *                 ( ( theta_n - theta_{n'} ) + theta^Delta_{n,n'} / 2 )
  *                 - sin( theta^Delta_{n,n'} / 2 ) */
 v_def_beta_1.resize( nb_dc_lines );
 i_line = 0;
 for( auto & line_id : DC_lines ) {
  Index p = start_line[ line_id ];
  Index n = end_line[ line_id ];
  double delta_theta = PI * ( ( std::max )( v_line_max_angle[ line_id ] ,
                                            -v_line_min_angle[ line_id ] ) ) /
   180.0;
  double coeff_cos = cos( delta_theta / 2.0 );
  double coeff_sin = sin( delta_theta / 2.0 );

  auto lfunc = new LinearFunction();
  lfunc->add_variable( & v_beta[ i_line ] , 1.0 );
  lfunc->add_variable( & v_theta[ p ] , -coeff_cos );
  lfunc->add_variable( & v_theta[ n ] , coeff_cos );
  v_def_beta_1[ i_line ].set_function( lfunc );
  v_def_beta_1[ i_line ].set_lhs( -Inf< double >() );
  v_def_beta_1[ i_line ].set_rhs( coeff_sin - coeff_cos * delta_theta / 2.0 );
  ++i_line;
  }
 add_static_constraint( v_def_beta_1 , "v_def_beta_1" );

 v_def_beta_2.resize( nb_dc_lines );
 i_line = 0;
 for( auto & line_id : DC_lines ) {
  Index p = start_line[ line_id ];
  Index n = end_line[ line_id ];
  double delta_theta = PI * ( ( std::max )( v_line_max_angle[ line_id ] ,
                                            -v_line_min_angle[ line_id ] ) ) /
   180.0;
  double coeff_cos = cos( delta_theta / 2.0 );
  double coeff_sin = sin( delta_theta / 2.0 );

  auto lfunc = new LinearFunction();
  lfunc->add_variable( & v_beta[ i_line ] , 1.0 );
  lfunc->add_variable( & v_theta[ p ] , -coeff_cos );
  lfunc->add_variable( & v_theta[ n ] , coeff_cos );
  v_def_beta_2[ i_line ].set_function( lfunc );
  v_def_beta_2[ i_line ].set_lhs( -coeff_sin + coeff_cos * delta_theta / 2.0 );
  v_def_beta_2[ i_line ].set_rhs( Inf< double >() );
  ++i_line;
  }
 add_static_constraint( v_def_beta_2 , "v_def_beta_2" );

 /* McCormick relaxations of s_{n,n'} (product involving the sine):
  *
  *  s_{n,n'} >= underline{v}_n underline{v}_{n'} beta_{n,n'}
  *              - sin( theta^Delta_{n,n'} ) z_{n,n'}
  *              + underline{v}_n underline{v}_{n'} sin( theta^Delta_{n,n'} )
  *  s_{n,n'} >= overline{v}_n overline{v}_{n'} beta_{n,n'}
  *              + sin( theta^Delta_{n,n'} ) z_{n,n'}
  *              - overline{v}_n overline{v}_{n'} sin( theta^Delta_{n,n'} )
  *  s_{n,n'} <= underline{v}_n underline{v}_{n'} beta_{n,n'}
  *              + sin( theta^Delta_{n,n'} ) z_{n,n'}
  *              - underline{v}_n underline{v}_{n'} sin( theta^Delta_{n,n'} )
  *  s_{n,n'} <= overline{v}_n overline{v}_{n'} beta_{n,n'}
  *              - sin( theta^Delta_{n,n'} ) z_{n,n'}
  *              + overline{v}_n overline{v}_{n'} sin( theta^Delta_{n,n'} ) */
 v_def_s_1.resize( nb_dc_lines );
 i_line = 0;
 for( auto & line_id : DC_lines ) {
  Index p = start_line[ line_id ];
  Index n = end_line[ line_id ];
  double delta_theta = PI * ( ( std::max )( v_line_max_angle[ line_id ] ,
                                            -v_line_min_angle[ line_id ] ) ) /
   180.0;
  double coeff_sin = sin( delta_theta );

  auto lfunc = new LinearFunction();
  lfunc->add_variable( & v_diff_product_voltages[ line_id ] , 1.0 );
  lfunc->add_variable( & v_beta[ i_line ] ,
                       -min_voltage[ n ] * min_voltage[ p ] *
                       pow( f_C_v_scal , 2 ) );
  lfunc->add_variable( & v_z[ i_line ] , coeff_sin );
  v_def_s_1[ i_line ].set_function( lfunc );
  v_def_s_1[ i_line ].set_lhs( coeff_sin * min_voltage[ n ] *
                               min_voltage[ p ] *
                               pow( f_C_v_scal , 2 ) );
  v_def_s_1[ i_line ].set_rhs( Inf< double >() );
  ++i_line;
  }
 add_static_constraint( v_def_s_1 , "v_def_s_1" );

 v_def_s_2.resize( nb_dc_lines );
 i_line = 0;
 for( auto & line_id : DC_lines ) {
  Index p = start_line[ line_id ];
  Index n = end_line[ line_id ];
  double delta_theta = PI * ( ( std::max )( v_line_max_angle[ line_id ] ,
                                            -v_line_min_angle[ line_id ] ) ) /
   180.0;
  double coeff_sin = sin( delta_theta );

  auto lfunc = new LinearFunction();
  lfunc->add_variable( & v_diff_product_voltages[ line_id ] , 1.0 );
  lfunc->add_variable( & v_beta[ i_line ] ,
                       -max_voltage[ n ] * max_voltage[ p ] *
                       pow( f_C_v_scal , 2 ) );
  lfunc->add_variable( & v_z[ i_line ] , -coeff_sin );
  v_def_s_2[ i_line ].set_function( lfunc );
  v_def_s_2[ i_line ].set_lhs( -coeff_sin * max_voltage[ n ] *
                               max_voltage[ p ] *
                               pow( f_C_v_scal , 2 ) );
  v_def_s_2[ i_line ].set_rhs( Inf< double >() );
  ++i_line;
  }
 add_static_constraint( v_def_s_2 , "v_def_s_2" );

 v_def_s_3.resize( nb_dc_lines );
 i_line = 0;
 for( auto & line_id : DC_lines ) {
  Index p = start_line[ line_id ];
  Index n = end_line[ line_id ];
  double delta_theta = PI * ( ( std::max )( v_line_max_angle[ line_id ] ,
                                            -v_line_min_angle[ line_id ] ) ) /
   180.0;
  double coeff_sin = sin( delta_theta );

  auto lfunc = new LinearFunction();
  lfunc->add_variable( & v_diff_product_voltages[ line_id ] , 1.0 );
  lfunc->add_variable( & v_beta[ i_line ] ,
                       -min_voltage[ n ] * min_voltage[ p ] *
                       pow( f_C_v_scal , 2 ) );
  lfunc->add_variable( & v_z[ i_line ] , -coeff_sin );
  v_def_s_3[ i_line ].set_function( lfunc );
  v_def_s_3[ i_line ].set_lhs( -Inf< double >() );
  v_def_s_3[ i_line ].set_rhs( -coeff_sin * min_voltage[ n ] *
                               min_voltage[ p ] *
                               pow( f_C_v_scal , 2 ) );
  ++i_line;
  }
 add_static_constraint( v_def_s_3 , "v_def_s_3" );

 v_def_s_4.resize( nb_dc_lines );
 i_line = 0;
 for( auto & line_id : DC_lines ) {
  Index p = start_line[ line_id ];
  Index n = end_line[ line_id ];
  double delta_theta = PI * ( ( std::max )( v_line_max_angle[ line_id ] ,
                                            -v_line_min_angle[ line_id ] ) ) /
   180.0;
  double coeff_sin = sin( delta_theta );

  auto lfunc = new LinearFunction();
  lfunc->add_variable( & v_diff_product_voltages[ line_id ] , 1.0 );
  lfunc->add_variable( & v_beta[ i_line ] ,
                       -max_voltage[ n ] * max_voltage[ p ] *
                       pow( f_C_v_scal , 2 ) );
  lfunc->add_variable( & v_z[ i_line ] , coeff_sin );
  v_def_s_4[ i_line ].set_function( lfunc );
  v_def_s_4[ i_line ].set_lhs( -Inf< double >() );
  v_def_s_4[ i_line ].set_rhs( coeff_sin * max_voltage[ n ] *
                               max_voltage[ p ] *
                               pow( f_C_v_scal , 2 ) );
  ++i_line;
  }
 add_static_constraint( v_def_s_4 , "v_def_s_4" );

 }  // end( ACNetworkBlock::strengthen_SOCP_relaxation )

/*--------------------------------------------------------------------------*/

std::vector< std::pair< double , double > >
 ACNetworkBlock::recover_feasible_solution( void )
{
 /* Since the solution provided from the AC OPF relaxation problem is not
  * necessarily feasible, we implement a feasibility-recovery algorithm. */

 std::vector< std::pair< double , double > > v_feasible_sol;
 // each pair is the real and imaginary part

 // 1) first, get the solution of the relaxation problem
 std::vector< double > relaxed_power_flow;
 std::vector< double > relaxed_reactive_power_flow;

 std::transform( v_power_flow.begin() , v_power_flow.end() ,
                 relaxed_power_flow.begin() ,
                 []( const ColVariable & v ) { return( v.get_value() ); }
                 );
 std::transform( v_reactive_power_flow.begin() ,
                 v_reactive_power_flow.end() ,
                 relaxed_reactive_power_flow.begin() ,
                 []( const ColVariable & v ) { return( v.get_value() ); }
                 );

 // 1b) multiply back the obtained solutions by the earlier scale factor
 //     since indeed we have computed v_power_flow_tilde, and we care for
 //     v_power_flow -> the relation is v_power_flow_tilde = C * v_power_flow

 // 2) then compute spanning tree
 auto result = f_NetworkData->get_cycle_basis();

 // 3) do some magic (TODO)

 return( v_feasible_sol );
 }  // end( ACNetworkBlock::recover_feasible_solution )

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

NetworkBlockSolution * ACNetworkBlock::new_Solution( void ) const
{
 return( new ACNetworkBlockSolution() );
 }  // end( ACNetworkBlock::new_Solution )

/*--------------------------------------------------------------------------*/
/*------------------ METHODS OF ACNetworkBlockSolution ---------------------*/
/*--------------------------------------------------------------------------*/

void ACNetworkBlockSolution::deserialize( const netCDF::NcGroup & group )
{
 // call the method of the base class
 DCNetworkBlockSolution::deserialize( group );

 // deserialize the (reactive) flow injection - - - - - - - - - - - - - - -
 ::deserialize< double >( group , "NodeInjectionReactive" ,
                          v_node_injection_reactive , false );

 // deserialize the (reactive) flow variables - - - - - - - - - - - - - - -
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

 // deserialize the (reactive) flow injection - - - - - - - - - - - - - - -
 auto ncVar = group.getVar( "NodeInjectionReactive" );
 if( ncVar.isNull() )
  v_node_injection_reactive.clear();
 else {
  std::vector< size_t > cnt = { 1 , f_number_nodes };
  v_node_injection_reactive.resize( f_number_nodes );
  ncVar.getVar( strt , cnt , v_node_injection_reactive.data() );
  }

 // deserialize the (reactive) flow variables - - - - - - - - - - - - - - -
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
  // read the (reactive) node injection variables - - - - - - - - - - - - -
  auto RNI = ACNB->get_const_reactive_node_injection();
  for( Index n = 0 ; n < f_number_nodes ; ++n )
   v_node_injection_reactive[ n ] = ( *( RNI++ ) ).get_value();
  }

 if( ! v_reactive_flow_from.empty() ) {
  // read the (reactive) flow power variables - - - - - - - - - - - - - - -
  auto RFli = ACNB->get_const_reactive_power_flow().begin();
  for( Index l = 0 ; l < f_number_lines ; ++l )
   v_reactive_flow_from[ l ] = ( *( RFli++ ) ).get_value();
  for( Index l = 0 ; l < f_number_lines ; ++l )
   v_reactive_flow_to[ l ] = ( *( RFli++ ) ).get_value();
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
  // write the (reactive) node injection variables- - - - - - - - - - - - -
  auto RNI = ACNB->get_reactive_node_injection();
  for( Index n = 0 ; n < f_number_nodes ; ++n )
   ( *( RNI++ ) ).set_value( v_node_injection_reactive[ n ] );
  }

 if( ! v_reactive_flow_from.empty() ) {
  // write the (reactive) flow power variables- - - - - - - - - - - - - - -
  auto RFli = ACNB->get_reactive_power_flow().begin();
  for( Index l = 0 ; l < f_number_lines ; ++l )
   ( *( RFli++ ) ).set_value( v_reactive_flow_from[ l ] );
  for( Index l = 0 ; l < f_number_lines ; ++l )
   ( *( RFli++ ) ).set_value( v_reactive_flow_to[ l ] );
  }
 }  // end( ACNetworkBlockSolution::write )

/*--------------------------------------------------------------------------*/

void ACNetworkBlockSolution::serialize( netCDF::NcGroup & group ) const
{
 // call the method of the base class
 DCNetworkBlockSolution::serialize( group );

 // serialize the (reactive) flow injection - - - - - - - - - - - - - - - -
 if( ! v_node_injection_reactive.empty() ) {
  auto nn = group.getDim( "NumberNodes" );
  ::serialize< double >( group , "NodeInjectionReactive" ,
                         netCDF::NcDouble() , nn ,
                         v_node_injection_reactive );
  }

 // serialize the (reactive) flow variables - - - - - - - - - - - - - - - -
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
                                        size_t idx ) const
{
 // call the method of the base class
 DCNetworkBlockSolution::serialize( group , idx );

 // now serialize the data structures - - - - - - - - - - - - - - - - - - -

 netCDF::NcVar RNI;   // NodeInjectionReactive
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
   RNI = group.getVar( "NodeInjectionReactive" );

  if( ! v_reactive_flow_from.empty() ) {
   RFFV = group.getVar( "ReactiveFlowFromValue" );
   RFTV = group.getVar( "ReactiveFlowToValue" );
   }
  }

 std::vector< size_t > strt = { idx , 0 };

 if( ! RFFV.isNull() ) {  // reactive power injection has to be serialised
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
 // NetworkBlockSolution which calls clone() and therefore returns an
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
   v_reactive_flow_from[ l ] += ACNBS->v_reactive_flow_from[ l ] *
                                multiplier;
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
