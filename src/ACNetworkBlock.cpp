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


#ifndef PI
    #define PI 3.14159265358979323846
#endif

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;
using namespace std::complex_literals;

typedef Eigen::SparseMatrix< std::complex<double> > SpCMat;
typedef Eigen::SparseVector< std::complex<double> > SpCVec;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register ACNetworkBlock to the Block factory

SMSpp_insert_in_factory_cpp_1( ACNetworkBlock );

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS OF ACNetworkBlock ------------------------*/
/*--------------------------------------------------------------------------*/

void ACNetworkBlock::generate_abstract_variables( Configuration * stvv )
{
  // QJ: print for debug
  std::cout << "LineResistance size : "  << f_NetworkData->get_line_resistance().size() << std::endl;
  std::cout << "LineReactance size : "   << f_NetworkData->get_line_reactance().size() << std::endl;
  std::cout << "LineSusceptance size : " << f_NetworkData->get_line_susceptance().size() << std::endl;
  std::cout << "NodeSuceptance size : "  << f_NetworkData->get_node_susceptance().size() << std::endl;
  std::cout << "NodeConductance size : " << f_NetworkData->get_node_conductance().size() << std::endl;
  std::cout << "NodeMinVoltage size : "  << f_NetworkData->get_node_min_voltage().size() << std::endl;
  std::cout << "NodeMaxVoltage size : "  << f_NetworkData->get_node_max_voltage().size() << std::endl;
  std::cout << "LineMinAngle size : "    << f_NetworkData->get_line_min_angle().size() << std::endl;
  std::cout << "LineMaxAngle size : "    << f_NetworkData->get_line_max_angle().size() << std::endl;

  // generate first the same variables as in the DCNetwork
  DCNetworkBlock::generate_abstract_variables(stvv);


  const auto number_nodes = get_number_nodes();
  const auto number_lines = get_number_lines();

  // ----- complex power flow (real and imaginary part for both directions)
  /*
  S denotes the AC power for each line, therefore it is a 2-dimensional vector 
  (real and imaginary part).
  */
  S_power_flow.resize( boost::extents[ 2 ][ 2*number_lines ] );
  for( int i = 0; i < 2; ++i ) {
    for( Index line_id = 0 ; line_id < 2*number_lines ; ++line_id ) {
      S_power_flow[ i ][ line_id ].set_type( ColVariable::kContinuous );
    }
  }
  add_static_variable( S_power_flow , "S_power_flow" );

  // -----
  v_sum_product_voltages.resize(number_lines);
  for( Index line_id = 0 ; line_id < number_lines ; ++line_id ) {
    v_sum_product_voltages[ line_id ].set_type( ColVariable::kContinuous );
  }
  v_diff_product_voltages.resize(number_lines);
  for( Index line_id = 0 ; line_id < number_lines ; ++line_id ) {
    v_diff_product_voltages[ line_id ].set_type( ColVariable::kContinuous );
  }
  v_sqrt_voltages.resize(number_nodes);
  for( Index node_id = 0 ; node_id < number_nodes; ++node_id ) {
    v_sqrt_voltages[ node_id ].set_type( ColVariable::kContinuous );
  }
  add_static_variable(v_sum_product_voltages, "v_sum_product_voltages");
  add_static_variable(v_diff_product_voltages, "v_diff_product_voltages");
  add_static_variable(v_sqrt_voltages, "v_sqrt_voltages");

};

// ---------------------------------------
void ACNetworkBlock::generate_objective( Configuration * objc){
  // TODO : use the coefficient of the matpower instance
  DCNetworkBlock::generate_objective(objc);
};


// ---------------------------------------
void ACNetworkBlock::generate_abstract_constraints( Configuration * stcc ){
  if( constraints_generated() )  // constraints have already been generated
  return;                       // nothing to do

  const auto number_nodes = get_number_nodes();

  if( number_nodes <= 1 )
    return;

  const auto number_lines = get_number_lines();

  if( number_lines <= 0 )
    throw( std::logic_error( "ACNetworkBlock::generate_abstract_constraints: "
                           "number of lines of DCNetworkBlock is not set" ) );

  const auto & start_line = f_NetworkData->get_start_line();
  const auto & end_line = f_NetworkData->get_end_line();

  std::vector<Index> AC_lines = get_AC_lines();
  double base_mva = f_NetworkData->get_baseMVA();


  // ----- Voltage bounds
  /*
  We ensure that the voltage magnitude is bounded between min_voltage and max_voltage.
  To this aim, W_{n,n} = (V_n).(V_n)^H = |V_n|^2 so we impose the bounds directly on W_{n,n}.
  */
  v_voltage_bounds_const.resize(number_nodes);
  const auto & min_voltage = f_NetworkData->get_node_min_voltage();
  const auto & max_voltage = f_NetworkData->get_node_max_voltage();
  for( Index n = 0 ; n < number_nodes ; ++n ) {
    v_voltage_bounds_const[ n ].set_lhs( pow(min_voltage[n],2) );
    v_voltage_bounds_const[ n ].set_rhs( pow(max_voltage[n],2) );
    v_voltage_bounds_const[ n ].set_variable( & v_sqrt_voltages[n] );
  }
  add_static_constraint( v_voltage_bounds_const, "AC_voltage_bounds_limit" );

  // ----- Angle bounds
  /*
  We aim to incorporate Phase Angle Differnce (PAD) constraints for each line = (start,end):
    min_angle_{line} <= angle(V_{end}) - angle(V_{start}) <= max_angle_{line}
  This constraint can be taken into account using:
     tan(min_angle_{line}) Real(W_{start,end})<= Imag(W_{start,end}) <= max_angle_{line} Real(W_{start,end})
  */
  v_angle_bounds_const.resize(boost::multi_array< FRowConstraint , 2 >::extent_gen()[ 2 ][ number_lines ] );
  const auto & min_angle = f_NetworkData->get_line_min_angle();
  const auto & max_angle = f_NetworkData->get_line_max_angle();
  for (Index line_id = 0; line_id < number_lines; ++line_id) {
    double phi_min = PI * min_angle[line_id] / 180.;
    double phi_max = PI * max_angle[line_id] / 180.;
    Index p = start_line[line_id];
    Index n = end_line[line_id];
    auto lfunc_1 = new LinearFunction();
    lfunc_1->add_variable( & v_diff_product_voltages[ line_id ], 1.0);
    lfunc_1->add_variable( & v_sum_product_voltages[ line_id ], -tan(phi_min));
    v_angle_bounds_const[0][ line_id ].set_lhs( 0.0 );
    v_angle_bounds_const[0][ line_id ].set_rhs( Inf< double >() );
    v_angle_bounds_const[0][ line_id ].set_function( lfunc_1 );
    auto lfunc_2 = new LinearFunction();
    lfunc_2->add_variable( & v_diff_product_voltages[ line_id ], 1.0);
    lfunc_2->add_variable( & v_sum_product_voltages[ line_id ], -tan(phi_max));
    v_angle_bounds_const[1][ line_id ].set_lhs( -Inf< double >() );
    v_angle_bounds_const[1][ line_id ].set_rhs( 0.0 );
    v_angle_bounds_const[1][ line_id ].set_function( lfunc_2 );
  }
  add_static_constraint( v_angle_bounds_const, "AC_angle_bounds_limit" );

  // ----- Active and Reactive Power conservation:
  /*
  We define 
    M_{b} = Ys_b E_{bb} + \sum_{a:(b,a) in lines} (Yff_{ba}E_{bb} + Yft_{ba}E_{ba}) 
          + \sum_{a:(a,b) in lines} (Ytt_{ba}E_{aa} + Ytf_{ba}E_{ab})
  where Ys is the shunt admittance, Yff, Yft, Ytt, Ytf are the line impedances and E_{ab} is the classical basis matrices.

  Then, the power conservation constraint can be written
    Supply - Demand = <M,W>_F

  The complex matrix product <M,W>_F is then decomposed into a real part and an imaginary part.
  */
  v_power_flow_injection_const.resize(1*number_nodes); // QJ: TODO reactive power conservation

  // real part of the power flow conservation
  for( Index p = 0 ; p < number_nodes ; ++p ) {
    auto lfunc = new LinearFunction();
    lfunc->add_variable( & v_node_injection[ 0 ][ p ] , -1.0 );
    lfunc->add_variable( & v_sqrt_voltages[ p ] , -ACdata.Ys.coeff(p).real() / base_mva );

    for (Index line_id = 0; line_id < number_lines; ++line_id) {
      Index i = start_line[line_id];
      Index j = end_line[line_id];

      if (i == p) lfunc->add_variable( & S_power_flow[0][ line_id ], 1.);
      if (j == p) lfunc->add_variable( & S_power_flow[0][ number_lines + line_id ], 1.);
    }
    v_power_flow_injection_const[ p ].set_both( -v_ActiveDemand[ p ] / base_mva );
    v_power_flow_injection_const[ p ].set_function( lfunc );
  }

  add_static_constraint( v_power_flow_injection_const, "AC_power_flow_injection" );

  // ----- Definition of complex power flow
  /*
  We define the complex power flow S_{line} for each line = (start,end) as 
    S_{start,end} = Yff_{start,end} W_{start,start} + Yft_{start,end}W_{start,end}
    S_{end,start} = Ytt_{start,end} W_{start,start} + Ytf_{start,end}W_{end,end}
  Onece, again, we then split into two constraints (one for real and one for imaginary part)
  */
  v_voltage_definition_const.resize(boost::multi_array< FRowConstraint , 2 >::extent_gen()[ 2 ][ 2*number_lines ] );
  for (Index line_id = 0; line_id < number_lines; ++line_id) {
    Index p = start_line[line_id];
    Index n = end_line[line_id];
    auto lfunc_1 = new LinearFunction();
    lfunc_1->add_variable( & v_sqrt_voltages[ p ], ACdata.v_admittance[ line_id ].real());
    lfunc_1->add_variable( & v_sum_product_voltages[ line_id ], - ACdata.v_admittance[ line_id ].real());
    lfunc_1->add_variable( & v_diff_product_voltages[ line_id ], - ACdata.v_admittance[ line_id ].imag());
    lfunc_1->add_variable( & S_power_flow[0][line_id], -1.0);
    v_voltage_definition_const[0][ line_id ].set_both(0.0);
    v_voltage_definition_const[0][ line_id ].set_function( lfunc_1 );
    auto lfunc_2 = new LinearFunction();
    lfunc_2->add_variable( & v_sqrt_voltages[ p ],  - ACdata.v_admittance[ line_id ].imag());
    lfunc_2->add_variable( & v_sum_product_voltages[ line_id ],  ACdata.v_admittance[ line_id ].imag());
    lfunc_2->add_variable( & v_diff_product_voltages[ line_id ], -ACdata.v_admittance[ line_id ].real());
    lfunc_2->add_variable( & S_power_flow[1][line_id], -1.0);
    v_voltage_definition_const[1][ line_id ].set_both(0.0);
    v_voltage_definition_const[1][ line_id ].set_function( lfunc_2 );
  }
  for (Index line_id = 0; line_id < number_lines; ++line_id) {
    Index p = start_line[line_id];
    Index n = end_line[line_id];
    auto lfunc_1 = new LinearFunction();
    lfunc_1->add_variable( & v_sqrt_voltages[ n ], ACdata.v_admittance[ line_id ].real());
    lfunc_1->add_variable( & v_sum_product_voltages[ line_id ], - ACdata.v_admittance[ line_id ].real());
    lfunc_1->add_variable( & v_diff_product_voltages[ line_id ], ACdata.v_admittance[ line_id ].imag());
    lfunc_1->add_variable( & S_power_flow[0][number_lines + line_id], -1.0);
    v_voltage_definition_const[0][ number_lines + line_id ].set_both(0.0);
    v_voltage_definition_const[0][ number_lines + line_id ].set_function( lfunc_1 );
    auto lfunc_2 = new LinearFunction();
    lfunc_2->add_variable( & v_sqrt_voltages[ n ],  - ACdata.v_admittance[ line_id ].imag());
    lfunc_2->add_variable( & v_sum_product_voltages[ line_id ],  ACdata.v_admittance[ line_id ].imag());
    lfunc_2->add_variable( & v_diff_product_voltages[ line_id ], ACdata.v_admittance[ line_id ].real());
    lfunc_2->add_variable( & S_power_flow[1][number_lines + line_id], -1.0);
    v_voltage_definition_const[1][ number_lines + line_id ].set_both(0.0);
    v_voltage_definition_const[1][ number_lines + line_id ].set_function( lfunc_2 );
  }
  add_static_constraint( v_voltage_definition_const, "AC_voltage_defintion_const" );

  // ----- Thermal limit on lines
  /*
  We impose that |S_{line}| <= rateA_{line}, which corresponds to a thermal limitation.
  Then, to take into account this constraint, we use a DQuadfunction:
    Real(S_{line})^2 + Imag(S_{line})^2 <= rateA_{line}^2
  */
  v_thermal_limit.resize(2*number_lines);
  const auto & rate_A = f_NetworkData->get_line_rate_A();
  for (Index line_id = 0; line_id < number_lines; ++line_id) {
    Index p = start_line[line_id];
    Index n = end_line[line_id];
    auto qfunc_1 = new DQuadFunction();
    qfunc_1->add_variable( & S_power_flow[0][line_id], 0.0, 1.0);
    qfunc_1->add_variable( & S_power_flow[1][line_id], 0.0, 1.0);
    v_thermal_limit[ line_id ].set_lhs( -Inf< double >() );
    v_thermal_limit[ line_id ].set_rhs( pow(rate_A[line_id]/base_mva, 2) );
    v_thermal_limit[ line_id ].set_function( qfunc_1 );
    auto qfunc_2 = new DQuadFunction();
    qfunc_2->add_variable( & S_power_flow[0][number_lines + line_id], 0.0, 1.0);
    qfunc_2->add_variable( & S_power_flow[1][number_lines + line_id], 0.0, 1.0);
    v_thermal_limit[ number_lines + line_id ].set_lhs( -Inf< double >() );
    v_thermal_limit[ number_lines + line_id ].set_rhs( pow(rate_A[line_id]/base_mva, 2) );
    v_thermal_limit[ number_lines + line_id ].set_function( qfunc_2 );
  }
  add_static_constraint( v_thermal_limit, "AC_thermal_limit_const" );

  // -----QJ: up to now, only SOCP relaxation is available but it could be replaced by something else
  generate_SOCP_relaxation();
 };

// ---------------------------------------
 /*
 Links between generic variables 
    v_sum_product_voltages,
    v_diff_product_voltages,
    v_sqrt_voltages
  using a SOCP relaxation.
 */
void ACNetworkBlock::generate_SOCP_relaxation(){

  const auto number_nodes = get_number_nodes();
  const auto number_lines = get_number_lines();
  const auto & start_line = f_NetworkData->get_start_line();
  const auto & end_line = f_NetworkData->get_end_line();

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

  v_socp_const.resize( number_lines );
  for (Index line_id = 0; line_id < number_lines; ++line_id) {
    Index p = start_line[line_id];
    Index n = end_line[line_id];
    auto qfunc = new QuadFunction();
    qfunc->add_variable( & v_sqrt_voltages[ p ], 0.0, 0.0);
    qfunc->add_variable( & v_sqrt_voltages[ n ], 0.0, 0.0); 
    qfunc->add_nd_term( & v_sqrt_voltages[ p ], & v_sqrt_voltages[ n ], -1.0);
    qfunc->add_variable( & v_sum_product_voltages[ line_id ], 0.0, 1.0);
    v_socp_const[ line_id ].set_lhs( -Inf< double >() );
    v_socp_const[ line_id ].set_rhs( 0.0 );
    v_socp_const[ line_id ].set_function( qfunc );
  }
  add_static_constraint(v_socp_const, "AC_socp_const" );
  
 };



// ---------------------------------------
void ACNetworkBlock::add_ACdata(Index interval, Index node, UnitBlock* unit_block, Index t, Index g ) { 
  const auto & start_line = f_NetworkData->get_start_line();
  const auto & end_line = f_NetworkData->get_end_line();
  const auto number_nodes = get_number_nodes();
  const auto number_lines = get_number_lines();

  // Line impedances for AC network (old quantities, not used in optimization)
  ACdata.Yff = SpCMat(number_nodes,number_nodes);
  ACdata.Yft = SpCMat(number_nodes,number_nodes);
  ACdata.Ytf = SpCMat(number_nodes,number_nodes);
  ACdata.Ytt = SpCMat(number_nodes,number_nodes);
  
  for( auto& line_id : get_AC_lines() ) {
    Index i = start_line[line_id];
    Index j = end_line[line_id];
    double r = f_NetworkData->get_line_resistance().at(line_id);
    double x = f_NetworkData->get_line_reactance().at(line_id);
    double b = f_NetworkData->get_line_susceptance().at(line_id);
    double tau = f_NetworkData->get_line_ratio().at(line_id);
    if (abs(tau) < 1e-4) tau = 1.0;
    double theta = PI * f_NetworkData->get_line_angle().at(line_id) / 180;
    ACdata.Yff.insert(i,j) = (1./(r+1i*x) + 1i*b/2.)/pow(tau,2);
    ACdata.Yft.insert(i,j) = -1./((r+1i*x)*tau*exp(-1i*theta));
    ACdata.Ytf.insert(i,j) = -1./((r+1i*x)*tau*exp(1i*theta));
    ACdata.Ytt.insert(i,j) = 1./(r+1i*x) + 1i*b/2.;
  }

  // Shunt admittance
  ACdata.Ys = SpCVec(number_nodes);
  double base_mva = f_NetworkData->get_baseMVA();
  for(Index n = 0; n < number_nodes; ++n){
    double Gs = f_NetworkData->get_node_conductance().at(n)/base_mva;
    double Bs = f_NetworkData->get_node_susceptance().at(n)/base_mva;
    ACdata.Ys.insert(n) = Gs + 1i*Bs;
  }

  // New version without explicit definition of Yff, Yft, Ytf, Ytt
  ACdata.v_admittance = std::vector< std::complex<double> >(number_lines, 0.);
  ACdata.v_transformer = std::vector< std::complex<double> >(number_lines, 0.);
  for( Index line_id = 0 ; line_id < number_lines ; ++line_id ) {
    Index i = start_line[line_id];
    Index j = end_line[line_id];
    double r = f_NetworkData->get_line_resistance().at(line_id);
    double x = f_NetworkData->get_line_reactance().at(line_id);
    double b = f_NetworkData->get_line_susceptance().at(line_id);
    double tau = 1.;//f_NetworkData->get_line_ratio().at(line_id);
    if (abs(tau) < 1e-4) tau = 1.0;
    double theta = PI * f_NetworkData->get_line_angle().at(line_id) / 180;
    ACdata.v_admittance[ line_id ] = (1./(r+1i*x) + 1i*b/2.)/pow(tau,2);
    ACdata.v_transformer[ line_id ] = exp(1i*theta);
  }

};

/*--------------------------------------------------------------------------*/
/*--------------------- End File ACNetworkBlock.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
