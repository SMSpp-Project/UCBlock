/*--------------------------------------------------------------------------*/
/*--------------------- File ACNetworkBlock.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/

#include <map>

#include <complex>

#include "NetworkBlock.h"

#include "DCNetworkBlock.h"

#include "ACNetworkBlock.h"

#include "LinearFunction.h"

#include "OneVarConstraint.h"

#include "FRealObjective.h"

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
  std::cout << "LineResistance size : "  << f_NetworkData->get_line_resistance().size() << std::endl;
  std::cout << "LineReactance size : "   << f_NetworkData->get_line_reactance().size() << std::endl;
  std::cout << "LineSusceptance size : " << f_NetworkData->get_line_susceptance().size() << std::endl;
  std::cout << "NodeSuceptance size : "  << f_NetworkData->get_node_susceptance().size() << std::endl;
  std::cout << "NodeConductance size : " << f_NetworkData->get_node_conductance().size() << std::endl;

  DCNetworkBlock::generate_abstract_variables(stvv);


  const auto number_nodes = get_number_nodes();
  const auto number_lines = get_number_lines();

  // voltage relaxation matrix W = V.V^H
  W_voltage_real.resize( boost::extents[ number_nodes ][ number_nodes ] );
  W_voltage_imag.resize( boost::extents[ number_nodes ][ number_nodes ] );

  for( Index p = 0 ; p < number_nodes ; ++p ) {
    for( Index n = 0 ; n < number_nodes ; ++n ) {
      W_voltage_real[ p ][ n ].set_type( ColVariable::kContinuous );
      W_voltage_imag[ p ][ n ].set_type( ColVariable::kContinuous );
    }
  }
  add_static_variable( W_voltage_real , "W_voltage_real" );
  add_static_variable( W_voltage_imag , "W_voltage_imag" );
}


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


  // ----- Voltage bounds
  v_voltage_bounds_const.resize(number_nodes);
  const auto & min_voltage = f_NetworkData->get_node_min_voltage();
  const auto & max_voltage = f_NetworkData->get_node_max_voltage();
  for( Index n = 0 ; n < number_nodes ; ++n ) {
    v_voltage_bounds_const[ n ].set_lhs( pow(min_voltage[n],2) );
    v_voltage_bounds_const[ n ].set_rhs( pow(max_voltage[n],2) );
    v_voltage_bounds_const[ n ].set_variable( & W_voltage_real[n][n] );
  }


  // ----- Active and Reactive Power conservation: "Supply - Demand = <M,W>_F"
  v_power_flow_injection_const.resize(2*number_nodes);

  for( Index p = 0 ; p < number_nodes ; ++p ) {
    auto lfunc = new LinearFunction();
    lfunc->add_variable( & v_node_injection[ 0 ][ p ] , -1.0 );

    for( Index n = 0 ; n < number_nodes ; ++n ) {

      lfunc->add_variable( & W_voltage_real[p][n] , ACdata.M.coeff(p,n).real() );
      lfunc->add_variable( & W_voltage_imag[p][n] , ACdata.M.coeff(p,n).imag() );
    }
    v_power_flow_injection_const[ p ].set_both( -v_ActiveDemand[ p ] );
    v_power_flow_injection_const[ p ].set_function( lfunc );
  }
  for( Index p = 0 ; p < number_nodes ; ++p ) {
    auto lfunc = new LinearFunction();
    lfunc->add_variable( & v_node_injection[ 1 ][ p ] , -1.0 );

    for( Index n = 0 ; n < number_nodes ; ++n ) {

      lfunc->add_variable( & W_voltage_real[p][n] , -ACdata.M.coeff(p,n).imag() );
      lfunc->add_variable( & W_voltage_imag[p][n] , ACdata.M.coeff(p,n).real() );
    }
    v_power_flow_injection_const[ p ].set_both( -v_ActiveDemand[ p ] ); // TODO must be ReactiveDemand
    v_power_flow_injection_const[ p ].set_function( lfunc );
  }


 };



// ---------------------------------------
void ACNetworkBlock::add_ACdata(Index interval, Index node, UnitBlock* unit_block, Index t, Index g ) { 
  std::cout << "Add ACdata to node " << node << std::endl;
  const auto & start_line = f_NetworkData->get_start_line();
  const auto & end_line = f_NetworkData->get_end_line();
  const auto number_nodes = get_number_nodes();

  // classical quantities for AC network
  ACdata.Yff = SpCMat(number_nodes,number_nodes);
  ACdata.Yft = SpCMat(number_nodes,number_nodes);
  ACdata.Ytf = SpCMat(number_nodes,number_nodes);
  ACdata.Ytt = SpCMat(number_nodes,number_nodes);
  ACdata.M   = SpCMat(number_nodes,number_nodes);
  
  for( auto& line_id : get_AC_lines() ) {
    Index i = start_line[line_id];
    Index j = end_line[line_id];
    double r = f_NetworkData->get_line_resistance().at(line_id);
    double x = f_NetworkData->get_line_reactance().at(line_id);
    double b = f_NetworkData->get_line_susceptance().at(line_id);
    double tau = f_NetworkData->get_line_ratio().at(line_id);
    double theta = PI * f_NetworkData->get_line_angle().at(line_id) / 180;
    ACdata.Yff.insert(i,j) = (1./(r+1i*x) + 1i*b/2.)/pow(tau,2);
    ACdata.Yft.insert(i,j) = -1./((r+1i*x)*tau*exp(-1i*theta));
    ACdata.Ytf.insert(i,j) = -1./((r+1i*x)*tau*exp(1i*theta));
    ACdata.Ytt.insert(i,j) = 1./(r+1i*x) + 1i*b/2.;
  }

  // line shunt
  for(Index n = 0; n < number_nodes; ++n){
    double Gs = f_NetworkData->get_node_conductance().at(n);
    double Bs = f_NetworkData->get_node_susceptance().at(n);
    ACdata.Ys.insert(n) = Gs + 1i*Bs;
  }

  // Construct matrix of constraints
  for(Index i = 0; i < number_nodes; ++i){
    for (SpCMat::InnerIterator itY(ACdata.Yff,i); itY; ++itY) {
      ACdata.M.coeffRef(i,i) += itY.value();
    }
    for (SpCMat::InnerIterator itY(ACdata.Yft,i); itY; ++itY) {
      Index j = itY.row();
      ACdata.M.coeffRef(i,j) += itY.value();
    }
    for (SpCMat::InnerIterator itY(ACdata.Ytt,i); itY; ++itY) {
      Index j = itY.row();
      ACdata.M.coeffRef(j,j) += itY.value();
    }
    for (SpCMat::InnerIterator itY(ACdata.Ytf,i); itY; ++itY) {
      Index j = itY.row();
      ACdata.M.coeffRef(j,i) += itY.value();
    }
    ACdata.M.coeffRef(i,i) += ACdata.Ys.coeff(i);
  }


  // matrices HM and ZM (real part and imaginary part)
  ACdata.HM = 0.5*(ACdata.M + ACdata.M.conjugate());
  ACdata.ZM = 0.5*(ACdata.M - ACdata.M.conjugate());



  std::cout << "ACdata added" << std::endl;
};

/*--------------------------------------------------------------------------*/
/*--------------------- End File ACNetworkBlock.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
