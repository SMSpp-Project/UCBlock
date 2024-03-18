/*--------------------------------------------------------------------------*/
/*--------------------------- File ACNetworkBlock.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 *
 * Header file for class ACNetworkBlock, which derives from DCNetworkBlock and
 * defines the standard SOCP relaxation corresponding to the "AC model"
 * of the transmission network in the Unit Commitment problem.
 *
 * @ Quentin J.
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __ACNetworkBlock
 #define __ACNetworkBlock
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Block.h"

#include "LinearFunction.h"

#include "FRowConstraint.h"

#include "OneVarConstraint.h"

#include "DCNetworkBlock.h"

#include "FRealObjective.h"

#include <iostream>
/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)

namespace SMSpp_di_unipi_it
{



class ACNetworkBlock : public DCNetworkBlock
{

 public:

  struct ACNetworkData 
  {
    std::map< std::pair<Index,Index>, std::complex<double> > Yff;
    std::map< std::pair<Index,Index>, std::complex<double> > Yft;
    std::map< std::pair<Index,Index>, std::complex<double> > Ytf;
    std::map< std::pair<Index,Index>, std::complex<double> > Ytt;

    std::map<Index, std::complex<double> > Ys;

    std::map< std::pair<Index,Index>, std::complex<double> > M;
  };

explicit ACNetworkBlock( Block * f_block = nullptr )
  : DCNetworkBlock( f_block ) {}

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

virtual void add_ACdata(Index interval, Index node, UnitBlock* unit_block, Index t, Index g ) override 
{
  const double PI  =3.141592653589793238463;
  std::complex<double> jj;
  const auto & start_line = f_NetworkData->get_start_line();
  const auto & end_line = f_NetworkData->get_end_line();

  for( auto& line_id : get_AC_lines() ) {
    Index i = start_line[line_id];
    Index j = end_line[line_id];
    double r = f_NetworkData->get_line_resistance(line_id);
    double b = f_NetworkData->get_line_reactance(line_id);
    double bb = f_NetworkData->get_line_susceptance(line_id);
    double tau = f_NetworkData->get_line_ratio(line_id);
    double theta = PI * f_NetworkData->get_line_angle(line_id) / 180;
    ACdata.Yff[std::make_tuple(i,j)] = (1./(r+jj*x) + jj*bb/2)/tau**2;
    ACdata.Yft[std::make_tuple(i,j)] = -1./((r+jj*x)*tau*exp(-jj*theta));
    ACdata.Ytf[std::make_tuple(i,j)] = -1./((r+jj*x)*tau*exp(jj*theta));
    ACdata.Ytt[std::make_tuple(i,j)] = 1./(r+jj*x) + jj*bb/2;
  }
  for(int n = 0; n < get_number_nodes(); ++n){
    double Gs = f_NetworkData->get_node_conductance();
    double Bs = f_NetworkData->get_node_susceptance();
    ACdata.Ys[n] = Gs + jj*Bs;
  }

  // Construct matrix of constraints
  std::complex<double> Mval;
  for( auto it = ACdata.Yff.begin(); it != ACdata.Yff.end(); ++it ) {
    std::complex<double> y = it->second;
    Index i = get<0>(it->first);
    Index j = get<1>(it->first);
    auto t = std::make_tuple(i,i);
    auto itM = ACdata.M.find(t);
    if (itM != ACdata.end()) Mval = itM->second;
    else Mval = 0.;
    ACdata.M[t] = Mval + y;
  }
  for( auto it = ACdata.Yft.begin(); it != ACdata.Yft.end(); ++it ) {
    std::complex<double> y = it->second;
    Index i = get<0>(it->first);
    Index j = get<1>(it->first);
    auto t = std::make_tuple(i,j);
    auto itM = ACdata.M.find(t);
    if (itM != ACdata.end()) Mval = itM->second;
    else Mval = 0.;
    ACdata.M[t] = Mval + y;
  }
  for( auto it = ACdata.Ytt.begin(); it != ACdata.Ytt.end(); ++it ) {
    std::complex<double> y = it->second;
    Index i = get<0>(it->first);
    Index j = get<1>(it->first);
    auto t = std::make_tuple(j,j);
    auto itM = ACdata.M.find(t);
    if (itM != ACdata.end()) Mval = itM->second;
    else Mval = 0.;
    ACdata.M[t] = Mval + y;
  }
  for( auto it = ACdata.Ytf.begin(); it != ACdata.Ytf.end(); ++it ) {
    std::complex<double> y = it->second;
    Index i = get<0>(it->first);
    Index j = get<1>(it->first);
    auto t = std::make_tuple(j,i);
    auto itM = ACdata.M.find(t);
    if (itM != ACdata.end()) Mval = itM->second;
    else Mval = 0.;
    ACdata.M[t] = Mval + y;
  }
  for(int n = 0; n < get_number_nodes(); ++n){
    auto t = std::make_tuple(n,n);
    auto itM = ACdata.M.find(t);
    if (itM != ACdata.end()) Mval = itM->second;
    else Mval = 0.;
    ACdata.M[t] = Mval + y;
  }

};

 
 protected:

 ACNetworkData ACdata;
 std::map<std::pair<Index,Index>, std::pair<ColVariable,ColVariable>> m_voltage;
 
 private:

 SMSpp_insert_in_factory_h;


 static void static_initialization( void ) {

  register_method< ACNetworkBlock , MF_dbl_it , Subset && , bool >(
   "DCNetworkBlock::set_active_demand" , &ACNetworkBlock::set_active_demand );

  register_method< ACNetworkBlock , MF_dbl_it , Range >(
   "DCNetworkBlock::set_active_demand" , &ACNetworkBlock::set_active_demand );
 }

};  // end( class( DCNetworkBlock ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* ACNetworkBlock.h included */

/*--------------------------------------------------------------------------*/
/*-------------------- End File ACNetworkBlock.h ---------------------------*/
/*--------------------------------------------------------------------------*/
