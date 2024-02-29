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

  struct GeneratorData 
  {
    double f_MinPower;
    double f_MaxPower;
    double f_MinReactivePower;
    double f_MaxReactivePower;
    double f_VoltageMagnitude;
    std::vector<double> v_PowerCostCoeffs;
    Index f_CostModel;
    Index f_NumberCostCoeffs;
  };

explicit ACNetworkBlock( Block * f_block = nullptr )
  : DCNetworkBlock( f_block ) {}

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

virtual void add_generator_data(Index interval, Index node, UnitBlock* unit_block, Index t, Index g ) override 
{
  GeneratorData gdata;
  gdata.f_MinPower = unit_block->get_min_power(t, g);
  gdata.f_MaxPower = unit_block->get_max_power(t, g);
  gdata.f_MinReactivePower = unit_block->get_min_reactive_power(t, g);
  gdata.f_MaxReactivePower = unit_block->get_max_reactive_power(t, g);
  gdata.f_VoltageMagnitude = unit_block->get_voltage_magnitude(t, g);
  gdata.f_NumberCostCoeffs = unit_block->get_number_cost_coeffs();
  gdata.f_CostModel = unit_block->get_cost_model();
  for (int i = 0; i < gdata.f_NumberCostCoeffs; ++i){
    gdata.v_PowerCostCoeffs.push_back(unit_block->get_cost_coeff(i,g));
  }

  auto it = m_generators.find(std::make_pair(node,interval));
  if (it != m_generators.end()){
    it->second.push_back(gdata);
  }
  else {
    std::vector<GeneratorData> v_gdata = {gdata};
    m_generators[std::make_pair(node,interval)] = v_gdata; 
  }
};

 
 protected:

 std::map<std::pair<Index,Index>, std::vector<GeneratorData>> m_generators;
 
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
