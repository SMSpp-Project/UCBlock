/*--------------------------------------------------------------------------*/
/*--------------------- File ACNetworkBlock.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/

#include <map>

#include "NetworkBlock.h"

#include "DCNetworkBlock.h"

#include "ACNetworkBlock.h"

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

  for (auto it_map = m_generators.begin(); it_map != m_generators.end(); ++it_map){
    std::cout << "Generator data for node " << it_map->first.first << ", interval " << it_map->first.second << " is :" << std::endl;
    for (auto& generator : it_map->second){
      std::cout << "\tMinPower: "           << generator.f_MinPower << std::endl;
      std::cout << "\tMaxPower: "           << generator.f_MaxPower << std::endl;
      std::cout << "\tMinReactivePower: "   << generator.f_MinReactivePower << std::endl;
      std::cout << "\tMaxReactivePower: "   << generator.f_MaxReactivePower << std::endl;
      std::cout << "\tCost Model: "         << generator.f_CostModel << std::endl;
      std::cout << "\tNumber cost coeffs: " << generator.f_NumberCostCoeffs << std::endl;
      int i = 0;
      for (auto & coeff : generator.v_PowerCostCoeffs){
        std::cout << "\t\tPower cost coeffs " << i << ": " << coeff << std::endl;
        ++i;
      }
    }
  }

  DCNetworkBlock::generate_abstract_variables(stvv);
}

/*--------------------------------------------------------------------------*/
/*--------------------- End File ACNetworkBlock.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
