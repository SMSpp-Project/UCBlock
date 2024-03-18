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

  DCNetworkBlock::generate_abstract_variables(stvv);
}

 void generate_abstract_constraints( Configuration * stcc = nullptr ){
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


 // TODO


 };

/*--------------------------------------------------------------------------*/
/*--------------------- End File ACNetworkBlock.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
