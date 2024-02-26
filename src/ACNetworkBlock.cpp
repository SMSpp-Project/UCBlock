/*--------------------------------------------------------------------------*/
/*--------------------- File DCNetworkBlock.cpp ----------------------------*/
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

// register DCNetworkBlock to the Block factory

SMSpp_insert_in_factory_cpp_1( ACNetworkBlock );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

// register DCNetworkData to the NetworkData factory

typedef DCNetworkBlock::DCNetworkData DCNetworkData;

SMSpp_insert_in_factory_cpp_1( DCNetworkData );

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS OF DCNetworkBlock ------------------------*/
/*--------------------------------------------------------------------------*/

ACNetworkBlock::~ACNetworkBlock()
{
  //TODO
}

void ACNetworkBlock::generate_abstract_variables( Configuration * stvv )
{
  std::cout << "QJ test override for AC" << std::endl;
  DCNetworkBlock::generate_abstract_variables(stvv);
}

/*--------------------------------------------------------------------------*/
/*--------------------- End File DCNetworkBlock.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
