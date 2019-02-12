/*--------------------------------------------------------------------------*/
/*--------------------- File NetworkBlock.cpp ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the NetworkBlock class.
 *
 * \version 0.10
 *
 * \date 03 - 09 - 2016
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Kostas Tavlaridis-Gyparakis \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Antonio Frangioni, Kostas Tavlaridis-Gyparakis
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/
#include <map>
#include "NetWorkBlock.h"
#include "UCBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace std;
using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

NetWorkBlock::NetWorkBlock(UCBlock * flbock) : Block(flbock) {
	    
}

/*--------------------------------------------------------------------------*/

NetWorkBlock::~NetWorkBlock() {

}

/*--------------------------------------------------------------------------*/

std::map<std::string,NetWorkBlock::NetWorkFactory>& NetWorkBlock::f_factory()
{
 static std::map<std::string,NetWorkBlock::NetWorkFactory>  *  ans =
     new  std::map<std::string,NetWorkBlock::NetWorkFactory>(); 
     //initializing the NetWorkFactory
 return  * ans;
 }

