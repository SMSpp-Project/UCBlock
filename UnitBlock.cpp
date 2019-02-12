/*--------------------------------------------------------------------------*/
/*------------------------ File UnitBlock.cpp ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the UnitBlock class.
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
#include "UnitBlock.h"
#include "UCBlock.h"


/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

using namespace std;

/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

UnitBlock::UnitBlock(UCBlock * flbock) : Block(flbock) {


UCBlock * f_UC_Block = dynamic_cast <UCBlock *> (get_f_Block());
  if(f_UC_Block){
	 U.resize(f_UC_Block->get_t());
	 P.resize(f_UC_Block->get_t());

	 for(int i = 0 ; i<f_UC_Block->get_t() ; i++){
		 //set pointer of the father Block
		 U[i].set_Block(this); 
		 P[i].set_Block(this); 
		 //set lower and upper bounds
		 U[i].set_lb(0.0); 
		 U[i].set_ub(1.0); 
		 P[i].set_lb(0.0);
		 P[i].set_ub(Inf<double>()); 
		 //set type of variables
		 U[i].set_type(ColVariable::binary); 
		 P[i].set_type(ColVariable::continuous);
	 }
  
         add_static_variable(U); //adding the commitement variables to the Static Variables Vector
 
	 add_static_variable(P); //adding the power variables to the Static Variables Vector
  }

}

/*--------------------------------------------------------------------------*/

UnitBlock::~UnitBlock() {

}
void UnitBlock::set_Unit_Block ( UCBlock * flbock ){

set_f_Block(flbock);

UCBlock * f_UC_Block = dynamic_cast <UCBlock *> (get_f_Block());
  if(f_UC_Block){
	 U.resize(f_UC_Block->get_t());
	 P.resize(f_UC_Block->get_t());

	 for(int i = 0 ; i<f_UC_Block->get_t() ; i++){
		 //set pointer of the father Block
		 U[i].set_Block(this); 
		 P[i].set_Block(this); 
		 //set lower and upper bounds
		 U[i].set_lb(0.0); 
		 U[i].set_ub(1.0); 
		 P[i].set_lb(0.0);
		 P[i].set_ub(Inf<double>()); 
		 //set type of variables
		 U[i].set_type(ColVariable::binary); 
		 P[i].set_type(ColVariable::continuous);
	 }
  
         add_static_variable(U); //adding the commitement variables to the Static Variables Vector
 
	 add_static_variable(P); //adding the power variables to the Static Variables Vector
  }


}
/*--------------------------------------------------------------------------*/

std::map<std::string,UnitBlock::UnitFactory>& UnitBlock::U_factory()
{

 static std::map<std::string,UnitBlock::UnitFactory>*  ans =
    new std::map<std::string,UnitBlock::UnitFactory>();
//Initializing the UnitFactory

 return *ans;
 }
