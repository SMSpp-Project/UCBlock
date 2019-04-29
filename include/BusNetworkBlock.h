/*--------------------------------------------------------------------------*/
/*--------------------------- File BusNetWorkBlock.h -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *derived* class BusNetworkBlcok, which derives from the 
 * NetworkBlock, in order to define the Network of the Academic Version of the 
 * Unit Commitement Problem. The BusNetworkBlock class is characterized by the 
 * following ingredinets:
 *
 * - A global Demand Constraint that needs to be satisfied from all the diffe-
 *   rent Units of the Academic UC throughout all the time steps of the opti-
 *   misation horizon.
 - - A global Reserve Power Constraint that needs to be satisfied from all the 
 *   different Units of the Academic UC throughout all the time steps of the 
 *   optimisation horizon. 
 *
 * Based on the above description the class has been constructed having the
 * following elements:
 *
 * - A public method that reads and initializes all the data that descri-
 *   bes the BusNetwork instance.
 * - A public method that initializes and stores in the static vector of the
 *   Block the Demand and Reserve Power Constraints.
 * - A vector of doubles to store the rhs values of the Reserve Power Constra-
 *   ints
 * - Two vectors of Linear Constraint Objects that are used in order to store
 *   the information regarding the Demand and Reserve Power Constraints.
 * - A small private class with sole purpose of defining a static member, which 
 *   has to be initialized when main() is executed. Τhe class constructos is 
 *   set to register BusNetworkBlock in the Factory of BusNetwork.
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
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/
#ifndef BUSNETWORKBLOCK_H_
#define BUSNETWORKBLOCK_H_

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "NetWorkBlock.h"
#include "UnitBlock.h"
#include "LinearConstraint.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it {

class BusNetworkBlock : public NetWorkBlock {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:
/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */
 /** constructor of BusNetWorkBlock, taking possibly a pointer of its father 
     Block */
 BusNetworkBlock(UCBlock * flbock = nullptr);

/*--------------------------------------------------------------------------*/

 virtual ~BusNetworkBlock();
///< destructor of BusNetWorkBlock: it is virtual, and empty
/*@} -----------------------------------------------------------------------*/
/*------------ METHODS FOR READING THE DATA OF THE BusNetworkBlock ---------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the BusNetworkBlock
    @{ */

 void load(std::istream& inStream);
 ///< Method for receiving data and constructing the constraints of the Network

/*@}------------------------------------------------------------------------*/
/*--------- METHODS CONSTRUCTING THE CONSTRAINTS OF BusNetworkBlock --------*/
/*--------------------------------------------------------------------------*/
/** @name Methods constructing the constraints of BusNetworkBlock
 *  @{ */

 void generate_static_constraints( void );
 ///< Method for generating and adding the demand and reserve constraint

/*@} -----------------------------------------------------------------------*/
/*----------- METHODS FOR HANDLING THE DATA OF THE BusNetworkBlock ---------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the data of the UnitBlock
    @{ */

    std::vector<LinearConstraint>* get_Demand_Constraint() {return &Demand_Const;} 
    ///< Method for returning the vector of commitement variables
    
    LinearConstraint * get_D(int i) {return &Demand_Const[i];} 
    ///< Method for returning the i-th commitement variables

    std::vector<LinearConstraint> * get_Reserve_Constraint() {return &Reserve_Const;} 
    ///< Method for returning the vector of power variables

    LinearConstraint * get_R(int i) {return &Reserve_Const[i];} 
    ///< Method for returning the i-th commitement variables

    
/*@}------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

protected:

/*--------------------------------------------------------------------------*/
/*------------------------------ CONSTRAINTS  ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constraints of BusNetworkBlock
 *  @{ */
 std::vector<LinearConstraint > Demand_Const;
 ///< vector to store Demand Constraints
 
std::vector<LinearConstraint> Reserve_Const;
///< vector to store Reserve Power Constraints

/*@} -----------------------------------------------------------------------*/
/*----------------------- VARIABLES OF THE CLASS ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Variables of BusNetworkBlock
 *  @{ */

 std::vector<double> r; 
 ///< vector to store the spinning reserve values


/*@} -----------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/**< Very small, "fake" class _initializer. Its only meaning is to define a
   static member _init that is initialized at the very beginning of the
   main(). Hence the constructor is called, and the constructor registers
   the BusNetworkBlock class into the (static) NetWorkBlock::f_factory. */

 static class _init {
         public:
          _init();
         }_initializer;
 };
 } /* namespace SMSpp_di_unipi_it */

#endif /* BUSNETWORKBLOCK_H_ */
