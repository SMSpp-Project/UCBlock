/*--------------------------------------------------------------------------*/
/*--------------------------- File UCBlock.h -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *derived* class UCBlcok, which
 * implements the base class Block, in order to define a very basic Unit Com-
 * mitement Block, that will be able to fit and be used as a base for almost 
 * any different variation of a Unit Commitement Problem. As a result of this 
 * the UCBlock class is characterized by the following ingridients:
 *
 * - A Network (that can potentially refer to satisfy of demand, reserve ener-,
 *   gy and any other global/linking constraint, characteristic)
 * - A set of Units (that may be referring to any different kind of Units, such
 *   as thermal, hydro, etc) 
 * - A time-horizon of the optimization problem
 * - A size of the total number of units included in the optimization problem.
 * 
 * Based on the above description the classes has been constructed having the
 * following elemnts:
 *
 * - Some basic public methods to read the data and initialize the optimisa-
 *   tion problem 
 * - Two public integer variables to store the time-horizon and the number of
 *   units
 * - A Vector of Pointers to objects of UnitBlock class, which is the base 
 *   class for any possible derived type of Unit that may be referring to a
 *   UC Problem.
 * - A Pointer to an object of NetWorkBlock class, that is the base class of
 *   any possible derived type of Network that may be attached to a UC Prob-
 *   lem.
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
#ifndef UCBLOCK_H_
#define UCBLOCK_H_ /* self-identification: #endif at the end
				     * of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Block.h"
#include "DQuadObjectiveFunction.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it {

 class UnitBlock;  ///< forward declaration of UnitBlock
 class NetWorkBlock;  ///< forward declaration of NetworkBlock


class UCBlock : public Block {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */
 /// constructor of UCBlock, taking possibly a pointer of its fater Block
 UCBlock(Block *father = nullptr) : Block(father) {}

/*--------------------------------------------------------------------------*/
 virtual ~UCBlock();   ///< destructor of UCBlock: it is virtual, and empty

/*@} -----------------------------------------------------------------------*/
/*--------------- METHODS FOR READING THE DATA OF THE UCBlock --------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the UCBlock
    @{ */

void load( std::istream &input ) {};  
/**< load method of the Block, that needs to be declared, as it is a pure virtual
method of te Block, in this point is initialized with empty body */

void instance(std::istream& inStream); 
///< Method that is initializing the instance and passing all the needed data

/*@} -----------------------------------------------------------------------*/
/*--------------- METHODS FOR HANDLING THE DATA OF THE UCBlock -------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the data of the UCBlock
    @{ */

 int get_t() {return t;} 
 ///< returns the time horizon of the Problem
 void set_t(int taf) {t=taf;} 
 ///< sets the time horizon of the Problem

 int get_units_size() {return units_size;} 
 ///< returns the number of thermal units of the Problem

 void set_units_size(int units) {units_size = units;} 
 ///< sets the number of thermal units of the Problem

 std::vector<UnitBlock *> get_units() {return Units;}

/*@} -----------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

protected:

 int t; ///< the time horizon of the Problem
 int units_size; ///< the number of thermal units of the Problem

 std::vector<UnitBlock *> Units; 
/**< vector of pointers of the objects of UnitBlock that refer to any possible
different Unit that can be derived from UnitBlock */

 NetWorkBlock * Network; 
/**< Pointer of the object of NetworkBlock that refer to any possible different 
Unit that can be derived from NetworkBlock */

LinearObjectiveFunction l_of;
DQuadObjectiveFunction  q_of;  
///<Objective function of the Academic Thermal Block

};

} /* namespace SMSpp_di_unipi_it */


#endif /* UCBLOCK_H_ */
