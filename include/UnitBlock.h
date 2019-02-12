/*--------------------------------------------------------------------------*/
/*--------------------------- File UnitBlock.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *derived* class UnitkBlcok, which derives from the 
 * Block, in order to define a base for any possible ThermalUnit that can be 
 * attached to a UCBlock, having the very basic information that can characte-
 * rize almost any different kind of thermal Unit, which in the case of this 
 * base class are restricted to only two very basic sets of Variables, one for
 * the Power Variables and one for the Commitement Variables of each Unit.
 * Based on the above description the class has been constructed having the
 * following elements:
 *
 * - A virtual public method in used in order to initialize and read the data
 *   of any possible derived UnitBlock class.
 * - A factory that is used in order for any possible derived class to be able 
 *   to "automatically" register itself and be initialized. The Factory is de-
 *   fined by a static method that initializes the static map that is used in
 *   order to store all the different possible derived classes that are linked
 *   with a unique string.
 * - Two vector of ColVariable Objects, that are used in order to store the in-
 *   formation regarding the Power and Commitement Variables of each Unit.
 * - A pointer to the UCBlock to which any derived UnitBlock class is attached.
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

#ifndef UNITBLOCK_H_
#define UNITBLOCK_H_  /* self-identification: #endif at the end
				     * of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "boost/function.hpp"
#include "Block.h"
#include "ColVariable.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it {

class UCBlock; ///< forward declaration of UCBlock

class UnitBlock : public Block {

/*--------------------------------------------------------------------------*/
/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

protected:

/*--------------------------------------------------------------------------*/
/*---------------------- PROTECTED TYPES OF THE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

 typedef boost::function <UnitBlock*(UCBlock *)> UnitFactory;
/**< Definition of the NetworkFactory, used to properly initialize all the di-
ferent possible derived classes of UnitBlock */

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/


public:
/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

 UnitBlock( UCBlock * flbock = nullptr );
/**< Constructor of UnitBlock, taking possibly a pointer of its fater Block
and also initializing the vectors of Power and Commitement Variables */

/*--------------------------------------------------------------------------*/

 virtual ~UnitBlock();
///< destructor of NetWorkBlock: it is virtual, and empty

/*@} -----------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE DATA OF THE UnitBlock -------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the UnitBlock
    @{ */

 static std::map<std::string,UnitBlock::UnitFactory>& U_factory();
 /**< Static Member Method for initialization of the map of the Factory of
  UnitBlock */

void set_Unit_Block ( UCBlock * flbock );

 void load( std::istream &input ) =0;
 /**< load method of the Block, pure virtual method of the Block, 
 that will be used to load the data */

/*@} -----------------------------------------------------------------------*/
/*-------------- METHODS FOR HANDLING THE DATA OF THE UnitBlock ------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the data of the UnitBlock
    @{ */

    std::vector<ColVariable>* get_U() {return &U;} 
    ///< Method for returning the vector of commitement variables
    
    ColVariable * get_U(int i) {return &U[i];} 
    ///< Method for returning the i-th commitement variables

    int get_total_U() {return U.size();}
    ///< Method for returning the total number of commitement variables

    std::vector<ColVariable> * get_P() {return &P;} 
    ///< Method for returning the vector of power variables

    ColVariable * get_P(int i) {return &P[i];} 
    ///< Method for returning the i-th commitement variables

    int get_total_P() {return P.size();}
    ///< Method for returning the total number of power variables

/*@} -----------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

protected:

 std::vector<ColVariable> U;  ///< Vector of commitement variables
 std::vector<ColVariable> P;  ///< Vector of power variables

 };



} /* namespace SMSpp_di_unipi_it */


#endif /* UNITBLOCK_H_ */
