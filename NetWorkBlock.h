/*--------------------------------------------------------------------------*/
/*--------------------------- File NetWorkBlock.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *derived* class NetworkBlcok, which derives from the 
 * Block, in order to define a very basic interface for any possible derived 
 * type of Network of a Unit Commitement Block. Τhe basis of Network is consi-
 * dered to be very generic and thus with the minimum possible ingrients and
 * in the case of this base class is restricted to only a very basic Demand
 * Constraint.
 * Based on the above description the class has been constructed having the
 * following elements:
 *
 * - A virtual public method in used in order to initialize and read the data
 *   of any possible derived NetwrokBlock class.
 * - A factory that is used in order for any possible derived class to be able 
 *   to "automatically" register itself and be initialized. The Factory is de-
 *   fined by a static method that initializes the static map that is used in
 *   order to store all the different possible derived classes that are linked
 *   with a unique string.
 * - A vector of doubles used to store the values that refer to the rhs of a
 *   global demand constraints that links together all the units of the exa-
 *   mined UC Problem.
 * - A pointer to the UCBlock to which any derived NetwrokBlock class is at-
 *   tached.  
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

#ifndef NETWORKBLOCK_H_
#define NETWORKBLOCK_H_ /* self-identification: #endif at the end
				     * of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "boost/function.hpp"
#include "Block.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it {


class UCBlock; ///< forward declaration of UCBlock


class NetWorkBlock : public Block {

/*--------------------------------------------------------------------------*/
/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

protected:

/*--------------------------------------------------------------------------*/
/*---------------------- PROTECTED TYPES OF THE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

typedef boost::function<NetWorkBlock * (UCBlock *)> NetWorkFactory;
/**< Definition of the NetworkFactory, used to properly initialize all the di-
ferent possible derived classes of NetworkBlock */

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:
/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */
 /// constructor of NetWorkBlock, taking possibly a pointer of its fater Block
 NetWorkBlock(UCBlock * flbock = nullptr);

/*--------------------------------------------------------------------------*/

 virtual ~NetWorkBlock(); 
///< destructor of NetWorkBlock: it is virtual, and empty

/*@} -----------------------------------------------------------------------*/
/*------------ METHODS FOR READING THE DATA OF THE NetWorkBlock ------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the NetWorkBlock
    @{ */

 static std::map<std::string,NetWorkBlock::NetWorkFactory>& f_factory();
 /**< Static Member Method for initialization of the map of the Factory of
      NetworkBlock*/

 void load( std::istream &input ) =0;
/**< load method of the Block, pure virtual method of the Block, 
that will be used to load the data */




/*@} -----------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

protected:

 std::vector<double> d; ///<vector to store the demand of the network

};


}

#endif /* NETWORKBLOCK_H_ */
