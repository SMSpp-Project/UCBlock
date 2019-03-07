/*--------------------------------------------------------------------------*/
/*--------------------------- File NetworkBlock.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 *
 * Header file for the *derived* class NetworkBlock, which derives
 * from Block, in order to define a very basic interface for any
 * possible derived type of network of a UCBlock. Τhe basis of
 * NetworkBlock is considered to be very generic and thus with the
 * minimum possible ingredients. Based on the above description, the
 * class has been constructed having the following elements:
 *
 * - A virtual public method is used in order to initialize and read
 *   the data of any possible derived NetworkBlock class.
 *
 * - A factory that is used in order for any possible derived class to
 *   be able to "automatically" register itself and be
 *   initialized. The factory is defined by a static method that
 *   initializes the static map that is used in order to store all the
 *   different possible derived classes that are linked with a unique
 *   string.
 *
 * - An int that stores the time this NetworkBlock is associated with.
 *
 * - A pointer to a Network, which defines the network.
 *
 * - A vector of doubles used to store the values of the demand at
 *   each node of the network.
 *
 * - A vector of doubles used to store the values of the susceptance
 *   of each line of the network.
 *
 * - Two vectors of doubles to store the minimum and maximum power
 *   flow in each line of the network.
 *
 * - A pointer to the UCBlock to which any derived NetworkBlock class
 *   is attached.
 *
 * - Flow limit constraints.
 *
 * \version 0.11
 *
 * \date 05 - 03 - 2019
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Ali Ghezelsoflu \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Rafael Durbano Lobato \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Kostas Tavlaridis-Gyparakis \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Antonio Frangioni, Ali Ghezelsoflu, Rafael
 * Durbano Lobato, and Kostas Tavlaridis-Gyparakis
 */

/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __NetworkBlock
#define __NetworkBlock /* self-identification: #endif at the end of
			* the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "boost/function.hpp"
#include "Block.h"
#include "ColVariable.h"
#include "FRowConstraint.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it {

class UCBlock; ///< forward declaration of UCBlock

class NetworkBlock : public Block {

/*--------------------------------------------------------------------------*/
/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

protected:

/*--------------------------------------------------------------------------*/
/*---------------------- PROTECTED TYPES OF THE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

 /** Definition of the NetworkFactory, used to properly initialize all
  * the different possible derived classes of NetworkBlock. */
 typedef boost::function<NetworkBlock * (UCBlock *)> NetworkFactory;

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

 /// constructor of NetworkBlock, taking possibly a pointer to its father Block

 NetworkBlock( Block * fblock = nullptr );

/*--------------------------------------------------------------------------*/

 /// destructor of NetworkBlock

 virtual ~NetworkBlock();

/*@} -----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

 virtual void generate_abstract_variables( Configuration *stvv = nullptr )
   override;

 virtual void generate_abstract_constraints( Configuration *stcc = nullptr )
   override;

/*@} -----------------------------------------------------------------------*/
/*--------------- METHODS FOR MODIFYING THE NetworkBlock -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the NetworkBlock
 *  @{ */

 /// sets the time this NetworkBlock is associated with
 void set_time( int t ) {
   f_time = t;
 }

 /// sets the Network of this NetworkBlock
 void set_network( Network * network ) {
   f_network = network;
 }

/*@} -----------------------------------------------------------------------*/
/*----------- METHODS FOR READING THE DATA OF THE NetworkBlock -------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the NetworkBlock
    @{ */


/*@} -----------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

 /// method incapsulating the NetworkBlock factory
 /** This method returns the NetworkBlock factory, which is a static object.
  * The rationale for using a method is that this is the "Construct On First
  * Use Idiom" that solves the "static initialization order problem". */

 static std::map<std::string, NetworkBlock::NetworkFactory> & f_factory();

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 /// the time associated with this NetworkBlock
 int f_time;

 /// a pointer to the Network
 Network * f_network;

 /// vector to store the demand of each node of the network
 std::vector<double> v_demand;

 /// vector to store the susceptance of each line of the network
 std::vector<double> v_susceptance;

 /// vector to store the minimum power flow at each line
 std::vector<double> v_minimum_power_flow;

 /// vector to store the maximum power flow at each line
 std::vector<double> v_maximum_power_flow;

 /// flow limit constraints
 std::vector<FRowConstraint> v_flow_limit_constraints;

};   // end( class( NetworkBlock ) )

}  /* namespace SMSpp_di_unipi_it */

#endif /* NetworkBlock.h included */

/*--------------------------------------------------------------------------*/
/*------------------------ End File NetworkBlock.h -------------------------*/
/*--------------------------------------------------------------------------*/
