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
 * - The number of nodes and lines of the network.
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
 * - Flow limit constraints.
 *
 * \version 0.11
 *
 * \date 19 - 03 - 2019
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

#include "Block.h"
#include "ColVariable.h"
#include "FRowConstraint.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it {

class NetworkBlock : public Block {

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

 virtual void load( std::istream &input ) override { };

/*@} -----------------------------------------------------------------------*/
/*--------------- METHODS FOR MODIFYING THE NetworkBlock -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the NetworkBlock
 *  @{ */

 /// sets the number of nodes of the network
 void set_number_nodes( int number_nodes ) {
   f_number_nodes = number_nodes;
 }
 /// sets the number of lines of the network
 void set_number_lines( int number_lines ) {
   f_number_lines = number_lines;
 }

/*@} -----------------------------------------------------------------------*/
/*----------- METHODS FOR READING THE DATA OF THE NetworkBlock -------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the NetworkBlock
    @{ */

 /// returns the vector of node injection variables
 const std::vector<ColVariable> & get_node_injection() const {
   return v_node_injection;
 }

/*@} -----------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 /// number of nodes of the network
 int f_number_nodes = -1;

 /// number of lines of the network
 int f_number_lines = -1;

 /// vector to store the demand of each node of the network
 std::vector<double> v_demand;

 /// vector to store the susceptance of each line of the network
 std::vector<double> v_susceptance;

 /// vector to store the minimum power flow at each line
 std::vector<double> v_minimum_power_flow;

 /// vector to store the maximum power flow at each line
 std::vector<double> v_maximum_power_flow;

 /// power injection at each node
 std::vector<ColVariable> v_node_injection;

 /// flow limit constraints
 std::vector<FRowConstraint> v_flow_limit_constraints;

};   // end( class( NetworkBlock ) )

}  /* namespace SMSpp_di_unipi_it */

#endif /* NetworkBlock.h included */

/*--------------------------------------------------------------------------*/
/*------------------------ End File NetworkBlock.h -------------------------*/
/*--------------------------------------------------------------------------*/
