/*--------------------------------------------------------------------------*/
/*--------------------------- File DCNetworkBlock.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 *
 * Header file for class DCNetworkBlock, which derives from NetworkBlock and
 * defines the standard linear constraints corresponding to the "DC model"
 * of the transmission network in the Unit Commitment problem.
 *
 * \version 0.11
 *
 * \date 18 - 07 - 2019
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
 * Copyright &copy; by Antonio Frangioni, and Ali Ghezelsoflu
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __DCNetworkBlock
#define __DCNetworkBlock
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Block.h"
#include "FRowConstraint.h"
#include "NetworkBlock.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it {

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS DCNetworkBlock ---------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// NetworkBlock with more than one node, i.e., a "DC" transmission network
/** The DCNetworkBlock class derives from NetworkBlock, and defines the
 * standard linear constraints corresponding to the "DC model" of the
 * transmission network in the Unit Commitment problem. */

class DCNetworkBlock : public NetworkBlock {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
 *
 * DCNetworkBlock defines a main public type:
 *
 @{ */

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

 /// constructor of DCNetworkBlock
 /** Constructor of DCNetworkBlock, taking possibly a pointer of its
  * father Block. */

 DCNetworkBlock( Block * fblock = nullptr ) : NetworkBlock( fblock ) ,
  f_NetworkData( nullptr ) , f_local_NetworkData( false ) { }

/*--------------------------------------------------------------------------*/
/// destructor of DCNetworkBlock
 ~DCNetworkBlock() override {
  delete f_NetworkData;
 }

/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */
/// loads the DCNetworkBlock instance from memory
/** Like load( std::istream & ), if there is any Solver attached to this
 *  DCNetworkBlock then a NBModification (the "nuclear option") is issued.
 */
 void load( std::istream & input ) override {
  throw ( std::logic_error( "DCNetworkBlock::load() not implemented yet" ) );
 }
/*--------------------------------------------------------------------------*/
///generate abstract constraints of DCNetworkBlock
/** The topology of the transmission network is defined by a set of nodes
 *
 *   \f$ N \f$ and a set of lines connecting the nodes \f$ L \f$.
 *
 *   By considering a \f$ |L| \times |N| \f$ matrix
 *   \f$ B_t \f$, which constitutes the so-called Power Transfer Distribution
 *   Factor matrix which represents the linear relationship between power
 *   injections at each node of the grid and active power flows through the
 *   transmission lines.
 *
 * The flow limit equations can be written as follow:
 *
 * \f[
 *  P^{mn}_{\ell } \leq \sum_{ n' \in N} (B)_({\ell, n'})
 *  (S_{ n'} - D^{ac}_{n'}) \leq
 *  P^{mx}_{\ell}  \quad \ell \in \mathcal{L} \quad
 * \f]
 *   Where \f$ P^{mn}_{\ell }\f$ and \f$ P^{mx}_{\ell }\f$ are minimum and
 *   maximum power flow at each line \f$ \ell \in L \f$ and \f$ S_{n'} \f$ and
 *   \f$ D^{ac}_{n' } \f$ is the node injection variable and active power
 *   demand at node \f$ n \in N \f$ in the network respectively.
 */
 void generate_abstract_constraints( Configuration *stcc = nullptr )
 override;

/**@} ----------------------------------------------------------------------*/
/*--------------- METHODS FOR MODIFYING THE DCNetworkBlock -----------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the DCNetworkBlock
 *  @{ */

 void set_NetworkData( NetworkBlock::NetworkData * network_data = nullptr )
 override
 {
  // if there was a previous NetworkData and it was local, delete it
  if( f_NetworkData && f_local_NetworkData )
   delete f_NetworkData;

  f_NetworkData = network_data;
  f_local_NetworkData = false;
  }
/**@} ----------------------------------------------------------------------*/
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

 /// the NetworkData object
 NetworkBlock::NetworkData * f_NetworkData;

 /// true if the NetworkData object has not been passed from outside
 bool f_local_NetworkData;

 /// flow limit constraints
 std::vector<FRowConstraint> v_flow_limit_constraints;

  };   // end( class( DCNetworkBlock ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  /* namespace SMSpp_di_unipi_it */

/*--------------------------------------------------------------------------*/

#endif /* DCNetworkBlock.h included */

/*--------------------------------------------------------------------------*/
/*------------------------ End File DCNetworkBlock.h -----------------------*/
/*--------------------------------------------------------------------------*/
