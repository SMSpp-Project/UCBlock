/*--------------------------------------------------------------------------*/
/*--------------------------- File DCNetworkBlock.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 *
 * Header file for the *derived* class DCNetworkBlock, which derives from
 * NetworkBlock, in order to define a very basic interface for any possible
 * derived type of dc-network of a UCBlock. Τhe basis of DCNetworkBlock is
 * considered to be very generic and thus with the minimum possible
 * ingredients.
 *
 * \version 0.11
 *
 * \date 24 - 06 - 2019
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
 * Copyright &copy by Antonio Frangioni, and Ali Ghezelsoflu
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
#include "UCBlock.h"


/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it {

  class DCNetworkBlock : public NetworkBlock {


/*--------------------------------------------------------------------------*/
/*------------------------- CLASS DCNetworkBlock ---------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/

/** The DCNetworkBlock class derives from NetworkBlock, in order to define a
 *  very basic interface for any possible derived type of dc-network of a
 *  NetworkBlock. Τhe basis of DCNetworkBlock is considered to be very generic
 *  and thus with the minimum possible ingredients. Based on the above
 *  description the class has been constructed having the following elements:
 *
 * - A public method that reads and initializes all the data that describes
 *   the DCNetwork instance.
 *
 * - The Demand Constraint that needs to be satisfied from all the different
 *   units of the UC problem throughout all the time steps of the optimisation
 *   horizon.
 *
 *  A dc-network defined by a set of nodes \f$ \mathcal{N} \f$ and a set
 *  of lines connecting the nodes \f$ \mathcal{L} \f$.
 *
 *   By considering a \f$ |\mathcal{L}| \times |\mathcal{N}| \f$ matrix
 *   \f$ B_t \f$, which constitutes the so-called Power Transfer Distribution
 *   Factor matrix which represents the linear relationship between power
 *   injections at each node of the grid and active power flows through the
 *   transmission lines. The flow limit equations can be written as follow:
 *
 * \f[
 *  P^{mn}_{\ell , t} \leq \sum_{ n' \in \mathcal{N}} (B_t)_({\ell, n'})
 *  (S_{t, n'} - D^{ac}_{n' , t}) \leq
 *  P^{mx}_{\ell , t} \quad t \in \mathcal{T} \quad \ell \in \mathcal{L} \quad
 * \f]
 *   Where \f$ P^{mn}_{\ell , t}\f$ and \f$ P^{mx}_{\ell , t}\f$ are minimum
 *   and maximum power flow at each line \f$ \ell \in \mathcal{L}\f$ and
 *   \f$ S_{t, n'} \f$ and \f$ D^{ac}_{n' , t} \f$ is the node injection
 *   variable and active power demand at node \f$ n \in \mathcal{N} \f$ in the
 *   network respectively.
 */
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
 * - Index, the type of indices;
 @{ */

typedef unsigned int Index;

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

/** Constructor of DCNetworkBlock, taking possibly a pointer of its
 * father Block. */

  DCNetworkBlock( Block * fblock = nullptr ): NetworkBlock( fblock ) ,
     f_NetworkData( nullptr ) , f_local_NetworkData( false ) { }

/*--------------------------------------------------------------------------*/

/// destructor of DCNetworkBlock

 virtual ~DCNetworkBlock() {}

/*@} -----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

/// extends Block::deserialize( netCDF::NcGroup )
/** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
 * the DCNetworkBlock. Besides the mandatory "type" attribute of any :Block,
 * the group should contain the following:
 *
 * - the variable "ActiveDemand", of type double and indexed over the
 *   dimension "NumberNodes"; the i-th entry of the variable is assumed to
 *   contain the active power demand at node i in the network; the variable is
 *   optional and has to be defined when NumberNodes > 1; if NumberNodes == 1
 *   there is no DC-Network then this variable need not be defined, since is
 *   not loaded.
 */

virtual void deserialize( netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/


virtual void generate_abstract_constraints( Configuration *stcc = nullptr )
    override;
/*--------------------------------------------------------------------------*/

virtual void load( std::istream &input ) override { };

/*@} -----------------------------------------------------------------------*/
/*--------------- METHODS FOR MODIFYING THE DCNetworkBlock -----------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the DCNetworkBlock
 *  @{ */

 void set_NetworkData( UCBlock::NetworkData * network_data = nullptr ) override
 {
  // if there was a previous NetworkData and it was local, delete it
  if( f_NetworkData && f_local_NetworkData )
   delete f_NetworkData;

  f_NetworkData = network_data;
  f_local_NetworkData = false;
  }

/*@} -----------------------------------------------------------------------*/
/*--------------------- METHODS FOR SAVING THE DCNetworkBlock --------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the DCNetworkBlock
 *  @{ */

/// extends Block::serialize( netCDF::NcGroup )
/** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
 * DCNetworkBlock. See DCNetworkBlock::deserialize( netCDF::NcGroup ) for
 * details of the format of the created netCDF group.
 */

virtual void serialize( netCDF::NcGroup & group ) const override;

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

 /// the NetworkData object
 UCBlock::NetworkData * f_NetworkData;

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
