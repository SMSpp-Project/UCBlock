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
 * \date 01 - 06 - 2019
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
  class UCBlock;     // forward declaration of UCBlock

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS DCNetworkBlock ---------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
// TODO: short vesion is missing
/** The DCNetworkBlock class derives from NetworkBlock, and defines the
 * standard linear constraints corresponding to the "DC model" of the
 * transmission network in the Unit Commitment problem.
 *
 * TODO: move this part into the comments of
 *       generate_abstract_constraints()
 *
 *  The topology of the transmission network is defined by a set of nodes
 *   \f$ N \f$ and a set of lines connecting the nodes \f$ L \f$.
 *
 *   By considering a \f$ |L| \times |N| \f$ matrix
 *   \f$ B_t \f$, which constitutes the so-called Power Transfer Distribution
 *   Factor matrix which represents the linear relationship between power
 *   injections at each node of the grid and active power flows through the
 *   transmission lines.
 *
 * TODO: you have to drop the index t everywhere, the time instant is fixed 
 *
 * TODO: but how is B defined out of the data of the NetworkBlock??
 *
 * The flow limit equations can be written as follow:
 *
 * \f[
 *  P^{mn}_{\ell , t} \leq \sum_{ n' \in N} (B_t)_({\ell, n'})
 *  (S_{t, n'} - D^{ac}_{n' , t}) \leq
 *  P^{mx}_{\ell , t} \quad t \in \mathcal{T} \quad \ell \in \mathcal{L} \quad
 * \f]
 *   Where \f$ P^{mn}_{\ell , t}\f$ and \f$ P^{mx}_{\ell , t}\f$ are minimum
 *   and maximum power flow at each line \f$ \ell \in L \f$ and
 *   \f$ S_{t, n'} \f$ and \f$ D^{ac}_{n' , t} \f$ is the node injection
 *   variable and active power demand at node \f$ n \in N \f$ in the
 *   network respectively.
 */

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
 * TODO: no, erase. Index is already defined in NetworkBlock, BTW with a
 *       different definition. It makes no sense to have two.
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

/// TODO: short version
/** Constructor of DCNetworkBlock, taking possibly a pointer of its
 * father Block. */

 DCNetworkBlock( Block * fblock = nullptr ) : NetworkBlock( fblock ) ,
  f_NetworkData( nullptr ) , f_local_NetworkData( false ) { }

/*--------------------------------------------------------------------------*/

/// destructor of DCNetworkBlock
// TODO: it must destruct the NetworkData if it is local!!
 virtual ~DCNetworkBlock() {}

/*@} -----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

 // TODO: no, ActiveDemand is already there in NetworkBlock!
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

 void set_NetworkData( NetworkBlock::NetworkData * network_data = nullptr ) override
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
