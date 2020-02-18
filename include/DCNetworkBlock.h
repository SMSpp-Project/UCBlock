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
 * \date 18 - 02 - 2020
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
 * transmission network in the Unit Commitment problem. Generally, there exist
 * three different kinds of DCNetworkBlock:
 *
 *  - DCNetworkBlock with just HVDC lines; where the susceptance for all lines
 *    is equal to zero. It's also known as the Net Transfer Capacity (NTC)
 *    model.
 *
 *  - DCNetworkBlock with just AC lines; where the susceptance for all lines
 *    is strictly a positive value.
 *
 *  - DCNetworkBlock of an hybrid AC/HVDC grid (both AC and HVDC lines). This
 *    is a combination of first and second cases, where for some lines(not all
 *    of them) may have zero susceptance. */

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
/*--------------------------------------------------------------------------*/
/// generate the abstract variables of the DCNetworkBlock
/** The DCNetworkBlock class in general doesn't have any variable, but in a
 * special case when it corresponds to a model with a single connected grid
 * composed of HVDC lines only, the cass has the power flow variables. This
 * variable is optional, if it is created, its size will be the number of
 * lines. It is also possible to restrict which of the subsets are generated
 * with the parameter stvv. If stvv is not nullptr and it is a
 * SimpleConfiguration<int> or if
 * f_BlockConfig->f_static_variables_Configuration is not nullptr and it is a
 * SimpleConfiguration<int>, then the f_value (an int) indicates whether each
 * of the optional variables should be created. If the Configuration is not
 * available, the default value is taken to be 0. */

  void generate_abstract_variables( Configuration *stvv ) override;

/*--------------------------------------------------------------------------*/
///generate abstract constraints of DCNetworkBlock
/** Three different kinds of DCNetworkBlock constraints are defined as below:
 *  The topology of the transmission network is defined by a set of nodes
 *  \f$ N \f$ and a set of lines \f$ L \f$. Moreover, it's assumed that
 *  \f$ P^{mn}_l \f$ and \f$ P^{mx}_l \f$ are minimum and maximum power flows
 *  at each line \f$ l \in L \f$ and \f$ D^{ac}_{n} \f$ is active power demand
 *  at node \f$ n \in N \f$ in the network respectively. The node injection
 *  variable of each node \f$ n \in N \f$ in and the power flows variable of
 *  each line \f$ l \in L \f$  are defined as \f$S_{n}\f$ and \f$ F_l \f$
 *  respectively.
 *
 *  - The particular case of Net Transfer Capacity (NTC) model:
 *    In this special case the susceptance value for each line is equal to
 *    zero. In fact, this corresponds to a model with a single connected grid
 *    composed of HVDC lines only. In this case, the flow limit equations
 *    define as:
 *
 *    \f[
 *    P^{mn}_l \leq F_l  \leq P^{mx}_l
 *                                                     \quad l \in L \quad (1)
 *    \f]
 *
 *   where in each line \f$ l \f$, \f$ n \f$ and \f$ n' \f$ are supposed to
 *   be the start and the end point of that respectively. Besides, the
 *   following link between power flows and injected power at each node of the
 *   grid:
 *
 *    \f[
 *      \sum_{l=(n',n) } F_l - \sum_{l=(n,n')} F_l) = S_{n}
 *                                                   \quad n \in N   \quad (2)
 *    \f]
 *
 *  - DCNetworkBlock with just AC lines model:
 *    By considering a \f$ |L| \times |N| \f$ matrix
 *    \f$ B_t \f$, which constitutes the so-called Power Transfer Distribution
 *    Factor matrix (PTDF-matrix) which represents the linear relationship
 *    between power injections at each node of the grid and active power flows
 *    through the transmission lines.
 *
 *  The flow limit equations can be written as follow:
 *
 *  \f[
 *   P^{mn}_l\leq \sum_{ n \in N} (B)_({l , n})
 *   (S_{ n} - D^{ac}_{n}) \leq  P^{mx}_l
 *                                                     \quad l \in L \quad (3)
 *  \f]
 *
 *  - DCNetworkBlock of an hybrid AC/HVDC grid (both AC and HVDC lines):
 *    This is the case of an hybrid grid constituted of both AC and HVDC (High
 *    Voltage Direct Current) lines. The DC lines are characterized by the
 *    fact that the flow passing through those lines is fully controllable.
 *    However, this flow still has an impact on the flows passing through
 *    connected AC lines. To use the matrix formalism, we first introduce some
 *    additional notations:
 *
 *    - Lines of the grid are indexed by \f$ l = 1,···|L^{ac}| \f$ for AC
 *      lines, while indexes \f$ l = |L^{ac}| + 1,··· ,|L^{ac}|+|L^{dc}| \f$
 *      refer to DC lines.
 *
 *    - For any \f$ l \in \{1,···|L^{ac}| \}\f$ and
 *      \f$ k \in \{1,···|L^{dc}| \} \f$, and put \f$ \ell(l) \f$ the pair of
 *      nodes related by the AC line indexed by \f$ l \f$ and
 *      \f$ \ell(k+|L^{ac}|) \f$ denotes the pair of nodes related by the DC
 *      line indexed by \f$ k+|L^{dc}| \f$.
 *
 *    - \f$ A^{dc} \f$ denotes the \f$ |L^{dc}| \times |N| \f$ incidence
 *      matrix induced by DC lines of the grid and the
 *      \f$ |L| \times (|L^{dc}| + |N|) \f$ matrix of A, where obtained by
 *      concatenation of bloc matrices as follows:
 *
 *     \f[
 *
 *       A = \left[
 *       \begin{array}{cc}
 *       B & -B(A^{dc})^T \\
 *       0_{|L^{dc}| \times |N|} & I_{|L^{dc}| \times |L^{dc}|}
 *       \end{array}\right]                                          \quad (4)
 *
 *     \f]
 *
 *      Where \f$ 0_{|L^{dc}| \times |N|} \f$  denotes the
 *      \f$ |L^{dc}| \times |N| \f$ zero matrix and
 *      \f$ I_{|L^{dc}| \times |L^{dc}|}\f$ the \f$|L^{dc}| \times |L^{dc}|\f$
 *      identity matrix. Therefor, the flow limit equations can be transformed
 *      into:
 *
 *      \f[
 *
 *       P^{mn} \leq A \left[
 *       \begin{array}{c}
 *       a \\
 *       b
 *       \end{array}\right]
 *       \leq  P^{mx}                                                \quad (5)
 *
 *      \f]
 *      where the vector \f$ a = (a_n)_{n = 1, ... , |N| }\f$ and
 *      \f$ b = (b_m)_{m = 1, ... , |l^{dc}| }\f$ are such that for any
 *      \f$ n \in \{ 1, ... , |N|\}\f$ and \f$ m \in \{ 1, ... , |l^{dc}|\}\f$
 *      which \f$ a_n = \sum_{ i \in I_n} p^{ac}_i - D^{ac}_n \f$ and
 *      \f$ b_m = p_{m + |L^{ac}|} = p^{dc}_{\ell(m + |L^{ac}|)}\f$.
 *
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

  /// DC flow limit constraints
  std::vector<FRowConstraint> v_AC_flow_limit_constraints;

  /// HVDC power flow limit constraints
  std::vector<FRowConstraint> v_HVDC_flow_limit_constraints;

  /// HVDC power flow and node injection constraints
  std::vector<FRowConstraint> v_flow_injection_constraints;

  /// AC_HVDC power flow constraints
  std::vector<FRowConstraint> v_AC_HVDC_flow_constraints;

  };   // end( class( DCNetworkBlock ) )

/*-----------------------------variables------------------------------------*/
 /// the power flow variables
 std::vector< ColVariable > v_power_flow;

/*--------------------------------------------------------------------------*/

}  /* namespace SMSpp_di_unipi_it */

/*--------------------------------------------------------------------------*/

#endif /* DCNetworkBlock.h included */

/*--------------------------------------------------------------------------*/
/*------------------------ End File DCNetworkBlock.h -----------------------*/
/*--------------------------------------------------------------------------*/
