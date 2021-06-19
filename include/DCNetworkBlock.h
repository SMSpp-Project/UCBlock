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
 * \date 18 - 06 - 2021
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
#include "OneVarConstraint.h"
#include "NetworkBlock.h"
#include "FRealObjective.h"

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
 *    is a non-zero value.
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
 virtual ~DCNetworkBlock() override;

/*--------------------------------------------------------------------------*/
/// generate the abstract variables of the DCNetworkBlock
/** Depending on the susceptance for each line of the network, the
 * DCNetworkBlock class may have a power flow variable or not. In other word,
 * if the susceptance is equal to zero(or not defined), the corresponding line
 * is a HVDC line and it must have the power flow variable. It means, each
 * HVDC line correspond to a power flow variable, then for the Net Transfer
 * Capacity (NTC) model all lines must have a power flow variable. If the
 * susceptance value is a non-zero value, the corresponding line is called AC
 * and there is no needed to define the power flow variable for that line.
 * Therefor, in the case of pure AC line there is no needed to define power
 * flow variables. Consequently, for the mixed case AC-HVDC, the power flow
 * variable must define just for HVDC lines. Similarly, depending on the
 * NetworkCost for each line of the network, the DCNetworkBlock class may have
 * an auxiliary variable or not. In other word, if the NetworkCost is equal to
 * zero(or not defined), the auxiliary variable and corresponding constraints
 * will not be defined.
 *
 * Note that since in this class the configuration is ignored.*/

  void generate_abstract_variables( Configuration * stvv = nullptr ) override;

/*--------------------------------------------------------------------------*/
///generate abstract constraints of DCNetworkBlock
/** Three different kinds of DCNetworkBlock constraints are defined as below:
 * The topology of the transmission network is defined by a set of nodes
 * \f$ N \f$ and a set of lines \f$ L \f$. Moreover, it's assumed that
 * \f$ P^{mn}_l \f$ and \f$ P^{mx}_l \f$ are minimum and maximum power flows
 * at each line \f$ l \in L \f$ and \f$ D^{ac}_{n} \f$ is active power demand
 * at node \f$ n \in N \f$ in the network respectively. The node injection
 * variable of each node \f$ n \in N \f$ and the power flows variable and an
 * auxiliary variable(which is not be defined if there is no network cost), of
 * each line \f$ l \in L \f$ are defined as \f$S_{n}\f$, \f$ F_l \f$ and
 * \f$ V_l \f$ respectively.
 *
 *  - DCNetworkBlock with just HVCD lines or the Net Transfer Capacity (NTC)
 *    model:
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
 *      \sum_{l=(n,n') } F_l - \sum_{l=(n',n)} F_l = S_{n} - D_{n}
 *                                                   \quad n \in N   \quad (2)
 *    \f]
 *
 *    Moreover, when NetworkCost for each line is not equal to zero, DCNetwork
 *    will have an objective function which is equal to multiplying NetworkCost
 *    by absolut value of power flows variable. To relaxing the absolute value,
 *    an auxiliary variable and constraints as below are needed:
 *
 *    \f[
 *    F_l \leq V_l                                     \quad l \in L \quad (3)
 *    \f]
 *
 *    \f[
 *    -V_l \leq F_l                                    \quad l \in L \quad (4)
 *    \f]
 *
 *  - DCNetworkBlock with just AC lines model:
 *    By considering a \f$ |L| \times |N| \f$ matrix
 *    \f$ B \f$ which constitutes the so-called Power Transfer Distribution
 *    Factor matrix (PTDF-matrix) which represents the linear relationship
 *    between power injections at each node of the grid and active power flows
 *    through the transmission lines.
 *
 *  The flow limit equations can be written as follow:
 *
 *  \f[
 *   P^{mn}_l\leq \sum_{ n \in N} B_{(l , n)}
 *   (S_n - D^{ac}_n) \leq  P^{mx}_l
 *                                                     \quad l \in L \quad (5)
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
 *       \end{array}\right]                                          \quad (6)
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
 *       \leq  P^{mx}                                                \quad (7)
 *
 *      \f]
 *      where the vector \f$ a = (a_n)_{n = 1, ... , |N| }\f$ and
 *      \f$ b = (b_m)_{m = 1, ... , |L^{dc}| }\f$ are such that for any
 *      \f$ n \in \{ 1, ... , |N|\}\f$ and \f$ m \in \{ 1, ... , |L^{dc}|\}\f$
 *      which \f$ a_n = \sum_{ i \in I_n} p^{ac}_i - D^{ac}_n \f$ and
 *      \f$ b_m = p_{m + |L^{ac}|} = p^{dc}_{\ell(m + |L^{ac}|)}\f$.
 *
 * Note that since in this class the configuration is ignored.*/

 void generate_abstract_constraints( Configuration * stcc = nullptr )
 override;

  /*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// generate the objective of the DCNetworkBlock
/** Method that generates the objective of the DCNetworkBlock.
 *
 * - Objective function: the objective function of the DCNetworkBlock
 *   is given as below:
 *
 *   \f[
 *     \min ( \sum_{ l \in L } ( NC_l V_l )
 *   \f]
 *   where \f$ NC_l \f$, is a network cost and \f$ V_l \f$ is the auxiliary
 *   variable */

  void generate_objective( Configuration *objc ) override;

/**@} ----------------------------------------------------------------------*/
/*---------------- Methods for checking the DCNetworkBlock -----------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for checking solution information in the DCNetworkBlock
 *  @{ */

/*--------------------------------------------------------------------------*/
 /// returns true if the current solution is (approximately) feasible
 /** This function returns true if and only if the solution encoded in the
  * current value of the Variable of this DCNetworkBlock is approximately
  * feasible within the given tolerance. The tolerance can be provided by
  * either \p fsbc or by #f_BlockConfig->f_is_feasible_Congifuration. If \p
  * fsbc is a pointer to a SimpleConfiguration< double >, then the tolerance
  * is the value present in that SimpleConfiguration. If \p fsbc is not a
  * pointer to a SimpleConfiguration< double > but
  * #f_BlockConfig->f_is_feasible_Congifuration is, then the tolerance is the
  * value present in #f_BlockConfig->f_is_feasible_Congifuration. Otherwise,
  * the tolerance is considered 1e-8 by default.
  *
  * This function currently considers only the abstract constraints to
  * determine if the solution is feasible. So, the parameter \p useabstract is
  * currently ignored. Moreover, if no abstract Constraint has been generated,
  * then this method returns true.
  *
  * Each Constraint of this DCNetworkBlock is a RowConstraint and a solution
  * is considered feasible if and only if the relative violation of each
  * RowConstraint of this DCNetworkBlock is not greater than the
  * tolerance. See RowConstraint::rel_viol() for details about the relative
  * violation.
  *
  * @param useabstract This parameter is currently ignored.
  *
  * @param fsbc If it is a pointer to a SimpleConfiguration<double>, then the
  *        value stored in that SimpleConfiguration will be the tolerance that
  *        determines if a solution is feasible. */

 bool is_feasible( bool useabstract = false ,
                   Configuration * fsbc = nullptr ) override;

/**@} ----------------------------------------------------------------------*/
/*---------- METHODS FOR READING THE DATA OF THE DCNetworkBlock ------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the NetworkBlock
    @{ */

 /// returns the number of nodes
 /** Returns the number of nodes in the transmission network. If
  * get_NetworkData() returns nullptr, this is equivalent to
  * get_NetworkData()->get_number_nodes(). Otherwise, it returns zero. */

 Index get_number_nodes( void ) const override {
  if( f_NetworkData )
   return f_NetworkData->get_number_nodes();
  return 0;
 }

/*--------------------------------------------------------------------------*/

 /// returns a pointer to the NetworkData
 /** Return a pointer to the NetworkData. */

 NetworkData * get_NetworkData() const override {
  return f_NetworkData;
  }

/**@} ----------------------------------------------------------------------*/
/*---------- METHODS FOR READING THE Variable OF THE DCNetworkBlock --------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the Variable of the DCNetworkBlock
 *
 * @{ */

 /// returns the vector of power flow variables
 /** The returned std::vector< ColVariable >, say F, contains the power flow
  * variables and is indexed over the dimension number of lines. There are two
  * possible cases:
  *
  * - if F is empty(), then this variable is not defined;
  *
  * - otherwise, F must have f_number_lines rows and F[ l ] is the power flow
  *   variable for line l.
  *   */

 const std::vector< ColVariable > & get_power_flow( ) const {
  return v_power_flow;
 }
/*--------------------------------------------------------------------------*/
  /// returns the vector of auxiliary variables
  /** The returned std::vector< ColVariable >, say V, contains the auxiliary
   * variables and is indexed over the dimension number of lines. There are two
   * possible cases:
   *
   * - if V is empty(), then this variable is not defined;
   *
   * - otherwise, V must have f_number_lines rows and V[ l ] is the auxiliary
   *   variable for line l.
   *   */

  const std::vector< ColVariable > & get_auxiliary_variable( ) const {
   return v_auxiliary_variable;
  }
/**@} ----------------------------------------------------------------------*/
/*--------- METHODS FOR READING THE Constraint OF THE DCNetworkBlock -------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the Constraint of the DCNetworkBlock
 *
 * @{ */

 /// returns the vector of power flow limit constraints
 /** This function returns a const reference to the vector of power flow limit
  * constraints. The i-th element of this vector is a FRowConstraint for the
  * i-th line of the network. */

 const std::vector< FRowConstraint > &
 get_power_flow_limit_constraints( ) const {
  if( ! f_NetworkData )
   throw( std::logic_error( "DCNetworkBlock:get_power_flow_limit_constraints:"
                            " NetworkData has not been set." ) );

  switch( f_NetworkData->get_lines_type() ) {
   case( kAC ):
    return v_AC_power_flow_limit_constraints;
   case( kAC_HVDC ):
   default:
    return v_AC_HVDC_power_flow_limit_constraints;
  }
 }

  const std::vector< BoxConstraint > &
  get_power_flow_limit_HVDC_bounds( ) const {
   if( ! f_NetworkData )
    throw( std::logic_error("DCNetworkBlock:get_power_flow_limit_HVDC_bounds:"
                             " NetworkData has not been set." ) );
     return v_HVDC_power_flow_limit_constraints;
  }

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
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

  void set_active_demand( std::vector< double >::const_iterator values,
                          Subset && subset,
                          const bool ordered = false,
                          c_ModParam issuePMod = eNoBlck,
                          c_ModParam issueAMod = eNoBlck ) final;

  void set_active_demand( std::vector< double >::const_iterator values,
                          Range rng = Range( 0, Inf< Index >() ),
                          c_ModParam issuePMod = eNoBlck,
                          c_ModParam issueAMod = eNoBlck ) final;

  static void static_initialization() {
   /*!!
    * Not all C++ compilers enjoy the template wizardry behing the three-args
    * version of register_method<> with the compact MS_*_*::args(), so we just
    * use the slightly less compact one with the explicit argument and be done
    * with it. !!*/
   // register_method< DCNetworkBlock >( "DCNetworkBlock::set_active_demand",
   //                                    &DCNetworkBlock::set_active_demand,
   //                                    MS_dbl_sbst::args() );
   //
   // register_method< DCNetworkBlock >( "DCNetworkBlock::set_active_demand",
   //                                    &DCNetworkBlock::set_active_demand,
   //                                    MS_dbl_rngd::args() );

   register_method< DCNetworkBlock, MF_dbl_it, Subset &&, const bool >(
    "DCNetworkBlock::set_active_demand",
    &DCNetworkBlock::set_active_demand );

   register_method< DCNetworkBlock, MF_dbl_it, Range >(
    "DCNetworkBlock::set_active_demand",
    &DCNetworkBlock::set_active_demand );
  }

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------data--------------------------------------*/

 /// the NetworkData object
 NetworkBlock::NetworkData * f_NetworkData;

 /// true if the NetworkData object has not been passed from outside
 bool f_local_NetworkData;

/*-----------------------------variables------------------------------------*/

 /// the power flow variables
 std::vector< ColVariable > v_power_flow;

  /// the auxiliary network cost variable
  std::vector< ColVariable > v_auxiliary_variable;

/*----------------------------constraints-----------------------------------*/

 /// AC power flow limit constraints
 std::vector<FRowConstraint> v_AC_power_flow_limit_constraints;

 /// HVDC power flow limit constraints
 std::vector<BoxConstraint> v_HVDC_power_flow_limit_constraints;

 /// AC_HVDC power flow limit constraints
 std::vector<FRowConstraint> v_AC_HVDC_power_flow_limit_constraints;

 /// HVDC power flow and node injection constraints
 std::vector<FRowConstraint> v_power_flow_injection_constraints;

  /// HVDC power flow auxiliary variable 1 constraints
  std::vector<FRowConstraint> v_power_flow_auxiliary_variable_one_constraints;

  /// HVDC power flow auxiliary variable 2 constraints
  std::vector<FRowConstraint> v_power_flow_auxiliary_variable_two_constraints;

 /// AC_HVDC power flow constraints
 std::vector<FRowConstraint> v_AC_HVDC_power_flow_constraints;

  /// the objective function
  FRealObjective objective;

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
  private:

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE FIELDS -------------------------------*/
/*--------------------------------------------------------------------------*/

  SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/

 };  // end( class( DCNetworkBlock ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 }  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* DCNetworkBlock.h included */

/*--------------------------------------------------------------------------*/
/*-------------------- End File DCNetworkBlock.h ---------------------------*/
/*--------------------------------------------------------------------------*/
