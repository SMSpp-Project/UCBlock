/*--------------------------------------------------------------------------*/
/*-------------------------- File ECNetworkBlock.h -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 *
 * Header file for the class ECNetworkBlock, which derives from the
 * Block, in order to define the basic interface for the
 * constraints/optimization problems which describe the behaviour of the
 * transmission network in a specific time instant in the Unit Commitment
 * (UC) problem, as represented in UCBlock.
 *
 * Each user is connected to the public grid through each own
 * Point-of-Delivery (PoD), and each user is billed for the energy he
 * consumes and sells.
 * Each instance of this class refers to a specific time period, e.g., a
 * peak period, in the whole time horizon, i.e.,
 * \f$ \mathcal{w} \in \mathcal{W} \f$, and can span an arbitrary number of
 * sub time horizon or intervals, i.e.,
 * \f$ \mathcal{\hat{t}_w} \in \mathcal{\hat{T}_w} \subseteq \mathcal{T} \f$.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy; by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __ECNetworkBlock
 #define __ECNetworkBlock
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "NetworkBlock.h"
#include "FRealObjective.h"
#include "FRowConstraint.h"
#include "OneVarConstraint.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS ECNetworkBlock -----------------------*/
/*--------------------------------------------------------------------------*/
/*----------------------------- GENERAL NOTES ------------------------------*/
/*--------------------------------------------------------------------------*/

class ECNetworkBlock : public NetworkBlock
{

/*--------------------------------------------------------------------------*/
/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/**@} ----------------------------------------------------------------------*/
/*----------------------- CONSTRUCTOR AND DESTRUCTOR -----------------------*/
/*--------------------------------------------------------------------------*/
 /** @name Constructor and Destructor
 * @{ */

 /// Constructor, takes the father
 /** Constructor of ECNetworkBlock, taking possibly a pointer of its father
 * Block. */

 explicit ECNetworkBlock( Block * f_block = nullptr ) :
  NetworkBlock( f_block ) , f_NetworkData( nullptr ) {}

 /// Destructor of ECNetworkBlock

 virtual ~ECNetworkBlock() override;

/// generates the static variables of ECNetworkBlock
/** The base ECNetworkBlock class has just the node injection variables.
 * Since a "bus" network has just one node, and therefore a single value D for
 * the demand and a single injection variable s, which can hardly be called a
 * variable since the only possible way to satisfy the constraints is by
 * having s = D which in fact makes the variable a constant.*/

 void generate_abstract_variables( Configuration * stvv ) override;

/// Generate the static constraint of the ECNetworkBlock
/** This method generates the abstract constraints of the ECNetworkBlock.
 * Since the node injection variable is fixed to the active demand value, it
 * must be a BoxConstraint for that variable whose lower and upper bounds are
 * equal to the active demand value. */

 void generate_abstract_constraints( Configuration * stcc ) override;

/// generate the objective of the ECNetworkBlock
/** Method that generates the objective of the ECNetworkBlock.
 *
 * - Objective function: the objective function of the ECNetworkBlock
 *   is "empty" (a FRealObjective with a LinearFunction inside with no active
 *   variables) */

 void generate_objective( Configuration * objc ) override;

/**@} ----------------------------------------------------------------------*/
/*------- METHODS FOR READING THE DATA OF THE ECNetworkBlock --------*/
/*--------------------------------------------------------------------------*/
 /** @name Reading the data of the ECNetworkBlock
 * @{ */

/// returns the number of nodes
/** Returns the number of nodes in the transmission network. If
 * get_NetworkData() returns nullptr, this is equivalent to
 * get_NetworkData()->get_number_nodes(). Otherwise, it returns zero. */

 Index get_number_nodes() const override {
  if( f_NetworkData )
   return f_NetworkData->get_number_nodes();
  return 0;
 }

/// returns the number of intervals
/** Returns the number of intervals spanned by this NetworkBlock. If
 * get_NetworkData() returns nullptr, this is equivalent to
 * get_NetworkData()->get_number_intervals(). Otherwise, it returns zero. */

 Index get_number_intervals() const override {
  if( f_NetworkData )
   return f_NetworkData->get_number_intervals();
  return 0;
 }

/// returns a pointer to the NetworkData
/** Return a pointer to the NetworkData. */

 NetworkData * get_NetworkData() const override {
  return f_NetworkData;
 }

/// returns the matrix of active demands
/** Method for returning the active demand for the given interval, which is
 * assumed to have size get_number_intervals() per get_number_nodes().
 *
 * @param t The interval wrt the vector of demands for each user is
 *          returned. */

 const double * get_active_demand( Index t = 0 ) const override {
  if( v_active_demand.empty() )
   return nullptr;
  return &( v_active_demand.data()[ t * get_number_nodes() ] );
 }

/// returns the vector of sell prices
/** Method for returning the tariff that user gain to sell electricity to
 * the public market. */

 const std::vector< double > & get_sell_price() const {
  return v_sell_price;
 }

/// returns the vector of buy prices
/** Method for returning the tariff that user pay to buy electricity at each
 * time horizon from the public market. */

 const std::vector< double > & get_buy_price() const {
  return v_buy_price;
 }

/// returns the maximum tariff
/** Method for returning the tariff that user pay due to the peak power. */

 const double & get_max_tariff() const {
  return f_max_tariff;
 }

/**@} ----------------------------------------------------------------------*/
/*---------- METHODS FOR READING THE Variable OF THE ECNetworkBlock --------*/
/*--------------------------------------------------------------------------*/
 /** @name Reading the Variable of the ECNetworkBlock
  *
  * @{ */

 /// returns the vector of micro power injection variables
 /** Method for returning vector of micro power injection variables, which is
  * assumed to have size get_number_nodes(). */

 std::vector< ColVariable > & get_micro_power_injection() {
  return v_micro_power_injection;
 }

 /// returns the vector of micro power absorption variables
 /** Method for returning vector of micro public power absorption variables,
  * which is assumed to have size get_number_nodes(). */

 std::vector< ColVariable > & get_micro_power_absorption() {
  return v_micro_power_absorption;
 }

 /// returns the vector of public power injection variables
 /** Method for returning vector of public power injection variables, which is
  * assumed to have size get_number_nodes(). */

 std::vector< ColVariable > & get_public_power_injection() {
  return v_public_power_injection;
 }

 /// returns the vector of public power absorption variables
 /** Method for returning vector of public power absorbed variables, which is
  * assumed to have size get_number_nodes(). */

 std::vector< ColVariable > & get_public_power_absorption() {
  return v_public_power_absorption;
 }

 /// returns the matrix of node injection variables
 /** Method for returning the node injection for the given interval, which is
  * assumed to have size get_number_intervals() per get_number_nodes().
  *
  * @param t The interval wrt the vector of node injections for each user is
  *          returned. */

 ColVariable * get_node_injection( Index t = 0 ) override {
  if( v_node_injection.empty() )
   return nullptr;
  return &( v_node_injection.data()[ t * get_number_nodes() ] );
 }

/**@} ----------------------------------------------------------------------*/
/*--------- METHODS FOR READING THE Constraint OF THE ECNetworkBlock -------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the Constraint of the ECNetworkBlock
 * @{ */



/**@} ----------------------------------------------------------------------*/
/*--------------- METHODS FOR MODIFYING THE ECNetworkBlock -----------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the DCNetworkBlock
 * @{ */

/// method to set the ActiveDemand
/** This method can be called either before or after that deserialize() is
 * called to provide the NetworkBlock with the ActiveDemand data. This allows
 * all Active Power Demand data corresponding to some UC problem to be
 * "grouped" together (typically, in UCBlock) rather than "spread" among the
 * different NetworkBlock, which may be convenient for some user.
 *
 * If this method is called *before* deserialize(), the data is just copied.
 * However, when deserialize() is called, if ActiveDemand data is present in
 * the NcGroup then this data is used, replacing (and therefore ignoring) the
 * data set by this method.
 *
 * Similarly, if this method is called *after* deserialize(), but some the
 * ActiveDemand was already present in the NcGroup, then that data is kept and
 * the call to this method does nothing.
 *
 * When this method is called, if it is empty it is written into, otherwise
 * nothing happens. In deserialize(), if the data is there in the NcGroup then
 * it is written in v_active_demand (which therefore is no longer empty),
 * otherwise it is left empty so that it can be set by this method. */

 void set_ActiveDemand( const double * v ) override {
  if( v_active_demand.empty() ) {
   // TODO
  }
 }

/**@} ----------------------------------------------------------------------*/
/*------------------------- OTHER INITIALIZATIONS --------------------------*/
/*--------------------------------------------------------------------------*/
 /** @name Other initializations
 * @{ */

 /// deserialize a ECNetworkBlock out of a netCDF::NcGroup
 /** Deserialize a ECNetworkBlock out of a netCDF::NcGroup, which should
  * contain all the data necessary to describe a NetworkBlock (see
  * NetworkBlock::deserialize()).
  * In particular, we refer to that description for the dimension
  * "NumberIntervals". The netCDF::NcGroup must then also contain:
  *
  * - The variable "BuyPrice", of type double and either of size 1 or indexed
  *   over the dimension "NumberIntervals" (if "NumberIntervals" is not
  *   provided, then this variable must be of size 1). This is meant to
  *   represent the vector BuyP[ t ] that, for each time instant t, contains
  *   the tariff that user pay to buy electricity from the public market for
  *   the corresponding time step. If "BuyPrice" has length 1 then BuyP[ t ]
  *   contains the same value for all t.
  *
  * - The variable "SellPrice", of type double and either of size 1 or
  *   indexed over the dimension "NumberIntervals" (if "NumberIntervals" is not
  *   provided, then this variable must be of size 1). This is meant to
  *   represent the vector SellP[ t ] that, for each time instant t, contains
  *   the tariff that user gain to sell electricity to the public market
  *   for the corresponding time step. If "SellPrice" has length 1 then
  *   SellP[ t ] contains the same value for all t.
  * */

 void deserialize( const netCDF::NcGroup & group ) override;

/// loads the ECNetworkBlock instance from an input standard stream.
/** Like load( std::istream & ), if there is any Solver attached to this
 *  ECNetworkBlock then a NBModification (the "nuclear option") is issued.
 *  @warning this method is not implemented yet
 *  @param input an input stream
 *  @param frmt the verbosity level
 */

 void load( std::istream & input , char frmt = 0 ) override {
  throw ( std::logic_error( "ECNetworkBlock::load() not implemented yet" ) );
 }

/** @} ---------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

 void set_active_demand( std::vector< double >::const_iterator values ,
                         Subset && subset ,
                         bool ordered = false ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck ) override final;

 void set_active_demand( std::vector< double >::const_iterator values ,
                         Range rng = Range( 0 , Inf< Index >() ) ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck ) override final;

 static void static_initialization() {
  /*!!
   * Not all C++ compilers enjoy the template wizardry behind the three-args
   * version of register_method<> with the compact MS_*_*::args(), so we just
   * use the slightly less compact one with the explicit argument and be done
   * with it. !!*/
  // register_method< ECNetworkBlock >( "ECNetworkBlock::set_active_demand",
  //                                    &ECNetworkBlock::set_active_demand,
  //                                    MS_dbl_sbst::args() );
  //
  // register_method< ECNetworkBlock >( "ECNetworkBlock::set_active_demand",
  //                                    &ECNetworkBlock::set_active_demand,
  //                                    MS_dbl_rngd::args() );

  register_method< ECNetworkBlock , MF_dbl_it , Subset && , bool >(
   "ECNetworkBlock::set_active_demand" ,
   &ECNetworkBlock::set_active_demand );

  register_method< ECNetworkBlock , MF_dbl_it , Range >(
   "ECNetworkBlock::set_active_demand" ,
   &ECNetworkBlock::set_active_demand );
 }

/*--------------------------------------------------------------------------*/
/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED METHODS OF THE CLASS ---------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS ----------------------------*/
/*--------------------------------------------------------------------------*/

/*---------------------------------- data ----------------------------------*/

 /// the NetworkData object
 NetworkBlock::NetworkData * f_NetworkData;

 /// matrix to store, for each interval, the demand of each node of the network
 boost::multi_array< double , 2 > v_active_demand;

 // energy bought from the public market at the national
 // price /pi^{P-,V} + /pi^{P-,F}
 // (the second term, i.e., the fixed tariff, is given as part of the
 // constant term)

 /// tariff that user pay to buy electricity at each time horizon
 std::vector< double > v_buy_price; // /pi^{P-,V}

 /// tariff that user gain to sell electricity at each time horizon
 std::vector< double > v_sell_price; // /pi^{P+} where /pi^{P+} < /pi^{P-,V}

 /// tariff that user pay due to the peak power
 double f_max_tariff;

/*-------------------------------- variables -------------------------------*/

 /// power injection for each interval at each node
 boost::multi_array< ColVariable , 2 > v_node_injection;

 /// power injected (+) at each node of the network, i.e., at each user PoD,
 /// to the microgrid market / network
 std::vector< ColVariable > v_micro_power_injection; // P^{M^+}

 /// power absorbed (-) at each node of the network, i.e., at each user PoD,
 /// from the microgrid market / network
 std::vector< ColVariable > v_micro_power_absorption; // P^{M^-}

 /// power injected (+) at user PoD from the public market at each time
 /// horizon that is referred to a specific peak period, i.e., a specific
 /// interval in "NumberIntervals"
 std::vector< ColVariable > v_public_power_injection; // P^{P^+}

 /// power absorbed (-) at user PoD from the public market at each time
 /// horizon that is referred to a specific peak period, i.e., a specific
 /// interval in "NumberIntervals"
 std::vector< ColVariable > v_public_power_absorption; // P^{P^-}

 /// maximum power usage at user PoD of the corresponding peak power period,
 /// i.e., a specific interval in "NumberIntervals"
 std::vector< ColVariable > v_max_power; // P^{max}

/*------------------------------- constraints ------------------------------*/

 /// the power balance constraints within the microgrid market / network
 std::vector< FRowConstraint > micro_power_balance_constraints;

 /// the power balance constraints
 boost::multi_array< FRowConstraint , 2 > power_balance_constraints;

 /// the peak power flow limit constraints, i.e., the constraints
 /// on the peak power at user PoD
 boost::multi_array< FRowConstraint , 3 > power_flow_limit_constraints;


 /// the objective function
 FRealObjective objective;

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE FIELDS -------------------------------*/
/*--------------------------------------------------------------------------*/



 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/
/*---------------------------- PRIVATE METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/



}; // end( class( ECNetworkBlock ) )

} /* namespace SMSpp_di_unipi_it */

#endif /* ECNetworkBlock.h included */

/*--------------------------------------------------------------------------*/
/*------------------------ End File ECNetworkBlock.h -----------------------*/
/*--------------------------------------------------------------------------*/