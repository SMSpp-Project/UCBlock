/*--------------------------------------------------------------------------*/
/*---------------------- File ECNetworkBlock.h ----------------------*/
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
 * Each time period in the time horizon
 * \f$ \mathcal{\hat{t}_w} \in \mathcal{\hat{T}_w} \subseteq \mathcal{T} \f$
 * refers to an interval \f$ \mathcal{w} \in \mathcal{W} \f$ with a
 * corresponding tariff.
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
  *  @{ */

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
 * equal to the active demand value.
 */

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


/// returns a pointer to the NetworkData
/** Return a pointer to the NetworkData. */

 NetworkData * get_NetworkData() const override {
  return f_NetworkData;
 }

/// returns the matrix of active demands
/** Method for returning the active demand for the given node, which is
 * assumed to have size get_number_intervals() by get_number_nodes(). */

 const double * get_active_demand( Index t = 0 ) const override {
  // TODO
 }

/// returns the matrix of sell prices
/** Method for returning the tariff that user gain to sell electricity to
 * the public market. */

 const std::vector< double > & get_sell_price() const {
  return v_sell_price;
 }

/// returns the vector of buy prices
/** Method for returning the *variable* tariff that user pay to buy electricity
 * at each time horizon from the public market. */

 const std::vector< double > & get_buy_price() const {
  return v_buy_price;
 }

/// returns the matrix of consumption prices
/** Method for returning the *fixed* tariff that user pay to buy electricity
 * at each time horizon from the public market. */

 const std::vector< double > & get_consumption_price() const {
  return v_consumption_price;
 }

/**@} ----------------------------------------------------------------------*/
/*------ METHODS FOR READING THE Variable OF THE ECNetworkBlock -----*/
/*--------------------------------------------------------------------------*/
 /** @name Reading the Variable of the ECNetworkBlock
  *
  * These methods allow to read the just the one set of Variable (which are
  * power injection (+) and absorption (-) to/from the the public market or
  * microgrid market / network, ones for each node) that ECNetworkBlock
  * in necessarily has.
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

/**@} ----------------------------------------------------------------------*/
/*------------------------- OTHER INITIALIZATIONS --------------------------*/
/*--------------------------------------------------------------------------*/
 /** @name Other initializations
  *  @{ */

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
  *   the variable tariff that user pay to buy electricity from the public
  *   market for the corresponding time step. If "BuyPrice" has length 1 then
  *   BuyP[ t ] contains the same value for all t.
  *
  * - The variable "ConsumptionPrice", of type double and either of size 1 or
  *   indexed over the dimension "NumberIntervals" (if "NumberIntervals" is not
  *   provided, then this variable must be of size 1). This is meant to
  *   represent the vector ConsP[ t ] that, for each time instant t, contains
  *   the fixed tariff that user pay to buy electricity from the public market
  *   for the corresponding time step. If "ConsumptionPrice" has length 1 then
  *   ConsP[ t ] contains the same value for all t.
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

/*--------------------------------------------------------------------------*/
/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*---------------------------- PROTECTED FIELDS ----------------------------*/
/*--------------------------------------------------------------------------*/

/*---------------------------------- data ----------------------------------*/

 /// the NetworkData object
 NetworkBlock::NetworkData * f_NetworkData;


 /// vector to store the demand of each node of the network
 boost::multi_array< double , 2 > v_active_demand;

 // energy bought from the public market at the national
 // price /pi^{P-,V} + /pi^{P-,F}

 /// variable tariff that user pay to buy electricity, i.e.,
 /// the tariff on the withdrawing, in any time horizon
 std::vector< double > v_buy_price; // /pi^{P-,V}

 /// fixed tariff that user pay to buy electricity, i.e.,
 /// the tariff component on the demand consumption
 std::vector< double > v_consumption_price; // /pi^{P-,F}

 /// tariff that user gain to sell electricity
 std::vector< double > v_sell_price; // /pi^{P+} where /pi^{P+} < /pi^{P-,V}

/*------------------------------- variables --------------------------------*/

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

/*------------------------------ constraints -------------------------------*/

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

 /// Clear the constraints from any-sized boost::multi_array of FRowConstraint
 template< unsigned long T >
 void clear_constraints(
  boost::multi_array< FRowConstraint , T > & constraints );

}; // end( class( ECNetworkBlock ) )

} /* namespace SMSpp_di_unipi_it */

#endif /* ECNetworkBlock.h included */

/*--------------------------------------------------------------------------*/
/*-------------------- End File ECNetworkBlock.h --------------------*/
/*--------------------------------------------------------------------------*/