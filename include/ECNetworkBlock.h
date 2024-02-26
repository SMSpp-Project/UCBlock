/*--------------------------------------------------------------------------*/
/*-------------------------- File ECNetworkBlock.h -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class ECNetworkBlock, which derives from NetworkBlock
 * and describe the behaviour of the energy community network at a specific time
 * instant or in a time interval, i.e., a peak period that can span an
 * arbitrary number of sub time horizons, in the Unit Commitment problem.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Donato Meoli
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

/// namespace for the Structured Modeling System++ (SMS++)

namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS ECNetworkBlock --------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// an energy community NetworkBlock, i.e., an "EC" network
/**
 * The ECNetworkBlock class derives from NetworkBlock, and embeds the idea
 * that a number of users with no pre-installed electrical generators are
 * joining to create an energy community. Each user is connected to the
 * public grid through each own Point-of-Delivery (PoD), and each user is
 * billed for the energy he consumes and sells. The aggregation of the users,
 * under the umbrella of the formal entity Energy Community, is awarded with an
 * economic benefit that is proportional to the energy shared among the
 * users, which is the energy that in each time step is produced and sold by
 * users within the community. */

class ECNetworkBlock : public NetworkBlock
{

/*--------------------------------------------------------------------------*/
/*------------------------ PUBLIC PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
 *
 * NetworkBlock defines the ECNetworkData data type, a small auxiliary class
 * to bunch the basic data of the community network.
 * @{ */

/*--------------------------------------------------------------------------*/
/*-------------------- CLASS NetworkBlock::NetworkData ---------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
 /// auxiliary class holding basic data about the community network
 /** The ECNetworkData class is a nested sub-class which only serves to have a
  * quick way to load all the basic data that describe the energy community.
  * The rationale is that while often the network does not change during the
  * (short) time horizon of UC, it makes sense to allow for this to happen.
  * This means that individual NetworkBlock objects may in principle have
  * different ECNetworkData, but most often they can share the same. By
  * bunching all the information together we make it easy for this sharing to
  * happen. */

 class ECNetworkData : public NetworkBlock::NetworkData
 {

/*--------------------------------------------------------------------------*/
/*---------------- PUBLIC PART OF THE ECNetworkData CLASS ------------------*/
/*--------------------------------------------------------------------------*/

  public:

/**@} ----------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 * @{ */

  /// constructor of ECNetworkData, does nothing
  ECNetworkData( void ) : NetworkBlock::NetworkData() {}

  /// copy constructor of ECNetworkData, does nothing
  ECNetworkData( NetworkData * ec_network_data ) {}

  /// destructor of ECNetworkData: it is virtual, and empty
  virtual ~ECNetworkData() = default;

/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 * @{ */

  /// deserialize an ECNetworkData out of a netCDF::NcGroup
  /** Deserialize an ECNetworkData out of a netCDF::NcGroup in case the
   * following variables are the same for each ECNetworkBlock of the
   * problem, so they were given just one time in the netCDF, at the head of
   * the hierarchy, which should contain the following:
   *
   * - The dimension "NumberNodes" containing the number of nodes in the
   *   problem; this dimension is mandatory and it cannot be equals to 1
   *   since cannot exists an Energy Community with just one user;
   *
   * - The variable "BuyPrice", of type netCDF::NcDouble and containing the
   *   tariff that the user pays to buy electricity from the public market;
   *
   * - The variable "SellPrice", of type netCDF::NcDouble and containing the
   *   tariff that the user gains to sell electricity to the public market;
   *
   * - The variable "RewardPrice", of type netCDF::NcDouble and containing the
   *   reward benefit awarded to the energy community for the energy consumed
   *   within the community itself;
   *
   * - The variable "PeakTariff", of type netCDF::NcDouble and containing the
   *   tariff that the user pays due to the peak power. */

  virtual void deserialize( const netCDF::NcGroup & group ) override;

/**@} ----------------------------------------------------------------------*/
/*------------ METHODS FOR READING THE DATA OF THE ECNetworkData -----------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the ECNetworkData
 * @{ */

  /// returns the energy sell price
  /** Returns the tariff that the user gains to sell electricity to the
   * public market. */

  double get_sell_price( void ) const { return( f_SellPrice ); }

/*--------------------------------------------------------------------------*/
  /// returns the energy buy price
  /** Returns the tariff that the user pays to buy electricity from the public
   * market. */

  double get_buy_price( void ) const { return( f_BuyPrice ); }

/*--------------------------------------------------------------------------*/
  /// returns the energy reward price
  /** Returns the tariff that the user gains when it absorbs power from the
   * microgrid market / network (instead of from the public grid). */

  double get_reward_price( void ) const { return( f_RewardPrice ); }

/*--------------------------------------------------------------------------*/
  /// returns the peak tariff
  /** Returns the tariff that the user pays due to the peak power. */

  double get_peak_tariff( void ) const { return( f_PeakTariff ); }

/**@} ----------------------------------------------------------------------*/
/*------------------ METHODS FOR SAVING THE ECNetworkData ------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the ECNetworkData
 * @{ */

  /// serialize an ECNetworkData out of a netCDF::NcGroup
  /** Serialize an ECNetworkData out of a netCDF::NcGroup to the specific
   * format of an ECNetworkData. See NetworkBlock::deserialize( netCDF::NcGroup
   * ) for details of the format of the created netCDF group. */

  virtual void serialize( netCDF::NcGroup & group ) const override;

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

  protected:

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED METHODS OF THE CLASS ---------------------*/
/*--------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

  /// tariff that the user pays to buy electricity at each time horizon
  double f_BuyPrice{};

  /// tariff that the user gains to sell electricity at each time horizon
  double f_SellPrice{};

  /// tariff that the user gains when it absorbs power from the microgrid
  /// (instead of from the public grid) at each time horizon
  double f_RewardPrice{};

  /// tariff that the user pays due to the peak power
  double f_PeakTariff{};

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

  private:

/*--------------------------------------------------------------------------*/
/*-------------------- PRIVATE FIELDS OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/



  SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/
/*---------------------- PRIVATE METHODS OF THE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/



 };  // end( class( ECNetworkData ) )

/**@} ----------------------------------------------------------------------*/
/*----------------------- CONSTRUCTOR AND DESTRUCTOR -----------------------*/
/*--------------------------------------------------------------------------*/
 /** @name Constructor and Destructor
 * @{ */

 /// constructor, takes the father
 /** Constructor of ECNetworkBlock, taking possibly a pointer of its father
 * Block. */

 explicit ECNetworkBlock( Block * f_block = nullptr )
  : NetworkBlock( f_block ) , f_NetworkData( nullptr ) {}

/*--------------------------------------------------------------------------*/
 /// destructor of ECNetworkBlock

 virtual ~ECNetworkBlock() override;

/*--------------------------------------------------------------------------*/
 /// generates the static variables of ECNetworkBlock
 /** The size of node injection variable is the number of intervals spanned
  * by this ECNetworkBlock by the number of nodes, which can be read via the
  * NetworkData object (either in the NcGroup or because it has been passed
  * and NumberNodes > 1), so this variable has size "NumberNodes", which can
  * be read via NetworkData::get_number_nodes(). The community scenario also
  * brings with it variables to represent the energy injected (+) or absorbed
  * (-) from both the public grid and the microgrid within the community, and
  * the peak power variables. */

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the static constraint of the ECNetworkBlock
 /** The constraints of an ECNetworkBlock are defined as below:
  *
  * - The power balance w.r.t. the load demand:
  *
  *   \f[
  *    P_{n,t}^{P+} - P_{n,t}^{P-} = S_n - D_n
  *                        \quad n \in \mathcal{N}, t \in \mathcal{T} \quad (1)
  *   \f]
  *
  * - The power dispatch cannot go beyond the peak power at the user PoD, and
  *   is calculated as:
  *
  *   \f[
  *    P_n^{mx} \geq P_{n,t}^{P+} - P_{n,t}^{P-}
  *                       \quad n \in \mathcal{N}, t \in \mathcal{T} \quad (2a)
  *   \f]
  *
  *   \f[
  *    P_n^{mx} \geq - ( P_{n,t}^{P+} - P_{n,t}^{P-} )
  *                       \quad n \in \mathcal{N}, t \in \mathcal{T} \quad (2b)
  *   \f]
  *
  * - The max power shared within the microgrid w.r.t. the public market:
  *
  *   \f[
  *    P_{n,t}^{M} \leq P_{n,t}^{P+}
  *                       \quad n \in \mathcal{N}, t \in \mathcal{T} \quad (3a)
  *   \f]
  *
  *   \f[
  *    P_{n,t}^{M} \leq P_{n,t}^{P-}
  *                       \quad n \in \mathcal{N}, t \in \mathcal{T} \quad (3b)
  *   \f]
  */

 void generate_abstract_constraints( Configuration * stcc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the objective of the ECNetworkBlock
 /** Method that generates the objective of the ECNetworkBlock.
  *
  * - Objective function: the objective function of the ECNetworkBlock
  *   is given as below:
  *
  *   \f[
  *     \min ( \sum_{ n \in \mathcal{N} } \pi^{mx} P_n^{mx} +
  *         \sum_{ t \in \mathcal{T} }
  *         ( \pi_t^{-,v} P_{n,t}^{P-} -
  *         \pi_t^+ P_{n,t}^{P+} -
  *         \pi_t^r P_{n,t}^{M} ) + \pi_t^{-,f} ) )
  *   \f]
  *
  *   where \f$ \pi^{mx} \f$ is the cost due to peak power and \f$ P_n^{mx}
  *   \f$ is the peak power variable; \f$ \pi_t^{-,v} \f$ and
  *   \f$ \pi_t^{-,f} \f$ are the buy prices of the energy bought from the
  *   public market, the variable and fixed costs, i.e., the constant term,
  *   respectively, and \f$ \pi_t^r \f$ is the tariff that user gains
  *   when it absorbs power from the microgrid instead of from the
  *   public market, while \f$ P_{n,t}^{P-} \f$ and \f$ P_{n,t}^{P+} \f$
  *   are the absorption and injection variables form the public market
  *   respectively; \f$ \pi_t^+ \f$ is the sell price of the energy, and
  *   finally \f$ P_{n,t}^{M} \f$ are the energy shard variables into the
  *   microgrid market. */

 void generate_objective( Configuration * objc = nullptr ) override;

/**@} ----------------------------------------------------------------------*/
/*---------------- Methods for checking the ECNetworkBlock -----------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for checking solution information in the ECNetworkBlock
 * @{ */

 /// returns true if the current solution is (approximately) feasible
 /** This function returns true if and only if the solution encoded in the
  * current value of the Variable of this ECNetworkBlock is approximately
  * feasible within the given tolerance. That is, a solution is considered
  * feasible if and only if
  *
  * -# each ColVariable is feasible; and
  *
  * -# the violation of each Constraint of this ECNetworkBlock is not
  *    greater than the tolerance.
  *
  * Every Constraint of this ECNetworkBlock is a RowConstraint and its
  * violation is given by either the relative (see RowConstraint::rel_viol())
  * or the absolute violation (see RowConstraint::abs_viol()), depending on
  * the Configuration that is provided.
  *
  * The tolerance and the type of violation can be provided by either \p fsbc
  * or #f_BlockConfig->f_is_feasible_Configuration and they are determined as
  * follows:
  *
  * - If \p fsbc is not a nullptr and it is a pointer to a
  *   SimpleConfiguration< double >, then the tolerance is the value present
  *   in that SimpleConfiguration and the relative violation is considered.
  *
  * - If \p fsbc is not nullptr and it is a
  *   SimpleConfiguration< std::pair< double , int > >, then the tolerance is
  *   fsbc->f_value.first and the type of violation is determined by
  *   fsbc->f_value.second (any nonzero number for relative violation and
  *   zero for absolute violation);
  *
  * - Otherwise, if both #f_BlockConfig and
  *   f_BlockConfig->f_is_feasible_Configuration are not nullptr and the
  *   latter is a pointer to either a SimpleConfiguration< double > or to a
  *   SimpleConfiguration< std::pair< double , int > >, then the values of the
  *   parameters are obtained analogously as above;
  *
  * - Otherwise, by default, the tolerance is 0 and the relative violation
  *   is considered.
  *
  * This function currently considers only the abstract representation to
  * determine if the solution is feasible. So, the parameter \p useabstract is
  * currently ignored. If no abstract Variable has been generated, then this
  * function returns true. Moreover, if no abstract Constraint has been
  * generated, the solution is considered to be feasible with respect to the
  * set of Variable only. Notice also that, before checking if the solution
  * satisfies a Constraint, the Constraint is computed (Constraint::compute()).
  *
  * @param useabstract This parameter is currently ignored.
  *
  * @param fsbc The pointer to a Configuration that specifies the tolerance
  *             and the type of violation that must be considered. */

 bool is_feasible( bool useabstract = false ,
                   Configuration * fsbc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// returns true if the the energy is shared between users in the community

 bool is_cooperative( void ) {
  if( ! f_NetworkData )
   return( std::any_of( v_RewardPrice.begin() , v_RewardPrice.end() ,
                        []( double cst ) { return( cst != 0 ); } ) );
  return( f_NetworkData->get_reward_price() != 0 );
 }

/**@} ----------------------------------------------------------------------*/
/*----------- METHODS FOR READING THE DATA OF THE ECNetworkBlock -----------*/
/*--------------------------------------------------------------------------*/
 /** @name Reading the data of the ECNetworkBlock
 * @{ */

 /// returns the number of nodes of the network
 /** Returns the number of nodes in the community network. If
  * get_NetworkData() returns nullptr, this is equivalent to
  * get_NetworkData()->get_number_nodes(). Otherwise, it throws an exception
  * since cannot exists an Energy Community with just one user. */

 Index get_number_nodes( void ) const override {
  if( ! f_NetworkData )
   throw( std::invalid_argument( "ECNetworkBlock::get_number_nodes: cannot "
                                 "create an Energy Community with just one "
                                 "user" ) );
  return( f_NetworkData->get_number_nodes() );
 }

/*--------------------------------------------------------------------------*/
 /// returns a pointer to the ECNetworkData
 /** Return a pointer to the ECNetworkData. */

 NetworkData * get_NetworkData( void ) const override {
  return( f_NetworkData );
 }

/*--------------------------------------------------------------------------*/
 /// returns the matrix of active demands
 /** Returns the active demand for the given interval, which is assumed to
  * have size get_number_intervals() by get_number_nodes().
  *
  * @param interval The interval wrt the vector of demands for each user is
  *                 returned.
  */

 const double * get_active_demand( Index interval = 0 ) const override {
  if( v_ActiveDemand.empty() )
   return( nullptr );
  return( &( v_ActiveDemand.data()[ interval * get_number_nodes() ] ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the energy sell price at the given interval
 /** Returns the tariff that the user gains to sell electricity to the
  * public market at the given interval.
  *
  * @param interval The interval wrt the sell price of the energy is returned.
  */

 double get_sell_price( Index interval ) const {
  if( ! f_NetworkData )
   return( v_SellPrice[ interval ] );
  return( f_NetworkData->get_sell_price() );
 }

/*--------------------------------------------------------------------------*/
 /// returns the energy buy price at the given interval
 /** Returns the tariff that the user pays to buy electricity from the public
  * market at the given interval.
  *
  * @param interval The interval wrt the buy price of the energy is returned. */

 double get_buy_price( Index interval ) const {
  if( ! f_NetworkData )
   return( v_BuyPrice[ interval ] );
  return( f_NetworkData->get_buy_price() );
 }

/*--------------------------------------------------------------------------*/

 /// returns the energy reward price at the given interval
 /** Returns the tariff that the user gains when it absorbs power from the
  * microgrid market / network (instead of from the public grid) at the given
  * interval.
  *
  * @param interval The interval wrt the reward price of the energy is returned.
  */

 double get_reward_price( Index interval ) const {
  if( ! f_NetworkData )
   return( v_RewardPrice[ interval ] );
  return( f_NetworkData->get_reward_price() );
 }

/*--------------------------------------------------------------------------*/
 /// returns the peak tariff
 /** Returns the tariff that the user pays due to the peak power. */

 double get_peak_tariff( void ) const {
  if( ! f_NetworkData )
   return( f_PeakTariff );
  return( f_NetworkData->get_peak_tariff() );
 }

/**@} ----------------------------------------------------------------------*/
/*---------- METHODS FOR READING THE Variable OF THE ECNetworkBlock --------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the Variable of the ECNetworkBlock
 * @{ */

 /// returns the vector of public power injection variables
 /** Returns the vector of public power injection variables, which is assumed
  * to have size get_number_nodes(). */

 ColVariable * get_power_injection( Index interval = 0 ) {
  if( v_power_injection.empty() )
   return( nullptr );
  return( &( v_power_injection.data()[ interval * get_number_nodes() ] ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of public power absorption variables
 /** Returns the vector of public power absorbed variables, which is assumed
  * to have size get_number_nodes(). */

 ColVariable * get_power_absorption( Index interval = 0 ) {
  if( v_power_absorption.empty() )
   return( nullptr );
  return( &( v_power_absorption.data()[ interval * get_number_nodes() ] ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of microgrid power variables
 /** The returned std::vector< ColVariable >, say S, contains the energy shared
  * variables and is indexed over get_number_intervals(). There are two
  * possible cases:
  *
  * - if V is empty(), then this variable is not defined;
  *
  * - otherwise, V must have f_numer_intervals rows and S[ i ] is the
  * energy shared within the network at the given interval. */

 const std::vector< ColVariable > & get_shared_power( void ) const {
  return( v_shared_power );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of maximum powers
 /** The returned std::vector< ColVariable >, say V, contains the maximum
  * powers of the corresponding peak power period and is indexed over the
  * dimension f_number_nodes. There are two possible cases:
  *
  * - if V is empty(), then this variable is not defined;
  *
  * - otherwise, V must have f_number_nodes rows and V[ u ] is the
  * maximum peak power for user u. */

 const std::vector< ColVariable > & get_peak_power( void ) const {
  return( v_peak_power );
 }

/**@} ----------------------------------------------------------------------*/
/*--------------- METHODS FOR MODIFYING THE ECNetworkBlock -----------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the ECNetworkBlock
 * @{ */

 void set_NetworkData( NetworkBlock::NetworkData * nd = nullptr ) override {
  // if there was a previous ECNetworkData, and it was local, delete it
  if( f_NetworkData && f_local_NetworkData )
   delete( f_NetworkData );

  f_NetworkData = static_cast< ECNetworkData * >( nd );
  f_local_NetworkData = false;
 }

/*--------------------------------------------------------------------------*/
 /// method to set the number of intervals

 void set_number_intervals( const Index i ) override {
  f_number_intervals = i;
 }

/*--------------------------------------------------------------------------*/
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
  * it is written in v_ActiveDemand (which therefore is no longer empty),
  * otherwise it is left empty so that it can be set by this method. */

 void set_ActiveDemand(
  const std::vector< std::vector< double > > & v ) override {
  if( v_ActiveDemand.empty() ) {
   v_ActiveDemand.resize( boost::multi_array< double , 2 >::extent_gen()
                          [ get_number_intervals() ][ get_number_nodes() ] );
   auto demand = v_ActiveDemand.data();
   for( Index i = 0 ; i < get_number_intervals() ; i++ )
    for( Index j = 0 ; j < get_number_nodes() ; j++ )
     *(demand++) = v[ i ][ j ];
  }
 }

/**@} ----------------------------------------------------------------------*/
/*------------------------- OTHER INITIALIZATIONS --------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
* @{ */

 /// deserialize an ECNetworkBlock out of a netCDF::NcGroup
 /** Deserialize an ECNetworkBlock out of a netCDF::NcGroup in case the
  * following variables are different for each ECNetworkBlock of the problem,
  * so they were explicitly given in each netCDF, each of which should contain
  * the following:
  *
  * - The dimension "NumberIntervals" containing the number of intervals
  *   spanned by this network block; this dimension is optional, if it is
  *   not provided then it is taken to be equal to 1;
  *
  * - The variable "ActiveDemand", of type netCDF::NcDouble and indexed over
  *   the dimensions "NumberIntervals" and "NumberNodes".
  *   If the NetworkData object description is present in the NcGroup this is
  *   the dimension "NumberNodes", but the NetworkData object is optional and it
  *   may not be there. Thus, if "NumberNodes" is not there and "ActiveDemand"
  *   is, then the NetworkData object must have been passed by set_NetworkData(),
  *   and the number of nodes can be read via NetworkData::get_number_nodes().
  *   However, "ActiveDemand" itself is optional. If it is not found in the
  *   NcGroup, then it *must* be passed (either before or after the call to
  *   deserialize()) by calling set_active_demand(). Since both groups of data
  *   are optional, the NcGroup  can actually be empty which implies that all
  *   the data will be (or have been) passed by the in-memory interface. In
  *   this case, it would clearly be preferable to *entirely avoid the NcGroup
  *   to be there*, and in fact UCBlock has provisions for the NcGroup
  *   describing the NetworkBlock to be optional [see the comments to
  *   UCBlock::deserialize()];
  *
  * - The variable "BuyPrice", of type netCDF::NcDouble and either of size 1
  *   or indexed over the dimension "NumberIntervals" (if "NumberIntervals"
  *   is not provided, then this variable must be of size 1). This is meant to
  *   represent the vector BuyP[ t ] that, for each time instant t, contains
  *   the tariff that the user pays to buy electricity from the public market
  *   for the corresponding time step. If "BuyPrice" has length 1 then
  *   BuyP[ t ] contains the same value for all t;
  *
  * - The variable "SellPrice", of type netCDF::NcDouble and either of size
  *   1 or indexed over the dimension "NumberIntervals" (if
  *   "NumberIntervals" is not provided, then this variable must be of size
  *   1). This is meant to represent the vector SellP[ t ] that, for each
  *   time instant t, contains the tariff that the user gains to sell
  *   electricity to the public market for the corresponding time step. If
  *   "SellPrice" has length 1 then SellP[ t ] contains the same value for
  *   all t;
  *
  * - The variable "RewardPrice", of type netCDF::NcDouble and either of size
  *   1 or indexed over the dimension "NumberIntervals" (if
  *   "NumberIntervals" is not provided, then this variable must be of size
  *   1). This is meant to represent the vector RewardP[ t ] that, for each
  *   time instant t, contains the reward benefit awarded to the energy
  *   community for the energy consumed within the community itself, for the
  *   corresponding time step. If "RewardPrice" has length 1 then
  *   RewardP[ t ] contains the same value for all t;
  *
  * - The variable "PeakTariff", of type netCDF::NcDouble and containing the
  *   tariff that the user pays due to the peak power;
  *
  * - The variable "ConstantTerm", of type netCDF::NcDouble and containing the
  *   constant term, i.e., typically the fixed costs. */

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
 /// loads the ECNetworkBlock instance from an input standard stream.
 /** Like load( std::istream & ), if there is any Solver attached to this
  * ECNetworkBlock then a NBModification (the "nuclear option") is issued. */

 void load( std::istream & input , char frmt = 0 ) override {
  throw( std::logic_error( "ECNetworkBlock::load() not implemented yet" ) );
 }

/**@} ----------------------------------------------------------------------*/
/*-------------------- METHODS FOR SAVING THE ECNetworkBlock ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the ECNetworkBlock
 * @{ */

 /// extends Block::serialize( netCDF::NcGroup )
 /** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
  * NetworkBlock. See NetworkBlock::deserialize( netCDF::NcGroup ) for
  * details of the format of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/** @} ---------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

 /// set the active demand at the nodes specified by \p subset
 /** This function sets the active demand at each node in the given \p
  * subset. The active demand at the node whose index is specified by the i-th
  * element in \p subset is given by the i-th element of the vector pointed by
  * \p values, i.e., it is given by the value pointed by (values + i). The
  * parameter \p ordered indicates whether the \p subset is ordered.
  *
  * @param values An iterator to a vector containing the active demand.
  *
  * @param subset The indices of the nodes at which the active demand is being
  *               modified.
  *
  * @param ordered It indicates whether \p subset is ordered.
  *
  * @param issuePMod It controls how physical Modification are issued.
  *
  * @param issueAMod It controls how abstract Modification are issued. */

 void set_active_demand( MF_dbl_it values ,
                         Subset && subset ,
                         const bool ordered = false ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck ) override final;

/*--------------------------------------------------------------------------*/
 /// set the active demand at the nodes specified by \p rng
 /** This function sets the active demand at each node in the given Range \p
  * rng. For each i in the given Range (up to the number of nodes minus 1),
  * the active demand at node i is given by the element of the vector pointed
  * by \p values whose index is (i - rng.first), i.e., it is given by the
  * value pointed by (values + i - rng.first).
  *
  * @param values An iterator to a vector containing the active demand.
  *
  * @param rng A Range containing the indices of the nodes at which the active
  *        demand is being modified.
  *
  * @param issuePMod It controls how physical Modification are issued.
  *
  * @param issueAMod It controls how abstract Modification are issued. */
 void set_active_demand( MF_dbl_it values ,
                         Range rng = Range( 0 , Inf< Index >() ) ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck ) override final;

/*--------------------------------------------------------------------------*/
/*---------------------- PROTECTED PART OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED METHODS OF THE CLASS ---------------------*/
/*--------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

/*---------------------------------- data ----------------------------------*/

 /// the ECNetworkData object
 ECNetworkData * f_NetworkData;

 /// matrix to store, for each interval, the demand of each node of the network
 boost::multi_array< double , 2 > v_ActiveDemand;

 /// tariff that the user pays to buy electricity at each time horizon
 std::vector< double > v_BuyPrice;

 /// tariff that the user gains to sell electricity at each time horizon
 std::vector< double > v_SellPrice;

 /// tariff that the user gains when it absorbs power from the microgrid
 /// (instead of from the public grid) at each time horizon
 std::vector< double > v_RewardPrice;

 /// tariff that the user pays due to the peak power
 double f_PeakTariff{};

/*-------------------------------- variables -------------------------------*/

 /// power injected (+) by the user into the public market at each time
 /// horizon that is referred to a specific peak period, i.e., a specific
 /// interval in "NumberIntervals"
 boost::multi_array< ColVariable , 2 > v_power_injection;

 /// power absorbed (-) by the user from the public market at each time
 /// horizon that is referred to a specific peak period, i.e., a specific
 /// interval in "NumberIntervals"
 boost::multi_array< ColVariable , 2 > v_power_absorption;


 /// power shared into the network
 std::vector< ColVariable > v_shared_power;

 /// maximum power usage by the user at the corresponding peak period, i.e.,
 /// a specific interval in "NumberIntervals"
 std::vector< ColVariable > v_peak_power;

/*------------------------------- constraints ------------------------------*/

 /// the power balance constraints
 boost::multi_array< FRowConstraint , 2 > power_balance_const;

 /// the shared power constraints
 boost::multi_array< FRowConstraint , 2 > power_shared_const;

 /// the peak power flow limit constraints, i.e., the constraints
 /// on the peak power at user PoD
 boost::multi_array< FRowConstraint , 3 > power_flow_limit_const;


 /// the node injection bound constraints
 boost::multi_array< BoxConstraint , 2 > node_injection_bounds_const;


 /// the objective function
 FRealObjective objective;

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------- PRIVATE FIELDS OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/



 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/
/*---------------------- PRIVATE METHODS OF THE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

 static void static_initialization( void ) {

  /* Warning: Not all C++ compilers enjoy the template wizardry behind the
   * three-args version of register_method<> with the compact MS_*_*::args(),
   *
   * register_method< ECNetworkBlock >( "ECNetworkBlock::set_active_demand",
   *                                    &ECNetworkBlock::set_active_demand,
   *                                    MS_dbl_sbst::args() );
   *
   * so we just use the slightly less compact one with the explicit argument
   * and be done with it. */

  register_method< ECNetworkBlock , MF_dbl_it , Subset && , bool >(
   "ECNetworkBlock::set_active_demand" , &ECNetworkBlock::set_active_demand );

  register_method< ECNetworkBlock , MF_dbl_it , Range >(
   "ECNetworkBlock::set_active_demand" , &ECNetworkBlock::set_active_demand );
 }

};  // end( class( ECNetworkBlock ) )

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS ECNetworkBlockMod -------------------------*/
/*--------------------------------------------------------------------------*/

/// derived class from NetworkBlockMod for modifications to an ECNetworkBlock
class ECNetworkBlockMod : public NetworkBlockMod
{

 public:

 /// public enum for the types of NetworkBlockMod
 enum ECNetB_mod_type
 {
  eSetActD = 0 , ///< set active demand values
  eNetBModLastParam  ///< first allowed parameter value for derived classes
  /**< Convenience value to easily allow derived classes to extend the set of
   * types of ECNetworkBlockMod. */
 };

 /// constructor, takes the ECNetworkBlock and the type
 ECNetworkBlockMod( ECNetworkBlock * const fblock , const int type )
  : NetworkBlockMod( fblock , type ) {}

 /// destructor, does nothing
 virtual ~ECNetworkBlockMod() override = default;

 /// returns the Block to which the Modification refers
 Block * get_Block( void ) const override { return( f_Block ); }

 protected:

 /// prints the ECNetworkBlockMod
 void print( std::ostream & output ) const override {
  output << "ECNetworkBlockMod[" << this << "]: ";
  switch( f_type ) {
   default:
    output << "Set active demand values ";
  }
 }

};  // end( class( ECNetworkBlockMod ) )

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS ECNetworkBlockRngdMod -----------------------*/
/*--------------------------------------------------------------------------*/

/// derived from ECNetworkBlockMod for "ranged" modifications
class ECNetworkBlockRngdMod : public ECNetworkBlockMod
{

 public:

 /// constructor: takes the ECNetworkBlock, the type, and the range
 ECNetworkBlockRngdMod( ECNetworkBlock * const fblock , const int type ,
                        Block::Range rng )
  : ECNetworkBlockMod( fblock , type ) , f_rng( rng ) {}

 /// destructor, does nothing
 virtual ~ECNetworkBlockRngdMod() override = default;

 /// accessor to the range
 Block::c_Range & rng( void ) { return( f_rng ); }

 protected:

 /// prints the ECNetworkBlockRngdMod
 void print( std::ostream & output ) const override {
  ECNetworkBlockMod::print( output );
  output << "[ " << f_rng.first << ", " << f_rng.second << " )" << std::endl;
 }

 Block::Range f_rng;  ///< the range

};  // end( class( ECNetworkBlockRngdMod ) )

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS ECNetworkBlockSbstMod ----------------------*/
/*--------------------------------------------------------------------------*/

/// derived from ECNetworkBlockMod for "subset" modifications
class ECNetworkBlockSbstMod : public ECNetworkBlockMod
{

 public:

 /// constructor: takes the ECNetworkBlock, the type, and the subset
 ECNetworkBlockSbstMod( ECNetworkBlock * const fblock , const int type ,
                        Block::Subset && nms )
  : ECNetworkBlockMod( fblock , type ) , f_nms( std::move( nms ) ) {}

 /// destructor, does nothing
 virtual ~ECNetworkBlockSbstMod() override = default;

 /// accessor to the subset
 Block::c_Subset & nms( void ) { return( f_nms ); }

 protected:

 /// prints the ECNetworkBlockSbstMod
 void print( std::ostream & output ) const override {
  ECNetworkBlockMod::print( output );
  output << "(# " << f_nms.size() << ")" << std::endl;
 }

 Block::Subset f_nms;  ///< the subset

};  // end( class( ECNetworkBlockSbstMod ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* ECNetworkBlock.h included */

/*--------------------------------------------------------------------------*/
/*------------------------ End File ECNetworkBlock.h -----------------------*/
/*--------------------------------------------------------------------------*/