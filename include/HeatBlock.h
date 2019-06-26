/*--------------------------------------------------------------------------*/
/*---------------------------- File HeatBlock.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 *
 * Header file for cass HeatBlock, which derives from Block, in order to
 * define a class representing set of "nearby" units (and a storage) that
 * can be used to satisfy a demand for (a single type of) heat.
 *
 * \version 0.11
 *
 * \date 19 - 06 - 2019
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
 * Copyright &copy; by Antonio Frangioni, Ali Ghezelsoflu, and Rafael
 * Durbano Lobato
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __HeatBlock
 #define __HeatBlock  /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Block.h"
#include "ColVariable.h"
#include "FRowConstraint.h"
#include "FRealObjective.h"

/*--------------------------------------------------------------------------*/
/*------------------------------ NAMESPACE ---------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)

namespace SMSpp_di_unipi_it {

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS HeatBlock -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/

/// implementation of the Block concept for a set of "heat units"
/** The HeatBlock class implements the Block concept [see Block.h] for a
 * set of "heat units" as required in the plan4res project.
 *
 * The constraints regarding heat management in the plan4res project can be
 * described separately for each HeatBlock. In the parlance of plan4res, a HB
 * is the intersection of an Energy Cell and a heat-ID: a set of
 * heat-producing units that are geographically close enough to exchange heat
 * of the same type, together with possibly a (single) storage for this type
 * of heat, in order to satisfy a given heat demand. HBs are only connected to
 * each other because some heat-producing units are also electricity-producing
 * ones, where heat is a by-product or a co-product of electricity. However,
 * this linking “happens at the level of the main UC model”, and therefore it
 * is not required for the description of the internal constraints of a HB.
 * All this on a given time horizon, as in the UC problem. The operations of
 * the heat generating unit are described on a discrete time horizon which in
 * this description we indicate it with \f$ T \f$.
 *
 * The HB therefore has as primary data the description of a set \f$ I \f$
 * of heat-producing units, possibly of a single heat storage, and of the
 * demand that has to be satisfied.
 *
 * TODO: move the above in the comments togenerate_abstract_variables(),
 *       generate_abstract_constraints() and generate_objective()
 *
 * The variables of the HB are the following:
 *
 * - \f$ p^{he}_t \f$ : representing the power produced by heat unit
 *   \f$ i \in I \f$ at time \f$ t \in T \f$;
 *
 * - if the heat storage is defined, \f$ s_{t,+} \geq 0 \f$ and \f$ s_{t,-}
 *   \geq 0 \f$ representing respectively the amount of heat added to and
 *   removed from the storage at time instant \f$ t \in T \f$;
 *
 * - if the heat storage is defined, \f$ v_t \f$ representing the amount of
 *   heat available in the storage at time instant \f$ t \in T \f$.
 *
 *  The constraints in the HB write as follow:
 *
 * - Demand Constraints: with \f$ D_t \f$ denoting the heat demand of the HB
 *   at time period \f$ t \in T \f$:
 *   \f[
 *     \sum_{ i \in I } ( p^{he}_{t,i} - s^{h}_{t,+} + s^{h}_{t,-}
       \geq D_t                     \quad t \in T     \quad     (1)
 *   \f]
 *
 * - Heat production bounds Constraints: with \f$ P^{mn}_{t,i} \f$ and
 *   \f$ P^{mx}_{t,i} \f$ denoting respectively the minimum and maximum heat
 *   production of unit \f$ i \in I \f$ at time \f$ t \in T \f$, the
 *   heat production bounds are
 *   \f[
 *     P^{mn}_{t,i} \leq p^{he}_{t,i} \leq P^{mx}_{t,i}
 *         \quad i \in I    \quad t \in T   \quad     (2)
 *   \f]
 *
 * - Heat storage bounds Constraints: with \f$ V^{mn}_t \f$ and
 *   \f$ V^{mx}_t \f$ denoting respectively the minimum and maximum heat
 *   storage . For each heat block at time \f$ t \in T \f$, the heat
 *   storage bounds are
 *   \f[
 *     v^{mn}_{t} v_t \leq V^{mx}_t    \quad t \in T    \quad   (3)
 *   \f]
 *
 * - Evolution in the stored heat Constraints. Let three constants
 *   \f$ \rho_+ \f$, \f$ \rho_- \f$, and \f$ \rho \f$ be given representing
 *   inefficiencies in, respectively, storing heat in the heat storage,
 *   extracting heat from the heat storage, and keeping heat in the heat
 *   storage; then the evolution in the stored heat can be written as
 *   \f[
 *    v_t = \rho v_{t-1} + \rho_+ s_{t,+} - \rho_- s^{h}_{t,-}
 *                                   \quad t \in T      \quad  (4)
 *   \f]
 *
 * - Objective Function: the objective function of HB simply reads
 *   \f[
 *     \min \sum_{ i \in I} \sum_{ t \in T}
 *          C_{t,i} p^{he}_{t,i}
 *   \f]
 *   where \f$  C_{t,i} \f$ is the cost of producing one heat unit by unit
 *   \f$ i \in I \f$ at time \f$ t \in T \f$. Note that
 *   storing heat has no cost.
 */

class HeatBlock : public Block {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------- PUBLIC TYPES OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
 *
 * HeatBlock defines the following main public types:
 *
 * - Index, the type of parameters indices;
 *
 * @{ */

/*--------------------------------------------------------------------------*/

 typedef std::size_t Index;

/**@} ----------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

/// constructor, takes the father and the time horizon
/** Constructor of HeatBlock, taking possibly a pointer of its father
 * Block and the time horizon. */

 HeatBlock( Block * father_block = nullptr , Index t = 0 )
  : Block( father_block )c, f_time_horizon( t ) { }

/*--------------------------------------------------------------------------*/
/// destructor of HeatBlock: it is virtual, and empty

 virtual ~HeatBlock() { }

/*@}------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

/// loads the HeatBlock instance from a stream

 virtual void load( std::istream &input ) override {
  throw( std::logic_error( "HeatBlock::load() not implemented yet" ) );
  };

/*--------------------------------------------------------------------------*/
/// extends Block::deserialize( netCDF::NcGroup )
/** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
 * the HeatBlock. Besides the mandatory "type" attribute of any :Block,
 * the group should contain the following:
 *
 * - The dimension "NumberHeatUnits" containing the number heat-producing
 *   units. This is not optional, and has to be >= 1.
 *
 * - The dimension "TimeHorizon" containing the time horizon. The dimension
 *   is optional because the same information may be passed via the method
 *   set_time_horizon(), or directly retrieved from the father if it is a
 *   UCBlock; see the comments to set_time_horizon() for details.
 *
 * - The dimension "NumberIntervals", that is provided to allow that all
 *   time-dependent data in the HeatBlock can only change at a subset of
 *   the time time instants of the time interval, being therefore
 *   piecewise-constant (possibly, constant). "NumberIntervals" should
 *   therefore be <= "TimeHorizon", with three distinct cases:
 *  
 *    i)  "NumberIntervals" <= 1, which is taken to mean "NumberIntervals"
 *        == 1; this is what is assumed if the dimension, that is optional,
 *        is not there. This means that the value of each  relevant data in
 *        the HeatBlock (see e.g. "CostHeatUnit", "MinHeatProduction" and
 *        "MaxHeatProduction" below) is the same for each time instant
 *        0, ..., "TimeHorizon" - 1 in the time horizon. In this case, also
 *        the variable "ChangeIntervals" (see below) is ignored.
 *
 *   ii)  1 < "NumberIntervals" < "TimeHorizon", which means that in some
 *        time instants, *but not all of them*, the values of some of the
 *        relevant data are changing; the intervals are then described in
 *        variable "ChangeIntervals".
 *
 *   iii) "NumberIntervals" == "TimeHorizon",  which means that values of 
 *        the relevant data changes at every time interval (in principle;
 *	  of course there is nothing preventing the same value to be
 *        repeated in the netCDF input). Also in this case the variable
 *        "ChangeIntervals" is ignored, since it is useless.
 *
 *   Note that this (together with "ChangeIntervals", if defined) obviously
 *   sets the "maximum frequency" at which data can change; is some data
 *   changes less frequently (say, it is constant), then the same value
 *   will have to be repeated. Individual data can also have specific
 *   provisions for the case where the data is all equal despite
 *   "NumberIntervals" saying differently, see "CostHeatUnit" etc.
 *   for instances.
 *
 * - The variable "ChangeIntervals", of type integer and indexed over the
 *   dimension "NumberIntervals". The time horizon is subdivided into
 *   NumberIntervals = k of the form [ 0 , i_1 ], [ i_1 + 1 , i_2 ], ...
 *   [ i_{k-1} + 1 , "TimeHorizon" - 1 ]; "ChangeIntervals" then has to
 *   contain [ i_1 , i_2 , ... i_{k-1} ]. Note that, therefore,
 *   "ChangeIntervals" has one significant value less than
 *   "NumberIntervals", which means that
 *   ChangeIntervals[ NumberIntervals - 1 ] is ignored. Anyway, the whole
 *   variable is ignored if either "NumberIntervals" <= 1 (such as if it
 *   is not defined), or "NumberIntervals" >= "TimeHorizon".
 * 
 * - The variable "TotalHeatDemand", of type double and indexed over the
 *   dimension "TimeHorizon": entry HeatDemand[ t ] is assumed to contain
 *   the total heat demand of this heat block to be satisfied for the
 *   corresponding time instant t. Note that this variable does *not*
 *   use the NumberIntervals / ChangeIntervals system, as demand is
 *   usually changing from one time instant to the next.
 *
 * - The variable "CostHeatUnit", of type double double and indexed over two
 *   dimensions. The first dimension can have size 1 or "NumberIntervals".
 *   The second dimension has size "NumberHeatUnits". This is meant to
 *   represent the matrix CHU[ t , i ] which is assumed to contain the
 *   unitary cost of heat production of heat-producing unit i for time
 *   instant t. If the first dimension has size 1, then the cost is the
 *   same for all time instants (for the same unit). Otherwise,
 *   CostHeatUnit[ h , i ] is the fixed value of FC[ t , i ] for all t in
 *   the interval [ ChangeIntervals[ h - 1 ], ChangeIntervals[ h ] ],
 *   with the assumption that ChangeIntervals[ - 1 ] = 0.
 *
 * - The variable "MinHeatProduction", of type double double and indexed over
 *   two dimensions. The first dimension can have size 1 or
 *   "NumberIntervals". The second dimension has size "NumberHeatUnits". This
 *   is meant to represent the matrix MnHP[ t , i ] which is assumed to
 *   contain the minimum heat production of heat-producing unit i for time
 *   instant t. The variable is optional: if it is not provided at all, it
 *   is intended MnHP[ t , i ] == 0 for all i and t. If the first dimension
 *   has size 1, then the minimum heat production is the same for all time
 *   instants (for the same unit). Otherwise, MinHeatProduction[ h , i ] is
 *   the fixed value of MnHP[ t , i ] for all t in the interval
 *   [ ChangeIntervals[ h - 1 ], ChangeIntervals[ h ] ], with the assumption
 *   that ChangeIntervals[ - 1 ] = 0.
 *
 * - The variable "MaxHeatProduction", of type double double and indexed over
 *   two dimensions. The first dimension can have size 1 or
 *   "NumberIntervals". The second dimension has size "NumberHeatUnits". This
 *   is meant to represent the matrix MxHP[ t , i ] which is assumed to
 *   contain the maximum heat production of heat-producing unit i for time
 *   instant t. It is assumed MxHP[ t , i ] >= MnHP[ t , i ] >= 0 for all
 *   i and t, with strict inequality holding for at least some t for each
 *   unit i (otherwise the production of unit i is fixed and there is
 *   nothing to decide).
 *
 *   NOTE: I don't think it makes sense that the maximum heat production is
 *         0, it would mean no heat ever, what would be the point?
 *
 *   If the first dimension has size 1, then the maximum heat production is
 *   the same for all time instants (for the same unit). Otherwise,
 *   MaxHeatProduction[ h , i ] is  the fixed value of MxHP[ t , i ] for all
 *   t in the interval [ ChangeIntervals[ h - 1 ], ChangeIntervals[ h ] ],
 *   with the assumption that ChangeIntervals[ - 1 ] = 0.
 *
 * - The variable "MinHeatStorage", of type double and either indexed over
 *   the dimension "NumberIntervals" or has size 1. This is meant to represent
 *   the vector MnHS[ t ] which, for each time instant t, contains the
 *   minimum heat storage of this HB. The variable is optional, if it is not
 *   provided at all it is intended that MnHS[ t ] == 0 for all t. If the
 *   variable has size 1, then the minimum heat storage is the same for all
 *   time instants. Otherwise, MinHeatStorage[ h ] is the fixed value of
 *   MnHS[ t ] for all t in the interval [ ChangeIntervals[ h - 1 ],
 *   ChangeIntervals[ h ] ], with the assumption that
 *   ChangeIntervals[ - 1 ] = 0.
 *
 * - The variable "MaxHeatStorage", of type double and either indexed over
 *   the dimension "NumberIntervals" or has size 1. This is meant to represent
 *   the vector MxHS[ t ] which, for each time instant t, contains the
 *   maximum heat storage of this HB. The variable is optional, if it is not
 *   provided at all it is intended that MxHS[ t ] == 0 for all t. Since it
 *   is assumed that MxHS[ t ] >= MnHS[ t ] >= 0 for all t, this means that
 *   there is no heat storage in this HB. If the variable has size 1, then
 *   the maximum heat storage is the same for all time instants. Otherwise,
 *   MaxHeatStorage[ h ] is the fixed value of MxHS[ t ] for all t in the
 *   interval [ ChangeIntervals[ h - 1 ], ChangeIntervals[ h ] ], with the
 *   assumption that ChangeIntervals[ - 1 ] = 0.
 *
 * - The scalar variable "InitialHeatAvailable", of type double and not
 *   indexed over any dimension, which indicates the the amount of heat in
 *   the storage at the beginning of the first time instant in this HB. It
 *   is assumed that MxHS[ 0 ] >= InitialHeatAvailable >= MnHS[ 0 ]. This
 *   variable is optional, if it is not provided it is taken to be
 *   InitialHeatAvailable == MnHS[ 0 ]. If there is no heat storage (say,
 *   MaxHeatStorage is not defined) then this variable is not read, because
 *   it is not used.
 *
 * - The scalar variable "StoringHeatRho", of type double and not indexed
 *   over any dimension, which indicates the inefficiency of storing heat
 *   in the heat storage (if any) in this HB. This variable is optional and
 *   it must always be StoringHeatRho <= 1, if it is not provided it is taken
 *   to be StoringHeatRho == 1. If there is no heat storage (say,
 *   MaxHeatStorage is not defined) then this variable is not read, because
 *   it is not used.
 *
 * - The scalar variable "ExtractingHeatRho", of type double and not indexed
 *   over any dimension, and which indicates the inefficiency of extracting
 *   heat from the heat storage (if any) in this HB. This variable is
 *   optional and it must always be ExtractingHeatRho >= 1, if it is not
 *   provided it is taken to be ExtractingHeatRho == 1. If there is no heat
 *   storage (say, MaxHeatStorage is not defined) then this variable is not
 *   read, because it is not used.
 *
 * - The scalar variable "KeepingHeatRho", of type double and not indexed over
 *   any dimension, which indicates the double of keeping heat in the heat
 *   storage (if any) in this HB. This variable is optional and it must always
 *   be KeepingHeatRho <= 1, if it is not provided it is taken to be
 *   KeepingHeatRho == 1. If there is no heat storage (say, MaxHeatStorage is
 *   not defined) then this variable is not read, because it is not used. */

 virtual void deserialize( netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
/// generate the abstract variables of the HeatUnit
/** Method that generates the abstract variables of the HeatBlock.
 *
 * HeatBlock class has four different "groups" of variables:
 * 
 * - the heat added variables
 *
 * - the heat removed variables
 *
 * - the heat available variables
 *
 * - the heat variables
 *
 * All of these variables are optional, except the heat variables, in the
 * sense that the model may just not have them (say, because the problem may
 * have not heat added, or it cannot generate heat available value). However,
 * it is also possible to restrict which of the subsets are generated with the
 * parameter stvv.
 *
 * If stvv is not nullptr and it is a SimpleConfiguration<int>, or if
 * f_BlockConfig->f_static_variables_Configuration is not nullptr and it is a
 * SimpleConfiguration<int>, then the f_value (an int) indicates whether each
 * of the optional variables should be created. If the Configuration is not
 * available, the default value is taken to be 0. The value of the int is
 * interpreted bit-wise, with heat added variables being bit 0, heat removed
 * reserve variables being bit 1, heat available reserve variables being bit
 * 2, and heat variables being bit 3. If the bit associated with a variable is
 * 0 then the variable (assuming the model actually has it) *is* created,
 * otherwise it is *not*; hence, the default value of 0 means that all the
 * variables (that the model has) are created.
 *
 * Whenever a group of variables is created, its size will be the time
 * horizon.
 *
 * - if f_time_horizon > 0, a std::vector< ColVariable > with
 *   exactly t entries, the entry t = 0, ..., f_time_horizon - 1 corresponding
 *   to the heat added variables;
 *
 * - if f_time_horizon > 0, a std::vector< ColVariable > with
 *   exactly t entries, the entry t = 0, ..., f_time_horizon - 1 corresponding
 *   to the heat removed variables;
 *
 * - if f_time_horizon > 0, a std::vector< ColVariable > with
 *   exactly t entries, the entry t = 0, ..., f_time_horizon - 1 corresponding
 *   to the heat available variables;
 *
 * - if f_time_horizon > 0, a std::vector< ColVariable > with
 *   exactly t entries, the entry t = 0, ..., f_time_horizon - 1 corresponding
 *   to the heat variables;
 *
 * Note that derived classes are free to use the other bits of the int to
 * similarly encode for creation of their own specific groups of variables.*/

 virtual void generate_abstract_variables( Configuration *stvv = nullptr )
    override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// generate the static constraint of the HeatBlock
/** Method that generates the static constraint of the HeatBlock.
 *
 * TODO: add comments moving out from the comments to the class */

 virtual void generate_abstract_constraints( Configuration *stcc = nullptr )
  override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// generate the objective of the HeatBlock
/** Method that generates the objective of the HeatBlock.
 *
 * TODO: add comments moving out from the comments to the class */

 virtual void generate_objective( Configuration *objc = nullptr )  override;

/*@} -----------------------------------------------------------------------*/
/*-------------- Methods for reading the data of the HeatBlock -------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for reading the data of the HeatBlock
 *  @{  */

/// returns the time horizon of the problem
 Index get_time_horizon( void ) const { return( f_time_horizon ); }

/// returns the number of units in this HeatBlock
 Index get_number_heat_units( void ) const { return( f_number_heat_units ); }

 // TODO: I do not think these implementations are correct, because the
 //       vectors may have 1st dimension == 1. In this case you should
 //       not duplicate the vectors for all time instants
 //
 //       Anyway, you should not have "interval" as first argument, but
 //       rather "time". This is different.
 //
 //       Finally, you are using std::vector<> to represent matrices, but
 //       we have boost::multi_array<> to do the same thing in a more
 //       natural way, you should consider using these (although the final
 //       choice is yours)

 /** returns the minimum heat production of the given block for interval t of
 * unit i */
 inline double get_min_heat_production( Index interval , Index unit ) const {
  return v_min_heat_production[ interval * f_number_heat_units + unit ];
  }

/** returns the cost of heat unit of the given block for interval t of
 * unit i */
    inline double get_cost_heat_unit(  Index interval , Index unit) const {
        return v_cost_heat_unit[ interval * f_number_heat_units + unit ];
    }

/** returns the maximum heat production of the given block for interval t of
 * unit i */
    inline double get_max_heat_production(  Index interval , Index unit) const {
        return v_max_heat_production[ interval * f_number_heat_units + unit ];
    }

/**@} ----------------------------------------------------------------------*/
/*------------- METHODS FOR READING THE Variable OF THE HeatBlock ----------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the Variable of the HeatBlock
 *
 * These methods allow to read the only Variable of the HeatBlock that must
 * be "known" outside of it, i.e., the heat variables.
 *
 * TODO: the second method is surely wrong, for each t there is a separate
 *       heat variable for each unit i
 *
 * TODO: you should use boost::multi_array< ColVariable > for
 *       multi-dimensional groups of Variable like the heat ones, which
 *       would change the signature of the first method
 *
 * @{ */

 /// Method for returning the vector of heat variables

 const std::vector<ColVariable> & get_heat( void ) const {
  return( v_heat );
  }

 /// Method for returning the pointer to the heat variable at time t

 ColVariable * get_heat( int t ) { return & ( v_heat[ t ] ); }

/**@} ----------------------------------------------------------------------*/
/*--------------------- METHODS FOR SAVING THE HeatBlock -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the HeatBlock
 *  @{ */

/// extends Block::serialize( netCDF::NcGroup )
/** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
 * HeatBlock. See HeatBlock::deserialize( netCDF::NcGroup ) for
 * details of the format of the created netCDF group. */

 virtual void serialize( netCDF::NcGroup & group ) const override;

/*@} -----------------------------------------------------------------------*/
/*----------------- METHODS FOR MODIFYING THE HeatBlock --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the HeatBlock
 *  @{ */

 /// set the time horizon method
 /** Set the time horizon method.
  *
  * This method can be called *before* that deserialize() is called to
  * provide the HeatBlock with the time horizon. This allows the information
  * not to be duplicated in the netCDF group that describes the HB, since
  * usually (bit not necessarily) a HeatBlock is deserialized inside a
  * UCBlock, and all HB have the same time horizon, that can therefore be
  * read once and for all by the father UCBlock.
  *
  * If this method is *not* called, which means that f_time_horizon is at its
  * initial value of 0 (not initialized), then when deserialize() is called
  * the information has to be available by other means, i.e.:
  *
  * (i)  If there is no dimension TimeHorizon in netCDF input, then the
  *      HeatBlock must have a father, which must be a UCBlock: the
  *      time horizon is then taken to be that of the father. If the
  *      HeatBlock does not have a father (or it is not a UCBlock), then
  *      exception is thrown.
  *
  * (ii) If the dimension TimeHorizon is present in the netCDF input of
  *      HeatBlock, the value provided there is used with no check that
  *      the UCBlock has a father at all, that the father is a UCBlock,
  *      or that the two time horizon agree.
  *
  * If this method *is* called, which has to happen before that deserialize()
  * is called, and f_time_horizon is set at a value != 0, then if the 
  * dimension TimeHorizon is present in netCDF input, then the two values must
  * agree. If the dimension TimeHorizon is not present, then the value set by
  * this method is used. Note that, of course, the data in the netCDF file (if
  * the unit has any data indiced over the time horizon) has to agree with the
  * value set by this method.
  *
  * If this method is called *after* that deserialize() is called, this is
  * taken to mean that the UnitBlock is being "reset", and that immediatley
  * after deserialize() will be called again. The same rules as above are to
  * be followed for that subsequent call to deserialize(). */

 void set_time_horizon( Index t ) { f_time_horizon = t; }

/*@} -----------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/

 Index f_time_horizon;  ///< the time horizon of the HB

 Index f_number_heat_units;  ///< The number of units of the HB

 Index f_number_intervals;   ///< the number of intervals

 std::vector< int >  v_change_intervals;  ///< the vector of change interval

 // TODO: no, HeatDemand is indexed over TimeHorizon!!

 /// the vector of HeatDemand indexed over the dimensions NumberIntervals
 std::vector< double > v_heat_demand;

 // TODO: if it's a matrix, why not using a boost::multi_array<>?

 /** the matrix of MinHeatProduction indexed over the dimensions
  * NumberIntervals and NumberHeatUnits */
 std::vector< double > v_min_heat_production;

 /** the matrix of MaxHeatProduction indexed over the dimensions
  * NumberIntervals and NumberHeatUnits */
 std::vector< double > v_max_heat_production;

 /// the vector of MinHeatStorage indexed over the dimensions NumberIntervals
 std::vector< double > v_min_heat_storage;

 /// the vector of MaxHeatStorage indexed over the dimensions NumberIntervals
 std::vector< double > v_max_heat_storage;

 /** the matrix of CostHeatUnit indexed over the dimensions
  * NumberIntervals and NumberHeatUnits */
 std::vector< double > v_cost_heat_unit;

 /// Value of storing heat rho
 double f_storing_heat_rho;

 /// Value of extracting heat rho
 double f_extracting_heat_rho;

 /// Value of keeping heat rho
 double f_keeping_heat_rho;

 /// the initial amount of heat in the storage at the beginning of the time t
 double f_initial_heat_storage;

/*-----------------------------variables------------------------------------*/

 // TODO: that's a matrix, you should be using a boost::multi_array<>
 /// Vector of Heat variables
 std::vector< ColVariable > v_heat;

 /// Vector of HeatAdded variables
 std::vector< ColVariable > v_heat_added;

 /// Vector of HeatRemoved variables
 std::vector< ColVariable > v_heat_removed;

 /// Vector of HeatAvailable variables
 std::vector< ColVariable > v_heat_available;

/*----------------------------constraints-----------------------------------*/

 /// the heat demand satisfaction constraints
 std::vector< FRowConstraint > v_HeatDemand_Constraints;

 /// the heat bound satisfaction constraints
 boost::multi_array<FRowConstraint, 2> v_HeatBounds_Constraints;

 /// the heat storage bound satisfaction constraints
 std::vector< FRowConstraint > v_HeatStorageBounds_Constraints;

 /// the evolution in the  stored heat constraints
 std::vector< FRowConstraint > v_EvolutionStoredHeat_Constraints;

 /// the objective function
 FRealObjective objective;

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

// returns which variables must be generated
/* This method returns an int that indicates which variables of HeatBlock
 * must be generated by the generate_abstract_variables() method. This value
 * may be given in stvv as explained in generate_abstract_variables(). If this
 * value is not given in stvv, then this method returns the appropriate value
 * according to what is specified in the generate_abstract_variables() method.
 */

 int get_variables_to_be_generated( Configuration *stvv );

 void deserialize_time_horizon( netCDF::NcGroup & group );

 void deserialize_change_intervals( netCDF::NcGroup & group );

/*--------------------------------------------------------------------------*/

  }; // end( class( HeatBlock ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 }; // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/

#endif /* HeatBlock.h included */

/*--------------------------------------------------------------------------*/
/*-------------------------- End File HeatBlock.h --------------------------*/
/*--------------------------------------------------------------------------*/
