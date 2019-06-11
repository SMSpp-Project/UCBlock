/*--------------------------------------------------------------------------*/
/*---------------------------- File HeatBlock.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *derived* class HeatBlock, which derives from
 * the Block, in order to define a base for any possible unit that produces
 * heat and can be attached to a UCBlock. It has very basic information that
 * can characterize almost any different kind of heat-unit, which includes a
 * set of heat variables. This class has thus been constructed having the
 * following elements:
 *
 * - A virtual public method that is used to initialize and read the
 *   data of any possible derived HeatBlock class.
 *
 * - The time horizon of the problem.
 *
 * - a vector of ColVariable objects, that are used to store the
 *   information regarding:
 *
 *     (i)   the heat produced by the unit.
 *
 *   This vector either have size equal to the time horizon
 *   or is empty, in which case the corresponding variables simply do
 *   not exist (for instance, the unit may not produce heat).
 *
 * \version 0.11
 *
 * \date 11 - 06 - 2019
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
 * Copyright &copy by Antonio Frangioni, Ali Ghezelsoflu, and Rafael
 * Durbano Lobato
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __HeatBlock
#define __HeatBlock
/* self-identification: #endif at the end of the file */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Block.h"
#include "ColVariable.h"
#include <map>
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

/// implementation of the Block concept for the heat unit problem
/** The HeatBlock class implements the Block concept [see Block.h]
 * for the EDF Unit Commitment Problem.
 *
 *  The constraints regarding heat management can be described separately for
 *  each “heat block” (HB). In SMS++ parlance, a HB is the intersection of an
 *  Energy Cell and a heat-ID: a set of heat-producing units that are
 *  geographically near enough to exchange heat of the same type, together
 *  with demand and possibly storage of this heat. HBs are only connected to
 *  each other because some heat-producing units are also electricity
 *  producing ones, where heat is a by-product or a co-product of electricity.
 *  However, this linking “happens at the level of the main UC model”, and
 *  therefore it is not required for the description of the internal
 *  constraints of a HB.
 *
 *  The variable to the HB portion of the EUC model are the following:
 *
 * - \f$ p^{he}_{t,i} \f$ : representing the power produced by heat unit
 *   \f$ i \in \mathcal{I}(h) \f$ at time \f$ t \in \mathcal{T} \f$;
 *
 * - if the heat storage is defined, \f$ s_{t,+} \geq 0 \f$ and \f$ s_{t,-}
 *   \geq 0 \f$ representing respectively the amount of heat added to and
 *   removed from the storage at time instant \f$ t \in \mathcal{T} \f$ ;
 *
 * - if the heat storage is defined, \f$ v_t \f$ representing respectively the
 *   amount of heat available in the storage at time instant
 *   \f$ t \in \mathcal{T} \f$;
 *
 *  The heat constraints of unit commitment problem, on the time horizon
 *  \f$ \mathcal{T} \f$ write as follow:
 *
 * - Demand Constraints:
 *   In the unit commitment problem, \f$ D^h_{t} \f$ denotes the heat demand
 *   of the each heat block \f$ h \in \mathcal{H} \f$ for each time period
 *   \f$ t \in \mathcal{T} \f$ will be satisfied as follow:
 *
 * \f[
 *  \sum_{ i \in \mathcal{I}_h } (p^{h , he}_{t,i} - s^{h}_{t , +} +
 *  s^{h}_{t , -} \geq D^h_{t} \quad t \in \mathcal{T}           \quad     (1)
 * \f]
 *
 * - Heat production bounds Constraints:
 *   Let \f$ P^{h , mn}_{t , i} \f$ and \f$ P^{h , mx}_{t , i} \f$ denotes the
 *   minimum and maximum heat production respectively. For each heat unit
 *   \f$ i \in \mathcal{I}(h) \f$ in heat block \f$ h \in \mathcal{H}\f$ and
 *   at time \f$ t \in \mathcal{T} \f$, the heat production bounds will be
 *   satisfied as follow:
 *
 * \f[
 *   p^{he}_{t,i} \in [ P^{h , mn}_{t , i} ,  P^{h , mx}_{t , i}]
 *    \quad i \in \mathcal{I}(h) \quad t \in \mathcal{T}         \quad     (2)
 * \f]
 *
 * - Heat storage bounds Constraints:
 *   Let \f$ V^{h , mn}_{t} \f$ and \f$ V^{h , mx}_{t} \f$ denotes the
 *   minimum and maximum heat storage respectively. For each heat block
 *   \f$ h \in \mathcal{H} \f$  and at time \f$ t \in \mathcal{T} \f$, the heat
 *   storage bounds will be presented as follow:
 *
 * \f[
 *      v^h_{t} \in [v^{h , mn}_{t} , V^{h , mx}_{t}]
 *                          \quad t \in \mathcal{T}           \quad        (3)
 * \f]
 *
 * - Evolution in the stored heat Constraints:
 *   Let three constants \f$ \rho^{h}_{+} \f$, \f$ \rho^{h}_{-} \f$, and
 *   \f$ \rho^{h} \f$ representing inefficiencies in, respectively, storing
 *   heat in the heat storage, extracting heat from the heat storage, and
 *   keeping heat in the heat storage. Then the evolution in the stored heat
 *   constraints will be presented as follow:
 *
 * \f[
 *   v^{h}_{t} = \rho^{h} v^{h}_{t-1} + \rho^{h}_{+} s^{h}_{t , +} +
 *     \rho^{h}_{-} s^{h}_{t , -}  \quad t \in \mathcal{T}          \quad  (4)
 * \f]
 *
 * - Objective Function:
 *   We will now consider one specific fixed HB \f$ h \in \mathcal{H} \f$. The
 *   portion of the objective function of the EUC corresponding to this HB
 *   simply reads:
 *
 * \f[
 *   min \sum_{ i \in \mathcal{I}_h } \sum_{ t \in \mathcal{T} }
 *                       C^h_{t,i}p^{h , he}_{t,i}
 * \f]
 *   Where \f$  C^h_{t,i} \f$ is the cost of producing one heat unit
 *   \f$ i \in \mathcal{I}(h) \f$ at time \f$ t \in \mathcal{T} \f$.
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

    typedef unsigned int Index;

/*@}------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

/** Constructor of HeatBlock, taking possibly a pointer of its
 * father Block. */

    HeatBlock( Block * father_block = nullptr , Index t = 0 )
        : Block( father_block ), f_time_horizon( t ) {}

/*--------------------------------------------------------------------------*/

/// destructor of HeatBlock

    virtual ~HeatBlock() { };


/*@}------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

/// loads the HeatBlock instance from memory
/** Loads the HeatBlock instance from memory.
 * Like load( std::istream & ), if there is any Solver attached to this
 * HeatBlock then a NBModification (the "nuclear option") is issued.
 * */

    virtual void load( std::istream &input ) override { };

/*--------------------------------------------------------------------------*/
/// extends Block::deserialize( netCDF::NcGroup )
/** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
 * the HeatBlock. Besides the mandatory "type" attribute of any :Block,
 * the group should contain the following:
 *
 * - the dimension "NumberHeatUnits" containing the number heat-producing
 *   units;
 *
 * - the dimension "TimeHorizon" containing the time horizon;
 *
 * Note 1: consider set time horizon \f$\{0, \dots, "TimeHorizon-1"\}\f$ with
 * dimension "TimeHorizon". It can be presented as the union of some
 * intervals.
 *
 * Note 2: let's suppose the values of each variable may change independently,
 * and may have different values for some intervals along the set time horizon
 * \f$\{0, \dots, "TimeHorizon-1"\}\f$. It means each variable has its own
 * changes along some intervals independently. Without loss of generality,
 * let's take the union of the intersection of all the intervals between each
 * two variables separately. It means, at the end we may have a set of
 * intervals such as: \f$ [0 , a], [a+1 , b], \dots, [j+1 , k], [k+1 ,
 * TimeHorizon-1] \f$which where they cover all the changes for all the
 * variables. Then in each interval, one variable may change or not(if not,
 * copy the corresponding value of its' previous interval).
 *
 * - the dimension "NumberIntervals" which is a subset of \f$ \{1, ...,
 *   "TimeHorizon"\}\f$ and indicates the number of above intervals \f$([0
 *   , a] , [a+1 , b], ... , [j+1 , k] , [k+1, TimeHorizon-1])\f$ where the
 *   variables change. This dimension is optional. If it is not
 *   provided then it is taken to be 1.
 *
 *   Three scenarios may happen:
 *
 *    i). In the simplest case scenario, "NumberIntervals = 1" which
 *        means that the value of each variable does not change, i.e.,
 *        it is the same for each period in \f$ \{1, \dots ,
 *        "TimeHorizon"\}\f$.
 *
 *   ii). In the average case scenario, "1 < NumberIntervals <
 *        TimeHorizon", which means that in some time steps the values
 *        of some of the variables are changing.
 *
 *  iii). In the worst case scenario, "NumberIntervals = TimeHorizon",
 *        which means that in every time step, the value of each
 *        variable may change.
 *
 * - the variable "ChangeIntervals", of type integer and indexed over
 *   the dimension "NumberIntervals"; the \f$t_{th}\f$ entry of the
 *   variable indicates the positive number of \f$ a, b, ..., k,
 *   TimeHorizon-1\f$ on the above example.
 *
 * - the variable "TotalHeatDemand", of type double and indexed over the
 *   dimensions "NumberIntervals": entry HeatDemand[ t ] is assumed to contain
 *   the total heat demand of this heat block to be satisfied for the
 *   corresponding interval t;
 *
 * - the variable "CostHeatUnit", of type double and indexed both over the
 *   dimensions "NumberIntervals" and "NumberHeatUnits"; the entry
 *   CostHeatUnit[ t , i ] indicates the cost of producing value for each
 *   interval t of heat unit i in this heat block; if NumberHeatUnits == 0
 *   (say, it is not provided at all) then this variable need not be defined,
 *   since it is not loaded;
 *
 * - the variable "MinHeatProduction", of type double and indexed over the
 *   dimensions "NumberIntervals" and "NumberHeatUnits"; the entry
 *   MinHeatProduction[ t , i ] indicates the minimum heat production value
 *   for each  interval t of heat unit i in this heat block; if
 *   NumberHeatUnits == 0 (say, it is not provided at all) then this variable
 *   need not be defined, since it is not loaded;
 *
 * - the variable "MaxHeatProduction", of type double and indexed over the
 *   dimensions "NumberIntervals" and "NumberHeatUnits"; the entry
 *   MaxHeatProduction[ t , i ] indicates the maximum heat production value
 *   for each  interval t of heat unit i in this heat block; if
 *   NumberHeatUnits == 0 (say, it is not provided at all) then this variable
 *   need not be defined, since it is not loaded;
 *
 * - the variable "MinHeatStorage", of type double and indexed over the
 *   dimensions "NumberIntervals": entry MinHeatStorage[ t ] is assumed to
 *   contain the minimum heat storage of this heat block for each interval t;
 *
 * - the variable "MaxHeatStorage", of type double and indexed over the
 *   dimensions "NumberIntervals": entry MaxHeatStorage[ t ] is assumed to
 *   contain the maximum heat storage of this heat block for each interval t;
 *
 * - the scalar variable "InitialHeatAvailable", of type double and not
 *   indexed over any dimension and indicates the the initial amount of heat
 *   in the storage at the time -1 in this heat block; this variable is
 *   optional, if it is not provided it is taken to be
 *   InitialHeatAvailable == MinHeatStorage [ 0 ];
 *
 * - the scalar variable "StoringHeatRho", of type double and not indexed over
 *   any dimension and indicates the storing heat in the heat storage(if any)
 *   in this heat block; this variable is optional and always
 *   StoringHeatRho <= 1, if it is not provided it is taken to be
 *   StoringHeatRho == 0; if for all intervals t,
 *   MaxHeatStorage[ t ] == MinHeatStorage[ t ] then the heat block has no
 *   heat storage then this variable need not be defined, since they are not
 *   loaded;
 *
 * - the scalar variable "ExtractingHeatRho", of type double and not indexed
 *   over any dimension and indicates the extracting heat in the heat storage(
 *   if any) in this heat block; this variable is optional and always
 *   ExtractingHeatRho >= 1, if it is not provided it is taken to be
 *   ExtractingHeatRho == 0; if for all intervals t,
 *   MaxHeatStorage[ t ] == MinHeatStorage[ t ] then the heat block has no
 *   heat storage then this variable need not be defined, since they are not
 *   loaded;
 *
 * - the scalar variable "KeepingHeatRho", of type double and not indexed over
 *   any dimension and indicates the keeping heat in the heat storage(if any)
 *   in this heat block; this variable is optional and always
 *   KeepingHeatRho <= 1, if it is not provided it is taken to be
 *   KeepingHeatRho == 0; if for all intervals t,
 *   MaxHeatStorage[ t ] == MinHeatStorage[ t ] then the heat block has no
 *   heat storage then this variable need not be defined, since they are not
 *   loaded;
 */

    virtual void deserialize( netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/

    virtual void generate_abstract_variables( Configuration *stvv = nullptr )
    override;

/// generate the abstract variables of the HeatUnit
/** Method that generates the abstract variables of the HeatBlock.
 *
 * - the heat added variables [bit 0]
 *
 * - the heat removed variables [bit 1]
 *
 * - the heat available variables [bit 2]
 *
 * - the heat variables [bit 3]
 *
 * Note: all of these variables are optional, except the heat variables. The
 * parameter stvv is used to decide which of the optional variables should be
 * used. If stvv is not nullptr and it is a SimpleConfiguration<int> or if
 * f_BlockConfig->f_static_constraints_Configuration is not nullptr and it is
 * a SimpleConfiguration<int>, then the f_value (an int of this (the first
 * possible) Configuration indicates whether each of the optional variables
 * should be used. This is done according to the corresponding bit of f_value.
 * If the bit associated with a variable is 1, then the variable should be
 * used; otherwise, it should not. The first bit is associated with the
 * heat added variables, the second one with the heat removed variables and so
 * on according to the order the variables are listed above. Whenever a group
 * of variables should be used, its size will be the time horizon. In any
 * other case (i.e., if a SimpleConfiguration<int> is not present), none of
 * the variables above is considered.
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
 * */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// generate the static constraint of the HeatBlock
/** Method that generates the static constraint of the HeatBlock.
 * These are the:
 * //TODO
 */

    virtual void generate_abstract_constraints( Configuration *stcc = nullptr )
    override ;



/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// generate the objective of the HeatBlock
/** Method that generates the objective of the HeatBlock.
 *
*/
    virtual void generate_objective( Configuration *objc = nullptr )
    override ;
/*@} -----------------------------------------------------------------------*/
/*-------------- Methods for reading the data of the HeatBlock -------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for reading the data of the HeatBlock
 *  @{  */

/// returns the time horizon of the problem
    Index get_time_horizon( void ) const { return f_time_horizon; }

/// returns the number of units in this HeatBlock
    Index get_number_heat_units( void ) const { return f_number_heat_units; }

/** returns the minimum heat production of the given block for interval t of
 * unit i */
    inline double get_min_heat_production(  Index interval , Index unit) const {
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

/// Method for returning the vector of heat variables
    const std::vector<ColVariable> & get_heat( void ) const {
        return v_heat;
    }
/// Method for returning the pointer to the heat variable at time t
    ColVariable * get_heat( int t ) { return & ( v_heat[t] ); }


/*@} -----------------------------------------------------------------------*/
/*--------------------- METHODS FOR SAVING THE HeatBlock -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the HeatBlock
 *  @{ */

/// extends Block::serialize( netCDF::NcGroup )
/** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
 * HeatBlock. See HeatBlock::deserialize( netCDF::NcGroup ) for
 * details of the format of the created netCDF group.
 *
 * */

    virtual void serialize( netCDF::NcGroup & group ) const override;


/*@} -----------------------------------------------------------------------*/
/*----------------- METHODS FOR MODIFYING THE HeatBlock --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the HeatBlock
 *  @{ */

/// Set the time horizon
    void set_time_horizon( Index t );

/*@} -----------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

  protected:

/// the time horizon of the problem
    Index f_time_horizon;

/// The number of units of the problem
    Index f_number_heat_units;

/// the number of intervals
    Index f_number_intervals;

/// the vector of change interval
    std::vector< int >  v_change_intervals;

/// the vector of HeatDemand indexed over the dimensions NumberIntervals
    std::vector< double > v_heat_demand;

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

    SMSpp_insert_in_factory_h;
/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

/// returns which variables must be generated
/** This method returns an int that indicates which variables of HeatBlock
 * must be generated by the generate_abstract_variables() method. This value
 * may be given in stvv as explained in generate_abstract_variables(). If this
 * value is not given in stvv, then this method returns the appropriate value
 * according to what is specified in the generate_abstract_variables() method.
 * */
    int get_variables_to_be_generated( Configuration *stvv );

    void deserialize_time_horizon( netCDF::NcGroup & group );

    void deserialize_change_intervals( netCDF::NcGroup & group );


  }; // end( class( HeatBlock ) )

/*@}  end( class( HeatBlock ) ) --------------------------------------------*/
/*--------------------------------------------------------------------------*/
}; // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* HeatBlock.h included */

/*--------------------------------------------------------------------------*/
/*-------------------------- End File HeatBlock.h --------------------------*/
/*--------------------------------------------------------------------------*/
