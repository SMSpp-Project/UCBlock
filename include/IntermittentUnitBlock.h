/*--------------------------------------------------------------------------*/
/*----------------------- File IntermittentUnitBlock.h ---------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class IntermittentUnitBlock, which derives from
 * UnitBlock [see UnitBlock.h], in order to define a Unit representing
 * Intermittent Generation in the Unit Commitment Problem.
 *
 * \version 0.11
 *
 * \date 02 - 08 - 2019
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
 *
 * Copyright &copy by Antonio Frangioni, Ali Ghezelsoflu
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __IntermittentUnitBlock
#define __IntermittentUnitBlock
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "ColVariable.h"
#include "FRowConstraint.h"
#include "OneVarConstraint.h"
#include "UnitBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------------ NAMESPACE ---------------------------------*/
/*--------------------------------------------------------------------------*/

/// Namespace for the Structured Modeling System++ (SMS++)

namespace SMSpp_di_unipi_it {

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS IntermittentUnitBlock -----------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// implementation of the Block concept for the Intermittent Generation unit
/** The IntermittentUnitBlock class implements the Block concept [see Block.h]
 * for a units representing generation (be ir centralized or distributed) by
 * intermittent (= unreliable) sources in the unit commitment problem, such
 * as wind farms, solar parks and run-of-the-river hydroelectricity. Each
 * unit is supposed to be connected to a specific node of the clustered
 * network (which means that the "distributed" case refers to "distributed in
 * a small region", where of course "small" depends on the granularity of the
 * network description. The model relies mainly on historical data of local
 * generation of wind and solar at each node of the grid; these data are
 * used to develop normalized generation profiles associated with wind and
 * solar generators. Intermittent generators are supposed to be able to
 * contribute to primary and secondary reserves. Contribution to the system
 * inertia concerns more specifically run of river generators. The potential
 * contribution of solar or wind generation to inertia is still the subject
 * of active research. Reserve requirements are specified in order to be
 * symmetrically available to increase or decrease power injected into the
 * grid. Then the technical and physical constraints are mainly divided in
 * three different categories:
 *
 * - the active power bounds;
 *
 * - the maximum power output constraints according to primary and secondary
 *   spinning reserves;
 *
 * - the minimum power output constraints according to primary and secondary
 *   spinning reserves. */

class IntermittentUnitBlock : public UnitBlock {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:
/*--------------------------------------------------------------------------*/
/*---------------------- PUBLIC TYPES OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name constructor and destructor
 *  @{ */

/// constructor, takes the father and the time horizon
/** Constructor of IntermittentUnitBlock, taking possibly a pointer of its
 * father Block. */

 explicit IntermittentUnitBlock( Block * f_block = nullptr , Index t = 0 )
  : UnitBlock( f_block ) {}

/*--------------------------------------------------------------------------*/
/// destructor of IntermittentUnitBlock

 ~IntermittentUnitBlock() override = default;

/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

/// extends Block::deserialize( netCDF::NcGroup )
/** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
 * the IntermittentUnitBlock. Besides the mandatory "type" attribute of any
 * :Block, the group must contain all the data required by the base UnitBlock,
 * as described in the comments to UnitBlock::deserialize( netCDF::NcGroup ).
 * In particular, we refer to that description for the crucial dimensions
 * "TimeHorizon", "NumberIntervals" and "ChangeIntervals". The netCDF::NcGroup
 * must then also contain:
 *
 * - The variable "MinPower", of type double and either of size 1 or indexed
 *   over the dimension "NumberIntervals". This is meant to represent the
 *   vector MinP[ t ] that, for each time instant t, contains the minimum
 *   potential production value of the unit for the corresponding time step.
 *   If "MinPower" has length 1 then MinP[ t ] contains the same value for all
 *   t. Otherwise, MinPower[ i ] is the fixed value of MinP[ t ] for all t in
 *   the interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with
 *   the assumption that ChangeIntervals[ - 1 ] = 0. If NumberIntervals <= 1
 *   or NumberIntervals >= TimeHorizon, then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded. Note that it must
 *   be MnP[ t ] >= 0 for all t.
 *
 * - The variable "MaxPower", of type double and either of size 1 or indexed
 *   over the dimension "NumberIntervals". This is meant to represent the
 *   vector MaxP[ t ] that, for each time instant t, contains the maximum
 *   potential production value of the unit for the corresponding time step.
 *   If "MaxPower" has length 1 then MaxP[ t ] contains the same value for all
 *   t. Otherwise, MaxPower[ i ] is the fixed value of MaxP[ t ] for all t in
 *   the interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with
 *   the assumption that ChangeIntervals[ - 1 ] = 0. If NumberIntervals <= 1
 *   or NumberIntervals >= TimeHorizon, then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded. Note that it must
 *   be MxP[ t ] >= MnP[ t ] [>= 0] for all t. Yet, MxP[ t ] == MnP[ t ] is
 *   possible: it means that (at time instant t) the unit cannot be curtailed
 *   and cannot provide any reserve.
 *
 * - The scalar variable "Gamma", of type double and not indexed over any
 *   dimension. This variable is used to take into account an uncertainty on
 *   the maximal potential production. Note that it must be 0 <= Gamma <= 1;
 *   when Gamma == 0, the unit does not provide any reserve.
 *
 * - The variable "InertiaPower", of type double and either indexed over the
 *   dimension "NumberIntervals" or has size 1. This is meant to represent the
 *   vector IP[ t ] which, for each time instant t, contains the contribution
 *   that the unit can give to the inertia constraint which depends on the
 *   active power that it is currently generating (basically, the constant to
 *   be multiplied to the active power variable) at time t for this unit. The
 *   variable is optional; if it is not defined, IP[ t ] == 0 for each time
 *   instants t. If it has size 1 then the entry IP[ 0 ] is assumed to contain
 *   the the inertia power value for this unit and all time instants t.
 *   Otherwise, InertiaPower[ i ] is the fixed value of IP[ t ] for all t in
 *   the interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with
 *   the assumption that ChangeIntervals[ - 1 ] = 0. If
 *   NumberIntervals <= 1 or NumberIntervals >= TimeHorizon then the mapping
 *   clearly does not require "ChangeIntervals", which in fact is not loaded.
 *
 *   */

 void deserialize( netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
/// generate the abstract variables of the IntermittentUnitBlock
/** The IntermittentUnitBlock class use get_variable() method to access to
 *  each "group" of variable that may create in UnitBlock class which are:
 *
 *  - the primary spinning reserve variables;
 *
 *  - the secondary spinning reserve variables;
 *
 *  - the active power variables.
 *
 *  All of those variables are optional except the active power variables in
 *  the sense that the model may just not have them and whenever a group of
 *  above variables is created, its size will be the time horizon. It is
 *  possible to restrict which of the subsets are generated with the parameter
 *  stvv. If stvv is not nullptr and it is a SimpleConfiguration<int>, or if
 *  f_BlockConfig->f_static_variables_Configuration is not nullptr and it is a
 *  SimpleConfiguration<int>, then the f_value (an int) indicates whether each
 *  of the optional variables should be created. If the Configuration is not
 *  available, the default value is taken to be 0.
 * */
 void generate_abstract_variables( Configuration *stvv ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// generate the static constraint of the IntermittentUnitBlock
/** Method that generates the static constraint of the IntermittentUnitBlock.
 * These are the:
 * - maximum and minimum power output constraints according to primary and
 *   secondary spinning reserves are presented in (1)-(2). Each of them is a
 *   std::vector<FRowConstraint>; with the dimension of f_time_horizon, where
 *   the entry t = 0, ...,f_time_horizon - 1 being the maximum and minimum
 *   power output value according to the primary and the secondary spinning
 *   reserves at time t. these ensure the maximum(or minimum) amount of energy
 *   that unit can produce(or use) when it is on(or off).
 *   \f[
 *       p^{pr}_{t} + p^{sc}_{t} \leq \gamma(P^{mx}_{t} - p^{ac}_{t} )
 *          \quad t \in \mathcal{T}                              \quad (1)
 *   \f]
 *
 *   \f[
 *       p^{pr}_{t} + p^{sc}_{t} \leq  p^{ac}_{t} - P^{mn}_{t}
 *          \quad t \in \mathcal{T}                              \quad (2)
 *   \f]
 *   where \f$ P^{mx}_{t} \f$ and \f$ P^{mn}_{t} \f$ are the maximum and
 *   minimum power output parameters for each time t of the time horizon
 *   \f$ \mathcal{T} \f$ respectively.
 *
 * - the active power bounds.
 *
 *   \f[
 *    p^{ac}_{t} \in [ P^{mn}_{t} , P^{mx}_{t}]
 *                              \quad t \in \mathcal{T}          \quad (3)
 *   \f]
 *   */
 void generate_abstract_constraints( Configuration *stcc ) override;

/**@} ----------------------------------------------------------------------*/
/*------- METHODS FOR READING THE DATA OF THE IntermittentUnitBlock --------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the IntermittentUnitBlock
 *
 * These methods allow to read data that must be common to (in principle) all
 * the kind of Intermittent Generation units
 * @{ */

 /// Returns the gamma value
 double get_gamma() const { return f_gamma; }
/*--------------------------------------------------------------------------*/
/// returns the vector of minimum power
/** The method returned a std::vector< double > V and each element of V
 * contains the minimum power at time t. There are three possible cases:
 *
 * - if the vector is empty, then the minimum power of the unit is 0;
 *
 * - if the vector has only one element, then V[ 0 ] is the minimum power of
 *   the unit for all time horizon;
 *
 * - otherwise, the std::vector< double > V must have size get_time_horizon()
 *   and each V[ t ] represents the minimum power value at time t. */

 const std::vector< double > & get_minimum_power() const {
  return( v_minimum_power );
 }
/*--------------------------------------------------------------------------*/
/// returns the vector of maximum power
/** The method returned a std::vector< double > V and each element of V
 * contains the maximum power at time t. There are three possible cases:
 *
 * - if the vector is empty, then the maximum power of the unit is 0;
 *
 * - if the vector has only one element, then V[ 0 ] is the maximum power of
 *   the unit for all time horizon;
 *
 * - otherwise, the std::vector< double > V must have size get_time_horizon()
 *   and each V[ t ] represents the maximum power value at time t. */

 const std::vector< double > & get_maximum_power() const {
  return( v_maximum_power );
 }

/*--------------------------------------------------------------------------*/
/// returns the matrix of inertia power
/** The returned value U = get_inertia_power() contains the contribution
 *  to inertia (basically, the constants to be multiplied by the active power
 *  variables returned by get_active_power()) of all the generators at all
 *  time instants. There are four possible cases:
 *
 * - if the matrix is empty, then the inertia power is always 0;
 *
 * - if the matrix only has one row (i.e., the first dimension has size 1),
 *   then the inertia power for each generator g is U[ 0 , g ] for all t
 *   which means that the second dimension has size get_number_generators();
 *
 * - if the matrix only has one column with size get_time_horizon() (i.e., the
 *   second dimension has size 1), then the InertiaPower[ t , 0 ] gives the
 *   inertia power for the problem at time t. Since in this unit there is only
 *   one electrical generator, this case should happen by assumption;
 *
 * - otherwise, the matrix has size get_time_horizon() per
 *   get_number_generators(), then the InertiaPower[ t , g ] represents
 *   the inertia power for the problem at time t for each electrical generator
 *   g. */

 const boost::multi_array< double , 2 > & get_inertia_power()
 const override {
  return( v_inertia_power );
 }
/**@} ----------------------------------------------------------------------*/
/*-------------- METHODS FOR SAVING THE IntermittentUnitBlock---------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the IntermittentUnitBlock
 *  @{ */

/// extends Block::serialize( netCDF::NcGroup )
/** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
 * IntermittentGenerationUnitBlock. See
 * IntermittentGenerationUnitBlock::deserialize( netCDF::NcGroup ) for details of the
 * format of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/**@} ----------------------------------------------------------------------*/
/*----------- METHODS FOR INITIALIZING THE IntermittentUnitBlock -----------*/
/*--------------------------------------------------------------------------*/

/** @name Handling the data of the IntermittentUnitBlock
    @{ */

 void load( std::istream & input ) override {
  throw ( std::logic_error( "IntermittentUnitBlock::load() not "
                            "implemented yet") );
 };

/**@} ----------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------data--------------------------------------*/
 /// the vector of MinPower
 std::vector< double >  v_minimum_power;

 /// the vector of MaxPower
 std::vector< double >  v_maximum_power;

 /// the gamma value
 double f_gamma;

 /// the matrix of inertia power of generators
 boost::multi_array< double , 2 > v_inertia_power;

/*----------------------------constraints-----------------------------------*/
/// the active power upper bound constraints
 std::vector< FRowConstraint > active_power_upper_bound_Constraints;

/// the active power lower bound constraints
 std::vector< FRowConstraint > active_power_lower_bound_Constraints;

/// the active power bounds constraints
 std::vector< FRowConstraint > active_power_bounds_Constraints;

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

};  // end( class( IntermittentUnitBlock ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* IntermittentUnitBlock.h included */

/*--------------------------------------------------------------------------*/
/*------------------ End File IntermittentUnitBlock.h ----------------------*/
/*--------------------------------------------------------------------------*/
