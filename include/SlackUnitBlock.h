/*--------------------------------------------------------------------------*/
/*------------------------- File SlackUnitBlock.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class SlackUnitBlock, which derives from UnitBlock
 * [see UnitBlock.h].
 *
 * \version 0.11
 *
 * \date 10 - 12 - 2019
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
 * Copyright &copy by Antonio Frangioni, and Ali Ghezelsoflu
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __SlackUnitBlock
#define __SlackUnitBlock
/* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "FRowConstraint.h"
#include "OneVarConstraint.h"
#include "FRealObjective.h"
#include "DQuadFunction.h"
#include "UnitBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------------ NAMESPACE ---------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)

namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS SlackUnitBlock --------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ GENERAL NOTES -----------------------------*/
/*--------------------------------------------------------------------------*/
/// implementation of the Block concept for a slack unit
/** The SlackUnitBlock class derives from UnitBlock and implements a
 */

class SlackUnitBlock : public UnitBlock {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

 /// constructor, takes the father and the time horizon
 /** Constructor of SlackUnitBlock, taking possibly a pointer of its
  * father Block and the time horizon. */

 explicit SlackUnitBlock( Block * f_block = nullptr, Index t = 0 ) :
         UnitBlock( f_block ) {

  v_MaxInertia.resize( boost::extents[ 0 ][ 0 ]);
 }

/*--------------------------------------------------------------------------*/
 /// destructor of SlackUnitBlock, it is empty

 ~SlackUnitBlock() override = default;

/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

/// Extends Block::deserialize( netCDF::NcGroup )
/** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
 * the SlackUnitBlock. Besides the mandatory "type" attribute of any :Block,
 * the group must contain all the data required by the base UnitBlock, as
 * described in the comments to UnitBlock::deserialize( netCDF::NcGroup ).
 * In particular, we refer to that description for the crucial dimensions
 * "TimeHorizon", "NumberIntervals" and "ChangeIntervals". The netCDF::NcGroup
 * must then also contain:
 *
 * - The variable "MaxPower", of type double and either of size 1 or indexed
 *   over the dimension "NumberIntervals". This is meant to represent the
 *   vector MxP[ t ] that, for each time instant t, contains the maximum
 *   active power output value of the unit for the corresponding time step.
 *   If "MaxPower" has length 1 then MxP[ t ] contains the same value for all
 *   t. Otherwise, MaxPower[ i ] is the fixed value of MxP[ t ] for all t in
 *   the interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with
 *   the assumption that ChangeIntervals[ - 1 ] = 0. This variable is optional,
 *   if is not provided then MxP[ t ] == 0 for all t. Note that it must be
 *   MxP[ t ] >= 0 for all t. If NumberIntervals <= 1 or
 *   NumberIntervals >= TimeHorizon, then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded.
 *
 * - The variable "MaxPrimaryPower", of type double and either of size 1 or
 *   indexed over the dimension "NumberIntervals". This is meant to represent
 *   the vector MaxPP[ t ] that, for each time instant t, contains the maximum
 *   active power that can be used as primary reserve of the unit for the
 *   corresponding time step. If "MaxPrimaryPower" has length 1 then
 *   MaxPP[ t ] contains the same value for all t. Otherwise,
 *   MaxPrimaryPower[ i ] is the fixed value of MaxPP[ t ] for all t in the
 *   interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the
 *   assumption that ChangeIntervals[ - 1 ] = 0. This variable is optional,
 *   if is not provided then MaxPP[ t ] == 0 for all t. If NumberIntervals
 *   <= 1 or NumberIntervals >= TimeHorizon, then the mapping clearly does
 *   not require "ChangeIntervals", which in fact is not loaded.
 *
 * - The variable "MaxSecondaryPower", of type double and either of size 1 or
 *   indexed over the dimension "NumberIntervals". This is meant to represent
 *   the vector MaxSP[ t ] that, for each time instant t, contains the maximum
 *   active power that can be used as secondary reserve of the unit for the
 *   corresponding time step. If "MaxSecondaryPower" has length 1 then
 *   MaxSP[ t ] contains the same value for all t. Otherwise,
 *   MaxSecondaryPower[ i ] is the fixed value of MaxSP[ t ] for all t in the
 *   interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the
 *   assumption that ChangeIntervals[ - 1 ] = 0. This variable is optional,
 *   if is not provided then MaxSP[ t ] == 0 for all t. Note that
 *   MaxPP[ t ] == 0 implies MaxSP[ t ] == 0 (that is, if MaxPrimaryPower is
 *   not defined then neither should MaxSecondaryPower). If NumberIntervals
 *   <= 1 or NumberIntervals >= TimeHorizon, then the mapping clearly does
 *   not require "ChangeIntervals", which in fact is not loaded.
 *
 * - The variable "MaxInertia", of type double and either indexed over the
 *   dimension "NumberIntervals" or has size 1. This is meant to represent the
 *   vector MaxI[ t ] which, for each time instant t, contains the
 *   contribution that the unit can give to the inertia constraint for the
 *   sole fact that is is on (basically, the constant to be multiplied to the
 *   commitment variable) at time t. The variable is optional; if it is not
 *   defined, MaxI[ t ] == 0 for all time instants. If it has size 1, then
 *   MaxI[ t ] == MaxInertia[ 0 ] for all t, regardless to what
 *   "NumberIntervals" says. Otherwise, MaxInertia[ i ] is the fixed value of
 *   MaxI[ t ] for all t in the interval [ ChangeIntervals[ i - 1 ] ,
 *   ChangeIntervals[ i ] ], with the assumption that ChangeIntervals[ - 1 ]
 *   = 0.
 *
 * - The variable "ActivePowerCost", of type double and either indexed over
 *   the dimension "NumberIntervals" or has size 1. This is meant to represent
 *   the vector APC[ t ] that, for each time instant t, contains the cost of
 *   producing active power of the unit at the corresponding time step. This
 *   variable is optional; if it is not provided then it's taken to be zero.
 *   If "ActivePowerCost" has length 1 then APC[ t ] contains the same value
 *   for t. Otherwise, ActivePowerCost[ i ] is the fixed value of APC[ t ] for
 *   all t in the interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ]
 *   with the assumption that ChangeIntervals[ - 1 ] = 0. If NumberIntervals
 *   <= 1 or NumberIntervals >= TimeHorizon, then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded.
 *
 * - The variable "PrimaryCost", of type double and either indexed over the
 *   dimension "NumberIntervals" or has size 1. This is meant to represent the
 *   vector PC[ t ] that, for each time instant t, contains the cost of
 *   producing power that can be used as primary reserve of the unit at the
 *   corresponding time step. This variable is optional; if it is not provided
 *   then it's taken to be zero. If "PrimaryCost" has length 1 then PC[ t ]
 *   contains the same value for t. Otherwise, PrimaryCost[ i ] is the fixed
 *   value of PC[ t ] for all t in the interval [ ChangeIntervals[ i - 1 ] ,
 *   ChangeIntervals[ i ] ], with the assumption that ChangeIntervals[ - 1 ]
 *   = 0. If NumberIntervals <= 1 or NumberIntervals >= TimeHorizon, then the
 *   mapping clearly does not require "ChangeIntervals", which in fact is not
 *   loaded.
 *
 * - The variable "SecondaryCost", of type double and either indexed over the
 *   dimension "NumberIntervals" or has size 1. This is meant to represent the
 *   vector SC[ t ] that, for each time instant t, contains the cost of
 *   producing power that can be used as secondary reserve of the unit at the
 *   corresponding time step. This variable is optional; if it is not provided
 *   then it's taken to be zero. If "SecondaryCost" has length 1 then SC[ t ]
 *   contains the same value for t. Otherwise, SecondaryCost[ i ] is the fixed
 *   value of SC[ t ] for all t in the interval [ ChangeIntervals[ i - 1 ] ,
 *   ChangeIntervals[ i ] ], with the assumption that ChangeIntervals[ - 1 ]
 *   = 0. If NumberIntervals <= 1 or NumberIntervals >= TimeHorizon, then the
 *   mapping clearly does not require "ChangeIntervals", which in fact is not
 *   loaded.
 *   */

 void deserialize( netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
/// generate the abstract variables of the SlackUnitBlock
/** The SlackUnitBlock class use get_variable() method to access to each
 *  "group" of variable that may create in UnitBlock class which are:
 *
 *  - the binary commitment variables which takes the continues values between
 *    1 and zero.
 *
 *  - the primary spinning reserve variables;
 *
 *  - the secondary spinning reserve variables;
 *
 *  - the active power variables.
 */

 void generate_abstract_variables( Configuration *stvv ) override;

/*--------------------------------------------------------------------------*/
/// Generate the static constraint of the SlackUnitBlock
/** This method generates the abstract constraints of the SlackUnitBlock.
 *
 * The operations of the slack generating unit are described on a discrete
 * time horizon as dictated by the UnitBlock interface. In this description
 * we indicate it with \f$ \mathcal{T}=\{ 0, \dots , \mathcal{|T|} - 1\} \f$.
 * This unit just contains the bounds constraint on the ActivePower, Primary
 * and Secondary spinning reserve variables as below:
 *
 *   \f[
 *      0 \leq p^{ac}_{t} \leq P^{mx}_{t} \quad t \in \mathcal{T}    \quad (1)
 *   \f]
 *
 *   \f[
 *      0 \leq p^{pr}_{t} \leq P^{mxP}_{t} \quad t \in \mathcal{T}   \quad (2)
 *   \f]
 *
 *   \f[
 *      0 \leq p^{sc}_{t} \leq P^{mxS}_{t} \quad t \in \mathcal{T}   \quad (3)
 *   \f]
*/
 void generate_abstract_constraints( Configuration *stcc ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// generate the objective of the SlackUnitBlock
/** Method that generates the objective of the SlackUnitBlock.
 * - Objective function: the objective function of the SlackUnitBlock
 *   representing the total power production cost to be minimized has the
 *   form:
 *
 *   \f[
 *     \min ( \sum_{ t \in  [t_0 , \mathcal{T}]  } P^{mx}_{t} p^{ac}_{t} +
 *     P^{mxP}_{t} p^{pr}_{t} +P^{mxS}_{t} p^{sc}_{t} + P^{MaxI}_t v_t)
 *   \f]
 *
 *   where \f$ P^{mx}_{t} \f$, \f$ P^{mxP}_{t} \f$ , \f$ P^{mxS}_{t} \f$, and
 *   P^{MaxI} is the MaxPower, MaxPrimaryPower, and MaxSecondaryPower
 *   respectively.
 */

 void generate_objective( Configuration *objc ) override;

/**@} ----------------------------------------------------------------------*/
/*----------- METHODS FOR READING THE DATA OF THE SlackUnitBlock -----------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the SlackUnitBlock
 *
 * @{ */
/// returns the vector of maximum power
/** The returned vector contains to maximum power at time t. There are three
 * possible cases:
 *
 * - if the vector is empty, then the maximum power of the unit is 0;
 *
 * - if the vector has only one element, then the maximum power of the unit
 *   for all time horizon;
 *
 * - otherwise, the vector must have size get_time_horizon() and each element
 *   of vector represents the maximum power value at time t.   */

 const std::vector< double > & get_max_power() const {
  return( v_MaxPower );
 }
/*--------------------------------------------------------------------------*/
/// returns the vector of maximum primary power
/** The returned vector contains to maximum primary power at time t. There are
 * three possible cases:
 *
 * - if the vector is empty, then the maximum primary power of the unit is 0;
 *
 * - if the vector has only one element, then the maximum primary power of
 *   the unit for all time horizon;
 *
 * - otherwise, the vector must have size get_time_horizon() and each element
 *   of vector represents the maximum primary power value at time t.   */

 const std::vector< double > & get_max_primary_power() const {
  return( v_MaxPrimaryPower );
 }
/*--------------------------------------------------------------------------*/
/// returns the vector of maximum secondary power
/** The returned vector contains to maximum secondary power at time t. There
 * are three possible cases:
 *
 * - if the vector is empty, then the maximum secondary power of the unit is
 *   0;
 *
 * - if the vector has only one element, then the maximum secondary power of
 *   the unit for all time horizon;
 *
 * - otherwise, the vector must have size get_time_horizon() and each element
 *   of vector represents the maximum secondary power value at time t.   */

 const std::vector< double > & get_max_secondary_power() const {
  return( v_MaxSecondaryPower );
 }
/*--------------------------------------------------------------------------*/
/// returns the vector of maximum inertia
/** The returned value U = get_inertia_commitment() contains the contribution
 *  to inertia (basically, the constants to be multiplied by the commitment
 *  variables returned by get_commitment()) of all the generators at all time
 *  instants. There are four possible cases:
 *
 * - if the matrix is empty, then the maximum inertia is always 0;
 *
 * - if the matrix only has one row (i.e., the first dimension has size 1),
 *   then the maximum inertia for each generator g is U[ 0 , g ] for all t
 *   which means that the second dimension has size get_number_generators();
 *
 * - if the matrix only has one column with size get_time_horizon() (i.e., the
 *   second dimension has size 1), then the MaxInertia[ t , 0 ] gives the
 *   maximum inertia for the problem at time t. Since in this unit  there is
 *   only one electrical generator, this case should happen by  assumption;
 *
 * - otherwise, the matrix has size get_time_horizon() per
 *   get_number_generators(), then the MaxInertia[ t , g ] represents
 *   the maximum inertia for the problem at time t for each electrical
 *   generator g. */

 const boost::multi_array< double , 2 > & get_inertia_commitment()
 const override {
  return( v_MaxInertia );
 }
/**@} ----------------------------------------------------------------------*/
/*------------------ METHODS FOR SAVING THE SlackUnitBlock ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the SlackUnitBlock
 *  @{ */

/// Extends Block::serialize( netCDF::NcGroup )
/** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
 * SlackUnitBlock. See SlackUnitBlock::deserialize( netCDF::NcGroup ) for
 * details of the format of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/**@} ----------------------------------------------------------------------*/
/*--------------- METHODS FOR INITIALIZING THE SlackUnitBlock --------------*/
/*--------------------------------------------------------------------------*/

/** @name Handling the data of the SlackUnitBlock
    @{ */

 void load( std::istream & input ) override {
  throw ( std::logic_error( "SlackUnitBlock::load() not implemented yet") );
 }

/**@} ----------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------data--------------------------------------*/
 /// the vector of MaxPower
 std::vector< double > v_MaxPower;

 /// the vector of MaxPrimaryPower
 std::vector< double > v_MaxPrimaryPower;

 /// the vector of MaxSecondaryPower
 std::vector< double > v_MaxSecondaryPower;

 /// the matrix of MaxInertia
 boost::multi_array< double, 2 > v_MaxInertia;
/*----------------------------constraints-----------------------------------*/
 /// the active power bounds constraints
 std::vector< FRowConstraint > active_power_bounds_Constraints;

 /// the maximum primary bounds constraints
 std::vector< FRowConstraint > max_primary_bounds_Constraints;

 /// the maximum secondary bounds constraints
 std::vector< FRowConstraint > max_secondary_bounds_Constraints;
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

};  // end( class( SlackUnitBlock ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* SlackUnitBlock.h included */

/*--------------------------------------------------------------------------*/
/*------------------------ End File SlackUnitBlock.h -----------------------*/
/*--------------------------------------------------------------------------*/
