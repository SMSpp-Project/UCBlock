/*--------------------------------------------------------------------------*/
/*------------------------- File HydroUnitBlock.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class HydroUnitBlock, which derives from UnitBlock
 * [see UnitBlock.h], in order to define a "reasonably standard" hydro unit
 * of a Unit Commitment Problem.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Ali Ghezelsoflu \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Rafael Durbano Lobato \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Ali Ghezelsoflu,
 *                      Rafael Durbano Lobato
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __HydroUnitBlock
 #define __HydroUnitBlock
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "ColVariable.h"

#include "FRowConstraint.h"

#include "DQuadFunction.h"

#include "UnitBlock.h"

#include "OneVarConstraint.h"

#include "FRealObjective.h"

/*--------------------------------------------------------------------------*/
/*------------------------------ NAMESPACE ---------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)

namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS HydroUnitBlock ----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// implementation of the Block concept for the hydro unit problem
/** The HydroUnitBlock class implements the Block concept [see Block.h] for a
 * "reasonably standard" hydro unit of a Unit Commitment Problem. That is, the
 * class is designed in order to give mathematical formulation to describe the
 * operation of large set of hydro storage. To model complex reservoir systems
 * several technical parameters have to be considered. These are divided into
 * reservoir-specific parameters, the hydro links connecting the reservoirs
 * and finally the turbine/pump parameters. The values are collected within a
 * reservoir database, a hydro-link database and a turbine/pump-database. The
 * technical and physical constraints are mainly divided in several different
 * categories as:
 *
 * - maximum and minimum power output constraints according to primary and
 *   secondary spinning reserves;
 *
 * - primary and secondary spinning reserves relation with active power for
 *   turbines;
 *
 * - primary and secondary spinning reserves value for pumps ( == 0 );
 *
 * - flow-to-active-power function;
 *
 * - ramp-up and ramp-down constraints;
 *
 * - flow rate variable bounds;
 *
 * - final volumes of each reservoir constraints;
 *
 * - final volumes variable bounds. */

class HydroUnitBlock : public UnitBlock
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 * @{ */

/// constructor, takes the father block
/** Constructor of HydroUnitBlock, taking possibly a pointer of its father
 * Block. */

 explicit HydroUnitBlock( Block * f_block = nullptr )
  : UnitBlock( f_block ) {}

/*--------------------------------------------------------------------------*/
 /// destructor of HydroUnitBlock

 virtual ~HydroUnitBlock() override;

/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 * @{ */

/// extends Block::deserialize( netCDF::NcGroup )
/** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
 * the HydroUnitBlock. Besides the mandatory "type" attribute of any :Block,
 * the group must contain all the data required by the base UnitBlock, as
 * described in the comments to UnitBlock::deserialize( netCDF::NcGroup ).
 * In particular, we refer to that description for the crucial dimensions
 * "TimeHorizon", "NumberIntervals" and "ChangeIntervals". The netCDF::NcGroup
 * must then also contain:
 *
 * - The dimension "NumberReservoirs" containing the number of all reservoirs
 *   (or nodes) in the HydroUnitBlock. The dimension is optional, if it is not
 *   provided then it is taken to be == 1, which means that the (in principle)
 *   cascading hydro system is actually single hydro reservoir. Note, however,
 *   that a single reservoir can still have multiple hydro generating units
 *   (see NumberArcs below).
 *
 * - The dimension "NumberArcs" containing the set of arcs (or units)
 *   connecting the reservoirs in cascading system. Each arc represents either
 *   a turbine generating electricity by converting potential energy of water
 *   going downhill, or a pump consuming electricity for moving water uphill.
 *
 * - The variable "StartArc", of type netCDF::NcUint and indexed over the
 *   dimension "NumberArcs"; the r-th entry of the variable is the starting
 *   point of the arc (a number in 0, ..., NumberReservoirs - 1). Note that
 *   arcs are oriented; that is, a positive flow along arc r (turbine) means
 *   that water is being taken away from StartArc[ r ] and delivered to
 *   EndArc[ r ] (see next), a negative flow (pump) means vice-versa. Note
 *   that reservoir names here go from 0 to NumberReservoirs - 1.
 *
 * - The variable "EndArc", of type netCDF::NcUint and indexed over the
 *   dimension "NumberArcs"; the r-th entry of the variable is the ending
 *   point of the arc; this is a number in 0, ..., NumberReservoirs. Note:
 *   this is NumberReservoirs and *not* NumberReservoirs - 1, because arcs can
 *   end in the "fake" reservoir NumberReservoirs. This indicates that water
 *   that flows along that arc "goes away from the system" and it is no longer
 *   counted, because it can no longer be used to produce electricity further
 *   down the river, or pumped back into one of its reservoirs. Indeed, there
 *   will be something like "the most downstream turbines": after water has
 *   been used there, it just goes away down some river and does not go to any
 *   other reservoir. Arcs are oriented (see above); StartArc[ r ] == EndArc[
 *   r ] (a self-loop) is not allowed, but multiple arcs between the same pair
 *   of reservoirs are. Indeed, often the same physical equipment can be used
 *   both as a turbine and as a pump; in our model these are represented as
 *   two parallel arcs (but with different upper and lower flow capacity, see
 *   "MinFlow" and "MaxFlow" below).
 *
 * - The variable "MinFlow", of type netCDF::NcDouble and indexed over both
 *   dimensions "NumberIntervals" and "NumberArcs". The first dimension may
 *   have either size 1 or size "NumberIntervals" (if "NumberIntervals" is not
 *   provided, then the size can also be "TimeHorizon") whereas the second one
 *   always has size NumberArcs (if the variable is provided at all). This is
 *   meant to represent the matrix MinF[ t , l ] which, for each time instant
 *   t and each arc l, contains the minimum flow value of the unit. This
 *   variable is optional; if it is not provided then it is assumed that MinF[
 *   t , l ] == 0, i.e., the minimum flow of the unit is zero. If the first
 *   dimension has size 1 then the entry MinF[ 0 , l ] gives the fixed minimum
 *   flow value of the unit for all time steps and each arc l. Otherwise,
 *   MinFlow[ i , l ] is the fixed value of MinF[ t , l ] for all time t and
 *   arc l in the interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ]
 *   ], with the assumption that ChangeIntervals[ - 1 ] = 0. If
 *   NumberIntervals <= 1 or NumberIntervals >= TimeHorizon, then the mapping
 *   clearly does not require "ChangeIntervals", which in fact is not loaded.
 *
 * - The variable "MaxFlow", of type netCDF::NcDouble and indexed over both
 *   dimensions "NumberIntervals" and "NumberArcs". The first dimension may
 *   have either size 1 or size "NumberIntervals" (if "NumberIntervals" is not
 *   provided, then the size can also be "TimeHorizon") whereas the second one
 *   always has size NumberArcs (if the variable is provided at all). This is
 *   meant to represent the matrix MaxF[ t , l ] which, for each time instant
 *   t and each arc l, contains the maximum flow value of the unit. This
 *   variable is optional; if it is not provided then it is assumed that MaxF[
 *   t , l ] == 0, i.e., the maximum flow of the unit is zero. If the first
 *   dimension has size 1 then the entry MaxF[ 0 , l ] gives the fixed maximum
 *   flow value of the unit for all time steps and each arc l. Otherwise,
 *   MaxFlow[ i , l ] is the fixed value of MaxF[ t , l ] for all time t and
 *   arc l in the interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ]
 *   ], with the assumption that ChangeIntervals[ - 1 ] = 0. If
 *   NumberIntervals <= 1 or NumberIntervals >= TimeHorizon, then the mapping
 *   clearly does not require "ChangeIntervals", which in fact is not loaded.
 *
 * Note: MinFlow and MaxFlow values can be either positive or negative (or
 * zero); whenever MinF[ t , l ] <= MaxF[ t , l ] <= 0 for each t and l,
 * the unit is considered a pump and whenever 0 <= MinF[ t , l ] <=
 * MaxF[ t , l ], the unit is considered a turbine. Note that an arc must
 * *always* be the same kind for *all* instants, i.e., it is not allowed
 * that a unit suddenly changes between a turbine and a pump, or vice-versa.
 * This is because the flow-to-active-power function of turbines is a convex
 * piecewise function with possibly many pieces (see "NumberPieces",
 * "LinearTerm", "ConstantTerm" below) , whereas the flow-to-active-power
 * function of a pump is a simple linear function. In other words, the
 * "number of pieces" of a turbine is >= 1, whereas the "number of pieces"
 * of a pump is necessarily equal to 1. In reality, the same equipment can
 * sometimes be used both as a pump and as a turbine. In our model this is be
 * accounted for by artificially splitting the unit into “two units”, a pump
 * one and a turbine one, which must be done at the data processing stage.
 * This causes the possible problem that at some time instant both the pump
 * and the turbine be active, which is not possible in practice. This is
 * unlikely to happen (because pumps consume more than turbines produce for
 * the same amount of water, so this would be uneconomical), but should it
 * ever happen, this occurrence is not handled in our model (which lets it
 * happen).
 *
 * - The variable "MinVolumetric", of type netCDF::NcDouble and indexed over
 *   both dimensions "NumberReservoirs" and "NumberIntervals". The first
 *   dimension always has size "NumberReservoirs" (if it is provided at all),
 *   whereas the second one may have size one or size "NumberIntervals" (if
 *   "NumberIntervals" is not provided, then the size can also be
 *   "TimeHorizon"). This is meant to represent the matrix MinV[ r , t ]
 *   which, for each reservoir r at each time instant t contains the minimum
 *   volumetric value of the unit for each reservoir and corresponding time
 *   step. It must be that 0 <= MinV[ r , t ] <= MaxV[ r , t ] for all r and
 *   t. This variable is optional; if it's not provided then it is assumed
 *   that MinV[ r , t ] == 0, i.e., the minimum volumetric of the unit is
 *   zero. If the second dimension has size 1 then the entry MinV[ r , 0 ]
 *   gives the fixed minimum volumetric value of the unit for each reservoir r
 *   along all the time horizon. Otherwise, MinVolumetric[ r , i ] is the
 *   fixed value of MinV[ r , t ] for reservoir r and all t in the interval [
 *   ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the assumption
 *   that ChangeIntervals[ - 1 ] = 0. If NumberIntervals <= 1 or
 *   NumberIntervals >= TimeHorizon, then the mapping clearly does not require
 *   "ChangeIntervals", which in fact is not loaded.
 *
 * - The variable "MaxVolumetric", of type netCDF::NcDouble and indexed over
 *   both dimensions "NumberReservoirs" and "NumberIntervals". The first
 *   dimension always has size "NumberReservoirs" (if it is provided at
 *   all), whereas the second one may have size one or size "NumberIntervals"
 *   (if "NumberIntervals" is not provided, then the size can also be
 *   "TimeHorizon"). This is meant to represent the matrix MaxV[ r , t ]
 *   which, for each reservoir r at each time instant t contains the maximum
 *   volumetric value of the unit for each reservoir and corresponding time
 *   step. It must be that 0 <= MinV[ r , t ] <= MaxV[ r , t ] for all r and
 *   t. This variable is not optional (a reservoir must have some available
 *   volume). If the second dimension has size 1 then the entry MaxV[ r , 0 ]
 *   gives the fixed maximum volumetric value of the unit for each reservoir r
 *   during the all the time horizon. Otherwise, the MaxVolumetric[ r , i ] is
 *   the fixed value of MaxV[ r , t ] for reservoir r and all t in the
 *   interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the
 *   assumption that ChangeIntervals[ - 1 ] = 0. If NumberIntervals <= 1 or
 *   NumberIntervals >= TimeHorizon then the mapping clearly does not require
 *   "ChangeIntervals", which in fact is not loaded.
 *
 * Note: it may happen that MinV[ r , t ] == MaxV[ r , t ], but only *in a
 * subset of the time instants*. For instance, the user may want to fix the
 * final value of the reservoir, for whatever reason. So, if MinV and MaxV are
 * independent of t, then MinV[ r ] < MaxV[ r ] must surely happen. If,
 * instead, they depend on t, then equality can be accepted at some instants
 * (but not all of them).
 *
 * - The variable "Inflows", of type netCDF::NcDouble and indexed over both
 *   dimensions "NumberReservoirs" and "TimeHorizon". This is meant to
 *   represent the matrix InF[ r , t ] which, for each reservoir r, contains
 *   the amount of water that "naturally" enters into reservoir r (because of
 *   rain, ice melting, non-controlled rivers flowing, and of course net of
 *   water leaving by evaporation, human consumption etc.) during time all the
 *   time interval t, and therefore that is available in the reservoir at the
 *   end of time step t (hence, the beginning of time step t + 1, if
 *   any). This variable is optional; if it isn't defined, it is taken to be
 *   zero. Inflows can be either positive or negative.
 *
 * - The variable "MinPower", of type netCDF::NcDouble and indexed over both
 *   dimensions "NumberIntervals" and "NumberArcs". The first dimension may
 *   have either size 1 or size "NumberIntervals" (if "NumberIntervals" is not
 *   provided, then the size can also be "TimeHorizon"), whereas the second
 *   one always has size NumberArcs (if it is provided at all). This is meant
 *   to represent the matrix MinP[ t , l ] which, for each time instant t at
 *   each arc l contains the minimum power value of the unit; it must be that
 *   MinP[ t , l ] <= MaxP[ t , l ] for each time instant t and each arc l.
 *   This variable is optional; if it is not provided then it is assumed that
 *   MinP[ t , l ] == 0, i.e., the minimum power of all units is zero (which
 *   means, each unit is a turbine). If the first dimension has size 1 then
 *   the entry MinP[ 0 , l ] is assumed to contain the minimum power of arc l
 *   for all time instants. Otherwise, MinPower[ i , l ] is the fixed value of
 *   MinP[ t , l ] for arc l and all time t in the interval [ ChangeIntervals[
 *   i - 1 ] , ChangeIntervals[ i ] ], with the assumption that
 *   ChangeIntervals[ - 1 ] = 0. If NumberIntervals <= 1 or NumberIntervals >=
 *   TimeHorizon, then the mapping clearly does not require "ChangeIntervals",
 *   which in fact is not loaded.
 *
 * - The variable "MaxPower", of type netCDF::NcDouble and indexed over both
 *   dimensions "NumberIntervals" and "NumberArcs". The first dimension may
 *   have either size 1 or size "NumberIntervals" (if "NumberIntervals" is not
 *   provided, then the size can also be "TimeHorizon") whereas the second one
 *   always has size NumberArcs (if it is provided at all). This is meant to
 *   represent the matrix MaxP[ t , l ] which, for each time instant t at each
 *   arc l contains the maximum power value of the unit; it must be that MinP[
 *   t , l ] <= MaxP[ t , l ] for each time instant t and each arc l. This
 *   variable is optional; if it is not provided then it is assumed that MaxP[
 *   t , l ] == 0, i.e., the maximum power of all units is zero (i.e., all
 *   units are pumps). If the first dimension has size 1 then the entry MaxP[
 *   0 , l ] is assumed to contain the maximum power of arc l for all time
 *   instants. Otherwise, MaxPower[ i , l ] is the fixed value of MaxP[ t , l
 *   ] for arc l and all time t in the interval [ ChangeIntervals[ i - 1 ] ,
 *   ChangeIntervals[ i ] ], with the assumption that ChangeIntervals[ - 1 ] =
 *   0. If NumberIntervals <= 1 or NumberIntervals >= TimeHorizon, then the
 *   mapping clearly does not require "ChangeIntervals", which in fact is not
 *   loaded.
 *
 * - The variable "DeltaRampUp", of type netCDF::NcDouble and indexed over
 *   both dimensions "NumberIntervals" and "NumberArcs". The first dimension
 *   may have either size 1 or size "NumberIntervals" (if "NumberIntervals" is
 *   not provided, then the size can also be "TimeHorizon"), whereas the
 *   second one always has size NumberArcs (if it is provided at all). This is
 *   meant to represent the matrix DP[ t , l ] which contains the maximum
 *   possible increase of the flow rate at time instant t for arc l. This
 *   variable is optional; if it is not provided then it is assumed that DP[ t
 *   , l ] == MaxP[ t , l ] - MinP[ t , l ], i.e., all units can ramp up by an
 *   arbitrary amount, i.e., there are no ramp-up constraints. If the first
 *   dimension has size 1 then the entry DP[ 0 , l ] is assumed to contain the
 *   maximum possible increase of the flow rate of arc l for all time
 *   instants. Otherwise, DeltaRampUp[ i , l ] is the fixed value of DP[ t , l
 *   ] for arc l and all time t in the interval [ ChangeIntervals[ i - 1 ] ,
 *   ChangeIntervals[ i ] ], with the assumption that ChangeIntervals[ - 1 ] =
 *   0. If NumberIntervals <= 1 or NumberIntervals >= TimeHorizon, then the
 *   mapping clearly does not require "ChangeIntervals", which in fact is not
 *   loaded.
 *
 * - The variable "DeltaRampDown", of type netCDF::NcDouble and indexed over
 *   both dimensions "NumberIntervals" and "NumberArcs". The first dimension
 *   may have either size 1 or size "NumberIntervals" (if "NumberIntervals" is
 *   not provided, then the size can also be "TimeHorizon"), whereas the
 *   second one always has size NumberArcs (if it is provided at all). This is
 *   meant to represent the matrix DM[ t , l ] which contains the maximum
 *   possible decrease of the flow rate at each time instant t of each arc
 *   l. This variable is optional; if it is not provided then it is assumed
 *   that DM[ t , l ] == MaxP[ t , l ] - MinP[ t , l ], i.e., the unit can
 *   ramp down by an arbitrary amount, i.e., there are no ramp-down
 *   constraints. If first dimension has size 1 then the entry DM[ 0 , l ] is
 *   assumed to contain the maximum possible decrease of the flow rate of arc
 *   l for all time instants. Otherwise, DeltaRampDown[ i , l ] is the fixed
 *   value of DM[ t , l ] for arc l and all time t in the interval [
 *   ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the assumption
 *   that ChangeIntervals[ - 1 ] = 0. If NumberIntervals <= 1 or
 *   NumberIntervals >= TimeHorizon, then the mapping clearly does not require
 *   "ChangeIntervals", which in fact is not loaded.
 *
 * - The variable "PrimaryRho", of type netCDF::NcDouble and indexed both over
 *   the dimensions "NumberIntervals" and "NumberArcs". The first dimension
 *   may have either size 1 or size "NumberIntervals" (if "NumberIntervals" is
 *   not provided, then the size can also be "TimeHorizon"), whereas the
 *   second one always has size NumberArcs (if it is provided at all). This is
 *   meant to represent the matrix PR[ t , l ] which, for each time instant t
 *   and arc l, contains the maximum possible fraction of active power that
 *   can be used as primary reserve. This variable is optional, when it's not
 *   present then PR[ t , l ] == 0 for all t and l, i.e., the unit is not
 *   capable of producing any primary reserve. Note that only turbines can
 *   produce primary reserve, i.e., PR[ t , l ] > 0 ==> MaxP[ t , l ] > 0. If
 *   the first dimension has size 1 then the entry PR[ 0 , l ] is assumed to
 *   contain the maximum possible fraction of active power that can be used as
 *   primary reserve by arc l for all time instants. Otherwise, PrimaryRho[ i
 *   , l ] is the fixed value of PR[ t , l ] for arc l and all t in the
 *   interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the
 *   assumption that ChangeIntervals[ - 1 ] = 0 and all l. If NumberIntervals
 *   <= 1 or NumberIntervals >= TimeHorizon, then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded.
 *
 * - The variable "SecondaryRho", of type netCDF::NcDouble and indexed both
 *   over the dimensions "NumberIntervals" and "NumberArcs". The first
 *   dimension may have either size 1 or size "NumberIntervals" (if
 *   "NumberIntervals" is not provided, then the size can also be
 *   "TimeHorizon"), whereas the second one always has size NumberArcs (if it
 *   is provided at all). This is meant to represent the matrix SR[ t , l ]
 *   which, for each time instant t and arc l contains the maximum possible
 *   fraction of active power that can be used as secondary reserve. This
 *   variable is optional, when it's not present then SR[ t , l ] == 0 for all
 *   t and l, i.e., the unit is not capable of producing any secondary
 *   reserve. Note that only turbines can produce secondary reserve, i.e., SR[
 *   t , l ] > 0 ==> MaxP[ t , l ] > 0. If the first dimension has size 1
 *   then the entry SR[ 0 , l ] is assumed to contain the maximum possible
 *   fraction of active power that can use as secondary reserve by arc l for
 *   all time instant. Otherwise, SecondaryRho[ i , l ] is the fixed value of
 *   SR[ t , l ] for arc l and all t in the interval [ ChangeIntervals[ i - 1
 *   ] , ChangeIntervals[ i ] ], with the assumption that ChangeIntervals[ - 1
 *   ] = 0. If NumberIntervals <= 1 or NumberIntervals >= TimeHorizon, then
 *   the mapping clearly does not require "ChangeIntervals" which in fact is
 *   not loaded.
 *
 * - The variable "NumberPieces", of type netCDF::NcUint and indexed over the
 *   dimension "NumberArcs". NumberPieces[ l ] tells how many pieces the
 *   concave flow-to-active-power function has for unit (arc) l. Note that
 *   pumps must necessarily have exactly one piece. The sum over all i of
 *   NumberPieces[ i ] is the total number of pieces (say,
 *   "TotalNumberPieces"). Clearly, TotalNumberPieces >= NumberArcs; if the
 *   flow-to-active-power function for all turbines only have one piece (those
 *   of pumps necessarily are so), then TotalNumberPieces == NumberArcs and
 *   there is no need to define this variable. If, instead, if it is defined,
 *   then should always be such that TotalNumberPieces >= NumberArcs.
 *
 * - The variable "LinearTerm", of type netCDF::NcDouble and indexed over the
 *   set { 0 , ..., TotalNumberPieces - 1 } (see "NumberPieces"). LinearTerm[
 *   h ] gives the linear term a_h of the linear function a_h * f + b_h that
 *   defines the concave flow-to-active-power function for some unit; the
 *   total function if F2AP( t ) = min { a_h * f + b_h , h \in H } for some
 *   finite set H that depends on the individual unit. It is then necessary to
 *   be able to assign a unique index h = 0, 1, ..., TotalNumberPieces - 1 to
 *   each pair ( unit , linear function a_h * f + b_h ). When
 *   TotalNumberPieces == NumberArcs, the index is the same as i = 0, 1, ...,
 *   NumberArcs - 1 (there is a one-to-one correspondence between each (unit)
 *   arc and each piece). When, instead, TotalNumberPieces > NumberArcs, a
 *   mapping must be defined. The mapping is the obvious one: each index of i
 *   = 0, 1, ..., NumberArcs - 1, corresponds with a unit (arc), and the
 *   linear functions for each unit (arc) also have some natural ordering.
 *   Thus, in general the mapping is: piece 0 = first piece of unit (arc) 0
 *   piece 1 = second piece of unit (arc) 0 ... piece NumberPieces[ 0 ] - 1 =
 *   last piece of unit (arc) 0 piece NumberPieces[ 0 ] = first piece of unit
 *   (arc) 1 piece NumberPieces[ 0 ] + 1 = second piece of unit (arc) 1 ...
 *   which of course boils down to "h = i" when each arc has exactly one
 *   piece.
 *
 * - The variable "ConstantTerm", of type netCDF::NcDouble and indexed over
 *   the set { 0 , ..., TotalNumberPieces" - 1 }. ConstantTerm[ h ] gives the
 *   constant term b_h of the linear function a_h * f + b_h that defines the
 *   concave flow-to-active-power function for some unit; see the comments to
 *   "LinearTerm" for details.
 *
 * - The variable "InertiaPower", of type netCDF::NcDouble and indexed both
 *   over the dimensions "NumberIntervals" and "NumberArcs". The first
 *   dimension may have either size 1 or size "NumberIntervals" (if
 *   "NumberIntervals" is not provided, then the size can also be
 *   "TimeHorizon") whereas the second one always has size NumberArcs (if it
 *   is provided at all). This is meant to represent the matrix IP[ t , l ]
 *   which, for each time instant t and arc l, contains the contribution that
 *   the unit can give to the inertia constraint which depends on the active
 *   power that it is currently generating (basically, the constant to be
 *   multiplied to the active power variable) at time t for arc l. The
 *   variable is optional; if it is not defined, IP[ t , l ] == 0 for each
 *   time instants t and arc l. If the first dimension has size 1 then the
 *   entry IP[ 0 , l ] is assumed to contain the inertia power value for
 *   arc l and all time instants t. Otherwise, InertiaPower[ i , l ] is the
 *   fixed value of IP[ t , l ] for all t in the interval [ ChangeIntervals[ i
 *   - 1 ] , ChangeIntervals[ i ] ], with the assumption that ChangeIntervals[
 *   - 1 ] = 0 and all l. If NumberIntervals <= 1 or NumberIntervals >=
 *   TimeHorizon then the mapping clearly does not require "ChangeIntervals",
 *   which in fact is not loaded.
 *
 * - The variable "InitialFlowRate", of type netCDF::NcDouble and indexed over
 *   the dimension "NumberArcs". Each entry InFR[ i ] indicates the amount of
 *   the flow that was going along arc i at time instant -1. This is necessary
 *   to compute ramp-up and ramp-down limits (cf. "DeltaRampUp" and
 *   "DeltaRampDown"), and therefore it is useless if there are no ramp
 *   constraints on *any* unit (arc), in which case it is not loaded.
 *
 * - The variable "InitialVolumetric", of type netCDF::NcDouble and indexed
 *   over the dimension "NumberReservoirs". Each entry InV[ r ] indicates the
 *   volumes of water in reservoir r at time instant -1.
 *
 * - The negative or positive scalar variable "UphillFlow", of type
 *   netCDF::NcInt and indexed over the dimension "NumberArcs". Each entry
 *   UpF[ l ] indicates the uphill flow delay for each unit (arc) l; the
 *   nontrivial concept is detailed below. This variable is optional, if it is
 *   not provided it is taken to be UpF[ l ] == 0.
 *
 * - The positive scalar variable "DownhillFlow", of type netCDF::NcUint and
 *   indexed over the dimension "NumberArcs". Each entry DnF[ l ] indicates
 *   the downhill flow delay for each arc (unit) l. This variable is optional,
 *   if it is not provided it is taken to be DnF[ l ] == 0.
 *
 * The last two quantities require some comment. Let us assume that we have
 * any arc l, with ( StartArc[ l ] = n , EndArc[ l ] = n' ), which
 * corresponds to (say) a turbine. This means that a positive flow along l
 * implies that water is being taken away from n and delivered to n',
 * passing through the turbine to produce flow, as graphically depicted
 * below:
 *
 *   \ n / >==============> [ TURBINE ] >================> \ n' /
 *              UpF[ l ]                     DnF[ l ]
 *
 * The issue here is that the turbine can be geographically far enough from
 * both n and n' that the water can take a long time (one or more time
 * instant, especially if these are "short" such as 15 or 5 minutes) to
 * reach the turbine from n, and n' from the turbine.
 *
 * A particular note of caution has to be mentioned regarding the fact that
 * UpF[ l ] can be *negative*: this means that the water is used in the
 * turbine *before* it goes out of reservoir n. This is counter-intuitive,
 * but can be explained by the fact that the pipe between n and the turbine
 * can be full, and therefore works as a "mini reservoir" in itself. When
 * the turbine is started, a "bubble" (depression) is created uphill the
 * turbine; this "bubbles up" the n ==> TURBINE pipe until it reaches n, and
 * it is only at that point that the water in n starts flowing away. Hence,
 * there is a negative temporal delay between the water starting flowing in
 * the turbine and it starting flowing away from n. Note that if the
 * n ==> TURBINE pipe is rather empty, the delay is positive in that one
 * have to start sending the water, which may take some time before filling
 * the pipe and therefore starting the turbine. Of course these are all
 * somewhat crude approximations of the true physical behavior, but they
 * are accurate enough for this setting. Yet, the case UpF[ l ] < 0 cannot
 * be disregarded. */

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
 /// generate the abstract variables of the HydroUnitBlock
 /** The HydroUnitBlock class has five boost::multi_array< ColVariable , 2 >
  * variables where the first four are:
  *
  * - the primary spinning reserve variables;
  *
  * - the secondary spinning reserve variables;
  *
  * - the active power variables;
  *
  * - the flow rate variables
  *
  * Each of the boost::multi_array< ColVariable , 2 > has as first dimension
  * the time horizon and as second dimension the number of arcs (or
  * generators, which is returned by get_number_generators()). The last
  * boost::multi_array< ColVariable , 2 > variable is:
  *
  * - the volumetric variables
  *
  * Where its' first dimension is the number of reservoirs and the second
  * dimension is the time horizon. */

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the static constraint of the HydroUnit
 /** This method generates the static constraint of the HydroUnitBlock.
  * In order to describe a hydro generating unit system, it will be convenient
  * to see a cascading system as a graph. Let \f$ \mathcal{N}^{hy}\f$ be the
  * set of reservoirs (nodes) and \f$ \mathcal{L}^{hy}\f$ be the set of arcs
  * connecting these reservoirs respectively. Attached to each
  * \f$ l \in \mathcal{L}^{hy}\f$ are one or several plants (turbines or pumps).
  * This system is described on a discrete time horizon as dictated by the
  * UnitBlock interface. In this description we indicate it with
  * \f$ \mathcal{T}=\{ 0, \dots , \mathcal{|T|} - 1\} \f$. Each reservoir
  * \f$ n \in \mathcal{N}^{hy}\f$ has a continuous volumetric variables
  * \f$ v^{hy}_{n,t}\f$ in \f$ m^3 \f$ for \f$ t \in \mathcal{T}\f$ with
  * associated lower and upper bounds \f$ V^{hy,mn}_{n,t}\f$,
  * \f$ V^{hy,mx}_{n,t}\f$ and inflows \f$ A_{n,t}\f$ in \f$ m^3 /s \f$. The
  * uphill and downhill flow rate are defined as \f$ \tau^{up} \f$ and
  * \f$ \tau^{dn}\f$ respectively. For each time \f$ t \in \mathcal{T}\f$ and
  * each arc \f$ l \in \mathcal{L}^{hy}\f$ the continuous flow rate variable
  * \f$ f_{t,l} \f$ in \f$ m^3 /s \f$ and ramping conditions
  * \f$ \Delta^{up}_{t,l} \f$ and \f$ \Delta^{dn}_{t,l} \f$ in
  * \f$ (m^3 /s)/h \f$ are disposed. The flow rate variable will be subject to
  * bounds \f$ F^{mn}_{t,l} \f$ and \f$ F^{mx}_{t,l} \f$ and it's assumed
  * moreover given a cutting plane model describing power as a function of flow
  * rate as below:
  *
  * \f[
  *  p^{ac}_{t,l}(f) := min \{ P_l + \rho^{hy}_{l}f_{t,l}\}
  * \f]
  *
  * where \f$ P_l \f$ and \f$ \rho^{hy}_{l} \f$ are considered as the constant
  * and linear multipliers of the linear function (flow-to-active-power)
  * respectively. Power generated by the hydro unit in each time and for each
  * arc \f$ p^{ac}_{t,l}, p^{pr}_{t,l}, p^{sc}_{t,l} \f$ in MW will be subject
  * to bounds \f$ P^{mn}_{t,l} \f$ and \f$ P^{mx}_{t,l} \f$ respectively.
  * Besides, we emphasize that reserve requirements are specified in order to
  * be symmetrically available to increase or decrease power injected into the
  * grid. For some of the constraints we will need to distinguish between pumps
  * and turbines. The distinction is made by considering the set of feasible
  * flow rates. Whenever \f$ [ F^{mn}_{t,l} , F^{mx}_{t,l}] \subseteq R_- \f$
  * for each arc and each time, the unit is considered a pump, and whenever
  * \f$ [ F^{mn}_{t,l} , F^{mx}_{t,l}] \subseteq R_+ \f$ the unit is considered
  * a turbine. Any possible mixed situation can be accounted for by
  * artificially splitting the unit into “two units”, which should be done at
  * the data processing stage (see deserialize() comments). With above
  * description the mathematical constraint of hydro unit may present as below:
  *
  * - maximum and minimum power output constraints according to primary and
  *   secondary spinning reserves are are presented in (1)-(2). Each of them
  *   is a boost::multi_array< FRowConstraint , 2 >; with two dimensions which
  *   are
  *   f_time_horizon, and f_NumberArcs entries, where the entry
  *   t = 0, ...,f_time_horizon - 1 and the entry l = 0, ...,f_NumberArcs - 1
  *   being the maximum and minimum power output value according to the primary
  *   and the secondary spinning reserves at time t and arc l. these ensure the
  *   maximum (or minimum) amount of energy that unit can produce (or use) when
  *   it is on (or off).
  *
  *   \f[
  *      p^{ac}_{t,l} + p^{pr}_{t,l} + p^{sc}_{t,l} \leq P^{mx}_{t,l}
  *                   \quad t \in \mathcal{T}, l \in \mathcal{L}^{hy} \quad (1)
  *   \f]
  *
  *   \f[
  *     P^{mn}_{t,l} \leq p^{ac}_{t,l} - p^{pr}_{t,l} - p^{sc}_{t,l}
  *                   \quad t \in \mathcal{T}, l \in \mathcal{L}^{hy} \quad (2)
  *   \f]
  *
  * - primary and secondary spinning reserves relation with active power at
  *   each time and for each turbine: the same as inequalities (1)-(2), the
  *   inequalities (3)-(4) ensure that maximum amount of primary and secondary
  *   spinning reserve in the problem. Each of them is a
  *   boost::multi_array< FRowConstraint , 2 >; with two dimensions which are
  *   f_time_horizon, and f_NumberArcs entries, where
  *   t = 0, ...,f_time_horizon - 1 and l = 0, ...,f_NumberArcs - 1
  *
  *   \f[
  *      p^{pr}_{t,l} \leq \rho^{pr}_{t,l}p^{ac}_{t,l} \quad t \in \mathcal{T},
  *        l \in \mathcal{L}^{hy} \quad with
  *        \quad [ F^{mn}_{t,l} , F^{mx}_{t,l}] \subseteq R_+         \quad (3)
  *   \f]
  *
  *   \f[
  *      p^{sc}_{t,l} \leq \rho^{sc}_{t,l}p^{ac}_{t,l} \quad t \in \mathcal{T},
  *        l \in \mathcal{L}^{hy} \quad with
  *        \quad [ F^{mn}_{t,l} , F^{mx}_{t,l}] \subseteq R_+         \quad (4)
  *   \f]
  *
  * where \f$ \rho^{pr}_{t,l} \f$ and \f$ \rho^{sc}_{t,l}\f$ are the maximum
  * possible fraction of active power at each time and each arc that can be
  * used as primary and secondary reserve respectively.
  *
  * - primary and secondary spinning reserves at each time and for each pump:
  *   these equalities (5)-(6) ensure that the primary and secondary spinning
  *   reserve value for each pump is equal to zero. Each of them is a
  *   boost::multi_array< FRowConstraint , 2 >; with two dimensions which are
  *   f_time_horizon, and f_NumberArcs entries, where
  *   t = 0, ...,f_time_horizon - 1 and l = 0, ...,f_NumberArcs - 1
  *
  *   \f[
  *      p^{pr}_{t,l} = 0 \quad t \in \mathcal{T},
  *        l \in \mathcal{L}^{hy} \quad with
  *        \quad [ F^{mn}_{t,l} , F^{mx}_{t,l}] \subseteq R_-         \quad (5)
  *   \f]
  *
  *   \f[
  *      p^{sc}_{t,l} = 0 \quad t \in \mathcal{T},
  *        l \in \mathcal{L}^{hy} \quad with
  *        \quad [ F^{mn}_{t,l} , F^{mx}_{t,l}] \subseteq R_-         \quad (6)
  *   \f]
  *
  * - flow to active power function at each time and for each pump: this
  *   equality (7) gives the active power relation with flow rate for each
  *   pump at time t. This is a boost::multi_array< FRowConstraint , 2 >; with
  *   two dimensions which are f_time_horizon, and f_NumberArcs entries, where
  *   t = 0, ...,f_time_horizon - 1 and l = 0, ...,f_NumberArcs - 1
  *
  *   \f[
  *      p^{ac}_{t,l} = \rho^{hy}_{l}f_{t,l} \quad t \in \mathcal{T},
  *        l \in \mathcal{L}^{hy} \quad with
  *        \quad [ F^{mn}_{t,l} , F^{mx}_{t,l}] \subseteq R_-         \quad (7)
  *   \f]
  *
  * - flow-to-active-power function at each time and for each turbine ;
  *
  *   \f[
  *      p^{ac}_{t,l} \leq P_j + \rho^{hy}_{j}f_{t,l} \quad j \in \mathcal{J}_l
  *        \quad t \in \mathcal{T}, l \in \mathcal{L}^{hy} \quad with
  *        \quad [ F^{mn}_{t,l} , F^{mx}_{t,l}] \subseteq R_+         \quad (8)
  *   \f]
  *
  * - flow rate variable bounds: This inequality (9) indicates upper and lower
  *   bound of flow rate at time t and for ach arc l, thus that is a
  *   boost::multi_array< FRowConstraint , 2 >; with two dimensions which are
  *   f_time_horizon, and f_NumberArcs entries, where
  *   t = 0, ...,f_time_horizon - 1 and l = 0, ...,f_NumberArcs - 1
  *
  *   \f[
  *      f_{t,l} \in [ F^{mn}_{t,l} , F^{mx}_{t,l}]  \quad t \in \mathcal{T},
  *           l \in \mathcal{L}^{hy}                                  \quad (9)
  *   \f]
  *
  * - ramp-up and ramp-down constraints: These inequality (10)-(11) indicate
  *   ramp-up and ramp-down constraints at time t and for ach arc l, so each of
  *   them is a boost::multi_array< FRowConstraint , 2 >; with two dimensions
  *   which are f_time_horizon, and f_NumberArcs entries, where
  *   t = 0, ...,f_time_horizon - 1 and l = 0, ...,f_NumberArcs - 1
  *
  *   \f[
  *      f_{t,l} - f_{t-1,l} \leq \Delta^{up}_{t,l} \quad t \in \mathcal{T},
  *           l \in \mathcal{L}^{hy}                                 \quad (10)
  *   \f]
  *
  *   \f[
  *      f_{t-1,l} - f_{t,l} \leq \Delta^{dn}_{t,l} \quad t \in \mathcal{T},
  *           l \in \mathcal{L}^{hy}                                 \quad (11)
  *   \f]
  *
  * - final volumes of each reservoir constraints: this equality (12) gives the
  *   final volumes of each reservoir r at time t. This is a
  *   boost::multi_array< FRowConstraint , 2 >; with two dimensions which are
  *   f_NumberReservoirs, and f_time_horizon entries, where
  *   n = 0, ...,f_NumberReservoirs - 1 and t = 0, ...,f_time_horizon - 1
  *
  *   \f[
  *      v^{hy}_{n,t} = v^{hy}_{n,t-1} + A_{n,t} + (\sum_{l=(n',n) \in
  *      \mathcal{L}^{hy} } f_{t - \tau^{dn}_l , l } - \sum_{l=(n,n') \in
  *      \mathcal{L}^{hy} } f_{t - \tau^{up}_l , l }) \quad t \in \mathcal{T},
  *      \quad n \in \mathcal{N}^{hy}                                \quad (12)
  *   \f]
  *
  *   where in each arc \f$ l=(n,n') \in \mathcal{L}^{hy} \f$, \f$ n \f$ and
  *   \f$ n' \f$ are supposed to be the start and the end point of that
  *   respectively.
  *
  * - final volumes variable bounds: This inequality (13) indicates upper and
  *   lower bound of volumetric variables of each reservoir for each time t,
  *   thus that is a boost::multi_array< FRowConstraint , 2 >; with two
  *   dimensions which are f_NumberReservoirs, and f_time_horizon entries,
  *   where n = 0, ...,f_NumberReservoirs - 1 and t = 0, ...,f_time_horizon - 1
  *
  *   \f[
  *      v^{hy}_{n,t} \in [ V^{hy,mn}_{n,t} , V^{hy,mx}_{n,t}]  \quad
  *        n \in \mathcal{N}^{hy}, t \in \mathcal{T}                 \quad (13)
  *   \f]
  */

 void generate_abstract_constraints( Configuration * stcc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the objective of the HydroUnitBlock
 /** Method that generates the objective of the HydroUnitBlock.
  *
  * - Objective function: the objective function of the HydroUnitBlock is
  *   "empty" (a FRealObjective with a LinearFunction inside with no active
  *   variables) */

 void generate_objective( Configuration * objc = nullptr ) override;

/**@} ----------------------------------------------------------------------*/
/*---------------- Methods for checking the HydroUnitBlock -----------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for checking solution information in the HydroUnitBlock
 * @{ */

 /// returns true if the current solution is (approximately) feasible
 /** This function returns true if and only if the solution encoded in the
  * current value of the Variable of this HydroUnitBlock is approximately
  * feasible within the given tolerance. That is, a solution is considered
  * feasible if and only if
  *
  * -# each ColVariable is feasible; and
  *
  * -# the violation of each Constraint of this HydroUnitBlock is not
  *      greater than the tolerance.
  *
  * Every Constraint of this HydroUnitBlock is a RowConstraint and its
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

/**@} ----------------------------------------------------------------------*/
/*--------- METHODS FOR READING THE DATA OF THE HydroUnitBlock -------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the HydroUnitBlock
 *
 * These methods allow to read data that must be common to (in principle) all
 * the kind of hydro units
 * @{ */

 /// returns the number of reservoirs
 Index get_number_reservoirs( void ) const {
  return( f_NumberReservoirs ? f_NumberReservoirs : 1 );
 }

/*--------------------------------------------------------------------------*/
 /// returns the number of arcs/generators

 Index get_number_generators( void ) const override {
  return( f_NumberArcs ? f_NumberArcs : 1 );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of start arcs
 /** Method for returning the vector of starting point of each arc. This vector
  * may have size of 1 (single hydro unit with just one arc between two
  * reservoirs) or the size of number of reservoirs, then there are two
  * possible cases:
  *
  * - if f_NumberReservoirs == 2, this vector has size of 1 which means there
  *   is one arc at the system (single hydro case).
  *
  * - if f_NumberReservoirs > 2, this vector have size of f_NumberReservoirs
  *   and each element of the vectors gives starting point of each arc in the
  *   network. */

 const std::vector< Index > & get_start_arc( void ) const {
  return( v_StartArc );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of end arcs
 /** Method for returning the vector of ending point of each arc. This
  * vector may have size of 1 (single hydro unit with just one arc between two
  * reservoirs) or the size of number of reservoirs, then there are two
  * possible cases:
  *
  * - if f_NumberReservoirs == 2, this vector has size of 1 which means there
  *   is one arc at the system (single hydro case).
  *
  * - if f_NumberReservoirs > 2, this vector have size of f_NumberReservoirs
  *   and each element of the vectors gives ending point of each arc in the
  *   network. */

 const std::vector< Index > & get_end_arc( void ) const { return( v_EndArc ); }

/*--------------------------------------------------------------------------*/
 /// returns the matrix of inertia power
 /** The returned value U = get_inertia_power() contains the contribution to
  * inertia (basically, the constants to be multiplied by the active power
  * variables returned by get_active_power()) of eac arcs (generators) at each
  * time instants. There are three possible cases:
  *
  * - if the matrix is empty, then the inertia power is 0 and this function
  *   returns a nullptr;
  *
  * - if the matrix only has one row (i.e., the first dimension has size 1),
  *   then the inertia power for arc (generator) l is U[ 0 , l ] for all t
  *   which means that the second dimension has size get_number_arcs();
  *
  * - otherwise, the matrix has size the time horizon per
  *   get_number_arcs(), and U[ t , l ] contains the contribution to
  *   inertia power of arc (generator) l at time instant t. */

 const double * get_inertia_power( Index generator ) const override {
  if( v_InertiaPower.empty() )
   return( nullptr );
  return( v_InertiaPower.data() + generator * f_time_horizon );
 }

/*--------------------------------------------------------------------------*/
 /// returns the matrix of minimum volumetric
 /** The method returned a two-dimensional boost::multi_array<> M such that
  * M[ n , t ] gives the minimum volumetric of the reservoir n at the time
  * instant t. This two-dimensional boost::multi_array<> M considers three
  * possible cases:
  *
  * - if the boost::multi_array<> M is empty() then no minimum volumetric are
  *   defined, and there are no minimum volumetric constraints;
  *
  * - if the boost::multi_array<> M has only one column which in this case the
  *   boost::multi_array<> M is a (transpose of) vector with size
  *   get_number_reservoirs(). Each element of M[ n , 0 ] gives the minimum
  *   volumetric of the each reservoir n for all time instant t;
  *
  * - otherwise the two-dimensional boost::multi_array<> M must have
  *   get_number_reservoirs() row where each row must have size of
  *   get_time_horizon() and each element of M[ n , t ] gives the minimum
  *   volumetric of reservoir n at time instant t. */

 const boost::multi_array< double , 2 > & get_min_volumetric( void ) const {
  return( v_MinVolumetric );
 }

/*--------------------------------------------------------------------------*/
 /// returns the matrix of maximum volumetric
 /** The method returned a two-dimensional boost::multi_array<> M such that
  * M[ n , t ] gives the maximum volumetric of the reservoir n at the time
  * instant t. This two-dimensional boost::multi_array<> M considers three
  * possible cases:
  *
  * - if the boost::multi_array<> M is empty() then no maximum volumetric are
  *   defined, and there are no minimum volumetric constraints;
  *
  * - if the boost::multi_array<> M has only one column which in this case the
  *   boost::multi_array<> M is a (transpose of) vector with size
  *   get_number_reservoirs(). Each element of M[ n , 0 ] gives the maximum
  *   volumetric of the each reservoir n for all time instant t;
  *
  * - otherwise the two-dimensional boost::multi_array<> M must have
  *   get_number_reservoirs() row where each row must have size of
  *   get_time_horizon() and each element of M[ n , t ] gives the maximum
  *   volumetric of reservoir n at time instant t. */

 const boost::multi_array< double , 2 > & get_max_volumetric( void ) const {
  return( v_MaxVolumetric );
 }

/*--------------------------------------------------------------------------*/
 /// returns the matrix of inflows
 /** The method returned a two-dimensional boost::multi_array<> M such that
  * M[ n , t ] gives the inflows of the reservoir n at the time
  * instant t. This two-dimensional boost::multi_array<> M considers two
  * possible cases:
  *
  * - if the boost::multi_array<> M is empty() then no inflows are defined;
  *
  * - otherwise the two-dimensional boost::multi_array<> M must have
  *   get_number_reservoirs() row where each row must have size of
  *   get_time_horizon() and each element of M[ n , t ] represents the inflows
  *   of reservoir n at time instant t. */

 const boost::multi_array< double , 2 > & get_inflows( void ) const {
  return( v_inflows );
 }

/*--------------------------------------------------------------------------*/
 /// returns the matrix of minimum power
 /** The method returned a two-dimensional boost::multi_array<> M such that
  * M[ t , i ] gives the minimum power at each time t associated with unit (arc)
  * i. This two-dimensional boost::multi_array<> M considers three possible
  * cases:
  *
  * - if the boost::multi_array<> M is empty() then no minimum power are
  *   defined, and there are no minimum power constraints;
  *
  * - if the boost::multi_array<> M has only one row which in this case the
  *   boost::multi_array<> M is a vector with size get_number_arcs(). Each
  *   element of M[ 0 , i ] gives the minimum power for all time instant t of
  *   each unit i;
  *
  * - otherwise the two-dimensional boost::multi_array<> M must have
  *   get_time_horizon() rows and get_number_arcs() columns and each element
  *   of M[ t , i ] gives the minimum power at time t and unit i. */

 double get_min_power( Index t , Index generator = 0 ) const override {
  return( *( v_MinPower.data() + generator * f_time_horizon + t ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the matrix of maximum power
 /** The method returned a two-dimensional boost::multi_array<> M such that
  * M[ t , i ] gives the maximum power at each time t associated with unit (arc)
  * i. This two-dimensional boost::multi_array<> M considers three possible
  * cases:
  *
  * - if the boost::multi_array<> M is empty() then no maximum power are
  *   defined, and there are no maximum power constraints;
  *
  * - if the boost::multi_array<> M has only one row which in this case the
  *   boost::multi_array<> M is a vector with size get_number_arcs(). Each
  *   element of M[ 0 , i ] gives the maximum power for all time instant t of
  *   each unit i;
  *
  * - otherwise the two-dimensional boost::multi_array<> M must have
  *   get_time_horizon() rows and get_number_arcs() columns and each element
  *   of M[ t , i ] gives the maximum power at time t and unit i. */

 double get_max_power( Index t , Index generator = 0 ) const override {
  return( *(v_MaxPower.data() + generator * f_time_horizon + t) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the matrix of minimum flow
 /** The method returned a two-dimensional boost::multi_array<> M such that
  * M[ t , i ] gives the minimum flow at each time t associated with unit (arc)
  * i. This two-dimensional boost::multi_array<> M considers three possible
  * cases:
  *
  * - if the boost::multi_array<> M is empty() then no minimum flow are
  *   defined, and there are no minimum flow constraints;
  *
  * - if the boost::multi_array<> M has only one row which in this case the
  *   boost::multi_array<> M is a vector with size get_number_arcs(). Each
  *   element of M[ 0 , i ] gives the minimum flow for all time instant t of
  *   each unit i;
  *
  * - otherwise the two-dimensional boost::multi_array<> M must have
  *   get_time_horizon() rows and get_number_arcs() columns and each element
  *   of M[ t , i ] gives the minimum flow at time t and unit i. */

 const boost::multi_array< double , 2 > & get_min_flow( void ) const {
  return( v_MinFlow );
 }

/*--------------------------------------------------------------------------*/
 /// returns the matrix of maximum flow
 /** The method returned a two-dimensional boost::multi_array<> M such that
  * M[ t , i ] gives the maximum flow at each time t associated with unit (arc)
  * i. This two-dimensional boost::multi_array<> M considers three possible
  * cases:
  *
  * - if the boost::multi_array<> M is empty() then no maximum flow are
  *   defined, and there are no maximum flow constraints;
  *
  * - if the boost::multi_array<> M has only one row which in this case the
  *   boost::multi_array<> M is a vector with size get_number_arcs(). Each
  *   element of M[ 0 , i ] gives the maximum flow for all time instant t of
  *   each unit i;
  *
  * - otherwise the two-dimensional boost::multi_array<> M must have
  *   get_time_horizon() rows and get_number_arcs() columns and each element
  *   of M[ t , i ] gives the maximum flow at time t and unit i. */

 const boost::multi_array< double , 2 > & get_max_flow( void ) const {
  return( v_MaxFlow );
 }

/*--------------------------------------------------------------------------*/
 /// returns the matrix of delta ramp up
 /** The method returned a two-dimensional boost::multi_array<> M such that
  * M[ t , i ] gives the delta ramp up at each time t associated with unit (arc)
  * i. This two-dimensional boost::multi_array<> M considers three possible
  * cases:
  *
  * - if the boost::multi_array<> M is empty() then no ramping constraints;
  *
  * - if the boost::multi_array<> M has only one row which in this case the
  *   boost::multi_array<> M is a vector with size get_number_arcs(). Each
  *   element of M[ 0 , i ] gives the delta ramp up value for all time instant
  *   t of each unit i;
  *
  * - otherwise the two-dimensional boost::multi_array<> M must have
  *   get_time_horizon() rows and get_number_arcs() columns and each element
  *   of M[ t , i ] gives the delta ramp up value at time t and unit i. */

 const boost::multi_array< double , 2 > & get_delta_ramp_up( void ) const {
  return( v_DeltaRampUp );
 }

/*--------------------------------------------------------------------------*/
 /// returns the matrix of delta ramp down
 /** The method returned a two-dimensional boost::multi_array<> M such that
  * M[ t , i ] gives the delta ramp down at each time t associated with unit
  * (arc) i. This two-dimensional boost::multi_array<> M considers three
  * possible cases:
  *
  * - if the boost::multi_array<> M is empty() then no ramping constraints;
  *
  * - if the boost::multi_array<> M has only one row which in this case the
  *   boost::multi_array<> M is a vector with size get_number_arcs(). Each
  *   element of M[ 0 , i ] gives the delta ramp down value for all time
  *   instant t of each unit i;
  *
  * - otherwise the two-dimensional boost::multi_array<> M must have
  *   get_time_horizon() rows and get_number_arcs() columns and each element
  *   of M[ t , i ] gives the delta ramp down value at time t and unit i. */

 const boost::multi_array< double , 2 > & get_delta_ramp_down( void ) const {
  return( v_DeltaRampDown );
 }

/*--------------------------------------------------------------------------*/
 /// returns the matrix of primary rho
 /** The method returned a two-dimensional boost::multi_array<> M such that
  * M[ t , i ] gives the primary rho at each time t associated with unit
  * (arc) i. This two-dimensional boost::multi_array<> M considers three
  * possible cases:
  *
  * - if the boost::multi_array<> M is empty() then no primary reserve
  *   constraints;
  *
  * - if the boost::multi_array<> M has only one row which in this case the
  *   boost::multi_array<> M is a vector with size get_number_arcs(). Each
  *   element of M[ 0 , i ] gives the primary rho value for all time instant t
  *   of each unit i;
  *
  * - otherwise the two-dimensional boost::multi_array<> M must have
  *   get_time_horizon() rows and get_number_arcs() columns and each element
  *   of M[ t , i ] gives the primary rho value at time t and unit i. */

 const boost::multi_array< double , 2 > & get_primary_rho( void ) const {
  return( v_PrimaryRho );
 }

/*--------------------------------------------------------------------------*/
 /// returns the matrix of secondary rho
 /** The method returned a two-dimensional boost::multi_array<> M such that
  * M[ t , i ] gives the secondary rho at each time t associated with unit
  * (arc) i. This two-dimensional boost::multi_array<> M considers three
  * possible cases:
  *
  * - if the boost::multi_array<> M is empty() then no secondary reserve
  *   constraints;
  *
  * - if the boost::multi_array<> M has only one row which in this case the
  *   boost::multi_array<> M is a vector with size get_number_arcs(). Each
  *   element of M[ 0 , i ] gives the secondary rho value for all time instant
  *   t of each unit i;
  *
  * - otherwise the two-dimensional boost::multi_array<> M must have
  *   get_time_horizon() rows and get_number_arcs() columns and each element
  *   of M[ t , i ] gives the secondary rho value at time t and unit i. */

 const boost::multi_array< double , 2 > & get_secondary_rho( void ) const {
  return( v_SecondaryRho );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of number pieces
 /** The returned vector contains the number of pieces for each unit (arc) i.
  * There are three possible cases:
  *
  * - if the vector is empty, then the number of pieces of each unit is 0;
  *
  * - if the vector has only one element, then V[ 0 ] presents the number
  *   of pieces for all units;
  *
  * - otherwise, the vector must have size get_number_arcs() and each element
  *   of V[ i ] represents the number of pieces of each unit i. */

 const std::vector< Index > & get_number_pieces( void ) const {
  return( v_NumberPieces );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of constant term
 /** The returned vector contains the constant term value of each available
  * pieces h in the set of {0, ..., TotalNumberPieces}(see NumberPieces
  * deserialize() comments). There are three possible cases:
  *
  * - if the vector is empty, then the constant term value is 0;
  *
  * - if the vector has only one element, then V[ 0 ] presents the const term
  *   value for all the pieces;
  *
  * - otherwise, the returned V is a std::vector < double > and
  *   V.sized == TotalNumberPieces and each element of V[ h ] represents the
  *   const term value of each piece h. */

 const std::vector< double > & get_const_term( void ) const {
  return( v_ConstTerm );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of linear term
 /** The returned vector contains the linear term value of each available
  * pieces h in the set of {0, ..., TotalNumberPieces}(see NumberPieces
  * deserialize() comments). There are three possible cases:
  *
  * - if the vector is empty, then the linear term value is 0;
  *
  * - if the vector has only one element, then V[ 0 ] presents the linear term
  *   value for all the pieces;
  *
  * - otherwise, the returned V is a std::vector < double > and
  *   V.sized == TotalNumberPieces and each element of V[ h ] represents the
  *   linear term value of each piece h. */

 const std::vector< double > & get_linear_term( void ) const {
  return( v_LinearTerm );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of uphill delay
 /** The returned vector contains the uphill delay for each unit (arc) i.
  * There are three possible cases:
  *
  * - if the vector is empty, then the uphill delay for each unit is 0;
  *
  * - if the vector has only one element, then V[ 0 ] presents the uphill delay
  *   for all units;
  *
  * - otherwise, the vector must have size get_number_arcs() and each element
  *   of V[ i ] represents the uphill delay for each unit i. */

 const std::vector< int > & get_uphill_delay( void ) const {
  return( v_UphillDelay );
 }

/*--------------------------------------------------------------------------*/
 /// returns the uphill delay for the given \p arc
 /** This function returns the uphill delay for the given \p arc.
  *
  * @return The uphill delay for the given \p arc. */

 int get_uphill_delay( Index arc ) const {
  if( v_UphillDelay.empty() )
   return( 0 );
  assert( arc < v_UphillDelay.size() );
  return( v_UphillDelay[ arc ] );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of downhill delay
 /** The returned vector contains the downhill delay for each unit (arc) i.
  * There are three possible cases:
  *
  * - if the vector is empty, then the downhill delay for each unit is 0;
  *
  * - if the vector has only one element, then V[ 0 ] presents the downhill
  *   delay for all units;
  *
  * - otherwise, the vector must have size get_number_arcs() and each element
  *   of V[ i ] represents the downhill delay for each unit i. */

 const std::vector< Index > & get_downhill_delay( void ) const {
  return( v_DownhillDelay );
 }

/*--------------------------------------------------------------------------*/
 /// returns the downhill delay for the given \p arc
 /** This function returns the downhill delay for the given \p arc.
  *
  * @return The downhill delay for the given \p arc. */

 Index get_downhill_delay( Index arc ) const {
  if( v_DownhillDelay.empty() )
   return( 0 );
  assert( arc < v_DownhillDelay.size() );
  return( v_DownhillDelay[ arc ] );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of initial volumetric
 /** The returned vector contains the initial volumetric for each reservoir n.
  * There are three possible cases:
  *
  * - if the vector is empty, then the initial volumetric for each reservoir is
  *   0;
  *
  * - if the vector has only one element, then V[ 0 ] presents the initial
  *   volumetric for all reservoirs;
  *
  * - otherwise, the vector must have size get_number_reservoirs() and each
  *   element of V[ i ] represents the initial volumetric for each reservoir n.
  */

 const std::vector< double > & get_initial_volumetric( void ) const {
  return( v_InitialVolumetric );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of initial flow rate
 /** The returned vector contains the initial flow rate for each unit (arc) i.
  * There are three possible cases:
  *
  * - if the vector is empty, then the initial flow rate for each unit is 0;
  *
  * - if the vector has only one element, then V[ 0 ] presents the initial flow
  *   rate for all units;
  *
  * - otherwise, the vector must have size get_number_arcs() and each element
  *   of V[ i ] represents the initial flow rate for each unit i. */

 const std::vector< double > & get_initial_flow_rate( void ) const {
  return( v_InitialFlowRate );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of initial flow rate at the given \p arc
 /** This method returns the initial flow rate at the given \p arc.
  *
  * @param arc The index of the arc whose initial flow rate is desired.
  *
  * @return The initial flow rate at the given \p arc. */

 double get_initial_flow_rate( Index arc ) const {
  assert( arc < f_NumberArcs );
  return( v_InitialFlowRate.empty() ? 0.0 :
          ( ( v_InitialFlowRate.size() == 1 ) ? v_InitialFlowRate.front() :
            v_InitialFlowRate[ arc ] ) );
 }

/**@} ----------------------------------------------------------------------*/
/*---------- METHODS FOR READING THE Variable OF THE HydroUnitBlock --------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the Variable of the HydroUnitBlock
 *
 * These methods allow to read the five groups of Variable that any
 * HydroUnitBlock in principle has (although some may not):
 *
 * - the volumetric variables;
 *
 * - the flow rate variables;
 *
 * - the active power variables;
 *
 * - the primary spinning reserve variables;
 *
 * - the secondary spinning reserve variables;
 *
 * All these five groups of variables are (if not empty)
 * boost::multi_array< ColVariable , 2 >. The volumetric variables have the
 * number reservoirs as the first dimension and time horizon as the second
 * dimension; all other variables have number of arcs (generators) as the
 * first dimension and time horizon as the second one.
 * @{ */

 /// returns the array of volume variables of the given \p reservoir
 /** This method returns the array of ColVariable representing the volumes of
  * the given \p reservoir at all time instants t in {0, ..., time_horizon -
  * 1}.
  *
  * @param reservoir The index of a reservoir (a number between 0 and
  *                  get_number_reservoirs() - 1).
  *
  * @return The array of ColVariable representing the volumes of the given \p
  *         reservoir. */

 ColVariable * get_volumetric( Index reservoir ) {
  if( reservoir < get_number_reservoirs() ) {
   const auto offset = reservoir * f_time_horizon;
   if( offset < v_volumetric.num_elements() )
    return( v_volumetric.data() + offset );
  }
  return( nullptr );
 }

/*--------------------------------------------------------------------------*/
 /// returns the volume of the given reservoir at the given time
 /** Returns a pointer to the ColVariable representing the volume of the given
  * \p reservoir at the given \p time.
  *
  * @param reservoir The index of the reservoir whose volume is desired (a
  *                  number between 0 and get_number_reservoirs() - 1).
  *
  * @param time The time at which the volume is desired (a number between 0
  *             and get_time_horizon() - 1).
  *
  * @return A pointer to the ColVariable representing the volume of the given
  *         \p reservoir at the given \p time. */

 ColVariable * get_volume( Index reservoir , Index time ) {
  if( reservoir < get_number_reservoirs() && time < f_time_horizon ) {
   const auto offset = reservoir * f_time_horizon + time;
   if( offset < v_volumetric.num_elements() )
    return( v_volumetric.data() + offset );
  }
  return( nullptr );
 }

/*--------------------------------------------------------------------------*/
 /// returns the array of active power variables of the given \p generator
 /** This method returns the array of ColVariable representing the active
  * power of the given \p generator at all time instants t in {0, ...,
  * time_horizon - 1}.
  *
  * @param generator The index of a generator (a number between 0 and
  *                  get_number_generators() - 1).
  *
  * @return The array of ColVariable representing the active power of the
  *         given \p generator. */

 ColVariable * get_active_power( Index generator ) override {
  if( generator < get_number_generators() ) {
   const auto offset = generator * f_time_horizon;
   if( offset < v_active_power.num_elements() )
    return( v_active_power.data() + offset );
  }
  return( nullptr );
 }

/*--------------------------------------------------------------------------*/
 /// returns the active power of the given generator at the given time
 /** Returns a pointer to the ColVariable representing the active power of the
  * given \p generator at the given \p time.
  *
  * @param generator The index of the generator whose active power is desired
  *                  (a number between 0 and get_number_generators() - 1).
  *
  * @param time The time at which the active power is desired (a number
  *             between 0 and get_time_horizon() - 1).
  *
  * @return A pointer to the ColVariable representing the active power of the
  *         given \p generator at the given \p time. */

 ColVariable * get_active_power( Index generator , Index time ) {
  if( generator < get_number_generators() && time < f_time_horizon ) {
   const auto offset = generator * f_time_horizon + time;
   if( offset < v_active_power.num_elements() )
    return( v_active_power.data() + offset );
  }
  return( nullptr );
 }

/*--------------------------------------------------------------------------*/
 /// returns the array of flow rate variables along the given \p arc
 /** This method returns the array of ColVariable representing the flow rate
  * along the given \p arc at all time instants t in {0, ..., time_horizon -
  * 1}.
  *
  * @param arc The index of an arc (a number between 0 and
  *            get_number_generators() - 1).
  *
  * @return The array of ColVariable representing the flow rate along the
  *         given \p arc. */

 ColVariable * get_flow_rate( Index arc ) {
  if( arc < get_number_generators() ) {
   const auto offset = arc * f_time_horizon;
   if( offset < v_flow_rate.num_elements() )
    return( v_flow_rate.data() + offset );
  }
  return( nullptr );
 }

/*--------------------------------------------------------------------------*/
 /// returns the flow rate along the given arc at the given time
 /** Returns a pointer to the ColVariable representing the flow rate along the
  * given \p arc at the given \p time.
  *
  * @param arc The index of the arc whose flow rate is desired (a number
  *            between 0 and get_number_generators() - 1).
  *
  * @param time The time at which the flow rate is desired (a number
  *             between 0 and get_time_horizon() - 1).
  *
  * @return A pointer to the ColVariable representing the flow rate along the
  *         given \p arc at the given \p time. */

 ColVariable * get_flow_rate( Index arc , Index time ) {
  if( arc < get_number_generators() && time < f_time_horizon ) {
   const auto offset = arc * f_time_horizon + time;
   if( offset < v_flow_rate.num_elements() )
    return( v_flow_rate.data() + offset );
  }
  return( nullptr );
 }

/*--------------------------------------------------------------------------*/
 /// returns the array of primary spinning reserve of the given \p generator
 /** This method returns the array of ColVariable representing the primary
  * spinning reserve of the given \p generator at all time instants t in {0,
  * ..., time_horizon - 1}.
  *
  * @param generator The index of a generator (a number between 0 and
  *                  get_number_generators() - 1).
  *
  * @return The array of ColVariable representing the primary spinning
  *         reserve of the given \p generator. */

 ColVariable * get_primary_spinning_reserve( Index generator ) override {
  if( generator < get_number_generators() ) {
   const auto offset = generator * f_time_horizon;
   if( offset < v_primary_spinning_reserve.num_elements() )
    return( v_primary_spinning_reserve.data() + offset );
  }
  return( nullptr );
 }

/*--------------------------------------------------------------------------*/
 /// returns the primary spinning reserve of a generator at the given time
 /** Returns a pointer to the ColVariable representing the primary spinning
  * reserve of the given \p generator at the given \p time.
  *
  * @param generator The index of the generator whose primary spinning reserve
  *                  is desired (a number between 0 and
  *                  get_number_generators() - 1).
  *
  * @param time The time at which the primary spinning reserve is desired (a
  *             number between 0 and get_time_horizon() - 1).
  *
  * @return A pointer to the ColVariable representing the primary spinning
  *         reserve of the given \p generator at the given \p time. */

 ColVariable * get_primary_spinning_reserve( Index generator , Index time ) {
  if( generator < get_number_generators() && time < f_time_horizon ) {
   const auto offset = generator * f_time_horizon + time;
   if( offset < v_primary_spinning_reserve.num_elements() )
    return( v_primary_spinning_reserve.data() + offset );
  }
  return( nullptr );
 }

/*--------------------------------------------------------------------------*/
 /// returns the array of secondary spinning reserve of the given \p generator
 /** This method returns the array of ColVariable representing the secondary
  * spinning reserve of the given \p generator at all time instants t in {0,
  * ..., time_horizon - 1}.
  *
  * @param generator The index of a generator (a number between 0 and
  *        get_number_generators() - 1).
  *
  * @return The array of ColVariable representing the secondary spinning
  *         reserve of the given \p generator. */

 ColVariable * get_secondary_spinning_reserve( Index generator ) override {
  if( generator < get_number_generators() ) {
   const auto offset = generator * f_time_horizon;
   if( offset < v_secondary_spinning_reserve.num_elements() )
    return( v_secondary_spinning_reserve.data() + offset );
  }
  return( nullptr );
 }

/*--------------------------------------------------------------------------*/
 /// returns the secondary spinning reserve of a generator at the given time
 /** Returns a pointer to the ColVariable representing the secondary spinning
  * reserve of the given \p generator at the given \p time.
  *
  * @param generator The index of the generator whose secondary spinning
  *        reserve is desired (a number between 0 and get_number_generators()
  *        - 1).
  *
  * @param time The time at which the secondary spinning reserve is desired (a
  *        number between 0 and get_time_horizon() - 1).
  *
  * @return A pointer to the ColVariable representing the secondary spinning
  *         reserve of the given \p generator at the given \p time. */

 ColVariable * get_secondary_spinning_reserve( Index generator , Index time ) {
  if( generator < get_number_generators() && time < f_time_horizon ) {
   const auto offset = generator * f_time_horizon + time;
   if( offset < v_secondary_spinning_reserve.num_elements() )
    return( v_secondary_spinning_reserve.data() + offset );
  }
  return( nullptr );
 }

/** @} ---------------------------------------------------------------------*/
/*------------------ METHODS FOR SAVING THE HydroUnitBlock------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for printing & saving the HydroUnitBlock
 * @{ */

 /// extends Block::serialize( netCDF::NcGroup )
 /** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
  * HydroUnitBlock. See HydroUnitBlock::deserialize( netCDF::NcGroup ) for
  * details of the format of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/** @} ---------------------------------------------------------------------*/
/*--------------- METHODS FOR INITIALIZING THE HydroUnitBlock --------------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the data of the HydroUnitBlock
 * @{ */

 void load( std::istream & input , char frmt = 0 ) override {
  throw( std::logic_error( "HydroUnitBlock::load() not implemented yet" ) );
 }

/** @} ---------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

 void set_inflow( MF_dbl_it values ,
                  Subset && subset ,
                  const bool ordered = false ,
                  ModParam issuePMod = eNoBlck ,
                  ModParam issueAMod = eNoBlck );

 void set_inflow( MF_dbl_it values ,
                  Range rng = Range( 0 , Inf< Index >() ) ,
                  ModParam issuePMod = eNoBlck ,
                  ModParam issueAMod = eNoBlck );

 void set_inertia_power( MF_dbl_it values ,
                         Subset && subset ,
                         const bool ordered = false ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck );

 void set_inertia_power( MF_dbl_it values ,
                         Range rng = Range( 0 , Inf< Index >() ) ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck );

 void set_initial_volume( MF_dbl_it values ,
                          Subset && subset ,
                          const bool ordered = false ,
                          ModParam issuePMod = eNoBlck ,
                          ModParam issueAMod = eNoBlck );

 void set_initial_volume( MF_dbl_it values ,
                          Range rng = Range( 0 , Inf< Index >() ) ,
                          ModParam issuePMod = eNoBlck ,
                          ModParam issueAMod = eNoBlck );

 void set_initial_flow_rate( MF_dbl_it values ,
                             Subset && subset ,
                             const bool ordered = false ,
                             ModParam issuePMod = eNoBlck ,
                             ModParam issueAMod = eNoBlck );

 void set_initial_flow_rate( MF_dbl_it values ,
                             Range rng = Range( 0 , Inf< Index >() ) ,
                             ModParam issuePMod = eNoBlck ,
                             ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED METHODS OF THE CLASS ---------------------*/
/*--------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

/*---------------------------------- data ----------------------------------*/

 /// the number of reservoirs (nodes) of the problem
 Index f_NumberReservoirs{};

 /// the number of connecting arcs which are connecting the reservoirs in
 /// cascading system
 Index f_NumberArcs{};

 /// the total number of pieces
 Index f_TotalNumberPieces{};

 /// the vector of UphillDelay
 std::vector< int > v_UphillDelay;

 /// the vector of DownhillDelay
 std::vector< Index > v_DownhillDelay;

 /// the vector of starting arc
 std::vector< Index > v_StartArc;

 /// the vector of ending arcs
 std::vector< Index > v_EndArc;

 /// the vector of initial volumetric
 std::vector< double > v_InitialVolumetric;

 /// the vector of initial flow rate
 std::vector< double > v_InitialFlowRate;

 /// the vector of NumberPieces
 std::vector< Index > v_NumberPieces;

 /// the vector of LinearTerm
 std::vector< double > v_LinearTerm;

 /// the vector of ConstTerm
 std::vector< double > v_ConstTerm;

 /// the matrix of inertia power of generators
 boost::multi_array< double , 2 > v_InertiaPower;

 /// the matrix of MinVolumetric
 /** Indexed over the dimensions NumberReservoirs and NumberIntervals. */
 boost::multi_array< double , 2 > v_MinVolumetric;

 /// the matrix of MaxVolumetric
 /** Indexed over the dimensions NumberReservoirs and NumberIntervals. */
 boost::multi_array< double , 2 > v_MaxVolumetric;

 /// the matrix of Inflows
 /** Indexed over the dimensions NumberReservoirs and NumberIntervals. */
 boost::multi_array< double , 2 > v_inflows;

 /// the matrix of MinPower
 /** Indexed over the dimensions NumberIntervals and NumberArcs. */
 boost::multi_array< double , 2 > v_MinPower;

 /// the matrix of MaxPower
 /** Indexed over the dimensions NumberIntervals and NumberArcs. */
 boost::multi_array< double , 2 > v_MaxPower;

 /// the matrix of MinFlow
 /** Indexed over the dimensions NumberIntervals and NumberArcs. */
 boost::multi_array< double , 2 > v_MinFlow;

 /// the matrix of MaxFlow
 /** Indexed over the dimensions NumberIntervals and NumberArcs. */
 boost::multi_array< double , 2 > v_MaxFlow;

 /// the matrix of DeltaRampUp
 /** Indexed over the dimensions NumberIntervals and NumberGenerators. */
 boost::multi_array< double , 2 > v_DeltaRampUp;

 /// the matrix of DeltaRampDown
 /** Indexed over the dimensions NumberIntervals and NumberGenerators. */
 boost::multi_array< double , 2 > v_DeltaRampDown;

 /// the matrix of PrimaryRho
 /** Indexed over the dimensions NumberIntervals and NumberArcs. */
 boost::multi_array< double , 2 > v_PrimaryRho;

 /// the matrix of SecondaryRho
 /** Indexed over the dimensions NumberIntervals and NumberArcs. */
 boost::multi_array< double , 2 > v_SecondaryRho;

/*-------------------------------- variables -------------------------------*/

 /// the matrix of volumetric variables
 boost::multi_array< ColVariable , 2 > v_volumetric;

 /// the matrix of flow rate variables
 boost::multi_array< ColVariable , 2 > v_flow_rate;

 /// the active power variables
 boost::multi_array< ColVariable , 2 > v_active_power;

 /// the primary spinning reserve variables
 boost::multi_array< ColVariable , 2 > v_primary_spinning_reserve;

 /// the secondary spinning reserve variables
 boost::multi_array< ColVariable , 2 > v_secondary_spinning_reserve;

/*------------------------------- constraints ------------------------------*/

 /// maximum power output according to primary-secondary reserves constraints
 boost::multi_array< FRowConstraint , 2 > MaxPowerPrimarySecondary_Const;

 /// minimum power output according to primary-secondary reserves constraints
 boost::multi_array< FRowConstraint , 2 > MinPowerPrimarySecondary_Const;

 /// power output relation with to primary reserves constraints
 boost::multi_array< FRowConstraint , 2 > ActivePowerPrimary_Const;

 /// power output relation with to secondary reserves constraints
 boost::multi_array< FRowConstraint , 2 > ActivePowerSecondary_Const;

 /// flow to active power function constraints
 boost::multi_array< FRowConstraint , 2 > FlowActivePower_Const;

 /// active power bounds
 boost::multi_array< FRowConstraint , 2 > ActivePowerBounds_Const;

 /// ramp-up constraints
 boost::multi_array< FRowConstraint , 2 > RampUp_Const;

 /// ramp-down constraints
 boost::multi_array< FRowConstraint , 2 > RampDown_Const;

 /// final volumes fo each reservoir constraints
 boost::multi_array< FRowConstraint , 2 > FinalVolumeReservoir_Const;


 /// flow rate bounds constraints
 boost::multi_array< BoxConstraint , 2 > FlowRateBounds_Const;

 /// volumetric bounds constraints
 boost::multi_array< BoxConstraint , 2 > VolumetricBounds_Const;


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

 /// transposes a deserialized multi-array if needed
 /** We deal with two-dimensional arrays that have dimensions ( time horizon per
  * number of arcs ). When provided by a netCDF variable, the size of the
  * dimension associated with the time horizon is allowed to be 1 (even if the
  * time horizon is greater than 1). This means that the given data does not
  * change over time. Therefore, the dimensions of a given array could be ( 1
  * x number of arcs ). This can be viewed as a "row vector". This being a
  * vector, the user may decide to provide a one-dimensional array whose size
  * is the number of arcs. This, however, is translated into a two-dimensional
  * array whose dimensions are ( number of arcs x 1 ), i.e., a "column
  * vector". In this case, the array must be transposed, so that its second
  * dimension becomes the number of arcs (and therefore compatible with our
  * data structure).
  *
  * @tparam T The type of the boost::multi_array.
  *
  * @param a A boost::multi_array that has been just deserialized.
  */
 template< typename T >
 void transpose( boost::multi_array< T , 2 > & a );

 /// decompress a multi_array using the change intervals
 void decompress_array( boost::multi_array< double , 2 > & a );

 /// decompress a max/min volumetric multi_array using the change intervals
 void decompress_vol( boost::multi_array< double , 2 > & a );

/*--------------------------------------------------------------------------*/
 /// updates the constraints for the given arcs at time 0
 /** This function updates the right-hand side of the ramp-up constraints and
  * the left-hand side of the ramp-down constraints associated with the given
  * \p arcs at time 0 (which are the ramp constraints that depend on the
  * initial flow rate).
  *
  * @param arcs The indices of the arcs whose constraints must be updated.
  */

 void update_initial_flow_rate_in_cnstrs( Range arcs ,
                                          c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// updates the constraints for the given arcs at time 0
 /** This function updates the right-hand side of the ramp-up constraints and
  * the left-hand side of the ramp-down constraints associated with the given
  * \p arcs at time 0 (which are the ramp constraints that depend on the
  * initial flow rate).
  *
  * @param arcs The indices of the arcs whose associated constraints must
  *        be updated.
  */

 void update_initial_flow_rate_in_cnstrs( const Block::Subset & arcs ,
                                          c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/

 static void static_initialization( void ) {

  /* Warning: Not all C++ compilers enjoy the template wizardry behind the
   * three-args version of register_method<> with the compact MS_*_*::args()
   *
   * register_method< HydroUnitBlock >( "HydroUnitBlock::set_inflow",
   *                                    &HydroUnitBlock::set_inflow,
   *                                    MS_dbl_sbst::args() );
   *
   * so we just use the slightly less compact one with the explicit argument
   * and be done with it. */

  register_method< HydroUnitBlock , MF_dbl_it , Subset && , bool >(
   "HydroUnitBlock::set_inflow" , &HydroUnitBlock::set_inflow );

  register_method< HydroUnitBlock , MF_dbl_it , Range >(
   "HydroUnitBlock::set_inflow" , &HydroUnitBlock::set_inflow );

  register_method< HydroUnitBlock , MF_dbl_it , Subset && , bool >(
   "HydroUnitBlock::set_inertia_power" , &HydroUnitBlock::set_inertia_power );

  register_method< HydroUnitBlock , MF_dbl_it , Range >(
   "HydroUnitBlock::set_inertia_power" , &HydroUnitBlock::set_inertia_power );

  register_method< HydroUnitBlock , MF_dbl_it , Subset && , bool >(
   "HydroUnitBlock::set_initial_volume" , &HydroUnitBlock::set_initial_volume );

  register_method< HydroUnitBlock , MF_dbl_it , Range >(
   "HydroUnitBlock::set_initial_volume" , &HydroUnitBlock::set_initial_volume );
 }

};  // end( class( HydroUnitBlock ) )

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS HydroUnitBlockMod ------------------------*/
/*--------------------------------------------------------------------------*/

/// derived class from Modification for modifications to a HydroUnitBlock
class HydroUnitBlockMod : public UnitBlockMod
{

 public:

 /// public enum for the types of HydroUnitBlockMod
 enum HUB_mod_type
 {
  eSetInf = eUBModLastParam , ///< set inflow values
  eSetInerP ,                 ///< set inertia power values
  eSetInitF ,                 ///< set initial flow rate values
  eSetInitV                   ///< set initial volumetric values
 };

 /// constructor, takes the HydroUnitBlock and the type
 HydroUnitBlockMod( HydroUnitBlock * const fblock ,
                    const int type ) :
  UnitBlockMod( fblock , type ) {}

 /// destructor, does nothing
 virtual ~HydroUnitBlockMod( void ) override = default;

 /// returns the Block to which the Modification refers
 Block * get_Block( void ) const override { return( f_Block ); }

 protected:

 /// prints the HydroUnitBlockMod
 void print( std::ostream & output ) const override {
  output << "HydroUnitBlockMod[" << this << "]: ";
  switch( f_type ) {
   case( eSetInf ):
    output << "set inflow values ";
    break;
   case( eSetInerP ):
    output << "set inertia power values ";
    break;
   case( eSetInitF ):
    output << "set initial flow rate values ";
    break;
   default:
    output << "set initial volumetric values ";
  }
 }

 HydroUnitBlock * f_Block{};
 ///< pointer to the Block to which the Modification refers

};  // end( class( HydroUnitBlockMod ) )

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS HydroUnitBlockRngdMod -----------------------*/
/*--------------------------------------------------------------------------*/

/// derived from HydroUnitBlockMod for "ranged" modifications
class HydroUnitBlockRngdMod : public HydroUnitBlockMod
{

 public:

 /// constructor: takes the HydroUnitBlock, the type, and the range
 HydroUnitBlockRngdMod( HydroUnitBlock * const fblock ,
                        const int type ,
                        Block::Range rng )
  : HydroUnitBlockMod( fblock , type ) , f_rng( rng ) {}

 /// destructor, does nothing
 virtual ~HydroUnitBlockRngdMod() override = default;

 /// accessor to the range
 Block::c_Range & rng( void ) { return( f_rng ); }

 protected:

 /// prints the HydroUnitBlockRngdMod
 void print( std::ostream & output ) const override {
  HydroUnitBlockMod::print( output );
  output << "[ " << f_rng.first << ", " << f_rng.second << " )" << std::endl;
 }

 Block::Range f_rng;  ///< the range

};  // end( class( HydroUnitBlockRngdMod ) )

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS HydroUnitBlockSbstMod ---------------------*/
/*--------------------------------------------------------------------------*/

/// derived from HydroUnitBlockMod for "subset" modifications
class HydroUnitBlockSbstMod : public HydroUnitBlockMod
{

 public:

 /// constructor: takes the HydroUnitBlock, the type, and the subset
 HydroUnitBlockSbstMod( HydroUnitBlock * const fblock ,
                        const int type ,
                        Block::Subset && nms )
  : HydroUnitBlockMod( fblock , type ) , f_nms( std::move( nms ) ) {}

 /// destructor, does nothing
 virtual ~HydroUnitBlockSbstMod() override = default;

 /// accessor to the subset
 Block::c_Subset & nms( void ) { return( f_nms ); }

 protected:

 /// prints the HydroUnitBlockSbstMod
 void print( std::ostream & output ) const override {
  HydroUnitBlockMod::print( output );
  output << "(# " << f_nms.size() << ")" << std::endl;
 }

 Block::Subset f_nms;  ///< the subset

};  // end( class( HydroUnitBlockSbstMod ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* HydroUnitBlock.h included */

/*--------------------------------------------------------------------------*/
/*---------------------- End File HydroUnitBlock.h -------------------------*/
/*--------------------------------------------------------------------------*/
