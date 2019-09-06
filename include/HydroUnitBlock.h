/*--------------------------------------------------------------------------*/
/*------------------------- File HydroUnitBlock.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class HydroUnitBlock, which derives from UnitBlock
 * [see UnitBlock.h], in order to define a "reasonably standard" hydro unit
 * of a Unit Commitment Problem.
 *
 * \version 0.11
 *
 * \date 02 - 09 - 2019
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

#ifndef __HydroUnitBlock
#define __HydroUnitBlock
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "ColVariable.h"
#include "FRowConstraint.h"
#include "OneVarConstraint.h"
#include "FRealObjective.h"
#include "DQuadFunction.h"
#include "UnitBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------------ NAMESPACE ---------------------------------*/
/*--------------------------------------------------------------------------*/

/// Namespace for the Structured Modeling System++ (SMS++)

namespace SMSpp_di_unipi_it {

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

class HydroUnitBlock : public UnitBlock {

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
/** @name Constructor and Destructor
 *  @{ */

/// constructor, takes the father and the time horizon
/** Constructor of HydroUnitBlock, taking possibly a pointer of its
 * father Block.
 */

 explicit HydroUnitBlock( Block * f_block = nullptr , Index t = 0):
         UnitBlock( f_block ) {}

/*--------------------------------------------------------------------------*/

/// destructor of HydroUnitBlock

 ~HydroUnitBlock() override = default;

/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

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
 * - The variable "StartArc", of type int and indexed over the dimension
 *   "NumberArcs"; the r-th entry of the variable is the starting point of
 *   the arc (a number in 0, ..., NumberReservoirs - 1). Note that arcs are
 *   oriented; that is, a positive flow along arc r (turbine) means that water
 *   is being taken away from StartArc[ r ] and delivered to EndArc[ r ]
 *   (see next), a negative flow (pump) means vice-versa. Note that reservoir
 *   names here go from 0 to NumberReservoirs - 1.
 *
 * - The variable "EndArc", of type int and indexed over the dimension
 *   "NumberArcs"; the r-th entry of the variable is the ending point of the
 *   arc; this is a number in 0, ..., NumberReservoirs. Note: this is
 *   NumberReservoirs and *not* NumberReservoirs - 1, because arcs can end in
 *   the "fake" reservoir NumberReservoirs. This indicates that water that
 *   flows along that arc "goes away from the system" and it is no longer
 *   counted, because it can no longer be used to produce electricity further
 *   down the river, or pumped back into one of its reservoirs. Indeed, there
 *   will be something like "the most downstream turbines": after water has
 *   been used there, it just goes away down some river and does not go to
 *   any other reservoir. Arcs are oriented (see above); StartArc[ r ] ==
 *   EndArc[ r ] (a self-loop) is not allowed, but multiple arcs between the
 *   same pair of reservoirs are. Indeed, often the same physical equipment
 *   can be used both as a turbine and as a pump; in our model these are
 *   represented as two parallel arcs (but with different upper and lower
 *   flow capacity, see "MinFlow" and "MaxFlow" below).
 *
 * - The variable "MinFlow", of type double and indexed over both dimensions
 *   "NumberIntervals" and "NumberArcs". The first dimension may have either
 *   size 1 or size NumberIntervals whereas the second one always has size
 *   NumberArcs (if the variable is provided at all). This is meant to
 *   represent the matrix MinF[ t , l ] which, for each time instant t and
 *   each arc l, contains the minimum flow value of the unit. This variable
 *   is optional; if it is not provided then it is assumed that MinF[ t , l ]
 *   == 0, i.e., the minimum flow of the unit is zero. If the first dimension
 *   has size 1 then the entry MinF[ 0 , l ] gives the fixed minimum flow
 *   value of the unit for all time steps and each arc l. Otherwise,
 *   MinFlow[ i , l ] is the fixed value of MinF[ t , l ] for all time t and
 *   arc l in the interval [ ChangeIntervals[ i - 1 ] ,
 *   ChangeIntervals[ i ] ], with the assumption that ChangeIntervals[ - 1 ]
 *   = 0. If NumberIntervals <= 1 or NumberIntervals >= TimeHorizon, then the
 *   mapping clearly does not require "ChangeIntervals", which in fact is not
 *   loaded.
 *
 * - The variable "MaxFlow", of type double and indexed over both dimensions
 *   "NumberIntervals" and "NumberArcs". The first dimension may have either
 *   size 1 or size NumberIntervals whereas the second one always has size
 *   NumberArcs (if the variable is provided at all). This is meant to
 *   represent the matrix MaxF[ t , l ] which, for each time instant t and
 *   each arc l, contains the maximum flow value of the unit. This variable
 *   is optional; if it is not provided then it is assumed that MaxF[ t , l ]
 *   == 0, i.e., the maximum flow of the unit is zero. If the first dimension
 *   has size 1 then the entry MaxF[ 0 , l ] gives the fixed maximum flow
 *   value of the unit for all time steps and each arc l. Otherwise,
 *   MaxFlow[ i , l ] is the fixed value of MaxF[ t , l ] for all time t and
 *   arc l in the interval [ ChangeIntervals[ i - 1 ] ,
 *   ChangeIntervals[ i ] ], with the assumption that ChangeIntervals[ - 1 ]
 *   = 0. If NumberIntervals <= 1 or NumberIntervals >= TimeHorizon, then the
 *   mapping clearly does not require "ChangeIntervals", which in fact is not
 *   loaded.
 *
 * Note: MinFlow and MaxFlow values can be either positive or negative (or
 * zero); whenever MinF[ t , l ] < MaxF[ t , l ] <= 0 for each t and l,
 * the unit is considered a pump and whenever 0 <= MinF[ t , l ] <
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
 * - The variable "MinVolumetric", of type double and indexed over both
 *   dimensions "NumberReservoirs" and "NumberIntervals". The first dimension
 *   always has size NumberReservoirs (if it is provided at all), whereas the
 *   second one may have size one or size NumberIntervals. This is meant to
 *   represent the matrix MinV[ r , t ] which, for each reservoir r at each
 *   time instant t contains the minimum volumetric value of the unit for 
 *   each reservoir and corresponding time step. It must be that
 *   0 <= MinV[ r , t ] < MaxV[ r , t ] for all r and t and . This variable
 *   is optional; if it's not provided then it is assumed that
 *   MinV[ r , t ] == 0, i.e., the minimum volumetric of the unit is zero. If
 *   the second dimension has size 1 then the entry MinV[ r , 0 ] gives the
 *   fixed minimum volumetric value of the unit for each reservoir r along
 *   all the time horizon. Otherwise, MinVolumetric[ r , i ] is the fixed
 *   value of MinV[ r , t ] for reservoir r and all t in the interval
 *   [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the assumption
 *   that ChangeIntervals[ - 1 ] = 0. If NumberIntervals <= 1 or
 *   NumberIntervals >= TimeHorizon, then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded.
 *
 * - The variable "MaxVolumetric", of type double and indexed over both
 *   dimensions "NumberReservoirs" and "NumberIntervals". The first dimension
 *   is always has size NumberReservoirs (if it is provided at all), whereas
 *   the second one may have size one or size NumberIntervals. This is meant
 *   to represent the matrix MaxV[ r , t ] which, for each reservoir r at
 *   each time instant t contains the maximum volumetric value of the unit
 *   for each reservoir and corresponding time step. It must be that
 *   0 <= MinV[ r , t ] < MaxV[ r , t ] for all r and t. This variable is not
 *   optional (a reservoir must have some available volume). If the second
 *   dimension has size 1 then the entry MaxV[ r , 0 ] gives the fixed
 *   maximum volumetric value of the unit for each reservoir r during the all
 *   the time horizon. Otherwise, the MaxVolumetric[ r , i ] is the fixed
 *   value of MaxV[ r , t ] for reservoir r and all t in the interval
 *   [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the assumption
 *   that ChangeIntervals[ - 1 ] = 0. If NumberIntervals <= 1 or
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
 * - The variable "Inflows", of type double and indexed over both dimensions
 *   "NumberReservoirs" and "TimeHorizon". This is meant to represent the
 *   matrix InF[ r , t ] which, for each reservoir r, contains the amount of
 *   water that "naturally" enters into reservoir r (because of rain, ice
 *   melting, non-controlled rivers flowing, and of course net of water
 *   leaving by evaporation, human consumption etc.) during time al the time
 *   interval t, and therefore that is available in the reservoir at the end
 *   of time step t (hence, the beginning of time step t + 1, if any). This
 *   variable is optional; if it isn't defined, it is taken to be zero.
 *   Inflows can be either positive or negative.
 *
 * - The variable "MinPower", of type double and indexed over both dimensions
 *   "NumberIntervals" and "NumberArcs". The first dimension may have either
 *   size 1 or size NumberIntervals, whereas the second one always has size
 *   NumberArcs (if it is provided at all). This is meant to represent the
 *   matrix MinP[ t , l ] which, for each time instant t at each arc l
 *   contains the minimum power value of the unit; it must be that
 *   MinP[ t , l ] < MaxP[ t , l ] for each time instant t and each arc l.
 *   This variable is optional; if it is not provided then it is assumed
 *   that MinP[ t , l ] == 0, i.e., the minimum power of all units is zero
 *   (which means, each unit is a turbine). If the first dimension has size 1
 *   then the entry MinP[ 0 , l ] is assumed to contain the minimum power of
 *   arc l for all time instants. Otherwise, MinPower[ i , l ] is the fixed
 *   value of MinP[ t , l ] for arc l and all time t in the interval
 *   [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the
 *   assumption that ChangeIntervals[ - 1 ] = 0. If NumberIntervals <= 1 or
 *   NumberIntervals >= TimeHorizon, then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded.
 *
 * - The variable "MaxPower", of type double and indexed over both dimensions
 *   "NumberIntervals" and "NumberArcs". The first dimension may have either
 *   size 1 or size NumberIntervals whereas the second one always has size
 *   NumberArcs (if it is provided at all). This is meant to represent the
 *   matrix MaxP[ t , l ] which, for each time instant t at each arc l
 *   contains the maximum power value of the unit; it must be that
 *   MinP[ t , l ] < MaxP[ t , l ] for each time instant t and each arc l.
 *   This variable is optional; if it is not provided then it is assumed that
 *   MaxP[ t , l ] == 0, i.e., the maximum power of all units is zero (i.e.,
 *   all units are pumps). If the first dimension has size 1 then the entry
 *   MaxP[ 0 , l ] is assumed to contain the maximum power of arc l for all
 *   time instants. Otherwise, MaxPower[ i , l ] is the fixed value of
 *   MaxP[ t , l ] for arc l and all time t in the interval
 *   [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the assumption
 *   that ChangeIntervals[ - 1 ] = 0. If NumberIntervals <= 1 or
 *   NumberIntervals >= TimeHorizon, then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded.
 *
 * - The variable "DeltaRampUp", of type double and indexed over both
 *   dimensions "NumberIntervals" and "NumberArcs". The first dimension may
 *   have either size 1 or size NumberIntervals, whereas the second one
 *   always has size NumberArcs (if it is provided at all). This is meant
 *   to represent the matrix DP[ t , l ] which contains the maximum possible
 *   increase of the flow rate at time instant t for arc l. This variable is
 *   optional; if it is not provided then it is assumed that DP[ t , l ] ==
 *   MaxP[ t , l ] - MinP[ t , l ], i.e., all units can ramp up by an
 *   arbitrary amount, i.e., there are no ramp-up constraints. If the first
 *   dimension has size 1 then the entry DP[ 0 , l ] is assumed to contain
 *   the maximum possible increase of the flow rate of arc l for all time
 *   instants. Otherwise, DeltaRampUp[ i , l ] is the fixed value of
 *   DP[ t , l ] for arc l and all time t in the interval
 *   [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the assumption
 *   that ChangeIntervals[ - 1 ] = 0. If NumberIntervals <= 1 or
 *   NumberIntervals >= TimeHorizon, then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded.
 *
 * - The variable "DeltaRampDown", of type double and indexed over both
 *   dimensions "NumberIntervals" and "NumberArcs". The first dimension may
 *   have either size 1 or size NumberIntervals, whereas the second one
 *   always has size NumberArcs (if it is provided at all). This is meant
 *   to represent the matrix DM[ t , l ] which contains the maximum possible
 *   decrease of the flow rate at each time instant t of each arc l. This
 *   variable is optional; if it is not provided then it is assumed that
 *   DM[ t , l ] == MaxP[ t , l ] - MinP[ t , l ], i.e., the unit can ramp
 *   down by an arbitrary amount, i.e., there are no ramp-down constraints.
 *   If first dimension has size 1 then the entry DM[ 0 , l ] is assumed to
 *   contain the maximum possible decrease of the flow rate of arc l for all
 *   time instants. Otherwise, DeltaRampDown[ i , l ] is the fixed value of
 *   DM[ t , l ] for arc l and all time t  in the interval
 *   [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the assumption
 *   that ChangeIntervals[ - 1 ] = 0. If NumberIntervals <= 1 or
 *   NumberIntervals >= TimeHorizon, then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded.
 *
 * - The variable "PrimaryRho", of type double and indexed both over the
 *   dimensions "NumberIntervals" and "NumberArcs". The first dimension may
 *   have either size 1 or size NumberIntervals, whereas the second one
 *   always has size NumberArcs (if it is provided at all). This is meant
 *   to represent the matrix PR[ t , l ] which, for each time instant t and
 *   arc l, contains the maximum possible fraction of active power that can
 *   be used as primary reserve. This variable is optional, when it's not
 *   present then PR[ t , l ] == 0 for all t and l, i.e., the unit is not
 *   capable of producing any primary reserve. Note that only turbines can
 *   produce primary reserve, i.e., PR[ t , l ] > 0 ==> MaxP[ t , l ] > 0.
 *   If the first dimension has size 1 then the entry PR[ 0 , l ] is assumed
 *   to contain the maximum possible fraction of active power that can be
 *   used as primary reserve by arc l for all time instants. Otherwise,
 *   PrimaryRho[ i , l ] is the fixed value of PR[ t , l ] for arc l and all
 *   t in the interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ],
 *   with the assumption that ChangeIntervals[ - 1 ] = 0 and all l. If
 *   NumberIntervals <= 1 or NumberIntervals >= TimeHorizon, then the mapping
 *   clearly does not require "ChangeIntervals", which in fact is not loaded.
 *
 * - The variable "SecondaryRho", of type double and indexed both over the
 *   dimensions "NumberIntervals" and "NumberArcs". The first dimension may
 *   have either size 1 or size NumberIntervals, whereas the second one
 *   always has size NumberArcs (if it is provided at all). This is meant to
 *   represent the matrix SR[ t , l ] which, for each time instant t and arc
 *   l contains the maximum possible fraction of active power that can be
 *   used as secondary reserve. This variable is optional, when it's not
 *   present then SR[ t , l ] == 0 for all t and l, i.e., the unit is not
 *   capable of producing any secondary reserve. Note that only turbines can
 *   produce secondary reserve, i.e., SR[ t , l ] > 0 ==> MaxP[ t , l ] > 0.
 *   If the first dimension has size 1 then the entry SR[ 0 , l ] is assumed
 *   to contain the maximum possible fraction of active power that can use as
 *   secondary reserve by arc l for all time instant. Otherwise,
 *   SecondaryRho[ i , l ] is the fixed value of SR[ t , l ] for arc l and
 *   all t in the interval [ ChangeIntervals[ i - 1 ] ,
 *   ChangeIntervals[ i ] ], with the assumption that ChangeIntervals[ - 1 ]
 *   = 0. If NumberIntervals <= 1 or NumberIntervals >= TimeHorizon, then the
 *   mapping clearly does not require "ChangeIntervals" which in fact is not
 *   loaded.
 *
 * - The variable "NumberPieces", indexed over the dimension "NumberArcs".
 *   NumberPieces[ l ] tells how many pieces the concave flow-to-active-power
 *   function has for unit (arc) l. Note that pumps must necessarily have
 *   exactly one piece. The sum over all i of NumberPieces[ i ] is the total
 *   number of pieces (say, "TotalNumberPieces"). Clearly, TotalNumberPieces
 *   >= NumberArcs; if the flow-to-active-power function for all turbines only
 *   have one piece (those of pumps necessarily are so), then 
 *   TotalNumberPieces == NumberArcs and there is no need to define this
 *   variable. If, instead, if it is defined, then should always be such that
 *   TotalNumberPieces >= NumberArcs.
 *
 * - The variable "LinearTerm", of type double and indexed over the set
 *   { 0 , ..., TotalNumberPieces - 1 } (see "NumberPieces"). LinearTerm[ h ]
 *   gives the linear term a_h of the linear function a_h * f + b_h that
 *   defines the concave flow-to-active-power function for some unit; the
 *   total function if F2AP( t ) = min { a_h * f + b_h , h \in H } for some
 *   finite set H that depends on the individual unit. It is then necessary
 *   to be able to assign a unique index h = 0, 1, ..., TotalNumberPieces - 1
 *   to each pair ( unit , linear function a_h * f + b_h ). When
 *   TotalNumberPieces == NumberArcs, the index is the same as i = 0, 1, ...,
 *   NumberArcs - 1 (there is a one-to-one correspondence between each (unit)
 *   arc and each piece). When, instead, TotalNumberPieces > NumberArcs, a
 *   mapping must be defined. The mapping is the obvious one: each index of
 *   i = 0, 1, ..., NumberArcs - 1, corresponds with a unit (arc), and the
 *   linear functions for each unit (arc) also have some natural ordering.
 *   Thus, in general the mapping is:
 *     piece 0 = first piece of unit (arc) 0
 *     piece 1 = second piece of unit (arc) 0
 *     ...
 *     piece NumberPieces[ 0 ] - 1 = last piece of unit (arc) 0
 *     piece NumberPieces[ 0 ] = first piece of unit (arc) 1
 *     piece NumberPieces[ 0 ] + 1 = second piece of unit (arc) 1
 *     ...
 *   which of course boils down to "h = i" when each arc has exactly one
 *   piece.
 *
 * - The variable "ConstantTerm", of type double and indexed over the set
 *   { 0 , ..., TotalNumberPieces" - 1 }. ConstantTerm[ h ] gives the
 *   constant term b_h of the linear function a_h * f + b_h that defines the
 *   concave flow-to-active-power function for some unit; see the comments to
 *   "LinearTerm" for details.
 *
 * - The variable "InertiaPower", of type double and indexed both over the
 *   dimensions "NumberIntervals" and "NumberArcs". The first dimension may
 *   have either size 1 or size NumberIntervals whereas the second one
 *   always has size NumberArcs (if it is provided at all). This is meant to
 *   represent the matrix IP[ t , l ] which, for each time instant t and arc
 *   l, contains the contribution that the unit can give to the inertia
 *   constraint which depends on the active power that it is currently
 *   generating (basically, the constant to be multiplied to the active
 *   power variable) at time t for arc l. The variable is optional; if it is
 *   not defined, IP[ t , l ] == 0 for each time instants t and arc l. If the
 *   first dimension has size 1 then the entry IP[ 0 , l ] is assumed to
 *   contain the the inertia power value for arc l and all time instants t.
 *   Otherwise, InertiaPower[ i , l ] is the fixed value of IP[ t , l ] for
 *   all t in the interval [ ChangeIntervals[ i - 1 ] ,
 *   ChangeIntervals[ i ] ], with the assumption that ChangeIntervals[ - 1 ]
 *    = 0 and all l. If  NumberIntervals <= 1 or NumberIntervals >= 
 *   TimeHorizon then the mapping clearly does not require "ChangeIntervals",
 *   which in fact is not loaded.
 *
 * - The variable "InitialFlowRate", of type double and indexed over the
 *   dimension "NumberArcs". Each entry InFR[ i ] indicates the amount of the
 *   flow that was going along arc i at time instant -1. This is necessary to
 *   compute ramp-up and ramp-down limits (cf. "DeltaRampUp" and
 *   "DeltaRampDown"), and therefore it is useless if there are no ramp
 *   constraints on *any* unit (arc), in which case it is not loaded.
 *
 * - The variable "InitialVolumetric", of type double and indexed over the
 *   dimension "NumberReservoirs". Each entry InV[ r ] indicates the volumes
 *   of water in reservoir r at time instant -1.
 *
 * - The negative or positive scalar variable "UphillFlow", of type Int64 and
 *   indexed over the dimension "NumberArcs". Each entry UpF[ l ] indicates
 *   the uphill flow delay for each unit (arc) l; the nontrivial concept is
 *   detailed below. This variable is optional, if it is not provided it is
 *   taken to be UpF[ l ] == 0.
 *
 * - The positive scalar variable "DownhillFlow", of type UInt64 and indexed
 *   over the dimension "NumberArcs". Each entry DnF[ l ] indicates the
 *   downhill flow delay for each arc (unit) l. This variable is optional, if
 *   it is not provided it is taken to be DnF[ l ] == 0.
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

 void deserialize( netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
/// generate the abstract variables of the HydroUnitBlock
/** The HydroUnitBlock class use get_variable() method to access to each
 *  "group" of variable that may create in UnitBlock class which are:
 *
 *  - the primary spinning reserve variables;
 *
 *  - the secondary spinning reserve variables;
 *
 *  - the active power variables;
 *
 *  All of those variables are optional except the active power variables in
 *  the sense that the model may just not have them and whenever a group of
 *  above variables is created, its size will be the time horizon. Moreover,
 *  HydroUnitBlock is defined more groups of variables as follow:
 *
 *  - the volumetric variables
 *
 *  - the flow rate variables
 *
 *  These two groups of variables may have size f_time_horizon or empty size.
 *  All of these variables are optional,and it is also possible to restrict
 *  which of the subsets are generated with the parameter stvv. If stvv is not
 *  nullptr and it is a SimpleConfiguration<int>, or if
 *  f_BlockConfig->f_static_variables_Configuration is not nullptr and it is a
 *  SimpleConfiguration<int>, then the f_value (an int) indicates whether each
 *  of the optional variables should be created. If the Configuration is not
 *  available, the default value is taken to be 0.
 * */
 void generate_abstract_variables( Configuration *stvv ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// generate the static constraint of the HydroUnit
/** This method generates the static constraint of the HydroUnitBlock.
 * In order to describe a hydro generating unit system, it will be convenient
 * to see a cascading system as a graph. Let \f$ \mathcal{N}^{hy}\f$ be the
 * set of reservoirs (nodes) and \f$ \mathcal{L}^{hy}\f$ be the set of arcs
 * connecting these reservoirs respectively. Attached to each
 * \f$ l \in \mathcal{L}^{hy}\f$ are one or several plants(turbines or pumps).
 * This system is described on a discrete time horizon as dictated by the
 * UnitBlock interface. In this description we indicate it with
 * \f$ \mathcal{T}=\{ 0, \dots , \mathcal{|T|} - 1\} \f$. Each reservoir
 * \f$ n \in \mathcal{N}^{hy}\f$ has a continuous volumetric variables
 * \f$ v^{hy}_{n,t}\f$ in \f$ m^3 \f$ for \f$ t \in \mathcal{T}\f$ with
 * associated lower and upper bounds \f$ V^{hy,mn}_{n,t}\f$,
 * \f$ V^{hy,mx}_{n,t}\f$ and inflows \f$ A_{n,t}\f$ in \f$ m^3 /s \f$. The
 * uphill and downhill flow rate are defined as \f$ \tau^{up} \f$ and
 * \f$ \tau^{dn}\f$ respectively. For each arc \f$ l \in \mathcal{L}^{hy}\f$
 * in each time \f$ t \in \mathcal{T}\f$ the continuous flow rate variable
 * \f$ f_{l,t} \f$ in \f$ m^3 /s \f$ and ramping conditions
 * \f$ \Delta^{up}_{l,t} \f$ and \f$ \Delta^{dn}_{l,t} \f$ in
 * \f$ (m^3 /s)/h \f$ are disposed. The flow rate variable will be subject to
 * bounds \f$ F^{mn}_{l,t} \f$ and \f$ F^{mx}_{l,t} \f$ and it's assumed
 * moreover given a cutting plane model describing power as a function of flow
 * rate as below:
 *   \f[
 *      p^{ac}_{t,l}(f) := min \{ P_l + \rho^{hy}_{l}f_{t,l}\}
 *   \f]
 * where \f$ P_l\f$ and \f$ \rho^{hy}_{l} \f$ are considered as the constant
 * and linear multipliers of the linear function(flow-to-active-power)
 * respectively. Power generated by the hydro unit in each time and for each
 * arc \f$ p^{ac}_{t,l}, p^{pr}_{t,l}, p^{sc}_{t,l} \f$ in MW will be subject
 * to bounds \f$ P^{mn}_{t,l} \f$ and \f$ P^{mx}_{t,l} \f$ respectively.
 * Besides, we emphasize that reserve requirements are specified in order to
 * be symmetrically available to increase or decrease power injected into the
 * grid. For some of the constraints we will need to distinguish between pumps
 * and turbines. The distinction is made by considering the set of feasible
 * flow rates. Whenever \f$ [ F^{mn}_{l,t} , F^{mx}_{l,t}] \subseteq R_- \f$
 * for each arc and each time, the unit is considered a pump, and whenever
 * \f$ [ F^{mn}_{l,t} , F^{mx}_{l,t}] \subseteq R_+ \f$ the unit is considered
 * a turbine. Any possible mixed situation can be accounted for by
 * artificially splitting the unit into “two units”, which should be done at
 * the data processing stage (see deserialize() comments). With above
 * description the mathematical constraint of hydro unit may present as below:
 *
 * - maximum and minimum power output constraints according to primary and
 *   secondary spinning reserves are are presented in (1)-(2). Each of them
 *   is a boost::multi_array<FRowConstraint, 2>; with two dimensions which are
 *   f_time_horizon, and f_number_arcs entries, where the entry
 *   t = 0, ...,f_time_horizon - 1 and the entry z = 0, ...,f_number_arcs - 1
 *   being the maximum and minimum power output value according to the primary
 *   and the secondary spinning reserves at time t and arc l. these ensure the
 *   maximum(or minimum) amount of energy that unit can produce(or use) when
 *   it is on(or off).
 *
 *   \f[
 *
 *      p^{ac}_{t,l} + p^{pr}_{t,l} + p^{sc}_{t,l} \leq P^{mx}_{t,l}
 *          \quad t \in \mathcal{T}, l \in \mathcal{L}^{hy}          \quad (1)
 *
 *   \f]
 *
 *   \f[
 *
 *     P^{mn}_{t,l} \leq p^{ac}_{t,l} - p^{pr}_{t,l} - p^{sc}_{t,l}
 *         \quad t \in \mathcal{T}, l \in \mathcal{L}^{hy}           \quad (2)
 *
 *   \f]
 *
 * - primary and secondary spinning reserves relation with active power at
 *   each time and for each turbine: the same as inequalities(1)-(2), the
 *   inequalities(3)-(4) ensure that maximum amount of primary and secondary
 *   spinning reserve in the problem. Each of them is a
 *   boost::multi_array<FRowConstraint, 2>; with two dimensions which are
 *   f_time_horizon, and f_number_arcs entries, where
 *   t = 0, ...,f_time_horizon - 1 and z = 0, ...,f_number_arcs - 1
 *
 *   \f[
 *
 *      p^{pr}_{t,l} \leq \rho^{pr}_{t,l}p^{ac}_{t,l} \quad t \in \mathcal{T},
 *        l \in \mathcal{L}^{hy} \quad with
 *        \quad  [ F^{mn}_{l,t} , F^{mx}_{l,t}] \subseteq R_+       \quad (3)
 *
 *   \f]
 *
 *   \f[
 *
 *      p^{sc}_{t,l} \leq \rho^{sc}_{t,l}p^{ac}_{t,l} \quad t \in \mathcal{T},
 *        l \in \mathcal{L}^{hy} \quad with
 *        \quad  [ F^{mn}_{l,t} , F^{mx}_{l,t}] \subseteq R_+       \quad (4)
 *
 *   \f]
 *  where \f$ \rho^{pr}_{t,l} \f$ and \f$ \rho^{sc}_{t,l}\f$ are the maximum
 *  possible fraction of active power at each time and each arc that can be
 *  used as primary and secondary reserve respectively.
 *
 * - primary and secondary spinning reserves at each time and for each pump:
 *   these equalities(5)-(6) ensure that the primary and secondary spinning
 *   reserve value for each pump is equal to zero. Each of them is a
 *   boost::multi_array<FRowConstraint, 2>; with two dimensions which are
 *   f_time_horizon, and f_number_arcs entries, where
 *   t = 0, ...,f_time_horizon - 1 and z = 0, ...,f_number_arcs - 1
 *
 *   \f[
 *
 *      p^{pr}_{t,l} = 0 \quad t \in \mathcal{T},
 *        l \in \mathcal{L}^{hy} \quad with
 *        \quad  [ F^{mn}_{l,t} , F^{mx}_{l,t}] \subseteq R_-      \quad (5)
 *
 *   \f]
 *
 *   \f[
 *
 *      p^{sc}_{t,l} = 0 \quad t \in \mathcal{T},
 *        l \in \mathcal{L}^{hy} \quad with
 *        \quad  [ F^{mn}_{l,t} , F^{mx}_{l,t}] \subseteq R_-     \quad (6)
 *
 *   \f]
 *
 * - flow to active power function at each time and for each pump: this
 *   equality(7) gives the active power relation with flow rate for each
 *   pump at time t. This is a boost::multi_array<FRowConstraint, 2>; with two
 *   dimensions which are f_time_horizon, and f_number_arcs entries, where
 *   t = 0, ...,f_time_horizon - 1 and z = 0, ...,f_number_arcs - 1
 *
 *   \f[
 *
 *      p^{ac}_{t,l} = \rho^{hy}_{l}f_{t,l} \quad t \in \mathcal{T},
 *        l \in \mathcal{L}^{hy} \quad with
 *        \quad  [ F^{mn}_{l,t} , F^{mx}_{l,t}] \subseteq R_-       \quad (7)
 *
 *   \f]
 *
 * - flow-to-active-power function at each time and for each turbine ;
 *
 *   \f[
 *
 *      p^{ac}_{t,l} \leq P_j + \rho^{hy}_{j}f_{t,l} \quad j \in \mathcal{J}_l
 *        \quad t \in \mathcal{T}, l \in \mathcal{L}^{hy} \quad with
 *        \quad  [ F^{mn}_{l,t} , F^{mx}_{l,t}] \subseteq R_+       \quad (8)
 *
 *   \f]
 *
 * - flow rate variable bounds: This inequality(9) indicates upper and lower
 *   bound of flow rate at time t and for ach arc l, thus that is a
 *   boost::multi_array<FRowConstraint, 2>; with two dimensions which are
 *   f_time_horizon, and f_number_arcs entries, where
 *   t = 0, ...,f_time_horizon - 1 and z = 0, ...,f_number_arcs - 1
 *
 *   \f[
 *
 *      f_{t,l} \in [ F^{mn}_{t,l} , F^{mx}_{t,l}]  \quad t \in \mathcal{T},
 *           l \in \mathcal{L}^{hy}                             \quad (9)
 *
 *   \f]
 *
 * - ramp-up and ramp-down constraints: These inequality(10)-(11) indicate
 *   ramp-up and ramp-down constraints at time t and for ach arc l, so each of
 *   them is a boost::multi_array<FRowConstraint, 2>; with two dimensions
 *   which are f_time_horizon, and f_number_arcs entries, where
 *   t = 0, ...,f_time_horizon - 1 and z = 0, ...,f_number_arcs - 1
 *
 *   \f[
 *
 *      f_{t,l} - f_{t-1,l} \leq \Delta^{up}_{t,l} \quad t \in \mathcal{T},
 *           l \in \mathcal{L}^{hy}                             \quad (10)
 *
 *   \f]
 *   \f[
 *
 *      f_{t-1,l} - f_{t,l} \leq \Delta^{dn}_{t,l} \quad t \in \mathcal{T},
 *           l \in \mathcal{L}^{hy}                             \quad (11)
 *
 *   \f]
 *
 * - final volumes of each reservoir constraints: this equality(12) gives the
 *   final volumes of each reservoir r at time t. This is a
 *   boost::multi_array<FRowConstraint, 2>; with two dimensions which are
 *   f_number_reservoirs, and f_time_horizon entries, where
 *   r = 0, ...,f_number_reservoirs - 1 and z = 0, ...,f_time_horizon - 1
 *   \f[
 *
 *      v^{hy}_{n,t} = v^{hy}_{n,t-1} + 3600 A_{n,t} +
 *      3600 (\sum_{n' \in \mathcal{A}(n)}\sum_{ l \in \mathcal{L}^{hy} }
 *      f_{t - \tau^{dn}_l} - \sum_{n' \in \mathcal{F}(n)}
 *      \sum_{ l \in \mathcal{L}^{hy} } f_{t - \tau^{up}_l})
 *      \quad t \in \mathcal{T}, \quad n \in \mathcal{N}^{hy}    \quad (12)
 *
 *   \f]
 *
 * - final volumes variable bounds: This inequality(13) indicates upper and
 *   lower bound of volumetric variables of each reservoir for each time t,
 *   thus that is a boost::multi_array<FRowConstraint, 2>; with two dimensions
 *   which are f_number_reservoirs, and f_time_horizon entries, where
 *   r = 0, ...,f_number_reservoirs - 1 and z = 0, ...,f_time_horizon - 1
 *   \f[
 *
 *      v^{hy}_{n,t} \in [ V^{hy,mn}_{n,t} , V^{hy,mx}_{n,t}]  \quad
 *        n \in \mathcal{N}^{hy}, t \in \mathcal{T}            \quad (13)
 *
 *   \f]
 *
 */
 void generate_abstract_constraints( Configuration *stcc ) override;
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// generate the objective function of the HydroUnitBlock
/** Method that generates the objective function of the HydroUnitBlock.

 * - Objective function: the objective function of this unit is zero */
 void generate_objective( Configuration *objc ) override;

/**@} ----------------------------------------------------------------------*/
/*--------- METHODS FOR READING THE DATA OF THE HydroUnitBlock -------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the HydroUnitBlock
 *
 * These methods allow to read data that must be common to (in principle) all
 * the kind of hydro units
 * @{ */

/// returns the number of reservoirs
 Index get_number_reservoirs() const { return f_number_reservoirs; }

/*--------------------------------------------------------------------------*/
 /// returns the number of arcs
 Index get_number_arcs() const { return f_number_arcs; }

/*--------------------------------------------------------------------------*/
/// returns the vector of start arcs
/** Method for returning the vector of starting point of each arc. This vector
 * may have size of 1 (single hydro unit with just one arc between two
 * reservoirs) or the size of number of reservoirs, then there are two
 * possible cases:
 *
 *  - if f_number_reservoirs == 2, this vector has size of 1 which means there
 *    is one arc at the system (single hydro case).
 *
 *  - if f_number_reservoirs > 2, this vector have size of f_number_reservoirs
 *    and each element of the vectors gives starting point of each arc in the
 *    network.
 */
 const std::vector< Index > & get_start_arc() const { return v_start_arc; }

/*--------------------------------------------------------------------------*/
/// returns the vector of end arcs
/** Method for returning the vector of ending point of each arc. This
 *  vector may have size of 1 (single hydro unit with just one arc between two
 *  reservoirs) or the size of number of reservoirs, then there are two
 *  possible cases:
 *
 *  - if f_number_reservoirs == 2, this vector has size of 1 which means there
 *    is one arc at the system (single hydro case).
 *
 *  - if f_number_reservoirs > 2, this vector have size of f_number_reservoirs
 *    and each element of the vectors gives ending point of each arc in the
 *    network.
 */
 const std::vector< Index > & get_end_arc() const { return( v_end_arc ); }
 /*--------------------------------------------------------------------------*/

/// returns the matrix of inertia power
/** The returned value U = get_inertia_power() contains the contribution to
 *  inertia (basically, the constants to be multiplied by the active power
 *  variables returned by get_active_power()) of eac arcs (generators) at each
 *  time instants. There are three possible cases
 *
 * - if the matrix is empty, then the inertia power is 0;
 *
 * - if the matrix only has one row (i.e., the first dimension has size 1),
 *   then the inertia power for arc (generator) l is U[ 0 , l ] for all t
 *   which means that the second dimension has size get_number_arcs();
 *
 * - otherwise, the matrix has size the time horizon per
 *   get_number_arcs(), and U[ t , l ] contains the contribution to
 *   inertia power of arc(generator) l at time instant t. */
 const boost::multi_array< double , 2 > & get_inertia_power() const override {
  return ( v_inertia_power );
 }
/*--------------------------------------------------------------------------*/
/// returns the matrix of minimum volumetric
/** The method returned a two-dimensional boost::multi_array<> M such that
 * M[ n , t ] gives the minimum volumetric of the reservoir n at the time
 * instant t. This two-dimensional boost::multi_array<> M considers three
 * possible cases:
 *
 *  - if the boost::multi_array<> M is empty() then no minimum volumetric are
 *    defined, and there are no minimum volumetric constraints;
 *
 *  - if the boost::multi_array<> M has only one column which in this case the
 *    boost::multi_array<> M is a (transpose of) vector with size
 *    get_number_reservoirs(). Each element of M[ n , 0 ] gives the minimum
 *    volumetric of the each reservoir n for all time instant t;
 *
 *  - otherwise the two-dimensional boost::multi_array<> M must have
 *    get_number_reservoirs() row where each row must have size of
 *    get_time_horizon() and each element of M[ n , t ] gives the minimum
 *    volumetric of reservoir n at time instant t. */
 const boost::multi_array< double , 2 > & get_minimum_volumetric() const {
  return ( v_minimum_volumetric );
 }
/*--------------------------------------------------------------------------*/
/// returns the matrix of maximum volumetric
/** The method returned a two-dimensional boost::multi_array<> M such that
 * M[ n , t ] gives the maximum volumetric of the reservoir n at the time
 * instant t. This two-dimensional boost::multi_array<> M considers three
 * possible cases:
 *
 *  - if the boost::multi_array<> M is empty() then no maximum volumetric are
 *    defined, and there are no minimum volumetric constraints;
 *
 *  - if the boost::multi_array<> M has only one column which in this case the
 *    boost::multi_array<> M is a (transpose of) vector with size
 *    get_number_reservoirs(). Each element of M[ n , 0 ] gives the maximum
 *    volumetric of the each reservoir n for all time instant t;
 *
 *  - otherwise the two-dimensional boost::multi_array<> M must have
 *    get_number_reservoirs() row where each row must have size of
 *    get_time_horizon() and each element of M[ n , t ] gives the maximum
 *    volumetric of reservoir n at time instant t. */
 const boost::multi_array< double , 2 > & get_maximum_volumetric() const {
  return ( v_maximum_volumetric );
 }
/*--------------------------------------------------------------------------*/
/// returns the matrix of inflows
/** The method returned a two-dimensional boost::multi_array<> M such that
 * M[ n , t ] gives the inflows of the reservoir n at the time
 * instant t. This two-dimensional boost::multi_array<> M considers two
 * possible cases:
 *
 *  - if the boost::multi_array<> M is empty() then no inflows are defined;
 *
 *  - otherwise the two-dimensional boost::multi_array<> M must have
 *    get_number_reservoirs() row where each row must have size of
 *    get_time_horizon() and each element of M[ n , t ] represents the inflows
 *    of reservoir n at time instant t. */
 const boost::multi_array< double , 2 > & get_inflows() const {
  return ( v_inflows );
 }
/*--------------------------------------------------------------------------*/
/// returns the matrix of minimum power
/** The method returned a two-dimensional boost::multi_array<> M such that
 * M[ t , i ] gives the minimum power at each time t associated with unit(arc)
 * i. This two-dimensional boost::multi_array<> M considers three possible
 * cases:
 *
 *  - if the boost::multi_array<> M is empty() then no minimum power are
 *    defined, and there are no minimum power constraints;
 *
 *  - if the boost::multi_array<> M has only one row which in this case the
 *    boost::multi_array<> M is a vector with size get_number_arcs(). Each
 *    element of M[ 0 , i ] gives the minimum power for all time instant t of
 *    each unit i;
 *
 *  - otherwise the two-dimensional boost::multi_array<> M must have
 *    get_time_horizon() rows and get_number_arcs() columns and each element
 *    of M[ t , i ] gives the minimum power at time t and unit i. */
 const boost::multi_array< double , 2 > & get_minimum_power() const {
  return ( v_minimum_power );
 }
/*--------------------------------------------------------------------------*/
/// returns the matrix of maximum power
/** The method returned a two-dimensional boost::multi_array<> M such that
 * M[ t , i ] gives the maximum power at each time t associated with unit(arc)
 * i. This two-dimensional boost::multi_array<> M considers three possible
 * cases:
 *
 *  - if the boost::multi_array<> M is empty() then no maximum power are
 *    defined, and there are no maximum power constraints;
 *
 *  - if the boost::multi_array<> M has only one row which in this case the
 *    boost::multi_array<> M is a vector with size get_number_arcs(). Each
 *    element of M[ 0 , i ] gives the maximum power for all time instant t of
 *    each unit i;
 *
 *  - otherwise the two-dimensional boost::multi_array<> M must have
 *    get_time_horizon() rows and get_number_arcs() columns and each element
 *    of M[ t , i ] gives the maximum power at time t and unit i. */
 const boost::multi_array< double , 2 > & get_maximum_power() const {
  return ( v_maximum_power );
 }
/*--------------------------------------------------------------------------*/
/// returns the matrix of minimum flow
/** The method returned a two-dimensional boost::multi_array<> M such that
 * M[ t , i ] gives the minimum flow at each time t associated with unit(arc)
 * i. This two-dimensional boost::multi_array<> M considers three possible
 * cases:
 *
 *  - if the boost::multi_array<> M is empty() then no minimum flow are
 *    defined, and there are no minimum flow constraints;
 *
 *  - if the boost::multi_array<> M has only one row which in this case the
 *    boost::multi_array<> M is a vector with size get_number_arcs(). Each
 *    element of M[ 0 , i ] gives the minimum flow for all time instant t of
 *    each unit i;
 *
 *  - otherwise the two-dimensional boost::multi_array<> M must have
 *    get_time_horizon() rows and get_number_arcs() columns and each element
 *    of M[ t , i ] gives the minimum flow at time t and unit i. */
 const boost::multi_array< double , 2 > & get_minimum_flow() const {
  return ( v_minimum_flow );
 }
/*--------------------------------------------------------------------------*/
/// returns the matrix of maximum flow
/** The method returned a two-dimensional boost::multi_array<> M such that
 * M[ t , i ] gives the maximum flow at each time t associated with unit(arc)
 * i. This two-dimensional boost::multi_array<> M considers three possible
 * cases:
 *
 *  - if the boost::multi_array<> M is empty() then no maximum flow are
 *    defined, and there are no maximum flow constraints;
 *
 *  - if the boost::multi_array<> M has only one row which in this case the
 *    boost::multi_array<> M is a vector with size get_number_arcs(). Each
 *    element of M[ 0 , i ] gives the maximum flow for all time instant t of
 *    each unit i;
 *
 *  - otherwise the two-dimensional boost::multi_array<> M must have
 *    get_time_horizon() rows and get_number_arcs() columns and each element
 *    of M[ t , i ] gives the maximum flow at time t and unit i. */
 const boost::multi_array< double , 2 > & get_maximum_flow() const {
  return ( v_maximum_flow );
 }
/*--------------------------------------------------------------------------*/
/// returns the matrix of delta ramp up
/** The method returned a two-dimensional boost::multi_array<> M such that
 * M[ t , i ] gives the delta ramp up at each time t associated with unit(arc)
 * i. This two-dimensional boost::multi_array<> M considers three possible
 * cases:
 *
 *  - if the boost::multi_array<> M is empty() then no ramping constraints;
 *
 *  - if the boost::multi_array<> M has only one row which in this case the
 *    boost::multi_array<> M is a vector with size get_number_arcs(). Each
 *    element of M[ 0 , i ] gives the delta ramp up value for all time instant
 *    t of each unit i;
 *
 *  - otherwise the two-dimensional boost::multi_array<> M must have
 *    get_time_horizon() rows and get_number_arcs() columns and each element
 *    of M[ t , i ] gives the delta ramp up value at time t and unit i. */
 const boost::multi_array< double , 2 > & get_delta_ramp_up() const {
  return ( v_delta_ramp_up);
 }
/*--------------------------------------------------------------------------*/
/// returns the matrix of delta ramp down
/** The method returned a two-dimensional boost::multi_array<> M such that
 * M[ t , i ] gives the delta ramp down at each time t associated with unit
 * (arc) i. This two-dimensional boost::multi_array<> M considers three
 * possible cases:
 *
 *  - if the boost::multi_array<> M is empty() then no ramping constraints;
 *
 *  - if the boost::multi_array<> M has only one row which in this case the
 *    boost::multi_array<> M is a vector with size get_number_arcs(). Each
 *    element of M[ 0 , i ] gives the delta ramp down value for all time
 *    instant t of each unit i;
 *
 *  - otherwise the two-dimensional boost::multi_array<> M must have
 *    get_time_horizon() rows and get_number_arcs() columns and each element
 *    of M[ t , i ] gives the delta ramp down value at time t and unit i. */
 const boost::multi_array< double , 2 > & get_delta_ramp_down() const {
  return ( v_delta_ramp_down);
 }
/*--------------------------------------------------------------------------*/
/// returns the matrix of primary rho
/** The method returned a two-dimensional boost::multi_array<> M such that
 * M[ t , i ] gives the primary rho at each time t associated with unit
 * (arc) i. This two-dimensional boost::multi_array<> M considers three
 * possible cases:
 *
 *  - if the boost::multi_array<> M is empty() then no primary reserve
 *    constraints;
 *
 *  - if the boost::multi_array<> M has only one row which in this case the
 *    boost::multi_array<> M is a vector with size get_number_arcs(). Each
 *    element of M[ 0 , i ] gives the primary rho value for all time instant t
 *    of each unit i;
 *
 *  - otherwise the two-dimensional boost::multi_array<> M must have
 *    get_time_horizon() rows and get_number_arcs() columns and each element
 *    of M[ t , i ] gives the primary rho value at time t and unit i. */
 const boost::multi_array< double , 2 > & get_primary_rho() const {
  return ( v_primary_rho);
 }
/*--------------------------------------------------------------------------*/
/// returns the matrix of secondary rho
/** The method returned a two-dimensional boost::multi_array<> M such that
 * M[ t , i ] gives the secondary rho at each time t associated with unit
 * (arc) i. This two-dimensional boost::multi_array<> M considers three
 * possible cases:
 *
 *  - if the boost::multi_array<> M is empty() then no secondary reserve
 *    constraints;
 *
 *  - if the boost::multi_array<> M has only one row which in this case the
 *    boost::multi_array<> M is a vector with size get_number_arcs(). Each
 *    element of M[ 0 , i ] gives the secondary rho value for all time instant
 *    t of each unit i;
 *
 *  - otherwise the two-dimensional boost::multi_array<> M must have
 *    get_time_horizon() rows and get_number_arcs() columns and each element
 *    of M[ t , i ] gives the secondary rho value at time t and unit i. */
 const boost::multi_array< double , 2 > & get_secondary_rho() const {
  return ( v_secondary_rho);
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
 const std::vector< Index > & get_number_pieces() const {
  return ( v_number_pieces);
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
 const std::vector< double > & get_const_term() const {
  return ( v_const_term);
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
 const std::vector< double > & get_linear_term() const {
  return ( v_linear_term);
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
 const std::vector< Index > & get_uphill_delay() const {
  return ( v_uphill_delay);
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
 const std::vector< Index > & get_downhill_delay() const {
  return ( v_downhill_delay);
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
 *   */
 const std::vector< double > & get_initial_volumetric() const {
  return ( v_initial_volumetric);
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
 const std::vector< double > & get_initial_flow_rate() const {
  return ( v_initial_flow_rate);
 }
/**@} ----------------------------------------------------------------------*/
/*---------- METHODS FOR READING THE Variable OF THE HydroUnitBlock --------*/
/*--------------------------------------------------------------------------*/

/** @name Reading the Variable of the HydroUnitBlock
 *
 * These methods allow to read the two groups of Variable that any
 * HydroUnitBlock in principle has (although some may not):
 *
 * - the volumetric variables;
 *
 * - the flow rate variables;
 *
 * All these two groups of variables are (if not empty)
 * boost::multi_array< ColVariable , 2 > with first dimension time horizon
 * and second dimension number of arcs (generators).
 * @{ */

/// returns the matrix of volumetric variables
/** The returned boost::multi_array< ColVariable , 2 >, say V, contains the
 * volumetric variables and is indexed over the dimensions time horizon and
 * number of arcs (generators). There are two possible cases:
 *
 * - if V is empty(), then these variables are not defined;
 *
 * - otherwise, V must have f_time_horizon rows and f_number_arcs columns, and
 *   M[ t , a ] is the volumetric variable for time step t of arc (generator)
 *   a. */

 const boost::multi_array< ColVariable , 2 > & get_volumetric() const {
  return v_volumetric;
 }
 /*--------------------------------------------------------------------------*/
/// returns the matrix of flow rate variables
/** The returned boost::multi_array< ColVariable , 2 >, say F, contains the
 * flow rate variables and is indexed over the dimensions time horizon and
 * number of arcs (generators). There are two possible cases:
 *
 * - if V is empty(), then these variables are not defined;
 *
 * - otherwise, F must have f_time_horizon rows and f_number_arcs columns, and
 *   M[ t , a ] is the flow rate variable for time step t of arc (generator)
 *   a. */

 const boost::multi_array< ColVariable , 2 > & get_flow_rate() const {
  return v_flow_rate;
 }
/**@} ----------------------------------------------------------------------*/
/*------------------ METHODS FOR SAVING THE HydroUnitBlock------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the HydroUnitBlock
 *  @{ */

/// extends Block::serialize( netCDF::NcGroup )
/** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
 * HydroUnitBlock. See HydroUnitBlock::deserialize( netCDF::NcGroup ) for
 * details of the format of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/**@} ----------------------------------------------------------------------*/
/*--------------- METHODS FOR INITIALIZING THE HydroUnitBlock --------------*/
/*--------------------------------------------------------------------------*/

/** @name Handling the data of the HydroUnitBlock
    @{ */

 void load( std::istream & input ) override {
  throw ( std::logic_error( "HydroUnitBlock::load() not implemented yet") );
 };

/**@} ----------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------data--------------------------------------*/
 /// The number of reservoirs(nodes) of the problem
 Index f_number_reservoirs;

 /// The number of connecting arcs which are connecting the reservoirs in
 /// cascading system
 Index f_number_arcs;

 /// The vector of UphillDelay
 std::vector< Index > v_uphill_delay;

 /// The vector of DownhillDelay
 std::vector< Index > v_downhill_delay;

 /// The vector of starting arc
 std::vector< Index > v_start_arc;

 /// The vector of ending arcs
 std::vector< Index > v_end_arc;

 /// The vector of initial volumetric
 std::vector< double > v_initial_volumetric;

 /// The vector of initial flow rate
 std::vector< double > v_initial_flow_rate;

 /// The vector of NumberPieces
 std::vector< Index > v_number_pieces;

 /// The vector of LinearTerm
 std::vector< double > v_linear_term;

 /// The vector of ConstTerm
 std::vector< double > v_const_term;

 /// the vector of inertia power of generators
 boost::multi_array< double , 2 > v_inertia_power;

 /// The matrix of MinVolumetric
 /** Indexed over the dimensions NumberReservoirs and NumberIntervals. */
 boost::multi_array< double, 2 > v_minimum_volumetric;

 /// The matrix of MaxVolumetric
 /** Indexed over the dimensions NumberReservoirs and NumberIntervals. */
 boost::multi_array< double, 2 > v_maximum_volumetric;

 /// The matrix of Inflows
 /** Indexed over the dimensions NumberReservoirs and NumberIntervals. */
 boost::multi_array< double, 2 > v_inflows;

 /// The matrix of MinPower
 /** Indexed over the dimensions NumberIntervals and NumberArcs. */
 boost::multi_array< double, 2 > v_minimum_power;

 /// The matrix of MaxPower
 /** Indexed over the dimensions NumberIntervals and NumberArcs. */
 boost::multi_array< double, 2 > v_maximum_power;

 /// The matrix of MinFlow
 /** Indexed over the dimensions NumberIntervals and NumberArcs. */
 boost::multi_array< double, 2 > v_minimum_flow;

 /// The matrix of MaxFlow
 /** Indexed over the dimensions NumberIntervals and NumberArcs. */
 boost::multi_array< double, 2 > v_maximum_flow;

 /// The matrix of DeltaRampUp
 /** Indexed over the dimensions NumberIntervals and NumberGenerators. */
 boost::multi_array< double, 2 > v_delta_ramp_up;

 /// The matrix of DeltaRampDown
 /** Indexed over the dimensions NumberIntervals and NumberGenerators. */
 boost::multi_array< double, 2 > v_delta_ramp_down;

 /// The matrix of PrimaryRho
 /** Indexed over the dimensions NumberIntervals and NumberArcs. */
 boost::multi_array< double, 2 > v_primary_rho;

 /// The matrix of SecondaryRho
 /** Indexed over the dimensions NumberIntervals and NumberArcs. */
 boost::multi_array< double, 2 > v_secondary_rho;

/*-----------------------------variables------------------------------------*/

 /// the matrix of volumetric variables
 boost::multi_array< ColVariable , 2> v_volumetric;

 /// the matrix of flow rate variables
 boost::multi_array< ColVariable , 2> v_flow_rate;

/*----------------------------constraints-----------------------------------*/
 /// maximum power output according to primary-secondary reserves constraints
 boost::multi_array< FRowConstraint, 2 >  v_MaxPowerPrimarySecondary_Const;

 /// minimum power output according to primary-secondary reserves constraints
 boost::multi_array< FRowConstraint, 2 >  v_MinPowerPrimarySecondary_Const;

 /// power output relation with to primary reserves constraints
 boost::multi_array< FRowConstraint, 2 >  v_ActivePowerPrimary_Const;

 /// power output relation with to secondary reserves constraints
 boost::multi_array< FRowConstraint, 2 >  v_ActivePowerSecondary_Const;

 /// primary reserves constraints for pumps
 boost::multi_array< FRowConstraint, 2 >  v_PrimaryPumps_Const;

 /// secondary reserves constraints for pumps
 boost::multi_array< FRowConstraint, 2 >  v_SecondaryPumps_Const;

 /// flow to active power function constraints for pumps
 boost::multi_array< FRowConstraint, 2 >  v_FlowActivePowerPumps_Const;

 /// flow to active power function constraints for turbine
 boost::multi_array< FRowConstraint, 2 >  v_FlowActivePowerTurbines_Const;

 /// ramp-up constraints
 boost::multi_array< FRowConstraint, 2 >  v_RampUp_Const;

 /// ramp-down constraints
 boost::multi_array< FRowConstraint, 2 >  v_RampDown_Const;

 /// flow rate bounds constraints
 boost::multi_array< FRowConstraint, 2 >  v_FlowRateBounds_Const;

 /// final volumes fo each reservoir constraints
 boost::multi_array< FRowConstraint, 2 >  v_FinalVolumeReservoir_Const;

 /// volumetric bounds constraints
 boost::multi_array< FRowConstraint, 2 >  v_Volumetric_Const;

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

};  // end( class( HydroUnitBlock ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* HydroUnitBlock.h included */

/*--------------------------------------------------------------------------*/
/*---------------------- End File HydroUnitBlock.h -------------------------*/
/*--------------------------------------------------------------------------*/
