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
 * \date 11 - 07 - 2019
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
 * - maximum and minimum power output constraints according to primary and
 *   secondary spinning reserves;
 *
 * - active power relation with primary and secondary spinning reserves;
 *
 * - primary and secondary spinning reserves value for pumps( == 0 );
 *
 * - flow to active power function;
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
 *   (or nodes) in a hydro unit block. The dimension is optional, if it is not
 *   provided then it is taken to be == 1 and in this case the cascading
 *   system becomes to a single hydro unit.
 *
 * Note: The concept of "fake" reservoir, with no volumetric variable and no
 * volumetric constraint /todo
 *
 * - The dimension "NumberArcs" containing the set of arcs connecting the
 *   reservoirs in cascading system.
 *
 * - The variable "StartArc", of type int and indexed over the dimension
 *   "NumberReservoirs"; the r-th entry of the variable is the starting point
 *   of the arc (a number in 0, ..., NumberReservoirs - 1). Note that arcs are
 *   oriented; that is, a positive flow along arc r (turbine) means that water
 *   is being taken away from StartArc[ r ] and delivered to EndArc[ r ]
 *   (see next), a negative flow (pump) means vice-versa. Note that reservoir
 *   names here go from 0 to NumberReservoirs.getSize() - 1;
 *
 * - The variable "EndArc", of type int and indexed over the dimension
 *   "NumberReservoirs"; the r-th entry of the variable is the ending point
 *   of the arc (a number in 0, ..., NumberReservoirs - 1). Note that arcs are
 *   oriented (see above); StartArc[ r ] == EndArc[ r ] (a self-loop) is not
 *   allowed, but multiple arcs between the same pair of reservoirs are. Note
 *   that reservoir names here go from 0 to NumberReservoirs.getSize() - 1;
 *
 * - The variable "MinFlow", of type double and indexed over both dimensions
 *   "NumberIntervals" and "NumberArcs". Both dimensions may have either size
 *   1 or full size("NumberIntervals" and "NumberArcs", respectively). This is
 *   meant to represent the matrix MinF[ t , a ] which, for each time instant
 *   t at each arc a contains the minimum flow value of the unit; it must be
 *   that the entry MinF[ t , a ] <= MaxF[ t , a ] for each time instant t and
 *   each arc a. If both dimensions have size 1 then the entry MinF[ 0 , 0 ]
 *   gives the minimum flow value of the unit with only one existing arc for
 *   all the time steps. If firs dimension has size 1 then the entry
 *   MinF[ 0 , a ] is assumed to contain the minimum flow of each arc a for
 *   all time instant. Otherwise, two cases may happen such that both
 *   dimensions may have full size or second dimension could have size of 1
 *   then MinFlow[ i , a ] (or MinFlow[ i , 0 ]) is the fixed value of
 *   MinF[ t , a ] (or MinF[ t , 0 ]) for all time t and arc a (or the one
 *   available arc) in the interval
 *   [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the assumption
 *   that ChangeIntervals[ - 1 ] = 0. If "NumberIntervals" <= 1 or
 *   "NumberIntervals" >= "TimeHorizon" then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded.
 *
 * - The variable "MaxFlow", of type double and indexed over both dimensions
 *   "NumberIntervals" and "NumberArcs". Both dimensions may have either size
 *   1 or full size("NumberIntervals" and "NumberArcs", respectively). This is
 *   meant to represent the matrix MaxF[ t , a ] which, for each time instant
 *   t at each arc a contains the maximum flow value of the unit; it must be
 *   that the entry MinF[ t , a ] <= MaxF[ t , a ] for each time instant t and
 *   each arc a. If both dimensions have size 1 then the entry MaxF[ 0 , 0 ]
 *   gives the maximum flow value of the unit with only one existing arc for
 *   all the time steps. If firs dimension has size 1 then the entry
 *   MaxF[ 0 , a ] is assumed to contain the maximum flow of each arc a for
 *   all time instant. Otherwise, two cases may happen such that both
 *   dimensions may have full size or second dimension could have size of 1
 *   then MaxFlow[ i , a ] (or MaxFlow[ i , 0 ]) is the fixed value of
 *   MaxF[ t , a ] (or MaxF[ t , 0 ]) for all time t and arc a (or the one
 *   available arc) in the interval
 *   [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the assumption
 *   that ChangeIntervals[ - 1 ] = 0. If "NumberIntervals" <= 1 or
 *   "NumberIntervals" >= "TimeHorizon" then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded.
 *
 * Note: MinFlow and MaxFlow values can be either positive or negative (or
 * zero); whenever MinF[ t , a ] < 0 and MaxF[ t , a ] < 0 for each t and a,
 * the unit is considered a pump and whenever MinF[ t , a ] > 0 and
 * MaxF[ t , a ] > 0 the unit is considered a turbine. Any possible mixed
 * situation can be accounted for by artificially splitting the unit into “two
 * units”, which should be done at the data processing stage. When
 * MaxF[ t , a ] == 0 it means maximum flow rate of the unit (which is
 * considered as a pump) is zero, and when MinF[ t , a ] == 0 it means minimum
 * flow rate of the unit(which is considered as turbine) is zero.
 *
 * - The variable "MinVolumetric", of type double and indexed over both
 *   dimensions "NumberReservoirs" and "NumberIntervals". Both dimensions may
 *   have either size 1 or full size("NumberReservoirs" and "NumberIntervals",
 *   respectively). This is meant to represent the matrix MinV[ r , t ] which,
 *   for each reservoir r at each time instant t contains the minimum
 *   volumetric value of the unit for each reservoir and corresponding time
 *   step; it must be that the entry MinV[ r , t ] >= 0 for all r and t and
 *   MinV[ r , t ] <= MaxV[ r , t ]. If both dimensions have size 1 then the
 *   entry MinV[ 0 , 0 ] gives the minimum volumetric value of the unit with
 *   only one existing reservoir for all the time steps. If second dimension
 *   has size 1 then the entry MinV[ r , 0 ] is assumed to contain the minimum
 *   volumetric for each reservoir r for all the time instants. Otherwise, two
 *   cases may happen in which both dimensions may have full size or first
 *   dimension has size of 1 then MinVolumetric[ r , i ] (or
 *   MinVolumetric[ 0 , i ] ) is the fixed value of MinV[ r , t ] (or
 *   MinV[ 0 , t ]) for all reservoir r (the one available reservoir)and all
 *   t in the interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ],
 *   with the assumption that ChangeIntervals[ - 1 ] = 0. If
 *   "NumberIntervals" <= 1 or "NumberIntervals" >= "TimeHorizon" then the
 *   mapping clearly does not require "ChangeIntervals", which in fact is not
 *   loaded.
 *
 * - The variable "MaxVolumetric", of type double and indexed over both
 *   dimensions "NumberReservoirs" and "NumberIntervals". Both dimensions may
 *   have either size 1 or full size("NumberReservoirs" and "NumberIntervals",
 *   respectively). This is meant to represent the matrix MaxV[ r , t ] which,
 *   for each reservoir r at each time instant t contains the maximum
 *   volumetric value of the unit for each reservoir and corresponding time
 *   step; it must be that the entry MaxV[ r , t ] >= 0 for all r and t and
 *   MinV[ r , t ] <= MaxV[ r , t ]. If both dimensions have size 1 then the
 *   entry MaxV[ 0 , 0 ] gives the maximum volumetric value of the unit with
 *   only one existing reservoir for all the time steps. If second dimension
 *   has size 1 then the entry MaxV[ r , 0 ] is assumed to contain the maximum
 *   volumetric for each reservoir r for all the time instants. Otherwise, two
 *   cases may happen in which both dimensions may have full size or first
 *   dimension has size of 1 then MaxVolumetric[ r , i] (or
 *   MaxVolumetric[ 0 , i ] ) is the fixed value of MaxV[ r , t ] (or
 *   MaxV[ 0 , t ]) for all reservoir r (the one available reservoir)and all
 *   t in the interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ],
 *   with the assumption that ChangeIntervals[ - 1 ] = 0. If
 *   "NumberIntervals" <= 1 or "NumberIntervals" >= "TimeHorizon" then the
 *   mapping clearly does not require "ChangeIntervals", which in fact is not
 *   loaded.
 *
 * - The variable "Inflows", of type double and indexed over both dimensions
 *   "NumberReservoirs" and "NumberIntervals". Both dimensions may have either
 *   size 1 or full size("NumberReservoirs" and "NumberIntervals",
 *   respectively). This is meant to represent the matrix InF[ r , t ] which,
 *   for each reservoir r at each time instant t contains the amount of water
 *   that goes to each reservoir r at time t;; it must be that the entry
 *   InF[ r , t ] >= 0 for all r and t. If both dimensions have size 1 then the
 *   entry InF[ 0 , 0 ] gives the amount of water that goes to one existing
 *   reservoir for all the time steps. If second dimension has size 1 then the
 *   entry InF[ r , 0 ] is assumed to contain the amount of water that goes to
 *   each reservoir r for all the time instants. Otherwise, two cases may
 *   happen in which both dimensions may have full size or first
 *   dimension has size of 1 then Inflows[ r , i ] (or Inflows[ 0 , i ] ) is
 *   the fixed value of InF[ r , t ] (or InF[ 0 , t ]) for all reservoir n
 *   (the one available reservoir) and all t in the interval
 *   [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the assumption
 *   that ChangeIntervals[ - 1 ] = 0. If "NumberIntervals" <= 1 or
 *   "NumberIntervals" >= "TimeHorizon" then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded.
 *
 * - The variable "MinPower", of type double and indexed over both dimensions
 *   "NumberIntervals" and "NumberArcs". Both dimensions may have either size
 *   1 or full size("NumberIntervals" and "NumberArcs", respectively). This is
 *   meant to represent the matrix MinP[ t , a ] which, for each time instant
 *   t at each arc a contains the minimum power value of the unit; it must be
 *   that MinP[ t , a ] >= 0 and MinP[ t , a ] <= MaxP[ t , a ] for each time
 *   instant t and each arc a. If both dimensions have size 1 then the entry
 *   MinP[ 0 , 0 ] gives the minimum power value of the unit with only one
 *   existing arc for all the time steps. If firs dimension has size 1 then
 *   the entry MinP[ 0 , a ] is assumed to contain the minimum power of each
 *   arc a for all time instant. Otherwise, two cases may happen such that
 *   both dimensions may have full size or second dimension could have size of
 *   1 then MinPower[ i , a ] (or MinPower[ i , 0 ]) is the fixed value of
 *   MinP[ t , a ] (or MinP[ t , 0 ]) for all time t and arc a (or the one
 *   available arc) in the interval
 *   [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the assumption
 *   that ChangeIntervals[ - 1 ] = 0. If "NumberIntervals" <= 1 or
 *   "NumberIntervals" >= "TimeHorizon" then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded.
 *
 * - The variable "MaxPower", of type double and indexed over both dimensions
 *   "NumberIntervals" and "NumberArcs". Both dimensions may have either size
 *   1 or full size("NumberIntervals" and "NumberArcs", respectively). This is
 *   meant to represent the matrix MaxP[ t , a ] which, for each time instant
 *   t at each arc a contains the maximum power value of the unit; it must be
 *   that MaxP[ t , a ] >= 0 and MinP[ t , a ] <= MaxP[ t , a ] for each time
 *   instant t and each arc a. If both dimensions have size 1 then the entry
 *   MaxP[ 0 , 0 ] gives the maximum power value of the unit with only one
 *   existing arc for all the time steps. If firs dimension has size 1 then
 *   the entry MaxP[ 0 , a ] is assumed to contain the maximum power of each
 *   arc a for all time instant. Otherwise, two cases may happen such that
 *   both dimensions may have full size or second dimension could have size of
 *   1 then MaxPower[ i , a ] (or MaxPower[ i , 0 ]) is the fixed value of
 *   MaxP[ t , a ] (or MaxP[ t , 0 ]) for all time t and arc a (or the one
 *   available arc) in the interval
 *   [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the assumption
 *   that ChangeIntervals[ - 1 ] = 0. If "NumberIntervals" <= 1 or
 *   "NumberIntervals" >= "TimeHorizon" then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded.
 *
 * - The variable "DeltaRampUp", of type double and indexed over both
 *   dimensions "NumberIntervals" and "NumberArcs". Both dimensions may have
 *   either size 1 or full size("NumberIntervals" and "NumberArcs",
 *   respectively). This is meant to represent the matrix DP[ t , a ] which,
 *   contains the maximum possible increase of the flow rate at each time
 *   instant t of each arc a. If both dimensions have size 1 then the entry
 *   DP[ 0 , 0 ] gives the maximum possible increase of the flow rate value of
 *   the only one existing arc for all the time steps. If firs dimension has
 *   size 1 then the entry DP[ 0 , a ] is assumed to contain the maximum
 *   possible increase of the flow rate of each arc a for all time instant.
 *   Otherwise, two cases may happen such that both dimensions may have full
 *   size or second dimension could have size of 1 then DeltaRampUp[ i , a ]
 *   (or DeltaRampUp[ i , 0 ]) is the fixed value of DP[ t , a ] (or
 *   DP[ t , 0 ]) for all time t and arc a (or the one available arc) in the
 *   interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the
 *   assumption that ChangeIntervals[ - 1 ] = 0. If "NumberIntervals" <= 1 or
 *   "NumberIntervals" >= "TimeHorizon" then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded. This variable is
 *   optional; if it is not provided then it is assumed that
 *   DP[ t , a ] == MaxF[ t , a ], i.e., the unit can ramp up by an arbitrary
 *   amount, i.e., there are no ramp-up constraints.
 *
 * - The variable "DeltaRampDown", of type double and indexed over both
 *   dimensions "NumberIntervals" and "NumberArcs". Both dimensions may have
 *   either size 1 or full size("NumberIntervals" and "NumberArcs",
 *   respectively). This is meant to represent the matrix DM[ t , a ] which,
 *   contains the maximum possible decrease of the flow rate at each time
 *   instant t of each arc a. If both dimensions have size 1 then the entry
 *   DM[ 0 , 0 ] gives the maximum possible decrease of the flow rate value of
 *   the only one existing arc for all the time steps. If firs dimension has
 *   size 1 then the entry DM[ 0 , a ] is assumed to contain the maximum
 *   possible decrease of the flow rate of each arc a for all time instant.
 *   Otherwise, two cases may happen such that both dimensions may have full
 *   size or second dimension could have size of 1 then DeltaRampDown[ i , a ]
 *   (or DeltaRampDown[ i , 0 ]) is the fixed value of DM[ t , a ] (or
 *   DM[ t , 0 ]) for all time t and arc a (or the one available arc) in the
 *   interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the
 *   assumption that ChangeIntervals[ - 1 ] = 0. If "NumberIntervals" <= 1 or
 *   "NumberIntervals" >= "TimeHorizon" then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded. This variable is
 *   optional; if it is not provided then it is assumed that
 *   DM[ t , a ] == MaxF[ t , a ], i.e., the unit can ramp up by an arbitrary
 *   amount, i.e., there are no ramp-up constraints.
 *
 * - The variable "PrimaryRho", of type double and indexed both over the
 *   dimensions "NumberIntervals" and "NumberArcs". Both dimensions may have
 *   either size 1 or full size("NumberIntervals" and "NumberArcs",
 *   respectively). This is meant to represent the matrix PR[ t , a ] which,
 *   for each time instant t and arc a, contains the maximum possible fraction
 *   of active power that can be used as primary reserve. If both dimensions
 *   have size 1 then the entry PR[ 0 , 0 ] gives the maximum possible
 *   fraction of active power that can use as primary reserve for the only one
 *   existing arc for all the time steps. If firs dimension has size 1 then
 *   the entry PR[ 0 , a ] is assumed to contain the maximum
 *   possible fraction of active power that can use as primary reserve of each
 *   arc a for all time instant. Otherwise, two cases may happen such that
 *   both dimensions may have full size or second dimension could have size of
 *   1 then PrimaryRho[ i , a ] (or PrimaryRho[ i , 0 ]) is the fixed value of
 *   PR[ t , a ] (or PR[ t , 0 ]) for all t in the interval
 *   [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the assumption
 *   that ChangeIntervals[ - 1 ] = 0 and all a. If "NumberIntervals" <= 1 or
 *   "NumberIntervals" >= "TimeHorizon" then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded. This variable is
 *   optional, when it's not present it will not capable of producing any
 *   primary reserve, which correspond to PR[ t , a ] == 0 for all t and a.
 *
 * - The variable "SecondaryRho", of type double and indexed both over the
 *   dimensions "NumberIntervals" and "NumberArcs". Both dimensions may have
 *   either size 1 or full size("NumberIntervals" and "NumberArcs",
 *   respectively). This is meant to represent the matrix SR[ t , a ] which,
 *   for each time instant t and arc a, contains the maximum possible fraction
 *   of active power that can be used as secondary reserve. If both dimensions
 *   have size 1 then the entry SR[ 0 , 0 ] gives the maximum possible
 *   fraction of active power that can use as secondary reserve for the only
 *   one existing arc for all the time steps. If firs dimension has size 1
 *   then the entry SR[ 0 , a ] is assumed to contain the maximum possible
 *   fraction of active power that can use as secondary reserve of each arc a
 *   for all time instant. Otherwise, two cases may happen such that both
 *   dimensions may have full size or second dimension could have size of
 *   1 then SecondaryRho[ i , a ] (or SecondaryRho[ i , 0 ]) is the fixed
 *   value of SR[ t , a ] (or SR[ t , 0 ]) for all t in the interval
 *   [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the assumption
 *   that ChangeIntervals[ - 1 ] = 0 and all a. If "NumberIntervals" <= 1 or
 *   "NumberIntervals" >= "TimeHorizon" then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded. This variable is
 *   optional, when it's not present it will not capable of producing any
 *   secondary reserve, which correspond to SR[ t , a ] == 0 for all t and a.
 *
 * - The variable "PowerFlowRho", of type double and indexed both over the
 *   dimensions "NumberIntervals" and "NumberArcs". Both dimensions may have
 *   either size 1 or full size("NumberIntervals" and "NumberArcs",
 *   respectively). This is meant to represent the matrix PFR[ t , a ] which,
 *   for each time instant t and arc a, contains the exact fraction of active
 *   power that can be used as flow rate. If both dimensions have size 1 then
 *   the entry PFR[ 0 , 0 ] gives the exact fraction of active power that can
 *   be used as flow rate for the only existing arc of all the time steps. If
 *   firs dimension has size 1 then the entry PFR[ 0 , a ] is assumed to
 *   contain the the exact fraction of active power that can be used as flow
 *   rate of each arc a for all time instant. Otherwise, two cases may happen
 *   such that both dimensions may have full size or second dimension could
 *   have size of 1 then PowerFlowRho[ i , a ] (or PowerFlowRho[ i , 0 ]) is
 *   the fixed value of PFR[ t , a ] (or PFR[ t , 0 ]) for all t in the
 *   interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the
 *   assumption that ChangeIntervals[ - 1 ] = 0 and all a. If
 *   "NumberIntervals" <= 1 or "NumberIntervals" >= "TimeHorizon" then the
 *   mapping clearly does not require "ChangeIntervals", which in fact is not
 *   loaded. This variable is optional, when it's not present it will not
 *   capable of producing any flow rate, which correspond to PFR[ t , a ] == 0
 *   for all t and a.
 *
 * - The variable "InertiaPower", of type double and indexed both over the
 *   dimensions "NumberIntervals" and "NumberArcs". Both dimensions may have
 *   either size 1 or full size("NumberIntervals" and "NumberArcs",
 *   respectively). This is meant to represent the matrix IP[ t , a ] which,
 *   for each time instant t and arc a, contains the contribution that the
 *   unit can give to the inertia constraint which depends on the active power
 *   that it is currently generating (basically, the constant to be multiplied
 *   to the active power variable) at time t for arc a. If both dimensions
 *   have size 1 then the entry IP[ 0 , 0 ] gives the contribution that the
 *   unit can give to the inertia constraint which depends on the active power
 *   that it is currently generating for the only existing arc of all the time
 *   steps. If firs dimension has size 1 then the entry IP[ 0 , a ] is assumed
 *   to contain the the inertia power value for each arc a and all time t.
 *   Otherwise, two cases may happen such that both dimensions may have full
 *   size or second dimension could have size of 1 then InertiaPower[ i , a ]
 *   (or InertiaPower[ i , 0 ]) is the fixed value of IP[ t , a ] (or
 *   IP[ t , 0 ]) for all t in the interval
 *   [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the assumption
 *   that ChangeIntervals[ - 1 ] = 0 and all a. If  "NumberIntervals" <= 1 or
 *   "NumberIntervals" >= "TimeHorizon" then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded. The variable is
 *   optional; if it is not defined, IP[ t , a ] == 0 for each time instants t
 *   and arc a.
 *
 * - The positive scalar variable "UphillFlow", of type UInt64 and not indexed
 *   over any dimension, which indicates the uphill flow delay in this unit.
 *   This variable is optional, if it is not provided it is taken to be
 *   UphillFlow == 0, which means that ?? //todo.
 *
 * - The positive scalar variable "DownhillFlow", of type UInt64 and not
 *   indexed over any dimension, which indicates the downhill flow delay in
 *   this unit. This variable is optional, if it is not provided it is taken
 *   to be DownhillFlow == 0, which means that ?? //todo.
 * */

 void deserialize( netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
/// Generate the abstract variables of the HydroUnitBlock
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
/// Generate the static constraint of the HydroUnit
/** Method that generates the static constraint of the HydroUnitBlock.
 * These are the:
 * //TODO I should put all the mathematical constraints here
 *
*/
 void generate_abstract_constraints( Configuration *stcc ) override;
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// Generate the objective of the HydroUnitBlock
/** Method that generates the objective of the HydroUnitBlock.
 * //TODO I should put the objective function here
 *
*/
 void generate_objective( Configuration *objc ) override;

/**@} ----------------------------------------------------------------------*/
/*--------- METHODS FOR READING THE DATA OF THE HydroUnitBlock -------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the HydroUnitBlock
 *
 * These methods allow to read data that must be common to (in principle) all
 * the kind of hydro units
 * @{ */

/// returns the number of arcs
 Index get_number_arcs() const { return f_number_arcs; }
/*--------------------------------------------------------------------------*/
/// returns the vector of start arcs
/** Method for returning the vector of starting point of each arc. This
 *  vector may have size of 1 (single hydro unit with just one arc between two
 *  reservoirs) or the size of number of reservoirs, then there are two
 *  possible cases:
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
 *  time instants. There are four possible cases
 *
 * - if the matrix is empty, then the inertia power is 0;
 *
 * - if the matrix only has one row (i.e., the first dimension has size 1),
 *   then the inertia power for arc (generator) a is U[ 0 , a ] for all t
 *   which means that the second dimension has size get_number_arcs();
 *
 * - if the matrix only has one column (i.e., the second dimension has size
 *   1), then the inertia power for all arcs (generators) at time t is
 *   U[ t , 0 ]; which means that the first dimension has size = the time
 *   horizon;
 *
 * - otherwise, the matrix has size the time horizon per
 *   get_number_arcs(), and U[ t , a ] contains the contribution to
 *   inertia power of arc(generator) a at time instant t. */
 const boost::multi_array< double , 2 > & get_inertia_power() const override {
  return ( v_inertia_power );
 }
/*--------------------------------------------------------------------------*/

 const boost::multi_array< double , 2 > & get_minimum_volumetric() const {
  return ( v_minimum_volumetric );
 }
/*--------------------------------------------------------------------------*/

 const boost::multi_array< double , 2 > & get_maximum_volumetric() const {
  return ( v_maximum_volumetric );
 }
/*--------------------------------------------------------------------------*/
 const boost::multi_array< double , 2 > & get_inflows() const {
  return ( v_inflows );
 }
/*--------------------------------------------------------------------------*/

 const boost::multi_array< double , 2 > & get_minimum_power() const {
  return ( v_minimum_power );
 }
/*--------------------------------------------------------------------------*/
 const boost::multi_array< double , 2 > & get_maximum_power() const {
  return ( v_maximum_power );
 }
/*--------------------------------------------------------------------------*/
 const boost::multi_array< double , 2 > & get_minimum_flow() const {
  return ( v_minimum_flow );
 }
/*--------------------------------------------------------------------------*/
 const boost::multi_array< double , 2 > & get_maximum_flow() const {
  return ( v_maximum_flow );
 }
/*--------------------------------------------------------------------------*/
 const boost::multi_array< double , 2 > & get_delta_ramp_up() const {
  return ( v_delta_ramp_up);
 }
/*--------------------------------------------------------------------------*/
 const boost::multi_array< double , 2 > & get_delta_ramp_down() const {
  return ( v_delta_ramp_down);
 }
/*--------------------------------------------------------------------------*/
 const boost::multi_array< double , 2 > & get_primary_rho() const {
  return ( v_primary_rho);
 }
/*--------------------------------------------------------------------------*/
 const boost::multi_array< double , 2 > & get_secondary_rho() const {
  return ( v_secondary_rho);
 }
/*--------------------------------------------------------------------------*/
 const boost::multi_array< double , 2 > & get_power_flow_rho() const {
  return ( v_power_flow_rho);
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

 /// the UphillDelay value
 Index f_uphill_delay;

 /// the DownhillDelay value
 Index f_downhill_delay;

 /// The vector of starting arc
 std::vector< Index > v_start_arc;

 /// The vector of ending arcs
 std::vector< Index > v_end_arc;

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

 /// The matrix of PowerFlowRho
 /** Indexed over the dimensions NumberIntervals and NumberArcs. */
 boost::multi_array< double, 2 > v_power_flow_rho;

/*-----------------------------variables------------------------------------*/

 /// the matrix of volumetric variables
 boost::multi_array< ColVariable , 2> v_volumetric;

 /// the matrix of flow rate variables
 boost::multi_array< ColVariable , 2> v_flow_rate;

/*----------------------------constraints-----------------------------------*/
//TODO

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
