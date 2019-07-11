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
#include "Block.h"

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
/// Implementation of the Block concept for the hydro unit problem
/** The HydroUnitBlock class implements the Block concept [see Block.h] for a
 * "reasonably standard" hydro unit of a Unit Commitment Problem. That is, the
 * class is designed in order to give mathematical formulation to describe the
 * operation of large set of hydro storage. To model complex reservoir systems
 * several technical parameters have to be considered. These are divided into
 * reservoir-specific parameters, the hydro links connecting the reservoirs
 * and finally the turbine/pump parameters. The values are collected within a
 * reservoir database, a hydro-link database and a turbine/pump-database. The
 * technical and physical constraints are mainly divided in ?? different
 * categories:
 * - ??
 * -??
 * -??
 * -??
 */
class HydroUnitBlock : public UnitBlock {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:
/*--------------------------------------------------------------------------*/
/*---------------------- PUBLIC TYPES OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
 *
 * HydroUnitBlock defines a main public types:
 *
 * - Index, the type of parameters indices;
 *
 *
 @{ */

/*--------------------------------------------------------------------------*/

 typedef std::size_t Index;  ///< index of parameters

/**@} ----------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

/// Constructor, takes the father and the time horizon
/** Constructor of HydroUnitBlock, taking possibly a pointer of its
 * father Block.
 */

 explicit HydroUnitBlock( Block * flbock = nullptr , Index t = 0):
         UnitBlock( flbock ) {}

/*--------------------------------------------------------------------------*/

/// Destructor of HydroUnitBlock

 ~HydroUnitBlock() override = default;

/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */
/// Extends Block::deserialize( netCDF::NcGroup )
/** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
 * the HydroUnitBlock. Besides the mandatory "type" attribute of any :Block,
 * the group must contain all the data required by the base UnitBlock, as
 * described in the comments to UnitBlock::deserialize( netCDF::NcGroup ).
 * In particular, we refer to that description for the crucial dimensions
 * "TimeHorizon", "NumberIntervals" and "ChangeIntervals". The netCDF::NcGroup
 * must then also contain:
 *
 * - The dimension "NumberReservoirs" containing the number of reservoirs(or
 *   nodes) in a cascading system. The dimension is optional, if it is not
 *   provided then it is taken to be == 1 and in this case the cascading
 *   system becomes to a single hydro unit.
 *
 * - The dimension "NumberArcs" containing the set of arcs connecting the
 *   reservoirs in cascading system. If the "NumberReservoirs == 1" this
 *   dimension is not needed to be define since it is a single hydro system
 *   (there is't any arc).
 *
 * - The dimension "NumberGenerators" containing the number of generators
 *   (which could be turbines or pumps) attached to each arc of cascading
 *   hydro unit. The dimension is optional, if it is not provided then it is
 *   taken to be one.
 *
 * - The variable "GeneratorArc", of type int and indexed over the dimension
 *   "NumberGenerators"; the entry GeneratorArc[ i ] tells to which arc of the
 *   cascading system, generator g attached. //TODO Does we need it?
 *
 * - The variable "MinVolumetric", of type double and indexed over both
 *   dimensions "NumberReservoirs" and "NumberIntervals". This is meant to
 *   represent the matrix MinV[ n , t] which, for each reservoir n at each
 *   time instant t contains the minimum volumetric value of the unit for each
 *   reservoir and corresponding time steps; it must be that the entry
 *   MinV[ n , t] >= 0 for all n and t and MinV[ n , t] <= MaxV[ n , t].
 *   MinVolumetric[ n , f]  is the fixed value of MinV[ n , t] for each
 *   reservoir n and all t in the interval
 *   [ ChangeIntervals[ f - 1 ] , ChangeIntervals[ f ] ], with the assumption
 *   that ChangeIntervals[ - 1 ] = 0.  If
 *   "NumberIntervals" <= 1 or "NumberIntervals" >= "TimeHorizon" then the
 *   mapping clearly does not require "ChangeIntervals", which in fact is not
 *   loaded.
 *
 * - The variable "MaxVolumetric", of type double and indexed over both
 *   dimensions "NumberReservoirs" and "NumberIntervals". This is meant to
 *   represent the matrix MaxV[ n , t] which, for each reservoir n at each
 *   time instant t contains the maximum volumetric value of the unit for each
 *   reservoir and corresponding time steps; it must be that the entry
 *   MaxV[ n , t] >= 0 for all n and t and MinV[ n , t] <= MaxV[ n , t].
 *   MaxVolumetric[ n , f]  is the fixed value of MaxV[ n , t] for each
 *   reservoir n and all t in the interval
 *   [ ChangeIntervals[ f - 1 ] , ChangeIntervals[ f ] ], with the assumption
 *   that ChangeIntervals[ - 1 ] = 0.  If
 *   "NumberIntervals" <= 1 or "NumberIntervals" >= "TimeHorizon" then the
 *   mapping clearly does not require "ChangeIntervals", which in fact is not
 *   loaded.
 *
 * - The variable "MinPower", of type double and indexed over both dimensions
 *   "NumberIntervals" and "NumberGenerators". This is meant to represent the
 *   matrix MinP[ t , g ] which, for each time instant t and each generator g,
 *   contains the minimum power output value of the unit for the corresponding
 *   time steps and generator; it must be that MinP[ t , g ] >= 0 for all t and
 *   g, and MinP[ t ,  g ] <= MaxP[ t ,  g ]. MinPower[ i , g ] is the fixed
 *   value of MinP[ t , g ] for all t in the interval
 *   [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the assumption
 *   that ChangeIntervals[ - 1 ] = 0 and each generator g. If
 *   "NumberIntervals" <= 1 or "NumberIntervals" >= "TimeHorizon" then the
 *   mapping clearly does not require "ChangeIntervals", which in fact is not
 *   loaded.
 *
 * - The variable "MaxPower", of type double and indexed over both dimensions
 *   "NumberIntervals" and "NumberGenerators". This is meant to represent the
 *   matrix MinP[ t , g ] which, for each time instant t and each generator g,
 *   contains the maximum power output value of the unit for the corresponding
 *   time steps and generator; it must be that MaxP[ t , g ] >= 0 for all t
 *   and g, and MinP[ t ,  g ] <= MaxP[ t ,  g ]. MaxPower[ i , g ] is the
 *   fixed value of MaxP[ t , g ] for all t in the interval
 *   [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the assumption
 *   that ChangeIntervals[ - 1 ] = 0 and each generator g. If
 *   "NumberIntervals" <= 1 or "NumberIntervals" >= "TimeHorizon" then the
 *   mapping clearly does not require "ChangeIntervals", which in fact is not
 *   loaded.
 *
 * - The variable "MinFlow", of type double and indexed over both dimensions
 *   "NumberIntervals" and "NumberGenerators". This is meant to represent the
 *   matrix MinF[ t , g ] which, for each time instant t and each generator g,
 *   contains the minimum flow rate value of the unit for the corresponding
 *   time steps and generator; it must be that MinF[ t , g ] <= MaxF[ t , g ].
 *   MinFlow[ i , g ] is the fixed value of MinF[ t , g ] for all t in the
 *   interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the
 *   assumption that ChangeIntervals[ - 1 ] = 0 and each generator g. If
 *   "NumberIntervals" <= 1 or "NumberIntervals" >= "TimeHorizon" then the
 *   mapping clearly does not require "ChangeIntervals", which in fact is not
 *   loaded.
 *
 * - The variable "MaxFlow", of type double and indexed over both dimensions
 *   "NumberIntervals" and "NumberGenerators". This is meant to represent the
 *   matrix MaxF[ t , g ] which, for each time instant t and each generator g,
 *   contains the maximum flow rate value of the unit for the corresponding
 *   time steps and generator; it must be that MinF[ t , g ] <= MaxF[ t , g ].
 *   MaxFlow[ i , g ] is the fixed value of MaxF[ t , g ] for all t in the
 *   interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the
 *   assumption that ChangeIntervals[ - 1 ] = 0 and each generator g. If
 *   "NumberIntervals" <= 1 or "NumberIntervals" >= "TimeHorizon" then the
 *   mapping clearly does not require "ChangeIntervals", which in fact is not
 *   loaded.
 *
 * Note: MinFlow and MaxFlow values can be either positive or negative
 *       (or zero); whenever MinF[ t , g] < 0 and MaxF[ t , g] < 0 for each t
 *       and g, the unit is considered a pump and whenever MinF[ t , g] > 0
 *       and MaxF[ t , p] > 0 the unit is considered a turbine. Any possible
 *       mixed situation can be accounted for by artificially splitting the
 *       unit into “two units”, which should be done at the data processing
 *       stage. When MaxF[ t , g] == 0 it means maximum flow rate of the unit
 *       (which is considered as a pump) is zero, and when MinF[ t , g] == 0
 *       it means minimum flow rate of the unit(which is considered as
 *       turbine) is zero.
 *
 * - The variable "DeltaRampUp", of type double and indexed both over
 *   dimensions "NumberIntervals" and "NumberGenerators". This is meant to
 *   represent the matrix DP[ t , g ] which, for each time instant t and each
 *   generator g, contains the maximum possible increase of the flow rate.
 *   This variable is optional; if it is not provided then it is assumed that
 *   DP[ t , g ] == MxF[ t , g ], i.e., the unit can ramp up by an arbitrary
 *   amount, i.e., there are no ramp-up constraints. DeltaRampUp[ i , g ] is
 *   the fixed value of DP[ t , g ] for all t in the interval
 *   [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the assumption
 *   that ChangeIntervals[ - 1 ] = 0 and each generator g. If
 *   "NumberIntervals" <= 1 or "NumberIntervals" >= "TimeHorizon" then the
 *   mapping clearly does not require "ChangeIntervals", which in fact is not
 *   loaded.
 *
 * - The variable "DeltaRampDown", of type double and indexed both over
 *   dimensions "NumberIntervals" and "NumberGenerators". This is meant to
 *   represent the matrix DM[ t , g ] which, for each time instant t and each
 *   generator g, contains the minimum possible increase of the flow rate. This
 *   variable is optional; if it is not provided then it is assumed that
 *   DM[ t , g ] == MaxF[ t , g ], i.e., the unit can ramp down by an
 *   arbitrary amount, i.e., there are no ramp-down constraints.
 *   DeltaRampDown[ i , g ] the fixed value of DP[ t , g ] for all t in the
 *   interval is [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the
 *   assumption that ChangeIntervals[ - 1 ] = 0 and each generator g. If
 *   "NumberIntervals" <= 1 or "NumberIntervals" >= "TimeHorizon" then the
 *   mapping clearly does not require "ChangeIntervals", which in fact is not
 *   loaded.
 *
 * - The variable "PrimaryRho", of type double and indexed both over the
 *   dimensions "NumberIntervals" and "NumberGenerators". This is meant to
 *   represent the matrix PR[ t , g ] which, for each time instant t and
 *   generator g, contains the maximum possible fraction of active power that
 *   can be used as primary reserve. PrimaryRho[ i , g ] is the fixed value of
 *   PR[ t , g ] for all t in the interval
 *   [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the assumption
 *   that ChangeIntervals[ - 1 ] = 0 and all g. If "NumberIntervals" <= 1 or
 *   "NumberIntervals" >= "TimeHorizon" then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded. This variable is
 *   optional, when it's not present it will not capable of producing any
 *   primary reserve, which correspond to PR[ t , g ] == 0 for all t and g.
 *
 * - The variable "SecondaryRho", of type double and indexed both over the
 *   dimensions "NumberIntervals" and "NumberGenerators". This is meant to
 *   represent the matrix SR[ t , g ] which, for each time instant t and
 *   generator g, contains the maximum possible fraction of active power that
 *   can be used as secondary reserve. SecondaryRho[ i , g ] is the fixed
 *   value of SR[ t , g ] for all t in the interval
 *   [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the assumption
 *   that ChangeIntervals[ - 1 ] = 0 and all g. If "NumberIntervals" <= 1 or
 *   "NumberIntervals" >= "TimeHorizon" then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded. This variable is
 *   optional, when it's not present it will not capable of producing any
 *   secondary reserve, which correspond to SR[ t , g ] == 0 for all t and g.
 *
 * - The variable "PowerFlowRho", of type double and indexed both over the
 *   dimensions "NumberIntervals" and "NumberGenerators". This is meant to
 *   represent the matrix PFR[ t , g ] which, for each time instant t and
 *   generator g, contains the exact fraction of active power that can be used
 *   as flow rate. PowerFlowRho[ i , g ] is the fixed value of PFR[ t , g ]
 *   for all t in the interval
 *   [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with the assumption
 *   that ChangeIntervals[ - 1 ] = 0 and all g. If "NumberIntervals" <= 1 or
 *   "NumberIntervals" >= "TimeHorizon" then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded. This variable is
 *   optional, when it's not present it will not capable of producing any
 *   flow rate, which correspond to PFR[ t , g ] == 0 for all t and g.
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
/// generate the static constraint of the HydroUnit
/** Method that generates the static constraint of the HydroUnitBlock.
 * These are the:
 * //TODO I should put all the mathematical constraints here
 *
*/
 void generate_abstract_constraints( Configuration *stcc ) override;
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// generate the objective of the HydroUnitBlock
/** Method that generates the objective of the HydroUnitBlock.
 * //TODO I should put the objective function here
 *
*/
 void generate_objective( Configuration *objc ) override;

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
 * and second dimension number of generators.
 * @{ */

 /// returns the matrix of volumetric variables
 /** The returned boost::multi_array< ColVariable , 2 >, say V, contains the
  * volumetric variables and is indexed over the dimensions time horizon and
  * number of generators. There are two possible cases:
  *
  *  - if V is empty(), then these variables are not defined;
  *
  *  - otherwise, V must have f_time_horizon rows and f_number_generators
  *    columns, and M[ t , g ] is the volumetric variable for time step t of
  *    generator g. */

 const boost::multi_array< ColVariable , 2 > & get_volumetric() const {
  return v_volumetric;
 }
 /*--------------------------------------------------------------------------*/

 /// returns the matrix of flow rate variables
 /** The returned boost::multi_array< ColVariable , 2 >, say F, contains the
  * flow rate variables and is indexed over the dimensions time horizon and
  * number of generators. There are two possible cases:
  *
  *  - if V is empty(), then these variables are not defined;
  *
  *  - otherwise, F must have f_time_horizon rows and f_number_generators
  *    columns, and M[ t , g ] is the flow rate variable for time step t of
  *    generator g. */

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

 /// The number of turbines of the problem
 Index f_number_generators;

 /// v_turbine_arc[ i ] tells to which are turbine i belongs
 std::vector< Index > v_turbine_arc;

 /// The matrix of MinVolumetric
 /** Indexed over the dimensions NumberReservoirs and NumberIntervals. */
 boost::multi_array< double, 2 > v_minimum_volumetric;

 /// The matrix of MaxVolumetric
 /** Indexed over the dimensions NumberReservoirs and NumberIntervals. */
 boost::multi_array< double, 2 > v_maximum_volumetric;

 /// The matrix of MinPower
 /** Indexed over the dimensions NumberIntervals and NumberGenerators. */
 boost::multi_array< double, 2 > v_minimum_power;

 /// The matrix of MaxPower
 /** Indexed over the dimensions NumberIntervals and NumberGenerators. */
 boost::multi_array< double, 2 > v_maximum_power;

 /// The matrix of MinFlow
 /** Indexed over the dimensions NumberIntervals and NumberGenerators. */
 boost::multi_array< double, 2 > v_minimum_flow;

 /// The matrix of DeltaRampUp
 /** Indexed over the dimensions NumberIntervals and NumberGenerators. */
 boost::multi_array< double, 2 > v_delta_ramp_up;

 /// The matrix of DeltaRampDown
 /** Indexed over the dimensions NumberIntervals and NumberGenerators. */
 boost::multi_array< double, 2 > v_delta_ramp_down;

 /// The matrix of PrimaryRho
 /** Indexed over the dimensions NumberIntervals and NumberGenerators. */
 boost::multi_array< double, 2 > v_primary_rho;

 /// The matrix of SecondaryRho
 /** Indexed over the dimensions NumberIntervals and NumberGenerators. */
 boost::multi_array< double, 2 > v_secondary_rho;

 /// The matrix of PowerFlowRho
 /** Indexed over the dimensions NumberIntervals and NumberGenerators. */
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
