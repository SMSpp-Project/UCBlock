/*--------------------------------------------------------------------------*/
/*----------------- File IntermittentGenerationUnitBlock.h -----------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class IntermittentGenerationUnitBlock, which derives
 * from UnitBlock [see UnitBlock.h], in order to define a "reasonably
 * standard" Intermittent Generation unit at Unit Commitment Problem.
 *
 * \version 0.11
 *
 * \date 30 - 07 - 2019
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

#ifndef __IntermittentGenerationUnitBlock
#define __IntermittentGenerationUnitBlock
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "UnitBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------------ NAMESPACE ---------------------------------*/
/*--------------------------------------------------------------------------*/

/// Namespace for the Structured Modeling System++ (SMS++)

namespace SMSpp_di_unipi_it {

/*--------------------------------------------------------------------------*/
/*---------------- CLASS IntermittentGenerationUnitBlock -------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// implementation of the Block concept for the Intermittent Generation unit
/** The IntermittentGenerationUnitBlock class implements the Block concept
 * [see Block.h] for a "reasonably standard" Intermittent Generation and
 * Distributed generation units in a single unit of the unit commitment
 * Problem. That is, the class is designed in order to give mathematical
 * formulation to describe the operation of large set of Intermittent GenerationUnitBlock  //todo
 * The technical and physical constraints are mainly divided in ?? different
 * categories:
 * - ??
 * -??
 * -??
 * -??
 */
class IntermittentGenerationUnitBlock : public UnitBlock {

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

/// Constructor, takes the father and the time horizon
/** Constructor of IntermittentGenerationUnitBlock, taking possibly a pointer of its
 * father Block.
 */

 explicit IntermittentGenerationUnitBlock( Block * f_block = nullptr , Index t = 0):
         UnitBlock( f_block ) {}

/*--------------------------------------------------------------------------*/

/// Destructor of IntermittentGenerationUnitBlock

 ~IntermittentGenerationUnitBlock() override = default;

/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */
/// Extends Block::deserialize( netCDF::NcGroup )
/** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
 * the IntermittentGenerationUnitBlock. Besides the mandatory "type" attribute of any
 * :Block, the group must contain all the data required by the base UnitBlock,
 * as described in the comments to UnitBlock::deserialize( netCDF::NcGroup ).
 * In particular, we refer to that description for the crucial dimensions
 * "TimeHorizon", "NumberIntervals" and "ChangeIntervals". The netCDF::NcGroup
 * must then also contain:
 *
 * - The variable "MinPower", of type double and either of size 1 or indexed
 *   over the dimension "NumberIntervals". This is meant to represent the
 *   vector MinP[ t ] that, for each time instant t, contains the minimum
 *   active power output value of the unit for the corresponding time step.
 *   If "MinPower" has length 1 then MinP[ t ] contains the same value for all
 *   t. Otherwise, MinPower[ i ] is the fixed value of MinP[ t ] for all t in
 *   the interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with
 *   the assumption that ChangeIntervals[ - 1 ] = 0.  If NumberIntervals <= 1
 *   or NumberIntervals >= TimeHorizon, then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded. Note that it must
 *   be MxP[ t ] >= MnP[ t ] >= 0 for all t, and when MxP[ t ] == MnP[ t ] the
 *   unit cannot be curtailed and cannot provide any reserve.
 *
 * - The variable "MaxPower", of type double and either of size 1 or indexed
 *   over the dimension "NumberIntervals". This is meant to represent the
 *   vector MaxP[ t ] that, for each time instant t, contains the maximum
 *   active power output value of the unit for the corresponding time step.
 *   If "MaxPower" has length 1 then MaxP[ t ] contains the same value for all
 *   t. Otherwise, MaxPower[ i ] is the fixed value of MaxP[ t ] for all t in
 *   the interval [ ChangeIntervals[ i - 1 ] , ChangeIntervals[ i ] ], with
 *   the assumption that ChangeIntervals[ - 1 ] = 0. If NumberIntervals <= 1
 *   or NumberIntervals >= TimeHorizon, then the mapping clearly does not
 *   require "ChangeIntervals", which in fact is not loaded. Note that it must
 *   be MxP[ t ] >= MnP[ t ] >= 0 for all t, and when MxP[ t ] == MnP[ t ] the
 *   unit cannot be curtailed and cannot provide any reserve.
 *
 * - The scalar variable "Gama", of type double and not indexed over
 *   any dimension. This variable is used to take into account an uncertainty
 *   on the maximal potential production. Note that it must be 0 >= Gama >= 1;
 *   when Gama == 0, the unit does not provide any reserve.
 * */

 void deserialize( netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
/// Generate the abstract variables of the IntermittentGenerationUnitBlock
/** The IntermittentGenerationUnitBlockUnitBlock class use get_variable() method to access to
 *  each "group" of variable that may create in UnitBlock class which are:
 *
 *
 *  //todo
 * */
 void generate_abstract_variables( Configuration *stvv ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// Generate the static constraint of the IntermittentGenerationUnitBlock
/** Method that generates the static constraint of the IntermittentGenerationUnitBlock.
 * These are the:
 * //TODO I should put all the mathematical constraints here
 *
*/
 void generate_abstract_constraints( Configuration *stcc ) override;
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// Generate the objective of the IntermittentGenerationUnitBlock
/** Method that generates the objective of the IntermittentGenerationUnitBlock.
 * //TODO I should put the objective function here
 *
*/
 void generate_objective( Configuration *objc ) override;

/**@} ----------------------------------------------------------------------*/
/*-- METHODS FOR READING THE DATA OF THE IntermittentGenerationUnitBlock ---*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the IntermittentGenerationUnitBlock
 *
 * These methods allow to read data that must be common to (in principle) all
 * the kind of IntermittentGeneration units
 * @{ */
//todo


/**@} ----------------------------------------------------------------------*/
/*- METHODS FOR READING THE Variable OF THE IntermittentGenerationUnitBlock-*/
/*--------------------------------------------------------------------------*/

/** @name Reading the Variable of the IntermittentGenerationUnitBlock
 *
 * These methods allow to read the two groups of Variable that any
 * IntermittentGenerationUnitBlock in principle has (although some may not):
 *
 * - ??
 *
 * - ??
 *
 * All these two groups of variables are (if not empty) ???
 * boost::multi_array< ColVariable , 2 > with first dimension time horizon
 * and second dimension number of generators.
 * @{ */

/**@} ----------------------------------------------------------------------*/
/*-------- METHODS FOR SAVING THE IntermittentGenerationUnitBlock-----------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the IntermittentGenerationUnitBlock
 *  @{ */

/// Extends Block::serialize( netCDF::NcGroup )
/** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
 * IntermittentGenerationUnitBlock. See
 * IntermittentGenerationUnitBlock::deserialize( netCDF::NcGroup ) for details of the
 * format of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/**@} ----------------------------------------------------------------------*/
/*------ METHODS FOR INITIALIZING THE IntermittentGenerationUnitBlock ------*/
/*--------------------------------------------------------------------------*/

/** @name Handling the data of the IntermittentGenerationUnitBlock
    @{ */

 void load( std::istream & input ) override {
  throw ( std::logic_error( "IntermittentGenerationUnitBlock::load() not "
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


/*-----------------------------variables------------------------------------*/


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

};  // end( class( IntermittentGenerationUnitBlock ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* IntermittentGenerationUnitBlock.h included */

/*--------------------------------------------------------------------------*/
/*--------- End File IntermittentGenerationUnitBlock.h ---------------------*/
/*--------------------------------------------------------------------------*/
