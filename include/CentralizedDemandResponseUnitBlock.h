/*--------------------------------------------------------------------------*/
/*------------- File CentralizedDemandResponseUnitBlock.h ------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class CentralizedDemandResponseUnitBlock, which derives from
 * UnitBlock [see UnitBlock.h], in order to define a "reasonably standard"
 * Centralized Demand Response unit at Unit Commitment Problem.
 *
 * \version 0.11
 *
 * \date 25 - 07 - 2019
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

#ifndef __CentralizedDemandResponseUnitBlock
#define __CentralizedDemandResponseUnitBlock
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
/*------------- CLASS CentralizedDemandResponseUnitBlock -------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// Implementation of the Block concept for the CentralizedDemandResponseUnitBlock unit problem
/** The CentralizedDemandResponseUnitBlock class implements the Block concept
 * [see Block.h] for a "reasonably standard" Centralized Demand Response unit of a Unit
 * Commitment Problem. That is, the class is designed in order to give
 * mathematical formulation to describe the operation of large set of
 * CentralizedDemandResponseUnitBlock  //todo
 * The technical and physical constraints are mainly divided in ?? different
 * categories:
 * - ??
 * -??
 * -??
 * -??
 */
class CentralizedDemandResponseUnitBlock : public UnitBlock {

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
/** Constructor of CentralizedDemandResponseUnitBlock, taking possibly a pointer of its
 * father Block.
 */

 explicit CentralizedDemandResponseUnitBlock( Block * f_block = nullptr , Index t = 0):
         UnitBlock( f_block ) {}

/*--------------------------------------------------------------------------*/

/// Destructor of CentralizedDemandResponseUnitBlock

 ~CentralizedDemandResponseUnitBlock() override = default;

/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */
/// Extends Block::deserialize( netCDF::NcGroup )
/** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
 * the CentralizedDemandResponseUnitBlock. Besides the mandatory "type" attribute of any
 * :Block, the group must contain all the data required by the base UnitBlock,
 * as described in the comments to UnitBlock::deserialize( netCDF::NcGroup ).
 * In particular, we refer to that description for the crucial dimensions
 * "TimeHorizon", "NumberIntervals" and "ChangeIntervals". The netCDF::NcGroup
 * must then also contain:
 * //TODO Does
 *
 * */

 void deserialize( netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
/// Generate the abstract variables of the CentralizedDemandResponseUnitBlock
/** The CentralizedDemandResponseUnitBlock class use get_variable() method to access to
 *  each "group" of variable that may create in UnitBlock class which are:
 *
 *
 *  //todo
 * */
 void generate_abstract_variables( Configuration *stvv ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// Generate the static constraint of the CentralizedDemandResponseUnitBlock
/** Method that generates the static constraint of the CentralizedDemandResponseUnitBlock.
 * These are the:
 * //TODO I should put all the mathematical constraints here
 *
*/
 void generate_abstract_constraints( Configuration *stcc ) override;
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// Generate the objective of the CentralizedDemandResponseUnitBlock
/** Method that generates the objective of the CentralizedDemandResponseUnitBlock.
 * //TODO I should put the objective function here
 *
*/
 void generate_objective( Configuration *objc ) override;

/**@} ----------------------------------------------------------------------*/
/*- METHODS FOR READING THE DATA OF THE CentralizedDemandResponseUnitBlock -*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the CentralizedDemandResponseUnitBlock
 *
 * These methods allow to read data that must be common to (in principle) all
 * the kind of Centralized Demand Response units
 * @{ */
//todo


/**@} ----------------------------------------------------------------------*/
/*-- METHODS FOR READING THE Variable OF THE CentralizedDemandResponseUnit--*/
/*--------------------------------------------------------------------------*/

/** @name Reading the Variable of the CentralizedDemandResponseUnitBlock
 *
 * These methods allow to read the two groups of Variable that any
 * CentralizedDemandResponseUnitBlock in principle has (although some may
 * not):
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
/*-------- METHODS FOR SAVING THE CentralizedDemandResponseUnitBlock--------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the
 * CentralizedDemandResponseUnitBlock
 *  @{ */

/// Extends Block::serialize( netCDF::NcGroup )
/** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
 * CentralizedDemandResponseUnitBlock. See
 * CentralizedDemandResponseUnitBlock::deserialize( netCDF::NcGroup ) for
 * details of the
 * format of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/**@} ----------------------------------------------------------------------*/
/*----- METHODS FOR INITIALIZING THE CentralizedDemandResponseUnitBlock ----*/
/*--------------------------------------------------------------------------*/

/** @name Handling the data of the CentralizedDemandResponseUnitBlock
    @{ */

 void load( std::istream & input ) override {
  throw ( std::logic_error( "CentralizedDemandResponseUnitBlock::load() not "
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

};  // end( class( CentralizedDemandResponseUnitBlock ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* CentralizedDemandResponseUnitBlock.h included */

/*--------------------------------------------------------------------------*/
/*------------ End File CentralizedDemandResponseUnitBlock.h ---------------*/
/*--------------------------------------------------------------------------*/
