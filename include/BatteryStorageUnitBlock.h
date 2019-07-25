/*--------------------------------------------------------------------------*/
/*------------------- File BatteryStorageUnitBlock.h -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class BatteryStorageUnitBlock, which derives from
 * UnitBlock [see UnitBlock.h], in order to define a "reasonably standard"
 * battery storage unit at Unit Commitment Problem.
 *
 * \version 0.11
 *
 * \date 18 - 07 - 2019
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

#ifndef __BatteryStorageUnitBlock
#define __BatteryStorageUnitBlock
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
/*--------------------- CLASS BatteryStorageUnitBlock ----------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// Implementation of the Block concept for the battery storage unit problem
/** The BatteryStorageUnitBlock class implements the Block concept
 * [see Block.h] for a "reasonably standard" battery storage unit of a Unit
 * Commitment Problem. That is, the class is designed in order to give
 * mathematical formulation to describe the operation of large set of battery
 * storage. //todo To model complex reservoir systems
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
class BatteryStorageUnitBlock : public UnitBlock {

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
/** Constructor of BatteryStorageUnitBlock, taking possibly a pointer of its
 * father Block.
 */

 explicit BatteryStorageUnitBlock( Block * f_block = nullptr , Index t = 0):
         UnitBlock( f_block ) {}

/*--------------------------------------------------------------------------*/

/// Destructor of BatteryStorageUnitBlock

 ~BatteryStorageUnitBlock() override = default;

/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */
/// Extends Block::deserialize( netCDF::NcGroup )
/** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
 * the BatteryStorageUnitBlock. Besides the mandatory "type" attribute of any
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
/// Generate the abstract variables of the BatteryStorageUnitBlock
/** The BatteryStorageUnitBlock class use get_variable() method to access to
 *  each "group" of variable that may create in UnitBlock class which are:
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
 *  BatteryStorageUnitBlock is defined more groups of variables as follow:
 *
 *  -
 *
 *  -
 *
 *  //todo
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
/// Generate the static constraint of the BatteryStorageUnitBlock
/** Method that generates the static constraint of the BatteryStorageUnitBlock.
 * These are the:
 * //TODO I should put all the mathematical constraints here
 *
*/
 void generate_abstract_constraints( Configuration *stcc ) override;
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// Generate the objective of the BatteryStorageUnitBlock
/** Method that generates the objective of the BatteryStorageUnitBlock.
 * //TODO I should put the objective function here
 *
*/
 void generate_objective( Configuration *objc ) override;

/**@} ----------------------------------------------------------------------*/
/*------ METHODS FOR READING THE DATA OF THE BatteryStorageUnitBlock -------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the BatteryStorageUnitBlock
 *
 * These methods allow to read data that must be common to (in principle) all
 * the kind of battery storage units
 * @{ */
//todo


/**@} ----------------------------------------------------------------------*/
/*----- METHODS FOR READING THE Variable OF THE BatteryStorageUnitBlock ----*/
/*--------------------------------------------------------------------------*/

/** @name Reading the Variable of the BatteryStorageUnitBlock
 *
 * These methods allow to read the two groups of Variable that any
 * BatteryStorageUnitBlock in principle has (although some may not):
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
/*-------------- METHODS FOR SAVING THE BatteryStorageUnitBlock-------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the BatteryStorageUnitBlock
 *  @{ */

/// Extends Block::serialize( netCDF::NcGroup )
/** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
 * BatteryStorageUnitBlock. See
 * BatteryStorageUnitBlock::deserialize( netCDF::NcGroup ) for details of the
 * format of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/**@} ----------------------------------------------------------------------*/
/*--------- METHODS FOR INITIALIZING THE BatteryStorageUnitBlock -----------*/
/*--------------------------------------------------------------------------*/

/** @name Handling the data of the BatteryStorageUnitBlock
    @{ */

 void load( std::istream & input ) override {
  throw ( std::logic_error( "BatteryStorageUnitBlock::load() not "
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

};  // end( class( BatteryStorageUnitBlock ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* BatteryStorageUnitBlock.h included */

/*--------------------------------------------------------------------------*/
/*----------------- End File BatteryStorageUnitBlock.h ---------------------*/
/*--------------------------------------------------------------------------*/
