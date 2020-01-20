/*--------------------------------------------------------------------------*/
/*----------------------- File HydroSystemUnitBlock.h ----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 *
 * Header file for the class HydroSystemUnitBlock, which derives from the
 * Block, in order to define a base class for any possible "hydro unit" and
 * plus the linking PolyhedralFunctionBlock to describe the future value of
 * water function in a UCBlock.
 *
 * \version 0.11
 *
 * \date 16 - 01 - 2020
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
 * \copyright &copy; by Antonio Frangioni, Ali Ghezelsoflu
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __HydroSystemUnitBlock
#define __HydroSystemUnitBlock
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Block.h"
#include "PolyhedralFunctionBlock.h"
#include "UnitBlock.h"
/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it {

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS HydroSystemUnitBlock -----------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// Implementation of the Block concept for "a collection of hydro unit" in UC
/** The class HydroSystemUnitBlock, which derives from the Block, defines a
 * base class for any possible "hydro unit" and the linking
 * PolyhedralFunctionBlock that can be attached to a UCBlock to describe the
 * future value of water function. The base HydroSystemUnitBlock class only
 * has very basic information that can characterize almost any different kind
 * of hydro unit:
 *
 * - A set of hydro units, represented by derived classes of the base class
 * HydroUnitBlock.
 *
 * - A set of bellman values, represented by derived classes of the base class
 * HydroUnitBlock.
 * */

class HydroSystemUnitBlock : public Block {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------- PUBLIC TYPES OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
 *
 * HydroSystemUnitBlock defines the following main public type:
 *
 * @{ */

/**@} ----------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

 /// Constructor, takes the father
 /** Constructor of HydroSystemUnitBlock, taking possibly a pointer of its
  * father Block. */

 explicit HydroSystemUnitBlock( Block * father_block = nullptr );

/*--------------------------------------------------------------------------*/
 /// Destructor of HydroSystemUnitBlock: it is virtual, and empty

 ~HydroSystemUnitBlock() override = default;

/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

/// Extends Block::deserialize( netCDF::NcGroup )
/** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
 * the HydroSystemUnitBlock. Besides the mandatory "type" attribute of any
 * :Block, the group should contain the following:
 *
 * - The dimension "NumberHydroUnits" containing the number of hydro units
 *   (HydroUnitBlock) in the problem;
 *
 * - The groups "HydroUnitBlock_0", "HydroUnitBlock_1", ... ,
 *   "HydroUnitBlock_n" with n == NumberUnits - 1, containing each one
 *   HydroUnitBlock corresponding to one or more electricity generating unit
 *   (electrical generator). It is an error if the corresponding groups are
 *   not there. Each HydroUnitBlock can have more than one electrical
 *   generator.
 *
 */
 void deserialize( netCDF::NcGroup & group ) override;

/**@} ----------------------------------------------------------------------*/
/*-------- METHODS FOR READING THE DATA OF THE HydroSystemUnitBlock --------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the HydroSystemUnitBlock
 *
 * @{ */

/**@} ----------------------------------------------------------------------*/
/*------ METHODS FOR READING THE Variable OF THE HydroSystemUnitBlock ------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the Variable of the HydroSystemUnitBlock
 *
 * @{ */

/**@} ----------------------------------------------------------------------*/
/*--------------- METHODS FOR SAVING THE HydroSystemUnitBlock --------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the HydroSystemUnitBlock
 *  @{ */

 /// extends Block::serialize( netCDF::NcGroup )
 /** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
  *  HydroSystemUnitBlock. See HydroSystemUnitBlock::deserialize( netCDF::NcGroup ) for
  *  details of the format of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/**@} ----------------------------------------------------------------------*/
/*------------- METHODS FOR MODIFYING THE HydroSystemUnitBlock -------------*/
/*--------------------------------------------------------------------------*/
 /** @name Methods for modifying the HydroSystemUnitBlock
  *  @{ */

/**@} ----------------------------------------------------------------------*/
/*------------ METHODS FOR INITIALIZING THE HydroSystemUnitBlock -----------*/
/*--------------------------------------------------------------------------*/
 /** @name Handling the data of the HydroSystemUnitBlock
    @{ */

 void load( std::istream & input ) override {
  throw ( std::logic_error( "HydroSystemUnitBlock::load() not implemented "
                            "yet" ) );
 }

/**@} ----------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED METHODS OF THE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
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

};  // end( class( HydroSystemUnitBlock ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* HydroSystemUnitBlock.h included */

/*--------------------------------------------------------------------------*/
/*--------------------- End File HydroSystemUnitBlock.h --------------------*/
/*--------------------------------------------------------------------------*/
