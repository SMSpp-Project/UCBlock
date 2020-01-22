/*--------------------------------------------------------------------------*/
/*------------------- File PolyhedralFunctionBlock.h -----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class PolyhedralFunctionBlock, which derives from
 * AbstractBlock to define the class of Block who have the specific structure
 * of having a PolyhedralFunction as objective, but otherwise can contain any
 * kind of Variable and Constraint (provided these are handled by the base
 * AbstractBlock class).
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

#ifndef __PolyhedralFunctionBlock
#define __PolyhedralFunctionBlock
                     /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "HydroSystemUnitBlock.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*-------------------- CLASS PolyhedralFunctionBlock -----------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// ...
/**
 *
 * */
class PolyhedralFunctionBlock : public Block {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public Types
    @{ */

/**@} ----------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructing and destructing PolyhedralFunctionBlock
 *  @{ */

 explicit PolyhedralFunctionBlock( Block * f_block = nullptr):
 HydroSystemUnitBlock( f_block ) {}

/*--------------------------------------------------------------------------*/
/// destructor of PolyhedralFunctionBlock

 ~PolyhedralFunctionBlock() override = default;

/**@} ----------------------------------------------------------------------*/
/*------- Methods for reading the data of the PolyhedralFunctionBlock ------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for reading the data of the PolyhedralFunctionBlock
 *  @{ */


/**@} ----------------------------------------------------------------------*/
/*----- METHODS FOR READING THE Variable OF THE PolyhedralFunctionBlock ----*/
/*--------------------------------------------------------------------------*/
/** @name Reading the Variable of the PolyhedralFunctionBlock
 *
 * @{ */

/**@} ----------------------------------------------------------------------*/
/*---------- METHODS FOR MODIFYING THE PolyhedralFunctionBlock -------------*/
/*--------------------------------------------------------------------------*/
 /** @name Methods for modifying the PolyhedralFunctionBlock
  *  @{ */

 /// get PolyhedralFunction
 /**
  *
  */

/*--------------------------------------------------------------------------*/
 /// sets the set of active Variable of the PolyhedralFunction
 /** Sets the set of active Variable of the PolyhedralFunction.
  */
/**@} ----------------------------------------------------------------------*/
/*---------- METHODS FOR INITIALIZING THE PolyhedralFunctionBlock ----------*/
/*--------------------------------------------------------------------------*/
 /** @name Handling the data of the PolyhedralFunctionBlock
    @{ */

 /// load the PolyhedralFunctionBlock out of an istream
 /** Method to deserialize the PolyhedralFunctionBlock out of an istream.
  *
  *     IT IS CURRENTLY NOT IMPLEMENTED
  *
  * but it still have to be defined (throwing exception) to make the class
  * concrete. */

 void load ( std::istream &input ) override {
  throw( std::logic_error(
          "PolyhedralFunctionBlock::load not implemented yet" ) );
 }
/**@} ----------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:
/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED METHODS OF THE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED FIELDS  ----------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE METHODS ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

};  // end( class( PolyhedralFunctionBlock ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* PolyhedralFunctionBlock.h included */

/*--------------------------------------------------------------------------*/
/*------------------- End File PolyhedralFunctionBlock.h -------------------*/
/*--------------------------------------------------------------------------*/
