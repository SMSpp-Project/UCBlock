/*--------------------------------------------------------------------------*/
/*------------------------- File NetworkNode.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * The NetworkNode class represents a node in a network. It is
 * characterized by a set of units (UnitBlocks) that belong to this
 * node.
 *
 * \version 0.1
 *
 * \date 04 - 03 - 2019
 *
 * \author Ali Ghezelsoflu \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Rafael Durbano Lobato \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Ali Ghezelsoflu and Rafael Durbano Lobato
 */

/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __NetworkNode
#define __NetworkNode
/* self-identification: #endif at the end of the file */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <vector>

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

///< namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

class UnitBlock;  // forward definition of UnitBlock

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup NetworkNode_CLASSES Classes in NetworkNode.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS NetworkNode -----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// A node in a network
/** The NetworkNode class represents a node in a network. It is
 * characterized by a set of units (UnitBlock) that belong to this
 * node. */

class NetworkNode {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*--------------- METHODS FOR MODIFYING THE NetworkNode --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the NetworkNode
 *  @{ */

  /// adds a UnitBlock to this NetworkNode
  /** Adds a single (pointer to a) UnitBlock to this NetworkNode. */

  void add_unit_block( UnitBlock * block ) {
    v_unit_blocks.push_back( block );
  }

/*@} -----------------------------------------------------------------------*/
/*----------- METHODS FOR READING THE DATA OF THE NetworkNode --------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the NetworkNode
    @{ */

  /// returns the vector of pointers to UnitBlocks

  const std::vector<UnitBlock *> & get_unit_blocks( void ) {
    return v_unit_blocks;
  }

/*@} -----------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED FIELDS  ----------------------------*/
/*--------------------------------------------------------------------------*/

  std::vector<UnitBlock *> v_unit_blocks;

};   // end( class( NetworkNode ) )

/*@}  end( group( NetworkNode_CLASSES ) ) */

} /* namespace SMSpp_di_unipi_it */

#endif /* NetworkNode.h included */

/*--------------------------------------------------------------------------*/
/*------------------------ End File NetworkNode.h --------------------------*/
/*--------------------------------------------------------------------------*/
