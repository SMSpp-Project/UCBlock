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

class HydroUnitBlock;     // forward declaration of HydroUnitBlock

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
 *   (HydroUnitBlock) in the problem. The dimension is optional: if it is not
 *   provided then it is taken to be 0, which means that there is no hydro
 *   unit block in the problem.
 *
 * - The groups "HydroUnitBlock_0", "HydroUnitBlock_1", ... ,
 *   "HydroUnitBlock_(n-1)", "PolyhedralFunctionBlock_n" with
 *   n == NumberHydroUnits, containing each one HydroUnitBlock, and the last
 *   one is corresponding to the polyhedral function block of the UCBlock. When
 *   NumberHydroUnits == 0, these groups need not be there since they are not
 *   read. If, instead, NumberHydroUnits > 0, it is an error if the
 *   corresponding groups are not there.
 *
 * - The groups "BellmanValue_0", "BellmanValue_1", ... ,"BellmanValue_R" with
 *   R == NumberHydroUnits, containing the future value of water (Volumetric
 *   variable) of each reservoir in each HydroUnitBlock. Since each
 *   HydroUnitBlock can have more than one reservoir(cf. HydroUnitBlock::
 *   get_number_reservoirs()), a value that is useful in the following is the
 *   total number of reservoirs. We will refer to such number as
 *   "TotalNumberReservoirs", which is computed by just calling
 *   get_number_reservoirs() on each of the HydroUnitBlock and summing all the
 *   results. Clearly, TotalNumberReservoirs >= NumberHydroUnits. Some of the
 *   HydroUnitBlock may have just one reservoir; if this happens for all the
 *   hydro unit blocks (but this is not likely), then TotalNumberReservoirs ==
 *   NumberHydroUnits. It is then useful to be able to assign a unique index
 *   h = 0, 1, ..., NumberHydroUnits - 1 to each of the hydro unit block in
 *   the UCBlock. When TotalNumberReservoirs == NumberHydroUnits the index is
 *   the same as b = 0, 1, ..., NumberHydroUnits - 1 (there is a one-to-one
 *   correspondence between HydroUnitBlock and reservoir, but this is not
 *   likely to happen). When, instead, TotalNumberReservoirs >
 *   NumberHydroUnits, a mapping must be defined. The mapping is the obvious
 *   one: HydroUnitBlock have an ordering b = 0, 1, ..., NumberHydroUnits - 1
 *   (cf. the groups "HydroUnitBlock_0", "HydroUnitBlock_1", ... above), and
 *   the number of reservoirs into each HydroUnitBlock also have some natural
 *   ordering, Thus, in general the mapping is:
 *
 *     bellman value 0 = first reservoir of HydroUnitBlock_0
 *     bellman value 1 = second reservoir of HydroUnitBlock_0
 *     ...
 *     bellman value k = k-th reservoir of HydroUnitBlock_0
 *                    k = HydroUnitBlock_0->get_number_reservoirs()
 *     bellman value k + 1 = first reservoir of HydroUnitBlock_1
 *     bellman value k + 2 = second reservoir of HydroUnitBlock_1
 *     ...
 *
 *   which of course boils down to "h = b" when each HydroUnitBlock has
 *   exactly one reservoir (but this is not assumed to happen).
 */
 void deserialize( netCDF::NcGroup & group ) override;

/**@} ----------------------------------------------------------------------*/
/*-------- METHODS FOR READING THE DATA OF THE HydroSystemUnitBlock --------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the HydroSystemUnitBlock
 *
 * @{ */

 /// Returns the number of hydro units of the problem
 Index get_number_hydro_units() const { return f_number_hydro_units; }

/*--------------------------------------------------------------------------*/
/// Returns the vector of (pointers to) HydroUnitBlock elements.
/** The vector of hydro units in the problem. There are two possible
 * cases:
 *
 * - if the vector is empty, then the there is no hydro unit block;
 *
 * - otherwise the vector must have the size of the number of hydro unit
 *   blocks plus one, and the h-th entry gives the corresponding hydro unit
 *   block h and the last element is the PolyhedralFunctionBlock. */

 const std::vector< HydroUnitBlock * > & get_hydro_unit_blocks() const {
  return v_hydro_unit_blocks;
 }

/*--------------------------------------------------------------------------*/
/// Returns the vector of (pointers to) HydroUnitBlock elements.
/** The vector of bellman values in the problem. There are three possible
 * cases:
 *
 * - if the vector is empty, then the there is no hydro unit block;
 *
 * - otherwise the vector must have the size of the total number of reservoirs
 *   and the h-th entry gives the corresponding bellman value of reservoir h.
 *   */

 const std::vector< HydroUnitBlock * > & get_bellman_values() const {
  return v_bellman_values;
 }
/**@} ----------------------------------------------------------------------*/
/*--------------- METHODS FOR SAVING THE HydroSystemUnitBlock --------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the HydroSystemUnitBlock
 *  @{ */

/// extends Block::serialize( netCDF::NcGroup )
/** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
 *  HydroSystemUnitBlock. See HydroSystemUnitBlock::deserialize( netCDF::
 *  NcGroup ) for details of the format of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

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

/*--------------------------------data--------------------------------------*/
 /// The number of hydro units of the problem
 Index f_number_hydro_units;

 /// The set of HydroUnitBlock
 std::vector< HydroUnitBlock * > v_hydro_unit_blocks;

 /// The set of BellmanValues
 std::vector< HydroUnitBlock * > v_bellman_values;
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
