/*--------------------------------------------------------------------------*/
/*----------------------- File HydroSystemUnitBlock.h ----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 *
 * Header file for the class HydroSystemUnitBlock, which derives from the
 * Block, in order to define a base class for any possible "hydro unit" plus
 * the linking PolyhedralFunctionBlock to describe the future value of water
 * function in a UCBlock.
 *
 * \version 0.11
 *
 * \date 19 - 05 - 2020
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
#include "HydroUnitBlock.h"
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
 * - The number of HydroUnitBlock in the problem;
 *
 * - A set of hydro units, represented by derived classes of the base class
 *   HydroUnitBlock;
 *
 * - Possibly a PolyhedralFunctionBlock as sub-Block.
 *
 * The first sub-Block of this HydroSystemUnitBlock are the HydroUnitBlock. If
 * this HydroSystemUnitBlock also has a PolyhedralFunctionBlock, then the
 * PolyhedralFunctionBlock is the last sub-Block of this HydroSystemUnitBlock.
 */

class HydroSystemUnitBlock : public UnitBlock {

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

 explicit HydroSystemUnitBlock( Block * father_block = nullptr ) :
  UnitBlock( father_block ) {}

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
 *   (HydroUnitBlock) in the problem.
 *
 * - The groups "HydroUnitBlock_0", "HydroUnitBlock_1", ... ,
 *   "HydroUnitBlock_(n-1)", with n == NumberHydroUnits, containing each one
 *   HydroUnitBlock.
 *
 * - The group "PolyhedralFunctionBlock" which contains a
 *   PolyhedralFunctionBlock, whose PolyhedralFunction represents the
 *   future value of the water (a.k.a. "Bellman values") left at the end of
 *   the time horizon in all the reservoirs of all the HydroUnitBlock of the
 *   HydroSystemUnitBlock.
 *
 * The future value of water function is represented by the single
 * PolyhedralFunction which lives inside the PolyhedralFunctionBlock. The
 * vector of "active" variable of PolyhedralFunction is therefore in a
 * one-to-one correspondence with the set of ColVariable in the HydroUnitBlock
 * that represent the amount of water left in each reservoir at the end of
 * the time horizon. Thus, it is necessary to specify the order of the
 * active ColVariable of the PolyhedralFunction. Let us denote by X[ 0 ],
 * X[ 1 ], ... , X[ R - 1 ] the vector of active ColVariable (i.e.,
 * X[ i ] is the one returned by get_active_var( i ) and R =
 * get_num_active_var()). Since each HydroUnitBlock can have more than one
 * reservoir (cf. HydroUnitBlock::get_number_reservoirs()), R is just the
 * total number of reservoir, which is computed by just calling
 * get_number_reservoirs() on each of the HydroUnitBlock and summing all the
 * results. Clearly, R >= NumberHydroUnits. Some of the HydroUnitBlock may
 * have just one reservoir; if this happens for all the hydro unit blocks
 * (but this is not likely), then R == NumberHydroUnits. In this case the
 * mapping is obvious: X[ i ] is the ColVariable that represent the amount of
 * water left in the only reservoir of HydroUnitBlock_i at the end of the
 * time horizon. When, instead, R > NumberHydroUnits, a mapping must be
 * defined. The mapping is the obvious one: HydroUnitBlock have an ordering
 * n = 0, 1, ..., NumberHydroUnits - 1  (cf. the groups "HydroUnitBlock_0",
 * "HydroUnitBlock_1", ... above), and the reservoirs into each
 * HydroUnitBlock also have a natural ordering, Thus, in general the mapping
 * is:
 *
 *   X[ 0 ] = ColVariable representing the amount of water left in the first
 *            reservoir of HydroUnitBlock_0 at the end of the time horizon
 *
 *   X[ 1 ] = ColVariable representing the amount of water left in the second
 *            reservoir of HydroUnitBlock_0 at the end of the time horizon
 *
 *     ...
 *
 *   X[ k ] = ColVariable representing the amount of water left in the k-th
 *            reservoir of HydroUnitBlock_0 at the end of the time horizon,
 *            with k = HydroUnitBlock_0->get_number_reservoirs()
 *
 *   X[ k + 1 ] = ColVariable representing the amount of water left in the
 *                first reservoir of HydroUnitBlock_1 at the end of the time
 *                horizon
 *
 *   X[ k + 2 ] = ColVariable representing the amount of water left in the
 *                second reservoir of HydroUnitBlock_1 at the end of the time
 *                horizon
 *     ...
 *
 * This must be the format of the data (linear inequalities) that define the
 * PolyhedralFunction: the i-th entry of each vector is related to the
 * future value of the water stored in the reservoir identified by the
 * above mapping. See PolyhedralFunction::deserialize() for details about
 * how the data must be stored in the PolyhedralFunctionBlock group. */

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

 /// Returns the i-th HydroUnitBlock
 HydroUnitBlock * get_hydro_unit_block( Index i ) const;

/*--------------------------------------------------------------------------*/

 virtual Index get_number_generators( void ) const override {
  Index number_generators = 0;
  for( auto sub_block : get_nested_Blocks() ) {
   if( auto unit_block = dynamic_cast< UnitBlock * >( sub_block ) )
    number_generators += unit_block->get_number_generators();
  }
  return number_generators;
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

 /// Deserialize the sub-Blocks of HydroSystemUnitBlock
 void deserialize_sub_blocks( const netCDF::NcGroup & group );

/*--------------------------------------------------------------------------*/

 /// Deserialize the sub-Blocks of HydroSystemUnitBlock that have the given
 /// prefix name
 void deserialize_sub_blocks( const netCDF::NcGroup & group,
                              const std::string & sub_group_name_prefix,
                              Index num_sub_blocks );

/*--------------------------------------------------------------------------*/

 /// Deserialize the PolyhedralFunctionBlock
 void deserialize_polyhedral_function_block
 ( const netCDF::NcGroup & group , const std::string & sub_group_name );

/*--------------------------------------------------------------------------*/

 /// Compute the total number of reservoirs
 Index get_total_number_reservoirs() const;

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
