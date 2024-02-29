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
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Ali Ghezelsoflu \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Rafael Durbano Lobato \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Ali Ghezelsoflu,
 *                      Rafael Durbano Lobato
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

#include "FRealObjective.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)

namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS HydroSystemUnitBlock -----------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// implementation of the Block concept for "a collection of hydro unit" in UC
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

class HydroSystemUnitBlock : public UnitBlock
{

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
 * @{ */

 /// constructor, takes the father
 /** Constructor of HydroSystemUnitBlock, taking possibly a pointer of its
  * father Block. */

 explicit HydroSystemUnitBlock( Block * father_block = nullptr )
  : UnitBlock( father_block ) {}

/*--------------------------------------------------------------------------*/
 /// destructor of HydroSystemUnitBlock

 virtual ~HydroSystemUnitBlock() override;

/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 * @{ */

/// extends Block::deserialize( netCDF::NcGroup )
/** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
 * the HydroSystemUnitBlock. Besides the mandatory "type" attribute of any
 * :Block, the group should contain the following:
 *
 * - The dimension "NumberHydroUnits" containing the number of hydro units
 *   (HydroUnitBlock) in the problem.
 *
 * - The groups "HydroUnitBlock_0", "HydroUnitBlock_1", ...,
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
 * X[ 1 ], ..., X[ R - 1 ] the vector of active ColVariable (i.e.,
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

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the objective of the HydroSystemUnitBlock
 /** Method that generates the objective of the HydroSystemUnitBlock.
  *
  * - Objective function: the objective function of the HydroSystemUnitBlock is
  *   "empty" (a FRealObjective with a LinearFunction inside with no active
  *   variables). */

 void generate_objective( Configuration * objc = nullptr ) override;

/**@} ----------------------------------------------------------------------*/
/*-------- METHODS FOR READING THE DATA OF THE HydroSystemUnitBlock --------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the HydroSystemUnitBlock
 * @{ */

 /// returns the number of hydro units of the problem
 Index get_number_hydro_units( void ) const { return( f_number_hydro_units ); }

/*--------------------------------------------------------------------------*/
 /// returns the i-th HydroUnitBlock

 HydroUnitBlock * get_hydro_unit_block( Index i ) const;

/*--------------------------------------------------------------------------*/
 /// returns the PolyhedralFunctionBlock

 PolyhedralFunctionBlock * get_polyhedral_function_block( void ) const {
  assert( ! v_Block.empty() );
  return( static_cast< PolyhedralFunctionBlock * >( v_Block.back() ) );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of active power variables of each HydroUnitBlock

 ColVariable * get_active_power( Index generator ) override {
  auto temp = generator;
  for( auto sub_block : get_nested_Blocks() )
   if( auto unit_block = dynamic_cast< HydroUnitBlock * >( sub_block ) ) {
    if( temp < unit_block->get_number_generators() )
     return( unit_block->get_active_power( temp ) );
    else
     temp = temp - unit_block->get_number_generators();
   }
  return( nullptr );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of primary spinning reserve variables of each HydroUnitBlock

 ColVariable * get_primary_spinning_reserve( Index generator ) override {
  auto temp = generator;
  for( auto sub_block : get_nested_Blocks() )
   if( auto unit_block = dynamic_cast< HydroUnitBlock * >( sub_block ) ) {
    if( temp < unit_block->get_number_generators() )
     return( unit_block->get_primary_spinning_reserve( temp ) );
    else
     temp = temp - unit_block->get_number_generators();
   }
  return( nullptr );
 }

/*--------------------------------------------------------------------------*/
 /// returns the vector of secondary spinning reserve variables of each HydroUnitBlock

 ColVariable * get_secondary_spinning_reserve( Index generator ) override {
  auto temp = generator;
  for( auto sub_block : get_nested_Blocks() )
   if( auto unit_block = dynamic_cast< HydroUnitBlock * >( sub_block ) ) {
    if( temp < unit_block->get_number_generators() )
     return( unit_block->get_secondary_spinning_reserve( temp ) );
    else
     temp = temp - unit_block->get_number_generators();
   }
  return( nullptr );
 }

/*--------------------------------------------------------------------------*/

 virtual Index get_number_generators( void ) const override {
  Index number_generators = 0;
  for( auto sub_block : get_nested_Blocks() )
   if( auto unit_block = dynamic_cast< UnitBlock * >( sub_block ) )
    number_generators += unit_block->get_number_generators();
  return( number_generators );
 }

/*--------------------------------------------------------------------------*/

 const double * get_inertia_power( Index generator ) const override {
  auto temp = generator;
  for( auto sub_block : get_nested_Blocks() )
   if( auto unit_block = dynamic_cast< HydroUnitBlock * >( sub_block ) ) {
    if( temp < unit_block->get_number_generators() )
     return( unit_block->get_inertia_power( temp ) );
    else
     temp = temp - unit_block->get_number_generators();
   }
  return( nullptr );
 }

/*--------------------------------------------------------------------------*/

 double get_min_power( Index t , Index generator = 0 ) const override {
  auto temp = generator;
  for( auto sub_block : get_nested_Blocks() )
   if( auto unit_block = dynamic_cast< HydroUnitBlock * >( sub_block ) ) {
    if( temp < unit_block->get_number_generators() )
     return( unit_block->get_min_power( t, temp ) );
    else
     temp = temp - unit_block->get_number_generators();
   }
  return( 0 );
 }

/*--------------------------------------------------------------------------*/

 double get_max_power( Index t , Index generator = 0 ) const override {
  auto temp = generator;
  for( auto sub_block : get_nested_Blocks() )
   if( auto unit_block = dynamic_cast< HydroUnitBlock * >( sub_block ) ) {
    if( temp < unit_block->get_number_generators() )
     return( unit_block->get_max_power( t, temp ) );
    else
     temp = temp - unit_block->get_number_generators();
   }
  return( 0 );
 }

/**@} ----------------------------------------------------------------------*/
/*--------------- METHODS FOR SAVING THE HydroSystemUnitBlock --------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for printing & saving the HydroSystemUnitBlock
 * @{ */

/// extends Block::serialize( netCDF::NcGroup )
/** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
 * HydroSystemUnitBlock. See HydroSystemUnitBlock::deserialize( netCDF::
 * NcGroup ) for details of the format of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/**@} ----------------------------------------------------------------------*/
/*----------- METHODS FOR MODIFYING THE HydroSystemUnitBlock ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the HydroSystemUnitBlock
 * @{ */

 /// sets reserve vars method
 /** This method can be called *after* that deserialize() and before
  * generate_abstract_variables() and generate_abstract_constraints(). This is
  * called to provide the UCBlock with the reserve variables if it's needed.
  * The input parameter is a bitwise value that allows to specify which unit
  * could have the reserve variables:
  *
  * - 1 the unit could have primary spinning reserve variables
  * - 2 the unit could have secondary spinning reserve variables
  *
  * Note: this method is only to "destroy" the (primary, secondary and inertia)
  * reserve variables; it cannot create them if they are not there. */

 void set_reserve_vars( unsigned char what ) override {
  reserve_vars = what;
  for( auto * b : v_Block )
   if( auto ub = dynamic_cast< HydroUnitBlock * >( b ) )
    ub->set_reserve_vars( what );
 }

/** @} ---------------------------------------------------------------------*/
/*------------ METHODS FOR INITIALIZING THE HydroSystemUnitBlock -----------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the data of the HydroSystemUnitBlock
 * @{ */

 void load( std::istream & input , char frmt = 0 ) override {
  throw( std::logic_error(
   "HydroSystemUnitBlock::load() not implemented  yet" ) );
 }

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED METHODS OF THE CLASS ---------------------*/
/*--------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

/*---------------------------------- data ----------------------------------*/

 /// the number of hydro units of the problem
 Index f_number_hydro_units{};

/*-------------------------------- variables -------------------------------*/



/*------------------------------- constraints ------------------------------*/

 /// the objective function
 FRealObjective objective;

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------- PRIVATE FIELDS OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/



 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/
/*---------------------- PRIVATE METHODS OF THE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

 /// deserialize the sub-Blocks of HydroSystemUnitBlock

 void deserialize_sub_blocks( const netCDF::NcGroup & group );

/*--------------------------------------------------------------------------*/
 /// deserialize the sub-Blocks of HydroSystemUnitBlock that have the given
 /// prefix name

 void deserialize_sub_blocks( const netCDF::NcGroup & group ,
                              const std::string & sub_group_name_prefix ,
                              Index num_sub_blocks );

/*--------------------------------------------------------------------------*/
 /// deserialize the PolyhedralFunctionBlock

 void deserialize_polyhedral_function_block
  ( const netCDF::NcGroup & group , const std::string & sub_group_name );

/*--------------------------------------------------------------------------*/
 /// compute the total number of reservoirs

 Index get_total_number_reservoirs( void ) const;

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
