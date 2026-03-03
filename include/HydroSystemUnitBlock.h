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

/** @} ---------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 * @{ */

 /// constructor, takes the father
 /** Constructor of HydroSystemUnitBlock, taking possibly a pointer of its
  * father Block. */

 explicit HydroSystemUnitBlock( Block * father_block = nullptr )
  : UnitBlock( father_block ) , f_number_hydro_units( 0 ) {}

/*--------------------------------------------------------------------------*/
 /// destructor of HydroSystemUnitBlock

 virtual ~HydroSystemUnitBlock() override;

/** @} ---------------------------------------------------------------------*/
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
  * - Objective function: the objective function of the HydroSystemUnitBlock
  *   is "empty" (a FRealObjective with a LinearFunction inside with no
  *   active variables). */

 void generate_objective( Configuration * objc = nullptr ) override;

/**@} ----------------------------------------------------------------------*/
/*-------- METHODS FOR READING THE DATA OF THE HydroSystemUnitBlock --------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the HydroSystemUnitBlock
 * @{ */

 /// returns the number of hydro units of the problem
 Index get_number_hydro_units( void ) const {
  return( f_number_hydro_units );
  }

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
  if( generator >= v_gen_map.size() )
   return( nullptr );
  return( HUB( v_gen_map[ generator ].first )->get_active_power(
					   v_gen_map[ generator ].second ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of reactive power variables of each HydroUnitBlock

 ColVariable * get_reactive_power( Index generator ) override {
  if( ( ! f_reactive_power ) || ( generator >= v_gen_map.size() ) )
   return( nullptr );
  return( HUB( v_gen_map[ generator ].first )->get_reactive_power(
					   v_gen_map[ generator ].second ) );
 } 

/*--------------------------------------------------------------------------*/
 /// returns the vector of primary reserve variables of each HydroUnitBlock

 ColVariable * get_primary_spinning_reserve( Index generator ) override {
  if( generator >= v_gen_map.size() )
   return( nullptr );
  return( HUB( v_gen_map[ generator ].first )->get_primary_spinning_reserve(
					   v_gen_map[ generator ].second ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of secondary reserve variables of each HydroUnitBlock

 ColVariable * get_secondary_spinning_reserve( Index generator ) override {
  if( generator >= v_gen_map.size() )
   return( nullptr );
  return( HUB( v_gen_map[ generator ].first
	       )->get_secondary_spinning_reserve(
					   v_gen_map[ generator ].second ) );
 }

/*--------------------------------------------------------------------------*/

 Index get_number_generators( void ) const override {
  return( v_gen_map.size() );
  }

/*--------------------------------------------------------------------------*/

 const double * get_inertia_power( Index generator ) const override {
  if( generator >= v_gen_map.size() )
   return( nullptr );
  return( HUB( v_gen_map[ generator ].first )->get_inertia_power(
					   v_gen_map[ generator ].second ) );
  }

/*--------------------------------------------------------------------------*/

 double get_min_power( Index t , Index generator = 0 ) const override {
  if( generator >= v_gen_map.size() )
   return( nullptr );
  return( HUB( v_gen_map[ generator ].first )->get_min_power( t , 
					   v_gen_map[ generator ].second ) );
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

/*--------------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/
 /// returns a Solution storing the current one of this HydroSystemUnitBlock
 /** This method must construct and return a (pointer to a) Solution object
  * representing the current "solution state" of this HydroSystemUnitBlock.
  * This is a HydroSystemUnitBlockSolution extending UnitBlockSolution with
  * the specific extra solution information of HydroUnitBlock.
  *
  * This Solution object is basically a UnitBlockSolution containing in
  * addition the collection of HydroUnitBlockSolution, one for each of the
  * inner HydroUnitBlock of this HydroSystemUnitBlock. However, since the
  * "root" UnitBlockSolution already contains the active power and other
  * variables, the "inner" HydroUnitBlockSolution are configured not to
  * store the same information.
  *
  * The parameter for deciding which kind of Solution must be returned is a
  * single int value, coded bitwise:
  *
  * - the first seven bits (bit 0 to bit 6) are passed to the "inner"
  *   HydroUnitBlockSolution with the first four bits masked (set to 0)
  *   so as to avoid duplicating information
  *
  * This value is to be found as:
  *
  * - if solc is not nullptr and it is a SimpleConfiguration< int >, then it
  *   is solc->f_value;
  *
  * - otherwise, if f_BlockConfig is not nullptr,
  *   f_BlockConfig->f_solution_Configuration is not nullptr and it is a
  *   SimpleConfiguration< int >, then it is
  *   f_BlockConfig->f_solution_Configuration->f_value;
  *
  * - otherwise, it is 63 (save everything). */

 Solution * get_Solution( Configuration * solc = nullptr ,
			  bool emptys = true ) override;

/*--------------------------------------------------------------------------*/
 /// return the "appropriate" [HydroSystem]UnitBlockSolution

 UnitBlockSolution * new_Solution( void ) const override;

/** @} ---------------------------------------------------------------------*/
/*--------------- METHODS FOR SAVING THE HydroSystemUnitBlock --------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for printing & saving the HydroSystemUnitBlock
 * @{ */

/// extends Block::serialize( netCDF::NcGroup )
/** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
 * HydroSystemUnitBlock. See HydroSystemUnitBlock::deserialize( netCDF::
 * NcGroup ) for details of the format of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/** @} ---------------------------------------------------------------------*/
/*----------- METHODS FOR MODIFYING THE HydroSystemUnitBlock ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the HydroSystemUnitBlock
 * @{ */

 /// sets which reserve variables should be there
 /** Basically just calls the method in all the sub-Block. */

 void set_reserve_vars( unsigned char what = 0 ) override {
  UnitBlock:set_reserve_vars( what );
  for( auto * b : v_Block )
   if( auto ub = dynamic_cast< HydroUnitBlock * >( b ) )
    ub->set_reserve_vars( what );
  }

/*--------------------------------------------------------------------------*/
 /// sets whether or not reactive power variables should be there
 /** Basically just calls the method in all the sub-Block. */

 void set_reactive_power( bool reactive = false ) override {
  UnitBlock:set_reactive_power( reactive );
  for( auto * b : v_Block )
   if( auto ub = dynamic_cast< HydroUnitBlock * >( b ) )
    ub->set_reactive_power( reactive );
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

 Index f_number_hydro_units;  ///< the number of hydro units of the problem

/*-------------------------------- variables -------------------------------*/

/*------------------------------- constraints ------------------------------*/

 FRealObjective objective;  ///< the objective function

/*----------------------------- data structures ----------------------------*/

 std::vector< std::pair< Index , Index > > v_gen_map;
 ///< the map between total generators and those in the inner HydroUnitBlock
 /**< v_gen_map[ i ].first  = index of sub-Block in which generator i is
  *   v_gen_map[ i ].second = index of the generator in sub-Block
  *                           v_gen_map[ i ].first to which generator i
  *                           corresponds */

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

 void deserialize_polyhedral_function_block( const netCDF::NcGroup & group ,
					const std::string & sub_group_name );

/*--------------------------------------------------------------------------*/
 /// compute the total number of reservoirs

 Index get_total_number_reservoirs( void ) const;

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 };  // end( class( HydroSystemUnitBlock ) )

/*--------------------------------------------------------------------------*/
/*------------------ CLASS HydroSystemUnitBlockSolution --------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// a [UnitBlock]Solution of a HydroSystemUnitBlock
/** The HydroSystemUnitBlockSolution class derives from UnitBlockSolution and
 * adds the "standard" information stored in there (active power, possibly
 * commitment and primary/secondary reserve) the other information that is
 * typical of the HydroSystemUnitBlock, i.e.,
 *
 * - one HydroUnitBlockSolution for each of the "inner" HydroUnitBlock of
 *   the HydroSystemUnitBlock
 *
 * Note, however, that since the "root" UnitBlockSolution contains (if so
 * required) all the "standard" information, the same is not duplicated in
 * the "inner" HydroUnitBlockSolution */

class HydroSystemUnitBlockSolution : public UnitBlockSolution
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*------------------------------- FRIENDS ----------------------------------*/

 friend HydroSystemUnitBlock;  ///< make HydroSystemUnitBlock friend

/*------ CONSTRUCTING AND DESTRUCTING HydroSystemUnitBlockSolution ---------*/

 /// constructor, it has nothing to do
 explicit HydroSystemUnitBlockSolution( void ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void deserialize( const netCDF::NcGroup & group ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 ~HydroSystemUnitBlockSolution() override = default;
 ///< destructor: it is virtual, and empty

/*--- METHODS DESCRIBING THE BEHAVIOR OF A HydroSystemUnitBlockSolution ---*/

 void read( const Block * block ) override final;

 void write( Block * block ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// serialize a HydroSystemUnitBlockSolution into a netCDF::NcGroup
 /** Serialize a HydroSystemUnitBlockSolution into a netCDF::NcGroup. The
  * format is the one of UnitBlockSolution
  * [cf. UnitBlockSolution::serialize()], plus:
  *
  * - The dimension "NumberHydroUnits" containing the number of hydro units
  *   (HydroUnitBlock) in the problem and therefore the number of "inner"
  *   HydrUnitBlockSolution
  *
  * - The groups "HydroUnitBlockSolution_0", "HydroUnitBlockSolution_1", ...,
  *   "HydroUnitBlockSolution_(n-1)", with n == NumberHydroUnits, containing
  *   each one HydroUnitBlock.
  *
  * Note, however, that since the "root" UnitBlockSolution contains (if so
  * required) all the "standard" information, the same is not duplicated in
  * the "inner" HydroUnitBlockSolution */

 void serialize( netCDF::NcGroup & group ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 HydroSystemUnitBlockSolution * scale( double factor ) const override;

 void sum( const Solution * solution , double multiplier ) override;

 HydroSystemUnitBlockSolution * clone( bool empty = false ) const override;

/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/

 void print( std::ostream &output ) const override {
  output << "HydroSystemUnitBlockSolution [" << this << "]: " << std::endl;
  }

/*---------------------- PRIVATE PART OF THE CLASS -------------------------*/

 private:

/*---------------------------- PRIVATE FIELDS ------------------------------*/

 std::vector< HydroUnitBlockSolution * > v_innerSol;
 ///< v_innerSol[ i ] = HydroUnitBlockSolution of inner Block i

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class( HydroSystemUnitBlockSolution ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* __HydroSystemUnitBlock */

/*--------------------------------------------------------------------------*/
/*--------------------- End File HydroSystemUnitBlock.h --------------------*/
/*--------------------------------------------------------------------------*/
