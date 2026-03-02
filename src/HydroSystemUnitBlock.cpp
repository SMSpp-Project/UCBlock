/*--------------------------------------------------------------------------*/
/*---------------------- File HydroSystemUnitBlock.cpp ---------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the HydroSystemUnitBlock class.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Ali Ghezelsoflu \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Ali Ghezelsoflu
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "HydroSystemUnitBlock.h"

#include "FRealObjective.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register HydroSystemUnitBlock to the Block factory
SMSpp_insert_in_factory_cpp_1( HydroSystemUnitBlock );

// register HydroSystemUnitBlockSolution to the Solution factory
SMSpp_insert_in_factory_cpp_0( HydroSystemUnitBlockSolution );

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS OF HydroSystemUnitBlock --------------------*/
/*--------------------------------------------------------------------------*/

HydroSystemUnitBlock::~HydroSystemUnitBlock()
{
 for( auto block : v_Block )
  delete( block );
 v_Block.clear();

 objective.clear();
}

/*--------------------------------------------------------------------------*/

HydroUnitBlock * HydroSystemUnitBlock::get_hydro_unit_block( Index i ) const
{
 return( dynamic_cast< HydroUnitBlock * >( v_Block[ i ] ) );
}

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void HydroSystemUnitBlock::deserialize( const netCDF::NcGroup & group )
{
 #ifndef NDEBUG
  static std::vector< std::string > expected_dims = { "TimeHorizon" ,
                                                      "NumberIntervals" ,
                                                      "NumberHydroUnits" };
  check_dimensions( group , expected_dims , std::cerr );
 #endif

 deserialize_time_horizon( group );
 deserialize_dim( group , "NumberHydroUnits" , f_number_hydro_units );
 deserialize_sub_blocks( group );

 Block::deserialize( group );
 }

/*--------------------------------------------------------------------------*/

void HydroSystemUnitBlock::deserialize_sub_blocks(
					      const netCDF::NcGroup & group )
{
 for( auto block : v_Block )
  delete( block );
 v_Block.clear();

 v_Block.reserve( f_number_hydro_units + 1 );

 deserialize_sub_blocks( group , "HydroUnitBlock_" , f_number_hydro_units );
 deserialize_polyhedral_function_block( group , "PolyhedralFunctionBlock" );
 }

/*--------------------------------------------------------------------------*/

Block::Index HydroSystemUnitBlock::get_total_number_reservoirs( void ) const
{
 Index total_number_reservoirs = 0;
 assert( v_Block.size() >= f_number_hydro_units );
 for( Index i = 0 ; i < f_number_hydro_units ; ++i )
  if( auto hydro_unit_block = dynamic_cast< HydroUnitBlock * >( v_Block[ i ]
								) )
   total_number_reservoirs += hydro_unit_block->get_number_reservoirs();
 return( total_number_reservoirs );
 }

/*--------------------------------------------------------------------------*/

void HydroSystemUnitBlock::deserialize_polyhedral_function_block(
                                         const netCDF::NcGroup & group ,
                                         const std::string & sub_group_name )
{
 if( group.isNull() )
  return;

 auto sub_group = group.getGroup( sub_group_name );
 if( sub_group.isNull() )
  return;

 // Create the PolyhedralFunctionBlock
 std::string class_name = "PolyhedralFunctionBlock";
 auto class_name_attribute = sub_group.getAtt( "type" );
 if( ! class_name_attribute.isNull() )
  class_name_attribute.getValues( class_name );

 auto polyhedral_function_block =
  dynamic_cast< PolyhedralFunctionBlock * >( new_Block( class_name , this ) );

 if( ! polyhedral_function_block )
  throw( std::logic_error( "HydroSystemUnitBlock::deserialize: the type "
                           "attribute of group " + sub_group_name +
                           " must be either 'PolyhedralFunctionBlock' or the "
                           "name of a class derived from "
                           "PolyhedralFunctionBlock." ) );

 // Deserialize the PolyhedralFunctionBlock
 polyhedral_function_block->deserialize( sub_group );

 v_Block.push_back( polyhedral_function_block );
 }

/*--------------------------------------------------------------------------*/

void HydroSystemUnitBlock::deserialize_sub_blocks(
				 const netCDF::NcGroup & group ,
				 const std::string & sub_group_name_prefix ,
				 Index num_sub_blocks )
{
 for( Index i = 0 ; i < num_sub_blocks ; ++i ) {
  std::string sub_group_name = sub_group_name_prefix + std::to_string( i );
  auto sub_group = group.getGroup( sub_group_name );
  auto sub_block = new_Block( sub_group , this );

  if( ! sub_block )
   throw( std::invalid_argument( "HydroSystemUnitBlock::deserialize: error "
                                 "when creating Block from group " +
                                 sub_group_name + "." ) );

  v_Block.push_back( sub_block );
  }
 }

/*--------------------------------------------------------------------------*/

void HydroSystemUnitBlock::generate_abstract_variables( Configuration * stvv )
{
 if( variables_generated() )  // variables have already been generated
  return;                     // nothing to do

 for( auto block : v_Block )
  block->generate_abstract_variables();

 // Collect the active Variables of the PolyhedralFunction: these are the
 // variables representing the final volume of each reservoir.

 std::vector< ColVariable * > x;
 x.reserve( get_total_number_reservoirs() );

 for( Index h = 0 ; h < get_number_hydro_units() ; ++h ) {
  auto hydro_unit_block = get_hydro_unit_block( h );
  hydro_unit_block->generate_abstract_variables();

  const auto number_reservoirs = hydro_unit_block->get_number_reservoirs();
  for( Index r = 0 ; r < number_reservoirs ; ++r )
   x.push_back( hydro_unit_block->get_volume( r , f_time_horizon - 1 ) );
  }

 // Set the active Variable of the PolyhedralFunction
 get_polyhedral_function_block()->get_PolyhedralFunction().
  set_variables( std::move( x ) );

 set_variables_generated();

 }  // end( HydroSystemUnitBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void HydroSystemUnitBlock::generate_objective( Configuration * objc )
{
 if( objective_generated() )  // Objective has already been generated
  return;                     // nothing to do

 for( auto block : v_Block )
  block->generate_objective();

 objective.set_function( new LinearFunction() );

 // Set Block objective
 this->set_objective( &objective );

 set_objective_generated();

}  // end( HydroSystemUnitBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/

Solution * HydroSystemUnitBlock::get_Solution( Configuration * csolc ,
					       bool emptys )
{
 Index wsol = 63;
 if( ( ! csolc ) && f_BlockConfig )
  csolc = f_BlockConfig->f_solution_Configuration;

 if( auto config = dynamic_cast< SimpleConfiguration< int > * >( csolc ) )
  wsol = config->f_value;

 // call the method of the base class
 auto * sol = dynamic_cast< HydroSystemUnitBlockSolution * >(
		                UnitBlock::get_Solution( csolc , emptys ) );
 assert( sol );

 // build a SimpleConfiguration< int > containing the value of wsol with the
 // first four bits masked (zeroed)
 SimpleConfiguration< int > iC( wsol & ~15 );

 // build the empty "inner" HydroUnitBlockSolution
 sol->v_innerSol.resize( get_number_hydro_units() , nullptr );
 for( std::size_t i = 0 ; i < sol->v_innerSol.size() ; ++i )
  sol->v_innerSol[ i ] = static_cast< HydroUnitBlockSolution * >(
	    get_hydro_unit_block( Index( i ) )->get_Solution( &iC , true ) );

 if( ! emptys )
  sol->read( this );

 return( sol );
 }

/*--------------------------------------------------------------------------*/

UnitBlockSolution * HydroSystemUnitBlock::new_Solution( void ) const {
 return( new HydroSystemUnitBlockSolution() );
 }

/*--------------------------------------------------------------------------*/
/*--------------- METHODS FOR SAVING THE HydroSystemUnitBlock --------------*/
/*--------------------------------------------------------------------------*/

void HydroSystemUnitBlock::serialize( netCDF::NcGroup & group ) const
{
 Block::serialize( group );

 auto dim_number_hydro_units = group.addDim( "NumberHydroUnits" ,
                                             f_number_hydro_units );
 // Serialize sub-blocks
 for( Index i = 0 ; i < f_number_hydro_units ; ++i ) {
  auto sub_block = get_hydro_unit_block( i );
  auto sub_group = group.addGroup( "HydroUnitBlock_" + std::to_string( i ) );
  sub_block->serialize( sub_group );
 }

 if( v_Block.size() > f_number_hydro_units ) {
  auto sub_group = group.addGroup( "PolyhedralFunctionBlock" );
  v_Block.back()->serialize( sub_group );
  }
 }

/*--------------------------------------------------------------------------*/
/*-------------- METHODS OF HydroSystemUnitBlockSolution -------------------*/
/*--------------------------------------------------------------------------*/

void HydroSystemUnitBlockSolution::deserialize(
					     const netCDF::NcGroup & group )
{
 // call the method of the base class
 UnitBlockSolution::deserialize( group );

 Index n_units;
 deserialize_dim( group , "NumberHydroUnits" , n_units , false );

 v_innerSol.resize( n_units , nullptr );
 for( std::size_t i = 0 ; i < v_innerSol.size() ; ++i ) {
  std::string sub_group_name = "HydroSystemUnitSolution_" +
                                                        std::to_string( i );
  auto sub_group = group.getGroup( sub_group_name );
  auto HSUSi = dynamic_cast< HydroUnitBlockSolution * >(
				       Solution::new_Solution( sub_group ) );
  if( ! HSUSi )
    throw( std::invalid_argument(
		              "HydroSystemUnitBlockSolution::deserialize: " +
			      sub_group_name +
			      " not a valid HydroUnitBlockSolution" ) );

  v_innerSol[ i ] = HSUSi;
  }
 }  // end( HydroSystemUnitBlockSolution::deserialize )

/*--------------------------------------------------------------------------*/

void HydroSystemUnitBlockSolution::read( const Block * block )
{
 auto HSUB = dynamic_cast< const HydroSystemUnitBlock * >( block );
 if( ! HSUB )
  throw( std::invalid_argument( "HydroSystemUnitBlockSolution::read: block"
				" is not a HydroSystemUnitBlock" ) );

 UnitBlockSolution::read( HSUB );  // call the method of the base class

 for( std::size_t i = 0 ; i < v_innerSol.size() ; ++i )
  v_innerSol[ i ]->read( HSUB->get_hydro_unit_block( Index( i ) ) );

 }  // end( HydroSystemUnitBlockSolution::read )

/*--------------------------------------------------------------------------*/

void HydroSystemUnitBlockSolution::write( Block * block )
{
 UnitBlockSolution::write( block );  // call the method of the base class

 auto HSUB = dynamic_cast< const HydroSystemUnitBlock * >( block );
 if( ! HSUB )
  throw( std::invalid_argument( "HydroSystemUnitBlockSolution::write: block"
				" is not a HydroSystemUnitBlock" ) );

 if(  v_innerSol.size() != HSUB->get_number_hydro_units() )
  throw( std::invalid_argument( "HydroSystemUnitBlockSolution::write: "
				"inconsistent number of hydro units" ) );

 for( std::size_t i = 0 ; i < v_innerSol.size() ; ++i )
  v_innerSol[ i ]->write( HSUB->get_hydro_unit_block( Index( i ) ) );

 }  // end( HydroSystemUnitBlockSolution::write )

/*--------------------------------------------------------------------------*/

void HydroSystemUnitBlockSolution::serialize( netCDF::NcGroup & group ) const
{
 UnitBlockSolution::serialize( group );  // call the method of the base class

 auto nh = group.addDim( "NumberHydroUnits" , v_innerSol.size() );

 for( std::size_t i = 0 ; i < v_innerSol.size() ; ++i ) {
  std::string sub_group_name = "HydroSystemUnitSolution_" +
                                                        std::to_string( i );
  auto sub_group = group.addGroup( sub_group_name );
  v_innerSol[ i ]->serialize( sub_group );
  }
 }  // end( HydroSystemUnitBlockSolution::serialize )

/*--------------------------------------------------------------------------*/

HydroSystemUnitBlockSolution * HydroSystemUnitBlockSolution::scale(
						       double factor ) const
{
 auto sol = clone();

 if( factor != 1 )
  for( std::size_t i = 0 ; i < v_innerSol.size() ; ++i )
   v_innerSol[ i ]->scale( factor );

 return( sol );

 }  // end( HydroSystemUnitBlockSolution::scale )

/*--------------------------------------------------------------------------*/

void HydroSystemUnitBlockSolution::sum( const Solution * solution ,
					double multiplier )
{
 // call the method of the base class
 UnitBlockSolution::sum( solution , multiplier );

 auto HSUBS = dynamic_cast< const HydroSystemUnitBlockSolution * >(
								  solution );
 if( ! HSUBS )
  throw( std::invalid_argument( "HydroSystemUnitBlockSolution::sum: solution"
				" not a HydroSystemUnitBlockSolution" ) );

 if( v_innerSol.size() != HSUBS->v_innerSol.size() )
  throw( std::invalid_argument( "HydroSystemUnitBlockSolution::sum: "
				"inconsistent number of hydro units" ) );

 for( std::size_t i = 0 ; i < v_innerSol.size() ; ++i )
  v_innerSol[ i ]->sum( HSUBS->v_innerSol[ i ] , multiplier );

 }  // end( HydroSystemUnitBlockSolution::sum )

/*--------------------------------------------------------------------------*/

HydroSystemUnitBlockSolution * HydroSystemUnitBlockSolution::clone(
							  bool empty ) const
{
 auto * sol = new HydroSystemUnitBlockSolution();

 if( ! empty ) {
  guts_of_clone( sol );
  sol->v_innerSol.resize( v_innerSol.size() );
  for( std::size_t i = 0 ; i < v_innerSol.size() ; ++i )
   (sol->v_innerSol)[ i ] = v_innerSol[ i ]->clone( false );
  }

 return( sol );

 }  // end( HydroSystemUnitBlockSolution::clone )

/*--------------------------------------------------------------------------*/
/*------------------- End File HydroSystemUnitBlock.cpp --------------------*/
/*--------------------------------------------------------------------------*/
