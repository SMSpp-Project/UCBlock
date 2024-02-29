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

#include "UCBlock.h"

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

 UnitBlock::deserialize_time_horizon( group );
 ::deserialize_dim( group , "NumberHydroUnits" , f_number_hydro_units );
 deserialize_sub_blocks( group );

 Block::deserialize( group );
}

/*--------------------------------------------------------------------------*/

void HydroSystemUnitBlock::deserialize_sub_blocks( const netCDF::NcGroup & group )
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
  if( auto hydro_unit_block = dynamic_cast< HydroUnitBlock * >( v_Block[ i ] ) )
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

 auto polyhedral_function_block = dynamic_cast< PolyhedralFunctionBlock * >
 ( new_Block( class_name , this ) );

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
 const Index num_sub_blocks )
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
/*------------------- End File HydroSystemUnitBlock.cpp --------------------*/
/*--------------------------------------------------------------------------*/
