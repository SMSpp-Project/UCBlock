/*--------------------------------------------------------------------------*/
/*---------------------- File HydroSystemUnitBlock.cpp ---------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the HydroSystemUnitBlock class.
 *
 * \version 0.11
 *
 * \date 08 - 07 - 2019
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
 * \copyright &copy by Antonio Frangioni, Ali Ghezelsoflu
 */

/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "UCBlock.h"
#include "HydroSystemUnitBlock.h"

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
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

HydroUnitBlock * HydroSystemUnitBlock::get_hydro_unit_block( Index i ) const {
 return dynamic_cast<HydroUnitBlock *>( v_Block[ i ] );
}

/*--------------------------------------------------------------------------*/

void HydroSystemUnitBlock::deserialize( netCDF::NcGroup & group ) {

 ::deserialize_dim( group, "NumberHydroUnits", f_number_hydro_units, true );

 Block::deserialize( group );
}
/*--------------------------------------------------------------------------*/

void HydroSystemUnitBlock::deserialize_sub_hydro_blocks( const netCDF::NcGroup & group ) {

 for( auto block : v_Block )
  delete block;

 v_Block.clear();

 deserialize_sub_hydro_blocks( group, "HydroUnitBlock_", f_number_hydro_units );

}
/*--------------------------------------------------------------------------*/

void HydroSystemUnitBlock::deserialize_sub_hydro_blocks
        ( const netCDF::NcGroup & group, const std::string & sub_group_name_prefix,
          const int num_sub_blocks ) {

 for( int i = 0; i < num_sub_blocks; ++i ) {

  std::string sub_group_name = sub_group_name_prefix + std::to_string( i );
  auto sub_group = group.getGroup( sub_group_name );

  if( sub_group.isNull() ) {
   throw ( std::invalid_argument( "HydroSystemUnitBlock::deserialize: " +
                                  sub_group_name + " is not present" ) );
  }

  auto class_name_attribute = sub_group.getAtt( "type" );

  if( class_name_attribute.isNull() ) {
   throw ( std::invalid_argument
           ( "HydroSystemUnitBlock::deserialize: type attribute "
             "is not present in group " + sub_group_name ) );
  }

  std::string class_name;
  class_name_attribute.getValues( class_name );
  auto sub_block = new_Block( class_name, this );
  sub_block->deserialize( sub_group );
  v_Block.push_back( sub_block );
 }
}
/*--------------------------------------------------------------------------*/

void HydroSystemUnitBlock::deserialize_polyhedral_function_block
( const netCDF::NcGroup & group ) {

 deserialize_sub_hydro_blocks( group, "PolyhedralFunctionBlock", 1 );

}
/*--------------------------------------------------------------------------*/
/*------------- METHODS FOR MODIFYING THE HydroSystemUnitBlock -------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------- METHODS FOR SAVING THE HydroSystemUnitBlock --------------*/
/*--------------------------------------------------------------------------*/

void HydroSystemUnitBlock::serialize( netCDF::NcGroup & group ) const {

 Block::serialize( group );

 auto dim_number_hydro_units = group.addDim( "NumberHydroUnits", f_number_hydro_units );

 // Serialize sub-blocks

 for( Index i = 0; i < f_number_hydro_units; ++i ) {
  auto sub_block = get_hydro_unit_block( i );
  auto sub_group = group.addGroup( "HydroUnitBlock_" + std::to_string( i ) );
  sub_block->serialize( sub_group );
 }

}

/*--------------------------------------------------------------------------*/
/*------------------- End File HydroSystemUnitBlock.cpp --------------------*/
/*--------------------------------------------------------------------------*/
