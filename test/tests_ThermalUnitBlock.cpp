#include <gtest/gtest.h>

#include "ThermalUnitBlock.h"

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*------------------------------- GLOBALS ----------------------------------*/
/*--------------------------------------------------------------------------*/

struct TestParameters
{
 std::string test_file;
};

/*--------------------------------------------------------------------------*/
/*------------------------- PARAMETRIZED FIXTURE ---------------------------*/
/*--------------------------------------------------------------------------*/

class ThermalUnitBlockTest : public ::testing::TestWithParam< TestParameters >
{

 protected:
 ThermalUnitBlock * block{};

 ThermalUnitBlockTest( void ) = default;

 ~ThermalUnitBlockTest() override = default;

 void SetUp( void ) override {
  std::string filename( GetParam().test_file );
  netCDF::NcFile f( filename , netCDF::NcFile::read );
  ASSERT_FALSE( f.isNull() );

  netCDF::NcGroupAtt gtype = f.getAtt( "SMS++_file_type" );
  ASSERT_FALSE( gtype.isNull() );

  int type;
  gtype.getValues( &type );
  ASSERT_EQ( type , eBlockFile );

  netCDF::NcGroup bg = f.getGroup( "Block_0" );
  ASSERT_FALSE( bg.isNull() );

  block = static_cast< ThermalUnitBlock * >(
   Block::new_Block( "ThermalUnitBlock" ));
  block->deserialize( bg );
 }

 void TearDown( void ) override {
  delete block;
 }

};

/*--------------------------------------------------------------------------*/
/*------------------------ PARAMETRIZED TEST CASES -------------------------*/
/*--------------------------------------------------------------------------*/

TEST_P( ThermalUnitBlockTest , ChangeMaxPowerRange ) {

 std::vector< double > orig( block->get_max_power() );

 std::vector< double > new_values( 10 , 1 );

 block->set_maximum_power( new_values.begin() , Block::Range( 0 , 10 ) );

 auto max_power = block->get_max_power();

 for( int i = 0 ; i < 10 ; ++i ) {
  EXPECT_EQ( max_power[ i ] , 1 );
 }
}

TEST_P( ThermalUnitBlockTest , ChangeMaxPowerSubset ) {

 std::vector< double > orig( block->get_max_power() );

 std::vector< double > new_values( 10 , 1 );

 block->set_maximum_power( new_values.begin() , Block::Range( 0 , 10 ) );

 auto max_power = block->get_max_power();

 for( int i = 0 ; i < 10 ; ++i ) {
  EXPECT_EQ( max_power[ i ] , 1 );
 }
}

/*--------------------------------------------------------------------------*/
/*------------------------- TEST CASE INSTANCES ----------------------------*/
/*--------------------------------------------------------------------------*/

INSTANTIATE_TEST_SUITE_P( ThermalUnitBlockTests ,
                          ThermalUnitBlockTest ,
                          ::testing::Values(
                           TestParameters{
                            "netCDF_files/1UC_Data/24/S1ramp1_24.nc4" }
                          ) );

/*--------------------------------------------------------------------------*/
/*---------------------------------- MAIN ----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( int argc , char ** argv ) {
 ::testing::InitGoogleTest( &argc , argv );
 return( RUN_ALL_TESTS() );
}
