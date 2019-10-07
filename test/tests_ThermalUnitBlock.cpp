#include <gtest/gtest.h>
#include <UCBlock.h>
#include <ThermalUnitBlock.h>
#include <CPXMILPSolver.h>
// #include <BusNetworkBlock.h>

using namespace SMSpp_di_unipi_it;


struct TestParameters {
 std::string test_file;
};

/*--------------------------------------------------------------------------*/
/*------------------------- PARAMETRIZED FIXTURE ---------------------------*/
/*--------------------------------------------------------------------------*/

class ThermalUnitBlockTest :
 public ::testing::TestWithParam< TestParameters > {

 protected:

 ThermalUnitBlock * block{};

 ThermalUnitBlockTest() = default;

 ~ThermalUnitBlockTest() override = default;

 void SetUp() override {
  std::string filename( GetParam().test_file );
  netCDF::NcFile f( filename, netCDF::NcFile::read );
  ASSERT_FALSE( f.isNull() );

  netCDF::NcGroupAtt gtype = f.getAtt( "SMS++_file_type" );
  ASSERT_FALSE( gtype.isNull() );

  int type;
  gtype.getValues( &type );
  ASSERT_EQ( type, eBlockFile );

  netCDF::NcGroup bg = f.getGroup( "Block_0" );
  ASSERT_FALSE( bg.isNull() );

  block = dynamic_cast<ThermalUnitBlock *>(Block::new_Block( "ThermalUnitBlock" ));
  block->deserialize( bg );

  auto milpsolver = new CPXMILPSolver();
  block->register_Solver( milpsolver );
 }
};

/*--------------------------------------------------------------------------*/
/*------------------------ PARAMETRIZED TEST CASES -------------------------*/
/*--------------------------------------------------------------------------*/

TEST_P( ThermalUnitBlockTest, SimpleSolve ) {
 int tmp = 15;
 SimpleConfiguration< int > myconfig( tmp );

 ASSERT_NO_THROW( {
                   block->generate_abstract_variables( &myconfig );
                   block->generate_abstract_constraints( nullptr );
                   block->generate_objective( nullptr );
                  } );

 auto solver = block->get_registered_solvers().front();
 int status = solver->compute();
 ASSERT_EQ( status, Solver::kOK );
 // auto ub = solver->get_ub();
}

/*--------------------------------------------------------------------------*/
/*------------------------- TEST CASE INSTANCES ----------------------------*/
/*--------------------------------------------------------------------------*/

INSTANTIATE_TEST_CASE_P( ThermalUnitBlockTests,
                         ThermalUnitBlockTest,
                         ::testing::Values(
                          TestParameters{ "/Users/niccolo/Progetti/sms_plus_plus_project/UCBlock/netCDF_files/1UC_Data/24/S1ramp1_24.nc4" },
                          TestParameters{ "/Users/niccolo/Progetti/sms_plus_plus_project/UCBlock/netCDF_files/1UC_Data/24/S1ramp2_24.nc4" }
                          ));

/*--------------------------------------------------------------------------*/
/*---------------------------------- MAIN ----------------------------------*/
/*--------------------------------------------------------------------------*/
int main( int argc, char ** argv ) {
 ::testing::InitGoogleTest( &argc, argv );
 return RUN_ALL_TESTS();
}
