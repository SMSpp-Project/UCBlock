#include <iostream>
#include <fstream>

#include <AbstractBlock.h>
#include <UCBlock.h>
#include <ThermalUnitBlock.h>
#include <BusNetworkBlock.h>
#include <CPXMILPSolver.h>

using namespace SMSpp_di_unipi_it;

int main( int argc, char ** argv ) {

 std::string filename( argv[ 1 ] );

 netCDF::NcFile f( filename, netCDF::NcFile::read );
 if( f.isNull() ) {
  std::cerr << "cannot open nc4 file " << filename << std::endl;
  exit( 1 );
 }

 netCDF::NcGroupAtt gtype = f.getAtt( "SMS++_file_type" );
 if( gtype.isNull() ) {
  std::cerr << filename << " is not an SMS++ nc4 file" << std::endl;
  exit( 1 );
 }

 int type;
 gtype.getValues( &type );

 if( type != eBlockFile ) {
  std::cerr << filename << " is not an SMS++ nc4 Block file" << std::endl;
  exit( 1 );
 }

 netCDF::NcGroup bg = f.getGroup( "Block_0" );
 if( bg.isNull() ) {
  std::cerr << "Block_0 empty or undefined in " << filename << std::endl;
  exit( 1 );
 }

 // Deserialize
 auto tub = dynamic_cast<ThermalUnitBlock *>(Block::new_Block( "ThermalUnitBlock" ));
 tub->deserialize( bg );

 // netCDF::NcFile f1( "test.nc4", netCDF::NcFile::replace );
 // f1.putAtt( "SMS++_file_type", netCDF::NcInt(), eBlockFile );
 // auto bg1 = f1.addGroup( "Block_0" );
 // tub->serialize( bg1 );

 // Generate abstract representation
 int tmp = 15;
 SimpleConfiguration<int> myconfig(tmp);

 tub->set_verbosity(Block::high);
 tub->generate_abstract_variables( &myconfig );
 tub->generate_abstract_constraints( nullptr );
 tub->generate_objective( nullptr);
 std::cout << *tub;

 // Register solver
 // Solver * solver = Solver::new_Solver( "CPXMILPSolver" );
 Solver* solver = new CPXMILPSolver();
 tub->register_Solver( solver );

 // Write problem
 dynamic_cast<CPXMILPSolver*>(solver)->write_lp("output.lp");
 // Solve
 int status = solver->compute();
 auto ub = solver->get_ub();

 // Retrieve objective function
 auto obj = dynamic_cast<FRealObjective *>(tub->get_objective());
 auto obj_f = obj->get_function();

 std::cout << "Status = " << status << std::endl;
 std::cout << "Upper bound = " << ub << std::endl;
 std::cout << "Function value =  " << obj_f->get_value() << std::endl;

 return 0;
}
