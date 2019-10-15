#include <iostream>
#include <getopt.h>

#include <AbstractBlock.h>
#include <ThermalUnitBlock.h>
#include <CPXMILPSolver.h>

using namespace SMSpp_di_unipi_it;

std::string filename{};
std::string lp_file{};
std::string solver_name( "cplex" );

void print_help() {
 // http://docopt.org
 std::cout << "Usage: 1uc_solver [options] <file>" << std::endl
           << std::endl
           << "-s <solver>, --solver <solver>  Choose solver." << std::endl
           << "                                Available solvers are: cplex, dp." << std::endl
           << "-w <file>, --writelp <file>     Write LP problem on file." << std::endl
           << "-h, --help                      Print this help." << std::endl;
}

void process_args( int argc, char ** argv ) {

 if( argc < 2 ) {
  print_help();
  exit( 1 );
 }

 const char * const short_opts = "s:w:h";
 const option long_opts[] = {
  { "solver",  required_argument, nullptr, 's' },
  { "writelp", required_argument, nullptr, 'w' },
  { "help",    no_argument,       nullptr, 'h' },
  { nullptr,   no_argument,       nullptr, 0 }
 };

 // Options
 while( true ) {
  const auto opt = getopt_long( argc, argv, short_opts, long_opts, nullptr );

  if( -1 == opt ) {
   break;
  }

  switch( opt ) {
   case 's':
    solver_name = std::string( optarg );
   case 'w':
    lp_file = std::string( optarg );
    break;
   case 'h': // -h or --help
    print_help();
    exit( 0 );
   case '?': // Unrecognized option
   default:
    print_help();
    exit( 1 );
  }
 }

 // Last argument
 filename = std::string( argv[ optind ] );
}

int main( int argc, char ** argv ) {

 process_args( argc, argv );

 Solver * solver;
 if( solver_name == "cplex" ) {
  // Solver * solver = Solver::new_Solver( "CPXMILPSolver" );
  solver = new CPXMILPSolver();
 } else if( solver_name == "dp" ) {
  std::cerr << "Sorry, DP Solver is not available yet..." << std::endl;
  exit( 0 );
 } else {
  std::cerr << "Available solvers are: cplex, dp" << std::endl;
  exit( 1 );
 }

 netCDF::NcFile f;
 try {
  f.open( filename, netCDF::NcFile::read );
 } catch( netCDF::exceptions::NcException & e ) {
  std::cerr << "Cannot open nc4 file " << filename << std::endl;
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

 // Generate abstract representation
 int tmp = 15;
 SimpleConfiguration< int > myconfig( tmp );

 tub->set_verbosity( Block::high );
 tub->generate_abstract_variables( &myconfig );
 tub->generate_abstract_constraints( nullptr );
 tub->generate_objective( nullptr );
 std::cout << *tub;

 // Register solver
 tub->register_Solver( solver );

 // Write problem
 if( !lp_file.empty() ) {
  dynamic_cast<CPXMILPSolver *>(solver)->write_lp( lp_file );
 }
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
