#include <iostream>
#include <getopt.h>

#include <UCBlock.h>
#include <ThermalUnitBlock.h>
#include <BusNetworkBlock.h>
// #include <CPXMILPSolver.h>

using namespace SMSpp_di_unipi_it;

std::string filename{};
std::string lp_file{};
std::string solver_name{};

void print_help() {
 // http://docopt.org
 std::cout << "Usage: uc_solver [options] <nc4-file>" << std::endl
           << std::endl
           << "Options:" << std::endl
           << "  -s <solver>, --solver <solver>  Choose solver." << std::endl
           << "                                  Available solvers are: cplex, dp." << std::endl
           << "  -w <file>, --writelp <file>     Write LP problem(s) on file(s)." << std::endl
           << "  -h, --help                      Print this help." << std::endl;
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
 if( optind < argc ) {
  filename = std::string( argv[ optind ] );
 } else {
  print_help();
  exit( 1 );
 }
}

int main( int argc, char ** argv ) {

 solver_name = "cplex";
 process_args( argc, argv );

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

 // Deserialize block
 auto ucb = dynamic_cast<UCBlock *>(Block::new_Block( "UCBlock" ));
 ucb->deserialize( bg );

 // Configure blocks
 auto conf = new BlockConfig();
 for (int i = 0; i < 10; ++i) {
  auto subconf = new BlockConfig();
  subconf->f_static_variables_Configuration = new SimpleConfiguration< int >( 15 );
  conf->v_sub_BlockConfig.emplace_back(subconf);
 }

 // // Configure solver
 auto slv_conf = new BlockSolverConfig();

 if( solver_name == "cplex" ) {
  slv_conf->v_SolverNames.emplace_back( "CPXMILPSolver" );
  slv_conf->v_SolverConfigs.emplace_back( new ComputeConfig() );

 } else if( solver_name == "dp" ) {
  std::cerr << "Sorry, DP Solver is not available yet..." << std::endl;
  exit( 0 );
 } else {
  std::cerr << "Available solvers are: cplex, dp" << std::endl;
  exit( 1 );
 }

 // for (auto b : ucb->get_nested_Blocks()) {
 //  auto tub = dynamic_cast<ThermalUnitBlock*>(b);
 //  if (tub) {
 //   tub->set_BlockConfig( conf );
 //   tub->set_SolverConfig( slv_conf );
 //  } else {
 //   auto bnb = dynamic_cast<BusNetworkBlock *>(b);
 //   if(bnb ) {
 //    bnb->generate_abstract_variables( nullptr );
 //   }
 //  }
 // }

 ucb->set_BlockConfig( conf );
 ucb->set_SolverConfig( slv_conf );

 // int i = 0;
 // double acc = 0;
 // for (auto b : ucb->get_nested_Blocks()) {
 //  auto tub = dynamic_cast<ThermalUnitBlock*>(b);
 //  if (tub) {
 //   auto solver = tub->get_registered_solvers().front();
 //
 //   // // Write LP problem
 //   // // TODO: Use configuration instead, so no dependency from CPXMILPSolver
 //   // if( !lp_file.empty() ) {
 //   //  dynamic_cast<CPXMILPSolver *>(solver)->write_lp( lp_file + std::to_string(i) + ".lp" );
 //   // }
 //
 //   // Solve
 //   int status = solver->compute();
 //   auto ub = solver->get_ub();
 //   auto lb = solver->get_lb();
 //
 //   auto obj = dynamic_cast<FRealObjective *>(tub->get_objective());
 //   auto obj_f = obj->get_function();
 //   auto obj_value = obj_f->get_value();
 //   acc += obj_value;
 //
 //   std::cout << "Block " << i << std::endl;
 //   std::cout << "Status = " << status << std::endl;
 //   std::cout << "Upper bound = " << ub << std::endl;
 //   std::cout << "Lower bound = " << lb << std::endl;
 //   std::cout << std::endl;
 //  }
 //  ++i;
 // }
 // std::cout << "Sum of objective values = " << acc << std::endl;

 auto solver = ucb->get_registered_solvers().front();
 int status = solver->compute();
 auto ub = solver->get_ub();
 auto lb = solver->get_lb();
 std::cout << "Status = " << status << std::endl;
 std::cout << "Upper bound = " << ub << std::endl;
 std::cout << "Lower bound = " << lb << std::endl;

 return 0;
}
