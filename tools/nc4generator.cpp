/*--------------------------------------------------------------------------*/
/*------------------------- File nc4generator.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Small main() for constructing UCBlock netCDF files out of dat and mod ones.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Niccolo' Iardella \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Niccolo' Iardella
 */

#include <iomanip>
#include <getopt.h>

#include <SMSTypedefs.h>

/*--------------------------------------------------------------------------*/
/*------------------------------- ThermalUnit ------------------------------*/
/*--------------------------------------------------------------------------*/

/// A thermal unit as represented in a DAT or in a MOD file
struct ThermalUnit
{
 unsigned int index{};

 double QuadTerm{};
 double LinearTerm{};
 double ConstTerm{};
 double MinPower{};
 double MaxPower{};
 double InitialPower{};  // aka PZERO
 double InitUpDownTime{};
 double MinUpTime{};
 double MinDownTime{};
 double DeltaRampUp{};
 double DeltaRampDown{};

 // Only in MOD files
 double coolAndFuelCost{};
 double hotAndFuelCost{};
 double tau{};
 double tauMax{};
 double fixedCost{};
 double SUCC{};

 // Only in DAT files
 double StartUpCost{};
 double BoundOn{};
 double BoundDown{};

 /// Loads the data from a line of a MOD file
 void load( std::istream & in ) {
  std::string skip;
  in >> index
     >> QuadTerm
     >> LinearTerm
     >> ConstTerm
     >> MinPower
     >> MaxPower
     >> InitUpDownTime
     >> MinUpTime
     >> MinDownTime
     >> coolAndFuelCost
     >> hotAndFuelCost
     >> tau
     >> tauMax
     >> fixedCost
     >> SUCC
     >> InitialPower;
  // "RampConstraints"
  in >> skip >> DeltaRampUp >> DeltaRampDown;
 }

 /// Prints the data
 void print( std::ostream & out ) const {
  out << index << "\t"
      << QuadTerm << "\t"
      << LinearTerm << "\t"
      << ConstTerm << "\t"
      << MinPower << "\t"
      << MaxPower << "\t"
      << int( InitUpDownTime ) << "\t"
      << int( MinUpTime ) << "\t"
      << int( MinDownTime ) << "\t"
      << coolAndFuelCost << "\t"
      << hotAndFuelCost << "\t"
      << int( tau ) << "\t"
      << int( tauMax ) << "\t"
      << fixedCost << "\t"
      << SUCC << "\t"
      << InitialPower << "\n";
  out << "RampConstraints\t "
      << DeltaRampUp << " \t " << DeltaRampDown << "\n";
 }

 /// Generates StartUpCost from
 void generate_startup_cost( void ) {
  if( ( coolAndFuelCost != 0 ) || ( hotAndFuelCost != 0 ) ) {
   throw( std::invalid_argument(
    "Time-dependent start up costs are not allowed" ) );
  }

  StartUpCost = fixedCost;
  // double coolCost = coolAndFuelCost * ( 1 - exp( -MinDownTime / tau ) ) + fixedCost;
  // double hotCost = hotAndFuelCost * MinDownTime + fixedCost;
  // double cost = coolCost < hotCost ? coolCost : hotCost;
 }

 friend std::ostream &
 operator<<( std::ostream & out , const ThermalUnit & unit ) {
  unit.print( out );
  return( out );
 }

 friend std::istream & operator>>( std::istream & in , ThermalUnit & unit ) {
  unit.load( in );
  return( in );
 }
};

/*--------------------------------------------------------------------------*/
/*-------------------------------- HydroUnit -------------------------------*/
/*--------------------------------------------------------------------------*/

/// A hydro unit as represented in a MOD file
struct HydroUnit
{
 unsigned int index{};

 double volumeToPower{};
 double b_h{};
 double maxUsage{};
 double maxSpillage{};
 double initialFlood{};
 double minFlood{};
 double maxFlood{};

 std::vector< double > inflows;

 /// Loads the data from a line of a MOD file
 void load( std::istream & in ) {
  std::string skip;
  in >> index
     >> volumeToPower
     >> b_h
     >> maxUsage
     >> maxSpillage
     >> initialFlood
     >> minFlood
     >> maxFlood;
  for( double & inflow : inflows ) {
   in >> inflow;
  }
 }

 /// Prints the data
 void print( std::ostream & out ) const {
  out << index << "\t"
      << volumeToPower << "\t"
      << b_h << "\t"
      << maxUsage << "\t"
      << maxSpillage << "\t"
      << initialFlood << "\t"
      << minFlood << "\t"
      << maxFlood << "\n";
  for( double i : inflows ) {
   out << "\t" << i << "\t";
  }
  out << "\n";
 }

 friend std::ostream &
 operator<<( std::ostream & out , const HydroUnit & unit ) {
  unit.print( out );
  return( out );
 }

 friend std::istream & operator>>( std::istream & in , HydroUnit & unit ) {
  unit.load( in );
  return( in );
 }
};

/*--------------------------------------------------------------------------*/
/*----------------------------- HydroCascadeUnit ---------------------------*/
/*--------------------------------------------------------------------------*/

/// A hydro cascade unit as represented in a MOD file
struct HydroCascadeUnit
{
};

/*--------------------------------------------------------------------------*/
/*-------------------------------- LoadCurve -------------------------------*/
/*--------------------------------------------------------------------------*/

/// A load curve as represented in a MOD file
struct LoadCurve
{
 double MinSystemCapacity{};
 double MaxSystemCapacity{};
 double MaxThermalCapacity{};
 std::vector< std::vector< double > > Loads;
 std::vector< double > SpinningReserve;
};

/*--------------------------------------------------------------------------*/
/*-------------------------------- DAT file --------------------------------*/
/*--------------------------------------------------------------------------*/

/// A DAT file containing a single thermal unit
struct DatFile
{
 unsigned int TimeHorizon{};
 ThermalUnit thermal_unit;
 std::vector< double > Lambda;
 std::vector< double > Mu;

 /// Generates the linear and constant coefficients of the cost function
 void generate_bc( std::vector< double > & b , std::vector< double > & c ) {
  b.resize( TimeHorizon );
  c.resize( TimeHorizon );

  for( unsigned int t = 0 ; t < TimeHorizon ; ++t ) {
   b[ t ] = thermal_unit.LinearTerm - Lambda[ t ];
   c[ t ] = thermal_unit.ConstTerm - Mu[ t ] * thermal_unit.MaxPower;
  }

  // If all elements are identical, we use only one value
  if( std::adjacent_find( b.begin() ,
                          b.end() ,
                          std::not_equal_to<>() ) == b.end() ) {
   b.resize( 1 );
  }

  if( std::adjacent_find( c.begin() ,
                          c.end() ,
                          std::not_equal_to<>() ) == c.end() ) {
   c.resize( 1 );
  }
 }

 /// Loads the data from a DAT file
 void load( std::istream & in ) {
  std::string skip;
  in >> skip >> TimeHorizon
     >> skip >> thermal_unit.QuadTerm
     >> skip >> thermal_unit.LinearTerm
     >> skip >> thermal_unit.ConstTerm
     >> skip >> thermal_unit.MinPower
     >> skip >> thermal_unit.MaxPower
     >> skip >> thermal_unit.InitUpDownTime
     >> skip >> thermal_unit.MinUpTime
     >> skip >> thermal_unit.MinDownTime
     >> skip >> thermal_unit.StartUpCost
     >> skip >> thermal_unit.InitialPower
     >> skip >> thermal_unit.DeltaRampUp
     >> skip >> thermal_unit.DeltaRampDown
     >> skip >> thermal_unit.BoundOn
     >> skip >> thermal_unit.BoundDown;

  Lambda.resize( TimeHorizon );
  Mu.resize( TimeHorizon );
  in >> skip;
  for( unsigned int t = 0 ; t < TimeHorizon ; ++t ) {
   in >> Lambda[ t ];
  }
  in >> skip;
  for( unsigned int t = 0 ; t < TimeHorizon ; ++t ) {
   in >> Mu[ t ];
  }
 }

 /// Prints the data
 void print( std::ostream & out ) const {

  /*
   * "Term" is misspelled in the labels because it is misspelled in the
   * original input files and we want to compare the output with them.
   */
  out << "TimeHorizon\t" << TimeHorizon << "\n"
      << "QuadTherm\t" << thermal_unit.QuadTerm << "\n"
      << "LinearTherm\t" << thermal_unit.LinearTerm << "\n"
      << "ConstTherm\t" << thermal_unit.ConstTerm << "\n"
      << "MinPower\t" << thermal_unit.MinPower << "\n"
      << "MaxPower\t" << thermal_unit.MaxPower << "\n"
      << "InitUpDownTime\t" << thermal_unit.InitUpDownTime << "\n"
      << "MinUpTime\t" << thermal_unit.MinUpTime << "\n"
      << "MinDownTime\t" << thermal_unit.MinDownTime << "\n"
      << "StartUpCost\t" << thermal_unit.StartUpCost << "\n"
      << "PZERO\t\t" << thermal_unit.InitialPower << "\n"
      << "DeltaRampUp\t" << thermal_unit.DeltaRampUp << "\n"
      << "DeltaRampDown\t" << thermal_unit.DeltaRampDown << "\n"
      << "BoundOn\t\t" << thermal_unit.BoundOn << "\n"
      << "BoundDown\t" << thermal_unit.BoundDown << "\n";

  out << "Lambda" << "\n";
  for( unsigned int t = 0 ; t < TimeHorizon ; ++t ) {
   out << Lambda[ t ] << " ";
  }
  out << "\n";
  out << "Mu" << "\n";
  for( unsigned int t = 0 ; t < TimeHorizon ; ++t ) {
   out << Mu[ t ] << " ";
  }
  out << "\n";
 }

 friend std::ostream & operator<<( std::ostream & out , const DatFile & file ) {
  file.print( out );
  return( out );
 }

 friend std::istream & operator>>( std::istream & in , DatFile & file ) {
  file.load( in );
  return( in );
 }
};

/*--------------------------------------------------------------------------*/
/*-------------------------------- MOD file --------------------------------*/
/*--------------------------------------------------------------------------*/

/// A MOD file containing a load curve and multiple units
struct ModFile
{
 unsigned int ProblemNum{};
 unsigned int TimeHorizon{};
 unsigned int NumThermal{};
 unsigned int NumHydro{};
 unsigned int NumCascade{};

 LoadCurve load_curve;
 std::vector< ThermalUnit > thermal_units;
 std::vector< HydroUnit > hydro_units;
 std::vector< HydroCascadeUnit > hydro_cascade_units;

 /// Loads the data from a MOD file
 void load( std::istream & in ) {
  std::string skip;
  int rows , columns , elements;

  in >> skip >> ProblemNum;
  in >> skip >> TimeHorizon;
  in >> skip >> NumThermal;
  in >> skip >> NumHydro;
  in >> skip >> NumCascade;

  // "LoadCurve"
  in >> skip;
  in >> skip >> load_curve.MinSystemCapacity;
  in >> skip >> load_curve.MaxSystemCapacity;
  in >> skip >> load_curve.MaxThermalCapacity;

  // "Loads"
  in >> skip >> rows >> columns;
  load_curve.Loads.resize( rows );
  for( int r = 0 ; r < rows ; ++r ) {
   load_curve.Loads[ r ].resize( columns );
   for( int c = 0 ; c < columns ; ++c ) {
    in >> load_curve.Loads[ r ][ c ];
   }
  }

  // "SpinningReserve"
  in >> skip >> elements;
  load_curve.SpinningReserve.resize( elements );
  for( int e = 0 ; e < elements ; ++e ) {
   in >> load_curve.SpinningReserve[ e ];
  }

  // "ThermalSection"
  in >> skip;
  thermal_units.resize( NumThermal );
  for( unsigned int i = 0 ; i < NumThermal ; ++i ) {
   thermal_units[ i ].load( in );
  }

  // "HydroSection"
  in >> skip;
  hydro_units.resize( NumHydro );
  for( unsigned int i = 0 ; i < NumHydro ; ++i ) {
   hydro_units[ i ].inflows.resize( TimeHorizon );
   hydro_units[ i ].load( in );
  }

  // "HydroCascadeSection"
  in >> skip;
  hydro_cascade_units.resize( NumCascade );
  for( unsigned int i = 0 ; i < NumCascade ; ++i ) {}
 }

 /// Prints the data
 void print( std::ostream & out ) const {
  out << std::fixed << std::setprecision( 6 );

  out << "ProblemNum " << ProblemNum << "\n"
      << "HorizonLen " << TimeHorizon << "\n"
      << "NumThermal " << NumThermal << "\n"
      << "NumHydro " << NumHydro << "\n"
      << "NumCascade " << NumCascade << "\n";

  out << "LoadCurve\n"
      << "MinSystemCapacity \t " << load_curve.MinSystemCapacity << "\n"
      << "MaxSystemCapacity \t " << load_curve.MaxSystemCapacity << "\n"
      << "MaxThermalCapacity\t " << load_curve.MaxThermalCapacity << "\n";

  out << "Loads\t"
      << load_curve.Loads.size() << "\t"
      << load_curve.Loads[ 0 ].size() << "\n";
  for( auto & row : load_curve.Loads ) {
   for( double c : row ) {
    out << c << "\t";
   }
   out << "\n";
  }

  out << "SpinningReserve\t"
      << load_curve.SpinningReserve.size() << " \n";
  for( double c : load_curve.SpinningReserve ) {
   out << c << "\t";
  }
  out << "\n";

  out << "ThermalSection\n";
  for( auto & unit : thermal_units ) {
   unit.print( out );
  }

  out << "HydroSection\n";
  for( auto & unit : hydro_units ) {
   unit.print( out );
  }
  out << "HydroCascadeSection\n";
 }

 friend std::ostream & operator<<( std::ostream & out , const ModFile & file ) {
  file.print( out );
  return( out );
 }

 friend std::istream & operator>>( std::istream & in , ModFile & file ) {
  file.load( in );
  return( in );
 }
};

enum filetype
{
 ftDat = 0 ,
 ftMod = 1
};

/*--------------------------------------------------------------------------*/
/*------------------------------ Other stuff -------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

DatFile dat_file;
ModFile mod_file;
filetype type;

std::vector< double > b;
std::vector< double > c;

/*--------------------------------------------------------------------------*/

void serialize_thermalunit( netCDF::NcGroup & g , const ThermalUnit & unit ) {
 serialize( g , "MinPower" , netCDF::NcDouble() , unit.MinPower );
 serialize( g , "MaxPower" , netCDF::NcDouble() , unit.MaxPower );
 serialize( g , "DeltaRampUp" , netCDF::NcDouble() , unit.DeltaRampUp );
 serialize( g , "DeltaRampDown" , netCDF::NcDouble() , unit.DeltaRampDown );
 serialize( g , "QuadTerm" , netCDF::NcDouble() , unit.QuadTerm );
 serialize( g , "StartUpCost" , netCDF::NcDouble() , unit.StartUpCost );

 if( type == ftDat ) {
  auto NumberIntervals = g.getDim( "NumberIntervals" );

  if( b.size() == 1 ) {
   serialize( g , "LinearTerm" , netCDF::NcDouble() , b[ 0 ] );
  } else {
   serialize( g , "LinearTerm" , netCDF::NcDouble() , NumberIntervals , b );
  }

  if( c.size() == 1 ) {
   serialize( g , "ConstTerm" , netCDF::NcDouble() , c[ 0 ] );
  } else {
   serialize( g , "ConstTerm" , netCDF::NcDouble() , NumberIntervals , c );
  }

 } else {  // type == ftMod
  serialize( g , "LinearTerm" , netCDF::NcDouble() , unit.LinearTerm );
  serialize( g , "ConstTerm" , netCDF::NcDouble() , unit.ConstTerm );
 }

 serialize( g , "InitialPower" , netCDF::NcDouble() , unit.InitialPower );
 serialize( g , "InitUpDownTime" , netCDF::NcInt64() , unit.InitUpDownTime );
 serialize( g , "MinUpTime" , netCDF::NcUint64() , unit.MinUpTime );
 serialize( g , "MinDownTime" , netCDF::NcUint64() , unit.MinDownTime );
}

/*--------------------------------------------------------------------------*/

void serialize_hydrounit( netCDF::NcGroup & g , const HydroUnit & unit ) {
 serialize( g , "LinearTerm" , netCDF::NcDouble() , unit.volumeToPower );
 serialize( g , "MaxFlow" , netCDF::NcDouble() , unit.maxSpillage );
 serialize( g , "MaxPower" , netCDF::NcDouble() ,
            ( unit.maxUsage * unit.volumeToPower ) );
 serialize( g , "InitialVolumetric" , netCDF::NcDouble() , unit.initialFlood );
 serialize( g , "MinVolumetric" , netCDF::NcDouble() , unit.minFlood );
 serialize( g , "MaxVolumetric" , netCDF::NcDouble() , unit.maxFlood );

 auto TimeHorizon = g.getParentGroup().getDim( "TimeHorizon" );
 serialize( g , "Inflows" , netCDF::NcDouble() , TimeHorizon , unit.inflows );
}

/*--------------------------------------------------------------------------*/


std::string input_path{};     ///< input file name
std::string output_path{};    ///< input file name
bool verbose = false;         ///< if the tool should be verbose
std::string exe{};            ///< name of the executable file
std::string docopt_desc{};    ///< tool description

/*--------------------------------------------------------------------------*/

/// Gets the name of the executable from its full path
std::string get_filename( const std::string & fullpath ) {
 std::size_t found = fullpath.find_last_of( "/\\" );
 return( fullpath.substr( found + 1 ) );
}

/*--------------------------------------------------------------------------*/

/// Prints the tool description and usage
void docopt( void ) {
 // http://docopt.org
 std::cout << docopt_desc << std::endl;
 std::cout << "Usage:\n"
           << "  " << exe << " [-v] <input>\n"
           << "  " << exe << " -h | --help\n"
           << std::endl
           << "Options:\n"
           << "  -v, --verbose  Make the tool verbose.\n"
           << "  -h, --help     Print this help.\n";
}

/*--------------------------------------------------------------------------*/

/// Processes command line arguments
void process_args( int argc , char ** argv ) {

 const char * const short_opts = "vh";
 const option long_opts[] = {
  { "verbose" , no_argument , nullptr , 'v' } ,
  { "help" ,    no_argument , nullptr , 'h' } ,
  { nullptr ,   no_argument , nullptr , 0 }
 };

 // Options
 while( true ) {
  const auto opt = getopt_long( argc , argv , short_opts , long_opts ,
                                nullptr );

  if( -1 == opt ) {
   break;
  }
  switch( opt ) {
   case 'v':
    verbose = true;
    break;
   case 'h':
    docopt();
    exit( 0 );
   case '?':
   default:
    std::cout << "Try " << exe << "' --help' for more information.\n";
    exit( 1 );
  }
 }

 // Last argument
 if( optind < argc ) {
  input_path = std::string( argv[ optind ] );
 } else {
  std::cout << exe << ": no input file\n"
            << "Try " << exe << "' --help' for more information.\n";
  exit( 1 );
 }
}

/*--------------------------------------------------------------------------*/
/*---------------------------------- MAIN ----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( int argc , char ** argv ) {

 // Manage options and help
 docopt_desc = "NC4 Thermal Unit generator.\n";
 exe = get_filename( argv[ 0 ] );
 process_args( argc , argv );

 // // Check if input file exists
 // if( ! std::filesystem::exists( input_path ) ) {
 //  std::cerr << exe << ": cannot open file " << input_path << std::endl;
 //  exit( 1 );
 // }

 // Check if input file can be opened
 std::ifstream input_file( input_path );
 if( ! input_file.is_open() ) {
  std::cerr << exe << ": cannot open file " << input_path << std::endl;
  exit( 1 );
 }

 std::string ext = input_path.substr( input_path.size() - 4 , 4 );
 std::string dat( ".dat" );
 std::string mod( ".mod" );

 // Check input file type
 if( std::equal( ext.begin() , ext.end() , dat.begin() ,
                 []( auto a , auto b ) {
                  return( std::tolower( a ) == std::tolower( b ) );
                 } ) ) {
  type = ftDat;
 } else if( std::equal( ext.begin() , ext.end() , mod.begin() ,
                        []( auto a , auto b ) {
                         return( std::tolower( a ) == std::tolower( b ) );
                        } ) ) {
  type = ftMod;
 } else {
  std::cerr << "Error: Supported file formats are: dat, mod." << std::endl;
  input_file.close();
  return( 1 );
 }

 // Read input file
 if( type == ftDat ) {
  // Read DAT file
  input_file >> dat_file;
  dat_file.generate_bc( b , c );

  if( verbose ) {
   std::cout << dat_file;
  }

 } else {  // type == ftMod
  // Read MOD file
  input_file >> mod_file;

  if( verbose ) {
   std::cout << mod_file;
  }
 }

 // Generate output
 input_file.close();
 output_path = input_path;
 output_path.erase( output_path.size() - 4 , 4 );
 output_path.append( ".nc4" );

 netCDF::NcFile f( output_path , netCDF::NcFile::replace );
 f.putAtt( "SMS++_file_type" , netCDF::NcInt() , eBlockFile );

 if( type == ftDat ) {

  auto bg = f.addGroup( "Block_0" );
  bg.putAtt( "type" , "ThermalUnitBlock" );
  bg.addDim( "TimeHorizon" , dat_file.TimeHorizon );
  bg.addDim( "NumberIntervals" , dat_file.TimeHorizon );
  serialize_thermalunit( bg , dat_file.thermal_unit );

 } else {  // type == ftMod

  auto bg = f.addGroup( "Block_0" );
  bg.putAtt( "type" , "UCBlock" );
  bg.addDim( "TimeHorizon" , mod_file.TimeHorizon );
  auto time_h = bg.getDim( "TimeHorizon" );
  bg.addDim( "NumberUnits" , mod_file.NumThermal +
                             mod_file.NumHydro +
                             mod_file.NumCascade );
  bg.addDim( "NumberIntervals" , 1 );

  auto ng = bg.addGroup( "NetworkData" );
  ng.addDim( "NumberNodes" , 1 );

  serialize( bg ,
             "ActivePowerDemand" ,
             netCDF::NcDouble() ,
             time_h ,
             mod_file.load_curve.Loads[ 0 ] );

  unsigned int num_units = 0;
  for( unsigned int i = 0 ; i < mod_file.NumThermal ; ++i , ++num_units ) {
   auto ug = bg.addGroup( "UnitBlock_" + std::to_string( num_units ) );
   ug.putAtt( "type" , "ThermalUnitBlock" );
   mod_file.thermal_units[ i ].generate_startup_cost();
   serialize_thermalunit( ug , mod_file.thermal_units[ i ] );
  }

  for( unsigned int i = 0 ; i < mod_file.NumHydro ; ++i , ++num_units ) {
   auto ug = bg.addGroup( "UnitBlock_" + std::to_string( num_units ) );
   ug.putAtt( "type" , "HydroUnitBlock" );
   serialize_hydrounit( ug , mod_file.hydro_units[ i ] );
  }
 }

 std::cout << "Output written on " << output_path << std::endl;
 return( 0 );
}

/*--------------------------------------------------------------------------*/
/*----------------------- End File nc4generator.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
