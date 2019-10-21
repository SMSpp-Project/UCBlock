/*--------------------------------------------------------------------------*/
/*---------------------------- File dat2nc4.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Small main() for constructing UCBlock netCDF files out of dat and mod ones.
 *
 * \version 0.10
 *
 * \date 20 - 06 - 2019
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Niccolò Iardella \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Niccolò Iardella
 */

#include <iostream>
#include <iomanip>
#include <vector>
#include <getopt.h>

#include <netcdf>
#include <ncByte.h>

#include <UCBlock.h>

/// A thermal unit as represented in a DAT or in a MOD file
struct ThermalUnit {
 unsigned int index{};

 double QuadTerm{};
 double LinearTerm{};
 double ConstTerm{};
 double MinPower{};
 double MaxPower{};
 double InitialPower{}; // aka PZERO
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
 void generate_startupcost() {
  if( coolAndFuelCost != 0 || hotAndFuelCost != 0 ) {
   throw ( std::invalid_argument( "Time-dependent start up costs are not allowed" ) );
  }

  StartUpCost = fixedCost;
  // double coolCost = coolAndFuelCost * ( 1 - exp( -MinDownTime / tau ) ) + fixedCost;
  // double hotCost = hotAndFuelCost * MinDownTime + fixedCost;
  // double cost = coolCost < hotCost ? coolCost : hotCost;
 }

 friend std::ostream & operator<<( std::ostream & out, const ThermalUnit & unit ) {
  unit.print( out );
  return out;
 }

 friend std::istream & operator>>( std::istream & in, ThermalUnit & unit ) {
  unit.load( in );
  return in;
 }
};

/// A hydro unit as represented in a MOD file
struct HydroUnit {
};

/// A hydro cascade unit as represented in a MOD file
struct HydroCascadeUnit {
};

/// A load curve as represented in a MOD file
struct LoadCurve {
 double MinSystemCapacity{};
 double MaxSystemCapacity{};
 double MaxThermalCapacity{};
 std::vector< std::vector< double>> Loads;
 std::vector< double > SpinningReserve;
};

/// A DAT file containing a single thermal unit
struct DatFile {
 unsigned int TimeHorizon{};
 ThermalUnit thermal_unit;
 std::vector< double > Lambda;
 std::vector< double > Mu;

 /// Generates the linear and constant coefficients of the cost function
 void generate_bc( std::vector< double > & b, std::vector< double > & c ) {
  b.resize( TimeHorizon );
  c.resize( TimeHorizon );

  for( unsigned int t = 0; t < TimeHorizon; ++t ) {
   b[ t ] = thermal_unit.LinearTerm - Lambda[ t ];
   c[ t ] = thermal_unit.ConstTerm - Mu[ t ]*thermal_unit.MaxPower;
  }

  // If all elements are identical, we use only one value
  if( std::adjacent_find( b.begin(),
                          b.end(),
                          std::not_equal_to<>() ) == b.end() ) {
   b.resize( 1 );
  }

  if( std::adjacent_find( c.begin(),
                          c.end(),
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
  for( unsigned int t = 0; t < TimeHorizon; ++t ) {
   in >> Lambda[ t ];
  }
  in >> skip;
  for( unsigned int t = 0; t < TimeHorizon; ++t ) {
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
  for( unsigned int t = 0; t < TimeHorizon; ++t ) {
   out << Lambda[ t ] << " ";
  }
  out << "\n";
  out << "Mu" << "\n";
  for( unsigned int t = 0; t < TimeHorizon; ++t ) {
   out << Mu[ t ] << " ";
  }
  out << "\n";
 }

 friend std::ostream & operator<<( std::ostream & out, const DatFile & file ) {
  file.print( out );
  return out;
 }

 friend std::istream & operator>>( std::istream & in, DatFile & file ) {
  file.load( in );
  return in;
 }
};

/// A MOD file containing a load curve and multiple units
struct ModFile {
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
  int rows, columns, elements;

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
  for( int r = 0; r < rows; ++r ) {
   load_curve.Loads[ r ].resize( columns );
   for( int c = 0; c < columns; ++c ) {
    in >> load_curve.Loads[ r ][ c ];
   }
  }

  // "SpinningReserve"
  in >> skip >> elements;
  load_curve.SpinningReserve.resize( elements );
  for( int e = 0; e < elements; ++e ) {
   in >> load_curve.SpinningReserve[ e ];
  }

  // "ThermalSection"
  in >> skip;
  thermal_units.resize( NumThermal );
  for( unsigned int i = 0; i < NumThermal; ++i ) {
   thermal_units[ i ].load( in );
  }

  // "HydroSection"
  in >> skip;
  hydro_units.resize( NumHydro );
  for( unsigned int i = 0; i < NumHydro; ++i ) {}

  // "HydroCascadeSection"
  in >> skip;
  hydro_cascade_units.resize( NumCascade );
  for( unsigned int i = 0; i < NumCascade; ++i ) {}
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
  out << "HydroCascadeSection\n";
 }

 friend std::ostream & operator<<( std::ostream & out, const ModFile & file ) {
  file.print( out );
  return out;
 }

 friend std::istream & operator>>( std::istream & in, ModFile & file ) {
  file.load( in );
  return in;
 }
};

enum filetype {
 ftDat = 0,
 ftMod = 1
};

using namespace SMSpp_di_unipi_it;

DatFile dat_file;
ModFile mod_file;
filetype type;

std::vector< double > b;
std::vector< double > c;

void serialize_unit( netCDF::NcGroup & g, const ThermalUnit & unit ) {
 serialize( g, "MinPower", netCDF::NcDouble(), unit.MinPower );
 serialize( g, "MaxPower", netCDF::NcDouble(), unit.MaxPower );
 serialize( g, "DeltaRampUp", netCDF::NcDouble(), unit.DeltaRampUp );
 serialize( g, "DeltaRampDown", netCDF::NcDouble(), unit.DeltaRampDown );
 serialize( g, "QuadTerm", netCDF::NcDouble(), unit.QuadTerm );
 serialize( g, "StartUpCost", netCDF::NcDouble(), unit.StartUpCost );

 if( type == ftDat ) {
  auto NumberIntervals = g.getDim( "NumberIntervals" );

  if( b.size() == 1 ) {
   serialize( g, "LinearTerm", netCDF::NcDouble(), b[ 0 ] );
  } else {
   serialize( g, "LinearTerm", netCDF::NcDouble(), NumberIntervals, b );
  }

  if( c.size() == 1 ) {
   serialize( g, "ConstTerm", netCDF::NcDouble(), c[ 0 ] );
  } else {
   serialize( g, "ConstTerm", netCDF::NcDouble(), NumberIntervals, c );
  }

 } else { // type == ftMod
  serialize( g, "LinearTerm", netCDF::NcDouble(), unit.LinearTerm );
  serialize( g, "ConstTerm", netCDF::NcDouble(), unit.ConstTerm );
 }

 serialize( g, "InitialPower", netCDF::NcDouble(), unit.InitialPower );
 serialize( g, "InitUpDownTime", netCDF::NcInt64(), unit.InitUpDownTime );
 serialize( g, "MinUpTime", netCDF::NcUint64(), unit.MinUpTime );
 serialize( g, "MinDownTime", netCDF::NcUint64(), unit.MinDownTime );
}

/// Prints usage information
void print_help() {
 // http://docopt.org
 std::cout << "Usage: nc4generator <file>" << std::endl;
}

/// Input file name
std::string filename{};


/// Processes command line arguments
void process_args( int argc, char ** argv ) {
 // It doesn't do much, but it's extendable

 if( argc < 2 ) {
  print_help();
  exit( 1 );
 }

 const char * const short_opts = "h";
 const option long_opts[] = {
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
 if (optind < argc) {
  filename = std::string( argv[ optind ] );
 } else {
  print_help();
  exit( 1 );
 }
}

int main( int argc, char ** argv ) {

 process_args( argc, argv );

 std::ifstream inputFile( filename );
 if( !inputFile.is_open() ) {
  std::cerr << "Error: cannot open file " << filename << std::endl;
  return 1;
 }

 if( filename.size() < 5 ) {
  std::cerr << "Error: File name is too short" << std::endl;
  inputFile.close();
  return 1;
 }

 std::string ext = filename.substr( filename.size() - 4, 4 );
 std::string dat( ".dat" );
 std::string mod( ".mod" );

 if( std::equal( ext.begin(), ext.end(), dat.begin(),
                 []( auto a, auto b ) {
                  return ( std::tolower( a ) == std::tolower( b ) );
                 } ) ) {
  type = ftDat;
 } else if( std::equal( ext.begin(), ext.end(), mod.begin(),
                        []( auto a, auto b ) {
                         return ( std::tolower( a ) == std::tolower( b ) );
                        } ) ) {
  type = ftMod;
 } else {
  std::cerr << "Error: Supported file formats are: dat, mod." << std::endl;
  inputFile.close();
  return 1;
 }

 filename.erase( filename.size() - 4, 4 );
 filename.append( ".nc4" );

 if( type == ftDat ) {
  inputFile >> dat_file;
  dat_file.generate_bc( b, c );
  std::cout << dat_file;
 } else { // type == ftMod
  inputFile >> mod_file;
  std::cout << mod_file;
 }
 inputFile.close();

 netCDF::NcFile f( filename, netCDF::NcFile::replace );
 f.putAtt( "SMS++_file_type", netCDF::NcInt(), eBlockFile );

 if( type == ftDat ) {

  auto bg = f.addGroup( "Block_0" );
  bg.putAtt( "type", "ThermalUnitBlock" );
  bg.addDim( "TimeHorizon", dat_file.TimeHorizon );
  bg.addDim( "NumberIntervals", dat_file.TimeHorizon );
  serialize_unit( bg, dat_file.thermal_unit );

 } else { // type == ftMod

  auto bg = f.addGroup( "Block_0" );
  bg.putAtt( "type", "UCBlock" );
  bg.addDim( "TimeHorizon", mod_file.TimeHorizon );
  bg.addDim( "NumberUnits", mod_file.NumThermal );
  // ubg.addDim( "NumberIntervals", 1 );

  for( unsigned int i = 0; i < mod_file.NumThermal; ++i ) {
   auto ug = bg.addGroup( "UnitBlock_" + std::to_string( i ) );
   // auto ug = bg.addGroup( "Unit_" + std::to_string( bg.getGroupCount() ) );
   ug.putAtt( "type", "ThermalUnitBlock" );
   mod_file.thermal_units[ i ].generate_startupcost();
   serialize_unit( ug, mod_file.thermal_units[ i ] );
  }

  for( unsigned int i = 0; i < mod_file.TimeHorizon; ++i ) {
   auto ng = bg.addGroup( "NetworkBlock_" + std::to_string( i ) );
   ng.putAtt( "type", "BusNetworkBlock" );
   ng.addDim( "NumberNodes", 1 );
   // FIXME: Check Loads[][] bounds
   serialize( ng,
              "ActiveDemand",
              netCDF::NcDouble(),
              mod_file.load_curve.Loads[ 0 ][ i ] );
  }
 }
 return 0;
}
