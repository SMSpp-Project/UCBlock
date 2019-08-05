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
#include <UCBlock.h>
#include <netcdf>
#include <ncByte.h>

#include "helper.h"

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

 if (type == ftDat) {
  auto NumberIntervals = g.getDim( "NumberIntervals" );

  std::cout << "b" << "\n";
  for( unsigned int t = 0; t < b.size(); ++t ) {
   std::cout << b[ t ] << " ";
  }
  std::cout << "\n";

  std::cout << "c" << "\n";
  for( unsigned int t = 0; t < c.size(); ++t ) {
   std::cout << c[ t ] << " ";
  }
  std::cout << "\n";

  if (b.size() == 1) {
   serialize( g, "LinearTerm", netCDF::NcDouble(), b[0] );
  } else {
   serialize( g, "LinearTerm", netCDF::NcDouble(), NumberIntervals, b );
  }

  if (c.size() == 1) {
   serialize( g, "ConstTerm", netCDF::NcDouble(), c[0] );
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

int main( int argc, char ** argv ) {

 std::ifstream inputFile( argv[ 1 ] );
 if( !inputFile.is_open() ) {
  std::cerr << "Error: cannot open file " << argv[ 1 ] << std::endl;
  return 1;
 }

 std::string filename( argv[ 1 ] );
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
  std::cerr << "Error: Unsupported file type" << std::endl;
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
   auto ng = bg.addGroup( "Network_" + std::to_string( i ) );
   ng.putAtt( "type", "NetworkBlock" );
   ng.addDim( "NumberNodes", 1 );
  }
 }
 return 0;
}
