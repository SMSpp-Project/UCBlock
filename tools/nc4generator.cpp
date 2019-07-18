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
#include <fstream>
#include <UCBlock.h>
#include <ThermalUnitBlock.h>
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
 serialize( g, "InitialPower", netCDF::NcDouble(), unit.InitialPower );
 serialize( g, "MinUpTime", netCDF::NcUint64(), unit.MinUpTime );
 serialize( g, "MinDownTime", netCDF::NcUint64(), unit.MinDownTime );
 serialize( g, "InitUpDownTime", netCDF::NcInt64(), unit.InitUpDownTime );

 serialize( g, "MinPower", netCDF::NcDouble(), unit.MinPower );
 serialize( g, "MaxPower", netCDF::NcDouble(), unit.MaxPower );
 serialize( g, "DeltaRampUp", netCDF::NcDouble(), unit.DeltaRampUp );
 serialize( g, "DeltaRampDown", netCDF::NcDouble(), unit.DeltaRampDown );
 serialize( g, "QuadTerm", netCDF::NcDouble(), unit.QuadTerm );

 if (type == ftDat) {
  auto NumberIntervals = g.getParentGroup().getDim( "NumberIntervals" );
  serialize( g, "LinearTerm", netCDF::NcDouble(), NumberIntervals, b );
  serialize( g, "ConstTerm", netCDF::NcDouble(), NumberIntervals, c );
 } else {
  serialize( g, "LinearTerm", netCDF::NcDouble(), unit.LinearTerm );
  serialize( g, "ConstTerm", netCDF::NcDouble(), unit.ConstTerm );
 }

 serialize( g, "StartUpCost", netCDF::NcDouble(), unit.StartUpCost );
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
 } else {
  inputFile >> mod_file;
  std::cout << mod_file;
 }
 inputFile.close();

 netCDF::NcFile f( filename, netCDF::NcFile::replace );
 f.putAtt( "SMS++_file_type", netCDF::NcInt(), eBlockFile );


 if( type == ftDat ) {
  f.putAtt( "type", "UnitBlock" );
  f.addDim( "TimeHorizon", dat_file.TimeHorizon );
  f.addDim( "NumberIntervals", dat_file.TimeHorizon );

  auto bg = f.addGroup( "Block_0" );
  bg.putAtt( "type", "ThermalUnitBlock" );
  serialize_unit( bg, dat_file.thermal_unit );

 } else {
  f.putAtt( "type", "UnitBlock" );
  f.addDim( "TimeHorizon", mod_file.TimeHorizon );
  f.addDim( "NumberIntervals", 1 );

  for( int i = 0; i < mod_file.NumThermal; ++i ) {
   auto bg = f.addGroup( "Block_" + std::to_string( f.getGroupCount() ) );
   bg.putAtt( "type", "ThermalUnitBlock" );
   mod_file.thermal_units[ i ].generate_startupcost();
   serialize_unit( bg, mod_file.thermal_units[ i ] );
  }
 }
 return 0;
}
