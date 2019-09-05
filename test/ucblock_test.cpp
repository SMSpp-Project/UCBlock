#include <iostream>
#include <fstream>

#include <UCBlock.h>
#include <ThermalUnitBlock.h>
#include <BusNetworkBlock.h>

using namespace SMSpp_di_unipi_it;

int main( int argc, char ** argv ) {

 std::string filename( argv[ 1 ] );

 netCDF::NcFile f( filename, netCDF::NcFile::read );
 if (f.isNull()) {
  std::cerr << "cannot open nc4 file " << filename << std::endl;
  exit(1);
 }

 netCDF::NcGroupAtt gtype = f.getAtt("SMS++_file_type");
 if (gtype.isNull()) {
  std::cerr << filename << " is not an SMS++ nc4 file" << std::endl;
  exit(1);
 }

 int type;
 gtype.getValues(&type);

 if (type != eBlockFile) {
  std::cerr << filename << " is not an SMS++ nc4 Block file" << std::endl;
  exit(1);
 }

 netCDF::NcGroup bg = f.getGroup("Block_0");
 if (bg.isNull()) {
  std::cerr << "Block_0 empty or undefined in " << filename << std::endl;
  exit(1);
 }

 // UCBlock deserialize and serialize
 auto ucb = dynamic_cast<UCBlock *>(Block::new_Block( "UCBlock" ));
 auto tub = dynamic_cast<ThermalUnitBlock *>(Block::new_Block( "ThermalUnitBlock", ucb ));
 auto bnb = dynamic_cast<BusNetworkBlock *>(Block::new_Block( "BusNetworkBlock", ucb ));
 ucb->deserialize( bg );

 netCDF::NcFile f1( "test.nc4", netCDF::NcFile::replace );
 f1.putAtt( "SMS++_file_type", netCDF::NcInt(), eBlockFile );
 auto bg1 = f1.addGroup( "Block_0" );
 ucb->serialize( bg1 );

 return 0;
}
