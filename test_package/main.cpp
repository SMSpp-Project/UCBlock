#include <fstream>

#include <UCBlock.h>

using namespace SMSpp_di_unipi_it;

int main( int argc, char ** argv ) {

 auto filename = std::string( argv[ 1 ]);
 netCDF::NcFile f;
 f.open( filename, netCDF::NcFile::read );
 netCDF::NcGroup bg = f.getGroup( "Block_0" );
 auto ucb = static_cast< UCBlock * >(Block::new_Block( "UCBlock" ));
 ucb->deserialize( bg );

 return( 0 );
}
