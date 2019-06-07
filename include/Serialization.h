/*--------------------------------------------------------------------------*/
/*-------------------------- File Serialization.h --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Useful functions for serializing and deserializing objects with nectCDF.
 *
 * \version 0.1
 *
 * \date 07 - 06 - 2019
 *
 * \author Rafael Durbano Lobato \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 */

#ifndef __Serialization
#define __Serialization
/* self-identification: #endif at the end of the file */

#include "netcdf"
#include <exception>
#include <numeric>
#include <string>
#include <vector>

namespace SMSpp_di_unipi_it::Serialization {

/*--------------------------------------------------------------------------*/

template<class T>
inline void deserialize_dim( const netCDF::NcGroup & group,
                             const std::string & dim_name, T & data ,
                             bool optional = true ,
                             const std::string & class_name = "" ) {

  netCDF::NcDim ncDim = group.getDim( dim_name );

  if( ncDim.isNull() ) {
    if( optional )
      return;
    throw( std::invalid_argument( class_name + "::deserialize_dim: " +
                                  dim_name + " is not present" ) );
  }

  data = ncDim.getSize();
  if( data <= 0 )
    throw( std::invalid_argument( class_name + "deserialize_dim: " +
                                  dim_name + " must be positive" ) );
}

/*--------------------------------------------------------------------------*/

template<class T>
inline void deserialize( const netCDF::NcGroup & group,
                         const std::string & var_name,
                         std::vector<T> & data,
                         const std::vector<size_t> & sizes,
                         bool optional = true ,
                         const std::string & class_name = "" ) {

  auto total_size = std::accumulate( begin( sizes ), end( sizes ), 1,
                                     std::multiplies<size_t>() );

  if( total_size == 0 ) {
    data.resize( 0 );
    return;
  }

  auto ncVar = group.getVar( var_name );
  if( ncVar.isNull() ) {
    if( optional ) {
      data.resize( 0 );
      return;
    }
    throw( std::invalid_argument( class_name + "::deserialize: " +
                                  var_name + " is not present" ) );
  }

  data.resize( total_size );

  std::vector<size_t> start;
  start.assign( sizes.size(), 0 );

  ncVar.getVar( start, sizes, data.data() );
}

/*--------------------------------------------------------------------------*/

template<class T>
void deserialize( const netCDF::NcGroup & group,
                  const std::string & var_name,
                  T * data, const std::string & class_name = "" ) {

  auto ncVar = group.getVar( var_name );
  if( ncVar.isNull() )
    throw( std::invalid_argument
           ( class_name + "::deserialize: " + var_name + " is not present" ) );
  ncVar.getVar( data );
}

/*--------------------------------------------------------------------------*/

template<class T>
void deserialize( const netCDF::NcGroup & group,
                  const std::string & var_name,
                  const size_t & size, std::vector<T> & data ,
                  const std::string & class_name = "" ) {

  auto ncVar = group.getVar( var_name );
  if( ncVar.isNull() )
    throw( std::invalid_argument
           ( class_name + "::deserialize: " + var_name + " is not present" ) );

  data.resize( size );
  ncVar.getVar( { 0 }, { size }, data.data() );
}

/*--------------------------------------------------------------------------*/

template<class T>
inline void serialize( netCDF::NcGroup & group, const std::string & var_name,
                       const netCDF::NcType & ncType,
                       const std::vector<netCDF::NcDim> & ncDim,
                       const std::vector<T> & data ) {

  std::vector<size_t> start;
  start.assign( ncDim.size(), 0 );

  std::vector<size_t> sizes;
  sizes.resize( ncDim.size() );
  for( size_t i = 0; i < sizes.size(); ++i )
    sizes[i] = ncDim[i].getSize();

  auto total_size = std::accumulate( begin( sizes ), end( sizes ), 1,
                                     std::multiplies<size_t>() );

  if( total_size == 0 )
    return;

  group.addVar( var_name , ncType , ncDim )
    .putVar( start , sizes , data.data() );
}

/*--------------------------------------------------------------------------*/

template<class T>
void serialize( netCDF::NcGroup & group, const std::string & var_name,
                const netCDF::NcType ncType, T data ) {
  ( group.addVar( var_name , ncType ) ).putVar( & data );
}

/*--------------------------------------------------------------------------*/

template<class T>
void serialize( netCDF::NcGroup & group, const std::string & var_name,
                const netCDF::NcType & ncType, const netCDF::NcDim & ncDim,
                const std::vector<T> & data ) {

  group.addVar( var_name , ncType , ncDim )
    .putVar( { 0 } , { data.size() } , data.data() );
}

/*--------------------------------------------------------------------------*/

} // end( namespace SMS_di_unipi_it::Serialization )

#endif
