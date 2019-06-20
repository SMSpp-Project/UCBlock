/*--------------------------------------------------------------------------*/
/*------------------------- File NetworkData.h -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * This class represents all the network data needed. It has a set of nodes, a
 * set of lines (or arcs), the starting and ending point of the line (however,
 * lines are not oriented), contain the susceptance of line and the minimum
 * and maximum power flow at each line l.
 *
 * \version 0.1
 *
 * \date 20 - 06 - 2019
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Ali Ghezelsoflu \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Antonio Frangioni and Ali Ghezelsoflu
 */

/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __NetworkData
#define __NetworkData /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <vector>
#include "UCBlock.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

///< namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup Network_CLASSES Classes in NetworkData.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS NetworkData -----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// A network data class
/** The NetworkData class represents a network. It is characterized by a
 * set of nodes (NetworkNodes) and a set of lines (or arcs) and all the data
 * needed to construct the network (instead of the demand of each node).
 * */

  class NetworkData : public UCBlock {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

  public:

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

/// constructor, takes the father
/** Constructor of NetworkData, taking possibly a pointer of its
 * father Block.
 *
 */

NetworkData( Block * father_block = nullptr ): UCBlock( father_block ) { }
/*--------------------------------------------------------------------------*/

/// destructor of NetworkData

virtual ~NetworkData() { };

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
 *
 * NetworkData defines a main public type:
 *
 * - Index, the type of indices;
 @{ */

typedef std::size_t Index;                 ///< index of parameters

/*@} -----------------------------------------------------------------------*/
/*------------- METHODS FOR READING THE DATA OF THE NetworkData ------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the NetworkData
    @{ */

 static NetworkData * new_NetworkData( netCDF::NcGroup & group )
 {
  if( /* the dimension NumberNodes is *not* there or it is == 1 */ )
   return( nullptr );

  auto nd = new NetworkData();
  nd->deserialize( group );
  return( nd );
  }

/*@} -----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

/// extends Block::deserialize( netCDF::NcGroup )
/** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
 * the NetworkData. Besides the mandatory "type" attribute of any :Block, the
 * group should contain the following:
 *
 * - the dimension "NumberNodes" containing the number of nodes in the
 *   problem; this dimension is optional, if it is not provided then it is
 *   taken to be == 1;
 *
 * - the dimension "NumberLines" containing the number of arcs in the problem;
 *   if NumberNodes == 1 then this dimension need not to be present since it
 *   is not loaded;
 *
 * - the variable "StartLine", of type int and indexed over the dimension
 *   "NumberNodes"; the i-th entry of the variable is the starting point of
 *   the line (however, lines are not oriented)
 *
 * - the variable "EndLine", of type int and indexed over the dimension
 *   "NumberNodes"; the i-th entry of the variable is the ending point of the
 *   line (however, lines are not oriented)
 *
 * - the variable "MinPowerFlow", of type double and indexed over the
 *   dimension "NumberLines"; the i-th entry of the variable is assumed to
 *   contain the minimum power flow at line i; if NumberNodes == 1 then this
 *   variable need not to be present since it is not loaded;
 *
 * - the variable "MaxPowerFlow", of type double and indexed over the
 *   dimension "NumberLines"; the i-th entry of the variable is assumed to
 *   contain the maximum power flow at line i; if NumberNodes == 1 then this
 *   variable need not to be present since it is not loaded;
 *
 * - the variable "Susceptance", of type double and indexed over the
 *   "NumberLines"; the i-th entry of this variable is assumed to contain the
 *   susceptance of line i; if NumberNodes == 1 then this variable need not to
 *   be present since it is not loaded.
 */

virtual void deserialize( netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/


virtual void load( std::istream &input ) override {
  throw( std::logic_error( "NetworkData::load() not implemented yet" ) );
};

/*@} -----------------------------------------------------------------------*/
/*--------------------- METHODS FOR SAVING THE NetworkBlock ----------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the NetworkBlock
 *  @{ */

/// extends Block::serialize( netCDF::NcGroup )
/** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
 * NetworkBlock. See NetworkBlock::deserialize( netCDF::NcGroup ) for
 * details of the format of the created netCDF group.
 */

virtual void serialize( netCDF::NcGroup & group ) const override;

/*@} -----------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

  protected:

/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/

SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

/// number of nodes of the network
Index f_number_nodes;

/// number of lines of the network
Index f_number_lines;

/// set starting lines
std::vector< int > v_start_line;

/// set ending lines
std::vector< int > v_end_line;

/// vector to store the susceptance of each line of the network
std::vector< double > v_susceptance;

/// vector to store the minimum power flow at each line
std::vector< double > v_min_power_flow;

/// vector to store the maximum power flow at each line
std::vector< double > v_max_power_flow;



  };   // end( class( NetworkData) )

}  /* namespace SMSpp_di_unipi_it */

#endif /* NetworkData.h included */

/*--------------------------------------------------------------------------*/
/*------------------------ End File NetworkData.h --------------------------*/
/*--------------------------------------------------------------------------*/
