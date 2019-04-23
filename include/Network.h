/*--------------------------------------------------------------------------*/
/*--------------------------- File Network.h -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * This class represents a network. It has a set of nodes
 * (NetworkNodes) and a set of lines (or arcs). A line is defined as a
 * pair of NetworkNodes, which means that those NetworkNodes are
 * connected in the network.
 *
 * \version 0.1
 *
 * \date 04 - 03 - 2019
 *
 * \author Ali Ghezelsoflu \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Rafael Durbano Lobato \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Ali Ghezelsoflu and Rafael Durbano Lobato
 */

/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __Network
#define __Network
/* self-identification: #endif at the end of the file */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <vector>

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

///< namespace for the Structured Modeling System++ (SMS++)
namespace SMSpp_di_unipi_it
{

class NetworkNode;  // forward definition of NetworkNode

/*--------------------------------------------------------------------------*/
/*------------------------------- CLASSES ----------------------------------*/
/*--------------------------------------------------------------------------*/
/** @defgroup Network_CLASSES Classes in Network.h
 *  @{ */

/*--------------------------------------------------------------------------*/
/*---------------------------- CLASS Network -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// A network
/** The Network class represents a network. It is characterized by a
 * set of nodes (NetworkNodes) and a set of lines (or arcs). The lines
 * are represented by a vector of pairs, each pair being formed by two
 * (pointers to) NetworkNodes that correspond to connected nodes in
 * the network. */

class Network {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*----------------- METHODS FOR MODIFYING THE Network ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the Network
 *  @{ */

  /// adds a NetworkNode to this Network
  /** Adds a single (pointer to a) NetworkNode to this Network. */

  void add_network_node( NetworkNode * node ) {
    v_nodes.push_back( node );
  }

  /// adds a line to this Network
  /** Adds a line to this Network. A line (or arc) is defined as a
   * pair of (pointers to) NetworkNodes. */

  void add_line( std::pair<NetworkNode *, NetworkNode *> line ) {
    v_lines.push_back( line );
  }

/*@} -----------------------------------------------------------------------*/
/*------------- METHODS FOR READING THE DATA OF THE Network ----------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the Network
    @{ */

  /// returns the number of nodes in this Network

  int get_num_nodes() const {
    return v_nodes.size();
  }

  /// returns the number of lines in this Network

  int get_num_lines() const {
    return v_lines.size();
  }

  /// returns the vector of NetworkNodes

  const std::vector<NetworkNode *> & get_nodes() const {
    return v_nodes;
  }

  /// returns the i-th NetworkNode

  NetworkNode * get_node( int i ) const {
    return v_nodes[i];
  }

/*@} -----------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------------- PROTECTED FIELDS  ----------------------------*/
/*--------------------------------------------------------------------------*/

  /// vector of nodes (NetworkNodes) of this Network
  std::vector<NetworkNode *> v_nodes;

  /// lines (arcs) of this Network
  std::vector<std::pair<NetworkNode *, NetworkNode *>> v_lines;

};   // end( class( Network ) )

/*@}  end( group( Network_CLASSES ) ) */

} /* namespace SMSpp_di_unipi_it */

#endif /* Network.h included */

/*--------------------------------------------------------------------------*/
/*-------------------------- End File Network.h ----------------------------*/
/*--------------------------------------------------------------------------*/
