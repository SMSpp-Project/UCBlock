/*--------------------------------------------------------------------------*/
/*--------------------- File DCNetworkBlock.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the DCNetworkBlock class.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Ali Ghezelsoflu \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Wim van Ackooij \n
 *         EDF R&D OSIRIS \n
 *
 * \author Quentin Jacquet \n
 *         EDF R&D OSIRIS \n
 *
 * \author Rafael Durbano Lobato \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Ali Ghezelsoflu, Wim van Ackooij,
 *                      Quentin Jacquet, Rafael Durbano Lobato
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <map>

#include <queue>

#include <utility>

#include "NetworkBlock.h"

#include "DCNetworkBlock.h"

#include "LinearFunction.h"

#include "OneVarConstraint.h"

#include "FRealObjective.h"

#include <Eigen/Sparse>

#include <Eigen/SparseLU>

#include <chrono>

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

using SpMat = Eigen::SparseMatrix< double >;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register DCNetworkBlock to the Block factory
SMSpp_insert_in_factory_cpp_1( DCNetworkBlock );

// register DCNetworkBlockSolution to the Solution factory
SMSpp_insert_in_factory_cpp_0( DCNetworkBlockSolution );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
// register DCNetworkData to the NetworkData factory

using DCNetworkData = DCNetworkBlock::DCNetworkData;

SMSpp_insert_in_factory_cpp_0( DCNetworkData );

/*--------------------------------------------------------------------------*/
/*---------------------------- STATIC FUNCTIONS ----------------------------*/
/*--------------------------------------------------------------------------*/
/* Deserializes a 2-dim double multi_array where the first dimension is
 * optional: if the ncVar only have one dimension it is taken to be the
 * second one and a column-wise matrix (1 row, dim2 columns) is returned.
 * If, furthermore, the variable is scalar (both dimensions are 1) then
 * the column-wise matrix contains that single value repeated over. */

bool deserialize_opt( const netCDF::NcGroup & group ,
		      const std::string & name ,
		      std::size_t dim1 , std::size_t dim2 ,
		      boost::multi_array< double , 2 > & multi_array )
{
 using index = typename boost::multi_array< double , 2 >::index;
 static const std::vector< index > empty = { 0 , 0 };

 auto ncVar = group.getVar( name );
 if( ncVar.isNull() ) {
  multi_array.resize( empty );
  return( false );
  }

 if( ncVar.getDimCount() == 0 )
  throw( std::invalid_argument( "netCDF variable " + name + "is empty" ) );

 if( ncVar.getDimCount() > 2 )
  throw( std::invalid_argument( "netCDF variable " + name +
				"has too many dimensions" ) );

 std::vector< index > size = { dim1 , dim2 };
 if( ncVar.getDimCount() == 1 )  // 1-dimensional variable
  size[ 0 ] = 1;                 // ignore the first dimension

 multi_array.resize( size );

 if( size[ 0 ] == 1 ) {          // 1-dimensional variable
  if( (ncVar.getDims())[ 0 ] == 1 ) {  // in fact, a scalar variable
   ncVar.getVar( { 0 } , { 1 } , multi_array.data() );
   for( Index i = 1 ; i < multi_array.data().size() ; ++i )
    multi_array.data()[ i ] = multi_array.data()[ 0 ];  // replicate
   }
  else                                 // a vector
   ncVar.getVar( { 0 } , { dim2 } , multi_array.data() );
  }
 else                            // 2-dimensional variable
  ncVar.getVar( empty , size , multi_array.data() );

 return( true );
 }

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS OF DCNetworkData -------------------------*/
/*--------------------------------------------------------------------------*/

void DCNetworkData::deserialize( const netCDF::NcGroup & group )
{
 NetworkData::deserialize( group );

 if( ! deserialize_dim( group , "ReferenceNode" , f_reference_node , true ) )
  f_reference_node = 0;

 // Optional variables

 if( f_number_nodes == 1 )
  return;

 deserialize_dim( group , "NumberLines" , f_number_lines , false );

 if( ! deserialize_dim( group , "NumberBranches" , f_number_branches ,
			true ) )
  f_number_branches = f_number_lines;

 Index time_instants = 1;
 deserialize_dim( group , "NumberInstants" , time_instants , true );
 
 ::deserialize( group , "StartLine" , f_number_branches , v_start_line ,
                false , false );

 ::deserialize( group , "EndLine" , f_number_branches , v_end_line ,
                false , false );

 for( Index i = 0 ; i < f_number_branches ; ++i ) {
  if( ( v_start_line[ i ] < 0 ) || ( v_start_line[ i ] >= f_number_nodes ) )
   throw( std::invalid_argument( "DCNetworkData::deserialize: "
				 "wrong start node number " +
				 std::to_string( v_start_line[ i ] ) ) );

  if( ( v_end_line[ i ] < 0 ) || ( v_end_line[ i ] >= f_number_nodes ) )
   throw( std::invalid_argument( "DCNetworkData::deserialize: "
				 "wrong end node number " +
				 std::to_string( v_end_line[ i ] ) ) );

  if( v_start_line[ i ] == v_end_line[ i ] )
   throw( std::invalid_argument( "DCNetworkData::deserialize: "
				 "start node == end node for branch " +
				 std::to_string( v_end_line[ i ] ) ) );
  }

 deserialize_opt( group , "MaxPowerFlow" , time_instants ,
		  f_number_lines , v_max_power_flow );

 deserialize_opt( group , "MinPowerFlow" , time_instants ,
		  f_number_lines , v_min_power_flow );

 if( ! deserialize_opt( group , "Efficiency" , time_instants ,
			f_number_branches , v_efficiency ) ) {
  using index = boost::multi_array< double , 2 >::index;
  const std::vector< index > size = { 1 , f_number_branches };
  v_efficiency.resize( size );
  for( Index b = 0 ; b < f_number_branches ; ++b )
   v_efficiency[ 0 ][ b ] = 1;
  }

 ::deserialize( group , "NetworkCost" , f_number_lines , v_network_cost ,
                true , true );

 ::deserialize( group , "LineSusceptance" , f_number_lines ,
		v_line_susceptance , true , true );

 f_number_HVDC_lines = std::count_if( v_line_susceptance.begin() ,
				      v_line_susceptance.end() ,
				      []( auto s ) { return( s == 0 ); }
				      );
 if( f_number_HVDC_lines == f_number_lines )
  v_line_susceptance.clear();

 //!! TODO: in the hypergraph case the line names are garbled!!
 ::deserialize( group , "LineName" , f_number_lines , v_line_names );

 if( is_hypergraph() ) {
  std::vector< Index > id;
  ::deserialize( group , "HyperArcID" , f_number_branches , id ,
                 false , false );

  // < start node , end node , branch name >
  std::vector< std::vector< std::tuple< Index , Index , Index > > >
                                                      tmp( f_number_lines );
  for( Index i = 0 ; i < f_number_branches ; ++i ) {
   if( ( id[ i ] < 0 ) || ( id[ i ] >= f_number_lines ) )
    throw( std::invalid_argument( "DCNetworkData::deserialize: "
				  "wrong hyperarc id " +
				  std::to_string( id[ i ] ) ) );
   tmp[ id[ i ] ].push_back( std::make_tuple( v_start_line[ i ] ,
                                              v_end_line[ i ] , i ) );
   if( std::get< 0 >( tmp[ id[ i ] ].front() ) !=
       std::get< 0 >( tmp[ id[ i ] ].back() ) )
    throw( std::invalid_argument( "DCNetworkData::deserialize: "
				  "branches for line " +
				  std::to_string( id[ i ] ) +
				  " have different start" ) );
   }

  // sort all tmp[ i ] for increasing end node
  for( Index i = 0 ; i < f_number_lines ; ++i ) {
   if( tmp[ i ].empty() )
    throw( std::invalid_argument( "DCNetworkData::deserialize: "
				  "no branches for line " +
				  std::to_string( i ) ) );
   std::sort( tmp[ i ].begin() , tmp[ i ].end() ,
              []( auto & a , auto & b ) {
               return( std::get< 1 >( a ) < std::get< 1 >( b ) ); } );
   for( Index j = 1 ; j < tmp[ i ].size() ; ++j )
    if( std::get< 1 >( tmp[ i ][ j ] ) == std::get< 1 >( tmp[ i ][ j - 1 ] ) )
     throw( std::invalid_argument( "DCNetworkData::deserialize: "
				   "repeated end node for line " +
				   std::to_string( i ) ) );
   }

  // resize v_start_line and put there the right start nodes
  v_start_line.resize( f_number_lines );
  for( Index i = 0 ; i < f_number_lines ; ++i )
   v_start_line[ i ] = std::get< 0 >( tmp[ i ].front() );

  // clear v_end_line
  v_end_line.clear();

  // build v_end_lines and v_h_efficiency
  v_end_lines.resize( f_number_lines );
  time_instants = v_h_efficiency.shape()[ 0 ];
  v_h_efficiency.resize( time_instants );
  for( Index t = 0 ; t < time_instants ; ++t )
   v_h_efficiency[ t ].resize( f_number_lines );

  for( Index i = 0 ; i < f_number_lines ; ++i ) {
   if( ( ! v_line_susceptance.empty() ) &&
       ( v_line_susceptance[ i ] != 0.0 ) && ( tmp[ i ].size() > 1 ) )
    throw( std::invalid_argument( "DCNetworkData::deserialize: line " +
				  std::to_string( id[ i ] ) +
				  " is a hyperarc but has susceptance" ) );
   v_end_lines[ i ].resize( tmp[ i ].size() );
   for( Index j = 0 ; j < tmp[ i ].size() ; ++j )
    v_end_lines[ i ][ j ] = std::get< 1 >( tmp[ i ][ j ] );

   for( Index t = 0 ; t < time_instants ; ++t ) {
    v_h_efficiency[ t ][ i ].resize( tmp[ i ].size() );
    for( Index j = 0 ; j < tmp[ i ].size() ; ++j )
     v_h_efficiency[ t ][ i ][ j ] =
                        v_efficiency[ t ][ std::get< 2 >( tmp[ i ][ j ] ) ];
    }
   }

  // clear v_efficiency
  //!! using index = boost::multi_array< double , 2 >::index;
  //!! const std::vector< index > empty = { 0 , 0 };
  //!! v_efficiency.resize( empty );
  v_efficiency.resize( { 0 , 0 } );

  }  // end( if( is_hypergraph() ) )

 stored_B2 = SpMat( f_number_nodes - 1 , f_number_nodes - 1 );
 stored_B2_inv = SpMat( f_number_nodes - 1 , f_number_nodes - 1 );

 ::deserialize( group , "NodeName" , f_number_nodes , v_node_names );

 }  // end( DCNetworkData::deserialize )

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG

std::vector< std::string > DCNetworkData::expected_dims( void ) const {
 static const std::vector< std::string > ed =
 { "NumberLines" , "NumberBranches" , "ReferenceNode" };

 auto ret = NetworkData::expected_dims();
 ret.insert( ret.end() , ed.begin() , ed.end() );

 return( ret );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

std::vector< std::string > DCNetworkData::expected_vars( void ) const {
 static const std::vector< std::string > ev =
 { "StartLine" , "EndLine" , "HyperArcID" , "MinPowerFlow" , "MaxPowerFlow" ,
   "LineSusceptance" , "NetworkCost" , "LineName" , "Efficiency" };

 auto ret = NetworkData::expected_vars();
 ret.insert( ret.end() , ev.begin() , ev.end() );

 return( ret );
 }

#endif

/*--------------------------------------------------------------------------*/

int DCNetworkData::get_reducedIdx( int idx ) const {
 if( nb_components == 1 ) {
  const int ref = get_reference_node();
  if( idx == ref )
   return( -1 );
  return( idx > ref ? idx - 1 : idx );
  }
 return( v_reduced_idx[ idx ] );
 }

/*--------------------------------------------------------------------------*/

int DCNetworkData::get_originalIdx( int idx ) const {
 if( nb_components == 1 ) {
  const int ref = get_reference_node();
  return( idx >= ref ? idx + 1 : idx );
  }
 // In the multi-component case, we use the explicit inverse mapping
 // v_original_idx[red_idx] = original_node
 return( v_original_idx[ idx ] );
 }

/*--------------------------------------------------------------------------*/

void DCNetworkData::identify_connected_components( void )
{
 Index nb_nodes = get_number_nodes();
 const auto & DC_lines = get_DC_lines();
 const auto & start_line = get_start_line();
 const auto & end_line = get_end_line();
 std::vector< std::vector< Index > > adj( nb_nodes );

 for( auto line_id : DC_lines ) {
  Index i = start_line[ line_id ];
  Index j = end_line[ line_id ];
  adj[ i ].push_back( j );
  adj[ j ].push_back( i );
 }
 v_component.assign( nb_nodes , -1 );
 nb_components = 0;

 for( Index v = 0 ; v < nb_nodes ; ++v ) {
  if( v_component[ v ] != -1 ) continue;

  std::queue< Index > q;
  q.push( v );
  v_component[ v ] = nb_components;

  while( ! q.empty() ) {
   Index u = q.front();
   q.pop();
   for( Index w : adj[ u ] ) {
    if( v_component[ w ] == -1 ) {
     v_component[ w ] = nb_components;
     q.push( w );
    }
   }
  }

  ++nb_components;
 }

 // Now simply identify each collection of nodes
 v_nodes_in_component.resize( nb_components );
 for( Index v = 0 ; v < nb_nodes ; ++v )
  v_nodes_in_component[ v_component[ v ] ].push_back( v );
 }

/*--------------------------------------------------------------------------*/

void DCNetworkData::compute_DCDF( c_Subset & HVDC_lines ,
                                  const SpMat & PTDF_matrix )
{
 const auto number_nodes = get_number_nodes();
 const auto number_lines = get_number_lines();
 const auto & start_line = get_start_line();
 const auto & end_line = get_end_line();
 
 // Connected components must already be known
 const size_t nb_c = get_nb_connected_components();

 // linking constraints between DC and HVDC
 SpMat A_DC_transpose( number_nodes - nb_c , number_lines );
 double eta = 1.0;
 for( auto & line_id : HVDC_lines ) {
  // Start node contribution (+1)
  Index s = start_line[ line_id ];
  int rs = get_reducedIdx( s );
  if( rs >= 0 )
   A_DC_transpose.coeffRef( rs , line_id ) = 1.0;
  if( ! is_hypergraph() ) { // no hypergraph
   eta = get_line_efficiency( line_id ); // efficiency of the HVDC line
   Index e = end_line[ line_id ];
   int re = get_reducedIdx( e );
   if( re >= 0 )
    A_DC_transpose.coeffRef( re , line_id ) = -eta;
   }
  else { // with hypergraph
   for( Index i = 0 ; i < get_end_lines()[ line_id ].size() ; ++i ) {
    eta = get_line_efficiencies( line_id )[ i ]; // efficiency of the hyperarc
    Index e = get_end_lines()[ line_id ][ i ];
    int re = get_reducedIdx( e );
    if( re >= 0 )
     A_DC_transpose.coeffRef( re , line_id ) = -eta;
    }
   }
  }
 DCDF = -PTDF_matrix * A_DC_transpose;
 DCDF_was_computed = true;
 }

/*--------------------------------------------------------------------------*/

SpMat DCNetworkData::get_PTDF( c_Subset & DC_lines , double tikhonov_coeff )
{
 // get data
 const auto number_nodes = get_number_nodes();
 const auto number_lines = get_number_lines();
 if( number_lines <= 0 )
  throw( std::logic_error( "DCNetworkBlock::generate_abstract_constraints: "
    				                   "number of lines of DCNetworkBlock is not set" ) );

 const auto & start_line = get_start_line();
 const auto & end_line = get_end_line();

 // Compute any connected components now
 identify_connected_components();

 const size_t nb_c = get_nb_connected_components();
 const auto & v_comp = get_subgraphs();
 /*-- some possible printing for debugging if need be
  std::cout << " Found " << nb_c << " connected components\n";
 for( int i = 0 ; i < nb_c ; ++i ) {
  std::cout << " Component " << i << " has nodes = ";
  for( int j = 0 ; j < v_comp[ i ].size() ; ++j )
   std::cout << v_comp[ i ][ j ] << " ; ";
  std::cout << "\n";
 } */

 // construct the matrix using two sub-matrices B_bar and B_hat
 // note that v_line_susceptance may be empty but it is only used inside
 // the for( ... : DC_lines ), which do not execute if it is
 SpMat B_hat = SpMat( number_lines , number_nodes );
 for( auto line_id : DC_lines ) {
  auto s = v_line_susceptance[ line_id ];
  B_hat.insert( line_id , start_line[ line_id ] ) = s;
  B_hat.insert( line_id , end_line[ line_id ] ) = - s;
  }

 SpMat B_bar = SpMat( number_nodes , number_nodes );
 std::vector< double > B_bar_diag( number_nodes , 0.0 );
 for( auto line_id : DC_lines ) {
  auto s = v_line_susceptance[ line_id ];
  B_bar.insert( start_line[ line_id ] , end_line[ line_id ] ) = - s;
  B_bar.insert( end_line[ line_id ] , start_line[ line_id ] ) = - s;
  B_bar_diag[ start_line[ line_id ] ] += s;
  B_bar_diag[ end_line[ line_id ] ] += s;
  }

 for( Index node_id = 0 ; node_id < number_nodes ; ++node_id )
  B_bar.insert( node_id , node_id ) = B_bar_diag[ node_id ] + tikhonov_coeff;

 /*
 * When only a single connected component is there, the situation is simple
 *  we pick the reference node and form the Identity matrix from which we delete that row
 * When multiple connected components are there, potentially just having a single node altogether (isolated node)
 *  then the situation consists of picking a reference node per "AC block"
 *  this too is fairly simple and consists of deleting as many rows as need be 
 *  in fact one per connected component
 * In the latter case we will not (at least for the time being) allow the user to pick a reference node per component 
 *  and will do this for him. 
 * Should we opt to remove this freedom from the user alltogether then the code can be simplified
 */ 
 SpMat I_nref = SpMat( number_nodes , number_nodes - nb_c );
 if( nb_c == 1 ) {
  Index ref_node = get_reference_node();
  for( Index node_id = 0 ; node_id < ref_node ; ++node_id )
    I_nref.insert( node_id , node_id ) = 1.;
  for( Index node_id = ref_node ; node_id < number_nodes - 1 ; ++node_id )
    I_nref.insert( node_id + 1 , node_id ) = 1.;
 }
 else {
  // Identification of the reference nodes
  v_reduced_idx.resize( number_nodes , -1 );
  v_original_idx.resize( number_nodes - nb_c , -1 );
  std::vector< bool > is_ref( number_nodes , false );
  for( const auto & comp : v_comp )
   is_ref[ comp.front() ] = true;

  int col = 0;
  for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {
   if( ! is_ref[ node_id ] ) {
    I_nref.insert( node_id , col ) = 1.0;
    v_reduced_idx[ node_id ] = col;
    v_original_idx[ col ] = node_id;
    ++col;
   }
  }
 }
 I_nref.makeCompressed();

 // construction of B1 and B2 (see documentation)
 SpMat B1 = B_hat * I_nref;
 SpMat B2 = I_nref.transpose() * B_bar * I_nref;

 // compute inverse of B2 and deduce PTDF
 SpMat PTDF_matrix;
 SpMat B2_inv = SpMat( number_nodes - nb_c , number_nodes - nb_c );
 std::pair< SpMat , SpMat > t = get_stored_B2();
 if( (t.first.rows() == number_nodes - nb_c) && B2.isApprox( t.first ) )
  B2_inv = t.second;
 else {
  // Inversion of sparse matrix with eigen (solve B2*X = I)
  //Eigen::BiCGSTAB<SpMat> solver;
  Eigen::SparseLU< SpMat > solver;
  solver.compute( B2 );
  SpMat I( number_nodes - nb_c , number_nodes - nb_c );
  I.setIdentity();
  if( solver.info() != Eigen::Success )
   std::cout << "Inversion in PTDF not possible" << std::endl;
  else {
   B2_inv = solver.solve( I );
   set_stored_B2( B2 , B2_inv );
   }
  }

 PTDF_matrix = B1 * B2_inv;

 return( PTDF_matrix );
 }

/*--------------------------------------------------------------------------*/


/**
 * Compute a fundamental cycle basis of the DC network using
 * a depth‑first search (DFS)–based variant of Paton's algorithm.
 *
 * The algorithm operates on the DC subgraph only:
 *  - Nodes represent buses.
 *  - Edges represent DC transmission lines.
 *  - HVDC lines are explicitly excluded from the topology.
 *
 * Overview of the algorithm:
 * --------------------------
 * 1. Build an undirected adjacency list of the DC graph.
 *
 * 2. Traverse each connected component of the DC graph using
 *    an explicit (iterative) depth‑first search.
 *
 * 3. During DFS, build a spanning forest:
 *      - Each node stores its parent in the DFS tree.
 *      - Roots satisfy parent[root] == root.
 *
 * 4. Track the current DFS path explicitly using an `in_stack` flag.
 *    This is crucial to distinguish true back‑edges to ancestors
 *    from cross‑edges when using an iterative DFS.
 *
 * 5. Whenever a back‑edge (u → v) is encountered such that:
 *      - v is already on the current DFS path (v is an ancestor of u),
 *      - v is not the parent of u,
 *    a fundamental cycle is detected.
 *
 * 6. Construct the cycle by:
 *      - Starting from the back‑edge endpoint v,
 *      - Walking up the parent pointers from u until v is reached,
 *      - Closing the cycle at v.
 *
 *    The resulting cycle is stored as an ordered list of nodes
 *    with cycle.front() == cycle.back().
 *
 * 7. Repeat until all connected components have been explored.
 *
 * Properties of the computed cycle basis:
 * ---------------------------------------
 * - Each cycle corresponds to exactly one non‑tree (back) edge.
 * - The set of cycles forms a fundamental cycle basis:
 *      |cycles| = |E_DC| − |V_DC| + (number of DC connected components)
 * - Cycles are expressed in node form here and later converted
 *   to edge / line incidence form in a separate routine.
 *
 * Design notes:
 * -------------
 * - An explicit stack is used instead of recursion to avoid
 *   stack overflows on large networks.
 *
 * - The `in_stack` array replaces the recursion stack marker
 *   normally used in recursive DFS and is required for correctness.
 *
 * - Self‑loops are handled explicitly and added as trivial cycles.
 *
 * - The algorithm is root‑agnostic: each DC connected component
 *   has its own DFS root.
 *
 * References:
 * -----------
 * - K. Paton, "An algorithm for finding a fundamental set of
 *   cycles of a graph", Communications of the ACM, 1969.
 * - NetworkX implementation of cycle_basis (adapted to C++).
 */
void DCNetworkData::compute_cycle_basis( void ) {
 if( cycle_basis_was_computed )
  return;

 const Index number_nodes = get_number_nodes();
 const Index number_lines = get_number_lines();

 if( number_lines <= 0 )
  throw std::logic_error(
   "DCNetworkData::compute_cycle_basis: number of lines not set"
  );

 const auto & start_line = get_start_line();
 const auto & end_line = get_end_line();

 /*------------------------------------------------------------
  * Build adjacency list (undirected)
  *------------------------------------------------------------*/
 std::vector< std::vector< Index > > neighbors( number_nodes );

 v_cycle_basis.clear();
 const auto & DC_lines = get_DC_lines();
 for( Index line_id : DC_lines ) {
  Index u = start_line[ line_id ];
  Index v = end_line[ line_id ];

  if( u == v ) {
   // self-loop → trivial cycle
   v_cycle_basis.push_back( { u } );
   continue;
  }

  neighbors[ u ].push_back( v );
  neighbors[ v ].push_back( u );
 }

 /*------------------------------------------------------------
  * Initialise DFS / Paton state
  *------------------------------------------------------------*/

 m_spanning_parent.assign( number_nodes , -1 );

 std::vector< int > depth( number_nodes , -1 );

 // marks nodes on the current DFS path
 std::vector< bool > in_stack( number_nodes , false );

 // explicit DFS stack
 std::vector< Index > stack;

 /*------------------------------------------------------------
  * DFS over connected components
  *------------------------------------------------------------*/
 for( Index start = 0 ; start < number_nodes ; ++start ) {
  if( depth[ start ] != -1 || neighbors[ start ].empty() )
   continue;

  // new connected component
  depth[ start ] = 0;
  m_spanning_parent[ start ] = start;
  stack.push_back( start );
  in_stack[ start ] = true;

  while( ! stack.empty() ) {
   Index u = stack.back();
   stack.pop_back();

   bool pushed_child = false;
   for( Index v : neighbors[ u ] ) {
    // Tree edge
    if( depth[ v ] == -1 ) {
     depth[ v ] = depth[ u ] + 1;
     m_spanning_parent[ v ] = static_cast< int >( u );

     stack.push_back( u ); // resume u later
     stack.push_back( v ); // DFS into v
     in_stack[ v ] = true;

     pushed_child = true;
     break; // important: depth-first
    }

    // -------------------------------------------------
    // Back edge to ANCESTOR → fundamental cycle
    // -------------------------------------------------
    if( v != static_cast< Index >( m_spanning_parent[ u ] ) &&
     in_stack[ v ] ) {
     Subset cycle;
     cycle.push_back( v );

     int x = static_cast< int >( u );
#ifndef NDEBUG
     int guard = 0;
#endif
     while( x != static_cast< int >( v ) ) {
      cycle.push_back( static_cast< Index >( x ) );
      x = m_spanning_parent[ x ];
#ifndef NDEBUG
      // Safety guards
      assert( x >= 0 && x < static_cast< int >( number_nodes ) );
      assert( ++guard <= static_cast< int >( number_nodes ) );
#endif
     }

     cycle.push_back( v ); // close cycle
     v_cycle_basis.push_back( std::move( cycle ) );
    }
   }

   // -------------------------------------------------
   // Finished exploring u
   // -------------------------------------------------
   if( ! pushed_child ) {
    in_stack[ u ] = false;
   }
  }
 }
 cycle_basis_was_computed = true;
}

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS OF DCNetworkBlock ------------------------*/
/*--------------------------------------------------------------------------*/

DCNetworkBlock::~DCNetworkBlock()
{
 Constraint::clear( v_power_flow_injection_const );
 Constraint::clear( v_DC_HVDC_power_flow_const );
 Constraint::clear( v_power_flow_relax_abs );
 Constraint::clear( v_power_flow_def);
 Constraint::clear( v_power_flow_limit_const );
 Constraint::clear( v_power_flow_limit_design_const );
 overall_balanced_const.clear();
 Constraint::clear( v_CYCLE_def_flow_const );
 Constraint::clear( v_CYCLE_def_cycle_const );
 Constraint::clear( v_KIRCHHOFF_power_flow_def );
 Constraint::clear( v_KIRCHHOFF_node_balance_const );
 v_reference_angle_const.clear();

 objective.clear();

 // delete the DCNetworkData if it is local
 if( f_local_NetworkData )
  delete( f_NetworkData );
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void DCNetworkBlock::deserialize( const netCDF::NcGroup & group )
{
 // deal with the DCNetworkData first

 Index NumberNodes;
 if( deserialize_dim( group , "NumberNodes" , NumberNodes ) ) {
  // since the dimension "NumberNodes" has been provided, a DCNetworkData
  // is there: deserialize is and mark it as being local
  if( f_local_NetworkData ) {   // if the NetworkData has not been passed
   delete( f_NetworkData );     // from UCBlock, then delete it
   f_NetworkData = nullptr;
   }
  auto DCND = get_new_NetworkData();
  DCND->deserialize( group );
  if( f_NetworkData &&
      ( f_NetworkData->get_number_nodes() != DCND->get_number_nodes() ) )
   throw( std::logic_error( "DCNetworkBlock::deserialize: NumberNodes "
    				                    "not matching between NetworkData" ) );
  f_NetworkData = DCND;
  f_local_NetworkData = true;
  // a [DC]NetworkData has been provided, so the size of the given vector
  // of active demand must be equal to the number of nodes.
  ::deserialize( group , "ActiveDemand" , NumberNodes , v_ActiveDemand );
  }
 else {
  // a DCNetworkData has not been provided, but the active demand may still
  // have been provided

  auto ActiveDemand = group.getVar( "ActiveDemand" );

  if( ! ActiveDemand.isNull() ) {
   // the active demand has indeed been provided.

   if( ActiveDemand.getDimCount() != 1 )
    // the active demand must be a one-dimensional array
    throw( std::invalid_argument(
     "DCNetworkBlock::deserialize(): ActiveDemand should have one dimension, "
     "but it has " + std::to_string( ActiveDemand.getDimCount() ) ) );

   // retrieve the number of nodes from the size of the given netCDF variable
   NumberNodes = ActiveDemand.getDim( 0 ).getSize();

   // resize the vector of active demand.
   v_ActiveDemand.resize( NumberNodes );

   // retrieve the active demand from the netCDF variable.
   ActiveDemand.getVar( v_ActiveDemand.data() );
   }
  }

 // note: NetworkBlock::deserialize() is called after dealing with the
 //       DCNetworkData, so that it's there when the list of expected stuff
 //       is constructed and checked
 NetworkBlock::deserialize( group );

 }  // end( DCNetworkBlock::deserialize )

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG

std::vector< std::string > DCNetworkBlock::expected_vars( void ) const {
 auto ret = NetworkBlock::expected_vars();
 ret.push_back( "ActiveDemand" );
 ret.push_back( "Kappa" );

 return( ret );
 }

#endif

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::generate_abstract_variables( Configuration * stvv )
{
 if( variables_generated() )  // variables have already been generated
  return;                     // nothing to do

 // call the method of the base class to generate node injection variables
 NetworkBlock::generate_abstract_variables( stvv );

 // read Configuration to set the formulation
 Index wf = 0;  // 0: PTDF (default), 1: cycle, 2: Kirchhoff
 if( ( ! stvv ) && f_BlockConfig )
  stvv = f_BlockConfig->f_static_variables_Configuration;
 if( auto sci = dynamic_cast< SimpleConfiguration< int > * >( stvv ) )
  wf = sci->f_value;

 // flow and (if necessary) auxiliary cost variables are always there
 generate_PTDF_variables();
 
 switch( wf ) {
  case( 1 ): ftype = CYCLE;     generate_CYCLE_variables(); break;
  case( 2 ): ftype = KIRCHHOFF; generate_KIRCHHOFF_variables(); break;
  default:   ftype = PTDF; 
  }  

 set_variables_generated();

 }  // end( DCNetworkBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::generate_PTDF_variables( void )
{
 /** This formulation corresponds to the "PTDF + FLOW" formulation of
  * "Linear Optimal Power Flow Using Cycle Flows" of
  *  Jonas Horsch, Henrik Ronellenfitsch, Dirk Witthaut, Tom Brown */

 auto number_lines = get_number_lines();
 if( number_lines <= 0 )
  return;

 // the power flow Variable
 v_power_flow.resize( number_lines );
 for( auto & var : v_power_flow )
  var.set_type( ColVariable::kContinuous );
 add_static_variable( v_power_flow , "p_flow_network" );

 if( ! f_NetworkData->get_network_cost().empty() ) {
  // the auxiliary Variable
  v_auxiliary_variable.resize( number_lines );
  for( auto & var : v_auxiliary_variable )
   var.set_type( ColVariable::kContinuous );
  add_static_variable( v_auxiliary_variable , "aux_network" );
  }
 }  // end( DCNetworkBlock::generate_PTDF_variables )

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::generate_CYCLE_variables( void )
{
 /** Implementation of "Linear Optimal Power Flow Using Cycle Flows" of
  *    Jonas Horsch, Henrik Ronellenfitsch, Dirk Witthaut, Tom Brown
  *
  * Here, we opt for the "CYCLE + FLOW" formulation with
  *   - variables "v_power_flow" as in the PTDF formulation (f_l in the paper)
  *   - variables "v_cycle_flow" (h_c in the paper) */

 const auto number_nodes = get_number_nodes();
 if( number_nodes <= 1 )
  return;

 const auto number_lines = get_number_lines();
 if( number_lines <= 0 )
  return;
  
 // the power flow variable on cycle basis
 // we know the number of cycles by the graph theory, see the paper.
 // So, no reason to call get_lines_in_cycle()
 v_cycle_flow.resize( f_NetworkData->get_lines_in_cycles().size() );

 for( auto & var : v_cycle_flow )
  var.set_type( ColVariable::kContinuous );
 add_static_variable( v_cycle_flow , "cycle_flow_network" );

 }  // end( DCNetworkBlock::generate_CYCLE_variables )

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::generate_KIRCHHOFF_variables( void )
{
 /** The Kirchhoff formulation uses:
  *  - power flow variables F_l for each line (same as PTDF / CYCLE)
  *  - voltage angle variables theta_n for each node (new)
  *
  * For pure HVDC networks (all susceptances zero), no angle variables are
  * needed since the flows are fully controllable. */

 const auto number_nodes = get_number_nodes();
 if( number_nodes <= 1 )
  return;

 // voltage angle variables are needed only if there are DC lines
 if( ! f_NetworkData->is_HVDC() ) {
  v_voltage_angle.resize( number_nodes );
  for( auto & var : v_voltage_angle )
   var.set_type( ColVariable::kContinuous );
  add_static_variable( v_voltage_angle , "voltage_angle" );
  }
 }  // end( DCNetworkBlock::generate_KIRCHHOFF_variables )

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::generate_abstract_constraints( Configuration * stcc )
{
 if( constraints_generated() )  // constraints have already been generated
  return;                       // nothing to do

 if( ( ! stcc ) && f_BlockConfig )
  stcc = f_BlockConfig->f_static_constraints_Configuration;

 f_C_v_scal = 1;
 f_tikhonov_coeff = 1e-4;
 f_ptdf_round = 1e-16; // Default value doing nothing

 if( auto SCdd = dynamic_cast< SimpleConfiguration< double > * >( stcc ) )
  f_C_v_scal = SCdd->f_value;
 else
  if( auto SCdd = dynamic_cast< SimpleConfiguration< std::pair< double ,
                                                                double > >
                                                     * >( stcc ) ) {
   f_C_v_scal = SCdd->f_value.first;
   f_tikhonov_coeff = SCdd->f_value.second;
   }
  else
   if( auto SCdd = dynamic_cast< SimpleConfiguration< std::vector< double > >
                                                      * >( stcc ) ) {
    if( SCdd->f_value.size() > 0 )
      f_C_v_scal = SCdd->f_value[ 0 ];
    if( SCdd->f_value.size() > 1 )
      f_tikhonov_coeff = SCdd->f_value[ 1 ];
    if( SCdd->f_value.size() > 2 )
      f_ptdf_round = SCdd->f_value[ 2 ];
   }

 switch( ftype ) {
  case( PTDF ):      generate_PTDF_constraints( stcc ); break;
  case( CYCLE ):     generate_CYCLE_constraints( stcc ); break;
  case( KIRCHHOFF ): generate_KIRCHHOFF_constraints( stcc ); break;
  default :
   throw( std::logic_error( "DCNetworkBlock::generate_abstract_constraints: "
			                         "unknown formulation type" ) );
  }

 // generate common constraints: bounds and cost (if there)
 generate_bound_constraints();         // generate flow limits
 generate_network_cost_constraints();  // generate cost constraints

 set_constraints_generated();  // signal all done

 }  // end( DCNetworkBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::generate_CYCLE_constraints( Configuration * stcc ) {
 /** Implementation of "Linear Optimal Power Flow Using Cycle Flows" of
  *  Jonas Horsch, Henrik Ronellenfitsch, Dirk Witthaut, Tom Brown */

 const auto number_nodes = get_number_nodes();
 if( number_nodes <= 1 )
  return;

 const auto number_lines = get_number_lines();

 if( number_lines <= 0 )
  throw( std::logic_error( "DCNetworkBlock::generate_abstract_constraints: "
                           "number of lines of DCNetworkBlock is not set" ) );

 // ----- First step: compute the cycle basis and spanning tree
 // cycle incidence matrices C_{lc} in the paper
 auto basis = f_NetworkData->get_lines_in_cycles();

 auto lines_in_tree = f_NetworkData->get_lines_in_spanning_tree();
 const auto & parent = f_NetworkData->get_spanning_parent();
 // Construct the children from the tree
 std::vector< std::vector< Index > > tree_children( number_nodes );
 for( Index i = 0 ; i < number_nodes ; ++i ) {
  int p = parent[ i ];
  if( p >= 0 && p != static_cast< int >( i ) ) {
   tree_children[ p ].push_back( i );
  }
 }

 /* eq (25): for all lines l,  f_l = sum_i T_{li} p_i  +  sum_c C_{lc} h_c */
 const auto & start_line = f_NetworkData->get_start_line();
 const auto & end_line = f_NetworkData->get_end_line();

 /*--------------------------------------------------------------*/
 /* build constraints f_l = Σ_i T_{li} p_i + Σ_c C_{lc} h_c      */
 /*--------------------------------------------------------------*/
 auto & DC_lines = f_NetworkData->get_DC_lines();
 auto & HVDC_lines = f_NetworkData->get_HVDC_lines();
 int nb_dc_lines = DC_lines.size();

 v_CYCLE_def_flow_const.resize( nb_dc_lines );

 int id_dc_line = 0;
 for( auto & line_id : DC_lines ) {
  auto lfunc = new LinearFunction();

  /* -f_l term */
  double constant_term = 0.;
  lfunc->add_variable( &v_power_flow[ line_id ] , -1.0 );

  /* Σ_i T_{li} p_i : only if l is a tree edge */
  if( lines_in_tree.contains( line_id ) ) {
   int sign = -lines_in_tree[ line_id ];
   // be careful, path FROM node TO root, i.e, in the reverse contrary to the spanning tree
   // endpoints of the DC tree edge
   Index u = start_line[ line_id ];
   Index v = end_line[ line_id ];

   const Index number_nodes = get_number_nodes();
   std::vector< bool > in_S( number_nodes , false );
   // identify the child subtree S_l
   Index child;
   if( parent[ v ] == static_cast< int >( u ) )
    child = v;
   else if( parent[ u ] == static_cast< int >( v ) )
    child = u;
   else {
#ifndef NDEBUG
    assert(
     false && "Tree edge orientation inconsistent with spanning parent" );
#endif
   }

   std::queue< Index > q;
   q.push( child );
   in_S[ child ] = true;

   while( ! q.empty() ) {
    Index x = q.front();
    q.pop();
    for( Index c : tree_children[ x ] ) {
     in_S[ c ] = true;
     q.push( c );
    }
   }

   for( Index i = 0 ; i < number_nodes ; ++i ) {
    if( ! in_S[ i ] )
     continue;

    /* ---- expand p_i ---- */
    // nodal injection variable
    lfunc->add_variable( &v_node_injection[ 0 ][ i ] , sign );
    constant_term += sign * v_ActiveDemand[ i ];
   }

   /* HVDC contributions: only if the line crosses the cut */
   for( auto hvdc_line : HVDC_lines ) {
    Index a = start_line[ hvdc_line ];
    Index b = end_line[ hvdc_line ];

    bool a_in = in_S[ a ];
    bool b_in = in_S[ b ];

    if( a_in && ! b_in )
     lfunc->add_variable( &v_power_flow[ hvdc_line ] , -1.0 * sign );
    else if( b_in && ! a_in )
     lfunc->add_variable( &v_power_flow[ hvdc_line ] , 1.0 * sign );
    // else: does not cross cut → zero contribution
   }
  }

  /* Σ_c C_{lc} h_c term */
  for( int cycle_id = 0 ; cycle_id < ( int )basis.size() ; ++cycle_id ) {
   const auto & cycle = basis[ cycle_id ];
   auto it = cycle.find( line_id );
   if( it != cycle.end() ) {
    lfunc->add_variable( &v_cycle_flow[ cycle_id ] , it->second );
   }
  }

  v_CYCLE_def_flow_const[ id_dc_line ].set_both( constant_term );
  v_CYCLE_def_flow_const[ id_dc_line ].set_function( lfunc );
  ++id_dc_line;
 }
 if( nb_dc_lines > 0 )
  add_static_constraint( v_CYCLE_def_flow_const , "v_CYCLE_def_flow_const" );

 // eq (25): forall cycle c, sum_l C_{lc}x_lf_l = 0
 if( basis.size() > 0 ) {
  v_CYCLE_def_cycle_const.resize( basis.size() );
  for( size_t cycle_id = 0 ; cycle_id < basis.size() ; ++cycle_id ) {
   auto lfunc = new LinearFunction();
   const auto & cycle = basis[ cycle_id ];
   for( const auto & [ line_id , coeff ] : cycle ) {
    lfunc->add_variable( &v_power_flow[ line_id ] ,
                         coeff / get_line_susceptance( line_id ) );
   }

   v_CYCLE_def_cycle_const[ cycle_id ].set_both( 0.0 );
   v_CYCLE_def_cycle_const[ cycle_id ].set_function( lfunc );
  }
  add_static_constraint( v_CYCLE_def_cycle_const , "v_CYCLE_def_cycle_const" );
 }

 // eq (25): sum_i p_i = 0
 auto lfunc = new LinearFunction();
 double constant_term = 0.;
 for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {
  lfunc->add_variable( &v_node_injection[ 0 ][ node_id ] , 1.0 );
  constant_term += v_ActiveDemand[ node_id ];
 }
 // In case the hypergraph is specified (only HVDC lines) and if these have non-1
 // efficiency, these need to enter the overall balance since they can imply
 // a) Losses (when flows are in the sense of the line and efficiency < 1) or
 //         against the sense of the line and efficiency > 1
 // b) Additional Generation (when flows are against the sense of the line and efficiency < 1)
 //         with the flow of the line and efficiency > 1 
 // 
 double eta;
 for( auto & hvdc_l : HVDC_lines ) {
  if( ! f_NetworkData->is_hypergraph() ) {
   // if not hyperarc, losses are 1 - eta
   eta = get_line_efficiency( hvdc_l );
   lfunc->add_variable( &v_power_flow[ hvdc_l ] , eta - 1.0 );
  }
  else { // if hyperarc, losses are 1 - sum(etas)
   auto & etas = get_line_efficiencies( hvdc_l );
   double eta_sum = std::accumulate( etas.begin() , etas.end() , 0.0 );
   lfunc->add_variable( &v_power_flow[ hvdc_l ] , eta_sum - 1.0 );
  }
 }
 overall_balanced_const.set_function( lfunc );
 overall_balanced_const.set_lhs( constant_term );
 overall_balanced_const.set_rhs( constant_term );
 add_static_constraint( overall_balanced_const , "overall_balanced_const" );

 // Add the HVDC constraints for the cycle formulation

 // If desired an additional boolean can be intercepted from the Block config and plugged here
 if( f_NetworkData->is_HVDC() || f_NetworkData->is_DC_HVDC() )
  generate_HVDC_nodal_constraints();

 /* --- old stuff
 auto & HVDC_lines = f_NetworkData->get_HVDC_lines();
 int nb_hvdc_lines = HVDC_lines.size();

 v_CYCLE_def_HVDC_const.resize( 2 * nb_hvdc_lines );
 int i_hvdc_line = 0;
 for( auto & hvdc_l : HVDC_lines ) {
   auto lfunc_1 = new LinearFunction();
   lfunc_1->add_variable( &v_node_injection[ 0 ][ start_line[ hvdc_l ] ], -1.0 );
   lfunc_1->add_variable( &v_power_flow[ hvdc_l ], 1.0 );
   v_CYCLE_def_HVDC_const[ i_hvdc_line ].set_function( lfunc_1 );
   v_CYCLE_def_HVDC_const[ i_hvdc_line ].set_both( -v_ActiveDemand[ start_line[ hvdc_l ] ] );

   auto lfunc_2 = new LinearFunction();
   lfunc_2->add_variable( &v_node_injection[ 0 ][ end_line[ hvdc_l ] ], -1.0 );
   lfunc_2->add_variable( &v_power_flow[ hvdc_l ], -1.0 );
   v_CYCLE_def_HVDC_const[ nb_hvdc_lines + i_hvdc_line ].set_function( lfunc_2 );
   v_CYCLE_def_HVDC_const[ nb_hvdc_lines + i_hvdc_line ].set_both( -v_ActiveDemand[ end_line[ hvdc_l ] ] );

   ++i_hvdc_line;
 }
 add_static_constraint( v_CYCLE_def_HVDC_const , "v_CYCLE_def_HVDC_const" );
 */
} // end( DCNetworkBlock::generate_CYCLE_constraints )

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::generate_KIRCHHOFF_constraints( Configuration * stcc )
{
 /** The Kirchhoff formulation directly encodes:
  *
  *  1. Kirchhoff's Voltage Law (KVL) for DC lines:
  *     F_l - B_l * ( theta_{from(l)} - theta_{to(l)} ) = 0
  *
  *  2. Kirchhoff's Current Law (KCL) at every node:
  *     -S_n + sum_{l: start(l)=n} F_l - sum_{l: end(l)=n} eta_l F_l = -D_n
  *
  *  3. Reference node angle fixed to zero:
  *     theta_{ref} = 0
  *
  * For pure HVDC networks (all susceptances zero), no angle variables
  * exist, and only the KCL node balance and flow limits are generated. */

 const auto number_nodes = get_number_nodes();
 if( number_nodes <= 1 )
  return;

 const auto number_lines = get_number_lines();

 if( number_lines <= 0 )
  throw( std::logic_error( "DCNetworkBlock::generate_KIRCHHOFF_constraints: "
                           "number of lines of DCNetworkBlock is not set" ) );

 LinearFunction::v_coeff_pair vars;

 const auto & start_line = f_NetworkData->get_start_line();
 const auto & end_line = f_NetworkData->get_end_line();

 // --- KVL: flow-angle definition for DC lines ---
 // F_l - B_l * theta_{from(l)} + B_l * theta_{to(l)} = 0
 if( ! f_NetworkData->is_HVDC() ) {
  auto & DC_lines = f_NetworkData->get_DC_lines();
  v_KIRCHHOFF_power_flow_def.resize( DC_lines.size() );

  Index idx = 0;
  for( auto & line_id : DC_lines ) {
   double B_l = get_line_susceptance( line_id );
   vars.emplace_back( & v_power_flow[ line_id ] , 1.0 );
   vars.emplace_back( & v_voltage_angle[ start_line[ line_id ] ] , -B_l );
   vars.emplace_back( & v_voltage_angle[ end_line[ line_id ] ] , B_l );

   v_KIRCHHOFF_power_flow_def[ idx ].set_both( 0.0 );
   v_KIRCHHOFF_power_flow_def[ idx ].set_function(
                                    new LinearFunction( std::move( vars ) ) );
   ++idx;
   }

  add_static_constraint( v_KIRCHHOFF_power_flow_def ,
                         "KIRCHHOFF_power_flow_def" );
  }

 // --- Reference angle constraint ---
 generate_reference_angle_constraint();

 // --- KCL: node power balance at every node ---
 generate_node_balance_constraints();

 }  // end( DCNetworkBlock::generate_KIRCHHOFF_constraints )

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::generate_HVDC_nodal_constraints( bool full_formulation ) {
 auto number_nodes = get_number_nodes();
 auto & DC_lines = f_NetworkData->get_DC_lines();
 auto & HVDC_lines = f_NetworkData->get_HVDC_lines();

 const auto & start_line = f_NetworkData->get_start_line();
 const auto & end_line = f_NetworkData->get_end_line();

 // power flow node injection constraints for mixed DC - HVDC- - - - - - - -
 // if we have mixed lines, we have as many as nodes impacted and touched
 // by DC lines
 int nb_DCnodes = 0;
 // savagely setting all visited nodes to true will generate nodal
 // balances for all nodes
 // the default and subtle initialization should be with false
 std::vector< bool > nodes_vist( number_nodes , full_formulation );

 // Flip any visited nodes to true
 for( auto & line_id : HVDC_lines ) {
  nodes_vist[ start_line[ line_id ] ] = true;
  if( ! f_NetworkData->is_hypergraph() )
   nodes_vist[ end_line[ line_id ] ] = true;
  else
   for( auto e : f_NetworkData->get_end_lines()[ line_id ] )
    nodes_vist[ e ] = true;
 }

 for( Index n = 0 ; n < number_nodes ; ++n ) {
  if( nodes_vist[ n ] )
   ++nb_DCnodes;
 }

 v_DC_HVDC_power_flow_const.resize( nb_DCnodes );

 // Add power balance equations for impacted nodes
 int iDCnode = 0;
 for( Index n = 0 ; n < number_nodes ; ++n ) {
  if( nodes_vist[ n ] ) {
   auto lfunc = new LinearFunction();

   double nodal_slack = 0.0; // 0.05; -- Why is there a nodal slack ?
   lfunc->add_variable( &v_node_injection[ 0 ][ n ] , -1.0 );

   double eta = 1.0;
   for( auto & line_id : HVDC_lines ) {
    if( start_line[ line_id ] == n )
     lfunc->add_variable( &v_power_flow[ line_id ] , 1.0 );
    if( ! f_NetworkData->is_hypergraph() ) {
     eta = get_line_efficiency( line_id );
     if( end_line[ line_id ] == n )
      lfunc->add_variable( &v_power_flow[ line_id ] , -eta );
    }
    else {
     for( Index i = 0 ;
          i < f_NetworkData->get_end_lines()[ line_id ].size() ; ++i ) {
      if( f_NetworkData->get_end_lines()[ line_id ][ i ] == n ) {
       eta = get_line_efficiencies( line_id )[ i ];
       lfunc->add_variable( &v_power_flow[ line_id ] , -eta );
      }
     }
    }
   }
   for( auto & line_id : DC_lines ) {
    if( start_line[ line_id ] == n )
     lfunc->add_variable( &v_power_flow[ line_id ] , 1.0 );
    if( end_line[ line_id ] == n )
     lfunc->add_variable( &v_power_flow[ line_id ] , -1.0 );
   }

   v_DC_HVDC_power_flow_const[ iDCnode ].set_lhs(
    -v_ActiveDemand[ n ] - nodal_slack );
   // allow for a 0.01 MW deviation
   v_DC_HVDC_power_flow_const[ iDCnode ].set_rhs(
    -v_ActiveDemand[ n ] + nodal_slack );
   v_DC_HVDC_power_flow_const[ iDCnode ].set_function( lfunc );

   ++iDCnode; // update the index
  }
 }
 // Add the whole vector of constraints at once
 add_static_constraint( v_DC_HVDC_power_flow_const ,
                        "DCHVDC_power_flow_injection" );
}

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::generate_network_cost_constraints( void )
{
 // 0 <= V_l - F_l  and  0 <= V_l + F_l  (linearization of |F_l|)
 if( f_NetworkData->get_network_cost().empty() )
  return;

 const auto number_lines = get_number_lines();
 LinearFunction::v_coeff_pair vars;

 v_power_flow_relax_abs.resize( MAFRC_ext()[ 2 ][ number_lines ] );

 for( Index line_id = 0 ; line_id < number_lines ; ++line_id ) {
  vars.push_back( std::make_pair( &v_power_flow[ line_id ] , -1.0 ) );
  vars.push_back( std::make_pair( &v_auxiliary_variable[ line_id ] , 1.0 ) );
  v_power_flow_relax_abs[ 0 ][ line_id ].set_lhs( 0.0 );
  v_power_flow_relax_abs[ 0 ][ line_id ].set_rhs( Inf< double >() );
  v_power_flow_relax_abs[ 0 ][ line_id ].set_function(
                                   new LinearFunction( std::move( vars ) ) );

  vars.push_back( std::make_pair( &v_power_flow[ line_id ] , 1.0 ) );
  vars.push_back( std::make_pair( &v_auxiliary_variable[ line_id ] , 1.0 ) );
  v_power_flow_relax_abs[ 1 ][ line_id ].set_lhs( 0.0 );
  v_power_flow_relax_abs[ 1 ][ line_id ].set_rhs( Inf< double >() );
  v_power_flow_relax_abs[ 1 ][ line_id ].set_function(
                                   new LinearFunction( std::move( vars ) ) );
  }

 add_static_constraint( v_power_flow_relax_abs , "power_flow_relax_abs" );

 }  // end( DCNetworkBlock::generate_network_cost_constraints )

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::generate_reference_angle_constraint( void )
{
 if( f_NetworkData->is_HVDC() )
  return;  // no angle variables for pure HVDC

 const Index ref = f_NetworkData->get_reference_node();
 v_reference_angle_const.set_lhs( 0.0 );
 v_reference_angle_const.set_rhs( 0.0 );
 v_reference_angle_const.set_variable( & v_voltage_angle[ ref ] );
 add_static_constraint( v_reference_angle_const , "reference_angle" );

 }  // end( DCNetworkBlock::generate_reference_angle_constraint )

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::generate_node_balance_constraints( void )
{
 // KCL at every node:
 // -S_n + sum_{l: start(l)=n} F_l - sum_{l: end(l)=n} eta_l F_l = -D_n

 const auto number_nodes = get_number_nodes();
 const auto & start_line = f_NetworkData->get_start_line();
 const auto & end_line = f_NetworkData->get_end_line();
 LinearFunction::v_coeff_pair vars;

 v_KIRCHHOFF_node_balance_const.resize( number_nodes );

 for( Index n = 0 ; n < number_nodes ; ++n ) {
  vars.emplace_back( & v_node_injection[ 0 ][ n ] , -1.0 );

  // DC lines (regular arcs, never hyperarcs)
  if( ! f_NetworkData->is_HVDC() ) {
   for( auto & l : f_NetworkData->get_DC_lines() ) {
    if( start_line[ l ] == n )
     vars.emplace_back( & v_power_flow[ l ] , 1.0 );
    double eta = get_line_efficiency( l );
    if( end_line[ l ] == n )
     vars.emplace_back( & v_power_flow[ l ] , -eta );
    }
   }

  // HVDC lines (maybe hyperarcs)
  if( ! f_NetworkData->is_DC() ) {
   for( auto & l : f_NetworkData->get_HVDC_lines() ) {
    if( start_line[ l ] == n )
     vars.emplace_back( & v_power_flow[ l ] , 1.0 );

    if( ! f_NetworkData->is_hypergraph() ) {
     double eta = get_line_efficiency( l );
     if( end_line[ l ] == n )
      vars.emplace_back( & v_power_flow[ l ] , -eta );
     }
    else {
     for( Index i = 0 ;
          i < f_NetworkData->get_end_lines()[ l ].size() ; ++i ) {
      if( f_NetworkData->get_end_lines()[ l ][ i ] == n ) {
       double eta = get_line_efficiencies( l )[ i ];
       vars.emplace_back( & v_power_flow[ l ] , -eta );
       }
      }
     }
    }
   }

  v_KIRCHHOFF_node_balance_const[ n ].set_both( -v_ActiveDemand[ n ] );
  v_KIRCHHOFF_node_balance_const[ n ].set_function(
                                    new LinearFunction( std::move( vars ) ) );
  }

 add_static_constraint( v_KIRCHHOFF_node_balance_const ,
                        "KIRCHHOFF_node_balance" );

 }  // end( DCNetworkBlock::generate_node_balance_constraints )

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::generate_PTDF_constraints( Configuration * stcc )
{
 /* This formulation corresponds to the "PTDF + FLOW" formulation of
  * "Linear Optimal Power Flow Using Cycle Flows" of
  * Jonas Horsch, Henrik Ronellenfitsch, Dirk Witthaut, Tom Brown */

 auto number_nodes = get_number_nodes();
 if( number_nodes <= 1 )
  return;

 auto number_lines = get_number_lines();

 if( number_lines <= 0 )
  throw( std::logic_error( "DCNetworkBlock::generate_abstract_constraints: "
                           "number of lines of DCNetworkBlock is not set" ) );

 LinearFunction::v_coeff_pair vars;

 const auto & start_line = f_NetworkData->get_start_line();
 const auto & end_line = f_NetworkData->get_end_line();

 // Splitting DC and HVDC part- - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 auto & DC_lines = f_NetworkData->get_DC_lines();
 auto & HVDC_lines = f_NetworkData->get_HVDC_lines();

 // constraints on the HVDC part- - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! f_NetworkData->is_DC() ) {  // ... if any

  // power flow node injection constraints for pure HVDC- - - - - - - - - - -
  if( f_NetworkData->is_HVDC() ) {
   v_power_flow_injection_const.resize( number_nodes );

   for( Index n = 0 ; n < number_nodes ; ++n ) {
    vars.push_back( std::make_pair( &v_node_injection[ 0 ][ n ] , -1.0 ) );
    double eta = 1.0;
    for( auto & line_id : HVDC_lines ) {
     // start node
     if( start_line[ line_id ] == n )
      vars.push_back( std::make_pair( &v_power_flow[ line_id ] , 1.0 ) );

     if( ! f_NetworkData->is_hypergraph() ) {
      // if no hyperarch -> normal behavior from get_line_efficiency
      eta = get_line_efficiency( line_id );
      // efficiency of the HVDC line
      if( end_line[ line_id ] == n )
       vars.push_back( std::make_pair( &v_power_flow[ line_id ] , -eta ) );
      }
     else { // if hyperarch -> loop over v_end_lines and get_line_efficiencies
      auto & end_lines = f_NetworkData->get_end_lines()[ line_id ];
      for( Index i = 0 ; i < end_lines.size() ; ++i ) {
       if( end_lines[ i ] == n ) {
        eta = get_line_efficiencies( line_id )[ i ];
        vars.push_back( std::make_pair( &v_power_flow[ line_id ] , -eta ) );
        }
       }
      }
     }

    v_power_flow_injection_const[ n ].set_both( -v_ActiveDemand[ n ] );
    v_power_flow_injection_const[ n ].set_function(
				   new LinearFunction( std::move( vars ) ) );
    }

   add_static_constraint( v_power_flow_injection_const ,
			  "HVDC_power_flow_injection" );
   }

   // If desired an additional boolean can be intercepted from the Block config and plugged here
   if( f_NetworkData->is_DC_HVDC() ) {   
    generate_HVDC_nodal_constraints( );
   }
 }

 // Constraints on the DC part- - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! f_NetworkData->is_HVDC() ) {  // ... if any

  // compute PTDF and, if necessary, DCDF
  SpMat PTDF_matrix = f_NetworkData->get_PTDF( DC_lines , f_tikhonov_coeff);
  if( f_NetworkData->is_DC_HVDC() )
   if( ! f_NetworkData->was_DCDF_computed() )
    f_NetworkData->compute_DCDF( HVDC_lines , PTDF_matrix );

  // Flow limit constraints
  v_power_flow_def.resize( number_lines );

  for( auto & line_id : DC_lines ) {
   // TODO : verify if this does not entail a copy of the information
   // which would be inefficient
   Eigen::SparseMatrix< double > a_row =
                   PTDF_matrix.block( line_id , 0 , 1 , PTDF_matrix.cols() );

   vars.push_back( std::make_pair( &v_power_flow[ line_id ] , -1.0 ) );
   double constant_term = 0;

   // loop over all nodes
   // for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {
   for( int k = 0 ; k < a_row.outerSize() ; ++k ) {
    for( Eigen::SparseMatrix< double >::InnerIterator it( a_row , k ) ; it ;
	 ++it ) {
     int node_id = f_NetworkData->get_originalIdx( it.col() );
     if( node_id != f_NetworkData->get_reference_node() ) {
      double coefficient = round_to( it.value() , f_ptdf_round );
      // Distribution Factor Matrix
      vars.push_back( std::make_pair( &v_node_injection[ 0 ][ node_id ] ,
                                      coefficient ) );
      constant_term += coefficient * v_ActiveDemand[ node_id ];
      }
     }
    }  // for each node

   if( f_NetworkData->is_DC_HVDC() ) {
    SpMat DCDF_ = f_NetworkData->get_DCDF();
    // TODO : Also only loop over the non zero entries of DCDF only ...
    for( auto & dc_line_id : HVDC_lines ) {
     double coeff = round_to( DCDF_.coeff( line_id , dc_line_id ) ,
                              f_ptdf_round );
     vars.push_back( std::make_pair( &v_power_flow[ dc_line_id ] , coeff ) );
     }
    }

   // Set the constraint (AC)
   v_power_flow_def[ line_id ].set_function(
                                    new LinearFunction( std::move( vars ) ) );
   v_power_flow_def[ line_id ].set_lhs( constant_term );
   v_power_flow_def[ line_id ].set_rhs( constant_term );

   }  // end( for( DC lines ) )

  add_static_constraint( v_power_flow_def , "AC/HVDC_powerflow_def" );

  double constant_term = 0.;
  for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {
   vars.push_back( std::make_pair( &v_node_injection[ 0 ][ node_id ] , 1. ) );
   constant_term += v_ActiveDemand[ node_id ];
   }

  overall_balanced_const.set_function( new LinearFunction(
						       std::move( vars ) ) );
  overall_balanced_const.set_lhs( constant_term );
  overall_balanced_const.set_rhs( constant_term );
  add_static_constraint( overall_balanced_const , "overall_balanced_const" );
  }  // end( if( there are DC lines ) )
 }  // end( DCNetworkBlock::generate_PTDF_constraints )

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::generate_bound_constraints( void )
{
 /*-----------------------------------------------------------------------*/
 /*-------------------- flow limits with/without design ------------------*/
 /*-----------------------------------------------------------------------*/

 LinearFunction::v_coeff_pair vars;

 const auto number_lines = get_number_lines();

 /*-------------------------- with design --------------------------------*/
 /** For lines having a design variable x_l (get_design( l ) != nullptr),
  *  impose:
  *
  *   LOWER: F_l - kappa * MinP_l * x_l >= 0
  *   UPPER: F_l - kappa * MaxP_l * x_l <= 0 */

 if( has_design() ) {
  v_power_flow_limit_design_const.resize(  MAFRC_ext()[ 2 ][ number_lines ] );

  for( Index l = 0 ; l < number_lines ; ++l ) {
   ColVariable * x = get_design( l );
   if( ! x )   // no design on this line:
    continue;  // handled in the "without design" block

   const double kappa = get_kappa( l );
   const double Pmn   = get_min_power_flow( l );
   const double Pmx   = get_max_power_flow( l );

   // LOWER:  F_l - kappa * Pmn * x_l >= 0
   vars.emplace_back( & v_power_flow[ l ] , 1.0 );
   vars.emplace_back( x , - kappa * Pmn );
   v_power_flow_limit_design_const[ 0 ][ l ].set_lhs( 0.0 );
   v_power_flow_limit_design_const[ 0 ][ l ].set_rhs( Inf< double >() );
   v_power_flow_limit_design_const[ 0 ][ l ].set_function(
                                   new LinearFunction( std::move( vars ) ) );

   // UPPER:  F_l - kappa * Pmx * x_l <= 0
   vars.emplace_back( & v_power_flow[ l ] , 1.0 );
   vars.emplace_back( x , - kappa * Pmx );
   v_power_flow_limit_design_const[ 1 ][ l ].set_lhs( -Inf< double >() );
   v_power_flow_limit_design_const[ 1 ][ l ].set_rhs( 0.0 );
   v_power_flow_limit_design_const[ 1 ][ l ].set_function(
                                   new LinearFunction( std::move( vars ) ) );
   }

  add_static_constraint( v_power_flow_limit_design_const ,
			 "Power_flow_limit_design" );
  }

 /*------------------------- without design -------------------------------*/
 /** For lines with no design variable (get_design( l ) == nullptr),
  *  impose the standard box:
  *
  *      kappa * MinP_l  <=  F_l  <=  kappa * MaxP_l */

 if( ! all_design() ) {
  v_power_flow_limit_const.resize( number_lines );

  for( Index l = 0 ; l < number_lines ; ++l ) {
   if( get_design( l ) )  // already handled above
    continue;

   const double kappa = get_kappa( l );
   v_power_flow_limit_const[ l ].set_lhs( kappa * f_C_v_scal *
					  get_min_power_flow( l ) );
   v_power_flow_limit_const[ l ].set_rhs( kappa * f_C_v_scal *
					  get_max_power_flow( l ) );
   v_power_flow_limit_const[ l ].set_variable( &v_power_flow[ l ] );
   }

  add_static_constraint( v_power_flow_limit_const , "Power_flow_limit" );
  }
 }  // end( DCNetworkBlock::generate_bound_constraints )

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::generate_objective( Configuration * objc )
{
 if( objective_generated() )  // Objective has already been generated
  return;                     // nothing to do

 LinearFunction::v_coeff_pair vars;

 if( const auto & nc = f_NetworkData->get_network_cost() ; ! nc.empty() ) {
  auto nl = get_number_lines();
  vars.resize( nl );
 
  for( Index line_id = 0 ; line_id < nl ; ++line_id )
   vars[ line_id ] = std::make_pair( & v_auxiliary_variable[ line_id ] ,
				     nc[ line_id ] );
  }

 auto lf = new LinearFunction( std::move( vars ) );
 lf->set_constant_term( f_ConstTerm );

 objective.set_function( lf );
 objective.set_sense( Objective::eMin );

 this->set_objective( &objective , eNoMod );  // set Block objective

 set_objective_generated();

 }  // end( DCNetworkBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/

Solution * DCNetworkBlock::get_Solution( Configuration * csolc ,
					 bool emptys )
{
 Index wsol = 7;
 if( ( ! csolc ) && f_BlockConfig )
  csolc = f_BlockConfig->f_solution_Configuration;

 if( auto config = dynamic_cast< SimpleConfiguration< int > * >( csolc ) )
  wsol = config->f_value;

 // call the method of the base class
 auto * sol = dynamic_cast< DCNetworkBlockSolution * >(
			     NetworkBlock::get_Solution( csolc , emptys ) );
 assert( sol );

 if( wsol & 2 )
  sol->v_flow.resize( get_number_lines() );

 if( wsol & 4 )
  sol->v_cost.resize( get_number_lines() );

 if( ! emptys )
  sol->read( this );

 return( sol );

 }  // end( DCNetworkBlock::get_Solution )

/*--------------------------------------------------------------------------*/

NetworkBlockSolution * DCNetworkBlock::new_Solution( void ) const {
 return( new DCNetworkBlockSolution() );
 }

/*--------------------------------------------------------------------------*/
/*----------------- METHODS FOR CHECKING THE DCNetworkBlock ----------------*/
/*--------------------------------------------------------------------------*/

bool DCNetworkBlock::is_feasible( bool useabstract , Configuration * fsbc )
{
 // Retrieve the tolerance and the type of violation.
 double tol = 0;
 bool rel_viol = true;

 // Try to extract, from "c", the parameters that determine feasibility.
 // If it succeeds, it sets the values of the parameters and returns
 // true. Otherwise, it returns false.
 auto extract_parameters = [ & tol , & rel_viol ]( Configuration * c )
  -> bool {
  if( auto tc = dynamic_cast< SimpleConfiguration< double > * >( c ) ) {
   tol = tc->f_value;
   return( true );
   }
  if( auto tc = dynamic_cast< SimpleConfiguration<
                                     std::pair< double , int > > * >( c ) ) {
   tol = tc->f_value.first;
   rel_viol = tc->f_value.second;
   return( true );
   }
  return( false );
  };

 if( ( ! extract_parameters( fsbc ) ) && f_BlockConfig )
  // if the given Configuration is not valid, try the one from the BlockConfig
  extract_parameters( f_BlockConfig->f_is_feasible_Configuration );

 return(
  NetworkBlock::is_feasible( useabstract )
  // Variables
  && ColVariable::is_feasible( v_node_injection , tol )
  && ColVariable::is_feasible( v_power_flow , tol )
  && ColVariable::is_feasible( v_auxiliary_variable , tol )
  && ColVariable::is_feasible( v_cycle_flow , tol )
  && ColVariable::is_feasible( v_voltage_angle , tol )
  // Constraints
  && RowConstraint::is_feasible( v_power_flow_limit_const , tol , rel_viol )
  && RowConstraint::is_feasible( v_power_flow_limit_design_const , tol ,
         rel_viol )
  && RowConstraint::is_feasible( v_power_flow_injection_const , tol ,
         rel_viol )
  && RowConstraint::is_feasible( v_power_flow_def , tol , rel_viol )
  && RowConstraint::is_feasible( v_power_flow_relax_abs , tol , rel_viol )
  && RowConstraint::is_feasible( overall_balanced_const , tol , rel_viol )
  && RowConstraint::is_feasible( node_injection_bounds_const , tol ,
         rel_viol )
  && RowConstraint::is_feasible( v_CYCLE_def_flow_const , tol , rel_viol )
  && RowConstraint::is_feasible( v_CYCLE_def_cycle_const , tol , rel_viol )
  // Kirchhoff formulation
  && RowConstraint::is_feasible( v_KIRCHHOFF_power_flow_def , tol , rel_viol )
  && RowConstraint::is_feasible( v_KIRCHHOFF_node_balance_const , tol ,
         rel_viol )
  && RowConstraint::is_feasible( v_reference_angle_const , tol , rel_viol )
  );

 }  // end( DCNetworkBlock::is_feasible )

/*--------------------------------------------------------------------------*/
/*--------- METHODS FOR LOADING, PRINTING & SAVING THE DCNetworkBlock ------*/
/*--------------------------------------------------------------------------*/

void DCNetworkData::serialize( netCDF::NcGroup & group ) const
{
 NetworkData::serialize( group );

 if( f_number_nodes == 1 )
  return;

 auto NumberLines = group.addDim( "NumberLines" , f_number_lines );
 if( f_reference_node )
  group.addDim( "ReferenceNode" , f_reference_node );

 if( is_hypergraph() ) {  // an hypergraph
  auto NumberBranches = group.addDim( "NumberBranches" , f_number_branches );

  // remap v_start_line, v_end_lines and v_h_efficiency into vectors
  // f_number_branches-long, build the vector of branches id
  std::vector< Index > id( f_number_branches );
  std::vector< Index > sn( f_number_branches );
  std::vector< Index > en( f_number_branches );
  std::vector< double > eff( f_number_branches );

  Index curr = 0;
  for( Index i = 0 ; i < f_number_lines ; ++i )
   for( Index j = 0 ; j < v_end_lines[ i ].size() ; ++j , ++curr ) {
    id[ curr ] = i;
    sn[ curr ] = v_start_line[ i ];
    en[ curr ] = v_end_lines[ i ][ j ];
    eff[ curr ] = v_h_efficiency[ i ][ j ];
    }

  ::serialize( group , "StartLine" , netCDF::NcUint() , NumberBranches , sn );

  ::serialize( group , "EndLine" , netCDF::NcUint() , NumberBranches , en );

  ::serialize( group , "Efficiency" , netCDF::NcDouble() , NumberBranches ,
         eff );

  ::serialize( group , "HyperArcID" , netCDF::NcUint() , NumberBranches ,
         id );
  }
 else {  // a regular graph
  ::serialize( group , "StartLine" , netCDF::NcUint() , NumberLines ,
               v_start_line );

  ::serialize( group , "EndLine" , netCDF::NcUint() , NumberLines ,
               v_end_line );

  ::serialize( group , "Efficiency" , netCDF::NcDouble() , NumberLines ,
              v_efficiency );
  }

 ::serialize( group , "MaxPowerFlow" , netCDF::NcDouble() , NumberLines ,
              v_max_power_flow );

 ::serialize( group , "MinPowerFlow" , netCDF::NcDouble() , NumberLines ,
              v_min_power_flow );

 if( ! v_line_susceptance.empty() )
  ::serialize( group , "LineSusceptance" , netCDF::NcDouble() , NumberLines ,
	              v_line_susceptance );

 ::serialize( group , "NetworkCost" , netCDF::NcDouble() , NumberLines ,
              v_network_cost );

 if( ! v_line_names.empty() ) {
  auto LineName = group.addVar( "LineName" , netCDF::NcString() ,
        NumberLines );
  for( Index i = 0 ; i < v_line_names.size() ; ++i )
   LineName.putVar( { i } , v_line_names[ i ] );
  }
 }  // end( DCNetworkData::serialize )

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::serialize( netCDF::NcGroup & group ) const
{
 NetworkBlock::serialize( group );

 if( auto network_data = get_NetworkData() )
  // If a DCNetworkData is present, serialize it
  network_data->serialize( group );

 auto NumberNodes = group.getDim( "NumberNodes" );

 if( ! v_ActiveDemand.empty() ) {
  // This DCNetworkBlock has active demand, so it is serialized.

  if( NumberNodes.isNull() )
   /* The dimension "NumberNodes" is not present in the group (which means
    * that a DCNetworkData is not present). However, the number of nodes can
    * still be obtained from the size of the active demand vector. Notice that
    * the name "NumberNodes" is not used for this new dimension, because it
    * would indicate that a DCNetworkData is present (which is not the
    * case). Therefore, we create an alternative dimension in order to be able
    * to serialize the active demand. */
   NumberNodes = group.addDim( "__NumberNodes__" , v_ActiveDemand.size() );

  // Finally, serialize the active demand.
  ::serialize( group , "ActiveDemand" , netCDF::NcDouble() ,
               NumberNodes , v_ActiveDemand );
  }
 }  // end( DCNetworkBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

void DCNetworkBlock::set_active_demand( MF_dbl_it values , Subset && subset ,
                                        bool ordered ,
                                        c_ModParam issuePMod ,
                                        c_ModParam issueAMod )
{
 if( subset.empty() )
  return;

 if( v_ActiveDemand.empty() ) {
  if( std::all_of( values , values + subset.size() ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  v_ActiveDemand.resize( get_number_nodes() );
  }

 bool identical = true;
 for( auto i : subset ) {
  if( i >= v_ActiveDemand.size() )
   throw( std::invalid_argument( "DCNetworkBlock::set_active_demand: "
                                 "invalid value in subset" ) );

  auto demand = *( values++ );
  if( v_ActiveDemand[ i ] != demand ) {
   identical = false;

   if( not_dry_run( issuePMod ) )  // Change the physical representation
    v_ActiveDemand[ i ] = demand;
   }
  }

 if( identical )
  return;  // nothing changes; return

 if( not_dry_run( issuePMod ) && not_dry_run( issueAMod ) &&
     constraints_generated() ) {
  std::vector< Index > modified_nodes( subset.begin() , subset.end() );
  if( f_NetworkData->is_HVDC() )
   change_DC_power_flow_injection_constraints( modified_nodes , issueAMod );
  else
   change_active_demand_constraints( modified_nodes , issueAMod );
  }

 if( issue_pmod( issuePMod ) ) {  // Issue a Physical Modification
  if( ! ordered )
   std::sort( subset.begin() , subset.end() );

  Block::add_Modification( std::make_shared< NetworkBlockSbstMod >(
                            this , NetworkBlockMod::eSetActD ,
                            std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );
  }
 }  // end( DCNetworkBlock::set_active_demand( subset ) )

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::set_active_demand( MF_dbl_it values , Range rng ,
                                        c_ModParam issuePMod ,
                                        c_ModParam issueAMod )
{
 rng.second = std::min( rng.second , get_number_nodes() );
 if( rng.second <= rng.first )
  return;

 if( v_ActiveDemand.empty() ) {
  if( std::all_of( values , values + ( rng.second - rng.first ) ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  v_ActiveDemand.resize( get_number_nodes() );
  }

 // If nothing changes, return
 if( std::equal( values , values + ( rng.second - rng.first ) ,
                 v_ActiveDemand.begin() + rng.first ) )
  return;

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation
  std::copy( values , values + ( rng.second - rng.first ) ,
             v_ActiveDemand.begin() + rng.first );

  if( not_dry_run( issueAMod ) && constraints_generated() ) {
   std::vector< Index > modified_nodes( rng.second - rng.first );
   std::iota( modified_nodes.begin() , modified_nodes.end() , rng.first );
   if( f_NetworkData->is_HVDC() )
    change_DC_power_flow_injection_constraints( modified_nodes , issueAMod );
   else
    change_active_demand_constraints( modified_nodes , issueAMod );
   }
  }

 if( issue_pmod( issuePMod ) ) // issue a Physical Modification
  Block::add_Modification( std::make_shared< NetworkBlockRngdMod >( this ,
                            NetworkBlockMod::eSetActD , rng ) ,
                           Observer::par2chnl( issuePMod ) );
 } // end( DCNetworkBlock::set_active_demand( range ) )

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::change_active_demand_constraints(
                           c_Subset & modified_nodes ,
                           c_ModParam issueAMod )
{
 switch( ftype ) {
  case( PTDF ):
   change_active_demand_constraints_PTDF( modified_nodes , issueAMod );
   break;
  case( CYCLE ):
   change_active_demand_constraints_CYCLE( modified_nodes , issueAMod );
   break;
  case( KIRCHHOFF ):
   change_active_demand_constraints_KIRCHHOFF( modified_nodes , issueAMod );
   break;
  default:
   throw( std::logic_error(
    "DCNetworkBlock::change_active_demand_constraints: unknown formulation type"
   ) );
  }
 }  // end( DCNetworkBlock::change_active_demand_constraints )

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::change_active_demand_constraints_PTDF(
                           c_Subset & modified_nodes ,
                           c_ModParam issueAMod )
{
 auto & DC_lines = f_NetworkData->get_DC_lines();
 SpMat PTDF_matrix = f_NetworkData->get_PTDF( DC_lines , f_tikhonov_coeff );

 const Index ref = f_NetworkData->get_reference_node();

 for( auto & line_id : DC_lines ) {
  double constant_term = 0.0;

  Eigen::SparseMatrix< double > a_row =
      PTDF_matrix.block( line_id , 0 , 1 , PTDF_matrix.cols() );

  for( int k = 0 ; k < a_row.outerSize() ; ++k ) {
   for( Eigen::SparseMatrix< double >::InnerIterator it( a_row , k ) ; it ;
        ++it ) {
    int node_id = f_NetworkData->get_originalIdx( it.col() );
    if( node_id != ref ) {
     double coefficient = round_to( it.value() , f_ptdf_round );
     constant_term += coefficient * v_ActiveDemand[ node_id ];
    }
        }
  }

  v_power_flow_def[ line_id ].set_lhs( constant_term , issueAMod );
  v_power_flow_def[ line_id ].set_rhs( constant_term , issueAMod );
 }

 double balance_const = 0.0;
 for( Index node_id = 0 ; node_id < get_number_nodes() ; ++node_id )
  balance_const += v_ActiveDemand[ node_id ];

 overall_balanced_const.set_lhs( balance_const , issueAMod );
 overall_balanced_const.set_rhs( balance_const , issueAMod );

 if( f_NetworkData->is_DC_HVDC() )
  change_DC_HVDC_power_flow_injection_constraints( modified_nodes , issueAMod );
}

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::change_active_demand_constraints_CYCLE(
                           c_Subset & modified_nodes ,
                           c_ModParam issueAMod )
{
 const auto number_nodes = get_number_nodes();
 const auto & start_line = f_NetworkData->get_start_line();
 const auto & end_line = f_NetworkData->get_end_line();

 auto basis = f_NetworkData->get_lines_in_cycles();
 auto lines_in_tree = f_NetworkData->get_lines_in_spanning_tree();
 const auto & parent = f_NetworkData->get_spanning_parent();

 std::vector< std::vector< Index > > tree_children( number_nodes );
 for( Index i = 0 ; i < number_nodes ; ++i ) {
  int p = parent[ i ];
  if( p >= 0 && p != static_cast< int >( i ) )
   tree_children[ p ].push_back( i );
 }

 auto & DC_lines = f_NetworkData->get_DC_lines();
 auto & HVDC_lines = f_NetworkData->get_HVDC_lines();

 int id_dc_line = 0;
 for( auto & line_id : DC_lines ) {
  double constant_term = 0.0;

  if( lines_in_tree.contains( line_id ) ) {
   int sign = -lines_in_tree[ line_id ];

   Index u = start_line[ line_id ];
   Index v = end_line[ line_id ];

   std::vector< bool > in_S( number_nodes , false );
   Index child;

   if( parent[ v ] == static_cast< int >( u ) )
    child = v;
   else if( parent[ u ] == static_cast< int >( v ) )
    child = u;
   else
    throw( std::logic_error(
     "DCNetworkBlock::change_active_demand_constraints_CYCLE: "
     "tree edge orientation inconsistent with spanning parent" ) );

   std::queue< Index > q;
   q.push( child );
   in_S[ child ] = true;

   while( ! q.empty() ) {
    Index x = q.front();
    q.pop();
    for( Index c : tree_children[ x ] ) {
     in_S[ c ] = true;
     q.push( c );
    }
   }

   for( Index i = 0 ; i < number_nodes ; ++i ) {
    if( in_S[ i ] )
     constant_term += sign * v_ActiveDemand[ i ];
   }
  }

  v_CYCLE_def_flow_const[ id_dc_line ].set_lhs( constant_term , issueAMod );
  v_CYCLE_def_flow_const[ id_dc_line ].set_rhs( constant_term , issueAMod );
  ++id_dc_line;
 }

 double balance_const = 0.0;
 for( Index node_id = 0 ; node_id < number_nodes ; ++node_id )
  balance_const += v_ActiveDemand[ node_id ];

 overall_balanced_const.set_lhs( balance_const , issueAMod );
 overall_balanced_const.set_rhs( balance_const , issueAMod );

 if( f_NetworkData->is_DC_HVDC() )
  change_DC_HVDC_power_flow_injection_constraints( modified_nodes , issueAMod );
}

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::change_active_demand_constraints_KIRCHHOFF(
                           c_Subset & modified_nodes ,
                           c_ModParam issueAMod )
{
 for( auto & n : modified_nodes )
  v_KIRCHHOFF_node_balance_const[ n ].set_both(
   -v_ActiveDemand[ n ] , issueAMod );

 if( f_NetworkData->is_DC_HVDC() )
  change_DC_HVDC_power_flow_injection_constraints( modified_nodes , issueAMod );
}

/*--------------------------------------------------------------------------*/


void DCNetworkBlock::change_DC_HVDC_power_flow_injection_constraints(
                           c_Subset & modified_nodes ,
                           c_ModParam issueAMod )
{
 const auto number_nodes = get_number_nodes();
 auto & HVDC_lines = f_NetworkData->get_HVDC_lines();
 const auto & start_line = f_NetworkData->get_start_line();

 std::vector< bool > nodes_vist( number_nodes , false );

 for( auto & line_id : HVDC_lines ) {
  nodes_vist[ start_line[ line_id ] ] = true;
  if( ! f_NetworkData->is_hypergraph() )
   nodes_vist[ f_NetworkData->get_end_line()[ line_id ] ] = true;
  else
   for( auto e : f_NetworkData->get_end_lines()[ line_id ] )
    nodes_vist[ e ] = true;
 }

 int iDCnode = 0;
 for( Index n = 0 ; n < number_nodes ; ++n ) {
  if( ! nodes_vist[ n ] )
   continue;

  if( std::find( modified_nodes.begin() , modified_nodes.end() , n ) !=
      modified_nodes.end() ) {
   v_DC_HVDC_power_flow_const[ iDCnode ].set_lhs( -v_ActiveDemand[ n ] ,
                                                  issueAMod );
   v_DC_HVDC_power_flow_const[ iDCnode ].set_rhs( -v_ActiveDemand[ n ] ,
                                                  issueAMod );
      }

  ++iDCnode;
 }
}

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::set_network_cost( MF_dbl_it values ,
                                       Subset && subset ,
                                       bool ordered ,
                                       c_ModParam issuePMod ,
                                       c_ModParam issueAMod )
{
 if( subset.empty() )
  return;

 auto & network_cost = f_NetworkData->get_network_cost();

 if( network_cost.empty() ) {
  if( std::all_of( values , values + subset.size() ,
                   []( double cst ) { return( cst == 0.0 ); } ) )
   return;

  network_cost.assign( get_number_lines() , 0.0 );
  }

 bool identical = true;
 auto values_it = values;
 for( auto i : subset ) {
  if( i >= network_cost.size() )
   throw( std::invalid_argument(
    "DCNetworkBlock::set_network_cost: invalid value in subset: " +
    std::to_string( i ) ) );

  if( network_cost[ i ] != *( values_it++ ) ) {
   identical = false;
   break;
   }
  }

 if( identical )
  return;

 if( not_dry_run( issuePMod ) ) {
  values_it = values;
  for( auto i : subset )
   network_cost[ i ] = *( values_it++ );

  if( not_dry_run( issueAMod ) && objective_generated() ) {
   auto * lf = static_cast< LinearFunction * >( objective.get_function() );

   for( auto i : subset ) {
    const auto idx = lf->is_active( &v_auxiliary_variable[ i ] );

    if( idx == Inf< Index >() )
     throw( std::logic_error(
      "DCNetworkBlock::set_network_cost: expected Variable not found in "
      "objective." ) );

    lf->modify_coefficient( idx ,
                            network_cost[ i ] ,
                            issueAMod );
    }
   }
  }

 if( issue_pmod( issuePMod ) ) {
  if( ! ordered )
   std::sort( subset.begin() , subset.end() );

  Block::add_Modification( std::make_shared< DCNetworkBlockSbstMod >(
                            this , DCNetworkBlockMod::eSetNetCost ,
                            std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );
  }
}  // end( DCNetworkBlock::set_network_cost( subset ) )

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::set_network_cost( MF_dbl_it values ,
                                       Range rng ,
                                       c_ModParam issuePMod ,
                                       c_ModParam issueAMod )
{
 rng.second = std::min( rng.second , get_number_lines() );
 if( rng.second <= rng.first )
  return;

 c_Index sz = rng.second - rng.first;

 auto & network_cost = f_NetworkData->get_network_cost();

 if( network_cost.empty() ) {
  if( std::all_of( values , values + sz ,
                   []( double cst ) { return( cst == 0.0 ); } ) )
   return;

  network_cost.assign( get_number_lines() , 0.0 );
 }

 if( std::equal( values ,
                 values + sz ,
                 network_cost.begin() + rng.first ) )
  return;

 if( not_dry_run( issuePMod ) ) {
  std::copy( values ,
             values + sz ,
             network_cost.begin() + rng.first );

  if( not_dry_run( issueAMod ) && objective_generated() ) {
   auto * lf = static_cast< LinearFunction * >( objective.get_function() );

   for( Index i = rng.first ; i < rng.second ; ++i ) {
    const auto idx = lf->is_active( &v_auxiliary_variable[ i ] );

    if( idx == Inf< Index >() )
     throw( std::logic_error(
      "DCNetworkBlock::set_network_cost: expected Variable not found in "
      "objective." ) );

    lf->modify_coefficient( idx ,
                            network_cost[ i ] ,
                            issueAMod );
   }
  }
 }

 if( issue_pmod( issuePMod ) )
  Block::add_Modification( std::make_shared< DCNetworkBlockRngdMod >(
                            this , DCNetworkBlockMod::eSetNetCost , rng ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( DCNetworkBlock::set_network_cost( range ) )

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::set_kappa( MF_dbl_it values , Subset && subset ,
				bool ordered ,
                                c_ModParam issuePMod ,
                                c_ModParam issueAMod )
{
 if( subset.empty() )
  return;

 if( v_kappa.empty() ) {
  if( std::all_of( values , values + subset.size() ,
                   []( double cst ) { return( cst == 1 ); } ) )
   return;

  v_kappa.resize( get_number_lines() , 1 );
  }

 bool identical = true;
 for( auto i : subset ) {
  if( i >= v_kappa.size() )
   throw( std::invalid_argument( "DCNetworkBlock::set_kappa: invalid value in"
                                 " subset: " + std::to_string( i ) ) );

  const auto kappa = *( values++ );
  if( v_kappa[ i ] != kappa ) {
   identical = false;
   if( not_dry_run( issuePMod ) )
    // Change the physical representation
    v_kappa[ i ] = kappa;
   }
  }

 if( identical )
  return;  // nothing changes; return

 if( not_dry_run( issuePMod ) && not_dry_run( issueAMod ) &&
     constraints_generated() ) {

  // Change the abstract representation
  std::vector< Index > modified_lines( subset.begin() , subset.end() );
  change_power_flow_limit_constraints( modified_lines , issueAMod );
  // kappa only appears in limit constraints
  }

 if( issue_pmod( issuePMod ) ) {  // issue a Physical Modification
  if( ! ordered )
   std::sort( subset.begin() , subset.end() );

  Block::add_Modification( std::make_shared< DCNetworkBlockSbstMod >( this ,
           DCNetworkBlockMod::eSetKappa , std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );
  }
 }  // end( DCNetworkBlock::set_kappa( subset ) )

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::set_kappa( MF_dbl_it values , Range rng ,
                                c_ModParam issuePMod , c_ModParam issueAMod )
{
 rng.second = std::min( rng.second , get_number_lines() );
 if( rng.second <= rng.first )
  return;

 if( v_kappa.empty() ) {
  if( std::all_of( values , values + ( rng.second - rng.first ) ,
                   []( double cst ) { return( cst == 1 ); } ) )
   return;

  v_kappa.resize( get_number_lines() , 1 );
  }

 // If nothing changes, return
 if( std::equal( values , values + ( rng.second - rng.first ) ,
                 v_kappa.begin() + rng.first ) )
  return;

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation
  std::copy( values , values + ( rng.second - rng.first ) ,
             v_kappa.begin() + rng.first );

  if( not_dry_run( issueAMod ) && constraints_generated() ) {
   // Change the abstract representation
   std::vector< Index > modified_lines( rng.second - rng.first );
   std::iota( modified_lines.begin() , modified_lines.end() , rng.first );
   change_power_flow_limit_constraints( modified_lines , issueAMod );
   }
  }

 if( issue_pmod( issuePMod ) )  // issue a Physical Modification
  Block::add_Modification( std::make_shared< DCNetworkBlockRngdMod >( this ,
               DCNetworkBlockMod::eSetKappa , rng ) ,
                           Observer::par2chnl( issuePMod ) );

 }  // end( DCNetworkBlock::set_kappa( range ) )

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::change_power_flow_limit_constraints(
			   c_Subset & modified_lines , c_ModParam issueAMod )
{
 for( auto i : modified_lines ) {
  double kappa = get_kappa( i );

  if( get_design( i ) ) {
   // Constraints *with* design variables: two per line (LOWER / UPPER).
   // We only update the coefficient of x_i (which is the second variable in
   // the LF).

   // LOWER bound:  F_i - kappa * MinP_i * x_i >= 0
   double lower_coeff_x = -kappa * get_min_power_flow( i );
   auto * lf_low = static_cast< LinearFunction * >(
		   v_power_flow_limit_design_const[ 0 ][ i ].get_function() );
   lf_low->modify_coefficient( 1 , lower_coeff_x , issueAMod );

   // UPPER bound:  F_i - kappa * MaxP_i * x_i <= 0
   double upper_coeff_x = -kappa * get_max_power_flow( i );
   auto * lf_up = static_cast< LinearFunction * >(
		   v_power_flow_limit_design_const[ 1 ][ i ].get_function() );
   lf_up->modify_coefficient( 1 , upper_coeff_x , issueAMod );
   }
  else {
   // Constraints *without* design variable: simple bounds update
   v_power_flow_limit_const[ i ].set_lhs( kappa * get_min_power_flow( i ) ,
					  issueAMod );
   v_power_flow_limit_const[ i ].set_rhs( kappa * get_max_power_flow( i ) ,
					  issueAMod );
    }
  }
 }  // end( DCNetworkBlock::change_power_flow_limit_constraints )

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::change_DC_power_flow_injection_constraints(
                           c_Subset & modified_nodes , c_ModParam issueAMod )
{
 if( ! f_NetworkData->is_HVDC() )
  return;

 for( auto & n : modified_nodes )
  v_power_flow_injection_const[ n ].set_both( -v_ActiveDemand[ n ] );

 }  // end( DCNetworkBlock::change_DC_power_flow_injection_constraints )

/*--------------------------------------------------------------------------*/
/*------------------ METHODS OF DCNetworkBlockSolution ---------------------*/
/*--------------------------------------------------------------------------*/

void DCNetworkBlockSolution::deserialize( const netCDF::NcGroup & group )
{
 // call the method of the base class
 NetworkBlockSolution::deserialize( group );

 // "NumberLines" - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( get_number_nodes() > 1 )
  deserialize_dim( group , "NumberLines" , f_number_lines , false );

 // deserialize the Flow Variables- - - - - - - - - - - - - - - - - - - - - -
 ::deserialize< double >( group , "FlowValue" , v_flow , false );

 // deserialize the Dual Prices - - - - - - - - - - - - - - - - - - - - - - -
 ::deserialize< double >( group , "DualCost" , v_cost , false );

 }  // end( DCNetworkBlockSolution::deserialize( NcGroup & ) )

/*--------------------------------------------------------------------------*/

void DCNetworkBlockSolution::deserialize( const netCDF::NcGroup & group ,
                                          size_t idx )
{
 // call the method of the base class
 NetworkBlockSolution::deserialize( group , idx );

 // "NumberLines" - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( get_number_nodes() > 1 )
  deserialize_dim( group , "NumberLines" , f_number_lines , false );

 std::vector< size_t > strt = { idx , 0 };
 std::vector< size_t > cnt = { 1 , f_number_lines };

 // deserialize the Flow Variables - - - - - - - - - - - - - - - - - - - - -
 auto ncVar = group.getVar( "FlowValue" );
 if( ncVar.isNull() )
  v_flow.clear();
 else {
  v_flow.resize( f_number_lines );
  ncVar.getVar( strt , cnt , v_flow.data() );
  }

 // deserialize the Dual Prices- - - - - - - - - - - - - - - - - - - - - - -
 ncVar = group.getVar( "DualCost" );
 if( ncVar.isNull() )
  v_cost.clear();
 else {
  v_cost.resize( f_number_lines );
  ncVar.getVar( strt , cnt , v_cost.data() );
  }
 }  // end( DCNetworkBlockSolution::deserialize( NcGroup & , size_t ) )

/*--------------------------------------------------------------------------*/

void DCNetworkBlockSolution::read( const Block * block )
{
 // call the method of the base class
 NetworkBlockSolution::read( block );

 auto DCNB = dynamic_cast< const DCNetworkBlock * >( block );
 if( ! DCNB )
  throw( std::invalid_argument(
    "DCNetworkBlockSolution::read: block is not a DCNetworkBlock" ) );

 f_number_lines = DCNB->get_number_lines();

 if( ! v_flow.empty() ) {
  // read the flow power variables - - - - - - - - - - - - - - - - - - - - -
  auto Fl = DCNB->get_power_flow();
  for( Index l = 0 ; l < f_number_lines ; ++l )
   v_flow[ l ] = Fl[ l ].get_value();
  }

 if( ! v_cost.empty() )
  // read the dual prices- - - - - - - - - - - - - - - - - - - - - - - - - -
  DCNB->get_dual_prices( v_cost );

 }  // end( DCNetworkBlockSolution::read )

/*--------------------------------------------------------------------------*/

void DCNetworkBlockSolution::write( Block * block )
{
 // call the method of the base class
 NetworkBlockSolution::write( block );

 auto DCNB = dynamic_cast< DCNetworkBlock * >( block );
 if( ! DCNB )
  throw( std::invalid_argument(
    "DCNetworkBlockSolution::write: block is not a DCNetworkBlock" ) );

 if( f_number_lines != DCNB->get_number_lines() )
  throw( std::invalid_argument(
    "DCNetworkBlockSolution::write: inconsistent lines number" ) );

 // write the flow power variables - - - - - - - - - - - - - - - - - - - - -
 if( ! v_flow.empty() )
  DCNB->set_power_flow( v_flow );

 // write the dual prices- - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_cost.empty() )
  DCNB->set_dual_prices( v_cost );

 }  // end( DCNetworkBlockSolution::write )

/*--------------------------------------------------------------------------*/

void DCNetworkBlockSolution::serialize( netCDF::NcGroup & group ) const
{
 // call the method of the base class
 NetworkBlockSolution::serialize( group );

 // "NumberLines" - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 auto nl = group.addDim( "NumberLines" , f_number_lines );

 // serialize the Flow Variables- - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_flow.empty() )
  ::serialize< double >( group , "FlowValue" , netCDF::NcDouble() , nl ,
       v_flow );

 // serialize the Dual Prices - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_cost.empty() )
  ::serialize< double >( group , "DualCost" , netCDF::NcDouble() , nl ,
       v_cost );

 }  // end( DCNetworkBlockSolution::serialize( NcGroup & ) )

/*--------------------------------------------------------------------------*/

void DCNetworkBlockSolution::serialize( netCDF::NcGroup & group ,
                                        size_t idx ) const {
 // call the method of the base class
 NetworkBlockSolution::serialize( group , idx );

 // now serialize the data structures - - - - - - - - - - - - - - - - - - - -

 netCDF::NcVar FV;  // FlowValue
 netCDF::NcVar DC;  // DualCost

 if( idx == 0 ) {  // first call, have to initialize everything
  // "NumberLines - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  auto nl = group.addDim( "NumberLines" , f_number_lines );

  // "NumberNetworks" is mandatory, and it's checked in the base class
  auto nnw = group.getDim( "NumberNetworks" );

  if( ! v_flow.empty() )
   FV = group.addVar( "FlowValue" , netCDF::NcDouble() , { nnw , nl } );

  if( ! v_cost.empty() )
   DC = group.addVar( "DualCost" , netCDF::NcDouble() , { nnw , nl } );
  }
 else {  // subsequent call, read what is supposedly already there
  if( ! v_flow.empty() )
   FV = group.getVar( "FlowValue" );

  if( ! v_cost.empty() )
   DC = group.getVar( "DualCost" );
  }

 std::vector< size_t > strt = { idx , 0 };
 std::vector< size_t > cnt = { 1 , f_number_lines };

 if( ! FV.isNull() )  // if power flows have to be serialised
  FV.putVar( strt , cnt , v_flow.data() );

 if( ! DC.isNull() )  // if Flow Variables have to be serialised
  DC.putVar( strt , cnt , v_cost.data() );

 }  // end( DCNetworkBlockSolution::serialize( NcGroup & , size_t ) )

/*--------------------------------------------------------------------------*/

DCNetworkBlockSolution * DCNetworkBlockSolution::scale( double factor ) const
{
 // call the method of the base class, which calls clone() and therefore
 // returns a DCNetworkBlockSolution
 auto sol = dynamic_cast< DCNetworkBlockSolution * >(
                                     NetworkBlockSolution::scale( factor ) );
 assert( sol );

 if( factor == 1 )
  return( sol );

 if( ! v_flow.empty() )
  for( Index l = 0 ; l < f_number_lines ; ++l )
   sol->v_flow[ l ] *= factor;

 if( ! v_cost.empty() )
  for( Index l = 0 ; l < f_number_lines ; ++l )
   sol->v_cost[ l ] *= factor;

 return( sol );

 }  // end( DCNetworkBlockSolution::scale )

/*--------------------------------------------------------------------------*/

void DCNetworkBlockSolution::sum( const Solution * solution ,
                                  double multiplier )
{
 // call the method of the base class
 NetworkBlockSolution::sum( solution , multiplier );

 auto DCNBS = dynamic_cast< const DCNetworkBlockSolution * >( solution );
 if( ! DCNBS )
  throw( std::invalid_argument(
    "DCNetworkBlockSolution::sum: solution not a DCNetworkBlockSolution" ) );

 if( f_number_lines != DCNBS->f_number_lines )
  throw( std::invalid_argument(
    "DCNetworkBlockSolution::sum: inconsistent lines number" ) );

 if( ! v_flow.empty() )
  for( Index l = 0 ; l < f_number_lines ; ++l )
   v_flow[ l ] += DCNBS->v_flow[ l ] * multiplier;

 if( ! v_cost.empty() )
  for( Index l = 0 ; l < f_number_lines ; ++l )
   v_cost[ l ] += DCNBS->v_cost[ l ] * multiplier;

 }  // end( DCNetworkBlockSolution::sum )

/*--------------------------------------------------------------------------*/

DCNetworkBlockSolution * DCNetworkBlockSolution::clone( bool empty ) const
{
 auto sol = new DCNetworkBlockSolution();

 if( ! empty ) {
  NetworkBlockSolution::guts_of_clone( sol );

  sol->f_number_lines = f_number_lines;
  sol->v_flow = v_flow;
  sol->v_cost = v_cost;
  }

 return( sol );

 }  // end( DCNetworkBlockSolution::clone )

/*--------------------------------------------------------------------------*/
/*--------------------- End File DCNetworkBlock.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
