/*--------------------------------------------------------------------------*/
/*--------------------- File DCNetworkBlock.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the DCNetworkBlock class.
 *
 * \author Wim van Ackooij \n
 *         EDF R&D OSIRIS \n
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Ali Ghezelsoflu \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Quentin Jacquet \n
 *         EDF R&D OSIRIS \n
 *
 * \author Rafael Durbano Lobato \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Rafael Durbano Lobato
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <map>

#include <utility>

#include "NetworkBlock.h"

#include "DCNetworkBlock.h"

#include "LinearFunction.h"

#include "OneVarConstraint.h"

#include "FRealObjective.h"

#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/IterativeLinearSolvers>

#include <chrono>

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

typedef Eigen::SparseMatrix<double> SpMat;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register DCNetworkBlock to the Block factory
SMSpp_insert_in_factory_cpp_0( DCNetworkBlock );

// register DCNetworkBlockSolution to the Solution factory
SMSpp_insert_in_factory_cpp_0( DCNetworkBlockSolution );

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

// register DCNetworkData to the NetworkData factory

typedef DCNetworkBlock::DCNetworkData DCNetworkData;

SMSpp_insert_in_factory_cpp_0( DCNetworkData );

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS OF DCNetworkData -------------------------*/
/*--------------------------------------------------------------------------*/

void DCNetworkData::deserialize( const netCDF::NcGroup & group )
{
 #ifndef NDEBUG
  static std::vector< std::string > expected_dims = { "NumberNodes" ,
   "NumberLines" ,  // if called from UCBlock:
   "TimeHorizon" , "NumberUnits" , "NumberNetworks" ,
   "NumberElectricalGenerators" };
  check_dimensions( group , expected_dims , std::cerr );

  static std::vector< std::string > expected_vars = { "ActiveDemand" ,
   "StartLine" , "EndLine" , "MinPowerFlow" ,"MaxPowerFlow" ,
   "LineSusceptance" , "NodeSusceptance" , "NetworkCost" , "NodeName" , "LineName" ,
   "ConstantTerm" ,  // if called from UCBlock:
   "ActivePowerDemand" , "GeneratorNode" , "NetworkConstantTerms" ,
   "NetworkBlockClassname" , "NetworkDataClassname"};
  check_variables( group , expected_vars , std::cerr );
 #endif

 NetworkData::deserialize( group );

 // the baseMVA field is a simple scalar value specifying the system MVA base used for converting power into per unit quantities (see Matpower)
  auto gbaseMVA = group.getAtt( "baseMVA" );
  std::string tmp_base;
  if(gbaseMVA.isNull()) f_base_mva = 1.;
  else {
    gbaseMVA.getValues(tmp_base);
    try { f_base_mva = std::stod(tmp_base); }
    catch (...) { f_base_mva = 1.; }
  }
  

 // Optional variables

 if( f_number_nodes > 1 ) {

  deserialize_dim( group , "NumberLines" , f_number_lines , false );

  ::deserialize( group , "StartLine" , f_number_lines , v_start_line , false ,
                 true );

  for( Index i = 0 ; i < f_number_lines ; ++i ) {
   if( ( v_start_line[ i ] < 0 ) || ( v_start_line[ i ] >= f_number_nodes ) )
    throw( std::invalid_argument( "DCNetworkData::deserialize: "
                                  "wrong start node number " +
                                  std::to_string( v_start_line[ i ] ) ) );

  ::deserialize( group , "EndLine" , f_number_lines , v_end_line , false ,
                 true );

   if( ( v_end_line[ i ] < 0 ) || ( v_end_line[ i ] >= f_number_nodes ) )
    throw( std::invalid_argument( "DCNetworkData::deserialize: "
                                  "wrong end node number " +
                                  std::to_string( v_end_line[ i ] ) ) );

   if( v_start_line[ i ] == v_end_line[ i ] )
    throw( std::invalid_argument( "DCNetworkData::deserialize: "
                                  "start node == end node for line " +
                                  std::to_string( v_end_line[ i ] ) ) );
   }

  ::deserialize( group , "MinPowerFlow" , f_number_lines , v_min_power_flow ,
                 true , true );

  ::deserialize( group , "MaxPowerFlow" , f_number_lines , v_max_power_flow ,
                 true , true );

  ::deserialize( group , "NetworkCost" , f_number_lines , v_network_cost ,
                 true , true );

  if( ! deserialize_dim( group, "ReferenceNode", f_reference_node, true ) )
    f_reference_node = 0;
  }

  stored_B2 = SpMat(f_number_nodes - 1, f_number_nodes - 1);
  stored_B2_inv = SpMat(f_number_nodes - 1, f_number_nodes - 1);

  ::deserialize( group , "LineSusceptance" , f_number_lines ,
		 v_line_susceptance , true , true );

  ::deserialize( group , "NodeSusceptance" , f_number_nodes ,
		 v_node_susceptance , true , true );

 ::deserialize( group , "NodeName" , f_number_nodes , v_node_names );
 ::deserialize( group , "LineName" , f_number_lines , v_line_names );

 f_lines_type = -1;

 }  // end( DCNetworkData::deserialize )

 int DCNetworkData::get_reducedIdx( int idx ) {
 if( idx > get_reference_node() )
  return( idx - 1 );
 return( idx );
 }

int DCNetworkData::get_originalIdx( int idx ){
    if ( idx >= get_reference_node() )
      return( idx + 1 );
    return( idx );
 }

 void DCNetworkData::compute_DCDF( const std::vector< Index > & DC_lines, SpMat & PTDF_matrix ){
    const auto number_nodes = get_number_nodes();
    const auto number_lines = get_number_lines();
    const auto & start_line = get_start_line();
    const auto & end_line = get_end_line();

    // linking constraints between AC and HVDC
    SpMat A_DC_transpose( number_nodes - 1 , number_lines );
    for( auto & line_id : DC_lines ) {
      A_DC_transpose.coeffRef( get_reducedIdx( start_line[ line_id ] ) , line_id ) = 1.;
      A_DC_transpose.coeffRef( get_reducedIdx( end_line[ line_id ] ) , line_id ) = -1.; // QJ_TOCHECK 1 or -1 ?
    }
    DCDF = -PTDF_matrix * A_DC_transpose;
    DCDF_was_computed = true;
 }

SpMat DCNetworkData::get_PTDF(const std::vector<Index>& AC_lines, double tikhonov_coeff){
  // get data
  const auto & susceptance = get_line_susceptance();
  const auto number_nodes = get_number_nodes();
  const auto number_lines = get_number_lines();
  if( number_lines <= 0 ) {
    throw( std::logic_error( "DCNetworkBlock::generate_abstract_constraints: "
                             "number of lines of DCNetworkBlock is not set" ) );
  }
  const auto & start_line = get_start_line();
  const auto & end_line = get_end_line();
  
  // construct the matrix using two sub-matrices B_bar and B_hat
  SpMat B_hat = SpMat( number_lines , number_nodes );
  for( auto & line_id : AC_lines ) {
   B_hat.insert( line_id , start_line[ line_id ] ) = susceptance[ line_id ];
   B_hat.insert( line_id , end_line[ line_id ] ) = -susceptance[ line_id ];
  }

  SpMat B_bar = SpMat( number_nodes , number_nodes );
  std::vector< double > B_bar_diag( number_nodes , 0.0 );
  for( auto & line_id : AC_lines ) {
   B_bar.insert( start_line[ line_id ] , end_line[ line_id ] ) = -susceptance[ line_id ];
   B_bar.insert( end_line[ line_id ] , start_line[ line_id ] ) = -susceptance[ line_id ];
   B_bar_diag[ start_line[ line_id ] ] += susceptance[ line_id ];
   B_bar_diag[ end_line[ line_id ] ] += susceptance[ line_id ];
  }
  for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {
   B_bar.insert( node_id , node_id ) = B_bar_diag[ node_id ] + tikhonov_coeff;
  }

  Index ref_node = get_reference_node();
  SpMat I_nref = SpMat( number_nodes , number_nodes - 1 );
  for( Index node_id = 0 ; node_id < ref_node ; ++node_id ) {
   I_nref.insert( node_id , node_id ) = 1.;
  }
  for( Index node_id = ref_node ; node_id < number_nodes - 1 ; ++node_id ) {
   I_nref.insert( node_id + 1 , node_id ) = 1.;
  }

  // construction of B1 and B2 (see documentation)
  SpMat B1 = B_hat * I_nref;
  SpMat B2 = I_nref.transpose() * B_bar * I_nref;

  // compute inverse of B2 and deduce PTDF
  SpMat PTDF_matrix;
  SpMat B2_inv = SpMat( number_nodes - 1 , number_nodes - 1 );
  std::pair< SpMat , SpMat > t = get_stored_B2();
  if( B2.isApprox( t.first ) ) {
   B2_inv = t.second;
   //std::cout << "stored found" << std::endl;
  }
  else {
   //std::cout << "stored NOT found" << std::endl;
   // Inversion of sparse matrix with eigen (solve B2*X = I)
   //Eigen::BiCGSTAB<SpMat> solver;
   Eigen::SparseLU< SpMat > solver;
   solver.compute( B2 );
   SpMat I( number_nodes - 1 , number_nodes - 1 );
   I.setIdentity();
   if( solver.info() != Eigen::Success ) std::cout <<
    "Inversion in PTDF not possible" << std::endl;
   else {
    B2_inv = solver.solve( I );
    set_stored_B2( B2 , B2_inv );
   }
  }
  PTDF_matrix = B1 * B2_inv;
  
  return( PTDF_matrix );
}

/* -----------------------------------------------------------------------*/
void DCNetworkData::compute_cycle_basis(int opt_root){  // QJ: to move to parent class NetworkData ? As it does not require data of neither DC nor AC
   /*Compute a list of cycles which form a basis for cycles of G.

    A basis for cycles of a network is a minimal collection of
    cycles such that any cycle in the network can be written
    as a sum of cycles in the basis.  Here summation of cycles
    is defined as "exclusive or" of the edges. Cycle bases are
    useful, e.g. when deriving equations for electric circuits
    using Kirchhoff's Laws.


    Returns
    -------
    A list of cycle lists.  Each cycle list is a list of nodes
    which forms a cycle (loop) in G.

    Examples
    --------
    >>> G = nx.Graph()
    >>> nx.add_cycle(G, [0, 1, 2, 3])
    >>> nx.add_cycle(G, [0, 3, 4, 5])
    >>> nx.cycle_basis(G, 0)
    [[3, 4, 5, 0], [1, 2, 3, 0]]

    Notes
    -----
    This is adapted from algorithm CACM 491 [1]_.

    References
    ----------
    .. [1] Paton, K. An algorithm for finding a fundamental set of
       cycles of a graph. Comm. ACM 12, 9 (Sept 1969), 514-518.
  */

  if (cycle_basis_was_computed) {
    std::cout << "Cycle basis already computed" << std::endl;
    return;
  }

  // get data
  const auto number_nodes = get_number_nodes();
  const auto number_lines = get_number_lines();
  if( number_lines <= 0 ) {
    throw( std::logic_error( "DCNetworkData::compute_cycle_basis: "
                             "number of lines of DCNetworkBlock is not set" ) );
  }
  const auto & start_line = get_start_line();
  const auto & end_line = get_end_line();

  // First, compute neighbors // QJ: should be a method, if needed in other graph functions ?
  std::vector< std::set< Index > > neighbors(number_nodes, std::set<Index>());
  for (Index id_line = 0; id_line < number_lines; ++id_line){
    Index i = start_line[id_line];
    Index j = end_line[id_line];
    neighbors[i].insert(j);
    neighbors[j].insert(i);
  }

  Index root;
  bool use_root = true;
  if (opt_root < 0) root = get_reference_node();
  else root = opt_root;

  this->v_cycle_basis.clear();
  this->m_spanning_tree.clear();

  auto start_solve = std::chrono::high_resolution_clock::now();

  std::vector<Index> gnodes;
  for (Index i = 0; i < number_nodes; ++i) gnodes.push_back(i);
  while (gnodes.size()){ // loop over connected components
    if (use_root) {
      root = gnodes.back();
      gnodes.pop_back();
    }
    std::vector<Index> stack = {root};
    std::map<Index, Index> pred = {{root,root}};
    std::map<Index, std::set<Index>> used;
    used[root] = std::set<Index>();
    while (stack.size()){ // walk the spanning tree finding cycles
        Index z = stack.back();
        stack.pop_back(); // use last-in so cycles easier to find
        std::set<Index> zused = used[z];
        for (auto& nbr : neighbors[z]){
            if (used.find(nbr) == used.end()){  // new node
                pred[nbr] = z;
                stack.push_back(nbr);
                used[nbr] = std::set<Index>();
                used[nbr].insert(z);
            }
            else if (nbr == z){ // self loops
                this->v_cycle_basis.push_back(std::vector<Index>(1,z));
            }
            else if (zused.find(nbr) == zused.end()){  // found a cycle
                std::set<Index> pn = used[nbr];
                std::vector<Index> cycle = {nbr, z};
                Index p = pred[z];
                while (pn.find(p) == pn.end()){
                    cycle.push_back(p);
                    p = pred[p];
                }
                cycle.push_back(p);
                this->v_cycle_basis.push_back(cycle);
                used[nbr].insert(z);
            }
        }
    }
    for (auto it = pred.begin(); it != pred.end(); ++it){
        auto it_gnode = std::find(gnodes.begin(), gnodes.end(), it->first);
        if (it_gnode != gnodes.end()) gnodes.erase(it_gnode);
    }
    use_root = false; // reinit root
    this->m_spanning_tree.insert(pred.begin(), pred.end());
  }

  double time = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now()
            - start_solve).count()/1000.0;

  cycle_basis_was_computed = true;

  std::cout << "Time to compute cycle basis : " << time << " sec." << std::endl;
}

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS OF DCNetworkBlock ------------------------*/
/*--------------------------------------------------------------------------*/

DCNetworkBlock::~DCNetworkBlock()
{
 Constraint::clear( v_power_flow_limit_const );
 Constraint::clear( v_power_flow_injection_const );
 Constraint::clear( v_power_flow_relax_abs );
 Constraint::clear( v_power_flow_def);

 Constraint::clear( node_injection_bounds_const );

 objective.clear();

 // Delete the DCNetworkData if it is local.
 if( f_local_NetworkData )
  delete( f_NetworkData );
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void DCNetworkBlock::deserialize( const netCDF::NcGroup & group )
{

#ifndef NDEBUG
 static std::vector< std::string > expected_dims = { "NumberNodes" ,
                                                     "NumberLines" };
 check_dimensions( group , expected_dims , std::cerr );

 static std::vector< std::string > expected_vars = { "ActiveDemand" ,
                                                     "StartLine" ,
                                                     "EndLine" ,
                                                     "MinPowerFlow" ,
                                                     "MaxPowerFlow" ,
                                                     "LineSusceptance" ,
                                                     "NodeSusceptance" ,
                                                     "NetworkCost" ,
                                                     "NodeName" ,
                                                     "LineName" ,
                                                     "ConstantTerm" };
 check_variables( group , expected_vars , std::cerr );
#endif

 NetworkBlock::deserialize( group );

 // Optional variables

 Index NumberNodes;
 if( deserialize_dim( group , "NumberNodes" , NumberNodes ) ) {
  // Since the dimension "NumberNodes" has been provided, it means that a
  // DCNetworkData has been provided. Thus, the DCNetworkData is deserialized,
  // and it is marked as being local
  if( f_local_NetworkData )
   // if the NetworkData has not been passed from UCBlock, then delete it
   delete( f_NetworkData );
  auto DCND = new DCNetworkData();
  DCND->deserialize( group );
  if( f_NetworkData &&
    ( f_NetworkData->get_number_nodes() != DCND->get_number_nodes() ) )
   throw( std::logic_error(
    "DCNetworkBlock::deserialize: NumberNodes not matching between NetworkData" ) );
  f_NetworkData = DCND;
  f_local_NetworkData = true;
  // A DCNetworkData has been provided. So, the size of the given vector of
  // active demand must be equal to the number of nodes.
  ::deserialize( group , "ActiveDemand" , NumberNodes , v_ActiveDemand );
 } else {
  // A DCNetworkData has not been provided. However, the active demand may still
  // have been provided.

  auto ActiveDemand = group.getVar( "ActiveDemand" );

  if( ! ActiveDemand.isNull() ) {
   // The active demand has indeed been provided.

   if( ActiveDemand.getDimCount() != 1 )
    // The active demand must be a one-dimensional array.
    throw( std::invalid_argument(
     "DCNetworkBlock::deserialize(): ActiveDemand should have one dimension, "
     "but it has " + std::to_string( ActiveDemand.getDimCount() ) ) );

   // Retrieve the number of nodes from the size of the given netCDF variable.
   const auto number_nodes = ActiveDemand.getDim( 0 ).getSize();

   // Resize the vector of active demand.
   v_ActiveDemand.resize( number_nodes );

   // Retrieve the active demand from the netCDF variable.
   ActiveDemand.getVar( v_ActiveDemand.data() );
  }
 }
}  // end( DCNetworkBlock::deserialize )

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::generate_abstract_variables( Configuration * stvv )
{
 if( variables_generated() )  // variables have already been generated
  return;                     // nothing to do

  //TEMP: formulation choice
 ftype = CYCLE; // should be an option somewhere else

 NetworkBlock::generate_abstract_variables( stvv );

 if (ftype == PTDF) {
  generate_PTDF_variables(stvv); 
 }
 else if (ftype == CYCLE){
  generate_CYCLE_variables(stvv);
 }
 else {
  throw( std::logic_error( "Not Implemented yet" ) );
 }

 set_variables_generated();

}  // end( DCNetworkBlock::generate_abstract_variables )


/*--------------------------------------------------------------------------*/
void DCNetworkBlock::generate_PTDF_variables( Configuration * stvv )
{
  /**
   * This formulation corresponds to the "PTDF + FLOW" formulation of
   * "Linear Optimal Power Flow Using Cycle Flows" of
   *    Jonas Horsch, Henrik Ronellenfitsch, Dirk Witthaut, Tom Brown 
   * */
 const auto number_lines = get_number_lines();

 if( number_lines > 0 ) {
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
 }
}  // end( DCNetworkBlock::generate_PTDF_variables )

/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/
void DCNetworkBlock::generate_CYCLE_variables( Configuration * stvv )
{
  /**
  * Implementation of "Linear Optimal Power Flow Using Cycle Flows" of
  *    Jonas Horsch, Henrik Ronellenfitsch, Dirk Witthaut, Tom Brown 
  * 
  * Here, we opt for the "CYCLE + FLOW" formulation with
  *   - variables "v_power_flow" as in the PTDF formulation (f_l in the paper)
  *   - variables "v_cycle_flow" (h_c in the paper)
  * 
  */

  generate_PTDF_variables(stvv); // we have the same variables + others

  const auto number_nodes = get_number_nodes();
  if( number_nodes <= 1 )
   return;
  const auto number_lines = get_number_lines();
  
  if( number_lines > 0 && number_nodes > 0) {
   // the power flow variable on cycle basis
   v_cycle_flow.resize( number_lines - number_nodes + 1); 
      // we know the number of cycles by the graph theory, see the paper. 
      // So, no reason to call get_lines_in_cycle()
   for( auto & var : v_cycle_flow )
    var.set_type( ColVariable::kContinuous );
   add_static_variable( v_cycle_flow , "cycle_flow_network" );
  }
}  // end( DCNetworkBlock::generate_CYCLE_variables )

/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/

void DCNetworkBlock::generate_abstract_constraints( Configuration * stcc )
{
 if( constraints_generated() )  // constraints have already been generated
  return;                       // nothing to do

 if (ftype == PTDF) {
  generate_PTDF_constraints(stcc);
 }
 else if (ftype == CYCLE){
  generate_CYCLE_constraints(stcc);
 }
 else {
  throw( std::logic_error( "Not Implemented yet" ) );
 }

 set_constraints_generated();

}  // end( DCNetworkBlock::generate_abstract_constraints )


/*--------------------------------------------------------------------------*/
void DCNetworkBlock::generate_CYCLE_constraints( Configuration * stcc )
{
  /**
   * Implementation of "Linear Optimal Power Flow Using Cycle Flows" of
   *    Jonas Horsch, Henrik Ronellenfitsch, Dirk Witthaut, Tom Brown 
   */

 const auto number_nodes = get_number_nodes();
 if( number_nodes <= 1 )
  return;
 const auto number_lines = get_number_lines();

 if( number_lines <= 0 )
  throw( std::logic_error( "DCNetworkBlock::generate_abstract_constraints: "
                           "number of lines of DCNetworkBlock is not set" ) );

 const auto & start_line = f_NetworkData->get_start_line();
 const auto & end_line = f_NetworkData->get_end_line();
 const auto lines_type = f_NetworkData->get_lines_type();

  // First step: compute the cycle basis and spanning tree
  //std::cout << "Cycle basis:" << std::endl;
  auto basis = f_NetworkData->get_lines_in_cycles(); // cycle incidence matrices C_{lc} in the paper
  /* for Debug
  for (auto& cycle: basis){
    std::cout << "(" ;
    for (auto it = cycle.begin(); it != cycle.end(); ++it){
      std::cout << "line " << it->first << ": " << it->second << ",";
    }
    std::cout << ")" << std::endl;
  }*/

  //std::cout << "Spanning tree:" << std::endl;
  auto tree = f_NetworkData->get_lines_in_spanning_tree();
  /* for debug
  for (auto it = tree.begin(); it != tree.end(); ++it){
    std::cout << "(line " << it->first << ":" << it->second << "),";
  }
  std::cout << std::endl;
  */

  // Then, create the equations


  throw( std::logic_error( "Not Implemented yet" ) );
 }

/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
void DCNetworkBlock::generate_PTDF_constraints( Configuration * stcc )
{
 // In mixed mode we can write all nodal balances, instead of just strictly those needed;
 // set the following flag to true in that case
 double nodal_slack = 0.05;
 double ptdf_slack  = 1.0; 
 double ptdf_round  = 1e-7; 
 bool full_formulation = true;  
 if ( full_formulation )
    std::cout << "begin constraint extended" << std::endl;
 else
    std::cout << "begin constraint economic" << std::endl;

 const auto number_nodes = get_number_nodes();

 if( number_nodes <= 1 )
  return;

 const auto number_lines = get_number_lines();

 if( number_lines <= 0 )
  throw( std::logic_error( "DCNetworkBlock::generate_abstract_constraints: "
                           "number of lines of DCNetworkBlock is not set" ) );

 const auto & start_line = f_NetworkData->get_start_line();
 const auto & end_line = f_NetworkData->get_end_line();
 const auto lines_type = f_NetworkData->get_lines_type();

 // Splitting AC and DC part
 std::vector<Index> AC_lines = f_NetworkData->get_AC_lines();
 SpMat PTDF_matrix = f_NetworkData->get_PTDF(AC_lines);

 /*std::cout << "Debut ecriture" << std::endl;
 std::cout << "Line" << "\t" << "Node" << "Coefficient" << std::endl; 
 for( auto & line_id : AC_lines ) {
    Eigen::SparseMatrix<double> a_row = PTDF_matrix.block(line_id, 0, 1, PTDF_matrix.cols() );
    for (int k=0; k < a_row.outerSize(); ++k){
        for (Eigen::SparseMatrix<double>::InnerIterator it(a_row,k); it; ++it){
          int node_id = f_NetworkData->get_originalIdx( it.col() );
          if( node_id != f_NetworkData->get_reference_node() ) {
            double coefficient = round_to( it.value(), 1e-7 ); //PTDF_matrix.coeff( line_id , f_NetworkData->get_reducedIdx( node_id ) );
            std::cout << f_NetworkData->get_line_names()[line_id] << "\t" <<  f_NetworkData->get_node_names()[node_id] << "\t" << coefficient << std::endl;
          }
        }
      } // for each node
    }
 std::cout << "Fichier fermé" << std::endl;
 exit(0);*/
 std::vector<Index> DC_lines = f_NetworkData->get_DC_lines();

 // ===== auxiliary variables for nonempty cost
 if( ! f_NetworkData->get_network_cost().empty() ) {
   // 0 <= V_l - F_l && 0 <= V_l + F_l
   v_power_flow_relax_abs.resize(
    boost::multi_array< FRowConstraint , 2 >::extent_gen()[ 2 ][ number_lines ] );

   // Definition of absolute value of power flow
   for( Index line_id = 0 ; line_id < number_lines ; ++line_id ){
    auto lfunc_1 = new LinearFunction();
    lfunc_1->add_variable( &v_power_flow[ line_id ] , -1.0 );
    lfunc_1->add_variable( &v_auxiliary_variable[ line_id ] , 1.0 );
    v_power_flow_relax_abs[ 0 ][ line_id ].set_lhs( 0.0 );
    v_power_flow_relax_abs[ 0 ][ line_id ].set_rhs( Inf< double >() );
    v_power_flow_relax_abs[ 0 ][ line_id ].set_function( lfunc_1 );

    auto lfunc_2 = new LinearFunction();
    lfunc_2->add_variable( &v_power_flow[ line_id ] , 1.0 );
    lfunc_2->add_variable( &v_auxiliary_variable[ line_id ] , 1.0 );
    v_power_flow_relax_abs[ 1 ][ line_id ].set_lhs( 0.0 );
    v_power_flow_relax_abs[ 1 ][ line_id ].set_rhs( Inf< double >() );
    v_power_flow_relax_abs[ 1 ][ line_id ].set_function( lfunc_2 );
   }
   add_static_constraint( v_power_flow_relax_abs , "power_flow_relax_abs" );
 } // ===== end( cost not empty )

 // ===== constraints on the DC part
 if( lines_type == kHVDC || lines_type == kAC_HVDC ) {

  // Flow limit constraints
  v_power_flow_limit_const.resize( number_lines );

  for( Index line_id = 0 ; line_id < number_lines ; ++line_id  ) {
    const auto kappa = get_kappa( line_id );
    v_power_flow_limit_const[ line_id ].set_lhs( kappa * f_NetworkData->get_min_power_flow( line_id ) );
    v_power_flow_limit_const[ line_id ].set_rhs( kappa * f_NetworkData->get_max_power_flow( line_id ) );
    v_power_flow_limit_const[ line_id ].set_variable( & v_power_flow[ line_id ] );
    add_static_constraint( v_power_flow_limit_const, "Power_flow_limit" );
  }

  // Power flow and node injection constraints
  if( lines_type == kHVDC ){
   v_power_flow_injection_const.resize(number_nodes );

    for( Index n = 0 ; n < number_nodes ; ++n ) {
        auto lfunc = new LinearFunction();
        lfunc->add_variable( & v_node_injection[ 0 ][ n ] , -1.0 );

        for( auto & line_id : DC_lines ) {
          if( start_line[ line_id ] == n )  lfunc->add_variable( & v_power_flow[ line_id ] , 1.0 );
          if( end_line[ line_id ] == n )    lfunc->add_variable( & v_power_flow[ line_id ] , -1.0 );
        }
        if( lines_type == kHVDC ) {
          v_power_flow_injection_const[ n ].set_both( -v_ActiveDemand[ n ] );
          v_power_flow_injection_const[ n ].set_function( lfunc );
        }
      }
    add_static_constraint( v_power_flow_injection_const, "HVDC_power_flow_injection" );
  }
  // If we have mixed lines, we have as many as nodes impacted and touched by DC lines
  if( lines_type == kAC_HVDC ){
    int nb_DCnodes = 0;
    // Savagely setting all visited nodes to true will generate nodal balances for all nodes
    // the default and subtle initialization should be with false
    std::vector<bool> nodes_vist( number_nodes, full_formulation );
  
    // Flip any visited nodes to true
    for( auto & line_id : DC_lines ){
        nodes_vist[ start_line[ line_id ] ] = true;
        nodes_vist[ end_line[ line_id ] ] = true;
    }
  
    for( Index n = 0 ; n < number_nodes ; ++n ){
      if ( nodes_vist[n] )
        ++nb_DCnodes;
    }
    v_AC_HVDC_power_flow_const.resize( nb_DCnodes );

    // Add power balance equations for impacted nodes
    int iDCnode=0;
    for( Index n = 0 ; n < number_nodes ; ++n ){
      if ( nodes_vist[n] ){
        auto lfunc = new LinearFunction();
        lfunc->add_variable( & v_node_injection[ 0 ][ n ] , -1.0 );

        for( auto & line_id : DC_lines ) {
          if( start_line[ line_id ] == n )  lfunc->add_variable( & v_power_flow[ line_id ] , 1.0 );
          if( end_line[ line_id ] == n )    lfunc->add_variable( & v_power_flow[ line_id ] , -1.0 );
        }
        for( auto & line_id : AC_lines ) {
          if( start_line[ line_id ] == n )  lfunc->add_variable( & v_power_flow[ line_id ] , 1.0 );
          if( end_line[ line_id ] == n )    lfunc->add_variable( & v_power_flow[ line_id ] , -1.0 );
        }
        // v_AC_HVDC_power_flow_const[ iDCnode ].set_both( -1.0*v_ActiveDemand[ n ] );
        v_AC_HVDC_power_flow_const[ iDCnode ].set_lhs( -1.0*v_ActiveDemand[ n ] - nodal_slack ); // allow for a 0.01 MW deviation
        v_AC_HVDC_power_flow_const[ iDCnode ].set_rhs( -1.0*v_ActiveDemand[ n ] + nodal_slack );
        v_AC_HVDC_power_flow_const[ iDCnode ].set_function( lfunc );
       
        ++iDCnode; // update the index
      }
    }
    // Add the whole vector of constraints at once
    add_static_constraint(v_AC_HVDC_power_flow_const, "ACdHVDC_power_flow_injection");
  }  

 } // ===== end constraints on HVDC part

 // ===== constraints on AC Part
 if( lines_type == kAC || lines_type == kAC_HVDC ) {
  
  if( lines_type == kAC_HVDC ) {
   if ( !f_NetworkData->was_DCDF_computed() ) 
      f_NetworkData->compute_DCDF( DC_lines, PTDF_matrix );     
  }

  // 
  // const auto & l_names = f_NetworkData->get_line_names();
  // const auto & n_names = f_NetworkData->get_node_names();
  // std::cout << "Some output for line " << l_names[1045] << " ";
  // Eigen::SparseMatrix<double> a_row = PTDF_matrix.block(1045, 0, 1, PTDF_matrix.cols() );
  // for (int k=0; k < a_row.outerSize(); ++k){
	  //for (Eigen::SparseMatrix<double>::InnerIterator it(a_row,k); it; ++it){
		  //std::cout << "(" << it.row() << "," << it.col()+1 << " name=" << n_names[it.col()+1] << " ) = " << it.value();
//
 //   }
  //}
  //std::cout << "\n";

  // Flow limit constraints
  v_power_flow_def.resize( number_lines );
  for( auto & line_id : AC_lines ) {
    // TODO : verify if this does not entail a copy of the information which would be inefficient
    Eigen::SparseMatrix<double> a_row = PTDF_matrix.block(line_id, 0, 1, PTDF_matrix.cols() );
   
    auto lfunc = new LinearFunction();
    lfunc->add_variable( &v_power_flow[ line_id ] , -1.0 );
    double constant_term = 0;

    // loop over all nodes
    // for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {   
    for (int k=0; k < a_row.outerSize(); ++k){
      for (Eigen::SparseMatrix<double>::InnerIterator it(a_row,k); it; ++it){
        int node_id = f_NetworkData->get_originalIdx( it.col() );
        if( node_id != f_NetworkData->get_reference_node() ) {
          double coefficient = round_to( it.value(), ptdf_round ); //PTDF_matrix.coeff( line_id , f_NetworkData->get_reducedIdx( node_id ) );
          // Distribution Factor Matrix
          lfunc->add_variable( &v_node_injection[ 0 ][ node_id ] , coefficient );
          constant_term += coefficient * v_ActiveDemand[ node_id ];
        }
      }
    } // for each node
    
    if( lines_type == kAC_HVDC ) {
      SpMat DCDF_ = f_NetworkData->get_DCDF();
      // TODO : Also only loop over the non zero entries of DCDF only ...
      for( auto & dc_line_id : DC_lines ){
        double coeff = round_to(DCDF_.coeff( line_id , dc_line_id ), ptdf_round );
        lfunc->add_variable( &v_power_flow[ dc_line_id ] , coeff );
      }
    }

    // Set the constraint (AC)
    v_power_flow_def[ line_id ].set_function( lfunc );
    v_power_flow_def[ line_id ].set_lhs( constant_term - ptdf_slack );
    v_power_flow_def[ line_id ].set_rhs( constant_term + ptdf_slack );
  }

  add_static_constraint( v_power_flow_def, "AC/HVDC_powerflow_def" );
  
  auto lfunc = new LinearFunction();
  double constant_term = 0.;
  for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {
   lfunc->add_variable( &v_node_injection[ 0 ][ node_id ] , 1. );
   constant_term += v_ActiveDemand[ node_id ];
  }
  overall_balanced_const.set_function( lfunc );
  overall_balanced_const.set_lhs( constant_term );
  overall_balanced_const.set_rhs( constant_term );
  add_static_constraint( overall_balanced_const , "overall_balanced_const" );

 } // ===== end AC and AC-HVDC constraints

 set_constraints_generated();

}  // end( DCNetworkBlock::generate_PTDF_constraints )

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::generate_objective( Configuration * objc )
{
 if( objective_generated() )  // Objective has already been generated
  return;                     // nothing to do

 auto lf = new LinearFunction();

 if( ! f_NetworkData->get_network_cost().empty() )
  for( Index line_id = 0 ; line_id < get_number_lines() ; ++line_id )
   lf->add_variable( &v_auxiliary_variable[ line_id ] ,
                     f_NetworkData->get_network_cost()[ line_id ] ,
                     eDryRun );

 lf->set_constant_term( f_ConstTerm );

 objective.set_function( lf );
 objective.set_sense( Objective::eMin );

 // Set Block objective
 this->set_objective( &objective );

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
 }

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
  if( auto tc = dynamic_cast< SimpleConfiguration< std::pair< double , int > > * >( c ) ) {
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
  // Constraints
  && RowConstraint::is_feasible( v_power_flow_limit_const , tol , rel_viol )
  && RowConstraint::is_feasible( v_power_flow_injection_const , tol , rel_viol )
  && RowConstraint::is_feasible( v_power_flow_def , tol , rel_viol )
  && RowConstraint::is_feasible( v_power_flow_relax_abs , tol , rel_viol )
  && RowConstraint::is_feasible( node_injection_bounds_const , tol , rel_viol ) );
} // end( DCNetworkBlock::is_feasible )

/*--------------------------------------------------------------------------*/
/*--------- METHODS FOR LOADING, PRINTING & SAVING THE DCNetworkBlock ------*/
/*--------------------------------------------------------------------------*/

void DCNetworkData::serialize( netCDF::NcGroup & group ) const
{
 NetworkData::serialize( group );

 if( f_number_nodes > 1 ) {
  auto NumberNodes = group.getDim( "NumberNodes" );
  auto NumberLines = group.addDim( "NumberLines" );

  ::serialize( group , "StartLine" , netCDF::NcUint() , NumberLines ,
               v_start_line );

  ::serialize( group , "EndLine" , netCDF::NcUint() , NumberLines ,
               v_end_line );

  ::serialize( group , "MinPowerFlow" , netCDF::NcDouble() , NumberLines ,
               v_min_power_flow );

  ::serialize( group , "MaxPowerFlow" , netCDF::NcDouble() , NumberLines ,
               v_max_power_flow );

  ::serialize( group , "NetworkCost" , netCDF::NcDouble() , NumberLines ,
               v_network_cost );

  if( ! v_line_names.empty() ) {
   assert( v_line_names.size() == NumberLines.getSize() );
   auto LineName = group.addVar( "LineName" , netCDF::NcString() ,
				 NumberLines );
   for( Index i = 0 ; i < v_line_names.size() ; ++i )
    LineName.putVar( { i } , v_line_names[ i ] );
   }
  }

 if( ! v_node_names.empty() ) {
  auto NumberNodes = group.getDim( "NumberNodes" );
  assert( v_node_names.size() == NumberNodes.getSize() );
  auto NodeName = group.addVar( "NodeName" , netCDF::NcString() ,
				NumberNodes );
  for( Index i = 0 ; i < v_node_names.size() ; ++i )
   NodeName.putVar( { i } , v_node_names[ i ] );
  }
 }  // end( DCNetworkData::serialize )

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::serialize( netCDF::NcGroup & group ) const
{
 NetworkBlock::serialize( group );

 if( auto network_data = get_NetworkData() )
  // If a DCNetworkData is present, serialize it.
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
  if( std::all_of( values ,
                   values + subset.size() ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  v_ActiveDemand.resize( get_number_nodes() );
  }

 bool identical = true;
 for( auto i : subset ) {
  if( i >= v_ActiveDemand.size() )
   throw( std::invalid_argument( "DCNetworkBlock::set_active_demand: "
                                 "invalid value in subset." ) );

  auto demand = *(values++);
  if( v_ActiveDemand[ i ] != demand ) {
   identical = false;

   if( not_dry_run( issuePMod ) )
    // Change the physical representation
    v_ActiveDemand[ i ] = demand;
   }
  }

 if( identical )
  return;  // nothing changes; return

 if( not_dry_run( issuePMod ) &&
      not_dry_run( issueAMod ) &&
      constraints_generated() ) {
      // Change the abstract representation
      if( f_NetworkData->get_lines_type() == kHVDC ) {
          std::vector< Index > modified_nodes(subset.begin(), subset.end());
          change_DC_power_flow_injection_constraints( modified_nodes , issueAMod );
      }
      if( f_NetworkData->get_lines_type() == kAC || f_NetworkData->get_lines_type() == kAC_HVDC ) {
          std::vector< Index > modified_lines( f_NetworkData->get_number_lines());
          std::iota( modified_lines.begin(), modified_lines.end() , 0 );
          change_relax_abs_constraints(modified_lines, issueAMod ); // all lines are modified
          change_power_flow_limit_constraints(modified_lines, issueAMod );
      }
  }

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification
  if( ! ordered )
   std::sort( subset.begin() , subset.end() );

  Block::add_Modification( std::make_shared< NetworkBlockSbstMod >(
                            this , NetworkBlockMod::eSetActD , std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );
 }
}  // end( DCNetworkData::set_active_demand )

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
   // Change the physical representation

  std::copy( values , values + ( rng.second - rng.first ) ,
             v_ActiveDemand.begin() + rng.first );

  if( not_dry_run( issueAMod ) && constraints_generated() ) {
       // Change the abstract representation

       if( f_NetworkData->get_lines_type() == kHVDC ) {
          std::vector< Index > modified_nodes( rng.second - rng.first );
          std::iota( modified_nodes.begin() , modified_nodes.end() ,
		     rng.first );
          change_DC_power_flow_injection_constraints( modified_nodes ,
						      issueAMod );
        }
       if( f_NetworkData->get_lines_type() == kAC ||
	   f_NetworkData->get_lines_type() == kAC_HVDC ) {
	std::vector< Index > modified_lines( f_NetworkData->get_number_lines());
	std::iota( modified_lines.begin() , modified_lines.end(),0);
	change_relax_abs_constraints( modified_lines , issueAMod ); // all lines are modified
	change_power_flow_limit_constraints( modified_lines , issueAMod );
        }
    }
  }
 }

 if( issue_pmod( issuePMod ) )  // issue a Physical Modification
  Block::add_Modification( std::make_shared< NetworkBlockRngdMod >( this ,
					 NetworkBlockMod::eSetActD , rng ) ,
                           Observer::par2chnl( issuePMod ) );

 }  // end( DCNetworkData::set_active_demand )

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::set_kappa( MF_dbl_it values , Subset && subset ,
                                bool ordered ,
                                c_ModParam issuePMod , c_ModParam issueAMod )
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

  const auto kappa = *(values++);
  if( v_kappa[ i ] != kappa ) {
   identical = false;
   if( not_dry_run( issuePMod ) )
    // Change the physical representation
    v_kappa[ i ] = kappa;
   }
  }

 if( identical )
  return;  // nothing changes; return

 if( not_dry_run( issuePMod ) &&
     not_dry_run( issueAMod ) &&
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
 }  // end( DCNetworkData::set_kappa( subset ) )

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

 }  // end( DCNetworkData::set_kappa( range ) )

/*--------------------------------------------------------------------------*/
void DCNetworkBlock::change_power_flow_limit_constraints
(const std::vector<Index>& modified_lines, c_ModParam issueAMod){
  for( auto i : modified_lines ) {
   v_power_flow_limit_const[ i ].set_lhs
    ( v_kappa[ i ] * f_NetworkData->get_min_power_flow( i ) , issueAMod );
   v_power_flow_limit_const[ i ].set_rhs
    ( v_kappa[ i ] * f_NetworkData->get_max_power_flow( i ) , issueAMod );
  }
}

/*--------------------------------------------------------------------------*/
void DCNetworkBlock::change_relax_abs_constraints
(const std::vector< Index > & modified_lines , c_ModParam issueAMod ) {
 if( ( f_NetworkData->get_lines_type() == kAC ) ||
     ( f_NetworkData->get_lines_type() == kAC_HVDC ) ) {
  std::vector< Index > AC_lines = f_NetworkData->get_AC_lines();
  SpMat PTDF_matrix = f_NetworkData->get_PTDF( AC_lines );
  for( auto & i : modified_lines ) {
   double constant_term = 0;
   for( Index node_id = 0 ; node_id < get_number_nodes() ; ++node_id )
    constant_term -=
     PTDF_matrix.coeff( i , node_id ) * v_ActiveDemand[ node_id ];
   v_power_flow_relax_abs[ 0 ][ i ].set_lhs( constant_term );
   v_power_flow_relax_abs[ 1 ][ i ].set_lhs( -constant_term );
  }
 }
}

/*--------------------------------------------------------------------------*/

void DCNetworkBlock::change_DC_power_flow_injection_constraints(
         const std::vector< Index > & modified_nodes , c_ModParam issueAMod )
{
 if( f_NetworkData->get_lines_type() != kHVDC )
  return;

 for( auto & n : modified_nodes )
  v_power_flow_injection_const[ n ].set_both( -v_ActiveDemand[ n ] );
 }

/*--------------------------------------------------------------------------*/
/*------------------ METHODS OF DCNetworkBlockSolution ---------------------*/
/*--------------------------------------------------------------------------*/

void DCNetworkBlockSolution::deserialize( const netCDF::NcGroup & group )
{
 // call the method of the base class
 NetworkBlockSolution::deserialize( group );

 // "NumberLines" is mandatory- - - - - - - - - - - - - - - - - - - - - - - -
 deserialize_dim( group , "NumberLines" , f_number_lines , false );

 // deserialize the Flow Variables - - - - - - - - - - - - - - - - - - - - -
 ::deserialize< double >( group , "FlowValue" , v_flow , false );

 // deserialize the Dual Prices- - - - - - - - - - - - - - - - - - - - - - -
 ::deserialize< double >( group , "DualCost" , v_cost , false );

 }  // end( DCNetworkBlockSolution::deserialize )

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

 if( ! v_flow.empty() ) {
  // write the flow power variables- - - - - - - - - - - - - - - - - - - - -
  auto Fl = DCNB->get_power_flow();
  for( Index l = 0 ; l < f_number_lines ; ++l )
   Fl[ l ].set_value( v_flow[ l ] );
  }

 if( ! v_cost.empty() )
  // write the dual prices - - - - - - - - - - - - - - - - - - - - - - - - -
  DCNB->set_dual_prices( v_cost );

 }  // end( DCNetworkBlockSolution::write )

/*--------------------------------------------------------------------------*/

void DCNetworkBlockSolution::serialize( netCDF::NcGroup & group ) const
{
 // call the method of the base class
 NetworkBlockSolution::serialize( group );

 // "NumberLines" is mandatory- - - - - - - - - - - - - - - - - - - - - - - -
 auto nl = group.addDim( "NumberLines" , f_number_lines );

 // serialize the Flow Variables- - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_flow.empty() )
  ::serialize< double >( group , "FlowValue" , netCDF::NcDouble() , nl ,
			 v_flow );

 // serialize the Dual Prices - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_cost.empty() )
  ::serialize< double >( group , "DualCost" , netCDF::NcDouble() , nl ,
			 v_cost );

 }  // end( DCNetworkBlockSolution::serialize )

/*--------------------------------------------------------------------------*/

DCNetworkBlockSolution * DCNetworkBlockSolution::scale( double factor ) const
{
 // call the method of the base class
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
