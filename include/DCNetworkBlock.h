/*--------------------------------------------------------------------------*/
/*--------------------------- File DCNetworkBlock.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 *
 * Header file for class DCNetworkBlock, which derives from NetworkBlock and
 * defines the standard linear constraints corresponding to the "DC model"
 * of the transmission network in the Unit Commitment problem.
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
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __DCNetworkBlock
 #define __DCNetworkBlock
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Block.h"

#include "LinearFunction.h"

#include "FRowConstraint.h"

#include "OneVarConstraint.h"

#include "NetworkBlock.h"

#include "FRealObjective.h"

#include <Eigen/Sparse>

#include <algorithm>

#include <utility>

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)

namespace SMSpp_di_unipi_it
{
 using SpMat = Eigen::SparseMatrix< double >;

/*--------------------------------------------------------------------------*/
/*------------------------- CLASS DCNetworkBlock ---------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// a "DC" transmission NetworkBlock
/** DCNetworkBlock derives from NetworkBlock and defines the linear
 * constraints of the transmission network at one instant of a Unit Commitment
 * problem under the DC approximation of the power flow. In this approximation
 * the voltage magnitudes are at their nominal value, the lines have no
 * resistance and the differences of the voltage angles are small enough for
 * their sine to be replaced by the angle; hence, the flow on a line that
 * obeys Kirchhoff's laws is a linear function of the angles at its two ends,
 * and no power is lost on it. A DCNetworkBlock is the network of one instant
 * (or of the instants that share it) of the model described in
 * \ref ucblock_model.
 *
 * \par Notation
 * Let \f$ \mathcal{N} \f$ be the set of the nodes, \f$ N = |\mathcal{N}| \f$,
 * and \f$ \mathcal{L} \f$ the set of the lines. Each line \f$ l \f$ has a
 * start node \f$ s(l) \f$, an end node \f$ e(l) \f$ (several end nodes
 * \f$ e_j(l) \f$ for a hyperarc, see DCNetworkData::deserialize()) and a flow
 * \f$ F_l \f$, positive when power leaves \f$ s(l) \f$ towards the end
 * node(s). We split the lines into the set \f$ \mathcal{L}^{S} \f$ of those
 * with a nonzero susceptance \f$ \mathfrak{S}_l \f$, whose flow is fixed by
 * the angles of their end nodes, and the set \f$ \mathcal{L}^{H} \f$ of those
 * with zero susceptance, whose flow is fully controllable. In the methods of
 * this class the former are called "DC lines" (e.g.,
 * DCNetworkData::get_DC_lines()), since the DC approximation applies to them,
 * and the latter "HVDC lines" (DCNetworkData::get_HVDC_lines()). Of the flow
 * of a line in \f$ \mathcal{L}^{H} \f$ only the fraction \f$ \eta_l \f$ (the
 * efficiency, \f$ \eta_{l,j} \f$ for the branch \f$ j \f$ of a hyperarc)
 * reaches the end node, while the efficiency of a line in
 * \f$ \mathcal{L}^{S} \f$ is 1 whatever the data say. For each node
 * \f$ n \f$, \f$ S_n \f$ is the node injection variable that UCBlock links to
 * the units at the node, \f$ D^{ac}_n \f$ the active demand ("ActiveDemand")
 * and \f$ a_n = S_n - D^{ac}_n \f$ the net injection. Furthermore, the HVDC
 * incidence matrix \f$ A^{H} \f$, of size \f$ |\mathcal{L}^{H}| \times N \f$,
 * has \f$ A^{H}_{l,s(l)} = 1 \f$, \f$ A^{H}_{l,e_j(l)} = -\eta_{l,j} \f$ and
 * 0 elsewhere; with \f$ F^{H} \f$ the vector of the flows of the lines in
 * \f$ \mathcal{L}^{H} \f$, the vector \f$ a - ( A^{H} )^\top F^{H} \f$ is the
 * injection that the lines in \f$ \mathcal{L}^{S} \f$ have to carry. Hence,
 * the node balance of node \f$ n \f$ is
 * \f[
 *   - S_n + \sum_{ l \in \mathcal{L} : s(l) = n } F_l
 *         - \sum_{ l \in \mathcal{L} , j : e_j(l) = n } \eta_{l,j} F_l
 *   = - D^{ac}_n \; ,
 *   \tag{B}
 * \f]
 * the sums running over all the lines incident to \f$ n \f$: the net
 * injection of the node equals the power that leaves it along the lines
 * minus the power that reaches it. A network is pure HVDC when
 * \f$ \mathcal{L}^{S} = \emptyset \f$, pure DC when \f$ \mathcal{L}^{H} =
 * \emptyset \f$ and mixed otherwise (see DCNetworkData::is_HVDC(),
 * DCNetworkData::is_DC() and DCNetworkData::is_DC_HVDC()). A hyperarc is
 * always in \f$ \mathcal{L}^{H} \f$, since the voltage law relates the
 * angles of two nodes: its flow is a decision independent of the angles,
 * and it enters every formulation only through its row of \f$ A^{H} \f$,
 * i.e., as the injection \f$ - F_l \f$ at \f$ s(l) \f$ and
 * \f$ \eta_{l,j} F_l \f$ at each \f$ e_j(l) \f$, as a line in
 * \f$ \mathcal{L}^{H} \f$ with one end node does. Hence, hyperarcs and
 * lines with nonzero susceptance can be in the same network, the latter
 * making a meshed grid whose angles the hyperarcs do not touch.
 *
 * \par Components and references
 * Removing the lines in \f$ \mathcal{L}^{H} \f$ splits the network into the
 * components \f$ \mathcal{N}_1 , \ldots , \mathcal{N}_{n_c} \f$ that the
 * lines in \f$ \mathcal{L}^{S} \f$ connect, a node with no such line forming
 * a single-node component (see
 * DCNetworkData::identify_connected_components()). In particular, the network
 * need not be connected. One reference node \f$ r_k \f$ is taken in each
 * component, and \f$ R = \{ r_1 , \ldots , r_{n_c} \} \f$. If \f$ n_c = 1 \f$
 * the reference is "ReferenceNode"; otherwise \f$ r_k \f$ is the node of
 * lowest index of \f$ \mathcal{N}_k \f$, except in the KIRCHHOFF formulation,
 * where the component that holds "ReferenceNode" takes it as its reference.
 * Note that the choice changes neither the flows nor the node injections that
 * the model admits, although it may change its numerical behavior.
 *
 * \par The PTDF matrix
 * Let \f$ \hat B \f$ be the \f$ |\mathcal{L}| \times N \f$ matrix with
 * \f$ \hat B_{l,s(l)} = \mathfrak{S}_l \f$ and \f$ \hat B_{l,e(l)} =
 * -\mathfrak{S}_l \f$ for \f$ l \in \mathcal{L}^{S} \f$ and 0 elsewhere.
 * The \f$ N \times N \f$ nodal susceptance matrix \f$ \bar B \f$ has
 * \f$ \bar B_{nn} = \sum_{ l \in \mathcal{L}^{S} : n \in \{ s(l) , e(l) \}
 * } \mathfrak{S}_l \f$ and \f$ \bar B_{nm} = - \sum_{ l \in
 * \mathcal{L}^{S} : \{ s(l) , e(l) \} = \{ n , m \} } \mathfrak{S}_l \f$
 * for \f$ n \neq m \f$ (parallel lines add up), i.e., \f$ \bar B =
 * \hat A^\top \hat B \f$ with \f$ \hat A \f$ the \f$ |\mathcal{L}| \times N
 * \f$ incidence matrix of the lines in \f$ \mathcal{L}^{S} \f$ (\f$ +1 \f$
 * at the start node, \f$ -1 \f$ at the end node, and a zero row for each
 * line in \f$ \mathcal{L}^{H} \f$). With \f$ I_R \f$ the
 * \f$ N \times ( N - n_c ) \f$ identity matrix without the columns of the
 * nodes in \f$ R \f$, the Power Transfer Distribution Factor (PTDF) matrix
 * is
 * \f[
 *   \Psi = \hat B I_R \bigl( I_R^\top \bar B I_R + \tau^T I \bigr)^{-1} ,
 *   \tag{P}
 * \f]
 * where \f$ \tau^T \ge 0 \f$ is the Tikhonov coefficient, 0 by default (see
 * generate_abstract_constraints()), and DCNetworkData::get_PTDF() computes it
 * with a sparse LU factorization. Let \f$ \tau^T = 0 \f$. For any vector
 * \f$ x \f$ of injections, the flows \f$ F = \Psi I_R^\top x \f$ satisfy
 * \f$ I_R^\top \hat A^\top F = I_R^\top x \f$, i.e., the balance of the lines
 * in \f$ \mathcal{L}^{S} \f$ at every node not in \f$ R \f$: indeed, since
 * \f$ \hat A^\top \hat B = \bar B \f$,
 * \f[
 *   I_R^\top \hat A^\top \Psi = I_R^\top \hat A^\top \hat B I_R
 *     ( I_R^\top \bar B I_R )^{-1} = ( I_R^\top \bar B I_R )
 *     ( I_R^\top \bar B I_R )^{-1} = I \; .
 * \f]
 * Each column of \f$ \hat A^\top \f$ (a line) has its two nonzeros in one
 * component, and therefore the entries of \f$ \hat A^\top F \f$ sum to zero
 * over each \f$ \mathcal{N}_k \f$; hence, the balance at \f$ r_k \f$, i.e.,
 * \f$ ( \hat A^\top F )_{r_k} = x_{r_k} \f$, holds if and only if
 * \f$ \sum_{ n \in \mathcal{N}_k } x_n = 0 \f$. In that case \f$ F \f$ are
 * the flows \f$ \hat B \theta \f$ given by the angles \f$ \theta \f$ with
 * \f$ \theta_{r_k} = 0 \f$ that solve \f$ \bar B \theta = x \f$, since
 * \f$ \theta = I_R ( I_R^\top \bar B I_R )^{-1} I_R^\top x \f$ satisfies the
 * rows of \f$ \bar B \theta = x \f$ of the nodes not in \f$ R \f$, and those
 * of \f$ R \f$ by the same sums. Since \f$ I_R^\top \bar B I_R \f$ is block
 * diagonal with one block per component (up to a permutation of the nodes),
 * \f$ \Psi \f$ is the juxtaposition of the PTDF matrices of the components,
 * each computed with the reference of its component; a component with a
 * single node has no column. The matrix \f$ I_R^\top \bar B I_R \f$ is
 * regular when each component is connected by lines with positive
 * susceptance. If it is singular (which may happen with susceptances of
 * opposite sign), the factorization fails and get_PTDF() throws
 * std::logic_error. With \f$ \tau^T > 0 \f$, instead, the rows that use
 * \f$ \Psi \f$ no longer imply the node balances exactly. As for the HVDC
 * lines, their distribution factors are \f$ \mathrm{DCDF} = - \Psi I_R^\top (
 * A^{H} )^\top \f$ (see DCNetworkData::compute_DCDF()), hence \f$ \Psi
 * I_R^\top ( a - ( A^{H} )^\top F^{H} ) = \Psi I_R^\top a + \mathrm{DCDF} \,
 * F^{H} \f$. This matrix depends on the instant through the efficiencies.
 *
 * \par Formulations
 * One can write the network in three formulations, which are chosen when the
 * variables are generated (see generate_abstract_variables()) and whose rows
 * are given in generate_abstract_constraints(). First, the KIRCHHOFF one (the
 * default) has the angle \f$ \theta_n \f$ of each node, the relation \f$ F_l
 * = \mathfrak{S}_l ( \theta_{s(l)} - \theta_{e(l)} ) \f$ for each
 * \f$ l \in \mathcal{L}^{S} \f$, the balance (B) at each node and one fixed
 * angle per component. Second, the PTDF one defines the flows of the lines in
 * \f$ \mathcal{L}^{S} \f$ through \f$ \Psi \f$ and DCDF, and completes them
 * with the balance of each component. Third, the CYCLE one writes the flows
 * of the lines in \f$ \mathcal{L}^{S} \f$ on a spanning forest and on the
 * fundamental cycles of the lines out of it, with the voltage law on each
 * cycle. The PTDF (with explicit flows) and CYCLE formulations are those
 * of
 *
 *   J. H&ouml;rsch, H. Ronellenfitsch, D. Witthaut and T. Brown, "Linear
 *   optimal power flow using cycle flows", Electric Power Systems Research
 *   158, 2018 (preprint arXiv:1704.01881).
 *
 * All three admit the same flows and node injections (when
 * \f$ \tau^T = 0 \f$ and the coefficients are not rounded), and therefore
 * give the same optimal value. In all of them the flows are explicit
 * variables, on which the bounds and the cost of the flows are written; the
 * classical PTDF form, in which the flows of the lines in
 * \f$ \mathcal{L}^{S} \f$ are eliminated and the bounds read \f$ P^{mn}_l \le
 * ( \Psi I_R^\top a )_l \le P^{mx}_l \f$ for a pure DC network, is obtained
 * by substitution.
 *
 * \par The net transfer capacity model
 * When \f$ \mathcal{L}^{S} = \emptyset \f$ the lines represent the commercial
 * exchanges between zones rather than physical flows. Then the network
 * reduces to the bounds on \f$ F^{H} \f$ and to the balance
 * \f$ a = ( A^{H} )^\top F^{H} \f$ at each node, i.e., to (B). Summing the
 * balances of all the nodes gives \f$ \sum_{ n \in \mathcal{N} } a_n = \sum_{
 * l \in \mathcal{L}^{H} } ( 1 - \sum_j \eta_{l,j} ) F_l \f$, which therefore
 * needs no row in the PTDF and KIRCHHOFF formulations. Also, the network need
 * not be connected, and the demand of a node with no line is then met by the
 * units at the node.
 *
 * \par Reduced networks
 * The network may be the reduction of a larger transmission grid, in which
 * sets of buses have been aggregated into zones and the lines joining two
 * zones replaced by one or more equivalent lines. The Block makes no
 * distinction between the two cases. Each node is then a zone, whose
 * injection \f$ S_n \f$ and demand \f$ D^{ac}_n \f$ are the sums of those of
 * the units and of the demands of the buses it contains; the data have to be
 * built in this way, since the Block has no information on the buses. Each
 * line carries the bounds "MinPowerFlow" and "MaxPowerFlow" and, if it obeys
 * the voltage law, the susceptance "LineSusceptance" given in the data, and
 * \f$ \Psi \f$ is always computed from these susceptances by (P). We do not
 * model a PTDF matrix estimated otherwise (e.g., fitted on the flows of the
 * original grid), a constant offset added to the flow of a line to account
 * for the flows that the aggregation hides, and the rules by which the buses
 * are clustered and the capacities of the equivalent lines are calibrated
 * (e.g., as the largest or as a percentile of the exchanges observed on the
 * original grid). Hence, a reduced network is represented exactly only when
 * each equivalent line has an equivalent susceptance; the capacities of the
 * equivalent lines are those given in the data. A line with zero susceptance
 * is a transport link between two zones, whose flow is bounded by its
 * capacities and appears only in the node balances.
 *
 * \par What is not modeled
 * Besides the reduction of a network above, we do not represent (i)
 * susceptances that change with the instant within one NetworkData (the flow
 * bounds and the efficiencies may change, through the dimension
 * "NumberInstants" of DCNetworkData::deserialize(), and a NetworkBlock with a
 * separate NetworkData may have different susceptances); (ii) the losses on
 * the lines in \f$ \mathcal{L}^{S} \f$, the reactive power and the voltage
 * magnitudes (see ACNetworkBlock); (iii) a hyperarc with nonzero
 * susceptance, i.e., the voltage law along a line with several end nodes,
 * which deserialize() rejects. */

class DCNetworkBlock : public NetworkBlock
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------------- PUBLIC TYPES --------------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
 *
 * DCNetworkBlock defines two public types:
 *
 * - formulation_type, the formulations of the network;
 *
 * - DCNetworkData, a small auxiliary class to bunch the basic electrical data
 *   of the transmission network.
 * @{ */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// the formulations of the network [see generate_abstract_constraints()]

 enum formulation_type
 {
  PTDF = 0 ,
  CYCLE ,
  KIRCHHOFF
  };

/*--------------------------------------------------------------------------*/
/*------------------- CLASS DCNetworkBlock::DCNetworkData ------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
 /// auxiliary class holding basic data about the (DC) transmission network
 /** The DCNetworkData class is a nested sub-class which only serves to have
  * a quick way to load all the basic data (topology and electrical
  * characteristics) that describe the transmission network. The rationale
  * is that, while often the network does not change during the (short) time
  * horizon of UC, it makes sense to allow for this to happen. This means
  * that individual NetworkBlock objects may in principle have different
  * DCNetworkData, but most often they can share the same. By bunching all
  * the information together, we make it easy for this sharing to happen. */

class DCNetworkData : public NetworkData
{
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

 public:

/** @} ---------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/** @name Constructor and Destructor
 * @{ */

 /// constructor of DCNetworkData, does nothing
 DCNetworkData( void ) : f_number_lines( 0 ) , f_number_HVDC_lines( 0 ) ,
  f_reference_node( 0 ) , DCDF_was_computed( Inf< Index >() ) ,
  f_number_branches( 0 ) , cycle_basis_was_computed( false ) ,
  nb_components( 0 ) {}

 /// destructor of DCNetworkData: it is virtual, and empty
 ~DCNetworkData() override = default;

/** @} --------------------- OTHER INITIALIZATIONS -------------------------*/
/** @name Other initializations
 * @{ */

 /// deserialize a DCNetworkData out of a netCDF::NcGroup
 /** Deserialize a DCNetworkData out of a netCDF::NcGroup, which should
  * contain the following:
  *
  * - The dimension "NumberNodes" containing the number of nodes in the
  *   problem; this dimension is optional, if it is not provided then it is
  *   taken to be equal to 1.
  *
  * If NumberNodes == 1 (equivalently, it is not provided), the network is
  * a "bus" formed of only one node, and therefore all the subsequent
  * information need not be present since it is not loaded. If
  * NumberNodes > 1, then all the subsequent information is considered:
  *
  * - The dimension "NumberLines" containing the number of lines in the
  *   transmission network. Each line can be either a "regular" line (one
  *   start node and one end node) or a hyperarc (one start node and
  *   several end nodes); see the (optional) dimension "NumberBranches"
  *   right next. The dimension is mandatory.
  *
  * - The dimension "NumberBranches" that is used to describe hyperarcs.
  *   The dimension is optional, if it is not defined then it is assumed that
  *   "NumberBranches" == "NumberLines", i.e., all lines are "regular".
  *   Otherwise, "NumberBranches" >= "NumberLines" (in fact, >) must hold
  *   since one single hyperarc is described by its multiple "branches", as
  *   detailed in "HyperArcID".
  *
  * - The dimension "NumberInstants" that specifies how many different time
  *   instants the possibly time-varying variables ("MaxPowerFlow",
  *   "MinPowerFlow" and "Efficiency") cover; the DCNetworkData may then be
  *   shared by the NetworkBlock of different instants, each reading the
  *   entry of its own instant. The dimension is optional, if it is not
  *   present it is taken to be 1, i.e., all the possibly time-varying
  *   variables are in fact time-static.
  *
  * - The variable "StartLine", of type netCDF::NcUint and indexed over the
  *   dimension "NumberBranches" (if it is defined, otherwise "NumberLines");
  *   the l-th entry of the variable is the starting point of the line (a
  *   number in 0, ..., NumberNodes - 1). Note that lines are not oriented,
  *   but the flow of energy is; that is, a positive flow along line l means
  *   that energy is being taken away from StartLine[ l ] and delivered to
  *   EndLine[ l ] (see next), a negative flow means vice versa. The variable
  *   is mandatory.
  *
  * - The variable "EndLine", of type netCDF::NcUint and indexed over the
  *   dimension "NumberBranches" (if it is defined, otherwise "NumberLines");
  *   the l-th entry of the variable is the ending point of the line (a number
  *   in 0, ..., NumberNodes - 1; lines are not oriented, but see above).
  *   StartLine[ l ] == EndLine[ l ] (a self-loop) is not allowed, but
  *   multiple lines between the same pair of nodes are. The variable is
  *   mandatory.
  *
  * - The variable "HyperArcID", of type netCDF::NcUint and indexed over the
  *   dimension "NumberBranches". The variable is mandatory if
  *   "NumberBranches" exists and therefore "NumberBranches" > "NumberLines",
  *   and ignored otherwise. The variable is used to specify which entries of
  *   "StartLine" and "EndLine" are different "branches" that correspond to
  *   the same hyperarc. The entries of the variable are a number in 0, ...,
  *   NumberLines - 1: HyperArcID[ i ] == l means that StartLine[ i ] and
  *   EndLine[ i ] describe one of the "branches" of the (hyper)line(arc) l.
  *   If a (hyper)line(arc) l has only one branch, i.e., HyperArcID[ i ] == l
  *   happens precisely for one index i, then l is a "regular" line. Note
  *   that, FOR EACH l = 0, ..., NumberLines - 1, THERE MUST BE AT LEAST ONE
  *   INDEX i such that HyperArcID[ i ] == l. If HyperArcID[ i ] == l happens
  *   for more than one index i, then l is a hyperarc (line). It is required
  *   that StartLine[ i ] == StartLine[ j ] and EndLine[ i ] != EndLine[ j ]
  *   for all pairs ( i , j ), i != j, such that HyperArcID[ i ] ==
  *   HyperArcID[ j ], i.e., all the branches of a hyperarc have the same
  *   start node and different end nodes, otherwise exception is thrown.
  *
  * - The variable "MaxPowerFlow", of type netCDF::NcDouble and indexed in
  *   principle over both dimensions "NumberLines" and "NumberInstants",
  *   whose entry [ l ][ t ] is the maximum flow \f$ P^{mx}_l \f$ of line
  *   \f$ l \f$ at the instant \f$ t \f$ of the NetworkBlock that uses it
  *   (see NetworkBlock::set_time_instant()). The variable can also be
  *   indexed over "NumberLines" only (and it must be so if
  *   "NumberInstants" is not defined), in which case the maximum flow of
  *   each line is the same at every instant. If line \f$ l \f$ is a
  *   hyperarc, \f$ P^{mx}_l \f$ is still one number, the maximum amount of
  *   flow leaving the start node, although then some flow (not necessarily
  *   the same amount, see "Efficiency") reaches more than one end node. The
  *   variable is optional; if it is not provided, \f$ P^{mx}_l = 0 \f$ for
  *   every line. The bound actually imposed on \f$ F_l \f$ is
  *   \f$ \kappa_l C^v P^{mx}_l \f$ (see DCNetworkBlock::deserialize() and
  *   DCNetworkBlock::generate_abstract_constraints()); when the network is
  *   the reduction of a larger one, it is the capacity of the equivalent
  *   line, which the data have to give.
  *
  * - The variable "MinPowerFlow", of type netCDF::NcDouble and indexed as
  *   "MaxPowerFlow", whose entry [ l ][ t ] is the minimum flow
  *   \f$ P^{mn}_l \f$ of line \f$ l \f$ at instant \f$ t \f$ (typically,
  *   but not necessarily, a negative number, since a line can be used in
  *   both directions). The variable is optional; if it is not provided,
  *   \f$ P^{mn}_l = 0 \f$ for every line, i.e., the flows only go from the
  *   start node to the end node(s).
  *
  * - The variable "LineSusceptance", of type netCDF::NcDouble and indexed
  *   over the dimension "NumberLines", whose entry \f$ l \f$ is the
  *   susceptance \f$ \mathfrak{S}_l \f$ of line \f$ l \f$, the same at
  *   every instant: \f$ l \in \mathcal{L}^{S} \f$ if
  *   \f$ \mathfrak{S}_l \neq 0 \f$ and \f$ l \in \mathcal{L}^{H} \f$
  *   otherwise. The variable is optional; if it is not provided, or if all
  *   its entries are 0, every line is in \f$ \mathcal{L}^{H} \f$ and the
  *   network is the net transfer capacity model of the class notes, which
  *   need not be connected. A hyperarc (see "HyperArcID") must have zero
  *   susceptance, otherwise exception is thrown, while the lines with one
  *   end node may have any susceptance also in a network with hyperarcs
  *   (see the class notes).
  *
  * - The dimension "ReferenceNode", the node whose voltage angle is fixed
  *   to 0, with the role described in the class notes: it is the reference
  *   of the PTDF matrix when the lines with nonzero susceptance connect all
  *   the nodes, and it fixes the angle of its component in the KIRCHHOFF
  *   formulation, while in the PTDF and CYCLE formulations each component
  *   of a network with several of them takes as reference its node of
  *   lowest index. The choice does not change the model, but it may change
  *   its numerical behavior, and it is irrelevant for a pure HVDC network.
  *   The dimension is optional, 0 by default; a value not smaller than
  *   "NumberNodes" makes deserialize() throw.
  *
  * - The variable "NetworkCost", of type netCDF::NcDouble and indexed over
  *   the dimension "NumberLines", whose entry \f$ l \f$ is the cost
  *   \f$ c^{net}_l \f$ of one unit of flow on line \f$ l \f$ in either
  *   direction, i.e., the objective has the term \f$ c^{net}_l | F_l | \f$
  *   (see DCNetworkBlock::generate_objective()). If line \f$ l \f$ is a
  *   hyperarc, the cost is still one number, the cost of a unit of flow
  *   leaving the start node. The variable is optional, all costs being 0
  *   by default; a negative cost would make the problem unbounded (the
  *   variable \f$ V_l \ge | F_l | \f$ that carries it being bounded only
  *   from below), and deserialize() throws std::logic_error.
  *
  * - The variable "Efficiency", of type netCDF::NcDouble and indexed in
  *   principle over both dimensions "NumberBranches" (if it is defined,
  *   otherwise "NumberLines") and "NumberInstants", whose entry [ l ][ t ]
  *   is the efficiency \f$ \eta_l \f$ of line \f$ l \f$ at instant
  *   \f$ t \f$: if \f$ F_l \f$ leaves the start node, \f$ \eta_l F_l \f$
  *   reaches the end node. The variable can also be indexed over
  *   "NumberBranches" ("NumberLines") only, in which case the efficiency is
  *   the same at every instant. The efficiency applies to the lines in
  *   \f$ \mathcal{L}^{H} \f$ only, that of a line with nonzero susceptance
  *   being 1 whatever the data say. If \f$ l \f$ is a hyperarc, each branch
  *   has its own efficiency \f$ \eta_{l,j} \f$: say, a hyperarc with the
  *   branches 1 \f$ \to \f$ 2 and 1 \f$ \to \f$ 3 of efficiency 0.5 each
  *   delivers half of the flow leaving node 1 to node 2 and the other half
  *   to node 3. The efficiencies of the branches of a hyperarc need not sum
  *   to 1, nor is an efficiency required to be at most 1. The variable is
  *   optional, every efficiency being 1 by default.
  *
  * - The variable "LineName", of type netCDF::NcString() and indexed over
  *   the dimension "NumberLines". Its i-th entry, namely LineName[ i ],
  *   contains the name of the i-th transmission line. This variable is
  *   optional. */

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG
 /// extends NetworkData::expected_dims()

 std::vector< std::string > expected_dims( void ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// extends NetworkData::expected_vars()

 std::vector< std::string > expected_vars( void ) const override;

#endif

/** @} ------- METHODS FOR READING THE DATA OF THE DCNetworkData -----------*/
/** @name Reading the data of the DCNetworkData
 * @{ */

 /// returns the number of lines of the network
 /** Method for returning the total number of lines of the network. When
  * get_number_nodes() == 1 (the network is a bus), get_number_lines() == 0
  * (no self-loops are allowed, hence there is no line to be made with a
  * single node). */

 Index get_number_lines( void ) const { return( f_number_lines ); }

/*--------------------------------------------------------------------------*/
 /// returns the number of HVDC lines of the network
 /** Method for returning the number of HVDC lines of the network, i.e.,
  * those with 0 susceptance. Clearly, get_number_HVDC_lines() <=
  * get_number_lines(): when the two are equal this is a HVDC grid, when
  * get_number_HVDC_lines() == 0 this is a DC grid, in between this is a
  * DC-HVDC grid. */

 Index get_number_HVDC_lines( void ) const { return( f_number_HVDC_lines ); }

/*--------------------------------------------------------------------------*/
 /// returns true if this is a pure HVDC grid

 bool is_HVDC( void ) const { return( v_line_susceptance.empty() ); }

/*--------------------------------------------------------------------------*/
 /// returns true if this is a pure DC grid

 bool is_DC( void ) const { return( f_number_HVDC_lines == 0 ); }

/*--------------------------------------------------------------------------*/
 /// returns true if this is a mixed DC - HVDC grid

 bool is_DC_HVDC( void ) const {
  return( ( ! v_line_susceptance.empty() ) && ( f_number_HVDC_lines > 0 ) );
  }
 
/*--------------------------------------------------------------------------*/
 /// returns the reference node of the network
 /** Method for returning the reference node of the network. */

 Index get_reference_node( void ) const { return( f_reference_node ); }

/*--------------------------------------------------------------------------*/
 /// returns true if the network is a hypergraph
 /** Method for returning true if the network is a hypergraph, i.e., if it has
  * at least one line with multiple head buses. When is_hypergraph() == false
  * the network is a "regular graph" and therefore get_end_line() has to be
  * used, while if is_hypergraph() == true the network is a hypergraph and
  * therefore get_end_lines() has to be used. */

 bool is_hypergraph( void ) const {
  return( f_number_branches > f_number_lines );
  }

/*--------------------------------------------------------------------------*/
 /// returns true if \p line is an hyperarc (more than one head bus)

 bool is_hyperarc( Index line ) const {
  if( is_hypergraph() )
   return( v_end_lines[ line ].size() > 1 );
  return( false );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of start buses for all lines
 /** Method for returning the vector of starting point of each line. This
  * vector may have empty size (bus network) or the size of number of lines,
  * then there are two possible cases:
  *
  * - if get_number_nodes() == 1, this vector has empty size which means there
  *   is no line at network (bus network), and this vector is not needed;
  *
  * - if get_number_nodes() > 1, this vector have size of f_number_lines and
  *    get_start_line()[ l ] gives starting (tail) bus of line l. */

 const std::vector< Index > & get_start_line( void ) const {
  return( v_start_line );
  }

/*--------------------------------------------------------------------------*/
 /// returns the start bus of line \p line

 Index get_start_line( Index line ) const {
  return( v_start_line[ line ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of end buses for all lines
 /** Method for returning the vector of ending point of each line. This
  * vector may have empty size (bus network) or the size of number of lines,
  * then there are three possible cases:
  *
  * - if get_number_nodes() == 1, this vector has empty size which means there
  *   is no line at network (bus network), and this vector does not need to
  *   be defined.
  *
  * - if get_number_nodes() > 1 and get_number_hyperarcs() == 0, then the
  *   network is a "regular graph", this vector have size of f_number_lines,
  *   and  get_end_line()[ l ] gives ending (head) bus of line l.
  *
  * - if get_number_nodes() > 1 and get_number_hyperarcs() > 0, then the
  *   network is a hypergraph and get_end_lines() must be used to get the
  *   set of end buses of the lines; this vector is empty if no line has a
  *   nonzero susceptance, and otherwise get_end_line()[ l ] is the first
  *   end bus of line l, i.e., the end bus of a line with nonzero
  *   susceptance, which has only one.
  */

 const std::vector< Index > & get_end_line( void ) const {
  return( v_end_line );
  }

/*--------------------------------------------------------------------------*/
 /// returns the end (first, in the hypergraph case) bus of line \p line

 Index get_end_line( Index line ) const {
  if( is_hypergraph() )
   return( v_end_lines[ line ].front() );
  return( v_end_line[ line ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of (sets of) end lines
 /** Method for returning the vector of sets of ending point of each line.
  * This vector is empty if get_number_hyperarcs() == 0, i.e., the network is
  * "regular graph" (which is true in particular if get_number_nodes() == 1),
  * otherwise  get_end_lines()[ l ] is a (const) std::vector< Index >
  * containing the end buses / nodes of line l. Line l is a "regular arc" if
  * get_end_lines()[ l ].size() == 1, and an hyperarc if
  * get_end_lines()[ l ].size() > 1 (it cannot obviously be 0). The number of
  * lines l such that get_end_lines()[ l ].size() > 1 is equal to
  * get_number_hyperarcs(). Each std::vector< Index > is ordered in increasing
  * sense and without repeated elements. */

 const std::vector< std::vector< Index > > & get_end_lines( void ) const {
  return( v_end_lines );
  }

/*--------------------------------------------------------------------------*/
 /// returns minimum power flow of the given \p line at the given \p time
 /** This method returns the minimum power flow of the given \p line at
  * the given \p time.
  *
  * @return the minimum power flow of the given \p line. */

 double get_min_power_flow( Index line , Index time ) const {
  assert( line < f_number_lines );
  if( ! v_min_power_flow.num_elements() )
   return( 0 );
  if( v_min_power_flow.shape()[ 1 ] == 1 )
   time = 0;
  return( v_min_power_flow[ line ][ time ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns maximum power flow of the given \p line
 /** This method returns the maximum power flow of the given \p line.
  *
  * @return the maximum power flow of the given \p line. */

 double get_max_power_flow( Index line , Index time ) const {
  assert( line < f_number_lines );
  if( ! v_max_power_flow.num_elements() )
   return( 0 );
  if( v_max_power_flow.shape()[ 1 ] == 1 )
   time = 0;
  return( v_max_power_flow[ line ][ time ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the efficiency of the given \p line at the given \p time
 /** Returns the efficiency of the given \p line at the given \p time. If no
  * efficiencies are specified or \p line has a nonzero susceptance, 1 is
  * returned. If the network has hyperarcs it can be called only for a line
  * with nonzero susceptance (see get_line_efficiencies()). */

  double get_line_efficiency( Index line , Index time ) const {
  assert( line < f_number_lines );
  if( ( ! v_line_susceptance.empty() ) && ( v_line_susceptance[ line ] != 0 ) )
   return( 1 );
  if( is_hypergraph() )
   throw( std::logic_error( "DCNetworkData::get_line_efficiency: called "
                            "for a line with zero susceptance in a "
                            "hypergraph" ) );
  if( v_efficiency.empty() )
   return( 1 );

  if( v_efficiency.shape()[ 1 ] == 1 )   // time-independent data
   time = 0;                             // just ignore the time
  return( v_efficiency[ line ][ time ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the efficiencies for all heads of hyperline \p line at \p time

 const std::vector< double > & get_line_efficiencies(
					 Index line , Index time  ) const {
  if( ! is_hypergraph() )
   throw( std::logic_error(
                      "get_line_efficiencies() called but no hypergraph" ) );
  if( v_h_efficiency.size() == 1 )  // time-independent data
   time = 0;                        // just ignore the time
  return( v_h_efficiency[ time ][ line ] );
  }
/*--------------------------------------------------------------------------*/
 /// returns the DC lines
 /** This function returns the DC lines in the transmission network, i.e.,
  * the lines in \f$ \mathcal{L}^{S} \f$ of the class notes, those with
  * nonzero susceptance, in increasing order.
  * @return the lines with nonzero susceptance. */

 const Subset & get_DC_lines( void ) {
  if( v_DC_lines.empty() && ( f_number_lines > f_number_HVDC_lines ) ) {
    v_DC_lines.reserve( f_number_lines - f_number_HVDC_lines );
    for( Index line_id = 0 ; line_id < f_number_lines ; ++line_id )
     if( v_line_susceptance[ line_id ] != 0 )
      v_DC_lines.push_back( line_id );
   }

  return( v_DC_lines );
  }

/*--------------------------------------------------------------------------*/
 /// returns the HVDC lines
 /** This function returns the HVDC lines in the transmission network, i.e.,
  * the lines in \f$ \mathcal{L}^{H} \f$ of the class notes, those with
  * zero susceptance, in increasing order.
  * @return the lines with zero susceptance. */

 const Subset & get_HVDC_lines( void ) {
  if( v_HVDC_lines.empty() && ( f_number_HVDC_lines > 0 ) ) {
   if( ! v_line_susceptance.empty() ) {
    v_HVDC_lines.reserve( f_number_HVDC_lines );
    for( Index line_id = 0 ; line_id < f_number_lines ; ++line_id )
     if( v_line_susceptance[ line_id ] == 0 )
      v_HVDC_lines.push_back( line_id );
    }
   else {
    v_HVDC_lines.resize( f_number_lines );
    std::iota( v_HVDC_lines.begin() , v_HVDC_lines.end() , 0 );
    }
   }
  return( v_HVDC_lines );
  }

/*--------------------------------------------------------------------------*/
 /// returns vector of the susceptances
 /** Method for returning the vector of the susceptances
  * \f$ \mathfrak{S}_l \f$ of the lines. The vector is empty if the network
  * has no line with nonzero susceptance (a bus, or a pure HVDC network,
  * even if "LineSusceptance" is given with all its entries 0); otherwise it
  * has get_number_lines() entries, and entry l is the susceptance of line
  * l. */

 const std::vector< double > & get_line_susceptance( void ) const {
  return( v_line_susceptance );
  }

/*--------------------------------------------------------------------------*/
 /// the column of node \p idx in the PTDF matrix, -1 for a reference node
 /** Returns the index of the column of node \p idx in the PTDF matrix
  * \f$ \Psi \f$ of (P) in the class notes, i.e., the position of the node
  * among those not in \f$ R \f$, or -1 if \p idx is a reference node. It
  * can be called only after get_PTDF(). */

 int get_reducedIdx( int idx ) const;

/*--------------------------------------------------------------------------*/
 /// the node of column \p idx of the PTDF matrix
 /** The inverse of get_reducedIdx(): returns the node whose column in the
  * PTDF matrix is \p idx. It can be called only after get_PTDF(). */

 int get_originalIdx( int idx ) const;

/*--------------------------------------------------------------------------*/
 /// compute the distribution factors of the HVDC lines at instant \p time
 /** Computes the \f$ |\mathcal{L}| \times |\mathcal{L}| \f$ matrix
  * \f$ \mathrm{DCDF} = - \Psi I_R^\top ( A^{H} )^\top \f$ of the class
  * notes, whose only nonzero columns are those of the lines in
  * \p HVDC_lines, from the PTDF matrix \p PTDF_matrix given by get_PTDF()
  * and from the efficiencies of the instant \p time; entry
  * \f$ ( l , h ) \f$ is the change of the flow of line \f$ l \f$ due to
  * a unit of flow on the HVDC line \f$ h \f$. The result is read with
  * get_DCDF(). */

 void compute_DCDF( c_Subset & HVDC_lines , const SpMat & PTDF_matrix ,
		    Index time );

/*--------------------------------------------------------------------------*/
 /// the matrix computed by the last call to compute_DCDF()

 const SpMat & get_DCDF( void ) const { return( DCDF ); }

/*--------------------------------------------------------------------------*/
 /// true if get_DCDF() holds the matrix of instant \p t

 bool was_DCDF_computed( Index t ) const {
  return( DCDF_was_computed == t );
  }

/*--------------------------------------------------------------------------*/
 /// compute the components of the lines with nonzero susceptance
 /** Computes the components \f$ \mathcal{N}_1 , \ldots ,
  * \mathcal{N}_{n_c} \f$ of the class notes, i.e., the connected
  * components of the graph whose nodes are those of the network and whose
  * edges are the lines with nonzero susceptance; a node with no such line
  * is a component of its own. The components are numbered in the order of
  * their node of lowest index, and the nodes of each are listed in
  * increasing order (see get_subgraphs()). */

 void identify_connected_components( void );

/*--------------------------------------------------------------------------*/
 /// the nodes of each component, in increasing order
 /** Returns the vector whose entry k holds the nodes of the component
  * \f$ \mathcal{N}_{k+1} \f$, in increasing order, so that its first entry
  * is the node of lowest index of the component; it is empty until
  * identify_connected_components() is called. */

 std::vector< std::vector< Index > > & get_subgraphs( void ) {
  return v_nodes_in_component;
  }

/*--------------------------------------------------------------------------*/
 /// the number of components, 0 until they are computed

 size_t get_nb_connected_components( void ) const {
  return nb_components;
  }

/*--------------------------------------------------------------------------*/
 /// returns the PTDF matrix
 /** Returns the PTDF matrix \f$ \Psi \f$ of (P) in the class notes, of
  * size \f$ |\mathcal{L}| \times ( N - n_c ) \f$, built on the lines in
  * \p DC_lines (which should be get_DC_lines(), since the components and
  * the references are always those of all the lines with nonzero
  * susceptance), with the Tikhonov coefficient \p tikhonov_coeff
  * (\f$ \tau^T \f$, 0 by default) added to the diagonal of the reduced
  * nodal susceptance matrix before it is factorized. The column of node
  * \f$ n \f$ is get_reducedIdx( n ), and the rows of the lines not in
  * \p DC_lines are 0. The reduced matrix and its inverse are kept, and
  * the inverse is reused as long as the matrix does not change. If the
  * network has no line with nonzero susceptance the matrix has no column;
  * if the reduced matrix cannot be factorized (it is singular, which may
  * happen with susceptances of opposite sign) std::logic_error is
  * thrown. */

 SpMat get_PTDF( c_Subset & DC_lines , double tikhonov_coeff = 0 );

/*--------------------------------------------------------------------------*/
 /// returns the PTDF matrix of all the lines with nonzero susceptance

 SpMat get_PTDF( void ) { return( get_PTDF( get_DC_lines() ) ); }

/*--------------------------------------------------------------------------*/
 /// the reduced nodal susceptance matrix and its inverse, as last computed

 std::pair< SpMat , SpMat > get_stored_B2( void ) {
  return( std::make_pair( stored_B2 , stored_B2_inv ) );
  }

/*--------------------------------------------------------------------------*/
 /// store the reduced nodal susceptance matrix and its inverse

 void set_stored_B2( const SpMat & B2 , const SpMat & B2_inv ) {
  // in case sizes mismatch
  stored_B2.resize( B2.rows() , B2.cols() );
  stored_B2_inv.resize( B2_inv.rows() , B2_inv.cols() );

  stored_B2 = B2;
  stored_B2_inv = B2_inv;
  }

/*--------------------------------------------------------------------------*/
 /// true if compute_cycle_basis() has been called

 bool was_cycle_basis_computed( void ) const {
  return( cycle_basis_was_computed );
  }

/*--------------------------------------------------------------------------*/
 /// mark the cycle basis as computed

 void set_cycle_basis_computed( void ) { cycle_basis_was_computed = true; }

/*--------------------------------------------------------------------------*/
 /// compute a spanning forest and the fundamental cycles of the DC lines
 /** Computes a spanning forest of the graph of the lines with nonzero
  * susceptance (the set \f$ \mathcal{L}^{S} \f$ of the class notes), with
  * one tree for each component that has such a line, rooted at the node of
  * lowest index of the component, recording for each node its parent and
  * the line joining it to the parent; then the fundamental cycle of each
  * line of \f$ \mathcal{L}^{S} \f$ that is not in the forest, i.e., the
  * line itself followed by the path in the forest from its end node back
  * to its start node. Since the forest is made of lines rather than of
  * pairs of nodes, parallel lines are handled: all of them but one are out
  * of the forest, each closing a cycle of two lines (a line from a node to
  * itself would be a cycle of its own, but deserialize() rejects it).
  * get_cycle_basis() and get_spanning_parent() return the results in terms
  * of nodes, get_lines_in_spanning_tree() and get_lines_in_cycles() in
  * terms of lines; these are the data of the CYCLE formulation of
  * DCNetworkBlock::generate_abstract_constraints(). */

 void compute_cycle_basis( void );

/*--------------------------------------------------------------------------*/
 /// the fundamental cycles as sequences of nodes
 /** Entry c is the sequence of the nodes of the fundamental cycle c: the
  * start node of the line out of the forest that closes it, its end node,
  * and then the nodes of the path in the forest back to the start node. */

 const std::vector< std::vector< Index > > & get_cycle_basis( void ) {
  if( ! cycle_basis_was_computed )
   this->compute_cycle_basis();
  return( v_cycle_basis );
  }

/*--------------------------------------------------------------------------*/
 /// the parent of each node in the spanning forest
 /** Entry n is the parent of node n in the spanning forest, n itself for a
  * root and -1 for a node with no line of nonzero susceptance. */

 const std::vector<int> & get_spanning_parent( void ) {
  if( ! cycle_basis_was_computed )
   this->compute_cycle_basis();
  return( m_spanning_parent );
  }

/*--------------------------------------------------------------------------*/
/// returns the DC lines of the spanning forest, with their orientation
/** Returns a map from each line of the spanning forest computed by
 * compute_cycle_basis() to +1 if its direction (from its start node to its
 * end node) goes from the parent to the child in the forest, and -1
 * otherwise, i.e., to \f$ - \epsilon_l \f$ in the flow rows of the CYCLE
 * formulation (see DCNetworkBlock::generate_abstract_constraints()). A
 * line parallel to one of the forest is not in the map, since it closes a
 * cycle (see get_lines_in_cycles()). */

 std::map< Index , int > get_lines_in_spanning_tree( void ) {
 if( ! cycle_basis_was_computed )
  compute_cycle_basis();

 const auto & start_line = get_start_line();

 // the tree edge of each non-root node is the very line the spanning tree
 // has been built with, so that a line parallel to it is not one
 std::map< Index , int > lines_in_spanning_tree;
 for( Index v = 0 ; v < m_spanning_line.size() ; ++v )
  if( m_spanning_line[ v ] >= 0 ) {
   const Index line_id = m_spanning_line[ v ];
   lines_in_spanning_tree[ line_id ] =
    ( start_line[ line_id ] == Index( m_spanning_parent[ v ] ) ) ? +1 : -1;
   }

 return( lines_in_spanning_tree );
 }

/*--------------------------------------------------------------------------*/
/// returns the fundamental cycles as DC lines, with their orientation
/** Returns, for each fundamental cycle \f$ c \f$ computed by
 * compute_cycle_basis(), a map from each of its lines \f$ l \f$ to
 * \f$ C_{lc} = +1 \f$ if the cycle traverses the line in its direction
 * (from its start node to its end node) and \f$ C_{lc} = -1 \f$ otherwise:
 * these are the coefficients of the cycle flows \f$ h_c \f$ in the flow
 * rows of the CYCLE formulation and of the flows in its voltage law (see
 * DCNetworkBlock::generate_abstract_constraints()). */

 std::vector< std::map< Index , int > > get_lines_in_cycles( void ) {
 if( ! cycle_basis_was_computed )
  compute_cycle_basis();

 return( v_line_cycles );
 }

/*--------------------------------------------------------------------------*/
 /// returns vector of the network cost
 /** Returns the vector of the costs \f$ c^{net}_l \f$ of a unit of flow
  * on each line ("NetworkCost"), which is empty if the costs are not given
  * (all 0) and has get_number_lines() entries otherwise. */

 std::vector< double > & get_network_cost( void ) {
  return( v_network_cost );
  }


/*--------------------------------------------------------------------------*/
 /// returns the vector containing the name of the lines

 const std::vector< std::string > & get_line_names( void ) const {
  return( v_line_names );
  }

/** @} --------------- METHODS FOR SAVING THE DCNetworkData ----------------*/
/** @name Methods for loading, printing & saving the DCNetworkData
 * @{ */

 /// serialize a DCNetworkData out of a netCDF::NcGroup
 /** Serialize a DCNetworkData out of a netCDF::NcGroup to the specific
  * format of a DCNetworkData. See
  * DCNetworkData::deserialize( netCDF::NcGroup ) for details of the format
  * of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/

 protected:

/*--------------------- PROTECTED METHODS OF THE CLASS ---------------------*/

/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/

 Index f_number_lines;           ///< number of lines of the network
 Index f_number_HVDC_lines;      ///< number of HVDC lines of the network

 /// the node "ReferenceNode" [see deserialize()]
 Index f_reference_node;

 /// the instant of the DCDF matrix computed last, to avoid recomputing it
 Index DCDF_was_computed;

 /// A boolean to avoid recomputing the cycle basis algorithm
 bool cycle_basis_was_computed;

 Index f_number_branches;     ///< the number of branches of all hyperarcs

 Subset v_start_line;         ///< vector of starting lines

 Subset v_end_line;           ///< vector of ending lines

 /// vector of (vector of) sets of ending lines for hyperarcs
 std::vector< Subset > v_end_lines;

 /// vector to store the susceptance of each line of the network
 std::vector< double > v_line_susceptance;

 /// matrix to store the minimum power flow at each line and time
 boost::multi_array< double , 2 > v_min_power_flow;

 /// matrix to store the maximum power flow at each line and time
 boost::multi_array< double , 2 > v_max_power_flow;

 /// vector to store the network cost at each line
 std::vector< double > v_network_cost;

 /** matrix to store the network efficiency of each line and time in the
  * graph case, effective only for HVDC lines and ignored otherwise;
  * v_efficiency[ l ][ t ] contains the efficiency of the head node of
  * line l at time t */
 boost::multi_array< double , 2 > v_efficiency;

 /** matrix of the network efficiency of each (hyper)line in the
  * hypergraph case, effective only for HVDC lines and ignored otherwise;
  * v_h_efficiency[ t ][ l ] is the std::vector< double > containing the
  * efficiency of all head nodes of (hyper)line l at time t */
 std::vector< std::vector< std::vector< double > > > v_h_efficiency;

 std::vector< std::string > v_line_names;  ///< Line names

 Subset v_DC_lines;      ///< the indices of DC lines (susceptance != 0)
 Subset v_HVDC_lines;    ///< the indices of HVDC lines (susceptance == 0)

 /// vector to store the cycle basis
 std::vector< Subset > v_cycle_basis;

 /// vector to store the spanning tree
 // -1 will be a default value for not belonging to the tree
 std::vector< int > m_spanning_parent;

 /// the DC line joining each node to its parent in the spanning tree
 // -1 for the roots and the nodes with no DC line
 std::vector< int > m_spanning_line;

 /// the fundamental cycles as DC lines, each with its direction
 std::vector< std::map< Index , int > > v_line_cycles;

 /// the component of each node [see identify_connected_components()]
 std::vector< int > v_component;
 size_t nb_components;  ///< the number of components, 0 until computed
 /// the nodes of each component, in increasing order
 std::vector< std::vector< Index > > v_nodes_in_component;
 std::vector< int > v_reduced_idx; // For each node the reduced index
 std::vector< int > v_original_idx; // For each reduced index the original node

 /// to not recompute each time the PTDF
 SpMat stored_B2;
 SpMat stored_B2_inv;

 /// the distribution factors of the HVDC lines [see compute_DCDF()]
 /** The product of minus the PTDF matrix and the transpose of the
  * incidence matrix of the HVDC lines with their efficiencies, restricted
  * to the nodes that are not references. */
 SpMat DCDF;

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

/*-------------------- PRIVATE FIELDS OF THE CLASS -------------------------*/

 SMSpp_insert_in_factory_h;

/*---------------------- PRIVATE METHODS OF THE CLASS ----------------------*/

 };  // end( class( DCNetworkData ) )

/** @} ---------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 * @{ */

 /// constructor of DCNetworkBlock
 /** Constructor of DCNetworkBlock, taking possibly a pointer of its
  * father Block. */

 explicit DCNetworkBlock( Block * f_block = nullptr )
  : NetworkBlock( f_block ) , f_NetworkData( nullptr ) , ftype( PTDF ) ,
    v_design( nullptr ) , f_C_v_scal( 1 ) , f_tikhonov_coeff( 0 ) ,
    f_ptdf_round( 1e-16 ) {}

/*--------------------------------------------------------------------------*/
 /// destructor of DCNetworkBlock

 virtual ~DCNetworkBlock() override;

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 * @{ */

 /// deserialize a DCNetworkBlock out of a netCDF::NcGroup
 /** Deserialize a DCNetworkBlock out of a netCDF::NcGroup, which should
  * contain all the data necessary to describe a NetworkBlock (see
  * NetworkBlock::deserialize(), which reads the constant "ConstantTerm" of
  * the objective) and possibly the following:
  *
  * - The dimensions and variables of a DCNetworkData (see
  *   DCNetworkData::deserialize()), which are read if the dimension
  *   "NumberNodes" is there; otherwise the DCNetworkData must have been
  *   passed by set_NetworkData() (typically by UCBlock, which shares one
  *   among all its NetworkBlock), and the number of nodes is read there.
  *
  * - The variable "ActiveDemand", of type netCDF::NcDouble and indexed over
  *   the dimension "NumberNodes" (or, if the DCNetworkData is not in the
  *   NcGroup, over any dimension of size equal to the number of nodes),
  *   whose entry n is the active demand \f$ D^{ac}_n \f$ of node n. The
  *   variable is optional; if it is not found in the NcGroup, then it must
  *   be passed (either before or after the call to deserialize()) by
  *   set_ActiveDemand() or set_active_demand(). Since both groups of data
  *   are optional, the NcGroup can actually be empty, all the data being
  *   passed by the in-memory interface; UCBlock has provisions for the
  *   NcGroup describing the NetworkBlock to be absent altogether (see
  *   UCBlock::deserialize()).
  *
  * The constants \f$ \kappa_l \f$ that multiply the flow bounds of the lines
  * (see generate_abstract_constraints()) are all 1 after deserialize(),
  * which does not read them: they are changed only by set_kappa(). */

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG
 // extends NetworkBlock::expected_dims()
 /* not necessary since DCNetworkBlock does not have any new dims save those
  * of the DCNetworkData that are automatically taken into account.

 std::vector< std::string > expected_dims( void ) const override;
 */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// extends NetworkBlock::expected_vars()

 std::vector< std::string > expected_vars( void ) const override;

#endif

/*--------------------------------------------------------------------------*/
 /// loads the DCNetworkBlock instance from a stream
 /** Like load( std::istream & ), if there is any Solver attached to this
  * DCNetworkBlock then a NBModification (the "nuclear option") is issued. */

 void load( std::istream & input , char frmt = 0 ) override {
  throw( std::logic_error( "DCNetworkBlock::load() not implemented yet" ) );
  }

/*--------------------------------------------------------------------------*/
 /// generate the abstract variables of the DCNetworkBlock
 /** Generates the variables of the formulation chosen by the Configuration
  * (see below). If the network has more than one node, they are:
  *
  * - the node injections \f$ S_n \f$ of the class notes, one for each node
  *   (the group "s_network" of NetworkBlock::generate_abstract_variables()),
  *   free variables through which UCBlock links the network to the units;
  *
  * - the flows \f$ F_l \f$, one for each line, in all the formulations
  *   (the group "p_flow_network", see get_power_flow());
  *
  * - the auxiliary variables \f$ V_l \ge | F_l | \f$, one for each line, if
  *   some line has a nonzero cost "NetworkCost" (the group "aux_network",
  *   see get_auxiliary_variable() and generate_network_cost_constraints());
  *
  * - in the CYCLE formulation, the cycle flows \f$ h_c \f$, one for each
  *   fundamental cycle of DCNetworkData::compute_cycle_basis() (the group
  *   "cycle_flow_network", see get_cycle_flow());
  *
  * - in the KIRCHHOFF formulation, if some line has a nonzero susceptance,
  *   the voltage angles \f$ \theta_n \f$, one for each node (the group
  *   "voltage_angle").
  *
  * All of them are continuous and free, their bounds being rows or bounds
  * of generate_abstract_constraints(). The formulation is given by an int
  * \f$ w \f$, the f_value of a SimpleConfiguration< int > that is either
  * \p stvv or, if \p stvv is nullptr and f_BlockConfig is not nullptr,
  * f_BlockConfig->f_static_variables_Configuration; if neither is a
  * SimpleConfiguration< int >, \f$ w = 2 \f$. The formulation is CYCLE if
  * \f$ w = 1 \f$, KIRCHHOFF if \f$ w = 2 \f$ and PTDF otherwise
  * (\f$ w = 0 \f$ among the others). */

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the flows and, if some line has a cost, the variables V
 /** Generates the flows \f$ F_l \f$ of all the lines and, if some line has
  * a nonzero "NetworkCost", the auxiliary variables \f$ V_l \f$ (see
  * generate_abstract_variables()): the variables of every formulation. */

 void generate_PTDF_variables( void );

/*--------------------------------------------------------------------------*/
 /// generate the cycle flows of the CYCLE formulation

 void generate_CYCLE_variables( void );

/*--------------------------------------------------------------------------*/
 /// generate the voltage angles of the KIRCHHOFF formulation

 void generate_KIRCHHOFF_variables( void );

/*--------------------------------------------------------------------------*/
 /// generate abstract constraints of DCNetworkBlock
 /** Generates the rows of the formulation chosen in
  * generate_abstract_variables(), in the notation of the class notes; none
  * if the network has a single node. With \f$ C^v \f$ the scaling factor of
  * the Configuration (see below) and \f$ \kappa_l \f$ the constant of line
  * \f$ l \f$ (1 unless set_kappa() changes it), every formulation has the
  * flow bounds
  * \f[
  *   \kappa_l C^v P^{mn}_l \le F_l \le \kappa_l C^v P^{mx}_l
  *   \qquad l \in \mathcal{L}
  *   \tag{1}
  * \f]
  * (the group "Power_flow_limit", see get_power_flow_limit_constraints()),
  * the same for the lines in \f$ \mathcal{L}^{S} \f$ and in
  * \f$ \mathcal{L}^{H} \f$. Here \f$ P^{mn}_l \f$ and \f$ P^{mx}_l \f$
  * are those of the instant of the Block (see
  * NetworkBlock::set_time_instant() and DCNetworkData::deserialize()). A
  * line with a design variable
  * \f$ x_l \f$ (see set_design_variables() and DesignNetworkBlock) has
  * instead the rows
  * \f[
  *   F_l - \kappa_l C^v P^{mx}_l x_l \le 0 \tag{2}
  * \f]
  * ("Power_flow_limit_design") and, only if \f$ P^{mn}_l \neq 0 \f$,
  * \f[
  *   F_l - \kappa_l C^v P^{mn}_l x_l \ge 0 \tag{3}
  * \f]
  * ("Power_flow_limit_design_min"), the lower half of (1), i.e.,
  * \f$ F_l \ge 0 \f$, being a bound when \f$ P^{mn}_l = 0 \f$. If some line
  * has a nonzero cost \f$ c^{net}_l \f$ ("NetworkCost"), each line has the
  * rows
  * \f[
  *   V_l - F_l \ge 0 \; , \qquad V_l + F_l \ge 0
  *   \tag{4}
  * \f]
  * ("power_flow_relax_abs", see generate_network_cost_constraints()).
  * Thus \f$ V_l \ge | F_l | \f$, and the term \f$ c^{net}_l V_l \f$ of the
  * objective (see generate_objective()) is \f$ c^{net}_l | F_l | \f$ at an
  * optimum, since \f$ c^{net}_l \ge 0 \f$ (a negative cost is refused).
  *
  * \par PTDF formulation
  * If \f$ \mathcal{L}^{S} = \emptyset \f$, the node balance (B) of the class
  * notes is written at every node,
  * \f[
  *   - S_n + \sum_{ l : s(l) = n } F_l
  *         - \sum_{ l , j : e_j(l) = n } \eta_{l,j} F_l = - D^{ac}_n
  *   \qquad n \in \mathcal{N}
  *   \tag{5}
  * \f]
  * ("HVDC_power_flow_injection"), and nothing else. Otherwise, with
  * \f$ \Psi \f$ and \f$ \mathrm{DCDF} \f$ of the class notes, the flow of
  * each line in \f$ \mathcal{L}^{S} \f$ is defined by
  * \f[
  *   F_l = \sum_{ n \in \mathcal{N} \setminus R } \Psi_{ln} \bigl( S_n -
  *         D^{ac}_n \bigr) + \sum_{ h \in \mathcal{L}^{H} }
  *         \mathrm{DCDF}_{lh} F_h
  *   \qquad l \in \mathcal{L}^{S}
  *   \tag{6}
  * \f]
  * (one row for each line in \f$ \mathcal{L}^{S} \f$, in the order of
  * DCNetworkData::get_DC_lines(), the group "AC/HVDC_powerflow_def"; the
  * column of \f$ \Psi \f$ of node \f$ n \f$ is
  * DCNetworkData::get_reducedIdx(), and the coefficients are rounded to
  * the nearest multiple of the precision \f$ \delta^{\Psi} \f$ of the
  * Configuration), the overall balance
  * \f[
  *   \sum_{ n \in \mathcal{N} } S_n
  *   + \sum_{ h \in \mathcal{L}^{H} } \Bigl( \sum_j \eta_{h,j} - 1 \Bigr)
  *     F_h = \sum_{ n \in \mathcal{N} } D^{ac}_n
  *   \tag{7}
  * \f]
  * ("overall_balanced_const", see overall_balance_function()), which is
  * the sum of (B) over all the nodes (the second sum is the loss of the
  * HVDC lines), and the node balance (B) at each node of a set
  * \f$ \mathcal{N}^{b} \f$,
  * \f[
  *   - S_n + \sum_{ l : s(l) = n } F_l
  *         - \sum_{ l , j : e_j(l) = n } \eta_{l,j} F_l = - D^{ac}_n
  *   \qquad n \in \mathcal{N}^{b}
  *   \tag{8}
  * \f]
  * ("DCHVDC_power_flow_injection", see generate_HVDC_nodal_constraints()),
  * the sums running over all the lines incident to \f$ n \f$. The set
  * \f$ \mathcal{N}^{b} \f$ holds the start and end nodes of the lines in
  * \f$ \mathcal{L}^{H} \f$ and, when there are several components, the
  * references \f$ r_k \f$ that are not among them, except the one of
  * lowest index. By the property of \f$ \Psi \f$ recalled in the class
  * notes, with \f$ \tau^T = 0 \f$ the rows (6) imply (B) at every node
  * not in \f$ R \f$, whatever the flows of the HVDC lines are. Then (8)
  * gives (B) at every reference but at most one, whose balance is the
  * difference between (7) and the sum of the others. Thus (6)-(8) amount
  * to (B) at every node, i.e., to the balance of each component
  * \f$ \mathcal{N}_k \f$ net of the flows of the HVDC lines that leave or
  * reach it, together with the flows of the lines in \f$ \mathcal{L}^{S}
  * \f$ that the injections determine. With one component and no line in
  * \f$ \mathcal{L}^{H} \f$, (8) is empty and (7) is \f$ \sum_n S_n =
  * \sum_n D^{ac}_n \f$; substituting (6) in (1) then gives the classical
  * form \f$ \kappa_l C^v P^{mn}_l \le ( \Psi I_R^\top ( S - D^{ac} ) )_l
  * \le \kappa_l C^v P^{mx}_l \f$.
  *
  * \par CYCLE formulation
  * The spanning forest of the lines in \f$ \mathcal{L}^{S} \f$ of
  * DCNetworkData::compute_cycle_basis() has one tree in each component that
  * has such a line, rooted at its node of lowest index, and a fundamental
  * cycle \f$ c \f$ for each line out of the forest (a line parallel to one
  * of the forest closes a cycle of two lines), with \f$ C_{lc} = \pm 1 \f$
  * when the cycle traverses \f$ l \f$ along or against its direction and 0
  * when it does not traverse it (DCNetworkData::get_lines_in_cycles()). For
  * a line \f$ l \f$ of the forest let \f$ \mathcal{N}^{T}_l \f$ be the
  * nodes of the subtree below it, and \f$ \epsilon_l = -1 \f$ if \f$ l \f$
  * goes from the parent to the child, \f$ \epsilon_l = +1 \f$ otherwise;
  * for a line out of the forest let \f$ \mathcal{N}^{T}_l = \emptyset \f$.
  * The flow of each line in \f$ \mathcal{L}^{S} \f$ is
  * \f[
  *   F_l = \epsilon_l \Bigl( \sum_{ n \in \mathcal{N}^{T}_l }
  *         \bigl( S_n - D^{ac}_n \bigr)
  *       + \sum_{ h \in \mathcal{L}^{H} } \bigl( \sum_j \eta_{h,j}
  *         [ e_j(h) \in \mathcal{N}^{T}_l ] - [ s(h) \in \mathcal{N}^{T}_l ]
  *         \bigr) F_h \Bigr) + \sum_c C_{lc} h_c
  *   \qquad l \in \mathcal{L}^{S}
  *   \tag{9}
  * \f]
  * ("v_CYCLE_def_flow_const", in the order of
  * DCNetworkData::get_DC_lines()), where \f$ [ \cdot ] \f$ is 1 if the
  * condition holds and 0 otherwise. In words, the flow of a line of the
  * forest carries the net injection of the subtree below it, including
  * what the HVDC lines bring in and take out of the subtree, plus the cycle
  * flows. The voltage law on each cycle is
  * \f[
  *   \sum_{ l \in \mathcal{L}^{S} } \frac{ C_{lc} }{ \mathfrak{S}_l } F_l = 0
  *   \qquad \text{for each cycle } c
  *   \tag{10}
  * \f]
  * ("v_CYCLE_def_cycle_const"), and the network is completed by (7) and
  * (8), both written also when \f$ \mathcal{L}^{S} = \emptyset \f$, in
  * which case every node is a component of its own. The rows (9) imply
  * (B) at every node that is not the node of lowest index of its component
  * (these nodes play here the role of \f$ R \f$ in the definition of
  * \f$ \mathcal{N}^{b} \f$): for a node \f$ n \f$ that is not the root, let
  * \f$ l \f$ be the line of the forest that joins \f$ n \f$ to its parent;
  * \f$ \epsilon_l F_l \f$ is the flow from \f$ n \f$ to its parent, and
  * the row (9) of \f$ l \f$ times \f$ \epsilon_l \f$, minus the rows of
  * the lines of the forest that join \f$ n \f$ to its children, each times
  * its own \f$ \epsilon \f$, leaves on the right the net injection of
  * \f$ n \f$ itself (and the flows of the HVDC lines at \f$ n \f$), since
  * the subtrees of the children make up that of \f$ n \f$ without
  * \f$ n \f$; on the left, with the lines out of the forest incident to
  * \f$ n \f$ (whose rows (9) are their cycle flows), it leaves the net flow
  * out of \f$ n \f$, each cycle flow \f$ h_c \f$ entering and leaving
  * \f$ n \f$ along the cycle, so that it cancels. This difference is (B)
  * at \f$ n \f$. Also, (10) is the voltage
  * law \f$ F_l =
  * \mathfrak{S}_l ( \theta_{s(l)} - \theta_{e(l)} ) \f$ written without
  * the angles; hence (9), (10), (7) and (8) admit the same flows and
  * injections as the other formulations.
  *
  * \par KIRCHHOFF formulation
  * The angles of the nodes define the flows of the lines in
  * \f$ \mathcal{L}^{S} \f$,
  * \f[
  *   F_l - \mathfrak{S}_l \theta_{s(l)} + \mathfrak{S}_l \theta_{e(l)} = 0
  *   \qquad l \in \mathcal{L}^{S}
  *   \tag{11}
  * \f]
  * ("KIRCHHOFF_power_flow_def"), the node balance (B) is written at every
  * node,
  * \f[
  *   - S_n + \sum_{ l : s(l) = n } F_l
  *         - \sum_{ l , j : e_j(l) = n } \eta_{l,j} F_l = - D^{ac}_n
  *   \qquad n \in \mathcal{N}
  *   \tag{12}
  * \f]
  * ("KIRCHHOFF_node_balance", see generate_node_balance_constraints()),
  * and one angle per component is fixed,
  * \f[
  *   \theta_{r_k} = 0 \qquad k = 1 , \ldots , n_c
  *   \tag{13}
  * \f]
  * ("reference_angle", see generate_reference_angle_constraint()), where
  * \f$ r_k \f$ is "ReferenceNode" in its component and the node of lowest
  * index in each of the others. The angles are therefore unique (that of a
  * node with no line of nonzero susceptance is 0). A pure HVDC network has
  * no angle and no row (11) nor (13), and (12) is (5).
  *
  * \par Dual values
  * By the convention of RowConstraint, for a minimization the dual value of
  * a row is the coefficient of its left-hand side in the Lagrangian, i.e.,
  * minus the derivative of the optimal value with respect to its active
  * bound. The nodal (or locational marginal) price of node \f$ n \f$ at
  * the instant of the Block is the derivative of the optimal value with
  * respect to \f$ D^{ac}_n \f$, which is minus the dual value
  * \f$ y^{ac}_{t,n} \f$ of the row of UCBlock that links \f$ S_n \f$ to the
  * units (see \ref ucbm_dual_net), as the stationarity with respect to
  * \f$ S_n \f$ requires. In the KIRCHHOFF formulation, and in the PTDF
  * formulation of a pure HVDC network, the demand of node \f$ n \f$
  * appears only in the right-hand side \f$ - D^{ac}_n \f$ of its node
  * balance (12), resp. (5), whose dual value is therefore
  * \f$ - y^{ac}_{t,n} \f$. In the PTDF and CYCLE formulations of the other
  * networks the demand of a node appears in several rows ((6) or (9), (7)
  * and possibly (8)), and the price is a combination of their dual values.
  * In the PTDF formulation, let \f$ y^{ov} \f$ be the dual value of (7),
  * \f$ y^{b}_n \f$ that of (8) at \f$ n \in \mathcal{N}^{b} \f$,
  * \f$ \mu_l \f$ that of the bounds (1) of line \f$ l \f$ (see
  * get_dual_prices(); \f$ \mu_l \ge 0 \f$ if the upper bound is active,
  * \f$ \mu_l \le 0 \f$ if the lower one is), \f$ \zeta_l \f$ the sum of
  * the dual values of the two rows (4) of the line times the coefficient
  * of \f$ F_l \f$ in them (0 if the line has no cost) and \f$ \pi_l \f$ the
  * dual value of (6). The stationarity with respect to the flow \f$ F_l \f$
  * of a line in \f$ \mathcal{L}^{S} \f$, which appears in (6) with
  * coefficient \f$ -1 \f$, in (1), in (4) and in the rows (8) at its ends,
  * gives
  * \f[
  *   \pi_l = \mu_l + \zeta_l + \sum_{ m \in \mathcal{N}^{b} } \hat A_{lm}
  *     y^{b}_m \; ,
  * \f]
  * and that with respect to \f$ S_n \f$ gives the price
  * \f[
  *   - y^{ac}_{t,n} = - y^{ov} + [ n \in \mathcal{N}^{b} ] \, y^{b}_n
  *     - [ n \notin R ] \sum_{ l \in \mathcal{L}^{S} } \Psi_{ln} \pi_l
  *   \qquad n \in \mathcal{N} \; ,
  * \f]
  * where \f$ [ \cdot ] \f$ is 1 if the condition holds and 0 otherwise.
  * The price of a node is thus the price \f$ - y^{ov} \f$ common to all the
  * nodes, plus the dual value of its node balance if it has one, minus the
  * congestion component, which sums over the lines the effect of an
  * injection at the node on the flow of the line times the marginal value
  * \f$ \pi_l \f$ of the line. For a network with one component, no line in
  * \f$ \mathcal{L}^{H} \f$ and no cost of the flows there is no node
  * balance, \f$ \pi_l = \mu_l \f$, and the price is
  * \f$ - y^{ov} - \sum_l \Psi_{ln} \mu_l \f$, the price \f$ - y^{ov} \f$ at
  * the reference node. The same relations hold for the multipliers
  * \f$ - y^{ac}_{t,n} \f$ when the NetworkBlock is a subproblem of the
  * Lagrangian dual of UCBlock, whose objective then has the term
  * \f$ - \sum_n y^{ac}_{t,n} S_n \f$; this is why relaxing also the rows
  * (7), (8) and the flow bounds gives the same dual bound (see
  * \ref ucbm_dual_net). The dual values are values per instant, as all the
  * costs (see \ref ucblock_model).
  *
  * \par Configuration
  * The parameters of the rows are given by \p stcc or, if \p stcc is
  * nullptr and f_BlockConfig is not nullptr, by
  * f_BlockConfig->f_static_constraints_Configuration, which can be:
  *
  * - a SimpleConfiguration< double > giving \f$ C^v \f$;
  *
  * - a SimpleConfiguration< std::pair< double , double > > giving
  *   \f$ C^v \f$ and the Tikhonov coefficient \f$ \tau^T \f$ of (P) in the
  *   class notes;
  *
  * - a SimpleConfiguration< std::vector< double > > giving, in this order
  *   and as far as the vector goes, \f$ C^v \f$, \f$ \tau^T \f$ and the
  *   precision \f$ \delta^{\Psi} \f$ to which the coefficients of
  *   \f$ \Psi \f$ and
  *   of DCDF in (6) are rounded.
  *
  * The defaults, used for what the Configuration does not give, are
  * \f$ C^v = 1 \f$, \f$ \tau^T = 0 \f$ and \f$ \delta^{\Psi} = 10^{-16} \f$.
  * A
  * positive \f$ \tau^T \f$ perturbs \f$ \Psi \f$, and then (6)-(8) no longer
  * imply (B) at every node; \f$ C^v \f$ scales the bounds (1)-(3) only. */

 void generate_abstract_constraints( Configuration * stcc = nullptr )
  override;

/*--------------------------------------------------------------------------*/
 /// generate the rows (5)-(8) of the PTDF formulation
 /** Generates the rows of the PTDF formulation of
  * generate_abstract_constraints(), but the bounds and the rows of the cost
  * of the flows; \p stcc is not used. */

 void generate_PTDF_constraints( Configuration * stcc = nullptr );

/*--------------------------------------------------------------------------*/
 /// generate the rows (7)-(10) of the CYCLE formulation
 /** Generates the rows of the CYCLE formulation of
  * generate_abstract_constraints(), but the bounds and the rows of the cost
  * of the flows; \p stcc is not used. */

 void generate_CYCLE_constraints( Configuration * stcc = nullptr );

/*--------------------------------------------------------------------------*/
 /// the overall balance of the PTDF and CYCLE formulations
 /** Returns the LinearFunction of the left-hand side of (7) of
  * generate_abstract_constraints(), the sum of the node injections less the
  * losses of the HVDC lines, and writes in \p constant_term its right-hand
  * side, the total active demand. */

 LinearFunction * overall_balance_function( double & constant_term );

/*--------------------------------------------------------------------------*/
 /// generate the rows (11)-(13) of the KIRCHHOFF formulation
 /** Generates the rows of the KIRCHHOFF formulation of
  * generate_abstract_constraints(), but the bounds and the rows of the cost
  * of the flows; \p stcc is not used. */

 void generate_KIRCHHOFF_constraints( Configuration * stcc = nullptr );

/*--------------------------------------------------------------------------*/
 /// generate the node balances (8) of the PTDF and CYCLE formulations
 /** Generates the node balance (8) of generate_abstract_constraints() at
  * each node of the set \f$ \mathcal{N}^{b} \f$ defined there (the ends
  * of the HVDC lines and the references of the components that the overall
  * balance (7) does not close), or at every node if \p full_formulation is
  * true; nothing is generated if the set is empty. */

 void generate_HVDC_nodal_constraints( bool full_formulation = false );

/*--------------------------------------------------------------------------*/
 /// generate the rows (4) that bound the auxiliary variables V
 /** Generates, for each line \f$ l \f$, the rows \f$ V_l - F_l \ge 0 \f$
  * and \f$ V_l + F_l \ge 0 \f$, i.e., (4) of
  * generate_abstract_constraints(), if some line has a nonzero
  * "NetworkCost", and nothing otherwise. It is called by
  * generate_abstract_constraints() in every formulation, and by the
  * derived classes that build their own rows. */

 void generate_network_cost_constraints( void );

/*--------------------------------------------------------------------------*/
 /// returns true if some line carries a non-zero network cost
 /** The auxiliary variables \f$ V_l \f$ and the rows (4) of
  * generate_abstract_constraints() exist only if this is true: an all-zero
  * "NetworkCost" vector states the same thing as a missing one. */

 bool has_network_cost( void ) const {
  if( ! f_NetworkData )
   return( false );
  const auto & nc = f_NetworkData->get_network_cost();
  return( std::any_of( nc.begin() , nc.end() ,
                       []( double cost ) { return( cost != 0. ); } ) );
  }

/*--------------------------------------------------------------------------*/
 /// generate the reference angles (13)
 /** Fixes to 0 one voltage angle in each component of the lines with
  * nonzero susceptance, i.e., (13) of generate_abstract_constraints():
  * that of "ReferenceNode" in its component and that of the node of lowest
  * index in each of the others. Nothing is generated for a pure HVDC
  * network, which has no angle. It is called by
  * generate_KIRCHHOFF_constraints() and by the derived classes that use
  * the angles. */

 void generate_reference_angle_constraint( void );

/*--------------------------------------------------------------------------*/
 /// generate the node balances (12) at every node
 /** Generates the node balance (B) of the class notes at every node, i.e.,
  * (12) of generate_abstract_constraints(), with all the lines incident to
  * the node, the efficiencies of the HVDC lines and the branches of the
  * hyperarcs. It is called by generate_KIRCHHOFF_constraints() and by the
  * derived classes that build their own rows. */

 void generate_node_balance_constraints( void );

/*--------------------------------------------------------------------------*/
 /// generate the flow bounds (1) and the design rows (2) and (3)
 /** Generates (1)-(3) of generate_abstract_constraints(). It is called by
  * generate_abstract_constraints() in every formulation, and by the
  * derived classes that build their own rows. */

 void generate_bound_constraints( void );

/*--------------------------------------------------------------------------*/
 /// rounds \p value to the nearest multiple of \p precision

 static double round_to( double value , double precision = 1.0 ) {
  return( std::round( value / precision ) * precision );
  }

/*--------------------------------------------------------------------------*/
 /// generate the objective of the DCNetworkBlock
 /** Generates the objective of the DCNetworkBlock, to be minimized,
  * \f[
  *   c^{0} + \sum_{ l \in \mathcal{L} } c^{net}_l V_l \; ,
  * \f]
  * where \f$ c^{0} \f$ is the constant "ConstantTerm" of
  * NetworkBlock::deserialize() and the sum is there only if some line has
  * a nonzero cost \f$ c^{net}_l \f$ ("NetworkCost"), in which case
  * \f$ c^{net}_l V_l = c^{net}_l | F_l | \f$ at an optimum by (4) of
  * generate_abstract_constraints() (the costs are nonnegative). The costs
  * are those of the data, in every formulation: no scale factor applies to
  * a NetworkBlock. \p objc is not used. */

 void generate_objective( Configuration * objc = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
/*---------------- Methods for checking the DCNetworkBlock -----------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for checking solution information in the DCNetworkBlock
 * @{ */

 /// returns true if the current solution is (approximately) feasible
 /** This function returns true if and only if the solution encoded in the
  * current value of the Variable of this DCNetworkBlock is approximately
  * feasible within the given tolerance. That is, a solution is considered
  * feasible if and only if
  *
  * -# each ColVariable is feasible; and
  *
  * -# the violation of each Constraint of this DCNetworkBlock is not
  *    greater than the tolerance.
  *
  * Every Constraint of this DCNetworkBlock is a RowConstraint and its
  * violation is given by either the relative (see RowConstraint::rel_viol())
  * or the absolute violation (see RowConstraint::abs_viol()), depending on
  * the Configuration that is provided.
  *
  * The tolerance and the type of violation can be provided by either \p fsbc
  * or f_BlockConfig->f_is_feasible_Configuration and they are determined as
  * follows:
  *
  * - If \p fsbc is not a nullptr and it is a pointer to a
  *   SimpleConfiguration< double >, then the tolerance is the value present
  *   in that SimpleConfiguration and the relative violation is considered.
  *
  * - If \p fsbc is not nullptr and it is a
  *   SimpleConfiguration< std::pair< double , int > >, then the tolerance is
  *   fsbc->f_value.first and the type of violation is determined by
  *   fsbc->f_value.second (any nonzero number for relative violation and
  *   zero for absolute violation);
  *
  * - Otherwise, if both f_BlockConfig and
  *   f_BlockConfig->f_is_feasible_Configuration are not nullptr and the
  *   latter is a pointer to either a SimpleConfiguration< double > or to a
  *   SimpleConfiguration< std::pair< double , int > >, then the values of the
  *   parameters are obtained analogously as above;
  *
  * - Otherwise, by default, the tolerance is 0 and the relative violation
  *   is considered.
  *
  * This function considers only the abstract representation to
  * determine if the solution is feasible. So, the parameter \p useabstract is
  * ignored. If no abstract Variable has been generated, then this
  * function returns true. Moreover, if no abstract Constraint has been
  * generated, the solution is considered to be feasible with respect to the
  * set of Variable only. Notice also that, before checking if the solution
  * satisfies a Constraint, the Constraint is computed
  * (Constraint::compute()).
  *
  * @param useabstract This parameter is ignored.
  *
  * @param fsbc The pointer to a Configuration that specifies the tolerance
  *             and the type of violation that must be considered. */

 bool is_feasible( bool useabstract = false ,
                   Configuration * fsbc = nullptr ) override;


/** @} ---------------------------------------------------------------------*/
/*---------- METHODS FOR READING THE DATA OF THE DCNetworkBlock ------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the NetworkBlock
 * @{ */

 /// returns the number of nodes
 /** Returns the number of nodes in the transmission network. If
  * get_NetworkData() does not return nullptr, this is
  * get_NetworkData()->get_number_nodes(). Otherwise, the network is a bus
  * and the method returns 1.
  *
  * @return the number of nodes in the network. */

 Index get_number_nodes( void ) const override {
  if( ! f_NetworkData )
   return( 1 );
  return( f_NetworkData->get_number_nodes() );
  }

/*--------------------------------------------------------------------------*/
 /// returns the number of lines of the network
 /** This function returns the number of lines in the transmission network.
  * If get_NetworkData() does not return nullptr, this is
  * get_NetworkData()->get_number_lines(). Otherwise, the network is a bus
  * and the method returns zero.
  *
  * @return the number of lines in the network. */

 Index get_number_lines( void ) const {
  if( ! f_NetworkData )
   return( 0 );
  return( f_NetworkData->get_number_lines() );
  }

/*--------------------------------------------------------------------------*/
 /// returns the susceptance of each lines of the network

 double get_line_susceptance( Index l ) const {
  if( ! ( f_NetworkData ) || f_NetworkData->get_line_susceptance().empty() )
   return( 0 );
  return( f_NetworkData->get_line_susceptance()[ l ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the kappa constant associated with the given \p line
 /** This function returns the constant \f$ \kappa_l \f$ of the given
  * \p line, which multiplies its minimum and maximum flow in (1)-(3) of
  * generate_abstract_constraints() (1 unless set_kappa() changes it).
  *
  * @param line The index of a line (between 0 and get_number_lines() - 1).
  *
  * @return The kappa constant associated with the given \p line. */

 double get_kappa( Index line ) const {
  if( v_kappa.empty() )
   return( 1 );
  assert( line < v_kappa.size() );
  return( v_kappa[ line ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the minimum power flow on the given \p line
 /** This function returns the minimum power flow on the given \p line. If
  * this DCNetworkBlock has no NetworkData, this function returns
  * 0. Otherwise, it returns the minimum power flow specified by the
  * NetworkData object.
  *
  * @param line The index of a line.
  *
  * @return The minimum power flow on the given \p line. */

 double get_min_power_flow( Index line ) const {
  if( ! f_NetworkData )
   return( 0 );
  return( f_NetworkData->get_min_power_flow( line , f_time_instant ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the maximum power flow on the given \p line
 /** This function returns the maximum power flow on the given \p line. If
  * this DCNetworkBlock has no NetworkData, thus function returns
  * 0. Otherwise, it returns the maximum power flow specified by the
  * NetworkData object.
  *
  * @param line The index of a line.
  *
  * @return The maximum power flow on the given \p line. */

 double get_max_power_flow( Index line ) const {
  if( ! f_NetworkData )
   return( 0 );
  return( f_NetworkData->get_max_power_flow( line , f_time_instant ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the factor the flow limits are scaled by
 /** The bounds on the power flow are stated as
  * \f$ \kappa C^{v} P^{mn} \leq F \leq \kappa C^{v} P^{mx} \f$, with
  * \f$ C^{v} \f$ this factor, 1 unless the Configuration of the static
  * Constraint says otherwise [see generate_abstract_constraints()]. Whoever
  * reads the duals of those bounds needs it: the derivative of a bound with
  * respect to the design is \f$ \kappa C^{v} P \f$ and not
  * \f$ \kappa P \f$. */

 double get_C_v_scal( void ) const { return( f_C_v_scal ); }

/*--------------------------------------------------------------------------*/
 /// returns the efficiency of the given \p line

 double get_line_efficiency( Index line ) const {
  if( ! f_NetworkData )
   return( 1 );
  return( f_NetworkData->get_line_efficiency( line , f_time_instant ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the efficiencies for all heads of hyperline \p line at \p time

 const std::vector< double > & get_line_efficiencies( Index line ) const {
  static const std::vector< double > _ret = { 1 };
  if( ! f_NetworkData )
   return( _ret );
  return( f_NetworkData->get_line_efficiencies( line , f_time_instant ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns a pointer to the DCNetworkData
 /** Return a pointer to the DCNetworkData. */

 NetworkData * get_NetworkData( void ) const override {
  return( f_NetworkData );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of active demands
 /** Returns the active demand for the given interval, which is assumed to
  * have size by get_number_nodes().
  *
  * @param interval The interval wrt the vector of demands for each user is
  *                 returned. */

 const double * get_active_demand( Index interval = 0 ) const override {
  if( v_ActiveDemand.empty() )
   return( nullptr );
  return( &( v_ActiveDemand.front() ) );
  }

/** @} ---------------------------------------------------------------------*/
/*---------- METHODS FOR READING THE Variable OF THE DCNetworkBlock --------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the Variable of the DCNetworkBlock
 * @{ */

 /// returns the vector of power flow variables
 /** The returned std::vector< ColVariable >, say F, contains the power flow
  * variables and is indexed over the dimension "NumberLines". There are two
  * possible cases:
  *
  * - if F is empty(), then this variable is not defined;
  *
  * - otherwise, F must have f_number_lines rows and F[ l ] is the power flow
  *   variable for line l. */

 const std::vector< ColVariable > & get_power_flow( void ) const {
  return( v_power_flow );
  }

/*--------------------------------------------------------------------------*/

 const std::vector< ColVariable > & get_cycle_flow( void ) const {
  return( v_cycle_flow );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of auxiliary variables
 /** The returned std::vector< ColVariable >, say V, contains the auxiliary
  * variables and is indexed over the dimension "NumberLines". There are two
  * possible cases:
  *
  * - if V is empty(), then this variable is not defined;
  *
  * - otherwise, V must have f_number_lines rows and V[ l ] is the auxiliary
  *   variable for line l. */

 const std::vector< ColVariable > & get_auxiliary_variable( void ) const {
  return( v_auxiliary_variable );
  }

/*--------------------------------------------------------------------------*/
 /// returns true if there is any design variable associated with some line

 bool has_design( void ) const {
  return( v_design && ( ! v_design->empty() ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns true if all lines have associated design variable

 bool all_design( void ) const {
  return( has_design() && ( v_design->size() == get_number_lines() ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the design variable associated with the given \p line
 /** This function returns a pointer to the design variable corresponding to
  * the physical transmission line indexed by \p line, or nullptr if there
  * is no such design variable.
  *
  * @param line The index of the line, between 0 and get_number_lines() - 1.
  *
  * @return A pointer to the design variable corresponding to \p line,
  *         or nullptr if that line has no design variable.
  *
  * \note The returned pointer refers to a variable owned externally by the
  *       corresponding DesignNetworkBlock; it must not be deleted or
  *       modified outside the intended modeling interface. */

 ColVariable * get_design( Index line ) const {
  if( ( ! v_design ) || v_design->empty() )
   return( nullptr );

  if( ! v_dense_design.empty() )
   return( v_dense_design[ line ] );

  if( line >= v_design->size() )
   return( nullptr );

  return( & ( *v_design )[ line ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the (const) design variable associated with the given \p line
 /** Const-qualified overload of get_design(). Returns a const pointer to the
  * design variable corresponding to the given physical line index \p line.
  * See get_design() for details. */

 const ColVariable * get_const_design( Index line ) const {
  return( get_design( line ) );
  }

/** @} ---------------------------------------------------------------------*/
/*--------- METHODS FOR READING THE Constraint OF THE DCNetworkBlock -------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the Constraint of the DCNetworkBlock
 * @{ */

 /// returns the vector of power flow limit constraints
 /** This function returns a const reference to the vector of power flow limit
  * constraints. The vector is empty if no line has the bound (1) of
  * generate_abstract_constraints(); otherwise element \f$ l \f$ is the
  * bound of line \f$ l \f$, empty (with no Variable) for a line with the
  * rows (2) and (3). */

 const std::vector< BoxConstraint > &
  get_power_flow_limit_constraints( void ) const {
  return( v_power_flow_limit_const );
  }

/*--------------------------------------------------------------------------*/
 /// returns the flow bounds, as get_power_flow_limit_constraints()

 const std::vector< BoxConstraint > &
  get_power_flow_limit_HVDC_bounds( void ) const {
  return( v_power_flow_limit_const );
  }

/*--------------------------------------------------------------------------*/
 /// return the vector of power losses on lines, all 0 in the DC approximation

 virtual std::vector< double > get_line_losses( void ) const {
  return( std::vector( get_number_lines() , 0. ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns the dual prices of power flow limits
 /** Writes in \p dp, for each line \f$ l \f$, the dual value \f$ \mu_l \f$
  * of its flow limits (1)-(3) of generate_abstract_constraints(), with the
  * sign convention of RowConstraint: \f$ \mu_l \ge 0 \f$ when the upper
  * limit is active, \f$ \mu_l \le 0 \f$ when the lower one is, and
  * \f$ - \mu_l \f$ is the derivative of the optimal value with respect to
  * the active limit. For a line with no design variable \f$ \mu_l \f$ is
  * the dual value of its bound (1); for a line with a design variable it
  * is the sum of the dual values of its rows (2) and (3), or of (2) and of
  * the lower half of (1) when the line has no row (3), the dual value of a
  * lower row being nonpositive as that of the lower side of a bound. The
  * vector is empty if the network has no line. */

 void get_dual_prices( std::vector< double > & dp ) const {
  auto nl = get_number_lines();
  if( ! nl ) {
   dp.clear();
   return;
   }

  dp.resize( nl );
  for( Index l = 0 ; l < nl ; ++l ) {
   if( ( l < v_design_row.size() ) &&
       ( v_design_row[ l ] < Inf< Index >() ) ) {  // design on this line
    // the upper side is always a row; the lower one is a row when the line
    // has a nonzero minimum flow, and the lower half of the bound otherwise.
    // The dual of the upper row is >= 0 and that of the lower one <= 0, as
    // the two sides of a BoxConstraint carry them in a single dual: they add
    dp[ l ] = v_power_flow_limit_design_const[ v_design_row[ l ] ].get_dual();

    if( v_design_min_row[ l ] < Inf< Index >() )
     dp[ l ] += v_power_flow_limit_design_min_const[ v_design_min_row[ l ]
                                                     ].get_dual();
    else
     if( ! v_power_flow_limit_const.empty() )
      dp[ l ] += v_power_flow_limit_const[ l ].get_dual();

    continue;
    }

   dp[ l ] = v_power_flow_limit_const.empty() ? 0 :
             v_power_flow_limit_const[ l ].get_dual();
   }
  }

/** @} ---------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Solution
 * @{ */

 /// returns a Solution representing the current solution of this NetworkBlock
 /** This method must construct and return a (pointer to a) Solution object
  * representing the current "solution state" of this NetworkBlock. This is
  * a DCNetworkBlockSolution extending NetworkBlockSolution with the specific
  * extra solution information of DCNetworkBlock.
  *
  * The parameter for deciding which kind of Solution must be returned is a
  * single int value, coded bitwise:
  *
  * - bit 0 (& 1) is "taken" by the base :NetworkBlock[Solution]
  *
  * - bit 1 (& 2) means "store the flow values"
  *
  * - bit 2 (& 4) means "store the dual prices"
  *
  * This value is to be found as:
  *
  * - if solc is not nullptr and it is a SimpleConfiguration< int >, then it
  *   is solc->f_value;
  *
  * - otherwise, if f_BlockConfig is not nullptr,
  *   f_BlockConfig->f_solution_Configuration is not nullptr and it is a
  *   SimpleConfiguration< int >, then it is
  *   f_BlockConfig->f_solution_Configuration->f_value;
  *
  * - otherwise, it is 7 (save everything).
  */

 Solution * get_Solution( Configuration * solc = nullptr ,
                          bool emptys = true ) override;

/*--------------------------------------------------------------------------*/
 /// return the "appropriate" [DC]NetworkBlockSolution

 NetworkBlockSolution * new_Solution( void ) const override;

/** @} ---------------------------------------------------------------------*/
/*--------------- METHODS FOR MODIFYING THE DCNetworkBlock -----------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the DCNetworkBlock
 * @{ */

 void set_NetworkData( NetworkData * nd = nullptr ) override {
  // if there was a previous DCNetworkData, and it was local, delete it
  if( f_NetworkData && f_local_NetworkData )
   delete( f_NetworkData );

  f_NetworkData = dynamic_cast< DCNetworkData * >( nd );
  f_local_NetworkData = false;
  }

/*--------------------------------------------------------------------------*/
 /// method to set the ActiveDemand
 /** The method is actually implemented since DCNetworkBlock is a concrete
  * class. Note that DCNetworkBlock always covers one interval only,
  * hence we expect v[] to contain just get_number_nodes() elements. */

 void set_ActiveDemand( const boost::multi_array< double , 2 > & v )
  override {
  if( v_ActiveDemand.empty() )
   v_ActiveDemand.assign( v[ 0 ].begin() , v[ 0 ].end() );
  }

/*--------------------------------------------------------------------------*/
 /// sets the power flows

 void set_power_flow( const std::vector< double > & pf ) {
  for( Index l = 0 ; l < v_power_flow.size() ; ++l )
   v_power_flow[ l ].set_value( pf[ l ] );
  }

/*--------------------------------------------------------------------------*/
 /// sets the dual prices of power flow limits
 /** Writes \p dp[ l ] as the dual value of the bound (1) of line \f$ l \f$
  * of generate_abstract_constraints(); nothing is written if no line has
  * such a bound (every line having the rows (2) and (3)), while the entry
  * of a line with the rows (2) and (3) is an empty BoxConstraint, whose
  * dual value is written as well but means nothing. */

 void set_dual_prices( const std::vector< double > & dp ) {
  auto nl = get_number_lines();
  if( ! nl )
   return;

  if( v_power_flow_limit_const.size() < nl )
   return;  // no bound to write the dual price on

  for( Index l = 0 ; l < nl ; ++l )
   v_power_flow_limit_const[ l ].set_dual( dp[ l ] );
  }

/*--------------------------------------------------------------------------*/
 /// set the (shared) design variables
 /** Sets the (shared) design variables, that are used to dimension all the
  * lines in the network. These are shared since they are typically decided
  * once and then used throughout all the (short-term) time horizon.
  *
  * Note that DCNetworkBlock retains the pointers to the two vectors \p DV
  * and \p Which, which therefore must not be changed by the caller for all
  * the lifetime of the object; in turn, DCNetworkBlock cannot change them.
  * (since they are const).
  *
  * If \p DV == nullptr or *DV.empty() then no line has design variables.
  * This is the default if this method is never called.
  *
  * If \p WDV == nullptr or WDV->empty(), it is assumed that DV[ i ]
  * refers to line i for all i = 0 , ... , DV.size() - 1. Otherwise,
  * DV[ i ] refers to line WDV[ i ]. If nonempty, \p WDV is supposed to
  * contain numbers in 0, ..., get_number_lines() - 1, be ordered in
  * increasing sense and without repeated elements.
  *
  * Must be called before DCNetworkBlock::generate_abstract_constraints(). */

 void set_design_variables( std::vector< ColVariable > * DV = nullptr ,
			    c_Subset * WDV = nullptr ) {
  if( constraints_generated() )
   throw( std::logic_error( "DCNetworkBlock::set_design_variables: called "
			    "when constraints are already generated" ) );
  v_design = DV;
  if( WDV && ( ! WDV->empty() ) ) {
   #ifndef NDEBUG
    for( Index i = 0 ; i < WDV->size() - 1 ; ++i )
     if( ( *WDV )[ i ] >= ( *WDV )[ i + 1 ] )
      throw( std::invalid_argument( "DCNetworkBlock::set_design_variables: "
				    "WDV not ordered" ) );
   #endif
   
   if( WDV->back() >= get_number_lines() )
    throw( std::invalid_argument( "DCNetworkBlock::set_design_variables: "
				  "invalid line number in WDV" ) );

   v_dense_design.resize( get_number_lines() , nullptr );
   for( Index i = 0 , j = 0 ; i < v_dense_design.size() ; ++i )
    if( ( j < WDV->size() ) && ( ( *WDV )[ j ] == i ) )
     v_dense_design[ i ] = & ( *v_design )[ j++ ];
   }
  else
   v_dense_design.clear();
  }

/** @} ---------------------------------------------------------------------*/
/*-------------------- METHODS FOR SAVING THE DCNetworkBlock ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the DCNetworkBlock
 * @{ */

 /// extends Block::serialize( netCDF::NcGroup )
 /** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
  * NetworkBlock. See NetworkBlock::deserialize( netCDF::NcGroup ) for
  * details of the format of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/** @} ---------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for changing the data of the DCNetworkBlock
 *  @{ */

 /// set the network cost values
 /** This function sets the network cost values of this DCNetworkBlock.
  *
  * @param values  Iterator to a vector containing the network cost values.
  * @param subset  If non-empty, the network cost values corresponding to the
  *                indices in \p subset are set to the values pointed by
  *                \p values. If empty, no operation is performed.
  * @param ordered It indicates whether \p subset is ordered.
  * @param issuePMod Controls how physical Modifications are issued.
  * @param issueAMod Controls how abstract Modifications are issued. */
 
 void set_network_cost( MF_dbl_it values ,
                        Subset && subset , bool ordered = false ,
                        c_ModParam issuePMod = eNoBlck ,
                        c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// set the network cost values
 /** This function sets the network cost values of this DCNetworkBlock.
  *
  * @param values Iterator to a vector containing the network cost values.
  * @param rng    If non-empty, the network cost values corresponding to the
  *               indices in \p rng are set to the values pointed by
  *               \p values. If empty, no operation is performed.
  * @param issuePMod Controls how physical Modifications are issued.
  * @param issueAMod Controls how abstract Modifications are issued. */

 void set_network_cost( MF_dbl_it values ,
                        Range rng = Range( 0 , Inf< Index >() ) ,
                        c_ModParam issuePMod = eNoBlck ,
                        c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// set the network cost values
 /** This function sets the network cost values of this DCNetworkBlock.
  *
  * @param value     The value of the network cost.
  * @param issuePMod Controls how physical Modifications are issued.
  * @param issueAMod Controls how abstract Modifications are issued. */

 void set_network_cost( double value , c_ModParam issuePMod = eNoBlck ,
                                       c_ModParam issueAMod = eNoBlck ) {
  std::vector< double > vector = { value };
  set_network_cost( vector.cbegin() , Range( 0 , Inf< Index >() ) ,
                    issuePMod , issueAMod );
  }

/*--------------------------------------------------------------------------*/
 /// set the kappa constants for the lines specified by \p subset
 /** This function sets the kappa constant of each line in the given \p
  * subset. The kappa constant of each line whose index is specified by the
  * i-th element in \p subset is given by the i-th element of the vector
  * pointed by \p values, i.e., it is given by the value pointed by (values +
  * i). The parameter \p ordered indicates whether the \p subset is ordered.
  *
  * The constant \f$ \kappa_l \f$ multiplies both flow limits of line
  * \f$ l \f$, \f$ \kappa_l C^v P^{mn}_l \le F_l \le \kappa_l C^v
  * P^{mx}_l \f$ in (1) of generate_abstract_constraints() (and the
  * coefficients of the design variable in (2) and (3)), and represents an
  * investment in the capacity of the line: with \f$ P^{mn}_l \f$ and
  * \f$ P^{mx}_l \f$ the limits of the largest line that can be built,
  * \f$ \kappa_l \in [ 0 , 1 ] \f$ is the fraction installed. For a line
  * with no design variable \f$ \kappa_l \f$ appears only in the
  * right-hand sides of (1), hence the optimal value of a convex problem
  * that contains the network (e.g., the continuous relaxation of a
  * UCBlock, or its Lagrangian dual) is a convex function of
  * \f$ \kappa_l \f$, and with \f$ \mu_l \f$ the dual value of the
  * limits of the line given by get_dual_prices() a subgradient of it is
  * \f[
  *   - C^v \mu_l P^{mx}_l \;\; \text{if } \mu_l \ge 0 \; , \qquad
  *   - C^v \mu_l P^{mn}_l \;\; \text{if } \mu_l \le 0 \; ;
  * \f]
  * summing it over the DCNetworkBlock of all the instants gives the
  * subgradient with respect to a capacity that is the same over the horizon
  * (see get_C_v_scal() and \ref ucbm_cap_sens). For a line with a design
  * variable \f$ x_l \f$, instead, \f$ \kappa_l \f$ multiplies \f$ x_l \f$ in
  * (2) and (3), a bilinear term, and the optimal value need not be convex in
  * \f$ \kappa_l \f$; its derivative where it exists is the one above times the
  * value of \f$ x_l \f$, but it is not a subgradient in general, and
  * InvestmentFunction refuses such a line. The set_kappa() of
  * ACNetworkBlock and of OTSNetworkBlock throws, since their lines have
  * limits in other rows as well.
  *
  * @param values An iterator to a vector containing the kappa constants.
  *
  * @param subset The indices of the lines whose kappa constants are being
  *        modified.
  *
  * @param ordered It indicates whether \p subset is ordered.
  *
  * @param issuePMod It controls how physical Modification are issued.
  *
  * @param issueAMod It controls how abstract Modification are issued. */

 virtual void set_kappa( MF_dbl_it values , Subset && subset ,
                         bool ordered = false ,
                         c_ModParam issuePMod = eNoBlck ,
                         c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// set the kappa constants for the lines specified by \p rng
 /** This function sets the kappa constants for the lines in the given Range
  * \p rng. For each i in the given Range (up to the number of lines minus 1),
  * the kappa constant for line i is given by the element of the vector
  * pointed by \p values whose index is (i - rng.first), i.e., it is given by
  * the value pointed by (values + i - rng.first).
  *
  * @param values An iterator to a vector containing the kappa constants.
  *
  * @param rng A Range containing the indices of the lines whose kappa
  *        constants are being modified.
  *
  * @param issuePMod It controls how physical Modification are issued.
  *
  * @param issueAMod It controls how abstract Modification are issued. */

 virtual void set_kappa( MF_dbl_it values ,
                         Range rng = Range( 0 , Inf< Index >() ) ,
                         c_ModParam issuePMod = eNoBlck ,
                         c_ModParam issueAMod = eNoBlck );

/*--------------------------------------------------------------------------*/
 /// change the abstract representation of the power flow limit constraints
 /** This function changes the abstract representation of the power flow limit
  * constraints for indices in \p modified_lines.
  *
  * @param modified_lines A vector of the indices of constraints that
  *                         have to be modified.
  *
  * @param issueAMod It controls how abstract Modification are issued. */

 void change_power_flow_limit_constraints( c_Subset & modified_lines,
					   c_ModParam issueAMod );

/*--------------------------------------------------------------------------*/
 /// change the abstract representation of the injection constraints
 /** This function changes the right-hand side of the node balance (5) of
  * generate_abstract_constraints() (PTDF formulation of a pure HVDC
  * network) of the nodes in \p modified_nodes; it does nothing if the
  * network has a line with nonzero susceptance.
  *
  * @param modified_nodes A vector of the indices of nodes that
  *                         have to be modified.
  *
  * @param issueAMod It controls how abstract Modification are issued. */

 void change_DC_power_flow_injection_constraints( c_Subset & modified_nodes ,
						  c_ModParam issueAMod );

/*--------------------------------------------------------------------------*/
 /// set the active demand at the nodes specified by \p subset
 /** This function sets the active demand at each node in the given \p
  * subset. The active demand at the node whose index is specified by the i-th
  * element in \p subset is given by the i-th element of the vector pointed by
  * \p values, i.e., it is given by the value pointed by (values + i). The
  * parameter \p ordered indicates whether the \p subset is ordered. If the
  * abstract representation has been generated, the rows that hold the
  * demand are changed by change_active_demand_constraints(), which a
  * derived class with rows of its own redefines.
  *
  * @param values An iterator to a vector containing the active demand.
  *
  * @param subset The indices of the nodes at which the active demand is being
  *        modified.
  *
  * @param ordered It indicates whether \p subset is ordered.
  *
  * @param issuePMod It controls how physical Modification are issued.
  *
  * @param issueAMod It controls how abstract Modification are issued. */

 void set_active_demand( MF_dbl_it values , Subset && subset ,
                         bool ordered = false ,
			 ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck ) override final;

/*--------------------------------------------------------------------------*/
 /// set the active demand at the nodes specified by \p rng
 /** This function sets the active demand at each node in the given Range \p
  * rng. For each i in the given Range (up to the number of nodes minus 1),
  * the active demand at node i is given by the element of the vector pointed
  * by \p values whose index is (i - rng.first), i.e., it is given by the
  * value pointed by (values + i - rng.first).
  *
  * @param values An iterator to a vector containing the active demand.
  *
  * @param rng A Range containing the indices of the nodes at which the active
  *        demand is being modified.
  *
  * @param issuePMod It controls how physical Modification are issued.
  *
  * @param issueAMod It controls how abstract Modification are issued. */

 void set_active_demand( MF_dbl_it values ,
                         Range rng = Range( 0 , Inf< Index >() ) ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck ) override final;

/*--------------------------------------------------------------------------*/
 /// change the abstract representation of the demand-dependent constraints
 /** This function updates the right-hand sides of all the rows of
  * generate_abstract_constraints() that depend on the active demand of the
  * nodes in \p modified_nodes, according to the formulation in #ftype: the
  * node balances (5) in the PTDF formulation of a pure HVDC network, the
  * rows (6), (7) and (8) in the PTDF formulation of the other networks, the
  * rows (9), (7) and (8) in the CYCLE formulation, and the node balances
  * (12) in the KIRCHHOFF formulation. Nothing is done if the network has a
  * single node, which has no such row. It is called by
  * set_active_demand(), and a derived class that builds rows of its own
  * that hold the demand (see ACNetworkBlock) redefines it.
  *
  * @param modified_nodes A vector containing the indices of the nodes whose
  *        active demand has been modified.
  *
  * @param issueAMod It controls how abstract Modification are issued. */

 virtual void change_active_demand_constraints( c_Subset & modified_nodes ,
                                                c_ModParam issueAMod );

/*--------------------------------------------------------------------------*/
 /// change the abstract representation of the PTDF constraints
 /** This function updates the abstract representation of the constraints
  * whose right-hand side depends on the active demand in the PTDF
  * formulation.
  *
  * In particular, it recomputes:
  *
  * - the constant terms of the power-flow definition constraints for all
  *   DC lines;
  *
  * - the constant term of the overall balance constraint;
  *
  * - the right-hand sides of the node balances (8) of
  *   generate_abstract_constraints(), if any.
  *
  * @param modified_nodes A vector containing the indices of the nodes whose
  *        active demand has been modified.
  *
  * @param issueAMod It controls how abstract Modification are issued. */

 void change_active_demand_constraints_PTDF( c_Subset & modified_nodes ,
                                             c_ModParam issueAMod );

/*--------------------------------------------------------------------------*/
 /// change the abstract representation of the CYCLE constraints
 /** This function updates the abstract representation of the constraints
  * whose right-hand side depends on the active demand in the CYCLE
  * formulation.
  *
  * In particular, it recomputes:
  *
  * - the constant terms of the flow-definition constraints associated with
  *   the spanning-tree part of the model;
  *
  * - the constant term of the overall balance constraint;
  *
  * - the right-hand sides of the node balances (8) of
  *   generate_abstract_constraints(), if any.
  *
  * @param modified_nodes A vector containing the indices of the nodes whose
  *        active demand has been modified.
  *
  * @param issueAMod It controls how abstract Modification are issued. */

 void change_active_demand_constraints_CYCLE( c_Subset & modified_nodes ,
                                              c_ModParam issueAMod );

/*--------------------------------------------------------------------------*/
 /// change the abstract representation of the KIRCHHOFF constraints
 /** This function updates the abstract representation of the constraints
  * whose right-hand side depends on the active demand in the KIRCHHOFF
  * formulation.
  *
  * In particular, it updates the node balances (12) of
  * generate_abstract_constraints() of the nodes listed in
  * \p modified_nodes, which are the only rows that hold the demand.
  *
  * @param modified_nodes A vector containing the indices of the nodes whose
  *        active demand has been modified.
  *
  * @param issueAMod It controls how abstract Modification are issued. */

 void change_active_demand_constraints_KIRCHHOFF( c_Subset & modified_nodes ,
						  c_ModParam issueAMod );

/*--------------------------------------------------------------------------*/
 /// change the right-hand sides of the node balances (8)
 /** This function updates the right-hand sides of the node balances (8) of
  * generate_abstract_constraints() of the nodes in \p modified_nodes that
  * have one, i.e., the rows generated by generate_HVDC_nodal_constraints()
  * and stored in #v_DC_HVDC_power_flow_const; it does nothing if there is
  * none (as in the KIRCHHOFF formulation).
  *
  * @param modified_nodes A vector containing the indices of the nodes whose
  *        active demand has been modified.
  *
  * @param issueAMod It controls how abstract Modification are issued. */

 void change_DC_HVDC_power_flow_injection_constraints(
                                        c_Subset & modified_nodes ,
                                        c_ModParam issueAMod );

/** @} ---------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED METHODS OF THE CLASS ---------------------*/
/*--------------------------------------------------------------------------*/

 DCNetworkData * get_new_NetworkData( void ) const override {
  return( new DCNetworkData() );
  }

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

/*---------------------------------- data ----------------------------------*/

 /// the DCNetworkData object
 DCNetworkData * f_NetworkData;

 /// vector to store the demand of each node of the network
 std::vector< double > v_ActiveDemand;

 std::vector< double > v_kappa;   ///< the kappa constant for each line

 formulation_type ftype;          ///< choice of model

 double f_C_v_scal;               ///< scaling factor C^v of the flow bounds

 double f_tikhonov_coeff;         ///< Tikhonov coefficient of the PTDF

 /// the precision to which the coefficients of the PTDF rows are rounded
 double f_ptdf_round;

/*-------------------------------- variables -------------------------------*/

 /// the power flow variables
 std::vector< ColVariable > v_power_flow;

 /// the power flow variables on cycle basis
 std::vector< ColVariable > v_cycle_flow;

 /// the voltage angle variables (Kirchhoff formulation)
 std::vector< ColVariable > v_voltage_angle;

 /// the auxiliary network cost variable
 std::vector< ColVariable > v_auxiliary_variable;

 /// the design variables for lines
 std::vector< ColVariable > * v_design;

 /// the "densified" version of v_design (if needed)
 std::vector< ColVariable * > v_dense_design;

/*------------------------------- constraints ------------------------------*/

 /// the node balances (5), PTDF formulation of a pure HVDC network
 /** One for each node; an ACNetworkBlock keeps here its own balances (see
  * ACNetworkBlock::generate_abstract_constraints()). */
 std::vector< FRowConstraint > v_power_flow_injection_const;

 /// the node balances (8) of the PTDF and CYCLE formulations
 /** One for each node of the set of generate_abstract_constraints(), in
  * increasing order of the node. */
 std::vector< FRowConstraint > v_DC_HVDC_power_flow_const;

 /// the rows (4), V_l - F_l >= 0 in [ 0 ] and V_l + F_l >= 0 in [ 1 ]
 boost::multi_array< FRowConstraint , 2 > v_power_flow_relax_abs;

 /// the rows (6) of the PTDF formulation
 /** One row for each DC line, in the order of DCNetworkData::get_DC_lines():
  * the HVDC lines have none. */
 std::vector< FRowConstraint > v_power_flow_def;

 /// Power flow limit constraints
 /** v_power_flow_limit_const[ l ] is the bound of line l, whose dual is the
  * price of its flow limit [see get_dual_prices()]; a line that has a design
  * variable and a lower row is entirely fenced by rows, and its entry is an
  * empty BoxConstraint, whose dual is 0. */
 std::vector< BoxConstraint > v_power_flow_limit_const;

 /// Upper power flow limit design constraints
 /** Only the lines that have a design variable have a row here, hence the
  * index is not the line: it is v_design_row[ line ]. */
 std::vector< FRowConstraint > v_power_flow_limit_design_const;

 /// Lower power flow limit design constraints
 /** Of the lines that have a design variable, only those whose minimum power
  * flow is not zero have a row here: with a zero minimum flow the design
  * variable has a zero coefficient in the row and what is left is the sign
  * of the flow, which reaches the solver as a bound (the lower half of
  * v_power_flow_limit_const) rather than as a row. The index is not the
  * line: it is v_design_min_row[ line ]. */
 std::vector< FRowConstraint > v_power_flow_limit_design_min_const;

 /// for each line, its row in v_power_flow_limit_design_const
 /** v_design_row[ line ] is the index of the line in
  * v_power_flow_limit_design_const, or Inf< Index >() if the line has no
  * design variable and therefore no row there. */
 std::vector< Index > v_design_row;

 /// for each line, its row in v_power_flow_limit_design_min_const
 /** v_design_min_row[ line ] is the index of the line in
  * v_power_flow_limit_design_min_const, or Inf< Index >() if the line has
  * no row there, be it because it has no design variable or because its
  * minimum power flow is zero. */
 std::vector< Index > v_design_min_row;

 /// the overall balance (7) of the PTDF and CYCLE formulations
 FRowConstraint overall_balanced_const;

 /// the rows (9) of the CYCLE formulation, in the order of the DC lines
 std::vector< FRowConstraint > v_CYCLE_def_flow_const;

 /// the rows (10) of the CYCLE formulation, one for each cycle
 std::vector< FRowConstraint > v_CYCLE_def_cycle_const;

 /// not used
 std::vector< FRowConstraint > v_CYCLE_def_HVDC_const;

 /// the rows (11) of the KIRCHHOFF formulation, in the order of the DC lines
 std::vector< FRowConstraint > v_KIRCHHOFF_power_flow_def;

 /// the node balances (12) of the KIRCHHOFF formulation, one for each node
 std::vector< FRowConstraint > v_KIRCHHOFF_node_balance_const;

 /// the reference angles (13) of the KIRCHHOFF formulation
 /** One for each component of the lines with nonzero susceptance, fixing
  * to 0 the angle of "ReferenceNode" in its component and that of the node
  * of lowest index in each of the others. */
 std::vector< BoxConstraint > v_reference_angle_const;

 /// the objective function
 FRealObjective objective;

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------- PRIVATE FIELDS OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/
/*---------------------- PRIVATE METHODS OF THE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

 static void static_initialization( void )
 {
  register_method< DCNetworkBlock , MF_dbl_it , Subset && , bool >(
   "DCNetworkBlock::set_network_cost" ,
   & DCNetworkBlock::set_network_cost );

  register_method< DCNetworkBlock , MF_dbl_it , Range >(
   "DCNetworkBlock::set_network_cost" ,
   & DCNetworkBlock::set_network_cost );

  register_method< DCNetworkBlock , MF_dbl_it , Subset && , bool >(
   "DCNetworkBlock::set_active_demand" ,
   & DCNetworkBlock::set_active_demand );

  register_method< DCNetworkBlock , MF_dbl_it , Range >(
   "DCNetworkBlock::set_active_demand" ,
   & DCNetworkBlock::set_active_demand );

  register_method< DCNetworkBlock , MF_dbl_it , Subset && , bool >(
   "DCNetworkBlock::set_kappa" ,
   & DCNetworkBlock::set_kappa );

  register_method< DCNetworkBlock , MF_dbl_it , Range >(
   "DCNetworkBlock::set_kappa" ,
   & DCNetworkBlock::set_kappa );
  }

/*--------------------------------------------------------------------------*/

 };  // end( class( DCNetworkBlock ) )

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS DCNetworkBlockMod -------------------------*/
/*--------------------------------------------------------------------------*/
/// derived class from NetworkBlockMod for modifications to a DCNetworkBlock

class DCNetworkBlockMod : public NetworkBlockMod
{
 public:

 /// public enum for the types of DCNetworkBlockMod
 enum DCNetB_mod_type
 {
  eSetKappa = eNetBModLastParam ,  ///< set the kappa constants
  eSetNetCost ,                    ///< set network cost values
  eDCNetBModLastParam  ///< first allowed parameter value for derived classes
  /**< Convenience value to easily allow derived classes to extend the set of
   * types of DCNetworkBlockMod. */
  };

 /// constructor, takes the DCNetworkBlock and the type
 DCNetworkBlockMod( DCNetworkBlock * const fblock , int type )
  : NetworkBlockMod( fblock , type ) {}

 /// destructor, does nothing
 virtual ~DCNetworkBlockMod() override = default;

 /// returns the Block to which the Modification refers
 Block * get_Block( void ) const override { return( f_Block ); }

 protected:

 /// prints the DCNetworkBlockMod
 void print( std::ostream & output ) const override {
  output << "DCNetworkBlockMod[" << this << "]: ";
  switch( f_type ) {
   case( eSetKappa ):
    output << "Set kappa values ";
    break;
   case( eSetNetCost ):
    output << "Set network cost values ";
    break;
   default:;
   }
  }
 };  // end( class( DCNetworkBlockMod ) )

/*--------------------------------------------------------------------------*/
/*---------------------- CLASS DCNetworkBlockRngdMod -----------------------*/
/*--------------------------------------------------------------------------*/
/// derived from DCNetworkBlockMod for "ranged" modifications

class DCNetworkBlockRngdMod : public DCNetworkBlockMod
{
 public:

 /// constructor: takes the DCNetworkBlock, the type, and the range
 DCNetworkBlockRngdMod( DCNetworkBlock * const fblock , int type ,
                        const Block::Range & rng )
  : DCNetworkBlockMod( fblock , type ) , f_rng( rng ) {}

 /// destructor, does nothing
 virtual ~DCNetworkBlockRngdMod() override = default;

 /// accessor to the range
 Block::c_Range & rng( void ) { return( f_rng ); }

 protected:

 /// prints the DCNetworkBlockRngdMod
 void print( std::ostream & output ) const override {
  DCNetworkBlockMod::print( output );
  output << "[ " << f_rng.first << ", " << f_rng.second << " )" << std::endl;
  }

 Block::Range f_rng;  ///< the range

 };  // end( class( DCNetworkBlockRngdMod ) )

/*--------------------------------------------------------------------------*/
/*----------------------- CLASS DCNetworkBlockSbstMod ----------------------*/
/*--------------------------------------------------------------------------*/
/// derived from DCNetworkBlockMod for "subset" modifications

class DCNetworkBlockSbstMod : public DCNetworkBlockMod
{
 public:

 /// constructor: takes the DCNetworkBlock, the type, and the subset
 DCNetworkBlockSbstMod( DCNetworkBlock * const fblock , int type ,
                        Block::Subset && nms )
  : DCNetworkBlockMod( fblock , type ) , f_nms( std::move( nms ) ) {}

 /// destructor, does nothing
 virtual ~DCNetworkBlockSbstMod() override = default;

 /// accessor to the subset
 Block::c_Subset & nms( void ) { return( f_nms ); }

 protected:

 /// prints the DCNetworkBlockSbstMod
 void print( std::ostream & output ) const override {
  DCNetworkBlockMod::print( output );
  output << "(# " << f_nms.size() << ")" << std::endl;
  }

 Block::Subset f_nms;  ///< the subset

 };  // end( class( DCNetworkBlockSbstMod ) )

/*--------------------------------------------------------------------------*/
/*-------------------- CLASS DCNetworkBlockSolution ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// a [NetworkBlock]Solution of a DCNetworkBlock
/** The DCNetworkBlockSolution class derives from NetworkBlockSolution and
 * adds the "standard" information stored in there (the node injection
 * variables) the other information that is typical of the DCNetworkBlock,
 * i.e.,
 *
 * - the flow variables on each link
 *
 * - [if available] the dual prices \f$ \mu_l \f$ of the flow limits of the
 *   lines given by DCNetworkBlock::get_dual_prices(), with their sign
 *   (nonnegative when the upper limit is active, nonpositive when the lower
 *   one is)
 *
 * Note that one DCNetworkBlock covers one time instant, so these variables
 * do not need to be indexed over time instants (unlike those of the base
 * NetworkBlockSolution). */

class DCNetworkBlockSolution : public NetworkBlockSolution
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*------------------------------- FRIENDS ----------------------------------*/

 friend DCNetworkBlock;  ///< make DCNetworkBlock friend

/*---------- CONSTRUCTING AND DESTRUCTING DCNetworkBlockSolution -----------*/

 /// constructor, does nothing
 explicit DCNetworkBlockSolution( void ) : NetworkBlockSolution() ,
  f_number_lines( 0 ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// deserialize a DCNetworkBlockSolution from a netCDF::NcGroup

 void deserialize( const netCDF::NcGroup & group ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// deserialize a DCNetworkBlockSolution from a "global" netCDF::NcGroup

 void deserialize( const netCDF::NcGroup & group , size_t idx ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 ~DCNetworkBlockSolution() override = default;
 ///< destructor: it is virtual, and empty

/*------ METHODS DESCRIBING THE BEHAVIOR OF A DCNetworkBlockSolution ------*/

 void read( const Block * block ) override;

 void write( Block * block ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// serialize a DCNetworkBlockSolution into a netCDF::NcGroup
 /** Serialize a DCNetworkBlockSolution into a netCDF::NcGroup. The format is
  * the one of NetworkBlockSolution, cf. the comments in
  * NetworkBlockSolution::serialize( netCDF::NcGroup & ), except that
  *
  *     "NumberNetworks" IS NOT REALLY NEEDED, BECAUSE "TotalNumberInstants"
  *     AND "EndInstant" ARE NOT REQUIRED SINCE DCNetworkBlock ALWAYS HAS
  *     DCNetworkBlock::get_number_intervals() == 1, AND ALL THE
  *     NetworkBlock IN \p group ARE SUPPOSED TO BE DCNetworkBlock
  *
  * In addition, \p group must contain:
  *
  * - The dimension "NumberLines" containing the number of lines in the
  *   transmission network. It is mandatory. Note that
  *
  *       ALL THE DCNetworkBlock MUST HAVE THE SAME NUMBER OF LINES
  *
  * - The variable "FlowValue", of type netCDF::NcDouble and indexed over
  *   the dimension "NumberLines"; FlowValue[ l ] is the optimal value of
  *   the power flow on line l. The variable is optional.
  *
  * - The variable "DualCost", of type netCDF::NcDouble and indexed over
  *   the dimension "NumberLines"; DualCost[ l ] is the dual price of the
  *   flow limits of line l given by DCNetworkBlock::get_dual_prices(). The
  *   variable is optional. */

 void serialize( netCDF::NcGroup & group ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// serialize a DCNetworkBlockSolution into a "global" netCDF::NcGroup
 /** "nonstandard" version of serialize() that loads a DCNetworkBlockSolution
  * from a "global" netCDF::NcGroup, i.e., one where the solution information
  * of multiple DCNetworkBlock are stored together (to avoid performance
  * issues due to the fact that netCDF is not structured to work with a large
  * number of sub-NcGroup in a file). The format is the  "nonstandard" one of
  * NetworkBlockSolution, cf. the comments in
  * NetworkBlockSolution::serialize( netCDF::NcGroup & , size_t ), except
  * that
  *
  *     "NumberNetworks" IS NOT REALLY NEEDED, BECAUSE "TotalNumberInstants"
  *     AND "EndInstant" ARE NOT REQUIRED SINCE DCNetworkBlock ALWAYS HAS
  *     DCNetworkBlock::get_number_intervals() == 1, AND ALL THE
  *     NetworkBlock IN \p group ARE SUPPOSED TO BE DCNetworkBlock
  *
  * In addition, \p group must contain:
  *
  * - The dimension "NumberLines" containing the number of lines in the
  *   transmission network. It is mandatory. Note that there is only one
  *   copy of the dimension, and as a consequence
  *
  *       ALL THE DCNetworkBlock MUST HAVE THE SAME NUMBER OF LINES
  *
  *   (which is of course necessary since they all take their data from
  *   the same variables where "NumberLines" is one of the dimensions)
  *
  * - The variable "FlowValue", of type netCDF::NcDouble and indexed both
  *   over the dimension "NumberNetworks" (which is the same as
  *   "TotalNumberInstants", that does not exist) and the dimension
  *   "NumberLines"; FlowValue[ idx ][ l ] is the optimal value of
  *   the power flow on line l for this DCNetworkBlock. The variable is
  *   optional.
  *
  * - The variable "DualCost", of type netCDF::NcDouble and indexed both
  *   over the dimension "NumberNetworks" (which is the same as
  *   "TotalNumberInstants", that does not exist) and the dimension
  *   "NumberLines"; DualCost[ idx ][ l ] is the dual price of the flow
  *   limits of line l of this DCNetworkBlock given by
  *   DCNetworkBlock::get_dual_prices(). The variable is optional.
  *
  * Note that the variables are constructed when \p idx == 0 according to
  * the fact that the corresponding DCNetworkBlockSolution has or not been
  * Configure-d to hold them, which means that
  *
  *       ALL THE DCNetworkBlockSolution MUST HAVE BEEN Configure-d IN THE
  *       SAME WAY
  *
  * (although, technically, if some of the DCNetworkBlockSolution that
  * appears when \p idx > 0 is Configure-d with less information than that
  * when idx == 0 the code will not break, but there will be uninitialised
  * values in the netCDF). */

 void serialize( netCDF::NcGroup & group , size_t idx ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 DCNetworkBlockSolution * scale( double factor ) const override;

 void sum( const Solution * solution , double multiplier ) override;

 DCNetworkBlockSolution * clone( bool empty = false ) const override;

/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/

 void print( std::ostream &output ) const override {
  output << "DCNetworkBlockSolution [" << this << "]: " << std::endl;
  }

/*-------------------------- PROTECTED FIELDS ------------------------------*/

 Index f_number_lines;          ///< the number of lines

 std::vector< double > v_flow;  ///< v_flow[ l ] = flow variable on line l

 std::vector< double > v_cost;  /**< v_cost[ l ] = dual price of the flow
                                 *   limits of line l [see
                                 *   DCNetworkBlock::get_dual_prices()] */

/*---------------------- PRIVATE PART OF THE CLASS -------------------------*/

 private:

/*---------------------------- PRIVATE FIELDS ------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class( DCNetworkBlockSolution ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 }  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* __DCNetworkBlock */

/*--------------------------------------------------------------------------*/
/*-------------------- End File DCNetworkBlock.h ---------------------------*/
/*--------------------------------------------------------------------------*/
