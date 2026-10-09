/*--------------------------------------------------------------------------*/
/*--------------------------- File ACNetworkBlock.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 *
 * Header file for class ACNetworkBlock, which derives from DCNetworkBlock and
 * defines the standard SOCP relaxation corresponding to the "AC model"
 * of the transmission network in the Unit Commitment problem.
 *
 * \author Quentin Jacquet \n
 *         EDF R&D OSIRIS \n
 *
 * \author Wim van Ackooij \n
 *         EDF R&D OSIRIS \n
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Quentin Jacquet, Wim van Ackooij,
 *            Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __ACNetworkBlock
 #define __ACNetworkBlock
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "DCNetworkBlock.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)

namespace SMSpp_di_unipi_it
{
/*--------------------------------------------------------------------------*/
/*------------------------- CLASS ACNetworkBlock ---------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// an "AC" transmission NetworkBlock (SOCP relaxation)
/** ACNetworkBlock derives from DCNetworkBlock and adds the numerous Variable
 * and Constraint necessary to represent the AC version of Kirchhoff's laws by
 * means of a Second-Order Cone Programming (SOCP) relaxation, possibly
 * strengthened with McCormick-like inequalities, i.e., the QC relaxation of
 * the Optimal Power Flow of
 *
 *   C. Coffrin, H. L. Hijazi and P. Van Hentenryck, "The QC Relaxation: A
 *   Theoretical and Computational Study on Optimal Power Flow", IEEE
 *   Transactions on Power Systems 31(4), 3008-3018, 2016,
 *   doi:10.1109/TPWRS.2015.2463111
 *
 * with the convex envelopes of the trigonometric terms of
 *
 *   H. L. Hijazi, C. Coffrin and P. Van Hentenryck, "Convex quadratic
 *   relaxations for mixed-integer nonlinear programs in power systems",
 *   Mathematical Programming Computation 9, 321-367, 2017,
 *   doi:10.1007/s12532-016-0112-z.
 *
 * In the notation of DCNetworkBlock, a line of an ACNetworkBlock is an HVDC
 * line if its susceptance, reactance and resistance are all 0, and an AC line
 * otherwise (see ACNetworkData::deserialize()); the AC lines are those that
 * DCNetworkData::get_DC_lines() returns. Each line \f$ l \f$ has an active
 * and a reactive flow at its start node, \f$ F^{fr}_l \f$ and
 * \f$ Q^{fr}_l \f$, and at its end node, \f$ F^{to}_l \f$ and
 * \f$ Q^{to}_l \f$, all leaving the node. Each node \f$ n \f$ has the active
 * and reactive injections \f$ S_n \f$ and \f$ R_n \f$ and the variable
 * \f$ w_n \f$ that stands for \f$ | V_n |^2 \f$. Each AC line has the
 * variables \f$ c_l \f$ and \f$ s_l \f$ that stand for the real and imaginary
 * parts of \f$ V_{s(l)} V_{e(l)}^* \f$, where \f$ V_n \f$ is the complex
 * voltage of node \f$ n \f$ (all in per unit and multiplied by the scaling
 * factors of generate_abstract_constraints()). The rows are the active and
 * reactive balances of each node, which hold the demands, and the definition
 * of the four flows of each AC line as linear functions of \f$ w \f$,
 * \f$ c \f$ and \f$ s \f$ through the admittances of the line. Then come the
 * thermal limits \f$ ( F^{fr}_l )^2 + ( Q^{fr}_l )^2 \le ( r^A_l )^2 \f$ and
 * \f$ ( F^{to}_l )^2 + ( Q^{to}_l )^2 \le ( r^A_l )^2 \f$ of each AC line,
 * the bounds on \f$ w_n \f$ given by the voltage limits, the bounds on the
 * angle difference, \f$ \tan( \phi^{mn}_l ) c_l \le s_l \le \tan( \phi^{mx}_l
 * ) c_l \f$, if the data give them (the angle differences are not bounded
 * otherwise, see ACNetworkData::deserialize()), and the cone
 * \f$ c_l^2 + s_l^2 \le w_{s(l)} w_{e(l)} \f$, i.e.,
 * the convex relaxation of \f$ c_l^2 + s_l^2 = w_{s(l)} w_{e(l)} \f$.
 * Finally, there are the bounds (1) of DCNetworkBlock on \f$ F^{fr}_l \f$
 * and, for an HVDC line, \f$ F^{fr}_l + F^{to}_l = 0 \f$ (no loss, since the
 * efficiency is ignored) and the bounds on its reactive flows, if given.
 * Therefore, the model is a convex relaxation of the AC optimal power flow,
 * whose solutions may have \f$ c_l^2 + s_l^2 < w_{s(l)} w_{e(l)} \f$, which
 * no voltage profile gives; no voltage profile is recovered from them [see
 * recover_feasible_solution()]. Note that the
 * formulations and the angles of DCNetworkBlock play no role, since the only
 * rows of DCNetworkBlock that an ACNetworkBlock writes are the bounds (1) of
 * DCNetworkBlock::generate_abstract_constraints(). */

class ACNetworkBlock : public DCNetworkBlock
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
 * ACNetworkBlock defines one main public type:
 *
 * - ACNetworkData, an auxiliary class extending DCNetworkData with all the
 *   data necessary to describe the AC version of the transmission network.
 * @{ */

/*--------------------------------------------------------------------------*/
/*------------------- CLASS ACNetworkBlock::ACNetworkData ------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
 /// auxiliary class holding basic data about the (AC) transmission network
 /** The ACNetworkData class is a nested sub-class which only serves to have
  * a quick way to load all the basic data (topology and electrical
  * characteristics) that describe the transmission network. It extends
  * DCNetworkBlock::DCNetworkData with the (numerous) data necessary to
  * represent the AC version of Kirchhoff's laws: the impedances of the
  * lines, their thermal limits and the bounds on the angles, the shunt
  * admittances of the nodes and the bounds on their voltages (see
  * deserialize()). */

class ACNetworkData : public DCNetworkData
{
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

 public:

/** @} ---------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/** @name Constructor and Destructor
 * @{ */

 /// constructor of ACNetworkData, does nothing
 ACNetworkData( void ) {}

 /// copy constructor of ACNetworkData, does nothing
 explicit ACNetworkData( const NetworkData * ) {}

 /// destructor of ACNetworkData: it is virtual, and empty
 ~ACNetworkData() override = default;

/** @} --------------------- OTHER INITIALIZATIONS -------------------------*/
/** @name Other initializations
 * @{ */

 /// deserialize an ACNetworkData out of a netCDF::NcGroup
 /** Deserialize an ACNetworkData out of a netCDF::NcGroup, which should
  * contain all the data of a DCNetworkData (see
  * DCNetworkData::deserialize()) and, if "NumberNodes" > 1, the following:
  *
  * - the attribute "baseMVA", a string holding the base power of the per
  *   unit system (in MVA), which divides the thermal limits and the shunt
  *   admittances of the nodes; if it is not there, or it does not hold a
  *   number, it is 1;
  *
  * - the variables "LineReactance" and "LineResistance", of type
  *   netCDF::NcDouble and indexed over "NumberLines", the reactance
  *   \f$ x_l \f$ and the resistance \f$ r_l \f$ of each line in per unit,
  *   which are mandatory, since they decide which lines are AC lines (those
  *   with a nonzero susceptance, reactance or resistance) and which are
  *   HVDC lines (the others);
  *
  * - the variables "LineChargingSusceptance", "LineRatio" and
  *   "LineShiftAngle", indexed over "NumberLines", the charging susceptance
  *   \f$ b_l \f$, the ratio \f$ \tau_l \f$ of a transformer and its phase
  *   shift \f$ \nu_l \f$ (in degrees) of each line, optional with defaults
  *   0, 1 and 0, which give the admittances of the line \f$ Y^{tt}_l =
  *   Y_l + i b_l / 2 \f$, \f$ Y^{ff}_l = Y^{tt}_l / \tau_l^2 \f$,
  *   \f$ Y^{ft}_l = - Y_l / ( \tau_l e^{- i \nu_l} ) \f$ and
  *   \f$ Y^{tf}_l = - Y_l / ( \tau_l e^{i \nu_l} ) \f$, with
  *   \f$ Y_l = 1 / ( r_l + i x_l ) \f$;
  *
  * - the variable "LineRATEA", indexed over "NumberLines", the thermal
  *   limit \f$ r^A_l \f$ of each AC line (in MVA, divided by "baseMVA");
  *
  * - the variables "LineMinAngle" and "LineMaxAngle", indexed over
  *   "NumberLines", the bounds \f$ \phi^{mn}_l \le \phi^{mx}_l \f$ on the
  *   difference of the angles of the voltages at the ends of each AC line,
  *   in degrees, optional but given together (one without the other is
  *   refused with a std::invalid_argument); if they are not there, the
  *   angle differences are not bounded: no row (8) and no bound derived
  *   from the angles is written [see
  *   ACNetworkBlock::generate_abstract_constraints()], and the strengthened
  *   relaxation, whose envelopes need a finite range of each angle
  *   difference, is not available [see
  *   ACNetworkBlock::generate_abstract_variables()];
  *
  * - the variables "NodeConductance" and "NodeSusceptance", indexed over
  *   "NumberNodes", the shunt conductance \f$ G^s_n \f$ and susceptance
  *   \f$ B^s_n \f$ of each node (in MW and MVAr at voltage 1, divided by
  *   "baseMVA");
  *
  * - the variables "NodeMinVoltage" and "NodeMaxVoltage", indexed over
  *   "NumberNodes", the bounds \f$ V^{mn}_n \f$ and \f$ V^{mx}_n \f$ on the
  *   voltage magnitude of each node;
  *
  * - the variables "MinReactivePowerFlow" and "MaxReactivePowerFlow",
  *   indexed over "NumberLines", the bounds on the reactive flows of the
  *   HVDC lines (those of the AC lines being ignored), optional.
  *
  * All the variables but "LineReactance" and "LineResistance" are read as
  * optional; those without a default and without a meaning when absent
  * ("LineRATEA" and the data of the nodes) are however used by
  * ACNetworkBlock::generate_abstract_constraints(), and therefore have to
  * be there if the abstract representation is generated. */

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG
 // extends [DC]NetworkData::expected_dims()
 /* not necessary since ACNetworkData does not have any new dimensions

 std::vector< std::string > expected_dims( void ) const override;
 */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// extends [DC]NetworkData::expected_vars()

 std::vector< std::string > expected_vars( void ) const override;

#endif

/** @} ------- METHODS FOR READING THE DATA OF THE ACNetworkData -----------*/
/** @name Reading the data of the ACNetworkData
 * @{ */

 /// returns the system MVA base of the instance

 double get_baseMVA( void ) const { return( f_base_mva ); }

/*--------------------------------------------------------------------------*/
 /// returns the vector of nodal shunt conductances

 const std::vector< double > & get_node_conductance( void ) const {
  return( v_node_conductance );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of nodal shunt susceptances

 const std::vector< double > & get_node_susceptance( void ) const {
  return( v_node_susceptance );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of nodal maximum voltages

 const std::vector< double > & get_node_max_voltage( void ) const {
  return( v_node_max_voltage );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of nodal minimum voltages

 const std::vector< double > & get_node_min_voltage( void ) const {
  return( v_node_min_voltage );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of line reactances

 const std::vector< double > & get_line_reactance( void ) const {
  return( v_line_reactance );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of line resistances

 const std::vector< double > & get_line_resistance( void ) const {
  return( v_line_resistance );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of line charging susceptances

 const std::vector< double > & get_line_chargingsusceptance( void ) const {
  return( v_line_chargingsusceptance );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of line (transformer) tap ratios

 const std::vector< double > & get_line_ratio( void ) const {
  return( v_line_ratio );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of line thermal limits (rate A)

 const std::vector< double > & get_line_rate_A( void ) const {
  return( v_line_rate_A );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of line phase-shift angles

 const std::vector< double > & get_line_angle( void ) const {
  return( v_line_angle );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of line minimum angle differences
 /** Returns the vector of the minimum angle differences of the lines, in
  * degrees, which is empty if the data do not bound the angle differences
  * (see has_angle_bounds()). */

 const std::vector< double > & get_line_min_angle( void ) const {
  return( v_line_min_angle );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of line maximum angle differences
 /** Returns the vector of the maximum angle differences of the lines, in
  * degrees, which is empty if the data do not bound the angle differences
  * (see has_angle_bounds()). */

 const std::vector< double > & get_line_max_angle( void ) const {
  return( v_line_max_angle );
  }

/*--------------------------------------------------------------------------*/
 /// returns true if the data bound the angle differences of the lines
 /** Returns true if "LineMinAngle" and "LineMaxAngle" have been read, in
  * which case get_line_min_angle() and get_line_max_angle() have one entry
  * per line; false if they were not there, in which case both vectors are
  * empty and the angle differences are not bounded. */

 bool has_angle_bounds( void ) const {
  return( ! v_line_min_angle.empty() );
  }

/*--------------------------------------------------------------------------*/
 /// returns true if minimum and maximum reactive flow bounds have been loaded

 bool has_reactive_bounds( void ) const {
  return( ( ! v_min_reac_power_flow.empty() ) &&
          ( ! v_max_reac_power_flow.empty() ) );
  }

/*--------------------------------------------------------------------------*/
 /// returns minimum reactive power flow of the given \p line
 /** This method returns the minimum reactive power flow of the given
  * \p line.
  *
  * @return the minimum reactive power flow of the given \p line. */

 double get_min_reac_power_flow( Index line ) const {
  assert( line < f_number_lines );
  if( v_min_reac_power_flow.empty() )
   return( 0 );
  return( v_min_reac_power_flow[ line ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns maximum reactive power flow of the given \p line
 /** This method returns the maximum reactive power flow of the given
  * \p line.
  *
  * @return the maximum reactive power flow of the given \p line. */

 double get_max_reac_power_flow( Index line ) const {
  assert( line < f_number_lines );
  if( v_max_reac_power_flow.empty() )
   return( 0 );
  return( v_max_reac_power_flow[ line ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the partition of AC lines per node into direct and reverse sets
 /** For each node n, returns a pair < direct , reverse > of sets of AC line
  * ids: the first set contains the lines for which n is the start bus, the
  * second the lines for which n is the end bus. */

 std::vector< std::pair< std::set< Index > , std::set< Index > > >
  get_direct_and_reverse_AClines( void ) {
  std::vector< std::pair< std::set< Index > , std::set< Index > > > v(
                                                       get_number_nodes() );
  const auto & start_line = get_start_line();
  const auto & end_line = get_end_line();
  for( auto & line_id : get_DC_lines() ) {
   Index p = start_line[ line_id ];
   Index n = end_line[ line_id ];
   v[ p ].first.insert( line_id );
   v[ n ].second.insert( line_id );
   }
  return( v );
  }

/** @} --------------- METHODS FOR SAVING THE ACNetworkData ----------------*/
/** @name Methods for loading, printing & saving the ACNetworkData
 * @{ */

 /// serialize an ACNetworkData out of a netCDF::NcGroup
 /** Serialize an ACNetworkData out of a netCDF::NcGroup to the specific
  * format of an ACNetworkData: the parent DCNetworkData first, then the
  * "baseMVA" attribute and the AC-specific line and node data. See
  * ACNetworkData::deserialize( netCDF::NcGroup ) for details of the format
  * of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/

 protected:

/*--------------------- PROTECTED METHODS OF THE CLASS ---------------------*/

/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/

 double f_base_mva;  ///< the reference MVA base of the instance

 /* AC networks come with a number of additional data:
  *
  * - each power line has a Resistance (r), Reactance (x) and
  *   Susceptance (B). Note that when these are given in per unit, they can
  *   be converted to physical values by computing f = V^2 / mbase and
  *   multiplying r, x by f, while dividing B by f.
  *
  * - power lines also have a "chargingSusceptance" with symbol (b). The
  *   formula for this quantity is 2 * pi * f * C, with C the capacitance of
  *   a power line and f the nominal frequency in Hz (e.g., 50). This
  *   charging susceptance plays a role in the "Reactive AC power flow
  *   equations".
  *
  * - the Susceptance is the datum that a DCNetworkBlock reads as well.
  *
  *   /!\ : we make the assumption that when all three (r, x, B) are zero,
  *         then the line is in fact HVDC.
  *
  * - some lines corresponding to Transformers also come with a "Tap" ratio.
  *   This is the ratio of nominal voltages on both ends of the line. The
  *   default value is therefore 1.0.
  *
  * - each line also has a thermal limit (rate_A), which is the bound on
  *   total flow, the equivalent of the classic MaxPowerFlow.
  *
  * - each line comes with a phase angle difference (nominally zero) and
  *   bounds on this difference. Typical default values for these bounds
  *   would be +/- 20 degrees, 30 being a sort of maximal value, indicating
  *   closeness to instability.
  *
  * - nodes that have shunts (typically not the case) have an extra
  *   conductance term Gs and susceptance term Bs. The default values are 0
  *   and 0.
  *
  * - finally each node has bounds on allowed voltages. In typical per unit
  *   style these are assumed to be 0.9 and 1.1. However, as SMS++ is unit
  *   agnostic, one could specify physical units in which case for a 220 kV
  *   node we would specify the values 198 and 242.
  *
  * - for power lines that are HVDC, one can specify additional maximal and
  *   minimal reactive power flow bounds; these are independent of the
  *   active bounds already available in DCNetworkData for HVDC lines only.
  *   For AC lines these are linked to the active flow through the thermal
  *   limit constraints, which are of the form
  *     active_flow^2 + reactive_flow^2 <= Thermal_limit^2.
  *
  * - in the mixed case, for simplicity, the bounds are also given for AC
  *   lines but not used (!!!) */

 std::vector< double > v_line_reactance;
 std::vector< double > v_line_resistance;
 std::vector< double > v_line_chargingsusceptance;
 std::vector< double > v_line_ratio;
 std::vector< double > v_line_rate_A;
 std::vector< double > v_line_angle;
 std::vector< double > v_line_min_angle;
 std::vector< double > v_line_max_angle;
 std::vector< double > v_node_conductance;
 std::vector< double > v_node_susceptance;
 std::vector< double > v_node_max_voltage;
 std::vector< double > v_node_min_voltage;

 /// vector to store the minimum reactive power flow at each line
 std::vector< double > v_min_reac_power_flow;

 /// vector to store the maximum reactive power flow at each line
 std::vector< double > v_max_reac_power_flow;

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

/*-------------------- PRIVATE FIELDS OF THE CLASS -------------------------*/

 SMSpp_insert_in_factory_h;

/*---------------------- PRIVATE METHODS OF THE CLASS ----------------------*/

 };  // end( class( ACNetworkData ) )

/** @} ---------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 * @{ */

 /// constructor of ACNetworkBlock
 /** Constructor of ACNetworkBlock, taking possibly a pointer of its
  * father Block. */

 explicit ACNetworkBlock( Block * f_block = nullptr )
  : DCNetworkBlock( f_block ) , b_strongSOCP( true ) {}

/*--------------------------------------------------------------------------*/
 /// destructor of ACNetworkBlock

 virtual ~ACNetworkBlock() override;

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 * @{ */

 /// deserialize an ACNetworkBlock out of a netCDF::NcGroup
 /** Deserialize an ACNetworkBlock out of a netCDF::NcGroup. All the
  * AC-specific data lives in the [AC]NetworkData object, which is
  * automatically deserialize()-d by the base class via
  * get_new_NetworkData(); hence this method dispatches to
  * DCNetworkBlock::deserialize() and then only reads the optional
  * "ReactiveDemand" variable, which mirrors the (one-dimensional)
  * "ActiveDemand" one. */

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
 /// extends DCNetworkBlock::serialize( netCDF::NcGroup )
 /** Serialize the ACNetworkBlock into the given netCDF::NcGroup: the base
  * DCNetworkBlock first (which takes care of the [AC]NetworkData), then the
  * "ReactiveDemand" variable (if any). */

 void serialize( netCDF::NcGroup & group ) const override;

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG
 // extends [DC]NetworkBlock::expected_dims()
 /* not necessary since ACNetworkBlock does not have any new dims save those
  * of the ACNetworkData that are automatically taken into account.

 std::vector< std::string > expected_dims( void ) const override;
 */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// extends [DC]NetworkBlock::expected_vars()

 std::vector< std::string > expected_vars( void ) const override;

#endif

/*--------------------------------------------------------------------------*/
 /// generate the abstract variables of the ACNetworkBlock
 /** This method generates the following variables:
  *
  *  - v_power_flow
  *  - v_reactive_power_flow
  *  - v_sum_product_voltages
  *  - v_diff_product_voltages
  *  - v_sqrd_voltages
  *
  * Observe that, unlike in a DCNetworkBlock, both v_power_flow (the real
  * part of flow through a line) and v_reactive_power_flow (the imaginary
  * part of flow through a line) have as dimension twice the total number of
  * lines. This is because we need to distinguish between flow to and from
  * buses. The flows are one group of static Variable, "v_power_flow_real",
  * which takes the place of the group "p_flow_network" of DCNetworkBlock,
  * so that each of them is given to the Solvers once.
  *
  * The further variables stand for products of the complex voltages
  * \f$ V_n \f$, with \f$ s = s(l) \f$ and \f$ e = e(l) \f$ for an AC line
  * \f$ l \f$:
  *
  *  - v_sum_product_voltages, \f$ c_l = \mathrm{Re}( V_s ) \mathrm{Re}( V_e )
  *    + \mathrm{Im}( V_s ) \mathrm{Im}( V_e ) \f$;
  *  - v_diff_product_voltages, \f$ s_l = \mathrm{Im}( V_s )
  *    \mathrm{Re}( V_e ) - \mathrm{Re}( V_s ) \mathrm{Im}( V_e ) \f$;
  *  - v_sqrd_voltages, \f$ w_n = | V_n |^2 \f$ for each node \f$ n \f$.
  *
  * These appear in the standard rotated second-order cones.
  *
  * The Configuration parameter \p stvv (or, if \p stvv is nullptr and
  * f_BlockConfig is not nullptr,
  * f_BlockConfig->f_static_variables_Configuration) may be either a
  * SimpleConfiguration< int > or a SimpleConfiguration< std::vector< int > >:
  * in either case the first value decides whether the variables needed for
  * the stronger SOCP relaxation are added (see
  * generate_strengthened_variables()), which is the case if it is positive
  * and also when neither Configuration is given. The strengthened
  * relaxation needs the bounds on the angle differences (see
  * ACNetworkData::has_angle_bounds()): if the data do not give them and the
  * network has an AC line, it is not generated when no Configuration asks
  * for it, while a Configuration that asks for it makes the method throw a
  * std::invalid_argument before any Variable is added. A
  * SimpleConfiguration< int >
  * is also read by DCNetworkBlock::generate_abstract_variables(), which adds
  * the variables of the formulation it selects (of KIRCHHOFF with any other
  * Configuration); they play no role in the rows of an ACNetworkBlock. */

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the additional variables for the stronger SOCP relaxation
 /** Adds the auxiliary variables (v_voltage, v_theta, v_alpha, v_beta,
  * v_z) needed by the strengthened SOCP relaxation, i.e., the QC
  * relaxation of the Optimal Power Flow of Coffrin, Hijazi and Van
  * Hentenryck (IEEE Transactions on Power Systems 31(4), 2016, see the
  * class description), which adds multiple McCormick inequalities to
  * strengthen the basic SOCP relaxation. The internal
  * boolean #b_strongSOCP (settable through a BlockConfig of
  * static_variables kind) toggles this on or off. */

 void generate_strengthened_variables( void );

/*--------------------------------------------------------------------------*/
 /// generate the abstract constraints of the ACNetworkBlock
 /** Generates the SOCP relaxation of the AC optimal power flow described in
  * the detailed description of the class; none if the network has a single
  * node. With \f$ C^v \f$, \f$ \epsilon^{AC} \f$, \f$ C^{sc} \f$ and
  * \f$ d \f$ the four parameters of the Configuration (see below),
  * \f$ P^{base} \f$ the base power "baseMVA",
  * \f$ Y^s_n = ( G^s_n + i B^s_n ) / P^{base} \f$ the shunt admittance of
  * node \f$ n \f$, and the other data of ACNetworkData::deserialize(), the
  * rows are the following.
  *
  * The active and reactive balances of each node \f$ n \f$
  * ("AC_power_flow_injection", the active ones first, see
  * change_active_demand_constraints()):
  * \f{align*}{
  *   - C^v S_n - \mathrm{Re}( Y^s_n ) w_n / C^v
  *     + \sum_{ l : s(l) = n } F^{fr}_l + \sum_{ l : e(l) = n } F^{to}_l
  *   &= - C^v D^{ac}_n \; , \tag{1} \\
  *   - C^v R_n - \mathrm{Im}( Y^s_n ) w_n / C^v
  *     + \sum_{ l : s(l) = n } Q^{fr}_l + \sum_{ l : e(l) = n } Q^{to}_l
  *   &= - C^v D^{re}_n \; , \tag{2}
  * \f}
  * the sums running over all the lines, \f$ D^{re}_n \f$ being the reactive
  * demand ("ReactiveDemand").
  *
  * The flows of each AC line \f$ l \f$, with \f$ s = s(l) \f$ and
  * \f$ e = e(l) \f$ ("AC_voltage_definition_const"), within
  * \f$ \pm \epsilon^{AC} \f$:
  * \f{align*}{
  *   C^{sc} \bigl( \mathrm{Re}( Y^{ff}_l ) w_s + \mathrm{Re}( Y^{ft}_l ) c_l
  *     + \mathrm{Im}( Y^{ft}_l ) s_l - C^v F^{fr}_l \bigr) &\approx 0 \; ,
  *     \tag{3} \\
  *   C^{sc} \bigl( - \mathrm{Im}( Y^{ff}_l ) w_s
  *     - \mathrm{Im}( Y^{ft}_l ) c_l + \mathrm{Re}( Y^{ft}_l ) s_l
  *     - C^v Q^{fr}_l \bigr) &\approx 0 \; , \tag{4} \\
  *   C^{sc} \bigl( \mathrm{Re}( Y^{tt}_l ) w_e + \mathrm{Re}( Y^{tf}_l ) c_l
  *     - \mathrm{Im}( Y^{tf}_l ) s_l - C^v F^{to}_l \bigr) &\approx 0 \; ,
  *     \tag{5} \\
  *   C^{sc} \bigl( - \mathrm{Im}( Y^{tt}_l ) w_e
  *     - \mathrm{Im}( Y^{tf}_l ) c_l - \mathrm{Re}( Y^{tf}_l ) s_l
  *     - C^v Q^{to}_l \bigr) &\approx 0 \; , \tag{6}
  * \f}
  * i.e., the real and imaginary parts of the complex powers
  * \f$ ( Y^{ff}_l )^* | V_s |^2 + ( Y^{ft}_l )^* V_s V_e^* \f$ and
  * \f$ ( Y^{tt}_l )^* | V_e |^2 + ( Y^{tf}_l )^* V_e V_s^* \f$, each
  * coefficient of the admittances (times \f$ C^{sc} \f$) being rounded to
  * \f$ d \f$ significant digits; here \f$ \approx 0 \f$ means that the
  * left-hand side lies in \f$ [ - \epsilon^{AC} , \epsilon^{AC} ] \f$, an
  * equality with the default \f$ \epsilon^{AC} = 0 \f$.
  *
  * The thermal limits of each AC line ("AC_thermal_limit_const"),
  * \f[
  *   ( F^{fr}_l )^2 + ( Q^{fr}_l )^2 \le ( C^v r^A_l / P^{base} )^2 \; ,
  *   \qquad
  *   ( F^{to}_l )^2 + ( Q^{to}_l )^2 \le ( C^v r^A_l / P^{base} )^2 \; ;
  *   \tag{7}
  * \f]
  * the bounds \f$ ( C^v V^{mn}_n )^2 \le w_n \le ( C^v V^{mx}_n )^2 \f$ of
  * each node ("AC_voltage_bounds_limit"); for each AC line, with
  * \f$ \phi^{mn}_l \f$ and \f$ \phi^{mx}_l \f$ in radians, the bounds on the
  * angle difference
  * \f[
  *   \tan( \phi^{mn}_l ) c_l \le s_l \le \tan( \phi^{mx}_l ) c_l
  *   \tag{8}
  * \f]
  * ("AC_angle_bounds_limit"), the bounds
  * \f$ \min\{ \cos | \phi^{mn}_l | , \cos | \phi^{mx}_l | \} ( C^v )^2
  * V^{mn}_s V^{mn}_e \le c_l \le ( C^v )^2 V^{mx}_s V^{mx}_e \f$ and
  * \f$ | s_l | \le \sin( \phi^{mx}_l - \phi^{mn}_l ) ( C^v )^2 V^{mx}_s
  * V^{mx}_e \f$ ("AC_elem_bounds"), both groups being absent if the data
  * do not bound the angle differences (see
  * ACNetworkData::has_angle_bounds()), the cone
  * \f[
  *   c_l^2 + s_l^2 - w_s w_e \le 0
  *   \tag{9}
  * \f]
  * ("AC_socp_const", see generate_SOCP_relaxation()) and the bounds
  * \f$ | c_l | , | s_l | \le ( C^v )^2 V^{mx}_s V^{mx}_e \f$ that (9) and
  * the voltage bounds imply. For each HVDC line, \f$ F^{fr}_l + F^{to}_l =
  * 0 \f$ ("HVDC_flow_links") and, if the data give them, the bounds of the
  * reactive flows \f$ Q^{fr}_l \f$ and \f$ Q^{to}_l \f$ times \f$ C^v \f$
  * ("Reactive_Flow_Bounds"). Finally the bounds (1) of
  * DCNetworkBlock::generate_abstract_constraints() on \f$ F^{fr}_l \f$
  * (see DCNetworkBlock::generate_bound_constraints()), and the bounds on
  * \f$ R_n \f$ given by set_min_reactive_node_injection() and
  * set_max_reactive_node_injection(). If #b_strongSOCP is true, also
  * strengthen_SOCP_relaxation() is called, which adds the bounds of the
  * strengthened relaxation; its McCormick inequalities are separated by
  * generate_dynamic_constraints().
  *
  * The parameters, which mainly serve the numerical stability of the
  * nonlinear rows, are given by the Configuration \p stcc or, if \p stcc is
  * nullptr and f_BlockConfig is not nullptr, by
  * f_BlockConfig->f_static_constraints_Configuration:
  *
  * - \f$ C^v \f$ (C_v_scal), 1 by default, which scales the voltages and
  *   the flows as a change of the per unit base would: \f$ w \f$, \f$ c \f$
  *   and \f$ s \f$ stand for \f$ ( C^v )^2 \f$ times the corresponding
  *   products of voltages, and the flows of the rows are \f$ C^v \f$ times
  *   the power flows;
  *
  * - \f$ \epsilon^{AC} \f$ (f_ACvS), 0 by default, the slack of (3)-(6);
  *
  * - \f$ C^{sc} \f$ (f_scale), 1 by default, a factor multiplying (3)-(6);
  *
  * - \f$ d \f$ (f_digits), 16 by default, the number of significant digits
  *   to which the coefficients of (3)-(6) are rounded (see round_sig()).
  *
  * The Configuration is:
  *
  * - a SimpleConfiguration< double > for setting C_v_scal alone;
  *
  * - a SimpleConfiguration< std::pair< double , double > > for setting
  *   C_v_scal and f_ACvS;
  *
  * - a SimpleConfiguration< std::vector< double > > of length up to 4 so
  *   that its first element is C_v_scal, the second f_ACvS, the third
  *   f_scale and the fourth f_digits (if the vector is shorter than 4 the
  *   non-present parameters are kept to their default value, if it is
  *   longer the extra numbers are ignored). */

 void generate_abstract_constraints( Configuration * stcc = nullptr )
  override;

/*--------------------------------------------------------------------------*/
 /// the size of a line of an ACNetworkBlock is not changed by a kappa
 /** Besides the bounds (1) of DCNetworkBlock, the limit of a line of an
  * ACNetworkBlock is in rows that a kappa does not reach (the thermal limit
  * and the bounds of the reactive flow), hence sizing a line this way is not
  * supported and these throw.
  *
  * Supporting it means saying what the size of such a line is and writing
  * the kappa into the rows that carry it, i.e., the thermal limit and the
  * bounds of the reactive flow, and answering whether the susceptance
  * follows the size, which it does not do linearly: a modelling choice
  * rather than a translation, left to whoever needs it. */

 void set_kappa( MF_dbl_it values , Subset && subset , bool ordered = false ,
                 c_ModParam issuePMod = eNoBlck ,
                 c_ModParam issueAMod = eNoBlck ) override {
  throw( std::logic_error( "ACNetworkBlock::set_kappa: sizing a line of an "
   "ACNetworkBlock is not supported" ) );
  }

/*--------------------------------------------------------------------------*/

 void set_kappa( MF_dbl_it values , Range rng = Range( 0 , Inf< Index >() ) ,
                 c_ModParam issuePMod = eNoBlck ,
                 c_ModParam issueAMod = eNoBlck ) override {
  throw( std::logic_error( "ACNetworkBlock::set_kappa: sizing a line of an "
   "ACNetworkBlock is not supported" ) );
  }

/*--------------------------------------------------------------------------*/
 /// change the active balance rows after a change of the active demand
 /** The active demand of node \f$ n \f$ is the right-hand side
  * \f$ - C^v D^{ac}_n \f$ of the active part of the power balance of
  * \f$ n \f$ (the row "AC_power_flow_injection" of index \f$ n \f$, see
  * generate_abstract_constraints()), which is the only row an
  * ACNetworkBlock writes the demand in: this method changes it for the
  * nodes in \p modified_nodes, as DCNetworkBlock::set_active_demand()
  * asks, and the rows of the formulations of DCNetworkBlock, which an
  * ACNetworkBlock does not have, are left alone. */

 void change_active_demand_constraints( c_Subset & modified_nodes ,
                                        c_ModParam issueAMod ) override;

/*--------------------------------------------------------------------------*/
 /// separate the McCormick strengthening inequalities as dynamic cuts
 /** Separates, at the current point, the McCormick valid inequalities that
  * strengthen the SOCP relaxation (the z / c / beta / s families), adding to
  * v_SOCP_cuts only those violated by more than the tolerance read from \p
  * dycc (a SimpleConfiguration< double >, or the .first of a
  * SimpleConfiguration< pair< double , int > >; default 1e-6). Does nothing
  * unless the strengthened relaxation is active (b_strongSOCP). */

 void generate_dynamic_constraints( Configuration * dycc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the SOCP relaxation constraints
 /** Generates the rotated second-order cone constraints linking
  * v_sum_product_voltages, v_diff_product_voltages and v_sqrd_voltages,
  * i.e., the SOCP relaxation of the non-convex equality
  *   c_{n,n'}^2 + s_{n,n'}^2 = c_{n,n} c_{n',n'} . */

 void generate_SOCP_relaxation( void );

/*--------------------------------------------------------------------------*/
 /// strengthen the SOCP relaxation with McCormick-like inequalities
 /** Adds the auxiliary McCormick constraints relating v_voltage, v_theta,
  * v_alpha, v_beta and v_z to v_sqrd_voltages, v_sum_product_voltages and
  * v_diff_product_voltages, following the QC relaxation of Coffrin, Hijazi
  * and Van Hentenryck (2016) and the envelopes of Hijazi, Coffrin and Van
  * Hentenryck (2017) [see the class description for the references]. */

 void strengthen_SOCP_relaxation( void );

/** @} ---------------------------------------------------------------------*/
/*---------- METHODS FOR READING THE DATA OF THE ACNetworkBlock ------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the ACNetworkBlock
 * @{ */

 /// return the vector of power losses on lines
 /** Returns, for each line, the sum of the values of its "from" and "to"
  * active flow variables, i.e., \f$ C^v \f$ times its active loss [see
  * generate_abstract_constraints()]. */

 std::vector< double > get_line_losses( void ) const override {
  std::vector< double > losses;
  Index number_lines = get_number_lines();
  for( int line_id = 0 ; line_id < number_lines ; ++line_id ) {
   losses.push_back( v_power_flow[ line_id ].get_value()
                     + v_power_flow[ number_lines + line_id ].get_value() );
   }
  return( losses );
  }

/*--------------------------------------------------------------------------*/
 /// returns true since ACNetworkBlock handles reactive power

 bool handles_reactive( void ) const override { return( true ); }

/** @} ---------------------------------------------------------------------*/
/*---------- METHODS FOR READING THE Variable OF THE ACNetworkBlock --------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the Variable of the ACNetworkBlock
 * @{ */

 /// returns the node injection reactive power variables
 /** ACNetworkBlock does handle reactive power, so the method is actually
  * implemented here. Note that \p interval is ignored since ACNetworkBlock
  * always covers a single interval only. */

 ColVariable * get_reactive_node_injection( Index interval = 0 ) override {
  if( v_reactive_node_injection.empty() )
   return( nullptr );
  return( v_reactive_node_injection.data() );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// const version of get_reactive_node_injection()

 const ColVariable * get_const_reactive_node_injection( Index interval = 0 )
  const override {
  if( v_reactive_node_injection.empty() )
   return( nullptr );
  return( v_reactive_node_injection.data() );
  }

/*--------------------------------------------------------------------------*/
 /// returns the reactive power flow variables
 /** Returns the reactive power flow variables. Since in the AC formulation
  * each line has *two* reactive power variables, the "from" and the "to"
  * ones, the method returns a (const reference to a) vector RPF of
  * 2 * get_number_lines() variables: for l < get_number_lines(), RPF[ l ]
  * is the "from" reactive power variable of line l, otherwise it is the
  * "to" reactive power variable of line l - get_number_lines(). */

 std::vector< ColVariable > & get_reactive_power_flow( void ) {
  return( v_reactive_power_flow );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// const version of get_reactive_power_flow()

 const std::vector< ColVariable > & get_const_reactive_power_flow( void )
  const {
  return( v_reactive_power_flow );
  }

/*--------------------------------------------------------------------------*/
 /// recover a feasible solution: not implemented
 /** Not implemented: it returns an empty vector, whatever the values of the
  * Variable. */

 std::vector< std::pair< double , double > > recover_feasible_solution( void );

/** @} ---------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Solution
 * @{ */

 /// returns a Solution representing the current solution of this NetworkBlock
 /** This method must construct and return a (pointer to a) Solution object
  * representing the current "solution state" of this NetworkBlock. This is
  * a ACNetworkBlockSolution extending DCNetworkBlockSolution (which in
  * turn extends NetworkBlockSolution) with the specific extra solution
  * information of ACNetworkBlock.
  *
  * The parameter for deciding which kind of Solution must be returned is a
  * single int value, coded bitwise. The format is the same as that of
  * DCNetworkBlock::get_Solution(), with the only relevant bits for
  * ACNetworkBlock being
  *
  * - bit 0 (& 1) means "store the node injection"
  *
  * - bit 1 (& 2) means "store the flow values"
  *
  * The first is actually managed in the base NetworkBlock, and the second
  * in DCNetworkBlock. In these, they are taken to mean "store the
  * *active* power node injection / flow values". In
  * ACNetworkBlockSolution, this is extended to "store the *reactive*
  * power node injection / flow values".
  *
  * This value is to be found as:
  *
  * - if solc is not nullptr and it is a SimpleConfiguration< int >, then
  *   it is solc->f_value;
  *
  * - otherwise, if f_BlockConfig is not nullptr,
  *   f_BlockConfig->f_solution_Configuration is not nullptr and it is a
  *   SimpleConfiguration< int >, then it is
  *   f_BlockConfig->f_solution_Configuration->f_value;
  *
  * - otherwise, it is 7 (save everything). */

 Solution * get_Solution( Configuration * solc = nullptr ,
                          bool emptys = true ) override;

/*--------------------------------------------------------------------------*/
 /// return the "appropriate" [AC]NetworkBlockSolution

 NetworkBlockSolution * new_Solution( void ) const override;

/** @} ---------------------------------------------------------------------*/
/*--------------- METHODS FOR MODIFYING THE ACNetworkBlock -----------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the ACNetworkBlock
 * @{ */

 void set_NetworkData( NetworkData * nd = nullptr ) override {
  if( nd && ( ! dynamic_cast< ACNetworkData * >( nd ) ) )
   throw( std::invalid_argument( "ACNetworkBlock::set_NetworkData: not "
                                 "an ACNetworkData" ) );

  DCNetworkBlock::set_NetworkData( nd );
  }

/*--------------------------------------------------------------------------*/
 /// method to set the ReactiveDemand
 /** ACNetworkBlock does handle reactive power, so the method is actually
  * implemented here. Note that ACNetworkBlock always covers one interval
  * only, hence we expect v[] to contain just get_number_nodes() elements. */

 void set_ReactiveDemand( const boost::multi_array< double , 2 > & v )
  override {
  if( v_ReactiveDemand.empty() )
   v_ReactiveDemand.assign( v[ 0 ].begin() , v[ 0 ].end() );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of reactive demands
 /** ACNetworkBlock does handle reactive power, so the method is actually
  * implemented here. Note that ACNetworkBlock always covers one interval
  * only, hence \p interval is ignored. */

 const double * get_reactive_demand( Index interval = 0 ) const override {
  if( v_ReactiveDemand.empty() )
   return( nullptr );
  return( v_ReactiveDemand.data() );
  }

/*--------------------------------------------------------------------------*/
 /// method to set the MinReactiveNodeInjection
 /** ACNetworkBlock does handle reactive power, so the method is actually
  * implemented here. Note that ACNetworkBlock always covers one interval
  * only, hence \p t is ignored. */

 void set_min_reactive_node_injection( double min_inj , Index node ,
                                       Index t ) override {
  if( v_MinReactiveNodeInjection.empty() )
   v_MinReactiveNodeInjection.resize( get_number_nodes() );
  v_MinReactiveNodeInjection[ node ] = min_inj;
  if( ! reactive_node_injection_bounds_const.empty() )  // the rows are there
   reactive_node_injection_bounds_const[ node ].set_lhs( min_inj );
  }

/*--------------------------------------------------------------------------*/
 /// method to set the MaxReactiveNodeInjection
 /** ACNetworkBlock does handle reactive power, so the method is actually
  * implemented here. Note that ACNetworkBlock always covers one interval
  * only, hence \p t is ignored. */

 void set_max_reactive_node_injection( double max_inj , Index node ,
                                       Index t ) override {
  if( v_MaxReactiveNodeInjection.empty() )
   v_MaxReactiveNodeInjection.resize( get_number_nodes() );
  v_MaxReactiveNodeInjection[ node ] = max_inj;
  if( ! reactive_node_injection_bounds_const.empty() )  // the rows are there
   reactive_node_injection_bounds_const[ node ].set_rhs( max_inj );
  }

/** @} ---------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED METHODS OF THE CLASS ---------------------*/
/*--------------------------------------------------------------------------*/

 ACNetworkData * get_new_NetworkData( void ) const override {
  return( new ACNetworkData() );
  }

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

/*---------------------------------- data ----------------------------------*/

 /// boolean: true if we wish to generate the stronger SOCP relaxation
 bool b_strongSOCP;

 /// minimum reactive production of the electrical generators
 std::vector< double > v_MinReactiveNodeInjection;

 /// maximum reactive production of the electrical generators
 std::vector< double > v_MaxReactiveNodeInjection;

 /// reactive demand
 std::vector< double > v_ReactiveDemand;

/*-------------------------------- variables -------------------------------*/

 /// reactive power injection for each interval at each node
 std::vector< ColVariable > v_reactive_node_injection;

 /** the reactive power flow variables (real part is the standard
  * v_power_flow variable).
  *
  * One counter-intuitive aspect of AC flow is that power lines will have
  * a flow in both directions that are not opposites of each other, hence
  * all power flow variables are twice the lines (flow "to" and "from"). */
 std::vector< ColVariable > v_reactive_power_flow;

 /* Generic variables for AC-OPF.
  *
  * The voltages in each node are complex numbers having a module
  * |V_n| which is the usual voltage on which bounds are imposed, and an
  * angle which will intervene in the equations as differences between
  * nodes connected by a power line. To this end the terms
  *
  *   c_{n,n'} = Re(V_n) Re(V_n') + Im(V_n) Im(V_n')
  *   s_{n,n'} = Im(V_n) Re(V_n') - Re(V_n) Im(V_n')
  *
  * appear, using the well-known trigonometric identities
  *
  *   cos(a) cos(b) = 0.5 ( cos(a+b) - cos(a-b) )
  *   sin(a) sin(b) = 0.5 ( cos(a-b) - cos(a+b) )
  *   sin(a) cos(b) = 0.5 ( sin(a+b) + sin(a-b) )
  *
  * Alternative representations can be derived making appear the angle
  * difference on the line (on which bounds are known, typically +/- 30
  * degrees). */
 std::vector< ColVariable > v_sum_product_voltages;
 ///< \f$ c_l = | V_{s(l)} | | V_{e(l)} | \cos( \theta_{s(l)} -
 ///< \theta_{e(l)} ) \f$, for the DC lines only [see v_dc_line_position]

 std::vector< ColVariable > v_diff_product_voltages;
 ///< \f$ s_l = | V_{s(l)} | | V_{e(l)} | \sin( \theta_{s(l)} -
 ///< \theta_{e(l)} ) \f$, for the DC lines only [see v_dc_line_position]

 std::vector< Index > v_dc_line_position;
 ///< the position of each line among the DC lines, Inf< Index >() for an
 ///< HVDC line, which has no v_sum_product_voltages and
 ///< v_diff_product_voltages

 std::vector< ColVariable > v_sqrd_voltages;
 ///< \f$ w_n = | V_n |^2 \f$ for each node \f$ n \f$

 /* Variables for strengthening the SOCP relaxation by adding McCormick-like
  * inequalities:
  *
  *   alpha   ~ cos( theta_i - theta_j )
  *   beta    ~ sin( theta_i - theta_j )
  *   theta   is the angle in each node
  *   voltage is the modulus of voltage in each node
  *   z       is the auxiliary variable used in the McCormick relaxation
  *           of V_n V_n'
  *
  * N.B.: the auxiliary variables v_alpha, v_beta and v_z are only defined
  *       for the AC lines. As a result, if any HVDC line is present, they
  *       will have a different indexing than the other terms in the
  *       equations, most notably v_sum_product_voltages and
  *       v_diff_product_voltages. */
 std::vector< ColVariable > v_voltage;
 std::vector< ColVariable > v_theta;
 std::vector< ColVariable > v_alpha;
 std::vector< ColVariable > v_beta;
 std::vector< ColVariable > v_z;

/*------------------------------- constraints ------------------------------*/

 /// the node injection reactive power bound constraints
 std::vector< BoxConstraint > reactive_node_injection_bounds_const;

 // ----- Generic constraints for AC-OPF -----
 std::vector< BoxConstraint > v_voltage_bounds_const;
 MAFRC v_angle_bounds_const;
 MAFRC v_voltage_definition_const;
 std::vector< FRowConstraint > v_thermal_limit;
 std::vector< FRowConstraint > v_flow_dc;

 MABC v_basic_bounds_const;

 // ----- Bounds on reactive flow in lines (HVDC only) -----
 std::vector< BoxConstraint > v_reactive_flow_bounds;

 // ----- Specific constraints for SOCP relaxation -----
 std::vector< FRowConstraint > v_socp_const;

 /* Finite box bounds on the SOCP variables v_sum_product_voltages (c) and
  * v_diff_product_voltages (s), implied by the cone v_socp_const
  * (c^2 + s^2 <= w_p w_n) together with the voltage bounds, hence valid at
  * every feasible point and non-cutting. */
 std::vector< BoxConstraint > v_sum_product_voltages_bounds;
 std::vector< BoxConstraint > v_diff_product_voltages_bounds;

 // ----- Constraints for the stronger SOCP relaxation -----
 std::vector< BoxConstraint > v_volt_bounds;
 std::vector< FRowConstraint > v_theta_bounds;
 /// absolute box bounds pinning the free angle gauge; v_theta only ever
 /// appears as differences theta_p - theta_n, so this cannot cut any profile
 std::vector< BoxConstraint > v_theta_box_bounds;
 std::vector< BoxConstraint > v_alpha_bounds;
 ///< alpha ~ cos( theta_i - theta_j )
 std::vector< BoxConstraint > v_beta_bounds;
 ///< beta  ~ sin( theta_i - theta_j )

 /// McCormick relaxation of the square term V_n^2
 std::vector< FRowConstraint > v_diag_const_1;
 std::vector< FRowConstraint > v_diag_const_2;

 /* Convex envelope of the cosine of the phase-angle difference x:
  *   alpha <= 1 - ( 1 - cos( xbar ) ) / xbar^2 * x^2
  *   alpha >= cos( xbar )
  * (H. L. Hijazi, C. Coffrin and P. Van Hentenryck, "Convex quadratic
  * relaxations for mixed-integer nonlinear programs in power systems",
  * Mathematical Programming Computation 9, 321-367, 2017,
  * doi:10.1007/s12532-016-0112-z). */
 std::vector< FRowConstraint > v_def_alpha_1;
 std::vector< FRowConstraint > v_def_alpha_2;

 /* The McCormick inequalities that strengthen the SOCP relaxation are valid
  * inequalities (cuts), not part of the core model: they are therefore
  * handled as dynamic constraints, separated on demand by
  * generate_dynamic_constraints() instead of being all materialised up
  * front (which on large multi-period instances would create hundreds of
  * thousands of rows). The families are:
  *   z_{n,n'}    : classic McCormick envelope of the product v_n v_n'
  *   c_{n,n'}    : McCormick relaxation of Re( W_{n,n'} ) from z and alpha
  *   beta_{n,n'} : convex envelope of sin( theta_n - theta_n' )
  *   s_{n,n'}    : McCormick relaxation of Im( W_{n,n'} ) from z and beta
  * (eqs. (23b)-(23c) of Coffrin, Hijazi and Van Hentenryck, 2016, see the
  * class description); all are linear in v_z,
  * v_alpha, v_beta, v_voltage, v_theta, v_sum_product_voltages and
  * v_diff_product_voltages. */
 std::list< FRowConstraint > v_SOCP_cuts;

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

 /// downcast accessor to the underlying ACNetworkData

 ACNetworkData * ND( void ) {
  return( static_cast< ACNetworkData * >( f_NetworkData ) );
  }

/*--------------------------------------------------------------------------*/

 /// nothing of its own to register in the methods factory
 /** The methods an ACNetworkBlock can be asked for by name are those of
  * DCNetworkBlock, registered by DCNetworkBlock::static_initialization():
  * they reach an ACNetworkBlock as well, being called on it as on the
  * DCNetworkBlock it is, and set_active_demand() is virtual. Registering
  * them again here, under the same names, would replace the adapter of
  * DCNetworkBlock with one casting to ACNetworkBlock, which is wrong on a
  * DCNetworkBlock that is not one, and which of the two survives would
  * depend on the order of the static initialization. This one is defined
  * all the same so that the factory does not call the inherited one, which
  * would register the same methods twice. */

 static void static_initialization( void ) {}

/*--------------------------------------------------------------------------*/

 };  // end( class( ACNetworkBlock ) )

/*--------------------------------------------------------------------------*/
/*-------------------- CLASS ACNetworkBlockSolution ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// a [DCNetworkBlock]Solution of an ACNetworkBlock
/** The ACNetworkBlockSolution class derives from DCNetworkBlockSolution,
 * and therefore from NetworkBlockSolution, and adds the other information
 * that is typical of the ACNetworkBlock, i.e.,
 *
 * - the "from" and "to" reactive power variables
 *
 * Note that one ACNetworkBlock covers one time instant, so these variables
 * do not need to be indexed over time instants, like those in
 * DCNetworkBlockSolution and unlike those of the base NetworkBlockSolution.
 */

class ACNetworkBlockSolution : public DCNetworkBlockSolution
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*------------------------------- FRIENDS ----------------------------------*/

 friend ACNetworkBlock;  ///< make ACNetworkBlock friend

/*---------- CONSTRUCTING AND DESTRUCTING ACNetworkBlockSolution -----------*/

 /// constructor, does nothing
 explicit ACNetworkBlockSolution( void ) : DCNetworkBlockSolution() {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// deserialize an ACNetworkBlockSolution from a netCDF::NcGroup

 void deserialize( const netCDF::NcGroup & group ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// deserialize an ACNetworkBlockSolution from a "global" netCDF::NcGroup

 void deserialize( const netCDF::NcGroup & group , size_t idx ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 ~ACNetworkBlockSolution() override = default;
 ///< destructor: it is virtual, and empty

/*------ METHODS DESCRIBING THE BEHAVIOR OF A ACNetworkBlockSolution ------*/

 void read( const Block * block ) override final;

 void write( Block * block ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// serialize an ACNetworkBlockSolution into a netCDF::NcGroup
 /** Serialize an ACNetworkBlockSolution into a netCDF::NcGroup. The format
  * is the one of DCNetworkBlockSolution, which includes the one of
  * NetworkBlockSolution, cf. the comments in those methods, which in
  * particular means that
  *
  *     "NumberNetworks" IS NOT REALLY NEEDED, BECAUSE "TotalNumberInstants"
  *     AND "EndInstant" ARE NOT REQUIRED SINCE ACNetworkBlock ALWAYS HAS
  *     DCNetworkBlock::get_number_intervals() == 1, AND ALL THE
  *     NetworkBlock IN \p group ARE SUPPOSED TO BE ACNetworkBlock
  *
  * Owing to DCNetworkBlockSolution, \p group must contain
  *
  * - The dimension "NumberLines" containing the number of lines in the
  *   transmission network. It is mandatory. Note that
  *
  *       ALL THE DCNetworkBlock MUST HAVE THE SAME NUMBER OF LINES
  *
  * Furthermore, \p group must contain the ACNetworkBlock-specific
  * information:
  *
  * - The variable "NodeInjectionReactive", of type netCDF::NcDouble and
  *   indexed over the dimension "NumberNodes"; NodeInjectionReactive[ n ]
  *   is the optimal value of the reactive node injection on node (bus) n.
  *   The variable is optional.
  *
  * - The variable "ReactiveFlowFromValue", of type netCDF::NcDouble and
  *   indexed over the dimension "NumberLines"; ReactiveFlowFromValue[ l ]
  *   is the optimal value of the reactive power "from" on line l. The
  *   variable is optional.
  *
  * - The variable "ReactiveFlowToValue", of type netCDF::NcDouble and
  *   indexed over the dimension "NumberLines"; ReactiveFlowToValue[ l ]
  *   is the optimal value of the reactive power "to" on line l. The
  *   variable is optional, but it must be there if ReactiveFlowFromValue
  *   is there. */

 void serialize( netCDF::NcGroup & group ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// serialize an ACNetworkBlockSolution into a "global" netCDF::NcGroup
 /** "nonstandard" version of serialize() that loads an ACNetworkBlockSolution
  * from a "global" netCDF::NcGroup, i.e., one where the solution information
  * of multiple ACNetworkBlock are stored together (to avoid performance
  * issues due to the fact that netCDF is not structured to work with a
  * large number of sub-NcGroup in a file). The format is the "nonstandard"
  * one of the corresponding DCNetworkBlockSolution and NetworkBlockSolution,
  * which in particular means that
  *
  *     "NumberNetworks" IS NOT REALLY NEEDED, BECAUSE "TotalNumberInstants"
  *     AND "EndInstant" ARE NOT REQUIRED SINCE ACNetworkBlock ALWAYS HAS
  *     ACNetworkBlock::get_number_intervals() == 1, AND ALL THE
  *     NetworkBlock IN \p group ARE SUPPOSED TO BE ACNetworkBlock
  *
  * Owing to DCNetworkBlockSolution, \p group must contain
  *
  * - The dimension "NumberLines" containing the number of lines in the
  *   transmission network. It is mandatory. Note that
  *
  *       ALL THE DCNetworkBlock MUST HAVE THE SAME NUMBER OF LINES
  *
  * Furthermore, \p group must contain the ACNetworkBlock-specific
  * information:
  *
  * - The variable "NodeInjectionReactive", of type netCDF::NcDouble and
  *   indexed over both the dimension "NumberNetworks" (which is the same
  *   as "TotalNumberInstants", that does not exist) and the dimension
  *   "NumberNodes"; NodeInjectionReactive[ idx ][ n ] is the optimal value
  *   of reactive node injection on node (bus) n for this ACNetworkBlock.
  *   The variable is optional.
  *
  * - The variable "ReactiveFlowFromValue", of type netCDF::NcDouble and
  *   indexed over both the dimension "NumberNetworks" (which is the same
  *   as "TotalNumberInstants", that does not exist) and the dimension
  *   "NumberLines"; ReactiveFlowFromValue[ idx ][ l ] is the optimal value
  *   of reactive power "from" on line l for this ACNetworkBlock. The
  *   variable is optional.
  *
  * - The variable "ReactiveFlowToValue", of type netCDF::NcDouble and
  *   indexed over both the dimension "NumberNetworks" (which is the same
  *   as "TotalNumberInstants", that does not exist) and the dimension
  *   "NumberLines"; ReactiveFlowToValue[ idx ][ l ] is the optimal value
  *   of reactive power "to" on line l for this ACNetworkBlock. The
  *   variable is optional, but it must be there if ReactiveFlowFromValue
  *   is there.
  *
  * Note that the variables are constructed when \p idx == 0 according to
  * the fact that the corresponding DCNetworkBlockSolution has or not been
  * Configure-d to hold them, which means that
  *
  *       ALL THE ACNetworkBlockSolution MUST HAVE BEEN Configure-d IN THE
  *       SAME WAY
  *
  * (although, technically, if some of the ACNetworkBlockSolution that
  * appears when \p idx > 0 is Configure-d with less information than that
  * when idx == 0 the code will not break, but there will be uninitialised
  * values in the netCDF). */

 void serialize( netCDF::NcGroup & group , size_t idx ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 ACNetworkBlockSolution * scale( double factor ) const override;

 void sum( const Solution * solution , double multiplier ) override;

 ACNetworkBlockSolution * clone( bool empty = false ) const override;

/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/

 void print( std::ostream & output ) const override {
  output << "ACNetworkBlockSolution [" << this << "]: " << std::endl;
  }

/*---------------------- PRIVATE PART OF THE CLASS -------------------------*/

 private:

/*---------------------------- PRIVATE FIELDS ------------------------------*/

 /// injection of reactive power at each node
 std::vector< double > v_node_injection_reactive;

 /// v_reactive_flow_from[ l ] = reactive power "from" on line l
 std::vector< double > v_reactive_flow_from;

 /// v_reactive_flow_to[ l ] = reactive power "to" on line l
 std::vector< double > v_reactive_flow_to;

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class( ACNetworkBlockSolution ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 }  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* __ACNetworkBlock */

/*--------------------------------------------------------------------------*/
/*-------------------- End File ACNetworkBlock.h ---------------------------*/
/*--------------------------------------------------------------------------*/
