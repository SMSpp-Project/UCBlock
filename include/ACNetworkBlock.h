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
 * \copyright &copy; by Antonio Frangioni, Quentin Jacquet, Wim van Ackooij
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
/** The ACNetworkBlock class derives from DCNetworkBlock and adds the
 * numerous Variable and Constraint necessary to represent the AC version
 * of Kirchhoff's laws by means of a Second-Order Cone Programming (SOCP)
 * relaxation, possibly strengthened with McCormick-like inequalities (see
 * Coffrin, Hijazi, Van Hentenryck, "The QC Relaxation: A Theoretical and
 * Computational Study on Optimal Power Flow", IEEE TPWRS 31(4), 2016).
 *
 * TO BE COMPLETED
 */

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
  * represent the AC version of Kirchhoff's laws ...
  *
  * TO BE COMPLETED
  */

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
  * contain the following:
  *
  * TO BE COMPLETED
  *
  * - the "baseMVA" scalar variable, of type netCDF::NcDouble,
  *                                          ^^^^^^^^^^^^^^^^
  *   THAT'S WHAT ONE WOULD EXPECT, BUT IT SEEMS IT'S RATHER A STRING???
  *   specifying the system MVA base used for converting power into per unit
  *   quantities (see Matpower) */

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

 const std::vector< double > & get_line_min_angle( void ) const {
  return( v_line_min_angle );
  }

/*--------------------------------------------------------------------------*/
 /// returns the vector of line maximum angle differences

 const std::vector< double > & get_line_max_angle( void ) const {
  return( v_line_max_angle );
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
  * - each power line now has a Resistance (r), Reactance (x) and
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
  * - the Susceptance was typically already specified when DCNetworks were
  *   used.
  *
  *   /!\ : we make the assumption that when all three (r, x, B) are zero,
  *         then the line is in fact HVDC.
  *
  * - some lines corresponding to Transformers also come with a "Tap" ratio.
  *   This is the ratio of nominal voltages on both ends of the line. The
  *   default value is therefore 1.0.
  *
  * - each line also has a thermal limit (rate_A), which is the bound on
  *   total flow - the equivalent of the classic MaxPowerFlow.
  *
  * - each line comes with a phase angle difference (nominally zero) and
  *   bounds on this difference. Typical default values for these bounds
  *   would be +/- 20°, 30 being a sort of maximal value, indicating close
  *   to instability.
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
  * Observe that contrary to before, both v_power_flow (the real part of
  * flow through a line) and v_reactive_power_flow (the imaginary part of
  * flow through a line) now have as dimension twice the total number of
  * lines. This is because we need to distinguish between flow to and from
  * buses.
  *
  * The further variables "correspond to" voltages in each node:
  *
  *  - v_sum_product_voltages  = c_{n,n'} = Re(V_n)Re(V_n') + Im(V_n)Im(V_n')
  *  - v_diff_product_voltages = s_{n,n'} = Im(V_n)Re(V_n') - Re(V_n)Im(V_n')
  *  - v_sqrd_voltages         = W_{n,n}  = c_{n,n} = |V_n|^2
  *
  * These appear in the standard rotated second-order cones.
  *
  * The Configuration parameter \p stvv (or, if \p stvv is nullptr and
  * f_BlockConfig is not nullptr, f_BlockConfig->f_static_variables_Configuration)
  * may be either a SimpleConfiguration< int > or a
  * SimpleConfiguration< std::vector< int > >: in either case a non-zero
  * first value toggles on the addition of the variables needed for the
  * stronger SOCP relaxation (see generate_strengthened_variables()). */

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the additional variables for the stronger SOCP relaxation
 /** Adds the auxiliary variables (v_voltage, v_theta, v_alpha, v_beta,
  * v_z) needed by the strengthened SOCP relaxation. We follow
  *
  *   C. Coffrin, H. L. Hijazi and P. Van Hentenryck, "The QC Relaxation:
  *   A Theoretical and Computational Study on Optimal Power Flow", IEEE
  *   Transactions on Power Systems, vol. 31, no. 4, pp. 3008-3018, 2016,
  *   doi: 10.1109/TPWRS.2015.2463111.
  *
  * This paper considers and adds multiple McCormick inequalities to
  * strengthen the basic SOCP relaxation. The internal boolean
  * #b_strongSOCP (settable through a BlockConfig of static_variables kind)
  * toggles this on or off. */

 void generate_strengthened_variables( void );

/*--------------------------------------------------------------------------*/
 /// generate the abstract constraints of the ACNetworkBlock
 /** Generates the SOCP relaxation of the AC OPF constraints, plus the
  * thermal limits, the (optional) reactive flow bounds for HVDC lines, the
  * voltage bounds, the angle bounds, the power flow conservation and the
  * SOCP cone constraints. If #b_strongSOCP is true, also calls
  * strengthen_SOCP_relaxation() to add the McCormick-like inequalities.
  *
  * The above nonlinear constraints in the formulation are way more
  * numerically unstable than standard linear constraints, and therefore
  * careful scaling is needed. This is accomplished by defining four
  * numerical quantities:
  *
  * TODO: COMMENT BETTER WHAT EACH OF THESE DOES
  *
  * - C_v_scal => default 1.0; this is the main variable, working in a
  *   similar fashion as the "usual" per unit transform. Essentially all
  *   voltages become tilde_Voltage = C_v_scal * Voltage and power flows
  *   become C_v_scal * v_power_flow;
  *
  * - f_ACvS => default 0.0 (slack for AC_voltage_definition_const);
  *
  * - f_scale => default 1.0 (scaling constant for the
  *   AC_voltage_definition_const equations);
  *
  * - f_digits => default 16, the number of digits in round_sig(). This
  *   helps round some of the admittance matrix data that appears in the
  *   constraint up to f_digits digits.
  *
  * Setting these to non-default values is possible with the Configuration
  * parameter, that is either \p stcc or, if f_BlockConfig is not nullptr,
  * f_BlockConfig->f_static_constraints_Configuration. If the result is
  * not nullptr, then it is a SimpleConfiguration< ... > which can contain
  * up to four numbers, i.e.,
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
  * v_diff_product_voltages, following the QC relaxation of Coffrin et al.
  * (2016) and Hijazi et al. (2017). */

 void strengthen_SOCP_relaxation( void );

/** @} ---------------------------------------------------------------------*/
/*---------- METHODS FOR READING THE DATA OF THE ACNetworkBlock ------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the ACNetworkBlock
 * @{ */

 /// return the vector of power losses on lines
 /** Returns the per-line active losses, computed as the sum of the "from"
  * and "to" active flow variables. */

 std::vector< double > get_line_losses( void ) {
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
 /// recover a (relaxed) feasible solution
 /** Warning: only a relaxed feasible solution is recovered. */

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
  * The voltages in each node are now complex numbers having a module
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
  * difference on the line (on which bounds are known, typically +/- 30°). */
 std::vector< ColVariable > v_sum_product_voltages;
 ///< c_{n,n'} = v_n v_n' cos( theta_n - theta_n' )

 std::vector< ColVariable > v_diff_product_voltages;
 ///< s_{n,n'} = v_n v_n' sin( theta_n - theta_n' )

 std::vector< ColVariable > v_sqrd_voltages;  ///< c_{n,n} = |V_n|^2

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

 // ----- Constraints for the stronger SOCP relaxation -----
 std::vector< BoxConstraint > v_volt_bounds;
 std::vector< FRowConstraint > v_theta_bounds;
 std::vector< BoxConstraint > v_alpha_bounds;
 ///< alpha ~ cos( theta_i - theta_j )
 std::vector< BoxConstraint > v_beta_bounds;
 ///< beta  ~ sin( theta_i - theta_j )

 /// McCormick relaxation of the square term V_n^2
 std::vector< FRowConstraint > v_diag_const_1;
 std::vector< FRowConstraint > v_diag_const_2;

 /* The terms V_n V_n' appear in later McCormick relaxations. Therefore
  * z_{n,n'} representing this term appears in the classic McCormick
  * relaxation. */
 std::vector< FRowConstraint > v_def_z_1;
 std::vector< FRowConstraint > v_def_z_2;
 std::vector< FRowConstraint > v_def_z_3;
 std::vector< FRowConstraint > v_def_z_4;

 /* The McCormick relaxation of some later terms requires the convex
  * envelope of the cosine and sine functions; see
  *
  *   Hijazi, H., Coffrin, C. & Hentenryck, P.V. "Convex quadratic
  *   relaxations for mixed-integer nonlinear programs in power systems",
  *   Math. Prog. Comp. 9, 321-367 (2017),
  *   https://doi.org/10.1007/s12532-016-0112-z.
  *
  * For the cosine of x (the phase angle difference) we get
  *
  *   alpha <= 1 - ( 1 - cos( xbar ) ) / xbar^2 * x^2
  *   alpha >= cos( xbar ) */
 std::vector< FRowConstraint > v_def_alpha_1;
 std::vector< FRowConstraint > v_def_alpha_2;

 /* The following appear in the classic McCormick relaxation of the
  * "double" convex relaxation of
  *   Re( V_n V_n' ) = < ( V_n V_n' )^M ( cos( theta_n - theta_n' )^C ) >^M
  * which is thus the McCormick relaxation of the product term
  *     z_{n,n'} alpha_{n,n'}
  * This product term is none other than c_{n,n'}; see eq. (23b) in
  * Coffrin (2016). */
 std::vector< FRowConstraint > v_def_c_1;
 std::vector< FRowConstraint > v_def_c_2;
 std::vector< FRowConstraint > v_def_c_3;
 std::vector< FRowConstraint > v_def_c_4;

 /* The convex envelope of the sine function (with x the phase angle
  * difference) is as follows:
  *
  *   beta <= cos( xbar/2 ) ( x - xbar/2 ) + sin( xbar/2 )
  *   beta >= cos( xbar/2 ) ( x + xbar/2 ) - sin( xbar/2 ) */
 std::vector< FRowConstraint > v_def_beta_1;
 std::vector< FRowConstraint > v_def_beta_2;

 /* The following appear in the classic McCormick relaxation of the
  * "double" convex relaxation of
  *   Im( V_n V_n' ) = < ( V_n V_n' )^M ( sin( theta_n - theta_n' )^S ) >^M
  * which is thus the McCormick relaxation of the product term
  *     z_{n,n'} beta_{n,n'}
  * This product term is none other than s_{n,n'}; see eq. (23c) in
  * Coffrin (2016). */
 std::vector< FRowConstraint > v_def_s_1;
 std::vector< FRowConstraint > v_def_s_2;
 std::vector< FRowConstraint > v_def_s_3;
 std::vector< FRowConstraint > v_def_s_4;

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

 static void static_initialization( void )
 {
  register_method< ACNetworkBlock , MF_dbl_it , Subset && , bool >(
   "DCNetworkBlock::set_active_demand" ,
   & ACNetworkBlock::set_active_demand );

  register_method< ACNetworkBlock , MF_dbl_it , Range >(
   "DCNetworkBlock::set_active_demand" ,
   & ACNetworkBlock::set_active_demand );
  }

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
