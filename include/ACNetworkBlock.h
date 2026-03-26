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
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Quentin Jacquet
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
/*------------------------ CLASS ACNetworkBlock ----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// implementation of the NetworkBlock concept for AC equations
/** The ACNetworkBlock class derives from DCNetworkBlock and adds it the
 * numerous Variable and Constraint necessary to represent the AC version
 * of Kirchhoff's laws ...
 *
 * TO BE COMPLETED
 *
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
/*----------------- CLASS ACNetworkBlock::ACNetworkData --------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
 /// auxiliary class holding basic data about the (AC) transmission network
 /** The DCNetworkData class is a nested sub-class which only serves to have a
  * quick way to load all the basic data (topology and electrical
  * characteristics) that describe the transmission network. It extends
  * DCNetworkBlock::DCNetworkData with the (numerous) data necessary to
  * epresent the AC version of Kirchhoff's laws ...
  *
  * TO BE COMPLETED
  *
  */

 class ACNetworkData : public DCNetworkData
 {
 /*----------------------- PUBLIC PART OF THE CLASS ------------------------*/

 public:

 /*---------------------- CONSTRUCTOR AND DESTRUCTOR -----------------------*/
 /** @name Constructor and Destructor
  * @{ */

 /// constructor of ACNetworkData, does nothing

 ACNetworkData( void ) {}

 /// copy constructor of ACNetworkData, does nothing

 explicit ACNetworkData( const NetworkData * ) {}

 /// destructor of ACNetworkData: it is virtual, and empty

 virtual ~ACNetworkData( ) override = default;

 /** @} -------------------- OTHER INITIALIZATIONS ------------------------*/

 /// deserialize a DCNetworkData out of a netCDF::NcGroup
 /** Deserialize a DCNetworkData out of a netCDF::NcGroup, which should
  * contain the following:
  *
  * TO BE COMPLETED
  *
  * - the "baseMVA" scalar variable, of type netCDF::NcDouble, 
  *                                          ^^^^^^^^^^^^^^^^
  *   THAT'S WHAT ONE WOULD EXPECT, BUT IS SEEMS IT'S RATHER A STRING???
  *   specifying the
  *   system MVA base used for converting power into per unit quantities 
  *   (see Matpower) */

 virtual void deserialize( const netCDF::NcGroup & group ) override;

 /*-------------------------------------------------------------------------*/

#ifndef NDEBUG
 // extends [DC]NetworkData::expected_dims()
 /* not necessary since ACNetworkData does not have any new dimensions

 std::vector< std::string > expected_dims( void ) const override;
 */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// extends [DC]NetworkData::expected_vars()

 std::vector< std::string > expected_vars( void ) const override;

#endif

 /** @} ------- METHODS FOR READING THE DATA OF THE ACNetworkData ----------*/
 /** @name Reading the data of the DCNetworkData
  * @{ */

 double get_baseMVA( void ) const { return( f_base_mva ); }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const std::vector< double > & get_node_conductance( void ) const {
  return( v_node_conductance );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const std::vector< double > & get_node_susceptance( void ) const {
  return( v_node_susceptance );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const std::vector< double > & get_node_max_voltage( void ) const {
  return( v_node_max_voltage );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const std::vector< double > & get_node_min_voltage( void ) const {
  return( v_node_min_voltage );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const std::vector< double > & get_line_reactance( void ) const {
  return( v_line_reactance );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const std::vector< double > & get_line_resistance( void ) const {
  return( v_line_resistance );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const std::vector< double > & get_line_ratio( void ) const {
  return( v_line_ratio );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const std::vector< double > & get_line_rate_A( void ) const {
  return( v_line_rate_A );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const std::vector< double > & get_line_angle( void ) const {
  return( v_line_angle );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const std::vector< double > & get_line_min_angle( void ) const {
  return( v_line_min_angle );
  }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 const std::vector< double > & get_line_max_angle( void ) const {
  return( v_line_max_angle );
  }

/*--------------------------------------------------------------------------*/
 /// returns true if minimum and maximum reactive flow bounds have been loaded
 bool has_reactive_bounds( ) const{
    return ( !v_min_reac_power_flow.empty() && !v_max_reac_power_flow.empty() );
 }


/*--------------------------------------------------------------------------*/
 /// returns minimum Reactive power flow of the given \p line
 /** This method returns the minimum Reactive power flow of the given \p line.
  *
  * @return the minimum Reactive power flow of the given \p line. */

 double get_min_reac_power_flow( Index line ) const {
  assert( line < f_number_lines );
  if( v_min_reac_power_flow.empty() )
   return( 0 );
  return( v_min_reac_power_flow[ line ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns maximum Reactive power flow of the given \p line
 /** This method returns the maximum Reactive power flow of the given \p line.
  *
  * @return the maximum Reactive power flow of the given \p line. */

 double get_max_reac_power_flow( Index line ) const {
  assert( line < f_number_lines );
  if( v_max_reac_power_flow.empty() )
   return( 0 );
  return( v_max_reac_power_flow[ line ] );
  }
  
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

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

/** @} -------------- METHODS FOR SAVING THE ACNetworkData ----------------*/
/** @name Methods for loading, printing & saving the ACNetworkData
 * @{ */

 /// serialize a ACNetworkData out of a netCDF::NcGroup
 /** Serialize a ACNetworkData out of a netCDF::NcGroup to the specific
  * format of a ACNetworkData. See
  * ACNetworkData::deserialize( netCDF::NcGroup ) for details of the format
  * of the created netCDF group.
  *
  * TODO: IMPLEMENT
  */

 void serialize( netCDF::NcGroup & group ) const override {
  throw( std::logic_error( "ACNetworkData::serialize() not implemented yet" )
	 );
  }

 /** @} ---------------- PROTECTED PART OF THE CLASS -----------------------*/

 protected:

 /*-------------------- PROTECTED METHODS OF THE CLASS ---------------------*/

 /*-------------------- PROTECTED FIELDS OF THE CLASS ----------------------*/

 double f_base_mva;  ///< the reference mva of the instance

 std::vector< double > v_line_reactance;
 std::vector< double > v_line_resistance;
 std::vector< double > v_line_ratio;
 std::vector< double > v_line_rate_A;
 std::vector< double > v_line_angle;
 std::vector< double > v_line_min_angle;
 std::vector< double > v_line_max_angle;
 std::vector< double > v_node_conductance;
 std::vector< double > v_node_susceptance;
 std::vector< double > v_node_max_voltage;
 std::vector< double > v_node_min_voltage;

 // Data on Reactive max and min power flow

  /// vector to store the minimum Reactive power flow at each line
 std::vector< double > v_min_reac_power_flow;

 /// vector to store the maximum Reactive power flow at each line
 std::vector< double > v_max_reac_power_flow;

 /*----------------------- PRIVATE PART OF THE CLASS -----------------------*/

 private:

 /*-------------------- PRIVATE METHODS OF THE CLASS -----------------------*/

 /*-------------------- PRIVATE FIELDS OF THE CLASS ------------------------*/

 SMSpp_insert_in_factory_h;

 /*-------------------------------------------------------------------------*/
 /*-------------------------------------------------------------------------*/

 }; // end( class( ACNetworkData ) )

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 * @{ */

 explicit ACNetworkBlock( Block * f_block = nullptr )
  : DCNetworkBlock( f_block ) , b_strongSOCP( true ) {}

/*--------------------------------------------------------------------------*/

 virtual ~ACNetworkBlock() override;

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 * @{ */

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG
 // extends [DC]NetworkBlock::expected_dims()
 /* not necessary since ACNetworkBlock does not have any new dims save those
  * of the ACNetworkData that are automatically taken into account.

 std::vector< std::string > expected_dims( void ) const override;
 */

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 // extends [DC]NetworkBlock::expected_vars()
 /* not necessary since ACNetworkBlock does not have any new vars save those
  * of the ACNetworkData that are automatically taken into account.

 std::vector< std::string > expected_vars( void ) const override;
 */
#endif

/*--------------------------------------------------------------------------*/
 /// generate the abstract variables of the ACNetworkBlock
 ///
 ///  We will generate the following variables
 ///    v_power_flow
 ///    v_reactive_power_flow
 ///    v_sum_product_voltages
 ///    v_diff_product_voltages
 ///    v_sqrd_voltages
 ///
 ///  Observe that contrary to before, both v_power_flow (the real part of flow through a line) and v_reactive_power_flow (the imaginary part of flow through a line)
 ///    Now have as dimension twice the total number of lines. This is because we need to distinguish between flow to and from buses.
 ///
 ///  The further variables "correspond to" voltages in each node
 ///    v_sum_product_voltages  = c_{n,n'} = Re(V_n)Re(V_n') + Im(V_n)Im(V_n')
 ///    v_diff_product_voltages = s_{n,n'} = Im(V_n)Re(V_n') - Re(V_n)Im(V_n')
 ///    v_sqrd_voltages = W_{n,n'} = c_{n,n} = |V_n|^2
 ///  These appear in the standard rotated second order cones.
 ///
 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

 /// The following function adds additional variables for the stronger SOCP relaxation.
 /// Herein we follow the following paper
 /// C. Coffrin, H. L. Hijazi and P. Van Hentenryck, "The QC Relaxation: A Theoretical and Computational Study on Optimal Power Flow," in IEEE Transactions on Power Systems, vol. 31, no. 4, pp. 3008-3018, July 2016, doi: 10.1109/TPWRS.2015.2463111.
 // this paper considers and adds multiple McCormick inequalities to strenghten the basic SOCP relaxation
 /// 
 /// the internal variable b_strongSOCP which can be set through a BlockConfig (static_variables kind) toggles this on or off
 void generate_strengthened_variables( void ); 

/*--------------------------------------------------------------------------*/
 /// generate the abstract constraints of the ACNetworkBlock
 /** TODO: comment
  *
  * The above nonlinear constraints in the formulation are way more
  * numerically instable than standard linear constraints, and therefore
  * careful scaling is needed. This is accomplished by defining
  * four numerical quantities:
  *
  * TODO: COMMENT BETTER WHAT EACH OF THESE DOES
  *
  * - C_v_scal => default 1.0
  *
  * - f_ACvS => default 0.0 (Slack for AC_voltage_definition_const)
  *
  * - f_scale => default 1.0 (scaling constant for the
  *              AC_voltage_definition_const equations)
  *
  * - f_digits => (default 16) the # of digits in round_sig (default?)
  *
  * Setting these to non-default values is possible with the Configuration
  * parameter, that is either \p stcc or, if f_BlockConfig is not nullptr,
  * f_BlockConfig->f_static_constraints_Configuration. If the result is not
  * nullptr, then is is a SimpleConfiguration< ... > which can contain up
  * to four numbers, i.e.,
  *
  * - a SimpleConfiguration< double > for setting C_v_scal alone;
  *
  * - a SimpleConfiguration< std::pair< double , double > > for setting
  *   C_v_scal and f_ACvS;
  *
  * - a SimpleConfiguration< std::vector< double > > of lenght up to 4 so
  *   that its first element is C_v_scal, the second f_ACvS, the third
  *   f_scale and the fourth f_digits (if the vector is shorter than 4
  *   the non-present parameters are kept to their default value, if it is
  *   longer the extra numbers are ignored). */

 void generate_abstract_constraints( Configuration * stcc = nullptr )
  override;

/*--------------------------------------------------------------------------*/

 void generate_SOCP_relaxation( void );

 void strengthen_SOCP_relaxation( void );

/** @} ---------------------------------------------------------------------*/
/*---------- METHODS FOR READING THE DATA OF THE ACNetworkBlock ------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the ACNetworkBlock
 * @{ */

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
  *  implemented here. Note that interval is ignored since ACNetworkBlock
  *  always covers a single interval only. */

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
  *  each line has *two* reactive power variables, the "from" and the "to"
  *  ones, the method returns a (const reference to a) vector RPF of 
  *  2 * get_number_lines() variables: for l <  get_number_lines(), RPF[ l ]
  *  is the "from" reactive power variable of line l, otherwise it is the
  *  "to" reactive power variable of line l - get_number_lines(). */

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
  // warning only a relaxed solution
  std::vector< std::pair< double , double > > recover_feasible_solution(
								      void );

/** @} ---------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for handling Solution
 * @{ */

 /// returns a Solution representing the current solution of this NetworkBlock
 /** This method must construct and return a (pointer to a) Solution object
  * representing the current "solution state" of this NetworkBlock. This is
  * a ACNetworkBlockSolution extending DCNetworkBlockSolution (which in turn
  * extends NetworkBlockSolution) with the specific extra solution
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
  * The first is actually managed in the base NetworkBlock, and the second in
  * DCNetworkBlock. In these, they are taken to mean "store the *active*
  * power node injection / flow values". In ACNetworkBlockSolution, this
  * is extended to "store the *reactive* power node injection / flow values". 
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
				  " an ACNetworkData" ) );

   DCNetworkBlock::set_NetworkData( nd );
   }

/*--------------------------------------------------------------------------*/
 /// method to set the ReactiveDemand
 /** ACNetworkBlock does handle reactive power, so the method is actually
  *  implemented here. Note that ACNetworkBlock always covers one interval
  *  only, hence we expect v[] to contain just get_number_nodes() elements. */

 void set_ReactiveDemand( const boost::multi_array< double , 2 > & v )
  override {
  if( v_ReactiveDemand.empty() )
   v_ReactiveDemand.assign( v[ 0 ].begin() , v[ 0 ].end() );
  }

/*--------------------------------------------------------------------------*/
 /// method to set the MinReactiveNodeInjection
 /** ACNetworkBlock does handle reactive power, so the method is actually
  *  implemented here. Note that ACNetworkBlock always covers one interval
  *  only, hence \p t is ignored. */

 void set_min_reactive_node_injection( double min_inj , Index node ,
				       Index t ) override {
  if( v_MinReactiveNodeInjection.empty() )
   v_MinReactiveNodeInjection.resize( get_number_nodes() );
  v_MinReactiveNodeInjection[ node ] = min_inj;
  }

/*--------------------------------------------------------------------------*/
 /// method to set the MaxReactiveNodeInjection
 /** ACNetworkBlock does handle reactive power, so the method is actually
  *  implemented here. Note that ACNetworkBlock always covers one interval
  *  only, hence \p t is ignored.  */

 void set_max_reactive_node_injection( double max_inj , Index node ,
				       Index t ) override {
  if( v_MaxReactiveNodeInjection.empty() )
   v_MaxReactiveNodeInjection.resize( get_number_nodes() );
  v_MaxReactiveNodeInjection[ node ] = max_inj;
  }

/*--------------------------------------------------------------------------*/
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

 /// Boolean on if we wish to generate the stronger SOCP relaxation
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

  /// real part is the standard "v_power_flow" variable
  std::vector< ColVariable > v_reactive_power_flow;

  // ----- Generic variables for AC-OPF
  std::vector< ColVariable > v_sum_product_voltages;
  std::vector< ColVariable > v_diff_product_voltages;
  std::vector< ColVariable > v_sqrd_voltages;

  /// ----- Variables for stenghtening the SOCP relaxation by adding McCormick like inequalities
  std::vector< ColVariable > v_voltage;
  std::vector< ColVariable > v_theta;
  std::vector< ColVariable > v_alpha;
  std::vector< ColVariable > v_beta;
  std::vector< ColVariable > v_z;

/*------------------------------- constraints ------------------------------*/

  /// the node injection reactive power bound constraints
  std::vector< BoxConstraint > reactive_node_injection_bounds_const;

  // ----- Generic constraints for AC-OPF
  std::vector< BoxConstraint > v_voltage_bounds_const;
  MAFRC v_angle_bounds_const;
  MAFRC v_voltage_definition_const;
  std::vector< FRowConstraint > v_thermal_limit;
  std::vector< FRowConstraint > v_flow_dc;

  MAFRC v_basic_bounds_const;

  // ----- Bounds on Reactive flow in lines (HVDC only)
  std::vector< FRowConstraint > v_reactive_flow_bounds;

  // ----- Specific constraints for SOCP relaxation
  std::vector< FRowConstraint > v_socp_const;

  /// ----- Variables for stenghtening the SOCP relaxation by adding McCormick like inequalities
  std::vector< FRowConstraint > v_volt_bounds;  
  std::vector< FRowConstraint > v_theta_bounds;  
  std::vector< FRowConstraint > v_alpha_bounds; // alpha ~ cos( theta_i - theta_j )
  std::vector< FRowConstraint > v_beta_bounds;  // beta  ~ sin( theta_i - theta_j )

  /// -----  Various constraints for the stronger SOCP relaxation

  std::vector< FRowConstraint > v_diag_const_1;
  std::vector< FRowConstraint > v_diag_const_2;
  std::vector< FRowConstraint > v_def_alpha_1;
  std::vector< FRowConstraint > v_def_alpha_2;
  std::vector< FRowConstraint > v_def_beta_1;
  std::vector< FRowConstraint > v_def_beta_2;
  std::vector< FRowConstraint > v_def_z_1;
  std::vector< FRowConstraint > v_def_z_2;
  std::vector< FRowConstraint > v_def_z_3;
  std::vector< FRowConstraint > v_def_z_4;
  std::vector< FRowConstraint > v_def_c_1;
  std::vector< FRowConstraint > v_def_c_2;
  std::vector< FRowConstraint > v_def_c_3;
  std::vector< FRowConstraint > v_def_c_4;
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

 ACNetworkData * ND( void ) {
  return( static_cast< ACNetworkData * >( f_NetworkData ) );
  }

/*--------------------------------------------------------------------------*/

 static void static_initialization( void )
 {
  register_method< ACNetworkBlock , MF_dbl_it , Subset && , bool >(
   "DCNetworkBlock::set_active_demand" , & ACNetworkBlock::set_active_demand );

  register_method< ACNetworkBlock , MF_dbl_it , Range >(
   "DCNetworkBlock::set_active_demand" , & ACNetworkBlock::set_active_demand );
 }

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 }; // end( class( ACNetworkBlock ) )

/*--------------------------------------------------------------------------*/
/*-------------------- CLASS ACNetworkBlockSolution ------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// a [DCNetworkBlock]Solution of a ACNetworkBlock
/** The ACNetworkBlockSolution class derives from DCNetworkBlockSolution,
 * and therefore from NetworkBlockSolution, and adds the other information
 * that is typical of the ACNetworkBlock, i.e.,
 *
 * - the "from" and "to" reactive power variables
 *
 * Note that one ACNetworkBlock covers one time instant, so these variables
 * do not need to be indexed over time instants, like these in
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
  * Furthermore,  \p group must contain the ACNetworkBlock-specific
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
 /** "nonstandard" version of serialize() that loads a ACNetworkBlockSolution
  * from a "global" netCDF::NcGroup, i.e., one where the solution information
  * of multiple ACNetworkBlock are stored together (to avoid performance
  * issues due to the fact that netCDF is not structured to work with a large
  * number of sub-NcGroup in a file). The format is the  "nonstandard" one of
  * the corresponding DCNetworkBlockSolution and NetworkBlockSolution, which
  * in particular means that
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
  * Furthermore,  \p group must contain the ACNetworkBlock-specific
  * information:
  *
  * - The variable "NodeInjectionReactive", of type netCDF::NcDouble and
  *   indexed over the both the dimension "NumberNetworks" (which is the same
  *   as "TotalNumberInstants", that does not exist) and the dimension
  *   "NumberNodes"; NodeInjectionReactive[ idx ][ n ] is the optimal value
  *   of reactivenode injection on node (bus) n for this ACNetworkBlock. The
  *   variable is optional.
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

 void print( std::ostream &output ) const override {
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

} // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* __ACNetworkBlock */

/*--------------------------------------------------------------------------*/
/*-------------------- End File ACNetworkBlock.h ---------------------------*/
/*--------------------------------------------------------------------------*/
