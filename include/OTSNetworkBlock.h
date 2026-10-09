/*--------------------------------------------------------------------------*/
/*----------------------- File OTSNetworkBlock.h ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 *
 * Header file for class OTSNetworkBlock, which derives from DCNetworkBlock
 * and adds Optimal Transmission Switching (OTS) constraints. The OTS model
 * allows the solver to open or close individual transmission lines in each
 * time period, which may give cheaper dispatches.
 *
 * Four OTS formulations are provided:
 *
 * - Standard BigM: a single binary switching variable \f$ \zeta_l \f$ per
 *   line, with classical BigM relaxation of Kirchhoff's Voltage Law.
 *
 * - Directional BigM: two binary variables
 *   \f$ \zeta_l^+, \zeta_l^- \f$ per line encoding the flow direction, giving
 *   tighter flow bounds.
 *
 * - Elastic BigM: a binary \f$ \zeta_l \f$ plus a continuous
 *   \f$ \zeta_{1,l} \f$ (elastic variable), yielding a two-stage
 *   BigM relaxation with coefficients \f$ \alpha, \beta \f$ that tighten
 *   the LP relaxation.
 *
 * - Elastic Directional BigM: combines directional and elastic ideas
 *   (\f$ \zeta_l^+, \zeta_l^-, \zeta_{1,l} \f$) for the tightest LP
 *   relaxation among the four.
 *
 * The BigM coefficients are computed, for each line with nonzero
 * susceptance, as
 * \f[
 *   M_l = | \mathfrak{S}_l | \sum_{ k \in \mathcal{L}^{S} }
 *         \frac{ P^{mx}_k }{ | \mathfrak{S}_k | }
 * \f]
 *
 * The class also defines OTSNetworkData (extending DCNetworkData) which
 * adds per-line switching costs, with support for a single scalar cost
 * (shared by all lines) or per-line vector costs.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __OTSNetworkBlock
 #define __OTSNetworkBlock
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
/*----------------------- CLASS OTSNetworkBlock ----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// a DCNetworkBlock with Optimal Transmission Switching
/** OTSNetworkBlock derives from DCNetworkBlock and extends the Kirchhoff
 * formulation with switching (on/off) variables for each transmission line,
 * which allow the solver to decide, in each time period, which lines are
 * active and which are open.
 *
 * This class always uses the KIRCHHOFF formulation as the base DC model
 * (voltage angles + flows), since OTS requires the BigM relaxation of
 * Kirchhoff's Voltage Law, which needs explicit angle variables.
 *
 * ## Formulations
 *
 * One selects the OTS formulation via a SimpleConfiguration< int > passed to
 * generate_abstract_variables(), whose integer value is interpreted bitwise
 * as follows.
 *
 * - Bits 0-1 (mask 0x3): select the formulation method:
 *   - 0 = Standard BigM (kOTS_Standard)
 *   - 1 = Directional BigM (kOTS_Directional)
 *   - 2 = Elastic BigM (kOTS_Elastic)
 *   - 3 = Elastic Directional BigM (kOTS_ElasticDirectional)
 *
 * - Bit 2 (mask 0x4): not used.
 *
 * ## Design-variable coupling
 *
 * If a DesignNetworkBlock has set design variables on this Block (via
 * set_design_variables()), the switching variables are coupled with them by
 * \f[
 *   \zeta_l \le x_l \qquad \forall\, l \text{ with design variable}
 * \f]
 * and therefore a line that has not been built (\f$ x_l = 0 \f$) is
 * necessarily open.
 *
 * ## Standard BigM formulation (kOTS_Standard)
 *
 * This formulation has the variables \f$ \zeta_l \in \{0,1\} \f$ (1 = line
 * closed / on) and the constraints
 * \f{align}{
 *   F_l - \mathfrak{S}_l \Delta\theta_l &\le  M_l (1 - \zeta_l) \\
 *   F_l - \mathfrak{S}_l \Delta\theta_l &\ge -M_l (1 - \zeta_l) \\
 *   -P^{mx}_l \zeta_l \le F_l &\le  P^{mx}_l \zeta_l
 * \f}
 *
 * ## Directional BigM formulation (kOTS_Directional)
 *
 * Here the variables are \f$ \zeta_l^+ , \zeta_l^- \in \{0,1\} \f$
 * (positive/negative flow direction), and the constraints are
 * \f{align}{
 *   \zeta_l^+ + \zeta_l^- &\le 1 \\
 *   F_l &\le  P^{mx}_l \zeta_l^+ \\
 *   F_l &\ge -P^{mx}_l \zeta_l^- \\
 *   F_l - \mathfrak{S}_l \Delta\theta_l &\le
 *     M_l (1 - \zeta_l^+ - \zeta_l^-) \\
 *   F_l - \mathfrak{S}_l \Delta\theta_l &\ge
 *     -M_l (1 - \zeta_l^+ - \zeta_l^-)
 * \f}
 *
 * ## Elastic BigM formulation (kOTS_Elastic)
 *
 * In this case the variables are \f$ \zeta_l \in \{0,1\} \f$ and the
 * continuous \f$ \zeta_{1,l} \f$, which has no bounds of its own: it is
 * bounded below by the first row, and a value above 1 only tightens the
 * other rows (since \f$ \beta_l \geq 0 \f$), so that it never helps. The
 * elastic coefficients of each line are \f$ \alpha_l = P^{mx}_l / M_l \f$
 * and \f$ \beta_l = 1 - \alpha_l \f$ (\f$ \alpha_l = 1 \f$ and
 * \f$ \beta_l = 0 \f$ when \f$ M_l = 0 \f$), and the constraints are
 * \f{align}{
 *   \zeta_l &\le \zeta_{1,l} \\
 *   F_l - \mathfrak{S}_l \Delta\theta_l &\le  \alpha_l M_l (1-\zeta_l)
 *                                  + \beta_l M_l (1-\zeta_{1,l}) \\
 *   F_l - \mathfrak{S}_l \Delta\theta_l &\ge -\alpha_l M_l (1-\zeta_l)
 *                                  - \beta_l M_l (1-\zeta_{1,l}) \\
 *   -P^{mx}_l \zeta_l \le F_l &\le  P^{mx}_l \zeta_l
 * \f}
 *
 * ## Elastic Directional BigM formulation (kOTS_ElasticDirectional)
 *
 * This last formulation has the variables
 * \f$ \zeta_l^+, \zeta_l^- \in \{0,1\} \f$ and the continuous
 * \f$ \zeta_{1,l} \f$ (with no bounds of its own, as above), and the
 * constraints
 * \f{align}{
 *   \zeta_l^+ + \zeta_l^- &\le 1 \\
 *   \zeta_l^+ &\le \zeta_{1,l},\quad \zeta_l^- \le \zeta_{1,l} \\
 *   F_l &\le  P^{mx}_l \zeta_l^+ \\
 *   F_l &\ge -P^{mx}_l \zeta_l^- \\
 *   F_l - \mathfrak{S}_l \Delta\theta_l &\le
 *     \alpha_l M_l (1 - \zeta_l^+ - \zeta_l^-)
 *     + \beta_l M_l (1-\zeta_{1,l}) \\
 *   F_l - \mathfrak{S}_l \Delta\theta_l &\ge
 *     -\alpha_l M_l (1 - \zeta_l^+ - \zeta_l^-)
 *     - \beta_l M_l (1-\zeta_{1,l})
 * \f}
 *
 * ## BigM computation
 *
 * We compute the BigM values as
 * \f[
 *   M_l = | \mathfrak{S}_l | \sum_{ k \in \mathcal{L}^{S} }
 *         \frac{ P^{mx}_k }{ | \mathfrak{S}_k | }
 * \f]
 * where the sum runs over all the lines with nonzero susceptance (the set
 * \f$ \mathcal{L}^{S} \f$ of DCNetworkBlock) and \f$ P^{mx}_k \f$ is
 * "MaxPowerFlow" of line \f$ k \f$ at the instant of the Block. In all the
 * formulations, \f$ \Delta\theta_l = \theta_{s(l)} - \theta_{e(l)} \f$ is the
 * difference of the angles at the ends of line \f$ l \f$, and the flow bounds
 * \f$ \pm P^{mx}_l \zeta_l \f$ are symmetric, built on "MaxPowerFlow" alone,
 * without the constant \f$ \kappa_l \f$ and the scaling factor \f$ C^v \f$ of
 * DCNetworkBlock. In addition, each line has the bounds (1) of
 * DCNetworkBlock::generate_abstract_constraints(). Note that the switching
 * variables and the rows of this class concern only the lines with nonzero
 * susceptance.
 *
 * ## Objective
 *
 * This class adds to the objective of the parent one the switching costs (if
 * any),
 * \f[
 *   \sum_{l \in \mathcal{L}^{S}} c^{\mathrm{sw}}_l (1 - \zeta_l)
 * \f]
 * where the sum runs over the lines with nonzero susceptance, which are the
 * only ones that have switching variables. For the directional formulations,
 * we use \f$ \zeta_l = \zeta_l^+ + \zeta_l^- \f$. */

class OTSNetworkBlock : public DCNetworkBlock
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*--------------------- PUBLIC TYPES OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
 * @{ */

/*--------------------------------------------------------------------------*/
/*------------------------ CLASS OTSNetworkData ----------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// DCNetworkData extension carrying per-line switching costs
/** OTSNetworkData extends DCNetworkData with a vector of per-line switching
 * costs. When a line \f$ l \f$ is opened (switched off) in a given time
 * period, a cost \f$ c^{\mathrm{sw}}_l \f$ is incurred in the objective.
 *
 * The switching cost can be specified as:
 *
 * - a scalar (one netCDF value), in which case every line gets the same
 *   cost;
 *
 * - a vector indexed over NumberLines, giving individual per-line costs;
 *
 * - omitted entirely, in which case the default is 0 for all lines.
 *
 * The netCDF variable is named "SwitchingCost" and is of type
 * netCDF::NcDouble. */

class OTSNetworkData : public DCNetworkData
{
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

 public:

/** @} ---------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/** @name Constructor and Destructor
 * @{ */

 /// constructor of OTSNetworkData, does nothing special
 OTSNetworkData( void ) : DCNetworkData() {}

 /// destructor of OTSNetworkData: it is virtual, and empty
 ~OTSNetworkData() override = default;

/** @} --------------------- OTHER INITIALIZATIONS -------------------------*/
/** @name Other initializations
 * @{ */

 /// deserialize an OTSNetworkData out of a netCDF::NcGroup
 /** Deserializes the parent DCNetworkData first, then reads the optional
  * netCDF variable "SwitchingCost":
  *
  * - If the variable has dimension 1 (a single scalar), that value is
  *   replicated for every line.
  *
  * - If the variable has dimension NumberLines, each entry gives the
  *   switching cost for the corresponding line.
  *
  * - If the variable is absent, all switching costs default to 0. */

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/
 /// serialize an OTSNetworkData into a netCDF::NcGroup
 /** Serialize an OTSNetworkData into a netCDF::NcGroup: the parent
  * DCNetworkData first, then the optional "SwitchingCost" variable (which
  * is not written when all switching costs are zero); see
  * OTSNetworkData::deserialize( netCDF::NcGroup ) for details of the
  * format. */

 void serialize( netCDF::NcGroup & group ) const override;

/** @} ----------------------- ACCESSOR METHODS ----------------------------*/
/** @name Accessor methods
 * @{ */

 /// returns the switching cost for line \p l
 /** Returns the switching cost for the given line index. If the switching
  * cost vector is empty (all zeros), returns 0. */

 double get_switching_cost( Index l ) const {
  if( v_switching_cost.empty() )
   return( 0.0 );
  return( v_switching_cost[ l ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the full switching cost vector (may be empty if all zero)

 const std::vector< double > & get_switching_cost( void ) const {
  return( v_switching_cost );
  }

/*--------------------------------------------------------------------------*/
 /// returns true if any line has a nonzero switching cost

 bool has_switching_cost( void ) const {
  return( ! v_switching_cost.empty() );
  }

/** @} ----------------- PROTECTED PART OF THE CLASS -----------------------*/

 protected:

 /// per-line switching costs; empty means all zero
 std::vector< double > v_switching_cost;

/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/

 private:

/*-------------------- PRIVATE FIELDS OF THE CLASS -------------------------*/

 SMSpp_insert_in_factory_h;

/*---------------------- PRIVATE METHODS OF THE CLASS ----------------------*/

 };  // end( class( OTSNetworkData ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
 /// enum for the OTS formulation types
 /** The four OTS formulations, selected by bits 0-1 of the configuration
  * integer. */

 enum ots_formulation_type {
  kOTS_Standard           = 0 , ///< Standard BigM
  kOTS_Directional        = 1 , ///< Directional BigM
  kOTS_Elastic            = 2 , ///< Elastic BigM (two-stage)
  kOTS_ElasticDirectional = 3   ///< Elastic + Directional BigM
  };

/** @} ---------------------------------------------------------------------*/
/*--------------- CONSTRUCTOR AND DESTRUCTOR OF OTSNetworkBlock ------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and destructor
 * @{ */

 /// constructor of OTSNetworkBlock, takes a Block father
 /** Constructor: \p father is the Block to which this OTSNetworkBlock
  * belongs; it can be nullptr (default). */

 OTSNetworkBlock( Block * father = nullptr )
  : DCNetworkBlock( father ) , f_ots_type( kOTS_Standard ) {}

/*--------------------------------------------------------------------------*/
 /// destructor of OTSNetworkBlock

 ~OTSNetworkBlock() override;

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 * @{ */

 /// extends DCNetworkBlock::deserialize( netCDF::NcGroup )
 /** Deserializes the OTSNetworkBlock from a netCDF group. The format is
  * identical to DCNetworkBlock, with the additional fields read by
  * OTSNetworkData (in particular, "SwitchingCost"). */

 void deserialize( const netCDF::NcGroup & group ) override;

/** @} ---------------------------------------------------------------------*/
/*------------ METHODS FOR GENERATING THE abstract representation ----------*/
/*--------------------------------------------------------------------------*/
/** @name Generating the abstract representation
 * @{ */

 /// generate the abstract variables of the OTSNetworkBlock
 /** Generates all variables for the OTS model. This method:
  *
  * 1. Forces the KIRCHHOFF base formulation (ignoring the parent's
  *    formulation selection), since OTS requires voltage angle variables.
  *
  * 2. Reads the OTS configuration from a SimpleConfiguration< int >.
  *    The integer is interpreted bitwise:
  *    - bits 0-1 (mask 0x3): OTS formulation (see ots_formulation_type).
  *    - bit 2 (mask 0x4): reserved.
  *
  * 3. Depending on the selected formulation, creates:
  *    - kOTS_Standard: \f$ \zeta_l \in \{0,1\} \f$ for each DC line.
  *    - kOTS_Directional: \f$ \zeta_l^+, \zeta_l^- \in \{0,1\} \f$ for each
  *      DC line.
  *    - kOTS_Elastic: \f$ \zeta_l \in \{0,1\} \f$ and the continuous
  *      \f$ \zeta_{1,l} \f$ for each DC line.
  *    - kOTS_ElasticDirectional: \f$ \zeta_l^+, \zeta_l^- \in \{0,1\} \f$ and
  *      the continuous \f$ \zeta_{1,l} \f$ for each DC line.
  *
  * The Configuration is found as for the parent class:
  * - if \p stvv != nullptr and is a SimpleConfiguration< int >, use it;
  * - otherwise use f_BlockConfig->f_static_variables_Configuration. */

 void generate_abstract_variables( Configuration * stvv = nullptr )
  override;

/*--------------------------------------------------------------------------*/
 /// generate the abstract constraints of the OTSNetworkBlock
 /** Generates all constraints for the OTS model. This method completely
  * overrides (does not call) DCNetworkBlock::generate_abstract_constraints,
  * because the Kirchhoff KVL constraints are replaced by their BigM
  * relaxation.
  *
  * The constraints generated are, in order:
  *
  * 1. NetworkCost auxiliary constraints (inherited helper
  *    generate_network_cost_constraints()):
  *    \f$ V_l \ge F_l,\; V_l \ge -F_l \f$ if NetworkCost is present.
  *
  * 2. BigM KVL relaxation for each DC line, formulation-dependent
  *    (see the class-level documentation for the exact constraints per
  *    formulation).
  *
  * 3. OTS flow bounds for each DC line, formulation-dependent:
  *    - Standard/Elastic:
  *      \f$ -P^{mx}_l \zeta_l \le F_l \le P^{mx}_l \zeta_l \f$
  *    - Directional/ElasticDirectional:
  *      \f$ -P^{mx}_l \zeta_l^- \le F_l \le P^{mx}_l \zeta_l^+ \f$
  *
  * 4. Switching exclusivity (Directional/ElasticDirectional only):
  *    \f$ \zeta_l^+ + \zeta_l^- \le 1 \f$
  *
  * 5. Elastic precedence (Elastic/ElasticDirectional only):
  *    - Elastic: \f$ \zeta_l \le \zeta_{1,l} \f$
  *    - ElasticDirectional:
  *      \f$ \zeta_l^+ \le \zeta_{1,l},\; \zeta_l^- \le \zeta_{1,l} \f$
  *
  * 6. Design coupling (if design variables exist):
  *    \f$ \zeta_l \le x_l \f$ (or \f$ \zeta_l^+ + \zeta_l^- \le x_l \f$ for
  *    directional formulations).
  *
  * 7. Reference angles (inherited helper
  *    generate_reference_angle_constraint()): one angle fixed to 0 in each
  *    component of the lines with nonzero susceptance, as (13) of
  *    DCNetworkBlock::generate_abstract_constraints().
  *
  * 8. KCL node balance (inherited helper
  *    generate_node_balance_constraints()).
  *
  * 9. The flow bounds (1) of DCNetworkBlock::generate_abstract_constraints()
  *    via generate_bound_constraints(), for all the lines (the HVDC lines
  *    have no switching variable and no other bound). */

 void generate_abstract_constraints( Configuration * stcc = nullptr )
  override;

/*--------------------------------------------------------------------------*/
 /// the size of a line of an OTSNetworkBlock is not changed by a kappa
 /** Besides the bounds (1) of DCNetworkBlock, the limit of a line of an
  * OTSNetworkBlock is in rows that a kappa does not reach (the flow bounds
  * coupled with the switching variables and the big-M), hence sizing a line
  * this way is not supported and these throw.
  *
  * Supporting it means saying what the size of such a line is and writing
  * the kappa into the rows that carry it, i.e., the flow bounds coupled
  * with the switching variables and the big-M of the Kirchhoff rows, which
  * is the limit itself and has to follow it: a modelling choice rather than
  * a translation, left to whoever needs it. */

 void set_kappa( MF_dbl_it values , Subset && subset , bool ordered = false ,
                 c_ModParam issuePMod = eNoBlck ,
                 c_ModParam issueAMod = eNoBlck ) override {
  throw( std::logic_error( "OTSNetworkBlock::set_kappa: sizing a line of an "
   "OTSNetworkBlock is not supported" ) );
  }

/*--------------------------------------------------------------------------*/

 void set_kappa( MF_dbl_it values , Range rng = Range( 0 , Inf< Index >() ) ,
                 c_ModParam issuePMod = eNoBlck ,
                 c_ModParam issueAMod = eNoBlck ) override {
  throw( std::logic_error( "OTSNetworkBlock::set_kappa: sizing a line of an "
   "OTSNetworkBlock is not supported" ) );
  }

/*--------------------------------------------------------------------------*/
 /// generate the objective of the OTSNetworkBlock
 /** Generates the objective function. Extends the parent objective with
  * switching costs:
  * \f[
  *   \sum_{l \in \mathcal{L}^{S}} c^{\mathrm{sw}}_l (1 - \zeta_l)
  * \f]
  * over the lines with nonzero susceptance. For directional formulations
  * the effective switching variable is
  * \f$ \zeta_l = \zeta_l^+ + \zeta_l^- \f$, so the cost becomes
  * \f$ c^{\mathrm{sw}}_l (1 - \zeta_l^+ - \zeta_l^-) \f$.
  *
  * If all switching costs are zero, no additional terms are added. */

 void generate_objective( Configuration * objc = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
/*-------------- METHODS FOR MODIFYING THE OTSNetworkBlock -----------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the OTSNetworkBlock
 * @{ */

 void set_NetworkData( NetworkData * nd = nullptr ) override {
  // if there was a previous OTSNetworkData, and it was local, delete it
  if( f_NetworkData && f_local_NetworkData )
   delete( f_NetworkData );

  f_NetworkData = dynamic_cast< OTSNetworkData * >( nd );
  f_local_NetworkData = false;
  }

/** @} ---------------------------------------------------------------------*/
/*-------------- Methods for checking the OTSNetworkBlock -----------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for checking solution information in the OTSNetworkBlock
 * @{ */

 /// returns true if the current solution is (approximately) feasible
 /** Extends DCNetworkBlock::is_feasible() by also checking feasibility
  * of all OTS-specific variables and constraints (switching variables,
  * BigM KVL, flow bounds, exclusivity, precedence, design coupling). */

 bool is_feasible( bool useabstract = false ,
                   Configuration * fsbc = nullptr ) override;


/** @} ---------------------------------------------------------------------*/
/*---------- METHODS FOR READING THE DATA OF THE OTSNetworkBlock ----------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the OTSNetworkBlock
 * @{ */

 /// returns the OTS formulation type currently in use
 ots_formulation_type get_ots_formulation( void ) const {
  return( f_ots_type );
  }

/*--------------------------------------------------------------------------*/
 /// returns the BigM value for line \p l

 double get_big_M( Index l ) const {
  if( v_big_M.empty() )
   return( 0.0 );
  return( v_big_M[ l ] );
  }

/*--------------------------------------------------------------------------*/
 /// returns the full BigM vector

 const std::vector< double > & get_big_M( void ) const {
  return( v_big_M );
  }

/*--------------------------------------------------------------------------*/
 /// returns the switching variable \f$ \zeta_l \f$ (Standard/Elastic)
 /** Returns a const reference to the vector of binary switching variables.
  * Only meaningful for Standard and Elastic formulations. */

 const std::vector< ColVariable > & get_switching( void ) const {
  return( v_switching );
  }

/*--------------------------------------------------------------------------*/
 /// returns the positive switching variables z+_l (Directional/ElasDir)

 const std::vector< ColVariable > & get_switching_pos( void ) const {
  return( v_switching_pos );
  }

/*--------------------------------------------------------------------------*/
 /// returns the negative switching variables z-_l (Directional/ElasDir)

 const std::vector< ColVariable > & get_switching_neg( void ) const {
  return( v_switching_neg );
  }

/*--------------------------------------------------------------------------*/
 /// returns the elastic variables z1_l (Elastic/ElasticDirectional)

 const std::vector< ColVariable > & get_elastic( void ) const {
  return( v_elastic );
  }

/** @} ---------------------------------------------------------------------*/
/*-------------- METHODS FOR READING THE Constraint -----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the Constraint of the OTSNetworkBlock
 * @{ */

 /// returns the upper BigM KVL constraints

 const std::vector< FRowConstraint > & get_OTS_KVL_upper( void ) const {
  return( v_OTS_KVL_upper );
  }

/*--------------------------------------------------------------------------*/
 /// returns the lower BigM KVL constraints

 const std::vector< FRowConstraint > & get_OTS_KVL_lower( void ) const {
  return( v_OTS_KVL_lower );
  }

/** @} ---------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED METHODS OF THE CLASS ---------------------*/
/*--------------------------------------------------------------------------*/

 DCNetworkData * get_new_NetworkData( void ) const override {
  return( new OTSNetworkData() );
  }

/*--------------------------------------------------------------------------*/
 /// compute the BigM values of the class notes
 /** Computes \f$ M_l = | \mathfrak{S}_l | \sum_k P^{mx}_k /
  * | \mathfrak{S}_k | \f$ for every DC line, the sum running over all the
  * lines with nonzero susceptance, and stores the result in v_big_M. */

 void compute_big_M( void );

/*--------------------------------------------------------------------------*/
 /// generate OTS variables for Standard formulation
 void generate_OTS_Standard_variables( void );

 /// generate OTS variables for Directional formulation
 void generate_OTS_Directional_variables( void );

 /// generate OTS variables for Elastic formulation
 void generate_OTS_Elastic_variables( void );

 /// generate OTS variables for ElasticDirectional formulation
 void generate_OTS_ElasticDirectional_variables( void );

/*--------------------------------------------------------------------------*/
 /// generate BigM KVL constraints for Standard formulation
 void generate_OTS_Standard_constraints( void );

 /// generate BigM KVL + exclusivity constraints for Directional
 void generate_OTS_Directional_constraints( void );

 /// generate BigM KVL + precedence constraints for Elastic
 void generate_OTS_Elastic_constraints( void );

 /// generate all constraints for ElasticDirectional
 void generate_OTS_ElasticDirectional_constraints( void );

/*--------------------------------------------------------------------------*/
 /// generate flow bounds with switching variables
 /** Generates the flow bounds coupled with the switching variables. For
  * Standard/Elastic: \f$ -P^{mx}_l \zeta_l \le F_l \le P^{mx}_l \zeta_l \f$.
  * For Directional/ElasticDirectional:
  * \f$ -P^{mx}_l \zeta_l^- \le F_l \le P^{mx}_l \zeta_l^+ \f$.
  * Only applies to DC lines; HVDC lines get standard box bounds. */

 void generate_OTS_flow_bounds( void );

/*--------------------------------------------------------------------------*/
 /// generate design-coupling constraints
 /** If design variables exist, generates \f$ \zeta_l \le x_l \f$ (Standard/
  * Elastic) or \f$ \zeta_l^+ + \zeta_l^- \le x_l \f$ (Directional). */

 void generate_design_coupling_constraints( void );

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

 ots_formulation_type f_ots_type;  ///< the OTS formulation in use

/*---------------------------------- data ---------------------------------*/

 /// BigM values, one per DC line
 std::vector< double > v_big_M;

 /// elastic alpha coefficients: alpha_l = f_max_l / M_l
 std::vector< double > v_alpha;

 /// elastic beta coefficients: beta_l = 1 - alpha_l
 std::vector< double > v_beta;

/*-------------------------------- variables ------------------------------*/

 /// binary switching variables \f$ \zeta_l \f$ (Standard, Elastic)
 std::vector< ColVariable > v_switching;

 /// positive direction switching variables z+_l (Directional, ElasDir)
 std::vector< ColVariable > v_switching_pos;

 /// negative direction switching variables z-_l (Directional, ElasDir)
 std::vector< ColVariable > v_switching_neg;

 /// continuous elastic variables z1_l, without bounds (Elastic, ElasDir)
 std::vector< ColVariable > v_elastic;

/*------------------------------- constraints -----------------------------*/

 /// the upper BigM rows of the voltage law, one for each DC line
 std::vector< FRowConstraint > v_OTS_KVL_upper;

 /// the lower BigM rows of the voltage law, one for each DC line
 std::vector< FRowConstraint > v_OTS_KVL_lower;

 /// OTS flow bounds (coupled with switching variables)
 std::vector< FRowConstraint > v_OTS_flow_upper;
 std::vector< FRowConstraint > v_OTS_flow_lower;

 /// switching exclusivity: z+ + z- <= 1 (Directional formulations)
 std::vector< FRowConstraint > v_switching_exclusivity;

 /// elastic precedence: z <= z1 (or z+ <= z1, z- <= z1)
 std::vector< FRowConstraint > v_elastic_precedence;

 /// design coupling: z <= x (or z+ + z- <= x)
 std::vector< FRowConstraint > v_design_coupling;

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
 // static_initialization() needs be defined, even if void, for otherwise the
 // method of the base class DCNetworkBlock is called, which is private

 static void static_initialization( void ) {}
 
/*--------------------------------------------------------------------------*/

 };  // end( class( OTSNetworkBlock ) )

/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif  /* OTSNetworkBlock.h included */

/*--------------------------------------------------------------------------*/
/*-------------------- End File OTSNetworkBlock.h --------------------------*/
/*--------------------------------------------------------------------------*/
