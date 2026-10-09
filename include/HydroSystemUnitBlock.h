/*--------------------------------------------------------------------------*/
/*----------------------- File HydroSystemUnitBlock.h ----------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 *
 * Header file for the class HydroSystemUnitBlock, which derives from the
 * Block, in order to define a base class for any possible "hydro unit" plus
 * the linking PolyhedralFunctionBlock to describe the future value of water
 * function in a UCBlock.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Ali Ghezelsoflu \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Rafael Durbano Lobato \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Ali Ghezelsoflu,
 *                      Rafael Durbano Lobato
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __HydroSystemUnitBlock
 #define __HydroSystemUnitBlock
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Block.h"

#include "PolyhedralFunctionBlock.h"

#include "HydroUnitBlock.h"

#include "FRealObjective.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)

namespace SMSpp_di_unipi_it
{
/*--------------------------------------------------------------------------*/
/*----------------------- CLASS HydroSystemUnitBlock -----------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// a UnitBlock made of HydroUnitBlock and of the future cost of their water
/** HydroSystemUnitBlock, which derives from UnitBlock, groups a set of
 * HydroUnitBlock (typically the valleys whose water is valued jointly)
 * together with the future cost of the water left in their reservoirs at the
 * end of the horizon, as a single UnitBlock of a UCBlock (see
 * \ref ucblock_model). Its sub-Block are
 *
 * - the HydroUnitBlock, in the order of the groups "HydroUnitBlock_0",
 *   "HydroUnitBlock_1", ..., whose generators are, in this order, the
 *   generators of the HydroSystemUnitBlock;
 *
 * - possibly a PolyhedralFunctionBlock, the last sub-Block, whose
 *   PolyhedralFunction is the future cost of the water.
 *
 * Note that the HydroSystemUnitBlock has no Variable and no Constraint
 * outside its sub-Block, its Objective is identically zero, and therefore its
 * cost is the sum of the Objective of its sub-Block, while the scale factor
 * of UnitBlock is not implemented, i.e., get_scale() is 1.
 *
 * Future cost of the water. If \f$ v^{f} \f$ is the vector of the final
 * volumes \f$ v^{hy}_{n,T-1} \f$ of all the reservoirs of all the
 * HydroUnitBlock, in the order described in deserialize(), the
 * PolyhedralFunction is
 * \f[
 *   \check\nu( v^{f} ) = \max_{ j \in \mathcal{J} } \bigl\{ \alpha_j +
 *     \beta_j^\top v^{f} \bigr\} \; ,
 * \f]
 * where \f$ \mathcal{J} \f$ indexes its rows (cuts). A PolyhedralFunction can
 * also have a constant lower bound and vertical rows, i.e., constraints
 * \f$ \alpha_j + \beta_j^\top v^{f} \leq 0 \f$ on its domain (see
 * PolyhedralFunction). Since \f$ \check\nu \f$ is added to the cost of the
 * UCBlock, which is minimized, it must be convex (a maximum, not a minimum,
 * of affine functions), and a concave one is rejected by deserialize(). It is
 * a cost, which is nonincreasing in the volumes when the water left in the
 * reservoirs has a value. In particular, a single row with
 * \f$ \beta_j = - \omega \f$, the vector of the values \f$ \omega \geq 0 \f$
 * of a unit of water in each reservoir, values the water left at the end of
 * the horizon linearly, i.e., \f$ \check\nu( v^{f} ) = \alpha_j - \omega^\top
 * v^{f} \f$, where the constant is irrelevant (e.g., minus the value of the
 * initial water). When the PolyhedralFunctionBlock is linearized,
 * \f$ \check\nu \f$ is represented by an epigraph variable \f$ \phi \f$ and
 * the rows \f$ \alpha_j + \beta_j^\top v^{f} \leq \phi \f$ (see
 * PolyhedralFunctionBlock).
 *
 * Volume-dependent values of the water. A common case is a value of the water
 * that decreases as the reservoir fills: the final volume
 * \f$ v_n = v^{hy}_{n,T-1} \f$ of reservoir \f$ n \f$ is split by the levels
 * \f$ Y_{n,0} < Y_{n,1} < \cdots
 * < Y_{n,M_n} \f$ into the compartments \f$ [ Y_{n,k-1} , Y_{n,k} ) \f$,
 * each with a value \f$ \omega_{n,k} \f$ per unit of volume and
 * \f$ \omega_{n,1} > \omega_{n,2} > \cdots > \omega_{n,M_n} \geq 0 \f$.
 * The value of the water is the maximum of \f$ \sum_k \omega_{n,k}
 * g_{n,k} \f$ over the volumes \f$ g_{n,k} \f$ in the compartments, subject
 * to \f$ \sum_{ k = 1 }^{ M_n } g_{n,k} = v_n - Y_{n,0} \f$ and
 * \f$ 0 \leq g_{n,k} \leq Y_{n,k} - Y_{n,k-1} \f$. Since the values
 * decrease, the compartments fill from the bottom, and the value is
 * \f[
 *   W_n( v_n ) = \sum_{ k = 1 }^{ M_n } \omega_{n,k} \min \bigl\{
 *     \max \{ v_n - Y_{n,k-1} , 0 \} , Y_{n,k} - Y_{n,k-1} \bigr\} \; ,
 * \f]
 * a concave piecewise-linear function on \f$ [ Y_{n,0} , Y_{n,M_n} ] \f$. On
 * compartment \f$ k \f$ it is \f$ \omega_{n,k} ( v_n - Y_{n,k-1} ) + \sum_{ i
 * < k } \omega_{n,i} ( Y_{n,i} - Y_{n,i-1} ) \f$, and each of these affine
 * functions is above \f$ W_n \f$ elsewhere, since the slopes decrease. Hence,
 * the future cost \f$ - W_n \f$ is the convex function
 * \f[
 *   - W_n( v_n ) = \max_{ k = 1 , \ldots , M_n } \bigl\{ \alpha_{n,k} -
 *     \omega_{n,k} v_n \bigr\} \; , \qquad \alpha_{n,k} = \omega_{n,k}
 *     Y_{n,k-1} - \sum_{ i = 1 }^{ k - 1 } \omega_{n,i} ( Y_{n,i} -
 *     Y_{n,i-1} ) \; ,
 * \f]
 * whose \f$ M_n \f$ pieces are rows of the PolyhedralFunction with
 * \f$ \alpha_j = \alpha_{n,k} \f$ and the single nonzero coefficient
 * \f$ - \omega_{n,k} \f$ in \f$ \beta_j \f$, that of the final volume of
 * \f$ n \f$. The range \f$ Y_{n,0} \leq v_n \leq Y_{n,M_n} \f$ is given by
 * the volume bounds of the HydroUnitBlock at the last instant, and the value
 * of the water at the beginning of the horizon is a constant. Since a
 * HydroSystemUnitBlock has a single PolyhedralFunction, the sum of such
 * functions over several reservoirs needs a row for each combination of their
 * pieces, i.e., \f$ \prod_n M_n \f$ rows (since the maximum of a sum of
 * independent maxima is the maximum over all the combinations).
 *
 * Multistage models. In a multistage model (see SDDPBlock) the UCBlock is the
 * problem of one stage, and \f$ \check\nu \f$ is the cutting-plane model of
 * the expected cost of the following stages as a function of the final
 * volumes. The initial volumes \f$ V^0 \f$ of the next stage are those final
 * volumes, which enter only the right-hand sides of its water balances (12)
 * of instant 0 (see HydroUnitBlock::generate_abstract_constraints() and
 * HydroUnitBlock::set_initial_volume()). Instead, the other initial
 * conditions of the valleys, i.e., the initial flows of the arcs, are data of
 * the stage and not part of the state. Hence, the cuts come from the dual
 * values of those rows. Let \f$ \mathcal{V}( V^0 ) \f$ be the optimal value
 * of the (convex) problem of the next stage for one realization of its random
 * data, as a function of its initial volumes, and \f$ y_n \f$ the dual value
 * of the water balance of instant 0 of reservoir \f$ n \f$ at the initial
 * volumes \f$ \hat V^0 \f$. Since the dual value is minus the derivative of
 * the optimal value with respect to the right-hand side (see
 * \ref ucbm_dual_sign), \f$ - y \f$ is a subgradient of \f$ \mathcal{V} \f$
 * at \f$ \hat V^0 \f$, i.e.,
 * \f[
 *   \mathcal{V}( v ) \geq \mathcal{V}( \hat V^0 ) - y^\top ( v - \hat V^0 )
 *     \qquad \text{for every } v \; .
 * \f]
 * The expectation of these affine minorants over the realizations, i.e.,
 * \f$ \alpha_j + \beta_j^\top v^{f} \f$ with
 * \f$ \beta_j = - \mathbb{E}[ y ] \f$ and \f$ \alpha_j = \mathbb{E}[
 * \mathcal{V}( \hat V^0 ) ] - \beta_j^\top \hat V^0 \f$, is then a new row of
 * \f$ \check\nu \f$. By linear programming duality, \f$ \alpha_j \f$ can
 * equivalently be computed from the dual values of the other rows of the next
 * stage and their right-hand sides, where the rows of the future cost of that
 * stage contribute the convex combination \f$ \sum_i y^{\phi}_i \alpha_i \f$
 * of their constants and \f$ y^{\phi} \f$ are the dual values of its epigraph
 * rows. Note that the cut is valid only if the problem of the stage is
 * convex, i.e., if it is solved as a continuous relaxation or through its
 * Lagrangian dual. In the latter case, the dual values of the water balances
 * are those of the last solution of the subproblem of this
 * HydroSystemUnitBlock, and they give a subgradient only if that solution is
 * at optimal multipliers and its Solver is a CDASolver (see \ref ucbm_multi
 * for the argument). We do not model the aggregation of the volumes of
 * several reservoirs into fewer state variables; one can still give a cut
 * computed on aggregated volumes \f$ \mathfrak{A} v^{f} \f$ (where
 * \f$ \mathfrak{A} \f$ is a matrix with fewer rows than columns) as the cut
 * \f$ \alpha_j + ( \mathfrak{A}^\top \beta_j )^\top v^{f} \f$ on the
 * individual volumes.
 *
 * Features not modeled. The future cost is part of the cost of this
 * UnitBlock: when the UCBlock is decomposed by a Lagrangian relaxation of its
 * linking constraints, it stays inside the subproblem of the
 * HydroSystemUnitBlock, which comprises all its HydroUnitBlock, and no
 * multiplier is associated with its rows. Hence, the future cost is separable
 * across the HydroSystemUnitBlock of a UCBlock, and the reservoirs whose
 * final volumes are valued jointly must belong to the same
 * HydroSystemUnitBlock. Neither the relaxation of the epigraph of
 * \f$ \check\nu \f$ (with a multiplier per cut) nor that of an equality
 * between the (aggregated) final volumes and copies of them that only
 * \f$ \check\nu \f$ uses, which would let the reservoirs be in different
 * subproblems, is available. A future value of the energy left in other
 * storages (batteries, load curtailment contracts) is not modeled either. */

class HydroSystemUnitBlock : public UnitBlock
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*---------------------- PUBLIC TYPES OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
 *
 * HydroSystemUnitBlock defines the following main public type:
 *
 * @{ */

/** @} ---------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 * @{ */

 /// constructor, takes the father
 /** Constructor of HydroSystemUnitBlock, taking possibly a pointer of its
  * father Block. */

 explicit HydroSystemUnitBlock( Block * father_block = nullptr )
  : UnitBlock( father_block ) , f_number_hydro_units( 0 ) {}

/*--------------------------------------------------------------------------*/
 /// destructor of HydroSystemUnitBlock

 virtual ~HydroSystemUnitBlock() override;

/** @} ---------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 * @{ */

/// extends Block::deserialize( netCDF::NcGroup )
/** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
 * the HydroSystemUnitBlock. Besides the mandatory "type" attribute of any
 * :Block, the group should contain the following:
 *
 * - The dimension "TimeHorizon", as for any UnitBlock.
 *
 * - The dimension "NumberHydroUnits" containing the number of hydro units
 *   (HydroUnitBlock) in the problem.
 *
 * - The groups "HydroUnitBlock_0", "HydroUnitBlock_1", ...,
 *   "HydroUnitBlock_(n-1)", with n == NumberHydroUnits, containing each one
 *   HydroUnitBlock.
 *
 * - The group "PolyhedralFunctionBlock", which contains a
 *   PolyhedralFunctionBlock (or a Block of a class derived from it, as given
 *   by its "type" attribute) whose PolyhedralFunction is the future cost
 *   \f$ \check\nu \f$ of the water left at the end of the horizon in all the
 *   reservoirs of all the HydroUnitBlock (see the description of the
 *   class). The group is optional: without it there is no future cost, and
 *   get_polyhedral_function_block() returns nullptr. The PolyhedralFunction
 *   must be convex (no dimension "PolyFunction_sign" of size 0), otherwise
 *   std::invalid_argument is thrown.
 *
 * The active Variable of the PolyhedralFunction are the volumes of all the
 * reservoirs at the last instant \f$ T - 1 \f$, and the columns of its
 * rows ("PolyFunction_A") refer to them in this order: the HydroUnitBlock
 * in the order of their groups and, within each HydroUnitBlock, its
 * reservoirs in their order. That is, with \f$ R_i \f$ the number of
 * reservoirs of HydroUnitBlock_i (see HydroUnitBlock::get_number_reservoirs())
 * and X[ 0 ], X[ 1 ], ... the active Variable of the PolyhedralFunction
 * (X[ k ] being the one returned by its get_active_var( k )):
 *
 * - X[ 0 ], ..., X[ \f$ R_0 - 1 \f$ ] are the final volumes of the
 *   reservoirs 0, ..., \f$ R_0 - 1 \f$ of HydroUnitBlock_0;
 *
 * - X[ \f$ R_0 \f$ ], ..., X[ \f$ R_0 + R_1 - 1 \f$ ] are the final volumes
 *   of the reservoirs 0, ..., \f$ R_1 - 1 \f$ of HydroUnitBlock_1;
 *
 * - and so on, the number of active Variable being \f$ \sum_i R_i \f$
 *   ("PolyFunction_NumVar"), which is NumberHydroUnits when every
 *   HydroUnitBlock has a single reservoir.
 *
 * See PolyhedralFunction::deserialize() for the format of the data in the
 * group "PolyhedralFunctionBlock". */

 void deserialize( const netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG
 /// extends UnitBlock::expected_dims()

 std::vector< std::string > expected_dims( void ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /* extends UnitBlock::expected_vars()
  * not necessary, no new variables

 std::vector< std::string > expected_vars( void ) const override;
 */
#endif

/*--------------------------------------------------------------------------*/
 /// generate the abstract variables of the HydroSystemUnitBlock
 /** Generates the Variable of all the sub-Block, and makes the final volumes
  * of the reservoirs, in the order described in deserialize(), the active
  * Variable of the PolyhedralFunction, if any. The HydroSystemUnitBlock has
  * no Variable of its own. */

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// generate the objective of the HydroSystemUnitBlock
 /** Method that generates the objective of the HydroSystemUnitBlock.
  *
  * - Objective function: the objective function of the HydroSystemUnitBlock
  *   is "empty" (a FRealObjective with a LinearFunction inside with no
  *   active variables), the cost of the HydroSystemUnitBlock being the sum
  *   of the Objective of its sub-Block, i.e., of the costs of the
  *   HydroUnitBlock and of the future cost \f$ \check\nu \f$ of the water
  *   (see the description of the class); the objectives of the sub-Block
  *   are generated as well. */

 void generate_objective( Configuration * objc = nullptr ) override;

/** @} ---------------------------------------------------------------------*/
/*-------- METHODS FOR READING THE DATA OF THE HydroSystemUnitBlock --------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the HydroSystemUnitBlock
 * @{ */

 /// returns the number of hydro units of the problem
 Index get_number_hydro_units( void ) const {
  return( f_number_hydro_units );
  }

/*--------------------------------------------------------------------------*/
 /// returns the i-th HydroUnitBlock, 0 <= i < get_number_hydro_units()
 /** Returns the i-th HydroUnitBlock; the index is checked only if NDEBUG is
  * not defined, in which case std::invalid_argument is thrown when
  * i >= get_number_hydro_units(). */

 HydroUnitBlock * get_hydro_unit_block( Index i ) const;

/*--------------------------------------------------------------------------*/
 /// returns the PolyhedralFunctionBlock, nullptr if there is none

 PolyhedralFunctionBlock * get_polyhedral_function_block( void ) const {
  if( v_Block.size() <= f_number_hydro_units )
   return( nullptr );
  return( static_cast< PolyhedralFunctionBlock * >( v_Block.back() ) );
  }

/*--------------------------------------------------------------------------*/
 /// the storages are the reservoirs of all the HydroUnitBlock, in order

 Index get_number_storages( void ) const override;

/*--------------------------------------------------------------------------*/
 /// returns the volumes of the given reservoir [see get_number_storages()]

 ColVariable * get_storage_level( Index storage ) override;

/*--------------------------------------------------------------------------*/
 /// returns the vector of active power variables of each HydroUnitBlock

 ColVariable * get_active_power( Index generator ) override;

/*--------------------------------------------------------------------------*/
 /// returns the vector of reactive power variables of each HydroUnitBlock

 ColVariable * get_reactive_power( Index generator ) override;

/*--------------------------------------------------------------------------*/
 /// returns the vector of primary reserve variables of each HydroUnitBlock

 ColVariable * get_primary_spinning_reserve( Index generator ) override;

/*--------------------------------------------------------------------------*/
 /// true if any of the HydroUnitBlock provides the primary reserve
 /** The generators of a HydroSystemUnitBlock come from different
  * HydroUnitBlock, which need not agree on providing the reserve: the
  * system has it if any of them has it, and the generators of those that do
  * not are skipped when a Solution is read and written. */

 bool has_primary_reserve( void ) const override;

/*--------------------------------------------------------------------------*/
 /// returns the vector of secondary reserve variables of each HydroUnitBlock

 ColVariable * get_secondary_spinning_reserve( Index generator ) override;

/*--------------------------------------------------------------------------*/
 /// true if any of the HydroUnitBlock provides the secondary reserve

 bool has_secondary_reserve( void ) const override;

/*--------------------------------------------------------------------------*/
 /// the generators are those of the HydroUnitBlock, unit after unit

 Index get_number_generators( void ) const override {
  return( v_gen_map.size() );
  }

/*--------------------------------------------------------------------------*/
 /// returns the inertia power of the generator in its HydroUnitBlock

 const double * get_inertia_power( Index generator ) const override;

/*--------------------------------------------------------------------------*/
 /// returns the minimum power of the generator in its HydroUnitBlock

 double get_min_power( Index t , Index generator = 0 ) const override;

/*--------------------------------------------------------------------------*/
 /// returns the maximum power of the generator in its HydroUnitBlock

 double get_max_power( Index t , Index generator = 0 ) const override;

/*--------------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/
 /// returns a Solution storing the current one of this HydroSystemUnitBlock
 /** This method must construct and return a (pointer to a) Solution object
  * representing the current "solution state" of this HydroSystemUnitBlock.
  * This is a HydroSystemUnitBlockSolution extending UnitBlockSolution with
  * the specific extra solution information of HydroUnitBlock.
  *
  * This Solution object is basically a UnitBlockSolution containing in
  * addition the collection of HydroUnitBlockSolution, one for each of the
  * inner HydroUnitBlock of this HydroSystemUnitBlock. However, since the
  * "root" UnitBlockSolution already contains the active power and other
  * variables, the "inner" HydroUnitBlockSolution are configured not to
  * store the same information.
  *
  * The parameter for deciding which kind of Solution must be returned is a
  * single int value, coded bitwise:
  *
  * - the first seven bits (bit 0 to bit 6) are passed to the "inner"
  *   HydroUnitBlockSolution with the first four bits masked (set to 0)
  *   so as to avoid duplicating information
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
  * - otherwise, it is 63 (save everything). */

 Solution * get_Solution( Configuration * solc = nullptr ,
			  bool emptys = true ) override;

/*--------------------------------------------------------------------------*/
 /// checks whether the current solution is feasible for the hydro system
 /** The solution is feasible if every sub-Block is. The tolerance and the
  * type of violation are taken from \p fsbc if it is a
  * SimpleConfiguration< double > (tolerance, relative violation) or a
  * SimpleConfiguration< std::pair< double , int > > (tolerance, relative
  * violation if the second is nonzero), otherwise from
  * f_BlockConfig->f_is_feasible_Configuration in the same way, otherwise
  * they are 0 and the relative violation; each sub-Block is checked with
  * them, unless its BlockConfig has its own f_is_feasible_Configuration. */

 bool is_feasible( bool useabstract = false ,
                   Configuration * fsbc = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// return the "appropriate" [HydroSystem]UnitBlockSolution

 UnitBlockSolution * new_Solution( void ) const override;

/** @} ---------------------------------------------------------------------*/
/*--------------- METHODS FOR SAVING THE HydroSystemUnitBlock --------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for printing & saving the HydroSystemUnitBlock
 * @{ */

/// extends Block::serialize( netCDF::NcGroup )
/** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
 * HydroSystemUnitBlock. See HydroSystemUnitBlock::deserialize( netCDF::
 * NcGroup ) for details of the format of the created netCDF group. */

 void serialize( netCDF::NcGroup & group ) const override;

/** @} ---------------------------------------------------------------------*/
/*----------- METHODS FOR MODIFYING THE HydroSystemUnitBlock ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the HydroSystemUnitBlock
 * @{ */

 /// sets which reserve variables should be there
 /** Basically just calls the method in all the sub-Block. */

 void set_reserve_vars( unsigned char what = 0 ) override {
  UnitBlock::set_reserve_vars( what );
  for( auto * b : v_Block )
   if( auto ub = dynamic_cast< HydroUnitBlock * >( b ) )
    ub->set_reserve_vars( what );
  }

/*--------------------------------------------------------------------------*/
 /// sets whether or not reactive power variables should be there
 /** Basically just calls the method in all the sub-Block. */

 void set_reactive_power( bool reactive = false ) override {
  UnitBlock::set_reactive_power( reactive );
  for( auto * b : v_Block )
   if( auto ub = dynamic_cast< HydroUnitBlock * >( b ) )
    ub->set_reactive_power( reactive );
  }

/** @} ---------------------------------------------------------------------*/
/*------------ METHODS FOR INITIALIZING THE HydroSystemUnitBlock -----------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the data of the HydroSystemUnitBlock
 * @{ */

 /// loading from a stream is not implemented: it throws

 void load( std::istream & input , char frmt = 0 ) override {
  throw( std::logic_error(
		    "HydroSystemUnitBlock::load() not implemented  yet" ) );
  }

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*--------------------- PROTECTED METHODS OF THE CLASS ---------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

/*---------------------------------- data ----------------------------------*/

 Index f_number_hydro_units;  ///< the number of hydro units of the problem

/*-------------------------------- variables -------------------------------*/

/*------------------------------- constraints ------------------------------*/

 FRealObjective objective;  ///< the objective function

/*----------------------------- data structures ----------------------------*/

 std::vector< std::pair< Index , Index > > v_gen_map;
 ///< the map between total generators and those in the inner HydroUnitBlock
 /**< v_gen_map[ i ].first  = index of sub-Block in which generator i is
  *   v_gen_map[ i ].second = index of the generator in sub-Block
  *                           v_gen_map[ i ].first to which generator i
  *                           corresponds */

/*--------------------------------------------------------------------------*/
/*--------------------- PRIVATE PART OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*-------------------- PRIVATE FIELDS OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/
/*---------------------- PRIVATE METHODS OF THE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/

 /// deserialize the sub-Blocks of HydroSystemUnitBlock

 void deserialize_sub_blocks( const netCDF::NcGroup & group );

/*--------------------------------------------------------------------------*/
 /// deserialize the sub-Blocks of HydroSystemUnitBlock that have the given
 /// prefix name

 void deserialize_sub_blocks( const netCDF::NcGroup & group ,
                              const std::string & sub_group_name_prefix ,
                              Index num_sub_blocks );

/*--------------------------------------------------------------------------*/
 /// deserialize the PolyhedralFunctionBlock

 void deserialize_polyhedral_function_block( const netCDF::NcGroup & group ,
					const std::string & sub_group_name );

/*--------------------------------------------------------------------------*/
 /// compute the total number of reservoirs

 Index get_total_number_reservoirs( void ) const;

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

 };  // end( class( HydroSystemUnitBlock ) )

/*--------------------------------------------------------------------------*/
/*------------------ CLASS HydroSystemUnitBlockSolution --------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// a [UnitBlock]Solution of a HydroSystemUnitBlock
/** The HydroSystemUnitBlockSolution class derives from UnitBlockSolution and
 * adds the "standard" information stored in there (active power, possibly
 * commitment and primary/secondary reserve) the other information that is
 * typical of the HydroSystemUnitBlock, i.e.,
 *
 * - one HydroUnitBlockSolution for each of the "inner" HydroUnitBlock of
 *   the HydroSystemUnitBlock
 *
 * Note, however, that since the "root" UnitBlockSolution contains (if so
 * required) all the "standard" information, the same is not duplicated in
 * the "inner" HydroUnitBlockSolution */

class HydroSystemUnitBlockSolution : public UnitBlockSolution
{

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*------------------------------- FRIENDS ----------------------------------*/

 friend HydroSystemUnitBlock;  ///< make HydroSystemUnitBlock friend

/*------ CONSTRUCTING AND DESTRUCTING HydroSystemUnitBlockSolution ---------*/

 /// constructor, it has nothing to do
 explicit HydroSystemUnitBlockSolution( void ) {}

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 void deserialize( const netCDF::NcGroup & group ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 ~HydroSystemUnitBlockSolution() override = default;
 ///< destructor: it is virtual, and empty

/*--- METHODS DESCRIBING THE BEHAVIOR OF A HydroSystemUnitBlockSolution ---*/

 void read( const Block * block ) override final;

 void write( Block * block ) override final;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// serialize a HydroSystemUnitBlockSolution into a netCDF::NcGroup
 /** Serialize a HydroSystemUnitBlockSolution into a netCDF::NcGroup. The
  * format is the one of UnitBlockSolution
  * [cf. UnitBlockSolution::serialize()], plus:
  *
  * - The dimension "NumberHydroUnits" containing the number of hydro units
  *   (HydroUnitBlock) in the problem and therefore the number of "inner"
  *   HydrUnitBlockSolution
  *
  * - The groups "HydroUnitBlockSolution_0", "HydroUnitBlockSolution_1", ...,
  *   "HydroUnitBlockSolution_(n-1)", with n == NumberHydroUnits, containing
  *   each one HydroUnitBlock.
  *
  * Note, however, that since the "root" UnitBlockSolution contains (if so
  * required) all the "standard" information, the same is not duplicated in
  * the "inner" HydroUnitBlockSolution */

 void serialize( netCDF::NcGroup & group ) const override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

 HydroSystemUnitBlockSolution * scale( double factor ) const override;

 void sum( const Solution * solution , double multiplier ) override;

 HydroSystemUnitBlockSolution * clone( bool empty = false ) const override;

/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/

 protected:

/*-------------------------- PROTECTED METHODS -----------------------------*/

 void print( std::ostream &output ) const override {
  output << "HydroSystemUnitBlockSolution [" << this << "]: " << std::endl;
  }

/*---------------------- PRIVATE PART OF THE CLASS -------------------------*/

 private:

/*---------------------------- PRIVATE FIELDS ------------------------------*/

 std::vector< HydroUnitBlockSolution * > v_innerSol;
 ///< v_innerSol[ i ] = HydroUnitBlockSolution of inner Block i

/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/

 };  // end( class( HydroSystemUnitBlockSolution ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* __HydroSystemUnitBlock */

/*--------------------------------------------------------------------------*/
/*--------------------- End File HydroSystemUnitBlock.h --------------------*/
/*--------------------------------------------------------------------------*/
