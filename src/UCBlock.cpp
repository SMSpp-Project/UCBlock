/*--------------------------------------------------------------------------*/
/*----------------------- File UCBlock.cpp ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the UCBlock class.
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
 * \author Kostas Tavlaridis-Gyparakis \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Ali Ghezelsoflu,
 *                      Rafael Durbano Lobato, Kostas Tavlaridis-Gyparakis,
 *                      Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "BlockInspection.h"

#include "AbstractBlock.h"

#include "BlockSolverConfig.h"

#include "FRealObjective.h"

#include "LinearFunction.h"

#include "Objective.h"

#include <array>
#include <functional>
#include <map>
#include <memory>
#include <numeric>

#include "UCBlock.h"

#include "BatteryUnitBlock.h"

#include "ACNetworkBlock.h"

#include "DCNetworkBlock.h"

#include "DesignNetworkBlock.h"

#include "HydroUnitBlock.h"

#include "IntermittentUnitBlock.h"

#include "ThermalUnitBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register UCBlock to the Block factory
SMSpp_insert_in_factory_cpp_1( UCBlock );

// register UCBlockSolution to the Solution factory
SMSpp_insert_in_factory_cpp_0( UCBlockSolution );

/*--------------------------------------------------------------------------*/
/*--------------------------- STATIC FUNCTIONS -----------------------------*/
/*--------------------------------------------------------------------------*/

template< class T , std::size_t K >
static void copy_multi_array( boost::multi_array< T , K > & to ,
			      const boost::multi_array< T , K > & from )
{
 std::vector< size_t > extent;
 auto shape = from.shape();
 extent.assign( shape , shape + from.num_dimensions() );
 to.resize( extent );
 to = from;
 }

/*--------------------------------------------------------------------------*/

static inline UnitBlock * UB( Block * b ) {
 return( static_cast< UnitBlock * >( b ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

static inline NetworkBlock * NB( Block * b ) {
 return( static_cast< NetworkBlock * >( b ) );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

static inline LinearFunction * LF( Function * f ) {
 return( static_cast< LinearFunction * >( f ) );
 }

/*--------------------------------------------------------------------------*/
/* The method that the class of a sub-Block registers in the methods factory
 * as "<classname>" + suffix, for the methods of UCBlock that hand a size to
 * the units or to the NetworkBlock [see resize_unit()]. It is looked up again
 * only when the class changes from one sub-Block to the next, so that a run
 * of sub-Block of the same class costs one lookup. */

template< class F >
class ClassMethod {
 public:

 ClassMethod( const char * caller , const char * suffix ,
              const char * form = "" )
  : f_caller( caller ) , f_suffix( suffix ) , f_form( form ) {}

 /// the method of the class of block, which the message calls the kind i
 F & operator()( const Block * block , const char * kind , Block::Index i ) {
  if( ( ! f_method ) || ( block->classname() != f_class ) ) {
   f_class = block->classname();
   f_method = Block::get_method< F >( f_class + f_suffix );
   if( ! f_method )
    throw( std::logic_error( std::string( "UCBlock::" ) + f_caller + ": " +
                             kind + " " + std::to_string( i ) + " is a " +
                             f_class + ", which does not register '" +
                             f_class + f_suffix + "'" + f_form ) );
   }
  return( *f_method );
  }

 private:

 const char * f_caller;
 const char * f_suffix;
 const char * f_form;  ///< the signature looked for, if not the plain one
 std::string f_class;
 F * f_method = nullptr;
 };

/*--------------------------------------------------------------------------*/
/* Hands each unit in units the value in the same position of values, through
 * the method that its class registers as "<classname>" + suffix, with the
 * one scalar a unit has [see UCBlock::resize_unit()]. */

static void size_units( UCBlock & ucb , const char * caller ,
                        const char * suffix , Block::MF_dbl_sp values ,
                        Block::c_Subset & units , ModParam issuePMod ,
                        ModParam issueAMod )
{
 if( units.size() > values.size() )
  throw( std::invalid_argument( std::string( "UCBlock::" ) + caller +
                                ": the span is shorter than the units it is "
                                "asked to write" ) );

 ClassMethod< Block::FunctionType< double > > method( caller , suffix );

 for( Block::Index k = 0 ; k < units.size() ; ++k ) {
  const auto u = units[ k ];
  if( u >= ucb.get_number_units() )
   throw( std::invalid_argument( std::string( "UCBlock::" ) + caller +
                                 ": there is no unit of index " +
                                 std::to_string( u ) ) );
  const auto unit = ucb.get_unit_block( u );
  std::invoke( method( unit , "unit" , u ) , unit , values[ k ] ,
               issuePMod , issueAMod );
  }
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/

static Block::Subset units_in( const UCBlock & ucb , Block::Range rng )
{
 rng.second = std::min( rng.second , ucb.get_number_units() );
 Block::Subset units;
 for( auto u = rng.first ; u < rng.second ; ++u )
  units.push_back( u );
 return( units );
 }

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/* The NetworkBlock holding the bounds of the node injections of the
 * interval i of nb, and the interval in it: a DesignNetworkBlock passes them
 * on to its sub-network of that interval, as interval 0 [see
 * DesignNetworkBlock::set_min_node_injection()]; nullptr if it has none. */

using BoundsHolder = std::pair< const NetworkBlock * , Block::Index >;

static BoundsHolder bounds_holder( const NetworkBlock * nb , Block::Index i )
{
 if( ! dynamic_cast< const DesignNetworkBlock * >( nb ) )
  return( BoundsHolder( nb , i ) );
 if( i >= nb->get_number_nested_Blocks() )
  return( BoundsHolder( nullptr , 0 ) );
 return( BoundsHolder(
          static_cast< const NetworkBlock * >( nb->get_nested_Block( i ) ) ,
          0 ) );
 }

/*--------------------------------------------------------------------------*/
/*--------------------------- METHODS OF UCBlock ---------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/

UCBlock::~UCBlock()
{
 Constraint::clear( v_node_injection_Const );
 Constraint::clear( v_reactive_node_injection_Const );
 Constraint::clear( v_PrimaryDemand_Const );
 Constraint::clear( v_SecondaryDemand_Const );
 Constraint::clear( v_InertiaDemand_Const );
 Constraint::clear( v_PollutantBudget_Const );
 for( auto & rows : v_carrier_balance_Const )
  Constraint::clear( rows );

 // the NetworkBlock of the carriers are among the sub-Block
 for( auto & block : v_Block )
  delete( block );
 v_Block.clear();

 delete( f_NetworkData );
 for( auto & cd : v_carrier )
  delete( cd.f_NetworkData );

 objective.clear();
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void UCBlock::deserialize_sub_blocks( const netCDF::NcGroup & group ,
                                      const std::string & prefix ,
                                      Index num_sub_blocks )
{
 auto sz = v_Block.size();
 v_Block.resize( sz + num_sub_blocks );
 // the sub-groups are listed once: NcGroup::getGroup() lists them all at
 // each call
 const auto sub_groups = group.getGroups();
 for( Index i = 0 ; i < num_sub_blocks ; ++i ) {
  std::string sub_group_name = prefix + std::to_string( i );
  auto it = sub_groups.find( sub_group_name );
  if( it == sub_groups.end() )
   throw( std::invalid_argument( "UCBlock::deserialize: " +
                                 sub_group_name + " not present" ) );

  if( auto bk = new_Block( it->second , this ) )
   v_Block[ sz++ ] = bk;
  else
   throw( std::invalid_argument( "UCBlock::deserialize: " +
                                 sub_group_name + " deserialize failed" ) );
  }
 }

/*--------------------------------------------------------------------------*/

void UCBlock::deserialize_network_blocks( const netCDF::NcGroup & group )
{
 Index cntr = 0;
 v_network_blocks.resize( f_number_networks , nullptr );

 // the sub-groups are listed once: NcGroup::getGroup() lists them all at
 // each call
 const auto sub_groups = group.getGroups();
 for( Index i = 0 ; i < f_number_networks ; ++i ) {
  std::string sub_group_name = "NetworkBlock_" + std::to_string( i );
  auto it = sub_groups.find( sub_group_name );
  if( it == sub_groups.end() )
   continue;

  auto res = almost_new_Block( it->second , this );
  if( auto nbi = dynamic_cast< NetworkBlock * >( res.second ) ) {
   v_network_blocks[ i ] = nbi;
   nbi->set_NetworkData( f_NetworkData );  // set the global NetworkData
   nbi->deserialize( res.first );          // deserialize the NetworkBlock
   ++cntr;
   }
  else {
   delete( nbi );
   throw( std::invalid_argument( "UCBlock::deserialize:" + sub_group_name +
                                 " not a valid NetworkBlock" ) );
  }
 }

 if( cntr ) {
  v_Block.resize( f_number_units + f_number_networks );
  std::copy( v_network_blocks.begin() , v_network_blocks.end() ,
             std::next( v_Block.begin() , f_number_units ) );
  }
 else
  v_network_blocks.clear();
 }

/*--------------------------------------------------------------------------*/

void UCBlock::deserialize( const netCDF::NcGroup & group )
{
 // mandatory dimensions- - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 deserialize_dim( group , "TimeHorizon" , f_time_horizon , false );
 deserialize_dim( group , "NumberUnits" , f_number_units , false );

 // optional dimensions - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! deserialize_dim( group , "NumberNetworks" , f_number_networks ) )
  f_number_networks = f_time_horizon;

 // optional: the "ConstantTerm" of a "NetworkBlock_n" group, if any,
 // prevails over NetworkConstantTerms[ n ]
 ::deserialize( group , "NetworkConstantTerms" , f_number_networks ,
                v_network_constant_terms );

 // true if "NetworkBlock_n" is given with its own "ConstantTerm"; the
 // sub-groups are listed once, NcGroup::getGroup() lists them all at each
 // call
 const auto sub_groups = group.getGroups();
 auto own_constant = [ & sub_groups ]( Index n ) {
  auto it = sub_groups.find( "NetworkBlock_" + std::to_string( n ) );
  return( ( it != sub_groups.end() ) &&
          ( ! it->second.getVar( "ConstantTerm" ).isNull() ) );
  };

 if( ! ::deserialize( group , network_block_classname ,
                      "NetworkBlockClassname" ) )
  network_block_classname = "DCNetworkBlock";

 if( ! ::deserialize( group , network_data_classname ,
		      "NetworkDataClassname" ) )
  network_data_classname = ( network_block_classname == "ACNetworkBlock" )
                             ? "ACNetworkData" : "DCNetworkData";

 Index number_nodes = 1;  // default: no networks
 if( deserialize_dim( group , "NumberNodes" , number_nodes ) ) {
  // if "NumberNodes" is defined, then a "global" NetworkData is there

  f_NetworkData = NetworkBlock::NetworkData::new_NetworkData(
                                                    network_data_classname );
  if( ! f_NetworkData )
   throw( std::invalid_argument(
                  "UCBlock::deserialize: NetworkData missing in UCBlock" ) );
  f_NetworkData->deserialize( group );
  }

 f_number_primary_zones = 0;
 deserialize_dim( group , "NumberPrimaryZones" , f_number_primary_zones );

 f_number_secondary_zones = 0;
 deserialize_dim( group , "NumberSecondaryZones" , f_number_secondary_zones );

 f_number_inertia_zones = 0;
 deserialize_dim( group , "NumberInertiaZones" , f_number_inertia_zones );

 f_number_pollutants = 0;
 deserialize_dim( group , "NumberPollutants" , f_number_pollutants );

 // variables (all optional)- - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 ::deserialize( group , "ActivePowerDemand" ,
                { number_nodes , f_time_horizon } , v_active_power_demand );

 ::deserialize( group , "ReactivePowerDemand" ,
                { number_nodes , f_time_horizon } , v_reactive_power_demand );

 ::deserialize( group , "PrimaryZones" , number_nodes ,
                v_primary_zones , true , true );

 ::deserialize( group , "PrimaryDemand" ,
                { f_number_primary_zones , f_time_horizon } ,
                v_primary_demand , true , false );

 ::deserialize( group , "SecondaryZones" , number_nodes ,
                v_secondary_zones , true , true );

 ::deserialize( group , "SecondaryDemand" ,
                { f_number_secondary_zones , f_time_horizon } ,
                v_secondary_demand , true , false );

 ::deserialize( group , "InertiaZones" , number_nodes ,
                v_inertia_zones , true , true );

 ::deserialize( group , "InertiaDemand" ,
                { f_number_inertia_zones , f_time_horizon } ,
                v_inertia_demand , true , false );

 // without the vector all the nodes are in the only zone, which cannot be
 // if there are more
 if( ( ( f_number_primary_zones > 1 ) && v_primary_zones.empty() ) ||
     ( ( f_number_secondary_zones > 1 ) && v_secondary_zones.empty() ) ||
     ( ( f_number_inertia_zones > 1 ) && v_inertia_zones.empty() ) )
  throw( std::invalid_argument( "UCBlock::deserialize: PrimaryZones, "
                                "SecondaryZones and InertiaZones are "
                                "mandatory if there is more than one zone "
                                "of that kind" ) );

 // the pollutant zones; PollutantRho depends on the number of electrical
 // generators, which is only known after the UnitBlock are loaded, hence it
 // is read further down together with PollutantBudget
 v_number_pollutant_zones.clear();
 f_total_number_pollutant_zones = 0;
 if( f_number_pollutants ) {
  if( ! ::deserialize( group , "NumberPollutantZones" , f_number_pollutants ,
                       v_number_pollutant_zones , true , false ) )
   v_number_pollutant_zones.assign( f_number_pollutants , 1 );

  f_total_number_pollutant_zones = std::accumulate(
		             v_number_pollutant_zones.begin() ,
		             v_number_pollutant_zones.end() , Index( 0 ) );

  Index tnpz = 0;
  if( deserialize_dim( group , "TotalNumberPollutantZones" , tnpz ) &&
      ( tnpz != f_total_number_pollutant_zones ) )
   throw( std::invalid_argument( "UCBlock::deserialize: "
				 "TotalNumberPollutantZones is not the sum "
				 "of NumberPollutantZones" ) );

  if( ! ::deserialize( group , "PollutantZones" ,
                       { f_number_pollutants , number_nodes } ,
                       v_pollutant_zones , true , true ) )
   if( std::any_of( v_number_pollutant_zones.begin() ,
                    v_number_pollutant_zones.end() ,
                    []( Index nz ) { return( nz > 1 ); } ) )
    throw( std::invalid_argument( "UCBlock::deserialize: PollutantZones "
				  "mandatory if some pollutant has more "
				  "than one zone" ) );
  }
 else
  v_pollutant_zones.resize( boost::extents[ 0 ][ 0 ] );

 // the carriers other than the carrier 0, save their NetworkBlock, which
 // are built once the units are there- - - - - - - - - - - - - - - - - - - -
 deserialize_carriers( group );

 // reset all existing sub-Block, if any- - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 for( auto block : v_Block )
  delete( block );
 v_Block.clear();

 // load and "configure" the UnitBlock- - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 deserialize_sub_blocks( group , "UnitBlock_" , f_number_units );

 // check which reserve variables (primary, secondary, inertia) need be
 // defined and "inform" the UnitBlock about this

 unsigned char what = 0;
 if( f_number_primary_zones > 0 )  // primary spinning reserve variables
  what += 1;

 if( f_number_secondary_zones > 0 )  // secondary spinning reserve variables
  what += 2;

 if( f_number_inertia_zones > 0 )  // inertia reserve variables
  what += 4;

 if( what > 0 )
  for( auto * b : v_Block )
   if( auto ub = dynamic_cast< UnitBlock * >( b ) )
    ub->set_reserve_vars( what );

 // load or create the NetworkBlock (if necessary)- - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 deserialize_network_blocks( group );  // try to load the NetworkBlock

 if( v_network_blocks.empty() && v_active_power_demand.empty() )
  throw( std::invalid_argument( "UCBlock::deserialize: ActivePowerDemand "
                                "mandatory if no NetworkBlocks" ) );

 // if number_nodes == 1, NetworkBlocks are useless and therefore removed
 if( number_nodes == 1 ) {
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // curios case: there is only one node and therefore NetworkBlock are
  // useless, save that the active power demand may be stored in there
  // rather than in v_active_power_demand; move it there and delete the
  // useless NetworkBlock. reactive power, if a ReactivePowerDemand was
  // provided, is still handled here directly by the single-node reactive
  // balance (see below and generate_reactive_node_injection_constraints())

  if( v_active_power_demand.empty() ) {
   // if active power demand is not defined, do it now and preload it with
   // zeros in case some NetworkBlock is not there
   v_active_power_demand.resize( MAdouble_ext()[ 1 ][ f_time_horizon ] );
   for( auto apdit = v_active_power_demand.data() ;
        apdit != v_active_power_demand.data() + f_time_horizon ; )
	*( apdit++ ) = 0;
   }

  if( ! v_network_blocks.empty() ) {
   Index t = 0;
   for( Index n = 0 ; n < f_number_networks ; ++n )
    if( v_network_blocks[ n ] ) {
     for( Index i = 0 ; i < v_network_blocks[ n ]->get_number_intervals() ;
          ++i , ++t ) {
      if( auto ad = v_network_blocks[ n ]->get_active_demand( i ) )
       v_active_power_demand[ 0 ][ t ] = ad[ 0 ];
      }

     // its constant term, if any, goes in the Objective of the UCBlock
     // [see generate_objective()]
     if( own_constant( n ) ) {
      v_network_constant_terms.resize( f_number_networks , 0 );
      v_network_constant_terms[ n ] =
                                   v_network_blocks[ n ]->get_const_term();
      }

     delete v_network_blocks[ n ];
     }

   v_network_blocks.clear();
   v_Block.resize( f_number_units );
   }

  // single-bus reactive power: if a ReactivePowerDemand has been provided,
  // enable reactive power so that the UnitBlock create the reactive variables
  // and the (single-node) reactive balance constraints are generated below
  if( ! v_reactive_power_demand.empty() )
   f_has_reactive = true;
  }
 else {
  // number_nodes > 1 - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // if (some of) the NetworkBlock don't exist, create them now and feed them
  // the "global" NetworkData
  if( v_network_blocks.empty() ) {
   v_network_blocks.resize( f_number_networks );
   v_Block.resize( f_number_units + f_number_networks );
   }

  // meanwhile compute the number of covered intervals
  Index sum_intervals = 0;
  // Index t = 0;
  for( Index n = 0 ; n < f_number_networks ; ++n ) {
   auto nbi = v_network_blocks[ n ];
   if( ! nbi ) {  // NetworkBlock n does not exist: create a NetworkBlock
    nbi = dynamic_cast< NetworkBlock * >(
			      new_Block( network_block_classname , this ) );
    v_network_blocks[ n ] = nbi;
    v_Block[ f_number_units + n ] = nbi;
    }

   // tell it which time instant is its starting one
   nbi->set_time_instant( sum_intervals );
			 
   // if the NetworkBlock does not have its own NetworkData...
   if( ! nbi->get_NetworkData() ) {
    if( ! f_NetworkData )
     throw( std::invalid_argument( "UCBlock::deserialize: NetworkData "
                                   "missing in NetworkBlock_" +
                                   std::to_string( n ) + " and in UCBlock" )
	    );

    // ... then set the UCBlock global one
    nbi->set_NetworkData( f_NetworkData );
    // assert that the global NetworkData passed is of the right type
    assert( nbi->get_NetworkData() != nullptr );
    }

   sum_intervals += nbi->get_number_intervals();

   if( ! n )  // for the first one only, check if reactive power is there
    f_has_reactive = nbi->handles_reactive();

   // the constant term: that of the group of the NetworkBlock, if any,
   // otherwise NetworkConstantTerms[ n ], if given
   if( ( ! v_network_constant_terms.empty() ) && ( ! own_constant( n ) ) )
    nbi->set_constant_term( v_network_constant_terms[ n ] );
   }

  // the constant terms are in the NetworkBlock, which serialize() writes
  v_network_constant_terms.clear();

  // it must be that sum_intervals == f_time_horizon, i.e., the NetworkBlock
  // we receive in input / create span by the whole time horizon of the
  // problem, throw exception otherwise
  if( sum_intervals != f_time_horizon )
   throw( std::invalid_argument( "UCBlock::deserialize: the intervals "
				 "spanned by the NetworkBlock do not cover "
				 "the time horizon" ) );

  // now pass the active power demand, if any
  if( v_active_power_demand.empty() ) {
   // if not, check that the demand was there in the NetworkBlock
   for( auto * nbi : v_network_blocks )
    for( Index i = 0 ; i < nbi->get_number_intervals() ; ++i )
     if( ! nbi->get_active_demand( i ) )
      throw( std::invalid_argument( "UCBlock::deserialize: ActivePowerDemand "
                                    "missing in UCBlock and in NetworkBlock"
                                    ) ); 
   }
  else {
   Index t = 0;
   boost::multi_array< double , 2 > ap_v;
   using range = boost::multi_array_types::index_range;
   for( auto * nbi : v_network_blocks ) {
    auto ni = nbi->get_number_intervals();
    ap_v.resize( boost::extents[ ni ][ number_nodes ] );
    for( Index i = 0 ; i < ni ; ++i , ++t ) {
     auto ap_c = v_active_power_demand[
		          boost::indices[ range( 0 , number_nodes ) ][ t ] ];
     std::copy( ap_c.begin() , ap_c.end() , ap_v[ i ].begin() );
     }

    nbi->set_ActiveDemand( ap_v );
    }

   // v_active_power_demand used up, disband it
   v_active_power_demand.resize( MAdouble_ext()[ 0 ][ 0 ] );
   }

  // if reactive power is considered, do the same for it
  if( f_has_reactive ) {
   if( v_reactive_power_demand.empty() ) {
    // demand is not defined, check it was there in the NetworkBlock
    for( auto * nbi : v_network_blocks )
     for( Index i = 0 ; i < nbi->get_number_intervals() ; ++i )
      if( ! nbi->get_reactive_demand( i ) )
       throw( std::invalid_argument(
			       "UCBlock::deserialize: ReactivePowerDemand "
                               "missing in UCBlock and in NetworkBlock" ) ); 
    }
   else {
    Index t = 0;
    boost::multi_array< double , 2 > rp_v;
    using range = boost::multi_array_types::index_range;
    for( auto * nbi : v_network_blocks ) {
     auto ni = nbi->get_number_intervals();
     rp_v.resize( boost::extents[ ni ][ number_nodes ] );
     for( Index i = 0 ; i < ni ; ++i , ++t ) {
      auto rp_c = v_reactive_power_demand[
		          boost::indices[ range( 0 , number_nodes ) ][ t ] ];
      std::copy( rp_c.begin() , rp_c.end() , rp_v[ i ].begin() );
      }

     nbi->set_ReactiveDemand( rp_v );
     }

    // v_reactive_power_demand used up, disband it
    v_reactive_power_demand.resize( MAdouble_ext()[ 0 ][ 0 ] );
    }
   }  // end( if( f_has_reactive ) )
  }  // end( else( number_nodes > 1 ) )- - - - - - - - - - - - - - - - - - -

 // if reactive power is there, inform the UnitBlock- - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( f_has_reactive )
  for( Index u = 0 ; u < f_number_units ; ++u )
   UB( v_Block[ u ] )->set_reactive_power( true );

 // recover the generator --> node mapping- - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 Index total_generators = 0;
 for( Index i = 0 ; i < f_number_units ; ++i )
  total_generators += UB( v_Block[ i ] )->get_number_generators();

 // "NumberGenerators", or by its former name "NumberElectricalGenerators"
 if( ( ! deserialize_dim( group , "NumberGenerators" ,
                          f_number_elc_generators ) ) &&
     ( ! deserialize_dim( group , "NumberElectricalGenerators" ,
                          f_number_elc_generators ) ) )
  f_number_elc_generators = total_generators;
 else
  // a number that the units do not have is refused: everything indexed over
  // the generators, from GeneratorNode to the emission rates, would then be
  // read over the wrong length, and the rows this Block builds over them are
  // sized with it
  if( f_number_elc_generators != total_generators )
   throw( std::invalid_argument( "UCBlock::deserialize: "
	  "NumberGenerators is " +
	  std::to_string( f_number_elc_generators ) + " while the " +
	  std::to_string( f_number_units ) + " units have " +
	  std::to_string( total_generators ) + " generators" ) );

 deserialize_generator_nodes( group );

 // the unit of each generator and the carrier of each of its outputs, and
 // the NetworkBlock of the carriers other than the carrier 0- - - - - - - - -
 map_generators();
 build_carrier_networks();

 // read the pollutant budget and emission rates- - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 v_pollutant_budget.clear();
 v_pollutant_min_budget.clear();
 f_number_storages = 0;
 if( f_number_pollutants ) {
  ::deserialize( group , "PollutantBudget" , f_total_number_pollutant_zones ,
                 v_pollutant_budget , false , false );

  if( ! ::deserialize( group , "PollutantMinBudget" ,
                       f_total_number_pollutant_zones ,
                       v_pollutant_min_budget , true , false ) )
   v_pollutant_min_budget.assign( f_total_number_pollutant_zones ,
                                  -Inf< double >() );

  // the size is checked here, since the first dimension may be a singleton
  ::deserialize( group , "PollutantRho" , {} , v_pollutant_rho , false ,
                 false );

  auto shp = v_pollutant_rho.shape();
  if( ( ( shp[ 0 ] != 1 ) && ( shp[ 0 ] != f_time_horizon ) ) ||
      ( shp[ 1 ] != f_number_pollutants ) ||
      ( shp[ 2 ] != f_number_elc_generators ) )
   throw( std::invalid_argument( "UCBlock::deserialize: PollutantRho must "
				 "have size [ 1 or TimeHorizon , "
				 "NumberPollutants , "
				 "NumberElectricalGenerators ]" ) );

  // the storages of all the units, one after the other
  for( Index i = 0 ; i < f_number_units ; ++i )
   f_number_storages += UB( v_Block[ i ] )->get_number_storages();

  Index nst = 0;
  if( deserialize_dim( group , "NumberStorages" , nst ) &&
      ( nst != f_number_storages ) )
   throw( std::invalid_argument( "UCBlock::deserialize: NumberStorages is "
				 "not the number of storages of the units"
				 ) );

  if( ::deserialize( group , "PollutantStorageRho" , {} ,
                     v_pollutant_storage_rho , true , false ) ) {
   auto shs = v_pollutant_storage_rho.shape();
   if( ( ( shs[ 0 ] != 1 ) && ( shs[ 0 ] != f_time_horizon ) ) ||
       ( shs[ 1 ] != f_number_pollutants ) ||
       ( shs[ 2 ] != f_number_storages ) )
    throw( std::invalid_argument( "UCBlock::deserialize: PollutantStorageRho "
				  "must have size [ 1 or TimeHorizon , "
				  "NumberPollutants , NumberStorages ]" ) );
   }
  else
   v_pollutant_storage_rho.resize( boost::extents[ 0 ][ 0 ][ 0 ] );
  }
 else {
  v_pollutant_rho.resize( boost::extents[ 0 ][ 0 ][ 0 ] );
  v_pollutant_storage_rho.resize( boost::extents[ 0 ][ 0 ][ 0 ] );
  }

 // compute and store the min and max node injection into each NetworkBlock -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 set_node_injection_bounds();

 // finally call the method of the base class - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note that by calling this last we ensure that the NetworkData (if any) is
 // properly initialised so that its expected stuff can be counted when it
 // is checked there inside
 Block::deserialize( group );

 }  // end( UCBlock::deserialize )

/*--------------------------------------------------------------------------*/

void UCBlock::deserialize_carriers( const netCDF::NcGroup & group )
{
 // the NetworkBlock of the carriers are among the sub-Block, which are
 // deleted with them
 for( auto & cd : v_carrier )
  delete( cd.f_NetworkData );
 v_carrier.clear();
 v_carrier_name.clear();

 Index nc = 1;
 deserialize_dim( group , "NumberCarriers" , nc );
 if( ! nc )
  throw( std::invalid_argument( "UCBlock::deserialize: NumberCarriers is 0"
				) );

 if( ( ! ::deserialize( group , "CarrierName" , nc , v_carrier_name ) ) &&
     ( nc > 1 ) )
  throw( std::invalid_argument( "UCBlock::deserialize: CarrierName "
				"mandatory if NumberCarriers > 1" ) );

 // the names are distinct, and the electricity can only be the carrier 0
 for( Index c = 0 ; c < v_carrier_name.size() ; ++c ) {
  if( ( c > 0 ) && ( v_carrier_name[ c ] == UnitBlock::kElectricity ) )
   throw( std::invalid_argument( "UCBlock::deserialize: the carrier " +
				 v_carrier_name[ c ] + " can only be the "
				 "carrier 0" ) );
  for( Index d = 0 ; d < c ; ++d )
   if( v_carrier_name[ d ] == v_carrier_name[ c ] )
    throw( std::invalid_argument( "UCBlock::deserialize: the carrier " +
				  v_carrier_name[ c ] + " is there twice" ) );
  }

 if( nc == 1 )
  return;

 v_carrier.resize( nc );

 // the sub-groups are listed once: NcGroup::getGroup() lists them all at
 // each call
 const auto sub_groups = group.getGroups();
 for( Index c = 1 ; c < nc ; ++c ) {
  auto & cd = v_carrier[ c ];
  cd.network_block_classname = "DCNetworkBlock";
  cd.network_data_classname = "DCNetworkData";

  // without its group, a carrier has one node, no demand and the unique
  // zone of each pollutant
  auto it = sub_groups.find( "Carrier_" + std::to_string( c ) );
  if( it == sub_groups.end() ) {
   if( std::any_of( v_number_pollutant_zones.begin() ,
		    v_number_pollutant_zones.end() ,
		    []( Index nz ) { return( nz > 1 ); } ) )
    throw( std::invalid_argument( "UCBlock::deserialize: Carrier_" +
				  std::to_string( c ) + " mandatory, since "
				  "some pollutant has more than one zone" ) );
   continue;
   }
  const auto & cg = it->second;

  if( ! ::deserialize( cg , cd.network_block_classname ,
		       "NetworkBlockClassname" ) )
   cd.network_block_classname = "DCNetworkBlock";

  if( ! ::deserialize( cg , cd.network_data_classname ,
		       "NetworkDataClassname" ) )
   cd.network_data_classname =
                       ( cd.network_block_classname == "ACNetworkBlock" ) ?
                       "ACNetworkData" : "DCNetworkData";

  // the dimensions are those of the group itself, not of the UCBlock [see
  // deserialize_local_dim()]
  Index nn = 1;
  if( deserialize_local_dim( cg , "NumberNodes" , nn ) && ( nn > 1 ) ) {
   cd.f_NetworkData = NetworkBlock::NetworkData::new_NetworkData(
					           cd.network_data_classname );
   if( ! cd.f_NetworkData )
    throw( std::invalid_argument( "UCBlock::deserialize: " +
				  cd.network_data_classname + " of Carrier_" +
				  std::to_string( c ) + " not a NetworkData"
				  ) );
   cd.f_NetworkData->deserialize( cg );
   }

  ::deserialize( cg , "EnergyDemand" , { nn , f_time_horizon } ,
		 cd.v_demand );

  // its NetworkBlock, if given, as those of the carrier 0: "NumberNetworks"
  // of them (one per time instant if it is not there), "NetworkBlock_n"
  // being the n-th, which gets the NetworkData of the carrier; those that
  // are not given are built by build_carrier_networks()
  if( cd.f_NetworkData ) {
   Index nw = 0;
   const bool given_nw = deserialize_local_dim( cg , "NumberNetworks" , nw );
   if( given_nw && ( ! nw ) )
    throw( std::invalid_argument( "UCBlock::deserialize: NumberNetworks of "
				  "Carrier_" + std::to_string( c ) +
				  " is 0" ) );
   const auto cg_groups = cg.getGroups();
   const auto nb_group = [ & cg_groups ]( Index n ) {
    return( cg_groups.find( "NetworkBlock_" + std::to_string( n ) ) );
    };
   if( ! given_nw )
    nw = f_time_horizon;
   for( Index n = 0 ; n < nw ; ++n ) {
    auto jt = nb_group( n );
    if( jt == cg_groups.end() )
     continue;
    if( ! cd.f_given_networks ) {
     cd.f_given_networks = true;
     cd.v_network_blocks.assign( nw , nullptr );
     }
    auto res = almost_new_Block( jt->second , this );
    auto nb = dynamic_cast< NetworkBlock * >( res.second );
    if( ! nb ) {
     delete( res.second );
     throw( std::invalid_argument( "UCBlock::deserialize: Carrier_" +
				   std::to_string( c ) + "/NetworkBlock_" +
				   std::to_string( n ) +
				   " not a valid NetworkBlock" ) );
     }
    cd.v_network_blocks[ n ] = nb;
    nb->set_NetworkData( cd.f_NetworkData );
    nb->deserialize( res.first );
    }
   if( given_nw && ( ! cd.f_given_networks ) ) {
    cd.f_given_networks = true;
    cd.v_network_blocks.assign( nw , nullptr );
    }
   }

  if( f_number_pollutants &&
      ( ! ::deserialize( cg , "PollutantZones" , { f_number_pollutants ,
						     nn } ,
			 cd.v_pollutant_zones , true , true ) ) &&
      std::any_of( v_number_pollutant_zones.begin() ,
		   v_number_pollutant_zones.end() ,
		   []( Index nz ) { return( nz > 1 ); } ) )
   throw( std::invalid_argument( "UCBlock::deserialize: PollutantZones of "
				 "Carrier_" + std::to_string( c ) +
				 " mandatory if some pollutant has more than "
				 "one zone" ) );
  }
 }  // end( UCBlock::deserialize_carriers )

/*--------------------------------------------------------------------------*/

void UCBlock::build_carrier_networks( void )
{
 for( Index c = 1 ; c < get_number_carriers() ; ++c ) {
  auto & cd = v_carrier[ c ];
  cd.v_network_of_instant.clear();
  if( ! cd.f_NetworkData ) {
   cd.v_network_blocks.clear();
   continue;  // one node: the demand is in the balance (1c) itself
   }

  // one NetworkBlock per time instant, unless they are given [see
  // deserialize_carriers()], each with the NetworkData of the carrier and
  // the demand of the instants it covers
  if( ! cd.f_given_networks )
   cd.v_network_blocks.assign( f_time_horizon , nullptr );
  const auto nn = cd.f_NetworkData->get_number_nodes();
  for( Index n = 0 ; n < cd.v_network_blocks.size() ; ++n ) {
   auto nb = cd.v_network_blocks[ n ];
   if( ! nb ) {
    nb = dynamic_cast< NetworkBlock * >(
			     new_Block( cd.network_block_classname , this ) );
    if( ! nb )
     throw( std::invalid_argument( "UCBlock::deserialize: " +
				   cd.network_block_classname + " of Carrier_"
				   + std::to_string( c ) + " not a NetworkBlock"
				   ) );
    cd.v_network_blocks[ n ] = nb;
    nb->set_NetworkData( cd.f_NetworkData );
    }
   v_Block.push_back( nb );

   nb->set_time_instant( cd.v_network_of_instant.size() );
   // the reactive power is a matter of the electricity alone
   if( nb->handles_reactive() )
    throw( std::invalid_argument( "UCBlock::deserialize: the NetworkBlock "
				  "of Carrier_" + std::to_string( c ) +
				  " handles reactive power, which only the "
				  "carrier 0 can" ) );

   const auto ni = nb->get_number_intervals();
   boost::multi_array< double , 2 > demand( boost::extents[ ni ][ nn ] );
   for( Index i = 0 ; i < ni ; ++i ) {
    const auto t = cd.v_network_of_instant.size();
    for( Index node = 0 ; node < nn ; ++node )
     demand[ i ][ node ] = ( cd.v_demand.empty() || ( t >= f_time_horizon ) )
                           ? 0 : cd.v_demand[ node ][ t ];
    cd.v_network_of_instant.emplace_back( n , i );
    }
   nb->set_ActiveDemand( demand );
   }

  if( cd.v_network_of_instant.size() != f_time_horizon )
   throw( std::invalid_argument( "UCBlock::deserialize: the NetworkBlock of "
				 "Carrier_" + std::to_string( c ) + " cover " +
				 std::to_string( cd.v_network_of_instant.size()
						 ) + " time instants instead "
				 "of " + std::to_string( f_time_horizon ) ) );

  // the demand is in the NetworkBlock, as that of the carrier 0
  cd.v_demand.resize( boost::extents[ 0 ][ 0 ] );
  }
 }  // end( UCBlock::build_carrier_networks )

/*--------------------------------------------------------------------------*/

void UCBlock::deserialize_generator_nodes( const netCDF::NcGroup & group )
{
 const auto nc = get_number_carriers();
 v_generator_node.clear();
 v_generator_carrier_node.resize( boost::extents[ 0 ][ 0 ] );

 auto gn_var = group.getVar( "GeneratorNode" );
 if( gn_var.isNull() ) {
  // every generator is then at the only node of each carrier, which with
  // no generator at all (a network of demand only) says nothing
  for( Index c = 0 ; f_number_elc_generators && ( c < nc ) ; ++c )
   if( get_number_nodes( c ) > 1 )
    throw( std::invalid_argument( "UCBlock::deserialize: GeneratorNode "
				  "mandatory if a carrier has more than one "
				  "node" ) );
  return;
  }

 // read as signed, the -1 of a generator that has no node on a carrier not
 // fitting an unsigned type
 const auto to_node = []( int n ) {
  return( n < 0 ? Inf< Index >() : Index( n ) );
  };

 if( gn_var.getDimCount() <= 1 ) {
  // the nodes of the carrier 0: the other carriers must have one node, where
  // every generator is
  for( Index c = 1 ; c < nc ; ++c )
   if( get_number_nodes( c ) > 1 )
    throw( std::invalid_argument( "UCBlock::deserialize: GeneratorNode "
				  "must be indexed over the carriers too if "
				  "a carrier other than 0 has more than one "
				  "node" ) );

  std::vector< int > gn;
  ::deserialize( group , "GeneratorNode" , f_number_elc_generators , gn ,
		 false , true );
  v_generator_node.resize( gn.size() );
  std::transform( gn.begin() , gn.end() , v_generator_node.begin() ,
		  to_node );
  }
 else {
  boost::multi_array< int , 2 > gn;
  ::deserialize( group , "GeneratorNode" , { f_number_elc_generators , nc } ,
		 gn , false , false );
  v_generator_carrier_node.resize(
			  boost::extents[ f_number_elc_generators ][ nc ] );
  v_generator_node.resize( f_number_elc_generators );
  for( Index g = 0 ; g < f_number_elc_generators ; ++g ) {
   for( Index c = 0 ; c < nc ; ++c )
    v_generator_carrier_node[ g ][ c ] = to_node( gn[ g ][ c ] );
   v_generator_node[ g ] = v_generator_carrier_node[ g ][ 0 ];
   }
  }

 for( Index g = 0 ; g < f_number_elc_generators ; ++g )
  for( Index c = 0 ; c < nc ; ++c ) {
   const auto node = get_generator_node( g , c );
   if( ( node != Inf< Index >() ) && ( node >= get_number_nodes( c ) ) )
    throw( std::invalid_argument( "UCBlock::deserialize: GeneratorNode of "
				  "generator " + std::to_string( g ) +
				  " on carrier " + std::to_string( c ) +
				  " is not a node" ) );
   }
 }  // end( UCBlock::deserialize_generator_nodes )

/*--------------------------------------------------------------------------*/

void UCBlock::map_generators( void )
{
 v_generator_unit.clear();
 v_generator_unit.reserve( f_number_elc_generators );
 v_output_carrier_index.clear();
 v_output_carrier_index.reserve( f_number_elc_generators );
 v_first_generator.assign( 1 , 0 );
 v_first_storage.assign( 1 , 0 );

 for( Index u = 0 ; u < f_number_units ; ++u ) {
  const auto unit_block = get_unit_block( u );
  v_first_generator.push_back( v_first_generator.back() +
                               unit_block->get_number_generators() );
  v_first_storage.push_back( v_first_storage.back() +
                             unit_block->get_number_storages() );
  // the outputs given per generator, over the generators of the unit
  unit_block->check_outputs();
  for( Index g = 0 ; g < unit_block->get_number_generators() ; ++g ) {
   v_generator_unit.emplace_back( u , g );
   std::vector< Index > carriers( unit_block->get_number_outputs( g ) );
   for( Index k = 0 ; k < carriers.size() ; ++k )
    carriers[ k ] = get_carrier_index(
				    unit_block->get_output_carrier( g , k ) );
   v_output_carrier_index.push_back( std::move( carriers ) );
   }
  }

 // the NetworkBlock of the carrier 0 that covers each instant
 v_network_of_instant.clear();
 for( Index n = 0 ; n < v_network_blocks.size() ; ++n )
  for( Index i = 0 ; i < v_network_blocks[ n ]->get_number_intervals() ; ++i )
   v_network_of_instant.emplace_back( n , i );
 }

/*--------------------------------------------------------------------------*/

UCBlock::Index UCBlock::get_carrier_index( const std::string & name ) const
{
 if( name == UnitBlock::kElectricity )
  return( 0 );
 for( Index c = 0 ; c < v_carrier_name.size() ; ++c )
  if( v_carrier_name[ c ] == name )
   return( c );
 return( Inf< Index >() );
 }

/*--------------------------------------------------------------------------*/

UCBlock::Index UCBlock::get_generator_node( Index g , Index c ) const
{
 if( c == 0 )
  return( node_of_generator( g ) );
 if( v_generator_carrier_node.num_elements() )
  return( v_generator_carrier_node[ g ][ c ] );
 // GeneratorNode in one dimension, or not given: the carrier has one node,
 // where every generator is
 return( 0 );
 }

/*--------------------------------------------------------------------------*/

double UCBlock::get_carrier_efficiency( Index t , Index g , Index c ) const
{
 const auto [ unit , generator ] = v_generator_unit[ g ];
 const auto unit_block = get_unit_block( unit );
 const auto & carriers = v_output_carrier_index[ g ];
 double eta = 0;
 for( Index k = 0 ; k < carriers.size() ; ++k )
  if( carriers[ k ] == c )
   eta += unit_block->get_output_efficiency( t , generator , k );
 return( eta );
 }

/*--------------------------------------------------------------------------*/

UCBlock::Index UCBlock::get_storage_pollutant_zone( Index p , Index u ) const
{
 return( get_unit_block( u )->get_number_generators() ?
         get_pollutant_zone( p , v_first_generator[ u ] ) :
         get_pollutant_zone_of_node( p , 0 ) );
 }

/*--------------------------------------------------------------------------*/

std::pair< UCBlock::Index , UCBlock::Index >
UCBlock::generator_location( Index g ) const
{
 for( Index c = 0 ; c < get_number_carriers() ; ++c ) {
  const auto node = get_generator_node( g , c );
  if( node != Inf< Index >() )
   return( std::make_pair( c , node ) );
  }
 return( std::make_pair( Inf< Index >() , Inf< Index >() ) );
 }

/*--------------------------------------------------------------------------*/

void UCBlock::set_node_injection_bounds( void )
{
 // for each NetworkBlock of each carrier and each interval it covers, the
 // minimum and the maximum of the power that the generators at each node
 // inject: a generator that is off injects 0, or the opposite of its fixed
 // consumption (on the carrier 0), and one that is on the efficiency
 // towards the carrier times a power between its minimum and its maximum
 auto bounds = [ this ]( Index c , Index t , NetworkBlock * nbi , Index i ) {
  const auto nodes = nbi->get_number_nodes();
  std::vector< double > min_node_injection( nodes , 0 );
  std::vector< double > max_node_injection( nodes , 0 );

  for( Index g = 0 ; g < f_number_elc_generators ; ++g ) {
   const auto node = get_generator_node( g , c );
   if( node == Inf< Index >() )
    continue;
   const auto [ unit , generator ] = v_generator_unit[ g ];
   const auto ub = get_unit_block( unit );
   // multiplier for fleet-level node injection bound: scale times kappa
   // times the upper bound of the design variable (1 if binary)
   const auto mult = ub->get_scale() * ub->get_design_ub() * ub->get_kappa();
   const auto eta = get_carrier_efficiency( t , g , c );
   const auto p_mn = eta * ub->get_min_power( t , generator );
   const auto p_mx = eta * ub->get_max_power( t , generator );
   auto fixed_consumption = c ? nullptr
                              : ub->get_fixed_consumption( generator );
   min_node_injection[ node ] += mult * std::min( { p_mn , p_mx ,
			   fixed_consumption ? - fixed_consumption[ t ] : 0 } );
   max_node_injection[ node ] += mult * std::max( { 0.0 , p_mn , p_mx } );
   }

  for( Index node = 0 ; node < nodes ; ++node ) {
   nbi->set_min_node_injection( min_node_injection[ node ] , node , i );
   nbi->set_max_node_injection( max_node_injection[ node ] , node , i );
   }
  };

 for( Index t = 0 ; t < v_network_of_instant.size() ; ++t ) {
  const auto [ n , i ] = v_network_of_instant[ t ];
  bounds( 0 , t , v_network_blocks[ n ] , i );
  }

 for( Index c = 1 ; c < get_number_carriers() ; ++c ) {
  const auto & cd = v_carrier[ c ];
  for( Index t = 0 ; t < cd.v_network_of_instant.size() ; ++t ) {
   const auto [ n , i ] = cd.v_network_of_instant[ t ];
   bounds( c , t , cd.v_network_blocks[ n ] , i );
   }
  }

 // if reactive power is considered, do the same for it, on the carrier 0
 if( f_has_reactive && ( ! v_network_blocks.empty() ) )
  for( Index t = 0 ; t < v_network_of_instant.size() ; ++t ) {
   const auto [ n , i ] = v_network_of_instant[ t ];
   auto nbi = v_network_blocks[ n ];
   // for each interval covered by the network, compute min and max
   // reactive power output by all generators
   std::vector< double > min_node_injection( nbi->get_number_nodes() , 0 );
   std::vector< double > max_node_injection( nbi->get_number_nodes() , 0 );

   for( Index g = 0 ; g < f_number_elc_generators ; ++g ) {
    const auto node = get_generator_node( g , 0 );
    if( node == Inf< Index >() )
     continue;
    const auto [ unit , generator ] = v_generator_unit[ g ];
    const auto ub = get_unit_block( unit );
    const auto mult = ub->get_scale() * ub->get_design_ub() *
                      ub->get_kappa();
    // the reactive power of a thermal unit lies in [ Q^mn , Q^mx ] when it
    // is off, and the bounds move by Q^{mn,on}, Q^{mx,on} when it is on: the
    // injection bounds cover both cases
    const auto tub = dynamic_cast< const ThermalUnitBlock * >( ub );
    double q_mn = ub->get_min_reactive_power( t , generator );
    double q_mx = ub->get_max_reactive_power( t , generator );
    if( tub ) {
     q_mn += std::min( 0.0 , tub->get_min_reactive_power_on( t , generator ) );
     q_mx += std::max( 0.0 , tub->get_max_reactive_power_on( t , generator ) );
     }
    min_node_injection[ node ] += mult * q_mn;
    max_node_injection[ node ] += mult * q_mx;
    }

   for( Index node = 0 ; node < nbi->get_number_nodes() ; ++node ) {
    nbi->set_min_reactive_node_injection( min_node_injection[ node ] ,
					  node , i );
    nbi->set_max_reactive_node_injection( max_node_injection[ node ] ,
					  node , i );
    }
   }
 }  // end( UCBlock::set_node_injection_bounds )

/*--------------------------------------------------------------------------*/

void UCBlock::balance_terms( Index c , Index t , Index node , bool reactive ,
			     std::vector< BalanceTerm > & terms ,
			     double & rhs )
{
 terms.clear();
 rhs = 0;

 const auto nodes = get_number_nodes( c );
 if( nodes > 1 ) {
  // the node injection of the NetworkBlock covering t comes first
  NetworkBlock * nbi;
  Index i = 0;
  if( c == 0 ) {
   const auto [ n , ii ] = v_network_of_instant[ t ];
   nbi = v_network_blocks[ n ];
   i = ii;
   }
  else {
   const auto [ n , ii ] = v_carrier[ c ].v_network_of_instant[ t ];
   nbi = v_carrier[ c ].v_network_blocks[ n ];
   i = ii;
   }

  auto injection = reactive ? nbi->get_reactive_node_injection( i )
                            : nbi->get_node_injection( i );
  terms.push_back( { &injection[ node ] , -1 , 0 , Inf< Index >() } );
  }
 else  // one node: the demand is the right-hand side
  if( c == 0 )
   rhs = reactive ? v_reactive_power_demand[ 0 ][ t ]
                  : v_active_power_demand[ 0 ][ t ];
  else
   if( ! v_carrier[ c ].v_demand.empty() )
    rhs = v_carrier[ c ].v_demand[ 0 ][ t ];

 for( Index g = 0 ; g < f_number_elc_generators ; ++g ) {
  const auto gnode = get_generator_node( g , c );
  if( ( gnode == Inf< Index >() ) || ( ( nodes > 1 ) && ( gnode != node ) ) )
   continue;

  const auto [ unit , generator ] = v_generator_unit[ g ];
  const auto unit_block = get_unit_block( unit );
  const auto scale = unit_block->get_scale();

  if( reactive ) {
   // the reactive power alone, the fixed consumption being an active power
   if( auto rp = unit_block->get_reactive_power( generator ) )
    terms.push_back( { &rp[ t ] , scale , 0 , unit } );
   continue;
   }

  // the active power, times the efficiency towards the carrier
  if( const auto eta = get_carrier_efficiency( t , g , c ) )
   if( auto ap = unit_block->get_active_power( generator ) )
    terms.push_back( { &ap[ t ] , scale * eta , 0 , unit } );

  // the fixed consumption fc[ t ] * ( 1 - u[ t ] ) of a generator that is off,
  // an electric power: fc[ t ] * u[ t ] on the left and fc[ t ] on the right
  // (a generator without commitment variables is always on)
  if( c == 0 )
   if( auto fc = unit_block->get_fixed_consumption( generator ) )
    if( fc[ t ] )
     if( auto u = unit_block->get_commitment( generator ) ) {
      const auto fixed_consumption = scale * fc[ t ];
      terms.push_back( { &u[ t ] , fixed_consumption , fixed_consumption ,
			 unit } );
      rhs += fixed_consumption;
      }
  }
 }  // end( UCBlock::balance_terms )

/*--------------------------------------------------------------------------*/

void UCBlock::fill_balance_constraints( MAFRC & rows , Index c ,
				        bool reactive )
{
 rows.resize( MAFRC_ext()[ f_time_horizon ][ get_number_nodes( c ) ] );

 std::vector< BalanceTerm > terms;
 double rhs;
 for( Index t = 0 ; t < f_time_horizon ; ++t )
  for( Index n = 0 ; n < rows.shape()[ 1 ] ; ++n ) {
   balance_terms( c , t , n , reactive , terms , rhs );
   LinearFunction::v_coeff_pair vc;
   vc.reserve( terms.size() );
   for( const auto & term : terms )
    vc.emplace_back( term.var , term.coeff );
   rows[ t ][ n ].set_both( rhs , eNoMod );
   rows[ t ][ n ].set_function( new LinearFunction( std::move( vc ) ) ,
				eNoMod );
   }
 }

/*--------------------------------------------------------------------------*/

std::vector< std::string > UCBlock::expected_dims( void ) const {
 static const std::vector< std::string > ed =
 { "TimeHorizon" , "NumberUnits" ,"NumberNetworks" , "NumberPrimaryZones" ,
   "NumberSecondaryZones" , "NumberInertiaZones" , "NumberPollutants" ,
   "NumberNodes" , "NumberLines" , "NumberBranches" ,
   "NumberElectricalGenerators" , "NumberGenerators" ,
   "TotalNumberPollutantZones" , "NumberStorages" , "NumberIntervals" ,
   "NumberCarriers" };

 auto ret = Block::expected_dims();
 ret.insert( ret.end() , ed.begin() , ed.end() );
 if( f_NetworkData ) {
  auto ndev = f_NetworkData->expected_dims();
  ret.insert( ret.end() , ndev.begin() , ndev.end() );
  }

 return( ret );
 }

/*--------------------------------------------------------------------------*/

std::vector< std::string > UCBlock::expected_vars( void ) const {
 static const std::vector< std::string > ev =
 { "ActivePowerDemand" , "ReactivePowerDemand" , "GeneratorNode" ,
   "PrimaryZones" , "PrimaryDemand" ,
   "SecondaryZones" , "SecondaryDemand" , "InertiaZones" , "InertiaDemand" ,
   "NumberPollutantZones" , "PollutantZones" , "PollutantBudget" ,
   "PollutantRho" , "PollutantMinBudget" , "PollutantStorageRho" ,
   "NetworkConstantTerms" , "NetworkBlockClassname" ,
   "NetworkDataClassname" , "CarrierName"
   };

 auto ret = Block::expected_vars();
 ret.insert( ret.end() , ev.begin() , ev.end() );
 if( f_NetworkData ) {
  auto nded = f_NetworkData->expected_vars();
  ret.insert( ret.end() , nded.begin() , nded.end() );
  }

 return( ret );
 }

/*--------------------------------------------------------------------------*/

void UCBlock::generate_abstract_constraints( Configuration * stcc )
{
 if( constraints_generated() )  // constraints have already been generated
  return;                       // nothing to do

 // generate abstract constraints in all the sub-Block
 Block::generate_abstract_constraints( stcc );

 // generate the abstract constraints of UCBlock
 generate_node_injection_constraints();
 if( f_has_reactive )
  generate_reactive_node_injection_constraints();
 generate_carrier_balance_constraints();

 generate_primary_demand_constraints();
 generate_secondary_demand_constraints();
 generate_inertia_demand_constraints();
 generate_pollutant_budget_constraints();

 set_constraints_generated();

}  // end( UCBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void UCBlock::generate_node_injection_constraints( void )
{
 // with one node there are no NetworkBlock, and the rows are those (1b) of
 // the active power demand [see balance_terms()]
 fill_balance_constraints( v_node_injection_Const , 0 , false );
 add_static_constraint( v_node_injection_Const , "node_injection_c" );
 }  // end( UCBlock::generate_node_injection_constraints )

/*--------------------------------------------------------------------------*/

void UCBlock::generate_reactive_node_injection_constraints( void )
{
 fill_balance_constraints( v_reactive_node_injection_Const , 0 , true );
 add_static_constraint( v_reactive_node_injection_Const ,
			"reactive_node_injection_c" );
 }  // end( UCBlock::generate_reactive_node_injection_constraints )

/*--------------------------------------------------------------------------*/

void UCBlock::generate_carrier_balance_constraints( void )
{
 // sized once, so that the rows of each carrier stay where they are
 v_carrier_balance_Const.resize( get_number_carriers() );
 for( Index c = 1 ; c < get_number_carriers() ; ++c ) {
  fill_balance_constraints( v_carrier_balance_Const[ c ] , c , false );
  add_static_constraint( v_carrier_balance_Const[ c ] ,
			 "carrier_balance_c_" + std::to_string( c ) );
  }
 }  // end( UCBlock::generate_carrier_balance_constraints )

/*--------------------------------------------------------------------------*/

void UCBlock::generate_primary_demand_constraints( void )
{
 if( f_number_primary_zones == 0 )
  return;

 v_PrimaryDemand_Const.resize(
		   MAFRC_ext()[ f_time_horizon ][ f_number_primary_zones ] );

 const auto number_nodes = get_number_nodes();

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {
  for( Index zone_id = 0 ; zone_id < f_number_primary_zones ; ++zone_id ) {

   auto lf = new LinearFunction();

   for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {

    if( ! node_belongs_to_primary_zone( node_id , zone_id ) )
     continue;

    Index elc_generator = 0;
    for( Index unit_id = 0 ; unit_id < f_number_units ; ++unit_id ) {

     const auto unit_block = get_unit_block( unit_id );
     const auto scale = unit_block->get_scale();

     for( Index generator = 0 ;
          generator < unit_block->get_number_generators() ;
          ++generator , ++elc_generator ) {

      if( ! generator_belongs_to_node( elc_generator , node_id ) )
       continue;

      if( auto primary_s_r =
       unit_block->get_primary_spinning_reserve( generator ) ) {
       // Now we add the primary reserve variable to the LinearFunction.
       auto primary_spinning_reserve = &primary_s_r[ t ];
       lf->add_variable( primary_spinning_reserve , scale );
      }

     }  // end( for( generator ) )
    }  // end( for( unit_id ) )
   }  // end( for( node_id ) )

   const auto demand = v_primary_demand[ zone_id ][ t ];
   v_PrimaryDemand_Const[ t ][ zone_id ].set_lhs( demand );
   v_PrimaryDemand_Const[ t ][ zone_id ].set_rhs( Inf< double >() );
   v_PrimaryDemand_Const[ t ][ zone_id ].set_function( lf );

  }  // end( for( zone_id ) )
 }  // end( for( t ) )

 add_static_constraint( v_PrimaryDemand_Const , "primary_demand_c" );

}  // end( UCBlock::generate_primary_demand_constraints )

/*--------------------------------------------------------------------------*/

void UCBlock::generate_secondary_demand_constraints( void )
{
 if( f_number_secondary_zones == 0 )
  return;

 v_SecondaryDemand_Const.resize(
  boost::multi_array< FRowConstraint , 2 >::extent_gen()
  [ f_time_horizon ][ f_number_secondary_zones ] );


 const auto number_nodes = get_number_nodes();

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {
  for( Index zone_id = 0 ; zone_id < f_number_secondary_zones ; ++zone_id ) {

   auto lf = new LinearFunction();

   for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {

    if( ! node_belongs_to_secondary_zone( node_id , zone_id ) )
     continue;

    Index elc_generator = 0;
    for( Index unit_id = 0 ; unit_id < f_number_units ; unit_id++ ) {

     const auto unit_block = get_unit_block( unit_id );
     const auto scale = unit_block->get_scale();

     for( Index generator = 0 ;
          generator < unit_block->get_number_generators() ;
          ++generator , ++elc_generator ) {

      if( ! generator_belongs_to_node( elc_generator , node_id ) )
       continue;

      if( auto secondary_s_r =
       unit_block->get_secondary_spinning_reserve( generator ) ) {
       // Now we add the secondary reserve variable to the LinearFunction.
       auto secondary_spinning_reserve = &secondary_s_r[ t ];
       lf->add_variable( secondary_spinning_reserve , scale );
      }

     }  // end( for( generator ) )
    }  // end( for( unit_id ) )
   }  // end( for( node_id ) )

   const auto demand = v_secondary_demand[ zone_id ][ t ];
   v_SecondaryDemand_Const[ t ][ zone_id ].set_lhs( demand );
   v_SecondaryDemand_Const[ t ][ zone_id ].set_rhs( Inf< double >() );
   v_SecondaryDemand_Const[ t ][ zone_id ].set_function( lf );

  }  // end( for( zone_id ) )
 }  // end( for( t ) )

 add_static_constraint( v_SecondaryDemand_Const , "secondary_demand_c" );

}  // end( UCBlock::generate_secondary_demand_constraints )

/*--------------------------------------------------------------------------*/

void UCBlock::generate_inertia_demand_constraints( void )
{
 if( f_number_inertia_zones == 0 )
  return;

 v_InertiaDemand_Const.resize(
  boost::multi_array< FRowConstraint , 2 >::extent_gen()
  [ f_time_horizon ][ f_number_inertia_zones ] );


 const auto number_nodes = get_number_nodes();

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {
  for( Index zone_id = 0 ; zone_id < f_number_inertia_zones ; ++zone_id ) {

   auto lf = new LinearFunction();

   for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {

    if( ! node_belongs_to_inertia_zone( node_id , zone_id ) )
     continue;

    Index elc_generator = 0;
    for( Index unit_id = 0 ; unit_id < f_number_units ; ++unit_id ) {

     const auto unit_block = get_unit_block( unit_id );
     const auto scale = unit_block->get_scale();

     for( Index generator = 0 ;
          generator < unit_block->get_number_generators() ;
          ++generator , ++elc_generator ) {

      if( ! generator_belongs_to_node( elc_generator , node_id ) )
       continue;

      auto commitment = unit_block->get_commitment( generator );
      auto inertia_commitment =
                         unit_block->get_inertia_commitment( generator );

      if( commitment && inertia_commitment ) {


       auto commitment_t = &commitment[ t ];
       auto coefficient = scale * inertia_commitment[ t ];
       lf->add_variable( commitment_t , coefficient );
      }

      auto active_power = unit_block->get_active_power( generator );
      auto inertia_power = unit_block->get_inertia_power( generator );

      if( active_power && inertia_power ) {

       auto active_power_t = &active_power[ t ];
       auto coefficient = scale * inertia_power[ t ];
       lf->add_variable( active_power_t , coefficient );
      }

     }  // end( for( generator ) )
    }  // end( for( unit_id ) )
   }  // end( for( node_id ) )

   const auto demand = v_inertia_demand[ zone_id ][ t ];
   v_InertiaDemand_Const[ t ][ zone_id ].set_lhs( demand );
   v_InertiaDemand_Const[ t ][ zone_id ].set_rhs( Inf< double >() );
   v_InertiaDemand_Const[ t ][ zone_id ].set_function( lf );

  }  // end( for( zone_id ) )
 }  // end( for( t ) )

 add_static_constraint( v_InertiaDemand_Const , "inertia_demand_c" );

}  // end( UCBlock::generate_inertia_demand_constraints )

/*--------------------------------------------------------------------------*/

template< class F >
void UCBlock::for_each_pollutant_term( Index p , F && visit ,
                                       Index unit ) const
{
 assert( ( unit == Inf< Index >() ) || ( unit < f_number_units ) );

 const auto number_zones = v_number_pollutant_zones[ p ];

 // all the units, or the given one alone
 const bool all = ( unit == Inf< Index >() );
 const Index first_unit = all ? 0 : unit;
 const Index end_unit = all ? f_number_units : unit + 1;
 Index elc_generator = v_first_generator[ first_unit ];
 Index storage = v_first_storage[ first_unit ];
 for( Index unit_id = first_unit ; unit_id < end_unit ; ++unit_id ) {

  const auto unit_block = get_unit_block( unit_id );

  // the storages of a unit are where its first generator is
  const auto storage_zone = get_storage_pollutant_zone( p , unit_id );

  for( Index generator = 0 ;
       generator < unit_block->get_number_generators() ;
       ++generator , ++elc_generator ) {

   const auto zone_id = get_pollutant_zone( p , elc_generator );
   if( zone_id >= number_zones )
    continue;  // the generator belongs to no zone of pollutant p

   auto active_power = unit_block->get_active_power( generator );
   if( ! active_power )
    continue;

   for( Index t = 0 ; t < f_time_horizon ; ++t )
    if( const auto rho = get_pollutant_rho( t , p , elc_generator ) )
     visit( unit_id , unit_block , zone_id , &active_power[ t ] , rho );
   }

  const auto number_storages = unit_block->get_number_storages();
  if( ( ! v_pollutant_storage_rho.empty() ) &&
      ( storage_zone < number_zones ) )
   for( Index s = 0 ; s < number_storages ; ++s )
    if( auto level = unit_block->get_storage_level( s ) )
     for( Index t = 0 ; t < f_time_horizon ; ++t )
      if( const auto sigma = get_pollutant_storage_rho( t , p ,
							storage + s ) )
       visit( unit_id , unit_block , storage_zone , &level[ t ] , sigma );

  storage += number_storages;
  }  // end( for( unit_id ) )

 }  // end( UCBlock::for_each_pollutant_term )

/*--------------------------------------------------------------------------*/

void UCBlock::generate_pollutant_budget_constraints( void )
{
 if( f_number_pollutants == 0 )
  return;

 // one constraint per zone of each pollutant: the pollutants may have a
 // different number of zones, so a two-dimensional array would have
 // useless rows
 v_PollutantBudget_Const.resize( f_number_pollutants );

 Index first_zone = 0;  // index in the budgets of zone 0 of pollutant p
 for( Index p = 0 ; p < f_number_pollutants ; ++p ) {
  const auto number_zones = v_number_pollutant_zones[ p ];
  v_PollutantBudget_Const[ p ].resize( number_zones );

  // the terms of the constraint of each zone of pollutant p: those of a
  // unit are consecutive, which is what update_pollutant_budget_constraints()
  // relies upon
  std::vector< LinearFunction::v_coeff_pair > vcp( number_zones );

  for_each_pollutant_term( p , [ & ]( Index , UnitBlock * unit_block ,
                                      Index zone_id , ColVariable * var ,
                                      double factor ) {
   vcp[ zone_id ].push_back( std::make_pair( var ,
                                             unit_block->get_scale() *
                                             factor ) );
   } );

  for( Index zone_id = 0 ; zone_id < number_zones ; ++zone_id ) {
   const auto k = first_zone + zone_id;
   auto & constraint = v_PollutantBudget_Const[ p ][ zone_id ];
   constraint.set_lhs( v_pollutant_min_budget[ k ] );
   constraint.set_rhs( v_pollutant_budget[ k ] );
   constraint.set_function(
                     new LinearFunction( std::move( vcp[ zone_id ] ) ) );
   }

  first_zone += number_zones;
  }  // end( for( p ) )

 add_static_constraint( v_PollutantBudget_Const , "pollutant_budget_c" );

 }  // end( UCBlock::generate_pollutant_budget_constraints )

/*--------------------------------------------------------------------------*/

void UCBlock::generate_objective( Configuration * objc )
{
 if( objective_generated() )  // Objective has already been generated
  return;                     // nothing to do

 for( auto block : v_Block )
  block->generate_objective();

 // the Objective of the UCBlock is the sum of the constant terms of the
 // NetworkBlock deleted in the bus case [see deserialize()], zero otherwise
 objective.set_function( new LinearFunction( LinearFunction::v_coeff_pair() ,
                         std::accumulate( v_network_constant_terms.begin() ,
                                          v_network_constant_terms.end() ,
                                          0.0 ) ) );

 // Set Block objective
 this->set_objective( &objective , eNoMod );

 set_objective_generated();

 } // end( UCBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*---------------- METHODS FOR CHECKING THE UCBlock ------------------------*/
/*--------------------------------------------------------------------------*/

bool UCBlock::is_feasible( bool useabstract , Configuration * fsbc )
{
 // Retrieve the tolerance and the type of violation.
 double tol = 1e-6;
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
  if( auto tc = dynamic_cast<
          SimpleConfiguration< std::pair< double , int > > * >( c ) ) {
   tol = tc->f_value.first;
   rel_viol = tc->f_value.second;
   return( true );
  }
  return( false );
 };

 if( ( ! extract_parameters( fsbc ) ) && f_BlockConfig )
  // if the given Configuration is not valid, try the one from the BlockConfig
  extract_parameters( f_BlockConfig->f_is_feasible_Configuration );

 // the sub-Block are checked with the same tolerance and type of violation,
 // unless they have their own in their BlockConfig
 SimpleConfiguration< std::pair< double , int > > subc(
                                 std::pair< double , int >( tol , rel_viol ) );
 for( const auto & sbi : this->get_nested_Blocks() ) {
  const auto bc = sbi->get_BlockConfig();
  if( ! sbi->is_feasible( useabstract ,
                          ( bc && bc->f_is_feasible_Configuration )
                          ? nullptr : & subc ) )
   return( false );
  }

 return(
  // Constraints: notice that the ZOConstraints are not checked, since the
  // corresponding check is made on the ColVariable
  RowConstraint::is_feasible( v_node_injection_Const , tol , rel_viol )
  && RowConstraint::is_feasible( v_PrimaryDemand_Const , tol , rel_viol )
  && RowConstraint::is_feasible( v_SecondaryDemand_Const , tol , rel_viol )
  && RowConstraint::is_feasible( v_InertiaDemand_Const , tol , rel_viol )
  && RowConstraint::is_feasible( v_PollutantBudget_Const , tol , rel_viol )
);

}  // end( UClUnitBlock::is_feasible )

/*--------------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/

Solution * UCBlock::get_Solution( Configuration *solc , bool emptys )
{
 int wsol = 7;
 if( ( ! solc ) && f_BlockConfig )
  solc = f_BlockConfig->f_solution_Configuration;

 if( auto tsolc = dynamic_cast< SimpleConfiguration< int > * >( solc ) )
  wsol = tsolc->f_value;

 auto *sol = new UCBlockSolution();

 if( wsol & 1 )
  sol->v_unit_Solution.resize( get_number_units() );

 if( ( wsol & 2 ) && ( get_number_nodes() > 1 ) )
  sol->v_network_Solution.resize( get_number_networks() );

 sol->f_compressed_network = wsol & 4;

 using mad2 = boost::multi_array< double , 2 >;

 if( wsol & 8 )
  sol->v_demand_duals.resize(
	    mad2::extent_gen()[ get_time_horizon() ][ get_number_nodes() ] );

 if( wsol & 16 )
  sol->v_primary_duals.resize(
    mad2::extent_gen()[ get_time_horizon() ][ get_number_primary_zones() ] );

 if( wsol & 32 )
  sol->v_secondary_duals.resize(
  mad2::extent_gen()[ get_time_horizon() ][ get_number_secondary_zones() ] );

 if( wsol & 64 )
  sol->v_inertia_duals.resize(
    mad2::extent_gen()[ get_time_horizon() ][ get_number_inertia_zones() ] );

 if( wsol & 128 )
  sol->v_pollutant_duals.resize( f_total_number_pollutant_zones );


 // the NetworkBlock and the balances (1c) of the carriers c > 0
 if( get_number_carriers() > 1 ) {
  if( wsol & 2 ) {
   sol->v_carrier_network_Solution.resize( get_number_carriers() );
   for( Index c = 1 ; c < get_number_carriers() ; ++c )
    sol->v_carrier_network_Solution[ c ].resize(
				      v_carrier[ c ].v_network_blocks.size() );
   }

  if( wsol & 8 ) {
   sol->v_carrier_duals.resize( get_number_carriers() );
   for( Index c = 1 ; c < get_number_carriers() ; ++c )
    sol->v_carrier_duals[ c ].resize(
		mad2::extent_gen()[ get_time_horizon() ][ get_number_nodes( c ) ] );
   }
  }

 if( ! emptys )
  sol->read( this );

 return( sol );

 }  // end( UCBlock::get_Solution )

/*--------------------------------------------------------------------------*/
/*------------ METHODS FOR LOADING, PRINTING & SAVING THE UCBlock ----------*/
/*--------------------------------------------------------------------------*/

void UCBlock::serialize( netCDF::NcGroup & group ) const
{
 Block::serialize( group );

 netCDF::NcDim NumberNodes;

 if( f_NetworkData ) {
  f_NetworkData->serialize( group );
  NumberNodes = group.getDim( "NumberNodes" );
  }
 else
  NumberNodes = group.addDim( "NumberNodes" , 1 );

 auto TimeHorizon = group.addDim( "TimeHorizon" , f_time_horizon );
 auto NumberUnits = group.addDim( "NumberUnits" , f_number_units );
 auto NumberNetworks = group.addDim( "NumberNetworks" , f_number_networks );
 auto NumberElectricalGenerators = group.addDim( "NumberGenerators" ,
                                                 f_number_elc_generators );

 auto TotalNumberPollutantZones = group.addDim( "TotalNumberPollutantZones" ,
                                            f_total_number_pollutant_zones );

 auto NumberPrimaryZones = group.addDim( "NumberPrimaryZones" ,
                                         f_number_primary_zones );
 auto NumberSecondaryZones = group.addDim( "NumberSecondaryZones" ,
                                           f_number_secondary_zones );
 auto NumberInertiaZones = group.addDim( "NumberInertiaZones" ,
                                         f_number_inertia_zones );
 auto NumberPollutants = group.addDim( "NumberPollutants" ,
                                       f_number_pollutants );

 if( ! v_carrier_name.empty() ) {
  auto NumberCarriers = group.addDim( "NumberCarriers" ,
				      get_number_carriers() );
  ::serialize( group , "CarrierName" , NumberCarriers , v_carrier_name );
  }

 ::serialize( group , "ActivePowerDemand" , netCDF::NcDouble() ,
              { NumberNodes , TimeHorizon } , v_active_power_demand );

 // note: in the multi-node case both demands have been handed over to the
 // NetworkBlock (which serialize them) at deserialize() time, so the two
 // multi_array are empty and nothing is written here
 ::serialize( group , "ReactivePowerDemand" , netCDF::NcDouble() ,
              { NumberNodes , TimeHorizon } , v_reactive_power_demand );

 ::serialize( group , "PrimaryZones" , netCDF::NcUint() ,
              NumberNodes , v_primary_zones );

 ::serialize( group , "PrimaryDemand" , netCDF::NcDouble() ,
              { NumberPrimaryZones , TimeHorizon } , v_primary_demand );

 ::serialize( group , "SecondaryZones" , netCDF::NcUint() ,
              NumberNodes , v_secondary_zones );

 ::serialize( group , "SecondaryDemand" , netCDF::NcDouble() ,
              { NumberSecondaryZones , TimeHorizon } , v_secondary_demand );

 ::serialize( group , "InertiaZones" , netCDF::NcUint() ,
              NumberNodes , v_inertia_zones );

 ::serialize( group , "InertiaDemand" , netCDF::NcDouble() ,
              { NumberInertiaZones , TimeHorizon } , v_inertia_demand );

 if( f_number_pollutants ) {
  ::serialize( group , "NumberPollutantZones" , netCDF::NcUint() ,
               NumberPollutants , v_number_pollutant_zones );

  ::serialize( group , "PollutantZones" , netCDF::NcUint() ,
               { NumberPollutants , NumberNodes } , v_pollutant_zones );

  ::serialize( group , "PollutantBudget" , netCDF::NcDouble() ,
               TotalNumberPollutantZones , v_pollutant_budget );

  if( std::any_of( v_pollutant_min_budget.begin() ,
                   v_pollutant_min_budget.end() ,
                   []( double b ) { return( b > -Inf< double >() ); } ) )
   ::serialize( group , "PollutantMinBudget" , netCDF::NcDouble() ,
                TotalNumberPollutantZones , v_pollutant_min_budget );

  // the first dimension is a singleton if the rates do not depend on time
  ::serialize( group , "PollutantRho" , netCDF::NcDouble() ,
               { TimeHorizon , NumberPollutants ,
                 NumberElectricalGenerators } , v_pollutant_rho ,
               false , true );

  if( ! v_pollutant_storage_rho.empty() ) {
   auto NumberStorages = group.addDim( "NumberStorages" , f_number_storages );
   ::serialize( group , "PollutantStorageRho" , netCDF::NcDouble() ,
                { TimeHorizon , NumberPollutants , NumberStorages } ,
                v_pollutant_storage_rho , false , true );
   }
  }

 if( std::any_of( v_network_constant_terms.begin() ,
                  v_network_constant_terms.end() ,
                  []( double cst ) { return( cst != 0 ); } ) )
  ::serialize( group , "NetworkConstantTerms" , netCDF::NcDouble() ,
               NumberNetworks , v_network_constant_terms );

 if( network_block_classname != "DCNetworkBlock" )
  ::serialize( group , "NetworkBlockClassname" , netCDF::NcString() ,
               network_block_classname );

 // the classname can be left implicit only when it coincides with the
 // default that deserialize() uses when the variable is absent
 if( network_data_classname != ( network_block_classname == "ACNetworkBlock"
                                 ? "ACNetworkData" : "DCNetworkData" ) )
  ::serialize( group , "NetworkDataClassname" , netCDF::NcString() ,
               network_data_classname );

 // GeneratorNode, with -1 where a generator has no node: in two dimensions
 // with more than one carrier, as a NcUint as long as no -1 is there
 const auto to_int = []( Index n ) {
  return( n == Inf< Index >() ? -1 : int( n ) );
  };
 if( v_generator_carrier_node.num_elements() ) {
  boost::multi_array< int , 2 > gn( boost::extents[ f_number_elc_generators ]
				                  [ get_number_carriers() ] );
  std::transform( v_generator_carrier_node.data() ,
		  v_generator_carrier_node.data() +
		  v_generator_carrier_node.num_elements() , gn.data() ,
		  to_int );
  ::serialize( group , "GeneratorNode" , netCDF::NcInt() ,
	       { NumberElectricalGenerators ,
		 group.getDim( "NumberCarriers" ) } , gn );
  }
 else
  if( std::find( v_generator_node.begin() , v_generator_node.end() ,
		 Inf< Index >() ) != v_generator_node.end() ) {
   std::vector< int > gn( v_generator_node.size() );
   std::transform( v_generator_node.begin() , v_generator_node.end() ,
		   gn.begin() , to_int );
   ::serialize( group , "GeneratorNode" , netCDF::NcInt() ,
		NumberElectricalGenerators , gn );
   }
  else
   ::serialize( group , "GeneratorNode" , netCDF::NcUint() ,
		NumberElectricalGenerators , v_generator_node );

 // Serialize sub-blocks

 for( Index i = 0 ; i < f_number_units ; ++i ) {
  auto sub_block = get_unit_block( i );
  auto sub_group = group.addGroup( "UnitBlock_" + std::to_string( i ) );
  sub_block->serialize( sub_group );
  }

 for( Index n = 0 ; n < f_number_networks ; ++n )
  if( auto sub_block = get_network_block( n ) ) {
   auto sub_group = group.addGroup( "NetworkBlock_" + std::to_string( n ) );
   sub_block->serialize( sub_group );
   }

 // the carriers other than the carrier 0, whose NetworkBlock are built by
 // deserialize() and written only if they were given in the group of the
 // carrier: their demand is
 for( Index c = 1 ; c < get_number_carriers() ; ++c ) {
  const auto & cd = v_carrier[ c ];
  auto cg = group.addGroup( "Carrier_" + std::to_string( c ) );

  netCDF::NcDim CNumberNodes;
  if( cd.f_NetworkData ) {
   cd.f_NetworkData->serialize( cg );
   CNumberNodes = cg.getDim( "NumberNodes" , netCDF::NcGroup::Current );
   }
  else
   CNumberNodes = cg.addDim( "NumberNodes" , 1 );

  if( cd.network_block_classname != "DCNetworkBlock" )
   ::serialize( cg , "NetworkBlockClassname" , netCDF::NcString() ,
		cd.network_block_classname );
  if( cd.network_data_classname != ( cd.network_block_classname ==
				     "ACNetworkBlock" ? "ACNetworkData"
				                      : "DCNetworkData" ) )
   ::serialize( cg , "NetworkDataClassname" , netCDF::NcString() ,
		cd.network_data_classname );

  const auto nn = get_number_nodes( c );
  if( ! cd.v_network_blocks.empty() ) {
   boost::multi_array< double , 2 > demand( boost::extents[ nn ]
					                    [ f_time_horizon ] );
   for( Index t = 0 ; t < f_time_horizon ; ++t ) {
    const auto [ nw , i ] = cd.v_network_of_instant[ t ];
    const auto d = cd.v_network_blocks[ nw ]->get_active_demand( i );
    for( Index n = 0 ; n < nn ; ++n )
     demand[ n ][ t ] = d ? d[ n ] : 0;
    }
   ::serialize( cg , "EnergyDemand" , netCDF::NcDouble() ,
		{ CNumberNodes , TimeHorizon } , demand );
   }
  else
   ::serialize( cg , "EnergyDemand" , netCDF::NcDouble() ,
		{ CNumberNodes , TimeHorizon } , cd.v_demand );

  if( f_number_pollutants )
   ::serialize( cg , "PollutantZones" , netCDF::NcUint() ,
		{ NumberPollutants , CNumberNodes } , cd.v_pollutant_zones ,
		false , true );

  if( cd.f_given_networks ) {
   cg.addDim( "NumberNetworks" , cd.v_network_blocks.size() );
   for( Index n = 0 ; n < cd.v_network_blocks.size() ; ++n ) {
    auto sub_group = cg.addGroup( "NetworkBlock_" + std::to_string( n ) );
    cd.v_network_blocks[ n ]->serialize( sub_group );
    }
   }
  }

 }  // end( UCBlock::serialize )

/*--------------------------------------------------------------------------*/
/*--------------- METHODS FOR READING THE DATA OF THE UCBlock --------------*/
/*--------------------------------------------------------------------------*/

int UCBlock::get_objective_sense( void ) const { return( Objective::eMin ); }

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

void UCBlock::add_Modification( sp_Mod mod , ChnlName chnl )
{
 // the scale of a unit going to 1 and back in get_zero_copy_linearization()
 // ends where it started: what the unit says of it is dropped
 if( const auto rtu = f_round_trip_unit.load() )
  if( const auto umod = dynamic_cast< const UnitBlockMod * >( mod.get() ) )
   if( ( umod->type() == UnitBlockMod::eScale ) &&
       ( umod->get_Block() == rtu ) )
    return;

 std::vector< Index > modified_units;
 std::vector< Index > inertia_units;
 bool injection = false;

 /* A scaled UnitBlock has its Variable multiplied by the scale wherever this
  * Block uses them [see UnitBlock::scale()], hence the rows that use them
  * have to be rewritten. The scaling of several units can arrive inside one
  * GroupModification, which is why the walk goes into the groups: looking
  * only at the top level would leave the rows of those units untouched. */

 std::function< void( const Modification * ) > collect =
  [ & ]( const Modification * m ) {
   if( const auto grp = dynamic_cast< const GroupModification * >( m ) ) {
    for( const auto & submod : grp->sub_Modifications() )
     collect( submod.get() );
    return;
    }

   if( const auto tmod = dynamic_cast< const UnitBlockMod * >( m ) ) {
    /* A change of the inertia power of a HydroUnitBlock changes the
     * coefficients of its active power in the inertia rows; the
     * HydroUnitBlock may be a unit of this UCBlock or one of the
     * HydroUnitBlock of a HydroSystemUnitBlock unit. */
    const bool inertia =
     ( tmod->type() == HydroUnitBlockMod::eSetInerP ) &&
     dynamic_cast< const HydroUnitBlockMod * >( m );
    // the bounds of the node injections depend on the power range of the
    // units, whatever their class [see UnitBlockMod::changes_power_range()]
    if( tmod->changes_power_range() )
     injection = true;
    if( ( tmod->type() == UnitBlockMod::eScale ) || inertia ) {
     /* The index is looked up among the units of this UCBlock rather than
      * among the sub-Block of the father of the unit: while a
      * LagrangianDualSolver is attached, that father is the LagBFunction
      * holding the unit alone, where the unit is always the 0-th. */
     const auto units_end = v_Block.begin() + f_number_units;
     auto it = std::find( v_Block.begin() , units_end , tmod->get_Block() );
     if( inertia && ( it == units_end ) && tmod->get_Block() )
      it = std::find( v_Block.begin() , units_end ,
                      tmod->get_Block()->get_f_Block() );
     if( it != units_end )
      ( inertia ? inertia_units : modified_units ).push_back(
                                      std::distance( v_Block.begin() , it ) );
     }
    }
   };

 collect( mod.get() );

 // the bounds of the node injections given to the NetworkBlocks follow the
 // power range of the units [see set_node_injection_bounds()]; replicate()
 // and resize_unit() compute them once, when they are done with the last
 // unit, and only if one has changed [see size_deferred()]
 if( injection || ( ! modified_units.empty() ) ) {
  if( f_defer_scaled_rows )
   f_deferred_injection_bounds = true;
  else
   set_node_injection_bounds();
  }

 if( ! modified_units.empty() ) {
  if( f_defer_scaled_rows )
   // a set of units is being sized, and their rows are rewritten once when
   // the last one is done [see size_deferred()]
   v_scaled_units.insert( v_scaled_units.end() , modified_units.begin() ,
                          modified_units.end() );
  else
   update_scaled_rows( modified_units );
  }

 if( ! inertia_units.empty() ) {
  std::sort( inertia_units.begin() , inertia_units.end() );
  inertia_units.erase( std::unique( inertia_units.begin() ,
                                    inertia_units.end() ) ,
                       inertia_units.end() );

  auto ichnl = open_channel();
  update_inertia_demand_constraints( inertia_units ,
                                     make_par( eNoBlck , ichnl ) );
  close_channel( ichnl );
  }

 Block::add_Modification( mod , chnl );
}

/*--------------------------------------------------------------------------*/

void UCBlock::update_scaled_rows( std::vector< Index > & scaled_units )
{
 if( scaled_units.empty() )
  return;

 // Sort the IDs of the modified units, and name each of them once: a group
 // can carry more than one scaling of the same unit
 std::sort( scaled_units.begin() , scaled_units.end() );
 scaled_units.erase( std::unique( scaled_units.begin() ,
                                  scaled_units.end() ) ,
                     scaled_units.end() );

 /* Each of these rewrites one row per time instant, and per node where there
  * are nodes: the whole reaction to the scaling travels in one channel, so
  * that a Solver able to write a set of coefficients, or of sides, in one
  * operation does that once instead of once per instant [see
  * MILPSolver::process_group_modification()]. */
 auto chnl = open_channel();
 const auto upar = make_par( eNoBlck , chnl );

 update_node_injection_constraints( scaled_units , upar );
 update_node_injection_constraints( scaled_units , upar , true );
 for( Index c = 1 ; c < get_number_carriers() ; ++c )
  update_balance_constraints( c , scaled_units , upar );
 update_primary_demand_constraints( scaled_units , upar );
 update_secondary_demand_constraints( scaled_units , upar );
 update_inertia_demand_constraints( scaled_units , upar );
 update_pollutant_budget_constraints( scaled_units , upar );

 close_channel( chnl );
}

/*--------------------------------------------------------------------------*/

void UCBlock::update_node_injection_constraints(
 const std::vector< Index > & modified_units , ModParam issueMod ,
 bool reactive )
{
 update_balance_constraints( 0 , modified_units , issueMod , reactive );
 }  // end( UCBlock::update_node_injection_constraints )

/*--------------------------------------------------------------------------*/

void UCBlock::update_balance_constraints(
			      Index c , const std::vector< Index > & modified_units ,
			      ModParam issueMod , bool reactive )
{
 auto & rows = reactive ? v_reactive_node_injection_Const :
               ( c == 0 ? v_node_injection_Const
                        : v_carrier_balance_Const[ c ] );

 if( ( ! constraints_generated() ) || rows.empty() || modified_units.empty() )
  return;

 // the terms are those of the rows, in the same order [see
 // fill_balance_constraints()]: the coefficients of those of the modified
 // units are written again, and so is the right-hand side, which the fixed
 // consumption of a scaled unit changes
 std::vector< BalanceTerm > terms;
 double rhs;
 for( Index t = 0 ; t < rows.shape()[ 0 ] ; ++t )
  for( Index n = 0 ; n < rows.shape()[ 1 ] ; ++n ) {
   auto & constraint = rows[ t ][ n ];
   balance_terms( c , t , n , reactive , terms , rhs );

   LinearFunction::Vec_FunctionValue coefficients;
   Subset subset;
   for( Index i = 0 ; i < terms.size() ; ++i )
    if( ( terms[ i ].unit != Inf< Index >() ) &&
        std::binary_search( modified_units.begin() , modified_units.end() ,
                            terms[ i ].unit ) ) {
     assert( i < constraint.get_num_active_var() );
     assert( constraint.get_active_var( i ) == terms[ i ].var );
     coefficients.push_back( terms[ i ].coeff );
     subset.push_back( i );
     }

   // the (abstract) Modifications that this issues do not concern this
   // UCBlock
   constraint.set_both( rhs , issueMod );
   if( ! subset.empty() )
    LF( constraint.get_function() )->modify_coefficients(
                std::move( coefficients ) , std::move( subset ) , true ,
                issueMod );
   }
 }  // end( UCBlock::update_balance_constraints )

/*--------------------------------------------------------------------------*/

void UCBlock::update_primary_demand_constraints(
                               const std::vector< Index > & modified_units ,
                               ModParam issueMod )
{
 if( ( ! constraints_generated() ) || ( v_PrimaryDemand_Const.empty() ) ||
     modified_units.empty() )
  return;  // there is nothing to be updated

 const auto number_nodes = get_number_nodes();

 /* The walk is the one of generate_primary_demand_constraints(), so that
  * the active Variable of each row are met in the order they have been
  * added: node by node, and within a node unit by unit. The Variable of a
  * unit with generators at several nodes of a zone are therefore not
  * consecutive, and the position of each one is counted along the walk. */

 for( Index t = 0 ; t < f_time_horizon ; ++t )
  for( Index zone_id = 0 ; zone_id < f_number_primary_zones ; ++zone_id ) {

   auto & constraint = v_PrimaryDemand_Const[ t ][ zone_id ];
   LinearFunction::Vec_FunctionValue coefficients;
   Subset subset;
   Index pos = 0;  // position of the next active Variable of the row

   for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {

    if( ! node_belongs_to_primary_zone( node_id , zone_id ) )
     continue;

    Index elc_generator = 0;
    for( Index unit_id = 0 ; unit_id < f_number_units ; ++unit_id ) {

     const auto unit_block = get_unit_block( unit_id );
     const bool modified = std::binary_search( modified_units.begin() ,
                                               modified_units.end() ,
                                               unit_id );

     for( Index generator = 0 ;
          generator < unit_block->get_number_generators() ;
          ++generator , ++elc_generator ) {

      if( ! generator_belongs_to_node( elc_generator , node_id ) )
       continue;

      if( auto primary_s_r =
           unit_block->get_primary_spinning_reserve( generator ) ) {
       assert( ( pos < constraint.get_num_active_var() ) &&
               ( constraint.get_active_var( pos ) == &primary_s_r[ t ] ) );
       if( modified ) {
        coefficients.push_back( unit_block->get_scale() );
        subset.push_back( pos );
        }
       ++pos;
       }
      }  // end( for( generator ) )
     }  // end( for( unit_id ) )
    }  // end( for( node_id ) )

   if( ! subset.empty() )
    LF( constraint.get_function() )->modify_coefficients(
                std::move( coefficients ) , std::move( subset ) , true ,
                issueMod );

   }  // end( for( zone_id ) )
 }  // end( UCBlock::update_primary_demand_constraints )

/*--------------------------------------------------------------------------*/

void UCBlock::update_secondary_demand_constraints(
                               const std::vector< Index > & modified_units ,
                               ModParam issueMod )
{
 if( ( ! constraints_generated() ) || ( v_SecondaryDemand_Const.empty() ) ||
     modified_units.empty() )
  return;  // there is nothing to be updated

 const auto number_nodes = get_number_nodes();

 /* The walk is the one of generate_secondary_demand_constraints(), so that
  * the active Variable of each row are met in the order they have been
  * added: node by node, and within a node unit by unit. The Variable of a
  * unit with generators at several nodes of a zone are therefore not
  * consecutive, and the position of each one is counted along the walk. */

 for( Index t = 0 ; t < f_time_horizon ; ++t )
  for( Index zone_id = 0 ; zone_id < f_number_secondary_zones ; ++zone_id ) {

   auto & constraint = v_SecondaryDemand_Const[ t ][ zone_id ];
   LinearFunction::Vec_FunctionValue coefficients;
   Subset subset;
   Index pos = 0;  // position of the next active Variable of the row

   for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {

    if( ! node_belongs_to_secondary_zone( node_id , zone_id ) )
     continue;

    Index elc_generator = 0;
    for( Index unit_id = 0 ; unit_id < f_number_units ; ++unit_id ) {

     const auto unit_block = get_unit_block( unit_id );
     const bool modified = std::binary_search( modified_units.begin() ,
                                               modified_units.end() ,
                                               unit_id );

     for( Index generator = 0 ;
          generator < unit_block->get_number_generators() ;
          ++generator , ++elc_generator ) {

      if( ! generator_belongs_to_node( elc_generator , node_id ) )
       continue;

      if( auto secondary_s_r =
           unit_block->get_secondary_spinning_reserve( generator ) ) {
       assert( ( pos < constraint.get_num_active_var() ) &&
               ( constraint.get_active_var( pos ) == &secondary_s_r[ t ] ) );
       if( modified ) {
        coefficients.push_back( unit_block->get_scale() );
        subset.push_back( pos );
        }
       ++pos;
       }
      }  // end( for( generator ) )
     }  // end( for( unit_id ) )
    }  // end( for( node_id ) )

   if( ! subset.empty() )
    LF( constraint.get_function() )->modify_coefficients(
                std::move( coefficients ) , std::move( subset ) , true ,
                issueMod );

   }  // end( for( zone_id ) )
 }  // end( UCBlock::update_secondary_demand_constraints )

/*--------------------------------------------------------------------------*/

void UCBlock::update_inertia_demand_constraints(
			        const std::vector< Index > & modified_units ,
			        ModParam issueMod )
{
 if( ( ! constraints_generated() ) || ( v_InertiaDemand_Const.empty() ) ||
     modified_units.empty() )
  return;  // there is nothing to be updated

 const auto number_nodes = get_number_nodes();

 /* The walk is the one of generate_inertia_demand_constraints(), so that
  * the active Variable of each row are met in the order they have been
  * added: node by node, and within a node unit by unit. The Variable of a
  * unit with generators at several nodes of a zone are therefore not
  * consecutive, and the position of each one is counted along the walk.
  * A term that the row does not have at its position (a term added to a
  * unit after the generation) is refused. */

 for( Index t = 0 ; t < f_time_horizon ; ++t )
  for( Index zone_id = 0 ; zone_id < f_number_inertia_zones ; ++zone_id ) {

   auto & constraint = v_InertiaDemand_Const[ t ][ zone_id ];
   LinearFunction::Vec_FunctionValue coefficients;
   Subset subset;
   Index pos = 0;  // position of the next active Variable of the row

   auto term = [ & ]( const ColVariable * var , Index unit_id ,
                      bool modified , double coefficient ) {
    if( ( pos >= constraint.get_num_active_var() ) ||
        ( constraint.get_active_var( pos ) != var ) )
     throw( std::logic_error( "UCBlock::update_inertia_demand_constraints:"
                              " unit " + std::to_string( unit_id ) +
                              " has an inertia term that the inertia rows, "
                              "generated without it, do not have" ) );
    if( modified ) {
     coefficients.push_back( coefficient );
     subset.push_back( pos );
     }
    ++pos;
    };

   for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {

    if( ! node_belongs_to_inertia_zone( node_id , zone_id ) )
     continue;

    Index elc_generator = 0;
    for( Index unit_id = 0 ; unit_id < f_number_units ; ++unit_id ) {

     const auto unit_block = get_unit_block( unit_id );
     const bool modified = std::binary_search( modified_units.begin() ,
                                               modified_units.end() ,
                                               unit_id );
     const auto scale = unit_block->get_scale();

     for( Index generator = 0 ;
          generator < unit_block->get_number_generators() ;
          ++generator , ++elc_generator ) {

      if( ! generator_belongs_to_node( elc_generator , node_id ) )
       continue;

      const auto commitment = unit_block->get_commitment( generator );
      const auto inertia_commitment =
                              unit_block->get_inertia_commitment( generator );
      if( commitment && inertia_commitment )
       term( &commitment[ t ] , unit_id , modified ,
             scale * inertia_commitment[ t ] );

      const auto active_power = unit_block->get_active_power( generator );
      const auto inertia_power = unit_block->get_inertia_power( generator );
      if( active_power && inertia_power )
       term( &active_power[ t ] , unit_id , modified ,
             scale * inertia_power[ t ] );

      }  // end( for( generator ) )
     }  // end( for( unit_id ) )
    }  // end( for( node_id ) )

   if( ! subset.empty() )
    LF( constraint.get_function() )->modify_coefficients(
                std::move( coefficients ) , std::move( subset ) , true ,
                issueMod );

   }  // end( for( zone_id ) )
 }  // end( UCBlock::update_inertia_demand_constraints )

/*--------------------------------------------------------------------------*/

void UCBlock::update_pollutant_budget_constraints(
			        const std::vector< Index > & modified_units ,
			        ModParam issueMod )
{
 if( ( ! constraints_generated() ) || v_PollutantBudget_Const.empty() ||
     modified_units.empty() )
  return;  // there is nothing to be updated

 for( Index p = 0 ; p < f_number_pollutants ; ++p ) {
  const auto number_zones = v_number_pollutant_zones[ p ];
  auto & constraints = v_PollutantBudget_Const[ p ];

  // for the constraint of each zone of pollutant p, the index of its next
  // active Variable and the coefficients that change, with their indices
  std::vector< Index > next_var_index( number_zones , 0 );
  std::vector< LinearFunction::Vec_FunctionValue > coefficients(
							     number_zones );
  std::vector< Subset > subset( number_zones );

  // the walk is the one of generate_pollutant_budget_constraints(), so that
  // the active Variables are met in the order they have been added
  for_each_pollutant_term( p , [ & ]( Index unit_id , UnitBlock * unit_block ,
                                      Index zone_id , ColVariable * ,
                                      double factor ) {
   if( std::binary_search( modified_units.begin() , modified_units.end() ,
                           unit_id ) ) {
    assert( next_var_index[ zone_id ] <
            constraints[ zone_id ].get_num_active_var() );
    assert( constraints[ zone_id ].get_active_var(
                      next_var_index[ zone_id ] )->get_Block() == unit_block );

    coefficients[ zone_id ].push_back( unit_block->get_scale() * factor );
    subset[ zone_id ].push_back( next_var_index[ zone_id ] );
    }

   ++next_var_index[ zone_id ];
   } );

  for( Index zone_id = 0 ; zone_id < number_zones ; ++zone_id )
   if( ! subset[ zone_id ].empty() )
    LF( constraints[ zone_id ].get_function() )->
     modify_coefficients( std::move( coefficients[ zone_id ] ) ,
                          std::move( subset[ zone_id ] ) , true , issueMod );
  }  // end( for( p ) )

 }  // end( UCBlock::update_pollutant_budget_constraints )

/*--------------------------------------------------------------------------*/

void UCBlock::update_node_injection_constraints( Index time ,
						 Index node_index ,
                                                 double demand )
{
 // the demand is already in v_active_power_demand, from which the
 // right-hand side is computed, net of the fixed consumption of the units
 // that are off [see balance_terms()]
 std::vector< BalanceTerm > terms;
 double rhs;
 balance_terms( 0 , time , node_index , false , terms , rhs );
 v_node_injection_Const[ time ][ node_index ].set_both( rhs , eNoBlck );
}

/*--------------------------------------------------------------------------*/

void UCBlock::set_active_power_demand( MF_dbl_it values ,
                                       Block::Subset && subset ,
                                       bool ordered ,
                                       c_ModParam issuePMod ,
                                       c_ModParam issueAMod )
{
 if( subset.empty() )
  return;

 if( ! v_network_blocks.empty() ) {
  set_network_active_demand( values , subset , issuePMod , issueAMod );
  return;
  }

 // Update the demand present in this UCBlock

 assert( ! v_active_power_demand.empty() );

 bool changed = false;

 for( auto index : subset ) {
  const auto node_index = index / f_time_horizon;
  const auto time = index % f_time_horizon;
  const auto demand = *( values++ );

  if( v_active_power_demand[ node_index ][ time ] != demand ) {
   changed = true;

   if( not_dry_run( issuePMod ) ) {
    // Change the physical representation
    v_active_power_demand[ node_index ][ time ] = demand;

    if( not_dry_run( issueAMod ) && constraints_generated() )
     // Change the abstract representation
     update_node_injection_constraints( time , node_index , demand );
   }
  }
 }

 if( ! changed ) // if nothing changes, return
  return;

 if( issue_pmod( issuePMod , UCBlockMod::changes_of(
                 UCBlockMod::eSetActD ) ) ) {
  if( ! ordered )
   std::sort( subset.begin() , subset.end() );

  Block::add_Modification( std::make_shared< UCBlockSbstMod >( this ,
                            UCBlockMod::eSetActD , std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );
 }
} // end( UCBlock::set_active_power_demand( subset ) )

/*--------------------------------------------------------------------------*/

template< class F >
void UCBlock::replicate_linking_terms( Index unit , F && term ,
                                       double & constant ) const
{
 const auto block = get_unit_block( unit );
 const auto time_horizon = get_time_horizon();
 const auto num_generators = block->get_number_generators();

 /* The index, among the electrical generators of this UCBlock, of the first
  * generator of this unit: the generators are numbered consecutively, unit
  * after unit. */

 const Index first_generator = v_first_generator[ unit ];

 // The balances of the carriers - - - - - - - - - - - - - - - - - - - - - - -

 // active and reactive [see for_each_balance_term()]: each term of the unit,
 // divided by its scale, is factor * var - constant in a row whose dual
 // value is dual, the factor multiplying the power as it multiplies the
 // fixed consumption
 for( const bool reactive : { false , true } )
  for_each_balance_term( unit , reactive ,
   [ & ]( const FRowConstraint & row , const ColVariable * var ,
          double factor , double cst ) {
    const auto dual = row.get_dual();
    term( var , dual * factor );
    constant -= dual * cst;
    } );

 // The reserve and inertia constraints - - - - - - - - - - - - - - - - - - -

 /* The rows of the primary and secondary reserves and those of the inertia
  * have the same shape: a row per zone and time instant, and a generator
  * contributes to it only if its node belongs to that zone. What the factor
  * multiplies is the reserve in those of the reserves, and the commitment
  * and the active power, weighted by their inertia, in those of the
  * inertia. The generators of the unit that are in a zone are found once,
  * and sorted by zone, then node, then generator, so that the terms of an
  * instant come zone by zone and node by node. */

 const auto zone_contribution =
  [ & ]( const auto & constraints , Index number_zones , auto zone_of ,
         auto contribute ) {

   if( constraints.empty() )
    return;

   std::vector< std::array< Index , 3 > > in_zone;
   for( Index g = 0 ; g < num_generators ; ++g ) {
    const auto zone = zone_of( first_generator + g );
    if( zone < number_zones )
     in_zone.push_back( { zone , node_of_generator( first_generator + g ) ,
                          g } );
    }
   std::sort( in_zone.begin() , in_zone.end() );

   for( Index t = 0 ; t < time_horizon ; ++t )
    for( const auto & [ zone , node , g ] : in_zone )
     contribute( g , t , constraints[ t ][ zone ].get_dual() );
   };

 zone_contribution( v_PrimaryDemand_Const , f_number_primary_zones ,
                    [ this ]( Index g ) { return( get_primary_zone( g ) ); } ,
                    [ & ]( Index g , Index t , double dual ) {
                     if( auto r = block->get_primary_spinning_reserve( g ) )
                      term( & r[ t ] , dual ); } );

 zone_contribution( v_SecondaryDemand_Const , f_number_secondary_zones ,
                    [ this ]( Index g ) {
                     return( get_secondary_zone( g ) ); } ,
                    [ & ]( Index g , Index t , double dual ) {
                     if( auto r = block->get_secondary_spinning_reserve( g ) )
                      term( & r[ t ] , dual ); } );

 zone_contribution( v_InertiaDemand_Const , f_number_inertia_zones ,
                    [ this ]( Index g ) { return( get_inertia_zone( g ) ); } ,
                    [ & ]( Index g , Index t , double dual ) {
                     const auto u = block->get_commitment( g );
                     const auto iu = block->get_inertia_commitment( g );
                     if( u && iu )
                      term( & u[ t ] , dual * iu[ t ] );
                     const auto p = block->get_active_power( g );
                     const auto ip = block->get_inertia_power( g );
                     if( p && ip )
                      term( & p[ t ] , dual * ip[ t ] ); } );

 // The bounds of the node injections - - - - - - - - - - - - - - - - - - - -

 constant += injection_bounds_linearization( unit , true );

 // The pollutant budget constraints - - - - - - - - - - - - - - - - - - - - -

 // in the constraint of each zone the factor multiplies rho * p for each
 // generator of the unit and sigma * v for each of its storages [see
 // for_each_pollutant_term()]
 for( Index p = 0 ; p < v_PollutantBudget_Const.size() ; ++p )
  for_each_pollutant_term( p , [ & ]( Index , UnitBlock * , Index zone_id ,
                                      ColVariable * var , double factor ) {
   term( var , v_PollutantBudget_Const[ p ][ zone_id ].get_dual() * factor );
   } , unit );

}  // end( UCBlock::replicate_linking_terms )

/*--------------------------------------------------------------------------*/

double UCBlock::injection_bounds_linearization( Index unit ,
                                                bool copies ) const
{
 if( ! constraints_generated() )
  return( 0 );

 /* A bound of a node injection is the sum, over the generators at the node,
  * of the scale times the upper bound of the design times the kappa of the
  * unit times what the generator can inject [see set_node_injection_bounds()]:
  * its derivative with respect to the scale, or to the kappa, is what
  * multiplies it there, on the side that the sign of the dual says, as for
  * the flow limits [see DCNetworkBlock::get_resize_linearization()]. */
 const auto block = get_unit_block( unit );
 const double factor = block->get_design_ub() *
                       ( copies ? block->get_kappa() : block->get_scale() );
 const auto obj_sign = ( get_objective_sense() == Objective::eMin ) ? - 1 : 1;

 const Index first_generator = v_first_generator[ unit ];

 double linearization = 0;
 const auto add = [ & ]( const BoxConstraint & row , double mn , double mx ) {
  const auto dual = row.get_dual();
  linearization -= dual * factor * ( ( obj_sign * dual > 0 ) ? mn : mx );
  };

 for( Index c = 0 ; c < get_number_carriers() ; ++c ) {
  const auto & networks = c ? v_carrier[ c ].v_network_blocks
                            : v_network_blocks;
  const auto & of_instant = c ? v_carrier[ c ].v_network_of_instant
                              : v_network_of_instant;
  for( Index t = 0 ; t < of_instant.size() ; ++t ) {
   const auto [ n , i ] = of_instant[ t ];
   const auto [ holder , j ] = bounds_holder( networks[ n ] , i );
   if( ! holder )
    continue;
   const auto & rows = holder->get_node_injection_bounds();
   if( ! rows.num_elements() )
    continue;
   for( Index g = 0 ; g < block->get_number_generators() ; ++g ) {
    const auto node = get_generator_node( first_generator + g , c );
    if( node == Inf< Index >() )
     continue;
    const auto eta = get_carrier_efficiency( t , first_generator + g , c );
    const auto p_mn = eta * block->get_min_power( t , g );
    const auto p_mx = eta * block->get_max_power( t , g );
    const auto fc = c ? nullptr : block->get_fixed_consumption( g );
    add( rows[ node ][ j ] ,
         std::min( { p_mn , p_mx , fc ? - fc[ t ] : 0 } ) ,
         std::max( { 0.0 , p_mn , p_mx } ) );
    }
   }
  }

 // the reactive ones, on the carrier 0, of an ACNetworkBlock, also as the
 // sub-network of a DesignNetworkBlock
 if( f_has_reactive ) {
  const auto thermal = dynamic_cast< const ThermalUnitBlock * >( block );
  for( Index t = 0 ; t < v_network_of_instant.size() ; ++t ) {
   const auto [ n , i ] = v_network_of_instant[ t ];
   const auto ac = dynamic_cast< const ACNetworkBlock * >(
                          bounds_holder( v_network_blocks[ n ] , i ).first );
   if( ( ! ac ) || ac->get_reactive_node_injection_bounds().empty() )
    continue;
   const auto & rows = ac->get_reactive_node_injection_bounds();
   for( Index g = 0 ; g < block->get_number_generators() ; ++g ) {
    const auto node = get_generator_node( first_generator + g , 0 );
    if( node == Inf< Index >() )
     continue;
    double q_mn = block->get_min_reactive_power( t , g );
    double q_mx = block->get_max_reactive_power( t , g );
    if( thermal ) {
     q_mn += std::min( 0.0 , thermal->get_min_reactive_power_on( t , g ) );
     q_mx += std::max( 0.0 , thermal->get_max_reactive_power_on( t , g ) );
     }
    add( rows[ node ] , q_mn , q_mx );
    }
   }
  }

 return( linearization );
 }  // end( UCBlock::injection_bounds_linearization )

bool UCBlock::has_design( void ) const
{
 for( Index u = 0 ; u < f_number_units ; ++u )
  if( get_unit_block( u )->has_design() )
   return( true );

 for( Index c = 0 ; c < get_number_carriers() ; ++c )
  for( auto nb : get_network_blocks( c ) )
   if( auto dnb = dynamic_cast< const DesignNetworkBlock * >( nb ) ;
       dnb && dnb->get_number_design_lines() )
    return( true );

 return( false );
 }

/*--------------------------------------------------------------------------*/

double UCBlock::get_replicate_linearization( Index unit ) const
{
 if( ! constraints_generated() )
  throw( std::logic_error( "UCBlock::get_replicate_linearization: the "
                           "abstract representation has not been generated, "
                           "so there are no rows to read the duals of" ) );

 const auto block = get_unit_block( unit );
 const auto scale = block->get_scale();

 // with no copies the power of the unit is 0 whatever the duals are
 if( scale == 0 )
  return( get_zero_copy_linearization( unit ) );

 // the rows above the unit, at the solution
 double linearization = 0;
 replicate_linking_terms( unit ,
                          [ & linearization ]( const ColVariable * var ,
                                               double coeff ) {
                           linearization += coeff * var->get_value(); } ,
                          linearization );

 /* The Objective of a UnitBlock carrying the factor has the form k f(x)
  * [see UnitBlock.h], so what it contributes here is f(x), the value divided
  * by k. */

 const auto objective =
  dynamic_cast< FRealObjective * >( block->get_objective() );
 if( ! objective )
  throw( std::logic_error( "UCBlock::get_replicate_linearization: unit " +
                           std::to_string( unit ) + " has no FRealObjective, "
                           "of which the derivative counts the value" ) );
 objective->compute();
 linearization += objective->value() / scale;

 return( linearization );

} // end( UCBlock::get_replicate_linearization )

/*--------------------------------------------------------------------------*/

double UCBlock::get_zero_copy_linearization( Index unit ) const
{
 /* The copy is the mirror of the abstract representation of the unit [see
  * AbstractBlock::mirror()], taken while the unit is scaled to 1 so that its
  * Objective is that of one copy, f^1; it then sits in an AbstractBlock whose
  * Objective holds the terms of the rows above the unit on the Variable of
  * the copy, and the Solver of the AbstractBlock sums that Objective and
  * the one of the copy. Whatever the mirror cannot reproduce is left out of
  * the copy, which is then a relaxation of one copy of the unit: the minimum
  * is not larger than the true one, and since the number of copies cannot
  * decrease below 0 the linearization stays below the value function. */

 const auto block = get_unit_block( unit );

 // the Solver of the copy: the one the extra part of the BlockConfig gives
 BlockSolverConfig * bsc = nullptr;
 if( f_BlockConfig )
  bsc = dynamic_cast< BlockSolverConfig * >(
                                        f_BlockConfig->f_extra_Configuration );

 if( ( ! bsc ) || bsc->get_SolverNames().empty() )
  throw( std::logic_error( "UCBlock::get_replicate_linearization: unit " +
                           std::to_string( unit ) + " has no copies, and "
                           "the derivative is then the minimum over one "
                           "copy of the unit, for which a Solver is needed: "
                           "the extra Configuration of the BlockConfig of "
                           "this UCBlock is not a BlockSolverConfig with a "
                           "Solver" ) );

 /* The unit is scaled to 1 and back to its scale, 0, which leaves it as it
  * was found. The unit tells this UCBlock of each of the two changes, even
  * if asked to issue no Modification [see UnitBlock::issue_unit_mod()]:
  * add_Modification() drops them while f_round_trip_unit is the unit, so
  * that the rows that carry the factor and the bounds of the node
  * injections are not rewritten, and nothing is passed on to the Solver or
  * the father of this UCBlock. The mutex makes two calls on the same
  * UCBlock take turns. */

 auto copy = new AbstractBlock();
 {
  std::lock_guard< std::mutex > guard( f_round_trip_mutex );
  f_round_trip_unit = block;
  const auto scale = block->get_scale();
  try {
   block->scale( 1.0 , eNoMod , eNoMod );
   copy->mirror( block );
   }
  catch( ... ) {
   block->scale( scale , eNoMod , eNoMod );
   f_round_trip_unit = nullptr;
   delete copy;
   throw;
   }
  block->scale( scale , eNoMod , eNoMod );
  f_round_trip_unit = nullptr;
  }

 auto wrapper = new AbstractBlock();
 wrapper->add_nested_Block( copy );

 // what is built for the copy goes when it has been solved, or when
 // anything throws: the Solver registered for this solution only, which a
 // cleared BlockSolverConfig removes [see BlockSolverConfig::clear()], the
 // AbstractBlock with the copy, and its Objective
 FRealObjective * objective = nullptr;
 BlockSolverConfig * cbsc = nullptr;
 const auto cleanup = [ & ]() {
  if( cbsc ) {
   cbsc->clear();
   cbsc->apply( wrapper );
   delete cbsc;
   }
  delete wrapper;   // and the copy it holds
  delete objective;
  };

 double constant = 0;
 double value = 0;
 try {
  // a Variable of the unit can be in several rows, say the active power in
  // the node injection, an inertia zone and a pollutant budget, while a
  // LinearFunction takes each Variable once: the coefficients are summed
  std::map< ColVariable * , double > coefficients;
  replicate_linking_terms( unit ,
                           [ & ]( const ColVariable * var , double coeff ) {
                            if( coeff == 0 )
                             return;
                            auto cvar = copy->mirror_of( var );
                            if( ! cvar )
                             throw( std::logic_error(
                              "UCBlock::get_replicate_linearization: a "
                              "Variable of the unit " + std::to_string( unit )
                              + " has no copy" ) );
                            coefficients[ cvar ] += coeff; } ,
                           constant );

  objective = new FRealObjective( wrapper , new LinearFunction(
               LinearFunction::v_coeff_pair( coefficients.begin() ,
                                             coefficients.end() ) ) );
  objective->set_sense( copy->get_objective_sense() , eNoMod );
  wrapper->set_objective( objective , eNoMod );

  cbsc = bsc->clone();
  cbsc->apply( wrapper );
  if( wrapper->get_registered_solvers().empty() )
   throw( std::logic_error( "UCBlock::get_replicate_linearization: the "
                            "BlockSolverConfig of the copy of unit " +
                            std::to_string( unit ) +
                            " registers no Solver" ) );

  // the value of a solution bounds the minimum only from above, while the
  // cut needs it bounded from below: the lower bound is one also when the
  // solve stopped before the optimum, and is not finite if nothing is
  // proven
  auto solver = wrapper->get_registered_solvers().front();
  solver->compute();
  value = solver->get_lb();
  if( ! ( std::abs( value ) < Inf< double >() ) )
   throw( std::logic_error( "UCBlock::get_replicate_linearization: the copy "
                            "of unit " + std::to_string( unit ) +
                            " has no finite lower bound" ) );
  }
 catch( ... ) {
  cleanup();
  throw;
  }

 cleanup();

 return( value + constant );

} // end( UCBlock::get_zero_copy_linearization )

/*--------------------------------------------------------------------------*/

void UCBlock::set_active_power_demand( MF_dbl_it values , Block::Range rng ,
                                       c_ModParam issuePMod ,
                                       c_ModParam issueAMod )
{
 rng.second = std::min( rng.second , get_number_nodes() * f_time_horizon );

 if( rng.first >= rng.second )
  return;

 if( ! v_network_blocks.empty() ) {
  Subset index( rng.second - rng.first );
  std::iota( index.begin() , index.end() , rng.first );
  set_network_active_demand( values , index , issuePMod , issueAMod );
  return;
  }

 // Update the demand present in this UCBlock

 assert( ! v_active_power_demand.empty() );

 bool changed = false;

 for( Index index = rng.first ; index < rng.second ; ++index ) {
  const auto node_index = index / f_time_horizon;
  const auto time = index % f_time_horizon;
  const auto demand = *( values++ );

  if( v_active_power_demand[ node_index ][ time ] != demand ) {
   changed = true;

   if( not_dry_run( issuePMod ) ) {
    // Change the physical representation
    v_active_power_demand[ node_index ][ time ] = demand;

    if( not_dry_run( issueAMod ) && constraints_generated() )
     // Change the abstract representation
     update_node_injection_constraints( time , node_index , demand );
   }
  }
 }

 if( ! changed )  // if nothing changes, return
  return;

 // issue a Physical Modification
 if( issue_pmod( issuePMod , UCBlockMod::changes_of(
                 UCBlockMod::eSetActD ) ) )
  Block::add_Modification( std::make_shared< UCBlockRngdMod >(
                            this , UCBlockMod::eSetActD , rng ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( UCBlock::set_active_power_demand( range ) )

/*--------------------------------------------------------------------------*/

void UCBlock::set_network_active_demand( MF_dbl_it values ,
					 const Subset & index ,
					 c_ModParam issuePMod ,
					 c_ModParam issueAMod )
{
 // the NetworkBlock covering each instant, and the interval of the instant
 // in it (f_number_networks if none does)
 std::vector< std::pair< Index , Index > > where(
				    f_time_horizon , { f_number_networks , 0 } );
 for( Index n = 0 , t = 0 ; n < f_number_networks ; ++n ) {
  auto nb = v_network_blocks[ n ];
  if( ! nb )
   continue;
  const auto ni = nb->get_number_intervals();
  for( Index i = 0 ; ( i < ni ) && ( t < f_time_horizon ) ; ++i )
   where[ t++ ] = { n , i };
  }

 // the entries and the values of each NetworkBlock, in the order of index
 std::vector< Subset > entries( f_number_networks );
 std::vector< std::vector< double > > demand( f_number_networks );
 for( auto idx : index ) {
  const auto node_index = idx / f_time_horizon;
  const auto time = idx % f_time_horizon;
  const auto [ n , i ] = where[ time ];
  if( n >= f_number_networks )
   throw( std::logic_error( "UCBlock::set_active_power_demand: time index "
			    "out of range of NetworkBlocks" ) );
  auto nb = v_network_blocks[ n ];
  const auto nb_nodes = nb->get_number_nodes();
  if( node_index >= nb_nodes )
   throw( std::out_of_range(
    "UCBlock::set_active_power_demand: node index out of range (node=" +
    std::to_string( node_index ) + ", nodes=" + std::to_string( nb_nodes ) +
    ")" ) );

  entries[ n ].push_back( i * nb_nodes + node_index );
  demand[ n ].push_back( *( values++ ) );
  }

 for( Index n = 0 ; n < f_number_networks ; ++n )
  if( ! entries[ n ].empty() )
   v_network_blocks[ n ]->set_active_demand( demand[ n ].cbegin() ,
					     std::move( entries[ n ] ) ,
					     false , issuePMod , issueAMod );

 }  // end( UCBlock::set_network_active_demand )

/*--------------------------------------------------------------------------*/

bool UCBlock::set_carrier_demand( Index c , Index n , Index t ,
				  MF_dbl_it value , c_ModParam issuePMod ,
				  c_ModParam issueAMod )
{
 auto & cd = v_carrier[ c ];

 if( ! cd.v_network_blocks.empty() ) {
  // more than one node: the demand is in the NetworkBlock covering the
  // instant, which issues its own Modification; the entries of its demand
  // are numbered interval after interval, as in set_active_power_demand()
  const auto [ nw , i ] = cd.v_network_of_instant[ t ];
  auto nbi = cd.v_network_blocks[ nw ];
  const auto old = nbi->get_active_demand( i ) ?
                   nbi->get_active_demand( i )[ n ] : 0;
  if( old == *value )
   return( false );
  const auto flat = i * nbi->get_number_nodes() + n;
  nbi->set_active_demand( value , Range( flat , flat + 1 ) , issuePMod ,
			  issueAMod );
  return( true );
  }

 // one node: the demand is the right-hand side of the balance (1c)
 if( cd.v_demand.empty() ) {
  if( *value == 0 )
   return( false );
  cd.v_demand.resize( boost::extents[ 1 ][ f_time_horizon ] );
  std::fill_n( cd.v_demand.data() , f_time_horizon , 0.0 );
  }

 if( cd.v_demand[ 0 ][ t ] == *value )
  return( false );

 if( not_dry_run( issuePMod ) ) {
  cd.v_demand[ 0 ][ t ] = *value;

  if( not_dry_run( issueAMod ) && constraints_generated() ) {
   std::vector< BalanceTerm > terms;
   double rhs;
   balance_terms( c , t , 0 , false , terms , rhs );
   v_carrier_balance_Const[ c ][ t ][ 0 ].set_both( rhs , issueAMod );
   }
  }

 return( true );
 }

/*--------------------------------------------------------------------------*/

void UCBlock::set_energy_demand( MF_dbl_it values , Block::Subset && subset ,
				 bool ordered , c_ModParam issuePMod ,
				 c_ModParam issueAMod )
{
 if( subset.empty() )
  return;

 // the indices of the carrier 0 are those of set_active_power_demand(),
 // which is called for them; those of the other carriers are handled here
 const Index first = get_number_nodes() * f_time_horizon;
 Subset active , others;
 std::vector< double > active_values;
 for( auto index : subset ) {
  const auto value = values++;
  if( index < first ) {
   active.push_back( index );
   active_values.push_back( *value );
   continue;
   }

  // the carrier and the node of index
  Index m = index / f_time_horizon - get_number_nodes();
  Index c = 1;
  while( ( c < get_number_carriers() ) && ( m >= get_number_nodes( c ) ) )
   m -= get_number_nodes( c++ );
  if( c >= get_number_carriers() )
   throw( std::out_of_range( "UCBlock::set_energy_demand: index " +
			     std::to_string( index ) + " out of range" ) );

  if( set_carrier_demand( c , m , index % f_time_horizon , value ,
			  issuePMod , issueAMod ) )
   others.push_back( index );
  }

 if( ! active.empty() )
  set_active_power_demand( active_values.begin() , std::move( active ) ,
			   ordered , issuePMod , issueAMod );

 if( others.empty() )
  return;

 if( issue_pmod( issuePMod , UCBlockMod::changes_of(
                 UCBlockMod::eSetEnD ) ) ) {
  if( ! ordered )
   std::sort( others.begin() , others.end() );

  Block::add_Modification( std::make_shared< UCBlockSbstMod >( this ,
                            UCBlockMod::eSetEnD , std::move( others ) ) ,
                           Observer::par2chnl( issuePMod ) );
  }
 }  // end( UCBlock::set_energy_demand( subset ) )

/*--------------------------------------------------------------------------*/

void UCBlock::set_energy_demand( MF_dbl_it values , Block::Range rng ,
				 c_ModParam issuePMod ,
				 c_ModParam issueAMod )
{
 // all the nodes of all the carriers
 Index nodes = 0;
 for( Index c = 0 ; c < get_number_carriers() ; ++c )
  nodes += get_number_nodes( c );
 rng.second = std::min( rng.second , nodes * f_time_horizon );

 if( rng.first >= rng.second )
  return;

 // the part of the carrier 0, if any, goes to set_active_power_demand()
 const Index first = get_number_nodes() * f_time_horizon;
 if( rng.first < first ) {
  const auto last = std::min( rng.second , first );
  set_active_power_demand( values , Range( rng.first , last ) ,
			   issuePMod , issueAMod );
  values += last - rng.first;
  rng.first = last;
  }

 bool changed = false;
 for( Index index = rng.first ; index < rng.second ; ++index ) {
  Index m = index / f_time_horizon - get_number_nodes();
  Index c = 1;
  while( m >= get_number_nodes( c ) )
   m -= get_number_nodes( c++ );
  if( set_carrier_demand( c , m , index % f_time_horizon , values++ ,
			  issuePMod , issueAMod ) )
   changed = true;
  }

 if( changed && issue_pmod( issuePMod , UCBlockMod::changes_of(
				              UCBlockMod::eSetEnD ) ) )
  Block::add_Modification( std::make_shared< UCBlockRngdMod >(
                            this , UCBlockMod::eSetEnD , rng ) ,
                           Observer::par2chnl( issuePMod ) );
 }  // end( UCBlock::set_energy_demand( range ) )

/*--------------------------------------------------------------------------*/

bool UCBlock::set_pollutant_budget_k( Index k , double budget , bool lower ,
                                      ModParam issuePMod , ModParam issueAMod )
{
 if( k >= f_total_number_pollutant_zones )
  throw( std::out_of_range( lower ?
                "UCBlock::set_pollutant_min_budget: index out of range" :
                "UCBlock::set_pollutant_budget: index out of range" ) );

 auto & current = lower ? v_pollutant_min_budget[ k ]
                        : v_pollutant_budget[ k ];
 if( current == budget )
  return( false );

 if( not_dry_run( issuePMod ) ) {
  current = budget;

  if( not_dry_run( issueAMod ) && constraints_generated() ) {
   // the pollutant p and the zone z of the index k
   Index p = 0;
   Index z = k;
   while( z >= v_number_pollutant_zones[ p ] )
    z -= v_number_pollutant_zones[ p++ ];

   if( lower )
    v_PollutantBudget_Const[ p ][ z ].set_lhs( budget , issueAMod );
   else
    v_PollutantBudget_Const[ p ][ z ].set_rhs( budget , issueAMod );
   }
  }

 return( true );

 }  // end( UCBlock::set_pollutant_budget_k )

/*--------------------------------------------------------------------------*/

void UCBlock::set_pollutant_bounds( MF_dbl_it values ,
                                    Block::Subset && subset ,
                                    bool ordered , bool lower ,
                                    c_ModParam issuePMod ,
                                    c_ModParam issueAMod )
{
 if( subset.empty() )
  return;

 // the abstract representation is changed right here, and all the changes
 // of the sides travel together in one channel
 auto amod = un_ModBlock( issueAMod );
 const bool grouped = not_dry_run( issuePMod ) && not_dry_run( amod ) &&
                      constraints_generated();
 if( grouped )
  amod = make_par( par2mod( amod ) , open_channel( par2chnl( amod ) ) );

 bool changed = false;
 for( auto k : subset )
  if( set_pollutant_budget_k( k , *( values++ ) , lower , issuePMod , amod ) )
   changed = true;

 if( grouped )
  close_channel( par2chnl( amod ) );

 if( changed && issue_pmod( issuePMod , UCBlockMod::changes_of(
                            UCBlockMod::eSetPolMinB ) ) ) {
  if( ! ordered )
   std::sort( subset.begin() , subset.end() );

  Block::add_Modification( std::make_shared< UCBlockSbstMod >( this ,
                            lower ? UCBlockMod::eSetPolMinB :
                                    UCBlockMod::eSetPolB ,
                            std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );
  }
 }  // end( UCBlock::set_pollutant_bounds( subset ) )

/*--------------------------------------------------------------------------*/

void UCBlock::set_pollutant_bounds( MF_dbl_it values , Block::Range rng ,
                                    bool lower , c_ModParam issuePMod ,
                                    c_ModParam issueAMod )
{
 rng.second = std::min( rng.second , f_total_number_pollutant_zones );
 if( rng.first >= rng.second )
  return;

 // the abstract representation is changed right here, and all the changes
 // of the sides travel together in one channel
 auto amod = un_ModBlock( issueAMod );
 const bool grouped = not_dry_run( issuePMod ) && not_dry_run( amod ) &&
                      constraints_generated();
 if( grouped )
  amod = make_par( par2mod( amod ) , open_channel( par2chnl( amod ) ) );

 bool changed = false;
 for( Index k = rng.first ; k < rng.second ; ++k )
  if( set_pollutant_budget_k( k , *( values++ ) , lower , issuePMod , amod ) )
   changed = true;

 if( grouped )
  close_channel( par2chnl( amod ) );

 if( changed && issue_pmod( issuePMod , UCBlockMod::changes_of(
                            UCBlockMod::eSetPolMinB ) ) )
  Block::add_Modification( std::make_shared< UCBlockRngdMod >( this ,
                            lower ? UCBlockMod::eSetPolMinB :
                                    UCBlockMod::eSetPolB , rng ) ,
                           Observer::par2chnl( issuePMod ) );

 }  // end( UCBlock::set_pollutant_bounds( range ) )

/*--------------------------------------------------------------------------*/

void UCBlock::set_pollutant_budget( MF_dbl_it values ,
                                    Block::Subset && subset , bool ordered ,
                                    c_ModParam issuePMod ,
                                    c_ModParam issueAMod )
{
 set_pollutant_bounds( values , std::move( subset ) , ordered , false ,
                       issuePMod , issueAMod );
 }

/*--------------------------------------------------------------------------*/

void UCBlock::set_pollutant_budget( MF_dbl_it values , Block::Range rng ,
                                    c_ModParam issuePMod ,
                                    c_ModParam issueAMod )
{
 set_pollutant_bounds( values , rng , false , issuePMod , issueAMod );
 }

/*--------------------------------------------------------------------------*/

void UCBlock::set_pollutant_min_budget( MF_dbl_it values ,
                                        Block::Subset && subset ,
                                        bool ordered ,
                                        c_ModParam issuePMod ,
                                        c_ModParam issueAMod )
{
 set_pollutant_bounds( values , std::move( subset ) , ordered , true ,
                       issuePMod , issueAMod );
 }

/*--------------------------------------------------------------------------*/

void UCBlock::set_pollutant_min_budget( MF_dbl_it values , Block::Range rng ,
                                        c_ModParam issuePMod ,
                                        c_ModParam issueAMod )
{
 set_pollutant_bounds( values , rng , true , issuePMod , issueAMod );
 }

/*--------------------------------------------------------------------------*/

void UCBlock::size_deferred( const std::function< void( void ) > & size )
{
 v_scaled_units.clear();
 f_deferred_injection_bounds = false;
 f_defer_scaled_rows = true;

 const auto rewrite = [ this ]() {
  f_defer_scaled_rows = false;
  update_scaled_rows( v_scaled_units );
  if( f_deferred_injection_bounds )
   set_node_injection_bounds();
  };

 try {
  size();
  }
 catch( ... ) {
  rewrite();
  throw;
  }

 rewrite();
 }

/*--------------------------------------------------------------------------*/

void UCBlock::resize_unit( MF_dbl_sp values , Subset && subset , bool ,
                           ModParam issuePMod , ModParam issueAMod )
{
 // each unit resized would have the bounds of the node injections computed
 // again [see add_Modification()]: they are computed once, after the last
 size_deferred( [ & ]() {
  size_units( *this , "resize_unit" , "::resize" , values , subset ,
              issuePMod , issueAMod );
  } );
 }

/*--------------------------------------------------------------------------*/

void UCBlock::resize_unit( MF_dbl_sp values , Range rng ,
                           ModParam issuePMod , ModParam issueAMod )
{
 resize_unit( values , units_in( *this , rng ) , true , issuePMod ,
              issueAMod );
 }

/*--------------------------------------------------------------------------*/

void UCBlock::resize_unit( MF_dbl_sp values , MF_dbl_sp second_values ,
                           Subset && subset , bool , ModParam issuePMod ,
                           ModParam issueAMod )
{
 if( ( subset.size() > values.size() ) ||
     ( subset.size() > second_values.size() ) )
  throw( std::invalid_argument( "UCBlock::resize_unit: a span is shorter "
                                "than the subset it is asked to write" ) );

 ClassMethod< FunctionType< double , double > >
  method( "resize_unit" , "::resize" , " with a pair of scalars" );

 // the bounds of the node injections are computed once, as in the other
 // resize_unit()
 size_deferred( [ & ]() {
  for( Index k = 0 ; k < subset.size() ; ++k ) {
   const auto u = subset[ k ];
   if( u >= f_number_units )
    throw( std::invalid_argument( "UCBlock::resize_unit: there is no unit "
                                  "of index " + std::to_string( u ) ) );
   const auto unit = get_unit_block( u );
   std::invoke( method( unit , "unit" , u ) , unit , values[ k ] ,
                second_values[ k ] , issuePMod , issueAMod );
   }
  } );
 }

/*--------------------------------------------------------------------------*/

void UCBlock::resize_unit( MF_dbl_sp values , MF_dbl_sp second_values ,
                           Range rng , ModParam issuePMod ,
                           ModParam issueAMod )
{
 resize_unit( values , second_values , units_in( *this , rng ) , true ,
              issuePMod , issueAMod );
 }

/*--------------------------------------------------------------------------*/

void UCBlock::replicate( MF_dbl_sp values , Subset && subset , bool ,
                         ModParam issuePMod , ModParam issueAMod )
{
 // each unit scaled would have its rows and the bounds of the node
 // injections rewritten [see add_Modification()]: they are rewritten once,
 // after the last
 size_deferred( [ & ]() {
  size_units( *this , "replicate" , "::replicate" , values , subset ,
              issuePMod , issueAMod );
  } );
 }

/*--------------------------------------------------------------------------*/

void UCBlock::replicate( MF_dbl_sp values , Range rng ,
                         ModParam issuePMod , ModParam issueAMod )
{
 replicate( values , units_in( *this , rng ) , true , issuePMod , issueAMod );
 }

/*--------------------------------------------------------------------------*/

/* The lines of a UCBlock are numbered carrier after carrier, those of the
 * carrier 0 first: first[ c ] is the index of the first line of carrier c,
 * first.back() the number of lines. A carrier with one node, hence no
 * NetworkBlock, has none. */

static std::vector< Block::Index > first_lines( const UCBlock & ucb )
{
 const auto nc = ucb.get_number_carriers();
 std::vector< Block::Index > first( nc + 1 , 0 );
 for( Block::Index c = 0 ; c < nc ; ++c ) {
  // the lines are those of the NetworkData the NetworkBlock share, which a
  // DesignNetworkBlock has too
  const auto & networks = ucb.get_network_blocks( c );
  const auto dc = networks.empty() ? nullptr :
                  dynamic_cast< const DCNetworkBlock::DCNetworkData * >(
                                         networks.front()->get_NetworkData() );
  first[ c + 1 ] = first[ c ] + ( dc ? dc->get_number_lines() : 0 );
  }
 return( first );
 }

/*--------------------------------------------------------------------------*/
/// the carrier of the line l of a UCBlock whose lines start at first

static Block::Index carrier_of( const std::vector< Block::Index > & first ,
                                Block::Index l , const std::string & who )
{
 const Block::Index c = std::upper_bound( first.begin() , first.end() , l ) -
                        first.begin() - 1;
 if( c + 1 >= first.size() )
  throw( std::invalid_argument( "UCBlock::" + who + ": there is no line of "
                                "index " + std::to_string( l ) ) );
 return( c );
 }

/*--------------------------------------------------------------------------*/

void UCBlock::resize_line( MF_dbl_sp values , Subset && subset , bool ordered ,
                           ModParam issuePMod , ModParam issueAMod )
{
 const auto first = first_lines( *this );
 if( ! first.back() )
  throw( std::logic_error( "UCBlock::resize_line: this UCBlock has no "
                           "lines" ) );
 if( subset.size() > values.size() )
  throw( std::invalid_argument( "UCBlock::resize_line: the span is shorter "
                                "than the lines it is asked to write" ) );

 ClassMethod< FunctionType< MF_dbl_sp , Subset && , bool > >
  resize( "resize_line" , "::resize" );

 // the lines of each carrier, numbered within it, and their values
 const auto nc = get_number_carriers();
 std::vector< Subset > lines( nc );
 std::vector< std::vector< double > > kappa( nc );
 for( Index k = 0 ; k < subset.size() ; ++k ) {
  const auto c = carrier_of( first , subset[ k ] , "resize_line" );
  lines[ c ].push_back( subset[ k ] - first[ c ] );
  kappa[ c ].push_back( values[ k ] );
  }

 for( Index c = 0 ; c < nc ; ++c ) {
  if( kappa[ c ].empty() )
   continue;
  // each NetworkBlock takes a Subset of its own, which the call consumes
  const auto & networks = get_network_blocks( c );
  for( Index n = 0 ; n < networks.size() ; ++n )
   std::invoke( resize( networks[ n ] , "NetworkBlock" , n ) , networks[ n ] ,
                MF_dbl_sp( kappa[ c ] ) , Subset( lines[ c ] ) , ordered ,
                issuePMod , issueAMod );
  }
 }

/*--------------------------------------------------------------------------*/

void UCBlock::resize_line( MF_dbl_sp values , Range rng ,
                           ModParam issuePMod , ModParam issueAMod )
{
 const auto first = first_lines( *this );
 if( ! first.back() )
  throw( std::logic_error( "UCBlock::resize_line: this UCBlock has no "
                           "lines" ) );
 const auto end = std::min( rng.second , first.back() );
 if( ( end > rng.first ) && ( end - rng.first > values.size() ) )
  throw( std::invalid_argument( "UCBlock::resize_line: the span is shorter "
                                "than the lines it is asked to write" ) );

 ClassMethod< FunctionType< MF_dbl_sp , Range > >
  resize( "resize_line" , "::resize" );

 // the part of the range on each carrier, numbered within it
 for( Index c = 0 ; c + 1 < first.size() ; ++c ) {
  const auto lo = std::max( rng.first , first[ c ] );
  const auto hi = std::min( rng.second , first[ c + 1 ] );
  if( lo >= hi )
   continue;
  const auto & networks = get_network_blocks( c );
  for( Index n = 0 ; n < networks.size() ; ++n )
   std::invoke( resize( networks[ n ] , "NetworkBlock" , n ) , networks[ n ] ,
                values.subspan( lo - rng.first , hi - lo ) ,
                Range( lo - first[ c ] , hi - first[ c ] ) , issuePMod ,
                issueAMod );
  }
 }

/*--------------------------------------------------------------------------*/

void UCBlock::get_resize_unit_linearization( MF_dbl_msp msp ,
                                             c_Subset & units ) const
{
 if( units.size() > msp.size() )
  throw( std::invalid_argument( "UCBlock::get_resize_unit_linearization: the "
                                "span is shorter than the subset it is asked "
                                "to answer for" ) );

 ClassMethod< QueryType< double & > >
  query( "get_resize_unit_linearization" , "::get_resize_linearization" );

 for( Index k = 0 ; k < units.size() ; ++k ) {
  if( units[ k ] >= f_number_units )
   throw( std::invalid_argument( "UCBlock::get_resize_unit_linearization: "
                                 "there is no unit of index " +
                                 std::to_string( units[ k ] ) ) );
  const auto unit = get_unit_block( units[ k ] );
  std::invoke( query( unit , "unit" , units[ k ] ) , unit , msp[ k ] );
  // the bounds of the node injections, which the unit does not see
  msp[ k ] += injection_bounds_linearization( units[ k ] , false );
  }
 }

/*--------------------------------------------------------------------------*/

void UCBlock::get_resize_unit_linearization( MF_dbl_msp msp ,
                                             MF_dbl_msp second_msp ,
                                             c_Subset & units ) const
{
 if( ( units.size() > msp.size() ) || ( units.size() > second_msp.size() ) )
  throw( std::invalid_argument( "UCBlock::get_resize_unit_linearization: a "
                                "span is shorter than the subset it is asked "
                                "to answer for" ) );

 ClassMethod< QueryType< double & , double & > >
  query( "get_resize_unit_linearization" , "::get_resize_linearization" ,
         " with a pair of scalars" );

 for( Index k = 0 ; k < units.size() ; ++k ) {
  if( units[ k ] >= f_number_units )
   throw( std::invalid_argument( "UCBlock::get_resize_unit_linearization: "
                                 "there is no unit of index " +
                                 std::to_string( units[ k ] ) ) );
  const auto unit = get_unit_block( units[ k ] );
  std::invoke( query( unit , "unit" , units[ k ] ) , unit , msp[ k ] ,
               second_msp[ k ] );
  // the bounds of the node injections, which only the first size enters
  msp[ k ] += injection_bounds_linearization( units[ k ] , false );
  }
 }

/*--------------------------------------------------------------------------*/

void UCBlock::get_resize_line_linearization( MF_dbl_msp msp ,
                                             c_Subset & lines ,
                                             bool ordered ) const
{
 if( lines.size() > msp.size() )
  throw( std::invalid_argument( "UCBlock::get_resize_line_linearization: the "
                                "span is shorter than the subset it is asked "
                                "to answer for" ) );
 const auto first = first_lines( *this );
 if( ! first.back() )
  throw( std::logic_error( "UCBlock::get_resize_line_linearization: this "
                           "UCBlock has no lines" ) );

 ClassMethod< QueryType< MF_dbl_msp , c_Subset & , bool > >
  query( "get_resize_line_linearization" , "::get_resize_linearization" );

 // the lines of each carrier, numbered within it, and where their answers go
 const auto nc = get_number_carriers();
 std::vector< Subset > local( nc );
 std::vector< std::vector< Index > > where( nc );
 for( Index k = 0 ; k < lines.size() ; ++k ) {
  const auto c = carrier_of( first , lines[ k ] ,
                             "get_resize_line_linearization" );
  local[ c ].push_back( lines[ k ] - first[ c ] );
  where[ c ].push_back( k );
  }

 // the answer for the lines of a carrier is the sum of those of its
 // NetworkBlock
 std::vector< double > sum;
 std::vector< double > answer;
 for( Index c = 0 ; c < nc ; ++c ) {
  if( local[ c ].empty() )
   continue;
  sum.assign( local[ c ].size() , 0.0 );
  answer.resize( local[ c ].size() );
  const auto & networks = get_network_blocks( c );
  for( Index n = 0 ; n < networks.size() ; ++n ) {
   std::invoke( query( networks[ n ] , "NetworkBlock" , n ) , networks[ n ] ,
                MF_dbl_msp( answer ) , local[ c ] , ordered );
   for( Index k = 0 ; k < answer.size() ; ++k )
    sum[ k ] += answer[ k ];
   }
  for( Index k = 0 ; k < sum.size() ; ++k )
   msp[ where[ c ][ k ] ] = sum[ k ];
  }
 }

/*--------------------------------------------------------------------------*/

void UCBlock::get_resize_line_linearization( MF_dbl_msp msp , Range rng )
 const
{
 const auto first = first_lines( *this );
 if( ! first.back() )
  throw( std::logic_error( "UCBlock::get_resize_line_linearization: this "
                           "UCBlock has no lines" ) );

 ClassMethod< QueryType< MF_dbl_msp , Range > >
  query( "get_resize_line_linearization" , "::get_resize_linearization" );

 /* The range stops at the number of lines, and the entries of msp past it
  * stay as the caller left them; the answer for the lines of a carrier is
  * the sum of those of its NetworkBlock. */
 rng.second = std::min( rng.second , first.back() );
 if( rng.second <= rng.first )
  return;
 if( msp.size() < rng.second - rng.first )
  throw( std::invalid_argument( "UCBlock::get_resize_line_linearization: the "
                                "span is shorter than the range it is asked "
                                "to answer for" ) );

 std::vector< double > answer;
 for( Index c = 0 ; c + 1 < first.size() ; ++c ) {
  const auto lo = std::max( rng.first , first[ c ] );
  const auto hi = std::min( rng.second , first[ c + 1 ] );
  if( lo >= hi )
   continue;
  const auto out = msp.subspan( lo - rng.first , hi - lo );
  std::fill( out.begin() , out.end() , 0.0 );
  answer.resize( hi - lo );
  const auto & networks = get_network_blocks( c );
  for( Index n = 0 ; n < networks.size() ; ++n ) {
   std::invoke( query( networks[ n ] , "NetworkBlock" , n ) , networks[ n ] ,
                MF_dbl_msp( answer ) ,
                Range( lo - first[ c ] , hi - first[ c ] ) );
   for( Index k = 0 ; k < answer.size() ; ++k )
    out[ k ] += answer[ k ];
   }
  }
 }

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS OF UCBlockSolution -------------------------*/
/*--------------------------------------------------------------------------*/

void UCBlockSolution::deserialize( const netCDF::NcGroup & group )
{
 // "TimeHorizon" is mandatory- - - - - - - - - - - - - - - - - - - - - - - -
 deserialize_dim( group , "TimeHorizon" , f_time_horizon , false );

 // the sub-groups are listed once: NcGroup::getGroup() lists them all at
 // each call, and there is one per NetworkBlock in the standard format
 const auto sub_groups = group.getGroups();
 auto sub_group = [ & sub_groups ]( const std::string & name ) {
  auto it = sub_groups.find( name );
  return( it == sub_groups.end() ? netCDF::NcGroup() : it->second );
  };

 // deserialize the UnitBlockSolution - - - - - - - - - - - - - - - - - - - -
 Index number_units = 0;
 if( deserialize_dim( group , "NumberUnits" , number_units ) ) {
  v_unit_Solution.resize( number_units );
  for( Index i = 0 ; i < number_units ; ++i ) {
   std::string sub_group_name = "UnitBlock_" + std::to_string( i );
   auto USi = dynamic_cast< UnitBlockSolution * >(
		       Solution::new_Solution( sub_group( sub_group_name ) ) );
   if( ! USi )
    throw( std::invalid_argument( "UCBlockSolution::deserialize: UnitBlock_"
				  + std::to_string( i ) +
				  " not a valid UnitBlockSolution" ) );
   v_unit_Solution[ i ] = USi;
   }
  }

 // deserialize the NetworkBlockSolution- - - - - - - - - - - - - - - - - - -
 Index number_networks = 0;
 deserialize_dim( group , "NumberNetworks" , number_networks );
 if( number_networks ) {
  v_network_Solution.resize( number_networks );

  // differently handle the standard format from the compressed one
  auto nb_group = sub_group( "NetworkBlock" );
  f_compressed_network = ! nb_group.isNull();
  if( f_compressed_network )  // compressed format
   NetworkBlockSolution::deserialize( nb_group , v_network_Solution );
  else  // standard format
   for( Index i = 0 ; i < number_networks ; ++i ) {
    std::string sub_group_name = "NetworkBlock_" + std::to_string( i );
    auto NBi = sub_group( sub_group_name );
    if( NBi.isNull() ) {  // serialize() skips the missing NetworkBlock
     v_network_Solution[ i ] = nullptr;
     continue;
     }
    auto Si = Solution::new_Solution( NBi );
    auto NSi = dynamic_cast< NetworkBlockSolution * >( Si );
    if( ! NSi ) {
     delete Si;
     throw( std::invalid_argument( "UCBlockSolution::deserialize: " +
				   sub_group_name +
				   " not a valid NetworkBlockSolution" ) );
     }
    v_network_Solution[ i ] = NSi;
    }
  }

 // deserialize the ActivePowerDuals- - - - - - - - - - - - - - - - - - - - -
 f_number_nodes = 0;
 if( deserialize_dim( group , "NumberNodes" , f_number_nodes ) )
  ::deserialize< double , 2 >( group , "ActivePowerDuals" ,
                               { f_time_horizon , f_number_nodes } ,
                               v_demand_duals , false , true );

 // deserialize the PrimaryDuals- - - - - - - - - - - - - - - - - - - - - - -
 f_number_primary_zones = 0;
 if( deserialize_dim( group , "NumberPrimaryZones" ,
			f_number_primary_zones ) )
  ::deserialize< double , 2 >( group , "PrimaryDuals" ,
                               { f_time_horizon , f_number_primary_zones } ,
                               v_primary_duals , false , true );

 // deserialize the SecondaryDuals- - - - - - - - - - - - - - - - - - - - - -
 f_number_secondary_zones = 0;
 if( deserialize_dim( group , "NumberSecondaryZones" ,
			f_number_secondary_zones ) )
  ::deserialize< double , 2 >( group , "SecondaryDuals" ,
                               { f_time_horizon , f_number_secondary_zones } ,
                               v_secondary_duals , false , true );

 // deserialize the InertiaDuals- - - - - - - - - - - - - - - - - - - - - - -
 f_number_inertia_zones = 0;
 if( deserialize_dim( group , "NumberInertiaZones" ,
			f_number_inertia_zones ) )
  ::deserialize< double , 2 >( group , "InertiaDuals" ,
                               { f_time_horizon , f_number_inertia_zones } ,
                               v_inertia_duals , false , true );

 // deserialize the PollutantDuals- - - - - - - - - - - - - - - - - - - - - -
 f_total_number_pollutant_zones = 0;
 if( deserialize_dim( group , "TotalNumberPollutantZones" ,
			f_total_number_pollutant_zones ) )
  ::deserialize( group , "PollutantDuals" , f_total_number_pollutant_zones ,
                 v_pollutant_duals , false , false );


 // the carriers c > 0, each in its group "Carrier_c", whose dimensions are
 // its own and not those above [see deserialize_local_dim()] - - - - - - - -
 for( auto & cns : v_carrier_network_Solution )
  for( auto nbs : cns )
   delete nbs;
 v_carrier_network_Solution.clear();
 v_carrier_duals.clear();

 Index nc = 1;
 deserialize_local_dim( group , "NumberCarriers" , nc );
 if( nc > 1 ) {
  v_carrier_network_Solution.resize( nc );
  v_carrier_duals.resize( nc );
  for( Index c = 1 ; c < nc ; ++c ) {
   auto cg = sub_group( "Carrier_" + std::to_string( c ) );
   if( cg.isNull() )
    continue;

   Index cnw = 0;
   if( deserialize_local_dim( cg , "NumberNetworks" , cnw ) && cnw ) {
    auto & cns = v_carrier_network_Solution[ c ];
    cns.resize( cnw , nullptr );
    const auto cg_groups = cg.getGroups();
    auto it = cg_groups.find( "NetworkBlock" );
    if( it != cg_groups.end() )  // compressed format
     NetworkBlockSolution::deserialize( it->second , cns );
    else
     for( Index t = 0 ; t < cnw ; ++t ) {
      auto jt = cg_groups.find( "NetworkBlock_" + std::to_string( t ) );
      if( jt == cg_groups.end() )
       continue;
      auto Si = Solution::new_Solution( jt->second );
      auto NSi = dynamic_cast< NetworkBlockSolution * >( Si );
      if( ! NSi ) {
       delete Si;
       throw( std::invalid_argument( "UCBlockSolution::deserialize: "
				     "Carrier_" + std::to_string( c ) +
				     "/NetworkBlock_" + std::to_string( t ) +
				     " not a valid NetworkBlockSolution" ) );
       }
      cns[ t ] = NSi;
      }
    }

   Index cnn = 0;
   if( deserialize_local_dim( cg , "NumberNodes" , cnn ) )
    ::deserialize< double , 2 >( cg , "EnergyDuals" ,
				 { f_time_horizon , cnn } ,
				 v_carrier_duals[ c ] , false , true );
   }
  }
 }  // end( UCBlockSolution::deserialize )

/*--------------------------------------------------------------------------*/

void UCBlockSolution::read( const Block * block )
{
 auto UCB = dynamic_cast< const UCBlock * >( block );
 if( ! UCB )
  throw( std::invalid_argument(
			 "UCBlockSolution::read: block is not a UCBlock" ) );

 f_time_horizon = UCB->get_time_horizon();
 f_number_nodes = UCB->get_number_nodes();
 f_number_primary_zones = UCB->get_number_primary_zones();
 f_number_secondary_zones = UCB->get_number_secondary_zones();
 f_number_inertia_zones = UCB->get_number_inertia_zones();
 f_total_number_pollutant_zones = UCB->get_total_number_pollutant_zones();

 if( ! v_unit_Solution.empty() ) {
  // read the UnitBlockSolution - - - - - - - - - - - - - - - - - - - - - - -
  for( Index i = 0 ; i <  v_unit_Solution.size() ; ++i ) {
   if( v_unit_Solution[ i ] )
    delete v_unit_Solution[ i ];
   v_unit_Solution[ i ] = static_cast< UnitBlockSolution * >(
		 UCB->get_unit_block( i )->get_Solution( nullptr , false ) );
   }
  }

 if( ! v_network_Solution.empty() ) {
  // read the NetworkBlockSolution- - - - - - - - - - - - - - - - - - - - - -
  for( Index i = 0 ; i < v_network_Solution.size() ; ++i ) {
   if( v_network_Solution[ i ] )
    delete v_network_Solution[ i ];
   if( auto NBi = UCB->get_network_block( i ) )
    v_network_Solution[ i ] = static_cast< NetworkBlockSolution * >(
		                     NBi->get_Solution( nullptr , false ) );
   else
    v_network_Solution[ i ] = nullptr;
   }
  }

 if( ! v_demand_duals.empty() ) {
  // read the dual variables of the node injection constraints- - - - - - - -
  auto & NIC = UCB->get_const_node_injection_constraints();
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   for( Index i = 0 ; i < f_number_nodes ; ++i )
    v_demand_duals[ t ][ i ] = NIC[ t ][ i ].get_dual();
  }

 if( ! v_primary_duals.empty() ) {
  // read the dual variables of the primary demand constraints- - - - - - - -
  auto & PDC = UCB->get_const_primary_demand_constraints();
  if( PDC.empty() )
   throw( std::invalid_argument(
        "UCBlockSolution::read-ing duals of non-existent primary demand" ) );

  for( Index t = 0 ; t < f_time_horizon  ; ++t )
   for( Index i = 0 ; i < f_number_primary_zones ; ++i )
    v_primary_duals[ t ][ i ] = PDC[ t ][ i ].get_dual();
  }

 if( ! v_secondary_duals.empty() ) {
  // read the dual variables of the secondary demand constraints- - - - - - -
  auto & SDC = UCB->get_const_secondary_demand_constraints();
  if( SDC.empty() )
   throw( std::invalid_argument(
      "UCBlockSolution::read-ing duals of non-existent secondary demand" ) );

  for( Index t = 0 ; t < f_time_horizon ; ++t )
   for( Index i = 0 ; i < f_number_secondary_zones ; ++i )
    v_secondary_duals[ t ][ i ] = SDC[ t ][ i ].get_dual();
  }

 if( ! v_inertia_duals.empty() ) {
  // read the dual variables of the inertia demand constraints - - - - - - -
  auto & IDC = UCB->get_const_inertia_demand_constraints();
  if( IDC.empty() )
   throw( std::invalid_argument(
      "UCBlockSolution::read-ing duals of non-existent secondary demand" ) );

  for( Index t = 0 ; t < f_time_horizon  ; ++t )
   for( Index i = 0 ; i < f_number_inertia_zones ; ++i )
    v_inertia_duals[ t ][ i ] = IDC[ t ][ i ].get_dual();
  }

 if( ! v_pollutant_duals.empty() ) {
  // read the dual variables of the pollutant budget constraints- - - - - - -
  auto & PBC = UCB->get_const_pollutant_constraints();
  if( PBC.empty() )
   throw( std::invalid_argument(
       "UCBlockSolution::read-ing duals of non-existent pollutant budget" ) );

  v_pollutant_duals.clear();
  for( const auto & zones : PBC )
   for( const auto & constraint : zones )
    v_pollutant_duals.push_back( constraint.get_dual() );
  }

 // the NetworkBlock and the balances (1c) of the carriers c > 0- - - - - - - -
 for( Index c = 1 ; c < v_carrier_network_Solution.size() ; ++c ) {
  const auto & nbs = UCB->get_network_blocks( c );
  auto & cns = v_carrier_network_Solution[ c ];
  for( Index t = 0 ; t < cns.size() ; ++t ) {
   delete cns[ t ];
   cns[ t ] = static_cast< NetworkBlockSolution * >(
		                nbs[ t ]->get_Solution( nullptr , false ) );
   }
  }

 for( Index c = 1 ; c < v_carrier_duals.size() ; ++c ) {
  const auto & rows = UCB->get_const_carrier_balance_constraints( c );
  if( rows.empty() )
   throw( std::invalid_argument( "UCBlockSolution::read-ing duals of "
				 "non-existent carrier balances" ) );
  auto & duals = v_carrier_duals[ c ];
  for( Index t = 0 ; t < duals.shape()[ 0 ] ; ++t )
   for( Index n = 0 ; n < duals.shape()[ 1 ] ; ++n )
    duals[ t ][ n ] = rows[ t ][ n ].get_dual();
  }
 }  // end( UCBlockSolution::read )

/*--------------------------------------------------------------------------*/

void UCBlockSolution::write( Block * block )
{
 auto UCB = dynamic_cast< UCBlock * >( block );
 if( ! UCB )
  throw( std::invalid_argument(
			"UCBlockSolution::write: block is not a UCBlock" ) );

 if( f_time_horizon != UCB->get_time_horizon() )
  throw( std::invalid_argument(
		     "UCBlockSolution::write: inconsistent time horizon" ) );

 if( ! v_unit_Solution.empty() ) {
  // write the UnitBlockSolution- - - - - - - - - - - - - - - - - - - - - - -
  for( Index i = 0 ; i <  v_unit_Solution.size() ; ++i )
   v_unit_Solution[ i ]->write( UCB->get_unit_block( i ) );
  }

 if( ! v_network_Solution.empty() ) {
  // write the NetworkBlockSolution - - - - - - - - - - - - - - - - - - - - -
  for( Index i = 0 ; i < v_network_Solution.size() ; ++i )
   if( v_network_Solution[ i ] )
    if( auto NBi = UCB->get_network_block( i ) )
     v_network_Solution[ i ]->write( NBi );
  }

 if( ! v_demand_duals.empty() ) {
  // write the dual variables of the node injection constraints - - - - - - -
  if( f_number_nodes != UCB->get_number_nodes() )
   throw( std::invalid_argument(
		      "UCBlockSolution::write: inconsistent node number" ) );
  auto & NIC = UCB->get_node_injection_constraints();
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   for( Index i = 0 ; i < f_number_nodes ; ++i )
    NIC[ t ][ i ].set_dual( v_demand_duals[ t ][ i ] );
  }

 if( ! v_primary_duals.empty() ) {
  // write the dual variables of the primary demand constraints - - - - - - -
  if( f_number_primary_zones != UCB->get_number_primary_zones() )
   throw( std::invalid_argument(
		    "UCBlockSolution::write: inconsistent primary zones" ) );
  auto & PDC = UCB->get_primary_demand_constraints();
  if( PDC.empty() )
   throw( std::invalid_argument(
       "UCBlockSolution::write-ing duals of non-existent primary demand" ) );

  for( Index t = 0 ; t < f_time_horizon  ; ++t )
   for( Index i = 0 ; i < f_number_primary_zones ; ++i )
    PDC[ t ][ i ].set_dual( v_primary_duals[ t ][ i ] );
  }

 if( ! v_secondary_duals.empty() ) {
  // write the dual variables of the secondary demand constraints - - - - - -
  if( f_number_secondary_zones != UCB->get_number_secondary_zones() )
   throw( std::invalid_argument(
		  "UCBlockSolution::write: inconsistent secondary zones" ) );
  auto & SDC = UCB->get_secondary_demand_constraints();
  if( SDC.empty() )
   throw( std::invalid_argument(
     "UCBlockSolution::write-ing duals of non-existent secondary demand" ) );

  for( Index t = 0 ; t < f_time_horizon  ; ++t )
   for( Index i = 0 ; i < f_number_secondary_zones ; ++i )
    SDC[ t ][ i ].set_dual( v_secondary_duals[ t ][ i ] );
  }

 if( ! v_inertia_duals.empty() ) {
  // write the dual variables of the inertia demand constraints- - - - - - -
  if( f_number_inertia_zones != UCB->get_number_inertia_zones() )
   throw( std::invalid_argument(
		    "UCBlockSolution::write: inconsistent inertia zones" ) );
  auto & IDC = UCB->get_inertia_demand_constraints();
  if( IDC.empty() )
   throw( std::invalid_argument(
     "UCBlockSolution::write-ing duals of non-existent secondary demand" ) );

  for( Index t = 0 ; t < f_time_horizon  ; ++t )
   for( Index i = 0 ; i < f_number_inertia_zones ; ++i )
    IDC[ t ][ i ].set_dual( v_inertia_duals[ t ][ i ] );
  }

 if( ! v_pollutant_duals.empty() ) {
  // write the dual variables of the pollutant budget constraints - - - - - -
  if( f_total_number_pollutant_zones !=
      UCB->get_total_number_pollutant_zones() )
   throw( std::invalid_argument(
		   "UCBlockSolution::write: inconsistent pollutant zones" ) );
  auto & PBC = UCB->get_pollutant_constraints();
  if( PBC.empty() )
   throw( std::invalid_argument(
     "UCBlockSolution::write-ing duals of non-existent pollutant budget" ) );

  auto dual = v_pollutant_duals.begin();
  for( auto & zones : PBC )
   for( auto & constraint : zones )
    constraint.set_dual( *(dual++) );
  }

 // the NetworkBlock and the balances (1c) of the carriers c > 0- - - - - - - -
 if( std::max( v_carrier_network_Solution.size() ,
	       v_carrier_duals.size() ) > UCB->get_number_carriers() )
  throw( std::invalid_argument(
		      "UCBlockSolution::write: inconsistent carriers" ) );

 // a carrier whose networks or duals the Solution does not have is skipped
 for( Index c = 1 ; c < v_carrier_network_Solution.size() ; ++c ) {
  const auto & nbs = UCB->get_network_blocks( c );
  const auto & cns = v_carrier_network_Solution[ c ];
  if( cns.empty() )
   continue;
  if( cns.size() != nbs.size() )
   throw( std::invalid_argument(
	       "UCBlockSolution::write: inconsistent carrier networks" ) );
  for( Index t = 0 ; t < cns.size() ; ++t )
   if( cns[ t ] )
    cns[ t ]->write( nbs[ t ] );
  }

 for( Index c = 1 ; c < v_carrier_duals.size() ; ++c ) {
  const auto & duals = v_carrier_duals[ c ];
  if( ! duals.num_elements() )
   continue;
  if( duals.shape()[ 1 ] != UCB->get_number_nodes( c ) )
   throw( std::invalid_argument(
	      "UCBlockSolution::write: inconsistent carrier node number" ) );
  auto & rows = UCB->get_carrier_balance_constraints( c );
  if( rows.empty() )
   throw( std::invalid_argument( "UCBlockSolution::write-ing duals of "
				 "non-existent carrier balances" ) );
  for( Index t = 0 ; t < duals.shape()[ 0 ] ; ++t )
   for( Index n = 0 ; n < duals.shape()[ 1 ] ; ++n )
    rows[ t ][ n ].set_dual( duals[ t ][ n ] );
  }
 }  // end( UCBlockSolution::write )

/*--------------------------------------------------------------------------*/

void UCBlockSolution::serialize( netCDF::NcGroup & group ) const
{
 Solution::serialize( group );

 // "TimeHorizon" is mandatory- - - - - - - - - - - - - - - - - - - - - - - -
 auto th = group.addDim( "TimeHorizon" , f_time_horizon );

 // serialize the UnitBlockSolution - - - - - - - - - - - - - - - - - - - - -
 if( ! v_unit_Solution.empty() ) {
  auto nu = group.addDim( "NumberUnits" , v_unit_Solution.size() );

  std::vector< std::string > names( v_unit_Solution.size() );
  for( Index i = 0 ; i < names.size() ; ++i )
   names[ i ] = "UnitBlock_" + std::to_string( i );
  write_via_memory( group , names ,
		    [ & ]( std::size_t i , netCDF::NcGroup & g ) {
		     v_unit_Solution[ i ]->serialize( g ); } );
  }

 // serialize the NetworkBlockSolution- - - - - - - - - - - - - - - - - - - -
 if( ! v_network_Solution.empty() ) {
  group.addDim( "NumberNetworks" , v_network_Solution.size() );

  if( f_compressed_network ) {  // compressed format
   auto sub_group = group.addGroup( "NetworkBlock" );
   NetworkBlockSolution::serialize( sub_group , v_network_Solution );
   }
  else {  // standard format
   std::vector< Index > which;
   std::vector< std::string > names;
   for( Index i = 0 ; i < v_network_Solution.size() ; ++i )
    if( v_network_Solution[ i ] ) {
     which.push_back( i );
     names.push_back( "NetworkBlock_" + std::to_string( i ) );
     }
   write_via_memory( group , names ,
		     [ & ]( std::size_t k , netCDF::NcGroup & g ) {
		      v_network_Solution[ which[ k ] ]->serialize( g ); } );
   }
  }

 // serialize the ActivePowerDuals- - - - - - - - - - - - - - - - - - - - - -
 if( ! v_demand_duals.empty() ) {
  auto nn = group.addDim( "NumberNodes" , v_demand_duals.shape()[ 1 ] );

  ::serialize< double , 2 >( group , "ActivePowerDuals" , netCDF::NcDouble() ,
			     { th , nn } , v_demand_duals );
  }

 // serialize the PrimaryDuals- - - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_primary_duals.empty() ) {
  auto npz = group.addDim( "NumberPrimaryZones" , f_number_primary_zones );

  ::serialize< double , 2 >( group , "PrimaryDuals" , netCDF::NcDouble() ,
			     { th , npz } , v_primary_duals );
  }

 // serialize the SecondaryDuals- - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_secondary_duals.empty() ) {
  auto nsz = group.addDim( "NumberSecondaryZones" ,
			   f_number_secondary_zones );

  ::serialize< double , 2 >( group , "SecondaryDuals" , netCDF::NcDouble() ,
			     { th , nsz } , v_secondary_duals );
  }

 // serialize the InertiaDuals- - - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_inertia_duals.empty() ) {
  auto niz = group.addDim( "NumberInertiaZones" , f_number_inertia_zones );

  ::serialize< double , 2 >( group , "InertiaDuals" , netCDF::NcDouble() ,
			     { th , niz } , v_inertia_duals );
  }

 // serialize the PollutantDuals- - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_pollutant_duals.empty() ) {
  auto tnpz = group.addDim( "TotalNumberPollutantZones" ,
			    f_total_number_pollutant_zones );

  ::serialize( group , "PollutantDuals" , netCDF::NcDouble() , tnpz ,
	       v_pollutant_duals );
  }

 // the carriers c > 0, each in its group "Carrier_c"- - - - - - - - - - - -
 const auto nc = std::max( v_carrier_network_Solution.size() ,
			   v_carrier_duals.size() );
 if( nc > 1 ) {
  group.addDim( "NumberCarriers" , nc );
  for( Index c = 1 ; c < nc ; ++c ) {
   auto cg = group.addGroup( "Carrier_" + std::to_string( c ) );
   auto cth = cg.addDim( "TimeHorizon" , f_time_horizon );

   if( ( c < v_carrier_network_Solution.size() ) &&
       ( ! v_carrier_network_Solution[ c ].empty() ) ) {
    const auto & cns = v_carrier_network_Solution[ c ];
    cg.addDim( "NumberNetworks" , cns.size() );
    if( f_compressed_network ) {
     auto sub_group = cg.addGroup( "NetworkBlock" );
     NetworkBlockSolution::serialize( sub_group , cns );
     }
    else
     for( Index t = 0 ; t < cns.size() ; ++t )
      if( cns[ t ] ) {
       auto sub_group = cg.addGroup( "NetworkBlock_" + std::to_string( t ) );
       cns[ t ]->serialize( sub_group );
       }
    }

   if( ( c < v_carrier_duals.size() ) && ( ! v_carrier_duals[ c ].empty() ) ) {
    auto cnn = cg.addDim( "NumberNodes" , v_carrier_duals[ c ].shape()[ 1 ] );
    ::serialize< double , 2 >( cg , "EnergyDuals" , netCDF::NcDouble() ,
			       { cth , cnn } , v_carrier_duals[ c ] );
    }
   }
  }
 }  // end( UCBlockSolution::serialize )

/*--------------------------------------------------------------------------*/

UCBlockSolution * UCBlockSolution::scale( double factor ) const
{
 auto sol = clone();

 if( factor == 1 )
  return( sol );

 // replace each unit and network Solution of the clone with its scaled
 // version
 for( auto & vi : sol->v_unit_Solution ) {
  auto unit_solution = vi;
  vi = unit_solution->scale( factor );
  delete unit_solution;
  }

 for( auto & ni : sol->v_network_Solution )
  if( ni ) {
   auto network_solution = ni;
   ni = network_solution->scale( factor );
   delete network_solution;
   }

 if( ! v_demand_duals.empty() )
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   for( Index i = 0 ; i < f_number_nodes ; ++i )
    sol->v_demand_duals[ t ][ i ] *= factor;

 if( ! v_primary_duals.empty() )
  for( Index t = 0 ; t < f_time_horizon  ; ++t )
   for( Index i = 0 ; i < f_number_primary_zones ; ++i )
    sol->v_primary_duals[ t ][ i ] *= factor;

 if( ! v_secondary_duals.empty() )
  for( Index t = 0 ; t < f_time_horizon  ; ++t )
   for( Index i = 0 ; i < f_number_secondary_zones ; ++i )
    sol->v_secondary_duals[ t ][ i ] *= factor;

 if( ! v_inertia_duals.empty() )
  for( Index t = 0 ; t < f_time_horizon  ; ++t )
   for( Index i = 0 ; i < f_number_inertia_zones ; ++i )
    sol->v_inertia_duals[ t ][ i ] *= factor;

 for( auto & cns : sol->v_carrier_network_Solution )
  for( auto & ni : cns )
   if( ni ) {
    auto network_solution = ni;
    ni = network_solution->scale( factor );
    delete network_solution;
    }

 for( auto & duals : sol->v_carrier_duals )
  for( auto it = duals.data() ; it != duals.data() + duals.num_elements() ;
       ++it )
   *it *= factor;

 for( auto & di : sol->v_pollutant_duals )
  di *= factor;

 return( sol );

 }  // end( UCBlockSolution::scale )

/*--------------------------------------------------------------------------*/

void UCBlockSolution::sum( const Solution * solution , double multiplier )
{
 auto UCBS = dynamic_cast< const UCBlockSolution * >( solution );
 if( ! UCBS )
  throw( std::invalid_argument(
	       "UCBlockSolution::sum: solution is not a UCBlockSolution" ) );

 if( f_time_horizon != UCBS->f_time_horizon )
  throw( std::invalid_argument(
		       "UCBlockSolution::sum: inconsistent time horizon" ) );
 if( v_unit_Solution.size() != UCBS->v_unit_Solution.size() )
  throw( std::invalid_argument(
		     "UCBlockSolution::read: inconsistent unit solution" ) );
 if( v_network_Solution.size() != UCBS->v_network_Solution.size() )
  throw( std::invalid_argument(
		  "UCBlockSolution::read: inconsistent network solution" ) );
 if( f_number_nodes != UCBS->f_number_nodes )
  throw( std::invalid_argument(
		        "UCBlockSolution::sum: inconsistent node number" ) );
 if( f_number_primary_zones != UCBS->f_number_primary_zones )
  throw( std::invalid_argument(
		      "UCBlockSolution::sum: inconsistent primary zones" ) );
 if( f_number_secondary_zones != UCBS->f_number_secondary_zones )
  throw( std::invalid_argument(
		    "UCBlockSolution::sum: inconsistent secondary zones" ) );
 if( f_number_inertia_zones != UCBS->f_number_inertia_zones )
  throw( std::invalid_argument(
		     "UCBlockSolution::read: inconsistent inertia zones" ) );
 if( v_pollutant_duals.size() != UCBS->v_pollutant_duals.size() )
  throw( std::invalid_argument(
		   "UCBlockSolution::sum: inconsistent pollutant zones" ) );

 for( Index i = 0 ; i < v_unit_Solution.size() ; ++i )
  v_unit_Solution[ i ]->sum( UCBS->v_unit_Solution[ i ] , multiplier );

 for( Index i = 0 ; i < v_network_Solution.size() ; ++i )
  v_network_Solution[ i ]->sum( UCBS->v_network_Solution[ i ] , multiplier );

 if( ! v_demand_duals.empty() )
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   for( Index i = 0 ; i < f_number_nodes ; ++i )
    v_demand_duals[ t ][ i ] += UCBS->v_demand_duals[ t ][ i ] * multiplier;

 if( ! v_primary_duals.empty() )
  for( Index t = 0 ; t < f_time_horizon  ; ++t )
   for( Index i = 0 ; i < f_number_primary_zones ; ++i )
    v_primary_duals[ t ][ i ] += UCBS->v_primary_duals[ t ][ i ]  * multiplier;

 if( ! v_secondary_duals.empty() )
  for( Index t = 0 ; t < f_time_horizon  ; ++t )
   for( Index i = 0 ; i < f_number_secondary_zones ; ++i )
    v_secondary_duals[ t ][ i ] +=
     UCBS->v_secondary_duals[ t ][ i ] * multiplier;

 if( ! v_inertia_duals.empty() )
  for( Index t = 0 ; t < f_time_horizon  ; ++t )
   for( Index i = 0 ; i < f_number_inertia_zones ; ++i )
    v_inertia_duals[ t ][ i ] += UCBS->v_inertia_duals[ t ][ i ] * multiplier;

 for( Index i = 0 ; i < v_pollutant_duals.size() ; ++i )
  v_pollutant_duals[ i ] += UCBS->v_pollutant_duals[ i ] * multiplier;

 if( ( v_carrier_network_Solution.size() !=
       UCBS->v_carrier_network_Solution.size() ) ||
     ( v_carrier_duals.size() != UCBS->v_carrier_duals.size() ) )
  throw( std::invalid_argument(
		         "UCBlockSolution::sum: inconsistent carriers" ) );

 for( Index c = 0 ; c < v_carrier_network_Solution.size() ; ++c ) {
  auto & cns = v_carrier_network_Solution[ c ];
  const auto & other = UCBS->v_carrier_network_Solution[ c ];
  if( cns.size() != other.size() )
   throw( std::invalid_argument(
	        "UCBlockSolution::sum: inconsistent carrier networks" ) );
  for( Index t = 0 ; t < cns.size() ; ++t )
   if( cns[ t ] && other[ t ] )
    cns[ t ]->sum( other[ t ] , multiplier );
  }

 for( Index c = 0 ; c < v_carrier_duals.size() ; ++c ) {
  auto & duals = v_carrier_duals[ c ];
  const auto & other = UCBS->v_carrier_duals[ c ];
  if( duals.num_elements() != other.num_elements() )
   throw( std::invalid_argument(
	     "UCBlockSolution::sum: inconsistent carrier node number" ) );
  for( Index i = 0 ; i < duals.num_elements() ; ++i )
   duals.data()[ i ] += other.data()[ i ] * multiplier;
  }

 }  // end( UCBlockSolution::sum )

/*--------------------------------------------------------------------------*/

UCBlockSolution * UCBlockSolution::clone( bool empty ) const
{
 auto sol = new UCBlockSolution();

 if( ! empty ) {
  sol->f_time_horizon = f_time_horizon;
  sol->f_number_nodes = f_number_nodes;
  sol->f_number_primary_zones = f_number_primary_zones;
  sol->f_number_secondary_zones = f_number_secondary_zones;
  sol->f_number_inertia_zones = f_number_inertia_zones;
  sol->f_total_number_pollutant_zones = f_total_number_pollutant_zones;

  if( ! v_unit_Solution.empty() ) {
   sol->v_unit_Solution.resize( v_unit_Solution.size() );
   for( Index i = 0 ; i < v_unit_Solution.size() ; ++i )
    sol->v_unit_Solution[ i ] = v_unit_Solution[ i ]->clone();
   }

  if( ! v_network_Solution.empty() ) {
   sol->v_network_Solution.resize( v_network_Solution.size() );
   for( Index i = 0 ; i < v_network_Solution.size() ; ++i )
    if( v_network_Solution[ i ] )
     sol->v_network_Solution[ i ] = v_network_Solution[ i ]->clone();
    else
     sol->v_network_Solution[ i ] = nullptr;
   }

  copy_multi_array( sol->v_demand_duals , v_demand_duals );
  copy_multi_array( sol->v_primary_duals , v_primary_duals );
  copy_multi_array( sol->v_secondary_duals , v_secondary_duals );
  copy_multi_array( sol->v_inertia_duals , v_inertia_duals );
  sol->v_pollutant_duals = v_pollutant_duals;

  sol->v_carrier_network_Solution.resize( v_carrier_network_Solution.size() );
  for( Index c = 0 ; c < v_carrier_network_Solution.size() ; ++c ) {
   const auto & cns = v_carrier_network_Solution[ c ];
   auto & sol_cns = sol->v_carrier_network_Solution[ c ];
   sol_cns.resize( cns.size() , nullptr );
   for( Index t = 0 ; t < cns.size() ; ++t )
    if( cns[ t ] )
     sol_cns[ t ] = cns[ t ]->clone();
   }

  sol->v_carrier_duals.resize( v_carrier_duals.size() );
  for( Index c = 0 ; c < v_carrier_duals.size() ; ++c )
   copy_multi_array( sol->v_carrier_duals[ c ] , v_carrier_duals[ c ] );
  }

 return( sol );

 }  // end( UCBlockSolution::clone )

/*--------------------------------------------------------------------------*/
/*------------------------ End File UCBlock.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
