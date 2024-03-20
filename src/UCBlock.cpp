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

// TODO commented away until HeatBlock are properly managed
// #include "HeatBlock.h"

#include "BlockInspection.h"

#include "LinearFunction.h"

#include "Objective.h"

#include "UCBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register UCBlock to the Block factory

SMSpp_insert_in_factory_cpp_1( UCBlock );

/*--------------------------------------------------------------------------*/
/*--------------------------- METHODS OF UCBlock ---------------------------*/
/*--------------------------------------------------------------------------*/

UCBlock::~UCBlock()
{
 Constraint::clear( v_node_injection_Const );
 Constraint::clear( v_PrimaryDemand_Const );
 Constraint::clear( v_SecondaryDemand_Const );
 Constraint::clear( v_InertiaDemand_Const );
 /* TODO commented away until HeatBlock are properly managed
   Constraint::clear( v_power_Heat_Rho_Const );
 */

 for( auto & v_constraints : v_PollutantBudget_Const )
  Constraint::clear( v_constraints );
 v_PollutantBudget_Const.clear();

 for( auto & block : v_Block )
  delete( block );

 delete( f_NetworkData );
}

/*--------------------------------------------------------------------------*/

void UCBlock::deserialize_sub_blocks( const netCDF::NcGroup & group ,
                                      const std::string & prefix ,
                                      Index num_sub_blocks )
{
 auto sz = v_Block.size();
 v_Block.resize( sz + num_sub_blocks );
 for( Index i = 0 ; i < num_sub_blocks ; ++i ) {
  std::string sub_group_name = prefix + std::to_string( i );
  auto sub_group = group.getGroup( sub_group_name );
  if( sub_group.isNull() )
   throw( std::invalid_argument( "UCBlock::deserialize: " +
                                 sub_group_name + " not present" ) );

  if( auto bk = new_Block( sub_group , this ) )
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
 v_network_blocks.resize( f_number_networks );

 for( Index i = 0 ; i < f_number_networks ; ++i ) {
  std::string sub_group_name = "NetworkBlock_" + std::to_string( i );
  auto sub_group = group.getGroup( sub_group_name );
  if( sub_group.isNull() )
   continue;

  if( auto nbi = dynamic_cast< NetworkBlock * >(
   new_Block( sub_group , this ) ) ) {
   v_network_blocks[ i ] = nbi;
   delete( f_NetworkData );
   f_NetworkData = static_cast< NetworkBlock::NetworkData * >(
    NetworkBlock::NetworkData::new_NetworkData( network_data_classname ) );
   f_NetworkData->deserialize( sub_group );
   nbi->set_NetworkData( f_NetworkData );
   ++cntr;
  } else {
   delete( nbi );
   throw( std::invalid_argument( "UCBlock::deserialize:" + sub_group_name +
                                 " not a valid NetworkBlock" ) );
  }
 }

 if( cntr ) {
  v_Block.resize( f_number_units + f_number_networks );
  std::copy( v_network_blocks.begin() , v_network_blocks.end() ,
             std::next( v_Block.begin() , f_number_units ) );
 } else
  v_network_blocks.clear();
}

/*--------------------------------------------------------------------------*/

void UCBlock::deserialize( const netCDF::NcGroup & group )
{

#ifndef NDEBUG
 static std::vector< std::string > expected_dims = { "TimeHorizon" ,
                                                     "NumberUnits" ,
                                                     "NumberNetworks" ,
                                                     "NumberHeatBlocks" ,
                                                     "NumberPrimaryZones" ,
                                                     "NumberSecondaryZones" ,
                                                     "NumberInertiaZones" ,
                                                     "NumberPollutants" ,
                                                     "NumberNodes" ,
                                                     "NumberLines" ,
                                                     "NumberElectricalGenerators" ,
                                                     "TotalNumberPollutantZones" ,
                                                     "NumberIntervals" };
 check_dimensions( group , expected_dims , std::cerr );

 static std::vector< std::string > expected_vars = { "ActivePowerDemand" ,
                                                     "GeneratorNode" ,
                                                     "HeatNode" , "HeatSet" ,
                                                     "PowerHeatRho" ,
                                                     "PrimaryZones" ,
                                                     "PrimaryDemand" ,
                                                     "SecondaryZones" ,
                                                     "SecondaryDemand" ,
                                                     "InertiaZones" ,
                                                     "InertiaDemand" ,
                                                     "NumberPollutantZones" ,
                                                     "PollutantZones" ,
                                                     "PollutantBudget" ,
                                                     "PollutantRho" ,
                                                     "StartNetworkIntervals" ,
                                                     "NetworkConstantTerms" ,
                                                     "NetworkBlockClassname" ,
                                                     "NetworkDataClassname" ,
                                                     // DCNetworkBlockData
                                                     "StartLine" , "EndLine" ,
                                                     "MinPowerFlow" ,
                                                     "MaxPowerFlow" ,
                                                     "Susceptance" ,
                                                     "NetworkCost" ,
                                                     // ECNetworkBlockData
                                                     "BuyPrice" ,
                                                     "SellPrice" ,
                                                     "RewardPrice" ,
                                                     "PeakTariff" };
 check_variables( group , expected_vars , std::cerr );
#endif

 // Mandatory variables

 ::deserialize_dim( group , "TimeHorizon" , f_time_horizon , false );
 ::deserialize_dim( group , "NumberUnits" , f_number_units , false );

 // Optional variables

 if( ! ::deserialize_dim( group , "NumberNetworks" , f_number_networks ) )
  f_number_networks = f_time_horizon;

 if( ! ::deserialize( group , "StartNetworkIntervals" , f_number_networks ,
                      v_start_network_intervals ) ) {
  v_start_network_intervals.resize( f_number_networks );
  std::iota( v_start_network_intervals.begin() ,
             v_start_network_intervals.end() , 0 );
 }
 v_start_network_intervals.push_back( f_time_horizon );

 if( ! ::deserialize( group , "NetworkConstantTerms" , f_number_networks ,
                      v_network_constant_terms ) )
  v_network_constant_terms.resize( f_number_networks );

 if( ! ::deserialize( group , network_block_classname ,
                      "NetworkBlockClassname" ) )
  network_block_classname = "DCNetworkBlock";
 if( ! ::deserialize( group , network_data_classname ,
                      "NetworkDataClassname" ) )
  network_data_classname = "DCNetworkData";

 Index number_nodes;
 if( ! ::deserialize_dim( group , "NumberNodes" , number_nodes ) )
  number_nodes = 1;
 if( ( network_block_classname == "ECNetworkBlock" ) &&
     ( number_nodes == 1 ) )
  throw( std::invalid_argument( "UCBlock::deserialize: cannot create "
                                "a community network with just one user" ) );

 /* TODO commented away until HeatBlock are properly managed
 if( ! ::deserialize_dim( group , "NumberHeatGenerators" ,
                          f_number_heat_generators ) )
  f_number_heat_generators = 0;
 */

 if( ::deserialize( group , "ActivePowerDemand" , v_active_power_demand ) ) {

  // when the network is a bus, the "ActivePowerDemand" variable in the nc4
  // input file could be provided as a simple 1D array of `f_time_horizon`
  // size, so we reshape `v_active_power_demand` in order to make it available
  // in the expected shape, i.e., `number_nodes` (= 1) x `f_time_horizon`
  if( ( number_nodes == 1 ) &&
      // ensure if is in fact provided as a 1D array, ignore if it is given
      // in the correct shape
      ( v_active_power_demand.shape()[ 0 ] == f_time_horizon ) ) {
   using index = decltype( v_active_power_demand )::index;
   std::vector< index > shape = { 1 , f_time_horizon };
   v_active_power_demand.reshape( shape );
  }

  // always check if the demand is given in the correct shape
  assert( ( v_active_power_demand.shape()[ 0 ] == number_nodes ) &&
          ( v_active_power_demand.shape()[ 1 ] == f_time_horizon ) );
 }

 // Optional dimensions

 /* TODO commented away until HeatBlock are properly managed
 f_number_heat_blocks = 0;
 ::deserialize_dim( group , "NumberHeatBlocks" , f_number_heat_blocks );
 */

 f_number_primary_zones = 0;
 ::deserialize_dim( group , "NumberPrimaryZones" , f_number_primary_zones );

 f_number_secondary_zones = 0;
 ::deserialize_dim( group , "NumberSecondaryZones" , f_number_secondary_zones );

 f_number_inertia_zones = 0;
 ::deserialize_dim( group , "NumberInertiaZones" , f_number_inertia_zones );

 f_number_pollutants = 0;
 ::deserialize_dim( group , "NumberPollutants" , f_number_pollutants );

 /* TODO commented away until HeatBlock are properly managed
 ::deserialize( group , "HeatSet" ,
                { f_number_units , f_number_heat_blocks } , v_heat_set ,
                true , false );
 */

 ::deserialize( group , "PrimaryZones" , number_nodes ,
                v_primary_zones , true , true );

 if( ::deserialize( group , "PrimaryDemand" ,
                    v_primary_demand , true , false ) )

  ::deserialize( group , "SecondaryZones" , number_nodes ,
                 v_secondary_zones , true , true );

 if( ::deserialize( group , "SecondaryDemand" ,
                    v_secondary_demand , true , false ) )

  ::deserialize( group , "InertiaZones" , number_nodes ,
                 v_inertia_zones , true , true );

 if( ::deserialize( group , "InertiaDemand" ,
                    v_inertia_demand , true , false ) )

  ::deserialize( group , "NumberPollutantZones" , f_number_pollutants ,
                 v_number_pollutant_zones , true , true );

 if( ! ::deserialize_dim( group , "TotalNumberPollutantZones" ,
                         f_total_number_pollutant_zones ) ) {
  f_total_number_pollutant_zones = 0;
  for( const auto & n : v_number_pollutant_zones )
   f_total_number_pollutant_zones += n;
 }

 if( ! f_total_number_pollutant_zones )
  f_total_number_pollutant_zones = f_number_pollutants;

 if( f_total_number_pollutant_zones ) {
  ::deserialize( group , "PollutantZones" , v_pollutant_zones , true , true );

  /* TODO commented away until this is properly managed
  ::deserialize( group , "PollutantBudget" ,
                 { f_total_number_pollutant_zones } ,
                 v_pollutant_budget , true , false );
  */

  ::deserialize( group , "PollutantRho" , v_pollutant_rho , true , true );

  /* TODO commented away until HeatBlock are properly managed
  ::deserialize( group , "PollutantHeatRho" ,
                 v_pollutant_heat_rho , true , true );
  */
 } else {
  v_pollutant_zones.resize(
   boost::multi_array< Index , 2 >::extent_gen()[ 0 ][ 0 ] );
  v_pollutant_budget.clear();
  v_pollutant_rho.resize(
   boost::multi_array< double , 3 >::extent_gen()[ 0 ][ 0 ][ 0 ] );
  /* TODO commented away until HeatBlock are properly managed
  v_pollutant_heat_rho.resize(
            boost::multi_array< double , 3 >::extent_gen()[ 0 ][ 0 ][ 0 ] );
  */
 }

 /* TODO commented away until HeatBlock are properly managed
 ::deserialize( group , "PowerHeatRho" , f_number_units ,
                v_power_heat_rho , true , true );

 if( f_number_heat_blocks )
  ::deserialize( group , "HeatNode" ,
                 f_number_heat_blocks , v_heat_node , true , true );
 */
 // TODO
 /* Notice that for units into a HeatBlock that also are electrical
  * units, the v_heat_node variable provides another time an
  * information that is already known, i.e., to which node they
  * belong to. Of course *the two information must agree*,
  * otherwise the input file is ill-defined and exception is thrown. */

 // reset all existing sub-Block, if any
 for( auto block : v_Block )
  delete( block );
 v_Block.clear();

 // load all UnitBlock
 deserialize_sub_blocks( group , "UnitBlock_" , f_number_units );

 // Generate UnitBlock primary spinning reserve variables
 unsigned int what = 0;
 if( f_number_primary_zones > 0 )
  what += 1;

 // Generate UnitBlock secondary spinning reserve variables
 if( f_number_secondary_zones > 0 )
  what += 2;

 // Generate UnitBlock inertia reserve variables
 if( f_number_inertia_zones > 0 )
  what += 4;

 if( what > 0 )
  for( auto * b : v_Block )
   if( auto ub = dynamic_cast< UnitBlock * >( b ) )
    ub->set_reserve_vars( what );

 // load all NetworkBlock, if any
 deserialize_network_blocks( group );

 if( v_network_blocks.empty() && ( ! v_active_power_demand.num_elements() ) )
  throw( std::invalid_argument( "UCBlock::deserialize: ActivePowerDemand "
                                "mandatory if no NetworkBlocks" ) );

 // if number_nodes == 1, NetworkBlocks are useless and therefore removed
 if( number_nodes == 1 ) {

  if( ! v_active_power_demand.num_elements() ) {

   // if active power demand is not defined, do it now and preload it with
   // zeros in case some NetworkBlock is not there
   v_active_power_demand.resize(
    boost::multi_array< double , 2 >::extent_gen()[ 1 ][ f_time_horizon ] );
   for( auto apdit = v_active_power_demand.data() ;
        apdit != v_active_power_demand.data() + f_time_horizon ; )
    *(apdit++) = 0;
  }

  if( ! v_network_blocks.empty() ) {

   Index t = 0;
   for( Index n = 0 ; n < f_number_networks ; ++n )
    if( v_network_blocks[ n ] ) {
     for( Index i = 0 ;
          i < v_network_blocks[ n ]->get_number_intervals() ;
          ++i , ++t ) {
      auto ad = v_network_blocks[ n ]->get_active_demand( i );
      if( ad )
       v_active_power_demand[ 0 ][ t ] = ad[ 0 ];
     }
    }
   v_network_blocks.clear();
   v_Block.resize( f_number_units );
  }

 } else {  // number_nodes > 1

  // if they don't exist, create them now as NetworkBlock
  if( v_network_blocks.empty() ) {
   v_network_blocks.resize( f_number_networks );
   v_Block.resize( f_number_units + f_number_networks );

   delete( f_NetworkData );
   f_NetworkData = static_cast< NetworkBlock::NetworkData * >(
    NetworkBlock::NetworkData::new_NetworkData( network_data_classname ) );
   f_NetworkData->deserialize( group );
  }

  Index t = 0;
  for( Index n = 0 ; n < f_number_networks ; ++n ) {

   auto nbi = v_network_blocks[ n ];
   if( ! nbi ) {  // NetworkBlock n does not exist: create a NetworkBlock
    nbi = dynamic_cast< NetworkBlock * >(
     new_Block( network_block_classname , this ) );
    v_network_blocks[ n ] = nbi;
    v_Block[ f_number_units + n ] = nbi;
    // since the ECNetworkBlock does not exist before, and we just created it,
    // we need to set its "NumberIntervals" and "ConstantTerm"
    nbi->set_number_intervals( v_start_network_intervals[ n + 1 ] -
                               v_start_network_intervals[ n ] );
    nbi->set_constant_term( v_network_constant_terms[ n ] );
   }

   if( ! nbi->get_NetworkData() ) {
    if( ! f_NetworkData )
     throw( std::invalid_argument( "UCBlock::deserialize: NetworkData "
                                   "missing in NetworkBlock " +
                                   std::to_string( n ) + " and in UCBlock" ) );
    nbi->set_NetworkData( f_NetworkData );
   }

   std::vector< std::vector< double > > ap_v;
   ap_v.resize( v_network_blocks[ n ]->get_number_intervals() ,
                std::vector< double >( number_nodes ) );

   for( Index i = 0 ;
        i < v_network_blocks[ n ]->get_number_intervals() ;
        ++i , ++t ) {

    if( ! nbi->get_active_demand( i ) ) {
     if( ! v_active_power_demand.num_elements() )
      throw( std::invalid_argument( "UCBlock::deserialize: ActivePowerDemand "
                                    "missing in UCBlock and in NetworkBlock "
                                    + std::to_string( n ) ) );
     typedef boost::multi_array_types::index_range range;
     auto ap_c = v_active_power_demand[
      boost::indices[ range( 0 , number_nodes ) ][ t ] ];
     std::copy( ap_c.begin() , ap_c.end() , ap_v[ i ].begin() );
    }
   }
   nbi->set_ActiveDemand( ap_v );
  }

  Index sum_intervals = std::accumulate(
   v_network_blocks.begin() , v_network_blocks.end() , 0 ,
   []( int init , const NetworkBlock * nb ) {
    return( init + nb->get_number_intervals() );
   } );

  // sum_intervals == f_time_horizon, i.e., we receive in input / we create n
  // `NetworkBlock`(s) and the sum of intervals spanned by each of
  // them is equal to the time horizon of the problem, throw exception otherwise
  if( sum_intervals != f_time_horizon )
   throw( std::invalid_argument( "UCBlock::deserialize: The sum of the number "
                                 "of intervals spanned by each NetworkBlock "
                                 "must be equal to the number of time horizon." ) );

  // v_active_power_demand used up, disband it
  v_active_power_demand.resize(
   boost::multi_array< double , 2 >::extent_gen()[ 0 ][ 0 ] );

 }  // end( else( number_nodes > 1 ) )

 /* TODO commented away until HeatBlock are properly managed
 if( f_number_heat_blocks )
  deserialize_sub_blocks( group , "HeatBlock_" , f_number_heat_blocks );
 */

 if( ! ::deserialize_dim( group , "NumberElectricalGenerators" ,
                          f_number_elc_generators ) ) {
  f_number_elc_generators = 0;
  for( Index i = 0 ; i < f_number_units ; ++i )
   f_number_elc_generators += static_cast< UnitBlock * >(
    v_Block[ i ] )->get_number_generators();
 }

 ::deserialize( group , "GeneratorNode" , f_number_elc_generators ,
                v_generator_node , true , true );

 // store the min and max node injection into each NetworkBlock

 if( ! v_network_blocks.empty() ) {

  Index t = 0;
  for( Index n = 0 ; n < f_number_networks ; ++n )

   for( Index i = 0 ;
        i < v_network_blocks[ n ]->get_number_intervals() ;
        ++i , ++t )

    for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {

     double min_node_injection = 0.0;
     double max_node_injection = 0.0;

     Index elc_generator = 0;
     for( Index unit_id = 0 ; unit_id < f_number_units ; unit_id++ ) {

      const auto unit_block = get_unit_block( unit_id );

      for( Index g = 0 ; g < unit_block->get_number_generators() ;
           ++g , ++elc_generator ) {

       if( node_id != v_generator_node[ elc_generator ] )
        continue;

       auto fixed_consumption = unit_block->get_fixed_consumption( g );
       min_node_injection += std::min( unit_block->get_min_power( t , g ) ,
                                       fixed_consumption
                                       ? -fixed_consumption[ t ] : 0.0 );
       max_node_injection += std::max( 0.0 ,
                                       unit_block->get_max_power( t , g ) );
      }
     }

     v_network_blocks[ n ]->set_min_node_injection( i , node_id ,
                                                    min_node_injection );
     v_network_blocks[ n ]->set_max_node_injection( i , node_id ,
                                                    max_node_injection );
    }
 }

 // finally call the method of the base class
 Block::deserialize( group );

}  // end( UCBlock::deserialize )

/*--------------------------------------------------------------------------*/

void UCBlock::generate_abstract_constraints( Configuration * stcc )
{
 if( constraints_generated() )  // constraints have already been generated
  return;                       // nothing to do

 // generate abstract constraints in all the sub-Block

 Block::generate_abstract_constraints( stcc );

 // generate the abstract constraints of UCBlock

 generate_node_injection_constraints();
 generate_primary_demand_constraints();
 generate_secondary_demand_constraints();
 generate_inertia_demand_constraints();
 generate_pollutant_budget_constraints();
 generate_heat_constraints();

 set_constraints_generated();

}  // end( UCBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void UCBlock::generate_node_injection_constraints( void )
{
 const auto number_nodes = get_number_nodes();

 v_node_injection_Const.resize(
  boost::multi_array< FRowConstraint , 2 >::extent_gen()
  [ f_time_horizon ][ number_nodes ] );

 if( number_nodes > 0 ) {  // well, that'd be curious, but ...

  if( number_nodes == 1 ) {
   // special case: in a bus network there are no NetworkBlocks and the node
   // injection constraints actually are active power demand constraints
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {  // for each time instant
    // initialise demand as active power
    auto rhs = v_active_power_demand[ 0 ][ t ];

    // each generator surely contributes with active power, but it may also
    // contribute with fixed consumption linked to commitment status, so
    // the number of nonzeros can be at most twice the number of generators
    LinearFunction::v_coeff_pair vc( 2 * f_number_elc_generators );
    auto vcit = vc.begin();

    for( Index i = 0 ; i < f_number_units ; ++i ) {  // for each unit
     const auto unit_block = get_unit_block( i );
     const auto scale = unit_block->get_scale();

     // for each electrical generator within the unit
     for( Index g = 0 ; g < unit_block->get_number_generators() ; ++g ) {

      // surely add the contribution of the corresponding active power
      *(vcit++) = std::pair( &unit_block->get_active_power( g )[ t ] ,
                             scale );

      // if the generator also has nonzero fixed consumption at t
      // fixed consumption happens when the generator is off, and it
      // therefore has the form fc[ t ] * ( 1 - u[ t ] ); thus, the
      // RHS of the constraint also has to be decreased by fc[ t ]. note
      // that a unit with no commitment is always on, and therefore the
      // fixed consumption is always 0
      if( auto fc = unit_block->get_fixed_consumption( g ) )
       if( fc[ t ] )
        if( auto u = unit_block->get_commitment( g ) ) {
         const auto fixed_consumption = fc[ t ] * scale;
         // add the contribution of the corresponding commitment variables
         *(vcit++) = std::pair( &u[ t ] , -fixed_consumption );
         rhs -= fixed_consumption;    // update the RHS
        }
     }  // end( for( g ) )
    }  // end( for( i ) )

    // set the final RHS of the constraint (equality constraint)
    v_node_injection_Const[ t ][ 0 ].set_both( rhs , eNoMod );
    // resize vc so that it's of the right length
    vc.resize( std::distance( vc.begin() , vcit ) );
    // construct and pass the LinearFunction to the FRowConstraint
    v_node_injection_Const[ t ][ 0 ].set_function(
     new LinearFunction( std::move( vc ) ) , eNoMod );
   }  // end( for( t ) )

  } else {  // number_nodes > 1

   // Network needs GeneratorNode

   Index t = 0;
   for( Index n = 0 ; n < f_number_networks ; ++n ) {

    for( Index i = 0 ;
         i < v_network_blocks[ n ]->get_number_intervals() ;
         ++i , ++t ) {

     auto node_injection = v_network_blocks[ n ]->get_node_injection( i );

     for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {

      auto lf = new LinearFunction();

      lf->add_variable( &node_injection[ node_id ] , -1.0 , eNoMod );

      double rhs = 0.0;

      Index elc_generator = 0;
      for( Index unit_id = 0 ; unit_id < f_number_units ; unit_id++ ) {

       const auto unit_block = get_unit_block( unit_id );
       const auto scale = unit_block->get_scale();

       for( Index generator = 0 ;
            generator < unit_block->get_number_generators() ;
            ++generator , ++elc_generator ) {

        if( node_id != v_generator_node[ elc_generator ] )
         continue;

        if( auto ap = unit_block->get_active_power( generator ) ) {
         auto active_power = &ap[ t ];
         lf->add_variable( active_power , scale , eNoMod );
        }

        if( auto fc = unit_block->get_fixed_consumption( generator ) ) {
         if( auto c = unit_block->get_commitment( generator ) ) {
          auto fixed_consumption = fc[ t ] * scale;
          auto commitment = &c[ t ];
          lf->add_variable( commitment , -fixed_consumption , eNoMod );
          rhs -= fixed_consumption;
         }
        }
       }
      }
      v_node_injection_Const[ t ][ node_id ].set_both( rhs , eNoMod );
      v_node_injection_Const[ t ][ node_id ].set_function( lf );
     }
    }
   }
  }

  add_static_constraint( v_node_injection_Const , "node_injection_c" );
 }

}  // end( UCBlock::generate_node_injection_constraints )

/*--------------------------------------------------------------------------*/

void UCBlock::generate_primary_demand_constraints( void )
{
 if( f_number_primary_zones == 0 )
  return;

 v_PrimaryDemand_Const.resize(
  boost::multi_array< FRowConstraint , 2 >::extent_gen()
  [ f_time_horizon ][ f_number_primary_zones ] );

 // We assume that, if a generator has primary spinning reserve for a time
 // instant, then it has primary spinning reserve for all time instants.
 primary_var_index.resize
  ( boost::multi_array< Range , 2 >::
    extent_gen()[ f_number_units ][ f_number_primary_zones ] );

 std::fill( primary_var_index.data() , primary_var_index.data() +
                                       primary_var_index.num_elements() ,
            std::pair{ Inf< Index >() , Inf< Index >() } );

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

       if( primary_var_index[ unit_id ][ zone_id ].first == Inf< Index >() ) {
        // This is the first Variable of this unit to be added to the
        // LinearFunction, so we store its index, which is given by the
        // current number of active Variables of the LinearFunction (right
        // before this Variable is added).

        // Since all time steps have the same structure, this must be the
        // first time step.
        assert( t == 0 );

        const auto num_active_var = lf->get_num_active_var();
        primary_var_index[ unit_id ][ zone_id ].first = num_active_var;
        primary_var_index[ unit_id ][ zone_id ].second = num_active_var;
       }

       if( t == 0 )
        // Increment the upper bound of the range.
        ++primary_var_index[ unit_id ][ zone_id ].second;

       // Now we add the primary reserve variable to the LinearFunction.
       auto primary_spinning_reserve = &primary_s_r[ t ];
       lf->add_variable( primary_spinning_reserve , scale );
      }

     }  // end( for( generator ) )
    }  // end( for( unit_id ) )
   }  // end( for( node_id ) )

   const auto demand = get_primary_demand()[ zone_id ][ t ];
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

 // We assume that, if a generator has secondary spinning reserve for a time
 // instant, then it has secondary spinning reserve for all time instants.
 secondary_var_index.resize
  ( boost::multi_array< Range , 2 >::
    extent_gen()[ f_number_units ][ f_number_secondary_zones ] );

 std::fill( secondary_var_index.data() , secondary_var_index.data() +
                                         secondary_var_index.num_elements() ,
            std::pair{ Inf< Index >() , Inf< Index >() } );

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

       if( secondary_var_index[ unit_id ][ zone_id ].first == Inf< Index >() ) {
        // This is the first Variable of this unit to be added to the
        // LinearFunction, so we store its index, which is given by the
        // current number of active Variables of the LinearFunction (right
        // before this Variable is added)

        // Since all time steps have the same structure, this must be the
        // first time step.
        assert( t == 0 );

        const auto num_active_var = lf->get_num_active_var();
        secondary_var_index[ unit_id ][ zone_id ].first = num_active_var;
        secondary_var_index[ unit_id ][ zone_id ].second = num_active_var;
       }

       if( t == 0 )
        // Increment the upper bound of the range.
        ++secondary_var_index[ unit_id ][ zone_id ].second;

       // Now we add the secondary reserve variable to the LinearFunction.
       auto secondary_spinning_reserve = &secondary_s_r[ t ];
       lf->add_variable( secondary_spinning_reserve , scale );
      }

     }  // end( for( generator ) )
    }  // end( for( unit_id ) )
   }  // end( for( node_id ) )

   const auto demand = get_secondary_demand()[ zone_id ][ t ];
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

 // We assume that, if a generator has commitment variable, inertia
 // commitment, inertia power, or active power variable for some time instant,
 // then it has the same thing for all time instants.
 inertia_var_index.resize
  ( boost::multi_array< Index , 2 >::
    extent_gen()[ f_number_units ][ f_number_inertia_zones ] );

 std::fill( inertia_var_index.data() , inertia_var_index.data() +
                                       inertia_var_index.num_elements() ,
            Inf< Index >() );

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
      auto inertia_commitment = unit_block->get_inertia_commitment( generator );

      if( commitment && inertia_commitment ) {

       // The term with the commitment variable will be added to the function.

       if( inertia_var_index[ unit_id ][ zone_id ] == Inf< Index >() ) {
        // This is the first Variable of this unit to be added to the
        // LinearFunction, so we store its index, which is given by the
        // current number of active Variables of the LinearFunction (right
        // before this Variable is added).

        // Since all time steps have the same structure, this must be the
        // first time step.
        assert( t == 0 );

        const auto num_active_var = lf->get_num_active_var();
        inertia_var_index[ unit_id ][ zone_id ] = num_active_var;
       }

       auto commitment_t = &commitment[ t ];
       auto coefficient = scale * inertia_commitment[ t ];
       lf->add_variable( commitment_t , coefficient );
      }

      auto active_power = unit_block->get_active_power( generator );
      auto inertia_power = unit_block->get_inertia_power( generator );

      if( active_power && inertia_power ) {
       // The term with the active power will be added to the function.

       if( inertia_var_index[ unit_id ][ zone_id ] == Inf< Index >() ) {
        // This is the first Variable of this unit to be added to the
        // LinearFunction, so we store its index, which is given by the
        // current number of active Variables of the LinearFunction (right
        // before this Variable is added).

        // Since all time steps have the same structure, this must be the
        // first time step.
        assert( t == 0 );

        const auto num_active_var = lf->get_num_active_var();
        inertia_var_index[ unit_id ][ zone_id ] = num_active_var;
       }

       auto active_power_t = &active_power[ t ];
       auto coefficient = scale * inertia_power[ t ];
       lf->add_variable( active_power_t , coefficient );
      }

     }  // end( for( generator ) )
    }  // end( for( unit_id ) )
   }  // end( for( node_id ) )

   const auto demand = get_inertia_demand()[ zone_id ][ t ];
   v_InertiaDemand_Const[ t ][ zone_id ].set_lhs( demand );
   v_InertiaDemand_Const[ t ][ zone_id ].set_rhs( Inf< double >() );
   v_InertiaDemand_Const[ t ][ zone_id ].set_function( lf );

  }  // end( for( zone_id ) )
 }  // end( for( t ) )

 add_static_constraint( v_InertiaDemand_Const , "inertia_demand_c" );

}  // end( UCBlock::generate_inertia_demand_constraints )

/*--------------------------------------------------------------------------*/

void UCBlock::generate_pollutant_budget_constraints( void )
{
 // TODO These constraints must be fixed

 const auto number_nodes = get_number_nodes();

 if( f_number_pollutants > 0 ) {

  v_PollutantBudget_Const.resize(
   v_number_pollutant_zones[ f_total_number_pollutant_zones ] );

  LinearFunction::v_coeff_pair vars;

  if( number_nodes == 1 ) {

   for( Index pollutant = 0 ; pollutant < f_number_pollutants ; ++pollutant ) {

    for( Index zone = 0 ; zone < v_number_pollutant_zones[ pollutant ] ;
         ++zone ) {

     for( Index t = 0 ; t < f_time_horizon ; ++t ) {

      for( Index unit_id = 0 ; unit_id < f_number_units ; ++unit_id ) {

       const auto unit_block = get_unit_block( unit_id );
       const auto scale = unit_block->get_scale();

       for( Index generator = 0 ;
            generator < unit_block->get_number_generators() ; ++generator ) {

        auto node_id = get_generator_node()[ generator ];
        auto zone_id = get_pollutant_zone()[ pollutant ][ node_id ];

        if( zone_id >= v_number_pollutant_zones[ pollutant ] )
         continue;  // this unit does not belong to any zone

        if( auto ap = unit_block->get_active_power( generator ) ) {
         auto active_power = &ap[ t ];
         auto rho = get_pollutant_rho()[ t ][ pollutant ][ generator ];
         auto coefficient = scale * rho;
         vars.push_back( std::make_pair( active_power , coefficient ) );
        }
       }
      }

      // Terms associated with heat-only generation units
      /* TODO commented away until HeatBlock are properly managed
      if( f_number_heat_blocks > 0 ) {  // TODO Do we have any HeatBlock?

       for( Index h = 0; h < f_number_heat_blocks; ++h ) {

        for( std::vector< Index >::size_type i = 0;
             i <= f_number_units; ++i ) {

         //auto unit_id = v_heat_only_units[ i ];
         auto heat_id = v_heat_set[ i ];

         auto zone_id = get_pollutant_zone()[pollutant][ h ];
         if( zone_id >= v_number_pollutant_zones[ pollutant ] )
          continue;  // this unit does not belong to any zone

         auto heat = get_heat_block()[ h ]->get_heat()[ t ][ i ];
         auto rho = get_pollutant_heat_rho()[ t ][pollutant][ h ];

         auto lf = dynamic_cast< LinearFunction * >
         ( v_PollutantBudget_Const[pollutant][zone_id].get_function() );

         lf->add_variable( &heat, rho );
        }
       }
      }
      */

      v_PollutantBudget_Const[ pollutant ][ zone ].set_rhs(
       v_pollutant_budget[ v_number_pollutant_zones[ pollutant ] ][ pollutant ] );
      v_PollutantBudget_Const[ pollutant ][ zone ].set_lhs( -Inf< double >() );
      v_PollutantBudget_Const[ pollutant ][ zone ].set_function(
       new LinearFunction( std::move( vars ) ) );
     }
    }
   }
  } else {

   for( Index pollutant = 0 ; pollutant < f_number_pollutants ; ++pollutant ) {

    for( Index zone = 0 ; zone < v_number_pollutant_zones[ pollutant ] ;
         ++zone ) {

     for( Index t = 0 ; t < f_time_horizon ; ++t ) {

      Index pollutant_zone = 0;
      for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {

       if( zone == v_pollutant_zones[ pollutant ][ node_id ] ) {

        Index generator_id = 0;
        for( Index elc_generator = 0 ; elc_generator < f_number_elc_generators ;
             ++elc_generator ) {

         if( node_id == v_generator_node[ elc_generator ] ) {

          auto block = get_nested_Blocks()[ generator_id ];
          auto unit_block = dynamic_cast< UnitBlock * >(block);
          if( ! unit_block )
           continue;

          const auto scale = unit_block->get_scale();

          for( Index generator = 0 ;
               generator < unit_block->get_number_generators() ; ++generator ) {

           auto node_id = get_generator_node()[ generator ];
           auto zone_id = get_pollutant_zone()[ pollutant ][ node_id ];

           if( zone_id >= v_number_pollutant_zones[ pollutant ] )
            continue;  // this unit does not belong to any zone

           if( auto ap = unit_block->get_active_power( generator ) ) {
            auto active_power = &ap[ t ];
            auto rho = get_pollutant_rho()[ t ][ pollutant ][ generator ];
            auto coefficient = scale * rho;
            vars.push_back( std::make_pair( active_power , coefficient ) );
           }
          }
         }
         generator_id++;
        }
       }
       pollutant_zone++;
      }

      // Terms associated with heat-only generation units
      /* TODO commented away until HeatBlock are properly managed
      if( f_number_heat_blocks > 0 ) {  // TODO Do we have any HeatBlock?

       for( Index h = 0; h < f_number_heat_blocks; ++h ) {

        for( std::vector< Index >::size_type i = 0;
             i <= f_number_units; ++i ) {

         //auto unit_id = v_heat_only_units[ i ];
         auto heat_id = v_heat_set[ i ];

         auto zone_id = get_pollutant_zone()[pollutant][h];
         if( zone_id >= v_number_pollutant_zones[pollutant] )
          continue;  // this unit does not belong to any zone

         auto heat = get_heat_block()[h]->get_heat()[ t ][ i ];
         auto rho = get_pollutant_heat_rho()[ t ][ pollutant ][ h ];

         auto lf = dynamic_cast< LinearFunction * >
         ( v_PollutantBudget_Const[pollutant][zone_id].get_function() );

         lf->add_variable( &heat, rho );
        }
       }
      }
      */

      v_PollutantBudget_Const[ pollutant ][ zone ].set_rhs(
       v_pollutant_budget[ v_number_pollutant_zones[ pollutant ] ][ pollutant ] );
      v_PollutantBudget_Const[ pollutant ][ zone ].set_lhs( -Inf< double >() );
      v_PollutantBudget_Const[ pollutant ][ zone ].set_function(
       new LinearFunction( std::move( vars ) ) );
     }
    }
   }
  }

  add_static_constraint(
   v_PollutantBudget_Const[ f_total_number_pollutant_zones ] );
 }

}  // end( UCBlock::generate_pollutant_budget_constraints )

/*--------------------------------------------------------------------------*/

void UCBlock::generate_heat_constraints( void )
{
 /* TODO commented away until HeatBlock are properly managed
 // Heat constraints.
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 // TODO Deal with the scaling of UnitBlock

 if( f_number_heat_blocks > 0 ) {

  Index num_constraints_per_time = 0;

  // List of units that produce electricity and belong to some
  // HeatBlock.
  std::vector< Index > electricity_generators_inside_a_heat_block;

  if( v_power_Heat_Rho_Const.size() != f_time_horizon ) {

   // this should only happen once
   assert( v_power_Heat_Rho_Const.empty() );

   auto is_electricity_producing_inside_a_heat_block =
    [ this ]( Index unit ) {
     for( Index heat_block_id = 0; heat_block_id < f_number_heat_blocks;
          ++heat_block_id ) {
      if( get_heat_set()[ unit ] <
          v_heat_blocks[ heat_block_id ]->get_number_heat_generators() )
       return( true );
     }
     return( false );
    };

   for( Index unit_id = 0; unit_id < f_number_units; ++unit_id ) {
    if( is_electricity_producing_inside_a_heat_block( unit_id ) ) {
     electricity_generators_inside_a_heat_block
     [ num_constraints_per_time++ ] = unit_id;
    }
   }

   v_power_Heat_Rho_Const.resize(
    boost::multi_array< FRowConstraint * , 2 >::extent_gen()
    [ f_time_horizon ][ num_constraints_per_time ] );
  }

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {

   for( std::vector< Index >::size_type constraint_id = 0;
        constraint_id < num_constraints_per_time; ++constraint_id ) {

    auto generator_id = electricity_generators_inside_a_heat_block[ constraint_id ];

    for( Index heat_block_id = 0; heat_block_id < f_number_heat_blocks;
         ++heat_block_id ) {

     auto heat_unit_id = get_heat_set()[ generator_id ];

     if( heat_unit_id >=
         v_heat_blocks[ heat_block_id ]->get_number_heat_generators() )
      continue;  // unit_id does not belong to heat_block_id

     if( ! v_power_Heat_Rho_Const[ t ][ constraint_id ].get_function() ) {
      v_power_Heat_Rho_Const[ t ][ constraint_id ].
       set_lhs( -Inf< double >() );
      v_power_Heat_Rho_Const[ t ][ constraint_id ].set_rhs( 0.0 );
      v_power_Heat_Rho_Const[ t ][ constraint_id ].set_function(
       new LinearFunction() );
     }

     auto active_power = get_unit_block( generator_id )->get_active_power();
     auto heat = get_heat_block() [ heat_unit_id ]->get_heat()[ t ][generator_id];
     auto power_heat_rho = get_power_heat_rho()[ generator_id ];

     auto lf = dynamic_cast< LinearFunction * >
     ( v_power_Heat_Rho_Const[ t ][ constraint_id ].get_function() );
     lf->add_variable( &heat, 1, 0 );
     lf->add_variable( &active_power[ t ][ generator_id ], -power_heat_rho );
    }
   }
  }

  add_static_constraint( v_power_Heat_Rho_Const );
 }
*/
}  // end( UCBlock::generate_heat_constraints )

/*--------------------------------------------------------------------------*/
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
 } else
  NumberNodes = group.addDim( "NumberNodes" , 1 );

 auto TimeHorizon = group.addDim( "TimeHorizon" , f_time_horizon );
 auto NumberUnits = group.addDim( "NumberUnits" , f_number_units );
 auto NumberNetworks = group.addDim( "NumberNetworks" , f_number_networks );
 auto NumberElectricalGenerators = group.addDim( "NumberElectricalGenerators" ,
                                                 f_number_elc_generators );

 auto TotalNumberPollutantZones = group.addDim( "TotalNumberPollutantZones" ,
                                                f_total_number_pollutant_zones );

 /* TODO commented away until HeatBlock are properly managed
 auto NumberHeatBlocks = group.addDim( "NumberHeatBlocks" ,
                                       f_number_heat_blocks );

 auto NumberHeatGenerators = group.addDim( "NumberHeatGenerators" ,
                                           f_number_heat_generators );
 */

 auto NumberPrimaryZones = group.addDim( "NumberPrimaryZones" ,
                                         f_number_primary_zones );
 auto NumberSecondaryZones = group.addDim( "NumberSecondaryZones" ,
                                           f_number_secondary_zones );
 auto NumberInertiaZones = group.addDim( "NumberInertiaZones" ,
                                         f_number_inertia_zones );
 auto NumberPollutants = group.addDim( "NumberPollutants" ,
                                       f_number_pollutants );

 ::serialize( group , "ActivePowerDemand" , netCDF::NcDouble() ,
              { NumberNodes , TimeHorizon } , v_active_power_demand );

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

 ::serialize( group , "NumberPollutantZones" , netCDF::NcUint() ,
              NumberPollutants , v_number_pollutant_zones );

 ::serialize( group , "PollutantZones" , netCDF::NcUint() ,
              { NumberPollutants , NumberNodes } , v_pollutant_zones );

 /* TODO commented away until this is properly managed
 ::serialize( group , "PollutantBudget" , netCDF::NcDouble() ,
              TotalNumberPollutantZones , v_pollutant_budget );
 */

 ::serialize( group , "PollutantRho" , netCDF::NcDouble() ,
              { TimeHorizon , NumberPollutants , NumberElectricalGenerators } ,
              v_pollutant_rho );

 std::vector< Index > number_intervals( v_start_network_intervals.size() , 0 );
 std::adjacent_difference( v_start_network_intervals.begin() ,
                           v_start_network_intervals.end() ,
                           number_intervals.begin() );
 // the first will be zero, so we start checking from the second
 if( std::any_of( std::next( number_intervals.begin() ) ,
                  number_intervals.end() ,
                  []( double cst ) { return( cst != 1 ); } ) )
  ::serialize( group , "StartNetworkIntervals" , netCDF::NcUint() ,
               NumberNetworks , v_start_network_intervals );

 if( std::any_of( v_network_constant_terms.begin() ,
                  v_network_constant_terms.end() ,
                  []( double cst ) { return( cst != 0 ); } ) )
  ::serialize( group , "NetworkConstantTerms" , netCDF::NcDouble() ,
               NumberNetworks , v_network_constant_terms );

 if( network_block_classname != "DCNetworkBlock" )
  ::serialize( group , "NetworkBlockClassname" , netCDF::NcString() ,
               network_block_classname );

 if( network_data_classname != "DCNetworkData" )
  ::serialize( group , "NetworkDataClassname" , netCDF::NcString() ,
               network_data_classname );

 /* TODO commented away until HeatBlock are properly managed
 ::serialize( group, "HeatSet", netCDF::NcUint(),
              NumberHeatGenerators, v_heat_set );

 ::serialize( group, "PollutantHeatRho", netCDF::NcDouble(),
              {TimeHorizon, NumberPollutants,
               NumberHeatBlocks},
              v_pollutant_heat_rho );

 ::serialize( group, "PowerHeatRho", netCDF::NcDouble(),
              NumberUnits, v_power_heat_rho );

 ::serialize( group, "HeatNode", netCDF::NcUint(),
              NumberHeatBlocks, v_heat_node );
 */

 ::serialize( group , "GeneratorNode" , netCDF::NcUint() ,
              NumberElectricalGenerators , v_generator_node );

 // Serialize sub-blocks

 for( Index i = 0 ; i < f_number_units ; ++i ) {
  auto sub_block = get_unit_block( i );
  auto sub_group = group.addGroup( "UnitBlock_" + std::to_string( i ) );
  sub_block->serialize( sub_group );
 }

 for( Index t = 0 ; t < f_time_horizon ; ++t )
  if( auto sub_block = get_network_block( t ) ) {
   auto sub_group = group.addGroup( "NetworkBlock_" + std::to_string( t ) );
   sub_block->serialize( sub_group );
  }

 /* TODO commented away until HeatBlock are properly managed
 for( Index i = 0; i < f_number_heat_blocks; ++i ) {
  auto sub_block = get_heat_block( i );
  auto sub_group = group.addGroup( "HeatBlock_" + std::to_string( i ) );
  sub_block->serialize( sub_group );
  }
 */

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
 std::vector< Index > modified_units;

 // TODO Handle GroupModification in order to deal with multiple UnitBlockMod
 // at the same time.

 if( const auto tmod = dynamic_cast< UnitBlockMod * >( mod.get() ) ) {
  if( tmod->type() == UnitBlockMod::eScale ) {
   auto unit_id = inspection::get_block_index( tmod->get_Block() );
   modified_units.push_back( unit_id );
  }
 }

 if( ! modified_units.empty() ) {
  // Sort the IDs of the modified units
  std::sort( modified_units.begin() , modified_units.end() );

  update_node_injection_constraints( modified_units );
  update_primary_demand_constraints( modified_units );
  update_secondary_demand_constraints( modified_units );
  update_inertia_demand_constraints( modified_units );

  // TODO Implement the following methods when their constraints have been
  // properly implemented.

  // update_pollutant_budget_constraints( modified_units );
  // update_heat_constraints( modified_units );
 }

 Block::add_Modification( mod , chnl );
}

/*--------------------------------------------------------------------------*/

void UCBlock::update_node_injection_constraints(
 const std::vector< Index > & modified_units )
{
 if( ( ! constraints_generated() ) ||
     ( v_node_injection_Const.empty() ) || modified_units.empty() )
  return;

 // Lambda for determining if some unit has been modified

 auto has_been_modified = [ &modified_units ]( Index unit_id ) {
  if( std::binary_search( modified_units.begin() ,
                          modified_units.end() , unit_id ) )
   return( true );
  return( false );
 };

 // Compute the total number of generators of the units that have been
 // modified

 auto total_num_generators = 0;
 for( auto unit_id : modified_units ) {
  const auto unit_block = get_unit_block( unit_id );
  total_num_generators += unit_block->get_number_generators();
 }

 const auto number_nodes = get_number_nodes();

 if( number_nodes > 0 ) {

  if( number_nodes == 1 ) {
   for( Index t = 0 ; t < f_time_horizon ; ++t ) {  // for each time instant

    auto & constraint = v_node_injection_Const[ t ][ 0 ];

    // This will store the coefficients that must be updated, i.e., those of
    // the active Variables that belong to the units that have been modified.
    LinearFunction::Vec_FunctionValue coefficients;
    coefficients.reserve( 2 * total_num_generators );

    // Subset that will store the indices of the active Variables whose
    // coefficients have changed.
    Subset subset;
    subset.reserve( 2 * total_num_generators );

    // Index of the current active Variable
    Index active_var_index = 0;

    // Initialise demand as active power
    auto rhs = v_active_power_demand[ 0 ][ t ];

    for( Index i = 0 ; i < f_number_units ; ++i ) {  // for each unit

     const auto unit_block = get_unit_block( i );
     const auto scale = unit_block->get_scale();
     const bool modified = has_been_modified( i );

     // for each electrical generator within the unit
     for( Index g = 0 ; g < unit_block->get_number_generators() ; ++g ) {

      if( modified ) {
       // update the coefficient of the active power variable
       coefficients.push_back( scale );
       subset.push_back( active_var_index );

       assert( active_var_index < constraint.get_num_active_var() );
       assert( constraint.get_active_var( active_var_index )->get_Block()
               == unit_block );
      }

      // increment due to the active power variable
      ++active_var_index;

      if( auto fc = unit_block->get_fixed_consumption( g ) )
       if( fc[ t ] )
        if( unit_block->get_commitment( g ) ) {
         const auto fixed_consumption = fc[ t ] * scale;
         rhs -= fixed_consumption;    // update the RHS

         if( modified ) {
          // update the coefficient of the commitment variable
          coefficients.push_back( scale );
          subset.push_back( active_var_index );

          assert( active_var_index < constraint.get_num_active_var() );
          assert( constraint.get_active_var( active_var_index )->get_Block()
                  == unit_block );
         }

         // increment due to the commitment variable
         ++active_var_index;
        }
     }  // end( for( g ) )
    }  // end( for( i ) )

    // Finally, we update the RHS of the constraint and the coefficients of
    // the active Variables that have been modified. Notice that the
    // (abstract) Modifications that will be issued as a result of this update
    // do not concern this UCBlock.

    // update the RHS of the constraint (equality constraint)
    constraint.set_both( rhs , eNoBlck );

    // update the coefficients
    static_cast< LinearFunction * >
    ( constraint.get_function() )->modify_coefficients
     ( std::move( coefficients ) , std::move( subset ) , true , eNoBlck );

   }  // end( for( t ) )

  } else {  // number_nodes > 1

   // Network needs GeneratorNode

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {

    for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {

     // This will store the coefficients that must be updated, i.e., those of
     // the active Variables that belong to the units that have been modified.
     LinearFunction::Vec_FunctionValue coefficients;
     coefficients.reserve( 2 * total_num_generators );

     // Subset that will store the indices of the active Variables whose
     // coefficients have changed.
     Subset subset;
     subset.reserve( 2 * total_num_generators );

     // Index of the current active Variable
     Index active_var_index = 0;

     auto & constraint = v_node_injection_Const[ t ][ node_id ];

     // increment due to the node injection variable
     ++active_var_index;

     double rhs = 0.0;

     Index elc_generator = 0;
     for( Index unit_id = 0 ; unit_id < f_number_units ; unit_id++ ) {

      const auto unit_block = get_unit_block( unit_id );
      const auto scale = unit_block->get_scale();
      const bool modified = has_been_modified( unit_id );

      for( Index generator = 0 ;
           generator < unit_block->get_number_generators() ;
           ++generator , ++elc_generator ) {

       if( node_id != v_generator_node[ elc_generator ] )
        continue;

       if( unit_block->get_active_power( generator ) ) {

        if( modified ) {
         // update the coefficient of the active power variable
         coefficients.push_back( scale );
         subset.push_back( active_var_index );

         assert( active_var_index < constraint.get_num_active_var() );
         assert( constraint.get_active_var( active_var_index )->get_Block()
                 == unit_block );
        }

        // increment due to the active power variable
        ++active_var_index;
       }

       if( auto fc = unit_block->get_fixed_consumption( generator ) ) {
        if( unit_block->get_commitment( generator ) ) {
         auto fixed_consumption = fc[ t ] * scale;

         if( modified ) {
          // update the coefficient of the commitment variable
          coefficients.push_back( - fixed_consumption );
          subset.push_back( active_var_index );

          assert( active_var_index < constraint.get_num_active_var() );
          assert( constraint.get_active_var( active_var_index )->get_Block()
                  == unit_block );
         }

         // increment due to the commitment variable
         ++active_var_index;

         rhs -= fixed_consumption;
        }
       }
      }
     }

     // Finally, we update the RHS of the constraint and the coefficients of
     // the active Variables that have been modified. Notice that the
     // (abstract) Modifications that will be issued as a result of this
     // update do not concern this UCBlock.

     // update the RHS of the constraint (equality constraint)
     constraint.set_both( rhs , eNoBlck );

     // update the coefficients
     static_cast< LinearFunction * >( constraint.get_function() )->
      modify_coefficients( std::move( coefficients ) , std::move( subset ) ,
                           true , eNoBlck );

    }  // end( for( node_id ) )
   }  // end( for( t ) )
  }
 }  // end( if( number_nodes > 0 ) )

}  // end( UCBlock::update_node_injection_constraints )

/*--------------------------------------------------------------------------*/

void UCBlock::update_primary_demand_constraints(
 const std::vector< Index > & modified_units )
{
 if( ( ! constraints_generated() ) || ( v_PrimaryDemand_Const.empty() ) ||
     modified_units.empty() )
  return;  // there is nothing to be updated

 // Indices of the zones that are affected by the modified units.
 std::set< Index > affected_zones;

 // Number of modified generators in each zone.
 std::vector< Index > num_generators_per_zone( f_number_primary_zones , 0 );

 // Collect the affected zones and count the number of affected generators in
 // each zone.

 Index elc_generator = 0;
 Index overall_unit_id = 0;
 for( const auto unit_id : modified_units ) {

  // Skip the units that have not been modified.
  while( overall_unit_id < unit_id ) {
   elc_generator += get_unit_block( overall_unit_id )->get_number_generators();
   ++overall_unit_id;
  }

  const auto unit_block = get_unit_block( unit_id );
  const auto num_generators = unit_block->get_number_generators();

  for( Index g = 0 ; g < num_generators ; ++g , ++elc_generator ) {
   const auto zone = get_primary_zone( elc_generator );
   affected_zones.insert( zone );
   ++num_generators_per_zone[ zone ];
  }

  ++overall_unit_id;
 }

 // Now loop over all affected constraints

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {
  for( const auto zone_id : affected_zones ) {

   auto & constraint = v_PrimaryDemand_Const[ t ][ zone_id ];

   // This will store the coefficients that must be updated, i.e., those of
   // the active Variables that belong to the units that have been modified.
   LinearFunction::Vec_FunctionValue coefficients;
   coefficients.reserve( num_generators_per_zone[ zone_id ] );

   // Subset that will store the indices of the active Variables whose
   // coefficients have changed.
   Subset subset;
   subset.reserve( num_generators_per_zone[ zone_id ] );

   for( const auto unit_id : modified_units ) {

    if( primary_var_index[ unit_id ][ zone_id ].first == Inf< Index >() ) {
     // This unit has no active Variable in the primary demand constraints
     // associated with zone "zone_id".
     continue;
    }

    const auto unit_block = get_unit_block( unit_id );
    const auto scale = unit_block->get_scale();

    // Indices of the active Variables of the current UnitBlock: the indices
    // are consecutive and are given by the open-closed interval
    // [ primary_var_index[ unit_id ][ zone_id ].first ,
    //   primary_var_index[ unit_id ][ zone_id ].second ).
    const auto num_variables = primary_var_index[ unit_id ][ zone_id ].second -
                               primary_var_index[ unit_id ][ zone_id ].first;
    std::vector< Index > var_indices( num_variables );
    std::iota( var_indices.begin() , var_indices.end() ,
               primary_var_index[ unit_id ][ zone_id ].first );

    for( const auto var_index : var_indices ) {
     assert( var_index < constraint.get_num_active_var() );
     assert( constraint.get_active_var( var_index )->get_Block()
             == unit_block );
    }

    subset.insert( subset.end() , var_indices.begin() , var_indices.end() );
    coefficients.insert( coefficients.end() , num_variables , scale );

   }  // end( for( modified_units ) )

   // Update the coefficients of the active variables
   static_cast< LinearFunction * >( constraint.get_function() )->
    modify_coefficients( std::move( coefficients ) , std::move( subset ) ,
                         false , eNoBlck );

  }  // end( for( zone_id ) )
 }  // end( for( t ) )

}  // end( UCBlock::update_primary_demand_constraints )

/*--------------------------------------------------------------------------*/

void UCBlock::update_secondary_demand_constraints(
 const std::vector< Index > & modified_units )
{
 if( ( ! constraints_generated() ) || ( v_SecondaryDemand_Const.empty() ) ||
     modified_units.empty() )
  return;  // there is nothing to be updated

 // Indices of the zones that are affected by the modified units.
 std::set< Index > affected_zones;

 // Number of modified generators in each zone.
 std::vector< Index > num_generators_per_zone( f_number_secondary_zones , 0 );

 // Collect the affected zones and count the number of affected generators in
 // each zone.

 Index elc_generator = 0;
 Index overall_unit_id = 0;
 for( const auto unit_id : modified_units ) {

  // Skip the units that have not been modified.
  while( overall_unit_id < unit_id ) {
   elc_generator += get_unit_block( overall_unit_id )->get_number_generators();
   ++overall_unit_id;
  }

  const auto unit_block = get_unit_block( unit_id );
  const auto num_generators = unit_block->get_number_generators();

  for( Index g = 0 ; g < num_generators ; ++g , ++elc_generator ) {
   const auto zone = get_secondary_zone( elc_generator );
   affected_zones.insert( zone );
   ++num_generators_per_zone[ zone ];
  }

  ++overall_unit_id;
 }

 // Now loop over all affected constraints

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {
  for( const auto zone_id : affected_zones ) {

   auto & constraint = v_SecondaryDemand_Const[ t ][ zone_id ];

   // This will store the coefficients that must be updated, i.e., those of
   // the active Variables that belong to the units that have been modified.
   LinearFunction::Vec_FunctionValue coefficients;
   coefficients.reserve( num_generators_per_zone[ zone_id ] );

   // Subset that will store the indices of the active Variables whose
   // coefficients have changed.
   Subset subset;
   subset.reserve( num_generators_per_zone[ zone_id ] );

   for( const auto unit_id : modified_units ) {

    if( secondary_var_index[ unit_id ][ zone_id ].first == Inf< Index >() ) {
     // This unit has no active Variable in the secondary demand constraints
     // associated with zone "zone_id".
     continue;
    }

    const auto unit_block = get_unit_block( unit_id );
    const auto scale = unit_block->get_scale();

    // Indices of the active Variables of the current UnitBlock: the indices
    // are consecutive and are given by the open-closed interval
    // [ secondary_var_index[ unit_id ][ zone_id ].first ,
    //   secondary_var_index[ unit_id ][ zone_id ].second ).
    const auto num_variables =
     secondary_var_index[ unit_id ][ zone_id ].second -
     secondary_var_index[ unit_id ][ zone_id ].first;
    std::vector< Index > var_indices( num_variables );
    std::iota( var_indices.begin() , var_indices.end() ,
               secondary_var_index[ unit_id ][ zone_id ].first );

    for( const auto var_index : var_indices ) {
     assert( var_index < constraint.get_num_active_var() );
     assert( constraint.get_active_var( var_index )->get_Block()
             == unit_block );
    }

    subset.insert( subset.end() , var_indices.begin() , var_indices.end() );
    coefficients.insert( coefficients.end() , num_variables , scale );

   }  // end( for( modified_units ) )

   // Update the coefficients of the active variables
   static_cast< LinearFunction * >( constraint.get_function() )->
    modify_coefficients( std::move( coefficients ) , std::move( subset ) ,
                         false , eNoBlck );

  }  // end( for( zone_id ) )
 }  // end( for( t ) )

}  // end( UCBlock::update_secondary_demand_constraints )

/*--------------------------------------------------------------------------*/

void UCBlock::update_inertia_demand_constraints(
 const std::vector< Index > & modified_units )
{
 if( ( ! constraints_generated() ) || ( v_InertiaDemand_Const.empty() ) ||
     modified_units.empty() )
  return;  // there is nothing to be updated

 // Indices of the zones that are affected by the modified units.
 std::set< Index > affected_zones;

 // Number of modified generators in each zone.
 std::vector< Index > num_generators_per_zone( f_number_inertia_zones , 0 );

 // Collect the affected zones and count the number of affected generators in
 // each zone.

 Index elc_generator = 0;
 Index overall_unit_id = 0;
 for( const auto unit_id : modified_units ) {

  // Skip the units that have not been modified.
  while( overall_unit_id < unit_id ) {
   elc_generator += get_unit_block( overall_unit_id )->get_number_generators();
   ++overall_unit_id;
  }

  const auto unit_block = get_unit_block( unit_id );
  const auto num_generators = unit_block->get_number_generators();

  for( Index g = 0 ; g < num_generators ; ++g , ++elc_generator ) {
   const auto zone = get_inertia_zone( elc_generator );
   affected_zones.insert( zone );
   ++num_generators_per_zone[ zone ];
  }

  ++overall_unit_id;
 }

 const auto number_nodes = get_number_nodes();

 // Now loop over all affected constraints

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {
  for( const auto zone_id : affected_zones ) {

   auto & constraint = v_InertiaDemand_Const[ t ][ zone_id ];

   // This will store the coefficients that must be updated, i.e., those of
   // the active Variables that belong to the units that have been modified.
   LinearFunction::Vec_FunctionValue coefficients;
   coefficients.reserve( num_generators_per_zone[ zone_id ] );

   // Subset that will store the indices of the active Variables whose
   // coefficients have changed.
   Subset subset;
   subset.reserve( num_generators_per_zone[ zone_id ] );

   for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {

    if( ! node_belongs_to_inertia_zone( node_id , zone_id ) )
     continue;

    Index elc_generator = 0;
    Index overall_unit_id = 0;
    for( const auto unit_id : modified_units ) {

     // Skip the units that have not been modified.
     while( overall_unit_id < unit_id ) {
      elc_generator += get_unit_block(
       overall_unit_id )->get_number_generators();
      ++overall_unit_id;
     }

     const auto unit_block = get_unit_block( unit_id );

     if( inertia_var_index[ unit_id ][ zone_id ] == Inf< Index >() ) {
      // This unit has no active Variable in the inertia demand constraints
      // associated with zone "zone_id".
      elc_generator += unit_block->get_number_generators();
      ++overall_unit_id;
      continue;
     }

     const auto scale = unit_block->get_scale();
     const auto num_generators = unit_block->get_number_generators();
     auto next_var_index = inertia_var_index[ unit_id ][ zone_id ];

     for( Index generator = 0 ;
          generator < num_generators ; ++generator , ++elc_generator ) {

      if( ! generator_belongs_to_node( elc_generator , node_id ) )
       continue;

      const auto commitment = unit_block->get_commitment( generator );
      auto inertia_commitment = unit_block->get_inertia_commitment( generator );

      if( commitment && inertia_commitment ) {

       assert( next_var_index < constraint.get_num_active_var() );
       assert( constraint.get_active_var( next_var_index )->get_Block()
               == unit_block );

       const auto coefficient = scale * inertia_commitment[ t ];
       coefficients.push_back( coefficient );
       subset.push_back( next_var_index++ );
      }

      const auto active_power = unit_block->get_active_power( generator );
      const auto inertia_power = unit_block->get_inertia_power( generator );

      if( active_power && inertia_power ) {
       assert( next_var_index < constraint.get_num_active_var() );
       assert( constraint.get_active_var( next_var_index )->get_Block()
               == unit_block );

       const auto coefficient = scale * inertia_power[ t ];
       coefficients.push_back( coefficient );
       subset.push_back( next_var_index++ );
      }

     }  // end( for( generator ) )

     ++overall_unit_id;

    }  // end( for( unit_id ) )

    // Update the coefficients of the active variables
    static_cast< LinearFunction * >( constraint.get_function() )->
     modify_coefficients( std::move( coefficients ) , std::move( subset ) ,
                          true , eNoBlck );

   }  // end( for( node_id ) )
  }  // end( for( zone_id ) )
 }  // end( for( t ) )

}  // end( UCBlock::update_inertia_demand_constraints )

/*--------------------------------------------------------------------------*/

void UCBlock::update_node_injection_constraints( Index time ,
                                                 Index node_index ,
                                                 double demand )
{
 auto rhs = demand;
 for( Index i = 0 ; i < f_number_units ; ++i ) {  // for each unit
  const auto unit_block = static_cast< UnitBlock * >( v_Block[ i ] );
  const auto scale = unit_block->get_scale();
  // for each electrical generator within the unit
  for( Index g = 0 ; g < unit_block->get_number_generators() ; ++g ) {
   if( auto fc = unit_block->get_fixed_consumption( g ) )
    if( fc[ time ] )
     if( auto u = unit_block->get_commitment( g ) ) {
      // add the contribution of the corresponding commitment variables
      rhs -= scale * fc[ time ];  // update the RHS
     }
  }  // end( for( g ) )
 }  // end( for( i ) )

 v_node_injection_Const[ time ][ node_index ].set_both( rhs , eNoBlck );
}

/*--------------------------------------------------------------------------*/

void UCBlock::set_active_power_demand( MF_dbl_it values ,
                                       Block::Subset && subset ,
                                       const bool ordered ,
                                       c_ModParam issuePMod ,
                                       c_ModParam issueAMod )
{
 if( subset.empty() )
  return;

 const auto number_nodes = get_number_nodes();

 if( ! v_network_blocks.empty() ) {
  // Update the demand of the NetworkBlocks
  // TODO Optimize
  for( auto index : subset ) {
   const auto node_index = index / f_time_horizon;
   const auto time = index % f_time_horizon;
   const auto demand = *values;
   v_network_blocks[ time ]->set_active_demand
    ( values++ , Range( node_index , node_index + 1 ) , issuePMod , issueAMod );

   if( number_nodes == 1 ) {
    assert( node_index == 0 );
    v_active_power_demand[ node_index ][ time ] = demand;
    update_node_injection_constraints( time , node_index , demand );
   }
  }
  return;
 }

 // Update the demand present in this UCBlock

 assert( ! v_active_power_demand.empty() );

 bool changed = false;

 for( auto index : subset ) {
  const auto node_index = index / f_time_horizon;
  const auto time = index % f_time_horizon;
  const auto demand = *(values++);

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

 // If nothing changes, return
 if( ! changed )
  return;

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< UCBlockSbstMod >(
                            this , UCBlockMod::eSetActD , std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( UCBlock::set_active_power_demand( subset ) )

/*--------------------------------------------------------------------------*/

void UCBlock::set_active_power_demand( MF_dbl_it values ,
                                       Block::Range rng ,
                                       c_ModParam issuePMod ,
                                       c_ModParam issueAMod )
{
 const auto number_nodes = get_number_nodes();

 rng.second = std::min( rng.second , number_nodes * f_time_horizon );

 if( rng.first >= rng.second )
  return;

 if( ! v_network_blocks.empty() ) {
  // Update the demand of the NetworkBlocks
  // TODO Optimize
  for( Index index = rng.first ; index < rng.second ; ++index ) {
   const auto node_index = index / f_time_horizon;
   const auto time = index % f_time_horizon;
   const auto demand = *values;
   v_network_blocks[ time ]->set_active_demand
    ( values++ , Range( node_index , node_index + 1 ) , issuePMod , issueAMod );

   if( number_nodes == 1 ) {
    assert( node_index == 0 );
    v_active_power_demand[ node_index ][ time ] = demand;
    update_node_injection_constraints( time , node_index , demand );
   }
  }
  return;
 }

 // Update the demand present in this UCBlock

 assert( ! v_active_power_demand.empty() );

 bool changed = false;

 for( Index index = rng.first ; index < rng.second ; ++index ) {
  const auto node_index = index / f_time_horizon;
  const auto time = index % f_time_horizon;
  const auto demand = *(values++);

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

 // If nothing changes, return
 if( ! changed )
  return;

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< UCBlockRngdMod >(
                            this , UCBlockMod::eSetActD , rng ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( UCBlock::set_active_power_demand( range ) )

/*--------------------------------------------------------------------------*/
/*------------------------ End File UCBlock.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
