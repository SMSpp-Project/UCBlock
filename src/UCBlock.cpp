/*--------------------------------------------------------------------------*/
/*----------------------- File UCBlock.cpp ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the UCBlock class.
 *
 * \version 0.20
 *
 * \date 28 - 01 - 2022
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Ali Ghezelsoflu \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Rafael Durbano Lobato \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Kostas Tavlaridis-Gyparakis \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy by Antonio Frangioni, Ali Ghezelsoflu, Rafael
 * Durbano Lobato, and Kostas Tavlaridis-Gyparakis
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

// TODO commented away until HeatBlock are properly managed
// #include "HeatBlock.h"

#include "LinearFunction.h"
#include "UCBlock.h"
#include "BusNetworkBlock.h"
#include "DCNetworkBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

SMSpp_insert_in_factory_cpp_1( UCBlock );

/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

UCBlock::~UCBlock()
{
 auto clear_constraints =
         []( boost::multi_array< FRowConstraint , 2 > & constraints ) {
          auto constraint = constraints.data();
          auto n = constraints.num_elements();
          for( decltype( n ) i = 0; i < n; ++i, ++constraint )
           constraint->clear();
         };

 clear_constraints( v_node_injection_constraints );
 clear_constraints( v_PrimaryDemand_Const );
 clear_constraints( v_SecondaryDemand_Const );
 clear_constraints( v_InertiaDemand_Const );
  /* TODO commented away until HeatBlock are properly managed
    clear_constraints( v_power_Heat_Rho_Const );
  */

 for( auto & v : v_PollutantBudget_Const )
  for( auto & constraint : v )
   constraint.clear();

 for( auto & block : v_Block )
  delete block;

 delete f_NetworkData;
 }

/*--------------------------------------------------------------------------*/

void UCBlock::deserialize_sub_blocks( const netCDF::NcGroup & group ,
                                      const std::string & prefix ,
                                      int num_sub_blocks )
{
 auto sz = v_Block.size();
 v_Block.resize( sz + num_sub_blocks , nullptr );
 for( int i = 0 ; i < num_sub_blocks ; ++i ) {
  std::string sub_group_name = prefix + std::to_string( i );
  auto sub_group = group.getGroup( sub_group_name );
  if( sub_group.isNull() )
   throw( std::invalid_argument( "UCBlock::deserialize: " +
                                 sub_group_name + " not present" ) );

  v_Block[ sz++ ] = new_Block( sub_group , this );
  }
 }

/*--------------------------------------------------------------------------*/

void UCBlock::deserialize_network_blocks( const netCDF::NcGroup & group )
{
 Index cntr = 0;
 v_network_blocks.resize( f_time_horizon , nullptr );

 for( Index i = 0 ; i < f_time_horizon ; ++i ) {
  std::string sub_group_name = "NetworkBlock_" + std::to_string( i );
  auto sub_group = group.getGroup( sub_group_name );
  if( sub_group.isNull() )
   continue;

  auto nbi = new_Block( sub_group , this );
  if( ( v_network_blocks[ i ] = dynamic_cast< NetworkBlock * >( nbi ) ) )
   ++cntr;
  else {
   delete nbi;
   throw( std::invalid_argument( sub_group_name +
                                 " not a valid NetworkBlock" ) );
   }
  }

 if( cntr ) {
  v_Block.resize( f_number_units + f_time_horizon );
  std::copy( v_network_blocks.begin() , v_network_blocks.end() ,
             std::next( v_Block.begin() , f_number_units ) );
  }
 else
  v_network_blocks.clear();
 }

/*--------------------------------------------------------------------------*/

void UCBlock::deserialize( const netCDF::NcGroup & group )
{
 #ifndef NDEBUG
  static std::vector< std::string > expected_dims = { "TimeHorizon" ,
   "NumberUnits" , "NumberHeatBlocks" , "NumberPrimaryZones" ,
   "NumberSecondaryZones" , "NumberInertiaZones" , "NumberPollutants" ,
   "NumberNodes" , "NumberLines" , "NumberElectricalGenerators" ,
   "TotalNumberPollutantZones" };
  check_dimensions( group , expected_dims , std::cerr );
  static std::vector< std::string > expected_vars = { "ActivePowerDemand" ,
   "GeneratorNode" , "HeatNode" , "HeatSet", "PowerHeatRho" , "PrimaryZones" ,
   "PrimaryDemand" , "SecondaryZones" , "SecondaryDemand" , "InertiaZones" ,
   "InertiaDemand" , "NumberPollutantZones" , "PollutantZones" ,
   "PollutantBudget" , "PollutantRho" , "StartLine" , "EndLine" ,
   "MinPowerFlow" , "MaxPowerFlow" , "Susceptance" , "NetworkCost" };
  check_variables( group, expected_vars, std::cerr );
 #endif

 ::deserialize_dim( group , "TimeHorizon" , f_time_horizon , false );
 ::deserialize_dim( group , "NumberUnits" , f_number_units , false );

 Index number_nodes;
 if( ! ::deserialize_dim( group , "NumberNodes" , number_nodes , true ) )
  number_nodes = 1;

 if( number_nodes > 1 ) {
  delete f_NetworkData;
  f_NetworkData = new NetworkBlock::NetworkData();
  f_NetworkData->deserialize( group );
  }

 /* TODO commented away until HeatBlock are properly managed
 if( ! ::deserialize_dim( group , "NumberHeatGenerators" ,
                          f_number_heat_generators , true ) )
  f_number_heat_generators = 0;
 */

 auto ActivePowerDemand = group.getVar( "ActivePowerDemand" );
 if( ! ActivePowerDemand.isNull() ) {
  // ActivePowerDemand has been provided.
  using index = decltype( v_active_power_demand )::index;
  std::vector< index > shape = { number_nodes , f_time_horizon };
  v_active_power_demand.resize( shape );
  ActivePowerDemand.getVar( v_active_power_demand.data() );
  }

 // optional dimensions
 /* TODO commented away until HeatBlock are properly managed
 f_number_heat_blocks = 0;
 ::deserialize_dim( group , "NumberHeatBlocks" ,
                    f_number_heat_blocks , true );
 */

 f_number_primary_zones = 0;
 ::deserialize_dim( group , "NumberPrimaryZones" ,
                    f_number_primary_zones , true );

 f_number_secondary_zones = 0;
 ::deserialize_dim( group , "NumberSecondaryZones" ,
                    f_number_secondary_zones , true );

 f_number_inertia_zones = 0;
 ::deserialize_dim( group , "NumberInertiaZones" ,
                    f_number_inertia_zones , true );

 f_number_pollutants = 0;
 ::deserialize_dim( group , "NumberPollutants" ,
                    f_number_pollutants , true );

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
                          f_total_number_pollutant_zones , true ) ) {
  f_total_number_pollutant_zones = 0;
  for( const auto & n : v_number_pollutant_zones )
   f_total_number_pollutant_zones += n;
  }

 if( ! f_total_number_pollutant_zones )
  f_total_number_pollutant_zones = f_number_pollutants;

 if( f_total_number_pollutant_zones ) {
  ::deserialize( group , "PollutantZones" ,
                 v_pollutant_zones , true , true );

  ::deserialize( group , "PollutantBudget" ,
                 { f_total_number_pollutant_zones } ,
                 v_pollutant_budget , true , false );

 ::deserialize( group , "PollutantRho" ,
                v_pollutant_rho , true , true );

  /* TODO commented away until HeatBlock are properly managed
  ::deserialize( group , "PollutantHeatRho" ,
                 v_pollutant_heat_rho , true , true );
  */
  }
 else {
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
  delete block;
 v_Block.clear();

 // load all UnitBlock
 deserialize_sub_blocks( group , "UnitBlock_" , f_number_units );

 // Generate UnitBlock primary spinning reserve variables
 unsigned int what = 0;
 if( f_number_primary_zones > 0 ) {
  what += 1;
 }

 // Generate UnitBlock secondary spinning reserve variables
 if( f_number_secondary_zones > 0 ) {
  what += 2;
 }

 // Generate UnitBlock inertia reserve variables
 if( f_number_inertia_zones > 0 ) {
  what += 4;
 }

 if( what > 0 ) {
  for( auto * b: v_Block ) {
   if( auto ub = dynamic_cast<UnitBlock *>(b) ) {
    ub->set_reserve_vars(what);
   }
  }
 }

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
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    if( v_network_blocks[ t ] ) {
     auto & ad = v_network_blocks[ t ]->get_active_demand();
     if( ! ad.empty() )
      v_active_power_demand[ 0 ][ t ] = ad.front();
     delete v_network_blocks[ t ];
     }
   v_network_blocks.clear();
   v_Block.resize( f_number_units );
   }
  }
 else {  // number_nodes > 1
  // if they don't exist, create them now as DCNetworkBlock
  if( v_network_blocks.empty() ) {
   v_network_blocks.resize( f_time_horizon , nullptr );
   v_Block.resize( f_number_units + f_time_horizon , nullptr );
   }

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   auto nbi = v_network_blocks[ t ];
   if( ! nbi ) {  // NetworkBlock i does not exist: create a DCNetworkBlock
    nbi = new DCNetworkBlock( this );
    v_network_blocks[ t ] = nbi;
    v_Block[ f_number_units + t ] = nbi;
    }

   if( ! nbi->get_NetworkData() ) {
    if( ! f_NetworkData )
     throw( std::invalid_argument( "UCBlock::deserialize: NetworkData "
                                   "missing in NetworkBlock " +
                                   std::to_string( t ) + " and in UCBlock" )
           );
     nbi->set_NetworkData( f_NetworkData );
    }

   if( nbi->get_active_demand().empty() ) {
    if( ! v_active_power_demand.num_elements() )
     throw( std::invalid_argument(
      "UCBlock::deserialize: ActivePowerDemand missing in UCBlock and in "
      "NetworkBlock " + std::to_string( t ) ) );
    typedef boost::multi_array_types::index_range range;
    auto ap_c = v_active_power_demand[
                        boost::indices[ range( 0 , number_nodes ) ][ t ] ];
    std::vector< double > ap_v( number_nodes );
    std::copy( ap_c.begin() , ap_c.end() , ap_v.begin() );
    nbi->set_ActiveDemand( ap_v );
    }
   }  // end( for( t ) )

  // v_active_power_demand used up, disband it
  v_active_power_demand.resize(
                 boost::multi_array< double , 2 >::extent_gen()[ 0 ][ 0 ] );

  }  // end( else( number_nodes > 1 ) )

 /* TODO commented away until HeatBlock are properly managed
 if( f_number_heat_blocks )
  deserialize_sub_blocks( group , "HeatBlock_" , f_number_heat_blocks );
 */

 if( !::deserialize_dim( group , "NumberElectricalGenerators" ,
                         f_number_elc_generators , true ) ) {
  f_number_elc_generators = 0;
  for( Index i = 0 ; i < f_number_units ; ++i )
   f_number_elc_generators +=
    static_cast< UnitBlock * >( v_Block[ i ] )->get_number_generators();
  }

 ::deserialize( group , "GeneratorNode" , f_number_elc_generators ,
                v_generator_node , true , true );

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

 // mark all done

 set_constraints_generated();

}  // end( UCBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void UCBlock::generate_node_injection_constraints() {

 const auto number_nodes = get_number_nodes();

 v_node_injection_constraints.resize(
  boost::multi_array< FRowConstraint , 2 >::extent_gen()[ f_time_horizon ]
                                                        [ number_nodes ] );

 if( number_nodes > 0 ) {  // well, that'd be curious, but ...
  if( number_nodes == 1 ) {
   // special case: in a BusNetwork there are no NetworkBlocks and the node
   // injection constraints actually are active power demand constraints
   //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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
      *(vcit++) = std::pair( & unit_block->get_active_power( g )[ t ] , scale );

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
         *(vcit++) = std::pair( & u[ t ] , - fixed_consumption );
         rhs -= fixed_consumption;    // update the RHS
         }
      }  // end( for( g ) )
     }  // end( for( i ) )

    // set the final RHS of the constraint (equality constraint)
    v_node_injection_constraints[ t ][ 0 ].set_both( rhs , eNoMod );
    // resize vc so that it's of the right length
    vc.resize( std::distance( vc.begin() , vcit ) );
    // construct and pass the LinearFunction to the FRowConstraint
    v_node_injection_constraints[ t ][ 0 ].set_function(
                          new LinearFunction( std::move( vc ) ) , eNoMod  );
    }  // end( for( t ) )
   }
  else {  // number_nodes > 1
   // DCNetwork needs GeneratorNode

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {

    auto & node_injection = v_network_blocks[ t ]->get_node_injection();

    for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {

     auto linear_function = new LinearFunction();

     linear_function->add_variable( &node_injection[node_id] , -1.0 , eNoMod );

     v_node_injection_constraints[ t ][ node_id ].set_both( 0.0 );

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
        linear_function->add_variable( active_power , scale , eNoMod );
       }

       double fixed_consumption = 0.0;
       if( auto fc = unit_block->get_fixed_consumption( generator ) )
        fixed_consumption = fc[ t ] * scale;

       if( auto c = unit_block->get_commitment( generator ) ) {
        auto commitment = &c[ t ];
        linear_function->add_variable( commitment , - fixed_consumption ,
                                       eNoMod );
       }

       v_node_injection_constraints[ t ][ node_id ].set_both
        ( v_node_injection_constraints[ t ][ node_id ].get_rhs()
          - fixed_consumption );
      }
     }
     v_node_injection_constraints[t][node_id].set_function( linear_function );
    }
   }
  }
  add_static_constraint( v_node_injection_constraints , "node_injection_c" );
 }
}  // end( UCBlock::generate_node_injection_constraints )

/*--------------------------------------------------------------------------*/

void UCBlock::generate_primary_demand_constraints() {

 if( f_number_primary_zones == 0 )
  return;

 v_PrimaryDemand_Const.resize
  ( boost::multi_array< FRowConstraint , 2 >::
    extent_gen()[ f_time_horizon ][ f_number_primary_zones ] );

 const auto number_nodes = get_number_nodes();

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {
  for( Index zone_id = 0 ; zone_id < f_number_primary_zones ; ++zone_id ) {

   auto linear_function = new LinearFunction();

   for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {

    if( ( f_number_primary_zones > 1 ) &&
        ( zone_id != v_primary_zones[ node_id ] ) )
     continue;

    Index elc_generator = 0;
    for( Index unit_id = 0 ; unit_id < f_number_units ; ++unit_id ) {

     const auto unit_block = get_unit_block( unit_id );
     const auto scale = unit_block->get_scale();

     for( Index generator = 0 ;
          generator < unit_block->get_number_generators() ;
          ++generator , ++elc_generator ) {

      if( ( number_nodes > 1 ) &&
          ( node_id != v_generator_node[ elc_generator ] ) )
       continue;

      if( auto primary_s_r =
          unit_block->get_primary_spinning_reserve( generator ) ) {
       auto primary_spinning_reserve = & primary_s_r[ t ];
       linear_function->add_variable( primary_spinning_reserve , scale );
      }

     }  // end( for( generator ) )
    }  // end( for( unit_id ) )
   }  // end( for( node_id ) )

   const auto demand = get_primary_demand()[ zone_id ][ t ];
   v_PrimaryDemand_Const[ t ][ zone_id ].set_lhs( demand );
   v_PrimaryDemand_Const[ t ][ zone_id ].set_rhs( Inf< double >() );
   v_PrimaryDemand_Const[ t ][ zone_id ].set_function( linear_function );
  }  // end( for( zone_id ) )
 }  // end( for( t ) )

 add_static_constraint( v_PrimaryDemand_Const , "primary_demand_c" );

}  // end( UCBlock::generate_primary_demand_constraints )

/*--------------------------------------------------------------------------*/

void UCBlock::generate_secondary_demand_constraints() {

 const auto number_nodes = get_number_nodes();

 if( f_number_secondary_zones > 0 ) {

  v_SecondaryDemand_Const.resize
   ( boost::multi_array< FRowConstraint, 2 >::
     extent_gen()[ f_time_horizon ][ f_number_secondary_zones ] );

  if( f_number_secondary_zones == 1 ) {  //no need to SecondaryZones

   if( number_nodes == 1 ) {  //BusNetwork no need to GeneratorNode

    for( Index t = 0 ; t < f_time_horizon ; ++t ) {
     auto linear_function = new LinearFunction();

     for( Index unit_id = 0 ; unit_id < f_number_units ; unit_id++ ) {

      const auto unit_block = get_unit_block( unit_id );
      const auto scale = unit_block->get_scale();

      for( Index generator = 0 ;
           generator < unit_block->get_number_generators() ; ++generator ) {

       if( auto secondary_s_r =
           unit_block->get_secondary_spinning_reserve( generator ) ) {
        auto secondary_spinning_reserve = & secondary_s_r[ t ];
        linear_function->add_variable( secondary_spinning_reserve , scale );
       }
      }
     }
     v_SecondaryDemand_Const[ t ][ 0 ].set_lhs( get_secondary_demand()[0][t] );
     v_SecondaryDemand_Const[ t ][ 0 ].set_rhs( Inf< double >());
     v_SecondaryDemand_Const[ t ][ 0 ].set_function( linear_function );
    }
   }
   else {  //DCNetwork needs GeneratorNode

    for( Index t = 0 ; t < f_time_horizon ; ++t ) {

     auto linear_function = new LinearFunction();

     for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {

      Index elc_generator = 0;
      for( Index unit_id = 0 ; unit_id < f_number_units ; ++unit_id ) {

       const auto unit_block = get_unit_block( unit_id );
       const auto scale = unit_block->get_scale();

       for( Index generator = 0 ;
            generator < unit_block->get_number_generators() ;
            ++generator , ++elc_generator ) {

        if( node_id != v_generator_node[ elc_generator ] )
         continue;

        if( auto secondary_s_r =
            unit_block->get_secondary_spinning_reserve( generator ) ) {
         auto secondary_spinning_reserve = & secondary_s_r[ t ];
         linear_function->add_variable( secondary_spinning_reserve , scale );
        }
       }
      }
     }
     v_SecondaryDemand_Const[ t ][ 0 ].set_lhs( get_secondary_demand()[0][t] );
     v_SecondaryDemand_Const[ t ][ 0 ].set_rhs( Inf< double >() );
     v_SecondaryDemand_Const[ t ][ 0 ].set_function( linear_function );
    }
   }
  }
  else if( f_number_secondary_zones > 1 ) {   // SecondaryZones is needed

   if( number_nodes == 1 ) {  //BusNetwork no need to GeneratorNode

    for( Index t = 0 ; t < f_time_horizon ; ++t ) {

     for( Index zone_id = 0 ; zone_id < f_number_secondary_zones ; ++zone_id ) {

      auto linear_function = new LinearFunction();

      if( zone_id == v_secondary_zones[ 0 ] ) {

       for( Index unit_id = 0 ; unit_id < f_number_units ; ++unit_id ) {

        const auto unit_block = get_unit_block( unit_id );
        const auto scale = unit_block->get_scale();

        for( Index generator = 0 ;
             generator < unit_block->get_number_generators() ; ++generator ) {

         if( auto secondary_s_r =
             unit_block->get_secondary_spinning_reserve( generator ) ) {
          auto secondary_spinning_reserve = & secondary_s_r[ t ];
          linear_function->add_variable( secondary_spinning_reserve , scale );
         }
        }
       }
      }

      v_SecondaryDemand_Const[ t ][ zone_id ].set_lhs
       ( get_secondary_demand()[ zone_id ][ t ] );
      v_SecondaryDemand_Const[ t ][ zone_id ].set_rhs( Inf< double >() );
      v_SecondaryDemand_Const[ t ][ zone_id ].set_function( linear_function );
     }
    }
   }
   else {  //DCNetwork needs GeneratorNode

    for( Index t = 0 ; t < f_time_horizon ; ++t ) {

     for( Index zone_id = 0 ; zone_id < f_number_secondary_zones ; ++zone_id ) {

      auto linear_function = new LinearFunction();

      for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {
       if( zone_id == v_secondary_zones[ node_id ] ) {

        Index elc_generator = 0;
        for( Index unit_id = 0 ; unit_id < f_number_units ; unit_id++ ) {

         const auto unit_block = get_unit_block( unit_id );
         const auto scale = unit_block->get_scale();

         for( Index generator = 0 ;
              generator < unit_block->get_number_generators() ;
              ++generator , ++elc_generator ) {

          if( node_id != v_generator_node[ elc_generator ] )
           continue;

          if( auto secondary_s_r =
              unit_block->get_secondary_spinning_reserve( generator ) ) {

           auto secondary_spinning_reserve = & secondary_s_r[ t ];
           linear_function->add_variable( secondary_spinning_reserve , scale );
          }
         }
        }
       }
      }
      v_SecondaryDemand_Const[ t ][ zone_id ].set_lhs
       ( get_secondary_demand()[ zone_id ][ t ] );
      v_SecondaryDemand_Const[ t ][ zone_id ].set_rhs( Inf< double >() );
      v_SecondaryDemand_Const[ t ][ zone_id ].set_function( linear_function );
     }
    }
   }
  }
  add_static_constraint( v_SecondaryDemand_Const , "secondary_demand_c" );
 }
}  // end( UCBlock::generate_secondary_demand_constraints )

/*--------------------------------------------------------------------------*/

void UCBlock::generate_inertia_demand_constraints() {

 const auto number_nodes = get_number_nodes();

 if( f_number_inertia_zones > 0 ) {

  v_InertiaDemand_Const.resize
   ( boost::multi_array< FRowConstraint, 2 >::
     extent_gen()[ f_time_horizon ][ f_number_inertia_zones ] );

  if( f_number_inertia_zones == 1 ) {  //no need to InertiaZones

   if( number_nodes == 1 ) {  //BusNetwork no need to GeneratorNode

    for( Index t = 0 ; t < f_time_horizon ; ++t ) {

     auto linear_function = new LinearFunction();

     v_InertiaDemand_Const[ t ][ 0 ].set_lhs( get_inertia_demand()[ 0 ][ t ] );
     v_InertiaDemand_Const[ t ][ 0 ].set_rhs( Inf< double >() );

     for( Index unit_id = 0 ; unit_id < f_number_units ; unit_id++ ) {

      const auto unit_block = get_unit_block( unit_id );
      const auto scale = unit_block->get_scale();

      for( Index generator = 0 ;
           generator < unit_block->get_number_generators() ; ++generator ) {

       auto c = unit_block->get_commitment( generator );
       auto inertia_commitment =
        unit_block->get_inertia_commitment( generator );

       if( c && inertia_commitment ) {
        auto commitment = & c[ t ];
        auto coefficient = scale * inertia_commitment[ t ];
        linear_function->add_variable( commitment , coefficient );
       }

       auto ap = unit_block->get_active_power( generator );
       auto inertia_power = unit_block->get_inertia_power( generator );

       if( ap && inertia_power ) {
        auto active_power = & ap[ t ];
        auto coefficient = scale * inertia_power[ t ];
        linear_function->add_variable( active_power , coefficient );
       }
      }
     }
     v_InertiaDemand_Const[ t ][ 0 ].set_function( linear_function );
    }
   }
   else {  //DCNetwork needs GeneratorNode

    for( Index t = 0 ; t < f_time_horizon ; ++t ) {

     auto linear_function = new LinearFunction();

     v_InertiaDemand_Const[ t ][ 0 ].set_lhs( get_inertia_demand()[ 0 ][ t ] );
     v_InertiaDemand_Const[ t ][ 0 ].set_rhs( Inf< double >() );

     for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {

      Index elc_generator = 0;
      for( Index unit_id = 0 ; unit_id < f_number_units ; ++unit_id ) {

       const auto unit_block = get_unit_block( unit_id );
       const auto scale = unit_block->get_scale();

       for( Index generator = 0 ;
            generator < unit_block->get_number_generators() ;
            ++generator , ++elc_generator ) {

        if( node_id != v_generator_node[ elc_generator ] )
         continue;

        auto c = unit_block->get_commitment( generator );
        auto inertia_commitment =
         unit_block->get_inertia_commitment( generator );

        if( c && inertia_commitment ) {
         auto commitment = & c[ t ];
         auto coefficient = scale * inertia_commitment[ t ];
         linear_function->add_variable( commitment , coefficient );
        }

        auto ap = unit_block->get_active_power( generator );
        auto inertia_power = unit_block->get_inertia_power( generator );
        if( ap && inertia_power ) {
         auto active_power = & ap[ t ];
         auto coefficient = scale * inertia_power[ t ];
         linear_function->add_variable( active_power, coefficient );
        }
       }
      }
     }
     v_InertiaDemand_Const[ t ][ 0 ].set_function( linear_function );
    }
   }
  }
  else if( f_number_inertia_zones > 1 ) {   // InertiaZones is needed

   if( number_nodes == 1 ) {  //BusNetwork no need to GeneratorNode

    for( Index t = 0 ; t < f_time_horizon ; ++t ) {

     for( Index zone_id = 0 ; zone_id < f_number_inertia_zones ; ++zone_id ) {

      auto linear_function = new LinearFunction();

      if( zone_id == v_inertia_zones[ 0 ] ) {

       for( Index unit_id = 0 ; unit_id < f_number_units ; ++unit_id ) {

        const auto unit_block = get_unit_block( unit_id );
        const auto scale = unit_block->get_scale();

        for( Index generator = 0 ;
             generator < unit_block->get_number_generators() ; ++generator ) {

         auto c = unit_block->get_commitment( generator );
         auto inertia_commitment =
          unit_block->get_inertia_commitment( generator );

         if( c && inertia_commitment ) {
          auto commitment = & c[ t ];
          auto coefficient = scale * inertia_commitment[ t ];
          linear_function->add_variable( commitment , coefficient );
         }

         auto ap = unit_block->get_active_power( generator );
         auto inertia_power = unit_block->get_inertia_power( generator );

         if( ap && inertia_power ) {
          auto active_power = & ap[ t ];
          auto coefficient = scale * inertia_power[ t ];
          linear_function->add_variable( active_power , coefficient );
         }
        }
       }
      }
      v_InertiaDemand_Const[ t ][ zone_id ].set_lhs
       ( get_inertia_demand()[ zone_id ][ t ] );
      v_InertiaDemand_Const[ t ][ zone_id ].set_rhs( Inf< double >() );
      v_InertiaDemand_Const[ t ][ zone_id ].set_function( linear_function );
     }
    }
   }
   else {  //DCNetwork needs GeneratorNode

    for( Index t = 0 ; t < f_time_horizon ; ++t ) {

     for( Index zone_id = 0 ; zone_id < f_number_inertia_zones ; ++zone_id ) {

      auto linear_function = new LinearFunction();

      for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {
       if( zone_id == v_inertia_zones[ node_id ] ) {

        Index elc_generator = 0;
        for( Index unit_id = 0 ; unit_id < f_number_units ; ++unit_id ) {

         const auto unit_block = get_unit_block( unit_id );
         const auto scale = unit_block->get_scale();

         for( Index generator = 0 ;
              generator < unit_block->get_number_generators() ;
              ++generator , ++elc_generator ) {

          if( node_id != v_generator_node[ elc_generator ] )
           continue;

          auto c = unit_block->get_commitment( generator );
          auto inertia_commitment =
           unit_block->get_inertia_commitment( generator );

          if( c && inertia_commitment ) {
           auto commitment = & c[ t ];
           auto coefficient = scale * inertia_commitment[ t ];
           linear_function->add_variable( commitment , coefficient );
          }

          auto ap = unit_block->get_active_power( generator );
          auto inertia_power = unit_block->get_inertia_power( generator );

          if( ap && inertia_power ) {
           auto active_power = & ap[ t ];
           auto coefficient = scale * inertia_power[ t ];
           linear_function->add_variable( active_power , coefficient );
          }
         }
        }
       }
      }
      v_InertiaDemand_Const[ t ][ zone_id ].set_lhs
       ( get_inertia_demand()[ zone_id ][ t ] );
      v_InertiaDemand_Const[ t ][ zone_id ].set_rhs( Inf< double >() );
      v_InertiaDemand_Const[ t ][ zone_id ].set_function( linear_function );
     }
    }
   }
  }
  add_static_constraint( v_InertiaDemand_Const , "inertia_demand_c" );
 }
}  // end( UCBlock::generate_inertia_demand_constraints )

/*--------------------------------------------------------------------------*/

void UCBlock::generate_pollutant_budget_constraints() {

 // TODO These constraints must be fixed

 const auto number_nodes = get_number_nodes();

 if( f_number_pollutants > 0 ) {

  v_PollutantBudget_Const.resize
   ( v_number_pollutant_zones[ f_total_number_pollutant_zones ] );

  if( number_nodes == 1 ) {

   for( Index pollutant = 0 ; pollutant < f_number_pollutants ; ++pollutant ) {

    for( Index zone = 0 ; zone < v_number_pollutant_zones[ pollutant ] ;
         ++zone ) {

     for( Index t = 0 ; t < f_time_horizon ; ++t ) {

      // Terms associated with active power
      auto linear_function = new LinearFunction();

      for( Index unit_id = 0 ; unit_id < f_number_units ; ++unit_id ) {

       const auto unit_block = get_unit_block( unit_id );
       const auto scale = unit_block->get_scale();

       for( Index generator = 0 ;
            generator < unit_block->get_number_generators() ; ++generator ) {

        auto node_id = get_generator_node()[ generator ];
        auto zone_id = get_pollutant_zone()[ pollutant ][ node_id ];

        if( zone_id >= v_number_pollutant_zones[ pollutant ] )
         continue; // this unit does not belong to any zone

        if( auto ap = unit_block->get_active_power( generator ) ) {
         auto active_power = & ap[ t ];
         auto rho = get_pollutant_rho()[ t ][ pollutant ][ generator ];
         auto coefficient =  scale * rho;
         linear_function->add_variable( active_power , coefficient );
        }
       }
      }

       // Terms associated with heat-only generation units
       /* TODO commented away until HeatBlock are properly managed
       if( f_number_heat_blocks > 0 ) { //TODO Do we have any HeatBlock?

        for( Index h = 0; h < f_number_heat_blocks; ++h ) {

         for( std::vector< Index >::size_type i = 0;
              i <= f_number_units; ++i ) {

          //auto unit_id = v_heat_only_units[ i ];
          auto heat_id = v_heat_set[i];

          auto zone_id = get_pollutant_zone()[pollutant][h];
          if( zone_id >= v_number_pollutant_zones[pollutant] )
           continue; // this unit does not belong to any zone

          auto heat = get_heat_block()[h]->get_heat()[t][i];
          auto rho = get_pollutant_heat_rho()[t][pollutant][h];

          auto linear_function = dynamic_cast<LinearFunction *>
          ( v_PollutantBudget_Const[pollutant][zone_id].get_function());

          linear_function->add_variable( &heat, rho );
         }
        }
       }
       */

      v_PollutantBudget_Const[ pollutant ][ zone ].set_rhs
       ( v_pollutant_budget[v_number_pollutant_zones[pollutant]][ pollutant ] );
      v_PollutantBudget_Const[ pollutant ][ zone ].set_lhs( -Inf< double >() );
      v_PollutantBudget_Const[ pollutant ][ zone ].set_function
       ( linear_function );
     }
    }
   }
  }
  else { // DCNetwork

   for( Index pollutant = 0 ; pollutant < f_number_pollutants ; ++pollutant ) {

    for( Index zone = 0 ; zone < v_number_pollutant_zones[ pollutant ] ;
         ++zone ) {

     for( Index t = 0 ; t < f_time_horizon ; ++t ) {

      // Terms associated with active power
      auto linear_function = new LinearFunction();

      Index pollutant_zone = 0;
      for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {

       if( zone == v_pollutant_zones[ pollutant ][ node_id ] ) {

        Index generator_id = 0;
        for( Index elc_generator = 0 ; elc_generator < f_number_elc_generators ;
             ++elc_generator ) {

         if( node_id == v_generator_node[ elc_generator ] ) {

          auto block = get_nested_Blocks()[ generator_id ];
          auto unit_block = dynamic_cast<UnitBlock *>(block);
          if( ! unit_block )
           continue;

          const auto scale = unit_block->get_scale();

          for( Index generator = 0;
               generator < unit_block->get_number_generators(); ++generator ) {

           auto node_id = get_generator_node()[generator];
           auto zone_id = get_pollutant_zone()[pollutant][node_id];

           if( zone_id >= v_number_pollutant_zones[pollutant] )
            continue; // this unit does not belong to any zone

           if( auto ap = unit_block->get_active_power( generator ) ) {
            auto active_power = &ap[t];
            auto rho = get_pollutant_rho()[t][pollutant][generator];
            auto coefficient = scale * rho;
            linear_function->add_variable( active_power , coefficient );
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
       if( f_number_heat_blocks > 0 ) { //TODO Do we have any HeatBlock?

        for( Index h = 0; h < f_number_heat_blocks; ++h ) {

         for( std::vector< Index >::size_type i = 0;
              i <= f_number_units; ++i ) {

          //auto unit_id = v_heat_only_units[ i ];
          auto heat_id = v_heat_set[i];

          auto zone_id = get_pollutant_zone()[pollutant][h];
          if( zone_id >= v_number_pollutant_zones[pollutant] )
           continue; // this unit does not belong to any zone

          auto heat = get_heat_block()[h]->get_heat()[t][i];
          auto rho = get_pollutant_heat_rho()[t][pollutant][h];

          auto linear_function = dynamic_cast<LinearFunction *>
          ( v_PollutantBudget_Const[pollutant][zone_id].get_function());

          linear_function->add_variable( &heat, rho );
         }
        }
       }
       */

      v_PollutantBudget_Const[ pollutant ][ zone ].set_rhs
       ( v_pollutant_budget[v_number_pollutant_zones[pollutant]][ pollutant ] );
      v_PollutantBudget_Const[ pollutant ][ zone ].set_lhs( -Inf< double >() );
      v_PollutantBudget_Const[ pollutant ][ zone ].set_function
       ( linear_function );
     }
    }
   }
  }
  add_static_constraint
   ( v_PollutantBudget_Const[ f_total_number_pollutant_zones ] );
 }
}  // end( UCBlock::generate_pollutant_budget_constraints )

/*--------------------------------------------------------------------------*/

void UCBlock::generate_heat_constraints() {

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
       return true;
     }
     return false;
    };

   for( Index unit_id = 0; unit_id < f_number_units; ++unit_id ) {
    if( is_electricity_producing_inside_a_heat_block( unit_id ) ) {
     electricity_generators_inside_a_heat_block
     [ num_constraints_per_time++ ] = unit_id;
    }
   }

   v_power_Heat_Rho_Const.resize
    ( boost::multi_array< FRowConstraint *, 2 >::
      extent_gen()[ f_time_horizon ][ num_constraints_per_time ] );
  }

  for( Index t = 0; t < f_time_horizon; ++t ) {

   for( std::vector< Index >::size_type constraint_id = 0;
        constraint_id < num_constraints_per_time; ++constraint_id ) {

    auto generator_id = electricity_generators_inside_a_heat_block[ constraint_id ];

    for( Index heat_block_id = 0; heat_block_id < f_number_heat_blocks;
         ++heat_block_id ) {

     auto heat_unit_id = get_heat_set()[ generator_id ];

     if( heat_unit_id >=
         v_heat_blocks[ heat_block_id ]->get_number_heat_generators() )
      continue; // unit_id does not belong to heat_block_id

     if( !v_power_Heat_Rho_Const[ t ][ constraint_id ].get_function() ) {
      v_power_Heat_Rho_Const[ t ][ constraint_id ].
       set_lhs( -Inf< double >() );
      v_power_Heat_Rho_Const[ t ][ constraint_id ].set_rhs( 0.0 );
      v_power_Heat_Rho_Const[ t ][ constraint_id ].set_function
       ( new LinearFunction() );
     }

     auto active_power = get_unit_block( generator_id )->get_active_power();
     auto heat = get_heat_block() [ heat_unit_id ]->get_heat()[ t ][generator_id];
     auto power_heat_rho = get_power_heat_rho()[ generator_id ];

     auto linear_function = dynamic_cast<LinearFunction *>
     ( v_power_Heat_Rho_Const[ t ][ constraint_id ].get_function());
     linear_function->add_variable( &heat, 1, 0 );
     linear_function->add_variable( &active_power[ t ][ generator_id ], -power_heat_rho );
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

void UCBlock::serialize( netCDF::NcGroup & group ) const {

 Block::serialize( group );

 netCDF::NcDim NumberNodes;

 if( f_NetworkData ) {
  f_NetworkData->serialize( group );
  NumberNodes = group.getDim( "NumberNodes" );
 }
 else {
  NumberNodes = group.addDim( "NumberNodes" , 1 );
 }

 auto TimeHorizon = group.addDim( "TimeHorizon" , f_time_horizon );
 auto NumberUnits = group.addDim( "NumberUnits" , f_number_units );
 auto NumberElectricalGenerators = group.addDim( "NumberElectricalGenerators" ,
                                                 f_number_elc_generators );

 auto TotalNumberPollutantZones = group.addDim
  ( "TotalNumberPollutantZones" , f_total_number_pollutant_zones );

 /* TODO commented away until HeatBlock are properly managed
 auto NumberHeatBlocks =
  group.addDim( "NumberHeatBlocks", f_number_heat_blocks );

 auto NumberHeatGenerators = group.addDim( "NumberHeatGenerators" ,
                                           f_number_heat_generators );
 */

 auto NumberPrimaryZones =
  group.addDim( "NumberPrimaryZones" , f_number_primary_zones );
 auto NumberSecondaryZones =
  group.addDim( "NumberSecondaryZones" , f_number_secondary_zones );
 auto NumberInertiaZones =
  group.addDim( "NumberInertiaZones" , f_number_inertia_zones );
 auto NumberPollutants =
  group.addDim( "NumberPollutants" , f_number_pollutants );

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

 ::serialize( group , "PollutantBudget" , netCDF::NcDouble() ,
              TotalNumberPollutantZones , v_pollutant_budget );

 ::serialize( group , "PollutantRho" , netCDF::NcDouble() ,
              { TimeHorizon , NumberPollutants , NumberElectricalGenerators } ,
              v_pollutant_rho );

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

 for( Index i = 0; i < f_number_units; ++i ) {
  auto sub_block = get_unit_block( i );
  auto sub_group = group.addGroup( "UnitBlock_" + std::to_string( i ));
  sub_block->serialize( sub_group );
 }

 for( Index t = 0; t < f_time_horizon; ++t )
  if( auto sub_block = get_network_block( t ) ) {
   auto sub_group = group.addGroup( "NetworkBlock_" + std::to_string( t ) );
   sub_block->serialize( sub_group );
   }

 /* TODO commented away until HeatBlock are properly managed
 for( Index i = 0; i < f_number_heat_blocks; ++i ) {
  auto sub_block = get_heat_block( i );
  auto sub_group = group.addGroup( "HeatBlock_" + std::to_string( i ));
  sub_block->serialize( sub_group );
  }
 */

 }  // end( UCBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

void UCBlock::add_Modification( sp_Mod mod , ChnlName chnl ) {
 if( mod->concerns_Block() ) {
  if( const auto tmod = dynamic_cast< UnitBlockMod * >( mod.get() ) ) {
   mod->concerns_Block( false );
   if( tmod->type() == UnitBlockMod::eScale ) {
    update_node_injection_constraints_scale( tmod->get_Block() );
   }
  }
 }

 Block::add_Modification( mod , chnl );
}

/*--------------------------------------------------------------------------*/

void UCBlock::update_node_injection_constraints_scale( Block * block ) {

 if( ( ! constraints_generated() ) ||
     ( v_node_injection_constraints.size() == 0 ) )
  return;

 const auto given_unit_block = dynamic_cast< UnitBlock * >( block );
 if( ! given_unit_block )
  throw( std::invalid_argument( "UCBlock::update_node_injection_constraints_"
                                "scale: given Block is not a UnitBlock." ) );

 const auto number_nodes = get_number_nodes();

 if( number_nodes > 0 ) {
  if( number_nodes == 1 ) { // BusNetwork
   for( Index t = 0 ; t < f_time_horizon ; ++t ) {  // for each time instant

    /* The active Variables of the LinearFunction defining the constraint are
     * grouped by UnitBlocks. That is, all active Variables of a given
     * UnitBlock have consecutive indices in this LinearFunction. The
     * following will store the Range of indices of the active Variables of
     * this LinearFunction that belong to the given Block. */
    Range range( Inf< Index >() , Inf< Index >() );

    // This will store the coefficients that must be updated, i.e., those of
    // the active Variables that belong to the given Block.
    LinearFunction::Vec_FunctionValue coefficients;

    Index active_var_index = 0;

    // initialise demand as active power
    auto rhs = v_active_power_demand[ 0 ][ t ];

    for( Index i = 0 ; i < f_number_units ; ++i ) {  // for each unit
     const auto unit_block = get_unit_block( i );
     const auto scale = unit_block->get_scale();

     if( unit_block == given_unit_block )
      coefficients.reserve( 2 * unit_block->get_number_generators() );

     // for each electrical generator within the unit
     for( Index g = 0 ; g < unit_block->get_number_generators() ; ++g ) {

      if( unit_block == given_unit_block ) {

       // update the Range
       if( range.first == Inf< Index >() ) {
        range.first = active_var_index;
        range.second = active_var_index;
       }
       range.second++;

       // update the coefficient of the active power variable
       coefficients.push_back( scale );
      }

      // increment due to the active power variable
      ++active_var_index;

      if( auto fc = unit_block->get_fixed_consumption( g ) )
       if( fc[ t ] )
        if( unit_block->get_commitment( g ) ) {
         const auto fixed_consumption = fc[ t ] * scale;
         rhs -= fixed_consumption;    // update the RHS
         if( unit_block == given_unit_block ) {
          // update the coefficient of the commitment variable
          coefficients.push_back( - fixed_consumption );

          // update the Range
          ++range.second;
         }

         // increment due to the commitment variable
         ++active_var_index;
        }
     }  // end( for( g ) )
    }  // end( for( i ) )

    // Finally, we update the RHS of the constraint and the coefficients of
    // the active Variables that belong to the given UnitBlock. Notice that
    // the (abstract) Modifications that will be issued as a result of this
    // update do not concern this UCBlock.

    // update the RHS of the constraint (equality constraint)
    v_node_injection_constraints[ t ][ 0 ].set_both( rhs , eNoBlck );

    // update the coefficients
    static_cast< LinearFunction * >
     ( v_node_injection_constraints[ t ][ 0 ].get_function() )->
     modify_coefficients( std::move( coefficients ) , range , eNoBlck );

   }  // end( for( t ) )
  }
  else {  // number_nodes > 1
   // DCNetwork needs GeneratorNode

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {

    /* The active Variables of the LinearFunction defining the constraint are
     * grouped by UnitBlocks. That is, all active Variables of a given
     * UnitBlock have consecutive indices in this LinearFunction. The
     * following will store the Range of indices of the active Variables of
     * this LinearFunction that belong to the given Block. */
    Range range( Inf< Index >() , Inf< Index >() );

    // This will store the coefficients that must be updated, i.e., those of
    // the active Variables that belong to the given Block.
    LinearFunction::Vec_FunctionValue coefficients;

    Index active_var_index = 0;

    for( Index node_id = 0 ; node_id < number_nodes ; ++node_id ) {

     // increment due to the node injection variable
     ++active_var_index;

     double rhs = 0.0;

     Index elc_generator = 0;
     for( Index unit_id = 0 ; unit_id < f_number_units ; unit_id++ ) {

      const auto unit_block = get_unit_block( unit_id );
      const auto scale = unit_block->get_scale();

      if( unit_block == given_unit_block )
       coefficients.reserve( 2 * unit_block->get_number_generators() );

      for( Index generator = 0 ;
           generator < unit_block->get_number_generators() ;
           ++generator , ++elc_generator ) {

       if( node_id != v_generator_node[ elc_generator ] )
        continue;

       if( unit_block->get_active_power( generator ) ) {

        if( unit_block == given_unit_block ) {
         // update the Range
         if( range.first == Inf< Index >() ) {
          range.first = active_var_index;
          range.second = active_var_index;
         }
         range.second++;
         // update the coefficient of the active power variable
         coefficients.push_back( scale );
        }

        // increment due to the active power variable
        ++active_var_index;
       }

       double fixed_consumption = 0.0;
       if( auto fc = unit_block->get_fixed_consumption( generator ) )
        fixed_consumption = fc[ t ] * scale;

       if( unit_block->get_commitment( generator ) ) {
        if( unit_block == given_unit_block ) {
         // update the Range
         if( range.first == Inf< Index >() ) {
          range.first = active_var_index;
          range.second = active_var_index;
         }
         range.second++;
         // update the coefficient of the commitment variable
         coefficients.push_back( - fixed_consumption );
        }

        // increment due to the commitment variable
        ++active_var_index;
       }

       rhs -= fixed_consumption;
      }
     }

     // Finally, we update the RHS of the constraint and the coefficients of
     // the active Variables that belong to the given UnitBlock. Notice that
     // the (abstract) Modifications that will be issued as a result of this
     // update do not concern this UCBlock.

     // update the RHS of the constraint (equality constraint)
     v_node_injection_constraints[ t ][ node_id ].set_both( rhs , eNoBlck );

     // update the coefficients
     static_cast< LinearFunction * >
      ( v_node_injection_constraints[ t ][ node_id ].get_function() )->
      modify_coefficients( std::move( coefficients ) , range , eNoBlck );

    }  // end( for( node_id ) )
   }  // end( for( t ) )
  }
 }   // end( if( number_nodes > 0 ) )
}  // end( UCBlock::update_node_injection_constraints_scale )

/*--------------------------------------------------------------------------*/

void UCBlock::update_node_injection_constraints( Index time , Index node_index ,
                                                 double demand ) {
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
      rhs -= scale * fc[ time ]; // update the RHS
     }
  }  // end( for( g ) )
 }  // end( for( i ) )

 v_node_injection_constraints[ time ][ node_index ].set_both( rhs );
}

/*--------------------------------------------------------------------------*/

void UCBlock::set_active_power_demand
( std::vector< double >::const_iterator values , Block::Subset && subset ,
  const bool ordered , c_ModParam issuePMod , c_ModParam issueAMod ) {

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
  const auto demand = *( values++ );

  if( v_active_power_demand[ node_index ][ time ] != demand ) {
   changed = true;

   if( not_dry_run( issuePMod ) ) {

    // Change the physical representation
    v_active_power_demand[ node_index ][ time ] = demand;

    if( not_dry_run( issueAMod ) && constraints_generated() ) {
     // Change the abstract representation
     update_node_injection_constraints( time , node_index , demand );
    }
   }
  }
 }

 // If nothing changes, return
 if( ! changed )
  return;

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification

  Block::add_Modification(
   std::make_shared< UCBlockSbstMod >( this , UCBlockMod::eSetActD ,
                                       std::move( subset ) ) ,
   Observer::par2chnl( issuePMod ) );
 }
}

/*--------------------------------------------------------------------------*/

void UCBlock::set_active_power_demand
( std::vector< double >::const_iterator values , Block::Range rng ,
  c_ModParam issuePMod , c_ModParam issueAMod ) {

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
  const auto demand = *( values++ );

  if( v_active_power_demand[ node_index ][ time ] != demand ) {
   changed = true;

   if( not_dry_run( issuePMod ) ) {
    // Change the physical representation
    v_active_power_demand[ node_index ][ time ] = demand;

    if( not_dry_run( issueAMod ) && constraints_generated() ) {
     // Change the abstract representation
     update_node_injection_constraints( time , node_index , demand );
    }
   }
  }
 }

 // If nothing changes, return
 if( ! changed )
  return;

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification

  Block::add_Modification(
   std::make_shared< UCBlockRngdMod >( this , UCBlockMod::eSetActD , rng ) ,
   Observer::par2chnl( issuePMod ) );
 }
}

/*--------------------------------------------------------------------------*/
/*------------------------ End File UCBlock.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
