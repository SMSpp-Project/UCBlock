/*--------------------------------------------------------------------------*/
/*------------- File IntermittentUnitBlock.cpp -------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the IntermittentUnitBlock class.
 *
 * \version 0.11
 *
 * \date 25 - 07 - 2019
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 *
 * \author Ali Ghezelsoflu \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 *
 * Copyright &copy by Antonio Frangioni, Ali Ghezelsoflu
 */

/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <iostream>
#include <random>
#include "IntermittentUnitBlock.h"
#include "LinearFunction.h"
#include <map>
#include "UnitBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;


/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register IntermittentUnitBlock to the Block factory


SMSpp_insert_in_factory_cpp_1( IntermittentUnitBlock );

/*--------------------------------------------------------------------------*/
/*--------------------- METHODS OF IntermittentUnitBlock -------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
void IntermittentUnitBlock::deserialize( netCDF::NcGroup & group ) {


#ifndef NDEBUG
 std::cerr << "[DEBUG] IntermittentUnitBlock::deserialize() - Checking Dims"
           << std::endl;
 std::vector< std::string > expected_dims = { "TimeHorizon",
                                              "NumberIntervals" };
 check_dimensions( group, expected_dims, std::cerr );

 std::cerr << "[DEBUG] IntermittentUnitBlock::deserialize() - Checking Vars"
           << std::endl;
 std::vector< std::string > expected_vars = { "MinPower",
                                              "MaxPower",
                                              "InertiaPower",
                                              "Gamma",
                                              "Kappa" };
 check_variables( group, expected_vars, std::cerr );
#endif



 UnitBlock::deserialize_time_horizon( group );
 UnitBlock::deserialize_change_intervals( group );

  ::deserialize( group, "MinPower", f_time_horizon, v_minimum_power, true, true );

  ::deserialize( group, "MaxPower",f_time_horizon,v_maximum_power, true, true );

 ::deserialize( group, "MaxPower",v_maximum_power, true );

 ::deserialize( group, "InertiaPower", v_inertia_power, true, true );

 ::deserialize( group, "Gamma", &f_gamma );

 ::deserialize( group, "Kappa", &f_kappa );

 decompress_vector( v_minimum_power );
 decompress_vector( v_maximum_power );

 UnitBlock::deserialize( group );
}// end( IntermittentUnitBlock::deserialize )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::generate_abstract_variables
        ( Configuration *stvv )
{
 UnitBlock::generate_abstract_variables( stvv );

 if ( f_time_horizon > 0 ){

  // Active Power Variable

  if( v_active_power.size() != f_time_horizon ) {
   assert( v_active_power.empty() ); // this should only happen once
   v_active_power.resize( f_time_horizon );
   int n = 0;
   for( auto & i : v_active_power ) {
    i.set_type( ColVariable::kNonNegative );
    add_static_variable( i, "p_" + std::to_string( n++ ) );
   }
  }

  // Primary Spinning Reserve Variable

  if( v_primary_spinning_reserve.size() != f_time_horizon ) {
   assert( v_primary_spinning_reserve.empty() ); // this should only happen once
   v_primary_spinning_reserve.resize( f_time_horizon );
   int n = 0;
   for( auto & i : v_primary_spinning_reserve ) {
    i.set_type( ColVariable::kNonNegative );
    add_static_variable( i, "pr_" + std::to_string( n++ ) );
   }
  }

  // Secondary Spinning Reserve Variable

  if( v_secondary_spinning_reserve.size() != f_time_horizon ) {
   assert( v_secondary_spinning_reserve.empty() ); // this should only happen once
   v_secondary_spinning_reserve.resize( f_time_horizon );
   int n = 0;
   for( auto & i : v_secondary_spinning_reserve ) {
    i.set_type( ColVariable::kNonNegative );
    add_static_variable( i, "sr_" + std::to_string( n++ ) );
   }
  }

 }

 AR |= HasVar;
} // end( IntermittentUnitBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::generate_abstract_constraints
        ( Configuration *stcc )
{
 // initial condition of each vector
 std::vector<double> min_power = v_minimum_power;
 if (min_power.size() == 1) {
  min_power.resize(f_time_horizon, f_kappa * min_power[0]);

 } else if (min_power.size() < f_time_horizon) {
  min_power.resize(f_time_horizon);
  int j = 0;
  for (unsigned long i = 0; i < v_change_intervals.size(); ++i) {
   Index sup;
   if (i == v_change_intervals.size() - 1 ) {
    sup = f_time_horizon;
   } else {
    sup = v_change_intervals[i];
   }
   for (; j < sup; ++j) {
    min_power[j] = f_kappa * v_minimum_power[i];
   }
  }
 }

 std::vector<double> max_power = v_maximum_power;
 if( max_power.size() == 1 ) {
  max_power.resize( f_time_horizon, f_kappa * max_power[0] );
 }else if (max_power.size() < f_time_horizon) {
  max_power.resize(f_time_horizon);
  int j = 0;
  for (unsigned long i = 0; i < v_change_intervals.size(); ++i) {
   Index sup;
   if (i == v_change_intervals.size() - 1 ) {
    sup = f_time_horizon;
   } else {
    sup = v_change_intervals[i];
   }
   for (; j < sup; ++j) {
    max_power[j] = f_kappa * v_maximum_power[i];
   }
  }
 }

/*--------------------------------------------------------------------------*/
 // Initializing maximum power constraints

 if( f_gamma != 0 ) {
  // INITIAL CONDITION
  std::vector< double > Max_power_kappa_gamma = v_maximum_power;
  if( Max_power_kappa_gamma.size() == 1 ) {
   Max_power_kappa_gamma.resize( f_time_horizon, f_kappa * f_gamma * Max_power_kappa_gamma[0] );
  } else if( Max_power_kappa_gamma.size() < f_time_horizon ) {
   Max_power_kappa_gamma.resize( f_time_horizon );
   int j = 0;
   for( unsigned long i = 0; i < v_change_intervals.size(); ++i ) {
    Index sup;
    if( i == v_change_intervals.size() - 1 ) {
     sup = f_time_horizon;
    } else {
     sup = v_change_intervals[i];
    }
    for( ; j < sup; ++j ) {
     Max_power_kappa_gamma[j] = f_kappa * f_gamma * v_maximum_power[i];
    }
   }
  }
  if( MaxPower_Constraints.size() != f_time_horizon ) {
   // this should only happen once
   assert( MaxPower_Constraints.empty());

   MaxPower_Constraints.resize( f_time_horizon );
  }


  for( Index t = 0; t < f_time_horizon; ++t ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_active_power[t], f_gamma );
   linear_function->add_variable( &v_primary_spinning_reserve[t], 1.0 );
   linear_function->add_variable( &v_secondary_spinning_reserve[t], 1.0 );

   MaxPower_Constraints[t].set_lhs( -Inf< double >());
   MaxPower_Constraints[t].set_rhs(( Max_power_kappa_gamma[t] ));
   MaxPower_Constraints[t].set_function( linear_function );
  }

  add_static_constraint( MaxPower_Constraints, "MaxPower_c" );


  // Initializing minimum power constraints

  if( MinPower_Constraints.size() != f_time_horizon ) {
   // this should only happen once
   assert( MinPower_Constraints.empty());

   MinPower_Constraints.resize( f_time_horizon );
  }

  for( Index t = 0; t < f_time_horizon; ++t ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_active_power[t], 1.0 );
   linear_function->add_variable( &v_primary_spinning_reserve[t], -1.0 );
   linear_function->add_variable( &v_secondary_spinning_reserve[t], -1.0 );

   MinPower_Constraints[t].set_lhs( min_power[t] );
   MinPower_Constraints[t].set_rhs( Inf< double >());
   MinPower_Constraints[t].set_function( linear_function );
  }

  add_static_constraint( MinPower_Constraints, "MinPower_c" );
 }

 // Initializing active power bounds constraints
 if( active_power_bounds_Constraints.size() != f_time_horizon ) {
  // this should only happen once
  assert( active_power_bounds_Constraints.empty());

  active_power_bounds_Constraints.resize( f_time_horizon );
 }

 for( Index t = 0; t < f_time_horizon; ++t ) {

  auto linear_function = new LinearFunction();

  linear_function->add_variable( &v_active_power[t], 1.0 );

  active_power_bounds_Constraints[t].set_lhs( min_power[ t ]);
  active_power_bounds_Constraints[t].set_rhs( max_power[ t ] );
  active_power_bounds_Constraints[t].set_function( linear_function );
 }

 add_static_constraint( active_power_bounds_Constraints, "ActivePowerBound_c" );

 AR |= HasCst;
} // end( IntermittentUnitBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/
/*--- METHODS FOR LOADING, PRINTING & SAVING THE IntermittentUnitBlock -----*/
/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::serialize( netCDF::NcGroup & group ) const {

 UnitBlock::serialize( group );

 auto NumberIntervals = group.getDim( "NumberIntervals" );

 ::serialize( group, "Gamma", netCDF::NcDouble(), f_gamma );

 ::serialize( group, "Kappa", netCDF::NcDouble(), f_kappa );

 ::serialize( group, "MinPower", netCDF::NcDouble(),
              NumberIntervals, v_minimum_power, true );

 ::serialize( group, "MaxPower", netCDF::NcDouble(),
              NumberIntervals, v_maximum_power, true );

 ::serialize( group, "InertiaPower", netCDF::NcDouble(),
              { NumberIntervals }, v_inertia_power, true );
}  // end( IntermittentUnitBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/
void IntermittentUnitBlock::set_maximum_power(
 std::vector< double >::const_iterator values,
 Block::Subset && subset,
 const bool ordered,
 c_ModParam issuePMod,
 c_ModParam issueAMod ) {
 if( subset.empty() ) {
  return;
 }

 if( v_maximum_power.empty() ) {
  if( std::all_of( values,
                   values + subset.size(),
                   []( double cst ) {
                    return ( cst == 0 );
                   } ) ) {
   return;
  }

  Index max_index = *max_element( std::begin( subset ), std::end( subset ) );
  v_maximum_power.assign( max_index, 0 );
 }

 // If nothing changes, return
 bool identical = true;
 for( auto i : subset ) {
  if( i >= v_maximum_power.size() ) {
   throw ( std::invalid_argument( "invalid value in subset" ) );
  }
  if( v_maximum_power[ i ] != *( values++ ) ) {
   identical = false;
  }
 }
 if( identical ) {
  return;
 }

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  for( auto i : subset ) {
   v_maximum_power[ i ] = *( values++ );
  }

  if( not_dry_run( issueAMod ) && AR & HasCst ) {
   // Change the abstract representation

   if( !MaxPower_Constraints.empty() ) {
    for( auto t : subset ) {
     MaxPower_Constraints[ t ]
      .set_rhs( f_kappa * f_gamma * v_maximum_power[ t ], issueAMod );
     // FIXME: use a GroupModification
    }
   }
   if( !active_power_bounds_Constraints.empty() ) {
    for( auto t : subset ) {
     active_power_bounds_Constraints[ t ]
      .set_rhs( v_maximum_power[ t ], issueAMod );
     // FIXME: use a GroupModification
    }
   }
  }
 }

 if( issue_pmod( issuePMod ) ) {
  // Issue a Physical Modification
  if( !ordered ) {
   std::sort( subset.begin(), subset.end() );
  }

  Block::add_Modification(
   std::make_shared< IntermittentUnitBlockSbstMod >( this,
                                                     IntermittentUnitBlockMod::eSetMaxP,
                                                     std::move( subset ) ),
   Observer::par2chnl( issuePMod ) );
 }
}

void IntermittentUnitBlock::set_maximum_power(
 std::vector< double >::const_iterator values,
 Block::Range rng,
 c_ModParam issuePMod,
 c_ModParam issueAMod ) {
 rng.second = std::min( rng.second, f_number_intervals );
 if( rng.second <= rng.first ) {
  return;
 }

 if( v_maximum_power.empty() ) {
  if( std::all_of( values,
                   values + ( rng.second - rng.first ),
                   []( double cst ) {
                    return ( cst == 0 );
                   } ) ) {
   return;
  }

  Index max_index = rng.second;
  v_maximum_power.assign( max_index, 0 );
 }

 // If nothing changes, return
 if( std::equal( values,
                 values + ( rng.second - rng.first ),
                 v_maximum_power.begin() + rng.first ) ) {
  return;
 }

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation

  std::copy( values,
             values + ( rng.second - rng.first ),
             v_maximum_power.begin() + rng.first );

  if( not_dry_run( issueAMod ) && AR & HasCst ) {
   // Change the abstract representation

   if( !MaxPower_Constraints.empty() ) {
    for( Index t = rng.first; t < rng.second; ++t ) {
     MaxPower_Constraints[ t ]
      .set_rhs( f_kappa * f_gamma * v_maximum_power[ t ], issueAMod );
     // FIXME: use a GroupModification
    }
   }
   if( !active_power_bounds_Constraints.empty() ) {
    for( Index t = rng.first; t < rng.second; ++t ) {
     active_power_bounds_Constraints[ t ]
      .set_rhs( v_maximum_power[ t ], issueAMod );
     // FIXME: use a GroupModification
    }
   }
  }
 }

 if( issue_pmod( issuePMod ) ) {
  Block::add_Modification(
   std::make_shared< IntermittentUnitBlockRngdMod >( this,
                                                     IntermittentUnitBlockMod::eSetMaxP,
                                                     rng ),
   Observer::par2chnl( issuePMod ) );
 }
}


template< typename T >
void IntermittentUnitBlock::decompress_vector( std::vector< T > & v ) {
 if( v.size() == 1 ) {
  v.resize( f_time_horizon, v[ 0 ] );
 } else if( v.size() < f_time_horizon ) {
  std::vector< T > temp = v;
  v.resize( f_time_horizon );
  int j = 0;
  for( unsigned long i = 0; i < v_change_intervals.size(); ++i ) {
   Index sup;
   if( i == v_change_intervals.size() - 1 ) {
    sup = f_time_horizon;
   } else {
    sup = v_change_intervals[ i ];
   }
   for( ; j < sup; ++j ) {
    v[ j ] = temp[ i ];
   }
  }
 }
}

/*--------------------------------------------------------------------------*/
/*------------------- End File IntermittentUnitBlock.cpp -------------------*/
/*--------------------------------------------------------------------------*/
