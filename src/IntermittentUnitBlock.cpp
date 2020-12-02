/*--------------------------------------------------------------------------*/
/*------------- File IntermittentUnitBlock.cpp -------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the IntermittentUnitBlock class.
 *
 * \version 0.11
 *
 * \date 30 - 09 - 2020
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

IntermittentUnitBlock::~IntermittentUnitBlock() {
 for( auto & constraint : MinPower_Constraints )
  constraint.clear();
 for( auto & constraint : MaxPower_Constraints )
  constraint.clear();
 for( auto & constraint : active_power_bounds_Constraints )
  constraint.clear();
}

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::deserialize( netCDF::NcGroup & group ) {


#ifndef NDEBUG
 std::vector< std::string > expected_dims = { "TimeHorizon",
                                              "NumberIntervals" };
 check_dimensions( group, expected_dims, std::cerr );
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

 ::deserialize( group, "InertiaPower", v_inertia_power, true );

 ::deserialize( group, "Gamma", &f_gamma, false );

if (! ::deserialize( group, "Kappa", &f_kappa, true )) {
 f_kappa = 1;
}

 decompress_vector( v_minimum_power );
 decompress_vector( v_maximum_power );

 UnitBlock::deserialize( group );
}// end( IntermittentUnitBlock::deserialize )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::generate_abstract_variables
        ( Configuration *stvv )
{

 if( variables_generated() )
  return; // variables have already been generated

 UnitBlock::generate_abstract_variables( stvv );

 if( f_time_horizon > 0 ) {

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

 set_variables_generated();
} // end( IntermittentUnitBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::generate_abstract_constraints
        ( Configuration *stcc )
{

 if( constraints_generated() )
  return; // constraints have already been generated

 std::vector<double> max_power = v_maximum_power;
 std::vector<double> min_power = v_minimum_power;

/*--------------------------------------------------------------------------*/
 // Initializing maximum power constraints


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
   MaxPower_Constraints[t].set_rhs(( f_gamma * f_kappa * (max_power[t]) ));
   MaxPower_Constraints[t].set_function( linear_function );
  }

  add_static_constraint( MaxPower_Constraints, "MaxPower_Intermittent" );


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

   MinPower_Constraints[t].set_lhs( f_kappa * min_power[t] );
   MinPower_Constraints[t].set_rhs( Inf< double >());
   MinPower_Constraints[t].set_function( linear_function );
  }

  add_static_constraint( MinPower_Constraints, "MinPower_Intermittent" );

 // Initializing active power bounds constraints
 if( active_power_bounds_Constraints.size() != f_time_horizon ) {
  // this should only happen once
  assert( active_power_bounds_Constraints.empty());

  active_power_bounds_Constraints.resize( f_time_horizon );
 }

 for( Index t = 0; t < f_time_horizon; ++t ) {

  active_power_bounds_Constraints[t].set_lhs( f_kappa * min_power[ t ]);
  active_power_bounds_Constraints[t].set_rhs( f_kappa * max_power[ t ] );
  active_power_bounds_Constraints[t].set_variable( &v_active_power[t] );
 }

 add_static_constraint( active_power_bounds_Constraints, "ActivePowerBound_Intermittent" );

 set_constraints_generated();
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

  if( not_dry_run( issueAMod ) && constraints_generated() ) {
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

  if( not_dry_run( issueAMod ) && constraints_generated() ) {
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
