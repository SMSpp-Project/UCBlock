/*--------------------------------------------------------------------------*/
/*------------------------ File SlackUnitBlock.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the SlackUnitBlock class, which derives from UnitBlock
 * [see UnitBlock.h] and implements a "slack" unit; a (typically, fictitious)
 * unit capable of producing (typically, a large amount of) active power
 * and/or primary/secondary reserve and/or inertia at any time period
 * completely indeopendently from each other and from all other time periods,
 * albeit at a (typically, huge) cost.
 *
 * \version 0.11
 *
 * \date 10 - 12 - 2019
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
 * Copyright &copy by Antonio Frangioni, Ali Ghezelsoflu
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "SlackUnitBlock.h"
#include "DQuadFunction.h"
#include "FRealObjective.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register SlackUnitBlock to the Block factory


SMSpp_insert_in_factory_cpp_1( SlackUnitBlock );

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS OF SlackUnitBlock ------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
void SlackUnitBlock::deserialize( netCDF::NcGroup & group ) {

 UnitBlock::deserialize_time_horizon( group );
 UnitBlock::deserialize_change_intervals( group );

 ::deserialize( group, "MaxPower", f_number_intervals,
                v_MaxPower, true, true );
 ::deserialize( group, "MaxPrimaryPower", f_number_intervals,
                v_MaxPrimaryPower, true, true );
 ::deserialize( group, "MaxSecondaryPower", f_number_intervals,
                v_MaxSecondaryPower, true, true );
 ::deserialize( group, "ActivePowerCost", f_number_intervals,
                v_active_power_cost, true, true );
 ::deserialize( group, "PrimaryCost", f_number_intervals,
                v_primary_cost, true, true );
 ::deserialize( group, "SecondaryCost", f_number_intervals,
                v_secondary_cost, true, true );
 ::deserialize( group, "InertiaCost", f_number_intervals,
                v_inertia_cost, true, true );
 ::deserialize( group, "MaxInertia", v_MaxInertia, true, true );

 UnitBlock::deserialize( group );
}// end( SlackUnitBlock::deserialize )

/*--------------------------------------------------------------------------*/

void SlackUnitBlock::generate_abstract_variables
        ( Configuration *stvv )
{
/*--------------------------------------------------------------------------*/
 if ( f_time_horizon > 0 ){

  // Commitment Variable

  if( v_commitment.size() != f_time_horizon ) {
   assert( v_commitment.empty() ); // this should only happen once
   v_commitment.resize( f_time_horizon );
   int n = 0;
   for( auto & i : v_commitment ) {
    i.set_type( ColVariable::kBinary );
    add_static_variable( i, "u_" + std::to_string( n++ ) );
   }
  }

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

} // end( SlackUnitBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void SlackUnitBlock::generate_abstract_constraints
        ( Configuration *stcc )
{

} // end( SlackUnitBlock::generate_abstract_constraints )


/*--------------------------------------------------------------------------*/

void SlackUnitBlock::generate_objective( Configuration *objc )
{
 // Initial condition of each vector
 std::vector<double> active_power_cost = this->v_active_power_cost;
 if( active_power_cost.size() == 1 ) {
  active_power_cost.resize( f_time_horizon, active_power_cost[ 0 ] );
 }else if (active_power_cost.size() < f_time_horizon ) {
  active_power_cost.resize(f_time_horizon);
  int j = 0;
  for (unsigned long i = 0; i < v_change_intervals.size(); ++i) {
   Index sup;
   if (i == v_change_intervals.size() - 1 ) {
    sup = f_time_horizon;
   } else {
    sup = v_change_intervals[i];
   }
   for (; j < sup; ++j) {
    active_power_cost[j] = v_active_power_cost[i];
   }
  }
 }


 std::vector<double> primary_cost = this->v_primary_cost;
 if( primary_cost.size() == 1 ) {
  primary_cost.resize( f_time_horizon, primary_cost[ 0 ] );
 }else if (primary_cost.size() < f_time_horizon ) {
  primary_cost.resize(f_time_horizon);
  int j = 0;
  for (unsigned long i = 0; i < v_change_intervals.size(); ++i) {
   Index sup;
   if (i == v_change_intervals.size() - 1 ) {
    sup = f_time_horizon;
   } else {
    sup = v_change_intervals[i];
   }
   for (; j < sup; ++j) {
    primary_cost[j] = v_primary_cost[i];
   }
  }
 }

 std::vector<double> secondary_cost = this->v_secondary_cost;
 if( secondary_cost.size() == 1 ) {
  secondary_cost.resize( f_time_horizon, secondary_cost[ 0 ] );
 }else if (secondary_cost.size() < f_time_horizon ) {
  secondary_cost.resize(f_time_horizon);
  int j = 0;
  for (unsigned long i = 0; i < v_change_intervals.size(); ++i) {
   Index sup;
   if (i == v_change_intervals.size() - 1 ) {
    sup = f_time_horizon;
   } else {
    sup = v_change_intervals[i];
   }
   for (; j < sup; ++j) {
    secondary_cost[j] = v_secondary_cost[i];
   }
  }
 }


 std::vector<double> inertia_cost = this->v_inertia_cost;
 if( inertia_cost.size() == 1 ) {
  inertia_cost.resize( f_time_horizon, inertia_cost[ 0 ] );
 }else if (inertia_cost.size() < f_time_horizon ) {
  inertia_cost.resize(f_time_horizon);
  int j = 0;
  for (unsigned long i = 0; i < v_change_intervals.size(); ++i) {
   Index sup;
   if (i == v_change_intervals.size() - 1 ) {
    sup = f_time_horizon;
   } else {
    sup = v_change_intervals[i];
   }
   for (; j < sup; ++j) {
    inertia_cost[j] = v_inertia_cost[i];
   }
  }
 }

 boost::multi_array< double, 2 > MaxInertia = v_MaxInertia;

 if( MaxInertia.size() == 1 ) {
  MaxInertia.resize( boost::extents[f_time_horizon][1] );
  for( Index t = 0; t < f_time_horizon; ++t ) {
   MaxInertia[t][0] = v_MaxInertia[0][0];

  }
 } else if( MaxInertia.size() < f_time_horizon ) {

  MaxInertia.resize( boost::extents[f_time_horizon][1] );

   int j = 0;
   for( unsigned long i = 0; i < v_change_intervals.size(); ++i ) {
    Index sup;
    if( i == v_change_intervals.size() - 1 ) {
     sup = f_time_horizon;
    } else {
     sup = v_change_intervals[i];
    }
    for( ; j < sup; ++j ) {

     MaxInertia[j][0] = v_MaxInertia[i][0];

   }
  }
 }
/*--------------------------------------------------------------------------*/
 if( get_objective() != nullptr )  // an objective is there already
  return;                        // cowardly (and silently) return

 // Initialize objective function

 if( v_commitment.size() != f_time_horizon ) {
  throw ( std::logic_error
          ( "SlackUnitBlock::generate_objective: v_commitment must have "
            "size equal to the time horizon." ) );
 }
 if( v_active_power.size() != f_time_horizon ) {
  throw ( std::logic_error
          ( "SlackUnitBlock::generate_objective: v_active_power must have "
            "size equal to the time horizon." ) );
 }

 if( v_primary_spinning_reserve.size() != f_time_horizon ) {
  throw ( std::logic_error
          ( "SlackUnitBlock::generate_objective: v_primary_spinning_reserve "
            "must have size equal to the time horizon." ) );
 }

 if( v_secondary_spinning_reserve.size() != f_time_horizon ) {
  throw ( std::logic_error
          ( "SlackUnitBlock::generate_objective: v_secondary_spinning_reserve"
            "must have size equal to the time horizon." ) );
 }


 auto dquad_function = new DQuadFunction();

 for( Index t = 0; t < f_time_horizon; ++t ) {
  dquad_function->add_variable( &v_active_power[ t ],
                                active_power_cost[ t ],
                                0.0 );
  dquad_function->add_variable( &v_primary_spinning_reserve[ t ],
                                primary_cost[ t ],
                                0.0 );

  dquad_function->add_variable( &v_secondary_spinning_reserve[ t ],
                                secondary_cost[ t ],
                                0.0 );
  dquad_function->add_variable( &v_commitment[ t ],
                                inertia_cost[ t ] * MaxInertia[ t ][0],
                                0.0 );
 }
 objective.set_function( dquad_function );
 objective.set_sense( Objective::eMin );

 // Set Block objective
 this->set_objective( &objective );
}  // end( SlackUnitBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*------- METHODS FOR LOADING, PRINTING & SAVING THE SlackUnitBlock --------*/
/*--------------------------------------------------------------------------*/

void SlackUnitBlock::serialize( netCDF::NcGroup & group ) const {

 UnitBlock::serialize( group );

 auto NumberIntervals = group.getDim( "NumberIntervals" );

 ::serialize( group, "MaxPower", netCDF::NcDouble(),
              NumberIntervals, v_MaxPower, true );

 ::serialize( group, "MaxPrimaryPower", netCDF::NcDouble(),
              NumberIntervals, v_MaxPrimaryPower, true );

 ::serialize( group, "MaxSecondaryPower", netCDF::NcDouble(),
              NumberIntervals, v_MaxSecondaryPower, true );

 ::serialize( group, "ActivePowerCost", netCDF::NcDouble(),
              NumberIntervals, v_active_power_cost, true );

 ::serialize( group, "PrimaryCost", netCDF::NcDouble(),
              NumberIntervals, v_primary_cost, true );

 ::serialize( group, "SecondaryCost", netCDF::NcDouble(),
              NumberIntervals, v_secondary_cost, true );

 ::serialize( group, "InertiaCost", netCDF::NcDouble(),
              NumberIntervals, v_inertia_cost, true );

 ::serialize( group, "MaxInertia", netCDF::NcDouble(),
              { NumberIntervals }, v_MaxInertia, true );
}  // end( SlackUnitBlock::serialize )

/*--------------------------------------------------------------------------*/
/*----------------------- End File SlackUnitBlock.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
