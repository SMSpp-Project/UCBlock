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

 UnitBlock::deserialize( group );

 ::deserialize( group, "MinPower", f_number_intervals, v_minimum_power, true, true );

 ::deserialize( group, "MaxPower", f_number_intervals, v_maximum_power, true, true );

 ::deserialize( group, "InertiaPower", v_inertia_power, true, true );

 ::deserialize( group, "Gamma", &f_gamma );

 ::deserialize( group, "Kappa", &f_kappa );


 if (v_inertia_power.empty()) {
  v_inertia_power.resize( boost::extents[ f_time_horizon ][ 1 ] );
 }

}// end( IntermittentUnitBlock::deserialize )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::generate_abstract_variables
        ( Configuration *stvv )
{
 UnitBlock::generate_abstract_variables( stvv );

} // end( IntermittentUnitBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void IntermittentUnitBlock::generate_abstract_constraints
        ( Configuration *stcc )
{

double kappa = f_kappa ? : 1;
double gamma = f_gamma ? : 0;

 // initial condition of each vector
 std::vector<double> min_power = v_minimum_power;
 if (min_power.size() == 1) {
  min_power.resize(f_time_horizon, min_power[0]);

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
    min_power[j] = v_minimum_power[i];
   }
  }
 }

 std::vector<double> max_power = v_maximum_power;
 if( max_power.size() == 1 ) {
  max_power.resize( f_time_horizon, max_power[0] );
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
    max_power[j] = v_maximum_power[i];
   }
  }
 }

/*--------------------------------------------------------------------------*/

 // Initializing maximum power constraints
  MaxPower_Constraints.resize( f_time_horizon );

  for( Index t = 0; t < f_time_horizon; ++t ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_active_power[t][0], f_gamma );
   linear_function->add_variable( &v_primary_spinning_reserve[t][0], 1.0 );
   linear_function->add_variable( &v_secondary_spinning_reserve[t][0], 1.0 );

   MaxPower_Constraints[t].set_lhs( -Inf< double >());
   MaxPower_Constraints[t].set_rhs( ( kappa * gamma * max_power[ t ]) );
   MaxPower_Constraints[t].set_function( linear_function );
  }

 add_static_constraint( MaxPower_Constraints, "MaxPower_c" );


 // Initializing minimum power constraints

 MinPower_Constraints.resize( f_time_horizon );

  for( Index t = 0; t < f_time_horizon; ++t ) {

   auto linear_function = new LinearFunction();

   linear_function->add_variable( &v_active_power[t][0], 1.0 );
   linear_function->add_variable( &v_primary_spinning_reserve[t][0], -1.0 );
   linear_function->add_variable( &v_secondary_spinning_reserve[t][0], -1.0 );

   MinPower_Constraints[t].set_lhs( kappa * min_power[ t ]);
   MinPower_Constraints[t].set_rhs( Inf< double >() );
   MinPower_Constraints[t].set_function( linear_function );
  }

 add_static_constraint( MinPower_Constraints, "MinPower_c" );


 // Initializing active power bounds constraints

 active_power_bounds_Constraints.resize( f_time_horizon );

 for( Index t = 0; t < f_time_horizon; ++t ) {

  auto linear_function = new LinearFunction();

  linear_function->add_variable( &v_active_power[t][0], 1.0 );

  active_power_bounds_Constraints[t].set_lhs( kappa *min_power[ t ]);
  active_power_bounds_Constraints[t].set_rhs( kappa *max_power[ t ] );
  active_power_bounds_Constraints[t].set_function( linear_function );
 }

 add_static_constraint( active_power_bounds_Constraints, "ActivePowerBound_c" );

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
              {NumberIntervals}, v_inertia_power, true );
}  // end( IntermittentUnitBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------------- End File IntermittentUnitBlock.cpp -------------------*/
/*--------------------------------------------------------------------------*/