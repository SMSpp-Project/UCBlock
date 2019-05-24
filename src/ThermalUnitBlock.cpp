/*--------------------------------------------------------------------------*/
/*--------------------- File ThermalUnitBlock.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the ThermalUnitBlock class.
 *
 * \version 0.11
 *
 * \date 24 - 05 - 2019
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
 * Copyright &copy by Antonio Frangioni, Ali Ghezelsoflu, and Rafael
 * Durbano Lobato
 */

/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <iostream>
#include <random>
#include "ThermalUnitBlock.h"
#include "LinearFunction.h"
#include <map>
#include "FRowConstraint.h"
#include "DQuadFunction.h"
#include "FRealObjective.h"
#include "UnitBlock.h"
#include "UCBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;


/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register ThermalUnitBlock to the Block factory


SMSpp_insert_in_factory_cpp_1( ThermalUnitBlock );

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS OF ThermalUnitBlock ----------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

//void ThermalUnitBlock::load( )
//{ }

/*--------------------------------------------------------------------------*/

template<class T>
void deserialize( const netCDF::NcGroup & group,
                  const std::string & var_name,
                  T * data ) {

  auto ncVar = group.getVar( var_name );
  if( ncVar.isNull() )
    throw( std::invalid_argument( "ThermalUnitBlock::deserialize: " +
                                  var_name + " is not present" ) );
  ncVar.getVar( data );
}

/*--------------------------------------------------------------------------*/

template<class T>
void deserialize( const netCDF::NcGroup & group, const std::string & var_name,
                  const size_t & size, std::vector<T> & data ) {

  auto ncVar = group.getVar( var_name );
  if( ncVar.isNull() )
    throw( std::invalid_argument( "ThermalUnitBlock::deserialize: " +
                                  var_name + " is not present" ) );

  data.resize( size );
  ncVar.getVar( { 0 }, { size }, data.data() );
}

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::deserialize( netCDF::NcGroup & group ) {

  // check the data- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  auto TH = group.getDim( "TimeHorizon" );
  if( TH.isNull() )
    throw( std::invalid_argument
           ( "ThermalUnitBlock::deserialize: TimeHorizon not present" ) );

  f_time_horizon = TH.getSize();
  if( f_time_horizon <= 0 )
    throw( std::invalid_argument
           ( "ThermalUnitBlock::deserialize: TimeHorizon must be positive" ) );

  auto NV = group.getDim( "NumberIntervals" );
  if( NV.isNull() )
    f_number_intervals = 0;
  else {
    f_number_intervals = NV.getSize();
    if( ( f_number_intervals < 1 ) || ( f_number_intervals > f_time_horizon ) )
      throw( std::invalid_argument
             ( "ThermalUnitBlock::deserialize: invalid f_number_intervals" ) );
  }

  // check problem data- - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( f_number_intervals > 1 ) && ( f_number_intervals < f_time_horizon ) ) {

    ::deserialize( group, "ChangeIntervals", f_number_intervals,
                   f_change_interval );

    // Check that all numbers are between 1 and f_time_horizon, that
    // the last number is == f_time_horizon, and that they are ordered
    // in increasing sense

    if( f_change_interval.back() != f_time_horizon ) {
      throw( std::invalid_argument
             ( "ThermalUnitBlock::deserialize: invalid value in "
               "ChangeIntervals: the last element must be TimeHorizon." ) );
    }

    int previous_t = 0;

    for( auto t : f_change_interval ) {
      if( ! ( t > previous_t && t < f_time_horizon - 1 ) )
        throw( std::invalid_argument
               ( "ThermalUnitBlock::deserialize: invalid value in "
                 "ChangeIntervals: " + std::to_string( t ) + ". All values "
                 "must be between 1 and TimeHorizon and in strictly "
                 "increasing order." ) );

      previous_t = t;
    }
  }

  // read problem data- - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ::deserialize( group, "MinPower",      f_number_intervals, f_MinPower );
  ::deserialize( group, "MaxPower",      f_number_intervals, f_MaxPower );
  ::deserialize( group, "DeltaRampUp",   f_number_intervals, f_DeltaRampUp );
  ::deserialize( group, "DeltaRampDown", f_number_intervals, f_DeltaRampDown );
  ::deserialize( group, "PrimaryRho",    f_number_intervals, f_PrimaryRho );
  ::deserialize( group, "SecondaryRho",  f_number_intervals, f_SecondaryRho );
  ::deserialize( group, "LinearTerm",    f_number_intervals, f_LinearTerm );
  ::deserialize( group, "QuadTerm",      f_number_intervals, f_QuadTerm );
  ::deserialize( group, "ConstTerm",     f_number_intervals, f_ConstTerm );

 // ::deserialize( group, "FixedConsPower",  & f_FixedConsPower );
  ::deserialize( group, "PZero",           & f_PZero );
  ::deserialize( group, "StartUpCost",     & f_StartUpCost );
  ::deserialize( group, "MinUpTime",       & f_MinUpTime );
  ::deserialize( group, "MinDownTime",     & f_MinDownTime );
  ::deserialize( group, "IntUpDownTime",   & f_InitUpDownTime );

}  // end( ThermalUnitBlock::deserialize )

/*--------------------------------------------------------------------------*/


void ThermalUnitBlock::generate_abstract_variables( Configuration *stvv ) {


    if ( f_InitUpDownTime > 0 )
      {

        init_t =  ( f_InitUpDownTime >= f_MinUpTime ? 0 :
                    f_MinUpTime   - f_InitUpDownTime );
    }
    else
    {
        init_t = ( -f_InitUpDownTime >= f_MinDownTime ? 0 :
                   f_MinDownTime + f_InitUpDownTime );

    }

/*--------Define-Binary-variables-start_up-&-shut_down------*/

  if( f_time_horizon < 0 ) {
    throw( std::logic_error( "ThermalUnitBlock::generate_abstract_variables: "
                             "time horizon of ThermalUnitBlock is not set" ) );
  }


    if( v_start_up.size() && v_shut_down.size() !=
                                                    f_time_horizon - init_t) {
    assert( v_start_up.size() == 0 ); // this should only happen once
        v_start_up.resize( f_time_horizon - init_t);

    assert( v_shut_down.size() == 0 ); // this should only happen once
        v_shut_down.resize( f_time_horizon - init_t);


      for(int t = init_t ; t < f_time_horizon ; t++){


            v_start_up[t].set_type(ColVariable::kBinary);
            v_shut_down[t].set_type(ColVariable::kBinary);



      }
        add_static_variable ( v_start_up );
        add_static_variable ( v_shut_down );
    }
} // end( ThermalUnitBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::generate_abstract_constraints( Configuration *stcc )
{


/*------------Fixing-the-commitment-variables-to-0-or-1---------------------*/

    double commitment_variable_value = - 1.0;

    if( f_InitUpDownTime < 0 && -f_InitUpDownTime < f_MinDownTime ) {

        // unit must remain off from time 0 to init_t - 1

        commitment_variable_value = 0.0;

        for (int t = 0; t < init_t; ++t) {
            v_active_power[t].set_value( 0.0 );
            v_active_power[t].is_fixed( true );
        }
    }

    else if( f_InitUpDownTime > 0 && f_InitUpDownTime < f_MinUpTime)
        // unit must remain on from time 0 to init_t - 1
        commitment_variable_value = 1.0;

    if ( commitment_variable_value >= 0.0 ) {

        //we need to fix v_commitment variable to 0 or 1 for init_t time steps
        for (int t = 0; t < init_t ; ++t) {
            v_commitment[t].set_value( commitment_variable_value );
            v_commitment[t].is_fixed( true );
        }
    }

/*-----------------Power-Output-Constraints---------------------------------*/

//PowerOutPut Constraints when v_commitment fix to 1 for init_t time steps

PowerFix_Const.resize(init_t);


if( f_InitUpDownTime > 0 && f_InitUpDownTime < f_MinUpTime)

    for (int t = 0; t < init_t ; ++t){

        PowerFix_Const[t].set_variable( & v_active_power[ t ] , 1.0 );
        PowerFix_Const[t].set_rhs (f_MaxPower[t]);
        PowerFix_Const[t].set_lhs (f_MinPower[t]);
        PowerFix_Const[t].set_Block( this );

    }
    add_static_constraint(PowerFix_Const);

//initializing PMin Inequalities with commitment variables

    PMin_Const.resize(f_time_horizon - init_t );

    for (int t= init_t  ; t < f_time_horizon ; t++) {

        auto linear_function = new LinearFunction();

        linear_function->add_variable( & v_active_power[t],             1.0 );

        linear_function->add_variable( & v_primary_spinning_reserve[t],
                                                                       -1.0 );

        linear_function->add_variable( & v_secondary_spinning_reserve[t],
                                                                       -1.0 );

        linear_function->add_variable( & v_commitment[t],     -f_MinPower[t]);

        PMin_Const[t].set_rhs( Inf<double>() );
        PMin_Const[t].set_lhs( 0.0 );
        PMin_Const[t].set_function(linear_function);

    }
    add_static_constraint(PMin_Const);


//initializing PMax Inequalities with commitment variables

    PMax_Const.resize(f_time_horizon - init_t );

    for (int t = init_t ; t < f_time_horizon ; t++) {

        auto linear_function = new LinearFunction();

        linear_function->add_variable( & v_active_power[t],             1.0 );
        linear_function->add_variable( & v_primary_spinning_reserve[t], 1.0 );

        linear_function->add_variable( & v_secondary_spinning_reserve[t],
                                                                        1.0 );

        linear_function->add_variable( & v_commitment[t],     -f_MaxPower[t]);

        PMax_Const[t].set_rhs( 0.0 );
        PMax_Const[t].set_lhs(-Inf<double>());
        PMax_Const[t].set_function(linear_function);

    }
    add_static_constraint(PMax_Const);

/*
//initializing power injected Equalities
    PowerInjected_Const.resize(f_time_horizon);
    for (int t = 0 ; t < f_time_horizon ; t++) {

        auto linear_function = new LinearFunction();
        linear_function->add_variable( & v_active_power[t],            -1.0 );
        linear_function->add_variable( & v_commitment[t],
                                                        f_FixedConsPower );
        linear_function->add_variable( & v_power_injected[t],           1.0 );

        PowerInjected_Const[t].set_both(-f_FixedConsPower);
        PowerInjected_Const[t].set_function(linear_function);

    }
    add_static_constraint(PowerInjected_Const);
*/
//initializing PrimaryRho fraction Inequalities
PrimaryRho_Const.resize(f_time_horizon);
    for (int t = 0 ; t < f_time_horizon ; t++){
        auto linear_function = new LinearFunction();
        linear_function->add_variable( & v_active_power[t],f_PrimaryRho[t] );
        linear_function->add_variable( & v_primary_spinning_reserve[t],
                                                                     - 1.0 );
        PrimaryRho_Const[t].set_rhs(Inf<double>());
        PrimaryRho_Const[t].set_lhs(0.0);
        PrimaryRho_Const[t].set_function(linear_function);

    }
    add_static_constraint(PrimaryRho_Const);

//initializing SecondaryRho fraction Inequalities
    SecondaryRho_Const.resize(f_time_horizon);
    for (int t = 0 ; t < f_time_horizon ; t++){
        auto linear_function = new LinearFunction();
        linear_function->add_variable( & v_active_power[t],
                                                         f_SecondaryRho[t] );
        linear_function->add_variable( & v_secondary_spinning_reserve[t],
                                                                     - 1.0 );
        SecondaryRho_Const[t].set_rhs(Inf<double>());
        SecondaryRho_Const[t].set_lhs(0.0);
        SecondaryRho_Const[t].set_function(linear_function);

    }
    add_static_constraint(SecondaryRho_Const);

//initializing PMax Inequalities with u, v, and w variables

//TODO

//initializing PMax Inequalities with u, v, and w variables

//TODO
/*-----------------Min-Up-Down-Constraints----------------------------------*/

//initializing UVW connection Constraints

UVW_Const.resize(f_time_horizon - init_t);

    for(int t = init_t + 1 ; t < f_time_horizon  ; t++) {

        auto linear_function = new LinearFunction();

        linear_function->add_variable( & v_commitment[t],      1.0 );
        linear_function->add_variable( & v_commitment[t - 1], -1.0 );
        linear_function->add_variable( & v_start_up[t],       -1.0 );
        linear_function->add_variable( & v_shut_down[t],       1.0 );

        UVW_Const[t].set_both( 0.0 );
        UVW_Const[t].set_function( linear_function );
    }
//initializing Turn ON Inequalities (UV Constraints)
UV_Const.resize( f_time_horizon - init_t );

    for(int t  = f_MinUpTime + init_t ; t < f_time_horizon  ; t++){

            auto linear_function = new LinearFunction();

            for (int s =  t - f_MinUpTime +1 ;  s < t ; s++){

                linear_function->add_variable( & v_start_up[s], 1.0 );
            }

            linear_function->add_variable( & v_commitment[t],          -1.0 );

            UV_Const[t].set_rhs( 0.0 );
            UV_Const[t].set_lhs( - Inf<double>() );
            UV_Const[t].set_function( linear_function );

    }


//initializing Turn OFF Inequalities (UW Constraints)
UW_Const.resize(f_time_horizon - init_t);

    for(int t = f_MinDownTime - init_t ; t < f_time_horizon  ; t++){


            auto linear_function = new LinearFunction();

            for (int s =  t - f_MinDownTime + 1 ;  s < t; s++){

               linear_function->add_variable( & v_shut_down[s], 1.0 );
            }

            linear_function->add_variable( & v_commitment[t], 1.0 );

            UW_Const[t].set_rhs( 1.0 );
            UW_Const[t].set_lhs( 0.0 );
            UW_Const[t].set_function( linear_function );
    }

    add_static_constraint(UVW_Const);
    add_static_constraint(UV_Const);
    add_static_constraint(UW_Const);

/*-----------------Ramp-Up-Down-Constraints---------------------------------*/

//initializing RampUp time Inequalities

RampUp_Const.resize(f_time_horizon - init_t);

    for(int t = init_t ; t < f_time_horizon - 1 ; t++){
        auto linear_function = new LinearFunction();
        linear_function->add_variable( & v_active_power[t+1], -1.0 );
        linear_function->add_variable( & v_active_power[t], +1.0 );
        linear_function->add_variable( & v_start_up[t+1],
                 - f_DeltaRampUp[t] );
        linear_function->add_variable( & v_commitment[t+1],
                                  (f_MinPower[t] + f_DeltaRampUp[t]));
        linear_function->add_variable( & v_commitment[t], -f_MinPower[t]);

        RampUp_Const[t].set_rhs( Inf<double>() );
        RampUp_Const[t].set_lhs( 0.0 );
        RampUp_Const[t].set_function(linear_function);

    }

    add_static_constraint(RampUp_Const);

//initializing RampDown time Inequalities

RampDown_Const.resize(f_time_horizon - init_t);

    for(int t = init_t ; t < f_time_horizon - 1 ; t++){
        auto linear_function = new LinearFunction();
        linear_function->add_variable( & v_active_power[t+1], 1.0 );
        linear_function->add_variable( & v_active_power[t], -1.0 );
        linear_function->add_variable( & v_shut_down[t+1],
                      - f_DeltaRampDown[t] );
        linear_function->add_variable( & v_commitment[t],
                                       (f_MinPower[t] +  f_DeltaRampDown[t]));
        linear_function->add_variable( & v_commitment[t+1],   -f_MinPower[t]);

        RampDown_Const[t].set_rhs( Inf<double>() );
        RampDown_Const[t].set_lhs( 0.0 );
        RampDown_Const[t].set_function(linear_function);

    }

    add_static_constraint(RampDown_Const);

/*-----------------Start-Up-Costs-Constraints-------------------------------*/

  //TODO


} // end( ThermalUnitBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::generate_objective( Configuration *objc ) {

  if( ! get_objective().empty() )  // an objective is there already
    return;                        // cowardly (and silently) return

  // initialize objective function - - - - - - - - - - - - - - - - - - - - - -

  if( v_commitment.size() != f_time_horizon ) {
    throw( std::logic_error
           ( "ThermalUnitBlock::generate_objective: v_commitment must have "
             "size equal to the time horizon." ) );
  }

  auto dquad_function = new DQuadFunction();

  for( int t = 0; t < f_time_horizon - init_t; ++t ) {
    dquad_function->add_variable( & v_start_up[t], f_StartUpCost, 0.0 );
  }

  for( int t = 0; t < f_time_horizon; ++t ) {
    dquad_function->add_variable( & v_active_power[t], f_LinearTerm[t],
                                  f_QuadTerm[t] );
    dquad_function->add_variable( & v_commitment[t], f_ConstTerm[t], 0.0 );
  }

  objective.set_function( dquad_function );
  objective.set_sense( Objective::eMin );
  objective.set_Block( this );

}  // end( ThermalUnitBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*-------- METHODS FOR LOADING, PRINTING & SAVING THE ThermalUnitBlock -----*/
/*--------------------------------------------------------------------------*/

template<class T>
void serialize( netCDF::NcGroup & group, const std::string & var_name,
                const netCDF::NcType ncType, T data ) {
  ( group.addVar( var_name , ncType ) ).putVar( & data );
}

/*--------------------------------------------------------------------------*/

template<class T>
void serialize( netCDF::NcGroup & group, const std::string & var_name,
                const netCDF::NcType & ncType, const netCDF::NcDim & ncDim,
                const std::vector<T> & data ) {

  group.addVar( var_name , ncType , ncDim )
    .putVar( { 0 } , { data.size() } , data.data() );
}

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::serialize( netCDF::NcGroup & group ) const {

  group.putAtt( "type" , "ThermalUnitBlock" );

  group.addDim( "TimeHorizon" , f_time_horizon );
  netCDF::NcDim ncdim_number_intervals = group.addDim( "NumberIntervals",
                                                       f_number_intervals );

  //::serialize( group, "FixedConsPower", netCDF::NcDouble(), f_FixedConsPower );
  ::serialize( group, "PZero",          netCDF::NcDouble(), f_PZero );
  ::serialize( group, "StartUpCost",    netCDF::NcDouble(), f_StartUpCost );
  ::serialize( group, "MinUpTime",      netCDF::NcUint64(), f_MinUpTime );
  ::serialize( group, "MinDownTime",    netCDF::NcUint64(), f_MinDownTime );
  ::serialize( group, "InitUpDownTime", netCDF::NcUint64(), f_InitUpDownTime );

  ::serialize( group, "ChangeInterval", netCDF::NcUint64(),
               ncdim_number_intervals, f_change_interval );

  ::serialize( group, "MinPower", netCDF::NcDouble(),
               ncdim_number_intervals, f_MinPower );

  ::serialize( group, "MaxPower", netCDF::NcDouble(),
               ncdim_number_intervals, f_MaxPower );

  ::serialize( group, "DeltaRampUp", netCDF::NcDouble(),
               ncdim_number_intervals, f_DeltaRampUp );

  ::serialize( group, "DeltaRampDown", netCDF::NcDouble(),
               ncdim_number_intervals, f_DeltaRampDown );

  ::serialize( group, "PrimaryRho", netCDF::NcDouble(),
               ncdim_number_intervals, f_PrimaryRho );

  ::serialize( group, "SecondaryRho", netCDF::NcDouble(),
               ncdim_number_intervals, f_SecondaryRho );

  ::serialize( group, "QuadTerm", netCDF::NcDouble(),
               ncdim_number_intervals, f_QuadTerm );

  ::serialize( group, "LinearTerm", netCDF::NcDouble(),
               ncdim_number_intervals, f_LinearTerm );

  ::serialize( group, "ConstTerm", netCDF::NcDouble(),
               ncdim_number_intervals, f_ConstTerm );

}  // end( ThermalUnitBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------------- End File ThermalUnitBlock.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
