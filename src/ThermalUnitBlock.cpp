/*--------------------------------------------------------------------------*/
/*--------------------- File ThermalUnitBlock.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the ThermalUnitBlock class.
 *
 * \version 0.11
 *
 * \date 02 - 05 - 2019
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
void ThermalUnitBlock::deserialize( const netCDF::NcGroup & group,
                                    const std::string & var_name,
                                    const std::vector<int> & Change_Interval,
                                    std::vector<T> & data ) {

  auto ncVar = group.getVar( var_name );
  if( ncVar.isNull() )
    throw( std::invalid_argument( "ThermalUnitBlock::deserialize: " +
                                  var_name + " is not present" ) );

  std::vector < size_t > start = { 0 };
  std::vector < size_t > count_number_intervals = { f_number_intervals };

  data.resize( f_time_horizon );

  if( f_number_intervals == f_time_horizon )

    ncVar.getVar( start, count_number_intervals, data.data() );

  else {

    std::vector<T> netCDF_data( f_number_intervals );
    ncVar.getVar( start, count_number_intervals, netCDF_data.data() );

    if( f_number_intervals == 1)

      data.assign( f_time_horizon, netCDF_data[ 0 ] );

    else {

      int i = 0;
      for( int j = 0; j < f_number_intervals; ) {
        data[ i ] = netCDF_data[ j ];
        if( ++i > Change_Interval[ j ] )
          ++j;
      }
    }
  }
}

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::deserialize( netCDF::NcGroup & group ) {

  // check the data- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  netCDF::NcDim TH = group.getDim( "TimeHorizon" );
  if( TH.isNull() )
    throw( std::invalid_argument
           ( "ThermalUnitBlock::deserialize: TimeHorizon not present" ) );

  f_time_horizon = TH.getSize();
  if( f_time_horizon <= 0 )
    throw( std::invalid_argument
           ( "ThermalUnitBlock::deserialize: TimeHorizon must be positive" ) );

  netCDF::NcDim NV = group.getDim( "NumberIntervals" );
  if( NV.isNull() )
    f_number_intervals = 0;
  else {
    f_number_intervals = NV.getSize();
    if( ( f_number_intervals < 1 ) || ( f_number_intervals > f_time_horizon ) )
      throw( std::invalid_argument
             ( "ThermalUnitBlock::deserialize: invalid f_number_intervals" ) );
  }

  std::vector < int > Change_Interval;

  // check problem data- - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( f_number_intervals > 1 ) && ( f_number_intervals < f_time_horizon ) ) {

    Change_Interval.resize( f_number_intervals );
    netCDF::NcVar CI = group.getVar( "ChangeIntervals" );
    CI.getVar( { 0 } , { f_number_intervals } , Change_Interval.data() );

    // Check that all numbers are between 1 and f_time_horizon, that
    // the last number is == f_time_horizon, and that they are ordered
    // in increasing sense

    if( Change_Interval.back() != f_time_horizon ) {
      throw( std::invalid_argument
             ( "ThermalUnitBlock::deserialize: invalid value in "
               "ChangeIntervals: the last element must be TimeHorizon." ) );
    }

    int previous_t = 0;

    for( auto t : Change_Interval ) {
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

  deserialize( group, "MinPower",      Change_Interval, f_MinPower );
  deserialize( group, "MaxPower",      Change_Interval, f_MaxPower );
  deserialize( group, "DeltaRampUp",   Change_Interval, f_DeltaRampUp );
  deserialize( group, "DeltaRampDown", Change_Interval, f_DeltaRampDown );
  deserialize( group, "PrimaryRho",    Change_Interval, f_PrimaryRho );
  deserialize( group, "SecondaryRho",  Change_Interval, f_SecondaryRho );
  deserialize( group, "LinearTerm",    Change_Interval, f_LinearTerm );
  deserialize( group, "QuadTerm",      Change_Interval, f_QuadTerm );
  deserialize( group, "ConstTerm",     Change_Interval, f_ConstTerm );

  ::deserialize( group, "FixedConsPower",  & f_FixedConsPower );
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

    double variable_value = - 1.0;

    if( f_InitUpDownTime < 0 && -f_InitUpDownTime < f_MinDownTime )
        variable_value = 0.0;

    else if( f_InitUpDownTime > 0 && f_InitUpDownTime < f_MinUpTime)
        variable_value = 1.0;

    if ( variable_value >= 0.0 ) {

        //we need to fix v_commitment variable to 0 or 1 for init_t time steps
        for (int t = 0; t < init_t ; ++t) {
            v_commitment[t].set_value(variable_value);
            v_commitment[t].is_fixed(true);
        }
    }

/*----------------Fixing-the-power-variables-to-0---------------------------*/

    double power_value= -1.0;

    if( f_InitUpDownTime < 0 && -f_InitUpDownTime < f_MinDownTime )
        power_value = 0.0;

    else if( f_InitUpDownTime > 0 && f_InitUpDownTime < f_MinUpTime)
        power_value > 0;

    if (  power_value == 0.0 ) {


        for (int t = 0; t < init_t ; ++t) {
            //we just need to fix v_power variable to 0 for init_t time steps
            v_active_power[t].set_value(power_value);
            v_active_power[t].is_fixed(true);

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

void ThermalUnitBlock::serialize( netCDF::NcGroup & group ) const {

  group.putAtt( "type" , "ThermalUnitBlock" );

  netCDF::NcDim time_horizon = group.addDim( "TimeHorizon" , f_time_horizon );
  netCDF::NcDim number_values = group.addDim( "NumberIntervals" , f_number_intervals);

  ::serialize( group, "FixedConsPower", netCDF::NcDouble(), f_FixedConsPower );
  ::serialize( group, "PZero",          netCDF::NcDouble(), f_PZero );
  ::serialize( group, "StartUpCost",    netCDF::NcDouble(), f_StartUpCost );
  ::serialize( group, "MinUpTime",      netCDF::NcUint64(), f_MinUpTime );
  ::serialize( group, "MinDownTime",    netCDF::NcUint64(), f_MinDownTime );
  ::serialize( group, "InitUpDownTime", netCDF::NcUint64(), f_InitUpDownTime );

  // TODO serialize the vectors

  std::vector < size_t > startp = { 0 };
  std::vector < size_t > countpt = { (size_t)f_time_horizon };
  std::vector < size_t > countpu = { f_number_intervals };

/*---------------------MinPower-serialize-----------------------------------*/

    if( f_MinPower.size() )
        if( f_number_intervals == f_time_horizon )
            ( group.addVar( "MinPower" , netCDF::NcDouble() , time_horizon )
            ).putVar( startp , countpt , f_MinPower.data() );
        else {
            ( group.addVar( "MinPower" , netCDF::NcDouble() , number_values )
            ).putVar( startp , countpu , f_MinPower.data() );
        }

/*---------------------MaxPower-serialize-----------------------------------*/

    if( f_MaxPower.size() )
        if( f_number_intervals == f_time_horizon )
            ( group.addVar( "MaxPower" , netCDF::NcDouble() , time_horizon )
            ).putVar( startp , countpt , f_MaxPower.data() );
        else {
            ( group.addVar( "MaxPower" , netCDF::NcDouble() , number_values )
            ).putVar( startp , countpu , f_MaxPower.data() );
        }

/*---------------------DeltaRampUp-serialize--------------------------------*/


    if( f_DeltaRampUp.size() )
        if( f_number_intervals == f_time_horizon )
            ( group.addVar( "DeltaRampUp" , netCDF::NcDouble() , time_horizon )
            ).putVar( startp , countpt , f_DeltaRampUp.data() );
        else {
            ( group.addVar( "DeltaRampUp" , netCDF::NcDouble() , number_values )
            ).putVar( startp , countpu , f_DeltaRampUp.data() );
        }

/*---------------------DeltaRampDown-serialize------------------------------*/

    if( f_DeltaRampDown.size() )
        if( f_number_intervals == f_time_horizon )
            ( group.addVar( "DeltaRampDown" , netCDF::NcDouble() , time_horizon )
            ).putVar( startp , countpt , f_DeltaRampDown.data() );
        else{
            ( group.addVar( "DeltaRampDown" , netCDF::NcDouble() , number_values )
            ).putVar( startp , countpu , f_DeltaRampDown.data() );
        }

/*--------------------------PrimaryRho-serialize----------------------------*/

    if( f_PrimaryRho.size() )
        if( f_number_intervals == f_time_horizon )
            ( group.addVar( "PrimaryRho" , netCDF::NcDouble() , time_horizon )
        ).putVar( startp , countpt , f_PrimaryRho.data() );
        else{
            ( group.addVar( "PrimaryRho" , netCDF::NcDouble() , number_values )
            ).putVar( startp , countpu , f_PrimaryRho.data() );
        }

/*-------------------------SecondaryRho-serialize---------------------------*/

    if( f_SecondaryRho.size() )
        if( f_number_intervals == f_time_horizon )
            ( group.addVar( "SecondaryRho" , netCDF::NcDouble() , time_horizon )
        ).putVar( startp , countpt , f_SecondaryRho.data() );
        else{
            ( group.addVar( "SecondaryRho" , netCDF::NcDouble() , number_values )
            ).putVar( startp , countpu , f_SecondaryRho.data() );
        }
/*----------------Quadratic-Linear-Constant-Term-serialize------------------*/

    if( f_QuadTerm.size() )
        if( f_number_intervals == f_time_horizon )
            ( group.addVar( "QuadTerm" , netCDF::NcDouble() , time_horizon )
         ).putVar( startp , countpt , f_QuadTerm.data() );
        else{
            ( group.addVar( "QuadTerm" , netCDF::NcDouble() , number_values )
            ).putVar( startp , countpu , f_QuadTerm.data() );
        }

    if( f_LinearTerm.size() )
        if( f_number_intervals == f_time_horizon )
            ( group.addVar( "LinearTerm" , netCDF::NcDouble() , time_horizon )
          ).putVar( startp , countpt , f_LinearTerm.data() );
        else{
            ( group.addVar( "LinearTerm" , netCDF::NcDouble() , number_values )
            ).putVar( startp , countpu , f_LinearTerm.data() );
        }

    if( f_ConstTerm.size() )
        if( f_number_intervals == f_time_horizon )
            ( group.addVar( "ConstTerm" , netCDF::NcDouble() , time_horizon )
         ).putVar( startp , countpt , f_ConstTerm.data() );
        else{
            ( group.addVar( "ConstTerm" , netCDF::NcDouble() , number_values )
            ).putVar( startp , countpu , f_ConstTerm.data() );
        }

}  // end( ThermalUnitBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------------- End File ThermalUnitBlock.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
