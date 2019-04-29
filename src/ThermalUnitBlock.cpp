/*--------------------------------------------------------------------------*/
/*--------------------- File ThermalUnitBlock.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the ThermalUnitBlock class.
 *
 * \version 0.11
 *
 * \date 17 - 04 - 2019
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

/*@} -----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */



//void ThermalUnitBlock::load( )
//{

//generate_abstract_variables( );

//}


/*--------------------------------------------------------------------------*/



/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::deserialize( netCDF::NcGroup & group , Block * father )
{

// check the data- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 netCDF::NcDim TH = group.getDim( "TimeHorizon" );
 if( TH.isNull() )
  throw( std::invalid_argument( "TimeHorizon not present" ) );

 f_time_horizon = TH.getSize();
 if( f_time_horizon <= 0 )
  throw( std::invalid_argument( "TimeHorizon <= 0" ) );
  
 netCDF::NcDim NV = group.getDim( "NumberValues" );
 if( NV.isNull() )
  f_number_values = 0;
 else {
  f_number_values = NV.getSize();
  if( ( f_number_values < 1 ) || ( f_number_values > f_time_horizon ) )
   throw( std::invalid_argument( "invalid f_number_values" ) );
  }
// read problem data- - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 std::vector < size_t > start = { 0 };
 std::vector < size_t > countime = { (size_t)f_time_horizon };
 std::vector < size_t > countch = { f_number_values };
 std::vector < int > Change_Interval;

 if( ( f_number_values > 1 ) && ( f_number_values < f_time_horizon ) ) {
  Change_Interval.resize( f_number_values );
  netCDF::NcVar CI = group.getVar( "ChangeIntervals" );
  CI.getVar( start , countch , Change_Interval.data() );

  // TODO: check that all numbers are between 1 and f_time_horizon,
  // that the last number is == f_time_horizon, and that they are
  // ordered in increasing sense
  }
 
/*--------------------MinPower-deserialize----------------------------------*/

 netCDF::NcVar MinP = group.getVar( "MinPower" );
 if( MinP.isNull() )
  throw( std::invalid_argument( "Min Power not present" ) );

 f_MinPower.resize( f_time_horizon );
    
 if( f_number_values == f_time_horizon )
  MinP.getVar( start , countime , f_MinPower.data() );
 else {
  std::vector<double> tmpv( f_number_values );
  MinP.getVar( start , countch , f_MinPower.data() );

  if( f_number_values == 1 )
   f_MinPower.assign( f_time_horizon , tmpv[ 0 ] );
  else {
   int i = 0;
   for( int j = 0 ; j < f_number_values ; ) {
    f_MinPower[ i ] = tmpv[ j ];
    if( ++i > Change_Interval[ j ] )
     ++j;
    }
   }
  }

/*--------------------MaxPower-deserialize----------------------------------*/

    netCDF::NcVar MaxP = group.getVar( "MaxPower" );
    if( MaxP.isNull() )
        throw( std::invalid_argument( "Max Power not present" ) );

    f_MaxPower.resize( f_time_horizon );

    if( f_number_values == f_time_horizon )
        MaxP.getVar( start , countime , f_MaxPower.data() );
    else {
        std::vector<double> tMpv( f_number_values );
        MaxP.getVar( start , countch , f_MaxPower.data() );

        if( f_number_values == 1 )
            f_MaxPower.assign( f_time_horizon , tMpv[ 0 ] );
        else {
            int i = 0;
            for( int j = 0 ; j < f_number_values ; ) {
                f_MaxPower[ i ] = tMpv[ j ];
                if( ++i > Change_Interval[ j ] )
                    ++j;
            }
        }
    }

/*------------------FixedConsPower-deserialize------------------------------*/

 netCDF::NcVar FixCoPow = group.getVar( "FixedConsPower" );
 if( FixCoPow.isNull() )
  throw( std::invalid_argument( "FixedConsPower not present" ) );


/*--------------------DeltaRampUp-deserialize-------------------------------*/

    netCDF::NcVar RampUp = group.getVar( "DeltaRampUp" );
    if( RampUp.isNull() )
        throw( std::invalid_argument( "DeltaRampUp not present" ) );

    f_DeltaRampUp.resize( f_time_horizon );

    if( f_number_values == f_time_horizon )
        RampUp.getVar( start , countime , f_DeltaRampUp.data() );
    else {
        std::vector<double> trup( f_number_values );
        RampUp.getVar( start , countch , f_DeltaRampUp.data() );

        if( f_number_values == 1 )
            f_DeltaRampUp.assign( f_time_horizon , trup[ 0 ] );
        else {
            int i = 0;
            for( int j = 0 ; j < f_number_values ; ) {
                f_DeltaRampUp[ i ] = trup[ j ];
                if( ++i > Change_Interval[ j ] )
                    ++j;
            }
        }
    }

/*-------------------DeltaRampDown-deserialize------------------------------*/

    netCDF::NcVar RampDown = group.getVar( "DeltaRampDown" );
    if( RampDown.isNull() )
        throw( std::invalid_argument( "DeltaRampDown not present" ) );

    f_DeltaRampDown.resize( f_time_horizon );

    if( f_number_values == f_time_horizon )
        RampDown.getVar( start , countime , f_DeltaRampDown.data() );
    else {
        std::vector<double> trdn( f_number_values );
        RampDown.getVar( start , countch , f_DeltaRampDown.data() );

        if( f_number_values == 1 )
            f_DeltaRampDown.assign( f_time_horizon , trdn[ 0 ] );
        else {
            int i = 0;
            for( int j = 0 ; j < f_number_values ; ) {
                f_DeltaRampDown[ i ] = trdn[ j ];
                if( ++i > Change_Interval[ j ] )
                    ++j;
            }
        }
    }
/*-------------------------PrimaryRho-deserialize---------------------------*/

    netCDF::NcVar PrimRho = group.getVar( "PrimaryRho" );
    if( PrimRho.isNull() )
        throw( std::invalid_argument( "PrimaryRho not present" ) );

    f_PrimaryRho.resize( f_time_horizon );

    if( f_number_values == f_time_horizon )
        PrimRho.getVar( start , countime , f_PrimaryRho.data() );
    else {
        std::vector<double> tprho( f_number_values );
        PrimRho.getVar( start , countch , f_PrimaryRho.data() );

        if( f_number_values == 1 )
            f_PrimaryRho.assign( f_time_horizon , tprho[ 0 ] );
        else {
            int i = 0;
            for( int j = 0 ; j < f_number_values ; ) {
                f_PrimaryRho[ i ] = tprho[ j ];
                if( ++i > Change_Interval[ j ] )
                    ++j;
            }
        }
    }
/*
    netCDF::NcVar PrimRho = group.getVar("PrimaryRho");
    if ( ! PrimRho.isNull() ) {
        f_PrimaryRho.resize(f_time_horizon);
        PrimRho.getVar(start, countime, f_PrimaryRho.data());

    }
    */
/*------------------------SecondaryRho-deserialize--------------------------*/

    netCDF::NcVar SecondRho = group.getVar( "SecondaryRho" );
    if( SecondRho.isNull() )
        throw( std::invalid_argument( "SecondaryRho not present" ) );

    f_SecondaryRho.resize( f_time_horizon );

    if( f_number_values == f_time_horizon )
        SecondRho.getVar( start , countime , f_SecondaryRho.data() );
    else {
        std::vector<double> tsrho( f_number_values );
        SecondRho.getVar( start , countch , f_SecondaryRho.data() );

        if( f_number_values == 1 )
            f_SecondaryRho.assign( f_time_horizon , tsrho[ 0 ] );
        else {
            int i = 0;
            for( int j = 0 ; j < f_number_values ; ) {
                f_SecondaryRho[ i ] = tsrho[ j ];
                if( ++i > Change_Interval[ j ] )
                    ++j;
            }
        }
    }

/*---------------Quadratic-Linear-Constant-Term-deserialize-----------------*/

    netCDF::NcVar QTerm = group.getVar( "QuadTerm" );
    if( QTerm.isNull() )
        throw( std::invalid_argument( "QuadTerm not present" ) );

    f_QuadTerm.resize( f_time_horizon );

    if( f_number_values == f_time_horizon )
        QTerm.getVar( start , countime , f_QuadTerm.data() );
    else {
        std::vector<double> tqt( f_number_values );
        QTerm.getVar( start , countch , f_QuadTerm.data() );

        if( f_number_values == 1 )
            f_QuadTerm.assign( f_time_horizon , tqt[ 0 ] );
        else {
            int i = 0;
            for( int j = 0 ; j < f_number_values ; ) {
                f_QuadTerm[ i ] = tqt[ j ];
                if( ++i > Change_Interval[ j ] )
                    ++j;
            }
        }
    }
//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    netCDF::NcVar LTerm = group.getVar( "LinearTerm" );
    if( LTerm.isNull() )
        throw( std::invalid_argument( "LinearTerm not present" ) );

    f_LinearTerm.resize( f_time_horizon );

    if( f_number_values == f_time_horizon )
        LTerm.getVar( start , countime , f_LinearTerm.data() );
    else {
        std::vector<double> tlt( f_number_values );
        LTerm.getVar( start , countch , f_LinearTerm.data() );

        if( f_number_values == 1 )
            f_LinearTerm.assign( f_time_horizon , tlt[ 0 ] );
        else {
            int i = 0;
            for( int j = 0 ; j < f_number_values ; ) {
                f_LinearTerm[ i ] = tlt[ j ];
                if( ++i > Change_Interval[ j ] )
                    ++j;
            }
        }
    }

//- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    netCDF::NcVar CTerm = group.getVar( "ConstTerm" );
    if( CTerm.isNull() )
        throw( std::invalid_argument( "ConstTerm not present" ) );

    f_ConstTerm.resize( f_time_horizon );

    if( f_number_values == f_time_horizon )
        CTerm.getVar( start , countime , f_ConstTerm.data() );
    else {
        std::vector<double> tct( f_number_values );
        CTerm.getVar( start , countch , f_ConstTerm.data() );

        if( f_number_values == 1 )
            f_ConstTerm.assign( f_time_horizon , tct[ 0 ] );
        else {
            int i = 0;
            for( int j = 0 ; j < f_number_values ; ) {
                f_ConstTerm[ i ] = tct[ j ];
                if( ++i > Change_Interval[ j ] )
                    ++j;
            }
        }
    }

/*----------------------------PZero-deserialize-----------------------------*/

netCDF::NcVar PZero  = group.getVar( "PZero" );
    if( PZero.isNull() )
        throw( std::invalid_argument( "PZero not present" ) );

/*-------------------------StartUpCost-deserialize--------------------------*/

netCDF::NcVar StartUpC  = group.getVar( "StartUpCost" );
    if( StartUpC.isNull() )
        throw( std::invalid_argument( "StartUpC not present" ) );

/*-----------------MinUp--MinDown--InitUpDown-deserialize-------------------*/

netCDF::NcVar MinUp  = group.getVar( "MinUpTime" );
    if( MinUp.isNull() )
        throw( std::invalid_argument( "MinUpTime not present" ) );

netCDF::NcVar MinDown  = group.getVar( "MinDownTime" );
    if( MinDown.isNull() )
        throw( std::invalid_argument( "MinDownTime not present" ) );

netCDF::NcVar IntUpDown  = group.getVar( "IntUpDownTime" );
    if( IntUpDown.isNull() )
        throw( std::invalid_argument( "IntUpDownTime not present" ) );

}  // end( ThermalUnitBlock::deserialize )

/*--------------------------------------------------------------------------*/


void ThermalUnitBlock::generate_abstract_variables( Configuration *stvv )
{


    if ( InitUpDownTime_val > 0 )
    {

        init_t =  ( InitUpDownTime_val >= MinUpTime_val ? 0 :
                    MinUpTime_val   - InitUpDownTime_val );
    }
    else
    {
        init_t = ( -InitUpDownTime_val >= MinDownTime_val ? 0 :
                   MinDownTime_val + InitUpDownTime_val );

    }

/*--------Define-Binary-variables-start_up_thermal-&-shut_down_thermal------*/

  if( f_time_horizon < 0 ) {
    throw( std::logic_error( "ThermalUnitBlock::generate_abstract_variables: "
                             "time horizon of ThermalUnitBlock is not set" ) );
  }


    if( v_start_up_thermal.size() && v_shut_down_thermal.size() !=
                                                    f_time_horizon - init_t) {
    assert( v_start_up_thermal.size() == 0 ); // this should only happen once
        v_start_up_thermal.resize( f_time_horizon - init_t);

    assert( v_shut_down_thermal.size() == 0 ); // this should only happen once
        v_shut_down_thermal.resize( f_time_horizon - init_t);


      for(int t = init_t ; t < f_time_horizon ; t++){


            v_start_up_thermal[t].set_type(ColVariable::kBinary);
            v_shut_down_thermal[t].set_type(ColVariable::kBinary);



      }
        add_static_variable ( v_start_up_thermal );
        add_static_variable ( v_shut_down_thermal );
    }
} // end( ThermalUnitBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::generate_abstract_constraints( Configuration *stcc )
{


/*------------Fixing-the-commitment-variables-to-0-or-1---------------------*/

    double variable_value = - 1.0;

    if( InitUpDownTime_val < 0 && -InitUpDownTime_val < MinDownTime_val )
        variable_value = 0.0;

    else if( InitUpDownTime_val > 0 && InitUpDownTime_val < MinUpTime_val)
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

    if( InitUpDownTime_val < 0 && -InitUpDownTime_val < MinDownTime_val )
        power_value = 0.0;

    else if( InitUpDownTime_val > 0 && InitUpDownTime_val < MinUpTime_val)
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


if( InitUpDownTime_val > 0 && InitUpDownTime_val < MinUpTime_val)

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

        linear_function->add_variable( & v_commitment[t],               1.0 );
        linear_function->add_variable( & v_commitment[t - 1],          -1.0 );
        linear_function->add_variable( & v_start_up_thermal[t]  ,      -1.0 );
        linear_function->add_variable( & v_shut_down_thermal[t] ,       1.0 );

        UVW_Const[t].set_both( 0.0 );
        UVW_Const[t].set_function( linear_function );
    }
//initializing Turn ON Inequalities (UV Constraints)
UV_Const.resize( f_time_horizon - init_t );

    for(int t  = MinUpTime_val + init_t ; t < f_time_horizon  ; t++){

            auto linear_function = new LinearFunction();

            for (int s =  t - MinUpTime_val +1 ;  s < t ; s++){

                linear_function->add_variable( & v_start_up_thermal[s], 1.0 );
            }

            linear_function->add_variable( & v_commitment[t],          -1.0 );

            UV_Const[t].set_rhs( 0.0 );
            UV_Const[t].set_lhs( - Inf<double>() );
            UV_Const[t].set_function( linear_function );

    }


//initializing Turn OFF Inequalities (UW Constraints)
UW_Const.resize(f_time_horizon - init_t);

    for(int t = MinDownTime_val - init_t ; t < f_time_horizon  ; t++){


            auto linear_function = new LinearFunction();

            for (int s =  t - MinDownTime_val + 1 ;  s < t; s++){

               linear_function->add_variable( & v_shut_down_thermal[s], 1.0 );
            }

            linear_function->add_variable( & v_commitment[t],           1.0 );

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
        linear_function->add_variable( & v_start_up_thermal[t+1],
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
        linear_function->add_variable( & v_shut_down_thermal[t+1],
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

void ThermalUnitBlock::generate_objective( Configuration *objc )
{
    if( ! get_objective().empty() )  // an objective is there already
        return;                         // cowardly (and silently) return

// initialize objective function - - - - - - - - - - - - - - - - - - - - - -
 //TODO
}  // end( ThermalUnitBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*-------- METHODS FOR LOADING, PRINTING & SAVING THE ThermalUnitBlock -----*/
/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::serialize( netCDF::NcGroup & group ) const
{

group.putAtt( "type" , "ThermalUnitBlock" );

netCDF::NcDim time_horizon = group.addDim( "TimeHorizon" , f_time_horizon );
netCDF::NcDim number_values = group.addDim( "NumberValues" , f_number_values);


std::vector < size_t > startp = { 0 };
std::vector < size_t > countpt = { (size_t)f_time_horizon };
std::vector < size_t > countpu = { f_number_values };


/*---------------------MinPower-serialize-----------------------------------*/
    if( f_MinPower.size() )
        if( f_number_values == f_time_horizon )
            ( group.addVar( "MinPower" , netCDF::NcDouble() , time_horizon )
            ).putVar( startp , countpt , f_MinPower.data() );
        else {
            ( group.addVar( "MinPower" , netCDF::NcDouble() , number_values )
            ).putVar( startp , countpu , f_MinPower.data() );
        }


/*---------------------MinPower-serialize-----------------------------------*/


    if( f_MaxPower.size() )
        if( f_number_values == f_time_horizon )
            ( group.addVar( "MaxPower" , netCDF::NcDouble() , time_horizon )
            ).putVar( startp , countpt , f_MaxPower.data() );
        else {
            ( group.addVar( "MaxPower" , netCDF::NcDouble() , number_values )
            ).putVar( startp , countpu , f_MaxPower.data() );
        }


/*-------------------FixedConsPower-serialize-------------------------------*/

netCDF::NcVar FixCoPow  = group.addVar( "FixedConsPower" ,
                                                        netCDF::NcUint64() );

/*---------------------DeltaRampUp-serialize--------------------------------*/


    if( f_DeltaRampUp.size() )
        if( f_number_values == f_time_horizon )
            ( group.addVar( "DeltaRampUp" , netCDF::NcDouble() , time_horizon )
            ).putVar( startp , countpt , f_DeltaRampUp.data() );
        else {
            ( group.addVar( "DeltaRampUp" , netCDF::NcDouble() , number_values )
            ).putVar( startp , countpu , f_DeltaRampUp.data() );
        }

/*---------------------DeltaRampDown-serialize------------------------------*/

    if( f_DeltaRampDown.size() )
        if( f_number_values == f_time_horizon )
            ( group.addVar( "DeltaRampDown" , netCDF::NcDouble() , time_horizon )
            ).putVar( startp , countpt , f_DeltaRampDown.data() );
        else{
            ( group.addVar( "DeltaRampDown" , netCDF::NcDouble() , number_values )
            ).putVar( startp , countpu , f_DeltaRampDown.data() );
        }

/*--------------------------PrimaryRho-serialize----------------------------*/

    if( f_PrimaryRho.size() )
        if( f_number_values == f_time_horizon )
            ( group.addVar( "PrimaryRho" , netCDF::NcDouble() , time_horizon )
        ).putVar( startp , countpt , f_PrimaryRho.data() );
        else{
            ( group.addVar( "PrimaryRho" , netCDF::NcDouble() , number_values )
            ).putVar( startp , countpu , f_PrimaryRho.data() );
        }

/*-------------------------SecondaryRho-serialize---------------------------*/

    if( f_SecondaryRho.size() )
        if( f_number_values == f_time_horizon )
            ( group.addVar( "SecondaryRho" , netCDF::NcDouble() , time_horizon )
        ).putVar( startp , countpt , f_SecondaryRho.data() );
        else{
            ( group.addVar( "SecondaryRho" , netCDF::NcDouble() , number_values )
            ).putVar( startp , countpu , f_SecondaryRho.data() );
        }
/*-----------------------------PZero-serialize------------------------------*/

netCDF::NcVar PZero  = group.addVar( "PZero" ,       netCDF::NcUint64() );

/*----------------Quadratic-Linear-Constant-Term-serialize------------------*/

    if( f_QuadTerm.size() )
        if( f_number_values == f_time_horizon )
            ( group.addVar( "QuadTerm" , netCDF::NcDouble() , time_horizon )
         ).putVar( startp , countpt , f_QuadTerm.data() );
        else{
            ( group.addVar( "QuadTerm" , netCDF::NcDouble() , number_values )
            ).putVar( startp , countpu , f_QuadTerm.data() );
        }

    if( f_LinearTerm.size() )
        if( f_number_values == f_time_horizon )
            ( group.addVar( "LinearTerm" , netCDF::NcDouble() , time_horizon )
          ).putVar( startp , countpt , f_LinearTerm.data() );
        else{
            ( group.addVar( "LinearTerm" , netCDF::NcDouble() , number_values )
            ).putVar( startp , countpu , f_LinearTerm.data() );
        }

    if( f_ConstTerm.size() )
        if( f_number_values == f_time_horizon )
            ( group.addVar( "ConstTerm" , netCDF::NcDouble() , time_horizon )
         ).putVar( startp , countpt , f_ConstTerm.data() );
        else{
            ( group.addVar( "ConstTerm" , netCDF::NcDouble() , number_values )
            ).putVar( startp , countpu , f_ConstTerm.data() );
        }
/*--------------------------StartUpCost-serialize---------------------------*/

netCDF::NcVar StartUpC  = group.addVar( "StartUpCost" ,  netCDF::NcUint64() );

/*------------------MinUp--MinDown--InitUpDown-serialize--------------------*/


netCDF::NcVar MinUp  = group.addVar( "MinUpTime" ,       netCDF::NcUint64() );

netCDF::NcVar MinDown  = group.addVar( "MinDownTime" ,   netCDF::NcUint64() );

netCDF::NcVar InitUpDown  = group.addVar( "InitUpDownTime" ,
                                                         netCDF::NcUint64() );

}  // end( ThermalUnitBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------------- End File ThermalUnitBlock.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
