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


// read problem data- - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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

 std::vector < size_t > start = { 0 };
 std::vector < size_t > countime = { (size_t)f_time_horizon };
 std::vector < size_t > countch = { f_number_values };
 std::vector < ptrdiff_t > stridech = { f_number_values };

 std::vector<int> ci;
 if( ( f_number_values > 1 ) && ( f_number_values < f_time_horizon ) ) {
  ci.resize( f_number_values );
  netCDF::NcVar ChangeIntervals = group.getVar( "ChangeIntervals" );
  ci.getVar( start , countch , ci.data() );

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
    if( ++i > ci[ j ] )
     ++j;
    }
   }
  }

/*--------------------MaxPower-deserialize----------------------------------*/

 netCDF::NcVar MaxP = group.getVar("MaxPower");
 if ( ! MaxP.isNull() ) {
  f_MaxPower.resize(f_time_horizon);
  MaxP.getVar(start, countime, f_MaxPower.data());
 }

/*------------------FixedConsPower-deserialize------------------------------*/

 netCDF::NcVar FixCoPow = group.getVar( "FixedConsPower" );
 if( FixCoPow.isNull() )
  throw( std::invalid_argument( "FixedConsPower not present" ) );
 f_fix_cons_pow = ... ;

  
/*--------------------DeltaRampUp-deserialize-------------------------------*/

 netCDF::NcVar RampUp = group.getVar( "DeltaRampUp" );
 if ( ! RampUp.isNull() ) {
  f_DeltaRampUp.resize(f_time_horizon);
  RampUp.getVar(start, countime, f_DeltaRampUp.data());
  }

/*-------------------DeltaRampDown-deserialize------------------------------*/

    netCDF::NcVar RampDown = group.getVar("DeltaRampDown");
    if ( ! RampDown.isNull() ) {
            f_DeltaRampDown.resize(f_time_horizon);
            RampDown.getVar(start, countime, f_DeltaRampDown.data());
        }

/*-------------------------PrimaryRho-deserialize---------------------------*/

    netCDF::NcVar PrimRho = group.getVar("PrimaryRho");
    if ( ! PrimRho.isNull() ) {
        f_PrimaryRho.resize(f_time_horizon);
        PrimRho.getVar(start, countime, f_PrimaryRho.data());

    }
/*------------------------SecondaryRho-deserialize--------------------------*/

    netCDF::NcVar SecondRho = group.getVar("SecondaryRho");
    if ( ! SecondRho.isNull() ) {
        f_SecondaryRho.resize(f_time_horizon);
        SecondRho.getVar(start, countime, f_SecondaryRho.data());

    }

/*---------------Quadratic-Linear-Constant-Term-deserialize-----------------*/

    netCDF::NcVar QTerm = group.getVar("QuadTerm");
    if ( ! QTerm.isNull() ) {
        f_QuadTerm.resize(f_time_horizon);
        QTerm.getVar(start, countime, f_QuadTerm.data());

    }
        netCDF::NcVar LTerm = group.getVar("LinearTerm");
        if ( ! LTerm.isNull() ) {
            f_LinearTerm.resize(f_time_horizon);
            LTerm.getVar(start, countime, f_LinearTerm.data());

        }

            netCDF::NcVar CTerm = group.getVar("ConstTerm");
            if ( ! CTerm.isNull() ) {
                f_ConstTerm.resize(f_time_horizon);
                CTerm.getVar(start, countime, f_ConstTerm.data());
            }
/*----------------------------PZero-deserialize-----------------------------*/

netCDF::NcVar PZero  = group.getVar( "PZero" );

/*-------------------------StartUpCost-deserialize--------------------------*/

netCDF::NcVar StartUpC  = group.getVar( "StartUpCost" );

/*-----------------MinUp--MinDown--InitUpDown-deserialize-------------------*/

netCDF::NcVar MinUp  = group.getVar( "MinUpTime" );

netCDF::NcVar MinDown  = group.getVar( "MinDownTime" );

netCDF::NcVar IntUpDown  = group.getVar( "IntUpDownTime" );


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

    if (  power_value = 0.0 ) {


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
/*------------ METHODS FOR LOADING, PRINTING & SAVING THE MCFBlock ---------*/
/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::serialize( netCDF::NcGroup & group ) const
{

    group.putAtt( "type" , "ThermalUnitBlock" );
    group.putAtt( "infinity" , netCDF::NcDouble() , Inf<double>() );


netCDF::NcDim time_horizon = group.addDim( "TimeHorizon" , f_time_horizon );
netCDF::NcDim number_values = group.addDim( "NumberValues" , f_number_values);


std::vector < size_t > startp = { 0 };
std::vector < size_t > countpt = { (size_t)f_time_horizon };
std::vector < size_t > countpu = { f_number_values };
std::vector < ptrdiff_t > strideche = { f_number_values };


/*
    if( f_MinPower.size() )
        if ( f_time_horizon == f_number_values) {
            ( group.addVar( "MinPower" , netCDF::NcDouble() , time_horizon )
            ).putVar( startp , countpt , f_MinPower.data() );
        }
        else if ( f_number_values > 0 && f_number_values < f_time_horizon ) {
            ( group.addVar( "MinPower" , netCDF::NcDouble() , time_changes )
         ).putVar( startp , countpu , strideche , f_change_interval.data() );
        }
        else if ( f_number_values == 0 ) {
            netCDF::NcVar MinPow  = group.addVar( "MinPower" ,
                                                  netCDF::NcDouble() );
        }
        else {
            // There is something wrong
        }

*/
/*---------------------MinPower-serialize-----------------------------------*/

    if( f_MinPower.size() )
            ( group.addVar( "MinPower" , netCDF::NcDouble() , time_horizon )
            ).putVar( startp , countpt , f_MinPower.data() );


/*---------------------MinPower-serialize-----------------------------------*/


    if( f_MaxPower.size() )
        ( group.addVar( "MaxPower" , netCDF::NcDouble() , time_horizon )
            ).putVar( startp , countpt , f_MaxPower.data() );

/*-------------------FixedConsPower-serialize-------------------------------*/

netCDF::NcVar FixCoPow  = group.addVar( "FixedConsPower" ,
                                                        netCDF::NcUint64() );

/*---------------------DeltaRampUp-serialize--------------------------------*/


    if( f_DeltaRampUp.size() )
        ( group.addVar( "DeltaRampUp" , netCDF::NcDouble() , time_horizon )
            ).putVar( startp , countpt , f_DeltaRampUp.data() );

/*---------------------DeltaRampDown-serialize------------------------------*/

    if( f_DeltaRampDown.size() )
       ( group.addVar( "DeltaRampDown" , netCDF::NcDouble() , time_horizon )
            ).putVar( startp , countpt , f_DeltaRampDown.data() );

/*--------------------------PrimaryRho-serialize----------------------------*/

    if( f_PrimaryRho.size() )
        ( group.addVar( "PrimaryRho" , netCDF::NcDouble() , time_horizon )
        ).putVar( startp , countpt , f_PrimaryRho.data() );


/*-------------------------SecondaryRho-serialize---------------------------*/

    if( f_SecondaryRho.size() )
        ( group.addVar( "SecondaryRho" , netCDF::NcDouble() , time_horizon )
        ).putVar( startp , countpt , f_SecondaryRho.data() );

/*-----------------------------PZero-serialize------------------------------*/

netCDF::NcVar PZero  = group.addVar( "PZero" ,       netCDF::NcUint64() );

/*----------------Quadratic-Linear-Constant-Term-serialize------------------*/

    if( f_QuadTerm.size() )
        ( group.addVar( "QuadTerm" , netCDF::NcDouble() , time_horizon )
         ).putVar( startp , countpt , f_QuadTerm.data() );

    if( f_LinearTerm.size() )
        ( group.addVar( "LinearTerm" , netCDF::NcDouble() , time_horizon )
          ).putVar( startp , countpt , f_LinearTerm.data() );

    if( f_ConstTerm.size() )
        ( group.addVar( "ConstTerm" , netCDF::NcDouble() , time_horizon )
         ).putVar( startp , countpt , f_ConstTerm.data() );

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
