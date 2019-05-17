/*--------------------------------------------------------------------------*/
/*------------------------- File HeatBlock.cpp -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the HeatBlock class.
 *
 * \version 0.11
 *
 * \date 17 - 05 - 2019
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
#include "HeatBlock.h"
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

// register HeatBlock to the Block factory


SMSpp_insert_in_factory_cpp_1( HeatBlock );

/*--------------------------------------------------------------------------*/
/*--------------------------- METHODS OF HeatBlock -------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

//void HeatBlock::load( )
//{ }

/*--------------------------------------------------------------------------*/

template<class T>
void deserialize( const netCDF::NcGroup & group,
                  const std::string & var_name,
                  T * data ) {

    auto ncVar = group.getVar( var_name );
    if( ncVar.isNull() )
        throw( std::invalid_argument( "HeatBlock::deserialize: " +
                                      var_name + " is not present" ) );
    ncVar.getVar( data );
}

/*--------------------------------------------------------------------------*/

template<class T>
void deserialize( const netCDF::NcGroup & group, const std::string & var_name,
                  const size_t & size, std::vector<T> & data ) {

    auto ncVar = group.getVar( var_name );
    if( ncVar.isNull() )
        throw( std::invalid_argument( "HeatBlock::deserialize: " +
                                      var_name + " is not present" ) );

    data.resize( size );
    ncVar.getVar( { 0 }, { size }, data.data() );
}

/*--------------------------------------------------------------------------*/

void HeatBlock::deserialize( netCDF::NcGroup & group ) {

    // check the data- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    auto TH = group.getDim( "TimeHorizon" );
    if( TH.isNull() )
        throw( std::invalid_argument
                ( "HeatBlock::deserialize: TimeHorizon not present" ) );

    f_time_horizon = TH.getSize();
    if( f_time_horizon <= 0 )
        throw( std::invalid_argument
                ( "HeatBlock::deserialize: TimeHorizon must be positive" ) );


}  // end( HeatBlock::deserialize )

/*--------------------------------------------------------------------------*/

int HeatBlock::get_variables_to_be_generated( Configuration *stvv ) {

    if( ! stvv )
        return 0;

    // informs which variables must be generated
    int variables_to_be_generated = 0;

    auto tstvv = dynamic_cast<SimpleConfiguration<int> *>( stvv );

    if( ( ! tstvv ) && f_BlockConfig &&
        f_BlockConfig->f_static_variables_Configuration ) {

        tstvv = dynamic_cast<SimpleConfiguration<int> *>
        ( f_BlockConfig->f_static_constraints_Configuration );
    }

    if( tstvv )
        variables_to_be_generated = tstvv->f_value;

    return variables_to_be_generated;
}

/*--------------------------------------------------------------------------*/

void HeatBlock::generate_abstract_variables( Configuration *stvv ) {

    if( v_heat.size() != 0 ) {
        // the abstract variable should be generated only once
        return;
    }

    if( f_time_horizon == 0 ) {
        // there are no variable to be generated
        return;
    }

    typedef std::vector< std::pair< std::vector<ColVariable> * , int > > v_pairs;

    v_pairs variables_and_types = {
            std::make_pair( &v_heat,                       ColVariable::kNonNegative )
            // v_heat must be the last one in this list
    };

    auto variables_to_be_generated = get_variables_to_be_generated( stvv );

    // The heat variables must be always present
    variables_to_be_generated |=
            (int) std::pow( 2, variables_and_types.size() - 1 );

    int k = 1;
    for( auto [ variables, variable_type ] : variables_and_types ) {
        if( variables_to_be_generated & k ) {
            variables->resize( f_time_horizon );
            for( auto & variable : * variables )
                variable.set_type( variable_type );
            add_static_variable( * variables );
        }
        k *= 2;
    }
} // end( HeatBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void HeatBlock::generate_abstract_constraints( Configuration *stcc ) {



} // end( HeatBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void HeatBlock::generate_objective( Configuration *objc ) {



}  // end( HeatBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*---------- METHODS FOR LOADING, PRINTING & SAVING THE HeatBlock ----------*/
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

void HeatBlock::serialize( netCDF::NcGroup & group ) const {

    group.putAtt( "type" , "HeatBlock" );

    group.addDim( "TimeHorizon" , f_time_horizon );


}  // end( HeatBlock::serialize )

/*--------------------------------------------------------------------------*/
/*----------------------- End File HeatBlock.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
