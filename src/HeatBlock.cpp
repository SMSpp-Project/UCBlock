/*--------------------------------------------------------------------------*/
/*------------------------- File HeatBlock.cpp -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the HeatBlock class.
 *
 * \version 0.11
 *
 * \date 28 - 05 - 2019
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
/*
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
*/

template<class T>
inline void deserialize( const netCDF::NcGroup & group,
                         const std::string & var_name,
                         std::vector<T> & data,
                         const std::vector<size_t> & sizes,
                         bool optional = true ) {

    auto total_size = std::accumulate( begin( sizes ), end( sizes ), 1,
                                       std::multiplies<size_t>() );

    if( total_size == 0 ) {
        data.resize( 0 );
        return;
    }

    auto ncVar = group.getVar( var_name );
    if( ncVar.isNull() ) {
        if( optional ) {
            data.resize( 0 );
            return;
        }
        throw( std::invalid_argument( "HeatBlock::deserialize: " +
                                      var_name + " is not present" ) );
    }

    data.resize( total_size );

    std::vector<size_t> start;
    start.assign( sizes.size(), 0 );

    ncVar.getVar( start, sizes, data.data() );
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


    auto NHU = group.getDim( "NumberHeatUnits" );
    if( NHU.isNull() )
        throw( std::invalid_argument
                ( "HeatBlock::deserialize: NumberHeatUnits not present" ) );

    f_number_heat_units = NHU.getSize();
    if( f_number_heat_units <= 0 )
        throw( std::invalid_argument
                ( "HeatBlock::deserialize: NumberHeatUnits must be "
                                                             "positive" ) );


    auto NV = group.getDim( "NumberIntervals" );
    if( NV.isNull() )
        f_number_intervals = 0;
    else {
        f_number_intervals = NV.getSize();
        if( ( f_number_intervals < 1 ) || ( f_number_intervals >
                                                        f_time_horizon ) )
            throw( std::invalid_argument
               ( "HeatBlock::deserialize: invalid f_number_intervals" ) );
    }

// check problem data- - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if( ( f_number_intervals > 1 ) && ( f_number_intervals < f_time_horizon ))
                                                                           {
        ::deserialize( group, "ChangeIntervals", v_change_interval,
                                                     {f_number_intervals} );

        // Check that all numbers are between 1 and f_time_horizon, that
        // the last number is == f_time_horizon, and that they are ordered
        // in increasing sense

        if( v_change_interval.back() != f_time_horizon ) {
            throw( std::invalid_argument
                    ( "HeatBlock::deserialize: invalid value in "
                      "ChangeIntervals: the last element must be "
                      "TimeHorizon." ) );
        }

        int previous_t = 0;

        for( auto t : v_change_interval ) {
            if( ! ( t > previous_t && t < f_time_horizon - 1 ) )
                throw( std::invalid_argument
                        ( "HeatBlock::deserialize: invalid value in "
                          "ChangeIntervals: " + std::to_string( t ) + ". All "
                          "values must be between 1 and TimeHorizon and in "
                          "strictly increasing order." ) );

            previous_t = t;
        }
    }

// read problem data- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ::deserialize( group, "TotalHeatDemand", v_heat_demand,
                                                     {f_number_intervals} );
        ::deserialize( group, "CostHeatUnit", v_cost_heat_unit,
                                {f_number_intervals , f_number_heat_units});

        ::deserialize( group, "MinHeatProduction", v_min_heat_production,
                                {f_number_intervals , f_number_heat_units});

        ::deserialize( group, "MaxHeatProduction", v_max_heat_production,
                                {f_number_intervals , f_number_heat_units});

        ::deserialize( group, "MinHeatStorage", v_min_heat_storage,
                                                     {f_number_intervals} );

        ::deserialize( group, "MaxHeatStorage", v_max_heat_storage,
                                                     {f_number_intervals} );

        ::deserialize( group, "StoringHeatRho",      & f_storing_heat_rho );

        ::deserialize( group, "ExtractingHeatRho", & f_extracting_heat_rho);

        ::deserialize( group, "KeepingHeatRho",       & f_keeping_heat_rho);


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

/*
    if( v_heat.size() != 0 ) {
        // the abstract variable should be generated only once
        return;
    }

    if( f_time_horizon < 0 ) {
        throw( std::logic_error( "HeatBlock::generate_abstract_variables: "
                                 "time horizon of HeatBlock is not set" ) );
    }

    if( v_heat.size() != f_time_horizon ) {
        assert( v_heat.size() == 0 ); // this should only happen once
        v_heat.resize( f_time_horizon );

        for(int t = 0 ; t < f_time_horizon - 1 ; t++){

            v_heat[t].set_type(ColVariable::kNonNegative);

        }
        add_static_variable ( v_heat);
    }
*/
/*--------------------------------------------------------------------------*/

    if( v_heat.size() != 0 ||
        v_heat_added.size() != 0 ||
        v_heat_removed.size() != 0 ||
        v_heat_available.size() != 0 ) {
        // the abstract variables should be generated only once
        return;
    }

    if( f_time_horizon == 0 ) {
        // there are no variables to be generated
        return;
    }
    typedef std::vector< std::pair< std::vector<ColVariable> * , int > > v_pairs;

    v_pairs variables_and_types = {
     std::make_pair( &v_heat_added,               ColVariable::kNonNegative ),
     std::make_pair( &v_heat_removed,             ColVariable::kNonNegative ),
     std::make_pair( &v_heat_available,           ColVariable::kNonNegative ),
     std::make_pair( &v_heat,                     ColVariable::kNonNegative )
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

/*----------------Satisfaction-Heat-demand-constraints----------------------*/

if (f_number_heat_units > 0 ){

    if( v_HeatDemand_Const.size() != f_time_horizon ) {
        // this should only happen once
        assert( v_HeatDemand_Const.size() == 0 );

        v_HeatDemand_Const.resize( f_time_horizon );
    }

    for (int t = 0; t < f_time_horizon; ++t){

        v_HeatDemand_Const[t].set_rhs(Inf<double>() );
        v_HeatDemand_Const[t].set_lhs( v_heat_demand[t]);
        v_HeatDemand_Const[t].set_function( new LinearFunction() );

        auto linear_function = new LinearFunction();

        linear_function->add_variable( & v_heat_added[t],             -1.0 );
        linear_function->add_variable( & v_heat_removed[t],            1.0 );

        auto heat_variable = get_heat( t );

        for( Index unit_id = 0; unit_id < f_number_heat_units; ++unit_id ) {
            auto linear_function = static_cast<LinearFunction *>
            ( v_HeatDemand_Const[t].get_function() );
            linear_function->add_variable( & heat_variable[unit_id], 1.0 );

        }

    }

    add_static_constraint( v_HeatDemand_Const );

}
/*----------------Satisfaction-Heat-Bounds-constraints----------------------*/

    if (f_number_heat_units > 0 ) {

        if (v_HeatBounds_Const.size() != f_time_horizon) {
            // this should only happen once
            assert(v_HeatBounds_Const.size() == 0);

            v_HeatBounds_Const.resize
                    ( boost::multi_array<FRowConstraint, 2>::
                      extent_gen()[ f_time_horizon ][ f_number_heat_units ] );

        }
        for (int t = 0; t < f_time_horizon; ++t) {
          for( Index unit_id = 0; unit_id < f_number_heat_units; ++unit_id ) {

                v_HeatBounds_Const[t][unit_id].set_rhs( v_max_heat_production
                                                        [unit_id]);
                v_HeatBounds_Const[t][unit_id].set_lhs( v_min_heat_production
                                                        [unit_id]);
                v_HeatBounds_Const[t][unit_id].set_function
                                          ( new LinearFunction() );

                auto heat_variable = get_heat( unit_id );
                auto linear_function = new LinearFunction();

                linear_function->add_variable( & heat_variable[t],      1.0 );


            }

            }
        add_static_constraint( v_HeatBounds_Const );

}
/*------------Satisfaction-Heat-Storage-Bounds-constraints------------------*/

        if( v_HeatStorageBounds_Const.size() != f_time_horizon ) {
            // this should only happen once
            assert( v_HeatStorageBounds_Const.size() == 0 );

            v_HeatStorageBounds_Const.resize( f_time_horizon );
        }

        for (int t = 0; t < f_time_horizon; ++t){

            v_HeatStorageBounds_Const[t].set_rhs(v_max_heat_storage[t] );
            v_HeatStorageBounds_Const[t].set_lhs( v_min_heat_storage[t]);
            v_HeatStorageBounds_Const[t].set_function( new LinearFunction() );

            auto linear_function = new LinearFunction();

            linear_function->add_variable( & v_heat_available[t],       1.0 );


        }

        add_static_constraint( v_HeatStorageBounds_Const );

//TODO ADD other constraints

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
/*
template<class T>
void serialize( netCDF::NcGroup & group, const std::string & var_name,
                const netCDF::NcType & ncType, const netCDF::NcDim & ncDim,
                const std::vector<T> & data ) {

    group.addVar( var_name , ncType , ncDim )
            .putVar( { 0 } , { data.size() } , data.data() );
}
*/
template<class T>
inline void serialize( netCDF::NcGroup & group, const std::string & var_name,
                       const netCDF::NcType & ncType,
                       const std::vector<netCDF::NcDim> & ncDim,
                       const std::vector<T> & data ) {

    std::vector<size_t> start;
    start.assign( ncDim.size(), 0 );

    std::vector<size_t> sizes;
    sizes.resize( ncDim.size() );
    for( size_t i = 0; i < sizes.size(); ++i )
        sizes[i] = ncDim[i].getSize();

    auto total_size = std::accumulate( begin( sizes ), end( sizes ), 1,
                                       std::multiplies<size_t>() );

    if( total_size == 0 )
        return;

    group.addVar( var_name , ncType , ncDim )
            .putVar( start , sizes , data.data() );
}
/*--------------------------------------------------------------------------*/

void HeatBlock::serialize( netCDF::NcGroup & group ) const {

    group.putAtt( "type" , "HeatBlock" );

    auto dim_time_horizon = group.addDim( "TimeHorizon",     f_time_horizon );

    auto dim_number_units = group.addDim( "NumberHeatUnits",
                                                        f_number_heat_units );
    auto dim_number_intervals = group.addDim( "NumberIntervals",
                                                         f_number_intervals );

    ::serialize( group, "StoringHeatRho",   netCDF::NcDouble(),
                                                         f_storing_heat_rho );

    ::serialize( group, "ExtractingHeatRho",   netCDF::NcDouble(),
                                                      f_extracting_heat_rho );

    ::serialize( group, "StoringHeatRho",   netCDF::NcDouble(),
                                                         f_keeping_heat_rho );


    ::serialize( group, "ChangeInterval", netCDF::NcUint64(),
                 {dim_number_intervals}, v_change_interval );


    ::serialize( group, "TotalHeatDemand", netCDF::NcDouble(),
                 {dim_number_intervals}, v_heat_demand);

    ::serialize( group, "MinHeatStorage", netCDF::NcDouble(),
                 {dim_number_intervals},v_min_heat_storage);

    ::serialize( group, "MaxHeatStorage", netCDF::NcDouble(),
                 {dim_number_intervals}, v_max_heat_storage);

    ::serialize( group, "MinHeatProduction", netCDF::NcDouble(),
                 {dim_number_intervals , dim_number_units},
                                                       v_min_heat_production);

    ::serialize( group, "MaxHeatProduction",netCDF::NcDouble(),
                 {dim_number_intervals , dim_number_units},
                                                       v_max_heat_production);


}  // end( HeatBlock::serialize )

/*--------------------------------------------------------------------------*/
/*----------------------- End File HeatBlock.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
