/*--------------------------------------------------------------------------*/
/*------------------------- File HeatBlock.cpp -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the HeatBlock class.
 *
 * \version 0.11
 *
 * \date 19 - 06 - 2019
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

void HeatBlock::deserialize_time_horizon( netCDF::NcGroup & group )
{
 netCDF::NcDim TimeHorizon = group.getDim( "TimeHorizon" );
    if( TimeHorizon.isNull() ) {
        // dimension TimeHorizon is not present in the netCDF input

        if( f_time_horizon == 0 ) {
            auto f_B = dynamic_cast< UCBlock * >( get_f_Block() );
            if( f_B )
                // The father Block is available. Take time horizon from it.
                this->set_time_horizon( f_B->get_time_horizon() );
            else
                throw( std::invalid_argument(
                    "HeatBlock::deserialize: TimeHorizon is not present in the "
                    "netCDF input and HeatBlock does not have a father." ) );
        }
    }
    else {
        // dimension TimeHorizon is present in the netCDF input

        auto th = TimeHorizon.getSize();
        if( f_time_horizon == 0 )
            this->set_time_horizon( th );
        else
        if( f_time_horizon != th )
            throw( std::logic_error(
                "HeatBlock::deserialize: TimeHorizon is not present in the "
                "netCDF. The (nonzero) time horizon of HeatBlock is different "
                "from that of its father, but they should be equal." ) );
    }
}
/*--------------------------------------------------------------------------*/

void HeatBlock::deserialize_change_intervals( netCDF::NcGroup & group ) {

    auto NumberIntervals = group.getDim( "NumberIntervals" );
    if( NumberIntervals.isNull() )
        f_number_intervals = 0;
    else {
        f_number_intervals = NumberIntervals.getSize();
        if( ( f_number_intervals < 1 ) || ( f_number_intervals > f_time_horizon ))
            throw( std::invalid_argument
                ( "HeatBlock::deserialize: invalid NumberIntervals. "
                  "It must be between 1 and TimeHorizon." ) );
    }

    if( ( f_number_intervals > 1 ) && ( f_number_intervals < f_time_horizon ) ){

      ::deserialize( group, "ChangeIntervals", f_number_intervals,
                     v_change_intervals );

        // Check that all numbers are between 1 and f_time_horizon, that
        // the last number is == f_time_horizon, and that they are ordered
        // in increasing sense

        if( v_change_intervals.back() != f_time_horizon ) {
            throw( std::invalid_argument
                ( "HeatBlock::deserialize: invalid value in ChangeIntervals: "
                  "the last element must be TimeHorizon." ) );
        }

        Index previous_t = 0;

        for( auto t : v_change_intervals ) {
            if( ! ( t > previous_t && t < f_time_horizon - 1 ) )
                throw( std::invalid_argument
                    ( "HeatBlock::deserialize: invalid value in ChangeIntervals: " +
                      std::to_string( t ) + ". All values must be between 1 and "
                                            "TimeHorizon and in strictly increasing order." ) );

            previous_t = t;
        }
    }
}
/*--------------------------------------------------------------------------*/

void HeatBlock::deserialize( netCDF::NcGroup & group ) {

    deserialize_time_horizon( group );
    deserialize_change_intervals( group );

    ::deserialize( group, "TotalHeatDemand", f_time_horizon,
                   v_heat_demand );

    ::deserialize( group, "CostHeatUnit",
                   { f_number_intervals , f_number_heat_units },
                   v_cost_heat_unit );

    ::deserialize( group, "MinHeatProduction",
                   { f_number_intervals , f_number_heat_units },
                   v_min_heat_production );

    ::deserialize( group, "MaxHeatProduction",
                   { f_number_intervals , f_number_heat_units },
                   v_max_heat_production );

    ::deserialize( group, "MinHeatStorage", f_number_intervals,
                   v_min_heat_storage );

    ::deserialize( group, "MaxHeatStorage", f_number_intervals,
                   v_max_heat_storage );

    ::deserialize( group, "StoringHeatRho",       & f_storing_heat_rho );
    ::deserialize( group, "ExtractingHeatRho",    & f_extracting_heat_rho );
    ::deserialize( group, "KeepingHeatRho",       & f_keeping_heat_rho );
    ::deserialize( group, "InitialHeatAvailable", & f_initial_heat_storage );

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

    // Satisfaction Heat Demand constraints

    if ( f_number_heat_units > 0 ){

        if( v_HeatDemand_Constraints.size() != f_time_horizon ) {
            // this should only happen once
            assert( v_HeatDemand_Constraints.size() == 0 );

            v_HeatDemand_Constraints.resize( f_time_horizon );
        }

        for ( Index t = 0; t < f_time_horizon; ++t ){

            v_HeatDemand_Constraints[ t ].set_rhs( Inf<double>() );
            v_HeatDemand_Constraints[ t ].set_lhs( v_heat_demand[ t ]);
            v_HeatDemand_Constraints[ t ].set_function( new LinearFunction() );

            auto linear_function = new LinearFunction();

            linear_function->add_variable( & v_heat_added[ t ],               -1.0 );
            linear_function->add_variable( & v_heat_removed[ t ],              1.0 );

            auto heat_variable = get_heat( t );

            for( Index unit_id = 0; unit_id < f_number_heat_units; ++unit_id ) {
                auto linear_function = static_cast<LinearFunction *>
                ( v_HeatDemand_Constraints[ t ].get_function() );
                linear_function->add_variable( & heat_variable[ unit_id ],       1.0 );
            }
        }

        add_static_constraint( v_HeatDemand_Constraints );
    }

    // Satisfaction Heat Bounds constraints

    if ( f_number_heat_units > 0 ) {

        if ( v_HeatBounds_Constraints.size() != f_time_horizon ) {
            // this should only happen once
            assert( v_HeatBounds_Constraints.size() == 0 );

            v_HeatBounds_Constraints.resize
                ( boost::multi_array<FRowConstraint, 2>::
                  extent_gen()[ f_time_horizon ][ f_number_heat_units ] );

        }
        for ( Index t = 0; t < f_time_horizon; ++t ) {
            for( Index unit_id = 0; unit_id < f_number_heat_units; ++unit_id ) {

                v_HeatBounds_Constraints[ t ][ unit_id ].set_rhs
                    ( v_max_heat_production[ unit_id ]);
                v_HeatBounds_Constraints[ t ][ unit_id ].set_lhs
                    ( v_min_heat_production[ unit_id ]);
                v_HeatBounds_Constraints[ t ][ unit_id ].set_function
                    ( new LinearFunction() );

                auto heat_variable = get_heat( unit_id );
                auto linear_function = new LinearFunction();

                linear_function->add_variable( & heat_variable[ t ],           1.0 );
            }
        }
        add_static_constraint( v_HeatBounds_Constraints );
    }

    // Satisfaction Heat Storage Bounds constraints

    if( v_HeatStorageBounds_Constraints.size() != f_time_horizon ) {
        // this should only happen once
        assert( v_HeatStorageBounds_Constraints.size() == 0 );

        v_HeatStorageBounds_Constraints.resize( f_time_horizon );
    }

    for ( Index t = 0; t < f_time_horizon; ++t ){

        v_HeatStorageBounds_Constraints[ t ].set_rhs( v_max_heat_storage[ t ] );
        v_HeatStorageBounds_Constraints[ t ].set_lhs( v_min_heat_storage[ t ] );
        v_HeatStorageBounds_Constraints[ t ].set_function( new LinearFunction() );

        auto linear_function = new LinearFunction();

        linear_function->add_variable( & v_heat_available[ t ],            1.0 );
    }
    add_static_constraint( v_HeatStorageBounds_Constraints );

    // Satisfaction Heat Storage Bounds constraints
    if( v_EvolutionStoredHeat_Constraints.size() != f_time_horizon ) {
        // this should only happen once
        assert( v_EvolutionStoredHeat_Constraints.size() == 0 );

        v_EvolutionStoredHeat_Constraints.resize( f_time_horizon );
    }

    for ( Index t = 0; t < f_time_horizon; ++t ){

        if( t == 0 ) {
            v_EvolutionStoredHeat_Constraints[ 0 ].set_lhs ( 0.0 );
            v_EvolutionStoredHeat_Constraints[ 0 ].set_rhs
                ( f_initial_heat_storage * f_keeping_heat_rho );
            v_EvolutionStoredHeat_Constraints[ 0 ].set_function
                ( new LinearFunction() );

            auto linear_function = new LinearFunction();

            linear_function->add_variable( & v_heat_available[ 0 ],           1.0 );

            linear_function->add_variable( & v_heat_added[ 0 ],
                                           - f_storing_heat_rho );
            linear_function->add_variable( & v_heat_removed[ 0 ],
                                           f_extracting_heat_rho );
        }
        else {
            v_EvolutionStoredHeat_Constraints[ t ].set_both ( 0.0 );
            v_EvolutionStoredHeat_Constraints[ t ].set_function
                ( new LinearFunction() );

            auto linear_function = new LinearFunction();
            linear_function->add_variable( & v_heat_available[ t +1 ],        1.0 );
            linear_function->add_variable( & v_heat_available[ t ],
                                           - f_keeping_heat_rho );
            linear_function->add_variable( & v_heat_added[ t + 1 ],
                                           - f_storing_heat_rho );
            linear_function->add_variable( & v_heat_removed[ t  + 1],
                                           f_extracting_heat_rho );
        }
    }
    add_static_constraint( v_EvolutionStoredHeat_Constraints );

} // end( HeatBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void HeatBlock::generate_objective( Configuration *objc ) {

    if( ! get_objective().empty() )  // an objective is there already
        return;                        // cowardly (and silently) return

    // Initialize objective function

    if( v_heat.size() != f_time_horizon ) {
        throw( std::logic_error
            ( "HeatBlock::generate_objective: v_heat must have "
              "size equal to the time horizon." ) );
    }

    auto linear_function = new LinearFunction();

    for( Index unit_id = 0; unit_id < f_number_heat_units; ++unit_id ) {
        for ( Index t = 0; t < f_time_horizon; ++t ) {

            auto heat = get_heat( t );

            auto cost = get_cost_heat_unit( t , unit_id);

            linear_function->add_variable( & heat[ t ], cost );

        }
    }

    objective.set_function( linear_function );
    objective.set_sense( Objective::eMin );
    objective.set_Block( this );

}  // end( HeatBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*---------- METHODS FOR LOADING, PRINTING & SAVING THE HeatBlock ----------*/
/*--------------------------------------------------------------------------*/

void HeatBlock::serialize( netCDF::NcGroup & group ) const {

    group.putAtt( "type" , "HeatBlock" );
    group.addDim( "TimeHorizon" , f_time_horizon );

    auto NumberIntervals = group.addDim( "NumberIntervals", f_number_intervals);

    ::serialize( group, "ChangeInterval", netCDF::NcUint64(),
                 NumberIntervals, v_change_intervals );

    auto dim_number_units = group.addDim( "NumberHeatUnits",
                                          f_number_heat_units );

    ::serialize( group, "StoringHeatRho", netCDF::NcDouble(),
                 f_storing_heat_rho );

    ::serialize( group, "ExtractingHeatRho", netCDF::NcDouble(),
                 f_extracting_heat_rho );

    ::serialize( group, "StoringHeatRho", netCDF::NcDouble(),
                 f_keeping_heat_rho );

    ::serialize( group, "InitialHeatAvailable", netCDF::NcDouble(),
                 f_initial_heat_storage );

    ::serialize( group, "ChangeInterval", netCDF::NcUint64(),
                 {NumberIntervals}, v_change_intervals );

    ::serialize( group, "TotalHeatDemand", netCDF::NcDouble(),
                 {NumberIntervals}, v_heat_demand);

    ::serialize( group, "MinHeatStorage", netCDF::NcDouble(),
                 {NumberIntervals},v_min_heat_storage);

    ::serialize( group, "MaxHeatStorage", netCDF::NcDouble(),
                 {NumberIntervals}, v_max_heat_storage);

    ::serialize( group, "MinHeatProduction", netCDF::NcDouble(),
                 {NumberIntervals , dim_number_units},
                 v_min_heat_production);

    ::serialize( group, "MaxHeatProduction",netCDF::NcDouble(),
                 {NumberIntervals , dim_number_units},
                 v_max_heat_production);

}  // end( HeatBlock::serialize )

/*--------------------------------------------------------------------------*/
/*----------------------- End File HeatBlock.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
