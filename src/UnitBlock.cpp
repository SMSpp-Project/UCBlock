/*--------------------------------------------------------------------------*/
/*------------------------ File UnitBlock.cpp ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the UnitBlock class.
 *
 * \version 0.11
 *
 * \date 07 - 06 - 2019
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
 * \author Kostas Tavlaridis-Gyparakis \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Antonio Frangioni, Ali Ghezelsoflu, Rafael
 * Durbano Lobato, and Kostas Tavlaridis-Gyparakis
 */

/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Serialization.h"
#include "UCBlock.h"
#include "UnitBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register UnitBlock to the Block factory

SMSpp_insert_in_factory_cpp_1( UnitBlock );

/*--------------------------------------------------------------------------*/
/*--------------------------- METHODS OF UnitBlock -------------------------*/
/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void UnitBlock::deserialize_time_horizon( netCDF::NcGroup & group ) {

  netCDF::NcDim TimeHorizon = group.getDim( "TimeHorizon" );
  if( TimeHorizon.isNull() ) {

    // Dimension TimeHorizon is not present in the netCDF input.

    auto f_Block = get_f_Block();

    if( f_Block ) {

      // The father Block is available. Take time horizon from it.
      auto time_horizon_father =
        static_cast<UCBlock *>( get_f_Block() )->get_time_horizon();

      if( f_time_horizon == 0 )
        this->set_time_horizon( time_horizon_father );

      else if( f_time_horizon != time_horizon_father )
        throw( std::logic_error
               ( "UnitBlock::deserialize: TimeHorizon is not present in the "
                 "netCDF. The (nonzero) time horizon of UnitBlock is different "
                 "from that of its father, but they should be equal." ) );
    }
    else if( f_time_horizon != 0 )
      throw( std::invalid_argument
             ( "UnitBlock::deserialize: TimeHorizon is not present in the "
               "netCDF input and UnitBlock does not have a father." ) );
  }
  else {

    auto time_horizon_netcdf = TimeHorizon.getSize();

    if( true ) { // TODO The condition should be that this object was
                 // just created and only the time horizon may have
                 // been set so far.

      if( f_time_horizon == 0 )
        // Use the time horizon provided by the netCDF
        f_time_horizon = time_horizon_netcdf;

      else if( f_time_horizon != time_horizon_netcdf )
        throw( std::invalid_argument
               ( "UnitBlock::deserialize: TimeHorizon in netCDF is different "
                 "from that (nonzero) currently specified in UnitBlock." ) );
    }
    else {
      // Replace the current time horizon (and possibly reset this UnitBlock)
      f_time_horizon = time_horizon_netcdf;
    }
  }
}

/*--------------------------------------------------------------------------*/

void UnitBlock::deserialize_change_intervals( netCDF::NcGroup & group ) {

  auto NumberIntervals = group.getDim( "NumberIntervals" );
  if( NumberIntervals.isNull() )
    f_number_intervals = 0;
  else {
    f_number_intervals = NumberIntervals.getSize();
    if( ( f_number_intervals < 1 ) || ( f_number_intervals > f_time_horizon ) )
      throw( std::invalid_argument
             ( "UnitBlock::deserialize: invalid NumberIntervals. "
               "It must be between 1 and TimeHorizon." ) );
  }

  if( ( f_number_intervals > 1 ) && ( f_number_intervals < f_time_horizon ) ) {

    Serialization::deserialize( group, "ChangeIntervals", f_number_intervals,
                                v_change_intervals );

    // Check that all numbers are between 1 and f_time_horizon, that
    // the last number is == f_time_horizon, and that they are ordered
    // in increasing sense

    if( v_change_intervals.back() != f_time_horizon ) {
      throw( std::invalid_argument
             ( "UnitBlock::deserialize: invalid value in ChangeIntervals: "
               "the last element must be TimeHorizon." ) );
    }

    Index previous_t = 0;

    for( auto t : v_change_intervals ) {
      if( ! ( t > previous_t && t < f_time_horizon - 1 ) )
        throw( std::invalid_argument
               ( "UnitBlock::deserialize: invalid value in ChangeIntervals: " +
                 std::to_string( t ) + ". All values must be between 1 and "
                 "TimeHorizon and in strictly increasing order." ) );

      previous_t = t;
    }
  }
}

/*--------------------------------------------------------------------------*/

void UnitBlock::deserialize( netCDF::NcGroup & group ) {

  using SMSpp_di_unipi_it::Serialization::deserialize;

  guts_of_destructor();
  deserialize_time_horizon( group );
  deserialize_change_intervals( group );

  deserialize( group, "FixedConsumption", v_fixed_consumption,
               {f_number_intervals} );

  deserialize( group, "InertiaCommitment", v_inertia_commitment,
               {f_number_intervals} );

  deserialize( group, "InertiaPower", v_inertia_power,
               {f_number_intervals} );
}

/*--------------------------------------------------------------------------*/

int UnitBlock::get_variables_to_be_generated( Configuration *stvv ) {

  if( ! stvv )
    return 0;

  // informs which variables must be generated
  int variables_to_be_generated = 0;

  auto tstvv = dynamic_cast<SimpleConfiguration<int> *>( stvv );

  if( ( ! tstvv ) && f_BlockConfig &&
      f_BlockConfig->f_static_variables_Configuration ) {

    tstvv = dynamic_cast<SimpleConfiguration<int> *>
      ( f_BlockConfig->f_static_variables_Configuration );
  }

  if( tstvv )
    variables_to_be_generated = tstvv->f_value;

  return variables_to_be_generated;
}

/*--------------------------------------------------------------------------*/

void UnitBlock::generate_abstract_variables( Configuration *stvv ) {

  if( f_time_horizon == 0 ) {
    // there are no variables to be generated
    return;
  }

  if( v_commitment.size() != 0 ||
      v_primary_spinning_reserve.size() != 0 ||
      v_secondary_spinning_reserve.size() != 0 ||
      v_active_power.size() != 0 ) {
    // the abstract variables should be generated only once
    return;
  }

  typedef std::vector< std::pair< std::vector<ColVariable> * , int > > v_pairs;

  v_pairs variables_and_types = {
    std::make_pair( &v_commitment,                 ColVariable::kBinary ),
    std::make_pair( &v_primary_spinning_reserve,   ColVariable::kNonNegative ),
    std::make_pair( &v_secondary_spinning_reserve, ColVariable::kNonNegative ),
    std::make_pair( &v_active_power,               ColVariable::kNonNegative )
    // v_active_power must be the last one in this list
  };

  auto variables_to_be_generated = get_variables_to_be_generated( stvv );

  // The active power variables must be always present
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
}

/*--------------------------------------------------------------------------*/
/*------------------- METHODS FOR MODIFYING THE UnitBlock ------------------*/
/*--------------------------------------------------------------------------*/

void UnitBlock::set_time_horizon( Index t ) {
  if( f_time_horizon == t )
    return;

  if( f_time_horizon != 0 )
    throw std::logic_error( "UnitBlock::set_time_horizon: "
                            "time horizon has already been set.");

  f_time_horizon = t;
}

/*--------------------------------------------------------------------------*/
/*---------------------- METHODS FOR SAVING THE UnitBlock ------------------*/
/*--------------------------------------------------------------------------*/

void UnitBlock::serialize( netCDF::NcGroup & group ) const {
  group.putAtt( "type" , "UnitBlock" );
  group.addDim( "TimeHorizon" , f_time_horizon );

  using SMSpp_di_unipi_it::Serialization::serialize;


  auto NumberIntervals = group.addDim( "NumberIntervals", f_number_intervals );

  Serialization::serialize( group, "ChangeInterval", netCDF::NcUint64(),
                            NumberIntervals, v_change_intervals );

  //::serialize( group, "FixedConsPower", netCDF::NcDouble(), f_FixedConsPower );

  serialize( group, "FixedConsumption", netCDF::NcDouble(),
             {NumberIntervals}, v_fixed_consumption);

  serialize( group, "InertiaCommitment", netCDF::NcDouble(),
             {NumberIntervals}, v_inertia_commitment);

  serialize( group, "InertiaPower", netCDF::NcDouble(),
             {NumberIntervals}, v_inertia_power);
}

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

void UnitBlock::guts_of_destructor( void ) {

  // delete all Variables
  v_commitment.clear();
  v_active_power.clear();
  v_primary_spinning_reserve.clear();
  v_secondary_spinning_reserve.clear();

  // explicitly reset all Variables

  // this is done for the case where this method is called prior to
  // re-loading a new instance: if not, the new representation would
  // be added to the previous one
  reset_static_variables();
}

/*--------------------------------------------------------------------------*/
/*---------------------- End File UnitBlock.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
