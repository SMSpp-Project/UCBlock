/*--------------------------------------------------------------------------*/
/*------------------------ File UnitBlock.cpp ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the UnitBlock class.
 *
 * \version 0.11
 *
 * \date 21 - 05 - 2019
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

void UnitBlock::deserialize( netCDF::NcGroup & group ) {

  netCDF::NcDim TH = group.getDim( "TimeHorizon" );
  if( TH.isNull() ) {

    // Dimension TimeHorizon is not present in the netCDF input.

    auto f_Block = get_f_Block();

    if( f_Block ) {
      // The father Block is available. Take time horizon from it.
      auto time_horizon = static_cast<UCBlock *>( get_f_Block() )->
        get_time_horizon();
      this->set_time_horizon( time_horizon );
    }
    else
      throw( std::invalid_argument( "UnitBlock::deserialize: TimeHorizon "
                                    "is not present in the netCDF input" ) );
  }
  else {

    f_time_horizon = TH.getSize();
    auto f_Block = get_f_Block();

    if( f_Block ) {

      // Check whether the given time horizon is equal to that of the
      // father Block.
      
      auto time_horizon = static_cast<UCBlock *>( get_f_Block() )->
        get_time_horizon();

      if( f_time_horizon != time_horizon )
        throw( std::invalid_argument
               ( "UnitBlock::deserialize: TimeHorizon is different "
                 "from that specified in UCBlock" ) );
    }
    else if( f_time_horizon <= 0 )
      throw( std::invalid_argument
             ( "UnitBlock::deserialize: TimeHorizon must be positive" ) );
  }
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
      ( f_BlockConfig->f_static_constraints_Configuration );
  }

  if( tstvv )
    variables_to_be_generated = tstvv->f_value;

  return variables_to_be_generated;
}

/*--------------------------------------------------------------------------*/

void UnitBlock::generate_abstract_variables( Configuration *stvv ) {

  if( v_commitment.size() != 0 ||
      v_power_injected.size() != 0 ||
      v_active_power.size() != 0 ||
      v_primary_spinning_reserve.size() != 0 ||
      v_secondary_spinning_reserve.size() != 0) {
    // the abstract variables should be generated only once
    return;
  }

  if( f_time_horizon == 0 ) {
    // there are no variables to be generated
    return;
  }

  typedef std::vector< std::pair< std::vector<ColVariable> * , int > > v_pairs;

  v_pairs variables_and_types = {
    std::make_pair( &v_commitment,                 ColVariable::kBinary ),
    std::make_pair( &v_power_injected,             ColVariable::kNonNegative ),
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

void UnitBlock::set_time_horizon( int t ) {
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

  //TODO InertiaCommitment and InertiaPower
}

/*--------------------------------------------------------------------------*/
/*---------------------- End File UnitBlock.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
