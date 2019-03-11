/*--------------------------------------------------------------------------*/
/*------------------------ File UnitBlock.cpp ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the UnitBlock class.
 *
 * \version 0.11
 *
 * \date 05 - 03 - 2019
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
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

void UnitBlock::generate_abstract_variables( Configuration *stvv ) {

  if( v_power.size() != f_time_horizon ) { // TODO It works only if v_power is
                                           // always there

    assert( v_commitment.size() == 0 &&
            v_power.size() == 0 &&
            v_primary_spinning_reserve.size() == 0 &&
            v_secondary_spinning_reserve.size() == 0 &&
            v_heat.size() == 0 ); // this should only happen once

    auto tstvv = dynamic_cast<SimpleConfiguration<int> *>( stvv );

    if( ( ! tstvv ) && f_BlockConfig &&
        f_BlockConfig->f_static_variables_Configuration )

      tstvv = dynamic_cast<SimpleConfiguration<int> *>
        ( f_BlockConfig->f_static_constraints_Configuration );

    std::vector< std::vector<ColVariable> * > variables = {
      &v_commitment,
      &v_power,
      &v_primary_spinning_reserve,
      &v_secondary_spinning_reserve,
      &v_heat
    };

    int use_variables = tstvv->f_value;

    int k = 1;
    for( int i = 0; i < variables.size(); ++i, k *= 2 ) {
      if( use_variables & k ) {
        ( * variables[i] ).resize( f_time_horizon );
        add_static_variable( * variables[i] );
      }
    }
  }
}

/*--------------------------------------------------------------------------*/
/*---------------------- End File UnitBlock.cpp ----------------------------*/
/*--------------------------------------------------------------------------*/
