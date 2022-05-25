/*--------------------------------------------------------------------------*/
/*-------------------------- File BusNetworkBlock.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the class BusNetworkBlock, which derives from NetworkBlock
 * [see NetworkBlock.h] in order to define a "bus" transmission network in
 * the Unit Commitment problem. There is actually precious little that this
 * class has to do that is not done already by the base NetworkBlock class.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Ali Ghezelsoflu \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Rafael Durbano Lobato \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Kostas Tavlaridis-Gyparakis \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \copyright &copy; by Antonio Frangioni, Ali Ghezelsoflu,
 *                   Rafael Durbano Lobato
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __BusNetworkBlock
 #define __BusNetworkBlock
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "NetworkBlock.h"
#include "OneVarConstraint.h"
#include "FRealObjective.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it
{
/*--------------------------------------------------------------------------*/
/*------------------------ CLASS BusNetworkBlock ---------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/
/// A NetworkBlock with only one node, i.e., a "bus" transmission network
/** The BusNetworkBlock class, which derives from NetworkBlock [see
 * NetworkBlock.h] implements the Block concept [see Block.h] in order to
 * define a "bus" transmission network in the Unit Commitment problem for a
 * given instant in the time horizon. There is actually precious little that
 * this class has to do that is not done already by the base NetworkBlock
 * class.*/

class BusNetworkBlock : public NetworkBlock
{
/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 public:

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

 /// constructor of BusNetworkBlock, taking possibly a pointer to its father

 explicit BusNetworkBlock( Block * f_block = nullptr )
  : NetworkBlock( f_block ) {
  v_node_injection.resize( 1 );
 }

/*--------------------------------------------------------------------------*/
 /// destructor of BusNetworkBlock, (understandably) does nothing

 virtual ~BusNetworkBlock() override = default;

/**@} ----------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

 /// loads the BusNetworkBlock instance from file - not implemented yet

 void load( std::istream & input , char frmt = 0 ) override {
  throw( std::logic_error( "BusNetworkBlock::load() not implemented yet" ) );
 }

/*--------------------------------------------------------------------------*/
 /// generates the static variables of BusNetworkBlock
 /** The base BusNetworkBlock class has just the node injection variables.
  * Since a "bus" network has just one node, and therefore a single value D for
  * the demand and a single injection variable s, which can hardly be called a
  * variable since the only possible way to satisfy the constraints is by
  * having s = D which in fact makes the variable a constant. */

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;

/*--------------------------------------------------------------------------*/
 /// Generate the static constraint of the BusNetworkBlock
 /** This method generates the abstract constraints of the BusNetworkBlock.
  * Since the node injection variable is fixed to the active demand value, it
  * must be a BoxConstraint for that variable whose lower and upper bounds are
  * equal to the active demand value. */

 void generate_abstract_constraints( Configuration * stcc = nullptr ) override;

/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
 /// generate the objective of the BusNetworkBlock
 /** Method that generates the objective of the BusNetworkBlock.
  *
  * - Objective function: the objective function of the BusNetworkBlock
  *   is "empty" (a FRealObjective with a LinearFunction inside with no active
  *   variables) */

 void generate_objective( Configuration * objc = nullptr ) override;

/**@} ----------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

 void set_active_demand( std::vector< double >::const_iterator values ,
                         Subset && subset = { 0 } , bool ordered = false ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck ) override final;

 void set_active_demand( std::vector< double >::const_iterator values ,
                         Range rng = Range( 0 , 1 ) ,
                         ModParam issuePMod = eNoBlck ,
                         ModParam issueAMod = eNoBlck ) override final;

 static void static_initialization( void ) {

  /* Warning: Not all C++ compilers enjoy the template wizardry behind the
   * three-args version of register_method<> with the compact MS_*_*::args(),
   *
   * register_method< BusNetworkBlock >( "BusNetworkBlock::set_active_demand",
   *                                     &BusNetworkBlock::set_active_demand,
   *                                     MS_dbl_sbst::args() );
   *
   * so we just use the slightly less compact one with the explicit argument
   * and be done with it. */

  register_method< BusNetworkBlock , MF_dbl_it , Subset && , bool >(
   "BusNetworkBlock::set_active_demand" ,
   &BusNetworkBlock::set_active_demand );

  register_method< BusNetworkBlock , MF_dbl_it , Range >(
   "BusNetworkBlock::set_active_demand" ,
   &BusNetworkBlock::set_active_demand );
 }

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

 protected:

/*--------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------data--------------------------------------*/



/*-----------------------------variables------------------------------------*/



/*----------------------------constraints-----------------------------------*/

 /// the node injection bound constraints
 std::vector< BoxConstraint > NodeInjection_bound_Constraints;


 /// the objective function
 FRealObjective objective;

/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/
/*--------------------------- PRIVATE FIELDS -------------------------------*/
/*--------------------------------------------------------------------------*/

 SMSpp_insert_in_factory_h;

/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/



};  // end( class( BusNetworkBlock ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* BusNetworkBlock.h included */

/*--------------------------------------------------------------------------*/
/*----------------------- End File BusNetworkBlock.h -----------------------*/
/*--------------------------------------------------------------------------*/
