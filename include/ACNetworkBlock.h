/*--------------------------------------------------------------------------*/
/*--------------------------- File ACNetworkBlock.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 *
 * Header file for class ACNetworkBlock, which derives from DCNetworkBlock and
 * defines the standard SOCP relaxation corresponding to the "AC model"
 * of the transmission network in the Unit Commitment problem.
 *
 * @ Quentin J.
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __ACNetworkBlock
 #define __ACNetworkBlock
                      /* self-identification: #endif at the end of the file */

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Block.h"

#include "LinearFunction.h"

#include "FRowConstraint.h"

#include "OneVarConstraint.h"

#include "DCNetworkBlock.h"

#include "FRealObjective.h"

#include <iostream>

#include <Eigen/Sparse>
/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)

typedef Eigen::SparseMatrix< std::complex<double> > SpCMat;
typedef Eigen::SparseVector< std::complex<double> > SpCVec;

namespace SMSpp_di_unipi_it
{



class ACNetworkBlock : public DCNetworkBlock
{

 public:

  struct ACNetworkData 
  {
    // classical quantities for AC Network
    SpCMat Yff, Yft, Ytf, Ytt;

    // line shunt
    SpCVec Ys;

    // matrix for power flow equations
    SpCMat M, HM, ZM;
  };

  struct ComplexColVariable
  {
    ColVariable real;
    ColVariable imag;
  };

explicit ACNetworkBlock( Block * f_block = nullptr )
  : DCNetworkBlock( f_block ) {}

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;
 void generate_abstract_constraints( Configuration * stcc = nullptr ) override;

virtual void add_ACdata(Index interval, Index node, UnitBlock* unit_block, Index t, Index g ) override;

 
 protected:

 ACNetworkData ACdata;
 boost::multi_array< ColVariable , 2 > W_voltage_real;
 boost::multi_array< ColVariable , 2 > W_voltage_imag;

 std::vector< BoxConstraint > v_voltage_bounds_const;

 private:

 SMSpp_insert_in_factory_h;


 static void static_initialization( void ) {

  register_method< ACNetworkBlock , MF_dbl_it , Subset && , bool >(
   "DCNetworkBlock::set_active_demand" , &ACNetworkBlock::set_active_demand );

  register_method< ACNetworkBlock , MF_dbl_it , Range >(
   "DCNetworkBlock::set_active_demand" , &ACNetworkBlock::set_active_demand );
 }

};  // end( class( DCNetworkBlock ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* ACNetworkBlock.h included */

/*--------------------------------------------------------------------------*/
/*-------------------- End File ACNetworkBlock.h ---------------------------*/
/*--------------------------------------------------------------------------*/
