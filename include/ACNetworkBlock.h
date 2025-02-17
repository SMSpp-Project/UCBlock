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

explicit ACNetworkBlock( Block * f_block = nullptr )
  : DCNetworkBlock( f_block ) {}

 void generate_abstract_variables( Configuration * stvv = nullptr ) override;
 void generate_abstract_constraints( Configuration * stcc = nullptr ) override;
 void generate_objective( Configuration * objc = nullptr ) override;

 void generate_SOCP_relaxation();


 const std::vector< ColVariable > & get_power_flow_imag( void ) const { return( v_power_flow_imag ); };

 
 protected:

 // ----- Variables
 std::vector< ColVariable > v_power_flow_imag; // real part is the standard "v_power_flow" variable

 // ----- Generic variables for AC-OPF
 std::vector< ColVariable > v_sum_product_voltages;
 std::vector< ColVariable > v_diff_product_voltages;
 std::vector< ColVariable > v_sqrd_voltages;

 // ----- Specific variables for SOCP relaxation
 SpVarMat W_voltage; //< Sparse matrix (only defined for lines and reversed lines)

 // ----- Generic constraints for AC-OPF
 std::vector< BoxConstraint > v_voltage_bounds_const;
 boost::multi_array< FRowConstraint , 2 > v_angle_bounds_const;
 boost::multi_array< FRowConstraint , 2 > v_voltage_definition_const;
 std::vector< FRowConstraint > v_thermal_limit;

 // ----- Specific constraints for SOCP relaxation
 std::vector< FRowConstraint > v_socp_const;

 private:

 SMSpp_insert_in_factory_h;


 static void static_initialization( void ) {

  register_method< ACNetworkBlock , MF_dbl_it , Subset && , bool >(
   "DCNetworkBlock::set_active_demand" , &ACNetworkBlock::set_active_demand );

  register_method< ACNetworkBlock , MF_dbl_it , Range >(
   "DCNetworkBlock::set_active_demand" , &ACNetworkBlock::set_active_demand );
 }

};  // end( class( ACNetworkBlock ) )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

}  // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* ACNetworkBlock.h included */

/*--------------------------------------------------------------------------*/
/*-------------------- End File ACNetworkBlock.h ---------------------------*/
/*--------------------------------------------------------------------------*/
