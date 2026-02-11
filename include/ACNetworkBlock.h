/*--------------------------------------------------------------------------*/
/*--------------------------- File ACNetworkBlock.h ------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 *
 * Header file for class ACNetworkBlock, which derives from DCNetworkBlock and
 * defines the standard SOCP relaxation corresponding to the "AC model"
 * of the transmission network in the Unit Commitment problem.
 *
 * \author Quentin Jacquet \n
 *         EDF R&D OSIRIS \n
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

#include <Eigen/Sparse>

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)

typedef Eigen::SparseMatrix< std::complex< double > > SpCMat;
typedef Eigen::SparseVector< std::complex< double > > SpCVec;

namespace SMSpp_di_unipi_it {
 class ACNetworkBlock : public DCNetworkBlock {
 public:
  class ACNetworkData : public DCNetworkData {
   /*--------------------------------------------------------------------------*/
   /*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
   /*--------------------------------------------------------------------------*/

  public:
   /**@} ----------------------------------------------------------------------*/
   /*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
   /*--------------------------------------------------------------------------*/
   /** @name Constructor and Destructor
    * @{ */

   /// constructor of ACNetworkData, does nothing
   ACNetworkData( void ) {
   }

   /// copy constructor of ACNetworkData, does nothing
   explicit ACNetworkData( const NetworkData * ) {
   }

   /// destructor of ACNetworkData: it is virtual, and empty
   virtual ~ACNetworkData( ) override = default;

   virtual void deserialize( const netCDF::NcGroup & group ) override;

   const std::vector< double > & get_node_conductance( void ) const {
    return( v_node_conductance );
   }

   const std::vector< double > & get_node_max_voltage( void ) const {
    return( v_node_max_voltage );
   }

   const std::vector< double > & get_node_min_voltage( void ) const {
    return( v_node_min_voltage );
   }

   const std::vector< double > & get_line_reactance( void ) const {
    return( v_line_reactance );
   }

   const std::vector< double > & get_line_resistance( void ) const {
    return( v_line_resistance );
   }

   const std::vector< double > & get_line_ratio( void ) const {
    return( v_line_ratio );
   }

   const std::vector< double > & get_line_rate_A( void ) const {
    return( v_line_rate_A );
   }

   const std::vector< double > & get_line_angle( void ) const {
    return( v_line_angle );
   }

   const std::vector< double > & get_line_min_angle( void ) const {
    return( v_line_min_angle );
   }

   const std::vector< double > & get_line_max_angle( void ) const {
    return( v_line_max_angle );
   }

   std::vector< std::pair< std::set< Index >, std::set< Index > > > get_direct_and_reverse_AClines( void ) const {
        std::vector< std::pair< std::set< Index >, std::set< Index > > > v(get_number_nodes());
        const auto & start_line = get_start_line();
        const auto & end_line = get_end_line();
        for (auto& line_id: get_AC_lines()){
            Index p = start_line[ line_id ];
            Index n = end_line[ line_id ];
            v[p].first.insert( line_id );
            v[n].second.insert( line_id );
        }
        return v;
   }


   /*--------------------------------------------------------------------------*/
   /*--------------------- PROTECTED PART OF THE CLASS ------------------------*/
   /*--------------------------------------------------------------------------*/

  protected:
   /*--------------------------------------------------------------------------*/
   /*--------------------- PROTECTED METHODS OF THE CLASS ---------------------*/
   /*--------------------------------------------------------------------------*/

   /*--------------------------------------------------------------------------*/
   /*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
   /*--------------------------------------------------------------------------*/

   std::vector< double > v_line_reactance;
   std::vector< double > v_line_resistance;
   std::vector< double > v_line_ratio;
   std::vector< double > v_line_rate_A;
   std::vector< double > v_line_angle;
   std::vector< double > v_line_min_angle;
   std::vector< double > v_line_max_angle;
   std::vector< double > v_node_conductance;
   std::vector< double > v_node_max_voltage;
   std::vector< double > v_node_min_voltage;

   /*--------------------------------------------------------------------------*/
   /*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
   /*--------------------------------------------------------------------------*/

  private:
   /*--------------------------------------------------------------------------*/
   /*-------------------- PRIVATE FIELDS OF THE CLASS -------------------------*/
   /*--------------------------------------------------------------------------*/

   SMSpp_insert_in_factory_h;

   /*--------------------------------------------------------------------------*/
   /*---------------------- PRIVATE METHODS OF THE CLASS ----------------------*/
   /*--------------------------------------------------------------------------*/
  }; // end( class( ACNetworkData ) )

  explicit ACNetworkBlock( Block * f_block = nullptr )
   : DCNetworkBlock( f_block ), f_NetworkData( nullptr ) {
  }

  void deserialize( const netCDF::NcGroup & group ) override;

  void generate_abstract_variables( Configuration * stvv = nullptr ) override;

  void generate_abstract_constraints( Configuration * stcc = nullptr ) override;

  void generate_objective( Configuration * objc = nullptr ) override;

  // We need to override the functions, as the NetworkData is the one of AC and not the one of DC
  Index get_number_nodes( void ) const override {
   return( f_NetworkData ) ? ( f_NetworkData->get_number_nodes() ) : 1;
  }

  Index get_number_lines( void ) const {
   return( f_NetworkData ) ? ( f_NetworkData->get_number_lines() ) : 0;
  }

  // A rounding function
  inline double round_sig(double value, int digits = 16){
    if (value == 0.0) return 0.0;

    double abs_v = std::fabs(value);
    int exponent = static_cast<int>(std::floor(std::log10(abs_v)));
    double factor = std::pow(10.0, digits - 1 - exponent);

    return std::round(value * factor) / factor;
  }

  std::vector< double > get_line_losses( void ) {
   std::vector< double > losses;
   Index number_lines = get_number_lines();
   for( int line_id = 0 ; line_id < number_lines ; ++line_id ) {
    losses.push_back( v_power_flow[ line_id ].get_value()
     + v_power_flow[ number_lines + line_id ].get_value() );
   }
   return( losses );
  };

  void generate_SOCP_relaxation( );

  const std::vector< ColVariable > & get_reactive_power_flow( void ) const {
   return( v_reactive_power_flow );
  }; // warning only a relaxed solution
  std::vector< std::pair< double , double > > recover_feasible_solution( void );

  /*--------------------------------------------------------------------------*/
  /// returns a pointer to the DCNetworkData
  /** Return a pointer to the DCNetworkData. */

  NetworkData * get_NetworkData( void ) const override {
   return( f_NetworkData );
  }

  /** @} ---------------------------------------------------------------------*/
  /*--------------- METHODS FOR MODIFYING THE ACNetworkBlock -----------------*/
  /*--------------------------------------------------------------------------*/
  /** @name Methods for modifying the ACNetworkBlock
   * @{ */

  void set_NetworkData( NetworkData * nd = nullptr ) override {
   // if there was a previous ACNetworkData, and it was local, delete it
   if( f_NetworkData && f_local_NetworkData )
    delete( f_NetworkData );

   f_NetworkData = dynamic_cast< ACNetworkData * >( nd );
   DCNetworkBlock::set_NetworkData( nd );
   f_local_NetworkData = false;
  }

 protected:

  ACNetworkData * f_NetworkData; ///< the ACNetworkData object

  // ----- Variables
  std::vector< ColVariable > v_reactive_power_flow;
  // real part is the standard "v_power_flow" variable

  // ----- Generic variables for AC-OPF
  std::vector< ColVariable > v_sum_product_voltages;
  std::vector< ColVariable > v_diff_product_voltages;
  std::vector< ColVariable > v_sqrd_voltages;

  // ----- Generic constraints for AC-OPF
  std::vector< BoxConstraint > v_voltage_bounds_const;
  boost::multi_array< FRowConstraint , 2 > v_angle_bounds_const;
  boost::multi_array< FRowConstraint , 2 > v_voltage_definition_const;
  std::vector< FRowConstraint > v_thermal_limit;
  std::vector< FRowConstraint > v_flow_dc;

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
 }; // end( class( ACNetworkBlock ) )

 /*--------------------------------------------------------------------------*/
 /*--------------------------------------------------------------------------*/
} // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* __ACNetworkBlock */

/*--------------------------------------------------------------------------*/
/*-------------------- End File ACNetworkBlock.h ---------------------------*/
/*--------------------------------------------------------------------------*/
