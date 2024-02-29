/*--------------------------------------------------------------------------*/
/*--------------------- File ThermalUnitBlock.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the ThermalUnitBlock class.
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
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Tiziano Bacci \n
 *         Istituto di Analisi di Sistemi e Informatica "Antonio Ruberti" \n
 *         Consiglio Nazionale delle Ricerche \n
 *
 * \copyright &copy; by Antonio Frangioni, Ali Ghezelsoflu,
 *                      Rafael Durbano Lobato, Donato Meoli, Tiziano Bacci
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "LinearFunction.h"

#include "DQuadFunction.h"

#include "ThermalUnitBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*-------------------------------- CONSTANTS -------------------------------*/
/*--------------------------------------------------------------------------*/

static constexpr unsigned char FormMsk = 7;
// mask for removing all but the first three bits and only leaving the formulation


static constexpr unsigned char tbinForm = 0;
/// the "three binaries" (3bin) formulation is used

static constexpr unsigned char TForm = 1;
/// the T formulation is used

static constexpr unsigned char ptForm = 2;
/// the p_t formulation is used

static constexpr unsigned char DPForm = 3;
/// the "dynamic programming" (DP) formulation is used

static constexpr unsigned char SUForm = 4;
/// the "start-up" (SU) formulation is used

static constexpr unsigned char SDForm = 5;
/// the "shut-down" (SD) formulation is used


static constexpr unsigned char PCuts = 8;
/// 4th bit of AR == 1 if the perspective cuts are used


bool ThermalUnitBlock::f_ignore_netcdf_vars;
/// this variable indicates which netCDF variables must be ignored

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

template< typename T >
static bool identical( std::vector< T > & vec , const Block::Subset sbst ,
                       typename std::vector< T >::const_iterator it ) {
 // returns true if the sub-vector of vec[] corresponding to the indices
 // in sbst is identical to the vector starting at it
 for( auto t : sbst )
  if( vec[ t ] != *(it++) )
   return( false );

 return( true );
}

/*--------------------------------------------------------------------------*/

template< typename T >
static void assign( std::vector< T > & vec , const Block::Subset sbst ,
                    typename std::vector< T >::const_iterator it ) {
 // assign to the sub-vector of vec[] corresponding to the indices in sbst
 // the values found in vector starting at it
 for( auto t : sbst )
  vec[ t ] = *(it++);
}

/*--------------------------------------------------------------------------*/

Block::Subset subset_add( const Block::Subset & sbst , Block::Index dlt ) {
 Block::Subset ret = sbst;
 for( auto & t : ret )
  t += dlt;

 return( ret );
}

/*--------------------------------------------------------------------------*/

Block::Subset subset_sbtrct( const Block::Subset & sbst , Block::Index dlt ) {
 Block::Subset ret = sbst;
 for( auto & t : ret )
  t -= dlt;

 return( ret );
}

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register ThermalUnitBlock to the Block factory

SMSpp_insert_in_factory_cpp_1( ThermalUnitBlock );

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS OF ThermalUnitBlock ----------------------*/
/*--------------------------------------------------------------------------*/

ThermalUnitBlock::~ThermalUnitBlock()
{
 Constraint::clear( CommitmentDesign_Const );
 Constraint::clear( StartUp_ShutDown_Variables_Const );
 Constraint::clear( StartUp_Const );
 Constraint::clear( ShutDown_Const );
 Constraint::clear( RampUp_Const );
 Constraint::clear( RampDown_Const );
 Constraint::clear( MinPower_Const );
 Constraint::clear( MaxPower_Const );
 Constraint::clear( PrimaryRho_Const );
 Constraint::clear( SecondaryRho_Const );

 Constraint::clear( Eq_ActivePower_Const );
 Constraint::clear( Eq_Commitment_Const );
 Constraint::clear( Eq_StartUp_Const );
 Constraint::clear( Eq_ShutDown_Const );
 Constraint::clear( Network_Const );

 Constraint::clear( Init_PC_Const );
 Constraint::clear( Eq_PC_Const );

 Constraint::clear( PC_cuts );

 Constraint::clear( Commitment_bound_Const );
 Constraint::clear( StartUp_Binary_bound_Const );
 Constraint::clear( ShutDown_Binary_bound_Const );

 Constraint::clear( Commitment_fixed_to_One_Const );

 objective.clear();
}

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::deserialize( const netCDF::NcGroup & group )
{

#ifndef NDEBUG
 std::vector< std::string > expected_dims = { "TimeHorizon" ,
                                              "NumberIntervals" };
 check_dimensions( group , expected_dims , std::cerr );

 // we only check for unexpected fields if "this" is a "true"
 // ThermalUnitBlock, i.e., not any derived class. this is because derived
 // classes will likely *have* other fields that tha base class does not
 // know about, and therefore it would complain about them. the idea is that
 // derived classes will then have to check for all expected fields,
 // comprised those of the base class
 // we don't do the same for dimensions as it's unlikely that derived
 // classes will introduce entirely new dimensions
 if( typeid( ThermalUnitBlock ) == typeid( *this ) ) {
  std::vector< std::string > expected_vars = { "InvestmentCost" , "Capacity" ,
                                               "MinPower" , "MaxPower" ,
                                               "DeltaRampUp" , "DeltaRampDown" ,
                                               "PrimaryRho" , "SecondaryRho" ,
                                               "LinearTerm" , "QuadTerm" ,
                                               "ConstTerm" , "StartUpCost" ,
                                               "FixedConsumption" ,
                                               "InertiaCommitment" ,
                                               "InitialPower" , "MinUpTime" ,
                                               "MinDownTime" ,
                                               "InitUpDownTime" ,
                                               "Availability" ,
                                               "StartUpLimit" ,
                                               "ShutDownLimit" };
  check_variables( group , expected_vars , std::cerr );
 }
#endif

 UnitBlock::deserialize( group );

 // Mandatory variables

 ::deserialize( group , "MaxPower" , v_MaxPower , false );

 // Optional variables

 ::deserialize( group , f_InvestmentCost , "InvestmentCost" );

 ::deserialize( group , f_Capacity , "Capacity" );

 if( ::deserialize( group , f_MinUpTime , "MinUpTime" ) )
  f_MinUpTime = std::min( std::max( f_MinUpTime , ( Index ) 1 ) ,
                          f_time_horizon );

 if( ::deserialize( group , f_MinDownTime , "MinDownTime" ) )
  f_MinDownTime = std::min( std::max( f_MinDownTime , ( Index ) 1 ) ,
                            f_time_horizon );

 ::deserialize( group , f_InitialPower , "InitialPower" );

 if( ! ::deserialize( group , f_InitUpDownTime , "InitUpDownTime" ) ) {
  if( f_InitialPower == 0 )
   f_InitUpDownTime = -f_MinDownTime;
  else
   f_InitUpDownTime = f_MinUpTime;
 }

 if( ! ::deserialize( group , "MinPower" , v_MinPower ) )
  v_MinPower.resize( f_time_horizon );

 if( ! ::deserialize( group , "Availability" , v_Availability ) )
  v_Availability.resize( f_time_horizon , 1.0 );

 if( ! ::deserialize( group , "LinearTerm" , v_LinearTerm ) )
  v_LinearTerm.resize( f_time_horizon );

 if( ! ::deserialize( group , "QuadTerm" , v_QuadTerm ) )
  v_QuadTerm.resize( f_time_horizon );

 if( ! ::deserialize( group , "ConstTerm" , v_ConstTerm ) )
  v_ConstTerm.resize( f_time_horizon );

 if( ! ::deserialize( group , "StartUpCost" , v_StartUpCost ) )
  v_StartUpCost.resize( f_time_horizon );

 ::deserialize( group , "DeltaRampUp" , v_DeltaRampUp );

 ::deserialize( group , "DeltaRampDown" , v_DeltaRampDown );

 ::deserialize( group , "FixedConsumption" , v_FixedConsumption );

 ::deserialize( group , "InertiaCommitment" , v_InertiaCommitment );

 if( ! ( f_ignore_netcdf_vars & 1 ) ) {
  ::deserialize( group , "PrimaryRho" , v_PrimaryRho );
  ::deserialize( group , "SecondaryRho" , v_SecondaryRho );
 }

 // Decompress vectors
 decompress_vector( v_MinPower );
 decompress_vector( v_MaxPower );
 decompress_vector( v_Availability );
 decompress_vector( v_DeltaRampUp );
 decompress_vector( v_DeltaRampDown );
 decompress_vector( v_PrimaryRho );
 decompress_vector( v_SecondaryRho );
 decompress_vector( v_LinearTerm );
 decompress_vector( v_QuadTerm );
 decompress_vector( v_ConstTerm );
 decompress_vector( v_StartUpCost );
 decompress_vector( v_FixedConsumption );
 decompress_vector( v_InertiaCommitment );

 if( ! ::deserialize( group , "StartUpLimit" , v_StartUpLimit ) ) {
  v_StartUpLimit.resize( f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   v_StartUpLimit[ t ] = get_operational_min_power( t );
 }

 if( ! ::deserialize( group , "ShutDownLimit" , v_ShutDownLimit ) ) {
  v_ShutDownLimit.resize( f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   v_ShutDownLimit[ t ] = get_operational_min_power( t );
 }

 // Decompress vectors
 decompress_vector( v_StartUpLimit );
 decompress_vector( v_ShutDownLimit );

 check_data_consistency();

}  // end( ThermalUnitBlock::deserialize )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::check_data_consistency( void ) const
{
 // InvestmentCost- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ( f_InvestmentCost != 0 ) && ( f_InitUpDownTime >= 0 ) )
  throw( std::logic_error( "ThermalUnitBlock::check_data_consistency: the "
                           "presence of the investment cost of the thermal "
                           "allows the model to switch into the strategic "
                           "scenario mode, but the presence of a positive "
                           "initial up/down time, typical of the operative "
                           "scenario, is incompatible." ) );

 // Minimum and maximum power - - - - - - - - - - - - - - - - - - - - - - - -
 assert( v_MinPower.size() == f_time_horizon );
 assert( v_MaxPower.size() == f_time_horizon );

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {
  if( v_MinPower[ t ] > v_MaxPower[ t ] )
   throw( std::logic_error( "ThermalUnitBlock::check_data_consistency: "
                            "minimum power at time " + std::to_string( t ) +
                            " is " + std::to_string( v_MinPower[ t ] ) +
                            ", which is greater than the maximum power, which "
                            "is " + std::to_string( v_MaxPower[ t ] ) + "." ) );

  if( v_MinPower[ t ] < 0 )
   throw( std::logic_error( "ThermalUnitBlock::check_data_consistency: "
                            "minimum power for time step "
                            + std::to_string( t ) + " is " +
                            std::to_string( v_MinPower[ t ] ) +
                            ", but it must be nonnegative." ) );
 }

 // Availability- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_Availability.empty() ) {
  assert( v_Availability.size() == f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   if( ( v_Availability[ t ] < 0 ) || ( v_Availability[ t ] > 1 ) )
    throw( std::logic_error( "ThermalUnitBlock::check_data_consistency: "
                             "availability for time step " +
                             std::to_string( t ) + " is " +
                             std::to_string( v_Availability[ t ] ) +
                             ", but it must be between 0 and 1." ) );
 }

 // Delta ramp-up - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_DeltaRampUp.empty() ) {
  assert( v_DeltaRampUp.size() == f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   if( v_DeltaRampUp[ t ] < 0 )
    throw( std::logic_error( "ThermalUnitBlock::check_data_consistency: "
                             "delta ramp-up for time step " +
                             std::to_string( t ) + " is " +
                             std::to_string( v_DeltaRampUp[ t ] ) +
                             ", but it must be nonnegative." ) );
 }

 // Delta ramp-down - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_DeltaRampDown.empty() ) {
  assert( v_DeltaRampDown.size() == f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   if( v_DeltaRampDown[ t ] < 0 )
    throw( std::logic_error( "ThermalUnitBlock::check_data_consistency: "
                             "delta ramp-down for time step " +
                             std::to_string( t ) + " is " +
                             std::to_string( v_DeltaRampDown[ t ] ) +
                             ", but it must be nonnegative." ) );
 }

 // QuadTerm - - - - - - - - - - - - - - - -- - - - - - - - - - - - - - - - -
 if( ! v_QuadTerm.empty() ) {
  assert( v_QuadTerm.size() == f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   if( v_QuadTerm[ t ] < 0 )
    throw( std::logic_error( "ThermalUnitBlock::check_data_consistency: "
                             "quadratic term for time " +
                             std::to_string( t ) + " is " +
                             std::to_string( v_QuadTerm[ t ] ) +
                             ", but it must be nonnegative." ) );
 }

 // InitialPower- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( f_InitialPower < 0 )
  throw( std::logic_error( "ThermalUnitBlock::check_data_consistency: "
                           "initial power is " +
                           std::to_string( f_InitialPower ) +
                           ", but it must be nonnegative." ) );

 // StartUpLimit- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 for( Index t = 0 ; t < f_time_horizon ; ++t )
  if( ( v_StartUpLimit[ t ] < get_operational_min_power( t ) ) ||
      ( v_StartUpLimit[ t ] > get_operational_max_power( t ) ) )
   throw( std::logic_error( "ThermalUnitBlock::check_data_consistency: "
                            "start-up limit for time step " +
                            std::to_string( t ) + " is " +
                            std::to_string( v_StartUpLimit[ t ] ) +
                            ", but it must be between " +
                            std::to_string( get_operational_min_power( t ) )
                            + " and " +
                            std::to_string( get_operational_max_power( t ) )
                            + "." ) );

 // ShutDownLimit - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 for( Index t = 0 ; t < f_time_horizon ; ++t )
  if( ( v_ShutDownLimit[ t ] < get_operational_min_power( t ) ) ||
      ( v_ShutDownLimit[ t ] > get_operational_max_power( t ) ) )
   throw( std::logic_error( "ThermalUnitBlock::check_data_consistency: "
                            "shut-down limit for time step " +
                            std::to_string( t ) + " is " +
                            std::to_string( v_ShutDownLimit[ t ] ) +
                            ", but it must be between " +
                            std::to_string( get_operational_min_power( t ) )
                            + " and " +
                            std::to_string( get_operational_max_power( t ) )
                            + "." ) );

}  // end( ThermalUnitBlock::check_data_consistency )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::generate_abstract_variables( Configuration * stvv )
{
 if( variables_generated() )  // variables have already been generated
  return;                     // nothing to do

 UnitBlock::generate_abstract_variables( stvv );

 Index wf = 1;  // T formulation
 if( ( ! stvv ) && f_BlockConfig )
  stvv = f_BlockConfig->f_static_variables_Configuration;
 if( auto sci = dynamic_cast< SimpleConfiguration< int > * >( stvv ) )
  wf = sci->f_value;

 if( f_InitUpDownTime > 0 )
  init_t = ( f_InitUpDownTime >= f_MinUpTime ? 0 :
             f_MinUpTime - f_InitUpDownTime );
 else
  init_t = ( -f_InitUpDownTime >= f_MinDownTime ? 0 :
             f_MinDownTime + f_InitUpDownTime );

 // Design Binary Variable- - - - - - - - - - - - - - - - - - - - - - - - - -
 if( f_InvestmentCost != 0 ) {
  design.set_type( ColVariable::kBinary );
  add_static_variable( design , "x_thermal" );
 }

 // Commitment Variables- - - - - - - - - - - - - - - - - - - - - - - - - - -
 v_commitment.resize( f_time_horizon );
 for( auto & var : v_commitment )
  var.set_type( ColVariable::kBinary );
 add_static_variable( v_commitment , "u_thermal" );

 // Active Power Variables- - - - - - - - - - - - - - - - - - - - - - - - - -
 v_active_power.resize( f_time_horizon );
 for( auto & var : v_active_power )
  var.set_type( ColVariable::kNonNegative );
 add_static_variable( v_active_power , "p_thermal" );

 // Start-Up and Shut-Down Binary Variables - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 auto startup_shutdown_size = f_time_horizon - init_t;

 if( startup_shutdown_size > 0 ) {

  v_start_up.resize( startup_shutdown_size );
  for( auto & var : v_start_up )
   var.set_type( ColVariable::kBinary );
  add_static_variable( v_start_up , "v_thermal" );

  v_shut_down.resize( startup_shutdown_size );
  for( auto & var : v_shut_down )
   var.set_type( ColVariable::kBinary );
  add_static_variable( v_shut_down , "w_thermal" );
 }

 // Primary Spinning Reserve Variables- - - - - - - - - - - - - - - - - - - -
 if( reserve_vars & 1u )  // if UCBlock has primary demand variables
  if( ! v_PrimaryRho.empty() ) {  // if unit produces any primary reserve
   v_primary_spinning_reserve.resize( f_time_horizon );
   for( auto & var : v_primary_spinning_reserve )
    var.set_type( ColVariable::kNonNegative );
   add_static_variable( v_primary_spinning_reserve , "pr_thermal" );
  }

 // Secondary Spinning Reserve Variables- - - - - - - - - - - - - - - - - - -
 if( reserve_vars & 2u )  // if UCBlock has secondary demand variables
  if( ! v_SecondaryRho.empty() ) {  // if unit produces any secondary reserve
   v_secondary_spinning_reserve.resize( f_time_horizon );
   for( auto & var : v_secondary_spinning_reserve )
    var.set_type( ColVariable::kNonNegative );
   add_static_variable( v_secondary_spinning_reserve , "sc_thermal" );
  }

 // Possibly fixing the commitment variables to 0 or 1- - - - - - - - - - - -
 if( f_InitUpDownTime > 0 ) {

  for( Index t = 0 ; t < init_t ; ++t ) {
   if( ! v_commitment.empty() ) {
    v_commitment[ t ].set_value( 1.0 );
    v_commitment[ t ].is_fixed( true );
   }
  }

  for( Index t = init_t ;
       t < std::min( init_t + f_MinDownTime , f_time_horizon ) ; ++t ) {
   v_start_up[ t - init_t ].set_value( 0.0 );
   v_start_up[ t - init_t ].is_fixed( true );
  }

 } else {

  for( Index t = 0 ; t < init_t ; ++t ) {
   if( ! v_active_power.empty() ) {
    v_active_power[ t ].set_value( 0.0 );
    v_active_power[ t ].is_fixed( true );
   }
   if( ! v_commitment.empty() ) {
    v_commitment[ t ].set_value( 0.0 );
    v_commitment[ t ].is_fixed( true );
   }
   if( ! v_primary_spinning_reserve.empty() ) {
    v_primary_spinning_reserve[ t ].set_value( 0.0 );
    v_primary_spinning_reserve[ t ].is_fixed( true );
   }
   if( ! v_secondary_spinning_reserve.empty() ) {
    v_secondary_spinning_reserve[ t ].set_value( 0.0 );
    v_secondary_spinning_reserve[ t ].is_fixed( true );
   }
  }

  for( Index t = init_t ;
       t < std::min( init_t + f_MinUpTime , f_time_horizon ) ; ++t ) {
   v_shut_down[ t - init_t ].set_value( 0.0 );
   v_shut_down[ t - init_t ].is_fixed( true );
  }
 }

 bool f_cuts = wf & PCuts;

 // Prospective Cuts Variables- - - - - - - - - - - - - - - - - - - - - - - -
 if( f_cuts ) {

  AR |= PCuts;

  v_cut.resize( f_time_horizon );
  for( auto & var : v_cut )
   var.set_type( ColVariable::kNonNegative );
  add_static_variable( v_cut , "z_thermal" );
 }

 // Active Power & Prospective Cuts Variables for DP, SU and SD formulations-
 // (with auxiliary structures initialization)- - - - - - - - - - - - - - - -
 switch( wf & FormMsk ) {

  case( tbinForm ):  // 3bin formulation- - - - - - - - - - - - - - - - - - -

   // AR |= tbinForm;  // does nothing

   break;

  case( TForm ):  // T formulation- - - - - - - - - - - - - - - - - - - - - -

   AR |= TForm;

   break;

  case( ptForm ):  // pt formulation- - - - - - - - - - - - - - - - - - - - -
   // fall through
  case( DPForm ):  // DP formulation- - - - - - - - - - - - - - - - - - - - -
   // fall through
  case( SUForm ):  // SU formulation- - - - - - - - - - - - - - - - - - - - -
   // fall through
  case( SDForm ):  // SD formulation- - - - - - - - - - - - - - - - - - - - -

   if( f_InitUpDownTime > 0 ) {  // if initial committed, OFF_0

    for( Index k = init_t ; k <= f_time_horizon + 1 ; ++k ) {
     if( k == 0 )
      k++;
     v_Y_plus.push_back( std::make_pair( 0 , k ) );
    }

    for( Index k = init_t ; k <= f_time_horizon ; ++k ) {
     if( k <= f_time_horizon - f_MinDownTime - 1 )
      for( Index h = ( k + f_MinDownTime + 1 ) ;  // OFF_h
           h <= f_time_horizon ; ++h )
       v_Y_minus.push_back( std::make_pair( k , h ) );
     v_Y_minus.push_back( std::make_pair( k , f_time_horizon + 1 ) );
     v_nodes_minus.push_back( k );
    }

    for( Index h = ( init_t + f_MinDownTime + 1 ) ;  // OFF_h
         h <= f_time_horizon ; ++h ) {
     if( h <= f_time_horizon - f_MinUpTime + 1 )
      for( Index k = ( h + f_MinUpTime - 1 ) ;  // ON_k
           k <= f_time_horizon ; ++k )
       v_Y_plus.push_back( std::make_pair( h , k ) );
     v_Y_plus.push_back( std::make_pair( h , f_time_horizon + 1 ) );
     v_nodes_plus.push_back( h );
    }

   } else {

    for( Index h = init_t + 1 ; h <= f_time_horizon + 1 ; ++h )  // ON_k
     v_Y_minus.push_back( std::make_pair( 0 , h ) );

    for( Index h = init_t + 1 ; h <= f_time_horizon ; ++h ) {  // OFF_h
     if( h <= f_time_horizon - f_MinUpTime + 1 )
      for( Index k = ( h + f_MinUpTime - 1 ) ;  // ON_k
           k <= f_time_horizon ; ++k )
       v_Y_plus.push_back( std::make_pair( h , k ) );
     v_Y_plus.push_back( std::make_pair( h , f_time_horizon + 1 ) );
     v_nodes_plus.push_back( h );
    }

    for( Index k = init_t + f_MinUpTime ;
         k <= f_time_horizon ; ++k ) {  // ON_k
     if( k <= f_time_horizon - f_MinDownTime - 1 )
      for( Index h = ( k + f_MinDownTime + 1 ) ;  // OFF_h
           h <= f_time_horizon ; ++h )
       v_Y_minus.push_back( std::make_pair( k , h ) );
     v_Y_minus.push_back( std::make_pair( k , f_time_horizon + 1 ) );
     v_nodes_minus.push_back( k );
    }
   }

   v_commitment_plus.resize( v_Y_plus.size() );
   for( auto & var : v_commitment_plus )
    var.set_type( ColVariable::kBinary );
   add_static_variable( v_commitment_plus , "y_plus_thermal" );

   v_commitment_minus.resize( v_Y_minus.size() );
   for( auto & var : v_commitment_minus )
    var.set_type( ColVariable::kBinary );
   add_static_variable( v_commitment_minus , "y_minus_thermal" );

   for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
    for( Index t = 0 ; t < f_time_horizon ; ++t )
     if( ( v_Y_plus[ i ].first <= t + 1 ) &&
         ( t + 1 <= v_Y_plus[ i ].second ) )
      v_P_h_k.push_back(
       std::make_pair( t , std::make_pair( v_Y_plus[ i ].first ,
                                           v_Y_plus[ i ].second ) ) );

   if( ( wf & FormMsk ) == ptForm ) {  // pt formulation- - - - - - - - - - -

    AR |= ptForm;

   } else if( ( wf & FormMsk ) == DPForm ) {  // DP formulation - - - - - - -

    AR |= DPForm;

    v_active_power_h_k.resize( v_P_h_k.size() );
    for( auto & var : v_active_power_h_k )
     var.set_type( ColVariable::kNonNegative );
    add_static_variable( v_active_power_h_k , "p_h_k_thermal" );

    if( f_cuts ) {

     for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
      for( Index t = 0 ; t < f_time_horizon ; ++t )
       if( ( v_Y_plus[ i ].first <= t + 1 ) &&
           ( t + 1 <= v_Y_plus[ i ].second ) )
        v_Z_h_k.push_back(
         std::make_pair( t , std::make_pair( v_Y_plus[ i ].first ,
                                             v_Y_plus[ i ].second ) ) );

     v_cut_h_k.resize( v_Z_h_k.size() );
     for( auto & var : v_cut_h_k )
      var.set_type( ColVariable::kNonNegative );
     add_static_variable( v_cut_h_k , "z_h_k_thermal" );
    }

   } else if( ( wf & FormMsk ) == SUForm ) {  // SU formulation - - - - - - -

    AR |= SUForm;

    bool check_var;
    for( Index i = 0 ; i < v_Y_plus.size() ; ++i ) {
     check_var = true;
     if( i > 0 )
      for( Index j = 0 ; j < v_P_h.size() ; ++j )
       if( v_P_h[ j ].second == v_Y_plus[ i ].first )
        check_var = false;
     if( check_var )
      for( Index t = 0 ; t < f_time_horizon ; ++t )
       if( v_Y_plus[ i ].first <= t + 1 ) {
        v_P_h.push_back( std::make_pair( t , v_Y_plus[ i ].first ) );

        if( f_cuts )
         v_Z_h.push_back( std::make_pair( t , v_Y_plus[ i ].first ) );
       }
    }

    v_active_power_h.resize( v_P_h.size() );
    for( auto & var : v_active_power_h )
     var.set_type( ColVariable::kNonNegative );
    add_static_variable( v_active_power_h , "p_h_thermal" );

    if( f_cuts ) {

     v_cut_h.resize( v_Z_h.size() );
     for( auto & var : v_cut_h )
      var.set_type( ColVariable::kNonNegative );
     add_static_variable( v_cut_h , "z_h_thermal" );
    }

   } else if( ( wf & FormMsk ) == SDForm ) {  // SD formulation - - - - - - -

    AR |= SDForm;

    bool check_var;
    for( Index i = 0 ; i < v_Y_plus.size() ; ++i ) {
     check_var = true;
     if( i > 0 )
      for( Index j = 0 ; j < v_P_k.size() ; ++j )
       if( v_P_k[ j ].second == v_Y_plus[ i ].second )
        check_var = false;
     if( check_var )
      for( Index t = 0 ; t < f_time_horizon ; ++t )
       if( v_Y_plus[ i ].second >= t + 1 ) {
        v_P_k.push_back( std::make_pair( t , v_Y_plus[ i ].second ) );

        if( f_cuts )
         v_Z_k.push_back( std::make_pair( t , v_Y_plus[ i ].second ) );
       }
    }

    v_active_power_k.resize( v_P_k.size() );
    for( auto & var : v_active_power_k )
     var.set_type( ColVariable::kNonNegative );
    add_static_variable( v_active_power_k , "p_k_thermal" );

    if( f_cuts ) {

     v_cut_k.resize( v_Z_k.size() );
     for( auto & var : v_cut_k )
      var.set_type( ColVariable::kNonNegative );
     add_static_variable( v_cut_k , "z_k_thermal" );
    }
   }

   break;

  default:

   throw( std::invalid_argument(
    "ThermalUnitBlock::generate_abstract_variables: invalid formulation" ) );

 }  // end( switch )

 set_variables_generated();

}  // end( ThermalUnitBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::generate_abstract_constraints( Configuration * stcc )
{
 if( constraints_generated() )  // constraints have already been generated
  return;                       // nothing to do

 bool generate_ZOConstraints = false;
 if( ( ! stcc ) && f_BlockConfig )
  stcc = f_BlockConfig->f_static_constraints_Configuration;
 if( auto sci = dynamic_cast< SimpleConfiguration< int > * >( stcc ) )
  generate_ZOConstraints = sci->f_value;

 LinearFunction::v_coeff_pair vars;

 // Initializing commitment design binary variable constraints- - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( f_InvestmentCost != 0 ) {

  CommitmentDesign_Const.resize( f_time_horizon );

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {

   vars.push_back( std::make_pair( &v_commitment[ t ] , 1.0 ) );
   vars.push_back( std::make_pair( &design , -1.0 ) );

   CommitmentDesign_Const[ t ].set_lhs( -Inf< double >() );
   CommitmentDesign_Const[ t ].set_rhs( 0.0 );
   CommitmentDesign_Const[ t ].set_function(
    new LinearFunction( std::move( vars ) ) );
  }

  add_static_constraint( CommitmentDesign_Const ,
                         "CommitmentDesign_Const_Thermal" );
 }

 switch( AR & FormMsk ) {

  case( tbinForm ):  // 3bin formulation- - - - - - - - - - - - - - - - - - -
   // fall through
  case( TForm ): {  // T formulation- - - - - - - - - - - - - - - - - - - - -

   // Initializing start-up and shut-down variables connection constraints- -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   auto startup_shutdown_const_size = f_time_horizon - init_t;

   if( startup_shutdown_const_size > 0 ) {

    StartUp_ShutDown_Variables_Const.resize( startup_shutdown_const_size );

    for( Index t = init_t , cnstr_idx = 0 ; t < f_time_horizon ;
         ++t , ++cnstr_idx ) {

     vars.push_back( std::make_pair( &v_commitment[ t ] , 1.0 ) );
     vars.push_back( std::make_pair( &v_start_up[ t - init_t ] , -1.0 ) );
     vars.push_back( std::make_pair( &v_shut_down[ t - init_t ] , 1.0 ) );

     if( t > init_t )
      vars.push_back( std::make_pair( &v_commitment[ t - 1 ] , -1.0 ) );

     if( ( t == init_t ) && ( f_InitUpDownTime > 0 ) )
      StartUp_ShutDown_Variables_Const[ cnstr_idx ].set_both( 1.0 );
     else
      StartUp_ShutDown_Variables_Const[ cnstr_idx ].set_both( 0.0 );
     StartUp_ShutDown_Variables_Const[ cnstr_idx ].set_function(
      new LinearFunction( std::move( vars ) ) );
    }

    add_static_constraint( StartUp_ShutDown_Variables_Const ,
                           "StartUp_ShutDown_Variables_Const_Thermal" );
   }

   // Initializing turn on constraints (start-up constraints) - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   auto startup_const_size =
    static_cast< int >( f_time_horizon - ( init_t + f_MinUpTime - 1 ) );

   if( startup_const_size > 0 ) {

    StartUp_Const.resize( startup_const_size );

    for( Index t = ( init_t + f_MinUpTime - 1 ) , cnstr_idx = 0 ;
         t < f_time_horizon ; ++t , ++cnstr_idx ) {

     for( Index s = t - ( init_t + f_MinUpTime - 1 ) ; s <= t - init_t ; ++s )
      vars.push_back( std::make_pair( &v_start_up[ s ] , -1.0 ) );

     vars.push_back( std::make_pair( &v_commitment[ t ] , 1.0 ) );

     StartUp_Const[ cnstr_idx ].set_lhs( 0.0 );
     StartUp_Const[ cnstr_idx ].set_rhs( Inf< double >() );
     StartUp_Const[ cnstr_idx ].set_function(
      new LinearFunction( std::move( vars ) ) );
    }

    add_static_constraint( StartUp_Const , "StartUp_Commitment_Const_Thermal" );
   }

   // Initializing turn off constraints (shut-down constraints) - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   auto shutdown_const_size =
    static_cast< int >( f_time_horizon - ( init_t + f_MinDownTime - 1 ) );

   if( shutdown_const_size > 0 ) {

    ShutDown_Const.resize( shutdown_const_size );

    for( Index t = ( init_t + f_MinDownTime - 1 ) , cnstr_idx = 0 ;
         t < f_time_horizon ; ++t , ++cnstr_idx ) {

     for( Index s = t - ( init_t + f_MinDownTime - 1 ) ; s <= t - init_t ; ++s )
      vars.push_back( std::make_pair( &v_shut_down[ s ] , 1.0 ) );

     vars.push_back( std::make_pair( &v_commitment[ t ] , 1.0 ) );

     ShutDown_Const[ cnstr_idx ].set_lhs( -Inf< double >() );
     ShutDown_Const[ cnstr_idx ].set_rhs( 1.0 );
     ShutDown_Const[ cnstr_idx ].set_function(
      new LinearFunction( std::move( vars ) ) );
    }

    add_static_constraint( ShutDown_Const ,
                           "ShutDown_Commitment_Const_Thermal" );
   }

   break;
  }

  case( ptForm ):  // pt formulation- - - - - - - - - - - - - - - - - - - - -
   // fall through
  case( DPForm ):  // DP formulation- - - - - - - - - - - - - - - - - - - - -
   // fall through
  case( SUForm ):  // SU formulation- - - - - - - - - - - - - - - - - - - - -
   // fall through
  case( SDForm ): {  // SD formulation- - - - - - - - - - - - - - - - - - - -

   // Constraints connecting power variables of 3bin with those of DP, SU and
   // SD formulations - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if( ( AR & FormMsk ) == ptForm ) {  // pt formulation- - - - - - - - - - -
    ;  // does nothing
   } else {

    Eq_ActivePower_Const.resize( f_time_horizon );

    if( ( AR & FormMsk ) == DPForm ) {  // DP formulation - - - - - - - - - -

     for( Index t = 0 ; t < f_time_horizon ; ++t ) {

      vars.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );

      for( Index j = 0 ; j < v_P_h_k.size() ; ++j )
       if( v_P_h_k[ j ].first == t )
        vars.push_back( std::make_pair( &v_active_power_h_k[ j ] , -1.0 ) );

      Eq_ActivePower_Const[ t ].set_both( 0.0 );
      Eq_ActivePower_Const[ t ].set_function(
       new LinearFunction( std::move( vars ) ) );
     }

    } else if( ( AR & FormMsk ) == SUForm ) {  // SU formulation- - - - - - -

     for( Index t = 0 ; t < f_time_horizon ; ++t ) {

      vars.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );

      for( Index j = 0 ; j < v_P_h.size() ; ++j )
       if( v_P_h[ j ].first == t )
        vars.push_back( std::make_pair( &v_active_power_h[ j ] , -1.0 ) );

      Eq_ActivePower_Const[ t ].set_both( 0.0 );
      Eq_ActivePower_Const[ t ].set_function(
       new LinearFunction( std::move( vars ) ) );
     }

    } else if( ( AR & FormMsk ) == SDForm ) {  // SD formulation- - - - - - -

     for( Index t = 0 ; t < f_time_horizon ; ++t ) {

      vars.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );

      for( Index j = 0 ; j < v_P_k.size() ; ++j )
       if( v_P_k[ j ].first == t )
        vars.push_back( std::make_pair( &v_active_power_k[ j ] , -1.0 ) );

      Eq_ActivePower_Const[ t ].set_both( 0.0 );
      Eq_ActivePower_Const[ t ].set_function(
       new LinearFunction( std::move( vars ) ) );
     }
    }

    add_static_constraint( Eq_ActivePower_Const ,
                           "Eq_ActivePower_Const_Thermal" );
   }

   // Constraints connecting commitment variables of 3bin with those of pt,
   // DP, SU and SD formulations- - - - - - - - - - - - - - - - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Eq_Commitment_Const.resize( f_time_horizon );

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {

    vars.push_back( std::make_pair( &v_commitment[ t ] , 1.0 ) );

    for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
     if( ( v_Y_plus[ i ].first <= t + 1 ) && ( t + 1 <= v_Y_plus[ i ].second ) )
      vars.push_back( std::make_pair( &v_commitment_plus[ i ] , -1.0 ) );

    Eq_Commitment_Const[ t ].set_both( 0.0 );
    Eq_Commitment_Const[ t ].set_function(
     new LinearFunction( std::move( vars ) ) );
   }

   add_static_constraint( Eq_Commitment_Const , "Eq_Commitment_Const_Thermal" );

   // Constraints connecting start-up variables of 3bin with those of pt,
   // DP, SU and SD formulations- - - - - - - - - - - - - - - - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Eq_StartUp_Const.resize( f_time_horizon - init_t );

   for( Index t = init_t , cnstr_idx = 0 ; t < f_time_horizon ;
        ++t , ++cnstr_idx ) {

    vars.push_back( std::make_pair( &v_start_up[ t - init_t ] , 1.0 ) );

    for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
     if( ( v_Y_plus[ i ].first == t + 1 ) && ( t + 1 <= v_Y_plus[ i ].second ) )
      vars.push_back( std::make_pair( &v_commitment_plus[ i ] , -1.0 ) );

    Eq_StartUp_Const[ cnstr_idx ].set_both( 0.0 );
    Eq_StartUp_Const[ cnstr_idx ].set_function(
     new LinearFunction( std::move( vars ) ) );
   }

   add_static_constraint( Eq_StartUp_Const , "Eq_StartUp_Const_Thermal" );

   // Constraints connecting shut-down variables of 3bin with those of pt,
   // DP, SU and SD formulations- - - - - - - - - - - - - - - - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Eq_ShutDown_Const.resize( f_time_horizon - init_t );

   for( Index t = init_t , cnstr_idx = 0 ; t < f_time_horizon ;
        ++t , ++cnstr_idx ) {

    vars.push_back( std::make_pair( &v_shut_down[ t - init_t ] , 1.0 ) );

    for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
     if( ( v_Y_plus[ i ].first <= t ) && ( t == v_Y_plus[ i ].second ) )
      vars.push_back( std::make_pair( &v_commitment_plus[ i ] , -1.0 ) );

    Eq_ShutDown_Const[ cnstr_idx ].set_both( 0.0 );
    Eq_ShutDown_Const[ cnstr_idx ].set_function(
     new LinearFunction( std::move( vars ) ) );
   }

   add_static_constraint( Eq_ShutDown_Const , "Eq_ShutDown_Const_Thermal" );

   // Network Constraints - - - - - - - - - - - - - - - - - - - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Network_Const.resize( v_nodes_plus.size() + v_nodes_minus.size() + 2 );

   auto cnstr_idx = 0;

   for( Index t = 0 ; t < 1 ; ++t , ++cnstr_idx ) {

    if( f_InitUpDownTime > 0 )
     for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
      if( v_Y_plus[ i ].first == t )
       vars.push_back( std::make_pair( &v_commitment_plus[ i ] , -1.0 ) );

    if( f_InitUpDownTime <= 0 )
     for( Index i = 0 ; i < v_Y_minus.size() ; ++i )
      if( v_Y_minus[ i ].first == t )
       vars.push_back( std::make_pair( &v_commitment_minus[ i ] , -1.0 ) );

    Network_Const[ cnstr_idx ].set_both( -1.0 );
    Network_Const[ cnstr_idx ].set_function(
     new LinearFunction( std::move( vars ) ) );
   }

   for( Index t = 0 ; t < v_nodes_plus.size() ; ++t , ++cnstr_idx ) {

    for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
     if( v_Y_plus[ i ].first == v_nodes_plus[ t ] )
      vars.push_back( std::make_pair( &v_commitment_plus[ i ] , -1.0 ) );

    for( Index i = 0 ; i < v_Y_minus.size() ; ++i )
     if( v_Y_minus[ i ].second == v_nodes_plus[ t ] )
      vars.push_back( std::make_pair( &v_commitment_minus[ i ] , 1.0 ) );

    Network_Const[ cnstr_idx ].set_both( 0.0 );
    Network_Const[ cnstr_idx ].set_function(
     new LinearFunction( std::move( vars ) ) );
   }

   for( Index t = 0 ; t < v_nodes_minus.size() ; ++t , ++cnstr_idx ) {

    for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
     if( v_Y_plus[ i ].second == v_nodes_minus[ t ] )
      vars.push_back( std::make_pair( &v_commitment_plus[ i ] , 1.0 ) );

    for( Index i = 0 ; i < v_Y_minus.size() ; ++i )
     if( v_Y_minus[ i ].first == v_nodes_minus[ t ] )
      vars.push_back( std::make_pair( &v_commitment_minus[ i ] , -1.0 ) );

    Network_Const[ cnstr_idx ].set_both( 0.0 );
    Network_Const[ cnstr_idx ].set_function(
     new LinearFunction( std::move( vars ) ) );
   }

   for( Index t = f_time_horizon + 1 ;
        t < f_time_horizon + 2 ; ++t , ++cnstr_idx ) {

    for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
     if( v_Y_plus[ i ].second == t )
      vars.push_back( std::make_pair( &v_commitment_plus[ i ] , 1.0 ) );

    for( Index i = 0 ; i < v_Y_minus.size() ; ++i )
     if( v_Y_minus[ i ].second == t )
      vars.push_back( std::make_pair( &v_commitment_minus[ i ] , 1.0 ) );

    Network_Const[ cnstr_idx ].set_both( 1.0 );
    Network_Const[ cnstr_idx ].set_function(
     new LinearFunction( std::move( vars ) ) );
   }

   add_static_constraint( Network_Const , "Network_Const_Thermal" );

   break;
  }

  default:

   exit( 1 );

 }  // end( switch )

 // Initializing ramp-up constraints- - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! v_DeltaRampUp.empty() )
  if( f_InitUpDownTime > 0 )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    if( f_InitialPower + v_DeltaRampUp[ 0 ] < get_operational_min_power( 0 ) )
     throw( std::logic_error(
      "ThermalUnitBlock::RampUpConstraints: when f_InitUpDownTime > 0, "
      "it must be that f_InitialPower + v_DeltaRampUp[ 0 ] >= "
      "get_operational_min_power( 0 )." ) );

 if( ! v_DeltaRampDown.empty() )
  if( f_InitUpDownTime > 0 )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    if( f_InitialPower - v_DeltaRampDown[ 0 ] > get_operational_max_power( 0 ) )
     throw( std::logic_error(
      "ThermalUnitBlock::RampDownConstraints: when f_InitUpDownTime > 0, "
      "it must be that f_InitialPower - v_DeltaRampDown[ 0 ] <= "
      "get_operational_max_power( 0 )." ) );

 if( ! v_DeltaRampUp.empty() ) {

  if( ( AR & FormMsk ) == tbinForm ) {  // 3bin formulation - - - - - - - - -

   RampUp_Const.resize( f_time_horizon );

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {

    vars.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );

    if( t > 0 ) {
     vars.push_back( std::make_pair( &v_active_power[ t - 1 ] , -1.0 ) );
     vars.push_back( std::make_pair( &v_commitment[ t - 1 ] ,
                                     -v_DeltaRampUp[ t - 1 ] ) );
    }

    if( t >= init_t )
     vars.push_back( std::make_pair( &v_start_up[ t - init_t ] ,
                                     -v_StartUpLimit[ t ] ) );

    RampUp_Const[ t ].set_lhs( -Inf< double >() );
    if( ( t == 0 ) && ( f_InitUpDownTime > 0 ) )
     RampUp_Const[ t ].set_rhs( f_InitialPower + v_DeltaRampUp[ t ] );
    else if( ( t > 0 ) || ( ( t == 0 ) && ( f_InitUpDownTime <= 0 ) ) )
     RampUp_Const[ t ].set_rhs( 0.0 );
    RampUp_Const[ t ].set_function( new LinearFunction( std::move( vars ) ) );
   }

  } else if( ( AR & FormMsk ) == TForm ) {  // T formulation- - - - - - - - -

   RampUp_Const.resize( f_time_horizon );

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {

    vars.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );

    if( t > 0 ) {
     vars.push_back( std::make_pair( &v_active_power[ t - 1 ] , -1.0 ) );
     vars.push_back( std::make_pair( &v_commitment[ t - 1 ] ,
                                     get_operational_min_power( t - 1 ) ) );
    }

    if( t >= init_t )
     vars.push_back( std::make_pair( &v_start_up[ t - init_t ] ,
                                     -( v_StartUpLimit[ t ] -
                                        get_operational_min_power( t ) -
                                        v_DeltaRampUp[ t ] ) ) );

    vars.push_back( std::make_pair( &v_commitment[ t ] ,
                                    -( v_DeltaRampUp[ t ] +
                                       get_operational_min_power( t ) ) ) );

    RampUp_Const[ t ].set_lhs( -Inf< double >() );
    if( ( t == 0 ) && ( f_InitUpDownTime > 0 ) )
     RampUp_Const[ t ].set_rhs(
      f_InitialPower - get_operational_min_power( t ) );
    else if( ( t > 0 ) || ( ( t == 0 ) && ( f_InitUpDownTime <= 0 ) ) )
     RampUp_Const[ t ].set_rhs( 0.0 );
    RampUp_Const[ t ].set_function( new LinearFunction( std::move( vars ) ) );
   }

  } else if( ( AR & FormMsk ) == ptForm ) {  // pt formulation- - - - - - - -

   if( f_InitUpDownTime > 0 )
    RampUp_Const.resize( f_time_horizon );
   if( f_InitUpDownTime <= 0 )
    RampUp_Const.resize( f_time_horizon - 1 );

   auto cnstr_idx = 0;

   for( Index t = 0 ; t < f_time_horizon ; ++t )
    if( ( t > 0 ) || ( ( t == 0 ) && ( f_InitUpDownTime > 0 ) ) ) {

     vars.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );

     if( t > 0 ) {
      vars.push_back( std::make_pair( &v_active_power[ t - 1 ] , -1.0 ) );
      for( Index i = 0 ; i < v_Y_plus.size() ; ++i ) {
       if( ( v_Y_plus[ i ].first <= t ) && ( t < v_Y_plus[ i ].second ) )
        vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                        -v_DeltaRampUp[ t - 1 ] ) );
       if( ( v_Y_plus[ i ].first <= t ) && ( t == v_Y_plus[ i ].second ) )
        vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                        get_operational_min_power( t - 1 ) ) );
       if( ( v_Y_plus[ i ].first == t + 1 ) &&
           ( t + 1 <= v_Y_plus[ i ].second ) )
        vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                        -v_StartUpLimit[ t ] ) );
      }
     }

     if( ( t == 0 ) && ( f_InitUpDownTime > 0 ) )
      for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
       if( ( v_Y_plus[ i ].first == 0 ) && ( t + 1 <= v_Y_plus[ i ].second ) )
        vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                        -f_InitialPower -
                                        v_DeltaRampUp[ t ] ) );

     RampUp_Const[ cnstr_idx ].set_lhs( -Inf< double >() );
     RampUp_Const[ cnstr_idx ].set_rhs( 0.0 );
     RampUp_Const[ cnstr_idx ].set_function(
      new LinearFunction( std::move( vars ) ) );

     cnstr_idx++;
    }

  } else if( ( AR & FormMsk ) == DPForm ) {  // DP formulation- - - - - - - -

   auto ramp_up_cnstrs_size = 0;

   for( Index j = 0 ; j < v_P_h_k.size() ; ++j ) {
    auto t = v_P_h_k[ j ].first;
    if( ( ( t == 0 ) && ( f_InitUpDownTime > 0 ) ) ||
        ( ( t > 0 ) && ( v_P_h_k[ j ].second.first + 1 <= t + 1 ) &&
          ( v_P_h_k[ j ].second.second >= t + 1 ) ) )
     ramp_up_cnstrs_size++;
   }

   if( ramp_up_cnstrs_size > 0 ) {

    RampUp_Const.resize( ramp_up_cnstrs_size );

    auto cnstr_idx = 0;

    for( Index j = 0 ; j < v_P_h_k.size() ; ++j ) {

     auto t = v_P_h_k[ j ].first;
     if( ( ( t == 0 ) && ( f_InitUpDownTime > 0 ) ) ||
         ( ( t > 0 ) && ( v_P_h_k[ j ].second.first + 1 <= t + 1 ) &&
           ( v_P_h_k[ j ].second.second >= t + 1 ) ) ) {

      vars.push_back( std::make_pair( &v_active_power_h_k[ j ] , 1.0 ) );

      if( t > 0 )
       for( Index s = 0 ; s < v_P_h_k.size() ; ++s )
        if( ( v_P_h_k[ j ].second.first == v_P_h_k[ s ].second.first ) &&
            ( v_P_h_k[ j ].second.second == v_P_h_k[ s ].second.second ) )
         if( v_P_h_k[ s ].first == t - 1 )
          vars.push_back( std::make_pair( &v_active_power_h_k[ s ] , -1.0 ) );

      for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
       if( ( v_P_h_k[ j ].second.first == v_Y_plus[ i ].first ) &&
           ( v_P_h_k[ j ].second.second == v_Y_plus[ i ].second ) ) {
        if( t == 0 )
         vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                         -v_DeltaRampUp[ t ] -
                                         f_InitialPower ) );
        else if( t > 0 )
         vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                         -v_DeltaRampUp[ t ] ) );
       }

      RampUp_Const[ cnstr_idx ].set_lhs( -Inf< double >() );
      RampUp_Const[ cnstr_idx ].set_rhs( 0.0 );
      RampUp_Const[ cnstr_idx ].set_function(
       new LinearFunction( std::move( vars ) ) );

      cnstr_idx++;
     }
    }
   }

  } else if( ( AR & FormMsk ) == SUForm ) {  // SU formulation- - - - - - - -

   auto ramp_up_cnstrs_size = 0;

   for( Index t = 0 ; t < f_time_horizon ; ++t )
    for( Index j = 0 ; j < v_P_h.size() ; ++j )
     if( v_P_h[ j ].first == t ) {
      if( ( f_InitUpDownTime > 0 ) && ( t == 0 ) )
       ramp_up_cnstrs_size++;
      if( ( t > 0 ) && ( t + 1 > v_P_h[ j ].second ) )
       ramp_up_cnstrs_size++;
     }

   RampUp_Const.resize( ramp_up_cnstrs_size );

   auto cnstr_idx = 0;

   for( Index t = 0 ; t < f_time_horizon ; ++t )
    for( Index j = 0 ; j < v_P_h.size() ; ++j )
     if( v_P_h[ j ].first == t ) {
      if( ( f_InitUpDownTime > 0 ) && ( t == 0 ) ) {

       vars.push_back( std::make_pair( &v_active_power_h[ j ] , 1.0 ) );

       for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
        if( v_P_h[ j ].second == v_Y_plus[ i ].first )
         if( t + 1 <= v_Y_plus[ i ].second )
          vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                          -v_DeltaRampUp[ t ] -
                                          f_InitialPower ) );

       RampUp_Const[ cnstr_idx ].set_lhs( -Inf< double >() );
       RampUp_Const[ cnstr_idx ].set_rhs( 0.0 );
       RampUp_Const[ cnstr_idx ].set_function(
        new LinearFunction( std::move( vars ) ) );

       cnstr_idx++;
      }
      if( ( t > 0 ) && ( t + 1 > v_P_h[ j ].second ) ) {

       vars.push_back( std::make_pair( &v_active_power_h[ j ] , 1.0 ) );

       for( Index s = 0 ; s < v_P_h.size() ; ++s )
        if( ( v_P_h[ j ].first - 1 == v_P_h[ s ].first ) &&
            ( v_P_h[ j ].second == v_P_h[ s ].second ) )
         vars.push_back( std::make_pair( &v_active_power_h[ s ] , -1.0 ) );
       for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
        if( v_P_h[ j ].second == v_Y_plus[ i ].first ) {
         if( t + 1 <= v_Y_plus[ i ].second )
          vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                          -v_DeltaRampUp[ t ] ) );
         if( t == v_Y_plus[ i ].second )
          vars.push_back( std::make_pair(
           &v_commitment_plus[ i ] ,
           get_operational_min_power( t - 1 ) ) );
        }

       RampUp_Const[ cnstr_idx ].set_lhs( -Inf< double >() );
       RampUp_Const[ cnstr_idx ].set_rhs( 0.0 );
       RampUp_Const[ cnstr_idx ].set_function(
        new LinearFunction( std::move( vars ) ) );

       cnstr_idx++;
      }
     }

  } else if( ( AR & FormMsk ) == SDForm ) {  // SD formulation- - - - - - - -

   auto ramp_up_cnstrs_size = 0;

   for( Index t = 0 ; t < f_time_horizon ; ++t )
    for( Index j = 0 ; j < v_P_k.size() ; ++j )
     if( v_P_k[ j ].first == t ) {
      if( ( f_InitUpDownTime > 0 ) && ( t == 0 ) )
       ramp_up_cnstrs_size++;
      if( ( t > 0 ) && ( t + 1 <= v_P_k[ j ].second ) )
       ramp_up_cnstrs_size++;
     }

   RampUp_Const.resize( ramp_up_cnstrs_size );

   auto cnstr_idx = 0;

   for( Index t = 0 ; t < f_time_horizon ; ++t )
    for( Index j = 0 ; j < v_P_k.size() ; ++j )
     if( v_P_k[ j ].first == t ) {
      if( ( f_InitUpDownTime > 0 ) && ( t == 0 ) ) {

       vars.push_back( std::make_pair( &v_active_power_k[ j ] , 1.0 ) );

       for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
        if( v_P_k[ j ].second == v_Y_plus[ i ].second )
         if( v_Y_plus[ i ].first <= t )
          vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                          -v_DeltaRampUp[ t ] -
                                          f_InitialPower ) );

       RampUp_Const[ cnstr_idx ].set_lhs( -Inf< double >() );
       RampUp_Const[ cnstr_idx ].set_rhs( 0.0 );
       RampUp_Const[ cnstr_idx ].set_function(
        new LinearFunction( std::move( vars ) ) );

       cnstr_idx++;
      }
      if( ( t > 0 ) && ( t + 1 <= v_P_k[ j ].second ) ) {

       vars.push_back( std::make_pair( &v_active_power_k[ j ] , 1.0 ) );

       for( Index s = 0 ; s < v_P_k.size() ; ++s )
        if( ( v_P_k[ j ].first - 1 == v_P_k[ s ].first ) &&
            ( v_P_k[ j ].second == v_P_k[ s ].second ) )
         vars.push_back( std::make_pair( &v_active_power_k[ s ] , -1.0 ) );
       for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
        if( v_P_k[ j ].second == v_Y_plus[ i ].second ) {
         if( v_Y_plus[ i ].first <= t )
          vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                          -v_DeltaRampUp[ t ] ) );
         if( t + 1 == v_Y_plus[ i ].first )
          vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                          -v_StartUpLimit[ t ] ) );
        }

       RampUp_Const[ cnstr_idx ].set_lhs( -Inf< double >() );
       RampUp_Const[ cnstr_idx ].set_rhs( 0.0 );
       RampUp_Const[ cnstr_idx ].set_function(
        new LinearFunction( std::move( vars ) ) );

       cnstr_idx++;
      }
     }
  }

  add_static_constraint( RampUp_Const , "RampUp_Const_Thermal" );
 }

 // Initializing ramp-down constraints- - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! v_DeltaRampDown.empty() ) {

  if( ( AR & FormMsk ) == tbinForm ) {  // 3bin formulation - - - - - - - - -

   RampDown_Const.resize( f_time_horizon );

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {

    vars.push_back( std::make_pair( &v_active_power[ t ] , -1.0 ) );
    vars.push_back( std::make_pair( &v_commitment[ t ] ,
                                    -v_DeltaRampDown[ t ] ) );

    if( t > 0 )
     vars.push_back( std::make_pair( &v_active_power[ t - 1 ] , 1.0 ) );

    if( t >= init_t )
     vars.push_back( std::make_pair( &v_shut_down[ t - init_t ] ,
                                     -v_ShutDownLimit[ t ] ) );

    RampDown_Const[ t ].set_lhs( -Inf< double >() );
    if( ( t == 0 ) && ( f_InitUpDownTime > 0 ) )
     RampDown_Const[ t ].set_rhs( -f_InitialPower );
    else if( ( t > 0 ) || ( ( t == 0 ) && ( f_InitUpDownTime <= 0 ) ) )
     RampDown_Const[ t ].set_rhs( 0.0 );
    RampDown_Const[ t ].set_function( new LinearFunction( std::move( vars ) ) );
   }

  } else if( ( AR & FormMsk ) == TForm ) {  // T formulation- - - - - - - - -

   RampDown_Const.resize( f_time_horizon );

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {

    vars.push_back( std::make_pair( &v_active_power[ t ] , -1.0 ) );

    if( t > 0 ) {
     vars.push_back( std::make_pair( &v_active_power[ t - 1 ] , 1.0 ) );
     vars.push_back( std::make_pair(
      &v_commitment[ t - 1 ] ,
      -( v_DeltaRampDown[ t - 1 ] +
         get_operational_min_power( t - 1 ) ) ) );
    }

    if( t >= init_t )
     vars.push_back( std::make_pair( &v_shut_down[ t - init_t ] ,
                                     -( v_ShutDownLimit[ t ] -
                                        get_operational_min_power( t ) -
                                        v_DeltaRampDown[ t ] ) ) );

    vars.push_back( std::make_pair( &v_commitment[ t ] ,
                                    get_operational_min_power( t ) ) );

    RampDown_Const[ t ].set_lhs( -Inf< double >() );
    if( ( t == 0 ) && ( f_InitUpDownTime > 0 ) )
     RampDown_Const[ t ].set_rhs(
      -( f_InitialPower - v_DeltaRampDown[ t ] -
         get_operational_min_power( t ) ) );
    else if( ( t > 0 ) || ( ( t == 0 ) && ( f_InitUpDownTime <= 0 ) ) )
     RampDown_Const[ t ].set_rhs( 0.0 );
    RampDown_Const[ t ].set_function( new LinearFunction( std::move( vars ) ) );
   }

  } else if( ( AR & FormMsk ) == ptForm ) {  // pt formulation- - - - - - - -

   if( f_InitUpDownTime > 0 )
    RampDown_Const.resize( f_time_horizon );
   if( f_InitUpDownTime <= 0 )
    RampDown_Const.resize( f_time_horizon - 1 );

   auto cnstr_idx = 0;

   for( Index t = 0 ; t < f_time_horizon ; ++t )
    if( ( t > 0 ) || ( ( t == 0 ) && ( f_InitUpDownTime > 0 ) ) ) {

     vars.push_back( std::make_pair( &v_active_power[ t ] , -1.0 ) );

     if( t > 0 ) {
      vars.push_back( std::make_pair( &v_active_power[ t - 1 ] , 1.0 ) );
      for( Index i = 0 ; i < v_Y_plus.size() ; ++i ) {
       if( ( v_Y_plus[ i ].first <= t ) && ( t < v_Y_plus[ i ].second ) )
        vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                        -v_DeltaRampDown[ t - 1 ] ) );
       if( ( v_Y_plus[ i ].first <= t ) && ( t == v_Y_plus[ i ].second ) )
        vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                        -v_ShutDownLimit[ t - 1 ] ) );
       if( ( v_Y_plus[ i ].first == t + 1 ) &&
           ( t + 1 <= v_Y_plus[ i ].second ) )
        vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                        get_operational_min_power( t ) ) );
      }
     }

     if( ( t == 0 ) && ( f_InitUpDownTime > 0 ) )
      for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
       if( ( v_Y_plus[ i ].first == 0 ) && ( t + 1 <= v_Y_plus[ i ].second ) )
        vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                        f_InitialPower -
                                        v_DeltaRampDown[ t ] ) );

     RampDown_Const[ cnstr_idx ].set_lhs( -Inf< double >() );
     RampDown_Const[ cnstr_idx ].set_rhs( 0.0 );
     RampDown_Const[ cnstr_idx ].set_function(
      new LinearFunction( std::move( vars ) ) );

     cnstr_idx++;
    }

  } else if( ( AR & FormMsk ) == DPForm ) {  // DP formulation- - - - - - - -

   auto ramp_up_cnstrs_size = 0;

   for( Index j = 0 ; j < v_P_h_k.size() ; ++j ) {
    auto t = v_P_h_k[ j ].first;
    if( ( ( t == 0 ) && ( f_InitUpDownTime > 0 ) ) ||
        ( ( t > 0 ) && ( v_P_h_k[ j ].second.first + 1 <= t + 1 ) &&
          ( v_P_h_k[ j ].second.second >= t + 1 ) ) )
     ramp_up_cnstrs_size++;
   }

   if( ramp_up_cnstrs_size > 0 ) {

    RampDown_Const.resize( ramp_up_cnstrs_size );

    auto cnstr_idx = 0;

    for( Index j = 0 ; j < v_P_h_k.size() ; ++j ) {

     auto t = v_P_h_k[ j ].first;
     if( ( ( t == 0 ) && ( f_InitUpDownTime > 0 ) ) ||
         ( ( t > 0 ) && ( v_P_h_k[ j ].second.first + 1 <= t + 1 ) &&
           ( v_P_h_k[ j ].second.second >= t + 1 ) ) ) {

      vars.push_back( std::make_pair( &v_active_power_h_k[ j ] , -1.0 ) );

      if( t > 0 )
       for( Index s = 0 ; s < v_P_h_k.size() ; ++s )
        if( ( v_P_h_k[ j ].second.first == v_P_h_k[ s ].second.first ) &&
            ( v_P_h_k[ j ].second.second == v_P_h_k[ s ].second.second ) )
         if( v_P_h_k[ s ].first == t - 1 )
          vars.push_back( std::make_pair( &v_active_power_h_k[ s ] , 1.0 ) );

      for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
       if( ( v_P_h_k[ j ].second.first == v_Y_plus[ i ].first ) &&
           ( v_P_h_k[ j ].second.second == v_Y_plus[ i ].second ) ) {
        if( t == 0 )
         vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                         -v_DeltaRampDown[ t ] +
                                         f_InitialPower ) );
        else if( t > 0 )
         vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                         -v_DeltaRampDown[ t ] ) );
       }

      RampDown_Const[ cnstr_idx ].set_lhs( -Inf< double >() );
      RampDown_Const[ cnstr_idx ].set_rhs( 0.0 );
      RampDown_Const[ cnstr_idx ].set_function(
       new LinearFunction( std::move( vars ) ) );

      cnstr_idx++;
     }
    }
   }

  } else if( ( AR & FormMsk ) == SUForm ) {  // SU formulation- - - - - - - -

   auto ramp_down_cnstrs_size = 0;

   for( Index t = 0 ; t < f_time_horizon ; ++t )
    for( Index j = 0 ; j < v_P_h.size() ; ++j )
     if( v_P_h[ j ].first == t ) {
      if( ( f_InitUpDownTime > 0 ) && ( t == 0 ) )
       ramp_down_cnstrs_size++;
      if( ( t > 0 ) && ( t + 1 > v_P_h[ j ].second ) )
       ramp_down_cnstrs_size++;
     }

   RampDown_Const.resize( ramp_down_cnstrs_size );

   auto cnstr_idx = 0;

   for( Index t = 0 ; t < f_time_horizon ; ++t )
    for( Index j = 0 ; j < v_P_h.size() ; ++j )
     if( v_P_h[ j ].first == t ) {
      if( ( f_InitUpDownTime > 0 ) && ( t == 0 ) ) {

       vars.push_back( std::make_pair( &v_active_power_h[ j ] , -1.0 ) );

       for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
        if( v_P_h[ j ].second == v_Y_plus[ i ].first )
         if( t + 1 <= v_Y_plus[ i ].second )
          vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                          -v_DeltaRampDown[ t ] +
                                          f_InitialPower ) );

       RampDown_Const[ cnstr_idx ].set_lhs( -Inf< double >() );
       RampDown_Const[ cnstr_idx ].set_rhs( 0.0 );
       RampDown_Const[ cnstr_idx ].set_function(
        new LinearFunction( std::move( vars ) ) );

       cnstr_idx++;
      }

      if( ( t > 0 ) && ( t + 1 > v_P_h[ j ].second ) ) {

       vars.push_back( std::make_pair( &v_active_power_h[ j ] , -1.0 ) );

       for( Index s = 0 ; s < v_P_h.size() ; ++s )
        if( ( v_P_h[ j ].first - 1 == v_P_h[ s ].first ) &&
            ( v_P_h[ j ].second == v_P_h[ s ].second ) )
         vars.push_back( std::make_pair( &v_active_power_h[ s ] , 1.0 ) );
       for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
        if( v_P_h[ j ].second == v_Y_plus[ i ].first ) {
         if( t + 1 <= v_Y_plus[ i ].second )
          vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                          -v_DeltaRampDown[ t ] ) );
         if( t == v_Y_plus[ i ].second )
          vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                          -v_ShutDownLimit[ t - 1 ] ) );
        }

       RampDown_Const[ cnstr_idx ].set_lhs( -Inf< double >() );
       RampDown_Const[ cnstr_idx ].set_rhs( 0.0 );
       RampDown_Const[ cnstr_idx ].set_function(
        new LinearFunction( std::move( vars ) ) );

       cnstr_idx++;
      }
     }

  } else if( ( AR & FormMsk ) == SDForm ) {  // SD formulation- - - - - - - -

   auto ramp_down_cnstrs_size = 0;

   for( Index t = 0 ; t < f_time_horizon ; ++t )
    for( Index j = 0 ; j < v_P_k.size() ; ++j )
     if( v_P_k[ j ].first == t ) {
      if( ( f_InitUpDownTime > 0 ) && ( t == 0 ) )
       ramp_down_cnstrs_size++;
      if( ( t > 0 ) && ( t + 1 <= v_P_k[ j ].second ) )
       ramp_down_cnstrs_size++;
     }

   RampDown_Const.resize( ramp_down_cnstrs_size );

   auto cnstr_idx = 0;

   for( Index t = 0 ; t < f_time_horizon ; ++t )
    for( Index j = 0 ; j < v_P_k.size() ; ++j )
     if( v_P_k[ j ].first == t ) {
      if( ( f_InitUpDownTime > 0 ) && ( t == 0 ) ) {

       vars.push_back( std::make_pair( &v_active_power_k[ j ] , -1.0 ) );
       for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
        if( v_P_k[ j ].second == v_Y_plus[ i ].second )
         if( v_Y_plus[ i ].first <= t )
          vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                          -v_DeltaRampDown[ t ] +
                                          f_InitialPower ) );

       RampDown_Const[ cnstr_idx ].set_lhs( -Inf< double >() );
       RampDown_Const[ cnstr_idx ].set_rhs( 0.0 );
       RampDown_Const[ cnstr_idx ].set_function(
        new LinearFunction( std::move( vars ) ) );

       cnstr_idx++;
      }
      if( ( t > 0 ) && ( t + 1 <= v_P_k[ j ].second ) ) {

       vars.push_back( std::make_pair( &v_active_power_k[ j ] , -1.0 ) );

       for( Index s = 0 ; s < v_P_k.size() ; ++s )
        if( ( v_P_k[ j ].first - 1 == v_P_k[ s ].first ) &&
            ( v_P_k[ j ].second == v_P_k[ s ].second ) )
         vars.push_back( std::make_pair( &v_active_power_k[ s ] , 1.0 ) );
       for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
        if( v_P_k[ j ].second == v_Y_plus[ i ].second ) {
         if( v_Y_plus[ i ].first <= t )
          vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                          -v_DeltaRampDown[ t ] ) );
         if( t + 1 == v_Y_plus[ i ].first )
          vars.push_back( std::make_pair(
           &v_commitment_plus[ i ] ,
           get_operational_min_power( t - 1 ) ) );
        }

       RampDown_Const[ cnstr_idx ].set_lhs( -Inf< double >() );
       RampDown_Const[ cnstr_idx ].set_rhs( 0.0 );
       RampDown_Const[ cnstr_idx ].set_function(
        new LinearFunction( std::move( vars ) ) );

       cnstr_idx++;
      }
     }
  }

  add_static_constraint( RampDown_Const , "RampDown_Const_Thermal" );
 }

 // Initializing minimum power constraints- - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( ( AR & FormMsk ) == tbinForm ) ||  // 3bin formulation - - - - - - - -
     ( ( AR & FormMsk ) == TForm ) ) {  // T formulation- - - - - - - - - - -

  MinPower_Const.resize( f_time_horizon );

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {

   vars.push_back( std::make_pair( &v_commitment[ t ] ,
                                   -get_operational_min_power( t ) ) );
   vars.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );

   // if UCBlock has primary demand variables
   if( ( reserve_vars & 1u ) && ( ! v_PrimaryRho.empty() ) )
    vars.push_back( std::make_pair( &v_primary_spinning_reserve[ t ] ,
                                    -1.0 ) );

   // if UCBlock has secondary reserve variables
   if( ( reserve_vars & 2u ) && ( ! v_SecondaryRho.empty() ) )
    vars.push_back( std::make_pair( &v_secondary_spinning_reserve[ t ] ,
                                    -1.0 ) );

   MinPower_Const[ t ].set_lhs( 0.0 );
   MinPower_Const[ t ].set_rhs( Inf< double >() );
   MinPower_Const[ t ].set_function( new LinearFunction( std::move( vars ) ) );
  }

 } else if( ( AR & FormMsk ) == ptForm ) {  // pt formulation - - - - - - - -

  MinPower_Const.resize( f_time_horizon );

  for( Index t = 0 , constraint_index = 0 ; t < f_time_horizon ;
       ++t , ++constraint_index ) {

   vars.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );

   // if UCBlock has primary demand variables
   if( ( reserve_vars & 1u ) && ( ! v_PrimaryRho.empty() ) )
    vars.push_back( std::make_pair( &v_primary_spinning_reserve[ t ] ,
                                    -1.0 ) );

   // if UCBlock has secondary reserve variables
   if( ( reserve_vars & 2u ) && ( ! v_SecondaryRho.empty() ) )
    vars.push_back( std::make_pair( &v_secondary_spinning_reserve[ t ] ,
                                    -1.0 ) );

   for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
    if( ( v_Y_plus[ i ].first <= t + 1 ) &&
        ( t + 1 <= v_Y_plus[ i ].second ) )
     vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                     -get_operational_min_power( t ) ) );

   MinPower_Const[ constraint_index ].set_lhs( 0.0 );
   MinPower_Const[ constraint_index ].set_rhs( Inf< double >() );
   MinPower_Const[ constraint_index ].set_function(
    new LinearFunction( std::move( vars ) ) );
  }

 } else if( ( AR & FormMsk ) == DPForm ) {  // DP formulation - - - - - - - -

  MinPower_Const.resize( v_P_h_k.size() );

  for( Index j = 0 ; j < v_P_h_k.size() ; ++j ) {

   auto t = v_P_h_k[ j ].first;
   for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
    if( ( v_P_h_k[ j ].second.first == v_Y_plus[ i ].first ) &&
        ( v_P_h_k[ j ].second.second == v_Y_plus[ i ].second ) )
     vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                     -get_operational_min_power( t ) ) );

   vars.push_back( std::make_pair( &v_active_power_h_k[ j ] , 1.0 ) );

   MinPower_Const[ j ].set_lhs( 0.0 );
   MinPower_Const[ j ].set_rhs( Inf< double >() );
   MinPower_Const[ j ].set_function( new LinearFunction( std::move( vars ) ) );
  }

 } else if( ( AR & FormMsk ) == SUForm ) {  // SU formulation - - - - - - - -

  MinPower_Const.resize( v_P_h.size() );

  for( Index j = 0 ; j < v_P_h.size() ; ++j ) {

   auto t = v_P_h[ j ].first;
   for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
    if( ( v_P_h[ j ].second == v_Y_plus[ i ].first ) &&
        ( t + 1 <= v_Y_plus[ i ].second ) )
     vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                     -get_operational_min_power( t ) ) );

   vars.push_back( std::make_pair( &v_active_power_h[ j ] , 1.0 ) );

   MinPower_Const[ j ].set_lhs( 0.0 );
   MinPower_Const[ j ].set_rhs( Inf< double >() );
   MinPower_Const[ j ].set_function( new LinearFunction( std::move( vars ) ) );
  }

 } else if( ( AR & FormMsk ) == SDForm ) {  // SD formulation - - - - - - - -

  MinPower_Const.resize( v_P_k.size() );

  for( Index j = 0 ; j < v_P_k.size() ; ++j ) {

   auto t = v_P_k[ j ].first;
   for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
    if( ( v_P_k[ j ].second == v_Y_plus[ i ].second ) &&
        ( v_Y_plus[ i ].first <= t + 1 ) )
     vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                     -get_operational_min_power( t ) ) );

   vars.push_back( std::make_pair( &v_active_power_k[ j ] , 1.0 ) );

   MinPower_Const[ j ].set_lhs( 0.0 );
   MinPower_Const[ j ].set_rhs( Inf< double >() );
   MinPower_Const[ j ].set_function( new LinearFunction( std::move( vars ) ) );
  }
 }

 add_static_constraint( MinPower_Const , "MinPower_Const_Thermal" );

 // Initializing maximum power constraints- - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( AR & FormMsk ) == tbinForm ) {  // 3bin formulation- - - - - - - - - -

  MaxPower_Const.resize( f_MinUpTime == 1 ?
                         ( init_t == 0 ?
                           2 * ( f_time_horizon - init_t ) - 2 + init_t :
                           2 * ( f_time_horizon - init_t ) - 1 + init_t ) :
                         f_time_horizon );

  for( Index t = 0 , cnstr_idx = 0 ; t < f_time_horizon ; ++t , ++cnstr_idx ) {

   if( t >= init_t ) {
    if( t == 0 )
     vars.push_back( std::make_pair( &v_shut_down[ t + 1 - init_t ] ,
                                     v_ShutDownLimit[ t ] -
                                     get_operational_max_power( t ) ) );
    if( t == f_time_horizon - 1 )
     vars.push_back( std::make_pair( &v_start_up[ t - init_t ] ,
                                     v_StartUpLimit[ t ] -
                                     get_operational_max_power( t ) ) );

    if( f_MinUpTime == 1 ) {
     if( ( t > 0 ) && ( t < f_time_horizon - 1 ) ) {
      vars.push_back( std::make_pair( &v_shut_down[ t + 1 - init_t ] ,
                                      v_ShutDownLimit[ t ] -
                                      get_operational_max_power( t ) ) );
      vars.push_back( std::make_pair( &v_start_up[ t - init_t ] ,
                                      std::max( 0.0 ,
                                                v_ShutDownLimit[ t ] -
                                                v_StartUpLimit[ t ] ) ) );
     }
    } else {
     if( ( t > 0 ) && ( t < f_time_horizon - 1 ) ) {
      vars.push_back( std::make_pair( &v_shut_down[ t + 1 - init_t ] ,
                                      v_ShutDownLimit[ t ] -
                                      get_operational_max_power( t ) ) );
      vars.push_back( std::make_pair( &v_start_up[ t - init_t ] ,
                                      v_StartUpLimit[ t ] -
                                      get_operational_max_power( t ) ) );
     }
    }
   }

   vars.push_back( std::make_pair( &v_commitment[ t ] ,
                                   get_operational_max_power( t ) ) );
   vars.push_back( std::make_pair( &v_active_power[ t ] , -1.0 ) );

   // if UCBlock has primary demand variables
   if( ( reserve_vars & 1u ) && ( ! v_PrimaryRho.empty() ) )
    vars.push_back( std::make_pair( &v_primary_spinning_reserve[ t ] ,
                                    -1.0 ) );

   // if UCBlock has secondary reserve variables
   if( ( reserve_vars & 2u ) && ( ! v_SecondaryRho.empty() ) )
    vars.push_back( std::make_pair( &v_secondary_spinning_reserve[ t ] ,
                                    -1.0 ) );

   MaxPower_Const[ cnstr_idx ].set_lhs( 0.0 );
   MaxPower_Const[ cnstr_idx ].set_rhs( Inf< double >() );
   MaxPower_Const[ cnstr_idx ].set_function(
    new LinearFunction( std::move( vars ) ) );

   if( t >= init_t ) {
    if( f_MinUpTime == 1 ) {
     if( ( t > 0 ) && ( t < f_time_horizon - 1 ) ) {

      vars.push_back( std::make_pair( &v_shut_down[ t + 1 - init_t ] ,
                                      std::max( 0.0 ,
                                                -v_ShutDownLimit[ t ] +
                                                v_StartUpLimit[ t ] ) ) );
      vars.push_back( std::make_pair( &v_start_up[ t - init_t ] ,
                                      v_StartUpLimit[ t ] -
                                      get_operational_max_power( t ) ) );

      vars.push_back( std::make_pair( &v_commitment[ t ] ,
                                      get_operational_max_power( t ) ) );
      vars.push_back( std::make_pair( &v_active_power[ t ] , -1.0 ) );

      // if UCBlock has primary demand variables
      if( ( reserve_vars & 1u ) && ( ! v_PrimaryRho.empty() ) )
       vars.push_back( std::make_pair( &v_primary_spinning_reserve[ t ] ,
                                       -1.0 ) );

      // if UCBlock has secondary reserve variables
      if( ( reserve_vars & 2u ) && ( ! v_SecondaryRho.empty() ) )
       vars.push_back( std::make_pair( &v_secondary_spinning_reserve[ t ] ,
                                       -1.0 ) );

      cnstr_idx++;

      MaxPower_Const[ cnstr_idx ].set_lhs( 0.0 );
      MaxPower_Const[ cnstr_idx ].set_rhs( Inf< double >() );
      MaxPower_Const[ cnstr_idx ].set_function(
       new LinearFunction( std::move( vars ) ) );
     }
    }
   }
  }

 } else if( ( AR & FormMsk ) == TForm ) {  // T formulation - - - - - - - - -

  Index max_power_cnstrs_size = 0;

  std::vector< int > v_T_RU;
  std::vector< int > v_T_RD;
  std::vector< int > v_K_SD;
  std::vector< int > v_K_SU;

  // TODO split ramp-up and ramp-down cnstrs
  if( ( ! v_DeltaRampUp.empty() ) && ( ! v_DeltaRampDown.empty() ) ) {

   v_T_RU.resize( f_time_horizon );
   v_T_RD.resize( f_time_horizon );
   v_K_SD.resize( f_time_horizon );
   v_K_SU.resize( f_time_horizon );

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {
    v_T_RU[ t ] = std::floor( ( get_operational_max_power( t ) -
                                v_ShutDownLimit[ t ] ) / v_DeltaRampUp[ t ] );
    v_T_RD[ t ] = std::floor( ( get_operational_max_power( t ) -
                                v_StartUpLimit[ t ] ) / v_DeltaRampDown[ t ] );
    v_K_SD[ t ] = std::min( f_InitUpDownTime - 1 , v_T_RD[ t ] );
    v_K_SD[ t ] = std::min( ( int ) ( f_time_horizon - t ) - 1 , v_K_SD[ t ] );
    v_K_SU[ t ] = std::min( f_InitUpDownTime - 2 -
                            std::max( 0 , v_K_SD[ t ] ) , v_T_RU[ t ] );
    v_K_SU[ t ] = std::min( ( int ) ( t ) - 1 , v_K_SU[ t ] );
   }

   max_power_cnstrs_size += f_time_horizon - init_t;
   // size bound constraints 4 + 5

   for( Index t = init_t ; t < f_time_horizon ; ++t )
    if( f_MinUpTime - 2 < v_T_RU[ t ] )
     max_power_cnstrs_size++;
   // size bound constraints 6

   for( Index t = init_t ; t < f_time_horizon ; ++t )
    if( v_K_SD[ t ] > 0 )
     max_power_cnstrs_size++;
  }

  // size bound constraints 0
  max_power_cnstrs_size += f_time_horizon;

  if( f_MinUpTime > 1 )
   max_power_cnstrs_size += f_time_horizon - init_t;
   // size bound constraints 1
  else
   max_power_cnstrs_size += 2 * f_time_horizon - 2 * init_t;
   // size bound constraints 2 + 3

  MaxPower_Const.resize( max_power_cnstrs_size );

  auto cnstr_idx = 0;

  // Bound constraints 0- - - - - - - - - - - - - - - - - - - - - - - - - - -

  for( Index t = 0 ; t < f_time_horizon ; ++t , ++cnstr_idx ) {

   vars.push_back( std::make_pair( &v_commitment[ t ] ,
                                   get_operational_max_power( t ) ) );
   vars.push_back( std::make_pair( &v_active_power[ t ] , -1.0 ) );

   // if UCBlock has primary demand variables
   if( ( reserve_vars & 1u ) && ( ! v_PrimaryRho.empty() ) )
    vars.push_back( std::make_pair( &v_primary_spinning_reserve[ t ] ,
                                    -1.0 ) );

   // if UCBlock has secondary reserve variables
   if( ( reserve_vars & 2u ) && ( ! v_SecondaryRho.empty() ) )
    vars.push_back( std::make_pair( &v_secondary_spinning_reserve[ t ] ,
                                    -1.0 ) );

   MaxPower_Const[ cnstr_idx ].set_lhs( 0.0 );
   MaxPower_Const[ cnstr_idx ].set_rhs( Inf< double >() );
   MaxPower_Const[ cnstr_idx ].set_function(
    new LinearFunction( std::move( vars ) ) );
  }

  // Bound constraints 1- - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( f_MinUpTime > 1 )
   for( Index t = init_t ; t < f_time_horizon ; ++t , ++cnstr_idx ) {

    vars.push_back( std::make_pair( &v_commitment[ t ] ,
                                    get_operational_max_power( t ) ) );
    vars.push_back( std::make_pair( &v_active_power[ t ] , -1.0 ) );

    // if UCBlock has primary demand variables
    if( ( reserve_vars & 1u ) && ( ! v_PrimaryRho.empty() ) )
     vars.push_back( std::make_pair( &v_primary_spinning_reserve[ t ] ,
                                     -1.0 ) );

    // if UCBlock has secondary reserve variables
    if( ( reserve_vars & 2u ) && ( ! v_SecondaryRho.empty() ) )
     vars.push_back( std::make_pair( &v_secondary_spinning_reserve[ t ] ,
                                     -1.0 ) );

    if( t >= init_t ) {
     vars.push_back( std::make_pair( &v_start_up[ t - init_t ] ,
                                     -( get_operational_max_power( t ) -
                                        v_StartUpLimit[ t ] ) ) );
     if( t < ( f_time_horizon - 1 ) )
      vars.push_back( std::make_pair( &v_shut_down[ t + 1 - init_t ] ,
                                      -( get_operational_max_power( t + 1 ) -
                                         v_ShutDownLimit[ t + 1 ] ) ) );
    }

    MaxPower_Const[ cnstr_idx ].set_lhs( 0.0 );
    MaxPower_Const[ cnstr_idx ].set_rhs( Inf< double >() );
    MaxPower_Const[ cnstr_idx ].set_function(
     new LinearFunction( std::move( vars ) ) );
   }

  if( f_MinUpTime == 1 ) {

   // Bound constraints 2 - - - - - - - - - - - - - - - - - - - - - - - - - -

   for( Index t = init_t ; t < f_time_horizon ; ++t , ++cnstr_idx ) {

    vars.push_back( std::make_pair( &v_commitment[ t ] ,
                                    get_operational_max_power( t ) ) );
    vars.push_back( std::make_pair( &v_active_power[ t ] , -1.0 ) );

    // if UCBlock has primary demand variables
    if( ( reserve_vars & 1u ) && ( ! v_PrimaryRho.empty() ) )
     vars.push_back( std::make_pair( &v_primary_spinning_reserve[ t ] ,
                                     -1.0 ) );

    // if UCBlock has secondary reserve variables
    if( ( reserve_vars & 2u ) && ( ! v_SecondaryRho.empty() ) )
     vars.push_back( std::make_pair( &v_secondary_spinning_reserve[ t ] ,
                                     -1.0 ) );

    if( t >= init_t ) {
     vars.push_back( std::make_pair( &v_start_up[ t - init_t ] ,
                                     -( get_operational_max_power( t ) -
                                        v_StartUpLimit[ t ] ) ) );
     if( v_ShutDownLimit[ t ] != v_StartUpLimit[ t ] )
      if( t < ( f_time_horizon - 1 ) )
       vars.push_back( std::make_pair(
        &v_shut_down[ t + 1 - init_t ] ,
        -( v_StartUpLimit[ t + 1 ] - v_ShutDownLimit[ t + 1 ] ) > 0 ?
        -( v_StartUpLimit[ t + 1 ] - v_ShutDownLimit[ t + 1 ] ) : 0.0 ) );
    }

    MaxPower_Const[ cnstr_idx ].set_lhs( 0.0 );
    MaxPower_Const[ cnstr_idx ].set_rhs( Inf< double >() );
    MaxPower_Const[ cnstr_idx ].set_function(
     new LinearFunction( std::move( vars ) ) );
   }

   // Bound constraints 3 - - - - - - - - - - - - - - - - - - - - - - - - - -

   for( Index t = init_t ; t < f_time_horizon ; ++t , ++cnstr_idx ) {

    vars.push_back( std::make_pair( &v_commitment[ t ] ,
                                    get_operational_max_power( t ) ) );
    vars.push_back( std::make_pair( &v_active_power[ t ] , -1.0 ) );

    // if UCBlock has primary demand variables
    if( ( reserve_vars & 1u ) && ( ! v_PrimaryRho.empty() ) )
     vars.push_back( std::make_pair( &v_primary_spinning_reserve[ t ] ,
                                     -1.0 ) );

    // if UCBlock has secondary reserve variables
    if( ( reserve_vars & 2u ) && ( ! v_SecondaryRho.empty() ) )
     vars.push_back( std::make_pair( &v_secondary_spinning_reserve[ t ] ,
                                     -1.0 ) );

    if( t >= init_t ) {
     if( t < ( f_time_horizon - 1 ) )
      vars.push_back( std::make_pair( &v_shut_down[ t + 1 - init_t ] ,
                                      -( get_operational_max_power( t + 1 ) -
                                         v_ShutDownLimit[ t + 1 ] ) ) );
     if( v_ShutDownLimit[ t ] != v_StartUpLimit[ t ] )
      vars.push_back( std::make_pair(
       &v_start_up[ t - init_t ] ,
       -( v_ShutDownLimit[ t ] - v_StartUpLimit[ t ] ) > 0 ?
       -( v_ShutDownLimit[ t ] - v_StartUpLimit[ t ] ) : 0.0 ) );
    }

    MaxPower_Const[ cnstr_idx ].set_lhs( 0.0 );
    MaxPower_Const[ cnstr_idx ].set_rhs( Inf< double >() );
    MaxPower_Const[ cnstr_idx ].set_function(
     new LinearFunction( std::move( vars ) ) );
   }
  }

  // TODO split ramp-up and ramp-down cnstrs
  if( ( ! v_DeltaRampUp.empty() ) && ( ! v_DeltaRampDown.empty() ) ) {

   // Bound constraints 4 - - - - - - - - - - - - - - - - - - - - - - - - - -

   for( Index t = init_t ; t < f_time_horizon ; ++t , ++cnstr_idx ) {

    vars.push_back( std::make_pair( &v_commitment[ t ] ,
                                    get_operational_max_power( t ) ) );
    vars.push_back( std::make_pair( &v_active_power[ t ] , -1.0 ) );

    // if UCBlock has primary demand variables
    if( ( reserve_vars & 1u ) && ( ! v_PrimaryRho.empty() ) )
     vars.push_back( std::make_pair( &v_primary_spinning_reserve[ t ] ,
                                     -1.0 ) );

    // if UCBlock has secondary reserve variables
    if( ( reserve_vars & 2u ) && ( ! v_SecondaryRho.empty() ) )
     vars.push_back( std::make_pair( &v_secondary_spinning_reserve[ t ] ,
                                     -1.0 ) );

    int min_RU = std::min( ( int ) ( f_MinUpTime ) - 2 , v_T_RU[ t ] );

    if( t >= init_t ) {
     if( t < ( f_time_horizon - 1 ) )
      vars.push_back( std::make_pair( &v_shut_down[ t + 1 - init_t ] ,
                                      -( get_operational_max_power( t + 1 ) -
                                         v_ShutDownLimit[ t + 1 ] ) ) );
     for( int s = 0 ; s < min_RU ; ++s )
      if( t - init_t >= s )
       vars.push_back( std::make_pair(
        &v_start_up[ t - s - init_t ] ,
        -( get_operational_max_power( t - s ) -
           v_StartUpLimit[ t - s ] -
           ( double ) ( s + 1 ) * v_DeltaRampUp[ t - s ] ) ) );
    }

    MaxPower_Const[ cnstr_idx ].set_lhs( 0.0 );
    MaxPower_Const[ cnstr_idx ].set_rhs( Inf< double >() );
    MaxPower_Const[ cnstr_idx ].set_function(
     new LinearFunction( std::move( vars ) ) );
   }

   // Bound constraints 5 - - - - - - - - - - - - - - - - - - - - - - - - - -

   for( Index t = init_t ; t < f_time_horizon ; ++t )
    if( ( f_MinUpTime - 2 ) < v_T_RU[ t ] ) {

     vars.push_back( std::make_pair( &v_commitment[ t ] ,
                                     get_operational_max_power( t ) ) );
     vars.push_back( std::make_pair( &v_active_power[ t ] , -1.0 ) );

     // if UCBlock has primary demand variables
     if( ( reserve_vars & 1u ) && ( ! v_PrimaryRho.empty() ) )
      vars.push_back( std::make_pair( &v_primary_spinning_reserve[ t ] ,
                                      -1.0 ) );

     // if UCBlock has secondary reserve variables
     if( ( reserve_vars & 2u ) && ( ! v_SecondaryRho.empty() ) )
      vars.push_back( std::make_pair( &v_secondary_spinning_reserve[ t ] ,
                                      -1.0 ) );

     int min_RU = std::min( ( int ) ( f_MinUpTime ) - 1 , v_T_RU[ t ] );

     if( t >= init_t ) {
      for( int s = 0 ; s < min_RU ; ++s )
       if( t - init_t >= s )
        vars.push_back( std::make_pair(
         &v_start_up[ t - s - init_t ] ,
         -( get_operational_max_power( t - s ) -
            v_StartUpLimit[ t - s ] -
            ( double ) ( s + 1 ) * v_DeltaRampUp[ t - s ] ) ) );
     }

     MaxPower_Const[ cnstr_idx ].set_lhs( 0.0 );
     MaxPower_Const[ cnstr_idx ].set_rhs( Inf< double >() );
     MaxPower_Const[ cnstr_idx ].set_function(
      new LinearFunction( std::move( vars ) ) );

     cnstr_idx++;
    }

   // Bound constraints 6 - - - - - - - - - - - - - - - - - - - - - - - - - -

   for( Index t = init_t ; t < f_time_horizon ; ++t )
    if( v_K_SD[ t ] > 0 ) {

     vars.push_back( std::make_pair( &v_commitment[ t ] ,
                                     get_operational_max_power( t ) ) );
     vars.push_back( std::make_pair( &v_active_power[ t ] , -1.0 ) );

     // if UCBlock has primary demand variables
     if( ( reserve_vars & 1u ) && ( ! v_PrimaryRho.empty() ) )
      vars.push_back( std::make_pair( &v_primary_spinning_reserve[ t ] ,
                                      -1.0 ) );

     // if UCBlock has secondary reserve variables
     if( ( reserve_vars & 2u ) && ( ! v_SecondaryRho.empty() ) )
      vars.push_back( std::make_pair( &v_secondary_spinning_reserve[ t ] ,
                                      -1.0 ) );

     if( t >= init_t ) {

      for( int s = 0 ; s < v_K_SD[ t ] ; ++s )
       if( ( t + 1 + s ) < f_time_horizon )
        vars.push_back( std::make_pair(
         &v_shut_down[ t + 1 + s - init_t ] ,
         -( get_operational_max_power( t + 1 + s ) -
            v_ShutDownLimit[ t + 1 + s ] -
            ( int ) ( s + 1 ) * v_DeltaRampDown[ t + 1 + s ] ) ) );

      for( int s = 0 ; s < v_K_SU[ t ] ; ++s )
       if( t - init_t >= s )
        vars.push_back( std::make_pair(
         &v_start_up[ t - s - init_t ] ,
         -( get_operational_max_power( t - s ) -
            v_StartUpLimit[ t - s ] -
            ( int ) ( s + 1 ) * v_DeltaRampUp[ t - s ] ) ) );
     }

     MaxPower_Const[ cnstr_idx ].set_lhs( 0.0 );
     MaxPower_Const[ cnstr_idx ].set_rhs( Inf< double >() );
     MaxPower_Const[ cnstr_idx ].set_function(
      new LinearFunction( std::move( vars ) ) );

     cnstr_idx++;
    }
  }

 } else if( ( AR & FormMsk ) == DPForm ) {  // DP formulation - - - - - - - -

  MaxPower_Const.resize( v_P_h_k.size() );

  for( Index j = 0 ; j < v_P_h_k.size() ; ++j ) {

   auto t = v_P_h_k[ j ].first;
   for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
    if( ( v_P_h_k[ j ].second.first == v_Y_plus[ i ].first ) &&
        ( v_P_h_k[ j ].second.second == v_Y_plus[ i ].second ) ) {
     if( v_Y_plus[ i ].first == t + 1 )
      vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                      v_StartUpLimit[ t ] ) );
     else if( v_Y_plus[ i ].second == t + 1 )
      vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                      v_ShutDownLimit[ t ] ) );
     else
      vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                      get_operational_max_power( t ) ) );
    }

   vars.push_back( std::make_pair( &v_active_power_h_k[ j ] , -1.0 ) );

   MaxPower_Const[ j ].set_lhs( 0.0 );
   MaxPower_Const[ j ].set_rhs( Inf< double >() );
   MaxPower_Const[ j ].set_function( new LinearFunction( std::move( vars ) ) );
  }

 } else {  // pt, SU or SD formulations - - - - - - - - - - - - - - - - - - -

  std::vector< double > v_psi( v_P_h_k.size() );

  for( Index j = 0 ; j < v_P_h_k.size() ; ++j ) {

   Index t = v_P_h_k[ j ].first;

   v_psi[ j ] = get_operational_max_power( t );

   if( v_P_h_k[ j ].second.second <= f_time_horizon )

    if( ! v_DeltaRampDown.empty() )
     if( v_P_h_k[ j ].second.first > 0 )
      v_psi[ j ] = std::min( v_psi[ j ] ,
                             v_ShutDownLimit[ t ] +
                             v_DeltaRampDown[ t ] *
                             ( v_P_h_k[ j ].second.second - ( t + 1 ) ) );

   if( ! v_DeltaRampUp.empty() ) {
    if( ( v_P_h_k[ j ].second.first == 0 ) && ( f_InitUpDownTime > 0 ) )
     v_psi[ j ] = std::min( v_psi[ j ] ,
                            f_InitialPower + v_DeltaRampUp[ t ] * ( t + 1 ) );

    if( v_P_h_k[ j ].second.first > 0 )
     v_psi[ j ] = std::min( v_psi[ j ] ,
                            v_StartUpLimit[ t ] +
                            v_DeltaRampUp[ t ] *
                            ( ( t + 1 ) - v_P_h_k[ j ].second.first ) );
   }

   v_psi[ j ] = std::max( v_psi[ j ] , get_operational_min_power( t ) );
  }

  if( ( AR & FormMsk ) == ptForm ) {  // pt formulation - - - - - - - - - - -

   MaxPower_Const.resize( f_time_horizon );

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {

    vars.push_back( std::make_pair( &v_active_power[ t ] , -1.0 ) );

    // if UCBlock has primary demand variables
    if( ( reserve_vars & 1u ) && ( ! v_PrimaryRho.empty() ) )
     vars.push_back( std::make_pair( &v_primary_spinning_reserve[ t ] ,
                                     -1.0 ) );

    // if UCBlock has secondary reserve variables
    if( ( reserve_vars & 2u ) && ( ! v_SecondaryRho.empty() ) )
     vars.push_back( std::make_pair( &v_secondary_spinning_reserve[ t ] ,
                                     -1.0 ) );

    for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
     for( Index j = 0 ; j < v_P_h_k.size() ; ++j )
      if( v_P_h_k[ j ].first == t )
       if( ( v_P_h_k[ j ].second.first == v_Y_plus[ i ].first ) &&
           ( v_P_h_k[ j ].second.second == v_Y_plus[ i ].second ) ) {

        if( ( v_Y_plus[ i ].first < t + 1 ) &&
            ( t + 1 < v_Y_plus[ i ].second ) )
         vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                         v_psi[ j ] ) );
        if( f_MinUpTime >= 2 ) {
         if( ( v_Y_plus[ i ].first == t + 1 ) &&
             ( t + 1 <= v_Y_plus[ i ].second ) )
          vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                          v_StartUpLimit[ t ] ) );
         if( ( v_Y_plus[ i ].first <= t + 1 ) &&
             ( t + 1 == v_Y_plus[ i ].second ) )
          vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                          v_ShutDownLimit[ t ] ) );
        }
        if( f_MinUpTime == 1 ) {
         if( ( v_Y_plus[ i ].first == t + 1 ) &&
             ( t + 1 < v_Y_plus[ i ].second ) )
          vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                          v_StartUpLimit[ t ] ) );
         if( ( v_Y_plus[ i ].first < t + 1 ) &&
             ( t + 1 == v_Y_plus[ i ].second ) )
          vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                          v_ShutDownLimit[ t ] ) );
         if( ( v_Y_plus[ i ].first == t + 1 ) &&
             ( t + 1 == v_Y_plus[ i ].second ) )
          vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                          std::min( v_StartUpLimit[ t ] ,
                                                    v_ShutDownLimit[ t ] ) ) );
        }
       }

    MaxPower_Const[ t ].set_lhs( 0.0 );
    MaxPower_Const[ t ].set_rhs( Inf< double >() );
    MaxPower_Const[ t ].set_function( new LinearFunction( std::move( vars ) ) );
   }

  } else if( ( AR & FormMsk ) == SUForm ) {  // SU formulation- - - - - - - -

   MaxPower_Const.resize( v_P_h.size() );

   for( Index j = 0 ; j < v_P_h.size() ; ++j ) {

    auto t = v_P_h[ j ].first;
    for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
     if( ( v_P_h[ j ].second == v_Y_plus[ i ].first ) &&
         ( t + 1 <= v_Y_plus[ i ].second ) ) {

      if( t + 1 == v_Y_plus[ i ].first ) {

       if( v_Y_plus[ i ].first < v_Y_plus[ i ].second )
        vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                        v_StartUpLimit[ t ] ) );

       if( v_Y_plus[ i ].first == v_Y_plus[ i ].second )
        vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                        std::min( v_StartUpLimit[ t ] ,
                                                  v_ShutDownLimit[ t ] ) ) );

      } else {

       if( t + 1 == v_Y_plus[ i ].second )
        vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                        v_ShutDownLimit[ t ] ) );

       if( t + 1 < v_Y_plus[ i ].second )
        for( Index s = 0 ; s < v_P_h_k.size() ; ++s )
         if( ( t == v_P_h_k[ s ].first ) &&
             ( v_P_h_k[ s ].second.first == v_Y_plus[ i ].first ) &&
             ( v_P_h_k[ s ].second.second == v_Y_plus[ i ].second ) )
          vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                          v_psi[ s ] ) );
      }
     }

    vars.push_back( std::make_pair( &v_active_power_h[ j ] , -1.0 ) );

    MaxPower_Const[ j ].set_lhs( 0.0 );
    MaxPower_Const[ j ].set_rhs( Inf< double >() );
    MaxPower_Const[ j ].set_function( new LinearFunction( std::move( vars ) ) );
   }

  } else if( ( AR & FormMsk ) == SDForm ) {  // SD formulation- - - - - - - -

   MaxPower_Const.resize( v_P_k.size() );

   for( Index j = 0 ; j < v_P_k.size() ; ++j ) {

    auto t = v_P_k[ j ].first;
    for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
     if( ( v_P_k[ j ].second == v_Y_plus[ i ].second ) &&
         ( v_Y_plus[ i ].first <= t + 1 ) ) {

      if( t + 1 == v_Y_plus[ i ].second ) {

       if( v_Y_plus[ i ].first < v_Y_plus[ i ].second )
        vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                        v_ShutDownLimit[ t ] ) );

       if( v_Y_plus[ i ].first == v_Y_plus[ i ].second )
        vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                        std::min( v_StartUpLimit[ t ] ,
                                                  v_ShutDownLimit[ t ] ) ) );

      } else {

       if( t + 1 == v_Y_plus[ i ].first )
        vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                        v_StartUpLimit[ t ] ) );

       if( v_Y_plus[ i ].first < t + 1 )
        for( Index s = 0 ; s < v_P_h_k.size() ; ++s )
         if( ( t == v_P_h_k[ s ].first ) &&
             ( v_P_h_k[ s ].second.first == v_Y_plus[ i ].first ) &&
             ( v_P_h_k[ s ].second.second == v_Y_plus[ i ].second ) )
          vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                          v_psi[ s ] ) );
      }
     }

    vars.push_back( std::make_pair( &v_active_power_k[ j ] , -1.0 ) );

    MaxPower_Const[ j ].set_lhs( 0.0 );
    MaxPower_Const[ j ].set_rhs( Inf< double >() );
    MaxPower_Const[ j ].set_function( new LinearFunction( std::move( vars ) ) );
   }
  }
 }

 add_static_constraint( MaxPower_Const , "MaxPower_Const_Thermal" );

 if( reserve_vars & 1u )  // if UCBlock has primary demand variables
  if( ! v_PrimaryRho.empty() ) {  // if unit produces any primary reserve

   // Initializing primary rho fraction constraints - - - - - - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   PrimaryRho_Const.resize( f_time_horizon );

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {

    if( ! v_PrimaryRho.empty() )
     vars.push_back( std::make_pair( &v_active_power[ t ] ,
                                     v_PrimaryRho[ t ] ) );
    else
     vars.push_back( std::make_pair( &v_active_power[ t ] , 0.0 ) );

    vars.push_back( std::make_pair( &v_primary_spinning_reserve[ t ] , -1.0 ) );

    PrimaryRho_Const[ t ].set_lhs( 0.0 );
    PrimaryRho_Const[ t ].set_rhs( Inf< double >() );
    PrimaryRho_Const[ t ].set_function(
     new LinearFunction( std::move( vars ) ) );
   }

   add_static_constraint( PrimaryRho_Const , "PrimaryRho_Const_Thermal" );
  }

 if( reserve_vars & 2u )  // if UCBlock has secondary demand variables
  if( ! v_SecondaryRho.empty() ) {  // if unit produces any secondary reserve

   // Initializing secondary rho fraction constraints - - - - - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   SecondaryRho_Const.resize( f_time_horizon );

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {

    if( ! v_SecondaryRho.empty() )
     vars.push_back( std::make_pair( &v_active_power[ t ] ,
                                     v_SecondaryRho[ t ] ) );
    else
     vars.push_back( std::make_pair( &v_active_power[ t ] , 0.0 ) );

    vars.push_back( std::make_pair( &v_secondary_spinning_reserve[ t ] ,
                                    -1.0 ) );

    SecondaryRho_Const[ t ].set_lhs( 0.0 );
    SecondaryRho_Const[ t ].set_rhs( Inf< double >() );
    SecondaryRho_Const[ t ].set_function(
     new LinearFunction( std::move( vars ) ) );
   }

   add_static_constraint( SecondaryRho_Const , "SecondaryRho_Const_Thermal" );
  }

 // ZOConstraints - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( generate_ZOConstraints ) {

  // the commitment bound constraints
  Commitment_bound_Const.resize( f_time_horizon );

  for( Index t = 0 ; t < f_time_horizon ; ++t )
   Commitment_bound_Const[ t ].set_variable( &v_commitment[ t ] );

  add_static_constraint( Commitment_bound_Const , "Commitment_bound_Thermal" );

  auto startup_shutdown_size = f_time_horizon - init_t;

  // the startup binary bound constraints
  StartUp_Binary_bound_Const.resize( startup_shutdown_size );

  for( Index t = 0 ; t < startup_shutdown_size ; ++t )
   StartUp_Binary_bound_Const[ t ].set_variable( &v_start_up[ t ] );

  add_static_constraint( StartUp_Binary_bound_Const ,
                         "StartUp_binary_bound_Thermal" );

  // the shut-down binary bound constraints
  ShutDown_Binary_bound_Const.resize( startup_shutdown_size );

  for( Index t = 0 ; t < startup_shutdown_size ; ++t )
   ShutDown_Binary_bound_Const[ t ].set_variable( &v_shut_down[ t ] );

  add_static_constraint( ShutDown_Binary_bound_Const ,
                         "ShoutDown_binary_bound_Thermal" );
 }

 // BoxConstraint - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( init_t > 0 ) && ( f_InitUpDownTime > 0 ) &&
     ( f_InitUpDownTime < f_MinUpTime ) ) {

  // the commitment fixed to one BoxConstraints
  Commitment_fixed_to_One_Const.resize( f_time_horizon );

  for( Index t = 0 ; t < init_t ; ++t ) {
   Commitment_fixed_to_One_Const[ t ].set_both( 1 );
   Commitment_fixed_to_One_Const[ t ].set_variable( &v_commitment[ t ] );
  }

  add_static_constraint( Commitment_fixed_to_One_Const ,
                         "Commitment_fixed_to_one_Thermal" );
 }

 if( AR & PCuts ) {

  // Initial perspective cuts constraints - - - - - - - - - - - - - - - - - -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  auto cnstr_idx = 0;

  if( ( ( AR & FormMsk ) == tbinForm ) ||  // 3bin formulation- - - - - - - -
      ( ( AR & FormMsk ) == TForm ) ) {  // T formulation - - - - - - - - - -

   Init_PC_Const.resize( 2 * f_time_horizon );

   for( Index t = 0 ; t < f_time_horizon ; ++t )
    for( Index k = 0 ; k <= 1 ; ++k ) {

     auto value = ( k == 0 ? get_operational_min_power( t )
                           : get_operational_max_power( t ) );

     vars.push_back( std::make_pair( &v_active_power[ t ] , 2 * value ) );
     vars.push_back( std::make_pair( &v_cut[ t ] , -1.0 ) );
     vars.push_back( std::make_pair( &v_commitment[ t ] ,
                                     -std::pow( value , 2 ) ) );

     Init_PC_Const[ cnstr_idx ].set_lhs( -Inf< double >() );
     Init_PC_Const[ cnstr_idx ].set_rhs( 0.0 );
     Init_PC_Const[ cnstr_idx ].set_function(
      new LinearFunction( std::move( vars ) ) );

     cnstr_idx++;
    }

  } else if( ( AR & FormMsk ) == ptForm ) {  // pt formulation- - - - - - - -

   Init_PC_Const.resize( 2 * f_time_horizon );

   for( Index t = 0 ; t < f_time_horizon ; ++t )
    for( Index k = 0 ; k <= 1 ; ++k ) {

     auto value = ( k == 0 ? get_operational_min_power( t )
                           : get_operational_max_power( t ) );

     vars.push_back( std::make_pair( &v_active_power[ t ] , 2 * value ) );
     vars.push_back( std::make_pair( &v_cut[ t ] , -1.0 ) );

     for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
      if( ( v_Y_plus[ i ].first <= t + 1 ) &&
          ( t + 1 <= v_Y_plus[ i ].second ) )
       vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                       -std::pow( value , 2 ) ) );

     Init_PC_Const[ cnstr_idx ].set_lhs( -Inf< double >() );
     Init_PC_Const[ cnstr_idx ].set_rhs( 0.0 );
     Init_PC_Const[ cnstr_idx ].set_function(
      new LinearFunction( std::move( vars ) ) );

     cnstr_idx++;
    }

  } else if( ( AR & FormMsk ) == DPForm ) {  // DP formulation- - - - - - - -

   Init_PC_Const.resize( 2 * v_P_h_k.size() );

   for( Index j = 0 ; j < v_P_h_k.size() ; ++j )
    for( Index k = 0 ; k <= 1 ; ++k ) {

     auto t = v_P_h_k[ j ].first;
     auto value = ( k == 0 ? get_operational_min_power( t )
                           : get_operational_max_power( t ) );

     vars.push_back( std::make_pair( &v_active_power_h_k[ j ] , 2 * value ) );

     for( Index s = 0 ; s < v_P_h_k.size() ; ++s )
      if( ( v_P_h_k[ j ].second.first == v_Z_h_k[ s ].second.first ) &&
          ( v_P_h_k[ j ].second.second == v_Z_h_k[ s ].second.second ) )
       if( v_Z_h_k[ s ].first == t )
        vars.push_back( std::make_pair( &v_cut_h_k[ s ] , -1.0 ) );

     for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
      if( ( v_Y_plus[ i ].first == v_P_h_k[ j ].second.first ) &&
          ( v_Y_plus[ i ].second == v_P_h_k[ j ].second.second ) )
       vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                       -std::pow( value , 2 ) ) );

     Init_PC_Const[ cnstr_idx ].set_lhs( -Inf< double >() );
     Init_PC_Const[ cnstr_idx ].set_rhs( 0.0 );
     Init_PC_Const[ cnstr_idx ].set_function(
      new LinearFunction( std::move( vars ) ) );

     cnstr_idx++;
    }

  } else if( ( AR & FormMsk ) == SUForm ) {  // SU formulation- - - - - - - -

   Init_PC_Const.resize( 2 * v_P_h.size() );

   for( Index j = 0 ; j < v_P_h.size() ; ++j )
    for( Index k = 0 ; k <= 1 ; ++k ) {

     auto t = v_P_h[ j ].first;
     auto value = ( k == 0 ? get_operational_min_power( t )
                           : get_operational_max_power( t ) );

     vars.push_back( std::make_pair( &v_active_power_h[ j ] , 2 * value ) );

     for( Index s = 0 ; s < v_P_h.size() ; ++s )
      if( ( v_P_h[ j ].second == v_Z_h[ s ].second ) &&
          ( v_P_h[ j ].first == v_Z_h[ s ].first ) )
       vars.push_back( std::make_pair( &v_cut_h[ s ] , -1.0 ) );

     for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
      if( v_P_h[ j ].second == v_Y_plus[ i ].first )
       if( t + 1 <= v_Y_plus[ i ].second )
        vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                        -std::pow( value , 2 ) ) );

     Init_PC_Const[ cnstr_idx ].set_lhs( -Inf< double >() );
     Init_PC_Const[ cnstr_idx ].set_rhs( 0.0 );
     Init_PC_Const[ cnstr_idx ].set_function(
      new LinearFunction( std::move( vars ) ) );

     cnstr_idx++;
    }

  } else if( ( AR & FormMsk ) == SDForm ) {  // SD formulation- - - - - - - -

   Init_PC_Const.resize( 2 * v_P_k.size() );

   for( Index j = 0 ; j < v_P_k.size() ; ++j )
    for( Index k = 0 ; k <= 1 ; ++k ) {

     auto t = v_P_k[ j ].first;
     auto value = ( k == 0 ? get_operational_min_power( t )
                           : get_operational_max_power( t ) );

     vars.push_back( std::make_pair( &v_active_power_k[ j ] , 2 * value ) );

     for( Index s = 0 ; s < v_P_k.size() ; ++s )
      if( ( v_P_k[ j ].second == v_Z_k[ s ].second ) &&
          ( v_P_k[ j ].first == v_Z_k[ s ].first ) )
       vars.push_back( std::make_pair( &v_cut_k[ s ] , -1.0 ) );

     for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
      if( v_P_k[ j ].second == v_Y_plus[ i ].second )
       if( v_Y_plus[ i ].first <= t + 1 )
        vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                        -std::pow( value , 2 ) ) );

     Init_PC_Const[ cnstr_idx ].set_lhs( -Inf< double >() );
     Init_PC_Const[ cnstr_idx ].set_rhs( 0.0 );
     Init_PC_Const[ cnstr_idx ].set_function(
      new LinearFunction( std::move( vars ) ) );

     cnstr_idx++;
    }
  }

  add_static_constraint( Init_PC_Const , "Init_PC_Const_Thermal" );

  // Constraints connecting perspective cuts variables of 3bin with those of
  // DP, SU and SD formulations- - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( ( AR & FormMsk ) == tbinForm ) ||  // 3bin formulation- - - - - - - -
      ( ( AR & FormMsk ) == TForm ) ||  // T formulation- - - - - - - - - - -
      ( ( AR & FormMsk ) == ptForm ) ) {  // pt formulation - - - - - - - - -
   ;  // does nothing
  } else {  // DP, SU or SD formulations- - - - - - - - - - - - - - - - - - -

   Eq_PC_Const.resize( f_time_horizon );

   if( ( AR & FormMsk ) == DPForm ) {  // DP formulation- - - - - - - - - - -

    for( Index t = 0 ; t < f_time_horizon ; ++t ) {

     vars.push_back( std::make_pair( &v_cut[ t ] , 1.0 ) );

     for( Index j = 0 ; j < v_Z_h_k.size() ; ++j )
      if( v_Z_h_k[ j ].first == t )
       vars.push_back( std::make_pair( &v_cut_h_k[ j ] , -1.0 ) );

     Eq_PC_Const[ t ].set_both( 0.0 );
     Eq_PC_Const[ t ].set_function( new LinearFunction( std::move( vars ) ) );
    }

   } else if( ( AR & FormMsk ) == SUForm ) {  // SU formulation - - - - - - -

    for( Index t = 0 ; t < f_time_horizon ; ++t ) {

     vars.push_back( std::make_pair( &v_cut[ t ] , 1.0 ) );

     for( Index j = 0 ; j < v_Z_h.size() ; ++j )
      if( v_Z_h[ j ].first == t )
       vars.push_back( std::make_pair( &v_cut_h[ j ] , -1.0 ) );

     Eq_PC_Const[ t ].set_both( 0.0 );
     Eq_PC_Const[ t ].set_function( new LinearFunction( std::move( vars ) ) );
    }

   } else if( ( AR & FormMsk ) == SDForm ) {  // SD formulation - - - - - - -

    for( Index t = 0 ; t < f_time_horizon ; ++t ) {

     vars.push_back( std::make_pair( &v_cut[ t ] , 1.0 ) );

     for( Index j = 0 ; j < v_Z_k.size() ; ++j )
      if( v_Z_k[ j ].first == t )
       vars.push_back( std::make_pair( &v_cut_k[ j ] , -1.0 ) );

     Eq_PC_Const[ t ].set_both( 0.0 );
     Eq_PC_Const[ t ].set_function( new LinearFunction( std::move( vars ) ) );
    }
   }
  }

  add_static_constraint( Eq_PC_Const , "Eq_PC_Const_Thermal" );
 }

 set_constraints_generated();

}  // end( ThermalUnitBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::generate_dynamic_constraints( Configuration * dycc )
{
 if( AR & PCuts ) {
  double tol = 1e-3;  // threshold parameter for P/C separation
  double eps = 1e-4;  // tolerance value to consider a binary variable

  auto extract_parameters = [ & tol , & eps ]( Configuration * c )
   -> bool {
   if( auto tc = dynamic_cast< SimpleConfiguration< double > * >( c ) ) {
    tol = tc->f_value;
    return( true );
   }
   if( auto tc = dynamic_cast<
    SimpleConfiguration< std::pair< double , double > > * >( c ) ) {
    tol = tc->f_value.first;
    eps = tc->f_value.second;
    return( true );
   }
   return( false );
  };

  if( ( ! extract_parameters( dycc ) ) && f_BlockConfig )
   // if the given Configuration is not valid, try the one from the BlockConfig
   extract_parameters( f_BlockConfig->f_dynamic_constraints_Configuration );

  LinearFunction::v_coeff_pair vars;

  for( Index t = 0 ; t < f_time_horizon ; ++t )

   if( v_commitment[ t ].get_value() > eps ) {

    if( v_cut[ t ].get_value() <
        ( std::pow( v_active_power[ t ].get_value() , 2 ) /
          v_commitment[ t ].get_value() ) - tol ) {

     std::list< FRowConstraint > cut( 1 );

     vars.push_back( std::make_pair( &v_active_power[ t ] ,
                                     2 * ( v_active_power[ t ].get_value() /
                                           v_commitment[ t ].get_value() ) ) );
     vars.push_back( std::make_pair( &v_cut[ t ] , -1.0 ) );
     vars.push_back(
      std::make_pair( &v_commitment[ t ] ,
                      -( std::pow( v_active_power[ t ].get_value() , 2 ) /
                         std::pow( v_commitment[ t ].get_value() , 2 ) ) ) );

     cut.front().set_lhs( -Inf< double >() );
     cut.front().set_rhs( 0.0 );
     cut.front().set_function(
      new LinearFunction( std::move( vars ) , eNoMod ) );

     add_dynamic_constraints( PC_cuts , cut , eNoBlck );
    }
   }

  add_dynamic_constraint( PC_cuts , "PC_cuts_Thermal" );
 }
}  // end( ThermalUnitBlock::generate_dynamic_constraints )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::generate_objective( Configuration * objc )
{
 if( objective_generated() )  // Objective has already been generated
  return;                     // nothing to do

 // initialize Objective
 //
 // the order of the variables in the Objective Function is:
 //
 // - first f_time_horizon - init_t start-up variables
 //
 // - then f_time_horizon active power variables (which may have the
 //   nonzero quadratic cost coefficient, while the others do not)
 //
 // - then f_time_horizon commitment variables
 //
 // - then possibly f_time_horizon primary reserve variables
 //
 // - then possibly f_time_horizon secondary reserve variables
 //
 // this arrangement is exploited in add_Modification to easily map
 // indices in the coefficients of the Objective Function back into
 // indices of the original variables (and figure out the kind of variable)

 if( v_commitment.size() != f_time_horizon )
  throw( std::logic_error(
   "ThermalUnitBlock::generate_objective: v_commitment must have "
   "size equal to the time horizon." ) );

 if( v_active_power.size() != f_time_horizon )
  throw( std::logic_error(
   "ThermalUnitBlock::generate_objective: v_active_power must have "
   "size equal to the time horizon." ) );

 if( v_start_up.size() != f_time_horizon - init_t )
  throw( std::logic_error(
   "ThermalUnitBlock::generate_objective: v_start_up must have "
   "size equal to the time horizon - init_t." ) );

 DQuadFunction::v_coeff_triple vars;

 if( f_InvestmentCost != 0 )
  vars.push_back( std::make_tuple( &design , f_InvestmentCost , 0.0 ) );

 // add the start-up variables- - - - - - - - - - - - - - - - - - - - - - - -
 for( Index t = init_t ; t < f_time_horizon ; ++t )
  vars.push_back( std::make_tuple( &v_start_up[ t - init_t ] ,
                                   f_scale * v_StartUpCost[ t ] , 0.0 ) );

 // add the active power variables- - - - - - - - - - - - - - - - - - - - - -
 for( Index t = 0 ; t < f_time_horizon ; ++t )
  vars.push_back( std::make_tuple( &v_active_power[ t ] ,
                                   f_scale * v_LinearTerm[ t ] ,
                                   AR & PCuts ? 0.0 : f_scale * v_QuadTerm[ t ] ) );

 // add the commitment variables- - - - - - - - - - - - - - - - - - - - - - -
 for( Index t = 0 ; t < f_time_horizon ; ++t )
  vars.push_back( std::make_tuple( &v_commitment[ t ] ,
                                   f_scale * v_ConstTerm[ t ] , 0.0 ) );

 if( ( reserve_vars & 1u ) && ( ! v_primary_spinning_reserve.empty() ) ) {
  // add the primary spinning reserve variables - - - - - - - - - - - - - - -
  if( v_primary_spinning_reserve.size() != f_time_horizon )
   throw( std::logic_error( "ThermalUnitBlock::generate_objective: v_primary_"
                            "spinning_reserve must have size equal to the "
                            "time horizon." ) );

  if( v_PrimarySpinningReserveCost.empty() )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    vars.push_back( std::make_tuple( &v_primary_spinning_reserve[ t ] ,
                                     0.0 , 0.0 ) );
  else
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    vars.push_back( std::make_tuple( &v_primary_spinning_reserve[ t ] ,
                                     f_scale *
                                     v_PrimarySpinningReserveCost[ t ] ,
                                     0.0 ) );
 }

 if( ( reserve_vars & 2u ) && ( ! v_secondary_spinning_reserve.empty() ) ) {
  // add the secondary spinning reserve variables - - - - - - - - - - - - - -
  if( v_secondary_spinning_reserve.size() != f_time_horizon )
   throw( std::logic_error( "ThermalUnitBlock::generate_objective: v_secondary"
                            "_spinning_reserve must have size equal to the "
                            "time horizon." ) );

  if( v_SecondarySpinningReserveCost.empty() )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    vars.push_back( std::make_tuple( &v_secondary_spinning_reserve[ t ] ,
                                     0.0 , 0.0 ) );
  else
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    vars.push_back( std::make_tuple( &v_secondary_spinning_reserve[ t ] ,
                                     f_scale *
                                     v_SecondarySpinningReserveCost[ t ] ,
                                     0.0 ) );
 }

 if( AR & PCuts )
  // add the perspective cuts variables - - - - - - - - - - - - - - - - - - -
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   vars.push_back( std::make_tuple( &v_cut[ t ] ,
                                    f_scale * v_QuadTerm[ t ] , 0.0 ) );

 objective.set_function( new DQuadFunction( std::move( vars ) ) );
 objective.set_sense( Objective::eMin );

 // set Block objective
 set_objective( &objective );

 set_objective_generated();

}  // end( ThermalUnitBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*---------------- METHODS FOR CHECKING THE ThermalUnitBlock ---------------*/
/*--------------------------------------------------------------------------*/

bool ThermalUnitBlock::is_feasible( bool useabstract , Configuration * fsbc )
{
 // Retrieve the tolerance and the type of violation.
 double tol = 0;
 bool rel_viol = true;

 // Try to extract, from "c", the parameters that determine feasibility.
 // If it succeeds, it sets the values of the parameters and returns
 // true. Otherwise, it returns false.
 auto extract_parameters = [ & tol , & rel_viol ]( Configuration * c )
  -> bool {
  if( auto tc = dynamic_cast< SimpleConfiguration< double > * >( c ) ) {
   tol = tc->f_value;
   return( true );
  }
  if( auto tc = dynamic_cast< SimpleConfiguration< std::pair< double , int > > * >( c ) ) {
   tol = tc->f_value.first;
   rel_viol = tc->f_value.second;
   return( true );
  }
  return( false );
 };

 if( ( ! extract_parameters( fsbc ) ) && f_BlockConfig )
  // if the given Configuration is not valid, try the one from the BlockConfig
  extract_parameters( f_BlockConfig->f_is_feasible_Configuration );

 return(
  UnitBlock::is_feasible( useabstract )
  // Variables
  && ColVariable::is_feasible( v_start_up , tol )
  && ColVariable::is_feasible( v_shut_down , tol )
  && ColVariable::is_feasible( v_primary_spinning_reserve , tol )
  && ColVariable::is_feasible( v_secondary_spinning_reserve , tol )
  && ColVariable::is_feasible( v_commitment , tol )
  && ColVariable::is_feasible( v_commitment_plus , tol )
  && ColVariable::is_feasible( v_commitment_minus , tol )
  && ColVariable::is_feasible( v_active_power , tol )
  && ColVariable::is_feasible( v_active_power_h_k , tol )
  && ColVariable::is_feasible( v_active_power_h , tol )
  && ColVariable::is_feasible( v_active_power_k , tol )
  && ColVariable::is_feasible( v_cut , tol )
  && ColVariable::is_feasible( v_cut_h_k , tol )
  && ColVariable::is_feasible( v_cut_h , tol )
  && ColVariable::is_feasible( v_cut_k , tol )
  // Constraints: notice that the ZOConstraints are not checked, since the
  // corresponding check is made on the ColVariable
  && RowConstraint::is_feasible( CommitmentDesign_Const , tol , rel_viol )
  && RowConstraint::is_feasible( StartUp_ShutDown_Variables_Const , tol , rel_viol )
  && RowConstraint::is_feasible( StartUp_Const , tol , rel_viol )
  && RowConstraint::is_feasible( ShutDown_Const , tol , rel_viol )
  && RowConstraint::is_feasible( RampUp_Const , tol , rel_viol )
  && RowConstraint::is_feasible( RampDown_Const , tol , rel_viol )
  && RowConstraint::is_feasible( MinPower_Const , tol , rel_viol )
  && RowConstraint::is_feasible( MaxPower_Const , tol , rel_viol )
  && RowConstraint::is_feasible( PrimaryRho_Const , tol , rel_viol )
  && RowConstraint::is_feasible( SecondaryRho_Const , tol , rel_viol )
  && RowConstraint::is_feasible( Eq_ActivePower_Const , tol , rel_viol )
  && RowConstraint::is_feasible( Eq_Commitment_Const , tol , rel_viol )
  && RowConstraint::is_feasible( Eq_StartUp_Const , tol , rel_viol )
  && RowConstraint::is_feasible( Eq_ShutDown_Const , tol , rel_viol )
  && RowConstraint::is_feasible( Network_Const , tol , rel_viol )
  && RowConstraint::is_feasible( Init_PC_Const , tol , rel_viol )
  && RowConstraint::is_feasible( Eq_PC_Const , tol , rel_viol )
  && RowConstraint::is_feasible( PC_cuts , tol , rel_viol )
  && RowConstraint::is_feasible( Commitment_fixed_to_One_Const , tol , rel_viol ) );

}  // end( ThermalUnitBlock::is_feasible )

/*--------------------------------------------------------------------------*/
/*-------- METHODS FOR LOADING, PRINTING & SAVING THE ThermalUnitBlock -----*/
/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::serialize( netCDF::NcGroup & group ) const
{
 UnitBlock::serialize( group );

 // Serialize scalar variables

 if( f_InvestmentCost != 0 )
  ::serialize( group , "InvestmentCost" , netCDF::NcDouble() ,
               f_InvestmentCost );

 if( f_Capacity != 0 )
  ::serialize( group , "Capacity" , netCDF::NcDouble() , f_Capacity );

 ::serialize( group , "InitialPower" , netCDF::NcDouble() , f_InitialPower );
 ::serialize( group , "MinUpTime" , netCDF::NcUint() , f_MinUpTime );
 ::serialize( group , "MinDownTime" , netCDF::NcUint() , f_MinDownTime );
 ::serialize( group , "InitUpDownTime" , netCDF::NcInt() , f_InitUpDownTime );

 // Serialize one-dimensional variables

 auto TimeHorizon = group.getDim( "TimeHorizon" );
 auto NumberIntervals = group.getDim( "NumberIntervals" );

 /* This lambda identifies the appropriate dimension for the given variable
  * (whose name is "var_name") and serializes the variable. The variable may
  * have any of the following dimensions: TimeHorizon, NumberIntervals,
  * 1. "allow_scalar_var" indicates whether the variable can be serialized as
  * a scalar variable (in which case the variable must have dimension 1). */
 auto serialize = [ &group , &TimeHorizon , &NumberIntervals ]
  ( const std::string & var_name , const std::vector< double > & data ,
    const netCDF::NcType & ncType = netCDF::NcDouble() ,
    bool allow_scalar_var = true ) {
  if( data.empty() )
   return;
  netCDF::NcDim dimension;
  if( data.size() == TimeHorizon.getSize() )
   dimension = TimeHorizon;
  else if( data.size() == NumberIntervals.getSize() )
   dimension = NumberIntervals;
  else if( data.size() != 1 )
   throw( std::logic_error( "ThermalUnitBlock::serialize: invalid dimension "
                            "for variable " + var_name + ": " +
                            std::to_string( data.size() ) +
                            ". Its dimension must be one of the following: "
                            "TimeHorizon, NumberIntervals, 1." ) );

  ::serialize( group , var_name , ncType , dimension , data ,
               allow_scalar_var );
 };

 serialize( "MinPower" , v_MinPower );
 serialize( "MaxPower" , v_MaxPower );
 serialize( "Availability" , v_Availability );
 serialize( "DeltaRampUp" , v_DeltaRampUp );
 serialize( "DeltaRampDown" , v_DeltaRampDown );
 serialize( "PrimaryRho" , v_PrimaryRho );
 serialize( "SecondaryRho" , v_SecondaryRho );
 serialize( "QuadTerm" , v_QuadTerm );
 serialize( "LinearTerm" , v_LinearTerm );
 serialize( "ConstTerm" , v_ConstTerm );
 serialize( "StartUpCost" , v_StartUpCost );
 serialize( "FixedConsumption" , v_FixedConsumption );
 serialize( "InertiaCommitment" , v_InertiaCommitment );
 serialize( "StartUpLimit" , v_StartUpLimit );
 serialize( "ShutDownLimit" , v_ShutDownLimit );

}  // end( ThermalUnitBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::add_Modification( sp_Mod mod , ChnlName chnl )
{
 if( mod->concerns_Block() ) {
  mod->concerns_Block( false );
  guts_of_add_Modification( mod.get() , chnl );
 }

 Block::add_Modification( mod , chnl );

}  // end( ThermalUnitBlock::add_Modification )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::update_availability_dependents( Index t ,
                                                       ModParam issueAMod )
{
 if( ! constraints_generated() )
  return;

 not_ModBlock( issueAMod );

 // MaxPower_Const: the commitment variable is in position 0
 static_cast< LinearFunction * >( MaxPower_Const[ t ].get_function()
 )->modify_coefficient( 0 , get_operational_max_power( t ) , issueAMod );

 // MinPower_Const: the commitment variable is in position 0
 static_cast< LinearFunction * >( MinPower_Const[ t ].get_function()
 )->modify_coefficient( 0 , -get_operational_min_power( t ) , issueAMod );

 // RampUp_Const
 if( init_t == 0 ) {

  double coefficient = get_operational_min_power( t );
  if( t == 0 )
   coefficient *= -1.0;

  auto f = static_cast< LinearFunction * >( RampUp_Const[ t ].get_function() );
  auto var_index = f->is_active( &v_start_up[ t ] );
  assert( var_index < f->get_num_active_var() );
  f->modify_coefficient( var_index , coefficient , issueAMod );

 } else if( init_t > 0 ) {

  auto depends_on_min_power = ( t > init_t );
  depends_on_min_power |= ( t == init_t ) &&
                          ( ( ( f_InitUpDownTime < 0 ) &&
                              ( -f_InitUpDownTime < f_MinDownTime ) ) ||
                            ( ( f_InitUpDownTime > 0 ) &&
                              ( f_InitUpDownTime < f_MinUpTime ) ) );

  if( depends_on_min_power ) {

   auto coefficient = get_operational_min_power( t );
   if( t == init_t )
    coefficient *= -1.0;

   auto f = static_cast< LinearFunction * >( RampUp_Const[ t ].get_function() );
   auto var_index = f->is_active( &v_start_up[ t - init_t ] );
   assert( var_index < f->get_num_active_var() );
   f->modify_coefficient( var_index , coefficient , issueAMod );
  }
 }

 // RampDown_Const
 if( ( ( init_t == 0 ) && ( t == 0 ) ) ||
     ( ( init_t > 0 ) && ( t >= init_t ) ) ) {

  auto f = static_cast< LinearFunction * >( RampDown_Const[ t ].get_function() );
  auto var_index = f->is_active( &v_shut_down[ t - init_t ] );
  assert( var_index < f->get_num_active_var() );
  auto coefficient = get_operational_min_power( t );
  f->modify_coefficient( var_index , coefficient , issueAMod );
 }
}  // end( ThermalUnitBlock::update_availability_dependents )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_availability( MF_dbl_it values ,
                                         Subset && subset ,
                                         const bool ordered ,
                                         ModParam issuePMod ,
                                         ModParam issueAMod )
{
 if( subset.empty() )
  return;

 if( v_Availability.empty() ) {
  if( std::all_of( values ,
                   values + subset.size() ,
                   []( double cst ) { return( cst == 1.0 ); } ) )
   return;

  v_Availability.assign( f_time_horizon , 1.0 );
 }

 if( ! ordered )
  std::sort( subset.begin() , subset.end() );

 if( subset.back() >= v_Availability.size() )
  throw( std::invalid_argument(
   "ThermalUnitBlock::set_availability: invalid index in subset." ) );

 // If nothing changes, return
 bool identical = true;
 auto availability = values;
 for( auto t : subset ) {
  if( v_Availability[ t ] != *availability )  // Check change
   identical = false;

  // Check consistency
  if( ! availability_is_consistent( t , *availability ) )
   throw( std::logic_error(
    "ThermalUnitBlock::set_availability: availability (" +
    std::to_string( *availability ) + ") at time " +
    std::to_string( t ) + " is not consistent." ) );

  std::advance( availability , 1 );
 }

 if( identical )
  return;  // nothing changes; return

 if( not_dry_run( issuePMod ) ) {
  // Change the physical representation
  availability = values;
  for( auto t : subset )
   v_Availability[ t ] = *(availability++);
 }

 if( not_dry_run( issueAMod ) && constraints_generated() )
  // Change the abstract representation
  for( auto t : subset )
   update_availability_dependents( t , un_ModBlock( issueAMod ) );

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockSbstMod >(
                            this , ThermalUnitBlockMod::eSetAv ,
                            std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_availability( subset ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_availability( MF_dbl_it values ,
                                         Range rng ,
                                         ModParam issuePMod ,
                                         ModParam issueAMod )
{
 rng.second = std::min( rng.second , f_time_horizon );
 if( rng.second <= rng.first )
  return;

 if( v_Availability.empty() ) {
  if( std::all_of( values ,
                   values + ( rng.second - rng.first ) ,
                   []( double cst ) { return( cst == 1.0 ); } ) )
   return;

  v_Availability.assign( f_time_horizon , 1.0 );
 }

 // If nothing changes, return
 if( std::equal( values ,
                 values + ( rng.second - rng.first ) ,
                 v_Availability.begin() + rng.first ) )
  return;

 // Check consistency
 auto availability = values;
 for( Index t = rng.first ; t < rng.second ; ++t ) {
  if( ! availability_is_consistent( t , *availability ) )
   throw( std::logic_error(
    "ThermalUnitBlock::set_availability: availability (" +
    std::to_string( *availability ) + ") at time " +
    std::to_string( t ) + " is not consistent." ) );

  std::advance( availability , 1 );
 }

 if( not_dry_run( issuePMod ) )
  // Change the physical representation
  std::copy( values ,
             values + ( rng.second - rng.first ) ,
             v_Availability.begin() + rng.first );

 if( not_dry_run( issueAMod ) && constraints_generated() )
  // Change the abstract representation
  for( Index t = rng.first ; t < rng.second ; ++t )
   update_availability_dependents( t , un_ModBlock( issueAMod ) );


 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockRngdMod >(
                            this , ThermalUnitBlockMod::eSetAv , rng ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_availability( range ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_maximum_power( MF_dbl_it values ,
                                          Subset && subset ,
                                          const bool ordered ,
                                          ModParam issuePMod ,
                                          ModParam issueAMod )
{
 if( subset.empty() )
  return;

 if( v_MaxPower.empty() ) {
  if( std::all_of( values ,
                   values + subset.size() ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  v_MaxPower.assign( f_time_horizon , 0 );
 }

 if( ! ordered )
  std::sort( subset.begin() , subset.end() );

 if( subset.back() >= v_MaxPower.size() )
  throw( std::invalid_argument(
   "ThermalUnitBlock::set_maximum_power: invalid index in subset." ) );

 if( identical( v_MaxPower , subset , values ) )  // if nothing changes
  return;                                         // return

 if( not_dry_run( issuePMod ) )
  // Change the physical representation
  assign( v_MaxPower , subset , values );

 if( not_dry_run( issueAMod ) && constraints_generated() )
  // Change the abstract representation
  for( auto t : subset )
   // the commitment variable is in position 0 in the LF
   static_cast< LinearFunction * >( MaxPower_Const[ t ].get_function()
   )->modify_coefficient( 0 , get_operational_max_power( t ) ,
                          un_ModBlock( issueAMod ) );

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockSbstMod >(
                            this , ThermalUnitBlockMod::eSetMaxP ,
                            std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_maximum_power( subset ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_maximum_power( MF_dbl_it values ,
                                          Range rng ,
                                          ModParam issuePMod ,
                                          ModParam issueAMod )
{
 rng.second = std::min( rng.second , f_time_horizon );
 if( rng.second <= rng.first )
  return;

 if( v_MaxPower.empty() ) {
  if( std::all_of( values ,
                   values + ( rng.second - rng.first ) ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  v_MaxPower.assign( f_time_horizon , 0 );
 }

 // If nothing changes, return
 if( std::equal( values ,
                 values + ( rng.second - rng.first ) ,
                 v_MaxPower.begin() + rng.first ) )
  return;


 if( not_dry_run( issuePMod ) )
  // Change the physical representation
  std::copy( values ,
             values + ( rng.second - rng.first ) ,
             v_MaxPower.begin() + rng.first );

 if( not_dry_run( issueAMod ) && constraints_generated() )
  // Change the abstract representation
  for( Index t = rng.first ; t < rng.second ; ++t )
   // the commitment variable is in position 0 in the LF
   static_cast< LinearFunction * >( MaxPower_Const[ t ].get_function()
   )->modify_coefficient( 0 , get_operational_max_power( t ) ,
                          un_ModBlock( issueAMod ) );

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockRngdMod >(
                            this , ThermalUnitBlockMod::eSetMaxP , rng ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_maximum_power( range ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::update_initial_power_in_cnstrs( ModParam issueAMod )
{
 if( ! ( RampUp_Const.empty() || v_DeltaRampUp.empty() ) )
  if( f_InitUpDownTime > 0 )
   RampUp_Const[ 0 ].set_rhs( v_DeltaRampUp[ 0 ] + f_InitialPower ,
                              issueAMod );

 if( ! RampDown_Const.empty() )
  if( f_InitUpDownTime > 0 )
   RampDown_Const[ 0 ].set_lhs( f_InitialPower , issueAMod );

}  // end( ThermalUnitBlock::update_initial_power_in_cnstrs )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_initial_power( MF_dbl_it values ,
                                          Subset && subset ,
                                          const bool ordered ,
                                          ModParam issuePMod ,
                                          ModParam issueAMod )
{
 if( subset.empty() )
  return;

 // Find the last index 0
 auto index_it = std::find( subset.rbegin() , subset.rend() , 0 );

 if( index_it == subset.rend() )
  return;  // 0 is not in subset; return

 std::advance( values , std::distance( index_it , subset.rend() ) - 1 );

 if( f_InitialPower == *values )
  return;  // nothing changes; return

 if( not_dry_run( issuePMod ) )
  // Change the physical representation
  f_InitialPower = *values;

 if( not_dry_run( issueAMod ) && constraints_generated() )
  // Change the abstract representation
  update_initial_power_in_cnstrs( un_ModBlock( issueAMod ) );

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockMod >(
                            this , ThermalUnitBlockMod::eSetInitP ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_initial_power( subset ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_initial_power( MF_dbl_it values ,
                                          Range rng ,
                                          ModParam issuePMod ,
                                          ModParam issueAMod )
{
 rng.second = std::min( rng.second , Index( 1 ) );
 if( ! ( ( rng.first <= 0 ) && ( 0 < rng.second ) ) )
  return;  // 0 does not belong to the range; return

 std::advance( values , -rng.first );

 if( f_InitialPower == *values )
  return;  // nothing changes; return

 if( not_dry_run( issuePMod ) )
  // Change the physical representation
  f_InitialPower = *values;

 if( not_dry_run( issueAMod ) && constraints_generated() )
  // Change the abstract representation
  update_initial_power_in_cnstrs( un_ModBlock( issueAMod ) );

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockMod >(
                            this , ThermalUnitBlockMod::eSetInitP ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_initial_power( range ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_startup_costs( MF_dbl_it values ,
                                          Subset && subset ,
                                          const bool ordered ,
                                          ModParam issuePMod ,
                                          ModParam issueAMod )
{
 if( subset.empty() )
  return;

 if( v_StartUpCost.empty() ) {
  if( std::all_of( values ,
                   values + subset.size() ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  v_StartUpCost.assign( f_time_horizon , 0 );
 }

 if( ! ordered )
  std::sort( subset.begin() , subset.end() );

 if( subset.back() >= v_StartUpCost.size() )
  throw( std::invalid_argument(
   "ThermalUnitBlock::set_startup_costs: invalid index in subset." ) );

 if( identical( v_StartUpCost , subset , values ) )  // if nothing changes
  return;                                            // return

 // the order of the variables in the Objective Function is:
 //
 // - first f_time_horizon - init_t start-up variables
 //
 // - then the rest ...
 //
 // this means that the start_up variable t is in position t - init_t
 // hence, those in the range [ 0 , init_t ) do not exist and their cost
 // cannot be changed
 if( subset.front() < init_t )
  throw( std::invalid_argument(
   "ThermalUnitBlock::set_startup_costs: invalid starting index in subset." ) );

 if( not_dry_run( issuePMod ) )
  // Change the physical representation
  assign( v_StartUpCost , subset , values );

 if( not_dry_run( issueAMod ) && objective_generated() ) {
  // Change the abstract representation
  Subset tmps = subset_sbtrct( subset , init_t );
  DQuadFunction::Vec_FunctionValue tmpv( values , values + subset.size() );
  static_cast< DQuadFunction * >( objective.get_function()
  )->modify_linear_coefficients( std::move( tmpv ) , std::move( tmps ) ,
                                 true , un_ModBlock( issueAMod ) );
 }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockSbstMod >(
                            this , ThermalUnitBlockMod::eSetSUC ,
                            std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_startup_costs( subset ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_startup_costs( MF_dbl_it values ,
                                          Range rng ,
                                          ModParam issuePMod ,
                                          ModParam issueAMod )
{
 rng.second = std::min( rng.second , f_time_horizon );
 if( rng.second <= rng.first )
  return;

 c_Index sz = rng.second - rng.first;
 if( v_StartUpCost.empty() ) {
  if( std::all_of( values ,
                   values + sz ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  v_StartUpCost.assign( f_time_horizon , 0 );
 }

 // If nothing changes, return
 if( std::equal( values , values + sz , v_StartUpCost.begin() + rng.first ) )
  return;

 // the order of the variables in the Objective Function is:
 //
 // - first f_time_horizon - init_t start-up variables
 //
 // - then the rest ...
 //
 // this means that the start_up variable t is in position t - init_t
 // hence, those in the range [ 0 , init_t ) do not exist and their cost
 // cannot be changed
 if( rng.first < init_t )
  throw( std::invalid_argument( "ThermalUnitBlock::set_startup_costs: invalid"
                                " starting index in range." ) );

 if( not_dry_run( issuePMod ) )
  // Change the physical representation
  std::copy( values ,
             values + sz ,
             v_StartUpCost.begin() + rng.first );

 if( not_dry_run( issueAMod ) && objective_generated() ) {
  // Change the abstract representation
  DQuadFunction::Vec_FunctionValue tmpv( values , values + sz );
  static_cast< DQuadFunction * >( objective.get_function()
  )->modify_linear_coefficients( std::move( tmpv ) ,
                                 Range( rng.first - init_t ,
                                        rng.second - init_t ) ,
                                 un_ModBlock( issueAMod ) );
 }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockRngdMod >(
                            this , ThermalUnitBlockMod::eSetSUC , rng ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_startup_costs( range ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_const_term( MF_dbl_it values ,
                                       Subset && subset ,
                                       const bool ordered ,
                                       ModParam issuePMod ,
                                       ModParam issueAMod )
{
 if( subset.empty() )
  return;

 if( v_ConstTerm.empty() ) {
  if( std::all_of( values ,
                   values + subset.size() ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  v_ConstTerm.assign( f_time_horizon , 0 );
 }

 if( ! ordered )
  std::sort( subset.begin() , subset.end() );

 if( subset.back() >= v_ConstTerm.size() )
  throw( std::invalid_argument(
   "ThermalUnitBlock::set_const_term: invalid index in subset." ) );

 if( identical( v_ConstTerm , subset , values ) )  // if nothing changes
  return;                                          // return

 if( not_dry_run( issuePMod ) )
  // Change the physical representation
  assign( v_ConstTerm , subset , values );

 if( not_dry_run( issueAMod ) && objective_generated() ) {
  // Change the abstract representation
  // the order of the variables in the Objective Function is:
  //
  // - first f_time_horizon - init_t start-up variables
  //
  // - then f_time_horizon active power variables
  //
  // - then f_time_horizon commitment variables
  //
  // - then possibly the rest
  //
  // hence, the commitment variables, whose coefficient is the fixed
  // cost, start from position 2 * f_time_horizon - init_t
  const Index dpos = 2 * f_time_horizon - init_t;

  Subset tmps = subset_add( subset , dpos );
  DQuadFunction::Vec_FunctionValue tmpv( values , values + subset.size() );
  static_cast< DQuadFunction * >( objective.get_function()
  )->modify_linear_coefficients( std::move( tmpv ) , std::move( tmps ) ,
                                 true , un_ModBlock( issueAMod ) );
 }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockSbstMod >(
                            this , ThermalUnitBlockMod::eSetConstT ,
                            std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_const_term( subset ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_const_term( MF_dbl_it values ,
                                       Range rng ,
                                       ModParam issuePMod ,
                                       ModParam issueAMod )
{
 rng.second = std::min( rng.second , f_time_horizon );
 if( rng.second <= rng.first )
  return;

 c_Index sz = rng.second - rng.first;
 if( v_ConstTerm.empty() ) {
  if( std::all_of( values ,
                   values + sz ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  v_ConstTerm.assign( f_time_horizon , 0 );
 }

 // If nothing changes, return
 if( std::equal( values , values + sz , v_ConstTerm.begin() + rng.first ) )
  return;

 if( not_dry_run( issuePMod ) )
  // Change the physical representation
  std::copy( values ,
             values + sz ,
             v_ConstTerm.begin() + rng.first );

 if( not_dry_run( issueAMod ) && objective_generated() ) {
  // Change the abstract representation
  // the order of the variables in the Objective Function is:
  //
  // - first f_time_horizon - init_t start-up variables
  //
  // - then f_time_horizon active power variables
  //
  // - then f_time_horizon commitment variables
  //
  // - then possibly the rest
  //
  // hence, the commitment variables, whose coefficient is the fixed
  // cost, start from position 2 * f_time_horizon - init_t
  const Index dpos = 2 * f_time_horizon - init_t;

  DQuadFunction::Vec_FunctionValue tmpv( values , values + sz );
  static_cast< DQuadFunction * >( objective.get_function()
  )->modify_linear_coefficients( std::move( tmpv ) ,
                                 Range( rng.first + dpos ,
                                        rng.second + dpos ) ,
                                 un_ModBlock( issueAMod ) );
 }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockRngdMod >(
                            this , ThermalUnitBlockMod::eSetConstT , rng ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_const_term( range ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_linear_term( MF_dbl_it values ,
                                        Subset && subset ,
                                        const bool ordered ,
                                        ModParam issuePMod ,
                                        ModParam issueAMod )
{
 if( subset.empty() )
  return;

 if( v_LinearTerm.empty() ) {
  if( std::all_of( values ,
                   values + subset.size() ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  v_LinearTerm.assign( f_time_horizon , 0 );
 }

 if( ! ordered )
  std::sort( subset.begin() , subset.end() );

 if( subset.back() >= v_LinearTerm.size() )
  throw( std::invalid_argument(
   "ThermalUnitBlock::set_linear_term: invalid index in subset." ) );

 if( identical( v_LinearTerm , subset , values ) )  // if nothing changes
  return;                                           // return

 if( not_dry_run( issuePMod ) )
  // Change the physical representation
  assign( v_LinearTerm , subset , values );

 if( not_dry_run( issueAMod ) && objective_generated() ) {
  // Change the abstract representation
  // the order of the variables in the Objective Function is:
  //
  // - first f_time_horizon - init_t start-up variables
  //
  // - then f_time_horizon active power variables
  //
  // - then possibly the rest
  //
  // hence, the active power  variables, whose coefficient is the linear
  // term of the cost, start from position f_time_horizon - init_t
  const Index dpos = f_time_horizon - init_t;

  Subset tmps = subset_add( subset , dpos );
  DQuadFunction::Vec_FunctionValue tmpv( values , values + subset.size() );
  static_cast< DQuadFunction * >( objective.get_function()
  )->modify_linear_coefficients( std::move( tmpv ) , std::move( tmps ) ,
                                 true , un_ModBlock( issueAMod ) );
 }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockSbstMod >(
                            this , ThermalUnitBlockMod::eSetLinT ,
                            std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_linear_term( subset ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_linear_term( MF_dbl_it values ,
                                        Range rng ,
                                        ModParam issuePMod ,
                                        ModParam issueAMod )
{
 rng.second = std::min( rng.second , f_time_horizon );
 if( rng.second <= rng.first )
  return;

 c_Index sz = rng.second - rng.first;
 if( v_LinearTerm.empty() ) {
  if( std::all_of( values ,
                   values + sz ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  v_LinearTerm.assign( f_time_horizon , 0 );
 }

 // If nothing changes, return
 if( std::equal( values , values + sz , v_LinearTerm.begin() + rng.first ) )
  return;

 if( not_dry_run( issuePMod ) )
  // Change the physical representation
  std::copy( values , values + sz , v_LinearTerm.begin() + rng.first );

 if( not_dry_run( issueAMod ) && objective_generated() ) {
  // Change the abstract representation
  // the order of the variables in the Objective Function is:
  //
  // - first f_time_horizon - init_t start-up variables
  //
  // - then f_time_horizon active power variables
  //
  // - then possibly the rest
  //
  // hence, the active power variables, whose coefficient is the linear
  // term of the cost, start from position f_time_horizon - init_t
  const Index dpos = f_time_horizon - init_t;

  DQuadFunction::Vec_FunctionValue tmpv( values , values + sz );
  static_cast< DQuadFunction * >( objective.get_function()
  )->modify_linear_coefficients( std::move( tmpv ) ,
                                 Range( rng.first + dpos ,
                                        rng.second + dpos ) ,
                                 un_ModBlock( issueAMod ) );
 }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockRngdMod >(
                            this , ThermalUnitBlockMod::eSetLinT , rng ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_linear_term( range ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_quad_term( MF_dbl_it values ,
                                      Subset && subset ,
                                      const bool ordered ,
                                      ModParam issuePMod ,
                                      ModParam issueAMod )
{
 if( subset.empty() )
  return;

 if( v_QuadTerm.empty() ) {
  if( std::all_of( values ,
                   values + subset.size() ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  v_QuadTerm.assign( f_time_horizon , 0 );
 }

 if( ! ordered )
  std::sort( subset.begin() , subset.end() );

 if( subset.back() >= v_QuadTerm.size() )
  throw( std::invalid_argument(
   "ThermalUnitBlock::set_quad_term: invalid index in subset." ) );

 if( identical( v_QuadTerm , subset , values ) )  // if nothing changes
  return;                                         // return

 if( not_dry_run( issuePMod ) )
  // Change the physical representation
  assign( v_QuadTerm , subset , values );

 if( not_dry_run( issueAMod ) && objective_generated() ) {
  // Change the abstract representation
  // the order of the variables in the Objective Function is:
  //
  // - first f_time_horizon - init_t start-up variables
  //
  // - then f_time_horizon active power variables (which may have the
  //   nonzero quadratic cost coefficient, while the others do not)
  //
  // - then possibly the rest
  //
  // hence, if no perspective cuts are used, then the active power variables,
  // whose quadratic coefficient is the quadratic term of the cost, start
  // from position f_time_horizon - init_t, else the quadratic coefficient
  // become the quadratic term of the perspective cut variables, start from
  // position
  // 5 * f_time_horizon - init_t if both primary and secondary reserve are
  // defined, and
  // 4 * f_time_horizon - init_t if just one between primary or secondary
  // reverse is defined, and
  // 3 * f_time_horizon - init_t otherwise
  const Index dpos = ! ( AR & PCuts ) ? f_time_horizon - init_t :
                     ( ( ( ( ( ! v_primary_spinning_reserve.empty() ) &&
                             ( reserve_vars & 1u ) ) &&
                           ( ( ! v_secondary_spinning_reserve.empty() ) &&
                             ( reserve_vars & 2u ) ) ) ? 5 :
                         ( ( ( ( ! v_primary_spinning_reserve.empty() ) &&
                               ( reserve_vars & 1u ) ) ||
                             ( ( ! v_secondary_spinning_reserve.empty() ) &&
                               ( reserve_vars & 2u ) ) ) ? 4 : 3 ) ) *
                       f_time_horizon - init_t );

  Subset tmps = subset_add( subset , dpos );

  if( ! ( AR & PCuts ) ) {

   DQuadFunction::Vec_FunctionValue tmplv( subset.size() , 0 );
   if( ! v_LinearTerm.empty() ) {
    auto tmplvit = tmplv.begin();
    for( auto t : subset )
     *(tmplvit++) = v_LinearTerm[ t ];
   }

   static_cast< DQuadFunction * >( objective.get_function()
   )->modify_terms( values , tmplv.begin() , std::move( tmps ) , true ,
                    un_ModBlock( issueAMod ) );

  } else {

   DQuadFunction::Vec_FunctionValue tmplv( values , values + subset.size() );
   static_cast< DQuadFunction * >( objective.get_function()
   )->modify_linear_coefficients( std::move( tmplv ) , std::move( tmps ) ,
                                  true , un_ModBlock( issueAMod ) );
  }
 }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockSbstMod >(
                            this , ThermalUnitBlockMod::eSetQuadT ,
                            std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_quad_term( subset ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_quad_term( MF_dbl_it values ,
                                      Range rng ,
                                      ModParam issuePMod ,
                                      ModParam issueAMod )
{
 rng.second = std::min( rng.second , f_time_horizon );
 if( rng.second <= rng.first )
  return;

 c_Index sz = rng.second - rng.first;
 if( v_QuadTerm.empty() ) {
  if( std::all_of( values ,
                   values + sz ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  v_QuadTerm.assign( f_time_horizon , 0 );
 }

 // If nothing changes, return
 if( std::equal( values , values + sz , v_QuadTerm.begin() + rng.first ) )
  return;

 if( not_dry_run( issuePMod ) )
  // Change the physical representation
  std::copy( values , values + sz , v_QuadTerm.begin() + rng.first );

 if( not_dry_run( issueAMod ) && objective_generated() ) {
  // Change the abstract representation
  // the order of the variables in the Objective Function is:
  //
  // - first f_time_horizon - init_t start-up variables
  //
  // - then f_time_horizon active power variables (which may have the
  //   nonzero quadratic cost coefficient, while the others do not)
  //
  // - then possibly the rest
  //
  // hence, if no perspective cuts are used, then the active power variables,
  // whose quadratic coefficient is the quadratic term of the cost, start
  // from position f_time_horizon - init_t, else the quadratic coefficient
  // become the quadratic term of the perspective cut variables, start from
  // position
  // 5 * f_time_horizon - init_t if both primary and secondary reserve are
  // defined, and
  // 4 * f_time_horizon - init_t if just one between primary or secondary
  // reverse is defined, and
  // 3 * f_time_horizon - init_t otherwise
  const Index dpos = ! ( AR & PCuts ) ? f_time_horizon - init_t :
                     ( ( ( ( ( ! v_primary_spinning_reserve.empty() ) &&
                             ( reserve_vars & 1u ) ) &&
                           ( ( ! v_secondary_spinning_reserve.empty() ) &&
                             ( reserve_vars & 2u ) ) ) ? 5 :
                         ( ( ( ( ! v_primary_spinning_reserve.empty() ) &&
                               ( reserve_vars & 1u ) ) ||
                             ( ( ! v_secondary_spinning_reserve.empty() ) &&
                               ( reserve_vars & 2u ) ) ) ? 4 : 3 ) ) *
                       f_time_horizon - init_t );

  if( ! ( AR & PCuts ) ) {

   DQuadFunction::Vec_FunctionValue tmplv( sz , 0 );
   if( ! v_LinearTerm.empty() )
    std::copy( v_LinearTerm.begin() + rng.first ,
               v_LinearTerm.begin() + rng.second , tmplv.begin() );

   static_cast< DQuadFunction * >( objective.get_function()
   )->modify_terms( values , tmplv.begin() ,
                    Range( rng.first + dpos , rng.second + dpos ) ,
                    un_ModBlock( issueAMod ) );

  } else {

   DQuadFunction::Vec_FunctionValue tmplv( values , values + sz );
   static_cast< DQuadFunction * >( objective.get_function()
   )->modify_linear_coefficients( std::move( tmplv ) ,
                                  Range( rng.first + dpos ,
                                         rng.second + dpos ) ,
                                  un_ModBlock( issueAMod ) );
  }
 }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockRngdMod >(
                            this , ThermalUnitBlockMod::eSetQuadT , rng ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_quad_term( range ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_primary_spinning_reserve_cost( MF_dbl_it values ,
                                                          Subset && subset ,
                                                          const bool ordered ,
                                                          ModParam issuePMod ,
                                                          ModParam issueAMod )
{
 if( v_primary_spinning_reserve.empty() || ( ! ( reserve_vars & 1u ) ) )
  return;  // primary reserve is not there, silently return

 if( subset.empty() )
  return;

 if( v_PrimarySpinningReserveCost.empty() ) {
  // The primary spinning reserve costs are currently all zero.
  if( std::all_of( values ,
                   values + subset.size() ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;  // The given values are zero: nothing to do

  v_PrimarySpinningReserveCost.assign( f_time_horizon , 0 );
 }

 if( ! ordered )
  std::sort( subset.begin() , subset.end() );

 if( subset.back() >= v_PrimarySpinningReserveCost.size() )
  throw( std::invalid_argument(
   "ThermalUnitBlock::set_primary_spinning_reserve_cost: "
   "invalid index in subset." ) );

 if( identical( v_PrimarySpinningReserveCost , subset , values ) )
  return;

 if( not_dry_run( issuePMod ) )
  // Change the physical representation
  assign( v_PrimarySpinningReserveCost , subset , values );

 if( not_dry_run( issueAMod ) && objective_generated() ) {
  // Change the abstract representation
  // the order of the variables in the Objective Function is:
  //
  // - first f_time_horizon - init_t start-up variables
  //
  // - then f_time_horizon active power variables
  //
  // - then f_time_horizon commitment variables
  //
  // - then possibly f_time_horizon primary reserve variables
  //
  // - then possibly the rest
  //
  // hence, the primary reserve variables, whose linear coefficient is the
  // primary spinning reserve cost, start from position
  // 3 * f_time_horizon - init_t
  const Index dpos = 3 * f_time_horizon - init_t;

  Subset tmps = subset_add( subset , dpos );
  DQuadFunction::Vec_FunctionValue tmpv( values , values + subset.size() );
  static_cast< DQuadFunction * >( objective.get_function()
  )->modify_linear_coefficients( std::move( tmpv ) , std::move( tmps ) ,
                                 true , un_ModBlock( issueAMod ) );
 }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockSbstMod >(
                            this , ThermalUnitBlockMod::eSetPrSpResCost ,
                            std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_primary_spinning_reserve_cost( subset ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_primary_spinning_reserve_cost( MF_dbl_it values ,
                                                          Range rng ,
                                                          ModParam issuePMod ,
                                                          ModParam issueAMod )
{
 if( v_primary_spinning_reserve.empty() || ( ! ( reserve_vars & 1u ) ) )
  return;  // primary reserve is not there, silently return

 rng.second = std::min( rng.second , f_time_horizon );
 if( rng.second <= rng.first )
  return;  // Empty range. Return.

 c_Index sz = rng.second - rng.first;
 if( v_PrimarySpinningReserveCost.empty() ) {
  // The primary spinning reserve costs are currently all zero.
  if( std::all_of( values ,
                   values + sz ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;  // The given values are zero. So, there is nothing to be changed.

  v_PrimarySpinningReserveCost.assign( f_time_horizon , 0 );
 }

 // If nothing changes, return
 if( std::equal( values ,
                 values + sz ,
                 v_PrimarySpinningReserveCost.begin() + rng.first ) )
  return;

 if( not_dry_run( issuePMod ) )
  // Change the physical representation
  std::copy( values ,
             values + sz ,
             v_PrimarySpinningReserveCost.begin() + rng.first );

 if( not_dry_run( issueAMod ) && objective_generated() ) {
  // Change the abstract representation
  // the order of the variables in the Objective Function is:
  //
  // - first f_time_horizon - init_t start-up variables
  //
  // - then f_time_horizon active power variables
  //
  // - then f_time_horizon commitment variables
  //
  // - then possibly f_time_horizon primary reserve variables
  //
  // - then possibly the rest
  //
  // hence, the primary reserve variables, whose linear coefficient is the
  // primary spinning reserve cost, start from position
  // 3 * f_time_horizon - init_t
  const Index dpos = 3 * f_time_horizon - init_t;

  DQuadFunction::Vec_FunctionValue tmpv( values , values + sz );
  static_cast< DQuadFunction * >( objective.get_function()
  )->modify_linear_coefficients( std::move( tmpv ) ,
                                 Range( rng.first + dpos ,
                                        rng.second + dpos ) ,
                                 un_ModBlock( issueAMod ) );
 }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockRngdMod >(
                            this , ThermalUnitBlockMod::eSetPrSpResCost , rng ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_primary_spinning_reserve_cost( range ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_secondary_spinning_reserve_cost( MF_dbl_it values ,
                                                            Subset && subset ,
                                                            const bool ordered ,
                                                            ModParam issuePMod ,
                                                            ModParam issueAMod )
{
 if( v_secondary_spinning_reserve.empty() || ( ! ( reserve_vars & 2u ) ) )
  return;  // secondary reserve is not there, silently return

 if( subset.empty() )
  return;

 if( v_SecondarySpinningReserveCost.empty() ) {
  // The secondary spinning reserve costs are currently all zero.
  if( std::all_of( values ,
                   values + subset.size() ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;  // The given values are zero: nothing to do

  v_SecondarySpinningReserveCost.assign( f_time_horizon , 0 );
 }

 if( ! ordered )
  std::sort( subset.begin() , subset.end() );

 if( subset.back() >= v_SecondarySpinningReserveCost.size() )
  throw( std::invalid_argument( "ThermalUnitBlock::set_secondary_spinning_"
                                 "reserve_cost: invalid index in subset." ) );

 if( identical( v_SecondarySpinningReserveCost , subset , values ) )
  return;

 if( not_dry_run( issuePMod ) )
  // Change the physical representation
  assign( v_SecondarySpinningReserveCost , subset , values );

 if( not_dry_run( issueAMod ) && objective_generated() ) {
  // Change the abstract representation
  // the order of the variables in the Objective Function is:
  //
  // - first f_time_horizon - init_t start-up variables
  //
  // - then f_time_horizon active power variables
  //
  // - then f_time_horizon commitment variables
  //
  // - then possibly f_time_horizon primary reserve variables
  //
  // - then possibly f_time_horizon secondary reserve variables
  //
  // - then possibly the rest
  //
  // hence, the secondary reserve variables, whose linear coefficient is the
  // secondary spinning reserve cost, start from position
  // 4 * f_time_horizon - init_t if primary reserve is defined, and
  // 3 * f_time_horizon - init_t otherwise
  const Index dpos = ( v_primary_spinning_reserve.empty() ||
                       ( ! ( reserve_vars & 1u ) ) ? 3 : 4 ) *
                     f_time_horizon - init_t;

  Subset tmps = subset_add( subset , dpos );
  DQuadFunction::Vec_FunctionValue tmpv( values , values + subset.size() );
  static_cast< DQuadFunction * >( objective.get_function()
  )->modify_linear_coefficients( std::move( tmpv ) , std::move( tmps ) ,
                                 true , un_ModBlock( issueAMod ) );
 }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockSbstMod >(
                            this , ThermalUnitBlockMod::eSetSecSpResCost ,
                            std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_secondary_spinning_reserve_cost( subset ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_secondary_spinning_reserve_cost( MF_dbl_it values ,
                                                            Range rng ,
                                                            ModParam issuePMod ,
                                                            ModParam issueAMod )
{
 if( v_secondary_spinning_reserve.empty() || ( ! ( reserve_vars & 2u ) ) )
  return;  // secondary reserve is not there, silently return

 rng.second = std::min( rng.second , f_time_horizon );
 if( rng.second <= rng.first )
  return;  // Empty range. Return.

 c_Index sz = rng.second - rng.first;
 if( v_SecondarySpinningReserveCost.empty() ) {
  // The secondary spinning reserve costs are currently all zero.
  if( std::all_of( values ,
                   values + sz ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;  // The given values are zero. So, there is nothing to be changed.

  v_SecondarySpinningReserveCost.assign( f_time_horizon , 0 );
 }

 // If nothing changes, return
 if( std::equal( values ,
                 values + sz ,
                 v_SecondarySpinningReserveCost.begin() + rng.first ) )
  return;

 if( not_dry_run( issuePMod ) )
  // Change the physical representation
  std::copy( values ,
             values + sz ,
             v_SecondarySpinningReserveCost.begin() + rng.first );

 if( not_dry_run( issueAMod ) && objective_generated() ) {
  // Change the abstract representation
  // the order of the variables in the Objective Function is:
  //
  // - first f_time_horizon - init_t start-up variables
  //
  // - then f_time_horizon active power variables
  //
  // - then f_time_horizon commitment variables
  //
  // - then possibly f_time_horizon primary reserve variables
  //
  // - then possibly f_time_horizon secondary reserve variables
  //
  // - then possibly the rest
  //
  // hence, the secondary reserve variables, whose linear coefficient is the
  // secondary spinning reserve cost, start from position
  // 4 * f_time_horizon - init_t if primary reserve is defined, and
  // 3 * f_time_horizon - init_t otherwise
  const Index dpos = ( v_primary_spinning_reserve.empty() ||
                       ( ! ( reserve_vars & 1u ) ) ? 3 : 4 ) *
                     f_time_horizon - init_t;

  DQuadFunction::Vec_FunctionValue tmpv( values , values + sz );
  static_cast< DQuadFunction * >( objective.get_function()
  )->modify_linear_coefficients( std::move( tmpv ) ,
                                 Range( rng.first + dpos ,
                                        rng.second + dpos ) ,
                                 un_ModBlock( issueAMod ) );
 }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockRngdMod >(
                            this , ThermalUnitBlockMod::eSetSecSpResCost , rng ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_secondary_spinning_reserve_cost( range ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_init_updown_time( MF_int_it values ,
                                             Subset && subset ,
                                             const bool ordered ,
                                             ModParam issuePMod ,
                                             ModParam issueAMod )
{
 if( subset.empty() )
  return;

 // Find the last index 0
 auto index_it = std::find( subset.rbegin() , subset.rend() , 0 );

 if( index_it == subset.rend() )
  return;  // 0 is not in subset; return

 std::advance( values , std::distance( index_it , subset.rend() ) - 1 );

 if( f_InitUpDownTime == *values )
  return;  // nothing changes; return

 if( not_dry_run( issuePMod ) )
  // Change the physical representation
  f_InitUpDownTime = *values;

 if( not_dry_run( issueAMod ) && variables_generated() )  // TODO
  throw( std::logic_error( "ThermalUnitBlock::set_init_updown_time: it is "
                            "currently not possible to update the abstract "
                            "representation." ) );

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockMod >(
                            this , ThermalUnitBlockMod::eSetInitUD ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_init_updown_time( subset ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_init_updown_time( MF_int_it values ,
                                             Range rng ,
                                             ModParam issuePMod ,
                                             ModParam issueAMod )
{
 rng.second = std::min( rng.second , decltype( rng.second )( 1 ) );
 if( ! ( ( rng.first <= 0 ) && ( 0 < rng.second ) ) )
  return;  // 0 does not belong to the range; return

 std::advance( values , -rng.first );

 if( f_InitUpDownTime == *values )
  return;  // nothing changes; return

 if( not_dry_run( issuePMod ) )
  // Change the physical representation
  f_InitUpDownTime = *values;

 if( not_dry_run( issueAMod ) && variables_generated() )  // TODO
  throw( std::logic_error( "ThermalUnitBlock::set_init_updown_time: it is "
                            "currently not possible to update the abstract "
                            "representation." ) );

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< ThermalUnitBlockMod >(
                            this , ThermalUnitBlockMod::eSetInitUD ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_init_updown_time( range ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::scale( MF_dbl_it values ,
                              Subset && subset ,
                              const bool ordered ,
                              c_ModParam issuePMod ,
                              c_ModParam issueAMod )
{
 if( subset.empty() )
  return;  // Since the given Subset is empty, no operation is performed

 if( f_scale == *values )
  return;  // The scale factor does not change: nothing to do

 if( not_dry_run( issuePMod ) ) {
  f_scale = *values;  // Update the scale factor

  if( not_dry_run( issueAMod ) ) {
   // Update the abstract representation
   if( objective_generated() )
    // Update the Objective
    update_objective( Range( 0 , Inf< Index >() ) , issueAMod );
  }
 }

 if( issue_pmod( issuePMod ) )
  // Issue a Physical Modification
  Block::add_Modification( std::make_shared< UnitBlockMod >(
                            this , UnitBlockMod::eScale ) ,
                           Observer::par2chnl( issuePMod ) );
 else if( auto f_Block = get_f_Block() )
  f_Block->add_Modification( std::make_shared< UnitBlockMod >(
                              this , UnitBlockMod::eScale ) ,
                             Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::scale )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::update_objective_start_up( const Subset & subset ,
                                                  c_ModParam issueAMod )
{
 if( ! objective_generated() )
  return;  // the Objective has not been generated: nothing to be done

 auto function = dynamic_cast< DQuadFunction * >( objective.get_function() );

 if( ! function )
  return;

 for( auto t : subset ) {
  if( t < init_t )
   continue;

  auto var_index = function->is_active( &v_start_up[ t - init_t ] );
  assert( var_index < function->get_num_active_var() );
  function->modify_linear_coefficient( var_index ,
                                       f_scale * v_StartUpCost[ t ] ,
                                       issueAMod );
 }
}  // end( ThermalUnitBlock::update_objective_start_up )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::update_objective_active_power( const Subset & subset ,
                                                      c_ModParam issueAMod )
{
 if( ! objective_generated() )
  return;  // the Objective has not been generated: nothing to be done

 auto function = dynamic_cast< DQuadFunction * >( objective.get_function() );

 if( ! function )
  return;

 for( auto t : subset ) {
  auto var_index = function->is_active( &v_active_power[ t ] );
  assert( var_index < function->get_num_active_var() );
  function->modify_term( var_index ,
                         f_scale * v_LinearTerm[ t ] ,
                         AR & PCuts ? 0.0 : f_scale * v_QuadTerm[ t ] ,
                         issueAMod );
 }
}  // end( ThermalUnitBlock::update_objective_active_power )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::update_objective_commitment( const Subset & subset ,
                                                    c_ModParam issueAMod )
{
 if( ! objective_generated() )
  return;  // the Objective has not been generated: nothing to be done

 auto function = dynamic_cast< DQuadFunction * >( objective.get_function() );

 if( ! function )
  return;

 for( auto t : subset ) {
  auto var_index = function->is_active( &v_commitment[ t ] );
  assert( var_index < function->get_num_active_var() );
  function->modify_linear_coefficient( var_index ,
                                       f_scale * v_ConstTerm[ t ] ,
                                       issueAMod );
 }
}  // end( ThermalUnitBlock::update_objective_commitment )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::update_objective( const Subset & subset ,
                                         c_ModParam issueAMod )
{
 update_objective_start_up( subset , issueAMod );
 update_objective_active_power( subset , issueAMod );
 update_objective_commitment( subset , issueAMod );
}  // end( ThermalUnitBlock::update_objective( subset ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::update_objective( Range rng ,
                                         c_ModParam issueAMod )
{
 rng.second = std::min( rng.second , f_time_horizon );
 if( rng.second <= rng.first )
  return;

 Subset subset( rng.second - rng.first );
 std::iota( subset.begin() , subset.end() , rng.first );

 update_objective( subset , issueAMod );
}  // end( ThermalUnitBlock::update_objective( range ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::guts_of_add_Modification( p_Mod mod , ChnlName chnl )
{
 // process abstract Modification - - - - - - - - - - - - - - - - - - - - - -
 /* This requires to patiently sift through the possible Modification types
  * to find what this Modification exactly is and appropriately mirror the
  * changes to the "abstract representation" to the "physical one".
  *
  * Note that since ThermalUnitBlock is a "leaf" Block (has no sub-Block),
  * this method does not have to deal with GroupModification since these
  * are produced by Block::add_Modification(), but this method is called
  * *before* that one is.
  *
  * As an important consequence,
  *
  *   THE STATE OF THE DATA STRUCTURE IN ThermalUnitBlock WHEN THIS METHOD
  *   IS EXECUTED IS PRECISELY THE ONE IN WHICH THE Modification WAS
  *   ISSUED: NO COMPLICATED OPERATIONS (Variable AND/OR Constraint BEING
  *   ADDED/REMOVED ...) CAN HAVE BEEN PERFORMED IN THE MEANTIME
  *
  * This assumption drastically simplifies some logic here. */

 // VariableMod - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( const auto tmod = dynamic_cast< VariableMod * >( mod ) ) {
  // changing the Variable is not supported, but the Modification is issued
  // when they are first generated, in which case it must be ignored
  if( ! variables_generated() )
   return;

  throw( std::logic_error( "ThermalUnitBlock - VariableMod not supported." ) );
  /*
  auto v = dynamic_cast< ColVariable * const >( tmod->variable() );

  if( v->is_fixed() ) {
   // TODO: Do something to the physical representation

   } else {
   // TODO: Do something to the physical representation
   }
  */
  return;
 }

 // BlockMod - Generic modification - - - - - - - - - - - - - - - - - - - - -
 if( const auto tmod = dynamic_cast< BlockMod * >( mod ) ) {
  // changing the Objective is not supported, but the Modification is issued
  // when it is first set, in which case it must be ignored
  if( ! objective_generated() )
   return;

  throw( std::logic_error( "ThermalUnitBlock - BlockMod not supported." ) );

  // TODO: BlockMod - obj changed
  return;
 }

 // FunctionMod - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( const auto tmod = dynamic_cast< FunctionMod * >( mod ) ) {
  auto f = tmod->function();
  if( f == static_cast< FRealObjective * >( get_objective() )->get_function() ) {
   handle_objective_change( tmod , chnl );
   return;
  }

  std::ostringstream em;
  em << *mod;
  throw( std::invalid_argument(
   "ThermalUnitBlock::guts_of_add_Modification: unsupported " + em.str() + "." ) );
  return;
 }

 // any other Modification is not supported - - - - - - - - - - - - - - - - -

 std::ostringstream em;
 em << *mod;
 throw( std::invalid_argument(
  "ThermalUnitBlock::guts_of_add_Modification: unsupported " + em.str() + "." ) );

}  // end( ThermalUnitBlock::guts_of_add_Modification )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::handle_objective_change( FunctionMod * mod ,
                                                ChnlName chnl )
{
 const auto * qf = static_cast< const DQuadFunction * >( mod->function() );
 auto par = make_par( eNoBlck , chnl );
 Index th = f_time_horizon;

 // C05FunctionModLinRngd - - - - - - - - - - - - - - - - - - - - - - - - - -
 // split the C05FunctionModLinRngd in up to 5 physical Modification by
 // calling the appropriate set_*() methods (ranged version) for those among
 // startup, power, commitment, primary/secondary reserve variables whose
 // coefficient change. This heavily relies on the fact that variables of
 // the same type are consecutive (and ordered in the obvious way) when
 // set as coefficients in the Objective

 if( const auto tmod = dynamic_cast< C05FunctionModLinRngd * >( mod ) ) {

  Index l = tmod->range().first;
  Index r = tmod->range().second;

  if( tmod->range().second > qf->get_num_active_var() )
   throw( std::invalid_argument(
    "ThermalUnitBlock::handle_objective_change: invalid Range [" +
    std::to_string( l ) + ", " + std::to_string( r ) +
    ") in C05FunctionModLinRngd." ) );

  std::vector< double > nv( r - l + 1 );
  Index gl = 0;
  Index gr = th - init_t;

  if( l < gr ) {  // startup variables
   Index r2 = std::min( r , gr );
   auto nvit = nv.begin();
   for( Index i = l ; i < r2 ; )
    *(nvit++) = qf->get_linear_coefficient( i++ );
   set_startup_costs( nv.begin() , Range( l , r2 ) , par , eDryRun );
   l = r2;
   if( l == r )
    return;
  }

  gl = gr;
  gr = 2 * th - init_t;

  if( l < gr ) {  // active power variables
   Index r2 = std::min( r , gr );
   auto nvit = nv.begin();
   for( Index i = l ; i < r2 ; )
    *(nvit++) = qf->get_linear_coefficient( i++ );
   set_linear_term( nv.begin() , Range( l - gl , r2 - gl ) , par , eDryRun );
   l = r2;
   if( l == r )
    return;
  }

  gl = gr;
  gr = 3 * th - init_t;

  if( l < gr ) {  // commitment variables
   Index r2 = std::min( r , gr );
   auto nvit = nv.begin();
   for( Index i = l ; i < r2 ; )
    *(nvit++) = qf->get_linear_coefficient( i++ );
   set_const_term( nv.begin() , Range( l - gl , r2 - gl ) , par , eDryRun );
   l = r2;
   if( l == r )
    return;
  }

  gl = gr;
  gr = 4 * th - init_t;

  if( l < gr ) {  // primary spinning reserve variables
   Index r2 = std::min( r , gr );
   auto nvit = nv.begin();
   for( Index i = l ; i < r2 ; )
    *(nvit++) = qf->get_linear_coefficient( i++ );
   set_primary_spinning_reserve_cost( nv.begin() , Range( l - gl , r2 - gl ) ,
                                      par , eDryRun );
   l = r2;
   if( l == r )
    return;
  }

  gl = gr;
  gr = 5 * th - init_t;

  if( l < gr ) {  // secondary spinning reserve variables
   Index r2 = std::min( r , gr );
   auto nvit = nv.begin();
   for( Index i = l ; i < r2 ; )
    *(nvit++) = qf->get_linear_coefficient( i++ );
   set_secondary_spinning_reserve_cost( nv.begin() , Range( l - gl , r2 - gl ) ,
                                        par , eDryRun );
   if( r2 == r )
    return;
  }

  throw( std::invalid_argument(
   "ThermalUnitBlock::handle_objective_change: invalid variable in "
   "C05FunctionModLinRngd." ) );
  return;

 }  // end( C05FunctionModLinRngd )

 // C05FunctionModLinSbst - - - - - - - - - - - - - - - - - - - - - - - - - -
 // split the C05FunctionModLinSbst in up to 5 physical Modification by
 // calling the appropriate set_*() methods (subset version) for those among
 // startup, power, commitment, primary/secondary reserve variables whose
 // coefficient change. This heavily relies on the fact that variables of
 // the same type are consecutive (and ordered in the obvious way) when
 // set as coefficients in the Objective

 if( const auto tmod = dynamic_cast< C05FunctionModLinSbst * >( mod ) ) {

  if( tmod->subset().back() > qf->get_num_active_var() )
   throw( std::invalid_argument(
    "ThermalUnitBlock::handle_objective_change: invalid Subset in "
    "C05FunctionModLinSbst." ) );

  std::vector< double > nv( tmod->subset().size() );
  auto l = tmod->subset().begin();
  Index gl = 0;
  Index gr = th - init_t;

  if( *l < gr ) {  // startup variables
   auto r = l;
   for( ++r ; ( r != tmod->subset().end() ) && ( *r < gr ) ; )
    ++r;
   Subset nms( std::distance( l , r ) );
   auto nvit = nv.begin();
   auto nmsit = nms.begin();
   while( l != r ) {
    *(nvit++) = qf->get_linear_coefficient( *l );
    *(nmsit++) = *(l++);
   }
   set_startup_costs( nv.begin() , std::move( nms ) , true , par , eDryRun );
   if( r == tmod->subset().end() )
    return;
  }

  gl = gr;
  gr = 2 * th - init_t;

  if( *l < gr ) {  // active power variables
   auto r = l;
   for( ++r ; ( r != tmod->subset().end() ) && ( *r < gr ) ; )
    ++r;
   Subset nms( std::distance( l , r ) );
   auto nvit = nv.begin();
   auto nmsit = nms.begin();
   while( l != r ) {
    *(nvit++) = qf->get_linear_coefficient( *l );
    *(nmsit++) = *(l++) - gl;
   }
   set_linear_term( nv.begin() , std::move( nms ) , true , par , eDryRun );
   if( r == tmod->subset().end() )
    return;
  }

  gl = gr;
  gr = 3 * th - init_t;

  if( *l < gr ) {  // commitment variables
   auto r = l;
   for( ++r ; ( r != tmod->subset().end() ) && ( *r < gr ) ; )
    ++r;
   Subset nms( std::distance( l , r ) );
   auto nvit = nv.begin();
   auto nmsit = nms.begin();
   while( l != r ) {
    *(nvit++) = qf->get_linear_coefficient( *l );
    *(nmsit++) = *(l++) - gl;
   }
   set_const_term( nv.begin() , std::move( nms ) , true , par , eDryRun );
   if( r == tmod->subset().end() )
    return;
  }

  gl = gr;
  gr = 4 * th - init_t;

  if( *l < gr ) {  // primary spinning reserve variables
   auto r = l;
   for( ++r ; ( r != tmod->subset().end() ) && ( *r < gr ) ; )
    ++r;
   Subset nms( std::distance( l , r ) );
   auto nvit = nv.begin();
   auto nmsit = nms.begin();
   while( l != r ) {
    *(nvit++) = qf->get_linear_coefficient( *l );
    *(nmsit++) = *(l++) - gl;
   }
   set_primary_spinning_reserve_cost( nv.begin() , std::move( nms ) ,
                                      true , par , eDryRun );
   if( r == tmod->subset().end() )
    return;
  }

  gl = gr;
  gr = 5 * th - init_t;

  if( *l < gr ) {  // secondary spinning reserve variables
   auto r = l;
   for( ++r ; ( r != tmod->subset().end() ) && ( *r < gr ) ; )
    ++r;
   Subset nms( std::distance( l , r ) );
   auto nvit = nv.begin();
   auto nmsit = nms.begin();
   while( l != r ) {
    *(nvit++) = qf->get_linear_coefficient( *l );
    *(nmsit++) = *(l++) - gl;
   }
   set_secondary_spinning_reserve_cost( nv.begin() , std::move( nms ) ,
                                        true , par , eDryRun );
  }

  if( l != tmod->subset().end() )
   throw( std::invalid_argument(
    "ThermalUnitBlock::handle_objective_change: invalid variable in "
    "C05FunctionModLinSbst." ) );
  return;

 }  // end( C05FunctionModLinSbst )

 throw( std::invalid_argument(
  "ThermalUnitBlock:: unsupported FunctionMod from Objective." ) );

}  // end( ThermalUnitBlock::handle_objective_change )

/*--------------------------------------------------------------------------*/
/*------------------- End File ThermalUnitBlock.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
