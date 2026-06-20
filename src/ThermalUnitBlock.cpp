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
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

template< typename T >
static bool identical( std::vector< T > & vec , const Block::Subset sbst ,
                       typename std::vector< T >::const_iterator it ) {
 // returns true if the sub-vector of vec[] corresponding to the indices
 // in sbst is identical to the vector starting at it
 for( auto t : sbst )
  if( vec[ t ] != *( it++ ) )
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
  vec[ t ] = *( it++ );
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

// register ThermalUnitBlockSolution to the Solution factory
SMSpp_insert_in_factory_cpp_0( ThermalUnitBlockSolution );

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
 // Deserialize data from the base class
 UnitBlock::deserialize( group );

 // Mandatory variables

 ::deserialize( group , "MaxPower" , f_time_horizon , v_MaxPower ,
                false , true , v_change_intervals );

 // Optional variables

 ::deserialize( group , f_InvestmentCost , "InvestmentCost" );

 ::deserialize( group , f_Capacity , "Capacity" );

 ::deserialize( group , f_scale , "Scale" );

 if( ::deserialize( group , f_MinUpTime , "MinUpTime" ) )
  f_MinUpTime = std::min( std::max( f_MinUpTime , static_cast< Index >( 1 ) ) ,
                          f_time_horizon );

 if( ::deserialize( group , f_MinDownTime , "MinDownTime" ) )
  f_MinDownTime = std::min( std::max( f_MinDownTime , static_cast< Index >( 1 ) ) ,
                            f_time_horizon );

 ::deserialize( group , f_InitialPower , "InitialPower" );

 if( ! ::deserialize( group , f_InitUpDownTime , "InitUpDownTime" ) ) {
  if( f_InitialPower == 0 )
   f_InitUpDownTime = -f_MinDownTime;
  else
   f_InitUpDownTime = f_MinUpTime;
 }

 if( ! ::deserialize( group , "MinPower" , f_time_horizon , v_MinPower ,
                      true , true , v_change_intervals ) )
  v_MinPower.resize( f_time_horizon );

 if( ! ::deserialize( group , "Availability" , f_time_horizon ,
		      v_Availability , true , true , v_change_intervals ) )
  v_Availability.resize( f_time_horizon , 1 );

 if( ! ::deserialize( group , "LinearTerm" , f_time_horizon , v_LinearTerm ,
                      true , true , v_change_intervals ) )
  v_LinearTerm.resize( f_time_horizon );

 // the spinning-reserve linear costs (the dual of the reserve demand pushed
 // onto the unit, the counterpart of LinearTerm); optional -- left empty when
 // absent, i.e. when the unit prices no reserve
 ::deserialize( group , "PrimarySpinningReserveCost" , f_time_horizon ,
                v_PrimarySpinningReserveCost , true , true , v_change_intervals );
 ::deserialize( group , "SecondarySpinningReserveCost" , f_time_horizon ,
                v_SecondarySpinningReserveCost , true , true ,
                v_change_intervals );

 if( ! ::deserialize( group , "QuadTerm" , f_time_horizon , v_QuadTerm ,
                      true , true , v_change_intervals ) )
  v_QuadTerm.resize( f_time_horizon );

 if( ! ::deserialize( group , "ConstTerm" , f_time_horizon , v_ConstTerm ,
                      true , true , v_change_intervals ) )
  v_ConstTerm.resize( f_time_horizon );

 if( ! ::deserialize( group , "StartUpCost" , f_time_horizon , v_StartUpCost ,
                      true , true , v_change_intervals ) )
  v_StartUpCost.resize( f_time_horizon );

 ::deserialize( group , "DeltaRampUp" , f_time_horizon , v_DeltaRampUp ,
                true , true , v_change_intervals );

 ::deserialize( group , "DeltaRampDown" , f_time_horizon , v_DeltaRampDown ,
                true , true , v_change_intervals );

 ::deserialize( group , "FixedConsumption" , f_time_horizon ,
                v_FixedConsumption , true , true , v_change_intervals );

 ::deserialize( group , "InertiaCommitment" , f_time_horizon ,
                v_InertiaCommitment , true , true , v_change_intervals );

 if( ! ( f_ignore_netcdf_vars & 1 ) ) {
  ::deserialize( group , "PrimaryRho" , f_time_horizon , v_PrimaryRho ,
		 true , true , v_change_intervals );
  if( std::all_of( v_PrimaryRho.begin() , v_PrimaryRho.end() ,
		   []( double i ) { return( i == 0 ); } ) )
   v_PrimaryRho.clear();

  ::deserialize( group , "SecondaryRho" , f_time_horizon , v_SecondaryRho ,
                 true , true , v_change_intervals );
  if( std::all_of( v_SecondaryRho.begin() , v_SecondaryRho.end() ,
		   []( double i ) { return( i == 0 ); } ) )
   v_SecondaryRho.clear();
  }

 if( ::deserialize( group , f_fixToMax , "FixToMaximum" ) )
  f_fixToMax = std::max( f_fixToMax , 0 );

 // variables for AC elements
 if( ::deserialize( group , "MaxReactivePower" , f_time_horizon ,
		    v_MaxReactivePower , true , true , v_change_intervals ) )
  if( std::all_of( v_MaxReactivePower.begin() , v_MaxReactivePower.end() ,
		   []( double i ) { return( i == 0 ); } ) )
   v_MaxReactivePower.clear();

 if( ::deserialize( group , "MinReactivePower" , f_time_horizon ,
		    v_MinReactivePower , true , true , v_change_intervals ) )
  if( std::all_of( v_MinReactivePower.begin() , v_MinReactivePower.end() ,
		   []( double i ) { return( i == 0 ); } ) )
   v_MinReactivePower.clear();

 // variables for the reference schedule
 ::deserialize( group , "ReferenceSchedule" , f_time_horizon ,
		v_RefSchedule , true , true , v_change_intervals );

 // start-up, shut-down limits
 if( ! ::deserialize( group , "StartUpLimit" , f_time_horizon ,
		      v_StartUpLimit , true , true , v_change_intervals ) ) {
  v_StartUpLimit.resize( f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   v_StartUpLimit[ t ] = get_operational_min_power( t );
  }

 if( ! ::deserialize( group , "ShutDownLimit" , f_time_horizon ,
		      v_ShutDownLimit , true , true , v_change_intervals ) ) {
  v_ShutDownLimit.resize( f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   v_ShutDownLimit[ t ] = get_operational_min_power( t );
  }

 if( ! ::deserialize( group , "MaxRampUpSteps" , f_time_horizon ,
                      v_MaxRampSteps , true , true , v_change_intervals ) ) {
  v_MaxRampSteps.resize( f_time_horizon + 1 );
  if( f_InitUpDownTime > 0 ) {
   const auto delta_ramp_up = get_delta_ramp_up( 0 );
   if( delta_ramp_up == 0 )
    v_MaxRampSteps[ 0 ] = static_cast< int >( f_time_horizon - 1 );
   else
    v_MaxRampSteps[ 0 ] = std::min( static_cast< int >(
     ( ( get_operational_max_power( 0 ) - f_InitialPower ) /
       delta_ramp_up ) ) , static_cast< int >( f_time_horizon - 1 ) );
  }
  else
   v_MaxRampSteps[ 0 ] = -1;
  for( Index t = 1 ; t <= f_time_horizon ; ++t ) {
   const auto delta_ramp_up = get_delta_ramp_up( t - 1 );
   if( delta_ramp_up == 0 )
    v_MaxRampSteps[ t ] = static_cast< int >( f_time_horizon - t );
   else
    v_MaxRampSteps[ t ] = std::min( static_cast< int >(
     ( ( get_operational_max_power( t - 1 ) - get_operational_min_power( t - 1 ) ) /
       delta_ramp_up ) ) , static_cast< int >( f_time_horizon - t ) );
  }
 }

 if( ! ::deserialize( group , "MaxRampDownSteps" , f_time_horizon ,
                      v_MaxRampDownSteps , true , true , v_change_intervals ) ) {
  v_MaxRampDownSteps.resize( f_time_horizon + 1 );
  if( f_InitUpDownTime > 0 ) {
   const auto delta_ramp_down = get_delta_ramp_down( 0 );
   if( delta_ramp_down == 0 )
    v_MaxRampDownSteps[ 0 ] = static_cast< int >( f_time_horizon - 1 );
   else
    v_MaxRampDownSteps[ 0 ] = std::min( static_cast< int >(
     ( ( f_InitialPower - get_operational_min_power( 0 ) ) /
       delta_ramp_down ) ) , static_cast< int >( f_time_horizon - 1 ) );
  }
  else
   v_MaxRampDownSteps[ 0 ] = -1;
  for( Index t = 1 ; t <= f_time_horizon ; ++t ) {
   const auto delta_ramp_down = get_delta_ramp_down( t - 1 );
   if( delta_ramp_down == 0 )
    v_MaxRampDownSteps[ t ] = f_time_horizon - t;
   else
    v_MaxRampDownSteps[ t ] = std::min( static_cast< int >(
     ( ( get_operational_max_power( t - 1 ) - get_operational_min_power( t - 1 ) ) /
       delta_ramp_down ) ) , static_cast< int >( f_time_horizon - t ) );
  }
 }

 // enforce the InitialPower >= MinPower contract for an initially-on unit:
 // if InitUpDownTime > 0 the unit was on at instant -1, hence InitialPower
 // must be >= MinPower (see the comments on "InitialPower" in the .h). Some
 // data sources (e.g. plan4res) violate this, encoding an initially-on unit
 // with InitialPower < MinPower; left as-is the ramp constraints of the
 // various formulations would disagree on the implied initial trajectory.
 // Clamp InitialPower up to MinPower (the least contract-satisfying value)
 // and warn, so that all formulations and the DP solvers stay consistent
 if( ( f_InitUpDownTime > 0 ) &&
     ( f_InitialPower < get_operational_min_power( 0 ) ) ) {
  std::cerr << "ThermalUnitBlock::deserialize: warning: initially-on unit "
               "(InitUpDownTime = " << f_InitUpDownTime << ") has InitialPower "
            << f_InitialPower << " < MinPower "
            << get_operational_min_power( 0 )
            << "; clamping InitialPower to MinPower to satisfy the contract"
            << std::endl;
  f_InitialPower = get_operational_min_power( 0 );
  }

 check_data_consistency();

}  // end( ThermalUnitBlock::deserialize )

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG

/*
std::vector< std::string > ThermalUnitBlock::expected_dims( void )
 const {
 static const std::vector< std::string > ed = { };

 auto ret = UnitBlock::expected_dims();
 ret.insert( ret.end() , ed.begin() , ed.end() );

 return( ret );
 }

----------------------------------------------------------------------------*/

std::vector< std::string > ThermalUnitBlock::expected_vars( void )
 const {
 static const std::vector< std::string > ev =
  { "InvestmentCost" , "Capacity" , "Scale" , "MinPower" , "MaxPower" ,
    "DeltaRampUp" , "DeltaRampDown" , "PrimaryRho" , "SecondaryRho" ,
    "PrimarySpinningReserveCost" , "SecondarySpinningReserveCost" ,
    "LinearTerm" , "QuadTerm" , "ConstTerm" , "StartUpCost" ,
    "FixedConsumption" , "InertiaCommitment" , "InitialPower" , "MinUpTime" ,
    "MinDownTime" , "InitUpDownTime" , "Availability" , "StartUpLimit" ,
    "ShutDownLimit" , "MaxRampUpSteps" , "MaxRampDownSteps" ,
    "InitialReactivePower", "MaxReactivePower" , "MinReactivePower" ,
    "ReferenceSchedule" , "FixToMaximum"
    };

 auto ret = UnitBlock::expected_vars();
 ret.insert( ret.end() , ev.begin() , ev.end() );

 return( ret );
 }

#endif

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
                            "is " + std::to_string( v_MaxPower[ t ] ) ) );

  if( v_MinPower[ t ] < 0 )
   throw( std::logic_error( "ThermalUnitBlock::check_data_consistency: "
                            "minimum power at time " + std::to_string( t ) +
			                         " is " + std::to_string( v_MinPower[ t ] ) +
                            ", but it must be nonnegative" ) );
  }

 // Availability- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_Availability.empty() ) {
  assert( v_Availability.size() == f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   if( ( v_Availability[ t ] < 0 ) || ( v_Availability[ t ] > 1 ) )
    throw( std::logic_error( "ThermalUnitBlock::check_data_consistency: "
                             "availability at time " + std::to_string( t ) +
			                          " is " + std::to_string( v_Availability[ t ] ) +
                             ", but it must be between 0 and 1" ) );
  }

 // Delta ramp-up - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_DeltaRampUp.empty() ) {
  assert( v_DeltaRampUp.size() == f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   if( v_DeltaRampUp[ t ] < 0 )
    throw( std::logic_error( "ThermalUnitBlock::check_data_consistency: "
                             "delta ramp-up at time " +
                             std::to_string( t ) + " is " +
                             std::to_string( v_DeltaRampUp[ t ] ) +
                             ", but it must be nonnegative" ) );
  }

 // Delta ramp-down - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_DeltaRampDown.empty() ) {
  assert( v_DeltaRampDown.size() == f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   if( v_DeltaRampDown[ t ] < 0 )
    throw( std::logic_error( "ThermalUnitBlock::check_data_consistency: "
                             "delta ramp-down at time " +
                             std::to_string( t ) + " is " +
                             std::to_string( v_DeltaRampDown[ t ] ) +
                             ", but it must be nonnegative" ) );
  }

 // QuadTerm - - - - - - - - - - - - - - - -- - - - - - - - - - - - - - - - -
 if( ! v_QuadTerm.empty() ) {
  assert( v_QuadTerm.size() == f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   if( v_QuadTerm[ t ] < 0 )
    throw( std::logic_error( "ThermalUnitBlock::check_data_consistency: "
                             "quadratic term at time " +
                             std::to_string( t ) + " is " +
                             std::to_string( v_QuadTerm[ t ] ) +
                             ", but it must be nonnegative" ) );
  }

 // InitialPower- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( f_InitialPower < 0 )
  throw( std::logic_error( "ThermalUnitBlock::check_data_consistency: "
                           "initial power is " +
                           std::to_string( f_InitialPower ) +
                           ", but it must be nonnegative" ) );

 // if the unit is initially committed, the initial power must lie within
 // the operational band of the unit: formulations and solvers (e.g. the
 // DP ones) rely on this documented contract, and violating it silently
 // produces spurious must-run behaviour rather than a clean error
 if( ( f_InitUpDownTime > 0 ) && ( f_InitialPower < v_MinPower.front() ) )
  throw( std::logic_error( "ThermalUnitBlock::check_data_consistency: the "
                           "unit is initially committed but the initial "
                           "power is " + std::to_string( f_InitialPower ) +
                           ", which is below the minimum power, which is " +
                           std::to_string( v_MinPower.front() ) ) );

 // StartUpLimit- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 for( Index t = 0 ; t < f_time_horizon ; ++t )
  if( ( v_StartUpLimit[ t ] < get_operational_min_power( t ) ) ||
      ( v_StartUpLimit[ t ] > get_operational_max_power( t ) ) )
   throw( std::logic_error( "ThermalUnitBlock::check_data_consistency: "
                            "start-up limit at time " +
                            std::to_string( t ) + " is " +
                            std::to_string( v_StartUpLimit[ t ] ) +
                            ", but it must be between " +
                            std::to_string( get_operational_min_power( t ) )
                            + " and " +
                            std::to_string( get_operational_max_power( t ) ) ) );

 // ShutDownLimit - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 for( Index t = 0 ; t < f_time_horizon ; ++t )
  if( ( v_ShutDownLimit[ t ] < get_operational_min_power( t ) ) ||
      ( v_ShutDownLimit[ t ] > get_operational_max_power( t ) ) )
   throw( std::logic_error( "ThermalUnitBlock::check_data_consistency: "
                            "shut-down limit at time " +
                            std::to_string( t ) + " is " +
                            std::to_string( v_ShutDownLimit[ t ] ) +
                            ", but it must be between " +
                            std::to_string( get_operational_min_power( t ) )
                            + " and " +
                            std::to_string( get_operational_max_power( t ) )
                            ) );

 }  // end( ThermalUnitBlock::check_data_consistency )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::generate_abstract_variables( Configuration * stvv )
{
 if( variables_generated() )  // variables have already been generated
  return;                     // nothing to do

 UnitBlock::generate_abstract_variables( stvv );

 // read Configuration to set the formulation - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
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
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( f_InvestmentCost != 0 ) {
  design.set_type( ColVariable::kBinary );
  add_static_variable( design , "x_thermal" );
  }
 else
  design.set_value( std::numeric_limits< double >::quiet_NaN() );

 // Commitment Variables- - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 v_commitment.resize( f_time_horizon );
 for( auto & var : v_commitment )
  var.set_type( ColVariable::kBinary );
 add_static_variable( v_commitment , "u_thermal" );

 // Active Power Variables- - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 v_active_power.resize( f_time_horizon );
 for( auto & var : v_active_power )
  var.set_type( ColVariable::kNonNegative );
 add_static_variable( v_active_power , "p_thermal" );

 // Reactive Power Variables, if any- - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( f_reactive_power ) {
  v_reactive_power.resize( f_time_horizon );
  for( auto & var : v_reactive_power )
   var.set_type( ColVariable::kNonNegative );
  add_static_variable( v_reactive_power , "q_thermal" );
  }

 // Start-Up and Shut-Down Binary Variables - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 auto startup_shutdown_size = f_time_horizon - init_t;

 if( startup_shutdown_size > 0 ) {
  auto vartype = ( wf & ZWCont ) ? ColVariable::kPosUnitary
                                 : ColVariable::kBinary;

  v_start_up.resize( startup_shutdown_size );
  for( auto & var : v_start_up )
   var.set_type( vartype );
  add_static_variable( v_start_up , "v_thermal" );

  v_shut_down.resize( startup_shutdown_size );
  for( auto & var : v_shut_down )
   var.set_type( vartype );
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

 // possibly fix the commitment variables to 0 or 1 - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( f_InitUpDownTime > 0 ) {  // InitUpDownTime > 0- - - - - - - - - - - - -
  for( Index t = 0 ; t < init_t ; ++t ) {
   if( ! v_commitment.empty() ) {
    v_commitment[ t ].set_value( 1 );
    v_commitment[ t ].is_fixed( true , eNoMod );
    }
   }

  for( Index t = init_t ;
       t < std::min( init_t + f_MinDownTime , f_time_horizon ) ; ++t ) {
   v_start_up[ t - init_t ].set_value( 0 );
   v_start_up[ t - init_t ].is_fixed( true , eNoMod );
   }
  }
 else {  // InitUpDownTime <= 0 - - - - - - - - - - - - - - - - - - - - - - -
  for( Index t = 0 ; t < init_t ; ++t ) {
   if( ! v_active_power.empty() ) {
    v_active_power[ t ].set_value( 0 );
    v_active_power[ t ].is_fixed( true , eNoMod );
    }

   if( ! v_commitment.empty() ) {
    v_commitment[ t ].set_value( 0 );
    v_commitment[ t ].is_fixed( true , eNoMod );
    }

   if( ! v_primary_spinning_reserve.empty() ) {
    v_primary_spinning_reserve[ t ].set_value( 0 );
    v_primary_spinning_reserve[ t ].is_fixed( true , eNoMod );
    }

   if( ! v_secondary_spinning_reserve.empty() ) {
    v_secondary_spinning_reserve[ t ].set_value( 0 );
    v_secondary_spinning_reserve[ t ].is_fixed( true , eNoMod );
    }
   }

  for( Index t = init_t ;
       t < std::min( init_t + f_MinUpTime , f_time_horizon ) ; ++t ) {
   v_shut_down[ t - init_t ].set_value( 0 );
   v_shut_down[ t - init_t ].is_fixed( true , eNoMod );
   }
  }

 bool f_cuts = wf & PCuts;

 // Prospective Cuts Variables- - - - - - - - - - - - - - - - - - - - - - - -

 if( f_cuts ) {
  AR |= PCuts;

  prevpbar.resize( f_time_horizon , 0 );
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

  case( SUSDForm ):  // SUSD formulation- - - - - - - - - - - - - - - - - - -

   if( f_InitUpDownTime > 0 ) {  // if initial committed, OFF_0

    // intervals ( 0 , k ) represent the pre-horizon ON run terminating (i.e.,
    // the unit shutting down) at k; k == init_t == 0 is the degenerate ( 0 , 0
    // ) interval covering no in-horizon period, which enables an immediate
    // shut-down at t == 0 (legal only when init_t == 0, the minimum up time is
    // already met). Since ( 0 , 0 ) covers no period, the ramp / shut-down-
    // limit constraints cannot gate it: include it only when that shut-down is
    // actually feasible, i.e. InitialPower <= ShutDownLimit[ 0 ] (the unit can
    // reach 0 in one step). This lets these interval-based formulations match
    // 3bin/T, which leave the t == 0 shut-down variable free and gate it via
    // their ramp-down constraint
    for( Index k = init_t ; k <= f_time_horizon + 1 ; ++k ) {
     if( ( k == 0 ) && ( f_InitialPower > v_ShutDownLimit[ 0 ] ) )
      continue;
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

   if( ( wf & FormMsk ) == ptForm )  // pt formulation- - - - - - - - - - - -
    AR |= ptForm;

   if( ( wf & FormMsk ) == DPForm ) {  // DP formulation- - - - - - - - - - -
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
   }

   if( ( wf & FormMsk ) == SUForm || ( wf & FormMsk ) == SUSDForm )
   {  // SU or SUSD formulation - - - - - - - - - - - - - - - - - - - - - - -

     if( ( wf & FormMsk ) == SUForm )
       AR |= SUForm;
     if( ( wf & FormMsk ) == SUSDForm )
       AR |= SUSDForm;

     for( Index i = 0 ; i < v_Y_plus.size() ; ++i ) {
     bool check_var = true;
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

     // add variables for SUSD formulation maximizing perspective cuts
     // variables of SU and SD formulations - - - - - - - - - - - - - - - - -
     if( ( wf & FormMsk ) == SUSDForm ) {
      v_cut_teta.resize( f_time_horizon );
      for( auto & var : v_cut_teta )
       var.set_type( ColVariable::kNonNegative );
      add_static_variable( v_cut_teta , "teta_thermal" );
     }
    }
   }

   if( ( wf & FormMsk ) == SDForm || ( wf & FormMsk ) == SUSDForm )
   {  // SD or SUSD formulation - - - - - - - - - - - - - - - - - - - - - - -

    if( ( wf & FormMsk ) == SDForm )
     AR |= SDForm;
    if( ( wf & FormMsk ) == SUSDForm )
     AR |= SUSDForm;
    for( Index i = 0 ; i < v_Y_plus.size() ; ++i ) {
     bool check_var = true;
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

 // the variables for the reference schedule, if there- - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! v_RefSchedule.empty() ) {
  v_abs_ref_schedule.resize( f_time_horizon );
  for( auto & var : v_abs_ref_schedule )
   var.set_type( ColVariable::kNonNegative );
  add_static_variable( v_abs_ref_schedule , "v_abs_refschd" );
  }

 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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

 // commitment design binary variable constraints - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( f_InvestmentCost != 0 ) {
  CommitmentDesign_Const.resize( f_time_horizon );

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   vars.push_back( std::make_pair( & v_commitment[ t ] , 1 ) );
   vars.push_back( std::make_pair( & design , -1 ) );

   CommitmentDesign_Const[ t ].set_lhs( -Inf< double >() );
   CommitmentDesign_Const[ t ].set_rhs( 0 );
   CommitmentDesign_Const[ t ].set_function(
                                  new LinearFunction( std::move( vars ) ) );
   }

  add_static_constraint( CommitmentDesign_Const ,
                         "CommitmentDesign_Const_Thermal" );
  }

 // compute the \psi constants for DP-related formulations- - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 std::vector< double > v_psi;

 if( ( ( AR & FormMsk ) == ptForm ) || ( ( AR & FormMsk ) == SUForm ) ||
     ( ( AR & FormMsk ) == SDForm ) || ( ( AR & FormMsk ) == SUSDForm ) ) {
  v_psi.resize( v_P_h_k.size() );

  for( Index j = 0 ; j < v_P_h_k.size() ; ++j ) {
   Index t = v_P_h_k[ j ].first;
   v_psi[ j ] = get_operational_max_power( t );
   if( v_P_h_k[ j ].second.second <= f_time_horizon )
    if( ! v_DeltaRampDown.empty() )
     if( v_P_h_k[ j ].second.first > 0 )
      v_psi[ j ] = std::min( v_psi[ j ] ,
			     v_ShutDownLimit[ t ] + v_DeltaRampDown[ t ] *
			     ( v_P_h_k[ j ].second.second - ( t + 1 ) ) );

   if( ! v_DeltaRampUp.empty() ) {
    if( ( v_P_h_k[ j ].second.first == 0 ) && ( f_InitUpDownTime > 0 ) )
     v_psi[ j ] = std::min( v_psi[ j ] ,
			    f_InitialPower + v_DeltaRampUp[ t ] * ( t + 1 ) );

    if( v_P_h_k[ j ].second.first > 0 )
     v_psi[ j ] = std::min( v_psi[ j ] ,
			    v_StartUpLimit[ t ] + v_DeltaRampUp[ t ] *
			    ( ( t + 1 ) - v_P_h_k[ j ].second.first ) );
    }

   v_psi[ j ] = std::max( v_psi[ j ] , get_operational_min_power( t ) );
   }
  }

 // add the "fixed to maximum generation" constraints - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( f_fixToMax > 0 ) {
  fixed_to_max_Power_Const.resize( f_time_horizon );

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   // P_t >= Pmax(t)
   auto lfunck = new LinearFunction();
   lfunck->add_variable( & v_active_power[ t ], 1.0 );
   fixed_to_max_Power_Const[ t ].set_lhs( get_operational_max_power( t ) );
   fixed_to_max_Power_Const[ t ].set_rhs( Inf< double >() );
   fixed_to_max_Power_Const[ t ].set_function( lfunck );
   }
  add_static_constraint( fixed_to_max_Power_Const, "FixedGeneration" );
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
        t < f_time_horizon ; ++t, ++cnstr_idx ) {

     for( Index s = t - ( init_t + f_MinUpTime - 1 ) ; s <= t - init_t ; ++s )
      vars.push_back( std::make_pair( &v_start_up[ s ] , -1.0 ) );

     vars.push_back( std::make_pair( &v_commitment[ t ] , 1.0 ) );

     StartUp_Const[ cnstr_idx ].set_lhs( 0.0 );
     StartUp_Const[ cnstr_idx ].set_rhs( Inf< double >() );
     StartUp_Const[ cnstr_idx ].set_function(
      new LinearFunction( std::move( vars ) ) );
    }

    add_static_constraint( StartUp_Const ,
                           "StartUp_Commitment_Const_Thermal" );
   }

   // Initializing turn off constraints (shut-down constraints) - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   auto shutdown_const_size =
    static_cast< int >( f_time_horizon - ( init_t + f_MinDownTime - 1 ) );

   if( shutdown_const_size > 0 ) {

    ShutDown_Const.resize( shutdown_const_size );

    for( Index t = ( init_t + f_MinDownTime - 1 ) , cnstr_idx = 0 ;
	 t < f_time_horizon ; ++t, ++cnstr_idx ) {
     for( Index s = t - ( init_t + f_MinDownTime - 1 ) ; s <= t - init_t ;
	  ++s )
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
  case( SDForm ):  // SD formulation- - - - - - - - - - - - - - - - - - - - -
   // fall through
  case ( SUSDForm ):  {  // SUSD formulation- - - - - - - - - - - - - - - - -

   // Constraints connecting power variables of 3bin with those of DP, SU and
   // SD formulations - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   if( ( AR & FormMsk ) == ptForm ) {  // pt formulation- - - - - - - - - - -
    ;  // does nothing
   } else {

    if( ( AR & FormMsk ) == SUSDForm )
     Eq_ActivePower_Const.resize( 2 * f_time_horizon );
    else
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
    }

    if( ( AR & FormMsk ) == SUForm || ( AR & FormMsk ) == SUSDForm )
    {  // SU or SUSD formulation- - - - - - - - - - - - - - - - - - - - - - -

     for( Index t = 0 ; t < f_time_horizon ; ++t ) {

      vars.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );

      for( Index j = 0 ; j < v_P_h.size() ; ++j )
       if( v_P_h[ j ].first == t )
        vars.push_back( std::make_pair( &v_active_power_h[ j ] , -1.0 ) );

      Eq_ActivePower_Const[ t ].set_both( 0.0 );
      Eq_ActivePower_Const[ t ].set_function(
       new LinearFunction( std::move( vars ) ) );
     }
    }

    if( ( AR & FormMsk ) == SDForm || ( AR & FormMsk ) == SUSDForm )
    {  // SD or SUSD formulation- - - - - - - - - - - - - - - - - - - - - - -

     for( Index t = 0 ; t < f_time_horizon ; ++t ) {

      vars.push_back( std::make_pair( &v_active_power[ t ] , 1.0 ) );

      for( Index j = 0 ; j < v_P_k.size() ; ++j )
       if( v_P_k[ j ].first == t )
        vars.push_back( std::make_pair( &v_active_power_k[ j ] , -1.0 ) );

      if( ( AR & FormMsk ) == SUSDForm ) {
       Eq_ActivePower_Const[ f_time_horizon + t ].set_both( 0.0 );
       Eq_ActivePower_Const[ f_time_horizon + t ].set_function(
        new LinearFunction( std::move( vars ) ) );
      } else {
       Eq_ActivePower_Const[ t ].set_both( 0.0 );
       Eq_ActivePower_Const[ t ].set_function(
        new LinearFunction( std::move( vars ) ) );
      }
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

   add_static_constraint( Eq_Commitment_Const ,
                          "Eq_Commitment_Const_Thermal" );

   // Constraints connecting start-up variables of 3bin with those of pt,
   // DP, SU and SD formulations- - - - - - - - - - - - - - - - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Eq_StartUp_Const.resize( f_time_horizon - init_t );

   for( Index t = init_t , cnstr_idx = 0 ; t < f_time_horizon ;
       ++t, ++cnstr_idx ) {

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
       ++t, ++cnstr_idx ) {

    vars.push_back( std::make_pair( &v_shut_down[ t - init_t ] , 1.0 ) );

    for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
     if( ( v_Y_plus[ i ].first <= t ) && ( t == v_Y_plus[ i ].second ) )
      vars.push_back( std::make_pair( &v_commitment_plus[ i ] , -1.0 ) );

    Eq_ShutDown_Const[ cnstr_idx ].set_both( 0.0 );
    Eq_ShutDown_Const[ cnstr_idx ].set_function(
     new LinearFunction( std::move( vars ) ) );
   }

   add_static_constraint( Eq_ShutDown_Const ,
                          "Eq_ShutDown_Const_Thermal" );

   // Network Constraints - - - - - - - - - - - - - - - - - - - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   Network_Const.resize( v_nodes_plus.size() + v_nodes_minus.size() + 2 );

   auto cnstr_idx = 0;

  for( Index t = 0 ; t < 1 ; ++t, ++cnstr_idx ) {

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

  for( Index t = 0 ; t < v_nodes_plus.size() ; ++t, ++cnstr_idx ) {

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
 // new_SU = 1 if the SU formulation includes the new ramp-up constraints of the
 //            SUSD formulation and the corresponding new ramp-down constraints
 int new_SU = 0;
 // new_SU = 1 if the SU formulation includes the new ramp-down constraints of
 //            the SUSD formulation and the corresponding new ramp-up constraints
 int new_SD = 0;

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

   auto ramp_up_cnstrs_size = f_time_horizon;
   RampUp_Const.resize( ramp_up_cnstrs_size );
   auto cnstr_idx = 0;

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

    RampUp_Const[ cnstr_idx ].set_lhs( -Inf< double >() );
    if( ( t == 0 ) && ( f_InitUpDownTime > 0 ) )
     RampUp_Const[ cnstr_idx ].set_rhs( f_InitialPower + v_DeltaRampUp[ t ] );
    else if( ( t > 0 ) || ( ( t == 0 ) && ( f_InitUpDownTime <= 0 ) ) )
     RampUp_Const[ cnstr_idx ].set_rhs( 0.0 );
    RampUp_Const[ cnstr_idx ].set_function( new LinearFunction( std::move( vars ) ) );
    cnstr_idx++;
   }

  }
  if( ( AR & FormMsk ) == TForm ) {  // T formulation - - - - - - - - - - - -

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

  }
  if( ( AR & FormMsk ) == ptForm ) {  // pt formulation - - - - - - - - - - -

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
                                        -f_InitialPower - v_DeltaRampUp[ t ] ) );

     RampUp_Const[ cnstr_idx ].set_lhs( -Inf< double >() );
     RampUp_Const[ cnstr_idx ].set_rhs( 0.0 );
     RampUp_Const[ cnstr_idx ].set_function(
      new LinearFunction( std::move( vars ) ) );

     cnstr_idx++;
    }

  }
  if( ( AR & FormMsk ) == DPForm ) {  // DP formulation - - - - - - - - - - -

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
                                         -v_DeltaRampUp[ t ] - f_InitialPower ) );
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

  }
  if( ( AR & FormMsk ) == SUForm && new_SU == 0 ) {  // SU formulation- - - -

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
                                          -v_DeltaRampUp[ t ] - f_InitialPower ) );

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
          vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                          get_operational_min_power( t - 1 ) ) );
        }

       RampUp_Const[ cnstr_idx ].set_lhs( -Inf< double >() );
       RampUp_Const[ cnstr_idx ].set_rhs( 0.0 );
       RampUp_Const[ cnstr_idx ].set_function(
        new LinearFunction( std::move( vars ) ) );

       cnstr_idx++;
      }
     }
  }

  if( ( ( ( AR & FormMsk ) == SUForm ) && ( new_SU == 1 ) ) ||
          ( AR & FormMsk ) == SUSDForm ) {  // SUSD formulation - - - - - - -

   auto ramp_up_cnstrs_size = 0;

   if( f_InitUpDownTime > 0 )
    for( Index k = 0 ; k < v_MaxRampSteps[ 0 ] ; ++k )
     for( Index j = 0 ; j < v_P_h.size() ; ++j )
      if( ( v_P_h[ j ].first == k ) && ( v_P_h[ j ].second  == 0 ) )
       ramp_up_cnstrs_size++;

   for( Index t = 0 ; t < f_time_horizon ; ++t )
    for( Index j = 0 ; j < v_P_h.size() ; ++j )
     if( v_P_h[ j ].first == t )
      for( Index s = 0 ; s < v_P_h.size() ; ++s )
       if( v_P_h[ j ].second == v_P_h[ s ].second )
        for( Index k = 1 ; k <= v_MaxRampSteps[ t + 1 ] ; ++k )
         if( v_P_h[ s ].first == t + k )
          ramp_up_cnstrs_size++;

   RampUp_Const.resize( ramp_up_cnstrs_size );

   auto cnstr_idx = 0;

   if( f_InitUpDownTime > 0 )
    for( Index k = 0 ; k < v_MaxRampSteps[ 0 ] ; ++k )
     for( Index j = 0 ; j < v_P_h.size() ; ++j )
      if( ( v_P_h[ j ].first == k ) && ( v_P_h[ j ].second == 0 ) ) {

       vars.push_back( std::make_pair( &v_active_power_h[ j ] , 1.0 ) );

       for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
        if( v_P_h[ j ].second == v_Y_plus[ i ].first )
         if( k + 1 <= v_Y_plus[ i ].second )
          vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                          -( static_cast< int >( k + 1 ) *
                                             v_DeltaRampUp[ k ] )
                                          - f_InitialPower ) );

       RampUp_Const[ cnstr_idx ].set_lhs( -Inf< double >() );
       RampUp_Const[ cnstr_idx ].set_rhs( 0.0 );
       RampUp_Const[ cnstr_idx ].set_function(
        new LinearFunction( std::move( vars ) ) );

       cnstr_idx++;
      }

   for( Index t = 0 ; t < f_time_horizon ; ++t )
    for( Index j = 0 ; j < v_P_h.size() ; ++j )
     if( v_P_h[ j ].first == t )
      for( Index s = 0 ; s < v_P_h.size() ; ++s )
       if( v_P_h[ j ].second == v_P_h[ s ].second )
        for( Index k = 1 ; k <= v_MaxRampSteps[ t + 1 ] ; ++k )
         if( v_P_h[ s ].first == t + k ) {

          vars.push_back( std::make_pair( &v_active_power_h[ s ] , 1.0 ) );
          vars.push_back( std::make_pair( &v_active_power_h[ j ] , -1.0 ) );

          for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
           if( v_P_h[ j ].second == v_Y_plus[ i ].first ) {
            if( t + k + 1 <= v_Y_plus[ i ].second )
             vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                             -( static_cast< int >( k ) *
                                                v_DeltaRampUp[ t ] ) ) );
            if( ( t + 1 <= v_Y_plus[ i ].second ) &&
                ( t + k + 1 > v_Y_plus[ i ].second ) )
             vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                             get_operational_min_power( t ) ) );
           }

          RampUp_Const[ cnstr_idx ].set_lhs( -Inf< double >() );
          RampUp_Const[ cnstr_idx ].set_rhs( 0.0 );
          RampUp_Const[ cnstr_idx ].set_function(
           new LinearFunction( std::move( vars ) ) );

          cnstr_idx++;
         }
  }
  if( ( ( AR & FormMsk ) == SDForm ) && ( new_SD == 0 ) ) {  // SD formulation

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
                                          -v_DeltaRampUp[ t ] - f_InitialPower ) );

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

  if( ( ( AR & FormMsk ) == SDForm ) && ( new_SD == 1 ) ) {  // SD formulation
   auto ramp_up_cnstrs_size = 0;

   if( f_InitUpDownTime > 0 )
    for( Index k = 0 ; k < v_MaxRampSteps[ 0 ] ; ++k )
     for( Index j = 0 ; j < v_P_k.size() ; ++j )
      if( v_P_k[ j ].first == k )
       ramp_up_cnstrs_size++;

   for( Index t = 0 ; t < f_time_horizon ; ++t )
    for( Index j = 0 ; j < v_P_k.size() ; ++j )
     if( v_P_k[ j ].first == t )
      for( Index s = 0 ; s < v_P_k.size() ; ++s )
       if( v_P_k[ j ].second == v_P_k[ s ].second )
        for( Index k = 1 ; k <= v_MaxRampSteps[ t + 1 ] ; ++k )
         if( v_P_k[ s ].first == t + k )
          ramp_up_cnstrs_size++;

   RampUp_Const.resize( ramp_up_cnstrs_size );

   auto cnstr_idx = 0;

   if( f_InitUpDownTime > 0 )
    for( Index k = 0 ; k < v_MaxRampSteps[ 0 ] ; ++k )
     for( Index j = 0 ; j < v_P_k.size() ; ++j )
      if( v_P_k[ j ].first == k ) {

       vars.push_back( std::make_pair( &v_active_power_k[ j ] , 1.0 ) );

       for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
        if( v_P_k[ j ].second == v_Y_plus[ i ].second )
         if( 0 == v_Y_plus[ i ].first )
          vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                          -( static_cast< int >( k + 1 ) *
                                             v_DeltaRampUp[ k ] )
                                          - f_InitialPower ) );

       RampUp_Const[ cnstr_idx ].set_lhs( -Inf< double >() );
       RampUp_Const[ cnstr_idx ].set_rhs( 0.0 );
       RampUp_Const[ cnstr_idx ].set_function(
        new LinearFunction( std::move( vars ) ) );

       cnstr_idx++;
      }

   for( Index t = 0 ; t < f_time_horizon ; ++t )
    for( Index j = 0 ; j < v_P_k.size() ; ++j )
     if( v_P_k[ j ].first == t )
      for( Index s = 0 ; s < v_P_k.size() ; ++s )
       if( v_P_k[ j ].second == v_P_k[ s ].second )
        for( Index k = 1 ; k <= v_MaxRampSteps[ t + 1 ] ; ++k )
         if( v_P_k[ s ].first == t + k ) {

          vars.push_back( std::make_pair( &v_active_power_k[ s ] , 1.0 ) );
          vars.push_back( std::make_pair( &v_active_power_k[ j ] , -1.0 ) );

          for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
           if( v_P_k[ j ].second == v_Y_plus[ i ].second ) {
            if( t + 1 >= v_Y_plus[ i ].first )
             vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                             -( static_cast< int >( k ) *
                                                v_DeltaRampUp[ t ] ) ) );
            if( ( t + 1 < v_Y_plus[ i ].first ) &&
                ( t + k + 1 > v_Y_plus[ i ].first ) )
        for( Index ss = 0 ; ss < v_P_h_k.size() ; ++ss )
    if( ( t + k == v_P_h_k[ ss ].first ) &&
        ( v_P_h_k[ ss ].second.first == v_Y_plus[ i ].first ) &&
        ( v_P_h_k[ ss ].second.second == v_Y_plus[ i ].second ) )
      vars.push_back( std::make_pair( &v_commitment_plus[ i ] , -v_psi[ ss ] ) );

            if( ( t + k + 1 == v_Y_plus[ i ].first ) &&
                ( t + k + 1 < v_Y_plus[ i ].second ) )
            vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                             -v_StartUpLimit[ t + k + 1 ] ) );

       if( ( t + k + 1 == v_Y_plus[ i ].first ) &&
          ( t + k + 1 == v_Y_plus[ i ].second ) )
            vars.push_back(
             std::make_pair( &v_commitment_plus[ i ] ,
                                    -std::min( v_StartUpLimit[ t + k + 1 ],
                                               v_ShutDownLimit[ t + k + 1 ] ) ) );
           }

          RampUp_Const[ cnstr_idx ].set_lhs( -Inf< double >() );
          RampUp_Const[ cnstr_idx ].set_rhs( 0.0 );
          RampUp_Const[ cnstr_idx ].set_function(
           new LinearFunction( std::move( vars ) ) );

          cnstr_idx++;
         }

  }

  add_static_constraint( RampUp_Const , "RampUp_Const_Thermal" );
  }

 // initializing ramp-down constraints- - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ! v_DeltaRampDown.empty() ) {
  if( ( AR & FormMsk ) == tbinForm ) {  // 3bin formulation - - - - - - - - -

   RampDown_Const.resize( f_time_horizon );

   auto cnstr_idx = 0;

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {

    vars.push_back( std::make_pair( &v_active_power[ t ] , -1.0 ) );
    vars.push_back( std::make_pair( &v_commitment[ t ] ,
                                    -v_DeltaRampDown[ t ] ) );

    if( t > 0 )
     vars.push_back( std::make_pair( &v_active_power[ t - 1 ] , 1.0 ) );

    if( t >= init_t )
     vars.push_back( std::make_pair( &v_shut_down[ t - init_t ] ,
                                     -v_ShutDownLimit[ t ] ) );

    RampDown_Const[ cnstr_idx ].set_lhs( -Inf< double >() );
    if( ( t == 0 ) && ( f_InitUpDownTime > 0 ) )
     RampDown_Const[ cnstr_idx ].set_rhs( -f_InitialPower );
    else if( ( t > 0 ) || ( ( t == 0 ) && ( f_InitUpDownTime <= 0 ) ) )
     RampDown_Const[ cnstr_idx ].set_rhs( 0.0 );
    RampDown_Const[ cnstr_idx ].set_function( new LinearFunction( std::move( vars ) ) );
    cnstr_idx++;
   }

  }
  if( ( AR & FormMsk ) == TForm ) {  // T formulation - - - - - - - - - - - -

   RampDown_Const.resize( f_time_horizon );

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {

    vars.push_back( std::make_pair( &v_active_power[ t ] , -1.0 ) );

    if( t > 0 ) {
     vars.push_back( std::make_pair( &v_active_power[ t - 1 ] , 1.0 ) );
     vars.push_back( std::make_pair( &v_commitment[ t - 1 ] ,
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

  }
  if( ( AR & FormMsk ) == ptForm ) {  // pt formulation - - - - - - - - - - -

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
                                        f_InitialPower - v_DeltaRampDown[ t ] ) );

     RampDown_Const[ cnstr_idx ].set_lhs( -Inf< double >() );
     RampDown_Const[ cnstr_idx ].set_rhs( 0.0 );
     RampDown_Const[ cnstr_idx ].set_function(
      new LinearFunction( std::move( vars ) ) );

     cnstr_idx++;
    }

  }
  if( ( AR & FormMsk ) == DPForm ) {  // DP formulation - - - - - - - - - - -

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
                                         -v_DeltaRampDown[ t ] + f_InitialPower ) );
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

  }
  if( ( AR & FormMsk ) == SUForm && new_SU == 0 ) {  // SU formulation- - - -

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
                                          -v_DeltaRampDown[ t ] + f_InitialPower ) );

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

  }
 if( ( ( AR & FormMsk ) == SUForm ) && ( new_SU == 1 ) ) {  // SU formulation

    auto ramp_down_cnstrs_size = 0;

    if( f_InitUpDownTime > 0 )
     for( Index k = 0 ; k < v_MaxRampDownSteps[ 0 ] ; ++k )
      for( Index j = 0 ; j < v_P_h.size() ; ++j )
       if( v_P_h[ j ].first == k )
        ramp_down_cnstrs_size++;

    for( Index t = 0 ; t < f_time_horizon ; ++t )
     for( Index j = 0 ; j < v_P_h.size() ; ++j )
      if( v_P_h[ j ].first == t )
       for( Index s = 0 ; s < v_P_h.size() ; ++s )
        if( v_P_h[ j ].second == v_P_h[ s ].second )
         for( Index k = 1 ; k <= v_MaxRampDownSteps[ t + 1 ] ; ++k )
          if( v_P_h[ s ].first == t + k )
           ramp_down_cnstrs_size++;

    RampDown_Const.resize( ramp_down_cnstrs_size );

    auto cnstr_idx = 0;

    if( f_InitUpDownTime > 0 )
     for( Index k = 0 ; k < v_MaxRampDownSteps[ 0 ] ; ++k )
      for( Index j = 0 ; j < v_P_h.size() ; ++j )
       if( v_P_h[ j ].first == k ) {

       vars.push_back( std::make_pair( &v_active_power_h[ j ] , -1.0 ) );

       for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
        if( v_P_h[ j ].second == v_Y_plus[ i ].first )
         if( k + 1 <= v_Y_plus[ i ].second )
          vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                         -( static_cast< int >( k + 1 ) *
                                            v_DeltaRampDown[ k ] )
                                         + f_InitialPower ) );

       RampDown_Const[ cnstr_idx ].set_lhs( -Inf< double >() );
       RampDown_Const[ cnstr_idx ].set_rhs( 0.0 );
       RampDown_Const[ cnstr_idx ].set_function(
        new LinearFunction( std::move( vars ) ) );

       cnstr_idx++;
      }

   for( Index t = 0 ; t < f_time_horizon ; ++t )
    for( Index j = 0 ; j < v_P_h.size() ; ++j )
     if( v_P_h[ j ].first == t )
      for( Index s = 0 ; s < v_P_h.size() ; ++s )
       if( v_P_h[ j ].second == v_P_h[ s ].second )
        for( Index k = 1 ; k <= v_MaxRampDownSteps[ t + 1 ] ; ++k )
         if( v_P_h[ s ].first == t + k ) {

          vars.push_back( std::make_pair( &v_active_power_h[ s ] , -1.0 ) );
          vars.push_back( std::make_pair( &v_active_power_h[ j ] , 1.0 ) );

          for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
           if( v_P_h[ j ].second == v_Y_plus[ i ].first )
           {
            if( t + k + 1 <= v_Y_plus[ i ].second )
             vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                             -( static_cast< int >( k ) *
                                                v_DeltaRampDown[ t ] ) ) );
            if( ( t + 1 < v_Y_plus[ i ].second ) &&
                ( t + k + 1 > v_Y_plus[ i ].second ) )
             for( Index ss = 0 ; ss < v_P_h_k.size() ; ++ss )
              if( ( t == v_P_h_k[ ss ].first ) &&
               ( v_P_h_k[ ss ].second.first == v_Y_plus[ i ].first ) &&
               ( v_P_h_k[ ss ].second.second == v_Y_plus[ i ].second ) )
               vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                               -v_psi[ ss ] ) );

            if( ( t + 1 == v_Y_plus[ i ].second ) &&
                ( t + 1 > v_Y_plus[ i ].first ) )
             vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                             -v_ShutDownLimit[ t ] ) );

            if( ( t + 1 == v_Y_plus[ i ].second ) &&
                ( t + 1 == v_Y_plus[ i ].first ) )
             vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                             -std::min( v_StartUpLimit[ t ] ,
                                              v_ShutDownLimit[ t ] ) ) );
           }

          RampDown_Const[ cnstr_idx ].set_lhs( -Inf< double >() );
          RampDown_Const[ cnstr_idx ].set_rhs( 0.0 );
          RampDown_Const[ cnstr_idx ].set_function(
           new LinearFunction( std::move( vars ) ) );

          cnstr_idx++;
   }
 }
  if( ( ( ( AR & FormMsk ) == SDForm ) && ( new_SD == 1 ) ) ||
          ( AR & FormMsk ) == SUSDForm ) {  // SUSD formulation - - - - - - -
   auto ramp_down_cnstrs_size = 0;

   if( f_InitUpDownTime > 0 )
    for( Index k = 0 ; k < v_MaxRampDownSteps[ 0 ] ; ++k )
     for( Index j = 0 ; j < v_P_k.size() ; ++j )
      if( v_P_k[ j ].first == k )
       ramp_down_cnstrs_size++;

   for( Index t = 0 ; t < f_time_horizon ; ++t )
    for( Index j = 0 ; j < v_P_k.size() ; ++j )
     if( v_P_k[ j ].first == t )
      for( Index s = 0 ; s < v_P_k.size() ; ++s )
       if( v_P_k[ j ].second == v_P_k[ s ].second )
        for( Index k = 1 ; k <= v_MaxRampDownSteps[ t + 1 ] ; ++k )
         if( v_P_k[ s ].first == t + k )
          ramp_down_cnstrs_size++;

   RampDown_Const.resize( ramp_down_cnstrs_size );

   auto cnstr_idx = 0;

   if( f_InitUpDownTime > 0 )
    for( Index k = 0 ; k < v_MaxRampDownSteps[ 0 ] ; ++k )
     for( Index j = 0 ; j < v_P_k.size() ; ++j )
      if( v_P_k[ j ].first == k ) {

       vars.push_back( std::make_pair( &v_active_power_k[ j ] , -1.0 ) );

       for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
   if( v_P_k[ j ].second == v_Y_plus[ i ].second )
         if( 0 == v_Y_plus[ i ].first )
          vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                          -( static_cast< int >( k + 1 ) *
                                             v_DeltaRampDown[ k ] )
                                          + f_InitialPower ) );

       RampDown_Const[ cnstr_idx ].set_lhs( -Inf< double >() );
       RampDown_Const[ cnstr_idx ].set_rhs( 0.0 );
       RampDown_Const[ cnstr_idx ].set_function(
        new LinearFunction( std::move( vars ) ) );

       cnstr_idx++;

      }

   for( Index t = 0 ; t < f_time_horizon ; ++t )
    for( Index j = 0 ; j < v_P_k.size() ; ++j )
     if( v_P_k[ j ].first == t )
      for( Index s = 0 ; s < v_P_k.size() ; ++s )
       if( v_P_k[ j ].second == v_P_k[ s ].second )
        for( Index k = 1 ; k <= v_MaxRampDownSteps[ t + 1 ] ; ++k )
         if( v_P_k[ s ].first == t + k ) {

          vars.push_back( std::make_pair( &v_active_power_k[ s ] , -1.0 ) );
          vars.push_back( std::make_pair( &v_active_power_k[ j ] , 1.0 ) );

          for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
           if( v_P_k[ j ].second == v_Y_plus[ i ].second ) {
            if( t + 1 >= v_Y_plus[ i ].first )
             vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                             -( static_cast< int >( k ) *
                                                v_DeltaRampDown[ t ] ) ) );
            if( ( t + 1 < v_Y_plus[ i ].first ) &&
                ( t + k + 1 >= v_Y_plus[ i ].first ) )
             vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                              get_operational_min_power( t ) ) );
           }

          RampDown_Const[ cnstr_idx ].set_lhs( -Inf< double >() );
          RampDown_Const[ cnstr_idx ].set_rhs( 0.0 );
          RampDown_Const[ cnstr_idx ].set_function(
           new LinearFunction( std::move( vars ) ) );

          cnstr_idx++;
         }
  }
  if( ( ( AR & FormMsk ) == SDForm ) && ( new_SD == 0 ) ) {  // SD formulation
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
                                          -v_DeltaRampDown[ t ] + f_InitialPower ) );

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
          vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
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

 }
 if( ( AR & FormMsk ) == ptForm ) {  // pt formulation - - - - - - - -

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

 }
 if( ( AR & FormMsk ) == DPForm ) {  // DP formulation - - - - - - - -

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
 }

 if( ( AR & FormMsk ) == SUForm || ( AR & FormMsk ) == SUSDForm )
 {  // SU and SUSD formulation- - - - - - - - - - - - - - - - - - - - - - - -

  if( ( AR & FormMsk ) == SUSDForm )
   MinPower_Const.resize( v_P_h.size() + v_P_k.size() );
  else
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

 }
 if( ( AR & FormMsk ) == SDForm || ( AR & FormMsk ) == SUSDForm )
 {  // SD and SUSD formulation- - - - - - - - - - - - - - - - - - - - - - - -

  if( ( AR & FormMsk ) == SDForm )
   MinPower_Const.resize( v_P_k.size() );

  for( Index j = 0 ; j < v_P_k.size() ; ++j ) {

   auto t = v_P_k[ j ].first;
   for( Index i = 0 ; i < v_Y_plus.size() ; ++i )
    if( ( v_P_k[ j ].second == v_Y_plus[ i ].second ) &&
        ( v_Y_plus[ i ].first <= t + 1 ) )
     vars.push_back( std::make_pair( &v_commitment_plus[ i ] ,
                                     -get_operational_min_power( t ) ) );

   vars.push_back( std::make_pair( &v_active_power_k[ j ] , 1.0 ) );

   if( ( AR & FormMsk ) == SUSDForm ) {
    MinPower_Const[ v_P_h.size() + j ].set_lhs( 0.0 );
    MinPower_Const[ v_P_h.size() + j ].set_rhs( Inf< double >() );
    MinPower_Const[ v_P_h.size() + j ].set_function(
     new LinearFunction( std::move( vars ) ) );
   } else {
    MinPower_Const[ j ].set_lhs( 0.0 );
    MinPower_Const[ j ].set_rhs( Inf< double >() );
    MinPower_Const[ j ].set_function( new LinearFunction( std::move( vars ) ) );
   }
  }
 }

 add_static_constraint( MinPower_Const , "MinPower_Const_Thermal" );

 // initializing maximum power constraints- - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( ( AR & FormMsk ) == tbinForm ) {  // 3bin formulation- - - - - - - - - -

  MaxPower_Const.resize( f_MinUpTime == 1 ?
                         ( init_t == 0 ?
                           2 * ( f_time_horizon - init_t ) - 2 + init_t :
                           2 * ( f_time_horizon - init_t ) - 1 + init_t ) :
                         f_time_horizon );

  for( Index t = 0 , cnstr_idx = 0 ; t < f_time_horizon ; ++t , ++cnstr_idx ) {

   if( t >= init_t ) {
    if( t == 0 ) {
     vars.push_back( std::make_pair( &v_shut_down[ t + 1 - init_t ] ,
                                     v_ShutDownLimit[ t ] -
                                     get_operational_max_power( t ) ) );
     // a start-up at t == 0 (the unit was off before the horizon, init_t == 0)
     // is bounded by the start-up cap exactly like an interior start-up: add
     // the start-up term so the reserve band p + pr + sr <= StartUpLimit is
     // enforced here too. Without it the band would use the full max_power at a
     // t == 0 start-up (looser than the T/pt formulations and the DP solvers).
     // Energy-only is unaffected: the ramp-up constraint already caps p <=
     // StartUpLimit at a start-up, so the term only binds the reserves.
     if( ( f_MinUpTime != 1 ) && ( t < f_time_horizon - 1 ) )
      vars.push_back( std::make_pair( &v_start_up[ t - init_t ] ,
                                      v_StartUpLimit[ t ] -
                                      get_operational_max_power( t ) ) );
     }
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

 }
 if( ( AR & FormMsk ) == TForm ) {  // T formulation - - - - - - - - -

  Index max_power_cnstrs_size = 0;

  std::vector< int > v_T_RU;
  std::vector< int > v_K_SU;
  std::vector< int > v_K_SD;

  if( ( ! v_DeltaRampUp.empty() ) || ( ! v_DeltaRampDown.empty() ) ) {
   std::vector< int > v_T_RD;

   if( ( ! v_DeltaRampUp.empty() ) ) {
    v_T_RU.resize( f_time_horizon );
    v_K_SU.resize( f_time_horizon );
   }
   if( ( ! v_DeltaRampDown.empty() ) ) {
    v_T_RD.resize( f_time_horizon );
    v_K_SD.resize( f_time_horizon );
   }

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {

    if( ( ! v_DeltaRampUp.empty() ) )
     v_T_RU[ t ] = std::floor( ( get_operational_max_power( t ) -
                                 v_ShutDownLimit[ t ] ) / v_DeltaRampUp[ t ] );
    if( ( ! v_DeltaRampDown.empty() ) ) {
     v_T_RD[ t ] = std::floor( ( get_operational_max_power( t ) -
                                 v_StartUpLimit[ t ] ) / v_DeltaRampDown[ t ] );
     v_K_SD[ t ] = std::min( f_InitUpDownTime - 1 , v_T_RD[ t ] );
     v_K_SD[ t ] = std::min( static_cast< int >( f_time_horizon - t ) - 1 ,
                             v_K_SD[ t ] );
    }
    if( ( ! v_DeltaRampUp.empty() ) && ( ! v_DeltaRampDown.empty() ) ) {
     v_K_SU[ t ] = std::min( f_InitUpDownTime - 2 -
                             std::max( 0 , v_K_SD[ t ] ) , v_T_RU[ t ] );
     v_K_SU[ t ] = std::min( static_cast< int >( t ) - 1 , v_K_SU[ t ] );
    }
    if( ( ! v_DeltaRampUp.empty() ) && ( v_DeltaRampDown.empty() ) )
     v_K_SU[ t ] = std::min( static_cast< int >( t ) - 1 , v_T_RU[ t ] );
   }

   if( ( ! v_DeltaRampUp.empty() ) ) {
    // size bound constraints 4
    max_power_cnstrs_size += f_time_horizon - init_t;

   // size bound constraints 5
    for( Index t = init_t ; t < f_time_horizon ; ++t )
     if( f_MinUpTime - 2 < v_T_RU[ t ] )
      max_power_cnstrs_size++;
   }
   // size bound constraints 6
   if( ( ! v_DeltaRampUp.empty() ) && ( ! v_DeltaRampDown.empty() ) )
    for( Index t = init_t ; t < f_time_horizon ; ++t )
     if( v_K_SD[ t ] > 0 )
      max_power_cnstrs_size++;
  }

  // size bound constraints 0
  max_power_cnstrs_size += f_time_horizon;

  // size bound constraints 1
  if( f_MinUpTime > 1 )
   max_power_cnstrs_size += f_time_horizon - init_t;
  // size bound constraints 2 + 3
  else
   max_power_cnstrs_size += 2 * f_time_horizon - 2 * init_t;


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

  if( ( ! v_DeltaRampUp.empty() ) || ( ! v_DeltaRampDown.empty() ) ) {

   // Bound constraints 4 - - - - - - - - - - - - - - - - - - - - - - - - - -
  if( ( ! v_DeltaRampUp.empty() ) )
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

    int min_RU = std::min( static_cast< int >( f_MinUpTime ) - 2 , v_T_RU[ t ] );

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
           static_cast< double >( s + 1 ) * v_DeltaRampUp[ t - s ] ) ) );
    }

    MaxPower_Const[ cnstr_idx ].set_lhs( 0.0 );
    MaxPower_Const[ cnstr_idx ].set_rhs( Inf< double >() );
    MaxPower_Const[ cnstr_idx ].set_function(
     new LinearFunction( std::move( vars ) ) );
   }

   // Bound constraints 5 - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ( ! v_DeltaRampUp.empty() ) )
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

     int min_RU = std::min( static_cast< int >( f_MinUpTime ) - 1 , v_T_RU[ t ] );

     if( t >= init_t ) {
      for( int s = 0 ; s < min_RU ; ++s )
       if( t - init_t >= s )
        vars.push_back( std::make_pair(
         &v_start_up[ t - s - init_t ] ,
         -( get_operational_max_power( t - s ) -
            v_StartUpLimit[ t - s ] -
            static_cast< double >( s + 1 ) * v_DeltaRampUp[ t - s ] ) ) );
     }

     MaxPower_Const[ cnstr_idx ].set_lhs( 0.0 );
     MaxPower_Const[ cnstr_idx ].set_rhs( Inf< double >() );
     MaxPower_Const[ cnstr_idx ].set_function(
      new LinearFunction( std::move( vars ) ) );
     cnstr_idx++;
    }

   // Bound constraints 6 - - - - - - - - - - - - - - - - - - - - - - - - - -
  if( ( ! v_DeltaRampUp.empty() ) && ( ! v_DeltaRampDown.empty() ) )
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
  }

  if( ( AR & FormMsk ) == SUForm || ( AR & FormMsk ) == SUSDForm )
  {  // SU and SUSD formulation - - - - - - - - - - - - - - - - - - - - - - -

   if( ( AR & FormMsk ) == SUSDForm )
    MaxPower_Const.resize( v_P_h.size() + v_P_k.size() );
   else
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
  }

  if( ( AR & FormMsk ) == SDForm || ( AR & FormMsk ) == SUSDForm )
  {  // SD and SUSD formulation - - - - - - - - - - - - - - - - - - - - - - -

   if( ( AR & FormMsk ) == SDForm )
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

    if( ( AR & FormMsk ) == SUSDForm ) {
     MaxPower_Const[ v_P_h.size() + j ].set_lhs( 0.0 );
     MaxPower_Const[ v_P_h.size() + j ].set_rhs( Inf< double >() );
     MaxPower_Const[ v_P_h.size() + j ].set_function(
      new LinearFunction( std::move( vars ) ) );
    } else {
     MaxPower_Const[ j ].set_lhs( 0.0 );
     MaxPower_Const[ j ].set_rhs( Inf< double >() );
     MaxPower_Const[ j ].set_function(
      new LinearFunction( std::move( vars ) ) );
    }
   }
  }
 }

 add_static_constraint( MaxPower_Const , "MaxPower_Const_Thermal" );

 if( reserve_vars & 1u )  // if UCBlock has primary demand variables
  if( ! v_PrimaryRho.empty() ) {  // if unit produces any primary reserve
   // initializing primary rho fraction constraints - - - - - - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   PrimaryRho_Const.resize( f_time_horizon );

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {
    vars.push_back( std::make_pair( & v_active_power[ t ] ,
				    v_PrimaryRho[ t ] ) );

    vars.push_back( std::make_pair( & v_primary_spinning_reserve[ t ] ,
				    -1.0 ) );

    PrimaryRho_Const[ t ].set_lhs( 0.0 );
    PrimaryRho_Const[ t ].set_rhs( Inf< double >() );
    PrimaryRho_Const[ t ].set_function(
				   new LinearFunction( std::move( vars ) ) );
    }

   add_static_constraint( PrimaryRho_Const ,
                          "PrimaryRho_Const_Thermal" );
   }

 if( reserve_vars & 2u )  // if UCBlock has secondary demand variables
  if( ! v_SecondaryRho.empty() ) {  // if unit produces any secondary reserve
   // initializing secondary rho fraction constraints - - - - - - - - - - - -
   // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   SecondaryRho_Const.resize( f_time_horizon );

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {
    vars.push_back( std::make_pair( & v_active_power[ t ] ,
				    v_SecondaryRho[ t ] ) );

    vars.push_back( std::make_pair( & v_secondary_spinning_reserve[ t ] ,
                                    -1.0 ) );

    SecondaryRho_Const[ t ].set_lhs( 0.0 );
    SecondaryRho_Const[ t ].set_rhs( Inf< double >() );
    SecondaryRho_Const[ t ].set_function(
				   new LinearFunction( std::move( vars ) ) );
   }

   add_static_constraint( SecondaryRho_Const ,
			  "SecondaryRho_Const_Thermal" );
   }

 // ZOConstraints - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( generate_ZOConstraints ) {
  // the commitment bound constraints
  Commitment_bound_Const.resize( f_time_horizon );

  for( Index t = 0 ; t < f_time_horizon ; ++t )
   Commitment_bound_Const[ t ].set_variable( & v_commitment[ t ] );

  add_static_constraint( Commitment_bound_Const ,
			 "Commitment_bound_Thermal" );

  auto startup_shutdown_size = f_time_horizon - init_t;

  // the startup binary bound constraints
  StartUp_Binary_bound_Const.resize( startup_shutdown_size );

  for( Index t = 0 ; t < startup_shutdown_size ; ++t )
   StartUp_Binary_bound_Const[ t ].set_variable( & v_start_up[ t ] );

  add_static_constraint( StartUp_Binary_bound_Const ,
                         "StartUp_binary_bound_Thermal" );

  // the shut-down binary bound constraints
  ShutDown_Binary_bound_Const.resize( startup_shutdown_size );

  for( Index t = 0 ; t < startup_shutdown_size ; ++t )
   ShutDown_Binary_bound_Const[ t ].set_variable( & v_shut_down[ t ] );

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
  }

  if( ( AR & FormMsk ) == SUForm || ( AR & FormMsk ) == SUSDForm )
  {  // SU and SUSD formulation - - - - - - - - - - - - - - - - - - - - - - -

   if( ( AR & FormMsk ) == SUSDForm )
    Init_PC_Const.resize( 2 * v_P_h.size() + 2 * v_P_k.size() );
   else
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
  }

  if( ( AR & FormMsk ) == SDForm || ( AR & FormMsk ) == SUSDForm )
  {  // SD and SUSD formulation - - - - - - - - - - - - - - - - - - - - - - -

   if( ( AR & FormMsk ) == SDForm )
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

  // Constraints connecting variables of the SUSD formulations with the
  // maximum of the perspective function of the SU and the SD formulations- -
  //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( AR & FormMsk ) == SUSDForm ) {

   Max_SUSD_PC_Const.resize( 2 * f_time_horizon );

   if( ( AR & FormMsk ) == SUSDForm )
    for( Index t = 0 ; t < f_time_horizon ; ++t ) {

     vars.push_back( std::make_pair( &v_cut_teta[ t ] , -1.0 ) );

     for( Index j = 0 ; j < v_Z_h.size() ; ++j )
      if( v_Z_h[ j ].first == t )
       vars.push_back( std::make_pair( &v_cut_h[ j ] , 1.0 ) );

     Max_SUSD_PC_Const[ t ].set_lhs( -Inf< double >() );
     Max_SUSD_PC_Const[ t ].set_rhs( 0.0 );
     Max_SUSD_PC_Const[ t ].set_function(
      new LinearFunction( std::move( vars ) ) );
    }

   if( ( AR & FormMsk ) == SUSDForm )
    for( Index t = 0 ; t < f_time_horizon ; ++t ) {

     vars.push_back( std::make_pair( &v_cut_teta[ t ] , -1.0 ) );

     for( Index j = 0 ; j < v_Z_k.size() ; ++j )
      if( v_Z_k[ j ].first == t )
       vars.push_back( std::make_pair( &v_cut_k[ j ] , 1.0 ) );

     Max_SUSD_PC_Const[ f_time_horizon + t ].set_lhs( -Inf< double >() );
     Max_SUSD_PC_Const[ f_time_horizon + t ].set_rhs( 0.0 );
     Max_SUSD_PC_Const[ f_time_horizon + t ].set_function(
      new LinearFunction( std::move( vars ) ) );
    }
   add_static_constraint( Max_SUSD_PC_Const , "Max_SUSD_PC_Const_Thermal" );
  }

  // Constraints connecting perspective cuts variables of 3bin with those of
  // DP, SU and SD formulations- - - - - - - - - - - - - - - - - - - - - - -
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  if( ( ( AR & FormMsk ) == tbinForm ) ||  // 3bin formulation- - - - - - - -
      ( ( AR & FormMsk ) == TForm ) ||  // T formulation- - - - - - - - - - -
      ( ( AR & FormMsk ) == ptForm ) ) {  // pt formulation - - - - - - - - -
   ;  // does nothing
  } else {  // DP, SU, SD and SUSD formulations - - - - - - - - - - - - - - -
   if( ( AR & FormMsk ) == SUSDForm )
    Eq_PC_Const.resize( 2 * f_time_horizon );
   else
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
   }

   if( ( AR & FormMsk ) == SUForm || ( AR & FormMsk ) == SUSDForm )
   {  // SU and SUSD formulation- - - - - - - - - - - - - - - - - - - - - - -

    for( Index t = 0 ; t < f_time_horizon ; ++t ) {

     vars.push_back( std::make_pair( &v_cut[ t ] , 1.0 ) );

     for( Index j = 0 ; j < v_Z_h.size() ; ++j )
      if( v_Z_h[ j ].first == t )
       vars.push_back( std::make_pair( &v_cut_h[ j ] , -1.0 ) );

     Eq_PC_Const[ t ].set_both( 0.0 );
     Eq_PC_Const[ t ].set_function( new LinearFunction( std::move( vars ) ) );
    }
   }

   if( ( AR & FormMsk ) == SDForm || ( AR & FormMsk ) == SUSDForm )
   {  // SU and SUSD formulation- - - - - - - - - - - - - - - - - - - - - - -

    for( Index t = 0 ; t < f_time_horizon ; ++t ) {

     vars.push_back( std::make_pair( &v_cut[ t ] , 1.0 ) );

     for( Index j = 0 ; j < v_Z_k.size() ; ++j )
      if( v_Z_k[ j ].first == t )
       vars.push_back( std::make_pair( &v_cut_k[ j ] , -1.0 ) );

     if( ( AR & FormMsk ) == SUSDForm ) {
      Eq_PC_Const[ f_time_horizon + t ].set_both( 0.0 );
      Eq_PC_Const[ f_time_horizon + t ].set_function(
       new LinearFunction( std::move( vars ) ) );
     } else {
      Eq_PC_Const[ t ].set_both( 0.0 );
      Eq_PC_Const[ t ].set_function( new LinearFunction( std::move( vars ) ) );
     }
    }
   }
  }

  add_static_constraint( Eq_PC_Const , "Eq_PC_Const_Thermal" );
  }

 if ( ! v_RefSchedule.empty() ) {
  Reference_Schedule_Const.resize( 2 * f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   // | P - Pref | <= v_abs_ref_schedule
   auto lfunc_1 = new LinearFunction();
   lfunc_1->add_variable( & v_active_power[ t ], 1.0 );
   lfunc_1->add_variable( & v_abs_ref_schedule[ t ] , - 1.0 );
   Reference_Schedule_Const[ t ].set_lhs( -Inf< double >() );
   Reference_Schedule_Const[ t ].set_rhs( v_RefSchedule[ t ] );
   Reference_Schedule_Const[ t ].set_function( lfunc_1 );
   //
   auto lfunc_2 = new LinearFunction();
   lfunc_2->add_variable( & v_active_power[ t ], -1.0 );
   lfunc_2->add_variable( & v_abs_ref_schedule[ t ], -1.0 );
   Reference_Schedule_Const[ f_time_horizon + t ].set_lhs( -Inf< double >() );
   Reference_Schedule_Const[ f_time_horizon + t ].set_rhs(
						       - v_RefSchedule[ t ] );
   Reference_Schedule_Const[ f_time_horizon + t ].set_function( lfunc_2 );
   }

  add_static_constraint( Reference_Schedule_Const ,
			 "Norm1_Reference_Schedule" );
  }

 // reactive power bounds constraints (if any) - - - - - - - - - - - - - - -
 //- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

 if( f_reactive_power && 
     ( ( ! v_MinReactivePower.empty() ) || ( ! v_MaxReactivePower.empty() ) )
     ) {
  if( ReactivePower_Bound_Const.empty() )
   ReactivePower_Bound_Const.resize( f_time_horizon );

  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   ReactivePower_Bound_Const[ t ].set_rhs( get_max_reactive_power( t ) );
   ReactivePower_Bound_Const[ t ].set_lhs( get_min_reactive_power( t ) );
   ReactivePower_Bound_Const[ t ].set_variable( & v_reactive_power[ t ] );
   }

  add_static_constraint( ReactivePower_Bound_Const ,
                         "ReactivePowerBound_thermal" );
  }

 set_constraints_generated();

 }  // end( ThermalUnitBlock::generate_abstract_constraints )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::generate_dynamic_constraints( Configuration * dycc )
{
 if( AR & PCuts ) {
  double tol = 1e-7;  // threshold parameter for P/C separation
  double eps = 1e-6;  // tolerance value to consider a binary variable

  bool check_loop = false;
  double part1;

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

  // tbin and T formulations- - - - - - - - - - - - - - - - - - - - - - - - -
  if( ( AR & FormMsk ) == tbinForm || ( AR & FormMsk ) == TForm ) {
   bool check_cut = false;
   double pbar;

   for( Index t = 0 ; t < f_time_horizon ; ++t )

    if( v_commitment[ t ].get_value() > eps ) {

     check_cut = false;

     if( check_loop ) {
      pbar = v_active_power[ t ].get_value() / v_commitment[ t ].get_value();
      double value = ( v_active_power[ t ].get_value() *
       v_active_power[ t ].get_value() ) / v_commitment[ t ].get_value();
      if( v_cut[ t ].get_value() < value - tol )
       if( prevpbar[ t ] == 0 || ( prevpbar[ t ] != 0 && std::abs(
         ( prevpbar[ t ] - pbar ) / prevpbar[ t ] ) > eps ) )
        check_cut = true;
     } else {
      if( v_cut[ t ].get_value() + std::pow(
          v_active_power[ t ].get_value() / v_commitment[ t ].get_value() , 2 ) *
          v_commitment[ t ].get_value() < 1 )
       part1 = 1.0;
      else
       part1 = v_cut[ t ].get_value() + std::pow(
             v_active_power[ t ].get_value() /
               v_commitment[ t ].get_value() , 2 ) * v_commitment[ t ].get_value();

      if( 2.0 * ( v_active_power[ t ].get_value() / v_commitment[ t ].get_value() ) *
          v_active_power[ t ].get_value() - part1 >= tol * part1 )
       check_cut = true;
     }

     if( check_cut ) {
      prevpbar[ t ] = pbar;

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
  }

  // pt formulation - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if( ( AR & FormMsk ) == ptForm ) {
   for( Index t = 0 ; t < f_time_horizon ; ++t ) {
    double sumy = 0.0;
    for( int j = 0 ; j < v_Y_plus.size() ; j++ )
     if( v_Y_plus[ j ].first <= t + 1 && t + 1 <= v_Y_plus[ j ].second ) {
      sumy += v_commitment_plus[ j ].get_value();
     }

    if( sumy > eps ) {
     if( v_cut[ t ].get_value() +
      std::pow( v_active_power[ t ].get_value() / sumy , 2 ) * sumy < 1 )
      part1 = 1.0;
     else
      part1 = v_cut[ t ].get_value() +
       std::pow( v_active_power[ t ].get_value() / sumy , 2 ) * sumy;

     if( 2.0 * ( v_active_power[ t ].get_value() / sumy ) *
                 v_active_power[ t ].get_value() - part1 >= tol * part1 ) {
      std::list< FRowConstraint > cut( 1 );

      vars.push_back( std::make_pair( &v_active_power[ t ] ,
                                      2 * ( v_active_power[ t ].get_value() / sumy ) ) );
      vars.push_back( std::make_pair( &v_cut[ t ] , -1.0 ) );

      for( int j = 0 ; j < v_Y_plus.size() ; j++ )
       if( v_Y_plus[ j ].first <= t + 1 && t + 1 <= v_Y_plus[ j ].second ) {
        vars.push_back(
         std::make_pair( &v_commitment_plus[ j ] ,
                         -( std::pow( v_active_power[ t ].get_value() , 2 ) /
                            std::pow( sumy , 2 ) ) ) );
       }

      cut.front().set_lhs( -Inf< double >() );
      cut.front().set_rhs( 0.0 );
      cut.front().set_function(
       new LinearFunction( std::move( vars ) , eNoMod ) );

      add_dynamic_constraints( PC_cuts , cut , eNoBlck );
     }
    }
   }
  }

  // DP formulation - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if( ( AR & FormMsk ) == DPForm ) {
   for( Index i = 0 ; i < v_P_h_k.size() ; ++i ) {
    for( int j = 0 ; j < v_Y_plus.size() ; j++ )
     if( ( v_Y_plus[ j ].first == v_P_h_k[ i ].second.first ) &&
         ( v_P_h_k[ i ].second.second == v_Y_plus[ j ].second ) )
      if( v_commitment_plus[ j ].get_value() > eps )
       for( Index s = 0 ; s < v_P_h_k.size() ; ++s )
        if( ( v_P_h_k[ i ].first == v_Z_h_k[ s ].first ) &&
            ( v_P_h_k[ i ].second.first == v_Z_h_k[ s ].second.first ) &&
            ( v_P_h_k[ i ].second.second == v_Z_h_k[ s ].second.second ) ) {
         if( v_cut_h_k[ s ].get_value() + std::pow(
           v_active_power_h_k[ i ].get_value() /
             v_commitment_plus[ j ].get_value() , 2 ) *
             v_commitment_plus[ j ].get_value() < 1 )
          part1 = 1.0;
         else
          part1 = v_cut_h_k[ s ].get_value() + std::pow(
            v_active_power_h_k[ i ].get_value() /
            v_commitment_plus[ j ].get_value() , 2 ) *
           v_commitment_plus[ j ].get_value();

         if( 2.0 * ( v_active_power_h_k[ i ].get_value() /
             v_commitment_plus[ j ].get_value() ) *
             v_active_power_h_k[ i ].get_value() - part1 >= tol * part1 ) {
          std::list< FRowConstraint > cut( 1 );

          vars.push_back(
           std::make_pair( &v_active_power_h_k[ i ] ,
            2 * ( v_active_power_h_k[ i ].get_value() /
                  v_commitment_plus[ j ].get_value() ) ) );
          vars.push_back( std::make_pair( &v_cut_h_k[ s ] , -1.0 ) );

          vars.push_back(
           std::make_pair(
            &v_commitment_plus[ j ] ,
            -( std::pow( v_active_power_h_k[ i ].get_value() , 2 ) /
               std::pow( v_commitment_plus[ j ].get_value() , 2 ) ) ) );

          cut.front().set_lhs( -Inf< double >() );
          cut.front().set_rhs( 0.0 );
          cut.front().set_function(
           new LinearFunction( std::move( vars ) , eNoMod ) );

          add_dynamic_constraints( PC_cuts , cut , eNoBlck );
         }
        }
   }
  }

  // SU, SD and SUSD formulations - - - - - - - - - - - - - - - - - - - - - -
  if( ( AR & FormMsk ) == SUForm || ( AR & FormMsk ) == SDForm ||
      ( AR & FormMsk ) == SUSDForm ) {
   if( ( AR & FormMsk ) == SUForm || ( AR & FormMsk ) == SUSDForm )
    for( Index i = 0 ; i < v_P_h.size() ; ++i ) {
     double sumy = 0.0;
     for( int j = 0 ; j < v_Y_plus.size() ; j++ )
      if( ( v_Y_plus[ j ].first == v_P_h[ i ].second ) &&
          ( v_P_h[ i ].first + 1 <= v_Y_plus[ j ].second ) )
       sumy += v_commitment_plus[ j ].get_value();

     if( sumy > eps )
      for( Index s = 0 ; s < v_P_h.size() ; ++s )
       if( ( v_P_h[ i ].first == v_Z_h[ s ].first ) &&
           ( v_P_h[ i ].second == v_Z_h[ s ].second ) )
        if( v_Z_h[ s ].first == v_P_h[ i ].first ) {
         if( v_cut_h[ s ].get_value() + std::pow(
           v_active_power_h[ i ].get_value() / sumy , 2 ) * sumy < 1 )
          part1 = 1.0;
         else
          part1 = v_cut_h[ s ].get_value() +
           std::pow( v_active_power_h[ i ].get_value() / sumy , 2 ) * sumy;

         if( 2.0 * ( v_active_power_h[ i ].get_value() / sumy ) *
             v_active_power_h[ i ].get_value() - part1 >= tol * part1 ) {
          std::list< FRowConstraint > cut( 1 );

          vars.push_back(
           std::make_pair( &v_active_power_h[ i ] ,
                           2 * ( v_active_power_h[ i ].get_value() / sumy ) ) );
          vars.push_back( std::make_pair( &v_cut_h[ s ] , -1.0 ) );

          for( int j = 0 ; j < v_Y_plus.size() ; j++ )
           if( ( v_Y_plus[ j ].first == v_P_h[ i ].second ) &&
               ( v_P_h[ i ].first + 1 <= v_Y_plus[ j ].second ) ) {
            vars.push_back(
             std::make_pair( &v_commitment_plus[ j ] ,
                             -( std::pow( v_active_power_h[ i ].get_value() , 2 ) /
                                std::pow( sumy , 2 ) ) ) );
           }

          cut.front().set_lhs( -Inf< double >() );
          cut.front().set_rhs( 0.0 );
          cut.front().set_function(
           new LinearFunction( std::move( vars ) , eNoMod ) );

          add_dynamic_constraints( PC_cuts , cut , eNoBlck );
         }
        }
    }

   // SD and SUSD formulations- - - - - - - - - - - - - - - - - - - - - - - -
   if( ( AR & FormMsk ) == SDForm || ( AR & FormMsk ) == SUSDForm )
    for( Index i = 0 ; i < v_P_k.size() ; ++i ) {
     double sumy = 0.0;
     for( int j = 0 ; j < v_Y_plus.size() ; j++ )
      if( ( v_Y_plus[ j ].second == v_P_k[ i ].second ) &&
          ( v_P_k[ i ].first + 1 >= v_Y_plus[ j ].first ) )
       sumy += v_commitment_plus[ j ].get_value();

     if( sumy > eps )
      for( Index s = 0 ; s < v_P_k.size() ; ++s )
       if( ( v_P_k[ i ].first == v_Z_k[ s ].first ) &&
           ( v_P_k[ i ].second == v_Z_k[ s ].second ) )
        if( v_Z_k[ s ].first == v_P_k[ i ].first ) {
         if( v_cut_k[ s ].get_value() +
          std::pow( v_active_power_k[ i ].get_value() / sumy , 2 ) * sumy < 1 )
          part1 = 1.0;
         else
          part1 = v_cut_k[ s ].get_value() +
           std::pow( v_active_power_k[ i ].get_value() / sumy , 2 ) * sumy;

         if( 2.0 * ( v_active_power_k[ i ].get_value() / sumy ) *
          v_active_power_k[ i ].get_value() - part1 >= tol * part1 ) {

          std::list< FRowConstraint > cut( 1 );

          vars.push_back(
           std::make_pair( &v_active_power_k[ i ] ,
                           2 * ( v_active_power_k[ i ].get_value() / sumy ) ) );
          vars.push_back( std::make_pair( &v_cut_k[ s ] , -1.0 ) );

          for( int j = 0 ; j < v_Y_plus.size() ; j++ )
           if( ( v_Y_plus[ j ].second == v_P_k[ i ].second ) &&
               ( v_P_k[ i ].first + 1 >= v_Y_plus[ j ].first ) ) {
            vars.push_back(
             std::make_pair( &v_commitment_plus[ j ] ,
                             -( std::pow( v_active_power_k[ i ].get_value() , 2 ) /
                                std::pow( sumy , 2 ) ) ) );
           }

          cut.front().set_lhs( -Inf< double >() );
          cut.front().set_rhs( 0.0 );
          cut.front().set_function(
           new LinearFunction( std::move( vars ) , eNoMod ) );

          add_dynamic_constraints( PC_cuts , cut , eNoBlck );
         }
        }
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
 // - then possibly the perspective-cut and reactive power variables
 //
 // - and FINALLY, if present, the single design (investment) variable
 //
 // this arrangement is exploited in add_Modification to easily map
 // indices in the coefficients of the Objective Function back into
 // indices of the original variables (and figure out the kind of variable).
 // The design variable is kept last so that the time-indexed sections above
 // start at index 0 with no leading offset

 if( v_commitment.size() != f_time_horizon )
  throw( std::logic_error( "ThermalUnitBlock::generate_objective: "
			   "v_commitment must have time horizon" ) );

 if( v_active_power.size() != f_time_horizon )
  throw( std::logic_error( "ThermalUnitBlock::generate_objective: "
			   "v_active_power must have time horizon" ) );

 if( v_start_up.size() != f_time_horizon - init_t )
  throw( std::logic_error( "ThermalUnitBlock::generate_objective: "
			   "v_start_up must have size horizon - init_t" ) );

 DQuadFunction::v_coeff_triple vars;

 // NOTE: the design (investment) variable, if present, is added LAST (after
 // the reactive power variables); see the comment there. This keeps the
 // index mapping done by handle_objective_change() for all the other
 // (time-indexed) variable sections free of any leading offset.

 // add the start-up variables- - - - - - - - - - - - - - - - - - - - - - - -
 // add start-up variables for tbin and T formulations
 //if( ( AR & FormMsk ) == tbinForm || ( AR & FormMsk ) == TForm )
 for( Index t = init_t ; t < f_time_horizon ; ++t )
  vars.push_back( std::make_tuple( & v_start_up[ t - init_t ] ,
                                   f_scale * v_StartUpCost[ t ] , 0 ) );
 // add start-up variables for pt, DP, SU, SD and SUSD formulations
 /*
 if( ( AR & FormMsk ) == ptForm || ( AR & FormMsk ) == DPForm ||
     ( AR & FormMsk ) == SUForm || ( AR & FormMsk ) == SDForm ||
     ( AR & FormMsk ) == SUSDForm )
   for( Index j = 0 ; j < v_Y_minus.size() ; ++j )
   if( v_Y_minus[ j ].second != f_time_horizon + 1 )
    vars.push_back( std::make_tuple( &v_commitment_minus[ j ] ,
                                     f_scale *
                                     v_StartUpCost[ v_Y_minus[ j ].second - 1 ] ,
                                     0.0 ) );
 */

 // add the active power variables- - - - - - - - - - - - - - - - - - - - - -
 for( Index t = 0 ; t < f_time_horizon ; ++t )
  vars.push_back( std::make_tuple( & v_active_power[ t ] ,
                                   f_scale * v_LinearTerm[ t ] ,
                                   AR & PCuts ? 0 : f_scale * v_QuadTerm[ t ] ) );

 // add the commitment variables- - - - - - - - - - - - - - - - - - - - - - -
 for( Index t = 0 ; t < f_time_horizon ; ++t )
  vars.push_back( std::make_tuple( &v_commitment[ t ] ,
                                   f_scale * v_ConstTerm[ t ] , 0 ) );

 if ( ! v_RefSchedule.empty() )
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   vars.push_back( std::make_tuple( & v_abs_ref_schedule[ t ] , 1 , 0 ) );

 if( ( reserve_vars & 1u ) && ( ! v_primary_spinning_reserve.empty() ) ) {
  // add the primary spinning reserve variables - - - - - - - - - - - - - - -
  if( v_primary_spinning_reserve.size() != f_time_horizon )
   throw( std::logic_error( "ThermalUnitBlock::generate_objective: v_primary_"
                            "spinning_reserve must have size equal to the "
                            "time horizon." ) );

  // the reserve cost coefficient defaults to the participation factor rho,
  // but lives in a separate, modifiable vector (see v_PrimaryRho)
  if( v_PrimarySpinningReserveCost.empty() )
   v_PrimarySpinningReserveCost = v_PrimaryRho;

  if( v_PrimarySpinningReserveCost.empty() )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    vars.push_back( std::make_tuple( &v_primary_spinning_reserve[ t ] ,
                                     0.0 , 0.0 ) );
  else
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    vars.push_back( std::make_tuple( &v_primary_spinning_reserve[ t ] ,
                                     f_scale * v_PrimarySpinningReserveCost[ t ] ,
                                     0.0 ) );
 }

 if( ( reserve_vars & 2u ) && ( ! v_secondary_spinning_reserve.empty() ) ) {
  // add the secondary spinning reserve variables - - - - - - - - - - - - - -
  if( v_secondary_spinning_reserve.size() != f_time_horizon )
   throw( std::logic_error( "ThermalUnitBlock::generate_objective: v_secondary"
                            "_spinning_reserve must have size equal to the "
                            "time horizon." ) );

  if( v_SecondarySpinningReserveCost.empty() )
   v_SecondarySpinningReserveCost = v_SecondaryRho;

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

 if( AR & PCuts ) {
  // add the perspective cuts variables - - - - - - - - - - - - - - - - - - -
  // tbin, T and pt formulations- - - - - - - - - - - - - - - - - - - - - - -
  if( ( AR & FormMsk ) == tbinForm || ( AR & FormMsk ) == TForm ||
   ( AR & FormMsk ) == ptForm )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    vars.push_back( std::make_tuple( &v_cut[ t ] ,
                                     f_scale * v_QuadTerm[ t ] , 0.0 ) );
  // DP formulation - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if( ( AR & FormMsk ) == DPForm )
   for( Index i = 0 ; i < v_Z_h_k.size() ; ++i )
    vars.push_back( std::make_tuple( &v_cut_h_k[ i ] ,
                                     f_scale *
                                     v_QuadTerm[ v_Z_h_k[ i ].first ] , 0.0 ) );
  // SU formulation - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if( ( AR & FormMsk ) == SUForm )
   for( Index i = 0 ; i < v_Z_h.size() ; ++i )
    vars.push_back( std::make_tuple( &v_cut_h[ i ] ,
                                     f_scale * v_QuadTerm[ v_Z_h[ i ].first ] ,
                                     0.0 ) );
  // SD formulation - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if( ( AR & FormMsk ) == SDForm )
   for( Index i = 0 ; i < v_Z_k.size() ; ++i )
    vars.push_back( std::make_tuple( &v_cut_k[ i ] ,
                                     f_scale * v_QuadTerm[ v_Z_k[ i ].first ] ,
                                     0.0 ) );
  // SUSD formulation - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if( ( AR & FormMsk ) == SUSDForm )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    vars.push_back( std::make_tuple( &v_cut_teta[ t ] ,
                                     f_scale * v_QuadTerm[ t ] , 0.0 ) );
  }

 // add the reactive power variables LAST - - - - - - - - - - - - - - - - - -
 // q[t] carries no cost in the unit objective, so its coefficient is 0 unless
 // a dualizing Solver sets v_ReactiveLinearTerm (see set_reactive_linear_term):
 // keeping them last lets handle_objective_change() recognise a reactive
 // coefficient change as the trailing [ num_active_var - th , num_active_var )
 // index block, independent of which optional sections precede it.
 if( f_reactive_power )
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   vars.push_back( std::make_tuple(
    & v_reactive_power[ t ] ,
    f_scale * ( v_ReactiveLinearTerm.empty() ? 0.0 : v_ReactiveLinearTerm[ t ] ) ,
    0.0 ) );

 // add the design (investment) variable LAST - - - - - - - - - - - - - - - -
 // x carries the (per-module, scaled) investment cost; a dualizing Solver may
 // change this coefficient (e.g. the non-anticipativity multiplier in a nested
 // Lagrangian). Keeping it as the trailing, single-variable block lets
 // handle_objective_change() recognise a design-coefficient change as the last
 // index num_active_var - 1, peel it off, and route it to
 // update_objective_investment() (which issues a eSetInvCost
 // ThermalUnitBlockMod for the DP Solvers), independent of all the optional
 // sections that may precede it.
 if( f_InvestmentCost != 0 ) {
  vars.push_back( std::make_tuple( & design ,
                                   f_scale * f_InvestmentCost , 0 ) );
  // seed the change-detection cache with the initial design cost, so that the
  // first (no-op) rewrite of the coefficient vector by a dualizing Solver does
  // not issue a useless eSetInvCost (the DP Solvers read the initial value
  // directly at setup via get_design_cost()); unscaled, to match get_design_cost()
  f_last_design_cost = f_InvestmentCost;
  }

 objective.set_function( new DQuadFunction( std::move( vars ) ) );
 objective.set_sense( Objective::eMin );

 // set Block objective
 set_objective( &objective , eNoMod );

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
  && ColVariable::is_feasible( v_reactive_power , tol )
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
  && RowConstraint::is_feasible( Commitment_fixed_to_One_Const , tol ,
				 rel_viol )
  && RowConstraint::is_feasible( Reference_Schedule_Const , tol , rel_viol )
  && RowConstraint::is_feasible( ReactivePower_Bound_Const , tol
				 , rel_viol ) );

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

 if( f_scale != 1 )
  ::serialize( group , "Scale" , netCDF::NcDouble() , f_scale );

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
 // the spinning-reserve linear costs are, like LinearTerm, the Lagrangian dual
 // of a relaxed system constraint (here the reserve demand) pushed onto the
 // unit; serialise them too so a dumped unit carries its reserve prices (empty,
 // hence skipped, when the unit offers no reserve or none is priced)
 serialize( "PrimarySpinningReserveCost" , v_PrimarySpinningReserveCost );
 serialize( "SecondarySpinningReserveCost" , v_SecondarySpinningReserveCost );
 serialize( "ConstTerm" , v_ConstTerm );
 serialize( "StartUpCost" , v_StartUpCost );
 serialize( "FixedConsumption" , v_FixedConsumption );
 serialize( "InertiaCommitment" , v_InertiaCommitment );
 serialize( "StartUpLimit" , v_StartUpLimit );
 serialize( "ShutDownLimit" , v_ShutDownLimit );

}  // end( ThermalUnitBlock::serialize )

/*--------------------------------------------------------------------------*/
/*----------------------- Methods for handling Solution --------------------*/
/*--------------------------------------------------------------------------*/

Solution * ThermalUnitBlock::get_Solution( Configuration * csolc ,
             bool emptys )
{
 Index wsol = 15;
 if( ( ! csolc ) && f_BlockConfig )
  csolc = f_BlockConfig->f_solution_Configuration;

 if( auto config = dynamic_cast< SimpleConfiguration< int > * >( csolc ) )
  wsol = config->f_value;

 // call the method of the base class
 auto * sol = dynamic_cast< ThermalUnitBlockSolution * >(
                     UnitBlock::get_Solution( csolc , emptys ) );
 assert( sol );

 if( ! emptys )
  sol->read( this );

 return( sol );
 }

/*--------------------------------------------------------------------------*/

UnitBlockSolution * ThermalUnitBlock::new_Solution( void ) const {
 return( new ThermalUnitBlockSolution() );
 }

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
   v_Availability[ t ] = *( availability++ );
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
  return;                                              // return

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
 rng.second = std::min( rng.second , static_cast< Index >( 1 ) );
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
  return;                                                 // return

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
  return;                                               // return

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
  return;                                                // return

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
// the reactive power variables q[t] are the LAST section of the Objective
// (see generate_objective): they are added only when f_reactive_power, after
// every other section, so they occupy the active-variable indices
// [ num_active_var - f_time_horizon , num_active_var ). The DP solvers price
// q[t] over [Qmin,Qmax] using this coefficient.

void ThermalUnitBlock::set_reactive_linear_term( MF_dbl_it values ,
                                                 Subset && subset ,
                                                 const bool ordered ,
                                                 ModParam issuePMod ,
                                                 ModParam issueAMod )
{
 if( subset.empty() )
  return;

 if( v_ReactiveLinearTerm.empty() ) {
  if( std::all_of( values ,
                   values + subset.size() ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  v_ReactiveLinearTerm.assign( f_time_horizon , 0 );
 }

 if( ! ordered )
  std::sort( subset.begin() , subset.end() );

 if( subset.back() >= v_ReactiveLinearTerm.size() )
  throw( std::invalid_argument(
   "ThermalUnitBlock::set_reactive_linear_term: invalid index in subset." ) );

 if( identical( v_ReactiveLinearTerm , subset , values ) )  // nothing changes
  return;

 if( not_dry_run( issuePMod ) )
  assign( v_ReactiveLinearTerm , subset , values );

 if( not_dry_run( issueAMod ) && objective_generated() && f_reactive_power ) {
  auto * qf = static_cast< DQuadFunction * >( objective.get_function() );
  const Index dpos = qf->get_num_active_var() - f_time_horizon;
  Subset tmps = subset_add( subset , dpos );
  DQuadFunction::Vec_FunctionValue tmpv( values , values + subset.size() );
  qf->modify_linear_coefficients( std::move( tmpv ) , std::move( tmps ) ,
                                  true , un_ModBlock( issueAMod ) );
 }

 if( issue_pmod( issuePMod ) )
  Block::add_Modification( std::make_shared< ThermalUnitBlockSbstMod >(
                            this , ThermalUnitBlockMod::eSetReactiveLinT ,
                            std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_reactive_linear_term( subset ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_reactive_linear_term( MF_dbl_it values ,
                                                 Range rng ,
                                                 ModParam issuePMod ,
                                                 ModParam issueAMod )
{
 rng.second = std::min( rng.second , f_time_horizon );
 if( rng.second <= rng.first )
  return;

 c_Index sz = rng.second - rng.first;
 if( v_ReactiveLinearTerm.empty() ) {
  if( std::all_of( values ,
                   values + sz ,
                   []( double cst ) { return( cst == 0 ); } ) )
   return;

  v_ReactiveLinearTerm.assign( f_time_horizon , 0 );
 }

 if( std::equal( values , values + sz ,
                 v_ReactiveLinearTerm.begin() + rng.first ) )
  return;

 if( not_dry_run( issuePMod ) )
  std::copy( values , values + sz , v_ReactiveLinearTerm.begin() + rng.first );

 if( not_dry_run( issueAMod ) && objective_generated() && f_reactive_power ) {
  auto * qf = static_cast< DQuadFunction * >( objective.get_function() );
  const Index dpos = qf->get_num_active_var() - f_time_horizon;
  DQuadFunction::Vec_FunctionValue tmpv( values , values + sz );
  qf->modify_linear_coefficients( std::move( tmpv ) ,
                                  Range( rng.first + dpos , rng.second + dpos ) ,
                                  un_ModBlock( issueAMod ) );
 }

 if( issue_pmod( issuePMod ) )
  Block::add_Modification( std::make_shared< ThermalUnitBlockRngdMod >(
                            this , ThermalUnitBlockMod::eSetReactiveLinT ,
                            rng ) ,
                           Observer::par2chnl( issuePMod ) );

}  // end( ThermalUnitBlock::set_reactive_linear_term( range ) )

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
  return;                                              // return

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
     *( tmplvit++ ) = v_LinearTerm[ t ];
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
 rng.second = std::min( rng.second , static_cast< decltype( rng.second ) >( 1 ) );
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
    update_objective( Range( 0 , Inf< Index >() ) , issuePMod , issueAMod );
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
                                                  c_ModParam issueAMod ) const
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
                                                      c_ModParam issueAMod ) const
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
                                                    c_ModParam issueAMod ) const
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

double ThermalUnitBlock::get_design_cost( void ) const
{
 if( f_InvestmentCost == 0 )
  return( 0 );  // no design variable: no design cost

 if( objective_generated() )
  if( auto function = dynamic_cast< DQuadFunction * >(
                                                  objective.get_function() ) ) {
   auto var_index = function->is_active( & design );
   if( var_index < function->get_num_active_var() )
    // the current (possibly dualized) coefficient of the design variable;
    // the Objective stores f_scale * cost, so divide it out to return the
    // cost in the same unscaled units as get_investment_cost()
    return( function->get_linear_coefficient( var_index ) / f_scale );
   }

 // the Objective is not available: fall back to the "original" cost
 return( f_InvestmentCost );

}  // end( ThermalUnitBlock::get_design_cost )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::update_objective_investment( ModParam issuePMod ,
                                                    ModParam issueAMod )
{
 if( ! objective_generated() )
  return;  // the Objective has not been generated: nothing to be done

 if( f_InvestmentCost == 0 )
  return;  // no design term in the Objective: nothing to be done

 if( not_dry_run( issueAMod ) ) {
  // change the abstract representation, i.e., the coefficient of the design
  // variable in the Objective
  auto function = dynamic_cast< DQuadFunction * >( objective.get_function() );
  if( ! function )
   return;
  auto var_index = function->is_active( & design );
  assert( var_index < function->get_num_active_var() );
  function->modify_linear_coefficient( var_index ,
                                       f_scale * f_InvestmentCost ,
                                       un_ModBlock( issueAMod ) );
  }

 // only issue the physical Modification if the design cost actually changed
 // since the last one was issued: a dualizing Solver rewrites the whole
 // coefficient vector (design included) at every iteration, but the design
 // coefficient itself changes only when the dual multiplier on it does. Issuing
 // a eSetInvCost every time would force the dualizing Bundle to invalidate this
 // component's linearizations at each iteration (preventing convergence) and
 // flood the (DP) Solvers with useless global-pool rechecks.
 const double cur = get_design_cost();
 if( ! ( cur == f_last_design_cost ) ) {  // ( cur != cached ), NaN-safe
  f_last_design_cost = cur;

  if( issue_pmod( issuePMod ) )
   // Issue a Physical Modification so that Solvers that consume the structural
   // data (the DP Solvers) refresh their copy of the design cost. Note that the
   // design cost lives only in the abstract Objective (there is no separate
   // physical field), so the Solvers re-read it via get_design_cost()
   Block::add_Modification( std::make_shared< ThermalUnitBlockMod >(
                             this , ThermalUnitBlockMod::eSetInvCost ) ,
                            Observer::par2chnl( issuePMod ) );
  }

}  // end( ThermalUnitBlock::update_objective_investment )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::update_objective( const Subset & subset ,
                                         ModParam issuePMod ,
                                         c_ModParam issueAMod ) {
 update_objective_investment( issuePMod , issueAMod );
 update_objective_start_up( subset , issueAMod );
 update_objective_active_power( subset , issueAMod );
 update_objective_commitment( subset , issueAMod );
}  // end( ThermalUnitBlock::update_objective( subset ) )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::update_objective( Range rng ,
                                         ModParam issuePMod ,
                                         c_ModParam issueAMod ) {
 rng.second = std::min( rng.second , f_time_horizon );
 if( rng.second <= rng.first )
  return;

 Subset subset( rng.second - rng.first );
 std::iota( subset.begin() , subset.end() , rng.first );

 update_objective( subset , issuePMod , issueAMod );
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
  // THE Modification SHOULD NOT BE ISSUED WHEN THEY ARE CREATED!
  //if( ! variables_generated() )
  // return;

  throw( std::logic_error( "ThermalUnitBlock - VariableMod not supported" ) );
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
  // THE Modification SHOULD NOT BE ISSUED WHEN IT IS CREATED!
  //if( ! objective_generated() )
  // return;

  throw( std::logic_error( "ThermalUnitBlock - BlockMod not supported." ) );

  // TODO: BlockMod - obj changed
  return;
  }

 // FunctionMod - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( const auto tmod = dynamic_cast< FunctionMod * >( mod ) ) {
  auto f = tmod->function();
  if( f == static_cast< FRealObjective * >( get_objective()
					    )->get_function() ) {
   handle_objective_change( tmod , chnl );
   return;
   }

  std::ostringstream em;
  em << *mod;
  throw( std::invalid_argument( "ThermalUnitBlock::add_Modification: "
				"unsupported " + em.str() ) );
  return;
  }

 // any other Modification is not supported - - - - - - - - - - - - - - - - -

 std::ostringstream em;
 em << *mod;
 throw( std::invalid_argument( "ThermalUnitBlock::add_Modification: "
			       "unsupported " + em.str() ) );

 }  // end( ThermalUnitBlock::guts_of_add_Modification )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::set_solution( void )
{
 // canonical part: the caller has already set v_active_power[t] and
 // v_commitment[t] for all t. We read them back to derive the
 // formulation-specific auxiliaries below.
 auto Pi = get_const_active_power( 0 );
 auto Ci = get_const_commitment( 0 );
 if( ( ! Pi ) || ( ! Ci ) )
  return;             // no canonical variables: nothing to derive from

 // start_up[ t ] = 1 iff commitment goes off->on at t. start_up is
 // indexed from init_t onwards (size = time_horizon - init_t); the
 // boundary case t == init_t == 0 is handled via the pre-horizon
 // state in f_InitUpDownTime: if the unit was off before t = 0
 // (f_InitUpDownTime <= 0) and is on at t = 0, that counts as a
 // start-up at t = 0; otherwise start_up[ 0 ] = 0
 if( auto sup_it = get_start_up() ) {
  if( init_t == 0 )
   sup_it[ 0 ].set_value(
    ( ( f_InitUpDownTime <= 0 ) && ( Ci[ 0 ].get_value() > 0.5 ) )
    ? 1.0 : 0.0 );
  for( Index t = std::max( init_t , Index( 1 ) ) ;
       t < f_time_horizon ; ++t )
   sup_it[ t - init_t ].set_value(
    ( ( Ci[ t ].get_value() > 0.5 ) && ( Ci[ t - 1 ].get_value() <= 0.5 ) )
    ? 1.0 : 0.0 );
  }

 // shut_down[ t ] = 1 iff commitment goes on->off at t (symmetric)
 if( auto sdn_it = get_shut_down() ) {
  if( init_t == 0 )
   sdn_it[ 0 ].set_value(
    ( ( f_InitUpDownTime > 0 ) && ( Ci[ 0 ].get_value() <= 0.5 ) )
    ? 1.0 : 0.0 );
  for( Index t = std::max( init_t , Index( 1 ) ) ;
       t < f_time_horizon ; ++t )
   sdn_it[ t - init_t ].set_value(
    ( ( Ci[ t ].get_value() <= 0.5 ) && ( Ci[ t - 1 ].get_value() > 0.5 ) )
    ? 1.0 : 0.0 );
  }

 // perspective-cut auxiliary variables (only when PCuts is active).
 // The cost coefficient of v_cut[t] in the Objective is alpha_t =
 // f_scale * v_QuadTerm[t]; the perspective constraint v_cut >= p^2 / u
 // is tight at the integer optimum (u in {0,1}), giving alpha_t * p_t^2
 // -- the original quadratic at integer u. Setting v_cut to the same
 // value here keeps LagBFunction's "original cost at x*" recomputation
 // (which reads the variable value via the saved CostMatrix) consistent
 // with what the formulation that uses the original quadratic produces.
 if( has_perspective_cuts() ) {
  const auto form = get_formulation();
  if( ( form == tbinForm ) || ( form == TForm ) || ( form == ptForm ) ) {
   if( auto cut_it = get_cut() )
    for( Index t = 0 ; t < f_time_horizon ; ++t ) {
     double pt = Pi[ t ].get_value();
     bool ut = Ci[ t ].get_value() > 0.5;
     cut_it[ t ].set_value( ut ? pt * pt : 0.0 );
     }
   }
  // TODO: DPForm / SUForm / SDForm / SUSDForm need to populate
  // v_cut_h_k / v_cut_h / v_cut_k / v_cut_teta indexed by the
  // disaggregated graph; the value at the active arc is p_t^2 and 0
  // elsewhere. Those formulations also need v_active_power_h_k,
  // v_commitment_plus, etc., which are not yet handled here.
  }

 }  // end( ThermalUnitBlock::set_solution )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlock::handle_objective_change( FunctionMod * mod ,
                                                ChnlName chnl )
{
 const auto * qf = static_cast< const DQuadFunction * >( mod->function() );
 auto par = make_par( eNoBlck , chnl );
 // the method is implemented by calling the physical change methods
 // set_startup_costs(), set_linear_term() etc.; these need be called with
 // eNoBlck for PMod (the physical representation need be changed, but the
 // corresponding Modification has to be ignored by the ThermalUnitBlock
 // since it has been self-inflicted), and with eDryRun for AMod (the
 // abstract representation has been changed already)
 Index th = f_time_horizon;

 // C05FunctionModLinRngd / DQuadFunctionModRngd - - - - - - - - - - - - - - -
 // split the Modification in up to 5 physical Modification by calling the
 // appropriate set_*() methods (ranged version) for those among startup,
 // power, commitment, primary/secondary reserve variables whose coefficient
 // change. This heavily relies on the fact that variables of the same type
 // are consecutive (and ordered in the obvious way) when set as coefficients
 // in the Objective. C05FunctionModLinRngd only reports changes to linear
 // coefficients, while DQuadFunctionModRngd also reports changes to quadratic
 // coefficients; since only the active power variables carry a non-zero
 // quadratic coefficient, set_quad_term() is called exclusively for that
 // section, and only when the Modification is a DQuadFunctionModRngd (the
 // with_quad flag).

 const Range * rng = nullptr;
 bool with_quad = false;
 if( const auto tmod = dynamic_cast< C05FunctionModLinRngd * >( mod ) )
  rng = & tmod->range();
 else
  if( const auto tmod = dynamic_cast< DQuadFunctionModRngd * >( mod ) ) {
   rng = & tmod->range();
   with_quad = true;
   }

 if( rng ) {
  Index l = rng->first;
  Index r = rng->second;

  if( r > qf->get_num_active_var() )
   throw( std::invalid_argument(
		  "ThermalUnitBlock::add_Modification: invalid Range [" +
		  std::to_string( l ) + ", " + std::to_string( r ) + ")" ) );

  // the design (investment) variable, if present, is the trailing single-var
  // block at index num_active_var - 1: peel off a change to it and route it to
  // update_objective_investment(), which (re)issues a eSetInvCost
  // ThermalUnitBlockMod for the DP Solvers (the abstract Objective has been
  // changed already, hence eDryRun). Afterwards work with the design-free count
  // nav, which makes the time-indexed sections below start at index 0.
  const Index has_design = ( f_InvestmentCost != 0 ) ? Index( 1 ) : Index( 0 );
  const Index nav = qf->get_num_active_var() - has_design;
  if( has_design && ( r > nav ) ) {
   update_objective_investment( par , eDryRun );
   r = nav;
   if( r <= l )
    return;  // only the design coefficient changed
   }

  // reactive power variables are the trailing section: peel off any index in
  // [ q_start , nav ) and route it to set_reactive_linear_term(),
  // then let the code below handle the remaining (non-reactive) prefix only.
  if( f_reactive_power ) {
   const Index q_start = nav - th;
   if( r > q_start ) {
    const Index rl = std::max( l , q_start );
    std::vector< double > qv( r - rl );
    auto qvit = qv.begin();
    for( Index i = rl ; i < r ; )
     *( qvit++ ) = qf->get_linear_coefficient( i++ );
    set_reactive_linear_term( qv.begin() ,
                              Range( rl - q_start , r - q_start ) ,
                              par , eDryRun );
    if( l >= q_start )
     return;  // the whole range was reactive
    r = q_start;
    }
   }

  std::vector< double > nv( r - l + 1 );
  std::vector< double > nvq;
  if( with_quad )
   nvq.resize( r - l + 1 );
  Index gl = 0;
  Index gr = th - init_t;

  if( l < gr ) {  // startup variables
   Index r2 = std::min( r , gr );
   auto nvit = nv.begin();
   for( Index i = l ; i < r2 ; )
    *( nvit++ ) = qf->get_linear_coefficient( i++ );
   // set_startup_costs( Range ) expects rng in time-space [ init_t ,
   // f_time_horizon ); active-var indices [ l , r2 ) correspond to time
   // indices [ l + init_t , r2 + init_t )
   set_startup_costs( nv.begin() ,
                      Range( l + init_t , r2 + init_t ) ,
                      par , eDryRun );
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
    *( nvit++ ) = qf->get_linear_coefficient( i++ );
   set_linear_term( nv.begin() , Range( l - gl , r2 - gl ) , par , eDryRun );
   if( with_quad ) {
    auto nvqit = nvq.begin();
    for( Index i = l ; i < r2 ; )
     *( nvqit++ ) = qf->get_quadratic_coefficient( i++ );
    set_quad_term( nvq.begin() , Range( l - gl , r2 - gl ) , par , eDryRun );
    }

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
    *( nvit++ ) = qf->get_linear_coefficient( i++ );
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
    *( nvit++ ) = qf->get_linear_coefficient( i++ );
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
    *( nvit++ ) = qf->get_linear_coefficient( i++ );
   set_secondary_spinning_reserve_cost( nv.begin() , Range( l - gl , r2 - gl ) ,
                                        par , eDryRun );
   if( r2 == r )
    return;
   }

  throw( std::invalid_argument( "ThermalUnitBlock::add_Modification: invalid "
				"variable in Modification" ) );
  return;

  }  // end( Rngd )

 // C05FunctionModLinSbst / DQuadFunctionModSbst - - - - - - - - - - - - - - -
 // split the Modification in up to 5 physical Modification by calling the
 // appropriate set_*() methods (subset version); see the comment on the
 // ranged case above: the same considerations apply, with the quadratic
 // coefficient being handled only for the active power section and only
 // when the Modification is a DQuadFunctionModSbst.

 const Subset * sbs = nullptr;
 with_quad = false;
 if( const auto tmod = dynamic_cast< C05FunctionModLinSbst * >( mod ) )
  sbs = & tmod->subset();
 else
  if( const auto tmod = dynamic_cast< DQuadFunctionModSbst * >( mod ) ) {
   sbs = & tmod->subset();
   with_quad = true;
   }

 Subset des_reduced;  // storage when the trailing design index is peeled off
 Subset reduced;  // storage when the reactive tail has to be peeled off
 if( sbs ) {
  if( sbs->back() > qf->get_num_active_var() )
   throw( std::invalid_argument( "ThermalUnitBlock::add_Modification: "
				 "invalid Subset" ) );

  // the design (investment) variable, if present, is the trailing single-var
  // block at index num_active_var - 1 (i.e., nav): peel a change to it off the
  // (sorted) Subset and route it to update_objective_investment() (see the
  // ranged case). Afterwards work with the design-free count nav.
  const Index has_design = ( f_InvestmentCost != 0 ) ? Index( 1 ) : Index( 0 );
  const Index nav = qf->get_num_active_var() - has_design;
  if( has_design && ( ! sbs->empty() ) && ( sbs->back() >= nav ) ) {
   update_objective_investment( par , eDryRun );
   des_reduced.assign( sbs->begin() , sbs->end() - 1 );
   sbs = & des_reduced;
   if( sbs->empty() )
    return;  // only the design coefficient changed
   }

  // peel off the trailing reactive section (see the ranged case): entries in
  // [ q_start , nav ) go to set_reactive_linear_term(); the
  // remaining prefix is handled by the section walk below unchanged.
  if( f_reactive_power ) {
   const Index q_start = nav - th;
   auto qit = std::lower_bound( sbs->begin() , sbs->end() , q_start );
   if( qit != sbs->end() ) {
    Subset qms( std::distance( qit , sbs->end() ) );
    std::vector< double > qv( qms.size() );
    auto qvit = qv.begin();
    auto qmsit = qms.begin();
    for( auto it = qit ; it != sbs->end() ; ++it ) {
     *( qvit++ ) = qf->get_linear_coefficient( *it );
     *( qmsit++ ) = *it - q_start;
     }
    set_reactive_linear_term( qv.begin() , std::move( qms ) , true ,
                              par , eDryRun );
    if( qit == sbs->begin() )
     return;  // the whole subset was reactive
    reduced.assign( sbs->begin() , qit );
    sbs = & reduced;
    }
   }

  std::vector< double > nv( sbs->size() );
  std::vector< double > nvq;
  if( with_quad )
   nvq.resize( sbs->size() );
  auto l = sbs->begin();
  Index gl = 0;
  Index gr = th - init_t;

  if( *l < gr ) {  // startup variables
   auto r = l;
   for( ++r ; ( r != sbs->end() ) && ( *r < gr ) ; )
    ++r;
   Subset nms( std::distance( l , r ) );
   auto nvit = nv.begin();
   auto nmsit = nms.begin();
   // set_startup_costs( Subset ) expects subset entries in time-space
   // [ init_t , f_time_horizon ); active-var indices [ 0 , th - init_t )
   // map to time indices by adding init_t
   while( l != r ) {
    *( nvit++ ) = qf->get_linear_coefficient( *l );
    *( nmsit++ ) = *( l++ ) + init_t;
    }
   set_startup_costs( nv.begin() , std::move( nms ) , true , par , eDryRun );
   if( r == sbs->end() )
    return;
   }

  gl = gr;
  gr = 2 * th - init_t;

  if( *l < gr ) {  // active power variables
   auto r = l;
   for( ++r ; ( r != sbs->end() ) && ( *r < gr ) ; )
    ++r;
   Subset nms( std::distance( l , r ) );
   auto nvit = nv.begin();
   auto nmsit = nms.begin();
   auto nvqit = nvq.begin();
   while( l != r ) {
    *( nvit++ ) = qf->get_linear_coefficient( *l );
    if( with_quad )
     *( nvqit++ ) = qf->get_quadratic_coefficient( *l );
    *( nmsit++ ) = *( l++ ) - gl;
    }
   if( with_quad ) {
    Subset nmsq = nms;  // copy, as it is moved
    set_quad_term( nvq.begin() , std::move( nmsq ) , true , par , eDryRun );
    }
   set_linear_term( nv.begin() , std::move( nms ) , true , par , eDryRun );
   if( r == sbs->end() )
    return;
   }

  gl = gr;
  gr = 3 * th - init_t;

  if( *l < gr ) {  // commitment variables
   auto r = l;
   for( ++r ; ( r != sbs->end() ) && ( *r < gr ) ; )
    ++r;
   Subset nms( std::distance( l , r ) );
   auto nvit = nv.begin();
   auto nmsit = nms.begin();
   while( l != r ) {
    *( nvit++ ) = qf->get_linear_coefficient( *l );
    *( nmsit++ ) = *( l++ ) - gl;
    }
   set_const_term( nv.begin() , std::move( nms ) , true , par , eDryRun );
   if( r == sbs->end() )
    return;
   }

  gl = gr;
  gr = 4 * th - init_t;

  if( *l < gr ) {  // primary spinning reserve variables
   auto r = l;
   for( ++r ; ( r != sbs->end() ) && ( *r < gr ) ; )
    ++r;
   Subset nms( std::distance( l , r ) );
   auto nvit = nv.begin();
   auto nmsit = nms.begin();
   while( l != r ) {
    *( nvit++ ) = qf->get_linear_coefficient( *l );
    *( nmsit++ ) = *( l++ ) - gl;
    }
   set_primary_spinning_reserve_cost( nv.begin() , std::move( nms ) ,
                                      true , par , eDryRun );
   if( r == sbs->end() )
    return;
   }

  gl = gr;
  gr = 5 * th - init_t;

  if( *l < gr ) {  // secondary spinning reserve variables
   auto r = l;
   for( ++r ; ( r != sbs->end() ) && ( *r < gr ) ; )
    ++r;
   Subset nms( std::distance( l , r ) );
   auto nvit = nv.begin();
   auto nmsit = nms.begin();
   while( l != r ) {
    *( nvit++ ) = qf->get_linear_coefficient( *l );
    *( nmsit++ ) = *( l++ ) - gl;
    }
   set_secondary_spinning_reserve_cost( nv.begin() , std::move( nms ) ,
                                        true , par , eDryRun );
   }

  if( l != sbs->end() )
   throw( std::invalid_argument( "ThermalUnitBlock::add_Modification: invalid "
				 "variable in Modification" ) );
  return;

  }  // end( Sbst )

 throw( std::invalid_argument( "ThermalUnitBlock:: unsupported FunctionMod "
			       "from Objective" ) );

 }  // end( ThermalUnitBlock::handle_objective_change )

/*--------------------------------------------------------------------------*/
/*------------------ METHODS OF ThermalUnitBlockSolution -------------------*/
/*--------------------------------------------------------------------------*/

void ThermalUnitBlockSolution::deserialize( const netCDF::NcGroup & group )
{
 // call the method of the base class
 UnitBlockSolution::deserialize( group );

 if( f_number_generators != 1 )
  throw( std::logic_error( "ThermalUnitBlockSolution::deserialize: "
         "thermals have only one generator" ) );

 // deserialize the design - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! ::deserialize< double >( group , f_design , "ThermalDesign" ) )
  f_design = dNaN;

 }  // end( ThermalUnitBlockSolution::deserialize )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlockSolution::read( const Block * block )
{
 auto TUB = dynamic_cast< const ThermalUnitBlock * >( block );
 if( ! TUB )
  throw( std::invalid_argument( "ThermalUnitBlockSolution::read: block "
        "is not a ThermalUnitBlock" ) );

 UnitBlockSolution::read( TUB );  // call the method of the base class

 // read the design- - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 f_design = TUB->get_const_design().get_value();

 }  // end( ThermalUnitBlockSolution::read )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlockSolution::write( Block * block )
{
 UnitBlockSolution::write( block );  // call the method of the base class

 auto TUB = dynamic_cast< ThermalUnitBlock * >( block );
 if( ! TUB )
  throw( std::invalid_argument( "ThermalUnitBlockSolution::write: block "
        "is not a ThermalUnitBlock" ) );

 // write the design - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 TUB->get_design().set_value( f_design );

 // (p, u) have just been restored into the Block by UnitBlockSolution::
 // write(); delegate all the formulation-specific bookkeeping (start_up /
 // shut_down, perspective-cut auxiliaries, ...) to the Block itself. This
 // shares the implementation with the inner DP Solvers and ensures that
 // LagBFunction::get_linearization_constant(), which reads variable
 // values to recompute f(x*) via its saved CostMatrix, sees a state
 // consistent with the saved (p, u).
 TUB->set_solution();

 }  // end( ThermalUnitBlockSolution::write )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlockSolution::serialize( netCDF::NcGroup & group ) const
{
 UnitBlockSolution::serialize( group );  // call the method of the base class

 // serialize the design- - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! std::isnan( f_design ) )
  ::serialize< double >( group , "ThermalDesign" , netCDF::NcDouble() ,
       f_design );

 }  // end( ThermalUnitBlockSolution::serialize )

/*--------------------------------------------------------------------------*/

ThermalUnitBlockSolution * ThermalUnitBlockSolution::scale( double factor )
 const
{
 auto sol = clone();

 if( factor == 1 )
  return( sol );

 guts_of_scale( sol , factor );

 if( ! std::isnan( f_design ) )
  sol->f_design *= factor;

 return( sol );

 }  // end( ThermalUnitBlockSolution::scale )

/*--------------------------------------------------------------------------*/

void ThermalUnitBlockSolution::sum( const Solution * solution ,
           double multiplier )
{
 // call the method of the base class
 UnitBlockSolution::sum( solution , multiplier );

 auto TUBS = dynamic_cast< const ThermalUnitBlockSolution * >( solution );
 if( ! TUBS )
  throw( std::invalid_argument( "ThermalUnitBlockSolution::sum: "
        "solution not a ThermalUnitBlockSolution" ) );

 if( ! std::isnan( f_design ) )
  f_design += TUBS->f_design * multiplier;

 }  // end( ThermalUnitBlockSolution::sum )

/*--------------------------------------------------------------------------*/

ThermalUnitBlockSolution * ThermalUnitBlockSolution::clone( bool empty ) const
{
 auto * sol = new ThermalUnitBlockSolution();

 if( ! empty ) {
  guts_of_clone( sol );
  sol->f_design = f_design;
  }

 return( sol );

 }  // end( ThermalUnitBlockSolution::clone )

/*--------------------------------------------------------------------------*/
/*------------------- End File ThermalUnitBlock.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
