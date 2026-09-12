/*--------------------------------------------------------------------------*/
/*--------------------- File NuclearUnitBlock.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the NuclearUnitBlock class.
 *
 * \author Antonio Frangioni \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Antonio Frangioni, Donato Meoli
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "LinearFunction.h"

#include "NuclearUnitBlock.h"

#include <algorithm>

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

using coeff_pair = LinearFunction::coeff_pair;

using v_coeff_pair = LinearFunction::v_coeff_pair;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register NuclearUnitBlock to the Block factory

SMSpp_insert_in_factory_cpp_1( NuclearUnitBlock );

/*--------------------------------------------------------------------------*/

// register NuclearUnitBlockSolution to the Solution factory

SMSpp_insert_in_factory_cpp_0( NuclearUnitBlockSolution );

/*--------------------------------------------------------------------------*/
/*------------------------------- FUNCTIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

static LinearFunction * LF( Function * f )
{
 return( static_cast< LinearFunction * >( f ) );
 }

/*--------------------------------------------------------------------------*/

template< typename T >
static bool identical( std::vector< T > & vec , const Block::Subset sbst ,
           typename std::vector< T >::const_iterator it )
{
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
        typename std::vector< T >::const_iterator it )
{
 // assign to the sub-vector of vec[] corresponding to the indices in sbst
 // the values found in vector starting at it
 for( auto t : sbst )
  vec[ t ] = *( it++ );
 }

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS OF NuclearUnitBlock ----------------------*/
/*--------------------------------------------------------------------------*/

NuclearUnitBlock::~NuclearUnitBlock()
{
 Constraint::clear( Nuclear_cuts );
 Constraint::clear( DeepDownLink );
 Constraint::clear( DeepLinkConst );
 Constraint::clear( DeepDropConst );
 Constraint::clear( DeepLowConst );
 Constraint::clear( DeepDecreasesPerDayConst );
 Constraint::clear( StartUpsPerDayConst );
 Constraint::clear( ModulationsPerDayConst );
 Constraint::clear( ModulationStartLink );
 Constraint::clear( ModulationMaxLength );
 Constraint::clear( BandMove );
 Constraint::clear( BandKeep );
 Constraint::clear( ModulationEndLink );
 Constraint::clear( BandPower );
 Constraint::clear( BandChoice );
 Constraint::clear( StartUpStability );
 Constraint::clear( ModulationStability );
 Constraint::clear( Modulation_FullRampDown );
 Constraint::clear( Modulation_FullRampUp );
 Constraint::clear( ModulationSameDirection );
 Constraint::clear( DownModulationLink );
 Constraint::clear( ModulationConst );
 Constraint::clear( NoStartUpModulation );
 Constraint::clear( NoDownModulation );
 Constraint::clear( Modulation_RampDown_Constraints );
 Constraint::clear( Modulation_RampUp_Constraints );
 }

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void NuclearUnitBlock::deserialize( const netCDF::NcGroup & group )
{
 // call the method of the base class
 ThermalUnitBlock::deserialize( group );

 // load optional variables ModulationTime and InitModulation or give them
 // default values

 if( ! ::deserialize( group , f_modulation_interval , "ModulationTime" ) )
  f_modulation_interval = 2;

 if( ! ::deserialize( group , f_initial_modulation , "InitModulation" ) )
  f_initial_modulation = f_modulation_interval;

 // load the optional variables ModulationDeltaRampUp and
 // ModulationDeltaRampDown, create the expanded vectors (if needed); they
 // are 0 if not provided

 if( ! ::deserialize( group , "ModulationDeltaRampUp" , f_time_horizon ,
                      v_modulation_ramp_up , true , true ,
                      v_change_intervals ) )
  v_modulation_ramp_up.assign( f_time_horizon , 0 );
 if( ! ::deserialize( group , "ModulationDeltaRampDown" , f_time_horizon ,
                      v_modulation_ramp_down , true , true ,
                      v_change_intervals ) )
  v_modulation_ramp_down.assign( f_time_horizon , 0 );

 // the operating rules: all optional, the defaults giving the original
 // model (single-instant modulations, no daily limit, no cost)

 int tmp;
 f_max_modulation_length = 1;
 if( ::deserialize( group , tmp , "MaxModulationLength" ) )
  f_max_modulation_length = Index( std::max( tmp , 0 ) );

 f_stability_after_start = 0;
 if( ::deserialize( group , tmp , "StabilityAfterStartUp" ) )
  f_stability_after_start = Index( std::max( tmp , 0 ) );

 f_modulations_per_day = -1;
 if( ::deserialize( group , tmp , "ModulationsPerDay" ) )
  f_modulations_per_day = std::max( tmp , 0 );

 f_start_ups_per_day = -1;
 if( ::deserialize( group , tmp , "StartUpsPerDay" ) )
  f_start_ups_per_day = std::max( tmp , 0 );

 f_deep_decreases_per_day = -1;
 if( ::deserialize( group , tmp , "DeepDecreasesPerDay" ) )
  f_deep_decreases_per_day = std::max( tmp , 0 );

 f_day_length = 0;
 if( ::deserialize( group , tmp , "DayLength" ) )
  f_day_length = Index( std::max( tmp , 0 ) );

 // a vector that is not provided, or is all 0 when 0 is the default, is
 // kept empty
 auto opt_vec = [ & ]( const std::string & name , std::vector< double > & v ,
                       bool keep_zero ) {
  if( ::deserialize( group , name , f_time_horizon , v , true , true ,
                     v_change_intervals ) && ( ! keep_zero ) &&
      std::all_of( v.begin() , v.end() ,
                   []( double x ) { return( x == 0 ); } ) )
   v.clear();
  };

 // the two breakpoints of the bands of the output, if any
 v_power_bands.clear();
 if( auto nc = group.getVar( "PowerBands" ) ; ! nc.isNull() ) {
  if( ( nc.getDimCount() != 1 ) || ( nc.getDim( 0 ).getSize() != 2 ) )
   throw( std::invalid_argument( "NuclearUnitBlock::deserialize: PowerBands "
                                 "must have size 2" ) );
  v_power_bands.resize( 2 );
  nc.getVar( { 0 } , { 2 } , v_power_bands.data() );
  if( v_power_bands[ 0 ] >= v_power_bands[ 1 ] )
   throw( std::invalid_argument( "NuclearUnitBlock::deserialize: the two "
                                 "PowerBands are not increasing" ) );
  }

 opt_vec( "DownModulationCost" , v_down_modulation_cost , false );
 opt_vec( "DeepDecreaseThreshold" , v_deep_threshold , true );
 opt_vec( "DeepDecreaseGradient" , v_deep_gradient , true );
 opt_vec( "DeepDecreaseCost" , v_deep_cost , false );
 if( v_deep_threshold.empty() || v_deep_gradient.empty() ) {
  v_deep_threshold.clear();
  v_deep_gradient.clear();
  v_deep_cost.clear();
  }

 check_data_consistency();

 }  // end( NuclearUnitBlock::deserialize )

/*--------------------------------------------------------------------------*/

#ifndef NDEBUG

/*
std::vector< std::string > NuclearUnitBlock::expected_dims( void )
 const {
 static const std::vector< std::string > ed = { };

 auto ret = UnitBlock::expected_dims();
 ret.insert( ret.end() , ed.begin() , ed.end() );

 return( ret );
 }

----------------------------------------------------------------------------*/

std::vector< std::string > NuclearUnitBlock::expected_vars( void )
 const {
 static const std::vector< std::string > ev =
 { "ModulationTime" , "InitModulation" , "ModulationDeltaRampUp" ,
   "ModulationDeltaRampDown" , "MaxModulationLength" ,
   "StabilityAfterStartUp" , "PowerBands" , "ModulationsPerDay" ,
   "StartUpsPerDay" , "DeepDecreasesPerDay" , "DayLength" ,
   "DownModulationCost" , "DeepDecreaseThreshold" , "DeepDecreaseGradient" ,
   "DeepDecreaseCost"
   };

 auto ret = ThermalUnitBlock::expected_vars();
 ret.insert( ret.end() , ev.begin() , ev.end() );

 return( ret );
 }

#endif

/*--------------------------------------------------------------------------*/

void NuclearUnitBlock::check_data_consistency( void ) const
{
 static const std::string fn = "NuclearUnitBlock::check_data_consistency";

 // DeltaRampUp/Down - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // both ramp vectors inherited from ThermalUnitBlock are mandatory here
 if( v_DeltaRampUp.empty() )
  throw( std::invalid_argument( fn + ": DeltaRampUp not present" ) );

 if( v_DeltaRampDown.empty() )
  throw( std::invalid_argument( fn + ": DeltaRampDown not present" ) );

 // ModulationTime and InitModulation - - - - - - - - - - - - - - - - - - - -
 if( f_modulation_interval < 2 )
  throw( std::invalid_argument( fn + ": ModulationTime is " +
                                std::to_string( f_modulation_interval ) +
                                ", but it must be at least 2" ) );

 if( f_initial_modulation < 1 )
  throw( std::invalid_argument( fn + ": InitModulation is " +
                                std::to_string( f_initial_modulation ) +
                                ", but it must be at least 1" ) );

 // ModulationDeltaRampUp/Down - - - - - - - - - - - - - - - - - - - - - - -
 // for each t: 0 \leq v_modulation_ramp_up[t]   \leq v_DeltaRampUp[t]
 //             0 \leq v_modulation_ramp_down[t] \leq v_DeltaRampUp[t]
 assert( v_modulation_ramp_up.size() == f_time_horizon );
 assert( v_modulation_ramp_down.size() == f_time_horizon );

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {
  if( v_modulation_ramp_up[ t ] < 0 )
   throw( std::logic_error( fn + ": modulation ramp up at time " +
                            std::to_string( t ) + " is " +
                            std::to_string( v_modulation_ramp_up[ t ] ) +
                            " < 0" ) );

  if( v_modulation_ramp_up[ t ] > v_DeltaRampUp[ t ] )
   throw( std::logic_error( fn + ": modulation ramp up at time " +
                            std::to_string( t ) + " is " +
                            std::to_string( v_modulation_ramp_up[ t ] ) +
                            " > ramp up = " +
                            std::to_string( v_DeltaRampUp[ t ] ) ) );

  if( v_modulation_ramp_down[ t ] < 0 )
   throw( std::logic_error( fn + ": modulation ramp down at time " +
                            std::to_string( t ) + " is " +
                            std::to_string( v_modulation_ramp_down[ t ] ) +
                            " < 0" ) );

  if( v_modulation_ramp_down[ t ] > v_DeltaRampDown[ t ] )
   throw( std::logic_error( fn + ": modulation ramp down at time " +
                            std::to_string( t ) + " is " +
                            std::to_string( v_modulation_ramp_down[ t ] ) +
                            " > ramp down = " +
                            std::to_string( v_DeltaRampDown[ t ] ) ) );

  }  // end( for t )

 // the operating rules- - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( f_max_modulation_length < 1 )
  throw( std::invalid_argument( fn + ": MaxModulationLength must be >= 1" ) );

 for( Index t = 0 ; t < v_down_modulation_cost.size() ; ++t )
  if( v_down_modulation_cost[ t ] < 0 )
   throw( std::invalid_argument( fn + ": DownModulationCost at time " +
                                 std::to_string( t ) + " is < 0" ) );

 for( Index t = 0 ; t < v_deep_gradient.size() ; ++t )
  if( v_deep_gradient[ t ] <= 0 )
   throw( std::invalid_argument( fn + ": DeepDecreaseGradient at time " +
                                 std::to_string( t ) + " is <= 0" ) );

 for( Index t = 0 ; t < v_deep_cost.size() ; ++t )
  if( v_deep_cost[ t ] < 0 )
   throw( std::invalid_argument( fn + ": DeepDecreaseCost at time " +
                                 std::to_string( t ) + " is < 0" ) );
 }  // end( NuclearUnitBlock::check_data_consistency )

/*--------------------------------------------------------------------------*/

void NuclearUnitBlock::generate_abstract_variables( Configuration * stvv ) {

 if( variables_generated() )
  return; // variables have already been generated

 // defined here inside:
 // - init_t = first instant in which commitment is free
 // - v_commitment[ f_time_horizon ]
 // - v_active_power[ f_time_horizon ]
 // - v_start_up[ f_time_horizon - init_t ]
 // - v_shut_down[ f_time_horizon - init_t ]
 ThermalUnitBlock::generate_abstract_variables( stvv );

 // the formulation of the operating rules, from the same Configuration that
 // selects that of the ThermalUnitBlock [see TightRules]
 if( ( ! stvv ) && f_BlockConfig )
  stvv = f_BlockConfig->f_static_variables_Configuration;
 if( auto sci = dynamic_cast< SimpleConfiguration< int > * >( stvv ) ) {
  f_tight_rules = sci->f_value & TightRules;
  f_tight_ramp = sci->f_value & TightRamp;
  f_tight_cuts = f_tight_rules && ( sci->f_value & TightCuts );
  }

 // Modulation Variable- - - - - - - - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // resize v_modulation to f_time_horizon, thereby creating the ColVariable
 v_modulation.resize( f_time_horizon );

 // set the type to binary
 for( auto & var : v_modulation )
  var.set_type( ColVariable::kBinary );

 // add the corresponding group to ThermalUnitBlock static Variable
 add_static_variable( v_modulation , "m_thermal" );

 // fixing the modulation variable to 0 - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // for all time instants between 0 and f_modulation_interval -
 // f_initial_modulation (right extreme excluded)

 const int lock0 = std::min( f_modulation_interval - f_initial_modulation ,
                             int( f_time_horizon ) );
 for( int t = 0 ; t < lock0 ; ) {
  v_modulation[ t ].set_value( 0.0 );
  v_modulation[ t++ ].is_fixed( true , eNoMod );
  }

 // and then, if the unit is off at time 0 (f_InitUpDownTime <= 0) then all
 // the u_t for t = 0, ..., init_t - 1 are fixed to 0 as well, which means
 // that the m_t must be fixed to 0 due to the constraint m_t leq u_y
 if( f_InitUpDownTime <= 0 )
  for( Index t = 0 ; t <  init_t ; ) {
   v_modulation[ t ].set_value( 0.0 );
   v_modulation[ t++ ].is_fixed( true , eNoMod );
   }

 // the downward modulation Variable, if the direction matters: fixed to 0
 // wherever the modulation is
 if( has_modulation_direction() ) {
  v_modulation_down.resize( f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   v_modulation_down[ t ].set_type( ColVariable::kBinary );
   if( v_modulation[ t ].is_fixed() ) {
    v_modulation_down[ t ].set_value( 0.0 );
    v_modulation_down[ t ].is_fixed( true , eNoMod );
    }
   }
  add_static_variable( v_modulation_down , "m_down_nuclear" );
  }

 // the modulation start Variable, if a modulation may last more than one
 // instant and either the modulations per day are limited or the tight rows
 // are used, which are written on it
 if( ( f_max_modulation_length > 1 ) &&
     ( f_tight_rules || f_tight_cuts || ( f_modulations_per_day >= 0 ) ) ) {
  v_modulation_start.resize( f_time_horizon );
  for( auto & var : v_modulation_start )
   var.set_type( ColVariable::kPosUnitary );
  add_static_variable( v_modulation_start , "m_start_nuclear" );
  }

 // the band of the output and the end of a modulation, if the output is
 // banded: three binaries per instant, all zero when the unit is off, and
 // one [ 0 , 1 ] variable per instant that is 1 exactly at the last step of
 // a modulation
 if( has_power_bands() ) {
  v_band.resize( 3 * f_time_horizon );
  for( auto & var : v_band )
   var.set_type( ColVariable::kBinary );
  add_static_variable( v_band , "band_nuclear" );
  v_modulation_end.resize( f_time_horizon );
  for( auto & var : v_modulation_end )
   var.set_type( ColVariable::kPosUnitary );
  add_static_variable( v_modulation_end , "m_end_nuclear" );
  }

 // the deep-decrease Variable: at t = 0 there is no deep decrease unless
 // the unit is on at the beginning
 if( has_deep_decrease() ) {
  v_deep.resize( f_time_horizon );
  v_deep_drop.resize( f_time_horizon );
  v_deep_low.resize( f_time_horizon );
  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   v_deep[ t ].set_type( ColVariable::kBinary );
   v_deep_drop[ t ].set_type( ColVariable::kBinary );
   v_deep_low[ t ].set_type( ColVariable::kBinary );
   }
  if( f_InitUpDownTime <= 0 )
   for( auto v : { & v_deep , & v_deep_drop , & v_deep_low } ) {
    ( *v )[ 0 ].set_value( 0.0 );
    ( *v )[ 0 ].is_fixed( true , eNoMod );
    }
  add_static_variable( v_deep , "deep_nuclear" );
  add_static_variable( v_deep_drop , "deep_drop_nuclear" );
  add_static_variable( v_deep_low , "deep_low_nuclear" );
  }
 } // end( NuclearUnitBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void NuclearUnitBlock::generate_abstract_constraints( Configuration * stcc )
{
 if( constraints_generated() )
  return; // constraints have already been generated

 // call the method of the base class - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // note that set_constraints_generated() is called there inside, so that it
 // does not need to be done again
 ThermalUnitBlock::generate_abstract_constraints( stcc );

 // important information from the base class:
 // - if f_InitUpDownTime > 0 then the unit was on before the initial time
 //   instant 0, i.e.,  u_{0 - 1} = 1, otherwise it was off, i.e.,
 //   u_{0 - 1} = 1
 // - if u_{0 - 1} = 1, then f_InitialPower = p_{0 - 1}
 // - v_StartUpLimit, the maximum power on startup
 // - v_ShutDownLimit, the maximum power on shutdown
 // - v_DeltaRampUp, the ramp-up delta
 // - v_DeltaRampDown, the ramp-down delta

 // construct the modulation ramp-up constraint - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // p_t - p_{t-1} - \Delta^M_{t+} u_{t-1} -
 // ( \Delta_{t+} - \Delta^M_{t+} ) m_t - \bar{l}_t v_t \leq 0

 // when the direction matters an upward step cannot decrease the output
 // nor a downward one increase it, i.e., the row gets + \Delta_{t+} d_t
 const bool dir = has_modulation_direction();

 Modulation_RampUp_Constraints.resize( f_time_horizon );

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {
  Index np = t ? 5 : 3;
  if( t < init_t )
   --np;
  if( dir )
   ++np;
  LinearFunction::v_coeff_pair cf( np );
  double RHS = 0;
  auto cfit = cf.begin();

  *( cfit++ ) = coeff_pair( & v_active_power[ t ] , 1.0 );
  *( cfit++ ) = coeff_pair( & v_modulation[ t ] ,
        - ( v_DeltaRampUp[ t ] - v_modulation_ramp_up[ t ] )
      );
  // the two terms "- p_{t-1}" and "- \Delta^M_{t+} u_{t-1}" only exist if
  // t > 0, as otherwise p_{t-1} and u_{t-1} are undefined
  if( t ) {
   *( cfit++ ) = coeff_pair( & v_commitment[ t - 1 ] ,
       - v_modulation_ramp_up[ t ] );

   *( cfit++ ) = coeff_pair( & v_active_power[ t - 1 ] , -1.0 );
   }
  else {
   // if t == 0, the "- p_{t-1}" term is fixed and equal to - f_InitialPower,
   // so there is no explicit term in the constraint (since the variable does
   // not exist) and the RHS becomes f_InitialPower
   RHS = f_InitialPower;
   // similarly, the "- \Delta^M_{t+} u_{t-1}" term is fixed, and it is
   // equal to - v_modulation_ramp_up[ t ] if u_{t-1} = 1 (i.e.,
   // f_InitUpDownTime > 0) and 0 otherwise, so this has to be added to RHS
   // (changing the sign)
   if( f_InitUpDownTime > 0 )
    RHS += v_modulation_ramp_up[ 0 ];
   }

  // the term - \bar{l}_t v_t only exist if t >= init_t, as for t < init_t
  // the commitment status if fixed and start-ups are not allowed, hence
  // the corresponding start-up variables are not even defined
  if( t >= init_t )
   *( cfit++ ) = coeff_pair( & v_start_up[ t - init_t ] ,
                             - v_StartUpLimit[ t ] );

  if( dir )
   *cfit = coeff_pair( & v_modulation_down[ t ] , v_DeltaRampUp[ t ] );

  Modulation_RampUp_Constraints[ t ].set_lhs( - Inf< double >() );
  Modulation_RampUp_Constraints[ t ].set_rhs( RHS );
  Modulation_RampUp_Constraints[ t ].set_function(
            new LinearFunction( std::move( cf ) ) );
  }

 add_static_constraint( Modulation_RampUp_Constraints ,
      "Modulation_RampUp_Constraints_Nuclear" );

 // construct the modulation ramp-down constraint - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // p_{t-1} - p_t - \Delta^M_{t-} u_t -
 // ( \Delta_{t-} - \Delta^M_{t-} ) m_t - \bar{u}_t w_t \leq 0

 // when the direction matters the row is
 // p_{t-1} - p_t - \Delta^M_{t-} u_t + \Delta^M_{t-} m_t - \Delta_{t-} d_t
 // - \bar{u}_t w_t \leq 0
 Modulation_RampDown_Constraints.resize( f_time_horizon );

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {
  Index np = t ? 5 : 4;
  if( t < init_t )
   --np;
  if( dir )
   ++np;
  LinearFunction::v_coeff_pair cf( np );
  auto cfit = cf.begin();

  *( cfit++ ) = coeff_pair( & v_active_power[ t ] , -1.0 );
  *( cfit++ ) = coeff_pair( & v_commitment[ t ] ,
      - v_modulation_ramp_down[ t ] );
  *( cfit++ ) = coeff_pair( & v_modulation[ t ] , dir ?
      v_modulation_ramp_down[ t ] :
      - ( v_DeltaRampDown[ t ] - v_modulation_ramp_down[ t ] ) );
  if( dir )
   *( cfit++ ) = coeff_pair( & v_modulation_down[ t ] ,
                             - v_DeltaRampDown[ t ] );

  // the terms "p_{t-1}" only exists if t > 0, as otherwise p_{t-1} is
  // undefined
  if( t )
   *( cfit++ ) = coeff_pair( & v_active_power[ t - 1 ] , 1.0 );

  // the term - \bar{l}_t v_t only exist if t >= init_t, as for t < init_t
  // the commitment status if fixed and shut-downs are not allowed, hence
  // the corresponding shut-down variables are not even defined
  if( t >= init_t )
   *cfit = coeff_pair( & v_shut_down[ t - init_t ] , - v_ShutDownLimit[ t ] );

  Modulation_RampDown_Constraints[ t ].set_lhs( - Inf< double >() );
  // if t == 0, the "p_{t-1}" term is fixed and equal to f_InitialPower, so
  // there is no explicit term in the constraint (since the variable does
  // not exist) and the RHS becomes - f_InitialPower
  Modulation_RampDown_Constraints[ t ].set_rhs( t ? 0 : - f_InitialPower );
  Modulation_RampDown_Constraints[ t ].set_function(
            new LinearFunction( std::move( cf ) ) );
  }

 add_static_constraint( Modulation_RampDown_Constraints ,
      "Modulation_RampDown_Constraints_Nuclear" );

 // construct the logical constraints - - - - - - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // m_t - u_t \leq 0  (modulation ==> unit up)
 // note: these only have to be constructed for t >= init_t, as for
 // t < init_t either u_t is fixed to 1, and the constraint is redundant, or
 // u_t is fixed to 0 and m_t has been fixed in generate_abstract_variables()

 NoDownModulation.resize( f_time_horizon - init_t );

 for( Index t = init_t ; t < f_time_horizon ; ++t ) {
  LinearFunction::v_coeff_pair cf( 2 );

  cf[ 0 ] = coeff_pair( & v_modulation[ t ] , 1.0 );
  cf[ 1 ] = coeff_pair( & v_commitment[ t ] , -1.0 );

  NoDownModulation[ t - init_t ].set_lhs( - Inf< double >() );
  NoDownModulation[ t - init_t ].set_rhs( 0 );
  NoDownModulation[ t - init_t ].set_function(
            new LinearFunction( std::move( cf ) ) );
  }

 add_static_constraint( NoDownModulation , "NoDownModulation_Nuclear" );

 // construct the logical constraints m_t + v_t \leq 1  - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // the unit is not modulating while starting up
 // note: these only have to be constructed for t >= init_t, as for
 // t < init_t u_t is fixed (no matter if to 0 or 1) and therefore no
 // start-up can ever occur; in fact, the start-up variables are not even
 // defined for t < init_t. it may also be that the m_t are fixed for those
 // t: this happens if u_t is fixed to 0, but not if u_t is fixed to 1, in
 // which case modulations can occur within the first init_t periods unless
 // forbidden by the initial state (f_initial_modulation), but the latter
 // case is already taken care of in generate_abstract_variables()

 NoStartUpModulation.resize( f_time_horizon - init_t );

 for( Index t = init_t ; t < f_time_horizon ; ++t ) {
  LinearFunction::v_coeff_pair cf( 2 );

  cf[ 0 ] = coeff_pair( & v_modulation[ t ] , 1.0 );
  cf[ 1 ] = coeff_pair( & v_start_up[ t - init_t ] , -1.0 );

  NoStartUpModulation[ t - init_t ].set_lhs( - Inf< double >() );
  NoStartUpModulation[ t - init_t ].set_rhs( 1.0 );
  NoStartUpModulation[ t - init_t ].set_function(
            new LinearFunction( std::move( cf ) ) );
  }

 add_static_constraint( NoStartUpModulation ,
      "NoStartUpModulation_Nuclear" );

 // construct the modulation constraint proper- - - - - - - - - - - - - - - -
 // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // sum_{h = \max\{ 0 , t - \tau^M + 1 \}}^t m_h \leq 1
 // recall that \tau^M >= 2: thus, for t = 0 one has t - \tau^M + 1 < 0 and
 // the sum would go for h = 0 to 0, i.e., it would be m[ 0 ]; but
 // m[ 0 ] <= 1, hence the first constraint is also redundant
 // more in general: for t < f_modulation_interval - f_initial_modulation
 // all m_t are fixed to 0, hence the constraint is useless until
 // t >= f_modulation_interval - f_initial_modulation + 1
 // similarly, if the unit is off at time 0 (f_InitUpDownTime <= 0) then all
 // the u_t for t = 0, ..., init_t - 1 are fixed to 0 as well, which means
 // that the m_t must be fixed to 0 due to the constraint m_t leq u_t;
 // hence the constraint is useless until t >= init_t + 1
 // (the window is the stability rule of single-instant modulations; for
 // longer ones it is replaced by the rows of generate_operating_rules())
 Index first_c = std::max( f_modulation_interval - f_initial_modulation ,
         int( 0 ) );
 if( f_InitUpDownTime <= 0 )
  first_c = std::max( first_c , init_t );
 ++first_c;
 if( ( f_max_modulation_length > 1 ) || ( first_c > f_time_horizon ) )
  first_c = f_time_horizon;

 ModulationConst.resize( f_time_horizon - first_c );

 for( Index t = first_c ; t < f_time_horizon ; ++t ) {
  Index h = std::max( int( 0 ) , int( t ) - f_modulation_interval + 1 );
  LinearFunction::v_coeff_pair cf( t - h + 1 );

  for( auto cfit = cf.begin() ; h <= t ; )
   *( cfit++ ) = coeff_pair( & v_modulation[ h++ ] , 1.0 );

  ModulationConst[ t - first_c ].set_lhs( - Inf< double >() );
  ModulationConst[ t - first_c ].set_rhs( 1.0 );
  ModulationConst[ t - first_c ].set_function(
            new LinearFunction( std::move( cf ) ) );
  }

 add_static_constraint( ModulationConst , "ModulationConst_Nuclear" );

 generate_operating_rules();

 } // end( NuclearUnitBlock::generate_abstract_constraints )
/*--------------------------------------------------------------------------*/

void NuclearUnitBlock::generate_operating_rules( void )
{
 const Index T = f_time_horizon;
 const bool dir = has_modulation_direction();
 const Index L = f_max_modulation_length;
 const double INF = Inf< double >();

 // append the row lhs <= sum coeff * var <= rhs
 auto row = [ & ]( std::vector< FRowConstraint > & rows ,
                   LinearFunction::v_coeff_pair && cf , double lhs ,
                   double rhs ) {
  rows.emplace_back();
  rows.back().set_lhs( lhs );
  rows.back().set_rhs( rhs );
  rows.back().set_function( new LinearFunction( std::move( cf ) ) );
  };

 // the days, as [ first , past-the-end ) intervals of instants
 std::vector< std::pair< Index , Index > > days;
 for( Index d0 = 0 ; d0 < T ; ) {
  const Index d1 = f_day_length ? std::min( d0 + f_day_length , T ) : T;
  days.emplace_back( d0 , d1 );
  d0 = d1;
  }

 // d_t <= m_t - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( dir ) {
  DownModulationLink.reserve( T );
  for( Index t = 0 ; t < T ; ++t )
   row( DownModulationLink , { coeff_pair( & v_modulation_down[ t ] , 1.0 ) ,
                               coeff_pair( & v_modulation[ t ] , -1.0 ) } ,
        -INF , 0 );
  add_static_constraint( DownModulationLink , "DownModulationLink_Nuclear" );
  }

 if( L > 1 ) {
  // same direction along a modulation- - - - - - - - - - - - - - - - - - -
  // d_{t+1} - d_t + m_t <= 1 , d_t - d_{t+1} + m_{t+1} <= 1
  ModulationSameDirection.reserve( 2 * T );
  for( Index t = 0 ; t + 1 < T ; ++t ) {
   row( ModulationSameDirection ,
        { coeff_pair( & v_modulation_down[ t + 1 ] , 1.0 ) ,
          coeff_pair( & v_modulation_down[ t ] , -1.0 ) ,
          coeff_pair( & v_modulation[ t ] , 1.0 ) } , -INF , 1 );
   row( ModulationSameDirection ,
        { coeff_pair( & v_modulation_down[ t ] , 1.0 ) ,
          coeff_pair( & v_modulation_down[ t + 1 ] , -1.0 ) ,
          coeff_pair( & v_modulation[ t + 1 ] , 1.0 ) } , -INF , 1 );
   }
  add_static_constraint( ModulationSameDirection ,
                         "ModulationSameDirection_Nuclear" );

  // full ramp but at the last step- - - - - - - - - - - - - - - - - - - - -
  // p_t - p_{t-1} - D+_t m_{t+1} - M_t m_t + M_t d_t >= - M_t
  // p_{t-1} - p_t - D-_t m_{t+1} - M'_t d_t >= - M'_t
  // with M_t = D+_t + max{ D-_t , SD_t } and M'_t = D-_t + max{ D+_t , SU_t },
  // since the modulation ramp constraints bound the decrease of the output
  // by max{ D-_t , SD_t } and its increase by max{ D+_t , SU_t }; at t = 0
  // p_{-1} is the (fixed) initial power
  Modulation_FullRampUp.reserve( T );
  Modulation_FullRampDown.reserve( T );
  for( Index t = 0 ; t < T ; ++t ) {
   const double M = full_ramp_up_M( t );
   const double Md = full_ramp_down_M( t );

   LinearFunction::v_coeff_pair cu , cd;
   cu.push_back( coeff_pair( & v_active_power[ t ] , 1.0 ) );
   cd.push_back( coeff_pair( & v_active_power[ t ] , -1.0 ) );
   if( t ) {
    cu.push_back( coeff_pair( & v_active_power[ t - 1 ] , -1.0 ) );
    cd.push_back( coeff_pair( & v_active_power[ t - 1 ] , 1.0 ) );
    }
   if( t + 1 < T ) {
    cu.push_back( coeff_pair( & v_modulation[ t + 1 ] ,
                              - v_DeltaRampUp[ t ] ) );
    cd.push_back( coeff_pair( & v_modulation[ t + 1 ] ,
                              - v_DeltaRampDown[ t ] ) );
    }
   if( f_tight_ramp ) {
    // the same rows with one coefficient per case rather than one for all:
    // a stable instant moves down by at most the stability ramp, a downward
    // step by the full ramp, and only a shut-down brings the output to 0
    // from as high as the shut-down limit, so that
    // A_t = D+_t + MD-_t , B_t = D+_t + D-_t ,
    // A'_t = D-_t + MD+_t , B'_t = ( D+_t - MD+_t )^+
    const double A = full_ramp_up_const( t );
    const double B = v_DeltaRampUp[ t ] + v_DeltaRampDown[ t ];
    const double Ad = full_ramp_down_const( t );
    const double Bd = std::max( v_DeltaRampUp[ t ] -
                                v_modulation_ramp_up[ t ] , 0.0 );
    cu.push_back( coeff_pair( & v_modulation[ t ] , - A ) );
    cu.push_back( coeff_pair( & v_modulation_down[ t ] , B ) );
    cd.push_back( coeff_pair( & v_modulation_down[ t ] , - ( Ad + Bd ) ) );
    cd.push_back( coeff_pair( & v_modulation[ t ] , Bd ) );
    if( t >= init_t ) {
     cu.push_back( coeff_pair( & v_shut_down[ t - init_t ] ,
                               v_ShutDownLimit[ t ] ) );
     cd.push_back( coeff_pair( & v_start_up[ t - init_t ] ,
                               v_StartUpLimit[ t ] ) );
     }
    row( Modulation_FullRampUp , std::move( cu ) ,
         - A + ( t ? 0 : f_InitialPower ) , INF );
    row( Modulation_FullRampDown , std::move( cd ) ,
         - Ad - ( t ? 0 : f_InitialPower ) , INF );
    }
   else {
    cu.push_back( coeff_pair( & v_modulation[ t ] , - M ) );
    cu.push_back( coeff_pair( & v_modulation_down[ t ] , M ) );
    cd.push_back( coeff_pair( & v_modulation_down[ t ] , - Md ) );
    row( Modulation_FullRampUp , std::move( cu ) ,
         - M + ( t ? 0 : f_InitialPower ) , INF );
    row( Modulation_FullRampDown , std::move( cd ) ,
         - Md - ( t ? 0 : f_InitialPower ) , INF );
    }
   }
  add_static_constraint( Modulation_FullRampUp ,
                         "Modulation_FullRampUp_Nuclear" );
  add_static_constraint( Modulation_FullRampDown ,
                         "Modulation_FullRampDown_Nuclear" );

  // s_t >= m_t - m_{t-1}, the start of a modulation, whenever it is there
  if( ! v_modulation_start.empty() ) {
   ModulationStartLink.reserve( T );
   for( Index t = 0 ; t < T ; ++t ) {
    LinearFunction::v_coeff_pair cf;
    cf.push_back( coeff_pair( & v_modulation_start[ t ] , 1.0 ) );
    cf.push_back( coeff_pair( & v_modulation[ t ] , -1.0 ) );
    if( t )
     cf.push_back( coeff_pair( & v_modulation[ t - 1 ] , 1.0 ) );
    row( ModulationStartLink , std::move( cf ) , 0 , INF );
    }
   add_static_constraint( ModulationStartLink ,
                          "ModulationStartLink_Nuclear" );
   }

  // stability: with the tight rows, two modulations cannot start less than
  // tau^M instants apart, i.e., sum_{h=t}^{t+tau^M-1} s_h <= 1, which is a
  // clique of the conflict graph of the starts; otherwise a modulation
  // ending at t forbids the steps up to t + tau^M - 1, i.e.,
  // sum_{h=2}^{K} m_{t+h} + ( K - 1 ) ( m_t - m_{t+1} ) <= K - 1 with
  // K = min{ tau^M - 1 , T - 1 - t }, whose relaxation is much weaker, the
  // end of the modulation being spread over m_t - m_{t+1}
  ModulationStability.reserve( 3 * T );
  for( Index t = 0 ; t + 2 < T ; ++t ) {
   const Index K = std::min( Index( f_modulation_interval - 1 ) , T - 1 - t );
   if( K < 2 )
    continue;
   LinearFunction::v_coeff_pair cf;
   for( Index h = 2 ; h <= K ; ++h )
    cf.push_back( coeff_pair( & v_modulation[ t + h ] , 1.0 ) );
   cf.push_back( coeff_pair( & v_modulation[ t ] , double( K - 1 ) ) );
   cf.push_back( coeff_pair( & v_modulation[ t + 1 ] , - double( K - 1 ) ) );
   row( ModulationStability , std::move( cf ) , -INF , double( K - 1 ) );
   }
  if( f_tight_rules && ( ! f_tight_cuts ) )
   for( Index t = 0 ; t + 1 < T ; ++t ) {
    // at most one start in any tau^M instants, two of them being at least
    // that far apart
    const Index K = std::min( Index( f_modulation_interval ) , T - t );
    if( K >= 2 ) {
     LinearFunction::v_coeff_pair cf;
     for( Index h = 0 ; h < K ; ++h )
      cf.push_back( coeff_pair( & v_modulation_start[ t + h ] , 1.0 ) );
     row( ModulationStability , std::move( cf ) , -INF , 1.0 );
     }
    // a modulation ending at t forbids the starts up to t + tau^M - 1:
    // m_t - m_{t+1} + sum_{h=t+1}^{t+tau^M-1} s_h <= 1
    const Index K1 = std::min( Index( f_modulation_interval - 1 ) ,
                               T - 1 - t );
    if( K1 >= 1 ) {
     LinearFunction::v_coeff_pair cf;
     cf.push_back( coeff_pair( & v_modulation[ t ] , 1.0 ) );
     cf.push_back( coeff_pair( & v_modulation[ t + 1 ] , -1.0 ) );
     for( Index h = 1 ; h <= K1 ; ++h )
      cf.push_back( coeff_pair( & v_modulation_start[ t + h ] , 1.0 ) );
     row( ModulationStability , std::move( cf ) , -INF , 1.0 );
     }
    }
  add_static_constraint( ModulationStability , "ModulationStability_Nuclear" );

  // maximum length: sum_{h=t}^{t+L} m_h <= L, and, with the tight rows, each
  // step belongs to a modulation started in the last L instants, i.e.,
  // m_t <= sum_{h=(t-L+1)^+}^{t} s_h
  ModulationMaxLength.reserve( 2 * T );
  for( Index t = 0 ; t + L < T ; ++t ) {
   LinearFunction::v_coeff_pair cf;
   for( Index h = t ; h <= t + L ; ++h )
    cf.push_back( coeff_pair( & v_modulation[ h ] , 1.0 ) );
   row( ModulationMaxLength , std::move( cf ) , -INF , double( L ) );
   }
  if( f_tight_rules && ( ! f_tight_cuts ) )
   for( Index t = 0 ; t < T ; ++t ) {
    LinearFunction::v_coeff_pair cf;
    cf.push_back( coeff_pair( & v_modulation[ t ] , 1.0 ) );
    for( Index h = ( t >= L ? t - L + 1 : 0 ) ; h <= t ; ++h )
     cf.push_back( coeff_pair( & v_modulation_start[ h ] , -1.0 ) );
    row( ModulationMaxLength , std::move( cf ) , -INF , 0.0 );
    }
  add_static_constraint( ModulationMaxLength , "ModulationMaxLength_Nuclear" );
  }

 // the bands of the output - - - - - - - - - - - - - - - - - - - - - - - -
 // b^1_t + b^2_t + b^3_t = u_t : one band per on instant, none when off;
 // the output is in its band, save at the instants in which a modulation is
 // in progress and does not end, where it travels between two of them:
 //   p_t >= Pmin b^1 + B_1 b^2 + B_2 b^3 - ( Pmax - Pmin )( m_t - e_t )
 //   p_t <= B_1 b^1 + B_2 b^2 + Pmax b^3 + ( Pmax - Pmin )( m_t - e_t )
 // with e_t = m_t ( 1 - m_{t+1} ) the last step of a modulation; the band
 // only changes there, and when it does it moves to an adjacent one:
 //   b^k_t - b^k_{t-1} <= e_t + ( 1 - u_{t-1} ) , and the other way round
 //   b^k_t + b^k_{t-1} <= 2 - e_t , b^1_t + b^3_{t-1} <= 1 , and vice versa
 if( has_power_bands() ) {
  const double B1 = v_power_bands[ 0 ] , B2 = v_power_bands[ 1 ];
  auto band = [ & ]( Index k , Index t ) { return( & v_band[ k * T + t ] ); };

  BandChoice.reserve( T );
  BandPower.reserve( 2 * T );
  ModulationEndLink.reserve( 3 * T );
  BandKeep.reserve( 6 * T );
  BandMove.reserve( 5 * T );

  for( Index t = 0 ; t < T ; ++t ) {
   const double pmin = get_min_power( t ) , pmax = get_max_power( t );
   const double M = pmax - pmin;

   row( BandChoice , { coeff_pair( band( 0 , t ) , 1.0 ) ,
                       coeff_pair( band( 1 , t ) , 1.0 ) ,
                       coeff_pair( band( 2 , t ) , 1.0 ) ,
                       coeff_pair( & v_commitment[ t ] , -1.0 ) } , 0 , 0 );

   row( BandPower , { coeff_pair( & v_active_power[ t ] , 1.0 ) ,
                      coeff_pair( band( 0 , t ) , - pmin ) ,
                      coeff_pair( band( 1 , t ) , - B1 ) ,
                      coeff_pair( band( 2 , t ) , - B2 ) ,
                      coeff_pair( & v_modulation[ t ] , M ) ,
                      coeff_pair( & v_modulation_end[ t ] , - M ) } , 0 ,
        INF );
   row( BandPower , { coeff_pair( & v_active_power[ t ] , 1.0 ) ,
                      coeff_pair( band( 0 , t ) , - B1 ) ,
                      coeff_pair( band( 1 , t ) , - B2 ) ,
                      coeff_pair( band( 2 , t ) , - pmax ) ,
                      coeff_pair( & v_modulation[ t ] , - M ) ,
                      coeff_pair( & v_modulation_end[ t ] , M ) } , -INF , 0 );

   // e_t = m_t ( 1 - m_{t+1} ), the last step of a modulation
   row( ModulationEndLink , { coeff_pair( & v_modulation_end[ t ] , 1.0 ) ,
                              coeff_pair( & v_modulation[ t ] , -1.0 ) } ,
        -INF , 0 );
   if( t + 1 < T ) {
    row( ModulationEndLink ,
         { coeff_pair( & v_modulation_end[ t ] , 1.0 ) ,
           coeff_pair( & v_modulation[ t + 1 ] , 1.0 ) } , -INF , 1 );
    row( ModulationEndLink ,
         { coeff_pair( & v_modulation_end[ t ] , 1.0 ) ,
           coeff_pair( & v_modulation[ t ] , -1.0 ) ,
           coeff_pair( & v_modulation[ t + 1 ] , 1.0 ) } , 0 , INF );
    }
   // at the last instant of the horizon a modulation may be still in
   // progress, exactly as it may in the DP: e_{T-1} is then free below
   // m_{T-1}, and the band changes there or does not. A modulation that
   // the horizon cuts, however, still owes its last step, hence it may
   // only have L^M - 1 of them:
   // sum_{h=T-L^M}^{T-1} m_h - e_{T-1} <= L^M - 1
   else if( L > 1 ) {
    LinearFunction::v_coeff_pair cf;
    for( Index h = ( T >= L ? T - L : 0 ) ; h < T ; ++h )
     cf.push_back( coeff_pair( & v_modulation[ h ] , 1.0 ) );
    cf.push_back( coeff_pair( & v_modulation_end[ T - 1 ] , -1.0 ) );
    row( ModulationEndLink , std::move( cf ) , -INF , double( L - 1 ) );
    }

   if( ! t ) {
    // the instant 0 has the (constant) band of the initial power, if the
    // unit is on at the beginning; if it is off, the band it restarts in
    // is free
    if( f_InitUpDownTime <= 0 )
     continue;
    const Index b0 = ( f_InitialPower <= B1 ) ? 0 :
                     ( ( f_InitialPower <= B2 ) ? 1 : 2 );
    // ... and if it shuts down at 0 it has no band at all, which is what
    // the term in the commitment leaves room for
    for( Index k = 0 ; k < 3 ; ++k ) {
     const double d0 = ( k == b0 ) ? 1.0 : 0.0;
     row( BandKeep , { coeff_pair( band( k , 0 ) , 1.0 ) ,
                       coeff_pair( & v_modulation_end[ 0 ] , -1.0 ) } ,
          -INF , d0 );
     // d0 - b^k_0 <= e_0 + ( 1 - u_0 ): the band is the initial one unless
     // a modulation ends at 0, and the unit that shuts down there has none
     row( BandKeep , { coeff_pair( band( k , 0 ) , -1.0 ) ,
                       coeff_pair( & v_modulation_end[ 0 ] , -1.0 ) ,
                       coeff_pair( & v_commitment[ 0 ] , 1.0 ) } ,
          -INF , 1 - d0 );
     row( BandMove , { coeff_pair( band( k , 0 ) , 1.0 ) ,
                       coeff_pair( & v_modulation_end[ 0 ] , 1.0 ) } ,
          -INF , 2 - d0 );
     }
    if( b0 == 2 )
     row( BandMove , { coeff_pair( band( 0 , 0 ) , 1.0 ) } , -INF , 0 );
    if( b0 == 0 )
     row( BandMove , { coeff_pair( band( 2 , 0 ) , 1.0 ) } , -INF , 0 );
    continue;
    }

   for( Index k = 0 ; k < 3 ; ++k ) {
    // the band does not change, unless a modulation ends or the unit was
    // off (in which case the band it restarts in is free)
    row( BandKeep , { coeff_pair( band( k , t ) , 1.0 ) ,
                      coeff_pair( band( k , t - 1 ) , -1.0 ) ,
                      coeff_pair( & v_modulation_end[ t ] , -1.0 ) ,
                      coeff_pair( & v_commitment[ t - 1 ] , 1.0 ) } ,
         -INF , 1 );
    row( BandKeep , { coeff_pair( band( k , t - 1 ) , 1.0 ) ,
                      coeff_pair( band( k , t ) , -1.0 ) ,
                      coeff_pair( & v_modulation_end[ t ] , -1.0 ) ,
                      coeff_pair( & v_commitment[ t ] , 1.0 ) } , -INF , 1 );
    // when a modulation ends the band does change
    row( BandMove , { coeff_pair( band( k , t ) , 1.0 ) ,
                      coeff_pair( band( k , t - 1 ) , 1.0 ) ,
                      coeff_pair( & v_modulation_end[ t ] , 1.0 ) } ,
         -INF , 2 );
    }
   // and it moves to an adjacent band, never across the whole range
   row( BandMove , { coeff_pair( band( 0 , t ) , 1.0 ) ,
                     coeff_pair( band( 2 , t - 1 ) , 1.0 ) } , -INF , 1 );
   row( BandMove , { coeff_pair( band( 2 , t ) , 1.0 ) ,
                     coeff_pair( band( 0 , t - 1 ) , 1.0 ) } , -INF , 1 );
   }

  add_static_constraint( BandChoice , "BandChoice_Nuclear" );
  add_static_constraint( BandPower , "BandPower_Nuclear" );
  add_static_constraint( ModulationEndLink , "ModulationEnd_Nuclear" );
  add_static_constraint( BandKeep , "BandKeep_Nuclear" );
  add_static_constraint( BandMove , "BandMove_Nuclear" );
  }

 // stability after a start-up - - - - - - - - - - - - - - - - - - - - - - -
 // a unit that starts up at t cannot begin a modulation before t + A, i.e.,
 // sum_{h=t}^{t+A-1} m_h <= K ( 1 - v_t ) with K the number of terms (the
 // term h = t being there only for uniformity, a start-up instant is not a
 // modulation step anyway); the tight form says the same one instant at a
 // time, m_h + v_t <= 1, which the relaxation cannot spread over the sum
 if( f_stability_after_start > 1 ) {
  StartUpStability.reserve( ( f_tight_rules ? f_stability_after_start : 1 )
                            * T );
  for( Index t = init_t ; t < T ; ++t ) {
   const Index h1 = std::min( t + f_stability_after_start , T );
   if( f_tight_rules )
    for( Index h = t + 1 ; h < h1 ; ++h )
     row( StartUpStability ,
          { coeff_pair( & v_modulation[ h ] , 1.0 ) ,
            coeff_pair( & v_start_up[ t - init_t ] , 1.0 ) } , -INF , 1 );
   else {
    LinearFunction::v_coeff_pair cf;
    for( Index h = t ; h < h1 ; ++h )
     cf.push_back( coeff_pair( & v_modulation[ h ] , 1.0 ) );
    cf.push_back( coeff_pair( & v_start_up[ t - init_t ] ,
                              double( h1 - t ) ) );
    row( StartUpStability , std::move( cf ) , -INF , double( h1 - t ) );
    }
   }
  add_static_constraint( StartUpStability , "StartUpStability_Nuclear" );
  }

 // modulations per day- - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( f_modulations_per_day >= 0 ) {
  ModulationsPerDayConst.reserve( days.size() );
  for( auto [ d0 , d1 ] : days ) {
   LinearFunction::v_coeff_pair cf;
   for( Index t = d0 ; t < d1 ; ++t )
    cf.push_back( coeff_pair( L > 1 ? & v_modulation_start[ t ]
                                    : & v_modulation[ t ] , 1.0 ) );
   row( ModulationsPerDayConst , std::move( cf ) , -INF ,
        double( f_modulations_per_day ) );
   }
  add_static_constraint( ModulationsPerDayConst ,
                         "ModulationsPerDay_Nuclear" );
  }

 // start-ups per day- - - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( f_start_ups_per_day >= 0 ) {
  StartUpsPerDayConst.reserve( days.size() );
  for( auto [ d0 , d1 ] : days ) {
   LinearFunction::v_coeff_pair cf;
   for( Index t = std::max( d0 , init_t ) ; t < d1 ; ++t )
    cf.push_back( coeff_pair( & v_start_up[ t - init_t ] , 1.0 ) );
   if( ! cf.empty() )
    row( StartUpsPerDayConst , std::move( cf ) , -INF ,
         double( f_start_ups_per_day ) );
   }
  add_static_constraint( StartUpsPerDayConst , "StartUpsPerDay_Nuclear" );
  }

 // deep decreases - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 // p_t + pd_t dd''_t >= pd_t ,
 // p_{t-1} - p_t - G_t dd'_t <= Dd_t with G_t = ( max{ D-_t , SD_t } -
 // Dd_t )^+ ,
 // dd_t - dd'_t - dd''_t - u_t >= -2
 // for the instants with an on predecessor (t = 0 only if on at the start)
 if( has_deep_decrease() ) {
  const Index t0 = ( f_InitUpDownTime > 0 ) ? 0 : 1;
  DeepLowConst.reserve( T );
  DeepDropConst.reserve( T );
  DeepLinkConst.reserve( T );
  for( Index t = t0 ; t < T ; ++t ) {
   const double pd = v_deep_threshold[ t ];
   const double Dd = v_deep_gradient[ t ];
   // the decrease is at most max{ D-_t , SD_t } [see above]
   const double Pm = std::max( v_DeltaRampDown[ t ] , v_ShutDownLimit[ t ] );
   // the output of an on unit is at least its minimum power, hence the tight
   // form p_t + ( pd - Pmin )^+ dd''_t - pd u_t >= 0, which the off unit
   // (p_t = 0) satisfies with any dd''_t; without the commitment the
   // coefficient of dd''_t has to be the whole pd, which is weaker
   if( f_tight_rules )
    row( DeepLowConst ,
         { coeff_pair( & v_active_power[ t ] , 1.0 ) ,
           coeff_pair( & v_deep_low[ t ] ,
                       std::max( pd - get_min_power( t ) , 0.0 ) ) ,
           coeff_pair( & v_commitment[ t ] , - pd ) } , 0 , INF );
   else
    row( DeepLowConst , { coeff_pair( & v_active_power[ t ] , 1.0 ) ,
                          coeff_pair( & v_deep_low[ t ] , pd ) } , pd , INF );
   LinearFunction::v_coeff_pair cf;
   if( t )
    cf.push_back( coeff_pair( & v_active_power[ t - 1 ] , 1.0 ) );
   cf.push_back( coeff_pair( & v_active_power[ t ] , -1.0 ) );
   cf.push_back( coeff_pair( & v_deep_drop[ t ] ,
                             - std::max( Pm - Dd , 0.0 ) ) );
   row( DeepDropConst , std::move( cf ) , -INF ,
        Dd - ( t ? 0 : f_InitialPower ) );
   row( DeepLinkConst , { coeff_pair( & v_deep[ t ] , 1.0 ) ,
                          coeff_pair( & v_deep_drop[ t ] , -1.0 ) ,
                          coeff_pair( & v_deep_low[ t ] , -1.0 ) ,
                          coeff_pair( & v_commitment[ t ] , -1.0 ) } ,
        -2 , INF );
   }
  add_static_constraint( DeepLowConst , "DeepLow_Nuclear" );
  add_static_constraint( DeepDropConst , "DeepDrop_Nuclear" );
  add_static_constraint( DeepLinkConst , "DeepLink_Nuclear" );

  // a decrease by at least the deep gradient is larger than what a stable
  // instant allows, hence a deep decrease is a downward modulation step:
  // dd_t <= d_t (<= m_t when the direction does not matter)
  if( f_tight_rules ) {
   DeepDownLink.reserve( T );
   for( Index t = t0 ; t < T ; ++t )
    if( v_deep_gradient[ t ] > v_modulation_ramp_down[ t ] )
     row( DeepDownLink ,
          { coeff_pair( & v_deep[ t ] , 1.0 ) ,
            coeff_pair( dir ? & v_modulation_down[ t ] : & v_modulation[ t ] ,
                        -1.0 ) } , -INF , 0 );
   if( ! DeepDownLink.empty() )
    add_static_constraint( DeepDownLink , "DeepDownLink_Nuclear" );
   }

  if( f_deep_decreases_per_day >= 0 ) {
   DeepDecreasesPerDayConst.reserve( days.size() );
   for( auto [ d0 , d1 ] : days ) {
    LinearFunction::v_coeff_pair cf;
    for( Index t = d0 ; t < d1 ; ++t )
     cf.push_back( coeff_pair( & v_deep[ t ] , 1.0 ) );
    row( DeepDecreasesPerDayConst , std::move( cf ) , -INF ,
         double( f_deep_decreases_per_day ) );
    }
   add_static_constraint( DeepDecreasesPerDayConst ,
                          "DeepDecreasesPerDay_Nuclear" );
   }
  }
 // the rows that are separated are those written on the start Variable,
 // which only exists if a modulation may last more than one instant
 if( f_tight_cuts && ( ! v_modulation_start.empty() ) ) {
  v_cut_done.assign( 3 * T , 0 );
  add_dynamic_constraint( Nuclear_cuts , "Nuclear_cuts" );
  }
 else
  f_tight_cuts = false;

 }  // end( NuclearUnitBlock::generate_operating_rules )

/*--------------------------------------------------------------------------*/

void NuclearUnitBlock::generate_dynamic_constraints( Configuration * dycc )
{
 // the Perspective Cuts of the thermal part
 ThermalUnitBlock::generate_dynamic_constraints( dycc );

 if( ! f_tight_cuts )
  return;

 const Index T = f_time_horizon;
 const Index L = f_max_modulation_length;
 const double tol = 1e-6;

 auto val = [ & ]( const ColVariable & v ) { return( v.get_value() ); };

 // add the row lhs <= sum coeff * var <= rhs to the separated ones
 auto cut = [ & ]( LinearFunction::v_coeff_pair && cf , double rhs ) {
  std::list< FRowConstraint > one( 1 );
  one.front().set_lhs( - Inf< double >() );
  one.front().set_rhs( rhs );
  one.front().set_function( new LinearFunction( std::move( cf ) , eNoMod ) );
  add_dynamic_constraints( Nuclear_cuts , one , eNoBlck );
  };

 for( Index t = 0 ; t < T ; ++t ) {
  // at most one start in any tau^M instants
  if( ( ! v_cut_done[ t ] ) && ( t + 1 < T ) ) {
   const Index K = std::min( Index( f_modulation_interval ) , T - t );
   if( K >= 2 ) {
    double lhs = 0;
    for( Index h = 0 ; h < K ; ++h )
     lhs += val( v_modulation_start[ t + h ] );
    if( lhs > 1 + tol ) {
     LinearFunction::v_coeff_pair cf;
     for( Index h = 0 ; h < K ; ++h )
      cf.push_back( coeff_pair( & v_modulation_start[ t + h ] , 1.0 ) );
     cut( std::move( cf ) , 1.0 );
     v_cut_done[ t ] = 1;
     }
    }
   }

  // a modulation ending at t forbids the starts up to t + tau^M - 1
  if( ( ! v_cut_done[ T + t ] ) && ( t + 1 < T ) ) {
   const Index K = std::min( Index( f_modulation_interval - 1 ) , T - 1 - t );
   if( K >= 1 ) {
    double lhs = val( v_modulation[ t ] ) - val( v_modulation[ t + 1 ] );
    for( Index h = 1 ; h <= K ; ++h )
     lhs += val( v_modulation_start[ t + h ] );
    if( lhs > 1 + tol ) {
     LinearFunction::v_coeff_pair cf;
     cf.push_back( coeff_pair( & v_modulation[ t ] , 1.0 ) );
     cf.push_back( coeff_pair( & v_modulation[ t + 1 ] , -1.0 ) );
     for( Index h = 1 ; h <= K ; ++h )
      cf.push_back( coeff_pair( & v_modulation_start[ t + h ] , 1.0 ) );
     cut( std::move( cf ) , 1.0 );
     v_cut_done[ T + t ] = 1;
     }
    }
   }

  // each step belongs to a modulation started in the last L instants
  if( ! v_cut_done[ 2 * T + t ] ) {
   const Index h0 = ( t >= L ? t - L + 1 : 0 );
   double lhs = val( v_modulation[ t ] );
   for( Index h = h0 ; h <= t ; ++h )
    lhs -= val( v_modulation_start[ h ] );
   if( lhs > tol ) {
    LinearFunction::v_coeff_pair cf;
    cf.push_back( coeff_pair( & v_modulation[ t ] , 1.0 ) );
    for( Index h = h0 ; h <= t ; ++h )
     cf.push_back( coeff_pair( & v_modulation_start[ h ] , -1.0 ) );
    cut( std::move( cf ) , 0.0 );
    v_cut_done[ 2 * T + t ] = 1;
    }
   }
  }
 }  // end( NuclearUnitBlock::generate_dynamic_constraints )

/*--------------------------------------------------------------------------*/

void NuclearUnitBlock::set_solution( void )
{
 // the thermal part: the start-up, the shut-down and the perspective-cut
 // auxiliaries out of ( p , u )
 ThermalUnitBlock::set_solution();

 auto Pi = get_const_active_power( 0 );
 auto Ci = get_const_commitment( 0 );
 if( ( ! Pi ) || ( ! Ci ) )
  return;             // no canonical variables: nothing to derive from

 // the start of a modulation: s_t = ( m_t - m_{t-1} )^+, the unit not being
 // in the middle of one at the beginning of the horizon
 if( ( ! v_modulation_start.empty() ) && ( ! v_modulation.empty() ) )
  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   const double m = v_modulation[ t ].get_value();
   const double mp = t ? v_modulation[ t - 1 ].get_value() : 0.0;
   v_modulation_start[ t ].set_value( std::max( m - mp , 0.0 ) );
   }

 // the deep decreases: the output is low, the drop is at least the deep
 // gradient, and the unit is on; the instant before the first one is the
 // initial power, and an initially off unit has no predecessor to drop from
 if( ! v_deep.empty() ) {
  const Index t0 = ( f_InitUpDownTime > 0 ) ? 0 : 1;
  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   if( t < t0 ) {
    v_deep[ t ].set_value( 0 );
    v_deep_drop[ t ].set_value( 0 );
    v_deep_low[ t ].set_value( 0 );
    continue;
    }
   const double p = Pi[ t ].get_value();
   const double pp = t ? Pi[ t - 1 ].get_value() : f_InitialPower;
   const bool on = Ci[ t ].get_value() > 0.5;
   const bool low = p <= v_deep_threshold[ t ] + 1e-9;
   const bool drop = pp - p >= v_deep_gradient[ t ] - 1e-9;
   v_deep_low[ t ].set_value( low ? 1 : 0 );
   v_deep_drop[ t ].set_value( drop ? 1 : 0 );
   v_deep[ t ].set_value( ( on && low && drop ) ? 1 : 0 );
   }
  }
 }  // end( NuclearUnitBlock::set_solution )

/*--------------------------------------------------------------------------*/

void NuclearUnitBlock::generate_objective( Configuration * objc )
{
 if( objective_generated() )
  return;

 ThermalUnitBlock::generate_objective( objc );

 // the costs of the downward modulation steps and of the deep decreases are
 // appended after all the Variable of the ThermalUnitBlock [see
 // objective_tail()], in this order
 DQuadFunction::v_coeff_triple vars;
 v_obj_tail.clear();
 if( ( ! v_modulation_down.empty() ) && ( ! v_down_modulation_cost.empty() ) )
  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   v_obj_tail.push_back( f_scale * v_down_modulation_cost[ t ] );
   vars.push_back( std::make_tuple( & v_modulation_down[ t ] ,
                                    v_obj_tail.back() , 0.0 ) );
   }
 if( ( ! v_deep.empty() ) && ( ! v_deep_cost.empty() ) )
  for( Index t = 0 ; t < f_time_horizon ; ++t ) {
   v_obj_tail.push_back( f_scale * v_deep_cost[ t ] );
   vars.push_back( std::make_tuple( & v_deep[ t ] , v_obj_tail.back() ,
                                    0.0 ) );
   }
 if( ! vars.empty() )
  static_cast< DQuadFunction * >( objective.get_function() )->add_variables(
                                                std::move( vars ) , eNoMod );

 }  // end( NuclearUnitBlock::generate_objective )

/*--------------------------------------------------------------------------*/

void NuclearUnitBlock::objective_tail_change( const DQuadFunction * qf ,
                                             Index first , Index last )
{
 // the costs of the operating rules cannot change via the abstract
 // representation: a change is fine as long as it leaves them as they are,
 // which is what a restore of the original costs does
 const Index base = qf->get_num_active_var() - v_obj_tail.size();
 for( Index i = first ; i < last ; ++i )
  if( ( qf->get_linear_coefficient( base + i ) != v_obj_tail[ i ] ) ||
      ( qf->get_quadratic_coefficient( base + i ) != 0 ) )
   throw( std::invalid_argument( "NuclearUnitBlock::add_Modification: the "
    "costs of the modulation and of the deep decreases cannot change" ) );

 }  // end( NuclearUnitBlock::objective_tail_change )

/*--------------------------------------------------------------------------*/

bool NuclearUnitBlock::is_feasible( bool useabstract , Configuration * fsbc )
{
 // retrieve the tolerance and the type of violation
 double tol = 0;
 bool rel_viol = true;

 // try to extract, from "c", the parameters that determine feasibility.
 // if it succeeds, it sets the values of the parameters and returns
 // true; otherwise, it returns false
 auto extract_parameters = [ & tol , & rel_viol ]( Configuration * c )
  -> bool {
  if( auto tc = dynamic_cast< SimpleConfiguration< double > * >( c ) ) {
   tol = tc->f_value;
   return( true );
   }
  if( auto tc = dynamic_cast< SimpleConfiguration<
                                    std::pair< double , int > > * >( c ) ) {
   tol = tc->f_value.first;
   rel_viol = tc->f_value.second;
   return( true );
   }
  return( false );
  };

 if( ( ! extract_parameters( fsbc ) ) && f_BlockConfig )
  // if the given Configuration is not valid, try the one from the BlockConfig
  extract_parameters( f_BlockConfig->f_is_feasible_Configuration );

 return( ThermalUnitBlock::is_feasible( useabstract )
   // Variable
   && ColVariable::is_feasible( v_modulation , tol )
   // Constraints
   && RowConstraint::is_feasible( Modulation_RampUp_Constraints , tol , rel_viol )
   && RowConstraint::is_feasible( Modulation_RampDown_Constraints , tol , rel_viol )
   && RowConstraint::is_feasible( NoDownModulation , tol , rel_viol )
   && RowConstraint::is_feasible( NoStartUpModulation , tol , rel_viol )
   && RowConstraint::is_feasible( ModulationConst , tol , rel_viol )
   && ColVariable::is_feasible( v_modulation_down , tol )
   && ColVariable::is_feasible( v_modulation_start , tol )
   && ColVariable::is_feasible( v_deep , tol )
   && ColVariable::is_feasible( v_deep_drop , tol )
   && ColVariable::is_feasible( v_deep_low , tol )
   && RowConstraint::is_feasible( DownModulationLink , tol , rel_viol )
   && RowConstraint::is_feasible( ModulationSameDirection , tol , rel_viol )
   && RowConstraint::is_feasible( Modulation_FullRampUp , tol , rel_viol )
   && RowConstraint::is_feasible( Modulation_FullRampDown , tol , rel_viol )
   && RowConstraint::is_feasible( ModulationStability , tol , rel_viol )
   && RowConstraint::is_feasible( StartUpStability , tol , rel_viol )
   && ColVariable::is_feasible( v_band , tol )
   && ColVariable::is_feasible( v_modulation_end , tol )
   && RowConstraint::is_feasible( BandChoice , tol , rel_viol )
   && RowConstraint::is_feasible( BandPower , tol , rel_viol )
   && RowConstraint::is_feasible( ModulationEndLink , tol , rel_viol )
   && RowConstraint::is_feasible( BandKeep , tol , rel_viol )
   && RowConstraint::is_feasible( BandMove , tol , rel_viol )
   && RowConstraint::is_feasible( ModulationMaxLength , tol , rel_viol )
   && RowConstraint::is_feasible( ModulationStartLink , tol , rel_viol )
   && RowConstraint::is_feasible( ModulationsPerDayConst , tol , rel_viol )
   && RowConstraint::is_feasible( StartUpsPerDayConst , tol , rel_viol )
   && RowConstraint::is_feasible( DeepDecreasesPerDayConst , tol ,
                                  rel_viol )
   && RowConstraint::is_feasible( DeepLowConst , tol , rel_viol )
   && RowConstraint::is_feasible( DeepDropConst , tol , rel_viol )
   && RowConstraint::is_feasible( DeepLinkConst , tol , rel_viol )
   && RowConstraint::is_feasible( DeepDownLink , tol , rel_viol )
   && RowConstraint::is_feasible( Nuclear_cuts , tol , rel_viol )
   );

 }  // end( NuclearUnitBlock::is_feasible )

/*--------------------------------------------------------------------------*/
/*--------------------- Methods for handling Solution ----------------------*/
/*--------------------------------------------------------------------------*/

Solution * NuclearUnitBlock::get_Solution( Configuration * csolc ,
					  bool emptys )
{
 Index wsol = 15;
 if( ( ! csolc ) && f_BlockConfig )
  csolc = f_BlockConfig->f_solution_Configuration;

 if( auto config = dynamic_cast< SimpleConfiguration< int > * >( csolc ) )
  wsol = config->f_value;

 // call the method of the base class
 auto sol = static_cast< NuclearUnitBlockSolution * >(
			      ThermalUnitBlock::get_Solution( csolc , true ) );

 // the modulation goes with the commitment: it is the same kind of
 // information, and there is no point in saving one without the other
 if( ( wsol & 2 ) && get_modulation() )
  sol->v_modulation.resize( get_time_horizon() );
 if( ( wsol & 2 ) && get_modulation_down() )
  sol->v_modulation_down.resize( get_time_horizon() );

 // the deep decreases are not implied by ( p , u , m , m^- ): the rows only
 // force them where the output really falls, so a Solver may leave one at 1
 // where it is not needed, and a Solution has to say which
 if( ( wsol & 2 ) && get_deep_decrease() ) {
  sol->v_deep.resize( get_time_horizon() );
  sol->v_deep_drop.resize( get_time_horizon() );
  sol->v_deep_low.resize( get_time_horizon() );
  }

 if( ! emptys )
  sol->read( this );

 return( sol );

 }  // end( NuclearUnitBlock::get_Solution )

/*--------------------------------------------------------------------------*/

UnitBlockSolution * NuclearUnitBlock::new_Solution( void ) const {
 return( new NuclearUnitBlockSolution() );
 }

/*--------------------------------------------------------------------------*/
/*-------- METHODS FOR LOADING, PRINTING & SAVING THE NuclearUnitBlock -----*/
/*--------------------------------------------------------------------------*/

void NuclearUnitBlock::serialize( netCDF::NcGroup & group ) const
{
 ThermalUnitBlock::serialize( group );

 // serialize scalar variables.
 ::serialize( group , "ModulationTime" , netCDF::NcUint() ,
        f_modulation_interval );
 ::serialize( group , "InitModulation" , netCDF::NcUint() ,
        f_initial_modulation );

 // serialize one-dimensional variables.
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
  else
   if( data.size() == NumberIntervals.getSize() )
    dimension = NumberIntervals;
   else
    if( data.size() != 1 ) {
     throw( std::logic_error(
    "NuclearUnitBlock::serialize: invalid dimension for variable " +
          var_name ) );
  }

  ::serialize( group , var_name , ncType , dimension , data ,
               allow_scalar_var );
  };

 serialize( "ModulationDeltaRampUp" , v_modulation_ramp_up );
 serialize( "ModulationDeltaRampDown" , v_modulation_ramp_down );

 // the operating rules, only when they are not the defaults
 if( f_max_modulation_length != 1 )
  ::serialize( group , "MaxModulationLength" , netCDF::NcUint() ,
               int( f_max_modulation_length ) );
 if( f_stability_after_start )
  ::serialize( group , "StabilityAfterStartUp" , netCDF::NcUint() ,
               int( f_stability_after_start ) );
 if( ! v_power_bands.empty() ) {
  auto dim = group.addDim( "NumberPowerBands" , 2 );
  group.addVar( "PowerBands" , netCDF::NcDouble() , dim ).putVar(
                                        { 0 } , { 2 } , v_power_bands.data() );
  }

 if( f_modulations_per_day >= 0 )
  ::serialize( group , "ModulationsPerDay" , netCDF::NcUint() ,
               f_modulations_per_day );
 if( f_start_ups_per_day >= 0 )
  ::serialize( group , "StartUpsPerDay" , netCDF::NcUint() ,
               f_start_ups_per_day );
 if( f_deep_decreases_per_day >= 0 )
  ::serialize( group , "DeepDecreasesPerDay" , netCDF::NcUint() ,
               f_deep_decreases_per_day );
 if( f_day_length )
  ::serialize( group , "DayLength" , netCDF::NcUint() ,
               int( f_day_length ) );
 serialize( "DownModulationCost" , v_down_modulation_cost );
 serialize( "DeepDecreaseThreshold" , v_deep_threshold );
 serialize( "DeepDecreaseGradient" , v_deep_gradient );
 serialize( "DeepDecreaseCost" , v_deep_cost );

 }  // end( NuclearUnitBlock::serialize )

/*--------------------------------------------------------------------------*/
/*------------------------ METHODS FOR CHANGING DATA -----------------------*/
/*----------------------------------------------------------------------------

void NuclearUnitBlock::add_Modification( sp_Mod mod, ChnlName chnl )
{
 if( mod->concerns_Block() ) {
  mod->concerns_Block( false );
  guts_of_add_Modification( mod.get() , chnl );
  }

 Block::add_Modification( mod, chnl );
 }

----------------------------------------------------------------------------*/

void NuclearUnitBlock::set_modulation_ramp_up( MF_dbl_it values ,
                 Subset && subset ,
                 bool ordered ,
                 ModParam issuePMod ,
                 ModParam issueAMod )
{
 static const std::string fn = "NuclearUnitBlock::set_modulation_ramp_up";

 if( subset.empty() )
  return;

 if( ! ordered )
  std::sort( subset.begin() , subset.end() );

 if( subset.back() >= f_time_horizon )
  throw( std::invalid_argument( fn + ": invalid index in subset" ) );

 if( identical( v_modulation_ramp_up , subset , values ) )  // no changes
  return;                                                   // return

 // check correctness of new values w.r.t. v_DeltaRampUp
 auto vit = values;
 for( auto t : subset ) {
  auto mrut = *( vit++ );
  if( mrut < 0 )
   throw( std::logic_error( fn + ": new modulation ramp up at time " +
                            std::to_string( t ) + " is " +
                            std::to_string( mrut ) + " < 0" ) );

  if( mrut > v_DeltaRampUp[ t ] )
   throw( std::logic_error( fn + ": new modulation ramp up at time " +
          std::to_string( t ) + " is " +
          std::to_string( mrut ) + "> ramp up = " +
          std::to_string( v_DeltaRampUp[ t ] ) ) );
  }

 if( not_dry_run( issuePMod ) )  // change the physical representation
  assign( v_modulation_ramp_up , subset , values );

 if( not_dry_run( issueAMod ) && constraints_generated() ) {
  // change the abstract representation
  // now change the corresponding Modulation_RampUp_Constraint[ t ]. note
  // that the \Delta^M_{t+} appears as the coefficient of u_{t-1} (if t > 0),
  // with opposite sign, and in the coefficient
  // ( \Delta_{t+} - \Delta^M_{t+} ) of m_t, again with opposite sign
  // these are respectively the coefficient 1 and 3 (the latter, only if
  // t > 0) of the LinearFunction in the FRowConstraint
  // if t == 0 then there is no term in u_{t-1} in the LinearFunction, but
  // \Delta^M_{t+} rather appears in the RHS, summed to f_InitialPower, if
  // u_{0 - 1} = 1, i.e., f_InitUpDownTime > 0

  // since several "abstract Modification" will be issued, pack them all into
  // a single GroupModification
  auto nAM = un_ModBlock( make_par( par2mod( issueAMod ) ,
            open_channel( par2chnl( issueAMod ) ) ) );

  // the instant 0, if it is there, is the only one whose \Delta^M_{t+} is in
  // the RHS rather than in a coefficient: it is dealt with here, out of the
  // loop, which then has no test of its own to do
  auto vt = values;
  if( subset.front() == 0 ) {
   const auto mru0 = *( vt++ );
   LF( Modulation_RampUp_Constraints[ 0 ].get_function()
       )->modify_coefficient( 1 , - ( v_DeltaRampUp[ 0 ] - mru0 ) , nAM );
   Modulation_RampUp_Constraints[ 0 ].set_rhs(
        f_InitialPower + ( f_InitUpDownTime > 0 ? mru0 : 0 ) , nAM );
   }

  // the two coefficients of each of the other instants, the modulation
  // variable being in position 1 of the LinearFunction and the commitment
  // one in position 2, go in a single Modification
  for( auto sit = subset.begin() + ( subset.front() == 0 ? 1 : 0 ) ;
       sit != subset.end() ; ++sit ) {
   const auto t = *sit;
   const auto mrut = *( vt++ );
   LF( Modulation_RampUp_Constraints[ t ].get_function()
       )->modify_coefficients( { - ( v_DeltaRampUp[ t ] - mrut ) , - mrut } ,
                               Range( 1 , 3 ) , nAM );
   }

  close_channel( par2chnl( nAM ) );  // at the end close the channel
  }

 if( issue_pmod( issuePMod ) )  // issue a physical Modification
  Block::add_Modification( std::make_shared< NuclearUnitBlockSbstMod >( this ,
                                              NuclearUnitBlockMod::eSetModDP ,
                                              std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );

 }  // end( NuclearUnitBlock::set_modulation_ramp_up )

/*--------------------------------------------------------------------------*/

void NuclearUnitBlock::set_modulation_ramp_up( MF_dbl_it values , Range rng ,
                 ModParam issuePMod ,
                 ModParam issueAMod )
{
 static const std::string fn = "NuclearUnitBlock::set_modulation_ramp_up";

 rng.second = std::min( rng.second , f_time_horizon );
 if( rng.second <= rng.first )
  return;

 // if nothing changes, return
 if( std::equal( values , values + ( rng.second - rng.first ) ,
                 v_modulation_ramp_up.begin() + rng.first ) )
  return;

 // check correctness of new values w.r.t. v_DeltaRampUp
 auto vit = values;
 for( auto t = rng.first ; t < rng.second ; ++t ) {
  auto mrut = *( vit++ );
  if( mrut < 0 )
   throw( std::logic_error( fn + ": new modulation ramp up at time " +
                            std::to_string( t ) + " is " +
                            std::to_string( mrut ) + " < 0" ) );

  if( mrut > v_DeltaRampUp[ t ] )
   throw( std::logic_error( fn + ": new modulation ramp up at time " +
          std::to_string( t ) + " is " +
          std::to_string( mrut ) + "> ramp up = " +
          std::to_string( v_DeltaRampUp[ t ] ) ) );
  }


 if( not_dry_run( issuePMod ) )  // change the physical representation
  std::copy( values , values + ( rng.second - rng.first ) ,
             v_modulation_ramp_up.begin() + rng.first );

 if( not_dry_run( issueAMod ) && constraints_generated() ) {
  // change the abstract representation
  // now change the corresponding Modulation_RampUp_Constraint[ t ]. note
  // that the \Delta^M_{t+} appears as the coefficient of u_{t-1} (if t > 0),
  // with opposite sign, and in the coefficient
  // ( \Delta_{t+} - \Delta^M_{t+} ) of m_t, again with opposite sign
  // these are respectively the coefficient 1 and 3 (the latter, only if
  // t > 0) of the LinearFunction in the FRowConstraint
  // if t == 0 then there is no term in u_{t-1} in the LinearFunction, but
  // \Delta^M_{t+} rather appears in the RHS, summed to f_InitialPower, if
  // u_{0 - 1} = 1, i.e., f_InitUpDownTime > 0

  // since several "abstract Modification" will be issued, pack them all into
  // a single GroupModification
  auto nAM = un_ModBlock( make_par( par2mod( issueAMod ) ,
            open_channel( par2chnl( issueAMod ) ) ) );

  // the instant 0, if it is in the range, is the only one whose
  // \Delta^M_{t+} is in the RHS rather than in a coefficient: it is dealt
  // with here, out of the loop, which then has no test of its own to do
  auto t = rng.first;
  if( t == 0 ) {
   const auto mru0 = *( values++ );
   LF( Modulation_RampUp_Constraints[ 0 ].get_function()
       )->modify_coefficient( 1 , - ( v_DeltaRampUp[ 0 ] - mru0 ) , nAM );
   Modulation_RampUp_Constraints[ 0 ].set_rhs(
        f_InitialPower + ( f_InitUpDownTime > 0 ? mru0 : 0 ) , nAM );
   ++t;
   }

  // the two coefficients of each of the other instants, the modulation
  // variable being in position 1 of the LinearFunction and the commitment
  // one in position 2, go in a single Modification
  for( ; t < rng.second ; ++t ) {
   const auto mrut = *( values++ );
   LF( Modulation_RampUp_Constraints[ t ].get_function()
       )->modify_coefficients( { - ( v_DeltaRampUp[ t ] - mrut ) , - mrut } ,
                               Range( 1 , 3 ) , nAM );
   }

  close_channel( par2chnl( nAM ) );  // at the end close the channel
  }

 if( issue_pmod( issuePMod ) )
  Block::add_Modification( std::make_shared< NuclearUnitBlockRngdMod >( this ,
                                      NuclearUnitBlockMod::eSetModDP , rng ) ,
                           Observer::par2chnl( issuePMod ) );

 }  // end( NuclearUnitBlock::set_modulation_ramp_up )

/*--------------------------------------------------------------------------*/

void NuclearUnitBlock::set_modulation_ramp_down( MF_dbl_it values ,
             Subset && subset ,
             bool ordered ,
             ModParam issuePMod ,
             ModParam issueAMod )
{
 static const std::string fn = "NuclearUnitBlock::set_modulation_ramp_down";

 if( subset.empty() )
  return;

 if( ! ordered )
  std::sort( subset.begin() , subset.end() );

 if( subset.back() >= f_time_horizon )
  throw( std::invalid_argument( fn + ": invalid index in subset" ) );

 if( identical( v_modulation_ramp_down , subset , values ) )  // no changes
  return;                                                     // return

 // check correctness of new values w.r.t. v_DeltaRampDown
 auto vit = values;
 for( auto t : subset ) {
  auto mrdt = *( vit++ );
  if( mrdt < 0 )
   throw( std::logic_error( fn + ": new modulation ramp down at time " +
                            std::to_string( t ) + " is " +
                            std::to_string( mrdt ) + " < 0" ) );

  if( mrdt > v_DeltaRampDown[ t ] )
   throw( std::logic_error( fn + ": new modulation ramp down at time " +
          std::to_string( t ) + " is " +
          std::to_string( mrdt ) + "> ramp down = " +
          std::to_string( v_DeltaRampDown[ t ] ) ) );
  }

 if( not_dry_run( issuePMod ) )  // change the physical representation
  assign( v_modulation_ramp_down , subset , values );

 if( not_dry_run( issueAMod ) && constraints_generated() ) {
  // change the abstract representation
  // now change the corresponding Modulation_RampDown_Constraint[ t ]. note
  // that the \Delta^M_{t-} appears as the coefficient of u_t, with opposite
  // sign, and in the coefficient ( \Delta_{t-} - \Delta^M_{t-} ) of m_t,
  // again with opposite sign. these are respectively the coefficient 1 and 2
  // of the LinearFunction in the FRowConstraint

  // since several "abstract Modification" will be issued, pack them all into
  // a single GroupModification
  auto nAM = un_ModBlock( make_par( par2mod( issueAMod ) ,
            open_channel( par2chnl( issueAMod ) ) ) );
  // the commitment variable is in position 1 of the LinearFunction and the
  // modulation one in position 2, so the two go in a single Modification
  for( auto t : subset ) {
   const auto mrdt = *( values++ );
   LF( Modulation_RampDown_Constraints[ t ].get_function()
       )->modify_coefficients( { - mrdt , - ( v_DeltaRampDown[ t ] - mrdt ) } ,
                               Range( 1 , 3 ) , nAM );
   }

  close_channel( par2chnl( nAM ) );  // at the end close the channel
  }

 if( issue_pmod( issuePMod ) )  // issue a physical Modification
  Block::add_Modification( std::make_shared< NuclearUnitBlockSbstMod >( this ,
                                              NuclearUnitBlockMod::eSetModDM ,
                                              std::move( subset ) ) ,
                           Observer::par2chnl( issuePMod ) );

 }  // end( NuclearUnitBlock::set_modulation_ramp_down )

/*--------------------------------------------------------------------------*/

void NuclearUnitBlock::set_modulation_ramp_down( MF_dbl_it values ,
             Range rng ,
             ModParam issuePMod ,
             ModParam issueAMod )
{
 static const std::string fn = "NuclearUnitBlock::set_modulation_ramp_down";

 rng.second = std::min( rng.second , f_time_horizon );
 if( rng.second <= rng.first )
  return;

 // if nothing changes, return
 if( std::equal( values , values + ( rng.second - rng.first ) ,
                 v_modulation_ramp_down.begin() + rng.first ) )
  return;

 // check correctness of new values w.r.t. v_DeltaRampDown
 auto vit = values;
 for( auto t = rng.first ; t < rng.second ; ++t ) {
  auto mrdt = *( vit++ );
  if( mrdt < 0 )
   throw( std::logic_error( fn + ": new modulation ramp down at time " +
                            std::to_string( t ) + " is " +
                            std::to_string( mrdt ) + " < 0" ) );

  if( mrdt > v_DeltaRampDown[ t ] )
   throw( std::logic_error( fn + ": new modulation ramp down at time " +
          std::to_string( t ) + " is " +
          std::to_string( mrdt ) + "> ramp up = " +
          std::to_string( v_DeltaRampDown[ t ] ) ) );
  }


 if( not_dry_run( issuePMod ) )  // change the physical representation
  std::copy( values , values + ( rng.second - rng.first ) ,
             v_modulation_ramp_down.begin() + rng.first );

 if( not_dry_run( issueAMod ) && constraints_generated() ) {
  // change the abstract representation
  // now change the corresponding Modulation_RampDown_Constraint[ t ]. note
  // that the \Delta^M_{t-} appears as the coefficient of u_t, with opposite
  // sign, and in the coefficient ( \Delta_{t-} - \Delta^M_{t-} ) of m_t,
  // again with opposite sign. these are respectively the coefficient 1 and 2
  // of the LinearFunction in the FRowConstraint

  // since several "abstract Modification" will be issued, pack them all into
  // a single GroupModification
  auto nAM = un_ModBlock( make_par( par2mod( issueAMod ) ,
            open_channel( par2chnl( issueAMod ) ) ) );

  // the commitment variable is in position 1 of the LinearFunction and the
  // modulation one in position 2, so the two go in a single Modification
  for( auto t = rng.first ; t < rng.second ; ++t ) {
   const auto mrdt = *( values++ );
   LF( Modulation_RampDown_Constraints[ t ].get_function()
       )->modify_coefficients( { - mrdt , - ( v_DeltaRampDown[ t ] - mrdt ) } ,
                               Range( 1 , 3 ) , nAM );
   }

  close_channel( par2chnl( nAM ) );  // at the end close the channel
  }

 if( issue_pmod( issuePMod ) )
  Block::add_Modification( std::make_shared< NuclearUnitBlockRngdMod >( this ,
                                      NuclearUnitBlockMod::eSetModDP , rng ) ,
                           Observer::par2chnl( issuePMod ) );

 }  // end( NuclearUnitBlock::set_modulation_ramp_up )

/*--------------------------------------------------------------------------*/

void NuclearUnitBlock::update_initial_power_in_cnstrs( ModParam issueAMod )
{
 // call the method of the base class to work on the original constraints
 ThermalUnitBlock::update_initial_power_in_cnstrs( issueAMod );

 // f_InitialPower influences the following new constraints of
 // NuclearUnitBlock, that have to be updated herein:
 // - the RHS of Modulation_RampUp_Constraints[ 0 ] is f_InitialPower +
 //   ( f_InitUpDownTime > 0 ? v_modulation_ramp_up[ 0 ] : 0 );
 // - the RHS of Modulation_RampDown_Constraints[ 0 ] is - f_InitialPower

 // several "abstract Modification" are issued, hence they are all packed
 // into a single GroupModification
 auto nAM = un_ModBlock( make_par( par2mod( issueAMod ) ,
                                   open_channel( par2chnl( issueAMod ) ) ) );

 Modulation_RampUp_Constraints[ 0 ].set_rhs( f_InitialPower +
      ( f_InitUpDownTime > 0 ? v_modulation_ramp_up[ 0 ] : 0 ) , nAM );

 Modulation_RampDown_Constraints[ 0 ].set_rhs( - f_InitialPower , nAM );

 // the rows of the operating rules that contain p_{-1}
 if( ! Modulation_FullRampUp.empty() ) {
  Modulation_FullRampUp[ 0 ].set_lhs( - full_ramp_up_const( 0 ) +
                                      f_InitialPower , nAM );
  Modulation_FullRampDown[ 0 ].set_lhs( - full_ramp_down_const( 0 ) -
                                        f_InitialPower , nAM );
  }
 if( ( ! DeepDropConst.empty() ) && ( f_InitUpDownTime > 0 ) )
  DeepDropConst[ 0 ].set_rhs( v_deep_gradient[ 0 ] - f_InitialPower , nAM );

 close_channel( par2chnl( nAM ) );  // at the end close the channel

 }  // end( NuclearUnitBlock::update_initial_power_in_constraints )

/*----------------------------------------------------------------------------

void NuclearUnitBlock::guts_of_add_Modification( p_Mod mod , ChnlName chnl )
{

 }  // end( NuclearUnitBlock::guts_of_add_Modification )

----------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------ METHODS OF NuclearUnitBlockSolution -------------------*/
/*--------------------------------------------------------------------------*/

void NuclearUnitBlockSolution::deserialize( const netCDF::NcGroup & group )
{
 // call the method of the base class
 ThermalUnitBlockSolution::deserialize( group );

 // deserialize the modulation - - - - - - - - - - - - - - - - - - - - - - -
 auto ncVar = group.getVar( "Modulation" );
 if( ncVar.isNull() )
  v_modulation.clear();
 else {
  v_modulation.resize( f_time_horizon );
  ncVar.getVar( { 0 } , { f_time_horizon } , v_modulation.data() );
  }

 ncVar = group.getVar( "ModulationDown" );
 if( ncVar.isNull() )
  v_modulation_down.clear();
 else {
  v_modulation_down.resize( f_time_horizon );
  ncVar.getVar( { 0 } , { f_time_horizon } , v_modulation_down.data() );
  }

 auto ds = [ & ]( const std::string & nm , std::vector< double > & v ) {
  auto nv = group.getVar( nm );
  if( nv.isNull() )
   v.clear();
  else {
   v.resize( f_time_horizon );
   nv.getVar( { 0 } , { f_time_horizon } , v.data() );
   }
  };
 ds( "DeepDecrease" , v_deep );
 ds( "DeepDrop" , v_deep_drop );
 ds( "DeepLow" , v_deep_low );

 }  // end( NuclearUnitBlockSolution::deserialize )

/*--------------------------------------------------------------------------*/

void NuclearUnitBlockSolution::read( const Block * block )
{
 auto NUB = dynamic_cast< const NuclearUnitBlock * >( block );
 if( ! NUB )
  throw( std::invalid_argument( "NuclearUnitBlockSolution::read: block "
				"is not a NuclearUnitBlock" ) );

 // call the method of the base class
 ThermalUnitBlockSolution::read( NUB );

 // read the modulation- - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_modulation.empty() )
  if( auto mi = NUB->get_const_modulation() )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    v_modulation[ t ] = mi[ t ].get_value();

 if( ! v_modulation_down.empty() )
  if( auto mi = NUB->get_const_modulation_down() )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    v_modulation_down[ t ] = mi[ t ].get_value();

 // the deep decreases
 auto rd = [ & ]( std::vector< double > & v , const ColVariable * x ) {
  if( v.empty() || ( ! x ) )
   return;
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   v[ t ] = x[ t ].get_value();
  };
 rd( v_deep , NUB->get_const_deep_decrease() );
 rd( v_deep_drop , NUB->get_const_deep_drop() );
 rd( v_deep_low , NUB->get_const_deep_low() );

 }  // end( NuclearUnitBlockSolution::read )

/*--------------------------------------------------------------------------*/

void NuclearUnitBlockSolution::write( Block * block )
{
 auto NUB = dynamic_cast< NuclearUnitBlock * >( block );
 if( ! NUB )
  throw( std::invalid_argument( "NuclearUnitBlockSolution::write: block "
				"is not a NuclearUnitBlock" ) );

 // call the method of the base class
 ThermalUnitBlockSolution::write( NUB );

 // write the modulation - - - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_modulation.empty() )
  if( auto mi = NUB->get_modulation() )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    ( mi++ )->set_value( v_modulation[ t ] );

 if( ! v_modulation_down.empty() )
  if( auto mi = NUB->get_modulation_down() )
   for( Index t = 0 ; t < f_time_horizon ; ++t )
    ( mi++ )->set_value( v_modulation_down[ t ] );

 // the start of the modulations and the deep decreases follow from what has
 // just been written, which the base class could not know when it derived
 // the thermal auxiliaries
 NUB->set_solution();

 // ... unless the deep decreases were saved, in which case they are the
 // ones to use: the rows only force them where the output really falls,
 // hence a Solver may have left one at 1 where deriving it gives 0
 auto wr = [ & ]( const std::vector< double > & v , ColVariable * x ) {
  if( v.empty() || ( ! x ) )
   return;
  for( Index t = 0 ; t < f_time_horizon ; ++t )
   ( x++ )->set_value( v[ t ] );
  };
 wr( v_deep , NUB->get_deep_decrease() );
 wr( v_deep_drop , NUB->get_deep_drop() );
 wr( v_deep_low , NUB->get_deep_low() );

 }  // end( NuclearUnitBlockSolution::write )

/*--------------------------------------------------------------------------*/

void NuclearUnitBlockSolution::serialize( netCDF::NcGroup & group ) const
{
 // call the method of the base class
 ThermalUnitBlockSolution::serialize( group );

 // serialize the modulation - - - - - - - - - - - - - - - - - - - - - - - -
 if( ! v_modulation.empty() )
  group.addVar( "Modulation" , netCDF::NcDouble() ,
		group.getDim( "TimeHorizon" ) ).putVar(
		     { 0 } , { f_time_horizon } , v_modulation.data() );

 if( ! v_modulation_down.empty() )
  group.addVar( "ModulationDown" , netCDF::NcDouble() ,
		group.getDim( "TimeHorizon" ) ).putVar(
		     { 0 } , { f_time_horizon } , v_modulation_down.data() );

 auto sv = [ & ]( const std::string & nm , const std::vector< double > & v ) {
  if( ! v.empty() )
   group.addVar( nm , netCDF::NcDouble() , group.getDim( "TimeHorizon" )
                 ).putVar( { 0 } , { f_time_horizon } , v.data() );
  };
 sv( "DeepDecrease" , v_deep );
 sv( "DeepDrop" , v_deep_drop );
 sv( "DeepLow" , v_deep_low );

 }  // end( NuclearUnitBlockSolution::serialize )

/*--------------------------------------------------------------------------*/

NuclearUnitBlockSolution * NuclearUnitBlockSolution::scale( double factor )
 const
{
 auto sol = clone();

 if( factor == 1 )
  return( sol );

 guts_of_scale( sol , factor );

 if( ! std::isnan( get_design() ) )
  sol->set_design( get_design() * factor );

 // every indicator scales like the rest of the Solution: the Objective pays
 // the downward steps and the deep decreases through their own Variable, so
 // that, as with the start-up ones of the base class, a convex combination
 // has to carry the average of the indicators and not those of one of its
 // constituents [see ThermalUnitBlockSolution::get_start_up()]
 auto times = [ factor ]( std::vector< double > v ) {
  for( auto & x : v )
   x *= factor;
  return( v );
  };
 sol->set_start_up( times( get_start_up() ) );
 sol->set_shut_down( times( get_shut_down() ) );
 for( auto * v : { & sol->v_modulation , & sol->v_modulation_down ,
                   & sol->v_deep , & sol->v_deep_drop , & sol->v_deep_low } )
  for( auto & x : *v )
   x *= factor;

 return( sol );

 }  // end( NuclearUnitBlockSolution::scale )

/*--------------------------------------------------------------------------*/

void NuclearUnitBlockSolution::sum( const Solution * solution ,
				    double multiplier )
{
 // call the method of the base class
 ThermalUnitBlockSolution::sum( solution , multiplier );

 auto NUBS = dynamic_cast< const NuclearUnitBlockSolution * >( solution );
 if( ! NUBS )
  throw( std::invalid_argument( "NuclearUnitBlockSolution::sum: solution "
                                "not a NuclearUnitBlockSolution" ) );

 // the indicators of the rules are summed as those of the base class [see
 // scale()]: a Solution that does not carry one of them is summed with one
 // that does only if the two agree on which it has
 auto add = [ & ]( std::vector< double > & v ,
                   const std::vector< double > & w , const char * what ) {
  if( v.size() != w.size() )
   throw( std::invalid_argument( std::string(
    "NuclearUnitBlockSolution::sum: inconsistent " ) + what + ", " +
    std::to_string( v.size() ) + " against " + std::to_string( w.size() ) ) );
  for( Index i = 0 ; i < v.size() ; ++i )
   v[ i ] += w[ i ] * multiplier;
  };
 add( v_modulation , NUBS->v_modulation , "modulation indicators" );
 add( v_modulation_down , NUBS->v_modulation_down ,
      "downward modulation indicators" );
 add( v_deep , NUBS->v_deep , "deep-decrease indicators" );
 add( v_deep_drop , NUBS->v_deep_drop , "deep-decrease drop indicators" );
 add( v_deep_low , NUBS->v_deep_low , "deep-decrease low indicators" );

 }  // end( NuclearUnitBlockSolution::sum )

/*--------------------------------------------------------------------------*/

NuclearUnitBlockSolution * NuclearUnitBlockSolution::clone( bool empty ) const
{
 auto sol = new NuclearUnitBlockSolution();

 if( ! empty ) {
  guts_of_clone( sol );
  sol->set_design( get_design() );
  // the start-up and shut-down indicators of the base class, which its own
  // clone() copies and guts_of_clone() does not: without them a scale(),
  // which is a clone(), would return a Solution of a different shape, and
  // the convex combinations of a Lagrangian approach would refuse to sum it
  sol->set_start_up( std::vector< double >( get_start_up() ) );
  sol->set_shut_down( std::vector< double >( get_shut_down() ) );
  sol->v_modulation = v_modulation;
  sol->v_modulation_down = v_modulation_down;
  sol->v_deep = v_deep;
  sol->v_deep_drop = v_deep_drop;
  sol->v_deep_low = v_deep_low;
  }

 return( sol );

 }  // end( NuclearUnitBlockSolution::clone )

/*--------------------------------------------------------------------------*/
/*------------------- End File NuclearUnitBlock.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
