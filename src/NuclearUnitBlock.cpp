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
 * Copyright &copy by Antonio Frangioni
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "LinearFunction.h"

#include "NuclearUnitBlock.h"

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
  vec[ t ] = *(it++);
 }

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS OF NuclearUnitBlock ----------------------*/
/*--------------------------------------------------------------------------*/

NuclearUnitBlock::~NuclearUnitBlock()
{
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
#ifndef NDEBUG
 // check all expected variables, comprised those of the base class: see
 // ThermalUnitBlock::deserialize() for the rationale
 const std::vector< std::string > expected_vars = { "MinPower" , "MaxPower" ,
  "DeltaRampUp" , "DeltaRampDown" , "PrimaryRho" , "SecondaryRho" ,
  "LinearTerm" , "QuadTerm" , "ConstTerm" , "StartUpCost" ,
  "FixedConsumption" , "InertiaCommitment" , "InitialPower" , "MinUpTime" ,
  "MinDownTime" , "InitUpDownTime" , "Availability" , "ModulationTime" ,
  "InitModulation" , "ModulationDeltaRampUp" , "ModulationDeltaRampDown" };

 check_variables( group , expected_vars , std::cerr );
#endif

 // call the method of the base class 
 ThermalUnitBlock::deserialize( group );

 // check that DeltaRampUp/Down are defined
 if( v_DeltaRampUp.empty() )
  throw( std::invalid_argument(
		"NuclearUnitBlock::deserialize: DeltaRampUp not present" ) );

 if( v_DeltaRampDown.empty() )
  throw( std::invalid_argument(
	      "NuclearUnitBlock::deserialize: DeltaRampDown not present" ) );

 // load optional variables ModulationTime and InitModulation or give them
 // default values, check that the values are logically correct

 if( ! ::deserialize( group , f_modulation_interval , "ModulationTime" ) )
  f_modulation_interval = 2;

 if( f_modulation_interval < 2 )
  throw( std::invalid_argument(
		    "NuclearUnitBlock::deserialize: ModulationTime < 2" ) );
 
 if( ! ::deserialize( group , f_initial_modulation , "InitModulation" ) )
  f_initial_modulation = f_modulation_interval;

 if( f_modulation_interval < 1 )
  throw( std::invalid_argument(
		    "NuclearUnitBlock::deserialize: InitModulation < 1" ) );

 // load mandatory variables ModulationDeltaRampUp and
 // ModulationDeltaRampDown, create the expanded vectors (if needed)
 
 ::deserialize( group , "ModulationDeltaRampUp" , v_modulation_ramp_up ,
		false );
 ::deserialize( group , "ModulationDeltaRampDown" , v_modulation_ramp_down ,
		false );

 // decompress vectors
 decompress_vector( v_modulation_ramp_up );
 decompress_vector( v_modulation_ramp_down );

 // check consistency of v_modulation_ramp_up and v_modulation_ramp_down
 check_modulation_consistency();

 }  // end( NuclearUnitBlock::deserialize )

/*--------------------------------------------------------------------------*/

void NuclearUnitBlock::check_modulation_consistency( void ) const
{
 static const std::string fn =
                             "NuclearUnitBlock::check_modulation_consistency";

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
                            "> ramp up = " +
			    std::to_string( v_DeltaRampUp[ t ] ) ) );

  if( v_modulation_ramp_down[ t ] < 0 )
   throw( std::logic_error( fn + ": modulation ramp down at time " +
                            std::to_string( t ) + " is " +
                            std::to_string( v_modulation_ramp_down[ t ] ) +
                            " < 0" ) );

  if( v_modulation_ramp_down[ t ] > v_DeltaRampUp[ t ] )
   throw( std::logic_error( fn + ": modulation ramp down at time " +
			    std::to_string( t ) + " is " +
			    std::to_string( v_modulation_ramp_down[ t ] ) +
                            "> ramp up = " +
			    std::to_string( v_DeltaRampDown[ t ] ) ) );

  }  // end( for t )
 }  // end( NuclearUnitBlock::check_modulation_consistency )

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

 for( int t = 0 ; t < f_modulation_interval - f_initial_modulation ; ) {
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
 
 } // end( NuclearUnitBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void NuclearUnitBlock::generate_abstract_constraints( Configuration * stcc ) {

 if( constraints_generated() )
  return; // constraints have already been generated

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

 Modulation_RampUp_Constraints.resize( f_time_horizon );

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {
  Index np = t ? 5 : 3;
  if( t < init_t )
   --np;
  LinearFunction::v_coeff_pair cf( np );
  double RHS = 0;
  auto cfit = cf.begin();

  *(cfit++) = coeff_pair( & v_active_power[ t ] , 1.0 );
  *(cfit++) = coeff_pair( & v_modulation[ t ] ,
			  - ( v_DeltaRampUp[ t ] - v_modulation_ramp_up[ t ] )
			);
  // the two terms "- p_{t-1}" and "- \Delta^M_{t+} u_{t-1}" only exist if
  // t > 0, as otherwise p_{t-1} and u_{t-1} are undefined
  if( t ) {
   *(cfit++) = coeff_pair( & v_commitment[ t - 1 ] ,
			 - v_modulation_ramp_up[ t ] );

   *(cfit++) = coeff_pair( & v_active_power[ t - 1 ] , -1.0 );
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
   *cfit = coeff_pair( & v_start_up[ t - init_t ] , v_StartUpLimit[ t ] );

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

 Modulation_RampDown_Constraints.resize( f_time_horizon );

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {
  Index np = t ? 5 : 4;
  if( t < init_t )
   --np;
  LinearFunction::v_coeff_pair cf( np );
  auto cfit = cf.begin();

  *(cfit++) = coeff_pair( & v_active_power[ t ] , -1.0 );
  *(cfit++) = coeff_pair( & v_commitment[ t ] ,
			- v_modulation_ramp_down[ t ] );
  *(cfit++) = coeff_pair( & v_modulation[ t ] ,
			- ( v_DeltaRampDown[ t ] -
			    v_modulation_ramp_down[ t ] ) );

  // the terms "p_{t-1}" only exists if t > 0, as otherwise p_{t-1} is
  // undefined
  if( t )
   *(cfit++) = coeff_pair( & v_active_power[ t - 1 ] , 1.0 );

  // the term - \bar{l}_t v_t only exist if t >= init_t, as for t < init_t
  // the commitment status if fixed and shut-downs are not allowed, hence
  // the corresponding shut-down variables are not even defined
  if( t >= init_t )
   *cfit = coeff_pair( & v_shut_down[ t - init_t ] , v_ShutDownLimit[ t ] );
 
  Modulation_RampDown_Constraints[ t ].set_lhs( - Inf< double >() );
  // if t == 0, the "p_{t-1}" term is fixed and equal to f_InitialPower, so
  // so there is no explicit term in the constraint (since the variable does
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
 // start-up can ever occur; in facy, the start-up variables are not even
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
 Index first_c = std::max( f_modulation_interval - f_initial_modulation ,
			   int( 0 ) );
 if( f_InitUpDownTime <= 0 )
  first_c = std::max( first_c , init_t );
 ++first_c;

 ModulationConst.resize( f_time_horizon - first_c );

 for( Index t = first_c ; t < f_time_horizon ; ++t ) {
  Index h = std::max( int( 0 ) , int( t ) - f_modulation_interval + 1 );
  LinearFunction::v_coeff_pair cf( t - h + 1 );

  for( auto cfit = cf.begin() ; h <= t ; )
   *(cfit++) = coeff_pair( & v_modulation[ h++ ] , 1.0 );

  ModulationConst[ t - first_c ].set_lhs( - Inf< double >() );
  ModulationConst[ t - first_c ].set_rhs( 1.0 );
  ModulationConst[ t - first_c ].set_function(
				    new LinearFunction( std::move( cf ) ) );
  }

 add_static_constraint( ModulationConst , "ModulationConst_Nuclear" );

 } // end( NuclearUnitBlock::generate_abstract_constraints )

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
	 && RowConstraint::is_feasible( Modulation_RampUp_Constraints ,
					tol , rel_viol )
	 && RowConstraint::is_feasible( Modulation_RampDown_Constraints ,
					tol , rel_viol )
	 && RowConstraint::is_feasible( NoDownModulation , tol , rel_viol )
	 && RowConstraint::is_feasible( NoStartUpModulation , tol , rel_viol )
	 && RowConstraint::is_feasible( ModulationConst , tol , rel_viol )
	 );

 }  // end( NuclearUnitBlock::is_feasible )

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
 auto serialize = [ & group , & TimeHorizon , & NumberIntervals ](
	const std::string & var_name , const std::vector< double > & data ,
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
  std::sort( subset.begin(), subset.end() );

 if( subset.back() >= f_time_horizon )
  throw( std::invalid_argument( fn + ": invalid index in subset" ) );

 if( identical( v_modulation_ramp_up , subset , values ) )  // no changes
  return;                                                   // return

 // check correctness of new values w.r.t. v_DeltaRampUp
 auto vit = values;
 for( auto t : subset ) {
  auto mrut = *(vit++);
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
  // that the \Delta^M_{t+} appears as the coeeficient of u_{t-1} (if t > 0),
  // with opposite sign, and in the coeeficient
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

  // TODO: make it more efficient by treating the t == 0 case offline
  for( auto t : subset ) {
   auto mrut = *(values++);
   auto lf = LF( Modulation_RampUp_Constraints[ t ].get_function() );

   // TODO: make it more efficient by calling modify_coefficients( subset )
   // the modulation variable is in position 1 in the lf
   lf->modify_coefficient( 1 , - ( v_DeltaRampUp[ t ] - mrut ) , nAM );

   if( t )  // for t > 0 the commitment variable is in position 2 in the lf
    lf->modify_coefficient( 2 , - mrut , nAM );
   else     // for t == 0, \Delta^M_{t+} is in the RHS
    Modulation_RampUp_Constraints[ t ].set_rhs(
	       f_InitialPower + ( f_InitUpDownTime > 0 ? mrut : 0 ) , nAM );
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
  auto mrut = *(vit++);
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

  // TODO: make it more efficient by treating the t == 0 case offline
  for( auto t = rng.first ; t < rng.second ; ++t ) {
   auto mrut = *(values++);
   auto lf = LF( Modulation_RampUp_Constraints[ t ].get_function() );

   // TODO: make it more efficient by calling modify_coefficients( subset )
   // the modulation variable is in position 1 in the lf
   lf->modify_coefficient( 1 , - ( v_DeltaRampUp[ t ] - mrut ) , nAM );

   if( t )  // for t > 0 the commitment variable is in position 2 in the lf
    lf->modify_coefficient( 2 , - mrut , nAM );
   else     // for t == 0, \Delta^M_{t+} is in the RHS
    Modulation_RampUp_Constraints[ t ].set_rhs(
	       f_InitialPower + ( f_InitUpDownTime > 0 ? mrut : 0 ) , nAM );
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
  std::sort( subset.begin(), subset.end() );

 if( subset.back() >= f_time_horizon )
  throw( std::invalid_argument( fn + ": invalid index in subset" ) );

 if( identical( v_modulation_ramp_down , subset , values ) )  // no changes
  return;                                                     // return

 // check correctness of new values w.r.t. v_DeltaRampDown
 auto vit = values;
 for( auto t : subset ) {
  auto mrdt = *(vit++);
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
  // that the \Delta^M_{t-} appears as the coeeficient of u_t, with opposite
  // sign, and in the coeeficient ( \Delta_{t-} - \Delta^M_{t-} ) of m_t,
  // again with opposite sign. these are respectively the coefficient 1 and 2
  // of the LinearFunction in the FRowConstraint

  // since several "abstract Modification" will be issued, pack them all into
  // a single GroupModification
  auto nAM = un_ModBlock( make_par( par2mod( issueAMod ) ,
				    open_channel( par2chnl( issueAMod ) ) ) );
  for( auto t : subset ) {
   auto mrdt = *(values++);
   auto lf = LF( Modulation_RampDown_Constraints[ t ].get_function() );

   // TODO: make it more efficient by calling modify_coefficients( range )
   // the commitment variable is in position 1 in the lf
   lf->modify_coefficient( 1 , - mrdt , nAM );
   // the modulation variable is in position 2 in the lf
   lf->modify_coefficient( 2 , - ( v_DeltaRampDown[ t ] - mrdt ) , nAM );
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
  auto mrdt = *(vit++);
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
  // that the \Delta^M_{t-} appears as the coeeficient of u_t, with opposite
  // sign, and in the coeeficient ( \Delta_{t-} - \Delta^M_{t-} ) of m_t,
  // again with opposite sign. these are respectively the coefficient 1 and 2
  // of the LinearFunction in the FRowConstraint

  // since several "abstract Modification" will be issued, pack them all into
  // a single GroupModification
  auto nAM = un_ModBlock( make_par( par2mod( issueAMod ) ,
				    open_channel( par2chnl( issueAMod ) ) ) );

  for( auto t = rng.first ; t < rng.second ; ++t ) {
   auto mrdt = *(values++);
   auto lf = LF( Modulation_RampUp_Constraints[ t ].get_function() );

   // TODO: make it more efficient by calling modify_coefficients( range )
   // the commitment variable is in position 1 in the lf
   lf->modify_coefficient( 1 , - mrdt , nAM );
   // the modulation variable is in position 2 in the lf
   lf->modify_coefficient( 2 , - ( v_DeltaRampDown[ t ] - mrdt ) , nAM );
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

 // TODO: pack the two Modificaton into a GroupModification
 Modulation_RampUp_Constraints[ 0 ].set_rhs( f_InitialPower +
      ( f_InitUpDownTime > 0 ? v_modulation_ramp_up[ 0 ] : 0 ) , issueAMod );

 Modulation_RampDown_Constraints[ 0 ].set_rhs( - f_InitialPower ,
					       issueAMod );

 }  // end( NuclearUnitBlock::update_initial_power_in_constraints )

/*----------------------------------------------------------------------------

void NuclearUnitBlock::guts_of_add_Modification( p_Mod mod , ChnlName chnl )
{

 }  // end( NuclearUnitBlock::guts_of_add_Modification )

----------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------- End File NuclearUnitBlock.cpp ------------------------*/
/*--------------------------------------------------------------------------*/
