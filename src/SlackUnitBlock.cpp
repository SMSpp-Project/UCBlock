/*--------------------------------------------------------------------------*/
/*------------------------ File SlackUnitBlock.cpp -------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the SlackUnitBlock class, which derives from UnitBlock
 * [see UnitBlock.h] and implements a "slack" unit; a (typically, fictitious)
 * unit capable of producing (typically, a large amount of) active power
 * and/or primary/secondary reserve and/or inertia at any time period
 * completely independently from each other and from all other time periods,
 * albeit at a (typically, huge) cost.
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
 * \copyright &copy; by Antonio Frangioni, Ali Ghezelsoflu,
 *                      Rafael Durbano Lobato
 */
/*--------------------------------------------------------------------------*/
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "SlackUnitBlock.h"

#include "LinearFunction.h"

#include "FRealObjective.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*----------------------------- STATIC MEMBERS -----------------------------*/
/*--------------------------------------------------------------------------*/

// register SlackUnitBlock to the Block factory

SMSpp_insert_in_factory_cpp_1( SlackUnitBlock );

/*--------------------------------------------------------------------------*/
/*----------------------- METHODS OF SlackUnitBlock ------------------------*/
/*--------------------------------------------------------------------------*/

SlackUnitBlock::~SlackUnitBlock()
{
 Constraint::clear( Secondary_Spinning_Reserve_Bound_Const );
 Constraint::clear( Primary_Spinning_Reserve_Bound_Const );
 Constraint::clear( ActivePower_Bound_Const );

 Constraint::clear( Inertia_Bound_Const );

 objective.clear();
}

/*--------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/

void SlackUnitBlock::deserialize( const netCDF::NcGroup & group )
{

#ifndef NDEBUG
 static std::vector< std::string > expected_dims = { "TimeHorizon" ,
                                                     "NumberIntervals" };
 check_dimensions( group , expected_dims , std::cerr );

 static std::vector< std::string > expected_vars = { "MaxPower" ,
                                                     "MaxPrimaryPower" ,
                                                     "MaxSecondaryPower" ,
                                                     "ActivePowerCost" ,
                                                     "PrimaryCost" ,
                                                     "SecondaryCost" ,
                                                     "InertiaCost" ,
                                                     "MaxInertia" };
 check_variables( group , expected_vars , std::cerr );
#endif

 // Optional variables
 ::deserialize( group , "MaxPower" , v_MaxPower );
 ::deserialize( group , "MaxPrimaryPower" , v_MaxPrimaryPower );
 ::deserialize( group , "MaxSecondaryPower" , v_MaxSecondaryPower );
 ::deserialize( group , "ActivePowerCost" , v_ActivePowerCost );
 ::deserialize( group , "PrimaryCost" , v_PrimaryCost );
 ::deserialize( group , "SecondaryCost" , v_SecondaryCost );
 ::deserialize( group , "InertiaCost" , v_InertiaCost );
 ::deserialize( group , "MaxInertia" , v_MaxInertia );

 // Deserialize data from the base class
 UnitBlock::deserialize( group );

 // Decompress vectors
 decompress_vector( v_MaxPower );
 decompress_vector( v_MaxPrimaryPower );
 decompress_vector( v_MaxSecondaryPower );
 decompress_vector( v_ActivePowerCost );
 decompress_vector( v_PrimaryCost );
 decompress_vector( v_SecondaryCost );
 decompress_vector( v_InertiaCost );
 decompress_vector( v_MaxInertia );

}  // end( SlackUnitBlock::deserialize )

/*--------------------------------------------------------------------------*/

void SlackUnitBlock::generate_abstract_variables( Configuration * stvv )
{
 if( variables_generated() )  // variables have already been generated
  return;                     // nothing to do

 // Commitment Variable
 if( reserve_vars & 4u ) {  // if UCBlock has inertia demand variables
  if( ! v_MaxInertia.empty() ) {  // if unit produces any inertia reserve
   v_commitment.resize( f_time_horizon );
   for( auto & i : v_commitment )
    i.set_type( ColVariable::kPosUnitary );
   add_static_variable( v_commitment , "u_slack" );
  }
 }

 // Active Power Variable
 v_active_power.resize( f_time_horizon );
 for( auto & var : v_active_power )
  var.set_type( ColVariable::kNonNegative );
 add_static_variable( v_active_power , "p_slack" );

 // Primary Spinning Reserve Variable
 if( reserve_vars & 1u ) {  // if UCBlock has primary demand variables
  if( ! v_MaxPrimaryPower.empty() ) {  // if unit produces any primary reserve
   v_primary_spinning_reserve.resize( f_time_horizon );
   for( auto & var : v_primary_spinning_reserve )
    var.set_type( ColVariable::kNonNegative );
   add_static_variable( v_primary_spinning_reserve , "pr_slack" );
  }
 }

 // Secondary Spinning Reserve Variable
 if( reserve_vars & 2u ) {  // if UCBlock has secondary demand variables
  if( ! v_MaxSecondaryPower.empty() ) {  // if unit produces any secondary reserve
   v_secondary_spinning_reserve.resize( f_time_horizon );
   for( auto & var : v_secondary_spinning_reserve )
    var.set_type( ColVariable::kNonNegative );
   add_static_variable( v_secondary_spinning_reserve , "sr_slack" );
  }
 }

 set_variables_generated();

}  // end( SlackUnitBlock::generate_abstract_variables )

/*--------------------------------------------------------------------------*/

void SlackUnitBlock::generate_abstract_constraints( Configuration * stcc )
{
 if( constraints_generated() )  // constraints have already been generated
  return;                       // nothing to do

 bool generate_ZOConstraints = false;
 if( ( ! stcc ) && f_BlockConfig )
  stcc = f_BlockConfig->f_static_constraints_Configuration;
 if( auto sci = dynamic_cast< SimpleConfiguration< int > * >( stcc ) )
  generate_ZOConstraints = sci->f_value;

 // Initializing active power bounds constraints
 if( ActivePower_Bound_Const.size() != f_time_horizon ) {
  // this should only happen once
  assert( ActivePower_Bound_Const.empty() );

  ActivePower_Bound_Const.resize( f_time_horizon );
 }

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {

  if( ! v_MaxPower.empty() )
   ActivePower_Bound_Const[ t ].set_rhs( v_MaxPower[ t ] );
  else
   ActivePower_Bound_Const[ t ].set_rhs( 0.0 );

  ActivePower_Bound_Const[ t ].set_variable( &v_active_power[ t ] );
 }

 add_static_constraint( ActivePower_Bound_Const ,
                        "ActivePowerBound_Slack" );

 /*--------------------------------------------------------------------------*/

 // Initializing primary spinning reserve bounds constraints
 if( reserve_vars & 1u ) {  // if UCBlock has primary demand variables
  if( ! v_MaxPrimaryPower.empty() ) {  // if unit produces any primary reserve

   if( Primary_Spinning_Reserve_Bound_Const.size() != f_time_horizon ) {
    // this should only happen once
    assert( Primary_Spinning_Reserve_Bound_Const.empty() );

    Primary_Spinning_Reserve_Bound_Const.resize( f_time_horizon );
   }

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {
    Primary_Spinning_Reserve_Bound_Const[ t ].set_rhs(
     v_MaxPrimaryPower[ t ] );
    Primary_Spinning_Reserve_Bound_Const[ t ].set_variable(
     &v_primary_spinning_reserve[ t ] );
   }

   add_static_constraint( Primary_Spinning_Reserve_Bound_Const ,
                          "PrimarySpinningReserveBound_Slack" );
  }
 }

 /*--------------------------------------------------------------------------*/

 // Initializing secondary spinning reserve bounds constraints
 if( reserve_vars & 2u ) {  // if UCBlock has secondary demand variables
  if( ! v_MaxSecondaryPower.empty() ) {  // if unit produces any secondary reserve
   if( Secondary_Spinning_Reserve_Bound_Const.size() != f_time_horizon ) {
    // this should only happen once
    assert( Secondary_Spinning_Reserve_Bound_Const.empty() );

    Secondary_Spinning_Reserve_Bound_Const.resize( f_time_horizon );
   }

   for( Index t = 0 ; t < f_time_horizon ; ++t ) {
    Secondary_Spinning_Reserve_Bound_Const[ t ].set_rhs(
     v_MaxSecondaryPower[ t ] );
    Secondary_Spinning_Reserve_Bound_Const[ t ].set_variable(
     &v_secondary_spinning_reserve[ t ] );
   }

   add_static_constraint( Secondary_Spinning_Reserve_Bound_Const ,
                          "SecondarySpinningReserveBound_Slack" );
  }
 }

 /*------------------------------ ZOConstraint ------------------------------*/

 if( generate_ZOConstraints ) {

  // the commitment bound constraints
  if( reserve_vars & 4u ) {
   Inertia_Bound_Const.resize( f_time_horizon );

   for( Index t = 0 ; t < f_time_horizon ; ++t )
    Inertia_Bound_Const[ t ].set_variable( &v_commitment[ t ] );

   add_static_constraint( Inertia_Bound_Const , "Inertia_bound_Slack" );
  }
 }

 set_constraints_generated();

}  // end( SlackUnitBlock::generate_abstract_constraints )


/*--------------------------------------------------------------------------*/

void SlackUnitBlock::generate_objective( Configuration * objc )
{
 if( objective_generated() )  // Objective has already been generated
  return;                     // nothing to do

 if( reserve_vars & 4u )
  if( v_commitment.size() != f_time_horizon )
   throw( std::logic_error(
    "SlackUnitBlock::generate_objective: v_commitment must have "
    "size equal to the time horizon." ) );

 if( v_active_power.size() != f_time_horizon )
  throw( std::logic_error( "SlackUnitBlock::generate_objective: "
                           "v_active_power must have size equal to the time "
                           "horizon." ) );

 if( reserve_vars & 1u )  // if UCBlock has primary demand variables
  if( ! v_MaxPrimaryPower.empty() )
   if( v_primary_spinning_reserve.size() != f_time_horizon )
    throw( std::logic_error( "SlackUnitBlock::generate_objective: "
                             "v_primary_spinning_reserve must have size equal"
                             " to the time horizon." ) );

 if( reserve_vars & 2u )  // if UCBlock has secondary demand variables
  if( ! v_MaxSecondaryPower.empty() )
   if( v_secondary_spinning_reserve.size() != f_time_horizon )
    throw( std::logic_error( "SlackUnitBlock::generate_objective: "
                             "v_secondary_spinning_reserve must have size "
                             "equal to the time horizon." ) );

 auto lf = new LinearFunction();

 for( Index t = 0 ; t < f_time_horizon ; ++t ) {

  if( ! v_ActivePowerCost.empty() )
   lf->add_variable( &v_active_power[ t ] , v_ActivePowerCost[ t ] , eDryRun );
  else
   lf->add_variable( &v_active_power[ t ] , 0.0 , eDryRun );

  if( reserve_vars & 1u ) {  // if UCBlock has primary demand variables
   if( ! v_MaxPrimaryPower.empty() ) {
    if( ! v_PrimaryCost.empty() )
     lf->add_variable( &v_primary_spinning_reserve[ t ] ,
                       v_PrimaryCost[ t ] , eDryRun );
    else
     lf->add_variable( &v_primary_spinning_reserve[ t ] , 0.0 , eDryRun );
   }
  }

  if( reserve_vars & 2u ) {  // if UCBlock has secondary demand variables
   if( ! v_MaxSecondaryPower.empty() ) {
    if( ! v_SecondaryCost.empty() )
     lf->add_variable( &v_secondary_spinning_reserve[ t ] ,
                       v_SecondaryCost[ t ] , eDryRun );
    else
     lf->add_variable( &v_secondary_spinning_reserve[ t ] , 0.0 , eDryRun );
   }
  }

  if( reserve_vars & 4u ) {
   if( ( ! v_InertiaCost.empty() ) && ( ! v_MaxInertia.empty() ) )
    lf->add_variable( &v_commitment[ t ] ,
                      v_InertiaCost[ t ] * v_MaxInertia[ t ] , eDryRun );
   else
    lf->add_variable( &v_commitment[ t ] , 0.0 , eDryRun );
  }
 }

 objective.set_function( lf );
 objective.set_sense( Objective::eMin );

 // Set Block objective
 this->set_objective( &objective );

 set_objective_generated();

}  // end( SlackUnitBlock::generate_objective )

/*--------------------------------------------------------------------------*/
/*----------------- METHODS FOR CHECKING THE SlackUnitBlock ----------------*/
/*--------------------------------------------------------------------------*/

bool SlackUnitBlock::is_feasible( bool useabstract , Configuration * fsbc )
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
  && ColVariable::is_feasible( v_commitment , tol )
  && ColVariable::is_feasible( v_active_power , tol )
  && ColVariable::is_feasible( v_primary_spinning_reserve , tol )
  && ColVariable::is_feasible( v_secondary_spinning_reserve , tol )
  // Constraints: notice that the ZOConstraint are not checked, since the
  // corresponding check is made on the ColVariable
  && RowConstraint::is_feasible( ActivePower_Bound_Const , tol , rel_viol )
  && RowConstraint::is_feasible( Primary_Spinning_Reserve_Bound_Const , tol , rel_viol )
  && RowConstraint::is_feasible( Secondary_Spinning_Reserve_Bound_Const , tol , rel_viol ) );

} // end( SlackUnitBlock::is_feasible )

/*--------------------------------------------------------------------------*/
/*------- METHODS FOR LOADING, PRINTING & SAVING THE SlackUnitBlock --------*/
/*--------------------------------------------------------------------------*/

void SlackUnitBlock::serialize( netCDF::NcGroup & group ) const
{
 UnitBlock::serialize( group );

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
   throw( std::logic_error(
    "SlackUnitBlock::serialize: invalid dimension for variable " + var_name +
    ": " + std::to_string( data.size() ) + ". Its dimension must be one of " +
    "the following: TimeHorizon, NumberIntervals, 1." ) );

  ::serialize( group , var_name , ncType , dimension , data ,
               allow_scalar_var );
 };

 serialize( "MaxPower" , v_MaxPower );
 serialize( "MaxPrimaryPower" , v_MaxPrimaryPower );
 serialize( "MaxSecondaryPower" , v_MaxSecondaryPower );
 serialize( "MaxInertia" , v_MaxInertia );
 serialize( "ActivePowerCost" , v_ActivePowerCost );
 serialize( "PrimaryCost" , v_PrimaryCost );
 serialize( "SecondaryCost" , v_SecondaryCost );
 serialize( "InertiaCost" , v_InertiaCost );

}  // end( SlackUnitBlock::serialize )

/*--------------------------------------------------------------------------*/
/*----------------------- End File SlackUnitBlock.cpp ----------------------*/
/*--------------------------------------------------------------------------*/
