/*--------------------------------------------------------------------------*/
/*---------------------------- File test.cpp -------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Unit test of the UCBlock module on instances built in memory.
 *
 * The dynamic programming Solver of the module, i.e., ThermalUnitDPSolver
 * and ThermalUnitExtDPSolver on a ThermalUnitBlock and
 * NuclearUnitExtDPSolver on a NuclearUnitBlock, are compared with each other
 * and with a brute force that enumerates the commitment and runs a dynamic
 * programming over the integer power levels. With linear costs and integer
 * data the vertices of the economic dispatch of a fixed commitment are
 * integer (its rows are differences of two consecutive powers), hence the
 * brute force is exact there; with quadratic costs it only bounds the
 * optimum from above, and the optimum of a few small cases is known in
 * closed form. The instances cover the horizon of one instant, the minimum
 * up and down times that the initial state or the end of the horizon cuts,
 * the ramps with the limits of the start-up and of the shut-down, the
 * commitment fixed ON or OFF in some instants (infeasibility included), the
 * scale factor, and the operating rules of a nuclear unit that modulates.
 *
 * Then every UnitBlock and NetworkBlock of the module, and a UCBlock with
 * no unit or with one unit of each kind on two nodes, is written in a
 * netCDF group, read back, written again and read back again: the two
 * groups written by the Block have to coincide, and the data read have to
 * be those written, with the degenerate cases of a single instant and of
 * the optional data left out.
 *
 * Finally the setters of the data of a ThermalUnitBlock and of a UCBlock
 * are called with a FakeSolver attached, which has to receive the
 * physical Modification, while the abstract representation has to follow
 * the change, the scale factor included; the DP Solver attached to the
 * unit has to follow it, too.
 *
 * The test needs nothing but the core SMS++, so that the CI of UCBlock
 * builds this module alone. The netCDF groups are written in a dataset in
 * memory, no file being written.
 *
 * \author Donato Meoli \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>
#include <memory>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include <netcdf>

#include "BatteryUnitBlock.h"
#include "DCNetworkBlock.h"
#include "ACNetworkBlock.h"
#include "ECNetworkBlock.h"
#include "FakeSolver.h"
#include "FRealObjective.h"
#include "DQuadFunction.h"
#include "HydroUnitBlock.h"
#include "IntermittentUnitBlock.h"
#include "NuclearUnitBlock.h"
#include "NuclearUnitExtDPSolver.h"
#include "SlackUnitBlock.h"
#include "ThermalUnitBlock.h"
#include "ThermalUnitDPSolver.h"
#include "ThermalUnitExtDPSolver.h"
#include "UCBlock.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

using Index = Block::Index;
using Range = Block::Range;
using Subset = Block::Subset;
using sp_Mod = std::shared_ptr< Modification >;

/*--------------------------------------------------------------------------*/
/*-------------------------------- GLOBALS ---------------------------------*/
/*--------------------------------------------------------------------------*/

static int failures = 0;  // number of failed checks

static const double INF = Inf< double >();

static int NCID = -1;  // the netCDF dataset in memory the groups go in

static unsigned int n_groups = 0;  // groups written so far

static std::mt19937 rg( 20260926 );  // fixed seed

static const std::vector< std::string > THERMAL_DP = {
 "ThermalUnitDPSolver" , "ThermalUnitExtDPSolver" };

/*--------------------------------------------------------------------------*/
/*------------------------------ FUNCTIONS ---------------------------------*/
/*--------------------------------------------------------------------------*/

static void check( bool ok , const std::string & what )
{
 if( ! ok ) {
  ++failures;
  std::cout << "FAILED: " << what << std::endl;
  }
 }

/*--------------------------------------------------------------------------*/

static bool close( double a , double b , double eps = 1e-6 )
{
 if( ( a == INF ) || ( b == INF ) )
  return( a == b );
 return( std::abs( a - b ) <= eps * std::max( 1.0 , std::abs( b ) ) );
 }

/*--------------------------------------------------------------------------*/

static std::string str( double v )
{
 std::ostringstream s;
 s << v;
 return( s.str() );
 }

/*--------------------------------------------------------------------------*/
/// a new group, with a name never used before, in the netCDF dataset in
/// memory; if fresh, the dataset is a new one, the groups of the previous
/// one being no longer usable (a dataset that holds all the groups of the
/// test gets slower and slower to write)

static netCDF::NcGroup new_group( const std::string & prefix ,
				  bool fresh = false )
{
 if( fresh || ( NCID < 0 ) ) {
  if( NCID >= 0 )
   nc_close( NCID );
  if( nc_create( "unit_test.nc4" , NC_DISKLESS | NC_NETCDF4 | NC_CLOBBER ,
		 & NCID ) != NC_NOERR )
   throw( std::runtime_error( "new_group: cannot create a netCDF dataset "
			      "in memory" ) );
  }
 return( netCDF::NcGroup( NCID ).addGroup( prefix + "_" +
					   std::to_string( n_groups++ ) ) );
 }

/*--------------------------------------------------------------------------*/

static int rnd( int lo , int hi )
{
 return( std::uniform_int_distribution< int >( lo , hi )( rg ) );
 }

/*--------------------------------------------------------------------------*/
/*------------------------ WRITING THE NETCDF GROUPS -----------------------*/
/*--------------------------------------------------------------------------*/

static void put( netCDF::NcGroup & g , const std::string & n , double v )
{
 g.addVar( n , netCDF::NcDouble() ).putVar( & v );
 }

static void put_int( netCDF::NcGroup & g , const std::string & n , int v )
{
 g.addVar( n , netCDF::NcInt() ).putVar( & v );
 }

static void put_uint( netCDF::NcGroup & g , const std::string & n ,
		      unsigned int v )
{
 g.addVar( n , netCDF::NcUint() ).putVar( & v );
 }

static void put( netCDF::NcGroup & g , const std::string & n ,
		 const netCDF::NcDim & d , const std::vector< double > & v )
{
 g.addVar( n , netCDF::NcDouble() , d ).putVar( v.data() );
 }

static void put( netCDF::NcGroup & g , const std::string & n ,
		 const std::vector< netCDF::NcDim > & d ,
		 const std::vector< double > & v )
{
 g.addVar( n , netCDF::NcDouble() , d ).putVar( v.data() );
 }

static void put_int( netCDF::NcGroup & g , const std::string & n ,
		     const netCDF::NcDim & d , const std::vector< int > & v )
{
 g.addVar( n , netCDF::NcInt() , d ).putVar( v.data() );
 }

/*--------------------------------------------------------------------------*/
/*-------------------- THE DATA OF A (NUCLEAR) THERMAL UNIT ----------------*/
/*--------------------------------------------------------------------------*/
/* The data of a ThermalUnitBlock, or of a NuclearUnitBlock, as the test
 * writes it. Bounds, ramps and limits are the same at every instant, since
 * the ramps and the limits are indexed differently by the two sides of a
 * step and the test only relies on what the documentation says; the costs
 * change with the instant. An empty cost vector is not written. */

struct TUData {
 Index T = 1;                 // time horizon
 double minP = 1;             // minimum power
 double maxP = 10;            // maximum power
 double ru = 10;              // ramp-up
 double rd = 10;              // ramp-down
 double su = 10;              // start-up limit
 double sd = 10;              // shut-down limit
 double initP = 0;            // initial power, if on before the horizon
 int initUD = -1;             // initial up (> 0) or down (<= 0) time
 unsigned int minUp = 1;      // minimum up time
 unsigned int minDown = 1;    // minimum down time
 std::vector< double > lin;   // linear term
 std::vector< double > quad;  // quadratic term
 std::vector< double > cnst;  // constant term
 std::vector< double > suc;   // start-up cost
 std::vector< double > sdc;   // shut-down cost
 double scale = 1;            // scale factor
 // the operating rules of a nuclear unit, used only if nuclear
 bool nuclear = false;
 unsigned int modT = 2;       // ModulationTime
 unsigned int initMod = 2;    // InitModulation
 double mru = 0;              // ModulationDeltaRampUp
 double mrd = 0;              // ModulationDeltaRampDown
 std::vector< double > bands; // PowerBands, the two breakpoints, if any
 };

/*--------------------------------------------------------------------------*/

static double at( const std::vector< double > & v , Index t )
{
 return( v.empty() ? 0 : v[ t ] );
 }

/*--------------------------------------------------------------------------*/

static void write_TU( netCDF::NcGroup g , const TUData & d )
{
 g.putAtt( "type" , d.nuclear ? "NuclearUnitBlock" : "ThermalUnitBlock" );
 g.addDim( "TimeHorizon" , d.T );
 auto NI = g.addDim( "NumberIntervals" , d.T );
 auto cst = [ & ]( double v ) { return( std::vector< double >( d.T , v ) ); };

 put( g , "MinPower" , NI , cst( d.minP ) );
 put( g , "MaxPower" , NI , cst( d.maxP ) );
 put( g , "DeltaRampUp" , NI , cst( d.ru ) );
 put( g , "DeltaRampDown" , NI , cst( d.rd ) );
 put( g , "StartUpLimit" , NI , cst( d.su ) );
 put( g , "ShutDownLimit" , NI , cst( d.sd ) );
 if( ! d.lin.empty() )
  put( g , "LinearTerm" , NI , d.lin );
 if( ! d.quad.empty() )
  put( g , "QuadTerm" , NI , d.quad );
 if( ! d.cnst.empty() )
  put( g , "ConstTerm" , NI , d.cnst );
 if( ! d.suc.empty() )
  put( g , "StartUpCost" , NI , d.suc );
 if( ! d.sdc.empty() )
  put( g , "ShutDownCost" , NI , d.sdc );
 if( d.initUD > 0 )
  put( g , "InitialPower" , d.initP );
 put_int( g , "InitUpDownTime" , d.initUD );
 put_uint( g , "MinUpTime" , d.minUp );
 put_uint( g , "MinDownTime" , d.minDown );
 if( d.scale != 1 )
  put( g , "Scale" , d.scale );

 if( d.nuclear ) {
  put_uint( g , "ModulationTime" , d.modT );
  put_uint( g , "InitModulation" , d.initMod );
  put( g , "ModulationDeltaRampUp" , NI , cst( d.mru ) );
  put( g , "ModulationDeltaRampDown" , NI , cst( d.mrd ) );
  if( ! d.bands.empty() )
   put( g , "PowerBands" , g.addDim( "NumberPowerBands" , 2 ) , d.bands );
  }
 }

/*--------------------------------------------------------------------------*/
/// a ThermalUnitBlock (or NuclearUnitBlock) with the data d, read from netCDF

static ThermalUnitBlock * new_TU( const TUData & d )
{
 auto g = new_group( "TU" , true );
 write_TU( g , d );
 auto b = Block::new_Block( g );
 auto tub = dynamic_cast< ThermalUnitBlock * >( b );
 if( ! tub )
  throw( std::logic_error( "new_TU: the ThermalUnitBlock was not built" ) );
 return( tub );
 }

/*--------------------------------------------------------------------------*/
/*------------------------------ BRUTE FORCE -------------------------------*/
/*--------------------------------------------------------------------------*/
/* The optimum of the unit over the integer power levels: every commitment
 * that respects the minimum up and down times (with the state before the
 * horizon, and a run that the end of the horizon cuts being allowed) and
 * the fixings in fix[] (-1 free, 0 OFF, 1 ON) is enumerated, and for each
 * a dynamic programming over the power (and over the lockout of the
 * modulation of a nuclear unit, see NuclearUnitExtDPSolver.h) finds the
 * cheapest dispatch. The steps are
 *
 * - ON -> ON: the power moves within the ramps; for a nuclear unit within
 *   the modulation ramps (and the ramps), or, if the lockout is 0, within
 *   the ramps by a modulation that sets the lockout to ModulationTime - 1;
 *
 * - OFF -> ON: the power is at most the start-up limit;
 *
 * - ON -> OFF: the power before is at most the shut-down limit;
 *
 * and each of them but a modulation decreases the lockout by one. Returns
 * +INF if nothing is feasible. Exact for linear costs and integer data. */

static double brute_force( const TUData & d , const std::vector< int > & fix )
{
 const int P = int( d.maxP );
 const int L = d.nuclear ? int( d.modT ) : 1;
 const int l0 = d.nuclear ? std::max( int( d.modT ) - int( d.initMod ) , 0 )
                          : 0;
 const bool on0 = d.initUD > 0;
 const int run0 = std::abs( d.initUD );
 const int up = std::max( 1u , d.minUp );
 const int down = std::max( 1u , d.minDown );
 const int nru = d.nuclear ? int( std::min( d.ru , d.mru ) ) : int( d.ru );
 const int nrd = d.nuclear ? int( std::min( d.rd , d.mrd ) ) : int( d.rd );

 auto cost_on = [ & ]( Index t , int p ) {
  return( at( d.cnst , t ) + at( d.lin , t ) * p +
	  at( d.quad , t ) * p * p );
  };

 double best = INF;
 const Index T = d.T;
 std::vector< double > f( ( P + 1 ) * L ) , g( ( P + 1 ) * L );
 auto id = [ L ]( int p , int l ) { return( p * L + l ); };

 for( unsigned long mask = 0 ; mask < ( 1ul << T ) ; ++mask ) {
  // the commitment, the fixings and the minimum up and down times
  bool ok = true;
  bool prev = on0;
  int run = run0;
  for( Index t = 0 ; ok && ( t < T ) ; ++t ) {
   const bool u = ( mask >> t ) & 1;
   if( ( fix[ t ] >= 0 ) && ( fix[ t ] != int( u ) ) )
    ok = false;
   else
    if( u != prev ) {
     if( prev ? ( run < up ) : ( run < down ) )
      ok = false;
     prev = u;
     run = 1;
     }
    else
     ++run;
   }
  if( ! ok )
   continue;

  // the dispatch
  std::fill( f.begin() , f.end() , INF );
  f[ id( on0 ? int( d.initP ) : 0 , l0 ) ] = 0;
  prev = on0;
  for( Index t = 0 ; t < T ; ++t ) {
   const bool u = ( mask >> t ) & 1;
   std::fill( g.begin() , g.end() , INF );
   for( int q = 0 ; q <= P ; ++q )
    for( int l = 0 ; l < L ; ++l ) {
     const double fq = f[ id( q , l ) ];
     if( fq == INF )
      continue;
     const int ld = std::max( l - 1 , 0 );
     auto relax = [ & ]( int p , int nl , double c ) {
      if( fq + c < g[ id( p , nl ) ] )
       g[ id( p , nl ) ] = fq + c;
      };
     if( ! u ) {
      if( ! prev )
       relax( 0 , ld , 0 );
      else
       if( q <= d.sd )
	relax( 0 , ld , at( d.sdc , t ) );
      continue;
      }
     for( int p = int( d.minP ) ; p <= P ; ++p ) {
      if( ! prev ) {
       if( p <= d.su )
	relax( p , ld , cost_on( t , p ) + at( d.suc , t ) );
       continue;
       }
      if( ( p - q <= nru ) && ( q - p <= nrd ) )
       relax( p , ld , cost_on( t , p ) );
      else
       if( d.nuclear && ( l == 0 ) && ( p - q <= d.ru ) && ( q - p <= d.rd ) )
	relax( p , L - 1 , cost_on( t , p ) );
      }
     }
   std::swap( f , g );
   prev = u;
   }
  for( auto v : f )
   best = std::min( best , v );
  }

 return( best * d.scale );
 }

/*--------------------------------------------------------------------------*/
/*------------------------- RUNNING THE DP SOLVERS -------------------------*/
/*--------------------------------------------------------------------------*/
/* Solves the unit with the given Solver, attached for the occasion, and
 * checks that it gives the expected value (+INF = infeasible), that its
 * lower and upper bounds coincide and that the schedule of its Solution is
 * feasible for the unit. If the unit has its abstract representation, the
 * schedule is also written in the Variable, where it has to satisfy every
 * Constraint and to cost what the Solver says. */

static void check_DP( ThermalUnitBlock * tub , const std::string & sname ,
		      double expected , const std::string & what )
{
 const auto who = what + ", " + sname;
 auto slv = Solver::new_Solver( sname );
 if( ! slv ) {
  check( false , who + ": the Solver is not in the factory" );
  return;
  }
 tub->register_Solver( slv );

 int status;
 try {
  status = slv->compute();
  }
 catch( std::exception & e ) {
  check( false , who + ": compute() throws " + e.what() );
  tub->unregister_Solver( slv , true );
  return;
  }

 if( expected == INF ) {
  check( status == Solver::kInfeasible ,
	 who + ": infeasible expected, status " + std::to_string( status ) );
  tub->unregister_Solver( slv , true );
  return;
  }

 check( status == Solver::kOK , who + ": status " + std::to_string( status ) );
 const double ub = slv->get_ub();
 check( close( ub , expected ) ,
	who + ": value " + str( ub ) + " instead of " + str( expected ) );
 check( close( slv->get_lb() , ub ) , who + ": lb differs from ub" );

 if( slv->has_var_solution() ) {
  SimpleConfiguration< double > tol( 1e-7 );
  // a unit whose schedule does not answer for it is checked through the
  // Variable, which have to be there
  if( tub->is_sol_feasible_physical() || tub->get_commitment( 0 ) ) {
   auto sol = slv->get_Solution();
   check( sol && tub->is_sol_feasible( sol , & tol ) ,
	  who + ": the schedule of the Solution is not feasible" );
   delete sol;
   }

  if( tub->get_commitment( 0 ) && tub->get_objective() ) {
   slv->get_var_solution();
   check( tub->is_feasible( true , & tol ) ,
	  who + ": the schedule in the Variable is not feasible" );
   auto obj = static_cast< FRealObjective * >( tub->get_objective() );
   obj->compute();
   check( close( obj->value() , ub ) ,
	  who + ": the Objective values the schedule " +
	  str( obj->value() ) + " instead of " + str( ub ) );
   }
  }

 tub->unregister_Solver( slv , true );
 }

/*--------------------------------------------------------------------------*/
/// every DP Solver of the unit against the expected value

static void check_all_DP( ThermalUnitBlock * tub , double expected ,
			  const std::string & what )
{
 if( dynamic_cast< NuclearUnitBlock * >( tub ) )
  check_DP( tub , "NuclearUnitExtDPSolver" , expected , what );
 else
  for( const auto & s : THERMAL_DP )
   check_DP( tub , s , expected , what );
 }

/*--------------------------------------------------------------------------*/
/// the unit with its abstract representation, the default T formulation

static void generate_all( Block * b )
{
 b->generate_abstract_variables();
 b->generate_abstract_constraints();
 b->generate_objective();
 }

/*--------------------------------------------------------------------------*/
/*------------------------------ DP: KNOWN CASES ---------------------------*/
/*--------------------------------------------------------------------------*/
/// a horizon of one instant, linear and quadratic cost

static void test_DP_one_instant( void )
{
 TUData d;  // off since long, starts if it pays
 d.initUD = -3;
 d.lin = { -2 };
 d.cnst = { 3 };
 d.suc = { 4 };
 // on at p = 10: 4 + 3 - 20 = -13
 auto tub = new_TU( d );
 check_all_DP( tub , -13 , "one instant, linear" );
 delete tub;

 d.quad = { 1 };
 d.lin = { -8 };
 // on at p = 4: 4 + 3 + 16 - 32 = -9
 tub = new_TU( d );
 check_all_DP( tub , -9 , "one instant, quadratic" );
 delete tub;

 d.cnst = { 20 };  // now starting does not pay
 tub = new_TU( d );
 check_all_DP( tub , 0 , "one instant, stays off" );
 delete tub;

 d.initUD = 1;  // on before, the minimum up time keeps it on
 d.initP = 4;
 d.minUp = 3;
 d.quad.clear();
 d.lin = { 5 };
 d.cnst = { 1 };
 d.rd = 1;
 // p >= 3 by the ramp-down: 1 + 15
 tub = new_TU( d );
 check_all_DP( tub , 16 , "one instant, kept on by the minimum up time" );
 delete tub;
 }

/*--------------------------------------------------------------------------*/
/// ramps binding a quadratic cost, and the scale factor

static void test_DP_ramps( void )
{
 TUData d;  // on since long at 2, can move by 1 per instant
 d.T = 2;
 d.initUD = 5;
 d.initP = 2;
 d.ru = d.rd = 1;
 d.quad = { 1 , 1 };
 d.lin = { -12 , -12 };
 // p = 3 , 4: ( 9 - 36 ) + ( 16 - 48 ) = -59, better than shutting down at
 // 0 and restarting at 1 at p = 6 (-36)
 auto tub = new_TU( d );
 check_all_DP( tub , -59 , "ramps" );
 delete tub;

 d.scale = 3;
 tub = new_TU( d );
 check_all_DP( tub , -177 , "ramps, scale 3" );
 delete tub;

 d.scale = 0;
 tub = new_TU( d );
 check_all_DP( tub , 0 , "ramps, scale 0" );
 delete tub;
 }

/*--------------------------------------------------------------------------*/
/// the minimum down time cut by the initial state, the minimum up time cut
/// by the end of the horizon

static void test_DP_boundary( void )
{
 TUData d;  // off since 1 instant, has to stay off for 2 more
 d.T = 3;
 d.maxP = d.su = d.sd = d.ru = d.rd = 5;
 d.initUD = -1;
 d.minDown = 3;
 d.minUp = 3;
 d.lin = { -10 , -10 , -10 };
 // on at t = 2 only, which the minimum up time allows at the end: -50
 auto tub = new_TU( d );
 check_all_DP( tub , -50 , "minimum down time from the initial state" );
 delete tub;

 d.initUD = 0;  // shut down right at the end of -1: off at all the 3
 tub = new_TU( d );
 check_all_DP( tub , 0 , "just shut down" );
 delete tub;

 d.initUD = -3;  // free at once
 tub = new_TU( d );
 check_all_DP( tub , -150 , "free at once" );
 delete tub;

 d.initUD = 5;  // on at 5, cannot shut down at 0 above the limit
 d.initP = 5;
 d.sd = 2;
 d.rd = 1;
 d.lin = { 10 , 10 , 10 };
 // shutting down at t needs the power at t - 1 within the limit 2, which
 // the ramp-down reaches only at 2: on at 4 , 3 , 2, i.e. 90
 tub = new_TU( d );
 check_all_DP( tub , 90 , "shut-down trajectory" );
 delete tub;
 }

/*--------------------------------------------------------------------------*/
/*---------------------------- DP: RANDOM CASES ----------------------------*/
/*--------------------------------------------------------------------------*/

static TUData random_TU( bool quadratic )
{
 TUData d;
 d.T = rnd( 1 , 6 );
 d.minP = rnd( 1 , 3 );
 d.maxP = d.minP + rnd( 1 , 8 );
 d.ru = rnd( 1 , int( d.maxP ) );
 d.rd = rnd( 1 , int( d.maxP ) );
 d.su = rnd( int( d.minP ) , int( d.maxP ) );
 d.sd = rnd( int( d.minP ) , int( d.maxP ) );
 d.initUD = rnd( -6 , 6 );
 if( d.initUD > 0 )
  d.initP = rnd( int( d.minP ) , int( d.maxP ) );
 // also longer than the horizon, see test_DP_clamped_min_up_down()
 d.minUp = rnd( 1 , 9 );
 d.minDown = rnd( 1 , 9 );
 d.lin.resize( d.T );
 d.cnst.resize( d.T );
 d.suc.resize( d.T );
 for( Index t = 0 ; t < d.T ; ++t ) {
  d.lin[ t ] = rnd( -6 , 3 );
  d.cnst[ t ] = rnd( 0 , 6 );
  d.suc[ t ] = rnd( 0 , 8 );
  }
 if( rnd( 0 , 1 ) ) {
  d.sdc.resize( d.T );
  for( auto & c : d.sdc )
   c = rnd( 0 , 3 );
  }
 if( quadratic ) {
  d.quad.resize( d.T );
  for( auto & q : d.quad )
   q = rnd( 0 , 4 ) * 0.25;
  }
 return( d );
 }

/*--------------------------------------------------------------------------*/

static std::string describe( const TUData & d )
{
 std::ostringstream s;
 s << "T " << d.T << " P [" << d.minP << "," << d.maxP << "] ramps "
   << d.ru << "/" << d.rd << " limits " << d.su << "/" << d.sd
   << " init " << d.initUD << "@" << d.initP << " up/down " << d.minUp
   << "/" << d.minDown;
 if( d.nuclear )
  s << " modT " << d.modT << " initMod " << d.initMod << " mramps "
    << d.mru << "/" << d.mrd;
 return( s.str() );
 }

/*--------------------------------------------------------------------------*/
/// the two thermal DP Solvers against the brute force, linear costs

static void test_DP_random_linear( void )
{
 for( int i = 0 ; i < 150 ; ++i ) {
  auto d = random_TU( false );
  auto tub = new_TU( d );
  const auto bf = brute_force( d , std::vector< int >( d.T , -1 ) );
  check_all_DP( tub , bf , "random linear " + std::to_string( i ) + " (" +
		describe( d ) + ")" );
  delete tub;
  }
 }

/*--------------------------------------------------------------------------*/
/// the two thermal DP Solvers against each other, quadratic costs

static void test_DP_random_quadratic( void )
{
 for( int i = 0 ; i < 100 ; ++i ) {
  auto d = random_TU( true );
  auto tub = new_TU( d );
  const auto what = "random quadratic " + std::to_string( i ) + " (" +
                    describe( d ) + ")";
  std::vector< double > val;
  for( const auto & s : THERMAL_DP ) {
   auto slv = Solver::new_Solver( s );
   tub->register_Solver( slv );
   const int status = slv->compute();
   val.push_back( status == Solver::kOK ? slv->get_ub() : INF );
   tub->unregister_Solver( slv , true );
   }
  check( close( val[ 0 ] , val[ 1 ] ) ,
	 what + ": " + str( val[ 0 ] ) + " vs " + str( val[ 1 ] ) );
  // the integer powers are a restriction
  const auto bf = brute_force( d , std::vector< int >( d.T , -1 ) );
  check( val[ 1 ] <= bf + 1e-6 ,
	 what + ": " + str( val[ 1 ] ) + " above the integer " + str( bf ) );
  delete tub;
  }
 }

/*--------------------------------------------------------------------------*/
/// the commitment fixed ON or OFF in some instants

static void test_DP_fixed( void )
{
 // a known case first: on every instant pays 50, OFF fixed at 1
 TUData d;
 d.T = 3;
 d.maxP = d.su = d.sd = d.ru = d.rd = 5;
 d.initUD = -10;
 d.lin = { -10 , -10 , -10 };

 auto run = [ & ]( const std::vector< int > & fix , double expected ,
		   const std::string & what ) {
  auto tub = new_TU( d );
  generate_all( tub );
  auto u = tub->get_commitment( 0 );
  for( Index t = 0 ; t < d.T ; ++t )
   if( ( fix[ t ] >= 0 ) && ( ! u[ t ].is_fixed() ) ) {
    u[ t ].set_value( fix[ t ] );
    u[ t ].is_fixed( true );
    }
  check( close( brute_force( d , fix ) , expected ) ,
	 what + ": the brute force gives " + str( brute_force( d , fix ) ) );
  check_all_DP( tub , expected , what );
  delete tub;
  };

 run( { -1 , -1 , -1 } , -150 , "fixed, none" );
 run( { -1 , 0 , -1 } , -100 , "fixed, OFF in the middle" );
 d.minDown = 2;
 run( { 1 , 0 , -1 } , -50 , "fixed, ON then OFF, minimum down time 2" );
 d.minDown = 1;
 d.minUp = 2;
 run( { -1 , 0 , -1 } , -50 , "fixed, OFF, minimum up time 2" );
 d.minUp = 1;
 d.lin = { 10 , 10 , 10 };
 run( { -1 , 1 , -1 } , 10 , "fixed ON where it does not pay" );

 // on since long above the shut-down limit: OFF at 0 is infeasible
 d.initUD = 5;
 d.initP = 5;
 d.sd = 2;
 d.rd = 1;
 run( { 0 , -1 , -1 } , INF , "fixed OFF, infeasible" );

 // random fixings against the brute force
 for( int i = 0 ; i < 60 ; ++i ) {
  auto r = random_TU( false );
  auto tub = new_TU( r );
  generate_all( tub );
  auto u = tub->get_commitment( 0 );
  std::vector< int > fix( r.T , -1 );
  for( Index t = 0 ; t < r.T ; ++t )
   if( u[ t ].is_fixed() )
    fix[ t ] = u[ t ].get_value() >= 0.5;  // fixed by the initial state
   else
    if( rnd( 0 , 2 ) == 0 ) {
     fix[ t ] = rnd( 0 , 1 );
     u[ t ].set_value( fix[ t ] );
     u[ t ].is_fixed( true );
     }
  check_all_DP( tub , brute_force( r , fix ) , "random fixed " +
		std::to_string( i ) + " (" + describe( r ) + ")" );
  delete tub;
  }
 }

/*--------------------------------------------------------------------------*/
/// the minimum up and down times set before the Variable are generated

static void test_DP_set_min_up_down( void )
{
 TUData d;
 d.T = 4;
 d.maxP = d.su = d.sd = d.ru = d.rd = 5;
 d.initUD = -10;
 d.lin = { -10 , 20 , -10 , -10 };  // off at 1 unless the up time forbids
 auto tub = new_TU( d );
 check_all_DP( tub , -150 , "before set_min_up_down_time" );

 // with a minimum up time of 2 it stays on at 1 at the minimum power
 tub->set_min_up_down_time( 2 , 1 );
 d.minUp = 2;
 check( close( brute_force( d , std::vector< int >( d.T , -1 ) ) , -130 ) ,
	"set_min_up_down_time: the brute force" );
 check( tub->get_min_up_time() == 2 , "set_min_up_down_time: up time" );
 check_all_DP( tub , -130 , "after set_min_up_down_time" );

 tub->generate_abstract_variables();
 bool thrown = false;
 try {
  tub->set_min_up_down_time( 1 , 1 );
  }
 catch( std::exception & ) {
  thrown = true;
  }
 check( thrown , "set_min_up_down_time does not throw after the Variable" );
 delete tub;
 }

/*--------------------------------------------------------------------------*/
/* A minimum up (down) time longer than the horizon is bounded by the
 * horizon plus the number of instants the unit has been on (off) before it
 * [see ThermalUnitBlock::deserialize()], which keeps the unit as it is for
 * the whole horizon if the minimum time says so: here the minimum time is 4,
 * the horizon 1 and the unit has been so for 3 instants, hence it has to
 * stay as it is at 0. Then the same through set_min_up_down_time(). */

static void test_DP_clamped_min_up_down( void )
{
 TUData d;
 d.maxP = d.su = d.sd = d.ru = d.rd = 5;
 d.initUD = -3;  // off since 3, starting at 0 would pay 50
 d.minDown = 4;
 d.lin = { -10 };
 check( brute_force( d , { -1 } ) == 0 ,
	"clamped down time: the brute force" );
 auto tub = new_TU( d );
 check( tub->get_min_down_time() == 4 ,
	"clamped down time: MinDownTime 4 is taken to be " +
	std::to_string( tub->get_min_down_time() ) );
 check_all_DP( tub , 0 , "clamped down time, off since 3 of 4" );
 tub->set_min_up_down_time( 1 , 9 );
 check( tub->get_min_down_time() == 4 ,
	"set_min_up_down_time: MinDownTime 9 is taken to be " +
	std::to_string( tub->get_min_down_time() ) );
 check_all_DP( tub , 0 , "clamped down time, set_min_up_down_time" );
 delete tub;

 d.initUD = 3;  // on since 3 at the minimum, staying on at 0 costs 10
 d.initP = 1;
 d.minDown = 1;
 d.minUp = 4;
 d.lin = { 10 };
 check( brute_force( d , { -1 } ) == 10 , "clamped up time: the brute force" );
 tub = new_TU( d );
 check( tub->get_min_up_time() == 4 ,
	"clamped up time: MinUpTime 4 is taken to be " +
	std::to_string( tub->get_min_up_time() ) );
 check_all_DP( tub , 10 , "clamped up time, on since 3 of 4" );
 delete tub;

 d.minUp = 9;  // larger than any bound: the same unit
 tub = new_TU( d );
 check( tub->get_min_up_time() == 4 ,
	"clamped up time: MinUpTime 9 is taken to be " +
	std::to_string( tub->get_min_up_time() ) );
 check_all_DP( tub , 10 , "clamped up time 9, on since 3" );
 delete tub;
 }

/*--------------------------------------------------------------------------*/
/*------------------------------- DP: NUCLEAR ------------------------------*/
/*--------------------------------------------------------------------------*/
/// NuclearUnitExtDPSolver against the brute force, and against the thermal
/// DP Solvers when the modulation ramps are the ramps

static void test_DP_nuclear( void )
{
 // a known case: on since long at 2, the price asks to go up by 3 at 0
 // and down by 3 at 1, while stable moves are 1 and a modulation is 3
 TUData d;
 d.nuclear = true;
 d.T = 2;
 d.minP = 2;
 d.maxP = d.su = d.sd = 5;
 d.ru = d.rd = 3;
 d.mru = d.mrd = 1;
 d.initUD = 10;
 d.initP = 2;
 d.modT = 3;
 d.initMod = 3;
 d.lin = { -10 , 10 };
 d.sdc = { 100 , 100 };  // shutting down does not pay
 // one modulation only in 3 instants: up to 5 at 0 (-50), then 4 at 1
 // (40), i.e. -10; or 3 at 0 (-30) and 2 at 1 (20), i.e. -10 as well,
 // while a thermal unit with the same ramps goes 5 , 2, i.e. -30
 auto tub = new_TU( d );
 check( close( brute_force( d , { -1 , -1 } ) , -10 ) ,
	"nuclear known: the brute force gives " +
	str( brute_force( d , { -1 , -1 } ) ) );
 check_all_DP( tub , -10 , "nuclear, one modulation" );
 delete tub;

 auto td = d;
 td.nuclear = false;
 tub = new_TU( td );
 check_all_DP( tub , -30 , "nuclear, the thermal twin" );
 delete tub;

 d.modT = 2;  // two consecutive modulations are forbidden still
 d.initMod = 2;
 tub = new_TU( d );
 check_all_DP( tub , brute_force( d , { -1 , -1 } ) , "nuclear, modT 2" );
 delete tub;

 d.initMod = 1;  // a modulation right before the horizon
 tub = new_TU( d );
 check_all_DP( tub , brute_force( d , { -1 , -1 } ) ,
	       "nuclear, locked at 0" );
 delete tub;

 for( int i = 0 ; i < 80 ; ++i ) {
  auto r = random_TU( false );
  r.nuclear = true;
  r.modT = rnd( 2 , 4 );
  r.initMod = rnd( 1 , int( r.modT ) );
  r.mru = rnd( 0 , int( r.ru ) );
  r.mrd = rnd( 0 , int( r.rd ) );
  const auto what = "random nuclear " + std::to_string( i ) + " (" +
                    describe( r ) + ")";
  tub = new_TU( r );
  check_all_DP( tub , brute_force( r , std::vector< int >( r.T , -1 ) ) ,
		what );
  delete tub;

  // with the full ramps the modulation says nothing
  r.mru = r.ru;
  r.mrd = r.rd;
  tub = new_TU( r );
  auto tr = r;
  tr.nuclear = false;
  const auto bf = brute_force( tr , std::vector< int >( r.T , -1 ) );
  check_all_DP( tub , bf , what + ", full modulation ramps" );
  delete tub;
  }
 }

/*--------------------------------------------------------------------------*/
/*------------------- THE CHECK OF A SCHEDULE IN A Solution ----------------*/
/*--------------------------------------------------------------------------*/
/// a schedule against the fixed Variable of the unit, and a nuclear one

static void test_sol_feasible( void )
{
 SimpleConfiguration< double > tol( 1e-7 );

 // off since long, then on at 1 and 2: the start-up is at 1
 TUData d;
 d.T = 3;
 d.maxP = d.su = d.sd = d.ru = d.rd = 5;
 d.initUD = -10;
 d.lin = { -10 , -10 , -10 };
 auto tub = new_TU( d );
 generate_all( tub );
 auto u = tub->get_commitment( 0 );
 auto p = tub->get_active_power( 0 );
 auto su = tub->get_start_up();
 const std::vector< double > U = { 0 , 1 , 1 };
 const std::vector< double > P = { 0 , 5 , 5 };
 for( Index t = 0 ; t < d.T ; ++t ) {
  u[ t ].set_value( U[ t ] );
  p[ t ].set_value( P[ t ] );
  }
 check( tub->get_number_start_up() == d.T ,
	"schedule: one start-up Variable per instant" );
 for( Index t = 0 ; t < tub->get_number_start_up() ; ++t )
  su[ t ].set_value( t == 1 ? 1 : 0 );

 auto sol = tub->get_Solution( nullptr , false );
 check( tub->is_sol_feasible_physical() ,
	"schedule: a thermal unit reads the Solution" );
 check( tub->is_sol_feasible( sol , & tol ) , "schedule: nothing fixed" );

 su[ 2 ].set_value( 1 );
 su[ 2 ].is_fixed( true );
 check( ! tub->is_sol_feasible( sol , & tol ) ,
	"schedule: a start-up fixed where the schedule has none" );
 su[ 2 ].is_fixed( false );
 su[ 2 ].set_value( 0 );
 su[ 1 ].is_fixed( true );
 check( tub->is_sol_feasible( sol , & tol ) ,
	"schedule: a start-up fixed where the schedule has it" );
 su[ 1 ].is_fixed( false );
 delete sol;
 delete tub;

 // on at 2, then 5 and 2: the ramps of a thermal unit allow it, while a
 // nuclear one modulates twice within its ModulationTime
 TUData n;
 n.nuclear = true;
 n.T = 2;
 n.minP = 2;
 n.maxP = n.su = n.sd = 5;
 n.ru = n.rd = 3;
 n.mru = n.mrd = 1;
 n.initUD = 10;
 n.initP = 2;
 n.modT = 3;
 n.initMod = 3;
 n.lin = { -10 , 10 };

 for( bool nuclear : { false , true } ) {
  n.nuclear = nuclear;
  const std::string who = nuclear ? "nuclear" : "the thermal twin";
  tub = new_TU( n );
  generate_all( tub );
  u = tub->get_commitment( 0 );
  p = tub->get_active_power( 0 );
  u[ 0 ].set_value( 1 );
  u[ 1 ].set_value( 1 );
  p[ 0 ].set_value( 5 );
  p[ 1 ].set_value( 2 );
  sol = tub->get_Solution( nullptr , false );
  check( tub->is_sol_feasible_physical() != nuclear ,
	 who + ": the schedule answers for a thermal unit only" );
  check( tub->is_sol_feasible( sol , & tol ) != nuclear ,
	 who + ": up by 3 and down by 3" );
  delete sol;
  delete tub;
  }

 // the nuclear DP does not read the bands of the output: a band that is
 // fixed is refused, rather than left out of the schedule it gives
 n.nuclear = true;
 n.bands = { 3 , 4 };
 for( bool fix : { false , true } ) {
  tub = new_TU( n );
  generate_all( tub );
  auto nub = static_cast< NuclearUnitBlock * >( tub );
  check( nub->get_band() , "bands: the band Variable are there" );
  if( fix && nub->get_band() ) {
   nub->get_band()[ 1 ].set_value( 1 );
   nub->get_band()[ 1 ].is_fixed( true );
   }
  auto slv = Solver::new_Solver( "NuclearUnitExtDPSolver" );
  tub->register_Solver( slv );
  bool refused = false;
  int status = Solver::kError;
  try { status = slv->compute(); }
  catch( std::logic_error & ) { refused = true; }
  if( fix )
   check( refused , "bands: a fixed band is refused by the nuclear DP" );
  else
   check( ( ! refused ) && ( status == Solver::kOK ) ,
	  "bands: with no band fixed the nuclear DP solves, status " +
	  std::to_string( status ) );
  tub->unregister_Solver( slv , true );
  delete tub;
  }
 }

/*--------------------------------------------------------------------------*/
/*----------------------------- ROUND TRIPS --------------------------------*/
/*--------------------------------------------------------------------------*/
/* Compares two netCDF groups: the same attributes, dimensions, variables
 * (over dimensions of the same names and sizes, with the same values) and
 * subgroups, recursively. The first difference goes in why. */

static std::string att_value( const netCDF::NcAtt & a )
{
 if( a.getType() == netCDF::NcChar() ) {
  std::string s;
  a.getValues( s );
  return( s );
  }
 std::vector< double > v( a.getAttLength() );
 a.getValues( v.data() );
 std::ostringstream s;
 for( auto x : v )
  s << x << " ";
 return( s.str() );
 }

static bool same_group( const netCDF::NcGroup & a , const netCDF::NcGroup & b ,
			std::string & why )
{
 const std::string where = a.getName() + " vs " + b.getName() + ": ";

 auto aa = a.getAtts();
 auto ba = b.getAtts();
 if( aa.size() != ba.size() ) {
  why = where + "different number of attributes";
  return( false );
  }
 for( auto & [ n , at ] : aa ) {
  auto it = ba.find( n );
  if( ( it == ba.end() ) || ( att_value( at ) != att_value( it->second ) ) ) {
   why = where + "attribute " + n;
   return( false );
   }
  }

 auto ad = a.getDims();
 auto bd = b.getDims();
 if( ad.size() != bd.size() ) {
  why = where + "different number of dimensions";
  return( false );
  }
 for( auto & [ n , dim ] : ad ) {
  auto it = bd.find( n );
  if( ( it == bd.end() ) || ( dim.getSize() != it->second.getSize() ) ) {
   why = where + "dimension " + n;
   return( false );
   }
  }

 auto av = a.getVars();
 auto bv = b.getVars();
 if( av.size() != bv.size() ) {
  why = where + "different number of variables";
  return( false );
  }
 for( auto & [ n , var ] : av ) {
  auto it = bv.find( n );
  if( it == bv.end() ) {
   why = where + "variable " + n + " missing";
   return( false );
   }
  const auto & wr = it->second;
  auto vd = var.getDims();
  auto wd = wr.getDims();
  if( vd.size() != wd.size() ) {
   why = where + "variable " + n + " has different dimensions";
   return( false );
   }
  std::size_t sz = 1;
  for( std::size_t i = 0 ; i < vd.size() ; ++i ) {
   if( ( vd[ i ].getName() != wd[ i ].getName() ) ||
       ( vd[ i ].getSize() != wd[ i ].getSize() ) ) {
    why = where + "variable " + n + " has different dimensions";
    return( false );
    }
   sz *= vd[ i ].getSize();
   }
  if( sz == 0 )
   continue;
  if( var.getType() == netCDF::NcString() ) {
   if( wr.getType() != netCDF::NcString() ) {
    why = where + "variable " + n + " has different types";
    return( false );
    }
   std::vector< char * > x( sz ) , y( sz );
   var.getVar( x.data() );
   wr.getVar( y.data() );
   bool eq = true;
   for( std::size_t i = 0 ; i < sz ; ++i ) {
    eq &= ( std::string( x[ i ] ) == std::string( y[ i ] ) );
    free( x[ i ] );
    free( y[ i ] );
    }
   if( ! eq ) {
    why = where + "variable " + n + " has different values";
    return( false );
    }
   continue;
   }
  std::vector< double > x( sz ) , y( sz );
  var.getVar( x.data() );
  wr.getVar( y.data() );
  if( x != y ) {
   why = where + "variable " + n + " has different values";
   return( false );
   }
  }

 auto ag = a.getGroups();
 auto bg = b.getGroups();
 if( ag.size() != bg.size() ) {
  why = where + "different number of groups";
  return( false );
  }
 for( auto & [ n , grp ] : ag ) {
  auto it = bg.find( n );
  if( it == bg.end() ) {
   why = where + "group " + n + " missing";
   return( false );
   }
  if( ! same_group( grp , it->second , why ) )
   return( false );
  }

 return( true );
 }

/*--------------------------------------------------------------------------*/
/* The round trip of the Block described by the group g: it is read, written,
 * read again and written again, and the two groups written have to
 * coincide. Returns the Block read the second time, which the caller checks
 * against the data written in g and deletes. */

static Block * round_trip( const netCDF::NcGroup & g ,
			   const std::string & what )
{
 Block * b1 = nullptr;
 Block * b2 = nullptr;
 try {
  b1 = Block::new_Block( g );
  if( ! b1 ) {
   check( false , what + ": the Block is not read" );
   return( nullptr );
   }
  auto g1 = new_group( "RT" );
  b1->serialize( g1 );
  b2 = Block::new_Block( g1 );
  if( ! b2 ) {
   check( false , what + ": the Block written is not read back" );
   delete b1;
   return( nullptr );
   }
  auto g2 = new_group( "RT" );
  b2->serialize( g2 );
  std::string why;
  check( same_group( g1 , g2 , why ) , what + ": " + why );
  check( b2->classname() == b1->classname() ,
	 what + ": " + b2->classname() + " read back as " +
	 b1->classname() );
  }
 catch( std::exception & e ) {
  check( false , what + ": throws " + e.what() );
  delete b1;
  delete b2;
  return( nullptr );
  }
 delete b1;
 return( b2 );
 }

/*--------------------------------------------------------------------------*/
/// the spanning forest and the cycles of the CYCLE formulation
/** The DC graph has three parallel lines between nodes 0 and 1, one of them
 * reversed, a triangle 0, 1, 2, a second component 3, 4 with two opposite
 * lines, a node with no DC line and HVDC lines (zero susceptance) between
 * the components and towards a node reached only by them. The forest must
 * have one DC line per non-root node, oriented as its parent, of the three
 * parallel lines exactly one, and no HVDC line; every cycle must be a
 * circulation and every DC line out of the forest in some cycle. Then, on
 * two parallel DC lines 0 -> 1 and an HVDC line 2 -> 1 with efficiency 0.5,
 * a point that satisfies the nodal balance with the losses must satisfy the
 * constraints of the CYCLE formulation, which a flow of a tree line written
 * without the efficiency would violate, while one off the balance must
 * not. */

static void test_cycle_basis( void )
{
 auto g = new_group( "CY" , true );
 g.putAtt( "type" , "DCNetworkBlock" );
 auto N = g.addDim( "NumberNodes" , 7 );
 auto L = g.addDim( "NumberLines" , 9 );
 //                  0  1  2  3  4  5  6  7  8
 put_int( g , "StartLine" , L , { 0 , 0 , 1 , 1 , 2 , 3 , 4 , 2 , 5 } );
 put_int( g , "EndLine" ,   L , { 1 , 1 , 0 , 2 , 0 , 4 , 3 , 3 , 6 } );
 put( g , "MinPowerFlow" , L , std::vector< double >( 9 , -100 ) );
 put( g , "MaxPowerFlow" , L , std::vector< double >( 9 , 100 ) );
 put( g , "LineSusceptance" , L , { 1 , 2 , 3 , 1 , 1 , 1 , 1 , 0 , 0 } );
 put( g , "ActiveDemand" , N , std::vector< double >( 7 , 0 ) );

 std::unique_ptr< Block > b;
 try {
  b.reset( Block::new_Block( g ) );
  }
 catch( std::exception & e ) {
  check( false , std::string( "cycle basis: reading throws " ) + e.what() );
  return;
  }
 auto nb = dynamic_cast< DCNetworkBlock * >( b.get() );
 auto nd = nb ? dynamic_cast< DCNetworkBlock::DCNetworkData * >(
                                         nb->get_NetworkData() ) : nullptr;
 if( ! nd ) {
  check( false , "cycle basis: the DCNetworkBlock is not read" );
  return;
  }

 const auto & st = nd->get_start_line();
 const auto & en = nd->get_end_line();
 const auto tree = nd->get_lines_in_spanning_tree();
 const auto cycles = nd->get_lines_in_cycles();
 const auto & parent = nd->get_spanning_parent();

 check( tree.size() == 3 , "cycle basis: 3 DC lines in the forest, " +
	std::to_string( tree.size() ) );
 check( cycles.size() == 4 , "cycle basis: 4 cycles, " +
	std::to_string( cycles.size() ) );

 int parallel = 0;
 for( const auto & [ l , sign ] : tree ) {
  const Index child = ( sign > 0 ) ? en[ l ] : st[ l ];
  const Index father = ( sign > 0 ) ? st[ l ] : en[ l ];
  check( parent[ child ] == int( father ) ,
	 "cycle basis: tree line " + std::to_string( l ) + " oriented as "
	 "the parent of its child" );
  check( l < 7 , "cycle basis: HVDC line " + std::to_string( l ) +
	 " in the forest" );
  if( l < 3 )
   ++parallel;
  }
 check( parallel == 1 , "cycle basis: " + std::to_string( parallel ) +
	" of the 3 parallel lines in the forest" );

 std::vector< bool > covered( 9 , false );
 for( const auto & [ l , sign ] : tree )
  covered[ l ] = true;
 for( std::size_t c = 0 ; c < cycles.size() ; ++c ) {
  std::vector< int > balance( 7 , 0 );
  for( const auto & [ l , sign ] : cycles[ c ] ) {
   check( l < 7 , "cycle basis: HVDC line in a cycle" );
   check( ( sign == 1 ) || ( sign == -1 ) ,
	  "cycle basis: a line twice in a cycle" );
   covered[ l ] = true;
   if( st[ l ] != en[ l ] ) {
    balance[ st[ l ] ] -= sign;
    balance[ en[ l ] ] += sign;
    }
   }
  check( std::all_of( balance.begin() , balance.end() ,
		      []( int x ) { return( x == 0 ); } ) ,
	 "cycle basis: cycle " + std::to_string( c ) +
	 " is not a circulation" );
  }
 for( Index l = 0 ; l < 7 ; ++l )
  check( covered[ l ] , "cycle basis: DC line " + std::to_string( l ) +
	 " neither in the forest nor in a cycle" );
 check( parent[ 5 ] == -1 && parent[ 6 ] == -1 ,
	"cycle basis: nodes with no DC line have no parent" );

 // the HVDC line with losses in the flows of the CYCLE formulation
 g = new_group( "CY" , true );
 g.putAtt( "type" , "DCNetworkBlock" );
 N = g.addDim( "NumberNodes" , 3 );
 L = g.addDim( "NumberLines" , 3 );
 put_int( g , "StartLine" , L , { 0 , 0 , 2 } );
 put_int( g , "EndLine" ,   L , { 1 , 1 , 1 } );
 put( g , "MinPowerFlow" , L , { -100 , -100 , 0 } );
 put( g , "MaxPowerFlow" , L , { 100 , 100 , 100 } );
 put( g , "LineSusceptance" , L , { 1 , 1 , 0 } );
 put( g , "Efficiency" , L , { 1 , 1 , 0.5 } );
 put( g , "ActiveDemand" , N , { 4 , 1 , 0 } );

 b.reset( Block::new_Block( g ) );
 nb = dynamic_cast< DCNetworkBlock * >( b.get() );
 if( ! nb ) {
  check( false , "cycle flows: the DCNetworkBlock is not read" );
  return;
  }
 SimpleConfiguration< int > cycle( 1 );
 nb->generate_abstract_variables( & cycle );
 nb->generate_abstract_constraints();

 // node 2 injects 10 on the HVDC line, which brings 5 to node 1, whose
 // demand is 1: the 4 left go to node 0 over the two parallel lines, 2 each
 // as their susceptances are the same, with the cycle flow that makes it so
 auto set_point = [ & ]( double hvdc ) {
  auto inj = nb->get_node_injection( 0 );
  inj[ 0 ].set_value( 0 );
  inj[ 1 ].set_value( 0 );
  inj[ 2 ].set_value( hvdc );
  auto & f = nb->get_power_flow();
  const_cast< ColVariable & >( f[ 0 ] ).set_value( -2 );
  const_cast< ColVariable & >( f[ 1 ] ).set_value( -2 );
  const_cast< ColVariable & >( f[ 2 ] ).set_value( hvdc );
  auto & h = nb->get_cycle_flow();
  auto nd2 = static_cast< DCNetworkBlock::DCNetworkData * >(
                                                  nb->get_NetworkData() );
  const auto tree2 = nd2->get_lines_in_spanning_tree();
  const auto cycles2 = nd2->get_lines_in_cycles();
  for( std::size_t c = 0 ; c < h.size() ; ++c )
   for( const auto & [ l , sign ] : cycles2[ c ] )
    if( ! tree2.contains( l ) )
     const_cast< ColVariable & >( h[ c ] ).set_value( -2.0 * sign );
  };

 set_point( 10 );
 check( nb->is_feasible() ,
	"cycle flows: the balance with the losses of the HVDC line" );
 set_point( 5 );  // node 1 gets 2.5 from the HVDC line instead of 5
 check( ! nb->is_feasible() , "cycle flows: a point off the balance" );
 }

/*--------------------------------------------------------------------------*/
/// ThermalUnitBlock and NuclearUnitBlock, with all the data and with the
/// fewest

static void test_RT_thermal( void )
{
 // all the data the test writes, over 3 instants
 TUData d;
 d.T = 3;
 d.minP = 2;
 d.maxP = 8;
 d.ru = 3;
 d.rd = 4;
 d.su = 5;
 d.sd = 6;
 d.initUD = 2;
 d.initP = 4;
 d.minUp = 3;
 d.minDown = 2;
 d.lin = { 1 , 2 , 3 };
 d.quad = { 0.1 , 0.2 , 0.3 };
 d.cnst = { 4 , 5 , 6 };
 d.suc = { 7 , 8 , 9 };
 d.sdc = { 1 , 0 , 1 };
 d.scale = 2;
 auto g = new_group( "TU" , true );
 write_TU( g , d );
 {
  auto nd = g.getDim( "NumberIntervals" );
  put( g , "PrimaryRho" , nd , { 0.1 , 0.1 , 0.2 } );
  put( g , "SecondaryRho" , nd , { 0.2 , 0.1 , 0.1 } );
  put( g , "Availability" , nd , { 1 , 1 , 0.75 } );
  put( g , "FixedConsumption" , nd , { 0.5 , 0.5 , 0.5 } );
  put( g , "InertiaCommitment" , nd , { 1 , 2 , 3 } );
  }
 auto b = dynamic_cast< ThermalUnitBlock * >( round_trip( g , "thermal" ) );
 if( b ) {
  check( b->get_time_horizon() == 3 , "thermal: time horizon" );
  check( b->get_min_up_time() == 3 , "thermal: MinUpTime" );
  check( b->get_min_down_time() == 2 , "thermal: MinDownTime" );
  check( b->get_init_up_down_time() == 2 , "thermal: InitUpDownTime" );
  check( b->get_initial_power() == 4 , "thermal: InitialPower" );
  check( b->get_scale() == 2 , "thermal: Scale" );
  check( b->get_delta_ramp_up( 1 ) == 3 , "thermal: DeltaRampUp" );
  check( b->get_delta_ramp_down( 2 ) == 4 , "thermal: DeltaRampDown" );
  check( b->get_linear_term( 2 ) == 3 , "thermal: LinearTerm" );
  check( b->get_quad_term( 1 ) == 0.2 , "thermal: QuadTerm" );
  check( b->get_const_term( 0 ) == 4 , "thermal: ConstTerm" );
  check( b->get_start_up_cost() == std::vector< double >( { 7 , 8 , 9 } ) ,
	 "thermal: StartUpCost" );
  check( b->get_shut_down_cost() == std::vector< double >( { 1 , 0 , 1 } ) ,
	 "thermal: ShutDownCost" );
  check( b->get_availability( 2 ) == 0.75 , "thermal: Availability" );
  check( b->get_max_power( 0 ) == 8 , "thermal: MaxPower" );
  check( b->get_min_power( 2 ) == 2 , "thermal: MinPower" );
  check( b->get_start_up_limit()[ 0 ] == 5 , "thermal: StartUpLimit" );
  check( b->get_shut_down_limit()[ 0 ] == 6 , "thermal: ShutDownLimit" );
  check( b->get_primary_rho()[ 2 ] == 0.2 , "thermal: PrimaryRho" );
  delete b;
  }

 // the fewest data over a single instant: the optional costs stay empty
 g = new_group( "TU" , true );
 g.putAtt( "type" , "ThermalUnitBlock" );
 g.addDim( "TimeHorizon" , 1 );
 put( g , "MinPower" , 0 );
 put( g , "MaxPower" , 5 );
 b = dynamic_cast< ThermalUnitBlock * >( round_trip( g , "thermal, fewest" ) );
 if( b ) {
  check( b->get_time_horizon() == 1 , "thermal, fewest: time horizon" );
  // the costs not given are 0, and a free shut-down keeps its vector empty
  check( ( b->get_linear_term( 0 ) == 0 ) && ( b->get_quad_term( 0 ) == 0 ) &&
	 ( b->get_const_term( 0 ) == 0 ) ,
	 "thermal, fewest: the costs not given are not 0" );
  check( std::all_of( b->get_start_up_cost().begin() ,
		      b->get_start_up_cost().end() ,
		      []( double c ) { return( c == 0 ); } ) ,
	 "thermal, fewest: StartUpCost" );
  check( b->get_shut_down_cost().empty() , "thermal, fewest: ShutDownCost" );
  check( b->get_max_power( 0 ) == 5 , "thermal, fewest: MaxPower" );
  check( b->get_scale() == 1 , "thermal, fewest: Scale" );
  delete b;
  }

 // a nuclear unit
 d.nuclear = true;
 d.scale = 1;
 d.modT = 3;
 d.initMod = 2;
 d.mru = 1;
 d.mrd = 2;
 g = new_group( "NU" , true );
 write_TU( g , d );
 auto nb = dynamic_cast< NuclearUnitBlock * >( round_trip( g , "nuclear" ) );
 if( nb ) {
  check( nb->get_min_up_time() == 3 , "nuclear: MinUpTime" );
  check( nb->get_linear_term( 1 ) == 2 , "nuclear: LinearTerm" );
  delete nb;
  }
 }

/*--------------------------------------------------------------------------*/
/// HydroUnitBlock, one reservoir and one arc, and with the fewest data

static void test_RT_hydro( void )
{
 auto g = new_group( "HU" , true );
 g.putAtt( "type" , "HydroUnitBlock" );
 auto T = g.addDim( "TimeHorizon" , 3 );
 auto R = g.addDim( "NumberReservoirs" , 1 );
 auto A = g.addDim( "NumberArcs" , 1 );
 put_int( g , "StartArc" , A , { 0 } );
 put_int( g , "EndArc" , A , { 1 } );
 put( g , "Inflows" , { R , T } , { 1 , 2 , 3 } );
 put( g , "MinFlow" , { T , A } , { 0 , 0 , 0 } );
 put( g , "MaxFlow" , { T , A } , { 10 , 10 , 10 } );
 put( g , "MinPower" , { T , A } , { 0 , 0 , 0 } );
 put( g , "MaxPower" , { T , A } , { 20 , 20 , 20 } );
 put( g , "LinearTerm" , A , { 2 } );
 put( g , "ConstantTerm" , A , { 0 } );
 put( g , "InitialVolumetric" , R , { 50 } );
 put( g , "MinVolumetric" , { R , T } , { 0 , 0 , 0 } );
 put( g , "MaxVolumetric" , { R , T } , { 100 , 100 , 100 } );
 auto b = dynamic_cast< HydroUnitBlock * >( round_trip( g , "hydro" ) );
 if( b ) {
  check( b->get_time_horizon() == 3 , "hydro: time horizon" );
  check( b->get_number_reservoirs() == 1 , "hydro: NumberReservoirs" );
  check( b->get_number_generators() == 1 , "hydro: NumberArcs" );
  check( b->get_initial_volumetric( 0 ) == 50 , "hydro: InitialVolumetric" );
  delete b;
  }
 }

/*--------------------------------------------------------------------------*/
/// BatteryUnitBlock, IntermittentUnitBlock and SlackUnitBlock

static void test_RT_other_units( void )
{
 auto g = new_group( "BU" , true );
 g.putAtt( "type" , "BatteryUnitBlock" );
 auto T = g.addDim( "TimeHorizon" , 2 );
 put( g , "MinStorage" , T , { 0 , 0 } );
 put( g , "MaxStorage" , T , { 10 , 10 } );
 put( g , "MaxPower" , T , { 5 , 5 } );
 put( g , "MinPower" , T , { -4 , -4 } );
 put( g , "InitialStorage" , 3 );
 put( g , "StoringBatteryRho" , T , { 0.8 , 0.8 } );
 put( g , "ExtractingBatteryRho" , T , { 0.9 , 0.9 } );
 put( g , "Cost" , T , { 1 , 2 } );
 auto b = round_trip( g , "battery" );
 if( auto bb = dynamic_cast< BatteryUnitBlock * >( b ) ) {
  check( bb->get_time_horizon() == 2 , "battery: time horizon" );
  check( bb->get_initial_storage() == 3 , "battery: InitialStorage" );
  check( bb->get_max_power( 1 ) == 5 , "battery: MaxPower" );
  check( bb->get_min_power( 0 ) == -4 , "battery: MinPower" );
  }
 delete b;

 g = new_group( "IU" , true );
 g.putAtt( "type" , "IntermittentUnitBlock" );
 T = g.addDim( "TimeHorizon" , 2 );
 put( g , "MaxPower" , T , { 7 , 3 } );
 b = round_trip( g , "intermittent" );
 if( auto ib = dynamic_cast< IntermittentUnitBlock * >( b ) ) {
  check( ib->get_max_power( 0 ) == 7 , "intermittent: MaxPower" );
  check( ib->get_max_power( 1 ) == 3 , "intermittent: MaxPower" );
  check( ib->get_min_power( 1 ) == 0 , "intermittent: MinPower" );
  }
 delete b;

 g = new_group( "IU" , true );
 g.putAtt( "type" , "IntermittentUnitBlock" );
 T = g.addDim( "TimeHorizon" , 1 );
 put( g , "MaxPower" , T , { 2 } );
 put( g , "InvestmentCost" , 5 );
 put( g , "MaxCapacityDesign" , 1 );
 put( g , "Scale" , 2 );
 b = round_trip( g , "intermittent, design" );
 if( auto ib = dynamic_cast< IntermittentUnitBlock * >( b ) )
  check( ib->get_scale() == 2 , "intermittent, design: Scale" );
 delete b;

 g = new_group( "SU" , true );
 g.putAtt( "type" , "SlackUnitBlock" );
 T = g.addDim( "TimeHorizon" , 2 );
 put( g , "MaxPower" , T , { 100 , 100 } );
 put( g , "ActivePowerCost" , T , { 1000 , 1000 } );
 b = round_trip( g , "slack" );
 if( auto sb = dynamic_cast< SlackUnitBlock * >( b ) )
  check( sb->get_max_power( 1 ) == 100 , "slack: MaxPower" );
 delete b;

 g = new_group( "SU" , true );  // with no data at all but the horizon
 g.putAtt( "type" , "SlackUnitBlock" );
 g.addDim( "TimeHorizon" , 1 );
 delete round_trip( g , "slack, fewest" );
 }

/*--------------------------------------------------------------------------*/
/// DCNetworkBlock and ACNetworkBlock with their own NetworkData, and
/// ECNetworkBlock

static void test_RT_networks( void )
{
 auto write_DC = [ & ]( netCDF::NcGroup g , const std::string & type ) {
  g.putAtt( "type" , type );
  auto N = g.addDim( "NumberNodes" , 3 );
  auto L = g.addDim( "NumberLines" , 2 );
  put_int( g , "StartLine" , L , { 0 , 1 } );
  put_int( g , "EndLine" , L , { 1 , 2 } );
  put( g , "MinPowerFlow" , L , { -5 , -6 } );
  put( g , "MaxPowerFlow" , L , { 5 , 6 } );
  put( g , "LineSusceptance" , L , { 0 , 2 } );
  put( g , "ActiveDemand" , N , { 1 , 2 , 3 } );
  return( N );
  };

 auto g = new_group( "DC" , true );
 write_DC( g , "DCNetworkBlock" );
 auto b = round_trip( g , "DC network" );
 if( auto nb = dynamic_cast< DCNetworkBlock * >( b ) ) {
  check( nb->get_NetworkData() &&
	 ( nb->get_NetworkData()->get_number_nodes() == 3 ) ,
	 "DC network: NumberNodes" );
  auto nd = dynamic_cast< DCNetworkBlock::DCNetworkData * >(
					    nb->get_NetworkData() );
  check( nd && ( nd->get_number_lines() == 2 ) ,
	 "DC network: NumberLines" );
  auto ad = nb->get_active_demand( 0 );
  check( ad && ( ad[ 2 ] == 3 ) , "DC network: ActiveDemand" );
  }
 delete b;

 g = new_group( "DC" , true );  // a single node: no line at all
 g.putAtt( "type" , "DCNetworkBlock" );
 auto N = g.addDim( "NumberNodes" , 1 );
 put( g , "ActiveDemand" , N , { 4 } );
 delete round_trip( g , "DC network, one node" );

 g = new_group( "AC" , true );
 N = write_DC( g , "ACNetworkBlock" );
 {
  auto L = g.getDim( "NumberLines" );
  put( g , "LineReactance" , L , { 0.1 , 0.2 } );
  put( g , "LineResistance" , L , { 0.01 , 0.02 } );
  put( g , "NodeMaxVoltage" , N , { 1.1 , 1.1 , 1.1 } );
  put( g , "NodeMinVoltage" , N , { 0.9 , 0.9 , 0.9 } );
  put( g , "ReactiveDemand" , N , { 0.5 , 0.5 , 0.5 } );
  }
 delete round_trip( g , "AC network" );

 g = new_group( "EC" , true );
 g.putAtt( "type" , "ECNetworkBlock" );
 auto I = g.addDim( "NumberIntervals" , 2 );
 N = g.addDim( "NumberNodes" , 2 );  // an energy community of two users
 put( g , "BuyPrice" , I , { 3 , 4 } );
 put( g , "SellPrice" , I , { 1 , 1 } );
 put( g , "PeakTariff" , 2 );
 put( g , "ActiveDemand" , { I , N } , { 5 , 6 , 7 , 8 } );
 delete round_trip( g , "EC network" );
 }

/*--------------------------------------------------------------------------*/
/// UCBlock with no unit, and with one unit of each kind on two nodes

static void test_RT_UCBlock( void )
{
 auto g = new_group( "UC" , true );
 g.putAtt( "type" , "UCBlock" );
 auto T = g.addDim( "TimeHorizon" , 2 );
 g.addDim( "NumberUnits" , 0 );
 auto one = g.addDim( "NumberNodes" , 1 );
 put( g , "ActivePowerDemand" , { one , T } , { 0 , 0 } );
 auto b = dynamic_cast< UCBlock * >( round_trip( g , "UCBlock, no unit" ) );
 if( b ) {
  check( b->get_number_units() == 0 , "UCBlock, no unit: NumberUnits" );
  check( b->get_time_horizon() == 2 , "UCBlock, no unit: TimeHorizon" );
  try {
   generate_all( b );
   }
  catch( std::exception & e ) {
   check( false , std::string( "UCBlock, no unit: the abstract "
			       "representation throws " ) + e.what() );
   }
  delete b;
  }

 g = new_group( "UC" , true );
 g.putAtt( "type" , "UCBlock" );
 T = g.addDim( "TimeHorizon" , 2 );
 g.addDim( "NumberUnits" , 5 );
 auto N = g.addDim( "NumberNodes" , 2 );
 auto L = g.addDim( "NumberLines" , 1 );
 auto G = g.addDim( "NumberElectricalGenerators" , 5 );
 put_int( g , "StartLine" , L , { 0 } );
 put_int( g , "EndLine" , L , { 1 } );
 put( g , "MinPowerFlow" , L , { -10 } );
 put( g , "MaxPowerFlow" , L , { 10 } );
 put( g , "ActivePowerDemand" , { N , T } , { 3 , 4 , 5 , 6 } );
 put_int( g , "GeneratorNode" , G , { 0 , 1 , 0 , 1 , 0 } );

 TUData d;
 d.T = 2;
 d.lin = { 1 , 2 };
 d.initUD = -1;
 auto u = g.addGroup( "UnitBlock_0" );
 write_TU( u , d );

 u = g.addGroup( "UnitBlock_1" );
 u.putAtt( "type" , "HydroUnitBlock" );
 auto R = u.addDim( "NumberReservoirs" , 1 );
 auto A = u.addDim( "NumberArcs" , 1 );
 put_int( u , "StartArc" , A , { 0 } );
 put_int( u , "EndArc" , A , { 1 } );
 put( u , "Inflows" , { R , T } , { 1 , 2 } );
 put( u , "MaxFlow" , { T , A } , { 10 , 10 } );
 put( u , "MaxPower" , { T , A } , { 20 , 20 } );
 put( u , "LinearTerm" , A , { 2 } );
 put( u , "InitialVolumetric" , R , { 50 } );
 put( u , "MaxVolumetric" , { R , T } , { 100 , 100 } );

 u = g.addGroup( "UnitBlock_2" );
 u.putAtt( "type" , "BatteryUnitBlock" );
 put( u , "MinStorage" , T , { 0 , 0 } );
 put( u , "MaxStorage" , T , { 10 , 10 } );
 put( u , "MaxPower" , T , { 5 , 5 } );

 u = g.addGroup( "UnitBlock_3" );
 u.putAtt( "type" , "IntermittentUnitBlock" );
 put( u , "MaxPower" , T , { 7 , 3 } );

 u = g.addGroup( "UnitBlock_4" );
 u.putAtt( "type" , "SlackUnitBlock" );
 put( u , "MaxPower" , T , { 100 , 100 } );
 put( u , "ActivePowerCost" , T , { 1000 , 1000 } );

 b = dynamic_cast< UCBlock * >( round_trip( g , "UCBlock, two nodes" ) );
 if( b ) {
  check( b->get_number_units() == 5 , "UCBlock: NumberUnits" );
  check( b->get_number_nodes() == 2 , "UCBlock: NumberNodes" );
  check( b->get_generator_node() ==
	 std::vector< Index >( { 0 , 1 , 0 , 1 , 0 } ) ,
	 "UCBlock: GeneratorNode" );
  check( dynamic_cast< ThermalUnitBlock * >( b->get_unit_block( 0 ) ) &&
	 dynamic_cast< HydroUnitBlock * >( b->get_unit_block( 1 ) ) &&
	 dynamic_cast< BatteryUnitBlock * >( b->get_unit_block( 2 ) ) &&
	 dynamic_cast< IntermittentUnitBlock * >( b->get_unit_block( 3 ) ) &&
	 dynamic_cast< SlackUnitBlock * >( b->get_unit_block( 4 ) ) ,
	 "UCBlock: the kinds of the units" );
  check( b->get_number_networks() == 2 , "UCBlock: NumberNetworks" );
  delete b;
  }
 }

/*--------------------------------------------------------------------------*/
/*-------------------------------- SETTERS ---------------------------------*/
/*--------------------------------------------------------------------------*/
/// the Modification the FakeSolver received, GroupModification unpacked

static void flatten( const sp_Mod & m ,
		     std::vector< sp_Mod > & out )
{
 if( auto gm = std::dynamic_pointer_cast< const GroupModification >( m ) ) {
  for( const auto & s : gm->sub_Modifications() )
   flatten( s , out );
  }
 else
  out.push_back( m );
 }

static std::vector< sp_Mod > received( FakeSolver * fs )
{
 std::vector< sp_Mod > out;
 for( const auto & m : fs->get_Modification_list() )
  flatten( m , out );
 return( out );
 }

/*--------------------------------------------------------------------------*/
/// the physical Modification of the given type among those received

static int count_TUB_mods( FakeSolver * fs , int type )
{
 int n = 0;
 for( const auto & m : received( fs ) )
  if( auto um = std::dynamic_pointer_cast< const UnitBlockMod >( m ) )
   if( um->type() == type )
    ++n;
 return( n );
 }

/*--------------------------------------------------------------------------*/
/// the coefficients of the Objective of the unit on the Variable v

static double lin_coef( ThermalUnitBlock * tub , const ColVariable * v )
{
 auto qf = static_cast< DQuadFunction * >(
	   static_cast< FRealObjective * >( tub->get_objective() )
                                                         ->get_function() );
 return( qf->get_linear_coefficient( qf->is_active( v ) ) );
 }

static double quad_coef( ThermalUnitBlock * tub , const ColVariable * v )
{
 auto qf = static_cast< DQuadFunction * >(
	   static_cast< FRealObjective * >( tub->get_objective() )
                                                         ->get_function() );
 return( qf->get_quadratic_coefficient( qf->is_active( v ) ) );
 }

/*--------------------------------------------------------------------------*/
/* The cost setters of a ThermalUnitBlock, both in the Range and in the
 * Subset form: the FakeSolver receives the physical Modification, the data
 * change, and so does the Objective, which carries the scale factor times
 * the cost. An unchanged value issues nothing, eNoMod changes without
 * issuing and eDryRun does not change. The DP Solver attached to the unit
 * follows every change. */

static void test_setters_thermal( double scale )
{
 const auto what = "setters, scale " + str( scale );
 TUData d;
 d.T = 3;
 d.maxP = d.su = d.sd = d.ru = d.rd = 5;
 d.initUD = -10;
 d.lin = { -10 , -10 , -10 };
 d.quad = { 0 , 0 , 0 };
 d.cnst = { 1 , 1 , 1 };
 d.suc = { 2 , 2 , 2 };
 d.sdc = { 1 , 1 , 1 };
 d.scale = scale;
 auto tub = new_TU( d );
 generate_all( tub );

 auto fs = new FakeSolver();
 tub->register_Solver( fs );
 auto dp = Solver::new_Solver( "ThermalUnitExtDPSolver" );
 tub->register_Solver( dp );
 fs->get_Modification_list().clear();

 auto p = tub->get_active_power( 0 );
 auto u = tub->get_commitment( 0 );
 auto v = tub->get_start_up();
 auto w = tub->get_shut_down();

 auto resolve = [ & ]( const std::string & when ) {
  const int status = dp->compute();
  const auto bf = brute_force( d , std::vector< int >( d.T , -1 ) );
  check( ( status == Solver::kOK ) && close( dp->get_ub() , bf ) ,
	 what + ", " + when + ": the DP gives " + str( dp->get_ub() ) +
	 " instead of " + str( bf ) );
  };
 resolve( "at the start" );

 // linear term, Range
 std::vector< double > val = { 3 , 4 };
 tub->set_linear_term( val.begin() , Range( 1 , 3 ) );
 d.lin = { -10 , 3 , 4 };
 check( count_TUB_mods( fs , ThermalUnitBlockMod::eSetLinT ) == 1 ,
	what + ": set_linear_term( Range ) issues no eSetLinT" );
 check( tub->get_linear_term( 2 ) == 4 ,
	what + ": set_linear_term( Range ) does not change the data" );
 check( close( lin_coef( tub , & p[ 1 ] ) , scale * 3 ) &&
	close( lin_coef( tub , & p[ 2 ] ) , scale * 4 ) ,
	what + ": set_linear_term( Range ) puts " +
	str( lin_coef( tub , & p[ 2 ] ) ) + " in the Objective, not " +
	str( scale * 4 ) );
 resolve( "linear term changed" );
 fs->get_Modification_list().clear();

 // the same values again: nothing
 tub->set_linear_term( val.begin() , Range( 1 , 3 ) );
 check( fs->get_Modification_list().empty() ,
	what + ": an unchanged linear term issues a Modification" );

 // linear term, Subset
 val = { -7 };
 tub->set_linear_term( val.begin() , Subset( { 0 } ) );
 d.lin[ 0 ] = -7;
 check( count_TUB_mods( fs , ThermalUnitBlockMod::eSetLinT ) == 1 ,
	what + ": set_linear_term( Subset ) issues no eSetLinT" );
 check( close( lin_coef( tub , & p[ 0 ] ) , scale * -7 ) ,
	what + ": set_linear_term( Subset ) puts " +
	str( lin_coef( tub , & p[ 0 ] ) ) + " in the Objective" );
 resolve( "linear term changed by Subset" );
 fs->get_Modification_list().clear();

 val = { -8 };
 tub->set_linear_term( val.begin() , Range( 0 , 1 ) , eNoMod , eNoMod );
 // eNoMod: the data and the Objective change, and the Objective issues no
 // Modification
 {
  int n = 0;
  for( const auto & m : received( fs ) )
   if( ! std::dynamic_pointer_cast< const UnitBlockMod >( m ) )
    ++n;
  check( n == 0 , what + ": eNoMod issues " + std::to_string( n ) +
	 " abstract Modification" );
  }
 check( tub->get_linear_term( 0 ) == -8 ,
	what + ": eNoMod does not change the data" );
 check( close( lin_coef( tub , & p[ 0 ] ) , scale * -8 ) ,
	what + ": eNoMod does not change the Objective" );
 d.lin[ 0 ] = -8;

 // eDryRun: nothing changes
 val = { -9 };
 tub->set_linear_term( val.begin() , Range( 0 , 1 ) , eDryRun , eDryRun );
 check( tub->get_linear_term( 0 ) == -8 ,
	what + ": eDryRun changes the data" );
 check( close( lin_coef( tub , & p[ 0 ] ) , scale * -8 ) ,
	what + ": eDryRun changes the Objective" );
 // the DP was not told of the eNoMod change: it is reloaded
 tub->unregister_Solver( dp , true );
 dp = Solver::new_Solver( "ThermalUnitExtDPSolver" );
 tub->register_Solver( dp );
 resolve( "reloaded" );
 fs->get_Modification_list().clear();

 // constant term
 val = { 5 , 6 };
 tub->set_const_term( val.begin() , Range( 0 , 2 ) );
 d.cnst = { 5 , 6 , 1 };
 check( count_TUB_mods( fs , ThermalUnitBlockMod::eSetConstT ) == 1 ,
	what + ": set_const_term issues no eSetConstT" );
 check( close( lin_coef( tub , & u[ 1 ] ) , scale * 6 ) ,
	what + ": set_const_term puts " + str( lin_coef( tub , & u[ 1 ] ) ) +
	" in the Objective" );
 resolve( "constant term changed" );
 fs->get_Modification_list().clear();

 // start-up cost; the unit is off since long, so every v_t is there
 val = { 9 };
 tub->set_startup_costs( val.begin() , Subset( { 2 } ) );
 d.suc = { 2 , 2 , 9 };
 check( count_TUB_mods( fs , ThermalUnitBlockMod::eSetSUC ) == 1 ,
	what + ": set_startup_costs issues no eSetSUC" );
 check( v && close( lin_coef( tub , & v[ 2 ] ) , scale * 9 ) ,
	what + ": set_startup_costs puts " +
	( v ? str( lin_coef( tub , & v[ 2 ] ) ) : "nothing" ) +
	" in the Objective" );
 resolve( "start-up cost changed" );
 fs->get_Modification_list().clear();

 // shut-down cost
 val = { 3 , 3 , 3 };
 tub->set_shutdown_costs( val.begin() , Range( 0 , 3 ) );
 d.sdc = { 3 , 3 , 3 };
 check( count_TUB_mods( fs , ThermalUnitBlockMod::eSetSDC ) == 1 ,
	what + ": set_shutdown_costs issues no eSetSDC" );
 check( w && close( lin_coef( tub , & w[ 1 ] ) , scale * 3 ) ,
	what + ": set_shutdown_costs puts " +
	( w ? str( lin_coef( tub , & w[ 1 ] ) ) : "nothing" ) +
	" in the Objective" );
 resolve( "shut-down cost changed" );
 fs->get_Modification_list().clear();

 // quadratic term: the linear coefficient of p_t stays
 val = { 0.5 };
 tub->set_quad_term( val.begin() , Range( 1 , 2 ) );
 check( count_TUB_mods( fs , ThermalUnitBlockMod::eSetQuadT ) == 1 ,
	what + ": set_quad_term issues no eSetQuadT" );
 check( close( quad_coef( tub , & p[ 1 ] ) , scale * 0.5 ) ,
	what + ": set_quad_term puts " + str( quad_coef( tub , & p[ 1 ] ) ) +
	" in the Objective" );
 check( close( lin_coef( tub , & p[ 1 ] ) , scale * 3 ) ,
	what + ": set_quad_term changes the linear coefficient to " +
	str( lin_coef( tub , & p[ 1 ] ) ) );
 val = { 0 };  // back, so that the brute force stays exact
 tub->set_quad_term( val.begin() , Range( 1 , 2 ) );
 resolve( "quadratic term changed back" );
 fs->get_Modification_list().clear();

 // the scale factor itself
 val = { 2 * scale + 1 };
 tub->scale( val.begin() , Subset( { 0 } ) );
 check( count_TUB_mods( fs , UnitBlockMod::eScale ) == 1 ,
	what + ": scale() issues no eScale" );
 check( tub->get_scale() == 2 * scale + 1 ,
	what + ": scale() does not change the scale factor" );
 check( close( lin_coef( tub , & p[ 2 ] ) , ( 2 * scale + 1 ) * 4 ) ,
	what + ": scale() puts " + str( lin_coef( tub , & p[ 2 ] ) ) +
	" in the Objective" );
 d.scale = 2 * scale + 1;
 resolve( "scale changed" );

 tub->unregister_Solver( dp , true );
 tub->unregister_Solver( fs , true );
 delete tub;
 }

/*--------------------------------------------------------------------------*/
/// the demand of a UCBlock, which the node balance has to follow

static void test_setters_UCBlock( void )
{
 auto g = new_group( "UC" , true );
 g.putAtt( "type" , "UCBlock" );
 auto T = g.addDim( "TimeHorizon" , 2 );
 g.addDim( "NumberUnits" , 1 );
 auto one = g.addDim( "NumberNodes" , 1 );
 put( g , "ActivePowerDemand" , { one , T } , { 3 , 4 } );
 TUData d;
 d.T = 2;
 d.lin = { 1 , 2 };
 auto u = g.addGroup( "UnitBlock_0" );
 write_TU( u , d );

 auto b = dynamic_cast< UCBlock * >( Block::new_Block( g ) );
 if( ! b ) {
  check( false , "UCBlock setters: the UCBlock is not read" );
  return;
  }
 generate_all( b );
 auto fs = new FakeSolver();
 b->register_Solver( fs );
 fs->get_Modification_list().clear();

 std::vector< double > val = { 7 };
 b->set_active_power_demand( val.begin() , Range( 1 , 2 ) );
 int n = 0;
 for( const auto & m : received( fs ) )
  if( auto um = std::dynamic_pointer_cast< const UCBlockMod >( m ) )
   if( um->type() == UCBlockMod::eSetActD )
    ++n;
 check( n == 1 , "UCBlock: set_active_power_demand issues no eSetActD" );
 check( b->get_active_power_demand()[ 0 ][ 1 ] == 7 ,
	"UCBlock: set_active_power_demand does not change the demand" );
 const auto & nic = b->get_const_node_injection_constraints();
 check( ( nic[ 1 ][ 0 ].get_lhs() == 7 ) && ( nic[ 1 ][ 0 ].get_rhs() == 7 ) ,
	"UCBlock: the node balance has [ " + str( nic[ 1 ][ 0 ].get_lhs() ) +
	" , " + str( nic[ 1 ][ 0 ].get_rhs() ) + " ], not the demand 7" );
 check( ( nic[ 0 ][ 0 ].get_lhs() == 3 ) && ( nic[ 0 ][ 0 ].get_rhs() == 3 ) ,
	"UCBlock: the node balance of the instant not changed moved" );

 b->unregister_Solver( fs , true );
 delete b;
 }

/*--------------------------------------------------------------------------*/
/*---------------------------------- MAIN ----------------------------------*/
/*--------------------------------------------------------------------------*/

int main( void )
{
 int ret = 0;
 try {
  test_DP_one_instant();
  test_DP_ramps();
  test_DP_boundary();
  test_DP_random_linear();
  test_DP_random_quadratic();
  test_DP_fixed();
  test_DP_set_min_up_down();
  test_DP_clamped_min_up_down();
  test_DP_nuclear();
  test_sol_feasible();

  test_RT_thermal();
  test_RT_hydro();
  test_RT_other_units();
  test_RT_networks();
  test_cycle_basis();
  test_RT_UCBlock();

  test_setters_thermal( 1 );
  test_setters_thermal( 3 );
  test_setters_UCBlock();
  }
 catch( std::exception & e ) {
  std::cout << "uncaught exception: " << e.what() << std::endl;
  ++failures;
  }

 if( NCID >= 0 )
  nc_close( NCID );

 if( failures ) {
  std::cout << failures << " check(s) FAILED" << std::endl;
  ret = 1;
  }
 else
  std::cout << "All tests passed!!" << std::endl;

 return( ret );
 }

/*--------------------------------------------------------------------------*/
/*------------------------- End File test.cpp ------------------------------*/
/*--------------------------------------------------------------------------*/
