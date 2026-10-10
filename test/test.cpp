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
 * The linking rows of a UCBlock are checked on two nodes: the zones of the
 * reserves (a node in no zone with one zone or with several, the zone
 * vector absent, the scaling of a unit in no zone), the constant terms of
 * the NetworkBlock (that of a group prevails, the bus case keeps them), and
 * the sign of the dual values of the node injection, reserve, inertia and
 * pollutant rows on instances whose duals are known.
 *
 * The changes of MaxPower, Availability, InitialPower and InitUpDownTime of
 * a ThermalUnitBlock and of a NuclearUnitBlock after the generation of the
 * abstract representation are compared, in every formulation, with the
 * unit read afresh from the changed data [see test_power_limits.cpp]: the
 * model that a :MILPSolver writes, row group by row group, and the values
 * of the :MILPSolver and of the dynamic programming Solver. By default the
 * reduced matrix of that test is run; with the option --power-limits the
 * tester runs only that test, the further options going to it (--full for
 * the whole matrix, see test_power_limits()).
 *
 * The test needs nothing but the core SMS++, so that the CI of UCBlock
 * builds this module alone; the signs of the duals and the test of the
 * power limits need a :MILPSolver in the factory, and are skipped when the
 * build has none. The netCDF groups are written in a dataset in memory, and
 * the test of the power limits writes the models in files of the working
 * directory, which it removes unless a case fails.
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
#include <numbers>
#include <cmath>
#include <cstdio>
#include <functional>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <vector>

#include <netcdf>

#include "BatteryUnitBlock.h"
#include "DCNetworkBlock.h"
#include "ACNetworkBlock.h"
#include "DesignNetworkBlock.h"
#include "BlockSolverConfig.h"
#include "CDASolver.h"
#include "ECNetworkBlock.h"
#include "FakeSolver.h"
#include "FRealObjective.h"
#include "DQuadFunction.h"
#include "HydroSystemUnitBlock.h"
#include "HydroUnitBlock.h"
#include "IntermittentUnitBlock.h"
#include "LinearFunction.h"
#include "NuclearUnitBlock.h"
#include "NuclearUnitExtDPSolver.h"
#include "SlackUnitBlock.h"
#include "ThermalUnitBlock.h"
#include "ThermalUnitDPSolver.h"
#include "ThermalUnitExtDPSolver.h"
#include "UCBlock.h"

// where the configuration files of test_power_limits() are
#ifndef UCBLOCK_TEST_DIR
#define UCBLOCK_TEST_DIR "."
#endif

// the changes of the power limits after the generation, compared with the
// unit read afresh [see test_power_limits.cpp]: the number of failed cases
int test_power_limits( const std::string & dir ,
                       const std::vector< std::string > & args );

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
 std::vector< double > prho;  // PrimaryRho
 std::vector< double > srho;  // SecondaryRho
 std::vector< double > prc;   // PrimarySpinningReserveCost
 std::vector< double > src;   // SecondarySpinningReserveCost
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
 if( ! d.prho.empty() )
  put( g , "PrimaryRho" , NI , d.prho );
 if( ! d.srho.empty() )
  put( g , "SecondaryRho" , NI , d.srho );
 if( ! d.prc.empty() )
  put( g , "PrimarySpinningReserveCost" , NI , d.prc );
 if( ! d.src.empty() )
  put( g , "SecondarySpinningReserveCost" , NI , d.src );
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
/*----------------- THE ROWS OF THE NETWORK FORMULATIONS -------------------*/
/*--------------------------------------------------------------------------*/
/* The rows of the three formulations of DCNetworkBlock are checked on
 * points that satisfy them or not, on networks with several components of
 * the lines with susceptance, a ReferenceNode that is not the lowest node
 * of its component, a singular matrix of the susceptances, and a change of
 * the demand after the rows are generated (for an ACNetworkBlock as well);
 * the Objective of an OTSNetworkBlock has to have the cost of each flow
 * once, and the dual prices of the flow limits one sign convention whether
 * the limits are bounds or rows. No Solver is needed.
 *
 * A line of a network written by write_net(): start and end node,
 * susceptance (0 for an HVDC line), flow bounds and efficiency. */

struct NetLine {
 Index s;
 Index e;
 double sus;
 double mn;
 double mx;
 double eta = 1;
 };

/*--------------------------------------------------------------------------*/
/// a NetworkBlock of the given type in a new group, with the given lines,
/// demand and (if positive) ReferenceNode

static netCDF::NcGroup write_net( const std::string & type , Index N ,
				  const std::vector< NetLine > & lines ,
				  const std::vector< double > & demand ,
				  Index ref = 0 )
{
 auto g = new_group( "NET" , true );
 g.putAtt( "type" , type );
 auto ND = g.addDim( "NumberNodes" , N );
 auto L = g.addDim( "NumberLines" , lines.size() );
 if( ref > 0 )  // a dimension of size 0 would be an unlimited one
  g.addDim( "ReferenceNode" , ref );
 std::vector< int > s , e;
 std::vector< double > sus , mn , mx , eta;
 for( const auto & l : lines ) {
  s.push_back( l.s );
  e.push_back( l.e );
  sus.push_back( l.sus );
  mn.push_back( l.mn );
  mx.push_back( l.mx );
  eta.push_back( l.eta );
  }
 put_int( g , "StartLine" , L , s );
 put_int( g , "EndLine" , L , e );
 put( g , "LineSusceptance" , L , sus );
 put( g , "MinPowerFlow" , L , mn );
 put( g , "MaxPowerFlow" , L , mx );
 put( g , "Efficiency" , L , eta );
 put( g , "ActiveDemand" , ND , demand );
 return( g );
 }

/*--------------------------------------------------------------------------*/
/// the DCNetworkBlock of the group, with its abstract representation in
/// the formulation wf (0 PTDF, 1 CYCLE, 2 KIRCHHOFF)

static DCNetworkBlock * new_net( const netCDF::NcGroup & g , int wf )
{
 auto nb = dynamic_cast< DCNetworkBlock * >( Block::new_Block( g ) );
 if( ! nb )
  throw( std::logic_error( "new_net: not a DCNetworkBlock" ) );
 SimpleConfiguration< int > f( wf );
 nb->generate_abstract_variables( & f );
 nb->generate_abstract_constraints();
 nb->generate_objective();
 return( nb );
 }

/*--------------------------------------------------------------------------*/
/// sets node injections, flows and (if any) voltage angles of the network

static void set_net_point( DCNetworkBlock * nb ,
			   const std::vector< double > & S ,
			   const std::vector< double > & F ,
			   const std::vector< double > & theta = {} )
{
 auto inj = nb->get_node_injection( 0 );
 for( Index n = 0 ; n < S.size() ; ++n )
  inj[ n ].set_value( S[ n ] );
 auto & f = nb->get_power_flow();
 for( Index l = 0 ; l < F.size() ; ++l )
  const_cast< ColVariable & >( f[ l ] ).set_value( F[ l ] );
 if( auto a = nb->get_static_variable_v< ColVariable >( "voltage_angle" ) )
  for( Index n = 0 ; n < a->size() ; ++n )
   ( *a )[ n ].set_value( n < theta.size() ? theta[ n ] : 0 );
 }

/*--------------------------------------------------------------------------*/

static const std::vector< std::string > FORMULATION = {
 "PTDF" , "CYCLE" , "KIRCHHOFF" };

/*--------------------------------------------------------------------------*/
/* The balance of every component of the lines with nonzero susceptance, in
 * every formulation: a point that moves power from a node to another of a
 * different component, with no line carrying it, satisfies all the rows of
 * the flows and the overall balance, and it has to be cut by the balance of
 * the components. The networks: two components joined by an HVDC line of
 * capacity 0, two components and no HVDC line, a component and a node with
 * no line, an HVDC line and two nodes with no line; the point that only
 * moves power along the lines is feasible in all of them. */

static void test_net_components( void )
{
 const double B = 1000;
 struct Case {
  std::string what;
  Index N;
  std::vector< NetLine > lines;
  std::vector< double > demand;
  std::vector< double > S_teleport;  // power from a component to another
  std::vector< double > S_right;     // power where the demand is
  };
 const std::vector< Case > cases = {
  { "two components and an HVDC line of capacity 0" , 4 ,
    { { 0 , 1 , 10 , -B , B } , { 2 , 3 , 10 , -B , B } ,
      { 1 , 3 , 0 , 0 , 0 } } ,
    { 0 , 0 , 100 , 0 } , { 100 , 0 , 0 , 0 } , { 0 , 0 , 100 , 0 } } ,
  { "two components and no HVDC line" , 4 ,
    { { 0 , 1 , 10 , -B , B } , { 2 , 3 , 10 , -B , B } } ,
    { 0 , 0 , 100 , 0 } , { 100 , 0 , 0 , 0 } , { 0 , 0 , 100 , 0 } } ,
  { "a node with no line" , 3 , { { 0 , 1 , 10 , -B , B } } ,
    { 0 , 0 , 100 } , { 100 , 0 , 0 } , { 0 , 0 , 100 } } ,
  { "an HVDC line and two nodes with no line" , 4 ,
    { { 0 , 1 , 0 , -B , B } } ,
    { 0 , 0 , 100 , 0 } , { 0 , 0 , 0 , 100 } , { 0 , 0 , 100 , 0 } } };

 for( const auto & c : cases )
  for( int wf = 0 ; wf < 3 ; ++wf ) {
   const auto what = "components, " + FORMULATION[ wf ] + ", " + c.what;
   try {
    std::unique_ptr< DCNetworkBlock > nb( new_net(
		  write_net( "DCNetworkBlock" , c.N , c.lines , c.demand ) ,
		  wf ) );
    std::vector< double > F( c.lines.size() , 0 );
    set_net_point( nb.get() , c.S_right , F );
    check( nb->is_feasible() , what + ": the balanced point is cut" );
    set_net_point( nb.get() , c.S_teleport , F );
    check( ! nb->is_feasible() , what + ": power moves with no line" );
    }
   catch( std::exception & e ) {
    check( false , what + ": throws " + e.what() );
    }
   }
 }

/*--------------------------------------------------------------------------*/
/* The reference nodes. A line 0 -> 1 with susceptance, an HVDC line 0 -> 2
 * and ReferenceNode = 1, which is not the lowest node of its component:
 * the flow 100 on the line from 0 to the demand at 1 is feasible in every
 * formulation, with the angles of KIRCHHOFF that fix to 0 that of node 1
 * and that of node 2, the only node of the other component; an angle of
 * node 2 other than 0 is not feasible, one reference angle being fixed in
 * each component. Then a ReferenceNode that is not a node, and the PTDF
 * matrix with its default Tikhonov coefficient, which is 0: the factor of
 * a single line is exactly -1. */

static void test_net_reference( void )
{
 const double B = 1000;
 const std::vector< NetLine > lines = { { 0 , 1 , 10 , -B , B } ,
					{ 0 , 2 , 0 , -B , B } };
 for( int wf = 0 ; wf < 3 ; ++wf ) {
  const auto what = "reference, " + FORMULATION[ wf ];
  try {
   std::unique_ptr< DCNetworkBlock > nb( new_net(
		 write_net( "DCNetworkBlock" , 3 , lines , { 0 , 100 , 0 } ,
			    1 ) , wf ) );
   set_net_point( nb.get() , { 100 , 0 , 0 } , { 100 , 0 } , { 10 , 0 , 0 } );
   check( nb->is_feasible() , what + ": ReferenceNode = 1 cuts the flow "
	  "of its own component" );
   if( wf == 2 ) {
    auto ra = nb->get_static_constraint_v< BoxConstraint >(
							 "reference_angle" );
    check( ra && ( ra->size() == 2 ) , what + ": not one reference angle "
	   "per component" );
    set_net_point( nb.get() , { 100 , 0 , 0 } , { 100 , 0 } ,
		   { 10 , 0 , 5 } );
    check( ! nb->is_feasible() , what + ": the angle of the other "
	   "component is free" );
    }
   }
  catch( std::exception & e ) {
   check( false , what + ": throws " + e.what() );
   }
  }

 // Block::new_Block() may report the exception and return nullptr
 bool thrown = false;
 try {
  auto b = Block::new_Block( write_net( "DCNetworkBlock" , 3 , lines ,
					{ 0 , 100 , 0 } , 3 ) );
  thrown = ! b;
  delete b;
  }
 catch( std::exception & ) { thrown = true; }
 check( thrown , "reference: ReferenceNode = NumberNodes is accepted" );

 std::unique_ptr< DCNetworkBlock > nb( dynamic_cast< DCNetworkBlock * >(
   Block::new_Block( write_net( "DCNetworkBlock" , 2 ,
				{ { 0 , 1 , 10 , -B , B } } , { 0 , 0 } ) ) ) );
 auto nd = static_cast< DCNetworkBlock::DCNetworkData * >(
					       nb->get_NetworkData() );
 auto ptdf = nd->get_PTDF( nd->get_DC_lines() );
 check( ( ptdf.cols() == 1 ) && close( ptdf.coeff( 0 , 0 ) , -1 , 1e-12 ) ,
	"reference: the PTDF factor of a single line is not -1" );
 }

/*--------------------------------------------------------------------------*/
/* The PTDF matrix: it cannot be computed when the matrix of the
 * susceptances of the nodes but the references is singular (susceptances
 * of opposite sign that cancel out), and the formulation throws instead of
 * using a zero matrix; with no line with susceptance it has no column,
 * also when asked with no argument. A network with hyperarcs and a line
 * with susceptance is not accepted. */

static void test_net_ptdf( void )
{
 const double B = 1000;
 bool thrown = false;
 try {
  std::unique_ptr< DCNetworkBlock > nb( new_net(
	   write_net( "DCNetworkBlock" , 3 ,
		      { { 0 , 1 , 1 , -B , B } , { 1 , 2 , 1 , -B , B } ,
			{ 0 , 2 , -0.5 , -B , B } } , { 0 , 0 , 100 } ) , 0 ) );
  }
 catch( std::logic_error & ) { thrown = true; }
 catch( std::exception & ) {}
 check( thrown , "PTDF: a singular matrix of the susceptances does not "
	"throw" );

 try {
  std::unique_ptr< DCNetworkBlock > nb( dynamic_cast< DCNetworkBlock * >(
       Block::new_Block( write_net( "DCNetworkBlock" , 3 ,
				    { { 0 , 1 , 0 , -B , B } ,
				      { 1 , 2 , 0 , -B , B } } ,
				    { 0 , 0 , 0 } ) ) ) );
  auto nd = static_cast< DCNetworkBlock::DCNetworkData * >(
					       nb->get_NetworkData() );
  auto ptdf = nd->get_PTDF();
  check( ptdf.cols() == 0 , "PTDF: a pure HVDC network has a column" );
  }
 catch( std::exception & e ) {
  check( false , std::string( "PTDF: a pure HVDC network throws " ) +
	 e.what() );
  }

 auto g = new_group( "NET" , true );
 g.putAtt( "type" , "DCNetworkBlock" );
 auto ND = g.addDim( "NumberNodes" , 3 );
 auto L = g.addDim( "NumberLines" , 2 );
 auto BR = g.addDim( "NumberBranches" , 3 );
 put_int( g , "StartLine" , BR , { 0 , 0 , 1 } );
 put_int( g , "EndLine" , BR , { 1 , 2 , 2 } );
 put_int( g , "HyperArcID" , BR , { 0 , 0 , 1 } );
 put( g , "LineSusceptance" , L , { 0 , 5 } );
 put( g , "MaxPowerFlow" , L , { B , B } );
 put( g , "MinPowerFlow" , L , { -B , -B } );
 put( g , "ActiveDemand" , ND , { 0 , 0 , 0 } );
 thrown = false;
 try {
  auto b = Block::new_Block( g );
  thrown = ! b;
  delete b;
  }
 catch( std::exception & ) { thrown = true; }
 check( thrown , "PTDF: hyperarcs and a line with susceptance accepted" );
 }

/*--------------------------------------------------------------------------*/
/* A change of the demand after the abstract representation is generated,
 * in every formulation, at a node at the end of an HVDC line: a pure HVDC
 * network, whose nodes all have a balance row, and a line with susceptance
 * 0 -> 1 with an HVDC line 0 -> 2; the point that meets the new demand has
 * to be feasible. Then an ACNetworkBlock, whose own active balance row of
 * the node changed has to hold the new demand. */

static void test_net_demand_change( void )
{
 const double B = 1000;
 for( int wf = 0 ; wf < 3 ; ++wf ) {
  auto what = "demand change, " + FORMULATION[ wf ] + ", pure HVDC";
  try {
   std::unique_ptr< DCNetworkBlock > nb( new_net(
	      write_net( "DCNetworkBlock" , 2 , { { 0 , 1 , 0 , -B , B } } ,
			 { 0 , 100 } ) , wf ) );
   std::vector< double > d = { 50 };
   nb->set_active_demand( d.cbegin() , Subset( { 1 } ) , true );
   set_net_point( nb.get() , { 50 , 0 } , { 50 } );
   check( nb->is_feasible() , what + ": the new demand is not met" );
   set_net_point( nb.get() , { 100 , 0 } , { 100 } );
   check( ! nb->is_feasible() , what + ": the old demand is met" );
   }
  catch( std::exception & e ) {
   check( false , what + ": throws " + e.what() );
   }

  what = "demand change, " + FORMULATION[ wf ] + ", mixed";
  try {
   std::unique_ptr< DCNetworkBlock > nb( new_net(
	      write_net( "DCNetworkBlock" , 3 ,
			 { { 0 , 1 , 10 , -B , B } , { 0 , 2 , 0 , -B , B } } ,
			 { 0 , 100 , 0 } ) , wf ) );
   std::vector< double > d = { 10 };
   nb->set_active_demand( d.cbegin() , Range( 2 , 3 ) );
   set_net_point( nb.get() , { 110 , 0 , 0 } , { 100 , 10 } ,
		  { 0 , -10 , 0 } );
   check( nb->is_feasible() , what + ": the new demand is not met" );
   }
  catch( std::exception & e ) {
   check( false , what + ": throws " + e.what() );
   }
  }

 // an ACNetworkBlock of three nodes on a line
 auto g = new_group( "AC" , true );
 g.putAtt( "type" , "ACNetworkBlock" );
 auto N = g.addDim( "NumberNodes" , 3 );
 auto L = g.addDim( "NumberLines" , 2 );
 put_int( g , "StartLine" , L , { 0 , 1 } );
 put_int( g , "EndLine" , L , { 1 , 2 } );
 put( g , "MinPowerFlow" , L , { -100 , -100 } );
 put( g , "MaxPowerFlow" , L , { 100 , 100 } );
 put( g , "LineSusceptance" , L , { 1 , 2 } );
 put( g , "LineRATEA" , L , { 100 , 100 } );
 put( g , "LineReactance" , L , { 0.1 , 0.2 } );
 put( g , "LineResistance" , L , { 0.01 , 0.02 } );
 put( g , "LineMinAngle" , L , { -30 , -30 } );
 put( g , "LineMaxAngle" , L , { 30 , 30 } );
 put( g , "NodeConductance" , N , { 0 , 0 , 0 } );
 put( g , "NodeSusceptance" , N , { 0 , 0 , 0 } );
 put( g , "NodeMaxVoltage" , N , { 1.1 , 1.1 , 1.1 } );
 put( g , "NodeMinVoltage" , N , { 0.9 , 0.9 , 0.9 } );
 put( g , "ReactiveDemand" , N , { 0.5 , 0.5 , 0.5 } );
 put( g , "ActiveDemand" , N , { 1 , 2 , 3 } );
 const std::string what = "demand change, ACNetworkBlock";
 try {
  std::unique_ptr< ACNetworkBlock > ab( dynamic_cast< ACNetworkBlock * >(
						 Block::new_Block( g ) ) );
  for( Index n = 0 ; n < 3 ; ++n ) {
   ab->set_min_node_injection( -10 , n );
   ab->set_max_node_injection( 10 , n );
   ab->set_min_reactive_node_injection( -10 , n , 0 );
   ab->set_max_reactive_node_injection( 10 , n , 0 );
   }
  ab->generate_abstract_variables();
  ab->generate_abstract_constraints();
  std::vector< double > d = { 7 };
  ab->set_active_demand( d.cbegin() , Subset( { 1 } ) , true );
  auto rows = ab->get_static_constraint_v< FRowConstraint >(
						   "AC_power_flow_injection" );
  check( rows && ( rows->size() == 6 ) &&
	 ( ( *rows )[ 1 ].get_lhs() == -7 ) &&
	 ( ( *rows )[ 1 ].get_rhs() == -7 ) &&
	 ( ( *rows )[ 0 ].get_rhs() == -1 ) ,
	 what + ": the active balance of the node does not hold the new "
	 "demand" );
  }
 catch( std::exception & e ) {
  check( false , what + ": throws " + e.what() );
  }
 }

/*--------------------------------------------------------------------------*/
/* The cost of the flows of an OTSNetworkBlock: each priced line has its
 * auxiliary Variable once in the Objective, with the cost of the data as
 * the coefficient, also if the Objective is asked for twice. */

static void test_net_OTS_cost( void )
{
 const double B = 1000;
 auto g = write_net( "OTSNetworkBlock" , 3 ,
		     { { 0 , 1 , 10 , -B , B } , { 1 , 2 , 5 , -B , B } } ,
		     { 0 , 0 , 1 } );
 put( g , "NetworkCost" , g.getDim( "NumberLines" ) , { 2 , 3 } );
 const std::string what = "OTS cost";
 try {
  std::unique_ptr< DCNetworkBlock > nb( dynamic_cast< DCNetworkBlock * >(
						 Block::new_Block( g ) ) );
  nb->generate_abstract_variables();
  nb->generate_abstract_constraints();
  nb->generate_objective();
  nb->generate_objective();
  auto lf = static_cast< LinearFunction * >( static_cast< FRealObjective * >(
				       nb->get_objective() )->get_function() );
  const auto & aux = nb->get_auxiliary_variable();
  const std::vector< double > cost = { 2 , 3 };
  for( Index l = 0 ; l < 2 ; ++l ) {
   int n = 0;
   double c = 0;
   for( const auto & [ v , coeff ] : lf->get_v_var() )
    if( v == & aux[ l ] ) {
     ++n;
     c += coeff;
     }
   check( ( n == 1 ) && ( c == cost[ l ] ) , what + ": line " +
	  std::to_string( l ) + " has its cost " + std::to_string( n ) +
	  " times, " + str( c ) + " in all" );
   }
  }
 catch( std::exception & e ) {
  check( false , what + ": throws " + e.what() );
  }
 }

/*--------------------------------------------------------------------------*/
/* A negative NetworkCost leaves V_l >= | F_l | unbounded above in the
 * Objective: it is refused when read and by both set_network_cost(), while
 * a zero cost is accepted. */

static void test_net_negative_cost( void )
{
 const double B = 1000;
 auto make = [ & ]( const std::vector< double > & cost ) {
  auto g = write_net( "DCNetworkBlock" , 3 ,
		      { { 0 , 1 , 10 , -B , B } , { 1 , 2 , 5 , -B , B } } ,
		      { 0 , 0 , 1 } );
  put( g , "NetworkCost" , g.getDim( "NumberLines" ) , cost );
  DCNetworkBlock * nb = nullptr;
  try {
   nb = dynamic_cast< DCNetworkBlock * >( Block::new_Block( g ) );
   }
  catch( std::exception & ) {}
  return( nb );
  };

 auto nb = make( { 2 , -1 } );
 check( ! nb , "network cost -1: accepted when read" );
 delete nb;

 nb = make( { 2 , 0 } );
 if( ! nb ) {
  check( false , "network cost 0: refused when read" );
  return;
  }
 nb->generate_abstract_variables();
 nb->generate_abstract_constraints();
 nb->generate_objective();
 auto throws = [ & ]( auto f ) {
  try {
   f();
   }
  catch( std::invalid_argument & ) {
   return( true );
   }
  return( false );
  };
 std::vector< double > val = { -1 };
 check( throws( [ & ]() { nb->set_network_cost( val.cbegin() ,
						Range( 1 , 2 ) ); } ) ,
	"network cost: set_network_cost( range ) accepts -1" );
 check( throws( [ & ]() { nb->set_network_cost( val.cbegin() ,
						Subset( { 1 } ) ); } ) ,
	"network cost: set_network_cost( subset ) accepts -1" );
 val = { 4 };
 check( ! throws( [ & ]() { nb->set_network_cost( val.cbegin() ,
						  Range( 1 , 2 ) ); } ) ,
	"network cost: set_network_cost( range ) refuses 4" );
 delete nb;
 }

/*--------------------------------------------------------------------------*/
/* The dual prices of the flow limits: for a line with a design variable
 * and a nonzero minimum flow the two sides are two rows, whose duals
 * (>= 0 for the upper one, <= 0 for the lower one) add up to the price of
 * the line, with the sign that the dual of the bound of a line with no
 * design variable has. */

static void test_net_dual_prices( void )
{
 auto g = write_net( "DCNetworkBlock" , 3 ,
		     { { 0 , 1 , 0 , -5 , 5 } , { 1 , 2 , 0 , 0 , 6 } ,
		       { 0 , 2 , 0 , -7 , 7 } } , { 0 , 0 , 1 } );
 const std::string what = "dual prices";
 try {
  std::unique_ptr< DCNetworkBlock > nb( dynamic_cast< DCNetworkBlock * >(
						 Block::new_Block( g ) ) );
  std::vector< ColVariable > x( 2 );
  Subset which = { 0 , 1 };
  nb->set_design_variables( & x , & which );
  SimpleConfiguration< int > kirchhoff( 2 );
  nb->generate_abstract_variables( & kirchhoff );
  nb->generate_abstract_constraints();
  auto up = nb->get_static_constraint_v< FRowConstraint >(
						 "Power_flow_limit_design" );
  auto lo = nb->get_static_constraint_v< FRowConstraint >(
					     "Power_flow_limit_design_min" );
  auto box = nb->get_static_constraint_v< BoxConstraint >(
							"Power_flow_limit" );
  if( ! ( up && lo && box && ( up->size() == 2 ) && ( lo->size() == 1 ) ) ) {
   check( false , what + ": the rows of the flow limits" );
   return;
   }
  ( *up )[ 0 ].set_dual( 0 );    // line 0: lower row active
  ( *lo )[ 0 ].set_dual( -3 );
  ( *up )[ 1 ].set_dual( 4 );    // line 1: upper row active
  ( *box )[ 1 ].set_dual( 0 );   // its lower half, a bound
  ( *box )[ 2 ].set_dual( -2 );  // line 2: no design, lower bound active
  std::vector< double > dp;
  nb->get_dual_prices( dp );
  check( ( dp.size() == 3 ) && ( dp[ 0 ] == -3 ) && ( dp[ 1 ] == 4 ) &&
	 ( dp[ 2 ] == -2 ) , what + ": " + str( dp[ 0 ] ) + " , " +
	 str( dp[ 1 ] ) + " , " + str( dp[ 2 ] ) + " instead of -3 , 4 , -2" );
  }
 catch( std::exception & e ) {
  check( false , what + ": throws " + e.what() );
  }
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
/*--------------------------------- HYDRO ----------------------------------*/
/*--------------------------------------------------------------------------*/
/* The data of a HydroUnitBlock as the test writes it: an empty vector is
 * not written. The matrices indexed by instant and arc are stored [ t ][ l ],
 * those indexed by reservoir and instant [ n ][ t ], InertiaPower
 * [ l ][ t ]. */

struct HUData {
 Index T = 1;                     // time horizon
 Index R = 1;                     // number of reservoirs
 Index A = 1;                     // number of arcs
 std::vector< int > start;        // StartArc
 std::vector< int > end;          // EndArc
 std::vector< int > up;           // UphillFlow
 std::vector< int > dn;           // DownhillFlow
 std::vector< double > minF;      // MinFlow
 std::vector< double > maxF;      // MaxFlow
 std::vector< double > minP;      // MinPower
 std::vector< double > maxP;      // MaxPower
 std::vector< double > rampUp;    // DeltaRampUp
 std::vector< double > lin;       // LinearTerm
 std::vector< double > cnst;      // ConstantTerm
 std::vector< double > inertia;   // InertiaPower
 std::vector< double > initV;     // InitialVolumetric
 std::vector< double > maxV;      // MaxVolumetric
 std::vector< double > inflows;   // Inflows
 };

/*--------------------------------------------------------------------------*/

static void write_HU( netCDF::NcGroup g , const HUData & d )
{
 g.putAtt( "type" , "HydroUnitBlock" );
 auto T = g.addDim( "TimeHorizon" , d.T );
 auto R = g.addDim( "NumberReservoirs" , d.R );
 auto A = g.addDim( "NumberArcs" , d.A );
 auto put_i = [ & ]( const std::string & n , const std::vector< int > & v ) {
  if( ! v.empty() )
   put_int( g , n , A , v );
  };
 auto put_ta = [ & ]( const std::string & n ,
		      const std::vector< double > & v ) {
  if( ! v.empty() )
   put( g , n , { T , A } , v );
  };
 auto put_rt = [ & ]( const std::string & n ,
		      const std::vector< double > & v ) {
  if( ! v.empty() )
   put( g , n , { R , T } , v );
  };
 put_i( "StartArc" , d.start );
 put_i( "EndArc" , d.end );
 put_i( "UphillFlow" , d.up );
 put_i( "DownhillFlow" , d.dn );
 put_ta( "MinFlow" , d.minF );
 put_ta( "MaxFlow" , d.maxF );
 put_ta( "MinPower" , d.minP );
 put_ta( "MaxPower" , d.maxP );
 put_ta( "DeltaRampUp" , d.rampUp );
 if( ! d.lin.empty() )
  put( g , "LinearTerm" , A , d.lin );
 if( ! d.cnst.empty() )
  put( g , "ConstantTerm" , A , d.cnst );
 if( ! d.inertia.empty() )
  put( g , "InertiaPower" , { A , T } , d.inertia );
 if( ! d.initV.empty() )
  put( g , "InitialVolumetric" , R , d.initV );
 put_rt( "MaxVolumetric" , d.maxV );
 put_rt( "Inflows" , d.inflows );
 }

/*--------------------------------------------------------------------------*/
/// the HydroUnitBlock of the data, nullptr (and the message) if it throws

static HydroUnitBlock * new_HU( const HUData & d , std::string & why )
{
 auto g = new_group( "HU" , true );
 write_HU( g , d );
 try {
  return( dynamic_cast< HydroUnitBlock * >( Block::new_Block( g ) ) );
  }
 catch( std::exception & e ) {
  why = e.what();
  }
 return( nullptr );
 }

/*--------------------------------------------------------------------------*/
/// a HydroUnitBlock whose data are a turbine from reservoir 0 to the sea

static HUData turbine_HU( Index T )
{
 HUData d;
 d.T = T;
 d.start = { 0 };
 d.end = { 1 };
 d.maxF.assign( T , 10 );
 d.maxP.assign( T , 10 );
 d.lin = { 1 };
 d.maxV.assign( T , 100 );
 return( d );
 }

/*--------------------------------------------------------------------------*/
/// whether the point of flows f, powers p and volumes v (reservoir by
/// reservoir) of the HydroUnitBlock is feasible, absolute tolerance 1e-9

static bool hydro_point( HydroUnitBlock * b ,
			 const std::vector< std::vector< double > > & f ,
			 const std::vector< std::vector< double > > & p ,
			 const std::vector< std::vector< double > > & v )
{
 for( Index l = 0 ; l < f.size() ; ++l )
  for( Index t = 0 ; t < f[ l ].size() ; ++t ) {
   b->get_flow_rate( l , t )->set_value( f[ l ][ t ] );
   b->get_active_power( l , t )->set_value( p[ l ][ t ] );
   }
 for( Index n = 0 ; n < v.size() ; ++n )
  for( Index t = 0 ; t < v[ n ].size() ; ++t )
   b->get_volume( n , t )->set_value( v[ n ][ t ] );
 SimpleConfiguration< std::pair< double , int > > tol( { 1e-9 , 0 } );
 return( b->is_feasible( false , & tol ) );
 }

/*--------------------------------------------------------------------------*/
/* The delays of the water balance at the borders of the horizon. With
 * UphillFlow = -1 the flow at instant s is withdrawn from the reservoir at
 * s - 1: the flow at 1 at 0, the flow at 0 never (it left before the
 * horizon). With UphillFlow = 1 the flow at 2 would be withdrawn at 3, after
 * the horizon, hence never. A negative delay must not drop the withdrawals
 * from the balance altogether. */

static void test_hydro_delays( void )
{
 std::string why;
 auto d = turbine_HU( 3 );
 d.initV = { 5 };
 d.up = { -1 };
 auto b = new_HU( d , why );
 if( ! b ) {
  check( false , "hydro, UphillFlow -1: throws " + why );
  return;
  }
 generate_all( b );
 check( hydro_point( b , { { 0 , 1 , 0 } } , { { 0 , 1 , 0 } } ,
		     { { 4 , 4 , 4 } } ) ,
	"hydro, UphillFlow -1: the flow at 1 is not withdrawn at 0" );
 check( ! hydro_point( b , { { 0 , 1 , 0 } } , { { 0 , 1 , 0 } } ,
		       { { 5 , 5 , 5 } } ) ,
	"hydro, UphillFlow -1: the flow at 1 is never withdrawn" );
 check( hydro_point( b , { { 1 , 0 , 0 } } , { { 1 , 0 , 0 } } ,
		     { { 5 , 5 , 5 } } ) ,
	"hydro, UphillFlow -1: the flow at 0, which left before the "
	"horizon, is withdrawn" );
 delete b;

 d.up = { 1 };
 b = new_HU( d , why );
 if( ! b ) {
  check( false , "hydro, UphillFlow 1: throws " + why );
  return;
  }
 generate_all( b );
 check( hydro_point( b , { { 1 , 0 , 0 } } , { { 1 , 0 , 0 } } ,
		     { { 5 , 4 , 4 } } ) ,
	"hydro, UphillFlow 1: the flow at 0 is not withdrawn at 1" );
 check( hydro_point( b , { { 0 , 0 , 1 } } , { { 0 , 0 , 1 } } ,
		     { { 5 , 5 , 5 } } ) ,
	"hydro, UphillFlow 1: the flow at 2, due after the horizon, is "
	"withdrawn" );
 delete b;
 }

/*--------------------------------------------------------------------------*/
/* A negative InitialVolumetric is the cyclic closure, which shapes the rows
 * of instant 0: once they are generated the sign of the value cannot
 * change, while a change that keeps it moves the right-hand side. */

static void test_hydro_initial_volume( void )
{
 std::string why;
 auto d = turbine_HU( 2 );
 d.initV = { 5 };
 auto b = new_HU( d , why );
 if( ! b ) {
  check( false , "hydro initial volume: throws " + why );
  return;
  }
 generate_all( b );
 std::vector< double > val = { 3 };
 b->set_initial_volume( val.begin() , Range( 0 , 1 ) );
 check( hydro_point( b , { { 0 , 0 } } , { { 0 , 0 } } , { { 3 , 3 } } ) ,
	"hydro initial volume: the change to 3 is not in the balance" );
 bool thrown = false;
 try {
  val = { -1 };
  b->set_initial_volume( val.begin() , Range( 0 , 1 ) );
  }
 catch( std::invalid_argument & ) {
  thrown = true;
  }
 check( thrown , "hydro initial volume: the cyclic closure set after the "
	"rows are generated" );
 delete b;
 }

/*--------------------------------------------------------------------------*/
/* The cyclic closure (InitialVolumetric < 0): the data check asks only the
 * total inflow to fit the total flows (an inflow above the capacity of the
 * arcs at one instant is buffered by the reservoir), also for one reservoir
 * without StartArc and EndArc; a change of the inflow at instant 0, or of
 * the (negative) initial volume, keeps the initial volume out of the
 * balance of instant 0, as the generation does. */

static void test_hydro_cyclic( void )
{
 std::string why;
 auto d = turbine_HU( 3 );
 d.initV = { -1 };
 d.inflows = { 9 , 3 , 3 };      // within MaxFlow 10 at every instant
 auto b = new_HU( d , why );
 if( ! b ) {
  check( false , "hydro cyclic: throws " + why );
  return;
  }
 generate_all( b );
 check( hydro_point( b , { { 5 , 5 , 5 } } , { { 5 , 5 , 5 } } ,
		     { { 4 , 2 , 0 } } ) ,
	"hydro cyclic: the point of the data is not feasible" );

 // set_inflow( range ) at instant 0: the row is v_0 - v_2 + f_0 = 12
 std::vector< double > val = { 12 };
 b->set_inflow( val.begin() , Range( 0 , 1 ) );
 check( hydro_point( b , { { 6 , 6 , 6 } } , { { 6 , 6 , 6 } } ,
		     { { 6 , 3 , 0 } } ) ,
	"hydro cyclic, set_inflow( range ): the initial volume is in the "
	"balance of instant 0" );

 // set_inflow( subset ) at instant 0: back to 9
 val = { 9 };
 b->set_inflow( val.begin() , Subset( { 0 } ) , true );
 check( hydro_point( b , { { 5 , 5 , 5 } } , { { 5 , 5 , 5 } } ,
		     { { 4 , 2 , 0 } } ) ,
	"hydro cyclic, set_inflow( subset ): the initial volume is in the "
	"balance of instant 0" );

 // a new negative initial volume keeps the closure: nothing moves
 val = { -7 };
 b->set_initial_volume( val.begin() , Range( 0 , 1 ) );
 check( hydro_point( b , { { 5 , 5 , 5 } } , { { 5 , 5 , 5 } } ,
		     { { 4 , 2 , 0 } } ) ,
	"hydro cyclic, set_initial_volume( range ): the initial volume is in "
	"the balance of instant 0" );
 val = { -3 };
 b->set_initial_volume( val.begin() , Subset( { 0 } ) , true );
 check( hydro_point( b , { { 5 , 5 , 5 } } , { { 5 , 5 , 5 } } ,
		     { { 4 , 2 , 0 } } ) ,
	"hydro cyclic, set_initial_volume( subset ): the initial volume is "
	"in the balance of instant 0" );
 delete b;

 // an inflow above the capacity of the arcs at one instant is buffered
 d.inflows = { 15 , 0 , 0 };
 b = new_HU( d , why );
 if( b ) {
  generate_all( b );
  check( hydro_point( b , { { 5 , 5 , 5 } } , { { 5 , 5 , 5 } } ,
		      { { 10 , 5 , 0 } } ) ,
	 "hydro cyclic, inflow 15 at one instant: the buffered point is not "
	 "feasible" );
  delete b;
  }
 else
  check( false , "hydro cyclic, inflow 15 > MaxFlow 10 at one instant: "
	 "refused " + why );

 // a total inflow equal to the total capacity is accepted, one more no
 d.inflows = { 15 , 15 , 0 };
 b = new_HU( d , why );
 check( b != nullptr , "hydro cyclic, total inflow 30 = capacity: "
	"throws " + why );
 delete b;
 d.inflows = { 15 , 15 , 1 };
 b = new_HU( d , why );
 check( ! b , "hydro cyclic, total inflow 31 > capacity 30: accepted" );
 delete b;

 // the same with one reservoir and no StartArc / EndArc
 d.start.clear();
 d.end.clear();
 b = new_HU( d , why );
 check( ! b , "hydro cyclic, one reservoir without arcs, total inflow "
	"31 > capacity 30: accepted" );
 delete b;
 d.inflows = { 30 , 0 , 0 };
 b = new_HU( d , why );
 check( b != nullptr , "hydro cyclic, one reservoir without arcs, total "
	"inflow 30: throws " + why );
 delete b;

 // an arc that enters the reservoir: its minimum flow must be absorbed
 // too; reservoir 1 receives arc 0 (MinFlow 4) and is emptied by arc 1
 HUData e;
 e.T = 2;
 e.R = 2;
 e.A = 2;
 e.start = { 0 , 1 };
 e.end = { 1 , 2 };
 e.minF = { 4 , 0 , 4 , 0 };
 e.maxF = { 10 , 5 , 10 , 5 };
 e.maxP = { 10 , 5 , 10 , 5 };
 e.lin = { 1 , 1 };
 e.initV = { 50 , -1 };
 e.maxV.assign( 4 , 100 );
 e.inflows = { 0 , 0 , 2 , 0 };  // reservoir 1: 2 + 8 in, at most 10 out
 b = new_HU( e , why );
 check( b != nullptr , "hydro cyclic, entering arc, net inflow 10 = "
	"capacity: throws " + why );
 delete b;
 e.inflows = { 0 , 0 , 3 , 0 };
 b = new_HU( e , why );
 check( ! b , "hydro cyclic, entering arc, net inflow 11 > capacity 10: "
	"accepted" );
 delete b;
 }

/*--------------------------------------------------------------------------*/
/* The data that the HydroUnitBlock rejects, and the limit cases it accepts:
 * the arcs, the kind of an arc (turbine, pump, idle) and the LinearTerm of a
 * turbine. */

static void test_hydro_data_checks( void )
{
 std::string why;
 auto expect = [ & ]( const HUData & d , bool ok , const std::string & what ) {
  why.clear();
  auto b = new_HU( d , why );
  if( ok )
   check( b != nullptr , "hydro data, " + what + ": throws " + why );
  else  // Block::new_Block() reports the exception and returns nullptr
   check( ! b , "hydro data, " + what + ": accepted" );
  delete b;
  };

 auto d = turbine_HU( 2 );
 d.lin.clear();
 expect( d , false , "turbine without LinearTerm" );

 d.maxF = { 0 , 0 };
 expect( d , true , "idle arc without LinearTerm" );

 d.minF = { -10 , -10 };
 d.minP = { -10 , -10 };
 d.maxP = { 0 , 0 };
 expect( d , true , "pump without LinearTerm" );

 d = turbine_HU( 2 );
 d.minF = { -5 , -5 };
 d.maxF = { 5 , 5 };
 expect( d , false , "MinFlow < 0 < MaxFlow" );

 d.minF = { 0 , -10 };
 d.maxF = { 10 , 0 };
 expect( d , false , "a turbine at 0 and a pump at 1" );

 d.minF = { 0 , 0 };
 d.maxF = { 10 , 0 };
 expect( d , true , "a turbine at 0 and idle at 1" );

 d = turbine_HU( 2 );
 d.R = 2;
 d.maxV.assign( 4 , 100 );
 expect( d , true , "two reservoirs" );

 d.start.clear();
 d.end.clear();
 expect( d , false , "two reservoirs without StartArc and EndArc" );

 d.R = 1;
 d.maxV.assign( 2 , 100 );
 expect( d , true , "one reservoir without StartArc and EndArc" );

 d.start = { 0 };
 expect( d , false , "StartArc without EndArc" );

 d = turbine_HU( 2 );
 d.R = 2;
 d.maxV.assign( 4 , 100 );
 d.start = { 2 };
 expect( d , false , "StartArc == NumberReservoirs" );

 d.start = { 0 };
 d.end = { 3 };
 expect( d , false , "EndArc > NumberReservoirs" );

 d.end = { 0 };
 expect( d , false , "a self-loop" );
 }

/*--------------------------------------------------------------------------*/
/* The flow-to-power row of a single-piece turbine with no constant term is
 * an equality only if a spillage arc can take the water that the turbine
 * would pass without producing: same StartArc, EndArc and delays, enough
 * MaxFlow at every instant, a turbine with MinFlow 0 and no ramp
 * constraint. The point "the turbine passes 5 and produces nothing" is
 * feasible exactly when the row is an inequality. */

static void test_hydro_spillage( void )
{
 // arc 0: turbine 0 -> 1; arc 1: spillage 0 -> 1; volumes 10 and 0
 auto base = []( void ) {
  HUData d;
  d.R = 2;
  d.A = 2;
  d.start = { 0 , 0 };
  d.end = { 1 , 1 };
  d.maxF = { 10 , 10 };
  d.maxP = { 10 , 0 };
  d.lin = { 1 , 0 };
  d.cnst = { 0 , 0 };
  d.initV = { 10 , 0 };
  d.maxV = { 100 , 100 };
  return( d );
  };

 // whether the row of the turbine is an equality
 auto equality = [ & ]( const HUData & d , const std::string & what ) {
  std::string why;
  auto b = new_HU( d , why );
  if( ! b ) {
   check( false , "hydro spillage, " + what + ": throws " + why );
   return( false );
   }
  generate_all( b );
  std::vector< std::vector< double > > f( d.A , { 0 } );
  std::vector< std::vector< double > > p( d.A , { 0 } );
  f[ 0 ][ 0 ] = 5;
  const bool eq = ! hydro_point( b , f , p , { { 5 } , { 5 } } );
  delete b;
  return( eq );
  };

 auto d = base();
 check( equality( d , "compatible" ) ,
	"hydro spillage: a compatible spillage arc gives no equality" );

 d.end = { 1 , 2 };
 check( ! equality( d , "other EndArc" ) ,
	"hydro spillage: equality with a spillage arc to another reservoir" );

 d = base();
 d.maxF = { 10 , 9 };
 check( ! equality( d , "small spillage" ) ,
	"hydro spillage: equality with a spillage arc smaller than the "
	"turbine" );

 d = base();
 d.dn = { 0 , 1 };
 check( ! equality( d , "other delay" ) ,
	"hydro spillage: equality with a spillage arc with another delay" );

 d = base();
 d.minF = { 1 , 0 };
 check( ! equality( d , "MinFlow" ) ,
	"hydro spillage: equality on a turbine with MinFlow > 0" );

 d = base();
 d.rampUp = { 100 , 100 };
 check( ! equality( d , "ramps" ) ,
	"hydro spillage: equality on a turbine with ramp constraints" );

 d = base();
 d.cnst = { 1 , 0 };
 check( ! equality( d , "ConstantTerm" ) ,
	"hydro spillage: equality on a turbine with a constant term" );

 // two turbines on one spillage arc, whose MaxFlow has to cover both
 d = base();
 d.A = 3;
 d.start = { 0 , 0 , 0 };
 d.end = { 1 , 1 , 1 };
 d.maxF = { 10 , 10 , 15 };
 d.maxP = { 10 , 10 , 0 };
 d.lin = { 1 , 1 , 0 };
 d.cnst = { 0 , 0 , 0 };
 check( ! equality( d , "two turbines" ) ,
	"hydro spillage: equality with a spillage arc smaller than the two "
	"turbines" );
 d.maxF = { 10 , 10 , 20 };
 check( equality( d , "two turbines, large spillage" ) ,
	"hydro spillage: no equality with a spillage arc as large as the "
	"two turbines" );
 }

/*--------------------------------------------------------------------------*/
/// a UCBlock with one node, one instant, one inertia zone requiring 5 and
/// the given unit (written by unit) as UnitBlock_0

static UCBlock * inertia_UC( std::function< void( netCDF::NcGroup ) > unit )
{
 auto g = new_group( "UC" , true );
 g.putAtt( "type" , "UCBlock" );
 auto T = g.addDim( "TimeHorizon" , 1 );
 g.addDim( "NumberUnits" , 1 );
 auto one = g.addDim( "NumberNodes" , 1 );
 auto Z = g.addDim( "NumberInertiaZones" , 1 );
 put( g , "ActivePowerDemand" , { one , T } , { 3 } );
 put( g , "InertiaDemand" , { Z , T } , { 5 } );
 unit( g.addGroup( "UnitBlock_0" ) );
 return( dynamic_cast< UCBlock * >( Block::new_Block( g ) ) );
 }

/*--------------------------------------------------------------------------*/
/// the coefficient of v in the inertia row of the UCBlock, NaN if absent

static double inertia_coef( UCBlock * b , const ColVariable * v )
{
 auto lf = static_cast< LinearFunction * >(
	    b->get_inertia_demand_constraints()[ 0 ][ 0 ].get_function() );
 const auto i = lf->is_active( v );
 if( i >= lf->get_num_active_var() )
  return( std::nan( "" ) );
 return( lf->get_coefficient( i ) );
 }

/*--------------------------------------------------------------------------*/
/* A change of the InertiaPower of a HydroUnitBlock has to reach the inertia
 * row of the UCBlock, both when the HydroUnitBlock is a unit of the UCBlock
 * and when it is inside a HydroSystemUnitBlock (with no
 * PolyhedralFunctionBlock, which is allowed); an InertiaPower given after
 * the rows were generated without it cannot be added to them, and is
 * refused. */

static void test_hydro_inertia( void )
{
 auto d = turbine_HU( 1 );
 d.inertia = { 1 };
 d.initV = { 50 };

 for( int sys = 0 ; sys < 2 ; ++sys ) {
  const std::string what = sys ? "hydro inertia, in a HydroSystemUnitBlock"
                               : "hydro inertia";
  auto b = inertia_UC( [ & ]( netCDF::NcGroup u ) {
   if( sys ) {
    u.putAtt( "type" , "HydroSystemUnitBlock" );
    u.addDim( "TimeHorizon" , 1 );
    u.addDim( "NumberHydroUnits" , 1 );
    write_HU( u.addGroup( "HydroUnitBlock_0" ) , d );
    }
   else
    write_HU( u , d );
   } );
  if( ! b ) {
   check( false , what + ": the UCBlock is not read" );
   continue;
   }
  HydroUnitBlock * hub = nullptr;
  if( sys ) {
   auto hsub = dynamic_cast< HydroSystemUnitBlock * >(
						  b->get_unit_block( 0 ) );
   if( hsub ) {
    check( hsub->get_polyhedral_function_block() == nullptr ,
	   what + ": a PolyhedralFunctionBlock out of nothing" );
    hub = hsub->get_hydro_unit_block( 0 );
    }
   }
  else
   hub = dynamic_cast< HydroUnitBlock * >( b->get_unit_block( 0 ) );
  if( ! hub ) {
   check( false , what + ": the HydroUnitBlock is not there" );
   delete b;
   continue;
   }
  try {
   generate_all( b );
   }
  catch( std::exception & e ) {
   check( false , what + ": the abstract representation throws " +
	  e.what() );
   delete b;
   continue;
   }
  auto fs = new FakeSolver();
  b->register_Solver( fs );
  const auto p = hub->get_active_power( 0 , 0 );
  check( close( inertia_coef( b , p ) , 1 ) ,
	 what + ": coefficient " + str( inertia_coef( b , p ) ) +
	 " instead of 1" );
  std::vector< double > val = { 2 };
  hub->set_inertia_power( val.begin() , Range( 0 , 1 ) );
  check( close( inertia_coef( b , p ) , 2 ) ,
	 what + ": coefficient " + str( inertia_coef( b , p ) ) +
	 " after set_inertia_power( 2 )" );
  b->unregister_Solver( fs , true );
  delete b;
  }

 d.inertia.clear();
 auto b = inertia_UC( [ & ]( netCDF::NcGroup u ) { write_HU( u , d ); } );
 if( ! b ) {
  check( false , "hydro inertia added: the UCBlock is not read" );
  return;
  }
 generate_all( b );
 auto fs = new FakeSolver();
 b->register_Solver( fs );
 bool thrown = false;
 try {
  std::vector< double > val = { 2 };
  static_cast< HydroUnitBlock * >( b->get_unit_block( 0 ) )->
   set_inertia_power( val.begin() , Range( 0 , 1 ) );
  }
 catch( std::logic_error & ) {
  thrown = true;
  }
 check( thrown , "hydro inertia added after the rows: accepted" );
 b->unregister_Solver( fs , true );
 delete b;
 }

/*--------------------------------------------------------------------------*/
/* HydroSystemUnitBlock: the active Variable of the PolyhedralFunction are
 * the final volumes, unit after unit and reservoir after reservoir; a
 * concave PolyhedralFunction is refused. */

static void test_hydro_system( void )
{
 auto d0 = turbine_HU( 2 );  // two reservoirs
 d0.R = 2;
 d0.maxV.assign( 4 , 100 );
 auto d1 = turbine_HU( 2 );  // one reservoir

 for( int concave = 0 ; concave < 2 ; ++concave ) {
  auto g = new_group( "HS" , true );
  g.putAtt( "type" , "HydroSystemUnitBlock" );
  g.addDim( "TimeHorizon" , 2 );
  g.addDim( "NumberHydroUnits" , 2 );
  write_HU( g.addGroup( "HydroUnitBlock_0" ) , d0 );
  write_HU( g.addGroup( "HydroUnitBlock_1" ) , d1 );
  auto pg = g.addGroup( "PolyhedralFunctionBlock" );
  pg.putAtt( "type" , "PolyhedralFunctionBlock" );
  auto nv = pg.addDim( "PolyFunction_NumVar" , 3 );
  auto nr = pg.addDim( "PolyFunction_NumRow" , 1 );
  put( pg , "PolyFunction_A" , { nr , nv } , { -1 , -2 , -3 } );
  put( pg , "PolyFunction_b" , nr , { 0 } );
  if( concave )
   pg.addDim( "PolyFunction_sign" , 0 );

  HydroSystemUnitBlock * b = nullptr;
  std::string why;
  try {
   b = dynamic_cast< HydroSystemUnitBlock * >( Block::new_Block( g ) );
   }
  catch( std::invalid_argument & e ) {
   why = e.what();
   }
  if( concave ) {
   check( ! b , "hydro system: a concave future cost accepted" );
   delete b;
   continue;
   }
  if( ! b ) {
   check( false , "hydro system: not read " + why );
   continue;
   }
  generate_all( b );
  auto & pf = b->get_polyhedral_function_block()->get_PolyhedralFunction();
  const auto h0 = b->get_hydro_unit_block( 0 );
  const auto h1 = b->get_hydro_unit_block( 1 );
  check( ( pf.get_num_active_var() == 3 ) &&
	 ( pf.get_active_var( 0 ) == h0->get_volume( 0 , 1 ) ) &&
	 ( pf.get_active_var( 1 ) == h0->get_volume( 1 , 1 ) ) &&
	 ( pf.get_active_var( 2 ) == h1->get_volume( 0 , 1 ) ) ,
	 "hydro system: the arguments of the future cost are not the final "
	 "volumes in order" );
  #ifndef NDEBUG
  bool thrown = false;
  try {
   b->get_hydro_unit_block( 2 );
   }
  catch( std::invalid_argument & ) {
   thrown = true;
   }
  check( thrown , "hydro system: get_hydro_unit_block( 2 ) of 2 accepted" );
  #endif
  delete b;
  }
 }

/*--------------------------------------------------------------------------*/
/*------------------- STORAGE, INTERMITTENT AND SLACK UNITS ----------------*/
/*--------------------------------------------------------------------------*/
/// the coefficient of v in the row c, NaN if v is not in it

static double row_coef( const FRowConstraint & c , const Variable * v )
{
 auto lf = dynamic_cast< const LinearFunction * >( c.get_function() );
 if( ! lf )
  return( std::numeric_limits< double >::quiet_NaN() );
 const auto i = lf->is_active( v );
 if( i >= lf->get_num_active_var() )
  return( std::numeric_limits< double >::quiet_NaN() );
 return( lf->get_coefficient( i ) );
 }

/*--------------------------------------------------------------------------*/
/// the coefficient of v in the (linear) Objective of the unit b, NaN if v is
/// not in it

static double obj_coef( Block * b , const Variable * v )
{
 auto lf = dynamic_cast< const LinearFunction * >(
	    static_cast< FRealObjective * >( b->get_objective() )
	                                                    ->get_function() );
 if( ! lf )
  return( std::numeric_limits< double >::quiet_NaN() );
 const auto i = lf->is_active( v );
 if( i >= lf->get_num_active_var() )
  return( std::numeric_limits< double >::quiet_NaN() );
 return( lf->get_coefficient( i ) );
 }

/*--------------------------------------------------------------------------*/
/// a BatteryUnitBlock of horizon T with the given data, read from a group;
/// nullptr (and the message in why) if reading it throws

static BatteryUnitBlock * new_BU( Index T ,
				  const std::vector< std::pair< std::string ,
				   std::vector< double > > > & data ,
				  std::string * why = nullptr )
{
 auto g = new_group( "BU" , true );
 g.putAtt( "type" , "BatteryUnitBlock" );
 auto TH = g.addDim( "TimeHorizon" , T );
 for( const auto & d : data )
  if( d.second.size() == 1 )
   put( g , d.first , d.second[ 0 ] );
  else
   put( g , d.first , TH , d.second );
 try {
  return( dynamic_cast< BatteryUnitBlock * >( Block::new_Block( g ) ) );
  }
 catch( std::exception & e ) {
  if( why )
   *why = e.what();
  return( nullptr );
  }
 }

/*--------------------------------------------------------------------------*/
/* The storage balance with a single efficiency given: the absent one is 1,
 * so that with only "ExtractingBatteryRho" the intake fills the battery
 * with coefficient 1, as with "StoringBatteryRho" = 1 written out, and with
 * only "StoringBatteryRho" the outtake empties it with coefficient 1; at
 * instant 0 (InitialStorage on the right-hand side) as at the others. With
 * the two efficiencies equal the pair is folded and the common value
 * multiplies the active power. */

static void test_battery_storage_balance( void )
{
 const std::vector< std::pair< std::string , std::vector< double > > > base =
  { { "MinStorage" , { 0 , 0 } } , { "MaxStorage" , { 20 , 20 } } ,
    { "MaxPower" , { 10 , 10 } } , { "MinPower" , { -10 , -10 } } ,
    { "InitialStorage" , { 5 } } };

 struct Case {
  std::string what;
  std::vector< double > sbr , ebr;
  double in , out;  // expected coefficients of intake and outtake
  };
 const std::vector< Case > cases = {
  { "only ExtractingBatteryRho" , { } , { 1.25 , 1.25 } , -1 , 1.25 } ,
  { "both, StoringBatteryRho = 1" , { 1 , 1 } , { 1.25 , 1.25 } , -1 , 1.25 } ,
  { "only StoringBatteryRho" , { 0.8 , 0.8 } , { } , -0.8 , 1 } };

 for( const auto & c : cases ) {
  auto data = base;
  if( ! c.sbr.empty() )
   data.push_back( { "StoringBatteryRho" , c.sbr } );
  if( ! c.ebr.empty() )
   data.push_back( { "ExtractingBatteryRho" , c.ebr } );
  auto bu = new_BU( 2 , data );
  if( ! bu ) {
   check( false , "battery balance, " + c.what + ": not read" );
   continue;
   }
  generate_all( bu );
  const auto & bal = bu->get_storage_balance_constraints();
  auto & il = bu->get_intake_level();
  auto & ol = bu->get_outtake_level();
  auto & sl = bu->get_storage_level();
  check( ( bal.size() == 2 ) && ( il.size() == 2 ) && ( ol.size() == 2 ) ,
	 "battery balance, " + c.what + ": intake and outtake not split" );
  if( ( bal.size() == 2 ) && ( il.size() == 2 ) && ( ol.size() == 2 ) )
   for( Index t = 0 ; t < 2 ; ++t ) {
    check( close( row_coef( bal[ t ] , & il[ t ] ) , c.in ) &&
	   close( row_coef( bal[ t ] , & ol[ t ] ) , c.out ) ,
	   "battery balance, " + c.what + ", t = " + std::to_string( t ) +
	   ": intake " + str( row_coef( bal[ t ] , & il[ t ] ) ) +
	   " and outtake " + str( row_coef( bal[ t ] , & ol[ t ] ) ) +
	   " instead of " + str( c.in ) + " and " + str( c.out ) );
    check( close( row_coef( bal[ t ] , & sl[ t ] ) , 1 ) ,
	   "battery balance, " + c.what + ": level coefficient" );
    }
  if( ! bal.empty() )
   check( ( bal[ 0 ].get_lhs() == 5 ) && ( bal[ 0 ].get_rhs() == 5 ) ,
	  "battery balance, " + c.what + ": instant 0 has InitialStorage " +
	  str( bal[ 0 ].get_rhs() ) + " on the right-hand side" );
  delete bu;
  }

 // the two efficiencies equal: folded onto the active power
 auto data = base;
 data.push_back( { "StoringBatteryRho" , { 2 , 2 } } );
 data.push_back( { "ExtractingBatteryRho" , { 2 , 2 } } );
 if( auto bu = new_BU( 2 , data ) ) {
  generate_all( bu );
  const auto & bal = bu->get_storage_balance_constraints();
  check( bu->get_intake_level().empty() ,
	 "battery balance, equal efficiencies: the pair is split" );
  auto p = bu->get_active_power( 0 );
  if( p && ( bal.size() == 2 ) )
   check( close( row_coef( bal[ 1 ] , & p[ 1 ] ) , 2 ) ,
	  "battery balance, equal efficiencies: active power coefficient " +
	  str( row_coef( bal[ 1 ] , & p[ 1 ] ) ) + " instead of 2" );
  delete bu;
  }
 else
  check( false , "battery balance, equal efficiencies: not read" );
 }

/*--------------------------------------------------------------------------*/
/* The sign of InitialStorage selects the cyclic balance at instant 0: once
 * the Constraints are generated a change of sign throws, a change that
 * keeps the sign only moves the right-hand side (a change between two
 * negative values changes nothing), and before generation any value is
 * accepted and decides the form of the row. */

static void test_battery_initial_storage_sign( void )
{
 const std::vector< std::pair< std::string , std::vector< double > > > base =
  { { "MinStorage" , { 0 , 0 , 0 } } , { "MaxStorage" , { 20 , 20 , 20 } } ,
    { "MaxPower" , { 10 , 10 , 10 } } };

 auto with_init = [ & ]( double v0 ) {
  auto data = base;
  data.push_back( { "InitialStorage" , { v0 } } );
  return( new_BU( 3 , data ) );
  };

 auto throws = [ & ]( BatteryUnitBlock * bu , double v ) {
  std::vector< double > val = { v };
  try {
   bu->set_initial_storage( val.cbegin() , Range( 0 , 1 ) );
   }
  catch( std::logic_error & ) {
   return( true );
   }
  return( false );
  };

 if( auto bu = with_init( 5 ) ) {
  generate_all( bu );
  check( throws( bu , -1 ) ,
	 "battery: InitialStorage from 5 to -1 after generation accepted" );
  check( ! throws( bu , 3 ) , "battery: InitialStorage from 5 to 3 throws" );
  const auto & bal = bu->get_storage_balance_constraints();
  check( ( bu->get_initial_storage() == 3 ) && ( bal[ 0 ].get_rhs() == 3 ) ,
	 "battery: InitialStorage 3 gives the right-hand side " +
	 str( bal[ 0 ].get_rhs() ) );
  // the Subset version too
  std::vector< double > val = { -2 };
  bool thrown = false;
  try {
   bu->set_initial_storage( val.cbegin() , Subset( { 0 } ) );
   }
  catch( std::logic_error & ) {
   thrown = true;
   }
  check( thrown , "battery: set_initial_storage( Subset ) accepts a change "
	 "of sign after generation" );
  delete bu;
  }
 else
  check( false , "battery, InitialStorage 5: not read" );

 if( auto bu = with_init( -1 ) ) {
  generate_all( bu );
  auto & sl = bu->get_storage_level();
  const auto & bal = bu->get_storage_balance_constraints();
  check( ! std::isnan( row_coef( bal[ 0 ] , & sl[ 2 ] ) ) ,
	 "battery, cyclic: the row of instant 0 lacks the last level" );
  check( throws( bu , 4 ) ,
	 "battery: InitialStorage from -1 to 4 after generation accepted" );
  check( ! throws( bu , -2 ) ,
	 "battery: InitialStorage from -1 to -2 throws" );
  check( bal[ 0 ].get_rhs() == 0 ,
	 "battery, cyclic: right-hand side " + str( bal[ 0 ].get_rhs() ) );
  delete bu;
  }
 else
  check( false , "battery, InitialStorage -1: not read" );

 if( auto bu = with_init( 5 ) ) {  // before generation: anything goes
  check( ! throws( bu , -1 ) ,
	 "battery: a change of sign before generation throws" );
  generate_all( bu );
  auto & sl = bu->get_storage_level();
  check( ! std::isnan( row_coef( bu->get_storage_balance_constraints()[ 0 ] ,
				 & sl[ 2 ] ) ) ,
	 "battery: made cyclic before generation, the row of instant 0 "
	 "lacks the last level" );
  delete bu;
  }
 }

/*--------------------------------------------------------------------------*/
/* The Cost of a battery set by an unsorted Subset: each value goes to the
 * index it is given with, in the data and in the Objective, where the price
 * on the net injection is - Cost times the active power. */

static void test_battery_set_cost_unsorted( void )
{
 auto bu = new_BU( 4 , { { "MinStorage" , { 0 , 0 , 0 , 0 } } ,
			 { "MaxStorage" , { 9 , 9 , 9 , 9 } } ,
			 { "MaxPower" , { 5 , 5 , 5 , 5 } } ,
			 { "Cost" , { 1 , 2 , 3 , 4 } } } );
 if( ! bu ) {
  check( false , "battery set_cost: not read" );
  return;
  }
 generate_all( bu );
 std::vector< double > val = { 7 , 9 };
 bu->set_cost( val.cbegin() , Subset( { 3 , 1 } ) , false );
 check( bu->get_cost() == std::vector< double >( { 1 , 9 , 3 , 7 } ) ,
	"battery set_cost( { 3 , 1 } , { 7 , 9 } ) gives Cost " +
	str( bu->get_cost()[ 1 ] ) + " at 1 and " + str( bu->get_cost()[ 3 ] ) +
	" at 3" );
 auto p = bu->get_active_power( 0 );
 check( close( obj_coef( bu , & p[ 3 ] ) , -7 ) &&
	close( obj_coef( bu , & p[ 1 ] ) , -9 ) ,
	"battery set_cost: Objective coefficients " +
	str( obj_coef( bu , & p[ 1 ] ) ) + " and " +
	str( obj_coef( bu , & p[ 3 ] ) ) + " instead of -9 and -7" );
 delete bu;
 }

/*--------------------------------------------------------------------------*/
/* A battery with a positive MinPower is refused when intake and outtake are
 * split (the intake would have a negative upper bound), accepted when they
 * are folded; a store with an unbounded MinStorage is accepted. */

static void test_battery_data_checks( void )
{
 std::string why;
 auto bu = new_BU( 2 , { { "MinStorage" , { 0 , 0 } } ,
			 { "MaxStorage" , { 9 , 9 } } ,
			 { "MaxPower" , { 5 , 5 } } , { "MinPower" , { 0 , 1 } } ,
			 { "StoringBatteryRho" , { 0.8 , 0.8 } } } , & why );
 check( ! bu , "battery: split with MinPower 1 accepted" );
 delete bu;

 bu = new_BU( 2 , { { "MinStorage" , { 0 , 0 } } ,
		    { "MaxStorage" , { 9 , 9 } } ,
		    { "MaxPower" , { 5 , 5 } } , { "MinPower" , { 0 , 1 } } } ,
	      & why );
 check( bu != nullptr , "battery: folded with MinPower 1 refused: " + why );
 delete bu;

 bu = new_BU( 2 , { { "MinStorage" , { -INF , -INF } } ,
		    { "MaxStorage" , { 9 , 9 } } ,
		    { "MaxPower" , { 5 , 5 } } ,
		    { "StoringBatteryRho" , { 0.8 , 0.8 } } } , & why );
 check( bu != nullptr , "battery: MinStorage -inf refused: " + why );
 if( bu ) {
  generate_all( bu );
  delete bu;
  }
 }

/*--------------------------------------------------------------------------*/
/* The converter fence of a battery with no design Variable: there only if
 * "ConverterMaxPower" is in the data (or the converter is designed), one
 * row on intake plus outtake or two on the active power, kappa times the
 * converter power on the right-hand side and non-binding where the latter
 * is 0; set_kappa() moves it; the default is not written back. */

static void test_battery_converter( void )
{
 const std::vector< std::pair< std::string , std::vector< double > > > base =
  { { "MinStorage" , { 0 , 0 } } , { "MaxStorage" , { 20 , 20 } } ,
    { "MaxPower" , { 10 , 10 } } };

 // no datum, no row
 if( auto bu = new_BU( 2 , base ) ) {
  generate_all( bu );
  check( bu->get_converter_bounds().num_elements() == 0 ,
	 "battery: converter rows with no ConverterMaxPower" );
  auto g = new_group( "BUs" );
  bu->serialize( g );
  check( g.getVar( "ConverterMaxPower" ).isNull() ,
	 "battery: the default ConverterMaxPower is written" );
  delete bu;
  }

 // split, the second converter power 0
 auto data = base;
 data.push_back( { "StoringBatteryRho" , { 0.8 , 0.8 } } );
 data.push_back( { "ConverterMaxPower" , { 2 , 0 } } );
 if( auto bu = new_BU( 2 , data ) ) {
  generate_all( bu );
  const auto & cb = bu->get_converter_bounds();
  auto & il = bu->get_intake_level();
  auto & ol = bu->get_outtake_level();
  check( ( cb.shape()[ 0 ] == 1 ) && ( cb.shape()[ 1 ] == 2 ) ,
	 "battery, split: not one converter row per instant" );
  if( cb.num_elements() == 2 ) {
   check( close( cb[ 0 ][ 0 ].get_rhs() , 2 ) &&
	  close( row_coef( cb[ 0 ][ 0 ] , & il[ 0 ] ) , 1 ) &&
	  close( row_coef( cb[ 0 ][ 0 ] , & ol[ 0 ] ) , 1 ) ,
	  "battery, split: converter row il + ol <= " +
	  str( cb[ 0 ][ 0 ].get_rhs() ) + " instead of 2" );
   check( cb[ 0 ][ 1 ].get_rhs() == INF ,
	  "battery, split: converter power 0 binds" );
   }
  auto g = new_group( "BUs" );
  bu->serialize( g );
  check( ! g.getVar( "ConverterMaxPower" ).isNull() ,
	 "battery: the ConverterMaxPower given is not written" );
  delete bu;
  }

 // folded, kappa 2, then set_kappa( 0.5 )
 data = base;
 data.push_back( { "ConverterMaxPower" , { 3 , 3 } } );
 data.push_back( { "Kappa" , { 2 } } );
 if( auto bu = new_BU( 2 , data ) ) {
  generate_all( bu );
  const auto & cb = bu->get_converter_bounds();
  auto p = bu->get_active_power( 0 );
  check( ( cb.shape()[ 0 ] == 2 ) && ( cb.shape()[ 1 ] == 2 ) ,
	 "battery, folded: not two converter rows per instant" );
  if( cb.num_elements() == 4 ) {
   check( close( cb[ 0 ][ 1 ].get_rhs() , 6 ) &&
	  close( cb[ 1 ][ 1 ].get_rhs() , 6 ) &&
	  close( row_coef( cb[ 0 ][ 1 ] , & p[ 1 ] ) , 1 ) &&
	  close( row_coef( cb[ 1 ][ 1 ] , & p[ 1 ] ) , -1 ) ,
	  "battery, folded: converter rows are not | p | <= 6" );
   bu->set_kappa( 0.5 );
   check( close( cb[ 0 ][ 0 ].get_rhs() , 1.5 ) &&
	  close( cb[ 1 ][ 1 ].get_rhs() , 1.5 ) ,
	  "battery, folded: set_kappa( 0.5 ) leaves the converter at " +
	  str( cb[ 0 ][ 0 ].get_rhs() ) );
   }
  delete bu;
  }

 // a designed converter on a battery with no design: the power is a
 // coefficient of the converter design Variable
 data = base;
 data.push_back( { "StoringBatteryRho" , { 0.8 , 0.8 } } );
 data.push_back( { "ConverterInvestmentCost" , { 1 } } );
 data.push_back( { "ConverterMaxCapacityDesign" , { 2 } } );
 if( auto bu = new_BU( 2 , data ) ) {
  generate_all( bu );
  const auto & cb = bu->get_converter_bounds();
  check( cb.num_elements() == 2 ,
	 "battery: a designed converter with no battery design has no row" );
  if( cb.num_elements() == 2 )
   check( close( row_coef( cb[ 0 ][ 0 ] , & bu->get_conv_design() ) , -10 ) &&
	  ( cb[ 0 ][ 0 ].get_rhs() == 0 ) ,
	  "battery: the designed converter row has coefficient " +
	  str( row_coef( cb[ 0 ][ 0 ] , & bu->get_conv_design() ) ) +
	  " instead of - kappa MaxPower = -10" );
  delete bu;
  }
 }

/*--------------------------------------------------------------------------*/
/* get_kappa_linearization() of a unit with a design Variable throws, also
 * when only the converter of a battery is designed; with none it does not
 * (the duals being all 0 it returns 0). */

static void test_kappa_linearization_design( void )
{
 auto throws = []( UnitBlock * ub ) {
  try {
   ub->get_kappa_linearization();
   }
  catch( std::logic_error & ) {
   return( true );
   }
  return( false );
  };

 const std::vector< std::pair< std::string , std::vector< double > > > base =
  { { "MinStorage" , { 0 , 0 } } , { "MaxStorage" , { 20 , 20 } } ,
    { "MaxPower" , { 10 , 10 } } };

 for( const std::string & cost : { "" , "BatteryInvestmentCost" ,
				   "ConverterInvestmentCost" } ) {
  auto data = base;
  if( ! cost.empty() )
   data.push_back( { cost , { 1 } } );
  auto bu = new_BU( 2 , data );
  if( ! bu ) {
   check( false , "kappa linearization, battery " + cost + ": not read" );
   continue;
   }
  generate_all( bu );
  check( throws( bu ) == ( ! cost.empty() ) ,
	 "kappa linearization, battery " + ( cost.empty() ? "with no design"
	 : "with " + cost ) + ": throws " + ( cost.empty() ? "" : "not " ) );
  delete bu;
  }

 for( bool design : { false , true } ) {
  auto g = new_group( "IU" , true );
  g.putAtt( "type" , "IntermittentUnitBlock" );
  auto T = g.addDim( "TimeHorizon" , 2 );
  put( g , "MaxPower" , T , { 7 , 3 } );
  if( design )
   put( g , "InvestmentCost" , 1 );
  auto ib = dynamic_cast< IntermittentUnitBlock * >( Block::new_Block( g ) );
  if( ! ib ) {
   check( false , "kappa linearization, intermittent: not read" );
   continue;
   }
  generate_all( ib );
  check( throws( ib ) == design ,
	 std::string( "kappa linearization, intermittent " ) +
	 ( design ? "with design: does not throw" : "with no design: throws" ) );
  delete ib;
  }
 }

/*--------------------------------------------------------------------------*/
/* The reserve rows (1)-(2) of an intermittent unit with a design Variable:
 * they exist, with gamma kappa MaxPower and kappa MinPower as coefficients
 * of the design Variable (0 where MinPower is), and follow set_kappa() and
 * set_maximum_power(); with Gamma = 0 there is no reserve and no row; with
 * no design the rows keep the data on the right-hand side. */

static void test_intermittent_design_reserve( void )
{
 auto make = []( double gamma , bool design ) {
  auto g = new_group( "IU" , true );
  g.putAtt( "type" , "IntermittentUnitBlock" );
  auto T = g.addDim( "TimeHorizon" , 2 );
  put( g , "MaxPower" , T , { 10 , 10 } );
  put( g , "MinPower" , T , { 1 , 0 } );
  put( g , "Gamma" , gamma );
  if( design ) {
   put( g , "InvestmentCost" , 1 );
   put( g , "MaxCapacityDesign" , 2 );
   }
  auto ib = dynamic_cast< IntermittentUnitBlock * >( Block::new_Block( g ) );
  if( ib ) {
   ib->set_reserve_vars( 3 );
   generate_all( ib );
   }
  return( ib );
  };

 if( auto ib = make( 0.5 , true ) ) {
  const auto & mn = ib->get_min_power_constraints();
  const auto & mx = ib->get_max_power_constraints();
  const auto & x = ib->get_const_design();
  check( ( mn.size() == 2 ) && ( mx.size() == 2 ) ,
	 "intermittent, design: no reserve rows" );
  if( ( mn.size() == 2 ) && ( mx.size() == 2 ) ) {
   auto p = ib->get_active_power( 0 );
   auto pr = ib->get_primary_spinning_reserve( 0 );
   check( close( row_coef( mx[ 0 ] , & x ) , -5 ) &&
	  ( mx[ 0 ].get_rhs() == 0 ) &&
	  close( row_coef( mx[ 0 ] , & p[ 0 ] ) , 0.5 ) &&
	  close( row_coef( mx[ 0 ] , & pr[ 0 ] ) , 1 ) ,
	  "intermittent, design: row (1) design coefficient " +
	  str( row_coef( mx[ 0 ] , & x ) ) + " instead of -5" );
   check( close( row_coef( mn[ 0 ] , & x ) , -1 ) &&
	  close( row_coef( mn[ 1 ] , & x ) , 0 ) &&
	  ( mn[ 0 ].get_lhs() == 0 ) &&
	  close( row_coef( mn[ 0 ] , & pr[ 0 ] ) , -1 ) ,
	  "intermittent, design: row (2) design coefficients " +
	  str( row_coef( mn[ 0 ] , & x ) ) + " and " +
	  str( row_coef( mn[ 1 ] , & x ) ) + " instead of -1 and 0" );
   ib->set_kappa( 2 );
   check( close( row_coef( mx[ 1 ] , & x ) , -10 ) &&
	  close( row_coef( mn[ 0 ] , & x ) , -2 ) ,
	  "intermittent, design: set_kappa( 2 ) gives " +
	  str( row_coef( mx[ 1 ] , & x ) ) + " in (1)" );
   std::vector< double > val = { 4 };
   ib->set_maximum_power( val.cbegin() , Subset( { 1 } ) );
   check( close( row_coef( mx[ 1 ] , & x ) , -4 ) ,
	  "intermittent, design: set_maximum_power gives " +
	  str( row_coef( mx[ 1 ] , & x ) ) + " in (1) instead of -4" );
   }
  delete ib;
  }

 if( auto ib = make( 0 , true ) ) {
  check( ib->get_min_power_constraints().empty() &&
	 ib->get_max_power_constraints().empty() &&
	 ( ib->get_primary_spinning_reserve( 0 ) == nullptr ) ,
	 "intermittent, design, Gamma 0: reserve rows or variables" );
  delete ib;
  }

 if( auto ib = make( 0.5 , false ) ) {
  const auto & mx = ib->get_max_power_constraints();
  const auto & mn = ib->get_min_power_constraints();
  check( ( mx.size() == 2 ) && close( mx[ 0 ].get_rhs() , 5 ) &&
	 ( mn.size() == 2 ) && close( mn[ 0 ].get_lhs() , 1 ) ,
	 "intermittent, no design: reserve rows not on the data" );
  delete ib;
  }
 }

/*--------------------------------------------------------------------------*/
// the :MILPSolver of the build and the value it finds (defined below)

static std::string milp_solver( void );
static double milp_value( Block * b , const std::string & sname );

/*--------------------------------------------------------------------------*/
/* A design variable whose two bounds are 1 is fixed to 1, i.e., the asset
 * is built and its investment cost paid, for a battery (battery and
 * converter design) and for an intermittent unit; with the bounds 0 and 1
 * it is not built. The Block alone, with no other cost, is solved by a
 * :MILPSolver. */

static void test_design_fixed_to_one( void )
{
 const auto sname = milp_solver();
 if( sname.empty() ) {
  std::cout << "design fixed to 1: no :MILPSolver, skipped" << std::endl;
  return;
  }

 for( const std::string & what : { "Battery" , "Converter" } )
  for( double lb : { 1.0 , 0.0 } ) {
   auto bu = new_BU( 2 , { { "MinStorage" , { 0 , 0 } } ,
			   { "MaxStorage" , { 20 , 20 } } ,
			   { "MaxPower" , { 10 , 10 } } ,
			   { what + "InvestmentCost" , { 7 } } ,
			   { what + "MinCapacityDesign" , { lb } } ,
			   { what + "MaxCapacityDesign" , { 1 } } } );
   if( ! bu ) {
    check( false , "design fixed to 1, " + what + ": not read" );
    continue;
    }
   generate_all( bu );
   const double v = milp_value( bu , sname );
   check( close( v , 7 * lb ) , "design fixed to 1, " + what +
	  ", bounds [ " + str( lb ) + " , 1 ]: value " + str( v ) +
	  " instead of " + str( 7 * lb ) );
   delete bu;
   }

 for( double lb : { 1.0 , 0.0 } )
  for( double ub : { 1.0 , -1.0 } ) {  // continuous and binary
   auto g = new_group( "IU" , true );
   g.putAtt( "type" , "IntermittentUnitBlock" );
   auto T = g.addDim( "TimeHorizon" , 2 );
   put( g , "MaxPower" , T , { 7 , 3 } );
   put( g , "InvestmentCost" , 5 );
   put( g , "MinCapacityDesign" , lb );
   put( g , "MaxCapacityDesign" , ub );
   auto ib = dynamic_cast< IntermittentUnitBlock * >( Block::new_Block( g ) );
   if( ! ib ) {
    check( false , "design fixed to 1, intermittent: not read" );
    continue;
    }
   generate_all( ib );
   const double v = milp_value( ib , sname );
   check( close( v , 5 * lb ) , "design fixed to 1, intermittent, bounds [ "
	  + str( lb ) + " , " + str( ub ) + " ]: value " + str( v ) +
	  " instead of " + str( 5 * lb ) );
   delete ib;
   }
 }

/*--------------------------------------------------------------------------*/
/* The storage and the converter of a battery as two assets, each with its
 * kappa. On a single bus over 4 instants with demand 0, 0, 150, 150, a
 * source of 100 at price 10 and one of 1000 at price 50, a battery moves
 * energy bought at 10 to the last two instants: q per instant saves 80 q,
 * and it needs a storage of 2 q and a converter of q. With the storage
 * priced 10 per unit of energy and the converter 25 per unit of power, the
 * continuous relaxation of the battery with the two design variables, i.e.,
 * of the whole problem, builds a storage of 100 and a converter of 50 for
 * 7000 - 4000 + 1000 + 1250 = 5250, and so do the two kappas at ( 100 , 50
 * ), while a single kappa, which resizes both, does best at 100 for 6500.
 * With a storage of 4 per unit, priced 40, the optimum is 25 units of
 * storage and 50 of converter, again 5250, and the single kappa does best at
 * 50 for 6250. The linearizations of the operational cost are checked at
 * ( 60 , 20 ), where the converter binds ( 0 , -80 ), and at ( 30 , 40 ),
 * where the storage does ( -40 , 0 ); at ( 0 , 0 ), where every fence of
 * the battery meets the other, they are checked to be a subgradient, as the
 * one of the single kappa is along the diagonal. A battery with no
 * "ConverterKappa" resizes its converter with set_kappa(), the datum is
 * written by serialize() only if given, and the kappa of the converter is
 * refused once the Constraint of a battery with no converter rows exist,
 * while before it gives them. */

static void test_battery_two_kappas( void )
{
 const auto sname = milp_solver();
 if( sname.empty() ) {
  std::cout << "battery, two kappas: no :MILPSolver, skipped" << std::endl;
  return;
  }

 const double NaN = std::numeric_limits< double >::quiet_NaN();

 // the UCBlock, with a battery of storage vmax per unit: with the design
 // variables (price cs per unit of storage, 25 per unit of converter) if
 // design, otherwise with the kappas ks and kc (kc not written if NaN)
 auto make = [ & ]( double vmax , double cs , bool design , double ks ,
		    double kc ) -> UCBlock * {
  auto g = new_group( "UC2K" , true );
  g.putAtt( "type" , "UCBlock" );
  auto th = g.addDim( "TimeHorizon" , 4 );
  g.addDim( "NumberUnits" , 3 );
  auto ng = g.addDim( "NumberElectricalGenerators" , 3 );
  auto nn = g.addDim( "NumberNodes" , 1 );
  put_int( g , "GeneratorNode" , ng , { 0 , 0 , 0 } );
  put( g , "ActivePowerDemand" , { nn , th } , { 0 , 0 , 150 , 150 } );
  for( Index u = 0 ; u < 2 ; ++u ) {
   auto s = g.addGroup( "UnitBlock_" + std::to_string( u ) );
   s.putAtt( "type" , "SlackUnitBlock" );
   put( s , "MaxPower" , u ? 1000 : 100 );
   put( s , "ActivePowerCost" , u ? 50 : 10 );
   }
  auto b = g.addGroup( "UnitBlock_2" );
  b.putAtt( "type" , "BatteryUnitBlock" );
  put( b , "MinStorage" , 0 );
  put( b , "MaxStorage" , vmax );
  put( b , "MaxPower" , 10 );
  put( b , "ConverterMaxPower" , 1 );
  if( design ) {
   put( b , "BatteryInvestmentCost" , cs );
   put( b , "BatteryMaxCapacityDesign" , 1000 );
   put( b , "ConverterInvestmentCost" , 25 );
   put( b , "ConverterMaxCapacityDesign" , 1000 );
   }
  else {
   put( b , "Kappa" , ks );
   if( ! std::isnan( kc ) )
    put( b , "ConverterKappa" , kc );
   }
  return( dynamic_cast< UCBlock * >( Block::new_Block( g ) ) );
  };

 auto battery = []( UCBlock * uc ) {
  return( static_cast< BatteryUnitBlock * >( uc->get_unit_block( 2 ) ) );
  };

 // the operational cost of uc, its dual values written in its Constraint;
 // the Solver stays registered, and receives the changes of kappa
 std::vector< std::pair< UCBlock * , BlockSolverConfig * > > registered;
 auto lp = [ & ]( UCBlock * uc ) {
  if( uc->get_registered_solvers().empty() ) {
   auto bsc = new BlockSolverConfig( 1 );
   bsc->add_ComputeConfig( std::string( sname ) , nullptr );
   bsc->apply( uc );
   bsc->clear();
   registered.emplace_back( uc , bsc );
   auto slv = uc->get_registered_solvers().front();
   const auto par = slv->int_par_str2idx( "intLogVerb" );
   if( par < Inf< Solver::idx_type >() )
    slv->set_par( par , 0 );
   }
  auto slv = static_cast< CDASolver * >(
                                     uc->get_registered_solvers().front() );
  if( ( slv->compute( false ) != Solver::kOK ) ||
      ( ! slv->has_dual_solution() ) )
   return( NaN );
  slv->get_var_solution();
  slv->get_dual_solution();
  return( slv->get_var_value() );
  };

 auto release = [ & ]( UCBlock * uc ) {
  for( auto & [ b , bsc ] : registered )
   if( b == uc ) {
    bsc->apply( uc );
    delete bsc;
    b = nullptr;
    }
  delete uc;
  };

 // the two instances: storage per unit, its price, the optimal kappas, and
 // the best single kappa with its value
 struct Case { double vmax , cs , ks , kc , k1 , v1; };
 for( const Case & c : { Case{ 1 , 10 , 100 , 50 , 100 , 6500 } ,
			 Case{ 4 , 40 , 25 , 50 , 50 , 6250 } } ) {
  const std::string what = "battery, two kappas, storage " + str( c.vmax ) +
                           " per unit";

  // the whole problem, with the design variables
  if( auto uc = make( c.vmax , c.cs , true , 1 , NaN ) ) {
   const double v = milp_value( uc , sname );
   check( close( v , 5250 ) , what + ": the design gives " + str( v ) +
	  " instead of 5250" );
   delete uc;
   }
  else
   check( false , what + ": the design UCBlock is not read" );

  // the two kappas at the optimum
  auto uc = make( c.vmax , c.cs , false , c.ks , c.kc );
  if( ! uc ) {
   check( false , what + ": the UCBlock is not read" );
   continue;
   }
  auto bu = battery( uc );
  check( ( ! bu->converter_follows_kappa() ) &&
	 close( bu->get_converter_kappa() , c.kc ) ,
	 what + ": the converter kappa read is " +
	 str( bu->get_converter_kappa() ) );
  const double v2 = lp( uc ) + c.cs * c.ks + 25 * c.kc;
  check( close( v2 , 5250 ) , what + ": the two kappas at ( " + str( c.ks )
	 + " , " + str( c.kc ) + " ) give " + str( v2 ) + " instead of 5250" );

  // the best single kappa, on a grid that contains it, is worse
  auto uc1 = make( c.vmax , c.cs , false , 0 , NaN );
  if( ! uc1 ) {
   check( false , what + ": the single kappa UCBlock is not read" );
   release( uc );
   continue;
   }
  auto b1 = battery( uc1 );
  double best = INF , kbest = -1;
  for( double k = 0 ; k <= 200 ; k += 5 ) {
   b1->set_kappa( k );
   check( close( b1->get_converter_kappa() , k ) , what +
	  ": the converter does not follow set_kappa( " + str( k ) + " )" );
   const double v = lp( uc1 ) + ( c.cs + 25 ) * k;
   if( v < best - 1e-9 ) {
    best = v;
    kbest = k;
    }
   }
  check( close( best , c.v1 ) && close( kbest , c.k1 ) ,
	 what + ": the single kappa does best at " + str( kbest ) + " with " +
	 str( best ) + " instead of " + str( c.k1 ) + " with " + str( c.v1 ) );

  // the linearization of the single kappa along the diagonal, from 0 too
  for( double k0 : { 0.0 , 40.0 } ) {
   b1->set_kappa( k0 );
   const double v0 = lp( uc1 );
   const double g = b1->get_kappa_linearization();
   for( double k : { 0.0 , 20.0 , 50.0 , 100.0 , 150.0 } ) {
    b1->set_kappa( k );
    const double v = lp( uc1 );
    check( v >= v0 + g * ( k - k0 ) - 1e-6 , what + ": the linearization "
	   + str( g ) + " of the single kappa at " + str( k0 ) +
	   " is above the value " + str( v ) + " at " + str( k ) );
    }
   }
  release( uc1 );
  release( uc );
  }

 // the coefficients of the two kappas where one of them binds, and their
 // validity at ( 0 , 0 )
 auto uc = make( 1 , 10 , false , 60 , 20 );
 if( ! uc ) {
  check( false , "battery, two kappas: the UCBlock is not read" );
  return;
  }
 auto bu = battery( uc );
 auto at = [ & ]( double ks , double kc , double & gs , double & gc ) {
  bu->set_kappa( ks );
  bu->set_converter_kappa( kc );
  const double v = lp( uc );
  gs = bu->get_kappa_linearization();
  gc = bu->get_converter_kappa_linearization();
  return( v );
  };
 double gs , gc;
 at( 60 , 20 , gs , gc );
 check( close( gs , 0 ) && close( gc , -80 ) , "battery, two kappas: at ( "
	"60 , 20 ) the linearization is ( " + str( gs ) + " , " + str( gc ) +
	" ) instead of ( 0 , -80 )" );
 at( 30 , 40 , gs , gc );
 check( close( gs , -40 ) && close( gc , 0 ) , "battery, two kappas: at ( "
	"30 , 40 ) the linearization is ( " + str( gs ) + " , " + str( gc ) +
	" ) instead of ( -40 , 0 )" );
 const double v0 = at( 0 , 0 , gs , gc );
 check( close( v0 , 7000 ) , "battery, two kappas: at ( 0 , 0 ) the value "
	"is " + str( v0 ) + " instead of 7000" );
 for( const auto & [ ks , kc ] : std::vector< std::pair< double , double > >
	{ { 100 , 50 } , { 30 , 40 } , { 60 , 20 } , { 10 , 80 } ,
	  { 0 , 50 } , { 100 , 0 } } ) {
  double hs , hc;
  const double v = at( ks , kc , hs , hc );
  check( v >= v0 + gs * ks + gc * kc - 1e-6 , "battery, two kappas: the "
	 "linearization ( " + str( gs ) + " , " + str( gc ) + " ) at ( 0 , 0 "
	 ") is above the value " + str( v ) + " at ( " + str( ks ) + " , " +
	 str( kc ) + " )" );
  }
 {
  auto g = new_group( "BUs" );
  bu->serialize( g );
  check( ! g.getVar( "ConverterKappa" ).isNull() ,
	 "battery, two kappas: the ConverterKappa given is not written" );
  }
 release( uc );

 // no "ConverterKappa": not written, and the converter follows the kappa
 const std::vector< std::pair< std::string , std::vector< double > > > base =
  { { "MinStorage" , { 0 , 0 } } , { "MaxStorage" , { 20 , 20 } } ,
    { "MaxPower" , { 10 , 10 } } };
 auto data = base;
 data.push_back( { "ConverterMaxPower" , { 3 , 3 } } );
 if( auto b = new_BU( 2 , data ) ) {
  generate_all( b );
  b->set_kappa( 2 );
  check( b->converter_follows_kappa() && close( b->get_converter_kappa() , 2 )
	 && close( b->get_converter_bounds()[ 0 ][ 0 ].get_rhs() , 6 ) ,
	 "battery, one kappa: the converter does not follow set_kappa( 2 )" );
  auto g = new_group( "BUs" );
  b->serialize( g );
  check( g.getVar( "ConverterKappa" ).isNull() ,
	 "battery, one kappa: a ConverterKappa is written" );
  // its own kappa from now on, and set_kappa() leaves it
  b->set_converter_kappa( 0.5 );
  b->set_kappa( 3 );
  check( close( b->get_converter_bounds()[ 0 ][ 0 ].get_rhs() , 1.5 ) ,
	 "battery: the converter kappa 0.5 gives the converter rows " +
	 str( b->get_converter_bounds()[ 0 ][ 0 ].get_rhs() ) );
  delete b;
  }

 // no converter rows: refused after the Constraint are generated, while
 // before they are given
 if( auto b = new_BU( 2 , base ) ) {
  generate_all( b );
  bool thrown = false;
  try {
   b->set_converter_kappa( 2 );
   }
  catch( std::logic_error & ) {
   thrown = true;
   }
  check( thrown , "battery: a converter kappa is accepted with no rows" );
  delete b;
  }
 if( auto b = new_BU( 2 , base ) ) {
  b->set_converter_kappa( 2 );
  generate_all( b );
  const auto & cb = b->get_converter_bounds();
  check( ( cb.num_elements() == 4 ) && close( cb[ 0 ][ 1 ].get_rhs() , 20 ) ,
	 "battery: a converter kappa 2 given before the Constraint does not "
	 "give the rows | p | <= 20" );
  delete b;
  }
 }

/*--------------------------------------------------------------------------*/
/* A battery of 3 modules (BatteryMinCapacityDesign = BatteryMaxCapacityDesign
 * = 3) with the binary variables (negative prices and a round-trip loss)
 * and the reserves: the binary rows (8) let the 3 modules (with a converter
 * of 30) charge and discharge 3 times the power of one, and the reserve
 * bounds (9) are those of the modules built, x_b kappa MaxPrimaryPower,
 * also after set_kappa(); with up to 1 module the rows are those of one
 * module. Integer designs with a minimum above 1 are accepted. */

static void test_battery_design_modules( void )
{
 auto make = []( double mn , double mx ) {
  auto bu = new_BU( 2 , { { "MinStorage" , { 0 , 0 } } ,
			  { "MaxStorage" , { 100 , 100 } } ,
			  { "MaxPower" , { 10 , 10 } } ,
			  { "MinPower" , { -10 , -10 } } ,
			  { "ConverterMaxPower" , { 30 , 30 } } ,
			  { "MaxPrimaryPower" , { 4 , 4 } } ,
			  { "StoringBatteryRho" , { 0.9 , 0.9 } } ,
			  { "Cost" , { -5 , 5 } } ,
			  { "BatteryInvestmentCost" , { 1 } } ,
			  { "BatteryMinCapacityDesign" , { mn } } ,
			  { "BatteryMaxCapacityDesign" , { mx } } } );
  if( bu ) {
   bu->set_reserve_vars( 1 );
   SimpleConfiguration< int > negative_prices( 1 );
   bu->generate_abstract_variables( & negative_prices );
   bu->generate_abstract_constraints();
   bu->generate_objective();
   }
  return( bu );
  };

 auto bu = make( 3 , 3 );
 if( ! bu ) {
  check( false , "battery modules: not read" );
  return;
  }
 const auto & x = bu->get_const_batt_design();
 const auto & lb0 = bu->get_primary_reserve_bounds();
 const auto & rd = bu->get_primary_reserve_design_rows();
 check( bu->get_max_intake_binary_constraints() != nullptr ,
	"battery modules: no binary rows" );
 check( ( lb0.size() == 2 ) && close( lb0[ 0 ].get_rhs() , 12 ) ,
	"battery modules: reserve bound " +
	( lb0.empty() ? std::string( "absent" ) : str( lb0[ 0 ].get_rhs() ) ) +
	" instead of 3 * 4" );
 check( ( rd.size() == 2 ) && close( row_coef( rd[ 1 ] , & x ) , -4 ) &&
	( rd[ 1 ].get_rhs() == 0 ) ,
	"battery modules: no reserve row of the design, or coefficient " +
	( rd.empty() ? std::string( "-" ) : str( row_coef( rd[ 1 ] , & x ) ) ) );
 const auto sname = milp_solver();
 if( ! sname.empty() ) {
  // charge 30 at 0 (3 modules of 10, paid 5 each), store 27, discharge
  // them at 1 (paid 5 each), investment 3
  const double v = milp_value( bu , sname );
  check( close( v , 3 - 5 * 30 - 5 * 27 ) , "battery modules: value " +
	 str( v ) + " instead of " + str( 3 - 5 * 30 - 5 * 27 ) + " (the "
	 "binary rows cap the power at one module)" );
  }
 bu->set_kappa( 2 );
 check( ( rd.size() == 2 ) && close( row_coef( rd[ 0 ] , & x ) , -8 ) &&
	close( lb0[ 0 ].get_rhs() , 24 ) ,
	"battery modules, set_kappa( 2 ): reserve rows not updated" );
 if( auto r = bu->get_max_outtake_binary_constraint( 1 ) )
  check( close( r->get_rhs() , 60 ) , "battery modules, set_kappa( 2 ): "
	 "binary row with right-hand side " + str( r->get_rhs() ) +
	 " instead of 2 * 3 * 10" );
 delete bu;

 bu = make( 2 , -3 );  // integer, at least 2 modules: accepted
 check( bu != nullptr , "battery, integer design in { 2 , 3 }: refused" );
 delete bu;

 bu = make( 0 , -1 );  // one module at most: the rows of one module
 if( bu ) {
  const auto & b0 = bu->get_primary_reserve_bounds();
  auto r = bu->get_max_outtake_binary_constraint( 0 );
  check( ( b0.size() == 2 ) && close( b0[ 0 ].get_rhs() , 4 ) && r &&
	 close( r->get_rhs() , 10 ) ,
	 "battery, one module: the rows are not those of one module" );
  delete bu;
  }
 }

/*--------------------------------------------------------------------------*/
/* A slack unit with the reactive power: the absolute value of the latter
 * costs 0.7 times ActivePowerCost, nothing when that is absent; a slack
 * unit in a UCBlock with inertia zones but no MaxInertia has no commitment
 * and its Objective is generated, while with MaxInertia the commitment
 * costs InertiaCost times MaxInertia. A negative MaxPower is a sink. */

static void test_slack_unit( void )
{
 auto make = []( const std::vector< std::pair< std::string ,
		  std::vector< double > > > & data , bool reactive ,
		 unsigned char reserve ) {
  auto g = new_group( "SU" , true );
  g.putAtt( "type" , "SlackUnitBlock" );
  auto T = g.addDim( "TimeHorizon" , 2 );
  for( const auto & d : data )
   put( g , d.first , T , d.second );
  auto sb = dynamic_cast< SlackUnitBlock * >( Block::new_Block( g ) );
  if( sb ) {
   sb->set_reactive_power( reactive );
   sb->set_reserve_vars( reserve );
   }
  return( sb );
  };

 auto all_coef = []( SlackUnitBlock * sb ) {
  auto lf = static_cast< const LinearFunction * >(
	    static_cast< FRealObjective * >( sb->get_objective() )
	                                                   ->get_function() );
  std::vector< double > c;
  for( const auto & vp : lf->get_v_var() )
   c.push_back( vp.second );
  return( c );
  };

 if( auto sb = make( { { "MaxPower" , { 100 , 100 } } } , true , 0 ) ) {
  try {
   generate_all( sb );
   const auto c = all_coef( sb );
   check( std::all_of( c.begin() , c.end() ,
		       []( double v ) { return( v == 0 ); } ) ,
	  "slack, reactive, no ActivePowerCost: a nonzero cost" );
   }
  catch( std::exception & e ) {
   check( false , std::string( "slack, reactive, no ActivePowerCost: " )
	  + e.what() );
   }
  delete sb;
  }

 if( auto sb = make( { { "MaxPower" , { 100 , 100 } } ,
		       { "ActivePowerCost" , { 10 , 20 } } } , true , 0 ) ) {
  generate_all( sb );
  const auto c = all_coef( sb );
  check( ( std::count_if( c.begin() , c.end() , []( double v ) {
	   return( close( v , 7 ) ); } ) == 1 ) &&
	 ( std::count_if( c.begin() , c.end() , []( double v ) {
	   return( close( v , 14 ) ); } ) == 1 ) ,
	 "slack, reactive: the absolute reactive power does not cost 0.7 "
	 "times ActivePowerCost" );
  delete sb;
  }

 if( auto sb = make( { { "MaxPower" , { -5 , 3 } } } , false , 4 ) ) {
  try {
   generate_all( sb );
   check( ( ! sb->has_commitment() ) &&
	  ( sb->get_commitment( 0 ) == nullptr ) ,
	  "slack, inertia zones, no MaxInertia: a commitment" );
   }
  catch( std::exception & e ) {
   check( false , std::string( "slack, inertia zones, no MaxInertia: " )
	  + e.what() );
   }
  check( ( sb->get_min_power( 0 ) == -5 ) && ( sb->get_max_power( 0 ) == 0 )
	 && ( sb->get_min_power( 1 ) == 0 ) && ( sb->get_max_power( 1 ) == 3 ) ,
	 "slack: the bounds of a negative and a positive MaxPower" );
  delete sb;
  }

 if( auto sb = make( { { "MaxPower" , { 5 , 5 } } ,
		       { "MaxInertia" , { 5 , 5 } } ,
		       { "InertiaCost" , { 2 , 2 } } } , false , 4 ) ) {
  generate_all( sb );
  auto u = sb->get_commitment( 0 );
  check( u && close( obj_coef( sb , & u[ 1 ] ) , 10 ) ,
	 "slack, MaxInertia: the commitment does not cost InertiaCost times "
	 "MaxInertia" );
  delete sb;
  }
 }

/*--------------------------------------------------------------------------*/
/*------------------------- THE LINKING ROWS OF UCBlock --------------------*/
/*--------------------------------------------------------------------------*/
/// a UCBlock on two nodes joined by one line of capacity 30, with three
/// IntermittentUnitBlock (costs 10, 30 and 50, the first at node 0, the
/// others at node 1, the third the only one with an inertia) and a
/// SlackUnitBlock at node 1; the caller adds the data of the linking rows

static netCDF::NcGroup write_UC_two_nodes( Index T ,
					   const std::vector< double > & demand ,
					   bool slack_inertia = false )
{
 auto g = new_group( "UC" , true );
 g.putAtt( "type" , "UCBlock" );
 auto TH = g.addDim( "TimeHorizon" , T );
 g.addDim( "NumberUnits" , 4 );
 auto N = g.addDim( "NumberNodes" , 2 );
 auto L = g.addDim( "NumberLines" , 1 );
 auto G = g.addDim( "NumberElectricalGenerators" , 4 );
 put_int( g , "StartLine" , L , { 0 } );
 put_int( g , "EndLine" , L , { 1 } );
 put( g , "MinPowerFlow" , L , { -30 } );
 put( g , "MaxPowerFlow" , L , { 30 } );
 put( g , "ActivePowerDemand" , { N , TH } , demand );
 put_int( g , "GeneratorNode" , G , { 0 , 1 , 1 , 1 } );

 const std::vector< double > cost = { 10 , 30 , 50 };
 for( Index i = 0 ; i < 3 ; ++i ) {
  auto u = g.addGroup( "UnitBlock_" + std::to_string( i ) );
  u.putAtt( "type" , "IntermittentUnitBlock" );
  put( u , "MaxPower" , 100.0 );
  put( u , "ActivePowerCost" , cost[ i ] );
  if( i < 2 )
   put( u , "Gamma" , 1.0 );
  else
   put( u , "InertiaPower" , 1.0 );
  }

 auto u = g.addGroup( "UnitBlock_3" );
 u.putAtt( "type" , "SlackUnitBlock" );
 put( u , "MaxPower" , 1000.0 );
 put( u , "ActivePowerCost" , 1000.0 );
 if( slack_inertia )
  put( u , "MaxInertia" , 0.0 );

 return( g );
 }

/*--------------------------------------------------------------------------*/
/// the units whose Variable are in the given row, in increasing order

static std::vector< Index > units_in( UCBlock * uc ,
				      const FRowConstraint & row )
{
 std::vector< Index > units;
 for( Index k = 0 ; k < row.get_num_active_var() ; ++k ) {
  const auto bk = row.get_active_var( k )->get_Block();
  for( Index i = 0 ; i < uc->get_number_units() ; ++i )
   if( bk == uc->get_unit_block( i ) )
    if( std::find( units.begin() , units.end() , i ) == units.end() )
     units.push_back( i );
  }
 std::sort( units.begin() , units.end() );
 return( units );
 }

/*--------------------------------------------------------------------------*/
/// the coefficient in the LinearFunction of the row of the first Variable
/// of the given unit, NaN if there is none

static double coef_of( UCBlock * uc , const FRowConstraint & row ,
		       Index unit )
{
 const auto lf = static_cast< const LinearFunction * >( row.get_function() );
 for( Index k = 0 ; k < row.get_num_active_var() ; ++k )
  if( row.get_active_var( k )->get_Block() == uc->get_unit_block( unit ) )
   return( lf->get_coefficient( k ) );
 return( std::numeric_limits< double >::quiet_NaN() );
 }

/*--------------------------------------------------------------------------*/
/* The zones of the reserve and inertia rows: with one primary zone the
 * marker of a node in no zone (a zone index >= NumberPrimaryZones) keeps
 * the units of that node out of the row, as it does with several zones
 * (C08); a zone vector absent with more than one zone of that kind is
 * refused (C33); and scaling, with the rows generated, a unit at a node in
 * no zone leaves the rows as they are, while scaling a unit in a zone
 * changes its coefficient (N2). */

static void test_zones_UCBlock( void )
{
 // one primary zone, node 1 in no zone: only unit 0 is in the row - - - -
 {
  auto g = write_UC_two_nodes( 1 , { 10 , 60 } );
  g.addDim( "NumberPrimaryZones" , 1 );
  auto N = g.getDim( "NumberNodes" );
  put_int( g , "PrimaryZones" , N , { 0 , 1 } );
  put( g , "PrimaryDemand" , { g.getDim( "NumberPrimaryZones" ) ,
                               g.getDim( "TimeHorizon" ) } , { 50 } );
  auto uc = dynamic_cast< UCBlock * >( Block::new_Block( g ) );
  check( uc , "zones, one zone with a marker: the UCBlock is not read" );
  if( uc ) {
   generate_all( uc );
   const auto units = units_in( uc ,
                                uc->get_const_primary_demand_constraints()
                                [ 0 ][ 0 ] );
   check( units == std::vector< Index >( { 0 } ) ,
	  "zones, one zone with a marker: the primary row has " +
	  std::to_string( units.size() ) + " units, not unit 0 alone" );
   delete uc;
   }
  }

 // one primary zone and no vector: every unit with a reserve is in the row
 {
  auto g = write_UC_two_nodes( 1 , { 10 , 60 } );
  g.addDim( "NumberPrimaryZones" , 1 );
  put( g , "PrimaryDemand" , { g.getDim( "NumberPrimaryZones" ) ,
                               g.getDim( "TimeHorizon" ) } , { 50 } );
  auto uc = dynamic_cast< UCBlock * >( Block::new_Block( g ) );
  check( uc , "zones, one zone, no vector: the UCBlock is not read" );
  if( uc ) {
   generate_all( uc );
   const auto units = units_in( uc ,
                                uc->get_const_primary_demand_constraints()
                                [ 0 ][ 0 ] );
   // the two units with a reserve, at nodes 0 and 1
   check( units == std::vector< Index >( { 0 , 1 } ) ,
	  "zones, one zone, no vector: the primary row has " +
	  std::to_string( units.size() ) + " units, not units 0 and 1" );
   delete uc;
   }
  }

 // two zones of each kind and no vector: refused- - - - - - - - - - - - - -
 for( const std::string kind : { "Primary" , "Secondary" , "Inertia" } ) {
  auto g = write_UC_two_nodes( 1 , { 10 , 60 } , true );
  auto Z = g.addDim( "Number" + kind + "Zones" , 2 );
  put( g , kind + "Demand" , { Z , g.getDim( "TimeHorizon" ) } , { 5 , 5 } );
  Block * b = nullptr;
  bool refused = false;
  try {
   b = Block::new_Block( g );
   refused = ( b == nullptr );
   }
  catch( std::exception & ) {
   refused = true;
   }
  check( refused , "zones: two " + kind + " zones and no " + kind +
	 "Zones are not refused" );
  delete b;
  }

 // two primary zones, node 1 in no zone; scale the units- - - - - - - - - -
 {
  auto g = write_UC_two_nodes( 1 , { 10 , 60 } );
  auto Z = g.addDim( "NumberPrimaryZones" , 2 );
  put_int( g , "PrimaryZones" , g.getDim( "NumberNodes" ) , { 0 , 2 } );
  put( g , "PrimaryDemand" , { Z , g.getDim( "TimeHorizon" ) } , { 5 , 0 } );
  auto uc = dynamic_cast< UCBlock * >( Block::new_Block( g ) );
  check( uc , "zones, node in no zone: the UCBlock is not read" );
  if( uc ) {
   generate_all( uc );
   auto fs = new FakeSolver();
   uc->register_Solver( fs );
   const auto & rows = uc->get_const_primary_demand_constraints();
   check( ( units_in( uc , rows[ 0 ][ 0 ] ) ==
	    std::vector< Index >( { 0 } ) ) &&
	  units_in( uc , rows[ 0 ][ 1 ] ).empty() ,
	  "zones, node in no zone: the rows do not hold unit 0 alone" );

   uc->get_unit_block( 1 )->scale( 0.5 , eModBlck , eModBlck );
   check( ( units_in( uc , rows[ 0 ][ 0 ] ) ==
	    std::vector< Index >( { 0 } ) ) &&
	  units_in( uc , rows[ 0 ][ 1 ] ).empty() &&
	  close( coef_of( uc , rows[ 0 ][ 0 ] , 0 ) , 1 ) ,
	  "zones, node in no zone: scaling a unit in no zone changes the "
	  "rows" );

   uc->get_unit_block( 0 )->scale( 3 , eModBlck , eModBlck );
   check( close( coef_of( uc , rows[ 0 ][ 0 ] , 0 ) , 3 ) ,
	  "zones, node in no zone: scaling unit 0 puts " +
	  str( coef_of( uc , rows[ 0 ][ 0 ] , 0 ) ) +
	  " in its primary row, not 3" );

   uc->unregister_Solver( fs , true );
   delete uc;
   }
  }
 }

/*--------------------------------------------------------------------------*/
/* An inertia zone with a SlackUnitBlock that gives no inertia, i.e., has no
 * "MaxInertia" (N1): the abstract representation is generated, and the
 * slack is not in the inertia rows. */

static void test_slack_no_inertia_UCBlock( void )
{
 auto g = write_UC_two_nodes( 2 , { 10 , 10 , 60 , 20 } );
 auto Z = g.addDim( "NumberInertiaZones" , 1 );
 put( g , "InertiaDemand" , { Z , g.getDim( "TimeHorizon" ) } , { 15 , 0 } );
 auto uc = dynamic_cast< UCBlock * >( Block::new_Block( g ) );
 check( uc , "slack without inertia: the UCBlock is not read" );
 if( ! uc )
  return;
 try {
  generate_all( uc );
  const auto units = units_in( uc ,
			       uc->get_const_inertia_demand_constraints()
			       [ 0 ][ 0 ] );
  check( std::find( units.begin() , units.end() , 3 ) == units.end() ,
	 "slack without inertia: the slack is in the inertia row" );
  }
 catch( std::exception & e ) {
  check( false , std::string( "slack without inertia: the abstract "
			      "representation throws " ) + e.what() );
  }
 delete uc;
 }

/*--------------------------------------------------------------------------*/
/* The constant terms of the NetworkBlock (C09): the "ConstantTerm" of a
 * "NetworkBlock_n" group prevails over NetworkConstantTerms[ n ], which
 * applies to the NetworkBlock without it; with a single node, where the
 * NetworkBlock are deleted, the constants go in the Objective of the
 * UCBlock and survive a round trip. */

static void test_network_constants_UCBlock( void )
{
 auto network_group = []( netCDF::NcGroup & g , Index n , double c ) {
  auto nb = g.addGroup( "NetworkBlock_" + std::to_string( n ) );
  nb.putAtt( "type" , "DCNetworkBlock" );
  put( nb , "ConstantTerm" , c );
  };

 auto read = []( netCDF::NcGroup & g , const std::string & what ) {
  UCBlock * uc = nullptr;
  try {
   uc = dynamic_cast< UCBlock * >( Block::new_Block( g ) );
   }
  catch( std::exception & e ) {
   check( false , what + ": deserialize() throws " + e.what() );
   return( uc );
   }
  check( uc , what + ": the UCBlock is not read" );
  return( uc );
  };

 auto ucb_constant = []( UCBlock * uc ) {
  generate_all( uc );
  return( static_cast< FRealObjective * >( uc->get_objective() )
	  ->get_constant_term() );
  };

 // two nodes, the group with its constant, no NetworkConstantTerms - - - - -
 {
  auto g = write_UC_two_nodes( 2 , { 10 , 10 , 60 , 20 } );
  g.addDim( "NumberNetworks" , 2 );
  network_group( g , 0 , 1000 );
  if( auto uc = read( g , "network constants, own" ) ) {
   check( uc->get_network_block( 0 )->get_const_term() == 1000 ,
	  "network constants, own: the ConstantTerm of the group is lost" );
   check( uc->get_network_block( 1 )->get_const_term() == 0 ,
	  "network constants, own: a NetworkBlock without one has " +
	  str( uc->get_network_block( 1 )->get_const_term() ) );
   check( ucb_constant( uc ) == 0 ,
	  "network constants, own: the UCBlock has a constant" );
   delete uc;
   }
  }

 // two nodes, the group with its constant and NetworkConstantTerms- - - - -
 {
  auto g = write_UC_two_nodes( 2 , { 10 , 10 , 60 , 20 } );
  auto NN = g.addDim( "NumberNetworks" , 2 );
  put( g , "NetworkConstantTerms" , NN , { 7 , 8 } );
  network_group( g , 0 , 1000 );
  if( auto uc = read( g , "network constants, both" ) ) {
   check( uc->get_network_block( 0 )->get_const_term() == 1000 ,
	  "network constants, both: NetworkConstantTerms prevails over the "
	  "ConstantTerm of the group" );
   check( uc->get_network_block( 1 )->get_const_term() == 8 ,
	  "network constants, both: NetworkConstantTerms is not applied to "
	  "the NetworkBlock without a group" );
   delete uc;
   }
  }

 // one node, the group with its constant: in the UCBlock Objective- - - - -
 {
  auto g = new_group( "UC" , true );
  g.putAtt( "type" , "UCBlock" );
  auto T = g.addDim( "TimeHorizon" , 2 );
  g.addDim( "NumberUnits" , 1 );
  g.addDim( "NumberNetworks" , 2 );
  auto one = g.addDim( "NumberNodes" , 1 );
  put( g , "ActivePowerDemand" , { one , T } , { 3 , 4 } );
  auto u = g.addGroup( "UnitBlock_0" );
  u.putAtt( "type" , "SlackUnitBlock" );
  put( u , "MaxPower" , 100.0 );
  put( u , "ActivePowerCost" , 1.0 );
  network_group( g , 1 , 1000 );
  if( auto uc = read( g , "network constants, one node" ) ) {
   check( uc->get_network_blocks().empty() ,
	  "network constants, one node: NetworkBlock left" );
   check( ucb_constant( uc ) == 1000 ,
	  "network constants, one node: the UCBlock Objective has the "
	  "constant " + str( ucb_constant( uc ) ) + ", not 1000" );
   delete uc;
   }
  if( auto b = dynamic_cast< UCBlock * >( round_trip( g ,
					   "network constants, one node" ) ) ) {
   check( ucb_constant( b ) == 1000 ,
	  "network constants, one node: the constant read back is " +
	  str( ucb_constant( b ) ) );
   delete b;
   }
  }

 // one node, NetworkConstantTerms only- - - - - - - - - - - - - - - - - - -
 {
  auto g = new_group( "UC" , true );
  g.putAtt( "type" , "UCBlock" );
  auto T = g.addDim( "TimeHorizon" , 2 );
  g.addDim( "NumberUnits" , 1 );
  auto NN = g.addDim( "NumberNetworks" , 2 );
  auto one = g.addDim( "NumberNodes" , 1 );
  put( g , "ActivePowerDemand" , { one , T } , { 3 , 4 } );
  put( g , "NetworkConstantTerms" , NN , { 5 , 6 } );
  auto u = g.addGroup( "UnitBlock_0" );
  u.putAtt( "type" , "SlackUnitBlock" );
  put( u , "MaxPower" , 100.0 );
  put( u , "ActivePowerCost" , 1.0 );
  if( auto uc = read( g , "network constants, one node, NCT" ) ) {
   check( ucb_constant( uc ) == 11 ,
	  "network constants, one node, NCT: the UCBlock Objective has the "
	  "constant " + str( ucb_constant( uc ) ) + ", not 11" );
   delete uc;
   }
  }
 }

/*--------------------------------------------------------------------------*/
/* The sign of the dual values of the linking rows (C48), on instances whose
 * duals are known: the dual of a node injection row is minus the marginal
 * cost of the demand of the node, those of the >= rows of the reserves and
 * of the inertia are minus the marginal cost of the requirement, that of a
 * pollutant row is the marginal saving of a unit more of budget when the
 * budget binds, minus the marginal cost of a unit more of the lower bound
 * when the latter binds. The ActivePowerDuals of the UCBlockSolution are
 * indexed [ TimeHorizon , NumberNodes ]. Needs a :MILPSolver in the
 * factory, skipped otherwise. */

static void test_dual_signs_UCBlock( void )
{
 std::string sname;
 for( const auto & name : { "CPXMILPSolver" , "GRBMILPSolver" ,
                            "HiGHSMILPSolver" , "SCIPMILPSolver" } )
  if( Solver::has_Solver( name ) ) {
   sname = name;
   break;
   }

 if( sname.empty() ) {
  std::cout << "no :MILPSolver in this build, the signs of the duals are "
	    << "not checked" << std::endl;
  return;
  }

 // solves the UCBlock read from g, returns it with its duals, or nullptr
 auto solve = [ & ]( netCDF::NcGroup & g , double opt ,
		     const std::string & what ) -> UCBlock * {
  auto uc = dynamic_cast< UCBlock * >( Block::new_Block( g ) );
  if( ! uc ) {
   check( false , what + ": the UCBlock is not read" );
   return( nullptr );
   }
  BlockSolverConfig bsc( 1 );
  bsc.add_ComputeConfig( std::string( sname ) , nullptr );
  bsc.apply( uc );
  auto slv = static_cast< CDASolver * >(
                                     uc->get_registered_solvers().front() );
  const auto par = slv->int_par_str2idx( "intLogVerb" );
  if( par < Inf< Solver::idx_type >() )
   slv->set_par( par , 0 );
  const auto status = slv->compute( false );
  bool ok = ( status == Solver::kOK ) && slv->has_dual_solution();
  if( ok ) {
   slv->get_var_solution();
   slv->get_dual_solution();
   ok = close( slv->get_var_value() , opt );
   }
  check( ok , what + ": " + sname + " gives status " +
	 std::to_string( status ) + " and value " +
	 str( slv->get_var_value() ) + ", not " + str( opt ) );
  bsc.clear();
  bsc.apply( uc );
  if( ok )
   return( uc );
  delete uc;
  return( nullptr );
  };

 auto check_dual = []( const FRowConstraint & row , double expected ,
		       const std::string & what ) {
  check( close( row.get_dual() , expected , 1e-5 ) ,
	 what + ": the dual is " + str( row.get_dual() ) + ", not " +
	 str( expected ) );
  };

 // node injection, inertia and pollutant budget- - - - - - - - - - - - - -
 // t = 0: the line binds, the inertia 15 is met by unit 2 alone; t = 1:
 // the emission of unit 0, the only one emitting, binds at 25
 {
  auto g = write_UC_two_nodes( 2 , { 10 , 10 , 60 , 20 } , true );
  auto TH = g.getDim( "TimeHorizon" );
  auto Z = g.addDim( "NumberInertiaZones" , 1 );
  put( g , "InertiaDemand" , { Z , TH } , { 15 , 0 } );
  g.addDim( "NumberPollutants" , 1 );
  auto TZ = g.addDim( "TotalNumberPollutantZones" , 1 );
  put( g , "PollutantBudget" , TZ , { 25 } );
  put( g , "PollutantRho" , { TH , g.getDim( "NumberPollutants" ) ,
			      g.getDim( "NumberElectricalGenerators" ) } ,
       { 0 , 0 , 0 , 0 , 1 , 0 , 0 , 0 } );
  if( auto uc = solve( g , 2000 , "duals, base" ) ) {
   const auto & nic = uc->get_const_node_injection_constraints();
   check_dual( nic[ 0 ][ 0 ] , -10 , "duals, base, node 0 at t = 0" );
   check_dual( nic[ 0 ][ 1 ] , -30 , "duals, base, node 1 at t = 0" );
   check_dual( nic[ 1 ][ 0 ] , -30 , "duals, base, node 0 at t = 1" );
   check_dual( nic[ 1 ][ 1 ] , -30 , "duals, base, node 1 at t = 1" );
   check_dual( uc->get_const_inertia_demand_constraints()[ 0 ][ 0 ] , -20 ,
	       "duals, base, inertia at t = 0" );
   check_dual( uc->get_const_pollutant_constraints()[ 0 ][ 0 ] , 20 ,
	       "duals, base, pollutant budget" );

   SimpleConfiguration< int > wsol( 8 );
   auto sol = uc->get_Solution( & wsol , false );
   auto sg = new_group( "SOL" );
   sol->serialize( sg );
   auto apd = sg.getVar( "ActivePowerDuals" );
   check( ( ! apd.isNull() ) && ( apd.getDimCount() == 2 ) &&
	  ( apd.getDim( 0 ).getName() == "TimeHorizon" ) &&
	  ( apd.getDim( 1 ).getName() == "NumberNodes" ) ,
	  "duals, base: ActivePowerDuals is not indexed [ TimeHorizon , "
	  "NumberNodes ]" );
   if( ( ! apd.isNull() ) && ( apd.getDimCount() == 2 ) ) {
    std::vector< double > v( 4 );
    apd.getVar( v.data() );
    check( close( v[ 1 ] , -30 , 1e-5 ) && close( v[ 2 ] , -30 , 1e-5 ) ,
	   "duals, base: ActivePowerDuals[ 0 , 1 ] is " + str( v[ 1 ] ) +
	   ", not -30" );
    }
   delete sol;
   delete uc;
   }
  }

 // primary and secondary reserve - - - - - - - - - - - - - - - - - - - - -
 {
  auto g = write_UC_two_nodes( 1 , { 80 , 20 } );
  auto TH = g.getDim( "TimeHorizon" );
  auto Zp = g.addDim( "NumberPrimaryZones" , 1 );
  auto Zs = g.addDim( "NumberSecondaryZones" , 1 );
  put( g , "PrimaryDemand" , { Zp , TH } , { 4 } );
  put( g , "SecondaryDemand" , { Zs , TH } , { 6 } );
  if( auto uc = solve( g , 1100 , "duals, reserves" ) ) {
   const auto & nic = uc->get_const_node_injection_constraints();
   check_dual( nic[ 0 ][ 0 ] , -20 , "duals, reserves, node 0" );
   check_dual( nic[ 0 ][ 1 ] , -20 , "duals, reserves, node 1" );
   check_dual( uc->get_const_primary_demand_constraints()[ 0 ][ 0 ] , -10 ,
	       "duals, reserves, primary" );
   check_dual( uc->get_const_secondary_demand_constraints()[ 0 ][ 0 ] , -10 ,
	       "duals, reserves, secondary" );
   delete uc;
   }
  }

 // pollutant row whose lower bound binds- - - - - - - - - - - - - - - - - -
 {
  auto g = write_UC_two_nodes( 2 , { 10 , 10 , 60 , 20 } );
  g.addDim( "NumberPollutants" , 1 );
  auto TZ = g.addDim( "TotalNumberPollutantZones" , 1 );
  auto PT = g.addDim( "PollutantRhoTime" , 1 );
  put( g , "PollutantBudget" , TZ , { 1000 } );
  put( g , "PollutantMinBudget" , TZ , { 50 } );
  put( g , "PollutantRho" , { PT , g.getDim( "NumberPollutants" ) ,
			      g.getDim( "NumberElectricalGenerators" ) } ,
       { 0 , 0 , 1 , 0 } );
  if( auto uc = solve( g , 3000 , "duals, lower budget" ) ) {
   check_dual( uc->get_const_pollutant_constraints()[ 0 ][ 0 ] , -40 ,
	       "duals, lower budget" );
   delete uc;
   }
  }
 }

/*--------------------------------------------------------------------------*/
/*---------------------------------- MAIN ----------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------- THERMAL UNIT: DATA AND 3BIN EDGE CASES ---------------*/
/*--------------------------------------------------------------------------*/
/* A ThermalUnitBlock written from per-instant vectors, so that the data
 * that TUData keeps constant (or always writes, as the ramps) can change
 * with the instant or be left out; names and vectors are written as they
 * are, over "NumberIntervals" = T. */

static ThermalUnitBlock * new_TU_vec(
 Index T , const std::vector< std::pair< std::string ,
                                         std::vector< double > > > & vecs ,
 int initUD , double initP , unsigned int minUp , unsigned int minDown )
{
 auto g = new_group( "TUv" , true );
 g.putAtt( "type" , "ThermalUnitBlock" );
 g.addDim( "TimeHorizon" , T );
 auto NI = g.addDim( "NumberIntervals" , T );
 for( const auto & nv : vecs )
  put( g , nv.first , NI , nv.second );
 if( initUD > 0 )
  put( g , "InitialPower" , initP );
 put_int( g , "InitUpDownTime" , initUD );
 put_uint( g , "MinUpTime" , minUp );
 put_uint( g , "MinDownTime" , minDown );
 auto tub = dynamic_cast< ThermalUnitBlock * >( Block::new_Block( g ) );
 if( ! tub )
  throw( std::logic_error( "new_TU_vec: no ThermalUnitBlock built" ) );
 return( tub );
 }

/*--------------------------------------------------------------------------*/
/// the unit with the abstract representation of the formulation wf

static void generate_all_wf( Block * b , int wf )
{
 SimpleConfiguration< int > f( wf );
 b->generate_abstract_variables( & f );
 b->generate_abstract_constraints();
 b->generate_objective();
 }

/*--------------------------------------------------------------------------*/
/* The maximum power rows of the 3bin formulation: a horizon of one instant
 * with the minimum up time 1 or 2 and the commitment free at 0, and a unit
 * held off for the whole horizon (init_t == T) with the minimum up time 1,
 * the cases in which the rows were counted wrong or read past the end; the
 * DP Solvers write their schedule in the Variable, which has to satisfy
 * the rows. Then a start-up at 0 with the minimum up time 1 and no ramp,
 * whose power has to be capped by the start-up limit as in the T
 * formulation and in the DP Solvers. */

static void test_thermal_3bin_max_power( void )
{
 // one instant, on since long at 5, can shut down at 0
 for( unsigned int up : { 1u , 2u } ) {
  TUData d;
  d.minP = 1;
  d.maxP = d.su = d.sd = d.ru = d.rd = 10;
  d.initUD = 5;
  d.initP = 5;
  d.minUp = up;
  d.lin = { -1 };
  const auto what = "3bin, T = 1, MinUpTime " + std::to_string( up );
  auto tub = new_TU( d );
  try {
   generate_all_wf( tub , ThermalUnitBlock::tbinForm );
   check_all_DP( tub , brute_force( d , { -1 } ) , what );
   }
  catch( std::exception & e ) {
   check( false , what + ": throws " + e.what() );
   }
  delete tub;
  }

 // off for 1 instant, MinDownTime 5: init_t = min( 4 , T ) = T = 2
 {
  TUData d;
  d.T = 2;
  d.minP = 1;
  d.maxP = d.su = d.sd = d.ru = d.rd = 10;
  d.initUD = -1;
  d.minDown = 5;
  d.lin = { -1 , -1 };
  const std::string what = "3bin, init_t == T";
  auto tub = new_TU( d );
  try {
   generate_all_wf( tub , ThermalUnitBlock::tbinForm );
   check_all_DP( tub , brute_force( d , { 0 , 0 } ) , what );
   }
  catch( std::exception & e ) {
   check( false , what + ": throws " + e.what() );
   }
  delete tub;
  }

 // off since long, MinUpTime 1, StartUpLimit 5 < MaxPower 10, no ramp: a
 // start-up at 0 at full power is not feasible, in 3bin as in T
 for( int wf : { int( ThermalUnitBlock::tbinForm ) ,
                 int( ThermalUnitBlock::TForm ) } ) {
  const auto what = "start-up at 0 above StartUpLimit, wf " +
                    std::to_string( wf );
  auto tub = new_TU_vec( 2 , { { "MinPower" , { 1 , 1 } } ,
                               { "MaxPower" , { 10 , 10 } } ,
                               { "StartUpLimit" , { 5 , 5 } } ,
                               { "ShutDownLimit" , { 10 , 10 } } ,
                               { "LinearTerm" , { -1 , -1 } } } ,
                         -5 , 0 , 1 , 1 );
  generate_all_wf( tub , wf );
  auto u = tub->get_commitment( 0 );
  auto p = tub->get_active_power( 0 );
  auto v = tub->get_start_up();
  auto w = tub->get_shut_down();
  // on at 0 and 1, started at 0 at 10
  u[ 0 ].set_value( 1 );  u[ 1 ].set_value( 1 );
  v[ 0 ].set_value( 1 );  v[ 1 ].set_value( 0 );
  w[ 0 ].set_value( 0 );  w[ 1 ].set_value( 0 );
  p[ 0 ].set_value( 10 );  p[ 1 ].set_value( 10 );
  SimpleConfiguration< double > tol( 1e-7 );
  check( ! tub->is_feasible( true , & tol ) ,
         what + ": p_0 = 10 > StartUpLimit[ 0 ] = 5 is feasible" );
  // at the start-up limit it is
  p[ 0 ].set_value( 5 );
  check( tub->is_feasible( true , & tol ) ,
         what + ": p_0 = 5 = StartUpLimit[ 0 ] is not feasible" );
  // the optimum: 5 at 0 and 10 at 1, i.e. -15, in the DP Solvers as well
  TUData d;
  d.T = 2;
  d.minP = 1;
  d.maxP = d.sd = d.ru = d.rd = 10;
  d.su = 5;
  d.initUD = -5;
  d.lin = { -1 , -1 };
  check( close( brute_force( d , { -1 , -1 } ) , -15 ) ,
         what + ": the brute force gives " +
         str( brute_force( d , { -1 , -1 } ) ) );
  check_all_DP( tub , -15 , what );
  delete tub;
  }
 }

/*--------------------------------------------------------------------------*/
/* The checks of the data of a ThermalUnitBlock: an InitialPower below the
 * minimum power of a unit on before the horizon is raised to MinPower[ 0 ]
 * also when the unit is unavailable at 0 (where the operational minimum
 * power is 0), and one above MaxPower[ 0 ] is only warned about; a negative
 * FixedConsumption is rejected, a zero one is not; the check of the
 * initial power against the ramp-down of a unit on before the horizon
 * throws with the name of the method. */

static void test_thermal_data_checks( void )
{
 // InitialPower 0.5 < MinPower 1, available at 0 or not
 for( double av0 : { 1.0 , 0.0 } ) {
  const auto what = "InitialPower below MinPower, Availability[ 0 ] " +
                    str( av0 );
  try {
   auto tub = new_TU_vec( 2 , { { "MinPower" , { 1 , 1 } } ,
                                { "MaxPower" , { 10 , 10 } } ,
                                { "Availability" , { av0 , 1 } } ,
                                { "StartUpLimit" , { av0 , 10 } } ,
                                { "ShutDownLimit" , { av0 , 10 } } ,
                                { "DeltaRampDown" , { 10 , 10 } } } ,
                          3 , 0.5 , 1 , 1 );
   check( close( tub->get_initial_power() , 1 ) ,
          what + ": InitialPower is " + str( tub->get_initial_power() ) +
          " instead of 1" );
   delete tub;
   }
  catch( std::exception & e ) {
   check( false , what + ": throws " + e.what() );
   }
  }

 // InitialPower 12 > MaxPower 10: kept, with a warning
 try {
  auto tub = new_TU_vec( 1 , { { "MinPower" , { 1 } } ,
                               { "MaxPower" , { 10 } } } ,
                         3 , 12 , 1 , 1 );
  check( close( tub->get_initial_power() , 12 ) ,
         "InitialPower above MaxPower is changed to " +
         str( tub->get_initial_power() ) );
  delete tub;
  }
 catch( std::exception & e ) {
  check( false , std::string( "InitialPower above MaxPower: throws " ) +
         e.what() );
  }

 // FixedConsumption: -1 is rejected, 0 is not
 for( double fc : { -1.0 , 0.0 } ) {
  bool thrown = false;
  try {
   auto tub = new_TU_vec( 2 , { { "MinPower" , { 1 , 1 } } ,
                                { "MaxPower" , { 10 , 10 } } ,
                                { "FixedConsumption" , { 0 , fc } } } ,
                          -3 , 0 , 1 , 1 );
   delete tub;
   }
  catch( std::logic_error & ) {
   thrown = true;
   }
  check( thrown == ( fc < 0 ) ,
         "FixedConsumption " + str( fc ) + ( thrown ? " is" : " is not" ) +
         " rejected" );
  }

 // on before at 10 = MaxPower, half available at 0 (operational maximum 5),
 // DeltaRampDown 2: 10 - 2 > 5, which generate_abstract_constraints() rejects
 {
  auto tub = new_TU_vec( 2 , { { "MinPower" , { 1 , 1 } } ,
                               { "MaxPower" , { 10 , 10 } } ,
                               { "Availability" , { 0.5 , 1 } } ,
                               { "StartUpLimit" , { 5 , 10 } } ,
                               { "ShutDownLimit" , { 5 , 10 } } ,
                               { "DeltaRampDown" , { 2 , 2 } } } ,
                         3 , 10 , 1 , 1 );
  std::string msg;
  try {
   generate_all_wf( tub , ThermalUnitBlock::TForm );
   }
  catch( std::logic_error & e ) {
   msg = e.what();
   }
  const std::string who = "ThermalUnitBlock::generate_abstract_constraints";
  check( msg.compare( 0 , who.size() , who ) == 0 ,
         "infeasible initial ramp-down: the message is \"" + msg + "\"" );
  delete tub;
  }
 }

/*--------------------------------------------------------------------------*/
/* The DP formulation caps the power of a run of one instant by the smaller
 * of the start-up and of the shut-down limit, as the other formulations and
 * the DP Solvers do: a unit off since long, on at 0 only, StartUpLimit 10
 * and ShutDownLimit 3, cannot produce more than 3 at 0. The schedule is
 * written in every Variable of the formulation, the run being the interval
 * ( h , k ) = ( 1 , 1 ) of the state-space graph. */

static void test_thermal_DP_one_instant_run( void )
{
 const std::string what = "DP formulation, run of one instant";
 auto tub = new_TU_vec( 2 , { { "MinPower" , { 1 , 1 } } ,
                              { "MaxPower" , { 10 , 10 } } ,
                              { "StartUpLimit" , { 10 , 10 } } ,
                              { "ShutDownLimit" , { 3 , 3 } } ,
                              { "LinearTerm" , { -1 , 20 } } } ,
                        -5 , 0 , 1 , 1 );
 generate_all_wf( tub , ThermalUnitBlock::DPForm );

 auto u = tub->get_commitment( 0 );
 auto p = tub->get_active_power( 0 );
 auto v = tub->get_start_up();
 auto w = tub->get_shut_down();
 auto yp = tub->get_static_variable_v< ColVariable >( "y_plus_thermal" );
 auto ym = tub->get_static_variable_v< ColVariable >( "y_minus_thermal" );
 auto phk = tub->get_static_variable_v< ColVariable >( "p_h_k_thermal" );
 const auto & Y = tub->get_Y_plus();
 if( ( ! yp ) || ( ! ym ) || ( ! phk ) || ( ym->size() != 5 ) ) {
  check( false , what + ": unexpected Variable of the DP formulation" );
  delete tub;
  return;
  }

 // on at 0 only: started at 0, shut down at 1
 u[ 0 ].set_value( 1 );  u[ 1 ].set_value( 0 );
 v[ 0 ].set_value( 1 );  v[ 1 ].set_value( 0 );
 w[ 0 ].set_value( 0 );  w[ 1 ].set_value( 1 );
 p[ 1 ].set_value( 0 );
 for( auto & y : *yp ) y.set_value( 0 );
 for( auto & y : *ym ) y.set_value( 0 );
 for( auto & q : *phk ) q.set_value( 0 );
 // the arcs: y^-( 0 , 1 ) (off until the start at 0), y^+( 1 , 1 ), and
 // y^-( 1 , 3 ) (off from the shut-down at 1 to the end), in the order in
 // which generate_abstract_variables() makes them, ( 0 , 1 ) , ( 0 , 2 ) ,
 // ( 0 , 3 ) , ( 1 , 3 ) , ( 2 , 3 )
 ( *ym )[ 0 ].set_value( 1 );
 ( *ym )[ 3 ].set_value( 1 );
 // p^{hk}_t are made interval by interval, instant by instant
 ColVariable * p11 = nullptr;
 for( Index i = 0 , j = 0 ; i < Y.size() ; ++i )
  for( Index t = 0 ; t < 2 ; ++t )
   if( ( Y[ i ].first <= t + 1 ) && ( t + 1 <= Y[ i ].second ) ) {
    if( ( Y[ i ].first == 1 ) && ( Y[ i ].second == 1 ) ) {
     ( *yp )[ i ].set_value( 1 );
     p11 = & ( *phk )[ j ];
     }
    ++j;
    }
 if( ! p11 ) {
  check( false , what + ": no run ( 1 , 1 )" );
  delete tub;
  return;
  }

 SimpleConfiguration< double > tol( 1e-7 );
 for( double p0 : { 10.0 , 3.0 } ) {
  p[ 0 ].set_value( p0 );
  p11->set_value( p0 );
  check( tub->is_feasible( true , & tol ) == ( p0 <= 3 ) ,
         what + ": p_0 = " + str( p0 ) + ( p0 <= 3 ? " is not" : " is" ) +
         " feasible" );
  }

 delete tub;

 // the optimum, -3 (staying on at 1 costs at least 20), in the DP Solvers
 // as well, on a twin unit with no abstract representation, since the
 // Solution of a DP Solver does not set the Variable of the DP formulation
 tub = new_TU_vec( 2 , { { "MinPower" , { 1 , 1 } } ,
                         { "MaxPower" , { 10 , 10 } } ,
                         { "StartUpLimit" , { 10 , 10 } } ,
                         { "ShutDownLimit" , { 3 , 3 } } ,
                         { "LinearTerm" , { -1 , 20 } } } ,
                   -5 , 0 , 1 , 1 );
 check_all_DP( tub , -3 , what );
 delete tub;
 }

/*--------------------------------------------------------------------------*/
/* The initial power changed after the abstract representation is
 * generated: in every formulation the rows that contain it follow (in the
 * SUSD one some of them appear or disappear), and the schedule of the DP
 * Solvers (which read the new value) satisfies them and costs what a fresh
 * load of the unit with the new value costs, both when the initial power
 * goes down (ramp-down binding, the unit forced on) and up (ramp-up
 * binding); the initial up/down time cannot change after the generation,
 * and the change is refused with std::logic_error. Then the same changes,
 * and that of the initial up/down time, before the generation: the SUSD
 * formulation, whose number of ramp rows from the initial power depends on
 * both, has then as many rows as a fresh load. */

static void test_thermal_initial_power_change( void )
{
 for( double lin : { -1.0 , 1.0 } )
  for( double ip : { 20.0 , 90.0 } )
   for( int wf = 0 ; wf < 7 ; ++wf ) {
    TUData d;
    d.T = 4;
    d.minP = 10;
    d.maxP = d.su = d.sd = 100;
    d.ru = d.rd = 15;
    d.initUD = 3;
    d.initP = 50;
    d.minUp = lin > 0 ? 10 : 1;  // forced on with a positive price
    d.lin = std::vector< double >( d.T , lin );
    const auto what = "InitialPower 50 -> " + str( ip ) + ", price " +
                      str( lin ) + ", wf " + std::to_string( wf );
    auto tub = new_TU( d );
    generate_all_wf( tub , wf );
    std::vector< double > val = { ip };
    try {
     tub->set_initial_power( val.begin() , Range( 0 , 1 ) );
     }
    catch( std::exception & e ) {
     check( false , what + ": throws " + e.what() );
     delete tub;
     continue;
     }
    d.initP = ip;
    check_all_DP( tub , brute_force( d , std::vector< int >( d.T , -1 ) ) ,
                  what );
    delete tub;
    }

 // the initial up/down time cannot change after the generation
 {
  TUData d;
  d.T = 4;
  d.minP = 10;
  d.maxP = d.su = d.sd = 100;
  d.ru = d.rd = 15;
  d.initUD = 3;
  d.initP = 50;
  auto tub = new_TU( d );
  generate_all_wf( tub , ThermalUnitBlock::ptForm );
  std::vector< int > ud = { -2 };
  bool thrown = false;
  try {
   tub->set_init_updown_time( ud.begin() , Range( 0 , 1 ) );
   }
  catch( std::logic_error & ) {
   thrown = true;
   }
  check( thrown && ( tub->get_init_up_down_time() == 3 ) ,
         "InitUpDownTime after the generation: " +
         std::string( thrown ? "refused" : "accepted" ) + ", now " +
         std::to_string( tub->get_init_up_down_time() ) );
  delete tub;
  }

 // before the generation, SUSD: as many ramp rows as a fresh load (the
 // rows of the ramps of several steps are dynamic)
 auto rows = []( ThermalUnitBlock * b , const std::string & n ) -> Index {
  auto v = b->get_static_constraint_v< FRowConstraint >( n );
  auto l = b->get_dynamic_constraint< FRowConstraint >( n );
  return( ( v ? v->size() : 0 ) + ( l ? l->size() : 0 ) );
  };
 for( int c = 0 ; c < 3 ; ++c ) {
  TUData d;
  d.T = 6;
  d.minP = 10;
  d.maxP = d.su = d.sd = 100;
  d.ru = d.rd = 15;
  d.initUD = c == 0 ? -1 : 3;
  d.initP = c == 1 ? 20 : 90;
  auto tub = new_TU( d );
  std::string what;
  if( c == 0 ) {  // off -> on before the horizon
   std::vector< int > ud = { 3 };
   tub->set_init_updown_time( ud.begin() , Range( 0 , 1 ) );
   std::vector< double > val = { 50 };
   tub->set_initial_power( val.begin() , Range( 0 , 1 ) );
   d.initUD = 3;
   d.initP = 50;
   what = "SUSD, InitUpDownTime -1 -> 3 before the generation";
   }
  else {
   std::vector< double > val = { c == 1 ? 90.0 : 20.0 };
   tub->set_initial_power( val.begin() , Range( 0 , 1 ) );
   d.initP = val[ 0 ];
   what = "SUSD, InitialPower " + str( c == 1 ? 20 : 90 ) + " -> " +
          str( val[ 0 ] ) + " before the generation";
   }
  generate_all_wf( tub , ThermalUnitBlock::SUSDForm );
  auto fresh = new_TU( d );
  generate_all_wf( fresh , ThermalUnitBlock::SUSDForm );
  for( const std::string n : { "RampUp_Const_Thermal" ,
                               "RampDown_Const_Thermal" ,
                               "RampUp_SUSD_Const_Thermal" ,
                               "RampDown_SUSD_Const_Thermal" } )
   check( rows( tub , n ) == rows( fresh , n ) ,
          what + ": " + std::to_string( rows( tub , n ) ) + " rows of " + n +
          " instead of " + std::to_string( rows( fresh , n ) ) );
  delete fresh;
  delete tub;
  }
 }

/*--------------------------------------------------------------------------*/
/*---------------- NUCLEAR: BANDS AT THE END OF THE HORIZON ----------------*/
/*--------------------------------------------------------------------------*/
/* A NuclearUnitBlock whose output is split into the bands [ 10 , 20 ],
 * [ 20 , 40 ] and [ 40 , 50 ], whose modulations may last 2 instants and
 * move the output by the full ramp 5, while a stable output does not move.
 * A modulation that the horizon cuts does not land, and its steps keep the
 * output in the band of origin, but it still leaves that band towards an
 * adjacent one: none starts
 * upwards from the top band, nor downwards from the bottom one, not even at
 * the last instant of the horizon or at the first one. The schedules that
 * do are written in the Variable, which then have to violate the rows of
 * the bands, while the mirror schedules, which move towards the inside of
 * the range, have to satisfy all of them; the dynamic programming Solver,
 * which never starts such a modulation, has to find the value of the unit
 * that keeps its output, whose schedule has to satisfy the rows as well.
 *
 * The test needs nothing but the core SMS++ (no MILP Solver): the rows are
 * checked on given points. */

/// the nuclear unit of test_nuclear_band_horizon(), initial power initP and
/// linear cost lin, over the horizon lin.size()

static NuclearUnitBlock * new_band_NU( const TUData & d , unsigned int L ,
				       const std::vector< double > & bands )
{
 auto g = new_group( "NU_bands" , true );
 write_TU( g , d );
 put_uint( g , "MaxModulationLength" , L );
 auto nb = g.addDim( "NumberPowerBands" , 2 );
 put( g , "PowerBands" , nb , bands );

 auto nub = dynamic_cast< NuclearUnitBlock * >( Block::new_Block( g ) );
 if( ! nub )
  throw( std::logic_error( "new_band_NU: the NuclearUnitBlock was not "
			   "built" ) );
 return( nub );
 }

static NuclearUnitBlock * new_banded_NU( double initP ,
					 const std::vector< double > & lin )
{
 TUData d;
 d.nuclear = true;
 d.T = lin.size();
 d.minP = 10;
 d.maxP = 50;
 d.ru = d.rd = 5;
 d.su = d.sd = 10;       // a unit above 10 cannot shut down
 d.initUD = 10;
 d.initP = initP;
 d.modT = 2;
 d.initMod = 2;          // free to modulate at 0
 d.mru = d.mrd = 0;      // a stable output does not move
 d.lin = lin;
 return( new_band_NU( d , 2 , { 20 , 40 } ) );
 }

/*--------------------------------------------------------------------------*/
/// writes the schedule ( p , m , d ), the unit being always on, in the
/// Variable of the unit, derives the others, and checks the rows

static bool banded_feasible( NuclearUnitBlock * nub ,
			     const std::vector< double > & p ,
			     const std::vector< double > & m ,
			     const std::vector< double > & d )
{
 auto P = nub->get_active_power( 0 );
 auto U = nub->get_commitment( 0 );
 auto M = nub->get_modulation();
 auto D = nub->get_modulation_down();
 if( ( ! P ) || ( ! U ) || ( ! M ) || ( ! D ) )
  throw( std::logic_error( "banded_feasible: Variable missing" ) );
 for( Index t = 0 ; t < p.size() ; ++t ) {
  P[ t ].set_value( p[ t ] );
  U[ t ].set_value( 1 );
  M[ t ].set_value( m[ t ] );
  D[ t ].set_value( d[ t ] );
  }
 nub->set_solution();
 SimpleConfiguration< double > tol( 1e-7 );
 return( nub->is_feasible( true , & tol ) );
 }

/*--------------------------------------------------------------------------*/

static void test_nuclear_band_horizon( void )
{
 // stable in the top band, a price at the last instant: going up by the
 // full ramp at 2 is a modulation out of the top band that the horizon
 // cuts (it stays in the top band, hence only the rows on the extreme
 // bands exclude it); going down by the full ramp, cut as well, is valid
 // as long as the output stays in the top band (from 46), and from 41 the
 // same step is the last one of a modulation that lands in the band below
 auto nub = new_banded_NU( 41 , { 0 , 0 , -10 } );
 generate_all( nub );
 check( ! banded_feasible( nub , { 41 , 41 , 46 } , { 0 , 0 , 1 } ,
			   { 0 , 0 , 0 } ) ,
	"nuclear bands: a step up out of the top band cut by the horizon is "
	"feasible" );
 check( banded_feasible( nub , { 41 , 41 , 36 } , { 0 , 0 , 1 } ,
			 { 0 , 0 , 1 } ) ,
	"nuclear bands: a full step down from the top band into the one below "
	"at the last instant is infeasible" );
 delete nub;
 nub = new_banded_NU( 46 , { 0 , 0 , -10 } );
 generate_all( nub );
 check( banded_feasible( nub , { 46 , 46 , 41 } , { 0 , 0 , 1 } ,
			 { 0 , 0 , 1 } ) ,
	"nuclear bands: a step down from the top band cut by the horizon is "
	"infeasible" );
 delete nub;
 // a modulation down to the breakpoint 40 at 0, which lands in the
 // intermediate band, and one up to 45 at 2: -450 (the MILP without the
 // rows on the extreme bands gets -460, going up from 41 at 2)
 nub = new_banded_NU( 41 , { 0 , 0 , -10 } );
 generate_all( nub );
 check_all_DP( nub , -450 , "nuclear bands, top at the end" );
 delete nub;

 // the mirror: stable in the bottom band, a cost at the last instant
 nub = new_banded_NU( 15 , { 0 , 0 , 10 } );
 generate_all( nub );
 check( ! banded_feasible( nub , { 15 , 15 , 10 } , { 0 , 0 , 1 } ,
			   { 0 , 0 , 1 } ) ,
	"nuclear bands: a step down out of the bottom band cut by the "
	"horizon is feasible" );
 check( banded_feasible( nub , { 15 , 15 , 20 } , { 0 , 0 , 1 } ,
			 { 0 , 0 , 0 } ) ,
	"nuclear bands: a step up from the bottom band cut by the horizon "
	"is infeasible" );
 delete nub;
 // the output cannot leave 15 downwards, nor shut down from above 10: 150
 // (the MILP without the rows gets 100)
 nub = new_banded_NU( 15 , { 0 , 0 , 10 } );
 generate_all( nub );
 check_all_DP( nub , 150 , "nuclear bands, bottom at the end" );
 delete nub;

 // a horizon of one instant, shorter than the longest modulation: the
 // band of the instant before is the band of the initial power, and the
 // modulation that starts at 0 is cut by the horizon
 nub = new_banded_NU( 41 , { -10 } );
 generate_all( nub );
 check( ! banded_feasible( nub , { 46 } , { 1 } , { 0 } ) ,
	"nuclear bands: a step up out of the top band at 0 is feasible" );
 check( banded_feasible( nub , { 36 } , { 1 } , { 1 } ) ,
	"nuclear bands: a full step down at 0 from the top band into the one "
	"below is infeasible" );
 delete nub;
 nub = new_banded_NU( 46 , { -10 } );
 generate_all( nub );
 check( banded_feasible( nub , { 41 } , { 1 } , { 1 } ) ,
	"nuclear bands: a step down from the top band at 0, cut by the "
	"horizon, is infeasible" );
 delete nub;
 nub = new_banded_NU( 41 , { -10 } );
 generate_all( nub );
 check_all_DP( nub , -410 , "nuclear bands, top at 0" );
 delete nub;

 nub = new_banded_NU( 15 , { 10 } );
 generate_all( nub );
 check( ! banded_feasible( nub , { 10 } , { 1 } , { 1 } ) ,
	"nuclear bands: a step down out of the bottom band at 0 is "
	"feasible" );
 check( banded_feasible( nub , { 20 } , { 1 } , { 0 } ) ,
	"nuclear bands: a step up from the bottom band at 0, cut by the "
	"horizon, is infeasible" );
 delete nub;
 }

/*--------------------------------------------------------------------------*/
/* The initial power of a banded NuclearUnitBlock changed after its abstract
 * representation is generated (T formulation): within the band of the
 * initial power the rows are updated, and the unit has to be the one a
 * fresh load gives (the dynamic programming Solver and the rows, which
 * check_DP() checks on its schedule); towards another band the rows of the
 * band at instant 0 cannot be updated, hence set_initial_power() has to
 * throw and leave the initial power as it was. */

static void test_nuclear_band_initial_power( void )
{
 // 41 -> 45, both in the top band [ 40 , 50 ]: the output stays at 45
 auto nub = new_banded_NU( 41 , { 0 , 0 , -10 } );
 generate_all( nub );
 std::vector< double > ip = { 45 };
 bool thrown = false;
 try {
  nub->set_initial_power( ip.cbegin() , Range( 0 , 1 ) );
  }
 catch( std::logic_error & ) { thrown = true; }
 check( ! thrown , "nuclear bands: set_initial_power() in the same band "
	"throws" );
 check( close( nub->get_initial_power() , 45 ) ,
	"nuclear bands: set_initial_power() in the same band is not applied" );
 check_all_DP( nub , -450 , "nuclear bands, initial power 41 -> 45" );
 delete nub;

 auto fresh = new_banded_NU( 45 , { 0 , 0 , -10 } );
 generate_all( fresh );
 check_all_DP( fresh , -450 , "nuclear bands, initial power 45" );
 delete fresh;

 // 41 -> 30, from the top band to the intermediate one
 nub = new_banded_NU( 41 , { 0 , 0 , -10 } );
 generate_all( nub );
 ip = { 30 };
 thrown = false;
 try {
  nub->set_initial_power( ip.cbegin() , Range( 0 , 1 ) );
  }
 catch( std::logic_error & ) { thrown = true; }
 check( thrown , "nuclear bands: set_initial_power() to another band does "
	"not throw" );
 check( close( nub->get_initial_power() , 41 ) ,
	"nuclear bands: a failed set_initial_power() changes the initial "
	"power" );
 check_all_DP( nub , -450 , "nuclear bands, failed initial power change" );
 delete nub;
 }

/*--------------------------------------------------------------------------*/
/*------------------ NUCLEAR: A MODULATION CROSSES ONE BAND ----------------*/
/*--------------------------------------------------------------------------*/
/* The optimum of a banded nuclear unit that is on before the horizon and
 * stays on, over the integer power levels, enumerated from the definition
 * of the operating rules: a stable instant moves the output by at most the
 * stability ramps within its band; a modulation starts from a stable
 * instant with no lockout, the steps that precede its last one move the
 * output by exactly the full ramp and keep it in the band of origin, the
 * last one moves it by at most the full ramp, in the same direction, into
 * the adjacent band, a modulation lasts at most L steps, and after its end
 * no step happens for ModulationTime - 1 instants. A modulation that the
 * horizon cuts is accepted. The band of the initial power is the lowest one
 * that contains it. Exact for linear costs and integer data (no deep
 * decrease, no daily limit, no cost of the downward steps). */

static double brute_force_bands( const TUData & d , unsigned int L ,
				 const std::vector< double > & bands )
{
 const int P0 = int( d.minP ) , P1 = int( d.maxP );
 const int np = P1 - P0 + 1;
 const int lo[ 3 ] = { P0 , int( bands[ 0 ] ) , int( bands[ 1 ] ) };
 const int hi[ 3 ] = { int( bands[ 0 ] ) , int( bands[ 1 ] ) , P1 };
 auto in = [ & ]( int p , int b ) { return( ( p >= lo[ b ] ) &&
					     ( p <= hi[ b ] ) ); };
 const int B = int( d.modT ) - 1;  // the lockout after the end
 // mode 0 stable, with lockout k in { 0 , ... , B }; mode 1 ( 2 ) in the
 // middle of an upward (downward) modulation, with k in { 1 , ... , L - 1 }
 // steps taken; b the band, or the band of origin in mode 1, 2
 const int nk = std::max( B + 1 , int( L ) );
 auto id = [ & ]( int p , int mode , int k , int b ) {
  return( ( ( ( p - P0 ) * 3 + mode ) * nk + k ) * 3 + b ); };
 const int ns = np * 3 * nk * 3;
 std::vector< double > f( ns , INF ) , g( ns );

 const int p_1 = int( d.initP );
 const int b_1 = ( p_1 <= hi[ 0 ] ) ? 0 : ( ( p_1 <= hi[ 1 ] ) ? 1 : 2 );
 const int l_1 = std::max( int( d.modT ) - int( d.initMod ) , 0 );
 f[ id( p_1 , 0 , l_1 , b_1 ) ] = 0;
 const int ru = int( d.ru ) , rd = int( d.rd );
 const int sru = int( std::min( d.ru , d.mru ) );
 const int srd = int( std::min( d.rd , d.mrd ) );

 for( Index t = 0 ; t < d.T ; ++t ) {
  std::fill( g.begin() , g.end() , INF );
  auto cost = [ & ]( int p ) {
   return( at( d.cnst , t ) + at( d.lin , t ) * p ); };
  auto relax = [ & ]( double fq , int p , int mode , int k , int b ) {
   if( ( p < P0 ) || ( p > P1 ) )
    return;
   const double v = fq + cost( p );
   if( v < g[ id( p , mode , k , b ) ] )
    g[ id( p , mode , k , b ) ] = v;
   };
  for( int q = P0 ; q <= P1 ; ++q )
   for( int mode = 0 ; mode < 3 ; ++mode )
    for( int k = 0 ; k < nk ; ++k )
     for( int b = 0 ; b < 3 ; ++b ) {
      const double fq = f[ id( q , mode , k , b ) ];
      if( fq == INF )
       continue;
      if( mode == 0 ) {
       // stable
       for( int p = q - srd ; p <= q + sru ; ++p )
	if( in( p , b ) )
	 relax( fq , p , 0 , std::max( k - 1 , 0 ) , b );
       if( k > 0 )
	continue;
       // a modulation of one step, up or down, landing next door
       if( b < 2 )
	for( int p = q ; p <= q + ru ; ++p )
	 if( in( p , b + 1 ) )
	  relax( fq , p , 0 , B , b + 1 );
       if( b > 0 )
	for( int p = q - rd ; p <= q ; ++p )
	 if( in( p , b - 1 ) )
	  relax( fq , p , 0 , B , b - 1 );
       // the first step of a longer one, staying in the band
       if( L > 1 ) {
	if( ( b < 2 ) && in( q + ru , b ) )
	 relax( fq , q + ru , 1 , 1 , b );
	if( ( b > 0 ) && in( q - rd , b ) )
	 relax( fq , q - rd , 2 , 1 , b );
	}
       continue;
       }
      // in the middle of a modulation: continue, or land next door
      const int dir = ( mode == 1 ) ? 1 : -1;
      const int r = ( mode == 1 ) ? ru : rd;
      if( ( k + 1 < int( L ) ) && in( q + dir * r , b ) )
       relax( fq , q + dir * r , mode , k + 1 , b );
      const int nb = b + dir;
      for( int s = 0 ; s <= r ; ++s )
       if( in( q + dir * s , nb ) )
	relax( fq , q + dir * s , 0 , B , nb );
      }
  std::swap( f , g );
  }

 double best = INF;
 for( auto v : f )
  best = std::min( best , v );
 return( best * d.scale );
 }

/*--------------------------------------------------------------------------*/
/* A modulation crosses exactly one boundary of its band: the steps before
 * the last one keep the output in the band of origin and the last one lands
 * in the adjacent band. Two cases: (a) the bands [ 230 , 248 ],
 * [ 248 , 432 ], [ 432 , 450 ] with the full ramps 15 down and 30 up, where
 * a modulation out of an extreme band has a single step; (b) the bands
 * [ 200 , 268 ], [ 268 , 382 ], [ 382 , 450 ] with the full ramp 30, halved
 * to [ 100 , 134 ], [ 134 , 191 ], [ 191 , 225 ] and 15 so that the brute
 * force stays small, where a modulation down from 150 is 150 -> 135 and a
 * last step into [ 120 , 134 ], while 150 -> 135 -> 120 -> [ 105 , 120 ],
 * whose second step leaves the band of origin, is not a modulation even if
 * MaxModulationLength allows 3 steps. Then random instances of both units,
 * with a large constant cost that keeps them on, against the brute force
 * of the operating rules. */

static void test_nuclear_band_rule( void )
{
 TUData a;
 a.nuclear = true;
 a.minP = 230;
 a.maxP = 450;
 a.ru = 30;
 a.rd = 15;
 a.su = a.sd = 230;
 a.initUD = 10;
 a.modT = 2;
 a.initMod = 2;
 a.mru = a.mrd = 0;
 const std::vector< double > ba = { 248 , 432 };

 // (a) from 440 in the top band, down: one step into [ 425 , 432 ]
 a.T = 2;
 a.initP = 440;
 a.lin = { 1 , 1 };
 auto nub = new_band_NU( a , 2 , ba );
 generate_all( nub );
 check( banded_feasible( nub , { 425 , 425 } , { 1 , 0 } , { 1 , 0 } ) ,
	"nuclear band rule (a): a one-step modulation 440 -> 425 is "
	"infeasible" );
 check( ! banded_feasible( nub , { 425 , 410 } , { 1 , 1 } , { 1 , 1 } ) ,
	"nuclear band rule (a): a modulation whose first step leaves the top "
	"band is feasible" );
 delete nub;
 nub = new_band_NU( a , 2 , ba );
 generate_all( nub );
 // 440 -> 425 at 0, locked at 1: 850
 check_all_DP( nub , brute_force_bands( a , 2 , ba ) ,
	       "nuclear band rule (a), down from 440" );
 check( close( brute_force_bands( a , 2 , ba ) , 850 ) ,
	"nuclear band rule (a): the brute force gives " +
	str( brute_force_bands( a , 2 , ba ) ) );
 delete nub;

 // (b) from 150 in the intermediate band, down, with L^M = 3
 TUData b = a;
 b.minP = 100;
 b.maxP = 225;
 b.ru = b.rd = 15;
 b.su = b.sd = 100;
 b.T = 3;
 b.initP = 150;
 b.lin = { 1 , 1 , 1 };
 const std::vector< double > bb = { 134 , 191 };
 nub = new_band_NU( b , 3 , bb );
 generate_all( nub );
 check( banded_feasible( nub , { 135 , 125 , 125 } , { 1 , 1 , 0 } ,
			 { 1 , 1 , 0 } ) ,
	"nuclear band rule (b): 150 -> 135 -> 125 is infeasible" );
 check( ! banded_feasible( nub , { 135 , 120 , 110 } , { 1 , 1 , 1 } ,
			   { 1 , 1 , 1 } ) ,
	"nuclear band rule (b): 150 -> 135 -> 120 -> 110 is feasible" );
 delete nub;
 nub = new_band_NU( b , 3 , bb );
 generate_all( nub );
 // 150 -> 135 -> 120 (landing at the end of the full ramp), locked at 2:
 // 135 + 120 + 120 = 375
 check( close( brute_force_bands( b , 3 , bb ) , 375 ) ,
	"nuclear band rule (b): the brute force gives " +
	str( brute_force_bands( b , 3 , bb ) ) );
 check_all_DP( nub , 375 , "nuclear band rule (b), down from 150" );
 delete nub;

 // random instances of the two units
 for( int i = 0 ; i < 60 ; ++i ) {
  TUData r = ( i % 2 ) ? b : a;
  const auto & bands = ( i % 2 ) ? bb : ba;
  r.T = rnd( 2 , 5 );
  const unsigned int L = rnd( 2 , 4 );
  r.modT = rnd( 2 , 3 );
  r.initMod = rnd( 1 , int( r.modT ) );
  r.mru = rnd( 0 , 2 );
  r.mrd = rnd( 0 , 2 );
  r.initP = rnd( int( r.minP ) , int( r.maxP ) );
  r.lin.resize( r.T );
  r.cnst.assign( r.T , -1e5 );       // shutting down never pays
  for( auto & c : r.lin )
   c = rnd( -3 , 3 );
  const auto what = "random band rule " + std::to_string( i ) + " (L = " +
                    std::to_string( L ) + ", " + describe( r ) + ")";
  nub = new_band_NU( r , L , bands );
  generate_all( nub );
  check_all_DP( nub , brute_force_bands( r , L , bands ) , what );
  delete nub;
  }
 }

/*--------------------------------------------------------------------------*/
/*------------- RESERVES, AVAILABILITY AND THE CHANGES OF DATA -------------*/
/*--------------------------------------------------------------------------*/
/// the first :MILPSolver of the build, empty if there is none

static std::string milp_solver( void )
{
 for( const auto & name : { "CPXMILPSolver" , "GRBMILPSolver" ,
                            "HiGHSMILPSolver" , "SCIPMILPSolver" } )
  if( Solver::has_Solver( name ) )
   return( name );
 return( "" );
 }

/*--------------------------------------------------------------------------*/
/// the value of the Block, with its abstract representation, that the
/// :MILPSolver sname finds: INF if infeasible, NaN if it finds none

static double milp_value( Block * b , const std::string & sname )
{
 auto slv = Solver::new_Solver( sname );
 if( ! slv )
  return( std::nan( "" ) );
 b->register_Solver( slv );
 double v = std::nan( "" );
 try {
  const auto status = slv->compute();
  if( status == Solver::kOK )
   v = slv->get_var_value();
  else
   if( status == Solver::kInfeasible )
    v = INF;
  }
 catch( std::exception & e ) {
  std::cout << sname << " throws " << e.what() << std::endl;
  }
 b->unregister_Solver( slv , true );
 return( v );
 }

/*--------------------------------------------------------------------------*/
/// whether the Variable v is in the Objective of the unit

static bool in_obj( ThermalUnitBlock * tub , const ColVariable * v )
{
 auto qf = static_cast< DQuadFunction * >(
	   static_cast< FRealObjective * >( tub->get_objective() )
                                                         ->get_function() );
 return( qf->is_active( v ) < qf->get_num_active_var() );
 }

/*--------------------------------------------------------------------------*/
/// every DP Solver of the unit, generated in any formulation (the last
/// argument), against the expected value, with the schedule they write in
/// the Variable, which ThermalUnitBlock::set_solution() completes in every
/// formulation

static void check_all_DP_wf( ThermalUnitBlock * tub , double expected ,
                             const std::string & what , int )
{
 check_all_DP( tub , expected , what );
 }

/*--------------------------------------------------------------------------*/
/// every DP Solver of the unit against the expected value, the value alone
/// (for a unit whose rows are not those of its data, e.g., after a change
/// whose abstract part is a dry run)

static void check_all_DP_value( ThermalUnitBlock * tub , double expected ,
                                const std::string & what )
{
 std::vector< std::string > names = THERMAL_DP;
 if( dynamic_cast< NuclearUnitBlock * >( tub ) )
  names = { "NuclearUnitExtDPSolver" };
 for( const auto & s : names ) {
  auto slv = Solver::new_Solver( s );
  tub->register_Solver( slv );
  const auto status = slv->compute();
  check( ( status == Solver::kOK ) && close( slv->get_ub() , expected ) ,
         what + ", " + s + ": status " + std::to_string( status ) +
         ", value " + str( slv->get_ub() ) + " instead of " +
         str( expected ) );
  tub->unregister_Solver( slv , true );
  }
 }

/*--------------------------------------------------------------------------*/
/* The cost of the reserves of a ThermalUnitBlock is 0 when the data do not
 * give it, PrimaryRho and SecondaryRho being only the fractions: a reserve
 * whose cost is 0 is in the Objective only if the Configuration of the
 * Objective (the parameter of generate_objective(), or else the one of the
 * BlockConfig) asks for it, bit 0 the primary and bit 1 the secondary, and
 * with a nonzero cost it is there anyway. A cost changed after the
 * Objective is generated without the term is refused, unless it is 0;
 * zeros set before and after the generation give the same Objective; and a
 * cost of the secondary reserve lands on the secondary Variable also when
 * the shut-down section moves the reserve sections. */

static void test_thermal_reserve_cost( void )
{
 const Index T = 3;
 auto unit = [ & ]( const std::vector< double > & prc ) {
  std::vector< std::pair< std::string , std::vector< double > > > v = {
   { "MinPower" , { 10 , 10 , 10 } } , { "MaxPower" , { 100 , 100 , 100 } } ,
   { "LinearTerm" , { 1 , 1 , 1 } } , { "ShutDownCost" , { 5 , 5 , 5 } } ,
   { "PrimaryRho" , { 0.1 , 0.1 , 0.1 } } ,
   { "SecondaryRho" , { 0.2 , 0.2 , 0.2 } } };
  if( ! prc.empty() )
   v.push_back( { "PrimarySpinningReserveCost" , prc } );
  auto tub = new_TU_vec( T , v , 2 , 50 , 1 , 1 );
  tub->set_reserve_vars( 3 );
  tub->generate_abstract_variables();
  tub->generate_abstract_constraints();
  return( tub );
  };

 auto coefs = []( ThermalUnitBlock * tub , const ColVariable * v ) {
  std::vector< double > c;
  for( Index t = 0 ; t < T ; ++t )
   c.push_back( in_obj( tub , v + t ) ? lin_coef( tub , v + t )
                                      : std::nan( "" ) );
  return( c );
  };

 // no Configuration, no cost: neither reserve in the Objective - - - - - -
 {
  auto tub = unit( { } );
  tub->generate_objective();
  check( ( ! in_obj( tub , tub->get_primary_spinning_reserve( 0 ) ) ) &&
         ( ! in_obj( tub , tub->get_secondary_spinning_reserve( 0 ) ) ) ,
         "reserve cost, no cost and no Configuration: a reserve term in "
         "the Objective" );
  check( tub->get_primary_spinning_reserve_cost().empty() ,
         "reserve cost, no cost: the cost is not 0" );
  // zeros are accepted, and a nonzero cost refused with nothing changed
  std::vector< double > z( T , 0 ) , c( T , 2 );
  bool thrown = false;
  try {
   tub->set_primary_spinning_reserve_cost( z.begin() , Range( 0 , T ) );
   }
  catch( std::exception & ) {
   thrown = true;
   }
  check( ! thrown , "reserve cost, no term: zeros are refused" );
  thrown = false;
  try {
   tub->set_primary_spinning_reserve_cost( c.begin() , Range( 0 , T ) );
   }
  catch( std::logic_error & ) {
   thrown = true;
   }
  const auto & pc = tub->get_primary_spinning_reserve_cost();
  check( thrown && std::all_of( pc.begin() , pc.end() ,
                                []( double v ) { return( v == 0 ); } ) ,
         "reserve cost, no term: a nonzero cost is not refused, or it "
         "changes the data" );
  delete tub;
  }

 // the parameter asks for the primary alone, then for both - - - - - - - -
 for( int bits : { 1 , 3 } ) {
  auto tub = unit( { } );
  SimpleConfiguration< int > objc( bits );
  tub->generate_objective( & objc );
  const auto pr = coefs( tub , tub->get_primary_spinning_reserve( 0 ) );
  const auto sr = coefs( tub , tub->get_secondary_spinning_reserve( 0 ) );
  check( ( pr == std::vector< double >( T , 0 ) ) &&
         ( ( bits == 3 ) ? ( sr == std::vector< double >( T , 0 ) )
                         : std::isnan( sr[ 0 ] ) ) ,
         "reserve cost, Configuration " + std::to_string( bits ) +
         ": the reserve terms are not those asked for, with cost 0" );
  if( bits == 3 ) {
   // the secondary cost after the shut-down and the primary sections
   std::vector< double > c2( T , 4 ) , c1 = { 0 , 2 , 0 };
   tub->set_secondary_spinning_reserve_cost( c2.begin() , Range( 0 , T ) );
   tub->set_primary_spinning_reserve_cost( c1.begin() , Range( 0 , T ) );
   const auto pr2 = coefs( tub , tub->get_primary_spinning_reserve( 0 ) );
   const auto sr2 = coefs( tub , tub->get_secondary_spinning_reserve( 0 ) );
   const auto sd = coefs( tub , tub->get_shut_down() );
   const auto p = coefs( tub , tub->get_active_power( 0 ) );
   check( ( pr2 == c1 ) && ( sr2 == c2 ) &&
          ( sd == std::vector< double >( T , 5 ) ) &&
          ( p == std::vector< double >( T , 1 ) ) ,
          "reserve cost, both terms: the costs set after the generation "
          "are not on the reserve Variable, or they move other costs" );
   }
  delete tub;
  }

 // the BlockConfig asks for both - - - - - - - - - - - - - - - - - - - - -
 {
  auto tub = unit( { } );
  std::istringstream text( "BlockConfig 1 2 * * * * * "
                           "SimpleConfiguration<int> 3 * * * * " );
  auto bc = dynamic_cast< BlockConfig * >(
                                     Configuration::deserialize( text ) );
  check( bc , "reserve cost, BlockConfig: not read" );
  if( bc ) {
   bc->apply( tub );
   delete bc;
   tub->generate_objective();
   check( in_obj( tub , tub->get_primary_spinning_reserve( 0 ) ) &&
          in_obj( tub , tub->get_secondary_spinning_reserve( 0 ) ) ,
          "reserve cost, BlockConfig: the reserve terms are not there" );
   }
  delete tub;
  }

 // a nonzero cost in the data: the primary is there, with its cost - - - -
 {
  auto tub = unit( { 0 , 2 , 0 } );
  tub->generate_objective();
  check( ( coefs( tub , tub->get_primary_spinning_reserve( 0 ) ) ==
           std::vector< double >( { 0 , 2 , 0 } ) ) &&
         ( ! in_obj( tub , tub->get_secondary_spinning_reserve( 0 ) ) ) ,
         "reserve cost, cost given: the terms are not the primary alone "
         "with its cost" );
  delete tub;
  }

 // zeros before and after the generation give the same Objective - - - - -
 {
  std::vector< std::vector< double > > got;
  for( int after = 0 ; after < 2 ; ++after ) {
   auto tub = unit( { } );
   SimpleConfiguration< int > objc( 1 );
   std::vector< double > z( T , 0 );
   if( ! after )
    tub->set_primary_spinning_reserve_cost( z.begin() , Range( 0 , T ) );
   tub->generate_objective( & objc );
   if( after )
    tub->set_primary_spinning_reserve_cost( z.begin() , Range( 0 , T ) );
   got.push_back( coefs( tub , tub->get_primary_spinning_reserve( 0 ) ) );
   delete tub;
   }
  check( got[ 0 ] == got[ 1 ] ,
         "reserve cost, zeros before and after the generation: different "
         "Objective" );
  }
 }

/*--------------------------------------------------------------------------*/
/// the data of a unit with reserves and ramps whose instant 0 the ramps
/// from InitialPower bind; random if i > 0, the case below if i == 0

static TUData reserve_TU( int i )
{
 TUData d;
 if( i == 0 ) {  // on at 20 MW, ramps 2, every MW of output worth 1
  d.T = 3;
  d.minP = 13;
  d.maxP = 54;
  d.ru = d.rd = 2;
  d.su = 43;
  d.sd = 44;
  d.initP = 20;
  d.initUD = 100;
  d.lin.assign( d.T , -1 );
  d.suc.assign( d.T , 1000 );  // no shut-down and start-up at 43 MW
  d.prho.assign( d.T , 0.05 );
  d.srho.assign( d.T , 0.08 );
  d.prc.assign( d.T , -0.05 );
  d.src.assign( d.T , -0.05 );
  return( d );
  }

 d.T = rnd( 2 , 5 );
 d.minP = 10;
 d.maxP = rnd( 30 , 60 );
 d.ru = rnd( 2 , 15 );
 d.rd = rnd( 2 , 15 );
 d.su = rnd( 10 , int( d.maxP ) );
 d.sd = rnd( 10 , int( d.maxP ) );
 d.initUD = rnd( -2 , 3 );
 if( d.initUD == 0 )
  d.initUD = 1;
 d.initP = rnd( 10 , int( d.maxP ) );
 d.minUp = rnd( 1 , 2 );
 d.minDown = rnd( 1 , 2 );
 for( Index t = 0 ; t < d.T ; ++t ) {
  d.lin.push_back( rnd( -30 , 30 ) / 10.0 );
  d.cnst.push_back( rnd( 0 , 5 ) );
  d.suc.push_back( rnd( 0 , 20 ) );
  d.prho.push_back( rnd( 5 , 30 ) / 100.0 );
  d.srho.push_back( rnd( 5 , 30 ) / 100.0 );
  // a reserve worth more than the output, often
  d.prc.push_back( - rnd( 0 , 40 ) / 10.0 );
  d.src.push_back( - rnd( 0 , 40 ) / 10.0 );
  }
 return( d );
 }

/*--------------------------------------------------------------------------*/
/* The rows of the reserve deliverability at instant 0 of a unit on before
 * the horizon: they bound the reserve by the ramp left after the move from
 * InitialPower, in every formulation, and the dynamic programming Solvers
 * do the same. The instance on at 20 MW with ramps 2, whose output goes up
 * by the full ramp, has no reserve at 0 (value -72, while without the rows
 * the reserve at 0 earns 0.143); random instances with reserves worth more
 * than the output (on and off before the horizon, ramps that bind at 0) give
 * the same value to the :MILPSolver in the seven formulations and to the
 * dynamic programming Solvers, whose schedule satisfies the rows; a unit
 * that shuts down at 0 leaves the row slack; the rows at 0 follow a
 * change of InitialPower; and a NuclearUnitBlock with reserves gives the
 * same value to the :MILPSolver and to NuclearUnitExtDPSolver when its
 * primary reserve is rewarded. */

static void test_thermal_reserve_t0( void )
{
 const auto sname = milp_solver();
 auto rows = []( ThermalUnitBlock * tub ) {
  return( tub->get_static_constraint_v< FRowConstraint >(
                                              "Reserve_Const_Thermal" ) );
  };
 auto make = []( const TUData & d , int wf ) {
  auto tub = new_TU( d );
  tub->set_reserve_vars( 3 );
  generate_all_wf( tub , wf );
  return( tub );
  };

 // the instance on at 20 MW: no reserve at 0, both DPs and the MILP - - - -
 for( int wf = 0 ; wf < 7 ; ++wf ) {
  const auto d = reserve_TU( 0 );
  auto tub = make( d , wf );
  const auto what = "reserve at 0, on at 20 MW, wf " + std::to_string( wf );
  auto r = rows( tub );
  check( r && ( r->size() == 5 * d.T ) ,
         what + ": not 5 T reserve rows" );
  if( ! sname.empty() ) {
   const auto v = milp_value( tub , sname );
   check( close( v , -72 ) , what + ": " + sname + " gives " + str( v ) +
          " instead of -72" );
   }
  check_all_DP_wf( tub , -72 , what , wf );
  delete tub;
  }

 // off before the horizon: no row at 0 - - - - - - - - - - - - - - - - - -
 {
  auto d = reserve_TU( 0 );
  d.initUD = -1;
  auto tub = make( d , ThermalUnitBlock::TForm );
  auto r = rows( tub );
  check( r && ( r->size() == 5 * d.T - 2 ) ,
         "reserve at 0, off before: not 5 T - 2 reserve rows" );
  delete tub;
  }

 // on before, positive prices: shut down at 0, the row (5) at 0 slack - - -
 {
  auto d = reserve_TU( 0 );
  d.initUD = 3;
  d.initP = 40;
  d.lin.assign( d.T , 1 );
  for( int wf = 0 ; wf < 7 ; ++wf ) {
   auto tub = make( d , wf );
   const auto what = "reserve at 0, shut-down at 0, wf " +
                     std::to_string( wf );
   if( ! sname.empty() ) {
    const auto v = milp_value( tub , sname );
    check( close( v , 0 ) , what + ": " + sname + " gives " + str( v ) +
           " instead of 0" );
    }
   check_all_DP_wf( tub , 0 , what , wf );
   delete tub;
   }
  }

 // InitialPower changed after the generation, 3bin and T formulations- - -
 for( int wf : { int( ThermalUnitBlock::tbinForm ) ,
                 int( ThermalUnitBlock::TForm ) } ) {
  auto d = reserve_TU( 0 );
  auto tub = make( d , wf );
  std::vector< double > ip = { 30 };
  tub->set_initial_power( ip.begin() , Range( 0 , 1 ) );
  const auto what = "reserve at 0, InitialPower 20 -> 30, wf " +
                    std::to_string( wf );
  auto r = rows( tub );
  if( r && ( r->size() == 5 * d.T ) ) {
   auto & r4 = ( *r )[ 3 * d.T ];
   auto & r5 = ( *r )[ 4 * d.T ];
   check( close( r4.get_rhs() , 32 ) ,
          what + ": the right-hand side of the ramp-up row at 0 is " +
          str( r4.get_rhs() ) );
   check( close( row_coef( r5 , tub->get_commitment( 0 ) ) , 28 ) ,
          what + ": the commitment at 0 has " +
          str( row_coef( r5 , tub->get_commitment( 0 ) ) ) +
          " in the ramp-down row at 0" );
   }
  else
   check( false , what + ": not 5 T reserve rows" );
  d.initP = 30;
  auto fresh = make( d , wf );
  if( ! sname.empty() ) {
   const auto v = milp_value( tub , sname );
   const auto vf = milp_value( fresh , sname );
   check( close( v , vf ) , what + ": " + sname + " gives " + str( v ) +
          ", a fresh load " + str( vf ) );
   check_all_DP( tub , vf , what );
   }
  delete fresh;
  delete tub;
  }

 if( sname.empty() ) {
  std::cout << "no :MILPSolver in this build, the reserves at 0 of random "
            << "instances are not checked" << std::endl;
  return;
  }

 // random instances, all the formulations in turn- - - - - - - - - - - - - -
 for( int i = 1 ; i <= 60 ; ++i ) {
  const auto d = reserve_TU( i );
  for( int wf : { 1 , i % 7 } ) {
   auto tub = make( d , wf );
   const auto what = "reserve at 0, random " + std::to_string( i ) +
                     ", wf " + std::to_string( wf ) + " (" + describe( d ) +
                     ")";
   const auto v = milp_value( tub , sname );
   if( std::isnan( v ) )
    check( false , what + ": " + sname + " finds nothing" );
   else
    check_all_DP_wf( tub , v , what , wf );
   delete tub;
   }
  }

 // a NuclearUnitBlock with reserves, on before the horizon - - - - - - - - -
 for( int i = 0 ; i < 10 ; ++i ) {
  auto d = reserve_TU( i + 1 );
  d.nuclear = true;
  d.initUD = 3;
  d.mru = d.ru;
  d.mrd = d.rd;
  d.src.assign( d.T , 0 );  // the primary reserve alone is rewarded
  auto tub = make( d , ThermalUnitBlock::TForm );
  const auto what = "reserve at 0, nuclear " + std::to_string( i ) + " (" +
                    describe( d ) + ")";
  const auto v = milp_value( tub , sname );
  if( std::isnan( v ) )
   check( false , what + ": " + sname + " finds nothing" );
  else
   check_all_DP( tub , v , what );
  delete tub;
  }
 }

/*--------------------------------------------------------------------------*/
/* The on-to-on transition of the dynamic programming Solvers when the
 * rewarded reserves fill a fraction kappa = rho^pr + rho^sc (or a
 * single rho) of the output larger than 1/2 and the ramp-down
 * deliverability row binds: the minimiser q of (4) of
 * ThermalUnitDPSolverBase.h then moves with the landing power p by
 * 1 - kappa < 1/2. The unit on at 27 MW before a horizon of 2 instants,
 * with no energy cost and both reserves rewarded (kappa = 0.51 at 1), has
 * its optimum at the outputs 25 and 25.5, where the rows (30), (31) and
 * (32) of ThermalUnitBlock.h all leave 12.5 MW of reserve at 1: the value
 * is 7 - 2.8 * 7 - 0.2 * 1.75 - 3.7 * 6.63 - 2.2 * 5.87 = -50.395. The same
 * transition away from instant 0 (an instant with no reserve price in
 * front, -48.395), with the secondary reserve the more rewarded one
 * (ordering swapped), with the primary reserve alone (rho^pr = 0.51 at 1, cost
 * -2.2: 7 - 19.95 - 2.2 * 12.5 = -40.45), and for a NuclearUnitBlock with
 * the modulation window of the ramps and a narrower one (-50.395, the move
 * being 0.5); then random instances whose rewarded reserves are often
 * above one half of the output, against the :MILPSolver. */

static void test_thermal_DP_two_rewarded_reserves( void )
{
 const auto sname = milp_solver();
 const auto rg0 = rg;  // the random instances of the other tests unchanged
 auto t2 = []( void ) {
  TUData d;
  d.T = 2;
  d.minP = 10;
  d.maxP = 38;
  d.ru = 13;
  d.rd = 12;
  d.su = 22;
  d.sd = 20;
  d.initP = 27;
  d.initUD = 3;
  d.minUp = 1;
  d.minDown = 2;
  d.cnst = { 2 , 5 };
  d.suc = { 16 , 20 };
  d.prho = { 0.28 , 0.26 };
  d.srho = { 0.07 , 0.25 };
  d.prc = { -2.8 , -3.7 };
  d.src = { -0.2 , -2.2 };
  return( d );
  };
 auto run = [ & ]( const TUData & d , double expected ,
                   const std::string & what , int wf ) {
  auto tub = new_TU( d );
  tub->set_reserve_vars( 3 );
  generate_all_wf( tub , wf );
  const auto who = what + ", wf " + std::to_string( wf );
  double v = expected;
  if( ! sname.empty() ) {
   v = milp_value( tub , sname );
   if( ! std::isnan( expected ) )
    check( close( v , expected ) , who + ": " + sname + " gives " +
           str( v ) + " instead of " + str( expected ) );
   }
  if( ! std::isnan( v ) )
   check_all_DP_wf( tub , v , who , wf );
  else
   if( ! sname.empty() )
    check( false , who + ": " + sname + " finds nothing" );
  delete tub;
  };

 for( int wf : { int( ThermalUnitBlock::tbinForm ) ,
                 int( ThermalUnitBlock::TForm ) } ) {
  run( t2() , -50.395 , "two rewarded reserves, T = 2" , wf );

  auto d = t2();  // an instant with no reserve price in front
  d.T = 3;
  d.cnst = { 2 , 2 , 5 };
  d.suc = { 16 , 16 , 20 };
  d.prho = { 0.28 , 0.28 , 0.26 };
  d.srho = { 0.07 , 0.07 , 0.25 };
  d.prc = { 0 , -2.8 , -3.7 };
  d.src = { 0 , -0.2 , -2.2 };
  run( d , -48.395 , "two rewarded reserves, T = 3" , wf );

  d = t2();       // the secondary reserve is the more rewarded one at 1
  d.prc[ 1 ] = -1.5;
  run( d , std::nan( "" ) , "two rewarded reserves, swapped" , wf );

  d = t2();       // the primary reserve alone, half of the output
  d.prho[ 1 ] = 0.51;
  d.prc[ 1 ] = -2.2;
  d.src[ 1 ] = 0;
  run( d , -40.45 , "one rewarded reserve, rho 0.51" , wf );

  if( wf != int( ThermalUnitBlock::TForm ) )
   continue;
  d = t2();       // nuclear, the window of the ramps and a narrower one
  d.nuclear = true;
  d.modT = 2;
  d.initMod = 2;
  d.mru = 13;
  d.mrd = 12;
  run( d , -50.395 , "two rewarded reserves, nuclear" , wf );
  d.mrd = 6;
  run( d , -50.395 , "two rewarded reserves, nuclear, window 6" , wf );
  }

 if( sname.empty() ) {
  std::cout << "no :MILPSolver in this build, random instances with two "
            << "rewarded reserves are not checked" << std::endl;
  return;
  }

 for( int i = 1 ; i <= 40 ; ++i ) {
  auto d = reserve_TU( i );
  for( Index t = 0 ; t < d.T ; ++t ) {
   d.prho[ t ] = rnd( 20 , 60 ) / 100.0;
   d.srho[ t ] = rnd( 5 , 40 ) / 100.0;
   }
  if( i > 30 ) {  // nuclear, on before the horizon
   d.nuclear = true;
   d.initUD = 3;
   d.mru = d.ru;
   d.mrd = ( i % 2 ) ? d.rd : std::max( 1.0 , d.rd / 2 );
   }
  run( d , std::nan( "" ) , "two rewarded reserves, random " +
       std::to_string( i ) + " (" + describe( d ) + ")" ,
       d.nuclear ? int( ThermalUnitBlock::TForm ) : 1 + i % 6 );
  }
 rg = rg0;
 }

/*--------------------------------------------------------------------------*/
/* The dynamic programming Solvers use the operational bounds of the rows,
 * Availability included, and, without ramp data, a ramp that no move
 * reaches: an availability 0.6 at 1 with ramps 30 (-3300), an availability
 * 0 at 1 that keeps the unit on at no output (-1997), and a maximum power
 * 30 at 2 with no ramp data, which the output reaches from 100 at once
 * (-3300); the :MILPSolver gives the same values in the seven
 * formulations. A change of the availability after a solve that is only
 * physical (the abstract one being a dry run) reaches the Solvers, which
 * then give the value of a fresh load. */

static void test_thermal_DP_availability( void )
{
 const auto sname = milp_solver();
 struct Case {
  std::string what;
  std::vector< std::pair< std::string , std::vector< double > > > data;
  double value;
  };
 const std::vector< Case > cases = {
  { "availability 0.6 at 1, ramps 30" ,
    { { "MinPower" , { 10 , 10 , 10 , 10 } } ,
      { "MaxPower" , { 100 , 100 , 100 , 100 } } ,
      { "DeltaRampUp" , { 30 , 30 , 30 , 30 } } ,
      { "DeltaRampDown" , { 30 , 30 , 30 , 30 } } ,
      { "LinearTerm" , { -10 , -10 , -10 , -10 } } ,
      { "Availability" , { 1 , 0.6 , 1 , 1 } } } , -3300 } ,
  { "availability 0 at 1" ,
    { { "MinPower" , { 10 , 10 , 10 } } ,
      { "MaxPower" , { 100 , 100 , 100 } } ,
      { "LinearTerm" , { -10 , -10 , -10 } } ,
      { "ConstTerm" , { 1 , 1 , 1 } } ,
      { "StartUpCost" , { 200 , 200 , 200 } } ,
      { "Availability" , { 1 , 0 , 1 } } } , -1997 } ,
  { "maximum power 30 at 2, no ramp data" ,
    { { "MinPower" , { 10 , 10 , 10 , 10 } } ,
      { "MaxPower" , { 100 , 100 , 30 , 100 } } ,
      { "LinearTerm" , { -10 , -10 , -10 , -10 } } } , -3300 } };

 for( const auto & c : cases ) {
  const Index T = c.data.front().second.size();
  const double ip = ( c.what[ 0 ] == 'm' ) ? 100 : 50;
  for( int wf = 0 ; wf < 7 ; ++wf ) {
   auto tub = new_TU_vec( T , c.data , 2 , ip , 1 , 1 );
   generate_all_wf( tub , wf );
   const auto what = c.what + ", wf " + std::to_string( wf );
   if( ! sname.empty() ) {
    const auto v = milp_value( tub , sname );
    check( close( v , c.value ) , what + ": " + sname + " gives " +
           str( v ) + " instead of " + str( c.value ) );
    }
   check_all_DP_wf( tub , c.value , what , wf );
   delete tub;
   }
  }

 // the availability changed, the abstract change being a dry run - - - - -
 auto data = cases[ 0 ].data;
 data.pop_back();  // no Availability
 auto tub = new_TU_vec( 4 , data , 2 , 50 , 1 , 1 );
 generate_all_wf( tub , ThermalUnitBlock::TForm );
 check_all_DP( tub , -3800 , "availability set by a dry run, before" );
 std::vector< double > av = { 0.6 };
 tub->set_availability( av.begin() , Range( 1 , 2 ) , eModBlck , eDryRun );
 check_all_DP( tub , -3300 , "availability set by a dry run" );
 delete tub;
 }

/*--------------------------------------------------------------------------*/
/* The maximum power and the availability of a ThermalUnitBlock changed once
 * the Constraint are generated, in every formulation, with and without ramp
 * data: the change is made, and the dynamic programming Solver give the
 * value of a fresh load with a schedule that satisfies the rows as they
 * are after the change [see check_DP()], while a change whose abstract
 * part is a dry run leaves the rows as they were, and the value alone is
 * that of a fresh load. Values that a fresh load refuses (a maximum power
 * below the minimum one, an availability outside [ 0 , 1 ]) are refused,
 * and the data stay as they were. MaxRampUpSteps and MaxRampDownSteps
 * given in the data are refused. Every formulation, compared row by row
 * with the unit read afresh, is in test_power_limits(). */

static void test_thermal_power_limits_change( void )
{
 for( int ramps = 0 ; ramps < 2 ; ++ramps )
  for( int wf = 0 ; wf < 7 ; ++wf ) {
   std::vector< std::pair< std::string , std::vector< double > > > data = {
    { "MinPower" , { 10 , 10 , 10 , 10 } } ,
    { "MaxPower" , { 100 , 100 , 100 , 100 } } ,
    { "LinearTerm" , { -10 , -10 , -10 , -10 } } };
   if( ramps ) {
    data.push_back( { "DeltaRampUp" , { 30 , 30 , 30 , 30 } } );
    data.push_back( { "DeltaRampDown" , { 30 , 30 , 30 , 30 } } );
    }
   auto tub = new_TU_vec( 4 , data , 2 , 100 , 1 , 1 );
   generate_all_wf( tub , wf );
   const auto what = "maximum power change, wf " + std::to_string( wf ) +
                     ( ramps ? ", ramps" : ", no ramps" );

   // 150 at 1: the ramps let it reach 130 only
   std::vector< double > mp = { 150 } , av = { 0.6 } , bad = { 5 };
   tub->set_maximum_power( mp.begin() , Range( 1 , 2 ) );
   check( tub->get_max_power()[ 1 ] == 150 ,
          what + ": MaxPower not changed after the generation" );
   check_all_DP_wf( tub , ramps ? -4300 : -4500 , what , wf );

   // below the minimum power: refused
   bool thrown = false;
   try {
    tub->set_maximum_power( bad.begin() , Range( 2 , 3 ) );
    }
   catch( std::logic_error & ) {
    thrown = true;
    }
   check( thrown && ( tub->get_max_power()[ 2 ] == 100 ) ,
          what + ": MaxPower below MinPower accepted, or the datum changed" );

   // the availability 0.6 at 1, i.e., 90 there, without ramps (with them
   // the unit can neither go down to 90 nor shut down there)
   if( ! ramps ) {
    tub->set_availability( av.begin() , Range( 1 , 2 ) );
    check( ( ! tub->get_availability().empty() ) &&
           ( tub->get_availability()[ 1 ] == 0.6 ) ,
           what + ": Availability not changed after the generation" );
    check_all_DP_wf( tub , -3900 , what + ", availability" , wf );
    }
   std::vector< double > av_bad = { 1.5 };
   thrown = false;
   try {
    tub->set_availability( av_bad.begin() , Range( 2 , 3 ) );
    }
   catch( std::logic_error & ) {
    thrown = true;
    }
   check( thrown && ( tub->get_availability().empty() ||
                      ( tub->get_availability()[ 2 ] == 1 ) ) ,
          what + ": Availability 1.5 accepted, or the datum changed" );

   if( ! ramps ) {
    // the rows are those of the old MaxPower, which a dry run leaves:
    // the value alone
    std::vector< double > mp2 = { 200 };
    tub->set_maximum_power( mp2.begin() , Range( 2 , 3 ) , eModBlck ,
                            eDryRun );
    check( tub->get_max_power()[ 2 ] == 200 ,
           what + ": MaxPower not changed by a dry run" );
    check_all_DP_value( tub , -4900 , what + ", dry run" );
    }
   delete tub;
   }

 for( const std::string n : { "MaxRampUpSteps" , "MaxRampDownSteps" } ) {
  auto g = new_group( "TUs" , true );
  g.putAtt( "type" , "ThermalUnitBlock" );
  g.addDim( "TimeHorizon" , 3 );
  auto NI = g.addDim( "NumberIntervals" , 3 );
  put( g , "MinPower" , NI , { 10 , 10 , 10 } );
  put( g , "MaxPower" , NI , { 100 , 100 , 100 } );
  put( g , "DeltaRampUp" , NI , { 30 , 30 , 30 } );
  put( g , "DeltaRampDown" , NI , { 30 , 30 , 30 } );
  put_int( g , n , NI , { 1 , 1 , 1 } );
  put_int( g , "InitUpDownTime" , -1 );
  bool refused = false;
  Block * b = nullptr;
  try {
   b = Block::new_Block( g );
   refused = ( b == nullptr );
   }
  catch( std::exception & ) {
   refused = true;
   }
  check( refused , n + " given in the data: accepted" );
  delete b;
  }
 }

/*--------------------------------------------------------------------------*/
/* The reactive power of a ThermalUnitBlock may be negative (absorbed), and
 * an absent bound of it is 0: with MinReactivePower -20, MaxReactivePower
 * 30 and a reactive term 1 the reactive power is -20 (value -3060), and with
 * no reactive data and a term -1 it is 0 (-3000, unbounded if the bounds
 * were not written), for the :MILPSolver in the seven formulations as for
 * the dynamic programming Solvers. A unit held on for 2 instants of 4 by
 * its minimum up time has 2 rows fixing its commitment, not 4. */

static void test_thermal_reactive_bounds( void )
{
 const auto sname = milp_solver();
 for( int given = 0 ; given < 2 ; ++given )
  for( int wf = 0 ; wf < 7 ; ++wf ) {
   std::vector< std::pair< std::string , std::vector< double > > > data = {
    { "MinPower" , { 10 , 10 , 10 } } , { "MaxPower" , { 100 , 100 , 100 } } ,
    { "LinearTerm" , { -10 , -10 , -10 } } };
   if( given ) {
    data.push_back( { "MinReactivePower" , { -20 , -20 , -20 } } );
    data.push_back( { "MaxReactivePower" , { 30 , 30 , 30 } } );
    }
   auto tub = new_TU_vec( 3 , data , 2 , 50 , 1 , 1 );
   tub->set_reactive_power( true );
   generate_all_wf( tub , wf );
   std::vector< double > rt( 3 , given ? 1 : -1 );
   tub->set_reactive_linear_term( rt.begin() , Range( 0 , 3 ) );
   const double value = given ? -3060 : -3000;
   const auto what = std::string( "reactive power, " ) +
                     ( given ? "bounds -20 and 30" : "no bounds" ) +
                     ", wf " + std::to_string( wf );
   if( ! sname.empty() ) {
    const auto v = milp_value( tub , sname );
    check( close( v , value ) , what + ": " + sname + " gives " + str( v ) +
           " instead of " + str( value ) );
    }
   check_all_DP_wf( tub , value , what , wf );
   delete tub;
   }

 auto tub = new_TU_vec( 4 , { { "MinPower" , { 10 , 10 , 10 , 10 } } ,
                              { "MaxPower" , { 100 , 100 , 100 , 100 } } } ,
                        1 , 50 , 3 , 1 );
 generate_all( tub );
 auto fx = tub->get_static_constraint_v< BoxConstraint >(
                                       "Commitment_fixed_to_one_Thermal" );
 check( fx && ( fx->size() == 2 ) ,
        "commitment fixed to one: " +
        std::to_string( fx ? fx->size() : 0 ) + " rows, not 2" );
 delete tub;
 }

/*--------------------------------------------------------------------------*/
/* The bounds of the node injections given to the NetworkBlock follow the
 * scale and the kappa of the units after the generation, as a fresh load
 * would compute them: in a UCBlock on two nodes, unit 0 at node 0 scaled
 * to 3 and the kappa of unit 1 at node 1 set to 2; and a NetworkBlock
 * whose rows of the bounds are generated (ECNetworkBlock) rewrites them. */

static void test_node_injection_bounds( void )
{
 auto g = write_UC_two_nodes( 1 , { 10 , 60 } );
 auto uc = dynamic_cast< UCBlock * >( Block::new_Block( g ) );
 check( uc , "node injection bounds: the UCBlock is not read" );
 if( uc ) {
  generate_all( uc );
  auto fs = new FakeSolver();
  uc->register_Solver( fs );
  auto nb = uc->get_network_block( 0 );
  // node 0: unit 0 (100); node 1: units 1 and 2 (100) and the slack (1000)
  check( close( nb->get_max_node_injection( 0 )[ 0 ] , 100 ) &&
         close( nb->get_max_node_injection( 0 )[ 1 ] , 1200 ) ,
         "node injection bounds: not 100 and 1200 at load" );
  uc->get_unit_block( 0 )->scale( 3 , eModBlck , eModBlck );
  check( close( nb->get_max_node_injection( 0 )[ 0 ] , 300 ) ,
         "node injection bounds: node 0 has " +
         str( nb->get_max_node_injection( 0 )[ 0 ] ) +
         " after the scale 3 of its unit, not 300" );
  static_cast< IntermittentUnitBlock * >( uc->get_unit_block( 1 ) )->
   set_kappa( 2 , eModBlck , eModBlck );
  check( close( nb->get_max_node_injection( 0 )[ 1 ] , 1300 ) ,
         "node injection bounds: node 1 has " +
         str( nb->get_max_node_injection( 0 )[ 1 ] ) +
         " after the kappa 2 of unit 1, not 1300" );
  uc->unregister_Solver( fs , true );
  delete uc;
  }

 auto e = new_group( "EC" , true );
 e.putAtt( "type" , "ECNetworkBlock" );
 auto I = e.addDim( "NumberIntervals" , 2 );
 auto N = e.addDim( "NumberNodes" , 2 );
 put( e , "BuyPrice" , I , { 3 , 4 } );
 put( e , "SellPrice" , I , { 1 , 1 } );
 put( e , "PeakTariff" , 2 );
 put( e , "ActiveDemand" , { I , N } , { 5 , 6 , 7 , 8 } );
 auto ec = dynamic_cast< NetworkBlock * >( Block::new_Block( e ) );
 check( ec , "node injection bounds, EC: not read" );
 if( ec ) {
  for( Index i = 0 ; i < 2 ; ++i )
   for( Index n = 0 ; n < 2 ; ++n ) {
    ec->set_min_node_injection( -10 , n , i );
    ec->set_max_node_injection( 10 , n , i );
    }
  generate_all( ec );
  ec->set_min_node_injection( -7 , 1 , 1 );
  ec->set_max_node_injection( 9 , 0 , 1 );
  auto b = ec->get_static_constraint< BoxConstraint , 2 >(
                                  "Node_Injection_Bound_Const_Network" );
  check( b && ( ( *b )[ 1 ][ 1 ].get_lhs() == -7 ) &&
         ( ( *b )[ 0 ][ 1 ].get_rhs() == 9 ) &&
         ( ( *b )[ 0 ][ 0 ].get_rhs() == 10 ) ,
         "node injection bounds, EC: the rows do not follow the bounds" );
  delete ec;
  }
 }

/*--------------------------------------------------------------------------*/
/* The rows of the inertia and of the reserves of a zone hold the Variable
 * node by node, so that those of a unit with generators at two nodes are
 * not consecutive: a HydroUnitBlock with arcs at nodes 0 and 1 (inertia 2
 * and 3) and an IntermittentUnitBlock at node 0 (inertia 1) in one zone,
 * the inertia of the arcs set to 5 and 7 after the generation, give 5, 7
 * and 1. A secondary zone over two nodes, its units scaled after the
 * generation, gives each unit its scale. */

static void test_zone_rows_two_nodes( void )
{
 auto g = new_group( "UC" , true );
 g.putAtt( "type" , "UCBlock" );
 auto TH = g.addDim( "TimeHorizon" , 1 );
 g.addDim( "NumberUnits" , 2 );
 auto N = g.addDim( "NumberNodes" , 2 );
 auto L = g.addDim( "NumberLines" , 1 );
 auto G = g.addDim( "NumberElectricalGenerators" , 3 );
 auto Z = g.addDim( "NumberInertiaZones" , 1 );
 put_int( g , "StartLine" , L , { 0 } );
 put_int( g , "EndLine" , L , { 1 } );
 put( g , "MinPowerFlow" , L , { -30 } );
 put( g , "MaxPowerFlow" , L , { 30 } );
 put( g , "ActivePowerDemand" , { N , TH } , { 5 , 5 } );
 put( g , "InertiaDemand" , { Z , TH } , { 40 } );
 put_int( g , "GeneratorNode" , G , { 0 , 1 , 0 } );
 HUData d;
 d.A = 2;
 d.start = { 0 , 0 };
 d.end = { 1 , 1 };
 d.maxF = { 10 , 10 };
 d.maxP = { 10 , 10 };
 d.lin = { 1 , 1 };
 d.inertia = { 2 , 3 };
 d.initV = { 100 };
 d.maxV = { 1000 };
 d.inflows = { 0 };
 write_HU( g.addGroup( "UnitBlock_0" ) , d );
 auto u = g.addGroup( "UnitBlock_1" );
 u.putAtt( "type" , "IntermittentUnitBlock" );
 put( u , "MaxPower" , 100.0 );
 put( u , "ActivePowerCost" , 50.0 );
 put( u , "InertiaPower" , 1.0 );

 auto uc = dynamic_cast< UCBlock * >( Block::new_Block( g ) );
 check( uc , "zone rows, two nodes: the UCBlock is not read" );
 if( uc ) {
  generate_all( uc );
  auto fs = new FakeSolver();
  uc->register_Solver( fs );
  auto hub = static_cast< HydroUnitBlock * >( uc->get_unit_block( 0 ) );
  std::vector< double > val = { 5 , 7 };
  hub->set_inertia_power( val.begin() , Range( 0 , 2 ) );
  auto ub = static_cast< UnitBlock * >( hub );
  const auto c0 = inertia_coef( uc , ub->get_active_power( 0 ) );
  const auto c1 = inertia_coef( uc , ub->get_active_power( 1 ) );
  const auto c2 = inertia_coef( uc ,
                          uc->get_unit_block( 1 )->get_active_power( 0 ) );
  check( close( c0 , 5 ) && close( c1 , 7 ) && close( c2 , 1 ) ,
         "zone rows, two nodes: inertia coefficients " + str( c0 ) + ", " +
         str( c1 ) + " and " + str( c2 ) + " instead of 5, 7 and 1" );
  uc->unregister_Solver( fs , true );
  delete uc;
  }

 auto s = write_UC_two_nodes( 1 , { 10 , 60 } );
 auto SZ = s.addDim( "NumberSecondaryZones" , 1 );
 put( s , "SecondaryDemand" , { SZ , s.getDim( "TimeHorizon" ) } , { 5 } );
 auto us = dynamic_cast< UCBlock * >( Block::new_Block( s ) );
 check( us , "zone rows, secondary: the UCBlock is not read" );
 if( us ) {
  generate_all( us );
  auto fs = new FakeSolver();
  us->register_Solver( fs );
  const auto & row = us->get_const_secondary_demand_constraints()[ 0 ][ 0 ];
  us->get_unit_block( 1 )->scale( 0.5 , eModBlck , eModBlck );
  us->get_unit_block( 0 )->scale( 3 , eModBlck , eModBlck );
  check( close( coef_of( us , row , 0 ) , 3 ) &&
         close( coef_of( us , row , 1 ) , 0.5 ) ,
         "zone rows, secondary: coefficients " +
         str( coef_of( us , row , 0 ) ) + " and " +
         str( coef_of( us , row , 1 ) ) + " instead of 3 and 0.5" );
  us->unregister_Solver( fs , true );
  delete us;
  }
 }

/*--------------------------------------------------------------------------*/
/* BatteryUnitBlock: with one instant and the cyclic balance the storage
 * level is once in the row of instant 0, with coefficient 1 minus the
 * standing loss; an infinite MinStorage stays infinite whatever kappa is
 * (also 0), and with a design variable it gives that variable no (infinite
 * or NaN) coefficient, the row being free on that side; the :MILPSolver
 * gives -50 and -497 to the first and the last. SlackUnitBlock with the
 * reactive power and no bound of it has the rows of its absolute value,
 * the reactive power being 0. */

static void test_storage_infinite_bounds( void )
{
 const auto sname = milp_solver();
 const std::vector< std::pair< std::string , std::vector< double > > > base =
  { { "MaxStorage" , { 100 } } , { "MaxPower" , { 50 } } ,
    { "MinPower" , { -50 } } , { "Cost" , { -5 } } };

 auto data = base;
 data.push_back( { "MinStorage" , { 0 } } );
 data.push_back( { "InitialStorage" , { -1 } } );
 data.push_back( { "StandingBatteryRho" , { 0.9 } } );
 if( auto bu = new_BU( 1 , data ) ) {
  generate_all( bu );
  const auto & bal = bu->get_storage_balance_constraints();
  const auto sl = & bu->get_storage_level()[ 0 ];
  Index n = 0;
  for( Index k = 0 ;
       ( ! bal.empty() ) && ( k < bal[ 0 ].get_num_active_var() ) ; ++k )
   if( bal[ 0 ].get_active_var( k ) == sl )
    ++n;
  check( ( n == 1 ) && close( row_coef( bal[ 0 ] , sl ) , 0.1 ) ,
         "battery, one instant, cyclic: the level is " + std::to_string( n ) +
         " times in the row, coefficient " +
         str( row_coef( bal[ 0 ] , sl ) ) );
  if( ! sname.empty() ) {
   const auto v = milp_value( bu , sname );
   check( close( v , -50 ) , "battery, one instant, cyclic: " + sname +
          " gives " + str( v ) + " instead of -50" );
   }
  delete bu;
  }
 else
  check( false , "battery, one instant, cyclic: not read" );

 for( int design = 0 ; design < 2 ; ++design ) {
  auto data2 = base;
  data2.push_back( { "MinStorage" , { -INF , -INF } } );
  data2.push_back( { "InitialStorage" , { 0 } } );
  data2.push_back( { "Kappa" , { design ? 1.0 : 0.0 } } );
  if( design )
   data2.push_back( { "BatteryInvestmentCost" , { 3 } } );
  for( auto & dd : data2 )  // the per-instant data over 2 instants
   if( ( dd.second.size() == 1 ) && ( dd.first != "InitialStorage" ) &&
       ( dd.first != "Kappa" ) && ( dd.first != "BatteryInvestmentCost" ) )
    dd.second.assign( 2 , dd.second[ 0 ] );
  const auto what = std::string( "battery, MinStorage -inf, " ) +
                    ( design ? "design" : "kappa 0" );
  auto bu = new_BU( 2 , data2 );
  if( ! bu ) {
   check( false , what + ": not read" );
   continue;
   }
  generate_all( bu );
  if( design ) {
   auto rows = bu->get_static_constraint< FRowConstraint , 2 >(
                                         "StorageLevel_Design_Battery" );
   bool ok = rows && ( rows->num_elements() > 0 );
   if( ok )
    for( auto it = rows->data() ; it != rows->data() + rows->num_elements() ;
         ++it ) {
     auto lf = static_cast< const LinearFunction * >( it->get_function() );
     for( const auto & vp : lf->get_v_var() )
      ok = ok && std::isfinite( vp.second );
     }
   ok = ok && ( ( *rows )[ 0 ][ 0 ].get_lhs() == -INF );
   check( ok , what + ": an infinite or NaN coefficient, or a bound row" );
   if( ! sname.empty() ) {
    const auto v = milp_value( bu , sname );
    check( close( v , -497 ) , what + ": " + sname + " gives " + str( v ) +
           " instead of -497" );
    }
   }
  else {
   const auto & sb = bu->get_storage_level_bounds();
   check( ( sb.size() == 2 ) && ( sb[ 0 ].get_lhs() == -INF ) ,
          what + ": the lower bound of the level is " +
          str( sb.empty() ? 0 : sb[ 0 ].get_lhs() ) );
   }
  delete bu;
  }

 auto g = new_group( "SU" , true );
 g.putAtt( "type" , "SlackUnitBlock" );
 auto T = g.addDim( "TimeHorizon" , 2 );
 put( g , "MaxPower" , T , { 50 , 50 } );
 put( g , "ActivePowerCost" , T , { 10 , 10 } );
 auto sb = dynamic_cast< SlackUnitBlock * >( Block::new_Block( g ) );
 if( sb ) {
  sb->set_reactive_power( true );
  generate_all( sb );
  auto ab = sb->get_static_constraint_v< FRowConstraint >(
                                                   "Lin_of_Abs_Reactive" );
  auto qb = sb->get_static_constraint_v< BoxConstraint >(
                                            "ReactivePowerBound_thermal" );
  check( ab && ( ab->size() == 4 ) && qb && ( qb->size() == 2 ) &&
         ( ( *qb )[ 0 ].get_lhs() == 0 ) && ( ( *qb )[ 0 ].get_rhs() == 0 ) ,
         "slack, reactive power with no bound: the rows are not there, or "
         "do not hold it at 0" );
  delete sb;
  }
 else
  check( false , "slack, reactive power with no bound: not read" );
 }

/*--------------------------------------------------------------------------*/
/* The flows of an ACNetworkBlock are one group of static Variable, the one
 * that holds all of them, and no other group holds them. */

static void test_ac_flow_group( void )
{
 auto g = new_group( "AC" , true );
 g.putAtt( "type" , "ACNetworkBlock" );
 auto N = g.addDim( "NumberNodes" , 3 );
 auto L = g.addDim( "NumberLines" , 2 );
 put_int( g , "StartLine" , L , { 0 , 1 } );
 put_int( g , "EndLine" , L , { 1 , 2 } );
 put( g , "MinPowerFlow" , L , { -5 , -6 } );
 put( g , "MaxPowerFlow" , L , { 5 , 6 } );
 put( g , "LineSusceptance" , L , { 0 , 2 } );
 put( g , "ActiveDemand" , N , { 1 , 2 , 3 } );
 put( g , "LineReactance" , L , { 0.1 , 0.2 } );
 put( g , "LineResistance" , L , { 0.01 , 0.02 } );
 put( g , "NodeMaxVoltage" , N , { 1.1 , 1.1 , 1.1 } );
 put( g , "NodeMinVoltage" , N , { 0.9 , 0.9 , 0.9 } );
 put( g , "ReactiveDemand" , N , { 0.5 , 0.5 , 0.5 } );
 auto ac = dynamic_cast< ACNetworkBlock * >( Block::new_Block( g ) );
 check( ac , "AC flows: not read" );
 if( ! ac )
  return;
 ac->generate_abstract_variables();
 auto f = ac->get_static_variable_v< ColVariable >( "v_power_flow_real" );
 Index holders = 0;
 for( Index i = 0 ; i < ac->get_number_static_variables() ; ++i )
  if( f && ( ac->get_static_variable_v< ColVariable >( i ) == f ) )
   ++holders;
 check( f && ( f->size() == 4 ) && ( holders == 1 ) ,
        "AC flows: " + std::to_string( holders ) + " groups hold the " +
        std::to_string( f ? f->size() : 0 ) + " flows" );
 delete ac;
 }

/*--------------------------------------------------------------------------*/
/*------------------- THERMAL: EDGE CASES OF THE ROWS ----------------------*/
/*--------------------------------------------------------------------------*/
/* A ThermalUnitBlock or a NuclearUnitBlock written from per-instant vectors
 * and named scalars, so that every datum can change with the instant or be
 * left out; the vectors are written over "NumberIntervals" = T, the bands
 * (if any) over "NumberPowerBands". */

static ThermalUnitBlock * new_unit_gen(
 bool nuclear , Index T ,
 const std::vector< std::pair< std::string , std::vector< double > > > & vecs ,
 const std::vector< std::pair< std::string , double > > & dbls ,
 const std::vector< std::pair< std::string , int > > & ints ,
 const std::vector< std::pair< std::string , unsigned int > > & uints ,
 const std::vector< double > & bands = { } )
{
 auto g = new_group( "TUg" , true );
 g.putAtt( "type" , nuclear ? "NuclearUnitBlock" : "ThermalUnitBlock" );
 g.addDim( "TimeHorizon" , T );
 auto NI = g.addDim( "NumberIntervals" , T );
 for( const auto & nv : vecs )
  put( g , nv.first , NI , nv.second );
 for( const auto & nv : dbls )
  put( g , nv.first , nv.second );
 for( const auto & nv : ints )
  put_int( g , nv.first , nv.second );
 for( const auto & nv : uints )
  put_uint( g , nv.first , nv.second );
 if( ! bands.empty() ) {
  auto nb = g.addDim( "NumberPowerBands" , bands.size() );
  put( g , "PowerBands" , nb , bands );
  }
 auto tub = dynamic_cast< ThermalUnitBlock * >( Block::new_Block( g ) );
 if( ! tub )
  throw( std::logic_error( "new_unit_gen: no ThermalUnitBlock built" ) );
 return( tub );
 }

/*--------------------------------------------------------------------------*/
/// the value of the unit in every formulation (MILP, if a :MILPSolver is
/// there) and with every DP Solver, against the expected one

static void check_all_forms( const std::function< ThermalUnitBlock *( void )
                             > & make , double expected ,
                             const std::string & what ,
                             std::vector< int > forms = { 0 , 1 , 2 , 3 ,
                                                          4 , 5 , 6 } )
{
 const auto sname = milp_solver();
 for( int wf : forms ) {
  auto tub = make();
  generate_all_wf( tub , wf );
  const auto who = what + ", wf " + std::to_string( wf );
  if( ! sname.empty() ) {
   const auto v = milp_value( tub , sname );
   check( close( v , expected ) , who + ": " + sname + " gives " + str( v ) +
          " instead of " + str( expected ) );
   }
  check_all_DP_wf( tub , expected , who , wf );
  delete tub;
  }
 }

/*--------------------------------------------------------------------------*/
/* The minimum up and down times over a horizon shorter than them: on since
 * 5 instants with the minimum up time 5 over T = 4, the unit may shut down
 * at 0 and start again at 1 (minimum down time 1), but then it has to stay
 * on to the end of the horizon, a run of 2 instants that ends inside it
 * being shorter than the minimum up time even if no window of 5 instants
 * fits in the horizon (the optimum -190, the run of 2 giving -200); the
 * same for the minimum down time with a unit off since 5 instants (-180,
 * the restart after 2 instants giving -200). Every formulation and both DP
 * Solvers, against the brute force. */

static void test_thermal_short_horizon_min_up_down( void )
{
 TUData d;
 d.T = 4;
 d.minP = 1;
 d.maxP = d.ru = d.rd = d.su = d.sd = 10;
 d.initUD = 5;
 d.initP = 5;
 d.minUp = 5;
 d.minDown = 1;
 d.lin = { 10 , -10 , -10 , 10 };
 const std::vector< int > free( d.T , -1 );
 check( close( brute_force( d , free ) , -190 ) ,
        "short horizon, minimum up time: the brute force gives " +
        str( brute_force( d , free ) ) + " instead of -190" );
 check_all_forms( [ & ]() { return( new_TU( d ) ); } , -190 ,
                  "short horizon, minimum up time 5 over T = 4" );

 d.initUD = -5;
 d.minUp = 1;
 d.minDown = 5;
 d.lin = { -10 , 10 , 10 , -10 };
 check( close( brute_force( d , free ) , -180 ) ,
        "short horizon, minimum down time: the brute force gives " +
        str( brute_force( d , free ) ) + " instead of -180" );
 check_all_forms( [ & ]() { return( new_TU( d ) ); } , -180 ,
                  "short horizon, minimum down time 5 over T = 4" );
 }

/*--------------------------------------------------------------------------*/
/* The multi-step ramp rows of the SUSD formulation with bounds that change
 * in time: an outage at 2 (Availability 0) of a unit on at 100 with the
 * ramp-up 20, which then restarts from 0 and reaches at most 20 at 3
 * (-200; the rows that start at the outage instant are not implied by the
 * bounds there); and a minimum power that drops from 90 to 10 at 2 with the
 * ramp-down 20, which keeps the output at least 70 at 2 and 50 at 3 (500,
 * a shut-down being too expensive).
 * Every formulation and both DP Solvers. */

static void test_thermal_SUSD_varying_bounds( void )
{
 check_all_forms( []() {
   return( new_unit_gen( false , 4 ,
           { { "MinPower" , { 0 , 0 , 0 , 0 } } ,
             { "MaxPower" , { 100 , 100 , 100 , 100 } } ,
             { "Availability" , { 1 , 1 , 0 , 1 } } ,
             { "DeltaRampUp" , { 20 , 20 , 20 , 20 } } ,
             { "DeltaRampDown" , { 100 , 100 , 100 , 100 } } ,
             { "StartUpLimit" , { 100 , 100 , 0 , 100 } } ,
             { "ShutDownLimit" , { 100 , 100 , 0 , 100 } } ,
             { "LinearTerm" , { 0 , 0 , 0 , -10 } } ,
             { "StartUpCost" , { 0 , 0 , 0 , 10000 } } } ,
           { { "InitialPower" , 100 } } , { { "InitUpDownTime" , 2 } } ,
           { { "MinUpTime" , 1 } , { "MinDownTime" , 1 } } ) ); } ,
   -200 , "SUSD ramps, outage at 2" );

 check_all_forms( []() {
   return( new_unit_gen( false , 4 ,
           { { "MinPower" , { 90 , 90 , 10 , 10 } } ,
             { "MaxPower" , { 100 , 100 , 100 , 100 } } ,
             { "DeltaRampUp" , { 100 , 100 , 100 , 100 } } ,
             { "DeltaRampDown" , { 20 , 20 , 20 , 20 } } ,
             { "LinearTerm" , { 0 , 0 , 0 , 10 } } ,
             { "ShutDownCost" , { 1e4 , 1e4 , 1e4 , 1e4 } } } ,
           { { "InitialPower" , 100 } } , { { "InitUpDownTime" , 2 } } ,
           { { "MinUpTime" , 1 } , { "MinDownTime" , 1 } } ) ); } ,
   500 , "SUSD ramps, minimum power from 90 to 10" );
 }

/*--------------------------------------------------------------------------*/
/* The ramps of the bound rows of the T formulation equal to 0: a unit that
 * cannot move from the start-up limit 30 (DeltaRampUp 0) is generated with
 * no division by 0 and its value is the one of every other formulation. */

static void test_thermal_zero_ramp( void )
{
 check_all_forms( []() {
   return( new_unit_gen( false , 3 ,
           { { "MinPower" , { 10 , 10 , 10 } } ,
             { "MaxPower" , { 100 , 100 , 100 } } ,
             { "DeltaRampUp" , { 0 , 0 , 0 } } ,
             { "DeltaRampDown" , { 0 , 0 , 0 } } ,
             { "StartUpLimit" , { 30 , 30 , 30 } } ,
             { "ShutDownLimit" , { 30 , 30 , 30 } } ,
             { "LinearTerm" , { -1 , -1 , -1 } } } ,
           { } , { { "InitUpDownTime" , -1 } } ,
           { { "MinUpTime" , 1 } , { "MinDownTime" , 1 } } ) ); } ,
   -90 , "zero ramps" );
 }

/*--------------------------------------------------------------------------*/
/* InitialPower changed after the generation without ramp data, across
 * ShutDownLimit[ 0 ] (50): the unit may shut down at 0 iff InitialPower <=
 * ShutDownLimit[ 0 ], in every formulation, and the change is accepted both
 * after the Variable alone and after the Constraint are generated; with
 * the rows, the dynamic programming Solver (which read the new value)
 * give the value of a fresh load with a schedule that satisfies them: 0
 * with InitialPower 45, the unit shutting down at 0, and 10 with 60, the
 * unit on at 0 at MinPower and off from 1 on. */

static void test_thermal_initial_power_no_ramps( void )
{
 auto make = []() {
  return( new_unit_gen( false , 3 ,
          { { "MinPower" , { 10 , 10 , 10 } } ,
            { "MaxPower" , { 100 , 100 , 100 } } ,
            { "ShutDownLimit" , { 50 , 50 , 50 } } ,
            { "LinearTerm" , { 1 , 1 , 1 } } } ,
          { { "InitialPower" , 40 } } , { { "InitUpDownTime" , 2 } } ,
          { { "MinUpTime" , 1 } , { "MinDownTime" , 1 } } ) );
  };
 for( int wf = 0 ; wf < 7 ; ++wf )
  for( int what = 0 ; what < 2 ; ++what ) {
   auto tub = make();
   SimpleConfiguration< int > f( wf );
   tub->generate_abstract_variables( & f );
   if( what ) {
    tub->generate_abstract_constraints();
    tub->generate_objective();
    }
   const auto who = std::string( "InitialPower, no ramps, wf " ) +
                    std::to_string( wf ) + ( what ? ", all" : ", Variable" );
   std::vector< double > ip = { 45 };
   bool thrown = false;
   try {
    tub->set_initial_power( ip.begin() , Range( 0 , 1 ) );
    }
   catch( std::logic_error & ) {
    thrown = true;
    }
   check( ( ! thrown ) && ( tub->get_initial_power() == 45 ) ,
          who + ": 40 -> 45 refused" );
   if( what )
    check_all_DP_wf( tub , 0 , who + ", 45" , wf );
   ip = { 60 };
   thrown = false;
   try {
    tub->set_initial_power( ip.begin() , Range( 0 , 1 ) );
    }
   catch( std::logic_error & ) {
    thrown = true;
    }
   check( ( ! thrown ) && ( tub->get_initial_power() == 60 ) ,
          who + ": 45 -> 60 refused" );
   if( what )
    check_all_DP_wf( tub , 10 , who + ", 60" , wf );
   delete tub;
   }

 // a NuclearUnitBlock with its Variable alone: the change is accepted
 auto nub = new_unit_gen( true , 3 ,
                          { { "MinPower" , { 10 , 10 , 10 } } ,
                            { "MaxPower" , { 100 , 100 , 100 } } ,
                            { "DeltaRampUp" , { 20 , 20 , 20 } } ,
                            { "DeltaRampDown" , { 20 , 20 , 20 } } ,
                            { "ModulationDeltaRampUp" , { 20 , 20 , 20 } } ,
                            { "ModulationDeltaRampDown" , { 20 , 20 , 20 } } ,
                            { "LinearTerm" , { 1 , 1 , 1 } } } ,
                          { { "InitialPower" , 40 } } ,
                          { { "InitUpDownTime" , 2 } } ,
                          { { "MinUpTime" , 1 } , { "MinDownTime" , 1 } ,
                            { "ModulationTime" , 2 } ,
                            { "InitModulation" , 2 } } );
 nub->generate_abstract_variables();
 std::vector< double > ip = { 45 };
 bool thrown = false;
 try {
  nub->set_initial_power( ip.begin() , Range( 0 , 1 ) );
  }
 catch( std::exception & ) {
  thrown = true;
  }
 check( ( ! thrown ) && ( nub->get_initial_power() == 45 ) ,
        "InitialPower of a NuclearUnitBlock with the Variable alone: "
        "refused" );
 delete nub;
 }

/*--------------------------------------------------------------------------*/
/* The shut-down at 0 of a unit on before the horizon, whose minimum up
 * time is over (InitUpDownTime 2, MinUpTime 1), with the default
 * ShutDownLimit (MinPower, 10) and a positive price, so that the unit wants
 * to be off: the output before a shut-down is at most the shut-down limit,
 * hence with InitialPower 80 the unit cannot shut down at 0 and it is on at
 * 0 at its minimum power, off from 1 on (100), without ramp data ("sd0")
 * and with a ramp-down that does not bind ("sd0r"); with InitialPower 10
 * it shuts down at 0 (0, "sd0e"); with a ramp-down of 30 it goes down
 * 80 -> 50 -> 20 -> 10 and shuts down at 3 (800, "sd0d"). The same value in
 * every formulation (with a :MILPSolver, if there is one) and with every
 * dynamic programming Solver. */

static void test_thermal_shut_down_at_zero( void )
{
 struct Case {
  const char * name;
  double ip;            // InitialPower
  double rd;            // DeltaRampDown, if > 0
  double value;
  };
 for( const auto & c : { Case{ "sd0" , 80 , 0 , 100 } ,
                         Case{ "sd0r" , 80 , 200 , 100 } ,
                         Case{ "sd0e" , 10 , 0 , 0 } ,
                         Case{ "sd0d" , 80 , 30 , 800 } } )
  check_all_forms( [ & ]() {
   std::vector< std::pair< std::string , std::vector< double > > > v = {
    { "MinPower" , { 10 , 10 , 10 , 10 } } ,
    { "MaxPower" , { 100 , 100 , 100 , 100 } } ,
    { "LinearTerm" , { 10 , 10 , 10 , 10 } } };
   if( c.rd > 0 )
    v.push_back( { "DeltaRampDown" , std::vector< double >( 4 , c.rd ) } );
   return( new_unit_gen( false , 4 , v , { { "InitialPower" , c.ip } } ,
                         { { "InitUpDownTime" , 2 } } ,
                         { { "MinUpTime" , 1 } , { "MinDownTime" , 1 } } ) );
   } , c.value , std::string( "shut-down at 0, " ) + c.name );
 }

/*--------------------------------------------------------------------------*/
/* A unit on before the horizon whose InitialPower (80) is above MaxPower at
 * every instant, without ramp data and with the default limits: no row
 * limits the move from InitialPower, the unit cannot shut down at 0 (80 is
 * above ShutDownLimit[ 0 ]), so it is on at 0 at its minimum power and off
 * from 1 on. Unavailable over the whole horizon (MinPower = MaxPower = 0,
 * "ip>mx0"): value 0; MaxPower 50 and MinPower 10 ("ip>mx"): value 100.
 * The same value in every formulation and with every dynamic programming
 * Solver, which must not take the largest MaxPower as the ramp the move
 * from InitialPower is limited by. */

static void test_thermal_initial_power_above_max( void )
{
 struct Case {
  const char * name;
  double mn;            // MinPower
  double mx;            // MaxPower
  double value;
  };
 for( const auto & c : { Case{ "ip>mx0" , 0 , 0 , 0 } ,
                         Case{ "ip>mx" , 10 , 50 , 100 } } )
  check_all_forms( [ & ]() {
   return( new_unit_gen( false , 4 ,
                         { { "MinPower" , std::vector< double >( 4 , c.mn ) } ,
                           { "MaxPower" , std::vector< double >( 4 , c.mx ) } ,
                           { "LinearTerm" , { 10 , 10 , 10 , 10 } } } ,
                         { { "InitialPower" , 80 } } ,
                         { { "InitUpDownTime" , 2 } } ,
                         { { "MinUpTime" , 1 } , { "MinDownTime" , 1 } } ) );
   } , c.value , std::string( "InitialPower above MaxPower, " ) + c.name );
 }

/*--------------------------------------------------------------------------*/
/* The unit of test_thermal_initial_power_above_max() ("ip>mx"), without
 * ramp data, with InitialPower 40 changed to 80 by set_initial_power()
 * between two compute() of the same dynamic programming Solver. Both are
 * above ShutDownLimit[ 0 ] (MinPower 10), so the unit is on at 0 at its
 * minimum power and off from 1 on, value 100, the move from InitialPower
 * being free; a Solver that keeps the ramp 50 it took for InitialPower 40
 * needs p_0 >= 30 at 80, and then p_1 = 10 to shut down, value 400. */

static void test_thermal_initial_power_above_max_change( void )
{
 for( const auto & s : THERMAL_DP ) {
  const auto what = "InitialPower 40 -> 80 above MaxPower, " + s;
  auto tub = new_unit_gen( false , 4 ,
                           { { "MinPower" , std::vector< double >( 4 , 10 ) } ,
                             { "MaxPower" , std::vector< double >( 4 , 50 ) } ,
                             { "LinearTerm" , { 10 , 10 , 10 , 10 } } } ,
                           { { "InitialPower" , 40 } } ,
                           { { "InitUpDownTime" , 2 } } ,
                           { { "MinUpTime" , 1 } , { "MinDownTime" , 1 } } );
  auto slv = Solver::new_Solver( s );
  tub->register_Solver( slv );
  try {
   auto status = slv->compute();
   check( ( status == Solver::kOK ) && close( slv->get_ub() , 100 ) ,
          what + ", at 40: status " + std::to_string( status ) +
          ", value " + str( slv->get_ub() ) + " instead of 100" );
   std::vector< double > val = { 80 };
   tub->set_initial_power( val.begin() , Range( 0 , 1 ) );
   status = slv->compute();
   check( ( status == Solver::kOK ) && close( slv->get_ub() , 100 ) ,
          what + ", at 80: status " + std::to_string( status ) +
          ", value " + str( slv->get_ub() ) + " instead of 100" );
   }
  catch( std::exception & e ) {
   check( false , what + ": throws " + e.what() );
   }
  tub->unregister_Solver( slv , true );
  delete tub;
  }
 }

/*--------------------------------------------------------------------------*/
/* The setters of the data of a ThermalUnitBlock with a Subset given out of
 * order (ordered = false): each value goes with the instant it is given
 * with, e.g., MaxPower 60 at 4 and 140 at 1 for the Subset { 4 , 1 } and
 * the values { 60 , 140 }, for the ten setters that take a Subset. */

static void test_thermal_unordered_subset( void )
{
 const Index T = 6;
 auto cst = [ & ]( double v ) { return( std::vector< double >( T , v ) ); };
 auto tub = new_unit_gen( false , T ,
                          { { "MinPower" , cst( 10 ) } ,
                            { "MaxPower" , cst( 100 ) } ,
                            { "Availability" , cst( 1 ) } ,
                            { "LinearTerm" , cst( 1 ) } ,
                            { "QuadTerm" , cst( 0.1 ) } ,
                            { "ConstTerm" , cst( 5 ) } ,
                            { "StartUpCost" , cst( 7 ) } ,
                            { "ShutDownCost" , cst( 3 ) } ,
                            { "ReactiveLinearTerm" , cst( 2 ) } ,
                            { "PrimaryRho" , cst( 0.1 ) } ,
                            { "SecondaryRho" , cst( 0.1 ) } ,
                            { "PrimarySpinningReserveCost" , cst( 1 ) } ,
                            { "SecondarySpinningReserveCost" ,
                              cst( 1 ) } } ,
                          { } , { { "InitUpDownTime" , -1 } } ,
                          { { "MinUpTime" , 1 } , { "MinDownTime" , 1 } } );

 using Setter = std::function< void( Block::MF_dbl_it , Subset && ) >;
 using Getter = std::function< const std::vector< double > & ( void ) >;
 struct S {
  const char * name;
  Setter set;
  Getter get;
  std::vector< double > val;  // the values for the instants 4 and 1
  };
 std::vector< S > setters = {
  { "MaxPower" ,
    [ & ]( Block::MF_dbl_it v , Subset && s ) {
     tub->set_maximum_power( v , std::move( s ) ); } ,
    [ & ]() -> const std::vector< double > & {
     return( tub->get_max_power() ); } , { 60 , 140 } } ,
  { "Availability" ,
    [ & ]( Block::MF_dbl_it v , Subset && s ) {
     tub->set_availability( v , std::move( s ) ); } ,
    [ & ]() -> const std::vector< double > & {
     return( tub->get_availability() ); } , { 0.5 , 0.8 } } ,
  { "StartUpCost" ,
    [ & ]( Block::MF_dbl_it v , Subset && s ) {
     tub->set_startup_costs( v , std::move( s ) ); } ,
    [ & ]() -> const std::vector< double > & {
     return( tub->get_start_up_cost() ); } , { 11 , 13 } } ,
  { "ShutDownCost" ,
    [ & ]( Block::MF_dbl_it v , Subset && s ) {
     tub->set_shutdown_costs( v , std::move( s ) ); } ,
    [ & ]() -> const std::vector< double > & {
     return( tub->get_shut_down_cost() ); } , { 4 , 6 } } ,
  { "ConstTerm" ,
    [ & ]( Block::MF_dbl_it v , Subset && s ) {
     tub->set_const_term( v , std::move( s ) ); } ,
    [ & ]() -> const std::vector< double > & {
     return( tub->get_const_term() ); } , { 8 , 9 } } ,
  { "LinearTerm" ,
    [ & ]( Block::MF_dbl_it v , Subset && s ) {
     tub->set_linear_term( v , std::move( s ) ); } ,
    [ & ]() -> const std::vector< double > & {
     return( tub->get_linear_term() ); } , { -2 , 3 } } ,
  { "ReactiveLinearTerm" ,
    [ & ]( Block::MF_dbl_it v , Subset && s ) {
     tub->set_reactive_linear_term( v , std::move( s ) ); } ,
    [ & ]() -> const std::vector< double > & {
     return( tub->get_reactive_linear_term() ); } , { -4 , 5 } } ,
  { "QuadTerm" ,
    [ & ]( Block::MF_dbl_it v , Subset && s ) {
     tub->set_quad_term( v , std::move( s ) ); } ,
    [ & ]() -> const std::vector< double > & {
     return( tub->get_quad_term() ); } , { 0.2 , 0.3 } } ,
  { "PrimarySpinningReserveCost" ,
    [ & ]( Block::MF_dbl_it v , Subset && s ) {
     tub->set_primary_spinning_reserve_cost( v , std::move( s ) ); } ,
    [ & ]() -> const std::vector< double > & {
     return( tub->get_primary_spinning_reserve_cost() ); } , { 2 , 4 } } ,
  { "SecondarySpinningReserveCost" ,
    [ & ]( Block::MF_dbl_it v , Subset && s ) {
     tub->set_secondary_spinning_reserve_cost( v , std::move( s ) ); } ,
    [ & ]() -> const std::vector< double > & {
     return( tub->get_secondary_spinning_reserve_cost() ); } ,
    { 6 , 8 } } };

 for( auto & st : setters ) {
  st.set( st.val.cbegin() , Subset( { 4 , 1 } ) );
  const auto & v = st.get();
  check( ( v.size() == T ) && ( v[ 4 ] == st.val[ 0 ] ) &&
         ( v[ 1 ] == st.val[ 1 ] ) ,
         std::string( "unordered Subset, " ) + st.name + ": [ 4 ] = " +
         ( v.size() == T ? str( v[ 4 ] ) : std::string( "?" ) ) +
         ", [ 1 ] = " + ( v.size() == T ? str( v[ 1 ] ) :
                                           std::string( "?" ) ) +
         " instead of " + str( st.val[ 0 ] ) + ", " + str( st.val[ 1 ] ) );
  }
 delete tub;
 }

/*--------------------------------------------------------------------------*/
/* The costs of the reserves and of the shut-down set before the Variable
 * are generated: a reserve cost set before set_reserve_vars() and before
 * the generation is in the Objective; after the generation, a unit with no
 * reserve Variable ignores it; a nonzero shut-down cost of a unit generated
 * without the shut-down term is refused and leaves the data as they are. */

static void test_thermal_cost_order( void )
{
 auto make = []() {
  return( new_unit_gen( false , 2 ,
          { { "MinPower" , { 10 , 10 } } , { "MaxPower" , { 100 , 100 } } ,
            { "LinearTerm" , { 1 , 1 } } ,
            { "PrimaryRho" , { 0.1 , 0.1 } } ,
            { "SecondaryRho" , { 0.1 , 0.1 } } } ,
          { { "InitialPower" , 50 } } , { { "InitUpDownTime" , 2 } } ,
          { { "MinUpTime" , 1 } , { "MinDownTime" , 1 } } ) );
  };

 auto tub = make();
 std::vector< double > c = { -3 , -4 };
 tub->set_primary_spinning_reserve_cost( c.begin() , Range( 0 , 2 ) );
 tub->set_secondary_spinning_reserve_cost( c.begin() , Range( 0 , 2 ) );
 tub->set_reserve_vars( 3 );
 generate_all( tub );
 auto pr = tub->get_primary_spinning_reserve( 0 );
 auto sc = tub->get_secondary_spinning_reserve( 0 );
 check( pr && sc && in_obj( tub , pr + 1 ) && in_obj( tub , sc + 1 ) &&
        close( lin_coef( tub , pr + 1 ) , -4 ) &&
        close( lin_coef( tub , sc + 1 ) , -4 ) ,
        "reserve costs set before the generation: not in the Objective" );
 delete tub;

 tub = make();  // no reserve Variable: the cost is ignored, no throw
 generate_all( tub );
 bool thrown = false;
 try {
  tub->set_primary_spinning_reserve_cost( c.begin() , Range( 0 , 2 ) );
  }
 catch( std::exception & ) {
  thrown = true;
  }
 check( ( ! thrown ) && tub->get_primary_spinning_reserve_cost().empty() ,
        "reserve cost of a unit with no reserve Variable: refused, or kept" );

 thrown = false;  // no ShutDownCost: a nonzero one refused, nothing kept
 try {
  tub->set_shutdown_costs( c.begin() , Range( 0 , 2 ) );
  }
 catch( std::logic_error & ) {
  thrown = true;
  }
 check( thrown && tub->get_shut_down_cost().empty() ,
        "shut-down cost of a unit with no shut-down term: not refused, or "
        "the data changed" );
 std::vector< double > z = { 0 , 0 };
 thrown = false;
 try {
  tub->set_shutdown_costs( z.begin() , Range( 0 , 2 ) );
  }
 catch( std::exception & ) {
  thrown = true;
  }
 check( ! thrown , "zero shut-down cost of a unit with no shut-down term: "
        "refused" );
 delete tub;
 }

/*--------------------------------------------------------------------------*/
/* What the DP Solvers take from the data must be what the rows have: a unit
 * with priced reserves (PrimaryRho, SecondaryRho and negative costs) whose
 * reserve Variable are not there, the enclosing UCBlock not asking for
 * them, has no reserve reward (the value of the brute force, which has no
 * reserve, i.e., 0 for a unit with no energy cost whose reserve would be
 * worth -1 per MW); a unit with FixToMaximum produces its operational
 * maximum power, being on also at an outage instant at p = 0 (154 instead
 * of 0); a NuclearUnitBlock with PowerBands that cannot shut down at an
 * outage is on there at p = 0, in the lowest band, as for the DP Solver
 * (20, the rows with the nominal MinPower being infeasible); a deep decrease
 * fixed to 1 is refused by NuclearUnitExtDPSolver, and one fixed to 0 is
 * not. */

static void test_thermal_DP_as_rows( void )
{
 // priced reserves with no Variable - - - - - - - - - - - - - - - - - - - -
 {
  auto d = reserve_TU( 0 );  // on at 20 MW, no energy cost: the reserve
  d.lin.assign( d.T , 0 );    // of the DP would be worth -1 per MW
  d.prc.assign( d.T , -1 );
  d.src.assign( d.T , -1 );
  auto d0 = d;
  d0.prho.clear();
  d0.srho.clear();
  d0.prc.clear();
  d0.src.clear();
  const double v = brute_force( d0 , std::vector< int >( d0.T , -1 ) );
  check_all_forms( [ & ]() { return( new_TU( d ) ); } , v ,
                   "priced reserves without reserve Variable" );
  }

 // FixToMaximum, with an outage at 1 - - - - - - - - - - - - - - - - - - -
 check_all_forms( []() {
   return( new_unit_gen( false , 4 ,
           { { "MinPower" , { 1 , 1 , 1 , 1 } } ,
             { "MaxPower" , { 10 , 10 , 10 , 10 } } ,
             { "Availability" , { 1 , 0 , 1 , 1 } } ,
             { "StartUpLimit" , { 10 , 0 , 10 , 10 } } ,
             { "ShutDownLimit" , { 10 , 0 , 10 , 10 } } ,
             { "LinearTerm" , { 5 , 5 , 5 , 5 } } ,
             { "ConstTerm" , { 1 , 1 , 1 , 1 } } } ,
           { { "InitialPower" , 10 } } ,
           { { "InitUpDownTime" , 2 } , { "FixToMaximum" , 1 } } ,
           { { "MinUpTime" , 1 } , { "MinDownTime" , 1 } } ) ); } ,
   154 , "FixToMaximum, outage at 1" );

 // nuclear, PowerBands, outage at 1 - - - - - - - - - - - - - - - - - - - -
 auto nuclear = []( bool deep ) {
  std::vector< std::pair< std::string , std::vector< double > > > v = {
   { "MinPower" , { 10 , 10 , 10 } } , { "MaxPower" , { 50 , 50 , 50 } } ,
   { "Availability" , { 1 , 0 , 1 } } ,
   { "DeltaRampUp" , { 50 , 50 , 50 } } ,
   { "DeltaRampDown" , { 50 , 50 , 50 } } ,
   { "ModulationDeltaRampUp" , { 50 , 50 , 50 } } ,
   { "ModulationDeltaRampDown" , { 50 , 50 , 50 } } ,
   { "StartUpLimit" , { 10 , 0 , 10 } } ,
   { "ShutDownLimit" , { 10 , 0 , 10 } } ,
   { "LinearTerm" , { 1 , 1 , 1 } } ,
   { "ShutDownCost" , { 100 , 100 , 100 } } };
  if( deep ) {
   v.push_back( { "DeepDecreaseThreshold" , { 15 , 15 , 15 } } );
   v.push_back( { "DeepDecreaseGradient" , { 3 , 3 , 3 } } );
   }
  return( new_unit_gen( true , 3 , v , { { "InitialPower" , 15 } } ,
                        { { "InitUpDownTime" , 10 } } ,
                        { { "MinUpTime" , 1 } , { "MinDownTime" , 1 } ,
                          { "ModulationTime" , 2 } ,
                          { "InitModulation" , 2 } } ,
                        deep ? std::vector< double >( )
                             : std::vector< double >( { 20 , 40 } ) ) );
  };
 check_all_forms( [ & ]() { return( nuclear( false ) ); } , 20 ,
                  "nuclear, PowerBands, outage at 1" ,
                  { int( ThermalUnitBlock::TForm ) } );

 // nuclear, deep decrease fixed - - - - - - - - - - - - - - - - - - - - - -
 for( int val : { 0 , 1 } ) {
  auto nub = static_cast< NuclearUnitBlock * >( nuclear( true ) );
  generate_all( nub );
  auto dd = nub->get_deep_decrease();
  if( ! dd ) {
   check( false , "nuclear, deep decrease: no Variable" );
   delete nub;
   continue;
   }
  dd[ 2 ].set_value( val );
  dd[ 2 ].is_fixed( true , eNoMod );
  auto slv = Solver::new_Solver( "NuclearUnitExtDPSolver" );
  nub->register_Solver( slv );
  bool thrown = false;
  try {
   slv->compute();
   }
  catch( std::logic_error & ) {
   thrown = true;
   }
  check( thrown == ( val == 1 ) , "nuclear, deep decrease fixed to " +
         std::to_string( val ) + ( val ? ": not refused" : ": refused" ) );
  nub->unregister_Solver( slv , true );
  delete nub;
  }
 }

/*--------------------------------------------------------------------------*/
/* A Solution of a ThermalUnitBlock saved and written back gives a feasible
 * point of the same value in every formulation: the MILP schedule of a unit
 * that shuts down, restarts and shuts down again (and of one with the same
 * data, on to the end), saved, every Variable of the formulation then set
 * to 0.37, and the Solution written back; ThermalUnitBlock::set_solution()
 * rebuilds the path of the schedule and the power of each run in the five
 * path formulations. With PCuts (bit 8) the Variable of the cuts take the
 * square of the power, i.e., the value is that of the quadratic cost. */

static void test_thermal_solution_restore( void )
{
 const auto sname = milp_solver();
 if( sname.empty() ) {
  std::cout << "no :MILPSolver in this build, the Solution of the "
            << "formulations is not checked" << std::endl;
  return;
  }

 static const std::vector< std::string > names = {
  "u_thermal" , "p_thermal" , "v_thermal" , "w_thermal" , "pr_thermal" ,
  "sc_thermal" , "y_plus_thermal" , "y_minus_thermal" , "p_h_k_thermal" ,
  "p_h_thermal" , "p_k_thermal" , "z_thermal" , "z_h_k_thermal" ,
  "z_h_thermal" , "z_k_thermal" , "teta_thermal" };

 TUData base;
 base.T = 6;
 base.minP = 2;
 base.maxP = 10;
 base.ru = base.rd = 4;
 base.su = base.sd = 6;
 base.initUD = 3;
 base.initP = 5;
 base.minUp = 2;
 base.minDown = 1;
 base.lin = { 4 , 4 , -6 , -6 , 4 , -3 };
 base.quad = { 0.1 , 0.1 , 0.1 , 0.1 , 0.1 , 0.1 };
 base.cnst = { 1 , 1 , 1 , 1 , 1 , 1 };
 base.suc = { 2 , 2 , 2 , 2 , 2 , 2 };
 auto ending = base;  // on to the end of the horizon
 ending.lin = { 4 , 4 , -6 , -6 , -6 , -6 };

 for( const auto * d : { & base , & ending } )
  for( int wf = 0 ; wf < 15 ; ++wf ) {
   if( wf == 7 )
    continue;
   const int form = ( wf < 7 ) ? wf : ( wf - 8 ) | 8;
   auto tub = new_TU( *d );
   tub->set_reserve_vars( 3 );
   generate_all_wf( tub , form );
   const auto who = std::string( "Solution restored, " ) +
                    ( d == & base ? "two runs" : "on to the end" ) +
                    ", wf " + std::to_string( form );
   auto slv = Solver::new_Solver( sname );
   tub->register_Solver( slv );
   const auto status = slv->compute();
   if( status != Solver::kOK ) {
    check( false , who + ": " + sname + " status " +
           std::to_string( status ) );
    tub->unregister_Solver( slv , true );
    delete tub;
    continue;
    }
   slv->get_var_solution();
   auto obj = static_cast< FRealObjective * >( tub->get_objective() );
   obj->compute();
   const double v0 = obj->value();
   auto sol = tub->get_Solution( nullptr , false );
   tub->unregister_Solver( slv , true );

   for( const auto & n : names )
    if( auto vv = tub->get_static_variable_v< ColVariable >( n ) )
     for( auto & var : *vv )
      if( ! var.is_fixed() )
       var.set_value( 0.37 );

   try {
    sol->write( tub );
    }
   catch( std::exception & e ) {
    check( false , who + ": writing the Solution throws " + e.what() );
    }
   delete sol;
   SimpleConfiguration< double > tol( 1e-6 );
   check( tub->is_feasible( true , & tol ) ,
          who + ": the restored point is not feasible" );
   obj->compute();
   if( form < 8 )
    check( close( obj->value() , v0 , 1e-5 ) , who + ": value " +
           str( obj->value() ) + " instead of " + str( v0 ) );
   else  // the cuts give at most the quadratic cost, the restore that
    check( obj->value() >= v0 - 1e-5 , who + ": value " +
           str( obj->value() ) + " below the one of the MILP " + str( v0 ) );
   delete tub;
   }
 }

/*--------------------------------------------------------------------------*/
/* A battery of modules whose converter has no size of its own (no
 * "ConverterMaxPower", no converter design variable): each module has a
 * converter of the power of the module, so that the design that chooses 3
 * modules has 3 times the converter power of one. With and without the
 * binary variables and a round-trip loss (the split of the active power)
 * and with neither (the folded one), 3 modules fixed and up to 3 modules
 * chosen, the value is that of 3 modules of 10 (charge 30, store 27 or 30,
 * discharge as much); the converter rows
 * have the coefficient -10 of the battery design, -20 after set_kappa( 2 ).
 * With "ConverterMaxPower" 10 given the converter stays of one module. */

static void test_battery_converter_modules( void )
{
 const auto sname = milp_solver();
 auto make = []( double mn , double mx , bool conv , bool split ) {
  std::vector< std::pair< std::string , std::vector< double > > > data = {
   { "MinStorage" , { 0 , 0 } } , { "MaxStorage" , { 100 , 100 } } ,
   { "MaxPower" , { 10 , 10 } } , { "MinPower" , { -10 , -10 } } ,
   { "Cost" , { -5 , 5 } } ,
   { "BatteryInvestmentCost" , { 1 } } ,
   { "BatteryMinCapacityDesign" , { mn } } ,
   { "BatteryMaxCapacityDesign" , { mx } } };
  if( conv )
   data.push_back( { "ConverterMaxPower" , { 10 , 10 } } );
  if( split )  // a round-trip loss: the pair is not folded
   data.push_back( { "StoringBatteryRho" , { 0.9 , 0.9 } } );
  auto bu = new_BU( 2 , data );
  if( bu ) {
   SimpleConfiguration< int > negative_prices( split ? 1 : 0 );
   bu->generate_abstract_variables( & negative_prices );
   bu->generate_abstract_constraints();
   bu->generate_objective();
   }
  return( bu );
  };

 for( int split = 0 ; split < 2 ; ++split )
  for( double mx : { 3.0 , -3.0 } ) {
   const double mn = ( mx > 0 ) ? 3 : 0;
   const auto who = std::string( "battery converter per module, " ) +
                    ( split ? "split" : "folded" ) +
                    ( mx > 0 ? ", 3 modules" : ", up to 3 modules" );
   auto bu = make( mn , mx , false , split );
   if( ! bu ) {
    check( false , who + ": not read" );
    continue;
    }
   const auto & rows = bu->get_intake_outtake_design_rows();
   const auto & x = bu->get_const_batt_design();
   const Index nconv = split ? 1 : 2;
   bool ok = ( rows.shape()[ 0 ] == 2 + nconv );
   for( Index r = 2 ; ok && ( r < 2 + nconv ) ; ++r )
    ok = close( row_coef( rows[ r ][ 0 ] , & x ) , -10 );
   check( ok , who + ": no converter row with the coefficient -10 of the "
          "battery design" );
   if( ! sname.empty() ) {
    const double v = milp_value( bu , sname );
    const double e = split ? 3 - 5 * 30 - 5 * 27 : 3 - 5 * 30 - 5 * 30;
    check( close( v , e ) , who + ": value " + str( v ) + " instead of " +
           str( e ) + " (the converter of one module)" );
    }
   bu->set_kappa( 2 );
   ok = ( rows.shape()[ 0 ] == 2 + nconv );
   for( Index r = 2 ; ok && ( r < 2 + nconv ) ; ++r )
    ok = close( row_coef( rows[ r ][ 1 ] , & x ) , -20 );
   check( ok , who + ", set_kappa( 2 ): converter rows not updated" );
   delete bu;
   }

 // ConverterMaxPower given: the converter is of fixed size, one module's
 if( ! sname.empty() ) {
  auto bu = make( 3 , 3 , true , true );
  if( bu ) {
   const double v = milp_value( bu , sname );
   check( close( v , 3 - 5 * 10 - 5 * 9 ) , "battery modules, converter "
          "10 given: value " + str( v ) + " instead of " +
          str( 3 - 5 * 10 - 5 * 9 ) );
   delete bu;
   }
  }
 }

/*--------------------------------------------------------------------------*/
/*--------------- THE MODULATION RAMPS OF A NuclearUnitBlock ---------------*/
/*--------------------------------------------------------------------------*/
/// a nuclear unit on since long at 30 over 6 instants, whose prices ask to
/// move at every instant, with the modulation ramps mru and mrd, modulations
/// of at most L instants and, if deep, deep decreases from 25 by 6

static NuclearUnitBlock * new_ramp_NU( unsigned int L ,
                                       const std::vector< double > & mru ,
                                       const std::vector< double > & mrd ,
                                       bool deep )
{
 std::vector< std::pair< std::string , std::vector< double > > > vecs = {
  { "MinPower" , std::vector< double >( 6 , 10 ) } ,
  { "MaxPower" , std::vector< double >( 6 , 50 ) } ,
  { "DeltaRampUp" , std::vector< double >( 6 , 20 ) } ,
  { "DeltaRampDown" , std::vector< double >( 6 , 20 ) } ,
  { "StartUpLimit" , std::vector< double >( 6 , 50 ) } ,
  { "ShutDownLimit" , std::vector< double >( 6 , 50 ) } ,
  { "ShutDownCost" , std::vector< double >( 6 , 1000 ) } ,
  { "LinearTerm" , { -10 , 10 , -10 , 10 , -10 , 10 } } ,
  { "ModulationDeltaRampUp" , mru } ,
  { "ModulationDeltaRampDown" , mrd } };
 if( deep ) {
  vecs.push_back( { "DeepDecreaseThreshold" ,
                    std::vector< double >( 6 , 25 ) } );
  vecs.push_back( { "DeepDecreaseGradient" ,
                    std::vector< double >( 6 , 6 ) } );
  vecs.push_back( { "DeepDecreaseCost" , std::vector< double >( 6 , 1 ) } );
  }
 auto tub = new_unit_gen( true , 6 , vecs , { { "InitialPower" , 30 } } ,
                          { { "InitUpDownTime" , 10 } } ,
                          { { "MinUpTime" , 1 } , { "MinDownTime" , 1 } ,
                            { "ModulationTime" , 2 } ,
                            { "InitModulation" , 2 } ,
                            { "MaxModulationLength" , L } } );
 auto nub = dynamic_cast< NuclearUnitBlock * >( tub );
 if( ! nub )
  throw( std::logic_error( "new_ramp_NU: no NuclearUnitBlock built" ) );
 return( nub );
 }

/*--------------------------------------------------------------------------*/
/// the value of the NuclearUnitExtDPSolver on the unit (NaN if not kOK)

static double nuclear_DP_value( NuclearUnitBlock * nub )
{
 auto slv = Solver::new_Solver( "NuclearUnitExtDPSolver" );
 nub->register_Solver( slv );
 const auto status = slv->compute();
 const double v = ( status == Solver::kOK ) ? slv->get_ub() : std::nan( "" );
 nub->unregister_Solver( slv , true );
 return( v );
 }

/*--------------------------------------------------------------------------*/
/// the first difference between the static groups of rows of two units
/// built alike, i.e., in the number of groups or of rows, in the sides or in
/// the coefficients, in order, of a row (the Variable being those of each
/// unit); "" if there is none

static std::string rows_differ( Block * a , Block * b )
{
 if( a->get_number_static_constraints() !=
     b->get_number_static_constraints() )
  return( "the numbers of groups differ" );
 for( Index i = 0 ; i < a->get_number_static_constraints() ; ++i ) {
  auto ra = a->get_static_constraint_v< FRowConstraint >( i );
  auto rb = b->get_static_constraint_v< FRowConstraint >( i );
  if( ( ! ra ) || ( ! rb ) ) {
   if( ra || rb )
    return( "group " + std::to_string( i ) + " is a row group in one only" );
   continue;
   }
  const auto g = "group " + std::to_string( i ) + " row ";
  if( ra->size() != rb->size() )
   return( "group " + std::to_string( i ) + " has " +
           std::to_string( ra->size() ) + " rows instead of " +
           std::to_string( rb->size() ) );
  for( Index r = 0 ; r < ra->size() ; ++r ) {
   const auto & x = ( *ra )[ r ];
   const auto & y = ( *rb )[ r ];
   if( ( x.get_lhs() != y.get_lhs() ) || ( x.get_rhs() != y.get_rhs() ) )
    return( g + std::to_string( r ) + ": the sides differ" );
   const auto & vx = static_cast< const LinearFunction * >(
                                          x.get_function() )->get_v_var();
   const auto & vy = static_cast< const LinearFunction * >(
                                          y.get_function() )->get_v_var();
   if( vx.size() != vy.size() )
    return( g + std::to_string( r ) + ": the numbers of Variable differ" );
   for( Index k = 0 ; k < vx.size() ; ++k )
    if( vx[ k ].second != vy[ k ].second )
     return( g + std::to_string( r ) + ": coefficient " +
             std::to_string( k ) + " is " + str( vx[ k ].second ) +
             " instead of " + str( vy[ k ].second ) );
   }
  }
 return( "" );
 }

/*--------------------------------------------------------------------------*/
/* The modulation ramps of a NuclearUnitBlock changed after the generation
 * of the abstract representation, in the T and 3bin formulations with and
 * without the tight operating rules, ramps and cuts, with single-instant
 * modulations and with modulations of 2 instants (the direction of a step
 * matters, and the ramp down is not in the coefficient of m_t as without
 * it): the ramp up at 0 (in the right-hand side of the unit on before the
 * horizon) and at 3 given as an unordered Subset, the ramp down on a Range,
 * the down ramp 0 at an instant included. The unit has then to have the
 * rows of a unit generated from the new data, coefficient by coefficient,
 * and its value, for the :MILPSolver (if
 * any) and for the DP Solver, both the one attached before the change and
 * a new one, whose schedule has to satisfy the changed rows; the old value
 * differs, so that the change is seen. Each setter issues the physical
 * Modification of its own type. A ramp larger than the full ramp is
 * refused, leaving the data and the rows; a dry run of the abstract part
 * changes the data and not the rows. With deep decreases and the tight
 * rules the down ramp decides which instants have the row dd_t <= d_t:
 * a change that adds or removes some of them, all of them included, is
 * refused with the data restored, one that keeps them is done; without
 * the tight rules it is done. */

static void test_nuclear_modulation_ramps_change( void )
{
 const auto sname = milp_solver();
 const std::vector< double > base( 6 , 5 );
 const std::vector< double > mru1 = { 2 , 5 , 5 , 8 , 5 , 5 };
 const std::vector< double > mrd1 = { 5 , 12 , 0 , 7 , 5 , 5 };

 for( unsigned int L : { 1u , 2u } )
  for( int wf : { 1 , 0 , 1 + 32 , 1 + 32 + 64 , 0 + 32 + 64 ,
                  1 + 32 + 64 + 128 } ) {
   const auto what = "modulation ramps change, L " + std::to_string( L ) +
                     ", wf " + std::to_string( wf );

   auto ref_unit = new_ramp_NU( L , mru1 , mrd1 , false );
   generate_all_wf( ref_unit , wf );
   const double ref = nuclear_DP_value( ref_unit );
   if( ! sname.empty() ) {
    const auto v = milp_value( ref_unit , sname );
    check( close( v , ref ) , what + ": the unit generated from the new "
           "data gives " + str( v ) + " with " + sname + " and " +
           str( ref ) + " with the DP" );
    }

   auto nub = new_ramp_NU( L , base , base , false );
   generate_all_wf( nub , wf );
   auto fs = new FakeSolver();
   nub->register_Solver( fs );
   auto dp = Solver::new_Solver( "NuclearUnitExtDPSolver" );
   nub->register_Solver( dp );
   const double old = ( dp->compute() == Solver::kOK ) ? dp->get_ub() :
                                                         std::nan( "" );
   check( ! close( old , ref ) , what + ": the change does not change the "
          "value (" + str( old ) + ")" );
   fs->get_Modification_list().clear();

   std::vector< double > up = { 8 , 2 } , down = { 12 , 0 , 7 };
   nub->set_modulation_ramp_up( up.cbegin() , Block::Subset{ 3 , 0 } ,
                                false );
   nub->set_modulation_ramp_down( down.cbegin() , Range( 1 , 4 ) );
   check( ( nub->get_modulation_ramp_up() == mru1 ) &&
          ( nub->get_modulation_ramp_down() == mrd1 ) ,
          what + ": the ramps are not the new ones" );
   check( ( count_TUB_mods( fs , NuclearUnitBlockMod::eSetModDP ) == 1 ) &&
          ( count_TUB_mods( fs , NuclearUnitBlockMod::eSetModDM ) == 1 ) ,
          what + ": not one physical Modification of each ramp" );
   const auto diff = rows_differ( nub , ref_unit );
   check( diff.empty() , what + ": the rows are not those of the new data, " +
          diff );

   const auto st = dp->compute();
   check( ( st == Solver::kOK ) && close( dp->get_ub() , ref ) ,
          what + ": the DP attached before gives " + str( dp->get_ub() ) +
          " instead of " + str( ref ) );
   nub->unregister_Solver( dp , true );
   nub->unregister_Solver( fs , true );
   if( ! sname.empty() ) {
    const auto v = milp_value( nub , sname );
    check( close( v , ref ) , what + ": " + sname + " gives " + str( v ) +
           " instead of " + str( ref ) );
    }
   check_all_DP( nub , ref , what );

   // larger than the full ramp: refused, nothing changes
   std::vector< double > bad = { 25 };
   bool thrown = false;
   try {
    nub->set_modulation_ramp_down( bad.cbegin() , Range( 2 , 3 ) );
    }
   catch( std::logic_error & ) {
    thrown = true;
    }
   check( thrown && ( nub->get_modulation_ramp_down()[ 2 ] == 0 ) ,
          what + ": a ramp down of 25 > 20 accepted, or the datum changed" );
   if( ! sname.empty() ) {
    const auto v = milp_value( nub , sname );
    check( close( v , ref ) , what + ": after the refusal " + sname +
           " gives " + str( v ) + " instead of " + str( ref ) );
    }

   // a dry run of the abstract part: the datum alone
   std::vector< double > up0 = { 5 };
   nub->set_modulation_ramp_up( up0.cbegin() , Range( 0 , 1 ) , eModBlck ,
                                eDryRun );
   check( nub->get_modulation_ramp_up()[ 0 ] == 5 ,
          what + ": the ramp up not changed by a dry run" );
   if( ! sname.empty() ) {
    const auto v = milp_value( nub , sname );
    check( close( v , ref ) , what + ": after a dry run " + sname +
           " gives " + str( v ) + " instead of " + str( ref ) );
    }
   delete nub;
   delete ref_unit;
   }

 // the rows dd_t <= d_t of the tight rules exist where the deep gradient
 // (6) is larger than the down ramp, i.e., at every instant but 3
 std::vector< double > mrd = { 5 , 5 , 5 , 7 , 5 , 5 };
 for( int wf : { 1 + 32 , 1 } ) {
  const auto what = "modulation ramps change, deep decreases, wf " +
                    std::to_string( wf );
  const bool tight = ( wf & 32 );
  auto nub = new_ramp_NU( 1 , base , mrd , true );
  generate_all_wf( nub , wf );

  auto refused = [ & ]( std::vector< double > v , Range rng ,
                        const std::string & which ) {
   const auto before = nub->get_modulation_ramp_down();
   bool thrown = false;
   try {
    nub->set_modulation_ramp_down( v.cbegin() , rng );
    }
   catch( std::logic_error & ) {
    thrown = true;
    }
   if( tight )
    check( thrown && ( nub->get_modulation_ramp_down() == before ) ,
           what + ": " + which + " accepted, or the data changed" );
   else
    check( ! thrown , what + ": " + which + " refused" );
   };
  refused( { 5 } , Range( 3 , 4 ) , "a row added at 3" );
  refused( { 6 } , Range( 1 , 2 ) , "the row at 1 removed" );
  refused( std::vector< double >( 6 , 6 ) , Range( 0 , 6 ) ,
           "every row removed" );

  // the same instants: done
  std::vector< double > v2 = { 4 };
  nub->set_modulation_ramp_down( v2.cbegin() , Range( 2 , 3 ) );
  auto ref_unit = new_ramp_NU( 1 , base , nub->get_modulation_ramp_down() ,
                               true );
  generate_all_wf( ref_unit , wf );
  const auto diff = rows_differ( nub , ref_unit );
  check( diff.empty() , what + ": the rows are not those of the new data, " +
         diff );
  const double ref = nuclear_DP_value( ref_unit );
  if( ! sname.empty() ) {
   const auto r = milp_value( ref_unit , sname );
   const auto v = milp_value( nub , sname );
   check( close( r , ref ) && close( v , ref ) ,
          what + ": " + sname + " gives " + str( v ) + " (new unit " +
          str( r ) + ") instead of " + str( ref ) );
   }
  check_all_DP( nub , ref , what );
  delete ref_unit;
  delete nub;
  }
 }

/*--------------------------------------------------------------------------*/
/* A fixing of a Variable of the operating rules that the
 * NuclearUnitExtDPSolver refuses (a modulation fixed to 0.5, a deep
 * decrease fixed to 1) makes compute() throw, and the read lock that
 * load_fixings() takes on the unit has to be released all the same. */

static void test_nuclear_DP_fixing_unlocks( void )
{
 const std::vector< double > base( 6 , 5 );
 for( int c = 0 ; c < 2 ; ++c ) {
  const std::string what = c ? "deep decrease fixed to 1" :
                               "modulation fixed to 0.5";
  auto nub = new_ramp_NU( 1 , base , base , c == 1 );
  generate_all_wf( nub , 1 );
  ColVariable * v = c ? nub->get_deep_decrease() : nub->get_modulation();
  if( ! v ) {
   check( false , what + ": Variable missing" );
   delete nub;
   continue;
   }
  v[ 2 ].set_value( c ? 1 : 0.5 );
  v[ 2 ].is_fixed( true , eNoMod );

  auto slv = Solver::new_Solver( "NuclearUnitExtDPSolver" );
  nub->register_Solver( slv );
  bool thrown = false;
  try {
   slv->compute();
   }
  catch( std::logic_error & ) {
   thrown = true;
   }
  check( thrown , what + ": compute() does not throw" );
  check( nub->is_owned_by( nullptr ) ,
         what + ": the unit is left locked after the throw" );
  if( nub->is_owned_by( nullptr ) )
   nub->unregister_Solver( slv , true );
  delete nub;
  }
 }

/*--------------------------------------------------------------------------*/
/*------------- THE ROWS OF A DERIVED CLASS IN update_rows() ---------------*/
/*--------------------------------------------------------------------------*/
/* A ThermalUnitBlock with rows of its own, written as a derived class writes
 * them [see ThermalUnitBlock::build_rows()]: p_t <= MaxPower[ t ] - 5 by
 * push_row() at the instants t < 2 where MaxPower[ t ] > 50, the group
 * registered only if it has rows, and the same rows in a dynamic group at
 * the instants t >= 2, registered only if it has rows when generated. */

class PushRowsTUB : public ThermalUnitBlock {
 public:
  PushRowsTUB( void ) : ThermalUnitBlock( nullptr ) {}

  ~PushRowsTUB() override {
   Constraint::clear( pushed );
   Constraint::clear( late.rows );
   }

  std::vector< FRowConstraint > pushed;
  DynRows late;

 protected:
  void build_rows( bool generate_ZOConstraints ) override {
   ThermalUnitBlock::build_rows( generate_ZOConstraints );
   auto p = get_active_power( 0 );
   pushed.reserve( 2 );  // Constraint cannot be copied
   for( Index t = 0 ; t < 2 ; ++t )
    if( get_max_power( t ) > 50 )
     push_row( pushed , { std::make_pair( p + t , 1.0 ) } ,
               -Inf< double >() , get_max_power( t ) - 5 );
   if( ! pushed.empty() )
    add_rows( pushed , "Pushed_Test" );
   for( Index t = 2 ; t < get_time_horizon() ; ++t )
    if( get_max_power( t ) > 50 )
     put_dyn_row( late , { int( t ) } , { std::make_pair( p + t , 1.0 ) } ,
                  -Inf< double >() , get_max_power( t ) - 5 );
   if( ! late.rows.empty() )
    add_dyn_rows( late , "Late_Test" );
   }
 };

/*--------------------------------------------------------------------------*/
/* MaxPower 60 60 40 40, i.e., two pushed rows and no dynamic group: 70 at 1
 * keeps the count and changes the row; 45 at 0 (one row less) and 45 at 0
 * and 1 (no row left) are refused, as 55 at 3 (a row of a dynamic group
 * that is not there) is, all with the data and the rows kept. */

static void test_thermal_update_rows_groups( void )
{
 auto g = new_group( "TUpush" , true );
 g.addDim( "TimeHorizon" , 4 );
 auto NI = g.addDim( "NumberIntervals" , 4 );
 put( g , "MinPower" , NI , { 10 , 10 , 10 , 10 } );
 put( g , "MaxPower" , NI , { 60 , 60 , 40 , 40 } );
 put( g , "LinearTerm" , NI , { -10 , -10 , -10 , -10 } );
 put_int( g , "InitUpDownTime" , 1 );
 put( g , "InitialPower" , 20 );
 PushRowsTUB tub;
 tub.deserialize( g );
 generate_all_wf( & tub , 0 );
 check( tub.pushed.size() == 2 , "derived rows: " +
        std::to_string( tub.pushed.size() ) + " pushed rows, not 2" );

 auto refused = [ & ]( std::vector< double > mp , Range rng ,
                       const std::string & what ) {
  const auto old = tub.get_max_power();
  const auto n = tub.pushed.size();
  std::vector< double > rhs;
  for( const auto & r : tub.pushed )
   rhs.push_back( r.get_rhs() );
  bool thrown = false;
  try {
   tub.set_maximum_power( mp.begin() , rng );
   }
  catch( std::logic_error & ) {
   thrown = true;
   }
  bool kept = ( tub.get_max_power() == old ) && ( tub.pushed.size() == n ) &&
              tub.late.rows.empty();
  for( Index i = 0 ; kept && ( i < n ) ; ++i )
   kept = ( tub.pushed[ i ].get_rhs() == rhs[ i ] );
  check( thrown && kept , "derived rows, " + what + ": " +
         ( thrown ? "the data or the rows changed" : "accepted" ) );
  };

 std::vector< double > mp = { 70 };
 tub.set_maximum_power( mp.begin() , Range( 1 , 2 ) );
 check( ( tub.pushed.size() == 2 ) && ( tub.pushed[ 1 ].get_rhs() == 65 ) &&
        ( tub.pushed[ 0 ].get_rhs() == 55 ) ,
        "derived rows: MaxPower 70 at 1 does not give the rhs 65 there" );

 refused( { 45 } , Range( 0 , 1 ) , "one pushed row less" );
 refused( { 45 , 45 } , Range( 0 , 2 ) , "no pushed row left" );
 refused( { 55 } , Range( 3 , 4 ) , "a row of a dynamic group not there" );
 }

/*--------------------------------------------------------------------------*/
/*------------------ InitUpDownTime AFTER THE GENERATION -------------------*/
/*--------------------------------------------------------------------------*/
/* With MinUpTime 3 and MinDownTime 2 over 6 instants and prices that ask to
 * switch, in the seven formulations: a unit on since 5 becomes on since 8
 * or since 3 (free from 0 in all three), and a unit off since 4 becomes off
 * since 2 (free from 0 in both): these are done, the value of a :MILPSolver
 * attached before the change and of the DP Solvers being that of the unit
 * read with the new value. On since 5 to on since 2 or off since 2, off
 * since 4 to off since 1 or on since 3 change the Variable and are refused
 * with the datum kept, as on since 1 to on since 2 is with MinUpTime 9,
 * which the old state cuts to 7 (fixed on for the whole horizon in both,
 * but the minimum time of the data, hence the one of the new unit, is not
 * known); a change whose abstract part is a dry run is done. */

static void test_thermal_init_updown_change( void )
{
 const auto sname = milp_solver();
 struct Case { int from; int to; unsigned int min_up; bool done; };
 const std::vector< Case > cases = {
  { 5 , 8 , 3 , true } , { 5 , 3 , 3 , true } , { -4 , -2 , 3 , true } ,
  { 5 , 2 , 3 , false } , { 5 , -2 , 3 , false } , { -4 , -1 , 3 , false } ,
  { -4 , 3 , 3 , false } , { 1 , 2 , 9 , false } };

 for( int wf = 0 ; wf < 7 ; ++wf )
  for( const auto & c : cases ) {
   const std::vector< std::pair< std::string , std::vector< double > > >
    data = { { "MinPower" , std::vector< double >( 6 , 10 ) } ,
             { "MaxPower" , std::vector< double >( 6 , 100 ) } ,
             { "DeltaRampUp" , std::vector< double >( 6 , 40 ) } ,
             { "DeltaRampDown" , std::vector< double >( 6 , 40 ) } ,
             { "LinearTerm" , { -10 , 20 , -10 , 20 , -10 , -10 } } };
   const auto what = "InitUpDownTime " + std::to_string( c.from ) + " to " +
                     std::to_string( c.to ) + ", MinUpTime " +
                     std::to_string( c.min_up ) + ", wf " +
                     std::to_string( wf );
   auto tub = new_TU_vec( 6 , data , c.from , 50 , c.min_up , 2 );
   generate_all_wf( tub , wf );

   Solver * ms = sname.empty() ? nullptr : Solver::new_Solver( sname );
   if( ms ) {
    tub->register_Solver( ms );
    ms->compute();
    }

   std::vector< int > v = { c.to };
   bool thrown = false;
   try {
    tub->set_init_updown_time( v.cbegin() , Range( 0 , 1 ) );
    }
   catch( std::logic_error & ) {
    thrown = true;
    }
   const int now = c.done ? c.to : c.from;
   check( ( thrown != c.done ) && ( tub->get_init_up_down_time() == now ) ,
          what + ": " + ( thrown ? "refused" : "done" ) +
          ", InitUpDownTime " +
          std::to_string( tub->get_init_up_down_time() ) );

   // the unit read with the value it should have
   auto fresh = new_TU_vec( 6 , data , now , 50 , c.min_up , 2 );
   generate_all_wf( fresh , wf );
   double expected = std::nan( "" );
   if( ! sname.empty() ) {
    expected = milp_value( fresh , sname );
    double value = std::nan( "" );
    if( ms->compute() == Solver::kOK )
     value = ms->get_var_value();
    check( close( value , expected ) , what + ": the attached " + sname +
           " finds " + str( value ) + " instead of " + str( expected ) );
    tub->unregister_Solver( ms , true );
    }
   else {
    auto dp = Solver::new_Solver( "ThermalUnitDPSolver" );
    fresh->register_Solver( dp );
    if( dp->compute() == Solver::kOK )
     expected = dp->get_ub();
    fresh->unregister_Solver( dp , true );
    }
   check_all_DP_wf( tub , expected , what , wf );
   delete fresh;

   // a dry run of the abstract part: the datum alone
   if( ( wf == 0 ) && ( ! c.done ) ) {
    tub->set_init_updown_time( v.cbegin() , Range( 0 , 1 ) , eModBlck ,
                               eDryRun );
    check( tub->get_init_up_down_time() == c.to ,
           what + ": not changed by a dry run of the abstract part" );
    }
   delete tub;
   }
 }

/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*---------------- THERMAL AND NUCLEAR: FINAL VERIFICATION -----------------*/
/*--------------------------------------------------------------------------*/
/* A change of the start-up or shut-down costs after a solve reaches a DP
 * Solver that stays attached: off before the horizon with a cheap start-up
 * the unit runs at every instant, and a start-up cost of 200 keeps it off;
 * on before the horizon with a cheap shut-down it shuts down at 0, and a
 * shut-down cost of 100 keeps it on at the minimum power. Each value is the
 * brute force of the data, after a first compute() of the same Solver. */

static void test_thermal_DP_cost_change( void )
{
 for( int sd = 0 ; sd < 2 ; ++sd )
  for( const auto & sname : THERMAL_DP ) {
   TUData d;
   d.T = 3;
   d.maxP = d.su = d.sd = d.ru = d.rd = 5;
   if( sd ) {
    d.initUD = 10;
    d.initP = 5;
    d.lin = { 10 , 10 , 10 };
    d.sdc = { 1 , 1 , 1 };
    }
   else {
    d.initUD = -10;
    d.lin = { -10 , -10 , -10 };
    d.suc = { 1 , 1 , 1 };
    }
   const std::string what = std::string( "DP cost change, " ) +
                            ( sd ? "shut-down" : "start-up" ) + ", " + sname;
   const std::vector< int > nofix( d.T , -1 );
   auto tub = new_TU( d );
   auto slv = Solver::new_Solver( sname );
   tub->register_Solver( slv );
   auto st = slv->compute();
   check( ( st == Solver::kOK ) &&
          close( slv->get_ub() , brute_force( d , nofix ) ) ,
          what + ": before the change, " + str( slv->get_ub() ) );

   const std::vector< double > c( d.T , sd ? 100 : 200 );
   if( sd ) {
    tub->set_shutdown_costs( c.cbegin() );
    d.sdc = c;
    }
   else {
    tub->set_startup_costs( c.cbegin() );
    d.suc = c;
    }
   const double expected = brute_force( d , nofix );
   check( ! close( expected , slv->get_ub() ) ,
          what + ": the change does not change the optimum" );
   st = slv->compute();
   check( ( st == Solver::kOK ) && close( slv->get_ub() , expected ) ,
          what + ": after the change " + str( slv->get_ub() ) +
          " instead of " + str( expected ) );
   tub->unregister_Solver( slv , true );
   delete tub;
   }
 }

/*--------------------------------------------------------------------------*/
/* A fixed Variable of an extended formulation, which the DP Solvers cannot
 * honour, is refused by them with std::logic_error, also for a nuclear
 * unit, while the same unit with no such fixing is solved; the MILP takes
 * the fixing, here an arc of the DP formulation fixed to 0, which leaves
 * the optimum as it is. */

static void test_thermal_DP_extended_fixing( void )
{
 const auto sname = milp_solver();
 for( int nuc = 0 ; nuc < 2 ; ++nuc ) {
  TUData d;
  d.T = 3;
  d.maxP = d.su = d.sd = d.ru = d.rd = 5;
  d.initUD = -10;
  d.lin = { -10 , -10 , -10 };
  d.nuclear = nuc;
  if( nuc )
   d.mru = d.mrd = 5;
  const std::string what = std::string( "DP extended fixing, " ) +
                           ( nuc ? "nuclear" : "thermal" );
  auto tub = new_TU( d );
  generate_all_wf( tub , ThermalUnitBlock::DPForm );
  auto y = tub->get_static_variable_v< ColVariable >( "y_plus_thermal" );
  check( y && ( ! y->empty() ) , what + ": no y_plus_thermal" );
  if( ( ! y ) || y->empty() ) {
   delete tub;
   continue;
   }
  check_all_DP( tub , brute_force( d , std::vector< int >( d.T , -1 ) ) ,
                what + ", not fixed" );

  // the last arc, a run that is never the optimal one
  y->back().set_value( 0 );
  y->back().is_fixed( true );
  std::vector< std::string > names = THERMAL_DP;
  if( nuc )
   names = { "NuclearUnitExtDPSolver" };
  for( const auto & s : names ) {
   auto slv = Solver::new_Solver( s );
   tub->register_Solver( slv );
   bool refused = false;
   try {
    slv->compute();
    }
   catch( std::logic_error & ) {
    refused = true;
    }
   check( refused , what + ", " + s + ": the fixing is not refused" );
   tub->unregister_Solver( slv , true );
   }
  if( ! sname.empty() ) {
   const auto v = milp_value( tub , sname );
   check( close( v , -150 ) , what + ": " + sname + " gives " + str( v ) +
          " instead of -150" );
   }
  delete tub;
  }
 }

/*--------------------------------------------------------------------------*/
/* set_initial_power() refuses, before any change, a negative value and,
 * for a unit on before the horizon, a value below MinPower[ 0 ], as
 * check_data_consistency() does: the initial power stays as it was, with
 * and without the abstract representation; for a unit off before the
 * horizon a value below the minimum power is taken. */

static void test_thermal_initial_power_refused( void )
{
 for( int gen = 0 ; gen < 2 ; ++gen )
  for( int on = 0 ; on < 2 ; ++on ) {
   TUData d;
   d.T = 3;
   d.minP = 2;
   d.maxP = d.su = d.sd = d.ru = d.rd = 5;
   d.initUD = on ? 3 : -3;
   d.initP = on ? 4 : 0;
   d.lin = { -10 , -10 , -10 };
   const std::string what = std::string( "initial power refused, " ) +
                            ( on ? "on" : "off" ) +
                            ( gen ? ", generated" : "" );
   auto tub = new_TU( d );
   if( gen )
    generate_all( tub );
   const double old = tub->get_initial_power();

   auto refused = [ & ]( double v ) {
    std::vector< double > val = { v };
    try {
     tub->set_initial_power( val.cbegin() );
     }
    catch( std::logic_error & ) {
     return( true );
     }
    return( false );
    };

   check( refused( -1 ) && ( tub->get_initial_power() == old ) ,
          what + ": a negative value is not refused" );
   if( on )
    check( refused( 1 ) && ( tub->get_initial_power() == old ) ,
           what + ": 1 < MinPower is not refused" );
   else
    check( ( ! refused( 1 ) ) && ( tub->get_initial_power() == 1 ) ,
           what + ": 1 is not taken" );
   check( ( ! refused( 3 ) ) && ( tub->get_initial_power() == 3 ) ,
          what + ": 3 is not taken" );
   delete tub;
   }
 }

/*--------------------------------------------------------------------------*/
/* is_feasible() checks every row, those of FixToMaximum (p_t at the
 * operational maximum power) included: a unit on before the horizon, fixed
 * to its maximum 10, at 10 is feasible and at 5 is not, in the 3bin and in
 * the T formulation. */

static void test_thermal_is_feasible_fix_to_max( void )
{
 for( int wf : { 0 , 1 } ) {
  auto tub = new_unit_gen( false , 1 ,
                           { { "MinPower" , { 1 } } , { "MaxPower" , { 10 } } ,
                             { "DeltaRampUp" , { 10 } } ,
                             { "DeltaRampDown" , { 10 } } ,
                             { "StartUpLimit" , { 10 } } ,
                             { "ShutDownLimit" , { 10 } } } ,
                           { { "InitialPower" , 10 } } ,
                           { { "InitUpDownTime" , 5 } ,
                             { "FixToMaximum" , 1 } } ,
                           { { "MinUpTime" , 1 } , { "MinDownTime" , 1 } } );
  generate_all_wf( tub , wf );
  const auto what = "is_feasible, FixToMaximum, wf " + std::to_string( wf );
  auto u = tub->get_commitment( 0 );
  auto p = tub->get_active_power( 0 );
  auto v = tub->get_start_up();
  auto w = tub->get_shut_down();
  u[ 0 ].set_value( 1 );
  if( v )
   v[ 0 ].set_value( 0 );
  if( w )
   w[ 0 ].set_value( 0 );
  p[ 0 ].set_value( 10 );
  check( tub->is_feasible( true ) , what + ": the maximum is not feasible" );
  p[ 0 ].set_value( 5 );
  check( ! tub->is_feasible( true ) , what + ": 5 is feasible" );
  delete tub;
  }
 }

/*--------------------------------------------------------------------------*/
/* The costs of the downward modulation steps and of the deep decreases set
 * through an unordered Subset go to the instants of the Subset: the values
 * 7 and 3 at the instants 3 and 0 give 3 at 0 and 7 at 3. */

static void test_nuclear_rule_costs_unordered( void )
{
 TUData d;
 d.T = 4;
 d.maxP = d.su = d.sd = d.ru = d.rd = 10;
 d.initUD = -10;
 d.lin = { -10 , -10 , -10 , -10 };
 d.nuclear = true;
 d.mru = d.mrd = 2;
 auto g = new_group( "NU" , true );
 write_TU( g , d );
 auto NI = g.getDim( "NumberIntervals" );
 put( g , "DownModulationCost" , NI , { 1 , 1 , 1 , 1 } );
 auto nub = dynamic_cast< NuclearUnitBlock * >( Block::new_Block( g ) );
 check( nub , "nuclear rule costs, unordered: not read" );
 if( ! nub )
  return;
 const std::vector< double > val = { 7 , 3 };
 nub->set_down_modulation_costs( val.cbegin() , Block::Subset{ 3 , 0 } ,
                                 false );
 const auto & c = nub->get_down_modulation_cost();
 check( ( c.size() == 4 ) && ( c[ 0 ] == 3 ) && ( c[ 3 ] == 7 ) &&
        ( c[ 1 ] == 1 ) && ( c[ 2 ] == 1 ) ,
        "nuclear rule costs, unordered: the costs are not 3 , 1 , 1 , 7" );
 delete nub;
 }

/*--------------------------------------------------------------------------*/
/* The bounds of the reactive node injections that a UCBlock gives to an
 * ACNetworkBlock cover a ThermalUnitBlock on as well as off: with
 * Q^mn = -5, Q^mx = 5, Q^{mn,on} = -10 and Q^{mx,on} = 20 the unit at node
 * 0 can inject from -15 to 25, which the rows of node 0 have to allow. */

static void test_reactive_injection_bounds_on( void )
{
 auto g = new_group( "UCQ" , true );
 g.putAtt( "type" , "UCBlock" );
 auto TH = g.addDim( "TimeHorizon" , 1 );
 g.addDim( "NumberUnits" , 2 );
 auto N = g.addDim( "NumberNodes" , 2 );
 auto L = g.addDim( "NumberLines" , 1 );
 auto G = g.addDim( "NumberElectricalGenerators" , 2 );
 {
  const char * cls = "ACNetworkBlock";
  g.addVar( "NetworkBlockClassname" , netCDF::NcString() ).putVar( & cls );
  }
 put_int( g , "StartLine" , L , { 0 } );
 put_int( g , "EndLine" , L , { 1 } );
 put( g , "MinPowerFlow" , L , { -30 } );
 put( g , "MaxPowerFlow" , L , { 30 } );
 put( g , "LineReactance" , L , { 0.1 } );
 put( g , "LineResistance" , L , { 0.01 } );
 put( g , "LineMinAngle" , L , { -30 } );
 put( g , "LineMaxAngle" , L , { 30 } );
 put( g , "LineRATEA" , L , { 30 } );
 put( g , "NodeConductance" , N , { 0 , 0 } );
 put( g , "NodeSusceptance" , N , { 0 , 0 } );
 put( g , "NodeMaxVoltage" , N , { 1.1 , 1.1 } );
 put( g , "NodeMinVoltage" , N , { 0.9 , 0.9 } );
 put( g , "ActivePowerDemand" , { N , TH } , { 1 , 2 } );
 put( g , "ReactivePowerDemand" , { N , TH } , { 0.5 , 0.5 } );
 put_int( g , "GeneratorNode" , G , { 0 , 1 } );

 auto u = g.addGroup( "UnitBlock_0" );
 u.putAtt( "type" , "ThermalUnitBlock" );
 put( u , "MinPower" , 1.0 );
 put( u , "MaxPower" , 10.0 );
 put( u , "MinReactivePower" , -5.0 );
 put( u , "MaxReactivePower" , 5.0 );
 put( u , "MinReactivePowerOn" , -10.0 );
 put( u , "MaxReactivePowerOn" , 20.0 );
 put_int( u , "InitUpDownTime" , -5 );
 put_uint( u , "MinUpTime" , 1 );
 put_uint( u , "MinDownTime" , 1 );

 auto s = g.addGroup( "UnitBlock_1" );
 s.putAtt( "type" , "SlackUnitBlock" );
 put( s , "MaxPower" , 1000.0 );
 put( s , "ActivePowerCost" , 1000.0 );

 std::unique_ptr< UCBlock > uc;
 try {
  uc.reset( dynamic_cast< UCBlock * >( Block::new_Block( g ) ) );
  }
 catch( std::exception & e ) {
  check( false , std::string( "reactive injection bounds: " ) + e.what() );
  return;
  }
 check( uc != nullptr , "reactive injection bounds: the UCBlock is not read" );
 if( ! uc )
  return;
 generate_all( uc.get() );
 auto nb = uc->get_network_block( 0 );
 auto b = nb->get_static_constraint_v< BoxConstraint >(
                              "Reactive_Node_Injection_Bound_Const_Network" );
 const bool ok = b && ( b->size() == 2 );
 check( ok && close( ( *b )[ 0 ].get_lhs() , -15 ) &&
        close( ( *b )[ 0 ].get_rhs() , 25 ) ,
        "reactive injection bounds: node 0 has [ " +
        str( ok ? ( *b )[ 0 ].get_lhs() : std::nan( "" ) ) + " , " +
        str( ok ? ( *b )[ 0 ].get_rhs() : std::nan( "" ) ) +
        " ], not [ -15 , 25 ]" );
 }

/*--------------------------------------------------------------------------*/
/* An ACNetworkBlock whose data do not have "LineMinAngle" and
 * "LineMaxAngle": the angle differences are not bounded, so the rows
 * AC_angle_bounds_limit and AC_elem_bounds are not there, the strengthened
 * relaxation is not generated by default and refused if asked, the data
 * written keep the bounds out, and one of the two without the other is
 * refused. A UCBlock on such a network, a ThermalUnitBlock and a
 * SlackUnitBlock on two nodes, is solved by the :MILPSolver of the build,
 * and its value is at most that with the bounds of +-30 degrees. */

static void test_ac_no_angle_bounds( void )
{
 auto net = [ & ]( bool min_angle , bool max_angle ) {
  auto g = new_group( "ACA" , true );
  g.putAtt( "type" , "ACNetworkBlock" );
  auto N = g.addDim( "NumberNodes" , 3 );
  auto L = g.addDim( "NumberLines" , 2 );
  put_int( g , "StartLine" , L , { 0 , 1 } );
  put_int( g , "EndLine" , L , { 1 , 2 } );
  put( g , "MinPowerFlow" , L , { -100 , -100 } );
  put( g , "MaxPowerFlow" , L , { 100 , 100 } );
  put( g , "LineRATEA" , L , { 100 , 100 } );
  put( g , "LineReactance" , L , { 0.1 , 0.2 } );
  put( g , "LineResistance" , L , { 0.01 , 0.02 } );
  if( min_angle )
   put( g , "LineMinAngle" , L , { -30 , -30 } );
  if( max_angle )
   put( g , "LineMaxAngle" , L , { 30 , 30 } );
  put( g , "NodeConductance" , N , { 0 , 0 , 0 } );
  put( g , "NodeSusceptance" , N , { 0 , 0 , 0 } );
  put( g , "NodeMaxVoltage" , N , { 1.1 , 1.1 , 1.1 } );
  put( g , "NodeMinVoltage" , N , { 0.9 , 0.9 , 0.9 } );
  put( g , "ReactiveDemand" , N , { 0.5 , 0.5 , 0.5 } );
  put( g , "ActiveDemand" , N , { 1 , 2 , 3 } );
  return( g );
  };

 // one of the two bounds without the other is refused
 for( bool min_angle : { true , false } ) {
  const std::string what = std::string( "AC angles, only " ) +
                           ( min_angle ? "LineMinAngle" : "LineMaxAngle" );
  // Block::new_Block() reports the exception and returns nullptr
  bool refused = false;
  try {
   std::unique_ptr< Block > b( Block::new_Block( net( min_angle ,
                                                      ! min_angle ) ) );
   refused = ! b;
   }
  catch( std::invalid_argument & ) {
   refused = true;
   }
  catch( std::exception & e ) {
   check( false , what + ": throws " + e.what() );
   continue;
   }
  check( refused , what + ": accepted" );
  }

 // neither: no bound on the angle differences
 std::string what = "AC angles absent";
 try {
  std::unique_ptr< ACNetworkBlock > ab( dynamic_cast< ACNetworkBlock * >(
                              Block::new_Block( net( false , false ) ) ) );
  check( ab != nullptr , what + ": not read" );
  if( ! ab )
   return;
  auto nd = static_cast< ACNetworkBlock::ACNetworkData * >(
                                                      ab->get_NetworkData() );
  check( ( ! nd->has_angle_bounds() ) && nd->get_line_min_angle().empty() &&
         nd->get_line_max_angle().empty() ,
         what + ": the data have bounds on the angle differences" );

  // the data written have no bound either, and are read back
  auto w = new_group( "ACW" );
  ab->serialize( w );
  check( w.getVar( "LineMinAngle" ).isNull() &&
         w.getVar( "LineMaxAngle" ).isNull() ,
         what + ": the bounds are written" );
  std::unique_ptr< ACNetworkBlock > rb( dynamic_cast< ACNetworkBlock * >(
                                                   Block::new_Block( w ) ) );
  check( rb && ( ! static_cast< ACNetworkBlock::ACNetworkData * >(
                            rb->get_NetworkData() )->has_angle_bounds() ) ,
         what + ": not read back without bounds" );

  for( Index n = 0 ; n < 3 ; ++n ) {
   ab->set_min_node_injection( -10 , n );
   ab->set_max_node_injection( 10 , n );
   ab->set_min_reactive_node_injection( -10 , n , 0 );
   ab->set_max_reactive_node_injection( 10 , n , 0 );
   }
  generate_all( ab.get() );
  check( ( ! ab->get_static_variable_v< ColVariable >( "v_theta" ) ) &&
         ( ! ab->get_static_constraint_v< FRowConstraint >(
                                              "v_theta_bounds" ) ) ,
         what + ": the strengthened relaxation is generated" );
  check( ! ab->get_static_constraint_v< FRowConstraint >(
                                             "AC_angle_bounds_limit" ) ,
         what + ": the rows AC_angle_bounds_limit are there" );
  check( ! ab->get_static_constraint_v< BoxConstraint >( "AC_elem_bounds" ) ,
         what + ": the bounds AC_elem_bounds are there" );
  auto socp = ab->get_static_constraint_v< FRowConstraint >(
                                                         "AC_socp_const" );
  check( socp && ( socp->size() == 2 ) ,
         what + ": no cone for each of the 2 lines" );
  }
 catch( std::exception & e ) {
  check( false , what + ": throws " + e.what() );
  }

 // the strengthened relaxation asked for is refused, before any Variable
 what = "AC angles absent, strengthened relaxation asked";
 try {
  std::unique_ptr< ACNetworkBlock > ab( dynamic_cast< ACNetworkBlock * >(
                              Block::new_Block( net( false , false ) ) ) );
  SimpleConfiguration< int > strong( 1 );
  bool thrown = false;
  try {
   ab->generate_abstract_variables( & strong );
   }
  catch( std::invalid_argument & ) {
   thrown = true;
   }
  check( thrown , what + ": accepted" );
  check( ab->get_number_static_variables() == 0 ,
         what + ": " + std::to_string( ab->get_number_static_variables() ) +
         " groups of Variable added before refusing" );
  }
 catch( std::exception & e ) {
  check( false , what + ": throws " + e.what() );
  }

 // the UCBlock, with and without the bounds
 const auto sname = milp_solver();
 if( sname.empty() ) {
  std::cout << "AC angles absent: no :MILPSolver, solve skipped"
            << std::endl;
  return;
  }
 auto uc_value = [ & ]( bool angles , const std::string & what ) {
  auto g = new_group( "UCA" , true );
  g.putAtt( "type" , "UCBlock" );
  auto TH = g.addDim( "TimeHorizon" , 1 );
  g.addDim( "NumberUnits" , 2 );
  auto N = g.addDim( "NumberNodes" , 2 );
  auto L = g.addDim( "NumberLines" , 1 );
  auto G = g.addDim( "NumberElectricalGenerators" , 2 );
  {
   const char * cls = "ACNetworkBlock";
   g.addVar( "NetworkBlockClassname" , netCDF::NcString() ).putVar( & cls );
   }
  put_int( g , "StartLine" , L , { 0 } );
  put_int( g , "EndLine" , L , { 1 } );
  put( g , "MinPowerFlow" , L , { -30 } );
  put( g , "MaxPowerFlow" , L , { 30 } );
  put( g , "LineReactance" , L , { 0.1 } );
  put( g , "LineResistance" , L , { 0.01 } );
  if( angles ) {
   put( g , "LineMinAngle" , L , { -30 } );
   put( g , "LineMaxAngle" , L , { 30 } );
   }
  put( g , "LineRATEA" , L , { 30 } );
  put( g , "NodeConductance" , N , { 0 , 0 } );
  put( g , "NodeSusceptance" , N , { 0 , 0 } );
  put( g , "NodeMaxVoltage" , N , { 1.1 , 1.1 } );
  put( g , "NodeMinVoltage" , N , { 0.9 , 0.9 } );
  put( g , "ActivePowerDemand" , { N , TH } , { 1 , 2 } );
  put( g , "ReactivePowerDemand" , { N , TH } , { 0.5 , 0.5 } );
  put_int( g , "GeneratorNode" , G , { 0 , 1 } );

  auto u = g.addGroup( "UnitBlock_0" );
  u.putAtt( "type" , "ThermalUnitBlock" );
  put( u , "MinPower" , 1.0 );
  put( u , "MaxPower" , 10.0 );
  put( u , "LinearTerm" , 1.0 );
  put( u , "MinReactivePower" , -5.0 );
  put( u , "MaxReactivePower" , 5.0 );
  put_int( u , "InitUpDownTime" , -5 );
  put_uint( u , "MinUpTime" , 1 );
  put_uint( u , "MinDownTime" , 1 );

  auto s = g.addGroup( "UnitBlock_1" );
  s.putAtt( "type" , "SlackUnitBlock" );
  put( s , "MaxPower" , 1000.0 );
  put( s , "ActivePowerCost" , 1000.0 );

  double v = std::nan( "" );
  try {
   std::unique_ptr< UCBlock > uc( dynamic_cast< UCBlock * >(
                                                  Block::new_Block( g ) ) );
   if( ! uc ) {
    check( false , what + ": the UCBlock is not read" );
    return( v );
    }
   generate_all( uc.get() );
   v = milp_value( uc.get() , sname );
   }
  catch( std::exception & e ) {
   check( false , what + ": throws " + e.what() );
   }
  return( v );
  };

 const double unb = uc_value( false , "AC angles absent, UCBlock" );
 const double bounded = uc_value( true , "AC angles of +-30, UCBlock" );
 check( std::isfinite( unb ) && std::isfinite( bounded ) &&
        ( unb <= bounded + 1e-6 * std::max( 1.0 , std::abs( bounded ) ) ) ,
        "AC angles absent, UCBlock: value " + str( unb ) +
        ", with the bounds " + str( bounded ) );
 }

/*--------------------------------------------------------------------------*/
/* Whether a maximum power row of the unit at instant t, among the static
 * rows of MaxPower_Const_Thermal and the dynamic ones of MaxPower5 and
 * MaxPower6, has exactly the given coefficients (u_t and p_t included). */

static bool has_max_row( ThermalUnitBlock * tub , const std::vector<
                         std::pair< const Variable * , double > > & coef ,
                         double eps = 1e-9 )
{
 auto match = [ & ]( const FRowConstraint & c ) {
  auto lf = dynamic_cast< const LinearFunction * >( c.get_function() );
  if( ( ! lf ) || ( lf->get_num_active_var() != coef.size() ) )
   return( false );
  for( const auto & vc : coef ) {
   const double a = row_coef( c , vc.first );
   if( std::isnan( a ) || ( std::abs( a - vc.second ) > eps ) )
    return( false );
   }
  return( true );
  };
 if( auto v = tub->get_static_constraint_v< FRowConstraint >(
                                              "MaxPower_Const_Thermal" ) )
  for( const auto & c : *v )
   if( match( c ) )
    return( true );
 for( const std::string n : { "MaxPower5_Const_Thermal" ,
                              "MaxPower6_Const_Thermal" } )
  if( auto l = tub->get_dynamic_constraint< FRowConstraint >( n ) )
   for( const auto & c : *l )
    if( match( c ) )
     return( true );
 return( false );
 }

/*--------------------------------------------------------------------------*/
/* An ACNetworkBlock that is a chain of AC lines, line l from node l to node
 * l + 1, with the given bounds on the angle differences (in degrees, NaN
 * written as such) and voltage bounds that change from node to node; no
 * "LineName" unless names is true. flows are the "MinPowerFlow" and
 * "MaxPowerFlow" of every line. */

static netCDF::NcGroup write_ac_chain(
                const std::vector< std::pair< double , double > > & angles ,
                bool names = false ,
                std::pair< double , double > flows = { -100 , 100 } )
{
 const Index nl = angles.size();
 auto g = new_group( "ACC" , true );
 g.putAtt( "type" , "ACNetworkBlock" );
 auto N = g.addDim( "NumberNodes" , nl + 1 );
 auto L = g.addDim( "NumberLines" , nl );
 std::vector< int > st( nl ) , en( nl );
 std::vector< double > mn( nl ) , mx( nl ) , fl( nl , flows.first ) ,
                       fu( nl , flows.second ) , ra( nl , 100 ) ,
                       x( nl ) , r( nl );
 for( Index l = 0 ; l < nl ; ++l ) {
  st[ l ] = l;
  en[ l ] = l + 1;
  mn[ l ] = angles[ l ].first;
  mx[ l ] = angles[ l ].second;
  x[ l ] = 0.1 + 0.01 * l;
  r[ l ] = 0.01 + 0.001 * l;
  }
 std::vector< double > vmn( nl + 1 ) , vmx( nl + 1 ) , zero( nl + 1 , 0 ) ,
                       dem( nl + 1 , 0.1 );
 for( Index n = 0 ; n <= nl ; ++n ) {
  vmn[ n ] = 0.85 + 0.01 * ( n % 5 );
  vmx[ n ] = 1.05 + 0.02 * ( n % 3 );
  }
 put_int( g , "StartLine" , L , st );
 put_int( g , "EndLine" , L , en );
 put( g , "MinPowerFlow" , L , fl );
 put( g , "MaxPowerFlow" , L , fu );
 put( g , "LineRATEA" , L , ra );
 put( g , "LineReactance" , L , x );
 put( g , "LineResistance" , L , r );
 put( g , "LineMinAngle" , L , mn );
 put( g , "LineMaxAngle" , L , mx );
 if( names ) {
  auto v = g.addVar( "LineName" , netCDF::NcString() , L );
  for( Index l = 0 ; l < nl ; ++l ) {
   const std::string s = "L" + std::to_string( l );
   const char * p = s.c_str();
   v.putVar( { l } , & p );
   }
  }
 put( g , "NodeConductance" , N , zero );
 put( g , "NodeSusceptance" , N , zero );
 put( g , "NodeMaxVoltage" , N , vmx );
 put( g , "NodeMinVoltage" , N , vmn );
 put( g , "ReactiveDemand" , N , dem );
 put( g , "ActiveDemand" , N , dem );
 return( g );
 }

/*--------------------------------------------------------------------------*/
/* The bounds on the angle differences of an ACNetworkBlock, read with the
 * convention of MATPOWER and turned into rows (8) and bounds (8a) line by
 * line [see ACNetworkBlock::generate_abstract_constraints()]: +-360 and
 * beyond, NaN and 0 / 0 leave a side (or both) unbounded, and such a line
 * has no row and no bound; a line with a range over 180 degrees has the
 * bounds and no row; all the others have both. Each row and each bound is
 * checked against a sampling of the box of the angle difference and of
 * the voltages, independent of the formulas: no sampled point violates
 * them, and every bound is the extreme of the sample (it is attained). The
 * strengthened relaxation covers the lines within +-90 degrees only. */

static void test_ac_angle_bounds_matpower( void )
{
 const double nan = std::nan( "" );
 const double pi = std::numbers::pi;
 // the cases, and whether the line has the bounds (8a) and the rows (8)
 struct Case { double mn , mx; bool box , rows , qc; };
 const std::vector< Case > cases = {
  { -360 , 360 , false , false , false } ,
  { -400 , 30 , false , false , false } ,
  { -30 , 360 , false , false , false } ,
  { 0 , 0 , false , false , false } ,
  { nan , 20 , false , false , false } ,
  { -80 , 80 , true , true , true } ,     // the rows divided by cos 80
  { 10 , 40 , true , true , true } ,      // not containing 0
  { -60 , -20 , true , true , true } ,    // not containing 0, negative
  { -120 , 120 , true , false , false } , // range over 180
  { 100 , 170 , true , true , false } ,   // c < 0 , s > 0
  { -30 , 30 , true , true , true } ,
  { -180 , 180 , false , false , false } ,// a whole turn
  { -170 , -100 , true , true , false } , // c < 0 , s < 0
  { 0 , 0.5 , true , true , true } ,      // one side 0
  { -10 , 200 , true , false , false } }; // over a maximizer of sin
 std::vector< std::pair< double , double > > angles;
 for( const auto & c : cases )
  angles.emplace_back( c.mn , c.mx );

 const std::string what = "AC angles, MATPOWER convention";
 std::unique_ptr< ACNetworkBlock > ab;
 try {
  ab.reset( dynamic_cast< ACNetworkBlock * >(
                              Block::new_Block( write_ac_chain( angles ) ) ) );
  }
 catch( std::exception & e ) {
  check( false , what + ": throws " + e.what() );
  return;
  }
 check( ab != nullptr , what + ": not read" );
 if( ! ab )
  return;
 auto nd = static_cast< ACNetworkBlock::ACNetworkData * >(
                                                      ab->get_NetworkData() );
 for( Index l = 0 ; l < cases.size() ; ++l ) {
  const auto [ mn , mx ] = nd->get_angle_difference_bounds( l );
  const bool lo = std::isfinite( mn ) , up = std::isfinite( mx );
  const bool exp_lo = ( cases[ l ].mn > -360 ) && ( ! ( ( cases[ l ].mn == 0 )
                      && ( cases[ l ].mx == 0 ) ) );
  const bool exp_up = ( cases[ l ].mx < 360 ) && ( ! ( ( cases[ l ].mn == 0 )
                      && ( cases[ l ].mx == 0 ) ) );
  check( ( lo == exp_lo ) && ( up == exp_up ) &&
         ( ( ! lo ) || ( mn == cases[ l ].mn ) ) &&
         ( ( ! up ) || ( mx == cases[ l ].mx ) ) ,
         what + ": line " + std::to_string( l ) + " [ " +
         str( cases[ l ].mn ) + " , " + str( cases[ l ].mx ) + " ] read as [ "
         + str( mn ) + " , " + str( mx ) + " ]" );
  }

 for( Index n = 0 ; n <= cases.size() ; ++n ) {
  ab->set_min_node_injection( -10 , n );
  ab->set_max_node_injection( 10 , n );
  ab->set_min_reactive_node_injection( -10 , n , 0 );
  ab->set_max_reactive_node_injection( 10 , n , 0 );
  }
 try {
  generate_all( ab.get() );
  }
 catch( std::exception & e ) {
  check( false , what + ": generating throws " + e.what() );
  return;
  }

 auto rows = ab->get_static_constraint< FRowConstraint , 2 >(
                                                    "AC_angle_bounds_limit" );
 auto box = ab->get_static_constraint< BoxConstraint , 2 >(
                                                           "AC_elem_bounds" );
 auto c = ab->get_static_variable_v< ColVariable >(
                                                  "v_sum_product_voltages" );
 auto s = ab->get_static_variable_v< ColVariable >(
                                                 "v_diff_product_voltages" );
 Index nrows = 0 , nbox = 0;
 for( const auto & cs : cases ) {
  nrows += cs.rows;
  nbox += cs.box;
  }
 check( rows && ( rows->shape()[ 0 ] == 2 ) &&
        ( rows->shape()[ 1 ] == nrows ) ,
        what + ": not " + std::to_string( nrows ) + " pairs of rows (8)" );
 check( box && ( box->shape()[ 0 ] == 2 ) && ( box->shape()[ 1 ] == nbox ) ,
        what + ": not " + std::to_string( nbox ) + " pairs of bounds (8a)" );
 check( c && s && ( c->size() == cases.size() ) &&
        ( s->size() == cases.size() ) , what + ": no c and s" );
 if( ( ! rows ) || ( ! box ) || ( ! c ) || ( ! s ) ||
     ( rows->shape()[ 1 ] != nrows ) || ( box->shape()[ 1 ] != nbox ) )
  return;

 const auto & vmn = nd->get_node_min_voltage();
 const auto & vmx = nd->get_node_max_voltage();
 Index ir = 0 , ib = 0;
 for( Index l = 0 ; l < cases.size() ; ++l ) {
  const auto & cs = cases[ l ];
  const std::string wl = what + ": line " + std::to_string( l ) + " [ " +
                         str( cs.mn ) + " , " + str( cs.mx ) + " ]";
  // no row and no bound of another line uses c_l and s_l
  auto uses = [ & ]( Index i ) {
   return( ( ! std::isnan( row_coef( ( *rows )[ 0 ][ i ] , & ( *s )[ l ] ) ) )
           || ( ( *box )[ 0 ][ i ].get_active_var( 0 ) == & ( *c )[ l ] ) );
   };
  if( cs.rows ) {
   check( uses( ir ) , wl + ": the row (8) is not that of the line" );
   if( ! uses( ir ) )
    return;
   }
  if( cs.box &&
      ( ( *box )[ 0 ][ ib ].get_active_var( 0 ) != & ( *c )[ l ] ) ) {
   check( false , wl + ": the bound (8a) is not that of the line" );
   return;
   }
  if( ! cs.box )
   continue;

  // the sample: the angle on a grid of 0.05 degrees, the voltages at
  // their bounds and in the middle
  const double a = pi * cs.mn / 180 , b = pi * cs.mx / 180;
  const Index ns = l , ne = l + 1;
  std::vector< double > z;
  for( double vs : { vmn[ ns ] , vmx[ ns ] , ( vmn[ ns ] + vmx[ ns ] ) / 2 } )
   for( double ve : { vmn[ ne ] , vmx[ ne ] , ( vmn[ ne ] + vmx[ ne ] ) / 2 } )
    z.push_back( vs * ve );
  const Index K = std::ceil( ( cs.mx - cs.mn ) / 0.05 );
  double cmn = INF , cmx = -INF , smn = INF , smx = -INF;
  bool row_ok = true , box_ok = true;
  const double tol = 1e-9;
  for( Index k = 0 ; k <= K ; ++k ) {
   const double th = a + ( b - a ) * k / K;
   for( double zz : z ) {
    const double cv = zz * std::cos( th ) , sv = zz * std::sin( th );
    cmn = std::min( cmn , cv );
    cmx = std::max( cmx , cv );
    smn = std::min( smn , sv );
    smx = std::max( smx , sv );
    for( int h = 0 ; h < 2 ; ++h ) {
     const auto & bc = ( *box )[ h ][ ib ];
     const double v = h ? sv : cv;
     if( ( v < bc.get_lhs() - tol ) || ( v > bc.get_rhs() + tol ) )
      box_ok = false;
     if( cs.rows ) {
      const auto & rc = ( *rows )[ h ][ ir ];
      const double rv = row_coef( rc , & ( *s )[ l ] ) * sv +
                        row_coef( rc , & ( *c )[ l ] ) * cv;
      if( ( rv < rc.get_lhs() - tol ) || ( rv > rc.get_rhs() + tol ) )
       row_ok = false;
      }
     }
    }
   }
  check( box_ok , wl + ": a feasible point violates the bounds (8a)" );
  check( row_ok , wl + ": a feasible point violates the rows (8)" );
  const double eps = 1e-5;
  check( ( std::abs( ( *box )[ 0 ][ ib ].get_lhs() - cmn ) <= eps ) &&
         ( std::abs( ( *box )[ 0 ][ ib ].get_rhs() - cmx ) <= eps ) &&
         ( std::abs( ( *box )[ 1 ][ ib ].get_lhs() - smn ) <= eps ) &&
         ( std::abs( ( *box )[ 1 ][ ib ].get_rhs() - smx ) <= eps ) ,
         wl + ": bounds c in [ " + str( ( *box )[ 0 ][ ib ].get_lhs() ) +
         " , " + str( ( *box )[ 0 ][ ib ].get_rhs() ) + " ], s in [ " +
         str( ( *box )[ 1 ][ ib ].get_lhs() ) + " , " +
         str( ( *box )[ 1 ][ ib ].get_rhs() ) + " ], sampled [ " + str( cmn ) +
         " , " + str( cmx ) + " ] and [ " + str( smn ) + " , " + str( smx ) +
         " ]" );
  if( cs.rows ) {
   // written s - tan( phi ) c if cos( phi ) >= 1 / 2, rotated otherwise
   for( int h = 0 ; h < 2 ; ++h ) {
    const double phi = h ? b : a;
    const double cs_exp = std::cos( phi ) >= 0.5 ? 1 : std::cos( phi );
    check( close( row_coef( ( *rows )[ h ][ ir ] , & ( *s )[ l ] ) ,
                  cs_exp , 1e-12 ) ,
           wl + ": the row (8) on side " + std::to_string( h ) +
           " has coefficient " +
           str( row_coef( ( *rows )[ h ][ ir ] , & ( *s )[ l ] ) ) +
           " of s, not " + str( cs_exp ) );
    }
   ++ir;
   }
  ++ib;
  }

 // the strengthened relaxation (on by default) only covers the lines within
 // +-90 degrees: v_theta_bounds has their bounds, in order
 auto tb = ab->get_static_constraint_v< FRowConstraint >( "v_theta_bounds" );
 auto al = ab->get_static_variable_v< ColVariable >( "v_alpha" );
 std::vector< Index > qc;
 for( Index l = 0 ; l < cases.size() ; ++l )
  if( cases[ l ].qc )
   qc.push_back( l );
 check( tb && al && ( tb->size() == qc.size() ) &&
        ( al->size() == qc.size() ) ,
        what + ": the strengthened relaxation does not cover " +
        std::to_string( qc.size() ) + " lines" );
 if( tb && ( tb->size() == qc.size() ) )
  for( Index i = 0 ; i < qc.size() ; ++i )
   check( close( ( *tb )[ i ].get_lhs() , pi * cases[ qc[ i ] ].mn / 180 ,
                 1e-12 ) &&
          close( ( *tb )[ i ].get_rhs() , pi * cases[ qc[ i ] ].mx / 180 ,
                 1e-12 ) ,
          what + ": v_theta_bounds " + std::to_string( i ) +
          " is not that of line " + std::to_string( qc[ i ] ) );

 // the dynamic cuts of the strengthened relaxation, at the point 0, use
 // only the Variable of the lines covered
 try {
  ab->generate_dynamic_constraints();
  }
 catch( std::exception & e ) {
  check( false , what + ": separating the cuts throws " + e.what() );
  }
 }

/*--------------------------------------------------------------------------*/
/* Bounds on the angle differences that are not read as such: a network
 * whose lines are all +-360, or 0 / 0, has no row (8), no bound (8a) and
 * no strengthened relaxation, which is refused if asked for; a finite pair
 * with LineMinAngle > LineMaxAngle is refused; a UCBlock on such a network
 * has the value of the same network without "LineMinAngle" and
 * "LineMaxAngle", and one with +-80 degrees a value between that and the
 * one with +-30 degrees. */

static void test_ac_angle_bounds_unbounded( void )
{
 for( const auto & p : std::vector< std::pair< double , double > >{
                       { -360 , 360 } , { 0 , 0 } , { -720 , 400 } } ) {
  const std::string what = "AC angles [ " + str( p.first ) + " , " +
                           str( p.second ) + " ]";
  try {
   std::unique_ptr< ACNetworkBlock > ab( dynamic_cast< ACNetworkBlock * >(
                         Block::new_Block( write_ac_chain( { p , p } ) ) ) );
   check( ab != nullptr , what + ": not read" );
   if( ! ab )
    continue;
   for( Index n = 0 ; n < 3 ; ++n ) {
    ab->set_min_node_injection( -10 , n );
    ab->set_max_node_injection( 10 , n );
    ab->set_min_reactive_node_injection( -10 , n , 0 );
    ab->set_max_reactive_node_injection( 10 , n , 0 );
    }
   generate_all( ab.get() );
   check( ! ab->get_static_constraint< FRowConstraint , 2 >(
                                              "AC_angle_bounds_limit" ) ,
          what + ": the rows (8) are there" );
   check( ! ab->get_static_constraint< BoxConstraint , 2 >(
                                                     "AC_elem_bounds" ) ,
          what + ": the bounds (8a) are there" );
   check( ! ab->get_static_variable_v< ColVariable >( "v_theta" ) ,
          what + ": the strengthened relaxation is generated" );
   }
  catch( std::exception & e ) {
   check( false , what + ": throws " + e.what() );
   }
  try {
   std::unique_ptr< ACNetworkBlock > ab( dynamic_cast< ACNetworkBlock * >(
                         Block::new_Block( write_ac_chain( { p , p } ) ) ) );
   SimpleConfiguration< int > strong( 1 );
   bool thrown = false;
   try {
    ab->generate_abstract_variables( & strong );
    }
   catch( std::invalid_argument & ) {
    thrown = true;
    }
   check( thrown , what + ": the strengthened relaxation asked is accepted" );
   }
  catch( std::exception & e ) {
   check( false , what + ", strengthened asked: throws " + e.what() );
   }
  }

 // LineMinAngle > LineMaxAngle, both finite
 {
  bool refused = false;
  try {
   std::unique_ptr< Block > b( Block::new_Block( write_ac_chain(
                                     { { -30 , 30 } , { 20 , 10 } } ) ) );
   refused = ! b;
   }
  catch( std::invalid_argument & ) {
   refused = true;
   }
  catch( std::exception & e ) {
   check( false , std::string( "AC angles [ 20 , 10 ]: throws " ) +
          e.what() );
   }
  check( refused , "AC angles [ 20 , 10 ]: accepted" );
  }

 const auto sname = milp_solver();
 if( sname.empty() ) {
  std::cout << "AC angles, unbounded: no :MILPSolver, solve skipped"
            << std::endl;
  return;
  }

 // the UCBlock of test_ac_no_angle_bounds(), the line given the angles
 auto uc_value = [ & ]( bool angles , double mn , double mx ) {
  const std::string what = "AC angles " + ( angles ? "[ " + str( mn ) +
                           " , " + str( mx ) + " ]" : std::string( "absent" ) )
                           + ", UCBlock";
  auto g = new_group( "UCM" , true );
  g.putAtt( "type" , "UCBlock" );
  auto TH = g.addDim( "TimeHorizon" , 1 );
  g.addDim( "NumberUnits" , 2 );
  auto N = g.addDim( "NumberNodes" , 2 );
  auto L = g.addDim( "NumberLines" , 1 );
  auto G = g.addDim( "NumberElectricalGenerators" , 2 );
  {
   const char * cls = "ACNetworkBlock";
   g.addVar( "NetworkBlockClassname" , netCDF::NcString() ).putVar( & cls );
   }
  put_int( g , "StartLine" , L , { 0 } );
  put_int( g , "EndLine" , L , { 1 } );
  put( g , "MinPowerFlow" , L , { -30 } );
  put( g , "MaxPowerFlow" , L , { 30 } );
  put( g , "LineReactance" , L , { 0.1 } );
  put( g , "LineResistance" , L , { 0.01 } );
  if( angles ) {
   put( g , "LineMinAngle" , L , { mn } );
   put( g , "LineMaxAngle" , L , { mx } );
   }
  put( g , "LineRATEA" , L , { 30 } );
  put( g , "NodeConductance" , N , { 0 , 0 } );
  put( g , "NodeSusceptance" , N , { 0 , 0 } );
  put( g , "NodeMaxVoltage" , N , { 1.1 , 1.1 } );
  put( g , "NodeMinVoltage" , N , { 0.9 , 0.9 } );
  put( g , "ActivePowerDemand" , { N , TH } , { 1 , 9 } );
  put( g , "ReactivePowerDemand" , { N , TH } , { 0.5 , 0.5 } );
  put_int( g , "GeneratorNode" , G , { 0 , 1 } );

  auto u = g.addGroup( "UnitBlock_0" );
  u.putAtt( "type" , "ThermalUnitBlock" );
  put( u , "MinPower" , 1.0 );
  put( u , "MaxPower" , 20.0 );
  put( u , "LinearTerm" , 1.0 );
  put( u , "MinReactivePower" , -5.0 );
  put( u , "MaxReactivePower" , 5.0 );
  put_int( u , "InitUpDownTime" , -5 );
  put_uint( u , "MinUpTime" , 1 );
  put_uint( u , "MinDownTime" , 1 );

  auto sl = g.addGroup( "UnitBlock_1" );
  sl.putAtt( "type" , "SlackUnitBlock" );
  put( sl , "MaxPower" , 1000.0 );
  put( sl , "ActivePowerCost" , 1000.0 );

  double v = std::nan( "" );
  try {
   std::unique_ptr< UCBlock > uc( dynamic_cast< UCBlock * >(
                                                  Block::new_Block( g ) ) );
   if( ! uc ) {
    check( false , what + ": the UCBlock is not read" );
    return( v );
    }
   generate_all( uc.get() );
   v = milp_value( uc.get() , sname );
   }
  catch( std::exception & e ) {
   check( false , what + ": throws " + e.what() );
   }
  check( std::isfinite( v ) , what + ": value " + str( v ) );
  return( v );
  };

 auto same = [ & ]( double a , double b ) {
  return( std::abs( a - b ) <= 1e-6 * std::max( 1.0 , std::abs( b ) ) );
  };
 const double absent = uc_value( false , 0 , 0 );
 const double m360 = uc_value( true , -360 , 360 );
 const double zero = uc_value( true , 0 , 0 );
 const double b80 = uc_value( true , -80 , 80 );
 const double b30 = uc_value( true , -30 , 30 );
 std::cout << "AC angles, UCBlock values: absent " << absent << ", +-360 "
           << m360 << ", 0 / 0 " << zero << ", +-80 " << b80 << ", +-30 "
           << b30 << std::endl;
 check( same( m360 , absent ) , "AC angles +-360, UCBlock: value " +
        str( m360 ) + ", without the bounds " + str( absent ) );
 check( same( zero , absent ) , "AC angles 0 / 0, UCBlock: value " +
        str( zero ) + ", without the bounds " + str( absent ) );
 check( ( absent <= b80 + 1e-6 * std::max( 1.0 , std::abs( b80 ) ) ) &&
        ( b80 <= b30 + 1e-6 * std::max( 1.0 , std::abs( b30 ) ) ) ,
        "AC angles +-80, UCBlock: value " + str( b80 ) + " not between " +
        str( absent ) + " (no bound) and " + str( b30 ) + " (+-30)" );
 }

/*--------------------------------------------------------------------------*/
/* The elementary check of the flows of ACNetworkBlock prints the name of a
 * line whose flow bounds cannot be met: with no "LineName" in the data it
 * prints an empty name, with the names the name of the line. */

static void test_ac_flow_check_no_names( void )
{
 for( bool names : { false , true } ) {
  const std::string what = std::string( "AC flow check, " ) +
                           ( names ? "with" : "without" ) + " line names";
  std::ostringstream out;
  auto old = std::cout.rdbuf( out.rdbuf() );
  try {
   // flows that the AC equations cannot give
   std::unique_ptr< ACNetworkBlock > ab( dynamic_cast< ACNetworkBlock * >(
       Block::new_Block( write_ac_chain( { { -30 , 30 } , { -360 , 360 } } ,
                                         names , { -2e6 , -1e6 } ) ) ) );
   if( ab ) {
    for( Index n = 0 ; n < 3 ; ++n ) {
     ab->set_min_node_injection( -10 , n );
     ab->set_max_node_injection( 10 , n );
     ab->set_min_reactive_node_injection( -10 , n , 0 );
     ab->set_max_reactive_node_injection( 10 , n , 0 );
     }
    generate_all( ab.get() );
    }
   std::cout.rdbuf( old );
   check( ab != nullptr , what + ": not read" );
   }
  catch( std::exception & e ) {
   std::cout.rdbuf( old );
   check( false , what + ": throws " + e.what() );
   continue;
   }
  const auto s = out.str();
  for( Index l = 0 ; l < 2 ; ++l ) {
   const std::string line = " The power line with index = " +
                            std::to_string( l ) + " and name " +
                            ( names ? "L" + std::to_string( l ) : "" ) +
                            " has induced bounds";
   check( s.find( line ) != std::string::npos ,
          what + ": the check of line " + std::to_string( l ) +
          " does not print \"" + line + "\"" );
   }
  }
 }

/*--------------------------------------------------------------------------*/
/* The maximum power rows (24)-(26) of the T formulation: a unit started at
 * t - s reaches at most StartUpLimit plus s ramps at t, and one shut down
 * at t + 1 + s at most ShutDownLimit plus s ramps. With P = 54.27,
 * SU = 43.41 and a ramp of 6.78 (MinUpTime 3) the row (24) at t = 2 has
 * the coefficients 10.86 of v_2 and 4.08 of v_1; with constant data
 * P = 100, SU = SD = 20 and ramps of 30 the rows are (38), (40) and (41) of
 * the T formulation of Knueven, Ostrowski and Watson, e.g., (24) at t = 3
 * with MinUpTime 4 is p_3 <= 100 u_3 - 80 w_4 - 80 v_3 - 50 v_2 - 20 v_1. */

static void test_thermal_T_ramp_bound_rows( void )
{
 {  // the hand values
  auto tub = new_unit_gen( false , 4 ,
                           { { "MinPower" , { 1 , 1 , 1 , 1 } } ,
                             { "MaxPower" , std::vector< double >( 4 ,
                                                                  54.27 ) } ,
                             { "StartUpLimit" , std::vector< double >( 4 ,
                                                                 43.41 ) } ,
                             { "ShutDownLimit" , std::vector< double >( 4 ,
                                                                 54.27 ) } ,
                             { "DeltaRampUp" , std::vector< double >( 4 ,
                                                                  6.78 ) } } ,
                           { } , { { "InitUpDownTime" , -5 } } ,
                           { { "MinUpTime" , 3 } , { "MinDownTime" , 1 } } );
  generate_all_wf( tub , ThermalUnitBlock::TForm );
  auto u = tub->get_commitment( 0 );
  auto p = tub->get_active_power( 0 );
  auto v = tub->get_start_up();
  auto w = tub->get_shut_down();
  check( has_max_row( tub , { { &u[ 2 ] , 54.27 } , { &p[ 2 ] , -1 } ,
                              { &w[ 3 ] , 0 } , { &v[ 2 ] , -10.86 } ,
                              { &v[ 1 ] , -4.08 } } , 1e-6 ) ,
         "T rows (24): no row p_2 <= 54.27 u_2 - 10.86 v_2 - 4.08 v_1" );
  delete tub;
  }

 for( unsigned int up : { 3u , 4u } ) {  // KOW (38), (40), (41)
  auto tub = new_unit_gen( false , 6 ,
                           { { "MinPower" , std::vector< double >( 6 , 1 ) } ,
                             { "MaxPower" , std::vector< double >( 6 ,
                                                                   100 ) } ,
                             { "StartUpLimit" , std::vector< double >( 6 ,
                                                                    20 ) } ,
                             { "ShutDownLimit" , std::vector< double >( 6 ,
                                                                     20 ) } ,
                             { "DeltaRampUp" , std::vector< double >( 6 ,
                                                                   30 ) } ,
                             { "DeltaRampDown" , std::vector< double >( 6 ,
                                                                     30 ) } } ,
                           { } , { { "InitUpDownTime" , -5 } } ,
                           { { "MinUpTime" , up } , { "MinDownTime" , 1 } } );
  generate_all_wf( tub , ThermalUnitBlock::TForm );
  auto u = tub->get_commitment( 0 );
  auto p = tub->get_active_power( 0 );
  auto v = tub->get_start_up();
  auto w = tub->get_shut_down();
  const auto what = "T rows, MinUpTime " + std::to_string( up );
  if( up == 4 ) {
   check( has_max_row( tub , { { &u[ 3 ] , 100 } , { &p[ 3 ] , -1 } ,
                               { &w[ 4 ] , -80 } , { &v[ 3 ] , -80 } ,
                               { &v[ 2 ] , -50 } , { &v[ 1 ] , -20 } } ) ,
          what + ": (24) at 3 is not KOW (38)" );
   check( has_max_row( tub , { { &u[ 1 ] , 100 } , { &p[ 1 ] , -1 } ,
                               { &w[ 2 ] , -80 } , { &w[ 3 ] , -50 } ,
                               { &w[ 4 ] , -20 } , { &v[ 1 ] , -80 } } ) ,
          what + ": (26) at 1 is not KOW (41)" );
   }
  else {
   check( has_max_row( tub , { { &u[ 3 ] , 100 } , { &p[ 3 ] , -1 } ,
                               { &w[ 4 ] , -80 } , { &v[ 3 ] , -80 } ,
                               { &v[ 2 ] , -50 } } ) ,
          what + ": (24) at 3 is not KOW (38)" );
   check( has_max_row( tub , { { &u[ 3 ] , 100 } , { &p[ 3 ] , -1 } ,
                               { &v[ 3 ] , -80 } , { &v[ 2 ] , -50 } ,
                               { &v[ 1 ] , -20 } } ) ,
          what + ": (25) at 3 is not KOW (40)" );
   check( has_max_row( tub , { { &u[ 1 ] , 100 } , { &p[ 1 ] , -1 } ,
                               { &w[ 2 ] , -80 } , { &w[ 3 ] , -50 } ,
                               { &w[ 4 ] , -20 } } ) ,
          what + ": (26) at 1 is not KOW (41)" );
   }
  delete tub;
  }
 }

/*--------------------------------------------------------------------------*/
/* The shut-down limit of the run that the state before the horizon
 * imposes: on since 1 instant with MinUpTime 3, the unit stays on at 0 and
 * 1 and may shut down at 2 only from at most ShutDownLimit[ 2 ] = 2 at 1;
 * with no ramp, a profit at 0 and 1 and a cost at 2, the optimum is -120
 * (10 at 0, 2 at 1, off at 2), in every formulation and with the DP
 * Solvers. */

static void test_thermal_shut_down_limit_at_t0( void )
{
 auto make = []( void ) {
  return( new_unit_gen( false , 3 ,
                        { { "MinPower" , { 1 , 1 , 1 } } ,
                          { "MaxPower" , { 10 , 10 , 10 } } ,
                          { "StartUpLimit" , { 10 , 10 , 10 } } ,
                          { "ShutDownLimit" , { 2 , 2 , 2 } } ,
                          { "LinearTerm" , { -10 , -10 , 100 } } } ,
                        { { "InitialPower" , 5 } } ,
                        { { "InitUpDownTime" , 1 } } ,
                        { { "MinUpTime" , 3 } , { "MinDownTime" , 1 } } ) );
  };
 check_all_forms( make , -120 , "shut-down limit at t0" );
 }


/*--------------------------------------------------------------------------*/
/* Tests of fvfix-HYD, to be pasted in UCBlock/test/test.cpp after
 * test_slack_unit() (they use write_HU(), new_HU(), turbine_HU(),
 * hydro_point(), new_BU(), obj_coef(), generate_all(), check(), close(),
 * str(), new_group(), put()). Extra #include: none (<cmath> is there).
 *
 * Call lines for main(), after test_slack_unit();
 *
 *  test_hydro_set_inflow_subset();
 *  test_hydro_initial_flow_absent();
 *  test_hydro_cost_without_datum();
 *  test_hydro_volume_bound_alone();
 *  test_single_value_costs();
 *  test_battery_design_ub();
 *  test_battery_design_zero_max();
 *  test_battery_kappa_infinite_storage();
 */
/*--------------------------------------------------------------------------*/
/* set_inflow( Subset ) takes the indices n T + t of all the instants, not
 * only those below the number of reservoirs, and refuses those from R T on;
 * the balance of the instant follows. */

static void test_hydro_set_inflow_subset( void )
{
 std::string why;
 auto d = turbine_HU( 3 );
 d.initV = { 5 };
 auto b = new_HU( d , why );
 if( ! b ) {
  check( false , "hydro set_inflow( Subset ): throws " + why );
  return;
  }
 generate_all( b );
 std::vector< double > val = { 4 };
 try {
  b->set_inflow( val.cbegin() , Subset( { 2 } ) );
  check( b->get_inflows()[ 0 ][ 2 ] == 4 ,
	 "hydro set_inflow( { 2 } ): inflow at 2 is " +
	 str( b->get_inflows()[ 0 ][ 2 ] ) + " instead of 4" );
  check( hydro_point( b , { { 0 , 0 , 0 } } , { { 0 , 0 , 0 } } ,
		      { { 5 , 5 , 9 } } ) ,
	 "hydro set_inflow( { 2 } ): the inflow 4 is not in the balance "
	 "of instant 2" );
  }
 catch( std::exception & e ) {
  check( false , std::string( "hydro set_inflow( { 2 } ) with R = 1, "
			      "T = 3: " ) + e.what() );
  }
 bool thrown = false;
 try {
  b->set_inflow( val.cbegin() , Subset( { 3 } ) );
  }
 catch( std::invalid_argument & ) {
  thrown = true;
  }
 check( thrown , "hydro set_inflow( { 3 } ) with R T = 3 accepted" );
 delete b;
 }

/*--------------------------------------------------------------------------*/
/* With "InitialFlowRate" absent, set_initial_flow_rate() on some of the
 * arcs gives the vector one entry per arc: the arcs not set keep 0 (a vector
 * of size 1 would be read as the value of every arc, a shorter one read out
 * of range). */

static void test_hydro_initial_flow_absent( void )
{
 HUData d;
 d.T = 2;
 d.A = 3;
 d.start = { 0 , 0 , 0 };
 d.end = { 1 , 1 , 1 };
 d.maxF.assign( 6 , 10 );
 d.maxP.assign( 6 , 10 );
 d.lin = { 1 , 1 , 1 };
 d.maxV.assign( 2 , 100 );
 std::string why;
 auto b = new_HU( d , why );
 if( ! b ) {
  check( false , "hydro initial flow, absent: throws " + why );
  return;
  }
 std::vector< double > val = { 2 };
 b->set_initial_flow_rate( val.cbegin() , Subset( { 0 } ) );
 check( ( b->get_initial_flow_rate( 0 ) == 2 ) &&
	( b->get_initial_flow_rate( 1 ) == 0 ) &&
	( b->get_initial_flow_rate( 2 ) == 0 ) ,
	"hydro set_initial_flow_rate( { 0 } , 2 ), datum absent: flows " +
	str( b->get_initial_flow_rate( 0 ) ) + ", " +
	str( b->get_initial_flow_rate( 1 ) ) + ", " +
	str( b->get_initial_flow_rate( 2 ) ) + " instead of 2, 0, 0" );
 delete b;

 b = new_HU( d , why );
 val = { 1 , 3 };
 b->set_initial_flow_rate( val.cbegin() , Range( 0 , 2 ) );
 check( ( b->get_initial_flow_rate( 0 ) == 1 ) &&
	( b->get_initial_flow_rate( 1 ) == 3 ) &&
	( b->get_initial_flow_rate( 2 ) == 0 ) ,
	"hydro set_initial_flow_rate( [ 0 , 2 ) ), datum absent: flow of "
	"arc 2 " + str( b->get_initial_flow_rate( 2 ) ) + " instead of 0" );
 delete b;
 }

/*--------------------------------------------------------------------------*/
/* An Objective generated without "ActivePowerCost" has no term in the
 * active power: set_active_power_cost() then adds the terms of the arcs it
 * sets, at every instant, instead of throwing after changing the data. */

static void test_hydro_cost_without_datum( void )
{
 for( bool range : { false , true } ) {
  const std::string w = range ? "( Range )" : "( Subset )";
  std::string why;
  auto b = new_HU( turbine_HU( 2 ) , why );
  if( ! b ) {
   check( false , "hydro cost without datum: throws " + why );
   return;
   }
  generate_all( b );
  std::vector< double > val = { 7 };
  try {
   if( range )
    b->set_active_power_cost( val.cbegin() , Range( 0 , 1 ) );
   else
    b->set_active_power_cost( val.cbegin() , Subset( { 0 } ) );
   check( close( obj_coef( b , b->get_active_power( 0 , 0 ) ) , 7 ) &&
	  close( obj_coef( b , b->get_active_power( 0 , 1 ) ) , 7 ) ,
	  "hydro set_active_power_cost" + w + " without the datum: "
	  "coefficients " + str( obj_coef( b , b->get_active_power( 0 , 0 ) ) )
	  + " and " + str( obj_coef( b , b->get_active_power( 0 , 1 ) ) ) +
	  " instead of 7" );
   val = { 2 };
   if( range )
    b->set_active_power_cost( val.cbegin() , Range( 0 , 1 ) );
   else
    b->set_active_power_cost( val.cbegin() , Subset( { 0 } ) );
   check( close( obj_coef( b , b->get_active_power( 0 , 1 ) ) , 2 ) ,
	  "hydro set_active_power_cost" + w + ", second change: "
	  "coefficient " + str( obj_coef( b , b->get_active_power( 0 , 1 ) ) )
	  + " instead of 2" );
   }
  catch( std::exception & e ) {
   check( false , "hydro set_active_power_cost" + w + " without the "
	  "datum: " + e.what() );
   }
  delete b;
  }
 }

/*--------------------------------------------------------------------------*/
/* A volume bound given alone is checked against the other one, absent and
 * hence 0: a positive MinVolumetric without MaxVolumetric (an empty
 * feasible set) and a negative MaxVolumetric are refused, a zero
 * MinVolumetric alone is accepted. */

static void test_hydro_volume_bound_alone( void )
{
 auto make = []( const std::string & name , double v ) -> HydroUnitBlock * {
  auto d = turbine_HU( 2 );
  d.maxV.clear();
  auto g = new_group( "HU" , true );
  write_HU( g , d );
  put( g , name , { g.getDim( "NumberReservoirs" ) ,
		    g.getDim( "TimeHorizon" ) } , { v , v } );
  try {
   return( dynamic_cast< HydroUnitBlock * >( Block::new_Block( g ) ) );
   }
  catch( std::exception & ) {
   return( nullptr );
   }
  };

 auto b = make( "MinVolumetric" , 3 );
 check( ! b , "hydro: MinVolumetric 3 without MaxVolumetric accepted" );
 delete b;
 b = make( "MaxVolumetric" , -1 );
 check( ! b , "hydro: MaxVolumetric -1 without MinVolumetric accepted" );
 delete b;
 b = make( "MinVolumetric" , 0 );
 check( b != nullptr , "hydro: MinVolumetric 0 without MaxVolumetric "
	"refused" );
 delete b;
 }

/*--------------------------------------------------------------------------*/
/* The overloads that set one cost value give it to every arc (hydro) or
 * every instant (battery, intermittent, slack), in the data and in the
 * Objective. */

static void test_single_value_costs( void )
{
 // hydro, two arcs
 {
  HUData d;
  d.T = 2;
  d.A = 2;
  d.start = { 0 , 0 };
  d.end = { 1 , 1 };
  d.maxF.assign( 4 , 10 );
  d.maxP.assign( 4 , 10 );
  d.lin = { 1 , 1 };
  d.maxV.assign( 2 , 100 );
  std::string why;
  if( auto b = new_HU( d , why ) ) {
   generate_all( b );
   b->set_active_power_cost( 5.0 );
   bool ok = true;
   for( Index l = 0 ; l < 2 ; ++l )
    for( Index t = 0 ; t < 2 ; ++t )
     ok = ok && close( obj_coef( b , b->get_active_power( l , t ) ) , 5 );
   check( ok , "hydro set_active_power_cost( 5 ): not 5 on every arc" );
   delete b;
   }
  else
   check( false , "hydro, two arcs: throws " + why );
  }

 // battery
 if( auto bu = new_BU( 4 , { { "MinStorage" , { 0 , 0 , 0 , 0 } } ,
			     { "MaxStorage" , { 9 , 9 , 9 , 9 } } ,
			     { "MaxPower" , { 5 , 5 , 5 , 5 } } } ) ) {
  generate_all( bu );
  bu->set_cost( 3.0 );
  check( bu->get_cost() == std::vector< double >( 4 , 3 ) ,
	 "battery set_cost( 3 ): the cost is not 3 at every instant" );
  auto p = bu->get_active_power( 0 );
  check( close( obj_coef( bu , & p[ 3 ] ) , -3 ) ,
	 "battery set_cost( 3 ): coefficient at 3 is " +
	 str( obj_coef( bu , & p[ 3 ] ) ) + " instead of -3" );
  delete bu;
  }
 else
  check( false , "battery for set_cost( double ): not read" );

 // intermittent
 {
  auto g = new_group( "IU" , true );
  g.putAtt( "type" , "IntermittentUnitBlock" );
  auto T = g.addDim( "TimeHorizon" , 3 );
  put( g , "MaxPower" , T , { 7 , 3 , 5 } );
  auto ib = dynamic_cast< IntermittentUnitBlock * >( Block::new_Block( g ) );
  if( ib ) {
   generate_all( ib );
   ib->set_active_power_cost( 4.0 );
   check( ib->get_active_power_cost() == std::vector< double >( 3 , 4 ) ,
	  "intermittent set_active_power_cost( 4 ): not 4 at every instant" );
   delete ib;
   }
  else
   check( false , "intermittent for set_active_power_cost( double ): not "
	  "read" );
  }

 // slack
 {
  auto g = new_group( "SU" , true );
  g.putAtt( "type" , "SlackUnitBlock" );
  auto T = g.addDim( "TimeHorizon" , 3 );
  put( g , "MaxPower" , T , { 100 , 100 , 100 } );
  auto sb = dynamic_cast< SlackUnitBlock * >( Block::new_Block( g ) );
  if( sb ) {
   generate_all( sb );
   sb->set_active_power_cost( 6.0 );
   check( sb->get_active_power_cost() == std::vector< double >( 3 , 6 ) ,
	  "slack set_active_power_cost( 6 ): not 6 at every instant" );
   delete sb;
   }
  else
   check( false , "slack for set_active_power_cost( double ): not read" );
  }
 }

/*--------------------------------------------------------------------------*/
/* The bound of the design used by the node injection bounds of a battery is
 * that of the battery design, which alone bounds the power: 5 with a
 * battery design up to 5 and no converter design (the default 1 of the
 * latter does not count), also with a converter design up to 2, and 1 with
 * a converter design alone. */

static void test_battery_design_ub( void )
{
 const std::vector< std::pair< std::string , std::vector< double > > > base =
  { { "MinStorage" , { 0 , 0 } } , { "MaxStorage" , { 20 , 20 } } ,
    { "MaxPower" , { 10 , 10 } } };

 struct Case { bool batt; bool conv; double ub; };
 for( const auto & c : { Case{ true , false , 5 } , Case{ true , true , 5 } ,
			 Case{ false , true , 1 } } ) {
  auto data = base;
  if( c.batt ) {
   data.push_back( { "BatteryInvestmentCost" , { 1 } } );
   data.push_back( { "BatteryMaxCapacityDesign" , { 5 } } );
   }
  if( c.conv ) {
   data.push_back( { "ConverterInvestmentCost" , { 1 } } );
   data.push_back( { "ConverterMaxCapacityDesign" , { 2 } } );
   }
  std::string why;
  auto bu = new_BU( 2 , data , & why );
  if( ! bu ) {
   check( false , "battery design ub: not read: " + why );
   continue;
   }
  check( bu->get_design_ub() == c.ub ,
	 std::string( "battery design ub with " ) +
	 ( c.batt ? "battery " : "" ) + ( c.conv ? "converter " : "" ) +
	 "design: " + str( bu->get_design_ub() ) + " instead of " +
	 str( c.ub ) );
  delete bu;
  }
 }

/*--------------------------------------------------------------------------*/
/* A continuous battery design with BatteryMaxCapacityDesign = 0 and a
 * positive BatteryMinCapacityDesign has an empty domain and is refused;
 * with both bounds 0 it is accepted. */

static void test_battery_design_zero_max( void )
{
 const std::vector< std::pair< std::string , std::vector< double > > > base =
  { { "MinStorage" , { 0 , 0 } } , { "MaxStorage" , { 20 , 20 } } ,
    { "MaxPower" , { 10 , 10 } } , { "BatteryInvestmentCost" , { 1 } } ,
    { "BatteryMaxCapacityDesign" , { 0 } } };
 auto data = base;
 data.push_back( { "BatteryMinCapacityDesign" , { 1 } } );
 auto bu = new_BU( 2 , data );
 check( ! bu , "battery: design in [ 1 , 0 ] accepted" );
 delete bu;
 std::string why;
 bu = new_BU( 2 , base , & why );
 check( bu != nullptr , "battery: design in [ 0 , 0 ] refused: " + why );
 delete bu;
 }

/*--------------------------------------------------------------------------*/
/* get_kappa_linearization() of a battery with an infinite MaxStorage (and
 * MinStorage) and all the duals 0 is 0, not the NaN of 0 times an infinite
 * bound. */

static void test_battery_kappa_infinite_storage( void )
{
 auto bu = new_BU( 2 , { { "MinStorage" , { -INF , -INF } } ,
			 { "MaxStorage" , { INF , INF } } ,
			 { "MaxPower" , { 10 , 10 } } } );
 if( ! bu ) {
  check( false , "battery with infinite storage bounds: not read" );
  return;
  }
 generate_all( bu );
 const double l = bu->get_kappa_linearization();
 check( l == 0 , "battery kappa linearization, infinite storage bounds and "
	"zero duals: " + str( l ) + " instead of 0" );
 delete bu;
 }

/* Tests of fvfix-NET, to be put in wt-ucb/test/test.cpp next to
 * test_ac_flow_group() (they use write_net(), new_net(), NetLine, put*(),
 * new_group(), check(), close(), str()).
 *
 * Extra #include (after "DCNetworkBlock.h"):
 *
 *   #include "DesignNetworkBlock.h"
 *
 * Call lines for main(), after test_ac_flow_group();
 *
 *   test_ac_line_losses_virtual();
 *   test_design_missing_subnetworks();
 */

/*--------------------------------------------------------------------------*/
/* The losses of an ACNetworkBlock read through a DCNetworkBlock pointer, as
 * UCBlock and InvestmentFunction do: the AC get_line_losses() overrides the
 * DC one (all 0), and gives for each line the sum of its "from" and "to"
 * flows; an HVDC line (all impedances 0) is among them. The stub
 * recover_feasible_solution() returns an empty vector, also with values on
 * the flows. */

static void test_ac_line_losses_virtual( void )
{
 auto g = new_group( "ACL" , true );
 g.putAtt( "type" , "ACNetworkBlock" );
 auto N = g.addDim( "NumberNodes" , 3 );
 auto L = g.addDim( "NumberLines" , 2 );
 put_int( g , "StartLine" , L , { 0 , 1 } );
 put_int( g , "EndLine" , L , { 1 , 2 } );
 put( g , "MinPowerFlow" , L , { -5 , -6 } );
 put( g , "MaxPowerFlow" , L , { 5 , 6 } );
 put( g , "LineSusceptance" , L , { 0 , 2 } );
 put( g , "ActiveDemand" , N , { 1 , 2 , 3 } );
 put( g , "LineReactance" , L , { 0 , 0.2 } );
 put( g , "LineResistance" , L , { 0 , 0.02 } );
 put( g , "NodeMaxVoltage" , N , { 1.1 , 1.1 , 1.1 } );
 put( g , "NodeMinVoltage" , N , { 0.9 , 0.9 , 0.9 } );
 put( g , "ReactiveDemand" , N , { 0.5 , 0.5 , 0.5 } );
 auto ac = dynamic_cast< ACNetworkBlock * >( Block::new_Block( g ) );
 check( ac , "AC losses: not read" );
 if( ! ac )
  return;
 ac->generate_abstract_variables();
 auto f = ac->get_static_variable_v< ColVariable >( "v_power_flow_real" );
 check( f && ( f->size() == 4 ) , "AC losses: no 4 flows" );
 if( ! ( f && ( f->size() == 4 ) ) ) {
  delete ac;
  return;
  }

 // from flows of lines 0, 1, then their to flows
 const std::vector< double > val = { 3 , -2 , -3 , 2.5 };
 for( Index i = 0 ; i < 4 ; ++i )
  ( *f )[ i ].set_value( val[ i ] );

 const DCNetworkBlock * dc = ac;
 const auto loss = dc->get_line_losses();
 check( loss.size() == 2 , "AC losses: " + std::to_string( loss.size() ) +
        " values through DCNetworkBlock instead of 2" );
 if( loss.size() == 2 ) {
  check( close( loss[ 0 ] , 0 ) , "AC losses: HVDC line " +
         str( loss[ 0 ] ) + " instead of 0" );
  check( close( loss[ 1 ] , 0.5 ) , "AC losses: AC line " +
         str( loss[ 1 ] ) + " instead of 0.5 (the DC version gives 0)" );
  }

 check( ac->recover_feasible_solution().empty() ,
        "AC recover_feasible_solution(): not empty" );

 delete ac;
 }

/*--------------------------------------------------------------------------*/
/* A DesignNetworkBlock with "NumberSubNetwork" = 3 and only the group
 * "NetworkBlock_0" (and then with none): the sub-Block 1 and 2 (all of
 * them, in the second case) are DCNetworkBlock with the NetworkData of the
 * DesignNetworkBlock, instead of being read out of range. */

static void test_design_missing_subnetworks( void )
{
 // the NetworkData of a two-node network with one line
 auto base = new_net( write_net( "DCNetworkBlock" , 2 , { { 0 , 1 , 1 ,
					      -5 , 5 } } , { 1 , 2 } ) , 2 );
 auto nd = base->get_NetworkData();

 for( bool first : { true , false } ) {
  const std::string who = std::string( "Design, " ) +
   ( first ? "only NetworkBlock_0" : "no NetworkBlock_<i>" );
  auto g = new_group( "DNB" , true );
  g.putAtt( "type" , "DesignNetworkBlock" );
  auto D = g.addDim( "NumberDesignLines" , 1 );
  g.addDim( "NumberSubNetwork" , 3 );
  put( g , "InvestmentCost" , D , { 10 } );
  put_int( g , "DesignLines" , D , { 0 } );
  if( first ) {
   auto s = g.addGroup( "NetworkBlock_0" );
   s.putAtt( "type" , "DCNetworkBlock" );
   auto SN = s.addDim( "NumberNodesD" , 2 );
   put( s , "ActiveDemand" , SN , { 1 , 2 } );
   }

  auto d = new DesignNetworkBlock();
  d->set_NetworkData( nd );
  try {
   d->deserialize( g );
   const auto & sub = d->get_nested_Blocks();
   check( sub.size() == 3 , who + ": " + std::to_string( sub.size() ) +
	  " sub-Block instead of 3" );
   for( Index n = 0 ; n < sub.size() ; ++n ) {
    auto nb = dynamic_cast< DCNetworkBlock * >( sub[ n ] );
    check( nb && ( nb->get_NetworkData() == nd ) , who + ": sub-Block " +
	   std::to_string( n ) + " not a DCNetworkBlock with the "
	   "NetworkData of the DesignNetworkBlock" );
    }
   }
  catch( std::exception & e ) {
   check( false , who + ": " + e.what() );
   }
  for( auto b : d->get_nested_Blocks() )
   delete b;
  delete d;
  }

 delete base;
 }

/*--------------------------------------------------------------------------*/
/* The bound psi of the run from before the horizon keeps the shut-down
 * term, as every other run: on since 5 instants at 40, with ShutDownLimit
 * 20 and a ramp down of 30, the run (0,2) (on at 0 and 1, off at 2) is at
 * most 20 + 30 = 50 at 0, which is the coefficient of its arc in the row
 * (27) of the pt formulation at 0 (100 without the term); the optimum, on
 * at 100 at every instant, is the same in every formulation. */

static void test_thermal_psi_initial_run( void )
{
 auto make = []( void ) {
  return( new_unit_gen( false , 3 ,
                        { { "MinPower" , { 10 , 10 , 10 } } ,
                          { "MaxPower" , { 100 , 100 , 100 } } ,
                          { "StartUpLimit" , { 100 , 100 , 100 } } ,
                          { "ShutDownLimit" , { 20 , 20 , 20 } } ,
                          { "DeltaRampUp" , { 100 , 100 , 100 } } ,
                          { "DeltaRampDown" , { 30 , 30 , 30 } } ,
                          { "LinearTerm" , { -10 , -10 , -10 } } } ,
                        { { "InitialPower" , 40 } } ,
                        { { "InitUpDownTime" , 5 } } ,
                        { { "MinUpTime" , 1 } , { "MinDownTime" , 1 } } ) );
  };
 auto tub = make();
 generate_all_wf( tub , ThermalUnitBlock::ptForm );
 auto rows = tub->get_static_constraint_v< FRowConstraint >(
                                                  "MaxPower_Const_Thermal" );
 std::vector< double > cf;
 if( rows && ( ! rows->empty() ) )
  if( auto lf = dynamic_cast< const LinearFunction * >(
                                       ( *rows )[ 0 ].get_function() ) )
   for( Index k = 0 ; k < lf->get_num_active_var() ; ++k )
    cf.push_back( lf->get_coefficient( k ) );
 std::string got;
 for( auto c : cf )
  got += " " + str( c );
 check( std::count( cf.begin() , cf.end() , 50.0 ) == 1 ,
        "psi of the initial run: the row (27) at 0 has" + got +
        ", no coefficient 50 of the arc ( 0 , 2 )" );
 delete tub;
 check_all_forms( make , -3000 , "psi of the initial run" );
 }

/*--------------------------------------------------------------------------*/
/// the Configuration in the file name of the directory of the test files
/// (UCBLOCK_TEST_DIR), nullptr if it cannot be read

static Configuration * read_test_cfg( const std::string & name )
{
 auto prefix = Configuration::get_filename_prefix();
 Configuration::set_filename_prefix( std::string( UCBLOCK_TEST_DIR ) + "/" );
 Configuration * c = nullptr;
 try {
  c = Configuration::deserialize( name );
  }
 catch( std::exception & e ) {
  std::cout << "cannot read " << name << ": " << e.what() << std::endl;
  }
 Configuration::set_filename_prefix( std::move( prefix ) );
 return( c );
 }

/*--------------------------------------------------------------------------*/
/// the unit generated in the formulation that the BlockConfig in the file
/// form selects [see TUBCfg-*.txt]; false if the file cannot be read

static bool generate_from_file( ThermalUnitBlock * tub ,
                                const std::string & form )
{
 auto bc = dynamic_cast< BlockConfig * >( read_test_cfg( form ) );
 if( ! bc )
  return( false );
 bc->apply( tub );
 delete bc;
 tub->generate_abstract_variables();
 tub->generate_abstract_constraints();
 tub->generate_objective();
 return( true );
 }

/*--------------------------------------------------------------------------*/
/// the value of the continuous relaxation of the unit, given by the
/// :MILPSolver of LPRelaxBSCfg.txt: NaN if that Solver is not in the build
/// or finds no value, INF if the relaxation is infeasible

static double relaxation_value( ThermalUnitBlock * tub )
{
 auto bsc = dynamic_cast< BlockSolverConfig * >(
                                     read_test_cfg( "LPRelaxBSCfg.txt" ) );
 if( ! bsc )
  return( std::nan( "" ) );
 for( const auto & s : bsc->get_SolverNames() )
  if( ! Solver::has_Solver( s ) ) {
   delete bsc;
   return( std::nan( "" ) );
   }
 bsc->apply( tub );
 bsc->clear();
 double v = std::nan( "" );
 if( ! tub->get_registered_solvers().empty() ) {
  auto slv = tub->get_registered_solvers().front();
  try {
   const auto status = slv->compute();
   if( status == Solver::kOK )
    v = slv->get_var_value();
   else
    if( status == Solver::kInfeasible )
     v = INF;
   }
  catch( std::exception & e ) {
   std::cout << "the relaxation throws " << e.what() << std::endl;
   }
  }
 bsc->apply( tub );
 delete bsc;
 return( v );
 }

/*--------------------------------------------------------------------------*/
/* The SUSD formulation has the ramp rows of the SU and SD formulations,
 * i.e., also those downwards on the powers p^h of the runs that start at h
 * and upwards on the powers of the runs that end at k, besides those over
 * several steps. Without them the continuous relaxation is weaker than the
 * DP one: over T = 4 with MinPower 10, MaxPower 100, StartUpLimit and
 * ShutDownLimit 40, both ramps 30 and the unit off since long, the half of
 * the run that starts at 0 and stays on to the end, aggregated with the
 * half that shuts down at 3, satisfies the other rows while the former
 * goes from 40 to 100 from 2 to 3, so that the relaxation is -90 instead
 * of -76, the integer optimum, which is also the value of the relaxation
 * of the DP formulation. The ramp groups are not empty in the SUSD
 * formulation, and with the unit on before the horizon (the rows at 0)
 * set_initial_power() updates them: the relaxation is then that of the
 * unit generated with the new initial power. */

static void test_thermal_SUSD_ramp_rows( void )
{
 auto make = []( double ip , int init ) {
  return( new_unit_gen( false , 4 ,
                        { { "MinPower" , { 10 , 10 , 10 , 10 } } ,
                          { "MaxPower" , { 100 , 100 , 100 , 100 } } ,
                          { "StartUpLimit" , { 40 , 40 , 40 , 40 } } ,
                          { "ShutDownLimit" , { 40 , 40 , 40 , 40 } } ,
                          { "DeltaRampUp" , { 30 , 30 , 30 , 30 } } ,
                          { "DeltaRampDown" , { 30 , 30 , 30 , 30 } } ,
                          { "QuadTerm" , { 0 , 0 , 0 , 0 } } ,
                          { "LinearTerm" , { -5 , -5 , 2 , -1 } } ,
                          { "ConstTerm" , { 68 , 68 , 68 , 68 } } ,
                          { "StartUpCost" , { 192 , 192 , 192 , 192 } } } ,
                        { { "InitialPower" , ip } } ,
                        { { "InitUpDownTime" , init } } ,
                        { { "MinUpTime" , 1 } , { "MinDownTime" , 1 } } ) );
  };

 // the integer optimum, every formulation and both DP Solvers
 check_all_forms( [ & ]() { return( make( 0 , -5 ) ); } , -76 ,
                  "SUSD ramp rows" );

 // the continuous relaxations of the SUSD and DP formulations
 for( const std::string form : { "TUBCfg-SUSD.txt" , "TUBCfg-DP.txt" } ) {
  auto tub = make( 0 , -5 );
  if( ! generate_from_file( tub , form ) ) {
   check( false , "SUSD ramp rows: cannot read " + form );
   delete tub;
   continue;
   }
  if( form == "TUBCfg-SUSD.txt" )
   for( const auto & g : { "RampUp_Const_Thermal" ,
                           "RampDown_Const_Thermal" } ) {
    auto rows = tub->get_static_constraint_v< FRowConstraint >( g );
    check( rows && ( ! rows->empty() ) ,
           std::string( "SUSD ramp rows: the group " ) + g + " is empty" );
    }
  const auto v = relaxation_value( tub );
  if( std::isnan( v ) )
   std::cout << "SUSD ramp rows, " << form << ": no relaxation value (no "
             << ":MILPSolver of LPRelaxBSCfg.txt?), skipped" << std::endl;
  else
   check( close( v , -76 ) , "SUSD ramp rows, " + form +
          ": the continuous relaxation is " + str( v ) +
          " instead of -76" );
  delete tub;
  }

 // on before the horizon at 40, then at 50 (above ShutDownLimit[ 0 ])
 auto tub = make( 40 , 5 );
 generate_from_file( tub , "TUBCfg-SUSD.txt" );
 std::vector< double > ip = { 50 };
 try {
  tub->set_initial_power( ip.cbegin() , Range( 0 , 1 ) );
  auto ref = make( 50 , 5 );
  generate_from_file( ref , "TUBCfg-SUSD.txt" );
  const auto v = relaxation_value( tub );
  const auto r = relaxation_value( ref );
  check( ( std::isnan( v ) && std::isnan( r ) ) || close( v , r ) ,
         "SUSD ramp rows, initial power 40 -> 50: the relaxation is " +
         str( v ) + " instead of " + str( r ) );
  delete ref;
  }
 catch( std::exception & e ) {
  check( false , std::string( "SUSD ramp rows, set_initial_power( 50 ) "
                              "throws " ) + e.what() );
  }
 delete tub;
 }

/*--------------------------------------------------------------------------*/
/* The bound psi of the maximum power rows (27) of the pt, SU, SD and SUSD
 * formulations holds at every instant of a run, the first and the last one
 * included, where the start-up and shut-down limits alone were used. With
 * MinPower 10 and MaxPower 100 at every instant, in the cases
 *
 * - A: T = 2, StartUpLimit 40, ShutDownLimit 30, both ramps 10, minimum up
 *   and down times 1, on since 3 instants at 20: the run from before the
 *   horizon that the end of the horizon cuts is at most 20 + 2 * 10 = 40
 *   at 1, both as the arc ( 0 , T ) and as the arc ( 0 , T + 1 ), which
 *   describe the same schedule (the former had 100, the maximum power);
 *
 * - B: T = 3, StartUpLimit 50, ShutDownLimit 30, DeltaRampUp 40 and no
 *   DeltaRampDown, minimum up and down times 1, off since 4 instants;
 *
 * - C: T = 3, StartUpLimit = ShutDownLimit = 50, DeltaRampUp 10,
 *   DeltaRampDown 30, minimum up time 2, on since 3 instants at 20: the
 *   run from before the horizon that shuts down at 2 is at most
 *   20 + 2 * 10 = 40 at its last instant 1, below the shut-down limit;
 *
 * - E: T = 3, StartUpLimit 40, ShutDownLimit 30, no ramps, on since 2
 *   instants at 60; F: T = 3, StartUpLimit 30, ShutDownLimit 40, only
 *   DeltaRampDown 20, minimum up time 2, off since 1 instant;
 *
 * with the costs below, the continuous relaxation of every formulation but
 * the 3bin and T ones is the integer optimum, which also the MILP of every
 * formulation and both DP Solvers give; before, those of the pt and SU
 * formulations were -137 in A and -134 in C, those of the pt and SD ones
 * -108 in B. The values come from an independent model of the rows. In A,
 * set_initial_power() to 30 gives the relaxation of the unit generated with
 * 30. A NuclearUnitBlock keeps the rows, which its operating rules do not
 * make invalid: A and C as nuclear units (modulation ramps 5, a modulation
 * every 2 instants, C also with the bands 30 and 60 and modulations of at
 * most 2 instants) and random small nuclear units have the same MILP value
 * in every formulation as the NuclearUnitExtDPSolver, and as the brute
 * force without bands. */

static void test_thermal_psi_first_last( void )
{
 struct Case {
  std::string name; Index T; double su , sd , ru , rd; unsigned int up;
  int init; double ip; std::vector< double > lin , cnst , suc;
  double opt;
  };
 const std::vector< Case > cases = {
  { "A" , 2 , 40 , 30 , 10 , 10 , 1 , 3 , 20 , { -5 , -2 } , { 51 , 48 } ,
    { 8 , 151 } , -131 } ,
  { "B" , 3 , 50 , 30 , 40 , -1 , 1 , -4 , 0 , { 1 , -1 , -3 } ,
    { 52 , 72 , 81 } , { 102 , 64 , 177 } , -103 } ,
  { "C" , 3 , 50 , 50 , 10 , 30 , 2 , 3 , 20 , { -4 , -2 , 2 } ,
    { 18 , 52 , 80 } , { 61 , 147 , 151 } , -130 } ,
  { "E" , 3 , 40 , 30 , -1 , -1 , 1 , 2 , 60 , { 0 , 1 , -6 } ,
    { 80 , 46 , 51 } , { 126 , 57 , 195 } , -413 } ,
  { "F" , 3 , 30 , 40 , -1 , 20 , 2 , -1 , 0 , { -3 , -5 , -2 } ,
    { 86 , 16 , 80 } , { 92 , 110 , 18 } , -516 } };

 auto make = []( const Case & c , double ip ) {
  auto cst = [ & ]( double v ) { return( std::vector< double >( c.T , v ) ); };
  std::vector< std::pair< std::string , std::vector< double > > > vecs = {
   { "MinPower" , cst( 10 ) } , { "MaxPower" , cst( 100 ) } ,
   { "StartUpLimit" , cst( c.su ) } , { "ShutDownLimit" , cst( c.sd ) } ,
   { "QuadTerm" , cst( 0 ) } , { "LinearTerm" , c.lin } ,
   { "ConstTerm" , c.cnst } , { "StartUpCost" , c.suc } };
  if( c.ru >= 0 )
   vecs.push_back( { "DeltaRampUp" , cst( c.ru ) } );
  if( c.rd >= 0 )
   vecs.push_back( { "DeltaRampDown" , cst( c.rd ) } );
  std::vector< std::pair< std::string , double > > dbls;
  if( c.init > 0 )
   dbls.push_back( { "InitialPower" , ip } );
  return( new_unit_gen( false , c.T , vecs , dbls ,
                        { { "InitUpDownTime" , c.init } } ,
                        { { "MinUpTime" , c.up } , { "MinDownTime" , 1 } } ) );
  };

 const std::vector< std::string > forms = { "TUBCfg-pt.txt" ,
  "TUBCfg-DP.txt" , "TUBCfg-SU.txt" , "TUBCfg-SD.txt" , "TUBCfg-SUSD.txt" };

 for( const auto & c : cases ) {
  const auto who = "psi at the first and last instant, case " + c.name;

  // the integer optimum, every formulation and both DP Solvers
  check_all_forms( [ & ]() { return( make( c , c.ip ) ); } , c.opt , who );

  // the continuous relaxations
  for( const auto & form : forms ) {
   auto tub = make( c , c.ip );
   if( ! generate_from_file( tub , form ) ) {
    check( false , who + ": cannot read " + form );
    delete tub;
    continue;
    }
   const auto v = relaxation_value( tub );
   if( std::isnan( v ) )
    std::cout << who << ", " << form << ": no relaxation value (no "
              << ":MILPSolver of LPRelaxBSCfg.txt?), skipped" << std::endl;
   else
    check( close( v , c.opt ) , who + ", " + form + ": the continuous "
           "relaxation is " + str( v ) + " instead of " + str( c.opt ) );
   delete tub;
   }
  }

 // A and C as nuclear units, and random small nuclear units
 auto nuclear = []( const Case & c ) {
  TUData d;
  d.nuclear = true;
  d.T = c.T;
  d.minP = 10;
  d.maxP = 100;
  d.su = c.su;
  d.sd = c.sd;
  d.ru = c.ru;
  d.rd = c.rd;
  d.initUD = c.init;
  d.initP = c.ip;
  d.minUp = c.up;
  d.lin = c.lin;
  d.cnst = c.cnst;
  d.suc = c.suc;
  d.modT = 2;
  d.initMod = 2;
  d.mru = d.mrd = 5;
  return( d );
  };
 for( Index i : { 0 , 2 } ) {
  const auto d = nuclear( cases[ i ] );
  const auto who = "psi at the first and last instant, nuclear case " +
                   cases[ i ].name;
  const double bf = brute_force( d , std::vector< int >( d.T , -1 ) );
  check_all_forms( [ & ]() { return( new_TU( d ) ); } , bf , who );
  auto nb = [ & ]() { return( new_band_NU( d , 2 , { 30 , 60 } ) ); };
  auto nub = nb();
  const double v = nuclear_DP_value( nub );
  delete nub;
  check_all_forms( nb , v , who + " with bands" );
  }
 for( int i = 0 ; i < 30 ; ++i ) {
  TUData d;
  d.nuclear = true;
  d.T = rnd( 2 , 4 );
  d.minP = 10;
  d.maxP = 100;
  d.su = rnd( 2 , 8 ) * 10;
  d.sd = rnd( 2 , 8 ) * 10;
  d.ru = rnd( 1 , 5 ) * 10;
  d.rd = rnd( 1 , 5 ) * 10;
  d.minUp = rnd( 1 , 3 );
  d.minDown = rnd( 1 , 2 );
  d.initUD = rnd( 0 , 1 ) ? rnd( 1 , 4 ) : -rnd( 1 , 4 );
  d.initP = rnd( 1 , 10 ) * 10;
  d.modT = rnd( 2 , 4 );
  d.initMod = rnd( 1 , int( d.modT ) );
  d.mru = rnd( 0 , int( d.ru ) );
  d.mrd = rnd( 0 , int( d.rd ) );
  for( Index t = 0 ; t < d.T ; ++t ) {
   d.lin.push_back( rnd( -6 , 2 ) );
   d.cnst.push_back( rnd( 0 , 99 ) );
   d.suc.push_back( rnd( 0 , 199 ) );
   }
  const double bf = brute_force( d , std::vector< int >( d.T , -1 ) );
  check_all_forms( [ & ]() { return( new_TU( d ) ); } , bf ,
                   "psi at the first and last instant, random nuclear " +
                   describe( d ) , { 2 , 3 , 4 , 5 , 6 } );
  }

 // A: the row (27) of the pt formulation at the last instant caps every
 // run on at it by 40, the arcs ( 0 , T ) and ( 0 , T + 1 ) included
 {
  auto tub = make( cases[ 0 ] , 20 );
  generate_from_file( tub , "TUBCfg-pt.txt" );
  auto rows = tub->get_static_constraint_v< FRowConstraint >(
                                                  "MaxPower_Const_Thermal" );
  std::string got;
  bool ok = rows && ( rows->size() == 2 );
  if( ok )
   if( auto lf = dynamic_cast< const LinearFunction * >(
                                       ( *rows )[ 1 ].get_function() ) )
    for( Index k = 0 ; k < lf->get_num_active_var() ; ++k ) {
     const auto cf = lf->get_coefficient( k );
     got += " " + str( cf );
     if( ( cf != -1 ) && ( cf != 40 ) )
      ok = false;
     }
  check( ok , "psi at the last instant, case A: the row (27) at 1 has" +
         got + " instead of 40 for every run" );
  delete tub;
  }

 // A: set_initial_power() to 30, the relaxation of the unit generated
 // with 30, in the pt and SU formulations
 for( const std::string form : { "TUBCfg-pt.txt" , "TUBCfg-SU.txt" } ) {
  auto tub = make( cases[ 0 ] , 20 );
  auto ref = make( cases[ 0 ] , 30 );
  generate_from_file( tub , form );
  generate_from_file( ref , form );
  std::vector< double > ip = { 30 };
  try {
   tub->set_initial_power( ip.cbegin() , Range( 0 , 1 ) );
   const auto v = relaxation_value( tub );
   const auto r = relaxation_value( ref );
   check( ( std::isnan( v ) && std::isnan( r ) ) || close( v , r ) ,
          "psi at the first and last instant, case A, " + form +
          ", initial power 20 -> 30: the relaxation is " + str( v ) +
          " instead of " + str( r ) );
   }
  catch( std::exception & e ) {
   check( false , "psi at the first and last instant, case A, " + form +
          ": set_initial_power( 30 ) throws " + e.what() );
   }
  delete ref;
  delete tub;
  }
 }

/*--------------------------------------------------------------------------*/

int main( int argc , char ** argv )
{
 // --power-limits: only the test of the power limits, with the options
 // that follow [see test_power_limits()]
 if( ( argc > 1 ) && ( std::string( argv[ 1 ] ) == "--power-limits" ) ) {
  try {
   failures = test_power_limits( UCBLOCK_TEST_DIR ,
                                 std::vector< std::string >( argv + 2 ,
                                                             argv + argc ) );
   }
  catch( std::exception & e ) {
   std::cout << "uncaught exception: " << e.what() << std::endl;
   ++failures;
   }
  std::cout << ( failures ? std::to_string( failures ) +
                            " case(s) FAILED" : "All tests passed!!" )
            << std::endl;
  return( failures ? 1 : 0 );
  }

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
  test_thermal_3bin_max_power();
  test_thermal_data_checks();
  test_thermal_DP_one_instant_run();
  test_thermal_initial_power_change();
  test_nuclear_band_horizon();
  test_nuclear_band_initial_power();
  test_nuclear_band_rule();
  test_thermal_reserve_cost();
  test_thermal_reserve_t0();
  test_thermal_DP_two_rewarded_reserves();
  test_thermal_DP_availability();
  test_thermal_power_limits_change();
  test_thermal_reactive_bounds();
  test_thermal_short_horizon_min_up_down();
  test_thermal_SUSD_varying_bounds();
  test_thermal_zero_ramp();
  test_thermal_initial_power_no_ramps();
  test_thermal_shut_down_at_zero();
  test_thermal_initial_power_above_max();
  test_thermal_initial_power_above_max_change();
  test_thermal_unordered_subset();
  test_thermal_cost_order();
  test_thermal_DP_as_rows();
  test_thermal_solution_restore();
  test_nuclear_modulation_ramps_change();
  test_nuclear_DP_fixing_unlocks();
  test_thermal_update_rows_groups();
  test_thermal_init_updown_change();
  failures += test_power_limits( UCBLOCK_TEST_DIR , { } );

  test_RT_thermal();
  test_RT_hydro();
  test_RT_other_units();
  test_RT_networks();
  test_cycle_basis();
  test_net_components();
  test_net_reference();
  test_net_ptdf();
  test_net_demand_change();
  test_net_OTS_cost();
  test_net_negative_cost();
  test_net_dual_prices();
  test_RT_UCBlock();

  test_setters_thermal( 1 );
  test_setters_thermal( 3 );
  test_setters_UCBlock();

  test_hydro_delays();
  test_hydro_initial_volume();
  test_hydro_cyclic();
  test_hydro_data_checks();
  test_hydro_spillage();
  test_hydro_inertia();
  test_hydro_system();

  test_battery_storage_balance();
  test_battery_initial_storage_sign();
  test_battery_set_cost_unsorted();
  test_battery_data_checks();
  test_battery_converter();
  test_kappa_linearization_design();
  test_intermittent_design_reserve();
  test_design_fixed_to_one();
  test_battery_two_kappas();
  test_battery_design_modules();
  test_battery_converter_modules();
  test_slack_unit();

  test_zones_UCBlock();
  test_slack_no_inertia_UCBlock();
  test_network_constants_UCBlock();
  test_dual_signs_UCBlock();
  test_node_injection_bounds();
  test_zone_rows_two_nodes();
  test_storage_infinite_bounds();
  test_ac_flow_group();
  test_ac_line_losses_virtual();
  test_design_missing_subnetworks();

  test_thermal_DP_cost_change();
  test_thermal_DP_extended_fixing();
  test_thermal_initial_power_refused();
  test_thermal_is_feasible_fix_to_max();
  test_nuclear_rule_costs_unordered();
  test_reactive_injection_bounds_on();
  test_ac_no_angle_bounds();
  test_ac_angle_bounds_matpower();
  test_ac_angle_bounds_unbounded();
  test_ac_flow_check_no_names();
  test_thermal_T_ramp_bound_rows();
  test_thermal_shut_down_limit_at_t0();
  test_thermal_psi_initial_run();
  test_thermal_SUSD_ramp_rows();
  test_thermal_psi_first_last();

  test_hydro_set_inflow_subset();
  test_hydro_initial_flow_absent();
  test_hydro_cost_without_datum();
  test_hydro_volume_bound_alone();
  test_single_value_costs();
  test_battery_design_ub();
  test_battery_design_zero_max();
  test_battery_kappa_infinite_storage();
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
