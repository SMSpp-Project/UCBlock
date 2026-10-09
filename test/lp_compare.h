/*--------------------------------------------------------------------------*/
/*--------------------------- File lp_compare.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Reading and comparing two models written in the LP format by a
 * :MILPSolver (strOutputFile), as needed to check that the abstract
 * representation of a Block changed in place equals the one generated
 * afresh from the changed data.
 *
 * A model is read into a canonical form in which nothing depends on the
 * order of the rows or of the terms: zero coefficients are dropped, the
 * terms of a row and of the Objective are sorted by the name of the
 * column, a row is a pair of sides [ lo , hi ] (one of them infinite for an
 * inequality, equal for an equality), and the slack column that a solver
 * adds for a ranged row (named "Rg" followed by the name of the row) is
 * folded back into the sides of its row. Columns are matched by their
 * name, which a :MILPSolver builds from the name of the group of the
 * Variable and the position in it (e.g., "p_thermal_0_3"), hence by role;
 * rows are matched by content, as a multiset, the row name only serving to
 * attribute a row to its group (the name stripped of the trailing
 * "_<number>" fields) and to pair a row that has no match with the row of
 * the same name in the other model, so that the report shows which
 * coefficients or sides differ.
 *
 * The header needs nothing but the standard library. */
/*--------------------------------------------------------------------------*/

#ifndef __LP_COMPARE_H
#define __LP_COMPARE_H

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

/*--------------------------------------------------------------------------*/

namespace lp_compare
{

/*--------------------------------------------------------------------------*/
/*------------------------------- TYPES ------------------------------------*/
/*--------------------------------------------------------------------------*/

using Terms = std::vector< std::pair< std::string , double > >;

/// a row lo <= sum terms <= hi, with sorted terms and no zero coefficient
struct Row {
 std::string name;   ///< the name in the file
 std::string group;  ///< the name without the trailing "_<number>" fields
 double lo = - std::numeric_limits< double >::infinity();
 double hi = std::numeric_limits< double >::infinity();
 Terms terms;
 };

/// a column: bounds and integrality
struct Column {
 double lo = 0;
 double hi = std::numeric_limits< double >::infinity();
 bool integer = false;
 };

/// the canonical form of a model
struct Model {
 bool minimize = true;
 std::map< std::string , double > obj;    ///< linear terms
 std::map< std::string , double > qobj;   ///< "x*y" (sorted) -> coefficient
 double obj_const = 0;
 std::vector< Row > rows;
 std::map< std::string , Column > cols;
 };

/// the differences in the rows of one group
struct GroupDiff {
 std::size_t in_a = 0;      ///< rows of the group in A
 std::size_t in_b = 0;      ///< rows of the group in B
 std::size_t only_a = 0;    ///< rows of A with no equal row in B
 std::size_t only_b = 0;    ///< rows of B with no equal row in A
 std::vector< std::string > examples;
 };

/// all the differences between two models
struct Diff {
 std::map< std::string , GroupDiff > rows;        ///< by row group
 std::map< std::string , std::size_t > cols;      ///< by column group
 std::vector< std::string > col_examples;
 std::vector< std::string > obj_examples;
 std::size_t n_obj = 0;
 std::size_t skipped_a = 0;  ///< rows of A in the skipped groups
 std::size_t skipped_b = 0;  ///< rows of B in the skipped groups

 bool empty( void ) const {
  if( n_obj || ( ! cols.empty() ) )
   return( false );
  for( const auto & g : rows )
   if( g.second.only_a || g.second.only_b )
    return( false );
  return( true );
  }

 /// the groups (rows, then columns, then "objective") that differ, comma
 /// separated
 std::string groups( void ) const {
  std::string s;
  auto add = [ & ]( const std::string & g ) {
   if( ! s.empty() )
    s += ",";
   s += g;
   };
  for( const auto & g : rows )
   if( g.second.only_a || g.second.only_b )
    add( g.first );
  for( const auto & c : cols )
   add( "col:" + c.first );
  if( n_obj )
   add( "objective" );
  return( s );
  }

 /// a readable report, at most max_ex examples per group
 std::string report( std::size_t max_ex = 3 ) const {
  std::ostringstream s;
  for( const auto & g : rows ) {
   if( ! ( g.second.only_a || g.second.only_b ) )
    continue;
   s << "  rows " << g.first << ": A " << g.second.in_a << ", B "
     << g.second.in_b << ", only in A " << g.second.only_a
     << ", only in B " << g.second.only_b << "\n";
   for( std::size_t i = 0 ;
        i < std::min( max_ex , g.second.examples.size() ) ; ++i )
    s << "    " << g.second.examples[ i ] << "\n";
   }
  for( const auto & c : cols )
   s << "  columns " << c.first << ": " << c.second << " differ\n";
  for( std::size_t i = 0 ; i < std::min( 3 * max_ex , col_examples.size() ) ;
       ++i )
   s << "    " << col_examples[ i ] << "\n";
  if( n_obj ) {
   s << "  objective: " << n_obj << " terms differ\n";
   for( std::size_t i = 0 ; i < std::min( max_ex , obj_examples.size() ) ;
        ++i )
    s << "    " << obj_examples[ i ] << "\n";
   }
  return( s.str() );
  }
 };

/*--------------------------------------------------------------------------*/
/*----------------------------- NUMBERS ------------------------------------*/
/*--------------------------------------------------------------------------*/

inline constexpr double INF = std::numeric_limits< double >::infinity();

/// true if s is a number (possibly signed, possibly "inf" / "infinity"),
/// and then its value in v
inline bool to_number( const std::string & s , double & v )
{
 if( s.empty() )
  return( false );
 std::string l( s );
 std::transform( l.begin() , l.end() , l.begin() ,
                 []( unsigned char c ) { return( std::tolower( c ) ); } );
 std::size_t i = ( ( l[ 0 ] == '+' ) || ( l[ 0 ] == '-' ) ) ? 1 : 0;
 const std::string body = l.substr( i );
 if( ( body == "inf" ) || ( body == "infinity" ) ) {
  v = ( l[ 0 ] == '-' ) ? - INF : INF;
  return( true );
  }
 if( body.empty() ||
     ! ( std::isdigit( static_cast< unsigned char >( body[ 0 ] ) ) ||
         ( body[ 0 ] == '.' ) ) )
  return( false );
 char * end = nullptr;
 v = std::strtod( s.c_str() , & end );
 return( end && ( * end == 0 ) );
 }

/// a == b within the relative tolerance tol (infinities equal if equal)
inline bool same( double a , double b , double tol )
{
 if( std::isinf( a ) || std::isinf( b ) )
  return( a == b );
 return( std::abs( a - b ) <=
         tol * std::max( { 1.0 , std::abs( a ) , std::abs( b ) } ) );
 }

inline std::string num( double v )
{
 if( std::isinf( v ) )
  return( v > 0 ? "+inf" : "-inf" );
 std::ostringstream s;
 s.precision( 10 );
 s << v;
 return( s.str() );
 }

/*--------------------------------------------------------------------------*/
/*------------------------------ NAMES -------------------------------------*/
/*--------------------------------------------------------------------------*/

/// the group of a row or column name: the name without the trailing
/// "_<number>" fields, e.g., "MaxPower_Const_Thermal_0_5" ->
/// "MaxPower_Const_Thermal"; a name made only of such fields, or a generic
/// solver name such as "R12" or "C7", is its own group "(unnamed)"
inline std::string group_of( const std::string & name )
{
 std::string g = name;
 for( ; ; ) {
  const auto pos = g.find_last_of( '_' );
  if( ( pos == std::string::npos ) || ( pos + 1 == g.size() ) )
   break;
  const auto tail = g.substr( pos + 1 );
  if( ! std::all_of( tail.begin() , tail.end() , []( unsigned char c ) {
        return( std::isdigit( c ) ); } ) )
   break;
  g.erase( pos );
  }
 if( ( g.size() > 1 ) && ( g == name ) &&
     ( ( g[ 0 ] == 'R' ) || ( g[ 0 ] == 'C' ) || ( g[ 0 ] == 'x' ) ||
       ( g[ 0 ] == 'c' ) ) &&
     std::all_of( g.begin() + 1 , g.end() , []( unsigned char c ) {
      return( std::isdigit( c ) ); } ) )
  return( "(unnamed)" );
 return( g );
 }

/*--------------------------------------------------------------------------*/
/*------------------------------ READING -----------------------------------*/
/*--------------------------------------------------------------------------*/

namespace detail
{

enum Section { eNone , eObj , eRows , eBounds , eInt , eBin , eIgnore ,
               eEnd };

/// the section a line opens, or eNone if it opens none
inline Section section_of( const std::string & line , bool & maximize )
{
 std::string l;
 for( char c : line )
  l += char( std::tolower( static_cast< unsigned char >( c ) ) );
 const auto b = l.find_first_not_of( " \t\r" );
 if( b == std::string::npos )
  return( eNone );
 const auto e = l.find_last_not_of( " \t\r" );
 l = l.substr( b , e - b + 1 );

 if( ( l == "minimize" ) || ( l == "minimum" ) || ( l == "min" ) ) {
  maximize = false;
  return( eObj );
  }
 if( ( l == "maximize" ) || ( l == "maximum" ) || ( l == "max" ) ) {
  maximize = true;
  return( eObj );
  }
 if( ( l == "subject to" ) || ( l == "such that" ) || ( l == "st" ) ||
     ( l == "s.t." ) )
  return( eRows );
 if( ( l == "bounds" ) || ( l == "bound" ) )
  return( eBounds );
 if( ( l == "generals" ) || ( l == "general" ) || ( l == "gen" ) ||
     ( l == "integers" ) || ( l == "integer" ) )
  return( eInt );
 if( ( l == "binaries" ) || ( l == "binary" ) || ( l == "bin" ) )
  return( eBin );
 if( ( l == "semi-continuous" ) || ( l == "semis" ) || ( l == "semi" ) ||
     ( l == "sos" ) || ( l == "sos1" ) || ( l == "sos2" ) ||
     ( l == "lazy constraints" ) || ( l == "user cuts" ) ||
     ( l == "general constraints" ) || ( l == "pwlobj" ) )
  return( eIgnore );
 if( l == "end" )
  return( eEnd );
 return( eNone );
 }

/// splits a piece of the file in tokens: blanks separate, the characters
/// [ ] ^ * / : are tokens of their own, and so are the relational
/// operators <= >= = < > =< =>; a sign glued to a name is split from it
inline void tokenize( const std::string & s , std::vector< std::string > & t )
{
 std::string cur;
 auto flush = [ & ]() {
  if( cur.empty() )
   return;
  double v;
  if( ( ( cur[ 0 ] == '+' ) || ( cur[ 0 ] == '-' ) ) && ( cur.size() > 1 ) &&
      ! to_number( cur , v ) ) {
   t.push_back( cur.substr( 0 , 1 ) );
   t.push_back( cur.substr( 1 ) );
   }
  else
   t.push_back( cur );
  cur.clear();
  };

 for( std::size_t i = 0 ; i < s.size() ; ++i ) {
  const char c = s[ i ];
  if( std::isspace( static_cast< unsigned char >( c ) ) ) {
   flush();
   continue;
   }
  if( ( c == '[' ) || ( c == ']' ) || ( c == '^' ) || ( c == '*' ) ||
      ( c == '/' ) || ( c == ':' ) ) {
   flush();
   t.push_back( std::string( 1 , c ) );
   continue;
   }
  if( ( c == '<' ) || ( c == '>' ) || ( c == '=' ) ) {
   flush();
   std::string op( 1 , c );
   if( ( i + 1 < s.size() ) &&
       ( ( s[ i + 1 ] == '=' ) || ( s[ i + 1 ] == '<' ) ||
         ( s[ i + 1 ] == '>' ) ) )
    op += s[ ++i ];
   t.push_back( op );
   continue;
   }
  // a sign starts a new token unless it is the sign of an exponent
  if( ( ( c == '+' ) || ( c == '-' ) ) && ( ! cur.empty() ) &&
      ! ( ( ( cur.back() == 'e' ) || ( cur.back() == 'E' ) ) &&
          ( cur.size() > 1 ) &&
          std::isdigit( static_cast< unsigned char >( cur[ cur.size() - 2 ] ) )
          ) )
   flush();
  cur += c;
  }
 flush();
 }

inline bool is_op( const std::string & s )
{
 return( ( s == "<=" ) || ( s == ">=" ) || ( s == "=" ) || ( s == "<" ) ||
         ( s == ">" ) || ( s == "=<" ) || ( s == "=>" ) || ( s == "==" ) );
 }

/// -1 for <=, +1 for >=, 0 for =
inline int op_dir( const std::string & s )
{
 if( ( s == "<=" ) || ( s == "<" ) || ( s == "=<" ) )
  return( -1 );
 if( ( s == ">=" ) || ( s == ">" ) || ( s == "=>" ) )
  return( 1 );
 return( 0 );
 }

/// reads a linear expression (with an optional quadratic part in brackets)
/// from t[ i ] on, until an operator, a name followed by ":", or the end;
/// adds the terms to lin and qd, the constant to cst
inline void read_expr( const std::vector< std::string > & t , std::size_t & i ,
                       std::map< std::string , double > & lin ,
                       std::map< std::string , double > & qd , double & cst ,
                       bool stop_at_label )
{
 double sign = 1;
 bool have_coef = false;
 double coef = 1;

 while( i < t.size() ) {
  const auto & tk = t[ i ];
  if( is_op( tk ) )
   return;
  if( stop_at_label && ( i + 1 < t.size() ) && ( t[ i + 1 ] == ":" ) )
   return;
  if( tk == "+" ) { sign = 1; ++i; continue; }
  if( tk == "-" ) { sign = -1; ++i; continue; }

  if( tk == "[" ) {  // quadratic part, up to "]" and an optional "/ 2"
   ++i;
   std::map< std::string , double > q;
   double qs = 1;
   bool qhave = false;
   double qc = 1;
   while( ( i < t.size() ) && ( t[ i ] != "]" ) ) {
    const auto & u = t[ i ];
    double v;
    if( u == "+" ) { qs = 1; ++i; continue; }
    if( u == "-" ) { qs = -1; ++i; continue; }
    if( to_number( u , v ) ) { qc = v; qhave = true; ++i; continue; }
    // a name: "x ^ 2" or "x * y"
    std::string a = u , b = u;
    ++i;
    if( ( i < t.size() ) && ( t[ i ] == "^" ) )
     i += 2;
    else
     if( ( i < t.size() ) && ( t[ i ] == "*" ) ) {
      b = t[ i + 1 ];
      i += 2;
      }
    if( b < a )
     std::swap( a , b );
    q[ a + "*" + b ] += qs * ( qhave ? qc : 1 );
    qs = 1;
    qhave = false;
    }
   ++i;  // "]"
   double div = 1;
   if( ( i + 1 < t.size() ) && ( t[ i ] == "/" ) ) {
    to_number( t[ i + 1 ] , div );
    i += 2;
    }
   for( const auto & e : q )
    qd[ e.first ] += sign * e.second / div;
   sign = 1;
   continue;
   }

  double v;
  if( to_number( tk , v ) ) {
   ++i;
   // a number followed by a name is a coefficient, otherwise a constant
   if( ( i < t.size() ) && ( ! is_op( t[ i ] ) ) && ( t[ i ] != "+" ) &&
       ( t[ i ] != "-" ) && ( t[ i ] != "[" ) &&
       ! ( stop_at_label && ( i + 1 < t.size() ) && ( t[ i + 1 ] == ":" ) ) ) {
    double w;
    if( ! to_number( t[ i ] , w ) ) {
     coef = v;
     have_coef = true;
     continue;
     }
    }
   cst += sign * v;
   sign = 1;
   continue;
   }

  // a name
  if( tk == "Constant" )   // the constant of the Objective, as a column
   cst += sign * ( have_coef ? coef : 1 );
  else
   lin[ tk ] += sign * ( have_coef ? coef : 1 );
  sign = 1;
  have_coef = false;
  ++i;
  }
 }

/// sets a bound of c from "name op value" (dir of the op as seen from name)
inline void set_bound( Column & c , int dir , double v )
{
 if( dir < 0 )
  c.hi = v;
 else
  if( dir > 0 )
   c.lo = v;
  else
   c.lo = c.hi = v;
 }

}  // end( namespace detail )

/*--------------------------------------------------------------------------*/
/// reads the LP file fname into its canonical form; throws
/// std::runtime_error if the file cannot be read or parsed

inline Model read_lp( const std::string & fname )
{
 using namespace detail;

 std::ifstream f( fname );
 if( ! f.is_open() )
  throw( std::runtime_error( "lp_compare::read_lp: cannot open " + fname ) );

 // split the file in sections, each a sequence of tokens
 std::map< Section , std::vector< std::string > > tok;
 Section cur = eNone;
 bool maximize = false;
 std::string line;
 while( std::getline( f , line ) ) {
  const auto b = line.find_first_not_of( " \t\r" );
  if( ( b == std::string::npos ) || ( line[ b ] == '\\' ) )
   continue;  // empty or comment
  const auto s = section_of( line , maximize );
  if( s != eNone ) {
   cur = s;
   if( cur == eEnd )
    break;
   continue;
   }
  if( ( cur == eNone ) || ( cur == eIgnore ) )
   continue;
  auto & v = tok[ cur ];
  tokenize( line , v );
  // a line of Bounds, Binaries or Generals ends there
  if( cur == eBounds )
   v.push_back( ";" );
  }

 Model m;
 m.minimize = ! maximize;

 auto col = [ & ]( const std::string & n ) -> Column & {
  return( m.cols[ n ] );
  };

 // Objective
 {
  auto & t = tok[ eObj ];
  std::size_t i = 0;
  if( ( t.size() > 1 ) && ( t[ 1 ] == ":" ) )
   i = 2;
  read_expr( t , i , m.obj , m.qobj , m.obj_const , false );
  for( const auto & e : m.obj )
   col( e.first );
  }

 // rows
 {
  auto & t = tok[ eRows ];
  std::size_t i = 0;
  std::size_t anon = 0;
  while( i < t.size() ) {
   Row r;
   if( ( i + 1 < t.size() ) && ( t[ i + 1 ] == ":" ) ) {
    r.name = t[ i ];
    i += 2;
    }
   else
    r.name = "_row" + std::to_string( anon++ );

   std::map< std::string , double > lin , qd;
   double cst = 0;
   double v;
   // "lo <= expr <= hi" (or with >=)
   double first = 0;
   int first_dir = 2;
   if( ( i + 1 < t.size() ) && to_number( t[ i ] , v ) &&
       is_op( t[ i + 1 ] ) ) {
    first = v;
    first_dir = op_dir( t[ i + 1 ] );
    i += 2;
    }
   read_expr( t , i , lin , qd , cst , true );
   if( ( i + 1 >= t.size() ) || ! is_op( t[ i ] ) )
    throw( std::runtime_error( "lp_compare::read_lp: row " + r.name +
                               " has no sense in " + fname ) );
   const int dir = op_dir( t[ i ] );
   double rhs;
   if( ! to_number( t[ i + 1 ] , rhs ) )
    throw( std::runtime_error( "lp_compare::read_lp: row " + r.name +
                               " has no right-hand side in " + fname ) );
   i += 2;
   rhs -= cst;
   if( dir < 0 )
    r.hi = rhs;
   else
    if( dir > 0 )
     r.lo = rhs;
    else
     r.lo = r.hi = rhs;
   if( first_dir != 2 ) {  // first <op> expr
    if( first_dir < 0 )
     r.lo = first - cst;
    else
     if( first_dir > 0 )
      r.hi = first - cst;
     else
      r.lo = r.hi = first - cst;
    }
   for( const auto & e : lin ) {
    col( e.first );
    if( e.second != 0 )
     r.terms.emplace_back( e.first , e.second );
    }
   if( ! qd.empty() )
    for( const auto & e : qd )
     r.terms.emplace_back( "[" + e.first + "]" , e.second );
   std::sort( r.terms.begin() , r.terms.end() );
   r.group = group_of( r.name );
   m.rows.push_back( std::move( r ) );
   }
  }

 // bounds, one per line
 {
  auto & t = tok[ eBounds ];
  std::size_t i = 0;
  while( i < t.size() ) {
   std::vector< std::string > l;
   while( ( i < t.size() ) && ( t[ i ] != ";" ) )
    l.push_back( t[ i++ ] );
   ++i;
   if( l.empty() )
    continue;
   double v , w;
   std::string lw = l.back();
   std::transform( lw.begin() , lw.end() , lw.begin() ,
                   []( unsigned char c ) { return( std::tolower( c ) ); } );
   if( ( l.size() == 2 ) && ( lw == "free" ) ) {
    auto & c = col( l[ 0 ] );
    c.lo = - INF;
    c.hi = INF;
    continue;
    }
   if( ( l.size() == 3 ) && is_op( l[ 1 ] ) ) {
    if( to_number( l[ 2 ] , v ) )         // name op v
     set_bound( col( l[ 0 ] ) , op_dir( l[ 1 ] ) , v );
    else
     if( to_number( l[ 0 ] , v ) )        // v op name
      set_bound( col( l[ 2 ] ) , - op_dir( l[ 1 ] ) , v );
    continue;
    }
   if( ( l.size() == 5 ) && is_op( l[ 1 ] ) && is_op( l[ 3 ] ) &&
       to_number( l[ 0 ] , v ) && to_number( l[ 4 ] , w ) ) {
    auto & c = col( l[ 2 ] );
    set_bound( c , - op_dir( l[ 1 ] ) , v );
    set_bound( c , op_dir( l[ 3 ] ) , w );
    continue;
    }
   std::string s;
   for( const auto & x : l )
    s += x + " ";
   throw( std::runtime_error( "lp_compare::read_lp: bound line \"" + s +
                              "\" not understood in " + fname ) );
   }
  }

 for( const auto & n : tok[ eInt ] )
  col( n ).integer = true;
 for( const auto & n : tok[ eBin ] ) {
  auto & c = col( n );
  c.integer = true;
  c.lo = std::max( c.lo , 0.0 );
  c.hi = std::min( c.hi , 1.0 );
  }

 // fold the slack column of a ranged row ("Rg<row name>") into the row
 std::map< std::string , std::size_t > by_name;
 for( std::size_t k = 0 ; k < m.rows.size() ; ++k )
  by_name[ m.rows[ k ].name ] = k;
 std::vector< std::string > folded;
 for( const auto & c : m.cols ) {
  if( c.first.compare( 0 , 2 , "Rg" ) || m.obj.count( c.first ) )
   continue;
  const auto it = by_name.find( c.first.substr( 2 ) );
  if( it == by_name.end() )
   continue;
  auto & r = m.rows[ it->second ];
  if( r.lo != r.hi )
   continue;
  auto pos = std::find_if( r.terms.begin() , r.terms.end() ,
                           [ & ]( const auto & e ) {
                            return( e.first == c.first ); } );
  if( pos == r.terms.end() )
   continue;
  // expr + a * s = rhs with s in [ l , u ]: expr in rhs - a * [ l , u ]
  const double a = pos->second;
  const double rhs = r.lo;
  const double x = rhs - a * c.second.lo , y = rhs - a * c.second.hi;
  r.lo = std::min( x , y );
  r.hi = std::max( x , y );
  r.terms.erase( pos );
  folded.push_back( c.first );
  }
 for( const auto & n : folded )
  m.cols.erase( n );

 return( m );
 }

/*--------------------------------------------------------------------------*/
/*----------------------------- COMPARING ----------------------------------*/
/*--------------------------------------------------------------------------*/

inline std::string row_string( const Row & r )
{
 std::ostringstream s;
 s << r.name << ": " << num( r.lo ) << " <= ";
 for( std::size_t k = 0 ; k < r.terms.size() ; ++k )
  s << ( k ? " + " : "" ) << num( r.terms[ k ].second ) << " "
    << r.terms[ k ].first;
 s << " <= " << num( r.hi );
 return( s.str() );
 }

/// what differs between a row of A and a row of B (usually the one with
/// the same name)
inline std::string row_delta( const Row & a , const Row & b , double tol )
{
 std::ostringstream s;
 s << a.name << ":";
 if( ! same( a.lo , b.lo , tol ) )
  s << " lo A " << num( a.lo ) << " B " << num( b.lo ) << ";";
 if( ! same( a.hi , b.hi , tol ) )
  s << " hi A " << num( a.hi ) << " B " << num( b.hi ) << ";";
 std::map< std::string , std::pair< double , double > > c;
 for( const auto & e : a.terms )
  c[ e.first ].first = e.second;
 for( const auto & e : b.terms )
  c[ e.first ].second = e.second;
 for( const auto & e : c )
  if( ! same( e.second.first , e.second.second , tol ) )
   s << " " << e.first << " A " << num( e.second.first ) << " B "
     << num( e.second.second ) << ";";
 return( s.str() );
 }

/*--------------------------------------------------------------------------*/
/// compares the models a and b; the rows of the groups in skip (e.g., the
/// dynamic groups of separated cuts) are left out of the comparison

inline Diff compare( const Model & a , const Model & b , double tol = 1e-7 ,
                     const std::set< std::string > & skip = {} )
{
 Diff d;

 // columns: existence, bounds, integrality
 for( const auto & c : a.cols ) {
  const auto it = b.cols.find( c.first );
  if( it == b.cols.end() ) {
   ++d.cols[ group_of( c.first ) ];
   d.col_examples.push_back( c.first + " only in A" );
   continue;
   }
  if( ( ! same( c.second.lo , it->second.lo , tol ) ) ||
      ( ! same( c.second.hi , it->second.hi , tol ) ) ||
      ( c.second.integer != it->second.integer ) ) {
   ++d.cols[ group_of( c.first ) ];
   d.col_examples.push_back( c.first + " A [" + num( c.second.lo ) + "," +
                             num( c.second.hi ) + "]" +
                             ( c.second.integer ? "i" : "" ) + " B [" +
                             num( it->second.lo ) + "," +
                             num( it->second.hi ) + "]" +
                             ( it->second.integer ? "i" : "" ) );
   }
  }
 for( const auto & c : b.cols )
  if( ! a.cols.count( c.first ) ) {
   ++d.cols[ group_of( c.first ) ];
   d.col_examples.push_back( c.first + " only in B" );
   }

 // Objective
 auto cmp_map = [ & ]( const std::map< std::string , double > & x ,
                       const std::map< std::string , double > & y ,
                       const std::string & what ) {
  std::set< std::string > keys;
  for( const auto & e : x ) keys.insert( e.first );
  for( const auto & e : y ) keys.insert( e.first );
  for( const auto & k : keys ) {
   const double u = x.count( k ) ? x.at( k ) : 0;
   const double v = y.count( k ) ? y.at( k ) : 0;
   if( ! same( u , v , tol ) ) {
    ++d.n_obj;
    d.obj_examples.push_back( what + k + " A " + num( u ) + " B " +
                              num( v ) );
    }
   }
  };
 cmp_map( a.obj , b.obj , "" );
 cmp_map( a.qobj , b.qobj , "quad " );
 if( ! same( a.obj_const , b.obj_const , tol ) ) {
  ++d.n_obj;
  d.obj_examples.push_back( "constant A " + num( a.obj_const ) + " B " +
                            num( b.obj_const ) );
  }
 if( a.minimize != b.minimize ) {
  ++d.n_obj;
  d.obj_examples.push_back( "sense differs" );
  }

 // rows: bucketed by support and finite sides, matched within tol
 auto key = []( const Row & r ) {
  std::string k = std::isinf( r.lo ) ? "-" : "L";
  k += std::isinf( r.hi ) ? "-" : "U";
  for( const auto & e : r.terms )
   k += "|" + e.first;
  return( k );
  };
 auto equal = [ & ]( const Row & x , const Row & y ) {
  if( ( ! same( x.lo , y.lo , tol ) ) || ( ! same( x.hi , y.hi , tol ) ) )
   return( false );
  for( std::size_t k = 0 ; k < x.terms.size() ; ++k )
   if( ! same( x.terms[ k ].second , y.terms[ k ].second , tol ) )
    return( false );
  return( true );
  };

 std::map< std::string , std::vector< const Row * > > bucket_b;
 std::map< std::string , const Row * > b_by_name , a_by_name;
 for( const auto & r : b.rows ) {
  if( skip.count( r.group ) ) { ++d.skipped_b; continue; }
  ++d.rows[ r.group ].in_b;
  bucket_b[ key( r ) ].push_back( & r );
  b_by_name[ r.name ] = & r;
  }
 for( const auto & r : a.rows )
  if( ! skip.count( r.group ) )
   a_by_name[ r.name ] = & r;

 std::set< const Row * > used;
 std::vector< const Row * > only_a;
 for( const auto & r : a.rows ) {
  if( skip.count( r.group ) ) { ++d.skipped_a; continue; }
  ++d.rows[ r.group ].in_a;
  bool found = false;
  const auto it = bucket_b.find( key( r ) );
  if( it != bucket_b.end() )
   for( const auto * q : it->second )
    if( ( ! used.count( q ) ) && equal( r , * q ) ) {
     used.insert( q );
     found = true;
     break;
     }
  if( ! found )
   only_a.push_back( & r );
  }

 for( const auto * r : only_a ) {
  auto & g = d.rows[ r->group ];
  ++g.only_a;
  const auto it = b_by_name.find( r->name );
  if( ( it != b_by_name.end() ) && ! used.count( it->second ) )
   g.examples.push_back( row_delta( * r , * it->second , tol ) );
  else
   g.examples.push_back( "only in A: " + row_string( * r ) );
  }
 for( const auto & r : b.rows ) {
  if( skip.count( r.group ) || used.count( & r ) )
   continue;
  auto & g = d.rows[ r.group ];
  ++g.only_b;
  const auto it = a_by_name.find( r.name );
  if( it == a_by_name.end() )
   g.examples.push_back( "only in B: " + row_string( r ) );
  }

 return( d );
 }

}  // end( namespace lp_compare )

/*--------------------------------------------------------------------------*/

#endif  /* lp_compare.h included */

/*--------------------------------------------------------------------------*/
/*----------------------- End File lp_compare.h ----------------------------*/
/*--------------------------------------------------------------------------*/
