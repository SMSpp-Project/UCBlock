
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>

#include "helper.h"

void DatFile::generate_bc( std::vector< double > & b, std::vector< double > & c ) {
 b.resize( TimeHorizon );
 c.resize( TimeHorizon );

 for( unsigned int t = 0; t < TimeHorizon; ++t ) {
  b[ t ] = thermal_unit.LinearTerm - Lambda[ t ];
  c[ t ] = thermal_unit.ConstTerm - Mu[ t ];
 }

 // If all elements are identical, we use only one value
 if( std::adjacent_find( b.begin(),
                         b.end(),
                         std::not_equal_to<>() ) == b.end() ) {
  b.resize(1);
 }

 if( std::adjacent_find( c.begin(),
                         c.end(),
                         std::not_equal_to<>() ) == c.end() ) {
  c.resize(1);
 }
}

void DatFile::load( std::istream & in ) {
 std::string skip;
 in >> skip >> TimeHorizon
    >> skip >> thermal_unit.QuadTerm
    >> skip >> thermal_unit.LinearTerm
    >> skip >> thermal_unit.ConstTerm
    >> skip >> thermal_unit.MinPower
    >> skip >> thermal_unit.MaxPower
    >> skip >> thermal_unit.InitUpDownTime
    >> skip >> thermal_unit.MinUpTime
    >> skip >> thermal_unit.MinDownTime
    >> skip >> thermal_unit.StartUpCost
    >> skip >> thermal_unit.InitialPower
    >> skip >> thermal_unit.DeltaRampUp
    >> skip >> thermal_unit.DeltaRampDown
    >> skip >> thermal_unit.BoundOn
    >> skip >> thermal_unit.BoundDown;

 Lambda.resize( TimeHorizon );
 Mu.resize( TimeHorizon );
 in >> skip;
 for( unsigned int t = 0; t < TimeHorizon; ++t ) {
  in >> Lambda[ t ];
 }
 in >> skip;
 for( unsigned int t = 0; t < TimeHorizon; ++t ) {
  in >> Mu[ t ];
 }
}

void DatFile::print( std::ostream & out ) const {

 /*
  * "Term" is misspelled in the labels because it is misspelled in the
  * original input files and we want to compare the output with them.
  */
 out << "TimeHorizon\t" << TimeHorizon << "\n"
     << "QuadTherm\t" << thermal_unit.QuadTerm << "\n"
     << "LinearTherm\t" << thermal_unit.LinearTerm << "\n"
     << "ConstTherm\t" << thermal_unit.ConstTerm << "\n"
     << "MinPower\t" << thermal_unit.MinPower << "\n"
     << "MaxPower\t" << thermal_unit.MaxPower << "\n"
     << "InitUpDownTime\t" << thermal_unit.InitUpDownTime << "\n"
     << "MinUpTime\t" << thermal_unit.MinUpTime << "\n"
     << "MinDownTime\t" << thermal_unit.MinDownTime << "\n"
     << "StartUpCost\t" << thermal_unit.StartUpCost << "\n"
     << "PZERO\t\t" << thermal_unit.InitialPower << "\n"
     << "DeltaRampUp\t" << thermal_unit.DeltaRampUp << "\n"
     << "DeltaRampDown\t" << thermal_unit.DeltaRampDown << "\n"
     << "BoundOn\t\t" << thermal_unit.BoundOn << "\n"
     << "BoundDown\t" << thermal_unit.BoundDown << "\n";

 out << "Lambda" << "\n";
 for( unsigned int t = 0; t < TimeHorizon; ++t ) {
  out << Lambda[ t ] << " ";
 }
 out << "\n";
 out << "Mu" << "\n";
 for( unsigned int t = 0; t < TimeHorizon; ++t ) {
  out << Mu[ t ] << " ";
 }
 out << "\n";
}

std::ostream & operator<<( std::ostream & out, const DatFile & file ) {
 file.print( out );
 return out;
}

std::istream & operator>>( std::istream & in, DatFile & file ) {
 file.load( in );
 return in;
}

void ThermalUnit::load( std::istream & in ) {
 std::string skip;
 in >> index
    >> QuadTerm
    >> LinearTerm
    >> ConstTerm
    >> MinPower
    >> MaxPower
    >> InitUpDownTime
    >> MinUpTime
    >> MinDownTime
    >> coolAndFuelCost
    >> hotAndFuelCost
    >> tau
    >> tauMax
    >> fixedCost
    >> SUCC
    >> InitialPower;
 // "RampConstraints"
 in >> skip >> DeltaRampUp >> DeltaRampDown;
}

void ThermalUnit::print( std::ostream & out ) const {
 out << index << "\t"
     << QuadTerm << "\t"
     << LinearTerm << "\t"
     << ConstTerm << "\t"
     << MinPower << "\t"
     << MaxPower << "\t"
     << int( InitUpDownTime ) << "\t"
     << int( MinUpTime ) << "\t"
     << int( MinDownTime ) << "\t"
     << coolAndFuelCost << "\t"
     << hotAndFuelCost << "\t"
     << int( tau ) << "\t"
     << int( tauMax ) << "\t"
     << fixedCost << "\t"
     << SUCC << "\t"
     << InitialPower << "\n";
 out << "RampConstraints\t "
     << DeltaRampUp << " \t " << DeltaRampDown << "\n";
}

void ThermalUnit::generate_startupcost() {
 if( coolAndFuelCost != 0 || hotAndFuelCost != 0 ) {
  throw ( std::invalid_argument( "Time-dependent start up costs are not allowed" ) );
 }

 StartUpCost = fixedCost;
 // double coolCost = coolAndFuelCost * ( 1 - exp( -MinDownTime / tau ) ) + fixedCost;
 // double hotCost = hotAndFuelCost * MinDownTime + fixedCost;
 // double cost = coolCost < hotCost ? coolCost : hotCost;
}

std::ostream & operator<<( std::ostream & out, const ThermalUnit & unit ) {
 unit.print( out );
 return out;
}

std::istream & operator>>( std::istream & in, ThermalUnit & unit ) {
 unit.load( in );
 return in;
}

void ModFile::load( std::istream & in ) {
 std::string skip;
 int rows, columns, elements;

 in >> skip >> ProblemNum;
 in >> skip >> TimeHorizon;
 in >> skip >> NumThermal;
 in >> skip >> NumHydro;
 in >> skip >> NumCascade;

 // "LoadCurve"
 in >> skip;
 in >> skip >> load_curve.MinSystemCapacity;
 in >> skip >> load_curve.MaxSystemCapacity;
 in >> skip >> load_curve.MaxThermalCapacity;

 // "Loads"
 in >> skip >> rows >> columns;
 load_curve.Loads.resize( rows );
 for( int r = 0; r < rows; ++r ) {
  load_curve.Loads[ r ].resize( columns );
  for( int c = 0; c < columns; ++c ) {
   in >> load_curve.Loads[ r ][ c ];
  }
 }

 // "SpinningReserve"
 in >> skip >> elements;
 load_curve.SpinningReserve.resize( elements );
 for( int e = 0; e < elements; ++e ) {
  in >> load_curve.SpinningReserve[ e ];
 }

 // "ThermalSection"
 in >> skip;
 thermal_units.resize( NumThermal );
 for( unsigned int i = 0; i < NumThermal; ++i ) {
  thermal_units[ i ].load( in );
 }

 // "HydroSection"
 in >> skip;
 hydro_units.resize( NumHydro );
 for( unsigned int i = 0; i < NumHydro; ++i ) {}

 // "HydroCascadeSection"
 in >> skip;
 hydro_cascade_units.resize( NumCascade );
 for( unsigned int i = 0; i < NumCascade; ++i ) {}
}

void ModFile::print( std::ostream & out ) const {
 out << std::fixed << std::setprecision( 6 );

 out << "ProblemNum " << ProblemNum << "\n"
     << "HorizonLen " << TimeHorizon << "\n"
     << "NumThermal " << NumThermal << "\n"
     << "NumHydro " << NumHydro << "\n"
     << "NumCascade " << NumCascade << "\n";

 out << "LoadCurve\n"
     << "MinSystemCapacity \t " << load_curve.MinSystemCapacity << "\n"
     << "MaxSystemCapacity \t " << load_curve.MaxSystemCapacity << "\n"
     << "MaxThermalCapacity\t " << load_curve.MaxThermalCapacity << "\n";

 out << "Loads\t"
     << load_curve.Loads.size() << "\t"
     << load_curve.Loads[ 0 ].size() << "\n";
 for( auto & row : load_curve.Loads ) {
  for( double c : row ) {
   out << c << "\t";
  }
  out << "\n";
 }

 out << "SpinningReserve\t"
     << load_curve.SpinningReserve.size() << " \n";
 for( double c : load_curve.SpinningReserve ) {
  out << c << "\t";
 }
 out << "\n";

 out << "ThermalSection\n";
 for( auto & unit : thermal_units ) {
  unit.print( out );
 }

 out << "HydroSection\n";
 out << "HydroCascadeSection\n";

}

std::ostream & operator<<( std::ostream & out, const ModFile & file ) {
 file.print( out );
 return out;
}

std::istream & operator>>( std::istream & in, ModFile & file ) {
 file.load( in );
 return in;
}
