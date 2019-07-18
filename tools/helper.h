#ifndef __DAT2NC4
#define __DAT2NC4

#include <vector>

/// A thermal unit as represented in a DAT or in a MOD file
struct ThermalUnit {
 unsigned int index{};

 double QuadTerm{};
 double LinearTerm{};
 double ConstTerm{};
 double MinPower{};
 double MaxPower{};
 double InitialPower{}; // aka PZERO
 double InitUpDownTime{};
 double MinUpTime{};
 double MinDownTime{};
 double DeltaRampUp{};
 double DeltaRampDown{};

 // Only in MOD files
 double coolAndFuelCost{};
 double hotAndFuelCost{};
 double tau{};
 double tauMax{};
 double fixedCost{};
 double SUCC{};

 // Only in DAT files
 double StartUpCost{};
 double BoundOn{};
 double BoundDown{};

 /// Loads the data from a line of a MOD file
 void load( std::istream & in );

 /// Prints the data
 void print( std::ostream & out ) const;

 /// Generates StartUpCost from
 void generate_startupcost();

 friend std::ostream & operator<<( std::ostream & out, const ThermalUnit & unit );

 friend std::istream & operator>>( std::istream & in, ThermalUnit & unit );
};

/// A hydro unit as represented in a MOD file
struct HydroUnit {
};

/// A hydro cascade unit as represented in a MOD file
struct HydroCascadeUnit {
};

/// A load curve as represented in a MOD file
struct LoadCurve {
 double MinSystemCapacity{};
 double MaxSystemCapacity{};
 double MaxThermalCapacity{};
 std::vector< std::vector< double>> Loads;
 std::vector< double > SpinningReserve;
};

/// A DAT file containing a single thermal unit
struct DatFile {
 unsigned int TimeHorizon{};
 ThermalUnit thermal_unit;
 std::vector< double > Lambda;
 std::vector< double > Mu;

 /// Generates the linear and constant coefficients of the cost function
 void generate_bc( std::vector< double > & b, std::vector< double > & c );

 /// Loads the data from a DAT file
 void load( std::istream & in );

 /// Prints the data
 void print( std::ostream & out ) const;

 friend std::ostream & operator<<( std::ostream & out, const DatFile & file );

 friend std::istream & operator>>( std::istream & in, DatFile & file );
};

/// A MOD file containing a load curve and multiple units
struct ModFile {
 unsigned int ProblemNum{};
 unsigned int TimeHorizon{};
 unsigned int NumThermal{};
 unsigned int NumHydro{};
 unsigned int NumCascade{};

 LoadCurve load_curve;
 std::vector< ThermalUnit > thermal_units;
 std::vector< HydroUnit > hydro_units;
 std::vector< HydroCascadeUnit > hydro_cascade_units;

 /// Loads the data from a MOD file
 void load( std::istream & in );

 /// Prints the data
 void print( std::ostream & out ) const;

 friend std::ostream & operator<<( std::ostream & out, const ModFile & file );

 friend std::istream & operator>>( std::istream & in, ModFile & file );
};

enum filetype {
 ftDat = 0,
 ftMod = 1
};

#endif //__DAT2NC4
