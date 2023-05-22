# Data Generator Tool

This is a tool to generate UC data to SMS++.  

## Getting started

The main file is UC2SMSpp.m which starts with two switchable categories as below:

1. % Spinning reserves: 

%   None                =====> There are no Primary and Secondary Spinning reserves

%   Prim                =====> There are just Primary Spinning reserves

%   PplusS              =====> There are both Primary and Secondary Spinning reserves


2. % Thermal:

%ThermOpt = 'noHydro';  =====> There are Just ThermalUnits

%ThermOpt = 'wHydro';   =====> There are ThermalUnits and HydroUnits

%ThermOpt = 'pHydro';   =====> Some subset of cascading systems -- to use with Caution


In each category, to generate each kind of instance, it's enough comment other options.
Thus, it gives a combinatorial option of picking "k" among "n".

The folder apogene_grace contains 7 different base data sets which are essential to generate the UC data to SMS++.
For that matter, the use should load them by changing the path (in lines 56-69) as his/him system.
Besides, user should give a path(in line 158) to where the SMSpp dataset gets written down.


### Requirements

MATLAB 

## Contributing

This section is not ready yet.

## Authors

### Current Lead Authors

- **Wim van Ackooij**  
   Expert Researcher at EDF R&D
  

## License

This section is not ready yet. 


