/*--------------------------------------------------------------------------*/
/*---------------------------- File HeatBlock.h ----------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *derived* class HeatBlock, which derives from
 * the Block, in order to define a base for any possible unit that produces
 * heat and can be attached to a UCBlock. It has very basic information that
 * can characterize almost any different kind of heat-unit, which includes a
 * set of heat variables. This class has thus been constructed having the
 * following elements:
 *
 * - A virtual public method that is used to initialize and read the
 *   data of any possible derived HeatBlock class.
 *
 * - The time horizon of the problem.
 *
 * - a vector of ColVariable objects, that are used to store the
 *   information regarding:
 *
 *     (i)   the heat produced by the unit.
 *
 *   This vector either have size equal to the time horizon
 *   or is empty, in which case the corresponding variables simply do
 *   not exist (for instance, the unit may not produce heat).
 *
 * \version 0.11
 *
 * \date 21 - 05 - 2019
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Ali Ghezelsoflu \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Rafael Durbano Lobato \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Antonio Frangioni, Ali Ghezelsoflu, and Rafael
 * Durbano Lobato
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef __HeatBlock
#define __HeatBlock
/* self-identification: #endif at the end of the file */
/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "Block.h"
#include "ColVariable.h"
#include <map>

/*--------------------------------------------------------------------------*/
/*------------------------------ NAMESPACE ---------------------------------*/
/*--------------------------------------------------------------------------*/

/// namespace for the Structured Modeling System++ (SMS++)

namespace SMSpp_di_unipi_it {


/*--------------------------------------------------------------------------*/
/*-------------------------- CLASS HeatBlock -------------------------------*/
/*--------------------------------------------------------------------------*/
/*--------------------------- GENERAL NOTES --------------------------------*/
/*--------------------------------------------------------------------------*/

/// implementation of the Block concept for the heat unit problem
/** The HeatBlock class implements the Block concept [see Block.h]
 * for the EDF Unit Commitment Problem.
 *
 *
 *
 *
 */

    class HeatBlock : public Block {


/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:

/*--------------------------------------------------------------------------*/
/*---------------------- PUBLIC TYPES OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Public types
 *
 * HeatBlock defines the following main public types:
 *
 * - Index, the type of parameters indices;
 *
 * @{ */

/*--------------------------------------------------------------------------*/

typedef unsigned int Index;

/*@}------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */

/** Constructor of HeatBlock, taking possibly a pointer of its
 * father Block. */

HeatBlock( Block * father_block = nullptr , int t = 0 )
        : Block( father_block ), f_time_horizon( t ) {}

/*--------------------------------------------------------------------------*/

/// destructor of HeatBlock

virtual ~HeatBlock() { };


/*@}------------------------------------------------------------------------*/
/*-------------------------- OTHER INITIALIZATIONS -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Other initializations
 *  @{ */

/// loads the HeatBlock instance from memory
/** Loads the HeatBlock instance from memory.
 * Like load( std::istream & ), if there is any Solver attached to this
 * HeatBlock then a NBModification (the "nuclear option") is issued.
 * */

virtual void load( std::istream &input ) override { };

/*--------------------------------------------------------------------------*/
/// extends Block::deserialize( netCDF::NcGroup )
/** Extends Block::deserialize( netCDF::NcGroup ) to the specific format of
 * the HeatBlock. Besides the mandatory "type" attribute of any :Block,
 * the group should contain the following:
 *
 * - the dimension "TimeHorizon" containing the time horizon;
 *
 * - the dimension "NumberHeatUnits" containing the number heat-producing
 *   units;
 *
 * NO, THIS WE DON'T NEED
 * - the dimension "NumberHeatOnlyUnits" containing the number of
 *   heat-only generation units in the problem; the dimension is
 *   optional: if it is not provided then it is taken to be 0, which
 *   means that there is no heat-only generation unit;
 *
 * - the dimension "NumberHeatID" containing the number of heat-ID in the
 *   problem; the dimension is optional: if it is not provided then it is
 *   taken to be 0, which means there is no heat-ID in the problem;
 *
 * - the dimension "NumberThermalStorage" containing the number of thermal
 *   storage units within an energy cell; the dimension is optional: if it is
 *   not provided then it is taken to be 0, which means no thermal storage
 *   units are presented in the problem TODO;
 *
 * - the variable "UnitEnergyCell", of type int and indexed over the dimension
 *   "NumberUnits"; the entry UnitEnergyCell[ i ] tells to which heat block unit
 *   i belongs; if UnitEnergyCell[ i ] >= NumberHeatBlocks, this means
 *   that unit i is not located in any heat block; if NumberEnergyCell == 0
 *   (say, it is not provided at all) then this variable need not be defined,
 *   since it is not loaded;
 *
 * - the variable "UnitHeatID", of type int and indexed over the dimension
 *   "NumberUnits"; the entry UnitHeatID[ i ] tells to which heat-ID unit i
 *   is assigned; if UnitHeatID[ i ] >= NumberHeatID, this means that unit i is
 *   not assigned to any heat-ID; if NumberHeatID == 0 (say, it is not
 *   provided at all), then this variable need not be defined, since it is not
 *   loaded;
 *
 * - the variable "HeatDemand", of type double and indexed both over the
 *   dimensions "NumberHeatID" and "TimeHorizon": entry HeatDemand[ i , t ] is
 *   assumed to contain the heat demand which are specified on the same heat-id
 *   i in the time t; if NumberHeatID == 0 (say, it is not provided at all),
 *   then this variable need not be defined, since it is not loaded;
 *
 * - the variable "HeatRho", of type double and indexed both over the
 *   dimensions "TimeHorizon" and "??": entry HeatRho[ t , i ] is
 *   assumed to contain the power to heat ratio at time t; TODO..;
 *
 * - the variable "ThermalStorage", of type int and indexed over the dimension
 *   "NumberUnits"; the entry ThermalStorage[ i ] tells to which thermal unit
 *   i belongs; if ThermalStorage[ i ] >= NumberThermalStorage, this means
 *   that unit i is not assigned to any thermal storage; if
 *   NumberThermalStorage == 0 (say, it is not provided at all), then this
 *   variable need not be defined, since it is not loaded TODO;
 *
 * - the scalar variable "MinHeat", of type double and not indexed over
 *   any dimension and indicates the minimum heat output value of the unit;
 *
 * - the scalar variable "MaxHeat", of type double and not indexed over
 *   any dimension and indicates the maximum heat output value of the unit;
 */

virtual void deserialize( netCDF::NcGroup & group ) override;

/*--------------------------------------------------------------------------*/

virtual void generate_abstract_variables( Configuration *stvv = nullptr )
override;

/// generate the abstract variables of the HeatUnit
/** Method that generates the abstract variables of the HeatBlock.
 * These are as std::vector< ColVariable >  with exactly :
 *
 *
 * */
/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// generate the static constraint of the HeatBlock
/** Method that generates the static constraint of the HeatBlock.
 * These are the:
 *
 */

virtual void generate_abstract_constraints( Configuration *stcc = nullptr )
override ;



/*- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -*/
/// generate the objective of the HeatBlock
/** Method that generates the objective of the HeatBlock.
 *
*/
virtual void generate_objective( Configuration *objc = nullptr )
    override ;
/*@} -----------------------------------------------------------------------*/
/*-------------- Methods for reading the data of the HeatBlock -------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for reading the data of the HeatBlock
 *  @{  */

/// returns the time horizon of the problem
Index get_time_horizon( void ) const { return f_time_horizon; }

/// returns the heat rho of the given unit at the given time
inline double get_heat_rho( Index unit, Index time ) const {
 return v_heat_rho[ unit * f_time_horizon + time ];
}

/// returns the heat demand of the given zone at the given time
inline double get_heat_demand( Index unit, Index time ) const {
  return v_heat_demand[ unit * f_time_horizon + time ];
}

/// Method for returning the vector of heat variables
const std::vector<ColVariable> & get_heat( void ) const {
     return v_heat;
}

/// Method for returning the pointer to the heat variable at time t
ColVariable * get_heat( int i ) { return & ( v_heat[i] ); }


/*@} -----------------------------------------------------------------------*/
/*--------------------- METHODS FOR SAVING THE HeatBlock -------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for loading, printing & saving the HeatBlock
 *  @{ */

/// extends Block::serialize( netCDF::NcGroup )
/** Extends Block::serialize( netCDF::NcGroup ) to the specific format of a
 * HeatBlock. See HeatBlock::deserialize( netCDF::NcGroup ) for
 * details of the format of the created netCDF group.
 *
 * */

virtual void serialize( netCDF::NcGroup & group ) const override;


/*@} -----------------------------------------------------------------------*/
/*----------------- METHODS FOR MODIFYING THE HeatBlock --------------------*/
/*--------------------------------------------------------------------------*/
/** @name Methods for modifying the HeatBlock
 *  @{ */

/// Set the time horizon
void set_time_horizon( int t );

/*@} -----------------------------------------------------------------------*/
/*-------------------- PROTECTED PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

protected:

/// The time horizon of the problem
Index f_time_horizon;

/** the matrix of HeatDemand indexed over the dimensions
* NumberHeatID and TimeHorizon */
std::vector<double> v_heat_demand;

/// the MinHeat value
int f_MinHeat;

/// the MaxHeat value
int f_MaxHeat;

/// Vector of heat rho
std::vector<double> v_heat_rho;

/// Vector of heat variables
std::vector<ColVariable> v_heat;
/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED FRIENDS -----------------------------*/
/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/
/*-------------------------- PROTECTED METHODS -----------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Protected methods for inserting and extracting
* @{ */
/// print the HeatBlock on an ostream with the given verbosity
/** Protected method to print information about the HeatBlock;
 * */

// virtual void print( std::ostream &output ) const override ;

/*--------------------------------------------------------------------------*/
/// loads the HeatBlock instance from standard .dat file format
/** Protected method for loading a HeatBlock out of a std::istream
 *
 *      //TODO
 */

//virtual void load( std::istream &input ) override final;


// void load(std::istream& inStream);

/*@}------------------------------------------------------------------------*/
/*-------------------- PROTECTED FIELDS OF THE CLASS -----------------------*/
/*--------------------------------------------------------------------------*/

/*-----------------------------variables------------------------------------*/




/*----------------------------constraints-----------------------------------*/


/*--------------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/
    private:

        SMSpp_insert_in_factory_h;
/*--------------------------------------------------------------------------*/
/*-------------------------- PRIVATE METHODS -------------------------------*/
/*--------------------------------------------------------------------------*/

/// returns which variables must be generated
/** This method returns an int that indicates which variables of
* HeatBlock must be generated by the generate_abstract_variables()
* method. This value may be given in stvv as explained in
* generate_abstract_variables(). If this value is not given in stvv,
* then this method returns the appropriate value according to what
* is specified in the generate_abstract_variables() method. */
int get_variables_to_be_generated( Configuration *stvv );

    }; // end( class( HeatBlock ) )

/*@}  end( class( HeatBlock ) ) --------------------------------------------*/
/*--------------------------------------------------------------------------*/
}; // end( namespace SMSpp_di_unipi_it )

/*--------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------*/

#endif /* HeatBlock.h included */

/*--------------------------------------------------------------------------*/
/*-------------------------- End File HeatBlock.h --------------------------*/
/*--------------------------------------------------------------------------*/
