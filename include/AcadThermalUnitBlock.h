/*--------------------------------------------------------------------------*/
/*------------------------ File AcadThermalUnitBlock.h ---------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *derived* class AcadThermalUnitBlock, which derives from 
 * UnitBlock, in order to define the Thermal Units of the Academic Version of 
 * Unit Commitement Problem. This class is designed in order to produce four dif-
 * ferent type of formulations of 1UC and for each one of them a p/c cuts formu-
 * lation is also considered.
 *
 * All different type of formulations are characterized by
 *
 * - A min & max power production of the Units (Pmin/Pmax Constraints)
 * - A min & max up and down duration of the Units (Min Up/Down Constraints)
 * - A min & max increase & decrease of power of the Units (Ramp Constraints)
 * - An initial state of the commitement Variable (Initial Status Constraint)
 * - A diag quadratic or linear cost objective function
 * - A time-depentent start-up cost
 *
 * Based on the above description the class has been constructed having the
 * following elements:
 *
 * - A public method that reads and initializes all the data that descri-
 *   bes the AcadThermalUnitBlock instance for each different formuation 
 *   (the selection being done via the variable named model that denotes 
 *   the specific model to be used and via the variable pc that denotes if
 *   perspective cuts will be used, both being tuned via main file) with
 *   the corresponding entities of the class being poppulated. 
 *
 * - A set of public methods, one for the initialization of each different
 *   type of constraints mentioned above for each different formulations with
 *   the methods also passing the constraints to the vector of static const-
 *   raints of the Block.
 * - Two public methods for the calculation of the start-up costs of the 
 *   thermal units.
 * - A public method for the initialization of the objective function, which
 *   can be linear or diagonial quadratic the selection of which is made via
 *   the variable of which is tuned in the main file.
 * - A number of different vectors of LinearConstraint Objects that are used
 *   to store all the information of the different sets of Constraints of
 *   each thermal unit.
 * - An object of DQuadObjectiveFunction and one of LinearObjectiveFunction 
 *   that are used to store the information of the quadratic or linear objec-
 *   tive function of each unit.
 * - Several different double variables that are used in order to store all
 *   the different values needed to descibe the above mentioned constraints
 *   and costs.
 * - A small private class with sole purpose of defining a static member, which 
 *   has to be initialized when main() is executed. Τhe class constructos is 
 *   set to register AcadThermalUnitBlock in the Factory of UnitBlock.

 * \version 0.10
 *
 * \date 03 - 09 - 2016
 *
 * \author Antonio Frangioni \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * \author Kostas Tavlaridis-Gyparakis \n
 *         Operations Research Group \n
 *         Dipartimento di Informatica \n
 *         Universita' di Pisa \n
 *
 * Copyright &copy by Antonio Frangioni, Kostas Tavlaridis-Gyparakis
 */
/*--------------------------------------------------------------------------*/
/*----------------------------- DEFINITIONS --------------------------------*/
/*--------------------------------------------------------------------------*/

#ifndef AcadThermalUnitBlock_H_
#define AcadThermalUnitBlock_H_

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "LinearConstraint.h"
#include "DQuadObjectiveFunction.h"
#include "UnitBlock.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it {

class UCBlock; ///< Fowrard Declaration of UCBlock Class

class AcadThermalUnitBlock : public UnitBlock {

/*--------------------------------------------------------------------------*/
/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

public:
/*--------------------------------------------------------------------------*/
/*---------------------- PUBLIC TYPES OF THE CLASS -------------------------*/
/*--------------------------------------------------------------------------*/

enum eShutStatus {eUnknow, eCooling, eBancking}; 
///< enum for different shutting status of the unit

/*--------------------------------------------------------------------------*/
/*--------------------- CONSTRUCTOR AND DESTRUCTOR -------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constructor and Destructor
 *  @{ */
	AcadThermalUnitBlock(UCBlock * flbock = nullptr): UnitBlock( flbock ) {
	
	rampconst = false;
}
/**< Constructor of AcadThermalUnitBlock, taking possibly a pointer of its father 
  *  Block, alongside with setting the Unit without rampconstraints */

	virtual ~AcadThermalUnitBlock(){};
///< destructor of BusNetWorkBlock: it is virtual, and empty

/*@} -----------------------------------------------------------------------*/
/*-------- METHODS FOR READING THE DATA OF THE AcadThermalUnitBlock --------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the AcadThermalUnitBlock
    @{ */

         void load(std::istream& inStream);
/**< Method for receiving data and constructing the constraints and costs 
     of the unit */

/*@} -----------------------------------------------------------------------*/
/*------- METHODS FOR HANDLING THE DATA OF THE AcadThermalUnitBlock --------*/
/*--------------------------------------------------------------------------*/
/** @name Handling the data of the AcadThermalUnitBlock
    @{ */
	double get_fMaxRampDown() {return fMaxRampDown;}
	double get_fMaxRampUp() {return fMaxRampUp;}
	double get_fInitPower() {return fInitPower;}
	double get_fMinPower() {return fMinPower;}
	double get_fMaxPower() {return fMaxPower;}
	double get_fInitUpDownTime() {return fInitUpDownTime;}
	double get_fMinDownTime() {return fMinDownTime;}
	double get_fMinUpTime() {return fMinUpTime;}
	double get_fBoundOn() {return fBoundOn;}
	double get_fBoundDown() {return fBoundDown;}

	double get_fQuadTherm() {return fQuadTherm;}
	double get_fLinearTherm() {return fLinearTherm;}
	double get_fConstTherm() {return fConstTherm;}
	double get_fStartupCosts() {return fStartupCosts;}

	std::vector<double> get_lambda() {return lambda;}
	std::vector<double> get_mew() {return mew;}

/**< Above methods are giving access to variables of the class that might be
     potentially needed elsewhere outside the Academic Thermal Unit Block */




/*@} -----------------------------------------------------------------------*/
/*------ METHODS CONSTRUCTING THE CONSTRAINTS OF AcadThermalUnitBlock ------*/
/*--------------------------------------------------------------------------*/
/** @name Methods constructing the constraints of AcadThermalUnitBlock
 *  @{ */

        void generate_static_constraints( void );
        ///< Method for generating and adding the model constraints

        void generate_dynamic_constraints( void );
        ///< Method used in order to generate perspective cuts

        void zt_Constraints( void );
        /**< Method used in order to generate the initial constraints refering to 
             perspective cuts */
	
	/// Methods for Constraints of model = 0
	void p_min_max_const(double pmin, double pmax);
	void min_up_down_const (double mindown, double minup);
        void init_up_down_const (double t0, double mindown, double minup);
    	void ramp_const(double D_down, double D_up, double pmin, double pmax);
    	
	/// Method applicable for Constraints of model = {1,2,3}
        void Network_Matrix();


	/// Methods for Constraints of model = 1
        void Convex_Constraints( void );


	/// Methods for Constraints of model = 2
        void Sum_Ramp_Constraints( void );
        void Sum_Pmin_Constraints( void );
        void Sum_Bound_Constraints( void );

	/// Methods for Constraints of model = 3
        void Semi_Ramp_Constraints( void );
        void Semi_Pmin_Constraints( void );
        void Semi_Bound_Constraints( void );

    	void obj_function( );
        /**< Method for constructing the Objective Function, common for all dif-
             ferent formulations
          */

        int mat_index (int i, int j); 
        ///< auxiliatry method for mapping between indexes and matrix positions

/*@} -----------------------------------------------------------------------*/
/*------ METHODS COMPUTING THE START UP COSTS OF AcadThermalUnitBlock ------*/
/*--------------------------------------------------------------------------*/
/** @name Methods computing the start-up costs of AcadThermalUnitBlock
 *  @{ */

        double ComputeStableSCT();
        ///<Method for Computing the Stable Start Up Cost

        void StartUpCost(int downTime, double & cost, eShutStatus& shutStatus);
        ///<Method for calculating the start up cost

/*@} -----------------------------------------------------------------------*/
/*-------------------- PUBLIC FIELDS OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*---------------------- CONSTRAINTS & OBJ FUNCTION ------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constraints and Objective Function of AcadThermalUnitBlock
 *  @{ */


    std::vector<LinearConstraint > Init_Const; 
 
    /// LinearConstraints of model = {0,1,2}

    std::vector<LinearConstraint > Pmin_Const; 
    std::vector<LinearConstraint > Pmax_Const; 
    std::vector<LinearConstraint > U_Const; 
    std::vector<LinearConstraint > P_Const; 
    std::vector<LinearConstraint > RampDown_Const; 
    std::vector<LinearConstraint > RampUp_Const; 

    /// LinearConstraints of model = {1,2,3}
    std::vector<LinearConstraint > Net_Matrix; 

    /// LinearConstraints of model = 0

    std::vector<LinearConstraint > UZ_Const; 
    std::vector<LinearConstraint > J_On_Const; 
    std::vector<LinearConstraint > J_Off_Const; 

    /// LinearConstraints of model = 1
    std::vector<LinearConstraint > Bound_on_Const;
    std::vector<LinearConstraint > Bound_Down_Const;  


    /// LinearConstraints of model = 3
    std::vector<LinearConstraint > Semi_Pmin_Const;
    std::vector<LinearConstraint > Semi_Pmax_Const;
    std::vector<LinearConstraint > Semi_Bound_on_Const; 
    std::vector<LinearConstraint > Semi_Ramp_up_Const;
    std::vector<LinearConstraint > Semi_Ramp_down_Const;

    /// LinearConstraints of model = {1,2,3}
    //std::vector<LinearConstraint > Net_Matrix; //Tiziano 


    /// LinearConstraints of Perspective Cuts
    std::vector<LinearConstraint > Zt_Const; 
    std::list<LinearConstraint > PC_Const; 

	LinearObjectiveFunction obj_fun; 
        ///<Linear Objective function of the Academic Thermal Block

        DQuadObjectiveFunction qbj_fun; 
        ///<Diag Quad Objective function of the Academic Thermal Block	



/*@} -----------------------------------------------------------------------*/
/*----------------------- UNIT VARIABLES OF THE CLASS ----------------------*/
/*--------------------------------------------------------------------------*/
/** @name Variables of AcadThermalUnitBlock
 *  @{ */


     std::vector<ColVariable> U_10;  
     std::vector<ColVariable> U_11;  
     std::vector<ColVariable> U_00;  
     std::vector<ColVariable> U_01;
     /**< Vectors that store the sets of auxiliary variables of the Academic 
          Thermal Unit (refering to model = 0)
     */

     std::vector<ColVariable> Y_plus;  //vector storing y^{h,k}_+ variables (model = {1,2,3})
     std::vector<ColVariable> Y_minus; //vector storing y^{h,k}_- variables (model = {1,2,3})
   
     std::vector<ColVariable> P_hk;   //vector storing p_{h,k} variables (model = 1)

     std::vector<ColVariable> P_h;    //vector storing p_{h} variables (model = 3)

     std::vector<ColVariable> z_t;  //vector storing z_t variables for perspective cuts



/*@} -----------------------------------------------------------------------*/
/*------------------ DATA & AUXILIARY VARIABLES OF THE CLASS ---------------*/
/*--------------------------------------------------------------------------*/
/** @name Variables of AcadThermalUnitBlock
 *  @{ */


     std::vector<pair<int,int>> Y_plus_pair;
     std::vector<pair<int,int>> Y_minus_pair;
      /**< The vectors Y_plus_pair and Y_minus_pair are used in order to define
       *   for each Colvariable in vectors Y_plus/minus what are their corersponding
       *   (h -> k ) indexes. Thus for each variable in Y_plus/minus there corres-
       *   ponds in the Y_plus/minus_pair vector in the exact same position the
       *   the pair of integers denoting the corresponding (h,k).
       *   As it is clear the vectors Y_plus and Y_plus_pair have the same size
       *   the same applies for the minus vectors as well.
       */
    
     std::vector<int> v_plus;  
     std::vector<int> v_minus;  

     std::vector<int> v_p_size;  
     std::vector<int> v_p_ind;  

     std::vector<int> v_ph_num;  // time step refering to prodution p_t 
     std::vector<int> v_ph_size; // total number of #h Variables p^t_h
     std::vector<int> v_ph_ind; //index of corresponding Y_plus 
     std::vector<int> v_ph_tot; //total number of variables until that time_Step 

      struct t {
        t()
         :index()
         ,sindex()
         ,b_on()
         ,b_dn()
         ,p_mx()
      {}

        std::vector<int> index;
        std::vector<int> sindex;
        std::vector<int> b_on;
        std::vector<int> b_dn;
        std::vector<int> p_mx;
  }; //struct that holds information needed to construct the variables of model = 2
    

    std::vector< t > v_t; //vector holding the struct

    std::vector<int> r_init; //vector denoting info for initial ramp up constraints

    std::vector<double> p_k; //varibale with the different \bar_{p}_k values

    double fQuadTherm, fLinearTherm, fConstTherm, fStartupCosts; 
    ///< Variables for the obj. function's cost

    double fMaxRampDown,fMaxRampUp,fInitPower,fMinPower, fMaxPower, 
           fInitUpDownTime, fMinUpTime, fMinDownTime, fBoundOn , fBoundDown; 
    ///< Variables for the constraints of Unit

    double coolAndFuelCost, hotAndFuelCost, fixedCost, tau; 
    ///< Variables for the start up costs

    double fGetStableSCT;
    /**< the maximum time duration after which, the start up cost function 
         is considered stable to its own maximum value. */

    bool rampconst; 
    ///< bool variable to denote if current problem has ramp constraints

    int init_t;
    ///< variable denoting the time-steps unit is subjected to initial conditions
    
    int data; 
    ///< variable used to denote the version of the input data
    
    double f_start_cost;
    /**< variable used to denote the start-up cost in the case of data input version
         with langrangian multipliers */

    int pc; //variable denoting the use of p/c

    int of; //variable for denote type of objective function

    int model; //variable for retrieving the selection of the formulation


protected: 
    std::vector<double> lambda;
    std::vector<double> mew; 
    ///< vectors that store teh dual information of the problem for the unit



/*@} -----------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

 private:

/*--------------------------------------------------------------------------*/

/* Very small, "fake" class _initializer. Its only meaning is to define a
   static member _init that is initialized at the very beginning of the
   main(). Hence the constructor is called, and the constructor registers
   the AcademicThermalUnitBlock class into the (static) UnitBlock::U_factory. */

 static class _init {
 public:
  _init(); 
  }_initializer; 
 };


class AcadMod : public BlockModificationAD
{

/*----------------------- PUBLIC PART OF THE CLASS -------------------------*/

 public:

/*---------------------------- PUBLIC TYPES --------------------------------*/

 AcadMod( int mod ) : BlockModificationAD( mod ) {};

 ~AcadMod() {};

std::list<LinearConstraint> newlist;
};

} /* namespace SMSpp_di_unipi_it */


//BOOST_CLASS_EXPORT_KEY(SMSpp_di_unipi_it::AcadThermalUnitBlock)

#endif /* AcadThermalUnitBlock_H_ */
