/*--------------------------------------------------------------------------*/
/*----------------------- File AcadThermalUnitMPBlock.h --------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Header file for the *derived* class UnitBlock, which describes the MIP For-
 * mulation that derives from the DP algorithms that are destriced in class 
 * AcadThemalUnitGraphSolver. More details soon to come... Stay tuned...
 * dients:
 *
 * \version 0.10
 *
 * \date 19 - 04 - 2016
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

#ifndef AcadThermalUnitMPBlock_H_
#define AcadThermalUnitMPBlock_H_

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include "LinearConstraint.h"
#include "DQuadObjectiveFunction.h"
//#include "UnitBlock.h"
#include "Block.h"

/*--------------------------------------------------------------------------*/
/*--------------------------- NAMESPACE ------------------------------------*/
/*--------------------------------------------------------------------------*/

namespace SMSpp_di_unipi_it {

//class UCBlock; ///< Fowrard Declaration of UCBlock Class

class UCBlock; ///< forward declaration of UCBlock

class AcadThermalUnitMPBlock : public Block {

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
	AcadThermalUnitMPBlock(UCBlock * fblock = nullptr);
/**< Constructor of AcadThermalUnitBlock, taking possibly a pointer of its father 
  *  Block, alongside with setting the Unit without rampconstraints */

	virtual ~AcadThermalUnitMPBlock();
///< destructor of BusNetWorkBlock: it is virtual, and empty

/*@} -----------------------------------------------------------------------*/
/*------ METHODS FOR READING THE DATA OF THE AcadThermalUnitMPBlock -------*/
/*--------------------------------------------------------------------------*/
/** @name Reading the data of the AcadThermalUnitMPBlock
    @{ */

         void load(std::istream& inStream);
/**< Method for receiving data and constructing the constraints and costs 
     of the unit */

/*@} -----------------------------------------------------------------------*/
/*------- METHODS FOR HANDLING THE DATA OF THE AcadThermalUnitMPBlock -----*/
/*--------------------------------------------------------------------------*/
/** @name Handling the data of the AcadThermalUnitMPBlock
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
     potentially needed elsewhere outside the Academic Thermal Unit Block
*/

/*@} -----------------------------------------------------------------------*/
/*---- METHODS CONSTRUCTING THE CONSTRAINTS OF AcadThermalUnitMPBlock -----*/
/*--------------------------------------------------------------------------*/
/** @name Methods constructing the constraints of AcadThermalUnitMPBlock
 *  @{ */

  //      void generate_static_constraints( void );
        ///< Method for generating and adding the demand and reserve constraint


/*@} -----------------------------------------------------------------------*/
/*---- METHODS COMPUTING THE START UP COSTS OF AcadThermalUnitMPBlock -----*/
/*--------------------------------------------------------------------------*/
/** @name Methods computing the start-up costs of AcadThermalUnitMPBlock
 *  @{ */

        double ComputeStableSCT();
        ///<Method for Computing the Stable Start Up Cost

        void StartUpCost(int downTime, double & cost, eShutStatus& shutStatus);
        ///<Method for calculating the start up cost

        void Network_Matrix();

        int mat_index (int i, int j);

        void obj_function( );
        ///<Method for constructing the Objective Function

        void generate_static_constraints( void );
        ///< Method for generating and adding the demand and reserve constraint
      

/*@} -----------------------------------------------------------------------*/
/*-------------------- PUBLIC FIELDS OF THE CLASS --------------------------*/
/*--------------------------------------------------------------------------*/
//protected: 
/*--------------------------------------------------------------------------*/
/*---------------------- CONSTRAINTS & OBJ FUNCTION ------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Constraints and Objective Function of AcadThermalUnitMPBlock
 *  @{ */


    std::vector<LinearConstraint > Net_Matrix; 
    ///<vector to store Equations of the network Matrix

    DQuadObjectiveFunction obj_fun; 
    ///<Objective function of the Academic Thermal Block

/*@} -----------------------------------------------------------------------*/
/*----------------------- VARIABLES OF THE CLASS ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @name Variables of AcadThermalUnitBlock
 *  @{ */

     std::vector<ColVariable> Y_plus;  //vector storing y^{h,k}_+ variables
     std::vector<ColVariable> Y_minus; //vector storing y^{h,k}_- variables

     /**< The vectors Y_plus_pair and Y_minus_pair are used in order to define
       *  for each Colvariable in vectors Y_plus/minus what are their corersponding
       *  (h -> k ) indexes. Thus for each variable in Y_plus/minus there corres-
       *  ponds in the Y_plus/minus_pair vector in the exact same position the
       *  the pair of integers denoting the corresponding (h,k).
       *  As it is clear the vectors Y_plus and Y_plus_pair have the same size
       *  the same applies for the minus vectors as well.
      */
     std::vector<pair<int,int>> Y_plus_pair;
     std::vector<pair<int,int>> Y_minus_pair;

    
     std::vector<int> v_plus;  
     std::vector<int> v_minus;  

     std::vector<int> v_p_size;  
     std::vector<int> v_p_ind;  

     std::vector<ColVariable> P;    //prepei na ginei boost::multi_array<ColVariable,2>

	std::vector<LinearConstraint > Pmax_Const;
	std::vector<LinearConstraint > Bound_on_Const;
	std::vector<LinearConstraint > Bound_Down_Const;  
	std::vector<LinearConstraint > Pmin_Const;


	std::vector<LinearConstraint > RampDown_Const;  
	std::vector<LinearConstraint > RampUp_Const;    

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
         with langrangian multipliers      
    */

    std::vector<double> lambda;
    std::vector<double> mew; 
    ///< vectors that store teh dual information of the problem for the unit


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
  };
    
  std::vector< t > v_t;

/*@} -----------------------------------------------------------------------*/
/*----------------------- PRIVATE PART OF THE CLASS ------------------------*/
/*--------------------------------------------------------------------------*/

};
 




} /* namespace SMSpp_di_unipi_it */


//BOOST_CLASS_EXPORT_KEY(SMSpp_di_unipi_it::AcadThermalUnitMIPBlock)

#endif /* AcadThermalUnitBlock_H_ */
