/*--------------------------------------------------------------------------*/
/*----------------------------- File test.cpp ------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Small main() for testing SMS++.
 *
 * \version 0.10
 *
 * \date 08 - 05 - 2016
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
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <experimental/filesystem>
#include <boost/lambda/bind.hpp>
#include <chrono>
namespace fs = std::experimental::filesystem;
using namespace boost::lambda;

#include "SMSTypedefs.h"
#include "ColVariable.h"
#include "LinearObjectiveFunction.h"
#include "LinearConstraint.h"
#include "UCBlock.h"
#include "UnitBlock.h"
#include "NetworkBlock.h"
#include "Network.h"
#include "NetworkNode.h"
//#include "BusNetworkBlock.h"
#include "ThermalUnitBlock.h"
//#include "AcadThemalUnitWRSolver.h"
//#include "AcadThemalUnitGraphSolver.h"
//#include "MILPSolver.h"
//#include "AcadThermalUnitMIPBlock.h"
//#include "AcadThermalUnitMPBlock.h"

/*--------------------------------------------------------------------------*/
/*-------------------------------- USING -----------------------------------*/
/*--------------------------------------------------------------------------*/

using namespace std;
using namespace boost;
using namespace SMSpp_di_unipi_it;



int dom; //external variable to choose the type of formulation that will be used for the thermal units
int choice; //external variable to choose if the whole UC or the 1UC will be solved
int t_of; //external variable to choose if the objective function will be linear (set to 0) or diagonial quadratic (set to 1)
int t_pc; //external variable to choose if perspective cuts will be used or not


/*--------------------------------------------------------------------------*/
/*--------------------------------- Main -----------------------------------*/
/*--------------------------------------------------------------------------*/

#define _GLIBCXX_USE_CXX11_ABI 0

int main() {

	choice = 0 ;

	
	if (choice == 0){

	  /*****Working and Testing UC for MILPSolver*******/
	  std::string path = "/home/tiziano/Scrivania/UC_K/UC/Data/UC_Data/T-Ramp/";


	  t_of=1;
	  t_pc=0;
	  for(int nunit = 0; nunit < 7; nunit++){
	    switch(nunit){
	    case 0:
	      path = "../Data/UC_Data/T-Ramp/10";
	      break;
	    case 1:
	      path = "../Data/UC_Data/T-Ramp/20";
	      break;
	    case 2:
	      path = "../Data/UC_Data/T-Ramp/50";
	      break;
	    case 3:
	      path = "../Data/UC_Data/T-Ramp/75";
	      break;
	    case 4:
	      path = "../Data/UC_Data/T-Ramp/100";
	      break;
	    case 5:
	      path = "../Data/UC_Data/T-Ramp/150";
	      break;
	    case 6:
	      path = "../Data/UC_Data/T-Ramp/200";
	      break;
	    }

	    for (auto & p : fs::directory_iterator(path)){


	      cout<<endl <<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
	      cout<<endl << p;
	      cout<<endl <<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;
		
	      fstream log;
	      log.open ("./log.txt", fstream::app);

	      log<< p << " : " << endl;

	      //for(int s = 0 ; s < 3; s++){ //selection of model formulation

	      dom = 0;
	      std::string instance;
	      instance = (p.path().string());
	      //instance = "../Data/UC_Data/T-Ramp/10/10_0_5_w.mod";
	      // Opening the instance file
	      ifstream inStream (instance, ios::binary);
	      // Give to myStreamPosition the position of the input sequence
	      ios::pos_type myStreamPosition = inStream.tellg();
	      // Set the proper position for the stream
	      inStream.seekg (myStreamPosition);

	      UCBlock ucblock;
	      ucblock.instance(inStream);     
	      cout<<endl<<"domanda: " << ucblock.get_units_size() ;
	      cout<<endl<<" // " ;
	      for(int i = 0 ; i < ucblock.get_units_size() +1 ; i++){
	
	
		ucblock.get_nested_Blocks()[i]->generate_static_constraints();
		 
	      }
	      cout<<endl<<"domanda_2:  " << ucblock.get_units_size() ;
	      cout<<endl<<" // " ;
	      double lb2;
	      double lb3;
	      int status;
		
	      //1-root_node without cuts
	      MILPSolver milpsol;  
	      //cout<<endl<<"test_before solver_of = " << ucblock.get_num_active_var << endl;
	      milpsol.set_Block( &ucblock);
	      milpsol.set_par(Solver::kMaxTime,14400.0); //setting the time limit 
	      milpsol.set_par(Solver::kLogVerb,1);   //setting CPLEX Verbosity
	      milpsol.set_par(Solver::kMaxIter, (long) (0.0) ); //initially equal to zero for rood nooodes
	      milpsol.set_par(Solver::kLastAlgPar +CPX_PARAM_PREIND,1);//set on the preprocessing
	      milpsol.set_par(Solver::kLastAlgPar +CPX_PARAM_STARTALG, 4 ); //set LP Solver for root node
	      milpsol.set_par(Solver::kLastAlgPar +CPX_PARAM_THREADS, 1 ); //set number of processors to 1

	      //log<<endl <<"For Solver no." << j <<"Calculating the root node lower bound without CPLEX cuts..." << endl;
	      cout<<endl<<"For Solver no."  << " && Model no. " << dom <<" Calculating the root node lower bound without CPLEX cuts..." << endl;
	      //setting the parameters
	      milpsol.set_par(Solver::kLastAlgPar + CPXPARAM_MIP_Limits_EachCutLimit,0); //set initially no cuts
	      milpsol.set_par(Solver::kLastAlgPar +CPXPARAM_MIP_Cuts_Gomory,-1); //set initially no cuts
		
	      //solving and getting infos
	      std::chrono::time_point<std::chrono::system_clock> t0a, t1a;
	      t0a = std::chrono::system_clock::now();
	      milpsol.solve();
	      t1a= std::chrono::system_clock::now();
	      std::chrono::duration<double> t_final3 = t1a-t0a;		
	      lb3 = milpsol.get_lb();
	      status = milpsol.sol_status;
	      /*log << " lb = " <<  lb3 << " && status = " << status << "time = " << t_final2.count() <<endl;
		cout << " lb = " << lb3 << " && status = " << status << "time = " << t_final2.count() <<endl;*/


	      //2-root_node with cuts
	      //log<<endl <<"For Solver no." << j <<"Calculating the root node lower bound with CPLEX cuts = " << endl;
	      cout<<endl<<"For Solver no." << " && Model no. " << dom <<" Calculating the root node lower bound with CPLEX cuts = " << endl;
	      //setting the parameters
	      milpsol.set_par(Solver::kLastAlgPar + CPXPARAM_MIP_Limits_EachCutLimit,2100000000); //set initially no cuts
	      milpsol.set_par(Solver::kLastAlgPar +CPXPARAM_MIP_Cuts_Gomory,0); //set initially no cuts

	      //solving and getting infos
	      std::chrono::time_point<std::chrono::system_clock> t2, t3;
	      t2= std::chrono::system_clock::now();
	      milpsol.solve();
	      t3= std::chrono::system_clock::now();
	      std::chrono::duration<double> t_final2 = t3-t2;
	      lb2 = milpsol.get_lb();
	      status = milpsol.sol_status;
	      /*log  << " lb = " << lb2 << " && status = " << status << "time = " << t_final3.count() <<endl;
		cout << " lb = " << lb2 << " && status = " << status << "time = " << t_final3.count() <<endl;
		
		log<<" no_cuts_sol_time = "   <<   t_final.count()       <<        " cuts_sol_time = "   << t_final2.count() 
		<< " root_lb_with_cuts = " << lb2 << " root_lb_no_cuts = " << lb3  << endl;

		cout<<" no_cuts_sol_time = "   <<   t_final.count()      <<        " cuts_sol_time = "   << t_final2.count() 
		<< " root_lb_with_cuts = " << lb2 << " root_lb_no_cuts = " << lb3  << endl;*/
	      //solving instance
		
		
	      std::chrono::time_point<std::chrono::system_clock> t0, t1;
	      //log <<endl<<"Solving the UC via CPLEX MIQP Solver..." << endl;
	      cout<<endl<<"For Solver no." << " && Model no. " << dom <<" Solving the UC via CPLEX MILP Solver..." << endl;
	      //setting the parameters
	      milpsol.set_par(Solver::kMaxIter, (long) (9223372036800000000.0) ); //set initially no cuts
	      //solveing and getting infos
	      t0 = std::chrono::system_clock::now();
	      milpsol.solve();
	      t1= std::chrono::system_clock::now();
	      std::chrono::duration<double> t_final = t1-t0;
	      status = milpsol.sol_status;

	      double value;
	      if(t_of==1){
		LinearObjectiveFunction * ucblock_of = boost::any_cast< LinearObjectiveFunction * >(ucblock.get_objective_function());
		value= ucblock_of->value();	
	      }
	      else if(t_of==0) {
		DQuadObjectiveFunction * ucblock_of = boost::any_cast< DQuadObjectiveFunction * >(ucblock.get_objective_function());
		value= ucblock_of->value();
	      }

		
		
	      log<<"MILP Model no. " << dom << " && LP Solver "  << " : sol_time = "   <<t_final.count() <<" sol_cost = " << value 
		 << " sol_status = " << status << " #_nodes = " << milpsol.nodes 
		 << " root_lb_with_cuts = " << lb2 << " ( " << t_final2.count() << " ) " << " root_lb_no_cuts = " << lb3 << " ( " << t_final3.count() << " ) "  << endl;

	      cout<<"MILP Model no." << dom << " && LP Solver "  << " : sol_time = "   <<t_final.count() <<" sol_cost = " << value 
		  << " sol_status = " << status << " #_nodes = " << milpsol.nodes 
		  << " root_lb_with_cuts = " << lb2 << " ( " << t_final2.count() << " ) " << " root_lb_no_cuts = " << lb3 << " ( " << t_final3.count() << " ) "  << endl; 
		
	      //}//for loop of different models	

	      log.close();
	
	    }//for loop of instances

	  }

	}//else for whole UC

	else {
	  /*****Working and Testing 1UC with lang costs for DP+MILPSolvers*******/


	  std::string path = "../Data/1UC_Data/24/"
	    /*
	      /home/kostas/Data_UC/Lang_Data/24
	      /home/kostas/Data_UC/Lang_Data/96
	      /home/kostas/Data_UC/Lang_Data/168

	      /home/kostas/Data_UC/fake_4
	    */;

	  fstream log;
	  log.open ("./log.txt", fstream::app);

    for (auto & p : fs::directory_iterator(path)){


	cout<<endl <<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
	cout<<endl <<"~~~~~~~~~instance: " << p  << "~~~~~~~~~~ " ; //&& t = " << uno.get_t();
	cout<<endl <<"~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"<<endl;

log<< endl << p << " : " << endl;
for(int s = 0 ; s < 1 ; s++){
	for(int j = 0 ; j < 2 ; j++){
	if(j==0) cout<<endl<<"Solving original MIQP model:";
	if(j==1) cout<<endl<<"Solving P/C model:";

	t_of=1;
	t_pc=0;

	dom = 2;
	std::string instance;
	instance = (p.path().string());
	//instance = "/home/tiziano/Scrivania/UC_K/UC/Data/1UC_Data/24/S12ramp10_24.dat";
	// Opening the instance file
	ifstream inStream (instance, ios::binary);
	// Give to myStreamPosition the position of the input sequence
	ios::pos_type myStreamPosition = inStream.tellg();
	// Set the proper position for the stream
	inStream.seekg (myStreamPosition);

	UCBlock uno;
	
	uno.instance(inStream);

	AcadThermalUnitBlock * acd_unit =  dynamic_cast <AcadThermalUnitBlock * > (uno.get_nested_Blocks()[0]);

	/*acd_uc*/uno.get_nested_Blocks()[0]->generate_static_constraints();


	//MILPSolver

	
	MILPSolver milpsol;  
	milpsol.set_Block(uno.get_nested_Blocks()[0]);
	double lb2,lb3;
	int status;


	milpsol.set_par(Solver::kMaxTime,3600.0); //setting the time limit 
	milpsol.set_par(Solver::kLogVerb,1);   //setting CPLEX Verbosity
	milpsol.set_par(Solver::kMaxIter, (long) (0.0) ); //initially equal to zero for rood nooodes

	milpsol.set_par(Solver::kLastAlgPar +CPX_PARAM_PREIND,1);//set on the preprocessing
	milpsol.set_par(Solver::kLastAlgPar +CPX_PARAM_STARTALG, 4 ); //set LP Solver for root node
	milpsol.set_par(Solver::kLastAlgPar +CPX_PARAM_THREADS, 1 ); //set number of processors to 1


	//lb root_node without cplex cuts
	milpsol.set_par(Solver::kLastAlgPar + CPXPARAM_MIP_Limits_EachCutLimit,0); //set initially no cuts
	milpsol.set_par(Solver::kLastAlgPar +CPXPARAM_MIP_Cuts_Gomory,-1); //set initially no cuts
	//solving and getting infos
	status = milpsol.solve();		
	lb3    = milpsol.get_lb();
	cout<<endl<<"For model no. " << acd_unit->model << " sol_status_root_no_cuts = " << status;
	//lb root_node with cplex cuts
	milpsol.set_par(Solver::kLastAlgPar + CPXPARAM_MIP_Limits_EachCutLimit,2100000000); //set initially no cuts
	milpsol.set_par(Solver::kLastAlgPar +CPXPARAM_MIP_Cuts_Gomory,0); //set initially no cuts
	//solving and getting infos
	status = milpsol.solve();		
	lb2    = milpsol.get_lb();
	cout<<endl<<"For model no. " << acd_unit->model << " sol_status_root_con_cuts = " << status;
	//solving the instance now
	milpsol.set_par(Solver::kMaxIter, (long) (9223372036800000000.0) ); //set initially no cuts

	std::chrono::time_point<std::chrono::system_clock> t0, t1;
	t0 = std::chrono::system_clock::now();
	status = milpsol.solve();
	t1= std::chrono::system_clock::now();
	std::chrono::duration<double> t1_fin = t1-t0;

	double value=0;
	if(acd_unit->of==1){
		LinearObjectiveFunction * obj_fun = boost::any_cast< LinearObjectiveFunction * >(uno.get_nested_Blocks()[0]->get_objective_function());
		value=obj_fun->value();
	
	
	}
	else if(acd_unit->of==0) {
	 	DQuadObjectiveFunction * obj_fun = boost::any_cast< DQuadObjectiveFunction * >(uno.get_nested_Blocks()[0]->get_objective_function());
		value=obj_fun->value();
	}


	//	cout<<endl <<"MILP Model no. " << acd_unit->model << " : sol_time = " <<t1_fin.count() <<" sol_cost = " << value << " sol_status = " << status << " #_nodes = " << milpsol.nodes 
	//       << " root_lb_with_cuts = " << lb2 << " root_lb_no_cuts = " << lb3 ;

	//log<<"MILP Model no. " << acd_unit->model << " sol_time = " <<t1_fin.count() <<" sol_cost = " << value << " sol_status = " << status << " #_nodes = " << milpsol.nodes 
	//       << " root_lb_with_cuts = " << lb2 << " root_lb_no_cuts = " << lb3  << endl;
	/*
if (acd_unit->model >= 0){

	cout<<endl <<" schedule: [" ;
	for(int i = 0 ; i < acd_unit->get_total_P() ; i++) cout<< " " << acd_unit->get_P(i)->get_value() << " ; ";
        cout << "]";

}
	*/

}// for pc & of

}//for model

}// for (auto & p : fs::directory_iterator(path))

log.close();

}//if case for 1UC or UC




	cout << endl << "Finitooo..?"; // prints

 return( 0 );
 }

/*--------------------------------------------------------------------------*/
/*------------------------- End File test.cpp ------------------------------*/
/*--------------------------------------------------------------------------*/
