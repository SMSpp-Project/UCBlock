/*--------------------------------------------------------------------------*/
/*----------------------- File UCBlock.cpp ---------------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the UCBlock class.
 *
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
/*---------------------------- IMPLEMENTATION ------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------------------ INCLUDES ----------------------------------*/
/*--------------------------------------------------------------------------*/

#include <iostream>
#include <vector>
#include "UCBlock.h"
#include "UnitBlock.h"
#include "NetWorkBlock.h"


/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/
//int datas = 0;
extern int choice;
extern int t_of;
using namespace std;
using namespace SMSpp_di_unipi_it; 


/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/
UCBlock::~UCBlock(){ 
  
    if( ! Units.empty() )
      //std::for_each( Units.begin(), Units.end(), []( UnitBlock* p ) { delete p; } );
        for ( auto p : Units ) delete p;
   
    if(choice == 0 && Network != nullptr)
     delete Network;
    

} 

void UCBlock::instance(std::istream& inStream) {

 if(choice == 0){
 	std::string skip;
 	std::string thermal, network;

 	inStream >> skip >> skip>> skip >> t; //setting the no. of timesteps
 	inStream >> thermal >> units_size; //setting the no. of thermal units and the string of therml

 	inStream >> skip >> skip; // skipping hydro units (if any)
 	inStream >> skip >> skip; // skipping hydro cascades (if any)

 	v_Block.resize( units_size + 1 );
 	inStream >> network; //setting the string of network the corresponding factory

 	Units.resize( units_size );


 	for( int i = 0 ; i < units_size ; i++) { 
        	/* Note that we initialize first the Units since we need to initialize the Variable U and P,
         	* since they are used in the Constraints of the BusNetwork and can not be included propely 
         	*without being initialized in advance */

    	Units[ i ] = UnitBlock::U_factory()[thermal](this); // set the corresponding object via factory

   	}



 	Network =  NetWorkBlock::f_factory()[network](this);//[network]; //new BusNetworkBlock();
 	// set the corresponding object via factory

 	Network->load( inStream );
 	//pass the data to the corresponding Network Block

  	v_Block[ units_size ] = Network;
 	//add the Network to the vector of nested Blocks

 
  	inStream >> skip;
 	for( int i = 0 ; i < units_size ; i++) {
  		Units[ i ]->load( inStream );
  		// pass the data to the corresponding thermal units
 
 		v_Block[ i ] = (Units[ i ]);
 		 // add each thermal unit Block to the vector of nested Blocks
  	}
	if ( t_of == 0 ){
		q_of.set_type(ObjectiveFunction::eMin);
		set_objective_function(q_of); // add quad objective funtion to the Block
	}
	else if ( t_of == 1 ){
		l_of.set_type(ObjectiveFunction::eMin);
		set_objective_function(l_of); // add linear objective funtion to the Block
	}
}



else if (choice == 1){

 std::string skip;
 std::string thermal;

 inStream >> skip >> t; //setting the no. of timesteps

	
 v_Block.resize( 1 );
 Units.resize( 1 );
 thermal = "NumThermal";


 Units[ 0 ] = UnitBlock::U_factory()[thermal](this); // set the corresponding object via factory
 Units[ 0 ]->load( inStream );
  // pass the data to the corresponding thermal units
  v_Block[ 0 ] = (Units[ 0 ]);

}


}
