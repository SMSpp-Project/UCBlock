/*--------------------------------------------------------------------------*/
/*--------------------- File BusNetworkBlock.cpp ---------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the BusNetworkBlock class.
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

#include "boost/bind.hpp"
#include "boost/functional/factory.hpp"

#include "BusNetworkBlock.h"
#include "UCBlock.h"
#include "AcadThermalUnitBlock.h"


/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/
extern int dom;
using namespace std;

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

BusNetworkBlock::BusNetworkBlock(UCBlock * flbock ) : NetWorkBlock(flbock) {
}

/*--------------------------------------------------------------------------*/

BusNetworkBlock::~BusNetworkBlock() {
}

/*--------------------------------------------------------------------------*/

void BusNetworkBlock::load(std::istream& inStream) {

 std::string skip;

UCBlock * f_UC_Block = dynamic_cast <UCBlock *> (get_f_Block());

	// Skip row of 	LoadCurve
	//inStream >> skip;
	// Skip row of 	MinSystemCapacity
	inStream >> skip >> skip;
	 // Skip row of MaxSystemCapacity
	inStream >> skip >> skip;
	// Skip row of	MaxThermalCapacity
	inStream >> skip >> skip;

	 int rows, elemPerRow; //local variables for reading the input data
	    inStream >> skip >> rows >> elemPerRow;

	 if (rows * elemPerRow != f_UC_Block->get_t())
	     throw std::invalid_argument( "Load Curve Dimensions not matching the time horizon" );

	 d.resize(f_UC_Block->get_t()); //set proper size to the demand vector


	 for (int i=0; i<f_UC_Block->get_t(); ++i){
		 inStream >> d[i]; //passing the demand values
	}

	 	 inStream >> skip >> elemPerRow;


	 	if (!(f_UC_Block->get_t()/rows == elemPerRow))
	 		throw std::invalid_argument( "Reserve Dimensions not matching the time horizon" );


	 	std::vector<double> percSpinningReserve(elemPerRow);

	     r.resize(f_UC_Block->get_t()); //set proper size to the reserve power vector
	     //storing the values of reserve vector

	     for (int i=0; i<elemPerRow; ++i)
	          inStream >> percSpinningReserve[i]; //storing the data values
	     //storing the values of reserve vector 

	     for (int i=0; i<f_UC_Block->get_t(); ++i)
	         r[i] = d[i] * (1+percSpinningReserve [i%elemPerRow]);
	         //storing the values of reserve vector
}

/*--------------------------------------------------------------------------*/

void BusNetworkBlock::generate_static_constraints() {

 //Cosntrution of Demand_Constraint of the form: Sum_{i \in I} p_it = D_t
 //Cosntrution of Reserve_Constraint of the form: Sum_{i \in I} p_it > R_t

UCBlock * f_UC_Block = dynamic_cast <UCBlock *> (get_f_Block());

 Demand_Const.resize( f_UC_Block->get_t() );
 Reserve_Const.resize( f_UC_Block->get_t() );
 if(dom==0 || dom == 2){
 for( int i = 0 ; i < f_UC_Block->get_t() ; i++ ) {
  // initilaize the vector of pairs of the variables and coefficients
  // of the corresponding constraints

  LinearConstraint::v_coeff_pair *v_demand_pair = new LinearConstraint::v_coeff_pair( f_UC_Block->get_units_size() );
  LinearConstraint::v_coeff_pair *v_reserve_pair = new LinearConstraint::v_coeff_pair( f_UC_Block->get_units_size() );
  LinearConstraint::v_coeff_pair::iterator it;
  for( int j = 0 ; j < f_UC_Block->get_units_size() ; j++ ) {

	  it = v_demand_pair->begin() + j;

       *it = LinearConstraint::coeff_pair((f_UC_Block->get_units()[ j ]->get_P(i)), 1.0);
	   // store p_it of the Sum
		AcadThermalUnitBlock * Acad_Unit = dynamic_cast <AcadThermalUnitBlock *> (f_UC_Block->get_units()[ j ]);
       it = v_reserve_pair->begin() + j;
       *it = LinearConstraint::coeff_pair((f_UC_Block->get_units()[ j ]->get_U(i)) ,Acad_Unit->get_fMaxPower() );
//cout<<endl<<" fMaxPower_reserve = " << Acad_Unit->get_fMaxPower();
	   // store p_it of the Sum
   }
  
  Demand_Const[i].add_variables(v_demand_pair , false );
  Demand_Const[i].set_Block(this);
  Demand_Const[i].set_lhs(d[ i ]);
  Demand_Const[i].set_rhs(d[ i ]);

  Reserve_Const[i].add_variables(v_reserve_pair , false ); 
  Reserve_Const[i].set_Block(this);
  Reserve_Const[i].set_lhs(r[ i ]);
  Reserve_Const[i].set_rhs(Inf<double>());
 
  
  }
}//if model 0/2

 else if(dom==1){
int u_count = 0;
for( int i = 0 ; i < f_UC_Block->get_t() ; i++ ) {
  // initilaize the vector of pairs of the variables and coefficients
  // of the corresponding constraints
	u_count = 0;
        for( int z = 0 ; z < f_UC_Block->get_units_size() ; z++ ) {
		AcadThermalUnitBlock * t_unit = dynamic_cast <AcadThermalUnitBlock *> (f_UC_Block->get_units()[ z ]);

    	   	for(int j = 0 ; j < t_unit->Y_plus_pair.size(); j++){
    	    		if (t_unit->Y_plus_pair[j].first <= i + 1 && i + 1 <= t_unit->Y_plus_pair[j].second ){
        	        	u_count++;
        	    	}       
        	} // for j-loop
	}// for z-loop

  	LinearConstraint::v_coeff_pair *v_demand_pair = new LinearConstraint::v_coeff_pair( u_count );
  	LinearConstraint::v_coeff_pair *v_reserve_pair = new LinearConstraint::v_coeff_pair( u_count );
  	LinearConstraint::v_coeff_pair::iterator it_d = v_demand_pair->begin();
  	LinearConstraint::v_coeff_pair::iterator it_r = v_reserve_pair->begin();

        for( int z = 0 ; z < f_UC_Block->get_units_size() ; z++ ) {
		AcadThermalUnitBlock * t_unit = dynamic_cast <AcadThermalUnitBlock *> (f_UC_Block->get_units()[ z ]);
	        for(int j = 0 ; j < t_unit->Y_plus_pair.size(); j++){
            		if (t_unit->Y_plus_pair[j].first <= i + 1 && i + 1 <= t_unit->Y_plus_pair[j].second ){
                		
		                *it_r = LinearConstraint::coeff_pair( &t_unit->Y_plus[j] , t_unit->get_fMaxPower()); 
				it_r++;
	
                		int k = i+1 - t_unit->Y_plus_pair[j].first;
		                if (t_unit->Y_plus_pair[j].first == 0 ) k--;
		                *it_d = LinearConstraint::coeff_pair( &t_unit->P_hk[t_unit->v_p_size[j]-t_unit->v_p_ind[j]+k], 1.0);                                             
				it_d++;
	               }//if
	        } // for j
	}// for z-loop

  Demand_Const[i].add_variables(v_demand_pair , false );
  Demand_Const[i].set_Block(this);
  Demand_Const[i].set_lhs(d[ i ]);
  Demand_Const[i].set_rhs(d[ i ]);

  Reserve_Const[i].add_variables(v_reserve_pair , false ); 
  Reserve_Const[i].set_Block(this);
  Reserve_Const[i].set_lhs(r[ i ]);
  Reserve_Const[i].set_rhs(Inf<double>());
  }

 }
 else if(dom==3){
  int p_count=0;
  int p_ind=0;
  int ind = 0;

  int u_count = 0;

for( int i = 0 ; i < f_UC_Block->get_t() ; i++ ) {
  // initilaize the vector of pairs of the variables and coefficients
  // of the corresponding constraints
    p_count=0;
    for( int z = 0 ; z < f_UC_Block->get_units_size() ; z++ ) {
		AcadThermalUnitBlock * t_unit = dynamic_cast <AcadThermalUnitBlock *> (f_UC_Block->get_units()[ z ]);
		ind = 0;
		for (int j = 0 ; j < t_unit->v_ph_num.size(); j++){
			ind = t_unit->Y_plus_pair[t_unit->v_ph_ind[j]].second;
			if (ind == f_UC_Block->get_t() +1) ind--;
			if ( i+1 <= ind ){
				for (int h = 0 ; h < t_unit->v_ph_size[j] ; h++){
					if ( i+1 == ind-h ) p_count++;
				} //for h-loop
			} // if

		} //for j-loop
	} //for z-loop


  LinearConstraint::v_coeff_pair *v_demand_pair  = new LinearConstraint::v_coeff_pair( p_count );
  LinearConstraint::v_coeff_pair::iterator it = v_demand_pair->begin();
  for( int z = 0 ; z < f_UC_Block->get_units_size() ; z++ ) {
	AcadThermalUnitBlock * t_unit = dynamic_cast <AcadThermalUnitBlock *> (f_UC_Block->get_units()[ z ]);
	p_ind=0;
	for (int j = 0 ; j < t_unit->v_ph_num.size(); j++){
		ind = t_unit->Y_plus_pair[t_unit->v_ph_ind[j]].second;
		if (ind == f_UC_Block->get_t() +1) ind--;
			if ( i+1 <= ind ){
				for (int h = 0 ; h < t_unit->v_ph_size[j] ; h++){
					if ( i+1 == ind - t_unit->v_ph_size[j] + 1 + h ){
						p_ind = t_unit->v_ph_tot[j] - t_unit->v_ph_size[j] + h;
			                		*it = LinearConstraint::coeff_pair( &t_unit->P_h[ p_ind ] , 1.0); //store		
							it++;
					} // if
				} //for h-loop
			} // if
	} //for j-loop
   }//for z-loop

  Demand_Const[i].add_variables(v_demand_pair , false );
  Demand_Const[i].set_Block(this);
  Demand_Const[i].set_lhs(d[ i ]);
  Demand_Const[i].set_rhs(d[ i ]);

  u_count = 0;
  for( int z = 0 ; z < f_UC_Block->get_units_size() ; z++ ) {
	AcadThermalUnitBlock * t_unit = dynamic_cast <AcadThermalUnitBlock *> (f_UC_Block->get_units()[ z ]);
        	for(int j = 0 ; j < t_unit->Y_plus_pair.size(); j++){
        	    if (t_unit->Y_plus_pair[j].first <= i + 1 && i + 1 <= t_unit->Y_plus_pair[j].second ){
        	        u_count++;
        	    }       
        	} // for j
  } // for z

  LinearConstraint::v_coeff_pair *v_reserve_pair = new LinearConstraint::v_coeff_pair( u_count );

  it = v_reserve_pair->begin();

  for( int z = 0 ; z < f_UC_Block->get_units_size() ; z++ ) {
	AcadThermalUnitBlock * t_unit = dynamic_cast <AcadThermalUnitBlock *> (f_UC_Block->get_units()[ z ]);
	  for(int j = 0 ; j < t_unit->Y_plus_pair.size(); j++){
            		if (t_unit->Y_plus_pair[j].first <= i + 1 && i + 1 <= t_unit->Y_plus_pair[j].second ){
				  *it = LinearConstraint::coeff_pair(&t_unit->Y_plus[j] ,t_unit->get_fMaxPower() );
				   it++;
			}
	  }//for j-loop
  }//for z-loop


  Reserve_Const[i].add_variables(v_reserve_pair , false ); 
  Reserve_Const[i].set_Block(this);
  Reserve_Const[i].set_lhs(r[ i ]);
  Reserve_Const[i].set_rhs(Inf<double>());
  }

 }
 // Adding them to the Static Constraints Vector
 add_static_constraint( Demand_Const );
 add_static_constraint( Reserve_Const );

}


BusNetworkBlock::_init::_init()
{
NetWorkBlock::f_factory()["LoadCurve"] =
 boost::bind(boost::factory<BusNetworkBlock*>(), _1);
 }


 BusNetworkBlock::_init BusNetworkBlock::_initializer;

 /* Ensure that the static member _initializer is initialized, so that
    BusNetworkBlock is registered into NetWorkBlock::f_factory. */



