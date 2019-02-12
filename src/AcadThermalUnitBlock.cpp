/*--------------------------------------------------------------------------*/
/*-------------------- File AcadThermalUnitBlock.cpp -----------------------*/
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
#include <iostream>
#include <random>
#include "boost/bind.hpp"
#include "boost/functional/factory.hpp"

#include "AcadThermalUnitBlock.h"
#include "UCBlock.h"
//extern int solo;
extern int dom;
extern int choice;
extern int t_of;
extern int t_pc;

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

using namespace std;

/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/

void AcadThermalUnitBlock::load(std::istream& inStream){

UCBlock * f_UC_Block = dynamic_cast <UCBlock *> (get_f_Block());

    string skip;
    data = choice;
    model = dom;
    pc=t_pc;
    of=t_of;
    //reading the input data from original files without duals
	if (data == 0){
    inStream >> skip >> fQuadTherm >> fLinearTherm >> fConstTherm >> fMinPower >> fMaxPower >> fInitUpDownTime >> fMinUpTime >> fMinDownTime
	>> coolAndFuelCost >> hotAndFuelCost >> tau >> skip >> fixedCost >> skip >> fInitPower;



    //check for rampconstraints
    const string kToken("RampConstraints");
    ios::pos_type myStreamPosition = inStream.tellg();
    string myToken;

    inStream >> myToken;

    if (myToken == kToken)
    {
    	rampconst = true;
        inStream >> fMaxRampUp >> fMaxRampDown;
    }
    else{
    	rampconst = false;
        inStream.seekg(myStreamPosition);
    }

fBoundDown = fBoundOn = fMinPower;

lambda.resize(f_UC_Block->get_t());

for(int i = 0 ; i < f_UC_Block->get_t() ; ++i)
/*inStream >> */ lambda[i] = 0;

mew.resize(f_UC_Block->get_t());
//inStream >>	skip;
for(int i = 0 ; i < f_UC_Block->get_t() ; ++i)
/*inStream >>*/ mew[i] = 0;

	}

//reading the input data from files with duals
	else{
	
   	rampconst = true;

	inStream >> 
    skip >> fQuadTherm >> skip >> fLinearTherm >> skip >> fConstTherm >> 
    skip >> fMinPower >> skip >> fMaxPower >> 
    skip >> fInitUpDownTime >> skip >> fMinUpTime >> skip >> fMinDownTime >> 
    skip >> f_start_cost >>
    skip >> fInitPower >>
    skip >> fMaxRampUp >> skip >> fMaxRampDown >>
	skip >> fBoundOn >> skip >> fBoundDown;

   
lambda.resize(f_UC_Block->get_t());
inStream >>	skip;
for(int i = 0 ; i < f_UC_Block->get_t() ; ++i)
inStream >>  lambda[i];
mew.resize(f_UC_Block->get_t());
inStream >>	skip;
for(int i = 0 ; i < f_UC_Block->get_t() ; ++i)
inStream >> mew[i];



}//else for reading instances

if(model == 0){

  if ( fInitUpDownTime > 0 )
    {
    
      init_t =  ( fInitUpDownTime >= fMinUpTime ? 0 : fMinUpTime   - fInitUpDownTime );  
    }
  else
    {
      init_t = ( -fInitUpDownTime >= fMinDownTime ? 0 : fMinDownTime + fInitUpDownTime ); 
    
    } 

//commented out for dyn formulation
U_10.resize(f_UC_Block->get_t() - init_t);
U_11.resize(f_UC_Block->get_t() - init_t);
U_00.resize(f_UC_Block->get_t() - init_t);
U_01.resize(f_UC_Block->get_t() - init_t);

 for(int i = 0 ; i<f_UC_Block->get_t() - init_t ; i++){
		 //set pointer of the father Block
		 U_10[i].set_Block(this); 
		 U_11[i].set_Block(this); 
		 U_01[i].set_Block(this); 
		 U_00[i].set_Block(this); 
		 //set lower and upper bounds
		 U_10[i].set_lb(0.0); 
		 U_11[i].set_lb(0.0); 
		 U_01[i].set_lb(0.0); 
		 U_00[i].set_lb(0.0); 
 		 U_10[i].set_ub(1.0); 
		 U_11[i].set_ub(1.0); 
		 U_01[i].set_ub(1.0); 
		 U_00[i].set_ub(1.0); 
		 //set type of variables
		 U_10[i].set_type(ColVariable::binary);
		 U_11[i].set_type(ColVariable::binary);  
		 U_01[i].set_type(ColVariable::binary);
		 U_10[i].set_type(ColVariable::binary);
		 
	 }

add_static_variable(U_11);  
add_static_variable(U_10);
add_static_variable(U_00);    
add_static_variable(U_01);

if (pc == 1){

	p_k.resize(2);
	p_k[0]= fMinPower;
	p_k[1]= fMaxPower;

	z_t.resize( f_UC_Block->get_t() );
        for(int i = 0 ; i<f_UC_Block->get_t(); i++){
		 z_t[i].set_Block(this); 
		 z_t[i].set_lb(-Inf<double>()); 
 		 z_t[i].set_ub(Inf<double>()); 
		 z_t[i].set_type(ColVariable::continuous);
	}
add_static_variable(z_t);  
}


}
else {

  if ( fInitUpDownTime > 0 )
    {
    
      init_t =  ( fInitUpDownTime >= fMinUpTime ? 0 : fMinUpTime   - fInitUpDownTime ); 
    }
  else
    {
      init_t = ( -fInitUpDownTime >= fMinDownTime ? 1 : fMinDownTime + fInitUpDownTime + 1 ); 
    
    } 

int count_plus=0;
int count_minus=0;


//constructing all arc connections frmo the starting node

if (fInitUpDownTime > 0 ){

    /*start with the arc connections from start node to all other arcs, where we
      have the following two scenarios:
      i)  in case init_t = 0 this means that from the start node we have arcs y^+_s,k
          going towards all the offline nodes and the sink nodes 
      ii) in case init_t != 0 this means that from the start node we have arcs y^+_s,k
          going towards all offline nodes for k > init_t and the sink nodes 
    */
    for(int i = init_t ; i <= f_UC_Block->get_t()+1 ; i++){
	if(i==0) i++; //added
        Y_plus_pair.push_back(make_pair(0,i));
        count_plus++;
    }   

   //proceed with constructing all the arcs that are connecting the nodes of the graph

   /* 
       Construcing the y^- arcs that we have the following:
       i)   all offline nodes that are prior to init_t are excluded (condition only valid,
            in case there it init_t != 0), where no y^- arcs are considered
       ii)  all offline nodes that are greater than init_t but yet again smaller than T-\tau^-
            for which there are arcs to all the online nodes t+\tau^- and the sink node
       iii) all offline nodes that are greater than T-\tau^-, for which there are arcs 
            connecting them only with the sink node
    */
    for(int i = init_t ; i <= f_UC_Block->get_t() /*- fMinDownTime*/; i++){
        //construct all y- arcs starting from offline node i
        if( i < f_UC_Block->get_t() - fMinDownTime - 1 ){
            //the case where the node has arc connecting it with online nodes
            for(int j = i + fMinDownTime + 1 ; j <= f_UC_Block->get_t() ; j++){
            Y_minus_pair.push_back(make_pair(i,j));
            count_minus++;
            }
        }
       // connection of offline arc with sink node
       Y_minus_pair.push_back(make_pair(i,f_UC_Block->get_t()+1));
       count_minus++;
       v_minus.push_back(i);
    }

   /* 
       Construcing the y^+ arcs that we have the following:
       i)   all online nodes that are prior to init_t + \tau^- are excluded (condition only valid,
            in case there it init_t != 0), where no y^+ arcs are considered
       ii)  all online nodes that are greater than init_t + \tau^- but yet again smaller than T-\tau^+
            for which there are arcs to all the offline nodes t+\tau^+ and the sink node
       iii) all online nodes that are greater than T-\tau^+, for which there are arcs 
            connecting them only with the sink node
    */
    for(int i = init_t+fMinDownTime+1 ; i <= f_UC_Block->get_t() ; i++){
        //construct all y+ arcs starting from online node i
        if( i < f_UC_Block->get_t() - fMinUpTime + 1 ){
            //the case where the node has arc connecting it with online nodes
            for(int j = i + fMinUpTime-1 ; j <= f_UC_Block->get_t() ; j++){
            Y_plus_pair.push_back(make_pair(i,j));
            count_plus++;
            }
        }
       // connection of offline arc with sink node
       Y_plus_pair.push_back(make_pair(i,f_UC_Block->get_t()+1));
       v_plus.push_back(i);
       count_plus++;
    }


   
}

else if ( fInitUpDownTime < 0 ) {

    /*start with the arc connections from start node to all other arcs, where we
      have the following two scenarios:
      i)  in case init_t = 0 this means that from the start node we have arcs y^-_s,k
          going towards all the online nodes and the sink nodes 
      ii) in case init_t != 0 this means that from the start node we have arcs y^-_s,k
          going towards all the online nodes for k > init_t and the sink nodes 
    */
    for(int i = init_t ; i <= f_UC_Block->get_t()+1 ; i++){
        Y_minus_pair.push_back(make_pair(0,i));
        count_minus++;
        }

  //proceed with constructing all the arcs that are connecting the nodes of the graph

   /* 
       Construcing the y^+ arcs that we have the following:
       i)   all online nodes that are prior to init_t are excluded (condition only valid,
            in case there it init_t != 0), where no y^+ arcs are considered
       ii)  all online nodes that are greater than init_t but yet again smaller than T-\tau^+
            for which there are arcs to all the offline nodes t+\tau^+ and the sink node
       iii) all online nodes that are greater than T-\tau^+, for which there are arcs 
            connecting them only with the sink node
    */
    for(int i = init_t ; i <= f_UC_Block->get_t() ; i++){
        //construct all y+ arcs starting from online node i
        if( i <= f_UC_Block->get_t() - fMinUpTime + 1 ){
            //the case where the node has arc connecting it with online nodes
            for(int j = i + fMinUpTime - 1 ; j <= f_UC_Block->get_t() ; j++){
            Y_plus_pair.push_back(make_pair(i,j));
            count_plus++;
            }
        }
       // connection of offline arc with sink node
       Y_plus_pair.push_back(make_pair(i,f_UC_Block->get_t()+1));
       count_plus++;
       v_plus.push_back(i);
    }

   /* 
       Construcing the y^- arcs that we have the following:
       i)   all offline nodes that are prior to init_t + \tau^+ are excluded (condition only valid,
            in case there it init_t != 0), where no y^- arcs are considered
       ii)  all offline nodes that are greater than init_t + \tau^+ but yet again smaller than T-\tau^-
            for which there are arcs to all the online nodes t+\tau^- and the sink node
       iii) all offline nodes that are greater than T-\tau^-, for which there are arcs 
            connecting them only with the sink node
    */
    for(int i = init_t+fMinUpTime - 1 ; i <= f_UC_Block->get_t() ; i++){
        //construct all y+ arcs starting from online node i
        if( i <= f_UC_Block->get_t() - fMinDownTime - 1){
            //the case where the node has arc connecting it with online nodes
            for(int j = i + fMinDownTime +1 ; j <= f_UC_Block->get_t() ; j++){
            Y_minus_pair.push_back(make_pair(i,j));
            count_minus++;
            }
        }
       // connection of offline arc with sink node
       Y_minus_pair.push_back(make_pair(i,f_UC_Block->get_t()+1));
       count_minus++;
       v_minus.push_back(i);
    }



 }

Y_minus.resize(count_minus);
Y_plus.resize(count_plus);



 for(int i = 0 ; i<Y_minus.size() ; i++){
		 //set pointer of the father Block
		 Y_minus[i].set_Block(this); 
		 //set lower and upper bounds
		 Y_minus[i].set_lb(0.0); 
 		 Y_minus[i].set_ub(1.0); 
		 //set type of variables
		 Y_minus[i].set_type(ColVariable::binary);
	 }

 for(int i = 0 ; i<Y_plus.size() ; i++){
		 //set pointer of the father Block
		 Y_plus[i].set_Block(this); 
		 //set lower and upper bounds
		 Y_plus[i].set_lb(0.0); 
 		 Y_plus[i].set_ub(1.0); 
		 //set type of variables
		 Y_plus[i].set_type(ColVariable::binary);
		 if(model == 2){ if (Y_plus_pair[i].first == 0 /*&& Y_plus_pair[i].second > 1*/) {r_init.push_back(i);} }
	 }

add_static_variable(Y_plus);  
add_static_variable(Y_minus);

 if (model == 1){
	v_p_size.resize(count_plus);
	v_p_ind.resize(count_plus);
	v_p_size[0]=Y_plus_pair[0].second - Y_plus_pair[0].first + 1;
	if (Y_plus_pair[0].first == 0) v_p_size[0] = v_p_size[0] - 1;
	if (Y_plus_pair[0].second == f_UC_Block->get_t()+1 ) v_p_size[0] = v_p_size[0] - 1;
	v_p_ind[0] = v_p_size[0];
	for(int i = 1 ; i < Y_plus.size() ; ++i){

		v_p_size[i] = (Y_plus_pair[i].second - Y_plus_pair[i].first + 1) + v_p_size[i-1];
		if (Y_plus_pair[i].first == 0 ) v_p_size[i] = v_p_size[i] - 1;
		if (Y_plus_pair[i].second == f_UC_Block->get_t()+1 ) v_p_size[i] = v_p_size[i] - 1;
		v_p_ind[i] = v_p_size[i] - v_p_size[i-1] ;
	}

	P_hk.resize(v_p_size[Y_plus.size()-1]);

	for(int i = 0 ; i<P_hk.size() ; i++){
		 //set pointer of the father Block
		 P_hk[i].set_Block(this); 
		 //set lower and upper bounds
		 P_hk[i].set_lb(0.0); 
 		 P_hk[i].set_ub(Inf<double>()); 
		 //set type of variables
		 P_hk[i].set_type(ColVariable::continuous);
	 }

	add_static_variable(P_hk); 


	U_Const.resize(f_UC_Block->get_t() );
	P_Const.resize(f_UC_Block->get_t() );

	int u_count;
	for(int i = 0 ; i < f_UC_Block->get_t() ; i++){
        	u_count = 0;
        	for(int j = 0 ; j < Y_plus_pair.size(); j++){
        	    if (Y_plus_pair[j].first <= i + 1 && i + 1 <= Y_plus_pair[j].second ){
        	            u_count++;
        	    }       
        	} // for j

    		LinearConstraint::v_coeff_pair  *v_u_pair = new LinearConstraint::v_coeff_pair(u_count+1);
		LinearConstraint::v_coeff_pair  *v_p_pair = new LinearConstraint::v_coeff_pair(u_count+1);
		//LinearConstraint::v_coeff_pair  *v_p2_pair = new LinearConstraint::v_coeff_pair(u_count+1);
		LinearConstraint::v_coeff_pair::iterator it;

		//Store the ramp down elements
		int k=0;

	        for(int j = 0 ; j < Y_plus_pair.size(); j++){
            		if (Y_plus_pair[j].first <= i + 1 && i + 1 <= Y_plus_pair[j].second ){
                		it = v_u_pair->begin()+k;
		                *it = LinearConstraint::coeff_pair( &Y_plus[j] , 1); //store

                		int z = i+1 - Y_plus_pair[j].first;
		                if (Y_plus_pair[j].first == 0 ) z--;
                		it = v_p_pair->begin()+k;k++;
		                *it = LinearConstraint::coeff_pair( &P_hk[v_p_size[j]-v_p_ind[j]+z], 1);                                             
	               }//if
	        } // for j
	it = v_u_pair->begin()+k;
	*it = LinearConstraint::coeff_pair( &U[i] , -1); //store U_t
	it = v_p_pair->begin()+k;
	*it = LinearConstraint::coeff_pair( &P[i] , -1); //store U_t
   

	U_Const[i].add_variables( v_u_pair, false);
	U_Const[i].set_Block(this);
	U_Const[i].set_lhs(0);
	U_Const[i].set_rhs(0);

	P_Const[i].add_variables( v_p_pair, false);
	P_Const[i].set_Block(this);
	P_Const[i].set_lhs(0);
	P_Const[i].set_rhs(0);

  	}//for i-time steps of U

	add_static_constraint(U_Const);
	add_static_constraint(P_Const);


}//if model == 1

else if (model == 2){

v_t.resize(f_UC_Block->get_t() );
 
int ind = 0;
for(int i = 0; i < v_t.size() ; i++){
    ind = i+1;
    for(int j = 0 ; j < Y_plus_pair.size(); j++){
        if (Y_plus_pair[j].first <= ind && ind <= Y_plus_pair[j].second ){
            v_t[i].index.push_back(j);
            if ( ind + 1 <= Y_plus_pair[j].second && ind != f_UC_Block->get_t() )
                v_t[i].sindex.push_back(j);
	    
	    if (Y_plus_pair[j].first == ind && Y_plus_pair[j].second != ind)
	      v_t[i].b_on.push_back(j);
            if (Y_plus_pair[j].second == ind && Y_plus_pair[j].first != ind)
	      v_t[i].b_dn.push_back(j);
            if (Y_plus_pair[j].second == ind && Y_plus_pair[j].first == ind){
	      if(fBoundDown <= fBoundOn)
		v_t[i].b_dn.push_back(j);
	      else
		v_t[i].b_on.push_back(j);
	    }
	    
            if (Y_plus_pair[j].first != ind &&  Y_plus_pair[j].second != ind)
                v_t[i].p_mx.push_back(j);
    }//first if
      
    }


}

U_Const.resize(f_UC_Block->get_t() );
int u_count;
for(int i = 0 ; i < f_UC_Block->get_t() ; i++){
        u_count = 0;
        for(int j = 0 ; j < Y_plus_pair.size(); j++){
            if (Y_plus_pair[j].first <= i + 1 && i + 1 <= Y_plus_pair[j].second ){
                    u_count++;
            }       
        } // for j

    LinearConstraint::v_coeff_pair  *v_u_pair = new LinearConstraint::v_coeff_pair(u_count+1);
    LinearConstraint::v_coeff_pair::iterator it;

    //Store the ramp down elements
    int k=0;
     for(int j = 0 ; j < Y_plus_pair.size(); j++){
            if (Y_plus_pair[j].first <= i + 1 && i + 1 <= Y_plus_pair[j].second ){
                	it = v_u_pair->begin()+k;k++;
                    *it = LinearConstraint::coeff_pair( &Y_plus[j] , 1); //store                     
            }       
        } // 
    it = v_u_pair->begin()+k;
    *it = LinearConstraint::coeff_pair( &U[i] , -1); //store U_t
   

	U_Const[i].add_variables( v_u_pair, false);
	U_Const[i].set_Block(this);
	U_Const[i].set_lhs(0);
	U_Const[i].set_rhs(0);

  }//for i-time steps of U

add_static_constraint(U_Const);


} //if model  == 2



 if (model == 3){

	int ph_count=0;
	int z=0;
	int k=0;
	for(int i = 1 ; i < Y_plus.size() ; ++i){
		k++;
		z=0;
		if(Y_plus_pair[i-1].first == 0 && Y_plus_pair[i-1].second == f_UC_Block->get_t()+1 ) z--;
		else if(Y_plus_pair[i-1].first == 0 || Y_plus_pair[i-1].second == f_UC_Block->get_t()+1 ) z=0;
		else z++;
		
		if(Y_plus_pair[i].first != Y_plus_pair[i-1].first){
			 ph_count += (Y_plus_pair[i-1].second - Y_plus_pair[i-1].first + z);
			 v_ph_size.push_back(Y_plus_pair[i-1].second - Y_plus_pair[i-1].first + z);
			 v_ph_num.push_back(k);
			 v_ph_ind.push_back(i-1);
			 v_ph_tot.push_back(ph_count);
			 k=0;
			 //cout<<endl<<"ph_count2 =" << ph_count;
			 if (i == Y_plus.size()-1){ //in case we are in the last check and is valid we need to also add the last element
				z = 0;
				if(Y_plus_pair[i].first == 0 && Y_plus_pair[i].second == f_UC_Block->get_t()+1 ) z--;
				else if(Y_plus_pair[i].first == 0 || Y_plus_pair[i].second == f_UC_Block->get_t()+1 ) z=0;
				else z++;
				ph_count += (Y_plus_pair[i].second - Y_plus_pair[i].first + z);
		         	//cout<<endl<<"bike && ph_count =" << ph_count;
			 	v_ph_size.push_back(Y_plus_pair[i].second - Y_plus_pair[i].first + z);
	      		        v_ph_num.push_back(1);
			        v_ph_ind.push_back(i);
			        v_ph_tot.push_back(ph_count);
			}// if last element of y^+
		}//if first/last
	}
	//cout<<endl<<"Init_check: ph_count = " << ph_count;
	P_h.resize(ph_count);

	for(int i = 0 ; i<P_h.size() ; i++){
		 //set pointer of the father Block
		 P_h[i].set_Block(this); 
		 //set lower and upper bounds
		 P_h[i].set_lb(0.0); 
 		 P_h[i].set_ub(Inf<double>()); 
		 //set type of variables
		 P_h[i].set_type(ColVariable::continuous);
	 }

	add_static_variable(P_h); 

	U_Const.resize(f_UC_Block->get_t() );
	P_Const.resize(f_UC_Block->get_t() );
	int u_count;
	int p_count;

	for(int i = 0 ; i < f_UC_Block->get_t() ; i++){
	
		//connecting y_hk with u variables
        	u_count = 0;
        	for(int j = 0 ; j < Y_plus_pair.size(); j++){
        	    if (Y_plus_pair[j].first <= i + 1 && i + 1 <= Y_plus_pair[j].second ){
        	            u_count++;
        	    }       
        	} // for j
		
    		LinearConstraint::v_coeff_pair  *v_u_pair = new LinearConstraint::v_coeff_pair(u_count+1);
		LinearConstraint::v_coeff_pair::iterator it;
		k=0;

	        for(int j = 0 ; j < Y_plus_pair.size(); j++){
            		if (Y_plus_pair[j].first <= i + 1 && i + 1 <= Y_plus_pair[j].second ){
                		it = v_u_pair->begin()+k;k++;
		                *it = LinearConstraint::coeff_pair( &Y_plus[j] , 1); //store
	               }//if
	        } // for j
		it = v_u_pair->begin()+k;
		*it = LinearConstraint::coeff_pair( &U[i] , -1); //store U_t
	
		U_Const[i].add_variables( v_u_pair, false);
		U_Const[i].set_Block(this);
		U_Const[i].set_lhs(0);
		U_Const[i].set_rhs(0);


		//connecting p_h with p variables
		p_count=0;
		int ind = 0;
		for (int j = 0 ; j < v_ph_num.size(); j++){
			ind = Y_plus_pair[v_ph_ind[j]].second;
			if (ind == f_UC_Block->get_t() +1) ind--;
			if ( i+1 <= ind ){
				for (int h = 0 ; h < v_ph_size[j] ; h++){
					if ( i+1 == ind-h ) p_count++;
				} //for h-loop
			} // if

		} //for j-loop

    		LinearConstraint::v_coeff_pair  *v_p_pair = new LinearConstraint::v_coeff_pair(p_count+1);
		auto it_p = v_p_pair->begin();
		//cout<<endl<<"Cheeeck[" << i << "]:  v_p_pair->size() = " << v_p_pair->size();
		//cout<<endl<<"//-//";
		int p_ind=0;
		if(p_count > 0){
			for (int j = 0 ; j < v_ph_num.size(); j++){
				ind = Y_plus_pair[v_ph_ind[j]].second;
				if (ind == f_UC_Block->get_t() +1) ind--;
				if ( i+1 <= ind ){
					for (int h = 0 ; h < v_ph_size[j] ; h++){
						if ( i+1 == ind - v_ph_size[j] + 1 + h ){
							p_ind = v_ph_tot[j] - v_ph_size[j] + h;
							*it_p = LinearConstraint::coeff_pair( &P_h[ p_ind ] , 1); //store		
							it_p++;
						} // if p_h
					} //for h-loop
				} // if

			} //for j-loop
		} // if p_count


		*it_p = LinearConstraint::coeff_pair( &P[i] , -1); //store U_t

		P_Const[i].add_variables( v_p_pair, false);
		P_Const[i].set_Block(this);
		P_Const[i].set_lhs(0);
		P_Const[i].set_rhs(0);
	

  	}//for i-loop

	add_static_constraint(U_Const);
	add_static_constraint(P_Const);

}//if model == 3

 } // if mod > 0


     //Constructing the Objective Function
    obj_function( );


} 




void AcadThermalUnitBlock::generate_static_constraints( void ){

   //Constructing the Pmin & Pmax Constraints
	if (model == 0){
    p_min_max_const(fMinPower, fMaxPower);

    //Constructing the Min Up/Down Constraints
    min_up_down_const (fMinDownTime, fMinUpTime);

    //Constructing the Initial Min Up/Down Constraints
    init_up_down_const (fInitUpDownTime, fMinDownTime, fMinUpTime);

    //Constructing the Ramp Min Up/Down Constraints
    if(rampconst){
    ramp_const(fMaxRampDown,fMaxRampUp,fBoundDown, fBoundOn);}

    if(pc==1) {zt_Constraints();}
}

else if (model == 1){

Network_Matrix();

Convex_Constraints();
}


else if (model == 2){
  
  Network_Matrix();
  
  Sum_Ramp_Constraints();
  Sum_Pmin_Constraints();
  Sum_Bound_Constraints();
}

else if (model == 3){

Network_Matrix();

Semi_Ramp_Constraints();
Semi_Pmin_Constraints();
Semi_Bound_Constraints();

}


}

void AcadThermalUnitBlock::zt_Constraints( void ){

//storing initial p/c constraint of the form:
// (2*a*\bar_{p}_k + b) *p_t + (c - a*(\bar_{p}_k)^2) *u_t - z_t <= 0

UCBlock * f_UC_Block = dynamic_cast <UCBlock *> (get_f_Block());

PC_Const.resize(p_k.size()*f_UC_Block->get_t());
Zt_Const.resize(p_k.size()*f_UC_Block->get_t());
double value;
auto zt = Zt_Const.begin();
//auto pt = PC_Const.begin();
for (int i = 0; i < f_UC_Block->get_t() ; i++){
	for(int j = 0 ; j < p_k.size() ; j++){

	LinearConstraint::v_coeff_pair  *v_zt_pair = new LinearConstraint::v_coeff_pair(3);

	auto it = v_zt_pair->begin();

	//linear + quad
	value = fLinearTherm - lambda[i];
	*it = LinearConstraint::coeff_pair( &P[i] , 2 * fQuadTherm); //store p_i,t*b

	value = 2*fQuadTherm*p_k[j] + fLinearTherm - lambda[i] ; 
	*it = LinearConstraint::coeff_pair( &P[i] , value); //store (2*a*\bar_{p}_k + b) *p_t
	it++;



	value = fConstTherm - mew[i] * fMaxPower;
	value = fConstTherm - mew[i] * fMaxPower - fQuadTherm*p_k[j]*p_k[j] ; 
	*it = LinearConstraint::coeff_pair( &U[i] , value); //store (c - a*\bar_{p}^2) *u_t
	it++;

	*it = LinearConstraint::coeff_pair( &z_t[i] , -1.0);  //store - z_t
	 
	zt->add_variables( v_zt_pair, false);
	zt->set_Block(this);
	zt->set_lhs(-Inf<double>());
	zt->set_rhs(0);
	zt++;

	} //for j

}

add_static_constraint(Zt_Const);

}

void AcadThermalUnitBlock::generate_dynamic_constraints( void ){

//after having retrieved solution of the current CPLEX node we need to make the corresponding test
UCBlock * f_UC_Block = dynamic_cast <UCBlock *> (get_f_Block());
double value;
int count=0;
double tol = 1e-5;

 for(int i = 0 ; i<f_UC_Block->get_t(); i++){
		 cout<< " [ " << U[i].get_value() << " && " << P[i].get_value() << " |||| " << z_t[i].get_value() <<  " ] "; 	 
	 }
for (int i = 0; i < f_UC_Block->get_t() ; i++){
	if (U[i].get_value() != 0 ){
		value = (fQuadTherm * P[i].get_value() * P[i].get_value())/U[i].get_value();
		if (z_t[i].get_value() < value + tol) count++;
	}
}

std::list<LinearConstraint> * new_cuts = new std::list<LinearConstraint>;
new_cuts->resize(count);

std::shared_ptr<BlockModificationAD> mod = std::make_shared<BlockModificationAD>(BlockModificationAD::eAddConst );

auto nc = new_cuts->begin();
int count2=0;
double value2;
for (int i = 0; i < f_UC_Block->get_t() ; i++){
	if (U[i].get_value() != 0 ){
		value = (fQuadTherm * P[i].get_value() * P[i].get_value())/U[i].get_value();

		if (z_t[i].get_value() < value + tol){ // we need to add the new constraints
			
			count2++;
			LinearConstraint::v_coeff_pair  *v_cut_pair = new LinearConstraint::v_coeff_pair(3);
			auto it = v_cut_pair->begin();

			value2 = 2 * fQuadTherm * ( P[i].get_value() / U[i].get_value() )  +  fLinearTherm -  lambda[i] ;
			*it = LinearConstraint::coeff_pair( &P[i] , value2); 
			it++;
		
			value2 = fConstTherm - mew[i] * fMaxPower - fQuadTherm * ( P[i].get_value() / U[i].get_value() );
			*it = LinearConstraint::coeff_pair( &U[i] , value2);
			it++;

			*it = LinearConstraint::coeff_pair( &z_t[i] , -1.0);  //store - z_t

			nc->add_variables( v_cut_pair, false);
			nc->set_lhs(-Inf<double>());
			nc->set_rhs(0);
			nc->set_Block(this);
			nc++;

		} // if new constraint

	} // if u_t != 0
}//for loop
if(count > 0) {
        mod->whc_list = new_cuts;
	add_modification( mod );
}

}




void AcadThermalUnitBlock::Semi_Pmin_Constraints( void ){

Semi_Pmin_Const.resize(P_h.size());
//Construction of Constraint with form: P_min* Sum_{t = k,...,T}Y^+_{h,t} - p^h_k <=0
int j = 0;
int sum = 0;
int m = 0;
int ind = 0;
for(int i  = 0 ; i < v_ph_ind.size() ; i++){
	m = v_ph_ind[i]; //get index of last Y^+
	ind = Y_plus_pair[m].first;	//get the starting h index
	if (Y_plus_pair[m].first == 0) ind = 1;
	for (int k = 0 ; k < v_ph_size[i] ; k++){ //for all the p^h_k \ k = h,...,h+v_ph_size[i]
		sum = 0;
		for(int h = 0 ; h < v_ph_num[i] ; h++){		
		    if ( ind + k <= Y_plus_pair[m-h].second ) sum++;
		}//for h-loop
		LinearConstraint::v_coeff_pair  *v_p_min_pair = new LinearConstraint::v_coeff_pair(sum+1);
		auto it = v_p_min_pair->begin();
		for(int h = 0 ; h < v_ph_num[i] ; h++){
			if ( ind + k <= Y_plus_pair[m-h].second ){
				*it = LinearConstraint::coeff_pair( &Y_plus[m-h] , fMinPower); //store Y^+_{h,t}*P_min		
				 it++;
			}
		}
		*it = LinearConstraint::coeff_pair( &P_h[j] , -1); //store  - p^h_k
		Semi_Pmin_Const[j].add_variables( v_p_min_pair, false);
		Semi_Pmin_Const[j].set_Block(this);
		Semi_Pmin_Const[j].set_lhs(-Inf<double>());
		Semi_Pmin_Const[j].set_rhs(0);
	j++;		
	}//for k
}//for i

add_static_constraint(Semi_Pmin_Const);

}

void AcadThermalUnitBlock::Semi_Bound_Constraints( void ){

Semi_Pmax_Const.resize(P_h.size());

// Construction of Constraint exceeds two forms:
// for k = h: p^h_h - Bound_on * Sum_{t = h,...,T}   Y^+_{h,t} <=0
// for k > h: p^h_k - P_max    * Sum_{t = k+1,...,T} Y^+_{h,t} - Bound_down * Y^+_{h,k} <=0

int j = 0;
int sum = 0;
int m = 0;
int ind = 0;
int p_ind = 0;
int s =0;
for(int i  = 0 ; i < v_ph_ind.size() ; i++){
	m = v_ph_ind[i]; //get index of last Y^+
	ind = Y_plus_pair[m].first;	//get the starting h index
	s=1;
	sum = 0;
	if (Y_plus_pair[m].first == 0){
		ind = 1; //there is no bound on since the unit was initially turned on thus we jump to rest of constraints
		s=0; //this means that all the elements will have to go for k >h type of constraints!
	}
	

	else { //in case we are not in initial start up
        // we store p^h_h - Bound_on * Sum_{t = h,...,T}   Y^+_{h,t} <=0

		p_ind = v_ph_tot[i] - v_ph_size[i]; //storing the index of first p_h in the ph vectors
		for(int h = 0 ; h < v_ph_num[i] ; h++){		
		    if ( ind <= Y_plus_pair[m-h].second ) sum++;
		}//for h-loop
		LinearConstraint::v_coeff_pair  *v_p_max_pair = new LinearConstraint::v_coeff_pair(sum+1);
		auto it = v_p_max_pair->begin();
		for(int h = 0 ; h < v_ph_num[i] ; h++){
			if ( ind <= Y_plus_pair[m-h].second ){
				*it = LinearConstraint::coeff_pair( &Y_plus[m-h] , -fBoundOn); //store - Bound_on * Sum_{t = h,...,T}   Y^+_{h,t}	
				 it++;
			}
		}
		/*cout<<endl<<"Check_Bound_on: P_h.size() = " << P_h.size() << " && p_ind = " << p_ind << " && sum  = " << sum << " && ind = " << ind << " && Y_plus_pair[m].first = " << Y_plus_pair[m].first
		          << " && Y_plus_pair[m].second = " << Y_plus_pair[m].second << " && m = " << m;
		cout<<endl<<" // ";*/
		*it = LinearConstraint::coeff_pair( &P_h[p_ind] , 1); //store  p^h_h
		Semi_Pmax_Const[j].add_variables( v_p_max_pair, false);
		Semi_Pmax_Const[j].set_Block(this);
		Semi_Pmax_Const[j].set_lhs(-Inf<double>());
		Semi_Pmax_Const[j].set_rhs(0);
		j++;
	}//else of p^h_h
	
	// storing p^h_k - P_max * Sum_{t = k+1,...,T} Y^+_{h,t} - Bound_down * Y^+_{h,k} <=0
	for (int k = s ; k < v_ph_size[i] ; k++){ //for all the p^h_k with k = h+1,...,h+v_ph_size[i]
		sum = 0;
		for(int h = 0 ; h < v_ph_num[i] ; h++){		
		    if ( ind + k <= Y_plus_pair[m-h].second ) sum++;
		}//for h-loop
		LinearConstraint::v_coeff_pair  *v_p_min_pair = new LinearConstraint::v_coeff_pair(sum+1);
		auto it = v_p_min_pair->begin();
		for(int h = 0 ; h < v_ph_num[i] ; h++){
			if ( ind + k <= Y_plus_pair[m-h].second ){
				if(ind + k < Y_plus_pair[m-h].second)
					*it = LinearConstraint::coeff_pair( &Y_plus[m-h] , -fMaxPower); //store - Y^+_{h,t} * P_max
				else if(ind + k == Y_plus_pair[m-h].second)
					*it = LinearConstraint::coeff_pair( &Y_plus[m-h] , -fMaxPower); //store - Y^+_{h,k} * Bound_down				
				 it++;
			}
		}
		*it = LinearConstraint::coeff_pair( &P_h[j] , 1); //store  p^h_k
		Semi_Pmax_Const[j].add_variables( v_p_min_pair, false);
		Semi_Pmax_Const[j].set_Block(this);
		Semi_Pmax_Const[j].set_lhs(-Inf<double>());
		Semi_Pmax_Const[j].set_rhs(0);
	j++;		
	}//for k
}//for i


add_static_constraint(Semi_Pmax_Const);


}

void AcadThermalUnitBlock::Semi_Ramp_Constraints( void ){

int count = 0;
int j = 0;
for(int i = 0 ; i < v_ph_ind.size() ; ++i){
j=0;
if (Y_plus_pair[v_ph_ind[i]].first == 0) j = 1;
if(v_ph_size[i] > 1) count +=  v_ph_size[i] - 1 + j;

}

Semi_Ramp_up_Const.resize(count);
Semi_Ramp_down_Const.resize(count);

int sum = 0;
int m = 0;
int ind = 0;
int p_ind = 0;
j = 0;
for(int i = 0 ; i < v_ph_ind.size() ; ++i){

	m = v_ph_ind[i]; //get index of last Y^+
	ind = Y_plus_pair[m].first;	//get the starting h index
	if (Y_plus_pair[m].first == 0) ind++;
	if(v_ph_size[i] > 1){
		if (Y_plus_pair[m].first == 0){
			sum = 0;
			for(int h = 0 ; h < v_ph_num[i] ; h++){		
		    		if ( ind <= Y_plus_pair[m-h].second ) sum++;
			}//for h-loop
	
			//Init Ramp-Up: - sum_{k} y^+_{0,k} * (D^+ + P_0) + p^1_{k} <= 0  
			LinearConstraint::v_coeff_pair  *v_init_ramp_up_pair = new LinearConstraint::v_coeff_pair(sum+1);
			auto it_up = v_init_ramp_up_pair->begin();

			//Init Ramp-Down: -sum_{k} y^+_{0,k} * (D^- - P_0) - p^1_{k} <= 0
			LinearConstraint::v_coeff_pair  *v_init_ramp_dn_pair = new LinearConstraint::v_coeff_pair(sum+1);
			auto it_dn = v_init_ramp_dn_pair->begin();

			for(int h = 0 ; h < v_ph_num[i] ; h++){
				if ( ind <= Y_plus_pair[m-h].second ){
					*it_up = LinearConstraint::coeff_pair( &Y_plus[m-h] , - (fMaxRampUp+fInitPower) ); //store - sum_{k} y^+_{0,k} * (D^+ + P_0)		
					 it_up++;
					*it_dn = LinearConstraint::coeff_pair( &Y_plus[m-h] , - (fMaxRampDown-fInitPower) ); //store -sum_{k} y^+_{0,k} * (D^- - P_0)		
					 it_dn++;
				}
			}

			*it_up = LinearConstraint::coeff_pair( &P_h[0] ,   1 ); //store p^1_{k} 		
			*it_dn = LinearConstraint::coeff_pair( &P_h[0] , - 1 ); //store - p^1_{k}	

			/*cout<<endl<<"Init_Ramp_Check: j = " << j << "&& v_init_ramp_up_pair->size() =  " << v_init_ramp_up_pair->size()<< " && v_init_ramp_dn_pair->size() =  " 
                                  << v_init_ramp_dn_pair->size() << " && m = " << m << " && P_h.size() = " << P_h.size() << " && sum = " << sum << " && ind = " << ind;
			cout<<endl<<"//";*/

			Semi_Ramp_up_Const[j].add_variables( v_init_ramp_up_pair, false);
			Semi_Ramp_up_Const[j].set_Block(this);
			Semi_Ramp_up_Const[j].set_lhs(-Inf<double>());
			Semi_Ramp_up_Const[j].set_rhs(0);

			Semi_Ramp_down_Const[j].add_variables( v_init_ramp_dn_pair, false);
			Semi_Ramp_down_Const[j].set_Block(this);
			Semi_Ramp_down_Const[j].set_lhs(-Inf<double>());
			Semi_Ramp_down_Const[j].set_rhs(0);

			j++;
		}// if init

		p_ind = v_ph_tot[i] - v_ph_size[i]; //storing the index of first attached p_h
		int z=0;
		for (int k = 0 ; k < v_ph_size[i]-1 ; k++){ //for all the p^h_k with k = h,...,h+v_ph_size[i]-1
		//Construct: 
		//Ramp-Up:   p^h_t+1 - p^h_t - D^+ * sum_{t+1,...,T} y^+_{h,k} + p_min * y_{h,t} <= 0
		//Ramp-Down: p^h_t - p^h_t+1 - D^- * sum_{t+1,...,T} y^+_{h,k} - Bound_down * y_{h,t} <= 0
		z=p_ind+k;
			sum = 0;
			for(int h = 0 ; h < v_ph_num[i] ; h++){		
		    		if ( ind + k <= Y_plus_pair[m-h].second ) sum++;
			}//for h-loop
		
			LinearConstraint::v_coeff_pair  *v_ramp_up_pair = new LinearConstraint::v_coeff_pair(sum+2);
			auto it_up = v_ramp_up_pair->begin();

			LinearConstraint::v_coeff_pair  *v_ramp_dn_pair = new LinearConstraint::v_coeff_pair(sum+2);
			auto it_dn = v_ramp_dn_pair->begin();
			
			for(int h = 0 ; h < v_ph_num[i] ; h++){
				if ( ind + k == Y_plus_pair[m-h].second ){
					//cout<<endl<<"adding y_{h,t}";
					*it_up = LinearConstraint::coeff_pair( &Y_plus[m-h] , fMinPower); //store p_min * y_{h,t}		
				 	it_up++;
					*it_dn = LinearConstraint::coeff_pair( &Y_plus[m-h] , -fBoundDown); //store -Bound_down * y_{h,t}		
				 	it_dn++;
				}
				else if ( ind + k < Y_plus_pair[m-h].second ){
					//cout<<endl<<"adding sum_{t+1,...,T} y^+_{h,k}";
					*it_up = LinearConstraint::coeff_pair( &Y_plus[m-h] , -fMaxRampUp); //store - D^+ * sum_{t+1,...,T} y^+_{h,k}		
				 	it_up++;
					*it_dn = LinearConstraint::coeff_pair( &Y_plus[m-h] , -fMaxRampDown); //store - D^- * sum_{t+1,...,T} y^+_{h,k}		
				 	it_dn++;
				}
			}

			*it_up = LinearConstraint::coeff_pair( &P_h[z+1] , 1); //store p^h_t+1		
			it_up++;
			*it_up = LinearConstraint::coeff_pair( &P_h[z] , -1); //store - p^h_t		
			

			*it_dn = LinearConstraint::coeff_pair( &P_h[z] , 1); //store p^h_t		
			it_dn++;
			*it_dn = LinearConstraint::coeff_pair( &P_h[z+1] , -1); //store - p^h_t+1		

			/*cout<<endl<<"Ramp_Check: j = " << j << " && Semi_Ramp_up_Const.size() = " << Semi_Ramp_up_Const.size() <<
                                    " && v_ramp_up_pair->size() =  " << v_ramp_up_pair->size()<< " && Semi_Ramp_down_Const.size() = " << Semi_Ramp_down_Const.size() <<
                                    " && v_ramp_dn_pair->size() =  " << v_ramp_dn_pair->size() << " && z = " << z << " && P_h.size() = " << P_h.size() << " && sum = " << sum;
			cout<<endl<<"//";*/
			Semi_Ramp_up_Const[j].add_variables( v_ramp_up_pair, false);
			Semi_Ramp_up_Const[j].set_Block(this);
			Semi_Ramp_up_Const[j].set_lhs(-Inf<double>());
			Semi_Ramp_up_Const[j].set_rhs(0);

			Semi_Ramp_down_Const[j].add_variables( v_ramp_dn_pair, false);
			Semi_Ramp_down_Const[j].set_Block(this);
			Semi_Ramp_down_Const[j].set_lhs(-Inf<double>());
			Semi_Ramp_down_Const[j].set_rhs(0);

			j++;
		} //for k-loop




} //for i-loop

}
			/*cout<<endl<<"Finally j = " << j;
			cout<<endl<<"//";*/
add_static_constraint(Semi_Ramp_up_Const);
add_static_constraint(Semi_Ramp_down_Const);
}


void AcadThermalUnitBlock::Convex_Constraints( void ){

UCBlock * f_UC_Block = dynamic_cast <UCBlock *> (get_f_Block());

//Construct the Bound_On Constraints for all different nodes: p^h_{h, k} - y^+_{h,k}*Bound_on <= 0
//First we need to check the size, since Bound_on constraints apply only for h->k routes that start
//within the optimisation horizon
int count=0;
for(int i = 0; i < Y_plus.size() ; i++){
if(Y_plus_pair[i].first != 0) count++;
}

Bound_on_Const.resize(count);
int j = 0;
for(int i = Y_plus.size() - count ; i < Y_plus.size() ; ++i){

    /********* Bound_On Constraints *********/
    LinearConstraint::v_coeff_pair  *v_bound_on_pair = new LinearConstraint::v_coeff_pair(2);
    LinearConstraint::v_coeff_pair::iterator it;
    //Store the ramp down elements
	it = v_bound_on_pair->begin();
    /*if (i == 0 ) *it = LinearConstraint::coeff_pair( &P_hk[0] , 1); //store p^h_{h, k}
    else*/
	*it = LinearConstraint::coeff_pair( &P_hk[v_p_size[i]-v_p_ind[i]] , 1); //store p^h_{h, k}
	it = v_bound_on_pair->begin()+1;
	*it = LinearConstraint::coeff_pair( &Y_plus[i] , -fBoundOn); //store - y^+_{h,k}*Bound_on


	Bound_on_Const[j].add_variables( v_bound_on_pair, false);
	Bound_on_Const[j].set_Block(this);
	Bound_on_Const[j].set_lhs(-Inf<double>());
	Bound_on_Const[j].set_rhs(0);
	j++;
}

add_static_constraint(Bound_on_Const);

//Construct the Bound_Down Constraints for all different nodes: p^k_{h, k} - y^+_{h,k}*Bound_down <= 0
//First we need to check the size, since Bound_Down constraints apply only for h->k routes that finish
//within the optimisation horizon
count=0;
for(int i = 0; i < Y_plus.size() ; i++){
if( Y_plus_pair[i].second != f_UC_Block->get_t()+1 ) count++;
}
Bound_Down_Const.resize(count);
j = 0;
for(int i = 0; i < Y_plus.size() ; i++){

  
    if( Y_plus_pair[i].second != f_UC_Block->get_t()+1){

    /********* Bound_Down Constraints *********/
        LinearConstraint::v_coeff_pair  *v_bound_down_pair = new LinearConstraint::v_coeff_pair(2);
        LinearConstraint::v_coeff_pair::iterator it2;
        //Store the ramp down elements
	    it2 = v_bound_down_pair->begin();
            *it2 = LinearConstraint::coeff_pair( &P_hk[v_p_size[i]-1] , 1); //store p^k_{h, k}
	    it2 = v_bound_down_pair->begin()+1;
	    *it2 = LinearConstraint::coeff_pair( &Y_plus[i] , -fBoundDown); //store - y^+_{h,k}*Bound_Down


	    Bound_Down_Const[j].add_variables( v_bound_down_pair, false);
	    Bound_Down_Const[j].set_Block(this);
	    Bound_Down_Const[j].set_lhs(-Inf<double>());
	    Bound_Down_Const[j].set_rhs(0);
		j++;  
	    }
}
//Adding them to the Static Constraints Vector
add_static_constraint(Bound_Down_Const);

//Construct the Min Power Constraints for all different nodes: p_min * y^+_{h,k} - p^t_{h, k} <= 0 \for t = { h, ... , k}
Pmin_Const.resize(P_hk.size());
j = 0;
for(int i = 0 ; i < Pmin_Const.size() ; ++i){

    LinearConstraint::v_coeff_pair  *v_min_pair = new LinearConstraint::v_coeff_pair(2);
    LinearConstraint::v_coeff_pair::iterator it;
    //Store the ramp down elements
	it = v_min_pair->begin();
    if (i == v_p_size[j]) j++;
    *it = LinearConstraint::coeff_pair( &Y_plus[j], fMinPower ); //store p_min * y^+_{h,k}
	it = v_min_pair->begin()+1;
    *it = LinearConstraint::coeff_pair( &P_hk[i] , -1); //store - p^t_{h, k}
	

	Pmin_Const[i].add_variables( v_min_pair, false);
	Pmin_Const[i].set_Block(this);
	Pmin_Const[i].set_lhs(-Inf<double>());
	Pmin_Const[i].set_rhs(0);

}

//Adding them to the Static Constraints Vector
add_static_constraint(Pmin_Const);

// Construct the Max Power Constraints for all different nodes: p^t_{h, k} - p_max * y^+_{h,k} <= 0 \for t = { h+1 , ... , k-1 }
// The pmax constraint applies only for all (h,k) where the route is atleast three time steps and on top of that for special oc-
// casions of where the route is two time steps in cases:
// i) where the route continues from the previous time-horizon and the unit is active for the two first timesteps, in this case
//    pmax applies for the first time-step
// ii)where the route proceeds to the next time-horizon and the unit is active for the last two time-steps continuing to the next
//    one, where in this case p_max applies for the last time-step
count = 0;
for (int i = 0 ; i < Y_plus.size() ; i++){
//cout<<endl<<"( " << Y_plus_pair[i].first << " ; " << Y_plus_pair[i].second << " ) && v_p_ind[i] = "<< v_p_ind[i] << " && v_p_size[i] = " << v_p_size[i] <<  " && count = " << count;
if(v_p_ind[i] > 2) {
    if      (Y_plus_pair[i].second == f_UC_Block->get_t()+1 && Y_plus_pair[i].first == 0) {count = count +  v_p_ind[i]; /*cout<<" i_1 = " << i;*/}
    else if (Y_plus_pair[i].second == f_UC_Block->get_t()+1 || Y_plus_pair[i].first == 0) {count = count +  v_p_ind[i] - 1; /*cout<<" i_2 = " << i;*/}
    else if (Y_plus_pair[i].second != f_UC_Block->get_t()+1 && Y_plus_pair[i].first != 0) {count = count +  v_p_ind[i] - 2; /*cout<<" i_3 = " << i;*/}
}
else if (v_p_ind[i] == 2 && (Y_plus_pair[i].first == 0 || Y_plus_pair[i].second == f_UC_Block->get_t()+1 )) count++;
//cout<<" && upgraded_count = " << count;
}
//cout<<endl<<"yolo && count = " << count;
Pmax_Const.resize(count);

j = 0;
for (int i = 0 ; i < Y_plus.size() ; ++i){

if(v_p_ind[i] > 2){
int z = 0;
int l = 0;
if      (Y_plus_pair[i].second == f_UC_Block->get_t()+1 && Y_plus_pair[i].first == 0) {z=0; l=0;}
else if (Y_plus_pair[i].second != f_UC_Block->get_t()+1 && Y_plus_pair[i].first == 0) {z=1; l=0;}
else if (Y_plus_pair[i].second == f_UC_Block->get_t()+1 && Y_plus_pair[i].first != 0) {z=1; l=1;}
else if (Y_plus_pair[i].second != f_UC_Block->get_t()+1 && Y_plus_pair[i].first != 0) {z=2; l=1;}
    for (int k = 0 ; k < v_p_ind[i]-z ; k++){

        LinearConstraint::v_coeff_pair  *v_max_pair = new LinearConstraint::v_coeff_pair(2);
        LinearConstraint::v_coeff_pair::iterator it;
        //Store the max power elements
	    it = v_max_pair->begin();
	int s = v_p_size[i]-v_p_ind[i]+l+k;
        *it = LinearConstraint::coeff_pair( &P_hk[s], 1); 
                                        //store p^t_{h, k} for all h<t<k and we have v_p_size[i]-v_p_ind[i] giving us position h, 
                                        //so we need to add plus 1 to arrive to the first element that pma applies and the we continue with the counter
    	it = v_max_pair->begin()+1;
        *it = LinearConstraint::coeff_pair( &Y_plus[i], -fMaxPower ); //store - p_max * y^+_{h,k}
/*	cout<<endl<<"checking pair ( " << Y_plus_pair[i].first << " ; " << Y_plus_pair[i].second << " ): s = " << s << " v_p_size[i] = " << v_p_size[i] << " && v_p_ind[i] = " << v_p_ind[i] << " && z = " << z  << " && l = "  << l << " && k = " << k << " && P_hk.size() = " << P_hk.size() << " && Pmax_Const.size() = " << Pmax_Const.size() << " && j = " << j;
cout<<endl<<" //" ;*/
        Pmax_Const[j].add_variables( v_max_pair, false);
    	Pmax_Const[j].set_Block(this);
    	Pmax_Const[j].set_lhs(-Inf<double>());
    	Pmax_Const[j].set_rhs(0);

        j++;


    }

}

else if (v_p_ind[i] == 2 && Y_plus_pair[i].first == 0) {

        LinearConstraint::v_coeff_pair  *v_max_pair = new LinearConstraint::v_coeff_pair(2);
        LinearConstraint::v_coeff_pair::iterator it;
        //Store the max power elements
	    it = v_max_pair->begin();
        *it = LinearConstraint::coeff_pair( &P_hk[v_p_size[i]-2] , 1); //store p^t_{h, k} we need the first element of the two
    	it = v_max_pair->begin()+1;
        *it = LinearConstraint::coeff_pair( &Y_plus[i], -fMaxPower ); //store - p_max * y^+_{h,k}
       // cout<<endl<<"yolo i = " << i << " && also index = " << v_p_size[i]-2 << " && j = " << j;
       // cout<<endl<<"yolo i = " << i;
        Pmax_Const[j].add_variables( v_max_pair, false);
    	Pmax_Const[j].set_Block(this);
    	Pmax_Const[j].set_lhs(-Inf<double>());
    	Pmax_Const[j].set_rhs(0);

        j++;

}

else if (v_p_ind[i] == 2 && Y_plus_pair[i].second == f_UC_Block->get_t()+1) {

        LinearConstraint::v_coeff_pair  *v_max_pair = new LinearConstraint::v_coeff_pair(2);
        LinearConstraint::v_coeff_pair::iterator it;
        //Store the max power elements
	    it = v_max_pair->begin();
        *it = LinearConstraint::coeff_pair( &P_hk[v_p_size[i]-1] , 1); //store p^t_{h, k} we need the second element of the two
    	it = v_max_pair->begin()+1;
        *it = LinearConstraint::coeff_pair( &Y_plus[i], -fMaxPower ); //store - p_max * y^+_{h,k}

        Pmax_Const[j].add_variables( v_max_pair, false);
    	Pmax_Const[j].set_Block(this);
    	Pmax_Const[j].set_lhs(-Inf<double>());
    	Pmax_Const[j].set_rhs(0);

        j++;


}   

} // for i
//Adding them to the Static Constraints Vector
add_static_constraint(Pmax_Const);

count = 0;

for(int i = 0 ; i < Y_plus.size() ; ++i){

if (Y_plus_pair[i].first == 0 /*&& v_p_ind[i] > 1*/ ) count = count + 1;
if(v_p_ind[i] > 1) count = count +  v_p_ind[i] - 1;

}

RampDown_Const.resize(count);
RampUp_Const.resize(count);
//cout<<endl<<"yo count = " << count;
j=0;
int ind=0;
for(int i = 0 ; i < Y_plus.size() ; ++i){

//create initial ramp consts
if (Y_plus_pair[i].first == 0 /*&& v_p_ind[i] > 1*/ ){

//Init Ramp-Up: -y^+_{h,k} * (D^+ + P_0) + p^1_{h,k} <= 0  
        LinearConstraint::v_coeff_pair  *v_ramp_up_pair = new LinearConstraint::v_coeff_pair(2);
        LinearConstraint::v_coeff_pair::iterator it;
	
	it = v_ramp_up_pair->begin();
        *it = LinearConstraint::coeff_pair( &P_hk[v_p_size[i]-v_p_ind[i]] , 1);
    	it = v_ramp_up_pair->begin()+1;
        *it = LinearConstraint::coeff_pair( &Y_plus[i], - (fMaxRampUp+fInitPower) );  //fBoundDown-fMaxRampDown

        RampUp_Const[j].add_variables( v_ramp_up_pair, false);
    	RampUp_Const[j].set_Block(this);
    	RampUp_Const[j].set_lhs(-Inf<double>());
    	RampUp_Const[j].set_rhs(0); // fBoundDown - fInitPower

//Init Ramp-Down: -y^+_{h,k} * (D^- - P_0) - p^1_{h,k} <= 0
        LinearConstraint::v_coeff_pair  *v_ramp_down_pair = new LinearConstraint::v_coeff_pair(2);
        LinearConstraint::v_coeff_pair::iterator it2;
	it2 = v_ramp_down_pair->begin();
        *it2 = LinearConstraint::coeff_pair( &P_hk[v_p_size[i]-v_p_ind[i]] , -1);
	it2 = v_ramp_down_pair->begin()+1;
        *it2 = LinearConstraint::coeff_pair( &Y_plus[i], - (fMaxRampDown-fInitPower) );
			
        RampDown_Const[j].add_variables( v_ramp_down_pair, false);
    	RampDown_Const[j].set_Block(this);
    	RampDown_Const[j].set_lhs(-Inf<double>());
    	RampDown_Const[j].set_rhs(0);
        
        j++;
}

// Ramp Up/Down Constraints are only valid for routes bigger than one time step
if(v_p_ind[i] > 1){

	
   for (int k = 0 ; k < v_p_ind[i] - 1 ; k++){
           
        /********* Ramp_Up Constraints *********/
        //Construct the Ramp Up   Constraints for all different nodes: p^{t+1}_{h, k} - p^{t}_{h, k}     - D^+ * y^+_{h,k}  <= 0 \for t = { h , ... , k-1 }

        LinearConstraint::v_coeff_pair  *v_ramp_up_pair = new LinearConstraint::v_coeff_pair(3);
        LinearConstraint::v_coeff_pair::iterator it;
        //Store the Ramp_Up elements

	it  = v_ramp_up_pair->begin();
	ind = v_p_size[i]-v_p_ind[i]+k+1;
        *it = LinearConstraint::coeff_pair( &P_hk[ ind ], 1); 
                                        //store p^t+1_{h, k} for all h<=t<k and we have v_p_size[i]-v_p_ind[i] giving us position h, 
                                        // and then we continue with the counter, since we want the t+1 element we add plus one to the 
    	it  = v_ramp_up_pair->begin()+1;
	ind = v_p_size[i]-v_p_ind[i]+k;
        *it = LinearConstraint::coeff_pair( &P_hk[ ind ] , - 1);
                                         //store p^t_{h, k} for all h<=t<k and we have v_p_size[i]-v_p_ind[i] giving us position h, 
                                        //so we start from there as first element that applies ramp constraints and the we continue with the counter
    	it = v_ramp_up_pair->begin()+2;
        *it = LinearConstraint::coeff_pair( &Y_plus[i], -fMaxRampUp ); //store - y^+_{h,k} * D^+
	/*cout<<endl<<"Yolo: ind = " << ind << " && k = " << k << " && v_p_size["<<i<<"] = " << v_p_size[i] << " && v_p_ind["<<i<<"] = " << v_p_ind[i] << " && P_hk.size() = " << P_hk.size()
                  <<" && Y_plus.size() = " << Y_plus.size();
	cout<<endl<<"//";*/
        RampUp_Const[j].add_variables( v_ramp_up_pair, false);
    	RampUp_Const[j].set_Block(this);
    	RampUp_Const[j].set_lhs(-Inf<double>());
    	RampUp_Const[j].set_rhs(0);


        /********* Ramp_Down Constraints *********/
	//Construct the Ramp Down Constraints for all different nodes: p^t_{h, k}     - p^{t+1}_{h, k}   - D^- * y^+_{h,k}  <= 0 \for t = { h , ... , k-1 }

        LinearConstraint::v_coeff_pair  *v_ramp_down_pair = new LinearConstraint::v_coeff_pair(3);
        LinearConstraint::v_coeff_pair::iterator it2;
        //Store the Ramp_Down elements
	    it2 = v_ramp_down_pair->begin();
        *it2 = LinearConstraint::coeff_pair( &P_hk[v_p_size[i]-v_p_ind[i]+k] , 1); //store p^{t}_{h, k}
    	it2 = v_ramp_down_pair->begin()+1;
        *it2 = LinearConstraint::coeff_pair( &P_hk[v_p_size[i]-v_p_ind[i]+k+1] , - 1); //store -p^{t+1}_{h, k}
    	it2 = v_ramp_down_pair->begin()+2;
        *it2 = LinearConstraint::coeff_pair( &Y_plus[i], -fMaxRampDown ); //store - y^+_{h,k} * D^-

        RampDown_Const[j].add_variables( v_ramp_down_pair, false);
        RampDown_Const[j].set_Block(this);
        RampDown_Const[j].set_lhs(-Inf<double>());
        RampDown_Const[j].set_rhs(0);

        j++;

        }   

  }



}//for i


//Adding them to the Static Constraints Vector
add_static_constraint(RampUp_Const);
add_static_constraint(RampDown_Const);



}

void AcadThermalUnitBlock::Sum_Bound_Constraints (){


//Construct Constraints that for each p_t take into account bound on/down and p_max having following form:
// p_t - y^+_{h,t}*Bound_down - y^+_{t,k}*Bound_on - sum_{t \in (h,k)} y^+_{h,k}*p_max <= 0 within the opt horizon

Bound_Down_Const.resize(v_t.size());

for(int i = 0; i < v_t.size() ; i++){

        LinearConstraint::v_coeff_pair  *v_bound_down_pair = new LinearConstraint::v_coeff_pair(1+v_t[i].b_dn.size() + v_t[i].b_on.size() + v_t[i].p_mx.size() );

    	/*cout<<endl<<"Check for i = " << i << " && v_t[i].b_dn.size() " << v_t[i].b_dn.size() <<"&& v_t[i].b_on.size() = " << v_t[i].b_on.size() 
              <<" v_t[i].p_mx.size() = " << v_t[i].p_mx.size() << " && v_bound_down_pair.size() = " << v_bound_down_pair->size();*/

        LinearConstraint::v_coeff_pair::iterator it2;
        it2 = v_bound_down_pair->begin();
        *it2 = LinearConstraint::coeff_pair( &P[i] , 1); //store p_t
    
        if (v_t[i].b_dn.size() > 0){
	   for(auto j:v_t[i].b_dn){ 
	    	it2++;
	    	*it2 = LinearConstraint::coeff_pair( &Y_plus[j] , -fBoundDown); //store - y^+_{h,t}*Bound_Down
	    }   
	}

        if (v_t[i].b_on.size() > 0){
             for(auto j : v_t[i].b_on){
	            it2++;
	            *it2 = LinearConstraint::coeff_pair( &Y_plus[j] , -fBoundOn); //store - y^+_{t,k}*Bound_On
	     }
	}

        if (v_t[i].p_mx.size() > 0){
            for(auto j : v_t[i].p_mx){
	            it2++;
	            *it2 = LinearConstraint::coeff_pair( &Y_plus[j] , -fMaxPower); //store - y^+_{h,k}*P_max for h < t < k 
            }
	}

	Bound_Down_Const[i].add_variables( v_bound_down_pair, false);
	Bound_Down_Const[i].set_Block(this);
	Bound_Down_Const[i].set_lhs(-Inf<double>());
	Bound_Down_Const[i].set_rhs(0);
              
} //for loop
//Adding them to the Static Constraints Vector
add_static_constraint(Bound_Down_Const);
 
}

void AcadThermalUnitBlock::Sum_Pmin_Constraints (){

//Construct the Min Power Constraints for all different nodes: p_min * y^+_{h,k} - p^t_{h, k} <= 0 \for t = { h, ... , k}
Pmin_Const.resize(v_t.size());
int j = 0;
for(int i = 0 ; i < v_t.size() ; ++i){
    LinearConstraint::v_coeff_pair  *v_min_pair = new LinearConstraint::v_coeff_pair(1+v_t[i].index.size());
    LinearConstraint::v_coeff_pair::iterator it;
    //Store the ramp down elements
	it = v_min_pair->begin();
    *it = LinearConstraint::coeff_pair( &P[i] , -1); //store - p^t_{h, k}
    for(j = 0 ; j < v_t[i].index.size() ; j++){
	    it = v_min_pair->begin()+1+j;
        *it = LinearConstraint::coeff_pair( &Y_plus[v_t[i].index[j]], fMinPower ); //store p_min * y^+_{h,k}
    }

	

	Pmin_Const[i].add_variables( v_min_pair, false);
	Pmin_Const[i].set_Block(this);
	Pmin_Const[i].set_lhs(-Inf<double>());
	Pmin_Const[i].set_rhs(0);
}

//Adding them to the Static Constraints Vector
add_static_constraint(Pmin_Const);

/*int count = 0;
for(int i = 0 ; i < v_t.size() ; i++){
    if(v_t[i].p_mx.size() > 0) count++;
}
Pmax_Const.resize(count);
int l = 0;
for(int i = 0 ; i < v_t.size() ; ++i){
   if(v_t[i].p_mx.size() > 0) {
    LinearConstraint::v_coeff_pair  *v_max_pair = new LinearConstraint::v_coeff_pair(1+v_t[i].p_mx.size());
    LinearConstraint::v_coeff_pair::iterator it;
    //Store the ramp down elements
	it = v_max_pair->begin();
    *it = LinearConstraint::coeff_pair( &P[i] , 1); 
    for(j = 0 ; j < v_t[i].p_mx.size() ; j++){
	    it = v_max_pair->begin()+1+j;
        *it = LinearConstraint::coeff_pair( &Y_plus[v_t[i].p_mx[j]], -fMaxPower ); 
    }

	

	Pmax_Const[l].add_variables( v_max_pair, false);
	Pmax_Const[l].set_Block(this);
	Pmax_Const[l].set_lhs(-Inf<double>());
	Pmax_Const[l].set_rhs(0);
    l++;
}

}

//Adding them to the Static Constraints Vector
add_static_constraint(Pmax_Const);*/
}

void AcadThermalUnitBlock::Sum_Ramp_Constraints (){

int count = 0;
int z = 0;
for (int  i = 0 ; i < v_t.size() ; i++) 
	if (v_t[i].sindex.size() > 0 ) count++;

if (r_init.size() > 0){
	count++;
	z++;
}
	RampDown_Const.resize(count);
	RampUp_Const.resize(count);



    if (r_init.size() > 0){ // this means we have initial ramp up constraints

	//Init Ramp-Up: - \sum_{ 1 \in [h,k] } y^+_{h,k} * (D^+ + P_0) + p_1 <= 0  
	LinearConstraint::v_coeff_pair  *v_ramp_init_up_pair   = new LinearConstraint::v_coeff_pair(1+ r_init.size() );
	LinearConstraint::v_coeff_pair::iterator it;
	it = v_ramp_init_up_pair->begin();
	*it = LinearConstraint::coeff_pair( &P[0] , 1); //store p_1
        for(int i = 0 ; i < r_init.size() ; i++  ){
    		it = v_ramp_init_up_pair->begin()+1+i;
            	*it = LinearConstraint::coeff_pair ( &Y_plus[r_init[i]] , - (fMaxRampUp+fInitPower) );//store - \sum_{ 1 \in [h,k] } y^+_{h,k}
	}		

	RampUp_Const[0].add_variables( v_ramp_init_up_pair, false);
	RampUp_Const[0].set_Block(this);
	RampUp_Const[0].set_lhs(-Inf<double>());
	RampUp_Const[0].set_rhs(0);



	//Init Ramp-Down: - \sum_{ 1 \in [h,k] } y^+_{h,k} * (D^- - P_0) - p_1 <= 0

	LinearConstraint::v_coeff_pair  *v_ramp_init_down_pair = new LinearConstraint::v_coeff_pair(1+ r_init.size() );
	LinearConstraint::v_coeff_pair::iterator it2;
	it2 = v_ramp_init_down_pair->begin();
	*it2 = LinearConstraint::coeff_pair( &P[0] , -1); //store - p_1	 
        for(int i = 0 ; i < r_init.size() ; i++  ){
    		it2 = v_ramp_init_down_pair->begin()+1+i;
            	*it2 = LinearConstraint::coeff_pair ( &Y_plus[r_init[i]] , - (fMaxRampDown-fInitPower) ); //store - \sum_{ 1 \in [h,k] } y^+_{h,k} * (D^- - P_0)
            }
  	
 	RampDown_Const[0].add_variables( v_ramp_init_down_pair, false); 
	RampDown_Const[0].set_Block(this);
        RampDown_Const[0].set_lhs(-Inf<double>());
        RampDown_Const[0].set_rhs(0);
}


       

			



//int j = 0;
int r = 0;
//int m = 0;

for(int i = z ; i < RampUp_Const.size() ; ++i){

while (v_t[r].sindex.size() == 0 ) r++;

        /********* Ramp_Up Constraints *********/
	//Construct the Ramp Up   Constraints for all different nodes: p^{t+1} - p^{t}     - D^+ * \sum_{ t \in [h,k-1] } y^+_{h,k}  <= 0 \for t = { 1 , ... , T-1 }

        LinearConstraint::v_coeff_pair  *v_ramp_up_pair = new LinearConstraint::v_coeff_pair(2 + v_t[r].sindex.size());
        auto it = v_ramp_up_pair->begin();
        *it = LinearConstraint::coeff_pair( &P[r+1], 1); 
    	it++;// = v_ramp_up_pair->begin()+1;
        *it = LinearConstraint::coeff_pair( &P[r] , - 1);

	for(auto j : v_t[r].sindex ){
    	it++;
        *it = LinearConstraint::coeff_pair( &Y_plus[j], -fMaxRampUp + fBoundOn ); //store - D^+ * \sum_{ t \in [h,k-1] } y^+_{h,k}
        }

        RampUp_Const[i].add_variables( v_ramp_up_pair, false);
    	RampUp_Const[i].set_Block(this);
    	RampUp_Const[i].set_lhs(-Inf<double>());
    	RampUp_Const[i].set_rhs(fBoundOn);


        /********* Ramp_Down Constraints *********/
	//Construct the Ramp Down Constraints for all different nodes: p^t     - p^{t+1}   - D^- * \sum_{ t \in [h,k-1] } y^+_{h,k}  <= 0 \for t = { 1 , ... , T-1 }

        LinearConstraint::v_coeff_pair  *v_ramp_down_pair = new LinearConstraint::v_coeff_pair(2 + v_t[r].sindex.size());
        //Store the Ramp_Down elements
	auto it2 = v_ramp_down_pair->begin();
        *it2 = LinearConstraint::coeff_pair( &P[r] , 1); //store p^{t}_{h, k}
    	it2++;// = v_ramp_down_pair->begin()+1;
        *it2 = LinearConstraint::coeff_pair( &P[r+1] , - 1); //store -p^{t+1}_{h, k}

	for( auto j : v_t[r].sindex ) {
    	it2++;
        *it2 = LinearConstraint::coeff_pair( &Y_plus[j], -fMaxRampDown + fBoundDown ); //store - y^+_{h,k} * D^-
        }

        RampDown_Const[i].add_variables( v_ramp_down_pair, false);
        RampDown_Const[i].set_Block(this);
        RampDown_Const[i].set_lhs(-Inf<double>());
        RampDown_Const[i].set_rhs(fBoundDown);

r++;
        }   

//Adding them to the Static Constraints Vector
add_static_constraint(RampUp_Const);
add_static_constraint(RampDown_Const);

}

void AcadThermalUnitBlock::Network_Matrix(){

//UCBlock * f_UC_Block = dynamic_cast <UCBlock *> (get_f_Block());

Net_Matrix.resize(v_minus.size() + v_plus.size() + 2); 

/*we scan all elements of Y_plus vector and for each one of the y^+_h,k variables we add a positive entry to the
  h-row of the coeff matrix frmo where the arc is leaving and then a negative entry to the k-row of the coeff mat-
  rix where the arc is arriving

*/
for(int i = 0; i < Y_plus.size() ; ++i ){

        LinearConstraint::v_coeff_pair * v_t_pair = new LinearConstraint::v_coeff_pair(1);
        LinearConstraint::v_coeff_pair::iterator it1;
	    it1 = v_t_pair->begin();
	    *it1 = LinearConstraint::coeff_pair( &Y_plus[i] , 1.0); 
        int k = mat_index(Y_plus_pair[i].first,0);


        Net_Matrix[k].add_variables( v_t_pair, false); 

        LinearConstraint::v_coeff_pair * v_s_pair = new LinearConstraint::v_coeff_pair(1);
        LinearConstraint::v_coeff_pair::iterator it2;
	    it2 = v_s_pair->begin();
	    *it2 = LinearConstraint::coeff_pair( &Y_plus[i] , -1.0); 
        int h = mat_index(Y_plus_pair[i].second,1);

        Net_Matrix[h].add_variables( v_s_pair, false); 


}

for(int i = 0; i < Y_minus.size() ; ++i ){

        LinearConstraint::v_coeff_pair * v_t_pair = new LinearConstraint::v_coeff_pair(1);
        LinearConstraint::v_coeff_pair::iterator it1;
	    it1 = v_t_pair->begin();
	    *it1 = LinearConstraint::coeff_pair( &Y_minus[i] , 1.0); 
        int k = mat_index(Y_minus_pair[i].first,1);
        Net_Matrix[k].add_variables( v_t_pair, false); 

        LinearConstraint::v_coeff_pair * v_s_pair = new LinearConstraint::v_coeff_pair(1);
        LinearConstraint::v_coeff_pair::iterator it2;
	    it2 = v_s_pair->begin();
	    *it2 = LinearConstraint::coeff_pair( &Y_minus[i] , -1.0); 
        int h = mat_index(Y_minus_pair[i].second,0);
        Net_Matrix[h].add_variables( v_s_pair, false); 


}

	Net_Matrix[0].set_Block(this);
    Net_Matrix[0].set_lhs(1.0);
    Net_Matrix[0].set_rhs(1.0);

//proceed with filling all the rest info of the constraints
for(int i = 1 ; i < Net_Matrix.size() -1 ; i++){

	Net_Matrix[i].set_Block(this);
    Net_Matrix[i].set_lhs(0.0);
    Net_Matrix[i].set_rhs(0.0);

    }

Net_Matrix[Net_Matrix.size() -1].set_Block(this);
Net_Matrix[Net_Matrix.size() -1].set_lhs(-1.0);
Net_Matrix[Net_Matrix.size() -1].set_rhs(-1.0);

add_static_constraint(Net_Matrix);

}


void AcadThermalUnitBlock::p_min_max_const(double pmin, double pmax){
		 

UCBlock * f_UC_Block = dynamic_cast <UCBlock *> (get_f_Block());

	//Cosntrution of Pmin_Constraint of the form: u_it P_i^min - p_it <= 0
	//Cosntrution of Pmax_Constraint of the form: p_it - u_it P_i^max <= 0

	Pmin_Const.resize( f_UC_Block->get_t());
	Pmax_Const.resize( f_UC_Block->get_t() );

		for(int i=0; i<f_UC_Block->get_t(); i++){
		
		//initilaize the vector of pairs of the variables and coefficients of the corresponding constraints
	        LinearConstraint::v_coeff_pair * v_min_pair = new LinearConstraint::v_coeff_pair(2);
		LinearConstraint::v_coeff_pair * v_max_pair = new LinearConstraint::v_coeff_pair(2);

		LinearConstraint::v_coeff_pair::iterator it;
		it = v_min_pair->begin();
		*it = LinearConstraint::coeff_pair( &U[i] , pmin); //store u_it P_i^min u_it
		it = v_min_pair->begin() + 1;
		*it = LinearConstraint::coeff_pair( &P[i] , -1.0); //store - p_it

		it = v_max_pair->begin();
		*it = LinearConstraint::coeff_pair( &U[i] , -pmax); // store - u_it P_i^max
		it = v_max_pair->begin()+1;
		*it = LinearConstraint::coeff_pair( &P[i] , 1.0); // store - p_it of Pmax Const

		//Initializing the Pmin & Pmax Constraints
					 

		 Pmin_Const[i].add_variables( v_min_pair, false); 
		 Pmin_Const[i].set_Block(this);
                 Pmin_Const[i].set_lhs(-Inf<double>());
		 Pmin_Const[i].set_rhs(0);
		 
		 Pmax_Const[i].add_variables( v_max_pair, false); 
		 Pmax_Const[i].set_Block(this);
                 Pmax_Const[i].set_lhs(-Inf<double>());
		 Pmax_Const[i].set_rhs(0);
				 

		}

	//Adding them to the Static Constraints Vector
	add_static_constraint(Pmin_Const);
	add_static_constraint(Pmax_Const);

}


void AcadThermalUnitBlock::min_up_down_const (double mindown, double minup){

UCBlock * f_UC_Block = dynamic_cast <UCBlock *> (get_f_Block());

//initializing u-Z connection Constraints
int t1 = f_UC_Block->get_t() - init_t;
UZ_Const.resize(t1);

for(int i=0; i < t1 ; i++){

	if(i+1 >= minup){

       LinearConstraint::v_coeff_pair * v_uz_pair = new LinearConstraint::v_coeff_pair(2+minup);
	   auto it = v_uz_pair->begin();
	   *it = LinearConstraint::coeff_pair( &U[i+init_t] , 1);
	   it = v_uz_pair->begin() + 1;
   	   *it = LinearConstraint::coeff_pair( &U_11[i] , -1);
	   for (int j = 0 ; j < minup; j++){
		   it = v_uz_pair->begin() + 2 + j;
   		   *it = LinearConstraint::coeff_pair( &U_01[i - j] , -1);
	   }
		
		UZ_Const[i].add_variables(v_uz_pair, false);
 		UZ_Const[i].set_Block(this);
		UZ_Const[i].set_lhs(0);
		UZ_Const[i].set_rhs(0);
	}

	else {

      LinearConstraint::v_coeff_pair * v_uz_pair = new LinearConstraint::v_coeff_pair(3+i);
	   auto it = v_uz_pair->begin();
	   *it = LinearConstraint::coeff_pair( &U[i+init_t] , 1);
	   it = v_uz_pair->begin() + 1;
   	   *it = LinearConstraint::coeff_pair( &U_11[i] , -1);
	   it = v_uz_pair->begin() + 2;
   	   *it = LinearConstraint::coeff_pair( &U_01[i] , -1);
	   for (int j = 0 ; j < i; j++){
		   it = v_uz_pair->begin() + 3 + j;
   		   *it = LinearConstraint::coeff_pair( &U_01[i -1 - j] , -1);		
	   }
		
		UZ_Const[i].add_variables(v_uz_pair, false);
 		UZ_Const[i].set_Block(this);
		UZ_Const[i].set_lhs(0);
		UZ_Const[i].set_rhs(0);

	}

	}
//initializing j_on Constraints
J_On_Const.resize(t1);

LinearConstraint::v_coeff_pair * v_j_on_pair = new LinearConstraint::v_coeff_pair(2);
auto it = v_j_on_pair->begin();
*it = LinearConstraint::coeff_pair( &U_11[0] , 1);
it = v_j_on_pair->begin() + 1;
*it = LinearConstraint::coeff_pair( &U_10[0] , 1);

		J_On_Const[0].add_variables(v_j_on_pair, false);
 		J_On_Const[0].set_Block(this);
		if(fInitUpDownTime < 0){
			J_On_Const[0].set_lhs(0);
			J_On_Const[0].set_rhs(0);
		}
		if(fInitUpDownTime > 0){
			J_On_Const[0].set_lhs(1);
			J_On_Const[0].set_rhs(1);
		}

for(int i=0; i < t1-1 ; i++){

	if(i+1 >= minup){

       LinearConstraint::v_coeff_pair * v_j_on_pair = new LinearConstraint::v_coeff_pair(4);
	   auto it = v_j_on_pair->begin();
	   *it = LinearConstraint::coeff_pair( &U_11[i] , 1);
	   it = v_j_on_pair->begin() + 1;
   	   *it = LinearConstraint::coeff_pair( &U_11[i+1] , -1);
	   it = v_j_on_pair->begin() + 2;
   	   *it = LinearConstraint::coeff_pair( &U_10[i+1] , -1);
	   it = v_j_on_pair->begin() + 3;
   	   *it = LinearConstraint::coeff_pair( &U_01[i+1-minup] , 1);	


		J_On_Const[i+1].add_variables(v_j_on_pair, false);
 		J_On_Const[i+1].set_Block(this);
		J_On_Const[i+1].set_lhs(0);
		J_On_Const[i+1].set_rhs(0);			
	 
	}

	else {

       LinearConstraint::v_coeff_pair * v_j_on_pair = new LinearConstraint::v_coeff_pair(3);
	   auto it = v_j_on_pair->begin();
	   *it = LinearConstraint::coeff_pair( &U_11[i] , 1);
       it = v_j_on_pair->begin() + 1;
   	   *it = LinearConstraint::coeff_pair( &U_11[i+1] , -1);
	   it = v_j_on_pair->begin() + 2;
   	   *it = LinearConstraint::coeff_pair( &U_10[i+1] , -1);

		J_On_Const[i+1].add_variables(v_j_on_pair, false);
 		J_On_Const[i+1].set_Block(this);
		J_On_Const[i+1].set_lhs(0);
		J_On_Const[i+1].set_rhs(0);

	}

}


//initializing J_Off_Const Constraints
J_Off_Const.resize(t1);

LinearConstraint::v_coeff_pair * v_j_off_pair = new LinearConstraint::v_coeff_pair(2);
auto it1 = v_j_off_pair->begin();
*it1 = LinearConstraint::coeff_pair( &U_00[0] , 1);
it1 = v_j_off_pair->begin() + 1;
*it1 = LinearConstraint::coeff_pair( &U_01[0] , 1);

		J_Off_Const[0].add_variables(v_j_off_pair, false);
 		J_Off_Const[0].set_Block(this);
		if(fInitUpDownTime < 0){
			J_Off_Const[0].set_lhs(1);
			J_Off_Const[0].set_rhs(1);
		}
		if(fInitUpDownTime > 0){
			J_Off_Const[0].set_lhs(0);
			J_Off_Const[0].set_rhs(0);
		}

for(int i=0; i < t1-1 ; i++){

	if(i+1 >= mindown){

       LinearConstraint::v_coeff_pair * v_j_off_pair = new LinearConstraint::v_coeff_pair(4);
	   auto it = v_j_off_pair->begin();
	   *it = LinearConstraint::coeff_pair( &U_00[i] , 1);
	   it = v_j_off_pair->begin() + 1;
   	   *it = LinearConstraint::coeff_pair( &U_00[i+1] , -1);
	   it = v_j_off_pair->begin() + 2;
   	   *it = LinearConstraint::coeff_pair( &U_01[i+1] , -1);
	   it = v_j_off_pair->begin() + 3;
   	   *it = LinearConstraint::coeff_pair( &U_10[i+1-mindown] , 1);	


		J_Off_Const[i+1].add_variables(v_j_off_pair, false);
 		J_Off_Const[i+1].set_Block(this);
		J_Off_Const[i+1].set_lhs(0);
		J_Off_Const[i+1].set_rhs(0);			
	 
	}

	else {

       LinearConstraint::v_coeff_pair * v_j_off_pair = new LinearConstraint::v_coeff_pair(3);
	   auto it = v_j_off_pair->begin();
	   *it = LinearConstraint::coeff_pair( &U_00[i] , 1);
	   it = v_j_off_pair->begin() + 1;
   	   *it = LinearConstraint::coeff_pair( &U_00[i+1] , -1);
	   it = v_j_off_pair->begin() + 2;
   	   *it = LinearConstraint::coeff_pair( &U_01[i+1] , -1);

		J_Off_Const[i+1].add_variables(v_j_off_pair, false);
 		J_Off_Const[i+1].set_Block(this);
		J_Off_Const[i+1].set_lhs(0);
		J_Off_Const[i+1].set_rhs(0);

	}

}

add_static_constraint(UZ_Const);
add_static_constraint(J_Off_Const);
add_static_constraint(J_On_Const);


}

void AcadThermalUnitBlock::init_up_down_const (double t0, double mindown, double minup){
			
//comment ouf for dynamic part
Init_Const.resize(init_t);
//Init_Const.resize(1);

	// Based on the given data the unit may have to be set on or off on the time step 0

	//initilaize the vector of pairs of the variables and coefficients of the corresponding constraints
	LinearConstraint::v_coeff_pair  *v_init_pair = new LinearConstraint::v_coeff_pair(1);
	LinearConstraint::v_coeff_pair::iterator it = v_init_pair->begin();
	*it = LinearConstraint::coeff_pair( &U[0] , 1); //store u_{i,0}

	if(t0<0 && -t0 < mindown){ //we need to set u_{i,0} = 0
		
		for (int i = 0; i < init_t ; ++i){
		
			LinearConstraint::v_coeff_pair  *v_init_pair = new LinearConstraint::v_coeff_pair(1);
			LinearConstraint::v_coeff_pair::iterator it = v_init_pair->begin();
			*it = LinearConstraint::coeff_pair( &U[i] , 1); //store u_{i,0}
		
			Init_Const[i].add_variables(v_init_pair, false); 
			Init_Const[i].set_Block(this);
       		Init_Const[i].set_lhs(0);
			Init_Const[i].set_rhs(0);
		}
			add_static_constraint(Init_Const);
		   //Add vector to Static Constraint Vector
	}
	else if (t0>0 && t0 < minup){ //we need to set u_{i,0} = 1

		for (int i = 0; i < init_t ; ++i){
		
			LinearConstraint::v_coeff_pair  *v_init_pair = new LinearConstraint::v_coeff_pair(1);
			LinearConstraint::v_coeff_pair::iterator it = v_init_pair->begin();
			*it = LinearConstraint::coeff_pair( &U[i] , 1); //store u_{i,0}
		
			Init_Const[i].add_variables(v_init_pair, false);
			Init_Const[i].set_Block(this);
          	Init_Const[i].set_lhs(1);
			Init_Const[i].set_rhs(1);
		}
			add_static_constraint(Init_Const); 
		
                //Add vectors to Static Constraint Vector
	}
	delete v_init_pair;//Tiziano
}
//                         ramp_const(fMaxRampDown,  fMaxRampUp,  fBoundDown,     fBoundOn);}
void AcadThermalUnitBlock::ramp_const(double D_down, double D_up, double pmin, double pmax){


UCBlock * f_UC_Block = dynamic_cast <UCBlock *> (get_f_Block());

	/*
	 * Cosntrution of Ramp_Down_Constraint of the form:
	 *  p_{i, t} - p_{i,t+1} +u_{i,t+1}*(p_min-D^down) <= p_min
	 *
	 * Cosntrution of Ramp_Up_Constraint of the form:
	 * p_{i, t+1} - p_{i,t} + u_{i,t}*(p_min-D^up) <= p_min
	 */

				 

	RampDown_Const.resize(f_UC_Block->get_t());
	RampUp_Const.resize(f_UC_Block->get_t());



			//initilaize the vector of pairs of the variables and coefficients of the corresponding constraints
			LinearConstraint::v_coeff_pair  *v_ramp_init_down_pair = new LinearConstraint::v_coeff_pair(1);
			LinearConstraint::v_coeff_pair  *v_ramp_init_up_pair = new LinearConstraint::v_coeff_pair(2);
			LinearConstraint::v_coeff_pair::iterator it;
			//Store the ramp down elements
			it = v_ramp_init_down_pair->begin();
			*it = LinearConstraint::coeff_pair( &P[0] , 1); //store p_{i, t}

			//Store the ramp up elements
			it = v_ramp_init_up_pair->begin();
			*it = LinearConstraint::coeff_pair ( &U[0] , (pmin-D_down) ); //store p_{i, t+1}
			it = v_ramp_init_up_pair->begin()+1;
			*it = LinearConstraint::coeff_pair( &P[0] , -1); //store - p_{i,t}
			
			//Initializing the Ramp Up & Down Constraints
  		    double initial_rhs = (fInitUpDownTime > 0 );
			double rhs;
			rhs = ( initial_rhs ?  - fInitPower - D_up : -  pmax );	

		 	RampDown_Const[0].add_variables( v_ramp_init_down_pair, false); 
			RampDown_Const[0].set_Block(this);
            RampDown_Const[0].set_lhs(-Inf<double>());
		    RampDown_Const[0].set_rhs(-rhs);
				 
			rhs = pmin - fInitPower ;
					
			RampUp_Const[0].add_variables( v_ramp_init_up_pair, false);
			RampUp_Const[0].set_Block(this);
			RampUp_Const[0].set_lhs(-Inf<double>());
			RampUp_Const[0].set_rhs(rhs);

	for(int i=0; i<f_UC_Block->get_t()-1; i++){

			//initilaize the vector of pairs of the variables and coefficients of the corresponding constraints
			LinearConstraint::v_coeff_pair  *v_ramp_down_pair = new LinearConstraint::v_coeff_pair(3);
			LinearConstraint::v_coeff_pair  *v_ramp_up_pair = new LinearConstraint::v_coeff_pair(3);
			LinearConstraint::v_coeff_pair::iterator it;
			//Store the ramp down elements
			it = v_ramp_down_pair->begin();
			*it = LinearConstraint::coeff_pair( &P[i] , 1); //store p_{i, t}
			it = v_ramp_down_pair->begin()+1;
			*it = LinearConstraint::coeff_pair( &P[i+1] , -1); //store - p_{i,t+1}
			it = v_ramp_down_pair->begin()+2;
			*it = LinearConstraint::coeff_pair( &U[i+1] , (pmin-D_down) ); //store u_{i,t+1}*(p_min-D^down)

			//Store the ramp up elements
			it = v_ramp_up_pair->begin();
			*it = LinearConstraint::coeff_pair ( &P[i+1] , 1); //store p_{i, t+1}
			it = v_ramp_up_pair->begin()+1;
			*it = LinearConstraint::coeff_pair( &P[i] , -1); //store - p_{i,t}
			it = v_ramp_up_pair->begin()+2;
			*it = LinearConstraint::coeff_pair( &U[i] , (pmax-D_up) ); //store u_{i,t}*(pmax-D^up)

			//Initializing the Ramp Up & Down Constraints
				 

		 	RampDown_Const[i+1].add_variables( v_ramp_down_pair, false); 
			RampDown_Const[i+1].set_Block(this);
            RampDown_Const[i+1].set_lhs(-Inf<double>());
		    RampDown_Const[i+1].set_rhs(pmin);
				 
			RampUp_Const[i+1].add_variables( v_ramp_up_pair, false);
			RampUp_Const[i+1].set_Block(this);
			RampUp_Const[i+1].set_lhs(-Inf<double>());
			RampUp_Const[i+1].set_rhs(pmax);
		  
		}
	//Adding them to the Static Constraints Vector
	add_static_constraint(RampDown_Const);
	add_static_constraint(RampUp_Const);
}

void AcadThermalUnitBlock::StartUpCost(int downTime, double& cost, eShutStatus& shutStatus){

  if(data == 0){
	if (downTime>fGetStableSCT) // checking if the downtime of unit is greater than the maximum downtime cost
			downTime = fGetStableSCT;

		double coolCost;
		coolCost = coolAndFuelCost * (1-exp(-downTime/tau)) + fixedCost; //computing the cold start up

		double hotCost;
		hotCost = hotAndFuelCost * downTime + fixedCost; //computing the hot start up

		if (coolCost<hotCost)
		{
			cost = coolCost;
			shutStatus = eCooling;
		}
		else
		{
			cost = hotCost;
			shutStatus = eBancking;
		}
}
else{
cost = f_start_cost;
}


}

double AcadThermalUnitBlock::ComputeStableSCT()
{
   double tollerance = 0.01;
   fGetStableSCT = fMinDownTime;
   double t1cost=0.0;
   double t2cost=0.0;
   eShutStatus shutStatus;

   StartUpCost(fGetStableSCT, t1cost, shutStatus);
   StartUpCost(fGetStableSCT+1, t2cost, shutStatus);

   while ((t2cost-t1cost)/t1cost > tollerance)
   {
      assert (t2cost-t1cost > 0);

      t1cost = t2cost;
      StartUpCost(++fGetStableSCT, t2cost, shutStatus);
   }
   return fGetStableSCT;
}

void AcadThermalUnitBlock::obj_function( ){

if(model == 0){
UCBlock * f_UC_Block = dynamic_cast <UCBlock *> (get_f_Block());

	/*The objective function of the Thermal unit is of the form: f_p(t) = (p_i,t)^2 * a + p_i,t * b + u_i,t * c
	  and on top of that there are the start up costs of the unit that are calculated in the above methods */

	/*Ιnitilaize the vector of pairs of the variables and coefficients of the corresponding
	linear and diagonial quadratic part of the objective function */

if (pc == 0){
		LinearConstraint::v_coeff_pair  *v_l_pair = new LinearConstraint::v_coeff_pair((f_UC_Block->get_t())*2 + f_UC_Block->get_t() - init_t);
		LinearConstraint::v_coeff_pair  *v_q_pair = new LinearConstraint::v_coeff_pair(f_UC_Block->get_t());
	
		/* Pass the corresponding values to the vectors*/
		LinearConstraint::v_coeff_pair::iterator it;
		int k=0;
		double value=0;
		for(int i=0; i< f_UC_Block->get_t() ; i++){
			it = v_l_pair->begin()+k;k++;
			value = fLinearTherm - lambda[i];
			*it = LinearConstraint::coeff_pair( &P[i] , value); //store p_i,t*b
			it = v_l_pair->begin()+k;k++;
			value = fConstTherm - mew[i] * fMaxPower;
			*it = LinearConstraint::coeff_pair( &U[i], value); //store u_i,t*c
			it = v_q_pair->begin()+i;
			*it = LinearConstraint::coeff_pair( &P[i] , 2 * fQuadTherm); //store p_i,t*b
		}


	 //start_up_costs
	  if(data == 0){
		f_start_cost =0;
	        int a = ComputeStableSCT();
	
		std::vector<double> v_StartupCosts;
		v_StartupCosts.resize( a - fMinDownTime +1);
		eShutStatus shutStatus;
		for ( int i=0; i < a - fMinDownTime +1; ++i) {
		      double cost;
		      StartUpCost(i+fMinDownTime, cost,shutStatus);
		      v_StartupCosts[i] = cost;
  		}

		for ( int i = fMinDownTime ; i <= fGetStableSCT ; i ++ )
		      f_start_cost = v_StartupCosts[i-fMinDownTime] ;
	} //data == 0


	for (int i =0 ; i < f_UC_Block->get_t() - init_t ; ++i){
		it = v_l_pair->begin()+k;k++;
		*it = LinearConstraint::coeff_pair( &U_01[i] , f_start_cost); 
		}
	
       if(of==1){	      
       		obj_fun.add_variables(v_l_pair,false); //add the vector of the linear part of the obj function
		delete v_q_pair;//Tiziano
       		obj_fun.set_type(ObjectiveFunction::eMin);
       		set_objective_function(obj_fun); // add the objective funtion to the Block
       }
       else if (of==0){
       		qbj_fun.add_variables(v_l_pair,false); //add the vector of the linear part of the obj function
       		qbj_fun.add_q_variables(v_q_pair,false); //add the vector of the linear part of the obj function
       		qbj_fun.set_type(ObjectiveFunction::eMin);
       		set_objective_function(qbj_fun); // add the objective funtion to the Block
       }


}//if pc == 0
else if (pc == 1){

	//objective function is only linear containing the start up costs and the z_t
	LinearConstraint::v_coeff_pair  *v_l_pair = new LinearConstraint::v_coeff_pair(f_UC_Block->get_t() + f_UC_Block->get_t() - init_t);

	LinearConstraint::v_coeff_pair::iterator it;
	int k=0;

	for(int i=0; i< f_UC_Block->get_t() ; i++){
		it = v_l_pair->begin()+k;k++;
		*it = LinearConstraint::coeff_pair( &z_t[i] , 1.0); //store z_t
	}

	 //start_up_costs
	  if(data == 0){
		f_start_cost =0;
	        int a = ComputeStableSCT();
	
		std::vector<double> v_StartupCosts;
		v_StartupCosts.resize( a - fMinDownTime +1);
		eShutStatus shutStatus;
		for ( int i=0; i < a - fMinDownTime +1; ++i) {
		      double cost;
		      StartUpCost(i+fMinDownTime, cost,shutStatus);
		      v_StartupCosts[i] = cost;
  		}

		for ( int i = fMinDownTime ; i <= fGetStableSCT ; i ++ )
		      f_start_cost = v_StartupCosts[i-fMinDownTime] ;
	} //data == 0


	for (int i =0 ; i < f_UC_Block->get_t() - init_t ; ++i){
		it = v_l_pair->begin()+k;k++;
		*it = LinearConstraint::coeff_pair( &U_01[i] , f_start_cost); 
		}

	obj_fun.add_variables(v_l_pair,false); //add the vector of the linear part of the obj function
        obj_fun.set_type(ObjectiveFunction::eMin);
        set_objective_function(obj_fun); // add the objective funtion to the Block


} // if pc == 1
}

else if (model == 1){



UCBlock * f_UC_Block = dynamic_cast <UCBlock *> (get_f_Block());

	/*The objective function of the Thermal unit is of the form: f_p(t) = (p_i,t)^2 * a + p_i,t * b + u_i,t * c
	  and on top of that there are the start up costs of the unit that are calculated in the above methods */

	/*Ιnitilaize the vector of pairs of the variables and coefficients of the corresponding
	linear and diagonial quadratic part of the objective function */
int count=0;
for(int i = 0 ; i < Y_minus_pair.size() ; i++){

if (Y_minus_pair[i].second != f_UC_Block->get_t()+1) count++;

}

	LinearConstraint::v_coeff_pair  *v_l_pair = new LinearConstraint::v_coeff_pair(P_hk.size() + Y_plus.size() + count);
	LinearConstraint::v_coeff_pair  *v_q_pair = new LinearConstraint::v_coeff_pair(P_hk.size());

int l = 0;
int m = 0;
int z = 0;
LinearConstraint::v_coeff_pair::iterator it;
//cout<<endl<<"Echek obj fun:";
for(int i = 0 ; i < Y_plus.size() ; ++i){
//cout<<endl<<"For Yplus[" << i << "] of interval : (" << Y_plus_pair[i].first << " ; " << Y_plus_pair[i].second  << ") ";
  //  we proceed with adding first the constant part where we have u_i,t * c transform to y^+_{h,k} * \sum_t ={h,...,k} c_i * mew[t]
    double value = 0;
    double value2 = 0;
    z = 0;
    for (int j = Y_plus_pair[i].first ; j <= Y_plus_pair[i].second ; j++){
    //cout<<endl<<"Examining j: " << j;
        if( j != f_UC_Block->get_t()+1 && j != 0) {
            //cout<<endl<<"Passed and we have z : " << z << " && mew[ " << j-1 << "] = " << mew[j-1] << " && lambda[ " << j-1 << "] = " << lambda[j-1] << " && P_hk index: " << v_p_size[i]-v_p_ind[i]+z;
            value  += fConstTherm - mew[j-1] * fMaxPower;
            value2 = fLinearTherm - lambda[j-1];
            it = v_l_pair->begin()+l;l++;
           	*it = LinearConstraint::coeff_pair( &P_hk[v_p_size[i]-v_p_ind[i]+z]  , value2); //store p_i,t*b
            it = v_q_pair->begin()+m;m++;
    		*it = LinearConstraint::coeff_pair( &P_hk[v_p_size[i]-v_p_ind[i]+z] , 2 * fQuadTherm); //store p_i,t*b
            z++;
            value2 = 0;
        }
    }

    it = v_l_pair->begin()+l;l++;
    *it = LinearConstraint::coeff_pair( &Y_plus[i], value); //store y^+_{h,k} * \sum_t ={h,...,k} c_i * mew[t]


    //now we need to add the linear part of the power production p_i,t*b that is transformed

}

//Finally we need to also include the start up costs in here
for(int i = 0 ; i < Y_minus_pair.size() ; i++){

    if (Y_minus_pair[i].second != f_UC_Block->get_t()+1) {

        LinearConstraint::v_coeff_pair::iterator it;
        it = v_l_pair->begin()+l;l++;
        *it = LinearConstraint::coeff_pair( &Y_minus[i]  , f_start_cost); 
    }


}
       if(of==1){	      
       		obj_fun.add_variables(v_l_pair,false); //add the vector of the linear part of the obj function
		delete v_q_pair; //Tiziano
       		obj_fun.set_type(ObjectiveFunction::eMin);
       		set_objective_function(obj_fun); // add the objective funtion to the Block
       }
       else if (of==0){
       		qbj_fun.add_variables(v_l_pair,false); //add the vector of the linear part of the obj function
       		qbj_fun.add_q_variables(v_q_pair,false); //add the vector of the linear part of the obj function
       		qbj_fun.set_type(ObjectiveFunction::eMin);
       		set_objective_function(qbj_fun); // add the objective funtion to the Block
       }

 
}


else if (model == 2){

UCBlock * f_UC_Block = dynamic_cast <UCBlock *> (get_f_Block());

	//The objective function of the Thermal unit is of the form: f_p(t) = (p_i,t)^2 * a + p_i,t * b + u_i,t * c
	//  and on top of that there are the start up costs of the unit that are calculated in the above methods 

	//Ιnitilaize the vector of pairs of the variables and coefficients of the corresponding
	//linear and diagonial quadratic part of the objective function 

int count=0;
for(int i = 0 ; i < Y_minus_pair.size() ; i++){

if (Y_minus_pair[i].second != f_UC_Block->get_t()+1) count++;

}


	LinearConstraint::v_coeff_pair  *v_l_pair = new LinearConstraint::v_coeff_pair(P.size() + count + Y_plus.size());
	LinearConstraint::v_coeff_pair  *v_q_pair = new LinearConstraint::v_coeff_pair(P.size());

	// Pass the corresponding values to the vectors
	LinearConstraint::v_coeff_pair::iterator it;
	int k=0;
	double value=0;
	for(int i=0; i< f_UC_Block->get_t() ; i++){

		it = v_l_pair->begin()+k;k++;
		value = fLinearTherm - lambda[i];
		*it = LinearConstraint::coeff_pair( &P[i] , value); //store p_i,t*b
    	it = v_q_pair->begin()+i;
		*it = LinearConstraint::coeff_pair( &P[i] , 2 * fQuadTherm); //store p_i,t*b
	}



for(int i = 0 ; i < Y_plus_pair.size() ; i++){
  	    value = 0;
        for (int z = Y_plus_pair[i].first ; z <= Y_plus_pair[i].second ; z++){
            if( z != f_UC_Block->get_t()+1 && z != 0) {
                value = value + fConstTherm - mew[z-1] * fMaxPower;
            }}
    	it = v_l_pair->begin()+k;k++;
		*it = LinearConstraint::coeff_pair( &Y_plus[i], value); //store u_i,t*c
}

//Finally we need to also include the start up costs in here
for(int i = 0 ; i < Y_minus_pair.size() ; i++){

    if (Y_minus_pair[i].second != f_UC_Block->get_t()+1) {

        LinearConstraint::v_coeff_pair::iterator it;
        it = v_l_pair->begin()+k;k++;
        *it = LinearConstraint::coeff_pair( &Y_minus[i]  , f_start_cost); 
    }


} //for int i


       if(of==1){	  

       		obj_fun.add_variables(v_l_pair,false); //add the vector of the linear part of the obj function
		delete v_q_pair; //Tiziano
		obj_fun.set_type(ObjectiveFunction::eMin);
       		set_objective_function(obj_fun); // add the objective funtion to the Block
       }
       else if (of==0){

		qbj_fun.add_variables(v_l_pair,false); //add the vector of the linear part of the obj function
       		qbj_fun.add_q_variables(v_q_pair,false); //add the vector of the linear part of the obj function
       		qbj_fun.set_type(ObjectiveFunction::eMin);
       		set_objective_function(qbj_fun); // add the objective funtion to the Block
       }

} //if model 2


else if (model == 3){



UCBlock * f_UC_Block = dynamic_cast <UCBlock *> (get_f_Block());

	/*The objective function of the Thermal unit is of the form: f_p(t) = (p_i,t)^2 * a + p_i,t * b + u_i,t * c
	  and on top of that there are the start up costs of the unit that are calculated in the above methods */

	/*Ιnitilaize the vector of pairs of the variables and coefficients of the corresponding
	linear and diagonial quadratic part of the objective function */
int count=0;
for(int i = 0 ; i < Y_minus_pair.size() ; i++){

if (Y_minus_pair[i].second != f_UC_Block->get_t()+1) count++;

}

	LinearConstraint::v_coeff_pair  *v_l_pair = new LinearConstraint::v_coeff_pair(P_h.size() + Y_plus.size() + count);
	LinearConstraint::v_coeff_pair  *v_q_pair = new LinearConstraint::v_coeff_pair(P_h.size());

int l = 0;
LinearConstraint::v_coeff_pair::iterator it;
//we proceed with adding first the constant part where we have u_i,t * c transform to y^+_{h,k} * \sum_t ={h,...,k} c_i * mew[t]
for(int i = 0 ; i < Y_plus.size() ; ++i){
    double value = 0;
    for (int j = Y_plus_pair[i].first ; j <= Y_plus_pair[i].second ; j++){
        if( j != f_UC_Block->get_t()+1 && j != 0) value  += fConstTherm - mew[j-1] * fMaxPower;
        }//for loop
    it = v_l_pair->begin()+l;l++;
    *it = LinearConstraint::coeff_pair( &Y_plus[i], value); //store y^+_{h,k} * \sum_t ={h,...,k} c_i * mew[t]


}

int q = 0;
int z = 0;
int t = 0;
double value2 = 0;
//now we need to add the linear and the quadratic parts of the power production p_i,t*(b-lamda[t]) and p_i,t^2*a respectively
for (int j  = 0 ; j < v_ph_num.size() ; j++){
    	t = Y_plus_pair[v_ph_ind[j]].first;
	if(Y_plus_pair[v_ph_ind[j]].first == 0) t++;
    for(int s = 0 ; s < v_ph_size[j]; s++){
    	value2 = fLinearTherm - lambda[t-1];
    	it = v_l_pair->begin()+l;l++;
    	*it = LinearConstraint::coeff_pair( &P_h[z]  , value2); //store p_i,t*b
    	it = v_q_pair->begin()+q;q++;
    	*it = LinearConstraint::coeff_pair( &P_h[z] , 2 * fQuadTherm); //store p_i,t*a
    	z++;
    	t++;
    }

}
//cout<<endl<<"of_check: z = " << z << " && v_q_pair->size() = " << v_q_pair->size() << " &&  P_h.size() = " << P_h.size();
//Finally we need to also include the start up costs in here
for(int i = 0 ; i < Y_minus_pair.size() ; i++){

    if (Y_minus_pair[i].second != f_UC_Block->get_t()+1) {

        LinearConstraint::v_coeff_pair::iterator it;
        it = v_l_pair->begin()+l;l++;
        *it = LinearConstraint::coeff_pair( &Y_minus[i]  , f_start_cost); 
    }


}


       if(of==1){	      
       		obj_fun.add_variables(v_l_pair,false); //add the vector of the linear part of the obj function
		delete v_q_pair; //Tiziano
       		obj_fun.set_type(ObjectiveFunction::eMin);
       		set_objective_function(obj_fun); // add the objective funtion to the Block
       }
       else if (of==0){
       		qbj_fun.add_variables(v_l_pair,false); //add the vector of the linear part of the obj function
       		qbj_fun.add_q_variables(v_q_pair,false); //add the vector of the linear part of the obj function
       		qbj_fun.set_type(ObjectiveFunction::eMin);
       		set_objective_function(qbj_fun); // add the objective funtion to the Block
       }



}





}





int AcadThermalUnitBlock::mat_index (int i, int j){

int index=0;

UCBlock * f_UC_Block = dynamic_cast <UCBlock *> (get_f_Block());

if (i == f_UC_Block->get_t()+1 ) {
    index = Net_Matrix.size() -1;
    return index;
}

if (j==0){ // case of node (k,ON)
    for(int k = 0 ; k < v_plus.size() ; ++k)
        if (v_plus[k] == i) index = k+1;//since the first position is for the start node
}
else if (j==1){ // case of node (k,OFF)
    for(int k = 0 ; k < v_minus.size() ; ++k)
        if (v_minus[k] == i) index = k+1+v_plus.size();//since is after the start and the valid ON nodes
}


return index;
}





AcadThermalUnitBlock::_init::_init()
{
 UnitBlock::U_factory()["NumThermal"]
  = boost::bind(boost::factory<AcadThermalUnitBlock*>(), _1);
 }

AcadThermalUnitBlock::_init AcadThermalUnitBlock::_initializer;
 /* Ensure that the static member _initializer is initialized, so that
    AcadThermalUnitBlock is registered into UnitBlock::f_factory. */
