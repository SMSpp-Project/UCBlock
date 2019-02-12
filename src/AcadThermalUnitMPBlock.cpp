/*--------------------------------------------------------------------------*/
/*------------------ File AcadThermalUnitMPBlock.cpp ----------------------*/
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
#include "boost/multi_array.hpp"
#include <iostream>

#include "AcadThermalUnitMPBlock.h"
#include "UCBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

using namespace std;

/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

AcadThermalUnitMPBlock::AcadThermalUnitMPBlock(UCBlock * fblock) : Block(fblock) {
	
	rampconst = false;
}

/*--------------------------------------------------------------------------*/

AcadThermalUnitMPBlock::~AcadThermalUnitMPBlock() {
	// TODO Auto-generated destructor stub
}

/*--------------------------------------------------------------------------*/

void AcadThermalUnitMPBlock::load(std::istream& inStream){

UCBlock * f_UC_Block = dynamic_cast <UCBlock *> (get_f_Block());

	string skip;
	data = 0;

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
//inStream >>	skip;
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

  if ( fInitUpDownTime > 0 )
    {
    
      init_t =  ( fInitUpDownTime >= fMinUpTime ? /*0*/1 : fMinUpTime   - fInitUpDownTime );  // if ( a > b ? 1 : 2) => if a>b -> 1 || b<=a ->2
    }
  else
    {
      init_t = ( -fInitUpDownTime >= fMinDownTime ? /*0*/1 : fMinDownTime + fInitUpDownTime ); 
    
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
P.resize(f_UC_Block->get_t());

 for(int i = 0 ; i<P.size() ; i++){
		 //set pointer of the father Block
		 P[i].set_Block(this); 
		 //set lower and upper bounds
		 P[i].set_lb(0.0); 
 		 P[i].set_ub(Inf<double>()); 
		 //set type of variables
		 P[i].set_type(ColVariable::continuous);
	 }
//cout<<endl<<"Y_minus with size: " << Y_minus.size() << " && " << Y_minus_pair.size();
 for(int i = 0 ; i<Y_minus.size() ; i++){
		 //set pointer of the father Block
		 Y_minus[i].set_Block(this); 
		 //set lower and upper bounds
		 Y_minus[i].set_lb(0.0); 
 		 Y_minus[i].set_ub(1.0); 
		 //set type of variables
		 Y_minus[i].set_type(ColVariable::binary);
         //cout<<endl<<"Y_minus[ " << i << " ] =  ( " << Y_minus_pair[i].first << " ; " << Y_minus_pair[i].second << " )";
	 }
//cout<<endl<<"Y_plus with size: " << Y_plus.size() << " && " << Y_plus.size();
 for(int i = 0 ; i<Y_plus.size() ; i++){
		 //set pointer of the father Block
		 Y_plus[i].set_Block(this); 
		 //set lower and upper bounds
		 Y_plus[i].set_lb(0.0); 
 		 Y_plus[i].set_ub(1.0); 
		 //set type of variables
		 Y_plus[i].set_type(ColVariable::binary);
         //cout<<endl<<"Y_plus[ " << i << " ] =  ( " << Y_plus_pair[i].first << " ; " << Y_plus_pair[i].second << " )";
	 }
add_static_variable(P);    
add_static_variable(Y_plus);  
add_static_variable(Y_minus);


v_t.resize(f_UC_Block->get_t() /*- init_t*/);


for(int i = 0; i < v_t.size() ; i++){
    for(int j = 0 ; j < Y_plus_pair.size(); j++){
    
        if (Y_plus_pair[j].first <= i + 1 && i + 1 <= Y_plus_pair[j].second ){
            v_t[i].index.push_back(j);
            if ( i + 2 <= Y_plus_pair[j].second )
                v_t[i].sindex.push_back(j);
            if (Y_plus_pair[j].first == i + 1 /*&& Y_plus_pair[j].first != 0 */)
                v_t[i].b_on.push_back(j);
            if (Y_plus_pair[j].second == i + 1 /*&& Y_plus_pair[j].second != f_UC_Block->get_t() + 1 */)
                v_t[i].b_dn.push_back(j);
            if (Y_plus_pair[j].first != i + 1 &&  Y_plus_pair[j].second != i + 1)
                v_t[i].p_mx.push_back(j);
    }//first if
      
    }


}

//add objective function
obj_function( );

}

void AcadThermalUnitMPBlock::Network_Matrix(){

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
        /*cout<<endl <<"YOLOOOOOOO_k = " << k << " FOOOOR Y_plus_pair[i].first = " << Y_plus_pair[i].first;
        cout<<endl <<"//";*/

        Net_Matrix[k].add_variables( v_t_pair, false); 

        LinearConstraint::v_coeff_pair * v_s_pair = new LinearConstraint::v_coeff_pair(1);
        LinearConstraint::v_coeff_pair::iterator it2;
	    it2 = v_s_pair->begin();
	    *it2 = LinearConstraint::coeff_pair( &Y_plus[i] , -1.0); 
        int h = mat_index(Y_plus_pair[i].second,1);
        /*cout<<endl <<"YOLOOOOOOO_h = " << h << " FOOOOR Y_plus_pair[i].second = " << Y_plus_pair[i].second;
        cout<<endl <<"//";*/

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

void AcadThermalUnitMPBlock::generate_static_constraints ( ){

//call and run method for network matrix construction and usage
Network_Matrix();

//UCBlock * f_UC_Block = dynamic_cast <UCBlock *> (get_f_Block());
int j;

//Construct the Bound_Down Constraints for all different nodes: p^k_{h, k} - y^+_{h,k}*Bound_down <= 0
//First we need to check the size, since Bound_on constraints apply only for h->k routes that finish
//within the optimisation horizon

int count = 0;
for(int i = 0 ; i < v_t.size() ; i++){
    if(v_t[i].b_dn.size() > 0 || v_t[i].b_on.size() > 0 || v_t[i].p_mx.size() > 0) count++;
}
Bound_Down_Const.resize(count);
int l=0;
int k=0;
for(int i = 0; i < v_t.size() ; i++){
    k=0;
   if(v_t[i].b_dn.size() > 0 || v_t[i].b_on.size() > 0 || v_t[i].p_mx.size() > 0) {
       /********* Bound_Down Constraints *********/
        LinearConstraint::v_coeff_pair  *v_bound_down_pair = new LinearConstraint::v_coeff_pair(1+v_t[i].b_dn.size() + v_t[i].b_on.size() + v_t[i].p_mx.size() );
        LinearConstraint::v_coeff_pair::iterator it2;
        //Store the ramp down elements
	    it2 = v_bound_down_pair->begin();
        *it2 = LinearConstraint::coeff_pair( &P[i] , 1); //store p^k_{h, k}
        if (v_t[i].b_dn.size() > 0){
        for(j = 0 ; j < v_t[i].b_dn.size() ; j++){
	    it2 = v_bound_down_pair->begin()+1+k;k++;
	    *it2 = LinearConstraint::coeff_pair( &Y_plus[v_t[i].b_dn[j]] , -fBoundDown); //store - y^+_{h,k}*Bound_Down
        } }

        if (v_t[i].b_on.size() > 0){
             for(j = 0 ; j < v_t[i].b_on.size() ; j++){
	            it2 = v_bound_down_pair->begin()+1+k;k++;
	            *it2 = LinearConstraint::coeff_pair( &Y_plus[v_t[i].b_on[j]] , -fBoundOn); //store - y^+_{h,k}*Bound_Down
        }

        }

        if (v_t[i].p_mx.size() > 0){
            for(j = 0 ; j < v_t[i].p_mx.size() ; j++){
	            it2 = v_bound_down_pair->begin()+1+k;k++;
	            *it2 = LinearConstraint::coeff_pair( &Y_plus[v_t[i].p_mx[j]] , -fMaxPower); //store - y^+_{h,k}*Bound_Down
        }

        }

	    Bound_Down_Const[l].add_variables( v_bound_down_pair, false);
	    Bound_Down_Const[l].set_Block(this);
	    Bound_Down_Const[l].set_lhs(-Inf<double>());
	    Bound_Down_Const[l].set_rhs(0);
        l++;
        }
	
}
//Adding them to the Static Constraints Vector
add_static_constraint(Bound_Down_Const);

//Construct the Min Power Constraints for all different nodes: p_min * y^+_{h,k} - p^t_{h, k} <= 0 \for t = { h, ... , k}
Pmin_Const.resize(v_t.size());
j = 0;
for(int i = 0 ; i < Pmin_Const.size() ; ++i){

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




//Construct the Ramp Up   Constraints for all different nodes: p^{t+1}_{h, k} - p^{t}_{h, k}     - D^+ * y^+_{h,k}  <= 0 \for t = { h , ... , k-1 }
//Construct the Ramp Down Constraints for all different nodes: p^t_{h, k}     - p^{t+1}_{h, k}   - D^- * y^+_{h,k}  <= 0 \for t = { h , ... , k-1 }
// both of them are only valid for routes bigger than one time step
/*int count = 0;
for(int i = 0 ; i < v_t.size() ; i++){
    if(v_t[i].sindex.size() > 0) count++;
}
if(v_t[i].sindex.size() > 0)*/
RampDown_Const.resize(v_t.size());
RampUp_Const.resize(v_t.size());



			LinearConstraint::v_coeff_pair  *v_ramp_init_down_pair = new LinearConstraint::v_coeff_pair(1);
			LinearConstraint::v_coeff_pair  *v_ramp_init_up_pair = new LinearConstraint::v_coeff_pair(1+ v_t[0].index.size() );
			LinearConstraint::v_coeff_pair::iterator it;
			//Store the ramp down elements
			it = v_ramp_init_down_pair->begin();
			*it = LinearConstraint::coeff_pair( &P[0] , 1); //store p_{i, t}

			//Store the ramp up elements
			it = v_ramp_init_up_pair->begin();
			*it = LinearConstraint::coeff_pair( &P[0] , -1); //store - p_{i,t}
            for(int i = 0 ; i < v_t[0].index.size() ; i++  ){
    			it = v_ramp_init_up_pair->begin()+1+i;
            	*it = LinearConstraint::coeff_pair ( &Y_plus[v_t[0].index[i]] , (fBoundDown-fMaxRampDown) ); //store p_{i, t+1}
            }
			
			//Initializing the Ramp Up & Down Constraints
  		    double initial_rhs = (fInitUpDownTime > 0 );
			double rhs;

			rhs = ( initial_rhs ?  - fInitPower - fMaxRampUp : 0 );	//( initial_rhs ?  - fInitPower - fMaxRampUp : -  fBoundOn );	

		 	RampDown_Const[0].add_variables( v_ramp_init_down_pair, false); 
			RampDown_Const[0].set_Block(this);
            RampDown_Const[0].set_lhs(-Inf<double>());
		    RampDown_Const[0].set_rhs(-rhs);
				 
			rhs = fBoundDown - fInitPower ;
					
			RampUp_Const[0].add_variables( v_ramp_init_up_pair, false);
			RampUp_Const[0].set_Block(this);
			RampUp_Const[0].set_lhs(-Inf<double>());
			RampUp_Const[0].set_rhs(rhs);



j=0;
//LinearConstraint::v_coeff_pair::iterator it;
//LinearConstraint::v_coeff_pair::iterator it2;
for(int i = 0 ; i < RampUp_Const.size()-1 ; ++i){


           /********* Ramp_Up Constraints *********/
        LinearConstraint::v_coeff_pair  *v_ramp_up_pair = new LinearConstraint::v_coeff_pair(2 + v_t[i].sindex.size());


        //Store the Ramp_Up elements
	    auto it = v_ramp_up_pair->begin();
        *it = LinearConstraint::coeff_pair( &P[i+1], 1); 
                                        //store p^t+1_{h, k} for all h<=t<k and we have v_p_size[i]-v_p_ind[i] giving us position h, 
                                        // and then we continue with the counter, since we want the t+1 element we add plus one to the 
    	it = v_ramp_up_pair->begin()+1;
        *it = LinearConstraint::coeff_pair( &P[i] , - 1);
                                         //store p^t_{h, k} for all h<=t<k and we have v_p_size[i]-v_p_ind[i] giving us position h, 
                                        //so we start from there as first element that applies ramp constraints and the we continue with the counter
        int m = 2;
        for(j = 0 ; j < v_t[i].sindex.size() ;j++ ){
    	it = v_ramp_up_pair->begin()+m;
        m++;
        *it = LinearConstraint::coeff_pair( &Y_plus[v_t[i].sindex[j]], -fMaxRampUp + fBoundOn ); //store - y^+_{h,k} * D^+
    	   //    cout<<endl << "fMaxRampUp = " << fMaxRampUp << "&& j = " << j;
           //   cout<<endl << "fMaxRampUp = " << fMaxRampUp;
        }
  		rhs = ( m > 2 ?  fBoundOn : 0 );	
        RampUp_Const[i+1].add_variables( v_ramp_up_pair, false);
    	RampUp_Const[i+1].set_Block(this);
    	RampUp_Const[i+1].set_lhs(-Inf<double>());
    	RampUp_Const[i+1].set_rhs(rhs);


        /********* Ramp_Down Constraints *********/
        LinearConstraint::v_coeff_pair  *v_ramp_down_pair = new LinearConstraint::v_coeff_pair(2 + v_t[i].sindex.size());

        //Store the Ramp_Down elements
	    auto it2 = v_ramp_down_pair->begin();
        *it2 = LinearConstraint::coeff_pair( &P[i] , 1); //store p^{t}_{h, k}
    	it2 = v_ramp_down_pair->begin()+1;
        *it2 = LinearConstraint::coeff_pair( &P[i+1] , - 1); //store -p^{t+1}_{h, k}
        m=2;
        for(j = 0 ; j < v_t[i].sindex.size() ;j++ ){
    	it2 = v_ramp_down_pair->begin()+m;m++;
        *it2 = LinearConstraint::coeff_pair( &Y_plus[v_t[i].sindex[j]], -fMaxRampDown + fBoundDown ); //store - y^+_{h,k} * D^-
        }
        rhs = ( m > 2 ?  fBoundDown : 0 );	
        RampDown_Const[i+1].add_variables( v_ramp_down_pair, false);
        RampDown_Const[i+1].set_Block(this);
        RampDown_Const[i+1].set_lhs(-Inf<double>());
        RampDown_Const[i+1].set_rhs(rhs);

        }   

//Adding them to the Static Constraints Vector
add_static_constraint(RampUp_Const);
add_static_constraint(RampDown_Const);


}


void AcadThermalUnitMPBlock::obj_function( ){


UCBlock * f_UC_Block = dynamic_cast <UCBlock *> (get_f_Block());

	/*The objective function of the Thermal unit is of the form: f_p(t) = (p_i,t)^2 * a + p_i,t * b + u_i,t * c
	  and on top of that there are the start up costs of the unit that are calculated in the above methods */

	/*Ιnitilaize the vector of pairs of the variables and coefficients of the corresponding
	linear and diagonial quadratic part of the objective function */
int count=0;
for(int i = 0 ; i < Y_minus_pair.size() ; i++){

if (Y_minus_pair[i].second != f_UC_Block->get_t()+1) count++;

}


	LinearConstraint::v_coeff_pair  *v_l_pair = new LinearConstraint::v_coeff_pair(P.size() + count + Y_plus.size());
	LinearConstraint::v_coeff_pair  *v_q_pair = new LinearConstraint::v_coeff_pair(P.size());


	/* Pass the corresponding values to the vectors*/
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


}


 obj_fun.add_q_variables(v_q_pair,false); //add the vector of the linear part of the obj function
 obj_fun.add_variables(v_l_pair,false); //add the vector of the linear part of the obj function
 obj_fun.set_type(ObjectiveFunction::eMin);
 set_objective_function(obj_fun); // add the objective funtion to the Block

   

}

int AcadThermalUnitMPBlock::mat_index (int i, int j){

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
