/*--------------------------------------------------------------------------*/
/*------------------ File AcadThermalUnitMIPBlock.cpp ----------------------*/
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

#include "AcadThermalUnitMIPBlock.h"
#include "UCBlock.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

using namespace std;

/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

AcadThermalUnitMIPBlock::AcadThermalUnitMIPBlock(UCBlock * fblock) : Block(fblock) {
	
	rampconst = false;
}

/*--------------------------------------------------------------------------*/

AcadThermalUnitMIPBlock::~AcadThermalUnitMIPBlock() {
	// TODO Auto-generated destructor stub
}

/*--------------------------------------------------------------------------*/

void AcadThermalUnitMIPBlock::load(std::istream& inStream){

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
	 }

v_p_size.resize(count_plus);
v_p_ind.resize(count_plus);
//int coun=0;
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

add_static_variable(Y_plus);  
add_static_variable(Y_minus);
add_static_variable(P_hk);    

/*cout<<endl<<"So we have Y_plus.size(): " << Y_plus.size() << "with the plus pairs: " << endl;

for(int i = 0; i < Y_plus_pair.size() ; i++){

cout<< " [ " << Y_plus_pair[i].first << " ; " << Y_plus_pair[i].second << " ]";

}

cout<<endl<<"And then we have the Y_minus.size(): " << Y_minus.size() <<" with the minus pairs: " << endl;

for(int i = 0; i < Y_minus_pair.size() ; i++){

cout<< " [ " << Y_minus_pair[i].first << " ; " << Y_minus_pair[i].second << " ]";

}*/



/*cout<<endl<<"So we have v_minus: ";
for(int i = 0; i < v_minus.size() ; i++){

cout<< v_minus[i] <<  "  ";

}

cout<<endl<<"So we have v_plus: ";
for(int i = 0; i < v_plus.size() ; i++){

cout<< v_plus[i] <<  "  ";

}
cout<<endl;
cout<<endl<<"Also we have v_p_size: ";
for(int i = 0; i < v_p_size.size() ; i++){

cout<< v_p_size[i] <<  "  ";

}
cout<<endl<<"Also we have v_p_ind: ";
for(int i = 0; i < v_p_ind.size() ; i++){

cout<< v_p_ind[i] <<  "  ";

}
cout<<endl;*/
//add objective function
obj_function( );

}

void AcadThermalUnitMIPBlock::Network_Matrix(){

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

void AcadThermalUnitMIPBlock::generate_static_constraints ( ){

//call and run method for network matrix construction and usage
Network_Matrix();

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
    if (i == 0 ) *it = LinearConstraint::coeff_pair( &P_hk[0] , 1); //store p^h_{h, k}
    else *it = LinearConstraint::coeff_pair( &P_hk[v_p_size[i]-v_p_ind[i]] , 1); //store p^h_{h, k}
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
//First we need to check the size, since Bound_on constraints apply only for h->k routes that finish
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

if(v_p_ind[i] > 2) {
    if (Y_plus_pair[i].second == f_UC_Block->get_t()+1) count = count +  v_p_ind[i] - 1;
    if (Y_plus_pair[i].second != f_UC_Block->get_t()+1) count = count +  v_p_ind[i] - 2;
}
else if (v_p_ind[i] == 2 && (Y_plus_pair[i].first == 0 || Y_plus_pair[i].second == f_UC_Block->get_t()+1 )) count++;

}
//cout<<endl<<"yolo && count = " << count;
Pmax_Const.resize(count);

j = 0;
for (int i = 0 ; i < Y_plus.size() ; ++i){

if(v_p_ind[i] > 2){
int z = 0;
if (Y_plus_pair[i].second != f_UC_Block->get_t()+1) z=2;
else if (Y_plus_pair[i].second == f_UC_Block->get_t()+1) z=1;
    for (int k = 0 ; k < v_p_ind[i]-z ; k++){

        LinearConstraint::v_coeff_pair  *v_max_pair = new LinearConstraint::v_coeff_pair(2);
        LinearConstraint::v_coeff_pair::iterator it;
        //Store the max power elements
	    it = v_max_pair->begin();
        *it = LinearConstraint::coeff_pair( &P_hk[v_p_size[i]-v_p_ind[i]+1+k], 1); 
                                        //store p^t_{h, k} for all h<t<k and we have v_p_size[i]-v_p_ind[i] giving us position h, 
                                        //so we need to add plus 1 to arrive to the first element that pma applies and the we continue with the counter
    	it = v_max_pair->begin()+1;
        *it = LinearConstraint::coeff_pair( &Y_plus[i], -fMaxPower ); //store - p_max * y^+_{h,k}

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

//Construct the Ramp Up   Constraints for all different nodes: p^{t+1}_{h, k} - p^{t}_{h, k}     - D^+ * y^+_{h,k}  <= 0 \for t = { h , ... , k-1 }
//Construct the Ramp Down Constraints for all different nodes: p^t_{h, k}     - p^{t+1}_{h, k}   - D^- * y^+_{h,k}  <= 0 \for t = { h , ... , k-1 }
// both of them are only valid for routes bigger than one time step
count = 0;

for(int i = 0 ; i < Y_plus.size() ; ++i){

if (Y_plus_pair[i].first == 0) count = count + 1;
if(v_p_ind[i] > 1) count = count +  v_p_ind[i] - 1;

}

RampDown_Const.resize(count);
RampUp_Const.resize(count);
//cout<<endl<<"yo count = " << count;
j=0;
for(int i = 0 ; i < Y_plus.size() ; ++i){

//create initial ramp consts
if (Y_plus_pair[i].first == 0){

//cout<<endl<< " v_p_size[i] =  " << v_p_size[i] << " && v_p_ind[i] =  " << v_p_ind[i] << " && index = " << v_p_size[i]-v_p_ind[i] << " && P_hk.size() = " <<P_hk.size();  
        LinearConstraint::v_coeff_pair  *v_ramp_up_pair = new LinearConstraint::v_coeff_pair(2);
        LinearConstraint::v_coeff_pair::iterator it;
	    it = v_ramp_up_pair->begin();
        *it = LinearConstraint::coeff_pair( &P_hk[v_p_size[i]-v_p_ind[i]] , - 1);
    	it = v_ramp_up_pair->begin()+1;
        *it = LinearConstraint::coeff_pair( &Y_plus[i], (-fMaxRampDown+fInitPower) );  //fBoundDown-fMaxRampDown

        RampUp_Const[j].add_variables( v_ramp_up_pair, false);
    	RampUp_Const[j].set_Block(this);
    	RampUp_Const[j].set_lhs(-Inf<double>());
    	RampUp_Const[j].set_rhs(0); // fBoundDown - fInitPower

        LinearConstraint::v_coeff_pair  *v_ramp_down_pair = new LinearConstraint::v_coeff_pair(1);
        LinearConstraint::v_coeff_pair::iterator it2;
	    it2 = v_ramp_down_pair->begin();
        *it2 = LinearConstraint::coeff_pair( &P_hk[v_p_size[i]-v_p_ind[i]] , 1);
    	
  		    double initial_rhs = (fInitUpDownTime > 0 );
			double rhs;
			rhs = ( initial_rhs ?  - fInitPower - fMaxRampUp :  0 );	


        RampDown_Const[j].add_variables( v_ramp_down_pair, false);
    	RampDown_Const[j].set_Block(this);
    	RampDown_Const[j].set_lhs(-Inf<double>());
    	RampDown_Const[j].set_rhs(-rhs);
        
        j++;
}

if(v_p_ind[i] > 1){


   for (int k = 0 ; k < v_p_ind[i] - 1 ; k++){
            /********* Ramp_Up Constraints *********/
        LinearConstraint::v_coeff_pair  *v_ramp_up_pair = new LinearConstraint::v_coeff_pair(3);
        LinearConstraint::v_coeff_pair::iterator it;
  //      cout<<endl<< " v_p_size[i] =  " << v_p_size[i] << " && v_p_ind[i] =  " << v_p_ind[i] << " && k = " << k << " && index = " << v_p_size[i]-v_p_ind[i]+k+1 << " && P_hk.size() = " <<P_hk.size();
  //      cout<<endl<< " ... " ;
        //Store the Ramp_Up elements
	    it = v_ramp_up_pair->begin();
        *it = LinearConstraint::coeff_pair( &P_hk[v_p_size[i]-v_p_ind[i]+k+1], 1); 
                                        //store p^t+1_{h, k} for all h<=t<k and we have v_p_size[i]-v_p_ind[i] giving us position h, 
                                        // and then we continue with the counter, since we want the t+1 element we add plus one to the 
    	it = v_ramp_up_pair->begin()+1;
        *it = LinearConstraint::coeff_pair( &P_hk[v_p_size[i]-v_p_ind[i]+k] , - 1);
                                         //store p^t_{h, k} for all h<=t<k and we have v_p_size[i]-v_p_ind[i] giving us position h, 
                                        //so we start from there as first element that applies ramp constraints and the we continue with the counter
    	it = v_ramp_up_pair->begin()+2;
        *it = LinearConstraint::coeff_pair( &Y_plus[i], -fMaxRampUp ); //store - y^+_{h,k} * D^+
    //    cout<<endl << "fMaxRampUp = " << fMaxRampUp << "&& j = " << j;
     //   cout<<endl << "fMaxRampUp = " << fMaxRampUp;
        RampUp_Const[j].add_variables( v_ramp_up_pair, false);
    	RampUp_Const[j].set_Block(this);
    	RampUp_Const[j].set_lhs(-Inf<double>());
    	RampUp_Const[j].set_rhs(0);


        /********* Ramp_Down Constraints *********/
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


void AcadThermalUnitMIPBlock::obj_function( ){


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


 obj_fun.add_q_variables(v_q_pair,false); //add the vector of the linear part of the obj function
 obj_fun.add_variables(v_l_pair,false); //add the vector of the linear part of the obj function
 obj_fun.set_type(ObjectiveFunction::eMin);
 set_objective_function(obj_fun); // add the objective funtion to the Block

   

}

int AcadThermalUnitMIPBlock::mat_index (int i, int j){

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


