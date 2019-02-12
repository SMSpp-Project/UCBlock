/*--------------------------------------------------------------------------*/
/*-------------------- File AcadThemalUnitGraphSolver.cpp ------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the AcadThemalUnitGraphSolver class.
 *
 * \version 0.10
 *
 * \date 24 - 09 - 2016
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

#include "AcadThemalUnitGraphSolver.h"
#include "AcadThermalUnitBlock.h"
#include "UCBlock.h"
#include "Block.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

void AcadThemalUnitGraphSolver::set_Block( Block *block ){


  if( f_Block )                         // was attached to some other Block
  f_Block->unregister_Solver( this );  // no more so
 f_Block = block;                      // this is the new block now
 f_Block->register_Solver( this );     // register to it


 AcadThermalUnitBlock * thermalUnit = dynamic_cast <AcadThermalUnitBlock *> (f_Block);

 fMaxStartupLevel = thermalUnit->ComputeStableSCT();
 

startup_costs.resize( fMaxStartupLevel - thermalUnit->get_fMinDownTime() + 1);

for (int i = 0; i < startup_costs.size(); i++ ){

     double cost;
     AcadThermalUnitBlock::eShutStatus shutStatus;
     thermalUnit->StartUpCost(i+thermalUnit->get_fMinDownTime(), cost, shutStatus);
     startup_costs[i] = cost;
}

     build_graph();


}

/*--------------------------------------------------------------------------*/

void AcadThemalUnitGraphSolver::build_graph( void )
{


 double eps   = 1e-10 ;  
 AcadThermalUnitBlock * thermalUnit = dynamic_cast <AcadThermalUnitBlock *> (f_Block);
 UCBlock *f_UCBlock = static_cast<UCBlock*>( thermalUnit->get_f_Block() );


    //setting the proper size of the different vectors used
    v_cost.resize(f_UCBlock->get_t());
    v_route.resize(f_UCBlock->get_t()+1);
    v_nodes.resize(f_UCBlock->get_t() * /*f_UCBlock->get_t()*/2); //change_1

    /* In order to construct the graph properly we need to differntiate between
       the following three cases:
      
       1) The unit was initially on and needs to remain on for some time steps
          due to ramp and/or min_up time constratins
       2) The unit was initially off and needs to remain off for some time steps
          due to min_down time constraint
       3) The unit regardless is initial on/off state is not subjected to any con-
          straint to keep its initial state for any time-steps and as a result 
          can freely decide from the begin its status
      */


     kMin = 0 ;

        if ( thermalUnit->get_fInitPower() >= thermalUnit->get_fBoundDown() + eps )
           {
                double tmp = thermalUnit->get_fInitPower();
                tmp -= thermalUnit->get_fMaxRampDown() ;
                kMin ++ ;
                while ( tmp >= thermalUnit->get_fBoundDown() + eps )
                {
                    tmp -= thermalUnit->get_fMaxRampDown() ;
                    kMin ++ ;
                }
                kMin -- ;
            }
            if ( kMin >= f_UCBlock->get_t() ) kMin = f_UCBlock->get_t() - 1 ;

               /** In this case the unit is on in the first time step and thus we are
                in node (0,1) from where we need to construct all arc connections,
                which for this node are:
                - all nodes (k,0), (k+1,0) , ... , (t,0) where k is the first time 
                  instance in which the unit is able to be switched off
                - node (t,1), which is the case that the unit is not switched of at
                  all and remains on throughout the whole optimisation period  */

    /************************************Case_1**********************************/

    if ( thermalUnit->get_fInitUpDownTime() > 0 && ( thermalUnit->get_fInitUpDownTime() < thermalUnit->get_fMinUpTime() || kMin > 0 )  )
    {   // checking that unit needs to remain on

        //calculating the time-steps (if any) for which unit needs to
        // remain on due to ramp_constraints
       
            int h= hmin = 0;
            /* calculation of the time instances for which unit needs to remain on
               due to min_up constraint */
            int k= thermalUnit->get_fMinUpTime()-thermalUnit->get_fInitUpDownTime()-1;
            
            /* setting the total time-steps, for which the unit needs to remain on
               based on the min_up & ramp constraints */         
            if ( kMin < k ) kMin = k;
            else            k    = kMin ;

            
           {

                double c_i=0;
                for( int t = h ; t<k; ++t) c_i += thermalUnit->get_fConstTherm();
                
                v_nodes[0].v_arcs.resize(f_UCBlock->get_t()-k);
                int i=0;
                //arcs to connect node (0,1) with nodes (k,0), (k+1,0) , ... , (t,0)
                for (; k< f_UCBlock->get_t() ; k ++ )
                {
                    v_nodes[0].v_arcs[i].h= h;
                    v_nodes[0].v_arcs[i].k= k/*+f_UCBlock->get_t()*/;
                    c_i += thermalUnit->get_fConstTherm();
                    v_nodes[0].v_arcs[i].cost1= c_i ;
                    v_nodes[0].v_arcs[i].cost2= 0;
                    v_nodes[0].v_arcs[i].valid= -1;//from on->off
                    i++;
                    
                }
                //arc to connect node (0,1) with node (t,1)
             

            }
            

            /*After the constructions of all arc connections for the node of first time
              step we need to create all the different arc conncetions for all the rest
              of the nodes, where we have:

              1) The nodes that denote the unit being offline, in which case we start to
                 create arc connections for the nodes: (k,0), (k+1,0) , ... , (t,0). Sin-
                 ce the unit could not have been switched off before. For this connections 
                 we need to take into account the following two cases:
              
                 i)  For all the time instances that are smaller than t-min_down_time, where
                     each (j,0) node is connected with all the nodes (j+min_down_time,1)...(t,1)
                     to denote the possible starts of the unit with respect with the min_down
                     constraint and also the arc connection with (t,0) which represents the
                     case in which the unit remains off for all the rest of the optimisation
                     period.
               
                 ii) For all the time instances that are larger/equal than t-min_down_time, 
                     where the unit is not able to switch its status inside this optimisation
                     period due to min_down constraint and thus the unit will remain ofline
                     for the rest time steps and the only arc connection needed is (t,0) */
                 
            for ( h = kMin; h < f_UCBlock->get_t() ; h ++ ) {
                
                if(h > f_UCBlock->get_t()-1-thermalUnit->get_fMinDownTime()){
                    v_nodes[h+f_UCBlock->get_t()].v_arcs.resize(1);
                    v_nodes[h+f_UCBlock->get_t()].v_arcs[0].h= h;
                    v_nodes[h+f_UCBlock->get_t()].v_arcs[0].k= f_UCBlock->get_t()-1;
                    v_nodes[h+f_UCBlock->get_t()].v_arcs[0].cost1= 0 ;
                    v_nodes[h+f_UCBlock->get_t()].v_arcs[0].cost2= 0;
                    v_nodes[h+f_UCBlock->get_t()].v_arcs[0].valid= -2; //from off->off
                    
                }
                else{
                    k = h + thermalUnit->get_fMinDownTime() - 1 ;
                    v_nodes[h+f_UCBlock->get_t()].v_arcs.resize(f_UCBlock->get_t()-k+1);
                    int i = 0;
                    for (; k < f_UCBlock->get_t() ; k ++ )
                        {
                        v_nodes[h+f_UCBlock->get_t()].v_arcs[i].h= h/*+f_UCBlock->get_t()*/;
                        v_nodes[h+f_UCBlock->get_t()].v_arcs[i].k= k;
                        v_nodes[h+f_UCBlock->get_t()].v_arcs[i].cost1= 0;
                        v_nodes[h+f_UCBlock->get_t()].v_arcs[i].cost2= 0; //start-up cost
                        v_nodes[h+f_UCBlock->get_t()].v_arcs[i].valid= 1; //from off->on
                        i++;
                        }	    
                    v_nodes[h+f_UCBlock->get_t()].v_arcs[i].h= h/*+f_UCBlock->get_t()*/;
                    v_nodes[h+f_UCBlock->get_t()].v_arcs[i].k= /*f_UCBlock->get_t()+*/f_UCBlock->get_t()-1;
                    v_nodes[h+f_UCBlock->get_t()].v_arcs[i].cost1= 0 ;
                    v_nodes[h+f_UCBlock->get_t()].v_arcs[i].cost2= 0;
                    v_nodes[h+f_UCBlock->get_t()].v_arcs[i].valid= -2; //from off->off
                }//else
                

                

            } //for loop h

            /*
              2) The nodes that denote the unit being online, in which case we start to
                 create arc connections for the nodes:
                                     (k+mindowntime,0), (k+mindowntime+1,0) , ... , (t,0). 
                 Since the unit could not have been switched on again at any possible time 
                 before, after its first. 
                 For these connections we need as above to take into account the following 
                 two cases:
              
                 i)  For all the time instances that are smaller than t-min_up_time, where
                     each (j,1) node is connected with all the nodes (j+min_up_time,0)...(t,0)
                     to denote the possible shut downs of the unit with respect with the min-up
                     constraint and also the arc connection with (t,1) which represents the
                     case in which the unit remains on for all the rest of the optimisation
                     period.
               
                 ii) For all the time instances that are larger/equal than t-min_up_time, 
                     where the unit is not able to switch its status inside this optimisation
                     period due to min_up constraint and thus the unit will remain online
                     for the rest time steps and the only arc connection needed is (t,1) */

            for ( h = kMin + thermalUnit->get_fMinDownTime() - 1; h < f_UCBlock->get_t() ; h ++ ) {
                
              
                if(h > f_UCBlock->get_t()-1-thermalUnit->get_fMinUpTime()){
                    v_nodes[h].v_arcs.resize(1);
                    v_nodes[h].v_arcs[0].h= h;
                    v_nodes[h].v_arcs[0].k= f_UCBlock->get_t()-1;
                    v_nodes[h].v_arcs[0].cost1= thermalUnit->get_fConstTherm() * (f_UCBlock->get_t() - h ) ;
                    v_nodes[h].v_arcs[0].cost2= 0; 
                    v_nodes[h].v_arcs[0].valid= 2; //from on->on

                }
         
                else{

                     k = h + thermalUnit->get_fMinUpTime() - 1;
                    double c_i=0;
                    for( int t = h ; t<k; ++t) c_i += thermalUnit->get_fConstTherm();

                    int i = 0;
                    v_nodes[h].v_arcs.resize(f_UCBlock->get_t()-k/*+1*/);
                    for (; k < f_UCBlock->get_t() ; k ++ )
                    {
                        v_nodes[h].v_arcs[i].h= h;
                        v_nodes[h].v_arcs[i].k= k/*+f_UCBlock->get_t()*/;
                        c_i += thermalUnit->get_fConstTherm();
                        v_nodes[h].v_arcs[i].cost1= c_i;
                        v_nodes[h].v_arcs[i].cost2= 0; //start_up cost
                        v_nodes[h].v_arcs[i].valid= -1;//from on->off
                        i++;
                    }	      
                 
            }//else
         } //for loop

   } // if for case_1



     /************************************Case_2**********************************/

    else if ( thermalUnit->get_fInitUpDownTime() < 0  && -thermalUnit->get_fInitUpDownTime() < thermalUnit->get_fMinDownTime()){
            // case where the unit needs to remain offline

             /** In this case the unit is off in the first time step and thus we are
                in node (0,0) from where we need to construct all arc connections,
                which for this node are:
                - all nodes (k,1), (k+1,1) , ... , (t,1) where k is the first time 
                  instance in which the unit is able to be switched on
                - node (t,0), which is the case that the unit is not switched of at
                  all and remains on throughout the whole optimisation period  */

            int k = thermalUnit->get_fMinDownTime() + thermalUnit->get_fInitUpDownTime();
            hmin = k;
            kMin = hmin + thermalUnit->get_fMinUpTime() - 1 ;

            {
          
                v_nodes[f_UCBlock->get_t()].v_arcs.resize(f_UCBlock->get_t()-k+1);
                int i=0;
                for (; k< f_UCBlock->get_t() ; k ++ )
                {
                    v_nodes[f_UCBlock->get_t()].v_arcs[i].h= /*f_UCBlock->get_t()*/0 ;
                    v_nodes[f_UCBlock->get_t()].v_arcs[i].k= k;
                    v_nodes[f_UCBlock->get_t()].v_arcs[i].cost1= 0 ;
                    v_nodes[f_UCBlock->get_t()].v_arcs[i].cost2= 0; //start_up
                    v_nodes[f_UCBlock->get_t()].v_arcs[i].valid= 1;//from off->on
                    i++;
                    
                }
                    v_nodes[f_UCBlock->get_t()].v_arcs[i].h= /*f_UCBlock->get_t()*/0;
                    v_nodes[f_UCBlock->get_t()].v_arcs[i].k= /*2**/f_UCBlock->get_t()-1;
                    v_nodes[f_UCBlock->get_t()].v_arcs[i].cost1= 0;
                    v_nodes[f_UCBlock->get_t()].v_arcs[i].cost2= 0;
                    v_nodes[f_UCBlock->get_t()].v_arcs[i].valid= -2;//from off->off
               }
            /*After the constructions of all arc connections for the node of first time
              step we need to create all the different arc conncetions for all the rest
              of the nodes, where we have:

              1) The nodes that denote the unit being online, in which case we start to
                 create arc connections for the nodes: (k,1), (k+1,1) , ... , (t,1). Sin-
                 ce the unit could not have been switched on before. For this connections 
                 we need to take into account the following two cases:
              
                 i)  For all the time instances that are smaller than t-min_up_time, where
                     each (j,1) node is connected with all the nodes (j+min_up_time,0)...(t,0)
                     to denote the possible shut downs of the unit with respect with the min_down
                     constraint and also the arc connection with (t,1) which represents the
                     case in which the unit remains on for all the rest of the optimisation
                     period.
               
                 ii) For all the time instances that are larger/equal than t-min_up_time, 
                     where the unit is not able to switch its status inside this optimisation
                     period due to min_up constraint and thus the unit will remain ofline
                     for the rest time steps and the only arc connection needed is (t,1) */

            for (int h = hmin; h < f_UCBlock->get_t()/*-1*/ ; h ++ ) {
                
              
                if(h > f_UCBlock->get_t()-1-thermalUnit->get_fMinUpTime()){
                    v_nodes[h].v_arcs.resize(1);
                    v_nodes[h].v_arcs[0].h= h;
                    v_nodes[h].v_arcs[0].k= f_UCBlock->get_t()-1;
                    v_nodes[h].v_arcs[0].cost1= thermalUnit->get_fConstTherm() * (f_UCBlock->get_t() - h ) ;
                    v_nodes[h].v_arcs[0].cost2= 0; 
                    v_nodes[h].v_arcs[0].valid= 2; //from on->on

                }
         
                else{

                     k = h + thermalUnit->get_fMinUpTime() - 1 ;
                    double c_i=0;
                    for( int t = h ; t<k; ++t) c_i += thermalUnit->get_fConstTherm();

                    int i = 0;
                    v_nodes[h].v_arcs.resize(f_UCBlock->get_t()-k);
                    for (; k < f_UCBlock->get_t() ; k ++ )
                    {
                        v_nodes[h].v_arcs[i].h= h;
                        v_nodes[h].v_arcs[i].k= k;
                        c_i += thermalUnit->get_fConstTherm();
                        v_nodes[h].v_arcs[i].cost1= c_i;
                        v_nodes[h].v_arcs[i].cost2= 0; //start_up cost
                        v_nodes[h].v_arcs[i].valid= -1;//from on->off
                        i++;
                    }	      
            }//else
           
         }

    /*          2) The nodes that denote the unit being offline, in which case we start to
                 create arc connections for the nodes:
                                     (k+minuptime,0), (k+minuptime+1,0) , ... , (t,0). 
                 Since the unit could not have been switched off again at any possible time 
                 before, after its first. 
                 For these connections we need as above to take into account the following 
                 two cases:
              
                 i)  For all the time instances that are smaller than t-min_down_time, where
                     each (j,0) node is connected with all the nodes (j+min_down_time,1)...(t,1)
                     to denote the possible starts of the unit with respect with the min-down
                     constraint and also the arc connection with (t,0) which represents the
                     case in which the unit remains off for all the rest of the optimisation
                     period.
               
                 ii) For all the time instances that are larger/equal than t-min_down_time, 
                     where the unit is not able to switch its status inside this optimisation
                     period due to min_down constraint and thus the unit will remain online
                     for the rest time steps and the only arc connection needed is (t,0) */

            for (int h = hmin + thermalUnit->get_fMinUpTime() - 1; h < f_UCBlock->get_t()-1 ; h ++ ) {
                
                if(h > f_UCBlock->get_t()-1-thermalUnit->get_fMinDownTime()){
                    v_nodes[h+f_UCBlock->get_t()].v_arcs.resize(1);
                    v_nodes[h+f_UCBlock->get_t()].v_arcs[0].h= h;
                    v_nodes[h+f_UCBlock->get_t()].v_arcs[0].k= f_UCBlock->get_t()-1;
                    v_nodes[h+f_UCBlock->get_t()].v_arcs[0].cost1= 0 ;
                    v_nodes[h+f_UCBlock->get_t()].v_arcs[0].cost2= 0;
                    v_nodes[h+f_UCBlock->get_t()].v_arcs[0].valid= -2; //from off->off
                    
                }
                else{
                    k = h + thermalUnit->get_fMinDownTime() - 1 ;
                    v_nodes[h+f_UCBlock->get_t()].v_arcs.resize(f_UCBlock->get_t()-k+1);
                    int i = 0;
                    for (; k < f_UCBlock->get_t() ; k ++ )
                        {
                        v_nodes[h+f_UCBlock->get_t()].v_arcs[i].h= h/*+f_UCBlock->get_t()*/;
                        v_nodes[h+f_UCBlock->get_t()].v_arcs[i].k= k;
                        v_nodes[h+f_UCBlock->get_t()].v_arcs[i].cost1= 0;
                        v_nodes[h+f_UCBlock->get_t()].v_arcs[i].cost2= 0; //start-up cost
                        v_nodes[h+f_UCBlock->get_t()].v_arcs[i].valid= 1; //from off->on
                        i++;
                        }	    
                    v_nodes[h+f_UCBlock->get_t()].v_arcs[i].h= h/*+f_UCBlock->get_t()*/;
                    v_nodes[h+f_UCBlock->get_t()].v_arcs[i].k= f_UCBlock->get_t()-1/*+f_UCBlock->get_t()*/;
                    v_nodes[h+f_UCBlock->get_t()].v_arcs[i].cost1= 0 ;
                    v_nodes[h+f_UCBlock->get_t()].v_arcs[i].cost2= 0;
                    v_nodes[h+f_UCBlock->get_t()].v_arcs[i].valid= -2; //from off->off
                }//else
                

                

            } //for loop h

   
           
}// if for case_2
     /************************************Case_3**********************************/

   else if (thermalUnit->get_fInitUpDownTime() != 0  )
        {   /* case where there are no restrictions for the unit from the min_up/down
               nor from the ramp constraints */

          hmin = 0;

          /* In this case we construct for all different time_steps and unit states
             all possible connections, where:
           
             1) For (h,0) nodes that denote all the timesteps where the unit is offline
                we have the following arc connections:
                i)  For 0<=h<t-mindown_time there are arcs to all (h+mindown,1), ... , (t,1)
                    nodes that denote that the unit can be switched on at any corresponding
                    forward time_step based on the min_dowb constraints plus an arc to the 
                    node (t,0) that denotes the case that the unit remains of until the end 
                    of the optimisation period.

                ii) For h>t-mindown_time there is only one arc to connect the node with (t,0),
                    since due to the min_down constraint the unit is not able to changes its
                    state inside this optimisation period and as a result will remain off
                    until the last time step.

               
             2) For (h,1) nodes that denote all the timesteps where the unit is online
                we have the following arc connections:
                i)  For 0<=h<t-minup_time there are arcs to all (h+minup,0), ... , (t,0)
                    nodes that denote that the unit can be switched of at any corresponding
                    forward time_step based on the min_up constraints plus an arc to the no-
                    de (t,1) that denote the case that the unit remains on until the end of
                    the optimisation period.

                ii) For h>t-minup_time there is only one arc to connect the node with (t,1),
                    since due to the min_up constraint the unit is not able to changes its
                    state inside this optimisation period and as a result will remain on
                    until the last time step.       
           */
         for (int h = 0; h < f_UCBlock->get_t()/*-1*/ ; h ++ ) {
                
                if(h > 0 && thermalUnit->get_fInitUpDownTime() > 0 && h < thermalUnit->get_fMinDownTime())        continue;
                //constructing arcs for all (h,0) nodes
                if(h > f_UCBlock->get_t()-1-thermalUnit->get_fMinDownTime()){
                    v_nodes[h+f_UCBlock->get_t()].v_arcs.resize(1);
                    v_nodes[h+f_UCBlock->get_t()].v_arcs[0].h= h;
                    v_nodes[h+f_UCBlock->get_t()].v_arcs[0].k= f_UCBlock->get_t()-1;
                    v_nodes[h+f_UCBlock->get_t()].v_arcs[0].cost1= 0 ;
                    v_nodes[h+f_UCBlock->get_t()].v_arcs[0].cost2= 0;
                    v_nodes[h+f_UCBlock->get_t()].v_arcs[0].valid= -2; //from off->off
                    
                }
                else{
                    int k=0;
                    if(h != 0) k = h + thermalUnit->get_fMinDownTime() - 1 ;
                   
                    v_nodes[h+f_UCBlock->get_t()].v_arcs.resize(f_UCBlock->get_t()-k+1);
                    int i = 0;
                    for (; k < f_UCBlock->get_t() ; k ++ )
                        {
                        v_nodes[h+f_UCBlock->get_t()].v_arcs[i].h= h;
                        v_nodes[h+f_UCBlock->get_t()].v_arcs[i].k= k;
                        v_nodes[h+f_UCBlock->get_t()].v_arcs[i].cost1= 0;
                        v_nodes[h+f_UCBlock->get_t()].v_arcs[i].cost2= 0; //start-up cost
                        v_nodes[h+f_UCBlock->get_t()].v_arcs[i].valid= 1; //from off->on
                        i++;
                        }	    
                    v_nodes[h+f_UCBlock->get_t()].v_arcs[i].h= h;
                    v_nodes[h+f_UCBlock->get_t()].v_arcs[i].k= f_UCBlock->get_t()-1;
                    v_nodes[h+f_UCBlock->get_t()].v_arcs[i].cost1= 0 ;
                    v_nodes[h+f_UCBlock->get_t()].v_arcs[i].cost2= 0;
                    v_nodes[h+f_UCBlock->get_t()].v_arcs[i].valid= -2; //from off->off
                }//else
                

                   //constructing arcs for all (h,1) nodes
                    if(h > f_UCBlock->get_t()-1-thermalUnit->get_fMinUpTime()){
                    v_nodes[h].v_arcs.resize(1);
                    v_nodes[h].v_arcs[0].h= h;
                    v_nodes[h].v_arcs[0].k= f_UCBlock->get_t()-1;
                    v_nodes[h].v_arcs[0].cost1= thermalUnit->get_fConstTherm() * (f_UCBlock->get_t() - h ) ;
                    v_nodes[h].v_arcs[0].cost2= 0; 
                    v_nodes[h].v_arcs[0].valid= 2; //from on->on

                }
         
                else{
                    int k=0;
                    if(h != 0) k = h + thermalUnit->get_fMinUpTime() - 1;
                    double c_i=0;
                    for( int t = h ; t<k; ++t) c_i += thermalUnit->get_fConstTherm();

                    int i = 0;
                    v_nodes[h].v_arcs.resize(f_UCBlock->get_t()-k);
                    for (; k < f_UCBlock->get_t() ; k ++ )
                    {
                        v_nodes[h].v_arcs[i].h= h;
                        v_nodes[h].v_arcs[i].k= k;
                        c_i += thermalUnit->get_fConstTherm();
                        v_nodes[h].v_arcs[i].cost1= c_i;
                        v_nodes[h].v_arcs[i].cost2= 0; //start_up cost
                        v_nodes[h].v_arcs[i].valid= -1;//from on->off
                        i++;
                    }	      
            }//else
                

            } //for loop h

       } //if case_3
     else   assert(0);




}

/*--------------------------------------------------------------------------*/

int AcadThemalUnitGraphSolver::solve(void){

 AcadThermalUnitBlock * thermalUnit = dynamic_cast <AcadThermalUnitBlock *> (f_Block);
 UCBlock *f_UCBlock = static_cast<UCBlock*>( thermalUnit->get_f_Block() );

//We need to compute dual contributions for each EDP
for(int i = 0 ; i < f_UCBlock->get_t() ; i++){
    if (v_nodes[i].v_arcs.size() > 0 ){
        for(int j = 0 ; j < v_nodes[i].v_arcs.size() ; j++){
            if( v_nodes[i].v_arcs[j].valid == -1 || v_nodes[i].v_arcs[j].valid == 2 ){
                double m = 0;
                for(int k = v_nodes[i].v_arcs[j].h ; k <= v_nodes[i].v_arcs[j].k ; ++k ){
                    m += -thermalUnit->get_mew()[k]*thermalUnit->get_fMaxPower();
                   // cout<<endl<<"EDP(" << v_nodes[i].v_arcs[j].h << ", " << v_nodes[i].v_arcs[j].k << ")";
                } //for k
                v_nodes[i].v_arcs[j].cost2 = m;
            }// if valid
        }//for j
        
    } // if v_arcs.size()

}// for i

//We need to compute all EDPs
int size = f_UCBlock->get_t() - hmin;
v_EDP.resize(size);
for(int i= hmin; i < f_UCBlock->get_t(); i++)
   {
      v_EDP[i-hmin].initialize(i,f_UCBlock);
      v_EDP[i-hmin].ComputeCosts(v_cost, thermalUnit, f_UCBlock);
      if (v_nodes[i].v_arcs.size() > 0 ){
            for(int j = 0 ; j < v_nodes[i].v_arcs.size() ; j++){
               if( v_nodes[i].v_arcs[j].valid == -1 || v_nodes[i].v_arcs[j].valid == 2 ){
                        v_nodes[i].v_arcs[j].cost2 += v_cost[v_nodes[i].v_arcs[j].k];
        } //if
     } //for j
  } // if
} // for i
MinPath();
get_var_solution();
return 0;
}

/*--------------------------------------------------------------------------*/

void AcadThemalUnitGraphSolver::get_var_solution( ){

 AcadThermalUnitBlock * thermalUnit = dynamic_cast <AcadThermalUnitBlock *> (f_Block);
 UCBlock *f_UCBlock = static_cast<UCBlock*>( thermalUnit->get_f_Block() );


 fPowerVector.resize(f_UCBlock->get_t());
 fStatusVector.resize(f_UCBlock->get_t());

    for( int t = 0; t < f_UCBlock->get_t(); ++t )
    {
        fPowerVector[t]= 0;
        fStatusVector[t]= 0;
    }


    int k= v_route[f_UCBlock->get_t()].pred;
    while( k != -1)
    {
        int h = v_route[k].h;
        v_EDP[h].ComputePowerVariables(k, fPowerVector, thermalUnit);

        for( int t = h ;  t <= k ; ++t) fStatusVector[t]= 1;
        k = v_route[k].pred;
    }
    fTotalCost= v_route[f_UCBlock->get_t()].lab;

    for( int t = 0; t < f_UCBlock->get_t(); ++t )
    {
    	thermalUnit->get_U(t)->set_value(fStatusVector[t]);
        thermalUnit->get_P(t)->set_value(fPowerVector[t]);
        
    }


}

/*--------------------------------------------------------------------------*/

void AcadThemalUnitGraphSolver::MinPath(){

 double eps   = 1e-10 ;  
 const double UCINF = 1e+20;  

   /* Initializing all the possible route-connections of the shortest 
      that start from the source node */

   /** Considering the case that the unit has been initially on, where 
       the start-up costs are zero */

 AcadThermalUnitBlock * thermalUnit = dynamic_cast <AcadThermalUnitBlock *> (f_Block);
 UCBlock *f_UCBlock = static_cast<UCBlock*>( thermalUnit->get_f_Block() );

    if (thermalUnit->get_fInitUpDownTime()>0)
    {

    	/** Examining the case that unit had been turned on for a period 
            shorter than min_up. */

        if (thermalUnit->get_fInitUpDownTime()<thermalUnit->get_fMinUpTime())
        {
            const int h = 0;

            /** Initializing all nodes that can be connected with the source */

            for ( int k = kMin  ; k < f_UCBlock->get_t() ; k ++ )
            {
                /* Initializing the structure */
                v_route[k].h= h;
                v_route[k].lab=  v_nodes[h].v_arcs[k-kMin].cost1 + v_nodes[h].v_arcs[k-kMin].cost2;
                v_route[k].pred= -1;
            }

               /* The target node can not be directly connected to the source node */
                v_route[f_UCBlock->get_t()].lab= UCINF;
        }

        /** In case unit was on for a period larger than the min_up time, 
             we need to take into account the ramp constraint as well */

        else //( thermalUnit->get_fInitUpDownTime() >= thermalUnit->get_fMinUpTime() > 0 )
        {     
            if ( thermalUnit->get_fInitPower() < thermalUnit->get_fBoundDown() + eps )
            {
                /* In case the ramp constraint is not violated we create the
                   connection between source and target node */

                v_route[f_UCBlock->get_t()].h    = -1 ;
                v_route[f_UCBlock->get_t()].lab  =  0 ;
                v_route[f_UCBlock->get_t()].pred = -1 ;

                //Setting initially the cost of all other conections to inf
                for ( int k = 0  ; k < f_UCBlock->get_t()  ; k ++ )
                    v_route[k].lab = UCINF ;

      		/* Initialize the possible connections from the case the unit
                   is being turned off at the beginning of the optimisation period */

                int h = thermalUnit->get_fMinDownTime() ;
                for (  ;  h < f_UCBlock->get_t() - thermalUnit->get_fMinUpTime() + 1 ; h ++ )
                /* checking all the possible starts of unit that allow the unit to be
                  turned off again inside the optimisation period */
                {
                    const double currentstartupcost = ComputeStartupCosts(h);          
                    int i=0;
                    for ( int k = h + thermalUnit->get_fMinUpTime() - 1 ; k < f_UCBlock->get_t() ; k ++ )
                    {

                        const double label = currentstartupcost+v_nodes[h].v_arcs[i].cost1 + v_nodes[h].v_arcs[i].cost2;
                        i++;
                        if ( label < v_route[k].lab )
                        {
                            v_route[k].h= h;
                            v_route[k].lab= label;
                            v_route[k].pred= -1;
                        }
                    }
                }
                for ( ;  h < f_UCBlock->get_t() ; h ++ )
                /* checking the possible starts of unit, for which the unit has to remain
                  on until the end of the optimisation period */
                {
                    const int k     = f_UCBlock->get_t() - 1 ;
                    int   i     = 0;
                    const double label = ComputeStartupCosts(h)+v_nodes[h].v_arcs[i].cost1 + v_nodes[h].v_arcs[i].cost2;
                    i++;
                    if ( label < v_route[k].lab )
                    {
                        v_route[k].h= h;
                        v_route[k].lab= label;
                        v_route[k].pred= -1;
                    }
                }

		/** Checking the case that unit remains on in the begin of optimisation period,
                    creating pairs of (0, kMin), (0, kMin + 1), ..., (0, n-1) */

                for ( int k = kMin ; k < f_UCBlock->get_t() ; k ++ )
                {
                    
                    const double label = v_nodes[0].v_arcs[k - kMin].cost1 + v_nodes[0].v_arcs[k - kMin].cost2;
                    if ( label < v_route[k].lab )
                    {
                        v_route[k].h= 0;
                        v_route[k].lab= label ;
                        v_route[k].pred= -1;
                    }
                }

            }
            else //case that the unit can't be immidiently turned off due to ramp constraint
            {
                 /** No direct connection between source and target nodes */

                v_route[f_UCBlock->get_t()].h= -1;
                v_route[f_UCBlock->get_t()].lab= UCINF;
                v_route[f_UCBlock->get_t()].pred= -1;


		/** Initialisation of all possible pairs, starting from the first time the unit
                    can be turned off: (0, Kmin), (0, Kmin + 1), ..., (0, n-1) */

                const int h= 0;
                for ( int k = kMin ; k < f_UCBlock->get_t() ; k ++ )
                {

                    v_route[k].h= h;
                    v_route[k].lab= v_nodes[0].v_arcs[k - kMin].cost1 + v_nodes[0].v_arcs[k - kMin].cost2;
                    v_route[k].pred= -1;
                }    
            }
        }
    }

    else //examining the case the unit was initially turned off
    {

       
        int idxcs = (-thermalUnit->get_fInitUpDownTime() < thermalUnit->get_fMinDownTime() ? thermalUnit->get_fMinDownTime() : -thermalUnit->get_fInitUpDownTime() ) ;
        //calculating the time-steps that the unit was off

        int h = hmin ;

 	/* Initialize all pairs that start from the first time step unit can be turned on
           and are connected with source node */

        {
            const double currentstartupcost = ComputeStartupCosts(idxcs) ;   


            for ( int k= kMin ; k < f_UCBlock->get_t() ; k ++ )
            {

                const double label = currentstartupcost +v_nodes[h].v_arcs[k - kMin].cost1 + v_nodes[h].v_arcs[k - kMin].cost2;

                v_route[k].h= h;
                v_route[k].lab= label;
                v_route[k].pred= -1;
            }
        }

        ++h;
        ++idxcs;

        /** Initialization of all other possible pairs */


        for (; h < f_UCBlock->get_t() - thermalUnit->get_fMinUpTime() + 1 ; h ++ , idxcs ++ )
        {
            const double currentstartupcost = ComputeStartupCosts(idxcs) ;
            int i = 0;
            for ( int k = h+thermalUnit->get_fMinUpTime()-1 ; k< f_UCBlock->get_t()-1 ; k ++ )
            {
                const double label = currentstartupcost +v_nodes[h].v_arcs[i].cost1 + v_nodes[h].v_arcs[i].cost2;

                i++;
                if(label<v_route[k].lab)
                {
                    v_route[k].h= h;
                    v_route[k].lab= label;
                    v_route[k].pred= -1;
                }
            }
	    // Nodes which terminate after the end of the interval of definition ...
            {
                const int k     = f_UCBlock->get_t()-1;
                const double label = ComputeStartupCosts(idxcs) +v_nodes[h].v_arcs[i].cost1 + v_nodes[h].v_arcs[i].cost2;

                if ( label<v_route[k].lab ) 
                {
                    v_route[k].h= h;
                    v_route[k].lab= label;
                    v_route[k].pred= -1;
                }   
            }
        }
	// Initializing the pair (s, d)
        v_route[f_UCBlock->get_t()].lab  = 0 ;
        v_route[f_UCBlock->get_t()].pred = -1;
    }


    /** Taking into consideration all the other arcs, by considering the following
        three iteration circles: 

         1) for (k = Kmin; k <f_UCBlock->get_t() - get_fMinDownTime() - thermalUnit->get_fMinUpTime(); k ++)

         2) for (; k <f_UCBlock->get_t() - get_fMinDownTime() - 1; k ++)

         3) for (; k <f_UCBlock->get_t(); k ++)

         In case 2) the first of two cycles of r currently present is always skipped.
         The division of the cycle eliminates so thermalUnit->get_fMinUpTime()-1 test operations.
         In case 3) the first two cycles could be removed (thus leaving only
         If the bow ((h, k), d)) thus saving '2 * (get_fMinDownTime() + 1) test operations.
         */


    for( int k= kMin ; k < f_UCBlock->get_t() ; ++k )
    {
        int r = k + thermalUnit->get_fMinDownTime() + 1 ;
        for ( ; r < f_UCBlock->get_t() - thermalUnit->get_fMinUpTime() + 1 ; ++r )
        {
            const double costbeforenode = v_route[k].lab+ComputeStartupCosts( r - 1 - k ) ;
            int i=0;
            for ( int q = r + thermalUnit->get_fMinUpTime() - 1 ; q < f_UCBlock->get_t() ; ++ q )
            {
                const double label = costbeforenode +v_nodes[r].v_arcs[i].cost1 + v_nodes[r].v_arcs[i].cost2;
                i++;
                if ( v_route[q].lab > label )
                {
                    v_route[q].h= r;
                    v_route[q].lab= label;
                    v_route[q].pred= k;
                } // if cost
            } // for q
        } // for r

        for(; r < f_UCBlock->get_t() ; ++r )
        {
            const int q     = f_UCBlock->get_t()-1;
            const double label = v_route[k].lab+ComputeStartupCosts(r-1-k) +v_nodes[r].v_arcs[0].cost1 + v_nodes[r].v_arcs[0].cost2;
            if ( v_route[q].lab > label )
            {
                v_route[q].h= r;
                v_route[q].lab= label;
                v_route[q].pred= k;
            }
        } // for r

        //arco ((h,k),d)

        {
            const double label= v_route[k].lab;
            if ( label < v_route[f_UCBlock->get_t()].lab)
            {
                v_route[f_UCBlock->get_t()].lab= label;
                v_route[f_UCBlock->get_t()].pred= k;
            }
        }
    }// for k

}

/*--------------------------------------------------------------------------*/

double AcadThemalUnitGraphSolver::ComputeStartupCosts( int t)
{
 AcadThermalUnitBlock * thermalUnit = dynamic_cast <AcadThermalUnitBlock *> (f_Block);
    /** se l'unita' e' spenta da piu' di fMaxStartupLevel istanti allora
    i costi di start-up sono costanti */

    if ( t > fMaxStartupLevel) t = fMaxStartupLevel;

    /** altrimenti restituisci il costo corrispondente al numero di ore
    di spegnimento, perche' il primo valore di t nella tabella
    coincide con get_fMinDownTime(), che e' proprio il minimo numero di
    istanti per cui l'unita' puo' essere spenta */

    t -= thermalUnit->get_fMinDownTime();
    return (startup_costs[t]); 
}   




/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*--------------- End File AcadThemalUnitGraphSolver.cpp -------------------*/
/*--------------------------------------------------------------------------*/
