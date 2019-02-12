/*--------------------------------------------------------------------------*/
/*--------------------- File AcadThemalUnitWRSolver.cpp --------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation of the AcadThemalUnitWRSolver class.
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

#include "AcadThemalUnitWRSolver.h"
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

void AcadThemalUnitWRSolver::set_Block( Block *block ){


  if( f_Block )                         // was attached to some other Block
  f_Block->unregister_Solver( this );  // no more so
 f_Block = block;                      // this is the new block now
 f_Block->register_Solver( this );     // register to it


 AcadThermalUnitBlock * thermalUnit = dynamic_cast <AcadThermalUnitBlock *> (f_Block);

//constructing the vector with the start-up costs

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

void AcadThemalUnitWRSolver::build_graph( void )
{

 double eps   = 1e-10 ;  
 AcadThermalUnitBlock * thermalUnit = dynamic_cast <AcadThermalUnitBlock *> (f_Block);
 UCBlock *f_UCBlock = static_cast<UCBlock*>( thermalUnit->get_f_Block() );

    //setting the proper size of the different vectors used
    v_cost.resize(f_UCBlock->get_t());
    v_route.resize(f_UCBlock->get_t()+1);
    v_nodes.resize(f_UCBlock->get_t() * f_UCBlock->get_t());

    for( int i = 0 ; i<v_nodes.size(); ++i) v_nodes[i].valid= 0; 
    // seting initially for all the nodes the valid field equal to zero


    /* In order to construct the graph properly we need to take into 
        account the initial status of the unit */

    if ( thermalUnit->get_fInitUpDownTime() > 0) //checking if it was initially turned on
    {

    /* For the case that the unit was initially turned on we have the
       following two cases:

       1) If the unit had been on for a period shorter than
          min_up time, then it should remain on for a total period
          atleast eqaul to min_up time

       2) Otherwise, if the unit had been on for a period at least
          equal to the min_up time, it can be turned off (with respect 
          to min-diwb ramp constraint)
     */

        if ( thermalUnit->get_fInitUpDownTime() < thermalUnit->get_fMinUpTime() ) // case of being on less than min_up
        {

            /** We should build arc (0, k), taking into account that
                the unit remains on for a period of at least thermalUnit->get_fMinUpTime() */

            int h= hmin = 0;
            int k= thermalUnit->get_fMinUpTime()-thermalUnit->get_fInitUpDownTime()-1;

            //caclulating the kMin taking into consideration min_up time & ramp constraints

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

            if ( kMin < k ) kMin = k;
            else            k    = kMin ;

            
            /** Construction of all possible nodes connected directly
                with the source (0, k), ..., (0, n-1).
                Then initialize the fields that make up the node structure */

            {

                double c_i=0;
                for( int t = h ; t<k; ++t) c_i += thermalUnit->get_fConstTherm();

                for (; k< f_UCBlock->get_t() ; k ++ )
                {
                    int i= position(h, k); //finding the position of the node in the vector

                    /* initializing the node*/

                    v_nodes[i].h= h;
                    v_nodes[i].k= k;
                    c_i += thermalUnit->get_fConstTherm();
                    v_nodes[i].cost1= c_i ;
                    v_nodes[i].cost2= 0;
                    v_nodes[i].valid= 1;
                }

            }
            /* Construction of the other nodes */

            for ( h = thermalUnit->get_fMinUpTime()-thermalUnit->get_fInitUpDownTime()+thermalUnit->get_fMinDownTime() ;
                  h < f_UCBlock->get_t() - thermalUnit->get_fMinUpTime() + 1 ; h ++ )
            	// constructing all the pairs that start after the unit has shut down from the initial time up
            	// until the last possible start that has a distance from sink larged that min up time
            {
                k = h + thermalUnit->get_fMinUpTime() - 1 ;

                double c_i=0;
                for( int t = h; t < k ; ++ t ) c_i += thermalUnit->get_fConstTherm();

                for (; k < f_UCBlock->get_t() ; k ++ )
                {
                    int i= position(h, k);
                    v_nodes[i].h= h;
                    v_nodes[i].k= k;
                    c_i += thermalUnit->get_fConstTherm();
                    v_nodes[i].cost1= c_i ;
                    v_nodes[i].cost2= 0;
                    v_nodes[i].valid= 1;
                }	      
            }

            //nodes which continue in the next planning period
            for (; h < f_UCBlock->get_t() ; h ++ )
            {
                k = f_UCBlock->get_t() - 1;
                double c_i=0;
                for( int t = h; t <= k ; ++ t ) c_i += thermalUnit->get_fConstTherm();

                int i= position(h, k);
                v_nodes[i].h= h;
                v_nodes[i].k= k;
                v_nodes[i].cost1= c_i ;
                v_nodes[i].cost2= 0;
                v_nodes[i].valid= 1;
            }       
        }

        else  //case where thermalUnit->get_fInitUpDownTime() >= thermalUnit->get_fMinUpTime()
        {   

       	 //Calculation of kMin, taking into account the ramp-down constraint

            kMin= 0;

            if ( thermalUnit->get_fInitPower() >= thermalUnit->get_fBoundDown() + eps )
            {
                double tmp = thermalUnit->get_fInitPower() ;
                tmp -= thermalUnit->get_fMaxRampDown();
                kMin ++ ;
                while ( tmp >= thermalUnit->get_fBoundDown() + eps )
                {
                    tmp -= thermalUnit->get_fMaxRampDown();
                    kMin ++ ;
                }
                kMin -- ;
            }
            if ( kMin >= f_UCBlock->get_t() ) kMin = f_UCBlock->get_t() - 1 ;


            /* constructing of all possible different nodes connected directly 
               with source (0,kMin),(0,kMin+1),...,(0,n-1) */

            int h = hmin = 0;
            int k = kMin ;
            double c_i = 0 ;
            for( int t = 0 ; t < k ; ++ t ) c_i += thermalUnit->get_fConstTherm();
            for (; k < f_UCBlock->get_t() ; k ++ )
            {
                int i= position(h, k);
                v_nodes[i].h= h;
                v_nodes[i].k= k;
                c_i += thermalUnit->get_fConstTherm();
                v_nodes[i].cost1= c_i ;
                v_nodes[i].cost2= 0;
                v_nodes[i].valid= 1;
            }   

            /* construction of all the (h,k) possible nodes */

            for ( h = thermalUnit->get_fMinDownTime() ; h < f_UCBlock->get_t() - thermalUnit->get_fMinUpTime() + 1 ; h ++ )
            {
                int k = h + thermalUnit->get_fMinUpTime() - 1 ;

                double c_i = 0 ;
                for ( int t = h ; t < k ; t ++ ) c_i += thermalUnit->get_fConstTherm() ;

                for (; k < f_UCBlock->get_t() ; k ++ )
                {
                    int i= position(h, k);
                    v_nodes[i].h= h;
                    v_nodes[i].k= k;
                    c_i += thermalUnit->get_fConstTherm();
                    v_nodes[i].cost1= c_i ;
                    v_nodes[i].cost2= 0;
                    v_nodes[i].valid= 1;
                }
            }
            //nodes which continue in the next planning period
            for (; h < f_UCBlock->get_t() ; h ++ )
            {
                int k = f_UCBlock->get_t() - 1 ;
                double c_i = 0 ;
                for ( int t = h ; t <= k ; t ++ ) c_i += thermalUnit->get_fConstTherm() ;

                int i= position(h, k);
                v_nodes[i].h= h;
                v_nodes[i].k= k;
                v_nodes[i].cost1= c_i ;
                v_nodes[i].cost2= 0;
                v_nodes[i].valid= 1;
            }
        } // finish of case thermalUnit->get_fInitUpDownTime() >= thermalUnit->get_fMinUpTime()
    } // finish of case InitUpDownTime >0

    /** If the unit was switched off, you can have two possibilities:

        1) the unit to be switched off for a period less than get_fMinDownTime(),
           in which case the unit needs to remain switched off until it reaches
           get_fMinDownTime()

        2) the unit was off for at least get_fMinDownTime() time steps, in which
           case it can be turned off (with respect to min-up ramp constraint).
    */

    else
    {
        if ( thermalUnit->get_fInitUpDownTime() < 0 )   {

            // compute the minimum time hmin such that the unit can be switchend on 



            int h = (-thermalUnit->get_fInitUpDownTime()< thermalUnit->get_fMinDownTime() ? thermalUnit->get_fMinDownTime() + thermalUnit->get_fInitUpDownTime() :  0 );

            hmin = h ;
            kMin = hmin + thermalUnit->get_fMinUpTime() - 1 ;
            if ( kMin >= f_UCBlock->get_t() ) kMin = f_UCBlock->get_t() - 1 ;

            // construction of all nodes
            for (; h < f_UCBlock->get_t() - thermalUnit->get_fMinUpTime() + 1 ; h ++ )
            {
                int k = h + thermalUnit->get_fMinUpTime() - 1 ;

                double c_i = 0 ;
                for ( int t = h ; t < k ; t ++ ) c_i += thermalUnit->get_fConstTherm();

                for ( ; k < f_UCBlock->get_t() ; k ++ )
                {
                    int i= position(h, k);
                    v_nodes[i].h= h;
                    v_nodes[i].k= k;
                    c_i += thermalUnit->get_fConstTherm();
                    v_nodes[i].cost1= c_i ;
                    v_nodes[i].cost2= 0;
                    v_nodes[i].valid= 1;
                }
            }
            // nodi che terminano nel prossimo periodo di pianiificazione ...
            for (; h < f_UCBlock->get_t() ; h ++ )
            {
                int k = f_UCBlock->get_t() - 1 ;
                double c_i = 0 ;
                for ( int t = h ; t <= k ; t ++ ) c_i += thermalUnit->get_fConstTherm();

                int i= position(h, k);
                v_nodes[i].h= h;
                v_nodes[i].k= k;
                v_nodes[i].cost1= c_i ;
                v_nodes[i].cost2= 0;
                v_nodes[i].valid= 1;
            }
        }

        /** if thermalUnit->get_fInitUpDownTime() equal to zero we terminate the program,
            perche' e' impossibile visto che l'unita', prima dell'istante 
            iniziale da noi considerato poteva essere accesa o spenta  */

        else   //
            assert(0);
    }

}

/*--------------------------------------------------------------------------*/

int AcadThemalUnitWRSolver::solve(void){

 AcadThermalUnitBlock * thermalUnit = dynamic_cast <AcadThermalUnitBlock *> (f_Block);
 UCBlock *f_UCBlock = static_cast<UCBlock*>( thermalUnit->get_f_Block() );

 //compute contributions
 int nodoSize= f_UCBlock->get_t()*f_UCBlock->get_t();

    for( int i = 0; i < nodoSize; ++i){
        if( v_nodes[i].valid ) 
        {
            double  m = 0;
            for ( int j = v_nodes[i].h ; j <= v_nodes[i].k ; ++j )
                m += -thermalUnit->get_mew()[j]*thermalUnit->get_fMaxPower();
                  //  cout<<endl<<"EDP(" << v_nodes[i].h << ", " << v_nodes[i].k << ")";
                v_nodes[i].cost2= m;
         }
}


int size = f_UCBlock->get_t() - hmin ;
v_EDP.resize(size);
for(int i= hmin; i < f_UCBlock->get_t(); i++)
   {
    v_EDP[i-hmin].initialize(i,f_UCBlock);
        v_EDP[i-hmin].ComputeCosts(v_cost, thermalUnit, f_UCBlock);
        for( int j = i; j < f_UCBlock->get_t(); j++)
        {
            int k= position(i,j);
            if(v_nodes[k].valid){
                v_nodes[k].cost2 += v_cost[j];
             }
        }
   }

MinPath();
get_var_solution();
return 0;
}

/*--------------------------------------------------------------------------*/

void AcadThemalUnitWRSolver::get_var_solution( ){

 AcadThermalUnitBlock * thermalUnit = dynamic_cast <AcadThermalUnitBlock *> (f_Block);
 UCBlock *f_UCBlock = static_cast<UCBlock*>( thermalUnit->get_f_Block() );

 fPowerVector.resize(f_UCBlock->get_t());
 fStatusVector.resize(f_UCBlock->get_t());
    for( int t = 0; t < f_UCBlock->get_t(); ++t )
    {
        fPowerVector[t]= 0;
        fStatusVector[t]= 0;
    }

    /* percorso all'indietro del cammino ottimo  */
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

void AcadThemalUnitWRSolver::MinPath(){

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
                const int i = position(h, k);   //find the proper position of node

                /* Initializing the structure */
                v_route[k].h= h;
                v_route[k].lab= v_nodes[i].cost1 + v_nodes[i].cost2;
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
                    for ( int k = h + thermalUnit->get_fMinUpTime() - 1 ; k < f_UCBlock->get_t() ; k ++ )
                    {
                        const int   i     = position(h, k);
                        const double label = currentstartupcost+v_nodes[i].cost1+v_nodes[i].cost2;
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
                    const int   i     = position(h, k);
                    const double label = ComputeStartupCosts(h)+v_nodes[i].cost1+v_nodes[i].cost2;
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
                    const int   i     = position(0,k);
                    const double label = v_nodes[i].cost1 + v_nodes[i].cost2; 
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
                    const int   i     = position(h,k);
                    v_route[k].h= h;
                    v_route[k].lab= v_nodes[i].cost1 + v_nodes[i].cost2; ;
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

            // kMin == min ( h+thermalUnit->get_fMinUpTime()-1 , f_UCBlock->get_t() - 1) ;
            for ( int k= kMin ; k < f_UCBlock->get_t() ; k ++ )
            {
                const int   i     = position(h, k);
                const double label = currentstartupcost +v_nodes[i].cost1+v_nodes[i].cost2;
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
            for ( int k = h+thermalUnit->get_fMinUpTime()-1 ; k< f_UCBlock->get_t()-1 ; k ++ )
            {
                const int   i     = position(h, k);
                const double label = currentstartupcost +v_nodes[i].cost1+v_nodes[i].cost2;
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
                const int   i     = position(h, k);
                const double label = ComputeStartupCosts(idxcs) +v_nodes[i].cost1+v_nodes[i].cost2;
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


    /** Taking into consideration all the other nodes, by considering the following
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
            for ( int q = r + thermalUnit->get_fMinUpTime() - 1 ; q < f_UCBlock->get_t() ; ++ q )
            {
                const int   i     = position(r,q);
                const double label = costbeforenode + v_nodes[i].cost1 + v_nodes[i].cost2;
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
            const int   i     = position(r,q);
            const double label = v_route[k].lab+ComputeStartupCosts(r-1-k) + v_nodes[i].cost1+v_nodes[i].cost2;
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

double AcadThemalUnitWRSolver::ComputeStartupCosts( int t)
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
}   // fine metodo ComputeStartupCosts 

/*--------------------------------------------------------------------------*/

int AcadThemalUnitWRSolver::position(int h, int k){
 AcadThermalUnitBlock * thermalUnit = dynamic_cast <AcadThermalUnitBlock *> (f_Block);
 UCBlock *f_UCBlock = static_cast<UCBlock*>( thermalUnit->get_f_Block() );
    return (h*f_UCBlock->get_t()+k);
  }



/*--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------*/
/*------------------ End File AcadThemalUnitWRSolver.cpp -------------------*/
/*--------------------------------------------------------------------------*/
