/*--------------------------------------------------------------------------*/
/*---------------------------- File EDPSolver.cpp --------------------------*/
/*--------------------------------------------------------------------------*/
/** @file
 * Implementation file for the Economic Dispatch Solver of Academic ThermalUnit 
 * with Ramp Constraints
 *
 *
 * \version 0.10
 *
 * \date 03 - 10 - 2016
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

#include "EDPSolver.h"

/*--------------------------------------------------------------------------*/
/*------------------------- NAMESPACE AND USING ----------------------------*/
/*--------------------------------------------------------------------------*/

using namespace SMSpp_di_unipi_it;

/*--------------------------------------------------------------------------*/
/*--------------------------------- METHODS --------------------------------*/
/*--------------------------------------------------------------------------*/

void EDPSolver::initialize ( int k, UCBlock *f_UCBlock ) 
{

h = k;
kMax = f_UCBlock->get_t() ;

  int coeffsize = (f_UCBlock->get_t())*(f_UCBlock->get_t())+h*h-2*h*(f_UCBlock->get_t());
//cout<<endl<<"Just checking...  ( " << h << " ; " << coeffsize << " ; " << Coeff.size() << " )";
//cout<<endl<<"//";
  if ( coeffsize != Coeff.size() )
   {
      Coeff.resize( coeffsize );
      
      int msize = coeffsize + (f_UCBlock->get_t()) - h ;
      m.resize( msize );

      v.resize((f_UCBlock->get_t()) );
      Position.resize((f_UCBlock->get_t()) );
      unconstrPowerOpt.resize((f_UCBlock->get_t()) );
      constrPowerOpt.resize((f_UCBlock->get_t()) );
  }
  eps = 1e-10 ;


}

/*--------------------------------------------------------------------------*/

void EDPSolver::ComputeCosts(std::vector< double > &costVector, AcadThermalUnitBlock *unit_block, UCBlock *f_UCBlock)
{

 int mcnt;
  int coeffcnt;

  int k= h;

  /* Alias per dati */

  //const tRowCost &QuadTherm = *fpQuadTherm ;

  /* inizializzazione dei costi */

  Coeff[0].alfa= unit_block->get_fQuadTherm();
  Coeff[0].beta= unit_block->get_fLinearTherm()- unit_block->get_lambda()[k];
  Coeff[0].gamma= 0;
  //cout<<endl<< "coeff[0].gamma = " << Coeff[0].gamma;
  coeffcnt= 1;   //indice della prossima posizione libera nell'array Coeff
  v[k]= 0;   //poiche' per k= h il numero di tratti e' pari a 1


  // Because 'for k = h the number of strokes and' equal to 1

     /* Initialize the array containing the m of the various extremes
        traits. Initially it contains the two extremes of the individual
        stretch. On power-up the power can not 'be
        above a certain value / bar {s} _k */

  /* inizializzazione dell'array contenente gli estremi m dei vari
     tratti. Inizialmente contiene i due estremi del singolo
     tratto. Al momento dell'accensione la potenza non potra' essere
     superiore ad un certo valore /bar{l}_k  */

  if((h==0)&&(unit_block->get_fInitUpDownTime()>0)){
    m[0]= max(unit_block->get_fMinPower(), unit_block->get_fInitPower() - unit_block->get_fMaxRampDown());
    m[1]= min(unit_block->get_fMaxPower(), unit_block->get_fInitPower() + unit_block->get_fMaxRampUp());
  }
  else{
    m[0]= unit_block->get_fMinPower();
    m[1]= unit_block->get_fBoundOn();   // /bar{l}_k;
  }
  mcnt= 2;   //indice della prossima posizione libera in m

  /* inizializzazione dei campi di Position contenenti gli indici
     iniziali dei vari tratti */

  Position[k].begm= 0;
  Position[k].begt= 0;

  /* inizializzazione delle potenze ottime non vincolate
     unconstrPowerOpt[k]. Inizializzazione quindi di
     unconstrPowerOpt[h] */

  /* unconstrPowerOpt[k] non e' vincolato ad essere <=
     get_fBoundDown()[k] */


  {
    /* calcolo della potenza minima del tratto iniziale */
    
    double tmp_pStar_hk= - Coeff[0].beta/(2*Coeff[0].alfa);
    if ( tmp_pStar_hk < m[0])
      unconstrPowerOpt[k]= m[0]; //unconstrPowerOpt[k] non e' vincolato ad essere <= get_fBoundDown()[k]
    else
      if (tmp_pStar_hk > m[1])
	unconstrPowerOpt[k]= m[1];
      else
	unconstrPowerOpt[k]= tmp_pStar_hk;



  //cout<< endl << "unconstrPowerOpt[k]_0 = " << unconstrPowerOpt[k] << " && m[0] = " << m[0] << " && m[1] =" << m[1] << " && tmp_pStar_hk = " << tmp_pStar_hk;
  }

  /* inizializzazione delle potenze ottime vincolate constrPowerOpt[k]
     che verrannno calcolate ad ogni iterazione */

  /* constrPowerOpt[k] e' vincolato ad essere <= (*unit_block->get_fMinPower())[k] */



  if (( k < (f_UCBlock->get_t())-1 ) && ( unconstrPowerOpt[k] > unit_block->get_fMinPower()))
    constrPowerOpt[k]= unit_block->get_fMinPower();
  else
    constrPowerOpt[k]= unconstrPowerOpt[k];
  
  costVector[k]= Coeff[0].alfa * constrPowerOpt[k]*constrPowerOpt[k]+ Coeff[0].beta*constrPowerOpt[k];

  /* ciclo piu' esterno */

  for(k= h+1; k < kMax ; ++k){

    /* costruzione dei tratti:/bar{m}_0 e' il primo estremo del primo
       tratto diz_hk. Tale elemento verra' caricato nel vettore degli
       m */

    Position[k].begm= mcnt;
    Position[k].begt= coeffcnt;

    /* calcolo del valore di/bar{m}_0: primo estremo del tratto
       considerato*/

    if(unit_block->get_fMinPower() >m[Position[k-1].begm]-unit_block->get_fMaxRampDown())
      m[mcnt]= unit_block->get_fMinPower();
    else
      m[mcnt]= m[Position[k-1].begm]-unit_block->get_fMaxRampDown();

    double p_bar= m[mcnt];   //inizializzo p_bar=/bar{m}_0 di z_hk(la funzione successiva)
    int    v_bar= 0;   //inizializzazione. Dopo il passo 3) ci dira' quanto vale v[k]


    /* calcolo di q: indice della posizione dell'intervallo m[i] del
       tratto a cui appartiene p^*(bar{p}) */
    
    double pStarDIp_bar;
        //cout<< endl << "unconstrPowerOpt[k-1]_2 = " << unconstrPowerOpt[k-1];

    if ( p_bar < unconstrPowerOpt[k-1]) 
      {
	pStarDIp_bar= p_bar + unit_block->get_fMaxRampDown();
	if ( pStarDIp_bar > unconstrPowerOpt[k-1] )
	  pStarDIp_bar = unconstrPowerOpt[k-1] ;
      }
    else 
      {
	pStarDIp_bar= p_bar -unit_block->get_fMaxRampUp();
	if ( pStarDIp_bar < unconstrPowerOpt[k-1] )
	  pStarDIp_bar = unconstrPowerOpt[k-1] ;	
      }
    
    int qm= Position[k-1].begm;
      
      while ((pStarDIp_bar>= m[qm+1])&&(qm < Position[k].begm-2)){
	++qm;
      
      /* if ( pStarDIp_bar > m[qm+1] ) ERRORE */
      
    }   //fine ciclo while per il calcolo di q di appartenenza

    int q = qm - Position[k-1].begm + Position[k-1].begt;

    /* calcola l'ultimo estremo del tratto, u_bar */

    double u_bar ; 
    if(unit_block->get_fMaxPower()< m[mcnt-1]+unit_block->get_fMaxRampUp())
      u_bar= unit_block->get_fMaxPower();
    else
      u_bar= m[mcnt-1] + unit_block->get_fMaxRampUp();
    ++mcnt;


    int firstTime= 1;
    /* caso 1) */
    //cout<< endl << "unconstrPowerOpt[k-1]_3 = " << unconstrPowerOpt[k-1];
    while((unconstrPowerOpt[k-1]) > (p_bar+unit_block->get_fMaxRampDown()) + eps ){

      /* setta i campi di Coeff per calcolare la
	 /bar{z}^{/bar{v}}(p) */

      Coeff[coeffcnt].alfa=  unit_block->get_fQuadTherm() + Coeff[q].alfa;
      Coeff[coeffcnt].beta= unit_block->get_fLinearTherm() - unit_block->get_lambda()[k] + Coeff[q].beta+ 2*unit_block->get_fMaxRampDown()* Coeff[q].alfa;  //(*fpLinearTherm)[k]- (*fpget_lambda())[k]+ Coeff[q].beta+ 2*(*fpMaxRampDown)[k-1]* Coeff[q].alfa;
      Coeff[coeffcnt].gamma= Coeff[q].gamma+ Coeff[q].alfa*(unit_block->get_fMaxRampDown())*(unit_block->get_fMaxRampDown())+ Coeff[q].beta*unit_block->get_fMaxRampDown();

      /* calcolo del massimo valore di p_bar tale che p*_k(p_bar)
	 rimanga nel q-esimo intervallo,unconstrPowerOpt sia fuori
	 dell'intervallo ammissibile e p_bar resti ammissibile */

      if( ( m[qm+1]-unit_block->get_fMaxRampDown() ) < (unconstrPowerOpt[k-1]- unit_block->get_fMaxRampDown()) - eps ) {
	p_bar= m[qm+1]- unit_block->get_fMaxRampDown();
	++q; ++qm ;
      }   //fine if
      else
	p_bar=unconstrPowerOpt[k-1]- unit_block->get_fMaxRampDown();
      if (p_bar> u_bar)
	p_bar= u_bar;
      ++v_bar;
      m[mcnt++]= p_bar;

      /* calcolo diunconstrPowerOpt,ottimo non vincolato di z_hk */

      if((firstTime)&&(2*Coeff[coeffcnt].alfa*p_bar+Coeff[coeffcnt].beta >0 )){
	unconstrPowerOpt[k]= - Coeff[coeffcnt].beta/(2*Coeff[coeffcnt].alfa);
	if(unconstrPowerOpt[k] < m[mcnt-2])
	  unconstrPowerOpt[k]= m[mcnt-2];
	firstTime= 0;
      }   //fine if per il calcolo dell'ottimo non vincolato

      ++coeffcnt;
    }   //fine ciclo while del caso 1)


    /* Caso 2): un semplice if */

    if(unconstrPowerOpt[k-1]>= p_bar- unit_block->get_fMaxRampUp()){

      /* setta i campi di Coeff per calcolare la z_barDIv_bar(p) */

      Coeff[coeffcnt].alfa= unit_block->get_fQuadTherm();
      Coeff[coeffcnt].beta= unit_block->get_fLinearTherm() - unit_block->get_lambda()[k]; //((*fpLinearTherm)[k])- ((*fpget_lambda())[k]);
    /*   cout<<endl<< "For coeff[coeffcnt].gamma_2 : Coeff[coeffcnt].gamma_old =" << Coeff[coeffcnt].gamma << " Coeff[q].alfa = " << Coeff[q].alfa
      << " && unconstrPowerOpt[k-1] = " << unconstrPowerOpt[k-1] << " && Coeff[q].beta = " << Coeff[q].beta;*/
      Coeff[coeffcnt].gamma= Coeff[q].alfa*unconstrPowerOpt[k-1]*unconstrPowerOpt[k-1]+ Coeff[q].beta*unconstrPowerOpt[k-1]+ Coeff[q].gamma;
      //cout<<endl<< "coeff[coeffcnt].gamma_2 = " << Coeff[coeffcnt].gamma << " ";
     /* calcola il massimo valore di p_bar tale cheunconstrPowerOpt rimane
	all'interno dell'intervallo ammissibile e p_bar resti
	ammissibile */

      if((unconstrPowerOpt[k-1]+unit_block->get_fMaxRampUp())< u_bar)
	p_bar=unconstrPowerOpt[k-1]+unit_block->get_fMaxRampUp();
      else
	p_bar= u_bar;
      ++v_bar;
      m[mcnt++]= p_bar;
      if(((2*Coeff[coeffcnt].alfa*p_bar+Coeff[coeffcnt].beta)>0)&& firstTime){
	unconstrPowerOpt[k]= - Coeff[coeffcnt].beta/(2*Coeff[coeffcnt].alfa);
	if(unconstrPowerOpt[k] < m[mcnt-2])
	  unconstrPowerOpt[k]= m[mcnt-2];
	firstTime= 0;
      }
      ++coeffcnt;
    }   //fine if del caso 2)


    /* Caso 3) */

    while(p_bar< u_bar){

      /* setta i campi di Coeff per calcolare la z_barDIv_bar(p) */

      Coeff[coeffcnt].alfa= unit_block->get_fQuadTherm()+ Coeff[q].alfa;
      Coeff[coeffcnt].beta= unit_block->get_fLinearTherm() - unit_block->get_lambda()[k] + Coeff[q].beta- 2*unit_block->get_fMaxRampUp()* Coeff[q].alfa;//- (*fpget_lambda())[k]+ Coeff[q].beta- 2*(*fpMaxRampUp)[k-1]* Coeff[q].alfa;
      Coeff[coeffcnt].gamma= Coeff[q].gamma+ Coeff[q].alfa*(unit_block->get_fMaxRampUp())*(unit_block->get_fMaxRampUp())- Coeff[q].beta*unit_block->get_fMaxRampUp();
        //      cout<<endl<< "coeff[coeffcnt].gamma_3 = " << Coeff[coeffcnt].gamma;
      /* calcolo del massimo valore di p_bar tale che p*_k(p_bar)
	 rimanga nel q-esimo intervallo e p_bar resti ammissibile */

      if(m[qm+1]+unit_block->get_fMaxRampUp()< u_bar)
	p_bar= m[qm+1]+unit_block->get_fMaxRampUp();
      else
	p_bar= u_bar;
      ++v_bar;
      m[mcnt++]= p_bar;
      ++q; ++qm ;
      if(((2*Coeff[coeffcnt].alfa*p_bar+Coeff[coeffcnt].beta)>0)&& firstTime){
	unconstrPowerOpt[k]= - Coeff[coeffcnt].beta/(2*Coeff[coeffcnt].alfa);
	if(unconstrPowerOpt[k] < m[mcnt-2])
	  unconstrPowerOpt[k]= m[mcnt-2];
	firstTime= 0;
      }
      ++coeffcnt;
    }   //fine ciclo while del caso 3)


    v[k]= v_bar-1;

    /* se la funzione e' sempre decrescente allora ... */

    if ( firstTime )
      unconstrPowerOpt[k] = u_bar ;

    /* calcolo di constrPowerOpt[k], ottimo vincolato relativo all'intera funzione */

    if (( k < (f_UCBlock->get_t())-1 ) && ( unconstrPowerOpt[k] > unit_block->get_fMinPower()))
      constrPowerOpt[k]= unit_block->get_fMinPower();
    else
      constrPowerOpt[k]= unconstrPowerOpt[k];
    
    /* calcola costo nodo (h,k) in costVector[k] */

    qm = Position[k].begm;

    while (constrPowerOpt[k]> m[qm+1]){
      ++qm;   
    }   //fine ciclo while per il calcolo di q di appartenenza

    q = qm - Position[k].begm + Position[k].begt;

    costVector[k]= Coeff[q].alfa*constrPowerOpt[k]*constrPowerOpt[k]+Coeff[q].beta*constrPowerOpt[k]+Coeff[q].gamma;

    /*cout<<endl<<"costVector[k] = " << costVector[k] << " && Coeff[q].alfa = " << Coeff[q].alfa
              <<" && constrPowerOpt[k] = " << constrPowerOpt[k] << " && Coeff[q].beta = " << Coeff[q].beta 
              <<" && constrPowerOpt[k] =" << constrPowerOpt[k] << " && Coeff[q].gamma = " << Coeff[q].gamma;*/
  }   //fine ciclo for piu' esterno

}


/*--------------------------------------------------------------------------*/


void EDPSolver::ComputePowerVariables(int k, vector<double> &powerVector, AcadThermalUnitBlock *unit_block){

      int t;


      /* calcolo delle potenze ottime p_h[h], p_h[h+1],...,
	 p_h[k-1]  */
      powerVector[k]= constrPowerOpt[k];
      for(t= k-1; t>= h; --t){

	/* proiezione dell'ottimo non vincolatounconstrPowerOpt[t]
	   sull'intervallo ottenuto considerando la potenza ottima
	   dell'istante successivo, quindi su:
	   [powerVector[t+1]-get_fMaxRampUp()[t],
	   powerVector[t+1]+fRampRampDown[t]] */

	/* se l'ottimo non vincolato si trova a sinstra
	   dell'intervallo, allora la potenza ottima corrispondera'
	   all'estremo sinistro della funzione */

	if (unconstrPowerOpt[t]< powerVector[t+1]-unit_block->get_fMaxRampUp())
	  powerVector[t]= powerVector[t+1]-unit_block->get_fMaxRampUp();

	/* se l'ottimo non vincolato si trova all'interno
	   dell'intervallo, allora la potenza ottima sara' proprio
	   uguale all'ottimo non vincolato */

	else if(unconstrPowerOpt[t]<= powerVector[t+1]+unit_block->get_fMaxRampDown())
	  powerVector[t]=unconstrPowerOpt[t];

	/*  se l'ottimo non vincolato si trova a destra
	    dell'intervallo, allora la potenza ottima corrispondera'
	    all'estremo destro della funzione */

	else
	  powerVector[t]= powerVector[t+1]+unit_block->get_fMaxRampDown();
      }   //fine ciclo for
    }
