#include<math.h>
#include<stdio.h>
#include<stdlib.h>

//****** Parametros do gerador aleatorio ********************//
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)
#define PI acos(-1.0)
#define n_n  1000000
#define aloca 1000
#define condi 1

//******   Parametros do sistema ********************//
#define NN 1000         // numero of IF
#define N 2             // numero de equacoes
#define h 0.01          // passo de integracao
#define transient 1000     // ms
#define NMD 20000 //numero maximo de disparo de um neuronio: 10^4
#define conexMAX  10.0*pow(10,3)   // each neuron max conections number
#define P_exc_ini 0.8  //porcentagem de conexões excitatórias

//*******gerador de numeros aleatorios******//
float ran1(long *idum);
float gasdev(long *idum);

#define NR_END 1
#define FREE_ARG char*

//*******ponteiros para alocar memoria*****///
void nrerror(char error_text[]);
int *vector(long nl,long nh);
int **imatrix(long nrl, long nrh, long ncl, long nch);
void free_vector(int *v, long nl, long nh);
void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
double *dvector(long nl,long nh);
void free_dvector(double *v, long nl, long nh);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);

void derivs(double y[],double df[],double *aa,double gexc_exc,
double gini_ini ,double gexc_auto,
double gini_auto,double gexc_inh,
double ginh_exc,double *Gexc_exc,
double *Gini_ini,
double *Gexc_inh,double *Ginh_exc,
double *Gexc_auto,double *Gini_auto);

FILE *o,*p,*p1,*p2,*q;  

int main(void)
{
double *df,*y,*a,*b,*c,*x,*aa,*tpeak,*cont_Iext,step,
marca_tempo[1001][aloca +1],delay_exc,delay_ini,delay,CV_ant,*Gexc_exc,
*Gini_ini,tauex,tauin,*Gexc_auto,*Gini_auto,*Gexc_inh,*Ginh_exc,gexc_inh,ginh_exc;
int i,j,cont,t,aux,contISI,n,conta_disparo[1001],atualiza_disparo[1001],auxi,rep,kkk_cont;
int **listaexc,*conextotalex,**listaini,*conextotalin,cont_hyper;
double tempo,b_IF,gexc_exc,gini_ini,V_reset,gexc_auto,gini_auto;
double kinicial,kfinal,R,real,compl,Rmedio,**nk,*phi;
int contR,*k,*kmax; 
double ISI,ISI2,CV,desviopadrao,KM,F,g,FM,RM;
int numerodeconexoes,sum,auxx,jj;
double C,gL,VT,DeltaT,tauw,Vr_exc,Vr_ini,EL,taum,IsynM,IsynMI,CVM,IM2,IMI;

long idum=-1234567897; //semente para gerar números aleatórios

int variable;
char output_filename[100];   

int N_aut_exc, N_aut_inh, N_exc_exc,N_inh_inh,N_exc_inh, N_inh_exc;
int Cont_aut_exc, Cont_aut_inh, Cont_exc_exc, Cont_inh_inh, Cont_exc_inh, 
Cont_inh_exc;

scanf("%d", &variable);

auxi=1;

if(variable==1) auxi=0;

sprintf(output_filename,"POK_CV_F_%d.dat",variable); //reo=256.2
o=fopen(output_filename,"wt");

//////////////Alocando memória //////////////////////////////////////
cont_Iext=dvector(1,NN+1);  aa=dvector(1,NN+1); y=dvector(1,N*NN+1);
df=dvector(1,N*NN+1); x=dvector(1,N*NN+1); a=dvector(1,N*NN+1);
b=dvector(1,N*NN+1); c=dvector(1,N*NN+1); k=vector(1,NN+1);
kmax=vector(1,NN+1); phi=dvector(1,NN+1); nk=dmatrix(1,NMD+2,1,NN+2);
tpeak=dvector(1,NN+1); Gexc_exc=dvector(1,NN+1); Gexc_auto=dvector(1,NN+1);
Gini_ini=dvector(1,NN+1); Gini_auto=dvector(1,NN+1);Ginh_exc=dvector(1,NN+1);
Gexc_inh=dvector(1,NN+1); conextotalex=vector(1,NN+2);
listaexc=imatrix(1,NN+1,1,conexMAX+2); conextotalin=vector(1,NN+2);
listaini=imatrix(1,NN+1,1,conexMAX+2);

////// Nomeando arquivos //////////////////////////////////
  
//o=fopen("POK_A_delay0_gexc_g_rand0.dat","wt");  
p=fopen("RP.dat","wt");
p1=fopen("RPexc.dat","wt");
p2=fopen("RPinh.dat","wt");
q=fopen("Line_delay0_gexc_g_rand0.dat","wt"); 
  
if(o==NULL) { puts("erro no arquivo"); getchar(); exit(1); }
if(p==NULL) { puts ("erro no arquivo"); getchar(); exit(1);}
if(q==NULL){puts ("erro no arquivo"); getchar();exit(1);}

///////// definindo percentagem de conexão para cada componente////
///////// Numero de conexões = probabilidade * (Numero total de conexões possível)
N_aut_exc = 0.25 * (0.8 * NN);
N_aut_inh = 0.25 * (0.2 * NN);
N_exc_exc = 0.05 * (0.8 * NN) * (0.8 * NN-1);
N_inh_inh = 0.05 * (0.2 * NN) * (0.2 * NN-1);
N_exc_inh = 0.05 * (0.2 * NN) * (0.8 * NN) ; 
N_inh_exc = 0.05 * (0.2 * NN) * (0.8 * NN) ;
printf("N_exc_exc=%d,N_inh_inh=%d, N_aut_exc=%d,N_aut_inh=%d,N_exc_inh=%d,N_inh_exc=%d\n",
N_exc_exc, N_inh_inh,N_aut_exc, N_aut_inh,N_exc_inh, N_inh_exc);
/////////////////////////////////////////////////////////////////////////
///// Numero de conexões entre as populações inicialmente igual a zero.
Cont_aut_exc=0; Cont_aut_inh=0; Cont_exc_exc=0;
Cont_inh_inh=0; Cont_exc_inh=0; Cont_inh_exc=0;

//****************** Criando lista com as conexoes *********************
//-----------zerando tudo------------//  
for(i=1;i<=NN;i++)
{
conextotalex[i]=0.0;
for(j=1;j<=conexMAX;j++) listaexc[i][j]=0.0;
}

////////////////////////////////////////////////////////////////////////
//////////////////////montando as redes/////////////////////////////////
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
////////////////////////Rede exc-exc////////////////////////////////////
/////////--------pelo menos uma conexao para cada neuronio------////////  
////////////////////////////////////////////////////////////////////////

for(i=1;i<=NN*0.8;i++)   //// 1 até 800
{
j=(int)(NN*0.8*ran1(& idum));
if(j==0) j=NN*0.8;
conextotalex[i]=1.0; listaexc[i][1]=j; Cont_exc_exc++;	
}

////////////////// encontrando conexões restantes //////////////////////
	
while(Cont_exc_exc<N_exc_exc)
{
i=(int)NN*0.8*ran1(& idum); //i neurônio pós sinaptico
j=(int)NN*0.8*ran1(& idum);// j neurônio pré sinaptico

if(i==0) i=NN*0.8; if(j==0) j=NN*0.8; auxx=0.0;
      
for(jj=1;jj<=conextotalex[j];jj++) //conta apenas novas conexoes
if(listaexc[j][jj]==i) auxx=1.0;

if(auxx==0.0) if(j!=i) 
{ conextotalex[i]++; listaexc[i][conextotalex[i]]=j; Cont_exc_exc++;}
//fprintf(o,"%d %d\n", i, listaexc[i][conextotalex[i]]);
//printf("%d %d %d %d\n",  N_exc_exc, Cont_exc_exc,i,j);
}

///////////////////////////fim_exc_exc//////////////////////////////////

////////////////////////////////////////////////////////////////////////
//////////////////////////Rede Inh_Inh//////////////////////////////////
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////
/////////--------pelo menos uma conexao para cada neuronio------////////
////////////////////////////////////////////////////////////////////////
  
for(i=NN*0.8;i<=NN;i++)
{ 
j=((int)NN*(0.8+0.2*ran1(& idum)));
if(j==NN*0.8) j=NN;
	
conextotalex[i]=1.0; 	listaexc[i][1]=j;
Cont_inh_inh++;	
}

////////////////// encontrando conexões restantes //////////////////////

while(Cont_inh_inh<N_inh_inh)
{
i=(int)(NN*0.8+NN*0.2*ran1(& idum)); //i neurônio pós sinaptico
j=(int)(NN*0.8+NN*0.2*ran1(& idum));// j neurônio pré sinaptico

if(i==NN*0.8) i=NN; if(j==NN*0.8) j=NN; auxx=0.0;
      
for(jj=1;jj<=conextotalex[j];jj++) //conta apenas novas conexoes
if(listaexc[j][jj]==i) auxx=1.0;

if(auxx==0.0) 
if(j!=i) 
{conextotalex[i]++; listaexc[i][conextotalex[i]]=j; Cont_inh_inh++;}
//printf("%d %d %d %d\n",  N_inh_inh, Cont_inh_inh,i,j);
}
//////////////////////////////////////////////////////Fim rede inh-inhi

////////////////////////////////////////////////////////////////////////
///////// Autapses Exc//////////////////////////////////////////////////
while(Cont_aut_exc<N_aut_exc)
{
for(i=601;i<=800;i++)
{	
j=i;// j neurônio pré sinaptico

auxx=0.0;
      
for(jj=1;jj<=conextotalex[j];jj++) //conta apenas novas conexoes
if(listaexc[j][jj]==i) auxx=1.0;

if(auxx==0.0) 
if(j==i) 
{conextotalex[i]++; listaexc[i][conextotalex[i]]=i; Cont_aut_exc++;}
//printf("%d %d %d %d\n",  N_aut_exc, Cont_aut_exc,i,j);
//
}
}
///////////////////////////////////////////////////Fim rede autapse exc
////////////////////////////////////////////////////////////////////////
///////// Autapses inh//////////////////////////////////////////////////
while(Cont_aut_inh<N_aut_inh) 
{	
for(i=801;i<=850;i++)
{
j=i;// j neurônio pré sinaptico
 auxx=0.0;
      
for(jj=1;jj<=conextotalex[j];jj++) //conta apenas novas conexoes
if(listaexc[j][jj]==i) auxx=1.0;

if(auxx==0.0) if(j==i) 
{conextotalex[i]++; listaexc[i][conextotalex[i]]=i; Cont_aut_inh++; }
//printf("%d %d %d %d\n", N_aut_inh, Cont_aut_inh,i,j);
//fprintf(o,"%d %d\n", i, listaexc[i][conextotalex[i]]);
}
}

//////////////////////////fim rede autapse inh//////////////////////////

////////////////////////////////////////////////////////////////////////
////////////////////////// Rede exc_inh////////// Excc ----> Inh ///////
////////////////////////////////////////////////////////////////////////

while(Cont_exc_inh<N_exc_inh)
{
i=(int)(NN*0.8+NN*0.2*ran1(& idum)); //i neurônio pós sinaptico inh
j=(int)NN*0.8*ran1(& idum);// j neurônio pré sinaptico exc

if(i==N*0.8) i=NN; if(j==0)	j=NN*0.8; auxx=0.0;
      
for(jj=1;jj<=conextotalex[j];jj++) //conta apenas novas conexoes
if(listaexc[j][jj]==i) auxx=1.0;

if(auxx==0.0) 
if(j!=i) 
{conextotalex[j]++; listaexc[j][conextotalex[j]]=i;Cont_exc_inh++;}
//printf("%d %d %d %d\n",  N_exc_inh, Cont_exc_inh,i,j);
}

////////////////////////////Fim rede exc-inh////////////////////////////

////////////////////////////////////////////////////////////////////////
////////////// Rede inh_exc ////////Inh -----> Exc//////////////////////
////////////////////////////////////////////////////////////////////////

while(Cont_inh_exc<N_inh_exc)
{
i=(int)(NN*0.8*ran1(& idum)); //i neurônio pós sinaptico exc
j=(int)(NN*0.8+NN*0.2*ran1(& idum));// j neurônio pré sinaptico inh

if(j==N*0.8) j=NN; if(i==0)	i=NN*0.8; auxx=0.0;
      
for(jj=1;jj<=conextotalex[j];jj++) //conta apenas novas conexoes
if(listaexc[j][jj]==i) auxx=1.0;

if(auxx==0.0) 
if(j!=i) 
{conextotalex[j]++; listaexc[j][conextotalex[j]]=i; Cont_inh_exc++;}
//printf("%d %d %d %d\n",  N_inh_exc, Cont_inh_exc,i,j);
}

////////////////////////////////////////////Fim rede exc-inh////////////
//*****************************************************************/////

printf("N_exc_exc=%d,N_inh_inh=%d, N_aut_exc=%d,N_aut_inh=%d,N_exc_inh=%d,N_inh_exc=%d\n",
Cont_exc_exc, Cont_inh_inh,Cont_aut_exc,Cont_aut_inh,Cont_exc_inh,Cont_inh_exc);

/////////////////////////Parameters/////////////////////////////////////
gexc_exc=0.1; // [0- 0.7] /// Exc-Exc Gera sozinho gera sincronização///
gexc_inh=1.0; /////Exc-Inh, sozinho não gera sincronização na pop inh
ginh_exc=1.0;  //// Inh-Exc, sozinho não gera sincronização na pop exc.
gini_ini=1.0; //// Inh-Inh, sozinho não gera sincronização na pop exc.
gexc_auto=10.0; //// Gera burst nos nurônios com autosinapses
gini_auto=1.0;

step=0.0;
//for(gexc_exc=0.1*(variable-1);gexc_exc<0.1*variable;gexc_exc=gexc_exc+0.05)
//for(gini_ini=1.0*(variable-1);gini_ini<1.0*variable;gini_ini=gini_ini+1.0)
for(gexc_inh=1.0*(variable-1);gexc_inh<1.0*variable;gexc_inh=gexc_inh+1.0)
{
//for(gexc_auto=0.0;gexc_auto<=40.0;gexc_auto=gexc_auto+2.0)
//for(gini_auto=0.0;gini_auto<=300.0;gini_auto=gini_auto+30.0)
for(ginh_exc=0.0;ginh_exc<=10.0;ginh_exc=ginh_exc+1.0)
{
//////// Fixed parameters////////
{b_IF=70.0;  V_reset=-58.0; delay_exc=1.5; delay_ini=0.8;
C=200.0; gL=12.0; VT=-50.0; DeltaT=2.0; tauw=300.0; //ms
EL=-70.0; taum=C/gL; Vr_exc=0.0; Vr_ini=-80.0;}

for(rep=1;rep<=condi;rep++)
   {
///////////////////////////////////////////////////////////	   
////////////////////// Condições iniciais/////////////////
{
CVM=0;RM=0; IM2=0; IMI=0;FM=0; idum=-1234567897; ///// Semente
n=n_n; // total number of steps n*h
aux=0; ISI=0;  ISI2=0; contISI=0; F=0.0; IsynM=0;  
  
for(i=1;i<=NN;i++) aa[i]=1.9+0.2*ran1(& idum); } 


//********* Mais Condicoes iniciais das variáveis **************************
for(i=1;i<=NN;i++)
{
x[1+(i-1)*N]=-70.6+20.0*ran1(& idum); //V   
x[2+(i-1)*N]= 0.0+80.0*ran1(& idum);  //w 

Gexc_exc[i]=0*ran1(& idum);
Gexc_auto[i]=0*ran1(& idum);
Gini_ini[i]=0*ran1(& idum); 
Gini_auto[i]=0*ran1(& idum);
Gexc_inh[i]=0*ran1(& idum);
Ginh_exc[i]=0*ran1(& idum);                 //s
	 
for(j=1;j<=aloca;j++) marca_tempo[i][j]=-1000;
	 
conta_disparo[i]=0;
	       
//~ aa[i]=1.9+0.2*ran1(& idum); 
tpeak[i]=-100.0;	//tempo do último disparo
cont_Iext[i]=0.0;
k[i]=0.0;               // contador do numero de disparos de um neuronio
kmax[i]=0.0;	         // contador do total de disparos de um neuronio
nk[1][i]=0.0;           // tempo em que ocorrem os disparos de um neuronio    	
atualiza_disparo[i]=1.0;
tauex=2.728; //ms synaptic conductances time constants  
tauin=tauex;

	}
//********************* LOOP DO TEMPO  ***************************
tempo=0.0;
aux=0;cont_hyper=0;
IsynM=0.0;IsynMI=0.0;kkk_cont=0;
for(t=1;t<=n;t++)  
	{                   
	  tempo=tempo+h;      //em milisegundos
	  
	  for(i=1;i<=NN;i++)  // tempo anterior
	    { 
Gexc_exc[i]=Gexc_exc[i]*exp(-h/tauex);  // condutancia recebida	
Gini_ini[i]=Gini_ini[i]*exp(-h/tauin);
Gexc_auto[i]=Gexc_auto[i]*exp(-h/tauex);
Gini_auto[i]=Gini_auto[i]*exp(-h/tauin);
Gexc_inh[i]=Gexc_inh[i]*exp(-h/tauex);
Ginh_exc[i]=Ginh_exc[i]*exp(-h/tauin);

////////////////////////////////////////////////////////////////////////
//////////////////////RESET////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

if(x[1+(i-1)*N]>-35.0) // modelo de aIEF
		{
conta_disparo[i]=conta_disparo[i]+1.0;
			  
if(conta_disparo[i]>=aloca) conta_disparo[i]=1.0;
			  
//  delay_exc=5.0+0.1*gasdev(& idum);  
//  delay_ini=5.0+0.1*gasdev(& idum);
			  
if(i<=P_exc_ini*NN) delay=delay_exc;
      		else delay=delay_ini;
			  
marca_tempo[i][conta_disparo[i]]=tempo+delay;
   //fprintf(p,"%f %d\n", tempo ,i);  /////  raster plot
		      
if(tempo>transient && k[i]<NMD)     //salvando os disparos e o tempo em que ocorrem para cada neuronio
{ k[i]=k[i]+1; nk[k[i]][i]=tempo; 
if(tempo<transient+500.0) if(k[i]==1.0)   kinicial=tempo; }
	      		      
/////fprintf(p,"%f %d\n",tempo,i); // raster plot	
//	 if(i<=0.8*NN) fprintf(p1,"%f %d\n",tempo,i); // raster plot	exc
//	 else fprintf(p2,"%f %d\n",tempo,i); // raster plot	inh

if(tempo>transient)
  	  {
	  ISI=ISI+tempo-tpeak[i];
      ISI2=ISI2+(tempo-tpeak[i])*(tempo-tpeak[i]);
      contISI=contISI+1; 
      }
tpeak[i]=tempo; 
if(k[i]>=NMD) t=n+1;  
                 
x[1+(i-1)*N]=V_reset; //V
x[2+(i-1)*N]=x[2+(i-1)*N]+b_IF; //w	       // if(i<20) x[1+(i-1)*N]=-45;
		}
		
//printf("%.3f %.3f \n",delay_exc,delay_ini);
if(atualiza_disparo[i]>=aloca) atualiza_disparo[i]=1.0;
			
if(tempo>=marca_tempo[i][atualiza_disparo[i]] && tempo<marca_tempo[i][atualiza_disparo[i]]+h) 
      {	
   	if(i<=P_exc_ini*NN) for(j=1;j<=conextotalex[i];j++)
		  {    //// exc to exc   
		  if(listaexc[i][j]<=P_exc_ini*NN)
				{   //// one exc to another exc
				if(i!=listaexc[i][j])
					{
					Gexc_exc[listaexc[i][j]]=Gexc_exc[listaexc[i][j]]+1.0; // condutancia exc recebida						
					}
		         
				else
					{ ////one exc to the same exc
					Gexc_auto[listaexc[i][j]]=Gexc_auto[listaexc[i][j]]+1.0; // condutancia exc recebida
					//printf("Exc-Time=%f, Neuron %d\n", t,i);
					}
				}
					
	      if(listaexc[i][j]>P_exc_ini*NN)
				{ //// exc to inh
				Gexc_inh[listaexc[i][j]]=Gexc_inh[listaexc[i][j]]+1.0; // condutancia exc recebida
				//printf("Exc-Time=%f, Neuron %d\n", t,i);
				}  //// se o neurônio pré j de i, é inhibitório.
		  } ////// Se o neurônio pos i é excitatório.
	 else 
	for(j=1;j<=conextotalex[i];j++)
		 { 
		 if(listaexc[i][j]>P_exc_ini*NN)
				{
				if(i!=listaexc[i][j])
					{
					Gini_ini[listaexc[i][j]]=Gini_ini[listaexc[i][j]]+1.0; // condutancia exc recebida
					}
				else
					{
					Gini_auto[listaexc[i][j]]=Gini_auto[listaexc[i][j]]+1.0; // condutancia exc recebida
					//printf("Exc-Time=%f, Neuron %d\n", t,i);
					}
				}
				
		if(listaexc[i][j]<=P_exc_ini*NN)
			{
Ginh_exc[listaexc[i][j]]=Ginh_exc[listaexc[i][j]]+1.0; // condutancia exc recebida
//printf("%d %d\n",listaexc[i][j],i );
//printf("Exc-Time=%f, Neuron %d\n", t,i);
			}  //// se o neurônio pré j de i, é inhibitório.
		  }
               atualiza_disparo[i]=atualiza_disparo[i]+1.0;
	    }
	      y[1+(i-1)*N]=x[1+(i-1)*N];	
	      y[2+(i-1)*N]=x[2+(i-1)*N];	  
	      
if(tempo>transient) IsynM = IsynM + (
      + gexc_exc*Gexc_exc[i]*(Vr_exc-y[1+(i-1)*N]) 
       + gini_ini*Gini_ini[i]*(Vr_ini-y[1+(i-1)*N])  
        + gexc_auto*Gexc_auto[i]*(Vr_exc-y[1+(i-1)*N]) 
			+ gini_auto*Gini_auto[i]*(Vr_ini-y[1+(i-1)*N])    
			  + gexc_inh*Gexc_inh[i]*(Vr_exc-y[1+(i-1)*N]) 
	+ ginh_exc*Ginh_exc[i]*(Vr_ini-y[1+(i-1)*N]) 
)/NN;           
	    }
// ------------ Integrador numérico Runge-Kutta 4ª ordem--------------- 
 derivs(y,df,aa,gexc_exc,gini_ini,gexc_auto,gini_auto,gexc_inh,ginh_exc,
 Gexc_exc, Gini_ini,Gexc_inh,Ginh_exc,Gexc_auto,Gini_auto);
 for(i=1;i<=N*NN;i++)
	    {
	     a[i]=h*df[i];
	     y[i]=x[i]+a[i]/2.0;
	    }
 derivs(y,df,aa,gexc_exc,gini_ini,gexc_auto,gini_auto,gexc_inh,ginh_exc,
 Gexc_exc,Gini_ini,Gexc_inh,Ginh_exc,Gexc_auto,Gini_auto);
  for(i=1;i<=N*NN;i++)
	    {
	     b[i]=h*df[i];
	     y[i]=x[i]+b[i]/2.0;
	    }
derivs(y,df,aa,gexc_exc,gini_ini,gexc_auto,gini_auto,gexc_inh,ginh_exc,
Gexc_exc,Gini_ini,Gexc_inh,Ginh_exc,Gexc_auto,Gini_auto);
for(i=1;i<=N*NN;i++)
	    {
	     c[i]=h*df[i];
	     y[i]=x[i]+c[i]; 
	    }   
derivs(y,df,aa,gexc_exc,gini_ini,gexc_auto,gini_auto,gexc_inh,ginh_exc,
Gexc_exc,Gini_ini,Gexc_inh,Ginh_exc,Gexc_auto,Gini_auto);
	   for(i=1;i<=N*NN;i++)
	    x[i]=x[i]+(a[i]+h*df[i])/6.0+(b[i]+c[i])/3.0; 	
	  //--------------------------------------------------------------------		
//imprime todos os 10% dos valores de V em função do tempo
//~ fprintf(p,"%.2f",tempo);
// for(i=1;i<=N*NN/20.0;i=i+3)
//   fprintf(p," %.3f",x[i]);
//  fprintf(p,"\n"); 
	  
if(cont_hyper>200) for(i=1;i<=N;i++) if(x[i]<-100)
{ printf("Hiperpolariza\n");
printf("Parâmetros, gexc=%f, ginh=%f, gexc_auto=%f,ginh_auto=%f,ginh_exc=%f,gexc_inh=%f",
 gexc_exc,gini_ini,gexc_auto,   gini_auto,   
ginh_exc,gexc_inh);	   break; }
else cont_hyper=0;
cont_hyper++;

//    i=500;
//    if(kkk_cont>=1)  
//     {
//     fprintf(o,"%.3f %.3f %.3f %.3f %.3f %.3f %.3f\n",tempo,
//     gexc_exc*Gexc_exc[i]*(Vr_exc-y[1+(i-1)*N]), 
//       gini_ini*Gini_ini[i]*(Vr_ini-y[1+(i-1)*N]),  
//       gexc_auto*Gexc_auto[i]*(Vr_exc-y[1+(i-1)*N]), 
//			gini_auto*Gini_auto[i]*(Vr_ini-y[1+(i-1)*N]),    
//			  gexc_inh*Gexc_inh[i]*(Vr_exc-y[1+(i-1)*N]), 
//			    ginh_exc*Ginh_exc[i]*(Vr_ini-y[1+(i-1)*N])); 
//   fprintf(o,"%.3f %.3f %.3f %.3f %.3f %.3f %.3f\n",tempo, 
//   Gexc_exc[i],Gini_ini[i],Gexc_auto[i],Gini_auto[i],Gexc_inh[i],Ginh_exc[i]);  
//			  kkk_cont=0;  
//			    }
//   kkk_cont++;

}  // fim loop do tempo

// CALCULO DO PARAMETRO DE ORDEM //
///////////////////////////////////

////////////////////////////////////////
///////Zera condições iniciais/////////
/////////////////////////////////////////

tempo=0.0; Rmedio=0.0; contR=0.0; kfinal=n*h;

for(i=1;i<=NN;i=i+1)
{
phi[i]=0.0;   kmax[i]=k[i];   k[i]=1.0;	
if(kmax[i]>1) if(kfinal>nk[kmax[i]][i] && kfinal>kinicial)
  kfinal=nk[kmax[i]][i];
}

for(tempo=transient;tempo<=n*h;tempo=tempo+50.0*h)
{
for(i=1;i<=NN;i=i+1) if(kmax[i]>1) if(tempo>nk[k[i]][i] && tempo<nk[kmax[i]][i])
   {
   if(tempo<=nk[k[i]+1][i])	  
   phi[i]=2*PI*(tempo-nk[k[i]][i])/(nk[k[i]+1][i]-nk[k[i]][i]); //+2*PI*k[i]
        
    if(tempo>=nk[k[i]+1][i]) k[i]=k[i]+1;		         
    }

  if(tempo>=kinicial && tempo<=kfinal)
  {
  real=0.0; compl=0.0; cont=0; 
  for(i=1;i<=NN;i++) if(kmax[i]>1)
      {
      real=real+cos(phi[i]);  compl=compl+sin(phi[i]);  cont=cont+1;
      }
    
real=real/cont; compl=compl/cont; R=sqrt(real*real+compl*compl);
    Rmedio=Rmedio+R; contR=contR+1.0;
    //fprintf(o,"%f %f \n",tempo,R);
  }
}
     //fprintf(o,"\n");
///////////////////////////////////////////////////////////
// CALCULO DO CV 
desviopadrao=sqrt(ISI2/contISI - (ISI/contISI)*(ISI/contISI));
CV=desviopadrao/(ISI/contISI);  CV_ant=CV;  F=((1000.0)/(ISI/contISI));
  
  aux=1.0;     
  // if(desviopadrao>6.0 && CV>0.5)
  //aux=-1;
  
  CVM=CVM+CV*1.0/condi;  FM=FM+ F*1.0/condi;
  RM= RM+aux*Rmedio/(contR*condi);
  IM2=IM2 + IsynM*0.01/((tempo-transient)*condi);
  IMI=IMI + IsynMI*0.01/((tempo-transient)*condi);
} // fim loop rep

fprintf(o,"%.3f %.3f %.3f %.3f %.3f %.3f %.6f %.3f %.3f %.3f\n",
gexc_exc,gexc_inh,gexc_auto,gini_ini,ginh_exc,gini_auto,FM,CVM,aux*Rmedio/contR,IM2);

  
	
printf("g_e=%.3f, g_ei=%.3f, g_e_auta=%.3f, g_i=%.3f, g_ie=%.3f, g_i_auta=%.3f, F=%.3f, CV=%.3f, Rmedio= %.3f, IM= %.3f \n",
gexc_exc,gexc_inh,gexc_auto,gini_ini,ginh_exc,gini_auto,FM,CVM,aux*Rmedio/contR,IM2);	
} /////// fim segundo loop
fprintf(o,"\n");
} //////// fim primeiro loop

////Dispensa memórias////////
free_dvector(y,1,N*NN+1);   free_dvector(df,1,N*NN+1);  
free_dvector(x,1,N*NN+1);   free_dvector(a,1,N*NN+1);
free_dvector(b,1,N*NN+1);  free_dvector(c,1,N*NN+1);
free_dvector(cont_Iext,1,NN+1);   free_dvector(aa,1,NN+1);
free_vector(k,1,NN+1);   free_vector(kmax,1,NN+1);
free_dvector(phi,1,NN+1);   free_dmatrix(nk,1,NMD+2,1,NN+2);
free_dvector(tpeak,1,NN+1); free_dvector(Gexc_exc,1,NN+1); 
free_dvector(Gini_ini,1,NN+1); free_dvector(Gexc_auto,1,NN+1); 
free_dvector(Gini_auto,1,NN+1); free_dvector(Gexc_inh,1,NN+1); 
free_dvector(Ginh_exc,1,NN+1);  free_vector(conextotalex,1,NN+2);
free_imatrix(listaexc,1,NN+1,1,conexMAX+2);

///// Fecha arquivos .dat ///
fclose(o); fclose(p); fclose(p1); fclose(p2); fclose(q);
return 0;
}

void derivs(double y[],double df[],double *aa, double gexc_exc,
double gini_ini, double gexc_auto,double gini_auto,double gexc_inh,
double ginh_exc,double *Gexc_exc,double *Gini_ini, double *Gexc_inh,
double *Ginh_exc, double *Gexc_auto,double *Gini_auto) // Equacoes diferenciais acopladas
{
int i;
double C,gL,VT,DeltaT,tauw,Vr_exc,Vr_ini,tauex,tauin,EL,taum,I;
I=270;  C=200.0;  gL=12.0;  VT=-50.0;   
DeltaT=2.0; tauw=300.0;  EL=-70.0;
taum=C/gL;  Vr_exc=0.0;  Vr_ini=-80.0;  
//////////////////////////EQUAÇOES ACOPLADAS////////////////////

for(i=1;i<=NN;i++)
{
 //printf("%d %f\n", i , I[i]);

  df[1+(i-1)*N]=(1.0/C)*(-gL*(y[1+(i-1)*N]-EL)
  + gL*DeltaT*exp((y[1+(i-1)*N]-VT)/DeltaT) 
    - y[2+(i-1)*N] +I 
      + gexc_exc*Gexc_exc[i]*(Vr_exc-y[1+(i-1)*N]) 
       + gini_ini*Gini_ini[i]*(Vr_ini-y[1+(i-1)*N])  
        + gexc_auto*Gexc_auto[i]*(Vr_exc-y[1+(i-1)*N]) 
			+ gini_auto*Gini_auto[i]*(Vr_ini-y[1+(i-1)*N])    
			  + gexc_inh*Gexc_inh[i]*(Vr_exc-y[1+(i-1)*N]) 
	+ ginh_exc*Ginh_exc[i]*(Vr_ini-y[1+(i-1)*N])                                                          
					);
                                                                        
   df[2+(i-1)*N]=(1.0/tauw)*(aa[i]*(y[1+(i-1)*N]-EL) - y[2+(i-1)*N]);
   } 
}

float ran1(long *idum)
{
 int j;
 long k;
 static long iy=0;
 static long iv[NTAB];
 float temp;
 
 if(*idum<=0 || !iy)
   {
     if(-(*idum)<1) *idum=1;
     else *idum = -(*idum);
    for(j=NTAB+7;j>=0;j--)
      {
       k=(*idum)/IQ;
       *idum=IA*(*idum-k*IQ)-IR*k;
       if(*idum<0) *idum +=IM;
       if(j<NTAB) iv[j]=*idum;
      }
      iy=iv[0];
   }
   k=(*idum)/IQ;
   *idum=IA*(*idum-k*IQ)-IR*k;
   if(*idum<0) *idum += IM;
   j=iy/NDIV;
   iy=iv[j];
   iv[j]=*idum;
   if((temp=AM*iy)>RNMX) return RNMX;
   else return temp;
}

float gasdev(long *idum)
{
  float ran1(long *idum);
  static int iset=0;
  static float gset;
  float fac,rsq,v1,v2;  
  
  if(*idum<0) iset=0;
  if(iset==0) {
    do {
      v1=2.0*ran1(idum)-1.0;
      v2=2.0*ran1(idum)-1.0;
      rsq=v1*v1+v2*v2;
    } while (rsq>=1.0 || rsq==0.0);
    fac=sqrt(-2.0*log(rsq)/rsq);
    gset=v1*fac;
    iset=1;
    return v2*fac;
  } else {
    iset=0;
    return gset;
  }
} 

double *dvector(long nl,long nh)
{
   double *v;
   
   v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
   if (!v) nrerror("allocation failure in dvector()");
   return v-nl+NR_END;
}

void free_dvector(double *v, long nl, long nh)
{
   free((FREE_ARG) (v+nl-NR_END));
}

void nrerror(char error_text[])
{
fprintf(stderr,"Numerical Recipes run-time error...\n");
fprintf(stderr,"%s\n",error_text);
fprintf(stderr,"...now exiting to system...\n");
exit(1);
}

int *vector(long nl,long nh)
{
int *v;
   
v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
if (!v) nrerror("allocation failure in dvector()");
return v-nl+NR_END;
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
{
long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
int **m;

m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
if (!m) nrerror("allocation failure 1 in imatrix()");
m += NR_END;
m -= nrl;

m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
if (!m[nrl]) nrerror("allocation failure 2 in imatrix()");
m[nrl] += NR_END;
m[nrl] -= ncl;

for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

return m;
}

void free_vector(int *v, long nl, long nh)
{
free((FREE_ARG) (v+nl-NR_END));
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
{
free((FREE_ARG) (m[nrl]+ncl-NR_END));
free((FREE_ARG) (m+nrl-NR_END));
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
{
long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
double **m;

m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
if (!m) nrerror("allocation failure 1 in matrix()");
m += NR_END;
m -= nrl;

m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
m[nrl] += NR_END;
m[nrl] -= ncl;

for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;
return m;
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
{
free((FREE_ARG) (m[nrl]+ncl-NR_END));
free((FREE_ARG) (m+nrl-NR_END));
}
