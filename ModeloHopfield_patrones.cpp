#include <iostream>
#include <math.h>
#include <fstream>
#include <ostream>
#include <string>
#include <stdlib.h>
#include <cstdlib>
#include <ctime>
#include "gsl_rng.h"

#define N 20// dimension red de neuronas.
#define num 3//Numero de patrones almacenados
gsl_rng *tau; //declaración de puntero número aleatorio de forma global 





double calcular_energia(int red[N][N],int n,int m,double pesos[N][N][N][N],double theta[N][N]);
void iniciar(int red[N][N]);
void iniciarpatron(int psi[N][N][num]);
void iniciarpatronf(int psi[N][N][num]);
void iniciarf(int red[N][N]);
void pintar_red (int red[N][N]);
void calcular_pesos(double pesos[N][N][N][N],int psi[N][N][num],double a[num]);
void calcular_theta(double pesos[N][N][N][N],double theta[N][N]);
void calcular_solapamiento(double solapamiento[num], int red[N][N],int psi[N][N][num],double a[num]);
double calcular_unsolapamiento(int red[N][N], int imagenfinal[N][N], double a_final);

using namespace std;

int main()
{
	double T= 0.0001; //Temperatura
	int n,m;
	double z;
	int i,j,k,l,u,o,h,f,v;
	int red[N][N]; //Matriz de neuronas
	double pesos[N][N][N][N]; //Para almacenar los pesos sinápticos
	double theta[N][N]; //Umbral de disparo
	extern gsl_rng *tau; //declaracion del random
	//8732358 //1493270 //7912864 //3472479
	int semilla= 2475675; //2475675	 //semilla necesaria para generar números
	double energia,numpatrones,alpha;
	double p,e,aux;
	double solapamiento[num],unsolap;
	int aux2;
	int psi[N][N][num],imagenfinal[N][N];//para almacenar el fichero en binario
	double a[num],a_final;
	unsigned t0,t1;

	//Inicializo los ficheros
   
    ofstream psis;
    psis.open("psi.dat");

    ofstream imagen; // Fichero con la imagen dada por la red neuronal
    imagen.open("Hopfieldmult_data.dat");

    

	

	ofstream patrones; // Fichero número de patrones recordados
    patrones.open("patrones.dat");

	
	
	t0=clock();

	tau = gsl_rng_alloc(gsl_rng_taus); //necesario para lo de los números aleatorios
	gsl_rng_set(tau,semilla); // Damos la semilla al puntero de números aleatorios
	
	
	iniciar(red); //Iniciamos nuestra red neuronal con 0 y 1 aleatoriamente
	
	numpatrones=0;
    alpha,0.0;
	i=0;
	j=0;
	a[0]=0;
	a[1]=0;
	a[2]=0;
	a[3]=0;


  
// Cargamos el numero de patrones en psi

    iniciarpatron(psi);
    iniciarpatronf(psi);
 	

    for ( int v = 0; v < num; v++)
    {
		

        for ( int i = 0; i < N; i++)
        {
            for ( int j = 0; j < N; j++)
            {
                psis<<psi[i][j][v];
                psis<<" ";

				a[v]=a[v]+psi[i][j][v];

            }
            psis<<endl;
            
        }
        psis<<endl;

		 a[v]= 1.0*a[v]/(N*N);
		 cout<<"El valor de a["<<v<<"] es: "<<a[v]<<endl;
    }
   
    
	



	

	
 	
 	calcular_pesos(pesos,psi,a);
	calcular_theta(pesos,theta);


	//Vamos a hacer la evolucion del sistema usando el algoritmo de metropoli
	for (j=0;j<100;j++) // Pasos Montecarlo
	{

		//Imprimimos la imagen en el fichero Hopfield_data.dat
		/*if(j==99){
		for ( u = 0; u < N; u++)
		{
			for ( o = 0; o < N; o++)
			{
				if (o<N-1)
                    {
                        imagen<<red[u][o]<<", ";
                    }
                    else imagen<<red[u][o];
			}
			imagen<<endl;
		}
		imagen<<endl;
		}
        */
	
	  for(i=0;i<N*N;i++) //N*N = N^2 es un paso montecarlo,número de veces a realizar el algoritmo para dar 1 paso en 1 fotograma 				     //del gift
	  {
	   
	    n = gsl_rng_uniform_int(tau,N); 
	    m = gsl_rng_uniform_int(tau,N); //Escogemos aleatoriamente una neurona en la posición n,m
	    

	    energia= calcular_energia(red,n,m,pesos,theta); //Calculamos su energia
	    

	    if (1.0< exp(-1.0*(energia/T)))
	    {
	   	p = 1.0;
	    }
	    else p= exp(-1.0*(energia/T)); //Escogemos entre el mínimo de 1 y e^-deltaE/T
	    

	    z= gsl_rng_uniform(tau); //Generamos un número aleatorio entre 0 y 1

	    
	    if (z<p) // Si ese número es menor que nuestra p entonces cambiamos el estado de la neurona
	    {

	   	 red[n][m] = 1-red[n][m];

	    }

	  }
	 
		
		

		

		
		



	 
	}

    calcular_solapamiento(solapamiento, red, psi, a);

        

		for (int u = 0; u < num; u++)
		{
			if (solapamiento[u]>0.75 || solapamiento[u]<-0.75)
			{
				numpatrones= numpatrones+1;
			}
			cout<<"El solapamiento ["<<u<<"] es: "<<solapamiento[u]<<endl;
		}
		
		alpha=numpatrones/(N*N);


        //imprimimos alpha en el fichero patones.dat

		cout<<"El numero de patrones recordados es : "<<numpatrones<<endl;
		cout<<"El alpha es : "<<alpha<<endl;
    
	

    //Cerramos los ficheros utilizados

    
	
	patrones.close();
    imagen.close();


	t1=clock();

	double time = (double(t1-t0)/CLOCKS_PER_SEC);
	cout << "Execution Time: " << time << endl;



	return 0;
}




//Función para calcular delta de Energia
double calcular_energia(int red[N][N],int n,int m,double pesos[N][N][N][N],double theta[N][N])
{
	double energia;
	double aux; 
	int i,j;
	aux= 0;
	
	for (i=0;i<N;i++)
	{
	 for(j=0;j<N;j++)
	 {

	     aux = pesos[n][m][i][j]*red[i][j]+aux;
	   
	 }
	}
	energia = (1.0-2.0*red[n][m])*(theta[n][m]-0.5*pesos[n][m][n][m]-aux);
	return energia;
	
}


void iniciar(int red[N][N]) //Funcion para dar un valor aleatorio a mi red neuronal
{
	extern gsl_rng *tau;
	int i,j;

	
	for (i=0;i<N;i++)
	{
	  for (j=0;j<N;j++)
	   {
	     red[i][j] = gsl_rng_uniform_int(tau,2); //Le doy un valor aleatorio o bien 0 o 1
	    }
	}
	return;
}


void iniciarf(int red[N][N]) //Funcion para dar un valor a la imagen inicial pero con un error
{
	extern gsl_rng *tau;
	int i,j;
	double m;

	



	
	for (i=0;i<N;i++)
	{
	  for (j=0;j<N;j++)
	   {
		
		m= gsl_rng_uniform(tau);

			if (m<0.50){
				if (red[i][j]==1){

					red[i][j]=0;
					}
				else red[i][j]=1;
			}

	     
	     
		}
	}
	return;
}




void calcular_pesos(double pesos[N][N][N][N],int psi[N][N][num],double a[num]) //Función para calcular los pesos sinápticos
{

	
 double  aux;
 int i,j,k,l,o;
 
//Calculo los pesos:

for(i=0;i<N;i++)
{
  for(j=0;j<N;j++)
    {
	for(k=0;k<N;k++)
	{
	  for(l=0;l<N;l++)
	  {
	    if (i==k && j ==l)
	    {
	    	pesos[i][j][k][l] = 0.0;
	    }
	    else
	    {
			pesos[i][j][k][l]=0.0;

			for ( o = 0; o < num; o++)
			{
				pesos[i][j][k][l] = pesos[i][j][k][l] + (psi[i][j][o]-a[o])*(psi[k][l][o]-a[o]);
			}
			
			pesos[i][j][k][l]= pesos[i][j][k][l]/(N*N);
	      
	    }
	  }
	}
    }
}
return;
}

void calcular_theta(double pesos[N][N][N][N],double theta[N][N]) //Función para calcular el umbral de disparo
{

 int i,j,k,l;
 
for(i=0;i<N;i++)
{
  for(j=0;j<N;j++)
    {
	for(k=0;k<N;k++)
	{
	  for(l=0;l<N;l++)
	  {
	    theta[i][j] = 0.5*(pesos[i][j][k][l])+theta[i][j];
	  }
	}
    }
}
return;
}

void calcular_solapamiento(double solapamiento[num],int red[N][N],int psi[N][N][num],double a[num])
{
 
 int i,j,u;
 double aux4;
 

for ( u = 0; u < num; u++)
{
	aux4=0;
	


 	for(i=0;i<N;i++)
 	{
   		for(j=0;j<N;j++)
   		{
     		aux4 = aux4+(psi[i][j][u]-a[u])*(red[i][j]-a[u]);
   		}
 	}

	solapamiento[u]= 1.0*aux4*1.0/(N*N*a[u]*(1.0-a[u]));

    

}
 
 

return;
}


void iniciarpatron(int psi[N][N][num]) //Funcion para dar un valor aleatorio a mi red neuronal
{
	extern gsl_rng *tau;
	int i,j,l;

	for ( l = 0; l < num; l++)
    {
        
    
	for (i=0;i<N;i++)
	{
	  for (j=0;j<N;j++)
	   {
	     psi[i][j][l] = gsl_rng_uniform_int(tau,2); //Le doy un valor aleatorio o bien 0 o 1
	    }
	}
    }

	return;
}


void iniciarpatronf(int psi[N][N][num]) //Funcion para dar un valor a la imagen inicial pero con un error
{
	extern gsl_rng *tau;
	int i,j,l;
	double m;

	



	for (int  l = 0; l < num; l++)
    {
    
	for (i=0;i<N;i++)
	{
	  for (j=0;j<N;j++)
	   {
		
		m= gsl_rng_uniform(tau);

			if (m<0.65){
				if (psi[i][j][l]==1){

					psi[i][j][l]=0;
					}
				else psi[i][j][l]=1;
			}

            


            

	     
	     
		}
            
	}
    
    
    
    }

	return;
}