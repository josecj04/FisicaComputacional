#include <iostream>
#include <math.h>
#include <fstream>
#include <ostream>
#include <string>
#include <stdlib.h>
#include <cstdlib>
#include <ctime>
#include "gsl_rng.h"

#define N 30 // dimension red de neuronas
gsl_rng *tau; //declaración de puntero número aleatorio de forma global 

double calcular_energia(int red[N][N],int n,int m,double pesos[N][N][N][N],double theta[N][N]);
void iniciar(int red[N][N]);
void iniciarf(int red[N][N]);
void pintar_red (int red[N][N]);
void calcular_pesos(double pesos[N][N][N][N],int psi[N][N],double a);
void calcular_theta(double pesos[N][N][N][N],double theta[N][N]);
double calcular_solapamiento(int red[N][N],int psi[N][N],double a);

using namespace std;

int main()
{
	double T= 0.0001; //Temperatura
	int n,m;
	double z;
	int i,j,k,l,u,o;
	int red[N][N]; //Matriz de neuronas
	double pesos[N][N][N][N]; //Para almacenar los pesos sinápticos
	double theta[N][N]; //Umbral de disparo
	extern gsl_rng *tau; //declaracion del random
	//8732358 //1493270 //7912864 //3472479
	int semilla=2475675; //2475675 //1493270	 //semilla necesaria para generar números
	double energia;
	double p,e,aux,a;
	double solapamiento;
	int aux2;
	int psi[N][N];//para almacenar el fichero en binario
    double psi1[N][N];
	unsigned t0,t1;
	//Inicializo los ficheros
    ifstream fich; // Fichero con la imagen querida en binario
    fich.open("imagen.dat");

    ofstream imagen; // Fichero con la imagen dada por la red neuronal
    imagen.open("Hopfield_data.dat");

    ofstream solap;  // Fichero del solapamiento en funcion de pasos montecarlo
    solap.open("solapamiento_data0,0001.dat");

	ofstream animacion;//FIchero donde guardar la imagen final tras los pasos montecarlos
	animacion.open("imagenfinalHopfield.dat");



	t0=clock();


	tau = gsl_rng_alloc(gsl_rng_taus); //necesario para lo de los números aleatorios
	gsl_rng_set(tau,semilla); // Damos la semilla al puntero de números aleatorios
	
	
	
	
	
	i=0;
	j=0;
	a=0;



  

    for ( i =0; i < N; i++)
    {
        for ( j = 0; j < N; j++)
        {
            fich>>psi[i][j];
			red[i][j]=psi[i][j];
            cout<<psi[i][j];
			
            a=a+psi[i][j];
        }
        
        cout<<endl;
    }
    
 	a = a*1.0/(N*N);
	cout<<"EL valor de a es "<<a<<endl;
 
	//iniciarf(red); //Iniciamos la red a un valor de 25% cambio

	iniciar(red); //Iniciamos nuestra red neuronal con 0 y 1 aleatoriamente

	//Mostramos la red inicial por pantalla

	for ( i = 0; i < N; i++)
	{
		for (j = 0; j < N; j++)
		{
			cout<<red[i][j];
		}
		cout<<endl;
		
	}
	
for ( u = 0; u < N; u++)
		{
			for ( o = 0; o < N; o++)
			{
				if (o<N-1)
                    {
                        animacion<<red[u][o]<<", ";
                    }
                    else animacion<<red[u][o];
			}
			animacion<<endl;
		}
		animacion<<endl;
	
 	
 	calcular_pesos(pesos,psi,a);
	calcular_theta(pesos,theta);


	//Vamos a hacer la evolucion del sistema usando el algoritmo de metropoli
	for (j=0;j<100;j++) // Pasos Montecarlo
	{

		//Imprimimos la imagen en el fichero Hopfield_data.dat



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


		/*if(j==99){
		for ( u = 0; u < N; u++)
		{
			for ( o = 0; o < N; o++)
			{
				if (o<N-1)
                    {
                        animacion<<red[u][o]<<", ";
                    }
                    else animacion<<red[u][o];
			}
			animacion<<endl;
		}
		animacion<<endl;
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
	 
     solapamiento=calcular_solapamiento(red,psi,a);
	 
	 solap<<solapamiento<<endl;
	 
	}
    
	

    //Cerramos los ficheros utilizados

    fich.close();
	imagen.close();
    solap.close();

	t1 = clock();


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

			if (m<0.25){
				if (red[i][j]==1){

					red[i][j]=0;
					}
				else red[i][j]=1;
			}

	     
	     
		}
	}
	return;
}


void pintar_red (int red[N][N]) //función para ver como es mi red de spins
{
 int i,j;
 for(i=0;i<N;i++)
 {
   for(j=0;j<N;j++)
   {
    printf("%i",red[i][j]);	
   }
   printf("\n");
 } 
  
 return;
}


void calcular_pesos(double pesos[N][N][N][N],int psi[N][N],double a) //Función para calcular los pesos sinápticos
{

	
 double  aux;
 int i,j,k,l;
 
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
	      pesos[i][j][k][l] = 1.0*(psi[i][j]-a)*(psi[k][l]-a)*(1.0/(N*N));
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

double calcular_solapamiento(int red[N][N],int psi[N][N],double a)
{
 
 int i,j;
 double aux4,solapamiento;
 aux4=0;
 for(i=0;i<N;i++)
 {
   for(j=0;j<N;j++)
   {
     aux4 = aux4+(psi[i][j]-a)*(red[i][j]-a);
   }
 }
 
 
 solapamiento = 1.0*aux4*1.0/(N*N*a*(1.0-a));
return solapamiento;
}
