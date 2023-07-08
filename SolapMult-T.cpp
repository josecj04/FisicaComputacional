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
void calcular_pesos(double pesos[N][N][N][N],int psi[N][N][4],double a[4]);
void calcular_theta(double pesos[N][N][N][N],double theta[N][N]);
void calcular_solapamiento(double solapamiento[4], int red[N][N],int psi[N][N][4],double a[4]);
double calcular_unsolapamiento(int red[N][N], int imagenfinal[N][N], double a_final);

using namespace std;

int main()
{
	double T= 0.0001; //Temperatura
	int n,m;
	double z;
	int i,j,k,l,u,o,h,f;
	int red[N][N]; //Matriz de neuronas
	double pesos[N][N][N][N]; //Para almacenar los pesos sinápticos
	double theta[N][N]; //Umbral de disparo
	extern gsl_rng *tau; //declaracion del random
	//8732358 //1493270 //7912864 //3472479
	int semilla= 1493270; //2475675	 //semilla necesaria para generar números
	double energia;
	double p,e,aux;
	double solapamiento[4],unsolap;
	int aux2;
	int psi[N][N][4],imagenfinal[N][N];//para almacenar el fichero en binario
	double a[4],a_final;

	//Inicializo los ficheros
    ifstream fich1; // Fichero con la imagen arriba izquierda
    fich1.open("uno.dat");

    ifstream fich2; // Fichero con la imagen arriba derecha
    fich2.open("tres.dat");

    ifstream fich3; // Fichero con la imagen querida abajo izquierda
    fich3.open("cinco.dat");

    ifstream fich4; // Fichero con la imagen querida abajo derecha
    fich4.open("diez.dat");


    ofstream imagen; // Fichero con la imagen dada por la red neuronal
    imagen.open("Hopfieldmult_data.dat");

    ofstream solap1;  // Fichero del solapamiento en funcion de pasos montecarlo
    solap1.open("solapamientomult-T(I).dat");

	ofstream solap2;  // Fichero del solapamiento en funcion de pasos montecarlo
    solap2.open("solapamientomult-T(III).dat");

	ofstream solap3;  // Fichero del solapamiento en funcion de pasos montecarlo
    solap3.open("solapamientomult-T(V).dat");

	ofstream solap4;  // Fichero del solapamiento en funcion de pasos montecarlo
    solap4.open("solapamientomult-T(X).dat");

	ifstream inicial; // Fichero con la imagen querida en binario
    inicial.open("imagenmult.dat");

	ofstream solapimagenin; // Fcihero con el solapamiento de la imgen entera
	solapimagenin.open("solapamientomult-T(completo).dat");

	ofstream animacion;
	animacion.open("imagenfinalHopfield.dat");

	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			
			inicial >> red[i][j];
			imagenfinal[i][j]= red[i][j];
			cout<<imagenfinal[i][j];
			a_final= a_final + imagenfinal[i][j];
		}
		cout<<endl;
	}
	a_final= a_final/(N*N);
	cout<<a_final<<endl;
	


	tau = gsl_rng_alloc(gsl_rng_taus); //necesario para lo de los números aleatorios
	gsl_rng_set(tau,semilla); // Damos la semilla al puntero de números aleatorios
	
	
	//iniciar(red); //Iniciamos nuestra red neuronal con 0 y 1 aleatoriamente
	
	
	i=0;
	j=0;
	a[0]=0;
	a[1]=0;
	a[2]=0;
	a[3]=0;

  
// Cargamos la primera imagen en psi1
    for ( i =0; i < N; i++)
    {
        for ( j = 0; j < N; j++)
        {
            fich1>>psi[i][j][0];
			
			a[0]=a[0]+ psi[i][j][0];
            
        }
        
        
    }
    
 	a[0] = a[0]*1.0/(N*N);
	cout<<"EL valor de a[1] es "<<a[0]<<endl;


	// Cargamos la segunda imagen en psi2
    for ( i =0; i < N; i++)
    {
        for ( j = 0; j < N; j++)
        {
            fich2>>psi[i][j][1];
			
			a[1]=a[1]+ psi[i][j][1];            
        }
        
        
    }
    
 	a[1] = a[1]*1.0/(N*N);
	cout<<"EL valor de a[2] es "<<a[1]<<endl;


	// Cargamos la tercera imagen en psi3
    for ( i =0; i < N; i++)
    {
        for ( j = 0; j < N; j++)
        {
            fich3>>psi[i][j][2];
			
			
			a[2]=a[2]+ psi[i][j][2];
            
        }
		
        
        
    }
    
 	a[2] = a[2]*1.0/(N*N);
	cout<<"EL valor de a[3] es "<<a[2]<<endl;

	// Cargamos la cuarta imagen en psi1
    for ( i =0; i < N; i++)
    {
        for ( j = 0; j < N; j++)
        {
            fich4>>psi[i][j][3];

			

			
			a[3]=a[3]+ psi[i][j][3];
            
        }

		
        
        
    }
    
 	a[3] = a[3]*1.0/(N*N);
	cout<<"EL valor de a[4] es "<<a[3]<<endl;
 
	iniciarf(red); //Iniciamos la red a un valor de 15% cambio

	
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


	while (T <= 0.5000) 
	{


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
                        animacion<<red[u][o]<<", ";
                    }
                    else animacion<<red[u][o];
			}
			animacion<<endl;
		}
		animacion<<endl;
		}*/

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
	 
		calcular_solapamiento(solapamiento, red, psi, a);

		//Imprimimos los solapamientos en un fichero 

		
		if(j=99){
	 solap1<<T<<" ";
	 solap1<<solapamiento[0]<<endl;

	 solap2<<T<<" ";
	 solap2<<solapamiento[1]<<endl;

	 solap3<<T<<" ";
	 solap3<<solapamiento[2]<<endl;

	 solap4<<T<<" ";
	 solap4<<solapamiento[3]<<endl;
	 }
		

		//Calculamos el solapamiento con la imagen final
		unsolap=calcular_unsolapamiento(red, imagenfinal, a_final);

		
		if(j=99){
	 solapimagenin<<T<<" ";
	 solapimagenin<<unsolap<<endl;
	 }
		



	 
	}

	T=T+0.001;
    
	}

    //Cerramos los ficheros utilizados

    fich1.close();
	fich2.close();
	fich3.close();
	fich4.close();
	imagen.close();
    solap1.close();
	solap2.close();
	solap3.close();
	solap4.close();
	inicial.close();
	solapimagenin.close();

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

			if (m<0.20){
				if (red[i][j]==1){

					red[i][j]=0;
					}
				else red[i][j]=1;
			}

	     
	     
		}
	}
	return;
}




void calcular_pesos(double pesos[N][N][N][N],int psi[N][N][4],double a[4]) //Función para calcular los pesos sinápticos
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

			for ( o = 0; o < 4; o++)
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

void calcular_solapamiento(double solapamiento[4],int red[N][N],int psi[N][N][4],double a[4])
{
 
 int i,j,u;
 double aux4;
 

for ( u = 0; u < 4; u++)
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


double calcular_unsolapamiento(int red[N][N], int imagenfinal[N][N], double a_final){

	int i,j;
 double aux4,solapa;
 aux4=0;
 for(i=0;i<N;i++)
 {
   for(j=0;j<N;j++)
   {
     aux4 = aux4+(imagenfinal[i][j]-a_final)*(red[i][j]-a_final);
   }
 }
 
 
 solapa = 1.0*aux4*1.0/(N*N*a_final*(1.0-a_final));





	return solapa; 
}