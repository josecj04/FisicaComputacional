//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////// CALCULO ECUACION DE SCHRODINGER EN EC. DIFERENCIALES ///////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <math.h>
#include <fstream>
#include <ostream>
#include <string>
#include <stdlib.h>
#include <conio.h>
#include <cstdlib>
#include <ctime>
#include <complex>
using namespace std;

//PARAMETROS
#define PI 3.141592653
#define N 100
#define nciclos 1
#define lambda 5



int main(){

int i,j;
	double k_0; 
	double aux,auxiliarcoseno,auxiliarseno;
	double auxi;
	double Amenos, Amas,compro;
	double V[N+1]; // Potencial reescalado
	double s; // "discretizacion espacial reescalada"
	double modulo[N+1]; //modulo de phi
	double norma; //norma
	std::complex<double> auxiliaralpha;
	std::complex<double> alpha[N+1];
    std::complex<double> phi[N+1]; //Funcion de onda
	std::complex<double> auxb;
    std::complex<double> b[N+1];
	std::complex<double> auxiliarbeta;
	std::complex<double> beta[N+1];
	std::complex<double> aux2;
	std::complex<double> aux3;
	std::complex<double> aux4;
	std::complex<double> Acero;
	std::complex<double> chi[N+1];

    //Abrimos ficheros de inicio ofstream
    ofstream fich;
    fich.open("schrodinger_data.dat");
   

	k_0=(2.0*PI*nciclos)/(1.0*N);
	cout<<k_0<<endl;

    for (i=0;i<=N;i++)
	{
		if(i>=0.4*N && i<=0.6*N)
		{
			
			
			V[i] = 1.0*lambda*k_0*k_0; 

				
			
		}
		else V[i]=0;
		//cout<<"El valor de V[ "<< i<<" ]: "<<V[i]<<endl;
		
	}
	



    //Aplicamos las condiciones de contorno para la funcion de onda en el instante inicial
    phi[0]=complex(0.0,0.0);
	phi[N]=complex(0.0,0.0);

	norma=0.0;


    
    s = 1.0/(4.0*k_0*k_0);

	//Iniciar la función de onda según la fórmula para tiempo inicial y su norma 
	for (i=1;i<= N-1;i++) 
	{
		
		aux2=complex(cos(k_0*i),sin(k_0*i));
		//cout<<aux2<<endl;
		aux=exp(-8.0*pow((4.0*i-N), 2.0)/(N*N));
		//cout<<aux<<endl;
		
		phi[i]=complex(aux*aux2);

		//cout<<"El valor de phi[ "<< i<<" ]: "<<phi[i]<<endl;
		
		norma+=norm(phi[i]);

		
	
		
	}
	cout<<"La norma total para este estado es "<<norma<<endl;

	

	//Normalizamos el estado
	for ( i = 0; i <= N-1; i++)
	{
		phi[i]= complex(((phi[i])/ sqrt(norma)));
		//cout<<"El valor de phi[ "<< i<<" ]: "<<phi[i]<<endl;
		
		

	}
	

	//Condiciones a imponer sobre las alphas
	alpha[N]=complex(0.0,0.0);
	alpha[N-1]=complex(0.0,0.0);

	//Calculamos los alphas correspondientes independientes del tiempo
	for ( i = N-2; i >= 0; i--)
	{
		Acero=complex(0.0,0.0);
		Acero=complex(-2.0-V[i],2.0/s);
		//cout<<"El valor de Acero es en la iteracion "<<i<<" : "<<Acero<<endl;
		auxiliaralpha =complex(Acero + alpha[i+1]);
		//cout<<"El valor de auxiliaralpha es en la iteracion "<<i<<" : "<<auxiliaralpha<<endl;
		alpha[i]=complex(-1.0/auxiliaralpha);
	}
	

	


	//Avanzamos en el bucle del tiempo 
	for ( j = 0; j < 1000; j++)
	{

		//Calculamos la norma total de la onda y el modulo de la onda en cada uno de sus nodos
		norma=0.0;

		for (i=0;i<=N;i++)
		{
			modulo[i] = abs(phi[i]);

			

			norma += norm(phi[i]);
		}

		cout<< "La norma para el instante t "<<j<<" es igual a : "<<norma<<endl;

		//Imprimo los modulo de cada nodo de la onda y el potencial 

		for (i=0;i<=N;i++)
		{
			fich<<i<<", ";
			fich<<modulo[i]<<", ";
			fich<<V[i]<<endl;	
		}
		fich<<endl;
		
		
		
		


		//Calculamos los valores de b_j
		for (i=0;i<=N;i++)
		{
			auxb = complex(0.0,4.0/s);

			b[i]= complex(auxb*phi[i]);
			
		}


		//	Calculo los beta
		beta[N-1]= complex(0.0,0.0); //Aplicamos las condiciones de contorno

		for ( i = N-2; i >= 0; i--)
		{
			Acero=complex(0.0,0.0);
			Acero=complex(-2.0-V[i],2.0/s);
			auxiliarbeta =complex(Acero + alpha[i+1]);

			beta[i]= complex((b[i+1]-beta[i+1])/auxiliarbeta);

		}
	

		//Calculamos los chi de la ecuación

		chi[0]= complex(0.0,0.0); //Condicion de contorno de chi

		for (i=0;i<=N-2;i++)
		{
			chi[i+1]= complex((alpha[i]*chi[i])+beta[i]);
		}

		chi[N]= complex(0.0,0.0); //Otra condicion de contorno de chi


		//Calculamos las phi para el nuevo t

		for (i = 0; i < N; i++)
		{
			phi[i]= complex(chi[i]-phi[i]);
			//cout<<"Valor de phi para j: "<<j<<" e i: "<<i<<" es: "<<phi[i]<<endl;
		}

		

	}

	fich.close();


    return 0;
}