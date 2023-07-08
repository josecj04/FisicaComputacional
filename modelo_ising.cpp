//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////// CALCULO MODELO DE ISING CON ALGORITMO DE METROPOLI /////////////////////////////
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
using namespace std;

#define N 100

double CalcularEnergia(int spin[N][N], int n,int m);



int main(){

    int n, m;
    int i, j, k, l, numpos, numneg;
    int spin[N][N];//Matrices de spin 1 o -1
    double Energia, p, z, aux, T;
    T=3;

    


//Abrimos el fichero "ising.dat"
ofstream fich;
fich.open("ising_data.dat");

//Iniciamos la matriz de spin con valores -1 y 1

    for ( i = 0; i < N; i++)
    {
        for ( j = 0; j < N; j++)
        {
            spin[i][j]= rand() % 2;
            
            if (spin[i][j]==0)
            {
                spin[i][j]=-1;
            }
            
        }
        
    }

//Ponemos la matriz todo de 1

    for ( i = 0; i < N; i++)
    {
        for ( j = 0; j < N; j++)
        {
            spin[i][j]=1;
        }
        
    }




//Realizamos los calculos de algoritmo de metropoli para el cálculo de la matriz de spin

    for ( k = 0; k < 1000; k++)
    {
        numpos=0;
        numneg=0;

        //Imprimimos la matriz de spin el el fichero
            for ( i = 0; i < N; i++)
            {
                for ( j = 0; j < N; j++)
                {    
                    if (j<N-1)
                    {
                        fich<<spin[i][j]<<", ";
                    }
                    else fich<<spin[i][j];
                      
                     
                }
                fich<<endl;
            } 
            
            fich<<endl;
        
        //Iniciamos un paso de metropoli
        for ( l = 0; l < N*N; l++)
        {
            //Calculamos la energia de un punto random de la red cristalina
            n= rand() % N;
            m= rand() % N;
            Energia=CalcularEnergia(spin , n, m);
            
            

            if (1.0< exp(-1.0*(Energia/T)))
            {
                p=1.0;
            }
            else p= exp(-1.0*(Energia/T)); 

         

            //Generamos un numero aleatorio entre 0 y 1
            z= rand() / (1.0*RAND_MAX);


            if (z<p) // Si ese número es menor que nuestra p entonces cambiamos el sentido del spin
	        {
	   	        spin[n][m] = -spin[n][m];
	        }

        }


        for ( i = 0; i < N; i++)
        {
                for (j = 0; j < N; j++)
                {
                    if (spin[i][j]==1)
                    {
                        numpos += 1;
                    }
                    else numneg += 1;
                    
                }
                
        }
             

    cout<<"El numero total de 1 es la red en el paso de Metropoli "<<k<<" es: "<<numpos<<endl;
    cout<<"El numero total de -1 es la red en el paso de Metropoli "<<k<<" es: "<<numneg<<endl;
    cout<<endl;

        
            
        
        
    }
    
    




    fich.close();





}







/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////FUNCIONES USADAS///////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////


double CalcularEnergia(int spin[N][N], int n,int m){
    double energia;
	int aux,aux2,aux3,aux4; 
	//Lo siguiente son condiciones de contorno periódicas
	aux= n+1;
	aux2=n-1;
	aux3 = m+1;
	aux4 = m-1;
    if (aux== N)
	{
	  aux = 0;	
	}
	if (aux2== -1)
	{
	  aux2 = N-1;	
	}
	if (aux3== N)
	{
	  aux3 = 0;	
	}
	if (aux4==-1)
	{
	  aux4 = N-1;	
	}

    //Aplicamos la formula de la energia en el algoritmo de metropoli
    energia = 2.0*(spin[n][m]*(spin[aux][m]+spin[aux2][m]+spin[n][aux3]+spin[n][aux4]));


	return energia;



}
 


