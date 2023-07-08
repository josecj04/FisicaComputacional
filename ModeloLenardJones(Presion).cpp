//EJERCICIO VOLUNTARIO: POTENCIAL DE LENNARD-JONES

#include <iostream>
#include <math.h>
#include <fstream>
#include <ostream>
#include <string>
#include <stdlib.h>
#include <cstdlib>
#include <ctime>
#include "gsl_rng.h"


#define N 16
#define L 10.0
#define h 0.002
#define PI 3.141592654
#define iteraciones 10000

//Para reescalar la velocidad

#define incremento 7.0

gsl_rng *tau;

using namespace std;

int main()
{
    int i,j,k;
    double x[N],y[N],vx[N],vy[N],ax[N],ay[N],wx[N],wy[N],phi[N],delta[N];
    double auxx[N],auxy[N],dx,dy,mod,F,F0,TC,p,px,py,PR,auxPR,auxTC;
    int naux;
    double t=0.0;
    int semilla=18233247;

    ofstream f1;
    ofstream f2;
    
    f1.open("lj.txt");
    f2.open("presionytemp.txt");

    tau=gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(tau,semilla);

    cout<<endl;

    //CONDICIONES INICIALES

    //===========================================================================================

    auxPR=0.0;
    auxTC=0.0;
    naux=0;

    for (i = 0; i < N; i++)
    {
        //Auxiliar aleatorio entre 0 y 1.5

        delta[i]=(gsl_rng_uniform(tau))*1.5; 

        //Ángulo aleatorio entre 0 y 2PI

        phi[i]=(gsl_rng_uniform(tau))*2*PI; 

        //Velocidad aleatoria con módulo igual a 1

        vx[i]=cos(phi[i])*incremento;
        vy[i]=sin(phi[i])*incremento;

        //Inicializo a cero 

        ax[i]=0.0;
        ay[i]=0.0;

        auxx[N]=0.0;
        auxy[N]=0.0;

    }

    //Posiciones iniciales

    x[0]=1.0+delta[0];
    x[1]=4.0+delta[1];
    x[2]=5.0+delta[2];
    x[3]=8.0+delta[3];
    x[4]=1.0+delta[4];
    x[5]=4.0+delta[5];
    x[6]=5.0+delta[6];
    x[7]=8.0+delta[7];
    x[8]=1.0+delta[8];
    x[9]=4.0+delta[9];
    x[10]=5.0+delta[10];
    x[11]=8.0+delta[11];
    x[12]=1.0+delta[12];
    x[13]=4.0+delta[13];
    x[14]=5.0+delta[14];
    x[15]=8.0+delta[15];

    y[0]=1.0+delta[0];
    y[1]=1.0+delta[1];
    y[2]=1.0+delta[2];
    y[3]=1.0+delta[3];
    y[4]=3.0+delta[4];
    y[5]=3.0+delta[5];
    y[6]=3.0+delta[6];
    y[7]=3.0+delta[7];
    y[8]=5.0+delta[8];
    y[9]=5.0+delta[9];
    y[10]=5.0+delta[10];
    y[11]=5.0+delta[11];
    y[12]=7.0+delta[12];
    y[13]=7.0+delta[13];
    y[14]=7.0+delta[14];
    y[15]=7.0+delta[15];

    //Evalúo la aceleración y energías iniciales

    F0=(24.0/pow(3.0,7.0))*((2.0/pow(3.0,6.0))-1); 

    for (i=0; i<N-1; i++) 
        {
            for(j=i+1; j<N; j++) 
            {
                dx = (x[i]-x[j]);
                dy = (y[i]-y[j]);
                
                if (fabs(dx)>0.5*L) dx=(L-fabs(dx))*((-dx)/fabs(dx));
                if (fabs(dy)>0.5*L) dy=(L-fabs(dy))*((-dy)/fabs(dy));

                mod=sqrt(pow(dx,2.0)+pow(dy,2.0));
                
                if(mod<3.0)
                {
                    F=(24.0/pow(mod,7.0))*((2.0/pow(mod,6.0))-1);  

                    ax[i]=ax[i]+(F-F0)*(dx/mod);
                    ax[j]=ax[j]-(F-F0)*(dx/mod);
                    ay[i]=ay[i]+(F-F0)*(dy/mod);
                    ay[j]=ay[j]-(F-F0)*(dy/mod);
                }
            }
        }    

    //===========================================================================================

    //ALGORITMO DE VERDET

    //===========================================================================================

    for (k = 0; k < iteraciones; k++)
    {  
        t=t+h;

        p=0.0;
        px=0.0;
        py=0.0;

        PR=0.0;

        TC=0.0;

        //Primero evaluo x,y,wx,wy//

        for (i = 0; i < N; i++)
        {
            wx[i]=(vx[i]+h*ax[i]/2.0);
            wy[i]=(vy[i]+h*ay[i]/2.0);
        }

        for (i = 0; i < N; i++)
        {
            x[i]=x[i]+h*wx[i];
            y[i]=y[i]+h*wy[i];  
        }

        //Verifico las condiciones de contorno periodicas

        for (i = 0; i < N; i++)
        {
            if (x[i]>L)
            {
                x[i]=x[i]-L;
                px=px+2.0*fabs(vx[i]);
            } 
            if (x[i]<0.0)
            {
                x[i]=x[i]+L;
                px=px+2.0*fabs(vx[i]);
            } 

            if (y[i]>L)
            {
                y[i]=y[i]-L;
                py=py+2.0*fabs(vy[i]);
            }
            if (y[i]<0.0)
            {
                y[i]=y[i]+L;
                py=py+2.0*fabs(vy[i]);
            } 
        }

        p=sqrt(px*px+py*py);

        PR=(p/(L*L*h));

        auxPR=auxPR+PR;

        for (i = 0; i < N; i++)
        {
            ax[i]=0.0;
            ay[i]=0.0;
        }

        //Evaluo las nuevas aceleraciones

        for (i=0; i<N-1; i++) 
        {
            for(j=i+1; j<N; j++) 
            {
                dx = (x[i]-x[j]);
                dy = (y[i]-y[j]);
                
                if (fabs(dx)>0.5*L) dx=(L-fabs(dx))*((-dx)/fabs(dx));
                if (fabs(dy)>0.5*L) dy=(L-fabs(dy))*((-dy)/fabs(dy));

                mod=sqrt(pow(dx,2.0)+pow(dy,2.0));
                
                if(mod<3.0)
                {
                    F=(24.0/pow(mod,7.0))*((2.0/pow(mod,6.0))-1);  

                    ax[i]=ax[i]+(F-F0)*(dx/mod);
                    ax[j]=ax[j]-(F-F0)*(dx/mod);
                    ay[i]=ay[i]+(F-F0)*(dy/mod);
                    ay[j]=ay[j]-(F-F0)*(dy/mod);
                }
            }
        }
        
        //Por último evalúo la velocidad//

        for (i = 0; i < N; i++)
        {
            vx[i]=wx[i]+(h/2.0)*ax[i];
            vy[i]=wy[i]+(h/2.0)*ay[i];
        }

        //Temperatura y presión del sistema 


        for (i = 0; i < N; i++)
        {
            TC=TC+(vx[i]*vx[i]+vy[i]*vy[i]);
        }

        TC=TC/(2.0*N);

        auxTC=auxTC+TC;

        f2<<TC<<" "<<PR<<endl;
    
        //Repetimos
    }

    //===========================================================================================

    auxPR=(auxPR/(1.0*iteraciones));
    auxTC=(auxTC/(1.0*iteraciones));

    cout<<auxTC<<" "<<auxPR<<endl;


    cout<<endl;
    cout<<endl;
    cout<<"Terminado";
    cout<<endl;
    cout<<endl;

    f1.close();
    f2.close();

    return 0;
}

