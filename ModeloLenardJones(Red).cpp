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
#define L 4.0
#define h 0.002
#define PI 3.141592654

gsl_rng *tau;

using namespace std;

int main()
{
    int i,j,k;
    double x[N],y[N],vx[N],vy[N],ax[N],ay[N],wx[N],wy[N],phi[N],delta[N],V,T,F0,V0,E,TC;
    double aux,aux2,aux3,aux4,aux4x,aux4y,auxv[N],auxvx[N],auxvy[N],dx,dy,mod,F,max,maxx,maxy,min,minx,miny;
    double t=0.0;
    int nh,nhx,nhy;
    int semilla=18233247;

    ofstream f1;
    ofstream f2;
    ofstream f3;
    ofstream f4;
    ofstream f5;
    ofstream f6;
    ofstream f7;
    ofstream f8;
    ofstream f9;

    f1.open("lj.txt");
    f2.open("cinetica.txt");
    f3.open("temperatura.txt");
    f4.open("potencial.txt");
    f5.open("total.txt");
    f6.open("histograma.txt");
    f7.open("histogramax.txt");
    f8.open("histogramay.txt");

    tau=gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(tau,semilla);

    printf("\n");

    //CONDICIONES INICIALES

    //===========================================================================================

    //Inicializo a cero 

    V=0.0;
    T=0.0;
    E=0.0;
    TC=0.0;
    aux=0.0;
    aux2=0.0;
    aux3=0.0;
    aux4=0.0;
    aux4x=0.0;
    aux4y=0.0;
    max=0.0;
    min=10.0;
    maxx=0.0;
    minx=10.0;
    maxy=0.0;
    miny=10.0;
    T=0.0;
    nh=0;
    nhx=0;
    nhy=0;

    for (i = 0; i < N; i++)
    {
        

        //Ángulo aleatorio entre 0 y 2PI

        phi[i]=(gsl_rng_uniform(tau))*2*PI; 

        //Velocidad aleatoria con módulo igual a 1

        vx[i]=0;
        vy[i]=0;

        //Inicializo a cero 

        ax[i]=0.0;
        ay[i]=0.0;

        auxv[i]=0.0;
        auxvx[i]=0.0;
        auxvy[i]=0.0;
    }

    //Posiciones iniciales

    for (i = 0; i < N; i++)
    {
        delta[i]=(gsl_rng_uniform(tau))*0.1; 
    }
    

    x[0]=0.5+delta[0];
    x[1]=1.5+delta[1];
    x[2]=2.5+delta[2];
    x[3]=3.5+delta[3];
    x[4]=0.5+delta[4];
    x[5]=1.5+delta[5];
    x[6]=2.5+delta[6];
    x[7]=3.5+delta[7];
    x[8]=0.5+delta[8];
    x[9]=1.5+delta[9];
    x[10]=2.5+delta[10];
    x[11]=3.5+delta[11];
    x[12]=0.5+delta[12];
    x[13]=1.5+delta[13];
    x[14]=2.5+delta[14];
    x[15]=3.5+delta[15];


    y[0]=0.5+delta[0];
    y[1]=0.5+delta[1];
    y[2]=0.5+delta[2];
    y[3]=0.5+delta[3];
    y[4]=1.5+delta[4];
    y[5]=1.5+delta[5];
    y[6]=1.5+delta[6];
    y[7]=1.5+delta[7];
    y[8]=2.5+delta[8];
    y[9]=2.5+delta[9];
    y[10]=2.5+delta[10];
    y[11]=2.5+delta[11];
    y[12]=3.5+delta[12];
    y[13]=3.5+delta[13];
    y[14]=3.5+delta[14];
    y[15]=3.5+delta[15];

    //Evalúo la aceleración y energías iniciales

    F0=(24.0/pow(3.0,7.0))*((2.0/pow(3.0,6.0))-1); 
    V0=(4.0/pow(3.0,6.0))*((1.0/pow(3.0,6.0))-1); 

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
                    V=(4.0/pow(mod,6.0))*((1.0/pow(mod,6.0))-1)-V0+(mod*F0)-(3.0*F0);

                    aux=aux+V;    

                    ax[i]=ax[i]+(F-F0)*(dx/mod);
                    ax[j]=ax[j]-(F-F0)*(dx/mod);
                    ay[i]=ay[i]+(F-F0)*(dy/mod);
                    ay[j]=ay[j]-(F-F0)*(dy/mod);
                }
            }
        }
    
    for (i = 0; i < N; i++)
    {
        aux2=aux2+0.5*(vx[i]*vx[i]+vy[i]*vy[i]);
    }

    E=aux+aux2;


    //Imprimo las posiciones y energía inicial en el fichero

    for (j = 0; j < N; j++)
    {
        f1<<x[j]<<" "<<y[i]<<endl;
    }

    f1<<endl;



    f2<<t<<" "<<aux2<<endl;
    f4<<t<<" "<<aux<<endl;
    f5<<t<<" "<<E<<endl;



    //===========================================================================================

    //ALGORITMO DE VERDET

    //===========================================================================================

    for (k = 0; k < 3000; k++)
    {  

        t=t+h;

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
            if (x[i]>L) x[i]=x[i]-L;
            if (x[i]<0.0) x[i]=x[i]+L;
            if (y[i]>L) y[i]=y[i]-L;
            if (y[i]<0.0) y[i]=y[i]+L;
        }

        //Vuelvo a inicializar a cero

        V=0.0;
        T=0.0;
        E=0.0;
        TC=0.0;
        aux=0.0;
        aux2=0.0;
        aux3=0.0;
        aux4=0.0;
        aux4x=0.0;
        aux4y=0.0;
        T=0.0;

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
                    V=(4.0/pow(mod,6.0))*((1.0/pow(mod,6.0))-1)-V0+(mod*F0)-(3.0*F0);;

                    aux=aux+V;

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

        //Ahora evalúo las nuevas energías

        for (i = 0; i < N; i++)
        {
            aux2=aux2+0.5*(vx[i]*vx[i]+vy[i]*vy[i]);
        }

        E=aux+aux2;
        

        //Imprimo los resultados en el fichero

        for (i = 0; i < N; i++)
        {
            f1<<x[i]<<" "<<y[i]<<endl;
        }

        f1<<endl;

        f2<<t<<" "<<aux2<<endl;
        f4<<t<<" "<<aux<<endl;
    f5<<t<<" "<<E<<endl;

        //Repetimos
    }

    //===========================================================================================

     cout<<endl;
    cout<<endl;
    cout<<"Terminado"<<endl;
    cout<<endl;
    cout<<endl;
    

    f1.close();
    f2.close();
    f3.close();
    f4.close();
    f5.close();
    f6.close();
    f7.close();
    f8.close();

    return 0;
}
