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



#define N 20
#define L 10.0
#define h 0.002
#define PI 3.141592654
#define iteraciones 1500

//Para reescalar la velocidad

#define incremento 1.0

gsl_rng *tau;

using namespace std;

//En este programa no utilizaremos funciones ya que reciclaremos el código del sistema solar

int main()
{
    int i,j,k;
    double x[N],y[N],vx[N],vy[N],ax[N],ay[N],wx[N],wy[N],phi[N],delta[N],V,T,F0,V0,E,TC;
    double aux,aux2,aux3,aux4,aux4x,aux4y,auxv[N],auxvx[N],auxvy[N],dx,dy,mod,F;
    double t=0.0;
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
    f9.open("presionytemp.txt");

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

    for (i = 0; i < N; i++)
    {
        //Auxiliar aleatorio entre 0 y 1.5

        delta[i]=(gsl_rng_uniform(tau))*1.5; 

        //Ángulo aleatorio entre 0 y 2PI

        phi[i]=(gsl_rng_uniform(tau));

        //Velocidad aleatoria con módulo igual a 1

        vx[i]=(phi[i])*incremento;
        vy[i]=0.0;

        //Inicializo a cero 

        ax[i]=0.0;
        ay[i]=0.0;

        auxv[i]=0.0;
        auxvx[i]=0.0;
        auxvy[i]=0.0;
    }

    for (i = 0; i < N; i++)
    {
        auxv[i]=sqrt((vx[i]*vx[i]+vy[i]*vy[i]));
        auxvx[i]=sqrt((vx[i]*vx[i]));
        auxvy[i]=sqrt((vy[i]*vy[i]));
        f6 << auxv[i] << endl;
        f7 << auxvx[i] << endl;
        f8 << auxvy[i] << endl;
        
        
    }

    f6 << endl;
    f7 << endl;
    f8 << endl;

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
    x[16]=1.0+delta[16];
    x[17]=4.0+delta[17];
    x[18]=5.0+delta[18];
    x[19]=8.0+delta[19];

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
    y[16]=9.0+delta[16];
    y[17]=9.0+delta[17];
    y[18]=9.0+delta[18];
    y[19]=9.0+delta[19];

    //Evalúo la aceleración y energías iniciales

    F0=(24.0/pow(3.0,7.0))*((2.0/pow(3.0,6.0))-1); 
    V0=(4.0/pow(3.0,6.0))*((1.0/pow(3.0,6.0))-1); 

    for (i = 0; i < N-1; i++) 
        {
            for(j = i+1; j < N; j++) 
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

    //Temperatura inicial

    for (i = 0; i < N; i++)
    {
        T=T+(vx[i]*vx[i]+vy[i]*vy[i]);
    }

    T=T/(2.0*N);

    f3 << t <<" ";
    f3 << T<<endl;
  
        
    //Imprimo las posiciones y energía inicial en el fichero

    for (j = 0; j < N; j++)
    {
        f1<<x[j]<<" "<<y[j]<<endl;
        
    }

    f1<<endl;

    f2<<t<<" "<<aux2;
    f4<<t<<" "<<aux;
    f2<<t<<" "<<E;
    

    for (i = 0; i < N; i++)
    {
        auxv[i]=(vx[i]*vx[i]+vy[i]*vy[i]);
    }

    //===========================================================================================

    //ALGORITMO DE VERDET

    //===========================================================================================

    for (k = 0; k < iteraciones; k++)
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

        //Temperatura del sistema para t=20 y 50

        if (k==10000 || k==25000)
        {
            for (i = 0; i < N; i++)
            {
                TC=TC+(vx[i]*vx[i]+vy[i]*vy[i]);
                auxv[i]=sqrt((vx[i]*vx[i]+vy[i]*vy[i]));
                auxvx[i]=sqrt((vx[i]*vx[i]));
                auxvy[i]=sqrt((vy[i]*vy[i]));
                f6<<auxv[i]<<endl;
                f7<<auxvx[i]<<endl;
                f8<<auxvy[i]<<endl;
                
            }


            f6<<endl;
            f7<<endl;
            f8<<endl;
            

            TC=TC/(2.0*N);

            f3<<t<<" "<<TC<<endl;
            

        }
        

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
    f9.close();

    return 0;
}

