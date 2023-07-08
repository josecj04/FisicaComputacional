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
    double dxaux,dyaux,modaux;
    double x0,y0;
    int nh,nhx,nhy;
    int semilla=18237240;

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
    f9.open("fluctuaciondosparticulas.txt");

    tau=gsl_rng_alloc(gsl_rng_taus);
    gsl_rng_set(tau,semilla);

    cout<<endl;

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

        vx[i]=0.0;
        vy[i]=0.0;

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
        delta[i]=(gsl_rng_uniform(tau))*0.6;
    }

    x[0]=0.5;
    x0=0.5;
    x[1]=1.5;
    x[2]=2.5;
    x[3]=3.5;
    x[4]=0.5;
    x[5]=1.5;
    x[6]=2.5;
    x[7]=3.5;
    x[8]=0.5;
    x[9]=1.5;
    x[10]=2.5;
    x[11]=3.5;
    x[12]=0.5;
    x[13]=1.5;
    x[14]=2.5;
    x[15]=3.5;

    y[0]=0.5;
    y0=0.5;
    y[1]=0.5;
    y[2]=0.5;
    y[3]=0.5;
    y[4]=1.5;
    y[5]=1.5;
    y[6]=1.5;
    y[7]=1.5;
    y[8]=2.5;
    y[9]=2.5;
    y[10]=2.5;
    y[11]=2.5;
    y[12]=3.5;
    y[13]=3.5;
    y[14]=3.5;
    y[15]=3.5;

    dxaux = fabs(x[0]-x[1]);
    dyaux = fabs(y[0]-y[1]);

    if (fabs(dxaux)>0.5*L) dxaux=(L-fabs(dxaux));
    if (fabs(dyaux)>0.5*L) dyaux=(L-fabs(dyaux));

    modaux=(pow(dxaux,2.0)+pow(dyaux,2.0));

    f9<<t<<" "<<modaux<<endl;
   

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

                    ax[i]=ax[i]+(F-F0)*(dx/mod);
                    ax[j]=ax[j]-(F-F0)*(dx/mod);
                    ay[i]=ay[i]+(F-F0)*(dy/mod);
                    ay[j]=ay[j]-(F-F0)*(dy/mod);
                }
            }
        }

    //Temperatura inicial

    for (i = 0; i < N; i++)
    {
        T=T+(vx[i]*vx[i]+vy[i]*vy[i]);
    }

    T=T/(2.0*N);

    f3<<t<<" "<<T<<endl;



    //Imprimo las posiciones y energía inicial en el fichero

    for (j = 0; j < N; j++)
    {
        f1<<x[j]<<" "<<y[j]<<endl;
    }

    f1<<endl;


    //===========================================================================================

    //ALGORITMO DE VERDET

    //===========================================================================================

    for (k = 0; k < 130000; k++)
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

        if (k==30000 || k==60000 || k==90000)
        {
            for (i = 0; i < N; i++)
            {
                vx[i]=1.1*vx[i];
                vy[i]=1.1*vy[i];
            }
        }
        
        dxaux = fabs(x[0]-x[1]);
        dyaux = fabs(y[0]-y[1]);

        if (fabs(dxaux)>0.5*L) dxaux=(L-fabs(dxaux));
        if (fabs(dyaux)>0.5*L) dyaux=(L-fabs(dyaux));

        modaux=sqrt(pow(dxaux,2.0)+pow(dyaux,2.0));

        f9<<t<<" "<<pow(modaux,2.0)<<endl;
        

        //Imprimo los resultados en el fichero

        for (i = 0; i < N; i++)
        {
            f1<<x[i]<<" "<<y[i];
        }

        f1<<endl;
    
        //Repetimos

        for (i = 0; i < N; i++)
        {
            T=T+(vx[i]*vx[i]+vy[i]*vy[i]);
        }

        T=T/(2.0*N);

        f3<<t<<" "<<T<<endl;
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
    f9.close();

    return 0;
}