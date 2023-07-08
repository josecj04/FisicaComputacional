//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////// CALCULO TRAYECTORIA COHETE - METODO DE RUNGE KUTTA  ////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <math.h>
#include <fstream>
#include <ostream>
#include <string>
#include <stdlib.h>
#include <conio.h>

//CONSTANTES 
#define PI 3.141592653
#define G 6.67*pow(10,-11.0)
#define Mt 5.9736*pow(10,24.0)
#define Ml 0.07349*pow(10,24.0)
#define dtl 3.844*pow(10,8.0)
#define w 2.6617*pow(10,-6.0)
#define Rt 6.378160*pow(10,6.0)
#define Rl 1.7374*pow(10,6.0)
#define m 20*pow(10,3.0)


using namespace std;

double CalcularH(double, double, double, double);

int main(){

    double h;
	h = 30;

	double r,phi,pr,pphi; //Coordenadas
	double raux, praux, pphiaux, rprimaaux; //Coordenadas auxiliarles para calcular el hamiltoniano
	double k1r,k1phi,k1pr,k1pphi;
	double k2r,k2phi,k2pr,k2pphi;
	double k3r,k3phi,k3pr,k3pphi;
	double k4r,k4phi,k4pr,k4pphi; //ki para todas las coordenadas
	double aux,aux1;
	int i,j;
	double delta,mu,rprima,rprimaH; // constantes
	double t; //tiempo
	
	double xTierra,yTierra; //Posciones x e y de la Tierra,cohete y Luna.  
	double xcohete,ycohete;
	double xLuna, yLuna;
	double v0,theta;
	double H, H_prima;
	
	ofstream fich;
	ofstream f;

	fich.open("cohete_data.dat");
	f.open("Hamiltoniano_data.dat");


    //Valores iniciales
    t =0;
	r = Rt*pow(dtl,-1.0); //posicion reescalado
	cout<<"El radio de la tierra es "<<Rt<<endl;
	cout<<"La distancia de la tierra a la luna es "<<dtl<<endl;
	cout<<"La r inicial es "<<r<<endl;
	phi =PI/2.0;

	v0= 11.0835*1000*pow(dtl,-1.0); //velocidad reescalada
	cout<<"El valor inicial de la velocidad es "<<v0<<endl;
	theta = 90*PI/180.0;
	delta=  (G*Mt)*pow(dtl,-3.0);
	cout<<"El valor inicial de la delta es "<<delta<<endl;
    mu= Ml*pow(Mt,-1.0);
	cout<<"El valor inicial de la mu es "<<mu<<endl;
	pr= v0*cos(theta-phi);
	pphi = r*v0*sin(theta-phi);


    xLuna = 1;
	yLuna = 0;
	xcohete = r;
	ycohete = 0;

    cout<<endl;
	cout<<endl;

    j=0;

    for ( i = 0; i < 10000; i++)
    { 
       //Calculamos K1
	    k1r = h*pr;
		k1phi = h*pphi/(r*r);
	
		rprima = pow(1+r*r-2*r*cos(phi-w*t),0.5);
		k1pr = h*((pphi*pphi)/pow(r,3)-delta*(1/(r*r)+(mu/pow(rprima,3.0))*(r-cos(phi-w*t))));
	
		k1pphi = h*(-1.0)*delta*mu*(r/pow(rprima,3.0))*sin(phi-w*t);


	   //Calculamos K2

	    k2r = h*(pr+k1pr*0.5);
		k2phi = h*(pphi+k1pphi*0.5)/((r+k1r*0.5)*(r+k1r*0.5));

		rprima = pow(1+(r+k1r*0.5)*(r+k1r*0.5)-2*(r+k1r*0.5)*cos((phi+k1phi*0.5)-w*(t+h*0.5)),0.5);
			
		aux = (pphi+k1pphi*0.5)*(pphi+k1pphi*0.5)/(pow((r+k1r*0.5),3.0)) ;
		aux1 = 1.0/((r+k1r*0.5)*(r+k1r*0.5));
		
		k2pr = h*(aux-delta*(aux1+(mu/pow(rprima,3))*((r+k1r*0.5)-cos((phi+k1phi*0.5)-w*(t+h*0.5)))));
	
		k2pphi = h*(-1.0)*delta*mu*((r+k1r*0.5)/pow(rprima,3.0))*sin((phi+k1phi*0.5)-w*(t+h*0.5));

	   //Calculamos K3

	   	k3r = h*(pr+k2pr*0.5);
		k3phi = h*(pphi+k2pphi*0.5)/((r+k2r*0.5)*(r+k2r*0.5));
	
		rprima = pow(1+(r+k2r*0.5)*(r+k2r*0.5)-2*(r+k2r*0.5)*cos((phi+k2phi*0.5)-w*(t+h*0.5)),0.5);
		
		aux = (pphi+k2pphi*0.5)*(pphi+k2pphi*0.5)/pow((r+k2r*0.5),3.0) ;
		aux1 = 1/((r+k2r*0.5)*(r+k2r*0.5));
		k3pr = h*(aux-delta*(aux1+(mu/(pow(rprima,3.0)))*((r+k2r*0.5)-cos((phi+k2phi*0.5)-w*(t+h*0.5)))));
	
		k3pphi = h*(-1.0)*delta*mu*((r+k2r*0.5)/pow(rprima,3.0))*sin((phi+k2phi*0.5)-w*(t+h*0.5));

	   //Calculamos K4

	   k4r = h*(pr+k3pr);
		k4phi = h*(pphi+k3pphi)/((r+k3r)*(r+k3r));
	
		rprima = pow(1+(r+k3r)*(r+k3r)-2.0*(r+k3r)*cos((phi+k3phi)-w*(t+h)),0.5);
			
		aux = (pphi+k3pphi)*(pphi+k3pphi)/pow((r+k3r),3.0) ;
		aux1 = 1/((r+k3r)*(r+k3r));
		k4pr = h*(aux-delta*(aux1+(mu/pow(rprima,3.0))*((r+k3r)-cos((phi+k3phi)-w*(t+h)))));
	
		k4pphi = h*(-1.0)*delta*mu*((r+k3r)/pow(rprima,3.0))*sin((phi+k3phi)-w*(t+h));

	   //Calculamos los parametros t+h

	    r = r+(1/6.0)*(k1r+2*k2r+2*k3r+k4r);
		phi = phi+(1/6.0)*(k1phi+2*k2phi+2*k3phi+k4phi);
		pr = pr+(1/6.0)*(k1pr+2*k2pr+2*k3pr+k4pr);
		pphi =pphi+(1/6.0)*(k1pphi+2*k2pphi+2*k3pphi+k4pphi);


		//Desescalamos las variables para calcular el hamiltoniano

		praux= pr*m*dtl;

		pphiaux= pphi*m*dtl*dtl;

		raux= r*dtl; 
	


		// Calculamos el hamiltoniano con el cambio de variable

		rprimaH= sqrt(dtl*dtl+raux*raux-2*raux*dtl*cos(phi-w*(t)));


		H =(praux*praux)/(2*m)+(pphiaux*pphiaux)/(2*m*raux*raux)-(G*Mt*m)/(raux)-(G*m*Ml)/(rprimaH)-w*pphiaux;

		cout<<"El valor de H' en el instante "<<i<<" es igual a "<<H<<endl;

		//Lo metemos en un fichero

		f<<t<<" ";
		f<<H<<endl;


		//Aumentamos el tiempo

		t=t+h;

	   //Representamos el movimiento del sistema
	    xTierra = 0;
		yTierra = 0;
		xcohete = r*cos(phi);
		ycohete = r*sin(phi);
		xLuna = cos(w*t);
		yLuna = sin(w*t);

		
		if(j==30){
		fich<<xTierra<<", ";
		fich<<yTierra<<endl;
		fich<<xLuna<<", ";
		fich<<yLuna<<endl;
		fich<<xcohete<<", ";
		fich<<ycohete<<endl;
		fich<<endl;

		j=0;
		}
		else j=j+1;
		
		
		
		

		

    }
    


	fich.close();
	f.close();

    return 0;
}



double CalcularH(double pradio, double ptheta, double ra, double rprima){
	double valor,aux1,aux2,aux3,aux4;

	aux1=(pradio*pradio)/(2*m);
	aux2=(ptheta*ptheta)/(2*m*ra*ra);
	aux3=(G*m*Mt)/(ra);
	aux4=(G*m*Ml)/(rprima);

	valor=aux1+aux2-aux3-aux4;

	return valor;
}