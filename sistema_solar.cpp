//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////// CALCULO DEL SISTEMA SOLAR CON ALGORITMO DE VERLET //////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////


//Librerias implicadas en el codigo
#include <iostream>
#include <math.h>
#include <fstream>
#include <ostream>
#include <string>
#include <stdlib.h>
#include <conio.h>
using namespace std;

//Declaracion de las funciones usadas
double RescalarPosicion(double);
double RescalarMasa(double);
double RescalarTiempo(double);
double DesescalarTiempo(double);
double CalcularModulo(double,double,double,double);
double CalcularAceleracion(double, double, double,double,double);
double CalcularPotencial(double, double, double, double, double, double);



int main(){

    double r_x[9], r_y[9], m[9], v_x[9], v_y[9], a_x[9], a_y[9], w_x[9], w_y[9],aux_y[9];
    int i,j, k;
    double h, t, Ep, E[9], valor, E_total,P[9];

    //Inicializamos los valores de  a_x y a_y a 0
    for (i = 0; i <9; i++)
    {
        a_x[i]=0;
        a_y[i]=0;
        r_y[i]=0;
        r_x[i]=0;
        v_x[i]=0;
        v_y[i]=0;
        w_x[i]=0;
        w_y[i]=0;
        m[i]=0;
        E[i]=0;
        P[i]=0;
    }

    //Rellenamos los vectores r_x,v_y y masa con los ficheros

    ifstream fich;

    fich.open("masas.dat");

    if (fich.is_open() == true)
    {
        for ( i = 0; i < 9; i++)
        {
            fich >> m[i];
            m[i]=m[i]*pow(10,24);
            m[i]=RescalarMasa(m[i]);
            
        }
        fich.close();
        m[0]=1;
        
    }
    else cout<<"Error al abrir el fichero" <<endl;



    fich.open("radio.dat");
    
    if (fich.is_open() == true)
    {
        for ( i = 0; i < 9; i++)
        {
            fich >> r_x[i];
            r_x[i]=r_x[i]*pow(10,9);
            r_x[i]=RescalarPosicion(r_x[i]);
            
            
        }
        fich.close();
        r_x[0]=0;
    }
    else cout<<"Error al abrir el fichero" <<endl;


        fich.open("velocidades.dat");
    
    if (fich.is_open() == true)
    {
        for ( i = 0; i < 9; i++)
        {
            fich >> v_y[i];            

            //RESCALAMOS LA VELOCIDAD
            v_y[i]=v_y[i]*1000;
            v_y[i]= RescalarPosicion(v_y[i]);
            v_y[i]= RescalarTiempo(v_y[i]);
            
            
        
        }
        fich.close();
        v_y[0]=0;
    }
    else cout<<"Error al abrir el fichero" <<endl;


    //Calculamos la aceleracion inicial
    for ( i = 0; i < 9; i++)
    {
        for ( j = 0; j < 9; j++)
        {
            if (j!=i)
            {
                

                a_x[i]+=CalcularAceleracion(m[j],r_x[i],r_x[j],r_y[i],r_y[j]);
                a_y[i]+=CalcularAceleracion(m[j],r_y[i],r_y[j],r_x[i],r_x[j]);
                
            }
            
        }    
    }
    



    //INICIAMOS LAS K ITERACIONES REALIZADAS EN EL ALGORITMO DE VERLET
    h=0.05;
    t=0;

    ofstream fichero,f,dat;

    fichero.open("planets_data.dat");
    f.open("Valores_Energia.dat");
    dat.open("Periodos.dat");

    for ( k = 0; k < 10000; k++)
    {
        cout<<"Intento numero "<<k<<endl;
        //Imprimimos los valores en un fichero ofstream "planets_data"

        if (fichero.is_open() == true)
        {
            for ( i = 0; i < 9; i++)
            {
                fichero <<r_x[i]<< ", ";
            
                fichero<<r_y[i]<<endl;
            

            }
            fichero<<endl;
        }
        else cout<<"Error al abrir el fichero" <<endl;

        
        //Energia total del sistema 
        for(i = 0; i < 9; i++){

	        for (j = 0; j < 9; j++){
                
                Ep=0;
                if (i!=j){

		 	        Ep+=CalcularPotencial(m[i],m[j],r_x[i],r_x[j],r_y[i],r_y[j]);

		        }
		 	
		    }
            E[i]=-Ep*0.5;
            
	    }


        //Calculamos la posicion r(t+h) con la aceleracion a(t)

        for ( i = 0; i < 9; i++)
        {
            aux_y[i]=r_y[i]; //Almaceno la posicion y anterior (util para calculo del periodo)
            r_x[i]+= h*v_x[i]+(0.5*h*h)*a_x[i];
            r_y[i]+= h*v_y[i]+(0.5*h*h)*a_y[i];

        }
        
        //Pasamos la informacion de a(t) a un vector auxiliar w(i)

        for ( i = 0; i < 9; i++)
        {
            w_x[i]=v_x[i]+(0.5*h)*a_x[i];
            w_y[i]=v_y[i]+(0.5*h)*a_y[i];
            
        }
        


        //Calculamos la nueva aceleracion a(t+h) con la posicion r(t+h)
        for ( i = 0; i < 9; i++)
        {
            a_x[i]=0;
            a_y[i]=0;
            for ( j = 0; j < 9; j++)
            {
                if (j!=i)
                {

                    a_x[i]+=CalcularAceleracion(m[j],r_x[i],r_x[j],r_y[i],r_y[j]);
                    a_y[i]+=CalcularAceleracion(m[j],r_y[i],r_y[j],r_x[i],r_x[j]);
                
                }
            
            }

        }
        

        //Calculamos la nueva velocidad v(t+h) con a(t+h) e w(i)

        for ( i = 0; i < 9; i++)
        {
            v_x[i]=w_x[i]+(h*0.5)*a_x[i];
            v_y[i]=w_y[i]+(h*0.5)*a_y[i];
        }


        //Calculo del periodo de orbita de los planetas
        for(i = 1; i < 9; i++){
		    //cout<<"La y antigua es: "<<aux_y[i]<<endl;
            //cout<<"La y nueva es: "<<r_y[i]<<endl;
            if ( aux_y[i]<0 && r_y[i]>0 && P[i]==0){

			    P[i] =t;
                dat<<P[i]<<endl;
                cout<<"El periodo del planeta "<<i<<" es: "<<P[i]<<endl;

		    }
	    }
        
        //Evaluamos el tiempo
        t+=h;
        

        if (f.is_open() == true)
        {
            for ( i = 0; i < 9; i++)
            {

                f << E[i] << ", " ;
            

            }

            f<<endl;
        }
        else cout<<"Error al abrir el fichero" <<endl;

    }
    
    fichero.close();
    f.close();
    dat.close();
    


    return 0;
}






/////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////FUNCIONES USADAS///////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////

//Función que reescala de metros a unidades astronómicas
double RescalarPosicion(double r_i){
    double rescalado,c;
    c=1.496*pow(10,11);

    rescalado=r_i/c;

    return rescalado;
}



//Función que rescalar la masa de los planetas a masas solares
double RescalarMasa(double m_i){
    double rescalar,M_s;
    M_s=1.99*pow(10,30);

    rescalar=m_i/M_s;

    return rescalar;

}



//Función que reescala el tiempo a unidades sistema solar
double RescalarTiempo(double t_i){
    double rescalar,G,c,M_s,cte;
    c=1.496*pow(10,11);
    G=6.67*pow(10,-11);
    M_s=1.99*pow(10,30);
    cte=sqrt((G*M_s)/(c*c*c));

    rescalar=t_i/cte;

    return rescalar;
}





//Calculo del modulo de dos vectores r_i y r_j
double CalcularModulo(double rx_i,double rx_j, double ry_i, double ry_j){
    double mod;

    mod=sqrt((rx_i-rx_j)*(rx_i-rx_j)+(ry_i-ry_j)*(ry_i-ry_j));

     return mod;
}

//Calculo de la aceleracion en un eje para una componente
double CalcularAceleracion(double m,double rx_i,double rx_j,double ry_i, double ry_j){
    double aceleracion;

    

    aceleracion=(-1.0)*m*(rx_i-rx_j)*pow(pow(rx_i-rx_j,2)+pow(ry_i-ry_j,2),-1.5);

    return aceleracion;
}


double CalcularPotencial(double mi, double mj, double rx_i, double rx_j, double ry_i, double ry_j){
    double potencial,modulo;

    modulo = sqrt(  (rx_i-rx_j)*(rx_i-rx_j)+ (ry_i-ry_j)*(ry_i-ry_j));

    potencial= (mi*mj)/modulo;

    return potencial;
}