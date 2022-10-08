#include <iostream>
#include <cmath>
using namespace std;

 const double Beta=0.35; //Declaración constantes globales
 const double Gamma=0.08; 
 
 //Escribimos las ecuaciones acopladas como funciones
double f1(double t,double s, double i){ //ds/dt
	return (-Beta*s*i);//Lo que nos devuelve la función 1 Susceptibles
}
double f2(double t,double s, double i){ //di/dt
	return (Beta*s*i-Gamma*i);//Lo que nos devuelve la función 2 Infectados
}

//Programa Runge Kutta
void UnPasodeRungeKutta(double & t, double &s0,double &i0, double &dt){ 
	double ds1,ds2,ds3,ds4,di1,di2,di3,di4;
	ds1=dt*f1(t,s0,i0);
	di1=dt*f2(t,s0,i0);	
	
	ds2=dt*f1(t+dt/2,s0+ds1/2,i0+di1/2);
	di2=dt*f2(t+dt/2,s0+ds1/2,i0+di1/2);
	
	ds3=dt*f1(t+dt/2,s0+ds2/2,i0+di2/2);
	di3=dt*f2(t+dt/2,s0+ds2/2,i0+di2/2);
	
	ds4=dt*f1(t+dt,s0+ds3,i0+di3);
	di4=dt*f2(t+dt,s0+ds3,i0+di3);			
			
	s0+=(ds1+2*(ds2+ds3)+ds4)/6;
	i0+=(di1+2*(di2+di3)+di4)/6; 
	t+=dt;
}


//Función principal
int main(){
	
	double t, s, i, dt=0.1; //Declaración de variables
	s=0.999; i=0.001; //Condiciones iniciales 
	
	for(t=0;t<=80;){//ciclo for dentro de un tiempo de 80 dias
		cout<<t<<"	"<<s<<"	"<<i<<"	"<<-(s+i)+1<<endl;//imprime los valores de t, s(t), i(t) y r(y)
		UnPasodeRungeKutta(t,s,i,dt);
	}
	 	
	return 0;
}