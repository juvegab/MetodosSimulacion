#include <iostream>
#include <cmath>
using namespace std;
 
const double Gamma=0.08; //Ahora Gamma es nuestra unica constante global
 
double f1(double t,double s,double i, double Beta){ //Tanto la función 1 como la 2 van a 
	return (-Beta*s*i);                             //depender de Beta, ya que esta será variable
}
double f2(double t,double s,double i,double Beta){ 
	return (Beta*s*i-Gamma*i);
}


void UnPasodeRungeKutta(double & t, double & s0, double & i0, double dt, double & Beta){ // La función Runge Kutta también dependerá de Beta
	double ds1,ds2,ds3,ds4,di1,di2,di3,di4;
	ds1=dt*f1(t,s0,i0,Beta);
	di1=dt*f2(t,s0,i0,Beta);	
	
	ds2=dt*f1(t+dt/2,s0+ds1/2,i0+di1/2,Beta);
	di2=dt*f2(t+dt/2,s0+ds1/2,i0+di1/2,Beta);
	
	ds3=dt*f1(t+dt/2,s0+ds2/2,i0+di2/2,Beta);
	di3=dt*f2(t+dt/2,s0+ds2/2,i0+di2/2,Beta);
	
	ds4=dt*f1(t+dt,s0+ds3,i0+di3,Beta);
	di4=dt*f2(t+dt,s0+ds3,i0+di3,Beta);			
			
	s0+=(ds1+2*(ds2+ds3)+ds4)/6;
	i0+=(di1+2*(di2+di3)+di4)/6; 
	t+=dt;
}



int main(){
	double t, s, i, Beta, R, L, dt=0.1; //Decalarción de variables
		
	
	for(Beta=0.35;Beta<=8;Beta+=0.01){//Primer ciclo for que varia a Beta
		
		s=0.999; i=0.001; //Importante: las condiciones iniciales deben estar dentro del ciclo for de Beta
		L=200; //L será el valor de t cuando tiende a un valor muy grande
		
		for(t=0;t<=L;){//Ciclo for que corre a Runge Kutta hasta el tiempo L
			UnPasodeRungeKutta(t,s,i,dt,Beta);
		}
		
		R=Beta/Gamma; //Valor de R
		cout<<R<<"	"<<s<<endl; //Imprime el valor de R y S_infinito (s(L))
	}	
	return 0;
}