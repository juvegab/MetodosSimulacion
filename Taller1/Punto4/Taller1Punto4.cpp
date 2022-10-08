#include <iostream>
#include <cmath>
#include "../vector.h"
using namespace std;

//Constantes globales

const int N=3;
const double g=980;
const double ord=1e10;
const double K=10.0;

//const double K=0.2;
//const double K=0.5;
//const double K=1.0;
//const double K=2.5;
//const double K=10.0;
//constantes de PEFRL
const double Zeta=0.1786178958448091e00;
const double Lambda=-0.2123418310626054e0;
const double Chi=-0.6626458266981849e-1;
const double Coeficiente1=(1-2*Lambda)/2;
const double Coeficiente2=1-2*(Chi+Zeta);

//Declaración de las clases
class Cuerpo;
class Colisionador;

//---------- Clase Cuerpo --------------
class Cuerpo{
private:
  double theta, omega, tau,m,R,l,x0,I;
  
public:
  void Inicie(double theta0,double omega0, double m0,
	     double R0,double l0,double x00);
  void BorreTorque(void){tau=0;};
  void SumeTorque(double tau0){tau+=tau0;};
  void Mueva_theta(double dt,double coeficiente);
  void Mueva_omega(double dt,double coeficiente);
  void Dibujese(void);
  double Getx(void){return x0+l*sin(theta);}; //Inline
  double Gety(void){return -l*cos(theta);}; //Inline
  double Gettau(void){return tau;}; //Inline
  friend class Colisionador;
};
void Cuerpo::Inicie(double theta0,double omega0, double m0,double R0,double l0,double x00){
  theta=theta0; omega=omega0; m=m0; R=R0;l=l0; x0=x00;I=m*l*l;
}
void Cuerpo::Mueva_theta(double dt,double coeficiente){
  theta+=omega*(dt*coeficiente);
}
void Cuerpo::Mueva_omega(double dt,double coeficiente){
  omega+=tau*(dt*coeficiente/I);
}
void Cuerpo::Dibujese(void){
  cout<<" , "<<Getx()<<"+"<<R<<"*cos(t),"<<Gety()<<"+"<<R<<"*sin(t)";
    cout<<" , "<<x0<<"+"<<l/7<<"*t*sin("<<theta<<"),"<<"-"<<l/7<<"*t*cos("<<theta<<")";
}
//---------- Clase Colisionador --------------
class Colisionador{
private:
public:
  void CalculeTorques(Cuerpo * Pendulo);
  void CalculeTorqueEntre(Cuerpo & Pendulo1,Cuerpo & Pendulo2);    
};
void Colisionador::CalculeTorques(Cuerpo * Pendulo){
  int i,j;double tau0;
  //Borrar fuerzas
  for(i=0;i<N;i++){
  Pendulo[i].BorreTorque();
  tau0=-Pendulo[i].m*Pendulo[i].l*g*sin(Pendulo[i].theta);
  Pendulo[i].SumeTorque(tau0);
  }
  
  //Calcular las fuerzas entre todas las parejas de Pendulos
   for(i=N-1;i>0;i--)
     
     CalculeTorqueEntre(Pendulo[i],Pendulo[i-1]);
}
void Colisionador::CalculeTorqueEntre(Cuerpo & Pendulo1,Cuerpo & Pendulo2){
  double s=Pendulo2.Getx()+Pendulo2.R-(Pendulo1.Getx()-Pendulo1.R);double F=0;
  if(s>0) F=K*ord*pow(s,1.5);
  Pendulo1.SumeTorque(F*Pendulo1.l);
  Pendulo2.SumeTorque(-F*Pendulo2.l);

}

//----------- Funciones Globales -----------

void InicieAnimacion(void){
  //  cout<<"set terminal gif animate"<<endl; 
  //  cout<<"set output 'DosPendulos.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange[-5:23]"<<endl;
  cout<<"set yrange[-15:2]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
    cout<<"plot 0,0 ";
}
void TermineCuadro(void){
    cout<<endl;
}


int main(){
  Cuerpo Pendulo[N];
  Colisionador Newton;
  double m0=100, l0=12, R0=1.5;
  double T=2*M_PI*sqrt(l0/g);
  double t,tmax=1.1*T,dt=0.00001;
  
  double tdibujo,tcuadro=T/100;
  int i;
  
  //---------------(theta0, omega0, m0, R0, l0,x00)
  Pendulo[0].Inicie(-0.261799,0,m0,R0,l0, 0);
  for(i=1;i<N;i++){Pendulo[i].Inicie(-0,0,m0,R0,l0, 2*R0*i);}
  
  //InicieAnimacion();
  
  for(t=0,tdibujo=0; t<tmax; t+=dt,tdibujo+=dt){
    //Dibujar
  //  if(tdibujo>tcuadro){
    //  InicieCuadro();
    //  for(i=0;i<N;i++) Pendulo[i].Dibujese();
     // TermineCuadro();
     // tdibujo=0;
    //}         
    
    if(t>0.174){if(t<0.178){cout<<t<<" "<<Pendulo[1].Gettau()<<endl;};};
    // Mover por PEFRL
    //for(i=0;i<N;i++)
    for(i=0;i<N;i++)  Pendulo[i].Mueva_theta(dt,Zeta);
    Newton.CalculeTorques(Pendulo);
    for(i=0;i<N;i++) Pendulo[i].Mueva_omega(dt,Coeficiente1);
    for(i=0;i<N;i++) Pendulo[i].Mueva_theta(dt,Chi);
    Newton.CalculeTorques(Pendulo);
    for(i=0;i<N;i++) Pendulo[i].Mueva_omega(dt,Lambda);
    for(i=0;i<N;i++) Pendulo[i].Mueva_theta(dt,Coeficiente2);
    Newton.CalculeTorques(Pendulo);
    for(i=0;i<N;i++) Pendulo[i].Mueva_omega(dt,Lambda);
    for(i=0;i<N;i++) Pendulo[i].Mueva_theta(dt,Chi);
    Newton.CalculeTorques(Pendulo);
    for(i=0;i<N;i++) Pendulo[i].Mueva_omega(dt,Coeficiente1);
    for(i=0;i<N;i++) Pendulo[i].Mueva_theta(dt,Zeta);   
  }
  
  return 0;
}
