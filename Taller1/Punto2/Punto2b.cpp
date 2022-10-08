#include <iostream>
#include <cmath>
using namespace std;

double f1(double t, double x1, double x2, double Lam){
  return (-1/t)*x1-Lam*Lam*x2;
};
  
double f2(double t, double x1, double x2){
  return x1;
};

void UnPasoRungeKutta4(double &t0, double &x10, double &x20, double dt,double Lamb){
    double dx11,dx21,dx31,dx41; 
    double dx12,dx22,dx32,dx42;
    
    dx11=f1(t0,x10,x20,Lamb)*dt;                     dx12=f2(t0,x10,x20)*dt;
    dx21=f1(t0+dt/2,x10+dx11/2,x20+dx12/2,Lamb)*dt;  dx22=f2(t0+dt/2,x10+dx11/2,x20+dx12/2)*dt;
    dx31=f1(t0+dt/2,x10+dx21/2,x20+dx22/2,Lamb)*dt;  dx32=f2(t0+dt/2,x10+dx21/2,x20+dx22/2)*dt;
    dx41=f1(t0+dt,x10+dx31,x20+dx32,Lamb)*dt;        dx42=f2(t0+dt,x10+dx31,x20+dx32)*dt;
    x10+=(dx11+2*(dx21+dx31)+dx41)/6;           x20+=(dx12+2*(dx22+dx32)+dx42)/6;
    t0+=dt;
}

int main(){
  double t,x1,x2,Lm;double dt=0.01;
  for(Lm=0.1;Lm<15;Lm+=dt){
    for(t=0.001,x1=0,x2=1;t<1;  ){
      UnPasoRungeKutta4(t,x1,x2,dt,Lm);
    }
    cout<<Lm<<" "<<x2<<" "<<endl;
  }
  return 0;
}