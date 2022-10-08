//Tambor RungeKutta
#include <iostream>
#include <cmath>
using namespace std;

const double Lambda=1.0; 
double f1(double t, double x1, double x2){
  return (-1/t)*x1-Lambda*Lambda*x2;
}

double f2(double t, double x1, double x2){
  return x1;
}

void UnPasoRungeKutta4(double &t0, double &x10, double &x20, double dt){
  double dx11,dx21,dx31,dx41; 
  double dx12,dx22,dx32,dx42;

  dx11=f1(t0,x10,x20)*dt;                     dx12=f2(t0,x10,x20)*dt;
  dx21=f1(t0+dt/2,x10+dx11/2,x20+dx12/2)*dt;  dx22=f2(t0+dt/2,x10+dx11/2,x20+dx12/2)*dt;
  dx31=f1(t0+dt/2,x10+dx21/2,x20+dx22/2)*dt;  dx32=f2(t0+dt/2,x10+dx21/2,x20+dx22/2)*dt;
  dx41=f1(t0+dt,x10+dx31,x20+dx32)*dt;        dx42=f2(t0+dt,x10+dx31,x20+dx32)*dt;
  x10+=(dx11+2*(dx21+dx31)+dx41)/6;           x20+=(dx12+2*(dx22+dx32)+dx42)/6;
  t0+=dt;
  
}

int main(){
  double t,x1,x2;double dt=0.01;

  for(t=0.01,x1=0,x2=1;t<10;  ){
    cout<<t<<" "<<x2<<" "<<x1<<endl;
    UnPasoRungeKutta4(t,x1,x2,dt);
  }
  
  return 0;
}
