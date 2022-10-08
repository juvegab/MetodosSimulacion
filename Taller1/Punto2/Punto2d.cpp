//Ceros por Biseccion de la Funcion de Bessel
#include <iostream>
#include <cmath>
using namespace std;

const double ErrMax=1e-7;

int main(){
  double x,t,dt=0.01, N=10;
  double a=13, b=15, m, fa, fm;

  /*for(t=0;t<N;t+=dt){
    cout<<" "<<t<<" "<<cyl_bessel_j(0,t)<<endl;
    }*/
  
  fa=cyl_bessel_j(0,a);
  while(b-a >= ErrMax){
    m = (b+a)/2; fm=cyl_bessel_j(0,m);
    if(fa*fm>0)
      {a=m; fa=fm;}
    else
      b=m;
  }
  cout<<"El cero es "<<(a+b)/2<<endl;
  return 0;
}


