#include <iostream>
#include<cmath>
#include<iomanip>
#include <math.h>
using namespace std;

int main(){
    double beta, gamma, b, d, a, v, I_0, N, S_0, t, E_0, R_0;
    beta=0.2;
    gamma=0.1;
    b=pow(10,-4);
    d=pow(10,-4);
    a=0.1;
    v=pow(10,-3);
    I_0=pow(10,-6);
    N=1;
    S_0=N-I_0;
    t=0.01;
    E_0=0;
    R_0=0;
    //cout<<beta<< gamma<< b<< d<< a<< v<< I_0<< N<< S_0<< t<< E_0<< R_0;
    double h, l, c, x1, x2, f, g;
    f=d*t+gamma*t+1;
    g=1+d*t+a*t;
    h=(gamma*v*pow(t,3)*a)/(1+d*t+v*t)-f*g*t*beta;
    l=b*pow(t,2)*a+t*a*S_0+(R_0*pow(t,2)*v*a)/(1+d*t+v*t)-f*g-f*g*d*t-I_0*g*t*beta;
    c=I_0*(g+g*d*t);
    double x = (l*l) - (4*h*c); 
    
    if (x <= 0){
        x = x*(-1);
        cout<<"Solución solo en números complejos"<<endl;
        cout<<"Solución en numeros complejos: " <<(-l/(2*h))<<" + "<<(sqrt(x)/(2*h))<<"i y "<<(-l/(2*h))<<" - "<<(sqrt(x)/(2*h))<<"i"<<endl;
    }else{
        x1 = (double)((-l + sqrt(x)) / (2*h));
        x2 = (double)((-l - sqrt(x)) / (2*h));
        
        cout<<"x1 = "<<x1<<endl;
        cout<<"x2 = "<<x2<<endl;
    }
    double I_1=x2;
    double R_1=(R_0+I_1*gamma*t)/(1+d*t+v*t);
    double S_1=(b*t+S_0+R_1*v*t)/(1+d*t+I_1*t*beta);
    double E_1=(E_0+I_1*t*S_1*beta)/(1+d*t+a*t);
    cout<< "I = "<< I_1<<" R = " <<R_1<<" S = "<< S_1<<" E =" << E_1<<"\n";

    for(int i=0; i<=1000; i++){
    f=d*t+gamma*t+1;
    g=1+d*t+a*t;
    h=(gamma*v*pow(t,3)*a)/(1+d*t+v*t)-f*g*t*beta;
    l=b*pow(t,2)*a+t*a*S_0+(R_0*pow(t,2)*v*a)/(1+d*t+v*t)-f*g-f*g*d*t-I_0*g*t*beta;
    c=I_0*(g+g*d*t);
    double x = (l*l) - (4*h*c);
    
    if (x <= 0){
        x = x*(-1);
        cout<<"Solución solo en números complejos"<<endl;
        cout<<"Solución en numeros complejos: " <<(-l/(2*h))<<" + "<<(sqrt(x)/(2*h))<<"i y "<<(-l/(2*h))<<" - "<<(sqrt(x)/(2*h))<<"i"<<endl;
    }else{
        x1 = (double)((-l + sqrt(x)) / (2*h));
        x2 = (double)((-l - sqrt(x)) / (2*h));
        
        cout<<"x1 = "<<x1<<endl;
        cout<<"x2 = "<<x2<<endl;
    }
    double I_1=x2;
    double R_1=(R_0+I_1*gamma*t)/(1+d*t+v*t);
    double S_1=(b*t+S_0+R_1*v*t)/(1+d*t+I_1*t*beta);
    double E_1=(E_0+I_1*t*S_1*beta)/(1+d*t+a*t);
    cout<< "I = "<< I_1<<" R = " <<R_1<<" S = "<< S_1<<" E =" << E_1<<"\n";
    I_0=I_1;
    R_0=R_1;
    S_0=S_1;
    E_0=E_1;
    t=t+0.01;
    }
}