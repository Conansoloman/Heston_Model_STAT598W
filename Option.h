
#ifndef Option_h
#define Option_h

#include <iostream>
#include <vector>
#include <random>
#include <math.h>
#include <stack>
#include <queue>

using namespace std;

const int M = 1000; //iteration for Monte Carlo Simulation
const double tick = 0.00001; //time tick

class Option{
    protected:
        double S0, V0, mu, kappa, theta, xi, rho, T;

    public:
    Option(double S0, double V0, double mu, double kappa, double theta, double xi, double rho, double T): S0(S0), V0(V0), mu(mu), kappa(kappa), theta(theta), xi(xi), rho(rho), T(T){
        };

        // Destructor
        ~Option() {
            //Not necessary because new is not used
        };

        virtual double getPrice(){
  
            return 0.0;
        };
};

//Monte Carlo using Heston Process
double milstenScheme(double S0, double V0, double mu, double kappa, double theta, double xi, double rho, double T){
    double sPrice=0.0;
    stack<double> St, Vt;
    St.push(S0);
    Vt.push(V0);
    double Zv,Zs;
    default_random_engine gen1, gen2;
    unsigned seed1 = unsigned(clock());
    unsigned seed2 = unsigned(clock() + 10000);
    gen1.seed(seed1);
    gen2.seed(seed2);
    normal_distribution<double> dist(0.0,1.0);
    //double sum=0.0;
    int NMAX = T/tick;
    for (int j = 0; j < NMAX; j ++){
        //count = count + 1;
        Zv = dist(gen1);
        //cout << Zv << endl;
        Zs = rho*Zv+sqrt(1-rho*rho)*dist(gen2);
        St.push(St.top()*(1+mu*tick+sqrt(Vt.top()*tick)*Zs+0.25*St.top()*tick*(Zs*Zs-1)));
        Vt.push(Vt.top()+kappa*(theta-Vt.top())*tick+xi*sqrt(Vt.top()*tick)*Zv+0.25*xi*xi*tick*(Zv*Zv-1));
    }
    //cout << "St: " << St.top()<<endl;
    sPrice = St.top();
    return sPrice;
}

//Monte Carlo using Advanced QE-Scheme
double QEScheme(double S0, double V0, double mu, double kappa, double theta, double xi, double rho, double T){
    double sPrice;
    queue<double> Vt;
    stack<double> St;
    St.push(S0);
    Vt.push(V0);
    double Zv,Zs;
    double Uv;
    default_random_engine gen1, gen2;
    unsigned seed1 = unsigned(clock());
    unsigned seed2 = unsigned(clock() + 10000);
    gen1.seed(seed1);
    gen2.seed(seed2);
    normal_distribution<double> dist(0.0,1.0);
    uniform_real_distribution<double> unidist(0.0,1.0);
    int NMAX = T/tick;
    double phi_c = 1.5;
    double gamma_1 = 0.5;
    double gamma_2 = 0.5;
    for (int j = 0; j < NMAX; j++){
        double m = theta + (Vt.front()-theta)*exp(-kappa*tick);
        double s2 = ((Vt.front()*xi*xi*exp(-kappa*tick))/kappa)*(1-exp(-kappa*tick))+(theta*xi*xi/(2*kappa))*(1-exp(-kappa*tick))*(1-exp(-kappa*tick));
        double phi = s2/(m*m);
        if (phi <= phi_c) {
            double b2 = 2/phi-1+2*sqrt(2/phi)*sqrt(2/phi-1);
            double b = sqrt(b2);
            double a = m/(1+b2);    // phi <= 2
            Zv = dist(gen1);
            Vt.push(a*(b+Zv)*(b+Zv));
        }
        else if (phi > phi_c){
            double p = (phi-1)/(phi+1);
            double beta = (1-p)/m;
            Uv = unidist(gen1);
            if (Uv <= p)
                Vt.push(0.0);
            else
                Vt.push(log((1-p)/(1-Uv))/beta);
        }
        Zs = dist(gen2);
        double K0 = -rho*kappa*theta*tick/xi;
        double K1 = gamma_1*tick*(kappa*rho/xi-0.5)-rho/xi;
        double K2 = gamma_2*tick*(kappa*rho/xi-0.5)+rho/xi;
        double K3 = gamma_1*tick*(1-rho*rho);
        double K4 = gamma_2*tick*(1-rho*rho);
        double tmp = Vt.front();
        Vt.pop();
        St.push(St.top()+mu*tick+K0+K1*tmp+K2*Vt.front()+sqrt(K3*tmp+K4*Vt.front())*Zs);
    }
    sPrice = St.top();
    return sPrice;
}


class EuropeanCall : public Option{
    //using Option::Option;
    private:
        double K;
    public:
    EuropeanCall(double S0, double V0, double mu, double kappa, double theta, double xi, double rho, double T, double K):Option(S0, V0, mu, kappa, theta, xi, rho, T), K(K){
        };
    double getPrice(){
        stack<double> optionPrice;
        double optionP = 0.0;
        double stockPrice;
        for (int i = 0; i < M; i ++){
            //stockPrice = milstenScheme(S0, V0, mu, kappa, theta, xi, rho, T);
            stockPrice = QEScheme(S0, V0, mu, kappa, theta, xi, rho, T);
            if (stockPrice-K > 0.0)
                optionPrice.push(stockPrice-K);
            else
                optionPrice.push (0.0);
        }
        while (!optionPrice.empty()) {
            optionP = optionP + optionPrice.top();
            optionPrice.pop();
        }
        optionP = optionP/M;
        return optionP;
    };
};

class EuropeanPut : public Option{
    private:
        double K;
    public:
        EuropeanPut(double S0, double V0, double mu, double kappa, double theta, double xi, double rho, double T, double K):Option(S0, V0, mu, kappa, theta, xi, rho, T), K(K){
    };
    double getPrice(){
        stack<double> optionPrice;
        double optionP = 0.0;
        double stockPrice;
        for (int i = 0; i < M; i++){
            stockPrice = QEScheme(S0, V0, mu, kappa, theta, xi, rho, T);
            if (K - stockPrice > 0.0)
                optionPrice.push(K - stockPrice);
            else
                optionPrice.push(0.0);
        }
        while (!optionPrice.empty()) {
            optionP = optionP + optionPrice.top();
            optionPrice.pop();
        }
        optionP = optionP/M;
        return optionP;
    };
};

#endif

