#include <iostream>
#include <vector>
#include <math.h>
#include <fstream>

#include <iomanip>
#include <chrono>
#include <thread>

#include <utility>
#include <algorithm> // std::min_element


typedef struct{
    double t;
    double y;
}Point;

typedef double (*fun)(double);

typedef double (*function)(double t, std::vector<double> y);


typedef struct {
    function *f;
    int n;
}ODE_system;



typedef struct{
    fun *f;
    int n;
}ODE_system_sol;


typedef std::vector<double> (*Method)(double t, std::vector<double> y, ODE_system *system , double h);

typedef std::pair<double, std::vector<double>>(*adaptive_method)(
    double t, std::vector<double> y, ODE_system *system , double h, double tolernace);

typedef struct{
    std::vector<double> v;
    double t;
}ODE_solution;

std::vector<double> euler(double t, std::vector<double> y, ODE_system *system, double h);
std::vector<double> mid_point(double t, std::vector<double> y, ODE_system *system, double h);
std::vector<double> heun(double t, std::vector<double> y, ODE_system *system, double h);

// std::vector<double> heun_nowy(double t, std::vector<double> y, ODE_system *system, double h){
//     std::vector<double> result(system->n);
//     std::vector<double> k1(y.size());  
//     std::vector<double> k1p(y.size());
//     std::vector<double> k2(y.size());    

//     for(int i = 0; i < y.size(); i++){
//         k1[i] = system->f[i](t, y);
//     }
//     for(int i = 0; i < y.size(); i++){
//         k1p[i] = y[i] + h * k1[i];
//     }
//     for(int i = 0; i < y.size(); i++){
//         k2[i] = system->f[i](t + h, k1p);
//     }
//     for(int i = 0; i < y.size(); i++){
//         result[i] = y[i] + (h/2)*(k1[i] + k2[i]);
//     }
    
//     return result;
// }

std::vector<double> rk4(double t, std::vector<double> y, ODE_system *system, double h);

std::vector<ODE_solution> ode(
    double a, double b, std::vector<double> y0, ODE_system *system , double h,
    Method method
);

double weighted_sum(int stop, std::vector<double> b, std::vector<std::vector<double>>& k_values, int col);

std::pair<double, std::vector<double>> Fehlberg(
    double t, std::vector<double> y, ODE_system *system, double h, double tolerance
);

std::pair<double, std::vector<double>> Cash_Karp(
    double t, std::vector<double> y, ODE_system *system, double h, double tolerance
);

std::pair<std::vector<ODE_solution>, std::vector<ODE_solution>> adaptive_ode(
    double a, double b, double t, std::vector<double> y0, ODE_system *system , double h, double tolerance,
    adaptive_method method
);

void print_ode_solution(std::vector<ODE_solution> s);
std::vector<Point> analitic_solution(double a, double b, double h ,fun f);
void save_data(std::string nazwa, std::vector<ODE_solution> s);
void save_data(std::string nazwa, std::vector<Point> dane);
// fun f <- dokładne rozwiązanie
void save_dif(std::string nazwa, const std::vector<ODE_solution>& przyblizone, const ODE_system_sol& dokladne);

void pokazPostep(float aktualny, float maksymalny);