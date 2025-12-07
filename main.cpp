#include "ode_solver.h"

//y(0) = 1
double df(double t, std::vector<double> y){
    return 2*y[0];
}

double f(double t){
    return 1 * exp(2*t);
}

//[1, 0, 0]

double dx1(double t, std::vector<double> y){
    return y[1];
}

double dx2(double t, std::vector<double> y){
    return y[2];
}

double dx3(double t, std::vector<double> y){
    return 6*y[2] - 11*y[1] + 6*y[0];
}

double dx_solution(double t){
    return exp(t)*(-3*exp(t) + exp(2*t) + 3);
}


double dy(double t, std::vector<double> y){
    return sin(y[0]) + pow(t, 1.0 / 10.0);
}


#define sigma 10.0
#define ro 28.0
#define beta 8.0/3.0

double l1(double t, std::vector<double> y){
    return sigma*(y[1] - y[0]);
}

double l2(double t, std::vector<double> y){
    return y[0]*(ro - y[2]) - y[1];
}

double l3(double t, std::vector<double> y){
    return y[0]*y[1] - beta*y[2];
}



int main(void){
    { // y' = 2y
        // function *functions = (function*)malloc(1 * sizeof(function));
        // functions[0] = df;
        // ODE_system system = {functions, 1};
        // std::vector<double> y0 = {1};
        // double start = 0;
        // double stop = 10;
        // double step = 0.1;
        // double tolerance = 1e-6;
        
        // std::vector<ODE_solution> euler_result_2y = ode(start, stop, y0, &system, step, euler);
        // save_data("wynik/euler_result_2y.txt", euler_result_2y);
        
        // std::vector<ODE_solution> mid_point_result_2y = ode(start, stop, y0, &system, step, mid_point);
        // save_data("wynik/mid_point_result_2y.txt", mid_point_result_2y);
        
        // std::vector<ODE_solution> heun_result_2y = ode(start, stop, y0, &system, step, heun);
        // save_data("wynik/heun_result_2y.txt", heun_result_2y);
        
        // std::vector<ODE_solution> rk4_result_2y = ode(start, stop, y0, &system, step, rk4);
        // save_data("wynik/rk4_result_2y.txt", rk4_result_2y);

        // //adaptacyjna metoda 
        
        // auto fehlberg = adaptive_ode(start, stop, start ,y0, &system, step, tolerance, Fehlberg);
        // std::vector<ODE_solution> fehlberg_result_2y = fehlberg.first;
        // save_data("wynik/fehlberg_result_2y.txt", fehlberg_result_2y);
        // std::cout<<"adas\n";


        // auto cash_karp = adaptive_ode(start, stop, start ,y0, &system, step, tolerance, Cash_Karp);
        // std::vector<ODE_solution> cash_karp_result_2y = cash_karp.first;
        // save_data("wynik/cash_karp_result_2y.txt", cash_karp_result_2y);
        // std::cout<<"xddd\n";

        // std::vector<Point> solution_2y = analitic_solution(start, stop, 0.01, f);
        // save_data("wynik/solution_2y.txt", solution_2y);

        
        // fun *sols = (fun*)malloc(1 * sizeof(fun));
        // sols[0] = f;
        // ODE_system_sol system_sol = {sols, 1};
        
        // save_dif("dif_euler.txt", euler_result_2y, system_sol);
        
        // save_dif("dif_mid_point.txt", mid_point_result_2y, system_sol);
        
        // save_dif("dif_heun.txt", heun_result_2y, system_sol);
        
        // save_dif("dif_rk4.txt", rk4_result_2y, system_sol);

        // save_dif("cash_karp_dif.txt", cash_karp_result_2y, system_sol);

        // save_dif("fehlberg_dif.txt", fehlberg_result_2y, system_sol);
    }
    



    { // y''' = ...
        // function *functions = (function*)malloc(3 * sizeof(function));
        // functions[0] = dx1;
        // functions[1] = dx2;
        // functions[2] = dx3;
        // ODE_system system = {functions, 3};
        // std::vector<double> y0 = {1, 0, 0};
        // double start = 0;
        // double stop = 1;
        // double step = 0.1;
        // double tolerance = 1e-6;

        // std::vector<ODE_solution> euler_result_sin = ode(start, stop, y0, &system, step, euler);
        // save_data("wynik/euler_result_sin.txt", euler_result_sin);
        
        // std::vector<ODE_solution> mid_point_result_sin = ode(start, stop, y0, &system, step, mid_point);
        // save_data("wynik/mid_point_result_sin.txt", mid_point_result_sin);
        
        // std::vector<ODE_solution> heun_result_sin = ode(start, stop, y0, &system, step, heun);
        // save_data("wynik/heun_result_sin.txt", heun_result_sin);
        
        // std::vector<ODE_solution> rk4_result_sin = ode(start, stop, y0, &system, step, rk4);
        // save_data("wynik/rk4_result_sin.txt", rk4_result_sin);

        // // Adaptacyjne metody
        // auto cash_karp = adaptive_ode(start, stop, start ,y0, &system, step, tolerance, Cash_Karp);
        // std::vector<ODE_solution> cash_karp_result_sin = cash_karp.first;
        // std::vector<ODE_solution> huj = cash_karp.second;
        // save_data("wynik/cash_karp_result_y'''.txt", cash_karp_result_sin);
        // std::cout<<"xdasda\n";

        // auto fehlberg = adaptive_ode(start, stop, start ,y0, &system, step, tolerance, Fehlberg);
        // std::vector<ODE_solution> fehlberg_result_sin = fehlberg.first;
        // save_data("wynik/fehlberg_result_y'''.txt", fehlberg_result_sin);
        // std::cout<<"dsadasdas\n";


        // std::vector<Point> solution_sin = analitic_solution(start, stop, 0.01, dx_solution);
        // save_data("wynik/solution_sin.txt", solution_sin);

        // fun *sols_sin = (fun*)malloc(1 * sizeof(fun));
        // sols_sin[0] = dx_solution;
        // ODE_system_sol system_sol_sin = {sols_sin, 1};
        
        // save_dif("wynik/dif_dyy_euler.txt", euler_result_sin, system_sol_sin);
        // save_dif("wynik/dif_dyy_mid_point.txt", mid_point_result_sin, system_sol_sin);
        // save_dif("wynik/dif_dyy_heun.txt", heun_result_sin, system_sol_sin);        
        // save_dif("wynik/dif_dyy_rk4.txt", rk4_result_sin, system_sol_sin);
        // save_dif("wynik/dif_dyy_RKF_CK.txt", cash_karp_result_sin, system_sol_sin);
        // save_dif("wynik/dif_dyy_RKF.txt", fehlberg_result_sin, system_sol_sin);
        // free(functions);
    }
    {
        // function *functions = (function*)malloc(1 * sizeof(function));
        // functions[0] = dy;
        // ODE_system system = {functions, 1};
        // std::vector<double> y0 = {1};
        // double start = 0;
        // double stop = 10;
        // double step = 0.5;
        // double tolerance = 1e-6;
        
        // std::vector<ODE_solution> euler_result_2y = ode(start, stop, y0, &system, step, euler);
        // save_data("wynik/euler_result_bez.txt", euler_result_2y);
        
        // std::vector<ODE_solution> mid_point_result_2y = ode(start, stop, y0, &system, step, mid_point);
        // save_data("wynik/mid_point_result_bez.txt", mid_point_result_2y);
        
        // std::vector<ODE_solution> heun_result_2y = ode(start, stop, y0, &system, step, heun);
        // save_data("wynik/heun_result_bez.txt", heun_result_2y);
        
        // std::vector<ODE_solution> rk4_result_2y = ode(start, stop, y0, &system, step, rk4);
        // save_data("wynik/rk4_result_bez.txt", rk4_result_2y);

        // //adaptacyjna metoda 

        // auto cash_karp = adaptive_ode(start, stop, start ,y0, &system, step, tolerance, Cash_Karp);
        // std::vector<ODE_solution> cash_karp_result_2y = cash_karp.first;
        // // std::vector<ODE_solution> cash_karp_result_failed_2y = cash_karp.second;
        // save_data("wynik/cash_karp_result_bez.txt", cash_karp_result_2y);
        // // save_data("cash_karp_result_failed_bez.txt", cash_karp_result_failed_2y);


        // auto fehlberg = adaptive_ode(start, stop, start ,y0, &system, step, tolerance, Fehlberg);
        // std::vector<ODE_solution> fehlberg_result_bez = fehlberg.first;
        // // std::vector<ODE_solution> fehlberg_result_failed_bez = fehlberg.second;
        // save_data("wynik/fehlberg_result_bez.txt", fehlberg_result_bez);
        // // save_data("cash_karp_result_failed_bez.txt", fehlberg_result_failed_bez);
        // free(functions);
    }

        { // lorenz = ...
        function *functions = (function*)malloc(3 * sizeof(function));
        functions[0] = l1;
        functions[1] = l2;
        functions[2] = l3;
        ODE_system system = {functions, 3};
        std::vector<double> y0 = {1, 1, 1};
        double start = 0;
        double stop = 50;
        double step = 0.1;
        double tolerance = 1e-6;

        // std::vector<ODE_solution> euler_lorenz = ode(start, stop, y0, &system, step, euler);
        // save_data("wynik/euler_lorenz.txt", euler_lorenz);
        
        // std::vector<ODE_solution> mid_point_lorenz = ode(start, stop, y0, &system, step, mid_point);
        // save_data("wynik/mid_point_lorenz.txt", mid_point_lorenz);
        
        // std::vector<ODE_solution> heun_lorenz = ode(start, stop, y0, &system, step, heun);
        // save_data("wynik/heun_lorenz.txt", heun_lorenz);
        
        // std::vector<ODE_solution> rk4_lorenz = ode(start, stop, y0, &system, step, rk4);
        // save_data("wynik/rk4_lorenz.txt", rk4_lorenz);

        // Adaptacyjne metody
        // auto cash_karp = adaptive_ode(start, stop, start ,y0, &system, step, tolerance, Cash_Karp);
        // std::vector<ODE_solution> cash_karp_lorenz = cash_karp.first;
        // std::vector<ODE_solution> huj = cash_karp.second;
        // save_data("wynik/cash_karp_lorenz.txt", cash_karp_lorenz);
        // std::cout<<"xdasda\n";

        auto fehlberg = adaptive_ode(start, stop, start ,y0, &system, step, tolerance, Fehlberg);
        std::vector<ODE_solution> fehlberg_lorenz = fehlberg.first;
        save_data("wynik/fehlberg_lorenz.txt", fehlberg_lorenz);
        std::cout<<"dsadasdas\n";        
        free(functions);
    }

    printf("KONIEC\n");
}