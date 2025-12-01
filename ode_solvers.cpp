#include "ode_solver.h"


void pokazPostep(float aktualny, float maksymalny) {
    float procent = (100 * aktualny) / maksymalny;
    std::cout << "\rPostęp: " << std::setw(7) << procent << "%, t:"<< std::setw(7) << aktualny;
    std::cout.flush();
}



std::vector<double> euler(double t, std::vector<double> y, ODE_system *system, double h){
    std::vector<double> result(system->n);
    for(int i = 0; i < system->n; i++){
        result[i] = y[i] + h * system->f[i](t, y); 
    }
    return result;
}

std::vector<double> mid_point(double t, std::vector<double> y, ODE_system *system, double h){
    std::vector<double> result(system->n);
    std::vector<double> y_half = y;    
    double t_half = t + (h/2);

    for(int j = 0; j < y.size(); j++){
        y_half[j] = y[j] + (h/2) * system->f[j](t, y);
    }
    for(int i = 0; i < system->n; i++){
        result[i] = y[i] + h * system->f[i](t_half, y_half); 
    }
    return result;
}


std::vector<double> heun(double t, std::vector<double> y, ODE_system *system, double h){
    std::vector<double> result(system->n);
    std::vector<double> k1(y.size());  
    std::vector<double> k1p(y.size());
    std::vector<double> k2(y.size());    

    for(int i = 0; i < y.size(); i++){
        k1[i] = system->f[i](t, y);
        k1p[i] = y[i] + h * k1[i];
    }
    for(int i = 0; i < y.size(); i++){
        k2[i] = system->f[i](t + h, k1p);
        result[i] = y[i] + (h/2)*(k1[i] + k2[i]);
    }
    return result;
}

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



std::vector<double> rk4(double t, std::vector<double> y, ODE_system *system, double h){
    std::vector<double> result(system->n);
    std::vector<double> k1(y.size()), k2(y.size()), k3(y.size()), k4(y.size()), k_temp(y.size());  

    for(int i = 0; i < y.size(); i++){
        k1[i] = system->f[i](t, y);
    }
    for(int i = 0; i < y.size(); i++){
        k_temp[i] = y[i] + (h/2.0) * k1[i];
    }

    for(int i = 0; i < y.size(); i++){
        k2[i] = system->f[i](t + (h/2.0), k_temp);
    }
    for(int i = 0; i < y.size(); i++){
        k_temp[i] = y[i] + (h/2.0) * k2[i];
    }

    for(int i = 0; i < y.size(); i++){
        k3[i] = system->f[i](t + (h/2.0), k_temp);
        
    }
    for(int i = 0; i < y.size(); i++){
        k_temp[i] = y[i] + h * k3[i];
    }

    for(int i = 0; i < y.size(); i++){
        k4[i] = system->f[i](t + h, k_temp);   
    }
    for(int i = 0; i < y.size(); i++){
        result[i] = y[i] + (h/6.0)*(k1[i] + 2.0*k2[i] + 2.0*k3[i] + k4[i]);
    }

    return result;
}



std::vector<ODE_solution> ode(
    double a, double b, std::vector<double> y0, ODE_system *system , double h,
    Method method
){
    double t = a;
    int n = (b-a)/h;
    std::vector<ODE_solution> values;
    std::vector<double> y = y0;    
    double end_h = 0;
    for(int i = 0; i <= n; i++){
        values.push_back(ODE_solution{y, t});    
        y = method(t, y, system, h);
        t += h;
    }
    return values;
}


double weighted_sum(int stop, std::vector<double> b, std::vector<std::vector<double>>& k_values, int col){
    double sum = 0;
    for(int i = 0; i < stop; i++){
        sum += b[i]*k_values[i][col];
    }
    return sum;
}


std::pair<double, std::vector<double>> Fehlberg(
    double t, std::vector<double> y, ODE_system *system, double h, double tolerance
){
    std::vector<double> a0 = {0, 0, 0, 0, 0};
    std::vector<double> a1 = {1.0/4 , 0, 0, 0, 0};
    std::vector<double> a2 = {3.0/32, 9.0/32 , 0, 0, 0};
    std::vector<double> a3 = {1932.0/2197, -7200.0/2197, 7296.0/2197, 0.0, 0.0};
    std::vector<double> a4 = {439.0/216, -8.0, 3680.0/513, -845.0/4104, 0.0};
    std::vector<double> a5 = {-8.0/27, 2.0, -3544.0/2565, 1859.0/4104, -11.0/40};

    std::vector<double> e = {0, 1.0/4, 3.0/8, 12.0/13, 1.0, 1.0/2};

    std::vector<double> c = {16.0/135, 0.0, 6656.0/12825, 28561.0/56430, -9.0/50, 2.0/55}; // 5 stopien
    std::vector<double> c_star = {25.0/216, 0.0, 1408.0/2565, 2197.0/4104, -1.0/5, 0.0}; // 4 stopien

    std::vector<std::vector<double>> k_values;


    std::vector<double> k1(y.size()), k2(y.size()), k3(y.size()), 
                        k4(y.size()), k5(y.size()), k6(y.size()),
                        k_temp(y.size());

    for(int i = 0; i < y.size(); i++){
        k1[i] = h * system->f[i](t, y);
    }k_values.push_back(k1);
    
    for(int i = 0; i < y.size(); i++){
        k_temp[i] = y[i] + weighted_sum(1, a1, k_values, i);
    }

    for(int i = 0; i < y.size(); i++){
        k2[i] = h * system->f[i](t + e[1]*h, k_temp);
    }k_values.push_back(k2);

    for(int i = 0; i < y.size(); i++){
        k_temp[i] = y[i] + weighted_sum(2, a2, k_values, i);
    }

    for(int i = 0; i < y.size(); i++){
        k3[i] = h * system->f[i](t + e[2]*h, k_temp);
    }k_values.push_back(k3);

    for(int i = 0; i < y.size(); i++){
        k_temp[i] = y[i] + weighted_sum(3, a3, k_values, i);
    }
    
    for(int i = 0; i < y.size(); i++){
        k4[i] = h * system->f[i](t + e[3]*h, k_temp);
    }k_values.push_back(k4);

    for(int i = 0; i < y.size(); i++){
        k_temp[i] = y[i] + weighted_sum(4, a4, k_values, i);
    }

    for(int i = 0; i < y.size(); i++){
        k5[i] = h * system->f[i](t + e[4]*h, k_temp);
    }k_values.push_back(k5);

    for(int i = 0; i < y.size(); i++){
        k_temp[i] = y[i] + weighted_sum(5, a5, k_values, i);
    }

    for(int i = 0; i < y.size(); i++){
        k6[i] = h * system->f[i](t + e[5]*h, k_temp);
    }k_values.push_back(k6);


    std::vector<double> dy(y.size(), 0.0);
    std::vector<double> eh(y.size(), 0.0);

    for(int i = 0; i < y.size(); i++){
        for(int j = 0; j < 6 ; j++){
            dy[i] += c[j] * k_values[j][i];
            double err = (c[j] - c_star[j]) * k_values[j][i];            
            eh[i] += err;
        }
        eh[i] = std::sqrt(eh[i] * eh[i]);
    }
    double error = *std::max_element(eh.begin(), eh.end());
    return std::make_pair(error, dy);
}


std::pair<double, std::vector<double>> Cash_Karp(
    double t, std::vector<double> y, ODE_system *system, double h, double tolerance
){
    std::vector<double> a0 = {0, 0, 0, 0, 0};
    std::vector<double> a1 = {1.0/5, 0, 0, 0, 0};
    std::vector<double> a2 = {3.0/40, 9.0/40, 0, 0, 0};
    std::vector<double> a3 = {3.0/10, -9.0/10, 6.0/5, 0, 0};
    std::vector<double> a4 = {-11.0/54, 5.0/2, -70.0/27, 35.0/27, 0};
    std::vector<double> a5 = {1631.0/55296, 175.0/512, 575.0/13824, 44275.0/110592, 253.0/4096};

    std::vector<double> e = {0, 1.0/5, 3.0/10, 3.0/5, 1.0, 7.0/8};

    std::vector<double> c = {37.0/378, 0.0, 250.0/621, 125.0/594, 0.0, 512.0/1771}; // 5 stopien
    std::vector<double> c_star = {2825.0/27648, 0, 18575.0/48384, 13525.0/55296, 277.0/14336, 1.0/4}; // 4 stopien

    std::vector<std::vector<double>> k_values;


    std::vector<double> k1(y.size()), k2(y.size()), k3(y.size()), 
                        k4(y.size()), k5(y.size()), k6(y.size()),
                        k_temp(y.size());

    for(int i = 0; i < y.size(); i++){
        k1[i] = h * system->f[i](t, y);
    }k_values.push_back(k1);
    
    for(int i = 0; i < y.size(); i++){
        k_temp[i] = y[i] + weighted_sum(1, a1, k_values, i);
    }

    for(int i = 0; i < y.size(); i++){
        k2[i] = h * system->f[i](t + e[1]*h, k_temp);
    }k_values.push_back(k2);

    for(int i = 0; i < y.size(); i++){
        k_temp[i] = y[i] + weighted_sum(2, a2, k_values, i);
    }

    for(int i = 0; i < y.size(); i++){
        k3[i] = h * system->f[i](t + e[2]*h, k_temp);
    }k_values.push_back(k3);

    for(int i = 0; i < y.size(); i++){
        k_temp[i] = y[i] + weighted_sum(3, a3, k_values, i);
    }
    
    for(int i = 0; i < y.size(); i++){
        k4[i] = h * system->f[i](t + e[3]*h, k_temp);
    }k_values.push_back(k4);

    for(int i = 0; i < y.size(); i++){
        k_temp[i] = y[i] + weighted_sum(4, a4, k_values, i);
    }

    for(int i = 0; i < y.size(); i++){
        k5[i] = h * system->f[i](t + e[4]*h, k_temp);
    }k_values.push_back(k5);

    for(int i = 0; i < y.size(); i++){
        k_temp[i] = y[i] + weighted_sum(5, a5, k_values, i);
    }

    for(int i = 0; i < y.size(); i++){
        k6[i] = h * system->f[i](t + e[5]*h, k_temp);
    }k_values.push_back(k6);


    std::vector<double> dy(y.size(), 0.0);
    std::vector<double> eh(y.size(), 0.0);

    for(int i = 0; i < y.size(); i++){
        for(int j = 0; j < 6 ; j++){
            dy[i] += c[j] * k_values[j][i];
            double err = (c[j] - c_star[j]) * k_values[j][i];            
            eh[i] += err;
        }
        eh[i] = std::sqrt(eh[i] * eh[i]);
    }
    double error = *std::max_element(eh.begin(), eh.end());
    return std::make_pair(error, dy);
}


std::pair<std::vector<ODE_solution>, std::vector<ODE_solution>> adaptive_ode(
    double a, double b, double t, std::vector<double> y0, ODE_system *system , double h, double tolerance,
    adaptive_method method
){
    int n = (b-a)/h;
    std::vector<ODE_solution> values;
    std::vector<ODE_solution> rejected_values;
    std::vector<double> y = y0;
    std::vector<double> dy(y.size(), 0);    
    std::pair<double, std::vector<double>> pair = std::make_pair(1, dy);
    // values.push_back(ODE_solution{y, t});   
    bool stop = false;
    double new_h = 0;
    while(true){
        pokazPostep(t, b);
        // std::cout<<"t: "<<t<<" h:"<< h <<"\n";
        pair = method(t, y, system, h, tolerance);
        if(pair.first <= tolerance){
            values.push_back(ODE_solution{y, t});   
            for (int j = 0; j < y.size(); j++){
                y[j] += pair.second[j];
            }
            t += h;
        }else{
            rejected_values.push_back(ODE_solution{y, t}); 
        }
        if(stop){
            break;
        }
        if(pair.first != 0){
            new_h = 0.9 * h * pow((tolerance/pair.first), 0.2);
        }else{
            new_h = h;
        }
        if ((h >= 0.0) == ((t + h) >= b)){
            new_h = b - t;
        }
        if (t >= b) {
            stop = true;
        }
        h = new_h;
    }
    return std::make_pair(values, rejected_values);
}

void print_ode_solution(std::vector<ODE_solution> s){
    for(int i = 0; i < s.size(); i++){
        printf("\n");
        printf("t = %lf, ", s[i].t);
        for (int j = 0; j < s[i].v.size(); j++){
        }
    }
}

std::vector<Point> analitic_solution(double a, double b, double h ,fun f){
    std::vector<Point> values;
    double n = (b-a)/h;
    double tn = a;
    values.push_back({tn, f(tn)});
    for (int i = 0; i < n; i++){
        tn += h;
        values.push_back({tn, f(tn)});
    }
    return values;
}

void save_data(std::string nazwa, std::vector<ODE_solution> s){
    std::ofstream file(nazwa);
    // file << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10);
    file << std::fixed << std::setprecision(6);

    if(!file.is_open()){
        std::cerr << "Nie można otworzyć pliku!" << std::endl;
        return;
    }
    file <<"t; ";
    for(int i = 0; i < s[i].v.size(); i++){
        file<<"y"<<i<<"; ";
    }
    file <<"\n";

    for(int i = 0; i < s.size(); i++){
        file << s[i].t << "; ";
        for (int j = 0; j < s[i].v.size(); j++){
            file << s[i].v[j] << "; ";
        }
        file <<"\n";
    }
    file.close();
}

void save_data(std::string nazwa, std::vector<Point> dane){
    std::ofstream file(nazwa);
    file << std::fixed << std::setprecision(std::numeric_limits<double>::max_digits10);
    if(!file.is_open()){
        std::cerr << "Nie można otworzyć pliku!" << std::endl;
        return;
    }
    file << "t" << "; " << "y0;" << "\n";
    for (const auto& point : dane) {
        file << point.t << "; " << point.y << ";\n";
    }
    file.close();
}


// fun f <- dokładne rozwiązanie
void save_dif(std::string nazwa, const std::vector<ODE_solution>& przyblizone, const ODE_system_sol& dokladne) {
    std::ofstream file(nazwa);
    if (!file.is_open()) {
        std::cerr << "Nie można otworzyć pliku!" << std::endl;
        return;
    }

    file << std::fixed << std::setprecision(6);

    file << "index & t";
    for (int i = 0; i < dokladne.n; ++i) {
        file << "&y" << i << "_approx&y" << i << "_exact&y" << i << "_diff(%)";
    }
    file << "\n";

    for (size_t i = 0; i < przyblizone.size(); ++i) {
        const double t = przyblizone[i].t;
        file << i << " & "<<t;

        for (int j = 0; j < dokladne.n; ++j) {
            double approx = przyblizone[i].v[j];
            double exact = dokladne.f[j](t);
            double diff = std::abs(approx - exact);
            double diff_percent = (std::abs(exact) > 1e-12) ? (std::abs(approx - exact) / std::abs(exact) * 100.0) : 0.0;
            // double diff_percent = (exact != 0.0) ? (std::abs(approx - exact) / std::abs(exact) * 100.0) : 0.0;

            file<<" & " << approx << " & " << exact << " & " << diff_percent << "\\\\";
        }
        file << "\n";
    }
    file.close();
}