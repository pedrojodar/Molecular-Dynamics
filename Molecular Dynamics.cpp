// Molecular Dynamics.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <time.h>
using namespace std;
double random() {
    return rand() / (RAND_MAX * 1.0);
}
class particle {
    public:
        double x, y;
        double length_x, length_y;
        double ax, ay;
        double deltat;
        double commom, distx, disty;
        double dist2_inv, dist6_inv, dist12_inv;
        particle() {
            
        }

        void initialize(double x_0, double y_0, double deltat, double length_x_val, double length_y_val) {
            x = x_0;
            y = y_0;
            this->deltat = deltat;
            this->length_x = length_x_val;
            this->length_y = length_y_val;
        }
        void evolve(){
            x = x + deltat * ax;
            y = y + deltat * ay;
            x = x - floor(x / length_x) * length_x;
            y = y - floor(y / length_y) * length_y;

        }
        void new_step() {
            evolve();
            ax = 0;
            ay = 0;

        }
        //Falta poner las condiciones de contorno
        vector<double> calculate_a(double r2_x, double r2_y) {
            distx = x - r2_x;  
            disty = y - r2_y;
            if (distx > length_x) {
                distx = distx - floor(distx / length_x) * length_x;
            }
            if (disty > length_y) {
                distx = distx - floor(distx / length_y) * length_y;

            }
            double dist_2 = distx * distx + disty * disty;
            if (dist_2 < pow(2, 1 / 3.0)) {
                dist2_inv = 1 / (dist_2);
                dist6_inv = dist2_inv * dist2_inv * dist2_inv;
                dist12_inv = dist6_inv * dist6_inv;
                commom = 4 * (12 * dist12_inv - 6 * dist6_inv) * (dist2_inv);
                ax = commom * distx;
                ay = commom * disty;
                return {commom * distx, commom*disty};
            }
            else {
                return { 0.0, 0.0 };
            }
        }
        void plot_msd() {

        }
};
class active_particle : public particle {
    public:
        double amplitude, omega;
        void force_magnitude(double ax_val, double omega_val){
            amplitude = ax_val;
            omega = omega_val;
        }
        void active_force(double iter) {
            ax = ax + amplitude * sin(omega * iter * deltat);
        }
};
class systemage {
    public:
        int iter;
        vector <particle> pasive_particles;
        vector <active_particle> active_particles;
        systemage(int number_particles, double deltat, double lenx, double leny, double frac_pass, double omega, double force_magnitude ) {
            
            pasive_particles.resize((int )frac_pass* number_particles);
            active_particles.resize(number_particles- pasive_particles.size());
            int total_particles = pasive_particles.size() + active_particles.size();
            int part_x = floor (sqrt(total_particles)+1);
            int counter = 0;
            iter = 0;
            for (int i = 0 ; i < part_x ; i ++){
                for (int j = 0; j < part_x; j++){
                    double x_0 = lenx * (i + 1) /(1.0* (part_x + 2));
                    double y_0 = leny * (j + 1) /(1.0 * (part_x + 2));
                    if (counter >= total_particles-1) break;
                    if (counter < pasive_particles.size()) {
                        pasive_particles[counter].initialize(x_0, y_0, deltat, lenx, leny);
                    }
                    else {
                        active_particles[counter-pasive_particles.size()].initialize(x_0, y_0, deltat, lenx, leny);
                        active_particles[counter - pasive_particles.size()].force_magnitude(force_magnitude, omega);
                    }
                    counter++;
                }
            }
        }
        void evolve() {
            vector <double> aux; 
            aux.resize(2);
            for (int i = 0; i < pasive_particles.size(); i++) {
                for (int j = i+1; j < pasive_particles.size(); j++) {
                    aux = pasive_particles[i].calculate_a(pasive_particles[j].x, pasive_particles[j].y);
                    pasive_particles[j].ax = pasive_particles[j].ax - aux[0];
                    pasive_particles[j].ay = pasive_particles[j].ay - aux[1];

                }
                for (int j = 0; j < active_particles.size(); j++) {
                    aux = pasive_particles[i].calculate_a(active_particles[j].x, active_particles[j].y);
                    active_particles[j].ax = active_particles[j].ax - aux[0];
                    active_particles[j].ay = active_particles[j].ay - aux[1];
                }
            }
            for (int i = 0; i < active_particles.size(); i++) {
                for (int j = i + 1; j < active_particles.size(); j++) {
                    aux = active_particles[i].calculate_a(active_particles[j].x, active_particles[j].y);
                    active_particles[j].ax = active_particles[j].ax - aux[0];
                    active_particles[j].ay = active_particles[j].ay - aux[1];
                }
                active_particles[i].active_force(iter);
            }            
            #pragma omp parallel for
            for (int i = 0; i < pasive_particles.size(); i++) {
                pasive_particles[i].new_step();
            }
            #pragma omp parallel for 
            for (int i = 0; i < active_particles.size(); i++) {
                active_particles[i].new_step();
            }
            iter++;
        }
        void check() {

        }
        void print_system() {
            ofstream fich; 
            fich.open("image");
            for (int i = 0; i < pasive_particles.size(); i++) {
                fich << pasive_particles[i].x << "\t" << pasive_particles[i].y << endl;
            }
            fich << endl << endl;
            fich.close();
            
        }

};
int main()
{
    systemage s(10, 0.0001, 20, 20, 0.5, 1, 10);
    for (int i = 0; i < 1000; i++) {
        cout << i << endl;
        s.evolve();
        s.print_system();
    }

}
