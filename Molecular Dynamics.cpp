// Molecular Dynamics.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <time.h>
#include <random>
using namespace std;
#define PI 3.141592653589793238
double random() {
    return rand() / (RAND_MAX * 1.0);
}
double gaussian() {
    double u1, u2, S;
    do {
        u1 = random();
        u2 = random();
        u1 = 2 * u1 - 1;
        u2 = 2 * u2 - 1;
        S = u1 * u1 + u2 * u2;
    } while (S >= 1);
    return sqrt(-2 * log(S)/S) *u1;
}
void checkgaussian() {
    ofstream fich;
    fich.open("gaussian_test.txt");
    for (int i = 0; i < 500; i++) {
        fich << gaussian() << endl;
    }
    fich.close();
    system("python3 check_gaussian.py");
}
class particle {
    private:
        double x0, y0;
    public:
        double x, y;

        double length_x, length_y;
        double ax, ay;
        double deltat;
        double Temp;
        double commom, distx, disty;
        double dist2_inv, dist6_inv, dist12_inv;
        particle() {
            
        }
        void position_thermalized() {
            x0 = x;
            y0 = y;
        }
        void initialize(double x_0, double y_0, double deltat, double length_x_val, double length_y_val, double Temp) {
            this->x0 = x_0;
            this->y0 = y_0;
            x = x_0;
            y = y_0;
            this->deltat = deltat;
            this->length_x = length_x_val;
            this->length_y = length_y_val;
            this->Temp = Temp;
            ax = 0;
            ay = 0;
        }
        void evolve(){
            x = x + deltat * ax + sqrt( 2*Temp * deltat)  *gaussian();
            y = y + deltat * ay + sqrt( 2*Temp * deltat) * gaussian();
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
                cout << distx << "\t" << disty << endl;
                cout << dist_2 << endl;

                cout << ax << "\t" << ay << endl;
                return { ax, ay};
            }
            else {
                return { 0.0, 0.0 };
            }
        }
         double msd() {
             return  (x - x0) * (x - x0) + (y - y0) * (y - y0);
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
        systemage(int number_particles, double deltat, double lenx, double leny, double frac_pass, double omega, double force_magnitude, double temp ) {
            
            pasive_particles.resize((int )(frac_pass* number_particles));
            active_particles.resize(number_particles- pasive_particles.size());
            cout << pasive_particles.size() << "\t" << active_particles.size() << endl;
            int total_particles = pasive_particles.size() + active_particles.size();
            int part_x = floor (sqrt(total_particles)+1);
            cout << part_x << endl;
            int counter = 0;
            iter = 0;
            for (int i = 0 ; i < part_x ; i ++){
                for (int j = 0; j < part_x; j++){
                    double x_0 = lenx * (i + 1) /(1.0* (part_x)+0.5);
                    double y_0 = leny * (j + 1) /(1.0 * (part_x)+0.5);
                    cout << x_0 << "\t" << y_0 << endl;
                    if (counter < pasive_particles.size()) {
                        pasive_particles[counter].initialize(x_0, y_0, deltat, lenx, leny, temp);
                    }
                    else {
                        
                        active_particles[counter-pasive_particles.size()].initialize(x_0, y_0, deltat, lenx, leny, temp);
                        active_particles[counter - pasive_particles.size()].force_magnitude(force_magnitude, omega);
                    }
                    if (counter >= total_particles - 1) break;
                    counter++;
                }
            }
            /*
            if (number_particles != pasive_particles.size())
                active_particles[active_particles.size() - 1].initialize(52.7, 5, deltat, lenx, leny, temp);
            else
                pasive_particles[pasive_particles.size() - 1].initialize(52.7, 5, deltat, lenx, leny, temp);
            */
            for (int i = 0; i < active_particles.size(); i++) {
                int index = random() * number_particles;
                    if (index < pasive_particles.size()) {
                        double aux_x = pasive_particles[index].x;
                        double aux_y = pasive_particles[index].y;
                        pasive_particles[index].x = active_particles[i].x;
                        pasive_particles[index].y = active_particles[i].y;
                        active_particles[i].x=aux_x;
                        active_particles[i].y= aux_y;

                    }
            }
            

            ofstream fich;
            fich.open("position.txt");
            for (int i = 0; i < pasive_particles.size(); i++) {
                fich << pasive_particles[i].x << "\t" << pasive_particles[i].y << endl;
            }

            fich.close();
            fich.open("position_ac.txt");
            for (int i = 0; i < active_particles.size(); i++) {
                fich << active_particles[i].x << "\t" << active_particles[i].y << endl;
            }
            system("Python3 initial.py");
        }
        void evolve() {
            vector <double> aux; 
            aux.resize(2);
            for (int i = 0; i < pasive_particles.size(); i++) {
                for (int j = i+1; j < pasive_particles.size(); j++) {
                    aux = pasive_particles[i].calculate_a(pasive_particles[j].x, pasive_particles[j].y);
                    pasive_particles[i].ax = pasive_particles[i].ax + aux[0];
                    pasive_particles[i].ay = pasive_particles[i].ay + aux[1];
                    pasive_particles[j].ax = pasive_particles[j].ax - aux[0];
                    pasive_particles[j].ay = pasive_particles[j].ay - aux[1];
                    

                }
                for (int j = 0; j < active_particles.size(); j++) {
                    aux = pasive_particles[i].calculate_a(active_particles[j].x, active_particles[j].y);
                    pasive_particles[i].ax = pasive_particles[i].ax + aux[0];
                    pasive_particles[i].ay = pasive_particles[i].ay + aux[1];
                    active_particles[j].ax = active_particles[j].ax - aux[0];
                    active_particles[j].ay = active_particles[j].ay - aux[1];
                }
            }
            for (int i = 0; i < active_particles.size(); i++) {
                for (int j = i + 1; j < active_particles.size(); j++) {
                    aux = active_particles[i].calculate_a(active_particles[j].x, active_particles[j].y);
                    active_particles[i].ax = active_particles[i].ax + aux[0];
                    active_particles[i].ay = active_particles[i].ay + aux[1];
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
            vector<double> mean_msq;
            int n_steps = 10;
            int evolution = 100; 
            int n_steps_thermalized = 1000;
            for (int j = 0 ; j < n_steps_thermalized; j++) evolve();
            for (int j = 0; j < pasive_particles.size(); j++) pasive_particles[j].position_thermalized();
            for (int j = 0 ; j < n_steps; j ++){
                //Termati
                for (int j = 0; j < evolution; j++) evolve();
                double mean_msq_aux=0;
                for (int i = 0; i < pasive_particles.size(); i++) {
                    mean_msq_aux += pasive_particles[i].msd();
                }
                mean_msq_aux = mean_msq_aux/(1.0*pasive_particles.size()* pasive_particles[0].deltat*evolution);
                mean_msq.push_back(mean_msq_aux);
            }
            cout << "Finito calculo" << endl;
            ofstream fich;
            fich.open("msd.txt");
            for (int i = 0; i < mean_msq.size(); i++) {
                fich <<i+1<<"\t" <<mean_msq[i] << endl;
            }
            fich.close();
            system("python3 plotmsd.py");
        }
        void print_system(int i) {
            ofstream fich; 
            string str = "video2\\" + to_string(i) + ".txt";
            fich.open(str);
            for (int i = 0; i < pasive_particles.size(); i++) {
                fich << pasive_particles[i].x << "\t" << pasive_particles[i].y << endl;
            }
            fich << endl << endl;
            fich.close();
            str = "video2\\" + to_string(i) + "ac.txt";
            fich.open(str);
            for (int i = 0; i < active_particles.size(); i++) {
                fich << active_particles[i].x << "\t" << active_particles[i].y << endl;
            }
            fich << endl << endl;
            fich.close();
        }

};
int main()
{
    checkgaussian();
    srand(time(NULL));
    systemage s(500, 0.01, 63.24, 63.24, 0.8, 15, 5,0.05);
    //s.check();
    
    for (int i = 0; i < 10000; i++) {
        cout << i << endl;
        s.evolve();
        if (i%100==0)  s.print_system(i/200);
    }
    if (s.active_particles.size()!=0)
        system("python video_ac.py 50");
    else
        system("python video.py 50");

    return 0;
}
