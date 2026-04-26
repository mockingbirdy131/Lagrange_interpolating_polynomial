#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <random>
#include <cmath>
using namespace std;

const double a = 0.0;       
const double b = 6.28;   
const int K = 3;                            // количество подотрезков      
const int N = 6;                            // степень многочленов
const int M = N + (N-1)*(K-1);             // всего узлов на отрезке
const int M_vis = 1000;                  

double func(double x);
void print(double value, int size, int space, ofstream &fout);
double interpol(double **p, double *y, int k, int segment);
void product(double **p, double *x, double step, int m);
int deltas(double *delta, double *absolut, double d);

double func(double x){
    return pow(x, 2);           
}

void print(double value, int size, int space, ofstream &fout){
    fout << scientific << setprecision(size) << setw (space) << right << value;
}

double interpol(double **p, double *y, int k, int segment){
    double sum = 0.0;
    int i, j;
    for (i = 0; i < N; i ++){
        j = i + (N-1)*(segment-1);
        sum += y[j]*p[k][i];
    }
    return sum;
}

void product(double **p, double *x, double step, int m){          
    int i, j, k = 0;
    double tochka = a;
    while (1){
        for (i = 0; i < N; i ++){
            p[k][i] = 1.0;
            for (j = 0; j < N; j ++)
                if (j != i)
                    p[k][i] *= (tochka-x[j])/(x[i]-x[j]);
        }
        tochka += step;
        k ++;
        if (k == m)
            break;
    }
}

int deltas(double *delta, double *absolut, double d){
    int i;
    double otn_l1 = 0.0, otn_l2 = 0.0, otn_oo = absolut[0];
    double abs_l1 = 0.0, abs_l2 = 0.0, abs_oo = delta[0];

    ofstream fout("delta.txt");

    for (i = 0; i < d-1; i ++){           // в b погрешность не смотрим
        abs_l1 += delta[i];
        abs_l2 += delta[i] * delta[i];
        if (delta[i] > abs_oo)
            abs_oo = delta[i];
        otn_l1 += absolut[i];
        otn_l2 += absolut[i] * absolut[i];
        if (absolut[i] > otn_oo)
            otn_oo = absolut[i];
    } 
    abs_l2 = sqrt(abs_l2);
    otn_l2 = abs_l2 / sqrt(otn_l2);
    otn_l1 = abs_l1 / otn_l1;
    otn_oo = abs_oo / otn_oo;
    
    fout << "otn_l1 = "; print(otn_l1, 3, 9, fout); fout << "      abs_l1 = "; print(abs_l1, 3, 9, fout); fout << "\n"; 
    fout << "otn_l2 = "; print(otn_l2, 3, 9, fout); fout << "      abs_l2 = "; print(abs_l2, 3, 9, fout); fout << "\n"; 
    fout << "otn_oo = "; print(otn_oo, 3, 9, fout); fout << "      abs_oo = "; print(abs_oo, 3, 9, fout); fout << "\n"; 

    fout.close();
    return 0;
}

int main(){                        
    int i, j, k, m = M_vis/K;;   
    double x[M] = {0}, y[M] = {0};
    double res, h = (b-a)/(M-1), step = (b-a)/(M_vis-1);
    double **p;
    ofstream fout("data.txt");

    for (i = 0; i < M; i ++)
        x[i] = a + i*h;                                              
    for (i = 0; i < M; i ++)
        y[i] = func(x[i]);                                              
    fout << M << " " << M_vis << " " << N << " " << K << "\n";
    for (i = 0; i < M; i ++){
        print(x[i], 6, 15, fout);
        print(y[i], 6, 15, fout);
        fout << "\n";
    }

    // визуализация -- произведение для b тоже считаем
    p = new double* [m+1];
    for (i = 0; i < m+1; i++)
        p[i] = new double [N];            
    product (p, x, step, m+1);
    // график строим -- доходим в каждом отрезке до m-1 точки, b отдельно рисуем
    k = 1;                           
    for (j = 0; j < M_vis-1; j ++) {
        res = a + j * step;
        if (j % m == 0 && j != 0)
            k ++;
        print(res, 6, 15, fout);
        print(func(res), 6, 15, fout);
        print(interpol(p, y, j % m, k), 6, 15, fout);
        fout << "\n";
    }
    print(res+step, 6, 15, fout);
    print(func(res+step), 6, 15, fout);
    print(interpol(p, y, m, k), 6, 15, fout);
    fout << "\n";

    for (i = 0; i < M_vis/K+1; i++)
        delete[] p[i]; 
    delete[] p;

    // Погрешности -- произведение для b можно не считать
    int d = M+99*(M-1);   // количество точек погрешностей на всем отрезке
    double *delta = new double [d];
    double *absolut = new double [d];
    m = d/K;
    p = new double* [m];
    for (i = 0; i < m; i++)
        p[i] = new double [N];                         
    step = h/100;                                          
    product(p, x, step, m);
    // Считаем погрешности
    k = 1;                   
    for (j = 0; j < d-1; j ++){
        res = a + j * step;
        if (j % m == 0 && j != 0)
            k ++;
        double inter = interpol(p, y, j % m, k);
        double abs = func(res);
        delta[j] = fabs(inter - abs);
        absolut[j] = fabs(abs);
    }
    deltas(delta, absolut, d);

    
    for (i = 0; i < m; i++)
        delete[] p[i]; 
    delete[] p;
    fout.close();
    system("python 2.py");      // python3 !!!!
    return 0;
}


