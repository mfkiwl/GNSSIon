#include "mainwindow.h"
#include "ui_mainwindow.h"
#include "math.h"
#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <ctime>
#include <cmath>
#include <fstream>
#include <string>

using namespace std;

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    connect(ui->pushButton, SIGNAL(clicked()), this, SLOT(OnPushButtonClicked()));
}

MainWindow::~MainWindow()
{
    delete ui;
}

const double R0 = 6370; 	// Earth's radius, km
const double PI = 3.14159265;
const double SoL = 3;	// Speed of light in 10^8 m/s
double X, Z, PSY, Path, Eps, EpsDiffX, EpsDiffZ, *a, *b, *c, *xk, f0, d_f0, step, TECv;
int n; //количество узлов интерполяции
bool withinTheBorders = true;

void interpolation () {
    string line, set, str = "";
    int numOfRows = 0, j = 0;

    ifstream ReadFile("data.txt");

    while (!ReadFile.eof()) {
        line = "";
        getline (ReadFile, line);
        set += line;
        set += " ";
        numOfRows++;
    }
    ReadFile.close();

    double arr[numOfRows], y[numOfRows + 1];
    n = numOfRows - 1;

    xk = new double [n + 2];
    a = new double [n + 2];
    b = new double [n + 2];
    c = new double [n + 2];

    for (int i = 0; i < set.length(); i++) {
        if (set[i] == ' ') {
            arr[j] = stod(str);
            str = "";
            j++;
            continue;
        }
        str += set[i];
    }
    step = arr[0];

    for (int k = 1; k <= n; k++) {           //значения узлов + расширение
        xk[k] = step * (k - 1);
    }
    xk[0] = xk[1];
    xk[n + 1] = xk[n];

    for (int k = 1; k <= n; k++) {   //значения функции критической частоты в узлах + расширение
        y[k] = arr[k];
    }
    y[0] = (y[2] - y[1]) / (xk[2] - xk[1]);
    y[n + 1] = (y[n] - y[n - 1]) / (xk[n] - xk[n - 1]);

    for (int k = 1; k <= n; k++) {           //значения коэффициентов a[i]
        a[k] = y[k];
    }

    b[1] = y[0];
    for (int k = 2; k <= n; k++) {       //значения коэффициентов b[i] + расширение
        b[k] = 2 * (y[k] - y[k - 1]) / step - b[k - 1];
    }
    b[n + 1] = y[n + 1];


    for (int k = 1; k <= n; k++) {       // значения коэффициентов c[i] + расширение
        c[k] = (b[k + 1] - b[k]) / 2 / step;
    }
}

void criticalFrequency (double x) {
    int interval = x / step + 1;
    if (interval > n) {
        withinTheBorders = false;
    }
    f0 = a[interval] + b[interval] * (x - xk[interval]) + c[interval] * (x - xk[interval]) * (x - xk[interval]);
    d_f0 = b[interval] + 2 * c[interval] * (x - xk[interval]);
}

void verticalTEC (double HasE, double endOfIon, double Z0, double Zm, double Wm, double fr, double Z_01, double S, double h) {
    bool notAssigned = true;
    double z = 0, eps; TECv = 0;

    while (z < HasE) {
        if (z < endOfIon) {
            if (z <= Z0) {
                eps = 1 - exp(-(z - Zm) * (z - Zm) / Wm / Wm) * f0 * f0 / fr / fr;
            }
            else {
                eps = 1 - 0.1 * exp(- (z - Z_01) / S) * f0 * f0 / fr / fr;
            }
        }
        else if (notAssigned) {
                eps = 1;
                notAssigned = false;
        }
        TECv = TECv + h * fr * fr * (1. - eps) / 808.;
        z += h;
    }
}
//******************************//
//  dZ/dX, dPsy/dX, dP/dX		//
//  functions definition		//
//******************************//

double f (double Z, double PSY) {	//dZ/dX
    return (1 + Z / R0)/ tan(PSY);
}

double g (double Z, double PSY) {	//dPSY/dX
    return EpsDiffX / 2 / Eps / tan(PSY) - (1 + Z / R0) * EpsDiffZ / 2 / Eps - 1 / R0;
}

double Q (double Z,double PSY) {     //dP/dX
    return (1 + Z / R0)/ sin(PSY)/pow(Eps, 0.5) / 0.3;
}

void MainWindow :: OnPushButtonClicked() {
    interpolation();
    ui -> textEdit -> clear();
    double fr, Ang1, Ang2, Ang_step, Zm, Wm, Z_01, HasE, Ang0, endOfIon;
    double k1, k2, k3, k4, m1, m2, m3, m4, q1, q2, q3, q4, h;   //Runge-Kutta coefficients
    double C, D, Z0, S, Path = 0, dz = 0, dPSY = 0, dp = 0, bx = 0, PathASE = 0, GPath = 0, TEC, TECt, Delay, DelayT;
    bool notAssigned = true;
    fr = ui -> lineEdit -> text().toDouble();
    Ang1 = ui -> lineEdit_2 -> text().toDouble();
    Ang2 = ui-> lineEdit_3 -> text().toDouble();
    Ang_step = ui -> lineEdit_4 -> text().toDouble();
    h = ui -> lineEdit_12 -> text().toDouble();
    HasE = ui -> lineEdit_13 -> text().toDouble(); // Satellite's altitude
    Zm = ui -> lineEdit_7 -> text().toDouble();
    Wm = ui -> lineEdit_8 -> text().toDouble();
    Z_01 = ui -> lineEdit_9 -> text().toDouble();
    endOfIon = ui -> lineEdit_5 -> text().toDouble();
    double d1;
    bool q = true;
    double B1 = Ang1 / Ang_step, B2 = Ang2 / Ang_step;

    //Calculation of profile exp2 + exp
    C = 2 * Zm * Z_01 - Zm * Zm - Wm * Wm * log(0.1); D = Z_01 * Z_01 - C; Z0 = Z_01 - pow(D, 0.5); S = Wm * Wm / (Z0 - Zm) / 2;

    //******************************//
    //	Calculation of PSY and Z	//
    //******************************//


    for (double B0 = B1; B0 <= B2; B0++) {

        X = 0; Z = 0; PSY = (90 - B0 * Ang_step) * PI / 180.; Ang0 = PSY; Path = 0; TEC = 0; notAssigned = true; q = true;
        ui -> textEdit -> insertPlainText(QString("\n\n\t\t***\tPSY=%1\t***\t\n\n").arg(90 - PSY * 180 / PI));

        while (Z > -h and Z < HasE) {

            if (withinTheBorders == false) {
                ui -> textEdit -> insertPlainText(QString("Выход за пределы интерполяционной сетки!\n\n"));
                break;
            }

            if (Z < endOfIon) {
                criticalFrequency(X);
                if (Z <= Z0) {
                    Eps = 1 - (f0 * f0 / fr / fr) * exp(- (Z - Zm) * (Z - Zm) / Wm / Wm);
                    EpsDiffX = - 2 * (1 / fr / fr) * exp(- (Z - Zm) * (Z - Zm) / Wm / Wm) * d_f0 * f0;
                    EpsDiffZ = 2 * (f0 * f0 / fr / fr) * exp(- (Z - Zm) * (Z - Zm) / Wm / Wm) * ((Z - Zm) / Wm / Wm);
                    //printf("%f %f\n", EpsDiffX, EpsDiffZ);
                }
                else {
                    Eps = 1 - 0.1 * (f0 * f0 / fr / fr) * exp(-((Z - Z_01) / S));
                    EpsDiffX = - (0.2 * exp(-(Z - Z_01) / S) * d_f0 * f0 / fr / fr);
                    EpsDiffZ = 0.1 * exp(-(Z - Z_01) / S) * f0 * f0 / fr / fr / S;
                }
            }
            else if (notAssigned) {
                    Eps = 1;
                    EpsDiffX = 0;
                    EpsDiffZ = 0;
                    notAssigned = false;
            }

            k1 = f(Z, PSY) * h;
            m1 = g(Z, PSY) * h;
            q1 = Q(Z, PSY) * h;
            //printf("\tk1=%5.2f\tm1=%5.2f\t",k1,m1);
            k2 = f(Z + k1 / 2, PSY + m1 / 2) * h;
            m2 = g(Z + k1 / 2, PSY + m1 / 2) * h;
            q2 = Q(Z + k1 / 2, PSY + m1 / 2) * h;
            //printf("\tk2=%5.2f\tm2=%5.2f\t",k2,m2);
            k3 = f(Z + k2 / 2, PSY + m2 / 2) * h;
            m3 = g(Z + k2 / 2, PSY + m2 / 2) * h;
            q3 = Q(Z + k2 / 2, PSY + m2 / 2) * h;
            //printf("\tk3=%5.2f\tm3=%5.2f\t",k3,m3);
            k4 = f(Z + k3, PSY + m3) * h;
            m4 = g(Z + k3, PSY + m3) * h;
            q4 = Q(Z + k3, PSY + m3) * h;
            //printf("\tk4=%5.2f\tm4=%5.2f\n",k4,m4);
            dz = (k1 + 2 * k2 + 2 * k3 + k4) / 6.;
            Z = Z + dz;
            dPSY = (m1 + 2 * m2 + 2 * m3 + m4) / 6.;
            PSY = PSY + dPSY;
            dp = (q1 + 2 * q2 + 2 * q3 + q4) / 6.;
            Path = Path + dp;

            if (Z >= Zm && q) {
                d1 = X;
                verticalTEC(HasE, endOfIon, Z0, Zm, Wm, fr, Z_01, S, h);
                q = false;
            }

            X += h;
            TEC = TEC + dp * fr * fr * (1. - Eps) * 0.3 / 808.;
        }

        if (withinTheBorders == false) {
            withinTheBorders = true;
        }
        else {
            bx = (Z - HasE) * h / dz;//(pow(1+h/dz,0.5))*(Z-HasE)*cos(PI/2-PSY);
            PathASE = Path - (Z - HasE) / 0.3 / cos(PSY - dPSY);//(pow(1+h*h/dz/dz,0.5))*(Z-HasE)/0.3;
            GPath = pow(R0 * R0 + (R0 + HasE) * ((R0 + HasE) - 2 * R0 * cos((X - bx) / R0)), 0.5) / 0.3;
            TECt = TECv * (1 + 16 * pow(0.53 - B0 * Ang_step / 180., 3));
            Delay = 40.3 * 100000 * TEC / fr / fr / 3.;
            DelayT = 40.3 * 100000 * TECt / fr / fr / 3.;

            ui -> textEdit -> insertPlainText(QString("Ang=%1\tXase=%2\tD1=%3\n\n").arg(90 - PSY * 180 / PI, 0, 'f', 2).arg(bx, 0, 'f', 3).arg(d1, 0, 'f', 3));

            ui -> textEdit -> insertPlainText(QString("P_Path=%1 us\tG_Path=%2 us\tDif_Path=%3 ns\n").arg(PathASE, 0, 'f', 3).arg(GPath, 0, 'f', 3).arg((PathASE - GPath) * 1000, 0, 'f', 2));
            ui -> textEdit -> insertPlainText(QString("P_Path=%1 km\tG_Path=%2 km\tDif_Path=%3 m\n\n").arg(PathASE * 0.3, 0,'f',3).arg(GPath * 0.3, 0, 'f', 3).arg((PathASE - GPath) * 300, 0, 'f', 2));

            ui -> textEdit-> insertPlainText(QString("TEC=%1 U \tTECt=%2 U\tDif_TEC=%3 mU\n").arg(TEC, 0, 'f', 3).arg(TECt, 0, 'f', 3).arg((TEC - TECt) * 1000, 0, 'f' ,3));

            ui -> textEdit -> insertPlainText(QString("\nDelay=%1 ns (%2 m)\tDelayT=%3 ns (%4 m)\t(D-DT)=%5 ns (%6 m) (D-DP)=%7 ns (%8 m)\n").arg(Delay, 0, 'f', 2).arg(Delay * 0.3, 0, 'f', 2).arg(DelayT, 0, 'f', 2).arg(DelayT * 0.3, 0, 'f', 2).arg(Delay - DelayT, 0, 'f', 2).arg((Delay - DelayT) * 0.3, 0, 'f', 2).arg((PathASE - GPath) * 1000 - Delay, 0, 'f', 2).arg((PathASE - GPath) * 300 - 0.3 * Delay, 0, 'f', 2));
            // ui->textEdit->insertPlainText(QString("P_Path=%1 km\tG_Path=%2 km\tDif_Path=%3 m\n").arg(PathASE*0.3,0,'f',3).arg(GPath*0.3,0,'f',3).arg((PathASE-GPath)*300,0,'f',3));

            ui -> progressBar -> setValue(100 * (B0 - B1) / (B2 - B1));
        }
    }
    delete [] xk;
    delete [] a;
    delete [] b;
    delete [] c;
}
