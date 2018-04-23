// Numeric_lesson_integrating.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <ctime>
using namespace std;

double intagratingByTrapeze(double x0, double xk, double dx);
double intagratingByTrapezeArr(double* ArrX, double* ArrY,int n);
double integratingMonteCarlo(double x0, double xk, double y0, double yk, int n);
double integratingUpgradedMonteCarlo(double x0, double xk, int n);
double f(double x);


int main()
{
	double x0 = 0, xk = 10;
	double dx = 0.001;
	double ArrX[] = { 0,2,3,5 };
	double ArrY[] = {0,4,9,25 };
	int n = 4;
	cout << "pole f kwadratowej: " << intagratingByTrapeze(x0, xk, dx) << endl;
	//cout << "pole f kwadratowej: " << intagratingByTrapezeArr(ArrX,ArrY,n);
	cout << "pole f kwadratowej: " << integratingMonteCarlo(x0, xk, -2, f(xk), 1000000) << endl;
	cout << "pole f kwadratowej: " << integratingUpgradedMonteCarlo(x0,xk,1000000)<< endl;
	
	system("pause");
    return 0;
}

double f(double x) {
	return x*x-2;
}

double intagratingByTrapezeArr(double* ArrX, double* ArrY, int n){
	double field = 0;
	for (int i = 1; i < n;i++) {
		double dx = ArrX[i] - ArrX[i - 1];
		field += ((ArrY[i] + ArrY[i-1])*dx) / 2.0;
	}
	return field;
}

double intagratingByTrapeze(double x0, double xk, double dx) {
	double xp = 0;
	double xn = x0;
	double field = 0;
	for (double i = 0; i <= xk; i += dx) {
		xp = xn;
		xn = xp + dx;
		field += ((f(xn) + f(xp))*dx) / 2.0;
	}
	return field;
}

double integratingMonteCarlo(double x0, double xk, double y0, double yk, int n) {
	double rectField = abs(xk - x0) * abs(yk - y0);
	double px, py;
	int c = 0;
	srand(time(NULL));
	for (int i = 0; i < n; i++) {
		px = ((double)rand() / (RAND_MAX)) * abs(xk - x0) + x0;
		py = ((double)rand() / (RAND_MAX)) * abs(yk - y0) + y0;

		if (py > 0 && py <= f(px)) {
			c++;
		}
		else if (py < 0 && f(px) < py) {
			c--;
		}
	}
	double field = rectField*(c / (double)n);
	return field;
}
double integratingUpgradedMonteCarlo(double x0, double xk, int n) {
	double px, py;
	int c = 0;
	srand(time(NULL));
	for (int i = 0; i < n; i++) {
		px = ((double)rand() / (RAND_MAX)) * abs(xk - x0) + x0;
		c += f(px);
	}
	double field = (c / (double)n)*abs(xk - x0);
	return field;
}

