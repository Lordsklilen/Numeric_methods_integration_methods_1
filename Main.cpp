// Numeric_lesson_integrating.cpp : Defines the entry point for the console application.
//

#include <iostream>
#include <ctime>
using namespace std;

double intagratingByTrapeze(double x0, double xk, double dx);
double intagratingByTrapezeArr(double* ArrX, double* ArrY, int n);
double integratingMonteCarlo(double x0, double xk, double y0, double yk, int n);
double integratingUpgradedMonteCarlo(double x0, double xk, int n);
double f(double x);
double ** ilorazy(double *X, double * Y, int n);
double LagrangePolynomial(double* xs, double* ys, int n, double x);
double intagratingByTrapezeCotesNewton(double x0, double xk, int n, double* xs, double* ys, int size);
double intagratingByMonteCarloCotesNewton(double x0, double xk, int n, double xs[], double ys[], int size);
double Simson(double x0,double xk,int n);
double Simson38(double xp, double xk, int n);
int main()
{
	/*double x0 = 0, xk = 10;
	double dx = 0.001;
	double ArrX[] = {0,2,3,5 };
	double ArrY[] = {0,4,9,25 };
	int n = 4;
	cout << "pole f kwadratowej: " << intagratingByTrapeze(x0, xk, dx) << endl;
	cout << "pole f kwadratowej: " << integratingMonteCarlo(x0, xk, -2, f(xk), 1000000) << endl;
	cout << "pole f kwadratowej: " << integratingUpgradedMonteCarlo(x0,xk,1000000)<< endl;*/

	double xs[] = { -4, -2, 0, 1,3 };
	double ys[] = { 16, 4, 0, 1, 9 };
	int n = 5;
	//double ** dif = ilorazy(X, Y, n);
	//cout << LagrangePolynomial( xs,  ys,  n,2);
	cout << intagratingByTrapezeCotesNewton(0, 4, 100, xs, ys, n) << endl;
	cout << intagratingByMonteCarloCotesNewton(0, 4, 100, xs, ys, n) << endl;
	cout << "Simson3/8: " << Simson38(0, 10, 1000000) << endl;
	cout << "Simson: " << Simson(0,10,1000000) << endl;

	system("pause");
	return 0;
}

double Simson38(double xp,double xk,int n) {
	
	double dx = (xk - xp) / (double)n;
	double sum = f(xp) + f(xk);
	for (int i = 1; i < n; i++) {
		if (i % 3 == 0)
			sum += 2 * (f(xp + i * dx));
		else 
			sum += 3 * (f(xp + i * dx));
	}
	return 3.0*dx/8.0 * sum;
}
double Simson(double xp, double xk, int n) {
	double sum = f(xp)+f(xk);
	double h = (xk - xp) / (double)n;
	double tmp = 0;
	for (int i = 0; i < n; i++) {
		double x = xp + i * h;
		tmp += f(x-h/2);
		sum += f(x);
	}
	tmp += f(xk - h / 2);
	sum = (h / 6.0)*(f(xp) + f(xk) + 2 * sum + 4 * tmp);
	return sum;

}


double f(double x) {
	return x * x - 2;
}
double** ilorazy(double *X, double * Y, int n) {
	double **dif = new double*[n];
	for (int i = 0; i < n; i++) {
		dif[i] = new double[n];
		dif[0][i] = Y[i];
	}
	for (int i = 1; i < n; i++) {
		for (int j = i; j < n; j++) {
			dif[i][j] = (dif[i - 1][j] - dif[i - 1][j - 1]) / (X[j] - X[j - i]);
		}
	}
	return dif;
}

double LagrangePolynomial(double* xs, double* ys, int n, double x) {
	double wynik = 0;
	for (int i = 0; i < n; i++) {
		double tmp = 1;
		for (int j = 0; j < n; j++) {
			if (j != i)
				tmp *= (x - xs[j]) / (xs[i] - xs[j]);
		}
		wynik += tmp * ys[i];
	}
	return wynik;
}
double intagratingByTrapezeCotesNewton(double x0, double xk, int n, double xs[], double ys[], int size) {
	double xp = 0;
	double xn = x0;
	double field = 0;
	double dx = (xk - xn) / (double)n;
	for (double i = 1; i < n; i++) {
		xp = x0 + (i - 1.)*dx;;
		xn = x0 + (i)*dx;
		field += ((LagrangePolynomial(xs, ys, size, xn) + LagrangePolynomial(xs, ys, size, xp))*dx) / 2.0;

	}
	return field;
}
double intagratingByMonteCarloCotesNewton(double x0, double xk, int n, double xs[], double ys[], int size) {
	double px, py;
	int c = 0;
	srand(time(NULL));
	for (int i = 0; i < n; i++) {
		px = (((double)rand() / (RAND_MAX)) * abs(xk - x0)) + x0;
		c += LagrangePolynomial(xs, ys, size, px);
		/*double t= LagrangePolynomial(xs, ys, size, px);
		cout << px << "," << t << endl;
		c += t;*/
	}
	double field = (c / (double)n)*abs(xk - x0);
	return field;
}

double intagratingByTrapezeArr(double* ArrX, double* ArrY, int n) {
	double field = 0;
	for (int i = 1; i < n; i++) {
		double dx = ArrX[i] - ArrX[i - 1];
		field += ((ArrY[i] + ArrY[i - 1])*dx) / 2.0;
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
	double field = rectField * (c / (double)n);
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




