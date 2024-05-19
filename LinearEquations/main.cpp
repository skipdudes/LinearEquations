#include <iostream>
#include <string>
#include <fstream>
#include <cmath>
#include <chrono>
#include "Matrix.h"
using namespace std;

// c = 2, d = 5, e = 4, f = 0
// N = 925
// a1 = 5 + 4 = 9
// a2 = a3 = -1
// b[n] = sin(n * (0 + 1)) = sin(n)

// Sta³e
const int N = 925;
const string outputName = "output.txt";

// Norma euklidesowa wektora residuum
double norm(int size, double* v)
{
	double norm = 0;

	for (int i = 0; i < size; i++)
		norm += v[i] * v[i];

	return norm;
}

// Metoda Jacobiego
void Jacobi(const int size, Matrix* A, std::vector<double>& x, std::vector<double>& b, double resi_norm, int& iterations, int& time, double& norm_final)
{
	std::vector<double> r(size);
	std::vector<double> x_tmp(size, 1);

	int i_counter = 0;
	double resi = 0;
	double sum = 0;

	chrono::steady_clock::time_point t_start = chrono::steady_clock::now();

	do
	{
		for (int i = 0; i < size; i++)
		{
			sum = 0;

			for (int j = 0; j < size; j++)
				if (j != i)
					sum += A->tab[i][j] * x_tmp[j];

			x[i] = (b[i] - sum) / A->tab[i][i];
		}

		x_tmp = x;

		// Residuum
		r = *A * x;
		for (int i = 0; i < size; i++)
			r[i] -= b[i];
		resi = norm(size, r.data());

		i_counter++;
	} while (resi > resi_norm);

	chrono::steady_clock::time_point t_end = chrono::steady_clock::now();
	time = (int)chrono::duration_cast<chrono::milliseconds>(t_end - t_start).count();
	iterations = i_counter;
	norm_final = resi;
}

// Metoda Gaussa-Seidla
void Gauss_Seidl(const int size, Matrix* A, std::vector<double>& x, std::vector<double>& b, double resi_norm, int& iterations, int& time, double& norm_final)
{
	std::vector<double> r(size);
	std::vector<double> x_tmp(size, 1);

	int i_counter = 0;
	double resi = 0;
	double sum = 0;

	chrono::steady_clock::time_point t_start = chrono::steady_clock::now();

	do
	{
		for (int i = 0; i < size; i++)
		{
			sum = 0;

			for (int j = 0; j < i; j++)
				sum += A->tab[i][j] * x[j];

			for (int j = i + 1; j < size; j++)
				sum += A->tab[i][j] * x_tmp[j];

			x[i] = (b[i] - sum) / A->tab[i][i];
		}

		x_tmp = x;

		// Residuum
		r = *A * x;
		for (int i = 0; i < size; i++)
			r[i] -= b[i];
		resi = norm(size, r.data());

		i_counter++;
	} while (resi > resi_norm);

	chrono::steady_clock::time_point t_end = chrono::steady_clock::now();
	time = (int)chrono::duration_cast<chrono::milliseconds>(t_end - t_start).count();
	iterations = i_counter;
	norm_final = resi;
}

// Metoda faktoryzacji LU
void LU(const int size, Matrix* A, std::vector<double>& x, std::vector<double>& b, int& time, double& norm_final)
{
	// Rozbicie macierzy
	Matrix* U = new Matrix(*A);
	Matrix* L = new Matrix(size, 1, 0, 0);

	chrono::steady_clock::time_point t_start = chrono::steady_clock::now();

	// A = L * U
	for (int i = 0; i < size - 1; i++)
	{
		for (int j = i + 1; j < size; j++)
		{
			L->tab[j][i] = U->tab[j][i] / U->tab[i][i];

			for (int k = i; k < size; k++)
			{
				U->tab[j][k] = U->tab[j][k] - L->tab[j][i] * U->tab[i][k];
			}
		}
	}

	// L * y = b
	std::vector<double> y(size, 0);
	for (int i = 0; i < size; i++)
	{
		double sum = 0;
		for (int j = 0; j < i; j++)
			sum += L->tab[i][j] * y[j];
		y[i] = (b[i] - sum) / L->tab[i][i];
	}

	// U * x = y
	for (int i = size - 1; i >= 0; i--)
	{
		double sum = 0;
		for (int j = i + 1; j < size; j++)
			sum += U->tab[i][j] * x[j];
		x[i] = (y[i] - sum) / U->tab[i][i];
	}

	chrono::steady_clock::time_point t_end = chrono::steady_clock::now();
	time = (int)chrono::duration_cast<chrono::milliseconds>(t_end - t_start).count();

	// Residuum
	std::vector<double> r = *A * x;
	for (int i = 0; i < size; i++)
		r[i] -= b[i];
	norm_final = norm(size, r.data());

	delete U;
	delete L;
}

int main()
{
	// Zadanie A
	int a1 = 9;
	int a2 = -1;
	int a3 = -1;

	Matrix* A = new Matrix(N, a1, a2, a3);
	std::vector<double> x(N, 0);
	std::vector<double> b(N);

	for (int i = 0; i < N; i++)
		b[i] = sin(i);

	// Zadanie B
	int iterations_j = 0, iterations_gs = 0;
	int time_j = 0, time_gs = 0;
	double norm_j = 0, norm_gs = 0;
	double resi_norm = 1e-9;

	Jacobi(N, A, x, b, resi_norm, iterations_j, time_j, norm_j);
	Gauss_Seidl(N, A, x, b, resi_norm, iterations_gs, time_gs, norm_gs);
	cout << "Zadanie B" << endl;
	cout << "Metoda Jacobiego: " << iterations_j << " iteracji, " << time_j << " ms, norma " << norm_j << endl;
	cout << "Metoda Gaussa-Seidla: " << iterations_gs << " iteracji, " << time_gs << " ms, norma " << norm_gs << endl << endl;

	// Zadanie C
	a1 = 3;
	delete A;
	A = new Matrix(N, a1, a2, a3);

	Jacobi(N, A, x, b, resi_norm, iterations_j, time_j, norm_j);
	Gauss_Seidl(N, A, x, b, resi_norm, iterations_gs, time_gs, norm_gs);
	cout << "Zadanie C" << endl;
	cout << "Metoda Jacobiego: " << iterations_j << " iteracji, " << time_j << " ms, norma " << norm_j << endl;
	cout << "Metoda Gaussa-Seidla: " << iterations_gs << " iteracji, " << time_gs << " ms, norma " << norm_gs << endl << endl;

	// Zadanie D
	int time_lu = 0;
	double norm_lu = 0;

	LU(N, A, x, b, time_lu, norm_lu);
	cout << "Zadanie D" << endl;
	cout << "Metoda faktoryzacji LU: " << time_lu << " ms, norma " << norm_lu << endl << endl;

	// Zadanie E
	a1 = 9;
	int N_values[] = { 100,500,1000,2000,3000,5000,7500 };

	// Plik z wynikami
	ofstream output;
	output.open(outputName.c_str());
	output << "n;jacobi;gauss-seidl;faktoryzacja-lu" << endl;

	for (int i = 0; i < 7; i++)
	{
		delete A;

		A = new Matrix(N_values[i], a1, a2, a3);
		x.resize(N_values[i], 0);
		b.resize(N_values[i]);

		for (int i = 0; i < N_values[i]; i++)
			b[i] = sin(i);

		Jacobi(N_values[i], A, x, b, resi_norm, iterations_j, time_j, norm_j);
		Gauss_Seidl(N_values[i], A, x, b, resi_norm, iterations_gs, time_gs, norm_gs);
		LU(N_values[i], A, x, b, time_lu, norm_lu);

		cout << "Zadanie E (N = " << N_values[i] << ")" << endl;
		cout << "Metoda Jacobiego: " << time_j << " ms" << endl;
		cout << "Metoda Gaussa-Seidla: " << time_gs << " ms" << endl;
		cout << "Metoda faktoryzacji LU: " << time_lu << " ms" << endl << endl;

		// Zapis wyników
		output << N_values[i] << ';' << time_j << ';' << time_gs << ';' << time_lu << endl;
	}

	output.close();

	// Dealokacja
	delete A;

	return 0;
}