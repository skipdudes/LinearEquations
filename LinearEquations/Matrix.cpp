#include "Matrix.h"
#include <cmath>

Matrix::Matrix(int N, int a1, int a2, int a3)
{
	tab = new double* [N];
	size = N;

	for (int i = 0; i < N; i++)
		tab[i] = new double[N] {};

	// Macierz pasmowa
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (i == j)
				tab[i][j] = a1;
			else if (j == i - 1 || j == i + 1)
				tab[i][j] = a2;
			else if (j == i - 2 || j == i + 2)
				tab[i][j] = a3;
		}
	}
}

Matrix::Matrix(const Matrix& A)
{
	tab = new double* [A.size];
	size = A.size;

	for (int i = 0; i < A.size; i++)
		tab[i] = new double[A.size] {};

	// Przepisanie wartoœci
	for (int i = 0; i < A.size; i++)
		for (int j = 0; j < A.size; j++)
			tab[i][j] = A.tab[i][j];
}

Matrix::~Matrix()
{
	for (int i = 0; i < size; i++)
		delete[] tab[i];
	delete[] tab;
}

std::vector<double> Matrix::operator*(const std::vector<double>& v)
{
	std::vector<double> u(size, 0);

	for (int i = 0; i < size; i++)
	{
		double sum = 0;
		for (int j = 0; j < size; j++)
			sum += tab[i][j] * v[j];
		u[i] = sum;
	}

	return u;
}