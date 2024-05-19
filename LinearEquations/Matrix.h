#pragma once
#include <vector>

class Matrix
{
public:
	// Konstruktory
	Matrix(int N, int a1, int a2, int a3);
	Matrix(const Matrix& A);

	// Destruktor (dealokacja zarezerwowanej pamiêci)
	~Matrix();

	// Operator mno¿enia (poprawiony: memory leak fix)
	std::vector<double> operator*(const std::vector<double>& v);

//private:
	// Pola
	double** tab = nullptr;
	int size = 0;
};