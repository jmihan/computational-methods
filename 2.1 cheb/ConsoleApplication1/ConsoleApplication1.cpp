#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#define Pi 3.1415926535897932384626433832795028841971
using namespace std;

void copy(double*, double*, int);
void proverka(double**, double*, int);
double* mul_b(double**, double*, int);
void mul(double*, int, double);
double Norma(double*, int);

int main()
{
	setlocale(LC_ALL, "RU");
	ifstream in("input.txt", ios::in);
	int size;
	in >> size;
	double** matrix = new double* [size];
	double* b = new double[size];
	double** result = new double* [size];
	for (int i = 0; i < size; i++)
	{
		matrix[i] = new double[size];
		result[i] = new double[size];
		for (int j = 0; j < size; j++)
		{
			in >> matrix[i][j];
		}
	}
	for (int i = 0; i < size; i++)
	{
		in >> b[i];
	}
	cout << "Матрица:" << endl;
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			cout << matrix[i][j] << " ";
		}
		cout << "    " << b[i] << endl;
	}

	//степенной метод нахождения собственных значений
	//НЕОБХОДИМОЕ УСЛОВИЕ A*=A>=0 (это нужно, чтобы собств. значения были неотрицательные и вещественные)
	//находим макс. собств. значение
	double* predX = new double[size];
	double* predXNorm = new double[size];
	double* x = new double[size];
	double* d = new double[size];
	for (int i = 0; i < size; i++)
	{
		predX[i] = 1;
	}
	double eps = 0.000001;
	double max = 0;

	do
	{
		double norma = Norma(predX, size); //считаем норму
		copy(predXNorm, predX, size);
		mul(predXNorm, size, 1.0 / norma);
		x = mul_b(matrix, predXNorm, size);
		for (int i = 0; i < size; i++)
		{
			d[i] = abs(predX[i] - x[i]);
		}
		max = d[0];
		for (int i = 1; i < size; i++)
		{
			if (d[i] > max)
			{
				max = d[i];
			}
		}
		copy(predX, x, size);
	} while (eps < max);
	double lambMax = Norma(x, size);
	cout << "Максимальное собственное значение: " << lambMax << endl;

	//находим мин. собств. знач. (находим макс. с.з. для (||A|| * E - A))
	double beta = 0.0;
	for (int i = 0; i < size; i++)
	{
		if (Norma(matrix[i], size) > beta)
		{
			beta = Norma(matrix[i], size);
		}
	}
	double** E = new double* [size];
	for (int i = 0; i < size; i++)
	{
		E[i] = new double[size];
		for (int j = 0; j < size; j++)
		{
			if (i != j)
			{
				E[i][j] = 0;
			}
			else
			{
				E[i][j] = 1;
			}
			result[i][j] = E[i][j] * beta - matrix[i][j];
		}
	}
	for (int i = 0; i < size; i++)
	{
		predX[i] = 1;
	}
	do
	{
		double norma = Norma(predX, size);
		copy(predXNorm, predX, size);
		mul(predXNorm, size, 1.0 / norma);
		x = mul_b(result, predXNorm, size);
		for (int i = 0; i < size; i++)
		{
			d[i] = abs(predX[i] - x[i]);
		}
		max = d[0];
		for (int i = 1; i < size; i++)
		{
			if (d[i] > max)
			{
				max = d[i];
			}
		}
		copy(predX, x, size);
	} while (eps < max);
	double p = Norma(x, size);
	double lambMin = beta - p;
	cout << "Минимальное собственное значение: " << lambMin << endl;

	//метод Ричардсона с Чебышевским набором параметрами
	for (int i = 0; i < size; i++)
	{
		predX[i] = 1;
	}
	double tauk;
	int n = 1000000;
	int k = 0;
	double* xa;
	do
	{
		k++;
		tauk = 2.0 / ((lambMax + lambMin) + (lambMax - lambMin) * cos((2.0 * k - 1) * Pi / (2.0 * n)));
		max = 0.0;
		xa = mul_b(matrix, predX, size);
		for (int i = 0; i < size; i++)
		{
			xa[i] -= b[i];
			xa[i] *= tauk;
		}
		for (int i = 0; i < size; i++)
		{
			x[i] = predX[i] - xa[i];
			d[i] = abs(x[i] - predX[i]);
			if (d[i] > max)
			{
				max = d[i];
			}
		}
		copy(predX, x, size);
	} while (max > eps);
	cout << endl << "Решение при помощи метода Ричардсона с Чебышевским набором параметрами:" << endl;
	for (int i = 0; i < size; i++)
	{
		cout << x[i] << endl;
	}
	cout << "Количество итераций: " << k;
	proverka(matrix, x, size);
	return 0;
}
void copy(double* a, double* b, int size)
{
	for (int i = 0; i < size; i++)
	{
		a[i] = b[i];
	}
}
void proverka(double** matrix, double* b, int size)
{
	double* c = new double[size];
	for (int i = 0; i < size; i++)
	{
		c[i] = 0;
	}

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			c[i] += matrix[i][j] * b[j];
		}
	}
	cout << endl << "VERIFICATION:" << endl;
	for (int i = 0; i < size; i++)
	{
		cout << c[i] << endl;
	}
}
double* mul_b(double** matrix, double* b, int size)
{
	double* c = new double[size];
	for (int i = 0; i < size; i++)
	{
		c[i] = 0;
	}
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			c[i] += matrix[i][j] * b[j];
		}
	}
	return c;
}
void mul(double* a, int size, double alpha)
{
	for (int j = 0; j < size; j++)
	{
		a[j] *= alpha;
	}
}
double Norma(double* x, int size)
{
	double res = 0;
	for (int i = 0; i < size; i++)
	{
		res += abs(x[i]);
	}
	//res = sqrt(res);
	return res;
}