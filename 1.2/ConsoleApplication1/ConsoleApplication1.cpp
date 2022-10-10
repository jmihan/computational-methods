#include <iostream>
#include <fstream>
#include <iomanip>
using namespace std;
ifstream input("input.txt", ios::in);
void read(double*, int);
int main()
{
    int n = 0;
    input >> n;
    double** matr = new double*[n];
    for (int i = 0; i < n; i++)
    {
        matr[i] = new double[n];
        for (int j = 0; j < n; j++)
        {
            input >> matr[i][j];
        }
    }
    double* f = new double[n];
    read(f, n);

    double* d = new double[n];
    double* u = new double[n - 1];
    d[0] = matr[0][0];
    u[0] = (1 / d[0]) * matr[0][1];
    for (int i = 1; i < (n - 1); i++) //LU разложение
    {
        d[i] = matr[i][i] - matr[i][i - 1] * u[i - 1];
        u[i] = (1 / d[i]) * matr[i][i + 1];
    }
    d[n - 1] = matr[n - 1][n - 1] - matr[n - 1][n - 2] * u[n - 2];

    double* y = new double[n];
    y[0] = (1 / d[0]) * f[0];
    for (int i = 1; i < (n - 1); i++) //Ly = f
    {
        y[i] = (1 / d[i]) * (f[i] - matr[i][i - 1] * y[i - 1]);
    }
    y[n - 1] = (1 / d[n - 1]) * (f[n - 1] - matr[n - 1][n - 2] * y[n - 2]);

    double* x = new double[n];
    x[n - 1] = y[n - 1];
    for (int i = (n - 2); i >= 0; i--) //Ux = y
    {
        x[i] = y[i] - u[i] * x[i + 1];
    }
    
    for (int i = 0; i < n; i++)
    {
        cout << fixed << setprecision(8) << setw(12) << "X" << i <<" = " << x[i] << endl;
    }
    
    input.close();
}
void read(double* mas, int n)
{
    for (int i = 0; i < n; i++)
    {
        input >> mas[i];
    }
}