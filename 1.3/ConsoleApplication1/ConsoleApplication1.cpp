#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;
ifstream in("input.txt", ios::in);
int main()
{
    int n = 0;
    in >> n;
    double** source = new double* [n];
    double** matr = new double* [n];
    for (int i = 0; i < n; i++)
    {
        source[i] = new double[n + 1];
        matr[i] = new double[n + 1];
        for (int j = 0; j < (n + 1); j++)
        {
            in >> source[i][j];
            matr[i][j] = source[i][j];
        }
    }
    double* w = new double[n];
    double** otr = new double* [n];
    for (int i = 0; i < n; i++)
    {
        otr[i] = new double[n];
        for (int j = 0; j < n; j++)
        {
            otr[i][j] = 0;
        }
    }
    double** res = new double* [n];
    for (int i = 0; i < n; i++)
    {
        res[i] = new double[n + 1];
        for (int j = 0; j < (n + 1); j++)
        {
            res[i][j] = 0;
        }
    }
    for (int k = 0; k < (n - 1); k++) //цикл по шагам отражения
    {
        double d = 0;
        for (int i = k; i < n; i++) //вычисляем скалярное произведение
        {
            d += matr[i][k] * matr[i][k];
        }

        double b = sqrt(d);
        if (matr[k][k] < 0)
        {
            b = -b;
        }
        double p = sqrt(2 * (d + b * matr[k][k]));
        for (int i = 0; i < n; i++) //вычисление вектора w
        {
            if (i < k)
            {
                w[i] = 0;
            }
            else
            {
                if (i > k)
                {
                    w[i] = matr[i][k] / p;
                }
                else
                {
                    w[i] = (matr[k][k] + b) / p;
                }
            }
        }

        for (int i = 0; i < n; i++) //вычисляем матрицу отражений
        {
            for (int j = 0; j < n; j++)
            {
                if (i != j)
                {
                    otr[i][j] = -2 * w[i] * w[j];
                }
                else
                {
                    otr[i][j] = 1 - 2 * w[i] * w[j];
                }
            }
        }

        for (int i = 0; i < n; i++) //вычисляем otr*matr
        {
            for (int j = 0; j < (n + 1); j++)
            {
                res[i][j] = 0;
                for (int t = 0; t < n; t++)
                {
                    res[i][j] = res[i][j] + otr[i][t] * matr[t][j];
                }
            }
        }
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < (n + 1); j++)
            {
                matr[i][j] = res[i][j];
            }
        }
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < (n + 1); j++)
        {
            cout << fixed << setprecision(8) << setw(12) << matr[i][j] << " ";
        }
        cout << endl;
    }
    double* x = new double[n];
    for (int i = (n - 1); i >= 0; i--)//находим x обратным ходом метода Гаусса
    {
        x[i] = matr[i][n];
        for (int j = (i + 1); j < n; j++)
        {
            x[i] = x[i] - x[j] * matr[i][j];
        }
        x[i] = x[i] / matr[i][i];
    }

    for (int i = 0; i < n; i++)
    {
        cout << fixed << setprecision(8) << setw(12) << "X" << i << " = " << x[i] << endl;
    }

    for (int i = 0; i < n; i++)
    {
        delete[] source[i];
        delete[] matr[i];
        delete[] res[i];
        delete[] otr[i];
    }
    delete[] source;
    delete[] matr;
    delete[] res;
    delete[] otr;
    delete[] w;
    delete[] x;
    return 0;
}
