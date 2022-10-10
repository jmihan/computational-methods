#include <iostream>
#include <iomanip>
#include <fstream>
using namespace std;
void inversion(double**, int);
void switchstr(double**, int, int, double*);
int main()
{
    ifstream fin;
    fin.open("input.txt");

    int n;//размер матрицы
    fin >> n;
    
    cout << "Input matrix: " << endl;
    double **matr = new double*[n];//матрица A
    for (int i = 0; i < n; i++)
    {
        cout << "(";
        matr[i] = new double[n];
        for (int j = 0; j < n; j++)
        {
            fin >> matr[i][j];
            cout << fixed << setprecision(8) << setw(12) << matr[i][j];
        }
        cout << "  )" << endl;
    }
    double** reverse = new double* [n]; //обратная матрица, пока что копия исходной
    for (int i = 0; i < n; i++)
    {
        reverse[i] = new double[n];
        for (int j = 0; j < n; j++)
        {
            reverse[i][j] = matr[i][j];
        }
    }

    cout << "Column b: " << endl;
    double* b = new double[n]; //столбец b
    for (int i = 0; i < n; i++)
    {
        fin >> b[i];
        cout << b[i] << endl;
    }

    
    for (int k = 0; k < n - 1; k++) //схема единственного деления(метод Гаусса)
    {
        double ved = matr[k][k]; //ведущий элемент
        if (ved == 0)
        {
            switchstr(matr, k, n, b);
        }
        ved = matr[k][k];
        for (int i = (k + 1); i < n; i++)
        {
            double ma = matr[i][k]; 
            for (int j = k; j < n; j++)
            {
                matr[i][j] = matr[i][j] - (ma * matr[k][j] / ved);
            }
            b[i] = b[i] - (ma * b[k] / ved);
        }
    }
    cout << "Upper triangular matrix: " << endl;
    for (int i = 0; i < n; i++)
    {
        cout << "(";
        for (int j = 0; j < n; j++)
        {
            cout << fixed << setprecision(8) << setw(12) << matr[i][j];
        }
        cout<<"  )\t"<<b[i] << endl;
    }
    
    double* x = new double[n];

    for (int k = n - 1; k > 0; k--)
    {
        double temp = matr[k][k];
        for (int i = k - 1; i >= 0; i--)
        {
            double temp1 = matr[i][k];
            for (int j = k; j >= 0; j--)
            {
                matr[i][j] -= (matr[k][j] / temp) * temp1;
            }
            
            b[i] = b[i] - b[k] * (temp1 / temp);
        }
    }
    for (int i = 0; i < n; i++)
    {
        x[i] = b[i] / matr[i][i];
    }

    cout << "Solutions: " << endl;
    for (int i = 0; i < n; i++)
    {
        cout << "X"<<i<<" = " << x[i] << endl;
    }

    double det=1;
    for (int i = 0; i < n; i++) //находим дискриминант
    {
        det *= matr[i][i];
    }
    cout << "Determinant: " << det << endl;

    inversion(reverse, n); //нахождение обратной матрицы
    cout << "Reverse matrix: " << endl;
    for (int i = 0; i < n; i++)
    {
        cout << "(";
        for (int j = 0; j < n; j++)
        {
            cout << fixed << setprecision(8) << setw(12) << reverse[i][j];
        }
        cout << "  )" << endl;
        delete[] matr[i];
        delete[] reverse[i];
    }
    delete[] matr;
    delete[] reverse;
    delete[] b;
    delete[] x;
    fin.close();
    return 0;
}

void inversion(double** a, int n) //функция вычисления обратной матрицы
{
    double temp;

    double** e = new double*[n]; //единичная матрица

    for (int i = 0; i < n; i++)
    {
        e[i] = new double[n];
    }

    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            if (i == j)
            {
                e[i][j] = 1;
            }
            else
            {
                e[i][j] = 0;
            }
        }
    }

    for (int k = 0; k < n; k++)
    {
        temp = a[k][k];
        if (temp == 0)
        {
            switchstr(a, k, n, NULL);
        }
        temp = a[k][k];
        for (int j = 0; j < n; j++)
        {
            a[k][j] /= temp;
            e[k][j] /= temp;
        }

        for (int i = k + 1; i < n; i++)
        {
            temp = a[i][k];

            for (int j = 0; j < n; j++)
            {
                a[i][j] -= a[k][j] * temp;
                e[i][j] -= e[k][j] * temp;
            }
        }
    }

    for (int k = n - 1; k > 0; k--)
    {
        for (int i = k - 1; i >= 0; i--)
        {
            temp = a[i][k];

            for (int j = 0; j < n; j++)
            {
                a[i][j] -= a[k][j] * temp;
                e[i][j] -= e[k][j] * temp;
            }
            
        }
    }
    
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            a[i][j] = e[i][j];
        }
    }

    for (int i = 0; i < n; i++)
    {
        delete[] e[i];
    }

    delete[] e;
}
void switchstr(double** m, int k, int n, double* st)//функция смены строк
{
    for (int i = (k + 1); i < n; i++)
    {
        double buf;
        if (m[i][k] != 0)
        {
            for (int j = 0; j < n; j++)
            {
                buf = m[k][j];
                m[k][j] = m[i][j];
                m[i][j] = buf;

                buf = st[k];
                st[k] = st[i];
                st[i] = buf;
            }
            
            return;
        }
    }
}