#include<iostream>
#include<iomanip>
#include<math.h>
#include<stdlib.h>
#include<Windows.h>

#define f1(x,y,z,w)  (11.2+0.9*y-1.2*z-0.4*w)/14.4
#define f2(x,y,z,w)  (-20.1+0.9*x-0.8*y-0.9*w)/20.6
#define f3(x,y,z,w)  (13.9-1.2*x-0.8*y-1.3*w)/19.6
#define f4(x,y,z,w)  (10.7-0.4*x-0.4*y-1.3*z)/17.6

#define   SIZE   10

using namespace std;

int countt=0;

void setcolor(unsigned short color)
{
	HANDLE hcon = GetStdHandle(STD_OUTPUT_HANDLE);
	SetConsoleTextAttribute(hcon, color);
}


void Gauss()
{
	float a[SIZE][SIZE], x[SIZE], ratio;
	int i, j, k, n;

	cout << setprecision(3) << fixed;

	setcolor(3);
	cout << "\nEnter the number of unknowns: ";
	setcolor(15);
	cin >> n;
	cout << "\nEnter the matrix:"; setcolor(4); cout << "(including the coefficients of b) : " << endl <<endl; setcolor(15);
	for (i = 1; i <= n; i++)
	{
		for (j = 1; j <= n + 1; j++)
		{
			cout << "a[" << i << "]" << j << "] = ";
			cin >> a[i][j];
		}
	}

	for (i = 1; i <= n - 1; i++)
	{
		for (j = i + 1; j <= n; j++)
		{
			ratio = a[j][i] / a[i][i];

			for (k = 1; k <= n + 1; k++)
			{
				countt++;
				a[j][k] = a[j][k] - ratio * a[i][k];
			}
		}
	}

	x[n] = a[n][n + 1] / a[n][n];

	for (i = n - 1; i >= 1; i--)
	{
		x[i] = a[i][n + 1];
		for (j = i + 1; j <= n; j++)
		{
			countt++;
			x[i] = x[i] - a[i][j] * x[j];
		}
		x[i] = x[i] / a[i][i];
	}

	setcolor(4);
	cout << endl << "\nSolutions: " << endl;
	setcolor(15);
	for (i = 1; i <= n; i++)
	{
		cout << "x" << i << " = " << x[i] << endl;
	}
	setcolor(4);
	cout << "\nNumber of iterations: "; 
	setcolor(15); 
	cout << countt;
}
void Jacobi()
{
	double x0 = 0, y0 = 0, z0 = 0, w0 = 0, x1, y1, z1, w1, e1, e2, e3, e4, e;
	int pas = 1;

	cout << setprecision(3) << fixed;
	cout << "_______________________________________________________________";
	setcolor(9);  cout << endl << "k"; setcolor(15); cout << " | \tx1\t\tx2\t\tx3\t\tx4" <<setw(5) <<"|" << endl;
	cout << "__|___________________________________________________________|" <<endl;

	int i = 0;
	do
	{
		x1 = f1(x0, y0, z0, w0);
		y1 = f2(x0, y0, z0, w0);
		z1 = f3(x0, y0, z0, w0);
		w1 = f4(x0, y0, z0, w0);
		setcolor(3);
		cout << pas; setcolor(15); cout << " |\t"; 
		cout << x1 << "\t\t" << y1 << "\t\t" << z1 << "\t\t" << w1 << " |" << endl;
	
		e1 = fabs(x0 - x1);
		e2 = fabs(y0 - y1);
		e3 = fabs(z0 - z1);
		e4 = fabs(w0 - w1);

		pas++;

		x0 = x1;
		y0 = y1;
		z0 = z1;
		w0 = w1;

		i++;
		countt++;

	} while (i < 3);
	cout << "__|___________________________________________________________|";
	cout << endl << "Solutions: x1 = " << x1 << ", x2 = " << y1 << ", x3 = " << z1 << ", x4 = " << w1 << endl;
	setcolor(4);
	cout << "\nNumber of iterations: ";
	setcolor(15);
	cout << countt;
}

void GaussSeidel(int n, double* A, double* b, double* x)
{
	int max_iteration = 20;
	double epsilon = 0.001;
	double max_err = 0;

	for (int k = 1; k <= max_iteration; k++)
	{
		max_err = 0.0;

		for (int i = 0; i < n; i++)
		{
			double sum = 0.0;

			for (int j = 0; j < n; j++)
			{
				if (i != j)
					sum = sum + A[i * n + j] * x[j];
				    countt++;
			}

			double aux = (b[i] - sum) / A[n * i + i];

			double relative_err = abs((x[i] - aux) / aux);

			if (relative_err > max_err)
				max_err = relative_err;
			    countt++;

			x[i] = aux;
		}

		if (max_err < epsilon)
		{
			countt++;
			return;
		}
	}
}

void Gauss_Seidel()
{

	int n;
	setcolor(3);
	cout << "\nEnter the number of unknowns: ";
	setcolor(15);
	cin >> n;

	double* A = new double[n * n];
	double* b = new double[n];
	double* x = new double[n];

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			A[i * n + j] = 0;
		}
	}

	for (int i = 0; i < n; i++)
	{
		b[i] = 0;
		x[i] = 0;
	}

	setcolor(3);
	cout << "\nEnter the matrix that already respects the convergence condition:\n";
	setcolor(15);

	int indexi = 0;

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << "a[" << i + 1 << "][" << j + 1 << "] = ";
			cin >> A[indexi];
			indexi++;
		}
	}

	setcolor(3);
	cout << "\n\nEnter the b's: \n\n";
	setcolor(15);

	for (int i = 0; i < n; i++)
	{
		cout << "b" << i + 1 << " = ";
		cin >> b[i];
	}

	setcolor(3);
	cout << "\nMatrix A:\n\n";
	setcolor(15);
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			cout << A[i * n + j] << "\t";
		}
		cout << endl;
	}

	setcolor(3);
	cout << "\nMatrix B: \n\n";
	setcolor(15);
	for (int i = 0; i < n; i++)
	{
		cout << b[i] << endl;
	}

	GaussSeidel(n, A, b, x);

	setcolor(4);
	cout << "\nSolutions: ";
	setcolor(15);
	for (int i = 0; i < n; i++)
	{
		cout << "\nx" << (i + 1) << " = " << x[i];
	}
	cout << "\n";
	setcolor(4);
	cout << "\nNumber of iterations: ";
	setcolor(15);
	cout << countt ;

	delete[]A;
	delete[]b;
	delete[]x;
}

int main()
{
	int option;

	setcolor(13);
	cout << "Matrix calculation methods:"; setcolor(15);
	cout << "\n1.Gauss\n2.Jacobi\n3.Gauss - Seidel\n";
	setcolor(13);
	cout << "\n\nChoose the method: "; setcolor(15);
	cin >> option;

	switch (option)
	{
	case 1:
		Gauss();
		break;
	case 2:
		Jacobi();
		break;
	case 3:
		Gauss_Seidel();
		break;
	}

	return 0;
}


