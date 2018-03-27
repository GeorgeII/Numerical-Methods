#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

void printMatrix(const double* matrix, int n)
{
	for (int i = 0; i < n; i++)	{
		for (int j = 0; j < n + 1; j++)	{
			cout << setw(8) << matrix[i * (n + 1) + j];
		}
		cout << endl;
	}
	cout << endl;
}

//check element for being max in a column
bool maxElem (const double* matrix, int &maxRow, int i, int n)
{
	int j = i;
	bool flag = true;
	double max = fabs(matrix[(n + 1) * i + j]);
	for (++i; i < n; i++)
	{
		if ( fabs(matrix[(n + 1) * i + j]) > max )
		{
			maxRow = i;
			max = fabs(matrix[(n + 1) * i + j]);
			flag = false;
		}
	}

	return flag;
}

//swap 'i' row and 'maxRow'
void swapRows(double* matrix, int maxRow, int i, int n)
{
	double temp;
	for (int j = 0; j < n + 1; j++)
	{
		temp = matrix[i * (n + 1) + j];
		matrix[i * (n + 1) + j] = matrix[maxRow * (n + 1) + j];
		matrix[maxRow * (n + 1) + j] = temp;
	}
	cout << "swapping  " << i << endl;
	printMatrix(&matrix[0], n);
}

int main()
{
	const int n = 3;

	double matrix[n][n + 1] = { 
				{ 1, 0, 1, 7 },
				{ 2, -1, -1, 0 },
				{ 1, -2, -1, 2 }
	};

	printMatrix(&matrix[0][0], n);

	
	for (int i = 0; i < n; i++) {
		int maxRow = i;

		//if an element is not max in a column
		if (!maxElem(&matrix[0][0], maxRow, i, n))	
		{
			swapRows(&matrix[0][0], maxRow, i, n);
		}

		//copy i row so we can overwrite i-row in matrix with no losses
		double copiedRow[n + 1];
		for (int j = i; j < n + 1; j++)
		{
			copiedRow[j] = matrix[i][j];
		}

		for (int k = i; k < n; k++)
		{
			for (int j = i; j < n + 1; j++) {
				if (k == i)
				{
					matrix[k][j] /= copiedRow[i];
				}
				else
				{
					matrix[k][j] -= copiedRow[j] * matrix[k][i] / copiedRow[i];
				}
			}
		}
	}

	cout << "1st step is done:" << endl;
	printMatrix(&matrix[0][0], n);


	//2nd step begins

	//solution vector
	double x[n];

	for (int i = n - 1; i >= 0; i--)
	{
		double sum = 0;
		for (int j = n - 1; j > i; j--)
		{
			sum += matrix[i][j] * x[j];
		}
		x[i] = matrix[i][n] - sum;
	}

	for (int i = 0; i < n; i++)
		cout << setw(6) << x[i];
	
	return 0;
}