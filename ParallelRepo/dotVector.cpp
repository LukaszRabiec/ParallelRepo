#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <cassert>
#include <omp.h>

double** AllocateMemoryAndZerosMatrix(int size)
{
	double** matrix = new double*[size];

	for (int i = 0; i < size; i++)
	{
		matrix[i] = new double[size];

		for (int j = 0; j < size; j++)
		{
			matrix[i][j] = 0;
		}
	}

	return matrix;
}

double* AllocateMemoryAndZerosVector(int size)
{
	double* vector = new double[size];

	for (int i = 0; i < size; i++)
	{
		vector[i] = 0;
	}

	return vector;
}

double* MultiplyMatrixByVector(double** matrix, double* vector, int size, double& duration)
{
	double* resultVector = AllocateMemoryAndZerosVector(size);
	double startOmp = omp_get_wtime();

#pragma omp parallel for shared(resultVector, matrix, vector, size)	
	for (int i = 0; i < size; i++)
	{
		double sum = 0;
		for (int j = 0; j < size; j++)
		{
			sum += matrix[i][j] * vector[j];
		}
		resultVector[i] = sum;
	}
	duration = omp_get_wtime() - startOmp;

	return resultVector;
}

void MultiplyMatrixByVectorAndReturn(double* resultVector, double** matrix, double* vector, int size, double& duration)
{
	double startOmp = omp_get_wtime();

#pragma omp parallel for shared(resultVector, matrix, vector, size)	
	for (int i = 0; i < size; i++)
	{
		double sum = 0;
		for (int j = 0;j < size; j++)
		{
			sum += matrix[i][j] * vector[j];
		}
		resultVector[i] = sum;
	}

	duration = omp_get_wtime() - startOmp;
}

double** ReadMatrixFromFile(const char* fileName, int &size)
{
	FILE* file = 0;
	file = fopen(fileName, "r");

	if (!file)
	{
		return 0;
	}

	char ignore = ' ';

	while (ignore != '\n')
	{
		fscanf(file, "%c", &ignore);
	}

	int n, m, length;
	fscanf(file, "%d", &n);
	fscanf(file, "%d", &m);
	if (n == m)
	{
		size = n;
	}
	else
	{
		std::cout << "The matrix must be square!" << std::endl;
		getchar();
		return 0;
	}
	fscanf(file, "%d", &length);

	double** matrix = AllocateMemoryAndZerosMatrix(size);

	int row, col;
	double value;
	while (!feof(file))
	{
		fscanf(file, "%d %d", &row, &col);
		fscanf(file, "%lf", &value);
		matrix[row - 1][col - 1] = value;
		if (row != col)
		{
			matrix[col - 1][row - 1] = value;
		}
	}

	fclose(file);
	return matrix;
}

double* ReadVectorFromFile(const char* fileName, int size)
{
	FILE* file = 0;
	file = fopen(fileName, "r");

	if (!file)
	{
		return 0;
	}

	double* vector = AllocateMemoryAndZerosVector(size);

	for (int i = 0; i < size; i++)
	{
		fscanf(file, "%lf", &vector[i]);
	}

	fclose(file);
	return vector;
}

double** FillMatrix(int size)
{
	double** matrix = AllocateMemoryAndZerosMatrix(size);

	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			matrix[i][j] = 0.1;

			if (i == j)
			{
				matrix[i][j] = 1.0;
			}
		}
	}

	return matrix;
}

double* FillVector(int size, bool inSequenceVector)
{
	double* vector = AllocateMemoryAndZerosVector(size);

	if (inSequenceVector)
	{
		for (int i = 0; i < size; i++)
		{
			vector[i] = ((double)i) + 1.0;
		}
	}
	else
	{
		srand(time(0));

		for (int i = 0; i < size; i++)
		{
			vector[i] = (double)(rand() % 10);
		}
	}

	return vector;
}

double VectorLength(double* vector, int size)
{
	double result = 0;

	for (int i = 0; i < size; i++)
	{
		result += vector[i] * vector[i];
	}
	result = sqrt(result);

	return result;
}

void NormalizeVector(double* vector, int size)
{
	double length = VectorLength(vector, size);

	for (int i = 0; i < size; i++)
	{
		vector[i] = vector[i] / length;
	}
}

void SwapVectors(double* vector1, double* vector2, int size)
{
	double tmp = 0;

	for (int i = 0; i < size; i++)
	{
		tmp = vector1[i];
		vector1[i] = vector2[i];
		vector2[i] = tmp;
	}
}

double DotVectors(double* vector1, double* vector2, int size)
{
	double result = 0;

	for (int i = 0; i < size; i++)
	{
		result += vector1[i] * vector2[i];
	}

	return result;
}

void TestDot()
{
	double tab1[2] = { 1.0, 2.0 };
	double tab2[2] = { 3.0, 4.0 };

	assert(DotVectors(tab1, tab2, 2) == 11.0);
}

int main(int argc, char** argv)
{
	if (argc != 4)
	{
		std::cout << "Invalid number of arguments!\nCorrect values: 'programName' 'matrixSize' 'iterationsNumber' 'procNumber'!" << std::endl;
		getchar();
		return 0;
	}

	int size, iterations, proc;

	size = atoi(argv[1]);
	iterations = atoi(argv[2]);
	proc = atoi(argv[3]);

	double** matrix = AllocateMemoryAndZerosMatrix(size);
	double* vector = AllocateMemoryAndZerosVector(size);
	matrix = FillMatrix(size);
	vector = FillVector(size, true);
	NormalizeVector(vector, size);

	omp_set_num_threads(proc);

	double time = 0.0;
	double* resultVector = AllocateMemoryAndZerosVector(size);

	for (int i = 0; i < iterations; i++)
	{
		MultiplyMatrixByVectorAndReturn(resultVector, matrix, vector, size, time);
		NormalizeVector(resultVector, size);
		SwapVectors(resultVector, vector, size);
	}
	
	printf("Matrix size: %dx%d\n", size, size);
	printf("Number of iterations: %d\n", iterations);
	printf("\nTime for %d threads: %.13fs\n", proc, time);
	//printf("\nVector length: %f\n", VectorLength(resultVector, size));
	TestDot();

	MultiplyMatrixByVectorAndReturn(vector, matrix, resultVector, size, time);

	printf("DOT: %.13e\n", DotVectors(vector, resultVector, size)); 

	for (int i = 0; i < size; i++)
	{
		delete[] matrix[i];
	}
	delete[] matrix;
	delete[] vector;
	delete[] resultVector;

	getchar();
	return EXIT_SUCCESS;
}