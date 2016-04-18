#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <cassert>
#include <omp.h>
#include <fstream>

using namespace std;

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

	return resultVector;
}

void MultiplyMatrixByVectorAndReturn(double* resultVector, double** matrix, double* vector, int size)
{
	//double startOmp = omp_get_wtime();

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

	//duration = omp_get_wtime() - startOmp;
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
	/*else
	{
		srand(time(0));

		for (int i = 0; i < size; i++)
		{
			vector[i] = (double)(rand() % 10);
		}
	}*/

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

bool ValidateArguments(int argc, char** argv)
{
	if (argc != 7)
	{
		cout << "Invalid number of arguments!\nUsage: 'programName' 'iterations' 'outFileName' 'procNum' 'sizeMin' 'sizeMax' 'sizeStep'!" << endl;
		getchar();
		return false;
	}

	if (atoi(argv[4]) > atoi(argv[5]))
	{
		cout << "sizeMax(arg 6) must be greater or equal than sizeMin(arg 5)!" << endl;
		getchar();
		return false;
	}

	if (atoi(argv[3]) < 1)
	{
		cout << "procNum(arg 4) must be greater or equal 1" << endl;
		getchar();
		return false;
	}

	return true;
}

int main(int argc, char** argv)
{
	/*if (!ValidateArguments(argc, argv))
	{
		return 0;
	}*/

	const char* fileName = "bcsstk19.mtx";
	ofstream outFile;
	int iterations, sizeMin, sizeMax, sizeStep, procNum;
	double timeOmp = 0.0, startOmp;
	double dot = 0.0;
	double** matrix = 0;
	double* vector = 0;
	double* resultVector = 0;
	char* outFileName;

	iterations = atoi(argv[1]);
	outFileName = argv[2];
	procNum = atoi(argv[3]);
	//sizeMin = atoi(argv[4]);
	//sizeMax = atoi(argv[5]);
	//sizeStep = atoi(argv[6]);

	outFile.open(outFileName);

	cout << "Starting calculations... please wait.\n" << endl;
	cout << "Size:\t\tCores:\t\tDot:\t\tTime:" << endl;

	//for (int size = sizeMin; size <= sizeMax; size += sizeStep)
	//{
		for (int proc = 1; proc <= procNum; proc++)
		{
			timeOmp = 0.0;
			omp_set_num_threads(proc);
			int size;
			//matrix = AllocateMemoryAndZerosMatrix(size);
			matrix = ReadMatrixFromFile(fileName, size);
			vector = AllocateMemoryAndZerosVector(size);
			resultVector = AllocateMemoryAndZerosVector(size);
			matrix = FillMatrix(size);
			vector = FillVector(size, true);

			startOmp = omp_get_wtime();

			for (int i = 0; i < iterations; i++)
			{
				MultiplyMatrixByVectorAndReturn(resultVector, matrix, vector, size);
				NormalizeVector(resultVector, size);
				SwapVectors(resultVector, vector, size);
			}

			MultiplyMatrixByVectorAndReturn(vector, matrix, resultVector, size);
			dot = DotVectors(vector, resultVector, size);

			timeOmp = omp_get_wtime() - startOmp;

			cout << size << "\t\t" << proc << "\t\t" << dot << "\t\t" << timeOmp << endl;
			outFile << size << "\t" << proc << "\t" << dot << "\t" << timeOmp << endl;
		}
	//}

	outFile.close();

	//for (int i = 0; i < size; i++)
	//{
	//	delete[] matrix[i];
	//}
	//delete[] matrix;
	//delete[] vector;
	//delete[] resultVector;

	cout << "\nDone. Check data file. Press key to close." << endl;
	getchar();
	return EXIT_SUCCESS;
}