/** 
 * Quicksort implementation for practical course
 **/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstring>
#include <omp.h>
#include <algorithm>
#include <iterator>

#include "Stopwatch.h"

//#define NUM_THREAD 32
#define PARALLEL 1 //Comment out if Serial version is needed
#define LENGTH 10000000
//#define FINAL (LENGTH / 48)

void print_list(double *data, int length)
{
	int i;
	for (i = 0; i < length; i++)
		printf("%e\t", data[i]);
	printf("\n");
}

void quicksort(double *data, int length, int num_thread)
{
	if (length <= 1)
		return;

// #pragma omp critical
// 	{
// 		print_list(data, length);
// 	}

	double pivot = data[0];
	double temp;
	int left = 1;
	int right = length - 1;

	do
	{
		while (left < (length - 1) && data[left] <= pivot)
			left++;
		while (right > 0 && data[right] >= pivot)
			right--;

		/* swap elements */
		if (left < right)
		{
			temp = data[left];
			data[left] = data[right];
			data[right] = temp;
		}

	} while (left < right);

	if (data[right] < pivot)
	{
		data[0] = data[right];
		data[right] = pivot;
	}

//print_list(data, length);

/* recursion */
#ifdef PARALLEL

#pragma omp task default(none) shared(data, num_thread) firstprivate(right) final(right < LENGTH/num_thread)
	{
		quicksort(data, right,num_thread);
	}
#pragma omp task default(none) shared(data,num_thread) firstprivate(length, left) final(length - left < LENGTH/num_thread)
	{
		quicksort(&(data[left]), length - left,num_thread);
	}
#pragma omp taskwait
#else
	quicksort(data, right,1);
	quicksort(&(data[left]), length - left,1);
#endif
}

int check(double *data, int length)
{
	int i;
	for (i = 1; i < length; i++)
		if (data[i] < data[i - 1])
			return 1;
	return 0;
}

int main(int argc, char **argv)
{
	int length;
	double *data;

	int mem_size;

	int i, j, k;

	length = LENGTH;
	if (argc > 1)
	{
		length = atoi(argv[1]);
	}

	data = (double *)malloc(length * sizeof(double));
	if (0 == data)
	{
		printf("memory allocation failed");
		return 0;
	}

	/* initialisation */
	srand(0);
	for (i = 0; i < length; i++)
	{
		data[i] = (double)rand() / (double)RAND_MAX;
	}

	Stopwatch stopwatch;
	double time = 0;
	double *data_cpy = (double *)malloc(length * sizeof(double));

	for (size_t i = 0; i < 10; i++)
	{
		std::copy(data, data+length, data_cpy);
		stopwatch.start();
		//print_list(data, length);
		#ifdef PARALLEL
		//omp_set_num_threads(NUM_THREAD);
		#pragma omp parallel
		{
			#pragma omp single
			quicksort(data_cpy, length,omp_get_max_threads());
		}
		#else
		quicksort(data_cpy, length,1);
		#endif
		time += stopwatch.stop();
		if (check(data_cpy, length) != 0)
		{
			printf("Quicksort incorrect.\n");
		}
	}

	//print_list(data, length);

	printf("Size of dataset: %d, elapsed time[s]: %e \n", length, time/10.);


	return (0);
}
