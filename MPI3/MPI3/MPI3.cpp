#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <windows.h>
#include "mpi.h"



void create_random_array(int **array, int n)
{
	int *arr;
	int i;
	arr = (int*)malloc(n * sizeof(int));

	if (arr == NULL)
	{
		fprintf(stderr, "Could not malloc, exiting.\n");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}
	for (i = 0; i < n; i++) arr[i] = rand() % n;
	*array = arr;
}

int main(int argc, char *argv[])
{
	int i, j, k, check, n;
	int my_rank, num_procs, my_n;
	int compare_to, divider;
	int *a = nullptr, *my_a = nullptr, *tmp = nullptr, *tmp2 = nullptr, tmp_no;
	double t1, t2;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

	int N = 4;
	if (argc >= 2) N = atoi(argv[1]);

	check = num_procs;
	while (check > 1)
	{
		if (check % 2 != 0)
		{
			fprintf(stderr, "Number of processes must power of 2.");
			MPI_Abort(MPI_COMM_WORLD, 1);
		}
		check /= 2;
	}
	my_n = N / num_procs;               

	if (num_procs > 1)
	{
		my_a = (int*)malloc(my_n * sizeof(int));
		tmp = (int*)malloc(my_n * sizeof(int));
		tmp2 = (int*)malloc(my_n * sizeof(int));

	}

	if (my_rank == 0) create_random_array(&a, my_n * num_procs);

	/* start timing */
	if (my_rank == 0) t1 = MPI_Wtime();

	if (num_procs > 1) MPI_Scatter(a, my_n, MPI_INT, my_a, my_n, MPI_INT, 0, MPI_COMM_WORLD);
	else my_a = a;

	/* phase 1: do shellsort on own part */

	for (divider = my_n / 2; divider > 0; divider /= 2)
	{
		for (i = divider; i < my_n; i++)
		{
			tmp_no = my_a[i];
			for (j = i; j >= divider && tmp_no < my_a[j - divider]; j -= divider)
			{
				my_a[j] = my_a[j - divider];
			}
			my_a[j] = tmp_no;
		}
	}
	/* phase 2: compare-split with processes far away,
	divide and move closer */

	divider = num_procs;
	compare_to = divider - my_rank - 1;
	divider /= 2;

	while (divider > 0)
	{

		/* compare-split */
		MPI_Sendrecv(my_a, my_n, MPI_INT, compare_to, divider, tmp, my_n, MPI_INT, compare_to, divider, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		/* merge sort if different */
		if (my_rank < compare_to && tmp[0] < my_a[my_n - 1])
		{
			j = k = i = 0;
			while (k < my_n)
			{
				if (tmp[j] < my_a[i]) tmp2[k++] = tmp[j++];
				else tmp2[k++] = my_a[i++];
			}
			memcpy(my_a, tmp2, my_n * sizeof(int));
		}
		else if (my_rank > compare_to && tmp[my_n - 1] > my_a[0])
		{
			j = k = i = my_n - 1;
			while (k >= 0)
			{
				if (tmp[j] > my_a[i]) tmp2[k--] = tmp[j--];
				else tmp2[k--] = my_a[i--];
			}
			memcpy(my_a, tmp2, my_n * sizeof(int));
		}
		n = my_rank / divider;
		compare_to = (n + 1) * divider - my_rank % divider - 1;
		divider /= 2;
	}

	/* phase 3: odd even compare split */
	int flag = num_procs;
	while (--flag) {
		if (flag % 2) {
			compare_to = (my_rank % 2 == 0) ? my_rank + 1 : my_rank - 1;
		}
		else {
			compare_to = (my_rank % 2 == 0) ? my_rank - 1 : my_rank + 1;
			if (compare_to == -1 || compare_to == num_procs)
				compare_to = MPI_PROC_NULL;
		}

		/* compare-split */
		MPI_Sendrecv(my_a, my_n, MPI_INT, compare_to, 10 + flag,
			tmp, my_n, MPI_INT, compare_to, 10 + flag,
			MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		/* merge sort if different */
		if (my_rank < compare_to && tmp[0] < my_a[my_n - 1]) {
			j = k = i = 0;
			while (k < my_n) {
				if (tmp[j] < my_a[i])
					tmp2[k++] = tmp[j++];
				else
					tmp2[k++] = my_a[i++];
			}
			memcpy(my_a, tmp2, my_n * sizeof(int));
		}
		else if (my_rank > compare_to && tmp[my_n - 1] > my_a[0]) {
			j = k = i = my_n - 1;
			while (k >= 0) {
				if (tmp[j] > my_a[i])
					tmp2[k--] = tmp[j--];
				else
					tmp2[k--] = my_a[i--];
			}
			memcpy(my_a, tmp2, my_n * sizeof(int));
		}
	}

	if (num_procs > 1) {
		MPI_Gather(my_a, my_n, MPI_INT, a, my_n, MPI_INT,
			0, MPI_COMM_WORLD);
	}

	if (my_rank == 0) {
		t2 = MPI_Wtime();
		if (num_procs > 1)
		{
			if (N < 200) {
				for (int i = 0; i < N; i++)
					std::cout << a[i] << " ";
			}
		}
		printf("\n\nDone sorting; elapsed time is %f\n", t2 - t1);

		free(a);
	}

	if (num_procs>1) {
		free(my_a);
		free(tmp);
		free(tmp2);
	}

	MPI_Finalize();
	return 0;
}