#pragma once
#include <string.h>
#include <limits.h>

/*
* C qsort() is very slow, much slower than C++ std::sort().
* This is because qsort() doesn't utilize data-type information at compile time,
* and it has redundant pointer dereference since it requires a compare function.
* For projects that use old C, it's impossible to convert to C++/newer C.
* 
* So we implement a simple quicksort that is ~~25% faster than std::sort() 
* with mostly random data, and much faster with structured/sorted data 
*/

static const int INF = INT_MAX; // 2^31 - 1

/*	we use macro instead of function to avoid extra pointer dereference. */
static double TEMP_SWAP_F64;
#define swap_f64(x,y) {TEMP_SWAP_F64 = (x); (x) = (y); (y) = TEMP_SWAP_F64;}

static int TEMP_SWAP_INT;
#define swap_int(x,y) {TEMP_SWAP_INT = (x); (x) = (y); (y) = TEMP_SWAP_INT;}

/* For quicksort, we need to choose a random pivot. Any random function should work. Even bad ones.
   cycle size = 25000004
*/
static int
cheap_random()
{
	const int base = 2 * 100 * 1000 * 1000 + 33;
	static int seed = 0;
	seed = seed * 7 + 13;
	if (seed > base) seed %= base;
	return seed;
}

/* Another way to choose pivot is using median of left, middle, and right values*/
inline int
median_index_of3(const double arr[], const int a, const int b, const int c)
{
	if (arr[a] < arr[b]) {
		if (arr[b] < arr[c]) return b;
		else if (arr[a] < arr[c]) return c;
		else return a;
	}
	else {
		if (arr[a] < arr[c]) return a;
		else if (arr[b] < arr[c]) return c;
		else return b;
	}
}

//****************************
//****************************	MAIN SORTING PART

/* Insertion sort is best when the array is small*/
inline void
insertion_sort(double a[], int l, int r)
{
	if (r <= l) return;
	int i, j;
	double value;

	/* Performing 1 quicksort loops increase performance by 20-30% on random data,
	 * and completely eliminate worst-case scenario. It's slower when the array is
	 * sorted, but in that case the sort time is already tiny.
	 */
	i = l; j = r;
	value = a[(l + r) >> 1];
	while (i <= j) {
		while (a[i] < value) i++;
		while (a[j] > value) j--;

		if (i <= j) {
			swap_f64(a[i], a[j]);
			i++;
			j--;
		}
	}

	for (i = l + 1; i <= r; i++) {
		j = i - 1;
		value = a[i];
		while (j >= l && a[j] > value) {
			a[j + 1] = a[j];
			j--;
		}
		a[j + 1] = value;
	}
}

/* This function choose a pivot, move all smaller elements to the left,
 * larger ones to the right, and put the pivot in its correct position.
 * It also checks for sorted/reverse-sorted array to speed-up special cases.
 */
static void
fastsort_partition(double a[], const int left, const int right, int* first_end_ptr, int* second_start_ptr) {
	int low, high, i, pivot, mid;
	double value;
	int increasing = 1, decreasing = 1;

	/*******/
	/* choose a random way to choose pivot, to prevent all possible worst-cases*/
	if ((right - left) & 1) pivot = left + cheap_random() % (right - left);
	else pivot = median_index_of3(a, left, (left + right) >> 1, right);
	value = a[pivot];

	/*******/
	/* Skip through smaller values on left and larger values on right*/
	low = left; high = right;
	while (a[low] < value) {
		low++;
		decreasing = 0;
		if (a[low] < a[low - 1]) increasing = 0;
	}

	while (a[high] > value) {
		high--;
		decreasing = 0;
		if (a[high] > a[high + 1]) increasing = 0;
	}

	increasing &= a[high] >= a[low];
	decreasing &= a[high] <= a[low];

	/*******/
	/* Resolve degenerate input cases */
	if (increasing) {
		/* Choose randomly looping from left->right or right->left. Prevents all possible worst-cases */
		if ((right - left) & 1) {
			for (i = low + 1; i <= high; i++) if (a[i] < a[i - 1]) {
				increasing = 0;
				break;
			}
		}
		else {
			for (i = high; i >= low + 1; i--) if (a[i] < a[i - 1]) {
				increasing = 0;
				break;
			}
		}
		if (increasing) {	/* sorted */
			*first_end_ptr = INF;
			return;
		}
	}

	if (decreasing) {
		if ((right - left) & 1) {
			for (i = low + 1; i <= high; i++) if (a[i] > a[i - 1]) {
				decreasing = 0;
				break;
			}
		}
		else {
			for (i = high; i >= low + 1; i--) if (a[i] > a[i - 1]) {
				decreasing = 0;
				break;
			}
		}
		if (decreasing) {
			mid = (right - left + 1) >> 1;
			for (i = 0; i < mid; i++) swap_f64(a[left + i], a[right - i]);
			*first_end_ptr = INF;
			return;
		}
	}

	/******/
	while (low <= high) {
		while (a[low] < value) low++;
		while (a[high] > value) high--;

		if (low <= high) {
			swap_f64(a[low], a[high]);
			low++;
			high--;
		}
	}

	*first_end_ptr = high;
	*second_start_ptr = low;
}

void
fastsort_recursive(double a[], int l, int r)
{
	int pivot, first_end, second_start;
	while (l < r) {
		if (r - l <= 70) {	/* determined through experiments and benchmarks, not randomly. 70-150 works fine on random/mixed (hard) data */
			insertion_sort(a, l, r);
			return;
		}

		fastsort_partition(a, l, r, &first_end, &second_start);
		if (first_end == INF) return; /* sorted */

		/* Recurse into smaller branch to avoid stack overflow */
		if (first_end - l < r - second_start) {
			fastsort_recursive(a, l, first_end);
			l = second_start;
		}
		else {
			fastsort_recursive(a, second_start, r);
			r = first_end;
		}
	}
}

/*
We separate this from the main sort function so that
if we ever need to initialize/shuffle/etc, we can do
it here without changing the core sorting code.
*/
void
fastsort(int n, double a[])
{
	fastsort_recursive(a, 0, n - 1);
}

//****************************************//
//****************************************//
//****************************************//
//****************************************//
//****************************************//	
//******** SORTING INDEX *****************//
/* Return array index[] such that a[index[]] is sorted */

inline int
median_index_of3_index(const double arr[], int index[], const int a, const int b, const int c)
{
	if (arr[index[a]] < arr[index[b]]) {
		if (arr[index[b]] < arr[index[c]]) return b;
		else if (arr[index[a]] < arr[index[c]]) return c;
		else return a;
	}
	else {
		if (arr[index[a]] < arr[index[c]]) return a;
		else if (arr[index[b]] < arr[index[c]]) return c;
		else return b;
	}
}


/* Insertion sort is best when the array is small*/
inline void
insertion_sort_index(const double a[], int index[], int l, int r)
{
	if (r <= l) return;
	int i, j, current_index;
	double value;

	/* Performing 1 quicksort loops increase performance by 20-30% on random data,
	 * and completely eliminate worst-case scenario. It's slower when the array is
	 * sorted, but in that case the sort time is already tiny.
	 */
	i = l; j = r;
	value = a[index[(l + r) >> 1]];
	while (i <= j) {
		while (a[index[i]] < value) i++;
		while (a[index[j]] > value) j--;

		if (i <= j) {
			swap_int(index[i], index[j]);
			i++;
			j--;
		}
	}

	for (i = l + 1; i <= r; i++) {
		j = i - 1;
		value = a[index[i]];
		current_index = index[i];

		while (j >= l && a[index[j]] > value) {
			index[j + 1] = index[j];
			j--;
		}
		index[j + 1] = current_index;
	}
}

//***************
static void
fastsort_partition_index(const double a[], int index[], const int left, const int right, int* first_end_ptr, int* second_start_ptr) {
	int low, high, i, pivot, mid;
	double value;
	int increasing = 1, decreasing = 1;

	/*******/
	/* choose a random way to choose pivot, to prevent all possible worst-cases*/
	if ((right - left) & 1) pivot = left + cheap_random() % (right - left);
	else pivot = median_index_of3_index(a, index, left, (left + right) >> 1, right);
	value = a[index[pivot]];

	/*******/
	/* Skip through smaller values on left and larger values on right*/
	low = left; high = right;
	while (a[index[low]] < value) {
		low++;
		decreasing = 0;
		if (a[index[low]] < a[index[low - 1]]) increasing = 0;
	}

	while (a[index[high]] > value) {
		high--;
		decreasing = 0;
		if (a[index[high]] > a[index[high + 1]]) increasing = 0;
	}

	increasing &= a[index[high]] >= a[index[low]];
	decreasing &= a[index[high]] <= a[index[low]];

	/*******/
	/* Resolve degenerate input cases */
	if (increasing) {
		if ((right - left) & 1) {
			for (i = low + 1; i <= high; i++) if (a[index[i]] < a[index[i - 1]]) {
				increasing = 0;
				break;
			}
		}
		else {
			for (i = high; i >= low + 1; i--) if (a[index[i]] < a[index[i - 1]]) {
				increasing = 0;
				break;
			}
		}
		if (increasing) {	/* sorted */
			*first_end_ptr = INF;
			return;
		}
	}

	if (decreasing) {
		if ((right - left) & 1) {
			for (i = low + 1; i <= high; i++) if (a[index[i]] > a[index[i - 1]]) {
				decreasing = 0;
				break;
			}
		}
		else {
			for (i = high; i >= low + 1; i--) if (a[index[i]] > a[index[i - 1]]) {
				decreasing = 0;
				break;
			}
		}
		if (decreasing) {
			mid = (right - left + 1) >> 1;
			for (i = 0; i < mid; i++) swap_int(index[left + i], index[right - i]);
			*first_end_ptr = INF;
			return;
		}
	}

	/******/
	while (low <= high) {
		while (a[index[low]] < value) low++;
		while (a[index[high]] > value) high--;

		if (low <= high) {
			swap_int(index[low], index[high]);
			low++;
			high--;
		}
	}

	*first_end_ptr = high;
	*second_start_ptr = low;
}


//***************
void
fastsort_recursive_index(const double a[], int index[], int l, int r)
{
	int pivot, first_end, second_start;
	while (l < r) {
		if (r - l <= 70) {	/* determined through experiments and benchmarks, not randomly. 70-150 works fine on random/mixed (hard) data */
			insertion_sort_index(a, index, l, r);
			return;
		}

		fastsort_partition_index(a, index, l, r, &first_end, &second_start);
		if (first_end == INF) return; /* sorted */

		/* Recurse into smaller branch to avoid stack overflow */
		if (first_end - l < r - second_start) {
			fastsort_recursive_index(a, index, l, first_end);
			l = second_start;
		}
		else {
			fastsort_recursive_index(a, index, second_start, r);
			r = first_end;
		}
	}
}    

//***************
void
fastsort_index(int n, const double a[], int index[])
{
	int i;
	for (i = 0; i < n; i++) index[i] = i;
	fastsort_recursive_index(a, index, 0, n - 1);
}

