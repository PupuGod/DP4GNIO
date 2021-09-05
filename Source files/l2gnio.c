#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <limits.h>
#include <cfloat>



void l2gnio(double* data, double* w, double* l_read, double* m_read, const int n, double* solution )   {
	double* bp = NULL;
	double* df_a = NULL;
	double* df_b = NULL;
	double dleft[2] = { 0.0,0.0 };
	double dright[2] = { 0.0,0.0 };
	double* bp_i = (double*)malloc((2 * n + 1) * sizeof(double));
	double* df_a_i = (double*)malloc((2 * n + 1) * sizeof(double));
	double* df_b_i = (double*)malloc((2 * n + 1) * sizeof(double));
	
	bp = bp_i + n;
	df_a = df_a_i + n ;
	df_b = df_b_i + n ;  // set the pointers to the middle of the pre-required memory
	
	*bp = 0.0;
	*df_a = 0.0;
	*df_b = 0.0;

	double* bp_start = bp;
	double* bp_end = bp;
	double* df_a_s = df_a;
	double* df_a_e = df_a;
	double* df_b_s = df_b;
	double* df_b_e = df_b;

	double* leftposi = (double*)malloc( (n - 1) * sizeof(double));
	double* rightposi = (double*)malloc((n - 1) * sizeof(double));
		
	double cur = 0.0;
	double bposi = 0.0;
	double bnega = 0.0;
	double a = 0.0;
	double b = 0.0;
	double lbd = 0.0;
	double mu = 0.0;

	for (int i = 0; i <= (n - 2); i++) {

		// III.a The sum procedure
		dleft[0] += (*w);
		dleft[1] -= (2 * (*w) * (*data));
		dright[0] += (*w);
		dright[1] -= (2 * (*w) * (*data));

		/// III.b The truncate procedure
		lbd = *l_read;
		mu = *m_read;
		if (lbd >= 0 && lbd< (double)DBL_MAX) {
			a = dleft[0];
			b = dleft[1]; // read leftmost coefficients

			while (bp_start <= bp_end) {
				cur = *bp_start; // read the leftmost bp
				if (-lbd <= (2 * cur * a + b)) {
					bnega = -(b + lbd) / (2 * a);
					bp_start--;
					*bp_start = bnega;
					leftposi[i] = bnega;
					dleft[0] = 0;
					dleft[1] = -lbd; // reviese the leftmost coefficients
					df_a_s--;
					df_b_s--; // left-move difference pointers
					*df_a_s = a;
					*df_b_s = b + lbd; // add new values of differences
					break;
				}

				a += (*df_a_s);
				b += (*df_b_s);
				bp_start++;
				df_a_s++;
				df_b_s++;
			}
			if (bp_start > bp_end) {
				bnega = -(b + lbd) / (2 * a);
				*bp_start = bnega;
				leftposi[i] = bnega;
				dleft[0] = 0;
				dleft[1] = -lbd; // reviese the leftmost coefficients
				bp_end = bp_start;
				df_a_e = df_a_s;
				df_b_e = df_b_s;
				*df_a_s = a;
				*df_b_s = b + lbd;

			}

		}
		else {
			bnega = -(double)DBL_MAX;
			leftposi[i] = bnega;
		}


		if (mu >= 0 && mu < (double)DBL_MAX) {
			a = dright[0];
			b = dright[1]; // read rightmost coefficients

			while (bp_end >= bp_start) {
				cur = *bp_end; // read the rightmost bp
				if (mu >= (2 * cur * a + b)) {
					bposi = (mu - b) / (2 * a);
					bp_end++;
					*bp_end = bposi;
					rightposi[i] = bposi;
					dright[0] = 0;
					dright[1] = mu;

					df_a_e++;
					df_b_e++;

					*df_a_e = -a;
					*df_b_e = mu - b;
					break;
				}
				a -= (*df_a_e);
				b -= (*df_b_e);
				bp_end--;
				df_a_e--;
				df_b_e--;
			}
			if (bp_end < bp_start) {
				bposi = (mu - b) / (2 * a);
				*bp_end = bposi;

				rightposi[i] = bposi;

				dright[0] = 0;
				dright[1] = mu; // reviese the leftmost coefficients

				bp_start = bp_end;
				df_a_s = df_a_e;
				df_b_s = df_b_e;
				*df_a_e = -a;
				*df_b_e = mu - b;

			}

		}
		else {
			bposi = (double)DBL_MAX;
			rightposi[i] = bposi;
		}

		data++;
		w++;
		l_read++;
		m_read++;

	}
	dleft[0] += (*w);
	dleft[1] -= (2 * (*w) * (*data));
	dright[0] += (*w);
	dright[1] -= (2 * (*w) * (*data));


	double xmin = 0.0;
	a = dleft[0];
	b = dleft[1];

	while (bp_start <= bp_end) {
		cur = *bp_start;
		if (0 <= (2 * cur * a + b)) {
			xmin = -(b) / (2 * a);
			break;
		}
		a += (*df_a_s);
		b += (*df_b_s);
		bp_start++;
		df_a_s++;
		df_b_s++;
	}
	if (bp_start > bp_end) {
		xmin = -(b) / (2 * a);
	}

	solution[n - 1] = xmin;
	double xold = xmin;
	double xnew = 0.0;
	// The recover procedure
	for (int i = 1; i <= (n - 1); i++) {
		if (xold > rightposi[n - 1 - i]) {
			xnew = rightposi[n - 1 - i];
		}
		else if (xold < leftposi[n - 1 - i]) {
			xnew = leftposi[n - 1 - i];
		}
		else {
			xnew = xold;
		}
		solution[n - 1 - i] = xnew;
		xold = xnew;
	}
	free(leftposi); free(rightposi); free(bp_i); free(df_a_i); free(df_b_i);
}


int main() {
/*Section 1 -- Initialization*/
	printf("The L2-GNIO computation software (C version)\nVersion 1.1 \n");
	printf("============================================================\n");
	printf("please check path and the size of the problem before use \n");
	

	char* filepath = "toy_example.txt"; //please input the path of the data here (with required format)
	const int n = 10; // please input the problem size here
	double* data = NULL;
	double* w = NULL;
	double* lbd = NULL;
	double* mu = NULL;
	

	printf("The file path of data is:   %s \n", filepath);
	printf("The problem size is:   %ld \n", n);
	printf("============================================================\n");
	printf("Reading the data from files, please wait... \n");


/*Section 2 -- Read the data from sources  */

	FILE* fp;
	fp = fopen( filepath, "r");

	double* data_all = (double*) malloc( (4 * n - 2) * sizeof(double));
	int i = 0;
	while (1) {
		fscanf(fp, "%lf", &data_all[i]);
		i = i + 1;
		if (i >= (4 * n - 2))
			break;
	}
	fclose(fp);

	data = data_all;
	w = data + n;
	lbd = w + n;
	mu = lbd + n - 1;
	printf("Data has been successfully read, start computation: \n");
	printf("============================================================\n");

/*Section 3 -- Computing the GNIO problem*/
	double* solution = (double*)malloc(n * sizeof(double));
	clock_t comp_start, comp_end;
	comp_start = clock();
	l2gnio(data, w, lbd, mu, n, solution);
	comp_end = clock();
	double comp_time = (double)(comp_end - comp_start) / CLOCKS_PER_SEC;


	printf("Solution succssful generated! \n");
	printf("Computation time for l2gnio is: %.16f s \n", comp_time);
	printf("============================================================\n");
	
/*Section 4 -- Output the solution (optional)
In default setting, solution won't be output, if you want to
save the solution for further use, please switch `save_solution` to 1
*/
	int save_solution = 0;
	FILE* fout = NULL;
	double number = 0;
	if (save_solution) {
		fout = fopen( "solution.txt", "w");
		for (int i = 0; i < n; i++) {
			number = solution1[i];
			fprintf(fout, "%lf\n", number);
		}
		printf("The solution is saved in solution.txt\n");
	}
	fclose(fout);
	
	
	
	
	
	
	
	
	
	
	
	
	printf("Thanks for using our software!\n");

	return 2021;
}



