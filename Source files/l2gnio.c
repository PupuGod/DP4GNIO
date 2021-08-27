#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <cfloat>



void l2gnio(double* data, double* w, double* l_read, double* m_read, const int n, double* solution )   {
	double* bp = NULL;
	double* df_a = NULL;
	double* df_b = NULL;
	double dleft[2] = { 0.0,0.0 };
	double dright[2] = { 0.0,0.0 };
	bp = (double*)malloc((2 * n - 1) * sizeof(double));
	df_a = (double*)malloc((2 * n - 1) * sizeof(double));
	df_b = (double*)malloc((2 * n - 1) * sizeof(double));
	bp += (n - 1);
	df_a += (n - 1);
	df_b += (n - 1); // set the pointers to the middle of the pre-required memory
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
				*df_a_e = a;
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

}


int main() {
/*Section 1 -- Initialization*/
	printf("The L2-GNIO computation software (C version)\nVersion 1.1\nAuthor: Xuyu Chen \n");
	printf("============================================================\n");
	printf("please check path and the size of the problem before use \n");
	

	char* filepath = "please input the path of the data here (with required format)";
	const int n = 2; // please input the problem size here
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
	errno_t read_state = fopen_s(&fp, filepath, "r");
	if (read_state != 0) {
		printf("Error: can't read the data \n");
		exit(0);
	}	
	double* data_all = (double*) malloc( (4 * n - 2) * sizeof(double));
	int i = 0;
	while (1) {
		fscanf_s(fp, "%lf", &data_all[i]);
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
	printf("Thanks for using our software!");

	return 2021;
}


