#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <cfloat>
#define M 1000000

int main(void) {
	printf("The L2-GNIO testfile (in C), Version 1.0, author: Xuyu Chen \n");
	// Part I -- Reading data from python --  Be Care full of the STACK OVERFLOW

	printf("Reading the data from python ....... \n");
	clock_t read_start, read_end;
	FILE* fp;
	read_start = clock();
	// Ensure the PATH and LENGTH befor use 
	errno_t read_state = fopen_s(&fp, "C:\\Users\\AI\\Desktop\\zzy\\random_e4.txt", "r");

	if (read_state != 0) {
		printf("Error: can't read the data \n");
		exit(0);
	}

	double data_all[(5 * M - 1)];
	//double* data_all= calloc((5 * M - 1), sizeof(double));

	int i = 0;
	while (1) {
		fscanf_s(fp, "%lf", &data_all[i] );
		i = i + 1;
		if (i >= (5*M - 1))
			break;
	}
	fclose(fp);
	read_end = clock();	
	printf("Data has been successfully read, start computing: \n");

	clock_t comp_start, comp_end ;
	comp_start = clock();
	// Part II -- Initialization -- Reading data, construct piecewise functions
	// II.a inintialize the piecewise function
	double* bp = NULL;
	double* df_a = NULL;
	double* df_b = NULL;
	double dleft[2] = { 0.0,0.0 };
	double dright[2] = { 0.0,0.0 };
	bp = (double*)malloc((2 * M - 1) * sizeof(double));
	df_a = (double*)malloc((2 * M - 1) * sizeof(double));
	df_b = (double*)malloc((2 * M - 1) * sizeof(double));
	bp += (M - 1);
	df_a += (M - 1);
	df_b += (M - 1); // set the pointers to the middle of the pre-required memory
	*bp = 0.0;
	*df_a = 0.0;
	*df_b = 0.0;

	double* bp_start = bp;
	double* bp_end = bp;
	double* df_a_s = df_a;
	double* df_a_e = df_a;
	double* df_b_s = df_b;
	double* df_b_e = df_b;

	double leftposi[M - 1];
	double rightposi[M - 1];

	// II.b  Reading the data via pointers => (data,weight,lbd,mu,solution,time)
	double* d_read = NULL;
	double* w_read = NULL;
	double* l_read = NULL;
	double* m_read = NULL;
	double* s_read = NULL;

	d_read = data_all;
	w_read = d_read + M;
	l_read = w_read + M;
	m_read = l_read + M - 1;
	s_read = m_read + M - 1;

    // II.c The variables in iterations
	double cur = 0.0;
	double bposi = 0.0;
	double bnega = 0.0;
	double a = 0.0;
	double b = 0.0;
	double lbd = 0.0;
	double mu = 0.0;

	// Part III -- Main Loop of GNIO 
	for (int i = 0; i <= (M - 2); i++) {

		// III.a The sum procedure
		dleft[0] += (*w_read);
		dleft[1] -= (2 * (*w_read) * (*d_read));
		dright[0] += (*w_read);
		dright[1] -= (2 * (*w_read) * (*d_read));
		
		/// III.b The truncate procedure
		lbd = *l_read;
		mu = *m_read;
		if (lbd >= 0) {
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


		if (mu >= 0) {
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

		d_read++;
		w_read++;
		l_read++;
		m_read++;

	}
	dleft[0] += (*w_read);
	dleft[1] -= (2 * (*w_read) * (*d_read));
	dright[0] += (*w_read);
	dright[1] -= (2 * (*w_read) * (*d_read));  // Updating the piecewise functions and store breakpoints;

    //Part IV -- Find xmin and the Recover procedure (left => right method)
	// May be optimized by (right => left method)

	// IV.a Find xmin 
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

	// IV.b The recover procedure

	double solution[M];
	solution[M - 1] = xmin;
	double xold = xmin;
	double xnew = 0.0;
	// The recover procedure
	for (int i = 1; i <= (M - 1); i++) {
		if (xold > rightposi[M - 1 - i]) {
			xnew = rightposi[M - 1 - i];
		}
		else if (xold < leftposi[M - 1 - i]) {
			xnew = leftposi[M - 1 - i];
		}
		else {
			xnew = xold;
		}
		solution[M - 1 - i] = xnew;
		xold = xnew;
	}
	comp_end = clock();
	printf("Solution succssful generated! Turn to test! \n");
	printf("Computing time: %.16f s \n", (double) ( (comp_end - comp_start)*0.001) );


	// Part V -- Check the solution and out put logfile
	double check_tol = 1e-6;
	double diff = 0.0;
	printf("Start testing with torrence to be %.8f .... \n", check_tol);
	for (int i = 0; i <= (M - 1); i++) {
		diff = fabs(*s_read - solution[i]);
		s_read++;
		if (diff > check_tol) {
			printf("Caution: the difference is larger than the torrence ! \n");
			printf("The difference reaches %.8f \n", diff);
			printf("The index is %d", i);
		}
	}
	printf("Generating the logfile, please wait! \n");
	// Output solution file & logfile
	FILE* log_f = NULL;
	errno_t log_state = fopen_s(&log_f, "logfile.txt", "a");

	if (log_state != 0) {
		printf("Error: can't open the file");
		exit(0);
	}
	else {
		fprintf(log_f, "======================================\n");
		fprintf(log_f, "Solution successfully generated! \n");
		fprintf(log_f, "Data size: ");
		fprintf(log_f, "%d \n", M);
		fprintf(log_f, "Data reading time:");
		fprintf(log_f, " %lf s \n ", (double)((read_end - read_start) * 0.001));
		fprintf(log_f, "Computing time:");
		fprintf(log_f, " %lf s \n ", (double)((comp_end - comp_start) * 0.001));
		fprintf(log_f, "======================================\n");
	}

	printf("Logfile successfully generated! \n");
		/*
		FILE * fout = NULL;
		errno_t out_state = fopen_s(&log_f, "csolution.txt", "w");
		for (int i = 0; i <= (M - 1); i++) {
			fprintf(fout, "%lf \n", solution[i]);
		}
		fclose(fout);
		*/

	return 3;



}