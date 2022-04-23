#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <set>
#include <utility>      // std::pair, std::make_pair
#include <float.h>
#include <limits.h>
#include <map>

using namespace std;
typedef pair<double, double> pairs;

extern "C" {
	int HelloWorld();

	void  l1gnio(double* data, double* w, double* l_read, double* m_read, const int n, double* solution);

	void  l2gnio(double* data, double* w, double* l_read, double* m_read, const int n, double* solution);

	void obj_l1(double* data, double* w, double* lbd, double* mu, const int n, double* solution, double* obj);

	void obj_l2(double* data, double* w, double* lbd, double* mu, const int n, double* solution, double* obj);

}

int HelloWorld() {
	printf("Thanks for using our software \n");
	printf("Copyright: Xuyu Chen and Xudong Li \n");
	printf("For any bug, please contact: chenxy18 AT fudan.edu.cn");
	printf("\n");
	return 30;
}
void l1gnio(double* data, double* w, double* l_read, double* m_read, const int n, double* solution) {

	map<double, double> bp;
	double lslope = 0.0;
	double rslope = 0.0;
	double left = 0.0;
	double right = 0.0;
	bp.insert(pairs(0.0, 0.0));
	double yi = 0.0;
	double lbd = 0.0;
	double mu = 0.0;
	double wi = 0.0;
	double bnega = 0.0;
	double bposi = 0.0;
	double* leftposi = (double*)malloc((n - 1) * sizeof(double));
	double* rightposi = (double*)malloc((n - 1) * sizeof(double));

	map<double, double>::iterator sum_ptr;
	map<double, double>::iterator l_ptr;
	map<double, double>::iterator r_ptr;


	for (int j = 0; j <= (n - 2); j++) {
		yi = *data;
		wi = *w;
		lbd = *l_read;
		mu = *m_read;


		sum_ptr = bp.find(yi);
		if (sum_ptr == bp.end())
			bp.insert(pairs(yi, 2 * wi));
		else
		{
			sum_ptr->second += 2 * wi;
		};
		rslope += wi;
		lslope -= wi;



		if (-lbd <= lslope && mu >= rslope) {
			bnega = -(double)DBL_MAX;
			bposi = (double)DBL_MAX;
		}
		else if (-lbd > lslope && mu < rslope) {
			left = lslope;
			while (left < -lbd) {
				l_ptr = bp.begin();
				bnega = l_ptr->first;
				left += (l_ptr->second);
				bp.erase(bnega);
			}
			bp.insert(pairs(bnega, left + lbd));
			lslope = -lbd;

			right = rslope;
			while (right > mu)
			{
				r_ptr = bp.end();
				r_ptr--;
				bposi = r_ptr->first;
				right -= (r_ptr->second);
				bp.erase(bposi);
			}
			bp.insert(pairs(bposi, mu - right));
			rslope = mu;
		}
		else if (-lbd <= lslope && mu < rslope) {
			bnega = -(double)DBL_MAX;

			right = rslope;
			while (right > mu)
			{
				r_ptr = bp.end();
				r_ptr--;
				bposi = r_ptr->first;
				right -= (r_ptr->second);
				bp.erase(bposi);
			}
			bp.insert(pairs(bposi, mu - right));
			rslope = mu;
		}
		else {
			left = lslope;
			while (left < -lbd) {
				l_ptr = bp.begin();
				bnega = l_ptr->first;
				left += (l_ptr->second);
				bp.erase(bnega);
			}
			bp.insert(pairs(bnega, left + lbd));
			lslope = -lbd;

			bposi = (double)DBL_MAX;
		}

		leftposi[j] = bnega;
		rightposi[j] = bposi;

		data++;
		w++;
		l_read++;
		m_read++;

	};
	yi = *data;
	wi = *w;
	sum_ptr = bp.find(yi);
	if (sum_ptr == bp.end())
		bp.insert(pairs(yi, 2 * wi));
	else
	{
		sum_ptr->second += 2 * wi;
	};
	rslope += wi;
	lslope -= wi;


	double xmin = 0.0;
	double s = 0.0;
	s = lslope;
	if (s >= 0) {
		l_ptr = bp.begin();
		xmin = l_ptr->first;
	}
	else {
		while (s < 0) {
			l_ptr = bp.begin();
			xmin = (l_ptr->first);
			s += (l_ptr->second);
			bp.erase(xmin);
		}
	}

	solution[n - 1] = xmin;
	double xold = xmin;
	double xnew = 0.0;

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
	free(leftposi);
	free(rightposi);

}
void l2gnio(double* data, double* w, double* l_read, double* m_read, const int n, double* solution) {
	double* bp = NULL;
	double* df_a = NULL;
	double* df_b = NULL;
	double dleft[2] = { 0.0,0.0 };
	double dright[2] = { 0.0,0.0 };
	double* bp_i = (double*)malloc((2 * n + 1) * sizeof(double));
	double* df_a_i = (double*)malloc((2 * n + 1) * sizeof(double));
	double* df_b_i = (double*)malloc((2 * n + 1) * sizeof(double));

	bp = bp_i + n;
	df_a = df_a_i + n;
	df_b = df_b_i + n; // set the pointers to the middle of the pre-required memory


	*bp = 0.0;
	*df_a = 0.0;
	*df_b = 0.0;

	double* bp_start = bp;
	double* bp_end = bp;
	double* df_a_s = df_a;
	double* df_a_e = df_a;
	double* df_b_s = df_b;
	double* df_b_e = df_b;

	double* leftposi = (double*)malloc((n - 1) * sizeof(double));
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
		if (lbd >= 0 && lbd < (double)DBL_MAX) {
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
void obj_l1(double* data, double* w, double* lbd, double* mu, const int n, double* solution, double* obj) {

	double objective = 0.0;

	for (int i = 0; i < n; i++) {
		objective = objective + w[i] * fabs(data[i] - solution[i]);
	}

	for (int i = 0; i < (n - 1); i++) {
		if (lbd[i] < (double)DBL_MAX) {
			objective = objective + lbd[i] * fmax(0, (solution[i] - solution[i + 1]));
		}
		if (mu[i] < (double)DBL_MAX) {
			objective = objective + mu[i] * fmax(0, (solution[i + 1] - solution[i]));
		}
	}

	*obj = objective;
}
void obj_l2(double* data, double* w, double* lbd, double* mu, const int n, double* solution, double* obj) {
	double objective = 0.0;

	for (int i = 0; i < n; i++) {
		objective = objective + w[i] * (data[i] - solution[i]) * (data[i] - solution[i]);
	}

	for (int i = 0; i < (n - 1); i++) {
		if (lbd[i] < (double)DBL_MAX) {
			objective = objective + lbd[i] * fmax(0, (solution[i] - solution[i + 1]));
		}
		if (mu[i] < (double)DBL_MAX) {
			objective = objective + mu[i] * fmax(0, (solution[i + 1] - solution[i]));
		}
	}

	*obj = objective;
}
