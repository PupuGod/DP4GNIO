#include <stdio.h>
#include <vector>
#include <set>
#include <utility>      // std::pair, std::make_pair
#include <string>       // std::string
#include <iostream>     // std::cout
#include <time.h>
#include <limits.h>
#include <map>
#define M 1000000

using namespace std;
typedef pair<double, double> pairs;

int main(){
	// 0 Input the data
	cout <<"The L1-GNIO testfile (in C++), Version 1.0, author: Xuyu Chen \n" << endl;
	// Part I -- Reading data from python --  Be Care full of the STACK OVERFLOW

	cout << "Reading the data from python ....... \n" << endl;
	clock_t read_start, read_end;
	FILE* fp;
	read_start = clock();
	// Ensure the PATH and LENGTH befor use 
	errno_t read_state = fopen_s(&fp, "C:\\Users\\AI\\Desktop\\zzy\\random_e4_l1.txt", "r");

	if (read_state != 0) {
		cout << "Error: can't read the data \n" << endl;
		exit(0);
	}

	double data_all[(5 * M - 1)];
	//double* data_all= calloc((5 * M - 1), sizeof(double));

	int i = 0;
	while (1) {
		fscanf_s(fp, "%lf", &data_all[i]);
		i = i + 1;
		if (i >= (5 * M - 1))
			break;
	}
	fclose(fp);
	read_end = clock();
	cout << "Data has been successfully read, start computing: \n" << endl;

	clock_t comp_start, comp_end;
	comp_start = clock();
	
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



	// 1 initialization breakpoint pair
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
	double leftposi[M - 1];
	double rightposi[M - 1];


	map<double, double>::iterator sum_ptr;
	map<double, double>::iterator l_ptr;
	map<double, double>::iterator r_ptr;

	comp_start = clock();

	for (int j = 0; j <= (M - 2); j++) {
		yi = *d_read;
		wi = *w_read;
		lbd = *l_read;
		mu = *m_read;
	
		// The ``sum procedure''
		sum_ptr = bp.find( yi );
		if (sum_ptr == bp.end())
			bp.insert(pairs(yi, 2 * wi));
		else
		{
			sum_ptr->second += 2*wi;
		};
		rslope += wi;
		lslope  -= wi;


		// The ``Truncate procedure''
		if (-lbd <= lslope && mu >= rslope) {
			bnega = -(double) DBL_MAX;
			bposi = (double) DBL_MAX;
		}else if( -lbd > lslope && mu< rslope ){
			left = lslope;
			while ( left  < -lbd) {
				l_ptr = bp.begin();
				bnega = l_ptr->first;
				left += (l_ptr->second);
				bp.erase(bnega);
			}
			bp.insert(pairs(bnega, left + lbd));
			lslope = -lbd;

			right = rslope;
			while ( right > mu)
			{	
				r_ptr = bp.end();
				r_ptr--;
				bposi = r_ptr->first;
				right -= (r_ptr->second);
				bp.erase(bposi);
			}
			bp.insert(pairs(bposi, mu - right));
			rslope = mu;
		}else if (-lbd <= lslope && mu < rslope) {
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

		leftposi[j]  = bnega;
		rightposi[j] = bposi;

		d_read++;
		w_read++;
		l_read++;
		m_read++;

	}; 

	

	// Last ``sum'' procedure
	yi = *d_read;
	wi = *w_read;
	sum_ptr = bp.find(yi);
	if (sum_ptr == bp.end())
		bp.insert(pairs(yi, 2 * wi));
	else
	{
		sum_ptr->second += 2 * wi;
	};
	rslope += wi;
	lslope -= wi;



	// 3 Find optimal solution and solve xn
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

	// 4 Recover procedure

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

	cout << "Solution succssful generated! Turn to test! \n" <<
	"Computing time:" << (double)((comp_end - comp_start) * 0.001) << "s \n" << endl;
	
	
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


	return 0;



};