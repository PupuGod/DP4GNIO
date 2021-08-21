#include <stdio.h>
#include <vector>
#include <set>
#include <utility>      // std::pair, std::make_pair
#include <string>       // std::string
#include <iostream>     // std::cout
#include <fstream>
#include <time.h>
#include <limits.h>
#include <map>


using namespace std;
typedef pair<double, double> pairs;

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
		yi =  *data;
		wi =  *w;
		lbd = *l_read;
		mu =  *m_read;

		
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

}

int main() {
	/*Section 1 -- Initialization*/
	printf("The L1-GNIO computation software (C++ version)\nVersion 1.1\nAuthor: Xuyu Chen \n");
	printf("============================================================\n");
	printf("please check path and the size of the problem before use \n");


	const char* filepath = "Please input the path of the data here (with given format)";
	double* data = NULL;
	double* w = NULL;
	double* lbd = NULL;
	double* mu = NULL;
	const int n = 0; // please input the problem size here

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
	double* data_all = (double*)malloc((4 * n - 2) * sizeof(double));
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
	l1gnio(data, w, lbd, mu, n, solution);
	comp_end = clock();
	double comp_time = (double)(comp_end - comp_start) / CLOCKS_PER_SEC;
	printf("Solution succssful generated! \n");
	printf("Computation time for l1gnio is: %.16f s \n", comp_time );
	printf("============================================================\n");
	printf("Thanks for using our software!");

	return 2021;
}


