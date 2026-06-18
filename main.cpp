/*Author: Chaoyi Long <chaoyilong@mail.sdu.edu.cn>,Jan 16th,2025*/
 
#include "structure.h"
using namespace std;
#ifndef APP_NAME
#define APP_NAME "noibic"
#endif
#ifndef APP_VERSION
#define APP_VERSION "dev"
#endif
#ifndef APP_BUILD_DATE
#define APP_BUILD_DATE "unknown"
#endif
#ifndef APP_HOST
#define APP_HOST "unknown"
#endif

size_t n, m;

vector<vector<Array> >A, B;
vector<vector<double>> NoiseMatrix;
vector<string> conds, genes;

int main(int argc,char* argv[])
{
	if (argc > 1 && (string(argv[1]) == "--version" || string(argv[1]) == "-V")) {
        cout << APP_NAME << " " << APP_VERSION
                  << " (built " << APP_BUILD_DATE << " on " << APP_HOST << ")\n";
        return 0;
    }
	//start the timer
	uglyTime(NULL);
	get_options(argc,argv);
	read_file(po->FN);

	cout << "Input matrix: rows=" << n << " cols=" << m << endl;
	
	if (po->FN_NOISE[0] != '\0')
		read_noise_matrix(po->FN_NOISE);
	
	data_preprocessing();
	seed_generation();
	cluster();
	cluster_expand();
	
	char filename[MX];
	strcpy(filename,po->FP);
	strcat(filename,".blocks");
	result_output(filename);
	return 0;
}
