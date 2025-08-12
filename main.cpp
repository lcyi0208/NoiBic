/*Author: Chaoyi Long <chaoyilong@mail.sdu.edu.cn>,Jan 16th,2025*/

#include "structure.h"

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
vector<string> conds, genes;


int main(int argc,char* argv[])
{
	if (argc > 1 && (std::string(argv[1]) == "--version" || std::string(argv[1]) == "-V")) {
        std::cout << APP_NAME << " " << APP_VERSION
                  << " (built " << APP_BUILD_DATE << " on " << APP_HOST << ")\n";
        return 0;
    }
	//start the timer
	uglyTime(NULL);
	get_options(argc,argv);
	
	ifstream input(po->FN);
	string nums, tmp;
	n = 0, m = 0;
	//get the expression matrix
	while (getline(input, nums))
	{
		istringstream ci(nums);
		if (n == 0)
		{
			ci >> tmp;
			if(tmp!="o")
			{
				m++;
				conds.emplace_back(tmp);
			}	
			while (ci >> tmp)
			{
				m++;
				conds.emplace_back(tmp);
			}
			
		}
		else
		{
			vector<Array> temp;
			temp.reserve(m);
			ci >> tmp;
			genes.emplace_back(tmp);
			for (size_t j = 0; j < m; j++)
			{
				Array s;
				s.con = j;
				ci >> s.val;
				s.old_val = s.val;
				temp.emplace_back(s);
			}
			A.emplace_back(temp);
		}
		n++;                       
	}
	n--;
	cout<<"col_num: "<<m<<"  row_num: "<<n<<endl;
	
	
	B= vector<vector<Array> >(n, vector<Array>(m));
	B = A;
	
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
