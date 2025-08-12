#include "structure.h"

bool CMP_node(const Node &a, const Node &b)
{
	if (a.area == b.area)
	{
		size_t x = a.in.size() + a.de.size(), y = b.in.size() + b.de.size();
		return x > y;
	}
	return a.area > b.area;
}
PIN OverlapCal(size_t x1,size_t x2)
{
    PIN p;
    size_t de_r1 = BiCluster[x1].de.size(), de_r2 = BiCluster[x2].de.size();
	size_t c1 = BiCluster[x1].index_column.size(), c2 = BiCluster[x2].index_column.size();
	size_t in_r1 = BiCluster[x1].in.size(), in_r2 = BiCluster[x2].in.size();
	vector<uint8_t> flag1(n + 5, 0), flag2(m + 5, 0);
    p.second.de.reserve(de_r2),p.second.in.reserve(in_r2);
    p.second.index_column.reserve(c2);

	for (size_t i = 0; i < in_r1; i++) flag1[BiCluster[x1].in.at(i)] = 1;
	for (size_t i = 0; i < de_r1; i++) flag1[BiCluster[x1].de.at(i)] = 1;
	for (size_t i = 0; i < c1; i++) flag2[BiCluster[x1].index_column.at(i)] = 1;
	size_t r = 0, c = 0;
	for (size_t i = 0; i < de_r2; i++) 
        if (flag1[BiCluster[x2].de.at(i)]) r++;
        else p.second.de.emplace_back(BiCluster[x2].de.at(i));
    
	for (size_t i = 0; i < in_r2; i++) 
        if (flag1[BiCluster[x2].in.at(i)]) r++;
        else p.second.in.emplace_back(BiCluster[x2].in.at(i));
	for (size_t i = 0; i < c2; i++) 
        if (flag2[BiCluster[x2].index_column.at(i)]) c++;
        else p.second.index_column.emplace_back(BiCluster[x2].index_column.at(i));
	
    p.second.area=(p.second.de.size()+p.second.in.size())*p.second.index_column.size();

	if(c==0&&r==0)
	{
		p.first=0;
		return p;
	}
	if(c==0) c=1;
	if(r==0) r=1;
	p.first = c*r;
	return p;
}
void result_output(char* _out)
{
	ofstream output(_out);
	vector<uint8_t> flag(BiCluster_num,0);

	size_t cnt = 0;

    size_t il1,il2,icol,jl1,jl2,jcol,area1,area2;
	for (size_t i = 0; i < BiCluster_num; i++)
	{
		il1=BiCluster[i].in.size();
		il2=BiCluster[i].de.size();
		icol=BiCluster[i].index_column.size();
		area1=(il1+il2)*icol;
		BiCluster[i].area=area1;

		for (size_t j = i + 1; j < BiCluster_num; j++)
		{
            il1=BiCluster[i].in.size();
			il2=BiCluster[i].de.size();
			icol=BiCluster[i].index_column.size();
			area1=(il1+il2)*icol;
			BiCluster[i].area=area1;

            jl1=BiCluster[j].in.size();
            jl2=BiCluster[j].de.size();
            jcol=BiCluster[j].index_column.size();
			area2=(jl1+jl2)*jcol;
			BiCluster[j].area=area2;

			size_t min_id;PIN p;
			if(area1>=area2)
			{
				p = OverlapCal(i, j);
				min_id=j;
			}
			else 
			{
				p = OverlapCal(j, i);
				min_id=i;
			}

			if(po->FILTER<eps&&(p.first*1.0 / min(area1,area2)) > (po->FILTER))
			{
				BiCluster[min_id]=p.second;
			}
			else if (po->FILTER>eps&&(p.first*1.0 / min(area1,area2)) >= (po->FILTER))
			{
				
				BiCluster[min_id]=p.second;
				
			}
		}
	}
	
	sort(BiCluster.begin(), BiCluster.begin()+BiCluster_num, CMP_node);
    size_t row_len,col_len;
	for (size_t i = 0; i < BiCluster_num; i++)
	{
		il1=BiCluster[i].in.size()+BiCluster[i].de.size();
        icol= BiCluster[i].index_column.size();
		if(po->MIN_LENGTH>0)
			row_len=ceil(n*po->MIN_LENGTH),col_len=ceil(m*po->MIN_LENGTH);
		else
			row_len=col_len=0;
		row_len=max(po->ROW_WIDTH,row_len);
		col_len=max(po->COL_WIDTH,col_len);
		if (il1 <row_len || icol < col_len||(po->IS_SINGLE_CELL_DATA&&(il1<po->CLUSTER_WIDTH||icol<po->CLUSTER_WIDTH))) flag[i] = 1;
		else cnt++;
	}
	
	
	cnt = min(po->BLOCK_NUM, cnt);
	
	output<<"#Parameters: "<<"q:"<<po->QUANTILE<<" a:"<<po->ABSOLUTE_QUANTILE<<" s:"<<po->SINGLE_CELL_PROCESSING<<" z:"<<po->ZERO_CLEAN;
	output<<" l:"<<po->LCS_TOLERANCE<<" n:"<<po->CLUSTER_WIDTH<<" d:"<<po->DISCRETIZATION<<" t:"<<po->THREADS_NUM;
	output<<" e:"<<po->EXPAND_TOLERANCE<<" b:"<<po->DICHOTOMY_TOLERANCE;
	output<<" c:"<<po->COL_WIDTH<<" r:"<<po->ROW_WIDTH<<" m:"<<po->MIN_LENGTH<<" f:"<<po->FILTER;
	output<<" S:"<<po->IS_SINGLE_CELL_DATA<<" o:"<<po->BLOCK_NUM<<"\n";
	
	output << "BiCluster_Num:"<<cnt << "\n";
	size_t j=0;
	for (size_t i = 0; j<cnt; i++)
	{
		if(flag[i]) continue;
		il1 = BiCluster[i].in.size(),il2= BiCluster[i].de.size(),  icol= BiCluster[i].index_column.size();
		output<<"BC: "<<i<<"\n";
		output <<"PC_Genes ["<< il1<<"]: ";
		sort(BiCluster[i].in.begin(), BiCluster[i].in.end());
		for (size_t k = 0; k < il1; k++) output << genes.at(BiCluster[i].in.at(k)) << " ";
		output << "\n";
        output <<"NC_Genes ["<< il2<<"]: ";
        sort(BiCluster[i].de.begin(), BiCluster[i].de.end());
		for (size_t k = 0; k < il2; k++) output << genes.at(BiCluster[i].de.at(k)) << " ";
        output << "\n";

		output <<"Conds [" <<icol << "]: ";
		sort(BiCluster[i].index_column.begin(), BiCluster[i].index_column.end());
		for (size_t k = 0; k < icol; k++) output << conds.at(BiCluster[i].index_column.at(k)) << " ";
		output << "\n\n";
		j++;
	}
	uglyTime("%d bicluster are written to %s", cnt,_out);
	cout<<"All jobs completed"<<endl;
}
