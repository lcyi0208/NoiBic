#include "structure.h"
#include <random>
#include<cmath>

struct min_max
{
	double Min, Max;
};
vector<min_max> mm;
vector<size_t> Start_pos;

size_t Upper_bound_Array(size_t l,size_t r,double val,const vector<Array> &a)
{
    size_t mid;
    while (l < r)
    {
        mid = (l + r) / 2;
        if (a[mid].val <= val) l = mid + 1;
        else r = mid;
    }
    if (a[l].val <= val) return (l + 1);
    return l;
}
size_t Lower_bound_Array(size_t l,size_t r,double val,const vector<Array> &a)
{
	size_t mid;
	while (l < r)
	{
		mid = (l + r) / 2;
		if (a[mid].val < val) l = mid + 1;
		else r = mid;
	}
	
	if (a[l].val < val) return (l + 1);
	return l;
}
inline double mapToRange(double value,double mean,double stddev,double minValue,double maxValue)
{
    double scaledValue = mean + stddev * value;
	return max(minValue, min(maxValue, scaledValue));
}
bool check_val(const double& val,const double& d,const double& mid,const size_t& j,const size_t& s,const size_t& t)
{
	if(val==-INF)
		return true;
	if(po->ABSOLUTE_QUANTILE==1)
	{
		if(j>s&&j<t) return true;
		else return false;
	}
	if(val < mid+d && val > (mid - d))
		return true; 
	if(j>s&&j<t)
	{
		if(fabs(val-(mid+d))<eps||fabs(val-(mid-d))<eps)
			return true;
	}
	return false;
}
void data_preprocessing()
{
	mm = vector<min_max>(n + 5);
	
	for (size_t i = 0; i < n; i++)
	{
		sort(A[i].begin(), A[i].end());
		mm[i].Min = A[i][0].val;
		mm[i].Max = A[i][m - 1].val;
	}

	if(po->SINGLE_CELL_PROCESSING)
	{
		random_device rd;
    	mt19937 gen(rd());  
		double stddev,mean_val,random_value,minValue,maxValue;
		for(size_t i=0;i<n;i++)
		{
			size_t up_index=Upper_bound_Array(0,m-1,0,A[i]);
			size_t low_index=Lower_bound_Array(0,m-1,0,A[i]);
			if(up_index<m) maxValue=A[i][up_index].val;
			else maxValue=A[i][m-1].val;
			if(low_index>0) minValue=A[i][low_index-1].val;
			else minValue=A[i][0].val;

			normal_distribution<double> normalDist(0.0, 1.0);
			stddev=(maxValue-minValue)/8.0;
			mean_val=(maxValue+minValue)/2.0;
			for(size_t j=low_index;j<up_index;j++)
			{
				random_value=normalDist(gen);
				A[i][j].val=mapToRange(random_value, mean_val, stddev, minValue, maxValue);
				A[i][j].old_val=A[i][j].val;
			}
			
			sort(A[i].begin(), A[i].end());
			mm[i].Min=A[i][0].val;
			mm[i].Max=A[i][m-1].val;
		}
	}
	if(po->ZERO_CLEAN)
	{
		for(size_t i=0;i<n;i++)
		{
			size_t up_index=Upper_bound_Array(0,m-1,0,A[i]);
			size_t low_index=Lower_bound_Array(0,m-1,0,A[i]);
			for(size_t j=low_index;j<up_index;j++)
			{
				A[i][j].val=-INF;
				A[i][j].old_val=A[i][j].val;
			}
			sort(A[i].begin(), A[i].end());
		}
	}
	
	Start_pos=vector<size_t> (n+5);
	size_t count=0;
	if ((po->QUANTILE) < 0.5)
	{
		size_t l = m / 2;
		size_t s = (size_t)((po->QUANTILE)*m), t = (size_t)((1 - (po->QUANTILE))*m);
		double left, right, d, mid;
		
		for (size_t i = 0; i < n; i++)
		{
			left = A[i][s].val, right = A[i][t].val;
			d = min((A[i][l].val - A[i][s].val), (A[i][t].val - A[i][l].val));
			mid = A[i][l].val;
			
			for (size_t j = 0; j < m; j++)
			{
				if (check_val(A[i][j].val,d,mid,j,s,t))
				{
					count++;
					
					A[i][j].val = -2;
					A[i][j].old_val = -2;
					B[i][A[i][j].con].val = -2;
					B[i][A[i][j].con].old_val = -2;
				
				}
				else
				{
					double tmp_val;
					if(fabs(mm[i].Max - mm[i].Min)<eps)
						tmp_val=1;
					else tmp_val=mm[i].Max - mm[i].Min;
					A[i][j].val = (A[i][j].val - mm[i].Min) / tmp_val;
					A[i][j].old_val = A[i][j].val;
					B[i][A[i][j].con].val = A[i][j].val;
					B[i][A[i][j].con].old_val = A[i][j].val;
				}
			}
			Start_pos[i]=count;
			count=0;
			sort(A[i].begin(), A[i].end());
		}
	}
	else
	{
		for (size_t i = 0; i < n; i++)
		{
			for (size_t j = 0; j < m; j++)
			{
				if (A[i][j].val==-INF)
				{
					count++;
					A[i][j].val = -2;
					A[i][j].old_val = -2;
					B[i][A[i][j].con].val = -2;
					B[i][A[i][j].con].old_val = -2;
				
				}
				else
				{
					double tmp_val;
					if(fabs(mm[i].Max - mm[i].Min)<eps)
						tmp_val=1;
					else tmp_val=mm[i].Max - mm[i].Min;
					A[i][j].val = (A[i][j].val - mm[i].Min) / tmp_val;
					A[i][j].old_val = A[i][j].val;
					B[i][A[i][j].con].val = A[i][j].val;
					B[i][A[i][j].con].old_val = A[i][j].val;
				}
			}
			Start_pos[i]=count;
			count=0;
			sort(A[i].begin(), A[i].end());
		}
	}
	
	if(m>2500)
		po->CLUSTER_WIDTH=max((size_t)(m*po->QUANTILE*2)/20+3,po->CLUSTER_WIDTH);
	else
		po->CLUSTER_WIDTH=max(m/35+3,po->CLUSTER_WIDTH);
	
	//cout<<po->CLUSTER_WIDTH<<endl;

}
// Initial commit (no-op)
