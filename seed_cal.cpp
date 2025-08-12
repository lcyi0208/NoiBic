#include "structure.h"


size_t Upper_bound_in(size_t l, size_t r, double val, vector<double> tmp)
{
	size_t mid;
	while (l < r)
	{
		mid = (l + r) / 2;
		if (tmp[mid] <= val) l = mid + 1;
		else r = mid;
	}
	if (tmp[l] <= val) return (l + 1);
	return l;
}
size_t Lower_bound_in(size_t l, size_t r, double val, vector<double> tmp)
{
	size_t mid;
	while (l < r)
	{
		mid = (l + r) / 2;
		if (tmp[mid] < val) l = mid + 1;
		else r = mid;
	}
	if (tmp[l] < val) return (l + 1);
	return l;
}
int Upper_bound_de(size_t l, size_t r, double val, vector<double> tmp)
{
	size_t mid;
	while (l < r)
	{
		mid = ((l + r) >> 1) + 1;
		if (tmp[mid] >= val) l = mid;
		else r = mid - 1;
	}
	if (tmp[l] < val) return -1;
	return l;
}
int Lower_bound_de(size_t l, size_t r, double val, vector<double> tmp)
{
	size_t mid;
	while (l < r)
	{
		mid = ((l + r) >> 1) + 1;
		if (tmp[mid] > val) l = mid;
		else r = mid - 1;
	}
	if (tmp[l] <= val) return -1;
	return l;
}
size_t seed_cal(Node S, size_t col)
{
	size_t len1 = S.de.size(), len2 = S.in.size(), len = S.index_column.size();
	
	vector<size_t> sum(len+5,0);
	size_t M = 0;
	size_t low,up;
	
	for (size_t i = 0; i < len1; i++)
	{
		if (B[S.de[i]][col].val == -2) continue;
		vector<double> tmp;
		tmp.reserve(len);
		double Sum = 0, var = 0;
		size_t ll = floor(len / 4.0), rr = ceil(len*3 / 4.0);
		for (size_t j = 0; j < len; j++)
		{
			tmp.emplace_back(B[S.de.at(i)][S.index_column.at(j)].val);
			if (j >= ll && j < rr) Sum += B[S.de.at(i)][S.index_column.at(j)].val;
		}
		Sum = Sum / (rr - ll);
		for (size_t j = ll; j < rr; j++)
		{
			var += (tmp[j] - Sum)*(tmp[j] - Sum);
		}
		
		var = sqrt(var / (rr - ll));
		sort(tmp.begin(), tmp.end(), greater<double>());
		
		double change = fabs(var/Sum) * po->DICHOTOMY_TOLERANCE;
		double val_l, val_r;
		val_l = B[S.de.at(i)][col].val*(1 - change);
		val_r = B[S.de.at(i)][col].val*(1 + change);
		up = Upper_bound_de(0, S.index_column.size() - 1,val_l, tmp) + 1;
		low = Lower_bound_de(0, S.index_column.size() - 1, val_r, tmp) + 1;
		
		for (size_t j = low; j <= up; j++)
		{
			sum[j]++;
			if (sum[j] > M) M = sum[j];
		}
	}
	
	for (size_t i = 0; i < len2; i++)
	{
		if (B[S.in[i]][col].val == -2) continue;
		vector<double> tmp;
		double Sum = 0, var = 0;
		tmp.reserve(len);
		size_t ll = floor(len / 4.0), rr = ceil(len*3 / 4.0);
		for (size_t j = 0; j < len; j++)
		{
			tmp.emplace_back(B[S.in.at(i)][S.index_column.at(j)].val);
			if (j >= ll && j < rr) Sum += B[S.in.at(i)][S.index_column.at(j)].val;
		}
		Sum = Sum / (rr - ll);
		for (size_t j = ll; j < rr; j++)
		{
			var += (tmp[j] - Sum)*(tmp[j] - Sum);
		}
		var = sqrt(var / (rr - ll));
		sort(tmp.begin(), tmp.end());
		
		double change = fabs(var/Sum) * po->DICHOTOMY_TOLERANCE;
		double val_l, val_r;
		val_l = B[S.in.at(i)][col].val*(1 - change);
		val_r = B[S.in.at(i)][col].val*(1 + change);
		up = Upper_bound_in(0, S.index_column.size() - 1,val_r, tmp);
		low = Lower_bound_in(0, S.index_column.size() - 1, val_l, tmp);
		
		for (size_t j = low; j <= up; j++)
		{
			sum[j]++;
			if (sum[j] > M) M = sum[j];
		}
	}
	return M;
}
