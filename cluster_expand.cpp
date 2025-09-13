#include "structure.h"

void column_expand()
{
	for (size_t i = 0; i < BiCluster_num; i++)
	{
		vector<uint8_t> flag(m + 5, 0);
		size_t len = BiCluster[i].index_column.size();
		size_t l1 = BiCluster[i].in.size(), l2 = BiCluster[i].de.size();
		if (len < po->COL_WIDTH / 2 || (l1 + l2) < po->COL_WIDTH/2)
			continue;


		for (size_t j = 0; j < len; j++)
			flag[BiCluster[i].index_column.at(j)] = 1;
		for (size_t j = 0; j < m; j++)
		{
			if (flag[j])
				continue;
			size_t ans = seed_cal(BiCluster[i], j);
			if (ans >= (l1 + l2) * (po->EXPAND_TOLERANCE))
			{
				BiCluster[i].index_column.emplace_back(j);
			}
		}
	}
}
PSS CLUSTER_EXPAND(const size_t &j, const vector<vector<Array> > &a)
{
	size_t len = a[0].size();
	vector<Array> b(len);
	for (size_t k = 0; k < len; k++)
	{
		b[k].con = a[0][k].con;
		b[k].val = B[j][b[k].con].val;
		b[k].old_val = b[k].val;
	}

	sort(b.begin(), b.end());
	Ans ans = LCS(a[0], b, len, len, 0);
	
	SEED tmp;
	tmp.x1 = j, tmp.x2 = ans.sig, tmp.len = ans.len;
	PSS P;
	P.first = tmp;

	ans = LCS(a[1], b, len, len, 0);
	tmp.x1 = j, tmp.x2 = ans.sig, tmp.len = ans.len;
	P.second = tmp;
	return P;
}

void row_expand()
{
	size_t tmp = BiCluster_num;
	
	for (size_t i = 0; i < tmp; i++)
	{
		if (BiCluster[i].index_column.size() < po->COL_WIDTH / 2)
			continue;
		size_t len2 = BiCluster[i].index_column.size();
		size_t l1 = BiCluster[i].in.size(), l2 = BiCluster[i].de.size();
		size_t t = BiCluster_num;
		vector<vector<Array>> a(2, vector<Array>(len2));
		BiCluster_num++;
		

		for (size_t j = 0; j < len2; j++)
		{
			a[0][j].con = BiCluster[i].index_column.at(j);
			a[0][j].val = B[BiCluster[i].x1][a[0][j].con].val;
			a[0][j].old_val = a[0][j].val;
			a[1][j].con = BiCluster[i].index_column.at(j);
			a[1][j].val = B[BiCluster[i].x2][a[1][j].con].val;
			a[1][j].old_val = a[1][j].val;
		}

		sort(a[0].begin(), a[0].end());
		sort(a[1].begin(), a[1].end());
	
		vector<uint8_t> flag(n + 5, 0);
		if (po->SINGLE_CELL_PROCESSING || po->ZERO_CLEAN||po->QUANTILE<0.2)
		{
			for (size_t j = 0; j < l1; j++)
				flag[BiCluster[i].in[j]] = 1;
			for (size_t j = 0; j < l2; j++)
				flag[BiCluster[i].de[j]] = 1;
		}
		else
		{
			BiCluster[i].de.clear();
			BiCluster[i].in.clear();
		}
		BiCluster[t] = BiCluster[i];
		ThreadPool pool(po->THREADS_NUM);
		vector<future<PSS> > results;
		results.reserve(n);
		for (size_t j = 0; j < n; j++)
		{
			if (flag[j])
				continue;
			results.emplace_back(pool.enqueue(CLUSTER_EXPAND, j, a));
		}
		
		for (auto &&result : results)
		{
			PSS tmp = result.get();
			if (tmp.first.len >= (po->EXPAND_TOLERANCE) * len2)
			{
				if (tmp.first.x2)
					BiCluster[i].in.emplace_back(tmp.first.x1);
				else
					BiCluster[i].de.emplace_back(tmp.first.x1);
			}

			if (tmp.second.len >= (po->EXPAND_TOLERANCE) * len2)
			{
				if (tmp.second.x2)
					BiCluster[t].in.emplace_back(tmp.second.x1);
				else
					BiCluster[t].de.emplace_back(tmp.second.x1);
			}
		}
		
		if (BiCluster[i].de.size() + BiCluster[i].in.size() < 2)
			BiCluster[i] = BiCluster[t];
		if (BiCluster[t].de.size() + BiCluster[t].in.size() < 2)
		{
			BiCluster_num--;
			BiCluster[t].de.clear();
			BiCluster[t].in.clear();
			BiCluster[t].index_column.clear();
		}
		else 
		{
			l1 = BiCluster[t].in.size(), l2 = BiCluster[t].de.size();
			len2 = BiCluster[t].index_column.size();
		}
	}
}
void cluster_expand()
{
	
	if(po->DICHOTOMY_TOLERANCE!=2)
		column_expand();
	row_expand();
	
	uglyTime("%d bicluster generated after expanding", BiCluster_num);
}
// Initial commit (no-op)
