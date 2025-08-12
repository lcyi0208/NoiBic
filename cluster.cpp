#include "structure.h"
#include<boost/dynamic_bitset/dynamic_bitset.hpp>
using namespace boost;
#include <thread>
#include <mutex>

PIA CLUSTER(const size_t &i,const vector<Array> &a)
{
	size_t len=a.size();
	vector<Array> b(len);
	for (size_t j = 0; j < len; j++)
	{
		b[j].con = a[j].con;
		b[j].val = B[i][b[j].con].val;
		b[j].old_val = b[j].val;
	}

	sort(b.begin(), b.end());
	PIA tmp;
	tmp.id=i;
	tmp.ans = LCS(a, b, len, len, 1);
	
	return tmp;
		
}
bool cmp(const PIA &a, const PIA &b)
{
	if (a.ans.len == b.ans.len) return a.id < b.id;
	return a.ans.len > b.ans.len;
}

size_t BiCluster_num = 0;
vector<Node> BiCluster;

PII overlapping_check(size_t x1, size_t x2)
{
	PII p;
	size_t de_r1 = BiCluster[x1].de.size(), de_r2 = BiCluster[x2].de.size(), c1 = BiCluster[x1].index_column.size(), c2 = BiCluster[x2].index_column.size();
	size_t in_r1 = BiCluster[x1].in.size(), in_r2 = BiCluster[x2].in.size();
	vector<uint8_t> flag1(n + 5, 0), flag2(m + 5, 0);
	for (size_t i = 0; i < in_r1; i++) flag1[BiCluster[x1].in.at(i)] = 1;
	for (size_t i = 0; i < de_r1; i++) flag1[BiCluster[x1].de.at(i)] = 1;
	for (size_t i = 0; i < c1; i++) flag2[BiCluster[x1].index_column.at(i)] = 1;
	size_t r = 0, c = 0;
	for (size_t i = 0; i < de_r2; i++) if (flag1[BiCluster[x2].de.at(i)]) r++;
	for (size_t i = 0; i < in_r2; i++) if (flag1[BiCluster[x2].in.at(i)]) r++;
	for (size_t i = 0; i < c2; i++) if (flag2[BiCluster[x2].index_column.at(i)]) c++;
	p.first = c, p.second = r;
	return p;
}
bool bicluster_check(size_t x1, size_t x2,const vector<Gene> &row_cluster)
{
	size_t len1 = row_cluster[x1].c.size(), len2 = row_cluster[x2].c.size();
	for (size_t i = 0; i < len1; i++)
	{
		for (size_t j = 0; j < len2; j++)
		{
			if (row_cluster[x1].c[i] == row_cluster[x2].c[j]) return false;
			PII p = overlapping_check(row_cluster[x1].c[i], row_cluster[x2].c[j]);
			if (p.second >= 2) return false;
		}
	}
	return true;
}
size_t Max(size_t x1, size_t x2,const vector<Gene> &row_cluster)
{
	size_t len1 = row_cluster[x1].c.size(), len2 = row_cluster[x2].c.size();
	size_t ans = 0;
	for (size_t i = 0; i < len1; i++)
	{
		size_t len = BiCluster[row_cluster[x1].c[i]].in.size() + BiCluster[row_cluster[x1].c[i]].de.size();
		ans = max(len, ans);
	}
	for (size_t j = 0; j < len2; j++)
	{
		size_t len = BiCluster[row_cluster[x2].c[j]].in.size() + BiCluster[row_cluster[x2].c[j]].de.size();
		ans = max(len, ans);
	}
	return ans;
}
bool is_seed(const SEED &S,const vector<Gene> &row_cluster)
{
	size_t x1 = row_cluster[S.x1].c.size(), x2 = row_cluster[S.x2].c.size(), len = S.len;
	if ((!x1) || (!x2))
	{
		return true;
	}
	else
	{
		if (bicluster_check(S.x1, S.x2,row_cluster) && len >= Max(S.x1, S.x2,row_cluster))
			return true;
	}
	return false;
}
bool check(size_t a, size_t b, size_t c, size_t d,size_t Length)
{ 
	if(c<po->CLUSTER_WIDTH||c<Length*po->CLUSTER_SIZE) return false;
	if (a>=Length*po->CLUSTER_SIZE*3||min(a, b) <= min(c, d)) return true;
	return false;
}
void cluster()
{
	vector<Gene> row_cluster(n);
	for (size_t i = 0; i < n; i++)
	{
		row_cluster[i].cluster_index = vector<uint8_t>(po->SEED_NUM * 2, 0);
	}
	BiCluster = vector<Node>(Seed.size() * 2);
	while (!Seed.empty())
	{
		SEED S = Seed.top();
		Seed.pop();
		if (!is_seed(S,row_cluster)) continue;
		Ans ans=LCS(A[S.x1],A[S.x2],m,m,1);
		
		size_t t = BiCluster_num;
		BiCluster_num++;
		size_t x1,x2;
		x1 = S.x1, x2 = S.x2;
		BiCluster[t].x1=x1,BiCluster[t].x2=x2;
		BiCluster[t].in.emplace_back(x1);
		if(ans.sig) BiCluster[t].in.emplace_back(x2);
		else BiCluster[t].de.emplace_back(x2);

		BiCluster[t].index_column.reserve(ans.len);
		for (size_t h = 1; h <= ans.len; h++)
		{
			BiCluster[t].index_column.emplace_back(ans.p[h]);
		}
		size_t len1 = ans.len, len2;
		size_t l1 = BiCluster[t].in.size(), l2 = BiCluster[t].de.size();
		len2=l1+l2;

		vector<Array> a(len1);
		row_cluster[x1].cluster_index[t] = 1;
		row_cluster[x2].cluster_index[t] = 1;
		row_cluster[x1].c.emplace_back(t);
		row_cluster[x2].c.emplace_back(t);
		for (size_t i = 0; i < len1; i++)
		{
			a[i].con = BiCluster[t].index_column[i];
			a[i].val = B[x1][BiCluster[t].index_column[i]].val;
			a[i].old_val = a[i].val;
		}
		sort(a.begin(), a.end());

		vector<dynamic_bitset<> >  bs(n, dynamic_bitset<>(m));
		ThreadPool pool(po->THREADS_NUM);
		vector<future<PIA> > results;
		results.reserve(n-2);
		vector<PIA> LCS_Matrix;
		LCS_Matrix.reserve(n-1);
		PIA tmp;
		tmp.ans.len=0;
		LCS_Matrix.emplace_back(tmp);
		

		for (size_t i = 0; i < n; i++)
		{
			if (row_cluster[i].cluster_index[t]) continue;
			results.emplace_back(pool.enqueue(CLUSTER,i,a));
			
		}
		
		for (auto&& result : results)
		{
			PIA tmp=result.get();
			LCS_Matrix.emplace_back(tmp);
			string s(m,'0');
			for (size_t j = 1; j <= tmp.ans.len; j++)
			{
				s[tmp.ans.p[j]]='1';
			}
			dynamic_bitset<> db1(s);
			bs.at(tmp.id)= db1;
		}
		
		
		sort(LCS_Matrix.begin(),LCS_Matrix.end(),cmp);
		
		dynamic_bitset<> temp(m);
		while (check(len1, len2, LCS_Matrix[0].ans.len, len2 + 1,ans.len))
		{
			if (LCS_Matrix[0].ans.sig) BiCluster[t].in.emplace_back(LCS_Matrix[0].id);
			else BiCluster[t].de.emplace_back(LCS_Matrix[0].id);
			row_cluster[LCS_Matrix[0].id].c.emplace_back(t);
			row_cluster[LCS_Matrix[0].id].cluster_index[t] = 1;

			temp=bs[LCS_Matrix[0].id];
			len1=LCS_Matrix[0].ans.len;
			size_t i = 1;
			while (LCS_Matrix[i].ans.len)
			{
				bs[LCS_Matrix[i].id] &= bs[LCS_Matrix[0].id];
				LCS_Matrix[i].ans.len = bs[LCS_Matrix[i].id].count();
				i++;
			}
			LCS_Matrix[0].ans.len = 0;
			sort(LCS_Matrix.begin(),LCS_Matrix.end(), cmp);
			l1 = BiCluster[t].in.size(), l2 = BiCluster[t].de.size();
			len2 = l1 + l2;
		}
		if((l1+l2)>2)
		{
			BiCluster[t].index_column.clear();
			for (int i = m - 1; i >= 0; i--)
			{
				if (temp[i] & 1) BiCluster[t].index_column.emplace_back(m - 1 - i);
			}
		}
	}
	
	uglyTime("%d bicluster generated after clustering",BiCluster_num);
}
