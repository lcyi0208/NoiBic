#include "structure.h"
#include <boost/dynamic_bitset.hpp>
using namespace boost;
using namespace std;
#include <thread>
#include <mutex>

PIA CLUSTER(const size_t& i, const vector<Array>& a, const size_t& x1)
{
	const size_t len = a.size();
	vector<Array> b;
	b.reserve(len);
	for (size_t j = 0; j < len; ++j)
	{
		Array t;
		t.con = a[j].con;
		t.val = B[i][t.con].val;
		t.old_val = t.val;
		b.emplace_back(t);
	}

	std::sort(b.begin(), b.end());
	PIA tmp;
	tmp.id = i;
	tmp.ans = LCS(a, b, len, len, x1, i, 1);
	return tmp;
}
bool cmp(const PIA& a, const PIA& b)
{
	if (a.ans.len == b.ans.len)
		return a.id < b.id;
	return a.ans.len > b.ans.len;
}

size_t BiCluster_num = 0;
vector<Node> BiCluster;
namespace
{
	struct MaterializedSeed
	{
		size_t x1 = 0;
		size_t x2 = 0;
		size_t sig = 1;
		vector<size_t> cols;
	};

	struct RepeatInfo
	{
		double max_ratio = 0.0;
		vector<size_t> cols;
	};

	inline bool cluster_check_val(double a, double b, double ratio)
	{
		if (min(fabs(a), fabs(b)) < eps)
			return fabs(a - b) >= ratio * po->LCS_TOLERANCE;
		return (fabs(a - b) / min(fabs(b), fabs(a))) >= ratio * po->LCS_TOLERANCE;
	}

	static RepeatInfo repeat_info_on_sorted_cols(const vector<size_t>& cols, size_t row)
	{
		RepeatInfo info;
		const size_t len = cols.size();
		if (len == 0)
			return info;

		vector<Array> sorted_cols(len);
		for (size_t i = 0; i < len; ++i)
		{
			size_t col = cols[i];
			sorted_cols[i].con = col;
			sorted_cols[i].val = B[row][col].val;
			sorted_cols[i].old_val = sorted_cols[i].val;
		}
		sort(sorted_cols.begin(), sorted_cols.end());

		info.cols.reserve(len);
		for (const Array& item : sorted_cols)
			info.cols.emplace_back(item.con);

		size_t cnt = 1;
		for (size_t i = 1; i < len; ++i)
		{
			const double prev = sorted_cols[i - 1].val;
			const double cur = sorted_cols[i].val;
			const double ts = NoiseMatrix[row][sorted_cols[i].con];

			if (fabs(cur - prev) < eps || (!isnan(ts) && !cluster_check_val(cur, prev, ts*po->DISCRETIZATION)))
			{
				++cnt;
			}
			else
			{
				double ratio = cnt * 1.0 / len;
				info.max_ratio = max(info.max_ratio, ratio);
				cnt = 1;
			} 
		}

		double ratio = cnt * 1.0 / len;
		info.max_ratio = max(info.max_ratio, ratio);

		return info;
	}

	static MaterializedSeed materialize_seed(const SEED& S)
	{
		MaterializedSeed seed;
		seed.x1 = S.x1;
		seed.x2 = S.x2;

		Ans ans = LCS(A[S.x1], A[S.x2], m, m, S.x1, S.x2, 1);
		seed.sig = ans.sig;
		seed.cols.reserve(ans.len);
		for (size_t i = 1; i <= ans.len; ++i)
			seed.cols.emplace_back(ans.p[i]);

		return seed;
	}

	static void ensure_bicluster_slot(vector<Gene>& row_cluster)
	{
		if (BiCluster_num >= BiCluster.size())
			BiCluster.resize(max(BiCluster.size() * 2 + 1, BiCluster_num + 1));

		size_t need = BiCluster_num + 1;
		if (!row_cluster.empty() && row_cluster[0].cluster_index.size() < need)
		{
			size_t new_size = max(row_cluster[0].cluster_index.size() * 2 + 1, need);
			for (Gene& gene : row_cluster)
				gene.cluster_index.resize(new_size, 0);
		}
	}

	static void print_bicluster(size_t t)
	{
		cout << "BC:" << BiCluster_num << ": rowin:" << BiCluster[t].in.size();
		cout << " rowde:" << BiCluster[t].de.size() << " col:" << BiCluster[t].index_column.size() << endl;
		for (size_t i = 0; i < BiCluster[t].in.size(); i++)
			cout << BiCluster[t].in[i] << " ";
		cout << endl;
		for (size_t i = 0; i < BiCluster[t].de.size(); i++)
			cout << BiCluster[t].de[i] << " ";
		cout << endl;
		for (size_t i = 0; i < BiCluster[t].index_column.size(); i++)
			cout << BiCluster[t].index_column[i] << " ";
		cout << endl;

		for (size_t i = 0; i < BiCluster[t].in.size(); i++)
		{
			cout << genes[BiCluster[t].in[i]] << " ";
			sort(BiCluster[t].index_column.begin(), BiCluster[t].index_column.end(),
				[&](size_t c1, size_t c2)
				{
					double v1 = B[BiCluster[t].in[i]][c1].val;
					double v2 = B[BiCluster[t].in[i]][c2].val;
					if (v1 != v2)
						return v1 < v2;
					return c1 < c2;
				});
			for (size_t j = 0; j < BiCluster[t].index_column.size(); j++)
			{
				size_t col = BiCluster[t].index_column[j];
				cout << B[BiCluster[t].in[i]][col].con << "-" << B[BiCluster[t].in[i]][col].val << " ";
			}
			cout << endl;
		}
		for (size_t i = 0; i < BiCluster[t].de.size(); i++)
		{
			cout << genes[BiCluster[t].de[i]] << " ";
			sort(BiCluster[t].index_column.begin(), BiCluster[t].index_column.end(),
				[&](size_t c1, size_t c2)
				{
					double v1 = B[BiCluster[t].de[i]][c1].val;
					double v2 = B[BiCluster[t].de[i]][c2].val;
					if (v1 != v2)
						return v1 < v2;
					return c1 < c2;
				});
			for (size_t j = 0; j < BiCluster[t].index_column.size(); j++)
			{
				size_t col = BiCluster[t].index_column[j];
				cout << B[BiCluster[t].de[i]][col].con << "-" << B[BiCluster[t].de[i]][col].val << " ";
			}
			cout << endl;
		}
	}

	static void cluster_seed_part(const MaterializedSeed& seed, size_t left, size_t right,
		vector<Gene>& row_cluster, ThreadPool& pool)
	{
		if (right < left || right >= seed.cols.size())
			return;

		const size_t len1 = right - left + 1;
		if (len1 < po->COL_WIDTH)
			return;

		ensure_bicluster_slot(row_cluster);
		size_t t = BiCluster_num++;
		BiCluster[t] = Node{};
		BiCluster[t].index_column.assign(seed.cols.begin() + left, seed.cols.begin() + right + 1);

		size_t x1 = seed.x1;
		size_t x2 = seed.x2;
		vector<Array> a(len1);
		for (size_t i = 0; i < len1; i++)
		{
			size_t col = seed.cols[left + i];
			a[i].con = col;
			a[i].val = B[x1][col].val;
			a[i].old_val = a[i].val;
		}

		BiCluster[t].x1 = x1;
		BiCluster[t].x2 = x2;
		BiCluster[t].sig = seed.sig;
		BiCluster[t].in.emplace_back(x1);
		if (seed.sig)
			BiCluster[t].in.emplace_back(x2);
		else
			BiCluster[t].de.emplace_back(x2);

		row_cluster[x1].cluster_index[t] = 1;
		row_cluster[x2].cluster_index[t] = 1;
		row_cluster[x1].c.emplace_back(t);
		row_cluster[x2].c.emplace_back(t);

		size_t l1 = BiCluster[t].in.size();
		size_t l2 = BiCluster[t].de.size();
		size_t len2 = l1 + l2;

		vector<dynamic_bitset<>> bs(n, dynamic_bitset<>(m));
		vector<future<PIA>> results;
		results.reserve(n > 2 ? n - 2 : 0);
		vector<PIA> LCS_Matrix;
		LCS_Matrix.reserve(n);
		PIA empty{};
		empty.id = n;
		empty.ans.len = 0;
		LCS_Matrix.emplace_back(empty);

		for (size_t i = 0; i < n; i++)
		{
			if (row_cluster[i].cluster_index[t])
				continue;
			results.emplace_back(pool.enqueue([i, &a, x1]()
				{
					return CLUSTER(i, a, x1);
				}));
		}

		for (auto&& result : results)
		{
			PIA tmp = result.get();
			LCS_Matrix.emplace_back(tmp);
			dynamic_bitset<> db1(m);
			for (size_t j = 1; j <= tmp.ans.len; ++j)
				db1.set(tmp.ans.p[j]);
			bs[tmp.id] = std::move(db1);
		}

		std::sort(LCS_Matrix.begin(), LCS_Matrix.end(), cmp);

		dynamic_bitset<> temp(m);
		size_t current_col_len = len1;
		while (check(current_col_len, len2, LCS_Matrix[0].ans.len, len2 + 1, len1))
		{
			if (LCS_Matrix[0].ans.sig)
				BiCluster[t].in.emplace_back(LCS_Matrix[0].id);
			else
				BiCluster[t].de.emplace_back(LCS_Matrix[0].id);
			row_cluster[LCS_Matrix[0].id].c.emplace_back(t);
			row_cluster[LCS_Matrix[0].id].cluster_index[t] = 1;

			temp = bs[LCS_Matrix[0].id];
			current_col_len = LCS_Matrix[0].ans.len;
			size_t i = 1;
			while (i < LCS_Matrix.size() && LCS_Matrix[i].ans.len)
			{
				bs[LCS_Matrix[i].id] &= bs[LCS_Matrix[0].id];
				LCS_Matrix[i].ans.len = bs[LCS_Matrix[i].id].count();
				i++;
			}
			LCS_Matrix[0].ans.len = 0;
			std::sort(LCS_Matrix.begin(), LCS_Matrix.end(), cmp);
			l1 = BiCluster[t].in.size();
			l2 = BiCluster[t].de.size();
			len2 = l1 + l2;
		}

		if ((l1 + l2) > 2)
		{
			BiCluster[t].index_column.clear();
			for (size_t i = 0; i < m; ++i)
			{
				if (temp[i])
					BiCluster[t].index_column.emplace_back(i);
			}
		}

		//print_bicluster(t);
	}
	
	static void print_repeat_info(const char* label, const RepeatInfo& info, size_t row)
	{
		cout << label << " repeat:" << info.max_ratio << endl;
		cout << label << " cols:";
		for (size_t col : info.cols)
			cout << " " << col;
		cout << endl;
		cout << label << " vals:";
		for (size_t col : info.cols)
			cout << " " << B[row][col].val;
		cout << endl;
	}

	static void process_seed_parts(MaterializedSeed seed, vector<Gene>& row_cluster, ThreadPool& pool)
	{
		const size_t len = seed.cols.size();
		if (len < po->COL_WIDTH)
			return;

		RepeatInfo r1 = repeat_info_on_sorted_cols(seed.cols, seed.x1);
		RepeatInfo r2 = repeat_info_on_sorted_cols(seed.cols, seed.x2);
		//cout << "r1p: " << r1.max_ratio << " r2p: " << r2.max_ratio << " t: " << threshold << " "<<len*threshold<<endl;
		//print_repeat_info("r1", r1, seed.x1);
		//print_repeat_info("r2", r2, seed.x2);

		if (r2.max_ratio < r1.max_ratio)
		{
			swap(seed.x1, seed.x2);
			seed.cols = std::move(r2.cols);
		}
		else
		{
			seed.cols = std::move(r1.cols);
		}

		cluster_seed_part(seed, 0, len - 1, row_cluster, pool);
	}
}
PII overlapping_check(size_t x1, size_t x2)
{
	PII p;
	size_t de_r1 = BiCluster[x1].de.size(), de_r2 = BiCluster[x2].de.size(), c1 = BiCluster[x1].index_column.size(), c2 = BiCluster[x2].index_column.size();
	size_t in_r1 = BiCluster[x1].in.size(), in_r2 = BiCluster[x2].in.size();
	vector<uint8_t> flag1(n, 0), flag2(m, 0);
	for (size_t i = 0; i < in_r1; i++)
		flag1[BiCluster[x1].in.at(i)] = 1;
	for (size_t i = 0; i < de_r1; i++)
		flag1[BiCluster[x1].de.at(i)] = 1;
	for (size_t i = 0; i < c1; i++)
		flag2[BiCluster[x1].index_column.at(i)] = 1;
	size_t r = 0, c = 0;
	for (size_t i = 0; i < de_r2; i++)
		if (flag1[BiCluster[x2].de.at(i)])
			r++;
	for (size_t i = 0; i < in_r2; i++)
		if (flag1[BiCluster[x2].in.at(i)])
			r++;
	for (size_t i = 0; i < c2; i++)
		if (flag2[BiCluster[x2].index_column.at(i)])
			c++;
	p.first = c, p.second = r;
	return p;
}
bool bicluster_check(size_t x1, size_t x2, const vector<Gene>& row_cluster)
{
	size_t len1 = row_cluster[x1].c.size(), len2 = row_cluster[x2].c.size();
	for (size_t i = 0; i < len1; i++)
	{
		for (size_t j = 0; j < len2; j++)
		{
			if (row_cluster[x1].c[i] == row_cluster[x2].c[j])
				return false;
			PII p = overlapping_check(row_cluster[x1].c[i], row_cluster[x2].c[j]);
			if (p.second >= 2)
				return false;
		}
	}
	return true;
}
size_t Max(size_t x1, size_t x2, const vector<Gene>& row_cluster)
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
bool is_seed(const SEED& S, const vector<Gene>& row_cluster)
{
	size_t x1 = row_cluster[S.x1].c.size(), x2 = row_cluster[S.x2].c.size(), len = S.len;
	if ((!x1) || (!x2))
	{
		return true;
	}
	else
	{
		if (bicluster_check(S.x1, S.x2, row_cluster) && len >= Max(S.x1, S.x2, row_cluster))
			return true;
	}
	return false;
}
bool check(size_t a, size_t b, size_t c, size_t d, size_t Length)
{
	if (c < po->CLUSTER_WIDTH || c < Length * po->CLUSTER_SIZE)
		return false;
	if ((a >= Length * po->CLUSTER_SIZE * 3) && po->IS_SINGLE_CELL_DATA)
		return true;
	if (min(a, b) <= min(c, d))
		return true;
	return false;
}
void cluster()
{
	vector<Gene> row_cluster(n);
	size_t initial_slots = max<size_t>(Seed.size() * 2, 1);
	for (size_t i = 0; i < n; i++)
		row_cluster[i].cluster_index = vector<uint8_t>(initial_slots, 0);

	BiCluster_num = 0;
	BiCluster.assign(initial_slots, Node{});
	//cout << "seed num:" << Seed.size() << endl;
	ThreadPool pool(max<size_t>(po->THREADS_NUM, 1));

	while (!Seed.empty())
	{
		SEED S = Seed.top();
		//cout << "seed:" << S.x1 << " " << S.x2 << " " << S.len << endl;
		Seed.pop();
		if (!is_seed(S, row_cluster))
			continue;

		MaterializedSeed seed = materialize_seed(S);
		process_seed_parts(std::move(seed), row_cluster, pool);
	}

	uglyTime("%zu bicluster generated after clustering", BiCluster_num);
}
