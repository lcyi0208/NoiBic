#include "structure.h"
using namespace std;
namespace
{
	static vector<double> collect_column_values(const Node &bic, size_t col)
	{
		vector<double> vals;
		vals.reserve(bic.in.size() + bic.de.size());

		for (size_t i = 0; i < bic.in.size(); ++i)
		{
			double x = B[bic.in[i]][col].val;
			if (!isnan(x))
				vals.push_back(x);
		}
		for (size_t i = 0; i < bic.de.size(); ++i)
		{
			double x = B[bic.de[i]][col].val;
			if (!isnan(x))
				vals.push_back(x);
		}
		return vals;
	}

	static double column_cv_meanAbs_all(const Node &bic, size_t col)
	{
		vector<double> vals = collect_column_values(bic, col);
		if (vals.empty())
			return 1e100; 
		if (vals.size() == 1)
			return 0.0;

		double mean = 0.0;
		double meanAbs = 0.0;
		for (size_t i = 0; i < vals.size(); ++i)
		{
			mean += vals[i];
			meanAbs += fabs(vals[i]);
		}
		mean /= (double)vals.size();
		meanAbs /= (double)vals.size();

		double var = 0.0;
		for (size_t i = 0; i < vals.size(); ++i)
		{
			double d = vals[i] - mean;
			var += d * d;
		}
		var /= (double)vals.size();

		double sd = sqrt(var);
		double eps = 1e-12;
		return sd / (meanAbs + eps);
	}
	static void rebuild_existing_columns(Node &bic)
	{
		const size_t row_count = bic.in.size() + bic.de.size();
		if (bic.index_column.size() <= 1 || row_count == 0)
			return;

		vector<size_t> ordered_cols = bic.index_column;
		sort(ordered_cols.begin(), ordered_cols.end(),
			 [&](size_t c1, size_t c2)
			 {
				 double v1 = B[bic.x1][c1].val;
				 double v2 = B[bic.x1][c2].val;
				 const bool n1 = isnan(v1);
				 const bool n2 = isnan(v2);
				 if (n1 != n2)
					 return !n1;
				 if (!n1 && v1 != v2)
					 return v1 < v2;
				 return c1 < c2;
			 });

		const size_t missing_pos = numeric_limits<size_t>::max();
		vector<size_t> active_cols = ordered_cols;
		vector<size_t> active_pos(m, missing_pos);
		for (size_t i = 0; i < active_cols.size(); ++i)
			active_pos[active_cols[i]] = i;

		Node probe = bic;
		const double threshold = row_count * po->EXPAND_TOLERANCE;
		for (size_t candidate : ordered_cols)
		{
			size_t pos = active_pos[candidate];
			if (pos == missing_pos)
				continue;

			size_t last_col = active_cols.back();
			active_cols[pos] = last_col;
			active_pos[last_col] = pos;
			active_cols.pop_back();
			active_pos[candidate] = missing_pos;

			if (active_cols.empty())
			{
				active_pos[candidate] = active_cols.size();
				active_cols.emplace_back(candidate);
				continue;
			}

			probe.index_column = std::move(active_cols);
			size_t ans = seed_cal(probe, candidate);
			active_cols = std::move(probe.index_column);
			if (ans >= threshold)
			{
				active_pos[candidate] = active_cols.size();
				active_cols.emplace_back(candidate);
			}
		}

		Node rebuilt_bic = bic;
		rebuilt_bic.index_column.clear();
		rebuilt_bic.index_column.reserve(ordered_cols.size());
		for (size_t col : ordered_cols)
		{
			if (active_pos[col] != missing_pos)
				rebuilt_bic.index_column.emplace_back(col);
		}

		bic = std::move(rebuilt_bic);
	}
}

void column_expand()
{
	for (size_t i = 0; i < BiCluster_num; i++)
	{
		size_t len = BiCluster[i].index_column.size();
		size_t l1 = BiCluster[i].in.size();
		size_t l2 = BiCluster[i].de.size();

		if (len < po->COL_WIDTH / 2 || (l1 + l2) < po->COL_WIDTH / 2)
			continue;

		vector<uint8_t> flag(m, 0);
		for (size_t col : BiCluster[i].index_column)
			flag[col] = 1;

		rebuild_existing_columns(BiCluster[i]);
		len = BiCluster[i].index_column.size();
		if (len < po->COL_WIDTH / 2)
			continue;

		/*cout << "before col expand:" << endl;
		cout << "Bicluster:" << i
			 << ": rowin:" << l1
			 << " rowde:" << l2
			 << " column:" << len << endl;

		for (int j = 0; j < (int)l1; j++)
		{
			sort(BiCluster[i].index_column.begin(), BiCluster[i].index_column.end(),
				 [&](size_t c1, size_t c2)
				 {
					 double v1 = B[BiCluster[i].in[j]][c1].val;
					 double v2 = B[BiCluster[i].in[j]][c2].val;
					 if (v1 != v2)
						 return v1 < v2;
					 return c1 < c2;
				 });

			cout << genes[BiCluster[i].in[j]] << ":\t";
			for (int k = 0; k < (int)len; k++)
			{
				cout << B[BiCluster[i].in[j]][BiCluster[i].index_column[k]].con << "-";
				cout << B[BiCluster[i].in[j]][BiCluster[i].index_column[k]].val << "\t";
			}
			cout << endl;
		}
		for (int j = 0; j < (int)l2; j++)
		{
			sort(BiCluster[i].index_column.begin(), BiCluster[i].index_column.end(),
				 [&](size_t c1, size_t c2)
				 {
					 double v1 = B[BiCluster[i].de[j]][c1].val;
					 double v2 = B[BiCluster[i].de[j]][c2].val;
					 if (v1 != v2)
						 return v1 < v2;
					 return c1 < c2;
				 });
			cout << genes[BiCluster[i].de[j]] << ":\t";
			for (int k = 0; k < (int)len; k++)
			{
				cout << B[BiCluster[i].de[j]][BiCluster[i].index_column[k]].con << "-";
				cout << B[BiCluster[i].de[j]][BiCluster[i].index_column[k]].val << "\t";
			}
			cout << endl;
		}
		cout << endl;*/

		vector<pair<double, size_t>> cand_cols;
		cand_cols.reserve(m);

		for (size_t j = 0; j < m; ++j)
		{
			if (flag[j])
				continue;

			double cv = column_cv_meanAbs_all(BiCluster[i], j);
			cand_cols.push_back(make_pair(cv, j));
		}

		sort(cand_cols.begin(), cand_cols.end(),
			 [](const pair<double, size_t> &a, const pair<double, size_t> &b)
			 {
				 if (a.first != b.first)
					 return a.first < b.first; 
				 return a.second < b.second;
			 });

		/*cout << "candidate columns sorted by cv:" << endl;
		for (size_t t = 0; t < cand_cols.size(); ++t)
		{
			cout << "rank " << t
				 << " col " << cand_cols[t].second
				 << " " << conds[cand_cols[t].second]
				 << " cv " << cand_cols[t].first << endl;
		}*/

		for (size_t t = 0; t < cand_cols.size(); ++t)
		{
			size_t col = cand_cols[t].second;
			double cv = cand_cols[t].first;

			size_t ans = seed_cal(BiCluster[i], col);

			/*cout << "try col: " << col
				 << " " << conds[col]
				 << " cv: " << cv
				 << " ans: " << ans
				 << " threshold: " << (l1 + l2) * (po->EXPAND_TOLERANCE)
				 << endl;*/

			if (ans >= (l1 + l2) * (po->EXPAND_TOLERANCE))
			{
				// cout << "add col " << col << endl;
				BiCluster[i].index_column.emplace_back(col);
			}
			else
			{
				// cout << "reject col " << col << endl;
			}
		}
		/*cout << "after col expand: " << BiCluster[i].in.size() << " ";
		cout << BiCluster[i].de.size() << " " << BiCluster[i].index_column.size() << endl;
		for (int j = 0; j < BiCluster[i].in.size(); j++)
			cout << BiCluster[i].in[j] << " ";
		cout << endl;
		for (int j = 0; j < BiCluster[i].de.size(); j++)
			cout << BiCluster[i].de[j] << " ";
		cout << endl;
		for (int j = 0; j < BiCluster[i].index_column.size(); j++)
			cout << BiCluster[i].index_column[j] << " ";
		cout << endl;*/
	}
}
inline bool is_right_val(double a, double b, double ratio)
{
	if (min(fabs(a), fabs(b)) < eps)
		return fabs(a - b) >= ratio * po->LCS_TOLERANCE;
	return (fabs(a - b) / min(fabs(b), fabs(a))) >= ratio * po->LCS_TOLERANCE;
}
static bool anchor_repeat_ok(const vector<Array> &a, size_t row_id)
{
	const size_t len = a.size();
	if (len == 0)
		return false;
	if (len == 1)
		return true;

	size_t max_cnt = 1, cnt = 1;
	double pre = a[0].val;

	for (size_t i = 1; i < len; ++i)
	{
		const double cur = a[i].val;
		const double ts = NoiseMatrix[row_id][a[i].con];

		if (isnan(ts) || is_right_val(cur, pre, ts * po->DISCRETIZATION))
		{
			max_cnt = max(max_cnt, cnt);
			cnt = 1;
			pre = cur;
		}
		else
		{
			++cnt;
		}
	}

	max_cnt = max(max_cnt, cnt);
	// cout << (max_cnt * 1.0 / len) << endl;
	return (max_cnt * 1.0 / len) <= po->DISCRETIZATION;
}

PSS CLUSTER_EXPAND(const size_t &j,
				   const vector<vector<Array>> &a,
				   const size_t &x1,
				   const size_t &x2,
				   bool use_x2)
{
	const size_t len = a[0].size();
	vector<Array> b(len);

	for (size_t k = 0; k < len; ++k)
	{
		b[k].con = a[0][k].con;
		b[k].val = B[j][b[k].con].val;
		b[k].old_val = b[k].val;
	}

	sort(b.begin(), b.end());

	PSS P;
	P.first.len = 0;
	P.first.x1 = j;
	P.first.x2 = 0;
	P.second = P.first;

	Ans ans = LCS(a[0], b, len, len, x1, j, 0);
	P.first.x1 = j;
	P.first.x2 = ans.sig;
	P.first.len = ans.len;

	if (use_x2)
	{
		ans = LCS(a[1], b, len, len, x2, j, 0);
		P.second.x1 = j;
		P.second.x2 = ans.sig;
		P.second.len = ans.len;
	}

	return P;
}

void row_expand(ThreadPool &pool)
{
	size_t old_num = BiCluster_num;

	for (size_t i = 0; i < old_num; ++i)
	{
		if (BiCluster[i].index_column.size() < po->COL_WIDTH / 2)
			continue;

		size_t len2 = BiCluster[i].index_column.size();
		size_t l1 = BiCluster[i].in.size();
		size_t l2 = BiCluster[i].de.size();
		/*cout << "before row expand " << "i:" << i << " rowin: " << l1 << "rowde: " << l2 << " col: " << len2 << endl;
		for (int k = 0; k < l1; k++)
		{
			cout << BiCluster[i].in[k] << " ";
		}
		cout << endl;
		for (int k = 0; k < l2; k++)
			cout << BiCluster[i].de[k] << " ";
		cout << endl;
		for (int k = 0; k < len2; k++)
			cout << BiCluster[i].index_column[k] << " ";
		cout << endl;*/
		vector<vector<Array>> a(2, vector<Array>(len2));
		for (size_t j = 0; j < len2; ++j)
		{
			size_t col = BiCluster[i].index_column[j];

			a[0][j].con = col;
			a[0][j].val = B[BiCluster[i].x1][col].val;
			a[0][j].old_val = a[0][j].val;

			a[1][j].con = col;
			a[1][j].val = B[BiCluster[i].x2][col].val;
			a[1][j].old_val = a[1][j].val;
		}

		sort(a[0].begin(), a[0].end());
		sort(a[1].begin(), a[1].end());

		bool use_x2 = anchor_repeat_ok(a[1], BiCluster[i].x2);
		// cout << "use x2 check:"<<use_x2<<endl;
		size_t t = (size_t)-1;
		if (use_x2)
		{
			t = BiCluster_num;
			if (BiCluster_num >= BiCluster.size())
				BiCluster.resize(max(BiCluster.size() * 2 + 1, BiCluster_num + 1));
			BiCluster_num++;
			BiCluster[t].index_column = BiCluster[i].index_column;
			BiCluster[t].area = BiCluster[i].area;
			BiCluster[t].x1 = BiCluster[i].x2;
			BiCluster[t].x2 = BiCluster[i].x1;
			BiCluster[t].sig = BiCluster[i].sig;
			if(BiCluster[i].sig)
			{
				BiCluster[t].in = BiCluster[i].in;
				BiCluster[t].de = BiCluster[i].de;
			}
			else
			{
				BiCluster[t].in = BiCluster[i].de;
				BiCluster[t].de = BiCluster[i].in;
			}
		}

		vector<uint8_t> flag(n, 0);
		for (size_t j = 0; j < l1; ++j)
			flag[BiCluster[i].in[j]] = 1;
		for (size_t j = 0; j < l2; ++j)
			flag[BiCluster[i].de[j]] = 1;

		vector<future<PSS>> results;
		results.reserve(n);

		for (size_t j = 0; j < n; ++j)
		{
			if (flag[j])
				continue;
			results.emplace_back(
				pool.enqueue([j, &a, x1 = BiCluster[i].x1, x2 = BiCluster[i].x2, use_x2]()
							 { return CLUSTER_EXPAND(j, a, x1, x2, use_x2); }));
		}

		for (auto &&result : results)
		{
			PSS tmp = result.get();
			//cout << "try row: " << tmp.first.x1 << " len: " << tmp.first.len << " threshold: " << po->EXPAND_TOLERANCE * len2 << endl;
			if (tmp.first.len >= po->EXPAND_TOLERANCE * len2)
			{
				if (tmp.first.x2)
				{
					//cout << "add row: " << tmp.first.x1 << " to " << tmp.first.x2 << endl;
					BiCluster[i].in.emplace_back(tmp.first.x1);
				}
				else
				{
					//cout << "add row: " << tmp.first.x1 << " to " << tmp.first.x2 << endl;
					BiCluster[i].de.emplace_back(tmp.first.x1);
				}
			}
			//if (use_x2)
			{
				//cout << "try row: " << tmp.second.x1 << " len: " << tmp.second.len << " threshold: " << po->EXPAND_TOLERANCE * len2 << endl;
			}
			if (use_x2 && tmp.second.len >= po->EXPAND_TOLERANCE * len2)
			{
				if (tmp.second.x2)
				{
					//cout << "add row: " << tmp.second.x1 << " to " << tmp.second.x2 << endl;
					BiCluster[t].in.emplace_back(tmp.second.x1);
				}
				else
				{
					//cout << "add row: " << tmp.second.x1 << " to " << tmp.second.x2 << endl;
					BiCluster[t].de.emplace_back(tmp.second.x1);
				}
			}
		}

		size_t size_i = BiCluster[i].in.size() + BiCluster[i].de.size();
		/*cout << "i after row expand:" << i;
		cout << "seed:" << BiCluster[i].x1 << " rowin:" << BiCluster[i].in.size() << " rowde: " << BiCluster[i].de.size() << " ";
		cout << "col: " << BiCluster[i].index_column.size() << endl;*/
		/*for (int k = 0; k < BiCluster[i].in.size(); k++)
			cout << BiCluster[i].in[k] << " ";
		cout << endl;
		for (int k = 0; k < BiCluster[i].de.size(); k++)
			cout << BiCluster[i].de[k] << " ";
		cout << endl;
		for (int k = 0; k < BiCluster[i].index_column.size(); k++)
			cout << BiCluster[i].index_column[k] << " ";
		cout << endl;*/

		/*for (int j = 0; j < BiCluster[i].in.size(); j++)
		{
			sort(BiCluster[i].index_column.begin(), BiCluster[i].index_column.end(),
				 [&](size_t c1, size_t c2)
				 {
					 double v1 = B[BiCluster[i].in[j]][c1].val;
					 double v2 = B[BiCluster[i].in[j]][c2].val;
					 if (v1 != v2)
						 return v1 < v2;
					 return c1 < c2;
				 });
			cout << genes[BiCluster[i].in[j]] << ":\t";
			for (int k = 0; k < BiCluster[i].index_column.size(); k++)
			{
				cout << B[BiCluster[i].in[j]][BiCluster[i].index_column[k]].con << "-";
				cout << B[BiCluster[i].in[j]][BiCluster[i].index_column[k]].val << "\t";
			}
			cout << endl;
		}

		for (int j = 0; j < BiCluster[i].de.size(); j++)
		{
			sort(BiCluster[i].index_column.begin(), BiCluster[i].index_column.end(),
				 [&](size_t c1, size_t c2)
				 {
					 double v1 = B[BiCluster[i].de[j]][c1].val;
					 double v2 = B[BiCluster[i].de[j]][c2].val;
					 if (v1 != v2)
						 return v1 > v2;
					 return c1 < c2;
				 });
			cout << genes[BiCluster[i].de[j]] << ":\t";
			for (int k = 0; k < BiCluster[i].index_column.size(); k++)
			{
				cout << B[BiCluster[i].de[j]][BiCluster[i].index_column[k]].con << "-";
				cout << B[BiCluster[i].de[j]][BiCluster[i].index_column[k]].val << "\t";
			}
			cout << endl;
		}
		cout << endl;*/

		if (!use_x2)
			continue;

		size_t size_t_ = BiCluster[t].in.size() + BiCluster[t].de.size();

		//cout << "t after row expand:" << t;
		//cout << "seed:" << BiCluster[i].x2 << " rowin:" << BiCluster[t].in.size() << " rowde: " << BiCluster[t].de.size() << " ";
		//cout << "col: " << BiCluster[t].index_column.size() << endl;

		/*for (int j = 0; j < BiCluster[t].in.size(); j++)
		{
			sort(BiCluster[t].index_column.begin(), BiCluster[t].index_column.end(),
				 [&](size_t c1, size_t c2)
				 {
					 double v1 = B[BiCluster[t].in[j]][c1].val;
					 double v2 = B[BiCluster[t].in[j]][c2].val;
					 if (v1 != v2)
						 return v1 < v2;
					 return c1 < c2;
				 });
			cout << genes[BiCluster[t].in[j]] << ":\t";
			for (int k = 0; k < BiCluster[t].index_column.size(); k++)
			{
				cout << B[BiCluster[t].in[j]][BiCluster[t].index_column[k]].con << "-";
				cout << B[BiCluster[t].in[j]][BiCluster[t].index_column[k]].val << "\t";
			}
			cout << endl;
		}

		for (int j = 0; j < BiCluster[t].de.size(); j++)
		{
			sort(BiCluster[t].index_column.begin(), BiCluster[t].index_column.end(),
				 [&](size_t c1, size_t c2)
				 {
					 double v1 = B[BiCluster[t].de[j]][c1].val;
					 double v2 = B[BiCluster[t].de[j]][c2].val;
					 if (v1 != v2)
						 return v1 > v2;
					 return c1 < c2;
				 });
			cout << genes[BiCluster[t].de[j]] << ":\t";
			for (int k = 0; k < BiCluster[t].index_column.size(); k++)
			{
				cout << B[BiCluster[t].de[j]][BiCluster[t].index_column[k]].con << "-";
				cout << B[BiCluster[t].de[j]][BiCluster[t].index_column[k]].val << "\t";
			}
			cout << endl;
		}
		cout << endl;*/
		if (size_i < 2 && size_t_ >= 2)
		{
			BiCluster[i] = std::move(BiCluster[t]);

			BiCluster_num--;
			BiCluster[t].in.clear();
			BiCluster[t].de.clear();
			BiCluster[t].index_column.clear();
		}
		else if (size_t_ < 2)
		{
			BiCluster_num--;
			BiCluster[t].in.clear();
			BiCluster[t].de.clear();
			BiCluster[t].index_column.clear();
		}
	}
}
void cluster_expand()
{
	if (po->DICHOTOMY_TOLERANCE != 2 && po->DICHOTOMY_TOLERANCE != 4)
		column_expand();
	if (po->DICHOTOMY_TOLERANCE == 3 || po->DICHOTOMY_TOLERANCE == 4)
	{
		uglyTime("%d bicluster generated after expanding", BiCluster_num);
		return;
	}
	ThreadPool pool(max<size_t>(po->THREADS_NUM, 1));
	row_expand(pool);
	// column_expand();
	uglyTime("%d bicluster generated after expanding", BiCluster_num);
}
