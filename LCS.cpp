#include "structure.h"
using namespace std;
void print_path(size_t i, size_t j, size_t k, size_t sig, const LCS_NODE &T, Ans &ans)
{
	if (k <= 0)
		return;
	if (T.path[sig][i][j] == 1)
	{
		ans.p[k] = T.a1[sig][i - 1].con;
		print_path(i - 1, j - 1, k - 1, sig, T, ans);
	}
	else if (T.path[sig][i][j] == 0)
	{
		print_path(i, j - 1, k, sig, T, ans);
	}
	else
	{
		print_path(i - 1, j, k, sig, T, ans);
	}
}

inline bool check_val(double a, double b, double ratio)
{
	if (min(fabs(a), fabs(b)) < eps)
		return fabs(a - b) >= ratio * po->LCS_TOLERANCE;
	return (fabs(a - b) / min(fabs(b), fabs(a))) >= ratio * po->LCS_TOLERANCE;
}
double check1(size_t x1, size_t x2, size_t sig, double ts, const LCS_NODE &T)
{
	if (isnan(ts))
		return -1;
	const size_t missing_pos = numeric_limits<size_t>::max();
	if (x1 >= T.sig1[sig].size() || x2 >= T.sig1[sig].size() ||
		x1 >= T.sig2[sig].size() || x2 >= T.sig2[sig].size())
		return -1;
	size_t id1 = T.sig1[sig][x1];
	size_t id11 = T.sig1[sig][x2];
	size_t id2 = T.sig2[sig][x2];
	size_t id22 = T.sig2[sig][x1];
	if (id1 == missing_pos || id11 == missing_pos ||
		id2 == missing_pos || id22 == missing_pos)
		return -1;
	if (!(id11 > id1 && id22 > id2))
		return -1;
	if (T.constant_check[1] && fabs(T.a2[sig][id22].val - T.a2[sig][id2].val) < eps)
		return -1;
	if (T.repeat_check[1] && check_val(T.a2[sig][id22].val, T.a2[sig][id2].val, ts * po->DISCRETIZATION))
		return -1;
	if (check_val(T.a2[sig][id22].val, T.a2[sig][id2].val, ts) || check_val(T.a2[sig][id22].val, T.a2[sig][id2].old_val, ts))
		return -1;
	if (check_val(T.a2[sig][id2].val, T.a2[sig][id22].val, ts) || check_val(T.a2[sig][id2].val, T.a2[sig][id22].old_val, ts))
		return -1;
	ts = fabs(ts);
	if (ts < eps)
		ts = 1;
	double tmp = min(fabs(T.a2[sig][id22].val), fabs(T.a2[sig][id2].val));
	if (tmp < eps)
		tmp = 1;
	return fabs(T.a2[sig][id22].val - T.a2[sig][id2].val) / tmp / ts;
}
double check2(size_t x1, size_t x2, size_t sig, double ts, const LCS_NODE &T)
{
	if (isnan(ts))
		return -1;
	const size_t missing_pos = numeric_limits<size_t>::max();
	if (x1 >= T.sig1[sig].size() || x2 >= T.sig1[sig].size() ||
		x1 >= T.sig2[sig].size() || x2 >= T.sig2[sig].size())
		return -1;
	size_t id1 = T.sig1[sig][x1];
	size_t id11 = T.sig1[sig][x2];
	size_t id2 = T.sig2[sig][x2];
	size_t id22 = T.sig2[sig][x1];
	if (id1 == missing_pos || id11 == missing_pos ||
		id2 == missing_pos || id22 == missing_pos)
		return -1;

	if (!(id11 > id1 && id22 > id2))
		return -1;
	if (T.constant_check[0] && fabs(T.a1[sig][id11].val - T.a1[sig][id1].val) < eps)
		return -1;
	if (T.repeat_check[0] && check_val(T.a1[sig][id11].val, T.a1[sig][id1].val, ts * po->DISCRETIZATION))
		return -1;
	if (check_val(T.a1[sig][id11].val, T.a1[sig][id1].val, ts) || check_val(T.a1[sig][id11].val, T.a1[sig][id1].old_val, ts))
		return -1;
	if (check_val(T.a1[sig][id1].val, T.a1[sig][id11].val, ts) || check_val(T.a1[sig][id1].val, T.a1[sig][id11].old_val, ts))
		return -1;
	ts = fabs(ts);
	if (ts < eps)
		ts = 1;
	double tmp = min(fabs(T.a1[sig][id11].val), fabs(T.a1[sig][id1].val));
	if (tmp < eps)
		tmp = 1;
	return (fabs(T.a1[sig][id11].val - T.a1[sig][id1].val) / tmp) / ts;
}

static double check_seed1(size_t x1, size_t x2, size_t sig, double ts, const SeedLCSWorkspace &T)
{
	if (isnan(ts))
		return -1;
	const size_t missing_pos = numeric_limits<size_t>::max();
	if (x1 >= T.sig1[sig].size() || x2 >= T.sig1[sig].size() ||
		x1 >= T.sig2[sig].size() || x2 >= T.sig2[sig].size())
		return -1;
	size_t id1 = T.sig1[sig][x1];
	size_t id11 = T.sig1[sig][x2];
	size_t id2 = T.sig2[sig][x2];
	size_t id22 = T.sig2[sig][x1];
	if (id1 == missing_pos || id11 == missing_pos ||
		id2 == missing_pos || id22 == missing_pos)
		return -1;
	if (!(id11 > id1 && id22 > id2))
		return -1;
	if (T.constant_check[1] && fabs(T.a2[sig][id22].val - T.a2[sig][id2].val) < eps)
		return -1;
	if (T.repeat_check[1] && check_val(T.a2[sig][id22].val, T.a2[sig][id2].val, ts * po->DISCRETIZATION))
		return -1;
	if (check_val(T.a2[sig][id22].val, T.a2[sig][id2].val, ts) || check_val(T.a2[sig][id22].val, T.a2[sig][id2].old_val, ts))
		return -1;
	if (check_val(T.a2[sig][id2].val, T.a2[sig][id22].val, ts) || check_val(T.a2[sig][id2].val, T.a2[sig][id22].old_val, ts))
		return -1;
	ts = fabs(ts);
	if (ts < eps)
		ts = 1;
	double tmp = min(fabs(T.a2[sig][id22].val), fabs(T.a2[sig][id2].val));
	if (tmp < eps)
		tmp = 1;
	return fabs(T.a2[sig][id22].val - T.a2[sig][id2].val) / tmp / ts;
}

static double check_seed2(size_t x1, size_t x2, size_t sig, double ts, const SeedLCSWorkspace &T)
{
	if (isnan(ts))
		return -1;
	const size_t missing_pos = numeric_limits<size_t>::max();
	if (x1 >= T.sig1[sig].size() || x2 >= T.sig1[sig].size() ||
		x1 >= T.sig2[sig].size() || x2 >= T.sig2[sig].size())
		return -1;
	size_t id1 = T.sig1[sig][x1];
	size_t id11 = T.sig1[sig][x2];
	size_t id2 = T.sig2[sig][x2];
	size_t id22 = T.sig2[sig][x1];
	if (id1 == missing_pos || id11 == missing_pos ||
		id2 == missing_pos || id22 == missing_pos)
		return -1;

	if (!(id11 > id1 && id22 > id2))
		return -1;
	if (T.constant_check[0] && fabs(T.a1[sig][id11].val - T.a1[sig][id1].val) < eps)
		return -1;
	if (T.repeat_check[0] && check_val(T.a1[sig][id11].val, T.a1[sig][id1].val, ts * po->DISCRETIZATION))
		return -1;
	if (check_val(T.a1[sig][id11].val, T.a1[sig][id1].val, ts) || check_val(T.a1[sig][id11].val, T.a1[sig][id1].old_val, ts))
		return -1;
	if (check_val(T.a1[sig][id1].val, T.a1[sig][id11].val, ts) || check_val(T.a1[sig][id1].val, T.a1[sig][id11].old_val, ts))
		return -1;
	ts = fabs(ts);
	if (ts < eps)
		ts = 1;
	double tmp = min(fabs(T.a1[sig][id11].val), fabs(T.a1[sig][id1].val));
	if (tmp < eps)
		tmp = 1;
	return (fabs(T.a1[sig][id11].val - T.a1[sig][id1].val) / tmp) / ts;
}

static size_t solve_seed(size_t id1, size_t id2, size_t sig,
						 size_t j,
						 double ts1, double ts2,
						 SeedLCSWorkspace &T)
{
	const size_t diag = T.prev[sig][j - 1];
	const size_t up = T.prev[sig][j];
	const size_t left = T.cur[sig][j - 1];

	if (id1 == id2)
		return diag + 1;

	double tmp1 = check_seed1(id1, id2, sig, ts2, T);
	double tmp2 = check_seed2(id1, id2, sig, ts1, T);

	if (tmp1 == -1 && tmp2 == -1)
		return (up > left) ? up : left;

	if (tmp1 != -1 && (tmp2 == -1 || tmp2 >= tmp1))
	{
		T.a2[sig][T.sig2[sig][id2]].con = id1;
		T.a2[sig][T.sig2[sig][id1]].con = id2;

		double temp = T.a2[sig][T.sig2[sig][id1]].old_val;
		T.a2[sig][T.sig2[sig][id1]].old_val = T.a2[sig][T.sig2[sig][id2]].old_val;
		T.a2[sig][T.sig2[sig][id2]].old_val = temp;
		swap(T.sig2[sig][id2], T.sig2[sig][id1]);
		return diag + 1;
	}

	T.a1[sig][T.sig1[sig][id2]].con = id1;
	T.a1[sig][T.sig1[sig][id1]].con = id2;

	double temp = T.a1[sig][T.sig1[sig][id1]].old_val;
	T.a1[sig][T.sig1[sig][id1]].old_val = T.a1[sig][T.sig1[sig][id2]].old_val;
	T.a1[sig][T.sig1[sig][id2]].old_val = temp;

	swap(T.sig1[sig][id2], T.sig1[sig][id1]);
	return diag + 1;
}

static void prepare_seed_workspace(SeedLCSWorkspace &T, size_t len1, size_t len2)
{
	const size_t missing_pos = numeric_limits<size_t>::max();

	T.sig1.resize(2);
	T.sig2.resize(2);
	T.a1.resize(2);
	T.a2.resize(2);

	for (size_t sig = 0; sig < 2; ++sig)
	{
		T.sig1[sig].assign(m, missing_pos);
		T.sig2[sig].assign(m, missing_pos);
		T.a1[sig].resize(len1);
		T.a2[sig].resize(len2);
		T.prev[sig].assign(len2 + 1, 0);
		T.cur[sig].assign(len2 + 1, 0);
	}
}

void solve(size_t id1, size_t id2, size_t sig, size_t i, size_t j, double ts1, double ts2, LCS_NODE &T)
{
	if (id1 == id2)
	{
		T.dp[sig][i][j] = T.dp[sig][i - 1][j - 1] + 1;
		if (T.store_path)
			T.path[sig][i][j] = 1;
		return;
	}
	double tmp1 = check1(id1, id2, sig, ts2, T), tmp2 = check2(id1, id2, sig, ts1, T);
	// if(tmp1!=-1||tmp2!=-1)
	//	cout<<"check1:"<<tmp1<<" check2: "<<tmp2<<endl;
	if (tmp1 == -1 && tmp2 == -1)
	{
		if (T.dp[sig][i - 1][j] > T.dp[sig][i][j - 1])
		{
			T.dp[sig][i][j] = T.dp[sig][i - 1][j];
			if (T.store_path)
				T.path[sig][i][j] = 2;
		}
		else
		{
			T.dp[sig][i][j] = T.dp[sig][i][j - 1];
			if (T.store_path)
				T.path[sig][i][j] = 0;
		}
		return;
	}
	else if (tmp1 != -1 && (tmp2 == -1 || tmp2 >= tmp1))
	{
		// cout<<"case 1"<<endl;
		T.dp[sig][i][j] = T.dp[sig][i - 1][j - 1] + 1;
		if (T.store_path)
			T.path[sig][i][j] = 1;
		double temp;

		T.a2[sig][T.sig2[sig][id2]].con = id1;
		T.a2[sig][T.sig2[sig][id1]].con = id2;

		temp = T.a2[sig][T.sig2[sig][id1]].old_val;
		T.a2[sig][T.sig2[sig][id1]].old_val = T.a2[sig][T.sig2[sig][id2]].old_val;
		T.a2[sig][T.sig2[sig][id2]].old_val = temp;
		swap(T.sig2[sig][id2], T.sig2[sig][id1]);
		return;
	}
	else
	{
		// cout<<"case 2"<<endl;
		T.dp[sig][i][j] = T.dp[sig][i - 1][j - 1] + 1;
		if (T.store_path)
			T.path[sig][i][j] = 1;
		double temp;

		T.a1[sig][T.sig1[sig][id2]].con = id1;
		T.a1[sig][T.sig1[sig][id1]].con = id2;

		temp = T.a1[sig][T.sig1[sig][id1]].old_val;
		T.a1[sig][T.sig1[sig][id1]].old_val = T.a1[sig][T.sig1[sig][id2]].old_val;
		T.a1[sig][T.sig1[sig][id2]].old_val = temp;

		swap(T.sig1[sig][id2], T.sig1[sig][id1]);
		return;
	}
}

SequenceRepeatStats sequence_repeat_stats(const vector<Array> &a, size_t len, size_t row)
{
	SequenceRepeatStats stats;
	if (len == 0)
		return stats;

	size_t max_cnt = 1;
	size_t cnt = 1;
	double pre = a[0].val;
	size_t exact_max_cnt = 1;
	size_t exact_cnt = 1;
	double exact_pre = a[0].val;

	for (size_t i = 1; i < len; ++i)
	{
		const double cur = a[i].val;
		const double ts = NoiseMatrix[row][a[i].con];

		//if (fabs(cur - pre) < eps)
		if (fabs(cur - pre) < eps || (!isnan(ts) && !check_val(cur, pre, ts)))
			++cnt; 
		else
		{
			max_cnt = max(max_cnt, cnt);
			cnt = 1;
			pre = cur;
		}

		if (fabs(cur - exact_pre) < eps)
			++exact_cnt;
		else
		{
			exact_max_cnt = max(exact_max_cnt, exact_cnt);
			exact_cnt = 1;
			exact_pre = cur;
		}
	}

	max_cnt = max(max_cnt, cnt);
	exact_max_cnt = max(exact_max_cnt, exact_cnt);
	stats.repeat = max_cnt * 1.0 / len > po->DISCRETIZATION;
	stats.constant = exact_max_cnt * 1.0 / len > po->DISCRETIZATION || exact_max_cnt >= po->COL_WIDTH;
	return stats;
}

Ans LCS(const vector<Array> &A1, const vector<Array> &A2,
		size_t n1, size_t n2, size_t r1, size_t r2,
		size_t is_print)
{
	Ans ans;
	ans.len = 0;
	ans.sig = 1;

	size_t len1 = 0, len2 = 0;
	while (len1 < n1 && !isnan(A1[len1].val))
		++len1;
	while (len2 < n2 && !isnan(A2[len2].val))
		++len2;

	if (len1 < po->COL_WIDTH || len2 < po->COL_WIDTH)
	{
		ans.len = 0;
		ans.sig = 1;
		return ans;
	}

	LCS_NODE T;
	T.store_path = (is_print != 0);
	if (T.store_path)
	{
		T.path = vector<vector<vector<size_t>>>(
			2, vector<vector<size_t>>(len1 + 5, vector<size_t>(len2 + 5, 0)));
	}

	T.dp = vector<vector<vector<size_t>>>(
		2, vector<vector<size_t>>(len1 + 5, vector<size_t>(len2 + 5, 0)));

	const size_t missing_pos = numeric_limits<size_t>::max();
	T.sig1 = vector<vector<size_t>>(2, vector<size_t>(m, missing_pos));
	T.sig2 = vector<vector<size_t>>(2, vector<size_t>(m, missing_pos));

	T.a1 = vector<vector<Array>>(2, vector<Array>(len1));
	T.a2 = vector<vector<Array>>(2, vector<Array>(len2));

	for (size_t i = 0; i < len1; i++)
	{
		T.a1[0][i].con = A1[i].con;
		T.a1[0][i].val = A1[i].val;
		T.a1[0][i].old_val = A1[i].old_val;

		T.a1[1][i].con = A1[i].con;
		T.a1[1][i].val = A1[i].val;
		T.a1[1][i].old_val = A1[i].old_val;

		T.sig1[0][A1[i].con] = i;
		T.sig1[1][A1[i].con] = i;
	}

	for (size_t i = 0; i < len2; i++)
	{

		T.a2[0][i].con = A2[i].con;
		T.a2[0][i].val = A2[i].val;
		T.a2[0][i].old_val = A2[i].old_val;

		T.a2[1][i].con = A2[i].con;
		T.a2[1][i].val = -A2[i].val;
		T.a2[1][i].old_val = -A2[i].old_val;

		T.sig2[0][A2[i].con] = i;
	}

	sort(T.a2[1].begin(), T.a2[1].end());
	for (size_t i = 0; i < len2; i++)
	{
		T.sig2[1][T.a2[1][i].con] = i;
	}

	// ---------- DP ----------
	
	SequenceRepeatStats stats1 = sequence_repeat_stats(T.a1[0], len1, r1);
	SequenceRepeatStats stats2 = sequence_repeat_stats(T.a2[0], len2, r2);
	T.repeat_check[0] = stats1.repeat;
	T.repeat_check[1] = stats2.repeat;
	T.constant_check[0] = stats1.constant;
	T.constant_check[1] = stats2.constant;
	
	for (size_t i = 1; i <= len1; i++)
	{
		for (size_t j = 1; j <= len2; j++)
		{
			double ts1 = NoiseMatrix[r1][T.a1[0][i - 1].con], ts2 = NoiseMatrix[r2][T.a2[0][j - 1].con];
			solve(T.a1[0][i - 1].con, T.a2[0][j - 1].con, 0, i, j, ts1, ts2, T);
			ts1 = NoiseMatrix[r1][T.a1[1][i - 1].con], ts2 = NoiseMatrix[r2][T.a2[1][j - 1].con];
			solve(T.a1[1][i - 1].con, T.a2[1][j - 1].con, 1, i, j, ts1, ts2, T);
		}
	}

	if (T.dp[0][len1][len2] >= T.dp[1][len1][len2])
	{
		ans.len = T.dp[0][len1][len2];
		ans.sig = 1;
		ans.p = vector<size_t>(ans.len + 1);

		if (is_print)
		{
			print_path(len1, len2, T.dp[0][len1][len2], 0, T, ans);
		}
	}
	else
	{
		ans.len = T.dp[1][len1][len2];
		ans.sig = 0;
		ans.p = vector<size_t>(ans.len + 1);

		if (is_print)
		{
			print_path(len1, len2, T.dp[1][len1][len2], 1, T, ans);
		}
	}

	return ans;
}

Ans LCS_seed(const vector<Array> &A1, const vector<Array> &A2,
			 size_t n1, size_t n2, size_t r1, size_t r2,
			 SeedLCSWorkspace &T,
			 const vector<uint8_t> &row_repeat,
			 const vector<uint8_t> &row_constant)
{
	Ans ans;
	ans.len = 0;
	ans.sig = 1;

	size_t len1 = 0, len2 = 0;
	while (len1 < n1 && !isnan(A1[len1].val))
		++len1;
	while (len2 < n2 && !isnan(A2[len2].val))
		++len2;

	if (len1 < po->COL_WIDTH || len2 < po->COL_WIDTH)
		return ans;

	prepare_seed_workspace(T, len1, len2);

	for (size_t i = 0; i < len1; i++)
	{
		T.a1[0][i].con = A1[i].con;
		T.a1[0][i].val = A1[i].val;
		T.a1[0][i].old_val = A1[i].old_val;

		T.a1[1][i].con = A1[i].con;
		T.a1[1][i].val = A1[i].val;
		T.a1[1][i].old_val = A1[i].old_val;

		T.sig1[0][A1[i].con] = i;
		T.sig1[1][A1[i].con] = i;
	}

	for (size_t i = 0; i < len2; i++)
	{
		T.a2[0][i].con = A2[i].con;
		T.a2[0][i].val = A2[i].val;
		T.a2[0][i].old_val = A2[i].old_val;

		T.a2[1][i].con = A2[i].con;
		T.a2[1][i].val = -A2[i].val;
		T.a2[1][i].old_val = -A2[i].old_val;

		T.sig2[0][A2[i].con] = i;
	}

	sort(T.a2[1].begin(), T.a2[1].end());
	for (size_t i = 0; i < len2; i++)
		T.sig2[1][T.a2[1][i].con] = i;

	if (row_repeat.size() > max(r1, r2))
	{
		T.repeat_check[0] = (row_repeat[r1] != 0);
		T.repeat_check[1] = (row_repeat[r2] != 0);
	}
	else
	{
		SequenceRepeatStats stats1 = sequence_repeat_stats(T.a1[0], len1, r1);
		SequenceRepeatStats stats2 = sequence_repeat_stats(T.a2[0], len2, r2);
		T.repeat_check[0] = stats1.repeat;
		T.repeat_check[1] = stats2.repeat;
		T.constant_check[0] = stats1.constant;
		T.constant_check[1] = stats2.constant;
	}
	if (row_constant.size() > max(r1, r2))
	{
		T.constant_check[0] = (row_constant[r1] != 0);
		T.constant_check[1] = (row_constant[r2] != 0);
	}

	for (size_t i = 1; i <= len1; i++)
	{
		fill(T.cur[0].begin(), T.cur[0].end(), 0);
		fill(T.cur[1].begin(), T.cur[1].end(), 0);

		for (size_t j = 1; j <= len2; j++)
		{
			double ts1 = NoiseMatrix[r1][T.a1[0][i - 1].con];
			double ts2 = NoiseMatrix[r2][T.a2[0][j - 1].con];
			T.cur[0][j] = solve_seed(T.a1[0][i - 1].con, T.a2[0][j - 1].con, 0, j, ts1, ts2, T);

			ts1 = NoiseMatrix[r1][T.a1[1][i - 1].con];
			ts2 = NoiseMatrix[r2][T.a2[1][j - 1].con];
			T.cur[1][j] = solve_seed(T.a1[1][i - 1].con, T.a2[1][j - 1].con, 1, j, ts1, ts2, T);
		}

		swap(T.prev[0], T.cur[0]);
		swap(T.prev[1], T.cur[1]);
	}

	if (T.prev[0][len2] >= T.prev[1][len2])
	{
		ans.len = T.prev[0][len2];
		ans.sig = 1;
	}
	else
	{
		ans.len = T.prev[1][len2];
		ans.sig = 0;
	}

	return ans;
}
