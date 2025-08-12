#include "structure.h"


void print_path(size_t i, size_t j, size_t k, size_t sig, const LCS_NODE &T, Ans &ans)
{  
	if (k <= 0) return;
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
	if (min(fabs(a),fabs(b)) < eps) return fabs(a - b) >= ratio*po->LCS_TOLERANCE;
	return (fabs(a - b) / min(fabs(b),fabs(a))) >= ratio*po->LCS_TOLERANCE;
}
double check1(size_t x1, size_t x2, size_t sig, double ts,const LCS_NODE &T)
{
	auto it11 = T.sig1[sig].find(x2);
	auto it22 = T.sig2[sig].find(x1);
	if(it11 == T.sig1[sig].end()||it22 == T.sig2[sig].end()) return -1;
	auto it1 = T.sig1[sig].find(x1);
	auto it2 = T.sig2[sig].find(x2);
	size_t id1 = it1->second;
	size_t id11 = it11->second;
	size_t id2 = it2->second;
	size_t id22 = it22->second;
	if(!(id11>id1&&id22>id2)) return -1;
	if (check_val(T.a2[sig][id22].val, T.a2[sig][id2].val, ts) || check_val(T.a2[sig][id22].val, T.a2[sig][id2].old_val, ts)) return -1;
	if (check_val(T.a2[sig][id2].val, T.a2[sig][id22].val, ts) || check_val(T.a2[sig][id2].val, T.a2[sig][id22].old_val, ts)) return -1;
	ts = fabs(ts);
	if (ts < eps) ts = 1;
	double tmp = min(fabs(T.a2[sig][id22].val), fabs(T.a2[sig][id2].val));
	if (tmp < eps) tmp = 1;
	return fabs(T.a2[sig][id22].val - T.a2[sig][id2].val) / tmp / ts;
}
double check2(size_t x1, size_t x2, size_t sig, double ts, const LCS_NODE &T)
{
	auto it11 = T.sig1[sig].find(x2);
	
	auto it22 = T.sig2[sig].find(x1);
	if(it11 == T.sig1[sig].end()||it22 == T.sig2[sig].end()) return -1;
	auto it1 = T.sig1[sig].find(x1);
	auto it2 = T.sig2[sig].find(x2);
	size_t id1 = it1->second;
	size_t id11 = it11->second;
	size_t id2 = it2->second;
	size_t id22 = it22->second;
	
	if(!(id11>id1&&id22>id2)) return -1;
	if (check_val(T.a1[sig][id11].val, T.a1[sig][id1].val, ts) || check_val(T.a1[sig][id11].val, T.a1[sig][id1].old_val, ts)) return -1;//第一个check限制交换两个数差值，第二个限制该数和原本的差值
	if (check_val(T.a1[sig][id1].val, T.a1[sig][id11].val, ts) || check_val(T.a1[sig][id1].val, T.a1[sig][id11].old_val, ts)) return -1;
	ts = fabs(ts);
	if (ts < eps) ts = 1;
	double tmp = min(fabs(T.a1[sig][id11].val), fabs(T.a1[sig][id1].val));
	if (tmp < eps) tmp = 1;
	return (fabs(T.a1[sig][id11].val - T.a1[sig][id1].val) / tmp) / ts;
}

void solve(size_t id1, size_t id2, size_t sig, size_t i, size_t j, double ts1, double ts2, LCS_NODE &T)
{
	if (id1 == id2)
	{
		T.dp[sig][i][j] = T.dp[sig][i - 1][j - 1] + 1;
		T.path[sig][i][j] = 1;
		return ;
	}
	double tmp1 = check1(id1, id2, sig, ts2, T), tmp2 = check2(id1, id2, sig, ts1, T);

	if (tmp1 == -1 && tmp2 == -1)
	{
		if (T.dp[sig][i - 1][j] > T.dp[sig][i][j - 1])
		{
			T.dp[sig][i][j] = T.dp[sig][i - 1][j];
			T.path[sig][i][j] = 2;
		}
		else
		{
			T.dp[sig][i][j] = T.dp[sig][i][j - 1];
			T.path[sig][i][j] = 0;
		}
		return ;
	}
	else if (tmp1 != -1 && (tmp2 == -1 || tmp2 >= tmp1))
	{
		T.dp[sig][i][j] = T.dp[sig][i - 1][j - 1] + 1;
		T.path[sig][i][j] = 1;
		double temp;

		T.a2[sig][T.sig2[sig][id2]].con = id1;
		T.a2[sig][T.sig2[sig][id1]].con = id2;

		temp = T.a2[sig][T.sig2[sig][id1]].old_val;
		T.a2[sig][T.sig2[sig][id1]].old_val = T.a2[sig][T.sig2[sig][id2]].old_val;
		T.a2[sig][T.sig2[sig][id2]].old_val = temp;
		swap(T.sig2[sig][id2],T.sig2[sig][id1]);
		return ;
	}
	else
	{
		T.dp[sig][i][j] = T.dp[sig][i - 1][j - 1] + 1;
		T.path[sig][i][j] = 1;
		double temp;

		T.a1[sig][T.sig1[sig][id2]].con = id1;
		T.a1[sig][T.sig1[sig][id1]].con = id2;

		temp = T.a1[sig][T.sig1[sig][id1]].old_val;
		T.a1[sig][T.sig1[sig][id1]].old_val = T.a1[sig][T.sig1[sig][id2]].old_val;
		T.a1[sig][T.sig1[sig][id2]].old_val = temp;

		swap(T.sig1[sig][id2],T.sig1[sig][id1]);
		return ;
	}
}

size_t repeat_check(const LCS_NODE& T,const size_t& sig,const Ans& ans,const double& ts1,const double& ts2)
{
	size_t max_cnt1=0,max_cnt2=0;
	size_t cnt1=0,cnt2=0;
	
	double pre1=-2,pre2=-2;
	
	for(size_t i=1;i<=ans.len;i++)
	{
		if(fabs(T.a1[sig][T.sig1[sig].at(ans.p[i])].val-pre1)/min(fabs(T.a1[sig][T.sig1[sig].at(ans.p[i])].val),fabs(pre1))>=ts1*po->LCS_TOLERANCE)
		{
			max_cnt1=max(cnt1,max_cnt1);
			cnt1=1;
			pre1=T.a1[sig][T.sig1[sig].at(ans.p[i])].val;
		}
		else cnt1++;
		if(fabs(T.a2[sig][T.sig2[sig].at(ans.p[i])].val-pre2)/min(fabs(T.a2[sig][T.sig2[sig].at(ans.p[i])].val),fabs(pre2))>=ts2*po->LCS_TOLERANCE)
		{
			max_cnt2=max(cnt2,max_cnt2);
			cnt2=1;
			pre2=T.a2[sig][T.sig2[sig].at(ans.p[i])].val;
		}
		else cnt2++;
	}
	max_cnt1=max(cnt1,max_cnt1);
	max_cnt2=max(cnt2,max_cnt2);
	if(max_cnt1*1.0/ans.len>po->DISCRETIZATION&&max_cnt2*1.0/ans.len>po->DISCRETIZATION)
	{
		return 0;
	}
	if(max_cnt1*1.0/ans.len<=max_cnt2*1.0/ans.len||max_cnt2*1.0/ans.len>po->DISCRETIZATION)
	{
		return 1;
	}
	else
	{
		return 2;
	}
}
Ans LCS(const vector<Array> &A1, const vector<Array> &A2, size_t n1, size_t n2, size_t is_print,size_t is_seed_gen,size_t st1,size_t st2,size_t is_repeat_check)
{
	Ans ans;
	size_t l1 = 0, l2 = 0;
	
	if(is_seed_gen)
	{
		l1=Start_pos[st1];
		l2=Start_pos[st2];
	}
	else
	{
		while ((l1 < n1)&&(fabs(A1[l1].val + 2) < eps)) l1++;
		while ((l2 < n2)&&(fabs(A2[l2].val + 2) < eps)) l2++;
	}

	if((n1-l1<po->COL_WIDTH)||(n2-l2<po->COL_WIDTH))
	{
		ans.len=0;
		ans.sig=1;
		return ans;
	}
	
	LCS_NODE T;
	T.path = vector<vector<vector<size_t> > >(2, vector<vector<size_t> >(n1-l1+5, vector<size_t>(n2-l2+5)));
	T.dp = vector<vector<vector<size_t> > >(2, vector<vector<size_t> >(n1-l1+5, vector<size_t>(n2-l2+5, 0)));
	
	T.sig1=vector<unordered_map<size_t,size_t> >(2,unordered_map<size_t,size_t>(n1-l1));
	T.sig2=vector<unordered_map<size_t,size_t> >(2,unordered_map<size_t,size_t>(n2-l2));
	T.a1 = vector<vector<Array> >(2, vector<Array>(n1-l1));
	
	T.a2 = vector<vector<Array> >(2, vector<Array>(n2 - l2));
	double sum1 = 0, sum2 = 0, var1 = 0, var2 = 0;
	size_t ll1 = floor((n1 - l1) / 4.0), rr1 = ceil((n1 - l1)*3 / 4.0);
	for (size_t i = 0; i < n1 - l1; i++)
	{
		T.a1[0][i].con = A1[i + l1].con;
		T.a1[1][i].val = A1[i + l1].val;
		T.a1[1][i].con = A1[i + l1].con;
		T.a1[0][i].val = A1[i + l1].val;
		T.a1[1][i].old_val = A1[i + l1].val;
		T.a1[0][i].old_val = A1[i + l1].val;
		T.sig1[0][A1[i + l1].con] = i;
		T.sig1[1][A1[i + l1].con] = i;
		if (i >= ll1 && i < rr1) 
		{
			sum1 += A1[i + l1].val;
		}
	}
	
	size_t ll2 = floor((n2 - l2) / 4.0), rr2 = ceil((n2 - l2)*3 / 4.0);
	
	for (size_t i = 0; i < n2 - l2; i++)
	{
		T.a2[0][i].con = A2[i + l2].con;
		T.a2[0][i].old_val = A2[i + l2].val;
		T.a2[0][i].val = A2[i + l2].val;
		T.a2[1][i].con = A2[i + l2].con;
		if (A2[i + l2].val == -2)
		{
			T.a2[1][i].val = -2;
			T.a2[1][i].old_val = -2;
		}
		else
		{
			T.a2[1][i].val = -A2[i + l2].val;
			T.a2[1][i].old_val = -A2[i + l2].old_val;
		}
		T.sig2[0][A2[i + l2].con] = i;
		if (i >= ll2 && i < rr2) 
		{
			sum2 += A2[i + l2].val;
		}
	}
	
	sort(T.a2[1].begin(), T.a2[1].end());
	for (size_t i = 0; i < n2 - l2; i++)
	{
		T.sig2[1][T.a2[1][i].con] = i;
	}
	

	sum1 = sum1 / (rr1 - ll1);
	for (size_t i = ll1; i < rr1; i++)
	{
		var1 += (A1[i + l1].val - sum1)*(A1[i + l1].val - sum1);
	}
	var1 = sqrt(var1 / (rr1 - ll1));
	
	sum2 = sum2 / (rr2 - ll2);
	
	for (size_t i = ll2; i < rr2; i++)
	{
		var2 += (A2[i + l2].val - sum2)*(A2[i + l2].val - sum2);
	}
	var2 = sqrt(var2 / (rr2 - ll2));

	double ts1, ts2;
	if(fabs(sum1)<eps) ts1=0;
	else ts1=fabs(var1/sum1);
	if(fabs(sum2)<eps) ts2=0;
	else ts2=fabs(var2/sum2);
	
	for (size_t i = 1; i <= n1 - l1; i++)
	{
		
		for (size_t j = 1; j <= n2 - l2; j++)
		{
			solve(T.a1[0][i - 1].con, T.a2[0][j - 1].con, 0, i, j, ts1, ts2, T);
			solve(T.a1[1][i - 1].con, T.a2[1][j - 1].con, 1, i, j, ts1, ts2, T);
		}
	}

	if (T.dp[0][n1 - l1][n2 - l2] >= T.dp[1][n1 - l1][n2 - l2])
	{
		ans.len = T.dp[0][n1 - l1][n2 - l2];
		ans.sig = 1;
		ans.p = vector<size_t>(ans.len + 1);
		if (is_print)
		{
			T.path_cnt=0;
			print_path(n1 - l1, n2 - l2, T.dp[0][n1 - l1][n2 - l2], 0, T, ans);
			if(is_repeat_check) ans.repeat_check=repeat_check(T,0,ans,ts1,ts2);
		}
	}
	else
	{
		ans.len = T. dp[1][n1 - l1][n2 - l2];
		ans.sig = 0;
		ans.p = vector<size_t>(ans.len + 1);
		if (is_print)
		{
			T.path_cnt=0;
			print_path(n1 - l1, n2 - l2, T.dp[1][n1 - l1][n2 - l2], 1, T, ans);
			if(is_repeat_check) ans.repeat_check=repeat_check(T,1,ans,ts1,ts2);
		}
	}
	return ans;
}
