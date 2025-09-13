#include "structure.h"
#include <thread>
#include <mutex>

priority_queue<SEED, vector<SEED>, cmp2> Seed;
SEED SEED_GEN(const size_t &k, const size_t &j)
{
	Ans ans;
	ans = LCS(A[k], A[j], m, m, 1, 1, k, j,1); 

	SEED seed;
	seed.len = ans.len;
	seed.x1 = k, seed.x2 = j;
	seed.flag=false;
	if(ans.len<po->COL_WIDTH) return seed;
	
	if(ans.repeat_check==0)
	{
		seed.len=0;
		return seed;
	}
	
	if(ans.repeat_check==2)
	{
		seed.flag=true;
		seed.x1 = j, seed.x2 = k;
	}
	return seed;
}

void seed_generation()
{
	size_t len = (size_t)ceil(n / 4.0);
	len = (len > 2 ? len : n);
	int step, k = -1, j = 0, end;
	step = len;
	ThreadPool pool(po->THREADS_NUM);
	vector<future<SEED> > results;
	priority_queue<SEED, vector<SEED>, cmp1> minseed;
	while (++k < n)
	{
		j = k;
		if (n - (int)(k / step) * step < step)
			end = n;
		else end = (int)(k / step + 1) * step;
		
		while (++j < end)
		{
			results.emplace_back(pool.enqueue(SEED_GEN, k, j));
		}
	}

	 
	for (auto &&result : results)
	{
		SEED tmp=result.get();
		if (tmp.len > po->COL_WIDTH) minseed.emplace(tmp);
		while (minseed.size() > po->SEED_NUM)
			minseed.pop(); 
	}

	while (!minseed.empty())
	{
		Seed.emplace(minseed.top());
		minseed.pop();
	}

	uglyTime("%d seeds generated", Seed.size());
}
// Initial commit (no-op)
