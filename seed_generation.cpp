#include "structure.h"
#include <atomic>
#include <cstdint>
#include <iterator>
#include <memory>
#include <random>
using namespace std;

priority_queue<SEED, vector<SEED>, cmp2> Seed;

static void build_seed_row_cache(vector<uint8_t>& row_repeat, vector<uint8_t>& row_constant)
{
	row_repeat.assign(n, 0);
	row_constant.assign(n, 0);

	for (size_t i = 0; i < n; ++i)
	{
		size_t len = 0;
		while (len < m && !isnan(A[i][len].val))
			++len;
		SequenceRepeatStats stats = sequence_repeat_stats(A[i], len, i);
		row_repeat[i] = stats.repeat;
		row_constant[i] = stats.constant;
	}
}

SEED SEED_GEN(size_t k, size_t j,
			  SeedLCSWorkspace& workspace,
			  const vector<uint8_t>& row_repeat,
			  const vector<uint8_t>& row_constant)
{
	Ans ans = LCS_seed(A[k], A[j], m, m, k, j, workspace, row_repeat, row_constant);

	SEED seed{};
	seed.len = ans.len;
	seed.x1 = k;
	seed.x2 = j;
	return seed;
}

static inline void keep_top_seed(priority_queue<SEED, vector<SEED>, cmp1>& top_seed, SEED&& seed, size_t limit)
{
	if (limit == 0 || seed.len < po->COL_WIDTH)
		return;

	if (top_seed.size() == limit)
	{
		const SEED& worst = top_seed.top();
		if (seed.len < worst.len ||
			(seed.len == worst.len &&
			 (seed.x1 > worst.x1 || (seed.x1 == worst.x1 && seed.x2 >= worst.x2))))
			return;
		top_seed.pop();
	}

	top_seed.emplace(std::move(seed));
}

static size_t seed_row_cap()
{
	return max<size_t>(1, po->BLOCK_NUM);
}

static size_t seed_row_gap()
{
	if (po->SEED_NUM == 0)
		return 0;
	return n / po->SEED_NUM;
}

static size_t seed_pool_limit()
{
	if (po->SEED_NUM != 0 &&
		po->SEED_POOL_MULTIPLIER > numeric_limits<size_t>::max() / po->SEED_NUM)
		errAbort("-u * seed_num overflows seed pool limit");
	return po->SEED_NUM * po->SEED_POOL_MULTIPLIER;
}

static bool seed_better(const SEED& a, const SEED& b)
{
	if (a.len != b.len)
		return a.len > b.len;
	if (a.x1 != b.x1)
		return a.x1 < b.x1;
	return a.x2 < b.x2;
}

static mt19937_64 make_seed_rng()
{
	static atomic<uint64_t> counter{0};
	random_device rd;
	const uint64_t c = counter.fetch_add(1, memory_order_relaxed);
	const uint64_t tid = static_cast<uint64_t>(hash<thread::id>{}(this_thread::get_id()));
	seed_seq seq{
		rd(), rd(),
		static_cast<unsigned>(clock1000()),
		static_cast<unsigned>(c),
		static_cast<unsigned>(c >> 32),
		static_cast<unsigned>(tid),
		static_cast<unsigned>(tid >> 32)
	};
	return mt19937_64(seq);
}

static void keep_random_seed(vector<SEED>& random_seeds, SEED&& seed,
							 size_t limit, size_t& seen,
							 mt19937_64& rng)
{
	if (limit == 0 || seed.len < po->COL_WIDTH)
		return;

	++seen;
	if (random_seeds.size() < limit)
	{
		random_seeds.emplace_back(std::move(seed));
		return;
	}

	uniform_int_distribution<size_t> dist(0, seen - 1);
	const size_t pos = dist(rng);
	if (pos < limit)
		random_seeds[pos] = std::move(seed);
}

static vector<SEED> select_random_seeds(vector<SEED>& candidates, size_t limit)
{
	if (limit == 0 || candidates.empty())
		return {};

	mt19937_64 rng = make_seed_rng();
	shuffle(candidates.begin(), candidates.end(), rng);

	const size_t selected_size = min(limit, candidates.size());
	vector<SEED> selected;
	selected.reserve(selected_size);
	for (size_t i = 0; i < selected_size; ++i)
		selected.emplace_back(std::move(candidates[i]));

	return selected;
}

class RowGapIndex
{
public:
	explicit RowGapIndex(size_t row_count)
		: blocked(row_count, 0)
	{
	}

	bool can_anchor(size_t row) const
	{
		return row < blocked.size() && !blocked[row];
	}

	void anchor(size_t row, size_t gap)
	{
		if (row >= blocked.size() || blocked[row])
			return;

		const size_t last = blocked.size() - 1;
		const size_t left = row > gap ? row - gap : 0;
		const size_t right = row + min(gap, last - row);
		fill(blocked.begin() + left, blocked.begin() + right + 1, uint8_t{1});
	}

private:
	vector<uint8_t> blocked;
};

static bool seed_respects_row_cap(const SEED& seed, const vector<size_t>& row_counts)
{
	const size_t row_cap = seed_row_cap();
	if (row_counts[seed.x1] >= row_cap || row_counts[seed.x2] >= row_cap)
		return false;
	return true;
}

static vector<SEED> select_diverse_seeds(vector<SEED>& candidates, size_t limit)
{
	sort(candidates.begin(), candidates.end(), seed_better);

	vector<SEED> selected;
	selected.reserve(min(limit, candidates.size()));
	vector<size_t> row_counts(n, 0);
	RowGapIndex anchors(n);
	const size_t gap = seed_row_gap();

	for (const SEED& seed : candidates)
	{
		if (selected.size() >= limit)
			break;
		if (!seed_respects_row_cap(seed, row_counts))
			continue;

		bool x1_anchor = anchors.can_anchor(seed.x1);
		bool x2_anchor = anchors.can_anchor(seed.x2);
		if (!x1_anchor && !x2_anchor)
			continue;

		selected.emplace_back(seed);
		++row_counts[seed.x1];
		++row_counts[seed.x2];
		if (x1_anchor)
			anchors.anchor(seed.x1, gap);
		if (x2_anchor && anchors.can_anchor(seed.x2))
			anchors.anchor(seed.x2, gap);
	}

	return selected;
}

static vector<SEED> drain_seed_queue(priority_queue<SEED, vector<SEED>, cmp1>& seeds)
{
	vector<SEED> result;
	result.reserve(seeds.size());
	while (!seeds.empty())
	{
		result.emplace_back(std::move(seeds.top()));
		seeds.pop();
	}
	return result;
}

static vector<SEED> seed_generation_worker(shared_ptr<atomic<size_t>> next_k,
										   size_t group_end,
										   size_t seed_limit,
										   shared_ptr<const vector<uint8_t>> row_repeat,
										   shared_ptr<const vector<uint8_t>> row_constant)
{
	priority_queue<SEED, vector<SEED>, cmp1> local_seed;
	vector<SEED> local_random_seeds;
	vector<SEED> local_candidates;
	SeedLCSWorkspace workspace;
	size_t random_seen = 0;
	mt19937_64 rng;
	if (po->SEED_SELECT_MODE == 2)
		rng = make_seed_rng();

	while (true)
	{
		size_t k = next_k->fetch_add(1, memory_order_relaxed);
		if (k >= group_end)
			break;

		for (size_t j = k + 1; j < group_end; ++j)
		{
			SEED seed = SEED_GEN(k, j, workspace, *row_repeat, *row_constant);
			if (seed.len < po->COL_WIDTH)
				continue;

			if (po->SEED_SELECT_MODE == 1 || po->SEED_SELECT_MODE == 3)
				local_candidates.emplace_back(std::move(seed));
			else if (po->SEED_SELECT_MODE == 2)
				keep_random_seed(local_random_seeds, std::move(seed), seed_limit, random_seen, rng);
			else
				keep_top_seed(local_seed, std::move(seed), seed_limit);
		}
	}

	if (po->SEED_SELECT_MODE == 1)
		return select_diverse_seeds(local_candidates, seed_limit);
	if (po->SEED_SELECT_MODE == 2)
		return local_random_seeds;
	if (po->SEED_SELECT_MODE == 3)
		return local_candidates;

	vector<SEED> seeds = drain_seed_queue(local_seed);
	return seeds;
}

void seed_generation()
{
	Seed = priority_queue<SEED, vector<SEED>, cmp2>();

	size_t group_size = static_cast<size_t>(ceil(n / 4.0));
	group_size = max<size_t>(group_size, 2);

	size_t thread_num = max<size_t>(po->THREADS_NUM, 1);
	ThreadPool pool(thread_num);
	vector<future<vector<SEED>>> results;
	results.reserve(thread_num * 4);
	auto row_repeat = make_shared<vector<uint8_t>>();
	auto row_constant = make_shared<vector<uint8_t>>();
	build_seed_row_cache(*row_repeat, *row_constant);

	const size_t seed_limit = po->SEED_NUM;
	const size_t candidate_limit = seed_pool_limit();

	priority_queue<SEED, vector<SEED>, cmp1> minseed;
	vector<SEED> global_candidates;

	for (size_t start = 0; start < n; start += group_size)
	{
		size_t end = min(n, start + group_size);
		size_t rows = end - start;
		if (rows < 2)
			continue;

		size_t worker_num = min(thread_num, rows);
		auto next_k = make_shared<atomic<size_t>>(start);
		for (size_t i = 0; i < worker_num; ++i)
			results.emplace_back(pool.enqueue(seed_generation_worker, next_k, end, candidate_limit, row_repeat, row_constant));
	}

	for (auto& fut : results)
	{
		vector<SEED> local_seeds = fut.get();
		if (po->SEED_SELECT_MODE == 1 || po->SEED_SELECT_MODE == 2 || po->SEED_SELECT_MODE == 3)
		{
			global_candidates.insert(global_candidates.end(),
									 make_move_iterator(local_seeds.begin()),
									 make_move_iterator(local_seeds.end()));
			continue;
		}

		for (SEED& seed : local_seeds)
			keep_top_seed(minseed, std::move(seed), seed_limit);
	}

	if (po->SEED_SELECT_MODE == 1)
	{
		vector<SEED> selected = select_diverse_seeds(global_candidates, po->SEED_NUM);
		for (SEED& seed : selected)
			Seed.emplace(std::move(seed));

		uglyTime("%zu seeds generated", Seed.size());
		return;
	}

	if (po->SEED_SELECT_MODE == 2)
	{
		vector<SEED> selected = select_random_seeds(global_candidates, po->SEED_NUM);
		for (SEED& seed : selected)
			Seed.emplace(std::move(seed));

		uglyTime("%zu seeds generated", Seed.size());
		return;
	}

	if (po->SEED_SELECT_MODE == 3)
	{
		for (SEED& seed : global_candidates)
			Seed.emplace(std::move(seed));

		uglyTime("%zu seeds generated", Seed.size());
		return;
	}

	while (!minseed.empty())
	{
		Seed.emplace(std::move(minseed.top()));
		minseed.pop();
	}

	uglyTime("%zu seeds generated", Seed.size());
}
