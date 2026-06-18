#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <limits>
#include <memory>
#include <mutex>
#include <queue>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>
#include <iostream>
#include <cstdarg>

#include <getopt.h>
#include <sys/time.h>
#include <unistd.h>
#include <sys/sysinfo.h>
#include <condition_variable>
#include <utility>
#include <future>
#include <functional>
#include <type_traits>
#include <fstream>
#include <sstream>
#include <cerrno>
#include <climits>
#include <cstdlib>

constexpr int MX = 10000;
constexpr double eps = 1e-6;
constexpr double MISSING_VALUE = std::numeric_limits<double>::quiet_NaN();
#define sameString(a, b) (strcmp((a), (b)) == 0)

using PII = std::pair<size_t, size_t>;

struct Array
{
	size_t con = 0;
	double val = MISSING_VALUE;
	double old_val = MISSING_VALUE;

	bool operator<(const Array& a) const
	{
		const bool this_nan = std::isnan(val);
		const bool other_nan = std::isnan(a.val);

		if (this_nan && other_nan) return con < a.con;
		if (this_nan) return false;
		if (other_nan) return true;

		if (std::fabs(val - a.val) < eps) return con < a.con;
		return val < a.val;
	}
};

struct Node
{
	size_t area = 0;
	size_t x1 = 0;
	size_t x2 = 0;
	size_t sig = 1;
	std::vector<size_t> in;
	std::vector<size_t> de;
	std::vector<size_t> index_column;
};

struct NoiseStatRow
{
	std::vector<double> all_vals;
	std::vector<double> down_vals;
	std::vector<double> up_vals;

	void clear()
	{
		all_vals.clear();
		down_vals.clear();
		up_vals.clear();
	}
};

using PIN = std::pair<PII, Node>;

struct SEED
{
	//bool flag = false;
	size_t x1 = 0;
	size_t x2 = 0;
	size_t len = 0;
};

using PSS = std::pair<SEED, SEED>;

struct cmp1
{
	bool operator()(const SEED& a, const SEED& b) const
	{
		if (a.len != b.len)
			return a.len > b.len;
		if (a.x1 != b.x1)
			return a.x1 < b.x1;
		return a.x2 < b.x2;
	}
};

struct cmp2
{
	bool operator()(const SEED& a, const SEED& b) const
	{
		if (a.len != b.len)
			return a.len < b.len;
		if (a.x1 != b.x1)
			return a.x1 > b.x1;
		return a.x2 > b.x2;
	}
};

struct Ans
{
	size_t len = 0;
	size_t sig = 0;
	std::vector<size_t> p;
};

struct PIA
{
	size_t id = 0;
	Ans ans;
};

struct Gene
{
	std::vector<size_t> c;
	std::vector<uint8_t> cluster_index;
};

struct SUM
{
	size_t index = 0;
	size_t Sum = 0;
};

struct LCS_NODE
{
	bool store_path = false;
	bool repeat_check[2] = {false, false};
	bool constant_check[2] = {false, false};
	std::vector<std::vector<std::vector<size_t>>> path;
	std::vector<std::vector<std::vector<size_t>>> dp;
	std::vector<std::vector<size_t>> sig1, sig2;
	std::vector<std::vector<Array>> a1, a2;
};

struct SeedLCSWorkspace
{
	bool repeat_check[2] = {false, false};
	bool constant_check[2] = {false, false};
	std::vector<std::vector<size_t>> sig1, sig2;
	std::vector<std::vector<Array>> a1, a2;
	std::vector<size_t> prev[2], cur[2];
};

struct SequenceRepeatStats
{
	bool repeat = false;
	bool constant = false;
};

struct Prog_options
{
	char FN[MX] = {0};        // input matrix
	char FP[MX] = {0};        // output prefix
	char FN_NOISE[MX] = {0};  // optional noise matrix file

	double LCS_TOLERANCE = 0.0;
	double DICHOTOMY_TOLERANCE = 0.0;
	double EXPAND_TOLERANCE = 0.0;
	double FILTER = 0.0;

	size_t COL_WIDTH = 0;
	size_t ROW_WIDTH = 0;
	size_t SEED_NUM = 0;
	size_t BLOCK_NUM = 0;
	size_t SEED_NUM_MULTIPLIER = 20;
	size_t SEED_SELECT_MODE = 0;
	size_t SEED_POOL_MULTIPLIER = 1;

	double QUANTILE = 0.0;
	size_t ABSOLUTE_QUANTILE = 0; 
	double MIN_LENGTH = 0.0;
	size_t THREADS_NUM = 0;
	size_t CLUSTER_WIDTH = 0;
	double CLUSTER_SIZE = 0.0;
	double DISCRETIZATION = 0.0;
	size_t IS_SINGLE_CELL_DATA = 0;
	size_t NORMALIZATION = 0;
};

// ==================== global variables ====================
extern std::priority_queue<SEED, std::vector<SEED>, cmp2> Seed;
extern std::vector<std::vector<Array>> A, B;
extern std::vector<std::vector<double>> NoiseMatrix;
extern size_t m, n;
extern std::vector<Node> BiCluster;
extern size_t BiCluster_num;
extern std::vector<std::string> conds, genes;
extern Prog_options* po;

// ==================== ThreadPool ====================
class ThreadPool
{
public:
	explicit ThreadPool(size_t numThreads) : stop(false)
	{
		for (size_t i = 0; i < numThreads; ++i)
		{
			workers.emplace_back([this]()
			{
				while (true)
				{
					std::function<void()> task;

					{
						std::unique_lock<std::mutex> lock(queueMutex);
						condition.wait(lock, [this]()
						{
							return stop || !taskQueue.empty();
						});

						if (stop && taskQueue.empty())
							return;

						task = std::move(taskQueue.front());
						taskQueue.pop();
					}

					task();
				}
			});
		}
	}

	template <class F, class... Args>
	auto enqueue(F&& f, Args&&... args)
		-> std::future<std::invoke_result_t<F, Args...>>
	{
		using return_type = std::invoke_result_t<F, Args...>;

		auto task = std::make_shared<std::packaged_task<return_type()>>(
			std::bind(std::forward<F>(f), std::forward<Args>(args)...)
		);

		std::future<return_type> result = task->get_future();

		{
			std::unique_lock<std::mutex> lock(queueMutex);
			if (stop)
				throw std::runtime_error("enqueue on stopped ThreadPool");

			taskQueue.emplace([task]() { (*task)(); });
		}

		condition.notify_one();
		return result;
	}

	~ThreadPool()
	{
		{
			std::unique_lock<std::mutex> lock(queueMutex);
			stop = true;
		}

		condition.notify_all();

		for (std::thread& worker : workers)
		{
			if (worker.joinable())
				worker.join();
		}
	}

private:
	std::vector<std::thread> workers;
	std::queue<std::function<void()>> taskQueue;
	std::mutex queueMutex;
	std::condition_variable condition;
	bool stop;
};

// ==================== utility ====================
void progress(const char* format, ...);
void err(const char* format, ...);
void errAbort(const char* format, ...);
long clock1000();
void uglyTime(const char* label, ...);
FILE* mustOpen(const char* filename, const char* mode);

// ==================== get_options ====================
void get_options(int argc, char* argv[]);

// ==================== data_processing ====================
void read_file(const char* filename);
void read_noise_matrix(const char* filename);
void data_preprocessing();
uint8_t which_side(const double& val, const double& d, const double& mid,
                  const size_t& j, const size_t& s, const size_t& t);


// ==================== seed_generation ====================
void seed_generation();
SEED SEED_GEN(size_t k, size_t j, SeedLCSWorkspace& workspace,
              const std::vector<uint8_t>& row_repeat,
              const std::vector<uint8_t>& row_constant);

// ==================== LCS ====================
Ans LCS(const std::vector<Array>& A1, const std::vector<Array>& A2,
        size_t n1, size_t n2,
        size_t r1, size_t r2,
        size_t is_print);
Ans LCS_seed(const std::vector<Array>& A1, const std::vector<Array>& A2,
             size_t n1, size_t n2,
             size_t r1, size_t r2,
             SeedLCSWorkspace& workspace,
             const std::vector<uint8_t>& row_repeat,
             const std::vector<uint8_t>& row_constant);
SequenceRepeatStats sequence_repeat_stats(const std::vector<Array> &a, size_t len, size_t row);

void solve(size_t id1, size_t id2, size_t sig,
           size_t i, size_t j,
           double ts1, double ts2,
           LCS_NODE& T);

double check2(size_t x1, size_t x2, size_t sig, double ts, const LCS_NODE& T);
double check1(size_t x1, size_t x2, size_t sig, double ts, const LCS_NODE& T);
void print_path(size_t i, size_t j, size_t k, size_t sig, const LCS_NODE& T, Ans& ans);

// ==================== cluster ====================
void cluster();
bool check(size_t a, size_t b, size_t c, size_t d, size_t Length);
PII overlapping_check(size_t x1, size_t x2);
bool cmp(const PIA& a, const PIA& b);
bool bicluster_check(size_t x1, size_t x2, const std::vector<Gene>& row_cluster);
size_t Max(size_t x1, size_t x2, const std::vector<Gene>& row_cluster);
bool is_seed(const SEED& S, const std::vector<Gene>& row_cluster);
PIA CLUSTER(const size_t& i, const std::vector<Array>& a, const size_t& x1);

// ==================== cluster_expand ====================
void cluster_expand();
void column_expand();
PSS CLUSTER_EXPAND(const size_t& j,
                   const std::vector<std::vector<Array>>& a,
                   const size_t& x1,
                   const size_t& x2,
                   bool use_x2);
size_t seed_cal(const Node& S, size_t col);

// ==================== result_output ====================
void result_output(char* _out);
PIN OverlapCal(size_t x1, size_t x2);
