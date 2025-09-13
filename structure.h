#pragma once
#include <bits/stdc++.h>
#include <sys/time.h>
#include <unistd.h>
#include <getopt.h>
#include <sys/sysinfo.h>
#include <condition_variable>
#include <utility>
#include <future>
#include <functional>

using namespace std;
const int MX=1e4;
const double eps = 1e-6;

typedef pair<size_t, size_t> PII;

#define sameString(a, b) (strcmp((a), (b)) == 0)
#define INF std::numeric_limits<double>::infinity()
struct Array
{
	size_t con;
	double val;
	double old_val;
	bool operator<(const Array &a) const
	{
		if (fabs(val - a.val) < eps)
			return con < a.con;
		return val < a.val;
	}
};

struct Node
{
	size_t area;
	size_t x1, x2;
	vector<size_t> in, de;
	vector<size_t> index_column;
};

typedef pair<size_t, Node> PIN;

struct SEED
{
	bool flag;
	size_t x1, x2;
	size_t len;
};

typedef pair<SEED,SEED> PSS;

struct cmp1
{
	bool operator()(SEED a, SEED b)
	{
		return a.len > b.len;
	}
};

struct cmp2
{
	bool operator()(SEED a, SEED b)
	{
		return a.len < b.len;
	}
};

struct Ans
{
	size_t repeat_check;
	size_t len;
	size_t sig;
	vector<size_t> p;
};

struct PIA
{
	size_t id;
	Ans ans;
};

struct Gene
{
	vector<size_t> c;
	vector<uint8_t> cluster_index;
};

struct SUM
{
	size_t index;
	size_t Sum;
};

struct LCS_NODE
{
	
	vector<vector<vector<size_t> > > path;
	vector<vector<vector<size_t> > > dp;
	vector<unordered_map<size_t, size_t> > sig1, sig2;
	vector<vector<Array> > a1, a2;
	size_t path_cnt;
};

struct Prog_options
{
	char FN[MX];
	char FP[MX];
	double LCS_TOLERANCE;
	double DICHOTOMY_TOLERANCE;
	double EXPAND_TOLERANCE;
	double FILTER;
	size_t COL_WIDTH;
	size_t ROW_WIDTH;
	size_t SEED_NUM;
	size_t BLOCK_NUM;
	double QUANTILE;
	size_t ABSOLUTE_QUANTILE;
	double MIN_LENGTH;
	size_t THREADS_NUM;
	size_t SINGLE_CELL_PROCESSING;
	size_t ZERO_CLEAN;
	size_t CLUSTER_WIDTH;
	double CLUSTER_SIZE;
	double DISCRETIZATION;
	size_t IS_SINGLE_CELL_DATA;
};

// global var
extern priority_queue<SEED, vector<SEED>, cmp2> Seed;
extern vector<vector<Array>> A, B;
extern size_t m, n;
extern vector<Node> BiCluster;
extern size_t BiCluster_num;
extern vector<string> conds, genes;
extern Prog_options *po;
extern vector<size_t> Start_pos;

// pthreadpool
class ThreadPool
{
public:
	ThreadPool(size_t numThreads) : stop(false)
	{
		for (size_t i = 0; i < numThreads; ++i)
		{
			workers.emplace_back([this]
								 {
                while (true) {
                    function<void()> task;

                    {
                        unique_lock<mutex> lock(queueMutex);
                        condition.wait(lock, [this] { return stop || !taskQueue.empty(); });

                        if (stop && taskQueue.empty()) {
                            return;
                        }

                        task = std::move(taskQueue.front());
                        taskQueue.pop();
                    }

                    task();
                } });
		}
	}

	template <class F, class... Args>
	auto enqueue(F &&f, Args &&...args) -> future<typename result_of<F(Args...)>::type>
	{
		using return_type = typename result_of<F(Args...)>::type;

		auto task = make_shared<packaged_task<return_type()>>(bind(forward<F>(f), forward<Args>(args)...));

		future<return_type> result = task->get_future();

		{
			unique_lock<mutex> lock(queueMutex);

			if (stop)
			{
				throw runtime_error("enqueue on stopped ThreadPool");
			}

			taskQueue.emplace([task]()
							  { (*task)(); });
		}

		condition.notify_one();
		return result;
	}

	~ThreadPool()
	{
		{
			unique_lock<mutex> lock(queueMutex);
			stop = true;
		}

		condition.notify_all();

		for (thread &worker : workers)
		{
			worker.join();
		}
	}

private:
	vector<thread> workers;
	queue<function<void()>> taskQueue;

	mutex queueMutex;
	condition_variable condition;
	bool stop;
};

// struct
void progerss(const char *format, ...);
void err(const char *format, ...);
void errAbort(const char *format, ...);
long clock1000();
void uglyTime(const char *label, ...);
FILE *mustOpen(const char *filename, const char *mode);

// get_options
void get_options(int argc, char *argv[]);

// data_processing
void data_preprocessing();
size_t Lower_bound_Array(size_t l, size_t r, double val, const vector<Array> &a);
inline double mapToRange(double value, double mean, double stddev, double minValue, double maxValue);
size_t Upper_bound_Array(size_t l, size_t r, double val, const vector<Array> &a);
bool check_val(const double& val,const double& d,const double& mid,const size_t& j,const size_t& s,const size_t& t);

// seed_generation
void seed_generation();
SEED SEED_GEN(const size_t &k, const size_t &j);

// LCS
Ans LCS(const vector<Array> &A1, const vector<Array> &A2, size_t n1, size_t n2, size_t is_print, size_t is_seed_gen = 0, size_t st1 = 0, size_t st2 = 0,size_t is_repeat_check=0);
size_t solve(size_t id1, size_t id2, size_t sig, size_t i, size_t j, LCS_NODE &T);
double check2(size_t x1, size_t x2, size_t sig, const LCS_NODE &T);
double check1(size_t x1, size_t x2, size_t sig, const LCS_NODE &T);
inline bool check_val(double a, double b, double ratio);
void print_path(size_t i, size_t j, size_t k, size_t sig, const LCS_NODE &T, Ans &ans);

// cluster
void cluster();
bool check(size_t a, size_t b, size_t c, size_t d,size_t Length);
PII overlapping_check(size_t x1, size_t x2);
bool cmp(const PIA &a, const PIA &b);
bool bicluster_check(size_t x1, size_t x2,const vector<Gene> &row_cluster);
size_t Max(size_t x1, size_t x2,const vector<Gene> &row_cluster);
bool is_seed(const SEED &S,const vector<Gene> &row_cluster);
PIA CLUSTER(const size_t &i,const vector<Array> &a);
//int repeat_check(const SEED &S,const Ans &ans);

// cluster_expand
void cluster_expand();
void column_expand();
void row_expand();
PSS CLUSTER_EXPAND(const size_t &j, const vector<vector<Array>> &a);

// result_output
void result_output(char *_out);
PIN OverlapCal(size_t x1, size_t x2);

size_t seed_cal(Node S, size_t col);
// Initial commit (no-op)
