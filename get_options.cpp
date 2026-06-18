
/***************************************************************************/

#include "structure.h"
using namespace std;
/***************************************************************************/
static const char USAGE[] =
    "===================================================================\n\
[Usage]\n\
$ ./noibic -i filename [argument list]\n\
===================================================================\n\
-V, --version :    Show program version and exit. \n\
===================================================================\n\
[Input]\n\
-i : The input file must be in one of two tab-delimited formats:\n\
     o        cond1    cond2    cond3\n\
     gene1      2.4      3.5     -2.4\n\
     gene2     -2.1      0.0      1.2\n\
-w : Optional noise file.\n\
     A tab-delimited numeric matrix without row/column labels, with the same\n\
     dimensions as the input. Each entry is the allowed relative fluctuation of\n\
     the corresponding expression value. If omitted, row-wise CV is used.\n\
     Example:\n\
          0.10    0.15    0.12\n\
          0.08    0.10    0.09\n\
===================================================================\n\
[Seeding]\n\
-M : Seed selection mode. 0 keeps the original top-k strategy; 1 enables row-diverse selection; 2 enables random selection; 3 keeps all qualified seeds.  \n\
	     Integer value 0, 1, 2, or 3; default: 1.  \n\
-k : Seed number multiplier. seed_num is derived as k * -o.  \n\
	     Positive integer, default: 20.  \n\
-u : Candidate seed pool multiplier. Each worker keeps up to u * seed_num candidates.  \n\
	     Positive integer, default: 5.  \n\
	     row-diverse gap uses rows / seed_num.  \n\
===================================================================\n\
[Data Preprocessing]\n\
-q : Remove non-expressed data based on numerical values  \n\
     Floating-point value in the range (0, 0.5], default: 0.1  \n\
     (see details in the Methods section of the paper).  \n\
-a : Remove non-expressed data by index position.  \n\
     Binary variable (0 or 1), default: 0  \n\
	 Recommended: set to 1 for matrices with more than 2000 rows.  \n\
     (see details in the Methods section of the paper).  \n\
-N : Normalize input data before biclustering.\n\
     Binary variable (0 or 1), default: 0.\n\
===================================================================\n\
[Biclustering]\n\
-l : Permissible mismatch rate during element swapping when searching  \n\
     for the Longest Common Subsequence (LCS) between two sequences.  \n\
     Floating-point value in the range [0, 1], default: 0.1.  \n\
-n : Minimum cluster width as a fraction of the original seed length.  \n\
     Floating-point value in the range [0, 1], default: 0.12.  \n\
-d : Threshold for filtering candidate seed sequences based on maximum redundancy.  \n\
     Floating-point value in the range [0, 1], default: 0.3.  \n\
-t : Number of threads for multi-threaded execution.  \n\
     Positive integer, default: 16.  \n\
===================================================================\n\
[Expansion]\n\
-e : Permissible error rate during bicluster expansion.  \n\
     Floating-point value in the range (0.5, 1], default: 0.85.  \n\
-b : Column/row expansion mode.  \n\
     Number in [0, 1] for column-expansion tolerance, 2 to disable column expansion, 3 to disable row expansion, or 4 to disable both.  \n\
     Default: 0.05.  \n\
===================================================================\n\
[Output]\n\
-c : Minimum number of columns in a block.  \n\
     Integer ≥ 3, default: 6.  \n\
-r : Minimum number of rows in a block.  \n\
     Integer ≥ 3, default: 4.  \n\
-m : Minimum block size as a fraction of the original matrix dimensions.\n\
    Floating-point value in the range [0, 1], default:0.0.\n\
-f : Overlap filtering threshold for biclusters.  \n\
     Floating-point value in the range [0, 1.0], default: 0.8.  \n\
     Smaller values remove more overlapping biclusters; 1.0 keeps almost all.  \n\
-S : Specify whether the input data is single-cell data.  \n\
     Binary variable (0 or 1), default: 0.  \n\
-o : Number of biclusters to report.  \n\
     Positive integer, default: 30.  \n\
-p : Output file name.  \n\
     Defaults to the input file name.  \n\
===================================================================\n";

static bool parse_double(const char *s, double &out)
{
  if (s == nullptr || *s == '\0')
    return false;

  char *end = nullptr;
  errno = 0;
  out = std::strtod(s, &end);
  return errno == 0 && end != s && *end == '\0';
}

static bool parse_size_t(const char *s, size_t &out)
{
  if (s == nullptr || *s == '\0' || *s == '-')
    return false;

  char *end = nullptr;
  errno = 0;
  unsigned long long v = std::strtoull(s, &end, 10);
  if (errno != 0 || end == s || *end != '\0')
    return false;
  if (v > static_cast<unsigned long long>(std::numeric_limits<size_t>::max()))
  {
    return false;
  }

  out = static_cast<size_t>(v);
  return true;
}

static bool parse_binary(const char *s, size_t &out)
{
  size_t v = 0;
  if (!parse_size_t(s, v))
    return false;
  if (v != 0 && v != 1)
    return false;
  out = v;
  return true;
}

template <size_t N>
static void copy_cstr(char (&dst)[N], const char *src)
{
  if (src == nullptr)
  {
    dst[0] = '\0';
    return;
  }
  std::snprintf(dst, N, "%s", src);
}

static void init_options()
{
  /* default parameters */

  po->FN[0] = '\0';
  po->FP[0] = '\0';
  po->FN_NOISE[0] = '\0';

  po->COL_WIDTH = 6;
  po->ROW_WIDTH = 4;

  po->QUANTILE = 0.1;
  po->ABSOLUTE_QUANTILE = 0;
  po->LCS_TOLERANCE = 0.1;
  po->DICHOTOMY_TOLERANCE = 0.05;
  po->EXPAND_TOLERANCE = 0.8;
  po->NORMALIZATION = 0;

  po->BLOCK_NUM = 30;
  po->SEED_NUM_MULTIPLIER = 20;
  po->SEED_NUM = po->SEED_NUM_MULTIPLIER * po->BLOCK_NUM;
  po->SEED_SELECT_MODE = 1;
  po->SEED_POOL_MULTIPLIER = 5;
  po->FILTER = 0.8;
  po->MIN_LENGTH = 0.0;
  po->THREADS_NUM = 16;
  po->IS_SINGLE_CELL_DATA = 0;
  po->CLUSTER_WIDTH = 12;

  po->CLUSTER_SIZE = 0.12;
  po->DISCRETIZATION = 0.3;
}

/*argc is a count of the arguments supplied to the program and argc[] is an
 * array of pointers to the strings which are those arguments-its type is array
 * of pointer to char
 */
void get_options(int argc, char *argv[])
{
  int op;
  bool is_valid = true;

  po = new Prog_options;
  /*Initialize the point*/
  // Prog_options *po;
  init_options();

  /*The getopt function gets the next option argument from the argument list
   *specified by the argv and argc arguments. Normally these values come
   *directly from the arguments received by main
   */
  /*An option character in this string can be followed by a colon (:) to
   *indicate that it takes a required argument. If an option character is
   *followed by two colons (::), its argument is optional if an option character
   *is followed by no colons, it does not need argument
   */
  while ((op = getopt(argc, argv, "i:w:p:l:b:d:q:a:e:f:c:r:o:m:t:S:n:N:M:k:u:h")) > 0)
  {
    switch (op)
    {
    /*optarg is set by getopt to point at the value of the option argument, for
     * those options that accept arguments*/
    case 'i':
      copy_cstr(po->FN, optarg);
      break;

    case 'w':
      copy_cstr(po->FN_NOISE, optarg);
      break;

    case 'p':
      copy_cstr(po->FP, optarg);
      break;
    case 'l':
      if (!parse_double(optarg, po->LCS_TOLERANCE))
      {
        err("-l should be a floating-point number");
        is_valid = false;
      }
      break;
    case 'b':
      if (!parse_double(optarg, po->DICHOTOMY_TOLERANCE))
      {
        err("-b should be a number: [0,1], 2, 3, or 4");
        is_valid = false;
      }
      break;
    case 'd':
      if (!parse_double(optarg, po->DISCRETIZATION))
      {
        err("-d should be a floating-point number");
        is_valid = false;
      }
      break;

    case 'M':
      if (!parse_size_t(optarg, po->SEED_SELECT_MODE))
      {
        err("-M should be an integer value of 0/1/2/3");
        is_valid = false;
      }
      break;

    case 'k':
      if (!parse_size_t(optarg, po->SEED_NUM_MULTIPLIER))
      {
        err("-k should be a positive integer");
        is_valid = false;
      }
      break;

    case 'u':
      if (!parse_size_t(optarg, po->SEED_POOL_MULTIPLIER))
      {
        err("-u should be a positive integer");
        is_valid = false;
      }
      break;

    case 'q':
      if (!parse_double(optarg, po->QUANTILE))
      {
        err("-q should be a floating-point number");
        is_valid = false;
      }
      break;
    case 'a':
      if (!parse_binary(optarg, po->ABSOLUTE_QUANTILE))
      {
        err("-a should be a binary value of 0/1");
        is_valid = false;
      }
      break;

    case 'e':
      if (!parse_double(optarg, po->EXPAND_TOLERANCE))
      {
        err("-e should be a floating-point number");
        is_valid = false;
      }
      break;
    case 'f':
      if (!parse_double(optarg, po->FILTER))
      {
        err("-f should be a floating-point number");
        is_valid = false;
      }
      break;

    case 'c':
      if (!parse_size_t(optarg, po->COL_WIDTH))
      {
        err("-c should be a positive integer");
        is_valid = false;
      }
      break;

    case 'r':
      if (!parse_size_t(optarg, po->ROW_WIDTH))
      {
        err("-r should be a positive integer");
        is_valid = false;
      }
      break;

    case 'o':
      if (!parse_size_t(optarg, po->BLOCK_NUM))
      {
        err("-o should be a positive integer");
        is_valid = false;
      }
      break;

    case 'm':
      if (!parse_double(optarg, po->MIN_LENGTH))
      {
        err("-m should be a floating-point number");
        is_valid = false;
      }
      break;

    case 't':
      if (!parse_size_t(optarg, po->THREADS_NUM))
      {
        err("-t should be a positive integer");
        is_valid = false;
      }
      break;

    case 'S':
      if (!parse_binary(optarg, po->IS_SINGLE_CELL_DATA))
      {
        err("-S should be a binary value of 0/1");
        is_valid = false;
      }
      break;

    case 'n':
      if (!parse_double(optarg, po->CLUSTER_SIZE))
      {
        err("-n should be a floating-point number");
        is_valid = false;
      }
      break;

    case 'N':
      if (!parse_binary(optarg, po->NORMALIZATION))
      {
        err("-N should be a binary value of 0/1");
        is_valid = false;
      }
      break;
    case 'h':
      puts(USAGE);
      exit(0);
    /*if expression does not match any constant-expression, control is
     * transferred to the statement(s) that follow the optional default label*/
    default:
      is_valid = false;
    }
  }
  /* basic sanity check */
  if (is_valid && po->FN[0] == '\0')
  {
    puts(USAGE);
    exit(0);
  }

  if (po->FP[0] == '\0')
  {
    copy_cstr(po->FP, po->FN);
  }

  if ((po->LCS_TOLERANCE > 1) || (po->LCS_TOLERANCE < 0))
  {
    err("-l noise ratio should be [0,1]");
    is_valid = false;
  }

  if (((po->DICHOTOMY_TOLERANCE > 1) || (po->DICHOTOMY_TOLERANCE < 0)) &&
      po->DICHOTOMY_TOLERANCE != 2 && po->DICHOTOMY_TOLERANCE != 3 &&
      po->DICHOTOMY_TOLERANCE != 4)
  {
    err("-b expansion mode should be [0,1], 2, 3, or 4");
    is_valid = false;
  }

  if ((po->DISCRETIZATION > 1) || (po->DISCRETIZATION < 0))
  {
    err("-d discretization ratio should be [0,1]");
    is_valid = false;
  }

  if (po->SEED_SELECT_MODE > 3)
  {
    err("-M seed selection mode should be 0, 1, 2, or 3");
    is_valid = false;
  }

  if (po->SEED_POOL_MULTIPLIER == 0)
  {
    err("-u seed pool multiplier should be >0");
    is_valid = false;
  }

  if (po->SEED_NUM_MULTIPLIER == 0)
  {
    err("-k seed number multiplier should be >0");
    is_valid = false;
  }

  if ((po->QUANTILE > 0.5) || (po->QUANTILE <= 0))
  {
    err("-q quantile discretization should be (0,.5]");
    is_valid = false;
  }

  if (po->ABSOLUTE_QUANTILE != 0 && po->ABSOLUTE_QUANTILE != 1)
  {
    err("-a should be a bool value of 0/1");
    is_valid = false;
  }

  if ((po->EXPAND_TOLERANCE > 1) || (po->EXPAND_TOLERANCE <= 0.5))
  {
    err("-e noise ratio should be (0.5,1]");
    is_valid = false;
  }

  if ((po->FILTER > 1) || (po->FILTER < 0))
  {
    err("-f overlapping filtering should be [0,1.]");
    is_valid = false;
  }

  if (po->COL_WIDTH < 3)
  {
    err("-c minimum column width should be >=3");
    is_valid = false;
  }

  if (po->ROW_WIDTH < 3)
  {
    err("-r minimum row width should be >=3");
    is_valid = false;
  }

  if (po->BLOCK_NUM <= 0)
  {
    err("-o number of blocks to report should be >0");
    is_valid = false;
  }

  if (po->BLOCK_NUM > 0 && po->SEED_NUM_MULTIPLIER > 0)
  {
    if (po->BLOCK_NUM > numeric_limits<size_t>::max() / po->SEED_NUM_MULTIPLIER)
    {
      err("-o * -k overflows seed_num");
      is_valid = false;
    }
    else
    {
      po->SEED_NUM = po->BLOCK_NUM * po->SEED_NUM_MULTIPLIER;
    }
  }

  if ((po->MIN_LENGTH > 1) || (po->MIN_LENGTH < 0))
  {
    err("-m should be [0,1]");
    is_valid = false;
  }

  if (po->THREADS_NUM <= 0)
  {
    err("-t number of threads should be a positive integer");
    is_valid = false;
  }

  if (po->IS_SINGLE_CELL_DATA != 0 && po->IS_SINGLE_CELL_DATA != 1)
  {
    err("-S should be a bool value of 0/1");
    is_valid = false;
  }

  if ((po->CLUSTER_SIZE > 1) || (po->CLUSTER_SIZE < 0))
  {
    err("-n cluster column ratio should be [0,1]");
    is_valid = false;
  }

  if (po->NORMALIZATION != 0 && po->NORMALIZATION != 1)
  {
    err("-N should be a bool value of 0/1");
    is_valid = false;
  }

  if (!is_valid)
    errAbort("Type -h to view possible options");
}
/***************************************************************************/
// Initial commit (no-op)
