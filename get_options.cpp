
/***************************************************************************/

#include "structure.h"

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
===================================================================\n\
[Data Preprocessing]\n\
-q : Remove non-expressed data based on numerical values  \n\
     Floating-point value in the range (0, 0.5], default: 0.06  \n\
     (see details in the Methods section of the paper).  \n\
-a : Remove non-expressed data by index position.  \n\
     Binary variable (0 or 1), default: 0  \n\
     (see details in the Methods section of the paper).  \n\
-s : Replace zero values with random numbers drawn from N(0,1).  \n\
     Binary variable (0 or 1), default: 0.  \n\
-z : Exclude zero values from clustering.  \n\
     Binary variable (0 or 1), default: 0.  \n\
===================================================================\n\
[Biclustering]\n\
-l : Permissible mismatch rate during element swapping when searching  \n\
     for the Longest Common Subsequence (LCS) between two sequences.  \n\
     Floating-point value in the range [0, 1], default: 0.2.  \n\
-n : Minimum cluster width as a fraction of the original seed length.  \n\
     Floating-point value in the range [0, 1], default: 0.12.  \n\
-d : Threshold for filtering candidate seed sequences based on maximum redundancy.  \n\
     Floating-point value in the range [0, 1], default: 0.9.  \n\
-t : Number of threads for multi-threaded execution.  \n\
     Positive integer, default: 16.  \n\
===================================================================\n\
[Expansion]\n\
-e : Permissible error rate during bicluster expansion.  \n\
     Floating-point value in the range (0.5, 1], default: 0.85.  \n\
-b : Permissible numerical error for column expansion using binary search.  \n\
     Floating-point value in the range [0, 1], or 2 to disable column expansion.  \n\
     Default: 0.2.  \n\
===================================================================\n\
[Output]\n\
-c : Minimum number of columns in a block.  \n\
     Integer ≥ 3, default: 6.  \n\
-r : Minimum number of rows in a block.  \n\
     Integer ≥ 3, default: 8.  \n\
-m : Minimum block size as a fraction of the original matrix dimensions.\n\
    Floating-point value in the range [0, 1], default:0.0.\n\
-f : Overlap filtering threshold for biclusters.  \n\
     Floating-point value in the range [0, 1.0], default: 1 (no filtering).  \n\
-S : Specify whether the input data is single-cell data.  \n\
     Binary variable (0 or 1), default: 0.  \n\
-o : Number of biclusters to report.  \n\
     Positive integer, default: 30.  \n\
-p : Output file name.  \n\
     Defaults to the input file name.  \n\
===================================================================\n";

static void init_options() {
  /* default parameters */
  
  strcpy(po->FN, " ");
  strcpy(po->FP, " ");
  po->COL_WIDTH = 6;
  po->ROW_WIDTH=4;
  
  po->QUANTILE = 0.1;
  po->ABSOLUTE_QUANTILE=0;
  po->LCS_TOLERANCE = 0.2;
  po->DICHOTOMY_TOLERANCE=0.2;
  po->EXPAND_TOLERANCE=0.85;
  
  po->BLOCK_NUM = 30;
  po->SEED_NUM = 500;
  po->FILTER = 0.8;
  po->MIN_LENGTH=0.0;	  
  po->THREADS_NUM=16;
  po->IS_SINGLE_CELL_DATA=0;
  po->SINGLE_CELL_PROCESSING=0;
  po->ZERO_CLEAN=0;
  po->CLUSTER_WIDTH=12;
  po->CLUSTER_SIZE=0.12;
  po->DISCRETIZATION=0.9;
}

/*argc is a count of the arguments supplied to the program and argc[] is an
 * array of pointers to the strings which are those arguments-its type is array
 * of pointer to char
 */
void get_options(int argc, char *argv[]) {
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
  while ((op = getopt(argc, argv, "i:p:l:b:d:q:a:e:f:c:r:o:m:t:s:S:z:n:h")) > 0) {
    switch (op) {
    /*optarg is set by getopt to point at the value of the option argument, for
     * those options that accept arguments*/
    case 'i':
      strcpy(po->FN, optarg);
      break;
    case 'p':
      strcpy(po->FP, optarg);
      break;
    case 'l':
      po->LCS_TOLERANCE=atof(optarg);
      break;
    case 'b':
      po->DICHOTOMY_TOLERANCE=atof(optarg);
      break;
    case 'd':
      po->DISCRETIZATION=atof(optarg);
      break;
    /*atof can convert string to double*/
    case 'q':
      po->QUANTILE = atof(optarg);
      break;
    case 'a':
      po->ABSOLUTE_QUANTILE = atoi(optarg);
      break;
    /*atoi can convert string to integer*/
    case 'e':
      po->EXPAND_TOLERANCE = atof(optarg);
      break;
    case 'f':
      po->FILTER = atof(optarg);
      break;
    case 'c':
      po->COL_WIDTH = atoi(optarg);
      break;
    case 'r':
      po->ROW_WIDTH = atoi(optarg);
      break;
    case 'o':
      po->BLOCK_NUM = atoi(optarg);
      po->SEED_NUM = 20*po->BLOCK_NUM;
      break;
    case 'm':
      po->MIN_LENGTH = atof(optarg);
      break;
    case 't':
      po->THREADS_NUM=atoi(optarg);
      break;
    case 's':
      po->SINGLE_CELL_PROCESSING=atoi(optarg);
      break;
    case 'S':
      po->IS_SINGLE_CELL_DATA=atoi(optarg);
      break;
    case 'z':
      po->ZERO_CLEAN=atoi(optarg);
      break;
    case 'n':
      po->CLUSTER_SIZE=atof(optarg);
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
  if (is_valid && po->FN[0] == ' ') {
    puts(USAGE);
    exit(0);
  }
  /* option value range check */
  if(po->FP[0] == ' ')
    strcpy(po->FP,po->FN);

  if ((po->LCS_TOLERANCE > 1) || (po->LCS_TOLERANCE < 0)) {
    err("-l noise ratio should be [0,1]");
    is_valid = false;
  }

  if (((po->DICHOTOMY_TOLERANCE > 1) || (po->DICHOTOMY_TOLERANCE < 0))&&po->DICHOTOMY_TOLERANCE!=2) {
    err("-b dichotomy noise ratio should be [0,1] or 2");
    is_valid = false;
  }

   if ((po->DISCRETIZATION > 1) || (po->DISCRETIZATION < 0)) {
    err("-d discretization ratio should be [0,1]");
    is_valid = false;
  }

  if ((po->QUANTILE > .5) || (po->QUANTILE <= 0)) {
    err("-q quantile discretization should be (0,.5]");
    is_valid = false;
  }

  if (po->ABSOLUTE_QUANTILE!=0&&po->ABSOLUTE_QUANTILE!=1) {
    err("-a should be a bool value of 0/1");
    is_valid = false;
  }

  if ((po->EXPAND_TOLERANCE > 1) || (po->EXPAND_TOLERANCE <= 0.5)) {
    err("-e noise ratio should be (0.5,1]");
    is_valid = false;
  }

  if ((po->FILTER > 1) || (po->FILTER < 0)) {
    err("-f overlapping filtering should be [0,1.]");
    is_valid = false;
  }
 
  if (po->COL_WIDTH < 3 && po->COL_WIDTH != -1) {
    err("-c minimum column width should be >=3");
    is_valid = false;
  }

  if (po->ROW_WIDTH < 3 && po->ROW_WIDTH != -1) {
    err("-r minimum row width should be >=3");
    is_valid = false;
  }
  
  if (po->BLOCK_NUM <= 0) {
    err("-o number of blocks to report should be >0");
    is_valid = false;
  }
  
  if ((po->MIN_LENGTH > 1) || (po->MIN_LENGTH < 0)) {
    err("-m should be [0,1]");
    is_valid = false;
  }

  if (po->THREADS_NUM !=static_cast<int>(po->THREADS_NUM)) {
    err("-t number of threads should be a positive integer");
    is_valid = false;
  }
  
  if (po->SINGLE_CELL_PROCESSING!=0&&po->SINGLE_CELL_PROCESSING!=1) {
    err("-s should be a bool value of 0/1");
    is_valid = false;
  }

  if (po->IS_SINGLE_CELL_DATA!=0&&po->IS_SINGLE_CELL_DATA!=1) {
    err("-S should be a bool value of 0/1");
    is_valid = false;
  }

  if (po->ZERO_CLEAN!=0&&po->ZERO_CLEAN!=1) {
    err("-z should be a bool value of 0/1");
    is_valid = false;
  }
  
  if ((po->CLUSTER_SIZE > 1) || (po->CLUSTER_SIZE < 0)) {
    err("-n cluster column ratio should be [0,1]");
    is_valid = false;
  }

  if (!is_valid)
    errAbort("Type -h to view possible options");
}
/***************************************************************************/
// Initial commit (no-op)
