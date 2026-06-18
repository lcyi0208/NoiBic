#include "structure.h"
#include <random>
#include <cmath>
using namespace std;

struct min_max
{
    double Min, Max;
};
static inline bool is_missing(double x)
{
    return isnan(x);
}

static double cv_mean_abs(const vector<double> &vals)
{
    if (vals.empty())
        return MISSING_VALUE;
    if (vals.size() == 1)
        return 0.0;

    double mean = 0.0, meanAbs = 0.0;
    for (double x : vals)
    {
        mean += x;
        meanAbs += fabs(x);
    }
    mean /= (double)vals.size();
    meanAbs /= (double)vals.size();

    double var = 0.0;
    for (double x : vals)
    {
        double d = x - mean;
        var += d * d;
    }
    var /= (double)vals.size();

    double sd = sqrt(var);
    double eps_local = 1e-12;
    return sd / (meanAbs + eps_local);
}

static void record_noise_tag_and_value(size_t row,
                                       size_t col,
                                       double val,
                                       uint8_t tag,
                                       NoiseStatRow &stat)
{
    if (is_missing(val))
    {
        NoiseMatrix[row][col] = MISSING_VALUE;
        return;
    }

    if (po->QUANTILE >= 0.5)
    {
        stat.all_vals.push_back(val);
        return;
    }

    if (tag == 1)
    {
        NoiseMatrix[row][col] = 1; 
        stat.down_vals.push_back(val);
    }
    else if (tag == 2)
    {
        NoiseMatrix[row][col] = 2; 
        stat.up_vals.push_back(val);
    }
    else
    {
        NoiseMatrix[row][col] = MISSING_VALUE;
    }
}

static void finalize_noise_row_cv(size_t row, const NoiseStatRow &stat)
{
    if (po->QUANTILE >= 0.5)
    {
        double cv_all = cv_mean_abs(stat.all_vals);
        for (size_t j = 0; j < m; ++j)
        {
            if (is_missing(B[row][j].val))
                NoiseMatrix[row][j] = MISSING_VALUE;
            else
                NoiseMatrix[row][j] = cv_all;
        }
    }
    else
    {
        double cv_down = cv_mean_abs(stat.down_vals);
        double cv_up = cv_mean_abs(stat.up_vals);

        for (size_t j = 0; j < m; ++j)
        {
            if (is_missing(B[row][j].val) || isnan(NoiseMatrix[row][j]))
            {
                NoiseMatrix[row][j] = MISSING_VALUE;
            }
            else if (NoiseMatrix[row][j] == 1)
            {
                NoiseMatrix[row][j] = cv_down;
            }
            else if (NoiseMatrix[row][j] == 2)
            {
                NoiseMatrix[row][j] = cv_up;
            }
            else
            {
                NoiseMatrix[row][j] = MISSING_VALUE;
            }
        }
    }
}

uint8_t which_side(const double &val, const double &d, const double &mid,
                   const size_t &j, const size_t &s, const size_t &t)
{
    if (po->ABSOLUTE_QUANTILE == 1)
    {
        if (j > s && j < t)
            return 0;
        if (j <= s)
            return 1;
        if (j >= t)
            return 2;
    }

    if (val < mid + d && val > mid - d)
        return 0;

    if (j > s && j < t)
    {
        if (fabs(val - (mid + d)) < eps || fabs(val - (mid - d)) < eps)
            return 0;
    }
    if (val <= mid - d)
        return 1;
    if (val >= mid + d)
        return 2;
    return 0;
}

void data_preprocessing()
{
    vector<min_max> mm(n);

    for (size_t i = 0; i < n; i++)
    {
        sort(A[i].begin(), A[i].end());

        mm[i].Min = A[i][0].val;
        mm[i].Max = A[i][m - 1].val;
    }

    if (po->QUANTILE < 0.5)
    {
        size_t l = m / 2;
        size_t s = (size_t)(po->QUANTILE * m);
        size_t t = (size_t)((1 - po->QUANTILE) * m);

        double d, mid;

        for (size_t i = 0; i < n; i++)
        {
            NoiseStatRow stat;

            mid = A[i][l].val;
            d = min(A[i][l].val - A[i][s].val, A[i][t].val - A[i][l].val);

            for (size_t j = 0; j < m; j++)
            {
                uint8_t ckv = which_side(A[i][j].val, d, mid, j, s, t);
                size_t col = A[i][j].con;

                if (!ckv)
                {
                    A[i][j].val = MISSING_VALUE;
                    A[i][j].old_val = MISSING_VALUE;
                }
                else if (po->NORMALIZATION)
                {
                    double tmp_val = (fabs(mm[i].Max - mm[i].Min) < eps)
                                         ? 1.0
                                         : (mm[i].Max - mm[i].Min);

                    A[i][j].val = (A[i][j].val - mm[i].Min) / tmp_val;
                    A[i][j].old_val = A[i][j].val;
                }

                B[i][col].con = col;
                B[i][col].val = A[i][j].val;
                B[i][col].old_val = A[i][j].old_val;

                if (po->FN_NOISE[0] == '\0')
                    record_noise_tag_and_value(i, col, B[i][col].val, ckv, stat);
            }

            if (po->FN_NOISE[0] == '\0')
                finalize_noise_row_cv(i, stat);

            sort(A[i].begin(), A[i].end());
        }
    }
    else
    {
        for (size_t i = 0; i < n; i++)
        {
            NoiseStatRow stat;

            for (size_t j = 0; j < m; j++)
            {
                size_t col = A[i][j].con;

                if (po->NORMALIZATION)
                {
                    double tmp_val = (fabs(mm[i].Max - mm[i].Min) < eps)
                                         ? 1.0
                                         : (mm[i].Max - mm[i].Min);

                    A[i][j].val = (A[i][j].val - mm[i].Min) / tmp_val;
                    A[i][j].old_val = A[i][j].val;
                }

                B[i][col].con = col;
                B[i][col].val = A[i][j].val;
                B[i][col].old_val = A[i][j].old_val;

                if (po->FN_NOISE[0] == '\0')
                    record_noise_tag_and_value(i, col, B[i][col].val, 3, stat);
            }

            if (po->FN_NOISE[0] == '\0')
                finalize_noise_row_cv(i, stat);

            sort(A[i].begin(), A[i].end());
        }
    }

    /*cout << "data_preprcessing" << endl;
    for (int i = 0; i < (int)n; i++)
    {
        cout << genes[i] << ":\t";
        for (int j = 0; j < (int)m; j++)
        {
            if (isnan(A[i][j].val))
                cout << "NaN";
            else
                cout << A[i][j].val;
            cout << "-" << A[i][j].con << "\t";
        }
        cout << endl;
    }
    cout << endl;
    cout<<"B matrix:"<<endl;
    for (int i = 0; i < (int)n; i++)
    {
        cout << genes[i] << ":\t";
        for (int j = 0; j < (int)m; j++)
        {
            if (isnan(B[i][j].val))
                cout << "NaN";
            else
                cout << B[i][j].val<<"\t";
        }
        cout << endl;
    }
    cout<<"noise matrix:"<<endl;
    for(int i=0;i<n;i++)
    {
        cout << genes[i] << ":\t";
        for(int j=0;j<m;j++)
        cout<<NoiseMatrix[i][j]<<"\t";
        cout<<endl;
    }
    cout<<endl;*/
    if (m > 2500)
        po->CLUSTER_WIDTH = max((size_t)(m * po->QUANTILE * 2) / 20 + 3, po->CLUSTER_WIDTH);
    else
        po->CLUSTER_WIDTH = max(m / 35 + 3, po->CLUSTER_WIDTH);

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < m; ++j)
        {
            if (is_missing(NoiseMatrix[i][j]))
                NoiseMatrix[i][j] = 0.0;
        }
    }
    // cout << po->CLUSTER_WIDTH << endl;
}

void read_noise_matrix(const char *filename)
{
    ifstream fin(filename);
    if (!fin)
        errAbort("cannot open noise matrix file");

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < m; ++j)
        {
            string tok;
            if (!(fin >> tok))
                errAbort("noise matrix size does not match input matrix n*m");

            if (tok == "NaN" || tok == "nan" || tok == "NAN")
                NoiseMatrix[i][j] = MISSING_VALUE;
            else
                NoiseMatrix[i][j] = stod(tok);
        }
    }
    
    double extra;
    if (fin >> extra)
    {
        err("warning: noise matrix file contains extra values after n*m entries");
    }
}
void read_file(const char *filename)
{
    ifstream input(filename);
    if (!input.is_open())
    {
        cerr << "Cannot open input file: " << filename << endl;
        exit(1);
    }

    string line, tmp;

    
    string header_line;
    while (getline(input, header_line))
    {
        if (header_line.find_first_not_of(" \t\r\n") != string::npos)
            break;
    }
    if (header_line.empty())
    {
        cerr << "Input file is empty." << endl;
        exit(1);
    }

    
    string first_data_line;
    while (getline(input, first_data_line))
    {
        if (first_data_line.find_first_not_of(" \t\r\n") != string::npos)
            break;
    }
    if (first_data_line.empty())
    {
        cerr << "No data row found in input file." << endl;
        exit(1);
    }

    auto split_tokens = [](const string &s) -> vector<string>
    {
        vector<string> tokens;
        istringstream iss(s);
        string tok;
        while (iss >> tok)
            tokens.emplace_back(tok);
        return tokens;
    };

    vector<string> header_tokens = split_tokens(header_line);
    vector<string> first_tokens = split_tokens(first_data_line);

    if (header_tokens.empty())
    {
        cerr << "Header line is invalid." << endl;
        exit(1);
    }
    if (first_tokens.size() < 2)
    {
        cerr << "First data row is invalid." << endl;
        exit(1);
    }

    
    size_t header_offset = 0;

    if (first_tokens.size() == header_tokens.size() + 1)
    {
        m = header_tokens.size();
        header_offset = 0;
    }
    else if (first_tokens.size() == header_tokens.size())
    {
        if (header_tokens.size() < 2)
        {
            cerr << "Header format is invalid." << endl;
            exit(1);
        }
        m = header_tokens.size() - 1;
        header_offset = 1;
    }
    else
    {
        cerr << "Header/data format mismatch." << endl;
        cerr << "Header token count = " << header_tokens.size()
             << ", first data row token count = " << first_tokens.size() << endl;
        exit(1);
    }

    conds.clear();
    conds.reserve(m);
    for (size_t j = 0; j < m; ++j)
        conds.emplace_back(header_tokens[j + header_offset]);

    A.clear();
    genes.clear();
    n = 0;

    auto parse_data_row = [&](const vector<string> &tokens, size_t line_no)
    {
        if (tokens.size() != m + 1)
        {
            cerr << "Invalid data row at line " << line_no
                 << ": expected " << (m + 1)
                 << " tokens, got " << tokens.size() << endl;
            exit(1);
        }

        genes.emplace_back(tokens[0]);

        vector<Array> temp;
        temp.reserve(m);

        for (size_t j = 0; j < m; ++j)
        {
            Array s;
            s.con = j;
            try
            {
                s.val = stod(tokens[j + 1]);
            }
            catch (const exception &)
            {
                cerr << "Invalid numeric value at line " << line_no
                     << ", column " << (j + 2)
                     << ": " << tokens[j + 1] << endl;
                exit(1);
            }
            s.old_val = s.val;
            temp.emplace_back(s);
        }

        A.emplace_back(std::move(temp));
        ++n;
    };

    parse_data_row(first_tokens, 2);

    size_t line_no = 2;
    while (getline(input, line))
    {
        ++line_no;
        if (line.find_first_not_of(" \t\r\n") == string::npos)
            continue;

        vector<string> tokens = split_tokens(line);
        parse_data_row(tokens, line_no);
    }

    B = vector<vector<Array>>(n, vector<Array>(m));
    NoiseMatrix.assign(n, vector<double>(m, MISSING_VALUE));
    /*cout << "col_num: " << m << "  row_num: " << n << endl;
    for(int i=0;i<m;i++)
    {
        cout<<conds[i]<<" ";
    }
    cout<<endl;
    for(int i=0;i<n;i++)
    {
        cout<<genes[i]<<": ";
        for(int j=0;j<m;j++)
        {
            cout<<A[i][j].val<<" ";
        }
        cout<<endl;
    }
    cout<<endl;*/
}