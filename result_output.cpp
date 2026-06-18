#include "structure.h"
using namespace std;
bool CMP_node(const Node &a, const Node &b)
{
    if (a.area == b.area)
    {
        size_t x = a.in.size() + a.de.size(), y = b.in.size() + b.de.size();
        return x > y;
    }
    return a.area > b.area;
}
PIN OverlapCal(size_t x1, size_t x2)
{
    PIN p;
    size_t de_r1 = BiCluster[x1].de.size(), de_r2 = BiCluster[x2].de.size();
    size_t c1 = BiCluster[x1].index_column.size(), c2 = BiCluster[x2].index_column.size();
    size_t in_r1 = BiCluster[x1].in.size(), in_r2 = BiCluster[x2].in.size();
    vector<uint8_t> flag1(n, 0), flag2(m, 0);
    p.second.de.reserve(de_r2), p.second.in.reserve(in_r2);
    p.second.index_column.reserve(c2);

    for (size_t i = 0; i < in_r1; i++)
        flag1[BiCluster[x1].in.at(i)] = 1;
    for (size_t i = 0; i < de_r1; i++)
        flag1[BiCluster[x1].de.at(i)] = 1;
    for (size_t i = 0; i < c1; i++)
        flag2[BiCluster[x1].index_column.at(i)] = 1;
    size_t r = 0, c = 0;
    for (size_t i = 0; i < de_r2; i++)
        if (flag1[BiCluster[x2].de.at(i)])
            r++;
        else
            p.second.de.emplace_back(BiCluster[x2].de.at(i));

    for (size_t i = 0; i < in_r2; i++)
        if (flag1[BiCluster[x2].in.at(i)])
            r++;
        else
            p.second.in.emplace_back(BiCluster[x2].in.at(i));
    for (size_t i = 0; i < c2; i++)
        if (flag2[BiCluster[x2].index_column.at(i)])
            c++;
        else
            p.second.index_column.emplace_back(BiCluster[x2].index_column.at(i));

    p.second.area = (p.second.de.size() + p.second.in.size()) * p.second.index_column.size();

    p.first = {r, c};
    return p;
}
static inline size_t calc_area(const Node &x)
{
    return (x.in.size() + x.de.size()) * x.index_column.size();
}

void result_output(char *_out)
{
    ofstream output(_out);
    if (!output)
    {
        errAbort("cannot open output file");
    }

    size_t cnt = 0;

    for (size_t i = 0; i < BiCluster_num; ++i)
        BiCluster[i].area = calc_area(BiCluster[i]);

    for (size_t i = 0; i < BiCluster_num; ++i)
    {
        for (size_t j = i + 1; j < BiCluster_num; ++j)
        {
            size_t area1 = BiCluster[i].area;
            size_t area2 = BiCluster[j].area;
            if (area1 == 0 || area2 == 0)
                continue;

            size_t min_id;
            PIN p;
            if (area1 >= area2)
            {
                p = OverlapCal(i, j);
                min_id = j;
            }
            else
            {
                p = OverlapCal(j, i);
                min_id = i;
            }

            size_t rows1 = BiCluster[i].in.size() + BiCluster[i].de.size();
            size_t rows2 = BiCluster[j].in.size() + BiCluster[j].de.size();
            size_t cols1 = BiCluster[i].index_column.size();
            size_t cols2 = BiCluster[j].index_column.size();
            size_t row_denom = min(rows1, rows2);
            size_t col_denom = min(cols1, cols2);
            if (row_denom == 0 || col_denom == 0)
                continue;

            double row_overlap = p.first.first * 1.0 / row_denom;
            double col_overlap = p.first.second * 1.0 / col_denom;
            bool should_filter = (po->FILTER < eps) ? (row_overlap > po->FILTER && col_overlap > po->FILTER)
                                                    : (row_overlap >= po->FILTER && col_overlap >= po->FILTER);

            if (should_filter)
            {
                BiCluster[min_id] = std::move(p.second);
                BiCluster[min_id].area = calc_area(BiCluster[min_id]);
            }
        }
    }

    sort(BiCluster.begin(), BiCluster.begin() + BiCluster_num, CMP_node);

    vector<size_t> valid_ids;
    valid_ids.reserve(BiCluster_num);

    size_t row_len, col_len;
    for (size_t i = 0; i < BiCluster_num; ++i)
    {
        size_t rows = BiCluster[i].in.size() + BiCluster[i].de.size();
        size_t cols = BiCluster[i].index_column.size();

        if (po->MIN_LENGTH > 0)
        {
            row_len = ceil(n * po->MIN_LENGTH);
            col_len = ceil(m * po->MIN_LENGTH);
        }
        else
        {
            row_len = col_len = 0;
        }

        row_len = max(po->ROW_WIDTH, row_len);
        col_len = max(po->COL_WIDTH, col_len);

        if (rows < row_len || cols < col_len ||
            (po->IS_SINGLE_CELL_DATA && (rows < po->CLUSTER_WIDTH || cols < po->CLUSTER_WIDTH)))
            continue;

        valid_ids.push_back(i);
    }

    cnt = min(po->BLOCK_NUM, valid_ids.size());
    //cnt=10;
    output << "#Parameters: "
           << "q:" << po->QUANTILE
           << " a:" << po->ABSOLUTE_QUANTILE
           << " N:" << (int)po->NORMALIZATION
           << " l:" << po->LCS_TOLERANCE
           << " n:" << po->CLUSTER_SIZE
           << " d:" << po->DISCRETIZATION
           << " M:" << po->SEED_SELECT_MODE
           << " k:" << po->SEED_NUM_MULTIPLIER
           << " u:" << po->SEED_POOL_MULTIPLIER
           << " t:" << po->THREADS_NUM
           << " e:" << po->EXPAND_TOLERANCE
           << " b:" << po->DICHOTOMY_TOLERANCE
           << " c:" << po->COL_WIDTH
           << " r:" << po->ROW_WIDTH
           << " m:" << po->MIN_LENGTH
           << " f:" << po->FILTER
           << " S:" << po->IS_SINGLE_CELL_DATA
           << " o:" << po->BLOCK_NUM << "\n";

    output << "BiCluster_Num:" << cnt << "\n";

    for (size_t t = 0; t < cnt; ++t)
    {
        size_t i = valid_ids[t];

        auto in = BiCluster[i].in;
        auto de = BiCluster[i].de;
        auto cols = BiCluster[i].index_column;

        sort(in.begin(), in.end());
        sort(de.begin(), de.end());
        sort(cols.begin(), cols.end());

        output << "BC: " << t << "\n";

        output << "PC_Genes [" << in.size() << "]: ";
        for (size_t x : in)
            output << genes.at(x) << " ";
        output << "\n";

        output << "NC_Genes [" << de.size() << "]: ";
        for (size_t x : de)
            output << genes.at(x) << " ";
        output << "\n";

        output << "Conds [" << cols.size() << "]: ";
        for (size_t c : cols)
            output << conds.at(c) << " ";
        output << "\n\n";
    }

    uglyTime("%zu bicluster are written to %s", cnt, _out);
    cout << "All jobs completed" << endl;
}
