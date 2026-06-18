#include "structure.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include <functional>

using namespace std;

namespace
{
	static vector<double> collect_row_values(size_t row, const vector<size_t> &cols)
	{
		vector<double> vals;
		vals.reserve(cols.size());
		for (size_t c : cols)
		{
			double x = B[row][c].val;
			if (!isnan(x))
				vals.push_back(x);
		}
		return vals;
	}

	static pair<size_t, size_t> insertion_range_increasing(
		const vector<double> &sorted_vals,
		double val_l,
		double val_r)
	{
		auto it_low = lower_bound(sorted_vals.begin(), sorted_vals.end(), val_l);
		auto it_up = upper_bound(sorted_vals.begin(), sorted_vals.end(), val_r);

		return {
			static_cast<size_t>(it_low - sorted_vals.begin()),
			static_cast<size_t>(it_up - sorted_vals.begin())};
	}

	static pair<size_t, size_t> insertion_range_decreasing(
		const vector<double> &sorted_vals,
		double val_l,
		double val_r)
	{
		auto comp = greater<double>();

		auto it_low = lower_bound(sorted_vals.begin(), sorted_vals.end(), val_r, comp);
		auto it_up = upper_bound(sorted_vals.begin(), sorted_vals.end(), val_l, comp);

		return {
			static_cast<size_t>(it_low - sorted_vals.begin()),
			static_cast<size_t>(it_up - sorted_vals.begin())};
	}

	static void vote_one_row(
		size_t row,
		size_t candidate_col,
		const vector<size_t> &bic_cols,
		bool increasing, 
		vector<size_t> &sum,
		size_t &M)
	{
		const double x = B[row][candidate_col].val;
		if (isnan(x))
			return; 
		if (bic_cols.empty())
			return;

		vector<double> vals = collect_row_values(row, bic_cols);
		if (vals.empty())
			return;
			
		if (increasing)
		{
			sort(vals.begin(), vals.end());
		}
		else
		{
			sort(vals.begin(), vals.end(), greater<double>());
		}

		double change = (NoiseMatrix[row][candidate_col]) * po->DICHOTOMY_TOLERANCE;
		/*cout<<"row: "<<row<<" ";
		for(int i=0;i<vals.size();i++)
			cout<<vals[i]<<" ";
			cout<<endl;*/
		
		double a = x * (1.0 - change);
		double b = x * (1.0 + change);
		double val_l = min(a, b);
		double val_r = max(a, b);

		size_t low = 0, up = 0;
		if (increasing)
		{
			tie(low, up) = insertion_range_increasing(vals, val_l, val_r);
		}
		else
		{
			tie(low, up) = insertion_range_decreasing(vals, val_l, val_r);
		}

		const size_t max_slot = vals.size(); 
		if (low > max_slot)
			low = max_slot;
		if (up > max_slot)
			up = max_slot;

		/*cout << (increasing ? "in row: " : "de row: ")
				  << row
				  << " x=" << x
				  << " noise=" <<NoiseMatrix[row][candidate_col]
				  << " change=" << change
				  << " range=[" << val_l << "," << val_r << "]"
				  << " slot=[" << low << "," << up << "]\n";*/
		
		if (low > up)
			return;
			
		for (size_t j = low; j <= up; ++j)
		{
			++sum[j];
			if (sum[j] > M)
				M = sum[j];
		}
	}

} // namespace

size_t seed_cal(const Node &S, size_t col)
{
	const size_t len = S.index_column.size();
	vector<size_t> sum(len + 1, 0); 
	size_t M = 0;
	
	/*cout << "candidate column:"<<col<<" "<<S.index_column.size()<<endl;
	for (size_t i = 0; i < S.de.size(); ++i)
	{
		cout<<genes[S.de[i]]<<"\t";
		for(int j=0;j<S.index_column.size();++j)
		cout << B[S.de[i]][S.index_column[j]].con << "-" << B[S.de[i]][S.index_column[j]].val << "\t";
		cout<<endl;
	}
	for (size_t i = 0; i < S.in.size(); ++i)
	{
		cout<<genes[S.in[i]]<<"\t";
		for(int j=0;j<S.index_column.size();++j)
		cout << B[S.in[i]][S.index_column[j]].con << "-" << B[S.in[i]][S.index_column[j]].val << "\t";
		cout<<endl;
	}
	cout << '\n';*/

	for (size_t i = 0; i < S.de.size(); ++i)
	{
		/*cout<<"de row "<<genes[S.de[i]]<<": ";
		for(int j=0;j<S.index_column.size();j++)
		{
			cout<<B[S.de[i]][S.index_column[j]].con<<"-"<<B[S.de[i]][S.index_column[j]].val<<"\t";
		}
		cout<<endl;*/
		vote_one_row(S.de[i], col, S.index_column, false, sum, M);
	}

	for (size_t i = 0; i < S.in.size(); ++i)
	{
		/*cout<<"in row "<<genes[S.in[i]]<<": ";
		for(int j=0;j<S.index_column.size();j++)
		{
			cout<<B[S.in[i]][S.index_column[j]].con<<"-"<<B[S.in[i]][S.index_column[j]].val<<"\t";
		}
		cout<<endl;*/
		vote_one_row(S.in[i], col, S.index_column, true, sum, M);
	}

	/*cout << "max sum: " << M << '\n';
	for (size_t i = 0; i < sum.size(); ++i)
	{
		cout << "slot " << i << " count: " << sum[i] << '\n';
	}*/

	return M;
}