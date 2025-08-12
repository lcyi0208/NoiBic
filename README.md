# NoiBic
A noise-tolerant biclustering algorithm for analyzing gene expression data at various scales.

## Table of contents

- [Installation](#installation)
- [Usage](#usage)
- [Example](#example)
- [Bicluster Similarity Metrics](#bicluster-similarity-metrics)
  
## Installation:
 g++ with support for C++11 (e.g. 8.4.0)

 Boost >= 1.47.0

 ### Installing Boost
 
 ```bash
 
 a) download a version of boost from https://www.boost.org/releases/latest/ and unpack it

    $ tar zxvf boost_1_88_0.tar.gz

 b) change to the boost directory and run ./bootstrap.sh

    $ cd boost_1_88_0
    $ ./bootstrap.sh

 c) run

    $ ./b2 install --prefix=/your/boost/dir
```
### Installing NoiBic

```bash
git clone https://github.com/lcyi0208/NoiBic.git
cd NoiBic
make BOOST_INCLUDE=/your/boost/dir
```
## Usage

```bash
./noibic -i filename [argument list]
```
The specific parameter list for each step of the algorithm is as follows.

### Input

```bash
-i : The input file must be in one of two tab-delimited formats: 
     o        cond1    cond2    cond3 
     gene1      2.4      3.5     -2.4 
     gene2     -2.1      0.0      1.2
```
### Data Preprocessing

```bash
-q : Remove non-expressed data based on numerical values   
     Floating-point value in the range (0, 0.5], default: 0.06   
     (see details in the Methods section of the paper).   
-a : Remove non-expressed data by index position.   
     Binary variable (0 or 1), default: 0      
-s : Replace zero values with random numbers drawn from N(0,1).   
     Binary variable (0 or 1), default: 0.   
-z : Exclude zero values from clustering.   
     Binary variable (0 or 1), default: 0.
```
### Biclustering

```bash
-l : Permissible mismatch rate during element swapping when searching   
     for the Longest Common Subsequence (LCS) between two sequences.   
     Floating-point value in the range [0, 1], default: 0.2.   
-n : Minimum cluster width as a fraction of the original seed length.   
     Floating-point value in the range [0, 1], default: 0.12.   
-d : Threshold for filtering candidate seed sequences based on maximum redundancy.   
     Floating-point value in the range [0, 1], default: 0.9.   
-t : Number of threads for multi-threaded execution.   
     Positive integer, default: 16. 
```
### Expansion

```bash
-e : Permissible error rate during bicluster expansion.   
     Floating-point value in the range (0.5, 1], default: 0.85.   
-b : Permissible numerical error for column expansion using binary search.   
     Floating-point value in the range [0, 1], or 2 to disable column expansion.   
     Default: 0.2. 
```
### Output

```bash
-c : Minimum number of columns in a block.   
     Integer ≥ 3, default: 6.   
-r : Minimum number of rows in a block.   
     Integer ≥ 3, default: 8.   
-m : Minimum block size as a fraction of the original matrix dimensions. 
    Floating-point value in the range [0, 1], default:0.0. 
-f : Overlap filtering threshold for biclusters.   
     Floating-point value in the range [0, 1.0], default: 1 (no filtering).   
-S : Specify whether the input data is single-cell data.   
     Binary variable (0 or 1), default: 0.   
-o : Number of biclusters to report.   
     Positive integer, default: 30.   
-p : Output file name.   
     Defaults to the input file name.  
```
## Example
```bash
./noibic -i ./test/toy.txt -p output -q 0.5 -o 3
```
Running the above command produces output.blocks.
In the .blocks file, the first line contains the algorithm’s runtime parameters, the second line shows the number of output biclusters, and the subsequent lines present the detailed information for each BC.
We classify the genes into positively and negatively regulated genes based on their correlations, as shown below.
```bash
#Parameters: q:0.5 a:0 s:0 z:0 l:0.2 n:12 d:0.9 t:16 e:0.85 b:0.2 c:6 r:4 m:0 f:0.8 S:0 o:3
BiCluster_Num:3
BC: 0
PC_Genes [12]: gene40 gene42 ... 
NC_Genes [8]: gene41 gene47 ... 
Conds [20]: cond40 cond41 cond42 ...
...
```
## Bicluster Similarity Metrics
When the biclusters embedded in the input matrix are known, the clustering results are assessed using the following criteria. We computes list-wise bicluster similarity between two sets of biclusters:

- G (Expected): Ground-truth biclusters

- D (Found): Biclusters discovered by an algorithm

It implements two evaluation standards:

---

### (i). Jaccard-based Recovery & Relevance

For a pair of biclusters $\( e \in G \)$ and $\( f \in D \)$:

$$J(e, f) = \frac{|e \cap f|}{|e \cup f|}$$

Where:

- $\lvert e \cap f\rvert$ = **intersection area**=number of matrix cells that belong to both $e$ and $f$=
  $\lvert R_e \cap R_f\rvert \times \lvert C_e \cap C_f\rvert$.

- $\lvert e \cup f\rvert$ = **union area**= number of cells that belong to either $e$ or $f$.


**List-level Recovery $(S(G,D))$:**

$$\text{Recovery}(G, D) = \frac{1}{|G|} \sum_{e \in G} \max_{f \in D} J(e, f)$$

**List-level Relevance $(S(D,G))$ :**

$$\text{Relevance}(G, D) = \frac{1}{|D|} \sum_{f \in D} \max_{e \in G} J(e, f)$$

---

### (ii). Asymmetric Recovery / Relevance
These metrics measure overlap relative to only one bicluster’s area:
- **Pairwise Recovery:**

$$\text{recovery}(e, f) = \frac{|e \cap f|}{|e|}$$

Fraction of the ground-truth bicluster $e$ that is covered by predicted bicluster $f$.

- **Pairwise Relevance:**

$$\text{relevance}(e, f) = \frac{|e \cap f|}{|f|}$$

Fraction of the predicted bicluster $f$ that overlaps with ground-truth bicluster $e$.
List-level definitions follow the same “average best-match” idea:
- **Recovery $(S(G, D))$** = Average of the best recovery $(e,f)$ values for each $e$ in $G$
- **Relevance $(S(D, G))$** = Average of the best relevance $(e,f)$ values for each $f$ in $D$
---

### Input File Format

Each bicluster is described using **two lines**:
1. Row indices (space-separated integers)
2. Column indices (space-separated integers)

Biclusters are separated by a **blank line**.  

---

### Usage
Command-line example:
```bash
python bicluster_metrics.py --expected expected.txt --found predicted.txt --metric jaccard --show-matrix
```
Parameters:

- --expected / -g → File containing ground-truth biclusters

- --found / -d → File containing predicted biclusters

- --metric / -m

  - jaccard = Jaccard-based Recovery & Relevance

  - rr = Asymmetric Recovery/Relevance

- --show-matrix → If set, prints a pairwise similarity matrix (rows=expected, cols=found)

---
### Output
Example:
```bash
# Metric: jaccard Recovery (S(G,D))=0.875000 Relevance (S(D,G))=0.842105
# Pairwise score matrix (rows=expected, cols=found):
0.000000 0.800000 0.000000
1.000000 0.000000 0.000000
0.000000 0.000000 0.750000
```
- **Recovery / Relevance values:** Overall summary scores

- **Matrix:** Each entry $(i,j)$ is the chosen similarity metric between expected $[i]$ and found $[j]$ .
