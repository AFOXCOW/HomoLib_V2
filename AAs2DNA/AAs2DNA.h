#ifndef AA2DNA
#define AA2DNA
#include <vector>
#include <string>
#include <ctime>
#include "../Amino_Acid_Table/Amino_Acid_Table.h"
#include "../Remove_RestrictionSite/Remove_ResSite.h"
using namespace std;

//Iteration to find all Permutation and combination and score.
void Iteration_Score(Amino_Acid_Table table, vector<int>indexs, int iteration);

//choose the best score codon for every site.
vector<string> reverseProteinToDNA(vector<string> Aligned_ProteinSequences, Amino_Acid_Table table);
#endif
