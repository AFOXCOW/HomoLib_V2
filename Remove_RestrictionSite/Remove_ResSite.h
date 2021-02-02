#ifndef Remove_ResSite
#define Remove_ResSite
#include<regex>
#include"../Amino_Acid_Table/Amino_Acid_Table.h"
#include<vector>
#include<string>
#include<iostream>
using namespace std;

string SiteNumberParse(string RawSite);
void FindResIndex(string res, string seq, vector<int>& index, vector<int>& length);
bool NoResSite(string CheckString, vector<string> res);

static string temp_tcga;
static string best_tcga;
void Iteration(string Protein, string ChangeString, string CheckString, int iter_level, Amino_Acid_Table table, vector<string> res, int start_Amino_Acid_inCheckStr,int end_Amino_Acid_inCheckStr,vector<vector<string>>Options);
string Choose_BestOne(string CheckString, int start_Amino_Acid_inCheckStr, int end_Amino_Acid_inCheckStr, Amino_Acid_Table table, vector<string> res,int tol);
int ReplaceResSites(vector<string> res, string& seq, Amino_Acid_Table table,int tol);
#endif
