#ifndef UTILS_H
#define UTILS_H
#include<string>
#include<vector>
#include<fstream>
#include<iostream>
using namespace std;

typedef struct seq {
	vector<string> sequences;
	vector<string> seq_names;
}FASTA_Seq;
void erase_all_space(string& s);
void SplitString(const string& s, vector<string>& tokens, const string& delimiters);
vector<vector<double>> Read_Similar_Table(string filepath, string final_Element_List);
FASTA_Seq Read_FASTA(string filepath);
string replace(const string& str, const string& src, const string& dest);
string complementary(string source);
int Find_length_index_forward_all(vector<string> Aligned_Sequences, int start, int length);
int Find_length_index_forward_one(vector<string> Aligned_Sequences, int start, int length);
int Find_length_index_backward_all(vector<string> Aligned_Sequences, int start, int length);
int Find_length_index_backward_one(vector<string> Aligned_Sequences, int start, int length);
#endif