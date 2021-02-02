#ifndef Amino_Acid_Table_H
#define Amino_Acid_Table_H
#include<string>
#include<vector>
using namespace std;

class Amino_Acid_Table
{
private:
	
public:
	double Codon_Frequency[20][6];
	char AAOrder[20] = { 'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V' };
	vector<vector<string>>Input_option = {
	{"GCT","GCC","GCA","GCG"},
	{"CGT","CGC","CGA","CGG","AGA","AGG"},
	{"AAT","AAC"," "," "},
	{"GAT","GAC"," "," "},
	{"TGT","TGC"," "," "},
	{" "," ","CAA","CAG"},
	{" "," ","GAA","GAG"},
	{"GGT","GGC","GGA","GGG"},
	{"CAT","CAC"," "," "},
	{"ATT","ATC","ATA"," "},
	{"CTT","CTC","CTA","CTG","TTA","TTG"},
	{" "," ","AAG","AAA"},
	{" "," "," ","ATG"},
	{"TTC","TTT"," "," "},
	{"CCT","CCA","CCG","CCC"},
	{"TCT","TCC","AGC","TCA","TCG","AGT"},
	{"ACT","ACC","ACA","ACG"},
	{" "," "," ","TGG"},
	{"TAT","TAC"," "," "},
	{"GTT","GTC","GTA","GTG"} };
	//only keep those option above codeFrequency_threshold.

	void Read_Codon_Table(const char* filename);
	Amino_Acid_Table(char* filename);
	~Amino_Acid_Table();
	string Final_option[20][6] = { {""} };
	//remove those codon options below codeFrequency_threshold
	void scan_input_option(double codeFrequency_threshold);
};
#endif

