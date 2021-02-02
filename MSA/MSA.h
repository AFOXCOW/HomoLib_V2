#ifndef HOMOLIB_MSA_H
#define HOMOLIB_MSA_H

#include<vector>
#include<string>
#include<set>
#include<fstream>
using namespace std;

// Label is a flag. Sequences with same label are aligned already.
typedef struct Sequence_{
	string sequence;
	unsigned int label=0;
} Sequence_labeled;

class MSA
{
private:

public:
	const int sequenceMaxNumber = 100;

	// For Protein sequence, Type of Amino Acid(20 sequence element). For DNA sequence, TCGA(4 seq element)
	int Number_of_Seq_Element = 0;
	int Number_of_Seq_Element_Plus_Spaceholder = Number_of_Seq_Element + 1;//one more letter(placeholder '-')

	// store the index.
	// For example, the input sequence element order is "ARND....."
	// so the for amino acid 'A',alphaBetHashTable[int('A')] = A.index in sequence element(0).
	int alphaBetHashTable[128] = { 0 };
	int number_of_Sequences = 0;

	// store the most similar groups` sequence label 
	int most_similar_group_one = 0;
	int most_similar_group_two = 1;

	// if only one sequence label,means all sequences in one group, means alignment finished. and Number_of_not_aligned_group=0
	int Number_of_not_aligned_group = 0;
	int Sequence_MaxLength = 2000;
	vector<Sequence_labeled> Aligning_Sequence;

	//Update Number_of_not_aligned_group, Once it is zero, End MSA.
	//void Update_Number_of_not_aligned_group();

	//choose most_similar_group
	void ChooseAlignedPair();

	//align two group, and apply the align pattern to both group.
	int AlignAndUpdate(int firstSet, int secondSet);

	//score two group
	double ScoreAlignAlignments(unsigned int firstSet, unsigned int secondSet);
	double AlignColumnToSpace(int setLabel, int columnLocation);
	double AlignColumnToColumn(int setLabelVer, int verColumnLocation, int setLabelHor, int horColumnLocation);
	double MaxofThree(double A, double B, double C);

	// It scores the similar between 2 elements.eg: 
	// for blosum62,A matches A will get score 4.
	// A matches spaceholder'-' will get score -4. 
	// The higher, the more similar
	double Similar_Table[21][21] = { 0 };
	MSA(vector<string> Input_Sequence, vector<vector<double>> Input_Similar_Table, string Seq_Element);
	~MSA(void);
	int IterationAlignProcess();
	vector<string>output();
};
#endif