#include "AAs2DNA.h"
int seq_index;
vector<int>number;
vector<int>index;
vector<string>temp_3char;
vector<string>best_3char;
double best_score = 0;
int times;
int count_atgc_1[4] = { 0 };
int count_atgc_2[4] = { 0 };
int count_atgc_3[4] = { 0 };
void Iteration_Score(Amino_Acid_Table table, vector<int>indexs, int iteration) {
	if (iteration == 0) {
		//calculate temp score
		for (int i = 0; i < number.size(); i++) {
			//site 1
			times = number[i];
			switch (temp_3char[i][0]) {
			case 'A':
				count_atgc_1[0] = count_atgc_1[0] + times;
				break;
			case 'T':
				count_atgc_1[1] = count_atgc_1[1] + times;
				break;
			case 'G':
				count_atgc_1[2] = count_atgc_1[2] + times;
				break;
			case 'C':
				count_atgc_1[3] = count_atgc_1[3] + times;
				break;
			}
			//site 2
			switch (temp_3char[i][1]) {
			case 'A':
				count_atgc_2[0] = count_atgc_2[0] + times;
				break;
			case 'T':
				count_atgc_2[1] = count_atgc_2[1] + times;
				break;
			case 'G':
				count_atgc_2[2] = count_atgc_2[2] + times;
				break;
			case 'C':
				count_atgc_2[3] = count_atgc_2[3] + times;
				break;
			}
			//site 3
			switch (temp_3char[i][2]) {
			case 'A':
				count_atgc_3[0] = count_atgc_3[0] + times;
				break;
			case 'T':
				count_atgc_3[1] = count_atgc_3[1] + times;
				break;
			case 'G':
				count_atgc_3[2] = count_atgc_3[2] + times;
				break;
			case 'C':
				count_atgc_3[3] = count_atgc_3[3] + times;
				break;
			}
		}
		int total = count_atgc_1[0] + count_atgc_1[1] + count_atgc_1[2] + count_atgc_1[3];
		double max = 0;
		double score = 0;
		for (int i = 0; i < 4; i++) {
			if (count_atgc_1[i] > max) {
				max = count_atgc_1[i];
			}
		}
		score = score + (max / total);
		max = 0;
		for (int i = 0; i < 4; i++) {
			if (count_atgc_2[i] > max) {
				max = count_atgc_2[i];
			}
		}
		score = score + (max / total);
		max = 0;
		for (int i = 0; i < 4; i++) {
			if (count_atgc_3[i] > max) {
				max = count_atgc_3[i];
			}
		}
		score = score + (max / total);
		if (score > best_score) {
			best_score = score;
			best_3char = temp_3char;
		}
		for (int i = 0; i < 4; i++) {
			count_atgc_1[i] = 0;
			count_atgc_2[i] = 0;
			count_atgc_3[i] = 0;
		}
	}
	else {
		vector<int>new_indexs;
		int length = 0;
		for (int i = 0; i < 6; i++) {
			if (table.Final_option[indexs[0]][i].compare("") != 0) {
				length++;
			}
		}
		for (int i = 1; i < indexs.size(); i++) {
			new_indexs.push_back(indexs[i]);
		}
		for (int i = 0; i < length; i++) {
			temp_3char.push_back(table.Final_option[indexs[0]][i]);
			Iteration_Score(table, new_indexs, iteration - 1);
			temp_3char.pop_back();
		}
	}
}
vector<string> reverseProteinToDNA(vector<string> Aligned_ProteinSequences, Amino_Acid_Table table)
{
	vector<string> DNA_sequences(Aligned_ProteinSequences.size());
	for (int i = 0; i < Aligned_ProteinSequences[0].size(); i++)
	{
		//{'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'}
		int number_amino_acid[21] = { 0 };
		for (int j = 0; j < Aligned_ProteinSequences.size(); j++) {
			switch (Aligned_ProteinSequences[j][i])
			{
			case 'A':
				number_amino_acid[0]++;
				break;
			case 'R':
				number_amino_acid[1]++;
				break;
			case 'N':
				number_amino_acid[2]++;
				break;
			case 'D':
				number_amino_acid[3]++;
				break;
			case 'C':
				number_amino_acid[4]++;
				break;
			case 'Q':
				number_amino_acid[5]++;
				break;
			case 'E':
				number_amino_acid[6]++;
				break;
			case 'G':
				number_amino_acid[7]++;
				break;
			case 'H':
				number_amino_acid[8]++;
				break;
			case 'I':
				number_amino_acid[9]++;
				break;
			case 'L':
				number_amino_acid[10]++;
				break;
			case 'K':
				number_amino_acid[11]++;
				break;
			case 'M':
				number_amino_acid[12]++;
				break;
			case 'F':
				number_amino_acid[13]++;
				break;
			case 'P':
				number_amino_acid[14]++;
				break;
			case 'S':
				number_amino_acid[15]++;
				break;
			case 'T':
				number_amino_acid[16]++;
				break;
			case 'W':
				number_amino_acid[17]++;
				break;
			case 'Y':
				number_amino_acid[18]++;
				break;
			case 'V':
				number_amino_acid[19]++;
				break;
			case '-':
				number_amino_acid[20]++;
				break;
			}
		}
		for (int i = 0; i < 20; i++) {
			if (number_amino_acid[i] != 0) {
				index.push_back(i);
				number.push_back(number_amino_acid[i]);
			}
		}
		if (index.size() > 1) {
			Iteration_Score(table, index, index.size());
		}
		else {
			int number_choice = 0;
			for (int i = 0; i < 6; ++i) {
				if (table.Final_option[index[0]][i] != "") {
					number_choice++;
				}
			}
			srand((unsigned)time(NULL));
			int n = rand() % (number_choice);
			best_3char.push_back(table.Final_option[index[0]][n]);
		}

		for (int j = 0; j < Aligned_ProteinSequences.size(); j++)
		{
			switch (Aligned_ProteinSequences[j][i])
			{
			case 'A':
				for (int i = 0; i < index.size(); i++) {
					if (index[i] == 0) {
						seq_index = i;
						break;
					}
				}
				DNA_sequences[j].append(best_3char[seq_index]);
				break;
			case 'R':
				for (int i = 0; i < index.size(); i++) {
					if (index[i] == 1) {
						seq_index = i;
						break;
					}
				}
				DNA_sequences[j].append(best_3char[seq_index]);
				break;
			case 'N':
				for (int i = 0; i < index.size(); i++) {
					if (index[i] == 2) {
						seq_index = i;
						break;
					}
				}
				DNA_sequences[j].append(best_3char[seq_index]);
				break;
			case 'D':
				for (int i = 0; i < index.size(); i++) {
					if (index[i] == 3) {
						seq_index = i;
						break;
					}
				}
				DNA_sequences[j].append(best_3char[seq_index]);
				break;
			case 'C':
				for (int i = 0; i < index.size(); i++) {
					if (index[i] == 4) {
						seq_index = i;
						break;
					}
				}
				DNA_sequences[j].append(best_3char[seq_index]);
				break;
			case 'Q':
				for (int i = 0; i < index.size(); i++) {
					if (index[i] == 5) {
						seq_index = i;
						break;
					}
				}
				DNA_sequences[j].append(best_3char[seq_index]);
				break;
			case 'E':
				for (int i = 0; i < index.size(); i++) {
					if (index[i] == 6) {
						seq_index = i;
						break;
					}
				}
				DNA_sequences[j].append(best_3char[seq_index]);
				break;
			case 'G':
				for (int i = 0; i < index.size(); i++) {
					if (index[i] == 7) {
						seq_index = i;
						break;
					}
				}
				DNA_sequences[j].append(best_3char[seq_index]);
				break;
			case 'H':
				for (int i = 0; i < index.size(); i++) {
					if (index[i] == 8) {
						seq_index = i;
						break;
					}
				}
				DNA_sequences[j].append(best_3char[seq_index]);
				break;
			case 'I':
				for (int i = 0; i < index.size(); i++) {
					if (index[i] == 9) {
						seq_index = i;
						break;
					}
				}
				DNA_sequences[j].append(best_3char[seq_index]);
				break;
			case 'L':
				for (int i = 0; i < index.size(); i++) {
					if (index[i] == 10) {
						seq_index = i;
						break;
					}
				}
				DNA_sequences[j].append(best_3char[seq_index]);
				break;
			case 'K':
				for (int i = 0; i < index.size(); i++) {
					if (index[i] == 11) {
						seq_index = i;
						break;
					}
				}
				DNA_sequences[j].append(best_3char[seq_index]);
				break;
			case 'M':
				for (int i = 0; i < index.size(); i++) {
					if (index[i] == 12) {
						seq_index = i;
						break;
					}
				}
				DNA_sequences[j].append(best_3char[seq_index]);
				break;
			case 'F':
				for (int i = 0; i < index.size(); i++) {
					if (index[i] == 13) {
						seq_index = i;
						break;
					}
				}
				DNA_sequences[j].append(best_3char[seq_index]);
				break;
			case 'P':
				for (int i = 0; i < index.size(); i++) {
					if (index[i] == 14) {
						seq_index = i;
						break;
					}
				}
				DNA_sequences[j].append(best_3char[seq_index]);
				break;
			case 'S':
				for (int i = 0; i < index.size(); i++) {
					if (index[i] == 15) {
						seq_index = i;
						break;
					}
				}
				DNA_sequences[j].append(best_3char[seq_index]);
				break;
			case 'T':
				for (int i = 0; i < index.size(); i++) {
					if (index[i] == 16) {
						seq_index = i;
						break;
					}
				}
				DNA_sequences[j].append(best_3char[seq_index]);
				break;
			case 'W':
				for (int i = 0; i < index.size(); i++) {
					if (index[i] == 17) {
						seq_index = i;
						break;
					}
				}
				DNA_sequences[j].append(best_3char[seq_index]);
				break;
			case 'Y':
				for (int i = 0; i < index.size(); i++) {
					if (index[i] == 18) {
						seq_index = i;
						break;
					}
				}
				DNA_sequences[j].append(best_3char[seq_index]);
				break;
			case 'V':
				for (int i = 0; i < index.size(); i++) {
					if (index[i] == 19) {
						seq_index = i;
						break;
					}
				}
				DNA_sequences[j].append(best_3char[seq_index]);
				break;
			case '-':
				break;
			}
		}
		//clear all pre infomation
		index.clear();
		number.clear();
		temp_3char.clear();
		best_3char.clear();
		best_score = 0;
		times = 0;
		for (int i = 0; i < 4; i++) {
			count_atgc_1[i] = 0;
			count_atgc_2[i] = 0;
			count_atgc_3[i] = 0;
		}
	}
	return DNA_sequences;
}