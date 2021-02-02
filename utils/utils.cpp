#include "utils.h"
#include <string>
using namespace std;

void erase_all_space(string& s) {
	int index = 0;
	if (!s.empty()) {
		while ((index = s.find(' ', index)) != string::npos) {
			s.erase(index, 1);
		}
	}
}
void SplitString(const string& s, vector<string>& tokens, const string& delimiters)
{
	string::size_type lastPos = s.find_first_not_of(delimiters, 0);
	string::size_type pos = s.find_first_of(delimiters, lastPos);
	while (string::npos != pos || string::npos != lastPos) {
		tokens.push_back(s.substr(lastPos, pos - lastPos));
		lastPos = s.find_first_not_of(delimiters, pos);
		pos = s.find_first_of(delimiters, lastPos);
	}
}
vector<vector<double>> Read_Similar_Table(string filepath, string final_Element_List) {
	/*
	Demo Blosum62
			C  S  T  P  A  G  N  D  E  Q  H  R  K  M  I  L  V  F  Y  W
		C   9
		S - 1  4
		T - 1  1  5
		P - 3 - 1 - 1  7
		A   0  1  0 - 1  4
		G - 3  0 - 2 - 2  0  6
		N - 3  1  0 - 2 - 2  0  6
		D - 3  0 - 1 - 1 - 2 - 1  1  6
		E - 4  0 - 1 - 1 - 1 - 2  0  2  5
		Q - 3  0 - 1 - 1 - 1 - 2  0  0  2  5
		H - 3 - 1 - 2 - 2 - 2 - 2  1 - 1  0  0  8
		R - 3 - 1 - 1 - 2 - 1 - 2  0 - 2  0  1  0  5
		K - 3  0 - 1 - 1 - 1 - 2  0 - 1  1  1 - 1  2  5
		M - 1 - 1 - 1 - 2 - 1 - 3 - 2 - 3 - 2  0 - 2 - 1 - 1  5
		I - 1 - 2 - 1 - 3 - 1 - 4 - 3 - 3 - 3 - 3 - 3 - 3 - 3  1  4
		L - 1 - 2 - 1 - 3 - 1 - 4 - 3 - 4 - 3 - 2 - 3 - 2 - 2  2  2  4
		V - 1 - 2  0 - 2  0 - 3 - 3 - 3 - 2 - 2 - 3 - 3 - 2  1  3  1  4
		F - 2 - 2 - 2 - 4 - 2 - 3 - 3 - 3 - 3 - 3 - 1 - 3 - 3  0  0  0 - 1  6
		Y - 2 - 2 - 2 - 3 - 2 - 3 - 2 - 3 - 2 - 1  2 - 2 - 2 - 1 - 1 - 1 - 1  3  7
		W - 2 - 3 - 2 - 4 - 3 - 2 - 4 - 4 - 3 - 2 - 2 - 3 - 3 - 1 - 3 - 2 - 3  1  2  11
		- - 4 - 4 - 4 - 4 - 4 - 4 - 4 - 4 - 4 - 4 - 4 - 4 - 4 - 4 - 4 - 4 - 4 - 4 - 4 -4 -4 -4 -4 0
	*/
	vector<vector<double>>similar_table;
	int Number_of_row = final_Element_List.size();
	int Number_of_column = final_Element_List.size();
	int Line_index = 0;
	ifstream inFile(filepath, iostream::in);
	string inLine;
	char row_current_Protein;
	char column_current_Protein;
	string input_Protein_List;
	int row_index_in_final_Protein_List;
	int column_index_in_final_Protein_List;
	vector<string> split_vector;
	similar_table.resize(Number_of_row);
	for (int i = 0; i < Number_of_row; i++) {
		similar_table[i].resize(Number_of_column);
	}
	if (!inFile) {
		cout << "Open Similar Score Table File Error!" << endl;
		exit(-1);
	}
	else {
		while (getline(inFile, inLine)) {
			if (Line_index == 0) {
				erase_all_space(inLine);
				input_Protein_List = inLine;
				Line_index++;
			}
			else {
				row_current_Protein = inLine[0];
				row_index_in_final_Protein_List = final_Element_List.find(row_current_Protein);
				SplitString(inLine, split_vector, " ");
				for (int i = 1; i < split_vector.size(); i++) {
					column_current_Protein = input_Protein_List[i];
					column_index_in_final_Protein_List = final_Element_List.find(column_current_Protein);
					similar_table[row_index_in_final_Protein_List][column_index_in_final_Protein_List] = stoi(split_vector[i]);
				}
				split_vector.clear();
				Line_index++;
			}
		}
	}
	return similar_table;
}
FASTA_Seq Read_FASTA(string filepath) {
	ifstream inFile(filepath, iostream::in);
	FASTA_Seq Sequences;
	string inLine;
	if (!inFile) {
		cout << "Open File Error!" << endl;
		exit(-1);
	}
	else {
		while (getline(inFile, inLine)) {
			if (inLine.find('>') == 0) {
				Sequences.seq_names.push_back(inLine);
				getline(inFile, inLine);
				Sequences.sequences.push_back(inLine);
			}
			else {
				cout << "FASTA Format Error!" << endl;
				exit(-1);
			}
		}
	}
	return Sequences;
}

string replace(const string& str, const string& src, const string& dest)
{
	string ret;
	string::size_type pos_begin = 0;
	string::size_type pos = str.find(src);
	while (pos != string::npos)
	{
		ret.append(str.data() + pos_begin, pos - pos_begin);
		ret += dest;
		pos_begin = pos + 1;
		pos = str.find(src, pos_begin);
	}
	if (pos_begin < str.length())
	{
		ret.append(str.begin() + pos_begin, str.end());
	}
	return ret;
}


string complementary(string source) {
	string complementary_string(source);
	int strLen = source.size();
	for (int i = 0; i < strLen; i++) {
		switch (source[strLen - 1 - i])
		{
		case 'A':
			complementary_string[i] = 'T'; 
			break;
		case 'a':
			complementary_string[i] = 'T'; 
			break;
		case 'T':
			complementary_string[i] = 'A'; 
			break;
		case 't':
			complementary_string[i] = 'A'; 
			break;
		case 'G':
			complementary_string[i] = 'C'; 
			break;
		case 'g':
			complementary_string[i] = 'C'; 
			break;
		case 'C':
			complementary_string[i] = 'G'; 
			break;
		case 'c':
			complementary_string[i] = 'G'; 
			break;
		default:
			complementary_string[i] = 'U'; //Unknown
		} // end switch
	} // end for
	return complementary_string;
}


//find the index which means all aligned sequences` substring [start,end) have length size(gap not include) 
int Find_length_index_forward_all(vector<string> Aligned_Sequences, int start, int length) {
	int end = start + length; //begin with start(included), end with 'end'(not included)
	if (end > Aligned_Sequences[0].size()) {
		return -1; //too short. 
	}
	bool exist_short = true;
	while (true) {
		for (int i = 0; i < Aligned_Sequences.size(); ++i) {
			string substring = Aligned_Sequences[i].substr(start, end - start);
			string rp_substring = replace(substring, "-", "");
			if (rp_substring.size() < length) {
				//exist at least one sequence, its length is shorter than requirement.
				end++;
				exist_short = true;
				break;
			}
			exist_short = false;
		}
		if (end > Aligned_Sequences[0].size()) {
			return -1; //can not find the index.
		}
		if (exist_short) {
			continue;
		}
		else {
			break;
		}
	}
	return end;
}

//find the index that all sequences length below maxlength.
int Find_length_index_forward_one(vector<string> Aligned_Sequences, int start, int length) {
	int end = start + length;
	if (end > Aligned_Sequences[0].size()) {
		return Aligned_Sequences[0].size();
	}
	bool exist_equal = false;
	while (true) {
		for (int i = 0; i < Aligned_Sequences.size(); ++i) {
			string substring = Aligned_Sequences[i].substr(start, end - start);
			string rp_substring = replace(substring, "-", "");
			if (rp_substring.size() == length) {
				//exist at least one sequence, its length is equal to requirements.
				return end;
			}
		}
		end++;
		if (end > Aligned_Sequences[0].size()) {
			return Aligned_Sequences[0].size();
		}
	}
}


//find the index which means all aligned sequences` substring [end,start] have at least length size(gap not include) 
int Find_length_index_backward_all(vector<string> Aligned_Sequences, int start, int length) {
	//start included, end included.
	int end = start - length + 1;
	if (end < 0) {
		return -1; //can not find the index.
	}
	bool exist_short = true;
	while (true) {
		for (int i = 0; i < Aligned_Sequences.size(); ++i) {
			string substring = Aligned_Sequences[i].substr(end, start - end + 1);
			string rp_substring = replace(substring, "-", "");
			if (rp_substring.size() < length) {
				//exist at least one sequence, its length is shorter than requirement.
				end--;
				exist_short = true;
				break;
			}
			exist_short = false;
		}
		if (end < 0) {
			return -1; //can not find the index.
		}
		if (exist_short) {
			continue;
		}
		else {
			break;
		}
	}
	return end;
}

//find the index that all sequences length below maxlength.
int Find_length_index_backward_one(vector<string> Aligned_Sequences, int start, int length) {
	int end = start - length + 1;
	if (end < 0) {
		return 0;
	}
	bool exist_equal = false;
	while (true) {
		for (int i = 0; i < Aligned_Sequences.size(); ++i) {
			string substring = Aligned_Sequences[i].substr(end, start - end + 1);
			string rp_substring = replace(substring, "-", "");
			if (rp_substring.size() == length) {
				return end;
			}
		}
		end--;
		if (end < 0) {
			return 0;
		}
	}
}