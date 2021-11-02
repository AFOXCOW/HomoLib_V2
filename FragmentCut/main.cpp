#define _CRT_SECURE_NO_WARNINGS
#include<iostream>
#include<string>
#include<vector>
#include<fstream>
#include"../utils/utils.h"
#include"../MSA/MSA.h"
#include<time.h>
using namespace std;
//find the index which means all aligned sequences` substring [start,index] have length size(gap not include) 
int Find_length_index_forward(vector<string> Aligned_Sequences, int start, int length) {
	int end = start + length;
	bool exist_short = true;
	while (true) {
		for (int i = 0; i < Aligned_Sequences.size(); ++i) {
			string substring = Aligned_Sequences[i].substr(start, end - start);
			string rp_substring = replace(substring, "-", "");
			if (rp_substring.size() < length) {
				//exist at least one sequence, its length is shorter than requirement.
				end++;
				exist_short = true;
				if (end > Aligned_Sequences[0].size()) {
					return -1; //can not find the index.
				}
				break;
			}
			exist_short = false;
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
int Find_length_index_backward(vector<string> Aligned_Sequences, int start, int length) {
	int end = start - length;
	if (end < 0) {
		return -1; //can not find the index.
	}
	bool exist_short = true;
	while (true) {
		for (int i = 0; i < Aligned_Sequences.size(); ++i) {
			string substring = Aligned_Sequences[i].substr(end, start-end);
			string rp_substring = replace(substring, "-", "");
			if (rp_substring.size() < length) {
				//exist at least one sequence, its length is shorter than requirement.
				end--;
				exist_short = true;
				if (end < 0) {
					return -1; //can not find the index.
				}
				break;
			}
			exist_short = false;
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
double ScoreAlignedSite(vector<string>Aligned_Sequence, int site, vector<vector<double>>similar_table) {
	string final_tcga_list = "ATCG-";
	double score = 0;
	char first;
	char second;
	for (int i = 0; i < Aligned_Sequence.size(); ++i) {
		first = Aligned_Sequence[i][site];
		for (int j = 0; j < Aligned_Sequence.size(); ++j) {
			second = Aligned_Sequence[j][site];
			score += similar_table[final_tcga_list.find(first)][final_tcga_list.find(second)];
		}
	}
	return score;
}
vector<int> Cut_into_Fragments(vector<string>Aligned_Sequence, int frag_number, int overlap_Length, vector<vector<double>>similar_table, int min_length, int max_length)
{
	int avg_length=0;
	for (int i = 0; i < Aligned_Sequence.size(); ++i) {
		avg_length += replace(Aligned_Sequence[i], "-", "").size();
	}
	avg_length /= Aligned_Sequence.size();
	int frag_Length = Aligned_Sequence[0].size() / frag_number;
	int aligned_Length = Aligned_Sequence[0].size();
	int countA = 0, countT = 0, countC = 0, countG = 0;
	int sequenceMaxLength = 2500;
	int temp_length = 0;
	int reverse_temp_length = 0;
	vector<int>cut;//store cut site
	vector<int>temp_cutsite;
	cut.push_back(0);
	int start, end, best_cutsite;
	int score, max_score = 0;

	temp_cutsite.push_back(0);
	//at first, cut the sequence on average length.
	for (int i = 1; i < frag_number; i++) {
		int tmp_cut = Find_length_index_forward(Aligned_Sequence, 0, avg_length/ frag_number * i);
		temp_cutsite.push_back(tmp_cut);
	}

	for (int k = 1; k < temp_cutsite.size(); k++) {
		//forward and backward 75 site range. to find a better cutsite.
		int min_start = Find_length_index_forward(Aligned_Sequence, temp_cutsite[k - 1], min_length);
		start = temp_cutsite[k] - 75;
		start = start > min_start ? start : min_start;
		end = temp_cutsite[k] + 75;
		max_score = 0;
		score = 0;
		for (; start < end - overlap_Length; start++) {
			for (int j = start; j < start + overlap_Length; j++) {
				score += ScoreAlignedSite(Aligned_Sequence,j,similar_table);
			}
			if (score > max_score) {
				max_score = score;
				best_cutsite = start + overlap_Length;//not include best_cutsite
			}
			score = 0;
		}
		cut.push_back(best_cutsite);
	}
	cut.push_back(Aligned_Sequence[0].size());
	//every cut site is the end of high homogours region(not include). the start of next frag(include).
	return cut;
}

int main(int argc, char* argv[]) {
	time_t t1 = clock();

	//Pre Info
	cout << "Loading user`s Parameters" << endl;
	string Email = argv[1];
	string UserName = argv[2];
	int min_length = atoi(argv[3]);
	int max_length = atoi(argv[4]);
	string  FASTA_path= argv[5];
	int overlapvalue = atoi(argv[6]);
	string dirname = argv[7];
	string PrimerF = argv[8];
	string PrimerR = argv[9];
	int mail_or_not = atoi(argv[10]); //1 for mail. others for not mail.
	cout << "Finished" << endl;

	//Multiple Sequences Alignment
	cout << "Processing Multiple Sequences Alignment" << endl;
	FASTA_Seq Sequences;
	vector<string> Aligned_Sequences;
	vector<vector<double>>similar_table;
	string final_Protein_List = "ATCG-";
	int Number_of_row = 5;
	int Number_of_column = 5;
	similar_table.resize(Number_of_row);
	for (int i = 0; i < Number_of_row; i++) {
		similar_table[i].resize(Number_of_column);
	}

	string similar_table_filePath = "Similar_Table.txt";
	Sequences = Read_FASTA(FASTA_path);
	string final_tcga_list = "ATCG-";
	similar_table = Read_Similar_Table(similar_table_filePath, final_tcga_list);
	MSA MSA_object(Sequences.sequences, similar_table, final_Protein_List);
	MSA_object.IterationAlignProcess();
	Aligned_Sequences = MSA_object.output();
	cout << "Finished" << endl;
	time_t t2 = clock();
	cout << "MSA costs " << (double)(t2 - t1) / CLOCKS_PER_SEC << " s" << endl;


	//estimate the Fragment number.
	cout << "Estimate Fragments number." << endl;
	vector<int> cut;
	double Avg_length = 0.0;
	double sum_length = 0.0;
	for (int i = 0; i < Sequences.sequences.size(); i++) {
		sum_length += Sequences.sequences[i].size();
	}
	Avg_length = sum_length / Aligned_Sequences.size();
	double expected_frag_length = (double(min_length) + double(max_length)) / 2;
	int fragmentnumber = Avg_length / expected_frag_length + 0.5;
	cout << "Finished" << endl;

	if (fragmentnumber < 2) {
		string error_file = "./" + dirname + "/" + UserName + "_error.txt";
		fstream out_error(error_file, iostream::out);
		out_error << "The input sequences are too short to be cut into fragments" << endl;
		exit(-1);
	}

	//cut into Fragments
	cout << "Cutting DNA sequences into fragments." << endl;
	cut = Cut_into_Fragments(Aligned_Sequences, fragmentnumber, overlapvalue, similar_table,min_length,max_length);
	cout << "Finished" << endl;

	//save to file 
	cout << "Save to file." << endl;
	string file = "./"+dirname + "/" + UserName + "_FragCut_fragmentsresult.txt";
	fstream outfragfile(file, iostream::out);
	outfragfile << "HomoLibrary Output - " << UserName << endl;
	outfragfile << endl;
	time_t nowtime;
	nowtime = time(NULL);
	struct tm* local, * gm;
	local = localtime(&nowtime);
	outfragfile << "Job started on " << local->tm_mday << "/" << local->tm_mon + 1 << "/" << local->tm_year + 1900 << " at " << local->tm_hour << ":" << local->tm_min << ":" << local->tm_sec << endl << endl;
	outfragfile << "User Name: " << UserName << endl << endl;
	outfragfile << "Email: " << Email << endl << endl;
	outfragfile << "Fragment Length: " << min_length << "--" << max_length << endl << endl;
	outfragfile << "Overlap Length: " << overlapvalue << endl << endl;
	outfragfile << "Primer-F: " << PrimerF << endl << endl;
	outfragfile << "Primer-R: " << PrimerR << endl << endl;
	outfragfile << "Input DNA Sequences: " << endl;

	for (int j = 0; j < Sequences.sequences.size(); j++)
	{
		outfragfile << Sequences.seq_names[j]<< "(Length = " << Sequences.sequences[j].size() << "): " << Sequences.sequences[j];
		outfragfile << endl;
	}
	outfragfile << endl << endl;
	outfragfile << "The Result of Aligned DNA Sequences :" << endl;
	for (int j = 0; j < Aligned_Sequences.size(); j++)
	{
		outfragfile << Sequences.seq_names[j] << "(Length = " << Aligned_Sequences[j].size() << "): " << Aligned_Sequences[j];
		outfragfile << endl;
	}
	outfragfile << endl << endl;

	outfragfile << "The Result of Designed Fragments Sequences :" << endl;
	for (int point = 0; point < cut.size() - 1; point++)
	{
		outfragfile << "Fragment " << point + 1 << endl;
		string fasta_file = "./"+dirname + "/" + UserName + "_FragCut_fragment" + to_string(point + 1) + ".fasta";
		fstream fasta_out(fasta_file, iostream::out);
		for (int i = 0; i < Aligned_Sequences.size(); i++) {
			outfragfile << Sequences.seq_names[i];
			fasta_out << Sequences.seq_names[i] <<"_frag"<<point+1<< endl;
			if (point == 0) {//First Fragment
				outfragfile <<PrimerF <<replace(Aligned_Sequences[i].substr(cut[point], cut[point + 1] - cut[point]), "-", "") << PrimerR<< "  Length without primers: " << replace(Aligned_Sequences[i].substr(cut[point], cut[point + 1] - cut[point]), "-", "").length() << endl;
				if (i == Aligned_Sequences.size() - 1) {
					fasta_out << PrimerF <<replace(Aligned_Sequences[i].substr(cut[point], cut[point + 1] - cut[point]), "-", "")<< PrimerR;
				}
				else {
					fasta_out << PrimerF <<replace(Aligned_Sequences[i].substr(cut[point], cut[point + 1] - cut[point]), "-", "") << PrimerR<< endl;
				}
			}
			else {
				outfragfile << PrimerF <<replace(Aligned_Sequences[i].substr(cut[point] - overlapvalue, cut[point + 1] - cut[point] + overlapvalue), "-", "") << PrimerR << "  Length without primers: " << replace(Aligned_Sequences[i].substr(cut[point] - overlapvalue, cut[point + 1] - cut[point] + overlapvalue), "-", "").length() << endl;
				if (i == Aligned_Sequences.size() - 1) {
					fasta_out << PrimerF <<replace(Aligned_Sequences[i].substr(cut[point] - overlapvalue, cut[point + 1] - cut[point] + overlapvalue), "-", "")<< PrimerR;
				}
				else {
					fasta_out << PrimerF <<replace(Aligned_Sequences[i].substr(cut[point] - overlapvalue, cut[point + 1] - cut[point] + overlapvalue), "-", "") << PrimerR << endl;
				}
			}
		}
		fasta_out.close();
		if (point == cut.size() - 2) {//write the fragment number for afterwards operation(oligo design).
			string fasta_num = dirname + "/" + UserName + "_FragCut_fragment.number";
			fstream number_out(fasta_num, iostream::out);
			number_out << point + 1;
			number_out.close();
		}
	}
	outfragfile.close();
	cout << "Finished" << endl;

	//mail to user
	cout << "Mail to user" << endl;
	if (mail_or_not == 1) {
		string align = "mail_dir_file.exe " + dirname + " " + Email + " " + UserName;
		char al[100];
		strcpy(al, align.c_str());
		system(al);
	}
	cout << "Finished" << endl;
	time_t t3 = clock();
	cout << "Total program costs " << (double)(t3 - t1) / CLOCKS_PER_SEC << " s" << endl;
	return 0;
}