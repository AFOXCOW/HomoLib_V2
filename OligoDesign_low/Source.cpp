#define _CRT_SECURE_NO_WARNINGS
#include<string>
#include"../utils/utils.h"
#include"../MSA/MSA.h"
#include<time.h>
#include <algorithm>

using namespace std;

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




bool FindCutSite(vector<string> Aligned_Sequences, int min_length, int max_length, int high_homo_length, vector<vector<double>> similar_table, vector<int>& cut, vector<int>& reversecut)
{	//[cut[0],cut[1])
	//the min_length should be greater than high homo length
	bool warn_flag=false;
	cut.push_back(0);
	reversecut.push_back(Aligned_Sequences[0].size());
	//more length restriction on minlength
	vector<int> cutRange = { Find_length_index_forward_all(Aligned_Sequences,0,min_length + 10) ,Find_length_index_forward_one(Aligned_Sequences,0,max_length) };
	
	if (cutRange[1] - cutRange[0] < 1) {
		warn_flag = true;
		cutRange[0] = cutRange[1] - (max_length - min_length);
	}

	int score = 0;
	int best_score = -1000000;
	int best_site = -1;
	//find forward oligo cut site.
	while (true) {
		for (int end = cutRange[0]; end < cutRange[1]; ++end) {
			//calculate homo score.
			for (int site = end - high_homo_length; site < end; ++site) {
				score += ScoreAlignedSite(Aligned_Sequences, site, similar_table);
			}
			if (score > best_score) {
				best_score = score;
				best_site = end;
			}
			score = 0;
		}
		cut.push_back(best_site);
		/*
		if (best_site == -1) {
			return warn_flag;
		}
		*/
		cutRange[0] = Find_length_index_forward_all(Aligned_Sequences, best_site, min_length + 10);
		cutRange[1] = Find_length_index_forward_one(Aligned_Sequences, best_site, max_length);
		//clear
		score = 0;
		best_score = -1000000;
		best_site = -1;
		if (cutRange[0] == -1) {
			break;
		}
		if (cutRange[1] - cutRange[0] < 1) {
			warn_flag = true;
			cutRange[0] = cutRange[1] - (max_length - min_length);
		}
	}

	//find reverse oligo cut site.
	for (int i = cut.size() - 1; i > 0;--i) {
		reversecut.push_back(cut[i] - high_homo_length);
	}
	return warn_flag;
}


void Cut2Oligo(vector<string> Aligned_Sequences, vector<int> cut, vector<int> reversecut, vector<vector<string>>& seq, vector<vector<string>>& reverseseq) {
	// 根据以上的切割位点，对每个蛋白质序列作正向切割和反向
	seq.resize(Aligned_Sequences.size());
	reverseseq.resize(Aligned_Sequences.size());
	int start = 0;
	for (int j = 0; j < Aligned_Sequences.size(); ++j) {
		for (int i = 0; i < cut.size() - 1; ++i) {
			string substring = Aligned_Sequences[j].substr(cut[i], cut[i + 1] - cut[i]);
			seq[j].push_back(replace(substring, "-", ""));
		}
		string lastString = Aligned_Sequences[j].substr(cut[cut.size() - 1], Aligned_Sequences[j].size() - cut[cut.size() - 1]);
		seq[j].push_back(replace(lastString, "-", ""));
	}
	sort(reversecut.begin(), reversecut.end(), less<int>());
	for (int j = 0; j < Aligned_Sequences.size(); ++j) {
		string firstString = Aligned_Sequences[j].substr(0, reversecut[0]);
		firstString = replace(firstString, "-", "");
		string cp_firstString = complementary(firstString);
		reverseseq[j].push_back(cp_firstString);
		for (int i = 0; i < cut.size() - 1; ++i) {
			string substring = Aligned_Sequences[j].substr(reversecut[i], reversecut[i + 1] - reversecut[i]);
			substring = replace(substring, "-", "");
			string cp_substring = complementary(substring);
			reverseseq[j].push_back(cp_substring);
		}
	}
}

int main(int argc, char* argv[]) {
	time_t t1 = clock();

	//pre info
	cout << "Loading Parameters" << endl;
	string Email = argv[1];
	string userName = argv[2];
	int min_length = atoi(argv[3]);
	int max_length = atoi(argv[4]);
	string FASTA_path = argv[5];
	string dirName = argv[6];
	int frag_Index = atoi(argv[7]);
	int high_homo_length = atoi(argv[8]);
	int homo_both_side = atoi(argv[9]);//1 for both side. 0 for 3` side. default 0.
	string PrimerF = argv[10];
	string PrimerR = argv[11];
	int mail_or_not = atoi(argv[12]);
	cout << "Finished" << endl;

	//MSA
	cout << "Processing MSA" << endl;
	FASTA_Seq Sequences;
	vector<string> Aligned_Sequences;
	vector<vector<double>>similar_table;
	double Avg_length = 0.0;
	double sum_length = 0.0;
	for (int i = 0; i < Aligned_Sequences.size(); i++) {
		sum_length += Aligned_Sequences[i].size();
	}
	Avg_length = sum_length / Aligned_Sequences.size();
	double expected_frag_length = (double(min_length) + double(max_length)) / 2;
	int fragmentnumber = Avg_length / expected_frag_length + 0.5;
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


	//Oligo Design
	cout << "Finding High Homologous Regions" << endl;
	vector<int> cut;
	vector<int> reversecut;
	vector<vector<string>> Seq;
	vector<vector<string>> reverseSeq;
	bool warn_flag = FindCutSite(Aligned_Sequences, min_length, max_length, high_homo_length, similar_table, cut, reversecut);
	cout << "Finished" << endl;
	if (homo_both_side == 1) {
		//once you choose this, plz make sure a greater min_length and high homo length.
		//cut site shift to the high homo region`s middle thus making 5` 3` both side high homo.	
		for (int i = 1; i < cut.size(); ++i) {
			cut[i] = cut[i] - high_homo_length / 2;
		}
		for (int i = 2; i < reversecut.size(); ++i) {
			reversecut[i] = reversecut[i] + high_homo_length / 2;
		}
	}

	cout << "Designing Oligos" << endl;
	Cut2Oligo(Aligned_Sequences, cut, reversecut, Seq, reverseSeq);
	cout << "Finished" << endl;


	//Save to File
	cout << "Save to file" << endl;
	string file = "./" + dirName + "/" + userName + "_OligoDesign_result_" + to_string(frag_Index) + ".txt";
	fstream outoligofile(file, iostream::out);
	outoligofile << "HomoLibrary Output - " << userName << endl;
	outoligofile << endl;
	time_t nowtime;
	nowtime = time(NULL);
	struct tm* local, * gm;
	local = localtime(&nowtime);
	outoligofile << "Job started on " << local->tm_mday << "/" << local->tm_mon + 1 << "/" << local->tm_year + 1900 << " at " << local->tm_hour << ":" << local->tm_min << ":" << local->tm_sec << endl << endl;
	outoligofile << "Job name: " << userName << endl << endl;
	outoligofile << "Email: " << Email << endl << endl;
	outoligofile << "Oligo Length: " << min_length << "--" << max_length << endl << endl;
	outoligofile << "High homologous Length: " << high_homo_length << endl << endl;
	outoligofile << "Fragment Index: " << frag_Index << endl << endl;
	outoligofile << "Primer-F: " << PrimerF << endl << endl;
	outoligofile << "Primer-R: " << PrimerR << endl << endl;
	if (homo_both_side == 0) {
		outoligofile << "Only Require 3` high homologous" << endl << endl;
	}
	else {
		outoligofile << "Require 3` and 5` high homologous" << endl << endl;
		outoligofile << "---------notice---------------" << endl;
		outoligofile << "The first forward oligo or last reverse oligo may be very short" << endl;
		outoligofile << "you may refine the short oligos manually" << endl;
	}
	if (warn_flag) {
		outoligofile << "----------------------------WARNING------------------------------------------" << endl;
		outoligofile << "Some Sequences are too short, thus making Oligos length constraint can not be satisfied!!" << endl;
		outoligofile << "But we still try to design oligos for other sequences." << endl;
		outoligofile << "Therefore for some short sequences, or some large sequences, you need to refine their oligos manully!!" << endl;
		outoligofile << "----------------------------WARNING------------------------------------------" << endl;
	}
	outoligofile << "Input DNA sequences: " << endl;
	for (int j = 0; j < Aligned_Sequences.size(); j++)
	{
		outoligofile << Sequences.seq_names[j] << "(Length = " << replace(Aligned_Sequences[j], "-", "").size() << "): ";
		outoligofile << replace(Aligned_Sequences[j], "-", "") << endl;
	}
	outoligofile << endl << endl;
	outoligofile << "The Result of Aligned Sequences :" << endl;
	for (int j = 0; j < Aligned_Sequences.size(); j++)
	{
		outoligofile << Sequences.seq_names[j] << "(Length = " << Aligned_Sequences[j].size() << "): ";
		outoligofile << Aligned_Sequences[j] << endl;
	}
	outoligofile << endl << endl;
	outoligofile << "The Result of Designed Oligos :" << endl;
	for (int i = 0; i < Aligned_Sequences.size(); i++) {
		outoligofile << "Fragment" << endl;
		outoligofile << replace(Aligned_Sequences[i], "-", "") << endl << endl;
		outoligofile << "Forward Oligo" << endl;
		outoligofile << "5`----seq----3`" << endl;
		for (int j = 0; j < Seq[i].size(); j++) {
			outoligofile << PrimerF << Seq[i][j] << PrimerR << endl;;
		}
		outoligofile << endl;
		outoligofile << "Reverse Oligo" << endl;
		outoligofile << "5`----seq----3`" << endl;
		for (int j = 0; j < reverseSeq[i].size(); j++) {
			outoligofile << PrimerF << reverseSeq[i][j] << PrimerR << endl;
		}
		outoligofile << endl;
	}
	outoligofile.close();
	cout << "Finished" << endl;

	//mail
	cout << "Mail to user" << endl;
	if (mail_or_not == 1) {
		string align = "mail_dir_file.exe " + dirName + " " + Email + " " + userName;
		char al[100];
		strcpy(al, align.c_str());
		system(al);
	}
	cout << "Finished" << endl;
	time_t t3 = clock();
	cout << "Total Program costs " << (double)(t3 - t1) / CLOCKS_PER_SEC << " s" << endl;
	return 0;
}
