#define _CRT_SECURE_NO_WARNINGS
#include "../MSA/MSA.h"
#include<iostream>
#include<string>
#include<time.h>
#include "../Amino_Acid_Table/Amino_Acid_Table.h"
#include "mysql.h"
#include "../DataBase/DataBase.h"
#include "../AAs2DNA/AAs2DNA.h"
#include "../utils/utils.h"
#include "../Remove_RestrictionSite/Remove_ResSite.h"
#pragma comment(lib,"libmysql.lib")

using namespace std;
int main(int argc, char* argv[]) {
	time_t t1 = clock();
	cout << "Loading user`s Parameters" << endl;
	//Pre Info
	string FASTA_filePath;
	int pos;
	char* Codon_FrequecyTable_Path;
	string Email = argv[1];
	string userName = argv[2];
	string restrictionName = argv[3];
	FASTA_filePath = argv[4];
	string dirName = argv[5];
	string customResSite = argv[6];
	double codonFrequency_threshold = atof(argv[7]);
	Codon_FrequecyTable_Path = argv[8];
	int tol = atoi(argv[9]);//1 for using rare codon when removing restriction site, others for not.
	string PrimerF = argv[10];
	string PrimerR = argv[11];
	int mail_or_not = atoi(argv[12]); //1 for mail.
	string host = "temp";
	if (strcmp(Codon_FrequecyTable_Path, "./E_coli.txt") == 0) {
		host = "E.coli";
	}
	else if (strcmp(Codon_FrequecyTable_Path, "./P_pastoris.txt") == 0) {
		host = "P.pastoris";
	}
	else if (strcmp(Codon_FrequecyTable_Path, "./Bacillus_subtilis.txt") == 0) {
		host = "Bacillus_subtilis";
	}
	else if (strcmp(Codon_FrequecyTable_Path, "./S_cerevisiae.txt") == 0) {
		host = "S.cerevisiae";
	}
	else {
		host = "User-Defined";
	}
	cout << "Finished" << endl;
	//Multiple Sequences Alignment
	cout << "Processing Multiple Sequences Alignment" << endl;
	FASTA_Seq Sequences;
	vector<string> Aligned_Sequence;
	vector<vector<double>>similar_table;
	string final_Protein_List = "CSTPAGNDEQHRKMILVFYW-";
	int Number_of_row = final_Protein_List.size();
	int Number_of_column = final_Protein_List.size();
	similar_table.resize(Number_of_row);
	for (int i = 0; i < Number_of_row; i++) {
		similar_table[i].resize(Number_of_column);
	}

	string Blosum_filePath = "Blosum.txt";
	Sequences = Read_FASTA(FASTA_filePath);
	similar_table = Read_Similar_Table(Blosum_filePath, final_Protein_List);
	MSA MSA_object(Sequences.sequences, similar_table, final_Protein_List);
	MSA_object.IterationAlignProcess();
	Aligned_Sequence = MSA_object.output();
	cout << "Finished" << endl;
	time_t t2 = clock();
	cout << "MSA costs " << (double)(t2 - t1) / CLOCKS_PER_SEC << " s" << endl;

	//Read Host Codon Frequency Table From TxT File
	cout << "Reading Host Codon Frequency Table" << endl;
	Amino_Acid_Table aa_table(Codon_FrequecyTable_Path);
	aa_table.scan_input_option(codonFrequency_threshold * 100);// Means threshold *100%
	cout << "Finished" << endl;
	
	//Convert Protein Sequences To DNA Sequences
	cout << "Convert Protein to DNA." << endl;
	vector<string>DNA_seqs;
	DNA_seqs = reverseProteinToDNA(Aligned_Sequence, aa_table);
	cout << "Finished." << endl;
	time_t t3 = clock();
	cout << "Reverse Transltation costs " << (double)(t3 - t2) / CLOCKS_PER_SEC << " s" << endl;
	
	//Read Restriction Site Both From Database And also Custom Site
	cout << "Read Restriction Sites Both From Database And also Custom Sites." << endl;
	vector<string> Restriction_Sites;
	if (restrictionName.length() != 0)
	{
		int pos = restrictionName.find('|', 0);
		int star = 0, resnum = 0;
		string ResSite_inDB;
		while (pos != string::npos) {
			string ResName = restrictionName.substr(star, pos - star);
			ResSite_inDB = mysqlres(ResName);
			SplitString(ResSite_inDB, Restriction_Sites, ",");
			star = pos + 1;
			pos = restrictionName.find('|', star);
			resnum++;
		}
		ResSite_inDB = mysqlres(restrictionName.substr(star, restrictionName.size() - star));
		SplitString(ResSite_inDB, Restriction_Sites, ",");
		resnum++;
	}
	SplitString(customResSite, Restriction_Sites, "|");
	cout << "Finished." << endl;

	//Remove Restriction Site In DNA Sequences(optimize)
	cout << "Removing Restriction Sites" << endl;
	for (int i = 0; i < DNA_seqs.size(); i++) {
		cout<<"Restriction Site Count of Seq "<< i <<": "<< ReplaceResSites(Restriction_Sites, DNA_seqs[i],aa_table,tol)<<endl;
	}
	cout << "Finished." << endl;
	time_t t4 = clock();
	cout << "Removing Restriction Site costs " << (double)(t4 - t3) / CLOCKS_PER_SEC << " s" << endl;

	//Save to file
	cout << "Save to file" << endl;
	string file = "./"+dirName +"/"+ userName + "_ReverseTranslate_resultGeneSequences.txt";
	string fasta_file = "./" + dirName + "/" + userName + "_ReverseTranslate_resultGeneSequences.fasta";
	fstream output(file, iostream::out);
	fstream fasta_output(fasta_file, iostream::out);
	output << "HomoLibrary Output - " << userName << endl << endl;
	time_t nowtime;
	nowtime = time(NULL);
	struct tm* local, * gm;
	local = localtime(&nowtime);
	output << "Job started on " << local->tm_mday << "/" << local->tm_mon + 1 << "/" << local->tm_year + 1900 << " at " << local->tm_hour << ":" << local->tm_min << ":" << local->tm_sec << endl << endl;
	output << "Job name: " << userName << endl << endl;
	output << "Email: " << Email << endl << endl;
	output << "Restriction Site Name: " << restrictionName << endl << endl;
	output << "Custom Site: " << customResSite << endl << endl;
	output << "code frequency threshold: " << codonFrequency_threshold << endl << endl;
	output << "host: " << host << endl << endl;
	output << "Primer-F: " << PrimerF << endl << endl;
	output << "Primer-R: " << PrimerR << endl << endl;
	if (tol == 1) {
		output << "Using rare codon when necessary"<< endl << endl;
	}
	else {
		output << "Do not use rare codon in any condition" << endl << endl;
	}
	output << "Input Protein sequences: " << endl;
	for (int j = 0; j < Sequences.sequences.size(); j++)
	{
		output << Sequences.seq_names[j] << "(Length = " << Sequences.sequences[j].size() << "): ";
		output << Sequences.sequences[j];
		output << endl;
	}
	output << endl << endl;
	output << "The Result of Aligned Protein Sequences :" << endl;
	for (int j = 0; j < Aligned_Sequence.size(); j++)
	{
		output << Sequences.seq_names[j] << "(Length = " << Aligned_Sequence[j].size() << "): ";
		output << Aligned_Sequence[j];
		output << endl;
	}
	output << endl << endl;
	output << "The Result of designed DNA Sequences (with Primer) :" << endl;
	for (int k = 0; k < DNA_seqs.size(); k++)
	{
		output  << Sequences.seq_names[k] << "_DNA_(Length without primer = " << DNA_seqs[k].size() << "): " << endl;
		fasta_output << Sequences.seq_names[k] << "_DNA:" << endl;
		if (k == DNA_seqs.size() - 1) {
			fasta_output << PrimerF <<DNA_seqs[k]<<PrimerR;
			output << PrimerF <<DNA_seqs[k]<< PrimerR;
		}
		else {
			fasta_output << PrimerF <<DNA_seqs[k] << PrimerR << endl;
			output <<PrimerF<< DNA_seqs[k] << PrimerR << endl;
		}
	}
	output.close();
	fasta_output.close();
	cout << "Finished" << endl;

	//mail to user
	cout << "Mail to user" << endl;
	if (mail_or_not == 1) {
		string align = "mail_dir_file.exe " + dirName + " " + Email + " " + userName;
		char al[100];
		strcpy(al, align.c_str());
		system(al);
	}
	cout << "Finished" << endl;
	time_t t5 = clock();
	cout << "Total program costs " << (double)(t5 - t1) / CLOCKS_PER_SEC << " s" << endl;
	return 0;
}



















/*
void SplitString(const std::string& s, std::vector<std::string>& v, const std::string& c)
{
	std::string::size_type pos1, pos2;
	pos2 = s.find(c);
	pos1 = 0;
	while (std::string::npos != pos2)
	{
		v.push_back(s.substr(pos1, pos2 - pos1));

		pos1 = pos2 + c.size();
		pos2 = s.find(c, pos1);
	}
	if (pos1 != s.length())
		v.push_back(s.substr(pos1));
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
*/
/* An unknow error for this function*/
/*
int Save_to_File(FASTA_Seq Sequences,vector<string> Aligned_Protein_Seqs, vector<string> DNA_seqs,baseInfo info) {

	string file = info.dirName + info.userName + "_ReverseTranslate_resultGeneSequences.txt";
	string fasta_file = info.dirName + info.userName + "_ReverseTranslate_resultGeneSequences.fasta";

	fstream output(file, iostream::out);
	fstream fasta_output(fasta_file, iostream::out);
	output << "HomoLibrary Output - " << info.userName << endl<<endl;
	time_t nowtime;
	nowtime = time(NULL);
	struct tm* local, * gm;
	local = localtime(&nowtime);
	output << "Job started on " << local->tm_mday << "/" << local->tm_mon + 1 << "/" << local->tm_year + 1900 << " at " << local->tm_hour << ":" << local->tm_min << ":" << local->tm_sec << endl << endl;
	output << "Job name: " << info.userName << endl << endl;
	output << "Email: " << info.Email << endl << endl;
	output << "Restriction Site Name: " << info.restrictionName << endl << endl;
	output << "Custom Site: " << info.customResSite << endl << endl;
	output << "code frequency: " << info.codonFrequency_threshold << endl << endl;
	output << "host: " << info.host << endl << endl;
	output << "Input Protein sequences: " << endl;

	for (int j = 0; j < Aligned_Protein_Seqs.size(); j++)
	{
		output << "Protein Sequence " << j << "(Length = " << replace(Aligned_Protein_Seqs[j], "-", "").size() << "): " << replace(Aligned_Protein_Seqs[j], "-", "");
		output << endl;
	}
	output << endl << endl;
	output << "The result of designed DNA DNA_sequences :" << endl;

	for (int k = 0; k < DNA_seqs.size(); k++)
	{
		output << k << "-" << Sequences.seq_names[k] << "(Length = " << DNA_seqs[k].size() << "): " << endl;
		fasta_output << '>' << k << "-" << Sequences.seq_names[k] << "(Length = " << DNA_seqs[k].size() << "): " << endl;
		if (k == Aligned_Protein_Seqs.size() - 1) {
			fasta_output << DNA_seqs[k];
			output << DNA_seqs[k];
		}
		else {
			fasta_output << DNA_seqs[k] << endl;
			output << DNA_seqs[k] << endl;
		}
	}
	output.close();
	fasta_output.close();
	return 0;
}
*/