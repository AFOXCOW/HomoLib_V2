#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <string>
#include <fstream>
using namespace std;

string quota = "\"";
int main(int argc, char* argv[])
{
	string cmd;
	string email = argv[1];
	string username = argv[2];
	string resname = argv[3];
	string filename = argv[4];
	string dirname = argv[5];
	string customsite = argv[6];
	float codeFrequency = atof(argv[7]);
	string codon_frequency_file = argv[8];
	string rare_codon = argv[9];
	string reverse_translation_Primer_F = argv[10];
	string reverse_translation_Primer_R = argv[11];
	string mail_or_not = "0";//not mail util all progress finished.
	cmd = "Reverse_Translate.exe " + quota + email + quota + " " + quota + username + quota + " " + quota + resname + quota + " " + quota + filename + quota + " " + quota + dirname + quota + " " + quota + customsite + quota + " " + quota + to_string(codeFrequency) + quota + " " + quota + codon_frequency_file + quota+" "\
		+quota+ rare_codon+quota+" "+quota+ reverse_translation_Primer_F+quota+" "+quota+ reverse_translation_Primer_R+quota+" "+quota+mail_or_not+quota+" >./"+dirname+"/Gene_Design.log";
	char al[500];
	strcpy(al, cmd.c_str());
	system(al);
	

	int min_length = atoi(argv[12]);
	int max_length = atoi(argv[13]);
	string temp = "";
	temp = username;
	temp = "./" + dirname + "/" + temp + "_ReverseTranslate_resultGeneSequences.fasta";
	char tmp[100];
	strcpy(tmp, temp.c_str());
	filename = tmp; //逆翻译后的geneSeq文件名
	int overlapvalue = atoi(argv[14]);
	string fc_Primer_F = argv[15];
	string fc_Primer_R = argv[16];
	cmd = "FragmentCut.exe " + quota+email+quota + " " + quota+username+quota + " " + quota+to_string(min_length)+quota + " " + quota+to_string(max_length) +quota+ " " + quota+filename +quota+ " " +\
		quota+to_string(overlapvalue)+quota + " " +quota+ dirname+quota+" "+quota+fc_Primer_F+quota+" "+quota+fc_Primer_R+quota+" "+quota+mail_or_not+quota + " >./" + dirname + "/Fragments_Design.log";
	strcpy(al, cmd.c_str());
	system(al);

	temp = username;
	temp = "./" + dirname + "/" + temp + "_FragCut_fragment.number";
	strcpy(tmp, temp.c_str());
	filename = tmp;
	int frag_number;
	fstream in_number(filename, iostream::in);
	in_number >> frag_number;

	int oligo_min_length = atoi(argv[17]);
	int oligo_max_length = atoi(argv[18]);
	int high_homo_len = atoi(argv[19]);
	int both_side = atoi(argv[20]);//0 for 3` side.
	string oligo_Primer_F = argv[21];
	string oligo_Primer_R = argv[22];

	for (int i = 0; i < frag_number; i++) {
		temp = username;
		temp = "./" + dirname + "/" + temp + "_FragCut_fragment" + to_string(i + 1) + ".fasta";
		strcpy(tmp, temp.c_str());
		filename = tmp;
		cmd = "OligoDesign.exe " + quota+email+quota + " " +quota+ username+quota + " " + quota+to_string(oligo_min_length) +quota+ " " + quota+to_string(oligo_max_length) +quota+ " " + quota+filename+quota + " " + quota+dirname +quota+ " " +\
			quota+to_string(i + 1)+quota+" "+quota+to_string(high_homo_len)+quota+" "+quota+to_string(both_side)+quota+" "+quota+ oligo_Primer_F+quota+" "+quota+ oligo_Primer_R+quota+" "+quota+mail_or_not+quota + " >./" + dirname + "/Oligos_Design.log";
		strcpy(al, cmd.c_str());
		system(al);
	}

	//mail files in dir
	string align = "mail_dir_file.exe " + dirname + " " + email + " " + username;
	strcpy(al, align.c_str());
	system(al);
}