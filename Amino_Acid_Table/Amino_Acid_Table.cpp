#define _CRT_SECURE_NO_WARNINGS
#include "Amino_Acid_Table.h"
#include <iostream>

void Amino_Acid_Table::scan_input_option(double codeFrequency_threshold) {
	for (int i = 0; i < 20; i++) {
		int k = 0;
		for (int j = 0; j < 6; j++) {
			if (Codon_Frequency[i][j] > codeFrequency_threshold) {
				Final_option[i][k] = Input_option[i][j];
				k++;
			}
		}
	}
}
Amino_Acid_Table::Amino_Acid_Table(char *filename)
{
	Read_Codon_Table(filename);
}

Amino_Acid_Table::~Amino_Acid_Table()
{
	//cout << "Delete Amino_Acid_Table Object" << endl;
}

void Amino_Acid_Table::Read_Codon_Table(const char* filename)
{
	FILE* fp;
	fp = fopen(filename, "r");
	float codon_freq_thousand[64];
	char str[30];
	char subbuff[6];
	float sum;
	int i = 0;
	while (1)
	{
		if (i % 16 == 0 && i != 0) {
			fgets(str, 2, fp);
		}
		fgets(str, 19, fp);
		memcpy(subbuff, &str[3], 5);
		subbuff[5] = '\0';
		codon_freq_thousand[i] = atof(subbuff);
		if (feof(fp))
		{
			break;
		}
		i++;
	}
	fclose(fp);
	//A 
	float gcu = codon_freq_thousand[49];
	float gcc = codon_freq_thousand[53];
	float gca = codon_freq_thousand[57];
	float gcg = codon_freq_thousand[61];
	sum = gcu + gcc + gca + gcg;
	Codon_Frequency[0][0] = gcu / sum;
	Codon_Frequency[0][1] = gcc / sum;
	Codon_Frequency[0][2] = gca / sum;
	Codon_Frequency[0][3] = gcg / sum;
	//R
	float cgu = codon_freq_thousand[19];
	float cgc = codon_freq_thousand[23];
	float cga = codon_freq_thousand[27];
	float cgg = codon_freq_thousand[31];
	float aga = codon_freq_thousand[43];
	int agg = codon_freq_thousand[47];
	sum = cgu + cgc + cga + cgg + aga + agg;
	Codon_Frequency[1][0] = cgu / sum;
	Codon_Frequency[1][1] = cgc / sum;
	Codon_Frequency[1][2] = cga / sum;
	Codon_Frequency[1][3] = cgg / sum;
	Codon_Frequency[1][4] = aga / sum;
	Codon_Frequency[1][5] = agg / sum;
	//N
	float aau = codon_freq_thousand[34];
	float aac = codon_freq_thousand[38];
	sum = aau + aac;
	Codon_Frequency[2][0] = aau / sum;
	Codon_Frequency[2][1] = aac / sum;
	//D
	float gau = codon_freq_thousand[50];
	float gac = codon_freq_thousand[54];
	sum = gau + gac;
	Codon_Frequency[3][0] = gau / sum;
	Codon_Frequency[3][1] = gac / sum;
	//C
	float ugu = codon_freq_thousand[3];
	float ugc = codon_freq_thousand[7];
	sum = ugu + ugc;
	Codon_Frequency[4][0] = ugu / sum;
	Codon_Frequency[4][1] = ugc / sum;
	//Q
	float caa = codon_freq_thousand[26];
	float cag = codon_freq_thousand[30];
	sum = caa + cag;
	Codon_Frequency[5][2] = caa / sum;
	Codon_Frequency[5][3] = cag / sum;
	//E
	float gaa = codon_freq_thousand[58];
	float gag = codon_freq_thousand[62];
	sum = gaa + gag;
	Codon_Frequency[6][2] = gaa / sum;
	Codon_Frequency[6][3] = gag / sum;
	//G
	float ggu = codon_freq_thousand[51];
	float ggc = codon_freq_thousand[55];
	float gga = codon_freq_thousand[59];
	float ggg = codon_freq_thousand[63];
	sum = ggu + ggc + gga + ggg;
	Codon_Frequency[7][0] = ggu / sum;
	Codon_Frequency[7][1] = ggc / sum;
	Codon_Frequency[7][2] = gga / sum;
	Codon_Frequency[7][3] = ggg / sum;

	//H
	float cau = codon_freq_thousand[18];
	float cac = codon_freq_thousand[22];
	sum = cau + cac;
	Codon_Frequency[8][0] = cau / sum;
	Codon_Frequency[8][1] = cac / sum;
	//I
	float auu = codon_freq_thousand[32];
	float auc = codon_freq_thousand[36];
	float aua = codon_freq_thousand[40];
	sum = auu + auc + aua;
	Codon_Frequency[9][0] = auu / sum;
	Codon_Frequency[9][1] = auc / sum;
	Codon_Frequency[9][2] = aua / sum;
	//L
	float uua = codon_freq_thousand[8];
	float uug = codon_freq_thousand[12];
	float cuu = codon_freq_thousand[16];
	float cuc = codon_freq_thousand[20];
	float cua = codon_freq_thousand[24];
	float cug = codon_freq_thousand[28];
	sum = uua + uug + cuu + cuc + cua + cug;
	Codon_Frequency[10][0] = cuu / sum;
	Codon_Frequency[10][1] = cuc / sum;
	Codon_Frequency[10][2] = cua / sum;
	Codon_Frequency[10][3] = cug / sum;
	Codon_Frequency[10][4] = uua / sum;
	Codon_Frequency[10][5] = uug / sum;

	//K
	float aaa = codon_freq_thousand[42];
	float aag = codon_freq_thousand[46];
	sum = aaa + aag;
	Codon_Frequency[11][2] = aag / sum;
	Codon_Frequency[11][3] = aaa / sum;

	//M
	Codon_Frequency[12][3] = 1;
	//F
	float uuu = codon_freq_thousand[0];
	float uuc = codon_freq_thousand[4];
	sum = uuu + uuc;
	Codon_Frequency[13][0] = uuc / sum;
	Codon_Frequency[13][1] = uuu / sum;

	//P
	float ccu = codon_freq_thousand[17];
	float ccc = codon_freq_thousand[21];
	float cca = codon_freq_thousand[25];
	float ccg = codon_freq_thousand[29];
	sum = ccu + ccc + cca + ccg;
	Codon_Frequency[14][0] = ccu / sum;
	Codon_Frequency[14][1] = cca / sum;
	Codon_Frequency[14][2] = ccg / sum;
	Codon_Frequency[14][3] = ccc / sum;
	//S
	float ucu = codon_freq_thousand[1];
	float ucc = codon_freq_thousand[5];
	float uca = codon_freq_thousand[9];
	float ucg = codon_freq_thousand[13];
	float agu = codon_freq_thousand[35];
	float agc = codon_freq_thousand[39];
	sum = ucu + ucc + uca + ucg + agu + agc;
	Codon_Frequency[15][0] = ucu / sum;
	Codon_Frequency[15][1] = ucc / sum;
	Codon_Frequency[15][2] = agc / sum;
	Codon_Frequency[15][3] = uca / sum;
	Codon_Frequency[15][4] = ucg / sum;
	Codon_Frequency[15][5] = agu / sum;
	//T
	float acu = codon_freq_thousand[33];
	float acc = codon_freq_thousand[37];
	float aca = codon_freq_thousand[41];
	float acg = codon_freq_thousand[45];
	sum = acu + acc + aca + acg;
	Codon_Frequency[16][0] = acu / sum;
	Codon_Frequency[16][1] = acc / sum;
	Codon_Frequency[16][2] = aca / sum;
	Codon_Frequency[16][3] = acg / sum;
	//W
	Codon_Frequency[17][3] = 1;
	//Y
	float uau = codon_freq_thousand[2];
	float uac = codon_freq_thousand[6];
	sum = uau + uac;
	Codon_Frequency[18][0] = uau / sum;
	Codon_Frequency[18][1] = uac / sum;
	//V
	float guu = codon_freq_thousand[48];
	float guc = codon_freq_thousand[52];
	float gua = codon_freq_thousand[56];
	float gug = codon_freq_thousand[60];
	sum = guu + guc + gua + gug;
	Codon_Frequency[19][0] = guu / sum;
	Codon_Frequency[19][1] = guc / sum;
	Codon_Frequency[19][2] = gua / sum;
	Codon_Frequency[19][3] = gug / sum;
	for (int i = 0; i < 20; i++) {
		for (int j = 0; j < 6; j++) {
			Codon_Frequency[i][j] = Codon_Frequency[i][j] * 100;
		}
	}
}