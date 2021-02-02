#include "Remove_ResSite.h"
using namespace std;
/*
R=A+G
Y=C+T
M=A+C
K=G+T
S=G+C
W=A+T
H=A+T+C
B=G+T+C
D=G+A+T
V=G+A+C
N=A+C+G+T
*/

string Degenerate_Bases = "RYMKSWHBDVN";
string   str[] = { "AG", "CT", "AC", "GT", "GC", "AT", "ATC", "GTC", "GAT", "GAC","ACGT" };
vector<string> Replace_Degenerate(str, str + sizeof(str) / sizeof(str[0]));

//AC(N)5CTCC ==> ACNNNNNCTCC
//in Database, there is no restriction site which have more than one (); 
//if there is more than one ,you can run this function several times.
string SiteNumberParse(string RawSite) {
	string numberChar = "0123456789";
	string out;
	int start = 0, number = 0, end = 0;
	char Degenerate;
	start = RawSite.find('(');
	Degenerate = RawSite[start + 1];
	if (start == -1) {
		return RawSite;
	}
	end = start + 2;
	int k = end + 1;
	string strNumber;
	while (numberChar.find(RawSite[k]) != numberChar.npos) {
		strNumber.push_back(RawSite[k]);
		++k;
	}
	number = atoi(strNumber.c_str());
	for (int i = 0; i < start; i++) {
		out.push_back(RawSite[i]);
	}
	for (int i = 0; i < number; i++) {
		out.push_back(Degenerate);
	}
	for (int i = k; i < RawSite.size(); i++) {
		out.push_back(RawSite[i]);
	}
	return out;
}


void FindResIndex(string res, string seq, vector<int> &index,vector<int> &length) {
	string pattern;
	for (int i = 0; i < res.size(); ++i) {
		if (Degenerate_Bases.find(res[i]) != Degenerate_Bases.npos) {
			//found Degenerate Bases
			pattern = pattern + '[' + str[Degenerate_Bases.find(res[i])] + ']';
			
		}
		else {
			pattern.push_back(res[i]);
		}
	}
	regex word_regex(pattern);
	auto words_begin = sregex_iterator(seq.begin(), seq.end(), word_regex);
	auto words_end = sregex_iterator();
	for (sregex_iterator i = words_begin; i != words_end; ++i) {
		smatch match = *i;
		index.push_back(match.position(0));
		length.push_back(match.str().size());
	}
}

bool NoResSite(string CheckString, vector<string> res) {
	vector<int>index,length;
	for (int i = 0; i < res.size(); ++i) {
		FindResIndex(res[i], CheckString, index, length);
		if (!index.empty()) {
			return false;
		}
	}
	return true;
}

static int best_score = 0;
void Iteration(string Protein,string ChangeString, string CheckString,int iter_level, Amino_Acid_Table table, vector<string> NoNumberRes, int start_Amino_Acid_inCheckStr, int end_Amino_Acid_inCheckStr,vector<vector<string>>Options) {
	if (iter_level == 0) {
		//calculate temp score
		int curr_score = 0;
		for (int i = 0; i < ChangeString.size(); ++i) {
			if (temp_tcga[i] == ChangeString[i]) {
				curr_score++;
			}
		}
		string ReplacedString = CheckString.replace(start_Amino_Acid_inCheckStr, end_Amino_Acid_inCheckStr - start_Amino_Acid_inCheckStr + 1, temp_tcga);
		if (!NoResSite(CheckString, NoNumberRes)) {
			curr_score = 0;
		}
		if (curr_score > best_score) {
			best_score = curr_score;
			best_tcga = temp_tcga;
		}
	}
	else {
		int Protein_idx;
		for (int i = 0; i < 20; ++i) {
			if (Protein[Protein.size() - iter_level] == table.AAOrder[i]) {
				Protein_idx = i;
				break;
			}
		}
		for (int i = 0; i < 6; ++i) {
			if (Options[Protein_idx][i] == "") {
				break;
			}
			temp_tcga.append(Options[Protein_idx][i]);
			Iteration(Protein, ChangeString,CheckString, iter_level - 1,table, NoNumberRes, start_Amino_Acid_inCheckStr, end_Amino_Acid_inCheckStr,Options);
			temp_tcga.pop_back();
			temp_tcga.pop_back();
			temp_tcga.pop_back();
		}
	}
}
string Choose_BestOne(string CheckString, int start_Amino_Acid_inCheckStr,int end_Amino_Acid_inCheckStr, Amino_Acid_Table table, vector<string> NoNumberRes,int tol) {
	string ChangeString = CheckString.substr(start_Amino_Acid_inCheckStr, end_Amino_Acid_inCheckStr - start_Amino_Acid_inCheckStr + 1);
	string Protein;
	string codon;
	bool Found = false;
	vector<vector<string>>Options;
	Options.resize(20);
	for (int i = 0; i < 20; ++i) {
		Options[i].resize(6);
	}
	if (tol == 1) {
		for (int i = 0; i < 20; ++i) {
			int k = 0;
			for (int j = 0; j < table.Input_option[i].size(); ++j) {
				if (table.Input_option[i][j] != " ") {
					Options[i][k] = table.Input_option[i][j];
					k++;
				}
			}
		}
	}
	else {
		for (int i = 0; i < 20; ++i) {
			for (int j = 0; j < 6; ++j) {
				Options[i][j] = table.Final_option[i][j];
			}
		}
	}
	//translate
	for (int i = 0; i < ChangeString.size(); i = i + 3) {
		codon = ChangeString.substr(i, 3);
		for (int k = 0; k < 20; ++k) {
			for (int j = 0; j < 6; ++j) {
				if (codon == Options[k][j]) {
					Protein.push_back(table.AAOrder[k]);
					Found = true;
					break;
				}
			}
			if (Found) {
				Found = false;
				break;
			}
		}
	}

	//Permutation and combination, and find best one.
	Iteration(Protein, ChangeString, CheckString, Protein.size(),table, NoNumberRes, start_Amino_Acid_inCheckStr, end_Amino_Acid_inCheckStr,Options);
	CheckString.replace(start_Amino_Acid_inCheckStr, end_Amino_Acid_inCheckStr - start_Amino_Acid_inCheckStr + 1, best_tcga);
	best_tcga.clear();
	temp_tcga.clear();
	best_score = 0;
	return CheckString;
}


int ReplaceResSites(vector<string> res, string& seq, Amino_Acid_Table table,int tol)
{
	int max_length = 0,curr_length = 0;
	vector<int> res_index;
	vector<string> NoNumberRes;
	for (int i = 0; i < res.size(); ++i) {
		NoNumberRes.push_back(SiteNumberParse(res[i]));
	}

	for (int i = 0; i < NoNumberRes.size(); ++i) {
		curr_length = NoNumberRes[i].size();
		if (curr_length > max_length) {
			max_length = curr_length;
		}
	}


	//find res site and store in index and length.
	vector<int> index, length;
	for (int i = 0; i < NoNumberRes.size(); ++i) {
		FindResIndex(NoNumberRes[i], seq, index, length);
	}

	//sort by index
	int swap_temp;
	for (int i = 0; i < index.size(); ++i) {
		for (int j = i+1; j < index.size(); ++j) {
			if (index[i] > index[j]) {
				swap_temp = index[i];
				index[i] = index[j];
				index[j] = swap_temp;

				swap_temp = length[i];
				length[i] = length[j];
				length[j] = swap_temp;
			}
		}
	}

	//replace one by one
	for (int i = 0; i < index.size(); ++i) {
		int start = index[i];
		int start_Amino_Acid = index[i] - index[i] % 3;
		while (i+1<index.size()-1 && index[i + 1] - index[i] - length[i] < max_length - 1) {
			//two resSite too close. take it as one.
			i++;
		}
		int end = index[i] + length[i];
		int end_Amino_Acid = index[i] + length[i] +3 - (index[i] + length[i]) % 3;
		int sub_start = 0, sub_end = 0;
		if (start - max_length - 1 < 0) {
			sub_start = 0;
		}
		else {
			sub_start = start - max_length - 1;
		}
		if (end + max_length -1 >seq.size()-1) {
			sub_end = seq.size()-1;
		}
		else {
			sub_end = end + max_length - 1;
		}

		//start at sub_start,and end at sub_end(include)
		string BestOne_Check;
		string CheckString = seq.substr(sub_start, sub_end- sub_start+1);
		int before_replace_len = CheckString.size();
		int start_Amino_Acid_inCheckStr = start_Amino_Acid - sub_start;
		int end_Amino_Acid_inCheckStr = end_Amino_Acid - sub_start-1;
		BestOne_Check = Choose_BestOne(CheckString, start_Amino_Acid_inCheckStr, end_Amino_Acid_inCheckStr, table, NoNumberRes,tol);
		int after_replace_len = BestOne_Check.size();
		if (before_replace_len != after_replace_len) {
			cout << "Can Not Remove Restriction Site while keeping the Protein Sequences the same!" << endl;
			exit(-1);
		}
		seq.replace(sub_start, sub_end - sub_start + 1, BestOne_Check);
	}
	return index.size();
}

/*
string Degenerate_Bases = "RYMKSWHBDVN";
int marka = 0, markr = 0, markn = 0, markd = 0, markc = 0, markq = 0, marke = 0, markg = 0, markh = 0, marki = 0, markl = 0,\
			markk = 0, markm = 0, markf = 0, markp = 0, marks = 0, markt = 0, markw = 0, marky = 0, markv = 0, markz = 0;



string   str[] = { "AG", "CT", "AC", "GT", "GC", "AT", "ATC", "GTC", "GAT", "GAC","ACGT" };
vector<string> Replace_Degenerate(str, str + sizeof(str) / sizeof(str[0]));
vector<string> ResConvert(string In_Site) {
	vector<string> out_sites;
	string ParsedSite = SiteNumberParse(In_Site);
	for (int i = 0; i < In_Site.size(); ++i) {
		if (Degenerate_Bases.find(In_Site[i]) == Degenerate_Bases.npos) {
			//Not found
			continue;
		}
		else {

		}
	}
}

//some input restriction site contain M,Y etc.

R=A+G
Y=C+T
M=A+C
K=G+T
S=G+C
W=A+T
H=A+T+C
B=G+T+C
D=G+A+T
V=G+A+C
N=A+C+G+T

vector<string> ResSitesConvert(vector<string> Input_ResSites) {
	vector<string> Output_ResSites;
	string In_ResSite;
	for (int i = 0; i < Input_ResSites.size(); i++) {
		In_ResSite = Input_ResSites[i];
	}
}











int toint(string s) {
	int res = 0;
	for (int i = 0; i < s.size(); i++){
		res *= 10;
		res += s[i] - '0';
	}
	return res;
}
int change_res_coden(string& seq, int sta) {
	string s;
	string cod = seq.substr(sta, 3);//²»Ó¦¸Ã¼õ1
	if (strcmp(cod.c_str(), "GCT") == 0 || strcmp(cod.c_str(), "GCC") == 0 || strcmp(cod.c_str(), "GCA") == 0 || strcmp(cod.c_str(), "GCG") == 0) {
		if (cod[2] == 'T') { seq.replace(sta + 2, 1, "A"); }
		else if (cod[2] == 'C') { seq.replace(sta + 2, 1, "G"); }
		else if (cod[2] == 'A') { seq.replace(sta + 2, 1, "T"); }
		else { seq.replace(sta + 2, 1, "C"); }
	}
	else if (strcmp(cod.c_str(), "CGT") == 0 || strcmp(cod.c_str(), "CGC") == 0 || strcmp(cod.c_str(), "CGA") == 0 || strcmp(cod.c_str(), "CGG") == 0 || strcmp(cod.c_str(), "AGA") == 0 || strcmp(cod.c_str(), "AGG") == 0) {
		if (cod[0] == 'C') {
			if (cod[2] == 'C') { seq.replace(sta + 2, 1, "G"); }
			else if (cod[2] == 'A') { seq.replace(sta + 2, 1, "T"); }
			else if (cod[2] == 'G') { seq.replace(sta + 2, 1, "C"); }
			else { seq.replace(sta + 2, 1, "A"); }
		}
		else {
			if (cod[2] == 'A') { seq.replace(sta, 1, "C"); }//AGA->CAG
			else { seq.replace(sta + 2, 1, "C"); }//AGG->AGA
		}
		// string optionr[]={"CGT","CGC","CGA","CGG","AGA","AGG"};
	}
	else if (strcmp(cod.c_str(), "AAT") == 0 || strcmp(cod.c_str(), "AAC") == 0) {
		if (cod[2] == 'T') { seq.replace(sta + 2, 1, "C"); }
		else { seq.replace(sta + 2, 1, "T"); }
		// string optionn[]={"AAT","AAC"}; 
	}
	else if (strcmp(cod.c_str(), "GAT") == 0 || strcmp(cod.c_str(), "GAC") == 0) {
		if (cod[2] == 'T') { seq.replace(sta + 2, 1, "C"); }
		else { seq.replace(sta + 2, 1, "T"); }
		// string optiond[]={"GAT","GAC"};
	}
	else if (strcmp(cod.c_str(), "TGT") == 0 || strcmp(cod.c_str(), "TGC") == 0) {
		if (cod[2] == 'T') { seq.replace(sta + 2, 1, "C"); }
		else { seq.replace(sta + 2, 1, "T"); }
		//string optionc[]={"TGT","TGC"};
	}
	else if (strcmp(cod.c_str(), "CAA") == 0 || strcmp(cod.c_str(), "CAG") == 0) {
		if (cod[2] == 'A') { seq.replace(sta + 2, 1, "G"); }
		else { seq.replace(sta + 2, 1, "A"); }
		//string optionq[]={"CAA","CAG"};
	}
	else if (strcmp(cod.c_str(), "GAA") == 0 || strcmp(cod.c_str(), "GAG") == 0) {
		if (cod[2] == 'A') { seq.replace(sta + 2, 1, "G"); }
		else { seq.replace(sta + 2, 1, "A"); }
		// string optione[]={"GAA","GAG"};
	}
	else if (strcmp(cod.c_str(), "GGT") == 0 || strcmp(cod.c_str(), "GGC") == 0 || strcmp(cod.c_str(), "GGA") == 0 || strcmp(cod.c_str(), "GGG") == 0) {
		if (cod[2] == 'T') { seq.replace(sta + 2, 1, "A"); }
		else if (cod[2] == 'A') { seq.replace(sta + 2, 1, "T"); }
		else if (cod[2] == 'G') { seq.replace(sta + 2, 1, "C"); }
		else { seq.replace(sta + 2, 1, "G"); }
		//string optiong[]={"GGT","GGC","GGA","GGG"};
	}
	else if (strcmp(cod.c_str(), "CAT") == 0 || strcmp(cod.c_str(), "CAC") == 0) {
		if (cod[2] == 'T') { seq.replace(sta + 2, 1, "C"); }
		else { seq.replace(sta + 2, 1, "T"); }
		//string optionh[]={"CAT","CAC"};
	}
	else if (strcmp(cod.c_str(), "ATT") == 0 || strcmp(cod.c_str(), "ATC") == 0 || strcmp(cod.c_str(), "ATA") == 0) {
		if (cod[2] == 'T') { seq.replace(sta + 2, 1, "A"); }
		else if (cod[2] == 'C') {
			if (marki == 0) { seq.replace(sta + 2, 1, "T"); marki = 1; }
			else { seq.replace(sta + 2, 1, "A"); marki = 0; }
		}
		else { seq.replace(sta + 2, 1, "T"); }
		// string optioni[]={"ATT","ATC","ATA"};
	}
	else if (strcmp(cod.c_str(), "TTA") == 0 || strcmp(cod.c_str(), "TTG") == 0 || strcmp(cod.c_str(), "CTT") == 0 || strcmp(cod.c_str(), "CTC") == 0 || strcmp(cod.c_str(), "CTA") == 0 || strcmp(cod.c_str(), "CTG") == 0) {
		if (cod[0] == 'C') {
			if (cod[2] == 'T') { seq.replace(sta + 2, 1, "A"); }
			else if (cod[2] == 'A') { seq.replace(sta + 2, 1, "T"); }
			else if (cod[2] == 'G') { seq.replace(sta + 2, 1, "C"); }
			else { seq.replace(sta + 2, 1, "G"); }
		}
		else {
			if (cod[2] == 'A') {
				if (markl == 0) { seq.replace(sta, 1, "C"); markl = 1; }
				else { seq.replace(sta + 2, 1, "G"); markl = 0; }
			}
			else {
				if (markl == 0) { seq.replace(sta + 2, 1, "A"); markl = 1; }
				else { seq.replace(sta, 1, "C"); markl = 0; }
			}
		}
		//string optionl[]={"TTA","TTG","CTT","CTC","CTA","CTG"};
	}
	else if (strcmp(cod.c_str(), "AAA") == 0 || strcmp(cod.c_str(), "AAG") == 0) {
		if (cod[2] == 'A') { seq.replace(sta + 2, 1, "G"); }
		else { seq.replace(sta + 2, 1, "A"); }
		//string optionk[]={"AAA","AAG"};
	}
	else if (strcmp(cod.c_str(), "TTT") == 0 || strcmp(cod.c_str(), "TTC") == 0) {
		if (cod[2] == 'T') { seq.replace(sta + 2, 1, "C"); }
		else { seq.replace(sta + 2, 1, "T"); }
		// string optionf[]={"TTT","TTC"};
	}
	else if (strcmp(cod.c_str(), "CCT") == 0 || strcmp(cod.c_str(), "CCC") == 0 || strcmp(cod.c_str(), "CCA") == 0 || strcmp(cod.c_str(), "CCG") == 0) {
		if (cod[2] == 'A') { seq.replace(sta + 2, 1, "G"); }
		else { seq.replace(sta + 2, 1, "A"); }
		//string optionp[]={"CCT","CCC","CCA","CCG"};
	}
	else if (strcmp(cod.c_str(), "TCT") == 0 || strcmp(cod.c_str(), "TCC") == 0 || strcmp(cod.c_str(), "TCA") == 0 || strcmp(cod.c_str(), "TCG") == 0 || strcmp(cod.c_str(), "AGT") == 0 || strcmp(cod.c_str(), "AGC") == 0) {
		if (cod[0] == 'T') {
			if (cod[2] == 'T') { seq.replace(sta + 2, 1, "A"); }
			else if (cod[2] == 'A') { seq.replace(sta + 2, 1, "T"); }
			else if (cod[2] == 'G') { seq.replace(sta + 2, 1, "C"); }
			else { seq.replace(sta + 2, 1, "G"); }
		}
		else {
			if (cod[2] == 'T') {
				seq.replace(sta, 3, "TCT");
				//  seq.replace(i,3,"TCA");
			}
			else {
				seq.replace(sta, 3, "TCC");
				//  seq.replace(i,3,"TCG");
			}
		}
		//string options[]={"TCT","TCC","TCA","TCG","AGT","AGC"};
	}
	else if (strcmp(cod.c_str(), "ACT") == 0 || strcmp(cod.c_str(), "ACC") == 0 || strcmp(cod.c_str(), "ACA") == 0 || strcmp(cod.c_str(), "ACG") == 0) {
		if (cod[2] == 'T') { seq.replace(sta + 2, 1, "A"); }
		else if (cod[2] == 'C') { seq.replace(sta + 2, 1, "G"); }
		else if (cod[2] == 'G') { seq.replace(sta + 2, 1, "C"); }
		else { seq.replace(sta + 2, 1, "T"); }
		//string optiont[]={"ACT","ACC","ACA","ACG"};
	}
	else if (strcmp(cod.c_str(), "TAT") == 0 || strcmp(cod.c_str(), "TAC") == 0) {
		if (cod[2] == 'T') { seq.replace(sta + 2, 1, "C"); }
		else { seq.replace(sta + 2, 1, "T"); }
		// string optiony[]={"TAT","TAC"};
	}
	else if (strcmp(cod.c_str(), "GTT") == 0 || strcmp(cod.c_str(), "GTC") == 0 || strcmp(cod.c_str(), "GTA") == 0 || strcmp(cod.c_str(), "GTG") == 0) {
		if (cod[2] == 'T') { seq.replace(sta + 2, 1, "A"); }
		else if (cod[2] == 'C') { seq.replace(sta + 2, 1, "G"); }
		else if (cod[2] == 'G') { seq.replace(sta + 2, 1, "C"); }
		else { seq.replace(sta + 2, 1, "T"); }
		// string optionv[]={"GTT","GTC","GTA","GTG"};
	}
	else if (strcmp(cod.c_str(), "TAG") == 0 || strcmp(cod.c_str(), "TGA") == 0 || strcmp(cod.c_str(), "TAA") == 0) {
		if (cod[2] == 'G') { seq.replace(sta + 1, 2, "GA"); }
		else if (cod[1] == 'G') { seq.replace(sta + 1, 2, "AG"); }
		else {
			if (markz == 0) {
				seq.replace(sta + 2, 1, "G"); markz = 1;
			}
			else { seq.replace(sta + 1, 1, "G"); markz = 0; }
		}
		//string optionz[]={"TAG","TGA","TAA"};
	}
	return 1;
}
int findres(string& seq, string res) {
	int numN, markN = 0, numres = 0;
	int pos, posN;
	string resp1, resp2;
	pos = res.find('N');
	if (pos != string::npos) {
		posN = pos;
		markN = 1;
		resp1 = res.substr(0, pos - 1);
		numN = toint(res.substr(pos + 2, 1));
		if (res.size() > pos + 2) {
			resp2 = res.substr(pos + 3, res.size() - pos - 3);
		}
		else {
			resp2 = "";
		}
	}
	if (markN == 0) {
		vector<string> res1;
		vector<string> res2;
		for (int j = 0; j < res.size(); j++) {
			if (res[j] == 'R') {
				if (res1.size() == 0)
				{
					res1.push_back("A");
					res1.push_back("G");
				}
				else
				{
					for (int i = 0; i < res1.size(); i++)
					{
						res2.push_back(res1[i]);
					}
					res1.clear();
					for (int i = 0; i < res2.size(); i++)
					{
						res1.push_back(res2[i].append("A"));
						res1.push_back(res2[i].append("G"));
					}
				}
			}
			if (res[j] == 'Y') {
				if (res1.size() == 0)
				{
					res1.push_back("C");
					res1.push_back("T");
				}
				else
				{
					for (int i = 0; i < res1.size(); i++)
					{
						res2.push_back(res1[i]);
					}
					res1.clear();
					for (int i = 0; i < res2.size(); i++)
					{
						res1.push_back(res2[i].append("C"));
						res1.push_back(res2[i].append("T"));
					}
				}
			}
			if (res[j] == 'M') {
				if (res1.size() == 0)
				{
					res1.push_back("A");
					res1.push_back("C");
				}
				else
				{
					for (int i = 0; i < res1.size(); i++)
					{
						res2.push_back(res1[i]);
					}
					res1.clear();
					for (int i = 0; i < res2.size(); i++)
					{
						res1.push_back(res2[i].append("A"));
						res1.push_back(res2[i].append("C"));
					}

				}
			}
			if (res[j] == 'K') {
				if (res1.size() == 0)
				{
					res1.push_back("G");
					res1.push_back("T");
				}
				else
				{
					for (int i = 0; i < res1.size(); i++)
					{
						res2.push_back(res1[i]);
					}
					res1.clear();
					for (int i = 0; i < res2.size(); i++)
					{
						res1.push_back(res2[i].append("G"));
						res1.push_back(res2[i].append("T"));
					}
				}
			}
			if (res[j] == 'S') {
				if (res1.size() == 0)
				{
					res1.push_back("C");
					res1.push_back("G");
				}
				else
				{
					for (int i = 0; i < res1.size(); i++)
					{
						res2.push_back(res1[i]);
					}
					res1.clear();
					for (int i = 0; i < res2.size(); i++)
					{
						res1.push_back(res2[i].append("C"));
						res1.push_back(res2[i].append("G"));
					}
				}
			}
			if (res[j] == 'W') {
				if (res1.size() == 0)
				{
					res1.push_back("A");
					res1.push_back("T");
				}
				else
				{
					for (int i = 0; i < res1.size(); i++)
					{
						res2.push_back(res1[i]);
					}
					res1.clear();
					for (int i = 0; i < res2.size(); i++)
					{
						res1.push_back(res2[i].append("A"));
						res1.push_back(res2[i].append("T"));
					}
				}
			}
			if (res[j] == 'H') {
				if (res1.size() == 0)
				{
					res1.push_back("A");
					res1.push_back("C");
					res1.push_back("T");
				}
				else
				{
					for (int i = 0; i < res1.size(); i++)
					{
						res2.push_back(res1[i]);
					}
					res1.clear();
					for (int i = 0; i < res2.size(); i++)
					{
						res1.push_back(res2[i].append("A"));
						res1.push_back(res2[i].append("C"));
						res1.push_back(res2[i].append("T"));
					}
				}
			}
			if (res[j] == 'B') {
				if (res1.size() == 0)
				{
					res1.push_back("C");
					res1.push_back("G");
					res1.push_back("T");
				}
				else
				{
					for (int i = 0; i < res1.size(); i++)
					{
						res2.push_back(res1[i]);
					}
					res1.clear();
					for (int i = 0; i < res2.size(); i++)
					{
						res1.push_back(res2[i].append("C"));
						res1.push_back(res2[i].append("G"));
						res1.push_back(res2[i].append("T"));
					}
				}
			}
			if (res[j] == 'V') {
				if (res1.size() == 0)
				{
					res1.push_back("A");
					res1.push_back("C");
					res1.push_back("G");
				}
				else
				{
					for (int i = 0; i < res1.size(); i++)
					{
						res2.push_back(res1[i]);
					}
					res1.clear();
					for (int i = 0; i < res2.size(); i++)
					{
						res1.push_back(res2[i].append("A"));
						res1.push_back(res2[i].append("C"));
						res1.push_back(res2[i].append("G"));
					}
				}
			}
			if (res[j] == 'H') {
				if (res1.size() == 0)
				{
					res1.push_back("A");
					res1.push_back("G");
					res1.push_back("T");
				}
				else
				{
					for (int i = 0; i < res1.size(); i++)
					{
						res2.push_back(res1[i]);
					}
					res1.clear();
					for (int i = 0; i < res2.size(); i++)
					{
						res1.push_back(res2[i].append("A"));
						res1.push_back(res2[i].append("G"));
						res1.push_back(res2[i].append("T"));
					}
				}
			}
			if (res[j] == 'A' || res[j] == 'T' || res[j] == 'C' || res[j] == 'G')
			{
				if (res1.size() == 0)
				{
					if (res[j] == 'A')
						res1.push_back("A");
					if (res[j] == 'T')
						res1.push_back("T");
					if (res[j] == 'C')
						res1.push_back("C");
					if (res[j] == 'G')
						res1.push_back("G");
				}
				else
				{
					for (int i = 0; i < res1.size(); i++)
					{
						if (res[j] == 'A')
							res1[i].append("A");
						if (res[j] == 'T')
							res1[i].append("T");
						if (res[j] == 'C')
							res1[i].append("C");
						if (res[j] == 'G')
							res1[i].append("G");
					}
				}
			}
		}
		for (int i = 0; i < res1.size(); i++) {
			pos = seq.find(res1[i].c_str());
			while (pos != string::npos) {
				numres++;
				if ((pos + res1[i].size()) % 3 == 0)
				{
					//change_res_coden(seq,pos+res1[i].size()-3);
					change_res_coden(seq, pos + res1[i].size() - 6);
				}

				if ((pos + res1[i].size()) % 3 == 1)
				{
					// change_res_coden(seq,pos+res1[i].size()-1);
					change_res_coden(seq, pos + res1[i].size() - 4);
				}
				if ((pos + res1[i].size()) % 3 == 2)
				{
					//change_res_coden(seq,pos+res1[i].size()-2);
					change_res_coden(seq, pos + res1[i].size() - 5);
				}
				pos = seq.find(res1[i].c_str(), pos + 1);
			}
		}
	}
	if (markN = 1)
	{
		vector<string> res11;
		vector<string> res12;
		for (int j = 0; j < resp1.size(); j++) {
			if (resp1[j] == 'R') {
				if (res11.size() == 0)
				{
					res11.push_back("A");
					res11.push_back("G");
				}
				else
				{

					for (int i = 0; i < res11.size(); i++)
					{
						res12.push_back(res11[i]);
					}
					res11.clear();
					for (int i = 0; i < res12.size(); i++)
					{
						res11.push_back(res12[i].append("A"));
						res11.push_back(res12[i].append("G"));
					}

				}
			}
			if (resp1[j] == 'Y') {
				if (res11.size() == 0)
				{
					res11.push_back("C");
					res11.push_back("T");
				}
				else
				{

					for (int i = 0; i < res11.size(); i++)
					{
						res12.push_back(res11[i]);
					}
					res11.clear();
					for (int i = 0; i < res12.size(); i++)
					{
						res11.push_back(res12[i].append("C"));
						res11.push_back(res12[i].append("T"));
					}

				}
			}
			if (resp1[j] == 'M') {
				if (res11.size() == 0)
				{
					res11.push_back("A");
					res11.push_back("C");
				}
				else
				{

					for (int i = 0; i < res11.size(); i++)
					{
						res12.push_back(res11[i]);
					}
					res11.clear();
					for (int i = 0; i < res12.size(); i++)
					{
						res11.push_back(res12[i].append("A"));
						res11.push_back(res12[i].append("C"));
					}

				}
			}
			if (resp1[j] == 'K') {
				if (res11.size() == 0)
				{
					res11.push_back("G");
					res11.push_back("T");
				}
				else
				{

					for (int i = 0; i < res11.size(); i++)
					{
						res12.push_back(res11[i]);
					}
					res11.clear();
					for (int i = 0; i < res12.size(); i++)
					{
						res11.push_back(res12[i].append("G"));
						res11.push_back(res12[i].append("T"));
					}

				}
			}
			if (resp1[j] == 'S') {
				if (res11.size() == 0)
				{
					res11.push_back("C");
					res11.push_back("G");
				}
				else
				{

					for (int i = 0; i < res11.size(); i++)
					{
						res12.push_back(res11[i]);
					}
					res11.clear();
					for (int i = 0; i < res12.size(); i++)
					{
						res11.push_back(res12[i].append("C"));
						res11.push_back(res12[i].append("G"));
					}

				}
			}
			if (resp1[j] == 'W') {
				if (res11.size() == 0)
				{
					res11.push_back("A");
					res11.push_back("T");
				}
				else
				{

					for (int i = 0; i < res11.size(); i++)
					{
						res12.push_back(res11[i]);
					}
					res11.clear();
					for (int i = 0; i < res12.size(); i++)
					{
						res11.push_back(res12[i].append("A"));
						res11.push_back(res12[i].append("T"));
					}

				}
			}
			if (resp1[j] == 'H') {
				if (res11.size() == 0)
				{
					res11.push_back("A");
					res11.push_back("C");
					res11.push_back("T");
				}
				else
				{

					for (int i = 0; i < res11.size(); i++)
					{
						res12.push_back(res11[i]);
					}
					res11.clear();
					for (int i = 0; i < res12.size(); i++)
					{
						res11.push_back(res12[i].append("A"));
						res11.push_back(res12[i].append("C"));
						res11.push_back(res12[i].append("T"));
					}

				}
			}
			if (resp1[j] == 'B') {
				if (res11.size() == 0)
				{
					res11.push_back("C");
					res11.push_back("G");
					res11.push_back("T");
				}
				else
				{

					for (int i = 0; i < res11.size(); i++)
					{
						res12.push_back(res11[i]);
					}
					res11.clear();
					for (int i = 0; i < res12.size(); i++)
					{
						res11.push_back(res12[i].append("C"));
						res11.push_back(res12[i].append("G"));
						res11.push_back(res12[i].append("T"));
					}

				}
			}
			if (resp1[j] == 'V') {
				if (res11.size() == 0)
				{
					res11.push_back("A");
					res11.push_back("C");
					res11.push_back("G");
				}
				else
				{

					for (int i = 0; i < res11.size(); i++)
					{
						res12.push_back(res11[i]);
					}
					res11.clear();
					for (int i = 0; i < res12.size(); i++)
					{
						res11.push_back(res12[i].append("A"));
						res11.push_back(res12[i].append("C"));
						res11.push_back(res12[i].append("G"));
					}

				}
			}
			if (resp1[j] == 'H') {
				if (res11.size() == 0)
				{
					res11.push_back("A");
					res11.push_back("G");
					res11.push_back("T");
				}
				else
				{

					for (int i = 0; i < res11.size(); i++)
					{
						res12.push_back(res11[i]);
					}
					res11.clear();
					for (int i = 0; i < res12.size(); i++)
					{
						res11.push_back(res12[i].append("A"));
						res11.push_back(res12[i].append("G"));
						res11.push_back(res12[i].append("T"));
					}

				}
			}
			if (resp1[j] == 'A' || resp1[j] == 'T' || resp1[j] == 'C' || resp1[j] == 'G')
			{
				if (res11.size() == 0)
				{
					if (resp1[j] == 'A')
						res11.push_back("A");
					if (resp1[j] == 'T')
						res11.push_back("T");
					if (resp1[j] == 'C')
						res11.push_back("C");
					if (resp1[j] == 'G')
						res11.push_back("G");
				}
				else
				{

					for (int i = 0; i < res11.size(); i++)
					{
						if (resp1[j] == 'A')
							res11[i].append("A");
						if (resp1[j] == 'T')
							res11[i].append("T");
						if (resp1[j] == 'C')
							res11[i].append("C");
						if (resp1[j] == 'G')
							res11[i].append("G");
					}

				}
			}
		}
		vector<string> res21;
		vector<string> res22;
		for (int j = 0; j < resp2.size(); j++) {
			if (resp2[j] == 'R') {
				if (res21.size() == 0)
				{
					res21.push_back("A");
					res21.push_back("G");
				}
				else
				{

					for (int i = 0; i < res21.size(); i++)
					{
						res22.push_back(res21[i]);
					}
					res21.clear();
					for (int i = 0; i < res22.size(); i++)
					{
						res21.push_back(res22[i].append("A"));
						res21.push_back(res22[i].append("G"));
					}

				}
			}
			if (resp2[j] == 'Y') {
				if (res21.size() == 0)
				{
					res21.push_back("C");
					res21.push_back("T");
				}
				else
				{

					for (int i = 0; i < res21.size(); i++)
					{
						res22.push_back(res21[i]);
					}
					res21.clear();
					for (int i = 0; i < res22.size(); i++)
					{
						res21.push_back(res22[i].append("C"));
						res21.push_back(res22[i].append("T"));
					}

				}
			}
			if (resp2[j] == 'M') {
				if (res21.size() == 0)
				{
					res21.push_back("A");
					res21.push_back("C");
				}
				else
				{

					for (int i = 0; i < res21.size(); i++)
					{
						res22.push_back(res21[i]);
					}
					res21.clear();
					for (int i = 0; i < res22.size(); i++)
					{
						res21.push_back(res22[i].append("A"));
						res21.push_back(res22[i].append("C"));
					}

				}
			}
			if (resp2[j] == 'K') {
				if (res21.size() == 0)
				{
					res21.push_back("G");
					res21.push_back("T");
				}
				else
				{

					for (int i = 0; i < res21.size(); i++)
					{
						res22.push_back(res21[i]);
					}
					res21.clear();
					for (int i = 0; i < res22.size(); i++)
					{
						res21.push_back(res22[i].append("G"));
						res21.push_back(res22[i].append("T"));
					}

				}
			}
			if (resp2[j] == 'S') {
				if (res21.size() == 0)
				{
					res21.push_back("C");
					res21.push_back("G");
				}
				else
				{

					for (int i = 0; i < res21.size(); i++)
					{
						res22.push_back(res21[i]);
					}
					res21.clear();
					for (int i = 0; i < res22.size(); i++)
					{
						res21.push_back(res22[i].append("C"));
						res21.push_back(res22[i].append("G"));
					}

				}
			}
			if (resp2[j] == 'W') {
				if (res21.size() == 0)
				{
					res21.push_back("A");
					res21.push_back("T");
				}
				else
				{
					for (int i = 0; i < res21.size(); i++)
					{
						res22.push_back(res21[i]);
					}
					res21.clear();
					for (int i = 0; i < res22.size(); i++)
					{
						res21.push_back(res22[i].append("A"));
						res21.push_back(res22[i].append("T"));
					}
				}
			}
			if (resp2[j] == 'H') {
				if (res21.size() == 0)
				{
					res21.push_back("A");
					res21.push_back("C");
					res21.push_back("T");
				}
				else
				{
					for (int i = 0; i < res21.size(); i++)
					{
						res22.push_back(res21[i]);
					}
					res21.clear();
					for (int i = 0; i < res22.size(); i++)
					{
						res21.push_back(res22[i].append("A"));
						res21.push_back(res22[i].append("C"));
						res21.push_back(res22[i].append("T"));
					}
				}
			}
			if (resp2[j] == 'B') {
				if (res21.size() == 0)
				{
					res21.push_back("C");
					res21.push_back("G");
					res21.push_back("T");
				}
				else
				{
					for (int i = 0; i < res21.size(); i++)
					{
						res22.push_back(res21[i]);
					}
					res21.clear();
					for (int i = 0; i < res22.size(); i++)
					{
						res21.push_back(res22[i].append("C"));
						res21.push_back(res22[i].append("G"));
						res21.push_back(res22[i].append("T"));
					}
				}
			}
			if (resp2[j] == 'V') {
				if (res21.size() == 0)
				{
					res21.push_back("A");
					res21.push_back("C");
					res21.push_back("G");
				}
				else
				{
					for (int i = 0; i < res21.size(); i++)
					{
						res22.push_back(res21[i]);
					}
					res21.clear();
					for (int i = 0; i < res22.size(); i++)
					{
						res21.push_back(res22[i].append("A"));
						res21.push_back(res22[i].append("C"));
						res21.push_back(res22[i].append("G"));
					}
				}
			}
			if (resp2[j] == 'H') {
				if (res21.size() == 0)
				{
					res21.push_back("A");
					res21.push_back("G");
					res21.push_back("T");
				}
				else
				{
					for (int i = 0; i < res21.size(); i++)
					{
						res22.push_back(res21[i]);
					}
					res21.clear();
					for (int i = 0; i < res22.size(); i++)
					{
						res21.push_back(res22[i].append("A"));
						res21.push_back(res22[i].append("G"));
						res21.push_back(res22[i].append("T"));
					}
				}
			}
			if (resp2[j] == 'A' || resp2[j] == 'T' || resp2[j] == 'C' || resp2[j] == 'G')
			{
				if (res21.size() == 0)
				{
					if (resp2[j] == 'A')
						res21.push_back("A");
					if (resp2[j] == 'T')
						res21.push_back("T");
					if (resp2[j] == 'C')
						res21.push_back("C");
					if (resp2[j] == 'G')
						res21.push_back("G");
				}
				else
				{

					for (int i = 0; i < res21.size(); i++)
					{
						if (resp2[j] == 'A')
							res21[i].append("A");
						if (resp2[j] == 'T')
							res21[i].append("T");
						if (resp2[j] == 'C')
							res21[i].append("C");
						if (resp2[j] == 'G')
							res21[i].append("G");
					}

				}
			}
		}
		for (int i = 0; i < res11.size(); i++) {
			//cout<<res11[i]<<endl;
			pos = seq.find(res11[i].c_str());
			while (pos != string::npos) {
				for (int j = 0; j < res21.size(); j++) {
					if (pos + res11[i].size() + numN > seq.size() - 1) { break; }
					if (res21[j] == seq.substr(pos + res11[i].size() + numN, res21[j].size())) {
						numres++;
						if ((pos + res11[i].size() + numN + res21[j].size()) % 3 == 0)
						{
							//change_res_coden(seq,pos+res11[i].size()+numN+res21[j].size());
							change_res_coden(seq, pos + res11[i].size() + numN + res21[j].size() - 6);
						}
						if ((pos + res11[i].size() + numN + res21[j].size()) % 3 == 1)
						{
							// change_res_coden(seq,pos+res11[i].size()+numN+res21[j].size()-1);
							change_res_coden(seq, pos + res11[i].size() + numN + res21[j].size() - 4);
						}
						if ((pos + res11[i].size() + numN + res21[j].size()) % 3 == 2)
						{
							//change_res_coden(seq,pos+res11[i].size()+numN+res21[j].size()-2);
							change_res_coden(seq, pos + res11[i].size() + numN + res21[j].size() - 5);
						}
						break;
					}
				}
				pos = seq.find(res11[i].c_str(), pos + 1);
			}
		}
	}
	return numres;
}
*/