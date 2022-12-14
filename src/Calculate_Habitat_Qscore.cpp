#pragma once
#include<map>
#include<math.h>
#include<string>
#include<vector>
#include<sstream>
#include<fstream>
#include<assert.h>
#include<stdlib.h>
#include<string.h>
#include<iostream>
#include<sys/stat.h>
#include<omp.h>
#include<string.h>
using namespace std;


char Ref_db;
string Mode;
double weight[3] = { 1,1,1 };
string taxonomy_list_path = "";
string Qscore_path = "";



class Score_table {
public:
	string otu_num;
	double senstivity[32];
	double Wprecision[32];

	void set_score_table(string otu_num, double s[])
	{
		this->otu_num = otu_num;
		for (int i = 0; i < 64; i++)
		{
			if (i % 2 == 0)
				this->senstivity[i / 2] = s[i];
			else if (i % 2 == 1)
				this->Wprecision[i / 2] = s[i];
		}
	}


};

// String splited
vector<string> split(const string& str, const string& delim) {
	vector<string> res;
	if ("" == str) return res;
	char* strs = new char[str.length() + 1]; 
	strcpy(strs, str.c_str());

	char* d = new char[delim.length() + 1];
	strcpy(d, delim.c_str());

	char* p = strtok(strs, d);
	while (p) {
		string s = p; 
		res.push_back(s); 
		p = strtok(NULL, d);
	}

	return res;
}


template <class Type>
Type stringToNum(const string& str)
{
	istringstream iss(str);
	Type num;  
	iss >> num;
	return num;
}


void printhelp() {
	printf("Comp-taxa version : 1.0\r\n");
	printf("Compute the best amplication stagety among habitat types\r\n");
	printf("Usage:\r\n");
	printf("Qscore [Option vaule]\r\n");
	printf("Options:\r\n");
	printf("\t-D (upper) ref database, default is R (NCBI-Refseq), or G (GreenGenes-13-8), or C (GreenGenes-13-8-99), or S (SILVA)");
	printf("\t[Input options, required]\r\n");
	printf("\t  -i taxonomy name and abundance\r\n\tor\r\n");


	//printf("\t  -R If the input table is reversed, T(rue) or F(alse), default is false [Optional for -T]");

	printf("\t[Output options]\r\n");
	printf("\t  -o Output file, default is to output on screen\r\n");

	printf("\t[Other options]\r\n");
	printf("\t  -w alignment, Wprecision and cost weighted, default is 1 1 1\r\n");
	printf("\t  -h Help\r\n");

	exit(0);

}

int Parse_Para(int argc, char* argv[]) {

	Ref_db = 'R';

	Qscore_path = "../example/Qscore.txt";



	int i = 1;

	if (argc == 1)
		printhelp();

	while (i < argc) {
		if (argv[i][0] != '-') {
			cerr << "Argument # " << i << " Error : Arguments must start with -" << endl;
			exit(0);
		};
		switch (argv[i][1]) {
		case 'D': Ref_db = argv[i + 1][0]; break;

		case 'i': taxonomy_list_path = argv[i + 1];  Mode = "queryByTaxonomy"; break;

		case 'w': weight[0] = stringToNum<double>(argv[i + 1]);  weight[1] = stringToNum<double>(argv[i + 2]);  weight[2] = stringToNum<double>(argv[i + 3]); i += 2; Mode = 1; break;

		case 'o': Qscore_path = argv[i + 1]; break;


		case 'h': printhelp(); break;
		default: cerr << "Error: Unrec argument " << argv[i] << endl; printhelp(); break;
		}
		i += 2;
	}

	//int max_core_number = sysconf(_SC_NPROCESSORS_CONF);

	//if ((Coren <= 0) || (Coren > max_core_number)) {
	//	//cerr << "Core number must be larger than 0, change to automatic mode" << endl;
	//	Coren = max_core_number;
	//}

}


void Read_Score_tables(map<string, Score_table>& score_table) {
	string filepath;
	switch (Ref_db) {
	case 'G': filepath = "../Database/" + Mode + "/GG_13.txt"; break;
	case 'R': filepath = "../Database/" + Mode + "/ncbi_refseq.txt"; break;
	case 'C': filepath = "../Database/" + Mode + "/GG_13_99.txt"; break;
	case 'S': filepath = "../Database/" + Mode + "/silva_16s.txt"; break;
	default:  filepath = "../Database/" + Mode + "/ncbi_refseq.txt"; break;
	}
	ifstream infile;
	infile.open(filepath.c_str(), ios::in);
	if (!infile.is_open())
		printf("Can't open %s\n", filepath.c_str());
	else
	{
		string read_row;
		Score_table sctemp;
		vector<string> row_temp;
		double s[64];
		getline(infile, read_row);
		while (getline(infile, read_row)) {
			row_temp = split(read_row, "\r\n");
			row_temp = split(row_temp[0], "\t");
			for (int i = 1; i <= 64; i++)
				s[i - 1] = stringToNum<double>(row_temp[i]);
			string otu_num = row_temp[0];
			sctemp.set_score_table(otu_num, s);
			score_table.insert(pair<string, Score_table>(sctemp.otu_num, sctemp));
		}
		printf("Read Score table done\n");
	}
	infile.close();
}




void Read_Taxonomy_list(string& filepath, map<string, double>& taxonomy_list) {

	ifstream infile;
	infile.open(filepath.c_str(), ios::in);
	if (!infile.is_open())
		printf("Can't open %s\n", filepath.c_str());
	else
	{
		string read_row;
		vector<string> row_temp;
		double count_sum = 0;
		while (getline(infile, read_row)) {
			if(read_row[0]!='#'){
				row_temp = split(read_row, "\r\n");
				row_temp = split(row_temp[0], "\t");
				if (row_temp.size() >= 2) {
					taxonomy_list.insert(pair<string, double>(row_temp[0], stringToNum<double>(row_temp[1])));
					count_sum += stringToNum<double>(row_temp[1]);
				}
			}

			
		}

		for (map<string, double>::iterator it = taxonomy_list.begin(); it != taxonomy_list.end(); it++) {
			it->second /= count_sum;
		}
		printf("Read taxonomy list done\n");
	}
	infile.close();
}

//Ð´ÈëÎÄµµ
void write_score(string filepath, double senstivity[32], double acc[32], double Qscore[32]) {
	FILE* fp;

	if ((fp = fopen(filepath.c_str(), "a")) == NULL) {
		printf("Open Failed.\n");
		return;
	}
	else {
		double max_value = 0, max_value2 = 0;
		string region[32] = { "V1 100bp S", "V1 150bp S", "V1 300bp S", "V1V2 600bp P", "V1V3 600bp P", "V3 100bp S",
		 	"V3 150bp S", "V3 300bp S", "V3 300bp P", "V3V4 600bp P", "V3V5 600bp P", "V4 100bp S",
			"V4 150bp S", "V4 300bp S", "V4 300bp P", "V4V5 600bp P", "V4V6 600bp P", "V5 100bp S",
		   	"V5 150bp S", "V5 200bp P", "V5 300bp S", "V5V6 300bp P", "V6 100bp S", "V6 150bp S",
		    "V6 200bp P", "V6 300bp S", "V6V8 600bp P", "V7 100bp S", "V7 150bp S", "V7 300bp S", "V7V8 600bp P", "V7V9 600bp P" };
		fprintf(fp, "Region\t");
		for (int i = 0; i < 32; i++) {
			fprintf(fp, "%s\t", region[i].c_str());
		}
		fprintf(fp, "\n");


		fprintf(fp, "senstivity\t");
		for (int i = 0; i < 32; i++) {
			fprintf(fp, "%f\t", senstivity[i]);
		}
		fprintf(fp, "\n");
		fprintf(fp, "Wprecision\t");
		for (int i = 0; i < 32; i++) {
			fprintf(fp, "%f\t", acc[i]);
		}
		fprintf(fp, "\n");

		fprintf(fp, "Qscore\t");
		for (int i = 0; i < 32; i++) {
			fprintf(fp, "%f\t", Qscore[i]);
			if (max_value < Qscore[i])
				max_value = Qscore[i];

		}



	}


	fclose(fp);
}

void Qscore_mode(int length[32], double senstivity[32], double Wprecision[32], double value[32], map<string, double>& taxonomy_list, map<string, double>& unmap_taxonomy_list, map<string, Score_table>& score_table) {
	double sum_abu = 0;
	for (map<string, double>::iterator tl_it = taxonomy_list.begin(); tl_it != taxonomy_list.end(); tl_it++) {
		map<string, Score_table>::iterator sco_it = score_table.find(tl_it->first);
		if (sco_it != score_table.end()) {
			sum_abu += tl_it->second;
			for (int i = 0; i < 32; i++) {
				if (sco_it->second.senstivity[i] != 0) {
					senstivity[i] += (sco_it->second.senstivity[i] * tl_it->second);
					Wprecision[i] += (sco_it->second.Wprecision[i] * tl_it->second);
				}
			}

		}
		else {
			unmap_taxonomy_list.insert(*tl_it);
		}
	}
	
		for (int i = 0; i < 32; i++) {			
			Wprecision[i] /= sum_abu;
			senstivity[i] /= sum_abu;
			
				if (senstivity[i] != 0) {
					Wprecision[i] /= senstivity[i];
				}
				else Wprecision[i] = 0;
			
			value[i] = (weight[0] * senstivity[i] + weight[1] * Wprecision[i] + weight[2] * (1000.0 / (1000.0 + length[i]))) / (weight[0] + weight[1] + weight[2]);

		}
	
	
		
		

}

void Qscore() {

	int length[32] = { 100,150,300,600,600,100,150,300,300,600,600,100,150,300,300,600,600,
				100,150,200,300,300,100,150,200,300,600,100,150,300,600,600 };
	double senstivity[32] = { 0 };
	double Wprecision[32] = { 0 };
	double value[32] = { 0 };

	map<string, double> taxonomy_list, unmap_taxonomy_list;
	Read_Taxonomy_list(taxonomy_list_path, taxonomy_list);

	map<string, Score_table> score_table;
	Read_Score_tables(score_table);

	Qscore_mode(length, senstivity, Wprecision, value, taxonomy_list, unmap_taxonomy_list, score_table);



	write_score(Qscore_path, senstivity, Wprecision, value);

	for (int l = 0; l < 32; l++) {
		senstivity[l] = 0;
		Wprecision[l] = 0;
		value[l] = 0;
	}

	taxonomy_list.clear();


}


int main(int argc, char* argv[])
{
	Parse_Para(argc, argv);
	if (taxonomy_list_path != "")
		Qscore();

	return 0;
}

