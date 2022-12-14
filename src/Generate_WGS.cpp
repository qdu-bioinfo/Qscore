#include <vector>
#include <algorithm>
#include <map>
#include <math.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

using namespace std;

int Sequence_length = 300;
int Sliding_window_max_step = 2000;
int Mode = 0;
string write_type = "T";
string read_path = "", write_path = "Generate_genome";


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




class OTUs {
public:
	string OTUs_id;
	string OTUs_seq;

	void clear() {
		this->OTUs_id = "";
		this->OTUs_seq = "";
	}

	void set_num(string num) {
		this->OTUs_id = num;
	}
	void set_OTUs_seq(string OTUs_seq) {
		this->OTUs_seq = OTUs_seq;
	}
	void del_seq(int pos) {
		this->OTUs_seq.erase(pos, 1);
	}
	void ins_seq(int pos, string bp) {
		this->OTUs_seq.insert(pos, bp);
	}
	void replace_seq(int& pos, string& bp) {
		this->OTUs_seq.replace(pos, 1, bp);
	}
	int get_seq_length() {
		return OTUs_seq.size();
	}

};

void otu_add_errs(OTUs& OTUs_temp) {
	/*	Qvalue = 37 - 1.8 * 10^(-4) * i^2
	*	90% replace error
	*	10%	insert or delete error
	*/
	srand(time(NULL));
	double error1 = 0.10;
	double error2 = 0.50;
	double Qvalue;

	for (int i = 0; i < OTUs_temp.get_seq_length(); i++) {

		Qvalue = pow(10, -(int)(37 - 0.000188 * (i) * (i)+0.5) / 10.0);

		double te_er = rand() * 1.0 / RAND_MAX;

		if (te_er <= Qvalue) {

			double er1 = rand() / double(RAND_MAX);
			if (er1 <= error1) {
				er1 = rand() / double(RAND_MAX);
				double er2 = rand() / double(RAND_MAX);
				if (er2 < error2) {
					int ra = rand() % 4;
					string bp;
					if (ra == 0)
						bp = "A";
					else if (ra == 1)
						bp = "T";
					else if (ra == 2)
						bp = "G";
					else bp = "C";
					OTUs_temp.ins_seq(i, bp);
				}
				else OTUs_temp.del_seq(i);
			}
			else {
				string base_temp = "A";
				switch (rand() % 4) {
				case 0:
					base_temp = "A";
					break;
				case 1:
					base_temp = "T";
					break;
				case 2:
					base_temp = "G";
					break;
				case 3:
					base_temp = "C";
					break;
				default:
					base_temp = "A";
				}
				OTUs_temp.replace_seq(i, base_temp);

			}
		}

	}

}

void Read_fasta(string& filepath, vector<OTUs>& OTUs_list)
{

	ifstream infile;
	infile.open(filepath.c_str(), ios::in);

	if (!infile.is_open()) {
		printf("Can't open %s\n", filepath.c_str());
	}
	else {
		string read_row = "", Seq_temp = "";
		OTUs OTUs_temp;
		int count_num = 1;


		while (getline(infile, read_row)) {
			if (read_row[0] == '>') {
				if (Seq_temp.size() >= 100) {
					OTUs_temp.OTUs_seq = Seq_temp;
					OTUs_list.push_back(OTUs_temp);
					OTUs_temp.clear();
					count_num++;
				}
				OTUs_temp.OTUs_id = read_row.substr(1) + "_" + to_string(count_num);
				Seq_temp = "";
			}
			else {
				Seq_temp += read_row;
			}
		}
		OTUs_temp.OTUs_seq = Seq_temp;
		OTUs_list.push_back(OTUs_temp);

	}
	infile.close();
}

void Write_OTU_list(string& filepath, vector<OTUs>& otu_list) {

	FILE* fp;
	string write_file_path = "";

	if (write_type == "F") {
		write_file_path = filepath + ".fq";
	}
	else {
		write_file_path = filepath + ".fa";
	}

	if ((fp = fopen(write_file_path.c_str(), "a")) == NULL)
	{
		printf("Open Failed.\n");
		return;
	}
	else {
		srand(time(NULL));
		for (vector<OTUs>::iterator it = otu_list.begin(); it != otu_list.end(); it++) {

			int rand_int = rand() % Sliding_window_max_step;

			for (int i = rand_int; i < it->OTUs_seq.size(); i += rand_int) {
				OTUs OTUs_temp;
				OTUs_temp.OTUs_id = it->OTUs_id + "_site:" + to_string(i);
				OTUs_temp.OTUs_seq = it->OTUs_seq.substr(i, Sequence_length + 50);
				otu_add_errs(OTUs_temp);
				fprintf(fp, ">%s\n%s\n", OTUs_temp.OTUs_id.c_str(), OTUs_temp.OTUs_seq.substr(0, Sequence_length).c_str());
				if (write_type == "F") {
					fprintf(fp, "+\r\n");
					for (int j = 0; j < Sequence_length; j++)
						fprintf(fp, "%c", (int)(70 - 0.000188 * (j) * (j)+0.5));
					fprintf(fp, "\r\n");

				}
				rand_int = rand() % Sliding_window_max_step + 1;
			}



		}
	}

	fclose(fp);


}

void printhelp() {
	printf("Generate_genome version : 1.0\r\n");
	printf("Generate simulated metagenomic sequencing sequences\r\n");
	printf("Usage:\r\n");
	printf("Generate_genome [Option vaule]\r\n");
	printf("Options:\r\n");

	printf("\t[Input options, required]\r\n");
	printf("\t  -i full-length genome fasta file path\r\n\tor\r\n");
	printf("\t  -l full-length genome fasta file paths list\r\n");

	printf("\t[Output options]\r\n");
	printf("\t  -o Output file, default is Generate_genome\r\n");
	printf("\t  -f Output file type, T is fasta,F is fastq default is T\r\n");

	printf("\t[Other options]\r\n");
	printf("\t  -w Random sliding window max step, default is 2,000\r\n");
	printf("\t  -s Sequences length, default is 300\r\n");
	printf("\t  -h Help\r\n");

	exit(0);

}


int Parse_Para(int argc, char* argv[]) {

	int i = 1;

	if (argc == 1)
		printhelp();

	while (i < argc) {
		if (argv[i][0] != '-') {
			cerr << "Argument # " << i << " Error : Arguments must start with -" << endl;
			exit(0);
		};
		switch (argv[i][1]) {

		case 'i': read_path = argv[i + 1];  Mode = 1; break;

		case 'l': read_path = argv[i + 1]; Mode = 2; break;

		case 'o': write_path = argv[i + 1]; break;

		case 'w': Sliding_window_max_step = stringToNum<int>(argv[i + 1]); break;

		case 's': Sequence_length = stringToNum<int>(argv[i + 1]); break;

		case 'f': write_type = argv[i + 1]; break;

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


void Generate_genome_specimen() {

	vector<OTUs> database;
	if (Mode == 1) {
		Read_fasta(read_path, database);
		Write_OTU_list(write_path, database);
		database.clear();
	}
	else if (Mode == 2) {

		vector<string> read_path_list;

		ifstream infile;
		infile.open(read_path.c_str(), ios::in);

		if (!infile.is_open()) {
			printf("Can't open %s\n", read_path.c_str());
		}
		else {
			string read_row = "";

			while (getline(infile, read_row)) {
				read_row = split(read_row, "\r\n")[0];
				read_path_list.push_back(read_row);
			}

		}
		for (vector<string>::iterator it = read_path_list.begin(); it != read_path_list.end(); it++) {
			Read_fasta(*it, database);
			Write_OTU_list(write_path, database);
			database.clear();
		}

	}


}

int main(int argc, char* argv[])
{
	Parse_Para(argc, argv);
	if (read_path != "" && write_path != "")
		Generate_genome_specimen();

	return 0;
}




