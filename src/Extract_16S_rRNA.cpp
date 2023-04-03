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
#include <regex>
using namespace std;

smatch m;

static int Sequence_length = 300;
static int Mode = 0;
string direction = "T";
string write_type = "T";
string read_path = "", write_path = "Extract_16S_rRNA/";




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

class Extract_16S {

public:
	string OTU_num;
	string seqs[13] = {};
	int site[13] = {};	

	void clear() {
		OTU_num = "";
		for (int i = 0; i < 13; i++) {
			seqs[i] = "";
			site[i] = 0;
		}

	}
};

void convert_seq(string& seqs) {
	/*
	*	Convert sequence for reverse complementarity	
	*/
	string temp = "";
	for (int i = seqs.size() - 1; i >= 0; i--) {
		if (seqs.at(i) == 'A') {
			temp += 'T';
		}
		else if (seqs.at(i) == 'T') {
			temp += 'A';
		}
		else if (seqs.at(i) == 'G') {
			temp += 'C';
		}
		else if (seqs.at(i) == 'C') {
			temp += 'G';
		}
	}
	seqs = temp;
}

void sequencing_error(OTUs& OTUs_temp, int& region) {
	/*	
	*	Forward: Qvalue = 37 - 1.8 * 10^(-4) * i^2
	*	Reverse: Qvalue = 33 - 1.8 * 10^(-4) * i^2
	*	90% replace error
	*	10%	insert or delete error
	*/
	srand(time(NULL));
	double error1 = 0.10;
	double error2 = 0.50;
	double Qvalue;

	for (int i = 0; i < OTUs_temp.get_seq_length(); i++) {
		if (region < 6) {
			Qvalue = pow(10, -(int)(37 - 0.000188 * (i) * (i)+0.5) / 10.0);
		}
		else {
			Qvalue = pow(10, -(int)(33 - 0.000188 * (i) * (i)+0.5) / 10.0);
		}
		double te_er = rand() * 1.0 / RAND_MAX;

		if (te_er <= Qvalue) {

			double er1 = rand() / double(RAND_MAX);
			if (er1 <= error1) {
				er1 = rand() / double(RAND_MAX);
				double er2 = rand() / double(RAND_MAX);
				if (er2 < error2) {
					int rand_base = rand() % 4;
					string bp;
					if (rand_base == 0)
						bp = "A";
					else if (rand_base == 1)
						bp = "T";
					else if (rand_base == 2)
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







void Write_OTU_list(string& filepath, vector<Extract_16S>& otu_list, multimap<string, int>& target_list) {

	string region[13] = { "8F","341F","515F","U789F","967F","1099F",
		"357R","518R","806R","926R","1064R","1406R","1492R" };
	string write_path = filepath + "/16S_rRNA_gene_list.txt";
	int primer_length[13] = { 20,17,19,18,19,16,15,19,20,20,19,17,16 };

	string command = "mkdir -p " + filepath;
	system(command.c_str());

	FILE* fp;
	int count_num = 0;
	if ((fp = fopen(write_path.c_str(), "a")) == NULL)
	{
		printf("Open Failed.\n");
		return;
	}
	else {
		for (multimap<string, int>::iterator it = target_list.begin(); it != target_list.end(); it++) {
			fprintf(fp, "%s_site:_%d\n", it->first.c_str(), it->second);
		}

	}

	for (int i = 0; i < 13; i++) {
		if (write_type == "F") {
			write_path = filepath + "/" + region[i] + "_" + to_string(Sequence_length) + "bp.fq";
		}
		else {
			write_path = filepath + "/" + region[i] + "_" + to_string(Sequence_length) + "bp.fa";
		}

		FILE* fp;
		int count_num = 0;
		if ((fp = fopen(write_path.c_str(), "a")) == NULL)
		{
			printf("Open Failed.\n");
			return;
		}
		else {
			srand(time(NULL));

			for (vector<Extract_16S>::iterator it = otu_list.begin(); it != otu_list.end(); it++) {

				if (it->seqs[i].size() >= 100) {
					OTUs OTUs_temp;
					OTUs_temp.OTUs_id = it->OTU_num;
					OTUs_temp.OTUs_seq = it->seqs[i];
					sequencing_error(OTUs_temp, i);
					fprintf(fp, ">%s\n%s\n", OTUs_temp.OTUs_id.c_str(), OTUs_temp.OTUs_seq.substr(primer_length[i], Sequence_length - 30).c_str());
					if (write_type == "F") {
						fprintf(fp, "+\r\n");
						for (int j = 0; j < Sequence_length; j++) {
							if (i < 6)
								fprintf(fp, "%c", (int)(70 - 0.000188 * (j) * (j)+0.5));
							else
								fprintf(fp, "%c", (int)(66 - 0.000188 * (j) * (j) + 0.5));
						}

						fprintf(fp, "\r\n");

					}




				}


			}
		}


		fclose(fp);

	}
}


void printhelp() {
	printf("Extraxt 16S rRNA version : 1.0\r\n");
	printf("Extraxt 16S rRNA sequences\r\n");
	printf("Usage:\r\n");
	printf("Extraxt_16S_rRNA [Option vaule]\r\n");
	printf("Options:\r\n");

	printf("\t[Input options, required]\r\n");
	printf("\t  -i full-length genome fasta file path\r\n\tor\r\n");
	printf("\t  -l full-length genome fasta file paths list\r\n");

	printf("\t[Output options]\r\n");
	printf("\t  -o Output file, default is Extraxt_16S_rRNA\r\n");
	printf("\t  -f Output file type, T is fasta,F is fastq default is T\r\n");

	printf("\t[Other options]\r\n");
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

void Extraxt_16S_rRNA() {
	srand(time(NULL));
	string seq_Pri[13] = { "AG[AG]GTT[CT]GAT[CT][ATCG]TGGCTCAG","CCTACGGG[ATCG]GGC[AT]GCAG","GTG[CT]CAGC[AC]GCCGCGGTAA","TAGATACCC[ACTG][CG]GTAGTCC","CAACGCGAAGAACCTTACC","G[CT]AACGAGCGCAACCC",
		"TACGG[AG]AGGCAGCAG","GTGCCAGC[CA]GCCGCGGTAA","ATTAGA[AT]ACCC[TCG][ATGC]GTAGTCC","AAACT[TC]AAA[TG][AG]AATTG[AG]CGG","AGGTG[ATCG]TGCATGG[TC][TC]GTCG","TG[TC]AC[AT]CAC[TC]GCCCGTC","AAGTCGTAACAAGGTA" };

	string region[13] = { "8F","341F","515F","U789F","967F","1099F",
		"357R","518R","806R","926R","1064R","1406R","1492R" };
	int site[13] = { 8,341,515,789,967,1099,357,518,806,926,1064,1406,1492 };

	vector<OTUs> database;
	multimap<string, int> target_list;

	Extract_16S otu_temp;
	vector<Extract_16S> OTU_list;

	if (Mode == 1) {
		Read_fasta(read_path, database);

		for (vector<OTUs>::iterator db_it = database.begin(); db_it != database.end(); db_it++) {
			for (int i = 0; i < 13; i++) {
				regex is_Pri(seq_Pri[i]);
				bool found;
				string seqs_temp = db_it->OTUs_seq;
				int otu_count = 1;
				while (found = regex_search(seqs_temp, m, is_Pri)) {
					otu_temp.site[i] = db_it->OTUs_seq.length() - m.suffix().str().length();
					if (i < 6) {

						otu_temp.seqs[i] = m.str() + m.suffix().str().substr(0, Sequence_length);
					}
					else {
						int pos = m.prefix().str().size() - (Sequence_length);
						if (pos < 0)
							pos = 0;
						otu_temp.seqs[i] = m.prefix().str().substr(pos, Sequence_length) + m.str();
						convert_seq(otu_temp.seqs[i]);
					}

					otu_temp.OTU_num = db_it->OTUs_id + "_site:_" + to_string(otu_temp.site[i]);
					OTU_list.push_back(otu_temp);
					multimap<string, int>::iterator target_it_beg = target_list.lower_bound(db_it->OTUs_id);
					multimap<string, int>::iterator target_it_end = target_list.upper_bound(db_it->OTUs_id);
					if (target_it_beg != target_list.end()) {
						bool flag = true;
						for (; target_it_beg != target_it_end; target_it_beg++) {
							if (abs(otu_temp.site[i] - target_it_beg->second) < 1500) {
								flag = false;
								break;
							}
						}
						if (flag == true)
							target_list.insert(pair<string, int>(db_it->OTUs_id, otu_temp.site[i] - site[i]));
					}
					else
						target_list.insert(pair<string, int>(db_it->OTUs_id, otu_temp.site[i] - site[i]));

					otu_temp.clear();
					otu_count++;
					seqs_temp = m.suffix().str();

				}
			}
		}
		Write_OTU_list(write_path, OTU_list, target_list);
		database.clear();
		OTU_list.clear();
		target_list.clear();
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

			for (vector<OTUs>::iterator db_it = database.begin(); db_it != database.end(); db_it++) {
				for (int i = 0; i < 13; i++) {
					regex is_Pri(seq_Pri[i]);
					bool found;
					string seqs_temp = db_it->OTUs_seq;
					int otu_count = 1;
					while (found = regex_search(seqs_temp, m, is_Pri)) {

						otu_temp.site[i] = db_it->OTUs_seq.length() - m.suffix().str().length();
						if (i < 6) {

							otu_temp.seqs[i] = m.str() + m.suffix().str().substr(0, Sequence_length + 50);
						}
						else {
							int pos = m.prefix().str().size() - (Sequence_length + 50);
							if (pos < 0)
								pos = 0;
							otu_temp.seqs[i] = m.prefix().str().substr(pos, Sequence_length + 50) + m.str();
						}

						otu_temp.OTU_num = db_it->OTUs_id + "_site:_" + to_string(otu_temp.site[i]);
						OTU_list.push_back(otu_temp);
						multimap<string, int>::iterator target_it_beg = target_list.lower_bound(db_it->OTUs_id);
						multimap<string, int>::iterator target_it_end = target_list.upper_bound(db_it->OTUs_id);
						if (target_it_beg != target_list.end()) {
							bool flag = true;
							for (; target_it_beg != target_it_end; target_it_beg++) {
								if (target_it_beg->second - otu_temp.site[i] < 1500) {
									flag = false;
									break;
								}
							}
							if (flag == true)
								target_list.insert(pair<string, int>(db_it->OTUs_id, otu_temp.site[i]));
						}
						else
							target_list.insert(pair<string, int>(db_it->OTUs_id, otu_temp.site[i]));
						otu_temp.clear();
						otu_count++;
						seqs_temp = m.suffix().str();

					}
				}
			}

			Write_OTU_list(write_path, OTU_list, target_list);
			database.clear();
			OTU_list.clear();
		}
	}




}





int main(int argc, char* argv[])
{
	Parse_Para(argc, argv);
	if (read_path != "" && write_path != "")
		Extraxt_16S_rRNA();

	return 0;
}
