#include <algorithm> 
#include <getopt.h>
#include <stdlib.h>
#include <unistd.h>

#include "BioSeq.h"
#include "Reader.h"
#include "Profile.h"
#include "PairwiseAlign.h"

using std::cout;
using std::endl;

// initial dimensions of the matrices; using in Profile.c++
int dim = 1024; 
// program arguments
char* input_fasta(NULL);
std::string NTsubst("NUC-45");
std::string AAsubst("BLOSUM62");
std::string algo("DEFAULT");
int gap_open(-10), gap_extension(-3), gap_frame(-15), stop_cost(-50), kmer(10);
int bonus(0); 

#define BORDER1 1000

#define BORDER2 15

void usage(char *name)
{
	cout << "usage: " << name << endl;
	cout << "short options (POSIX standart)" << endl;
	cout << "\t-h show this message" << endl;
	cout << "\t-i input file name (FASTA format expected)" << endl;
	cout << "\t-n filename with NT substitution matrix" << endl;
	cout << "\t\tDEFAULT: match +5 mismatch -4" << endl;
	cout << "\t-a filename with AA substitution matrix" << endl;
	cout << "\t\tDEFAULT: BLOSUM62 matrix" << endl;
	cout << "\t-g cost of creating a gap" << endl;
	cout << "\t\tDEFAULT: " << gap_open << endl;
	cout << "\t-e cost of a gap extension" << endl;
	cout << "\t\tDEFAULT: " << gap_extension << endl;
	cout << "\t-f cost of a gap frameshift " << endl;
	cout << "\t\tDEFAULT: " << gap_frame << endl;
	cout << "\t-s cost of a stop codon not at the end of the sequence" << endl;
	cout << "\t\tDEFAULT: " << stop_cost << endl;
	cout << "\t-b bonus for not break the frame" << endl;
	cout << "\t\tDEFAULT: " << bonus << endl;
	cout << "\t-d initial dimensions of the matrices" << endl;
	cout << "\t\tDEFAULT: 1024x1024" << endl;
	cout << "\t-k use k-mers in UPGMA" << endl;
	cout << "\t\tDEFAULT: k = " << kmer << endl;
	cout << "\t-p use pairwise aligner in UPGMA" << endl;
	cout << "long options" << endl;
	cout << "\t--help show this message" << endl;
	cout << "\t--input input file name (FASTA format expected)" << endl;
	cout << "\t--NT_subst filename with NT substitution matrix" << endl;
	cout << "\t\tDEFAULT: match +5 mismatch -4" << endl;
	cout << "\t--AA_subst filename with AA substitution matrix" << endl;
	cout << "\t\tDEFAULT: BLOSUM62 matrix" << endl;
	cout << "\t--gap_open cost of creating a gap" << endl;
	cout << "\t\tDEFAULT: " << gap_open << endl;
	cout << "\t--gap_extension cost of a gap extension" << endl;
	cout << "\t\tDEFAULT: " << gap_extension << endl;
	cout << "\t--gap_frame cost of a gap frameshift " << endl;
	cout << "\t\tDEFAULT: " << gap_frame << endl;
	cout << "\t--stop_cost cost of a stop codon not at the end of the sequence" << endl;
	cout << "\t\tDEFAULT: " << stop_cost << endl;
	cout << "\t--bonus bonus for not break the frame" << endl;
	cout << "\t\tDEFAULT: " << bonus << endl;
	cout << "\t--dimension initial dimensions of the matrices" << endl;
	cout << "\t\tDEFAULT: 1024x1024" << endl;
	cout << "\t--k-mers use k-mers in UPGMA" << endl;
	cout << "\t\tDEFAULT: k = " << kmer << endl;
	cout << "\t--pairwise use pairwise aligner in UPGMA" << endl;
}

void CheckSequences(std::vector<BioSeq*> &data) {
	for (auto index = data.begin(); index != data.end(); ) {
		BioSeq* seq = *index;
		bool flag = true;
		for (int i = 0; i < seq->Length(); i++) {
			char c = (*seq)[i];
			if (c == '\n' || c == '\r' || c == '\t' || c == ' ') {
				seq->nt_seq.erase(i, 1);
				i--;
			} else if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
				cout << "Unknown nuc '" << c << "' in sequence: " << seq->name << endl;
				cout << "Position: " << i << endl;
				cout << "Sequence removed" << endl;
				index = data.erase(index);
				flag = false;
				break;
			}
		}
		if (flag) index++;
	}
}

int main(int argc, char** argv) {
	static struct option opt[] = {
		{"help", 0, 0, 'h'},
		{"input", 1, 0, 'i'},
		{"NT_subst", 1, 0, 'n'},
		{"AA_subst", 1, 0, 'a'},
		{"gap_open", 1, 0, 'g'},
		{"gap_extension", 1, 0, 'e'},
		{"gap_frame", 1, 0, 'f'},
		{"stop_cost", 1, 0, 's'},
		{"k-mers", 2, 0, 'k'},
		{"pairwise", 0, 0, 'p'},
		{"dimension", 1, 0, 'd'},
		{"bonus", 1, 0, 'b'},
		{0,0,0,0}
	};
	int ind, code;
	opterr = 1;
	while ((code=getopt_long(argc,argv,"hpk::i:n:a:g:e:f:s:d:b:",opt,&ind)) != -1) 
	{
		switch (code) {
			case 'h':
				usage(argv[0]);
				return 0;
			case 'i':
				input_fasta = optarg;
				break;
			case 'n':
				NTsubst = optarg;
				break;
			case 'a':
				AAsubst = optarg;
				break;
			case 'g':
				gap_open = atoi(optarg);
				break;
			case 'e':
				gap_extension = atoi(optarg);
				break;
			case 'f':
				gap_frame = atoi(optarg);
				break;
			case 's':
				stop_cost = atoi(optarg);
				break;
			case 'k':
				algo = "KMERS";
				if (optarg) kmer = atoi(optarg);
				break;
			case 'p':
				algo = "PAIRWISE";
				break;
			case 'd':
				dim = atoi(optarg);
				if (dim < 1) dim = 1;
				break;
			case 'b':
				bonus = atoi(optarg);
				break;
			default:
				usage(argv[0]);
				return -1;
		}
	}
	if (!input_fasta) {
		cout << "nothing to align" << endl;
		return 0;
	}
	//===============
	cout << "Reading sequences..." << endl;
	std::vector<BioSeq*> data; 
	if (ReadFastaFile(input_fasta, data)) {
		cout << "can't open input file: " << input_fasta << endl;
		return -1;
	}
	CheckSequences(data);
	cout << data.size() << " sequences were obtained" << endl; 
	cout << "Input parameters:" << endl;
	cout << "NT substitution matrix\t" << NTsubst << endl;
	cout << "AA substitution matrix\t" << AAsubst << endl;
	cout << "Gap open cost\t\t" << gap_open << endl;
	cout << "Gap extension cost\t" << gap_extension << endl;
	cout << "Gap frame cost\t\t" << gap_frame << endl;
	cout << "Stop codon cost\t\t" << stop_cost << endl;
	// creating aligner tool
	PairwiseAlign aligner(NTsubst.c_str(), AAsubst.c_str(), stop_cost, gap_open, 
												gap_extension, gap_frame);
	// check substitution matrices
	if (!aligner.GetNTscoreMatrix()) {
		cout << "can't open NT substitution matrix file: " << NTsubst << endl;
		return -1;
	}
	if (!aligner.GetAAscoreMatrix()) {
		cout << "can't open AA substitution matrix file: " << AAsubst << endl;
		return -1;
	}
	// calculating MSA
	if (algo == "DEFAULT") {
		int average_length = 0;
		for_each(data.begin(), data.end(), [&average_length] (BioSeq* seq) {
			average_length += seq->Length();
		} );
		average_length /= data.size();
		if (average_length > BORDER1 || data.size() > BORDER2) 
			algo = "KMERS";
		else
			algo = "PAIRWISE";
	}
	func distance;
	if (algo == "PAIRWISE") 
		distance = [&aligner] (BioSeq* seq1, BioSeq* seq2) -> int { 
			return aligner.Align(seq1, seq2); 
		};
	else 
		distance = [] (BioSeq* seq1, BioSeq* seq2) -> int { 
			int score = 0;
			for (int i = 0; i < seq1->Length() - kmer + 1; i++) {
				for (int j = 0; j < seq2->Length() - kmer + 1; j++) {
					// test i-th k-mer in seq1 and j-th k-mer in seq2
					bool flag = true;
					for (int z = 0; flag && z < kmer; z++) {
						flag = ((*seq1)[i + z] == (*seq2)[j + z]);
					}
					if (flag) score++;
				}
			}
			return score;
		};
	// this call will change data!
	UPGMA(data, distance, gap_open, gap_extension, gap_frame, stop_cost, 
				aligner.GetNTscoreMatrix(), aligner.GetAAscoreMatrix()); 
	// print results
	cout << endl;
	for_each(data.begin(), data.end(), [] (BioSeq* seq) {seq->PrintNT(cout);} );
	cout << endl;
	for_each(data.begin(), data.end(), [] (BioSeq* seq) {seq->PrintAA(cout);} );
	// free memory
	for_each(data.begin(), data.end(), [] (BioSeq* seq) { delete seq; } );
	return 0;
}