#include "Reader.h"

#include <fstream>

int ReadFastaFile(char* input_fasta, std::vector<BioSeq*>& result) {
	std::filebuf fb;
	if (fb.open(input_fasta, std::ios::in)) {
		std::istream is(&fb);
		while (is) {
			std::string line;
			std::getline(is, line);
			if (line.size() && line[0] == '>') result.push_back(new BioSeq(line, ""));
			else if (line.size()) result.back()->nt_seq += line;
		}
		fb.close();
		return 0;
	} 
	return 1;
}