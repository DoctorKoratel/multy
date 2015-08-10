#ifndef _Reader_H_
#define _Reader_H_

#include <vector>

#include "BioSeq.h"

int ReadFastaFile(char* input_fasta, std::vector<BioSeq*>& buf);

#endif