#ifndef _Profile_H_
#define _Profile_H_

#include <vector>
#include <algorithm>

#include "BioSeq.h"
#include "PairwiseAlign.h"

struct Profile {
	// data
	std::vector<BioSeq*> sequences;
	int frequency[128];
	// constructor
	Profile() { };
	// destructor
	~Profile() {
		// everything will be deleted in main program
		/*
		for_each(sequences.begin(), sequences.end(), 
			[] (BioSeq* seq) { delete seq; } );
		*/
	};
	// alignment of two profiles
	Profile& operator + (Profile& another);
	// alignment of profile and sequence
	Profile& operator + (BioSeq* sequence);
	// get score in column
	float ColumnNTscore(BioSeq* seq, int index1, int index2);
	float ColumnAAscore(BioSeq* seq, int index1, int index2);
	float ColumnNTscore(Profile& another, int index1, int index2);
	float ColumnAAscore(Profile& another, int index1, int index2);
	// fill an array of frequency
	void CalcFrequenciesNT(int position);
	void CalcFrequenciesAA(int position);
	// insert gaps in sequences
	void InsertGap(int pos, int count) {
		for_each(sequences.begin(), sequences.end(), [pos, count] (BioSeq* s) {
			if (count) s->nt_seq.insert(pos, count, '-');
		});
	}
	void InsertGap(int pos) {
		for_each(sequences.begin(), sequences.end(), [pos] (BioSeq* s) {
			s->nt_seq.insert(pos, 1, '-');
		});
	}
};

// distance function type
typedef std::function<int (BioSeq* seq1, BioSeq* seq2)> func;

// UPGMA function change input sequences!
void UPGMA(std::vector<BioSeq*>& sequences, func distance, 
					 // input parameters
					 int go, int ge, int gf, int sc, const int* nt, const int* aa);
void UPGMAfree();

#endif