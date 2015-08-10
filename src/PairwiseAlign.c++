#include <fstream>
#include <string.h>
#include <algorithm>

#include "PairwiseAlign.h"

#define MAXINT 0x7FFFFFFF

PairwiseAlign :: PairwiseAlign(const char* nt_score_matrix, 
			const char* aa_score_matrix, int stop_cost, int gap_open, 
			int gap_extension, int gap_frame) {
	this->gap_open = gap_open;
	this->gap_extension = gap_extension;
	this->gap_frame = gap_frame;
	this->stop_cost = stop_cost;
	ChangeNTscoreMatrix(nt_score_matrix);
	ChangeAAscoreMatrix(aa_score_matrix);
}

void PairwiseAlign :: NewScoreMatrix(const char* file_name, int* score_matrix) {
	std::ifstream ifs;
	ifs.open (file_name, std::ifstream::in);
	if (!ifs.is_open()) return;
	// initialization
	for (int i = 0; i < 128; i++)
		for (int j = 0; j < 128; j++)
			score_matrix[i*128 + j] = -MAXINT;
	// skip comments
	char c;
	std::string comment_string;
	while (ifs.good() && (c = ifs.get()) == '#')
		std::getline(ifs, comment_string);
	ifs.putback(c);
	// reading alphabet
	std::string alphabet;
	std::getline(ifs, alphabet);
	alphabet.erase(std::remove_if(alphabet.begin(), alphabet.end(),
								[](char c){ return (c == ' ' || c == '\t'); }), alphabet.end());
	// reading matrix
	for (unsigned int i = 0; i < alphabet.length(); i++) {
		c = ifs.get();
		if (c != '*')
			for (unsigned int j = 0; j < alphabet.length(); j++) {
				if (alphabet[j] != '*') ifs >> score_matrix[c*128 + alphabet[j]];
				else {
					int value;
					ifs >> value;
					for (int z = 0; z < 128; z++)
						if (score_matrix[c*128 + z] == -MAXINT)
							score_matrix[c*128 + z] = value;
				}
			}
			else
				for (unsigned int j = 0; j < alphabet.length(); j++) {
					int value;
					ifs >> value;
					if (alphabet[j] != '*') {
						for (int z = 0; z < 128; z++)
							if (score_matrix[z*128 + alphabet[j]] == -MAXINT)
								score_matrix[z*128 + alphabet[j]] = value;
					}
					else {
						for (int z = 0; z < 128; z++)
							if (score_matrix[z*128 + z] == -MAXINT)
								score_matrix[z*128 + z] = value;
					}
				}
				std::getline(ifs, comment_string);
	}
	ifs.close();
	/*  DEBUG
	 *    for (int i = 0; i < 128; i++) {
	 *        for (int j = 0; j < 128; j++) {
	 *            std::cout << score_matrix[i*128 + j] << ' ';
	 *        }
	 *        std::cout << std::endl;
	 *    }
	 */
}

int PairwiseAlign :: Align(const BioSeq* s1, const BioSeq* s2) {
	n = s1->nt_seq.length() + 1;
	m = s2->nt_seq.length() + 1;
	// check borders
	int new_size = 1;
	while (n > new_size || m > new_size) {
		new_size *= 2;
	}
	// reallocate memory (if necessary)
	if (new_size > size) {
		if (F) delete[] F; if (W) delete[] W;
		F = new int[new_size * new_size];
		W = new int[new_size * new_size];
		size = new_size;
	}
	// initialization
	memset(F, 0, sizeof(int)*size*size);
	memset(W, 0, sizeof(int)*size*size);
	// calculating alignment score
	for (int i = 1; i < n; i++) {
		for (int j = 1; j < m; j++) {
			// translation
			char AA1 = '?', AA2 = '?';
			if (i > 2) AA1 = s1->TranslateNTtoAA(i-3);
			if (j > 2) AA2 = s2->TranslateNTtoAA(j-3);
			// stop-codon check
			int stopS1, stopS2;
			stopS1 = (AA1 == '*') ? stop_cost : 0;
			stopS2 = (AA2 == '*') ? stop_cost : 0;
			// AA substitution
			int subst_AA;
			if (AA1 == '*' || AA2 == '*') subst_AA = stopS1 + stopS2;
			else if (i - 3 >= 0 && j - 3 >= 0)
				subst_AA = aa_score_matrix[AA1*128 + AA2];
			else subst_AA = MAXINT;
			// calculating optimal pass
			// NT align
			int score = nt_score_matrix[s1->nt_seq[i-1]*128+s2->nt_seq[j-1]]
					+ F[(i-1)*size + j - 1];
			int way = 12;
			if (j - 2 >= 0 && score < F[(i-1)*size + j-2]) {
				score = F[(i-1)*size + j - 2];
				way = 13;
			}
			if (i - 2 >= 0 && score < F[(i-2)*size + j-1]) {
				score = F[(i-2)*size + j-1];
				way = 14;
			}
			if (i - 2 >= 0 && j-2 >= 0 && score < F[(i-2)*size + j-2]) {
				score = F[(i-2)*size + j-2];
				way = 15;
			}
			score += 2 * gap_frame;
			// AA align
			int tmp;
			if (i-3 >= 0) {
				tmp = F[(i-3)*size + j] + stopS1 + gap_extension;
				if (W[(i-3)*size + j] > 6) tmp += gap_open;
				if (score <= tmp) {
					score = tmp;
					way = 1;
				}
			}
			if (j-3 >= 0) {
				tmp = F[i*size + j-3] + stopS2 + gap_extension;
				if (W[i*size + j-3] > 6) tmp += gap_open;
				if (score <= tmp) {
					score = tmp;
					way = 2;
				}
			}
			if (i-3 >= 0 && j-3 >= 0 && score <= F[(i-3)*size+j-3] + subst_AA) {
				score = F[(i-3)*size+j-3] + subst_AA;
				way = 7;
			}
			// AA & NT
			if (i > 2 && j > 1 && score < F[(i-3)*size+j-2]+stopS1+gap_frame) {
				score = F[(i-3)*size+j-2]+stopS1+gap_frame;
				way = 8;
			}
			if (i > 2 && score < F[(i-3)*size + j-1] + stopS1 + gap_frame) {
				score = F[(i-3)*size + j-1] + stopS1 + gap_frame;
				way = 9;
			}
			if (i > 1 && j > 2 && score < F[(i-2)*size+j-3]+stopS2+gap_frame) {
				score = score < F[(i-2)*size+j-3]+stopS2+gap_frame;
				way = 10;
			}
			if (j > 2 && score < F[(i-1)*size + j-3] + stopS2 + gap_frame) {
				score = F[(i-1)*size+j-3]+stopS2+gap_frame;
				way = 11;
			}
			// NT gap
			tmp = F[i*size+j-1] + gap_extension + gap_frame;
			if (W[i*size+j-1] > 6) tmp += gap_open;
			if (score < tmp) {
				score = tmp;
				way = 3;
			}
			if (j > 1) {
				tmp = F[i*size+j-2] + gap_extension + gap_frame;
				if (W[i*size+j-2] > 6) tmp += gap_open;
				if (score < tmp) {
					score = tmp;
					way = 4;
				}
			}
			tmp = F[(i-1)*size+j] + gap_extension + gap_frame;
			if (W[(i-1)*size+j] > 6) tmp += gap_open;
			if (score < tmp) {
				score = tmp;
				way = 5;
			}
			if (i > 1) {
				tmp = F[(i-2)*size+j] + gap_extension + gap_frame;
				if (W[(i-2)*size+j] > 6) tmp += gap_open;
				if (score < tmp) {
					score = tmp;
					way = 6;
				}
			}
			// save the best variant
			F[i*size+j] = score;
			W[i*size+j] = way;\
		}
	}
	// searching the best score
	int result_score = F[n*size-1];
	int i = n-1, j = m-1;
	// ...search in last line
	for (int index = 0; index < m-1; index++) 
		if (F[(n-1)*size+index] > result_score) {
			j = index;
			result_score = F[(n-1)*size+index];
		}
	// ...search in last column
	for (int index = 0; index < n-1; index++) 
		if (F[index*size + m-1] > result_score) {
			i = index; 
			j = m-1;
			result_score = F[index*size + m-1];
		}
	// i, j - point of best score
	return F[i*size + j];
}