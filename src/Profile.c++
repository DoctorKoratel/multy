#include <string.h>

#include "Profile.h"
#include "PairwiseAlign.h"

#define MAXINT 0x7FFFFFFF

// debug flags
// compile program with -DDEBUG and -DPRINT_MATRIX to obtain debug info
//#define DEBUG
//#define PRINT_MATRIX

#ifdef TIME
	#include <ctime>
	#include <iostream>
	using std::cout;
	using std::endl;
#endif

#ifdef DEBUG
	#include <iostream>
	using std::cout;
	using std::endl;
#endif

// global tables for merge sequences
// free in UPGMAfree
extern int dim;
unsigned int* best_mvt = NULL;
float* scores = NULL;

// array of parameters 
// parameters[0] = gap open cost
// parameters[1] = gap extension cost
// parameters[2] = gap frame cost
// parameters[3] = stop cost
int parameters[5]; 
const int* nt_subst; // int array [128 * 128] nt_subst['A'*128+'C']<-score A & C
const int* aa_subst; // int array [128 * 128]
extern int bonus;

#define NUC 4
#define AMINO 23
const unsigned char nucleotides[] = {'A', 'C', 'G', 'T'};
const unsigned char amino[] = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 
			'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'X', 'Y', 'Z'};

void UPGMA(std::vector<BioSeq*>& sequences, func distance, 
						int go, int ge, int gf, int sc, const int* nt, const int* aa) 
{
	unsigned int n = sequences.size();
	if (!n) return; // zero check
	// memory allocation for global tables
	if (!best_mvt) best_mvt = new unsigned int [dim*dim];
	if (!scores) scores = new float [dim*dim];
	// fill parameters
	parameters[0] = go; // gap open
	parameters[1] = ge; // gap extension
	parameters[2] = gf; // gap frame
	parameters[3] = sc; // stop cost
	nt_subst = nt; // NT substitution matrix
	aa_subst = aa; // AA substitution matrix
	// initialization
	double* matrix = new double[n * n];
	double max = -((int)MAXINT);
	unsigned int max_i = 0, max_j = 0;
	
#ifdef TIME
	std::clock_t start;
	double duration;
	start = std::clock();
#endif
	
	// fill UPGMA table
	for (unsigned int i = 0; i < n; i++) {
		// calculating distances between sequences
		for (unsigned int j = 0; j < i; j++) {
			matrix[i*n + j] = distance(sequences[i],sequences[j]);
			// searching max score
			if (matrix[i*n + j] > max) {
				max = matrix[i*n + j];
				max_i = i;
				max_j = j;
			}
		}
	}
	
#ifdef TIME
	duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
	cout << "Building UPGMA matrix done! " << duration << "\n";
#endif

#ifdef DEBUG
	cout << "Building UPGMA matrix done!" << endl;
	for (unsigned int i = 0; i < n; i++) {
		for (unsigned int j = 0; j < i; j++) {
			cout << matrix[i*n + j] << '\t';
		}
		cout << '.' << endl;
	}
	cout << "Profile: (sequence #" << max_i << ')' << endl; 
#endif

	// create profile
	Profile profile;
	profile.sequences.push_back(sequences[max_i]);
	// add sequnces to profile
	int alive = n;
	bool* corpse = new bool[n];
	memset(corpse, false, sizeof(bool)*n); // everybody is alive
	while (alive > 1) {
		// add sequence max_j to profile

#ifdef DEBUG
	cout << "add sequence #" << max_j << endl;
#endif

		parameters[4] = bonus;  // bonus for not break the frame

#ifdef TIME
	start = std::clock();
#endif

		profile = profile + sequences[max_j];

#ifdef TIME
	duration = ( std::clock() - start ) / (double) CLOCKS_PER_SEC;
	cout << "add sequence #" << max_j << " " << duration << "\n";
#endif

		// killing max_j
		corpse[max_j] = true; 
		alive--;
		// updating distances && searching new max
		max = -MAXINT;
		int new_max_j = max_j;
		// ...in max_i column
		for (unsigned int k = max_i + 1; k < n; k++) {
			if (!corpse[k]) {
				if (k > max_j)
					matrix[k*n+max_i] += matrix[k*n+max_j] * matrix[max_j*(n+1)];
				else 
					matrix[k*n+max_i] += matrix[max_j*n+k] * matrix[max_j*(n+1)];
				if (max < matrix[k*n+max_i]) {
					max = matrix[k*n+max_i];
					new_max_j = k;
				}
			}
		}
		// ...in max_i row
		for (unsigned int k = 0; k < max_i; k++) {
			if (!corpse[k]) {
				if (k > max_j)
					matrix[max_i*n+k] += matrix[k*n+max_j] * matrix[max_j*(n+1)];
				else 
					matrix[max_i*n+k] += matrix[max_j*n+k] * matrix[max_j*(n+1)];
				if (max < matrix[max_i*n+k]) {
					max = matrix[max_i*n+k];
					new_max_j = k;
				}
			}
		}
		max_j = new_max_j;
	}
	// completion frames
	int l = sequences[0]->Length();
	int f = (3 - l % 3) % 3;
	for_each(sequences.begin(), sequences.end(), [l, f] (BioSeq* s) {
		s->nt_seq.insert(l, f, '-');
	});
	// deleting temporary matrices
	delete[] corpse;
	delete[] matrix;
	// comment next line if you want to use UPGMA again for other sequences
	// and call UPGMAfree in main program in the end
	UPGMAfree();
}

inline void UPGMAfree() { 
	if (best_mvt) delete[] best_mvt;
	if (scores) delete[] scores;
	best_mvt = NULL;
	scores = NULL;
	dim = 0;
}

Profile& Profile :: operator + (BioSeq* sequence) {
	// reallocate memory (if necessary)
	int new_size = (dim) ? dim : 1;
	int dim1 = sequences[0]->Length()+1, dim2 = sequence->Length()+1;
	while (dim1 > new_size || dim2 > new_size) {
		new_size *= 2;
	}
	if (new_size > dim) {
		UPGMAfree();
		dim = new_size;
		best_mvt = new unsigned int[dim*dim];
		scores = new float[dim*dim];
	}
	// initialization
	memset(best_mvt, 0xFF, sizeof(unsigned int)*dim*dim);
	memset(scores, 0, sizeof(float)*dim*dim);
	// calculating new profile
	//  * fill score matrix

#ifdef DEBUG
	cout << "profile length: " << dim1-1 << endl;
	cout << "sequence length: " << dim2-1 << endl;
	cout << "allocated matrix: " << dim << 'x' << dim << endl;
#endif

	for (int i = 1; i < dim1; i++) {
		for (int j = 1; j < dim2; j++) {
			/* 25 possible moves
				NT align (x10)
				================
				X | XX | -X | XX
				- | -X | XX | XX
				================
				- | XX | X- | X
				X | X- | XX | X
				================
				  | XX | -- |
				  | -- | XX |
			*/
			float score = ColumnNTscore(sequence, i-1, j-1) + scores[(i-1)*dim+j-1];
			int way = 15;
			if (best_mvt[i*dim+j-1] < 8 // gap extension
				&& score < scores[i*dim+j-1] + parameters[1]) {
				score = scores[i*dim+j-1] + parameters[1];
				way = 1;
			} else if (best_mvt[i*dim+j-1] > 7 // gap open
				&& score < scores[i*dim+j-1] + parameters[0]) {
				score = scores[i*dim+j-1] + parameters[0];
				way = 1;
			}
			if (best_mvt[(i-1)*dim+j] > 7 && best_mvt[(i-1)*dim+j] < 15 
				&& score < scores[(i-1)*dim+j] + parameters[1]) {
				score = scores[(i-1)*dim+j] + parameters[1];
				way = 8;
			} else if ((best_mvt[(i-1)*dim+j] > 14 || best_mvt[(i-1)*dim+j] < 8)
				&& score < scores[(i-1)*dim+j] + parameters[0]) {
				score = scores[(i-1)*dim+j] + parameters[0];
				way = 8;
			}
			if (i - 2 >= 0 && j - 2 >= 0) {
				if (score < ColumnNTscore(sequence, i-1, j-1) 
					+ ColumnNTscore(sequence, i-2, j-2) + scores[(i-2)*dim+j-2]) {
					score = ColumnNTscore(sequence, i-1, j-1) 
						+ ColumnNTscore(sequence, i-2, j-2) + scores[(i-2)*dim+j-2];
					way = 16;
				}
			}
			if (i - 2 >= 0) {
				if (score < scores[(i-2)*dim+j-1] + ColumnNTscore(sequence, i-1, j-1)) {
					score = scores[(i-2)*dim+j-1] + ColumnNTscore(sequence, i-1, j-1);
					way = 17;
				} 
				if (score < scores[(i-2)*dim+j-1] + ColumnNTscore(sequence, i-2, j-1)) {
					score = scores[(i-2)*dim+j-1] + ColumnNTscore(sequence, i-2, j-1);
					way = 2;
				}
				if (best_mvt[(i-2)*dim+j] < 8 
					&& score < scores[(i-2)*dim+j] + parameters[1]) {
					score = scores[(i-2)*dim+j] + parameters[1];
					way = 3;
				} else if (best_mvt[(i-2)*dim+j] > 7 
					&& score < scores[(i-2)*dim+j] + parameters[0]) {
					score = scores[(i-2)*dim+j] + parameters[0];
					way = 3;
				}
			}
			if (j - 2 >= 0) {
				if (score < scores[(i-1)*dim+j-2] + ColumnNTscore(sequence, i-1, j-1)) {
					score = scores[(i-1)*dim+j-2] + ColumnNTscore(sequence, i-1, j-1);
					way = 18;
				} 
				if (score < scores[(i-1)*dim+j-2] + ColumnNTscore(sequence, i-1, j-2)) {
					score = scores[(i-1)*dim+j-2] + ColumnNTscore(sequence, i-1, j-2);
					way = 9;
				}
				if (best_mvt[i*dim+j-2] > 7 && best_mvt[i*dim+j-2] < 15
					&& score < scores[i*dim+j-2] + parameters[1]) {
					score = scores[i*dim+j-2] + parameters[1];
					way = 10;
				} else if ((best_mvt[i*dim+j-2] > 14 || best_mvt[i*dim+j-2] < 8)
					&& score < scores[i*dim+j-2] + parameters[0]) {
					score = scores[i*dim+j-2] + parameters[0];
					way = 10;
				}
			}
			if (j > 2) score += parameters[2]; // and one more addition after AA check
			// AA align
			// with frame shift
			if (i - 3 >= 0 && j - 2 >= 0) {
				if (score < scores[(i-3)*dim+j-2] + ColumnNTscore(sequence, i-2, j-2) 
					+ ColumnNTscore(sequence, i-1, j-1)) {
					score = scores[(i-3)*dim+j-2] + ColumnNTscore(sequence, i-2, j-2) 
					+ ColumnNTscore(sequence, i-1, j-1);
					way = 20;
				}
				if (score < scores[(i-3)*dim+j-2] + ColumnNTscore(sequence,i-3, j-2) 
					+ ColumnNTscore(sequence,i-1, j-1)) {
					score = scores[(i-3)*dim+j-2] + ColumnNTscore(sequence,i-3, j-2) 
					+ ColumnNTscore(sequence,i-1, j-1);
					way = 22;
				}
				if (score < scores[(i-3)*dim+j-2] + ColumnNTscore(sequence,i-3, j-2) 
					+ ColumnNTscore(sequence,i-2, j-1)) {
					score = scores[(i-3)*dim+j-2] + ColumnNTscore(sequence,i-3, j-2) 
					+ ColumnNTscore(sequence,i-2, j-1);
					way = 4;
				}
			}
			if (i - 3 >= 0) {
				if (score < scores[(i-3)*dim+j-1] + ColumnNTscore(sequence,i-1, j-1)) {
					score = scores[(i-3)*dim+j-1] + ColumnNTscore(sequence,i-1, j-1);
					way = 23;
				}
				if (score < scores[(i-3)*dim+j-1] + ColumnNTscore(sequence,i-2, j-1)) {
					score = scores[(i-3)*dim+j-1] + ColumnNTscore(sequence,i-2, j-1);
					way = 5;
				}
				if (score < scores[(i-3)*dim+j-1] + ColumnNTscore(sequence,i-3, j-1)) {
					score = scores[(i-3)*dim+j-1] + ColumnNTscore(sequence,i-3, j-1);
					way = 6;
				}
			}
			if (i - 2 >= 0 && j - 3 >= 0) {
				if (score < scores[(i-2)*dim+j-3] + ColumnNTscore(sequence,i-2, j-2) 
					+ ColumnNTscore(sequence,i-1, j-1)) {
					score = scores[(i-2)*dim+j-3] + ColumnNTscore(sequence,i-2, j-2) 
					+ ColumnNTscore(sequence,i-1, j-1);
					way = 21;
				}
				if (score < scores[(i-2)*dim+j-3] + ColumnNTscore(sequence,i-2, j-3) 
					+ ColumnNTscore(sequence,i-1, j-1)) {
					score = scores[(i-2)*dim+j-3] + ColumnNTscore(sequence,i-2, j-3) 
					+ ColumnNTscore(sequence,i-1, j-1);
					way = 24;
				}
				if (score < scores[(i-2)*dim+j-3] + ColumnNTscore(sequence,i-2, j-3) 
					+ ColumnNTscore(sequence,i-1, j-2)) {
					score = scores[(i-2)*dim+j-3] + ColumnNTscore(sequence,i-2, j-3) 
					+ ColumnNTscore(sequence,i-2, j-1);
					way = 11;
				}
			}
			if (j - 3 >= 0) {
				if (score < scores[(i-1)*dim+j-3] + ColumnNTscore(sequence,i-1, j-1)) {
					score = scores[(i-1)*dim+j-3] + ColumnNTscore(sequence,i-1, j-1);
					way = 25;
				}
				if (score < scores[(i-1)*dim+j-3] + ColumnNTscore(sequence,i-1, j-2)) {
					score = scores[(i-1)*dim+j-3] + ColumnNTscore(sequence,i-1, j-2);
					way = 12;
				}
				if (score < scores[(i-1)*dim+j-3] + ColumnNTscore(sequence,i-1, j-3)) {
					score = scores[(i-1)*dim+j-3] + ColumnNTscore(sequence,i-1, j-3);
					way = 13;
				}
			}
			score += parameters[2];
			// AA match
			if (i - 3 >= 0) {
				if (best_mvt[(i-3)*dim+j-1] < 8) {
					if (score < scores[(i-3)*dim+j] + parameters[1]) {
						score = scores[(i-3)*dim+j] + parameters[1];
						way = 7;
					}
				} else {
					if (score < scores[(i-3)*dim+j] + parameters[0]) {
						score = scores[(i-3)*dim+j] + parameters[0];
						way = 7;
					}
				}
			}
			if (j - 3 >= 0) {
				if (best_mvt[i*dim+j-3] > 7 && best_mvt[i*dim+j-3] < 15) {
					if (score < scores[i*dim+j-3] + parameters[1]) {
						score = scores[i*dim+j-3] + parameters[1];
						way = 14;
					}
				} else {
					if (score < scores[i*dim+j-3] + parameters[0]) {
						score = scores[i*dim+j-3] + parameters[0];
						way = 14;
					}
				}
			}
			if (i - 3 >= 0 && j - 3 >= 0) { 
				if (best_mvt[(i-3)*dim+j-3] == 19) {
					if (score < ColumnAAscore(sequence, i-3,j-3) 
						+ ColumnNTscore(sequence, i-3,j-3)+ColumnNTscore(sequence, i-2,j-2) 
						+ ColumnNTscore(sequence, i-1, j-1) + scores[(i-3)*dim+j-3] 
						+ parameters[4]) {
						score = ColumnAAscore(sequence, i-3,j-3) 
						+ ColumnNTscore(sequence, i-3,j-3)+ColumnNTscore(sequence, i-2,j-2) 
						+ ColumnNTscore(sequence, i-1, j-1) + scores[(i-3)*dim+j-3] 
						+ parameters[4];
						way = 19;
						parameters[4]++;
						if (parameters[4] > bonus) parameters[4] = bonus;
					} 
				} else {
					if (score < ColumnAAscore(sequence, i-3,j-3) 
						+ ColumnNTscore(sequence, i-3,j-3)+ColumnNTscore(sequence, i-2,j-2) 
						+ ColumnNTscore(sequence, i-1, j-1) + scores[(i-3)*dim+j-3]) {
						score = ColumnAAscore(sequence, i-3,j-3) 
						+ ColumnNTscore(sequence, i-3,j-3)+ColumnNTscore(sequence, i-2,j-2) 
						+ ColumnNTscore(sequence, i-1, j-1)  + scores[(i-3)*dim+j-3];
						way = 19;
						parameters[4] = bonus;
					}
				}
			}
			// score & way contain best option
			scores[i*dim+j] = score;
			best_mvt[i*dim+j] = way;
		}
	}
	//  * update sequences
	// searching the best score
	float result_score = scores[(dim1-1)*dim+dim2-1];
	int i = dim1-1, j = dim2-1;
	// ...search in last line
	for (int index = dim2 - dim2 / 10; index < dim2-1; index++) 
		if (scores[(dim1-1)*dim+index] > result_score) {
			j = index;
			result_score = scores[(dim1-1)*dim+index];
		}
	// ...search in last column
	for (int index = dim1 - dim1 / 10; index < dim1-1; index++) 
		if (scores[index*dim + dim2-1] > result_score) {
			i = index; 
			j = dim2-1;
			result_score = scores[index*dim + dim2-1];
		}
	// i, j - point of best score
	// updating sequences
	this->InsertGap(dim1-1, dim2-1-j); // insert gaps in Profile 
	sequence->InsertGap(dim2-1, dim1-1-i); // insert gaps in sequence
	
#ifdef PRINT_MATRIX
	for (int j = 0; j < dim2; j++) { // print indexes
		cout << j << '\t';
	}
	cout << endl;
	for (int i = 0; i < dim1; i++) {
		for (int j = 0; j < dim2; j++) {
			cout << best_mvt[i*dim+j] << '\t';
		}
		cout << endl;
	}
	cout << endl << endl;
	for (int i = 0; i < dim1; i++) {
		for (int j = 0; j < dim2; j++) {
			cout << scores[i*dim+j] << '\t';
		}
		cout << endl;
	}
	cout << endl << "(max) " << i << " - " << j << ' ' << result_score << endl;
#endif 
	
	switch (best_mvt[i*dim + j]) {
		case 1:
			sequence->InsertGap(j); 
			i--;
			break;
		case 2:
			sequence->InsertGap(j);
			i -= 2;
			j--;
			break;
		case 3:
			sequence->InsertGap(j, 2);
			i -= 2;
			break;
		case 8:
			InsertGap(i);
			j--;
			break;
		case 9:
			InsertGap(i);
			i--;
			j -= 2;
			break;
		case 10:
			InsertGap(i, 2);
			j -= 2;
			break;
		case 15:
			i--;
			j--;
			break;
		case 16:
			i -= 2;
			j -= 2;
			break;
		case 17:
			i -= 2;
			j--;
			sequence->InsertGap(j);
			break;
		case 18:
			i--;
			j -= 2;
			InsertGap(i);
			break;
	}
	
	while (i && j) {
		switch (best_mvt[i*dim + j]) {
			case 1:
				sequence->InsertGap(j, 3); 
				i--;
				// save frame
				InsertGap(i, 2);
				break;
			case 2:
				sequence->InsertGap(j);
				i -= 2;
				j--;
				// save frame
				InsertGap(i);
				sequence->InsertGap(j);
				break;
			case 3:
				sequence->InsertGap(j, 3);
				i -= 2;
				// save frame
				InsertGap(i);
				break;
			case 4:
				sequence->InsertGap(j);
				i -= 3;
				j -= 2;
				break;
			case 5:
				sequence->InsertGap(j);
				j--;
				sequence->InsertGap(j);
				i -= 3;
				break;
			case 6:
				sequence->InsertGap(j, 2);
				i -= 3;
				j--;
				break;
			case 7:
				sequence->InsertGap(j, 3);
				i -= 3;
				break;
			case 8:
				InsertGap(i, 3);
				j--;
				// save frame
				sequence->InsertGap(j, 2);
				break;
			case 9:
				InsertGap(i);
				i--;
				j -= 2;
				// save frame
				InsertGap(i);
				sequence->InsertGap(j);
				break;
			case 10:
				InsertGap(i, 3);
				j -= 2;
				// save frame
				sequence->InsertGap(j);
				break;
			case 11:
				InsertGap(i);
				i -= 2;
				j -= 3;
				break;
			case 12:
				InsertGap(i);
				i--;
				InsertGap(i);
				j -= 3;
				break;
			case 13:
				InsertGap(i, 2);
				i--;
				j -= 3;
				break;
			case 14:
				InsertGap(i, 3);
				j -= 3;
				break;
			case 15:
				i--;
				j--;
				// save frame
				InsertGap(i, 2);
				sequence->InsertGap(j, 2);
				break;
			case 16:
				i -= 2;
				j -= 2;
				// save frame
				InsertGap(i);
				sequence->InsertGap(j);
				break;
			case 17:
				i -= 2;
				j--;
				// save frame
				InsertGap(i);
				sequence->InsertGap(j, 2);
				break;
			case 18:
				i--;
				j -= 2;
				// save frame
				InsertGap(i, 2);
				sequence->InsertGap(j);
				break;
			case 19:
				i -= 3;
				j -= 3;
				break;
			case 20:
				i -= 3;
				j -= 2;
				sequence->InsertGap(j);
				break;
			case 21:
				i -= 2;
				InsertGap(i);
				j -= 3;
				break;
			case 22:
				i -= 3;
				j--;
				sequence->InsertGap(j);
				j--;
				break;
			case 23:
				i -= 3;
				j--;
				sequence->InsertGap(j, 2);
				break;
			case 24:
				i--;
				InsertGap(i);
				i--;
				j -= 3;
				break;
			case 25:
				i--;
				InsertGap(i, 2);
				j -= 3;
				break;
			default:
				// some problems
				// return NULL; // DEBUG
				return *this;
		}
	}
	InsertGap(0, j);
	sequence->InsertGap(0, i);
	// add sequence to Profile
	sequences.push_back(sequence);
	
#ifdef DEBUG
	for_each(sequences.begin(), sequences.end(), [](BioSeq* s) {
		s->PrintNT(cout);
	});
	cout << endl;
#endif
	
	// return new Profile
	return *this;
}

Profile& Profile :: operator + (Profile& another) {
	// zero check
	if (another.sequences.size() == 0)
		return *this;
	if (sequences.size() == 0) {
		sequences = another.sequences;
		return *this;
	}
	// reallocate memory (if necessary)
	int new_size = (dim) ? dim : 1;
	int dim1 = sequences[0]->Length()+1, dim2 = another.sequences[0]->Length()+1;
	while (dim1 > new_size || dim2 > new_size) {
		new_size *= 2;
	}
	if (new_size > dim) {
		UPGMAfree();
		dim = new_size;
		best_mvt = new unsigned int[dim*dim];
		scores = new float[dim*dim];
	}
	// initialization
	memset(best_mvt, 0xFF, sizeof(unsigned int)*dim*dim);
	memset(scores, 0, sizeof(float)*dim*dim);
	// calculating new profile
	//  * fill score matrix
	
#ifdef DEBUG
	cout << "seq1 size: " << dim1-1 << endl;
	cout << "seq2 size: " << dim2-1 << endl;
	cout << "allocated matrix: " << dim << 'x' << dim << endl;
#endif
	
	for (int i = 1; i < dim1; i++) {
		for (int j = 1; j < dim2; j++) {
			/* 25 possible moves
				NT align (x10)
				================
				X | XX | -X | XX
				- | -X | XX | XX
				================
				- | XX | X- | X
				X | X- | XX | X
				================
				  | XX | -- |
				  | -- | XX |
			*/
			float score = ColumnNTscore(another, i-1, j-1) + scores[(i-1)*dim+j-1];
			int way = 15;
			if (j - 1 >= 0) {
				if (best_mvt[i*dim+j-1] < 8 // gap extension
					&& score < scores[i*dim+j-1] + parameters[1]) {
					score = scores[i*dim+j-1] + parameters[1];
					way = 1;
				} else if (best_mvt[i*dim+j-1] > 7 // gap open
					&& score < scores[i*dim+j-1] + parameters[0]) {
					score = scores[i*dim+j-1] + parameters[0];
					way = 1;
				}
			}
			if (i - 1 >= 0) {
				if (best_mvt[(i-1)*dim+j] > 7 && best_mvt[(i-1)*dim+j] < 15 
					&& score < scores[(i-1)*dim+j] + parameters[1]) {
					score = scores[(i-1)*dim+j] + parameters[1];
					way = 8;
				} else if ((best_mvt[(i-1)*dim+j] > 14 || best_mvt[(i-1)*dim+j] < 8)
					&& score < scores[(i-1)*dim+j] + parameters[0]) {
					score = scores[(i-1)*dim+j] + parameters[0];
					way = 8;
				}
			}
			if (i - 2 >= 0 && j - 2 >= 0) {
				if (score < ColumnNTscore(another, i-1, j-1) 
					+ ColumnNTscore(another, i-2, j-2) + scores[(i-2)*dim+j-2]) {
					score = ColumnNTscore(another, i-1, j-1) 
						+ ColumnNTscore(another, i-2, j-2) + scores[(i-2)*dim+j-2];
					way = 16;
				}
			}
			if (i - 2 >= 0 && j - 1 >= 0) {
				if (best_mvt[(i-2)*dim+j-1] < 8
					&& score < scores[(i-2)*dim+j-1] + parameters[1] 
					+ ColumnNTscore(another, i-1, j-1)) {
					score = scores[(i-2)*dim+j-1] + parameters[1] 
						+ ColumnNTscore(another, i-1, j-1);
					way = 17;
				} else if (best_mvt[(i-2)*dim+j-1] > 7
					&& score < scores[(i-2)*dim+j-1] + parameters[0] 
					+ ColumnNTscore(another, i-1, j-1)) {
					score = scores[(i-2)*dim+j-1] + parameters[0] 
						+ ColumnNTscore(another, i-1, j-1);
					way = 17;
				} 
				if (score < scores[(i-2)*dim+j-1] + ColumnNTscore(another, i-2, j-1) 
					+ parameters[0]) {
					score = scores[(i-2)*dim+j-1] + ColumnNTscore(another, i-2, j-1) 
						+ parameters[0];
					way = 2;
				}
			}
			if (i - 1 >= 0 && j - 2 >= 0) {
				if (best_mvt[(i-1)*dim+j-2] > 7 && best_mvt[(i-1)*dim+j-2] < 15
					&& score < scores[(i-1)*dim+j-2] + parameters[1] 
					+ ColumnNTscore(another, i-1, j-1)) {
					score = scores[(i-1)*dim+j-2] + parameters[1] 
						+ ColumnNTscore(another, i-1, j-1);
					way = 18;
				} else if ((best_mvt[(i-1)*dim+j-2] > 14 || best_mvt[(i-1)*dim+j-2] < 8)
					&& score < scores[(i-1)*dim+j-2] + parameters[0] 
					+ ColumnNTscore(another, i-1, j-1)) {
					score = scores[(i-1)*dim+j-2] + parameters[0] 
						+ ColumnNTscore(another, i-1, j-1);
					way = 18;
				} 
				if (score < scores[(i-1)*dim+j-2] + ColumnNTscore(another, i-1, j-2) 
					+ parameters[0]) {
					score = scores[(i-1)*dim+j-2] + ColumnNTscore(another, i-1, j-2) 
						+ parameters[0];
					way = 9;
				}
			}
			if (i - 2 >= 0) {
				if (best_mvt[(i-2)*dim+j] < 8 
					&& score < scores[(i-2)*dim+j] + parameters[1] * 2) {
					score = scores[(i-2)*dim+j] + parameters[1] * 2;
					way = 3;
				} else if (best_mvt[(i-2)*dim+j] > 7 
					&& score < scores[(i-2)*dim+j] + parameters[0] + parameters[1]) {
					score = scores[(i-2)*dim+j] + parameters[0] + parameters[1];
					way = 3;
				}
			}
			if (j - 2 >= 0) {
				if (best_mvt[i*dim+j-2] > 7 && best_mvt[i*dim+j-2] < 15
					&& score < scores[i*dim+j-2] + parameters[1] * 2) {
					score = scores[i*dim+j-2] + parameters[1] * 2;
					way = 10;
				} else if ((best_mvt[i*dim+j-2] > 14 || best_mvt[i*dim+j-2] < 8)
					&& score < scores[i*dim+j-2] + parameters[0] + parameters[1]) {
					score = scores[i*dim+j-2] + parameters[0] + parameters[1];
					way = 10;
				} 
			}
			score += parameters[2]; // and one more addition after AA check
			// AA align
			// with frame shift
			if (i - 3 >= 0 && j - 2 >= 0) {
				if (best_mvt[(i-3)*dim+j-2] < 8) {
					if (score < scores[(i-3)*dim+j-2] + ColumnNTscore(another, i-2, j-2) 
						+ ColumnNTscore(another, i-1, j-1) + parameters[1]) {
						score = scores[(i-3)*dim+j-2] + ColumnNTscore(another, i-2, j-2) 
						+ ColumnNTscore(another, i-1, j-1) + parameters[1];
						way = 20;
					}
				} else {
					if (score < scores[(i-3)*dim+j-2] + ColumnNTscore(another, i-2, j-2) 
						+ ColumnNTscore(another, i-1, j-1) + parameters[0]) {
						score = scores[(i-3)*dim+j-2] + ColumnNTscore(another, i-2, j-2) 
						+ ColumnNTscore(another, i-1, j-1) + parameters[0];
						way = 20;
					}
				}
				if (score < scores[(i-3)*dim+j-2] + ColumnNTscore(another,i-3, j-2) 
					+ ColumnNTscore(another,i-1, j-1) + parameters[0]) {
					score = scores[(i-3)*dim+j-2] + ColumnNTscore(another,i-3, j-2) 
					+ ColumnNTscore(another,i-1, j-1) + parameters[0];
					way = 22;
				}
				if (score < scores[(i-3)*dim+j-2] + ColumnNTscore(another,i-3, j-2) 
					+ ColumnNTscore(another,i-2, j-1) + parameters[0]) {
					score = scores[(i-3)*dim+j-2] + ColumnNTscore(another,i-3, j-2) 
					+ ColumnNTscore(another,i-2, j-1) + parameters[0];
					way = 4;
				}
			}
			if (i - 3 >= 0 && j - 1 >= 0) {
				int prev_gap = (best_mvt[(i-3)*dim+j-1] < 8) ? 1 : 0;
				if (score < scores[(i-3)*dim+j-1] + ColumnNTscore(another,i-1, j-1)
					+ parameters[prev_gap] + parameters[1]) {
					score = scores[(i-3)*dim+j-1] + ColumnNTscore(another,i-1, j-1)
					+ parameters[prev_gap] + parameters[1];
					way = 23;
				}
				if (score < scores[(i-3)*dim+j-1] + ColumnNTscore(another,i-2, j-1)
					+ parameters[prev_gap] + parameters[0]) {
					score = scores[(i-3)*dim+j-1] + ColumnNTscore(another,i-2, j-1)
					+ parameters[prev_gap] + parameters[0];
					way = 5;
				}
				if (score < scores[(i-3)*dim+j-1] + ColumnNTscore(another,i-3, j-1)
					+ parameters[0] + parameters[1]) {
					score = scores[(i-3)*dim+j-1] + ColumnNTscore(another,i-3, j-1)
					+ parameters[0] + parameters[1];
					way = 6;
				}
			}
			if (i - 2 >= 0 && j - 3 >= 0) {
				if (best_mvt[(i-2)*dim+j-3] > 7 && best_mvt[(i-2)*dim+j-3] < 15) {
					if (score < scores[(i-2)*dim+j-3] + ColumnNTscore(another,i-2, j-2) 
						+ ColumnNTscore(another,i-1, j-1) + parameters[1]) {
						score = scores[(i-2)*dim+j-3] + ColumnNTscore(another,i-2, j-2) 
						+ ColumnNTscore(another,i-1, j-1) + parameters[1];
						way = 21;
					}
				} else {
					if (score < scores[(i-2)*dim+j-3] + ColumnNTscore(another,i-2, j-2) 
						+ ColumnNTscore(another,i-1, j-1) + parameters[0]) {
						score = scores[(i-2)*dim+j-3] + ColumnNTscore(another,i-2, j-2) 
						+ ColumnNTscore(another,i-1, j-1) + parameters[0];
						way = 21;
					}
				}
				if (score < scores[(i-2)*dim+j-3] + ColumnNTscore(another,i-2, j-3) 
					+ ColumnNTscore(another,i-1, j-1) + parameters[0]) {
					score = scores[(i-2)*dim+j-3] + ColumnNTscore(another,i-2, j-3) 
					+ ColumnNTscore(another,i-1, j-1) + parameters[0];
					way = 24;
				}
				if (score < scores[(i-2)*dim+j-3] + ColumnNTscore(another,i-2, j-3) 
					+ ColumnNTscore(another,i-1, j-2) + parameters[0]) {
					score = scores[(i-2)*dim+j-3] + ColumnNTscore(another,i-2, j-3) 
					+ ColumnNTscore(another,i-2, j-1) + parameters[0];
					way = 11;
				}
			}
			if (i - 1 >= 0 && j - 3 >= 0) {
				int prev_gap = (best_mvt[(i-1)*dim+j-3] > 7 
					&& best_mvt[(i-1)*dim+j-3] < 15) ? 1 : 0;
				if (score < scores[(i-1)*dim+j-3] + ColumnNTscore(another,i-1, j-1)
					+ parameters[prev_gap] + parameters[1]) {
					score = scores[(i-1)*dim+j-3] + ColumnNTscore(another,i-1, j-1)
					+ parameters[prev_gap] + parameters[1];
					way = 25;
				}
				if (score < scores[(i-1)*dim+j-3] + ColumnNTscore(another,i-1, j-2)
					+ parameters[prev_gap] + parameters[0]) {
					score = scores[(i-1)*dim+j-3] + ColumnNTscore(another,i-1, j-2)
					+ parameters[prev_gap] + parameters[0];
					way = 12;
				}
				if (score < scores[(i-1)*dim+j-3] + ColumnNTscore(another,i-1, j-3)
					+ parameters[0] + parameters[1]) {
					score = scores[(i-1)*dim+j-3] + ColumnNTscore(another,i-1, j-3)
					+ parameters[0] + parameters[1];
					way = 13;
				}
			}
			if (i > 2 && j > 2)
				score += parameters[2];
			// AA match
			if (i - 3 >= 0) {
				if (best_mvt[(i-3)*dim+j-1] < 8) {
					if (score < scores[(i-3)*dim+j] + parameters[1] * 3) {
						score = scores[(i-3)*dim+j] + parameters[1] * 3;
						way = 7;
					}
				} else {
					if (score < scores[(i-3)*dim+j] + parameters[1] * 2 + parameters[0]) {
						score = scores[(i-3)*dim+j] + parameters[1] * 2 + parameters[0];
						way = 7;
					}
				}
			}
			if (j - 3 >= 0) {
				if (best_mvt[i*dim+j-3] > 7 && best_mvt[i*dim+j-3] < 15) {
					if (score < scores[i*dim+j-3] + parameters[1] * 3) {
						score = scores[i*dim+j-3] + parameters[1] * 3;
						way = 14;
					}
				} else {
					if (score < scores[i*dim+j-3] + parameters[1] * 2 + parameters[0]) {
						score = scores[i*dim+j-3] + parameters[1] * 2 + parameters[0];
						way = 14;
					}
				}
			}
			if (i - 3 >= 0 && j - 3 >= 0) { 
				if (score < ColumnAAscore(another, i-3,j-3) 
					+ ColumnNTscore(another, i-3,j-3) + ColumnNTscore(another, i-2,j-2) 
					+ ColumnNTscore(another, i-1, j-1) + scores[(i-3)*dim+j-3]) {
					score = ColumnAAscore(another, i-3,j-3) 
					+ ColumnNTscore(another, i-3,j-3) + ColumnNTscore(another, i-2,j-2) 
					+ ColumnNTscore(another, i-1, j-1)  + scores[(i-3)*dim+j-3];
					way = 19;
				}
			}
			// score & way contain best option
			scores[i*dim+j] = score;
			best_mvt[i*dim+j] = way;
		}
	}
	//  * update sequences
	// searching the best score
	float result_score = scores[(dim1-1)*dim+dim2-1];
	int i = dim1-1, j = dim2-1;
	// ...search in last line
	for (int index = 0; index < dim2-1; index++) 
		if (scores[(dim1-1)*dim+index] > result_score) {
			j = index;
			result_score = scores[(dim1-1)*dim+index];
		}
	// ...search in last column
	for (int index = 0; index < dim1-1; index++) 
		if (scores[index*dim + dim2-1] > result_score) {
			i = index; 
			j = dim2-1;
			result_score = scores[index*dim + dim2-1];
		}
	// i, j - point of best score
	// updating sequences
	this->InsertGap(dim1-1, dim2-1-j); // insert gaps in Profile 1
	another.InsertGap(dim2-1, dim1-1-i); // insert gaps in Profile 2
	
#ifdef PRINT_MATRIX
	for (int j = 0; j < dim2; j++) { // print indexes
		cout << j << '\t';
	}
	cout << endl;
	for (int i = 0; i < dim1; i++) {
		for (int j = 0; j < dim2; j++) {
			cout << best_mvt[i*dim+j] << '\t';
		}
		cout << endl;
	}
	cout << endl << endl;
	for (int i = 0; i < dim1; i++) {
		for (int j = 0; j < dim2; j++) {
			cout << scores[i*dim+j] << '\t';
		}
		cout << endl;
	}
	cout << endl << "(max) " << i << " - " << j << ' ' << result_score << endl;
#endif 
	
	while (i && j) {
		switch (best_mvt[i*dim + j]) {
			case 1:
				another.InsertGap(j, 3); 
				i--;
				// save frame
				InsertGap(i, 2);
				break;
			case 2:
				another.InsertGap(j);
				i -= 2;
				j--;
				// save frame
				InsertGap(i);
				another.InsertGap(j);
				break;
			case 3:
				another.InsertGap(j, 3);
				i -= 2;
				// save frame
				InsertGap(i);
				break;
			case 4:
				another.InsertGap(j);
				i -= 3;
				j -= 2;
				break;
			case 5:
				another.InsertGap(j);
				j--;
				another.InsertGap(j);
				i -= 3;
				break;
			case 6:
				another.InsertGap(j, 2);
				i -= 3;
				j--;
				break;
			case 7:
				another.InsertGap(j, 3);
				i -= 3;
				break;
			case 8:
				InsertGap(i, 3);
				j--;
				// save frame
				another.InsertGap(j, 2);
				break;
			case 9:
				InsertGap(i);
				i--;
				j -= 2;
				// save frame
				InsertGap(i);
				another.InsertGap(j);
				break;
			case 10:
				InsertGap(i, 3);
				j -= 2;
				// save frame
				another.InsertGap(j);
				break;
			case 11:
				InsertGap(i);
				i -= 2;
				j -= 3;
				break;
			case 12:
				InsertGap(i);
				i--;
				InsertGap(i);
				j -= 3;
				break;
			case 13:
				InsertGap(i, 2);
				i--;
				j -= 3;
				break;
			case 14:
				InsertGap(i, 3);
				j -= 3;
				break;
			case 15:
				i--;
				j--;
				// save frame
				InsertGap(i, 2);
				another.InsertGap(j, 2);
				break;
			case 16:
				i -= 2;
				j -= 2;
				// save frame
				InsertGap(i);
				another.InsertGap(j);
				break;
			case 17:
				i -= 2;
				j--;
				// save frame
				InsertGap(i);
				another.InsertGap(j, 2);
				break;
			case 18:
				i--;
				j -= 2;
				// save frame
				InsertGap(i, 2);
				another.InsertGap(j);
				break;
			case 19:
				i -= 3;
				j -= 3;
				break;
			case 20:
				i -= 3;
				j -= 2;
				another.InsertGap(j);
				break;
			case 21:
				i -= 2;
				InsertGap(i);
				j -= 3;
				break;
			case 22:
				i -= 3;
				j--;
				another.InsertGap(j);
				j--;
				break;
			case 23:
				i -= 3;
				j--;
				another.InsertGap(j, 2);
				break;
			case 24:
				i--;
				InsertGap(i);
				i--;
				j -= 3;
				break;
			case 25:
				i--;
				InsertGap(i, 2);
				j -= 3;
				break;
			default:
				// some problems
				// return NULL; // DEBUG
				return *this;
		}
	}
	InsertGap(0, j);
	another.InsertGap(0, i);
	// add sequences from Profile2 to Profile1
	for_each(another.sequences.begin(), another.sequences.end(), [&](BioSeq* s) {
		sequences.push_back(s);
	});
	
#ifdef DEBUG
	for_each(sequences.begin(), sequences.end(), [](BioSeq* s) {
		s->PrintNT(cout);
	});
	cout << endl;
#endif
	
	// return new Profile
	return *this;
}

float Profile :: ColumnNTscore(BioSeq* seq, int index1, int index2) {
	int score = 0;
	CalcFrequenciesNT(index1);
	for (int i = 0; i < NUC; i++) {
		unsigned char char1 = nucleotides[i];
		unsigned char char2 = (*seq)[index2];
		if (frequency[char1]) {
			score += nt_subst[char1*128+char2] * (frequency[char1] + 1); 
		}
	}
	return score / (sequences.size()+1.);
}

float Profile :: ColumnAAscore(BioSeq* seq, int index1, int index2) {
	int score = 0;
	CalcFrequenciesAA(index1);
	for (int i = 0; i < AMINO; i++) {
		unsigned char char1 = amino[i];
		unsigned char char2 = seq->TranslateNTtoAA(index2);
		if (frequency[char1]) {
			score += aa_subst[char1*128+char2] * (frequency[char1] + 1); 
		}
	}
	score += frequency['-'] * parameters[1]; // gap penalty
	score += frequency['!'] * parameters[2]; // gap frame
	if (index1+3 < sequences[0]->Length() && index2+3 < seq->Length())
		score += frequency['*'] * parameters[3]; // stop cost
	return score / (sequences.size()+1.);
}

float Profile :: ColumnNTscore(Profile& another, int index1, int index2) {
	int score = 0, steps = 0;
	float denominator = sequences.size() + another.sequences.size();
	CalcFrequenciesNT(index1); 
	another.CalcFrequenciesNT(index2);
	for (int i = 0; i < NUC; i++) {
		frequency[nucleotides[i]] += another.frequency[nucleotides[i]];
	}
	// frequency - array for index1+index2 columns
	for (int i = 0; i < NUC; i++) {
		unsigned char char1 = nucleotides[i];
		if (!frequency[char1]) continue;
		steps++;
		score += frequency[char1] * nt_subst[char1*128 + char1];
		for (int j = i+1; j < NUC; j++) {
			unsigned char char2 = nucleotides[j];
			if (!frequency[char2]) continue;
			score += (frequency[char1]+frequency[char2]) * nt_subst[char1*128+char2];
			steps++;
		}
	}
	if (frequency['-']) {
		score += frequency['-'] * parameters[1]; // gap penalty
		steps++;
	}
	return score / denominator / steps;
}

float Profile :: ColumnAAscore(Profile& another, int index1, int index2) {
	int score = 0, steps = 0;
	float denominator = sequences.size() + another.sequences.size();
 	CalcFrequenciesAA(index1); 
 	another.CalcFrequenciesAA(index2);
	for (int i = 0; i < AMINO; i++) {
		frequency[amino[i]] += another.frequency[amino[i]];
	}
	// frequency - array for index1+index2 columns
	for (int i = 0; i < AMINO; i++) {
		unsigned char char1 = amino[i];
		if (!frequency[char1]) continue;
		steps++;
		score += frequency[char1] * aa_subst[char1*128 + char1];
		for (int j = i+1; j < AMINO; j++) {
			unsigned char char2 = amino[j];
			if (!frequency[char2]) continue;
			score += (frequency[char1]+frequency[char2]) * aa_subst[char1*128+char2];
			steps++;
		}
	}
	if (frequency['-']) {
		score += frequency['-'] * parameters[1]; // gap penalty
		steps++;
	}
	if (frequency['!']) {
		score += frequency['!'] * parameters[2]; // gap frame
		steps++;
	}
	if (frequency['*']) {
		score += frequency['*'] * parameters[3]; // stop cost
		steps++;
	}
	return score / denominator / steps;
}

void Profile :: CalcFrequenciesAA(int position) {
	memset(frequency, 0, sizeof(int)*128);
	for_each(sequences.begin(), sequences.end(), 
		[&](BioSeq* s) { 
			frequency[(unsigned char)s->TranslateNTtoAA(position)]++; 
		} );
}

void Profile :: CalcFrequenciesNT(int position) {
	memset(frequency, 0, sizeof(int)*128);
	for_each(sequences.begin(), sequences.end(), 
		[&](BioSeq* s) { 
			frequency[(unsigned char)s->nt_seq[position]]++; 
		} );
}