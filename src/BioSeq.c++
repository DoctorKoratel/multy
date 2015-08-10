#include "BioSeq.h"

/* auto generated code */

std::string BioSeq :: GetAAseq() {
	std::string result = "";
	for (unsigned int i = 0; i+2 < nt_seq.length(); i += 3) {
		if (nt_seq[i] == '-' && nt_seq[i+1] == '-' && nt_seq[i+2] == '-') {
			result += '-';
			continue;
		}
		else if (nt_seq[i] == '-' || nt_seq[i+1] == '-' || nt_seq[i+2] == '-') {
			result += '!';
			continue;
		}
		switch (nt_seq[i]) {
			case 'A':
				switch (nt_seq[i+1]) {
					case 'A':
						switch (nt_seq[i+2]) {
							case 'A':
							case 'G':
								result += 'K';
								break;
							default:
								result += 'N';
						}
						break;
					case 'C':
						result += 'T';
						break;
					case 'G':
						switch (nt_seq[i+2]) {
							case 'A':
							case 'G':
								result += 'R';
								break;
							default:
								result += 'S';
						}
						break;
					default:
						switch (nt_seq[i+2]) {
							case 'G':
								result += 'M'; //start
								break;
							default:
								result += 'I';
						}
				}
				break;
			case 'C':
				switch (nt_seq[i+1]) {
					case 'A':
						switch (nt_seq[i+2]) {
							case 'A':
							case 'G':
								result += 'Q';
								break;
							default:
								result += 'H';
						}
						break;
					case 'C':
						result += 'P';
						break;
					case 'G':
						result += 'R';
						break;
					default:
						result += 'L';
				}
				break;
			case 'G':
				switch (nt_seq[i+1]) {
					case 'A':
						switch (nt_seq[i+2]) {
							case 'A':
							case 'G':
								result += 'E';
								break;
							default:
								result += 'D';
						}
						break;
					case 'C':
						result += 'A';
						break;
					case 'G':
						result += 'G';
						break;
					default:
						result += 'V';
				}
				break;
			default:
				switch (nt_seq[i+1]) {
					case 'A':
						switch (nt_seq[i+2]) {
							case 'A':
							case 'G':
								result += '*'; //stop
								break;
							default:
								result += 'Y';
						}
						break;
					case 'C':
						result += 'S';
						break;
					case 'G':
						switch (nt_seq[i+2]) {
							case 'A':
								result += '*'; //stop
								break;
							case 'G':
								result += 'W';
								break;
							default:
								result += 'C';
						}
						break;
					default:
						switch (nt_seq[i+2]) {
							case 'A':
							case 'G':
								result += 'L';
								break;
							default:
								result += 'F';
						}
				}
		}
	}
	if (nt_seq.length() % 3) result += '!';
	return result;
}

void BioSeq :: PrintAA(std::ostream& out) {
	out << name << std::endl;
	for (unsigned int i = 0; i+2 < nt_seq.length(); i += 3) {
		if (nt_seq[i] == '-' && nt_seq[i+1] == '-' && nt_seq[i+2] == '-') {
			out << '-';
			continue;
		}
		else if (nt_seq[i] == '-' || nt_seq[i+1] == '-' || nt_seq[i+2] == '-') {
			out << '!';
			continue;
		}
		switch (nt_seq[i]) {
			case 'A':
				switch (nt_seq[i+1]) {
					case 'A':
						switch (nt_seq[i+2]) {
							case 'A':
							case 'G':
								out << 'K';
								break;
							default:
								out << 'N';
						}
						break;
					case 'C':
						out << 'T';
						break;
					case 'G':
						switch (nt_seq[i+2]) {
							case 'A':
							case 'G':
								out << 'R';
								break;
							default:
								out << 'S';
						}
						break;
					default:
						switch (nt_seq[i+2]) {
							case 'G':
								out << 'M'; //start
								break;
							default:
								out << 'I';
						}
				}
				break;
			case 'C':
				switch (nt_seq[i+1]) {
					case 'A':
						switch (nt_seq[i+2]) {
							case 'A':
							case 'G':
								out << 'Q';
								break;
							default:
								out << 'H';
						}
						break;
					case 'C':
						out << 'P';
						break;
					case 'G':
						out << 'R';
						break;
					default:
						out << 'L';
				}
				break;
			case 'G':
				switch (nt_seq[i+1]) {
					case 'A':
						switch (nt_seq[i+2]) {
							case 'A':
							case 'G':
								out << 'E';
								break;
							default:
								out << 'D';
						}
						break;
							case 'C':
								out << 'A';
								break;
							case 'G':
								out << 'G';
								break;
							default:
								out << 'V';
				}
				break;
			default:
				switch (nt_seq[i+1]) {
					case 'A':
						switch (nt_seq[i+2]) {
							case 'A':
							case 'G':
								out << '*'; //stop
								break;
							default:
								out << 'Y';
						}
						break;
					case 'C':
						out << 'S';
						break;
					case 'G':
						switch (nt_seq[i+2]) {
							case 'A':
								out << '*'; //stop
								break;
							case 'G':
								out << 'W';
								break;
							default:
								out << 'C';
						}
						break;
					default:
						switch (nt_seq[i+2]) {
							case 'A':
							case 'G':
								out << 'L';
								break;
							default:
								out << 'F';
						}
				}
		}
	}
	if (nt_seq.length() % 3) out << '!';
	out << std::endl;
}

char BioSeq :: TranslateNTtoAA(int i) const {
	if (nt_seq.length() < (unsigned int)i+3) return '!'; // check sequence tail
	if (nt_seq[i] == '-' && nt_seq[i+1] == '-' && nt_seq[i+2] == '-') return '-';
	if (nt_seq[i] == '-' || nt_seq[i+1] == '-' || nt_seq[i+2] == '-') return '!';
	char result;
	switch (nt_seq[i]) {
		case 'A':
			switch (nt_seq[i+1]) {
				case 'A':
					switch (nt_seq[i+2]) {
						case 'A':
						case 'G':
							result = 'K';
							break;
						default:
							result = 'N';
					}
					break;
				case 'C':
					result = 'T';
					break;
				case 'G':
					switch (nt_seq[i+2]) {
						case 'A':
						case 'G':
							result = 'R';
							break;
						default:
							result = 'S';
					}
					break;
				default:
					switch (nt_seq[i+2]) {
						case 'G':
							result = 'M'; //start
							break;
						default:
							result = 'I';
					}
			}
			break;
		case 'C':
			switch (nt_seq[i+1]) {
				case 'A':
					switch (nt_seq[i+2]) {
						case 'A':
						case 'G':
							result = 'Q';
							break;
						default:
							result = 'H';
					}
					break;
				case 'C':
					result = 'P';
					break;
				case 'G':
					result = 'R';
					break;
				default:
					result = 'L';
			}
			break;
		case 'G':
			switch (nt_seq[i+1]) {
				case 'A':
					switch (nt_seq[i+2]) {
						case 'A':
						case 'G':
							result = 'E';
							break;
						default:
							result = 'D';
					}
					break;
				case 'C':
					result = 'A';
					break;
				case 'G':
					result = 'G';
					break;
				default:
					result = 'V';
			}
			break;
		default:
			switch (nt_seq[i+1]) {
				case 'A':
					switch (nt_seq[i+2]) {
						case 'A':
						case 'G':
							result = '*'; //stop
							break;
						default:
							result = 'Y';
					}
					break;
				case 'C':
					result = 'S';
					break;
				case 'G':
					switch (nt_seq[i+2]) {
						case 'A':
							result = '*'; //stop
							break;
						case 'G':
							result = 'W';
							break;
						default:
							result = 'C';
					}
					break;
				default:
					switch (nt_seq[i+2]) {
						case 'A':
						case 'G':
							result = 'L';
							break;
						default:
							result = 'F';
					}
			}
	}
	return result;
}