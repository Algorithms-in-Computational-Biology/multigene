#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "nnparams.h"
#include "thermalign.h"
#include "dinkelbach.h"
#include "multigene.h"
#include "stree_ukkonen.h"

// #define _debug
#define _log 	

using namespace std;

int minLength = 18; //2
int maxLength = 24; //3
float minGCContent = 40.0; //40
float maxGCContent = 60.0; //70
float minMeltingTemperature = 50.0; //50
float maxMeltingTemperature = 69.0; //65
float thresholdElongationEfficiency = 20.0;
int minProductSize = 550; //1
int maxProductSize = 1030; //2
float primerConcentration = 0.000001;
float templateConcentration = 0.000001;
float saltConcentration = 1; // K?
float t = 50.0f; // Temperature in C	

FILE *f = NULL;

// (primer/target)
float efficiency[4][4] = { { 0.0130, 0.3359, 0.0019, 1 }, { 0.1842, 0.0101, 1,
		0.1928 }, { 0.0206, 1, 0.0377, 0.1495 }, { 1, 0.3035, 0.4501, 0.0150 } };

float *newArray(int size) {
	float *array = (float*) malloc(size * sizeof(float));
	if (array == NULL) {
		printf("Memoria insuficiente!");

		exit(0);
	}
	return array;
}

char *newString(int length) {
	char* str = (char*) malloc((length + 1) * sizeof(char));
	if (str == NULL) {
		printf("Memoria insuficiente!");

		exit(0);
	}
	return str;
}

Primer *newPrimer(char* sequence, int length, int m) {
	Primer* primer = (Primer*) malloc(sizeof(Primer));
	if (primer == NULL) {
		printf("Memoria insuficiente!");

		exit(0);
	}
	primer->targetId = 0;
	primer->sequence = newString(length);
	strcpy(primer->sequence, sequence);
	primer->length = length;
	primer->forward = false;
	primer->reverse = false;
	primer->GCContent = 0.0;
	primer->temperature = 0.0;
	primer->freeEnergy = 0.0;
	primer->forwardElongationEfficiency = newArray(m);
	primer->reverseElongationEfficiency = newArray(m);
	primer->next = NULL;

	return primer;
}

void deletePrimer(Primer *primer) {
	free(primer->sequence);
	free(primer->forwardElongationEfficiency);
	free(primer->reverseElongationEfficiency);
	free(primer);
}

Primer *copyPrimer(Primer *primer, int m) {
	Primer* copy = newPrimer(primer->sequence, primer->length, m);
	copy->targetId = primer->targetId;
	copy->position = primer->position;
	copy->forward = primer->forward;
	copy->reverse = primer->reverse;
	copy->GCContent = primer->GCContent;
	copy->temperature = primer->temperature;
	copy->freeEnergy = primer->freeEnergy;
	memcpy(copy->forwardElongationEfficiency,
			primer->forwardElongationEfficiency, m * sizeof(float));
	memcpy(copy->reverseElongationEfficiency,
			primer->reverseElongationEfficiency, m * sizeof(float));

	return copy;
}

void addPrimer(Primer **head, Primer *primer) {
	primer->next = *head;
	*head = primer;
}

void deletePrimerList(Primer **primers) {
	Primer *primer = *primers;
	Primer *next;
	while (primer != NULL) {
		next = primer->next;
		deletePrimer(primer);
		primer = next;
	}
	*primers = NULL;
}

Pair *newPair(Primer* forward, Primer *reverse, int m) {
	Pair* pair = (Pair*) malloc(sizeof(Pair));
	if (pair == NULL) {
		printf("Memoria insuficiente!");

		exit(0);
	}
	pair->forward = copyPrimer(forward, m);
	pair->reverse = copyPrimer(reverse, m);
	pair->next = NULL;

	return pair;
}

Pair *copyPair(Pair *pair, int m) {
	Pair* copy = newPair(pair->forward, pair->reverse, m);
	copy->targetId = pair->targetId;
	copy->productSize = pair->productSize;
	copy->freeEnergy = pair->freeEnergy;

	return copy;
}

void deletePair(Pair *pair) {
	deletePrimer(pair->forward);
	deletePrimer(pair->reverse);
	free(pair);
}

void addPair(Pair **head, Pair *pair) {
	pair->next = *head;
	*head = pair;
}

void deletePairList(Pair **pairs) {
	Pair *pair = *pairs;
	Pair *next;
	while (pair != NULL) {
		next = pair->next;
		deletePair(pair);
		pair = next;
	}
	*pairs = NULL;
}

char *substring(char* str, int pos, int length) {
	char* sub = (char*) malloc((length + 1) * sizeof(char));
	if (sub == NULL) {
		printf("Memoria insuficiente!");

		exit(0);
	}
	strncpy(sub, str + pos, length);
	sub[length] = '\0';

	return sub;
}

bool isSubstring(SUFFIX_TREE tree, char *str, int n) {
	STREE_NODE node;
	int pos;
	return (stree_match(tree, str, n, &node, &pos) == n);
}

Primer *createCandidateList(int targetId, char* target, int m, SUFFIX_TREE tree) {
	register int i, j;
	Primer *head = NULL;
	int n = strlen(target);
	for (i = minLength; i <= maxLength; i++) {
		for (j = 0; j < n - i + 1; j++) {
			char* sequence = substring(target, j, i);
			// Search candidate in suffix tree
			if (!isSubstring(tree, sequence, i)) {
				Primer *primer = newPrimer(sequence, i, m);
				primer->position = j + 1;
				primer->targetId = targetId;

				addPrimer(&head, primer);
			}
			free(sequence);
		}
	}
	return head;
}

// Temp
int convertBase(char base) {
	if (base == 'A') {
		return 1;
	} else if (base == 'C') {
		return 2;
	} else if (base == 'G') {
		return 3;
	} else {
		return 4;
	}
}
//

char *addSigns(char *sequence) {
	int n = strlen(sequence);
	char *str = (char*) malloc((n + 3) * sizeof(char));
	str[0] = '$';
	for (int i = 0; i < n; i++) {
		str[i + 1] = sequence[i];
	}
	str[n + 1] = '$';
	str[n + 2] = '\0';

	return str;
}

char getComplementaryBase(char base) {
	if (base == 'A') {
		return 'T';
	} else if (base == 'G') {
		return 'C';
	} else if (base == 'C') {
		return 'G';
	} else if (base == 'T') {
		return 'A';
	}
	return 'A';
}

char *calculateReverseComplement(char *sequence, int n) {
	register int i;
	char* reverse = (char*) malloc((n + 1) * sizeof(char));
	for (i = n - 1; i >= 0; i--) {
		reverse[n - i - 1] = getComplementaryBase(sequence[i]);
	}
	reverse[n] = '\0';

	return reverse;
}

char *calculateComplement(char *sequence, int n) {
	register int i;
	char* complement = (char*) malloc((n + 1) * sizeof(char));
	for (i = 0; i < n; i++) {
		complement[i] = getComplementaryBase(sequence[i]);
	}
	complement[n] = '\0';

	return complement;
}

char *calculateReverse(char *sequence, int n) {
	register int i;
	char* reverse = (char*) malloc((n + 1) * sizeof(char));
	for (i = n - 1; i >= 0; i--) {
		reverse[n - i - 1] = sequence[i];
	}
	reverse[n] = '\0';

	return reverse;
}

int calculateProductSize(Primer *forward, Primer *reverse) {
	return reverse->position - (forward->position + forward->length);
}

float calculateGCContent(char *sequence, int n) {
	register int i;
	float count = 0;
	for (i = 0; i < n; i++) {
		if (sequence[i] == 'G' || sequence[i] == 'C') {
			count += 1;
		}
	}
	return 100 * count / (float) n;
}

bool isPerfectComplement(char *target, char *primer) {
	register int k;
	int primerLength = strlen(primer);
	int targetLength = strlen(target);
	if (targetLength != primerLength) {
		return false;
	}
	for (k = 0; k < primerLength; k++) {
		if (primer[k] != getComplementaryBase(target[k])) {
			return false;
		}
	}
	return true;
}

float calculateFreeEnergy(char *target1, char *primer1, PNNParams params,
		float t) {
	char *target = addSigns(target1);
	char *primer = addSigns(primer1);

	int primerLength = strlen(primer);
	int targetLength = strlen(target);

	CThermAlign tAlign(targetLength, primerLength, params);
	tAlign.InitStrings(target, primer, targetLength, primerLength);
	tAlign.CalculateTable();

	float tempK = tAlign.GetMeltingTempK(tAlign.maxloci, tAlign.maxlocj);
	GAlign gAlign(targetLength, primerLength, params);

	gAlign.InitStrings(target, primer, targetLength, primerLength);
	dinkelbach dinkel(params, &gAlign);
	dinkel.iteration(tempK);

	free(target);
	free(primer);

	return -gAlign.GetFreeEnergyK(gAlign.maxloci, gAlign.maxlocj, t + 273.0f);
}

float calculateMeltingTemperature(char *target1, char *primer1,
		PNNParams params) {
	char *target = addSigns(target1);
	char *primer = addSigns(primer1);

	int primerLength = strlen(primer);
	int targetLength = strlen(target);

	CThermAlign tAlign(targetLength, primerLength, params);
	tAlign.InitStrings(target, primer, targetLength, primerLength);
	tAlign.CalculateTable();
	float temperature;
	if (isPerfectComplement(target, primer)) {
		temperature = tAlign.GetMeltingTempC(tAlign.maxloci, tAlign.maxlocj);
	} else {
		float tempK = tAlign.GetMeltingTempK(tAlign.maxloci, tAlign.maxlocj);
		GAlign gAlign(targetLength, primerLength, params);

		gAlign.InitStrings(target, primer, targetLength, primerLength);
		dinkelbach dinkel(params, &gAlign);
		dinkel.iteration(tempK);
		temperature = gAlign.GetMeltingTempC(gAlign.maxloci, gAlign.maxlocj);
	}
	free(target);
	free(primer);

	return temperature;
}

BasePair *calculateThermodynamicAlignment(char *target1, char *primer1,
		PNNParams params) {
	char *target = addSigns(target1);
	char *primer = addSigns(primer1);

	int primerLength = strlen(primer);
	int targetLength = strlen(target);

	CThermAlign tAlign(targetLength, primerLength, params);
	tAlign.InitStrings(target, primer, targetLength, primerLength);
	tAlign.CalculateTable();
	float tempK = tAlign.GetMeltingTempK(tAlign.maxloci, tAlign.maxlocj);

	GAlign gAlign(targetLength, primerLength, params);
	gAlign.InitStrings(target, primer, targetLength, primerLength);
	dinkelbach dinkel(params, &gAlign);
	dinkel.iteration(tempK);
	gAlign.OutputLocalAlignment(cout);

	free(target);
	free(primer);

	BasePair *pairs = (BasePair*) malloc(2 * sizeof(BasePair));
	memcpy(pairs, gAlign.basePairs, 2 * sizeof(BasePair));

	return pairs;
}

// Validate GC content and temperature
bool validate(Primer *primer, PNNParams params) {
	float GCContent = calculateGCContent(primer->sequence, primer->length);
	if (GCContent < minGCContent || GCContent > maxGCContent) {
		return false;
	}
	primer->GCContent = GCContent;
	char* complement = calculateComplement(primer->sequence, primer->length);
	float temperature = calculateMeltingTemperature(primer->sequence,
			complement, params);
	// delete complement
	free(complement);
	if (temperature < minMeltingTemperature
			|| temperature > maxMeltingTemperature) {
		return false;
	}
	primer->temperature = temperature;

	return true;
}

bool validateElongationEfficiency(Primer *primer, int m) {
	register int i;
	if (primer->forward) {
		for (i = 0; i < m; i++) {
			float elongationEfficiency = primer->forwardElongationEfficiency[i];
			if (elongationEfficiency > thresholdElongationEfficiency) {
				return false;
			}
		}
	}
	if (primer->reverse) {
		for (i = 0; i < m; i++) {
			float elongationEfficiency = primer->reverseElongationEfficiency[i];
			if (elongationEfficiency > thresholdElongationEfficiency) {
				return false;

			}
		}
	}
	return true;
}

void printResult(Pair *result, int m) {
	register int i;

	fprintf(f,
			"%10s %12s %8s %25s %6s %8s %10s %11s %10s %21s %25s %6s %8s %10s %11s %10s %21s\n%10s %12s %8s %25s %6s %8s %10s %11s %10s %21s %25s %6s %8s %10s %11s %10s %21s\n",
			"Target id.", "Product size", "dG Dimer", "Forward primer",
			"Length", "Position", "GC content", "Temperature", "dG Hairpin",
			"Elongation efficiency", "Reverse primer", "Length", "Position",
			"GC content", "Temperature", "dG Hairpin", "Elongation efficiency",
			"----------", "------------", "--------", "--------------",
			"------", "--------", "----------", "-----------", "----------",
			"---------------------", "--------------", "------", "--------",
			"----------", "-----------", "----------", "---------------------");
	for (Pair *pair = result; pair != NULL; pair = pair->next) {
		Primer *forward = pair->forward;
		fprintf(f, "%10d %12d %8.2f ", pair->targetId, pair->productSize,
				pair->freeEnergy / 1000);

		fprintf(f, "%25s %6d %8d %10.2f %11.2f %10.2f ", forward->sequence,
				forward->length, forward->position, forward->GCContent,
				forward->temperature, forward->freeEnergy / 1000);

		for (i = 0; i < m; i++) {
			fprintf(f, "%21.2f ", forward->forwardElongationEfficiency[i]);
		}
		Primer *reverse = pair->reverse;
		fprintf(f, "%25s %6d %8d %10.2f %11.2f %10.2f ", reverse->sequence,
				reverse->length, reverse->position, reverse->GCContent,
				reverse->temperature, reverse->freeEnergy / 1000);
		for (i = 0; i < m; i++) {
			fprintf(f, "%21.2f", reverse->reverseElongationEfficiency[i]);
		}
		fprintf(f, "\n");
	}
}

void printPrimer(Primer *primer, int m) {
	register int i;

	fprintf(f, "%c %c %s %d %d %.2f%% %.2fC %.2f ", primer->forward ? 'F' : ' ',
			primer->reverse ? 'R' : ' ', primer->sequence, primer->length,
			primer->position, primer->GCContent, primer->temperature,
			primer->freeEnergy);

	if (primer->forward) {
		for (i = 0; i < m; i++) {
			fprintf(f, "(%.2f%%) ", primer->forwardElongationEfficiency[i]);
		}
	}
	if (primer->reverse) {
		for (i = 0; i < m; i++) {
			fprintf(f, "(%.2f%%) ", primer->reverseElongationEfficiency[i]);
		}
	}
	fprintf(f, "\n");
}

void printPair(Pair *pair, int m) {
	Primer *forward = pair->forward;
	Primer *reverse = pair->reverse;
	fprintf(f, "%d %d %.2f ", pair->targetId, pair->productSize,
			pair->freeEnergy);
	printPrimer(forward, m);
	fprintf(f, "%d %d %.2f ", pair->targetId, pair->productSize,
			pair->freeEnergy);
	printPrimer(reverse, m);
}

float calculateElongationEfficiency(char target, char primer) {
	int i = convertBase(target);
	int j = convertBase(primer);

	return efficiency[i - 1][j - 1] * 100;
}

// Temp
bool isMismatch(char target, char primer) {
	int num1 = convertBase(target);
	int num2 = convertBase(primer);

	return (num1 + num2) != 5;
}

SUFFIX_TREE buildGSTree(char **targets, int m, int j) {
	register int i;
	SUFFIX_TREE tree;

	if (targets == NULL || m == 0) {
		return NULL;
	}
	tree = stree_new_tree(4, 0, LINKED_LIST, 0);
	if (tree == NULL) {
		return NULL;
	}
	printf("Building suffix tree of non targets");
	for (i = 0; i < m; i++) {
		if (j != i) {
			printf(" %d", i + 1);
			int length = strlen(targets[i]);
			if (stree_ukkonen_add_string(tree, targets[i], NULL, length, i + 1)
					< 1) {
				stree_delete_tree(tree);
				return NULL;
			}
		}
	}
	printf("...\n");
	return tree;
}

Pair *design(char **targets, int m) {
	// Always use the register keyword in loops, because it is faster
	register int i, j, k;

	//## 1 ##
	PNNParams params = new CNNParams();
	params->InitParams(primerConcentration, templateConcentration,
			saltConcentration, SALT_METHOD_SANTALUCIA);
	//
	int non = 0;
	Pair *result = NULL;
	Primer *primers = NULL;
	for (i = 0; i < m; i++) { // Targets
#ifdef _log
		printf("Building candidate list of target %d...\n", (i + 1));
#endif
#ifdef _debug
		fprintf(f, "Building candidate list of target %d...\n", (i + 1));
#endif
		// Building suffix tree
		SUFFIX_TREE tree = buildGSTree(targets, m, non++);
		Primer *candidates = createCandidateList(i + 1, targets[i], m - 1,
				tree);
		stree_delete_tree(tree);
#ifdef _log
		printf(
				"Calculating the elongation efficiency of candidates against non targets sequences...\n");
#endif
#ifdef _debug
		fprintf(f, "Calculating the elongation efficiency of candidates against non targets sequences...\n");
		fprintf(f, "Primer|Sequence|Complement|Non target|Elongation efficiency\n");
#endif
		for (Primer *candidate = candidates; candidate != NULL; candidate =
				candidate->next) {
			Primer *primer = copyPrimer(candidate, m - 1);
			// Check all elongation efficiencies against non targets
			int index = 0;
			for (j = 0; j < m; j++) { // Non targets
				if (j != i) {
					char *nonTarget = targets[j];
					char *complement = calculateComplement(primer->sequence,
							primer->length);
					// Calculate thermodynamic alignment between non target and candidate
					//## 2 ##
					BasePair *basePairs = calculateThermodynamicAlignment(
							nonTarget, complement, params); // SLOW!!!
					//
					// Checks two ends of duplex
					// ################
					primer->forwardElongationEfficiency[index] = 100.0;
					primer->reverseElongationEfficiency[index] = 100.0;
					// ################
					for (k = 0; k < 2; k++) {
						// ## 4 ##
						char c1 = nonTarget[basePairs[k].i - 1];
						char c2 = complement[basePairs[k].j - 1];
						//
						if (isMismatch(c1, c2)) {
							float elongationEfficiency =
									calculateElongationEfficiency(c2, c1);
#ifdef _debug
							fprintf(f, "%c %s %s %s %c/%c(%.2f%%)\n", k?'F':'R', primer->sequence, complement, nonTarget, c2, c1, elongationEfficiency);
#endif
							if (k) {
								primer->forward = true;
								primer->forwardElongationEfficiency[index] =
										elongationEfficiency;
							} else {
								primer->reverse = true;
								primer->reverseElongationEfficiency[index] =
										elongationEfficiency;
							}
						}
					}
					free(basePairs);
					// delete complement
					free(complement);
					index++;
				}
			}
			// It is valid only if the primer has a terminal mismatch
			if (primer->forward || primer->reverse) {
				addPrimer(&primers, primer);
			}
		}
		// delete candidate list
		deletePrimerList(&candidates);
	}
	Pair* pairs = NULL;
	for (i = 0; i < m; i++) {
#ifdef _log	
		printf("Validating candidates of target %d...\n", (i + 1));
#endif
#ifdef _debug	
		fprintf(f, "Validating candidates of target %d...\n", (i + 1));
		fprintf(f, "Primer|Sequence|Length|Position|GC Content|Temperature|Elongation Efficiency\n");
#endif
		Primer *forwards = NULL;
		Primer *reverses = NULL;
		for (Primer *primer = primers; primer != NULL; primer = primer->next) {
			if (i == (primer->targetId - 1)) {
				if (validate(primer, params)) { // ## 5 ##
#ifdef _debug
						printPrimer(primer, m - 1);
#endif
					char *sequence = primer->sequence;
					char *reverse = calculateReverse(sequence,
							strlen(sequence));
					primer->freeEnergy = calculateFreeEnergy(sequence, reverse,
							params, t);
					free(reverse);
					if (primer->forward) {
						Primer *copy = copyPrimer(primer, m - 1);
						addPrimer(&forwards, copy);
					}
					if (primer->reverse) {
						Primer *copy = copyPrimer(primer, m - 1);
						addPrimer(&reverses, copy);
					}
				}
			}
		}
#ifdef _log
		printf(
				"Calculating the cartesian product between selected forward and reverse primers...\n");
#endif
#ifdef _debug
		fprintf(f, "Calculating the cartesian product between selected forward and reverse primers...\n");
		fprintf(f, "Forward primer|Reverse Primer|Product size\n");
#endif
		// Cartesian product
		for (Primer *forward = forwards; forward != NULL;
				forward = forward->next) {
			for (Primer *reverse = reverses; reverse != NULL;
					reverse = reverse->next) {
				int productSize = calculateProductSize(forward, reverse);
#ifdef _debug
				fprintf(f, "%s %s %d\n", forward->sequence, reverse->sequence, productSize);
#endif
				if (productSize >= minProductSize
						&& productSize <= maxProductSize) {
					forward->reverse = false;
					reverse->forward = false;

					Pair* pair = newPair(forward, reverse, m - 1);
					pair->targetId = i + 1;
					pair->productSize = productSize;
					pair->freeEnergy = calculateFreeEnergy(forward->sequence,
							reverse->sequence, params, t);
					addPair(&pairs, pair);
				}
			}
		}
		// delete list
		deletePrimerList(&forwards);
		deletePrimerList(&reverses);
	}
	// delete primer list
	deletePrimerList(&primers);
#ifdef _log
	printf("Saving the results...\n");
#endif
#ifdef _debug
	fprintf(f, "Saving the results...\n");
	fprintf(f, "Target Id.|Product size|Primer|Sequence|Length|Position|GC Content|Temperature|Elongation Efficiency \n");
#endif
	for (Pair *pair = pairs; pair != NULL; pair = pair->next) {
		// Check the dimer formation
		Primer *forward = pair->forward;
		Primer *reverse = pair->reverse;

		if (validateElongationEfficiency(forward, m - 1)
				&& validateElongationEfficiency(reverse, m - 1)) {
//#ifdef _debug
			Pair *copy = copyPair(pair, m - 1);
			addPair(&result, copy);
//#endif
		}
	}
	// delete pairs
	deletePairList(&pairs);

	return result;
}

int main(void) {
	int m = 2;
	char* targets[m] = {"ATGGCGGTGGCTTCGACCTCGCCGCTATCCGCCACGGCCCCCTCGCCGCCCGCTCCGGTGTCCGGGTTCCTCGCTCTCCCCGCCCGCCGCGGCTGCGCAACGCGCCTCGGCTCCGCCGCCGCGTGGAGGAGGCTTCGCGTGGAGGCGATCTGGAAGCAGCAGGAGAAGCAGCGGGCAGAGGTGTCCGTCGAGGAACCCGCCCCCGTCAGGGAGGCCGCCGCGCCCCTGGACGGAGTCGGAGCTGACGACCCCATGGTTCCTTCCTCGGACGAGAGCTGGGTGGTCAGGCTCGAGCAGTCGGTCAACATTTTCCTCACGGAATCGGTGATTATACTACTCAATACCGTGTACCGTGATCGGAACTACGCCAGGTTTTTTGTGCTGGAGACGATTGCCAGGGTGCCGTATTTCGCGTTCATATCGGTGCTTCACATGTATGAAACCTTTGGCTGGTGGAGACGAGCTGATTATCTAAAAGTTCACTTTGCGCAGAGCTTGAACGAGTTTCATCATCTCTTGATCATGGAAGAATTGGGTGGCAACGCTATATGGATTGATTGTTTCCTTGCTCGATTTATGGCGTTTTTTTACTACTTCATGACTGTTGCGATGTACATGTTGAGCCCACGAATGGCATATCACTTCTCTGAATGTGTGGAGAGACATGCGTACTCCACCTATGATAAGTTCCTCAAGCTCCATGAAGAGGAATTGAAAACACTACCAGCTCCAGAGGCAGCATTGAACTATTACCTGAATGAGGACCTTTACTTATTTGATGAGTTTCAGACAACAAGAATTCCATGTTCTAGGAGGCCTAAAATAGATAACTTGTATGATGTATTCGTCAATATACGAGATGACGAGGCAGAGCACTGCAAGACAATGAAGGCATGTCAAACACATGGAACTCTTCGTTCTCCTCACTCAATGCCGAACTGCTTAGAAGCTGCTACAGAATGTGTAATACCTGAAAACGATTGTGAAGGTATTGTGGACTGTGTCAAAAAGTCCCTTACAAAGTAA",
			    "ATGGCGGTGGCCTCGACCTCGCCGCTGTCCGCCAAGCCCGCCACGGCCCCTTCGCCGCCCGCTCCCGGATCCGGGCTCCTCGCTCTCGGCGTTCGCCGCGCCCCCGCCACTGCCGCGTGGAGGAGGCTCCGCGTGGAGGCGATCAGGACGCAGCGAACGGAGGTGCCCGTCGAGGAGTCCGCCCCCGCCAGGGACGCCGCCGCTGCCGCGCCCCTGGACGGAAACGGAGCCGGAGCGGACGGCTCCGTGGTTCCTTCCTCGGACGACAGCTGGGTTGTCAAGCTCGAGCAGTCGTTCAACATTTTCGCCACGGATTCGGTGATTATGGTACTCAAGGGCGTGTACGGTGATCGGTACTACGCCAGGTTCTTTGCGCTGGAGACGATTGCGAGGGTGCCGTACTTCGCATTCATATCGGTGCTTCACTTGTATGCGACCTTTGGATGGTGGAGACGAGCTGATTACATAAAGGTTCACTTTGCGCAGAGCTGGAACGAGTTCCATCACCTCTTGATCATGGAAGAATTGGGTGGCGACTCTTTGTGGTTTGACTGTTTTCTTGCTCGGTTTATGGCATTCTTTTACTACTTCATGACTGTTGCAATGTACATGCTGAGCCCACGAATGGCATATCACTTTTCCGAATGTGTGGAGAGACATGCATATTCCACCTATGATGAGTTCCTCAAGCTCCATGAAGAGGAATTGAAAAGACTACCAGCTCCAGAGGCAGCATTGAACTATTACATGAATGAGGACCTTTACTTATTCGATGAGTTTCAGGCATCAAGAACTCCAGGTTCTAGGAGGCCTAAAATAGATAACTTATACGATGTATTCGTTAATATACGAGAAGATGAGGCAGAGCACTGCAAGACAATGAAGACCTGTCAAACACATGGAAATCTTCGTTCTCCTCATTCAACGCCGAACTGCTTAGAAGATGATACGGAATGTGTAATACCTGAAAACGACTGTGAAGGTATTGTGGACTGTGTCAAAAAGTCCCTTACAAAGTAA"};

	//char* targets[m] = {"ATGGCGGTGGCTTCGACCTCGCCGCTATCCGCCACGGCCCCCTCGCCGCC",
	//	       "ATGGCGGTGGCCTCGACCTCGCCGCTGTCCGCCAAGCCCGCCACGGCCCCTTCGCCG"};

	//char* targets[m] = {"ATGAATCATAGTGCTGCAGCAAAATTATCCAGGTCGATTATATCTCAACTTAGTACTCGCGGCTTCTCAACAGCATCAATAAACTCTTCTGAAACCGCCCAAATCTTTGCCAGAGTGAGGCCTGCGTTCGGTACCCGGAATCTGAGTACTTCTGTTTCCCCCAATGGCTCGCCAAAGGACGAGAAAAATAACTCTGTCAGCTCTGACAAGTCGCCTGACGATAAGATCATTGGGAGCTATTGGGGCGTGGCACCCGCTAAGCTGACTAAAGAAGACGGTTCCGCCTGGAAATGGAATTGCTTCAGGCCCTGGGACACCTACACAGCGGACGTCTCCATTGATGTGACAAAGCACCACAAACCAGAAAACTTCAGGGATAAATTTGCTTACTGGACTGTCCAGACTCTCAAATACCCAACTCATTTATTCTTTCAGAGGCGCCACATGTGTCATGCAATGTTACTAGAGACAGTGGCGGCAGTTCCCGGAATGGTTGGAGGGATGCTTTTGCATTTCAAATCGCTGAGGCGATTCGAACACAGCGGCGGATGGATAAAAGCTCTACTTGAAGAAGCTGAAAACGAGCGTATGCATTTGATGACATTCATAGATTTGGCGAAGCCTGCCTGGTACGAACGTGCCCTTGTTTTTGCAGTCCAAGGAGTGTTTTTCAATGCATATTTTCTGGCCTACTTGGCTTCTCCAAAGCTCGCTCACCGCATAGTGGGCTACTTGGAAGAAGAAGCAGTGATATCTTACAGTGAATTCCTCAAAGATTTGGACAATGGTAGCTTTGAAAATGTCCCGGCTCCGGCAATCGCCATTGATTACTGGCGTTTGCCTCCGGATTCAACTCTCCGAGATGTTGTTGTGGTCATCAGAGCCGATGAGGCTCACCACCGTGACCTTAACCACTATGCATCGGATATTCAATGTCAAGGACATGAGCTCAAGGAAGCACCAGCGCCGATAGGATATCATTAA",
	//	"ATGTTGGCTGTGTTGGCTCCTCGATTGTTCTCCTCTGTAACCACTCGTGTGGTGACGGTGAGCCGATGGCAACCACGATGGTGACTGGTTACAAGCTTGGGATTGTTCACGTGAGAAACTGGAGCACTGTGGCTGTAGGTGAGAAGGAGCAGGAGGAGAAGAAACAGGCGACGGAAACAGCCGGTGTCGGTAACAACAAGGAAGAGAAAAGGATCGGAGTTACTGGGGTGTGGAAGTTCCGAAGGTCACTAAAGAAGATGGGACTGAATGGCGATGGAACTGCTTTAGGCCATGGGAGACTTACAAAGCTGACTTATCCATTGATCTGAAGAAACACCATGCGCCAGCAACATTTTTGGACAAAATGGCCTTTTGGACCGTGAAAGCTCTAAGATGGCCAACTGATTTGTTCTTCCAGAGGAGATATGGGTGCCGGGCAATGATGCTTGAGACGGTGGCAGCCGTGCCGGGAATGGTGGGAGGCTTGCTGTTGCACTGCAAGTCATTGAGGAAATTTGAGCACAGCGGGGGCTGGATCAAGGCGCTTTTGGAAGAAGCCGAAAACGAGAGAATGCATCTAATGACTTTCATGGAGGTGGCCAAGCCCAgGTGGTACGAGAGGGCTCTGGTTTTCGCAGTCCAAGGTGTATTCTTCAACGCCTACTTCCTGGGCTATTTGATCTCTCCGAAATTCGCTCACCGCATGGTCGGCTACCTGGAAGAAGAAGCAATTCACTCaTACACAGAATTCCTCAAAGAATTGGACAAAGGTAACATTGAAAACGTCCCAGCTCCTGCAATCGCCATAGACTACTGGCAAATGTCTCCGGACTCCACCTTgCGTGATGTTGTGATGGTGGTGAGAGCCGATGAGGCCCATCACCGaGATGTCAATCACTTCGCATCGGATGTACACTATCAAGGACGTGAACTGAGGGAGGCGCCAGCGCCAATTGGGTATCACTAA",
	//	"TTCTAATTTTCTCAGAGAAATGGCGTTTGTAAGGATTGTCGGTGATGCGAGGTCTTTTAAACGGCGGGAGGTACAGAAACCGGCACATTTGGACAGCGGTTTCCAGACGGCAGCTGGAGGTTTTGGAGAGAAACGGCTTGCGGTCTGCAGTTATGCAGCGTGGCGCTGAAGCGCAAGTGAAAGAGAAGAAAGAGGAAACGGAAAAGAAAGATGCCATGGTGTCCAGTTATTGGGGAATTGCCAGGCCAAAGATCACCAGAGAGGACGGCACTCTTGGAATTGCTTCATGCCATGGGAAACTGATCAGTTGGACTTATCTATTGATTTGAAGAAGCATCATGTTCCAAGGACATTTCTGGATAAATTTGCATACAAGACTGTCAAAATCCTTCGAGCTCCAACTGATATCTTTTTTCAGAGACATTATGGGTGTCGGGCAATGATGCTAGAAACTGTGGCTGCTGTGCCTGGAATAGTTGGGGGGATGCTGCTGCGTCTGAAGTCTCTCCGCAAGCTAGAGCAAAGTGGTGGCTGGGTCAAAGCCTTGCTCGAAGAAGCAGAGAATGAGAGGATGCATCTCATGACCATGGTGGAGCTTGTGCAGCCTAAATGGTATGAGAGGCTCCTGGTTCTTGCTGTGCAGGGAGTCCTTTTTAACTCTTTCTTTGTACTTTATGTACTCTCTCCCAAACTGGCACATAGAATTGTTGGGTACTTGGAGGAGGAAGCTATCCACTCGTATACAGAATATCTCAAGGATATTAGTAGCAGTGCAATTGAAAATGTTCCAGCCCCAGCTATCGCAATTGACTACTGGAGACTACCCAAGGATGCCACTCTTAAGGATGTTGTCACTGTGATCCGTGCTGGTAAGGCTCATCGCCGTGATGTCAACCATTTTGCTTCTAATATACAAGTTAAGGGGAAGGAATTGAGAGAAGCTCCGGCCCTCCCTTGGTTACATGGAGGCTTATTTTTCTGTTCAGGGTATAATTTCACTTTTGTATATGCTGCGGATGTACATAGGAGTCATAAATATGAACATAAATACTAAATAGTATACCGTACATGCATCTGTTGTAACAGAAAACAAATTATTGCAT",
	//	"TCTCAGGGAAATGGACTTTGTAACGATGACGATGATGCGAGGTATTTTAAATGGCGGGAGATCCGGCAACCGGTACATTTGGACGGCGGTTTCCAGACGGCAGCCCGAGGTTGTGGAGAGAAGCCGGTTGCAGTCTGCTGTTATGCAGTGGAGGAGGATGCTGAGCAGCAGCGGTGGAGCCGAATCGCAAGTGAAAGATGAGAAAGAGGAGAAGAAAGATGCGATGGTATCCAGTTACTGGGGAATTGCGAGGCCAAAGATCACCAGAGAGGACGGAACTGAGTGGCCTTGGAACTGCTTCATGCCATGGGAAACTTACCAGTCAAACTTATCTATTGATTTGAAGAAGCATCATGTACCAAGGACATTTCTGGATAAATTTGCATACAGGACTGTCAAAATCCTTAGAGTTCCAACTGATATCTTTTTTCAGAGACGTTATGGGTGTCGGGCAATGATGCTAGAAACTGTGGCTGCTGTGCCCGGAATGGTTGGAGGGATGTTGCAGCATCTGAAGTCTCTCCGCAAGATGGAACAAAGCGGTGGCTGGATCAAAGCCTTGCTTGAAGAAGCAGAGAATGAGAGGATGCATCTGATGACCATGGTGGAGCTTGTGCAGCCTCAATGGTATGAGAGGCTGCTGGTTCTTGCTGTGCAGGGAGTCTTCTTTAACTCATTCTTTGTGCTTTATGTACTCTCTCCCAAACTGGCACATAGAATTGTGGGGTATTTGGAGGAGGAAGCTATCCACTCGTATACAGAATATCTCAAGGATATCGATAGTGGTGCAATCGAAAATGTTCCAGCCCCAGCTATTGCAATTGACTACTGGAGACTTCCTAAGGATGCCACTCTCAAGGATGTTATCACTGTGATCCGTGCTGATGAGGCTCATCATCGCGATGTCAACCATTTTGCTTCTGATATACAAGTTCAGGGAAAGGAGTTGAGAGAGGCACCGGCCCCCCTTGGCTACCATTGAGGGCTATTTCACTGTTCAGGGTATAATTTCACTTTTGCATATGGCGAAAATAGTATATCAGTACTACTGTAGGAGTATGAAAATTGTATCTGCGGGAAACGAAAACAAATTATTGTATGTTGTGTATTATAACTTGATTGAATTTAACAAGCATTAAAGCTATAGAGATAATGCGTTTCTT",
	//	"TGGCCGTGCTTGGTCCTCGTTTGTTCTCCTCTGTTACCACTCGTGCGGTGTGCGGAGAGTCGGCGGCTACCTGGATTCTGACCGGTTACAATTACAAGCTTGGGATTCTTGGCGTGAGAAATAGGAGCACTTCGGCTTTTCATGTGAAGGAACAGAAGGAGGAGCAGAATGCGAAGGTAACCACCGGTGACGGTAACAATAACGAAGAGAAGAAGCTCACGAGTTACTGGGGTCTGGAAGCTCCTAAGTTCACCAAAGAAGATGGTACTTCATGGCCGTGGACCTGCTTTACGCCATGGGAGACTTACAGAGCTAACTTATCCATTGATCTGGAGAAACACCATGCGCCTGCAAAATTTATGGACAAAATGGCTTTTTGGATGGTCAAAATTCTCAGATGGCCTACTGATTTATTCTTTCAGAGGAGATTTGGCTGCCGGGCAATGATGCTCGAAACCGTGGCGGCAGTGCCTGGAATGGTGGGAGGCGTGCTGTTACACTGCAAGTCATTGAGGAAATTCGAGCAAAGCGGCGGCTGGATCAAGGCACTTTTGGAAGAAGCCGAGAATGAGAGAATGCACCTCATGACTTTCATGGAGGTCACAAAGCCCAGATGGTATGACAGAGCTCTGGTTTTGGTAGTCCAAGGTATATTCTTCAATGCATACTTTTTAGGCTATATGATTTCTCCGAAATTCGCTCACCGCGTGGTTGGGTACCTCGAAGAAGAAGCCATTCATTCTTACACAGAATTTCTCAAGGAATTGGACAACGGTAACATTGAAAATGTACCTGCTCCGGCGATTGCCCTTGATTACTGGCGTCTGCCTCCGGGATCAACTCTTCGTGATGTGGTGATGGTCGTCAGAGCCGATGAGGCCCACCACCGCGATGTCAATCACTTTGCATCGGATATCCATTGTCAAGGGCATGACTTGAGAGAGTGCGCAGCGCCGATTGGGTATCACTGACATGTTAAAA",
	//	"TGGCCGTGCTTGGTCCTCGTTTGTTCTCCTCTGTTACCACTCGTGCGGTGTGCGAAGAGTCGGCTGCTACCTGGATTCTGACCGGTTACAATTACAAGCTTGGGATTCTTGGCGTGAGAAATAGGAGCACTTCGGCTTTTCATGTGAAGGAACAGAAGGAGGAGCAGAATGCGAAGGTAACCACCGGTGACGGTAACAATAACGAAGAGAAGAAGCTCACGAGTTACTGGGGTCTGGAAGCTCCTAAGTTCACCAAAGAAGATGGTACTTCATGGCCGTGGACCTGCTTTACGCCATGGGAGACTTACAGAGCTAACTTATCCATTGATCTGGAGAAACACCATGCGCCTGCAAAATTTATGGACAAAATGGCTTTTTGGATGGTCAAAATTCTCAGATGGCCTACTGATTTATTCTTTCAGAGGAGATTTGGCTGCCGGGCAATGATGCTCGAAACCGTGGCGGCAGTGCCTGGAATGGTGGGAGGCGTGCTGTTACACTGCAAGTCATTGAGGAAATTCGAGCAAAGCGGCGGCTGGATCAAGGCACTTTTGGAAGAAGCCGAGAATGAGAGAATGCACCTCATGACTTTCATGGAGGTCACAAAGCCCAGATGGTATGACAGAGCTCTGGTTTTGGTAGTCCAAGGTATATTCTTCAATGCATACTTCTTAGGCTATATGATTTCTCCGAAATTCGCTCACCGCGTGGTAGGGTACCTCGAAGAAGAAGCCATTCATTCTTACACAGAATTTCTCAAGGAATTGGACAAAGGTAACATTGAAAATGTACCTGCTCCGGCAATTGCCCTTGATTACTGGCGCCTGCCTCCGGGATCAACTCTTCGTGATGTGGTGATGGTCGTCAGAGCCGATGAGGCCCACCACCGCGATGTCAATCATTTTGCATCGGATATCCATTGTCAAGGGCATGACTTGAGAGAGTGCGCAGCGCCGATTGGGTATCACTGACATGTTAAAA"};
	//

	//char* targets[m] = {"ATGGTAGATGGAGTGATGATTCTTCCTGTGCTTATCATGATTGCTCTCCCCTCCCCTAGTATGGAAGATGAGAAGCCCAAGGTCAACCCCAAACTCTACATGTGTGTGTGTGAAGGTCTCTCCTGCGGTAATGAGGACCACTGTGAAGGCCAGCAGTGCTTTTCCTCACTGAGCATCAACGATGGCTTCCACGTCTACCAGAAAGGCTGCTTCCAGGTTTATGAGCAGGGAAAGATGACCTGTAAGACCCCGCCGTCCCCTGGCCAAGCCGTGGAGTGCTGCCAAGGGGACTGGTGTAACAGGAACATCACAGCCCAGCTGCCCACTAAAGGAAAATCCTTCCCTGGAACACAGAATTTCCACTTGGAGGTTGGCCTCATTATTCTCTCTGTAGTGTTCGCAGTATGTCTTTTAGCCTGCCTGCTGGGAGTTGCTCTCCGAAAATTTAAAAGGCGCAACCAAGAACGCCTCAATCCCCGAGATGTGGAGTATGGCACTATCGAAGGGCTCATCACCACCAATGTTGGAGACAGCACTTTAGCAGATTTATTGGATCATTCGTGTACATCAGGAAGTGGCTCTGGTCTTCCTTTTCTGGTACAAAGAACAGTGGCTCGCCAGATTACACTGTTGGAGTGTGTTGGGAAAGGCAGGTATGGTGAGGTGTGGAGGGGAAGCTGGCAAGGGGAGAATGTTGCCGTGAAGATCTTCTCCTCCCGTGATGAGAAGTCATGGTTCAGGGAAACGGAATTGTACAACACTGTGATGCTGAGGCATGAAAATATCTTAGGTTTCATTGCTTCAGACATGACATCAAGACACTCCAGTACCCAGCTGTGGTTAATTACACATTATCATGAAATGGGATCGTTGTACGACTATCTTCAGCTTACTACTCTGGATACAGTTAGCTGCCTTCGAATAGTGCTGTCCATAGCTAGTGGTCTTGCACATTTGCACATAGAGATATTTGGGACCCAAGGGAAACCAGCCATTGCCCATCGAGATTTAAAGAGCAAAAATATTCTGGTTAAGAAGAATGGACAGTGTTGCATAGCAGATTTGGGCCTGGCAGTCATGCATTCCCAGAGCACCAATCAGCTTGATGTGGGGAACAATCCCCGTGTGGGCACCAAGCGCTACATGGCCCCCGAAGTTCTAGATGAAACCATCCAGGTGGATTGTTTCGATTCTTATAAAAGGGTCGATATTTGGGCCTTTGGACTTGTTTTGTGGGAAGTGGCCAGGCGGATGGTGAGCAATGGTATAGTGGAGGATTACAAGCCACCGTTCTACGATGTGGTTCCCAATGACCCAAGTTTTGAAGATATGAGGAAGGTAGTCTGTGTGGATCAACAAAGGCCAAACATACCCAACAGATGGTTCTCAGACCCGACATTAACCTCTCTGGCCAAGCTAATGAAAGAATGCTGGTATCAAAATCCATCCGCAAGACTCACAGCACTGCGTATCAAAAAGACTTTGACCAAAATTGATAATTCTCTCGACAAATTGAAAACTGACTGTTGA",
	//                    "ATGGTAGATGGAGTGATGATTCTTCCTGTGCTTATCATGATTGCTCTCCCCTCCCCTAGTATGGAAGATGAGAAGCCCAAGGTCAGCCCCAAACTCTACATGTGTGTGTGTGAAGGTCTCTCCTGCGGTAATGAGGACCACTGTGAAGGCCAGCAGTGCTTTTCCTCACTGAGCATCAACGATGGCTTCCACGTCTACCAGAAAGGCTGCTTCCAGGTTTATGAGCAGGGAAAGATGACCTGTAAGACCCCGCCGTCCCCTGGCCAAGCCGTGGAGTGCTGCCAAGGGGACTGGTGTAACAGGAACATCACTGCCCAGCTGCCCACTAAAGGAAAATCCTTCCCTGGAACACAGAATTTCCACTTGGAGGTTGGCCTCATTATTCTCTCTGTAGTGTTTGCAGTATGTCTTTTAGCCTGCCTGCTGGGAGTTGCTCTCCGAAAATTTAAAAGGCGCAACCAAGAACGCCTCAATCCCCGAGACGTGGAGTATGGCACTATCGAAGGGCTCATCACCACCAATGTTGGAGACAGCACTTTAGCAGATTTATTGGATCATTCGTGTACATCAGGAAGTGGCTCTGGTCTTCCTTTTCTGGTACAAAGAACGGTGGCTCGCCAGATTACACTGTTGGAGTGCGTCGGGAAAGGCAGGTATGGTGAGGTGTGGAGGGGCAGCTGGCAAGGGGAGAATGTTGCCGTGAAGATCTTCTCCTCCCGTGATGAGAAGTCATGGTTCAGGGAAACGGAATTGTACAACACTGTGATGCTGAGGCATGAAAATATCTTAGGTTTCATTGCTTCAGACATGACATCAAGACACTCCAGTACCCAGCTGTGGTTAATTACACATTATCATGAAATGGGATCATTGTACGACTATCTTCAGCTTACTACTCTGGATACAGTTAGCTGCCTTCGAATAGTGCTGTCCATAGCTAGTGGTCTTGCACATTTGCACATAGAGATATTTGGGACCCAAGGGAAACCAGCCATTGCCCATCGAGATTTAAAGAGCAAAAATATTCTGGTTAAGAAGAATGGACAGTGTTGCATAGCAGATTTGGGCCTGGCAGTCATGCATTCCCAGAGCACCAATCAGCTTGATGTGGGGAACAATCCCCGTGTGGGCACCAAGCGCTACATGGCCCCCGAAGTTCTAGATGAAACCATCCAGGTGGATTGTTTCGATTCTTATAAAAGGGTCGATATTTGGGCCTTTGGACTTGTTTTGTGGGAAGTGGCCAGGCGGATGGTGAGCAATGGTATAGTGGAGGATTACAAGCCACCGTTCTACGATGTGGTTCCCAATGACCCAAGTTTTGAAGATATGAGGAAGGTAGTCTGTGTGGATCAACAAAGGCCAAACATACCCAACAGATGGTTCTCAGACCCGACATTAACCTCTCTGGCCAAGCTAATGAAAGAATGCTGGTATCAAAATCCATCCGCAAGACTCACAGCACTGCGTATCAAAAAGACTTTGACCAAAATTGATAATTCCCTCGACAAATTGAAAACTGACTGTTGA",
	//                    "ATGGTAGATGGAGTGATGATTCTTCCTGTGCTTATCATGATTGCTCTCCCCTCCCCTAGTATGGAAGATGAGAAGCCCAAGGTCAACCCCAAACTCTACATGTGTGTGTGTGAAGGTCTCTCCTGCGGTAATGAGGACCACTGTGAAGGCCAGCAGTGCTTTTCCTCACTGAGCATCAACGATGGCTTCCACGTCTACCAGAAAGGCTGCTTCCAGGTTTATGAGCAGGGAAAGATGACCTGTAAGACCCCGCCGTCCCCTGGCCAAGCCGTGGAGTGCTGCCAAGGGGACTGGTGTAACAGGAACATCACGGCCCAGCTGCCCACTAAAGGAAAATCCTTCCCTGGAACACAGAATTTCCACTTGGAGGTTGGCCTCATTATTCTCTCTGTAGTGTTCGCAGTATGTCTTTTAGCCTGCCTGCTGGGAGTTGCTCTCCGAAAATTTAAAAGGCGCAACCAAGAACGCCTCAATCCCCGAGACGTGGAGTATGGCACTATCGAAGGGCTCATCACCACCAATGTTGGAGACAGCACTTTAGCAGATTTATTGGATCATTCGTGTACATCAGGAAGTGGCTCTGGTCTTCCTTTTCTGGTACAAAGAACAGTGGCTCGCCAGATTACACTGTTGGAGTGTGTCGGGAAAGGCAGGTATGGTGAGGTGTGGAGGGGCAGCTGGCAAGGGGAGAATGTTGCCGTGAAGATCTTCTCCTCCCGTGATGAGAAGTCATGGTTCAGGGAAACGGAATTGTACAACACTGTGATGCTGAGGCATGAAAATATCTTAGGTTTCATTGCTTCAGACATGACATCAAGACACTCCAGTACCCAGCTGTGGTTAATTACACATTATCATGAAATGGGATCGTTGTACGACTATCTTCAGCTTACTACTCTGGATACAGTTAGCTGCCTTCGAATAGTGCTGTCCATAGCTAGTGGTCTTGCACATTTGCACATAGAGATATTTGGGACCCAAGGGAAACCAGCCATTGCCCATCGAGATTTAAAGAGCAAAAATATTCTGGTTAAGAAGAATGGACAGTGTTGCATAGCAGATTTGGGCCTGGCAGTCATGCATTCCCAGAGCACCAATCAGCTTGATGTGGGGAACAATCCCCGTGTGGGCACCAAGCGCTACATGGCCCCCGAAGTTCTAGATGAAACCATCCAGGTGGATTGTTTCGATTCTTATAAAAGGGTCGATATTTGGGCCTTTGGACTTGTTTTGTGGGAAGTGGCCAGGCGGATGGTGAGCAATGGTATAGTGGAGGATTACAAGCCACCGTTCTACGATGTGGTTCCCAATGACCCAAGTTTTGAAGATATGAGGAAGGTAGTCTGTGTGGATCAACAAAGGCCAAACATACCCAACAGATGGTTCTCAGACCCGACATTAACCTCTCTGGCCAAGCTAATGAAAGAATGCTGGTATCAAAATCCATCCGCAAGACTCACAGCACTGCGTATCAAAAAGACTTTGACCAAAATTGATAATTCCCTCGACAAATTGAAAACTGACTGTTGA"};

#ifdef _debug
	fprintf(f, "Targets\n");
	for(int i = 0; i < m; i++)
	{
		fprintf(f, "%d. %s\n", i + 1, targets[i]);
	}
#endif
	f = fopen("resultado.txt", "w");
	if (f == NULL) {
		printf("Error opening file!\n");
		exit(1);
	}
	Pair *result = design(targets, m);
	printResult(result, m - 1);
	deletePairList(&result);
	fclose(f);

	return 0;
}

