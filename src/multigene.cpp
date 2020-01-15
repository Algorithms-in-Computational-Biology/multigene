#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include "nnparams.h"
#include "thermalign.h"
#include "dinkelbach.h"
#include "multigene.h"
#include "stree_ukkonen.h"

// #define _DEBUG
#define _LOG 	

using namespace std;

int min_length = 18; //2
int max_length = 24; //3
float min_gc_content = 40.0; //40
float max_gc_content = 60.0; //70
float min_melting_temperature = 50.0; //50
float max_melting_temperature = 69.0; //65
float threshold_elongation_efficiency = 20.0;
int min_product_size = 550; //1
int max_product_size = 1030; //2
float primer_concentration = 0.000001;
float template_concentration = 0.000001;
float salt_concentration = 1; // K?
float t = 50.0f; // Temperature in C	

FILE *f = NULL;

// Elongation efficiency (primer/target)
float efficiency[4][4] = { 
	{ 0.0130, 0.3359, 0.0019, 1 }, 
	{ 0.1842, 0.0101, 1, 0.1928 }, 
	{ 0.0206, 1, 0.0377, 0.1495 }, 
	{ 1, 0.3035, 0.4501, 0.0150 } 
	};

float *new_array(int size) 
{
	float *array = (float *) malloc(size * sizeof(float));
	if (array == NULL) 
	{
		printf("Memoria insuficiente!");

		exit(0);
	}
	return array;
}

char *new_string(int length) 
{
	char *str = (char *) malloc((length + 1) * sizeof(char));
	if (str == NULL) 
	{
		printf("Memoria insuficiente!");

		exit(0);
	}
	return str;
}

Primer *new_primer(char *sequence, int length, int m) 
{
	Primer *primer = (Primer *) malloc(sizeof(Primer));
	if (primer == NULL) 
	{
		printf("Memoria insuficiente!");

		exit(0);
	}
	primer->target_id = 0;
	primer->sequence = new_string(length);
	strcpy(primer->sequence, sequence);
	primer->length = length;
	primer->forward = false;
	primer->reverse = false;
	primer->gc_content = 0.0;
	primer->temperature = 0.0;
	primer->free_energy = 0.0;
	primer->forward_elongation_efficiency = new_array(m);
	primer->reverse_elongation_efficiency = new_array(m);
	primer->next = NULL;

	return primer;
}

void delete_primer(Primer *primer) 
{
	free(primer->sequence);
	free(primer->forward_elongation_efficiency);
	free(primer->reverse_elongation_efficiency);
	free(primer);
}

Primer *copy_primer(Primer *primer, int m) 
{
	Primer *copy = new_primer(primer->sequence, primer->length, m);
	copy->target_id = primer->target_id;
	copy->position = primer->position;
	copy->forward = primer->forward;
	copy->reverse = primer->reverse;
	copy->gc_content = primer->gc_content;
	copy->temperature = primer->temperature;
	copy->free_energy = primer->free_energy;
	memcpy(copy->forward_elongation_efficiency,
			primer->forward_elongation_efficiency, m * sizeof(float));
	memcpy(copy->reverse_elongation_efficiency,
			primer->reverse_elongation_efficiency, m * sizeof(float));

	return copy;
}

void add_primer(Primer **head, Primer *primer) 
{
	primer->next = *head;
	*head = primer;
}

void delete_primer_list(Primer **primers) 
{
	Primer *primer = *primers;
	Primer *next;
	while (primer != NULL) 
	{
		next = primer->next;
		delete_primer(primer);
		primer = next;
	}
	*primers = NULL;
}

Pair *new_pair(Primer *forward, Primer *reverse, int m) 
{
	Pair *pair = (Pair *) malloc(sizeof(Pair));
	if (pair == NULL) 
	{
		printf("Memoria insuficiente!");

		exit(0);
	}
	pair->forward = copy_primer(forward, m);
	pair->reverse = copy_primer(reverse, m);
	pair->next = NULL;

	return pair;
}

Pair *copy_pair(Pair *pair, int m) 
{
	Pair *copy = new_pair(pair->forward, pair->reverse, m);
	copy->target_id = pair->target_id;
	copy->product_size = pair->product_size;
	copy->free_energy = pair->free_energy;

	return copy;
}

void delete_pair(Pair *pair) 
{
	delete_primer(pair->forward);
	delete_primer(pair->reverse);
	free(pair);
}

void add_pair(Pair **head, Pair *pair) 
{
	pair->next = *head;
	*head = pair;
}

void delete_pair_list(Pair **pairs) 
{
	Pair *pair = *pairs;
	Pair *next;
	while (pair != NULL) 
	{
		next = pair->next;
		delete_pair(pair);
		pair = next;
	}
	*pairs = NULL;
}

char *substring(char *str, int pos, int length) 
{
	char *sub = (char *) malloc((length + 1) * sizeof(char));
	if (sub == NULL) 
	{
		printf("Memoria insuficiente!");

		exit(0);
	}
	strncpy(sub, str + pos, length);
	sub[length] = '\0';

	return sub;
}

bool is_substring(SUFFIX_TREE tree, char *str, int n) 
{
	STREE_NODE node;
	int pos;
	return (stree_match(tree, str, n, &node, &pos) == n);
}

Primer *create_candidate_list(int target_id, char *target, int m, 
	SUFFIX_TREE tree)
{
	register int i, j;
	Primer *head = NULL;
	int n = strlen(target);
	for (i = min_length; i <= max_length; i++) 
	{
		for (j = 0; j < n - i + 1; j++) 
		{
			char *sequence = substring(target, j, i);
			// Search candidate in suffix tree
			if (!is_substring(tree, sequence, i)) 
			{
				Primer *primer = new_primer(sequence, i, m);
				primer->position = j + 1;
				primer->target_id = target_id;

				add_primer(&head, primer);
			}
			free(sequence);
		}
	}
	return head;
}

// Temp
int convert_base(char base) 
{
	if (base == 'A') 
	{
		return 1;
	} 
	else if (base == 'C') 
	{
		return 2;
	} 
	else if (base == 'G') 
	{
		return 3;
	} 
	else 
	{
		return 4;
	}
}
//

char *add_signs(char *sequence)
{
	int n = strlen(sequence);
	char *str = (char *) malloc((n + 3) * sizeof(char));
	str[0] = '$';
	for (int i = 0; i < n; i++) 
	{
		str[i + 1] = sequence[i];
	}
	str[n + 1] = '$';
	str[n + 2] = '\0';

	return str;
}

char get_complementary_base(char base) 
{
	if (base == 'A') 
	{
		return 'T';
	} 
	else if (base == 'G') 
	{
		return 'C';
	} 
	else if (base == 'C') 
	{
		return 'G';
	} 
	else if (base == 'T') 
	{
		return 'A';
	}
	return 'A';
}

char *calculate_reverse_complement(char *sequence, int n) 
{
	register int i;
	char *reverse = (char *) malloc((n + 1) * sizeof(char));
	for (i = n - 1; i >= 0; i--) 
	{
		reverse[n - i - 1] = get_complementary_base(sequence[i]);
	}
	reverse[n] = '\0';

	return reverse;
}

char *calculate_complement(char *sequence, int n) 
{
	register int i;
	char *complement = (char *) malloc((n + 1) * sizeof(char));
	for (i = 0; i < n; i++) 
	{
		complement[i] = get_complementary_base(sequence[i]);
	}
	complement[n] = '\0';

	return complement;
}

char *calculate_reverse(char *sequence, int n) 
{
	register int i;
	char *reverse = (char *) malloc((n + 1) * sizeof(char));
	for (i = n - 1; i >= 0; i--) 
	{
		reverse[n - i - 1] = sequence[i];
	}
	reverse[n] = '\0';

	return reverse;
}

int calculate_product_size(Primer *forward, Primer *reverse) 
{
	return reverse->position - (forward->position + forward->length);
}

float calculate_gc_content(char *sequence, int n) 
{
	register int i;
	float count = 0;
	for (i = 0; i < n; i++) 
	{
		if (sequence[i] == 'G' || sequence[i] == 'C') {
			count += 1;
		}
	}
	return 100 * count / (float) n;
}

bool is_perfect_complement(char *target, char *primer) 
{
	register int k;
	int primer_length = strlen(primer);
	int target_length = strlen(target);
	if (target_length != primer_length) 
	{
		return false;
	}
	for (k = 0; k < primer_length; k++) 
	{
		if (primer[k] != get_complementary_base(target[k])) 
		{
			return false;
		}
	}
	return true;
}

float calculate_free_energy(char *target1, char *primer1, PNNParams params,
		float t) 
{
	char *target = add_signs(target1);
	char *primer = add_signs(primer1);

	int primer_length = strlen(primer);
	int target_length = strlen(target);

	CThermAlign tAlign(target_length, primer_length, params);
	tAlign.InitStrings(target, primer, target_length, primer_length);
	tAlign.CalculateTable();

	float tempK = tAlign.GetMeltingTempK(tAlign.maxloci, tAlign.maxlocj);
	GAlign gAlign(target_length, primer_length, params);

	gAlign.InitStrings(target, primer, target_length, primer_length);
	dinkelbach dinkel(params, &gAlign);
	dinkel.iteration(tempK);

	free(target);
	free(primer);

	return -gAlign.GetFreeEnergyK(gAlign.maxloci, gAlign.maxlocj, t + 273.0f);
}

float calculate_melting_temperature(char *target1, char *primer1,
		PNNParams params) 
{
	char *target = add_signs(target1);
	char *primer = add_signs(primer1);

	int primer_length = strlen(primer);
	int target_length = strlen(target);

	CThermAlign tAlign(target_length, primer_length, params);
	tAlign.InitStrings(target, primer, target_length, primer_length);
	tAlign.CalculateTable();
	float temperature;
	if (is_perfect_complement(target, primer)) 
	{
		temperature = tAlign.GetMeltingTempC(tAlign.maxloci, tAlign.maxlocj);
	} 
	else 
	{
		float tempK = tAlign.GetMeltingTempK(tAlign.maxloci, tAlign.maxlocj);
		GAlign gAlign(target_length, primer_length, params);

		gAlign.InitStrings(target, primer, target_length, primer_length);
		dinkelbach dinkel(params, &gAlign);
		dinkel.iteration(tempK);
		temperature = gAlign.GetMeltingTempC(gAlign.maxloci, gAlign.maxlocj);
	}
	free(target);
	free(primer);

	return temperature;
}

BasePair *calculate_thermodynamic_alignment(char *target1, char *primer1,
		PNNParams params) 
{
	char *target = add_signs(target1);
	char *primer = add_signs(primer1);

	int primer_length = strlen(primer);
	int target_length = strlen(target);

	CThermAlign tAlign(target_length, primer_length, params);
	tAlign.InitStrings(target, primer, target_length, primer_length);
	tAlign.CalculateTable();
	float tempK = tAlign.GetMeltingTempK(tAlign.maxloci, tAlign.maxlocj);

	GAlign gAlign(target_length, primer_length, params);
	gAlign.InitStrings(target, primer, target_length, primer_length);
	dinkelbach dinkel(params, &gAlign);
	dinkel.iteration(tempK);
	gAlign.OutputLocalAlignment(cout);

	free(target);
	free(primer);

	BasePair *pairs = (BasePair *) malloc(2 * sizeof(BasePair));
	memcpy(pairs, gAlign.basePairs, 2 * sizeof(BasePair));

	return pairs;
}

// Validate GC content and temperature
bool validate(Primer *primer, PNNParams params) 
{
	float gc_content = calculate_gc_content(primer->sequence, primer->length);
	if (gc_content < min_gc_content || gc_content > max_gc_content) 
	{
		return false;
	}
	primer->gc_content = gc_content;
	char *complement = calculate_complement(primer->sequence, primer->length);
	float temperature = calculate_melting_temperature(primer->sequence,
			complement, params);
	// delete complement
	free(complement);
	if (temperature < min_melting_temperature
			|| temperature > max_melting_temperature) 
			{
		return false;
	}
	primer->temperature = temperature;

	return true;
}

bool validate_elongation_efficiency(Primer *primer, int m) 
{
	register int i;
	if (primer->forward) 
	{
		for (i = 0; i < m; i++) 
		{
			float elongation_efficiency = 
				primer->forward_elongation_efficiency[i];
			if (elongation_efficiency > threshold_elongation_efficiency) 
			{
				return false;
			}
		}
	}
	if (primer->reverse) 
	{
		for (i = 0; i < m; i++) 
		{
			float elongation_efficiency = 
				primer->reverse_elongation_efficiency[i];
			if (elongation_efficiency > threshold_elongation_efficiency) 
			{
				return false;

			}
		}
	}
	return true;
}

void print_result(Pair *result, int m) 
{
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
	for (Pair *pair = result; pair != NULL; pair = pair->next) 
	{
		Primer *forward = pair->forward;
		fprintf(f, "%10d %12d %8.2f ", pair->target_id, pair->product_size,
				pair->free_energy / 1000);

		fprintf(f, "%25s %6d %8d %10.2f %11.2f %10.2f ", forward->sequence,
				forward->length, forward->position, forward->gc_content,
				forward->temperature, forward->free_energy / 1000);

		for (i = 0; i < m; i++) 
		{
			fprintf(f, "%21.2f ", forward->forward_elongation_efficiency[i]);
		}
		Primer *reverse = pair->reverse;
		fprintf(f, "%25s %6d %8d %10.2f %11.2f %10.2f ", reverse->sequence,
				reverse->length, reverse->position, reverse->gc_content,
				reverse->temperature, reverse->free_energy / 1000);
		for (i = 0; i < m; i++) 
		{
			fprintf(f, "%21.2f", reverse->reverse_elongation_efficiency[i]);
		}
		fprintf(f, "\n");
	}
}

void print_primer(Primer *primer, int m) 
{
	register int i;
	fprintf(f, "%c %c %s %d %d %.2f%% %.2fC %.2f ", primer->forward ? 'F' : ' ',
			primer->reverse ? 'R' : ' ', primer->sequence, primer->length,
			primer->position, primer->gc_content, primer->temperature,
			primer->free_energy);

	if (primer->forward) 
	{
		for (i = 0; i < m; i++) 
		{
			fprintf(f, "(%.2f%%) ", primer->forward_elongation_efficiency[i]);
		}
	}
	if (primer->reverse) 
	{
		for (i = 0; i < m; i++) 
		{
			fprintf(f, "(%.2f%%) ", primer->reverse_elongation_efficiency[i]);
		}
	}
	fprintf(f, "\n");
}

void print_pair(Pair *pair, int m) 
{
	Primer *forward = pair->forward;
	Primer *reverse = pair->reverse;
	fprintf(f, "%d %d %.2f ", pair->target_id, pair->product_size,
			pair->free_energy);
	print_primer(forward, m);
	fprintf(f, "%d %d %.2f ", pair->target_id, pair->product_size,
			pair->free_energy);
	print_primer(reverse, m);
}

float calculate_elongation_efficiency(char target, char primer) 
{
	int i = convert_base(target);
	int j = convert_base(primer);

	return efficiency[i - 1][j - 1] * 100;
}

// Temp
bool is_mismatch(char target, char primer) 
{
	int num1 = convert_base(target);
	int num2 = convert_base(primer);

	return (num1 + num2) != 5;
}

SUFFIX_TREE build_gstree(char **targets, int m, int j) 
{
	register int i;
	SUFFIX_TREE tree;

	if (targets == NULL || m == 0) 
	{
		return NULL;
	}
	tree = stree_new_tree(4, 0, LINKED_LIST, 0);
	if (tree == NULL) 
	{
		return NULL;
	}
	printf("Building suffix tree of non targets");
	for (i = 0; i < m; i++) 
	{
		if (j != i) 
		{
			printf(" %d", i + 1);
			int length = strlen(targets[i]);
			if (stree_ukkonen_add_string(tree, targets[i], NULL, 
				length, i + 1) < 1) 
				{
				stree_delete_tree(tree);
				return NULL;
			}
		}
	}
	printf("...\n");
	return tree;
}

// extern "C" Pair *design(char **targets, int m) 
Pair *design(char **targets, int m) 
{
	// Always use the register keyword in loops, because it is faster
	register int i, j, k;

	//## 1 ##
	CNNParams params;
	params.InitParams(primer_concentration, template_concentration,
			salt_concentration, SALT_METHOD_SANTALUCIA);
	//
	int non = 0;
	Pair *result = NULL;
	Primer *primers = NULL;
	for (i = 0; i < m; i++) 
	{ // Targets
#ifdef _LOG
		printf("Building candidate list of target %d...\n", (i + 1));
#endif
#ifdef _DEBUG
		fprintf(f, "Building candidate list of target %d...\n", (i + 1));
#endif
		// Building suffix tree
		SUFFIX_TREE tree = build_gstree(targets, m, non++);
		Primer *candidates = create_candidate_list(i + 1, targets[i], m - 1,
				tree);
		stree_delete_tree(tree);
#ifdef _LOG
		printf("Calculating the elongation efficiency of candidates against non targets sequences...\n");
#endif
#ifdef _DEBUG
		fprintf(f, "Calculating the elongation efficiency of candidates against non targets sequences...\n");
		fprintf(f, "Primer|Sequence|Complement|Non target|Elongation efficiency\n");
#endif
		for (Primer *candidate = candidates; candidate != NULL; candidate =
				candidate->next) 
				{
			Primer *primer = copy_primer(candidate, m - 1);
			// Check all elongation efficiencies against non targets
			int index = 0;
			for (j = 0; j < m; j++)
			{ // Non targets
				if (j != i) 
				{
					char *non_target = targets[j];
					char *complement = calculate_complement(primer->sequence,
							primer->length);
					// Calculate thermodynamic alignment between non target and candidate
					//## 2 ##
					BasePair *base_pairs = calculate_thermodynamic_alignment(
							non_target, complement, &params); // SLOW!!!
					//
					// Checks two ends of duplex
					// ################
					primer->forward_elongation_efficiency[index] = 100.0;
					primer->reverse_elongation_efficiency[index] = 100.0;
					// ################
					for (k = 0; k < 2; k++) 
					{
						// ## 4 ##
						char c1 = non_target[base_pairs[k].i - 1];
						char c2 = complement[base_pairs[k].j - 1];
						//
						if (is_mismatch(c1, c2)) 
						{
							float elongation_efficiency =
									calculate_elongation_efficiency(c2, c1);
#ifdef _DEBUG
							fprintf(f, "%c %s %s %s %c/%c(%.2f%%)\n", k?'F':'R', primer->sequence, complement, nonTarget, c2, c1, elongationEfficiency);
#endif
							if (k) 
							{
								primer->forward = true;
								primer->forward_elongation_efficiency[index] =
										elongation_efficiency;
							} 
							else 
							{
								primer->reverse = true;
								primer->reverse_elongation_efficiency[index] =
										elongation_efficiency;
							}
						}
					}
					free(base_pairs);
					// delete complement
					free(complement);
					index++;
				}
			}
			// It is valid only if the primer has a terminal mismatch
			if (primer->forward || primer->reverse) 
			{
				add_primer(&primers, primer);
			}
			else 
			{
				// delete copy
				delete_primer(primer);
			}
		}
		// delete candidate list
		delete_primer_list(&candidates);
	}
	Pair *pairs = NULL;
	for (i = 0; i < m; i++) 
	{
#ifdef _LOG	
		printf("Validating candidates of target %d...\n", (i + 1));
#endif
#ifdef _DEBUG	
		fprintf(f, "Validating candidates of target %d...\n", (i + 1));
		fprintf(f, "Primer|Sequence|Length|Position|GC Content|Temperature|Elongation Efficiency\n");
#endif
		Primer *forwards = NULL;
		Primer *reverses = NULL;
		for (Primer *primer = primers; primer != NULL; primer = primer->next) 
		{
			if (i == (primer->target_id - 1)) 
			{
				if (validate(primer, &params)) 
				{ // ## 5 ##
#ifdef _DEBUG
						print_primer(primer, m - 1);
#endif
					char *sequence = primer->sequence;
					char *reverse = calculate_reverse(sequence,
							strlen(sequence));
					primer->free_energy = calculate_free_energy(sequence, reverse,
							&params, t);
					free(reverse);
					if (primer->forward) 
					{
						Primer *copy = copy_primer(primer, m - 1);
						add_primer(&forwards, copy);
					}
					if (primer->reverse) 
					{
						Primer *copy = copy_primer(primer, m - 1);
						add_primer(&reverses, copy);
					}
				}
			}
		}
#ifdef _LOG
		printf("Calculating the cartesian product between selected forward and reverse primers...\n");
#endif
#ifdef _DEBUG
		fprintf(f, "Calculating the cartesian product between selected forward and reverse primers...\n");
		fprintf(f, "Forward primer|Reverse Primer|Product size\n");
#endif
		// Cartesian product
		for (Primer *forward = forwards; forward != NULL;
				forward = forward->next) 
				{
			for (Primer *reverse = reverses; reverse != NULL;
					reverse = reverse->next) 
					{
				int product_size = calculate_product_size(forward, reverse);
#ifdef _DEBUG
				fprintf(f, "%s %s %d\n", forward->sequence, reverse->sequence, productSize);
#endif
				if (product_size >= min_product_size
						&& product_size <= max_product_size) 
						{
					forward->reverse = false;
					reverse->forward = false;

					Pair *pair = new_pair(forward, reverse, m - 1);
					pair->target_id = i + 1;
					pair->product_size = product_size;
					pair->free_energy = calculate_free_energy(forward->sequence,
							reverse->sequence, &params, t);
					add_pair(&pairs, pair);
				}
			}
		}
		// delete list
		delete_primer_list(&forwards);
		delete_primer_list(&reverses);
	}
	// delete primer list
	delete_primer_list(&primers);
#ifdef _LOG
	printf("Saving the results...\n");
#endif
#ifdef _DEBUG
	fprintf(f, "Saving the results...\n");
	fprintf(f, "Target Id.|Product size|Primer|Sequence|Length|Position|GC Content|Temperature|Elongation Efficiency \n");
#endif
	for (Pair *pair = pairs; pair != NULL; pair = pair->next) 
	{
		// Check the dimer formation
		Primer *forward = pair->forward;
		Primer *reverse = pair->reverse;

		if (validate_elongation_efficiency(forward, m - 1)
				&& validate_elongation_efficiency(reverse, m - 1)) 
				{
//#ifdef _DEBUG
			Pair *copy = copy_pair(pair, m - 1);
			add_pair(&result, copy);
//#endif
		}
	}
	// delete pairs
	delete_pair_list(&pairs);

	return result;
}

int main(void) {
	int m = 6;

	// Rapid memory test
	// char* targets[m] = {"ATGGCGGTGGCTTCGACCTCGCCGCTATCCGCCACGGCCCCCTCGCCGCC",
	// 	       "ATGGCGGTGGCCTCGACCTCGCCGCTGTCCGCCAAGCCCGCCACGGCCCCTTCGCCG"};

	// PtoxZm
	// char* targets[m] = {"ATGGCGGTGGCTTCGACCTCGCCGCTATCCGCCACGGCCCCCTCGCCGCCCGCTCCGGTGTCCGGGTTCCTCGCTCTCCCCGCCCGCCGCGGCTGCGCAACGCGCCTCGGCTCCGCCGCCGCGTGGAGGAGGCTTCGCGTGGAGGCGATCTGGAAGCAGCAGGAGAAGCAGCGGGCAGAGGTGTCCGTCGAGGAACCCGCCCCCGTCAGGGAGGCCGCCGCGCCCCTGGACGGAGTCGGAGCTGACGACCCCATGGTTCCTTCCTCGGACGAGAGCTGGGTGGTCAGGCTCGAGCAGTCGGTCAACATTTTCCTCACGGAATCGGTGATTATACTACTCAATACCGTGTACCGTGATCGGAACTACGCCAGGTTTTTTGTGCTGGAGACGATTGCCAGGGTGCCGTATTTCGCGTTCATATCGGTGCTTCACATGTATGAAACCTTTGGCTGGTGGAGACGAGCTGATTATCTAAAAGTTCACTTTGCGCAGAGCTTGAACGAGTTTCATCATCTCTTGATCATGGAAGAATTGGGTGGCAACGCTATATGGATTGATTGTTTCCTTGCTCGATTTATGGCGTTTTTTTACTACTTCATGACTGTTGCGATGTACATGTTGAGCCCACGAATGGCATATCACTTCTCTGAATGTGTGGAGAGACATGCGTACTCCACCTATGATAAGTTCCTCAAGCTCCATGAAGAGGAATTGAAAACACTACCAGCTCCAGAGGCAGCATTGAACTATTACCTGAATGAGGACCTTTACTTATTTGATGAGTTTCAGACAACAAGAATTCCATGTTCTAGGAGGCCTAAAATAGATAACTTGTATGATGTATTCGTCAATATACGAGATGACGAGGCAGAGCACTGCAAGACAATGAAGGCATGTCAAACACATGGAACTCTTCGTTCTCCTCACTCAATGCCGAACTGCTTAGAAGCTGCTACAGAATGTGTAATACCTGAAAACGATTGTGAAGGTATTGTGGACTGTGTCAAAAAGTCCCTTACAAAGTAA",
	// 		    "ATGGCGGTGGCCTCGACCTCGCCGCTGTCCGCCAAGCCCGCCACGGCCCCTTCGCCGCCCGCTCCCGGATCCGGGCTCCTCGCTCTCGGCGTTCGCCGCGCCCCCGCCACTGCCGCGTGGAGGAGGCTCCGCGTGGAGGCGATCAGGACGCAGCGAACGGAGGTGCCCGTCGAGGAGTCCGCCCCCGCCAGGGACGCCGCCGCTGCCGCGCCCCTGGACGGAAACGGAGCCGGAGCGGACGGCTCCGTGGTTCCTTCCTCGGACGACAGCTGGGTTGTCAAGCTCGAGCAGTCGTTCAACATTTTCGCCACGGATTCGGTGATTATGGTACTCAAGGGCGTGTACGGTGATCGGTACTACGCCAGGTTCTTTGCGCTGGAGACGATTGCGAGGGTGCCGTACTTCGCATTCATATCGGTGCTTCACTTGTATGCGACCTTTGGATGGTGGAGACGAGCTGATTACATAAAGGTTCACTTTGCGCAGAGCTGGAACGAGTTCCATCACCTCTTGATCATGGAAGAATTGGGTGGCGACTCTTTGTGGTTTGACTGTTTTCTTGCTCGGTTTATGGCATTCTTTTACTACTTCATGACTGTTGCAATGTACATGCTGAGCCCACGAATGGCATATCACTTTTCCGAATGTGTGGAGAGACATGCATATTCCACCTATGATGAGTTCCTCAAGCTCCATGAAGAGGAATTGAAAAGACTACCAGCTCCAGAGGCAGCATTGAACTATTACATGAATGAGGACCTTTACTTATTCGATGAGTTTCAGGCATCAAGAACTCCAGGTTCTAGGAGGCCTAAAATAGATAACTTATACGATGTATTCGTTAATATACGAGAAGATGAGGCAGAGCACTGCAAGACAATGAAGACCTGTCAAACACATGGAAATCTTCGTTCTCCTCATTCAACGCCGAACTGCTTAGAAGATGATACGGAATGTGTAATACCTGAAAACGACTGTGAAGGTATTGTGGACTGTGTCAAAAAGTCCCTTACAAAGTAA"};

	// Anacardium
	char* targets[m] = {"ATGAATCATAGTGCTGCAGCAAAATTATCCAGGTCGATTATATCTCAACTTAGTACTCGCGGCTTCTCAACAGCATCAATAAACTCTTCTGAAACCGCCCAAATCTTTGCCAGAGTGAGGCCTGCGTTCGGTACCCGGAATCTGAGTACTTCTGTTTCCCCCAATGGCTCGCCAAAGGACGAGAAAAATAACTCTGTCAGCTCTGACAAGTCGCCTGACGATAAGATCATTGGGAGCTATTGGGGCGTGGCACCCGCTAAGCTGACTAAAGAAGACGGTTCCGCCTGGAAATGGAATTGCTTCAGGCCCTGGGACACCTACACAGCGGACGTCTCCATTGATGTGACAAAGCACCACAAACCAGAAAACTTCAGGGATAAATTTGCTTACTGGACTGTCCAGACTCTCAAATACCCAACTCATTTATTCTTTCAGAGGCGCCACATGTGTCATGCAATGTTACTAGAGACAGTGGCGGCAGTTCCCGGAATGGTTGGAGGGATGCTTTTGCATTTCAAATCGCTGAGGCGATTCGAACACAGCGGCGGATGGATAAAAGCTCTACTTGAAGAAGCTGAAAACGAGCGTATGCATTTGATGACATTCATAGATTTGGCGAAGCCTGCCTGGTACGAACGTGCCCTTGTTTTTGCAGTCCAAGGAGTGTTTTTCAATGCATATTTTCTGGCCTACTTGGCTTCTCCAAAGCTCGCTCACCGCATAGTGGGCTACTTGGAAGAAGAAGCAGTGATATCTTACAGTGAATTCCTCAAAGATTTGGACAATGGTAGCTTTGAAAATGTCCCGGCTCCGGCAATCGCCATTGATTACTGGCGTTTGCCTCCGGATTCAACTCTCCGAGATGTTGTTGTGGTCATCAGAGCCGATGAGGCTCACCACCGTGACCTTAACCACTATGCATCGGATATTCAATGTCAAGGACATGAGCTCAAGGAAGCACCAGCGCCGATAGGATATCATTAA",
		"ATGTTGGCTGTGTTGGCTCCTCGATTGTTCTCCTCTGTAACCACTCGTGTGGTGACGGTGAGCCGATGGCAACCACGATGGTGACTGGTTACAAGCTTGGGATTGTTCACGTGAGAAACTGGAGCACTGTGGCTGTAGGTGAGAAGGAGCAGGAGGAGAAGAAACAGGCGACGGAAACAGCCGGTGTCGGTAACAACAAGGAAGAGAAAAGGATCGGAGTTACTGGGGTGTGGAAGTTCCGAAGGTCACTAAAGAAGATGGGACTGAATGGCGATGGAACTGCTTTAGGCCATGGGAGACTTACAAAGCTGACTTATCCATTGATCTGAAGAAACACCATGCGCCAGCAACATTTTTGGACAAAATGGCCTTTTGGACCGTGAAAGCTCTAAGATGGCCAACTGATTTGTTCTTCCAGAGGAGATATGGGTGCCGGGCAATGATGCTTGAGACGGTGGCAGCCGTGCCGGGAATGGTGGGAGGCTTGCTGTTGCACTGCAAGTCATTGAGGAAATTTGAGCACAGCGGGGGCTGGATCAAGGCGCTTTTGGAAGAAGCCGAAAACGAGAGAATGCATCTAATGACTTTCATGGAGGTGGCCAAGCCCAgGTGGTACGAGAGGGCTCTGGTTTTCGCAGTCCAAGGTGTATTCTTCAACGCCTACTTCCTGGGCTATTTGATCTCTCCGAAATTCGCTCACCGCATGGTCGGCTACCTGGAAGAAGAAGCAATTCACTCaTACACAGAATTCCTCAAAGAATTGGACAAAGGTAACATTGAAAACGTCCCAGCTCCTGCAATCGCCATAGACTACTGGCAAATGTCTCCGGACTCCACCTTgCGTGATGTTGTGATGGTGGTGAGAGCCGATGAGGCCCATCACCGaGATGTCAATCACTTCGCATCGGATGTACACTATCAAGGACGTGAACTGAGGGAGGCGCCAGCGCCAATTGGGTATCACTAA",
		"TTCTAATTTTCTCAGAGAAATGGCGTTTGTAAGGATTGTCGGTGATGCGAGGTCTTTTAAACGGCGGGAGGTACAGAAACCGGCACATTTGGACAGCGGTTTCCAGACGGCAGCTGGAGGTTTTGGAGAGAAACGGCTTGCGGTCTGCAGTTATGCAGCGTGGCGCTGAAGCGCAAGTGAAAGAGAAGAAAGAGGAAACGGAAAAGAAAGATGCCATGGTGTCCAGTTATTGGGGAATTGCCAGGCCAAAGATCACCAGAGAGGACGGCACTCTTGGAATTGCTTCATGCCATGGGAAACTGATCAGTTGGACTTATCTATTGATTTGAAGAAGCATCATGTTCCAAGGACATTTCTGGATAAATTTGCATACAAGACTGTCAAAATCCTTCGAGCTCCAACTGATATCTTTTTTCAGAGACATTATGGGTGTCGGGCAATGATGCTAGAAACTGTGGCTGCTGTGCCTGGAATAGTTGGGGGGATGCTGCTGCGTCTGAAGTCTCTCCGCAAGCTAGAGCAAAGTGGTGGCTGGGTCAAAGCCTTGCTCGAAGAAGCAGAGAATGAGAGGATGCATCTCATGACCATGGTGGAGCTTGTGCAGCCTAAATGGTATGAGAGGCTCCTGGTTCTTGCTGTGCAGGGAGTCCTTTTTAACTCTTTCTTTGTACTTTATGTACTCTCTCCCAAACTGGCACATAGAATTGTTGGGTACTTGGAGGAGGAAGCTATCCACTCGTATACAGAATATCTCAAGGATATTAGTAGCAGTGCAATTGAAAATGTTCCAGCCCCAGCTATCGCAATTGACTACTGGAGACTACCCAAGGATGCCACTCTTAAGGATGTTGTCACTGTGATCCGTGCTGGTAAGGCTCATCGCCGTGATGTCAACCATTTTGCTTCTAATATACAAGTTAAGGGGAAGGAATTGAGAGAAGCTCCGGCCCTCCCTTGGTTACATGGAGGCTTATTTTTCTGTTCAGGGTATAATTTCACTTTTGTATATGCTGCGGATGTACATAGGAGTCATAAATATGAACATAAATACTAAATAGTATACCGTACATGCATCTGTTGTAACAGAAAACAAATTATTGCAT",
		"TCTCAGGGAAATGGACTTTGTAACGATGACGATGATGCGAGGTATTTTAAATGGCGGGAGATCCGGCAACCGGTACATTTGGACGGCGGTTTCCAGACGGCAGCCCGAGGTTGTGGAGAGAAGCCGGTTGCAGTCTGCTGTTATGCAGTGGAGGAGGATGCTGAGCAGCAGCGGTGGAGCCGAATCGCAAGTGAAAGATGAGAAAGAGGAGAAGAAAGATGCGATGGTATCCAGTTACTGGGGAATTGCGAGGCCAAAGATCACCAGAGAGGACGGAACTGAGTGGCCTTGGAACTGCTTCATGCCATGGGAAACTTACCAGTCAAACTTATCTATTGATTTGAAGAAGCATCATGTACCAAGGACATTTCTGGATAAATTTGCATACAGGACTGTCAAAATCCTTAGAGTTCCAACTGATATCTTTTTTCAGAGACGTTATGGGTGTCGGGCAATGATGCTAGAAACTGTGGCTGCTGTGCCCGGAATGGTTGGAGGGATGTTGCAGCATCTGAAGTCTCTCCGCAAGATGGAACAAAGCGGTGGCTGGATCAAAGCCTTGCTTGAAGAAGCAGAGAATGAGAGGATGCATCTGATGACCATGGTGGAGCTTGTGCAGCCTCAATGGTATGAGAGGCTGCTGGTTCTTGCTGTGCAGGGAGTCTTCTTTAACTCATTCTTTGTGCTTTATGTACTCTCTCCCAAACTGGCACATAGAATTGTGGGGTATTTGGAGGAGGAAGCTATCCACTCGTATACAGAATATCTCAAGGATATCGATAGTGGTGCAATCGAAAATGTTCCAGCCCCAGCTATTGCAATTGACTACTGGAGACTTCCTAAGGATGCCACTCTCAAGGATGTTATCACTGTGATCCGTGCTGATGAGGCTCATCATCGCGATGTCAACCATTTTGCTTCTGATATACAAGTTCAGGGAAAGGAGTTGAGAGAGGCACCGGCCCCCCTTGGCTACCATTGAGGGCTATTTCACTGTTCAGGGTATAATTTCACTTTTGCATATGGCGAAAATAGTATATCAGTACTACTGTAGGAGTATGAAAATTGTATCTGCGGGAAACGAAAACAAATTATTGTATGTTGTGTATTATAACTTGATTGAATTTAACAAGCATTAAAGCTATAGAGATAATGCGTTTCTT",
		"TGGCCGTGCTTGGTCCTCGTTTGTTCTCCTCTGTTACCACTCGTGCGGTGTGCGGAGAGTCGGCGGCTACCTGGATTCTGACCGGTTACAATTACAAGCTTGGGATTCTTGGCGTGAGAAATAGGAGCACTTCGGCTTTTCATGTGAAGGAACAGAAGGAGGAGCAGAATGCGAAGGTAACCACCGGTGACGGTAACAATAACGAAGAGAAGAAGCTCACGAGTTACTGGGGTCTGGAAGCTCCTAAGTTCACCAAAGAAGATGGTACTTCATGGCCGTGGACCTGCTTTACGCCATGGGAGACTTACAGAGCTAACTTATCCATTGATCTGGAGAAACACCATGCGCCTGCAAAATTTATGGACAAAATGGCTTTTTGGATGGTCAAAATTCTCAGATGGCCTACTGATTTATTCTTTCAGAGGAGATTTGGCTGCCGGGCAATGATGCTCGAAACCGTGGCGGCAGTGCCTGGAATGGTGGGAGGCGTGCTGTTACACTGCAAGTCATTGAGGAAATTCGAGCAAAGCGGCGGCTGGATCAAGGCACTTTTGGAAGAAGCCGAGAATGAGAGAATGCACCTCATGACTTTCATGGAGGTCACAAAGCCCAGATGGTATGACAGAGCTCTGGTTTTGGTAGTCCAAGGTATATTCTTCAATGCATACTTTTTAGGCTATATGATTTCTCCGAAATTCGCTCACCGCGTGGTTGGGTACCTCGAAGAAGAAGCCATTCATTCTTACACAGAATTTCTCAAGGAATTGGACAACGGTAACATTGAAAATGTACCTGCTCCGGCGATTGCCCTTGATTACTGGCGTCTGCCTCCGGGATCAACTCTTCGTGATGTGGTGATGGTCGTCAGAGCCGATGAGGCCCACCACCGCGATGTCAATCACTTTGCATCGGATATCCATTGTCAAGGGCATGACTTGAGAGAGTGCGCAGCGCCGATTGGGTATCACTGACATGTTAAAA",
		"TGGCCGTGCTTGGTCCTCGTTTGTTCTCCTCTGTTACCACTCGTGCGGTGTGCGAAGAGTCGGCTGCTACCTGGATTCTGACCGGTTACAATTACAAGCTTGGGATTCTTGGCGTGAGAAATAGGAGCACTTCGGCTTTTCATGTGAAGGAACAGAAGGAGGAGCAGAATGCGAAGGTAACCACCGGTGACGGTAACAATAACGAAGAGAAGAAGCTCACGAGTTACTGGGGTCTGGAAGCTCCTAAGTTCACCAAAGAAGATGGTACTTCATGGCCGTGGACCTGCTTTACGCCATGGGAGACTTACAGAGCTAACTTATCCATTGATCTGGAGAAACACCATGCGCCTGCAAAATTTATGGACAAAATGGCTTTTTGGATGGTCAAAATTCTCAGATGGCCTACTGATTTATTCTTTCAGAGGAGATTTGGCTGCCGGGCAATGATGCTCGAAACCGTGGCGGCAGTGCCTGGAATGGTGGGAGGCGTGCTGTTACACTGCAAGTCATTGAGGAAATTCGAGCAAAGCGGCGGCTGGATCAAGGCACTTTTGGAAGAAGCCGAGAATGAGAGAATGCACCTCATGACTTTCATGGAGGTCACAAAGCCCAGATGGTATGACAGAGCTCTGGTTTTGGTAGTCCAAGGTATATTCTTCAATGCATACTTCTTAGGCTATATGATTTCTCCGAAATTCGCTCACCGCGTGGTAGGGTACCTCGAAGAAGAAGCCATTCATTCTTACACAGAATTTCTCAAGGAATTGGACAAAGGTAACATTGAAAATGTACCTGCTCCGGCAATTGCCCTTGATTACTGGCGCCTGCCTCCGGGATCAACTCTTCGTGATGTGGTGATGGTCGTCAGAGCCGATGAGGCCCACCACCGCGATGTCAATCATTTTGCATCGGATATCCATTGTCAAGGGCATGACTTGAGAGAGTGCGCAGCGCCGATTGGGTATCACTGACATGTTAAAA"};
	
	//?
	//char* targets[m] = {"ATGGTAGATGGAGTGATGATTCTTCCTGTGCTTATCATGATTGCTCTCCCCTCCCCTAGTATGGAAGATGAGAAGCCCAAGGTCAACCCCAAACTCTACATGTGTGTGTGTGAAGGTCTCTCCTGCGGTAATGAGGACCACTGTGAAGGCCAGCAGTGCTTTTCCTCACTGAGCATCAACGATGGCTTCCACGTCTACCAGAAAGGCTGCTTCCAGGTTTATGAGCAGGGAAAGATGACCTGTAAGACCCCGCCGTCCCCTGGCCAAGCCGTGGAGTGCTGCCAAGGGGACTGGTGTAACAGGAACATCACAGCCCAGCTGCCCACTAAAGGAAAATCCTTCCCTGGAACACAGAATTTCCACTTGGAGGTTGGCCTCATTATTCTCTCTGTAGTGTTCGCAGTATGTCTTTTAGCCTGCCTGCTGGGAGTTGCTCTCCGAAAATTTAAAAGGCGCAACCAAGAACGCCTCAATCCCCGAGATGTGGAGTATGGCACTATCGAAGGGCTCATCACCACCAATGTTGGAGACAGCACTTTAGCAGATTTATTGGATCATTCGTGTACATCAGGAAGTGGCTCTGGTCTTCCTTTTCTGGTACAAAGAACAGTGGCTCGCCAGATTACACTGTTGGAGTGTGTTGGGAAAGGCAGGTATGGTGAGGTGTGGAGGGGAAGCTGGCAAGGGGAGAATGTTGCCGTGAAGATCTTCTCCTCCCGTGATGAGAAGTCATGGTTCAGGGAAACGGAATTGTACAACACTGTGATGCTGAGGCATGAAAATATCTTAGGTTTCATTGCTTCAGACATGACATCAAGACACTCCAGTACCCAGCTGTGGTTAATTACACATTATCATGAAATGGGATCGTTGTACGACTATCTTCAGCTTACTACTCTGGATACAGTTAGCTGCCTTCGAATAGTGCTGTCCATAGCTAGTGGTCTTGCACATTTGCACATAGAGATATTTGGGACCCAAGGGAAACCAGCCATTGCCCATCGAGATTTAAAGAGCAAAAATATTCTGGTTAAGAAGAATGGACAGTGTTGCATAGCAGATTTGGGCCTGGCAGTCATGCATTCCCAGAGCACCAATCAGCTTGATGTGGGGAACAATCCCCGTGTGGGCACCAAGCGCTACATGGCCCCCGAAGTTCTAGATGAAACCATCCAGGTGGATTGTTTCGATTCTTATAAAAGGGTCGATATTTGGGCCTTTGGACTTGTTTTGTGGGAAGTGGCCAGGCGGATGGTGAGCAATGGTATAGTGGAGGATTACAAGCCACCGTTCTACGATGTGGTTCCCAATGACCCAAGTTTTGAAGATATGAGGAAGGTAGTCTGTGTGGATCAACAAAGGCCAAACATACCCAACAGATGGTTCTCAGACCCGACATTAACCTCTCTGGCCAAGCTAATGAAAGAATGCTGGTATCAAAATCCATCCGCAAGACTCACAGCACTGCGTATCAAAAAGACTTTGACCAAAATTGATAATTCTCTCGACAAATTGAAAACTGACTGTTGA",
	//                    "ATGGTAGATGGAGTGATGATTCTTCCTGTGCTTATCATGATTGCTCTCCCCTCCCCTAGTATGGAAGATGAGAAGCCCAAGGTCAGCCCCAAACTCTACATGTGTGTGTGTGAAGGTCTCTCCTGCGGTAATGAGGACCACTGTGAAGGCCAGCAGTGCTTTTCCTCACTGAGCATCAACGATGGCTTCCACGTCTACCAGAAAGGCTGCTTCCAGGTTTATGAGCAGGGAAAGATGACCTGTAAGACCCCGCCGTCCCCTGGCCAAGCCGTGGAGTGCTGCCAAGGGGACTGGTGTAACAGGAACATCACTGCCCAGCTGCCCACTAAAGGAAAATCCTTCCCTGGAACACAGAATTTCCACTTGGAGGTTGGCCTCATTATTCTCTCTGTAGTGTTTGCAGTATGTCTTTTAGCCTGCCTGCTGGGAGTTGCTCTCCGAAAATTTAAAAGGCGCAACCAAGAACGCCTCAATCCCCGAGACGTGGAGTATGGCACTATCGAAGGGCTCATCACCACCAATGTTGGAGACAGCACTTTAGCAGATTTATTGGATCATTCGTGTACATCAGGAAGTGGCTCTGGTCTTCCTTTTCTGGTACAAAGAACGGTGGCTCGCCAGATTACACTGTTGGAGTGCGTCGGGAAAGGCAGGTATGGTGAGGTGTGGAGGGGCAGCTGGCAAGGGGAGAATGTTGCCGTGAAGATCTTCTCCTCCCGTGATGAGAAGTCATGGTTCAGGGAAACGGAATTGTACAACACTGTGATGCTGAGGCATGAAAATATCTTAGGTTTCATTGCTTCAGACATGACATCAAGACACTCCAGTACCCAGCTGTGGTTAATTACACATTATCATGAAATGGGATCATTGTACGACTATCTTCAGCTTACTACTCTGGATACAGTTAGCTGCCTTCGAATAGTGCTGTCCATAGCTAGTGGTCTTGCACATTTGCACATAGAGATATTTGGGACCCAAGGGAAACCAGCCATTGCCCATCGAGATTTAAAGAGCAAAAATATTCTGGTTAAGAAGAATGGACAGTGTTGCATAGCAGATTTGGGCCTGGCAGTCATGCATTCCCAGAGCACCAATCAGCTTGATGTGGGGAACAATCCCCGTGTGGGCACCAAGCGCTACATGGCCCCCGAAGTTCTAGATGAAACCATCCAGGTGGATTGTTTCGATTCTTATAAAAGGGTCGATATTTGGGCCTTTGGACTTGTTTTGTGGGAAGTGGCCAGGCGGATGGTGAGCAATGGTATAGTGGAGGATTACAAGCCACCGTTCTACGATGTGGTTCCCAATGACCCAAGTTTTGAAGATATGAGGAAGGTAGTCTGTGTGGATCAACAAAGGCCAAACATACCCAACAGATGGTTCTCAGACCCGACATTAACCTCTCTGGCCAAGCTAATGAAAGAATGCTGGTATCAAAATCCATCCGCAAGACTCACAGCACTGCGTATCAAAAAGACTTTGACCAAAATTGATAATTCCCTCGACAAATTGAAAACTGACTGTTGA",
	//                    "ATGGTAGATGGAGTGATGATTCTTCCTGTGCTTATCATGATTGCTCTCCCCTCCCCTAGTATGGAAGATGAGAAGCCCAAGGTCAACCCCAAACTCTACATGTGTGTGTGTGAAGGTCTCTCCTGCGGTAATGAGGACCACTGTGAAGGCCAGCAGTGCTTTTCCTCACTGAGCATCAACGATGGCTTCCACGTCTACCAGAAAGGCTGCTTCCAGGTTTATGAGCAGGGAAAGATGACCTGTAAGACCCCGCCGTCCCCTGGCCAAGCCGTGGAGTGCTGCCAAGGGGACTGGTGTAACAGGAACATCACGGCCCAGCTGCCCACTAAAGGAAAATCCTTCCCTGGAACACAGAATTTCCACTTGGAGGTTGGCCTCATTATTCTCTCTGTAGTGTTCGCAGTATGTCTTTTAGCCTGCCTGCTGGGAGTTGCTCTCCGAAAATTTAAAAGGCGCAACCAAGAACGCCTCAATCCCCGAGACGTGGAGTATGGCACTATCGAAGGGCTCATCACCACCAATGTTGGAGACAGCACTTTAGCAGATTTATTGGATCATTCGTGTACATCAGGAAGTGGCTCTGGTCTTCCTTTTCTGGTACAAAGAACAGTGGCTCGCCAGATTACACTGTTGGAGTGTGTCGGGAAAGGCAGGTATGGTGAGGTGTGGAGGGGCAGCTGGCAAGGGGAGAATGTTGCCGTGAAGATCTTCTCCTCCCGTGATGAGAAGTCATGGTTCAGGGAAACGGAATTGTACAACACTGTGATGCTGAGGCATGAAAATATCTTAGGTTTCATTGCTTCAGACATGACATCAAGACACTCCAGTACCCAGCTGTGGTTAATTACACATTATCATGAAATGGGATCGTTGTACGACTATCTTCAGCTTACTACTCTGGATACAGTTAGCTGCCTTCGAATAGTGCTGTCCATAGCTAGTGGTCTTGCACATTTGCACATAGAGATATTTGGGACCCAAGGGAAACCAGCCATTGCCCATCGAGATTTAAAGAGCAAAAATATTCTGGTTAAGAAGAATGGACAGTGTTGCATAGCAGATTTGGGCCTGGCAGTCATGCATTCCCAGAGCACCAATCAGCTTGATGTGGGGAACAATCCCCGTGTGGGCACCAAGCGCTACATGGCCCCCGAAGTTCTAGATGAAACCATCCAGGTGGATTGTTTCGATTCTTATAAAAGGGTCGATATTTGGGCCTTTGGACTTGTTTTGTGGGAAGTGGCCAGGCGGATGGTGAGCAATGGTATAGTGGAGGATTACAAGCCACCGTTCTACGATGTGGTTCCCAATGACCCAAGTTTTGAAGATATGAGGAAGGTAGTCTGTGTGGATCAACAAAGGCCAAACATACCCAACAGATGGTTCTCAGACCCGACATTAACCTCTCTGGCCAAGCTAATGAAAGAATGCTGGTATCAAAATCCATCCGCAAGACTCACAGCACTGCGTATCAAAAAGACTTTGACCAAAATTGATAATTCCCTCGACAAATTGAAAACTGACTGTTGA"};

#ifdef _DEBUG
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
	print_result(result, m - 1);
	delete_pair_list(&result);
	fclose(f);

	return 0;
}