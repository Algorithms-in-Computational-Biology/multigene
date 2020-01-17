#ifndef MULTIGENE_H
#define MULTIGENE_H

struct TPrimer 
{
    char *target_id;			// Target sequence identification
    char *sequence;				// Primer sequence
    int length;					// Primer length
    int position;				// Begin position in target
    bool forward;				// Indicate whether primer is forward
    bool reverse;				// Indicate whether primer is reverse
    float gc_content;			// GC content
    float temperature;			// Melting temperature in C
    // float free_energy;		// Free energy change (dG) of hairpin formation in cal/mol 
    // int number;				// Number of non target sequences
    float *forward_elongation_efficiency; // Forward elongation efficiency
    float *reverse_elongation_efficiency; // Reverse elongation efficiency
    struct TPrimer *next;
};
typedef TPrimer Primer;

struct TPair 
{
    char *target_id;			// Target identification
    int product_size;			// PCR Product size 
    // float free_energy;		// Free energy change (dG) of dimer formation in cal/mol
    Primer *forward;			// Forward primer
    Primer *reverse;			// Reverse primer
    struct TPair *next;
};
typedef TPair Pair;

struct TTarget
{
	char *id;                   // Target sequence identification
	char *sequence;             // Target sequence
};
typedef TTarget Target;

// extern "C" Pair *design(char **targets, int m);
extern "C" Pair *design(Target *targets, int m);
// Pair *design(Target *targets, int m);
#endif
