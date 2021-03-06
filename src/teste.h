struct TPrimer {
    int targetId;			// Target identification
    char *sequence;			// Primer sequence
    int length;				// Primer length
    int position;			// Begin position in target
    bool forward;			// Indicate whether primer is forward
    bool reverse;			// Indicate whether primer is reverse
    float GCContent;			// GC content
    float temperature;			// Melting temperature in C
    float freeEnergy;			// Free energy change (dG) of hairpin formation in cal/mol 
    // int number;			// Number of non target sequences
    float *forwardElongationEfficiency; // Forward elongation efficiency
    float *reverseElongationEfficiency; // Reverse elongation efficiency
    struct TPrimer *next;
};
typedef TPrimer Primer;

struct TPair {
    int targetId;			// Target identification 
    int productSize;			// PCR Product size 
    float freeEnergy;			// Free energy change (dG) of dimer formation in cal/mol
    Primer *forward;			// Forward primer
    Primer *reverse;			// Reverse primer
    struct TPair *next;
};
typedef TPair Pair;

extern "C" double getK(float dG, float T); 
extern "C" Pair *design(char **targets, int m); 
