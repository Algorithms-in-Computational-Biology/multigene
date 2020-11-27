#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "teste.h"

#define R 0.00199

extern "C" double getK(float dG, float T) 
{
	return exp(-dG/R * T);
}

double getKeff(double Kpt, double Kpd, double Kpf, double Ktf, float P) 
{
	double KpdP = Kpd * P;
	double sqr = sqrt(pow(Kpf + 1, 2) * 8 * KpdP);

	return (4 * Kpt * KpdP)/(-1 - Kpf + sqr) * (1 + Ktf);
}

double getEffhyb(float P, double Keff) 
{
	double PKeff = P * Keff; 

	return PKeff/(1 + PKeff);
}

double getEffamp(double Effhyd, double Effelong) 
{
	return Effhyd * Effelong;
}

float* newArray(int size) {
    float *array = (float*) malloc(size * sizeof(float));
    if (array == NULL) {
        printf("Memoria insuficiente!");

        exit(0);
    }
    return array;
}

char* newString(int length) {
   char* str = (char*) malloc((length + 1) * sizeof(char));
   if (str == NULL) {
        printf("Memoria insuficiente!");

        exit(0);
   } 
   return str;
}

Primer* newPrimer(char* sequence, int length, int m) {
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

Primer* copyPrimer(Primer *primer, int m) {
    Primer* copy = newPrimer(primer->sequence, primer->length, m);
    copy->targetId = primer->targetId;
    copy->position = primer->position;
    copy->forward = primer->forward;
    copy->reverse = primer->reverse;
    copy->GCContent = primer->GCContent;
    copy->temperature = primer->temperature;
    copy->freeEnergy = primer->freeEnergy;
    memcpy(copy->forwardElongationEfficiency, primer->forwardElongationEfficiency, m * sizeof(float));
    memcpy(copy->reverseElongationEfficiency, primer->reverseElongationEfficiency, m * sizeof(float));

    return copy;
}

void addPrimer(Primer **head, Primer *primer) {
    primer->next = *head;
    *head = primer;
}

Pair *newPair(Primer *forward, Primer *reverse, int m)
//Pair *newPair() 
{
    Pair *pair = (Pair*) malloc(sizeof(Pair));
    if (pair == NULL) 
    {
        printf("Memoria insuficiente!");

        exit(0);
    }
    pair->forward = copyPrimer(forward, m);
    pair->reverse = copyPrimer(reverse, m);
    pair->next = NULL;

    return pair;
}

void addPair(Pair **head, Pair *pair) 
{
    pair->next = *head;
    *head = pair;
}

extern "C" Pair *design(char **targets, int m) 
{
	for(int i = 0; i < m; i++) 
	{
		//printf("%s %d\n", targets[i], i);
	}
	Primer *forward = newPrimer(targets[0], 1, m);
    forward->position = 1;
    forward->targetId = 1;
	forward->forward = true;
	forward->reverse = false;
	forward->GCContent = 1.1;
	forward->temperature = 0.5;
	forward->freeEnergy = 0.5;
	forward->forwardElongationEfficiency[0] = 0.15;
	forward->reverseElongationEfficiency[0] = 0.25;

	Primer *reverse = newPrimer(targets[1], 1, m);
    reverse->position = 2;
    reverse->targetId = 1;
	reverse->forward = false;
	reverse->reverse = true;
	reverse->GCContent = 1.2;
	reverse->temperature = 1.5;
	reverse->freeEnergy = -0.5;
	reverse->forwardElongationEfficiency[0] = 0.35;
	reverse->reverseElongationEfficiency[0] = 0.45;

	Pair *pairs = NULL;
	Pair *pair = newPair(forward, reverse, m);

	pair->targetId = 1;
	pair->productSize = 2;
	pair->freeEnergy = -0.5f;
	
	addPair(&pairs, pair);

	return pairs;
}


/*int main() 
{
 	int m = 2;
	char *targets[] = {"A","B"};
	Pair *pairs = design(targets, m - 1);
        for (Pair *pair = pairs; pair != NULL; pair = pair->next) 
        {
	    printf("%d %d %.2f ", pair->targetId, pair->productSize, pair->freeEnergy);
	    Primer *forward = pair->forward;	
	    Primer *reverse = pair->reverse;
	    printf("%s %d %s %d\n", forward->sequence, forward->position, reverse->sequence, reverse->position);
        }
	
	//printf("%lf\n", getK(28.5, 154.6));
	return 0;
}*/
