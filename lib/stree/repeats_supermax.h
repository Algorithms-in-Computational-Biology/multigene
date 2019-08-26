#ifndef STRMAT
#define STRMAT
#endif

#ifndef _REPEATS_SUPERMAX_H_
#define _REPEATS_SUPERMAX_H_

#include "stree_strmat.h"

typedef struct supermax_node {
  char *S;
  int M, num_witness, num_leaves, percent;
  struct supermax_node *next;
  STREE_NODE node;
} STRUCT_SUPERMAX, *SUPERMAXIMALS;

SUPERMAXIMALS supermax_find(char *S, int M, int min_percent, int min_length);

SUPERMAXIMALS supermax_find(SUFFIX_TREE tree, int min_percent, int min_length);

#endif
