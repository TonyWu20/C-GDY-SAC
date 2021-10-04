#define PCRE2_CODE_UNIT_WIDTH 8
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include <string.h>
#include <pcre2.h>
typedef struct {
		double x,y,z;
		char *elm;
		int itemId;
} ATOM_BLOCK;
extern int atomCount;
long int countLines(FILE *file);
int returnLines(char *lineptr[], FILE *file);
void scanAtom(FILE *file);
//int reMatch(char *, PCRE2_SPTR , pcre2_match_data **, PCRE2_SIZE**);
