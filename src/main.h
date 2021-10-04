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
        int bCdSite; /* boolean for coordination site*/
} ATOM_BLOCK;
extern int atomCount;
long int countLines(FILE *file);
int returnLines(char *lineptr[], FILE *file);
/* parse_mol */
void scanAtom(FILE *file);
void resetXYZ(void);

/* parse_base */
//int reMatch(char *, PCRE2_SPTR , pcre2_match_data **, PCRE2_SIZE**);
