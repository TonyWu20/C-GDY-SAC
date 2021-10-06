#include "../main.h"
#define MAX_ELEMENTS 4
typedef struct {
		double coord[3];
		char *elm;
		int itemId;
        int bCdSite; /* boolean for coordination site*/
} ATOM_BLOCK;

typedef struct
{
		char latLines[8];
		double latVector[3][3];
		char elements[MAX_ELEMENTS];
		ATOM_BLOCK metalAtom;
		ATOM_BLOCK carbonSites[5];
		ATOM_BLOCK totalAtoms[];
} BASE_LATTICE;

/* parse_mol */
int scanAtom(FILE *file, ATOM_BLOCK *atom);
void resetXYZ(int atomCount, ATOM_BLOCK *atoms);
BASE_LATTICE parseBase(FILE *file);
