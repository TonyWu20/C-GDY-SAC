#include "../main.h"
#define MAX_ELEMENTS 4
typedef struct {
		double coord[3];
		char *elm;
		int elmId;
		int itemId;
        int bCdSite; /* boolean for coordination site*/
} ATOM_BLOCK;

typedef struct
{
		double latVector[3][3];
		//char elements[MAX_ELEMENTS];
		//ATOM_BLOCK metalAtom;
		//ATOM_BLOCK carbonSites[5];
		int atomNum;
		ATOM_BLOCK totalAtoms[];
} BASE_LATTICE;

/* parse_mol */
int scanAtom(FILE *file, ATOM_BLOCK *atom);
int countAtoms(FILE *file);
void resetXYZ(int atomCount, ATOM_BLOCK *atoms);
BASE_LATTICE* init_lattice(int atomNum);
BASE_LATTICE parseBase(FILE *file);
void appendAtoms(BASE_LATTICE source, ATOM_BLOCK add_atoms, BASE_LATTICE *target);

/* basic functions */

int saveItemId(char *, ATOM_BLOCK *);
int saveElmInfo(char *, ATOM_BLOCK *);
int checkCdSite(char *);
int saveCoord(char *, ATOM_BLOCK *);

int get_LatVector(FILE *file, double (*v)[3]);
