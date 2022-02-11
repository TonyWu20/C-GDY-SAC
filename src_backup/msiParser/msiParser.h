#define PCRE2_CODE_UNIT_WIDTH 8
#include "../main.h"
#include <pcre2.h>
#define MAX_ELEMENTS 4
#define MAXLINE 128
#define OUT 0
#define INBLOCK 1
#define GETXYZ 2
#define NEXT 3
#define LATTICE 10
#define ATOMS 11

/* basic pcre2 match */
pcre2_code *init_re(char *RegexStr);

/* basic walkthrough .msi */
int atomBlockWalk(int blockFlag, char *line, ATOM_BLOCK *atom);
int scanAtom(FILE *file, ATOM_BLOCK *atom);
void resetXYZ(int atomCount, ATOM_BLOCK *atoms);
/* parse_mol */
MOLECULE *init_mol(int atomNum);
int find_CdAtom(MOLECULE *);
int find_Stem(MOLECULE *);
MOLECULE *parseMol(FILE *file);
void load_StemVector(MOLECULE *, double *);
void init_xz_plane(MOLECULE *); /* init the stem direction and turn molecule
                                   from xy plane to xz plane */

/* parse_base */
int countAtoms(FILE *file);
int get_LatVector(FILE *file, double (*v)[3]);
BASE_LATTICE *init_lattice(int atomNum);
BASE_LATTICE *parseBase(FILE *file);
void assign_carbonVector(BASE_LATTICE *lat, int bAtomId, double *vec);

/* mod_msi */
BASE_LATTICE *init_adsorbed_lat(BASE_LATTICE *, MOLECULE *);
void appendMolAtoms(BASE_LATTICE *source, MOLECULE *add_mol,
                    BASE_LATTICE *target);
void placeMol(MOLECULE *mol, BASE_LATTICE *lat, int destId,
              BASE_LATTICE *target);
void rotMol(MOLECULE *, double angle, char axis);
void get_CoordMat(MOLECULE *s,
                  double[][3]); /* define the n-1 dimension when you want to
                                   pass the array as a pointer */
void assignCoordtoMol(double[][3], MOLECULE *);
void align_carbon_chain(MOLECULE *, double *chain_vec);
void attach_carbon_chain(MOLECULE *, BASE_LATTICE *, int carbonId);

/* basic functions */

int saveItemId(char *, ATOM_BLOCK *);
int saveElmInfo(char *, ATOM_BLOCK *);
int checkCdSite(char *);
int saveCoord(char *, ATOM_BLOCK *);
int checkStem(char *);
