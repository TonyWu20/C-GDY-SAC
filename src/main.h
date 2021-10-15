#ifndef MAIN_H
#define MAIN_H
#include <ctype.h>
#define __GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
typedef struct
{
    double coord[3];
    char elm[2];
    int elmId;
    int itemId;
    int bCdSite; /* boolean for coordination site*/
    int bStem;   /* boolean for stem */
} ATOM_BLOCK;

typedef struct
{
    int atomNum;
    int CdSiteId;
    int StemId;
    ATOM_BLOCK molAtoms[];
} MOLECULE;

typedef struct
{
    double latVector[3][3];
    // char elements[MAX_ELEMENTS];
    // ATOM_BLOCK metalAtom;
    // ATOM_BLOCK carbonSites[5];
    double carbon_chain_vec[3];
    double carbon_metal_vec[3];
    int atomNum;
    ATOM_BLOCK totalAtoms[];
} BASE_LATTICE;
#endif
