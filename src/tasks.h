#pragma once
#include "assemble.h"
#include "atom.h"
#include "castep_database.h"
#include "cell.h"
#include "lattice.h"
#include "misc.h"
#include "molecule.h"
#include "parser.h"
enum
{
    ETHYLENE,
    ACETIC_ACID,
    ETHANOL,
    ETHANOL_OTHER
};

/* Perform all tasks in one go */
void allocateTasks(int pathNameCode, int *progress);
/* return list of ads file paths based on pathway name*/
char **pathway_adsLists(int pathNameCode, int *adsListLen);
/* Helper of pathway_adsLists */
char **adsListBuild(char *pathName, char *adsList[], int listLen);
