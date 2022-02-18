#pragma once
#include "assemble.h"
#include "atom.h"
#include "cell.h"
#include "database/ads_database.h"
#include "database/database.h"
#include "database/lattice_database.h"
#include "lattice.h"
#include "molecule.h"
#include "parser.h"
#include <omp.h>
enum
{
    ETHYLENE,
    ACETIC_ACID,
    ETHANOL,
    ETHANOL_OTHER
};

/* Perform all tasks in one go */
void allocateTasks(ElmTableYAML *elmTableYAML, char **elements,
                   AdsTableYAML *AdsTableYAML, int *progress);
/* return list of ads file paths based on pathway name*/
char **pathway_adsLists(int pathNameCode, int *adsListLen);
/* Helper of pathway_adsLists */
char **adsListBuild(char *pathName, char *adsList[], int listLen);
void generator(ElmTableYAML *elmTableYAML, char **elements,
               AdsTableYAML *adsTableYAML, int *progress);
