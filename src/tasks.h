#pragma once
#include "assemble.h"
#include "atom.h"
#include "castep_database.h"
#include "cell.h"
#include "lattice.h"
#include "misc.h"
#include "molecule.h"
#include "parser.h"

/* Perform all tasks in one go */
void allocateTasks(char *pathName);
/* return list of ads file paths */
char **pathway_adsLists(char *pathName, int *adsListLen);
