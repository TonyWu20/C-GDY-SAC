#pragma once
#include "castep_database.h"
#include "cell.h"
#include "lattice.h"
#include "misc.h"
#include "my_maths.h"
#include "parser.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int getFinalCutoffEnergy(Cell *cell);

/* Parse fine cutoff energy from given potential file */
int parse_fineCutoffEnergy(const char *fileName);
