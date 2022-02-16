#pragma once
#include "cell.h"
#include "database/database.h"
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
/* Adjust necessary parameters according to the constitutional elements
 *, write to target directory
 */
void write_param(Cell *self);

/* Write kptaux content out to the target directory */
void write_kptaux(Cell *self);

void write_trjaux(Cell *self);

void copy_potentials(Cell *self);

void write_pbsScript(Cell *self);

void write_SMCastepExtension(Cell *self);
