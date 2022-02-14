#pragma once
/* #include "cell.h" */
#define PCRE2_CODE_UNIT_WIDTH 8
#include "atom.h"
#include "database/ads_database.h"
#include "lattice.h"
#include "molecule.h"
#include <pcre2.h>
#include <stdio.h>
#include <string.h>

/* Parse and construct the Adsorbate object */
Adsorbate *parse_adsorbate_from_file(char *fileName, char *name,
                                     AdsorbateInfo *adsInfo);
/* Parse and construct the Lattice object */
Lattice *parse_lattice_from_file(char *fileName, char *name);
