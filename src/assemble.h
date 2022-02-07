#pragma once
#include "atom.h"
#include "lattice.h"
#include "molecule.h"
#define NULLSITE 255

/* Append cd2 mol name to lattice name
 * Returns malloced string
 */
char *append_mol_name(Lattice *lat, Adsorbate *ads, int c1, int c2);

/* Align the adsorbate with the target sites */
int init_ads_direction(Lattice *lat, Adsorbate *ads, int id_from, int id_to);

/* Move the adsorbate to the target position (carbon chain ver)
 * The carbonId 1 and 2 must be the adjacent carbon atoms
 */
Lattice *Add_mol_to_lattice(Lattice *, Adsorbate *, int, int, char *pathName,
                            double height);
/* Move the adsorbate to the carbon metal position */
