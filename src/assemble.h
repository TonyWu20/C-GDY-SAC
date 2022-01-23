#pragma once
#include "atom.h"
#include "lattice.h"
#include "molecule.h"

// Align the adsorbate with the carbon chain
void ads_align_with_carbon_chain(Lattice *lat, Adsorbate *ads);
// Move the adsorbate to the target position
Lattice *Add_mol_to_carbon_chain(Lattice *, Adsorbate *, int, int);
