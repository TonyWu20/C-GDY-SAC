#pragma once
#include "atom.h"
#include "castep_database.h"
#include "lattice.h"
#include "molecule.h"
#include "uthash.h"
#include <stdbool.h>
#include <stdio.h>

/* Add CastepInfo item to hash table
 */

typedef struct Cell Cell;
struct Cell
{
    Lattice *lattice;
    CastepInfo *infoTab;
    bool atomSorted;
    int elmNums;
    char **elmLists;
    struct Cell_vtable *vtable;
    struct Cell_textFunc *textTable;
    void (*destroy)(Cell *self);
};

/* Virtual functions table for Cell */
struct Cell_vtable
{
    void (*sortAtoms)(Cell *self);
    char **(*sortElmList)(Cell *self, int *returnSize);
    void (*exportCell)(Cell *self);
};

/* Virtual functions table for dealing with text output for Cell */
struct Cell_textFunc
{
    char *(*blockWriter)(Cell *self, char *blockName,
                         char *(*blockTextWriter)(Cell *self));
};

/* Initialize a Cell object and construct a new hashtable for the existing
 * elements in the cell
 * Will automatically sort the atoms in atomic number (1st) and atom Id
 * (2nd) order
 */
Cell *createCell(Lattice *lat, CastepInfo *table);
/* Create a hashtable for the exisiting elements in cell
 */
void cellLoadInfoTab(Cell *self, CastepInfo *table);
/* Destroy and free memory of Cell object */
void destroyCell(Cell *self);
/* General method to produce BLOCK_XX in .cell
 */
char *cellWriteBlock(Cell *self, char *blockName,
                     char *(*blockTextWriter)(Cell *self));

void cellExport(Cell *self);
/* Methods to fill BLOCK_XX
 */

/* LATTICE_CART */
char *cell_latticeVector_writer(Cell *self);
/* POSITIONS_FRAC */
char *cell_fracCoord_writer(Cell *self);

/* BS_KPOINTS_LIST or KPOINTS_LIST */
char *cell_kPointsList_writer(Cell *self);

/* Non BLOCK wrapped options
 * FIX_ALL_CELL and FIX_COM
 */
char *cell_miscOptions_writer(Cell *self);
/* IONIC_CONSTRAINTS (empty) */
char *cell_ionicConstraints_writer(Cell *self);
/* EXTERNAL_EFIELD */
char *cell_externalEfield_writer(Cell *self);
/* EXTERNAL_PRESSURE */
char *cell_externalPressure_writer(Cell *self);
/* SPECIES_MASS */
char *cell_speciesMass_writer(Cell *self);
/* SPECIES_POT */
char *cell_speciesPot_writer(Cell *self);
/* SPECIES_LCAO_STATES */
char *cell_speciesLCAOstates_writer(Cell *self);

/* Sort array of [Atom *] by element atomic number and atom id
 */
void sortAtomsByElement(Cell *self);

/* Remove duplicate element name and generate a set of existing elements
 * in the lattice
 * params:
 *      Cell *self: self object of Cell;
 *      int *returnSize: get the size of return char **;
 * returns:
 * A malloced array of pointers to strings
 */
char **sortedElementList(Cell *self, int *returnSize);

/* Other required seed files */
