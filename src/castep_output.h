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
    struct Cell_vtable *vtable;
    struct Cell_textFunc *textTable;
    void (*destroy)(Cell *self);
};

struct Cell_vtable
{
    void (*sortAtoms)(Cell *self);
    char **(*sortElmList)(Cell *self, int *returnSize);
};

struct Cell_textFunc
{
    char *(*blockWriter)(Cell *self, char *blockName,
                         char *(*blockTextWriter)(Cell *self));
};

Cell *createCell(Lattice *lat, CastepInfo **table);
void cellLoadInfoTab(Cell *self, CastepInfo **table);
void destroyCell(Cell *self);
char *cellWriteBlock(Cell *self, char *blockName,
                     char *(*blockTextWriter)(Cell *self));

char *cell_latticeVector_writer(Cell *self);
char *cell_fracCoord_writer(Cell *self);
char *cell_kPointsList_writer(Cell *self);
char *cell_miscOptions_writer(Cell *self);
char *cell_ionicConstraints_writer(Cell *self);
char *cell_externalPressure_writer(Cell *self);
char *cell_speciesMass_writer(Cell *self);
char *cell_speciesPot_writer(Cell *self);
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

struct Cell_vtable cellVTable = {sortAtomsByElement, sortedElementList};
struct Cell_textFunc cellTextTable = {cellWriteBlock};
