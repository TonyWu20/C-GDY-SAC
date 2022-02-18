#pragma once
#include "my_maths.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct _Atom
{
    char *element;
    int elementId;
    vec_double3 coord;
    int atomId;
    struct Atom_vtable *vtable;
};
/* Only the member: element needs free()
 */
typedef struct _Atom Atom;
struct Atom_vtable
{
    vec_double3 (*get_coord)(Atom *);
    Atom *(*dupAtom)(Atom *self);
    char *(*export_text)(Atom *self);
    void (*destroy)(Atom *);
};

/* @abstract: Create an atom by giving its element symbol,
 * elementId, coord, and atomId
 *
 * The passed element will be duplicated by strdup(),
 * needs free() when destroy.
 */
Atom *createAtom(char *element, vec_double3 coord, int atomId, int elementId);

/* @abstract: Duplicate an atom */
Atom *dupAtom(Atom *self);
void destroyAtom(Atom *);

// Methods
vec_double3 Atom_get_coord(Atom *);
char *Atom_textblock(Atom *self);
