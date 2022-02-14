#include "atom.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Implementation of Atom methods

struct Atom_vtable atom_vtable = {
    Atom_get_coord, Atom_update_coord, Atom_get_atomId, Atom_set_atomId,
    dupAtom,        Atom_textblock,    destroyAtom};
Atom *createAtom(char *element, simd_double3 coord, int atomId)
{
    Atom *newAtom = malloc(sizeof(Atom));
    strncpy(newAtom->element, element, strlen(element));
    newAtom->coord = coord.xyz;
    newAtom->atomId = atomId;
    newAtom->vtable = &atom_vtable;
    return newAtom;
}

Atom *dupAtom(Atom *self)
{
    Atom *dup = createAtom(self->element, self->coord, self->atomId);
    return dup;
}
void destroyAtom(Atom *atomPtr)
{
    free(atomPtr);
}

simd_double3 Atom_get_coord(Atom *self)
{
    return self->coord;
}

void Atom_update_coord(Atom *self, double x, double y, double z)
{
    self->coord[0] = x;
    self->coord[1] = y;
    self->coord[2] = z;
}

int Atom_get_atomId(Atom *self)
{
    return self->atomId;
}

void Atom_set_atomId(Atom *self, int newId)
{
    self->atomId = newId;
}

void Atom_set_treeId(Atom *self, int newId)
{
    self->atomId = newId;
}

char *Atom_textblock(Atom *self)
{
    int textLen =
        1 + snprintf(NULL, 0,
                     "  (%d Atom\n    (A C ACL \"%d %s\")\n    (A C Label "
                     "\"%s\")\n    (A "
                     "D XYZ (%.12f %.12f %.12f))\n    (A I Id %d)\n  )\n",
                     self->atomId + 1, self->elementId, self->element,
                     self->element, self->coord.x, self->coord.y, self->coord.z,
                     self->atomId);
    char *buffer = malloc(textLen); // freed after export to msi
    snprintf(
        buffer, textLen,
        "  (%d Atom\n    (A C ACL \"%d %s\")\n    (A C Label \"%s\")\n    (A "
        "D XYZ (%.12f %.12f %.12f))\n    (A I Id %d)\n  )\n",
        self->atomId + 1, self->elementId, self->element, self->element,
        self->coord.x, self->coord.y, self->coord.z, self->atomId);
    return buffer;
}
// End of Atom
