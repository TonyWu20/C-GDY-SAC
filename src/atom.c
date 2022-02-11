#include "atom.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Implementation of Atom methods

struct Atom_vtable atom_vtable = {
    Atom_get_coord,  Atom_update_coord, Atom_get_atomId,
    Atom_set_atomId, Atom_get_treeId,   Atom_set_treeId,
    dupAtom,         Atom_textblock,    destroyAtom};
Atom *createAtom(char *element, char *label, double *coord, int atomId,
                 int treeId)
{
    Atom *newAtom = malloc(sizeof(Atom));
    newAtom->element = strdup(element);
    newAtom->ACL_Label = strdup(label);
    newAtom->coord = malloc(sizeof(double) * 4);
    memcpy(newAtom->coord, coord, sizeof(double) * 4);
    newAtom->atomId = atomId;
    newAtom->treeId = treeId;
    newAtom->vtable = &atom_vtable;
    return newAtom;
}

Atom *dupAtom(Atom *self)
{
    Atom *dup = createAtom(self->element, self->ACL_Label, self->coord,
                           self->atomId, self->treeId);
    return dup;
}
void destroyAtom(Atom *atomPtr)
{
    free(atomPtr->element);
    free(atomPtr->ACL_Label);
    free(atomPtr->coord);
    free(atomPtr);
}

double *Atom_get_coord(Atom *self)
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

int Atom_get_treeId(Atom *self)
{
    return self->treeId;
}

void Atom_set_treeId(Atom *self, int newId)
{
    self->atomId = newId;
}

char *Atom_textblock(Atom *self)
{
    int textLen = snprintf(
        NULL, 0,
        "  (%d Atom\n    (A C ACL \"%s\")\n    (A C Label \"%s\")\n    (A "
        "D XYZ (%.12f %.12f %.12f))\n    (A I Id %d)\n  )\n",
        self->treeId, self->ACL_Label, self->element, self->coord[0],
        self->coord[1], self->coord[2], self->atomId);
    textLen += 1;
    char *buffer = malloc(textLen); // freed after export to msi
    snprintf(buffer, textLen,
             "  (%d Atom\n    (A C ACL \"%s\")\n    (A C Label \"%s\")\n    (A "
             "D XYZ (%.12f %.12f %.12f))\n    (A I Id %d)\n  )\n",
             self->treeId, self->ACL_Label, self->element, self->coord[0],
             self->coord[1], self->coord[2], self->atomId);
    buffer = realloc(buffer, strlen(buffer) + 1);
    return buffer;
}
// End of Atom
