#include "atom.h"

// Implementation of Atom methods

struct Atom_vtable atom_vtable = {Atom_get_coord, dupAtom, Atom_textblock,
                                  destroyAtom};
Atom *createAtom(char *element, int elementId, vec_double3 coord, int atomId)
{
    Atom *newAtom = malloc(sizeof(Atom));
    newAtom->element = strdup(element);
    newAtom->elementId = elementId;
    newAtom->coord = coord;
    newAtom->atomId = atomId;
    newAtom->vtable = &atom_vtable;
    return newAtom;
}
Atom *dupAtom(Atom *self)
{
    Atom *dup =
        createAtom(self->element, self->elementId, self->coord, self->atomId);
    return dup;
}
void destroyAtom(Atom *atomPtr)
{
    free(atomPtr->element);
    free(atomPtr);
}

vec_double3 Atom_get_coord(Atom *self)
{
    return self->coord;
}

char *Atom_textblock(Atom *self)
{
    char *buffer;
    asprintf(
        &buffer,
        "  (%d Atom\n    (A C ACL \"%d %s\")\n    (A C Label \"%s\")\n    (A "
        "D XYZ (%.12f %.12f %.12f))\n    (A I Id %d)\n  )\n",
        self->atomId + 1, self->elementId, self->element, self->element,
        self->coord.x, self->coord.y, self->coord.z, self->atomId);
    return buffer;
}
// End of Atom
