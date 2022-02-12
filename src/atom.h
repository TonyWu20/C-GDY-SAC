#pragma once
#include <simd/simd.h>
struct _Atom
{
    char *element;
    char *ACL_Label;
    int elementId;
    simd_double3 coord;
    int atomId;
    struct Atom_vtable *vtable;
};
typedef struct _Atom Atom;
struct Atom_vtable
{
    simd_double3 (*get_coord)(Atom *);
    void (*update_coord)(Atom *, double, double, double);
    int (*get_atomId)(Atom *);
    void (*set_atomId)(Atom *, int);
    Atom *(*dupAtom)(Atom *self);
    char *(*export_text)(Atom *self);
    void (*destroy)(Atom *);
};

// Memory management
Atom *createAtom(char *element, char *label, simd_double3 coord, int atomId);

Atom *dupAtom(Atom *self);
void destroyAtom(Atom *);

// Methods
simd_double3 Atom_get_coord(Atom *);
void Atom_update_coord(Atom *, double x, double y, double z);
int Atom_get_atomId(Atom *);
void Atom_set_atomId(Atom *, int);
char *Atom_textblock(Atom *self);
