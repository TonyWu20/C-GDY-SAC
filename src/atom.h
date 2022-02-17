#pragma once
#include <simd/simd.h>
struct _Atom
{
    char *element;
    int elementId;
    simd_double3 coord;
    int atomId;
    struct Atom_vtable *vtable;
};
typedef struct _Atom Atom;
struct Atom_vtable
{
    simd_double3 (*get_coord)(Atom *);
    int (*get_atomId)(Atom *);
    void (*set_atomId)(Atom *, int);
    Atom *(*dupAtom)(Atom *self);
    char *(*export_text)(Atom *self);
    void (*destroy)(Atom *);
};

// Memory management

/* Create an Atom struct with element name, coord in simd_double3 and atomId */
Atom *createAtom(char *element, simd_double3 coord, int atomId, int elementId);

Atom *dupAtom(Atom *self);
void destroyAtom(Atom *);

// Methods
simd_double3 Atom_get_coord(Atom *);
void Atom_update_coord(Atom *, double x, double y, double z);
int Atom_get_atomId(Atom *);
void Atom_set_atomId(Atom *, int);
char *Atom_textblock(Atom *self);
