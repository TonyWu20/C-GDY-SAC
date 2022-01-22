#pragma once
#define PCRE2_CODE_UNIT_WIDTH 8
#include "atom.h"
#include "lattice.h"
#include "molecule.h"
#include <pcre2.h>
#include <stdio.h>
#include <string.h>

pcre2_code *init_re(char *RegexStr);
Atom **get_all_atoms(char *text, int *returnSize);
void get_cd_num(char *subject, int *cd_num);
void get_cd_info(char *subject, int *cd_num, int **cd_ids);
void get_stem_arr(char *subject, int *stem_arr);
void get_plane_arr(char *subject, int *plane_arr);
void get_lattice_vectors(char *subject, Matrix **result);
Atom *parse_atom(char *atom_block);
Adsorbate *parse_molecule_from_file(char *fileName, char *name);
Lattice *parse_lattice_from_file(char *fileName, char *name);
