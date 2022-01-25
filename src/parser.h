#pragma once
#define PCRE2_CODE_UNIT_WIDTH 8
#include "atom.h"
#include "lattice.h"
#include "molecule.h"
#include <pcre2.h>
#include <stdio.h>
#include <string.h>

/* Initialize re pattern */
pcre2_code *init_re(char *RegexStr);

/* Wrapper for one-time re match */
void re_match(pcre2_code *re_pattern, pcre2_match_data **match_data, int *rc,
              char *subject);

/* Match and parse all atoms
 * Returns: array of Atom *
 */
Atom **get_all_atoms(char *text, int *returnSize);

/* Parse number of coordinate atoms */
void get_cd_num(char *subject, int *cd_num);

/* Parse atom Ids of coordinate atoms and store in **cd_ids
 */
void get_cd_info(char *subject, int *cd_num, int **cd_ids);

/* Get if the molecule is symmetric */
void get_symmetric_info(char *subject, int *bSym);

/* Store the stem atomId in the pointer */
void get_stem_arr(char *subject, int *stem_arr);

/* Store the plane atomId in the pointer */
void get_plane_arr(char *subject, int *plane_arr);

/* Parse the lattice vectors */
void get_lattice_vectors(char *subject, Matrix **result);

/* Parse single atom block */
Atom *parse_atom(char *atom_block);

/* Parse and construct the Adsorbate object */
Adsorbate *parse_molecule_from_file(char *fileName, char *name);

/* Parse and construct the Lattice object */
Lattice *parse_lattice_from_file(char *fileName, char *name);
