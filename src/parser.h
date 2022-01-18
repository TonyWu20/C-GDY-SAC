#pragma once
#define PCRE2_CODE_UNIT_WIDTH 8
#include "atom.h"
#include <pcre2.h>
#include <stdio.h>
#include <string.h>

#define MAXLINE 1000
enum
{
    OUT,
    INBLOCK,
    GETXYZ,
    NEXT
};

pcre2_code *init_re(char *RegexStr);
Atom **atom_block(char *text, int *returnSize);
Atom *parse_atom(char *atom_block);
