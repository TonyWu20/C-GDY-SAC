#pragma once
#define PCRE2_CODE_UNIT_WIDTH 8
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
char **atom_block(char *text, int *returnSize);
