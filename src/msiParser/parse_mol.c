#include "msiParser.h"
#include <stdio.h>
#include <stdlib.h>
#define OUT 0
#define INBLOCK 1
#define NEXT 2
static int reMatch(char *RegexStr, PCRE2_SPTR subject,
                   pcre2_match_data **match_data, PCRE2_SIZE **ovector);
static int atomBlockWalk(int blockFlag, char *line, ATOM_BLOCK *atom);
extern ATOM_BLOCK ads[]; /* defined in main.c */

/** scan atoms in .msi file
 * args:
 * 		FILE *file
 * process:
 * 		ATOM_BLOCK ads
 */
static int atomBlockWalk(int blockFlag, char *line, ATOM_BLOCK *atom)
{
    /* Declaration and initialization variables for pcre2 */
    pcre2_match_data *match_data = NULL; /* pointer for match_data */
    PCRE2_SIZE *ovector = NULL;          /* pointer for ovector */
    PCRE2_UCHAR8 *buffer; /* pointer to buffer to get substring */
    PCRE2_SIZE size;      /* pointer to store size of substring */

    char *RegexStr[] = {"ACL \"([0-9]{1,}) ([a-zA-Z]{1,})\"",
                        "XYZ \\(([0-9.-]{1,}) ([0-9.-]{1,}) ([0-9.-]{1,})\\)"};

    char cdStr[] = "# cd";
    if (blockFlag == OUT &&
        (reMatch(RegexStr[0], (PCRE2_SPTR)line, &match_data, &ovector) >
         1)) /* Out of atom blocks or start of new atom block*/
    {
        /* Get itemId */
        pcre2_substring_get_bynumber(match_data, 1, &buffer, &size);
        atom->itemId = atoi((const char *)buffer);

        /* Get element */
        pcre2_substring_get_bynumber(match_data, 2, &buffer, &size);
        atom->elm = strdup((const char *)buffer);
        /* Free match_data memory */
        pcre2_match_data_free(match_data);
        /* Check if it is the coordination site atom */
        if ((reMatch(cdStr, (PCRE2_SPTR)line, &match_data, &ovector)) == 1)
        {
            printf("Found cd_atom at Atom %d\n", atom->itemId);
            atom->bCdSite = 1;
        }
        else
        {
            atom->bCdSite = 0;
        }
        return INBLOCK; /* inside the atom block now */
    }
    else if (blockFlag == INBLOCK && (reMatch(RegexStr[1], (PCRE2_SPTR)line,
                                              &match_data, &ovector) > 1))
    {
        /* Stored x, y, z*/
        for (int j = 0; j < 3; j++)
        {
            pcre2_substring_get_bynumber(match_data, j + 1, &buffer, &size);
            atom->coord[j] = atof((const char *)buffer);
        }
        pcre2_match_data_free(match_data);
        return NEXT; /* Complete a block, now next */
    }
    else
    {
        pcre2_match_data_free(match_data); /* free memory */
        return blockFlag; /* keep the state (OUT/INBLOCK) if the above condition
                             does not fully meet */
    }
}

int scanAtom(FILE *file, ATOM_BLOCK *atoms)
{
    char line[256];      /* buffer for fgets */
    int blockFlag = OUT; /* indicate whether in the same block of atom */
    int atomCount = 0;
    while (fgets(line, sizeof(line), file)) /* get lines one by one */
    {
        /* assign blockFlag to the return value itself */
        blockFlag = atomBlockWalk(blockFlag, line, &atoms[atomCount]);
        switch (blockFlag)
        {
        case (NEXT):         /* Next block*/
            blockFlag = OUT; /* Change to out of block */
            atomCount++;     /* next pointer */
            break;
        default: /* blockFlag == OUT|INBLOCK, keep it */
            break;
        }
    }
    return atomCount;
}

void resetXYZ(int atomCount, ATOM_BLOCK *atoms)
{
    int cd_atom;
    for (cd_atom = 0; cd_atom < atomCount && atoms[cd_atom].bCdSite == 0;
         ++cd_atom)
    {
        ;
    }
    if (!atoms[cd_atom].bCdSite)
    {
        printf("No cd_atom is commented. Default to first atom.\n");
        cd_atom = 0;
        atoms[cd_atom].bCdSite = 1;
    }
    for (int i = 0; i < atomCount; i++)
    {
        /* Set the cd_atom as the origin */
        for (int j = 0; j < 3; j++)
        {
            atoms[i].coord[j] -= atoms[cd_atom].coord[j];
        }
    }
    printf("The %d atom is coord site\n", atoms[cd_atom].itemId);
}

/**We want to direct the pointer of *match_data and *ovector
 * to the pointers created by pcre2 functions, so we need to pass
 * the address of our defined pointers into the function to direct them
 * to the pointers created by pcre2 functions.*/
static int reMatch(char *RegexStr, PCRE2_SPTR subject,
                   pcre2_match_data **match_data, PCRE2_SIZE **ovector)
{
    int errornumber;
    PCRE2_SIZE erroroffset;
    int rc;
    pcre2_code *re;
    PCRE2_SPTR pattern;
    pattern = (PCRE2_SPTR)RegexStr;
    re = pcre2_compile(pattern, PCRE2_ZERO_TERMINATED, 0, &errornumber,
                       &erroroffset, NULL);
    if (re == NULL) /* Failed*/
    {
        PCRE2_UCHAR buffer[256];
        pcre2_get_error_message(errornumber, buffer, sizeof(buffer));
        printf("PCRE2 compilation failed at offset %d: %s\n", (int)erroroffset,
               buffer);
    }
    /* direct the true match_data to the passed pointer to *match_data */
    *match_data = pcre2_match_data_create_from_pattern(re, NULL);
    rc = pcre2_match(re, subject, (PCRE2_SIZE)strlen((const char *)subject), 0,
                     0, *match_data, NULL);
    if (rc > 1)
    {
        /* direct the true ovector to the passed pointer to *ovector */
        *ovector = pcre2_get_ovector_pointer(*match_data);
    }
    pcre2_code_free(re);
    return rc;
}
