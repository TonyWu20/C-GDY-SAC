#include "main.h"
#include <stdio.h>
#include <stdlib.h>
static int reMatch(char *RegexStr, PCRE2_SPTR subject,
                   pcre2_match_data **match_data, PCRE2_SIZE **ovector);
int atomCount = 0;       /* extern variable defined in main.h */
extern ATOM_BLOCK ads[]; /* defined in main.c */

/** scan atoms in .msi file
 * args:
 * 		FILE *file
 * process:
 * 		ATOM_BLOCK ads
 */
void scanAtom(FILE *file)
{
    char line[256];    /* buffer for fgets */
    int blockFlag = 0; /* indicate whether in the same block of atom */
    while (fgets(line, sizeof(line), file)) /* get lines one by one */
    {
        /* Declaration and initialization variables for pcre2 */
        pcre2_match_data *match_data = NULL; /* pointer for match_data */
        PCRE2_SIZE *ovector = NULL;          /* pointer for ovector */
        PCRE2_UCHAR8 *buffer; /* pointer to buffer to get substring */
        PCRE2_SIZE size;      /* pointer to store size of substring */

        /* Pattern to match: 1. ItemId and element; 2. XYZ coords */
        char *RegexStr[] = {
            "ACL \"([0-9]{1,}) ([a-zA-Z]{1,})\"",
            "XYZ \\(([0-9.-]{1,}) ([0-9.-]{1,}) ([0-9.-]{1,})\\)"};

        char cdStr[] = "# cd";
        /* Walkthrough lines */
        if (blockFlag == 0) /* Out of atom blocks or start of new atom block*/
        {
            int rc =
                reMatch(RegexStr[0], (PCRE2_SPTR)line, &match_data, &ovector);
            if (rc > 1) /* Matched */
            {
                /* Get itemId */
                pcre2_substring_get_bynumber(match_data, 1, &buffer, &size);
                ads[atomCount].itemId = atoi((const char *)buffer);

                /* Get element */
                pcre2_substring_get_bynumber(match_data, 2, &buffer, &size);
                ads[atomCount].elm = strdup((const char *)buffer);
                /* Free match_data memory */
                pcre2_match_data_free(match_data);
                /* Check if it is the coordination site atom */
                if ((rc = reMatch(cdStr, (PCRE2_SPTR)line, &match_data,
                                  &ovector)) == 1)
                {
                    printf("Found cd_atom at Atom %d\n", ads[atomCount].itemId);
                    ads[atomCount].bCdSite = 1;
                }
                else
                {
                    ads[atomCount].bCdSite = 0;
                }
                blockFlag = 1; /* inside the atom block now */
            }
            else
            {
                pcre2_match_data_free(match_data); /* free memory */
            }
        }
        else if (blockFlag == 1)
        {
            int rc =
                reMatch(RegexStr[1], (PCRE2_SPTR8)line, &match_data, &ovector);
            if (rc > 1) /* Matched */
            {
                /* Stored x, y, z*/
                pcre2_substring_get_bynumber(match_data, 1, &buffer, &size);
                ads[atomCount].x = atof((const char *)buffer);
                pcre2_substring_get_bynumber(match_data, 2, &buffer, &size);
                ads[atomCount].y = atof((const char *)buffer);
                pcre2_substring_get_bynumber(match_data, 3, &buffer, &size);
                ads[atomCount].z = atof((const char *)buffer);
                pcre2_match_data_free(match_data);
                blockFlag = 0; /* Complete a block, now next */
                atomCount++;   /* Complete a block, now next */
            }
            else
            {
                pcre2_match_data_free(match_data); /* free memory */
            }
        }
    }
}

void resetXYZ()
{
    int cd_atom;
    for (cd_atom = 0; cd_atom < atomCount && ads[cd_atom].bCdSite == 0;
         ++cd_atom)
    {
        ;
    }
    if (!ads[cd_atom].bCdSite)
    {
        printf("No cd_atom is commented. Default to first atom.\n");
        cd_atom = 0;
        ads[cd_atom].bCdSite = 1;
    }
    for (int i = 0; i < atomCount; i++)
    {
        /* Set the cd_atom as the origin */
        ads[i].x -= ads[cd_atom].x;
        ads[i].y -= ads[cd_atom].y;
        ads[i].z -= ads[cd_atom].z;
    }
    printf("The %d atom is coord site\n", ads[cd_atom].itemId);
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
