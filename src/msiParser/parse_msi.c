#include "msiParser.h"

/** scan atoms in .msi file
 * args:
 * 		FILE *file
 * 		ATOM_BLOCK *atom
 * returns:
 * 		int atomCount
 */
int scanAtom(FILE *file, ATOM_BLOCK *atoms)
{
    char line[MAXLINE];  /* buffer for fgets */
    int blockFlag = OUT; /* indicate whether in the same block of atom */
    ATOM_BLOCK *ptr;
    ptr = atoms;
    while (fgets(line, sizeof(line), file)) /* get lines one by one */
    {
        /* assign blockFlag to the return value itself */
        blockFlag = atomBlockWalk(blockFlag, line, ptr);
        switch (blockFlag)
        {
        case (NEXT):         /* Next block*/
            blockFlag = OUT; /* Change to out of block */
            ptr++;           /* next pointer */
            break;
        default: /* blockFlag == OUT|INBLOCK, keep it */
            break;
        }
    }
    int atomCount = ptr - atoms;
    return atomCount;
}

/** Found the commented coord site atom, and then reset coords to set that atom
 * as (0,0,0)
 *  Args: int atomCount; number of atoms
 *  	  ATOM_BLOCK *atoms; pointer to the array of ATOM_BLOCK to be modified
 * return: void;
 */
void resetXYZ(int atomCount, ATOM_BLOCK *atoms)
{
    int cd_atom;
    for (cd_atom = 0; cd_atom < atomCount && atoms[cd_atom].bCdSite == 0;
         ++cd_atom)
    {
        ; /* increment cd_atom until the atom with bCdSite = 1 is found*/
    }
    if (!atoms[cd_atom].bCdSite) /* Not found */
    {
        printf("No cd_atom is commented. Default to first atom.\n");
        cd_atom = 0;
        atoms[cd_atom].bCdSite = 1;
    }
    for (int i = 0; i < atomCount; i++) /* reset each atom coordinate */
    {
        if (i != cd_atom)
        {
            /* Set the cd_atom as the origin */
            for (int j = 0; j < 3; j++)
            {
                atoms[i].coord[j] -= atoms[cd_atom].coord[j];
            }
        }
    }

    for (int j = 0; j < 3; j++)
    {
        atoms[cd_atom].coord[j] -= atoms[cd_atom].coord[j];
    }
}

/**We want to direct the pointer of *match_data and *ovector
 * to the pointers created by pcre2 functions, so we need to pass
 * the address of our defined pointers into the function to direct them
 * to the pointers created by pcre2 functions.*/
int reMatch(char *RegexStr, PCRE2_SPTR subject, pcre2_match_data **match_data,
            PCRE2_SIZE **ovector)
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

/** Walk every line of .msi file to parse atoms
 *  args:
 *  	 int blockFlag: OUT, INBLOCK, NEXT
 *  	 char *line: string pointer of the line from fgets
 *  	 ATOM_BLOCK *atom: pointer dereferenced from ATOM_BLOCK atom[atomCount]
 *  returns:
 *  	 int blockFlag
 */
int atomBlockWalk(int blockFlag, char *line, ATOM_BLOCK *atom)
{
    /* Based on status indicator blockFlag */
    switch (blockFlag)
    {
    case OUT:
        blockFlag = saveItemId(line, atom); /* OUT or INBLOCK */
        break;
    case INBLOCK:
        blockFlag = saveElmInfo(line, atom); /* INBLOCK or GETXYZ */
        break;
    case GETXYZ:
        blockFlag = saveCoord(line, atom); /* GETXYZ or OUT */
        break;
    default:
        break;
    }
    return blockFlag; /* repeat status flag input in the while-loop in
                         scanAtom() */
}

/** Save Item Id when found. Return INBLOCK or OUT to feedback the control case.
 *  Args:
 *  	char *line;
 *  	ATOM_BLOCK *atom;
 *  returns:
 *  	int blockFlag = OUT (when not found) or INBLOCK (when found, next stage)
 */
int saveItemId(char *line, ATOM_BLOCK *atom)
{
    /* Declaration and initialization variables for pcre2 */
    pcre2_match_data *match_data = NULL; /* pointer for match_data */
    PCRE2_SIZE *ovector = NULL;          /* pointer for ovector */
    PCRE2_UCHAR8 *buffer; /* pointer to buffer to get substring */
    PCRE2_SIZE size;      /* pointer to store size of substring */
    char *RegexStr = "([0-9]{1,}) Atom";
    if (reMatch(RegexStr, (PCRE2_SPTR8)line, &match_data, &ovector) == 2)
    {
        pcre2_substring_get_bynumber_8(match_data, 1, &buffer, &size);
        atom->itemId = atoi((const char *)buffer);
        pcre2_substring_free(buffer);
        pcre2_match_data_free(match_data);
        return INBLOCK;
    }
    pcre2_match_data_free(match_data);
    return OUT;
}

/** Save Element name and element id when found. Return INBLOCK or GETXYZ to
 * feedback the control case.
 * Args:
 * 		char *line;
 * 		ATOM_BLOCK *atom;
 * returns: int blockFlag = INBLOCK (when not found) or GETXYZ (when found, next
 * stage)
 */
int saveElmInfo(char *line, ATOM_BLOCK *atom)
{
    /* Declaration and initialization variables for pcre2 */
    pcre2_match_data *match_data = NULL; /* pointer for match_data */
    PCRE2_SIZE *ovector = NULL;          /* pointer for ovector */
    PCRE2_UCHAR8 *buffer; /* pointer to buffer to get substring */
    PCRE2_SIZE size;      /* pointer to store size of substring */
    char *RegexStr = "ACL \"([0-9]{1,}) ([a-zA-Z]{1,})\"";
    if (reMatch(RegexStr, (PCRE2_SPTR8)line, &match_data, &ovector) == 3)
    {
        /* Get elmId */
        pcre2_substring_get_bynumber(match_data, 1, &buffer, &size);
        atom->elmId = atoi((const char *)buffer);
        /* Get element */
        pcre2_substring_get_bynumber(match_data, 2, &buffer, &size);
        atom->elm = strdup((const char *)buffer);
        atom->bCdSite = checkCdSite(line);
        atom->bStem = checkStem(line);
        pcre2_match_data_free(match_data);
        pcre2_substring_free(buffer);
        return GETXYZ;
    }
    else
    {
        pcre2_match_data_free(match_data);
        return INBLOCK;
    }
}
/** Save coordinates of the atom when found. Return GETXYZ or OUT to feedback
 * the control case.
 *  Args:
 *  	char *line;
 *  	ATOM_BLOCK *atom;
 *  returns: int blockFlag = GETXYZ (when not found) or NEXT (when found, next
 * block)
 */
int saveCoord(char *line, ATOM_BLOCK *atom)
{
    /* Declaration and initialization variables for pcre2 */
    pcre2_match_data *match_data = NULL; /* pointer for match_data */
    PCRE2_SIZE *ovector = NULL;          /* pointer for ovector */
    char *RegexStr = "XYZ \\(([0-9.e-]{1,}) ([0-9.e-]{1,}) ([0-9.e-]{1,})\\)";
    int rc;
    rc = reMatch(RegexStr, (PCRE2_SPTR8)line, &match_data, &ovector);
    if (rc == 4)
    {
        /* Stored x, y, z*/
        for (int j = 0; j < 3; ++j)
        {
            PCRE2_UCHAR8 *buffer; /* pointer to buffer to get substring */
            PCRE2_SIZE size;      /* pointer to store size of substring */
            pcre2_substring_get_bynumber(match_data, j + 1, &buffer, &size);
            atom->coord[j] = atof((const char *)buffer);
            pcre2_substring_free(buffer);
        }
        pcre2_match_data_free(match_data);
        return NEXT;
    }
    else
    {
        pcre2_match_data_free(match_data);
        return GETXYZ;
    }
}

/** Check if the atom is commented as the coord site. Return 1 or 0.
 */
int checkCdSite(char *line)
{
    /* Declaration and initialization variables for pcre2 */
    pcre2_match_data *match_data = NULL; /* pointer for match_data */
    PCRE2_SIZE *ovector = NULL;          /* pointer for ovector */
    char cdStr[] = "# cd";
    if ((reMatch(cdStr, (PCRE2_SPTR)line, &match_data, &ovector)) == 1)
    {
        pcre2_match_data_free(match_data);
        return 1;
    }
    else
    {
        pcre2_match_data_free(match_data);
        return 0;
    }
}
int checkStem(char *line)
{
    /* Declaration and initialization variables for pcre2 */
    pcre2_match_data *match_data = NULL; /* pointer for match_data */
    PCRE2_SIZE *ovector = NULL;          /* pointer for ovector */
    char cdStr[] = "# stem";
    if ((reMatch(cdStr, (PCRE2_SPTR)line, &match_data, &ovector)) == 1)
    {
        pcre2_match_data_free(match_data);
        return 1;
    }
    else
    {
        pcre2_match_data_free(match_data);
        return 0;
    }
}
