#include "../maths/MyMaths.h"
#include "msiParser.h"
#include "pcre2.h"

static int matchAtomId(char *line);
static char
scanLatVector(char *line,
              double *vector); /* Pass latVector[i] is a dereference */
int countAtoms(FILE *file)
{
    int atomNum = 0;
    char line[128];
    while (fgets(line, sizeof(line), file))
    {
        if (matchAtomId(line) == 1)
        {
            ++atomNum;
        }
    }
    rewind(file);
    return atomNum;
}
BASE_LATTICE *init_lattice(int atomNum)
{
    BASE_LATTICE *lattice;
    lattice = malloc(sizeof(BASE_LATTICE) + atomNum * sizeof(ATOM_BLOCK));
    lattice->atomNum = atomNum;
    return lattice;
}
BASE_LATTICE *parseBase(FILE *file)
{
    int atomNum = 0;
    atomNum = countAtoms(file);
    BASE_LATTICE *lattice;
    lattice = init_lattice(atomNum);
    get_LatVector(file, lattice->latVector);
    scanAtom(file, lattice->totalAtoms);
    /*assign_carbonVector(lattice, 41, lattice->carbon_chain_vec);*/
    /*assign_carbonVector(lattice, 73, lattice->carbon_metal_vec);*/
    return lattice;
}
int get_LatVector(FILE *file, double (*v)[3])
{
    char line[MAXLINE];
    int DONE = 0;
    while (fgets(line, sizeof(line), file))
    {
        switch (scanLatVector(line, *v))
        {
        case 'A':
            v++;
            break;
        case 'B':
            v++;
            break;
        case 'C':
            DONE = 1;
            break;
        default:
            break;
        }
    }
    rewind(file);
    return DONE;
}
/** Scan if the passed line contains lattice vectors, and modify the passed
 * double pointer *vector, which points to an array of three double values.
 * Args:
 * 		char *line;
 * 		double *vector;
 * returns:
 * 		char 'A', 'B', 'C', or 'N' if not found;
 */
static char scanLatVector(char *line, double *vector)
{
    /* Declaration and initialization variables for pcre2 */
    int errornumber;
    PCRE2_SIZE erroroffset;
    int rc;
    pcre2_code *re;
    PCRE2_SIZE size; /* pointer to store size of substring */

    PCRE2_SPTR VecStr =
        (PCRE2_SPTR) "([A-C])3 \\(([0-9.-]{1,}) ([0-9.-]{1,}) ([0-9.-]{1,})\\)";
    char fail = 'N';
    re = pcre2_compile(VecStr, PCRE2_ZERO_TERMINATED, 0, &errornumber,
                       &erroroffset, NULL);
    pcre2_match_data *match_data =
        pcre2_match_data_create_from_pattern(re, NULL);
    rc = pcre2_match(re, (PCRE2_SPTR)line,
                     (PCRE2_SIZE)strlen((const char *)line), 0, 0, match_data,
                     NULL);

    if (rc > 1)
    {
        /* Get vector name */
        char VecName[1];
        pcre2_substring_copy_bynumber(match_data, 1, (PCRE2_UCHAR8 *)VecName,
                                      &size);
        /* convert groups of coords to double vector[3] */
        for (int i = 0; i < 3; i++)
        {
            PCRE2_UCHAR8 *buffer;
            pcre2_substring_get_bynumber(match_data, i + 2, &buffer, &size);
            vector[i] = atof((const char *)buffer);
            pcre2_substring_free(buffer);
        }
        pcre2_match_data_free(match_data);
        pcre2_code_free(re);
        return VecName[0];
    }
    else
    {
        pcre2_match_data_free(match_data);
        pcre2_code_free(re);
        return fail;
    }
}

/** Match "A I Id [0-9]{1.}" to get id of atom
 *  Args:
 *  	char *line;
 *  returns:
 *  	int rc; -1 or 1
 * */
static int matchAtomId(char *line)
{
    char RegexStr[] = "A I Id [0-9]{1,}";
    int rc = 0;
    pcre2_code *re = init_re(RegexStr);
    pcre2_match_data *match_data =
        pcre2_match_data_create_from_pattern(re, NULL);
    rc = pcre2_match(re, (PCRE2_SPTR)line,
                     (PCRE2_SIZE)strlen((const char *)line), 0, 0, match_data,
                     NULL);
    pcre2_match_data_free(match_data);
    pcre2_code_free(re);
    return rc; /* -1 or 1 */
}

void assign_carbonVector(BASE_LATTICE *lat, int bAtomId, double *vec)
{
    double *a, *b;
    a = lat->totalAtoms[40].coord;
    b = lat->totalAtoms[bAtomId].coord;
    initVector(a, b, vec);
}
