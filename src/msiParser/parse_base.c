#include "../maths/MyMaths.h"
#include "msiParser.h"

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
    assign_carbonVector(lattice, 41, lattice->carbon_chain_vec);
    assign_carbonVector(lattice, 73, lattice->carbon_metal_vec);
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
    pcre2_match_data *match_data = NULL; /* pointer for match_data */
    PCRE2_SIZE *ovector = NULL;          /* pointer for ovector */
    PCRE2_UCHAR8 *buffer; /* pointer to buffer to get substring */
    PCRE2_SIZE size;      /* pointer to store size of substring */

    char VecStr[] = "([A-C])3 \\(([0-9.-]{1,}) ([0-9.-]{1,}) ([0-9.-]{1,})\\)";
    char *VecName = NULL;
    char fail = 'N';

    if ((reMatch(VecStr, (PCRE2_SPTR)line, &match_data, &ovector)) > 1)
    {
        /* Get vector name */
        pcre2_substring_get_bynumber(match_data, 1, &buffer, &size);
        VecName = strdup((const char *)buffer);
        /* convert groups of coords to double vector[3] */
        for (int i = 0; i < 3; i++)
        {
            pcre2_substring_get_bynumber(match_data, i + 2, &buffer, &size);
            vector[i] = atof((const char *)buffer);
        }
        return *VecName;
    }
    else
        return fail;
}

/** Match "A I Id [0-9]{1.}" to get id of atom
 *  Args:
 *  	char *line;
 *  returns:
 *  	int rc; -1 or 1
 * */
static int matchAtomId(char *line)
{
    pcre2_match_data *match_data = NULL; /* pointer for match_data */
    PCRE2_SIZE *ovector = NULL;          /* pointer for ovector */
    char RegexStr[] = "A I Id [0-9]{1,}";
    int rc = 0;
    rc = reMatch(RegexStr, (PCRE2_SPTR)line, &match_data, &ovector);
    return rc; /* -1 or 1 */
}

void assign_carbonVector(BASE_LATTICE *lat, int bAtomId, double *vec)
{
    double *a, *b;
    a = lat->totalAtoms[40].coord;
    b = lat->totalAtoms[bAtomId].coord;
    initVector(a, b, vec);
}
