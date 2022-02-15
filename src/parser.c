#include "parser.h"
#include "misc.h"
#include "pcre2.h"
#include <stdio.h>
#include <string.h>

static pcre2_code *init_re(char *RegexStr);
static void re_match(pcre2_code *re_pattern, pcre2_match_data **match_data,
                     int *rc, char *subject);
static Atom *parse_atom(char *atom_block);
static Atom **get_all_atoms(char *subject, int *returnSize);

// [>regex<]
static pcre2_code *init_re(char *RegexStr)
{
    int errornumber;
    PCRE2_SIZE erroroffset;
    pcre2_code *re;
    re = pcre2_compile((PCRE2_SPTR)RegexStr, PCRE2_ZERO_TERMINATED,
                       PCRE2_MULTILINE, &errornumber, &erroroffset, NULL);
    if (re == NULL)
    {
        PCRE2_UCHAR buffer[256];
        pcre2_get_error_message(errornumber, buffer, sizeof(buffer));
        printf("PCRE2 compilation failed at offset %d: %s\n", (int)erroroffset,
               buffer);
    }
    return re;
}

static void re_match(pcre2_code *re_pattern, pcre2_match_data **match_data,
                     int *rc, char *subject)
{
    if (!re_pattern)
    {
        printf("Failed initializing re_pattern");
        return;
    }
    *match_data = pcre2_match_data_create_from_pattern(re_pattern, NULL);
    *rc = pcre2_match(re_pattern, (PCRE2_SPTR)subject, strlen(subject), 0, 0,
                      *match_data, NULL);
    if (*rc < 0)
    {
        switch (*rc)
        {
        case PCRE2_ERROR_NOMATCH:
            printf("No match.\n");
            break;
        default:
            printf("Matching error%d\n", *rc);
            break;
        }
        pcre2_match_data_free(*match_data);
        pcre2_code_free(re_pattern);
        return;
    }
}

/* Parse the lattice vectors */
static void get_lattice_vectors(char *subject, simd_double3x3 *result)
{
    char *RegexStr = ".*A3 \\(([0-9e. -]+)\\)\\)\r\n"
                     ".*B3 \\(([0-9e. -]+)\\)\\)\r\n"
                     ".*C3 \\(([0-9e. -]+)\\)\\)\r\n";
    pcre2_code *re = init_re(RegexStr);
    pcre2_match_data *match_data;
    int rc = 0;
    re_match(re, &match_data, &rc, subject);
    PCRE2_UCHAR8 *buffer = NULL;
    PCRE2_SIZE size = 0;
    simd_double3 xyz[3];
    for (int i = 0; i < 3; ++i)
    {
        pcre2_substring_get_bynumber(match_data, i + 1, &buffer, &size);
        char *rest = NULL;
        char *token;
        int j = 0;
        for (token = strtok_r((char *)buffer, " ", &rest); token && j < 3;
             token = strtok_r(NULL, " ", &rest), ++j)
        {
            xyz[i][j] = atof((const char *)token);
        }
        (*result).columns[i] = xyz[i];
        pcre2_substring_free(buffer);
    }
    pcre2_code_free(re);
    pcre2_match_data_free(match_data);
}

/* Match and parse all atoms
 * Returns: array of Atom *
 */
/* static Atom **get_all_atoms(char *subject, int *returnSize) */
static Atom **get_all_atoms(char *subject, int *returnSize)
{
    char RegexStr[] = "\\(([0-9]+) Atom\\R.*ACL \"([0-9]+) "
                      "([a-zA-Z]+).*\\R.*\\R.*"
                      "XYZ \\(([0-9.e-]+) ([0-9.e-]+) ([0-9.e-]+).*\\R"
                      ".*Id ([0-9]+)";
    Atom **ret = malloc(sizeof(Atom *));
    pcre2_code *re = init_re(RegexStr);
    pcre2_match_data *match_data;
    int rc = 0;
    re_match(re, &match_data, &rc, subject);
    PCRE2_SIZE *ovector = pcre2_get_ovector_pointer(match_data);
    PCRE2_SPTR substring_start = (PCRE2_SPTR)subject + ovector[0];
    PCRE2_SIZE substring_length = ovector[1] - ovector[0];
    char *substring = malloc(substring_length + 1);
    substring = memcpy(substring, substring_start, substring_length);
    substring[substring_length] = 0;
    (*returnSize)++;
    ret = realloc(ret, sizeof(Atom *) * (*returnSize));
    ret[*returnSize - 1] = parse_atom(substring);

    // Loop for second and subsequent matches
    for (;;)
    {
        uint32_t option = 0;
        PCRE2_SIZE start_offset = ovector[1];
        if (ovector[0] == ovector[1])
        {
            if (ovector[0] == strlen(subject))
                break;
        }
        //[ > Run the next matching operation < ]
        rc = pcre2_match(re, (PCRE2_SPTR)subject, strlen(subject), start_offset,
                         option, match_data, NULL);
        if (rc < 0)
        {
            pcre2_match_data_free(match_data);
            pcre2_code_free(re);
            free(substring);
            return ret;
        }

        if (rc == 0)
            printf("ovector was not big enough for all the captured "
                   "substrings\n");
        substring_start = (PCRE2_SPTR)subject + ovector[0];
        substring_length = ovector[1] - ovector[0];
        (*returnSize)++;
        substring = realloc(substring, substring_length + 1);
        memcpy(substring, substring_start, substring_length);
        substring[substring_length] = 0;
        ret = realloc(ret, sizeof(Atom *) * (*returnSize));
        ret[*returnSize - 1] = parse_atom(substring);
    }
    pcre2_match_data_free(match_data);
    pcre2_code_free(re);
    free(substring);
    return ret;
}

/* Parse single atom block */
static Atom *parse_atom(char *atom_block)
{
    char RegexStr[] = "\\([0-9]+ Atom\\R.*ACL \"([0-9]+) ([a-zA-Z"
                      "]+).*\\R.*Label \".*\r\n.*"
                      "XYZ \\(([0-9.e-]+) ([0-9.e-]+) ([0-9.e-]+).*\\R"
                      ".*Id ([0-9]+)";
    pcre2_code *re = init_re(RegexStr);
    pcre2_match_data *match_data;
    int rc = 0;
    re_match(re, &match_data, &rc, atom_block);
    PCRE2_UCHAR8 *buffer = NULL;
    PCRE2_SIZE size = 0;
    pcre2_substring_get_bynumber(match_data, 1, &buffer, &size);
    int elementId = atoi((const char *)buffer);
    pcre2_substring_free(buffer);
    // Element
    pcre2_substring_get_bynumber(match_data, 2, &buffer, &size);
    char *element = strdup((const char *)buffer);
    pcre2_substring_free(buffer);
    // XYZ matrix
    simd_double3 coord;
    for (int i = 0; i < 3; ++i)
    {
        pcre2_substring_get_bynumber(match_data, i + 3, &buffer, &size);
        coord[i] = atof((const char *)buffer);
        pcre2_substring_free(buffer);
    }
    // Atom Id
    pcre2_substring_get_bynumber(match_data, 6, &buffer, &size);
    int atomId = atoi((const char *)buffer);
    pcre2_substring_free(buffer);
    // Create new Atom object
    Atom *new = createAtom(element, coord, atomId, elementId);
    // Parse the element atomic number from ACL Label
    // Free memory
    free(element);
    pcre2_match_data_free(match_data);
    pcre2_code_free(re);
    return new;
}

Adsorbate *parse_adsorbate_from_file(char *fileName, char *name,
                                     AdsorbateInfo *adsInfo)
{
    char *body = readWholeFile(fileName);
    int atomNums = 0;
    Atom **atom_arr = get_all_atoms(body, &atomNums);
    free(body);
    Molecule *mol = createMolecule(name, atomNums, atom_arr);
    Adsorbate *ads =
        createAdsorbate(mol, adsInfo->coordAtomNum, adsInfo->coordAtomIds,
                        adsInfo->stemAtomIds, adsInfo->planeAtomIds,
                        adsInfo->vertical, adsInfo->bSym, adsInfo->upperAtomId);
    return ads;
}

Lattice *parse_lattice_from_file(char *fileName, char *name)
{
    char *body = readWholeFile((const char *)fileName);
    int atom_nums = 0;
    Atom **atom_arr = get_all_atoms(body, &atom_nums);
    Molecule *lat_mol = createMolecule(name, atom_nums, atom_arr);
    simd_double3x3 latVectors;
    get_lattice_vectors(body, &latVectors);
    Lattice *new = createLattice(lat_mol, latVectors);
    free(body);
    return new;
}
