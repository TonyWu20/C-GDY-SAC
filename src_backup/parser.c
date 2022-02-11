#include "parser.h"
#include "atom.h"
#include "lattice.h"
#include "matrix.h"
#include "misc.h"
#include "molecule.h"
#include "pcre2.h"
#include <stdio.h>
#include <string.h>

static pcre2_code *init_re(char *RegexStr);
static void re_match(pcre2_code *re_pattern, pcre2_match_data **match_data,
                     int *rc, char *subject);
static void get_cd_num(char *subject, int *cd_num);
static void get_cd_info(char *subject, int *cd_num, int **cd_ids);
static void get_symmetric_info(char *subject, int *bSym);
static void get_stem_arr(char *subject, int *stem_arr);
static void get_plane_arr(char *subject, int *plane_arr);
static void get_upperAtomId(char *subject, int *upperAtomId);
static Atom *parse_atom(char *atom_block);
static Atom **get_all_atoms(char *subject, int *returnSize);
// Table

// UThash

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

/* Parse number of coordinate atoms */
static void get_cd_num(char *subject, int *cd_num)
{
    char RegexStr[] = "# cd_num: ([0-9])";
    pcre2_code *re = init_re(RegexStr);
    pcre2_match_data *match_data;
    int rc = 0;
    re_match(re, &match_data, &rc, subject);
    PCRE2_UCHAR8 *buffer = NULL;
    PCRE2_SIZE size = 0;
    pcre2_substring_get_bynumber(match_data, 1, &buffer, &size);
    *cd_num = atoi((const char *)buffer);
    pcre2_substring_free(buffer);
    pcre2_match_data_free(match_data);
    pcre2_code_free(re);
}

/* Parse atom Ids of coordinate atoms and store in **cd_ids
 */
static void get_cd_info(char *subject, int *cd_num, int **cd_ids)
{
    get_cd_num(subject, cd_num);
    char RegexStr[] = "# cd_ids: ([0-9,]+)";
    pcre2_code *re = init_re(RegexStr);
    pcre2_match_data *match_data;
    int rc = 0;
    re_match(re, &match_data, &rc, subject);
    PCRE2_UCHAR8 *buffer = NULL;
    PCRE2_SIZE size = 0;
    pcre2_substring_get_bynumber(match_data, 1, &buffer, &size);
    char *rest = NULL;
    char *token;
    *cd_ids = realloc(*cd_ids, sizeof(int) * (*cd_num));
    int count = 0;
    for (token = strtok_r((char *)buffer, ",", &rest); token;
         token = strtok_r(NULL, ",", &rest))
    {
        (*cd_ids)[count++] = atoi(token);
    }
    pcre2_code_free(re);
    pcre2_substring_free(buffer);
    pcre2_match_data_free(match_data);
}

/* Get if the molecule is symmetric */
static void get_symmetric_info(char *subject, int *bSym)
{
    char RegexStr[] = "# Symmetric: ([0-1])";
    pcre2_code *re = init_re(RegexStr);
    pcre2_match_data *match_data;
    int rc = 0;
    re_match(re, &match_data, &rc, subject);
    PCRE2_UCHAR *buffer = NULL;
    PCRE2_SIZE size = 0;
    pcre2_substring_get_bynumber(match_data, 1, &buffer, &size);
    *bSym = atoi((const char *)buffer);
    pcre2_substring_free(buffer);
    pcre2_match_data_free(match_data);
    pcre2_code_free(re);
}

/* Store the stem atomId in the pointer */
static void get_stem_arr(char *subject, int *stem_arr)
{
    char rg[] = "# stem_ids: ([0-9,]+)";
    pcre2_code *re = init_re(rg);
    pcre2_match_data *match_data;
    int rc = 0;
    re_match(re, &match_data, &rc, subject);
    PCRE2_UCHAR8 *buffer = NULL;
    PCRE2_SIZE size = 0;
    pcre2_substring_get_bynumber(match_data, 1, &buffer, &size);
    char *rest = NULL;
    char *token;
    int i = 0;
    for (token = strtok_r((char *)buffer, ",", &rest); token && i < 2;
         token = strtok_r(NULL, ",", &rest), ++i)
    {
        stem_arr[i] = atoi((const char *)token);
    }
    pcre2_code_free(re);
    pcre2_substring_free(buffer);
    pcre2_match_data_free(match_data);
}

/* Store the plane atomId in the pointer */
static void get_plane_arr(char *subject, int *plane_arr)
{
    char rg[] = "# plane_ids: ([0-9,]+)";
    pcre2_code *re = init_re(rg);
    pcre2_match_data *match_data;
    int rc = 0;
    re_match(re, &match_data, &rc, subject);
    PCRE2_UCHAR8 *buffer = NULL;
    PCRE2_SIZE size = 0;
    pcre2_substring_get_bynumber(match_data, 1, &buffer, &size);
    char *rest = NULL;
    char *token;
    int i = 0;
    for (token = strtok_r((char *)buffer, ",", &rest); token && i < 3;
         token = strtok_r(NULL, ",", &rest), ++i)
    {
        plane_arr[i] = atoi((const char *)token);
    }
    pcre2_code_free(re);
    pcre2_substring_free(buffer);
    pcre2_match_data_free(match_data);
}

/* Parse the lattice vectors */
static void get_lattice_vectors(char *subject, Matrix **result)
{
    *result = create_matrix(4, 3);
    if (*result == NULL)
    {
        printf("Can't create matrix.\n");
        return;
    }
    char *RegexStr = ".*A3 \\(([0-9e. -]+)\\)\\)\r\n"
                     ".*B3 \\(([0-9e. -]+)\\)\\)\r\n"
                     ".*C3 \\(([0-9e. -]+)\\)\\)\r\n";
    pcre2_code *re = init_re(RegexStr);
    pcre2_match_data *match_data;
    int rc = 0;
    re_match(re, &match_data, &rc, subject);
    PCRE2_UCHAR8 *buffer = NULL;
    PCRE2_SIZE size = 0;
    for (int j = 0; j < 3; ++j)
    {
        pcre2_substring_get_bynumber(match_data, j + 1, &buffer, &size);
        char *rest = NULL;
        char *token;
        int i = 0;
        for (token = strtok_r((char *)buffer, " ", &rest); token && i < 3;
             token = strtok_r(NULL, " ", &rest), ++i)
        {
            (*result)->value[i][j] = atof((const char *)token);
        }
        pcre2_substring_free(buffer);
    }
    for (int i = 0; i < 3; ++i)
    {
        (*result)->value[3][i] = 1;
    }
    pcre2_code_free(re);
    pcre2_match_data_free(match_data);
}

/* Parse the upper atom Id
 * store in passed int pointer
 */
static void get_upperAtomId(char *subject, int *upperAtomId)
{
    char RegexStr[] = "# Up_atom: ([0-9])";
    pcre2_code *re = init_re(RegexStr);
    pcre2_match_data *match_data;
    int rc = 0;
    re_match(re, &match_data, &rc, subject);
    PCRE2_UCHAR8 *buffer = NULL;
    PCRE2_SIZE size = 0;
    pcre2_substring_get_bynumber_8(match_data, 1, &buffer, &size);
    *upperAtomId = atoi((const char *)buffer);
    pcre2_substring_free_8(buffer);
    pcre2_match_data_free_8(match_data);
    pcre2_code_free_8(re);
}

/* Match and parse all atoms
 * Returns: array of Atom *
 */
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
    char RegexStr[] = "\\(([0-9]+) Atom\\R.*ACL \"([0-9a-zA-Z "
                      "]+).*\\R.*Label \"([a-zA-Z]+).*\r\n.*"
                      "XYZ \\(([0-9.e-]+) ([0-9.e-]+) ([0-9.e-]+).*\\R"
                      ".*Id ([0-9]+)";
    pcre2_code *re = init_re(RegexStr);
    pcre2_match_data *match_data;
    int rc = 0;
    re_match(re, &match_data, &rc, atom_block);
    PCRE2_UCHAR8 *buffer = NULL;
    PCRE2_SIZE size = 0;
    // Tree Id
    pcre2_substring_get_bynumber(match_data, 1, &buffer, &size);
    int treeId = atoi((const char *)buffer);
    pcre2_substring_free(buffer);
    // Label
    pcre2_substring_get_bynumber(match_data, 2, &buffer, &size);
    char *label = strdup((const char *)buffer);
    pcre2_substring_free(buffer);
    // Element
    pcre2_substring_get_bynumber(match_data, 3, &buffer, &size);
    char *element = strdup((const char *)buffer);
    pcre2_substring_free(buffer);
    // XYZ matrix
    Matrix *coord = create_matrix(4, 1);
    for (int i = 0; i < 3; ++i)
    {
        pcre2_substring_get_bynumber(match_data, i + 4, &buffer, &size);
        coord->value[i][0] = atof((const char *)buffer);
        pcre2_substring_free(buffer);
    }
    coord->value[3][0] = 1;
    // Atom Id
    pcre2_substring_get_bynumber(match_data, 7, &buffer, &size);
    int atomId = atoi((const char *)buffer);
    // Create new Atom object
    char *labelDup = strdup(label);
    char *token, *rest;
    token = strtok_r(labelDup, " ", &rest);
    Atom *new = createAtom(element, label, coord, atomId, treeId);
    new->elementId = atoi(token);
    // Parse the element atomic number from ACL Label
    // Free memory
    pcre2_substring_free(buffer);
    pcre2_match_data_free(match_data);
    pcre2_code_free(re);
    free(labelDup);
    return new;
}

Adsorbate *parse_molecule_from_file(char *fileName, char *name)
{
    char *body = readWholeFile(fileName);
    int atom_nums = 0;
    Atom **atom_arr = get_all_atoms(body, &atom_nums);
    Molecule *mol = NULL;
    int cd_num = 0;
    int *cd_arr = NULL;
    int stem_arr[2];
    int plane_arr[3];
    int bSym = 0;
    int upperAtomId = 0;
    get_cd_info(body, &cd_num, &cd_arr);
    get_stem_arr(body, stem_arr);
    get_plane_arr(body, plane_arr);
    get_symmetric_info(body, &bSym);
    get_upperAtomId(body, &upperAtomId);
    mol = createMolecule(name, atom_nums, atom_arr);
    Adsorbate *ads = createAdsorbate(mol, cd_num, cd_arr, stem_arr, plane_arr,
                                     bSym, upperAtomId);
    free(cd_arr);
    free(body);
    return ads;
}

Lattice *parse_lattice_from_file(char *fileName, char *name)
{
    char *body = readWholeFile((const char *)fileName);
    int atom_nums = 0;
    Atom **atom_arr = get_all_atoms(body, &atom_nums);
    Molecule *lat_mol = createMolecule(name, atom_nums, atom_arr);
    Matrix *vectors;
    get_lattice_vectors(body, &vectors);
    Lattice *new = createLattice(lat_mol, vectors);
    free(body);
    return new;
}