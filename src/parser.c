#include "parser.h"
#include "atom.h"
#include "pcre2.h"
#include <stdint.h>

pcre2_code *init_re(char *RegexStr)
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

char **atom_block(char *subject, int *returnSize)
{
    char RegexStr[] =
        "\\(([0-9]+) Atom\r\n.*ACL \"([0-9]+) ([a-zA-Z]+).*\r\n.*\r\n.*"
        "XYZ \\(([0-9.e-]+) ([0-9.e-]+) ([0-9.e-]+).*\r\n"
        ".*Id ([0-9]+)";
    char **ret = malloc(sizeof(char *));
    pcre2_code *re = init_re(RegexStr);
    if (!re)
    {
        printf("Failed initializing re pattern");
        return NULL;
    }
    pcre2_match_data *match_data =
        pcre2_match_data_create_from_pattern(re, NULL);
    int rc = pcre2_match(re, (PCRE2_SPTR)subject, strlen(subject), 0, 0,
                         match_data, NULL);
    if (rc < 0)
    {
        switch (rc)
        {
        case PCRE2_ERROR_NOMATCH:
            printf("No match.\n");
            break;
        default:
            printf("Matching error%d\n", rc);
            break;
        }
        pcre2_match_data_free(match_data);
        pcre2_code_free(re);
        return NULL;
    }
    PCRE2_SIZE *ovector = pcre2_get_ovector_pointer(match_data);
    PCRE2_SPTR substring_start = (PCRE2_SPTR)subject + ovector[0];
    PCRE2_SIZE substring_length = ovector[1] - ovector[0];
    (*returnSize)++;
    ret = realloc(ret, sizeof(char *) * (*returnSize));
    ret[*returnSize - 1] = strndup((char *)substring_start, substring_length);
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
            return ret;
        }

        if (rc == 0)
            printf("ovector was not big enough for all the captured "
                   "substrings\n");
        PCRE2_SPTR substring_start = (PCRE2_SPTR)subject + ovector[0];
        size_t substring_length = ovector[1] - ovector[0];
        (*returnSize)++;
        ret = realloc(ret, sizeof(char *) * (*returnSize));
        ret[*returnSize - 1] =
            strndup((char *)substring_start, (int)substring_length);
    }
    pcre2_match_data_free(match_data);
    pcre2_code_free(re);
    return ret;
}
