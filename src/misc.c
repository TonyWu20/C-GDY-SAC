#include "misc.h"
#include "tasks.h"

void printProgress(int cur, int total, double percentage, char *name)
{
    int val = (int)(percentage * 100);
    int lpad = (int)(percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    int strPad = 32 - strlen(name);
    printf("\e[32m \rNow %.*s%*s %d/%d %3d%% [%.*s%*s]\e[m", 32, name, strPad,
           "", cur, total, val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}

int mkdir_p(char *dest)
{
    char *ret;
    int status;
    for (ret = strchr(dest, '/'); ret; ret= strchr(ret, '/'))
    {
        int pos = ret - dest;
        dest[pos] = '\0';
        struct stat s;
        int err = stat(dest, &s);
        if (err != -1)
        {
            dest[pos] = '/';
            ret++;
        }
        else {
            status = mkdir(dest, 0777);
            if (status)
            {
                printf("Fail to create %s\n", dest);
            }
            dest[pos] = '/';
            ret++;
        }
    }
    return status;
}

void createDirectory(char *dest)
{
    mkdir_p(dest);
}

char *extractStemName(char *filepath)
{
    char *buffer = strdup(filepath);
    char *last = strrchr(buffer, '/');
    char *stempath = NULL;
    if (last != NULL)
        stempath = last + 1;
    char *token, *rest;
    token = strtok_r(stempath, ".", &rest);
    char *ret = strdup(token);
    free(buffer);
    return ret;
}
