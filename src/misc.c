#include "misc.h"
#include <errno.h>

void printProgress(int cur, int total, double percentage, char *name)
{
    int val = (int)(percentage * 100);
    int lpad = (int)(percentage * PBWIDTH);
    int rpad = PBWIDTH - lpad;
    int strPad = 32 - strlen(name);
    printf("\e[32m \r%.*s%*s %d/%d %3d%% [%.*s%*s]\e[m", 32, name, strPad, "",
           cur, total, val, lpad, PBSTR, rpad, "");
    fflush(stdout);
}

void createDirectory(char *dest)
{
    struct stat s;
    if (stat(dest, &s) != -1)
    {
        return;
    }
    char *ret;
    int status;
    errno = 0;
    for (ret = strchr(dest, '/'); ret; ret = strchr(ret, '/'))
    {
        int pos = ret - dest;
        dest[pos] = '\0';
        struct stat s;
        int err = stat(dest, &s);
        if (err != -1)
        {
            dest[pos] = '/';
        }
        else
        {
            status = mkdir(dest, 0777);
            if (status && errno != EEXIST)
            {
                printf("Fail to create %s, ERROR:%d\n", dest, status);
            }
            dest[pos] = '/';
        }
        ret++;
    }
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

char *readWholeFile(const char *filePath)
{
    FILE *f = fopen(filePath, "r");
    if (!f)
    {
        printf("Open '%s' failed!\n", filePath);
    }
    fseek(f, 0, SEEK_END);
    long fSize = ftell(f);
    rewind(f);
    char *ret = malloc(fSize + 1);
    fread(ret, fSize, 1, f);
    fclose(f);
    ret[fSize] = 0;
    return ret;
}
