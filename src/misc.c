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

void createDirectory(char *dest)
{
    struct stat s;
    int err = stat(dest, &s);
    if (err == -1)
    {
        int destStringLen = strlen(dest);
        int commandLen = destStringLen + 9;
        char mkdir_command[commandLen + 1];
        snprintf(mkdir_command, commandLen + 1, "mkdir -p %s", dest);
        system(mkdir_command);
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
