#include "misc.h"

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
