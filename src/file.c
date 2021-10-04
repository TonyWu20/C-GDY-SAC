#include "main.h"

long int countLines(FILE *file)
{
    long int lines = 0;
    while (EOF != (fscanf(file, "%*[^\n]"), fscanf(file, "%*c")))
        ++lines;
    rewind(file);
    return lines;
}

int returnLines(char *lineptr[], FILE *file)
{
    char *line = NULL; /* line array */
    size_t len = 0;
    ssize_t read;
    int counter = 0;
    /* char content[MAX_NUM_LINES][MAX_LINE_LENGTH]; */

    while ((read = getline(&line, &len, file)) !=
           -1) /* loop thru lines and add them to strArray */
    {
        /* fgets(line, 100, file) */
        lineptr[counter++] = strdup(line);
    }
    return counter;
}
