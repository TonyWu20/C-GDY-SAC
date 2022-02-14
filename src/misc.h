#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#define PBSTR "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"
#define PBWIDTH 40

/* Simple progress bar */
void printProgress(int cur, int total, double percentage, char *name);

/* recursive mkdir */
int mkdir_p(char *dest);
/* Create directory if not existed */
void createDirectory(char *dest);

/* Return filename without suffix */
char *extractStemName(char *filepath);

char *readWholeFile(const char *filePath);
