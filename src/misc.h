#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#define TOTAL_ELEMENT_NUM 44

void createDirectory(char *dest);

void allocateTasks(char *pathName);
/* Return filename without suffix */
char *extractStemName(char *filepath);
/* Malloc file path of the base model */
char *findBaseByElementId(int i);

/* return list of ads file paths */
char **pathway_adsLists(char *pathName, int *adsListLen);
