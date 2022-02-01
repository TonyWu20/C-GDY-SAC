#include "cell.h"
#include "castep_database.h"
#include "misc.h"
#include "my_maths.h"
#include <string.h>

struct Cell_vtable cellVTable = {sortAtomsByElement, sortedElementList,
                                 cellExport};
struct Cell_textFunc cellTextTable = {cellWriteBlock};

/* Cell starts */
Cell *createCell(Lattice *lat, CastepInfo *table)
{
    Cell *new = malloc(sizeof(Cell));
    new->lattice = lat;
    new->atomSorted = false;
    new->destroy = destroyCell;
    new->vtable = &cellVTable;
    new->textTable = &cellTextTable;
    int elmNums = 0;
    new->elmLists = new->vtable->sortElmList(new, &elmNums);
    new->elmNums = elmNums;
    new->infoTab = NULL;
    cellLoadInfoTab(new, table);
    if (!new->infoTab)
    {
        printf("infoTab loading failed\n");
    }
    return new;
}

void cellLoadInfoTab(Cell *self, CastepInfo *table)
{
    for (int i = 0; i < self->elmNums; ++i)
    {
        CastepInfo *get = find_item(table, self->elmLists[i]);
        add_item(&self->infoTab, get->info);
    }
}

void destroyCell(Cell *self)
{
    self->lattice->vtable->destroy(self->lattice); // Destroy lattice here
    delete_all(&self->infoTab);
    for (int i = 0; i < self->elmNums; ++i)
        free(self->elmLists[i]);
    free(self->elmLists);
    free(self);
}

/* Returns a malloced string
 */
char *cellWriteBlock(Cell *self, char *blockName,
                     char *(*blockTextWriter)(Cell *self))
{
    char *content = blockTextWriter(self);
    int needed = snprintf(NULL, 0, "%%BLOCK %s\n%s%%ENDBLOCK %s\n\n", blockName,
                          content, blockName);
    char *output = malloc(needed + 1);
    snprintf(output, needed + 1, "%%BLOCK %s\n%s%%ENDBLOCK %s\n\n", blockName,
             content, blockName);
    free(content);
    return output;
}

char *cell_latticeVector_writer(Cell *self)
{
    double a[3], b[3], c[3];
    Lattice *lat = self->lattice;
    Matrix *latVectors = lat->lattice_vectors;
    for (int i = 0; i < 3; ++i)
    {
        a[i] = latVectors->value[i][0];
        b[i] = latVectors->value[i][1];
        c[i] = latVectors->value[i][2];
    }
    int bufSize =
        snprintf(NULL, 0,
                 "%24.18f%24.18f%24.18f\n%24.18f%24.18f%24.18f\n%24."
                 "18f%24.18f%24.18f\n",
                 a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2]) +
        1;
    char *resString = malloc(bufSize);
    snprintf(
        resString, bufSize,
        "%24.18f%24.18f%24.18f\n%24.18f%24.18f%24.18f\n%24.18f%24.18f%24.18f\n",
        a[0], a[1], a[2], b[0], b[1], b[2], c[0], c[1], c[2]);
    return resString;
}

char *cell_fracCoord_writer(Cell *self)
{
    if (self->atomSorted == false)
        self->vtable->sortAtoms(self);
    Molecule *mol = self->lattice->_mol;
    Matrix *molCoord = mol->vtable->get_mol_coords(mol);
    Matrix *xyzToFrac = fractionalCoordMatrix(self->lattice->lattice_vectors);
    Matrix *fracCoord;
    multiply_matrices(xyzToFrac, molCoord, &fracCoord);
    destroy_matrix(xyzToFrac);
    destroy_matrix(molCoord);
    free(xyzToFrac);
    free(molCoord);
    int lineLens[mol->atomNum];
    char **lines = malloc(sizeof(char *) * mol->atomNum);
    int totalLen = 0;
    for (int i = 0; i < mol->atomNum; ++i)
    {
        Atom *cur = mol->atom_arr[i];
        Matrix *cd = fracCoord;
        if (i < mol->atomNum - 1)
        {
            lineLens[i] = 1 + snprintf(NULL, 0, "%3s%20.16f%20.16f%20.16f\n",
                                       cur->element, cd->value[0][i],
                                       cd->value[1][i], cd->value[2][i]);
            lines[i] = malloc(sizeof(char) * lineLens[i]);
            snprintf(lines[i], lineLens[i], "%3s%20.16f%20.16f%20.16f\n",
                     cur->element, cd->value[0][i], cd->value[1][i],
                     cd->value[2][i]);
        }
        else
        {
            CastepInfo *metal =
                find_item(self->infoTab, self->lattice->metal_symbol);
            if (metal->info->spin)
            {
                lineLens[i] =
                    1 + snprintf(NULL, 0,
                                 "%3s%20.16f%20.16f%20.16f SPIN=%14.10f\n",
                                 cur->element, cd->value[0][i], cd->value[1][i],
                                 cd->value[2][i], (double)metal->info->spin);
                lines[i] = malloc(sizeof(char) * lineLens[i]);
                snprintf(lines[i], lineLens[i],
                         "%3s%20.16f%20.16f%20.16f SPIN=%14.10f\n",
                         cur->element, cd->value[0][i], cd->value[1][i],
                         cd->value[2][i], (double)metal->info->spin);
            }
            else
            {
                lineLens[i] =
                    1 + snprintf(NULL, 0, "%3s%20.16f%20.16f%20.16f\n",
                                 cur->element, cd->value[0][i], cd->value[1][i],
                                 cd->value[2][i]);
                lines[i] = malloc(sizeof(char) * lineLens[i]);
                snprintf(lines[i], lineLens[i], "%3s%20.16f%20.16f%20.16f\n",
                         cur->element, cd->value[0][i], cd->value[1][i],
                         cd->value[2][i]);
            }
        }
        totalLen += lineLens[i];
    }
    destroy_matrix(fracCoord);
    free(fracCoord);
    char *blockText = calloc(totalLen + 1, sizeof(char));
    for (int i = 0; i < mol->atomNum; ++i)
    {
        strncat(blockText, lines[i], lineLens[i]);
        free(lines[i]);
    }
    free(lines);
    return blockText;
}

char *cell_kPointsList_writer(Cell *self)
{
    char line[] = "   0.0000000000000000   0.0000000000000000   "
                  "0.0000000000000000       1.000000000000000\n";
    char *ret = strdup(line);
    return ret;
}

char *cell_miscOptions_writer(Cell *self)
{
    char line[] = "FIX_ALL_CELL : true\nFIX_COM : false\n";
    char *ret = strdup(line);
    return ret;
}

char *cell_ionicConstraints_writer(Cell *self)
{
    char *ret = strdup("");
    return ret;
}
char *cell_externalEfield_writer(Cell *self)
{
    char text[] = "    0.0000000000     0.0000000000     0.0000000000\n";
    char *ret = strdup(text);
    return ret;
}

char *cell_externalPressure_writer(Cell *self)
{
    char line[] = "    0.0000000000    0.0000000000    0.0000000000\n          "
                  "          0.0000000000    0.0000000000\n                    "
                  "                0.0000000000\n";
    char *ret = strdup(line);
    return ret;
}

char *cell_speciesMass_writer(Cell *self)
{
    int lineLens[self->elmNums];
    char **tmpLines = malloc(sizeof(char *) * self->elmNums);
    int totalLen = 0;
    for (int i = 0; i < self->elmNums; ++i)
    {
        const char format[] = "%8s%18.10f\n";
        CastepInfo *item = find_item(self->infoTab, self->elmLists[i]);
        lineLens[i] =
            1 + snprintf(NULL, 0, format, self->elmLists[i], item->info->mass);
        tmpLines[i] = malloc(lineLens[i]);
        snprintf(tmpLines[i], lineLens[i], format, self->elmLists[i],
                 item->info->mass);
        totalLen += lineLens[i];
    }
    char *ret = calloc(totalLen + 1, sizeof(char));
    for (int i = 0; i < self->elmNums; ++i)
    {
        strncat(ret, tmpLines[i], lineLens[i]);
        free(tmpLines[i]);
    }
    free(tmpLines);
    return ret;
}

char *cell_speciesPot_writer(Cell *self)
{
    int lineLens[self->elmNums];
    char **tmpLines = malloc(sizeof(char *) * self->elmNums);
    int totalLen = 0;
    const char format[] = "%8s  %s\n";
    for (int i = 0; i < self->elmNums; ++i)
    {
        CastepInfo *item = find_item(self->infoTab, self->elmLists[i]);
        char *potential_stem = strrchr(item->info->potential_file, '/') + 1;
        lineLens[i] =
            1 + snprintf(NULL, 0, format, self->elmLists[i], potential_stem);
        tmpLines[i] = malloc(lineLens[i]);
        snprintf(tmpLines[i], lineLens[i], format, self->elmLists[i],
                 potential_stem);
        totalLen += lineLens[i];
    }
    char *ret = calloc(totalLen + 1, sizeof(char));
    for (int i = 0; i < self->elmNums; ++i)
    {
        strncat(ret, tmpLines[i], lineLens[i]);
        free(tmpLines[i]);
    }
    free(tmpLines);
    return ret;
}

char *cell_speciesLCAOstates_writer(Cell *self)
{
    int lineLens[self->elmNums];
    char **tmpLines = malloc(sizeof(char *) * self->elmNums);
    int totalLen = 0;
    const char format[] = "%8s%10d\n";
    for (int i = 0; i < self->elmNums; ++i)
    {
        CastepInfo *item = find_item(self->infoTab, self->elmLists[i]);
        lineLens[i] =
            1 + snprintf(NULL, 0, format, self->elmLists[i], item->info->LCAO);
        tmpLines[i] = malloc(lineLens[i]);
        snprintf(tmpLines[i], lineLens[i], format, self->elmLists[i],
                 item->info->LCAO);
        totalLen += lineLens[i];
    }
    char *ret = calloc(totalLen + 1, sizeof(char));
    for (int i = 0; i < self->elmNums; ++i)
    {
        strncat(ret, tmpLines[i], lineLens[i]);
        free(tmpLines[i]);
    }
    free(tmpLines);
    return ret;
}

void cellExport(Cell *self, bool DOS)
{
    /* Filepath processing */
    char *stemName = self->lattice->_mol->name;
    char *exportDir = self->lattice->vtable->exportDir(self->lattice,
                                                       self->lattice->pathName);
    /* dest dir with "_opt" suffix */
    int subDirLen = 1 + snprintf(NULL, 0, "%s%s_opt/", exportDir, stemName);
    char *subDir = malloc(subDirLen);
    snprintf(subDir, subDirLen, "%s%s_opt/", exportDir, stemName);
    free(exportDir);
    /* Create directory routine in misc.h */
    createDirectory(subDir);
    char *fileName;
    /* Differ for DOS or not */
    if (DOS == false)
    {
        int fileNameLen = 1 + snprintf(NULL, 0, "%s%s.cell", subDir, stemName);
        fileName = malloc(fileNameLen);
        snprintf(fileName, fileNameLen, "%s%s.cell", subDir, stemName);
    }
    else
    {
        int fileNameLen =
            1 + snprintf(NULL, 0, "%s%s_DOS.cell", subDir, stemName);
        fileName = malloc(fileNameLen);
        snprintf(fileName, fileNameLen, "%s%s_DOS.cell", subDir, stemName);
    }
    /* Release malloc'd memory subDir */
    free(subDir);
    /* Ready to write */
    FILE *writeTo = fopen(fileName, "w");
    free(fileName);
    /* Get,Write and Free strings for each section */
    char *latVec = self->textTable->blockWriter(self, "LATTICE_CART",
                                                cell_latticeVector_writer);
    fputs(latVec, writeTo);
    free(latVec);
    char *fracCoord = self->textTable->blockWriter(self, "POSITIONS_FRAC",
                                                   cell_fracCoord_writer);
    fputs(fracCoord, writeTo);
    free(fracCoord);
    if (DOS)
    {
        char *BS_kPoints = self->textTable->blockWriter(
            self, "BS_KPOINTS_LIST", cell_kPointsList_writer);
        fputs(BS_kPoints, writeTo);
        free(BS_kPoints);
    }
    char *kPoints = self->textTable->blockWriter(self, "KPOINTS_LIST",
                                                 cell_kPointsList_writer);
    fputs(kPoints, writeTo);
    free(kPoints);
    char *misc = cell_miscOptions_writer(self);
    fputs(misc, writeTo);
    free(misc);
    char *ionicConstraint = self->textTable->blockWriter(
        self, "IONIC_CONSTRAINTS", cell_ionicConstraints_writer);
    fputs(ionicConstraint, writeTo);
    free(ionicConstraint);
    char *efield = self->textTable->blockWriter(self, "EXTERNAL_EFIELD",
                                                cell_externalEfield_writer);
    fputs(efield, writeTo);
    free(efield);
    char *ePressure = self->textTable->blockWriter(
        self, "EXTERNAL_PRESSURE", cell_externalPressure_writer);
    fputs(ePressure, writeTo);
    free(ePressure);
    char *spMass = self->textTable->blockWriter(self, "SPECIES_MASS",
                                                cell_speciesMass_writer);
    fputs(spMass, writeTo);
    free(spMass);
    char *spPot = self->textTable->blockWriter(self, "SPECIES_POT",
                                               cell_speciesPot_writer);
    fputs(spPot, writeTo);
    free(spPot);
    char *spLCAO = self->textTable->blockWriter(self, "SPECIES_LCAO_STATES",
                                                cell_speciesLCAOstates_writer);
    fputs(spLCAO, writeTo);
    free(spLCAO);
    /* Close file release pointer */
    fclose(writeTo);
}

static int atomCmp(const void *a, const void *b)
{
    Atom *atomA = *(Atom **)a;
    Atom *atomB = *(Atom **)b;
    int aNum = atomA->elementId;
    int bNum = atomB->elementId;
    if (aNum == bNum)
        return atomA->atomId - atomB->atomId;
    return aNum - bNum;
}

void sortAtomsByElement(Cell *self)
{
    Molecule *mol = self->lattice->_mol;
    Atom **atomArray = mol->atom_arr;
    qsort(atomArray, mol->atomNum, sizeof(Atom *), atomCmp);
    mol->atom_arr = atomArray;
    self->atomSorted = true;
}

char **sortedElementList(Cell *self, int *returnSize)
{
    if (self->atomSorted == false)
        self->vtable->sortAtoms(self);
    Molecule *mol = self->lattice->_mol;
    int size = 1;
    char **elementsSet = malloc(sizeof(char *));
    elementsSet[0] = strdup(mol->atom_arr[0]->element);
    for (int i = 1; i < mol->atomNum; ++i)
    {
        Atom *cur = mol->atom_arr[i];
        Atom *prev = mol->atom_arr[i - 1];
        if (prev->elementId == cur->elementId)
            continue;
        else
        {
            ++size;
            elementsSet = realloc(elementsSet, sizeof(char *) * (size));
            elementsSet[size - 1] = strdup(cur->element);
        };
    }
    *returnSize = size;
    return elementsSet;
}

/* Cell ends */
