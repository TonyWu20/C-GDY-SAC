#include "castep_output.h"
#include "castep_database.h"
#include "my_maths.h"
#include <string.h>

/* Cell starts */
Cell *createCell(Lattice *lat, CastepInfo **table)
{
    Cell *new = malloc(sizeof(Cell));
    new->lattice = lat;
    new->atomSorted = false;
    new->destroy = destroyCell;
    new->vtable = &cellVTable;
    new->textTable = &cellTextTable;
    new->infoTab = NULL;
    cellLoadInfoTab(new, table);
    return new;
}

void cellLoadInfoTab(Cell *self, CastepInfo **table)
{
    int elmNums = 0;
    char **elmList = self->vtable->sortElmList(self, &elmNums);
    for (int i = 0; i < elmNums; ++i)
    {
        CastepInfo *get = find_item(table, elmList[i]);
        add_item(&self->infoTab, get->info);
        free(elmList[i]);
    }
    free(elmList);
}

void destroyCell(Cell *self)
{
    self->lattice->vtable->destroy(self->lattice); // Destroy lattice here
    delete_all(&self->infoTab);
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
                find_item(&self->infoTab, self->lattice->metal_symbol);
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
    int elmNums = 0;
    char **elmList = sortedElementList(self, &elmNums);
    int lineLens[elmNums];
    char **tmpLines = malloc(sizeof(char *) * elmNums);
    int totalLen = 0;
    for (int i = 0; i < elmNums; ++i)
    {
        const char format[] = "%8s%18.10f\n";
        CastepInfo *item = find_item(&self->infoTab, elmList[i]);
        lineLens[i] =
            1 + snprintf(NULL, 0, format, elmList[i], item->info->mass);
        tmpLines[i] = malloc(lineLens[i]);
        snprintf(tmpLines[i], lineLens[i], format, elmList[i],
                 item->info->mass);
        totalLen += lineLens[i];
    }
    char *ret = calloc(totalLen + 1, sizeof(char));
    for (int i = 0; i < elmNums; ++i)
    {
        strncat(ret, tmpLines[i], lineLens[i]);
        free(tmpLines[i]);
    }
    free(tmpLines);
    return ret;
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
