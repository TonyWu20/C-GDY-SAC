#include "cell.h"
#include "misc.h"
#include "my_maths.h"
#include "param.h"
#include <string.h>

struct Cell_vtable cellVTable = {sortAtomsByElement, sortedElementList,
                                 cellExport, seedExport};
struct Cell_textFunc cellTextTable = {cellWriteBlock};

/* Cell starts */
Cell *createCell(Lattice *lat, HashNode *table)
{
    Cell *new = malloc(sizeof(Cell));
    new->lattice = lat;
    new->lattice->vtable->rotate_to_standard_orientation(new->lattice);
    new->atomSorted = false;
    new->destroy = destroyCell;
    new->vtable = &cellVTable;
    new->textTable = &cellTextTable;
    int elmNums = 0;
    new->elmLists = new->vtable->sortElmList(new, &elmNums);
    new->elmNums = elmNums;
    new->lookupTable = table;
    return new;
}

void destroyCell(Cell *self)
{
    self->lattice->vtable->destroy(self->lattice); // Destroy lattice here
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
    Lattice *lat = self->lattice;
    simd_double3x3 latVectors = lat->lattice_vectors;
    simd_double3 a, b, c;
    a = latVectors.columns[0];
    b = latVectors.columns[1];
    c = latVectors.columns[2];
    int bufSize = snprintf(NULL, 0,
                           "%24.18f%24.18f%24.18f\n%24.18f%24.18f%24.18f\n%24."
                           "18f%24.18f%24.18f\n",
                           a.x, a.y, a.z, b.x, b.y, b.z, c.x, c.y, c.z) +
                  1;
    char *resString = malloc(bufSize);
    snprintf(
        resString, bufSize,
        "%24.18f%24.18f%24.18f\n%24.18f%24.18f%24.18f\n%24.18f%24.18f%24.18f\n",
        a.x, a.y, a.z, b.x, b.y, b.z, c.x, c.y, c.z);
    return resString;
}

char *cell_fracCoord_writer(Cell *self)
{
    if (self->atomSorted == false)
        self->vtable->sortAtoms(self);
    Molecule *mol = self->lattice->mol;
    /* Convert */
    simd_double3x3 toFrac = fracCoordMat(self->lattice->lattice_vectors);
    int lineLens[mol->atomNum];
    char **lines = malloc(sizeof(char *) * mol->atomNum);
    int totalLen = 0;
    for (int i = 0; i < mol->atomNum; ++i)
    {
        Atom *cur = mol->atom_arr[i];
        simd_double3 fracCd = simd_mul(toFrac, cur->coord);
        if (i < mol->atomNum - 1)
        {
            lineLens[i] =
                1 + snprintf(NULL, 0, "%3s%20.16f%20.16f%20.16f\n",
                             cur->element, fracCd.x, fracCd.y, fracCd.z);
            lines[i] = malloc(sizeof(char) * lineLens[i]);
            snprintf(lines[i], lineLens[i], "%3s%20.16f%20.16f%20.16f\n",
                     cur->element, fracCd.x, fracCd.y, fracCd.z);
        }
        else
        {
            Atom *metalAtom = self->lattice->mol->vtable->get_atom_by_Id(
                self->lattice->mol, self->lattice->metal_site_id);
            HashNode *metalNode =
                find_item_by_str(self->lookupTable, metalAtom->element);
            ElmInfo *metal = (ElmInfo *)metalNode->val;
            if (metal->spin)
            {
                lineLens[i] =
                    1 + snprintf(NULL, 0,
                                 "%3s%20.16f%20.16f%20.16f SPIN=%14.10f\n",
                                 cur->element, fracCd.x, fracCd.y, fracCd.z,
                                 (double)metal->spin);
                lines[i] = malloc(sizeof(char) * lineLens[i]);
                snprintf(lines[i], lineLens[i],
                         "%3s%20.16f%20.16f%20.16f SPIN=%14.10f\n",
                         cur->element, fracCd.x, fracCd.y, fracCd.z,
                         (double)metal->spin);
            }
            else
            {
                lineLens[i] =
                    1 + snprintf(NULL, 0, "%3s%20.16f%20.16f%20.16f\n",
                                 cur->element, fracCd.x, fracCd.y, fracCd.z);
                lines[i] = malloc(sizeof(char) * lineLens[i]);
                snprintf(lines[i], lineLens[i], "%3s%20.16f%20.16f%20.16f\n",
                         cur->element, fracCd.x, fracCd.y, fracCd.z);
            }
        }
        totalLen += lineLens[i];
    }
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
    char line[] = "FIX_ALL_CELL : true\n\nFIX_COM : false\n";
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
        char *key = self->elmLists[i];
        if (strlen(key) > 2)
        {
            printf("Key error when SPECIES_MASS: %s\n", key);
        }
        const char format[] = "%8s%18.10f\n";
        HashNode *elmNode = find_item_by_str(self->lookupTable, key);
        ElmInfo *elmItem = (ElmInfo *)(elmNode->val);
        lineLens[i] =
            1 + snprintf(NULL, 0, format, self->elmLists[i], elmItem->mass);
        tmpLines[i] = malloc(lineLens[i]);
        snprintf(tmpLines[i], lineLens[i], format, self->elmLists[i],
                 elmItem->mass);
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
        char *key = self->elmLists[i];
        if (strlen(key) > 2)
        {
            printf("Key error in SPECIES_POT: %s\n", key);
        }
        HashNode *elmNode = find_item_by_str(self->lookupTable, key);
        ElmInfo *elmInfo = (ElmInfo *)(elmNode->val);
        lineLens[i] =
            1 + snprintf(NULL, 0, format, self->elmLists[i], elmInfo->potFile);
        tmpLines[i] = malloc(lineLens[i]);
        snprintf(tmpLines[i], lineLens[i], format, self->elmLists[i],
                 elmInfo->potFile);
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
        const char *key = self->elmLists[i];
        if (strlen(key) > 2)
        {
            printf("Key error in SPECIES_POT: %s\n", key);
        }
        HashNode *elmNode = find_item_by_str(self->lookupTable, key);
        ElmInfo *elmInfo = (ElmInfo *)(elmNode->val);
        lineLens[i] =
            1 + snprintf(NULL, 0, format, elmInfo->name, elmInfo->LCAO);
        tmpLines[i] = malloc(lineLens[i]);
        snprintf(tmpLines[i], lineLens[i], format, elmInfo->name,
                 elmInfo->LCAO);
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

void cellExport(Cell *self)
{
    /* Filepath processing */
    char *stemName = self->lattice->mol->name;
    char *exportDir = self->lattice->vtable->exportDir(self->lattice);
    createDirectory(exportDir);
    char *fileName;
    char *DOSfileName;
    /* Differ for DOS or not */
    int fileNameLen = 1 + snprintf(NULL, 0, "%s%s.cell", exportDir, stemName);
    fileName = malloc(fileNameLen);
    snprintf(fileName, fileNameLen, "%s%s.cell", exportDir, stemName);
    int DOSfileNameLen =
        1 + snprintf(NULL, 0, "%s%s_DOS.cell", exportDir, stemName);
    DOSfileName = malloc(DOSfileNameLen);
    snprintf(DOSfileName, DOSfileNameLen, "%s%s_DOS.cell", exportDir, stemName);
    /* Release malloc'd memory exportDir */
    free(exportDir);
    /* Ready to write */
    FILE *writeTo = fopen(fileName, "w");
    FILE *writeDOS = fopen(DOSfileName, "w");
    free(fileName);
    free(DOSfileName);
    /* Get,Write and Free strings for each section */
    char *latVec = self->textTable->blockWriter(self, "LATTICE_CART",
                                                cell_latticeVector_writer);
    fputs(latVec, writeTo);
    fputs(latVec, writeDOS);
    free(latVec);
    char *fracCoord = self->textTable->blockWriter(self, "POSITIONS_FRAC",
                                                   cell_fracCoord_writer);
    fputs(fracCoord, writeTo);
    fputs(fracCoord, writeDOS);
    free(fracCoord);
    char *BS_kPoints = self->textTable->blockWriter(self, "BS_KPOINTS_LIST",
                                                    cell_kPointsList_writer);
    fputs(BS_kPoints, writeDOS);
    free(BS_kPoints);
    char *kPoints = self->textTable->blockWriter(self, "KPOINTS_LIST",
                                                 cell_kPointsList_writer);
    fputs(kPoints, writeTo);
    fputs(kPoints, writeDOS);
    free(kPoints);
    char *misc = cell_miscOptions_writer(self);
    fputs(misc, writeTo);
    fputs(misc, writeDOS);
    free(misc);
    char *ionicConstraint = self->textTable->blockWriter(
        self, "IONIC_CONSTRAINTS", cell_ionicConstraints_writer);
    fputs(ionicConstraint, writeTo);
    fputs(ionicConstraint, writeDOS);
    free(ionicConstraint);
    char *efield = self->textTable->blockWriter(self, "EXTERNAL_EFIELD",
                                                cell_externalEfield_writer);
    fputs(efield, writeTo);
    fputs(efield, writeDOS);
    free(efield);
    char *ePressure = self->textTable->blockWriter(
        self, "EXTERNAL_PRESSURE", cell_externalPressure_writer);
    fputs(ePressure, writeTo);
    fputs(ePressure, writeDOS);
    free(ePressure);
    char *spMass = self->textTable->blockWriter(self, "SPECIES_MASS",
                                                cell_speciesMass_writer);
    fputs(spMass, writeTo);
    fputs(spMass, writeDOS);
    free(spMass);
    char *spPot = self->textTable->blockWriter(self, "SPECIES_POT",
                                               cell_speciesPot_writer);
    fputs(spPot, writeTo);
    fputs(spPot, writeDOS);
    free(spPot);
    char *spLCAO = self->textTable->blockWriter(self, "SPECIES_LCAO_STATES",
                                                cell_speciesLCAOstates_writer);
    fputs(spLCAO, writeTo);
    fputs(spLCAO, writeDOS);
    free(spLCAO);
    /* Close file release pointer */
    fclose(writeTo);
    fclose(writeDOS);
}

void seedExport(Cell *self)
{
    write_kptaux(self);
    write_param(self);
    write_trjaux(self);
    write_pbsScript(self);
    write_SMCastepExtension(self);
    /* copy_potentials(self); */
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
    Molecule *mol = self->lattice->mol;
    Atom **atomArray = mol->atom_arr;
    qsort(atomArray, mol->atomNum, sizeof(Atom *), atomCmp);
    mol->atom_arr = atomArray;
    self->atomSorted = true;
    self->lattice->metal_site_id = mol->atomNum;
}

char **sortedElementList(Cell *self, int *returnSize)
{
    if (self->atomSorted == false)
        self->vtable->sortAtoms(self);
    Molecule *mol = self->lattice->mol;
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
