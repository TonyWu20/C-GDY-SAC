#include "main.h"
#include "maths/MyMaths.h"
#include "msiParser/msiParser.h"
#include "msiWriter/msiWriter.h"
#include <stdio.h>

void init_mol_direction(FILE *file);
void test_latMsi_writer(FILE *file);
void test_attach_mol(FILE *mol, FILE *lat);
int main(int argc, char *argv[])
{
    char *fileName = argv[1];
    FILE *file = fopen(fileName, "r");
    if (file == NULL)
    {
        printf("Cannot open file.\n");
        return 1;
    }
    char *latName = argv[2];
    FILE *flat = fopen(latName, "r");
    if (file == NULL)
    {
        printf("Cannot open file.\n");
        return 1;
    }
    test_attach_mol(file, flat);
}

void init_mol_direction(FILE *file)
{
    MOLECULE *mol;
    mol = parseMol(file);
    fclose(file);
    init_xz_plane(mol);
    MSI_FILE *rotated_mol;
    rotated_mol = build_MolMsi(mol);
    free(mol);
    for (int i = 0; i < rotated_mol->ItemNum; i++)
    {
        printf("%s", rotated_mol->lines[i]);
    }
    free_MSI_MOL(rotated_mol);
}

void test_latMsi_writer(FILE *file)
{
    BASE_LATTICE *lat;
    lat = parseBase(file);
    fclose(file);
    MSI_FILE *lat_file;
    lat_file = build_LatMsi(lat);
    free(lat);
    free_MSI_LAT(lat_file);
}

void test_attach_mol(FILE *fmol, FILE *flat)
{
    MOLECULE *mol = parseMol(fmol);
    BASE_LATTICE *lat = parseBase(flat);
    BASE_LATTICE *added_lat = init_adsorbed_lat(lat, mol);
    fclose(fmol);
    fclose(flat);
    init_xz_plane(mol);
    align_carbon_chain(mol, lat->carbon_chain_vec);
    placeMol(mol, lat, 39, added_lat);
    MSI_FILE *new_flat = build_LatMsi(added_lat);
    for (int i = 0; i < new_flat->ItemNum; ++i)
    {
        printf("%s", new_flat->lines[i]);
    }
    free(mol);
    free(lat);
    free(added_lat);
    free_MSI_LAT(new_flat);
}
