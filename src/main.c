#include "assemble.h"
#include "atom.h"
#include "cell.h"
#include "database/ads_database.h"
#include "database/database.h"
#include "database/lattice_database.h"
#include "misc.h"
#include "my_maths.h"
#include "parser.h"
#include <time.h>

void test_hash()
{
    AdsTableYAML *adsTableYAML = load_adsTableYAML();
    HashNode *adsTable = init_adsInfoTable(adsTableYAML);
    HashNode *item = find_item_by_str(adsTable, "C2H4");
    AdsInfo *adsInfo = (AdsInfo *)(item->val);
    printf("%s, %d\n", adsInfo->name, ((AdsInfo *)item->val)->coordAtoms_count);
    Adsorbate *ads = parse_adsorbate_from_file(
        "./C2_pathways_ads/ethylene_path/C2H4.msi", "C2H4", item->val);
    ads->vtable->make_upright(ads);
    ads->vtable->destroy(ads);
    delete_all(&adsTable);
    adsTableYAML->destroy(&adsTableYAML);
}

void test_lat()
{
    Lattice *lat =
        parse_lattice_from_file("./msi_models/3d/SAC_GDY_V.msi", "SAC_GDY_V");
    for (int i = 0; i < 3; ++i)
    {
        printf("%f, %f, %f\n", lat->lattice_vectors.columns[i].x,
               lat->lattice_vectors.columns[i].y,
               lat->lattice_vectors.columns[i].z);
    }
    lat->vtable->rotate_to_standard_orientation(lat);
    printf("Rotate to standard orientation\n");
    for (int i = 0; i < 3; ++i)
    {
        printf("%f, %f, %f\n", lat->lattice_vectors.columns[i].x,
               lat->lattice_vectors.columns[i].y,
               lat->lattice_vectors.columns[i].z);
    }
    double alpha = simd_vector_angle(lat->lattice_vectors.columns[1],
                                     lat->lattice_vectors.columns[2]);
    printf("%f, %f\n", alpha, cos(alpha));
    ElmTableYAML *elmTableYAML = load_elmTableYAML();
    HashNode *elmTable = init_ElmInfoTable(elmTableYAML);
    HashNode *metal = find_item_by_str(elmTable, "Co");
    printf("%s\n", (char *)metal->key);
    ElmInfo *metal_info = (ElmInfo *)metal->val;
    printf("After find %p\n", metal->val);
    printf("%s,%d\n", ((ElmInfo *)metal->val)->name, metal_info->atomicNum);
    lat->vtable->modify_metal(lat, metal_info->name, metal_info->atomicNum);
    Cell *cell = createCell(lat, elmTable);
    char *cell_pot = cell->textTable->blockWriter(cell, "SPECIES_POT",
                                                  cell_speciesPot_writer);
    printf("%s\n", cell_pot);
    free(cell_pot);
    cell->destroy(cell);
    delete_all(&elmTable);
    elmTableYAML->destroy(&elmTableYAML);
}

void test_assemble()
{
    Lattice *lat =
        parse_lattice_from_file("./msi_models/3d/SAC_GDY_V.msi", "SAC_GDY_V");
    AdsTableYAML *adsTableYAML = load_adsTableYAML();
    HashNode *adsTable = init_adsInfoTable(adsTableYAML);
    AdsInfo *c2h4 = (AdsInfo *)find_item_by_str(adsTable, "C2H4")->val;
    Adsorbate *ads = parse_adsorbate_from_file(
        "./C2_pathways_ads/ethylene_path/C2H4.msi", "C2H4", c2h4);
    ElmTableYAML *elmTableYAML = load_elmTableYAML();
    HashNode *elmTable = init_ElmInfoTable(elmTableYAML);
    ElmInfo *Tm = find_item_by_str(elmTable, "Tm")->val;
    lat->vtable->modify_metal(lat, "Tm", Tm->atomicNum);
    Lattice *assembled =
        Add_mol_to_lattice(lat, ads, ads->taskLists->tasks[0][0],
                           ads->taskLists->tasks[0][1], 1.4);
    assembled->vtable->export_msi(assembled);
    Cell *assembled_cell = createCell(assembled, elmTable);
    assembled_cell->vtable->exportCell(assembled_cell);
    lat->vtable->destroy(lat);
    ads->vtable->destroy(ads);
    assembled_cell->destroy(assembled_cell);
    delete_all(&adsTable);
    delete_all(&elmTable);
    destroy_adsTableYAML(&adsTableYAML);
    elmTableYAML->destroy(&elmTableYAML);
}

int main(int argc, char *argv[])
{
    test_hash();
    test_assemble();
    return 0;
}
