#include "ads_database.h"
#include "atom.h"
#include "misc.h"
#include "my_maths.h"
#include "parser.h"
#include <time.h>

void test_parse_atom(char *msiPath, char *name)
{
    char *body = readWholeFile(msiPath);
    int nums = 0;
    Atom **atoms = get_all_atoms(body, &nums);
    free(body);
    simd_double3 atomCoords[nums];
    time_t t;
    for (int i = 0; i < nums; ++i)
        atomCoords[i] = atoms[i]->coord;
    t = clock();
    /* simd_double3 centroid = simd_centroid_of_points(atomCoords, nums); */
    t = clock() - t;
    double time_taken = ((double)t) / CLOCKS_PER_SEC;
    printf("SIMD centroid routine took %f\n", time_taken);
    printf("%lu\n", sizeof(atomCoords[0]));
    for (int i = 0; i < nums; ++i)
    {
        atoms[i]->vtable->destroy(atoms[i]);
    }
    free(atoms);
}

void test_ads(char *msiPath, char *name)
{
    Adsorbate *ads = parse_molecule_from_file(msiPath, name);
    ads->vtable->make_upright(ads);
    ads->vtable->export_msi(ads, NULL);
    ads->vtable->destroy(ads);
}

void test_hash()
{
    HashNode *adsTable = init_adsInfoTable();
    HashNode *item = find_item_by_str(adsTable, "HOHCCHO");
    printf("%s, %d\n", ((AdsorbateInfo *)item->val)->name,
           ((AdsorbateInfo *)item->val)->coordAtomNum);
    Adsorbate *ads = parse_adsorbate_from_file(
        "./C2_pathways_ads/ethanol_path/HOHCCHO.msi", "HOHCCHO", item->val);
    ads->vtable->make_upright(ads);
    for (int i = 0; i < ads->_mol->atomNum; ++i)
    {
        char *text =
            ads->_mol->atom_arr[i]->vtable->export_text(ads->_mol->atom_arr[i]);
        printf("%s\n", text);
        free(text);
    }
    ads->vtable->destroy(ads);
    delete_all(&adsTable);
}

int main(int argc, char *argv[])
{
    test_parse_atom("./msi_models/3d/SAC_GDY_V.msi", "SAC_GDY_V");
    test_ads("./C2_pathways_ads/ethanol_path/HOHCCHO.msi", "HOHCCHO");
    if (argc == 3)
        test_ads(argv[1], argv[2]);
    test_hash();
    return 0;
}
