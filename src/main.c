#include "atom.h"
#include "database/ads_database.h"
#include "database/database.h"
#include "database/lattice_database.h"
#include "misc.h"
#include "my_maths.h"
#include "parser.h"
#include <time.h>

void test_hash()
{
    HashNode *adsTable = init_adsInfoTable();
    HashNode *item = find_item_by_str(adsTable, "C2H4");
    printf("%s, %d\n", ((AdsorbateInfo *)item->val)->name,
           ((AdsorbateInfo *)item->val)->coordAtomNum);
    Adsorbate *ads = parse_adsorbate_from_file(
        "./C2_pathways_ads/ethylene_path/C2H4.msi", "C2H4", item->val);
    ads->vtable->make_upright(ads);
    for (int i = 0; i < ads->mol->atomNum; ++i)
    {
        char *text =
            ads->mol->atom_arr[i]->vtable->export_text(ads->mol->atom_arr[i]);
        printf("%s\n", text);
        free(text);
    }
    ads->vtable->destroy(ads);
    delete_all(&adsTable);
}

int main(int argc, char *argv[])
{
    test_hash();
    HashNode *elmTable = init_ElmInfoTable();
    delete_all(&elmTable);
    return 0;
}
