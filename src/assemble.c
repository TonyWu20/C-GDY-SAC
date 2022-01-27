#include "assemble.h"
#include "lattice.h"
#include "matrix.h"
#include "my_maths.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char *append_mol_name(Lattice *lat, Adsorbate *ads, int c1, int c2)
{
    char *suffix_1 = get_carbon_site_name(c1);
    char *suffix_2 = (c2 != NULLSITE) ? get_carbon_site_name(c2) : NULL;
    int suffix_len = (c2 != NULLSITE)
                         ? 1 + strlen(suffix_1) + 1 + strlen(suffix_2)
                         : 1 + strlen(suffix_1); // _ + suffix1 + _ + suffix2
    int ads_name_len = strlen(ads->_mol->name);
    int lat_name_len = strlen(lat->_mol->name);
    int total_new_name_len = lat_name_len + 1 + ads_name_len + suffix_len;
    char *buffer = malloc(total_new_name_len + 1);
    if (c2 != NULLSITE)
    {
        char namePattern[] = "%s_%s_%s_%s";
        snprintf(buffer, total_new_name_len + 1, namePattern, lat->_mol->name,
                 ads->_mol->name, suffix_1, suffix_2);
    }
    else
    {
        char namePattern[] = "%s_%s_%s";
        snprintf(buffer, total_new_name_len + 1, namePattern, lat->_mol->name,
                 ads->_mol->name, suffix_1);
    }
    return buffer;
}

int init_ads_direction(Lattice *lat, Adsorbate *ads, int id_from, int id_to)
{
    /* Init the molecule to be upright
     * Except C2H4, prefer it to be planar
     */
    if ((id_from == 73 && id_to != 41 && id_to != 42) ||
        (id_to == 73 && id_from != 41 && id_from != 42))
    {
        printf("When one site is metal, the other given site is not a valid "
               "carbon site\n");
        return -1;
    }
    if (strcmp(ads->_mol->name, "C2H4"))
        ads->ads_vtable->make_upright(ads);
    /* Get the stem vector to determine the direction of the adsorbate*/
    Matrix *ads_stem_vec = ads->ads_vtable->get_stem_vector(ads);
    /* Get the target direction vector */
    Matrix *direction_vec =
        lat->_mol->vtable->get_vector_ab(lat->_mol, id_from, id_to);
    Matrix *rot_mat = rotate_u_to_v(ads_stem_vec, direction_vec);
    ads->_mol->vtable->apply_transformation(ads->_mol, rot_mat,
                                            rotate_around_origin);
    /* Clean memory */
    destroy_matrix(ads_stem_vec);
    destroy_matrix(direction_vec);
    destroy_matrix(rot_mat);
    free(ads_stem_vec);
    free(direction_vec);
    free(rot_mat);
    return 0;
}

Lattice *Add_mol_to_lattice(Lattice *lat, Adsorbate *ads, int c1, int c2)
{
    // Init the adsorbate to be aligned with the carbon chain
    if (c2 != NULLSITE || ads->coordAtomNum == 2)
        init_ads_direction(lat, ads, c1, c2);
    else
        init_ads_direction(lat, ads, 41, 42);
    int cd_1 = ads->coordAtomIds[0];
    int cd_2 = (c2 != NULLSITE || ads->coordAtomNum == 2) ? ads->coordAtomIds[1]
                                                          : cd_1;
    int carbon_1 = c1;
    int carbon_2 = (c2 != NULLSITE || ads->coordAtomNum == 2) ? c2 : c1;
    // Get the center coordinates of the two carbon sites
    double *center_coords_of_carbon_sites = /* malloced */
        lat->_mol->vtable->get_centroid_ab(lat->_mol, carbon_1, carbon_2);
    double *center_coords_of_cd_sites =
        ads->_mol->vtable->get_centroid_ab(ads->_mol, cd_1, cd_2);
    Matrix *trans_m = translate_mat_a_to_b(center_coords_of_cd_sites,
                                           center_coords_of_carbon_sites);
    /* Set the hanging height on the GDY plane */
    trans_m->value[2][3] += 1.5;
    /* Apply transformation */
    ads->_mol->vtable->apply_transformation(ads->_mol, trans_m,
                                            translate_a_to_b);
    /* Clean memory */
    destroy_matrix(trans_m);
    free(trans_m);
    free(center_coords_of_carbon_sites);
    free(center_coords_of_cd_sites);
    char *newName = append_mol_name(lat, ads, c1, c2);
    Lattice *newModel = lat->vtable->attach_molecule(lat, ads, newName);
    free(newName);
    return newModel;
}
