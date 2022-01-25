#include "assemble.h"
#include "lattice.h"
#include "matrix.h"
#include "my_maths.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char *append_cd2_mol_name(Lattice *lat, Adsorbate *ads, int c1, int c2)
{
    char *suffix_1 = get_carbon_site_name(c1);
    char *suffix_2 = get_carbon_site_name(c2);
    int suffix_len = 1 + strlen(suffix_1) + 1 +
                     strlen(suffix_2); // _ + suffix1 + _ + suffix2
    int ads_name_len = strlen(ads->_mol->name);
    int lat_name_len = strlen(lat->_mol->name);
    int total_new_name_len = lat_name_len + 1 + ads_name_len + suffix_len;
    char *buffer = malloc(total_new_name_len + 1);
    snprintf(buffer, total_new_name_len + 1, "%s_%s_%s_%s", lat->_mol->name,
             ads->_mol->name, suffix_1, suffix_2);
    return buffer;
}

void ads_align_with_carbon_chain(Lattice *lat, Adsorbate *ads)
{
    /* Init the molecule to be upright
     * Except C2H4, prefer it to be planar
     */
    if (strcmp(ads->_mol->name, "C2H4"))
        ads->ads_vtable->make_upright(ads);
    /* Get the stem vector to determine the direction of the adsorbate*/
    Matrix *ads_stem_vec = ads->ads_vtable->get_stem_vector(ads);
    // Get the carbon chain vector
    Matrix *carbon_chain_vec = lat->vtable->get_carbon_chain_vector(lat);
    // Obtain the rotation matrix for the stem rotate to carbon chain, with the
    // axis being the cross product of the stem and chain vectors
    Matrix *rot_mat = rotate_u_to_v(ads_stem_vec, carbon_chain_vec);
    // Apply rotation
    ads->_mol->vtable->apply_transformation(ads->_mol, rot_mat,
                                            rotate_around_origin);
    // Destroy and Free used matrix
    destroy_matrix(ads_stem_vec);
    destroy_matrix(carbon_chain_vec);
    destroy_matrix(rot_mat);
    free(ads_stem_vec);
    free(carbon_chain_vec);
    free(rot_mat);
}

void ads_align_with_carbon_metal(Lattice *lat, Adsorbate *ads, int carbonId)
{
    if (carbonId != 41 && carbonId != 42)
    {
        printf("The target carbon site %d is not valid\n", carbonId);
        return;
    }
    /* Init the molecule to be upright
     * Except C2H4, prefer it to be planar
     */
    if (strcmp(ads->_mol->name, "C2H4"))
        ads->ads_vtable->make_upright(ads);
    /* Get the stem vector to determine the direction of the adsorbate */
    Matrix *ads_stem_vec = ads->ads_vtable->get_stem_vector(ads);
    /* Get the carbon to metal vector */
    Matrix *carbon_metal_vec =
        lat->vtable->get_carbon_metal_vector(lat, carbonId);
    /* Obtain the required transformation (rotation) matrix */
    Matrix *rot_mat = rotate_u_to_v(ads_stem_vec, carbon_metal_vec);
    /* Apply transformation */
    Molecule *mol = ads->_mol;
    mol->vtable->apply_transformation(mol, rot_mat, rotate_around_origin);

    /* Destroy and free used matrix */
    destroy_matrix(ads_stem_vec);
    destroy_matrix(carbon_metal_vec);
    destroy_matrix(rot_mat);
    free(ads_stem_vec);
    free(carbon_metal_vec);
    free(rot_mat);
}

Lattice *Add_cd2_mol_to_carbon_chain(Lattice *lat, Adsorbate *ads, int c1,
                                     int c2)
{
    // Init the adsorbate to be aligned with the carbon chain
    ads_align_with_carbon_chain(lat, ads);
    // Get the center coordinates of the two carbon sites
    double *center_coords_of_carbon_sites = /* malloced */
        lat->_mol->vtable->get_centroid_ab(lat->_mol, c1, c2);
    double *center_coords_of_cd_sites = ads->_mol->vtable->get_centroid_ab(
        ads->_mol, ads->coordAtomIds[0], ads->coordAtomIds[1]);
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
    char *newName = append_cd2_mol_name(lat, ads, c1, c2);
    Lattice *newModel = lat->vtable->attach_molecule(lat, ads, newName);
    free(newName);
    return newModel;
}

Lattice *Add_cd2_mol_to_carbon_metal(Lattice *lat, Adsorbate *ads, int carbonId)
{
    /* Init the adsorbate to the correct direction */
    ads_align_with_carbon_metal(lat, ads, carbonId);
    Molecule *latMol = lat->_mol, *adsMol = ads->_mol;
    double *center_of_carbon_metal =
        latMol->vtable->get_centroid_ab(latMol, carbonId, 73);
    double *center_of_cd_sites = adsMol->vtable->get_centroid_ab(
        adsMol, ads->coordAtomIds[0], ads->coordAtomIds[1]);
    Matrix *trans_m =
        translate_mat_a_to_b(center_of_cd_sites, center_of_carbon_metal);
    /* Set the hanging height  */
    trans_m->value[2][3] += 1.5;
    adsMol->vtable->apply_transformation(adsMol, trans_m, translate_a_to_b);

    /* Clean memory */
    destroy_matrix(trans_m);
    free(trans_m);
    free(center_of_carbon_metal);
    free(center_of_cd_sites);
    char *newName = append_cd2_mol_name(lat, ads, carbonId, 73);
    Lattice *newModel = lat->vtable->attach_molecule(lat, ads, newName);
    free(newName);
    return newModel;
}
