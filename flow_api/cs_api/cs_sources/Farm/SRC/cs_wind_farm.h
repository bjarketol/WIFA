#ifndef READ_SOURCES_H
#define READ_SOURCES_H

#include "cs_headers.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
  char         *WT_file_name;     /* file name */
  cs_lnum_t     n_WT;             /* number of turbines  */
  cs_lnum_t     n_WT_types;             /* number of turbines  */
  cs_real_t     total_power;
  cs_real_3_t  *WT_coords;
  cs_real_t    *WT_diameters;
  cs_lnum_t    *WT_types;
  cs_real_t    *WT_volume;
  cs_real_3_t  *WT_uvect;
  cs_real_t    *WT_u;
  cs_real_t    *WT_rho;
  cs_real_t    *WT_hub_u;
  cs_real_t    *WT_dir;
  cs_real_t    *WT_ra_dir;
  cs_real_t    *WT_ctstar;
  cs_real_t    *WT_cpstar;
  cs_real_t    *WT_thrust;
  cs_real_t    *WT_power;
} cs_wind_farm;

typedef struct _cs_ctstar_table cs_ctstar_table;
typedef struct _cs_cpstar_table cs_cpstar_table;

/* ============================================================================
* Static global variable
*============================================================================*/

/*------------Pointer to the emission structure------------------------------*/
extern cs_wind_farm *cs_glob_wind_farm;

/*----------------------------------------------------------------------------*/
/*!
 * \file  cs_wind_farm.c
 * \fn    create_wind_farm
 * \brief Allocate memory and create a global static variable cs_glob_wind_farm
 **
 * Allocate memory for a global static variable cs_glob_wind farm and fill it
 * with a values taken from WT_file_name.
 *
 * \param[in]   WT_file_name             file name of emissions
 * \return      Pointer to created structure or NULL in case some error
 */
/*----------------------------------------------------------------------------*/

cs_wind_farm *
cs_wind_farm_create_from_file(const char *WT_file_name);

/*----------------------------------------------------------------------------*/
/*!
 * \fn    cs_wind_farm_free
 * \brief Free memory for the variable of type cs_wind_farm
 *
 * \param[in]   cs_glob_wind_farm   Variable to free memory
 */
/*----------------------------------------------------------------------------*/
cs_wind_farm *
cs_wind_farm_free(cs_wind_farm *cs_glob_wind_farm);

/*----------------------------------------------------------------------------*/
/*!
 * \fn    cs_wind_farm_compute_source_term
 * \brief Computes source term for wind farm
 *
 */
/*----------------------------------------------------------------------------*/
void
cs_wind_farm_compute_source_term(void);

/*----------------------------------------------------------------------------*/
/*!
 * \fn    cs_wind_farm_write
 * \brief Writes output file for the wind farm
 *
 */
/*----------------------------------------------------------------------------*/
void
cs_wind_farm_write(void);

void
rotate2(cs_real_t *coords,
        cs_real_t tilt_angle,
        cs_real_t yaw_angle,
        cs_real_t roll_angle,
        cs_real_t *mean_coords,
        cs_real_t *rotated_coords);

cs_real_t
interp_ct_or_cp(cs_real_t *ct_or_cp_values, cs_real_t *ct_or_cp_speeds, cs_lnum_t number_of_values, cs_real_t disk_velocity);

#ifdef __cplusplus0
}
#endif

#endif // READ_SOURCES_H
