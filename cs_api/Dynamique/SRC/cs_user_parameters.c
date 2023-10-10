/*============================================================================
 * User functions for input of calculation parameters.
 *============================================================================*/

/* VERS */

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2021 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <math.h>
#include <string.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

/*----------------------------------------------------------------------------
 * PLE library headers
 *----------------------------------------------------------------------------*/

#include <ple_coupling.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Select physical model options, including user fields.
 *
 * This function is called at the earliest stages of the data setup,
 * so field ids are not available yet.
 */
/*----------------------------------------------------------------------------*/

void
cs_user_model(void)
{

  cs_wall_functions_t *wf = cs_get_glob_wall_functions();
  wf->iwalfs = CS_WALL_F_S_MONIN_OBUKHOV;
  wf->iwallf = CS_WALL_F_2SCALES_SMOOTH_ROUGH;

  /* Atmospheric module options */
  cs_glob_physical_model_flag[CS_ATMOSPHERIC] = CS_ATMO_DRY;

  /* Read the meteo file (1) or impose directly the input values to compute it
  * in Code_Saturne (2) */  
  cs_glob_atmo_option->meteo_profile = 2;

  /* Inverse LMO length */
  cs_glob_atmo_option->meteo_dlmo = cs_notebook_parameter_value_by_name("Lmoinv");
  /* Large scale roughness */
  cs_glob_atmo_option->meteo_z0 = cs_notebook_parameter_value_by_name("z0");
  /* Elevation for reference velocity */
  cs_glob_atmo_option->meteo_zref = cs_notebook_parameter_value_by_name("zref");
  /* Friction velocity */
  cs_glob_atmo_option->meteo_uref = cs_notebook_parameter_value_by_name("ureff");
  /* Velocity direction */
  cs_glob_atmo_option->meteo_angle = cs_notebook_parameter_value_by_name("teta");
  /* Ground temperature */ 
  cs_glob_atmo_option->meteo_t0 = cs_notebook_parameter_value_by_name("t0");
  /* X coordinate of the centre of the domain in Lambert 93 */
  cs_glob_atmo_option->x_l93 = cs_notebook_parameter_value_by_name("XL93");
  /* Y coordinate of the centre of the domain in Lambert 93 */
  cs_glob_atmo_option->y_l93 = cs_notebook_parameter_value_by_name("YL93");
  /* Option to compute ground elevation in the domain */
  cs_glob_atmo_option->compute_z_ground = false;

  /* Automatic open boundary conditions
   *   1: meteo mass flow rate is imposed with a constant large scale
   *      pressure gradient
   *   2: same plus velocity profile imposed at ingoing faces
   */
  cs_glob_atmo_option->open_bcs_treatment = 0;

  /* Time of the simulation (for radiative model or chemistry)
   * syear:  starting year
   * squant: starting quantile
   * shour:  starting hour (UTC)
   * smin:   starting minute
   * ssec:   starting second
   */
  cs_glob_atmo_option->syear = 2020;
  cs_glob_atmo_option->squant = 1;
  cs_glob_atmo_option->shour = 1;
  cs_glob_atmo_option->smin = 0;
  cs_glob_atmo_option->ssec = 0.;

  /* /\* Add boundary properties for roughness read in tif file *\/ */
  /* cs_parameters_add_property("gdalroughness", */
  /*                            1, */
  /*                            CS_MESH_LOCATION_BOUNDARY_FACES); */

  /* cs_parameters_add_property("roughness_zone", */
  /*                            1, */
  /*                            CS_MESH_LOCATION_BOUNDARY_FACES); */

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define or modify general numerical and physical user parameters.
 *
 * At the calling point of this function, most model-related most variables
 * and other fields have been defined, so specific settings related to those
 * fields may be set here.
 *
 * At this stage, the mesh is not built or read yet, so associated data
 * such as field values are not accessible yet, though pending mesh
 * operations and some fields may have been defined.
 *
 * \param[in, out]   domain    pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_parameters(cs_domain_t *domain)
{

  /* Example: activate mesh robustness options */
  /*-------------------------------------------*/

  cs_glob_mesh_quantities_flag |= CS_CELL_CENTER_CORRECTION;
  cs_glob_mesh_quantities_flag |= CS_CELL_FACE_CENTER_CORRECTION;
  cs_glob_mesh_quantities_flag |= CS_FACE_RECONSTRUCTION_CLIP;
  cs_glob_mesh_quantities_flag |= CS_CELL_VOLUME_RATIO_CORRECTION;
  cs_glob_mesh_quantities_flag |= CS_BAD_CELLS_WARPED_CORRECTION;
  cs_glob_mesh_quantities_flag |= CS_BAD_CELLS_REGULARISATION;
  cs_glob_mesh_quantities_flag |= CS_FACE_DISTANCE_CLIP;

  /* Uncomment following if global forcing term is applied */
  /*cs_velocity_pressure_param_t *vp_param = cs_get_glob_velocity_pressure_param();
   *vp_param->igpust=0;*/
  
  /* Modify variables properties */
  cs_equation_param_t *eqp;

  eqp = cs_field_get_equation_param(CS_F_(k));
  eqp->blencv = 0;
  /* Following option might be activated for calculus which do not converge */
  /* eqp->ircflu = 0; */

  eqp = cs_field_get_equation_param(CS_F_(eps));
  eqp->blencv = 0;
  /* Following option might be activated for calculus which do not converge */
  /* eqp->ircflu = 0; */
  
  /* Get the Key for the Inf for the convective scheme */
  int kccmin = cs_field_key_id("min_scalar");
  /* Set the Value for the Sup and Inf of the studied scalar
   * for the Gamma beta limiter of k and epsilon */
  cs_field_set_key_double(CS_F_(eps), kccmin, 0.);
  cs_field_set_key_double(CS_F_(k), kccmin, 0.);


}

/*----------------------------------------------------------------------------*/

END_C_DECLS
