/*============================================================================
 * Additional user-defined source terms for variable equations.
 *============================================================================*/

/* VERS */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2022 EDF S.A.

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
 * \brief Function called at each time step to define physical properties.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 * \param[in]       f_id     field id of the variable
 * \param[out]      st_exp   explicit source term
 * \param[out]      st_imp   implicit part of the source term
 */
/*----------------------------------------------------------------------------*/

void
cs_user_source_terms(cs_domain_t  *domain,
                     int           f_id,
                     cs_real_t    *st_exp,
                     cs_real_t    *st_imp)
{
  const cs_real_t *cell_vol = domain->mesh_quantities->cell_vol;
  const cs_real_3_t *cell_cen
    = (const cs_real_3_t *)domain->mesh_quantities->cell_cen;

  const cs_mesh_t *m = domain->mesh;
  
  const cs_real_t *cpro_rom = CS_F_(rho)->val;
  const cs_real_3_t *vel = CS_F_(vel)->val;

  const cs_field_t *fld = cs_field_by_id(f_id);
 
  /*physical variables for density variation */
  cs_real_t *source_term_x = cs_field_by_name("source_term_x")->val;
  cs_real_t *source_term_y = cs_field_by_name("source_term_y")->val;
  cs_real_t *st_coeff = cs_field_by_name("source_term_coeff")->val;
  cs_lnum_t   cell_id = 0;
  double domain_volume = 0.;
  
  cs_time_step_t *ts = cs_get_glob_time_step();


  /******************************************
   * READ Wind Turbines from coordinates file
   *******************************************/
  cs_real_t WT_d = 136.0;  //TODO: read as notebook param
  cs_real_t WT_radius = WT_d/2.;
  cs_real_t wind_dir=cs_glob_atmo_option->meteo_angle;
  
  //file name and tables
  char WT_file_name[1024];
  sprintf(WT_file_name,"placement_turbines.csv");
  //calculate number of lines in the file
  FILE* stream = fopen(WT_file_name, "r");
  char line[1024];
  char *tok0;
  const char sep0[4] = ";";
  size_t n_WT = 0;
  //read header
  //fgets(line, 1024, stream);
  while (fgets(line, 1024, stream) != NULL){
    // try to parse line
    tok0 = strtok(line, sep0);
    if (tok0 != NULL)
      n_WT++;
  }
  fclose(stream);

  //allocate memory for all members of emission structure
  cs_real_t *WT_x_coords=NULL;
  cs_real_t *WT_y_coords=NULL;
  cs_real_t *WT_z_coords=NULL;  
  BFT_MALLOC(WT_x_coords, n_WT, cs_real_t);
  BFT_MALLOC(WT_y_coords, n_WT, cs_real_t);
  BFT_MALLOC(WT_z_coords, n_WT, cs_real_t);

   //read and parse this file
  stream = fopen(WT_file_name, "r");
  char *tok;
  const char sep[4] = ",";
  size_t line_count = -1;
  size_t i=-1;

  while (fgets(line, 1024, stream) != NULL){
    i++;
    if (++line_count >= n_WT)
      break;

    // parse line
    tok = strtok(line, sep);
    if (tok == NULL)
      continue;
    WT_x_coords[i]=atof(tok);
    
    tok = strtok(NULL, sep);
    if (tok == NULL)
      continue;
    WT_y_coords[i]=atof(tok);

    tok = strtok(NULL, sep);
    if (tok == NULL)
      continue;
    WT_z_coords[i]=atof(tok);
  }
  fclose(stream);
  /******************************************
   * ACTUATOR DISK SOURCE TERMS
   *******************************************/
  char name[128];
  
  /* For velocity
   * ============ */

  if (fld == CS_F_(vel)) {

    const cs_equation_param_t *eqp = cs_field_get_equation_param(CS_F_(vel));
    if (eqp->verbosity == 1) {
      bft_printf(" User source terms for variable %s \n ",
                 cs_field_get_label(CS_F_(vel)));
    }
    const cs_real_t dpi = atan(1.0)*4.0;
    cs_real_t WT_surf = dpi*pow(WT_radius,2);
    cs_real_t WT_volume;
    cs_real_t ra_wind_dir= (270.0-wind_dir)*dpi/180.0; //radian absolute wind_dir
    cs_real_t Ct;
    //
    cs_real_t u_sonde, u_cell, u_x, u_y, u_z;
    cs_real_3_t *_st_exp = (cs_real_3_t *)st_exp;
    cs_real_33_t *_st_imp = (cs_real_33_t *)st_imp;
    
    cs_real_t rotor_rdr=1.0; //ratio for rotor radius relative to turbine radius
    cs_real_t rotor_thickness_rdr=0.1; //ratio for rotor thickness relative to turbine diameter
    cs_real_t AD_mesh_cell_size = 5.0;
    cs_real_t cancel_dist = 0.5*AD_mesh_cell_size;
    cs_real_t select_dist = AD_mesh_cell_size;
    /******************LOOP ON TURBINES*************************/
    for (cs_lnum_t WT_count=0; WT_count < n_WT; WT_count ++){
      //get Actuator Disk zone by name
      sprintf(name,"turbine_%d",WT_count+1);
      cs_zone_t *z = cs_volume_zone_by_name_try(name);
      
      /* u_sonde=0.0; */
      /* u_x=0.0; */
      /* u_y=0.0; */
      /* u_z=0.0; */
      /* u_cell=0.0;	 */
      /* for (cs_lnum_t nc =0; nc < z->n_elts; nc++){ */
      /* 	cell_id = z->elt_ids[nc]; */
      /* 	u_x = u_x + vel[cell_id][0]*cell_vol[cell_id]/domain_volume; */
      /* 	u_y = u_y + vel[cell_id][1]*cell_vol[cell_id]/domain_volume; */
      /* 	u_z = u_z + vel[cell_id][2]*cell_vol[cell_id]/domain_volume; */
      /* } */
      /* u_sonde=sqrt(pow(u_x,2)+pow(u_y,2)+pow(u_z,2)); */

      //for testing
      u_sonde=5.0;
      u_x = u_sonde*cos(ra_wind_dir);
      u_y = u_sonde*sin(ra_wind_dir);

      //TODO : read Ct from file
      Ct=0.762;
      //
      WT_volume=0.0;      
      for (cs_lnum_t nc = 0; nc < z->n_elts; nc++) {
	cs_lnum_t c_id = z->elt_ids[nc];
	//
	cs_real_t x_dist,y_dist,z_dist;
	x_dist = cell_cen[c_id][0] - WT_x_coords[WT_count];
	y_dist = cell_cen[c_id][1] - WT_y_coords[WT_count];
	z_dist = cell_cen[c_id][2] - WT_z_coords[WT_count];
	cs_real_t cell_radius = sqrt(pow(y_dist,2)+pow(z_dist,2));
	
	/* //calcul juste
	/* //1 = a*(rdr*WT) + b; */
	/* //0 = a*(rdr*WT + cancel_dist) + b; */
	/* //1 = a*(-cancel_dist); */
	/* //a = -1/cancel_dist; */
	/* //b = (x_WT + rdr*WT + cancel_dist)/cancel_dist; */
	//
	cs_real_t a_radius = -1.0/cancel_dist;
	cs_real_t b_radius = -a_radius*(WT_radius*rotor_rdr + cancel_dist);
	cs_real_t a_xdist = -1.0/cancel_dist;
	cs_real_t b_xdist = -a_xdist*(WT_radius*rotor_thickness_rdr + cancel_dist);

	st_coeff[c_id]=0.0;
	if (fabs(x_dist)<=rotor_thickness_rdr*WT_radius + select_dist)  {
	  st_coeff[c_id]=a_xdist*fabs(x_dist) + b_xdist;
	  if (cell_radius>rotor_rdr*WT_radius + select_dist) {
	    st_coeff[c_id]= 0.0;
	  }
	  else if (cell_radius>rotor_rdr*WT_radius) {
	    st_coeff[c_id]= a_radius*cell_radius + b_radius;
	  }
	}
	if (cell_radius<=rotor_rdr*WT_radius && fabs(x_dist)<=rotor_thickness_rdr*WT_radius) {
	  st_coeff[c_id]=1.0;
	}
	WT_volume=WT_volume+st_coeff[c_id]*cell_vol[c_id];
      }
	
      for (cs_lnum_t nc = 0; nc < z->n_elts; nc++) {
	cs_lnum_t c_id = z->elt_ids[nc];
	cs_real_t AD_coeff=0.5*cpro_rom[c_id]*WT_surf*Ct;
	// explicit source term
	_st_exp[c_id][0]    = -(cell_vol[c_id]*st_coeff[c_id])*AD_coeff*u_sonde*u_x/WT_volume;      
	_st_exp[c_id][1]    = -(cell_vol[c_id]*st_coeff[c_id])*AD_coeff*u_sonde*u_y/WT_volume;

	/* /\* // semi implicite source terms *\/ */
	/* _st_exp[c_id][0]    = -cell_vol[c_id]*AD_coeff*st_coeff[c_id]*(-1*u_x*u_x*u_x/u_sonde - u_x*u_y*u_y/u_sonde)/WT_volume; */
	/* _st_exp[c_id][1]    = -cell_vol[c_id]*AD_coeff*st_coeff[c_id]*(-1*u_y*u_y*u_y/u_sonde - u_y*u_x*u_x/u_sonde)/WT_volume; */
	
	/* _st_imp[c_id][0][0] = -cell_vol[c_id]*AD_coeff*st_coeff[c_id]*(u_sonde + u_x*u_x/u_sonde)/WT_volume; */
	/* _st_imp[c_id][0][1] = -cell_vol[c_id]*AD_coeff*st_coeff[c_id]*(u_x*u_y/u_sonde)/WT_volume; */
	
	/* _st_imp[c_id][1][0] = -cell_vol[c_id]*AD_coeff*st_coeff[c_id]*(u_y*u_x/u_sonde)/WT_volume; */
	/* _st_imp[c_id][1][1] = -cell_vol[c_id]*AD_coeff*st_coeff[c_id]*(u_sonde+u_y*u_y/u_sonde)/WT_volume; */
	//
	source_term_x[c_id] = -cell_vol[c_id]*AD_coeff*u_sonde*u_x*st_coeff[c_id]/WT_volume;
	source_term_y[c_id] = -cell_vol[c_id]*AD_coeff*u_sonde*u_y*st_coeff[c_id]/WT_volume;
	
      }
    }
    /************************************************************/

  }
  
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
