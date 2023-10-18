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
#include "prototypes.h"

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
yawing2(cs_real_t *coords, cs_real_t wind_dir, cs_real_t *mean_coords, cs_real_t *yawed_coords) {
  //
  const cs_real_t dpi = atan(1.0)*4.0;
  cs_real_t ra_wind_dir= (270.0-wind_dir)*dpi/180.0; //radian absolute wind_dir  
  //
  yawed_coords[0] = cos(ra_wind_dir)*(coords[0]-mean_coords[0]) - sin(ra_wind_dir)*(coords[1]-mean_coords[1]) + mean_coords[0];
  //
  yawed_coords[1] = sin(ra_wind_dir)*(coords[0]-mean_coords[0]) + cos(ra_wind_dir)*(coords[1]-mean_coords[1]) + mean_coords[1];
  //
  yawed_coords[2] = coords[2];
}

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
  cs_real_t WT_d = cs_notebook_parameter_value_by_name("WT_d");
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
  cs_real_t *power_values=NULL;
  BFT_MALLOC(WT_x_coords, n_WT, cs_real_t);
  BFT_MALLOC(WT_y_coords, n_WT, cs_real_t);
  BFT_MALLOC(WT_z_coords, n_WT, cs_real_t);
  BFT_MALLOC(power_values, n_WT, cs_real_t);

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
   * READ Cp and Ct from data files
   *******************************************/  
  //file name and tables
  char cp_file_name[1024];
  char ct_file_name[1024];
  sprintf(ct_file_name,"ct_table.csv");
  sprintf(cp_file_name,"cp_table.csv");
  //calculate number of lines in the file
  FILE* stream_ct = fopen(ct_file_name, "r");
  char *tok2;
  const char sep2[4] = ",";
  size_t n_ct = 0;
  //read header
  //fgets(line, 1024, stream);
  while (fgets(line, 1024, stream) != NULL){
    // try to parse line
    tok2 = strtok(line, sep2);
    if (tok2 != NULL)
      n_ct++;
  }
  fclose(stream_ct);

  //allocate memory for all members of emission structure
  cs_real_t *ct_speeds=NULL;
  cs_real_t *ct_values=NULL;
  BFT_MALLOC(ct_speeds, n_ct, cs_real_t);
  BFT_MALLOC(ct_values, n_ct, cs_real_t);

   //read and parse this file
  stream_ct = fopen(ct_file_name, "r");
  char *tok3;
  const char sep3[4] = ",";
  line_count = -1;
  i=-1;

  while (fgets(line, 1024, stream) != NULL){
    i++;
    if (++line_count >= n_ct)
      break;

    // parse line
    tok3 = strtok(line, sep3);
    if (tok3 == NULL)
      continue;
    ct_speeds[i]=atof(tok3)*pow(1.260/cs_glob_fluid_properties->ro0,1./3.);
    
    tok3 = strtok(NULL, sep3);
    if (tok3 == NULL)
      continue;
    ct_values[i]=atof(tok3);
  }
  fclose(stream_ct);

  //cp
  //calculate number of lines in the file
  FILE* stream_cp = fopen(cp_file_name, "r");
  char *tok4;
  const char sep4[4] = ",";
  size_t n_cp = 0;
  //read header
  //fgets(line, 1024, stream);
  while (fgets(line, 1024, stream) != NULL){
    // try to parse line
    tok4 = strtok(line, sep4);
    if (tok4 != NULL)
      n_cp++;
  }
  fclose(stream_cp);

  //allocate memory for all members of emission structure
  cs_real_t *cp_speeds=NULL;
  cs_real_t *cp_values=NULL;
  BFT_MALLOC(cp_speeds, n_cp, cs_real_t);
  BFT_MALLOC(cp_values, n_cp, cs_real_t);

   //read and parse this file
  stream_cp = fopen(cp_file_name, "r");
  char *tok5;
  const char sep5[4] = ",";
  line_count = -1;
  i=-1;

  while (fgets(line, 1024, stream) != NULL){
    i++;
    if (++line_count >= n_cp)
      break;

    // parse line
    tok5 = strtok(line, sep5);
    if (tok5 == NULL)
      continue;
    cp_speeds[i]=atof(tok5)*pow(1.260/cs_glob_fluid_properties->ro0,1./3.);
    
    tok5 = strtok(NULL, sep5);
    if (tok5 == NULL)
      continue;
    cp_values[i]=atof(tok5);
  }
  fclose(stream_cp);
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
    //
    cs_real_t WT_volume[n_WT];
    cs_real_t WT_ux[n_WT];
    cs_real_t WT_uy[n_WT];
    cs_real_t WT_uz[n_WT];
    cs_real_t WT_u[n_WT];
    cs_real_t u_sonde,u_x,u_y;
    cs_real_t WT_ct[n_WT];
    cs_real_t WT_ctstar[n_WT];
    cs_real_t WT_cp[n_WT];
    cs_real_t WT_cpstar[n_WT];
    //
    cs_real_t ra_wind_dir= (270.0-wind_dir)*dpi/180.0; //radian absolute wind_dir
    //
    cs_real_3_t *_st_exp = (cs_real_3_t *)st_exp;
    cs_real_33_t *_st_imp = (cs_real_33_t *)st_imp;
    
    cs_real_t rotor_rdr=1.0; //ratio for rotor radius relative to turbine radius
    cs_real_t rotor_thickness_rdr=0.05; //ratio for rotor thickness relative to turbine diameter
    cs_real_t AD_mesh_cell_size = 5.0;
    cs_real_t cancel_dist = 0.5*AD_mesh_cell_size;
    cs_real_t select_dist = AD_mesh_cell_size;
    /******************LOOP ON TURBINES*************************/
    for (cs_lnum_t WT_count=0; WT_count < n_WT; WT_count ++){
      //get Actuator Disk zone by name
      sprintf(name,"turbine_%d",WT_count+1);
      cs_zone_t *z = cs_volume_zone_by_name_try(name);
      //
      WT_volume[WT_count]=0.0;      
      for (cs_lnum_t nc = 0; nc < z->n_elts; nc++) {
	cs_lnum_t c_id = z->elt_ids[nc];
	//
	cs_real_t x_dist,y_dist,z_dist;
	x_dist = cell_cen[c_id][0] - WT_x_coords[WT_count];
	y_dist = cell_cen[c_id][1] - WT_y_coords[WT_count];
	z_dist = cell_cen[c_id][2] - WT_z_coords[WT_count];

	cs_real_t dxdydz[3];
	cs_real_t deyawed_dxdydz[3];
	cs_real_t mean_coords[3];
	dxdydz[0] = x_dist;
	dxdydz[1] = y_dist;
	dxdydz[2] = z_dist;
	mean_coords[0] = 0.0;
	mean_coords[1] = 0.0;
	mean_coords[2] = 0.0;	
	yawing2(dxdydz, 360.0 - wind_dir, mean_coords, deyawed_dxdydz);
	cs_real_t deyawed_x_dist,deyawed_y_dist,deyawed_z_dist;
	deyawed_x_dist = deyawed_dxdydz[0];
	deyawed_y_dist = deyawed_dxdydz[1];
	deyawed_z_dist = deyawed_dxdydz[2];
	
	cs_real_t cell_radius = sqrt(pow(deyawed_y_dist,2)+pow(deyawed_z_dist,2));
	//
	
	/* //calculate linear decay until 0 at cancel_dist
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
	if (cs_notebook_parameter_value_by_name("st_method")==0) {
	  if (fabs(deyawed_x_dist)<=rotor_thickness_rdr*WT_radius + select_dist)  {
	    st_coeff[c_id]=fmax(a_xdist*fabs(deyawed_x_dist) + b_xdist,0.0);
	    if (cell_radius>rotor_rdr*WT_radius + select_dist) {
	      st_coeff[c_id]= 0.0;
	    }
	    else if (cell_radius>rotor_rdr*WT_radius) {
	      st_coeff[c_id]= fmax(a_radius*cell_radius + b_radius,0.0);
	    }
	  }
	  if (cell_radius<=rotor_rdr*WT_radius && fabs(deyawed_x_dist)<=rotor_thickness_rdr*WT_radius) {
	    st_coeff[c_id]=1.0;
	  }
	}
	//Ren et al 2019
	else if (cs_notebook_parameter_value_by_name("st_method")==1) {
	  cs_real_t D0 = -4.0/3.0;
	  cs_real_t D1 = -2.0/9.0;
	  cs_real_t b=cell_radius/WT_radius;
	    
	  if (cell_radius<=rotor_rdr*WT_radius &&
	      fabs(deyawed_x_dist)<=rotor_thickness_rdr*WT_radius) {
	    st_coeff[c_id]= (D0+D1*pow(b,2)+pow(b,4))/(D0+0.5*D1+1./3.);
	  }
	  else {
	    st_coeff[c_id]= 0.0;
	  }
	}
	WT_volume[WT_count]=WT_volume[WT_count]+st_coeff[c_id]*cell_vol[c_id];
      }
      cs_parall_sum(1, CS_REAL_TYPE, &WT_volume[WT_count]);
      bft_printf(" Volume of WT %d is %.2f \n ",WT_count, WT_volume[WT_count]);
	      
      WT_u[WT_count]=0.0;
      WT_ux[WT_count]=0.0;
      WT_uy[WT_count]=0.0;
      WT_uz[WT_count]=0.0;
      for (cs_lnum_t nc =0; nc < z->n_elts; nc++){
      	cs_lnum_t c_id = z->elt_ids[nc];
      	WT_ux[WT_count] = WT_ux[WT_count] + vel[c_id][0]*st_coeff[c_id]*cell_vol[c_id]/WT_volume[WT_count];
      	WT_uy[WT_count] = WT_uy[WT_count] + vel[c_id][1]*st_coeff[c_id]*cell_vol[c_id]/WT_volume[WT_count];
      	WT_uz[WT_count] = WT_uz[WT_count] + vel[c_id][2]*st_coeff[c_id]*cell_vol[c_id]/WT_volume[WT_count];
      }
      cs_parall_sum(1, CS_REAL_TYPE, &WT_ux[WT_count]);
      cs_parall_sum(1, CS_REAL_TYPE, &WT_uy[WT_count]);
      cs_parall_sum(1, CS_REAL_TYPE, &WT_uz[WT_count]);      
      WT_u[WT_count]=sqrt(pow(WT_ux[WT_count],2)+pow(WT_uy[WT_count],2)+pow(WT_uz[WT_count],2));
      //bft_printf(" Velocities of WT %d are (ux,uy,umagn)=(%.4f,%.4f,%.4f) \n ",WT_count, WT_ux[WT_count],WT_uy[WT_count],WT_u[WT_count]);

      ////for testing
      //u_sonde=5.0;
      //u_x = u_sonde*cos(ra_wind_dir);
      //u_y = u_sonde*sin(ra_wind_dir);
      //
      cs_real_t ct_1, ct_2, ct_u_1, ct_u_2;
      ct_u_1 = ct_speeds[0];
      ct_u_2 = ct_speeds[n_ct-1];
      if(WT_u[WT_count]<ct_u_1 || WT_u[WT_count]>ct_u_2) {
	WT_ct[WT_count] = 0.0;
	WT_ctstar[WT_count] = 0.0;
      }
      else {
	ct_u_2 = ct_speeds[0];
	ct_2 = ct_values[0];
	for (cs_lnum_t ct_count=1; ct_count < n_ct; ct_count ++){
	  ct_u_1=ct_u_2;
	  ct_u_2=ct_speeds[ct_count];
	  ct_1=ct_2;
	  ct_2=ct_values[ct_count];
	  if(WT_u[WT_count]>=ct_u_1 && WT_u[WT_count]<=ct_u_2) {
	    /* //linear interpolation of ct between two values;
	       /* //ct_1 = a*ct_u_1 + b; */
	    /* //ct_2 = a*ct_u_2 + b; */
	    /* //a = (ct_2 - ct_1)/(ct_u_2-ct_u_1); */
	    /* //b = ct_1 - a*ct_u_1; */
	    //
	    cs_real_t a_ct = (ct_2 - ct_1)/(ct_u_2-ct_u_1);
	    cs_real_t b_ct = ct_1 - a_ct*ct_u_1;
	    WT_ct[WT_count] = a_ct*WT_u[WT_count]+b_ct;
	    WT_ctstar[WT_count] = 4.0*WT_ct[WT_count]/pow(1+sqrt(1-WT_ct[WT_count]),2);
	    break;
	  }
	}
      }
      cs_real_t cp_1, cp_2, cp_u_1, cp_u_2;
      cp_u_1 = cp_speeds[0];
      cp_u_2 = cp_speeds[n_cp-1];
      if(WT_u[WT_count]<cp_u_1 || WT_u[WT_count]>cp_u_2) {
	WT_cp[WT_count] = 0.0;
	WT_cpstar[WT_count] = 0.0;
      }
      else {
	cp_u_2 = cp_speeds[0];
	cp_2 = cp_values[0];
	for (cs_lnum_t cp_count=1; cp_count < n_cp; cp_count ++){
	  cp_u_1=cp_u_2;
	  cp_u_2=cp_speeds[cp_count];
	  cp_1=cp_2;
	  cp_2=cp_values[cp_count];
	  if(WT_u[WT_count]>=cp_u_1 && WT_u[WT_count]<=cp_u_2) {
	    /* //linear interpolation of cp between two values;
	       /* //cp_1 = a*cp_u_1 + b; */
	    /* //cp_2 = a*cp_u_2 + b; */
	    /* //a = (cp_2 - cp_1)/(cp_u_2-cp_u_1); */
	    /* //b = cp_1 - a*cp_u_1; */
	    //
	    cs_real_t a_cp = (cp_2 - cp_1)/(cp_u_2-cp_u_1);
	    cs_real_t b_cp = cp_1 - a_cp*cp_u_1;
	    WT_cp[WT_count] = a_cp*WT_u[WT_count]+b_cp;
	    WT_cpstar[WT_count] = 8.0*WT_cp[WT_count]/pow(1+sqrt(1-WT_ct[WT_count]),3);
	    break;
	  }
	}
      }

      if (cs_notebook_parameter_value_by_name("isol")>0) {
	//isolation
	if(WT_count != cs_notebook_parameter_value_by_name("isol")-1) {
	  WT_cpstar[WT_count] = 0.0;
	  WT_ctstar[WT_count] = 0.0;
	  WT_cp[WT_count] = 0.0;
	  WT_ct[WT_count] = 0.0;
	  WT_ux[WT_count] = 0.0;
	  WT_uy[WT_count] = 0.0;
	  WT_u[WT_count] = 0.0;
	}
      }
	
      //bft_printf(" Velocity, CT and CP of WT %d are %.4f, %.4f and %.4f\n ",WT_count, WT_u[WT_count], WT_ct[WT_count], WT_cp[WT_count]);
      //
      //WT_cpstar as u is taken at wind turbine position, not freestream
      power_values[WT_count] = 0.5*cs_glob_fluid_properties->ro0*WT_surf*WT_cpstar[WT_count]*pow(WT_u[WT_count],3);
      for (cs_lnum_t nc = 0; nc < z->n_elts; nc++) {
	cs_lnum_t c_id = z->elt_ids[nc];
	//WT_ctstar as u is taken at wind turbine position, not freestream
	cs_real_t AD_coeff=0.5*cpro_rom[c_id]*WT_surf*WT_ctstar[WT_count];
	// explicit source term
	_st_exp[c_id][0]    = -(cell_vol[c_id]*st_coeff[c_id])*AD_coeff*WT_u[WT_count]*WT_ux[WT_count]/WT_volume[WT_count];      
	_st_exp[c_id][1]    = -(cell_vol[c_id]*st_coeff[c_id])*AD_coeff*WT_u[WT_count]*WT_uy[WT_count]/WT_volume[WT_count];
	source_term_x[c_id] = _st_exp[c_id][0];
	source_term_y[c_id] = _st_exp[c_id][1];
	/* /\* // semi implicite source terms *\/ */
	/* _st_exp[c_id][0]    = -cell_vol[c_id]*AD_coeff*st_coeff[c_id]*(-1*u_x*u_x*u_x/u_sonde - u_x*u_y*u_y/u_sonde)/WT_volume[WT_count]; */
	/* _st_exp[c_id][1]    = -cell_vol[c_id]*AD_coeff*st_coeff[c_id]*(-1*u_y*u_y*u_y/u_sonde - u_y*u_x*u_x/u_sonde)/WT_volume[WT_count]; */
	
	/* _st_imp[c_id][0][0] = -cell_vol[c_id]*AD_coeff*st_coeff[c_id]*(u_sonde + u_x*u_x/u_sonde)/WT_volume[WT_count]; */
	/* _st_imp[c_id][0][1] = -cell_vol[c_id]*AD_coeff*st_coeff[c_id]*(u_x*u_y/u_sonde)/WT_volume[WT_count]; */
	
	/* _st_imp[c_id][1][0] = -cell_vol[c_id]*AD_coeff*st_coeff[c_id]*(u_y*u_x/u_sonde)/WT_volume[WT_count]; */
	/* _st_imp[c_id][1][1] = -cell_vol[c_id]*AD_coeff*st_coeff[c_id]*(u_sonde+u_y*u_y/u_sonde)/WT_volume[WT_count]; */
	/* source_term_x[c_id] = -(cell_vol[c_id]*st_coeff[c_id])*AD_coeff*u_sonde*u_x/WT_volume[WT_count]; */
	/* source_term_y[c_id] = -(cell_vol[c_id]*st_coeff[c_id])*AD_coeff*u_sonde*u_y/WT_volume[WT_count]; */
	
      }      
    }
    
    cs_real_t total_power;
    total_power=0.0;      
    for (cs_lnum_t WT_count=0; WT_count < n_WT; WT_count ++){
      total_power=total_power + power_values[WT_count];
    }
    if ((ts->nt_cur) == (ts->nt_max)-1){
      FILE* stream10 = fopen("power.txt", "w");
      fprintf(stream10,"#total wind farm power is %.4f\n",total_power);
      fprintf(stream10,"#wind_turbine, velocity, ct, ct*, cp, cp*, power\n");
      for (cs_lnum_t WT_count=0; WT_count < n_WT; WT_count ++){
	fprintf(stream10,"%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f\n",WT_count+1, WT_u[WT_count], WT_ct[WT_count],WT_ctstar[WT_count], WT_cp[WT_count], WT_cpstar[WT_count],power_values[WT_count]);
      }
      fclose(stream10);
    }
    /************************************************************/

  }
  
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
