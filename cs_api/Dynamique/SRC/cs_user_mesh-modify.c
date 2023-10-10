/*============================================================================
 * Definition of the calculation mesh.
 *
 * Mesh modification function examples.
 *============================================================================*/

/* VERS */

/*
  This file is part of code_saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2023 EDF S.A.

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
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"
#include "prototypes.h"
/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_mesh.c
 *
 * \brief Mesh modification example.
 *
 * See \ref cs_user_mesh for examples.
 */
/*----------------------------------------------------------------------------*/
/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Modify geometry and mesh.
 *
 * \param[in,out] mesh  pointer to a cs_mesh_t structure
*/
/*----------------------------------------------------------------------------*/

void
yawing(cs_real_t *coords, cs_real_t wind_dir, cs_real_t *mean_coords, cs_real_t *yawed_coords) {
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
cs_user_mesh_modify(cs_mesh_t  *mesh)
{
  /* Refine a selected portion of a mesh */

  /*! [mesh_modify_refine_1] */

  /******************************************
   * READ Wind Turbines from coordinates file
   *******************************************/
  
  cs_real_t WT_d = 136.0; //TODO: read as notebook param
  cs_real_t WT_radius = WT_d/2.; 
  cs_real_t wind_dir=cs_glob_atmo_option->meteo_angle;

  //file name and tables
  char WT_file_name[1024];
  sprintf(WT_file_name,"placement_turbines.csv");
  cs_real_t  *WT_x_coords = NULL;
  cs_real_t  *WT_y_coords = NULL;
  cs_real_t  *WT_z_coords = NULL;
  cs_real_t WT_min_x,WT_max_x,WT_min_y,WT_max_y,WT_min_z,WT_max_z;
  cs_real_t WF_length,WF_depth,WF_height,WF_size;
  cs_real_t WT_min_distance; //minimum distance between 2 turbines in the farm
  
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

    if(i==0) {
      WT_min_x=WT_x_coords[i];
      WT_max_x=WT_x_coords[i];
      WT_min_y=WT_y_coords[i];
      WT_max_y=WT_y_coords[i];      
      WT_min_z=WT_z_coords[i];
      WT_max_z=WT_z_coords[i];      
    }
    else {
      if(WT_x_coords[i]<WT_min_x) {
	WT_min_x=WT_x_coords[i];
      }
      if(WT_x_coords[i]>WT_max_x) {
	WT_max_x=WT_x_coords[i];
      }
      if(WT_y_coords[i]<WT_min_y) {
	WT_min_y=WT_y_coords[i];
      }
      if(WT_y_coords[i]>WT_max_y) {
	WT_max_y=WT_y_coords[i];
      }	
      if(WT_z_coords[i]<WT_min_z) {
	WT_min_z=WT_z_coords[i];
      }
      if(WT_z_coords[i]>WT_max_z) {
	WT_max_z=WT_z_coords[i];
      }	
    }
  }
  fclose(stream);
  //
  WF_length=WT_max_x-WT_min_x + WT_d;
  WF_depth=WT_max_y-WT_min_y + WT_d;
  WF_size=sqrt(pow(WF_length,2.0)+pow(WF_depth,2.0));
  
  /******************************************
   * MARK TURBINE LOCATIONS USING CYLINDERS
   *******************************************/ 
  cs_real_t mesh_coords[3];
  cs_real_t yawed_mesh_coords[3];
  cs_real_t mean_coords[3];
  cs_lnum_t  vtx_id;
  //mean=0 because rotation in global referential
  mean_coords[0] = 0.0;
  mean_coords[1] = 0.0;
  mean_coords[2] = 0.0;
  //
  char criteria[100];
  //Refinement Distance Ratio
  //relative to Wind Farm lenght, or Wind Turbine diameter
  cs_real_t rdr;
  //
  cs_lnum_t   n_selected_cells = 0;
  cs_lnum_t  *selected_cells = NULL;
  //
  cs_real_t base_corner_coords[3];
  cs_real_t xTranslate_base_corner_coords[3];
  cs_real_t yTranslate_base_corner_coords[3];
  cs_real_t zTranslate_base_corner_coords[3];
  //
  cs_real_t yawed_base_corner_coords[3];
  cs_real_t twice_yawed_base_corner_coords[3];
  cs_real_t yawed_xTranslate_base_corner_coords[3];
  cs_real_t yawed_yTranslate_base_corner_coords[3];
  cs_real_t yawed_zTranslate_base_corner_coords[3];
  //
  cs_real_t disk_1_coords[3];
  cs_real_t disk_2_coords[3];
  cs_real_t yawed_disk_1_coords[3];
  cs_real_t yawed_disk_2_coords[3];
  //
  cs_real_t thickness_rdr; //ratio to WT diameter
  
  rdr=1.2; //ratio for disk radius relative to turbine radius
  thickness_rdr=0.2; //ratio for disk thickness relative to turbine diameter 
  
  char name[128];//for Actuator Disk (AD) zone name
  
  for (cs_lnum_t WT_count=0; WT_count < n_WT; WT_count ++){    
    /***********************************************************/
    /* Cylinder corresponding to AD zone */
    disk_1_coords[0] = WT_x_coords[WT_count] - 0.5 * thickness_rdr*WT_d;
    disk_1_coords[1] = WT_y_coords[WT_count];
    disk_1_coords[2] = WT_z_coords[WT_count];
    //
    disk_2_coords[0] = WT_x_coords[WT_count] + 0.5 * thickness_rdr*WT_d;
    disk_2_coords[1] = WT_y_coords[WT_count];
    disk_2_coords[2] = WT_z_coords[WT_count];
    
    //yaw disks centers in the cylinder local referential 
    mean_coords[0] = 0.5*(disk_1_coords[0]+disk_2_coords[0]);
    mean_coords[1] = 0.5*(disk_1_coords[1]+disk_2_coords[1]);
    mean_coords[2] = 0.5*(disk_1_coords[2]+disk_2_coords[2]);
    yawing(disk_1_coords, wind_dir, mean_coords, yawed_disk_1_coords);
    yawing(disk_2_coords, wind_dir, mean_coords, yawed_disk_2_coords);

    //select zone
    sprintf(criteria, "cylinder[%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f]", \
    	    yawed_disk_1_coords[0],yawed_disk_1_coords[1],yawed_disk_1_coords[2],
    	    yawed_disk_2_coords[0],yawed_disk_2_coords[1],yawed_disk_2_coords[2],
    	    rdr*WT_radius);
    
    bft_printf("Turbine_%d criteria 3 : %s\n",WT_count+1,criteria);
    
    BFT_MALLOC(selected_cells, mesh->n_cells, cs_lnum_t);

    cs_selector_get_cell_list(criteria,
			      &n_selected_cells,
			      selected_cells);

    //tag AD zone as volume zone
    sprintf(name,"turbine_%d",WT_count+1);
    cs_volume_zone_define(name, criteria, CS_VOLUME_ZONE_SOURCE_TERM);    /* Cylinder corresponding to AD zone */
    
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
