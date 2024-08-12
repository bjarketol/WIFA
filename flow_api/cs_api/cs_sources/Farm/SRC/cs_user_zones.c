/*============================================================================
 * Boundary and volume zone definitions.
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/*----------------------------------------------------------------------------
 * Local headers
 *----------------------------------------------------------------------------*/

#include "cs_headers.h"
#include "prototypes.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Define volume and surface zones.
 *
 * See \ref sec_selection_criteria for details on selection criteria.
 */
/*----------------------------------------------------------------------------*/

#pragma weak cs_user_zones
void
cs_user_zones(void)
{

  /******************************************
   * READ Wind Turbines from coordinates file
   *******************************************/

  //file name and tables
  char WT_file_name[1024];
  sprintf(WT_file_name,"turbines_info.csv");
  cs_real_t *WT_x_coords = NULL;
  cs_real_t *WT_y_coords = NULL;
  cs_real_t *WT_z_coords = NULL;
  cs_real_t *WT_diameters=NULL;
  cs_lnum_t *WT_types=NULL;

  //calculate number of lines in the file
  FILE* stream = fopen(WT_file_name, "r");
  char line[1024];
  char *tok0;
  const char sep0[4] = ";";
  size_t n_WT = 0;
  size_t i=-1;
  //read header
  //fgets(line, 1024, stream);
  while (fgets(line, 1024, stream) != NULL){
    i++;
    if (i>0) {
      // try to parse line
      tok0 = strtok(line, sep0);
      if (tok0 != NULL) {
        n_WT++;
      }
    }
  }
  fclose(stream);

  //allocate memory for all members of emission structure
  BFT_MALLOC(WT_x_coords, n_WT, cs_real_t);
  BFT_MALLOC(WT_y_coords, n_WT, cs_real_t);
  BFT_MALLOC(WT_z_coords, n_WT, cs_real_t);
  BFT_MALLOC(WT_diameters, n_WT, cs_real_t);
  BFT_MALLOC(WT_types, n_WT, cs_lnum_t);

  //read and parse this file
  stream = fopen(WT_file_name, "r");
  char *tok;
  const char sep[4] = ",";
  size_t line_count = -1;
  i=-1;

  while (fgets(line, 1024, stream) != NULL){
    i++;
    if (++line_count >= n_WT+1)
      break;

    //skip header
    if(i>0) {
      // parse line
      tok = strtok(line, sep);
      if (tok == NULL)
  continue;
      WT_x_coords[i-1]=atof(tok);

      tok = strtok(NULL, sep);
      if (tok == NULL)
  continue;
      WT_y_coords[i-1]=atof(tok);

      tok = strtok(NULL, sep);
      if (tok == NULL)
  continue;
      WT_z_coords[i-1]=atof(tok);

      tok = strtok(NULL, sep);
      if (tok == NULL)
        continue;
      WT_diameters[i-1]=atof(tok);

      tok = strtok(NULL, sep);
      if (tok == NULL)
        continue;
      WT_types[i-1]=atoi(tok);
    }
  }
  fclose(stream);

  //
  /******************************************
   * MARK TURBINE LOCATIONS USING CYLINDERS
   *******************************************/
  //
  //
  //
  char name[128];//for Actuator Disk (AD) zone name
  char criteria[100];
  //
  cs_real_t base_corner_coords[3];
  cs_real_t xTranslate_base_corner_coords[3];
  cs_real_t yTranslate_base_corner_coords[3];
  cs_real_t zTranslate_base_corner_coords[3];

  cs_lnum_t start_WT_count, end_WT_count;
  if (cs_notebook_parameter_value_by_name("isol")>0) {
    start_WT_count = cs_notebook_parameter_value_by_name("isol")-1;
    end_WT_count = cs_notebook_parameter_value_by_name("isol");
  }
  else {
    start_WT_count = 0;
    end_WT_count = n_WT;
  }
  for (cs_lnum_t WT_count=start_WT_count; WT_count < end_WT_count; WT_count ++){
    sprintf(name,"turbine_%d",WT_count+1);
    cs_real_t WT_radius = WT_diameters[WT_count]/2;

    /***********************************************************/
    /* Cylinder corresponding to AD zone */
    base_corner_coords[0] =  WT_x_coords[WT_count]-2*WT_radius;
    base_corner_coords[1] =  WT_y_coords[WT_count]-2*WT_radius;
    base_corner_coords[2] =  WT_z_coords[WT_count]-2*WT_radius;
    //dx
    xTranslate_base_corner_coords[0] = 4*WT_radius;
    xTranslate_base_corner_coords[1] = 0.0;
    xTranslate_base_corner_coords[2] = 0.0;
    //dy
    yTranslate_base_corner_coords[0] = 0.0;
    yTranslate_base_corner_coords[1] = 4*WT_radius;
    yTranslate_base_corner_coords[2] = 0.0;
    //dz
    zTranslate_base_corner_coords[0] = 0.0;
    zTranslate_base_corner_coords[1] = 0.0;
    zTranslate_base_corner_coords[2] = 4*WT_radius;

    //select zone
    sprintf(criteria, "box[%.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f]", \
    base_corner_coords[0],base_corner_coords[1],base_corner_coords[2],
    xTranslate_base_corner_coords[0],xTranslate_base_corner_coords[1],
    xTranslate_base_corner_coords[2],yTranslate_base_corner_coords[0],
    yTranslate_base_corner_coords[1],yTranslate_base_corner_coords[2],
    zTranslate_base_corner_coords[0],zTranslate_base_corner_coords[1],
    zTranslate_base_corner_coords[2]);
    cs_volume_zone_define(name, criteria, CS_VOLUME_ZONE_SOURCE_TERM);

    bft_printf("Creation of zone %s \n",name);
  }
  BFT_FREE(WT_x_coords);
  BFT_FREE(WT_y_coords);
  BFT_FREE(WT_z_coords);
  BFT_FREE(WT_diameters);
  BFT_FREE(WT_types);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
