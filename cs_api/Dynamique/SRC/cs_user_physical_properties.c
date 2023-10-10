/*============================================================================
 * User definition of physical properties.
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

/*----------------------------------------------------------------------------*/
/*!
 * \file cs_user_physical_properties.c
 *
 * \brief User definition of physical properties.
 */
/*----------------------------------------------------------------------------*/

/*============================================================================
 * User function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Function called at each time step to define physical properties.
 *
 * \param[in, out]  domain   pointer to a cs_domain_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_user_physical_properties(cs_domain_t   *domain)
{
  /*physical variable */
  cs_real_t *actuator_disk = cs_field_by_name("actuator_disk")->val;
  cs_real_t *is_probe = cs_field_by_name("is_probe")->val;
  cs_lnum_t   cell_id = 0;
  double cell_vol = 0.;
  double domain_volume = 0.;


  /******************************************
   * READ Wind Turbines from coordinates file
   *******************************************/
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
  
  /******************************************
   * HIGHLIGHT ACTUATOR DISK ZONES
   *******************************************/
  char name[128];
   
  for (cs_lnum_t WT_count=0; WT_count < n_WT; WT_count ++){
    //get zone by name
    sprintf(name,"turbine_%d",WT_count+1);
    cs_zone_t *z = cs_volume_zone_by_name_try(name);
    
    //set "actuator_disk" marker to 1
    for (cs_lnum_t nc =0; nc < z->n_elts; nc++){
      cell_id = z->elt_ids[nc];
      cell_vol = cs_glob_mesh_quantities->cell_vol[cell_id];
      domain_volume = domain_volume + cell_vol;
      actuator_disk[cell_id]=1.0;
    }
  }
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
