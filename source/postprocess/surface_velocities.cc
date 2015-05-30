/*
  Copyright (C) 2011 - 2015 by the authors of the ASPECT code.

  This file is part of ASPECT.

  ASPECT is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2, or (at your option)
  any later version.

  ASPECT is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with ASPECT; see the file doc/COPYING.  If not see
  <http://www.gnu.org/licenses/>.
 */
/*  $Id: surface_velocities.cc $  */

//#include <aspect/postprocess/surface_velocities.h>
#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>
#include <deal.II/grid/grid_tools.h>
#include <boost/iostreams/filtering_stream.hpp>

#include <aspect/global.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>

// Figure out if these are needed
#include <aspect/geometry_model/interface.h>
#include <aspect/initial_conditions/interface.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>
#include <aspect/material_model/interface.h>
#include <fstream>
#include <iostream>
#include <string>



namespace aspect
{
  namespace Postprocess
  {
	template <int dim>
	class SurfaceVelocities : public Interface<dim>, public ::aspect::SimulatorAccess<dim>
	{
		public:
			/**
			 * Generate graphical output from the current solution.
			 */
			virtual
			std::pair<std::string,std::string>
			execute (TableHandler &statistics);
	};
  }
}

namespace aspect
{
	namespace Postprocess
	{
		template <int dim>
		std::pair<std::string,std::string>
		SurfaceVelocities<dim>::execute (TableHandler &statistics)
		{
			const QMidpoint<dim-1> quadrature_formula;

			FEFaceValues<dim> fe_face_values (this->get_mapping(),
					this->get_fe(),
					quadrature_formula,
					update_values |
					update_normal_vectors |
					update_q_points);

			std::vector<Tensor<1,dim> > velocity_values (quadrature_formula.size());

			std::string filename = this->get_output_directory() +
					"surface_velocities." +
					Utilities::int_to_string
					(this->get_triangulation().locally_owned_subdomain(), 4);
			std::ofstream output (filename.c_str());

			typename DoFHandler<dim>::active_cell_iterator
			cell = this->get_dof_handler().begin_active(),
			endc = this->get_dof_handler().end();

	        unsigned int cell_index = 0;
	        for (; cell!=endc; ++cell,++cell_index)
				if (cell->is_locally_owned())
		            if (cell->at_boundary())
		            {
		            	// check if cell is at the top boundary
		                unsigned int top_face_idx = numbers::invalid_unsigned_int;
		                {
		                	for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
		                        if (cell->at_boundary(f) && this->get_geometry_model().depth (cell->face(f)->center()) < cell->face(f)->minimum_vertex_distance()/3)
		                        {
		                            top_face_idx = f;
		                            break;
		                        }
		                }
		                if (top_face_idx == numbers::invalid_unsigned_int)
		                {
		                    (*return_value.second)(cell_index) = 0;
		                    continue;
		                }

		            }
	        			// extract velocity solution for cells at the top boundary
						{
							fe_face_values.reinit (cell, f);
							fe_face_values[this->introspection().extractors.velocities].get_function_values (this->get_solution(),
									velocity_values);

							for (unsigned int q=0; q<quadrature_formula.size(); ++q)
							{
								Point<3> p (fe_face_values.quadrature_point(q)[0],
										fe_face_values.quadrature_point(q)[1],
										fe_face_values.quadrature_point(q)[2]);

								/**
								 * Compute the latitude and longitude that corresponds to a point on a
								 * sphere.
								 *
								 * @param p The given point for which we'd like to know lat/long
								 * @return A pair of values where the first part corresponds to the latitude
								 *   and the second to the longitude. These values are returned in degrees,
								 *   where latitude is counted north from the equator and longitude is counted
								 *   east from the Greenwich meridian.
								 */

								std::pair<double,double> lat_long = make_pair(std::atan2(p[2],
										std::sqrt(p[0] * p[0] + p[1] * p[1]))
										* 180
										/ numbers::PI,
										std::atan2(p[1], p[0]) * 180 / numbers::PI);

								// output velocities in file with lat lon Vx Vy Vz in m/yr
								output << lat_long.second << ' ' << lat_long.first
										<< ' ' << (velocity_values[q] * year_in_seconds) / 1000
										<< std::endl;
							}
						}

			return std::pair<std::string,std::string>("Writing surface velocities...",
					filename);
		}
	}
}


// explicit instantiations
namespace aspect
{
	namespace Postprocess
	{
	ASPECT_REGISTER_POSTPROCESSOR(SurfaceVelocities,
			"surface velocities",
			"A postprocessor that outputs surface velocities at the centers "
			"of surface faces for a box.cc model with the surface assigned face 5.")
	}
}
