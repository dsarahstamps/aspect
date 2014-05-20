/*
 Copyright (C) 2011, 2012 by the authors of the ASPECT code.

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
/*  $Id: africa.cc 1538 2013-01-06 03:12:23Z bangerth $  */

#include <aspect/geometry_model/interface.h>
#include <aspect/initial_conditions/interface.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/tria_accessor.h>
#include <aspect/postprocess/interface.h>
#include <aspect/simulator_access.h>
#include <aspect/material_model/interface.h>
#include <aspect/simulator_access.h>
#include <fstream>
#include <iostream>


namespace
{
using namespace dealii;

/**
 * Compute the Cartesian coordinate corresponding to a point at given
 * latitude/longitude and a given height relative to the WGS84 reference
 * ellipsoid.
 *
 * @param latitude The latitude of the point, in degrees counted
 *   north from the equator.
 * @param longitude The latitude of the point, in degrees counted
 *   east from Greenwich meridian.
 * @param height The height (or, if negative, the depth) of a point
 *   relative to the WGS84 reference surface.
 * @return
 */
Point<3>
xyz_from_lat_long_height_in_wgs84(const double latitude,
		const double longitude,
		const double height)
		{
	/* WGS84 ellipsoid constants */
	const double radius = 6378137;
	const double ellipticity = 8.1819190842622e-2;

	/* convert to radians */
	const double latitude_rad = latitude * numbers::PI / 180.0;
	const double longitude_rad = longitude * numbers::PI / 180.0;

	/* intermediate calculation for vertical radius of curvature */
	const double N = radius
			/ std::sqrt(1
					- (ellipticity * ellipticity * std::sin(latitude_rad)
	* std::sin(latitude_rad)));

	/* results */
	const Point<3> p((N + height) * std::cos(latitude_rad)
	* std::cos(longitude_rad),
	(N + height) * std::cos(latitude_rad)
	* std::sin(longitude_rad),
	((1 - ellipticity * ellipticity) * N + height) * std::sin(latitude_rad));

	return p;
		}

/**
 * Like xyz_from_lat_long_height_wgs84, but with respect to a sphere that has the
 * same radius as the earth at the equator.
 */
Point<3>
xyz_from_lat_long_height(const double latitude,
		const double longitude,
		const double height)
		{
	const double radius = 6378137;

	/* convert to radians */
	const double latitude_rad = latitude * numbers::PI / 180.0;
	const double longitude_rad = longitude * numbers::PI / 180.0;

	return Point<3>((radius + height) * std::cos(latitude_rad)
	* std::cos(longitude_rad),
	(radius + height) * std::cos(latitude_rad)
	* std::sin(longitude_rad),
	(radius + height) * std::sin(latitude_rad));
		}

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
std::pair<double, double>
lat_long_from_xyz(const Point<3> &p)
{
	return std::make_pair(std::atan2(p[2],
			std::sqrt(p[0] * p[0] + p[1] * p[1]))
	* 180
	/ numbers::PI,
	std::atan2(p[1], p[0]) * 180 / numbers::PI);
}

/**
 * Compute the latitude and longitude that corresponds to a point on the
 * WGS84 surface.
 *
 * @param p The given point for which we'd like to know lat/long
 * @return A pair of values where the first part corresponds to the latitude
 *   and the second to the longitude. These values are returned in degrees,
 *   where latitude is counted north from the equator and longitude is counted
 *   east from the Greenwich meridian.
 */
std::pair<double, double>
lat_long_from_xyz_wgs84(const Point<3> &pos)
{
	/* WGS84 ellipsoid constants */
	const double radius = 6378137;
	const double ellipticity = 8.1819190842622e-2;

	/* calculations */
	const double b = std::sqrt(radius * radius
			* (1 - ellipticity * ellipticity));
	const double ep = std::sqrt((radius * radius - b * b) / (b * b));
	const double p = std::sqrt(pos(0) * pos(0) + pos(1) * pos(1));
	const double th = std::atan2(radius * pos(2), b * p);
	const double lon = std::atan2(pos(1), pos(0));
	const double lat = std::atan2((pos(2) + ep * ep * b * std::sin(th)
	* std::sin(th) * std::sin(th)),
			(p - (ellipticity
					* ellipticity
					* radius
					* (std::cos(th) * std::cos(th)
	* std::cos(th)))));
	const double N = radius
			/ (std::sqrt(1
					- ellipticity * ellipticity * std::sin(lat) * std::sin(lat)));
	const double alt = p / std::cos(lat) - N;

	/* convert to degrees */
	const double lon_degrees = lon * (180 / numbers::PI);
	const double lat_degrees = lat * (180 / numbers::PI);

	/* set all longitudes between [0,360] */
	if (lon_degrees < 0)
		return std::make_pair(lat_degrees, lon_degrees + 360);
	else if (lon_degrees > 360)
		return std::make_pair(lat_degrees, lon_degrees - 360);
	else
		return std::make_pair(lat_degrees, lon_degrees);
}
}

namespace aspect
{
namespace GeometryModel
{
using namespace dealii;

/**
 * A namespace in which we define the parameters that describe the
 * domain.
 */
namespace DomainData
{
double latitudes[6];
double longitudes[6];
double bottom_depth;
}

/**
 * A class that reads the topography in from a file and provides
 * a function that allows querying it at arbitrary latitudes/longitudes.
 */
class Topography
{
public:
	/**
	 * Constructor
	 *
	 * @param filename Name of a file from which to read the topography.
	 */
	double delta;
	int topo_flag;
	int number_coords;

	Topography(const std::string &filename)
	{

		std::ifstream input2(filename.c_str());
		if (input2.is_open())
		{
			std::string line;

			int count = 0;
			while (true)
			{
				double latitude_alt, longitude_alt, altitude;
				input2 >> latitude_alt >> longitude_alt >> altitude;
				if (input2.eof())
					break;

				latitudes_alt.push_back(latitude_alt);
				longitudes_alt.push_back(longitude_alt);
				altitudes.push_back(altitude);

				/** Find first 2 numbers that are different to use in
				 * calculating half the difference between each position as delta.
				 */
				if (count < 2)
				{
					count ++;
				}
				if (count == 2)
				{
					//	std::cout << "latitudes_alt[0]: "<< latitudes_alt[0] << std::endl;
					//	std::cout << "latitudes_alt[1]: "<< latitudes_alt[1] << std::endl;
					//	std::cout << "longitudes_alt[0]: "<< longitudes_alt[0] << std::endl;
					//	std::cout << "longitudes_alt[1]: "<< longitudes_alt[1] << std::endl;
					//	std::cout << "latitude_alt: "<< latitude_alt << std::endl;
					//	std::cout << "longitude_alt: "<< longitude_alt << std::endl;

					if ((std::fabs(latitudes_alt[0] - latitudes_alt[1]) > 1e-9) && (std::fabs(longitudes_alt[0] - longitudes_alt[1]) > 1e-9))
					{
						// If file not formatted correctly stop program
						std::cout << ""<< std::endl;
						throw std::ios_base::failure("Topography file not formatted correctly. " + filename + "Either longitude or latitudes must be increasing.");
						std::cout << ""<< std::endl;
					}
					if ((std::fabs(latitudes_alt[0] - latitudes_alt[1]) < 1e-9) && (std::fabs(longitudes_alt[0] - longitudes_alt[1]) < 1e-9))
					{
						// If file not formatted correctly stop program
						std::cout << ""<< std::endl;
						throw std::ios_base::failure("Topography file not formatted correctly. " + filename + "Either longitude or latitudes must be increasing.");
						std::cout << ""<< std::endl;
					}

					if (std::fabs(latitudes_alt[0] - latitudes_alt[1]) > 1e-9)
					{
						// Calculate half the distance between points for delta
						delta = (std::fabs((0.5)*(latitudes_alt[0] - latitudes_alt[1])));
						// Set flag for counting the number of latitudes later
						// If flag is 0 then longitudes grouped and we calculate delta from latitudes
						topo_flag = 0;
					}
					else
					{
						// Calculate half the distance between points for delta
						delta = std::fabs((0.5)*(longitudes_alt[0] - longitudes_alt[1]));
						// Set flag for counting number of longitudes later
						// If flag is 1 then latitudes are grouped and we calculate delta from longitudes
						topo_flag = 1;
					}
					std::cout << ""<< std::endl;
					std::cout<<"Topography file delta = "<< delta << std::endl;
					std::cout<<"Resolution of input topography in meters is approximately "<< delta*111*2 << std::endl;
					std::cout << ""<< std::endl;
					count++;
				}
				//end of section on calculating delta
			}
		}
		else
			throw std::ios_base::failure("Cannot open file topography file " + filename + "!!!");

		//Calculate the number of unique longitudes or latitudes
		double a,b;
		int count2;

		if ( topo_flag == 1 )
		{
			a = latitudes_alt[0];
			b = latitudes_alt[1];
			count2 = 2;

			while (a-b < 1e-9)
			{
				a = b;
				b = latitudes_alt[count2];
				count2++;
			}
		}
		if ( topo_flag == 0 )
		{
			a = longitudes_alt[0];
			b = longitudes_alt[1];
			count2 = 2;

			while (a-b < 1e-9)
			{
				a = b;
				b = longitudes_alt[count2];
				count2++;
			}
		}
		std::cout<<"number of unique latitudes or longitudes = "<< count2 - 1 << std::endl;
		number_coords = count2-1;
	}


	/**
	 * Look up the altitude for a given latitude/longitude.
	 *
	 * @param latitude The latitude in degrees, counted north from
	 *   the equator.
	 * @param longitude The longitude in degrees, counted east from
	 *   the Greenwich meridian.
	 * @return The altitude at this location, in meters.
	 */
	double
	get_altitude(const double latitude, const double longitude) const
	{
		// loop over the entire array and see if we find a point
		// that's within delta of what we're looking for. the data
		// is arranged in a way that keeps the latitude constant
		// while running over all longitudes, and when we're done
		// with that changes the latitude by one step. so if the
		// latitude is wrong, we can definitely skip ahead a whole
		// set of longitudes/latitudes, which we calculate above and
		// name number_coords
		for (unsigned int i = 0; i < latitudes_alt.size();)
			if (std::fabs(latitude - latitudes_alt[i]) <= delta)
			{

				if (std::fabs(longitude - longitudes_alt[i]) <= delta)
					return altitudes[i] * 2;
				else
					++i;
			}
			else
				i += number_coords;
		std::cout<<"Is your topography file ordered latitude, longitude, value? " << std::endl;
		AssertThrow(false, ExcInternalError());
		return 0;
	}
	// }
private:
	/**
	 * Read and access topography file topography.txt
	 */
	std::vector<double> latitudes_alt;
	std::vector<double> longitudes_alt;
	std::vector<double> altitudes;
};


/**
 * A class that describes the boundary of a mesh that is given by the
 * topology provided for the East African Rift region.
 */
class AfricaTopographyBoundary : public StraightBoundary<3>
{
public:
	/**
	 * Constructor.
	 * @param topography A pointer to a topography object.
	 */
	AfricaTopographyBoundary(const std_cxx1x::shared_ptr<const Topography> &topography)
:
	topography(topography)
{
}

	/**
	 * Refer to the general documentation of
	 * this class and the documentation of the
	 * base class.
	 */
	virtual Point<3>
	get_new_point_on_line(const Triangulation<3>::line_iterator &line) const
	{
		// find lat/long of the mid point of the line, then compute the new mid
		// point by projecting it out onto the real surface
		const std::pair<double, double> lat_long =
				lat_long_from_xyz_wgs84((line->vertex(0) + line->vertex(1)) / 2);

		return xyz_from_lat_long_height_in_wgs84(lat_long.first,
				lat_long.second,
				topography->get_altitude(lat_long.first,
						lat_long.second));
	}

	/**
	 * Refer to the general documentation of
	 * this class and the documentation of the
	 * base class.
	 */
	virtual Point<3>
	get_new_point_on_quad(const Triangulation<3>::quad_iterator &quad) const
	{
		// find lat/long of the mid point of the quad, then compute the new mid
		// point by projecting it out onto the real surface
		const std::pair<double, double> lat_long =
				lat_long_from_xyz_wgs84((quad->vertex(0) + quad->vertex(1)
						+ quad->vertex(2) + quad->vertex(3))
						/ 4);

		return xyz_from_lat_long_height_in_wgs84(lat_long.first,
				lat_long.second,
				topography->get_altitude(lat_long.first,
						lat_long.second));
	}

	/**
	 * Refer to the general
	 * documentation of this class
	 * and the documentation of the
	 * base class.
	 *
	 * Calls
	 * @p get_intermediate_points_between_points.
	 */
	virtual
	void
	get_intermediate_points_on_line(const Triangulation<3>::line_iterator &line,
			std::vector<Point<3> > &points) const
	{
		if (points.size() == 1)
			points[0] = get_new_point_on_line(line);
		else
			get_intermediate_points_between_points(line->vertex(0),
					line->vertex(1),
					points);
	}

	/**
	 * Refer to the general
	 * documentation of this class
	 * and the documentation of the
	 * base class.
	 *
	 * Only implemented for <tt>dim=3</tt>
	 * and for <tt>points.size()==1</tt>.
	 */
	virtual
	void
	get_intermediate_points_on_quad(const Triangulation<3>::quad_iterator &quad,
			std::vector<Point<3> > &points) const
	{
		if (points.size() == 1)
			points[0] = get_new_point_on_quad(quad);
		else
		{
			unsigned int m =
					static_cast<unsigned int>(std::sqrt(static_cast<double>(points.size())));
			Assert(points.size() == m * m, ExcInternalError());

			std::vector<Point<3> > lp0(m);
			std::vector<Point<3> > lp1(m);

			get_intermediate_points_on_line(quad->line(0), lp0);
			get_intermediate_points_on_line(quad->line(1), lp1);

			std::vector<Point<3> > lps(m);
			for (unsigned int i = 0; i < m; ++i)
			{
				get_intermediate_points_between_points(lp0[i], lp1[i], lps);

				for (unsigned int j = 0; j < m; ++j)
					points[i * m + j] = lps[j];
			}
		}
	}

	/**
	 * Implementation of the function
	 * declared in the base class.
	 *
	 * Refer to the general
	 * documentation of this class
	 * and the documentation of the
	 * base class.
	 */
	virtual Tensor<1, 3>
	normal_vector(const Triangulation<3>::face_iterator &face,
			const Point<3> &p) const
			{
		// assume radial normal vectors
		return p / p.norm();
			}

	/**
	 * Compute the normals to the
	 * boundary at the vertices of
	 * the given face.
	 *
	 * Refer to the general
	 * documentation of this class
	 * and the documentation of the
	 * base class.
	 */
	virtual
	void
	get_normals_at_vertices(const Triangulation<3>::face_iterator &face,
			Boundary<3,3>::FaceVertexNormals &face_vertex_normals) const
	{
		Assert(false, ExcNotImplemented());
	}

private:

	/**
	 * Called by
	 * @p get_intermediate_points_on_line
	 * and by
	 * @p get_intermediate_points_on_quad.
	 *
	 * Refer to the general
	 * documentation of
	 * @p get_intermediate_points_on_line
	 * in the documentation of the
	 * base class.
	 */
	void
	get_intermediate_points_between_points(const Point<3> &p0,
			const Point<3> &p1,
			std::vector<Point<3> > &points) const
	{
		const unsigned int n = points.size();
		Assert(n > 0, ExcInternalError());

		for (unsigned int i = 0; i < n; ++i)
		{
			// first compute the point along a straight line
			points[i] = (1 - 1. * (i + 1) / (n + 1)) * p0
					+ 1. * (i + 1) / (n + 1) * p1;
			// then pull it out to the surface
			const std::pair<double, double> lat_long =
					lat_long_from_xyz_wgs84(points[i]);

			points[i] =
					xyz_from_lat_long_height_in_wgs84(lat_long.first,
							lat_long.second,
							topography->get_altitude(lat_long.first,
									lat_long.second));
		}
	}

private:
	/**
	 * Pointer to an object that describes the topography.
	 */
	std_cxx1x::shared_ptr<const Topography> topography;
};

/**
 * A class that describes a geometry for the region below the East African Rift.
 */
template <int dim>
class Africa : public Interface<dim>
{
public:
	/**
	 * Generate a coarse mesh for the geometry described by this class.
	 */
	virtual
	void
	create_coarse_mesh(parallel::distributed::Triangulation<dim> &coarse_grid) const;

	/**
	 * Return the typical length scale one would expect of features in this geometry,
	 * assuming realistic parameters.
	 *
	 * We return 1/20th of the distance from the two endpoints (vertex_0 and vertex_7) of the cell.
	 */
	virtual
	double
	length_scale() const;

	virtual
	double
	depth(const Point<dim> &position) const;

	virtual Point<dim>
	representative_point(const double depth) const;

	virtual
	double
	maximal_depth() const;

	/**
	 * Return the set of boundary indicators that are used by this model. This
	 * information is used to determine what boundary indicators can be used in
	 * the input file.
	 *
	 * The box model uses boundary indicators zero through 2*dim-1, with the first
	 * two being the faces perpendicular to the x-axis, the next two perpendicular
	 * to the y-axis, etc.
	 */
	virtual std::set<types::boundary_id>
	get_used_boundary_indicators() const;

	/**
	 * Declare the parameters this class takes through input files.
	 */
	static
	void
	declare_parameters(ParameterHandler &prm);

	/**
	 * Read the parameters this class declares from the parameter
	 * file.
	 */
	virtual
	void
	parse_parameters(ParameterHandler &prm);

private:
	/**
	 * Pointer to an object that describes the topography.
	 */
	std_cxx1x::shared_ptr<const Topography> topography;
};

template <>
void
Africa<3>::create_coarse_mesh(parallel::distributed::Triangulation<3> &coarse_grid) const
{
	const int dim = 3;

	// create the locations of the vertices
	std::vector<Point<dim> > vertices(12);

	/* Define vertices at the bottom of the cell. these are located
	 * on a sphere */
	for (unsigned int vertex = 0; vertex < 6; ++vertex)
		vertices[vertex] =
				xyz_from_lat_long_height(DomainData::latitudes[vertex],
						DomainData::longitudes[vertex],
						-DomainData::bottom_depth);

	/* Define vertices at the top of the cell. these are located at
	 * certain depths relative to the WGS84 ellipsoid. */
	for (unsigned int vertex = 0; vertex < 6; ++vertex)
		vertices[vertex + 6] =
				xyz_from_lat_long_height_in_wgs84(DomainData::latitudes[vertex],
						DomainData::longitudes[vertex],
						topography->get_altitude(DomainData::latitudes[vertex],
								DomainData::longitudes[vertex]));

	// now define the two cells from it
	const unsigned int vertex_indices[2][8] =
	{
			{ 0, 1, 2, 3, 6, 7, 8, 9 },
			{ 2, 3, 4, 5, 8, 9, 10, 11 } };
	std::vector<CellData<dim> > cells(2);

	for (unsigned int c = 0; c < 2; ++c)
		for (unsigned int vertex = 0;
				vertex < GeometryInfo<dim>::vertices_per_cell; ++vertex)
			cells[c].vertices[vertex] = vertex_indices[c][vertex];

	cells[0].material_id = cells[1].material_id = 0;

	coarse_grid.create_triangulation(vertices, cells, SubCellData());

	// set all boundary indicators. we want the edges to be curved
	// as well. for this, it matters in which order we call
	// set_all_boundary_indicators() -- we have to do it last for
	// the inner and outer boundary, which conveniently is what
	// happens in the following loop
	for (Triangulation<3>::active_cell_iterator cell =
			coarse_grid.begin_active(); cell != coarse_grid.end(); ++cell)
		for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
			if (cell->face(f)->at_boundary())
				cell->face(f)->set_all_boundary_indicators(f);

	static AfricaTopographyBoundary top_boundary(topography);
	coarse_grid.set_boundary(5, top_boundary);

	static HyperShellBoundary<3> bottom_boundary;
	coarse_grid.set_boundary(4, bottom_boundary);
}

template <int dim>
void
Africa<dim>::create_coarse_mesh(parallel::distributed::Triangulation<dim> &coarse_grid) const
{
	Assert(false, ExcNotImplemented());
}

template <int dim>
std::set<types::boundary_id>
Africa<dim>::get_used_boundary_indicators() const
{
	// boundary indicators are zero through 2*dim-1
	std::set<types::boundary_id> s;
	for (unsigned int i = 0; i < 2 * dim; ++i)
		s.insert(i);
	return s;
}

template <int dim>
double
Africa<dim>::length_scale() const
{
	// diameter divided by 20
	return ((xyz_from_lat_long_height_in_wgs84(-10,
			26,
			-DomainData::bottom_depth)
			- xyz_from_lat_long_height_in_wgs84(5, 35, 0)).norm()
			/ 20);
}

template <int dim>
double
Africa<dim>::depth(const Point<dim> &position) const
{
	//FIXME: need to use actual ellipsoid
	const double radius = 6378137;
	const double d = radius - position.norm();

	return std::max(std::min(d, maximal_depth()), 0.0);
}

template <int dim>
double
Africa<dim>::maximal_depth() const
{
	return DomainData::bottom_depth;
}

// Declare variables for reading in coordinates of the region of interest.
double westLongitude;
double eastLongitude;
double northLatitude;
double southLatitude;
double bottomDepth;

template <int dim>
void
Africa<dim>::declare_parameters(ParameterHandler &prm)
{
	prm.enter_subsection("Geometry model");
	{
		prm.enter_subsection("Africa");
		{
			prm.declare_entry("Topography",
					"topography.txt",
					Patterns::FileName(),
					"Surface coordinates and topography. Units: degrees and meters.");
			prm.declare_entry("west longitude",
					"35.0",
					Patterns::Double(),
					"Western longitude of model region.");
			prm.declare_entry("east longitude",
					"26.0",
					Patterns::Double(),
					"Eastern longitude of model region");
			prm.declare_entry("north latitude",
					"5.0",
					Patterns::Double(),
					"Northern latitude of model region.");
			prm.declare_entry("south latitude",
					"-10.0",
					Patterns::Double(),
					"Southern latitude of model region.");
			prm.declare_entry("bottom depth",
					"50000.0",
					Patterns::Double(0),
					"Bottom depth of model region.");
		}
		prm.leave_subsection();
	}
	prm.leave_subsection();
}

template <int dim>
void
Africa<dim>::parse_parameters(ParameterHandler &prm)
{
	prm.enter_subsection("Geometry model");
	{
		prm.enter_subsection("Africa");
		{
			topography.reset(new Topography(prm.get("Topography")));

			/** Get latitude and longitudes defining region of interest from
			 * the parameter file.
			 */
			westLongitude = prm.get_double("west longitude");
			eastLongitude = prm.get_double("east longitude");
			northLatitude = prm.get_double("north latitude");
			southLatitude = prm.get_double("south latitude");
			bottomDepth = prm.get_double("bottom depth");

			//check if we can do calculations within this function
			// bottomDepth = bottomDepth + 1.0;
			// std::cout << "bottomDepth: " << bottomDepth << std::endl;
			//std::cout << "latitudes[3]: " << DomainData::latitudes[3] << std::endl;
			//

			/** For each vertex assign the input value or calculate the mid-point
			 * latitude.
			 */
			DomainData::latitudes[0] = southLatitude;
			DomainData::latitudes[1] = southLatitude;
			DomainData::latitudes[2] = (1.0/2.0)*(southLatitude + northLatitude);
			DomainData::latitudes[3] = (1.0/2.0)*(southLatitude + northLatitude);
			DomainData::latitudes[4] = northLatitude;
			DomainData::latitudes[5] = northLatitude;

			// Print the latitudes to verify your input coordinates
			std::cout << ""<< std::endl;
			std::cout << "Southern, middle, and northern latitudes = "<< std::endl;
			for (int i = 0; i<6; i++){
				string String = static_cast<std::ostringstream*>( &(std::ostringstream() << i))->str();
				std::cout << "latitudes["+ String +"]: " << DomainData::latitudes[i] << std::endl;
			}

			// For each vertex assign the input value
			DomainData::longitudes[0] = westLongitude;
			DomainData::longitudes[1] = eastLongitude;
			DomainData::longitudes[2] = westLongitude;
			DomainData::longitudes[3] = eastLongitude;
			DomainData::longitudes[4] = westLongitude;
			DomainData::longitudes[5] = eastLongitude;

			// Print the longitudes to verif your input coordinates
			std::cout << ""<< std::endl;
			std::cout << "Eastern and western longitudes = "<< std::endl;
			for (int i = 0; i<6; i++){
				string String = static_cast<std::ostringstream*>( &(std::ostringstream() << i))->str();
				std::cout << "longitudes["+ String +"]: " << DomainData::longitudes[i] << std::endl;
			}

			// Assign user input bottom depth of model
			DomainData::bottom_depth = bottomDepth;

			// Print to screen to verify input is read in correctly.
			std::cout << ""<< std::endl;
			std::cout << "Depth of region in kilometers = " << (DomainData::bottom_depth)/1000 << std::endl;
		}
		prm.leave_subsection();
	}
	prm.leave_subsection();
}

template <int dim>
Point<dim>
Africa<dim>::representative_point(const double depth) const
{
	return Point<dim>();
}

template <>
Point<3>
Africa<3>::representative_point(const double depth) const
{
	const int dim = 3;
	const double radius = 6378137;
	Assert(depth >= 0, ExcMessage("Given depth must be positive or zero."));
	Assert(depth <= maximal_depth(),
			ExcMessage("Given depth must be less than or equal to the maximal depth of this geometry."));

	// choose a point on the center axis of the domain
	Point<dim> p =
			//	(xyz_from_lat_long_height_in_wgs84(-10,
			//			26,
			//			-DomainData::bottom_depth)
			//			+ xyz_from_lat_long_height_in_wgs84(5, 35, 0))
			//			/ 2;
			(xyz_from_lat_long_height_in_wgs84(southLatitude,
					eastLongitude,
					-DomainData::bottom_depth)
					+ xyz_from_lat_long_height_in_wgs84(northLatitude, eastLongitude, 0))
					/ 2;
	p /= p.norm();
	p *= radius - depth;
	return p;
}

}
}

// explicit instantiations
namespace aspect
{
namespace GeometryModel
{
ASPECT_REGISTER_GEOMETRY_MODEL(Africa,
		"africa",
		"A 3D geometry that accounts to Earth's ellipticity assuming the WGS84 "
		"ellipsoid definition. The domain stretches from 10 degrees south to "
		"5 degrees north, and from 26 degree east to 35 degrees east, encompassing "
		"the central part of the East African Rift."
		"Faces of model are defined as 0, west; 1,east; 2, south; 3, north; 4, bottom; "
		"5, top.")
}
}

namespace aspect
{
namespace InitialConditions
{
using namespace dealii;

/**
 * A class that implements temperature initial
 * conditions.
 *
 * @ingroup InitialConditionsModels
 */
template <int dim>
class ModelRegions : public Interface<dim>
{
public:
	double delta;
	double delta_crust;
	int litho_flag;
	int crust_flag;
	int number_coords_litho;
	int number_coords_crust;

	/**
	 * Return the initial temperature as a function of position.
	 */
	virtual
	double
	initial_temperature(const Point<dim> &position) const;

	/**
	 * Return the true if within crust and false everywhere else.
	 */
	virtual
	bool
	crustal_region(const Point<dim> &position) const;

	/**
	 * Declare the parameters this class takes through input files.
	 * The default implementation of this function does not describe
	 * any parameters. Consequently, derived classes do not have to
	 * overload this function if they do not take any runtime parameters.
	 */
	static
	void
	declare_parameters(ParameterHandler &prm);

	/**
	 * Read the parameters this class declares from the parameter
	 * file. The default implementation of this function does not read
	 * any parameters. Consequently, derived classes do not have to
	 * overload this function if they do not take any runtime parameters.
	 */
	virtual
	void
	parse_parameters(ParameterHandler &prm);

private:

	/**
	 * Return the depth at which the temperature is 1673.15 K as a function of position.
	 */
	double
	get_lithosphere_isotherm(const double latitude,
			const double longitude) const;

	/**
	 * Find depth of crust.
	 */
	double
	get_crustal_depth(const double latitude,
			const double longitude) const;

	/**
	 * Read and access depth to isotherm file thickness.txt
	 */
	std::vector<double> latitudes_iso;
	std::vector<double> longitudes_iso;
	std::vector<double> depths_iso;

	/**
	 * Read and access depth to crustal file crust_thickness.txt
	 */
	std::vector<double> latitudes_crust;
	std::vector<double> longitudes_crust;
	std::vector<double> depths_crust;
};

template <>
double
ModelRegions<3>::get_lithosphere_isotherm(const double latitude,
		const double longitude) const
		{
	// loop over the entire array and see if we find a point
	// that's within delta of what we're looking for. the data
	// is arranged in a way that keeps the latitude constant
	// while running over all longitudes, and when we're done
	// with that changes the latitude by one step. so if the
	// latitude is wrong, we can definitely skip ahead a whole
	// set of longitudes. The number of values to skip is calculated.
	for (unsigned int i = 0; i <= latitudes_iso.size();)
		if (std::fabs(latitude - latitudes_iso[i]) <= delta)
		{
			if (std::fabs(longitude - longitudes_iso[i]) <= delta)
				return -depths_iso[i]*1000;
			else
				++i;
		}
		else
			i += number_coords_litho;

	Assert(false, ExcInternalError());
	return 0;
		}

template <>
double
ModelRegions<3>::get_crustal_depth(const double latitude,
		const double longitude) const
		{
	// loop over the entire array and see if we find a point
	// that's within delta of what we're looking for. the data
	// is arranged in a way that keeps the latitude constant
	// while running over all longitudes, and when we're done
	// with that changes the latitude by one step. so if the
	// latitude is wrong, we can definitely skip ahead a whole
	// set of longitudes. The number of values to skip is calculated.
	for (unsigned int i = 0; i <= latitudes_crust.size();)
		if (std::fabs(latitude - latitudes_crust[i]) <= delta_crust)
		{
			if (std::fabs(longitude - longitudes_crust[i]) <= delta_crust)
				return -depths_crust[i]*1000;
			else
				++i;
		}
		else
			i += number_coords_litho;

	Assert(false, ExcInternalError());
	return 0;
		}

template <int dim>
double
ModelRegions<dim>::initial_temperature (const Point<dim> &position) const
{
	Assert (false, ExcNotImplemented());
	return 0;
}

template <>
double
ModelRegions<3>::initial_temperature (const Point<3> &position) const
{
	// get the depth of the lithosphere isotherm for the current lat/long
	// position
	const std::pair<double, double> lat_long = lat_long_from_xyz_wgs84(position);

	//TODO: this is the depth with respect to the sphere; need WGS84 here
	const double radius = 6378137;
	const double depth = position.norm() - radius;

	// if above the isotherm, use a linear behavior up to the surface
	// below the isotherm, use an increase of .5 degrees per kilometer
	// (0.0005 degrees per meter)
	const double isotherm_temp = 1673.15;

	const double isotherm_depth = get_lithosphere_isotherm(lat_long.first, lat_long.second);
	if (depth < isotherm_depth)
		return isotherm_temp - (depth-isotherm_depth) * 0.0005;
	else
		return 273.15+depth/isotherm_depth*1400;
}

template <int dim>
bool
ModelRegions<dim>::crustal_region (const Point<dim> &position) const
{
	Assert (false, ExcNotImplemented());
	return 0;
}

template <>
bool
ModelRegions<3>::crustal_region (const Point<3> &position) const
{
	// get the depth of the Mohovorhic discontinuity for the current lat/long
	// position
	const std::pair<double, double> lat_long = lat_long_from_xyz_wgs84(position);

	//TODO: this is the depth with respect to the sphere; need WGS84 here
	const double radius = 6378137;
	const double depth = position.norm() - radius;

	// if above or equal to the Moho, assign region as true, else false

	const double crustal_depth = get_crustal_depth(lat_long.first, lat_long.second);
	if (depth <= crustal_depth)
		return true;
	else
		return false;
}

template <int dim>
void
ModelRegions<dim>::declare_parameters(ParameterHandler &prm)
{
	prm.enter_subsection("Initial conditions");
	{
		prm.enter_subsection("Model regions");
		{
			prm.declare_entry("Isotherm filename",
					"thickness.txt",
					Patterns::FileName(),
					"Surface coordinates and depths to the 1673.15 K isotherm. Units: degrees and kilometers.");

			prm.declare_entry("Crustal thickness filename",
					"crustal_thickness.txt",
					Patterns::FileName(),
					"Surface coordinates and depths to the Mohovoric discontinuity. Units: degrees and kilometers.");
		}
		prm.leave_subsection();
	}
	prm.leave_subsection();
}

template <int dim>
void
ModelRegions<dim>::parse_parameters(ParameterHandler &prm)

{
	std::string isotherm_file;
	std::string crustal_file;
	prm.enter_subsection("Initial conditions");
	{
		prm.enter_subsection("Model regions");
		{
			isotherm_file = prm.get("Isotherm filename");
			crustal_file = prm.get("Crustal thickness filename");
		}
		prm.leave_subsection();
	}
	prm.leave_subsection();

	std::ifstream input1(isotherm_file.c_str());
	AssertThrow (input1.is_open(),
			ExcMessage (std::string("Can't read from file <") + isotherm_file + ">"));
	std::string line;

	std::ifstream input2(crustal_file.c_str());
	AssertThrow (input2.is_open(),
			ExcMessage (std::string("Can't read from file <") + crustal_file + ">"));

	// Start loop with int count to calculate delta below for lithospheric thickness
	int count = 0;
	while (true)
	{
		double latitude_iso, longitude_iso, depth_iso;
		input1 >> latitude_iso >> longitude_iso >> depth_iso;
		if (input1.eof())
			break;

		latitudes_iso.push_back(latitude_iso);
		longitudes_iso.push_back(longitude_iso);
		depths_iso.push_back(depth_iso);

		/** Find first 2 numbers that are different to use in
		 * calculating half the difference between each position as delta.
		 */
		if (count < 2 )
		{
			count ++;
		}
		if (count == 2)
		{
			if ((std::fabs(latitudes_iso[0] - latitudes_iso[1]) > 1e-9) && (std::fabs(longitudes_iso[0] - longitudes_iso[1]) > 1e-9))
			{
				// Stop program if file formatted incorrectly.
				std::cout << ""<< std::endl;
				throw std::ios_base::failure("Lithospheric thickness file not formatted correctly. " + isotherm_file + "Make sure you have lat, lon, value with lat. or lon. varying.");
			}
			if ((std::fabs(latitudes_iso[0] - latitudes_iso[1]) < 1e-9) && (std::fabs(longitudes_iso[0] - longitudes_iso[1]) < 1e-9))
			{
				// Stop program if file formatted incorrectly.
				std::cout << ""<< std::endl;
				throw std::ios_base::failure("Lithospheric thickness file not formatted correctly. " + isotherm_file + "Make sure you have lat, lon, value with lat. or lon. varying.");
			}

			if (std::fabs(latitudes_iso[0] - latitudes_iso[1]) > 1e-9)
			{
				// Calculate delta as half the distance between points.
				delta = std::fabs((0.5)*(latitudes_iso[0] - latitudes_iso[1]));
				// If flag is 0 then longitudes grouped and we calculate delta from latitudes
				litho_flag = 0;
			}
			else
			{
				// Calculate delta as half the distance between points.
				delta = std::fabs((0.5)*(longitudes_iso[0] - longitudes_iso[1]));
				// If flag is 1 then latitudes are grouped and we calculate delta from longitudes
				litho_flag = 1;
			}
			std::cout << ""<< std::endl;
			std::cout<<"Lithosphere thickness delta = "<< delta << std::endl;
			std::cout<<"Resolution of input lithosphere thickness in meters is approximately = "<< delta*111*2 << std::endl;
			std::cout << ""<< std::endl;
			count++;
		}
	} // End of loop for calculating delta for the lithosphere thickness file


	// Start loop with int countc to calculate delta_crust_crust below
	int countc = 0;
	while (true)
	{
		double latitude_crust, longitude_crust, depth_crust;
		input2 >> latitude_crust >> longitude_crust >> depth_crust;
		if (input2.eof())
			break;

		latitudes_crust.push_back(latitude_crust);
		longitudes_crust.push_back(longitude_crust);
		depths_crust.push_back(depth_crust);


		/** Find first 2 numbers that are different to use in
		 * calculating half the difference between each position as delta_crust.
		 */
		if (countc < 2 )
		{
			countc ++;
		}
		if (countc == 2)
		{

			if ((std::fabs(latitudes_crust[0] - latitudes_crust[1]) > 1e-9) && (std::fabs(longitudes_crust[0] - longitudes_crust[1]) > 1e-9))
			{
				// Stop program if file formatted incorrectly.
				std::cout << ""<< std::endl;
				throw std::ios_base::failure("Lithospheric thickness file not formatted correctly. " + isotherm_file + "Make sure you have lat, lon, value with lat. or lon. varying.");
			}

			if ((std::fabs(latitudes_crust[0] - latitudes_crust[1]) < 1e-9) && (std::fabs(longitudes_crust[0] - longitudes_crust[1]) < 1e-9))
			{
				// Stop program if file formatted incorrectly.
				std::cout << ""<< std::endl;
				throw std::ios_base::failure("Lithospheric thickness file not formatted correctly. " + isotherm_file + "Make sure you have lat, lon, value with lat. or lon. varying.");
			}

			if (std::fabs(latitudes_crust[0] - latitudes_crust[1]) > 1e-9)
			{
				// Calculate delta_crust as half the distance between points.
				delta_crust = std::fabs((0.5)*(latitudes_crust[0] - latitudes_crust[1]));
				// If flag is 0 then longitudes grouped and we calculate delta_crust from latitudes
				crust_flag = 0;
			}
			else
			{
				// Calculate delta_crust as half the distance between points.
				delta_crust = std::fabs((0.5)*(longitudes_crust[0] - longitudes_crust[1]));
				// If flag is 1 then latitudes are grouped and we calculate delta_crust from longitudes
				crust_flag = 1;
			}

			std::cout << ""<< std::endl;
			std::cout<<"Crustal thickness delta_crust = "<< delta_crust << std::endl;
			std::cout<<"Resolution of input Crustal thickness in meters is approximately = "<< delta_crust*111*2 << std::endl;
			std::cout << ""<< std::endl;

			countc++;
		}
	} //End loop for calculate delta for crustal thickness.

	//Calculate the number of unique longitudes or latitudes from the lithosphere isotherm file and crustal thickness file.
	double c,d,r,s;
	int count3,count5;
	count3 = 0;
	count5 = 0;
	c = 0;
	d = 0;
	r = 0;
	s = 0;

	if ( litho_flag == 1 )
	{
		c = latitudes_iso[0];
		d = latitudes_iso[1];
		count3 = 2;

		while (c-d < 1e-9)
		{
			c = d;
			d = latitudes_iso[count3];
			count3++;
		}
	}

	if ( crust_flag == 1 )
	{
		r = latitudes_crust[0];
		s = latitudes_crust[1];
		count5 = 2;

		while (r-s < 1e-9)
		{
			r = s;
			s = latitudes_crust[count5];
			count5++;
		}
	}

	if ( litho_flag == 0 )
	{
		c = longitudes_iso[0];
		d = longitudes_iso[1];
		count3 = 2;

		while (c-d < 1e-9)
		{
			c = d;
			d = longitudes_iso[count3];
			count3++;
		}
	}

	if ( crust_flag == 0 )
	{
		r = longitudes_crust[0];
		s = longitudes_crust[1];
		count5 = 2;

		while (r-s < 1e-9)
		{
			r = s;
			s = longitudes_crust[count5];
			count5++;
		}
	}
	std::cout << ""<< std::endl;
	std::cout<<"number of unique latitudes or longitudes in lithosphere thickness file= "<< count3 - 1 << std::endl;
	number_coords_litho = count3-1;
	std::cout << ""<< std::endl;
	std::cout<<"number of unique latitudes or longitudes in crustal thickness file= "<< count5 - 1 << std::endl;
	number_coords_crust = count3-1;
	std::cout << ""<< std::endl;
}
}
}

// explicit instantiations
namespace aspect
{
namespace InitialConditions
{
ASPECT_REGISTER_INITIAL_CONDITIONS(ModelRegions,
		"model regions",
		"In subsection lithosphere isotherm we define "
		"the extent of the conductive heat "
		"equation is assumed to be 1400 C (1673.15 K) "
		"as previously used by Bird et al., 2008 "
		"and Stamps et al. (in prep). This assumption "
		"is consistent with the Schubert et al., 2004 "
		"definition of the mechanical lithosphere.")
}
}

namespace aspect
{
namespace MaterialModel
{
using namespace dealii;

/**
 * A material model that consists of globally constant values for all
 * material parameters except that the density decays linearly with the
 * temperature.
 *
 * The model is considered incompressible, following the definition
 * described in Interface::is_compressible. This is essentially
 * the material model used in the step-32 tutorial program.
 *
 * @ingroup MaterialModels
 */
template <int dim>
class Stamps : public MaterialModel::InterfaceCompatibility<dim>, public ::aspect::SimulatorAccess<dim>
{
public:
	/**
	 * @name Physical parameters used in the basic equations
	 * @{
	 */
	virtual double viscosity (const double                  temperature,
			const double                  pressure,
			const std::vector<double>    &compositional_fields,
			const SymmetricTensor<2,dim> &strain_rate,
			const Point<dim>             &position) const;

	virtual double density (const double temperature,
			const double pressure,
			const std::vector<double> &compositional_fields,
			const Point<dim> &position) const;

	virtual double compressibility (const double temperature,
			const double pressure,
			const std::vector<double> &compositional_fields,
			const Point<dim> &position) const;

	virtual double specific_heat (const double temperature,
			const double pressure,
			const std::vector<double> &compositional_fields,
			const Point<dim> &position) const;

	virtual double thermal_expansion_coefficient (const double      temperature,
			const double      pressure,
			const std::vector<double> &compositional_fields,
			const Point<dim> &position) const;

	virtual double thermal_conductivity (const double temperature,
			const double pressure,
			const std::vector<double> &compositional_fields,
			const Point<dim> &position) const;

	/**
	 * @}
	 */

	/**
	 * @name Qualitative properties one can ask a material model
	 * @{
	 */

	/**
	 * Return true if the viscosity() function returns something that
	 * may depend on the variable identifies by the argument.
	 */
	virtual bool
	viscosity_depends_on (const NonlinearDependence::Dependence dependence) const;

	/**
	 * Return true if the density() function returns something that
	 * may depend on the variable identifies by the argument.
	 */
	virtual bool
	density_depends_on (const NonlinearDependence::Dependence dependence) const;

	/**
	 * Return true if the compressibility() function returns something that
	 * may depend on the variable identifies by the argument.
	 *
	 * This function must return false for all possible arguments if the
	 * is_compressible() function returns false.
	 */
	virtual bool
	compressibility_depends_on (const NonlinearDependence::Dependence dependence) const;

	/**
	 * Return true if the specific_heat() function returns something that
	 * may depend on the variable identifies by the argument.
	 */
	virtual bool
	specific_heat_depends_on (const NonlinearDependence::Dependence dependence) const;

	/**
	 * Return true if the thermal_conductivity() function returns something that
	 * may depend on the variable identifies by the argument.
	 */
	virtual bool
	thermal_conductivity_depends_on (const NonlinearDependence::Dependence dependence) const;

	/**
	 * Return whether the model is compressible or not.  Incompressibility
	 * does not necessarily imply that the density is constant; rather, it
	 * may still depend on temperature or pressure. In the current
	 * context, compressibility means whether we should solve the continuity
	 * equation as $\nabla \cdot (\rho \mathbf u)=0$ (compressible Stokes)
	 * or as $\nabla \cdot \mathbf{u}=0$ (incompressible Stokes).:L
	 */
	virtual bool is_compressible () const;
	/**
	 * @}
	 */

	/**
	 * @name Reference quantities
	 * @{
	 */
	virtual double reference_viscosity () const;

	virtual double reference_density () const;

	virtual double reference_thermal_expansion_coefficient () const;

	double reference_thermal_diffusivity () const;

	double reference_cp () const;
	/**
	 * @}
	 */

	/** Variables used to read in crustal thickness file.
	 *
	 */
	double delta_crust;
	int crust_flag;
	int number_coords_crust;

	/**
	 * Return the true if within crust and false everywhere else.
	 */
	virtual bool crustal_region (const Point<dim> &position) const;

	/**
	 * Declare the parameters this class takes through input files.
	 */
	static
	void
	declare_parameters (ParameterHandler &prm);

	/**
	 * Read the parameters this class declares from the parameter
	 * file.
	 */
	virtual
	void
	parse_parameters (ParameterHandler &prm);
	/**
	 * @}
	 */

private:
	double reference_rho;
	double reference_T;
	double eta;
	double composition_viscosity_prefactor;
	double thermal_viscosity_exponent;
	double thermal_alpha;
	double reference_specific_heat;

	/**
	 * The thermal conductivity.
	 */
	double k_value;
	double compositional_delta_rho;

	/**
	 * Find depth of crust.
	 */
	double
	get_crustal_depth(const double latitude,
			const double longitude) const;

	/**
	 * Read and access depth to crustal file crust_thickness.txt
	 */
	std::vector<double> latitudes_crust;
	std::vector<double> longitudes_crust;
	std::vector<double> depths_crust;
};

}
}

namespace aspect
{
namespace MaterialModel
{

template <>
double
Stamps<3>::get_crustal_depth(const double latitude,
		const double longitude) const
		{
	// loop over the entire array and see if we find a point
	// that's within delta of what we're looking for. the data
	// is arranged in a way that keeps the latitude constant
	// while running over all longitudes, and when we're done
	// with that changes the latitude by one step. so if the
	// latitude is wrong, we can definitely skip ahead a whole
	// set of longitudes. The number of values to skip is calculated.
	for (unsigned int i = 0; i <= latitudes_crust.size();)
		if (std::fabs(latitude - latitudes_crust[i]) <= delta_crust)
		{
			if (std::fabs(longitude - longitudes_crust[i]) <= delta_crust)
				return -depths_crust[i]*1000;
			else
				++i;
		}
		else
			i += number_coords_crust;

	Assert(false, ExcInternalError());
	return 0;
		}

template <int dim>
bool
Stamps<dim>::crustal_region (const Point<dim> &position) const
{
	Assert (false, ExcNotImplemented());
	return 0;
}

template <>
bool
Stamps<3>::crustal_region (const Point<3> &position) const
{
	// get the depth of the Mohovorhic discontinuity for the current lat/long
	// position
	const std::pair<double, double> lat_long = lat_long_from_xyz_wgs84(position);

	//TODO: this is the depth with respect to the sphere; need WGS84 here
	const double radius = 6378137;
	const double depth = position.norm() - radius;

	// if above or equal to the Moho, assign region as true, else false

	const double crustal_depth = get_crustal_depth(lat_long.first, lat_long.second);
	if (depth <= crustal_depth)
		return true;
	else
		return false;
}

template <int dim>
double
Stamps<dim>::
viscosity (const double temperature,
		const double pressure,
		const std::vector<double> &composition,       /*composition*/
		const SymmetricTensor<2,dim> &,
		const Point<dim> &position) const
		{
	// Material parameters are from Karato and Wu (1993) for dry olivine)
	const double B_diff_m = 810000000000.0; 	// see p.248-240 Schubert et al., 2004; grain size = 3 mm
	const double R = 8.3144;                	// gas constant J/K.mol
	const double V_diff_m = 0.000006;         	// Activation volume m^3/mol
	const double E_diff_m = 300000;           	// J/mol
	const double Biot = 0.01;					// Biot's pore pressure

	if (temperature < 1673.15)
		return 1e25;
	else
		return (0.5 * B_diff_m * std::exp((E_diff_m+pressure*V_diff_m)/(R*temperature)));
		}
//	if (crustal_region(position) == true && temperature < 1673.15)
//
//		return (Coulomb friction law)
//
//				if (crustal_region(position) == false && temperature < 1673.15)
//
//					return (dislocation creep);
//
//				else
//					// Use diffusion creep flow law below the lithosphere
//					return (0.5 * B_diff_m * std::exp((E_diff_m+pressure*V_diff_m)/(R*temperature)));
//
//		}


template <int dim>
double
Stamps<dim>::
reference_viscosity () const
{
	return eta;
}

template <int dim>
double
Stamps<dim>::
reference_density () const
{
	return reference_rho;
}

template <int dim>
double
Stamps<dim>::
reference_thermal_expansion_coefficient () const
{
	return thermal_alpha;
}

template <int dim>
double
Stamps<dim>::
specific_heat (const double,
		const double,
		const std::vector<double> &, /*composition*/
		const Point<dim> &) const
		{
	return reference_specific_heat;
		}

template <int dim>
double
Stamps<dim>::
reference_cp () const
{
	return reference_specific_heat;
}

template <int dim>
double
Stamps<dim>::
thermal_conductivity (const double,
		const double,
		const std::vector<double> &, /*composition*/
		const Point<dim> &) const
		{
	return k_value;
		}

template <int dim>
double
Stamps<dim>::
reference_thermal_diffusivity () const
{
	return k_value/(reference_rho*reference_specific_heat);
}

template <int dim>
double
Stamps<dim>::
density (const double temperature,
		const double,
		const std::vector<double> &compositional_fields, /*composition*/
		const Point<dim> &) const
		{
	return (reference_rho * (1 - thermal_alpha * (temperature - reference_T))
			+
			(compositional_fields.size()>0
					?
							compositional_delta_rho * compositional_fields[0]
							                                               :
							0));
		}


template <int dim>
double
Stamps<dim>::
thermal_expansion_coefficient (const double temperature,
		const double,
		const std::vector<double> &, /*composition*/
		const Point<dim> &) const
		{
	return thermal_alpha;
		}


template <int dim>
double
Stamps<dim>::
compressibility (const double,
		const double,
		const std::vector<double> &, /*composition*/
		const Point<dim> &) const
		{
	return 0.0;
		}

template <int dim>
bool
Stamps<dim>::
viscosity_depends_on (const NonlinearDependence::Dependence dependence) const
{
	// compare this with the implementation of the viscosity() function
	// to see the dependencies
	if (((dependence & NonlinearDependence::temperature) != NonlinearDependence::none)
			&&
			(thermal_viscosity_exponent != 0))
		return true;
	else if (((dependence & NonlinearDependence::compositional_fields) != NonlinearDependence::none)
			&&
			(composition_viscosity_prefactor != 1.0))
		return true;
	else
		return false;
}


template <int dim>
bool
Stamps<dim>::
density_depends_on (const NonlinearDependence::Dependence dependence) const
{
	// compare this with the implementation of the density() function
	// to see the dependencies
	if (((dependence & NonlinearDependence::temperature) != NonlinearDependence::none)
			&&
			(thermal_alpha != 0))
		return true;
	else if (((dependence & NonlinearDependence::compositional_fields) != NonlinearDependence::none)
			&&
			(compositional_delta_rho != 0))
		return true;
	else
		return false;
}

template <int dim>
bool
Stamps<dim>::
compressibility_depends_on (const NonlinearDependence::Dependence) const
{
	return false;
}

template <int dim>
bool
Stamps<dim>::
specific_heat_depends_on (const NonlinearDependence::Dependence) const
{
	return false;
}

template <int dim>
bool
Stamps<dim>::
thermal_conductivity_depends_on (const NonlinearDependence::Dependence dependence) const
{
	return false;
}


template <int dim>
bool
Stamps<dim>::
is_compressible () const
{
	return false;
}

template <int dim>
void
Stamps<dim>::declare_parameters (ParameterHandler &prm)
{
	prm.enter_subsection("Material model");
	{
		prm.enter_subsection("Stamps model");
		{
			prm.declare_entry("Crustal thickness filename",
					"crustal_thickness.txt",
					Patterns::FileName(),
					"Surface coordinates and depths to the Mohovoric discontinuity. Units: degrees and kilometers.");
			prm.declare_entry ("Reference density", "3300",
					Patterns::Double (0),
					"Reference density $\\rho_0$. Units: $kg/m^3$.");
			prm.declare_entry ("Reference temperature", "293",
					Patterns::Double (0),
					"The reference temperature $T_0$. Units: $K$.");
			prm.declare_entry ("Viscosity", "5e24",
					Patterns::Double (0),
					"The value of the constant viscosity. Units: $kg/m/s$.");
			prm.declare_entry ("Composition viscosity prefactor", "1.0",
					Patterns::Double (0),
					"A linear dependency of viscosity on the first compositional field. "
					"Dimensionless prefactor. With a value of 1.0 (the default) the "
					"viscosity does not depend on the composition.");
			prm.declare_entry ("Thermal viscosity exponent", "0.0",
					Patterns::Double (0),
					"The temperature dependence of viscosity. Dimensionless exponent.");
			prm.declare_entry ("Thermal conductivity", "4.7",
					Patterns::Double (0),
					"The value of the thermal conductivity $k$. "
					"Units: $W/m/K$.");
			prm.declare_entry ("Reference specific heat", "1250",
					Patterns::Double (0),
					"The value of the specific heat $cp$. "
					"Units: $J/kg/K$.");
			prm.declare_entry ("Thermal expansion coefficient", "2e-5",
					Patterns::Double (0),
					"The value of the thermal expansion coefficient $\\beta$. "
					"Units: $1/K$.");
			prm.declare_entry ("Density differential for compositional field 1", "0",
					Patterns::Double(),
					"If compositional fields are used, then one would frequently want "
					"to make the density depend on these fields. In this simple material "
					"model, we make the following assumptions: if no compositional fields "
					"are used in the current simulation, then the density is simply the usual "
					"one with its linear dependence on the temperature. If there are compositional "
					"fields, then the density only depends on the first one in such a way that "
					"the density has an additional term of the kind $+\\Delta \\rho \\; c_1(\\mathbf x)$. "
					"This parameter describes the value of $\\Delta \\rho$. Units: $kg/m^3/\\textrm{unit "
					"change in composition}$.");
		}
		prm.leave_subsection();
	}
	prm.leave_subsection();
}



template <int dim>
void
Stamps<dim>::parse_parameters (ParameterHandler &prm)
{
	std::string crustal_file;
	prm.enter_subsection("Material model");
	{
		prm.enter_subsection("Stamps model");
		{
			crustal_file 			   = prm.get("Crustal thickness filename");
			reference_rho              = prm.get_double ("Reference density");
			reference_T                = prm.get_double ("Reference temperature");
			eta                        = prm.get_double ("Viscosity");
			composition_viscosity_prefactor = prm.get_double ("Composition viscosity prefactor");
			thermal_viscosity_exponent = prm.get_double ("Thermal viscosity exponent");
			k_value                    = prm.get_double ("Thermal conductivity");
			reference_specific_heat    = prm.get_double ("Reference specific heat");
			thermal_alpha              = prm.get_double ("Thermal expansion coefficient");
			compositional_delta_rho    = prm.get_double ("Density differential for compositional field 1");

			if (thermal_viscosity_exponent!=0.0 && reference_T == 0.0)
				AssertThrow(false, ExcMessage("Error: Material model Stamps with Thermal viscosity exponent can not have reference_T=0."));
		}
		prm.leave_subsection();
	}
	prm.leave_subsection();
}
}
}

// explicit instantiations
namespace aspect
{
namespace MaterialModel
{
ASPECT_REGISTER_MATERIAL_MODEL(Stamps,
		"stamps",
		"A material model that has constant values "
		"for all coefficients but the density and viscosity. "
		"This model uses the formulation that assumes an incompressible"
		" medium despite the fact that the density follows the law "
		"$\\rho(T)=\\rho_0(1-\\beta(T-T_{\\text{ref}})$. "
		"The temperature dependency of viscosity is "
		" switched off by default and follows the formula"
		"$\\eta(T)=\\eta_0*e^{\\eta_T*\\Delta T / T_{\\text{ref}})}$."
		"The value for the components of this formula and additional "
		"parameters are read from the parameter file in subsection "
		"'Stamps model'.")
}
}

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

#include <aspect/global.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_values.h>


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

	for (; cell!=endc; ++cell)
		if (cell->is_locally_owned())
			for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
				if (cell->at_boundary(f) &&
						cell->face(f)->boundary_indicator() == 5) // face is at top boundary
				{
					fe_face_values.reinit (cell, f);
					fe_face_values[this->introspection().extractors.velocities].get_function_values (this->get_solution(),
							velocity_values);

					for (unsigned int q=0; q<quadrature_formula.size(); ++q)
					{
						Point<3> p (fe_face_values.quadrature_point(q)[0],
								fe_face_values.quadrature_point(q)[1],
								fe_face_values.quadrature_point(q)[2]);
						std::pair<double,double> lat_long = lat_long_from_xyz (p);

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
		"Output surface velocities at the centers of surface faces.")
}
}
