/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#include "src/Geometry.h"

#include <netcdf.h>

#include <cmath>
#include <sstream>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/meshgenerator.h"
#include "atlas/util/Geometry.h"
#include "atlas/util/KDTree.h"
#include "atlas/util/Point.h"

#include "eckit/mpi/Comm.h"

#include "src/Variables.h"

#include "util/abor1_cpp.h"
#include "util/Logger.h"

#define ERR(e) {ABORT(nc_strerror(e));}

// -----------------------------------------------------------------------------
namespace quench {
// -----------------------------------------------------------------------------
Geometry::Geometry(const eckit::Configuration & conf)
  : comm_(eckit::mpi::comm()), gridType_("no_type"), groups_() {
  // Initialize eckit communicator for ATLAS
  eckit::mpi::setCommDefault(comm_.name().c_str());

  // Halo
  halo_ = conf.getInt("halo", 0);

  // Set regional grid flag
  regionalGrid_ = false;

  // Setup grid
  eckit::LocalConfiguration gridParams = conf.getSubConfiguration("grid");
  oops::Log::info() << "Info     : Grid config: " << gridParams << std::endl;
  if (gridParams.has("type")) {
    gridType_ = gridParams.getString("type");
    if (gridType_ == "regional") {
      regionalGrid_ = true;
    }
  }
  grid_ = atlas::Grid(gridParams);

  // Setup partitioner
  const bool noPointOnLastTask = conf.getBool("no point on last task", false);
  if (noPointOnLastTask && (comm_.size() > 1)) {
    partitioner_ = atlas::grid::Partitioner(conf.getString("partitioner", "equal_regions"),
      comm_.size()-1);
  } else {
    partitioner_ = atlas::grid::Partitioner(conf.getString("partitioner", "equal_regions"));
  }

  if (conf.getString("function space") == "StructuredColumns") {
    // StructuredColumns
    ASSERT(gridType_ != "unstructured");
    ASSERT(partitioner_);

    // Setup function space
    functionSpace_ = atlas::functionspace::StructuredColumns(grid_, partitioner_,
                     atlas::option::halo(halo_));

    // Bugfix for regional grids
    if (regionalGrid_) {
      auto lonlat = atlas::array::make_view<double, 2>(functionSpace_.lonlat());
      double lonlatPoint[] = {0, 0};
      atlas::functionspace::StructuredColumns fs(functionSpace_);
      atlas::StructuredGrid grid = fs.grid();
      auto view_i = atlas::array::make_view<int, 1>(fs.index_i());
      auto view_j = atlas::array::make_view<int, 1>(fs.index_j());
      for (int jj = 0; jj < fs.size(); ++jj) {
        grid.lonlat(view_i(jj)-1, view_j(jj)-1, lonlatPoint);
       lonlat(jj, 1) = lonlatPoint[1];
        lonlat(jj, 0) = lonlatPoint[0];
      }
    }

    // Setup mesh
    mesh_ = atlas::MeshGenerator("structured").generate(grid_, partitioner_);
  } else if (conf.getString("function space") == "NodeColumns") {
    // NodeColumns
    if (grid_.name().compare(0, 2, std::string{"CS"}) == 0) {
      // CubedSphere
      ASSERT(partitioner_);
      mesh_ = atlas::MeshGenerator("cubedsphere_dual").generate(grid_, partitioner_);
      functionSpace_ = atlas::functionspace::CubedSphereNodeColumns(mesh_);
    } else {
      // For now allow NodeColumns for unstructured triangulation meshes
      // Regular or Structured grids could be supported with extra code
      ASSERT(partitioner_);
      mesh_ = atlas::MeshGenerator("delaunay").generate(grid_, partitioner_);
      functionSpace_ = atlas::functionspace::NodeColumns(mesh_);
    }
  } else if (conf.getString("function space") == "PointCloud") {
    ABORT(conf.getString("function space") + " function space not supported");
  } else {
    ABORT(conf.getString("function space") + " function space not supported yet");
  }

  // Groups
  size_t groupIndex = 0;
  for (const auto & groupParams : conf.getSubConfigurations("groups")) {
    // Use this group index for all the group variables
    for (const auto & var : groupParams.getStringVector("variables")) {
      if (groupIndex_.find(var) != groupIndex_.end()) {
        ABORT("Same variable present in distinct groups");
      } else {
        groupIndex_[var] = groupIndex;
      }
    }

    // Define group
    groupData group;

    // Number of levels
    group.levels_ = groupParams.getInt("levels", 1);

    // Corresponding level for 2D variables (first or last)
    group.lev2d_ = groupParams.getString("lev2d", "first");

    // Vertical coordinate
    if (groupParams.has("vert_coord")) {
      std::vector<double> vert_coordParams = groupParams.getDoubleVector("vert_coord");
      if (vert_coordParams.size() != group.levels_) {
        ABORT("Wrong number of levels in the user-specified vertical coordinate");
      }
      for (size_t jlevel = 0; jlevel < group.levels_; ++jlevel) {
        group.vert_coord_.push_back(vert_coordParams[jlevel]);
      }
    } else {
      for (size_t jlevel = 0; jlevel < group.levels_; ++jlevel) {
        group.vert_coord_.push_back(static_cast<double>(jlevel+1));
      }
    }

    // Default mask, set to 1 (true)
    atlas::Field gmask = functionSpace_.createField<int>(
      atlas::option::name("gmask") | atlas::option::levels(group.levels_));
    auto maskView = atlas::array::make_view<int, 2>(gmask);
    maskView.assign(1);

    // Specific mask
    std::string maskType = groupParams.getString("mask type", "none");
    if (maskType == "none") {
      // No mask
    } else if (maskType == "sea") {
      // Read sea mask
      readSeaMask(groupParams.getString("mask path", "../quench/data/landsea.nc"),
        group.levels_, group.lev2d_, gmask);
    } else {
      ABORT("Wrong mask type");
    }

    // Fill geometry fields
    group.fields_ = atlas::FieldSet();

    if (regionalGrid_) {
      // 2D indices
      atlas::functionspace::StructuredColumns fs(functionSpace_);
      group.fields_->add(fs.index_i());
      group.fields_->add(fs.index_j());

      // Area
      atlas::Field area = functionSpace_.createField<double>(
      atlas::option::name("area") | atlas::option::levels(1));
      auto areaView = atlas::array::make_view<double, 2>(area);
      atlas::StructuredGrid grid = fs.grid();
      auto view_i = atlas::array::make_view<int, 1>(fs.index_i());
      auto view_j = atlas::array::make_view<int, 1>(fs.index_j());
      for (atlas::idx_t jnode = 0; jnode < area.shape(0); ++jnode) {
        // Initialization
        int i = view_i(jnode);
        int j = view_j(jnode);
        double dist_i = 0.0;
        double dist_j = 0.0;

        // i-direction component
        if (i == 1) {
          dist_i = atlas::util::Earth().distance(grid.lonlat(i-1, j-1), grid.lonlat(i, j-1));
        } else if (i == grid.nx(j-1)) {
          dist_i = atlas::util::Earth().distance(grid.lonlat(i-2, j-1), grid.lonlat(i-1, j-1));
        } else {
          dist_i = 0.5*atlas::util::Earth().distance(grid.lonlat(i-2, j-1), grid.lonlat(i, j-1));
        }

        // j-direction component
        if (j == 1) {
          dist_j = atlas::util::Earth().distance(grid.lonlat(i-1, j-1), grid.lonlat(i-1, j));
        } else if (j == grid.ny()) {
          dist_j = atlas::util::Earth().distance(grid.lonlat(i-1, j-2), grid.lonlat(i-1, j-1));
        } else {
          dist_j = 0.5*atlas::util::Earth().distance(grid.lonlat(i-1, j-2), grid.lonlat(i-1, j));
        }

        // Local scale
        areaView(jnode, 0) = dist_i*dist_j;
      }
      group.fields_->add(area);
    }

    // Vertical coordinate
    atlas::Field vert_coord = functionSpace_.createField<double>(
      atlas::option::name("vert_coord") | atlas::option::levels(group.levels_));
    auto vert_coordView = atlas::array::make_view<double, 2>(vert_coord);
    for (atlas::idx_t jnode = 0; jnode < vert_coord.shape(0); ++jnode) {
      for (size_t jlevel = 0; jlevel < group.levels_; ++jlevel) {
        vert_coordView(jnode, jlevel) = group.vert_coord_[jlevel];
      }
    }
    group.fields_->add(vert_coord);

    // Geographical mask
    group.fields_->add(gmask);

    // Owned points mask
    if (functionSpace_.type() == "StructuredColumns") {
      // Structured columns
      atlas::functionspace::StructuredColumns fs(functionSpace_);
      atlas::StructuredGrid grid = fs.grid();
      atlas::Field owned = fs.createField<int>(atlas::option::name("owned")
        | atlas::option::levels(1));
      auto ownedView = atlas::array::make_view<int, 2>(owned);
      auto ghostView = atlas::array::make_view<int, 1>(fs.ghost());
      auto view_i = atlas::array::make_view<int, 1>(fs.index_i());
      auto view_j = atlas::array::make_view<int, 1>(fs.index_j());
      for (atlas::idx_t j = fs.j_begin_halo(); j < fs.j_end_halo(); ++j) {
        for (atlas::idx_t i = fs.i_begin_halo(j); i < fs.i_end_halo(j); ++i) {
          atlas::idx_t jnode = fs.index(i, j);
          ownedView(jnode, 0) = ghostView(jnode) > 0 ? 0 : 1;
        }
      }

      // Special case for lon/lat grids
      if (grid_.name().compare(0, 1, std::string{"L"}) == 0) {
        for (atlas::idx_t j = fs.j_begin_halo(); j < fs.j_end_halo(); ++j) {
          for (atlas::idx_t i = fs.i_begin_halo(j); i < fs.i_end_halo(j); ++i) {
            atlas::idx_t jnode = fs.index(i, j);
            if (((view_j(jnode) == 1) || (view_j(jnode) == grid.ny())) && (view_i(jnode) != 1)) {
              ownedView(jnode, 0) = 0;
            }
          }
        }
      }

      // Add owned points mask
      group.fields_->add(owned);
    } else if (functionSpace_.type() == "NodeColumns") {
      // NodeColumns
      if (grid_.name().compare(0, 2, std::string{"CS"}) == 0) {
        // CubedSphere
        atlas::functionspace::NodeColumns fs(functionSpace_);
        atlas::Field owned = fs.createField<int>(atlas::option::name("owned")
          | atlas::option::levels(1));
        auto ownedView = atlas::array::make_view<int, 2>(owned);
        auto ghostView = atlas::array::make_view<int, 1>(fs.ghost());
        for (atlas::idx_t jnode = 0; jnode < owned.shape(0); ++jnode) {
          ownedView(jnode, 0) = ghostView(jnode) > 0 ? 0 : 1;
        }

        // Add owned points mask
        group.fields_->add(owned);
      } else {
        // Other NodeColumns
        atlas::functionspace::NodeColumns fs(functionSpace_);
        atlas::Field owned = fs.createField<int>(atlas::option::name("owned")
          | atlas::option::levels(1));
        auto ownedView = atlas::array::make_view<int, 2>(owned);
        auto ghostView = atlas::array::make_view<int, 1>(fs.ghost());
        for (atlas::idx_t jnode = 0; jnode < owned.shape(0); ++jnode) {
          ownedView(jnode, 0) = ghostView(jnode) > 0 ? 0 : 1;
        }

        // Add owned points mask
        group.fields_->add(owned);
      }
    }

    // Mask size
    group.gmaskSize_ = 0.0;
    size_t domainSize = 0.0;
    auto ghostView = atlas::array::make_view<int, 1>(functionSpace_.ghost());
    for (atlas::idx_t jnode = 0; jnode < gmask.shape(0); ++jnode) {
      for (atlas::idx_t jlevel = 0; jlevel < gmask.shape(1); ++jlevel) {
        if (ghostView(jnode) == 0) {
          if (maskView(jnode, jlevel) == 1) {
            group.gmaskSize_ += 1.0;
          }
          domainSize++;
        }
      }
    }
    comm_.allReduceInPlace(group.gmaskSize_, eckit::mpi::sum());
    comm_.allReduceInPlace(domainSize, eckit::mpi::sum());
    if (domainSize > 0) {
      group.gmaskSize_ = group.gmaskSize_/static_cast<double>(domainSize);
    }

    // Save group
    groups_.push_back(group);

    // Increment group index
    groupIndex++;
  }

  // Print summary
  this->print(oops::Log::info());
}
// -----------------------------------------------------------------------------
Geometry::Geometry(const Geometry & other) : comm_(other.comm_), halo_(other.halo_),
  grid_(other.grid_), gridType_(other.gridType_), regionalGrid_(other.regionalGrid_),
  partitioner_(other.partitioner_), mesh_(other.mesh_),
  groupIndex_(other.groupIndex_)  {
  // Copy function space
  if (other.functionSpace_.type() == "StructuredColumns") {
    // StructuredColumns
    functionSpace_ = atlas::functionspace::StructuredColumns(other.functionSpace_);
  } else if (other.functionSpace_.type() == "NodeColumns") {
    // NodeColumns
    if (grid_.name().compare(0, 2, std::string{"CS"}) == 0) {
      // CubedSphere
      functionSpace_ = atlas::functionspace::CubedSphereNodeColumns(other.functionSpace_);
    } else {
      // Other NodeColumns
      functionSpace_ = atlas::functionspace::NodeColumns(other.functionSpace_);
    }
  } else if (other.functionSpace_.type() == "PointCloud") {
    ABORT(other.functionSpace_.type() + " function space not supported");
  } else {
    ABORT(other.functionSpace_.type() + " function space not supported yet");
  }

  // Copy groups
  for (size_t groupIndex = 0; groupIndex < other.groups_.size(); ++groupIndex) {
    // Define group
    groupData group;

    // Copy number of levels
    group.levels_ = other.groups_[groupIndex].levels_;

    // Copy corresponding level for 2D variables (first or last)
    group.lev2d_ = other.groups_[groupIndex].lev2d_;

    // Copy vertical coordinate
    group.vert_coord_ = other.groups_[groupIndex].vert_coord_;

    // Copy geometry fields
    group.fields_ = atlas::FieldSet();
    group.fields_->add(other.groups_[groupIndex].fields_["vert_coord"]);
    group.fields_->add(other.groups_[groupIndex].fields_["gmask"]);
    if (other.groups_[groupIndex].fields_.has("owned")) {
      group.fields_->add(other.groups_[groupIndex].fields_["owned"]);
    }

    // Copy mask size
    group.gmaskSize_ = other.groups_[groupIndex].gmaskSize_;

    // Save group
    groups_.push_back(group);
  }
}
// -----------------------------------------------------------------------------
size_t Geometry::levels(const std::string & var) const {
  if (groupIndex_.count(var) == 0) {
    ABORT("Variable " + var + " not found in groupIndex_");
  }
  return groups_[groupIndex_.at(var)].levels_;
}
// -----------------------------------------------------------------------------
size_t Geometry::groupIndex(const std::string & var) const {
  if (groupIndex_.count(var) == 0) {
    ABORT("Variable " + var + " not found in groupIndex_");
  }
  return groupIndex_.at(var);
}
// -----------------------------------------------------------------------------
std::vector<size_t> Geometry::variableSizes(const Variables & vars) const {
  std::vector<size_t> sizes;
  for (const auto & var : vars.variablesList()) {
    sizes.push_back(levels(var));
  }
  return sizes;
}
// -----------------------------------------------------------------------------
void Geometry::latlon(std::vector<double> & lats, std::vector<double> & lons,
                      const bool includeHaloForRealLife) const {
  const auto lonlat = atlas::array::make_view<double, 2>(functionSpace_.lonlat());
  const auto ghost = atlas::array::make_view<int, 1>(functionSpace_.ghost());

  // TODO(Algo): Remove/fix the hack below when GeometryData local KD tree needs
  // to be set up correctly (e.g. when UnstructuredInterpolator is used).
  // For now never include halo in the latlon output because halo points from
  // some atlas grids (e.g. gaussian) can have unrealistic latitudes (e.g. more
  // than 90 degrees) and those latitudes can't be handled by KD trees.
  // Global KD trees created in GeometryData are used for communication and
  // don't need halo information.
  // Local KD trees in GeometryData need halo information but aren't used unless
  // UnstructuredInterpolator is used.
  bool includeHalo = false;
  const size_t npts = functionSpace_.size();
  const size_t nptsReturned = [&]() {
    if (includeHalo && comm_.size() > 1) {
      return npts;
    } else {
      size_t result = 0;
      for (atlas::idx_t i = 0; i < ghost.shape(0); ++i) {
        if (ghost(i) == 0) {
          result++;
        }
      }
      return result;
    }
  }();

  lats.resize(nptsReturned);
  lons.resize(nptsReturned);

  size_t count = 0;
  for (size_t jj = 0; jj < npts; ++jj) {
    // copy owned points, i.e. points with ghost==?
    if (ghost(jj) == 0 || (includeHalo && comm_.size() > 1)) {
      lats[count] = lonlat(jj, 1);
      lons[count] = lonlat(jj, 0);
      if (lons[count] < 0.0) lons[count] += 360.0;
      count++;
    }
  }
  ASSERT(count == nptsReturned);
}
// -----------------------------------------------------------------------------
void Geometry::print(std::ostream & os) const {
  std::string prefix;
  if (os.rdbuf() == oops::Log::info().rdbuf()) {
    prefix = "Info     : ";
  }
  os << prefix <<  "Quench geometry grid:" << std::endl;
  os << prefix << "- name: " << grid_.name() << std::endl;
  os << prefix << "- size: " << grid_.size() << std::endl;
  if (regionalGrid_) {
    os << prefix << "Regional grid detected" << std::endl;
  }
  if (partitioner_) {
    os << prefix << "Partitioner:" << std::endl;
    os << prefix << "- type: " << partitioner_.type() << std::endl;
  }
  os << prefix << "Function space:" << std::endl;
  os << prefix << "- type: " << functionSpace_.type() << std::endl;
  os << prefix << "- halo: " << halo_ << std::endl;
  os << prefix << "Groups: " << std::endl;
  for (size_t groupIndex = 0; groupIndex < groups_.size(); ++groupIndex) {
    os << prefix << "- Group " << groupIndex << ":" << std::endl;
    os << prefix << "  Vertical levels: " << std::endl;
    os << prefix << "  - number: " << levels(groupIndex) << std::endl;
    os << prefix << "  - vert_coord: " << groups_[groupIndex].vert_coord_ << std::endl;
    os << prefix << "  Mask size: " << static_cast<int>(groups_[groupIndex].gmaskSize_*100.0)
       << "%" << std::endl;
  }
}
// -----------------------------------------------------------------------------
void Geometry::readSeaMask(const std::string & maskPath,
                           const size_t & levels,
                           const std::string & lev2d,
                           atlas::Field & gmask) const {
  // Lon/lat sizes
  size_t nlon = 0;
  size_t nlat = 0;

  // NetCDF IDs
  int ncid, retval, nlon_id, nlat_id, lon_id, lat_id, lsm_id;

  if (comm_.rank() == 0) {
    // Open NetCDF file
    if ((retval = nc_open(maskPath.c_str(), NC_NOWRITE, &ncid))) ERR(retval);

    // Get lon/lat sizes
    if ((retval = nc_inq_dimid(ncid, "lon", &nlon_id))) ERR(retval);
    if ((retval = nc_inq_dimid(ncid, "lat", &nlat_id))) ERR(retval);
    if ((retval = nc_inq_dimlen(ncid, nlon_id, &nlon))) ERR(retval);
    if ((retval = nc_inq_dimlen(ncid, nlat_id, &nlat))) ERR(retval);
  }

  // Broadcast lon/lat sizes
  comm_.broadcast(nlon, 0);
  comm_.broadcast(nlat, 0);

  // Coordinates and land-sea mask
  std::vector<double> lon(nlon);
  std::vector<double> lat(nlat);
  std::vector<int> lsm(nlat*nlon);

  if (comm_.rank() == 0) {
    // Get lon/lat
    if ((retval = nc_inq_varid(ncid, "lon", &lon_id))) ERR(retval);
    if ((retval = nc_inq_varid(ncid, "lat", &lat_id))) ERR(retval);
    if ((retval = nc_inq_varid(ncid, "LSMASK", &lsm_id))) ERR(retval);

    // Read data
    float zlon[nlon][1];
    float zlat[nlat][1];
    uint8_t zlsm[nlat][nlon];
    if ((retval = nc_get_var_float(ncid, lon_id, &zlon[0][0]))) ERR(retval);
    if ((retval = nc_get_var_float(ncid, lat_id, &zlat[0][0]))) ERR(retval);
    if ((retval = nc_get_var_ubyte(ncid, lsm_id, &zlsm[0][0]))) ERR(retval);

    // Copy data
    for (size_t ilon = 0; ilon < nlon; ++ilon) {
      lon[ilon] = zlon[ilon][0];
    }
    for (size_t ilat = 0; ilat < nlat; ++ilat) {
      lat[ilat] = zlat[ilat][0];
    }
    for (size_t ilat = 0; ilat < nlat; ++ilat) {
     for (size_t ilon = 0; ilon < nlon; ++ilon) {
        lsm[ilat*nlon+ilon] = static_cast<int>(zlsm[ilat][ilon]);
      }
    }

    // Close file
    if ((retval = nc_close(ncid))) ERR(retval);
  }

  // Broadcast coordinates and land-sea mask
  comm_.broadcast(lon.begin(), lon.end(), 0);
  comm_.broadcast(lat.begin(), lat.end(), 0);
  comm_.broadcast(lsm.begin(), lsm.end(), 0);

  // Build KD-tree
  atlas::Geometry geometry(atlas::util::Earth::radius());
  atlas::util::IndexKDTree2D search(geometry);
  search.reserve(nlat*nlon);
  std::vector<double> lon2d;
  std::vector<double> lat2d;
  std::vector<size_t> payload2d;
  int jnode = 0;
  for (size_t ilat = 0; ilat < nlat; ++ilat) {
    for (size_t ilon = 0; ilon < nlon; ++ilon) {
      lon2d.push_back(lon[ilon]);
      lat2d.push_back(lat[ilat]);
      payload2d.push_back(jnode);
      ++jnode;
    }
  }
  search.build(lon2d, lat2d, payload2d);

  // Ghost points
  atlas::Field ghost = functionSpace_.ghost();
  auto ghostView = atlas::array::make_view<int, 1>(ghost);

  if (functionSpace_.type() == "StructuredColumns") {
    // StructuredColumns
    atlas::functionspace::StructuredColumns fs(functionSpace_);
    auto lonlatView = atlas::array::make_view<double, 2>(fs.xy());
    auto maskView = atlas::array::make_view<int, 2>(gmask);
    for (atlas::idx_t jnode = 0; jnode < fs.xy().shape(0); ++jnode) {
      if (ghostView(jnode) == 0) {
        // Find nearest neighbor
        size_t nn = search.closestPoint(atlas::PointLonLat{lonlatView(jnode, 0),
          lonlatView(jnode, 1)}).payload();

        // Ocean points for all levels
        for (size_t jlevel = 0; jlevel < levels; ++jlevel) {
          if (lsm[nn] == 0) {
             maskView(jnode, jlevel) = 1;
           } else {
             maskView(jnode, jlevel) = 0;
           }
         }

        // Ocean + small islands for:
        // - the first level of 3D fields,
        // - the 2D fields if lev2d = "first"
        if (lsm[nn] == 3) {
          if ((levels > 1) || (lev2d == "first")) {
            maskView(jnode, 0) = 1;
          }
        }
      }
    }
  } else {
    ABORT("Sea mask not supported for " + functionSpace_.type() + " yet");
  }
}
// -----------------------------------------------------------------------------
}  // namespace quench
