/*
 * (C) Copyright 2023 Meteorologisk Institutt
 * 
 */

#pragma once

#include <ostream>
#include <string>
#include <unordered_map>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"

#include "eckit/mpi/Comm.h"

#include "util/ObjectCounter.h"
#include "util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace quench {

// -----------------------------------------------------------------------------
/// Geometry handles geometry for quench model.

class Geometry : public util::Printable,
                 private util::ObjectCounter<Geometry> {
 public:
  static const std::string classname() {return "quench::Geometry";}

  explicit Geometry(const eckit::Configuration &);
  Geometry(const Geometry &);

  const eckit::mpi::Comm & getComm() const {return comm_;}
  const size_t halo() const {return halo_;}
  const atlas::Grid grid() const {return grid_;}
  const std::string gridType() const {return gridType_;}
  const atlas::grid::Partitioner partitioner() const {return partitioner_;}
  const atlas::Mesh mesh() const {return mesh_;}
  const atlas::FunctionSpace & functionSpace() const {return functionSpace_;}
  atlas::FunctionSpace & functionSpace() {return functionSpace_;}
  const atlas::FieldSet & fields() const {return groups_[0].fields_;}
  const atlas::FieldSet & fields(const size_t & groupIndex) const
    {return groups_[groupIndex].fields_;}
  size_t levels(const size_t & groupIndex) const {return groups_[groupIndex].levels_;}
  size_t levels(const std::string & var) const;
  size_t groups() const {return groups_.size();}
  size_t groupIndex(const std::string & var) const;

  size_t variableSize(const std::string &) const;
  size_t maskLevel(const std::string &, const size_t &) const;
  std::vector<size_t> variableSizes(const oops::Variables & vars) const;
  void latlon(std::vector<double> &, std::vector<double> &, const bool) const;
  bool levelsAreTopDown() const {return true;}

 private:
  void print(std::ostream &) const;
  void readSeaMask(const std::string &, const size_t &, const std::string &, atlas::Field &) const;
  const eckit::mpi::Comm & comm_;
  size_t halo_;
  atlas::Grid grid_;
  std::string gridType_;
  bool regionalGrid_;
  bool unstructuredGrid_;
  atlas::grid::Partitioner partitioner_;
  atlas::grid::Distribution distribution_;
  atlas::Mesh mesh_;
  atlas::FunctionSpace functionSpace_;
  std::unordered_map<std::string, size_t> groupIndex_;
  struct groupData {
    size_t levels_;
    std::string lev2d_;
    std::vector<double> vert_coord_;
    atlas::FieldSet fields_;
    double gmaskSize_;
  };
  std::vector<groupData> groups_;
};

// -----------------------------------------------------------------------------

}  // namespace quench
