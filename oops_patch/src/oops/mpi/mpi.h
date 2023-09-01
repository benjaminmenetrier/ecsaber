/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <Eigen/Dense>


#include <string>
#include <utility>
#include <vector>

#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"

#include "oops/util/Timer.h"

namespace util {
class DateTime;
}  // namespace util

namespace oops {
namespace mpi {

// ------------------------------------------------------------------------------------------------

/// Communicator with all MPI tasks (ie MPI_COMM_WORLD)
const eckit::mpi::Comm & world();

/// Communicator with each MPI task by itself
const eckit::mpi::Comm & myself();

const eckit::mpi::Comm & clone(const eckit::mpi::Comm &);

// ------------------------------------------------------------------------------------------------

void gather(const eckit::mpi::Comm & comm, const std::vector<double> & send,
            std::vector<double> & recv, const size_t root);

// ------------------------------------------------------------------------------------------------
// allGather for eigen vectors
// ------------------------------------------------------------------------------------------------
void allGather(const eckit::mpi::Comm & comm,
               const Eigen::VectorXd &, Eigen::MatrixXd &);

// ------------------------------------------------------------------------------------------------

// The following functions simplify the allGatherv operation on vectors (reducing it to a single
// function call).

// ------------------------------------------------------------------------------------------------

/// \brief Perform the MPI *all gather* operation on a vector of "plain old data".
///
/// This operation gathers data from all tasks and delivers the combined data to all tasks.
///
/// \tparam T must be a type for which there exists a specialization of eckit::mpi::Data::Type.
///
/// \param[in] comm
///   Communicator.
/// \param[inout] x
///   On input, data owned by this task that need to be delivered to all other tasks. On output,
///   combined data received from all tasks (concatenated in the order of increasing task ranks).
template <typename T>
void allGatherv(const eckit::mpi::Comm & comm, std::vector<T> &x) {
    eckit::mpi::Buffer<T> buffer(comm.size());
    comm.allGatherv(x.begin(), x.end(), buffer);
    x = std::move(buffer.buffer);
}

template <typename T>
void allGatherv(const eckit::mpi::Comm & comm, const std::vector<T> & send, std::vector<T> & recv) {
    eckit::mpi::Buffer<T> buffer(comm.size());
    comm.allGatherv(send.begin(), send.end(), buffer);
    recv = std::move(buffer.buffer);
}

// ------------------------------------------------------------------------------------------------

/// \brief Perform the MPI *all gather* operation on a vector of DateTime objects.
///
/// This operation gathers data from all tasks and delivers the combined data to all tasks.
///
/// \param[in] comm
///   Communicator.
/// \param[inout] x
///   On input, data owned by this task that need to be delivered to all other tasks. On output,
///   combined data received from all tasks (concatenated in the order of increasing task ranks).
void allGatherv(const eckit::mpi::Comm & comm, std::vector<util::DateTime> &x);

// ------------------------------------------------------------------------------------------------

/// \brief Perform the MPI *all gather* operation on a vector of DateTime objects.
///
/// This operation gathers data from all tasks and delivers the combined data to all tasks.
///
/// \param[in] comm
///   Communicator.
/// \param[inout] x
///   On input, data owned by this task that need to be delivered to all other tasks. On output,
///   combined data received from all tasks (concatenated in the order of increasing task ranks).
void allGatherv(const eckit::mpi::Comm & comm, std::vector<std::string> &x);

// ------------------------------------------------------------------------------------------------

/// \brief Perform the exclusive scan operation.
///
/// On output, `x` is set to the sum of the values of `x` passed to this function
/// on all ranks lower than the calling rank (and to 0 on rank 0).
void exclusiveScan(const eckit::mpi::Comm &comm, size_t &x);

// ------------------------------------------------------------------------------------------------
// MPI broadcast utilities based on eckit broadcast.

/// \brief broadcast a vector variable via the eckit broadcast
/// \param comm eckit communicator group
/// \param vectorVar vector for broadcasting
/// \param root root rank for broadcasting
template <typename VecType>
void broadcastVector(const eckit::mpi::Comm & comm, std::vector<VecType> & vectorVar,
                     const int root) {
    // eckit broadcast support vectors, but you need to have the vectors identically
    // sized on both sides before doing the broadcast. This routine will broadcast
    // the vector size so the receiving end can resize properly.
    int vecSize;
    if (comm.rank() == root) {
        vecSize = vectorVar.size();
        comm.broadcast(vecSize, root);
        comm.broadcast(vectorVar, root);
    } else {
        comm.broadcast(vecSize, root);
        vectorVar.resize(vecSize);
        comm.broadcast(vectorVar, root);
    }
}

/// @brief broadcast a bool variable via the eckit broadcast
/// @param comm eckit communicator group
/// @param boolVar variable for broadcasting
/// @param root root rank for broadcasting
void broadcastBool(const eckit::mpi::Comm & comm, bool & boolVar, const int root);

/// \brief broadcast a string variable via the eckit broadcast
/// \param comm eckit communicator group
/// \param stringVar string for broadcasting
/// \param root root rank for broadcasting
void broadcastString(const eckit::mpi::Comm & comm, std::string & stringVar, const int root);

// ------------------------------------------------------------------------------------------------
// MPI send/receive utilities based on eckit send/receive.

/// \brief send a string variable via the eckit send
/// \param comm eckit communicator group
/// \param stringVar string for sending
/// \param toRank rank we are sending to
void sendString(const eckit::mpi::Comm & comm, const std::string & stringVar, const int toRank);

/// \brief receive a string variable via the eckit receive
/// \param comm eckit communicator group
/// \param stringVar string for receiving
/// \param fromRank rank we are receiving from
void receiveString(const eckit::mpi::Comm & comm, std::string & stringVar, const int fromRank);

}  // namespace mpi
}  // namespace oops
