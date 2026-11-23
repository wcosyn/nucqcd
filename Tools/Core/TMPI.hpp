/*
 * TMPI.hpp
 *
 * Author: Pieter Vancraeyveld (pieter.vancraeyveld@ugent.be)
 *
 * Work started on: February,25 2009
 *
 */

#ifndef TMPI_H
#define TMPI_H
#define STRANGEDEUTERON_MPI

#define NROFRES 9

#if defined(STRANGEDEUTERON_MPI) && !defined(__MAKECINT__)
#include <mpi.h>
#else
typedef int MPI_Op;
typedef int MPI_Datatype;
#endif
#include <fstream>

/*! \brief Coordinates MPI calculations: one object per program (created once, it initializes MPI and
 * finalizes it in its destructor), with the rank, the master/slave roles, output silencing for the
 * slaves, and gathering of results of loops distributed over the processes (GatherResults) */
class TMPI /* FINAL */
{
 public:
  /// results to be gathered in GatherResults
  typedef struct { 
    double value[NROFRES]; ///< the values
    int    save;           ///< does this process hold the values?
  } ValueSave_t;

 public:
  TMPI(int*,char***); ///< initializes MPI with argc and argv; only one TMPI object can be created
  ~TMPI(); ///< finalizes MPI
  
  static bool IsMaster() { return fgRank==kMasterRank; } ///< is this the master process?
  static bool IsSlave() { return fgRank!=kMasterRank; } ///< is this a slave process?
  static int Rank() { return fgRank; } ///< rank of this process
  static int NumberOfProcesses() { return fgNproc; } ///< number of processes
  static int RankOfMaster() { return kMasterRank; } ///< rank of the master process (0)
 
  static std::ostream& Cout(); ///< standard output, /dev/null for the slaves
  static std::ostream& Cerr(); ///< standard error, /dev/null for the slaves
  static void SilenceSlaves(bool=true); ///< redirects cout and cerr of the slaves to /dev/null (false: restores them)

  static void GatherResults(int,double**); ///< collects in the master the results of a loop distributed as i = Rank(), Rank()+NumberOfProcesses(), ...
      
 private:
  TMPI(const TMPI&); ///< not copyable
  TMPI& operator=(const TMPI&); ///< not assignable
  static void CombineResults(void*,void*,int*,MPI_Datatype*); ///< MPI reduction operator used by GatherResults

 private:
  static const int kMasterRank=0; ///< rank of the master process
  static int fgRank; ///< rank of this process
  static int fgNproc; ///< number of processes
  static bool fgMpiInitialized; ///< has the TMPI object been created?
  static std::ofstream fgDevNull; ///< stream to /dev/null
  static MPI_Op MPI_COMBINE_RESULTS; ///< the reduction operator of GatherResults
  static MPI_Datatype datatype; ///< MPI data type of ValueSave_t

}; // class TMPI

#endif // TMPI_H
