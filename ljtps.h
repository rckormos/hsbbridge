#ifndef LJTPS_H
#define LJTPS_H
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <random>
#include <vector>
#include <fstream>
#include <cassert>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <ctime>
#include <omp.h>

using namespace std;

//*INITIAL FRAME GENERATION*//
// helper constructor to make suitable initial configurations in 1D
vector< vector<double> > tileStraight(int particles, int dimensions, double bondlength);
// helper constructor to make suitable initial configurations in 2D; tiles the first and
// second spatial dimensions with particles hexagonally close-packed at equilibrium length
vector< vector<double> > tileHexagonal(int particles, int dimensions, double bondlength);

//*VECTOR ARITHMETIC*//
vector<double> addvectors(vector<double> u, vector<double> v);
vector<double> scalevector(double c, vector<double> v);
vector<double> subvectors(vector<double> u, vector<double> v);
double dotvectors(vector<double> u, vector<double> v);
double squarevector(vector<double> u);

class Traj
{
public:
  // simulation parameters
    double kT; // kB * T
    double dt; // time step
    double epsilon; // potential well depth
    int trajlength; // length of trajectories
    int configdim; // dimension of configuration space
    int numparticles; // number of particles
    vector<double> gamma; // friction coefficients for each particle, units of reciprocal time
    vector<double> sigma; // distance of particles for which the potential is zero;
                          // for radii less than this, the potential is infinite
    bool calcvelocities; // should velocities at each step be computed or not?
    int phasedim; // dimension of phase space per particle = configdim *
                  // (static_cast<int>calcvelocities + 1)

  // the following vectors store copies of the trajectory, force, potential, and interparticle
  // distances; the copy with index 0 is operated upon by dynamical propagation; that with index 1
  // is the copy used for exchange at the end of a move; those with higher indices are used in annealing;
  // hierarchically, the indices proceed as: trajectory index, frame, particle, dimension
    vector< vector< vector< vector<double> > > > traj; // the complete trajectory
    vector< vector< vector< vector<double> > > > force; // the forces (potential derivatives) in each
                                                        // dimension at each frame
    vector< vector<double> > trajpot; // the potential of the system at each frame
    vector< vector< vector< vector<double> > > > distances; // matrix of interparticle distances
                                                            // for each frame of each trajectory
    int numcopies; // sets the number of copies of traj, force, trajpot, and distances (usually 2)

  // constructs helpful for hard sphere bridges
    vector< vector<double> > dr; // the differential change in configuration
    vector< vector<double> > pseudodistances; // the distances calculated when only some
                                              // of the particles have been moved
    vector< vector<int> > lincombs; // a matrix containing the linear combinations
                                    // to transform from Cartesian to Jacobi coordinates
    vector<int> clustersize;  // gives the size of a particle's cluster
    vector<int> cluster;  // denotes which Jacobi cluster a particle is a part of
    double Qrat; // contains the ratio of adjoint probabilities Q for a hard sphere bridge
    double prevQrat; // contrains the ratio from the previous trajectory
    bool wasted; // determines if a hard sphere bridge trajectory was rejected for overlaps

  // enhanced sampling parameters
    int shiftlength; // maximum length of a shifting move
    int anneallength; // maximum length of an annealing move
    int minlength; // minimum bridge length
    int maxlength; // maximum bridge length
    int numbridges; // number of bridge moves performed in an annealing move
    double Tmax; // maximum temperature of annealing moves, relative to kT

  // desired starting and ending states
    int startstate;
    int endstate;

  // rng for noises
    gsl_rng *generator;
    int seed; // will use time for now, could be modified later

  //*CONSTRUCTORS*//
  // default constructor
  Traj();
  // constructor with parameters
  Traj(double kTp, double dtp, double gammap, double epsilonp, double sigmap,
     int trajlengthp, int configdimp, int numparticlesp, bool calcvelocitiesp,
     int shiftlengthp, int anneallengthp, int minlengthp, int maxlengthp,
     int numbridgesp, double Tmaxp, int startstatep, int endstatep);
  // constructor given parameters and the integer index of a seed trajectory
  Traj(double kTp, double dtp, double gammap, double epsilonp, double sigmap,
     int trajlengthp, int configdimp, int numparticlesp, bool calcvelocitiesp,
     int shiftlengthp, int anneallengthp, int minlengthp, int maxlengthp,
     int numbridgesp, double Tmaxp, int startstatep, int endstatep, int trajseed);
  // destructor
  ~Traj();

  //*UTILITY FUNCTIONS*//
  // outputs the trajectory to .xyz format
  void printTrajPos(int iteration, int index, int startframe, int endframe);
  // outputs the force for each particle at each frame to .txt format
  void printTrajForce(int iteration, int index, int startframe, int endframe);
  // prints the current potential of each particle in as many columns as there
  // are relevant spatial dimensions for the spatial derivatives
  void printTrajPot(int iteration, int index, int startframe, int endframe);
  // outputs the satisfaction by a trajectory of a characteristic function at
  // each frame, taking as input the integer code for the characteristic function
  void printTrajhB(int iteration, int instate, int startframe, int endframe);
  // outputs the state of the trajectory as an integer index according to the
  // following: 0 for uncategorized, 1 for C_0^0, 2 for C_0^1, 3 for C_1^4, 4 for
  // C_1^1, 5 for C_1^3, 6 for C_2^1, 7 for C_2^4, 8 for C_2^0, 9 for C_1^0
  void printTrajState(int iteration, int startframe, int endframe);
  // outputs the position versus time of one particle for the whole trajectory
  void printPartPos(int iteration, int index, int particle, int startframe, int endframe);
  // prints a seed trajectory to as many files as there are particles in the system
  // in the format of printPartPos for each file
  void printSeedTraj(int iteration);

  //*DYNAMICS*//
  // helper to calculate all interparticle distances; if "particle" is -1, the
  // entire matrix is calculated; outputs 1 if any pair of particles are within
  // one another's hard sphere radii
  bool calcDistances(int frame, int particle);
  // helper to find the indices of the closest particle pairing, according to
  // the pseudodistances matrix; returns as a vector of two integers
  vector<int> findShortestDist(int frame);
  // helper to calculate the potential and its derivatives with respect to
  // single particle positions
  void updatePotential(int frame);
  // updates the trajectory and initial trajectory information after a TPS move;
  // copyindex denotes the storage trajectory to be updated with the new information
  void resetTraj(int startframe, int endframe, bool success, int copyindex);
  // updates the positions (and possibly velocities) of all the particles
  // according to overdamped Langevin equations; places hard walls at sigma
  // if reflecting is 1
  void propagateDynamics(int startframe, int endframe, bool velocities, bool reflecting);
  // creates a free Brownian bridge between two points
  void propagateBridgeDynamics(int startframe, int endframe, int particle,
    bool velocities, bool updatepot);
  // creates a Brownian bridge whereby no two particles pass within the sum
  // of their sigma values of one another; returns 1 in the rare case that they do
  bool propagateHardSphereDynamics(int startframe, int endframe, bool absorbing,
    bool velocities, bool updatepot);
  // propagates the two-body conditioned dynamics of a vector r given an unconditioned
  // displacement vector drrel, the presence of a mirror plane with vector representation rp,
  // a final configuration rf, a final time tf, and diffusion constant D, as well as a noise
  // vector drrel; if absorbing is 1, the plane is treated as an absorbing, not reflecting,
  // boundary condition
  vector<double> propagateMirrorDynamics(vector<double> r, vector<double> rp,
    vector<double> drrel, vector<double> rf, double tf, double t, double D, bool absorbing);

  //*CHARACTERISTIC FUNCTIONS*//
  // helper function; determines if two particles are within 1.2 rm of each other
  bool near(int frame, int particleA, int particleB);
  // helper function; determines if two particles are more than 1.8 rm away from each other
  bool far(int frame, int particleA, int particleB);
  // helper function; returns the number of nearest neighbors a given particle has
  int numNeighbors(int frame, int particle);
  // helper function; returns the number of particles that are far from a given particle
  int numFar(int frame, int particle);
  // determines whether a trajectory has a given particle with 6 neighbors
  bool inCenter(int frame, int particle);
  // determines whether a trajectory has a particle with no nearest neighbors
  bool scattered(int frame);
  // determines whether all particles are unphysically stacked together
  bool stacked(int frame);
  // helper function; returns a vector containing the number of particles with the index of neighbors
  vector<int> occupationVector(int frame);
  // gives an integer index corresponding to the state the system is in
  int state(int frame);
  // determines whether a trajectory is reactive from one state to another, taking
  // as input the integer state indices as specified by "state"
  bool reactive(int startstate, int endstate);
  // determines whether a segment of trajectory is reactive from one state to another, taking
  // as input the integer state indices as specified by "state"
  bool reactive(int startstate, int endstate, int startframe, int endframe);

  //*PATH ACTION CALCULATIONS*//
  // calculates the dimensionless path action for a section of trajectory
  double S(int startframe, int endframe, int particle, int dimension);
  // calculates the dimensionless path action for a section of trajectory,
  // as if no forces were acting
  double S0(int startframe, int endframe, int particle, int dimension);
  // more efficiently calculates the difference between S and S0
  double deltaS(int startframe, int endframe, int particle, int dimension);
  // calculates the effective action associated with including reflecting boundary
  // conditions between all close particle pairs and excluding the component of grad
  // lnQ perpendicular to the reflecting barrier in relative coordinate space;
  // works for a single step of dynamical propagation, and takes most of the same
  // inputs as "propagateMirrorDynamics"
  double Qratio(vector<double> r, vector<double> rnew, vector<double> rp,
    vector<double> rf, double tf, double t, double D, bool absorbing);

  //*TRANSITION PATH SAMPLING*//
  int moveframe; // returns the frame from which a TPS move began
  int endmoveframe; // only for bridges and annealing; returns the frame
                    // where the move ended
  int choice; // gives an integer index for the type of move employed

  vector< vector<double> > replicaprobs; // useful for statistical analysis of annealing

  // dynamically shoots forward from a randomly selected frame in the trajectory, overwriting
  // the previous states that were farther along
  bool shootForward();
  // dynamically shoots backward from a randomly selected frame in the trajectory, overwriting
  // the previous states farther back in the trajectory
  bool shootReverse();
  // shoots both ways from a randomly selected frame in the trajectory
  bool shootTwoWays();
  // shifts (reptates) forward, making a randomly selected state the new final state and dynamically
  // filling in new states at the start of the trajectory
  bool shiftForward();
  // shifts (reptates) backward, making a randomly selected state the new initial state and dynamically
  // filling in new states at the end of the trajectory
  bool shiftReverse();
  // connect two specified frames by a free Brownian bridge; the index of
  // initTraj to be updated after the move is specified as a parameter
  bool bridge(int startframe, int endframe, int copyindex);
  // connect two randomly selected frames by Brownian bridge
  bool bridge();
  // connect two specified frames by hard sphere Brownian bridge; the index of
  // initTraj to be updated after the move is specified as a parameter
  bool HSbridge(int startframe, int endframe, int copyindex);
  // connect two randomly selected frames by hard sphere Brownian bridge
  bool HSbridge();
  // employs an annealing strategy to generate a trajectory by raising and lowering
  // the temperature through a series of intermediate trajectories
  bool anneal();
  // randomly shoot, shift, anneal, or bridge with a given set of probabilities
  bool moveRandom(double shootprob, double shiftprob, double annealprob);

  // variables helpful for annealing in parallel
  vector<double> threadSvector;
  vector<double> bridgeSvector;

  //*TPS STATISTICS*//
  vector<long long int> successes; // keeps track of the number of successes of each move type
  vector<long long int> attempts; // keeps track of the number of attempts of each move type

  // output to screen success statistics for move types
  void TPSstats();
  // print success statistics for move types to a given file name
  void TPSstats(string filename);
};

#endif
