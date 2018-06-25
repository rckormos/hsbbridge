#include "ljtps.h"

// ADD FUNCTIONALITY FOR VELOCITY TRACKING LATER

//*INITIAL TRAJECTORY GENERATION*//
// helper constructor to make suitable initial configurations in 1D
vector< vector<double> > tileStraight(int particles, int dimensions, double sigma)
{
  double bondlength = pow(2,1/6) * 2.0 * sigma;
  vector< vector<double> > positions(particles, vector<double>(dimensions));
  double xposition = 0.0;
  double xdisplacement = bondlength;
  for(int partcounter = 0; partcounter < particles; partcounter++)
  {
    positions[partcounter][0] = xposition;
    xposition += xdisplacement;
    xdisplacement *= -(partcounter + 2);
    for(int dimcounter = 1; dimcounter < dimensions; dimcounter++)
    {
      positions[partcounter][dimcounter] = 0.0;
    }
  }
  return positions;
}

// helper constructor to make suitable initial configurations in 2D; tiles the first and
// second spatial dimensions with particles hexagonally close-packed at equilibrium length
vector< vector<double> > tileHexagonal(int particles, int dimensions, double sigma)
{
  assert(dimensions > 1);
  double rm = sigma * 2.0 * pow(2.0, 1.0 / 6.0);
  vector< vector<double> > positions(particles, vector<double>(dimensions));
  for(int dimcounter = 0; dimcounter < dimensions; dimcounter++)
    positions[0][dimcounter] = 0.0;
  short direction = 0; // 0 to 5, representing the six directions one can travel
                       // on a hexagonal lattice
  int ring = 1; // viewing the plane as a series of concentric hexagonal rings
                // radiating out from the origin, this is the current ring
  double xposition = 0.0; // the current x position the loop is at
  double yposition = 0.0; // the current y position the loop is at
  double xdisplacement = rm; // how far to go in x to place next particle
  double ydisplacement = 0.0; // how far to go in y
  int stepbound; // determines whether to step $ring$ or $ring - 1$
  int particlecounter = 1; // keeps track of which particle the loop is at
  while(particlecounter < particles)
  {
    // when spiraling out hexagonally, one side must be shorter than the others
    if(direction == 1)
    {
      stepbound = ring - 1;
    }
    else
    {
      stepbound = ring;
    }
    bool flag = 0;
    for(int stepcounter = 0; stepcounter < stepbound && flag == 0; stepcounter++)
    {
      xposition += xdisplacement;
      yposition += ydisplacement;
      positions[particlecounter][0] = xposition;
      positions[particlecounter][1] = yposition;
      // set positions in higher dimensions than the second equal to 0.0
      for(int dimcounter = 2; dimcounter < dimensions; dimcounter++)
      {
        positions[particlecounter][dimcounter] = 0.0;
      }
      particlecounter++;
      if(particlecounter >= particles)
      {
        flag = 1;
      }
    }
    direction = (direction + 1) % 6; // change direction (qualitative)
    switch(direction) // change direction (quantitative)
    {
      case 0:
          xdisplacement = rm;
          ydisplacement = 0.0;
          break;
      case 1:
          xdisplacement = rm * 0.5;
          ydisplacement = rm * 0.86602540378;
          break;
      case 2:
          xdisplacement = -rm * 0.5;
          ydisplacement = rm * 0.86602540378;
          break;
      case 3:
          xdisplacement = -rm;
          ydisplacement = 0.0;
          break;
      case 4:
          xdisplacement = -rm * 0.5;
          ydisplacement = -rm * 0.86602540378;
          break;
      case 5:
          xdisplacement = rm * 0.5;
          ydisplacement = -rm * 0.86602540378;
          break;
    }
    if(direction == 0)
    {
      ring++;
    }
  }
  return positions;
}

//*VECTOR ARITHMETIC*//
vector<double> addvectors(vector<double> u, vector<double> v)
{
  assert(u.size() == v.size());
  vector<double> sum(u.size(), 0.0);
  for(unsigned int i = 0; i < u.size(); i++)
  {
    sum[i] = u[i] + v[i];
  }
  return sum;
}

vector<double> scalevector(double c, vector<double> v)
{
  vector<double> product(v.size(), 0.0);
  for(unsigned int i = 0; i < v.size(); i++)
  {
    product[i] = c * v[i];
  }
  return product;
}

vector<double> subvectors(vector<double> u, vector<double> v)
{
  return addvectors(u,scalevector(-1.0,v));
}

double dotvectors(vector<double> u, vector<double> v)
{
  assert(u.size() == v.size());
  double dotproduct = 0.0;
  for(unsigned int i = 0; i < u.size(); i++)
  {
    dotproduct += u[i] * v[i];
  }
  return dotproduct;
}

double squarevector(vector<double> u)
{
  return dotvectors(u,u);
}

//*CONSTRUCTORS*//
Traj::Traj(): Traj(0.05, 0.001, 1.0, 1.0, 1.0, 5000, 2, 7, 0, 1000, 1000, 500, 1250,
  1000, 1.0, 0, 1, -1){}

Traj::Traj(double kTp, double dtp, double gammap, double epsilonp, double sigmap,
  int trajlengthp, int configdimp, int numparticlesp, bool calcvelocitiesp, int
  shiftlengthp, int anneallengthp, int minlengthp, int maxlengthp, int numbridgesp,
  double Tmaxp, int startstatep, int endstatep): Traj(kTp, dtp, gammap, epsilonp,
    sigmap, trajlengthp, configdimp, numparticlesp, calcvelocitiesp, shiftlengthp,
    anneallengthp, minlengthp, maxlengthp, numbridgesp, Tmaxp, startstatep,
    endstatep, -1){}

Traj::Traj(double kTp, double dtp, double gammap, double epsilonp, double sigmap,
  int trajlengthp, int configdimp, int numparticlesp, bool calcvelocitiesp, int
  shiftlengthp, int anneallengthp, int minlengthp, int maxlengthp, int numbridgesp,
  double Tmaxp, int startstatep, int endstatep, int trajseed)
{
  // define class parameters in terms of function arguments
  kT = kTp;
  dt = dtp;
  epsilon = epsilonp;
  trajlength = trajlengthp;
  configdim = configdimp;
  numparticles = numparticlesp;
  gamma = vector<double>(numparticles,gammap);
  sigma = vector<double>(numparticles,sigmap);
  calcvelocities = calcvelocitiesp;
  phasedim = configdim * (static_cast<int>(calcvelocities) + 1);

  shiftlength = shiftlengthp;
  anneallength = anneallengthp;
  minlength = minlengthp;
  maxlength = maxlengthp;
  numbridges = numbridgesp;
  Tmax = Tmaxp;

  startstate = startstatep;
  endstate = endstatep;

  // initialize replicaprobs vector (for statistical analysis of annealing)
  replicaprobs.resize(numbridges);
  for(int tempcounter = 0; tempcounter < numbridges; tempcounter++)
  {
    replicaprobs[tempcounter].resize(2);
  }

  // rng setup
  generator = gsl_rng_alloc(gsl_rng_mt19937);
  seed = rand()*time(NULL);
  gsl_rng_set(generator, seed);

  // initialize trajectory
  traj.resize(1);
  force.resize(1);
  trajpot.resize(1);
  distances.resize(1);
  // 3D dynamic array that stores phase space coordinates of each particle for
  // each step of the trajectory
  traj[0].resize(trajlength);
  // 3D dynamic array that stores the spatial derivatives of the potential at each
  // particle for each step of the trajectory
  force[0].resize(trajlength);
  // 1D dynamic array that stores the potential for each step of the trajectory
  trajpot[0].resize(trajlength);
  // 3D dynamic array that stores the distances between each pair of particles
  // for each step of the trajectory
  distances[0].resize(trajlength);
  // initialize substituent vectors
  for (int framecounter = 0; framecounter < trajlength; framecounter++)
  {
    traj[0][framecounter].resize(numparticles);
    force[0][framecounter].resize(numparticles);
    distances[0][framecounter].resize(numparticles);
    for (int partcounter = 0; partcounter < numparticles; partcounter++)
    {
      traj[0][framecounter][partcounter].resize(phasedim,0.0);
      force[0][framecounter][partcounter].resize(configdim,0.0);
      distances[0][framecounter][partcounter].resize(numparticles,0.0);
    }
  }

  // determine if a seed is intended to be used
  if (trajseed >= 0)
  {
    // read in trajectory from source files (one for each particle)
    string striteration = to_string(trajseed);
    for (int partcounter = 0; partcounter < numparticles; partcounter++)
    {
      string filename = "seed_" + striteration + "_particle_" +
        to_string(partcounter) + ".txt";
      ifstream sourcetrajfile;
      sourcetrajfile.open(filename);
      int rowcounter = -1; // starts at -1 so that it can be incremented to 0 for first row
      int columncounter = 0;
      double readnum;
      while (sourcetrajfile >> readnum)
      {
        if (columncounter == 0)
        {
          rowcounter++;
          columncounter = (columncounter + 1)%(configdim + 1);
        }
        else
        {
          traj[0][rowcounter][partcounter][columncounter - 1] = readnum;
          columncounter = (columncounter + 1)%(configdim + 1);
        }
      }
    }
    // calculate force and potential for seed trajectory
    for(int framecounter = 0; framecounter < trajlength; framecounter++)
    {
      updatePotential(framecounter);
    }
    // output relevant .xyz and .txt files for seed trajectory
    printTrajPos(0, 0, 0, trajlength - 1);
    printPartPos(0,0,0,0,1000);
    printPartPos(0,0,0,1,1000);
    printTrajState(0, 0, trajlength - 1);
    // determine if seed trajectory is reactive given the desired start and end states
    if (!reactive(startstate, endstate))
    {
      trajseed = -1;
    }
  }

  // generate an adequate seed trajectory using MD
  if (trajseed == -1)
  {
    cout << "No valid seed trajectory given.  Using MD to generate one." << '\n';
    if (configdim == 1)
    {
      traj[0][0] = tileStraight(numparticles, configdim, sigma[0]);
    }
    else
    {
      traj[0][0] = tileHexagonal(numparticles, configdim, sigma[0]);
    }
    updatePotential(0);

    // double trajectory length for performing MD and scanning segments
    traj[0].resize(2*trajlength);
    force[0].resize(2*trajlength);
    trajpot[0].resize(2*trajlength);
    distances[0].resize(2*trajlength);
    for (int framecounter = trajlength; framecounter < 2*trajlength; framecounter++)
    {
      traj[0][framecounter].resize(numparticles);
      force[0][framecounter].resize(numparticles);
      distances[0][framecounter].resize(numparticles);
      for (int partcounter = 0; partcounter < numparticles; partcounter++)
      {
        traj[0][framecounter][partcounter].resize(phasedim);
        force[0][framecounter][partcounter].resize(configdim);
        distances[0][framecounter][partcounter].resize(numparticles);
      }
    }

    // increase the temperature to facilitate more efficient seed trajectory searching
    double kTorig = kT;
    kT = 1.5 * kT;

    // run MD until a reactive segment of length trajlength is found
    propagateDynamics(0, 2*trajlength - 1, 0, 1);
    bool mdreactive = reactive(startstate, endstate);
    int mdstretchcounter = 0;
    while (!mdreactive)
    {
      mdstretchcounter++;
      // loop through trajectory segments of length trajlength
      for (int framecounter = 0; framecounter < trajlength
        && mdreactive == 0; framecounter++)
      {
        // check to see if a segment of trajectory might have another pair of
        // particles switch from center and edge of a C_0 configuration
        if (state(framecounter) == 2 && state(framecounter + trajlength) == 2)
        {
          int starttag = 0;
          int endtag = 0;
          for (int partcounter = 0; partcounter < numparticles; partcounter++)
          {
            if (numNeighbors(framecounter, partcounter) == 6)
            {
              starttag = partcounter;
            }
            if (numNeighbors(framecounter + trajlength, partcounter) == 6)
            {
              endtag = partcounter;
            }
          }
          // make the particle with 6 neighbors at framecounter the tagged particle
          for (int innercounter = 0; innercounter < 2*trajlength - 1
            && starttag - endtag != 0; innercounter++)
          {
            traj[0][innercounter][starttag].swap(traj[0][innercounter][0]);
            force[0][innercounter][starttag].swap(force[0][innercounter][0]);
            distances[0][innercounter][starttag].swap(distances[0][innercounter][0]);
          }
        }
        // determine if the current trajectory segment is reactive
        if (reactive(startstate, endstate, framecounter, framecounter + trajlength - 1))
        {
          for (int innercounter = framecounter; innercounter <
            framecounter + trajlength; innercounter++)
          {
            traj[0][innercounter - framecounter] = traj[0][innercounter];
            force[0][innercounter - framecounter] = force[0][innercounter];
            trajpot[0][innercounter - framecounter] = trajpot[0][innercounter];
            distances[0][innercounter - framecounter] = distances[0][innercounter];
          }
          mdreactive = 1;
        }
      }
      // if no segment is reactive, generate a new MD trajectory of length 2 * trajlength
      if (!mdreactive)
      {
        traj[0][0].swap(traj[0][2*trajlength - 1]);
        force[0][0].swap(force[0][2*trajlength - 1]);
        swap(trajpot[0][0], trajpot[0][2*trajlength - 1]);
        distances[0][0].swap(distances[0][2*trajlength - 1]);
        propagateDynamics(0, 2*trajlength - 1, 0, 1);
      }
    }
    // once a seed trajectory is found, reset the temperature and trajectory length
    kT = kTorig;
    traj[0].resize(trajlength);
    force[0].resize(trajlength);
    trajpot[0].resize(trajlength);
    distances[0].resize(trajlength);
    cout << "Seed trajectory found." << '\n';
    // ensure that no seed trajectories in the CWD are overwritten
    int iteration = 0;
    string filename = "seed_" + to_string(iteration) + "_particle_0.txt";
    ifstream ifile;
    ifile.open(filename);
    while ((bool)ifile)
    {
      iteration++;
      filename = "seed_" + to_string(iteration) + "_particle_0.txt";
      ifile.open(filename);
    }
    // print .xyz trajectory, .txt potential, and seed trajectory .txt files
    printTrajPos(0, 0, 0, trajlength - 1);
    printTrajPot(0, 0, 0, trajlength - 1);
    printSeedTraj(iteration);
  }

  // initialize constructs relevant to hard sphere propagation
  cluster = vector<int>(numparticles,0);
  clustersize = vector<int>(numparticles,0);
  prevQrat = 0.0;
  dr = vector< vector<double> >(numparticles,vector<double>(configdim,0.0));

  // specifies the number of copies (not counting the main trajectory) needed
  // usually two, one for standard moves and one extra for substituent bridges
  numcopies = 10;
  for(int threadcounter = 1; threadcounter <= numcopies; threadcounter++)
  {
    traj.push_back(traj[0]);
    force.push_back(force[0]);
    trajpot.push_back(trajpot[0]);
    distances.push_back(distances[0]);
  }
  threadSvector.resize(numcopies);

  // move success vector setup
  successes = {0, 0, 0, 0, 0, 0, 0};
  attempts = {0, 0, 0, 0, 0, 0, 0};
}

Traj::~Traj()
{
  gsl_rng_free(generator);  // freeing the generator prevents huge memory leaks
                            // whenever multiple traj classes are constructed
}

//*UTILITY FUNCTIONS*//
// outputs the trajectory to .xyz format
void Traj::printTrajPos(int iteration, int index, int startframe, int endframe)
{
  string title = "kT_" + to_string(kT) + "_dt_"
    + to_string(dt) + "_gam_" + to_string(gamma[0]) + "_epsilon_" +
    to_string(epsilon) + "_sigma_" + to_string(sigma[0]) +  "_len_" +
    to_string(trajlength) + "_part_" + to_string(numparticles) + "_dim_"
    + to_string(configdim) + "_" + to_string(iteration) + ".xyz";
  int step = 1;
  /*
  if ((endframe - startframe) > 1500)
  {
    step = (endframe - startframe) / 1500;
  }
  */
  ofstream positionoutput;
  positionoutput.open(title);
  for (int framecounter = startframe; framecounter <= endframe; framecounter++)
  {
    if (framecounter % step == 0)
    {
      positionoutput << numparticles << '\n';
      positionoutput << "frame " << to_string(framecounter) << '\n';
      for (int partcounter = 0; partcounter < numparticles; partcounter++)
      {
        if (partcounter == 0)
        {
          // positionoutput << "O" << " ";
          positionoutput << "H" << " ";
        }
        else
        {
          // positionoutput << "He" << " ";
          positionoutput << "H" << " ";
        }
        if (configdim == 1)
        {
          positionoutput << traj[index][framecounter][partcounter][0] << " " << 0 <<
            " " << 0 << '\n';
        }
        else if (configdim == 2)
        {
          positionoutput << traj[index][framecounter][partcounter][0] << " " <<
            traj[index][framecounter][partcounter][1] << " " << 0 << '\n';
        }
        else
        {
          positionoutput << traj[index][framecounter][partcounter][0] << " " <<
            traj[index][framecounter][partcounter][1] << " " <<
            traj[index][framecounter][partcounter][2] << '\n';
        }
      }
    }
  }
  positionoutput.close();
}

// outputs the force for each particle at each frame to .xyz format
void Traj::printTrajForce(int iteration, int index, int startframe, int endframe)
{
  string title = "kT_" + to_string(kT) + "_dt_"
    + to_string(dt) + "_gam_" + to_string(gamma[0]) + "_epsilon_" +
    to_string(epsilon) + "_sigma_" + to_string(sigma[0]) +  "_len_" +
    to_string(trajlength) + "_part_" + to_string(numparticles) + "_dim_"
    + to_string(configdim) + "_force_" + to_string(iteration) + ".txt";
  int step = 1;
  if ((endframe - startframe) > 3000)
  {
    step = (endframe - startframe) / 3000;
  }
  ofstream forceoutput;
  forceoutput.open(title);
  for (int framecounter = startframe; framecounter <= endframe; framecounter++)
  {
    if (framecounter % step == 0)
    {
      forceoutput << 1 << '\n';
      forceoutput << "frame " << to_string(framecounter) << '\n';
      forceoutput << "He" << '\t';
      for (int partcounter = 0; partcounter < numparticles; partcounter++)
      {
        for (int dimcounter = 0; dimcounter < configdim; dimcounter++)
        {
          forceoutput << force[index][framecounter][partcounter][dimcounter] << '\t';
        }
        if (configdim < 3)
        {
          for (int extracounter = configdim; extracounter < 3; extracounter++)
          {
            forceoutput << 0.0 << '\t';
          }
        }
        forceoutput << '\n';
      }
    }
  }
  forceoutput.close();
}

// outputs the potential at each frame to .txt format
void Traj::printTrajPot(int iteration, int index, int startframe, int endframe)
{
  string title = "kT_" + to_string(kT) + "_dt_"
    + to_string(dt) + "_gam_" + to_string(gamma[0]) + "_epsilon_" +
    to_string(epsilon) + "_sigma_" + to_string(sigma[0]) +  "_len_" +
    to_string(trajlength) + "_part_" + to_string(numparticles) + "_dim_"
    + to_string(configdim) + "_pots_" + to_string(iteration) + ".xyz";
  int step = 1;
  if ((endframe - startframe) > 3000)
  {
    step = (endframe - startframe) / 3000;
  }
  ofstream potentialoutput;
  potentialoutput.open(title);
  for (int framecounter = startframe; framecounter <= endframe; framecounter++)
  {
    if (framecounter % step == 0)
    {
      potentialoutput << framecounter*dt << '\t' << trajpot[index][framecounter];
      potentialoutput << '\n';
    }
  }
  potentialoutput.close();
}

// outputs the satisfaction by a trajectory of a characteristic function at
// each frame, taking as input the integer code for the characteristic function
void Traj::printTrajhB(int iteration, int instate, int startframe, int endframe)
{
  string title = "kT_" + to_string(kT) + "_dt_"
    + to_string(dt) + "_gam_" + to_string(gamma[0]) + "_epsilon_" +
    to_string(epsilon) + "_sigma_" + to_string(sigma[0]) +  "_len_" +
    to_string(trajlength) + "_part_" + to_string(numparticles) + "_dim_"
    + to_string(configdim) + "_pots_" + to_string(iteration) + ".xyz";
  int step = 1;
  if ((endframe - startframe) > 3000)
  {
    step = (endframe - startframe) / 3000;
  }
  ofstream potentialoutput;
  potentialoutput.open(title);
  for (int framecounter = startframe; framecounter <= endframe; framecounter++)
  {
    if (framecounter % step == 0)
    {
      potentialoutput << framecounter*dt << '\t' << (state(framecounter) == instate);
      potentialoutput << '\n';
    }
  }
  potentialoutput.close();
}

// outputs the state of the trajectory as an integer index according to the
// following: 0 for uncategorized, 1 for C_0^0, 2 for C_0^1, 3 for C_1^4, 4 for
// C_1^1, 5 for C_1^3, 6 for C_2^1, 7 for C_2^4, 8 for C_2^0, 9 for C_1^0
void Traj::printTrajState(int iteration, int startframe, int endframe)
{
  string title = "path_kT_" + to_string(kT) + "_dt_"
    + to_string(dt) + "_gam_" + to_string(gamma[0]) + "_epsilon_" +
    to_string(epsilon) + "_sigma_" + to_string(sigma[0]) +  "_len_" +
    to_string(trajlength) + "_part_" + to_string(numparticles) + "_dim_"
    + to_string(configdim) + "_states_" + to_string(iteration) + ".txt";
  int step = 1;
  if ((endframe - startframe) > 1500)
  {
    step = (endframe - startframe) / 1500;
  }
  ofstream stateoutput;
  stateoutput.open(title);
  for (int framecounter = startframe; framecounter <= endframe; framecounter++)
  {
    if (framecounter % step == 0)
    {
      stateoutput << " ";
      if (framecounter*dt < 10.0)
      {
        stateoutput << " ";
      }
      stateoutput << fixed << setprecision(4) << framecounter*dt << "        "
        << state(framecounter);
      stateoutput << '\n';
    }
  }
  stateoutput.close();
}

// outputs the position versus time of one particle for the whole trajectory
void Traj::printPartPos(int iteration, int index, int particle, int startframe, int endframe)
{
  string title = "kT_" + to_string(kT) + "_dt_"
    + to_string(dt) + "_gam_" + to_string(gamma[0]) + "_epsilon_" +
    to_string(epsilon) + "_sigma_" + to_string(sigma[0]) +  "_len_" +
    to_string(trajlength) + "_numpart_" + to_string(numparticles) + "_dim_"
    + to_string(configdim) + "_part_" + to_string(particle) + "_" +
    to_string(iteration) + ".txt";
  ofstream positionoutput;
  positionoutput.open(title);
  for (int framecounter = startframe; framecounter <= endframe; framecounter++)
  {
    positionoutput << dt*framecounter << '\t';
    for (int dimcounter = 0; dimcounter < configdim; dimcounter++)
    {
      positionoutput << traj[index][framecounter][particle][dimcounter] << '\t';
    }
    positionoutput << '\n';
  }
  positionoutput.close();
}

// prints a seed trajectory to as many files as there are particles in the system
// in the format of printPartPos for each file
void Traj::printSeedTraj(int iteration)
{
  for (int partcounter = 0; partcounter < numparticles; partcounter++)
  {
    string title = "seed_" + to_string(iteration) + "_particle_" + to_string(partcounter)
      + ".txt";
    ofstream positionoutput;
    positionoutput.open(title);
    for (int framecounter = 0; framecounter < trajlength; framecounter++)
    {
      positionoutput << dt*framecounter << '\t';
      for (int dimcounter = 0; dimcounter < configdim; dimcounter++)
      {
        positionoutput << traj[0][framecounter][partcounter][dimcounter] << '\t';
      }
      positionoutput << '\n';
    }
    positionoutput.close();
  }
}

//*DYNAMICS*//

// helper to calculate all interparticle distances; if "particle" is -1, the
// entire matrix is calculated; outputs 1 if any pair of particles are within
// one another's hard sphere radii
bool Traj::calcDistances(int frame, int particle)
{
  bool flag = 0;
  double distance = 0.0;
  for(int partcounter = 0; partcounter < numparticles - 1; partcounter++)
  {
    if(particle == -1 || particle == partcounter)
    {
      for(int intcounter = partcounter + 1; intcounter < numparticles; intcounter++)
      {
        distance = sqrt(squarevector(subvectors(traj[0][frame][partcounter],
          traj[0][frame][intcounter])));
        if(distance < sigma[partcounter] + sigma[intcounter])
        {
          flag = 1;
        }
        if(particle == -1)
        {
          distances[0][frame][partcounter][intcounter] = distance;
          distances[0][frame][intcounter][partcounter] = distance;
        }
        // the case of some particle distances but not others being calculates
        // only occurs in hard sphere propagation, where pseudodistances is referenced
        else
        {
          pseudodistances[partcounter][intcounter] = distance;
          pseudodistances[intcounter][partcounter] = distance;
        }
      }
    }
  }

  return flag;
}

// helper to find the indices of the closest particle pairing, according to
// the pseudodistances matrix; returns as a vector of two integers
vector<int> Traj::findShortestDist(int frame)
{
  vector<int> tuple(2,-1);
  double shortest = 1.7976931348623158e+308; // largest double; any distance is shorter
  for(int partcounter = 0; partcounter < numparticles - 1; partcounter++)
  {
    for(int intcounter = partcounter + 1; intcounter < numparticles; intcounter++)
    {
      if(pseudodistances[partcounter][intcounter] < shortest && (cluster[partcounter]
        != cluster[intcounter]))
      {
        shortest = pseudodistances[partcounter][intcounter];
        tuple[0] = partcounter;
        tuple[1] = intcounter;
      }
    }
  }
  return tuple;
}

// helper to calculate the potential and its derivatives with respect to
// single particle positions
void Traj::updatePotential(int frame)
// Pairwise Lennard-Jones potential
{
  calcDistances(frame,-1);
  double framePot = 0.0;
  vector< vector<double> > frameForce(numparticles, vector<double>(configdim, 0.0));
  for(int partcounter = 0; partcounter < numparticles - 1; partcounter++)
  {
    for(int intcounter = partcounter + 1; intcounter < numparticles; intcounter++)
    {
      vector<double> x_i = subvectors(traj[0][frame][partcounter],
        traj[0][frame][intcounter]);
      double r_square = pow(distances[0][frame][partcounter][intcounter],2);
      double sigsquare = pow(sigma[partcounter] + sigma[intcounter],2);
      double sixthpower = pow(sigsquare / r_square,3);
      framePot += 4 * epsilon * sixthpower * (sixthpower - 1);
      double dVdrOverr = 24 * epsilon * sixthpower / r_square * (1 - 2 * sixthpower);
      for(int dimcounter = 0; dimcounter < configdim; dimcounter++)
      {
        frameForce[partcounter][dimcounter] -= dVdrOverr * x_i[dimcounter];
        frameForce[intcounter][dimcounter] += dVdrOverr * x_i[dimcounter];
      }
    }
  }
  trajpot[0][frame] = framePot;
  force[0][frame] = frameForce;
}

// updates the trajectory and initial trajectory information after a TPS move;
// copyindex denotes the storage trajectory to be updated with the new information
void Traj::resetTraj(int startframe, int endframe, bool success, int copyindex)
{
  // if unsuccessful, swap out the failed trajectory of index 0 with the previous
  // trajectory, stored with index copyindex
  if (!success)
  {
    for(int framecounter = startframe; framecounter <= endframe; framecounter++)
    {
      traj[0][framecounter].swap(traj[copyindex][framecounter]);
      force[0][framecounter].swap(force[copyindex][framecounter]);
      swap(trajpot[0][framecounter], trajpot[copyindex][framecounter]);
      distances[0][framecounter].swap(distances[copyindex][framecounter]);
    }
  }

  // update copies of traj, force, trajpot, and distances
  for(int framecounter = startframe; framecounter <= endframe; framecounter++)
  {
    for(int partcounter = 0; partcounter < numparticles; partcounter++)
    {
      for(int dimcounter = 0; dimcounter < configdim; dimcounter++)
      {
        traj[copyindex][framecounter][partcounter][dimcounter] =
          traj[0][framecounter][partcounter][dimcounter];
        force[copyindex][framecounter][partcounter][dimcounter] =
          force[0][framecounter][partcounter][dimcounter];
        distances[copyindex][framecounter][partcounter][dimcounter] =
          distances[0][framecounter][partcounter][dimcounter];
      }
      for(int intcounter = configdim; intcounter < numparticles; intcounter++)
      {
        distances[copyindex][framecounter][partcounter][intcounter] =
          distances[0][framecounter][partcounter][intcounter];
      }
    }
    trajpot[copyindex][framecounter] = trajpot[0][framecounter];
  }
}

// updates the positions (and possibly velocities) of all the particles
// according to overdamped Langevin equations; places hard walls at sigma
// if reflecting is 1
void Traj::propagateDynamics(int startframe, int endframe, bool velocities, bool reflecting)
// Brownian (overdamped Langevin) dynamics
// Note: the terminology "startframe" means the first frame shot from, not the
// first frame chronologically.  This distinction is irrelevant for forward shooting.
{
  if(startframe < endframe)
  {
    for(int framecounter = startframe; framecounter < endframe; framecounter++)
    {
      for(int partcounter = 0; partcounter < numparticles; partcounter++)
      {
        for(int dimcounter = 0; dimcounter < configdim; dimcounter++)
        {
          double detpart = force[0][framecounter][partcounter][dimcounter]
            / gamma[partcounter];
          double probpart = sqrt(2 * kT / gamma[partcounter] / dt)
            * gsl_ran_gaussian(generator, 1.0);
          traj[0][framecounter + 1][partcounter][dimcounter] =
            traj[0][framecounter][partcounter][dimcounter] + (detpart + probpart) * dt;
        }
      }

      // enforces event-driven elastic collisions
      if(reflecting)
      {
        calcDistances(framecounter + 1, -1);
        for(int partcounter = 0; partcounter < numparticles - 1; partcounter++)
        {
          for(int intcounter = partcounter + 1; intcounter < numparticles; intcounter++)
          {
            if(distances[0][framecounter + 1][partcounter][intcounter] < sigma[partcounter]
              + sigma[intcounter])
            {
              vector<double> r = subvectors(traj[0][framecounter][partcounter],
                traj[0][framecounter][intcounter]);
              double rnorm = sqrt(squarevector(r));
              vector<double> rprime = subvectors(traj[0][framecounter + 1][partcounter],
                traj[0][framecounter + 1][intcounter]);
              double rprimeproj = dotvectors(r,rprime) / rnorm;
              vector<double> dr = scalevector((rprimeproj - sigma[partcounter]
                - sigma[intcounter]) / rnorm,r);
              traj[0][framecounter + 1][partcounter] = subvectors(traj[0][framecounter + 1]
                [partcounter], dr);
              traj[0][framecounter + 1][intcounter] = addvectors(traj[0][framecounter + 1]
                [intcounter], dr);
            }
          }
        }
      }

      updatePotential(framecounter + 1);
    }
  }
  else if (startframe > endframe)
  {
    for (int framecounter = startframe; framecounter > endframe; framecounter--)
    {
      for (int partcounter = 0; partcounter < numparticles; partcounter++)
      {
        for(int dimcounter = 0; dimcounter < configdim; dimcounter++)
        {
          double detpart = force[0][framecounter][partcounter][dimcounter]
            / gamma[partcounter];
          double probpart = sqrt(2*kT/dt/gamma[partcounter])
            * gsl_ran_gaussian(generator, 1.0);
          traj[0][framecounter - 1][partcounter][dimcounter] =
          traj[0][framecounter][partcounter][dimcounter] + (detpart + probpart)*dt;
        }
      }

      // enforces event-driven elastic collisions
      if(reflecting)
      {
        calcDistances(framecounter - 1, -1);
        for(int partcounter = 0; partcounter < numparticles - 1; partcounter++)
        {
          for(int intcounter = partcounter + 1; intcounter < numparticles; intcounter++)
          {
            if(distances[0][framecounter - 1][partcounter][intcounter] < sigma[partcounter]
              + sigma[intcounter])
            {
              vector<double> r = subvectors(traj[0][framecounter][partcounter],
                traj[0][framecounter][intcounter]);
              double rnorm = sqrt(squarevector(r));
              vector<double> rprime = subvectors(traj[0][framecounter - 1][partcounter],
                traj[0][framecounter - 1][intcounter]);
              double rprimeproj = dotvectors(r,rprime) / rnorm;
              vector<double> dr = scalevector((rprimeproj - sigma[partcounter]
                - sigma[intcounter]) / rnorm,r);
              traj[0][framecounter - 1][partcounter] = subvectors(traj[0][framecounter - 1]
                [partcounter], dr);
              traj[0][framecounter - 1][intcounter] = addvectors(traj[0][framecounter - 1]
                [intcounter], dr);
            }
          }
        }
      }

      updatePotential(framecounter - 1);
    }
  }
}

// creates a free Brownian bridge between two points
void Traj::propagateBridgeDynamics(int startframe, int endframe, int particle,
  bool velocities, bool updatepot)
// Brownian (overdamped Langevin) dynamics with fixed endpoints
// note: endframe and startframe should remain unchanged by the dynamics
{
  vector< vector<double> > R(numparticles, vector<double>(configdim, 0.0));
  double T = (endframe - startframe) * dt;
  for(int framecounter = startframe + 1; framecounter < endframe; framecounter++)
  {
    double t = (framecounter - startframe) * dt;
    double Rstd = sqrt(2 * kT / gamma[particle] * (1.0 / (T - t) - 1.0 / (T - t + dt)));
    for(int dimcounter = 0; dimcounter < configdim; dimcounter++)
    {
      R[particle][dimcounter] += Rstd*gsl_ran_gaussian(generator, 1.0);
      traj[0][framecounter][particle][dimcounter] =
            traj[0][startframe][particle][dimcounter] +
            (t / T) * (traj[0][endframe][particle][dimcounter]
            - traj[0][startframe][particle][dimcounter]) +
            (T - t) * R[particle][dimcounter];
    }
    if(updatepot)
    {
      updatePotential(framecounter);
    }
  }
}

// creates a Brownian bridge whereby no two particles pass within the sum
// of their sigma values of one another; returns 1 in the rare case that they do
bool Traj::propagateHardSphereDynamics(int startframe, int endframe, bool absorbing,
  bool velocities, bool updatepot)
// Brownian (overdamped Langevin) dynamics with fixed endpoints and the constraint
// that particles should not pass within a distance of sigmaA + sigmaB of one another
// note: currently only works when all gamma are equal, since the expansion coefficients
// are being treated as equal in the Jacobi vectors
{
  double tf = endframe * dt;
  vector< vector<double> > finalframe = traj[0][endframe];
  vector<double> COMf(configdim, 0.0);
  for(int partcounter = 0; partcounter < numparticles; partcounter++)
  {
    COMf = addvectors(COMf,scalevector(1.0 / numparticles,
      traj[0][endframe][partcounter]));
  }
  Qrat = 1.0;
  bool overlapping = 0;
  for(int frame = startframe; frame < endframe && !overlapping; frame++)
  {
    double t = frame * dt;
    pseudodistances = distances[0][frame];
    // prepare random noises
    for(int partcounter = 0; partcounter < numparticles; partcounter++)
    {
      for(int dimcounter = 0; dimcounter < configdim; dimcounter++)
      {
        dr[partcounter][dimcounter] = gsl_ran_gaussian(generator,
          sqrt(2.0 * kT / gamma[partcounter] * dt));
      }
    }
    // propagate COM Jacobi coordinate
    vector<double> COMnoise(configdim,0.0);
    vector<double> COM(configdim,0.0);
    for(int partcounter = 0; partcounter < numparticles; partcounter++)
    {
      COMnoise = addvectors(COMnoise,scalevector(1.0 / numparticles,dr[partcounter]));
      COM = addvectors(COM,scalevector(1.0 / numparticles,
        traj[0][frame][partcounter]));
    }
    vector<double> COMforce = scalevector(dt / (tf - t),subvectors(COMf,COM));
    for(int partcounter = 0; partcounter < numparticles; partcounter++)
    {
      if(frame != endframe - 1)
      {
        traj[0][frame + 1][partcounter] = addvectors(traj[0][frame][partcounter],
          addvectors(COMnoise,COMforce));
      }
      else
      {
        traj[0][frame + 1][partcounter] = traj[0][frame][partcounter];
      }
    }
    // propagate non-COM Jacobi coordinates
    lincombs.resize(0);
    for(int partcounter = 0; partcounter < numparticles; partcounter++)
    {
      cluster[partcounter] = partcounter + 1;
      clustersize[partcounter] = 1;
    }
    // ...by looping through the linearly combined vectors in the Jacobi coordinate system
    for(int vectorcounter = 0; vectorcounter < numparticles - 1; vectorcounter++)
    {
      vector<int> closestpair = findShortestDist(frame + 1);
      int partnerone = closestpair[0];
      int partnertwo = closestpair[1];
      lincombs.push_back(vector<int>(numparticles,0));
      // determine the coefficients of the current Jacobi vector
      for(int partcounter = 0; partcounter < numparticles; partcounter++)
      {
        if(cluster[partcounter] == cluster[partnerone])
        {
          lincombs[vectorcounter][partcounter] = - clustersize[partcounter];
          if(partcounter != partnerone)
          {
            clustersize[partcounter] += clustersize[partnertwo];
          }
        }
        if(cluster[partcounter] == cluster[partnertwo])
        {
          lincombs[vectorcounter][partcounter] = clustersize[partcounter];
          if(partcounter != partnertwo)
          {
            clustersize[partcounter] += clustersize[partnerone];
            cluster[partcounter] = cluster[partnerone];
          }
        }
      }
      int newsize = clustersize[partnerone] + clustersize[partnertwo];
      clustersize[partnerone] = newsize;
      clustersize[partnertwo] = newsize;
      cluster[partnertwo] = cluster[partnerone];
      // determine the inputs to propagateMirrorDynamics
      vector<double> r(configdim,0.0);
      vector<double> drrel(configdim,0.0);
      vector<double> rf(configdim,0.0);
      double D = 0;
      for(int partcounter = 0; partcounter < numparticles; partcounter++)
      {
        if (lincombs[vectorcounter][partcounter] != 0)
        {
          r = addvectors(r,scalevector(1.0 / lincombs[vectorcounter][partcounter],
            traj[0][frame + 1][partcounter]));
          drrel = addvectors(drrel,scalevector(1.0 / lincombs[vectorcounter][partcounter],
            dr[partcounter]));
          rf = addvectors(rf,scalevector(1.0 / lincombs[vectorcounter][partcounter],
            finalframe[partcounter]));
          D += fabs(1.0 / lincombs[vectorcounter][partcounter]) * kT / gamma[partcounter];
        }
      }
      vector<double> deltar = subvectors(traj[0][frame + 1][partnertwo],
        traj[0][frame + 1][partnerone]);
      double deltarnorm = sqrt(squarevector(deltar));
      vector<double> rp = scalevector(dotvectors(r,deltar) / pow(deltarnorm,2)
        + (sigma[partnerone] + sigma[partnertwo]) / deltarnorm - 1.0, deltar);
      // propagate the dynamics and calculate the ratio of adjoint probabilities Q
      drrel = propagateMirrorDynamics(r, drrel, rp, rf, tf, t, D, absorbing);
      Qrat *= Qratio(r, addvectors(r,drrel), rp, rf, tf, t, D, absorbing);
      for(int partcounter = 0; partcounter < numparticles; partcounter++)
      {
        if(lincombs[vectorcounter][partcounter] < 0 && frame != endframe - 1)
        {
          double factor = -1.0 - 1.0 * lincombs[vectorcounter][partcounter]
            / clustersize[partcounter];
          traj[0][frame + 1][partcounter] = addvectors(traj[0][frame + 1][partcounter],
            scalevector(factor,drrel));
        }
        if(lincombs[vectorcounter][partcounter] > 0 && frame != endframe - 1)
        {
          double factor = 1.0 - 1.0 * lincombs[vectorcounter][partcounter]
            / clustersize[partcounter];
          traj[0][frame + 1][partcounter] = addvectors(traj[0][frame + 1][partcounter],
            scalevector(factor,drrel));
        }
      }
      for(int partcounter = 0; partcounter < numparticles; partcounter++)
      {
        if(lincombs[vectorcounter][partcounter] != 0)
        {
          overlapping = calcDistances(frame + 1,partcounter);
        }
      }
    }
    overlapping = calcDistances(frame + 1,-1);
    if(updatepot)
    {
      updatePotential(frame);
    }
  }
  traj[0][endframe] = finalframe;
  return overlapping;
}

vector<double> Traj::propagateMirrorDynamics(vector<double> r, vector<double> drrel,
  vector<double> rp, vector<double> rf, double tf, double t, double D, bool absorbing)
{
  vector<double> F = subvectors(rf,r);
  vector<double> deltarf = scalevector(dotvectors(rp,rf) / squarevector(rp)
    - 1.0,rp);
  // determine if the final configuration lies on the allowed side of the mirror plane
  // if it is, propagate dynamics accordingly
  if(dotvectors(rp, deltarf) > 0)
  {
    F = scalevector(1.0 / (tf - t),F);
    double tanch = tanh((dotvectors(deltarf,r) - dotvectors(deltarf,rf) +
      squarevector(deltarf)) / (2 * D * (tf - t)));
    double factor;
    if(absorbing)
    {
      factor = - (1.0 - 1.0 / tanch) / (tf - t);
    }
    else
    {
      factor = - (1.0 - tanch) / (tf - t);
    }
    F = addvectors(F,scalevector(factor,deltarf));
  }
  // if the final configuration is inaccessible, define a fictitous final configuration
  // and final time and then propagate dynamics accordingly
  else
  {
    vector<double> rfprime = subvectors(rf,deltarf);
    F = subvectors(rfprime,r);
    double tfprime = t + (tf - t) * sqrt(squarevector(subvectors(rfprime, r))
      / squarevector(subvectors(rf, r)));
    if(tfprime < t + 5 * dt)
    {
      tfprime = t + 5 * dt;
    }
    F = scalevector(1.0 / (tfprime - t),F);
    if(absorbing)
    {
      F = addvectors(F,scalevector(- 2.0 * D /
        dotvectors(rp, subvectors(rfprime,r)),rp));
    }
  }
  drrel = addvectors(drrel,scalevector(dt,F));
  // enforce reflecting boundary conditions if necessary
  if(!absorbing)
  {
    double rpnorm = sqrt(squarevector(rp));
    double rproj = dotvectors(rp,r) / rpnorm;
    double drproj = dotvectors(rp,drrel) / rpnorm;
    if(rproj + drproj < rpnorm)
    {
      drrel = addvectors(drrel, scalevector(2.0 * (rpnorm - rproj - drproj)
        / rpnorm,rp));
    }
  }
  return drrel;
}

//*CHARACTERISTIC FUNCTIONS*//
// helper function; determines if two particles are within 1.2 rm of each other
bool Traj::near(int frame, int particleA, int particleB)
{
  bool near = (distances[0][frame][particleA][particleB] < 1.2 * pow(2.0,1.0 / 6.0)
    * (sigma[particleA] + sigma[particleB]));
  return near;
}

// helper function; determines if two particles are more than 1.8 rm away from each other
bool Traj::far(int frame, int particleA, int particleB)
{
  double r = sqrt(squarevector(subvectors(traj[0][frame][particleA],
    traj[0][frame][particleB])));
  if(r > 1.8 * pow(2.0, 1.0 / 6.0) * (sigma[particleA] + sigma[particleB]))
  {
    return 1;
  }
  return 0;
}

// helper function; returns the number of nearest neighbors a given particle has
int Traj::numNeighbors(int frame, int particle)
{
  int neighbors = 0;
  for(int partcounter = 0; partcounter < numparticles; partcounter++)
  {
    if(partcounter != particle)
    {
      if(near(frame, particle, partcounter))
      {
        neighbors++;
      }
    }
  }
  return neighbors;
}

// helper function; returns the number of particles that are far from a given particle
int Traj::numFar(int frame, int particle)
{
  int farparts = 0;
  for(int partcounter = 0; partcounter < numparticles; partcounter++)
  {
    if(partcounter != particle)
    {
      if(far(frame, particle, partcounter))
      {
        farparts++;
      }
    }
  }
  return farparts;
}

// determines whether a trajectory has a given particle with 6 neighbors
bool Traj::inCenter(int frame, int particle)
{
  return (numNeighbors(frame, particle) == 6);
}

// determines whether a trajectory has a particle with no nearest neighbors
bool Traj::scattered(int frame)
{
  bool scattered = 0;
  for(int partcounter = 0; partcounter < numparticles; partcounter++)
  {
    scattered = scattered || (numNeighbors(frame, partcounter) == 0);
  }
  return scattered;
}

// determines whether all particles are unphysically stacked together
bool Traj::stacked(int frame)
{
  bool stacked = 1;
  for(int partcounter = 0; partcounter < numparticles; partcounter++)
  {
      stacked = stacked && (numNeighbors(frame, partcounter) == 6);
  }
  return stacked;
}

// helper function; returns a vector containing the number of particles with the index of neighbors
vector<int> Traj::occupationVector(int frame)
{
  vector<int> occupations(numparticles, 0);
  for(int partcounter = 0; partcounter < numparticles; partcounter++)
  {
    int neighbors = numNeighbors(frame, partcounter);
    occupations[neighbors]++;
  }
  return occupations;
}

// gives an integer index corresponding to the state the system is in
int Traj::state(int frame)
{
  // once gave the following indices: -1 for uncategorized, 0 for C_0^0, 10 for C_0^1,
  // 41 for C_1^4, 31 for C_1^3, 11 for C_1^1, 12 for C_2^1

  // now gives Julia's indices: 0 for uncategorized, 1 for C_0^0, 2 for C_0^1,
  // 3 for C_1^4, 4 for C_1^1, 5 for C_1^3, 6 for C_2^1, 7 for C_2^4, 8 for C_2^0,
  // 9 for C_1^0
  vector<int> occupations = occupationVector(frame);
  vector<int> C0 = {0, 0, 0, 6, 0, 0, 1};
  vector<int> C1 = {0, 0, 3, 1, 2, 1, 0};
  vector<int> C2 = {0, 0, 2, 3, 1, 1, 0};
  vector<int> C3 = {0, 0, 2, 2, 3, 0, 0};
  int neighbors = numNeighbors(frame, 0);
  if(occupations == C0)
  {
    if(neighbors == 6)
    {
      return 1; // C_0^0
    }
    else
    {
      return 2; // C_0^1
    }
  }
  else if(occupations == C1)
  {
    if(neighbors == 5)
    {
      return 3; // C_1^4
    }
    else if(neighbors == 3)
    {
      return 9; // C_1^0
    }
    // find the unique particle with three neighbors
    int threefer = -1;
    for(int partcounter = 1; partcounter < numparticles; partcounter++)
    {
      if(numNeighbors(frame, partcounter) == 3)
      {
        threefer = partcounter;
      }
    }
    if(neighbors == 4 && near(frame, 0, threefer))
    {
      return 4; // C_1^1
    }
    else if(neighbors == 2 && near(frame, 0, threefer))
    {
      return 5; // C_1^3
    }
    else
    {
      return 0; // uncategorized
    }
  }
  else if(occupations == C2)
  {
    if(neighbors == 4)
    {
      return 6; // C_2^1
    }
    else if(neighbors == 5)
    {
      return 7; // C_2^4
    }
    else if(neighbors == 3)
    {
      bool neartwofer = 0;
      for(int partcounter = 1; partcounter < numparticles; partcounter++)
      {
        if(numNeighbors(frame, partcounter) == 2 && near(frame, 0, partcounter))
        {
          neartwofer = 1;
        }
      }
      if (!neartwofer)
      {
        return 8; // C_2^0
      }
    }
    // change if other C2 states are to be accounted for
    else
    {
      return 0; // uncategorized
    }
  }
  // change if C3 states are to be accounted for
  else if(occupations == C3)
  {
    return 0; // uncategorized
  }
  return 0; // uncategorized
}

// determines whether a trajectory is reactive from one state to another, taking
// as input the integer state indices as specified by "state"
bool Traj::reactive(int startstate, int endstate)
{
  return (state(0) == startstate && state(trajlength - 1) == endstate);
}

// determines whether a segment of trajectory is reactive from one state to another, taking
// as input the integer state indices as specified by "state"
bool Traj::reactive(int startstate, int endstate, int startframe, int endframe)
{
  return (state(startframe) == startstate && state(endframe) == endstate);
}

//*PATH ACTION CALCULATIONS*//
// calculates the dimensionless path action for a section of trajectory
double Traj::S(int startframe, int endframe, int particle, int dimension)
{
    double S = 0;
    for(int framecounter = startframe; framecounter < endframe; framecounter++)
    {
      double beta = 1.0 / kT;
      double deltax = traj[0][framecounter + 1][particle][dimension]
        - traj[0][framecounter][particle][dimension];
      S += 0.25 * beta * gamma[particle] * (pow(deltax,2) / dt
          - 2 * deltax * force[0][framecounter][particle][dimension] / gamma[particle]
          + dt * pow(force[0][framecounter][particle][dimension],2)
          / pow(gamma[particle],2));
    }
    return S;
}

// calculates the dimensionless path action for a section of trajectory, assuming
// no forces are acting
double Traj::S0(int startframe, int endframe, int particle, int dimension)
{
    double S0 = 0;
    for(int framecounter = startframe; framecounter < endframe; framecounter++)
    {
      double beta = 1.0 / kT;
      double deltax = traj[0][framecounter + 1][particle][dimension]
        - traj[0][framecounter][particle][dimension];
      S0 += 0.25 * beta * gamma[particle] * (pow(deltax,2) / dt);
    }
    return S0;
}

// more efficiently calculates the difference between S and S0
double Traj::deltaS(int startframe, int endframe, int particle, int dimension)
{
    double deltaS = 0;
    for(int framecounter = startframe; framecounter < endframe; framecounter++)
    {
      double beta = 1.0 / kT;
      double deltax = traj[0][framecounter + 1][particle][dimension]
        - traj[0][framecounter][particle][dimension];
      deltaS += 0.25*beta*gamma[particle]*(-2*deltax
          * force[0][framecounter][particle][dimension] / gamma[particle]
          + dt * pow(force[0][framecounter][particle][dimension],2)
          / pow(gamma[particle],2));
    }
    return deltaS;
}

// calculates the effective action associated with including reflecting boundary
// conditions between all close particle pairs and excluding the component of grad
// lnQ perpendicular to the reflecting barrier in relative coordinate space;
// works for a single step of dynamical propagation, and takes most of the same
// inputs as "propagateMirrorDynamics"
double Traj::Qratio(vector<double> r, vector<double> rnew, vector<double> rp,
  vector<double> rf, double tf, double t, double D, bool absorbing)
{
  vector<double> deltarf = scalevector(dotvectors(rp,rf) / squarevector(rp) - 1.0,rp);
  double Qrat = 1.0;
  // determine if the final configuration lies on the allowed side of the mirror plane
  if (dotvectors(rp, deltarf) > 0)
  {
    vector<double> rm = subvectors(rf, scalevector(2.0, deltarf));
    double sign = 1.0;
    if(absorbing)
    {
      sign = - 1.0;
      Qrat *= erf(sqrt(squarevector(deltarf) / (2 * D * (tf - t))));
      if (tf != t + dt)
      {
        Qrat /= erf(sqrt(squarevector(deltarf) / (2 * D * (tf - t - dt))));
      }
    }
    else
    {
      Qrat /= (exp(- squarevector(subvectors(rf, r)) / (4 * D * (tf - t)))
        + sign * exp(- squarevector(subvectors(rm, r)) / (4 * D * (tf - t))))
        / sqrt(tf - t);
      if (tf != t + dt)
      {
        Qrat *= (exp(- squarevector(subvectors(rf, rnew)) / (4 * D * (tf - t
          - dt))) + sign * exp(- squarevector(subvectors(rm, rnew)) / (4 * D
          * (tf - t - dt)))) / sqrt(tf - t - dt);
      }
    }
  }
  // if the final configuration does not lie on the allowed side of the mirror plane,
  // determine a fictitious final configuration and final time
  else
  {
    vector<double> rfprime = subvectors(rf,deltarf);
    double tfprime = t + (tf - t) * sqrt(squarevector(subvectors(rfprime, r))
      / squarevector(subvectors(rf, r)));
    if(tfprime < t + 5 * dt)
    {
      tfprime = t + 5 * dt;
    }
    if(absorbing)
    {
      Qrat /= - dotvectors(rp, subvectors(rfprime,r));
      Qrat *= - dotvectors(rp, subvectors(rfprime,rnew));
    }
    Qrat /= exp(- squarevector(subvectors(rfprime, r)) / (4 * D * (tfprime
      - t))) / sqrt(tfprime - t);
    Qrat *= exp(- squarevector(subvectors(rfprime, rnew)) / (4 * D * (tfprime
      - t - dt))) / sqrt(tfprime - t - dt);
  }
  return Qrat;
}

//*TRANSITION PATH SAMPLING*//

// dynamically shoots forward from a randomly selected frame in the trajectory, overwriting
// the previous states that were farther along
bool Traj::shootForward()
{
    choice = 0;
    double random = gsl_rng_uniform(generator);
    int frame = (trajlength - 1) * random;
    moveframe = frame;
    propagateDynamics(frame, trajlength - 1, calcvelocities, 1);
    bool selected = reactive(startstate,endstate);
    resetTraj(frame + 1, trajlength - 1, selected, 1);
    return selected;
}

// dynamically shoots backward from a randomly selected frame in the trajectory, overwriting
// the previous states farther back in the trajectory
bool Traj::shootReverse()
{
    choice = 1;
    double random = gsl_rng_uniform(generator);
    int frame = ((trajlength - 1) * random) + 1;
    moveframe = frame;
    propagateDynamics(frame, 0, calcvelocities, 1);
    bool selected = reactive(startstate,endstate);
    resetTraj(0, frame - 1, selected, 1);
    return selected;
}

// shoots both ways from a randomly selected frame in the trajectory
bool Traj::shootTwoWays()
{
    double random = gsl_rng_uniform(generator);
    int frame = trajlength * random;
    moveframe = frame;
    propagateDynamics(frame, 0, calcvelocities, 1);
    propagateDynamics(frame, trajlength - 1, calcvelocities, 1);
    bool selected = reactive(startstate,endstate);
    resetTraj(0, frame - 1, selected, 1);
    resetTraj(frame + 1, trajlength - 1, selected, 1);
    return selected;
}

// shifts (reptates) forward, making a randomly selected state the new final state and dynamically
// filling in new states back to the new start of the trajectory
bool Traj::shiftForward()
{
    choice = 2;
    double random = gsl_rng_uniform(generator);
    int frame = (trajlength - 1) - shiftlength * random;
    moveframe = frame;
    // perform the shifting
    for (int counter = 0; counter <= frame; counter++)
    {
      traj[0][trajlength-1-counter] = traj[0][frame-counter];
      updatePotential(trajlength-1-counter);
    }
    // fill in the gap
    propagateDynamics(trajlength - frame - 1, 0, calcvelocities, 1);
    bool selected = reactive(startstate,endstate);
    resetTraj(0, trajlength - 1, selected, 1);
    return selected;
}

// shifts (reptates) backward, making a randomly selected state the new initial state and dynamically
// filling in new states forward to the new end of the trajectory
bool Traj::shiftReverse()
{
    choice = 3;
    double random = gsl_rng_uniform(generator);
    int frame = shiftlength * random;
    moveframe = frame;
    // perform the shifting
    for (int counter = frame; counter < trajlength; counter++)
    {
      traj[0][counter-frame] = traj[0][counter];
      updatePotential(counter - frame);
    }
    // fill in the gap
    propagateDynamics(trajlength - frame - 1, trajlength - 1, calcvelocities, 1);
    bool selected = reactive(startstate,endstate);
    resetTraj(0, trajlength - 1, selected, 1);
    return selected;
}

// connects two specified points in the trajectory with a Brownian bridge
// copy index denotes the trajectory to be used as backup (DO NOT USE 0)
bool Traj::bridge(int startframe, int endframe, int copyindex)
{
  choice = 6;
  assert(copyindex != 0);
  assert(startframe < endframe);
  moveframe = startframe; // for reporting the frames outside the function
  endmoveframe = endframe; // for reporting the frames outside the function
  double oldDeltaS = 0;
  double oldS0 = 0;
  double oldS = 0;
  double newDeltaS = 0;
  double newS0 = 0;
  double newS = 0;
  // determine which particles have more than four nearest neighbors, and
  // calculate deltaS for those that do not and S for those that do
  vector<bool> constrained(numparticles, 0);
  int maxunconstrained = -1;
  for (int partcounter = 0; partcounter < numparticles; partcounter++)
  {
    constrained[partcounter] = (numNeighbors(startframe, partcounter) > 4);
    if (!constrained[partcounter])
    {
      maxunconstrained = partcounter;
    }
    for (int dimcounter = 0; dimcounter < configdim; dimcounter++)
    {
      if (!constrained[partcounter])
      {
        oldDeltaS += deltaS(startframe, endframe, partcounter, dimcounter);
        if (copyindex > 1)
        {
          oldS0 += S0(startframe, endframe, partcounter, dimcounter);
        }
      }
      else
      {
        oldS += S(startframe, endframe, partcounter, dimcounter);
      }
    }
  }

  // propagate bridge dynamics, calculating the potential only once all particles
  // have been propagated
  for (int partcounter = 0; partcounter < numparticles; partcounter++)
  {
      if (!constrained[partcounter] && partcounter != maxunconstrained)
      {
        propagateBridgeDynamics(startframe, endframe, partcounter, 0, 0);
      }
      else if (partcounter == maxunconstrained)
      {
        propagateBridgeDynamics(startframe, endframe, partcounter, 0, 1);
      }
  }

  // determine which particles have more than four nearest neighbors, and
  // calculate deltaS for those that do not and S for those that do
  for (int partcounter = 0; partcounter < numparticles; partcounter++)
  {
    for (int dimcounter = 0; dimcounter < configdim; dimcounter++)
    {
      if (!constrained[partcounter])
      {
        newDeltaS += deltaS(startframe, endframe, partcounter, dimcounter);
        if (copyindex > 1)
        {
          newS0 += S0(startframe, endframe, partcounter, dimcounter);
        }
      }
      else
      {
        newS += S(startframe, endframe, partcounter, dimcounter);
      }
    }
  }
  // accept or reject the bridge move
  bool selected = exp(oldDeltaS + oldS - newDeltaS - newS)
    > gsl_rng_uniform(generator);
  // if copyindex > 1, as is the case in annealing, determine the total action
  // of the new bridge segment
  if (copyindex > 1)
  {
    if (selected)
    {
      #pragma omp atomic
      successes[5]++;
      threadSvector[copyindex - 2] = newS + newDeltaS + newS0;
    }
    else
    {
      threadSvector[copyindex - 2] = oldS + oldDeltaS + oldS0;
    }
    #pragma omp atomic
    attempts[5]++;
  }

  resetTraj(startframe + 1, endframe - 1, selected, copyindex);

  return selected;
}

// connects two randomly selected points with a Brownian bridge
bool Traj::bridge()
{
  choice = 6;
  int length = minlength + gsl_rng_uniform(generator) * (maxlength - minlength);
  int startframe = (trajlength - 1 - length) * gsl_rng_uniform(generator);
  int endframe = startframe + length;
  return bridge(startframe, endframe, 1);
}

// connects two specified points in the trajectory with a hard sphere Brownian bridge
// copy index denotes the trajectory to be used as backup (DO NOT USE 0)
bool Traj::HSbridge(int startframe, int endframe, int copyindex)
{
  choice = 6;
  assert(copyindex != 0);
  assert(startframe < endframe);
  moveframe = startframe; // for reporting the frames outside the function
  endmoveframe = endframe; // for reporting the frames outside the function
  double oldDeltaS = 0;
  double oldS0 = 0;
  double newDeltaS = 0;
  double newS0 = 0;
  // calculate deltaS and S0 (for later calculating S) for the initial trajectory
  for (int partcounter = 0; partcounter < numparticles; partcounter++)
  {
    for (int dimcounter = 0; dimcounter < configdim; dimcounter++)
    {
      oldDeltaS += deltaS(startframe, endframe, partcounter, dimcounter);
      if (copyindex > 1)
      {
        oldS0 += S0(startframe, endframe, partcounter, dimcounter);
      }
    }
  }

  // determine Qrat for the previous trajectory
  prevQrat = 1.0;
  double tf = dt * endframe;
  vector< vector<double> > finalframe = traj[0][endframe];
  for(int frame = startframe; frame < endframe; frame++)
  {
    double t = dt * frame;
    pseudodistances = distances[0][frame];
    for(int partcounter = 0; partcounter < numparticles; partcounter++)
    {
      dr[partcounter] = subvectors(traj[0][frame + 1][partcounter],traj[0]
        [frame][partcounter]);
    }
    vector< vector<double> > nextframe = traj[0][frame + 1];
    traj[0][frame + 1] = traj[0][frame];
    lincombs.resize(0);
    for(int partcounter = 0; partcounter < numparticles; partcounter++)
    {
      cluster[partcounter] = partcounter + 1;
      clustersize[partcounter] = 1;
    }
    for(int vectorcounter = 0; vectorcounter < numparticles - 1; vectorcounter++)
    {
      vector<int> closestpair = findShortestDist(frame + 1);
      int partnerone = closestpair[0];
      int partnertwo = closestpair[1];
      lincombs.push_back(vector<int>(numparticles,0));
      // determine the coefficients of the current Jacobi vector
      for(int partcounter = 0; partcounter < numparticles; partcounter++)
      {
        if(cluster[partcounter] == cluster[partnerone])
        {
          lincombs[vectorcounter][partcounter] = - clustersize[partcounter];
          if(partcounter != partnerone)
          {
            clustersize[partcounter] += clustersize[partnertwo];
          }
        }
        if(cluster[partcounter] == cluster[partnertwo])
        {
          lincombs[vectorcounter][partcounter] = clustersize[partcounter];
          if(partcounter != partnertwo)
          {
            clustersize[partcounter] += clustersize[partnerone];
            cluster[partcounter] = cluster[partnerone];
          }
        }
      }
      int newsize = clustersize[partnerone] + clustersize[partnertwo];
      clustersize[partnerone] = newsize;
      clustersize[partnertwo] = newsize;
      cluster[partnertwo] = cluster[partnerone];
      // determine the inputs to propagateMirrorDynamics
      vector<double> r(configdim,0.0);
      vector<double> drrel(configdim,0.0);
      vector<double> rf(configdim,0.0);
      double D = 0;
      for(int partcounter = 0; partcounter < numparticles; partcounter++)
      {
        if (lincombs[vectorcounter][partcounter] != 0)
        {
          r = addvectors(r,scalevector(1.0 / lincombs[vectorcounter][partcounter],
            traj[0][frame + 1][partcounter]));
          drrel = addvectors(drrel,scalevector(1.0 / lincombs[vectorcounter][partcounter],
            dr[partcounter]));
          rf = addvectors(rf,scalevector(1.0 / lincombs[vectorcounter][partcounter],
            finalframe[partcounter]));
          D += fabs(1.0 / lincombs[vectorcounter][partcounter]) * kT / gamma[partcounter];
        }
      }
      vector<double> deltar = subvectors(traj[0][frame + 1][partnertwo],
        traj[0][frame + 1][partnerone]);
      vector<double> axis = subvectors(r,deltar);
      double deltarnorm = sqrt(squarevector(deltar));
      vector<double> rp = scalevector(dotvectors(axis,deltar) / pow(deltarnorm,2)
        + (sigma[partnerone] + sigma[partnertwo]) / deltarnorm, deltar);
      prevQrat *= Qratio(r, addvectors(r,drrel), rp, rf, tf, t, D, 0);
      for(int partcounter = 0; partcounter < numparticles; partcounter++)
      {
        if(lincombs[vectorcounter][partcounter] < 0 && frame != endframe - 1)
        {
          double factor = -1.0 - 1.0 * lincombs[vectorcounter][partcounter]
            / clustersize[partcounter];
          traj[0][frame + 1][partcounter] = addvectors(traj[0][frame + 1][partcounter],
            scalevector(factor,drrel));
        }
        if(lincombs[vectorcounter][partcounter] > 0 && frame != endframe - 1)
        {
          double factor = 1.0 - 1.0 * lincombs[vectorcounter][partcounter]
            / clustersize[partcounter];
          traj[0][frame + 1][partcounter] = addvectors(traj[0][frame + 1][partcounter],
            scalevector(factor,drrel));
        }
      }
      for(int partcounter = 0; partcounter < numparticles; partcounter++)
      {
        if(lincombs[vectorcounter][partcounter] != 0)
        {
          calcDistances(frame + 1,partcounter);
        }
      }
    }
    traj[0][frame + 1] = nextframe;
  }

  // propagate dynamics
  bool overlapping = propagateHardSphereDynamics(startframe, endframe, 0, 0, 1);
  wasted = overlapping;

  // determine deltaS and S0 (for calculating S) of the new trajectory
  for (int partcounter = 0; partcounter < numparticles; partcounter++)
  {
    for (int dimcounter = 0; dimcounter < configdim; dimcounter++)
    {
      newDeltaS += deltaS(startframe, endframe, partcounter, dimcounter);
      if (copyindex > 1)
      {
        newS0 += S0(startframe, endframe, partcounter, dimcounter);
      }
    }
  }
  // accept or reject the hard sphere bridge move
  bool selected = !overlapping && exp(oldDeltaS - newDeltaS) * prevQrat / Qrat
    > gsl_rng_uniform(generator);
  // if copyindex > 1, as is the case in annealing, determine the total action
  // of the new bridge segment
  if (copyindex > 1)
  {
    if (selected)
    {
      #pragma omp atomic
      successes[5]++;
      threadSvector[copyindex - 2] = newDeltaS + newS0;
    }
    else
    {
      threadSvector[copyindex - 2] = oldDeltaS + oldS0;
    }
    #pragma omp atomic
    attempts[5]++;
  }

  if(selected)
  {
    prevQrat = Qrat;
  }

  resetTraj(startframe + 1, endframe - 1, selected, 1);

  return selected;
}

// connects two randomly selected points with a hard sphere Brownian bridge
bool Traj::HSbridge()
{
  choice = 6;
  /*
  int forwardstate;
  int backwardstate;
  int minstartframe = 0;
  int maxendframe = trajlength - 1;
  for (int framecounter = 0; framecounter < trajlength; framecounter++)
  {
    forwardstate = state(framecounter);
    backwardstate = state(trajlength - 1 - framecounter);
    if (forwardstate == startstate)
    {
      minstartframe = framecounter;
    }
    if (backwardstate == endstate)
    {
      maxendframe = trajlength - 1 - framecounter;
    }
  }
  int startframe = minstartframe;
  int endframe = maxendframe;
  if (maxendframe - minstartframe > maxlength)
  {
    startframe = minstartframe + (maxendframe - minstartframe - maxlength)
      * gsl_rng_uniform(generator);
    endframe = startframe + maxlength;
  }
  */
  int length = minlength + gsl_rng_uniform(generator) * (maxlength - minlength);
  int startframe = (trajlength - 1 - length) * gsl_rng_uniform(generator);
  int endframe = startframe + length;
  return HSbridge(startframe, endframe, 1);
}

// employs an annealing strategy to generate a trajectory by raising and lowering
// the temperature through a series of intermediate trajectories
bool Traj::anneal()
{
  choice = 4;
  // select starting and ending frames
  int forwardstate;
  int backwardstate;
  int minstartframe = 0;
  int maxendframe = trajlength - 1;
  for(int framecounter = 0; framecounter < trajlength; framecounter++)
  {
    forwardstate = state(framecounter);
    backwardstate = state(trajlength - 1 - framecounter);
    if (forwardstate == startstate)
    {
      minstartframe = framecounter;
    }
    if (backwardstate == endstate)
    {
      maxendframe = trajlength - 1 - framecounter;
    }
  }
  int startframe = minstartframe;
  int endframe = maxendframe;
  if (maxendframe - minstartframe > anneallength)
  {
    startframe = (maxendframe - anneallength - minstartframe) * gsl_rng_uniform(generator);
    endframe = startframe + anneallength;
  }
  int length = endframe - startframe;

  double sbaths = 0.0;
  double kTorig = kT;
  double r = pow(Tmax / kTorig, 2.0 / (numbridges - 2 + numbridges % 2));
  double prevtemp = kT; // temperature before last temperature increase

  double currentS = 0;

  for(int tempcounter = 1; tempcounter <= numbridges; tempcounter++)
  {
    // determine the temperature of the step, employing a geometric progression
    if(tempcounter == 1){}
    else if(tempcounter < 0.5 * numbridges + 1)
    {
      kT *= r;
    }
    else if(tempcounter == 0.5 * numbridges + 1){}
    else
    {
      kT /= r;
    }
    sbaths += (1 - prevtemp / kT) * currentS;
    prevtemp = kT;
    int bridgelength = minlength + (maxlength - minlength) * gsl_rng_uniform(generator);
    int totalbridges = length / bridgelength + tempcounter % 2;
    bridgeSvector.resize(totalbridges + 1, 0);

    currentS = 0;

    //BEGIN PARALLELIZATION HERE//
    #pragma omp parallel for
    // propagate bridge dynamics for the full chain
    for(int bridgecounter = 0; bridgecounter <= totalbridges; bridgecounter++)
    {
      int threadnum = omp_get_thread_num();
      int bridgestartframe;
      int bridgeendframe;
      if(tempcounter % 2 == 1 && bridgecounter == 0)
      {
        bridgestartframe = startframe;
        bridgeendframe = startframe + bridgelength / 2;
      }
      else if(tempcounter % 2 == 1)
      {
        bridgestartframe = startframe + bridgelength / 2
          + bridgecounter * bridgelength;
        bridgeendframe = startframe + bridgelength / 2
          + (bridgecounter + 1) * bridgelength;;
      }
      else
      {
        bridgestartframe = startframe + bridgecounter * bridgelength;
        bridgeendframe = startframe + (bridgecounter + 1) * bridgelength;
      }
      if(bridgeendframe > endframe)
      {
        bridgeendframe = endframe;
      }
      // prepare storage vector
      resetTraj(bridgestartframe, bridgeendframe, 1, threadnum + 2);
      if(bridgestartframe < bridgeendframe)
      {
        bridge(bridgestartframe, bridgeendframe, threadnum + 2);
        bridgeSvector[bridgecounter] = threadSvector[threadnum];
      }
    }
    //END PARALLELIZATION HERE
    // reconstruct S from the actions of each individual segment
    for(int bridgecounter = 0; bridgecounter <= totalbridges; bridgecounter++)
    {
      currentS += bridgeSvector[bridgecounter];
    }
  }

  // reset the temperature and accept or reject the move
  kT = kTorig;
  bool selected = (exp(sbaths) > gsl_rng_uniform(generator));
  resetTraj(startframe + 1, endframe - 1, selected, 0);
  return selected;
}

// randomly shoot, shift, anneal, or bridge with a given set of probabilities
// note: probabilities that do not add to 1 will be complemented by remainder bridge moves
bool Traj::moveRandom(double shootprob, double shiftprob, double annealprob)
{
    wasted = 0;
    double randy = gsl_rng_uniform(generator);
    bool randbit = static_cast<bool>(rand()%2);
    bool selected;
    if (randy < shootprob)
    {
      if (randbit == 0)
      {
        selected = shootForward();
      }
      else
      {
        selected = shootReverse();
      }
    }
    else if (randy >= shootprob && randy < shootprob + shiftprob)
    {
      if (randbit == 0)
      {
        selected = shiftForward();
      }
      else
      {
        selected = shiftReverse();
      }
    }
    else if (randy >= shootprob + shiftprob && randy < shootprob + shiftprob + annealprob)
    {
      selected = anneal();
    }
    else
    {
      // selected = bridge();
      // selected = HSbridge(0, trajlength - 1, 1);
      selected = HSbridge();
    }
    if (selected)
    {
      successes[choice]++;
    }
    attempts[choice]++;
    return selected;
}

  //*TPS STATISTICS*//

  // output to screen success statistics for move types
  void Traj::TPSstats()
  {
    vector<double> rates(7, 0);
    for (int choicecounter = 0; choicecounter < 7; choicecounter++)
    {
      if (attempts[choicecounter] != 0)
      {
        rates[choicecounter] = 1.0 * successes[choicecounter] / attempts[choicecounter];
      }
      switch (choicecounter)
      {
        case 0: cout << "Forward shooting ";
                break;
        case 1: cout << "Reverse shooting ";
                break;
        case 2: cout << "Forward shifting ";
                break;
        case 3: cout << "Reverse shifting ";
                break;
        case 4: cout << "Annealing moves ";
                break;
        case 5: cout << "Substituent brownian bridges ";
                break;
        case 6: cout << "Independent brownian bridges ";
                break;
      }
      cout << "succeeded " << rates[choicecounter] * 100.0 << "% of the time.  ("
      << successes[choicecounter] << "/" << attempts[choicecounter] << ")" << '\n';
    }
  }

  // print success statistics for move types to a given file name
  void Traj::TPSstats(string filename)
  {
    ofstream ratefile;
    ratefile.open(filename);
    vector<double> rates(7, 0);
    for (int choicecounter = 0; choicecounter < 7; choicecounter++)
    {
      if (attempts[choicecounter] != 0)
      {
        rates[choicecounter] = 1.0 * successes[choicecounter] / attempts[choicecounter];
      }
      switch (choicecounter)
      {
        case 0: ratefile << "Forward shooting ";
                break;
        case 1: ratefile << "Reverse shooting ";
                break;
        case 2: ratefile << "Forward shifting ";
                break;
        case 3: ratefile << "Reverse shifting ";
                break;
        case 4: ratefile << "Annealing moves ";
                break;
        case 5: ratefile << "Substituent brownian bridges ";
                break;
        case 6: ratefile << "Independent brownian bridges ";
                break;
      }
      ratefile << "succeeded " << rates[choicecounter] * 100.0 << "% of the time.  ("
      << successes[choicecounter] << "/" << attempts[choicecounter] << ")" << '\n';
    }
    ratefile.close();
  }
