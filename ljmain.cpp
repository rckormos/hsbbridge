#include "ljtps.h"
#include <sys/types.h>
#include <sys/stat.h>

int main()
{
  /***PARAMETERS***/
  // physical
  double kT = 0.10;
  double dt = 0.001;
  double gamma = 1.0;
  double epsilon = 1.0;
  double sigma = 0.5;
  int trajlength = 5001;
  int configdim = 2;
  int numparticles = 7;
  bool calcvelocities = 0;

  // shifting-specific
  int shiftlength = trajlength * 0.1;

  // annealing-specific
  int anneallength = 10;
  int minlength = 50;
  int maxlength = 100;
  int numbridges = 100;
  double Tmax = 0.10;

  // sampling
  int nmoves = 10000;
  double probshoot = 0.45;
  double probshift = 0.45;
  double probanneal = 0.0;

  int startstate = 1;
  int endstate = 3;

  int trajseed = 1;
  /****************************************************************************/
  clock_t thistime;
  thistime = clock();

  ofstream accofminlength;
  accofminlength.open("accofminlength.txt");

  // iterate through hard sphere bridge lengths by multiples of ten
  for(int lengthcounter = 10; lengthcounter <= 150; lengthcounter += 10)
  {
    cout << "lengthcounter = " << lengthcounter << '\n';
    minlength = lengthcounter;
    maxlength = lengthcounter;
    Traj tpstrj(kT, dt, gamma, epsilon, sigma, trajlength, configdim,
      numparticles, calcvelocities, shiftlength, anneallength, minlength,
      maxlength, numbridges, Tmax, startstate, endstate, trajseed);
    cout << "Beginning TPS." << '\n';
    bool percentflag = 0;
    int numwasted = 0;
    for(int movecounter = 0; movecounter < nmoves; movecounter++)
    {
      bool success = tpstrj.moveRandom(probshoot, probshift, probanneal);
      if(success)
      {
        // tpstrj.printTrajPos(movecounter + 1, 0, 0, trajlength - 1);
      }
      if(tpstrj.wasted)
      {
        numwasted++;
      }

      // output percent completion
      if(movecounter % (nmoves / 100) == 0 && percentflag == 0)
      {
        cout << (movecounter / (nmoves / 100)) << "%" << '\n';
        /*
          tpstrj.printTrajPos(movecounter / (nmoves / 100) + 1, 0, 0,
            trajlength - 1);
        */
        percentflag = 1;
      }
      else if(movecounter % (nmoves / 100) != 0)
      {
        percentflag = 0;
      }
    }
    cout << "100%" << '\n';
    thistime = clock() - thistime;
    tpstrj.TPSstats();

    accofminlength << lengthcounter << '\t' << 1.0 * tpstrj.successes[6]
      / tpstrj.attempts[6] << '\n';

    cout << "Trajectories were wasted " << 100.0 * numwasted
      / tpstrj.attempts[6] << "% of the time." << '\n';
  }

  accofminlength.close();

  // output process time
  cout << "Process took " << static_cast<double>(thistime) /
    static_cast<double>(CLOCKS_PER_SEC) << " seconds." << '\n' << '\n';

  return 0;
}
