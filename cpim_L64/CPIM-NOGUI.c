// C code to expand the Contact Process code to incorporate 
// a basic Ising Model in a Fair fashion

#include <stdlib.h>
#include "mt64.h"    /* Pseudo-random number generation MT library (64 bit) */
#include <math.h>    /* Math to transform random n from continuous to discrete */
#include <time.h>    /* Used to seed pseudo-random number generator */
#include <stdio.h>

/* Lattice Size */
#define X_SIZE 64
#define Y_SIZE 64

/* Defaulfs */
// number of replicas for each parameter set 
#define REPLICAS 10
// default birth/colonization rate/probability
#define BETA   0.01
// default mortality/extinction rate/probability
#define DELTA  0.001
// default differentiation rate/probability
#define ALPHA  0.1
// strength of the coupling in positive terms (J =  -1*COUPLING kBT units)
// we should have 1/2 if we do not want to douple count pairs
// therefore
// use a positive number!
#define COUPLING (1)
// default Temperature
#define TEMPERATURE 2.269
// default interaction radius
#define RADIUS 1
// default initial condition chosen
#define INIT 1

/* Structure with the simulation data */
struct simulation
  {
  int lattice_configuration[X_SIZE][Y_SIZE]; /* Store latice configuration */
  int init_option;            /* Choice of initial condition*/
  int generation_time;        /* Generations simulated */
  int Ising_neighboorhood;    /* Ising Neighboorhood: r=1 (NN) vs r=2 (NNN)*/
  int occupancy;              /* Lattice occupancy */
  int vacancy;                /* Lattice vacancy*/
  int up;                     /* Number of spins in the up   (+1) state */
  int down;                   /* Number of spins in the down (-1) state */
  double birth_rate;          /* Contact Process' birth */
  double death_rate;          /* Contact Process' death */
  double differentiation_rate;/* Differentiation into spin state */
  double T;                   /* Ising's temperature */
  double J;                   /* Ising's coupling: ferro (-kB) or anti-ferro (+kB) */
  double lamda_rate;          /* Contact-Ising Monte Carlo biass*/
} s ;        // instance s of the structure to hold the simulation


double local_energy (int x, int y)
	{
  // Energy of site at coordinate (x,y)
  double energy; //in kB*T units
	int up = 0;   
	int down = 0; 
	if (s.Ising_neighboorhood == 1)  // Nearest Neighboorhood (NN) has 4 sites
    	{ 
      // we check the South (S) neighboor (#1)
      if (s.lattice_configuration[x][(int)((Y_SIZE + y+1)%Y_SIZE)] == 1){up++;}
    	 else if (s.lattice_configuration[x][(int)((Y_SIZE + y+1)%Y_SIZE)] == -1){down++;}
    	// we check the North (N) neighboor (#2)
      if (s.lattice_configuration[x][(int)((Y_SIZE + y-1)%Y_SIZE)] == 1){up++;}
    	 else if (s.lattice_configuration[x][(int)((Y_SIZE + y-1)%Y_SIZE)] == -1){down++;}
    	// we check the West (W) neighboor  (#3)
      if (s.lattice_configuration[(int)((X_SIZE + x-1)%X_SIZE)][y] == 1){up++;}
    	 else if (s.lattice_configuration[(int)((X_SIZE + x-1)%X_SIZE)][y] == -1){down++;}
    	// we check the East (E) neighboor  (#4)
      if (s.lattice_configuration[(int)((X_SIZE + x+1)%X_SIZE)][y] == 1){up++;}
    	 else if (s.lattice_configuration[(int)((X_SIZE + x+1)%X_SIZE)][y] == -1){down++;}
    	}
    else if (s.Ising_neighboorhood == 2) //Next Nearest Neighboorhood (NNN) has 12 sites
        { 
        // we check the South (S) neighboor (#1)
        if (s.lattice_configuration[x][(int)((Y_SIZE + y+1)%Y_SIZE)] == 1){up++;}
    	   else if (s.lattice_configuration[x][(int)((Y_SIZE + y+1)%Y_SIZE)] == -1){down++;}
    	  // we check the North (N) neighboor (#2)
        if (s.lattice_configuration[x][(int)((Y_SIZE + y-1)%Y_SIZE)] == 1){up++;}
    	   else if (s.lattice_configuration[x][(int)((Y_SIZE + y-1)%Y_SIZE)] == -1){down++;}
    	  // we check the West (W) neighboor (#3)
        if (s.lattice_configuration[(int)((X_SIZE + x-1)%X_SIZE)][y] == 1){up++;}
    	   else if (s.lattice_configuration[(int)((X_SIZE + x-1)%X_SIZE)][y] == -1){down++;}
    	  // we check the East (E) neighboor (#4)
        if (s.lattice_configuration[(int)((X_SIZE + x+1)%X_SIZE)][y] == 1){up++;}
    	   else if (s.lattice_configuration[(int)((X_SIZE + x+1)%X_SIZE)][y] == -1){down++;}
        // we check the South-South (SS) neighboor (#5)
        if (s.lattice_configuration[x][(int)((Y_SIZE + y+2)%Y_SIZE)] == 1){up++;}
    	   else if (s.lattice_configuration[x][(int)((Y_SIZE + y+2)%Y_SIZE)] == -1){down++;}
    	  // we check the North-North (NN) neighboor (#6)
        if (s.lattice_configuration[x][(int)((Y_SIZE + y-2)%Y_SIZE)] == 1){up++;}
    	   else if (s.lattice_configuration[x][(int)((Y_SIZE + y-2)%Y_SIZE)] == -1){down++;}
    	  // we check the West-West (WW) neighboor   (#7)
        if (s.lattice_configuration[(int)((X_SIZE + x-2)%X_SIZE)][y] == 1){up++;}
    	   else if (s.lattice_configuration[(int)((X_SIZE + x-2)%X_SIZE)][y] == -1){down++;}
    	  // we chack the East-East (EE) neighboor   (#8)
        if (s.lattice_configuration[(int)((X_SIZE + x+2)%X_SIZE)][y] == 1){up++;}
    	   else if (s.lattice_configuration[(int)((X_SIZE + x+2)%X_SIZE)][y] == -1){down++;}
    	  // we check the South-West (SW) neighboor  (#9)
        if (s.lattice_configuration[(int)((X_SIZE + x-1)%X_SIZE)][(int)((Y_SIZE + y+1)%Y_SIZE)] == 1){up++;}
    	   else if (s.lattice_configuration[(int)((X_SIZE + x-1)%X_SIZE)][(int)((Y_SIZE + y+1)%Y_SIZE)] == -1){down++;}
    	  // we check the North-East (NE) neighboor  (#10)
        if (s.lattice_configuration[(int)((X_SIZE + x+1)%X_SIZE)][(int)((Y_SIZE + y-1)%Y_SIZE)] == 1){up++;}
    	   else if (s.lattice_configuration[(int)((X_SIZE + x+1)%X_SIZE)][(int)((Y_SIZE + y-1)%Y_SIZE)] == -1){down++;}
    	  // we check the North-West (NW) neighboor  (#11)
        if (s.lattice_configuration[(int)((X_SIZE + x-1)%X_SIZE)][(int)((Y_SIZE + y-1)%Y_SIZE)] == 1){up++;}
    	   else if (s.lattice_configuration[(int)((X_SIZE + x-1)%X_SIZE)][(int)((Y_SIZE + y-1)%Y_SIZE)] == -1){down++;}
       	// we check the South-East (SE) neighboor  (#12)
        if (s.lattice_configuration[(int)((X_SIZE + x+1)%X_SIZE)][(int)((Y_SIZE + y+1)%Y_SIZE)] == 1){up++;}
    	   else if (s.lattice_configuration[(int)((X_SIZE + x+1)%X_SIZE)][(int)((Y_SIZE + y+1)%Y_SIZE)] == -1){down++;}
    	  }
	energy =  s.J * (double) (s.lattice_configuration[x][y] * (up-down));
  return energy;
	}


/* Update function */
int update_lattice (void)
  {
  // int random_neighbor;
  int random_neighbor_state, random_neighboor;
  double random_spin;
  // Energies
  double spin_energy, spin_energy_diff;
  // Probability of reactions
  double transition_probability;
  int random_x_coor, random_y_coor;
  // For the Contact Process we always consider NN interactions
  for (int site = 0; site < (int) (Y_SIZE*X_SIZE); site++)
    {
    /* Pick a random focal site */
    random_x_coor = (int) floor (genrand64_real1 ()* X_SIZE);
    random_y_coor = (int) floor (genrand64_real1 ()* Y_SIZE);
    switch (s.lattice_configuration[random_x_coor][random_y_coor])
      {
      case 0: /* Site is empty */
      /* Chose a random neighbor from the num_neighbors posible ones */
      random_neighboor = (int) floor (genrand64_real3()* 4);
			switch(random_neighboor)
					{
					case 0: // South
							random_neighbor_state = s.lattice_configuration[random_x_coor][(int) ((Y_SIZE + random_y_coor-1)%Y_SIZE)]; 
							break;
					case 1: // North
							random_neighbor_state =	s.lattice_configuration[random_x_coor][(int) ((Y_SIZE + random_y_coor+1)%Y_SIZE)]; 
							break;
					case 2: // East
							random_neighbor_state =	s.lattice_configuration[(int) ((X_SIZE + random_x_coor-1)%X_SIZE)][random_y_coor];
							break;
					case 3: // West
							random_neighbor_state =	s.lattice_configuration[(int) ((X_SIZE + random_x_coor+1)%X_SIZE)][random_y_coor];
							break;
					}
        /* If its random neighbor is occupied: put a copy at the focal site
           with probability brith_rate * dt */
        if (genrand64_real2 () < s.birth_rate)
           {
           switch(random_neighbor_state)
             {
              case 2: 
                s.lattice_configuration[random_x_coor][random_y_coor] = 2;
                s.occupancy ++; s.vacancy --;
               break;
              case 1: 
                s.lattice_configuration[random_x_coor][random_y_coor] = 1;
                s.occupancy ++; s.vacancy --;
                s.up ++;
               break;
              case -1:
                s.lattice_configuration[random_x_coor][random_y_coor] = -1;
                s.occupancy ++;s.vacancy --;
                s.down ++;
               break; 
              case 0:
                s.lattice_configuration[random_x_coor][random_y_coor] = 0;
                break;
             }
          }
        break; /* break case 0 */
      case 2: /* Focal point is in the occupied, undifferentiated state */
        // First we check if the site survives
        // No need for Gillespie as cells are macroscopic compare to its 
        // inner components which can undertake reactions only if the cell
        // indeed exists
        if (genrand64_real2 () < s.death_rate)
                       {
                        s.lattice_configuration[random_x_coor][random_y_coor] = 0;
                        s.occupancy --; s.vacancy ++;
                       }
             else if (genrand64_real2 () < s.differentiation_rate)
                      {
                       /* Set an occupied site in the middle of the lattice */
                       random_spin = (int) ((genrand64_int64 () % 2) * 2) - 1;
                      if (random_spin == 1)
                          {
                           s.lattice_configuration[random_x_coor][random_y_coor] = random_spin;
                           s.up ++;
                           }
                       else if (random_spin == -1)
                           {
                           s.lattice_configuration[random_x_coor][random_y_coor] = random_spin;
                           s.down ++;
                           }
                      }
        break;
      case 1: /* Focal point is in the up (+1) state */
        // We skip Gillespie because of separation of scales
        spin_energy = local_energy (random_x_coor, random_y_coor);
        spin_energy_diff = -(2) * spin_energy;
        transition_probability = exp (-spin_energy_diff/s.T);
        if (genrand64_real2 () < s.death_rate)
                        {
                        s.lattice_configuration[random_x_coor][random_y_coor] = 0;
                        s.occupancy --; s.vacancy ++;
                        s.up --;
                        }
                else if (spin_energy_diff < 0 ||
                                              genrand64_real2 () < transition_probability)
                        {
                        s.lattice_configuration[random_x_coor][random_y_coor] = -1;
                        s.up --;
                        s.down ++;
                        }
        break;
      case -1: /* Focal point is in the down (-1) state */
        // We skip Gillespie because of separation of scales
        spin_energy = local_energy (random_x_coor, random_y_coor);
        spin_energy_diff = -(2) * spin_energy;
        transition_probability = exp (-spin_energy_diff/s.T);
        if (genrand64_real2 () < s.death_rate)
                        {
                        s.lattice_configuration[random_x_coor][random_y_coor] = 0;
                        s.occupancy --; s.vacancy ++;
                        s.down --;
                        }
                else if (spin_energy_diff < 0 ||
                                              genrand64_real2 () < transition_probability)
                        {
                        s.lattice_configuration[random_x_coor][random_y_coor] = 1;
                        s.up ++;
                        s.down --;
                        }
        break;
      }
    }
  s.generation_time ++;
  return 0;
}


/* Callback to initialize lattice*/
static void init_lattice ()
  {
  int random_spin;
  int x,y;
  /* Fill the lattice with 0s (unoccupied state) */
  for (x = 0; x < X_SIZE; x++)
    {
    for (y = 0; y < Y_SIZE; y++)
      {
      s.lattice_configuration[x][y]= 0;
      }
    }
  s.occupancy = 0;
  s.up =  0;
  s.down = 0;
  s.vacancy = (int) X_SIZE*Y_SIZE;
  switch(s.init_option)
    {
      case 1:
            /* Set an occupied site in the middle of the lattice */
            random_spin = (int) ((genrand64_int64 () % 2) * 2) - 1;
            if (random_spin == 1)
              {
              s.lattice_configuration[(int) X_SIZE/2][(int) Y_SIZE/2] = random_spin;
              s.up ++; s.vacancy--; s.occupancy++;
              }
            else if (random_spin == -1)
              {
              s.lattice_configuration[(int) X_SIZE/2][(int) Y_SIZE/2] = random_spin;
              s.down ++; s.vacancy --; s.occupancy++;
              }
            break;
      case 2:
            /* Set an undifferentiated site in the middle of the lattice*/
              s.lattice_configuration[(int) X_SIZE/2][(int) Y_SIZE/2] = 2;
              s.vacancy--; s.occupancy++;

            break;
      case 3:
            // Set a small (r=2) cluster with undifferentiated sites in the middle of the lattice
           for (x = (int) X_SIZE/2 - 2 ; x < (int) X_SIZE/2 + 2; x++)
                                for (y = (int) X_SIZE/2 - 2; y < (int) X_SIZE/2 + 2; y++)
                                        {
                                        s.lattice_configuration[x][y]=2;
                                        s.occupancy ++; s.vacancy --;
                                        }
            break;
      case 4:
            for (x = (int) X_SIZE/2 - 2 ; x < (int) X_SIZE/2 + 2; x++)
                                for (y = (int) X_SIZE/2 - 2; y < (int) X_SIZE/2 + 2; y++)
                                    {
                                      random_spin = (int) ((genrand64_int64 () % 2) * 2) - 1;
                                      if (random_spin == 1)
                                        {
                                        s.lattice_configuration[x][y] = random_spin;
                                        s.up ++; s.vacancy--; s.occupancy++;
                                        }
                                        else if (random_spin == -1)
                                           {
                                           s.lattice_configuration[x][y] = random_spin;
                                           s.down ++; s.vacancy --; s.occupancy++;
                                           }
                                    }
         
            break;
      case 5:
            // Set a lattice fully occupied with undufferenciated particles
            for (x = 0; x < (int) X_SIZE; x++)
               for (y = 0; y < (int) Y_SIZE; y++)
                    {
                    s.lattice_configuration[x][y]=2;
                    s.occupancy ++; s.vacancy --;
                    }
            break;
    }
   s.generation_time = 0;
  }


static void initialize_simulation(void)
  {
  /* Initialize Mersenne Twister algorithm for random number genration */
  unsigned int seed = (unsigned int) time (NULL);
  init_genrand64 (seed);

  /* Set default parameters of the simulation */
  //initial condition option
  s.init_option = (int) INIT;

  // Contact Process
  s.birth_rate = (double) BETA;
  s.death_rate = (double) DELTA;

  // Cell differenciation
  s.differentiation_rate = (double) ALPHA;

  // Ising Model
  // interaction radius
  s.Ising_neighboorhood = (int) RADIUS;
  // Temperature
  // Line below is commented because temperature is being swapped in main ()
  // s.T = (double) TEMPERATURE;
  // Spin coupling
  s.J = -1 * (double) COUPLING;
  }


/* Main function spanning a Gtk Application object */
int main (int argc, char **argv)
  {
  /* Create a new text file and populate with the data headers */
  FILE *datafile;
	datafile = fopen(argv[1], "w");
	fprintf(datafile, "replica,lattice_size,temperature,generation_time,occupancy,up,down,magnetization\n");
	fclose(datafile);

  /* The following variables are used for calculating the standard deviation of 
    num_gens_averaged (e.g. 10) consecutive occupancy values */
  int num_gens_averaged = 100;
  double occupancy_arr[num_gens_averaged];
  double average_occupancy;
  double sum;
  double occupancy_variance;

  /* The general strategy is to perform REPLICAS times the whole swap of the 
    the parameters*/
	for (int N = 0; N < REPLICAS; N++)
	  {
		for (double temp = 1.00; temp < 3.6; temp += 0.01)
		  {
			initialize_simulation();
			init_lattice();
      s.T = temp;

      /* We run the simulation until close to the contact process' lattice
        occupancy convergence. */
      while ((double) s.occupancy / (double) (X_SIZE * Y_SIZE) < (double) (0.95 * (1 - (s.death_rate/s.birth_rate))))
        {
        update_lattice();
        }

      /* Criteria for stabilization: last num_gens_averaged (e.g. 100) lattice
        occupancy values' variance smaller or equal than 0.0001 */
      occupancy_variance = 1.0; /* arbitrary initial value to get inside the loop */
      while (occupancy_variance > 0.0001)
        {
        /* We run num_gens_averaged number of generations in order to compute
          the lattice occupancy values' variance over this epoch. */
        average_occupancy = 0.0;
        for (int i = 0; i < num_gens_averaged-1; i++)
          {
          occupancy_arr[i] = (double) s.occupancy / (double) (X_SIZE * Y_SIZE);
          average_occupancy += occupancy_arr[i];
          for (int j = 0; j < 10; j++)
            {update_lattice();}
          }
        sum = 0.0;
        average_occupancy = (double) average_occupancy / (double) num_gens_averaged;
        /* Now we compute the occupancy variance */
        for (int j = 0; j < num_gens_averaged-1; j++)
          {
          sum = sum + pow((occupancy_arr[j] - average_occupancy), 2);
          }
        occupancy_variance = (double) sum / (double) num_gens_averaged;
        }

      /* Lattice's occupacy stability criterium has been met, so we save the 
        next num_gens_averaged generation times' data. */
			for (int t = 0; t < num_gens_averaged; t++)
			  {
				datafile = fopen(argv[1], "a+");
				fprintf(datafile, "%d,%d,%f,%d,%f,%f,%f,%f\n", 
                N, 
                X_SIZE, 
                s.T,
								s.generation_time, 
                (double) s.occupancy / (double) (X_SIZE * Y_SIZE),
								(double) s.up / (double) s.occupancy,
								(double) s.down / (double) s.occupancy,
								(double) (s.up - s.down) / (double) (s.up + s.down));
				fclose(datafile);
        update_lattice();
			  }
		  }
	  }
  return 0;
  }
