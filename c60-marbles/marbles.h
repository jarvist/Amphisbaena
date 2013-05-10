/* Amphisbaena:- (pronounced am.fis.BEEN.uh)
 * A two-headed Slithering Snake implementation designed to simulate morphologies of Organic Solar Cells.
 * Code by Jarvist Frost, based on algorithm and physical setup by Dr Jenny Nelson & Felix Rickerman.
 * File begun June 16th 2005
 * 
 * Molecular Electronic Materials and Devices
 * Experimental Solid State Physics
 * Blackett Laboratory
 * Imperial College, London
 */

//Start of Simulation Configuration
#define X 1000  //lattice sizes
#define Y 1000 //40
#define Z 41
//In 2D console display, Z is the number of lines, X is the number of Characters. Y goes 'into' the screen & should be set to 1

#define MAX_SNAKES (7*X*Y) //total number of C60 bucky balls 
#define MAX_SEGMENTS 28 //cutoff of maximum number of segments in a snake 
                        //  - will cut off gaussian distribution at this value

#define l_0 10 //average length of snakes in gaussian distribution
#define sigma 0.005 //500 //snake lengths are gaussian distribution of this standard deviation
//set to a very small value to have snakes of exactly the 'average' length 

#define beta 0.02 //BETA //5 //0.01 //25.0		// B=1/T  T=temperature of the lattice, in units of k_B
/* From Felix's code; T is
 * very cold= 2.0
 * cold = 20.0
 * intermediate = 200.0 
 * hot= 2000.0
*/
#define DENSITY 0.20 //density of snake material as a fraction of the whole

#define SLITHER_PRINT 100
#define MC_BETWEEN_DROPS 5000
#define TOTAL_SLITHERS (((long long)MAX_SNAKES*(long long)MC_BETWEEN_DROPS) + (long long)(200 * 1000000))

float iE[2][2] =		//interaction energies between the different plastics
{
  {0.0, 0.0},
  {0.0, 500.0} //+ve energies are attractive, -ve energies are repulsive
};

#define SUBSTRATE_E 100.0

//End of Simulation Configuration Parameters

//Start of Global Variables + Types
int lattice[X][Y][Z];
int perc[X][Y][Z];
int num_snakes = 0; //start with zero initial snakes - global to allow selection of snakes using this
int num_segments=0; //total number of snake segments. Used for percolation stats.

struct coord{ int x;  int y;  int z; };

struct coord marbles[MAX_SNAKES];

//End of Global Variables + Types
//Start of Prototypes
float lattice_energy();
void empty_lattice ();
void print_lattice ();
void print_xmakemol ();
void print_povray ();
void print_lattice_pnm ();
void print_lattice_pnm_file(int i);
void print_snakes ();
void print_a_snake (int id);
void placemarble();

void wriggle (int wriggles);
void crawl(int x, int y, int z);
void percolate();  

//End of Prototypes

//long long rand_count=0;
long float_seed = -1;
double
rand_float ()
{
//   rand_count++; //count random numbers used
//  return (ran2 (&float_seed));
return(genrand_real2());
}

int
rand_int (int max)
{
// rand_count++; //count random numbers used
  return (int) (gen_rand32()%max);
  //return (int) (genrand_real2() * (double) max);
}
  
