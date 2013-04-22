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
#define X 100  //lattice sizes
#define Y 100 //40
#define Z 100 
//In 2D console display, Z is the number of lines, X is the number of Characters. Y goes 'into' the screen & should be set to 1

#define MAX_SNAKES 200000 //cutoff for maximum number of snakes 
                         // - will segfault if attempts to generate more than this
#define MAX_SEGMENTS 90 //cutoff of maximum number of segments in a snake 
                        //  - will cut off gaussian distribution at this value

#define l_0 5 //average length of snakes in gaussian distribution
#define sigma 0.00001 //500 //snake lengths are gaussian distribution of this standard deviation
//set to a very small value to have snakes of exactly the 'average' length 

#define beta 0.005 //BETA //5 //0.01 //25.0		// B=1/T  T=temperature of the lattice, in units of k_B
/* From Felix's code; T is
 * very cold= 2.0
 * cold = 20.0
 * intermediate = 200.0 
 * hot= 2000.0
*/
#define DENSITY 0.40 //density of snake material as a fraction of the whole

#define TOTAL_SLITHERS (400000) //units are now MegaSlithers Jarv 25-6-07
//given 1M slither per sec
//~50'000 M slithers in 12 hrs

#define SLITHER_PRINT 20000

#define iE_SUBSTRATE 50.0 //200.0  //substrate interaction energy - postive is sticky
#define iE_ALIGNED 50.0  //energy benefit for striping in Z axis
#define iE_PI 450.0 //was 50.0 //pi-pi inter-plane stacking energy, added onto usual iE_Interaction below
#define iE_INTERACTION 50.0 //bare cube-cube occupation interaction energy

float iE[2][2] =		//interaction energies between the different plastics
{
  {0.0, 0.0},
  {0.0, 50.0} //+ve energies are attractive, -ve energies are repulsive
};

//End of Simulation Configuration Parameters

//Start of Global Variables + Types
int lattice[X][Y][Z];
int perc[X][Y][Z];
int seg[X][Y][Z];

int num_snakes = 0; //start with zero initial snakes - global to allow selection of snakes using this
int num_segments=0; //total number of snake segments. Used for percolation stats.

struct coord{ int x;  int y;  int z; };

struct snake_struct
{
     int head; //initial head at oroborus[0]
     int segs; //number of segments in snake
     int id;   //identity - number of snake.int >1
     struct coord oroborus[MAX_SEGMENTS];
               //array of locations of snake segments, loops around on itself - hence oroborus
}
 snakes[MAX_SNAKES];
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
void fill_snakes (float target_density);
void wriggle ();
void crawl(int x, int y, int z);
void percolate();  
int fit_snake (struct snake_struct *snake, int x, int y, int z, int length,
	                      int segment);
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
  return (int) (genrand_real2() * (double) max);
  }
  
