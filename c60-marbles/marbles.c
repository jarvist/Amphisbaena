/* MARBLES - a C60 crystal sintering MC program; based on:
 * 
 * Amphisbaena:- (pronounced am.fis.BEEN.uh)
 * A two-headed Slithering Snake implementation designed to simulate morphologies of Organic Solar Cells.
 * Code by Jarvist Frost, based on algorithm and physical setup by Dr Jenny Nelson & Felix Rickerman.
 * File begun June 16th 2005
 * 
 * Molecular Electronic Materials and Devices
 * Experimental Solid State Physics
 * Blackett Laboratory
 * Imperial College, London
 */

#include <math.h>
#include <stdio.h>
#include <time.h>
//#include "mt19937ar-cok.c" //Mersenne twister random number generator; with Cok's speed optimisations
#include "SFMT/SFMT.c"
#include "marbles.h" //simulation parameters & prototypes
#include "lattice_util.c" //lattice print + test functions
//#include "gorgophone.c" //ToF simulator
//#include "ran2.c"		//Numerical Recipes < 100'000'000 sequence portable generator
//#include "irbit2.c"		//Numerical Recipes random bit generorator

//#include "df3.c" //requires perc lattice

int count_balls();

int main ()
{
   int i;
   int timestart,timenow,secleft;
//   char name[50];
   
   empty_lattice (); //clears lattice   
   
//  init_genrand(1); //seed MT random number 
  timestart=time(NULL); 
  init_gen_rand(timestart); //seed with current time in seconds since 1970 Unix Epoch
   

   fprintf(stderr,"Empty Lattice E: %f\n",lattice_energy());

   placemarble(); //MUST DO THIS AT LEAST ONCE BEFORE WRIGGLE!
   
   for (i = 0; i < MAX_SNAKES; i++)
    {
      wriggle (MC_BETWEEN_DROPS); //this function is the heart of the simulation - executed millions of times per second
       //might be worth unrolling this directly into this loop rather than have it as a subroutine...
      placemarble();
       
/*       print_lattice();
       
       if (count_balls()>num_snakes)
	 {
	    fprintf(stderr,"Oh no! Too many balls %d vs. %d...This ends now.\n",count_balls(),num_snakes);
	    print_lattice();
	    print_xmakemol("marbles.xyz");
	    return(-1);
	 }
  */     
       
       if (i>0 && i%SLITHER_PRINT==0) //display information for monitoring the progression of the simulation
	 {

//	    generate_df3(i/SLITHER_PRINT);
//

	    fprintf(stderr,"Slithers: %dM Of: %dM Marbles: %d Of: %d\n",(((int)i*(MC_BETWEEN_DROPS/1000))/1000),(MAX_SNAKES*(MC_BETWEEN_DROPS/1000))/1000,num_snakes,MAX_SNAKES);  

//	    print_lattice();
//	    print_xmakemol("marbles.xyz");
	    
	    timenow=time(NULL);
	    fprintf(stderr,"    %ds; Avg. %fM MC Moves per second; ",timenow-timestart,(double)((i*(MC_BETWEEN_DROPS/1000))/1000)/(double)(timenow-timestart));

	    secleft=(float)(MAX_SNAKES-i)*((float)(timenow-timestart)/(float)i); //we predict there will be this many seconds to finish
	    
	    fprintf(stderr,"Finish in: %dd %dhr %dm %ds\n",
		    (int)(secleft/(60*60*24)),
		    (int)(secleft/(60*60)%24),
		    (int)(secleft/(60))%(60),
		     (int)(secleft)%(60) );
	    
//	    timeend=(int)((float)MAX_SNAKES*((float)(timenow-timestart)/(float)i));	    
//	    fprintf(stderr,"At: %s\n",asctime(localtime((time_t) &timeend)));
//	    print_pdb ("traj.pdb", i/1000000);
		
	    
//	    print_povray ();
//	    fprintf(stderr,"Lattice Energy: %f\n",lattice_energy());
//        usleep(1000000);	    
//	print_lattice_pnm_file((int)(beta*1000.00));    
	 }
    }
//   fprintf(stderr,"Going for glory: 200M equil. MC steps...\n");
//   wriggle(200000000);
   
//  percolate();   //tests for percolation - and sets electric snakes as 'live' if part of percolating cluster
                 //its a recursive algorithm though - so will crash if too many lattice sites!
//  print_lattice();
 
  print_xmakemol("marbles.xyz");
//   print_povray("snakes.pov");
//  print_lattice_pnm_file((int)(beta*1000.00));    
//  
//  fprintf (stderr, "Num Snakes: %d Lattice Points: %d\n", num_snakes,
//	   X * Y * Z);
//   fprintf(stderr,"\nStats:\n\tRandom Numbers used: %lld\n\t Steps: %lld\n\tRand per Step: %f\n",
//	   rand_count,i,((double)rand_count/(double)i) );

//   scanf("%s",&name);
//   save_lattice_file(name);   
   
//   gorgophone();
//   printf("#Time Taken: %d seconds\n",time(NULL)-timestart);
}

int count_balls()
{
 int count=0,x,y,z;
   
   for (x=0;x<X;x++)
     for (y=0;y<Y;y++)
       for (z=0;z<Z;z++)
	 if (lattice[x][y][z]!=-1)
	   count++;
   return(count);
   
}

float lattice_energy() //calculate overall energy of lattice
{
  int x, y, z,i,dx,dy,dz,mat1,mat2; //m1,m2
   double E=0.0;
   int clumps=0;
   double clumpiness;
   
  for (x = 0; x < X; x++)
    for (y = 0; y < Y; y++)
      for (z = 0; z < Z; z++)
	 {
      for (i = 1; i <= 6; i++)
	              //all 6 directions
	      {
		 
		           dx = dy = dz = 0;
		           switch (i)
		               //choose which of 6 directions to attempt wriggle in
		   {
		      
		      
		    case 1:
		                    dx = 1;
		                    break;
		    case 2:
		                    dx = -1;
		                    break;
		    case 3:
		                    dy = 1;
		                    break;
		    case 4:
		                    dy = -1;
		                    break;
		    case 5:
		                    dz = 1;
		                    break;
		    case 6:
		                    dz = -1;
		                    break;
		   }
		 
		 
		           if (x + dx >= 0 && x + dx < X     // if prospective interaction site within lattice
			                     && y + dy >= 0 && y + dy < Y
			                     && z + dz >= 0 && z + dz < Z)
		   {		      
		      if (lattice[x][y][z]==-1)
			mat1=0;
		      else
			mat1=1;
		      if (lattice[x+dx][y+dy][z+dz]==-1)
			mat2=0;
		      else
			mat2=1;
		    E+=iE[mat1][mat2];
		      
		      if (mat1==mat2) clumps++; //tot up number of similar adjacent sites
//		      fprintf(stderr,"x: %d y: %d z: %d dx: %d dy: %d dz: %d E: %f\n",x,y,z,dx,dy,dz,E);
		   }
		   }
	 }
   //why does the following look so strange?
   //Because you need to divide by the number of three dimensional 'prison bars' that are formed between adjacent 
   //sites in different axes.
   clumpiness=(double)clumps/(2.0*(double)(
			     (X-1)*Y*Z
			     + X*(Y-1)*Z
			     + X*Y*(Z-1)
			     ));
   fprintf(stderr," E: %f Clumpiness: %f\n",E,clumpiness);
 return(E);  
}

int electricsegments;
int percolation;

/*void percolate()
{
   int x, y, z;
   
   for (x = 0; x < X; x++) //reset percolation / electrification lattice
         for (y = 0; y < Y; y++)
             for (z = 0; z < Z; z++)
	         perc[x][y][z] = 0;
   electricsegments=0;
   percolation=0;
   
   for (x = 0; x < X; x++) //electrifiy y=Y-1 plane
     for (y = 0; y < Y; y++)
       crawl(x,y,0); //if snake present
   
   fprintf(stderr,"Density: %f Electricsegments: %d Percolation: %d Fraction Electric Snakes: %f\n",DENSITY,electricsegments,percolation,(float)electricsegments/(float)num_segments);
}*/


void placemarble()
{
   int x,y,z;

 do 
     {
// x=rand_int(X/2)+X/4; y=rand_int(Y/2)+Y/4; //across aperture in core
 x=rand_int(X); y=rand_int(Y); //across whole substrate
 z=Z-2;
     } while (lattice[x][y][z]>=0);
   
// fprintf(stderr,"d");
   
 while (lattice[x][y][z]==-1 && z>0)
     z--;
 z++;
 if (lattice[x][y][z]>=0)
     fprintf(stderr,"Oh no! I've lost my marbles!\n");
   
// printf("X,Y,Z height of new marble: %d %d %d\n",x,y,z);
   
 marbles[num_snakes].x=x; marbles[num_snakes].y=y; marbles[num_snakes].z=z;
 lattice[x][y][z]=num_snakes;
   
 num_snakes++;
// fprintf(stderr, "Marble placed: %d %d %d num_snakes %d\n",x,y,z,num_snakes);  
}


/*void crawl(int x, int y, int z)
{
  if (x<0 || x>=X ||
      y<0 || y>=Y ||
      z<0 || z>=Z ||
      lattice[x][y][z]==-1 || 
      perc[x][y][z]==1) //no snake OR already been here by alternate route...
     return;
   
  perc[x][y][z]=1; 
  electricsegments++;
   
   if (z==0) //crossed lattice
     percolation=1;
   
   crawl(x+1,y,z);
   crawl(x-1,y,z);
   crawl(x,y+1,z);
   crawl(x,y-1,z);
   crawl(x,y,z+1);
   crawl(x,y,z-1);
}*/

/* The following function is the very core of the simulation.
 * When called, it does the following:-
 *  Pick a random snake
 *  Choose which end will be the 'head'
 *  Choose a random direction for the head to 'move in'
 *  If that location is not free, return immediately.
 *  Otherwise, calculated the Energy of the snakes in the new 'shifted' position versus the current configuration.
 *  Compare this energy shift with an energy randomly pulled from a Boltzmann distribution of temperature/energy
 *  If energetically favourable, move the snake.
 *  Otherwise, put the tail back in place
 *  Return!
 */
const int lp [2][12][3]= //look up tables ROXOR
  //lattice parameters: Given loc. [x][y][z], which are the offsets to the 12 nearest neighbours
  //considering the magic chequer board compression into a cubic matrix we use  (that's what the second half set of params are) :)
  //This used to be generated with a 'switch' twice within wriggle. Lookup tables should be much faster.
{
     {
	  {1,0,0}, //same z
	  {-1,0,0},  
	  {0,1,0},
	  {0,-1,0},
	
	  {0,0,1}, //+1 z
	  {0,1,1},
	  {1,0,1},
	  {1,1,1},
	
	  {0,0,-1},
	  {0,1,-1},
	  {1,0,-1},
	  {1,1,-1},
     },
   
     {
	  {1,0,0}, //same z
	  {-1,0,0},  
	  {0,1,0},
	  {0,-1,0},
	
	  {-1,-1,1}, //+1 z
	  {-1,0,1},
	  {0,-1,1},
	  {0,0,1},
	
	  {-1,-1,-1},
	  {-1,0,-1},
	  {0,-1,-1},
	  {0,0,-1},
	
     }
   
};


void
wriggle (int wriggles)
{
  int s;
  int dx, dy, dz;
  int x, y, z;
  double dE = 0.0;
  struct coord h, t;
  int mat, i,dir,count;

   for(count=0;count<wriggles;count++)
     {
	
  s = rand_int (num_snakes);
  //choose a snake at random

//  dx = dy = dz = 0;
/*  dx=rand_int(3)-1; //choose one of 9 pos. nearest on cube
  dy=rand_int(3)-1;
  dz=rand_int(3)-1;
*/


   dir=rand_int(12); //direction c.f. lookup table (lp) we are planning to move
   
   x=marbles[s].x; y=marbles[s].y; z=marbles[s].z;
   
   dx=lp[z%2][dir][0]; dy=lp[z%2][dir][1]; dz=lp[z%2][dir][2]; //such in values from lookup table 

   h.x=(x+dx + X)%X; h.y=(y+dy + Y)%Y; h.z=z+dz;

  //not allowed into edge of space 
   
//  fprintf(stderr,"I CAN HAS CHEESBISQUIT? (x,y,z) %d %d %d (dx,dy,dz) %d %d %d\n",x,y,z,dx,dy,dz);
  
  if (  z + dz <= 0 || z + dz >= Z)
    continue;
  //not allowed lattice location !
  //The majority of calls to this function will have returned at this point, 
  // so the above code is stripped to the bare essentials to reach here as quickly as possible
   
//   fprintf(stderr,"Try: (x,y,z) %d %d %d\n",x+dx,y+dy,z+dz);
   
  if (lattice[h.x][h.y][h.z] == -1)
    //gap where head is being forced
    // i.e.move is PHYSICALLY possible now need to see whether ENERGETICALLY FAVOURABLE
    {
//      fprintf(stderr,"Space at: (x,y,z) %d %d %d\n", x + dx, y + dy, z + dz);
      //printf("c");
       t.x=x; t.y=y; t.z=z;
       
       lattice[t.x][t.y][t.z]=-1; //LATTICE SET EMPTY
      
      dE = 0.0;
      for (i = 0; i < 12; i++)
	//all 6 directions
	{
	   
           dx=lp[z%2][i][0]; dy=lp[z%2][i][1]; dz=lp[z%2][i][2]; //suck in values from lookup table 
	   

	   //NB: Energy change as a result of where the new head is going affect the overall snake positively
	   //Enthalpy of where the tail used to be affects the overall energy negatively



//	      if (lattice[(h.x + dx + X)%X][(h.y + dy + Y)%Y][h.z + dz] == -1)
//		mat = 0;  //if empty lattice site material 0[non snake]
//	      else
//		mat = 1;  //if snake material type 1[snake type]

	      mat=lattice[(h.x + dx + X)%X][(h.y + dy + Y)%Y][h.z + dz] != -1;
//	      mat=(lattice[(h.x + dx + X)%X][(h.y + dy + Y)%Y][h.z + dz]+1)&1; //one line of arithmatic to replace 'if' statement selection rule above
	      dE += iE[1][mat] - iE[0][mat];

//	      if (lattice[(t.x + dx + X)%X][(t.y + dy +Y)%Y][t.z + dz] == -1)
//		mat = 0;   //if empty lattice site, material 0[non snake]
//	      else
//		mat = 1;   //if snake material type 1[snake type] 
//	      mat=(lattice[(t.x + dx + X)%X][(t.y + dy +Y)%Y][t.z + dz]+1)&1;

	      mat=lattice[(t.x + dx + X)%X][(t.y + dy +Y)%Y][t.z + dz] != -1;

	      dE += iE[0][mat] - iE[1][mat];

	   
	   //SUBSTRATE substrate interactions
	   if (t.z+dz==0) //added JMF 22-1-07
	     dE-=SUBSTRATE_E;
	   if (h.z+dz==0)
	     dE+=SUBSTRATE_E;
	}

//      fprintf(stderr, "dE: %f\t", dE);
//      fprintf(stderr,"From: (x,y,z) %d %d %d to %d %d %d with dE: %f", t.x,t.y,t.z,h.x,h.y,h.z,dE);
//      
      if (dE > 0.0 || exp(dE * beta) > rand_float ())
	//if dE exothermic, reaction progresses automatically
	//or if sufficient boltzmann energy to drive endothermic reaction
	{
//	   fprintf(stderr,"Going! (x,y,z) %d %d %d\n",h.x,h.y,h.z);
	   
	  lattice[h.x][h.y][h.z] = s;
//	  lattice[t.x][t.y][t.z] = -1;
	  marbles[s].x=h.x;
	  marbles[s].y=h.y;
	  marbles[s].z=h.z;
	  //mark new head location on lattice
//	   print_lattice();
//	   
//	   fprintf(stderr,"  I go...\n");
	   
	}
      else
	{
//	   fprintf(stderr,"Holding!\n");
	  //i.e.Tail put back in location
	      lattice[t.x][t.y][t.z]=s;
//		   fprintf(stderr,"  I stay...\n");
	}

    }
   
//  else
//    {
//      lattice[snakes[s].oroborus[tail].x]
//	[snakes[s].oroborus[tail].y]
//	[snakes[s].oroborus[tail].z] = snakes[s].id;
      //i.e.Tail put back in location

//             lattice[t.x][t.y][t.z]=s;

      // printf("f");
      //printf("No space at: (x,y,z) %d %d %d\n", x + dx, y + dy, z + dz);
//    }
     }
   
}

//the gaussian function
double
gauss (double x, double x_0, double mysigma)
{
  return exp (-0.5 * (x - x_0) * (x - x_0) / mysigma / mysigma);
}

