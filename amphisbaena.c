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

#include <math.h>
#include <stdio.h>
#include <time.h>
#include "mt19937ar-cok.c" //Mersenne twister random number generator; with Cok's speed optimisations
#include "amphisbaena.h" //simulation parameters & prototypes
#include "lattice_util.c" //lattice print + test functions
//#include "gorgophone.c" //ToF simulator
//#include "ran2.c"		//Numerical Recipes < 100'000'000 sequence portable generator
//#include "irbit2.c"		//Numerical Recipes random bit generorator

#include "df3.c"

main ()
{
  unsigned long long int i;
   int timestart;
   char name[50];
   
//  init_genrand(1); //seed MT random number 
  timestart=time(NULL); 
  init_genrand(timestart); //seed with current time in seconds since 1970 Unix Epoch
   
  empty_lattice (); //clears lattice
   fprintf(stderr,"Empty Lattice E: %f\n",lattice_energy());
  fill_snakes (DENSITY);
//  print_lattice();
//   sleep(1);
  //print_snakes();
  //print_a_snake(0);
   
  fprintf (stderr, "Snakes filled\nNum Snakes: %d Lattice Points: %d\n",
	   num_snakes, X * Y * Z);
//  fprintf(stderr,"Initial Snake Filled Lattice Energy: %f\n",lattice_energy());   
   
//   printf("#X,Y,Z: %d %d %d\n#Density: %f\n#Beta: %f\n#l_0: %d\n#sigma: %f\n#iE[1,1]: %f\n#Snakes: %d\n#Slithers: %d\n",
//	  X,Y,Z,DENSITY,beta,l_0,sigma, iE[1][1], num_snakes,TOTAL_SLITHERS);

   for (i = 0; i < (long long) TOTAL_SLITHERS; i++)
    {
      wriggle (); //this function is the heart of the simulation - executed millions of times per second
       //might be worth unrolling this directly into this loop rather than have it as a subroutine...
       
       if (i%SLITHER_PRINT==0) //display information for monitoring the progression of the simulation
	 {
	 fprintf(stderr,"Slithers: %lld Snakes: %d Lattice Points: %d \n",i,num_snakes,X*Y*Z);
//	    lattice_energy();
	    percolate();

//	    generate_df3(i/SLITHER_PRINT);
//
	    	    print_lattice();
	    
//	    print_povray ();
//	    fprintf(stderr,"Lattice Energy: %f\n",lattice_energy());
//        usleep(1000000);	    
//	print_lattice_pnm_file((int)(beta*1000.00));    
	 }
    }
   
  percolate();   //tests for percolation - and sets electric snakes as 'live' if part of percolating cluster
                 //its a recursive algorithm though - so will crash if too many lattice sites!
//  print_lattice();
//  print_xmakemol();
   print_povray("snakes.pov");
//  print_lattice_pnm_file((int)(beta*1000.00));    
   print_lattice_pnm_file(1);

  fprintf (stderr, "Num Snakes: %d Lattice Points: %d\n", num_snakes,
	   X * Y * Z);
//   fprintf(stderr,"\nStats:\n\tRandom Numbers used: %lld\n\t Steps: %lld\n\tRand per Step: %f\n",
//	   rand_count,i,((double)rand_count/(double)i) );

//   scanf("%s",&name);
//   save_lattice_file(name);   
   
//   gorgophone();
//   printf("#Time Taken: %d seconds\n",time(NULL)-timestart);
}

float lattice_energy() //calculate overall energy of lattice
{
  int x, y, z,m1,m2,i,dx,dy,dz,mat1,mat2;
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

void percolate()
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
}

void crawl(int x, int y, int z)
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
}

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
void
wriggle ()
{
  int s, head, heads, tail;
  int dx, dy, dz;
  int x, y, z;
  double dE = 0.0;
  struct coord h, t;
  int mat, i;

  s = rand_int (num_snakes);
  //choose a snake at random

  heads = -rand_int (2);
  //change this to use dedicated bit random generator
  // i.e.heads = 0, -1 with equal prob.

  dx = dy = dz = 0;
  switch (rand_int (6) + 1)
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

//  fprintf(stderr,"\nWriggle: s: %d\t head: %d\theads: %d\t (dx,dy,dz): %d %d %d\n", s, snakes[s].head, heads, dx, dy, dz);

  head = snakes[s].head + heads; //Either choose head1, where the .head points to, 
                                 //or the tail/head2 to act as head this time

  if (head >= snakes[s].segs) //The Oroborus datatype is circular, so wrap around...
    head = 0;
  if (head < 0)
    head = snakes[s].segs - 1;

//  printf("head: %d\t tail:%d\n", head, tail);
  x = snakes[s].oroborus[head].x;
  y = snakes[s].oroborus[head].y;
  z = snakes[s].oroborus[head].z;

  if (x + dx < 0 || x + dx >= X
      || y + dy < 0 || y + dy >= Y 
      || z + dz < 0 || z + dz >= Z)
    return;
  //not allowed lattice location !
  //The majority of calls to this function will have returned at this point, 
  // so the above code is stripped to the bare essentials to reach here as quickly as possible
  
  //Now we can calculate other variables necessary to actually move the snake - like the tail location;

  if (heads < 0)
    tail = head + 1; //tail location before head i.e.going forwards
  else
    tail = head - 1; //tail location after head i.e.going backwards

  if (tail >= snakes[s].segs) //The Oroborus datatype is circular, so wrap around...
    tail = 0;
  if (tail < 0)
    tail = snakes[s].segs - 1;
   
   
  lattice[snakes[s].oroborus[tail].x]
    [snakes[s].oroborus[tail].y][snakes[s].oroborus[tail].z] = -1;
  //i.e.Tail moves out of the way first
  // This is so that the 'head' can follow the tail around in a closed circle
  // I this allows enclosed / boxed-in snakes to 'cycle' around and then escape

  if (lattice[x + dx][y + dy][z + dz] == -1)
    //gap where head is being forced
    // i.e.move is PHYSICALLY possible now need to see whether ENERGETICALLY FAVOURABLE
    {
//      fprintf(stderr,"Space at: (x,y,z) %d %d %d\n", x + dx, y + dy, z + dz);
      //printf("c");

      h.x = x + dx; //h is new head location
      h.y = y + dy;
      h.z = z + dz;
      
      //maybe just copy the pointer to coord ?
      t.x = snakes[s].oroborus[tail].x; //t is tail location
      t.y = snakes[s].oroborus[tail].y;
      t.z = snakes[s].oroborus[tail].z;

      dE = 0.0;
      for (i = 1; i <= 6; i++)
	//all 6 directions
	{
	  dx = dy = dz = 0;
	  switch (i)
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

	   //NB: Energy change as a result of where the new head is going affect the overall snake positively
	   //Enthalpy of where the tail used to be affects the overall energy negatively
	  if (h.x + dx >= 0 && h.x + dx < X	// if prospective interaction site within lattice
	      && h.y + dy >= 0 && h.y + dy < Y
	      && h.z + dz >= 0 && h.z + dz < Z)
	    {
	      if (lattice[h.x + dx][h.y + dy][h.z + dz] == -1)
		mat = 0;  //if empty lattice site material 0[non snake]
	      else
		mat = 1;  //if snake material type 1[snake type]
	      dE += iE[1][mat] - iE[0][mat];

	      if (lattice[h.x + dx][h.y + dy][h.z + dz] == snakes[s].id) //if touching ourselves
		dE += iE[1][mat] - iE[0][mat];
	      //count energy change twice - once for head
	      //segment interaction and once for segment head interaction
	    }

	  if (t.x + dx >= 0 && t.x + dx < X	// if prospective interaction site within lattice
	      && t.y + dy >= 0 && t.y + dy < Y
	      && t.z + dz >= 0 && t.z + dz < Z)
	    {

	      if (lattice[t.x + dx][t.y + dy][t.z + dz] == -1)
		mat = 0;   //if empty lattice site, material 0[non snake]
	      else
		mat = 1;   //if snake material type 1[snake type] 
	      dE += iE[0][mat] - iE[1][mat];

	      if (lattice[t.x + dx][t.y + dy][t.z + dz] == snakes[s].id)
		dE += iE[0][mat] - iE[1][mat];
	      //if touching ourselves count energy change twice - 
	      //once for head-segment interaction
	      // and once for segment-head interaction}
	    }
	   
	   //SUBSTRATE substrate interactions
	   if (t.z+dz<0) //added JMF 22-1-07
	     dE-=500000.0;
	   if (h.z+dz<0)
	     dE+=500000.0;
	}

//      fprintf(stderr, "dE: %f\t", dE);

      if (dE > 0.0 || exp(dE * beta) > rand_float ())
	//if dE exothermic, reaction progresses automatically
	//or if sufficient boltzmann energy to drive endothermic reaction
	{
//	   fprintf(stderr,"Going! (x,y,z) %d %d %d\n",h.x,h.y,h.z);
	   
	  lattice[h.x][h.y][h.z] = snakes[s].id;
	  //mark new head location on lattice
//	   print_lattice();
	   
	  if (heads < 0)
	    {
	      //forwards 
	      snakes[s].head = (tail + 1) % snakes[s].segs;
	    }
	  else
	    snakes[s].head = tail;

	  snakes[s].oroborus[tail].x = h.x;	//tail now becomes new head...
	  snakes[s].oroborus[tail].y = h.y;
	  snakes[s].oroborus[tail].z = h.z;
	}
      else
	{
//	   fprintf(stderr,"Holding!\n");
	  lattice[snakes[s].oroborus[tail].x]
	    [snakes[s].oroborus[tail].y]
	    [snakes[s].oroborus[tail].z] = snakes[s].id;
	  //i.e.Tail put back in location
	}

    }
  else
    {
      lattice[snakes[s].oroborus[tail].x]
	[snakes[s].oroborus[tail].y]
	[snakes[s].oroborus[tail].z] = snakes[s].id;
      //i.e.Tail put back in location


      // printf("f");
      //printf("No space at: (x,y,z) %d %d %d\n", x + dx, y + dy, z + dz);
    }

}

//the gaussian function
double
gauss (double x, double x_0, double mysigma)
{
  return exp (-0.5 * (x - x_0) * (x - x_0) / mysigma / mysigma);
}

void
place_snake (struct snake_struct *snake, int length)
{
  int x, y, z, flag = 0;

  do
    {
      x = rand_int (X);
      y = rand_int (Y);
      z = rand_int (Z);
      //printf("Trying to place Snake: Length: %d\tat\atx: %d\ty: %d\tz: %d\n", length, x, y, z);
      flag = fit_snake (snake, x, y, z, length, 0);
    }
  while (flag == 0); //yes, this infinite loops if you ask it to do something impossible
  //printf("Succesfully Placed! Length:%d \tat\atx: %d\ty: %d\tz: %d\n", length, x, y, z);
  //print_lattice();

  snake->segs = length;
   num_segments+=length; //update total number of snake segments counter
  snake->head = 0;
  //setup snake structure with necessary info no that its created

}


void
fill_snakes (float target_density)
{
  //ok first we generate a table of probability of creating a snake of length l
  // drawn from the Gaussian distribution.
  // $p(l) = exp(-\frac {(l - l_0) ^ 2} {2 \ sigma ^ 2} $

  double p_l[MAX_SEGMENTS], s_l[MAX_SEGMENTS],
    p_total = 0.0, p_total_length = 0.0, weighting, tmprand;
  int l, segments = 0, target_segments = 20, size;

  target_segments =
    (int) (target_density * (float) (X * Y * Z));

  for (l = 1; l < MAX_SEGMENTS; l++)
    //minimum snake size is 2 segments
    {
      s_l[l] = 0;
      //make no snakes of this size
      p_l[l] = gauss ((double) l, l_0, sigma);
      //fill p_l with probability / length distribution

      p_total += p_l[l];
      //total 'probability' within our p range
      p_total_length += p_l[l] * (double) l;
      //total 'probability' times length.'. expectation

      // printf("%d\t%f\t%f\t%f\n", l, p_l[l], p_total, p_total_length);
    }
  //weighting = (target_density * (float) X * (float) Y * (float) Z) / p_total_length;
  //printf("Weighting: %f\t Target: %f\n", weighting, target_density);

  while (segments < target_segments)
    {
      size = rand_int (MAX_SEGMENTS - 2) + 1;
      //printf("Size: %d\n", size);
      tmprand = rand_float ();
      if (tmprand < (p_l[size] / p_total))
	{

	  //printf("Snake Choosen; size: %d\tp_l:%f\ttmprand:%f\n", size, p_l[size] / p_total, tmprand);
	  segments += size;
	  //segments++;
	  s_l[size]++;
	}
    }

  for (l = MAX_SEGMENTS - 1; l > 1; l--)
    //start filling in lattice with longest snakes
    while (s_l[l]-- > 0)
      //while number of snakes of this size is > 0
      {
	fprintf (stderr, "Placing Snake %d Size %d\n", num_snakes, l);
	if (num_snakes >= MAX_SNAKES)
	  {
	    fprintf (stderr, "Number snakes > MAX_snakes. Dieing!\n");
	    //exit(-1);
	  }
	snakes[num_snakes].id = num_snakes;
	//setup id
	snakes[num_snakes].segs = l;
	place_snake (&snakes[num_snakes], l);
	//randomly place snake of length l
	num_snakes++;
      }
}

int
fit_snake (struct snake_struct *snake, int x, int y, int z,
	   int length, int segment)
{
  int flag = 0;
  //returns 'snake fitted' flag down the tree of recursion
  int order[6] = { 1, 2, 3, 4, 5, 6 }; //order we do the directions in 
   int mix, site1, site2, dx, dy, dz, tmp;

  //printf("fit_snake: x:%d\ty:%d\tz:%d\tlength:%d\n", x, y, z, length);

  if (x < 0 || x >= X || y < 0 || y >= Y || z < 0 || z >= Z)
    return 0;
  //not allowed lattice location !
  //printf("Boundaries ok.\n");

  //printf("L: %d %d %d = %d\n", x, y, z, lattice[x][y][z]);

  //this location already occupied
  if (lattice[x][y][z] != -1)
    return 0;

  snake->oroborus[segment].x = x;
  //put location of this segment in oroborus array in snake
  snake->oroborus[segment].y = y;
  snake->oroborus[segment].z = z;
  lattice[x][y][z] = snake->id;
  //mark lattice location as occupied

  /*
   * printf("Segment: %d x,y,z %d %d %d =
   * %d\n",segment, snake->oroborus[segment].x,
   * snake->oroborus[segment].y,
   * snake->oroborus[segment].z,
   * lattice[x][y][z]);
   */

  //snake succesfully fitted !
  if (length - 1 == segment)
    return 1;
  //return flag to travel down recursion tree...
   
  // mix up order of trying different directions
  // try possible locations for next segment
  for (mix = 0; mix < 20; mix++)
    {

      site1 = rand_int (6);
      site2 = rand_int (6);
      tmp = order[site2];
      order[site2] = order[site1];
      order[site1] = tmp;
    }

  /*
   * for (mix=0;mix<6;mix++) //debugging
   * printer for order mixer
   * printf("%d\t",order[mix]); printf("\n");
   */

  for (mix = 0; mix < 6; mix++)
    {
      dx = dy = dz = 0;
      switch (order[mix])
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

      if (flag =
	  fit_snake (snake, x + dx, y + dy, z + dz, length, segment + 1) == 1)
	break;
    }

  //this location didn 't work out...
  if (flag == 0)
    lattice[x][y][z] = -1;
  //wipe lattice location - oroborus[segment] will get overwritten anyway

  return flag;
}
