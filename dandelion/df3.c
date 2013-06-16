/* generate_df3.c - generates .df3 volume density files for PovRay
 * Believe to be Public Domain & from http://astronomy.swin.edu.au/~pbourke/povray/df3/
 * Modified for Linear Combination of Atomic Orbitals visualisation by Jarvist Frost 2004
 */

#include <stdio.h>
#include <math.h>

#define SWAP_2(x) ( (((x) & 0xff) << 8) | ((unsigned short)(x) >> 8) )
#define SWAP_4(x) ( ((x) << 24) | \
		                (((x) << 8) & 0x00ff0000) | \
		                (((x) >> 8) & 0x0000ff00) | \
		                ((x) >> 24) )
#define FIX_SHORT(x) (*(unsigned short *)&(x) = SWAP_2(*(unsigned short *)&(x)))
#define FIX_INT(x)   (*(unsigned int *)&(x)   = SWAP_4(*(unsigned int *)&(x)))

int
generate_df3 (int count)
{
  int i, j, k, x, y, z, n;
   int ny=Y, nx=X, nz=Z;
   char name[20];
  short v;
  float themin = 1e32, themax = -1e32;
  FILE *fptr;
   
   fprintf(stderr,"Generating df3 file: %d\n",count);

  /* Calculate the bounds */
  for (i = 0; i < X; i++)
    {
      for (j = 0; j < Y; j++)
	{
	  for (k = 0; k < Z; k++)
	    {
	       if (lattice[i][j][k]<0.0) continue;
	       
	      if (themax < lattice[i][j][k])
		themax = lattice[i][j][k];
	      if (themin > lattice[i][j][k])
		themin = lattice[i][j][k];
	    }
	}
    }
  if (themin >= themax)
    {				/* There is no variation */
      themax = themin + 1;
      themin -= 1;
    }

//  printf ("themax: %f\tthemin: %f\n", themax, themin);
  /* Write it to a file */
   sprintf(name,"snakes_%.5d.df3",count);
  if ((fptr = fopen (name, "w")) == NULL)
    return (-1);

  fputc (nx >> 8, fptr);
  fputc (nx & 0xff, fptr);
  fputc (ny >> 8, fptr);
  fputc (ny & 0xff, fptr);
  fputc (nz >> 8, fptr);
  fputc (nz & 0xff, fptr);

  for (k = 0; k < nz; k++)
    {
      for (j = 0; j < ny; j++)
	{
	  for (i = 0; i < nx; i++)
	    {
	       if (lattice[i][j][k]>=0)
		 {
		    if (perc[i][j][k]>0)
		 //v = (int) ((float) 65535 * (lattice[i][j][k] - themin) / (themax - themin));
		      v=65535;
		    else
		      v=65535/4;
		 }
	       
	       else
		 v=0; //65535


	      FIX_SHORT (v);
	      fwrite (&v, 2, 1, fptr);
	    }
	}
    }
  fclose (fptr);
   fprintf(stderr,"Generated.\n");
}
