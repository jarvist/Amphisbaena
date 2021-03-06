

void
empty_lattice ()
{
  int x, y, z;
  for (x = 0; x < X; x++)
    for (y = 0; y < Y; y++)
      for (z = 0; z < Z; z++)
	lattice[x][y][z] = -1; //0;
   //-1 is unoccupied
   //0 and above indicate occupied, with the number being the ID of the snake
}

void
print_lattice ()
{
  int x, y=0, z = 0;
      for (z = Z-1; z >=0; z--)

    {
      fprintf (stderr,"\n");
   for (x = 0; x < X; x++)
	//for (z = 0; z < Z; z++)
	if (lattice[x][y][z] == -1)
	  fprintf (stderr,"%c[37m.%c[0m",27,27);
	else
	 {	    
//	  printf("o");
	  fprintf (stderr,"%c[%d",27,31+(lattice[x][y][z]%7));
	  if (perc[x][y][z]==1)		      
	      fprintf(stderr,";7");
fprintf(stderr,"m%c%c[0m",seg[x][y][z],27);
//printf("m%c%c[0m",(perc[x][y][z]%10)+'0',27);
	    
	 }
       
    }
  fprintf (stderr,"\n\n");
}

void calc_correlation()
{
 int x,y,z,tot,cx,cy,cz;

 tot=cx=cy=cz=0;

 for (z=0;z<Z-1;z++)
	 for (x=0;x<X-1;x++)
		 for (y=0;y<Y-1;y++)
			 if (lattice[x][y][z]!=-1)
			 {
				 tot++;
				 if (lattice[x][y][z+1]!=-1)
					 cz ++;
				 if (lattice[x][y+1][z]!=-1)
					 cy++;
				 if (lattice[x+1][y][z]!=-1)
					 cx++;
			 }
printf("Correlation in X: %f Y: %f Z: %f Sites: %d\n",(float)cx/(float)tot, (float)cy/(float)tot, (float)cz/(float)tot, tot);

}

void
  print_xmakemol ()
{
   int s,seg;
   char mole;
   printf("%d\n\n",num_segments);//total number of snake atoms
   for (s=0;s<num_snakes;s++)
     {
//	printf("molecule\n"); //each snake is new molecule
	if (perc[snakes[s].oroborus[0].x][snakes[s].oroborus[0].y][snakes[s].oroborus[0].z]==1) //electric snake?
	  mole='C';
	else
	  mole='O';
	
	for (seg=0;seg<snakes[s].segs;seg++)
	  printf("%c %.1f %.1f %.1f\n",
		 mole,
		 1.5*(float)snakes[s].oroborus[seg].x,
		 1.5*(float)snakes[s].oroborus[seg].y,
		 1.5*(float)snakes[s].oroborus[seg].z);
//	printf("\n");
	
     }   
}

void
    print_povray (char * name)
{
   
      int s,seg;
   float r,g,b;
      char mole;
      FILE *fo;
   fo=fopen(name,"w");

      fprintf(stderr,"%d\n\n",num_segments);//total number of snake atoms

      fprintf(fo,"camera {  location <-0.6,-0.6,0.5> up       <0,0,1> right    <1,0,0> sky <0,0,1> look_at  <0.5,0.5,0.5> angle    70 }         light_source { <-0.5,-2,0.5> color rgb <1,1,1>}  plane { z, -0.05 pigment { rgb 1 } }  #declare bbox = texture { pigment {rgb <0.4,0,0>} }   //marking X for substrate base cylinder {<0,0,0>,<1,1,0>,0.01 texture {bbox} } cylinder {<1,0,0>,<0,1,0>,0.01 texture {bbox} }\n");

      for (s=0;s<num_snakes;s++)
     {
	r=g=b=0;
	//      printf("molecule\n"); //each snake is new molecule
   //      if (perc[snakes[s].oroborus[0].x][snakes[s].oroborus[0].y][snakes[s].oroborus[0].z]==1) //electric snake?
     //      g=1;
       //  else
         //  r=1;

	g=2.0;
	r+=(float)snakes[s].oroborus[snakes[s].head].y/Z;

	fprintf(fo,"sphere{< %f, %f, %f>,%f texture{ pigment { rgb <%f,%f,%f>}}}\n",
	                        (float)snakes[s].oroborus[snakes[s].head].x/Z,
	                        (float)snakes[s].oroborus[snakes[s].head].y/Z,
	                        (float)snakes[s].oroborus[snakes[s].head].z/Z,
	                        0.25/Z,
	                        2.0,2.0,2.0);
        fprintf(fo,"sphere{< %f, %f, %f>,%f texture{ pigment { rgb <%f,%f,%f>}}}\n",
	                        (float)snakes[s].oroborus[(snakes[s].head-1)%snakes[s].segs].x/Z,
	                        (float)snakes[s].oroborus[(snakes[s].head-1)%snakes[s].segs].y/Z,
	                        (float)snakes[s].oroborus[(snakes[s].head-1)%snakes[s].segs].z/Z,
	                        0.20/Z,
	                        0.0,0.0,2.0);

         for (seg=snakes[s].head;seg<(snakes[s].segs+snakes[s].head-1);seg++)
	  {
//		  r=(float)(seg-snakes[s].head)/(float)snakes[s].segs;
//		  b=1-r;
//		  g=0;

           fprintf(fo,"cone{< %f, %f, %f>,%f <%f, %f, %f>,%f texture{ pigment { rgb <%f,%f,%f>}}}\n",
                  (float)snakes[s].oroborus[seg%snakes[s].segs].x/Z,
                  (float)snakes[s].oroborus[seg%snakes[s].segs].y/Z,
                  (float)snakes[s].oroborus[seg%snakes[s].segs].z/Z,
		  0.15/Z,
		  (float)snakes[s].oroborus[(seg+1)%snakes[s].segs].x/Z,
		  (float)snakes[s].oroborus[(seg+1)%snakes[s].segs].y/Z,
		  (float)snakes[s].oroborus[(seg+1)%snakes[s].segs].z/Z,
		  0.0,
		  r,g,b);
/*	     
	  fprintf(fo,"sphere{< %f, %f, %f>,%f texture{ pigment { rgb <%f,%f,%f>}}}\n",
		 (float)snakes[s].oroborus[seg%snakes[s].segs].x/Z,
		 (float)snakes[s].oroborus[seg%snakes[s].segs].y/Z,
		 (float)snakes[s].oroborus[seg%snakes[s].segs].z/Z,
		 0.15/Z,
		 r,g,b);
*/		 
 //      printf("\n");
	  } 
	
      }
   fclose(fo);
 } 


void
print_lattice_pnm ()
{
  printf ("P2\n%d %d\n%d\n", X, Y, num_snakes);

  int x, y, z = 0;
  for (x = 0; x < X; x++)
    {
      for (y = 0; y < Y; y++)
	{
	  if (lattice[x][y][z] == -1)
	    printf ("%d\t", num_snakes);
	  else
	    //printf("o");
	    printf ("%d\t", lattice[x][y][z]);
	}

      printf ("\n");
    }
  printf ("\n\n");
}

void save_lattice_file(char * name)
{
   int x,y,z;
   FILE *fo;
   fo=fopen(name,"w");
   fprintf(fo,"%d %d %d %d\n",X,Y,Z,num_snakes);
   
     for (x=0;x<X;x++)
       for (y=0;y<Y;y++)
	 {
	    fprintf(fo,"\n");
	 for (z=0;z<Z;z++)

	    fprintf (fo,"%d\t", lattice[x][y][z]);
	 }
   

     fprintf(fo,"\n\nMC: %d p: %f iE_int: %f iE_align: %f iE_PI: %f iE_SUBSTRATE: %f\n",
		     TOTAL_SLITHERS,DENSITY,
		     iE_INTERACTION,iE_ALIGNED,iE_PI,iE_SUBSTRATE);
   fclose(fo);
}

void load_lattice_file(char * name)
{
   int x,y,z,tmp;
   FILE *fo;
   fo=fopen(name,"r");
   fscanf(fo,"%d %d %d %d\n",&x,&y,&z,&num_snakes);
   if (x!=X || y!=Y || z!=Z )
     {
	
     fprintf(stderr,"CRASH! INPUT LATTICE PARAMS NOT COMPATIBLE. I READ %d %d %d, I EXPECTED %d %d %d!",x,y,z,X,Y,Z);
     exit(-1);
     }
   
     for (x=0;x<X;x++)
       for (y=0;y<Y;y++)
	 for (z=0;z<Z;z++)
	   fscanf(fo,"%d",&lattice[x][y][z]);
   fclose(fo);
}   

void turn_z_to_x()
{
 int x,y,z;
     for (x=0;x<X;x++)
	            for (y=0;y<Y;y++)
			             for (z=0;z<Z;z++)
					     perc[x][y][z]=lattice[x][z][y];
     for (x=0;x<X;x++)
	            for (y=0;y<Y;y++)
			             for (z=0;z<Z;z++)
					     lattice[x][y][z]=perc[x][y][z];
}

void print_lattice_pnm_file(int i)
{
 FILE *fo;
   char name[30]; //static buffers, security flaws, yada yada
     int x, y, z = 0;
   
   sprintf(name,"hydra_%d.pgm",i);
   fo=fopen(name,"w");   
   fprintf (fo,"P2\n%d %d\n%d\n", X, Y, num_snakes);
   

     for (y = 0; y < Y; y++)
     {
	      for (x = 0; x < X; x++)
	  {
	     
	               if (lattice[x][y][z] == -1)
	                   fprintf (fo,"%d\t", num_snakes);
	               else
	                   //printf("o");
	                   fprintf (fo,"%d\t", lattice[x][y][z]);
	  }
	
	
	      fprintf (fo,"\n");
     }
   
     fprintf (fo,"\n\n");
   fclose(fo);
   
}


void
print_snakes ()
{
  int i;
  for (i = 0; i < num_snakes; i++)
    {
      printf ("Id: %d \tSegments: %d \tHead: %d\n", snakes[i].id,
	      snakes[i].segs, snakes[i].head);

      printf ("Head Loc (x,y,z): %d %d %d\n",
	      snakes[i].oroborus[snakes[i].head].x,
	      snakes[i].oroborus[snakes[i].head].y,
	      snakes[i].oroborus[snakes[i].head].z);

    }
}

void
print_a_snake (int id)
{
  int i;
  printf ("Id: %d \tSegments: %d \tHead: %d\n", snakes[id].id,
	  snakes[id].segs, snakes[id].head);
  for (i = 0; i < snakes[id].segs; i++)
    printf ("Seg: %d (x,y,z) %d %d %d\n", i, snakes[id].oroborus[i].x,
	    snakes[id].oroborus[i].y, snakes[id].oroborus[i].z);
}



