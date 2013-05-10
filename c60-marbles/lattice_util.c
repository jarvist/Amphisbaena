int pdb_model=1;

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
  int x, y=0, z = 0,c=0;
   
      for (y = 0; y <Y; y++)
     {
	
      printf ("\n");
     for (x = 0; x < X; x++)
    {
       c=0;
	for (z = 0; z < Z; z++)
	   {if (lattice[x][y][z]>=0) c++;}
       
	if (c==0)
	  printf ("%c[37m.%c[0m",27,27);
	else
	 {	    
//	  printf("o");
	  printf ("%c[%d",27,31+(c%7));
//	  if (perc[x][y][z]>0)		      
//	      printf(";7");
	  printf("m%c%c[0m",(c%26)+'A'-1,27);
	    
	 }
       
       printf (" ");
       }
     }
   
  printf ("\n\n");
}

void
  print_xmakemol (char * name)
{
   int s;
   char mole;
   FILE *fo;
   fo=fopen(name,"w");
   
   //stash a bit of useful info in this file...
   fprintf(fo,"#Beta: %f NUM_SNAKES: %d (X Y Z): (%d %d %d) MC_BETWEEN_DROPS: %d\n",beta,MAX_SNAKES,X,Y,Z,MC_BETWEEN_DROPS);
   fprintf(fo,"#iE: %f SUBSTRATE_E: %f\n",iE[1][1],SUBSTRATE_E);
   
   fprintf(fo,"%d\n\n",num_snakes);//total number of snake atoms
   for (s=0;s<num_snakes;s++)
     {
//	printf("molecule\n"); //each snake is new molecule
//	if (perc[snakes[s].oroborus[0].x][snakes[s].oroborus[0].y][snakes[s].oroborus[0].z]==1) //electric snake?
//	  mole='C';
//	else
	  mole='O';
	
	if (marbles[s].z%2==1)
	  fprintf(fo,"%c %.1f %.1f %.1f\n",
		 mole,
		 5.0*(float)marbles[s].x,
		 5.0*(float)marbles[s].y,
		 3.16*(float)marbles[s].z);
	else //offset
	  fprintf(fo,"%c %.1f %.1f %.1f\n",
		 mole,
		 2.5+(5.0*(float)marbles[s].x),
		 2.5+(5.0*(float)marbles[s].y),
		 3.16*(float)marbles[s].z);
	  
//	printf("\n");
	
     }   
   fclose(fo);
}

void
  print_pdb (char * name, int t)
{
   //This kind of works;
   //the idea was automagically output trajectories
   //But I don't really care anymore. Jarv 11-06-07
   int s;
   char mole;
   FILE *fo;
   fo=fopen(name,"a");
   
//   fprintf(fo,"%d\n\n",num_snakes);//total number of snake atoms
   fprintf(fo,"TITLE C60_MARBLES t=%d.00\n",t);
   fprintf(fo,"CRYST1 %f %f %f 90.0 90.0 90.0 P 1           1\n",5.0*(float)X,5.0*(float)Y,5.0*(float)Z);
   fprintf(fo,"MODEL\t%d\n",pdb_model++);
   for (s=0;s<num_snakes;s++)
     {
//	printf("molecule\n"); //each snake is new molecule
//	if (perc[snakes[s].oroborus[0].x][snakes[s].oroborus[0].y][snakes[s].oroborus[0].z]==1) //electric snake?
//	  mole='C';
//	else
	  mole='C';
	
	if (marbles[s].z%2==1)
	  fprintf(fo,"HETATM %4d  %c%c%c CBB     1     %07.3f %07.3f %07.3f  1.00  0.00           C\n",
		 s+1,
		 mole,'A'+s%26,'A'+(s/26)%26,
		 5.0*(float)marbles[s].x+0.001,
		 5.0*(float)marbles[s].y+0.001,
		 3.16*(float)marbles[s].z+0.001
		  );
	else //offset
	  fprintf(fo,"HETATM %4d  %c%c%c CBB     1     %07.3f %07.3f %07.3f  1.00  0.00           C\n",
		 s+1,
		 mole,'A'+s%26,'A'+(s/26)%26,
		 2.5+(5.0*(float)marbles[s].x+0.001),
		 2.5+(5.0*(float)marbles[s].y+0.001),
		 3.16*(float)marbles[s].z+0.001
		  );
	  
//	printf("\n");
	
     }   
   
  fprintf(fo,"TER\nENDMDL\n");
  fclose(fo);
}


/*void
    print_povray (char * name)
{
   
      int s,seg;
   float r,g,b;
      char mole;
      FILE *fo;
   fo=fopen(name,"w");

      fprintf(stderr,"%d\n\n",num_segments);//total number of snake atoms
      for (s=0;s<num_snakes;s++)
     {
	r=g=b=0;
	//      printf("molecule\n"); //each snake is new molecule
         if (perc[snakes[s].oroborus[0].x][snakes[s].oroborus[0].y][snakes[s].oroborus[0].z]==1) //electric snake?
           g=1;
         else
           r=1;

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

           fprintf(fo,"cylinder{< %f, %f, %f>, <%f, %f, %f>,%f texture{ pigment { rgb <%f,%f,%f>}}}\n",
                  (float)snakes[s].oroborus[seg%snakes[s].segs].x/Z,
                  (float)snakes[s].oroborus[seg%snakes[s].segs].y/Z,
                  (float)snakes[s].oroborus[seg%snakes[s].segs].z/Z,
		  (float)snakes[s].oroborus[(seg+1)%snakes[s].segs].x/Z,
		  (float)snakes[s].oroborus[(seg+1)%snakes[s].segs].y/Z,
		  (float)snakes[s].oroborus[(seg+1)%snakes[s].segs].z/Z,
		  0.1/Z,
		  r,g,b);
	     
	  fprintf(fo,"sphere{< %f, %f, %f>,%f texture{ pigment { rgb <%f,%f,%f>}}}\n",
		 (float)snakes[s].oroborus[seg%snakes[s].segs].x/Z,
		 (float)snakes[s].oroborus[seg%snakes[s].segs].y/Z,
		 (float)snakes[s].oroborus[seg%snakes[s].segs].z/Z,
		 0.15/Z,
		 r,g,b);
		 
 //      printf("\n");
	  } 
	
      }
   fclose(fo);
 } 
*/

void
print_lattice_pnm ()
{
  printf ("P2\n%d %d\n%d\n", X, Y, Z);

  int x, y, z = 0,count=0;
  for (x = 0; x < X; x++)
    {
      for (y = 0; y < Y; y++)
	{
	   for (z=0;z<Z;z++)
	     if (lattice[x][z][z]!=-1)
	       count++;
	   
          printf ("%d\t", count);
	}

      printf ("\n");
    }
  printf ("\n\n");
}

void save_lattice_file(char * name)
{
   int x,y,z,count;
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
   

     fprintf(fo,"\n\n%d %f %f\n",TOTAL_SLITHERS,DENSITY,iE[1][1]);
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
     //exit(-1);
     }
   
     for (x=0;x<X;x++)
       for (y=0;y<Y;y++)
	 for (z=0;z<Z;z++)
	   fscanf(fo,"%d",&lattice[x][y][z]);
   fclose(fo);
}   

void print_lattice_pnm_file(int i)
{
 FILE *fo;
   char name[30]; //static buffers, security flaws, yada yada
     int x, y, z = 0,count=0;
   
   sprintf(name,"marbles_%d.pgm",i);
   fo=fopen(name,"w");   
   fprintf (fo,"P2\n%d %d\n%d\n", X, Y, Z/4);
   

     for (y = 0; y < Y; y++)
     {
	      for (x = 0; x < X; x++)
	  {
	   for (z=0;z<Z;z++)
	     if (lattice[x][y][z]!=-1)
	       count++;
	   if (count>Z/4)
	       {
		  fprintf(stderr,"High Stack; loss of contrast expected. Rescaling %d to %d.\n",count,Z/4);
		  count=Z/4;
		  
	       }
	     
          fprintf (fo,"%d\t", count);
	     
	     count=0;
	  }
	      fprintf (fo,"\n");
     }
   
     fprintf (fo,"\n\n");
   fclose(fo);
   
}


/*void
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
*/
