
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi/mpi.h>

#define MAXFN   200
#define MAXLINE 2000
#define MAXLOCS 200000

#define  ARRAYSIZE	400000
#define  MASTER		0
#define MTAG1 2
#define MTAG2 1

float  data[ARRAYSIZE];
MPI_Status status;
int numtasks, taskid,block, rc, offset, source;


int
load_locations (char *fn, float locs[][2], int first)
{
  char fullfn[MAXFN], line[MAXLINE], *kwd;
  int nlocs;
  FILE *f;

  sprintf (fullfn, "./data/%s", fn);
  if ((f = fopen (fullfn, "r")) == NULL) {
    fprintf (stderr, "Cannot open %s!\n", fullfn);
    exit (EXIT_FAILURE);
  }
  fgets (line, MAXLINE, f);
  fgets (line, MAXLINE, f);
  nlocs = first;
  while (fgets (line, MAXLINE, f) != NULL) {
    kwd = strtok (line, " ");
    locs[nlocs][1] = atof (kwd);
    kwd = strtok ((char *) NULL, " ");
    locs[nlocs][0] = atof (kwd);
    nlocs += 1;
    if (nlocs >= MAXLOCS) {
      fprintf (stderr, "Not enough space to hold all locations: %d!\n", nlocs);
      exit (EXIT_FAILURE);
    }
  }
  fclose (f);
  return nlocs;
}

//This is the function that we want to "distribute" to master and slaves. 
float
func(float myoffset, int fraction, int myid,float locs[][2],int i,float y1,float x1, int *nc,float *dsum,int *np,int *np_overall){
	int j, y2,x2,d;
	float hmsum = 0.0, hm;
	for (j = myoffset ; j < myoffset+fraction ; j++) {
	      if (i != j) {
		y2 = locs[j][0] - y1;
		x2 = locs[j][1] - x1;
		d = sqrt (y2 * y2 + x2 * x2);
		if (d <= 0.0) {
		  *nc += 1;
		} else {
		  *dsum += 1.0 / d;
		  *np += 1;
		  *np_overall += 1;
		}
	      }
	 }
	hm = *np/ *dsum;
	hmsum += 1.0 / hm;
	return hmsum;
}


void
coverage (float locs[][2], int nlocs, float *result, int *npaths, int *ncoin)
{

  
  int np_overall = 0, nc = 0, np, i, j;
  float hmsum = 0.0, dsum, hm, d, y1, x1, y2, x2,path,myhmsum;
  //For master node
  if(taskid == MASTER){
	//Here we portion the data in the arrays
	for(i=0; i<ARRAYSIZE; i++) {
	    data[i] =  i * 1.0;
	}	
	/*  

	size=numprocs;
	chunk=nlocs/size;
		
	for( i=0;i<size;i++){
		for(j=0;j<.....){
			//I should not have done it using arrays but the right way is this and in here i would assign the data to each slave beginning from the master like this
			master----> data(0) to chunk-1
			slave1---->data(chunk) to 2*chunk-
			slave2---->data(2*chunk) to 3*chunk-1
			slave3---->data(3*chunk) to 4*chunk-1
		After this master will send with MPI_Send the data to the slaves as the procedure continues down from here..
			}
			
		}
		
  */
		
		
	
	}
	offset = block;
	for (path=1; path<numtasks; path++) {
	    MPI_Send(&offset, 1, MPI_INT, path, MTAG1, MPI_COMM_WORLD);
	    MPI_Send(&data[offset], block, MPI_FLOAT, path, MTAG2, MPI_COMM_WORLD);
	    offset = offset + block;
	}
	offset = 0;
	for (i = 0; i < nlocs; i++) {
	    y1 = locs[i][0];
	    x1 = locs[i][1];
	    np = 0;
	    dsum = 0.0;
	    myhmsum = func(offset,block,taskid,locs,i,y1,x1,&nc,&dsum,&np,&np_overall);
	}

	// Receiving results from every single task
	for (i=1; i<numtasks; i++) {
	    source = i;
	    MPI_Recv(&offset, 1, MPI_INT, source, MTAG1, MPI_COMM_WORLD, &status);
	    MPI_Recv(&data[offset], block, MPI_FLOAT, source, MTAG2,
	      MPI_COMM_WORLD, &status);
	}

	// Obtaining the total sum  
  	MPI_Reduce(&myhmsum, &hmsum, 1, MPI_FLOAT, MPI_SUM, MASTER, MPI_COMM_WORLD);
	
	*result = nlocs / hmsum;
	*npaths = np_overall;
	*ncoin = nc;
	
	//Node of slave. 
  }if (taskid != MASTER){
	// Receiving chunk of array from the master task
  	source = MASTER;
  	MPI_Recv(&offset, 1, MPI_INT, source, MTAG1, MPI_COMM_WORLD, &status);
  	MPI_Recv(&data[offset], block, MPI_FLOAT, source, MTAG2, MPI_COMM_WORLD, &status);
	for (i = 0; i < nlocs; i++) {
	    y1 = locs[i][0];
	    x1 = locs[i][1];
	    np = 0;
	    dsum = 0.0;
	    myhmsum = func(offset,block,taskid,locs,i,y1,x1,&nc,&dsum,&np,&np_overall);
	   
	  }
	// Deliver the results to the master task back again
	  path = MASTER;
	  MPI_Send(&offset, 1, MPI_INT, path, MTAG1, MPI_COMM_WORLD);
	  MPI_Send(&data[offset], block, MPI_FLOAT, MASTER, MTAG2, MPI_COMM_WORLD);

	  MPI_Reduce(&myhmsum, &hmsum, 1, MPI_FLOAT, MPI_SUM, MASTER, MPI_COMM_WORLD);

	}
	 


}



float
mean_coverage1 (char *opname, char **imsets, int nimsets)
{
  int nfiles = 0, i, fno, nlocs, np, nc;
  float locs[MAXLOCS][2];
  float csum = 0.0, cov;
  char fn[MAXFN];

  for (i = 0; i < nimsets; i++) {
    for (fno = 1; fno < 7; fno++) {
      sprintf (fn, "%s_%s%d.txt", opname, imsets[i], fno);
      nlocs = load_locations (fn, locs, 0);
      coverage (locs, nlocs, &cov, &np, &nc);
      // printf ("  %s %s %f %d %d\n", opname, fn, cov, np, nc);
      csum += cov;
      nfiles += 1;
    }
  }
  return csum / nfiles;
}


float
mean_coverage2 (char *opname1, char *opname2, char **imsets, int nimsets)
{
  int nfiles = 0, i, fno, nlocs, np, nc;
  float locs[MAXLOCS][2];
  float csum = 0.0, cov;
  char fn1[MAXFN], fn2[MAXFN];

  for (i = 0; i < nimsets; i++) {
    for (fno = 1; fno < 7; fno++) {
      sprintf (fn1, "%s_%s%d.txt", opname1, imsets[i], fno);
      nlocs = load_locations (fn1, locs, 0);
      sprintf (fn2, "%s_%s%d.txt", opname2, imsets[i], fno);
      nlocs = load_locations (fn2, locs, nlocs);
      coverage (locs, nlocs, &cov, &np, &nc);
      // printf ("  %s %s %f %d %d\n", fn1, fn2, cov, np, nc);
      csum += cov;
      nfiles += 1;
    }
  }
  return csum / nfiles;
}


float
mean_coverage3 (char *opname1, char *opname2, char *opname3,
		char **imsets, int nimsets)
{
  int nfiles = 0, i, fno, nlocs, np, nc;
  float locs[MAXLOCS][2];
  float csum = 0.0, cov;
  char fn1[MAXFN], fn2[MAXFN], fn3[MAXFN];

  for (i = 0; i < nimsets; i++) {
    for (fno = 1; fno < 7; fno++) {
      sprintf (fn1, "%s_%s%d.txt", opname1, imsets[i], fno);
      nlocs = load_locations (fn1, locs, 0);
      sprintf (fn2, "%s_%s%d.txt", opname2, imsets[i], fno);
      nlocs = load_locations (fn2, locs, nlocs);
      sprintf (fn3, "%s_%s%d.txt", opname3, imsets[i], fno);
      nlocs = load_locations (fn3, locs, nlocs);
      coverage (locs, nlocs, &cov, &np, &nc);
      // printf ("  %s %s %s %f %d %d\n", fn1, fn2, fn3, cov, np, nc);
      csum += cov;
      nfiles += 1;
    }
  }
  return csum / nfiles;
}


float
mean_coverage4 (char *opname1, char *opname2, char *opname3, char *opname4,
		char **imsets, int nimsets)
{
  int nfiles = 0, i, fno, nlocs, np, nc;
  float locs[MAXLOCS][2];
  float csum = 0.0, cov;
  char fn1[MAXFN], fn2[MAXFN], fn3[MAXFN], fn4[MAXFN];

  for (i = 0; i < nimsets; i++) {
    for (fno = 1; fno < 7; fno++) {
      sprintf (fn1, "%s_%s%d.txt", opname1, imsets[i], fno);
      nlocs = load_locations (fn1, locs, 0);
      sprintf (fn2, "%s_%s%d.txt", opname2, imsets[i], fno);
      nlocs = load_locations (fn2, locs, nlocs);
      sprintf (fn3, "%s_%s%d.txt", opname3, imsets[i], fno);
      nlocs = load_locations (fn3, locs, nlocs);
      sprintf (fn4, "%s_%s%d.txt", opname4, imsets[i], fno);
      nlocs = load_locations (fn4, locs, nlocs);
      coverage (locs, nlocs, &cov, &np, &nc);
      // printf ("  %s %s %s %s %f %d %d\n", fn1, fn2, fn3, fn4, cov, np, nc);
      csum += cov;
      nfiles += 1;
    }
  }
  return csum / nfiles;
}


int
main (int argc, char **argv)
{
  char *opnames[] = {"ebr", "ibr", "mser", "sfop"};
  int nopnames = 4;
  char *imsets[] = {"bark", "bikes", "boat", "graf", "leuv", "trees", "ubc",
		    "wall"};
  int nimsets = 8;
  int i1, i2, i3, i4;
  float mc;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
  if (numtasks % 4 != 0) {
     MPI_Abort(MPI_COMM_WORLD, rc);
     exit(0);
  }
  MPI_Comm_rank(MPI_COMM_WORLD,&taskid);
  block = (ARRAYSIZE / numtasks);
  
 
  // Calculate the coverage for each individual operator.
  for (i1 = 0; i1 < nopnames; i1++) {
    mc = mean_coverage1 (opnames[i1], imsets, nimsets);
    printf ("%s: %f\n", opnames[i1], mc);
  }

  // Calculate the coverage for each combination of pairs of operators.
  for (i1 = 0; i1 < nopnames; i1++) {
    for (i2 = i1+1; i2 < nopnames; i2++) {
      mc = mean_coverage2 (opnames[i1], opnames[i2], imsets, nimsets);
      printf ("%s + %s: %f\n", opnames[i1], opnames[i2], mc);
    }
  }

  // Calculate the coverage for each combination of triplets of operators.
  for (i1 = 0; i1 < nopnames; i1++) {
    for (i2 = i1+1; i2 < nopnames; i2++) {
      for (i3 = i2+1; i3 < nopnames; i3++) {
	mc = mean_coverage3 (opnames[i1], opnames[i2], opnames[i3],
			     imsets, nimsets);
	printf ("%s + %s + %s: %f\n", opnames[i1], opnames[i2], opnames[i3], mc);
      }
    }
  }

  // Calculate the coverage for each combination of quadruplets of operators.
  for (i1 = 0; i1 < nopnames; i1++) {
    for (i2 = i1+1; i2 < nopnames; i2++) {
      for (i3 = i2+1; i3 < nopnames; i3++) {
	for (i4 = i3+1; i4 < nopnames; i4++) {
	  mc = mean_coverage4 (opnames[i1], opnames[i2], opnames[i3],
			       opnames[i4],  imsets, nimsets);
	  printf ("%s + %s + %s + %s: %f\n", opnames[i1], opnames[i2],
		  opnames[i3], opnames[i4], mc);
	}
      }
    }
  }
  MPI_Finalize();
  return EXIT_SUCCESS;
}

//-----------------------------------------------------------------------------
// End of coverage.c
//-----------------------------------------------------------------------------
