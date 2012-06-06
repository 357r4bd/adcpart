#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// may want to use this at some point, so keep it
typedef struct {
  //int degree, *adjNodes;
  double x,y,z;
} GraphNode;

GraphNode *nodes;
int numNodes, numElements;

void SkipToEndofLine(FILE *fpInp)
{ while (getc(fpInp) != '\n');}

void ReadInput() //read in the grid
{ int i, j, nodeNum, tmp;
  FILE *fp;
  fp = stdin;                                // modified code to take in graph def via stdin
  SkipToEndofLine(fp);                       // skip grid name
  fscanf(fp,"%d %d",&numElements,&numNodes); // read in numElements and numNodes
  // allocate memory to store list of vertices
  nodes = (GraphNode *)malloc(numNodes*sizeof(GraphNode));
  // read in vertices - note the extreme vertices for defining of box
  for (i=0; i<numNodes; i++)
  { fscanf(fp,"%d %lf %lf %lf",&tmp,&nodes[i].x,&nodes[i].y,&nodes[i].z);
     //printf("%lf %lf %lf\n",nodes[i].x,nodes[i].y,nodes[i].z);
    SkipToEndofLine(fp);
  }
  // close after reading in vertices - we don't care about anything else
  fclose(fp);
}

int main(int argc, char **argv)
{ int i; 
  ReadInput(); // suck in grid
  return 1;
}

