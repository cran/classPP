#include <R.h> 
#include <stdio.h> 
#include <string.h> 
#include <stdlib.h> 
#include <math.h> 
#include <time.h>


void cartgini(int *fn, int *fp, int *groups,double *fvals,int *group, int *ngroup,double *val); 
void cartentropy(int *fn, int *fp, int *groups,double *fvals,int *group, int *ngroup,double *val); 
void zero_int(int *mem, int size); 
void sort_group(double *vals,int *group,int *fp,int *fn,int left, int right) ;
void zero(double *mem, int size); 
void sort_data(double *x, int *index,int left, int right); 


void cartgini(int *fn, int *fp, int *groups,double *fvals,int *group, int *ngroup,double *val)
{
 int i, k,n,p,g,left,right,l,*xindex,j; 
 double dev, prob,index,minindex, *vals; 
 double *x;
 int *nright;

  n=*fn; p=*fp; g=*groups;

  vals = Calloc(p*n,double); 
  for(i=0; i<p; i++)
  { for(j=0; j<n; j++)   
        vals[i*n+j] = fvals[i*n+j];
  }     
  


  right = n-1; 
  left = 0;
  xindex = Calloc(n,int);zero_int(xindex,n);
  for (i=0; i<n; i++) 
    xindex[i] = group[i]; 
  sort_group(vals,xindex,fp,fn,left,right); 

/* data relocation and make index */ 
  x = Calloc(n,double); zero(x,n); 
  minindex=0; 
  nright = Calloc(g,int);
for(l=0; l<p; l++)
{ 
  for (i=0; i<n; i++) 
   x[i] = vals[l*n+i]; 


  left=0; 
  right=n-1; 
  sort_data(x, xindex,left,right) ; 

/* Calculate gini index */ 

  zero_int(nright,g); 
  index=1;
  for (i=0; i<g; i++) 
  {
     prob = ((double)ngroup[i])/((double)n); 
     if (prob>0) index -= prob*prob; 
   }
  for (i=0; i<n-1; i++) 
  {  (nright[(xindex[i]-1)])++; 
     dev=1;
     for (k=0; k<g; k++) 
     { prob = ((double) nright[k])/((double)(i+1)); 
       dev -= ((double)(i+1))*(prob*prob)/((double)n);  
       prob = ((double) (ngroup[k]-nright[k]))/((double)(n-i-1)); 
       dev -= ((double)(n-i-1))*(prob*prob)/((double)n);  

     } 
   if (dev<index) index = dev; 
  } 
   if(l==0 ) minindex= index; 
  else if(minindex < index) minindex = index;
}
*val = 1-minindex;

 Free(vals);
 Free(x);
 Free(xindex);
 Free(nright);
} 
/***********************
Entropy

*********************************************************************/ 

void cartentropy(int *fn, int *fp, int *groups,double *fvals,int *group, int *ngroup,double *val)
{
int i, k,n,p,g,left,right,l,*xindex,j; 
double dev, prob,index,minindex, *vals; 
 double *x;
 int *nright;

  n=*fn; p=*fp; g=*groups;

  vals = Calloc(p*n,double); 
  for(i=0; i<p; i++)
  { for(j=0; j<n; j++)   
        vals[i*n+j] = fvals[i*n+j];
  }     
  

  right = n-1; 
  left = 0;
  xindex = Calloc(n,int);zero_int(xindex,n); 
  for (i=0; i<n; i++) 
    xindex[i] = group[i]; 
  sort_group(vals,xindex,fp,fn,left,right); 

/* data relocation and make index */ 
  x = Calloc(n,double); zero(x,n); 

  minindex=0; 
  nright = Calloc(g,int);
for(l=0; l<p; l++)
{ 
  for (i=0; i<n; i++) 
    x[i] = vals[l*n+i]; 
   
  left=0; 
  right=n-1; 
  sort_data(x, xindex,left,right) ; 

  zero_int(nright,g); 
 index = 0; 
  for (i=0; i<g; i++) 
  {
     prob = ((double)ngroup[i])/((double)n); 
     if (prob>0) index -= prob*log(prob); 
   }
   for (i=0; i<n-1; i++) 
   {   
      (nright[(xindex[i]-1)])++; 
      dev=0; 
      for (k=0; k<g; k++) 
	{  prob = ((double) nright[k])/((double)(i+1));  
         if(prob>0) dev -= ((double)(i+1))*prob*log(prob)/((double)(n)); 
         prob = ((double) (ngroup[k]-nright[k]))/((double)(n-i-1)); 
         if(prob>0) dev -=((double)(n-i-1))* prob*log(prob)/((double)(n)); 
      } 
      if (dev<index)
      index = dev; 

   }
   if(l==0 ) minindex= index; 
  else if(minindex < index) minindex = index;

}
*val = 1-minindex/log(g);

 Free(x);
 Free(xindex);
 Free(nright);

} 
