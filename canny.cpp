/* source: http://marathon.csee.usf.edu/edge/edge_detection.html */
/* URL: ftp://figment.csee.usf.edu/pub/Edge_Comparison/source_code/canny.src */

/* ECPS 203 Assignment 1 solution, start for Assignment 2 */

/* off-by-one bugs fixed by Rainer Doemer, 10/05/23 */

/* Last modified by Abyukth Kumar on 10/17/2023	*/

/* Quick Note : This code has been indented with tabstop=3" */
/*******************************************************************************
* --------------------------------------------
*(c) 2001 University of South Florida, Tampa
* Use, or copying without permission prohibited.
* PERMISSION TO USE
* In transmitting this software, permission to use for research and
* educational purposes is hereby granted.  This software may be copied for
* archival and backup purposes only.  This software may not be transmitted
* to a third party without prior permission of the copyright holder. This
* permission may be granted only by Mike Heath or Prof. Sudeep Sarkar of
* University of South Florida (sarkar@csee.usf.edu). Acknowledgment as
* appropriate is respectfully requested.
* 
*  Heath, M., Sarkar, S., Sanocki, T., and Bowyer, K. Comparison of edge
*    detectors: a methodology and initial study, Computer Vision and Image
*    Understanding 69 (1), 38-54, January 1998.
*  Heath, M., Sarkar, S., Sanocki, T. and Bowyer, K.W. A Robust Visual
*    Method for Assessing the Relative Performance of Edge Detection
*    Algorithms, IEEE Transactions on Pattern Analysis and Machine
*    Intelligence 19 (12),  1338-1359, December 1997.
*  ------------------------------------------------------
*
* PROGRAM: canny_edge
* PURPOSE: This program implements a "Canny" edge detector. The processing
* steps are as follows:
*
*   1) Convolve the image with a separable gaussian filter.
*   2) Take the dx and dy the first derivatives using [-1,0,1] and [1,0,-1]'.
*   3) Compute the magnitude: sqrt(dx*dx+dy*dy).
*   4) Perform non-maximal suppression.
*   5) Perform hysteresis.
*
* We have three parameters required to perform this detection algorithm. 
* These are as follows:
*
*   SIGMA = The standard deviation of the gaussian smoothing filter.
*   T_LOW  = Specifies the low value to use in hysteresis. This is a 
*           fraction (0-1) of the computed high threshold edge strength value.
*   T_HIGH = Specifies the high value to use in hysteresis. This fraction (0-1)
*           specifies the percentage point in a histogram of the gradient of
*           the magnitude. Magnitude values of zero are not counted in the
*           histogram.
*
*	 The hard coded values of these parameters are as follows:
*	 SIGMA 	- 	0.6
*	 T_LOW 	- 	0.3
*	 T_HIGH 	-	0.8 
*
*******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define VERBOSE 	1		//For debugging

#define SUCCESS	0
#define FAILED		1

#define ROWS				240
#define COLUMNS			320
#define IMAGE_SIZE		ROWS*COLUMNS
#define MAX_GRAY_VAL		255

#define SIGMA 				0.6	
#define T_LOW				0.3	
#define T_HIGH				0.8

#define WINDOWSIZE	  	(int)(1 + 2*ceil(2.5 * SIGMA))

#define NOEDGE 			255
#define POSSIBLE_EDGE 	128
#define EDGE 				0
#define BOOSTBLURFACTOR 90.0
#ifndef M_PI
#define M_PI 3.14159265356789202346
#endif

//Functions used to read and write the image data in PGM format 
int read_pgm_image(char *infilename, unsigned char *image);
int write_pgm_image(char *outfilename, unsigned char *image);

//Functions used to perform canny edge detection
void canny(unsigned char *image, unsigned char *edge);
void gaussian_smooth(unsigned char *image, short int *smoothedim);
void make_gaussian_kernel(float *kernel);
void derivative_x_y(short int *smoothedim, short int *delta_x, short int *delta_y);
void magnitude_x_y(short int *delta_x, short int *delta_y, short int *magnitude);
void apply_hysteresis(short int *mag, unsigned char *nms, unsigned char *edge);
void non_max_supp(short *mag, short *gradx, short *grady, unsigned char *result);
void follow_edges(unsigned char *edgemapptr, short *edgemagptr, short lowval);

int main()
{
   char infilename[] = "golfcart.pgm";  						 /* Name of the i/p image */
   char outfilename[] = "golfcart_s_0.6_l_0.3_h_0.8.pgm"; /* Name of the o/p image */
   unsigned char image[IMAGE_SIZE];     						 /* The i/p image */
   unsigned char edge[IMAGE_SIZE];      						 /* The o/p edge image */

   /****************************************************************************
   * Read in the image. This read function allocates memory for the image.
   ****************************************************************************/
   
	if(VERBOSE) 
		printf("Reading the image %s.\n", infilename);

   if(read_pgm_image(infilename, image) == FAILED)
	{
      fprintf(stderr, "Error reading the input image, %s.\n", infilename);
      exit(EXIT_FAILURE);
   }

   /****************************************************************************
   * Perform the edge detection. All of the work takes place here.
   ****************************************************************************/
   
	if(VERBOSE) 
		printf("Starting Canny edge detection.\n");
   
	canny(image, edge);

   /****************************************************************************
   * Write out the edge image to a file.
   ****************************************************************************/

   if(VERBOSE) 
		printf("Writing the edge iname in the file %s.\n", outfilename);
   
	if(write_pgm_image(outfilename, edge) == FAILED)
	{
      fprintf(stderr, "Error writing the edge image, %s.\n", outfilename);
      exit(1);
   }

   return(0); /* exit cleanly */
}

/*******************************************************************************
* FUNCTION: 	canny
* PURPOSE: 		To perform canny edge detection.
* ARGUMENTS: 	unsigned char *image -	Buffer containing image data
* 					unsigned char *edge	-	Buffer to store o/p edge image data
* RETURN: 		None
*******************************************************************************/
void canny(unsigned char *image, unsigned char *edge)
{
   unsigned char nms[IMAGE_SIZE];    /* Points that are local maximal magnitude. */
   short int smoothedim[IMAGE_SIZE]; /* The image after gaussian smoothing.      */
   short int delta_x[IMAGE_SIZE];    /* The first derivative image, x-direction. */
   short int delta_y[IMAGE_SIZE];    /* The first derivative image, y-direction. */
   short int magnitude[IMAGE_SIZE];	 /* The magnitude of the gradient image.     */

   /****************************************************************************
   * Perform gaussian smoothing on the image using the input standard
   * deviation.
   ****************************************************************************/
   
	if(VERBOSE) 
		printf("Smoothing the image using a gaussian kernel.\n");
   
	gaussian_smooth(image, smoothedim);

   /****************************************************************************
   * Compute the first derivative in the x and y directions.
   ****************************************************************************/
   
	if(VERBOSE)
		printf("Computing the X and Y first derivatives.\n");
   
	derivative_x_y(smoothedim, delta_x, delta_y);

   /****************************************************************************
   * Compute the magnitude of the gradient.
   ****************************************************************************/
   
	if(VERBOSE) 
		printf("Computing the magnitude of the gradient.\n");
   
	magnitude_x_y(delta_x, delta_y, magnitude);

   /****************************************************************************
   * Perform non-maximal suppression.
   ****************************************************************************/
   
	if(VERBOSE) 
		printf("Doing the non-maximal suppression.\n");
   
	non_max_supp(magnitude, delta_x, delta_y, nms);

   /****************************************************************************
   * Use hysteresis to mark the edge pixels.
   ****************************************************************************/
   
	if(VERBOSE)
		printf("Doing hysteresis thresholding.\n");
   
	apply_hysteresis(magnitude, nms, edge);
}

/*******************************************************************************
* FUNCTION: 	magnitude_x_y
* PURPOSE: 		Compute the magnitude of the gradient. This is the square root of
* 					the sum of the squared derivative values.
* ARGUMENTS:	short int *delta_x	-	Derivative in x direction
* 					short int *delta_y	-	Derivative in y direction
* 					short int *magnitude	-	Variable to store magnitude of gradient 
* RETURN:		NONE
*******************************************************************************/
void magnitude_x_y(short int *delta_x, short int *delta_y, short int *magnitude)
{
   int r, c, pos, sq1, sq2;

   /****************************************************************************
   * Computing magnitude of the gradient.
   ****************************************************************************/
   
	for(r=0, pos=0; r < ROWS; r++)
	{
      for(c=0; c < COLUMNS; c++, pos++)
		{
         sq1 = (int)delta_x[pos] * (int)delta_x[pos];
         sq2 = (int)delta_y[pos] * (int)delta_y[pos];
         magnitude[pos] = (short)(0.5 + sqrt((float)sq1 + (float)sq2));
      }
   }
}

/*******************************************************************************
* FUNCTION: 	derivative_x_y
* PURPOSE: 		Compute the first derivative of the image in both the x any y
* 					directions. The differential filters that are used are:
*
*              		                            -1
*         		dx =  -1 0 +1     and       dy =   0
*                 		                         +1
*
* ARGUMENTS:	short int *smoothedim	-	Smoothened image
* 					short int *delta_x		-	Derivative in x direction
* 					short int *delta_y		-	Derivative in y direction				
* RETURN:		None
*******************************************************************************/
void derivative_x_y(short int *smoothedim, short int *delta_x,
		    short int *delta_y)
{
   int r, c, pos;

   /****************************************************************************
   * Compute the x-derivative. Adjust the derivative at the borders to avoid
   * losing pixels.
   ****************************************************************************/
   
	if(VERBOSE) 
		printf("   Computing the X-direction derivative.\n");
   
	for(r=0; r < ROWS; r++)
	{
      pos = r * COLUMNS;
      delta_x[pos] = smoothedim[pos+1] - smoothedim[pos];
      pos++;
      
		for(c=1; c < (COLUMNS-1); c++, pos++)
		{
         delta_x[pos] = smoothedim[pos+1] - smoothedim[pos-1];
      }
      delta_x[pos] = smoothedim[pos] - smoothedim[pos-1];
   }

   /****************************************************************************
   * Compute the y-derivative. Adjust the derivative at the borders to avoid
   * losing pixels.
   ****************************************************************************/
   
	if(VERBOSE)
		printf("   Computing the Y-direction derivative.\n");
   
	for(c=0; c < COLUMNS; c++)
	{
      pos = c;
      delta_y[pos] = smoothedim[pos+COLUMNS] - smoothedim[pos];
      pos += COLUMNS;
      
		for(r=1;	r < (ROWS-1); r++, pos+=COLUMNS)
		{
         delta_y[pos] = smoothedim[pos+COLUMNS] - smoothedim[pos-COLUMNS];
      }
      delta_y[pos] = smoothedim[pos] - smoothedim[pos-COLUMNS];
   }
}

/*******************************************************************************
* FUNCTION: 	gaussian_smooth
* PURPOSE: 		Blur an image with a gaussian filter.
* ARGUMENTS:	unsigned char *image		-	Buffer containing image data
* 					short int *smoothedim	-	Buffer to store smoothened image	
* RETURN:		None
*******************************************************************************/
void gaussian_smooth(unsigned char *image, short int *smoothedim)
{
   int r, c, rr, cc;     			/* Counter variables */
   int center = WINDOWSIZE/2; 	/* Half of the windowsize */
   float tempim[IMAGE_SIZE]; 		/* Buffer for separable filter gaussian smoothing */
   float kernel[WINDOWSIZE];  	/* A one dimensional gaussian kernel */
   float dot;            			/* Dot product summing variable */
   float sum;            			/* Sum of the kernel weights variable */

   /****************************************************************************
   * Create a 1-dimensional gaussian smoothing kernel.
   ****************************************************************************/
   if(VERBOSE) 
		printf("   Computing the gaussian smoothing kernel.\n");
   
	make_gaussian_kernel(kernel);

   /****************************************************************************
   * Blur in the x - direction.
   ****************************************************************************/
   
	if(VERBOSE)
		printf("   Blurring the image in the X-direction.\n");
  
	for(r=0; r < ROWS; r++)
	{
      for(c=0; c < COLUMNS; c++)
		{
         dot = 0.0;
         sum = 0.0;
         
			for(cc=(-center); cc <= center; cc++)
			{
            if(((c+cc) >= 0) && ((c+cc) < COLUMNS))
				{
               dot += (float)image[r*COLUMNS+(c+cc)] * kernel[center+cc];
               sum += kernel[center+cc];
            }
         }
         tempim[r*COLUMNS+c] = dot/sum;
      }
   }

   /****************************************************************************
   * Blur in the y - direction.
   ****************************************************************************/
   
	if(VERBOSE) 
		printf("   Blurring the image in the Y-direction.\n");
   
	for(c=0; c < COLUMNS; c++)
	{
      for(r=0; r < ROWS; r++)
		{
         sum = 0.0;
         dot = 0.0;
         
			for(rr=(-center); rr <= center; rr++)
			{
            if(((r+rr) >= 0) && ((r+rr) < ROWS))
				{
               dot += tempim[(r+rr)*COLUMNS+c] * kernel[center+rr];
               sum += kernel[center+rr];
            }
         }
         smoothedim[r*COLUMNS+c] = (short int)(dot*BOOSTBLURFACTOR/sum + 0.5);
      }
   }
}

/*******************************************************************************
* FUNCTION:		make_gaussian_kernel
* PURPOSE: 		Create a one dimensional gaussian kernel.
* ARGUMENTS:	float *kernel	-	Buffer to store kernel coefficients
* RETURN: 		None
*******************************************************************************/
void make_gaussian_kernel(float *kernel)
{
   int i; 
	int center = WINDOWSIZE/2;
   float x, fx, sum=0.0;

   if(VERBOSE) 
		printf("      The kernel has %d elements.\n", WINDOWSIZE);
   
	for(i=0; i < WINDOWSIZE; i++)
	{
      x = (float)(i - center);
      fx = pow(2.71828, -0.5*x*x/(SIGMA*SIGMA)) / (SIGMA * sqrt(6.2831853));
      kernel[i] = fx;
      sum += fx;
   }

   for(i=0; i < WINDOWSIZE; i++) 
		kernel[i] /= sum;

   if(VERBOSE)
	{
      printf("The filter coefficients are:\n");
      for(i=0; i < WINDOWSIZE; i++)
         printf("kernel[%d] = %f\n", i, kernel[i]);
   }
}


/*******************************************************************************
* FUNCTION: 	follow_edges
* PURPOSE: 		This procedure edges is a recursive routine that traces edges 
* 					along all paths whose magnitude values remain above a specified 
* 					lower threshhold.
* ARGUMENTS:	unsigned char *edgemapptr	-	Edge map pointer		
* 					short *edgemagptr				-	Edge magnitude pointer
* 					short lowval					-	Lower threshold value
* RETURN:		None
*******************************************************************************/
void follow_edges(unsigned char *edgemapptr, short *edgemagptr, short lowval)
{
   short *tempmagptr;
   unsigned char *tempmapptr;
   int i;
   int x[8] = {1,1,0,-1,-1,-1,0,1},
       y[8] = {0,1,1,1,0,-1,-1,-1};

   for(i=0; i < 8; i++)
	{
      tempmapptr = edgemapptr - y[i]*COLUMNS + x[i];
      tempmagptr = edgemagptr - y[i]*COLUMNS + x[i];

      if((*tempmapptr == POSSIBLE_EDGE) && (*tempmagptr > lowval))
		{
         *tempmapptr = (unsigned char) EDGE;
         follow_edges(tempmapptr,tempmagptr, lowval);
      }
   }
}

/*******************************************************************************
* FUNCTION: 	apply_hysteresis
* PURPOSE: 		This routine finds edges that are above some high threshhold or
* 					are connected to a high pixel by a path of pixels greater than a 
* 					low threshold.
* ARGUMENTS:	short int *mag				
* 					unsigned char *nms	
*					unsigned char *edge
* RETURN:		None
*******************************************************************************/
void apply_hysteresis(short int *mag, unsigned char *nms, unsigned char *edge)
{
   int r, c, pos, numedges, highcount, lowthreshold, highthreshold, hist[32768];
	short int maximum_mag = 0;

   /****************************************************************************
   * Initialize the edge map to possible edges everywhere the non-maximal
   * suppression suggested there could be an edge except for the border. At
   * the border we say there can not be an edge because it makes the
   * follow_edges algorithm more efficient to not worry about tracking an
   * edge off the side of the image.
   ****************************************************************************/
   
	for(r=0, pos=0; r < ROWS; r++)
	{
      for(c=0; c < COLUMNS; c++, pos++)
		{
	 		if(nms[pos] == POSSIBLE_EDGE) 
				edge[pos] = POSSIBLE_EDGE;
	 		else 
				edge[pos] = NOEDGE;
      }
   }

   for(r=0, pos=0; r < ROWS; r++, pos+=COLUMNS)
	{
      edge[pos] = NOEDGE;
      edge[pos+COLUMNS-1] = NOEDGE;
   }
   pos = (ROWS-1) * COLUMNS;
   
	for(c=0; c < COLUMNS; c++, pos++)
	{
      edge[c] = NOEDGE;
      edge[pos] = NOEDGE;
   }

   /****************************************************************************
   * Compute the histogram of the magnitude image. Then use the histogram to
   * compute hysteresis thresholds.
   ****************************************************************************/
   
	for(r=0; r < 32768; r++) 
		hist[r] = 0;

   for(r=0, pos=0; r < ROWS; r++)
	{
      for(c=0; c < COLUMNS; c++, pos++)
		{
	 		if(edge[pos] == POSSIBLE_EDGE) 
				hist[mag[pos]]++;
      }
   }

   /****************************************************************************
   * Compute the number of pixels that passed the nonmaximal suppression.
   ****************************************************************************/
   
	for(r=1, numedges=0; r < 32768; r++)
	{
      if(hist[r] != 0)
			 maximum_mag = r;
      
		numedges += hist[r];
   }

   highcount = (int)(numedges * T_HIGH + 0.5);

   /****************************************************************************
   * Compute the high threshold value as the (100 * T_HIGH) percentage point
   * in the magnitude of the gradient histogram of all the pixels that passes
   * non-maximal suppression. Then calculate the low threshold as a fraction
   * of the computed high threshold value. John Canny said in his paper
   * "A Computational Approach to Edge Detection" that "The ratio of the
   * high to low threshold in the implementation is in the range two or three
   * to one." That means that in terms of this implementation, we should
   * choose T_LOW ~= 0.5 or 0.33333.
   ****************************************************************************/
   
	r = 1;
   numedges = hist[1];
   
	while((r < (maximum_mag-1)) && (numedges < highcount))
	{
      r++;
      numedges += hist[r];
   }
	
   highthreshold = r;
   lowthreshold = (int)(highthreshold * T_LOW + 0.5);

   if(VERBOSE)
	{
      printf("The input low and high fractions of 0.3 and 0.8 computed to\n");
      printf("magnitude of the gradient threshold : %d %d\n", lowthreshold, highthreshold);
   }

   /****************************************************************************
   * This loop looks for pixels above the highthreshold to locate edges and
   * then calls follow_edges to continue the edge.
   ****************************************************************************/
   
	for(r=0, pos=0; r < ROWS; r++)
	{
      for(c=0; c < COLUMNS; c++, pos++)
		{
	 		if((edge[pos] == POSSIBLE_EDGE) && (mag[pos] >= highthreshold))
			{
            edge[pos] = EDGE;
            follow_edges((edge+pos), (mag+pos), lowthreshold);
	 		}
      }
   }

   /****************************************************************************
   * Set all the remaining possible edges to non-edges.
   ****************************************************************************/
   
	for(r=0, pos=0; r < ROWS; r++)
	{
      for(c=0; c < COLUMNS; c++,pos++) 
			if(edge[pos] != EDGE) 
				edge[pos] = NOEDGE;
   }
}

/*******************************************************************************
* FUNCTION:		non_max_supp
* PURPOSE: 		This routine applies non-maximal suppression to the magnitude of
* 					the gradient image.
* ARGUMENTS:	short *mag
* 					short *gradx
* 					short *grady
* 					unsigned char *result				
* RETURN:		None
*******************************************************************************/
void non_max_supp(short *mag, short *gradx, short *grady, unsigned char *result) 
{
	int rowcount, colcount, count;
   short *magrowptr, *magptr;
   short *gxrowptr, *gxptr;
   short *gyrowptr, *gyptr, z1 = 0, z2 = 0;
   short m00, gx = 0, gy = 0;
   float mag1, mag2, xperp = 0, yperp = 0;
   unsigned char *resultrowptr, *resultptr;
   

   /****************************************************************************
   * Zero the edges of the result image.
   ****************************************************************************/
   
	for(count=0, resultrowptr=result, resultptr=result+COLUMNS*(ROWS-1); 
       count < COLUMNS; resultptr++, resultrowptr++, count++)
	{
    	*resultrowptr = *resultptr = (unsigned char)0;
   }

   for(count=0, resultptr=result, resultrowptr=result+COLUMNS-1;
       count < ROWS; count++, resultptr+=COLUMNS, resultrowptr+=COLUMNS)
	{
     	*resultptr = *resultrowptr = (unsigned char)0;
   }

   /****************************************************************************
   * Suppress non-maximum points.
   ****************************************************************************/
   
	for(rowcount=1, magrowptr=mag+COLUMNS+1, gxrowptr=gradx+COLUMNS+1,
      gyrowptr=grady+COLUMNS+1, resultrowptr=result+COLUMNS+1;
      rowcount < ROWS-1; /* bug fix 10/05/23, RD */
      rowcount++, magrowptr+=COLUMNS, gyrowptr+=COLUMNS, gxrowptr+=COLUMNS,
      resultrowptr+=COLUMNS)
	{   
      for(colcount=1, magptr=magrowptr, gxptr=gxrowptr, gyptr=gyrowptr,
         resultptr=resultrowptr; colcount < COLUMNS-1; /* bug fix 10/05/23, RD */
         colcount++, magptr++, gxptr++, gyptr++, resultptr++)
		{
         m00 = *magptr;
         if(m00 == 0)
			{
            *resultptr = (unsigned char)NOEDGE;
         }
         else
			{
            xperp = -(gx = *gxptr)/((float)m00);
            yperp = (gy = *gyptr)/((float)m00);
         }

         if(gx >= 0)
			{
            if(gy >= 0)
				{
            	if (gx >= gy)
               {  
                	/* 111 */
                  /* Left point */
                  z1 = *(magptr - 1);
                  z2 = *(magptr - COLUMNS - 1);

                  mag1 = (m00 - z1)*xperp + (z2 - z1)*yperp;
                       
                  /* Right point */
                  z1 = *(magptr + 1);
                  z2 = *(magptr + COLUMNS + 1);

                  mag2 = (m00 - z1)*xperp + (z2 - z1)*yperp;
	            }
   	         else
               {    
                  /* 110 */
                  /* Left point */
                  z1 = *(magptr - COLUMNS);
                  z2 = *(magptr - COLUMNS - 1);

                  mag1 = (z1 - z2)*xperp + (z1 - m00)*yperp;

                  /* Right point */
               	z1 = *(magptr + COLUMNS);
                  z2 = *(magptr + COLUMNS + 1);

                  mag2 = (z1 - z2)*xperp + (z1 - m00)*yperp; 
               }
          	}
            else
            {
               if (gx >= -gy)
             	{
               	/* 101 */
                  /* Left point */
                  z1 = *(magptr - 1);
                  z2 = *(magptr + COLUMNS - 1);

                  mag1 = (m00 - z1)*xperp + (z1 - z2)*yperp;
            
                  /* Right point */
                  z1 = *(magptr + 1);
                  z2 = *(magptr - COLUMNS + 1);

                  mag2 = (m00 - z1)*xperp + (z1 - z2)*yperp;
					}	
               else
               {    
                  /* 100 */
                  /* Left point */
                  z1 = *(magptr + COLUMNS);
                  z2 = *(magptr + COLUMNS - 1);

                  mag1 = (z1 - z2)*xperp + (m00 - z1)*yperp;

                  /* Right point */
                  z1 = *(magptr - COLUMNS);
                  z2 = *(magptr - COLUMNS + 1);

                  mag2 = (z1 - z2)*xperp  + (m00 - z1)*yperp; 
               }
       		}
      	}
         else
         {
         	if ((gy = *gyptr) >= 0)
         	{
            	if (-gx >= gy)
               {          
               	/* 011 */
                  /* Left point */
                  z1 = *(magptr + 1);
                  z2 = *(magptr - COLUMNS + 1);

                  mag1 = (z1 - m00)*xperp + (z2 - z1)*yperp;

                  /* Right point */
                  z1 = *(magptr - 1);
                  z2 = *(magptr + COLUMNS - 1);

                  mag2 = (z1 - m00)*xperp + (z2 - z1)*yperp;
					}
               else
               {
               	/* 010 */
                  /* Left point */
                  z1 = *(magptr - COLUMNS);
                  z2 = *(magptr - COLUMNS + 1);

                  mag1 = (z2 - z1)*xperp + (z1 - m00)*yperp;

                  /* Right point */
                  z1 = *(magptr + COLUMNS);
                  z2 = *(magptr + COLUMNS - 1);

                  mag2 = (z2 - z1)*xperp + (z1 - m00)*yperp;
					}
        		}
            else
            {
            	if (-gx > -gy)
               {
               	/* 001 */
                  /* Left point */
                  z1 = *(magptr + 1);
                  z2 = *(magptr + COLUMNS + 1);

                  mag1 = (z1 - m00)*xperp + (z1 - z2)*yperp;

                  /* Right point */
                  z1 = *(magptr - 1);
                  z2 = *(magptr - COLUMNS - 1);

                  mag2 = (z1 - m00)*xperp + (z1 - z2)*yperp;
	            }
               else
        	      {
                  /* 000 */
                  /* Left point */
                  z1 = *(magptr + COLUMNS);
                  z2 = *(magptr + COLUMNS + 1);

                  mag1 = (z2 - z1)*xperp + (m00 - z1)*yperp;

                  /* Right point */
                  z1 = *(magptr - COLUMNS);
                  z2 = *(magptr - COLUMNS - 1);

                  mag2 = (z2 - z1)*xperp + (m00 - z1)*yperp;
					}
        		}
       	} 

         /* Now determine if the current point is a maximum point */

         if ((mag1 > 0.0) || (mag2 > 0.0))
         {
         	*resultptr = (unsigned char) NOEDGE;
         }
         else
         {    
            if (mag2 == 0.0)
        			*resultptr = (unsigned char) NOEDGE;
            else
            	*resultptr = (unsigned char) POSSIBLE_EDGE;
         }
   	} 
 	}
}

/******************************************************************************
* FUNCTION: 	read_pgm_image
* PURPOSE: 		This function reads in an image in PGM format.All comments in the 
* 					header are discarded in the process of reading the image
* ARGUMENTS:	char *infilename		-	File containing the image data
* 					unsigned char *image	-	Buffer to store image data
* RETURN:		0	-	Failure
* 					1	-	Success
******************************************************************************/
int read_pgm_image(char *infilename, unsigned char *image)
{
   FILE *fp;
   char buf[71];

   /***************************************************************************
   * Open the input image file for reading
   ***************************************************************************/
   
	if((fp = fopen(infilename, "r")) == NULL)
	{
      fprintf(stderr, "Error reading file in read_pgm_image()\n");
      return(0);
   }

   /***************************************************************************
   * Verify that the image is in PGM format, read in the number of columns
   * and ROWS in the image and scan past all of the header information.
   ***************************************************************************/
   
	fgets(buf, 70, fp);
   if(strncmp(buf,"P5",2) != 0)
	{
      fprintf(stderr, "The file %s is not in PGM format in ", infilename);
      fprintf(stderr, "read_pgm_image().\n");
      fclose(fp);
      return(FAILED);
   }

   do{ fgets(buf, 70, fp);}while(buf[0] == '#');  /* skip all comment lines */
   do{ fgets(buf, 70, fp);}while(buf[0] == '#');  /* skip all comment lines */

   if(ROWS != fread(image, COLUMNS, ROWS, fp))
	{
      fprintf(stderr, "Error reading the image data in read_pgm_image().\n");
      fclose(fp);
      return(FAILED);
 	}

   fclose(fp);
   return(0);
}

/******************************************************************************
* FUNCTION:		write_pgm_image
* PURPOSE:		This function writes an image in PGM format
* ARGUMENTS:	char *outfilename		-	File to write edge output data
* 					unsigned char *image	-	Buffer containing edge output data
* RETURN:		0	-	Failure
* 					1	-	Success
******************************************************************************/
int write_pgm_image(char *outfilename, unsigned char *image)
{
   FILE *fp;
   
   /***************************************************************************
   * Open the output image file for writing if a filename was given. If no
   * filename was provided, set fp to write to standard output.
   ***************************************************************************/
   
   if((fp = fopen(outfilename, "w")) == NULL)
	{
         fprintf(stderr, "Error writing the file inside write_pgm_image()\n");
         return(FAILED);
   }

   /***************************************************************************
   * Write the header information to the PGM file.
   ***************************************************************************/
   
	fprintf(fp, "P5\n%d %d\n", COLUMNS, ROWS);
  	
	fprintf(fp, "%d\n", MAX_GRAY_VAL);

   /***************************************************************************
   * Write the image data to the file.
   ***************************************************************************/
   if(ROWS != fwrite(image, COLUMNS, ROWS, fp))
	{
      fprintf(stderr, "Error writing the image data in write_pgm_image().\n");
      fclose(fp);
      return(FAILED);
   }

   fclose(fp);
   return(0);
}
