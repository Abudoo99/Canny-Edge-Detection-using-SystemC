/* source: http://marathon.csee.usf.edu/edge/edge_detection.html */
/* URL: ftp://figment.csee.usf.edu/pub/Edge_Comparison/source_code/canny.src */

/* ECPS 203 Assignment 6 */

/* Author : Abyukth Kumar */

/* off-by-one bugs fixed by Rainer Doemer, 10/05/23 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "systemc.h"

#define VERBOSE 0

#define NOEDGE 255
#define POSSIBLE_EDGE 128
#define EDGE 0
#define BOOSTBLURFACTOR 90.0
#ifndef M_PI
#define M_PI 3.14159265356789202346
#endif
#define SIGMA 0.6
#define TLOW  0.3
#define THIGH 0.8

#define WINDOWSIZE	21

#define COLS 2704
#define ROWS 1520
#define SIZE COLS*ROWS
#define VIDEONAME "Engineering"
#define IMG_IN    "video/" VIDEONAME "%03d.pgm"
#define IMG_OUT   VIDEONAME "%03d_edges.pgm"
#define IMG_NUM   30 /* number of images processed (1 or more) */
#define AVAIL_IMG 30 /* number of different image frames (1 or more) */
#define SET_STACK_SIZE set_stack_size(128*1024*1024);
/* upper bound for the size of the gaussian kernel
 * SIGMA must be less than 4.0
 * check for 'windowsize' below
 */

typedef struct Image_s
{
	unsigned char img[SIZE];

	Image_s(void)
	{
	   for (int i=0; i<SIZE; i++)
	   { 
	      img[i] = 0;
	   }
	}

	Image_s& operator=(const Image_s& copy)
	{
	   for (int i=0; i<SIZE; i++)
	   { 
	      img[i] = copy.img[i];
	   }            
	   return *this;
	}

	operator unsigned char*()
	{
	   return img;
	}

	unsigned char& operator[](const int index)
	{
	   return img[index];
	}
} IMAGE;

typedef struct SImage_s
{
   short int img[SIZE];

   SImage_s(void)
   {
      for (int i=0; i<SIZE; i++)
      {
         img[i] = 0;
      }
   }

   SImage_s& operator=(const SImage_s& copy)
   {
      for (int i=0; i<SIZE; i++)
      {
         img[i] = copy.img[i];
      }
      return *this;
   }

   operator short int*()
   {
      return img;
   }

   short int& operator[](const int index)
   {
      return img[index];
   }
}SIMAGE;

typedef struct Kernel_s
{
   float kernel[WINDOWSIZE];

   Kernel_s(void)
   {
      for (int i=0; i<WINDOWSIZE; i++)
      {
         kernel[i] = 0;
      }
   }

   Kernel_s& operator=(const Kernel_s& copy)
   {
      for (int i=0; i<WINDOWSIZE; i++)
      {
         kernel[i] = copy.kernel[i];
      }
      return *this;
   }

   operator float*()
   {
      return kernel;
   }

   float& operator[](const int index)
   {
      return kernel[index];
   }
}KERNEL;

typedef struct Tempim_s
{
   float tempim[SIZE];

   Tempim_s(void)
   {
      for (int i=0; i<SIZE; i++)
      {
         tempim[i] = 0;
      }
   }

   Tempim_s& operator=(const Tempim_s& copy)
   {
      for (int i=0; i<SIZE; i++)
      {
         tempim[i] = copy.tempim[i];
      }
      return *this;
   }

   operator float*()
   {
      return tempim;
   }

   float& operator[](const int index)
   {
      return tempim[index];
   }
}TEMPIM;



SC_MODULE(Stimulus){
	IMAGE imageout;
	sc_fifo_out<IMAGE> ImgOut;
	/******************************************************************************
	* Function: read_pgm_image
	* Purpose: This function reads in an image in PGM format. The image can be
	* read in from either a file or from standard input. The image is only read
	* from standard input when infilename = NULL. Because the PGM format includes
	* the number of columns and the number of rows in the image, these are read
	* from the file. Memory to store the image is allocated OUTSIDE this function.
	* The found image size is checked against the expected rows and cols.
	* All comments in the header are discarded in the process of reading the
	* image. Upon failure, this function returns 0, upon sucess it returns 1.
	******************************************************************************/
	int read_pgm_image(const char *infilename, unsigned char *image, int rows, int cols)
	{
	   FILE *fp;
	   char buf[71];
	   int r, c;

	   /***************************************************************************
	   * Open the input image file for reading if a filename was given. If no
	   * filename was provided, set fp to read from standard input.
	   ***************************************************************************/
	   if(infilename == NULL) fp = stdin;
	   else{
	      if((fp = fopen(infilename, "r")) == NULL){
	         fprintf(stderr, "Error reading the file %s in read_pgm_image().\n",
	            infilename);
	         return(0);
	      }
	   }

	   /***************************************************************************
	   * Verify that the image is in PGM format, read in the number of columns
	   * and rows in the image and scan past all of the header information.
	   ***************************************************************************/
	   fgets(buf, 70, fp);
	   if(strncmp(buf,"P5",2) != 0){
	      fprintf(stderr, "The file %s is not in PGM format in ", infilename);
	      fprintf(stderr, "read_pgm_image().\n");
	      if(fp != stdin) fclose(fp);
	      return(0);
	   }
	   do{ fgets(buf, 70, fp); }while(buf[0] == '#');  /* skip all comment lines */
	   sscanf(buf, "%d %d", &c, &r);
	   if(c != cols || r != rows){
	      fprintf(stderr, "The file %s is not a %d by %d image in ", infilename, cols, rows);
	      fprintf(stderr, "read_pgm_image().\n");
	      if(fp != stdin) fclose(fp);
	      return(0);
	   }
	   do{ fgets(buf, 70, fp); }while(buf[0] == '#');  /* skip all comment lines */

	   /***************************************************************************
	   * Read the image from the file.
	   ***************************************************************************/
	   if((unsigned)rows != fread(image, cols, rows, fp)){
	      fprintf(stderr, "Error reading the image data in read_pgm_image().\n");
	      if(fp != stdin) fclose(fp);
	      return(0);
	   }

	   if(fp != stdin) fclose(fp);
	   return(1);
	}

	void main(void)
	{
	   int i=0, n=0;
	   char infilename[40];

	   for(i=0; i<IMG_NUM; i++)
	   {
	      n = i % AVAIL_IMG;
	      sprintf(infilename, IMG_IN, n+1);

	      /****************************************************************************
	      * Read in the image.
	      ****************************************************************************/
	      if(VERBOSE) printf("Reading the image %s.\n", infilename);
	      if(read_pgm_image(infilename, imageout, ROWS, COLS) == 0){
	         fprintf(stderr, "Error reading the input image, %s.\n", infilename);
	         exit(1);
	      }
	      ImgOut.write(imageout);
	   }
	}

	SC_CTOR(Stimulus)
	{
	   SC_THREAD(main);
	   SET_STACK_SIZE
	}
};

SC_MODULE(Monitor){
	IMAGE imagein;
	sc_fifo_in<IMAGE>  ImgIn;

	/******************************************************************************
	* Function: write_pgm_image
	* Purpose: This function writes an image in PGM format. The file is either
	* written to the file specified by outfilename or to standard output if
	* outfilename = NULL. A comment can be written to the header if coment != NULL.
	******************************************************************************/
	int write_pgm_image(const char *outfilename, unsigned char *image, int rows,
	    int cols, const char *comment, int maxval)
	{
	   FILE *fp;

	   /***************************************************************************
	   * Open the output image file for writing if a filename was given. If no
	   * filename was provided, set fp to write to standard output.
	   ***************************************************************************/
	   if(outfilename == NULL) fp = stdout;
	   else{
	      if((fp = fopen(outfilename, "w")) == NULL){
	         fprintf(stderr, "Error writing the file %s in write_pgm_image().\n",
	            outfilename);
	         return(0);
	      }
	   }

	   /***************************************************************************
	   * Write the header information to the PGM file.
	   ***************************************************************************/
	   fprintf(fp, "P5\n%d %d\n", cols, rows);
	   if(comment != NULL)
	      if(strlen(comment) <= 70) fprintf(fp, "# %s\n", comment);
	   fprintf(fp, "%d\n", maxval);

	   /***************************************************************************
	   * Write the image data to the file.
	   ***************************************************************************/
	   if((unsigned)rows != fwrite(image, cols, rows, fp)){
	      fprintf(stderr, "Error writing the image data in write_pgm_image().\n");
	      if(fp != stdout) fclose(fp);
	      return(0);
	   }

	   if(fp != stdout) fclose(fp);
	   return(1);
	}

	void main(void)
	{
	   char outfilename[128];    /* Name of the output "edge" image */
	   int i, n;

	   for(i=0; i<IMG_NUM; i++)
	   {
	      ImgIn.read(imagein);

	      /****************************************************************************
	      * Write out the edge image to a file.
	      ****************************************************************************/
	      n = i % AVAIL_IMG;
	      sprintf(outfilename, IMG_OUT, n+1);
	      if(VERBOSE) printf("Writing the edge image in the file %s.\n", outfilename);
	      if(write_pgm_image(outfilename, imagein, ROWS, COLS,"", 255) == 0){
	         fprintf(stderr, "Error writing the edge image, %s.\n", outfilename);
	         exit(1);
	      }
			printf("Edge output generated %s\n", outfilename);
	   }

	   if(VERBOSE) printf("Monitor exits simulation.\n");
	      sc_stop();	// done testing, quit the simulation
	}
	
	SC_CTOR(Monitor)
	{
	   SC_THREAD(main);
	   SET_STACK_SIZE
	}
};

SC_MODULE(DataIn)
{
	sc_fifo_in<IMAGE> ImgIn;
	sc_fifo_out<IMAGE> ImgOut;
	IMAGE Image;

	void main()
	{
	   while(1)
	   {
	      ImgIn.read(Image);
	      ImgOut.write(Image);
	   }
	}
	
	SC_CTOR(DataIn)
	{
	   SC_THREAD(main);
	   SET_STACK_SIZE
	}
};

SC_MODULE(DataOut)
{
	sc_fifo_in<IMAGE> ImgIn;
	sc_fifo_out<IMAGE> ImgOut;
	IMAGE Image;

	void main()
	{
	   while(1)
	   {
	      ImgIn.read(Image);
	      ImgOut.write(Image);
	   }
	}
	
	SC_CTOR(DataOut)
	{
	   SC_THREAD(main);
	   SET_STACK_SIZE
	}
};

SC_MODULE(Gaussian_Kernel)
{
   KERNEL kernel;

	sc_fifo_out<KERNEL> G_Out1;
	sc_fifo_out<KERNEL> G_Out2;
	sc_fifo_out<int> C_Out1;
	sc_fifo_out<int> C_Out2;

	/*******************************************************************************
	* PROCEDURE: make_gaussian_kernel
	* PURPOSE: Create a one dimensional gaussian kernel.
	* NAME: Mike Heath
	* DATE: 2/15/96
	*******************************************************************************/
	void make_gaussian_kernel(float sigma, float *kernel, int *windowsize)
	{
	   int i, center;
	   float x, fx, sum=0.0;

	   *windowsize = 1 + 2 * ceil(2.5 * sigma);
	   center = (*windowsize) / 2;

	   if(VERBOSE) printf("      The kernel has %d elements.\n", *windowsize);

	   for(i=0;i<(*windowsize);i++){
	      x = (float)(i - center);
	      fx = pow(2.71828, -0.5*x*x/(sigma*sigma)) / (sigma * sqrt(6.2831853));
	      kernel[i] = fx;
	      sum += fx;
	   }

	   for(i=0;i<(*windowsize);i++) kernel[i] /= sum;

	   if(VERBOSE){
	      printf("The filter coefficients are:\n");
	      for(i=0;i<(*windowsize);i++)
	         printf("kernel[%d] = %f\n", i, kernel[i]);
	   }
	}

	void main()
	{
		int center, windowsize;

		make_gaussian_kernel(SIGMA, kernel, &windowsize);
		center = windowsize / 2;

		while(1)
		{
			G_Out1.write(kernel);
			G_Out2.write(kernel);
			C_Out1.write(center);	// Value of variable center
			C_Out2.write(center);	// Value of variable center
		}
	}

	SC_CTOR(Gaussian_Kernel)
	{
		SC_THREAD(main);
		SET_STACK_SIZE
	}
};

SC_MODULE(BlurX)
{
	IMAGE image;
	KERNEL kernel;
	TEMPIM tempim;

	sc_fifo_in<IMAGE> ImgIn;
	sc_fifo_in<KERNEL> KernelIn;
	sc_fifo_in<int> CenterIn; 
	sc_fifo_out<TEMPIM> TempImOut;

	/*******************************************************************************
	* PROCEDURE: blur_x_image
	* PURPOSE: Blur an image in X direction.
	* NAME: Mike Heath
	* DATE: 2/15/96
	*******************************************************************************/
	void blur_x_image(unsigned char *image, float *kernel, int center, int rows, 
		int cols, float *tempim)
	{
		int r, c, rr, cc;     			/* Counter variables */
		float dot;            			/* Dot product summing variable */
		float sum;            			/* Sum of the kernel weights variable */		

		if(VERBOSE) printf("   Bluring the image in the X-direction.\n");
		for(r=0;r<rows;r++){
			for(c=0;c<cols;c++){
				dot = 0.0;
				sum = 0.0;
				for(cc=(-center);cc<=center;cc++){
					if(((c+cc) >= 0) && ((c+cc) < cols)){
					dot += (float)image[r*cols+(c+cc)] * kernel[center+cc];
					sum += kernel[center+cc];
					}
				}
				tempim[r*cols+c] = dot/sum;
			}
		}		
	}

	void main()
	{
		int center;
		while(1)
		{
			ImgIn.read(image);
			KernelIn.read(kernel);
			CenterIn.read(center);

			blur_x_image(image, kernel, center, ROWS, COLS, tempim);
	
			TempImOut.write(tempim);
		}
	}

	SC_CTOR(BlurX)
	{
		SC_THREAD(main);
		SET_STACK_SIZE
	}
};

SC_MODULE(BlurY)
{
	KERNEL kernel;
	TEMPIM tempim;
	SIMAGE smoothedim;

	sc_fifo_in<KERNEL> KernelIn;
	sc_fifo_in<int> CenterIn;
	sc_fifo_in<TEMPIM> TempImIn;
	sc_fifo_out<SIMAGE> SmoothedImOut;

	/*******************************************************************************
	* PROCEDURE: blur_y_image
	* PURPOSE: Blur an image in Y direction.
	* NAME: Mike Heath
	* DATE: 2/15/96
	*******************************************************************************/
	void blur_y_image(float *kernel, int center, float *tempim, int rows, 
		int cols, short int *smoothedim)
	{
		int r, c, rr, cc;     			/* Counter variables */
		float dot;            			/* Dot product summing variable */
		float sum;            			/* Sum of the kernel weights variable */

		if(VERBOSE) printf("   Bluring the image in the Y-direction.\n");
		for(c=0;c<cols;c++){
			for(r=0;r<rows;r++){
				sum = 0.0;
				dot = 0.0;
				for(rr=(-center);rr<=center;rr++){
					if(((r+rr) >= 0) && ((r+rr) < rows)){
					dot += tempim[(r+rr)*cols+c] * kernel[center+rr];
					sum += kernel[center+rr];
					}
				}
				smoothedim[r*cols+c] = (short int)(dot*BOOSTBLURFACTOR/sum + 0.5);
			}
		}
	}

	void main()
	{
		int center;
		while(1)
		{
			KernelIn.read(kernel);
			CenterIn.read(center);
			TempImIn.read(tempim);
			blur_y_image(kernel, center, tempim, ROWS, COLS, smoothedim);

			SmoothedImOut.write(smoothedim);
		}
	}

	SC_CTOR(BlurY)
	{
		SC_THREAD(main);
		SET_STACK_SIZE
   	}
};

SC_MODULE(Gaussian_Smooth)
{
	IMAGE image;
   	SIMAGE smoothedim;

	sc_fifo_in<IMAGE> ImgIn;
	sc_fifo_out<SIMAGE> SmoothedImOut;

	sc_fifo<KERNEL> knl_q1;
	sc_fifo<KERNEL> knl_q2;
	sc_fifo<TEMPIM> tmp_q;
	sc_fifo<int> int_q1;
	sc_fifo<int> int_q2;

	Gaussian_Kernel g_kernel;
	BlurX blurx;
	BlurY blury;

	void before_end_of_elaboration(){
		g_kernel.G_Out1.bind(knl_q1);
		g_kernel.G_Out2.bind(knl_q2);
		g_kernel.C_Out1.bind(int_q1);
		g_kernel.C_Out2.bind(int_q2);

		blurx.ImgIn.bind(ImgIn);
		blurx.KernelIn.bind(knl_q1);
		blurx.CenterIn.bind(int_q1);
		blurx.TempImOut.bind(tmp_q);

		blury.KernelIn.bind(knl_q2);
		blury.CenterIn.bind(int_q2);
		blury.TempImIn.bind(tmp_q);
		blury.SmoothedImOut.bind(SmoothedImOut);
	}

	SC_CTOR(Gaussian_Smooth)
	:knl_q1("knl_q1", 2) 
	,knl_q2("knl_q2", 2) 
	,tmp_q("tmp_q", 2)
	,int_q1("int_q1", 2) 
	,int_q2("int_q2", 2) 
	,g_kernel("g_kernel")
	,blurx("blurx") 
	,blury("blury")
	{}
};

SC_MODULE(Derivative_X_Y)
{
	SIMAGE smoothedim, delta_x, delta_y;
	sc_fifo_in<SIMAGE> SmoothedImIn;
	sc_fifo_out<SIMAGE> XOut1;
	sc_fifo_out<SIMAGE> XOut2;
	sc_fifo_out<SIMAGE> YOut1;
	sc_fifo_out<SIMAGE> YOut2;

	/*******************************************************************************
	* PROCEDURE: derivative_x_y
	* PURPOSE: Compute the first derivative of the image in both the x any y
	* directions. The differential filters that are used are:
	*
	*                                          -1
	*         dx =  -1 0 +1     and       dy =  0
	*                                          +1
	*
	* NAME: Mike Heath
	* DATE: 2/15/96
	*******************************************************************************/
	void derivative_x_y(short int *smoothedim, int rows, int cols,
	        short int *delta_x, short int *delta_y)
	{
	   int r, c, pos;

		/****************************************************************************
		 * Compute the x-derivative. Adjust the derivative at the borders to avoid
		 * losing pixels.
		 ****************************************************************************/
		if(VERBOSE) printf("   Computing the X-direction derivative.\n");
		for(r=0;r<rows;r++){
			pos = r * cols;
			delta_x[pos] = smoothedim[pos+1] - smoothedim[pos];
			pos++;
			for(c=1;c<(cols-1);c++,pos++){
				delta_x[pos] = smoothedim[pos+1] - smoothedim[pos-1];
			}
			delta_x[pos] = smoothedim[pos] - smoothedim[pos-1];
		}

		/****************************************************************************
		 * Compute the y-derivative. Adjust the derivative at the borders to avoid
		 * losing pixels.
		 ****************************************************************************/
		if(VERBOSE) printf("   Computing the Y-direction derivative.\n");
		for(c=0;c<cols;c++){
			pos = c;
			delta_y[pos] = smoothedim[pos+cols] - smoothedim[pos];
			pos += cols;
			for(r=1;r<(rows-1);r++,pos+=cols){
				delta_y[pos] = smoothedim[pos+cols] - smoothedim[pos-cols];
			}
			delta_y[pos] = smoothedim[pos] - smoothedim[pos-cols];
		}
	}
	
	void main()
	{
		while(1)
		{
			SmoothedImIn.read(smoothedim);
			
			derivative_x_y(smoothedim, ROWS, COLS, delta_x, delta_y);
			
			XOut1.write(delta_x);
			YOut1.write(delta_y);
			XOut2.write(delta_x);
			YOut2.write(delta_y);
		}
	}

   SC_CTOR(Derivative_X_Y)
   {
      SC_THREAD(main);   
      SET_STACK_SIZE
   }
};

SC_MODULE(Magnitude_X_Y)
{
	SIMAGE delta_x, delta_y, magnitude;

	sc_fifo_in<SIMAGE> XIn;
	sc_fifo_in<SIMAGE> YIn;
	sc_fifo_out<SIMAGE> MagnitudeOut1;
	sc_fifo_out<SIMAGE> MagnitudeOut2;

	/*******************************************************************************
	* PROCEDURE: magnitude_x_y
	* PURPOSE: Compute the magnitude of the gradient. This is the square root of
	* the sum of the squared derivative values.
	* NAME: Mike Heath
	* DATE: 2/15/96
	*******************************************************************************/
	void magnitude_x_y(short int *delta_x, short int *delta_y, int rows, int cols,
	        short int *magnitude)
	{
	   int r, c, pos, sq1, sq2;

	   for(r=0,pos=0;r<rows;r++){
	      for(c=0;c<cols;c++,pos++){
	         sq1 = (int)delta_x[pos] * (int)delta_x[pos];
	         sq2 = (int)delta_y[pos] * (int)delta_y[pos];
	         magnitude[pos] = (short)(0.5 + sqrt((float)sq1 + (float)sq2));
	      }
	   }

	}

	void main()
	{
		while(1)
		{
			XIn.read(delta_x);
			YIn.read(delta_y);

			magnitude_x_y(delta_x, delta_y, ROWS, COLS, magnitude);

			MagnitudeOut1.write(magnitude);
			MagnitudeOut2.write(magnitude);
		}
	}

	SC_CTOR(Magnitude_X_Y)
	{
		SC_THREAD(main);
		SET_STACK_SIZE
	}
};

SC_MODULE(Non_Max_Supp)
{
	SIMAGE magnitude, delta_x, delta_y;
	IMAGE nms;

	sc_fifo_in<SIMAGE> GradxIn;
	sc_fifo_in<SIMAGE> GradyIn;
	sc_fifo_in<SIMAGE> MagIn;
	sc_fifo_out<IMAGE> NmsOut;

	/*******************************************************************************
	* PROCEDURE: non_max_supp
	* PURPOSE: This routine applies non-maximal suppression to the magnitude of
	* the gradient image.
	* NAME: Mike Heath
	* DATE: 2/15/96
	*******************************************************************************/
	void non_max_supp(short *mag, short *gradx, short *grady, int nrows, int ncols,
	    unsigned char *result)
	{
	    int rowcount, colcount,count;
	    short *magrowptr,*magptr;
	    short *gxrowptr,*gxptr;
	    short *gyrowptr,*gyptr,z1,z2;
	    short m00; short gx=0; short gy=0;
	    float mag1, mag2; float xperp=0;float yperp=0;
	    unsigned char *resultrowptr, *resultptr;

	   /****************************************************************************
	   * Zero the edges of the result image.
	   ****************************************************************************/
	    for(count=0,resultrowptr=result,resultptr=result+ncols*(nrows-1);
	        count<ncols; resultptr++,resultrowptr++,count++){
	        *resultrowptr = *resultptr = (unsigned char) 0;
	    }

	    for(count=0,resultptr=result,resultrowptr=result+ncols-1;
	        count<nrows; count++,resultptr+=ncols,resultrowptr+=ncols){
	        *resultptr = *resultrowptr = (unsigned char) 0;
	    }

	   /****************************************************************************
	   * Suppress non-maximum points.
	   ****************************************************************************/
	   for(rowcount=1,magrowptr=mag+ncols+1,gxrowptr=gradx+ncols+1,
	      gyrowptr=grady+ncols+1,resultrowptr=result+ncols+1;
      rowcount<nrows-1; /* bug fix 10/05/23, RD */
	      rowcount++,magrowptr+=ncols,gyrowptr+=ncols,gxrowptr+=ncols,
	      resultrowptr+=ncols){
	      for(colcount=1,magptr=magrowptr,gxptr=gxrowptr,gyptr=gyrowptr,
         resultptr=resultrowptr;colcount<ncols-1; /* bug fix 10/05/23, RD */
	         colcount++,magptr++,gxptr++,gyptr++,resultptr++){
	         m00 = *magptr;
	         if(m00 == 0){
	            *resultptr = (unsigned char) NOEDGE;
	         }
	         else{
	            xperp = -(gx = *gxptr)/((float)m00);
	            yperp = (gy = *gyptr)/((float)m00);
	         }

	         if(gx >= 0){
	            if(gy >= 0){
	                    if (gx >= gy)
	                    {
	                        /* 111 */
	                        /* Left point */
	                        z1 = *(magptr - 1);
	                        z2 = *(magptr - ncols - 1);

	                        mag1 = (m00 - z1)*xperp + (z2 - z1)*yperp;

	                        /* Right point */
	                        z1 = *(magptr + 1);
	                        z2 = *(magptr + ncols + 1);

	                        mag2 = (m00 - z1)*xperp + (z2 - z1)*yperp;
	                    }
	                    else
	                    {
	                        /* 110 */
	                        /* Left point */
	                        z1 = *(magptr - ncols);
	                        z2 = *(magptr - ncols - 1);

	                        mag1 = (z1 - z2)*xperp + (z1 - m00)*yperp;

	                        /* Right point */
	                        z1 = *(magptr + ncols);
	                        z2 = *(magptr + ncols + 1);

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
	                        z2 = *(magptr + ncols - 1);

	                        mag1 = (m00 - z1)*xperp + (z1 - z2)*yperp;

	                        /* Right point */
	                        z1 = *(magptr + 1);
	                        z2 = *(magptr - ncols + 1);

	                        mag2 = (m00 - z1)*xperp + (z1 - z2)*yperp;
	                    }
	                    else
	                    {
	                        /* 100 */
	                        /* Left point */
	                        z1 = *(magptr + ncols);
	                        z2 = *(magptr + ncols - 1);

	                        mag1 = (z1 - z2)*xperp + (m00 - z1)*yperp;

	                        /* Right point */
	                        z1 = *(magptr - ncols);
	                        z2 = *(magptr - ncols + 1);

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
	                        z2 = *(magptr - ncols + 1);

	                        mag1 = (z1 - m00)*xperp + (z2 - z1)*yperp;

	                        /* Right point */
	                        z1 = *(magptr - 1);
	                        z2 = *(magptr + ncols - 1);

	                        mag2 = (z1 - m00)*xperp + (z2 - z1)*yperp;
	                    }
	                    else
	                    {
	                        /* 010 */
	                        /* Left point */
	                        z1 = *(magptr - ncols);
	                        z2 = *(magptr - ncols + 1);

	                        mag1 = (z2 - z1)*xperp + (z1 - m00)*yperp;

	                        /* Right point */
	                        z1 = *(magptr + ncols);
	                        z2 = *(magptr + ncols - 1);

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
	                        z2 = *(magptr + ncols + 1);

	                        mag1 = (z1 - m00)*xperp + (z1 - z2)*yperp;

	                        /* Right point */
	                        z1 = *(magptr - 1);
	                        z2 = *(magptr - ncols - 1);

	                        mag2 = (z1 - m00)*xperp + (z1 - z2)*yperp;
	                    }
	                    else
	                    {
	                        /* 000 */
	                        /* Left point */
	                        z1 = *(magptr + ncols);
	                        z2 = *(magptr + ncols + 1);

	                        mag1 = (z2 - z1)*xperp + (m00 - z1)*yperp;

	                        /* Right point */
	                        z1 = *(magptr - ncols);
	                        z2 = *(magptr - ncols - 1);

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

   	void main()
   	{
		while(1)
		{
			GradxIn.read(delta_x);
			GradyIn.read(delta_y);
			MagIn.read(magnitude);

			non_max_supp(magnitude, delta_x, delta_y, ROWS, COLS, nms);

			NmsOut.write(nms);
		}
   	}

	SC_CTOR(Non_Max_Supp)
	{
		SC_THREAD(main);
		SET_STACK_SIZE	
	}
};

SC_MODULE(Apply_Hysteresis)
{
	IMAGE nms, edge;
	SIMAGE magnitude;

	sc_fifo_in<SIMAGE> MagIn;
	sc_fifo_in<IMAGE> NmsIn;
	sc_fifo_out<IMAGE> ImgOut;
	
	/*******************************************************************************
	* PROCEDURE: follow_edges
	* PURPOSE: This procedure edges is a recursive routine that traces edgs along
	* all paths whose magnitude values remain above some specifyable lower
	* threshhold.
	* NAME: Mike Heath
	* DATE: 2/15/96
	*******************************************************************************/
	void follow_edges(unsigned char *edgemapptr, short *edgemagptr, short lowval,
	   int cols)
	{
	   short *tempmagptr;
	   unsigned char *tempmapptr;
	   int i;
	   int x[8] = {1,1,0,-1,-1,-1,0,1},
	       y[8] = {0,1,1,1,0,-1,-1,-1};

	   for(i=0;i<8;i++){
	      tempmapptr = edgemapptr - y[i]*cols + x[i];
	      tempmagptr = edgemagptr - y[i]*cols + x[i];

	      if((*tempmapptr == POSSIBLE_EDGE) && (*tempmagptr > lowval)){
	         *tempmapptr = (unsigned char) EDGE;
	         follow_edges(tempmapptr,tempmagptr, lowval, cols);
	      }
	   }
	}

	

	/*******************************************************************************
	* PROCEDURE: apply_hysteresis
	* PURPOSE: This routine finds edges that are above some high threshhold or
	* are connected to a high pixel by a path of pixels greater than a low
	* threshold.
	* NAME: Mike Heath
	* DATE: 2/15/96
	*******************************************************************************/
	void apply_hysteresis(short int *mag, unsigned char *nms, int rows, int cols,
		float tlow, float thigh, unsigned char *edge)
	{
	   int r, c, pos, numedges, highcount, lowthreshold, highthreshold, hist[32768];
	   short int maximum_mag=0;

	   /****************************************************************************
	   * Initialize the edge map to possible edges everywhere the non-maximal
	   * suppression suggested there could be an edge except for the border. At
	   * the border we say there can not be an edge because it makes the
	   * follow_edges algorithm more efficient to not worry about tracking an
	   * edge off the side of the image.
	   ****************************************************************************/
	   for(r=0,pos=0;r<rows;r++){
	      for(c=0;c<cols;c++,pos++){
		 if(nms[pos] == POSSIBLE_EDGE) edge[pos] = POSSIBLE_EDGE;
		 else edge[pos] = NOEDGE;
	      }
	   }

	   for(r=0,pos=0;r<rows;r++,pos+=cols){
	      edge[pos] = NOEDGE;
	      edge[pos+cols-1] = NOEDGE;
	   }
	   pos = (rows-1) * cols;
	   for(c=0;c<cols;c++,pos++){
	      edge[c] = NOEDGE;
	      edge[pos] = NOEDGE;
	   }

	   /****************************************************************************
	   * Compute the histogram of the magnitude image. Then use the histogram to
	   * compute hysteresis thresholds.
	   ****************************************************************************/
	   for(r=0;r<32768;r++) hist[r] = 0;
	   for(r=0,pos=0;r<rows;r++){
	      for(c=0;c<cols;c++,pos++){
		 if(edge[pos] == POSSIBLE_EDGE) hist[mag[pos]]++;
	      }
	   }

	   /****************************************************************************
	   * Compute the number of pixels that passed the nonmaximal suppression.
	   ****************************************************************************/
	   for(r=1,numedges=0;r<32768;r++){
	      if(hist[r] != 0) maximum_mag = r;
	      numedges += hist[r];
	   }

	   highcount = (int)(numedges * thigh + 0.5);

	   /****************************************************************************
	   * Compute the high threshold value as the (100 * thigh) percentage point
	   * in the magnitude of the gradient histogram of all the pixels that passes
	   * non-maximal suppression. Then calculate the low threshold as a fraction
	   * of the computed high threshold value. John Canny said in his paper
	   * "A Computational Approach to Edge Detection" that "The ratio of the
	   * high to low threshold in the implementation is in the range two or three
	   * to one." That means that in terms of this implementation, we should
	   * choose tlow ~= 0.5 or 0.33333.
	   ****************************************************************************/
	   r = 1;
	   numedges = hist[1];
	   while((r<(maximum_mag-1)) && (numedges < highcount)){
	      r++;
	      numedges += hist[r];
	   }
	   highthreshold = r;
	   lowthreshold = (int)(highthreshold * tlow + 0.5);

	   if(VERBOSE){
	      printf("The input low and high fractions of %f and %f computed to\n",
		 tlow, thigh);
	      printf("magnitude of the gradient threshold values of: %d %d\n",
		 lowthreshold, highthreshold);
	   }

	   /****************************************************************************
	   * This loop looks for pixels above the highthreshold to locate edges and
	   * then calls follow_edges to continue the edge.
	   ****************************************************************************/
	   for(r=0,pos=0;r<rows;r++){
	      for(c=0;c<cols;c++,pos++){
		 		if((edge[pos] == POSSIBLE_EDGE) && (mag[pos] >= highthreshold)){
	            edge[pos] = EDGE;
	            follow_edges((edge+pos), (mag+pos), lowthreshold, cols);
				}
	   	}
	   }

	   /****************************************************************************
	   * Set all the remaining possible edges to non-edges.
	   ****************************************************************************/
	   for(r=0,pos=0;r<rows;r++){
	      for(c=0;c<cols;c++,pos++) if(edge[pos] != EDGE) edge[pos] = NOEDGE;
	   }
	}
	
	void main()
	{
		while(1)
		{
			NmsIn.read(nms);
			MagIn.read(magnitude);

			apply_hysteresis(magnitude, nms, ROWS, COLS, TLOW, THIGH, edge);

			ImgOut.write(edge);
		}
	}

	SC_CTOR(Apply_Hysteresis)
	{
		SC_THREAD(main);
		SET_STACK_SIZE
	}
};

SC_MODULE(DUT)
{
	IMAGE image, edge;

	sc_fifo_in<IMAGE> ImgIn;
	sc_fifo_out<IMAGE> ImgOut;

	sc_fifo<SIMAGE> q1;
	sc_fifo<SIMAGE> q2;
	sc_fifo<SIMAGE> q3;
	sc_fifo<SIMAGE> q4;
	sc_fifo<SIMAGE> q5;
	sc_fifo<SIMAGE> q6;
	sc_fifo<SIMAGE> q7;
	sc_fifo<IMAGE> nms;

   Gaussian_Smooth gaussian_smooth;
   Derivative_X_Y derivative_x_y;
   Magnitude_X_Y magnitude_x_y;
   Non_Max_Supp non_max_supp;
   Apply_Hysteresis apply_hysteresis;	

	void before_end_of_elaboration(){
		gaussian_smooth.ImgIn.bind(ImgIn);
     	gaussian_smooth.SmoothedImOut.bind(q1);

		derivative_x_y.SmoothedImIn.bind(q1);
		derivative_x_y.XOut1.bind(q2);
		derivative_x_y.XOut2.bind(q3);	
		derivative_x_y.YOut1.bind(q4);
		derivative_x_y.YOut2.bind(q5);

		magnitude_x_y.XIn.bind(q2);
		magnitude_x_y.YIn.bind(q4);
		magnitude_x_y.MagnitudeOut1.bind(q6);
		magnitude_x_y.MagnitudeOut2.bind(q7);

		non_max_supp.GradxIn.bind(q3);
		non_max_supp.GradyIn.bind(q5);
		non_max_supp.MagIn.bind(q6);
		non_max_supp.NmsOut.bind(nms);		

		apply_hysteresis.MagIn.bind(q7);
		apply_hysteresis.NmsIn.bind(nms);
		apply_hysteresis.ImgOut.bind(ImgOut);				
	}
	
	SC_CTOR(DUT)
	:q1("q1", 2)
	,q2("q2", 2)
	,q3("q3", 2)
	,q4("q4", 2)
	,q5("q5", 2)
	,q6("q6", 2)
	,q7("q7", 2)
	,nms("nms", 2)
	,gaussian_smooth("gaussian_smooth")
	,derivative_x_y("derivative_x_y")
	,magnitude_x_y("magnitude_x_y")
	,non_max_supp("non_max_supp")
	,apply_hysteresis("apply_hysteresis")
	{}
};

SC_MODULE(Platform)
{
	sc_fifo_in<IMAGE> ImgIn;
	sc_fifo_out<IMAGE> ImgOut;
	sc_fifo<IMAGE> q1;
	sc_fifo<IMAGE> q2;

	DataIn din;
	DUT canny;
	DataOut dout;

	void before_end_of_elaboration(){
	   din.ImgIn.bind(ImgIn);
	   din.ImgOut.bind(q1);
	   canny.ImgIn.bind(q1);
	   canny.ImgOut.bind(q2);
	   dout.ImgIn.bind(q2);
	   dout.ImgOut.bind(ImgOut);
	}

	SC_CTOR(Platform)
	:q1("q1",2)
	,q2("q2",2)
	,din("din")
	,canny("canny")
	,dout("dout")
	{}
};

SC_MODULE(Top)
{
	sc_fifo<IMAGE> q1;
	sc_fifo<IMAGE> q2;
	Stimulus stimulus;
	Platform platform;
	Monitor monitor;

	void before_end_of_elaboration(){
	   stimulus.ImgOut.bind(q1);
	   platform.ImgIn.bind(q1);
	   platform.ImgOut.bind(q2);
	   monitor.ImgIn.bind(q2);
	}

	SC_CTOR(Top)
	:q1("q1",2)
	,q2("q2",2)
	,stimulus("stimulus")
	,platform("platform")
	,monitor("monitor")
	{}
};

Top top("top");

int sc_main(int argc, char* argv[])
{
	sc_start();
	
	return 0;
}
