#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include <complex.h>
#include <fftw3.h>
#include <sndfile.h>

#include <math.h>

#include "gnuplot_i.h"

#define	FRAME_SIZE 1024
#define HOP_SIZE 1024

static gnuplot_ctrl *h;
static fftw_plan plan;

static void
print_usage (char *progname)
{	printf ("\nUsage : %s <input file> \n", progname) ;
	puts ("\n"
		) ;

} 
static void
fill_buffer(double *buffer, double *new_buffer)
{
  int i;
  double tmp[FRAME_SIZE-HOP_SIZE];
  
  /* save */
  for (i=0;i<FRAME_SIZE-HOP_SIZE;i++)
    tmp[i] = buffer[i+HOP_SIZE];
  
  /* save offset */
  for (i=0;i<(FRAME_SIZE-HOP_SIZE);i++)
    {
      buffer[i] = tmp[i];
    }
  
  for (i=0;i<HOP_SIZE;i++)
    {
      buffer[FRAME_SIZE-HOP_SIZE+i] = new_buffer[i];
    }
}

static int
read_n_samples (SNDFILE * infile, double * buffer, int channels, int n)
{

  if (channels == 1)
    {
      /* MONO */
      int readcount ;

      readcount = sf_readf_double (infile, buffer, n);

      return readcount==n;
    }
  else if (channels == 2)
    {
      /* STEREO */
      double buf [2 * n] ;
      int readcount, k ;
      readcount = sf_readf_double (infile, buf, n);
      for (k = 0 ; k < readcount ; k++)
	buffer[k] = (buf [k * 2]+buf [k * 2+1])/2.0 ;

      return readcount==n;
    }
  else
    {
      /* FORMAT ERROR */
      printf ("Channel format error.\n");
    }
  
  return 0;
} 

static int
read_samples (SNDFILE * infile, double * buffer, int channels)
{
  return read_n_samples (infile, buffer, channels, HOP_SIZE);
}

void
fft_init (complex in[FRAME_SIZE], complex spec[FRAME_SIZE])
{
  plan = fftw_plan_dft_1d(FRAME_SIZE, in, spec, FFTW_FORWARD, FFTW_ESTIMATE);
}

void
fft_exit (void)
{
  fftw_destroy_plan(plan);
}

void
fft_process (void)
{
  fftw_execute(plan);
}


int
main (int argc, char * argv [])
{	char 		*progname, *infilename;
	SNDFILE	 	*infile = NULL ;
	SF_INFO	 	sfinfo ;

	progname = strrchr (argv [0], '/') ;
	progname = progname ? progname + 1 : argv [0] ;

	if (argc != 2)
	{	print_usage (progname) ;
		return 1 ;
		} ;

	infilename = argv [1] ;

	if ((infile = sf_open (infilename, SFM_READ, &sfinfo)) == NULL)
	{	printf ("Not able to open input file %s.\n", infilename) ;
		puts (sf_strerror (NULL)) ;
		return 1 ;
		} ;

	/* Read WAV */
	int nb_frames = 0;
	double new_buffer[HOP_SIZE];
	double buffer[FRAME_SIZE];

	/* Plot Init */
	h=gnuplot_init();
	gnuplot_setstyle(h, "lines");

	
	int i;
	for (i=0;i<(FRAME_SIZE/HOP_SIZE-1);i++)
	  {
	    if (read_samples (infile, new_buffer, sfinfo.channels)==1)
	      fill_buffer(buffer, new_buffer);
	    else
	      {
		printf("not enough samples !!\n");
		return 1;
	      }
	  }



	complex samples[FRAME_SIZE];
	double amplitude[FRAME_SIZE];
	complex spectrum[FRAME_SIZE];

	/* FFT init */
	fft_init(samples, spectrum);
		
	while (read_samples (infile, new_buffer, sfinfo.channels)==1)
	  {
	    /* Process Samples */
	    printf("Processing frame %d\n",nb_frames);

	    /* hop size */
	    fill_buffer(buffer, new_buffer);


	    /* ----- */
	    /* Pitch Estimation */

	    double max=0;

	    /* Pitch analysis */

	    /* Difference Function */
	    double *df = malloc(sizeof(double)* (FRAME_SIZE/2));
	  
	    int tau,j, tauok=0, trouve=0;
	    double tmp;
  
	    for (tau=0; tau<FRAME_SIZE/2; tau++)
	      {
		df[tau] = 0.;
	      }
	    
	    /* TODO */
	    /* DIFF FUNCTION */
	    
	    
	    /* cumulative mean normalized difference function */
	    tmp = 0.;
	    df[0] = 1.;
  
	    for (tau=1; tau<FRAME_SIZE/2; tau++)
	      {
		/* TODO */
		/* CMNDF */
	      }

	    /* search pitch */
	    tau=1;
  
	    do 
	      {
		if(df[tau] < 0.225) { // valeur de settings à régler
		  while (df[tau+1] < df[tau]) {
		    tau++;
		  }
		  tauok=tau;
		  trouve=1;
		}
		tau++;
	      } while (tau<FRAME_SIZE/2 && !trouve);

	    /* TODO */
	    /* print pitch */


	    
	    

	    /* ----- */
	    /* Spectral */

	    // fft process
	    for (i=0; i < FRAME_SIZE; i++)
	      {
		// Fenetre Hann
		//samples[i] = buffer[i]*(0.5-0.5*cos(2.0*M_PI*(double)i/FRAME_SIZE));
		// Fenetre rect
		samples[i] = buffer[i];
	      }
	    for (i=FRAME_SIZE; i < FRAME_SIZE; i++)
	      {
		samples[i] = 0.0;
	      }

	    
	    fft_process();

	    // spectrum contient les complexes résultats de la fft
	    for (i=0; i < FRAME_SIZE; i++)
	      {
		amplitude[i] = cabs(spectrum[i]);
	      }


	    int imax=0;
	    max = 0.0;
	    
	    for (i=0; i < FRAME_SIZE/2; i++)
	      {
		if (amplitude[i] > max)
		  {
		    max = amplitude[i];
		    imax = i;
		  }    
	      }
	    printf("Spectral : max %d %f\n",imax, (double)imax*sfinfo.samplerate/FRAME_SIZE);
	    
	    /* TODO */  
	    
	    /* plot amplitude */
	    gnuplot_resetplot(h);
	    gnuplot_plot_x(h,amplitude,FRAME_SIZE/2,"amplitude");
	    sleep(1);
    
	    /* PLOT */
	    //gnuplot_resetplot(h);
	    //gnuplot_plot_x(h,buffer,FRAME_SIZE,"temporal frame");
	    //sleep(1);
    
    	    
	    nb_frames++;
	  }

	sf_close (infile) ;

	/* FFT exit */
	fft_exit();
	
	return 0 ;
} /* main */

