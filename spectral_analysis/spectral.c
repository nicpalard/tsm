#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include <complex.h>
#include <fftw3.h>
#include <sndfile.h>

#include <math.h>

#include "gnuplot_i.h"

#define SIZE 1024
#define	FRAME_SIZE 1024
#define HOP_SIZE 1024
#define N 1024
#define FECH 44100

static gnuplot_ctrl *h;
static fftw_plan plan;

static void
print_usage (char *progname)
{
  printf ("\nUsage : %s <input file> \n", progname);
  puts ("\n");
}

static void
fill_buffer(double *buffer, double *new_buffer)
{
  int i;
  double tmp[FRAME_SIZE-HOP_SIZE];

  /* save */
  for (i = 0; i < (FRAME_SIZE - HOP_SIZE); i++) {
    tmp[i] = buffer[i + HOP_SIZE];
  }

  /* save offset */
  for (i = 0; i < (FRAME_SIZE-HOP_SIZE); i++) {
    buffer[i] = tmp[i];
  }

  for (i = 0; i < HOP_SIZE; i++) {
    buffer[FRAME_SIZE-HOP_SIZE+i] = new_buffer[i];
  }
}

static int
read_n_samples (SNDFILE * infile, double * buffer, int channels, int n)
{

  if (channels == 1) {
    /* MONO */
    int readcount;

    readcount = sf_readf_double (infile, buffer, n);

    return readcount==n;
  }
  else if (channels == 2) {
    /* STEREO */
    double buf [2 * n];
    int readcount, k;
    readcount = sf_readf_double (infile, buf, n);
    for (k = 0; k < readcount; k++) {
      buffer[k] = (buf[k*2] + buf[k*2+1]) / 2.0;
    }
    return (readcount == n);
  }
  else {
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

double hann_frame(double value)
{
	return 0.5 - 0.5 * cos( (2 * M_PI * value) / N);
}

double barret_frame(double value)
{
	return 2 - (2 * value) / N;
}

void
fft_init (complex in[SIZE], complex spec[SIZE])
{
  plan = fftw_plan_dft_1d(SIZE, in, spec, FFTW_FORWARD, FFTW_ESTIMATE);
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

void
dft_process(double s[FRAME_SIZE], complex S[FRAME_SIZE])
{
  int i,j;

  for (i = 0; i < FRAME_SIZE; i++) {
    S[i] = 0;
    for (j = 0; j < FRAME_SIZE; j++) {
      S[i] += s[j]*cexp(2.0*M_PI*I*(double)i*j/FRAME_SIZE);
    }
  }
}

int
main (int argc, char * argv [])
{
  char *progname, *infilename;
  SNDFILE *infile = NULL;
  SF_INFO sfinfo;

  progname = strrchr (argv [0], '/');
  progname = progname ? progname + 1 : argv [0];

  if (argc != 2) {
    print_usage (progname);
    return 1;
  };

  infilename = argv [1];

  if ((infile = sf_open (infilename, SFM_READ, &sfinfo)) == NULL) {
    printf ("Not able to open input file %s.\n", infilename);
    puts (sf_strerror (NULL));
    return 1;
  };

  /* Read WAV */
  int nb_frames = 0;
  double new_buffer[HOP_SIZE];
  double buffer[FRAME_SIZE];

  /* Plot Init */
  h=gnuplot_init();
  gnuplot_setstyle(h, "lines");

  int i;
  for (i = 0; i < (FRAME_SIZE/HOP_SIZE-1); i++) {
    if (read_samples (infile, new_buffer, sfinfo.channels) == 1) {
      fill_buffer(buffer, new_buffer);
    }
    else {
      printf("not enough samples !!\n");
      return 1;
    }
  }

  complex samples[SIZE];
  double amplitude[FRAME_SIZE];
  double phase[FRAME_SIZE];
  complex spectrum[SIZE];

  /* FFT init */
  fft_init(samples, spectrum);


  while (read_samples (infile, new_buffer, sfinfo.channels) == 1) {
    /* Process Samples */
    printf("Processing frame %d\n", nb_frames);

    /* hop size */
	fill_buffer(buffer, new_buffer);
	//fill_buffer(zero_padd, new_buffer);

    for (int i = 0; i < SIZE; i++) {
      samples[i] = buffer[i] * hann_frame(i);
  	}

    /* FFT process */
    fftw_execute(plan);
	// Spectre d'amplitude symétrique -> moitié de frame
    double max_amplitude = cabs(spectrum[0]);
    int frequence = 0;
    for (int i = 0; i < FRAME_SIZE/2; i++) {
      amplitude[i] = cabs(spectrum[i]);
      if (amplitude[i] > max_amplitude)
	  {
		  max_amplitude = amplitude[i];
		  frequence = i;
      }
    }

    //max_amplitude /= (FRAME_SIZE/2);
    printf("  max amplitude : %f\n", max_amplitude);
    //frequence = (frequence * FECH) / FRAME_SIZE;
    printf("  fréquence : %d\n", frequence);

    float al = 20 * log10(amplitude[frequence - 1]);
    float ac = 20 * log10(amplitude[frequence]);
    float ar = 20 * log10(amplitude[frequence + 1]);

    float d = 0.5 * ((al - ar)/(al - 2 * ac + ar));

    float frequence2 = ((frequence + d) * FECH) / FRAME_SIZE;
    printf("  max amplitude : %f\n", max_amplitude);
    printf("  fréquence : %f\n", frequence2);
	printf("  note : %d\n", (int) round(57 + 12 * log2(frequence2/440.0))%12);


    /* plot amplitude */
    /* PLOT */
    gnuplot_resetplot(h);
    //gnuplot_plot_x(h, buffer, FRAME_SIZE, "temporal frame");
    //gnuplot_plot_x(h, amplitude, FRAME_SIZE/2, "spectral");
    //sleep(1);

    nb_frames++;
  }

  sf_close (infile);

  /* FFT exit */
  fft_exit();

  return 0;
} /* main */
