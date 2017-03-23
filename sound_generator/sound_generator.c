#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define M_PI       3.14159265358979323846
#define FE 44100

double oscillateur(double amplitude, double frequence, double period, double phi)
{
	return amplitude * sin(2 * M_PI * frequence * period) + phi;
}

double modulation_amplitude(double amplitude, double frequenceC, double frequenceM, double period, double c)
{
	return (c + oscillateur(amplitude, frequenceM, period, 0)) * oscillateur(1, frequenceC, period, 0);
}

double modulation_frequence(double amplitude, double frequenceC, double frequenceM, double period, double I)
{
	return amplitude * sin(2 * M_PI * frequenceC * period + I * oscillateur(amplitude, frequenceM, period, 0));
}

double create_harmonics(double amplitude, double frequence, int nb_harmonic, double period)
{
	double value = 0;
	for(int i = 0 ; i <= nb_harmonic ; i++)
	{
		value += oscillateur(amplitude/(double)(i+1), frequence * (double)(i+1), period, 0);
	}
	return value;
}

int main(int argc, char** argv)
{
	if(argc < 6)
	{
		printf("usage : %s <raw_file_name> <duree> <channel_number> <quantificiation_bits> <amplitude>\n", argv[0]);
		exit(1);
	}

	FILE* fd = fopen(argv[1], "wb");
	int duree = atoi(argv[2]);
	int channel_number = atoi(argv[3]);
	int quantification_bits = atoi(argv[4]);

	double amplitude = atof(argv[5]);

	int size = duree * FE;
	double* samples = calloc(size, sizeof(double));

	for(int i = 0 ; i < size ; i++)
	{
		// Ajout d'harmoniques
		//samples[i] = create_harmonics(amplitude, 1000.0, 0, i * 1.0/(double)FE);

		//Synthese additive
		//samples[i] = oscillateur(amplitude, 1000.0, (double)i/FE, 0);
		//samples[i] += oscillateur(amplitude, 2000.0, (double)i/FE, 0);

		// Modulation d'amplitude
		//samples[i] = modulation_amplitude(amplitude, 1000.0, 500.0, (double)i/FE, 0.5);

		// Modulation de fréquence
		//samples[i] = modulation_frequence(amplitude, 1000.0, 200.0, (double)i/FE, 3);
		//samples[i] += modulation_frequence(amplitude, 1200.0, 150.0, (double)i/FE, 3);
		samples[i] = modulation_frequence(amplitude, 100 + i%100, 10.0, (double)i/FE, 1);

		/*
		//Fondamentale
        samples[i] = oscillateur(0.4, 440.0, (i * 1.0/(double)FE), 0);
		//1ere harmonique
		samples[i] += oscillateur(0.4/2.0, 2*440.0, (i * 1.0/(double)FE), 0);
		//2eme harmonique
		samples[i] += oscillateur(1/3.0, 3*440.0, (i * 1.0/(double)FE), 0);
		//samples[i] = 14545721;
		//samples[i] = -1 + (double)rand() / ((double)RAND_MAX/2.0);
		*/
	}

	short* to_write = calloc(size, sizeof(short));
	for(int i = 0 ; i < size ; i++)
	{
		to_write[i] = (short) (pow(2, (quantification_bits - 1)) * samples[i]);
	}

	fwrite(to_write, sizeof(short), size, fd);

	// Executing Sox Command
	char cmd[256];
	// SoX command line arguments
	// -r : rate sampling rate
	// -c : channel nummber
	// -b : byte quantificiation
	// -e : encoding
	snprintf (cmd, 256,
		"sox -b%d -esigned-integer -c%d -r%d %s %s",
		quantification_bits,
		channel_number, (int)FE, argv[1], "out.wav");

	system (cmd);

	fclose(fd);
}
