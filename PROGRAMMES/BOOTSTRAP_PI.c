// BOOTSTRAP_PI.c

#include "libhdr"

int i, r, REPS, rnd;
double w, sum_pi,  pii[1000];
struct acc PII[100000];

FILE *fin, *fout;

main()
{
	getintandskip("REPS:",&REPS,1,infinity);

	fin = fopen ("r2_PB51","r");
	fout = fopen ("outfile_pi","w");

	for (i=1; i<=965; i++)
	{
		fscanf(fin,"%lf", &w);
		fscanf(fin,"%lf", &w);
		pii[i] = w;
	}
	fclose(fin);

	for(r=1;r<REPS;r++)
	{
		sum_pi = 0.0;
		for (i=1; i<=965; i++)
		{
			rnd = (int)(uniform()*965)+1;
			sum_pi += pii[rnd];
		}
		accum(&PII, sum_pi/965);
	}

	fprintf(fout, "mean = %f    sd = %f\n", accmean(&PII), sqrt(variance(&PII)));
	printf("mean = %f    sd = %f\n", accmean(&PII), sqrt(variance(&PII)));
}

