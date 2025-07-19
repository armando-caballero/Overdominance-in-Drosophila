/* SIMOVERDOM_ines17.c (28/03/2025) */

/* ***************************************************** */

#include "libhdr"
#define NN 10001  /* max number of NIND */
#define MM 8001  /* max number of NCRO */
#define NW 301  /* max number of 100 kb windows */
#define SS 100000  /* max number of SNPs */
#define GG 4000  /* max number of genes */

int NIND, NCRO, NLOCI, TOTLOCI, numSNPs, numSNPsneu, numSNPsdel, numSNPslet, numSNPsod, numSNPsqn, numSNPsadv, numSNPNP;
int i, j, k, l, rep, classes, replicates, g, b;
int window_size, nwindows;
int RM[NN], ran_i;
int gm[NN][MM][2];
int ss, cc, crom[SS], loc[SS], genes, bb;

int NINDNP, gmNP[NN][MM][2];
int chromNP[MM][31], chrom[MM][31];
unsigned long long int x, posNP[MM][31], pos[MM][31], posplink[SS], startgene[GG], endgene[GG], lengthgene[200000], posb[300000];
double sNP[MM][31], atNP[MM][31], hsNP[MM][31], hatNP[MM][31], qNP[MM][31];
double s[MM][31], s1[MM][31], s2[MM][31], s1w[MM][31], s2w[MM][31], at[MM][31], hs[MM][31], hat[MM][31], freqNP[MM][31], LEQ_NP;
double r2[SS], pii[SS], B[300000];

double pm_s[NN], pm_f[NN];

double d_a_del, alfa_a_del, va_del, vd_del, id_del;
double d_a_delw, alfa_a_delw, va_delw[NW], vd_delw[NW], id_delw[NW];
double d_a_let, alfa_a_let, va_let, vd_let, id_let;
double d_a_letw, alfa_a_letw, va_letw[NW], vd_letw[NW], id_letw[NW];
double d_a_od, alfa_a_od, va_od, vd_od, id_od;
double d_a_odw, alfa_a_odw, va_odw[NW], vd_odw[NW], id_odw[NW];
double d_a_qn, alfa_a_qn, va_qn, vd_qn, id_qn;
double d_a_qnw, alfa_a_qnw, va_qnw[NW], vd_qnw[NW], id_qnw[NW];
double d_a_adv, alfa_a_adv, va_adv, vd_adv, id_adv;
double d_a_advw, alfa_a_advw, va_advw[NW], vd_advw[NW], id_advw[NW];

double AA, Aa, aa, q[MM][31];
double FITN[NN], FITN_F[NN];
double mean_FITN, mean_FITN_F;

double numTajDwin[NW], TajDwin[NW], cwin[NW];

struct acc gmean_s, gvar_s, AVE_VA, AVE_VD, AVE_ID, AVE_VA_del, AVE_VD_del, AVE_ID_del, AVE_VA_let, AVE_VD_let, AVE_ID_let, AVE_VA_od, AVE_VD_od, AVE_ID_od, AVE_VA_qn, AVE_VD_qn, AVE_ID_qn, AVE_VA_adv, AVE_VD_adv, AVE_ID_adv;

struct acc numr2win[NW], r2win[NW], numpiwin[NW], piwin[NW], neuwin[NW], delwin[NW], letwin[NW], odwin[NW], hsdelwin[NW], sdelwin[NW], hsletwin[NW], sletwin[NW], hsodwin[NW], sodwin[NW], qletwin[NW], qodwin[NW], qdelwin[NW], qneuwin[NW], numgenes[NW], length[NW], Bswin[NW], qnwin[NW], qqnwin[NW],sqnwin[NW], hsqnwin[NW], advwin[NW], qadvwin[NW], sadvwin[NW], hsadvwin[NW];

struct acc AVE_numr2win[NW], AVE_r2win[NW], AVE_numpiwin[NW], AVE_piwin[NW], AVE_neuwin[NW], AVE_delwin[NW], AVE_letwin[NW], AVE_odwin[NW], AVE_qnwin[NW], AVE_advwin[NW], AVE_va_delwin[NW], AVE_vd_delwin[NW], AVE_id_delwin[NW], AVE_va_letwin[NW], AVE_vd_letwin[NW], AVE_id_letwin[NW], AVE_va_odwin[NW], AVE_vd_odwin[NW], AVE_id_odwin[NW], AVE_va_qnwin[NW], AVE_vd_qnwin[NW], AVE_id_qnwin[NW], AVE_hsdelwin[NW], AVE_sdelwin[NW], AVE_hsletwin[NW], AVE_sletwin[NW], AVE_hsodwin[NW], AVE_sodwin[NW], AVE_hsqnwin[NW], AVE_sqnwin[NW], AVE_hsadvwin[NW], AVE_sadvwin[NW], AVE_qdelwin[NW], AVE_qletwin[NW], AVE_qodwin[NW], AVE_qneuwin[NW], AVE_qqnwin[NW], AVE_qadvwin[NW], AVE_va_advwin[NW], AVE_vd_advwin[NW], AVE_id_advwin[NW];

struct acc AVE_numSNPsneu[SS], AVE_numSNPsdel[SS], AVE_numSNPslet[SS], AVE_numSNPsod[SS], AVE_numSNPsqn[SS], AVE_numSNPsadv[SS];
    
struct acc AVE_numTajDwin[NW], AVE_TajDwin[NW], AVE_cwin[NW];

struct acc q_0[NN];
struct acc AVE_q_0[NN];

struct acc AVE_mean_FITN, AVE_mean_FITN_F, AVE_ID_FITN;

FILE *fptr, *fgen, *frep, *fdat, *fpop, *ffmap, *ffreq, *ffped, *fr2file, *fpifile, *fTajDfile, *fRec, *fgenes, *fB, *ftable;

/* ***************************************************** */

main()
{
	fptr = fopen ("dfilename.dat","w");
	fgen = fopen ("genfile.dat","w");
	ffreq = fopen ("frequencies.dat","w");
    	ftable = fopen ("table.txt","w");

	getinputs();
	recombination_masks();
	natural_population();
	read_genes();
	read_Bs();

	fprintf(ffreq, "0.0 ");
	for (j=1; j<classes; j++)	fprintf(ffreq, "%4.2f ", (double)j/100);

	for (rep=1; rep<=replicates; rep++)
	{
		numSNPs = 0.0;

		frep = fopen ("repfile.dat","a");
		fprintf (frep,"\nreplicate %d\n", rep);
		fclose(frep);

		fprintf(ffreq, "\n");

	//	if (tracelevel!=0) fprintf (fptr,"\n\nreplicate %d\n\n", rep);

		sample();
		frequency_genes();
		genotypic_values();
		phenotypeB();
//		if (tracelevel!=0) dumpphenotypes();
		plink_files();
		int status1 = system("bash shell_calculations");
		READ_FILES();
		ID_HOMOZYGOTES();
		settozero();
	 	printout();
		//int status2 = system("graphs");
		writeseed();
  	}
}

/* ***************************************************** */

getinputs()
{
	tracestart();
	getseed();

	getintandskip("NINDNP (max 10000):",&NINDNP,2,10000);
	getintandskip("NIND (max 1000):",&NIND,2,10000);
	NLOCI=30;
	getintandskip("Number of frequeny classes :",&classes,1,1000);
	getintandskip("NWINDOWS:",&nwindows,1,infinity);
	getintandskip("Window size :",&window_size,1,100000);
	getintandskip("Number of replicates :",&replicates,1,infinity);
}

/* **************************************************** */

recombination_masks ()
{
	for (l=0; l<NLOCI; l++)   RM[l]=pow(2.0,(double)l);
}

/* ***************************************************** */

natural_population ()
{
	int x, ds, g0, g1;
	unsigned long long int dpos;

	double dps, da, dh, dq;

	/* ***** take effects of genes ***** */

	fdat=fopen("list_allsnps","r");

	fscanf(fdat,"%d", &x);
	numSNPNP = x;

	NCRO = numSNPNP/NLOCI;
	TOTLOCI = NCRO * NLOCI;

	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	{
		fscanf(fdat,"%d%llu%lf%lf%lf%lf", &ds, &dpos, &dps, &da, &dh, &dq);
		chromNP[k][l] = ds;
		//fprintf(fptr,"chromNP[%d][%d]=%d\n", k, l, chromNP[k][l]);
		posNP[k][l] = dpos;
		if (dps < -1.0) dps=(-1.0);
		if (da == -99.0) da=0.0;
		sNP[k][l] = dps;
		atNP[k][l] = da;
		hsNP[k][l] = dh;
		if (sNP[k][l] == -1.0) hsNP[k][l]=0.02;
		hatNP[k][l] = dh;
		freqNP[k][l] = dq / (2.0*(double)NINDNP);
		if((tracelevel!=0)&&(k==0)&&(l==0)) fprintf(fptr,"k=%d l=%d posNP=%llu chromNP=%d sNP=%f hNP=%f freqNP=%f\n", k, l, posNP[k][l], chromNP[k][l], sNP[k][l], hsNP[k][l], freqNP[k][l]);
	}

	/* ***** take genotypic values of natural population ***** */

	fpop=fopen("dataBP.ped","r");

	for (i=0; i<NINDNP; i++)
	{
		lookfortext("IND");

		for (j=1; j<=5; j++)	fscanf(fpop,"%d", &x);

		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		{
			fscanf(fpop,"%d%d", &g0, &g1);

			if (g0 == 2)	gmNP[i][k][0]=(gmNP[i][k][0] | RM[l]);
			if (g1 == 2)	gmNP[i][k][1]=(gmNP[i][k][1] | RM[l]);
		}
	}

	fclose(fpop);

/*	if (tracelevel!=0)
	{
		fprintf(fptr, "\nNatural population genotypes\n");
		for (i=0; i<NIND; i++)
		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		if (i==0)
		{
			if ((gmNP[i][k][0] & RM[l])==RM[l])	fprintf(fptr, "1 ");
			else					fprintf(fptr, "0 ");
			if ((gmNP[i][k][1] & RM[l])==RM[l])	fprintf(fptr, "1 ");
			else					fprintf(fptr, "0 ");
		}
		fprintf(fptr, "\n");
	}
*/
	/* ***** estimate LEQ in the natural population ***** */

	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	{
		AA=0.0; Aa=0.0; aa=0.0;

		for (i=0; i<NINDNP; i++)
		{
			if (((gmNP[i][k][0] & RM[l])==RM[l])&&((gmNP[i][k][1] & RM[l])==RM[l]))		aa+=1.0;
	    		else if (((gmNP[i][k][0] & RM[l])!=RM[l])&&((gmNP[i][k][1] & RM[l])!=RM[l]))	AA+=1.0;
		     else	Aa+=1.0;
		}

		qNP[k][l] = (aa/(double)NINDNP)+(Aa/(2.0*(double)NINDNP));

//		if ((tracelevel!=0)&&(k==0))	fprintf(fptr, "k=%d l=%d AA=%f Aa=%f aa=%f q=%f\n", k, l, AA, Aa, aa, qNP[k][l]);

		LEQ_NP += 2.0 * ((sNP[k][l]*hsNP[k][l]) - (sNP[k][l]/2.0)) * qNP[k][l] * (1.0 - qNP[k][l]);
	}

/*	if (tracelevel!=0)
	{
		fprintf(fptr, "\n LEQ_NP = %f\n", LEQ_NP);
		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		if (k==0)
		fprintf(fptr, "\n k=%d l=%d   sNP=%f  hsNP=%f  qNP=%f", k, l, sNP[k][l], hsNP[k][l], qNP[k][l]);
	}
*/
	fclose(fdat);
}

/* ***************************************************** */

read_genes()
{
	if (nwindows == 23000000) fgenes=fopen("genes_2L.txt","r");
 	if (nwindows == 21100000) fgenes=fopen("genes_2R.txt","r");
	if (nwindows == 24500000) fgenes=fopen("genes_3L.txt","r");
	if (nwindows == 27900000) fgenes=fopen("genes_3R.txt","r");

 	while (!feof(fgenes))
 	{
	 	genes ++;
	 	fscanf(fgenes,"%llu", &x);
	 	startgene[genes] = x;
		fscanf(fgenes,"%llu", &x);
		endgene[genes] = x;
		fscanf(fgenes,"%llu", &x);
		lengthgene[genes] = x;
		//		if (tracelevel!=0) fprintf(fptr,"genes=%d\n", genes);
		//    if (tracelevel!=0) fprintf(fptr,"genes=%d, startgene=%llu, endgene=%llu, lengthgene=%llu\n", genes, startgene[genes], endgene[genes], lengthgene[genes]);
	}
	fclose(fgenes);

	for (g=1; g<genes; g++)
	{
 		for (j=0; j<(nwindows/window_size); j++)
 		{
			if ( (startgene[g] >= (j*window_size)) && (startgene[g] < ((j+1)*window_size)) )
			{
	 			accum (&numgenes[j], 1.0);
				accum (&length[j], lengthgene[g]);
	 			//if (tracelevel!=0)    fprintf(fptr,"j=%d numgenes=%f	lengthgene=%f\n", j, accsum(&numgenes[j]), accmean(&length[j]));
	 			break;
			}
 		}
	}

}

/* ***************************************************** */

read_Bs()
{
	double w;
    
	if (nwindows == 23000000) fB=fopen("B-2L-1kb.txt","r");
	if (nwindows == 21100000) fB=fopen("B-2R-1kb.txt","r");
	if (nwindows == 24500000) fB=fopen("B-3L-1kb.txt","r");
	if (nwindows == 27900000) fB=fopen("B-3R-1kb.txt","r");

	while (!feof(fB))
	{
		bb ++;
		fscanf(fB,"%llu", &x);
		posb[bb] = x;
		fscanf(fB,"%lf", &w);
		B[bb] = w;
		//if (tracelevel!=0) fprintf(fptr,"bb=%d\n", bb);
		//if (tracelevel!=0) fprintf(fptr,"bb=%d, posb=%llu, B=%f\n", bb, posb[bb], B[bb]);
	}
	fclose(fB);

	for (b=1; b<bb; b++)
	{
		for (j=0; j<(nwindows/100000); j++)
		{
			if ( (posb[b] >= (j*100000)) && (posb[b] < ((j+1)*100000)) )
			{
				accum (&Bswin[j], B[b]);
				//if (tracelevel!=0)    fprintf(fptr,"j=%d, b=%d,	B=%f\n", j, b, B[b]);
				//if (tracelevel!=0)    fprintf(fptr,"j=%d,	Bswin=%f\n", j, accmean(&Bswin[j]));
				break;
			}
		}
	}

}
/* ***************************************************** */

sample ()
{
	int g;

	/* ***** sample the first NIND individuals from the Base Population ***** */

	/* ** Randomise NINDNP individuals ** */

	for (i=0; i<NINDNP; i++)
	{
		ran_i = (int)(uniform() * NINDNP);

		for (k=0; k<NCRO; k++)
		{
			g=gmNP[i][k][0]; gmNP[i][k][0]=gmNP[ran_i][k][0]; gmNP[ran_i][k][0]=g;
			g=gmNP[i][k][1]; gmNP[i][k][1]=gmNP[ran_i][k][1]; gmNP[ran_i][k][1]=g;
		}
	}

/*	if (tracelevel!=0)
	{
		fprintf(fptr, "\nFirst individual natural population\n");
		for (i=0; i<NINDNP; i++)
		for (k=0; k<NCRO; k++)
			if ((i==0)&&(k>=0)&&(k<=10))	fprintf (fptr," i = %d  k = %d  gmNP0 = %d   gmNP1 = %d\n", i, k, gmNP[i][k][0], gmNP[i][k][1]);

		for (i=0; i<NINDNP; i++)
		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		if (i==21)
		{
			if ((gmNP[i][k][0] & RM[l])==RM[l])	fprintf(fptr, "1 ");
			else					fprintf(fptr, "0 ");
			if ((gmNP[i][k][1] & RM[l])==RM[l])	fprintf(fptr, "1 ");
			else					fprintf(fptr, "0 ");
		}
		fprintf(fptr, "\n");
	}
*/
	for (i=0; i<NIND; i++)
	for (k=0; k<NCRO; k++)
	{
		gm[i][k][0]=gmNP[i][k][0];
		gm[i][k][1]=gmNP[i][k][1];
	}

	if (tracelevel!=0)
	{
		fprintf(fptr, "\nIndividual sampled\n");
		for (i=0; i<NIND; i++)
		for (k=0; k<NCRO; k++)
			if ((i==0)&&(k>=0)&&(k<=10))	fprintf (fptr," i = %d  k = %d  gm0 = %d   gm1 = %d\n", i, k, gm[i][k][0], gm[i][k][1]);

		int ALL=0, homos=0;
		for (i=0; i<NIND; i++)
		{
			ALL=0; homos=0;

			for (k=0; k<NCRO; k++)
			for (l=0; l<NLOCI; l++)
			if (i==21)
			{
				if ((gm[i][k][0] & RM[l])==RM[l])	fprintf(fptr, "1 ");
				else					fprintf(fptr, "0 ");
				if ((gm[i][k][1] & RM[l])==RM[l])	fprintf(fptr, "1 ");
				else					fprintf(fptr, "0 ");

				ALL ++;
				if ( ( ((gm[i][k][0] & RM[l])==RM[l]) && ((gm[i][k][1] & RM[l])==RM[l]) ) || ( ((gm[i][k][0] & RM[l])!=RM[l]) && ((gm[i][k][1] & RM[l])!=RM[l]) ) )	homos ++;
			}
			if (i==21) fprintf(fptr, "\ni = %d   ALL = %d   homos = %d", i,  ALL, homos);
		}
	}

	/* ***** take effects of genes from Base Population ***** */

	if (rep == 1)
	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	{
		chrom[k][l] = chromNP[k][l];
		//fprintf(fptr,"chrom[%d][%d]=%d\n", k, l, chrom[k][l]);
		pos[k][l] = posNP[k][l];
		s[k][l] = sNP[k][l];
		at[k][l] = atNP[k][l];
		hs[k][l] = hsNP[k][l];
		hat[k][l] = hatNP[k][l];
//		if((k==4000)&&(l==0)) fprintf(fptr,"\n\nk=%d l=%d %llu chrom=%d\n\n", k, l, pos[k][l], chrom[k][l]);

//		if (tracelevel!=0)    if ((k<4)&&(l<5)) fprintf(fptr,"\n k=%d l=%d s=%f hs=%f", k, l, s[k][l], hs[k][l]);
	}
}

/* ***************************************************** */

frequency_genes ()
{
	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	{
		AA=0.0; Aa=0.0; aa=0.0;

		for (i=0; i<NIND; i++)
		{
			if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))		aa+=1.0;
	    		else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	AA+=1.0;
			else	Aa+=1.0;
		}

		q[k][l] = (aa/(double)NIND)+(Aa/(2.0*(double)NIND));
        
        //		if ((tracelevel!=0)&&(k==0))	fprintf(fptr, "k=%d l=%d s=%f  hs=%f  AA=%f Aa=%f aa=%f q=%f\n", k, l, s[k][l], hs[k][l], AA, Aa, aa, q[k][l]);

		if ((q[k][l] != 0.0)&&(q[k][l] != 1.0)&&(s[k][l]==0.0))	numSNPs ++;
        //		if ((q[k][l] != 0.0)&&(q[k][l] != 1.0))			numSNPs ++;
        //		if (tracelevel!=0)	if ((q[k][l] != 0.0)&&(q[k][l] != 1.0)) fprintf(fptr,"\nk=%d l=%d numSNPs=%d q=%f\n", k, l, numSNPs, q[k][l]);

		if (q[k][l] == 0.0)	accum (&q_0[0], 1.0);
		for (j=1; j<classes; j++)	if ((q[k][l] > ((double)(j-1)/classes))&&(q[k][l] <= ((double)j/classes))) accum (&q_0[j], 1.0);
	}

	fprintf(ffreq, "%f ", accsum(&q_0[0]));
	accum (&AVE_q_0[0], accsum(&q_0[0]));
	for (j=1; j<classes; j++)
	{
		fprintf(ffreq, "%f ", accsum(&q_0[j]));
		accum (&AVE_q_0[j], accsum(&q_0[j]));
	}

	/* ******************* ADDITIVE AND DOMINANCE VARIANCE AND INBREEDING DEPRESSION RATE ***************** */

	va_del = 0.0;
	vd_del = 0.0;
	id_del = 0.0;
    
	va_let = 0.0;
	vd_let = 0.0;
	id_let = 0.0;
    
	va_od = 0.0;
	vd_od = 0.0;
	id_od = 0.0;
    
	va_qn = 0.0;
	vd_qn = 0.0;
	id_qn = 0.0;
    
	va_adv = 0.0;
	vd_adv = 0.0;
	id_adv = 0.0;

	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	if ((q[k][l] != 0.0)&&(q[k][l] != 1.0)&&(s[k][l] != 0.0))
	{
		//lethals
		if (s[k][l]<= (-0.9))
		{
			d_a_let = (s[k][l]/2.0) * (2.0*hs[k][l] - 1.0);
			alfa_a_let = (-s[k][l]/2.0) + ( d_a_let * (2.0*q[k][l] - 1.0) );
			va_let += 2.0 * q[k][l] * (1.0-q[k][l]) * alfa_a_let * alfa_a_let;
			vd_let += (2.0 * d_a_let * q[k][l] * (1.0-q[k][l])) * (2.0 * d_a_let * q[k][l] * (1.0-q[k][l]));
			id_let += (2.0 * d_a_let * q[k][l] * (1.0-q[k][l]));
		}
        
		//overdominants
		else if ((s[k][l]>0.0) && (hs[k][l] == 1.5))
		{
 			s1[k][l]=(s[k][l]*hs[k][l])/(1+s[k][l]*hs[k][l]);
			 s2[k][l]=(s[k][l]*(1+hs[k][l]))/(1+s[k][l]*hs[k][l]);
            
			 d_a_od = (s[k][l]*(hs[k][l]-0.5))/(1+s[k][l]*hs[k][l]);
  			alfa_a_od = (q[k][l] * s2[k][l]) - ((1-q[k][l]) * s1[k][l]);
            
			va_od += 2.0 * q[k][l] * (1.0-q[k][l]) * alfa_a_od * alfa_a_od;
			vd_od += (2.0 * d_a_od * q[k][l] * (1.0-q[k][l])) * (2.0 * d_a_od * q[k][l] * (1.0-q[k][l]));
			id_od += (2.0 * d_a_od * q[k][l] * (1.0-q[k][l]));
		}
        
		//deleterious
		else if ((s[k][l]<0.0) && (s[k][l] > (-0.9)) && (s[k][l] != (-0.0001)))
		{
			d_a_del = (s[k][l]/2.0) * (2.0*hs[k][l] - 1.0);
			alfa_a_del = (-s[k][l]/2.0) + ( d_a_del * (2.0*q[k][l] - 1.0) );
 			va_del += 2.0 * q[k][l] * (1.0-q[k][l]) * alfa_a_del * alfa_a_del;
 			vd_del += (2.0 * d_a_del * q[k][l] * (1.0-q[k][l])) * (2.0 * d_a_del * q[k][l] * (1.0-q[k][l]));
			id_del += (2.0 * d_a_del * q[k][l] * (1.0-q[k][l]));
		}
        
        		//QN
        		else if (s[k][l] == (-0.0001))
		{
			d_a_qn = (s[k][l]/2.0) * (2.0*hs[k][l] - 1.0);
			alfa_a_qn = (-s[k][l]/2.0) + ( d_a_qn * (2.0*q[k][l] - 1.0) );
 			va_qn += 2.0 * q[k][l] * (1.0-q[k][l]) * alfa_a_qn * alfa_a_qn;
 			vd_qn += (2.0 * d_a_qn * q[k][l] * (1.0-q[k][l])) * (2.0 * d_a_qn * q[k][l] * (1.0-q[k][l]));
			id_qn += (2.0 * d_a_qn * q[k][l] * (1.0-q[k][l]));
		}
        
		//advantageous
		else if ((s[k][l]>0.0) && (hs[k][l] != 1.5))
		{
			d_a_adv = (s[k][l]/2.0) * (2.0*hs[k][l] - 1.0);
			alfa_a_adv = (-s[k][l]/2.0) + ( d_a_adv * (2.0*q[k][l] - 1.0) );
 			va_adv += 2.0 * q[k][l] * (1.0-q[k][l]) * alfa_a_adv * alfa_a_adv;
 			vd_adv += (2.0 * d_a_adv * q[k][l] * (1.0-q[k][l])) * (2.0 * d_a_adv * q[k][l] * (1.0-q[k][l]));
			id_adv += (2.0 * d_a_adv * q[k][l] * (1.0-q[k][l]));
 		}
	}

	accum (&AVE_VA_let, va_let);
	accum (&AVE_VD_let, vd_let);
	accum (&AVE_ID_let, id_let);

	accum (&AVE_VA_od, va_od);
	accum (&AVE_VD_od, vd_od);
	accum (&AVE_ID_od, id_od);

	accum (&AVE_VA_del, va_del);
	accum (&AVE_VD_del, vd_del);
	accum (&AVE_ID_del, id_del);
    
	accum (&AVE_VA_qn, va_qn);
	accum (&AVE_VD_qn, vd_qn);
	accum (&AVE_ID_qn, id_qn);
    
	accum (&AVE_VA_adv, va_adv);
	accum (&AVE_VD_adv, vd_adv);
	accum (&AVE_ID_adv, id_adv);

	accum (&AVE_VA, va_let + va_od + va_del + va_qn + va_adv);
	accum (&AVE_VD, vd_let + vd_od + vd_del + vd_qn + vd_adv);
	accum (&AVE_ID, id_let + id_od + id_del + id_qn + id_adv);
    
	/* ******************* windows ADDITIVE AND DOMINANCE VARIANCE AND INBREEDING DEPRESSION RATE ***************** */

	for (j=0; j<(nwindows/window_size); j++)
	{
		va_letw[j] = 0.0; vd_letw[j] = 0.0; id_letw[j] = 0.0;
		va_odw[j] = 0.0; vd_odw[j] = 0.0; id_odw[j] = 0.0;
		va_delw[j] = 0.0; vd_delw[j] = 0.0; id_delw[j] = 0.0;
		va_qnw[j] = 0.0; vd_qnw[j] = 0.0; id_qnw[j] = 0.0;
		va_advw[j] = 0.0; vd_advw[j] = 0.0; id_advw[j] = 0.0;
	}

	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	if ((q[k][l] != 0.0)&&(q[k][l] != 1.0)&&(s[k][l] != 0.0))
    	{

        for (j=0; j<(nwindows/window_size); j++)
        if ( (pos[k][l] >= (j*window_size)) && (pos[k][l] < ((j+1)*window_size)) )
        {
            //lethals
            if (s[k][l]<= (-0.9))
            {
                d_a_letw = (s[k][l]/2.0) * (2.0*hs[k][l] - 1.0);
                alfa_a_letw = (-s[k][l]/2.0) + ( d_a_letw * (2.0*q[k][l] - 1.0) );
                va_letw[j] += 2.0 * q[k][l] * (1.0-q[k][l]) * alfa_a_letw * alfa_a_letw;
                vd_letw[j] += (2.0 * d_a_letw * q[k][l] * (1.0-q[k][l])) * (2.0 * d_a_letw * q[k][l] * (1.0-q[k][l]));
                id_letw[j] += (2.0 * d_a_letw * q[k][l] * (1.0-q[k][l]));
            }
            
            //overdominants
            else if ((s[k][l]>0.0) && (hs[k][l] == 1.5))
            {  
                s1w[k][l]=(s[k][l]*hs[k][l])/(1+s[k][l]*hs[k][l]);
                s2w[k][l]=(s[k][l]*(1+hs[k][l]))/(1+s[k][l]*hs[k][l]);
            
                d_a_odw = (s[k][l]*(hs[k][l]-0.5))/(1+s[k][l]*hs[k][l]);
                alfa_a_odw = (q[k][l] * s2w[k][l]) - ((1-q[k][l]) * s1w[k][l]);
                va_odw[j] += 2.0 * q[k][l] * (1.0-q[k][l]) * alfa_a_odw * alfa_a_odw;
                vd_odw[j] += (2.0 * d_a_odw * q[k][l] * (1.0-q[k][l])) * (2.0 * d_a_odw * q[k][l] * (1.0-q[k][l]));
                id_odw[j] += (2.0 * d_a_odw * q[k][l] * (1.0-q[k][l]));
                                
                //if ((tracelevel!=0) && (j==10)) fprintf(fptr,"rep=%d j=%d va_odw=%10.8f va_odwin=%10.8f\n", rep, j, va_odw, accsum(&va_odwin[j]));
                
                //if ((tracelevel!=0) && (j==10)) fprintf(fptr,"rep=%d j=%d vd_odw=%10.8f vd_odwin=%10.8f\n", rep, j, vd_odw, accsum(&vd_odwin[j]));
                
                //if ((tracelevel!=0) && (j==10)) fprintf(fptr,"rep=%d j=%d id_odw=%10.8f id_odwin=%10.8f\n", rep, j, id_odw, accsum(&id_odwin[j]));
            } 
            
            //deleterious
            else if ((s[k][l]<0.0) && (s[k][l] > (-0.9)) && (s[k][l] != (-0.0001)))
            {
                d_a_delw = (s[k][l]/2.0) * (2.0*hs[k][l] - 1.0);
                alfa_a_delw = (-s[k][l]/2.0) + ( d_a_delw * (2.0*q[k][l] - 1.0) );
                va_delw[j] += 2.0 * q[k][l] * (1.0-q[k][l]) * alfa_a_delw * alfa_a_delw;
                vd_delw[j] += (2.0 * d_a_delw * q[k][l] * (1.0-q[k][l])) * (2.0 * d_a_delw * q[k][l] * (1.0-q[k][l]));
                id_delw[j] += (2.0 * d_a_delw * q[k][l] * (1.0-q[k][l]));
                
                //if ((tracelevel!=0) && (j==1)) fprintf(fptr,"rep=%d j=%d id_delw=%15.12f\n", rep, j, id_delw[j]);
                //if ((tracelevel!=0) && (j==1)) fprintf(fptr,"rep=%d j=%d va_delw=%15.12f\n", rep, j, va_delw[j]);
                //if ((tracelevel!=0) && (j==1)) fprintf(fptr,"rep=%d j=%d vd_delw=%15.12f\n", rep, j, vd_delw[j]);
            } 
            
            //QN
            else if (s[k][l] == (-0.0001))
            {
                d_a_qnw = (s[k][l]/2.0) * (2.0*hs[k][l] - 1.0);
                alfa_a_qnw = (-s[k][l]/2.0) + ( d_a_qnw * (2.0*q[k][l] - 1.0) );
                va_qnw[j] += 2.0 * q[k][l] * (1.0-q[k][l]) * alfa_a_qnw * alfa_a_qnw;
                vd_qnw[j] += (2.0 * d_a_qnw * q[k][l] * (1.0-q[k][l])) * (2.0 * d_a_qnw * q[k][l] * (1.0-q[k][l]));
                id_qnw[j] += (2.0 * d_a_qnw * q[k][l] * (1.0-q[k][l]));
            }
            
            //advantageous
            else if ((s[k][l]>0.0) && (hs[k][l] != 1.5))
            {
                d_a_advw = (s[k][l]/2.0) * (2.0*hs[k][l] - 1.0);
                alfa_a_advw = (-s[k][l]/2.0) + ( d_a_advw * (2.0*q[k][l] - 1.0) );
                va_advw[j] += 2.0 * q[k][l] * (1.0-q[k][l]) * alfa_a_advw * alfa_a_advw;
                vd_advw[j] += (2.0 * d_a_advw * q[k][l] * (1.0-q[k][l])) * (2.0 * d_a_advw * q[k][l] * (1.0-q[k][l]));
                id_advw[j] += (2.0 * d_a_advw * q[k][l] * (1.0-q[k][l]));
                
                //if ((tracelevel!=0) && (j==1)) fprintf(fptr,"rep=%d j=%d id_advw=%15.12f\n", rep, j, id_advw[j]);
                
                //if (tracelevel!=0)  fprintf(fptr,"rep=%d j=%d k=%d l=%d s=%f  h=%f d_a_advw=%f alfa_a_advw=%f va_advw=%15.12f\n", rep, j, k, l, s[k][l], hs[k][l], d_a_advw, alfa_a_advw, va_advw[j]);
                
                //if ((tracelevel!=0) && (j==1)) fprintf(fptr,"rep=%d j=%d vd_advw=%15.12f\n", rep, j, vd_advw[j]);
            }
        }
    }

 	for (j=0; j<(nwindows/window_size); j++)
	{
		accum (&AVE_va_letwin[j], va_letw[j]);
		accum (&AVE_vd_letwin[j], vd_letw[j]);
		accum (&AVE_id_letwin[j], id_letw[j]);
		accum (&AVE_va_odwin[j], va_odw[j]);
		accum (&AVE_vd_odwin[j], vd_odw[j]);
		accum (&AVE_id_odwin[j], id_odw[j]);
		accum (&AVE_va_delwin[j], va_delw[j]);
		accum (&AVE_vd_delwin[j], vd_delw[j]);
		accum (&AVE_id_delwin[j], id_delw[j]);
		accum (&AVE_va_qnwin[j], va_qnw[j]);
		accum (&AVE_vd_qnwin[j], vd_qnw[j]);
		accum (&AVE_id_qnwin[j], id_qnw[j]);
		accum (&AVE_va_advwin[j], va_advw[j]);
		accum (&AVE_vd_advwin[j], vd_advw[j]);
		accum (&AVE_id_advwin[j], id_advw[j]);
          
 //       if (tracelevel!=0)  fprintf(fptr,"rep=%d j=%d AVE_va_advw=%15.12f\n", rep, j, accmean(&AVE_va_advwin[j])); 
//        if (tracelevel!=0)  fprintf(fptr,"rep=%d j=%d AVE_vd_advw=%15.12f\n", rep, j, accmean(&AVE_vd_advwin[j]));
 //       if (tracelevel!=0)  fprintf(fptr,"rep=%d j=%d AVE_id_advw=%15.12f\n", rep, j, accmean(&AVE_id_advwin[j]));
        
        //if (tracelevel!=0)  fprintf(fptr,"rep=%d j=%d AVE_va_qnw=%15.12f\n", rep, j, accmean(&AVE_va_qnwin[j])); 
        //if (tracelevel!=0)  fprintf(fptr,"rep=%d j=%d va_qnw=%15.12f\n", rep, j, va_qnw[j]);  
	}
    
    /* ******************* windows s, h, q ***************** */

    numSNPslet = 0.0;
    numSNPsod = 0.0;
    numSNPsdel = 0.0;
    numSNPsneu = 0.0;
    numSNPsqn = 0.0;
    numSNPsadv = 0.0;
    
    for (k=0; k<NCRO; k++)
    for (l=0; l<NLOCI; l++)
    if ((q[k][l] != 0.0)&&(q[k][l] != 1.0)&&(s[k][l] != 0.0))
    {
        //lethals
        if (s[k][l]<= (-0.9))
        {
            numSNPslet++;
            //if (tracelevel!=0) fprintf(fptr,"rep=%d k=%d l=%d pos=%llu q=%f h=%f s=%f let=%d\n", rep, k, l, pos[k][l], q[k][l], hs[k][l], s[k][l], numSNPslet);
        }
        
        //overdominants
        else if ((s[k][l]>0.0) && (hs[k][l] == 1.5))
        {
            numSNPsod ++;
            //if (tracelevel!=0) fprintf(fptr,"rep=%d k=%d l=%d pos=%llu q=%f h=%f s=%f od=%d\n", rep, k, l, pos[k][l], q[k][l], hs[k][l], s[k][l], numSNPsod);
        }
        
        //deleterious
        else if ((s[k][l]<0.0) && (s[k][l] > (-0.9)) && (s[k][l] != (-0.0001)))
        {
            numSNPsdel ++;
            //if (tracelevel!=0) fprintf(fptr,"rep=%d k=%d l=%d pos=%llu q=%f h=%f s=%f del=%d\n", rep, k, l, pos[k][l], q[k][l], hs[k][l], s[k][l], numSNPsdel);
        } 
        
        //QN
        else if ((s[k][l] == (-0.0001)))
        {
            numSNPsqn ++;
            //if (tracelevel!=0) fprintf(fptr,"rep=%d k=%d l=%d pos=%llu q=%f h=%f s=%f qn=%d\n", rep, k, l, pos[k][l], q[k][l], hs[k][l], s[k][l], numSNPsqn);
        }
        
        //advantageous
         else if ((s[k][l]>0.0) && (hs[k][l] != 1.5))
        {
            numSNPsadv ++;
        }
    }
 
	//neutrals
	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	if ((q[k][l] != 0.0)&&(q[k][l] != 1.0)&&(s[k][l]==0.0))
    	{
		numSNPsneu ++;
	}
    
    accum (&AVE_numSNPslet, numSNPslet);
    accum (&AVE_numSNPsod, numSNPsod);
    accum (&AVE_numSNPsdel, numSNPsdel);
    accum (&AVE_numSNPsqn, numSNPsqn);
    accum (&AVE_numSNPsadv, numSNPsadv);
    accum (&AVE_numSNPsneu, numSNPsneu);

    
    for (k=0; k<NCRO; k++)
    for (l=0; l<NLOCI; l++)
    if ((q[k][l] != 0.0)&&(q[k][l] != 1.0)&&(s[k][l] != 0.0))
    {
        for (j=0; j<(nwindows/window_size); j++)
        if ( (pos[k][l] >= (j*window_size)) && (pos[k][l] < ((j+1)*window_size)) )
        {
            //lethals
            if (s[k][l]<= (-0.9))
            {
                accum (&letwin[j], 1.0); 
                accum (&hsletwin[j], hs[k][l]);
                accum (&sletwin[j], s[k][l]);
                accum (&qletwin[j], q[k][l]);
               //if ((tracelevel!=0) && (j==5))    fprintf(fptr,"rep=%d k=%d l=%d pos=%llu j=%d let=%f\n", rep, k, l, pos[k][l], j, accsum(&letwin[j]));
            }
            
            //overdominants
            else if ((s[k][l]>0.0) && (hs[k][l] == 1.5))
            {
                accum (&odwin[j], 1.0);
                accum (&hsodwin[j], hs[k][l]);
                accum (&sodwin[j], s[k][l]);
                accum (&qodwin[j], q[k][l]);
                //if (tracelevel!=0) fprintf(fptr,"rep=%d k=%d l=%d pos=%llu j=%d od=%f q=%f h=%f s=%f\n", rep, k, l, pos[k][l], j, accsum(&odwin[j]), q[k][l], hs[k][l], s[k][l]);
            }
            
            //deleterious
            else if ((s[k][l]<0.0) && (s[k][l] > (-0.9)) && (s[k][l] != (-0.0001)))
            {
                accum (&delwin[j], 1.0);
                accum (&hsdelwin[j], hs[k][l]);
                accum (&sdelwin[j], s[k][l]);
                accum (&qdelwin[j], q[k][l]);
                
                //if (tracelevel!=0) fprintf(fptr,"rep=%d k=%d l=%d pos=%llu j=%d del=%f q=%f h=%f s=%f\n", rep, k, l, pos[k][l], j, accsum(&delwin[j]), q[k][l], hs[k][l], s[k][l]);
            }
            
            //QN
            else if (s[k][l] == (-0.0001))
            {
                accum (&qnwin[j], 1.0);
                accum (&hsqnwin[j], hs[k][l]);
                accum (&sqnwin[j], s[k][l]);
                accum (&qqnwin[j], q[k][l]);
                //if ((tracelevel!=0) && (j==5))    fprintf(fptr,"rep=%d k=%d l=%d pos=%llu j=%d qn=%f q=%f h=%f s=%f\n", rep, k, l, pos[k][l], j, accsum(&qnwin[j]), q[k][l], hs[k][l], s[k][l]);
            }
            
            //advantageous
            else if ((s[k][l]>0.0) && (hs[k][l] != 1.5))
            {
                accum (&advwin[j], 1.0);
                accum (&hsadvwin[j], hs[k][l]);
                accum (&sadvwin[j], s[k][l]);
                accum (&qadvwin[j], q[k][l]);
            }
        }
    }
    
    for (k=0; k<NCRO; k++)
    for (l=0; l<NLOCI; l++)
    if ((q[k][l] != 0.0)&&(q[k][l] != 1.0)&&(s[k][l] == 0.0))
    //neutrals
    {
        for (j=0; j<(nwindows/window_size); j++)
        if ( (pos[k][l] >= (j*window_size)) && (pos[k][l] < ((j+1)*window_size)) )
        {
            accum (&neuwin[j], 1.0);
            accum (&qneuwin[j], q[k][l]);
        }
    }

    for (j=0; j<(nwindows/window_size); j++)
    {
        accum (&AVE_neuwin[j], accsum(&neuwin[j]));
        accum (&AVE_letwin[j], accsum(&letwin[j]));
        accum (&AVE_odwin[j], accsum(&odwin[j]));
        accum (&AVE_delwin[j], accsum(&delwin[j]));
        accum (&AVE_qnwin[j], accsum(&qnwin[j]));
        accum (&AVE_advwin[j], accsum(&advwin[j]));
        
        //neutrals
        if (accmean(&qneuwin[j]) >= 0.0) accum (&AVE_qneuwin[j], accmean(&qneuwin[j]));
        
        //deleterious
        if (accmean(&hsdelwin[j]) >= 0.0) accum (&AVE_hsdelwin[j], accmean(&hsdelwin[j]));
        if (accmean(&sdelwin[j]) <= 0.0) accum (&AVE_sdelwin[j], accmean(&sdelwin[j]));
        if (accmean(&qdelwin[j]) >= 0.0) accum (&AVE_qdelwin[j], accmean(&qdelwin[j]));
        
        ///if (tracelevel!=0)   fprintf(fptr,"rep=%d j=%d del=%f q=%f h=%f s=%f let=%f  od=%f\n", rep, j, accsum(&del[j]), accmean(&qdelwin[j]), accmean(&hsdelwin[j]), accmean(&sdelwin[j]), accsum(&letwin[j]), accsum(&odwin[j]));

        //overdominants
        if (accmean(&hsodwin[j]) >= 0.0) accum (&AVE_hsodwin[j], accmean(&hsodwin[j]));
        if (accmean(&sodwin[j]) >= 0.0) accum (&AVE_sodwin[j], accmean(&sodwin[j]));
        if (accmean(&qodwin[j]) >= 0.0) accum (&AVE_qodwin[j], accmean(&qodwin[j]));
        
        //lethals
        if (accmean(&hsletwin[j]) >= 0.0) accum (&AVE_hsletwin[j], accmean(&hsletwin[j]));
        if (accmean(&sletwin[j]) <= 0.0) accum (&AVE_sletwin[j], accmean(&sletwin[j]));
        if (accmean(&qletwin[j]) >= 0.0) accum (&AVE_qletwin[j], accmean(&qletwin[j]));
        
        //QN
        if (accmean(&hsqnwin[j]) >= 0.0) accum (&AVE_hsqnwin[j], accmean(&hsqnwin[j]));
        if (accmean(&sqnwin[j]) <= 0.0) accum (&AVE_sqnwin[j], accmean(&sqnwin[j]));
        if (accmean(&qqnwin[j]) >= 0.0) accum (&AVE_qqnwin[j], accmean(&qqnwin[j]));
        
        //advantageous
        if (accmean(&hsadvwin[j]) >= 0.0) accum (&AVE_hsadvwin[j], accmean(&hsadvwin[j]));
        if (accmean(&sadvwin[j]) >= 0.0) accum (&AVE_sadvwin[j], accmean(&sadvwin[j]));
        if (accmean(&qadvwin[j]) >= 0.0) accum (&AVE_qadvwin[j], accmean(&qadvwin[j]));
    }
}
/* ***************************************************** */

genotypic_values ()
{
//	if (tracelevel!=0)	fprintf(fptr,"\n\n Genotypic values \n");

	for (i=0; i<NIND; i++)
	{
		pm_s[i] = 1.0;
		pm_f[i] = 1.0;

		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		//overdominance
		if ((s[k][l]>0.0) && (hs[k][l] == 1.5))
			{
				if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))  	pm_s[i] *= 1.0 - s2[k][l] /* aa */;
				else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l])) pm_s[i] *= 1.0 - s1[k][l] /* AA */;
				else	pm_s[i] *= 1.0 /* Aa */;
			}
        
        else
			{
				if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))  	pm_s[i] *= (1.0 + s[k][l]) /* aa */;
				else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l])) /* AA */;
				else	pm_s[i] *= (1.0 + (s[k][l]*hs[k][l])) /* Aa */;
			}
		//if (tracelevel!=0)
		//{
			//if (i <= 10)
			//fprintf(fptr,"genotypic_values %d    pm_s = %f\n", i, pm_s[i]);
			//for (k=0; k<NCRO; k++) if ((i==0)&&(k<=2)&&(gm[i][k][0]!=0)) fprintf (fptr,"%d   gm0=%d   gm1=%d\n", i, gm[i][k][0], gm[i][k][1]);
		//}
		}

}

/* ***************************************************** */

phenotypeB ()
{
	int ii, it;
	double gsum_s=0.0, gsum2_s=0.0, gsum_f=0.0;

	for (i=0; i<NIND; i++)
	{
		gsum_s += pm_s[i];
		gsum2_s += (pm_s[i]*pm_s[i]);
	}

	accum (&gmean_s, gsum_s/(double)NIND);
	accum (&gvar_s, (gsum2_s - (gsum_s*gsum_s / (double)NIND)) / ((double)NIND - 1.0));

//	if (tracelevel!=0)   fprintf(fptr,"\ngmean_s = %f  gvar_s = %f C2 = %f\n",
//	gsum_s/(double)NIND, (gsum2_s - (gsum_s*gsum_s / (double)NIND)) / (double)NIND, ( (gsum2_s*(double)NIND) / (gsum_s*gsum_s) ) - 1.0);

}

/* ***************************************************** */

dumpphenotypes()
{
//	if (tracelevel==0)   return (0);

//	fprintf(fptr,"\n Fitness values\n");
//	for (i=0; i<NIND; i++)   fprintf(fptr,"i=%d pm_s=%f\n", i, pm_s[i]);
}

/* ***************************************************** */

plink_files()
{
	int last, lastpos;

	ffped = fopen ("data.ped","w");
	ffmap = fopen ("data.map","w");

	// data.map

	last = 1;
	lastpos = 1;
	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	if ((q[k][l] != 0.0)&&(q[k][l] != 1.0)&&(s[k][l]==0.0))
	{
		if (chrom[k][l] != last)
		{
			last = chrom[k][l];
			lastpos = pos[k][l];
		}
		fprintf(ffmap,"%d SNP%llu 0 %llu\n", chrom[k][l], pos[k][l], pos[k][l]-lastpos);
	}

	for (i=0; i<NIND; i++)
	{
		// data.ped

		fprintf(ffped,"1 IND%d 0 0 1 -9 ", i);
		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		if ((q[k][l] != 0.0)&&(q[k][l] != 1.0)&&(s[k][l]==0.0))
		{
			if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))		fprintf(ffped,"T T ");
	    		else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	fprintf(ffped,"A A ");
	    		else if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))	fprintf(ffped,"T A ");
	    		else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))	fprintf(ffped,"A T ");
		}
		fprintf(ffped,"\n");
	}

	fclose (ffped);
	fclose (ffmap);

	int ALL = 0, homos = 0;
	if (tracelevel!=0)
	{
		fprintf(fptr, "\nOnly neutral SNPS\n");
		for (i=0; i<NIND; i++)
		{
			ALL = 0; homos = 0;
			for (k=0; k<NCRO; k++)
			for (l=0; l<NLOCI; l++)
			if ((q[k][l] != 0.0)&&(q[k][l] != 1.0)&&(s[k][l]==0.0))
//			if (i==40)
			{
//				if ((gm[i][k][0] & RM[l])==RM[l])	fprintf(fptr, "1 ");
//				else					fprintf(fptr, "0 ");
//				if ((gm[i][k][1] & RM[l])==RM[l])	fprintf(fptr, "1 ");
//				else					fprintf(fptr, "0 ");

				ALL ++;
				if ( ( ((gm[i][k][0] & RM[l])==RM[l]) && ((gm[i][k][1] & RM[l])==RM[l]) ) || ( ((gm[i][k][0] & RM[l])!=RM[l]) && ((gm[i][k][1] & RM[l])!=RM[l]) ) )	homos ++;
			}
			fprintf(fptr, "\ni = %d   ALL = %d   homos = %d", i,  ALL, homos);
		}
	}
}

/* ***************************************************** */

READ_FILES()
{
	double w;

	fr2file=fopen("r2file","r");
	fpifile=fopen("pifile","r");
	fTajDfile=fopen("TajimaDfile","r");
    
    if (nwindows == 23000000) fRec=fopen("recombinationfile_2L","r");
    if (nwindows == 21100000) fRec=fopen("recombinationfile_2R","r");
    if (nwindows == 24500000) fRec=fopen("recombinationfile_3L","r");
    if (nwindows == 27900000) fRec=fopen("recombinationfile_3R","r");

//	if (tracelevel!=0)   fprintf(fptr,"numSNPs=%d\n", numSNPs);

	for (ss=1; ss<numSNPs; ss++)
	{
		fscanf(fr2file,"%llu", &x);
		posplink[ss] = x;
		fscanf(fr2file,"%lf", &w);
		r2[ss] = w;
		fscanf(fpifile,"%lf", &w);
		pii[ss] = w;

    for (j=0; j<(nwindows/window_size); j++)
		{
			if ( (posplink[ss] >= (j*window_size)) && (posplink[ss] < ((j+1)*window_size)) )
			{
				accum (&numr2win[j], 1.0);
				if (accmean(&numr2win[j]) > 0.0) accum (&r2win[j], r2[ss]);
				accum (&numpiwin[j], 1.0);
				if (accmean(&numpiwin[j]) > 0.0) accum (&piwin[j], pii[ss]);
                //if ((tracelevel!=0) && (j==5)) fprintf(fptr,"rep=%d snp=%d posplink=%d r2=%f pi=%f\n", rep, ss, posplink[ss], r2[ss], pii[ss]);
				break;
			 }
		 }
 	  }

	for (j=0; j<(nwindows/window_size); j++)
	{
			fscanf(fTajDfile,"%lf", &w);
			numTajDwin[j] = w;
			fscanf(fTajDfile,"%lf", &w);
			TajDwin[j] = w;
			//if ((tracelevel!=0) && (j==5)) fprintf(fptr,"rep=%d j=%d TajD=%f\n", rep, j, TajDwin[j]);

	}

	for (j=0; j<(nwindows/window_size); j++)
	{
			fscanf(fRec,"%lf", &w);
			cwin[j] = w;
	}

	for (j=0; j<(nwindows/window_size); j++)
	{
		accum (&AVE_numr2win[j], accsum(&numr2win[j]));
		accum (&AVE_r2win[j], accmean(&r2win[j]));
		accum (&AVE_numpiwin[j], accsum(&numpiwin[j]));
		accum (&AVE_piwin[j], accmean(&piwin[j]));
		accum (&AVE_numTajDwin[j], numTajDwin[j]);
		accum (&AVE_TajDwin[j], TajDwin[j]);
		accum (&AVE_cwin[j], cwin[j]);
		//if (tracelevel!=0)    fprintf(fptr,"rep=%d j=%d numr2=%f r2=%f numpi=%f pii=%f numTajD=%f TajD=%f\n", rep, j, accsum(&numr2win[j]), accmean(&r2win[j]), accsum(&numpiwin[j]), accmean(&piwin[j]), numTajDwin[j], TajDwin[j]);
	}

}

/* ***************************************************** */

ID_HOMOZYGOTES()
{
	int rnd;

	for (i=0; i<NIND; i++)
	{
		FITN[i] = 1.0;
		FITN_F[i] = 1.0;

		if (uniform() < 0.5) rnd = 1;
		else				 rnd = 0;

		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)
		{
	    		if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))
			{
				FITN[i] *= (1.0 + s[k][l]*hs[k][l]);

				if (rnd == 1)
				{
					FITN_F[i] *= (1.0 + s[k][l]);
				}
				else		/*11*/;
			}
	    		else if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))
			{
				FITN[i] *= (1.0 + s[k][l]*hs[k][l]);

				if (rnd == 1)	/*11*/;
				else
				{
					FITN_F[i] *= (1.0 + s[k][l]);
				}
			}
			else if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))
			{
				FITN[i] *= (1.0 + s[k][l]);
				FITN_F[i] *= (1.0 + s[k][l]);
			}
		}
	}

	mean_FITN = 0.0;
	mean_FITN_F = 0.0;

	for (i=0; i<NIND; i++)
	{
		mean_FITN += FITN[i]/(double)NIND;
		mean_FITN_F += FITN_F[i]/(double)NIND;
	}

	accum(&AVE_mean_FITN, mean_FITN);
	accum(&AVE_mean_FITN_F, mean_FITN_F);

	accum(&AVE_ID_FITN, -log(mean_FITN_F/mean_FITN));

	return(0);
}

/* ***************************************************** */

settozero()
{
	for (j=0; j<classes; j++)	initacc (&q_0[j]);
	for (j=0; j<(nwindows/window_size); j++)
	{
		initacc (&numr2win[j]);
		initacc (&r2win[j]);
        
		initacc (&numpiwin[j]);
		initacc (&piwin[j]);
        
		initacc (&numTajDwin[j]);
		initacc (&TajDwin[j]);
        
		initacc (&cwin[j]);
        
        initacc (&neuwin[j]);
        initacc (&qneuwin[j]);
        
		initacc (&delwin[j]);
		initacc (&hsdelwin[j]);
		initacc (&sdelwin[j]);
		initacc (&qdelwin[j]);
        
        initacc (&qnwin[j]);
		initacc (&hsqnwin[j]);
		initacc (&sqnwin[j]);
		initacc (&qqnwin[j]);
        
        initacc (&advwin[j]);
		initacc (&hsadvwin[j]);
		initacc (&sadvwin[j]);
		initacc (&qadvwin[j]);
        
        initacc (&odwin[j]);
        initacc (&hsodwin[j]);
		initacc (&sodwin[j]);
        initacc (&qodwin[j]);
        
        initacc (&letwin[j]);
        initacc (&hsletwin[j]);
		initacc (&sletwin[j]);
        initacc (&qletwin[j]);
	}
}

/* ***************************************************** */

printout()
{
	double TOT_VA[NW], TOT_VD[NW], TOT_ID[NW];
        
	fgen = fopen ("genfile.dat","w");

	fprintf(fgen, "NINDNP=%d  NIND=%d  NCRO=%d  TOTLOCI=%d  numSNPNP=%d  numSNPsneu=%f  numSNPsdel=%f  numSNPsqn=%f  numSNPslet=%f  numSNPod=%f  numSNPadv=%f  LEQ_NP=%f\n", NINDNP, NIND, NCRO, TOTLOCI, numSNPNP, accmean(&AVE_numSNPsneu), accmean(&AVE_numSNPsdel), accmean(&AVE_numSNPsqn), accmean(&AVE_numSNPslet), accmean(&AVE_numSNPsod), accmean(&AVE_numSNPsadv), LEQ_NP);

	fprintf(fgen, "\n*********** Mean fitness and phenotype ***********\n\n");

	fprintf(fgen, "W           VG          VA          VAdel       VAqn       VAlet       VAod        VAadv        VD          VDdel          VDqn       VDlet       VDod        VDadv        B          Bdel          Bqn        Blet        Bod        Badv\n");
		fprintf(fgen, "%10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f  %10.8f\n", accmean(&gmean_s), accmean(&gvar_s), accmean(&AVE_VA), accmean(&AVE_VA_del), accmean(&AVE_VA_qn), accmean(&AVE_VA_let), accmean(&AVE_VA_od), accmean(&AVE_VA_adv), accmean(&AVE_VD), accmean(&AVE_VD_del), accmean(&AVE_VD_qn),accmean(&AVE_VD_let), accmean(&AVE_VD_od), accmean(&AVE_VD_adv), accmean(&AVE_ID), accmean(&AVE_ID_del), accmean(&AVE_ID_qn), accmean(&AVE_ID_let), accmean(&AVE_ID_od), accmean(&AVE_ID_adv));
    
	/*fprintf(fgen, "SD_W      SD_VG     SD_VA     SD_VAdel  SD_VAqn  SD_VAlet  SD_VAod   SD_VD     SD_VDdel     SD_VDqn  SD_VDlet  SD_VDod   SD_B      SD_Bdel      SD_Bqn   SD_Blet   SD_Bod\n");
		fprintf(fgen,   "%f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f  %f\n", sqrt(variance(&gmean_s)), sqrt(variance(&gvar_s)), sqrt(variance(&AVE_VA)), sqrt(variance(&AVE_VA_del)), sqrt(variance(&AVE_VA_qn)), sqrt(variance(&AVE_VA_let)), sqrt(variance(&AVE_VA_od)), sqrt(variance(&AVE_VD)), sqrt(variance(&AVE_VD_del)), sqrt(variance(&AVE_VD_qn)), sqrt(variance(&AVE_VD_let)), sqrt(variance(&AVE_VD_od)), sqrt(variance(&AVE_ID)), sqrt(variance(&AVE_ID_del)), sqrt(variance(&AVE_ID_qn)), sqrt(variance(&AVE_ID_let)), sqrt(variance(&AVE_ID_od)));
*/
	fprintf(fgen, "\n*********** ID from homozygotes ***********\n\n");

	fprintf(fgen, "mean_FITN=%6.4f    mean_FITN_F=%6.4f    ID_F1=%6.4f\n", accmean(&AVE_mean_FITN), accmean(&AVE_mean_FITN_F), accmean(&AVE_ID_FITN));
  	fprintf(fgen, "SD_FITN=%6.4f    SD_FITN_F=%6.4f    SD_ID_F1=%6.4f\n", sqrt(variance(&AVE_mean_FITN)), sqrt(variance(&AVE_mean_FITN_F)), sqrt(variance(&AVE_ID_FITN)));

	fprintf(fgen, "\n*********** Frequencies ***********\n\n");

	fprintf(fgen, "0.00 ");
	for (j=1; j<classes; j++)	fprintf(fgen, "%4.2f ", (double)j/100);
	fprintf(fgen, "\n");
	fprintf(fgen, "%f ", accmean(&AVE_q_0[0]));
	for (j=1; j<classes; j++)	fprintf(fgen, "%f ", accmean(&AVE_q_0[j]));

	for (j=0; j<(nwindows/window_size); j++) 
	{
		TOT_VA[j] = 0.0;
		if (accmean(&AVE_va_delwin[j]) > 0.0) TOT_VA[j] += accmean(&AVE_va_delwin[j]); 
		if (accmean(&AVE_va_letwin[j]) > 0.0) TOT_VA[j] += accmean(&AVE_va_letwin[j]); 
		if (accmean(&AVE_va_odwin[j]) > 0.0) TOT_VA[j] += accmean(&AVE_va_odwin[j]); 
		if (accmean(&AVE_va_qnwin[j]) > 0.0) TOT_VA[j] += accmean(&AVE_va_qnwin[j]); 
		if (accmean(&AVE_va_advwin[j]) > 0.0) TOT_VA[j] += accmean(&AVE_va_advwin[j]); 
    
		TOT_VD[j] = 0.0;
		if (accmean(&AVE_vd_delwin[j]) > 0.0) TOT_VD[j] += accmean(&AVE_vd_delwin[j]); 
		if (accmean(&AVE_vd_letwin[j]) > 0.0) TOT_VD[j] += accmean(&AVE_vd_letwin[j]); 
		if (accmean(&AVE_vd_odwin[j]) > 0.0) TOT_VD[j] += accmean(&AVE_vd_odwin[j]);
		if (accmean(&AVE_vd_qnwin[j]) > 0.0) TOT_VD[j] += accmean(&AVE_vd_qnwin[j]); 
		if (accmean(&AVE_vd_advwin[j]) > 0.0) TOT_VD[j] += accmean(&AVE_vd_advwin[j]); 
    
		TOT_ID[j] = 0.0;
		if (accmean(&AVE_id_delwin[j]) > 0.0) TOT_ID[j] += accmean(&AVE_id_delwin[j]); 
		if (accmean(&AVE_id_letwin[j]) > 0.0) TOT_ID[j] += accmean(&AVE_id_letwin[j]); 
		if (accmean(&AVE_id_odwin[j]) > 0.0) TOT_ID[j] += accmean(&AVE_id_odwin[j]);
		if (accmean(&AVE_id_qnwin[j]) > 0.0) TOT_ID[j] += accmean(&AVE_id_qnwin[j]); 
		if (accmean(&AVE_id_advwin[j]) > 0.0) TOT_ID[j] += accmean(&AVE_id_advwin[j]); 
    }
    
	fprintf(fgen, "\n\n*********** windows 100kb ***********\n\n");
	fprintf(fgen, "kb      nsnps       qneu       r2         pi         TajD        c          Bs         numgenes   length       del         qdel          hdel              sdel       qn         qqn          hqn              sqn           let        qlet          hlet              slet        od         qod                   hod              sod              adv         qadv                   hadv              sadv              VA             VAdel         VAqn         VAlet            VAod          VAadv          VD             VDdel            VDqn            VDlet            VDod      VDadv      B              Bdel            Bqn            Blet            Bod            Badv\n");
	for (j=0; j<(nwindows/window_size); j++) fprintf(fgen, "%d     %6.4f     %6.4f     %6.4f     %6.4f     %6.4f     %6.4f     %6.4f     %6.4f     %6.4f     %6.4f     %6.4f     %6.4f     %10.8f     %6.4f     %6.4f     %6.4f     %6.4f     %6.4f     %6.4f     %6.4f     %6.4f     %6.4f     %6.4f     %6.4f     %6.4f     %6.4f     %6.4f     %6.4f     %6.4f     %10.8f     %10.8f     %10.8f     %10.8f     %10.8f     %10.8f     %10.8f     %10.8f     %10.8f     %10.8f     %10.8f     %10.8f     %10.8f     %10.8f     %10.8f     %10.8f     %10.8f     %10.8f\n", (j+1)*100, accmean(&AVE_neuwin[j]), accmean(&AVE_qneuwin[j]), accmean(&AVE_r2win[j]), accmean(&AVE_piwin[j]), accmean(&AVE_TajDwin[j]), accmean(&AVE_cwin[j]), accmean(&Bswin[j]), accsum(&numgenes[j]), accmean(&length[j]), accmean(&AVE_delwin[j]), accmean(&AVE_qdelwin[j]), accmean(&AVE_hsdelwin[j]), accmean(&AVE_sdelwin[j]), accmean(&AVE_qnwin[j]), accmean(&AVE_qqnwin[j]), accmean(&AVE_hsqnwin[j]), accmean(&AVE_sqnwin[j]), accmean(&AVE_letwin[j]), accmean(&AVE_qletwin[j]), accmean(&AVE_hsletwin[j]), accmean(&AVE_sletwin[j]), accmean(&AVE_odwin[j]), accmean(&AVE_qodwin[j]), accmean(&AVE_hsodwin[j]), accmean(&AVE_sodwin[j]), accmean(&AVE_advwin[j]), accmean(&AVE_qadvwin[j]), accmean(&AVE_hsadvwin[j]), accmean(&AVE_sadvwin[j]), TOT_VA[j], accmean(&AVE_va_delwin[j]), accmean(&AVE_va_qnwin[j]), accmean(&AVE_va_letwin[j]), accmean(&AVE_va_odwin[j]), accmean(&AVE_va_advwin[j]), TOT_VD[j], accmean(&AVE_vd_delwin[j]), accmean(&AVE_vd_qnwin[j]), accmean(&AVE_vd_letwin[j]), accmean(&AVE_vd_odwin[j]), accmean(&AVE_vd_advwin[j]), TOT_ID[j], accmean(&AVE_id_delwin[j]), accmean(&AVE_id_qnwin[j]), accmean(&AVE_id_letwin[j]), accmean(&AVE_id_odwin[j]), accmean(&AVE_id_advwin[j]));

	fclose(fgen);

	ftable = fopen ("table.txt","w");
	for (j=0; j<(nwindows/window_size); j++) fprintf(ftable, "%d     %6.4f     %6.4f     %6.4f     %6.4f     %6.4f     %6.4f     %6.4f     %6.4f     %6.4f     %6.4f     %6.4f     %6.4f     %10.8f     %6.4f     %6.4f     %6.4f     %6.4f     %6.4f     %6.4f     %6.4f     %6.4f     %6.4f     %6.4f     %6.4f     %6.4f     %6.4f     %6.4f     %6.4f     %6.4f     %10.8f     %10.8f     %10.8f     %10.8f     %10.8f     %10.8f     %10.8f     %10.8f     %10.8f     %10.8f     %10.8f     %10.8f     %10.8f     %10.8f     %10.8f     %10.8f     %10.8f     %10.8f\n", (j+1)*100, accmean(&AVE_neuwin[j]), accmean(&AVE_qneuwin[j]), accmean(&AVE_r2win[j]), accmean(&AVE_piwin[j]), accmean(&AVE_TajDwin[j]), accmean(&AVE_cwin[j]), accmean(&Bswin[j]), accsum(&numgenes[j]), accmean(&length[j]), accmean(&AVE_delwin[j]), accmean(&AVE_qdelwin[j]), accmean(&AVE_hsdelwin[j]), accmean(&AVE_sdelwin[j]), accmean(&AVE_qnwin[j]), accmean(&AVE_qqnwin[j]), accmean(&AVE_hsqnwin[j]), accmean(&AVE_sqnwin[j]), accmean(&AVE_letwin[j]), accmean(&AVE_qletwin[j]), accmean(&AVE_hsletwin[j]), accmean(&AVE_sletwin[j]), accmean(&AVE_odwin[j]), accmean(&AVE_qodwin[j]), accmean(&AVE_hsodwin[j]), accmean(&AVE_sodwin[j]), accmean(&AVE_advwin[j]), accmean(&AVE_qadvwin[j]), accmean(&AVE_hsadvwin[j]), accmean(&AVE_sadvwin[j]), TOT_VA[j], accmean(&AVE_va_delwin[j]), accmean(&AVE_va_qnwin[j]), accmean(&AVE_va_letwin[j]), accmean(&AVE_va_odwin[j]), accmean(&AVE_va_advwin[j]), TOT_VD[j], accmean(&AVE_vd_delwin[j]), accmean(&AVE_vd_qnwin[j]), accmean(&AVE_vd_letwin[j]), accmean(&AVE_vd_odwin[j]), accmean(&AVE_vd_advwin[j]), TOT_ID[j], accmean(&AVE_id_delwin[j]), accmean(&AVE_id_qnwin[j]), accmean(&AVE_id_letwin[j]), accmean(&AVE_id_odwin[j]), accmean(&AVE_id_advwin[j]));
                                                     
	fclose(ftable);
}

/* ***************************************************** */

lookfortext(s)
char *s;
{
   int len, i, curchar;
   char c;

   curchar = 0;
   len = 0;

   for (i=0; i<=100; i++)
   {
      if (s[i] == '\0') break;
      len++;
   }
   do
   {
      c = getc(fpop);

      if (c==s[curchar])
      {
         curchar++;
         if (curchar==len) return(0);
      }
      else curchar = 0;
   }
   while (c != EOF);
}

/* ********************************************************************************************* */
