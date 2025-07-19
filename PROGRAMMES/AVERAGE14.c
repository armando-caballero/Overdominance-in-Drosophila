// AVERAGE14.c

#include "libhdr"

int i, x, r, a, n, REPS, nwindows, window[300], sum_ns_pb;
double w, y, z, ppi[300], ppi_corr[300], c[300];
double ns[300], qns[300], ndel[300], qndel[300], nqn[300], qnqn[300], nlet[300], qnlet[300], nod[300], qnod[300], nadv[300], qnadv[300], nr2[300];
double tot_nsnps, tot_qneu, tot_r2, tot_ppi, tot_TajD, tot_del, tot_qdel, tot_hsdel, tot_sdel, tot_qn, tot_qqn, tot_hsqn, tot_sqn, tot_let, tot_qlet, tot_hslet, tot_slet, tot_od, tot_qod, tot_hsod, tot_sod, tot_adv, tot_qadv, tot_hsadv, tot_sadv, tot_VA, tot_VAdel, tot_VAqn, tot_VAlet, tot_VAod, tot_VAadv, tot_VD, tot_VDdel, tot_VDqn, tot_VDlet, tot_VDod, tot_VDadv, tot_B, tot_Bdel, tot_Bqn, tot_Blet, tot_Bod, tot_Badv;
char filename[20];

double sum_n, sum_c, sum2_c, sum_pi, sum_pi_c, cov, var_c;
double 	sum_ns, sum_ndel, sum_nqn, sum_nlet, sum_nadv, sum_nod, sum_r2, sum_r2025, sum_ns025;
double 	sum_qns, sum_qndel, sum_qnqn, sum_qnlet, sum_qnadv, sum_qnod;
double a1, a2, b1, b2, c1, c2, e1, e2;

struct acc nsnps[300], r2[300], ppiwin[300], TajD[300], cwin[300], Bs[300], numgenes[300], length[300], del[300], qdel[300], qlet[300], qod[300], qqn[300], qneu[300], hsdel[300], sdel[300], let[300], hslet[300], slet[300], od[300], hsod[300], sod[300], qn[300], hsqn[300], sqn[300], VA[300], VAdel[300], VAqn[300], VAlet[300], VAod[300], VAadv[300], VD[300], VDdel[300], VDqn[300], VDlet[300], VDod[300], VDadv[300], B[300], Bdel[300], Bqn[300], Blet[300], Bod[300], Badv[300], adv[300], hsadv[300], sadv[300], qadv[300];
struct acc bb[10001];

FILE *fout, *fse, *ftab, *fpicorr01, *fDTcorr01, *fpicorr025, *fDTcorr025, *fpicorr, *fnq, *fr2, *fr2c025, *fc025, *fc025avg, *ftabcpi;

main()
{
	getintandskip("NREPS:",&REPS,1,infinity);
	getintandskip("NWINDOWS:",&nwindows,1,infinity);

	if (nwindows == 230) sum_ns_pb = 442530;
	if (nwindows == 211) sum_ns_pb = 329499;
	if (nwindows == 245) sum_ns_pb = 383732;
	if (nwindows == 279) sum_ns_pb = 376041;

	n = 50;
	for(i=1;i<n;i++)
	{
		a1 += 1.0/(double)i;
		a2 += 1.0/(double)(i*i);
	}
	b1 = (n+1.0)/(3.0*(n-1.0));
	b2 = 2.0*((n*n)+n+3.0)/(9.0*n*(n-1.0));
	c1 = b1 - (1.0/a1);
	c2 = b2 - ((n+2.0)/(a1*n)) + (a2/(a1*a1));
	e1 = c1 / a1;
	e2 = c2 / ((a1*a1) + a2);

	fout = fopen ("tableavg.txt","w");
	fse = fopen ("se.txt","w");
    
	fpicorr = fopen ("tablepicorr.txt","w");
	fr2 = fopen ("tabler2.txt","w");
	fr2c025 = fopen ("tabler2c025.txt","w");
    
	fnq = fopen ("tablenq.txt","w");
    
	fpicorr01 = fopen ("tablepicorr01.txt","w");
	fDTcorr01 = fopen ("tableDTcorr01.txt","w");
    
	fpicorr025 = fopen ("tablepicorr025.txt","w");
	fDTcorr025 = fopen ("tableDTcorr025.txt","w");
    
	fc025 = fopen ("tablec025.txt","w");
	fc025avg = fopen ("tablec025avg.txt","w");
    
	fprintf(fc025, "c   nsnps   r2  pi_corr\n");

	for (r=1; r<=REPS; r++)
	{
		sprintf(filename, "table%d.txt", r);
		ftab = fopen(filename, "r");

		sprintf(filename, "table_C_PI%d.txt", r);
		ftabcpi = fopen(filename, "w");

		sum_ns = 0.0;
		sum_ndel = 0.0;
		sum_nqn = 0.0;	
		sum_nlet = 0.0;
		sum_nod = 0.0;
		sum_nadv = 0.0;
        
		sum_qns = 0.0;
		sum_qndel = 0.0;
		sum_qnqn = 0.0;
 		sum_qnlet = 0.0;
        		sum_qnod = 0.0;
        		sum_qnadv = 0.0;
        
       		 sum_r2 = 0.0;

 		for (i=1; i<=nwindows; i++)
		{
			fscanf(ftab,"%d", &x);
			window[i] = x;

 			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)	
			{
				accum(&nsnps[i], w);
				ns[i] = w;
				sum_ns += ns[i];
			}
            
			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)	
			{
				accum(&qneu[i], w);
				qns[i] = w;
				sum_qns += qns[i]*ns[i];
			}

			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)	
			{
				accum(&r2[i], w);
				nr2[i] = w;
				sum_r2 += ns[i]*nr2[i];
			}

			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a > -100000)
			{
				ppi[i]=w;
			}
			else ppi[i] = -999;
            
			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)	accum(&TajD[i], w);
            
			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)
			{
				accum(&cwin[i], w);
				c[i] = w;
			}
			else c[i] = -999;

			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)	accum(&Bs[i], w);
            
			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)	accum(&numgenes[i], w);
            
 			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)	accum(&length[i], w);
            
			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)
			{
				accum(&del[i], w);
				ndel[i] = w;
				sum_ndel += ndel[i];
			}
            
			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)
			{
				accum(&qdel[i], w);
				qndel[i] = w;
				sum_qndel += qndel[i]*ndel[i];
			}
            
 			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)	accum(&hsdel[i], w);
            
			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)	accum(&sdel[i], w);
//			if (a < 0)	printf("r=%d i=%d w=%f a=%d %d\n", r, i, w, a, isdigit(a));
            
			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)
			{
				accum(&qn[i], w);
 				nqn[i] = w;
				sum_nqn += nqn[i];
			}
                
			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)
 			{
				accum(&qqn[i], w);
				qnqn[i] = w;
				sum_qnqn += qnqn[i]*nqn[i];
			}
            
			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)	accum(&hsqn[i], w);
            
			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)	accum(&sqn[i], w);
            
			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)
			{
				accum(&let[i], w);
				nlet[i] = w;
				sum_nlet += nlet[i];
			}
            
			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)
			{
				accum(&qlet[i], w);
				qnlet[i] = w;
 				sum_qnlet += qnlet[i]*nlet[i];
			}	
            
			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)	accum(&hslet[i], w);
            
			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)	accum(&slet[i], w);
            
			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)	
			{
				accum(&od[i], w);
				nod[i] = w;
				sum_nod += nod[i];
			}
            
  			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)	
 			{
				accum(&qod[i], w);
				qnod[i] = w;
 				sum_qnod += qnod[i]*nod[i];
			}
            
			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)	accum(&hsod[i], w);
            
			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)	accum(&sod[i], w);
            
			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)	
			{
				accum(&adv[i], w);
				nadv[i] = w;
				sum_nadv += nadv[i];
			}
            
			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)	
			{
				accum(&qadv[i], w);
				qnadv[i] = w;
				sum_qnadv += qnadv[i]*nadv[i];
			}
            
			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)	accum(&hsadv[i], w);
            
			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)	accum(&sadv[i], w);
            
			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)	accum(&VA[i], w);
            
			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)	accum(&VAdel[i], w);
            
			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)	accum(&VAqn[i], w);
            
			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)	accum(&VAlet[i], w);
        
			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)	accum(&VAod[i], w);
            
			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)	accum(&VAadv[i], w);
            
			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)	accum(&VD[i], w);
            
			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)	accum(&VDdel[i], w);
            
			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)	accum(&VDqn[i], w);
            
			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)	accum(&VDlet[i], w);
            
 			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)	accum(&VDod[i], w);
            
			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)	accum(&VDadv[i], w);
            
			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)	accum(&B[i], w);
            
			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)	accum(&Bdel[i], w);
            
			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)	accum(&Bqn[i], w);
            
			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)	accum(&Blet[i], w);
            
			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)	accum(&Bod[i], w);
            
			fscanf(ftab,"%lf", &w);
			a = (int)w;
			if (a >= -100000)	accum(&Badv[i], w);
		}
        
		int sum_win = 0;
		int sum_win01 = 0;
		int sum_win025 = 0;
		int sum_winc025 = 0;
		int sum_winr2 = 0;
        
		double sum_ppi = 0.0;
		double sum_ppi01 = 0.0;
		double sum_ppi025 = 0.0;
        
		double sum_D01 = 0.0;
		double sum_D025 = 0.0;
        
		double sum_r2025 = 0.0;
		double sum_ns025 = 0.0;
        
		double S_corr[300], d_corr[300], D_corr[300];

		for (i=1; i<=nwindows; i++)
		{
			if (ppi[i] != -999)
			{
				ppi_corr[i] = ((ppi[i]*ns[i]*(sum_ns_pb/sum_ns))/100000);
				S_corr[i] = (ns[i]*(sum_ns_pb/sum_ns))/100000;
				d_corr[i] = ppi_corr[i] - (S_corr[i] / a1);
				D_corr[i] = d_corr[i] / sqrt((e1*S_corr[i]) + (e2*S_corr[i]*(S_corr[i]-1.0)));
				accum(&ppiwin[i], ppi_corr[i]);
			}
			else ppi_corr[i] = -999;


			fprintf(ftabcpi, "%f %10.8f\n", c[i], ppi_corr[i]);


			if ((c[i] <= 0.1)&&(ppi_corr[i] != -999))
			{
				sum_win01 ++;
				sum_ppi01 += ppi_corr[i];
				sum_D01 += D_corr[i];
			}
        
			if ((c[i] <= 0.25)&&(ppi_corr[i] != -999))
			{
				sum_win025 ++;
				sum_ppi025 += ppi_corr[i];
				sum_D025 += D_corr[i];
			}
        
			if (ppi_corr[i] != -999)
			{
				sum_win ++;
				sum_ppi += ppi_corr[i];
			}
        
			if (c[i] <= 0.25)
			{
				sum_winc025 ++;
				sum_ns025 += ns[i];
				sum_r2025 += nr2[i]*ns[i];
			}
        
			if (c[i] <= 0.25)
			{
				fprintf(fc025, "%f %f %f %f\n", c[i], ns[i], nr2[i], ppi_corr[i]);
			}
		}
        
		fprintf(fpicorr01, "%f %d\n", sum_ppi01/(double)sum_win01, sum_win01);
		fprintf(fDTcorr01, "%f\n", sum_D01/(double)sum_win01);
        
		fprintf(fpicorr025, "%f %d\n", sum_ppi025/(double)sum_win025, sum_win025);
		fprintf(fDTcorr025, "%f\n", sum_D025/(double)sum_win025);
        
  		fprintf(fpicorr, "%f %d\n", sum_ppi/(double)sum_win, sum_win);
        
		fprintf(fr2, "%f\n", sum_r2/(double)sum_ns);
        
		fprintf(fr2c025, "%f %d\n", sum_r2025/sum_ns025, sum_winc025);
        
		fprintf(fnq, "%f    %f    %f    %f    %f    %f    %f    %f    %f    %f    %f    %f\n", sum_ns, sum_qns/sum_ns, sum_ndel, sum_qndel/sum_ndel, sum_nqn, sum_qnqn/sum_nqn, sum_nlet, sum_qnlet/sum_nlet, sum_nod, sum_qnod/sum_nod, sum_nadv, sum_qnadv/sum_nadv);
    
		sum_n = 0.0;
		sum_c = 0.0;
		sum2_c = 0.0;
		sum_pi = 0.0;
		sum_pi_c = 0.0;
        
		for (i=1; i<=nwindows; i++)
		{
			if ((ppi_corr[i] != -999) && (c[i] != -999)&& (c[i] <= 0.25))
			{
				sum_n ++;
				sum_c += c[i];
				sum2_c += c[i]*c[i];
				sum_pi += ppi_corr[i];
				sum_pi_c += c[i]*ppi_corr[i];
			}
		}

		cov = (sum_pi_c/sum_n) - ((sum_pi/sum_n) * (sum_c/sum_n));
		var_c = (sum2_c/sum_n) - ((sum_c/sum_n) * (sum_c/sum_n));
		accum(&bb, (cov / var_c));

		fclose(ftab);
	}
    
	/* ***************** printout ******************** */

	fprintf(fc025avg, "c   nsnps   r2  pi_corr\n");
	for (i=1; i<=nwindows; i++)
	{
		if (accmean(&cwin[i]) <= 0.25)
		{
			fprintf(fc025avg, "%f %f %f\n", accmean(&cwin[i]), accmean(&r2[i]), accmean(&ppiwin[i]));
		}
	}
	fclose(fc025avg);

	fprintf(fout, "kb      nsnps       qneu         r2         pi         TajD        c          Bs         numgenes   length       del         qdel          hdel              sdel           qn         qqn          hqn              sqn           let        qlet          hlet              slet          od         qod          hod              sod          adv         qadv          hadv              sadv          VA             VAdel         VAqn         VAlet            VAod          VAadv          VD             VDdel             VDqn            VDlet            VDod      VDadv      B              Bdel              Bqn            Blet            Bod            Badv\n");
	for (i=1; i<=nwindows; i++)
	{
		fprintf(fout, "%d %6.4f   %6.4f   %6.4f   %10.8f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %10.8f   %10.8f   %10.8f    %10.8f  %10.8f   %10.8f    %10.8f  %10.8f   %10.8f   %10.8f   %10.8f   %10.8f    %10.8f  %10.8f   %10.8f    %10.8f    %10.8f    %10.8f    %10.8f\n",
		window[i], accmean(&nsnps[i]), accmean(&qneu[i]), accmean(&r2[i]), accmean(&ppiwin[i]), accmean(&TajD[i]), accmean(&cwin[i]), accmean(&Bs[i]), accmean(&numgenes[i]), accmean(&length[i]), accmean(&del[i]), accmean(&qdel[i]), accmean(&hsdel[i]), accmean(&sdel[i]), accmean(&qn[i]), accmean(&qqn[i]), accmean(&hsqn[i]), accmean(&sqn[i]),accmean(&let[i]), accmean(&qlet[i]), accmean(&hslet[i]), accmean(&slet[i]), accmean(&od[i]), accmean(&qod[i]), accmean(&hsod[i]), accmean(&sod[i]), accmean(&adv[i]), accmean(&qadv[i]), accmean(&hsadv[i]), accmean(&sadv[i]),accmean(&VA[i]), accmean(&VAdel[i]), accmean(&VAqn[i]), accmean(&VAlet[i]), accmean(&VAod[i]), accmean(&VAadv[i]), accmean(&VD[i]), accmean(&VDdel[i]), accmean(&VDqn[i]), accmean(&VDlet[i]), accmean(&VDod[i]), accmean(&VDadv[i]), accmean(&B[i]), accmean(&Bdel[i]), accmean(&Bqn[i]), accmean(&Blet[i]), accmean(&Bod[i]), accmean(&Badv[i]));
	}
   
	for (i=1; i<=nwindows; i++) 
	{
		if (accmean(&nsnps[i]) > 0.0) tot_nsnps += accmean(&nsnps[i]);
        
		if (accmean(&nsnps[i]) > 0.0) tot_qneu += (accmean(&qneu[i])*accmean(&nsnps[i]));
		if (accmean(&nsnps[i]) > 0.0) tot_r2 += (accmean(&r2[i])*accmean(&nsnps[i]));
		if (accmean(&nsnps[i]) > 0.0) tot_ppi += (accmean(&ppiwin[i]));
		if (accmean(&nsnps[i]) > 0.0) tot_TajD += (accmean(&TajD[i])*accmean(&nsnps[i]));
        
		if (accmean(&del[i]) > 0.0) tot_del += accmean(&del[i]);
		//printf("i=%d del=%f tot_del=%f\n", i, accmean(&del[i]), tot_del);
		if (accmean(&let[i]) > 0.0) tot_let += accmean(&let[i]);
		if (accmean(&od[i]) > 0.0) tot_od += accmean(&od[i]);
		if (accmean(&qn[i]) > 0.0) tot_qn += accmean(&qn[i]);
		if (accmean(&adv[i]) > 0.0) tot_adv += accmean(&adv[i]);
        
		if (accmean(&del[i]) > 0.0) 
		{
			tot_qdel += (accmean(&qdel[i])*accmean(&del[i]));
			//printf("i=%d del=%f qdel=%f tot_qdel=%f\n", i, accmean(&del[i]), accmean(&qdel[i]), tot_qdel);
			tot_hsdel += (accmean(&hsdel[i])*accmean(&del[i]));
			//printf("i=%d del=%f hsdel=%f tot_hsdel=%f\n", i, accmean(&del[i]), accmean(&hsdel[i]), tot_hsdel);
			tot_sdel += (accmean(&sdel[i])*accmean(&del[i]));
			//printf("i=%d del=%f s=%f tot_s=%f\n", i, accmean(&del[i]), accmean(&s[i]), tot_s);
            
			if (accmean(&VAdel[i]) > 0.0) tot_VAdel += (accmean(&VAdel[i]));
			if (accmean(&VDdel[i]) > 0.0) tot_VDdel += (accmean(&VDdel[i]));
			if (accmean(&Bdel[i]) > 0.0) tot_Bdel += (accmean(&Bdel[i]));
			//printf("i=%d del=%f Bdel=%f tot_Bdel=%f\n", i, accmean(&del[i]), accmean(&Bdel[i]), tot_Bdel);
		}
        
		if (accmean(&qn[i]) > 0.0) 
		{
			tot_qqn += (accmean(&qqn[i])*accmean(&qn[i]));
			tot_hsqn += (accmean(&hsqn[i])*accmean(&qn[i]));
			tot_sqn += (accmean(&sqn[i])*accmean(&qn[i]));
            
			if (accmean(&VAqn[i]) > 0.0) tot_VAqn += (accmean(&VAqn[i]));
			if (accmean(&VDqn[i]) > 0.0) tot_VDqn += (accmean(&VDqn[i]));
			if (accmean(&Bqn[i]) > 0.0) tot_Bqn += (accmean(&Bqn[i]));
		}
              
		if (accmean(&let[i]) > 0.0) 
		{
			tot_qlet += (accmean(&qlet[i])*accmean(&let[i]));
			tot_hslet += (accmean(&hslet[i])*accmean(&let[i]));
			tot_slet += (accmean(&slet[i])*accmean(&let[i]));
            
			if (accmean(&VAlet[i]) > 0.0) tot_VAlet += (accmean(&VAlet[i]));
			if (accmean(&VDlet[i]) > 0.0) tot_VDlet += (accmean(&VDlet[i]));
			if (accmean(&Blet[i]) > 0.0) tot_Blet += (accmean(&Blet[i]));
		}
        
		if (accmean(&od[i]) > 0.0) 
		{
			tot_qod += (accmean(&qod[i])*accmean(&od[i]));
			tot_hsod += (accmean(&hsod[i])*accmean(&od[i]));
			tot_sod += (accmean(&sod[i])*accmean(&od[i]));
            
			if (accmean(&VAod[i]) > 0.0) tot_VAod += (accmean(&VAod[i]));
			if (accmean(&VDod[i]) > 0.0) tot_VDod += (accmean(&VDod[i]));
			if (accmean(&Bod[i]) > 0.0) tot_Bod += (accmean(&Bod[i])); 
		}
        
		if (accmean(&adv[i]) > 0.0) 
		{
 			tot_qadv += (accmean(&qadv[i])*accmean(&adv[i]));
			tot_hsadv += (accmean(&hsadv[i])*accmean(&adv[i]));
			tot_sadv += (accmean(&sadv[i])*accmean(&adv[i]));
            
			if (accmean(&VAadv[i]) > 0.0) tot_VAadv += (accmean(&VAadv[i]));
			if (accmean(&VDadv[i]) > 0.0) tot_VDadv += (accmean(&VDadv[i]));
			if (accmean(&Badv[i]) > 0.0) tot_Badv += (accmean(&Badv[i])); 
		}
	}
    
	fprintf(fout, "\n********AVERAGE*******\n");
	fprintf(fout, "nsnps       qneu       r2       pi      b       sd(b)     se(b)     TajD       del        qdel        hdel        sdel         qn        qqn        hqn        sqn         let        qlet        hlet        slet                od        qod        hod        sod       adv        qadv        hadv        sadv       VA           VAdel        VAqn        VAlet         VAod        VAadv        VD           VDdel           VDqn         VDlet       VDod         VDadv         B             Bdel             Bqn        Blet         Bod         Badv\n");
	fprintf(fout, "%6.4f   %6.4f  %6.4f   %10.8f   %10.8f   %10.8f    %10.8f    %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %10.8f   %10.8f   %10.8f    %10.8f  %10.8f   %10.8f    %10.8f  %10.8f   %10.8f    %10.8f    %10.8f    %10.8f    %10.8f  %10.8f   %10.8f   %10.8f   %10.8f   %10.8f\n",tot_nsnps, (tot_qneu/tot_nsnps), (tot_r2/tot_nsnps), (tot_ppi/nwindows), accmean(&bb), sqrt(variance(&bb)), se(&bb), (tot_TajD/tot_nsnps), tot_del, (tot_qdel/tot_del), (tot_hsdel/tot_del), (tot_sdel/tot_del), tot_qn, (tot_qqn/tot_qn), (tot_hsqn/tot_qn), (tot_sqn/tot_qn), tot_let, (tot_qlet/tot_let), (tot_hslet/tot_let), (tot_slet/tot_let), tot_od, (tot_qod/tot_od), (tot_hsod/tot_od), (tot_sod/tot_od), tot_adv, (tot_qadv/tot_adv), (tot_hsadv/tot_adv), (tot_sadv/tot_adv), (tot_VAdel+tot_VAqn+tot_VAlet+tot_VAod+tot_VAadv), tot_VAdel, tot_VAqn, tot_VAlet, tot_VAod, tot_VAadv, (tot_VDdel+tot_VDqn+tot_VDlet+tot_VDod+tot_VDadv), tot_VDdel, tot_VDqn, tot_VDlet, tot_VDod, tot_VDadv, (tot_Bdel+tot_Bqn+tot_Blet+tot_Bod+tot_Badv), tot_Bdel, tot_Bqn, tot_Blet, tot_Bod, tot_Badv);
    
	fclose(fout);
    
	fprintf(fse, "kb      nsnps       qneu         r2         pi         TajD        c          Bs         numgenes   length       del         qdel          hdel              sdel           qn         qqn          hqn              sqn           let        qlet          hlet              slet          od         qod          hod              sod          adv         qadv          hadv              sadv          VA             VAdel         VAqn         VAlet            VAod          VAadv          VD             VDdel             VDqn            VDlet            VDod      VDadv      B              Bdel              Bqn            Blet            Bod            Badv\n");
	for (i=1; i<=nwindows; i++)
	{
		fprintf(fse, "%d %6.4f   %6.4f   %10.8f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %6.4f   %10.8f   %10.8f   %10.8f    %10.8f  %10.8f   %10.8f    %10.8f  %10.8f   %10.8f   %10.8f   %10.8f   %10.8f    %10.8f  %10.8f   %10.8f    %10.8f    %10.8f    %10.8f    %10.8f\n",
		window[i], se(&nsnps[i]), se(&qneu[i]), se(&r2[i]), se(&ppiwin[i]), se(&TajD[i]), se(&cwin[i]), se(&Bs[i]), se(&numgenes[i]), se(&length[i]), se(&del[i]), se(&qdel[i]), se(&hsdel[i]), se(&sdel[i]), se(&qn[i]), se(&qqn[i]), se(&hsqn[i]), se(&sqn[i]), se(&let[i]), se(&qlet[i]), se(&hslet[i]), se(&slet[i]), se(&od[i]), se(&qod[i]), se(&hsod[i]), se(&sod[i]), se(&adv[i]), se(&qadv[i]), se(&hsadv[i]), se(&sadv[i]), se(&VA[i]), se(&VAdel[i]), se(&VAqn[i]), se(&VAlet[i]), se(&VAod[i]), se(&VAadv[i]), se(&VD[i]), se(&VDdel[i]), se(&VDqn[i]), se(&VDlet[i]), se(&VDod[i]), se(&VDadv[i]), se(&B[i]), se(&Bdel[i]), se(&Bqn[i]), se(&Blet[i]), se(&Bod[i]), se(&Badv[i]));
	}
    
	fclose(fse);                                         
}

