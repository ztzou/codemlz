/* codemlz.c

   Maximum likelihood parameter estimation for codon sequences (seqtype=1)

         Modified from codeml (Copyright, Ziheng YANG, 1993-2003)
                         Zhengting Zou, 2016

               cc -o codeml -fast codeml.c tools.o -lm
                         codemlz <ControlFileName>

    This script tries to use the frame work of original codeml program to realize
    the ML estimation of a new parameter eta, which is the selection difference
    between non-synonymous transitions and non-synonymous transversions.
*/

#include "paml.h"

#define NS            7000
#define NBRANCH       (NS*2-2)
#define NNODE         (NS*2-1)
#define MAXNSONS      100
#define NGENE         2000
#define LSPNAME       50
#define NCODE         64
#define NCATG         40
#define NBTYPE        17

#define NP            (NBRANCH*2+NGENE-1+2+NCODE+2)
/*
#define NP            (NBRANCH+NGENE-1+189+2+NCODE+2)
*/
extern char BASEs[],AAs[];
extern int noisy, NFunCall, NEigenQ, NPMatUVRoot, *ancestor, GeneticCode[][64];
extern double *SeqDistance;
extern double SS,NN,Sd,Nd; /* kostas, SS=# of syn. sites, NN=# of non-syn. sites, Sd=# of syn. subs., Nd=# of non-syn. subs. as defined in DistanceMatNG86 in treesub.c */

//int  Forestry (FILE *fout);
int  GetMemPUVR(int nc, int nUVR);
int  sortwM3(double x[]);
void DetailOutput(FILE *fout, double x[], double var[]);
int  GetOptions (char *ctlf);
int  testx (double x[], int np);
int  SetxBound (int np, double xb[][2]);
int  SetxInitials (int np, double x[], double xb[][2]);
int  GetInitials (double x[], int*fromfile);
double *PointKappa (double xcom[], int igene);
double *PointOmega (double xcom[], int igene, int inode, int isiteclass);
int  GetCodonFreqs (void);
int  SetParameters (double x[]);
int  SetParametersNSsites (double x[]);
int  Set_UVR_BranchSite (int iclass, int branchlabel);
int  SetPGene (int igene, int _pi, int _UVRoot, int _alpha, double x[]);
int  SetPSiteClass(int iclass, double x[]);
int  PMatJC69like (double P[], double t, int n);
int  printfcode (FILE *fout, double fb61[], double space[]);
int  InitializeCodon (FILE *fout, double space[]);
int  AA2Codonf (double faa[20], double fcodon[]);
int  DistanceMatAA (FILE *fout);
int  GetDaa(FILE *fout, double daa[]);
void getpcodonClass(double x[], double pcodonClass[]);
int  SelectionCoefficients (FILE* fout, double kappa[], double ppi[], double omega);
int  eigenQcodon(int mode, double blength, double *S, double *dS, double *dN,
     double Root[], double U[], double V[], double *meanrate, double kappa[], double omega, double eta, double Q[]);
int  eigenQaa(FILE *fout, double Root[], double U[], double V[],double rate[]);
int  Qcodon2aa(double Qc[], double pic[], double Qaa[], double piaa[]);
int  SetAA1STEP(void);
int  GetOmegaAA(int OmegaAA[]);
int  TestModelQc(FILE *fout, double x[]);
double lfun2dSdN(double x[], int np);
int  VariancedSdN(double t, double omega, double vtw[2*2], double vdSdN[2*2]);
int  GetCodonFreqs2 (void);
int  PairwiseCodon(FILE *fout, FILE*fds, FILE*fdn, FILE*dt, double space[]);
int  PairwiseAA(FILE *fout, FILE *f2AA);
int  lfunNSsites_rate(FILE* fout, double x[], int np);
int  lfunNSsites_M2M8(FILE* frst, double x[], int np);
int  lfunNSsites_AC(FILE* frst, double x[], int np);
double GetBranchRate(int igene, int ibrate, double x[], int *ix);
int  GetPMatBranch(double Pt[], double x[], double t, int inode);
int  ConditionalPNode(int inode, int igene, double x[]);
double CDFdN_dS(double x,double par[]);
int  DiscreteNSsites(double par[]);
char GetAASiteSpecies(int species, int sitepatt);
void finishup(void);
int  mergeSeqs(FILE*fout);
void Get4foldSites(void);
int  AdHocRateSmoothing(FILE*fout, double x[NS*3], double xb[NS*3][2], double space[]);
void DatingHeteroData(FILE* fout);

int SlidingWindow(FILE*fout, FILE* fpair[], double space[]);

void SimulateData2s61(void);
void Ina(void);
void d4dSdN(FILE*fout);

//kostas functions
double logprior(double t, double w, double par[]);
double logistic_transformation(double point, double logmean, double stdlogpar);
double loglikelihood(double Ptmatrix[]);
int EstVariances(double *var);
int BayesPairwise(int is, int js, double x[], double var[], double maxlogl,
                    int npoints, double xb[][2], double space[]);
double logP(double x[], int np);
double CDFLogis( double x, double m, double s );
//end of kostas functions

struct common_info {
   unsigned char *z[NS];
   char *spname[NS], seqf[512],outf[512],treef[512],daafile[512], cleandata;
   char oldconP[NNODE];       /* update conP for nodes? to save computation */
   int seqtype, ns, ls, ngene, posG[NGENE+1], lgene[NGENE], npatt,*pose, readpattern;
   int runmode,clock, verbose,print, codonf,aaDist,model,NSsites;
   int nOmega, nbtype, nOmegaType;  /* branch partition, AA pair (w) partition */
   int method, icode, ncode, Mgene, ndata, bootstrap;
   int fix_rgene,fix_kappa,fix_omega,fix_alpha,fix_eta,fix_rho,nparK,fix_blength,getSE;
   int np, ntime, nrgene, nkappa, npi, nrate, nalpha, ncatG, hkyREV;
   size_t sconP, sspace;
   double *fpatt, *space, kappa,omega,alpha,eta,rho,rgene[NGENE], TipDate, TipDate_TimeUnit;
   double pi[NCODE], piG[NGENE][64], fb61[64];
   double f3x4[NGENE][12], *pf3x4, piAA[20];
   double freqK[NCATG], rK[NCATG], MK[NCATG*NCATG],daa[20*20], *conP, *fhK;
   double (*plfun)(double x[],int np);
   double hyperpar[4]; /* kostas, the hyperparameters for the prior distribution of distance & omega */
   double omega_fix;  /* fix the last w in the NSbranchB, NSbranch2 models
          for lineages.  Useful for testing whether w>1 for some lineages. */
   int     conPSiteClass; /* conPSiteClass=0 if (method==0) and =1 if (method==1)?? */
   int     NnodeScale;
   char   *nodeScale;        /* nScale[ns-1] for interior nodes */
   double *nodeScaleF;       /* nScaleF[npatt] for scale factors */
  /* pomega & pkappa are used to communicate between SetParameters & ConditionalPNode
     & eigenQcodon.  Try to remove them? */
   double *pomega, pkappa[5], *ppi,*peta;
}  com;
struct TREEB {
   int  nbranch, nnode, root, branches[NBRANCH][2];
   double lnL;
}  tree;
struct TREEN {
   int father, nson, sons[MAXNSONS], ibranch, ipop;
   double branch, age, omega, *conP, label;
   char *nodeStr, fossil, usefossil;
}  *nodes, **gnodes, nodes_t[2*NS-1];


/* for sptree.nodes[].fossil: lower, upper, bounds, gamma, inverse-gamma */
enum {LOWER_F=1, UPPER_F, BOUND_F} FOSSIL_FLAGS;
char *fossils[]={" ", "L", "U", "B"};

struct SPECIESTREE {
   int nbranch, nnode, root, nspecies, nfossil;
   struct TREESPN {
      char name[LSPNAME+1], fossil, usefossil;  /* fossil: 0, 1, 2, 3 */
      int father, nson, sons[2];
      double age, pfossil[7];   /* lower and upper bounds or alpha & beta */
      double *lnrates;          /* log rates for loci */
   } nodes[2*NS-1];
}  sptree;
/* all trees are binary & rooted, with ancestors unknown. */

struct DATA { /* locus-specific data and tree information */
   int ns[NGENE], ls[NGENE], npatt[NGENE], ngene, lgene[NGENE];
   int root[NGENE+1], BlengthMethod, fix_nu, nbrate[NGENE], icode[NGENE];
   int datatype[1];
   char   *z[NGENE][NS], cleandata[NGENE];
   char   idaafile[NGENE], daafile[NGENE][40];
   double *fpatt[NGENE], lnpT, lnpR, lnpDi[NGENE];
   double Qfactor[NGENE], pi[NGENE][NCODE];
   double rgene[NGENE], kappa[NGENE], alpha[NGENE], omega[NGENE];
   int NnodeScale[NGENE];
   char *nodeScale[NGENE];    /* nScale[data.ns[locus]-1] for interior nodes */
}  data;

extern double Small_Diff;
int Nsensecodon, FROM61[64], FROM64[64], FourFold[4][4];
int ChangedInIteration;  /* 1: t changed, update P(t); 2: paras changed, update UVRoot */
double *PMat, *U, *V, *Root, *_UU[NBTYPE+2], *_VV[NBTYPE+2], *_Root[NBTYPE+2];
/* 5 sets for branchsite models (YN2002); 6 sets for clade models */

double pcodon0[64],paa0[20], *pcodonClass;  /* for aaDist=FIT1 or FIT2 */

int BayesEB;  /* =1 for site models M2a & M8; =2 for branch-site models A & C */
int LASTROUND;
int IClass=-1;

int OmegaAA[190], AA1STEP[190];
enum {DNA, AA, CODON, MORPHC} DATATYPE;

double _rateSite=1;
double Qfactor_NS, Qfactor_NS_branch[NBTYPE];
int KGaussLegendreRule=16;

double AAchem[][20+1]={  /* last element is the max */
{8.1, 10.5, 11.6, 13, 5.5, 10.5, 12.3, 9, 10.4, 5.2,
 4.9, 11.3,  5.7, 5.2,  8,  9.2,  8.6, 5.4, 6.2, 5.9,    13}, /* p */
{ 31, 124,  56,  54,   55, 85, 83,   3, 96, 111,
 111, 119, 105, 132, 32.5, 32, 61, 170, 136, 84,        170}, /* v */
{0, 0.65, 1.33, 1.38, 2.75, 0.89, 0.92, 0.74, 0.58,
 0, 0, 0.33, 0, 0, 0.39, 1.42, 0.71, 0.13, 0.2, 0,      -999},/* c */
{-0.11, 0.079, -0.136, -0.285, -0.184, -0.067, -0.246, -0.073, 0.32, 0.001,
 -0.008, 0.049, -0.041, 0.438, -0.016, -0.153, -0.208, 0.493, 0.381, -0.155} /* a */
};   /* in the order p, v, c, a */


FILE *fout, *frub, *flnf, *frst, *frst1, *frst2=NULL, *finitials;
char *ratef="rates";
enum {Fequal, F1x4, F3x4, Fcodon, F1x4MG, F3x4MG, FMutSel0, FMutSel} CodonFreqs;
char *codonfreqs[]={"Fequal", "F1x4", "F3x4", "Fcodon", "F1x4MG", "F3x4MG", "FMutSel0", "FMutSel"};
enum {NSbranchB=1, NSbranch2, NSbranch3} NSBranchModels;
char *NSbranchmodels[]={"One dN/dS ratio",
     "free dN/dS Ratios for branches", "several dN/dS ratios for branches",
     "NSbranch3"};
enum {Poisson, EqualInput, Empirical, Empirical_F,
     FromCodon=6, REVaa_0=8, REVaa=9} AAModel;
char *aamodels[]={"Poisson", "EqualInput", "Empirical", "Empirical_F", "",
     "", "FromCodon", "", "REVaa_0", "REVaa"};
enum {NSnneutral=1, NSpselection, NSdiscrete, NSfreqs, NSgamma, NS2gamma,
     NSbeta, NSbetaw, NSbetagamma, NSbeta1gamma, NSbeta1normal, NS02normal,
     NS3normal, NSM2aRel=22, NSTgamma, NSTinvgamma, NSTgamma1, NSTinvgamma1} NSsitesModels;
char *NSsitesmodels[]={"one-ratio","NearlyNeutral", "PositiveSelection","discrete","freqs",
     "gamma","2gamma","beta","beta&w>1","beta&gamma", "beta&gamma+1",
     "beta&normal>1", "0&2normal>0", "3normal>0", "", "", "", "", "", "", "", "",
     "M2a_rel", "Tgamma", "Tinvgamma", "Tgamma+1", "Tinvgamma+1"};
int maxNSsitesModels=27;
enum {FIT1=11, FIT2=12} SiteClassModels;
enum {AAClasses=7 } aaDistModels;
char *clockstr[]={"", "Global clock", "Local clock", "ClockCombined"};
enum {GlobalClock=1, LocalClock, ClockCombined} ClockModels;

/* variables for batch run of site models */
int ncatG0=10, insmodel=0, nnsmodels=1, nsmodels[15]={0};
/* used for sliding windows analysis */
int windowsize0=20, offset0=1, npositive=0;
double lnLmodel;

#define CODEML 1
#include "treesub.c"
#include "treespace.c"
#include "codeml_utils.c"

int main (int argc, char *argv[])
{
   FILE *fseq=NULL, *fpair[6];
   char pairfs[6][32]={"2NG.dS","2NG.dN","2NG.t", "2ML.dS","2ML.dN","2ML.t"};
   char ctlf[96]="codeml.ctl", *pmodel, timestr[64];
   char *seqtypestr[3]={"CODONML", "AAML", "CODON2AAML"};
   char *Mgenestr[]={"diff. rate", "separate data", "diff. rate & pi",
                     "diff. rate & k&w", "diff. rate & pi & k&w"};
   int getdistance=1, i, k, s2=0, idata, nc, nUVR, cleandata0;

   starttimer();

   /*
   printf("KGaussLegendreRule? ");
   scanf("%d", &KGaussLegendreRule);
   */
   com.ndata=1;
   noisy=9;           com.runmode=0;
   com.clock=0;       com.fix_rgene=0; /* 0: estimate rate factors for genes */
   com.cleandata=0;  /* 1: delete; 0:use missing data */
   com.seqtype=AAseq;
   com.model=Empirical_F;
   strcpy(com.daafile, "jones.dat");
   com.icode=0;       com.nrate=0;
   com.fix_kappa=0;   com.kappa=1;    com.omega=2.1;
   com.fix_alpha=1;   com.alpha=0.;   com.ncatG=4;   /* alpha=0 := inf */
   com.fix_eta=0;     com.eta=1.;
   com.fix_rho=1;     com.rho=0.;
   com.getSE=0;       com.print=0;    com.verbose=1;  com.fix_blength=0;
   com.method=0;      com.space=NULL;

   frub=gfopen("rub","w");
   frst=gfopen("rst","w");
   frst1=gfopen("rst1","w");

   SetSeed(1, 0);

   if(argc>1) strncpy(ctlf, argv[1], 95);

   GetOptions(ctlf);
   cleandata0 = com.cleandata;
   getdistance = 1;

   fprintf(frst, "Supplemental results for CODEML (seqf: %s  treef: %s)\n",
         com.seqf, com.treef);
   printf("%s in %s\n", seqtypestr[com.seqtype-1], pamlVerStr);

   fout = gfopen(com.outf, "w");

   if(noisy && com.seqtype==CODONseq)
         { printcu(F0,NULL,com.icode); puts("Nice code, uuh?"); }

   nUVR=1; nc=20; /* nc: number of codons */
   if(com.seqtype==CODONseq) {
       nc = 64;
       //if(com.model>=1) nUVR = NBTYPE+2;
   }
   else if (com.seqtype==CODONseq || com.model==FromCodon)
       nc = 64;
   GetMemPUVR(nc, nUVR); /* Get memory for "nUVR sets of matrices" */

   if((fseq=fopen(com.seqf,"r"))==NULL || com.seqf[0]=='\0') { /* open sequence file */
      printf ("\n\nSequence file %s not found!\n", com.seqf);
      exit (-1);
   }

   if (com.aaDist==AAClasses) {
      SetAA1STEP(); /* store whether two aas can exchange within 1 step */
      GetOmegaAA(OmegaAA);
   }

   if(com.seqtype==1) {
         for(i=0; i<3; i++)
            fpair[i]=(FILE*)gfopen(pairfs[i],"w");
   }

   /* Begin processing ndata number of datasets */
   for (idata=0; idata<com.ndata; idata++) {
      if (com.ndata>1) {
         printf ("\nData set %d ", idata+1);
         fprintf(fout, "\n\nData set %d\n", idata+1);
         fprintf(frst,"\t%d",idata+1);
         fprintf(frst1, "%d", idata+1);
         fprintf(frub,"\nData set %2d\n",idata+1);
      }

      if(idata)
         GetOptions(ctlf); /* warning: ndata, cleandata etc. are read again. */

      com.cleandata = cleandata0;
      /* ReadSeq may change seqtype*/
      ReadSeq((com.verbose?fout:NULL), fseq, com.cleandata, 0);
      SetMapAmbiguity();

      fprintf(frst1,"\t%d\t%d\t%d", com.ns, com.ls, com.npatt);
      printf("\ncom.ns=%d\tcom.ls=%d\tcom.npatt=%d\n", com.ns, com.ls, com.npatt);

      if(com.ngene>1) {
         if(com.seqtype==1 && com.npi)
            error2("codon models (estFreq) not implemented for ngene > 1"); /* ngene not for codon model */
            if(com.runmode==-2 && com.Mgene!=1) error2("use Mgene=1 for runmode=-2?");
            if(com.runmode==-3 && com.Mgene!=1) error2("use Mgene=1 for runmode=-3?");
            if(com.model) error2("NSbranchsites with ngene.");
            if(com.NSsites) error2("NSsites with ngene.");
            if(com.aaDist>=FIT1)  /* because of pcodon0[] */
               { error2("ngene for amino acid fitness models"); }
      }

      if(com.ndata==1) fclose(fseq);

      /* Begin reading trees */
      i = (com.ns*2-1)*sizeof(struct TREEN);
      if((nodes=(struct TREEN*)malloc(i))==NULL)
         error2("oom nodes");

      pmodel=(com.seqtype==CODONseq?NSbranchmodels[com.model]:aamodels[com.model]);
      fprintf(fout,"%s (in %s)  %s\n",seqtypestr[com.seqtype-1], pamlVerStr, com.seqf);
      fprintf(fout,"Model: %s for branches, ", pmodel);
      printf("\n%s (in %s)  %s\n",seqtypestr[com.seqtype-1], pamlVerStr, com.seqf);
      printf("\nModel: %s for branches\n", pmodel);
      if(com.seqtype==CODONseq||com.model==FromCodon) {
         if(com.fix_kappa) fprintf(fout, " kappa = %.3f fixed\n", com.kappa);
         if(com.fix_omega) fprintf(fout, " omega = %.3f fixed\n", com.omega);
         if(com.fix_eta) fprintf(fout, " eta = %.3f fixed\n", com.eta);
      }

      if(com.alpha) fprintf (fout, "dGamma (ncatG=%d) ", com.ncatG);
      if(com.ngene>1)
         fprintf (fout, " (%d genes: %s)  ", com.ngene, Mgenestr[com.Mgene]);
      if(com.alpha==0)  com.nalpha=0;
      else              com.nalpha=(com.nalpha?com.ngene:!com.fix_alpha);
      if(com.Mgene==1) com.nalpha=!com.fix_alpha;
      if(com.nalpha>1 && (!com.alpha || com.ngene==1 || com.fix_alpha))
         error2("Malpha");
      if(com.nalpha>1) fprintf (fout,"(%d gamma)", com.nalpha);
      if(com.Mgene && com.ngene==1) error2("Mgene for one gene.");
      if(com.seqtype==CODONseq) {
         fprintf (fout, "\nCodon frequency model: %s\n", codonfreqs[com.codonf]);
         printf ("\nCodon frequency model: %s\n", codonfreqs[com.codonf]);
         if(com.alpha)
            fputs("Warning: Gamma model for codons.  See documentation.",fout);
      }
      if((com.seqtype==CODONseq||com.model==FromCodon)
         && (com.aaDist && com.aaDist<10 && com.aaDist!=AAClasses))
         fprintf(fout,"%s, %s\n",com.daafile,(com.aaDist>0?"geometric":"linear"));

      fprintf(fout,"ns = %3d  ls = %3d\n\n", com.ns, com.ls);

      com.sspace = max2(5000000,3*com.ncode*com.ncode*sizeof(double)); /* 61x61 matrix space */
      printf("com.ncode = %d, com.sspace = %d\n",com.ncode,com.sspace);
      k = com.ns*(com.ns-1)/2;
      if((com.space = (double*)realloc(com.space,com.sspace))==NULL) { /* reallocate to com.space */
         printf("\nfailed to get %9lu bytes for space", com.sspace);
         error2("oom space");
      }
      if(getdistance) { /* reallocate space for seq distance (?) and initiation */
         SeqDistance=(double*)realloc(SeqDistance, k*sizeof(double));
         ancestor=(int*)realloc(ancestor, k*sizeof(int));
         if(SeqDistance==NULL||ancestor==NULL) error2("oom distance&ancestor");
         for(i=0; i<k; i++) SeqDistance[i] = -1;
      }

      if(com.seqtype==AAseq) {
         InitializeBaseAA (fout);
         if (com.model==FromCodon /* ||com.aaDist==AAClasses */)
            AA2Codonf(com.pi, com.fb61);  /* get codon freqs from aa freqs */
      }
      else {  /* codon sequences */
         if(com.sspace < max2(com.ngene+1,com.ns)*(64+12+4)*sizeof(double)) {
            com.sspace = max2(com.ngene+1,com.ns)*(64+12+4)*sizeof(double);
            if((com.space = (double*)realloc(com.space,com.sspace))==NULL)
               error2("oom space for #c");
         }
         if (InitializeCodon(fout,com.space))
            error2("giving up on stop codons");

      } // if(com.seqtype==AAseq) else

      printf("\ncom.pi[:5]: %f %f %f %f %f\n",com.pi[0],com.pi[1],com.pi[2],com.pi[3],com.pi[4]);

      if(getdistance) {
         if(com.seqtype==CODONseq)
            DistanceMatNG86(fout,fpair[0],fpair[1],fpair[2],0);
         else
            DistanceMatAA(fout);
      }
      fflush(fout);

      if(com.alpha || com.NSsites) {
         s2=com.npatt*com.ncatG*sizeof(double);
         if((com.fhK=(double*)realloc(com.fhK,s2))==NULL) error2("oom fhK");
      }

      com.sconP = 2L *com.ncode*com.npatt*sizeof(double);
      /* to be increased later in GetInitials() */
      /* com.sconP = (com.ns-1)*com.ncode*com.npatt*sizeof(double); */
      com.conP = (double*)realloc(com.conP, com.sconP);

      printf("\n%9u bytes for distance",com.ns*(com.ns-1)/2*sizeof(double));
      printf("\n%9u bytes for conP\n", com.sconP);
      printf ("%9u bytes for fhK\n%9u bytes for space\n", s2, com.sspace);
      if(com.conP==NULL)
         error2("oom conP");

      if(com.Mgene) error2("Mgene > 0 not implemented");
      if(com.runmode) {
    	 error2("Runmode != 0 not implemented");
      }
      else {
    	 Forestry(fout);
         printf("\nTime used: %s\n", printtime(timestr));
         fprintf(fout,"\nTime used: %s\n", printtime(timestr));
      }

      FPN(frst);  fflush(frst);
      FPN(frst1); fflush(frst1);
      free(nodes);

      printf("\nDataset finished # %d.\n",idata);

   } // for (idata=0; idata<com.ndata; idata++)

   fclose(frst);
   k=0;
   if(com.seqtype==1)
      k = ((com.runmode==-2 || com.runmode==-3) ? 6 : 3);
   else if (com.runmode==-2)
      k=1;
   FOR(i,k) fclose(fpair[i]);
   if(com.ndata>1 && fseq) fclose(fseq);
   fclose(fout);  fclose(frub);
   if(finitials)  fclose(finitials);
   FreeMemPUVR();
   free(com.pose);
   for(i=0; i<com.ns; i++) free(com.z[i]);

   printf("\nEnd of codemlz.\n");
   return (0);

}

int Forestry (FILE *fout)
{
   printf("\nentering Forestry\n");
   static int times=0;
   FILE *ftree, *frate=NULL;
   int  status=0, i,j=0,k, itree, ntree, np, iteration=1;
   int pauptree=0, haslength;
   double x[NP],xb[NP][2], xcom[NP-NBRANCH], lnL=0,lnL0=0, e=1e-8, tl=0, nchange=-1;
   double *g=NULL, *H=NULL;
#ifdef NSSITESBandits
   FILE *fM0tree;
#endif

   if ((ftree=fopen(com.treef,"r"))==NULL) {
      printf("\ntree file %s not found.\n", com.treef);
      exit(-1);
   }
   GetTreeFileType(ftree, &ntree, &pauptree, 0);
   if (com.alpha)
      frate=(FILE*)gfopen(ratef,"w");
   if (ntree>10 && com.npatt>10000 && com.print)
      puts("\nlnf file may be large");
   flnf=gfopen("lnf","w+");
   fprintf(flnf,"%6d %6d %6d\n", ntree, com.ls, com.npatt);

   /* seems useless for codemlz */
   if(com.seqtype==1 && com.aaDist>=FIT1) {
      xtoy(com.pi,pcodon0,64);
      zero(paa0,20);
      FOR(i,com.ncode) paa0[GeneticCode[com.icode][FROM61[i]]]+=pcodon0[i];
      pcodonClass=(double*)malloc(com.ncatG*64*sizeof(double));
      if(pcodonClass==NULL) error2("oom pcodonClass");
   }
   /* --- */

   for(itree=0; ntree==-1||itree<ntree; itree++,iteration=1) {
      if(ReadTreeN(ftree,&haslength, &i,0,1))
            { puts("end of tree file."); break; }

      printf("\nTREE # %2d\n", itree+1);
      fprintf(fout,"\n\nTREE # %2d:  ", itree+1);
      fprintf(flnf,"\n\n%2d\n", itree+1);
      if(com.print) fprintf (frst,"\n\nTREE # %2d\n", itree+1);
      fprintf(frub,"\n\nTREE #%2d\n", itree+1);

      /* we use haslength=0 and com.fix_blength=0 */
      if (com.fix_blength==2 && !haslength) error2("no branch lengths in tree");
      if (com.fix_blength>0 && !haslength) com.fix_blength=0;
      if (times++==0 && com.fix_blength>0 && haslength) {
         if(com.clock) puts("\nBranch lengths in tree are ignored");
         else {
            if(com.fix_blength==2)
               puts("\nBranch lengths in tree are fixed.");
            else if(com.fix_blength==1)
               puts("\nBranch lengths in tree used as initials.");
            if(com.fix_blength==1) {
               FOR(i,tree.nnode)
                  if((x[nodes[i].ibranch]=nodes[i].branch)<0)
                     x[nodes[i].ibranch]=1e-5;
            }
         }
      }
      LASTROUND=0;
      if(com.cleandata)
         nchange = MPScore(com.space);
      if(com.ns<40) { OutTreeN(F0,0,0); printf("   MP score: %.0f",nchange); }
      OutTreeN(fout,0,0); fprintf(fout,"   MP score: %.0f",nchange);

      if(!com.clock && nodes[tree.root].nson<=2 && com.ns>2) {
         puts("\nThis is a rooted tree, without clock.  Check.");
         fputs("\nThis is a rooted tree.  Please check!",fout);
      }
      GetInitials(x, &i);

      np = com.np;
      if(noisy>=3 && np<100) matout(F0,x,1,np);
      if(i==-1) iteration = 0;
      if(np>NP || np-com.ntime>NP-NBRANCH) error2("raise NP");
      if(com.sspace < spaceming2(np)) {
         com.sspace = spaceming2(np);
         printf ("\nspace adjusted to %9u bytes\n",com.sspace);
         if((com.space=(double*)realloc(com.space,com.sspace))==NULL) {
            printf("\ntrying to get %d bytes for ming2", com.sspace);
            error2("oom space");
         }
      }
      printf("\nntime & nrate & np:%6d%6d%6d\n",com.ntime,com.nrate,com.np);

/*
      if(itree && !finitials)  for(i=0;i<np-com.ntime;i++) x[com.ntime+i] = xcom[i];
*/
      if(iteration && np) {
         SetxBound(np, xb);
         SetxInitials (np, x, xb); /* start within the feasible region */
      }
      PointconPnodes ();
/*
for(i=0; i<com.npatt; i++)
com.fpatt[i] /= (double)com.ls;
*/
      lnL = com.plfun (x,np);
      if(noisy) {
         printf("\nnp =%6d", np);
         printf("\nlnL0 = %12.6f\n",-lnL);
      }

      if(iteration && np) {
         if(com.method == 1)
            j = minB (noisy>2?frub:NULL, &lnL,x,xb, e, com.space);
         else if (com.method==3)
            j = minB2(noisy>2?frub:NULL, &lnL,x,xb, e, com.space);
         else
            j = ming2(noisy>2?frub:NULL,&lnL,com.plfun,NULL,x,xb, com.space,e,np);

         if (j==-1 || lnL<=0 || lnL>1e7) status=-1;
         else status=0;
         if(status) fprintf(fout,"\ncheck convergence..");

      }
      printf("Out..\nlnL  = %12.6f\n",-lnL);

      printf("%d lfun, %d eigenQcodon, %d P(t)\n",NFunCall, NEigenQ, NPMatUVRoot);
      if (itree==0)
         { lnL0=lnL;  FOR(i,np-com.ntime) xcom[i]=x[com.ntime+i]; }
      else if (!j)
         for (i=0; i<np-com.ntime; i++) xcom[i]=xcom[i]*.2+x[com.ntime+i]*0.8;

      if(!LASTROUND && (com.NSsites==NSpselection||com.NSsites==NSM2aRel||com.NSsites==NSdiscrete
        ||com.NSsites==NSfreqs||com.NSsites==NS3normal)) { /* Ignore */
         /* transform back to p0, p1,... */
         k=com.ntime+com.nrgene+com.nkappa+com.npi;

         if(com.nparK) {   /* HMM model for w */
            k += com.ncatG;
            for(i=0; i<com.ncatG; i++,k+=com.ncatG-1)
               f_and_x(x+k,x+k,com.ncatG,0,0);
         }
         else {
            j = (com.NSsites==NS3normal ? 3 : com.ncatG);
            if(com.model && com.model<=NSbranch2) j=3;
            f_and_x(x+k,x+k,j,0,0);
         }
      }
      LASTROUND=1;
      if(com.NSsites==NSdiscrete && com.aaDist==0 && com.model==0)
         sortwM3(x);
      if(com.clock) { /* move times into x[] */
         for(i=0,j=!nodes[tree.root].fossil; i<tree.nnode; i++)
            if(i!=tree.root && nodes[i].nson && !nodes[i].fossil)
               x[j++] = nodes[i].age;
      }

      fprintf (fout,"\nlnL(ntime:%3d  np:%3d): %13.6f %+14.6f\n",
         com.ntime, np, -lnL, -lnL+lnL0);

      if(com.fix_blength<2) {
         OutTreeB(fout);  FPN(fout);
      }
/*
      OutTreeB(fout);  FPN(fout);
      if(com.fix_blength==2) {
         for(i=0; i<tree.nbranch; i++) fprintf(fout, " %8.5f", nodes[tree.branches[i][1]].branch);
         FPN(fout);
      }
*/
      for(i=0; i<np; i++) fprintf(fout," %8.6f",x[i]);
      FPN(fout); fflush(fout);

      if (com.getSE) { /* Ignore */
         puts("Calculating SE's");
         if(com.sspace < np*(np+1)*sizeof(double)) {
            com.sspace = np*(np+1)*sizeof(double);
            if((com.space=(double*)realloc(com.space,com.sspace))==NULL)
               error2("oom space for SE");
         }

         g = com.space;
         H = g + com.np;
         HessianSKT2004 (x, lnL, g, H);
         if(com.getSE>=2 && com.clock==0 && nodes[tree.root].nson==3) {  /* g & H */
            fprintf(frst2,"\n %d\n\n", com.ns);
            OutTreeN(frst2, 1, 1);  fprintf(frst2,"\n\n");
            for(i=0; i<com.ntime; i++)
               if(x[i]>0.0004 && fabs(g[i])<0.005) g[i] = 0;
            for(i=0; i<com.ntime; i++) fprintf(frst2," %9.6f", x[i]);  fprintf(frst2, "\n\n");
            for(i=0; i<com.ntime; i++) fprintf(frst2," %9.6f", g[i]);  fprintf(frst2, "\n\n");
            fprintf(frst2, "\nHessian\n\n");
            for(i=0; i<com.ntime; i++,FPN(frst2))
               for(j=0; j<com.ntime; j++)
                  fprintf(frst2," %10.4g", H[i*np+j]);
            fflush(frst2);
         }

         for(i=0; i<np*np; i++)  H[i] *= -1;
         matinv(H, np, np, H+np*np);
         fprintf(fout,"SEs for parameters:\n");
         for(i=0; i<np; i++)
            fprintf(fout," %8.6f", (H[i*np+i]>0. ? sqrt(H[i*np+i]) : -1));
         FPN(fout);
      } /* if (com.getSE) */

      if(com.seqtype==1 && com.ntime && com.clock==0)
         fprintf(fout,"\nNote: Branch length is defined as number of nucleotide substitutions per codon (not per neucleotide site).\n");
      if(com.Mgene>1) {
         fprintf(fout,"Note: Branch length is defined for the first gene (site partition).\n");
         fprintf(fout,"For other genes, look at \"rates for genes\".\n");
      }

      /* if (com.clock) SetBranch (x); */
      if(com.clock && com.nbtype>1)
         fputs("\nWarning: branch rates are not yet applied in tree length and branch lengths",fout);
      if(AbsoluteRate)
         fputs("\nNote: mutation rate is not applied to tree length.  Tree has times, for TreeView",fout);
      for(i=0,tl=0; i<tree.nnode; i++)
         if(i!=tree.root) tl += nodes[i].branch;
      fprintf(fout,"\ntree length = %9.5f%s\n",tl,com.ngene>1?" (1st gene)":"");



#ifdef NSSITESBandits
      if(com.NSsites==0) {
         for(i=com.ntime; i<com.np; i++) fprintf(frst1,"\t%.3f", x[i]);
         fprintf(frst1,"\t%.2f\t%.3f", tl, -lnL);

         fM0tree=(FILE*)gfopen("M0tree", (insmodel==0?"w":"a"));
         fprintf(fM0tree, "%d  %d\n", com.ns, 1);
         OutTreeN(fM0tree,1,1);  FPN(fM0tree);
         fclose(fM0tree);
      }
      else {
         for(i=com.ntime; i<com.np; i++) fprintf(frst1,"\t%.3f",x[i]);
         fprintf(frst1,"\t%.3f",-lnL);
      }
#else

      for(i=0; i<com.np; i++) fprintf(frst1,"\t%.3f",x[i]);
      fprintf(frst1,"\t%.3f", -lnL);

/*
      fprintf(frst1,"\t%.4f", (com.ns==2 ? x[0]*2 : 0));
      for(i=0; i<com.nkappa; i++) fprintf(frst1,"\t%.3f",x[com.ntime+i]);
      fprintf(frst1,"\t%.4f", com.omega);
      fprintf(frst1,"\t%.3f", -lnL);
*/
#endif

      FPN(fout); OutTreeN(fout,0,1);  FPN(fout);
      FPN(fout); OutTreeN(fout,1,1);  FPN(fout);
      if(com.clock) {
         FPN(fout); OutTreeN(fout,1,PrNodeNum); FPN(fout);
      }

      if(com.np-com.ntime || com.clock)
         DetailOutput(fout,x, H);

      if (com.seqtype==AAseq && com.model>=REVaa_0)
         eigenQaa(fout, Root, U, V, x+com.ntime+com.nrgene);

      if (com.NSsites)
         lfunNSsites_rate(frst,x,np);
      if (com.print) {
         if(com.rho==0 && com.nparK==0 && com.clock<=1)
            AncestralSeqs(frst,x);
         if(!com.NSsites && com.plfun!=lfun)
            lfunRates(frate,x,np);
      }
      com.print -= 9;
      lnL = com.plfun(x,np);
      com.print += 9;

      fflush(fout);  fflush(flnf);  fflush(frst);  fflush(frst1);
   }     /* for(itree) */

   fclose(ftree);
   if(frate) fclose(frate);
   if (com.aaDist && com.aaDist<10 && com.aaDist!=AAClasses
      && (com.seqtype==CODONseq||com.model==FromCodon))
      printf("\n%s, %s.\n", com.daafile, (com.aaDist>0 ? "geometric" : "linear"));
   if(com.seqtype==1 && com.aaDist>=FIT1) free(pcodonClass);
   if(ntree==-1)  ntree=itree;

   if(ntree>1) {
      rewind(flnf);
      rell(flnf, fout, ntree);
   }

   fclose(flnf);
   return (0);
}
