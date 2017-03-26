#include "hmm.h"
#include "utils.h"
#include <cstring>
int main(int argc, char *argv[])
{
    char infilename[100];
    char mdfilename[100];
    FILE *infile, *mdfile, *mdfile2;
    int i,j,k,m,n;
    int dim=2;
    float *dat;
    float **u;
    double *wt=NULL;
    int nseq, numdat, onelen=0, *len, *stcls;
    HmmModel *md=NULL;
    int numst=2;
    double *loglikehd, lhsum;
    float epsilon=EPSILON;
    float tp1, tp2;
    
    /*----------------------------------------------------------------*/
    /*---------------- Read in parameters from command line-----------*/
    /*----------------------------------------------------------------*/
    
    i = 1;
    while (i <argc)
    {
        if (*(argv[i]) != '-')
        {
            printf("**ERROR** bad arguments\n");
            exit(1);
        }
        else
        {
            switch(*(argv[i++] + 1))
            {
                case 'i':
                    strcpy(infilename,argv[i]);
                    break;
                case 'm':
                    strcpy(mdfilename,argv[i]);
                    break;
                case 'd':
                    sscanf(argv[i],"%d",&dim);
                    break;
                case 'n':
                    sscanf(argv[i],"%d",&nseq);
                    break;
                case 'l':
                    sscanf(argv[i],"%d",&onelen);
                    break;
                case 's':
                    sscanf(argv[i],"%d",&numst);
                    break;
                case 'e':
                    sscanf(argv[i],"%f",&epsilon);
                    break;
                default:
                {
                    printf("**ERROR** bad arguments\n");
                    exit(1);
                }
            }
            i++;
        }
    }
    
    /*----------------------------------------------------------------*/
    /*--------------------- open files -------------------------------*/
    /*----------------------------------------------------------------*/
    
    infile = fopen(infilename, "r");
    if (infile == NULL)
    {
        printf("Couldn't open input data file \n");
        exit(1);
    }
    
    
    mdfile2 = fopen(mdfilename, "rb");
    if (mdfile2 == NULL)
    {
        printf("Couldn't open output model file \n");
        exit(1);
    }
    
    /*----------------------------------------------------------------*/
    /*----------------- Read in data ---------------------------------*/
    /*----------------------------------------------------------------*/
    
    // Assume the same length for all the sequences
    len = (int *)calloc(nseq, sizeof(int));
    if (onelen>0) {
        for (i=0;i<nseq;i++) len[i]=onelen;
    }
    else {
        for (i=0;i<nseq;i++) {     //read in len from stdin
            fscanf(stdin, "%d", len+i);
        }
    }
    
    for (i=0,numdat=0;i<nseq;i++) { numdat+=len[i];}
    dat=(float *)calloc(numdat*dim,sizeof(float));
    u=(float **)calloc(nseq,sizeof(float *));
    for (i=0,m=0;i<nseq;i++) {
        u[i]=dat+m*dim;
        m+=len[i];
    }
    
    // For testing purpose only
    //fprintf(stdout, "nseq=%d, numdat=%d, dim=%d\n",nseq,numdat,dim);
    //for (i=0;i<nseq;i++)
    // fprintf(stdout, "len[%d]=%d\n", i,len[i]);r
    int size[2],dimension,seq_len;
    fread(size,sizeof(int),2, infile);
    dimension=size[0];
    seq_len=size[1];
    printf("dimension=%d, n=%d\n",dimension,seq_len);
    print_vector(dat, dim, numdat);
    fread(dat,sizeof(float),numdat*dim, infile);
    /*for (m=0;m<numdat;m++) {*/
    /*if (feof(nfile)) {*/
    /*fprintf(stderr, "Error: not enough data in input file\n");*/
    /*exit(0);*/
    /*}*/
    /*for (j=0;j<dim;j++) {*/
    /*fscanf(infile, "%e",&tp1);*/
    /*dat[m*dim+j]=tp1;*/
    /*}
     */
    /*fscanf(infile, "\n");*/
    /*}*/
    
    /*----------------------------------------------------------------*/
    /*----------------- Estimate HMM  ---------------------------------*/
    /*----------------------------------------------------------------*/
    
    //fprintf(stderr, "numdat=%d, nseq=%d\n",numdat,nseq);
    
    loglikehd=(double *)calloc(nseq,sizeof(double));
    md=(HmmModel *)calloc(1,sizeof(HmmModel));
    
    //Output loglikehd from hmmfit() is not written out
    
    // Binary file for the output model
//    write_model(md, mdfile);
    read_model(md, mdfile2);
    
    //Ascii file for the model
    md.print_model("");
    
    /*----------------------------------------------------------------*/
    /*------------ Compute the likelihood of the sequence  -----------*/
    /*----------------------------------------------------------------*/
    double loglikelihood,oneseq_likelihood;
    double *thetalog;
    thetalog=(double *)calloc(seq_len*md->numst, sizeof(double));
    for (i=0, loglikelihood=0.0; i<nseq; i++) {
        forward(u[i], len[i], thetalog,  md, &oneseq_likelihood);
        printf("loglikelihood for %d seq = %f\n", i, oneseq_likelihood);
        loglikelihood+=oneseq_likelihood;
    }
    printf("seq_len= %d, loglikelihood = %f\n",seq_len, loglikelihood);
    
}
