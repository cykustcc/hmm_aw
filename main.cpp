#include "hmm.h"
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>

#include <gflags/gflags.h>
#include <glog/logging.h>

DEFINE_string(infilename, "",
              "input data file name");
DEFINE_string(mdfilename, "",
              "output model file name");
DEFINE_int32(dim, 10,
             "data dimension");
DEFINE_int32(num, 1,
             "Optional, number of sequences");
DEFINE_int32(statenum, 3,
             "umber of states in HMM");
DEFINE_int32(len, 100,
             "sequence length");
DEFINE_bool(forcediag, false,
             "sequence length");


int main(int argc, char *argv[])
{
  gflags::SetUsageMessage("train HMM with Gaussian (or GMMs) emission\n"
                          "usage: train <command> <args>\n\n"
                          "commands:\n"
                          "  infilename      input data file name\n"
                          "  mdfilename      output model file name\n"
                          "  dim             data dimension\n"
                          "  num             number of sequences\n"
                          "  statenum        number of states in HMM\n"
                          "  l               sequence length");
  gflags::ParseCommandLineFlags(&argc, &argv, true);
  FILE *infile, *mdfile;
  int i,j,k,m,n;
  int dim = FLAGS_dim;
  int nseq = FLAGS_num;
  int numdat, onelen = FLAGS_len;
  int numst = FLAGS_statenum;
  double lhsum;
  float epsilon=EPSILON;
  float tp1, tp2;

  std::string infilenamestr = FLAGS_infilename;
  std::string mdfilenamestr = FLAGS_mdfilename;
  std::string hmmmdfilenamestr = mdfilenamestr.substr(0,mdfilenamestr.size()-3) + "txt";
//  std::cout<<hmmmdfilenamestr<<std::endl;
  const char * infilename = infilenamestr.c_str();
  const char * mdfilename = mdfilenamestr.c_str();
  const char * hmmmdfilename = hmmmdfilenamestr.c_str();

  /*----------------------------------------------------------------*/
  /*--------------------- open files -------------------------------*/
  /*----------------------------------------------------------------*/
  infile = fopen(infilename, "r");
  if (infile == NULL) {
      printf("Couldn't open input data file \n");
      exit(1);
   }

  mdfile = fopen(mdfilename, "wb");
  if (mdfile == NULL) {
      printf("Couldn't open output model file \n");
      exit(1);
  }

  /*----------------------------------------------------------------*/
  /*----------------- Read in data ---------------------------------*/
  /*----------------------------------------------------------------*/

  // Assume the same length for all the sequences
  std::vector<int> len(nseq, 0);
  if (onelen>0) {
    for (i=0;i<nseq;i++) len[i]=onelen;
  }
  else {
    for (i=0;i<nseq;i++) {     //read in len from stdin
      fscanf(stdin, "%d", &len[i]);
    }
  }
  std::ifstream input(FLAGS_infilename, std::ios::in | std::ifstream::binary);
  std::vector<std::vector<float>> u;
  for (i=0,numdat=0;i<nseq;i++){
    numdat+=len[i];
    u.push_back(std::vector<float>(len[i] * dim, 0.0));
    std::cout<<len[i]<<std::endl;
  }

  // For testing purpose only
  // fprintf(stdout, "nseq=%d, numdat=%d, dim=%d\n",nseq,numdat,dim);
  // for (i=0;i<nseq;i++)
  //  fprintf(stdout, "len[%d]=%d\n", i,len[i]);

  std::vector<int> size(2, 0); // dimension, seq_len
	input.read(reinterpret_cast<char*>(&size[0]), size.size() * sizeof(int));
  printf("dimension=%d, n=%d\n", size[0], size[1]);
  for (i=0; i<nseq; i++) {
    input.read(reinterpret_cast<char*>(&u[i][0]), u[i].size() * sizeof(float));
    print_vector(u[i], len[i], dim);
  }

  /*----------------------------------------------------------------*/
  /*----------------- Estimate HMM  ---------------------------------*/
  /*----------------------------------------------------------------*/

  //fprintf(stderr, "numdat=%d, nseq=%d\n",numdat,nseq);

  std::vector<double> loglikehd(nseq,0.0);

  std::vector<int> stcls;
  std::vector<double> wt;
  HmmModel md(dim, numst, numst, stcls);
  md.hmmfit(u, nseq, len, loglikehd, lhsum,
	 (double)epsilon, wt, FLAGS_forcediag);

  //Output loglikehd from hmmfit() is not written out
  md.write_model("");
  md.print_model("");
  return 0;
}
