#include "FastxParser.hpp"
#include <iostream>
#include <thread>
#include <vector>
#include <string>
#include <list>
#include <fstream>
#include <sstream>
#include <math.h>


#define MIN(x,y) ((x) < (y) ? (x) : (y)) //calculate minimum between two values
// Taken from https://www.tutorialspoint.com/cplusplus-program-to-implement-levenshtein-distance-computing-algorithm
int Levenstein_distance(std::string s1, std::string s2){
   int i,j,l1,l2,t,track;
   int dist[50][50];

   //stores the lenght of strings s1 and s2
   l1 = s1.length();
   l2= s2.length();
   for(i=0;i<=l1;i++) {
      dist[0][i] = i;
   }
   for(j=0;j<=l2;j++) {
      dist[j][0] = j;
   }
   for (j=1;j<=l1;j++) {
      for(i=1;i<=l2;i++) {
         if(s1[i-1] == s2[j-1]) {
            track= 0;
         } else {
            track = 1;
         }
         t = MIN((dist[i-1][j]+1),(dist[i][j-1]+1));
         dist[i][j] = MIN(t,(dist[i-1][j-1]+track));
      }
   }
   return dist[l2][l1];
}

struct Bases {
  uint32_t A, C, G, T;
};

int main(int argc, char* argv[]) {
   if (argc == 1) {
       std::cerr << "usage: error_correct_bc fastqR1 fastqR2 whitelist";
       return 1;
   } 
   int numFiles = argc - 1;
   if (numFiles != 3) {
       std::cerr << "you must provide 3 files, two fastq and one whitelist file!\n";
       return 1;
   }
   
   size_t numPairs = numFiles / 2;
   std::vector<std::string> files;
   std::vector<std::string> files2;
   std::vector<std::string> files3;

   files3.push_back(argv[3]);
   for (size_t i = 1; i <= numPairs; ++i) {
       files.push_back(argv[i]);
       files2.push_back(argv[i + numPairs]);
   }

   // The whitelist needs to be parsed and list of strings generated


   std::vector<std::string> barcodes;
   std::ifstream ifs;
   ifs.open (argv[3], std::ifstream::in);
   std::string line;

   while(std::getline(ifs, line)){

     std::stringstream    linestream(line);
     std::string          barcode;
     int                  v1;

     std::getline(linestream, barcode, '\t');

     linestream >> v1;

     barcodes.push_back(barcode);
     
   }

   ifs.close();

   //Now start the fastq parsing and check to see if the barcode is within levenstein distance 


  size_t nt = 4;
  size_t np = 2;
  fastx_parser::FastxParser<fastx_parser::ReadQualPair> parser(files, files2, nt, np);
  parser.start();

  std::vector<std::thread> readers;
  std::vector<Bases> counters(nt, {0, 0, 0, 0});

  std::ofstream outf;
  std::ofstream outf1;


  outf.open("output.fastq.1");
  outf1.open("output.fastq.2");

  for (size_t i = 0; i < nt; ++i) {
    readers.emplace_back([&, i]() {
      auto rg = parser.getReadGroup();
      while (true) {
        if (parser.refill(rg)) {
          for (auto& seqPair : rg) {
            auto& seq = seqPair.first;
            auto& seq2 = seqPair.second;

	      std::string c = seq.seq;

	      std::string bc_fastq = c.substr(0,24);

	      std::string umi = c.substr(24,16);


	      for(std::vector<int>::size_type  i = 0; i != barcodes.size(); ++i){
		if (Levenstein_distance(barcodes[i], bc_fastq) <= 2){

		  std::string bc_umi = bc_fastq + umi;
		    outf << '@' << seq.name  << '\n'<< bc_umi << '\n'<< '+' << '\n'<< seq.qual << '\n';
		    outf1 << '@'<< seq2.name  << '\n'<< seq2.seq << '\n'<< '+' << '\n'<< seq2.qual << '\n';

		  break;
		}else{
		}
		
	      }
	      // Need a function to parse whitelist and then generate a list of barcodes 
              // need to select barcode str.substr(0,24)
	      // Then see if the barcode is in the whitelist 




	      std::string c2 = seq2.seq;
	      //for second function if the read is in barcode then the second read is output, otherwise it is discarded
	      

          }
        } else {
          break;
        }
      }
    });
  }

  for (auto& t : readers) {
    t.join();
  }

  //parser.stop();
  outf.close();
  outf1.close();
  return 0;
}
