#include "nucaln.h"
#include <iostream>
#include <cstdint>

using namespace std;

void test(const string &s1, const string &s2);

int main(int argc, char* argv[]) {
   string topseq("CCCATCAAATG-CATCCTGGCTCTC");
   string bottomseq("CCCATC-AATGCCATCCTGGCTCTC");

   Nucaln stemaln(topseq, bottomseq);
   cout << stemaln.getConsensus().first << endl;

   test("CGAAATCCTCTATGTCAACATCTC", "CGAAATCCTCTATGTCAGCATCTC");
   test("GCCTG-AAGATTTACAACTG-AACTC", "--CTGAAAGATTTACAACTGAAACTC");
   test("CCACATCAAGGAGTATTTCTACTC", "CCACATCAAGGAGTATTTCTACTC");
   test("--CGACGAGAACTTGGCATTGAACTC", "AACGACGAGAACTTGGCATTGAACTC");
   int8_t x = -15;
   cout << "8 bits integer: " << x*1 << endl;
   cout << "8 bits integer: " << x*2 << endl;
   cout << "casted to int: " << int(x) << endl;

   return 0;
}

void test(const string &s1, const string &s2) {
   Nucaln stemaln(s1, s2);
   cout << "\nOriginal input:\n" << stemaln << endl;
   pair<string,string> result = stemaln.getConsensus();
   if (result.second.empty()) {
      cout << stemaln.getConsensus().first << endl;
   }
   else {
      cout << "two different sequences\n"
         << result.first << endl << result.second << endl;
   }
}
