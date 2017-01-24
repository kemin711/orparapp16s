#include "stalactitegraph.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <tuple>
#include <strformat.h>

using namespace std;

/**
 * Read the data structure for testing just get the first one
 */
int processAll(const string &file);
void testdrawCompound();
void testdrawC();
void testdrawString();

int main(int argc, char* argv[]) {
   string infile="BEIeven_vertcov.txt";
   int i=1;
   while (i<argc) {
      infile = argv[i];
      ++i;
   }
   int numcl = processAll(infile);
   cout << numcl << " clusters processed\n";

   return 0;
}

// for simulating the original data type
int processAll(const string &file) {
   StalactiteGraph stalac(PNG);

   ifstream ifs(file);
   string line, clname;
   getline(ifs, clname); // cluster name
   getline(ifs, line); // header 
   getline(ifs, line); // first data line
   int cnt=0;
   while (!ifs.eof()) {
      cerr << "\nprocessing cluster: " << clname << endl;
      vector<tuple<int, int, int> > res;
      while (!ifs.eof() && line.substr(0,2) != "cl") {
         vector<string> row = split(line, '\t');
         res.push_back(tuple<int,int,int>(stoi(row[1]), stoi(row[2]), stoi(row[3])));
         getline(ifs, line);
      }
      stalac.plot(res, clname);
      ++cnt;
      if (!ifs.eof()) {
         clname = line;
         getline(ifs, line);
         getline(ifs, line);
      }
   }
   return cnt;
}

void draw_c_curve(Plotter& plotter, double dx, double dy, int order) {
   static const float ratio=0.51;
   if (order >= 18)
      plotter.fcontrel(dx,dy);
   else {
      draw_c_curve(plotter, ratio*(dx-dy), ratio*(dx+dy), order+1);
      draw_c_curve(plotter, ratio*(dy+dx), ratio*(dy-dx), order+1);
   }
}

void testdrawC() {
   PlotterParams params;
   params.setplparam("PAGESIZE", (char*)"letter");
   ofstream ouf("testC.png");
   PNGPlotter plt(cin, ouf, cerr, params);
   if (plt.openpl() < 0) {
      throw runtime_error("Failed to open plotter");
   }
   plt.fspace(0,0, 1000.0, 1300.0);
   plt.flinewidth(0.25);
   plt.pencolorname("red");
   plt.erase();
   plt.fmove(600.0, 300.0);
   draw_c_curve(plt, 0.0, 400.0, 0);
   plt.closepl();
   cerr << "graphics in test.png\n";
}
      
void testdrawCompound() {
   PlotterParams params;
   params.setplparam("PAGESIZE", (char*)"letter");
   ofstream ouf("testcompound.png");
   PNGPlotter plt(cin, ouf, cerr, params);
   if (plt.openpl() < 0) {
      throw runtime_error("Failed to open plotter");
   }
   plt.fspace(0,0, 1000.0, 1000.0);
   plt.flinewidth(5);
   plt.pencolorname("green");
   plt.filltype(1);
   plt.erase();
   plt.orientation(1);
   plt.fbox(50.0, 50.0, 950.0, 950.0);
   plt.endsubpath();
   size_t i,j;
   for (i=0; i<2; ++i) {
      for (j=9; j>=3; --j) {
         plt.orientation(j%2? -1 : 1);
         plt.fcircle(250.0 + 500*i, 500.0, j*20.0);
         plt.endsubpath();
      }
      plt.fmove(225.0+500*i, 475.0);
      plt.fcont(250.0+500*i, 525.0);
      plt.fcont(275.0+500*i, 475.0);
      plt.endsubpath();
   }
   plt.endpath();
   plt.closepl();
   cerr << "graphics in test.png\n";
}
      
void drawBoxedString(Plotter &plotter, const char *s, double size, double angle) {
   double true_size, width;
   static const double expand=1.9;
   plotter.ftextangle(angle);
   true_size = plotter.ffontsize(size);
   width=plotter.flabelwidth(s);
   cerr << "size truesize: " << size << " " << true_size << " " << width << endl;
   plotter.fellipserel(0,0,expand*0.5*width, expand*0.5*true_size, angle);
   plotter.alabel('c', 'c', s);
}

void testdrawString() {
   PlotterParams params;
   double SIZE=150.0;
   params.setplparam("PAGESIZE", (char*)"letter");
   ofstream ouf("testString.png");
   PNGPlotter plt(cin, ouf, cerr, params);
   if (plt.openpl() < 0) {
      throw runtime_error("Failed to open plotter");
   }
   plt.fspace(-SIZE, -SIZE, SIZE, SIZE);
   plt.pencolorname("blue");
   plt.fillcolorname("white");
   plt.filltype(1);
   plt.fontname("NewCenturySchlbk-Roman");
   plt.erase();
   const char* label="GNU libplot!";
   for (size_t i=205; i>1; --i) {
      string tmp = itos(i) + " " + label;
      double theta, radius;
      theta=0.5*i;
      radius = SIZE/pow(theta, 0.35);
      plt.fmove(radius*cos(theta), radius*sin(theta));
      drawBoxedString(plt, tmp.c_str(), 0.04*radius, (180.0*theta/M_PI)-90.0);
   }
   if (plt.closepl() < 0) {
      throw runtime_error("failed to close");
   }
}

