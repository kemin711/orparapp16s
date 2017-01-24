#include <iostream>
#include <fstream>
#include <cmath>

#include "stalactitegraph.h"

double StalactiteGraph::borderFactor = 1.07;
int StalactiteGraph::fontSize=28;

void StalactiteGraph::initParam() {
   if (fileType == PNG) {
      string drawDim=to_string(amplifyBorder(plotWidth)) + "x" + to_string(amplifyBorder(plotHeight));
      plpa.setplparam("BITMAPSIZE", (char*)drawDim.c_str());
   }
   else if (fileType == SVG) {
      plpa.setplparam("PAGESIZE", (char*)"letter");
   }
   // the rest is not implemented
}

void StalactiteGraph::plot(const vector<tuple<int,int,int> > &points, const string &cln)
{
   vector<pair<int, int> > p1, p2, p3;
   reducePoint(points, p1, p2, p3);
   int depth = getMaxDepth(p3);
   int width=points.size();
   if (depth < 5) {
      cerr << "cluster depth too low to draw a graph\n";
      return;
   }
   double hr = plotWidth/(double)width; // for transformation
   double vr = plotHeight/(double)depth;
   // set up parameter, done in initializer
   // only two file types now
   string outfile = cln + ".png";
   if (fileType == SVG)
      outfile = cln + ".svg";
   if (!outputdir.empty())
      outfile = outputdir + "/" + outfile;
   ofstream ouf(outfile);
   Plotter *plt = new PNGPlotter(cin, ouf, cerr, plpa);
   if (fileType == SVG) 
      plt = new SVGPlotter(cin, ouf, cerr, plpa);
   if (plt->openpl() < 0) {
      throw runtime_error("Failed to open plotter");
   }
   plt->erase();
   plt->fspace(-(0.1*plotWidth),-(0.1*plotHeight), amplifyBorder(plotWidth), amplifyBorder(plotHeight));
   plt->flinewidth(4);
   plt->pencolorname("black");
   // draw coordinates
   plt->fbox(0, 0, plotWidth+5, plotHeight*1.03);
   plt->flinewidth(1);
   // use plotfont -T png --help-fonts to get a list
   plt->fontname("HersheySans");
   int true_size = plt->fontsize(fontSize);
   plt->fmove(0, -true_size);
   plt->alabel('c', 'c', (char*)"0");
   int xx = tickUnit(width, 7);
   //cout << " xx " << xx << " hr,vr: " << hr << "," << vr 
   //   << " depth: " << depth << endl;
   for (int x=xx; x < width; x += xx) {
      //cerr << x << " ";
      plt->fline(x*hr, true_size, x*hr, 0);
      plt->endpath();
      plt->fmove(x*hr, -true_size);
      plt->alabel('c', 'c', to_string(x).c_str());
   }
   //cerr << "\nx label done\n";
   int yy = tickUnit(depth, 5);
   plt->fmove(-10,0);
   plt->alabel('r', 'c', (char*)"0");
   for (int y=yy; y < depth; y += yy) {
      //cerr << y << " ";
      plt->fline(0, y*vr, true_size, y*vr);
      plt->endpath();
      plt->move(-10, y*vr);
      plt->alabel('r', 'c', to_string(y).c_str());
   }
   //cerr << "\ny label done\n";
   // draw curves
   plt->pencolorname("red");
   plt->flinewidth(1);
   plt->fmove(p1[0].first*hr, p1[0].second*vr);
   size_t i;
   for (i=1; i<p1.size(); ++i) {
      plt->fcont(p1[i].first*hr, p1[i].second*vr);
   }
   plt->endpath();
   plt->fline(1540, 100, 1500, 100);
   plt->fmove(1495, 100);
   plt->alabel('r', 'c', (char*)"Consensus");
   plt->pencolorname("green");
   plt->flinewidth(0.5);
   plt->fmove(p2[0].first*hr, p2[0].second*vr);
   for (i=1; i<p2.size(); ++i) {
      plt->fcont(p2[i].first*hr, p2[i].second*vr);
   }
   plt->endpath();
   plt->fline(1540, 130, 1500, 130);
   plt->fmove(1495, 130);
   plt->alabel('r', 'c', (char*)"Base");
   plt->pencolorname("blue");
   plt->linewidth(0.25);
   plt->fmove(p3[0].first*hr, p3[0].second*vr);
   for (i=1; i<p3.size(); ++i) {
      plt->fcont(p3[i].first*hr, p3[i].second*vr);
   }
   plt->endpath();
   plt->fline(1540, 160, 1500, 160);
   plt->fmove(1495, 160);
   plt->alabel('r', 'c', (char*)"Total");
   plt->closepl();
   cerr << "graphics written to " << outfile << endl;
   delete plt;
}

void StalactiteGraph::reducePoint(const vector<tuple<int,int,int> > &data3,
      vector<pair<int,int> > &p1, vector<pair<int,int> > &p2, 
      vector<pair<int,int> > &p3) 
{
   //vector<pair<int, int> > p1, p2, p3;
   size_t i=0;
   int x=0;
   int y=get<0>(data3[0]);
   p1.push_back(make_pair(1, get<0>(data3[0])));
   p2.push_back(make_pair(1, get<1>(data3[0])));
   p3.push_back(make_pair(1, get<2>(data3[0])));
   ++i;
   // x,y is the last point
   while (i<data3.size()) {
      if (get<0>(data3[i]) != y) {
         if (i-1 > x) {
            p1.push_back(make_pair(i, y));
            x=i;
         }
         y=get<0>(data3[i]);
         p1.push_back(make_pair(i+1, y));
      }
      ++i;
   }
   // fill p2
   x=0; i=1;
   y=get<1>(data3[0]);
   while (i<data3.size()) {
      if (get<1>(data3[i]) != y) {
         if (i-1 > x) {
            p2.push_back(make_pair(i, y));
            x=i;
         }
         y=get<1>(data3[i]);
         p2.push_back(make_pair(i+1, y));
      }
      ++i;
   }
   // fill p3
   x=0; i=1;
   y=get<2>(data3[0]);
   while (i<data3.size()) {
      if (get<2>(data3[i]) != y) {
         if (i-1 > x) {
            p3.push_back(make_pair(i, y));
            x=i;
         }
         y=get<2>(data3[i]);
         p3.push_back(make_pair(i+1, y));
      }
      ++i;
   }
}

int StalactiteGraph::getMaxDepth(const vector<pair<int,int> > &pp) {
   int maxdepth = pp[0].second;
   size_t i=1;
   while (i<pp.size()) {
      if (maxdepth < pp[i].second)
         maxdepth = pp[i].second;
      ++i;
   }
   cerr << "max depth: " << maxdepth << endl;
   return maxdepth;
}

// numd: number of divisions
// len total length of the axis
int StalactiteGraph::tickUnit(int len, int numd) {
   // first try to round to 100
   int u = round(float(len)/numd);
   if (u > 100) {
      u = u/100*100;
   }
   else if (u>10) {
      u = u/10*10;
   }
   return u;
}


