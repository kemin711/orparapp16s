#ifndef STALACTITEGRAPH_H
#define STALACTITEGRAPH_H

// (C) 2016 Kemin Zhou at Roche
// Added to the consensus library for plotting graphics output
// This is a simple wrapper for the gnu plotutils library

#include <plotter.h>
#include <string>
#include <tuple>
#include <vector>

using namespace std;

/**
 * Needs to be consilidated with other ploting programs
 */
enum PlotFileType {
   PNG, SVG
};

/**
 * Draw base vertical coverage plot for each cluster.
 */
class StalactiteGraph {
   public:
      /**
       * Default file type PNG
       */
      StalactiteGraph() : fileType(PNG), plotWidth(1800), plotHeight(1200), outputdir() { initParam(); }
      /**
       * Default width 1800, height 1200
       */
      StalactiteGraph(PlotFileType ft) : fileType(ft), plotWidth(1800), plotHeight(1200), 
         outputdir() { initParam(); }
      StalactiteGraph(PlotFileType ft, int width, int height) : fileType(ft), plotWidth(1800), 
            plotHeight(1200), outputdir() { initParam(); }
      /**
       * Can do repeated work.
       * @param points data points from the consensus algorithm.
       */
      void plot(const vector<tuple<int,int,int> > &points, const string &cln);
      void setOutputDir(const string &dirname) { outputdir=dirname; }

   private:
      static int tickUnit(int len, int numd);
      static int getMaxDepth(const vector<pair<int,int> > &pp);
      /**
       * Convert the input into three vectors.
       * This will merge identical data points. In case of SVG, the
       * number of objects will be greatly reduced. 
       */
      static void reducePoint(const vector<tuple<int,int,int> > &data3,
                  vector<pair<int,int> > &p1, vector<pair<int,int> > &p2, 
                        vector<pair<int,int> > &p3);
      static int amplifyBorder(int original) { return int(double(original)*borderFactor); }
      void initParam();
      PlotFileType fileType;
      /**
       * Plotting area width default 1800
       */
      int plotWidth;
      /**
       * Plotting area height default 1200
       */
      int plotHeight;
      PlotterParams plpa;
      /**
       * Output directory default empty
       */
      string outputdir;

      /**
       * default 1.07
       */
      static double borderFactor;
      static int fontSize;
};

#endif
