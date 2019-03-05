/*  histogram.h

    Mark Woolrich, Matthew Webster and Emma Robinson, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

/*  CCOPYRIGHT  */

#if !defined(_histogram_h)
#define _histogram_h

#include <iostream>
#include <fstream>
#define WANT_STREAM
#define WANT_MATH

#include "newmatap.h"
#include "newmatio.h"
#include "miscmaths.h"

using namespace NEWMAT;

namespace MISCMATHS {
 
  class Histogram
    {
    public:
      Histogram(){};
      const Histogram& operator=(const Histogram& in){
	sourceData=in.sourceData; calcRange=in.calcRange; histMin=in.histMin; histMax=in.histMax; bins=in.bins; histogram=in.histogram; CDF=in.CDF; datapoints=in.datapoints; exclusion=in.exclusion;
	return *this;
      }

      Histogram(const Histogram& in){ *this=in;}

      Histogram(const ColumnVector& psourceData, int numBins)
	: sourceData(psourceData), calcRange(true), bins(numBins){}

      Histogram(const ColumnVector& psourceData, float phistMin, float phistMax, int numBins) 
	: sourceData(psourceData), calcRange(false), histMin(phistMin), histMax(phistMax), bins(numBins){}
      
      void set(const ColumnVector& psourceData, int numBins) {	
	sourceData=psourceData; calcRange=true; bins=numBins;
      }

      void set(const ColumnVector& psourceData, float phistMin, float phistMax, int numBins) {	
	sourceData=psourceData; calcRange=false; histMin=phistMin; histMax=phistMax; bins=numBins;
      }

      void generate();
      void generate(ColumnVector exclusion_values);
      ColumnVector generateCDF();

      float getHistMin() const {return histMin;}
      float getHistMax() const {return histMax;}
      void setHistMax(float phistMax) {histMax = phistMax;}
      void setHistMin(float phistMin) {histMin = phistMin;}
      void setexclusion(ColumnVector exclusion_values) {exclusion =exclusion_values;}
      void smooth();

      int integrateAll() {return sourceData.Nrows();}

      const ColumnVector& getData() {return histogram;}
      void setData(const ColumnVector& phist) { histogram=phist;}

      int integrateToInf(float value) const { return integrate(value, histMax); }
      int integrateFromInf(float value) const { return integrate(histMin, value); }
      int integrate(float value1, float value2) const;
    
      void match(Histogram &);
      
      float mode() const;

      int getBin(float value) const;
      float getValue(int bin) const;
      float getPercentile(float perc);

      inline int getNumBins() const {return bins;}
      inline ColumnVector getCDF() const {return CDF;}
      inline ColumnVector getsourceData()const {return sourceData;}
    protected:

    private:

      ColumnVector sourceData;
      ColumnVector histogram;
      ColumnVector exclusion;
      ColumnVector CDF;

      bool calcRange;

      float histMin;
      float histMax;

      int bins; // number of bins in histogram
      int datapoints;
    };

  inline int Histogram::getBin(float value) const
    {
      float binwidth=(histMax-histMin)/bins;
      return Max(1, Min((int)((((float)bins)*((float)(value-(histMin-binwidth))))/((float)(histMax-histMin))),bins));
    }
  
  inline float Histogram::getValue(int bin) const
    {
      return (bin*(histMax-histMin))/(float)bins + histMin;
    }

  inline ColumnVector Histogram::generateCDF() 
    { 
     
     
      CDF.ReSize(bins);
      
      
      CDF(1)=histogram(1)/datapoints;
  
      for (int i=2;i<=bins;i++)
	CDF(i)=CDF(i-1)+ histogram(i)/datapoints;

     
      return CDF;
    }
}

#endif





