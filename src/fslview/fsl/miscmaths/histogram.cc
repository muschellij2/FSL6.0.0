/*  histogram.cc

    Mark Woolrich, Matthew Webster and Emma Robinson, FMRIB Image Analysis Group

    Copyright (C) 1999-2000 University of Oxford  */

/*  CCOPYRIGHT  */

#include "miscmaths.h"
#include "histogram.h"

using namespace std;

#ifndef NO_NAMESPACE
namespace MISCMATHS {
#endif

  float Histogram::getPercentile(float perc) 
  {
    if(histogram.Nrows()==0) generate();
      
    generateCDF();
    float percentile=getValue(1);
     
    for (int i=2;i<=CDF.Nrows();i++){
      if(CDF(i)>perc){
	double diff=(CDF(i)-perc)/(CDF(i)-CDF(i-1));
	percentile=getValue(i-1)+diff*(getValue(i)-getValue(i-1));
	break;
      }

    }
    return percentile;
  }

  void Histogram::generate()
  {
    Tracer ts("Histogram::generate");
  
    int size = sourceData.Nrows();
  
    if(calcRange)
      {
	// calculate range automatically
	histMin=histMax=sourceData(1);
	for(int i=1; i<=size; i++)
	  {
	    if (sourceData(i)>histMax)
	      histMax=sourceData(i);
	    if (sourceData(i)<histMin)
	      histMin=sourceData(i);
	  }
      }
    
  // zero histogram
    histogram.ReSize(bins);
    histogram=0;
    
    // create histogram; the MIN is so that the maximum value falls in the
    // last valid bin, not the (last+1) bin
    for(int i=1; i<=size; i++)
      {      
      histogram(getBin(sourceData(i)))++;
      }

    datapoints=size;
  }
  
  void Histogram::generate(ColumnVector exclude)
  {
    Tracer ts("Histogram::generate");
    int size = sourceData.Nrows();
    if(calcRange)
      {
	// calculate range automatically
	histMin=histMax=sourceData(1);
	for(int i=1; i<=size; i++)
	  {
	    if (sourceData(i)>histMax)
	      histMax=sourceData(i);
	    if (sourceData(i)<histMin)
	      histMin=sourceData(i);
	  }
      }
    
  
    histogram.ReSize(bins);
    histogram=0;
  
    datapoints=size;
  
    for(int i=1; i<=size; i++)
      {
	if(exclude(i)>0){
	  histogram(getBin(sourceData(i)))++;
	}else datapoints--;

      }
    

  }


  void Histogram::smooth()
    {
      Tracer ts("Histogram::smooth");

      ColumnVector newhist=histogram;

      // smooth in i direction
      newhist=0;
      ColumnVector kernel(3); 
      // corresponds to Gaussian with sigma=0.8 voxels
      //       kernel(1)=0.5;
      //       kernel(2)=0.2283;      
      //       kernel(3)=0.0219;
      // corresponds to Gaussian with sigma=0.6 voxels
      //       kernel(1)=0.6638;
      //       kernel(2)=0.1655;      
      //       kernel(3)=0.0026;

      //gauss(0.5,5,1)
      kernel(1)=0.7866;
      kernel(2)=0.1065;      
      kernel(3)=0.0003;

      for(int i=1; i<=bins; i++)
	  {
	    float val=0.5*histogram(i);
	    float norm=kernel(1);

	    if(i>1)
	      {
		val+=kernel(2)*(histogram(i-1));
		norm+=kernel(2);
	      }
	    if(i>2)
	      {
		val+=kernel(3)*(histogram(i-2));
		norm+=kernel(3);		
	      }
	    if(i<bins)
	      {
		val+=kernel(2)*(histogram(i+1));
		norm+=kernel(2);
	      }
	    if(i<bins-1)
	      {
		val+=kernel(3)*(histogram(i+2));
		norm+=kernel(3);		
	      }
	    val/=norm;

	    newhist(i)=val;
	  }

      histogram=newhist;

    }

  int Histogram::integrate(float value1, float value2) const
    {
      int upperLimit = getBin(value2);
      int sum = 0;

      for(int i = getBin(value1)+1; i< upperLimit; i++)
	{
	  sum += (int)histogram(i);
	}
      return sum;
    }

  float Histogram::mode() const
    {
      int maxbin = 0;
      int maxnum = 0;

      for(int i = 1; i< bins; i++)
	{
	  if((int)histogram(i) > maxnum) {
	    maxnum = (int)histogram(i);
	    maxbin = i;
	  }
	}

      return getValue(maxbin);
    }



  void Histogram::match(Histogram &target){
   
    int bin, newbin;
    double cdfval, val,dist;
    ColumnVector CDF_ref;
    ColumnVector olddata;
    ColumnVector histnew(bins);
    vector<double> rangein,rangetarg;

    CDF_ref=target.getCDF();
    olddata=sourceData;
    histnew=0; dist=0;
    if(exclusion.Nrows()==0){cout << " Histogram::match no excl " << endl; exclusion.ReSize(sourceData.Nrows()); exclusion=1;}
   
    float binwidthtarget=(target.getHistMax() - target.getHistMin())/target.getNumBins();
    
    
    for (int i=1;i<=sourceData.Nrows();i++){
     
      if(exclusion(i)>0){

	bin=getBin(sourceData(i));
      
	cdfval=CDF(bin);
	newbin=1;
	
	val=sourceData(i)-histMin;
      
	if(bin==target.getNumBins()){
	  newbin=bin;
	}else if( (cdfval- CDF_ref(bin)) > 1e-20){
	
	 
	  for (int j=bin;j<target.getNumBins();j++){
	    
	    if((cdfval - CDF_ref(j)) > 1e-20 && (cdfval - CDF_ref(j+1)) <= 1e-20 ){
	      newbin=j+1;
	      dist=(cdfval-CDF_ref(j))/(CDF_ref(j+1)-CDF_ref(j)); // find distance between current bins as proportion of bin width
	     break;
	    }
	  }
	  
	}
	else{
	  
	  //search backwards
	  
	  
	
	  int j=bin-1;
	  while (j>0){
	    
	    //  if(sourceData(i)<1.2 && sourceData(i)>1)
	    if((cdfval - CDF_ref(1)) < -1e-20 &&  fabs(CDF_ref(bin) - CDF_ref(1)) < 1e-20){newbin=j+1; break; }
	    else if ((cdfval - CDF_ref(1)) > -1e-20 &&  fabs(CDF_ref(bin) - CDF_ref(1)) < 1e-20){newbin=j+1; break; }
	    else if((cdfval - CDF_ref(j)) > 1e-20 && (cdfval - CDF_ref(j+1)) <= 1e-20 ){
	      newbin=j+1;

	      
	      dist=(cdfval-CDF_ref(j))/(CDF_ref(j+1)-CDF_ref(j)); // find distance between current bins as proportion of bin width
	      break;
	    }
	    j--;
	  }
	  //search forwards
	}
	
      

	sourceData(i)=target.getHistMin() + (newbin-1)*binwidthtarget+dist*binwidthtarget;

	histnew(newbin)++;
	
	if(sourceData(i) < target.getHistMin()) sourceData(i) = target.getHistMin();
	if(sourceData(i) > target.getHistMax()) sourceData(i) = target.getHistMax();

      }
    }
  }
  
#ifndef NO_NAMESPACE
  }
#endif




































