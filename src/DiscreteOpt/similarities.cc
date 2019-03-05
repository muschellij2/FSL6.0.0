/*  similarities.cc

    Emma Robinson, FMRIB Image Analysis Group

    Copyright (C) 2012 University of Oxford  */

/*  Part of FSL - FMRIB's Software Library
    http://www.fmrib.ox.ac.uk/fsl
    fsl@fmrib.ox.ac.uk
    
    Developed at FMRIB (Oxford Centre for Functional Magnetic Resonance
    Imaging of the Brain), Department of Clinical Neurology, Oxford
    University, Oxford, UK
    
    
    LICENCE
    
    FMRIB Software Library, Release 5.0 (c) 2012, The University of
    Oxford (the "Software")
    
    The Software remains the property of the University of Oxford ("the
    University").
    
    The Software is distributed "AS IS" under this Licence solely for
    non-commercial use in the hope that it will be useful, but in order
    that the University as a charitable foundation protects its assets for
    the benefit of its educational and research purposes, the University
    makes clear that no condition is made or to be implied, nor is any
    warranty given or to be implied, as to the accuracy of the Software,
    or that it will be suitable for any particular purpose or for use
    under any specific conditions. Furthermore, the University disclaims
    all responsibility for the use which is made of the Software. It
    further disclaims any liability for the outcomes arising from using
    the Software.
    
    The Licensee agrees to indemnify the University and hold the
    University harmless from and against any and all claims, damages and
    liabilities asserted by third parties (including claims for
    negligence) which arise directly or indirectly from the use of the
    Software or the sale of any products based on the Software.
    
    No part of the Software may be reproduced, modified, transmitted or
    transferred in any form or by any means, electronic or mechanical,
    without the express permission of the University. The permission of
    the University is not required if the said reproduction, modification,
    transmission or transference is done without financial return, the
    conditions of this Licence are imposed upon the receiver of the
    product, and all original and amended source code is included in any
    transmitted product. You may be held legally responsible for any
    copyright infringement that is caused or encouraged by your failure to
    abide by these terms and conditions.
    
    You are not permitted under this Licence to use this Software
    commercially. Use for which any financial return is received shall be
    defined as commercial use, and includes (1) integration of all or part
    of the source code or the Software into a product for sale or license
    by or on behalf of Licensee to third parties or (2) use of the
    Software or any derivative of it for research with the final aim of
    developing software products for sale or license to a third party or
    (3) use of the Software or any derivative of it for research with the
    final aim of developing non-software products for sale or license to a
    third party, or (4) use of the Software to provide any service to an
    external organisation for which payment is received. If you are
    interested in using the Software commercially, please contact Oxford
    University Innovation ("OUI"), the technology transfer company of the
    University, to negotiate a licence. Contact details are:
    Innovation@innovation.ox.ac.uk quoting reference DE/9564. */
#include "similarities.h"
/// HISTOGRAM 2D FUNCTIONS
namespace DISCRETEOPT{ 

histogram2D::~histogram2D()
{
  m_nbinsx=0; 
  m_nbinsy=0;
  m_min_x   = 0;
  m_min_y   = 0;
  m_max_x   = 0;
  m_max_y   = 0;
  m_width_x = 0;
  m_width_y = 0;
  
  _nsamps=0;
}

void histogram2D::Initialize(int nbinsx, int nbinsy, double maxx, double maxy, double minx, double miny)
{
  
  m_nbinsx  = nbinsx; 
  m_nbinsy  = nbinsy;
  m_min_x   = minx;
  m_min_y   = miny;
  m_max_x   = maxx;
  m_max_y   = maxy;
  m_width_x = (m_max_x - m_min_x) / (double)m_nbinsx;
  m_width_y = (m_max_y - m_min_y) / (double)m_nbinsy;
 
  Matrix newmat(m_nbinsx, m_nbinsy);
 
  _bins=newmat;
  _weights=newmat;
  _bins=0;
  _nsamps=0;
  _weights=0;
}

void histogram2D::Zero()
{

  int i,j;

  for (j = 1; j <= m_nbinsy; j++) {
    for (i = 1; i <= m_nbinsx; i++) {
      _bins(i,j) = 0;
    }
  }
  _nsamps=0;
}


void histogram2D::AddSample(double x, double y)
{
  int i, j;
  

  if (x < m_min_x) return;
  if (x > m_max_x) return;
  if (y < m_min_y) return;
  if (y > m_max_y) return;

  i = (int) MISCMATHS::round(m_nbinsx * (x - m_min_x - 0.5*m_width_x) / (m_max_x - m_min_x)); 
  j = (int) MISCMATHS::round(m_nbinsy * (y - m_min_y - 0.5*m_width_y) / (m_max_y - m_min_y));
  
  if (i < 1) i = 1;
  if (j < 1) j = 1;
  if (i > m_nbinsx) i = m_nbinsx;
  if (j > m_nbinsy) j = m_nbinsy ;

 
  _bins(i,j) += 1;
  _weights(i,j)=1;
  _nsamps      += 1;
  
 
}

void histogram2D::AddSample(double weight, double x, double y)
{
  int i, j;
  
  if (x < m_min_x) return;
  if (x > m_max_x) return;
  if (y < m_min_y) return;
  if (y > m_max_y) return;
  
  i = (int) MISCMATHS::round(m_nbinsx * (x - m_min_x - 0.5*m_width_x) / (m_max_x - m_min_x)); 
  j = (int) MISCMATHS::round(m_nbinsy * (y - m_min_y - 0.5*m_width_y) / (m_max_y - m_min_y));
  
  if (i < 1) i = 1;
  if (j < 1) j = 1;
  if (i > m_nbinsx) i = m_nbinsx;
  if (j > m_nbinsy) j = m_nbinsy ;
  
  _bins(i,j) += 1;
  _weights(i,j)=1;// set as one for the time being until we can find correct formula weight;
  _nsamps      += 1;
  
  
}

void histogram2D::DeleteSample(double x, double y)
{
  int i, j;
  
 
  if (x < m_min_x) return;
  if (x > m_max_x) return;
  if (y < m_min_y) return;
  if (y > m_max_y) return;

  i = (int) MISCMATHS::round(m_nbinsx * (x - m_min_x - 0.5*m_width_x) / (m_max_x - m_min_x));
  j = (int) MISCMATHS::round(m_nbinsy * (y - m_min_y - 0.5*m_width_y) / (m_max_y - m_min_y));
  
  if (i < 0) i = 1;
  if (j < 0) j = 1;
  if (i > m_nbinsx) i = m_nbinsx;
  if (j > m_nbinsy) j = m_nbinsy;

  _bins(i,j) -= 1;
  _nsamps     -= 1;
 
}

double histogram2D::MarginalEntropyX(){
  int i,j;
  ColumnVector M(m_nbinsx);
  double E=0;
 
  M=0;


  for (i=1;i<=m_nbinsx;i++){
    for(j=1;j<=m_nbinsy;j++){
      M(i)=M(i)+_bins(i,j)/(double)_nsamps;
       }

  }

  for (i=1;i<=m_nbinsx;i++){
    if(M(i)>0)
      E-=M(i)*log(M(i));
  }
  
  return E;
}

double histogram2D::MarginalEntropyY(){
  int i,j;
  ColumnVector M(m_nbinsy);
  double E=0;

  M=0;

  for (j=1;j<=m_nbinsy;j++){
    for(i=1;i<=m_nbinsx;i++){
      M(j)=M(j)+_bins(i,j)/(double)_nsamps;
     
    }
  }

  for (i=1;i<=m_nbinsx;i++){
    if(M(i)>0)
    E-=M(i)*log(M(i));
  }
  return E;
}

double histogram2D::JointEntropy(){
  int i,j;
  double E=0;
  Matrix JP(m_nbinsx,m_nbinsy);

  for (i=1;i<=m_nbinsy;i++){
    for(j=1;j<=m_nbinsx;j++){
      JP(i,j)=_bins(i,j)/(double)_nsamps;
    }
  }

  for (i=1;i<=m_nbinsy;i++){
    for(j=1;j<=m_nbinsx;j++){
      if(JP(i,j)>0)
	E-=_weights(i,j)*JP(i,j)*log(JP(i,j));

   
    }
  }
 
  return E;
}




double histogram2D::normalisedmutualinformation(){
  if(this->JointEntropy() >1e-5 )return (this->MarginalEntropyX()+this->MarginalEntropyY())/this->JointEntropy();
  else return 0;
}



  void simkernel::set_input(string path_to_input){ // assume matrix is full??

    m_A =boost::shared_ptr<BFMatrix > (new FullBFMatrix (path_to_input));
  
  }

  void simkernel::set_reference(string path_to_ref){
  m_B =boost::shared_ptr<BFMatrix > (new FullBFMatrix (path_to_ref));
  }
  
  void simkernel::Initialize(int simval){
  
    _sim=simval;
    
    if(m_A==NULL){    throw  DISCRETEOPTHOCRException("SIMILARITIES:: Connectivity matrices have not been initliased");}
    
    if(_sim==1|| _sim==2){
      
      _rmeanA=meanvector(*m_A);
      if(m_A->Nrows()==1) _meanA=_rmeanA(1);
      if(m_B!=NULL)     
	_rmeanB=meanvector(*m_B);
      else
	_rmeanB=_rmeanA;
      
      if(m_B->Nrows()==1) _meanB=_rmeanB(1);
    }
    
    
    /// for NMI initialise histogram
    if(_sim==3){
      calc_range(m_A,maxx,minx);
      if(m_B==NULL){
	miny=minx;
	maxy=maxx;
      }
      else
	calc_range(m_B,maxy,miny);
      
      if(m_A->Nrows()<256)    
	hist.Initialize(m_A->Nrows()-1,m_A->Nrows()-1,maxx,maxy,minx,miny);
      else{
	hist.Initialize(256,256,maxx,maxy,minx,miny);
      }
    }
    
  }
  
  ///// for the case where we are working on vector data and we have no knowledge of the full data m_A (currently used in discrete opt)
  void simkernel::Initialize(int simval,const vector<double> &inputdata, const vector<double> &refdata,const vector<double> &weights){
  
    _sim=simval;
    _initialised=true;
    if(_sim==3){
      
      calc_range(inputdata,maxx,minx);
      calc_range(refdata,maxy,miny);
      if(inputdata.size()<256)    
	hist.Initialize(inputdata.size()-1,inputdata.size()-1,maxx,maxy,minx,miny);
      else{
	hist.Initialize(256,256,maxx,maxy,minx,miny);
      }
      
    }
    else{
    _meanA=0;_meanB=0;
    double sum=0;
    for(unsigned int i=0;i<inputdata.size();i++){
      double varwght=0;
      if(weights.size()==inputdata.size()) varwght=weights[i];
      else varwght=1;
     
      _meanA+=varwght*inputdata[i];
      _meanB+=varwght*refdata[i];
      sum+=varwght;
    }
    if(sum>0){
      _meanA/=sum;
      _meanB/=sum;
    }
    }
    

  }

  RowVector simkernel::meanvector(const BFMatrix  &fdt_matrix)
  {  
    RowVector mean(fdt_matrix.Ncols());
    mean=0;
    
    if(fdt_matrix.Nrows()==1){
      double sum=0;
      for (unsigned int i=0; i < fdt_matrix.Ncols(); i++){
	sum=sum+fdt_matrix.Peek(1,i+1);
      }
      for (unsigned int i=0; i < fdt_matrix.Ncols(); i++)
	mean(i+1)=sum/fdt_matrix.Ncols();
    }
    else{
      for (unsigned int i=0; i < fdt_matrix.Ncols(); i++){
	double sum=0;
	
	for (unsigned int j=0; j < fdt_matrix.Nrows(); j++){
	  sum=sum+fdt_matrix.Peek(j+1,i+1);
	  
	}
	mean(i+1)=sum/fdt_matrix.Nrows();
	
      }
    }
    
    return mean;
  }

  

  void simkernel::calc_range(boost::shared_ptr<BFMatrix > mat, double &max, double &min){


    for (unsigned int i=1; i <= mat->Nrows(); i++){
      for (unsigned int j=1; j <= mat->Ncols(); j++){
	if( mat->Peek(i,j)>max)
	  max=mat->Peek(i,j);
	if(mat->Peek(i,j)< min)
	  min=mat->Peek(i,j);
      }
    }
  }

  void simkernel::calc_range(const vector<double> &mat, double &max, double &min){

    max=0; min=1e7;
    for (unsigned int i=0; i < mat.size(); i++){
	if( mat[i]>max)
	  max=mat[i];
	if(mat[i]< min)
	  min=mat[i];
      }
    
  }

  

  double simkernel::corr(int i,int j){
    
    int num=0,numB=0;;
    double prod,varA,varB;
    prod=0; varA=0; varB=0;

    if(m_B.get()==0)   m_B=m_A;
  

    if(_rmeanA.Ncols()==0){
      _rmeanA=meanvector(*m_A);
      if(_rmeanB.Ncols()==0) {
	if(m_B==NULL)
	  _rmeanB=_rmeanA;
	else
	  _rmeanB=meanvector(*m_B);
      }
    }
    boost::shared_ptr<SparseBFMatrix<double> > ptr =boost::dynamic_pointer_cast<SparseBFMatrix<double> >(m_A);

    double Bzerooffset=(0.0 - _rmeanB(j));  // result for all zero values of sparse mat
    double Azerooffset=(0.0 - _rmeanA(i));
      
     for (BFMatrixColumnIterator it=m_A->begin(i);it!=m_A->end(i);it++){
      prod +=(*it-_rmeanA(i))*(m_B->Peek(it.Row(),j)-_rmeanB(j));
      varA+=(*it-_rmeanA(i))*(*it-_rmeanA(i));
      varB+=(m_B->Peek(it.Row(),j)-_rmeanB(j))*(m_B->Peek(it.Row(),j)-_rmeanB(j));
      num++;
    }


     if(ptr){ //if sparse run  all for all rows where A had zero values

      varA+=(m_A->Nrows()-num)*Azerooffset*Azerooffset; /// for rows where A has no values

      for (BFMatrixColumnIterator it=m_B->begin(j);it!=m_B->end(j);it++){
	if(m_A->Peek(it.Row(),i)==0){
	  prod +=Azerooffset*(*it-_rmeanB(j));
	  num++;
	  varB+=(*it-_rmeanB(j))*(*it-_rmeanB(j));
	}
	numB++;
      } 
    
      
      varB+=(m_A->Nrows()-numB)*Bzerooffset*Bzerooffset;    /// for rows where B has no values
      prod+=(2*m_A->Nrows()-num)*Azerooffset*Bzerooffset;  /// for rows where A&B have no values
    }  

    if (varA == 0.0 || varB == 0.0) return  0.0; else
      return prod/(sqrt(varA)*sqrt(varB));
    
  }


  double simkernel::SSD(int i, int j){
    
    double prod=0;
    if(m_B.get()==0)   m_B=m_A;
    
    boost::shared_ptr<SparseBFMatrix<double> > ptr =boost::dynamic_pointer_cast<SparseBFMatrix<double> >(m_A);

    for (BFMatrixColumnIterator it=m_A->begin(i);it!=m_A->end(i);it++)
      prod +=(*it-m_B->Peek(it.Row(),j))*(*it-m_B->Peek(it.Row(),j));
    
    if(ptr){
      
      for (BFMatrixColumnIterator it=m_B->begin(j);it!=m_B->end(j);it++){
	if(m_A->Peek(it.Row(),i)==0)
	  prod +=(*it)*(*it);
      }
      
    }
    
    return sqrt(prod)/m_A->Nrows();
  }
  
  
  /* NMI derived similarity kernel*/
  double simkernel::NMI(int i, int j){
    
    if(m_A->Nrows()<256)    
      hist.ReSize(m_A->Nrows()-1,m_A->Nrows()-1);
    else
      hist.ReSize(256,256);
    
    
    if(m_B.get()==0)   m_B=m_A;
    
    for (unsigned int s=1;s<=m_A->Nrows();s++){
     
	hist.AddSample(m_A->Peek(s,i),m_B->Peek(s,j));
     
    }  
   
    return  hist.normalisedmutualinformation();
  }
  
  //// for vectors  - weight function is optional
  
  double simkernel::corr(const vector<double> &A,const vector<double> &B,const vector<double> &weights){
    
    double prod,varA,varB;
    double sum=0;
    prod=0; varA=0; varB=0;
    
    if(A.size()!=B.size()){cout << A.size() << " " << B.size() << " " << weights.size() << endl; throw   DISCRETEOPTHOCRException("SIMILARITIES:: correlation, data dimensions do not match");}
    Initialize(2,A,B,weights);

    for (unsigned int s=0;s<A.size();s++){
      double varwght=0;
      if(weights.size()==A.size()){ varwght=weights[s]; }//cout << s << " varwght " << varwght <<  endl;}
      else varwght=1;

      prod = prod+varwght*(A[s] - _meanA)*(B[s] - _meanB); 
   
      varA = varA+varwght*(A[s] - _meanA)*(A[s] - _meanA);
      varB = varB+varwght*(B[s] - _meanB)*(B[s] - _meanB);
     
      sum+=varwght;
      
    }

    if(sum>0){
      prod/=sum;
      varA/=sum;
      varB/=sum;
    }

   
    if (varA == 0.0 || varB == 0.0) return  0.0; else
      return prod/(sqrt(varA)*sqrt(varB));
    
  }

  double simkernel::corr(const map<int,float> &A,const map<int,float> &B){
    
    double prod,varA,varB,val;
    double MAPA_MEAN=0,MAPB_MEAN=0;
    int ind=0;
    prod=0; varA=0; varB=0;
    
    /// estimate mean
    for (map<int, float>::const_iterator iter = A.begin(); iter != A.end(); ++iter){
      if(B.find(iter->first) != B.end()){
	MAPA_MEAN+=iter->second;
	MAPB_MEAN+=B.find(iter->first)->second;
	ind++;
      }
    }

    if(ind==0){ cout << " ind==0 << A.size() " << A.size() << " " << B.size() << endl; val=MAXSIM; }
    else{
   
      MAPA_MEAN/=ind;
      MAPB_MEAN/=ind;
     
      for (map<int, float>::const_iterator iter = A.begin(); iter != A.end(); ++iter){
	if(B.find(iter->first) != B.end()){
	  //if(iter->second!=iter->second) cout << " iter->second!=iter->second " << iter->second << endl;
	  //if(B.find(iter->first)->second!=B.find(iter->first)->second) cout << " iB.find(iter->first)->second!=B.find(iter->first)->second " << B.find(iter->first)->second << endl;

	  prod = prod+(iter->second - MAPA_MEAN)*(B.find(iter->first)->second - MAPB_MEAN); 
	  
	  varA = varA+(iter->second - MAPA_MEAN)*(iter->second - MAPA_MEAN);
	  varB = varB+(B.find(iter->first)->second - MAPB_MEAN)*(B.find(iter->first)->second - MAPB_MEAN);
     
	  //if(varA!=varA) cout << " varA " << varA << " MAPA_MEAN " << MAPA_MEAN << " iter->second " << iter->second << " ind " <<  ind << endl;
	  //if(varB!=varB) cout << " varB " << varA << " MAPB_MEAN " << MAPB_MEAN <<" B.find(iter->first)->second " << B.find(iter->first)->second << " ind " <<  ind << endl;

	}
      
      }

      
      prod/=ind;
      varA/=ind;
      varB/=ind;
    

      if (varA == 0.0 || varB == 0.0) val=0.0; 
      else {val=prod/(sqrt(varA)*sqrt(varB));

	val=1-(1+val)/2;
      }
      if(val!=val)       cout << "prod " << prod << " varA " << varA << " varB " << varB <<" ind " <<  ind << endl;

    }

    if(val!=val)       cout << " 2 prod " << prod << " varA " << varA << " varB " << varB <<" ind " <<  ind << endl;


    return val;
  }
  
  // used for calculating componentwise 'correlation'
  double simkernel::corrdebug(const int &i, Matrix &COMP, const vector<double> &A,const vector<double> &B,const vector<double> &weights){
    
    double prod,varA,varB;
    double sum=0;
    prod=0; varA=0; varB=0;
  
    if(A.size()!=B.size()){cout << A.size() << " " << B.size() << " " << weights.size() << endl; throw   DISCRETEOPTHOCRException("SIMILARITIES:: correlation, data dimensions do not match");}
    Initialize(2,A,B,weights);

    
    for (unsigned int s=0;s<A.size();s++){
      double varwght=0;
      if(weights.size()==A.size()) varwght=weights[s];
      else varwght=1;

      COMP(s+1,i+1)=varwght*(A[s] - _meanA)*(B[s] - _meanB);

      prod = prod+COMP(s+1,i+1); 
      

      varA = varA+varwght*(A[s] - _meanA)*(A[s] - _meanA);
      varB = varB+varwght*(B[s] - _meanB)*(B[s] - _meanB);
     
      sum+=varwght;
        
      }

    if(sum>0){
      prod/=sum;
      varA/=sum;
      varB/=sum;
    }
   
    for (unsigned int s=0;s<A.size();s++){
      COMP(s+1,i+1)/=(sum*sqrt(varA)*sqrt(varB));
      
    }
    
    if (varA == 0.0 || varB == 0.0) return  0.0; else
      return prod/(sqrt(varA)*sqrt(varB));
    
  }
  
  double simkernel::SSD(const vector<double> &A,const vector<double> &B,const vector<double> &weights){
  
    double prod=0,sum=0;

    if(A.size()!=B.size()){ throw   DISCRETEOPTHOCRException("SIMILARITIES:: SSD data dimensions do not match");}
    Initialize(3,A,B,weights);

    for (unsigned int s=0;s<A.size();s++){
      double varwght=0;
      if(weights.size()==A.size()) varwght=weights[s];
      else varwght=1;
      prod = prod+varwght*(A[s] - B[s])*(A[s] - B[s]); 
      sum+=varwght;
      
    }						  
    if(sum>0) prod=prod/sum;
    return prod;
  }

  
  double simkernel::nSSD(const vector<double> &A,const vector<double> &B,const vector<double> &weights){ /// normalized by variance - was hoped this would make SSD better for discrete opt but doesn't seem to be the case
    
    double prod=0,sum=0;
    double varA=0.0,varB=0.0;
    if(A.size()!=B.size()){ throw   DISCRETEOPTHOCRException("SIMILARITIES:: SSD data dimensions do not match");}
    
    // only perform calculations for non zero values
     for (unsigned int s=0;s<A.size();s++){
      // 
       double varwght=0;
      if(weights.size()==A.size()) varwght=weights[s];
      else varwght=1;
     
      prod = prod+varwght*(A[s] - B[s])*(A[s] - B[s]); 
      varA = varA+varwght*(A[s] - _meanA)*(A[s] - _meanA);
      varB = varB+varwght*(B[s] - _meanB)*(B[s] - _meanB);
      sum+=varwght;

    }

    if(sum>0){
      prod/=sum;
      varA/=sum;
      varB/=sum;
    }

    if((sqrt(varA)*sqrt(varB))>0) prod=prod/(sqrt(varA)*sqrt(varB));

    return prod;
  }

  

  double simkernel::NMI(const vector<double> &A,const vector<double> &B,const vector<double> &weights){
    double val;
    if(A.size()!=B.size()){ throw   DISCRETEOPTHOCRException("SIMILARITIES:: NMI data dimensions do not match");}

    Initialize(4,A,B,weights);

    if(A.size()<256)    
      hist.ReSize(A.size()-1,A.size()-1);
    else
      hist.ReSize(256,256);
  
    for (unsigned int s=0;s<A.size();s++){
      double varwght=0;
      if(weights.size()==A.size()) varwght=1;//weights[s];  // Haven't verified weighted NMI yet
      else varwght=1;
      hist.AddSample(varwght,A[s],B[s]);
    }  
    
    val= hist.normalisedmutualinformation();

  
    return  val;
  }


  double simkernel::get_sim(int input, int reference){

    double val=0.0;
 
    switch (_sim){
    case 1:
      val=SSD(input,reference);
      break;
    case 2:
      val=corr(input,reference);
      break;
    case 3:
      val=NMI(input,reference);
      break;
    }
    // }

  return val;
}

 
  double simkernel::get_sim_for_min(int input, int reference){

    double val=0.0;
  
    switch (_sim){
    case 1:
      val=-SSD(input,reference);
      break;
    case 2:
      val=corr(input,reference);
      val=1-(1+val)/2;
      break;
    case 3:
      val=NMI(input,reference);
      break;
    }
    // }

  return val;
}
 
 
  double simkernel::get_sim_for_min(const vector<double> &input,const vector<double> &reference,const vector<double> &weights){
    
    double val=0.0;
   
    switch (_sim){
    case 1:
      val=SSD(input,reference,weights);
    case 2:
      val=corr(input,reference,weights);
      val=1-(1+val)/2; // scale between 0 and 1
      break;
    case 3:
      val=NMI(input,reference);
      break;
    }
    

    return val;
  }

  
  double simkernel::get_sim_for_mindebug(const int & i, Matrix &COMP, const vector<double> &input,const vector<double> &reference,const vector<double> &weights){
    
    double val=0.0;
  
    switch (_sim){
    case 1:
      val=SSD(input,reference,weights);
      break;
    case 2:
      val=corrdebug(i,COMP,input,reference,weights);
      val=1-(1+val)/2;// scale between 0 and 1
      break;
    case 3:
      val=NMI(input,reference);
      break;
    }
    

    return val;
  }

/////ALPHA ENTROPY **NOT CURRENTLY AVAILABLE **** 


/* vector<double> simkernel::alpha_entropy(vector<vector<double> > &M, vector<double> &weights, int &kNN, const double & eps,vector<vector<int> > & neighbours)
  {
    int dim=M[0].size();
    int maxPts=M.size();
    ANNpointArray		dataPts;				// data points
    ANNpoint			queryPt;				// query point
    ANNkd_tree*			kdTree;					// search structure
    ANNdistArray dists;					// near neighbor distances
    ANNidxArray			nnIdx;					// near neighbor indices
    vector<double> SSD;
    double tot=0.0;
    double X,Y;
    X=Y=0;
    queryPt = annAllocPt(dim);					// allocate query point
    dataPts = annAllocPts(maxPts, dim);			// allocate data points
    nnIdx = new ANNidx[kNN];						// allocate near neigh indices
    dists = new ANNdist[kNN];						// allocate near neighbor dists
    int zeros;
    
    if(kNN>(int)M.size()){ cout << kNN << " datasize " << M.size() << endl; kNN= M.size()-1;}
    neighbours.clear();

    // cout << "Data Points:\n";
    for (unsigned int i=0;i<M.size();i++){
      for (unsigned int j=0;j<M[i].size();j++){
	dataPts[i][j]=weights[j]*M[i][j];  /// M will hold vertex feature intensities; data combines feature intensities 
	//Weights influence contributions of dfferent features (contributions are assumed fixed for each control point dataset)
      }
    }
   
    kdTree	 = new ANNbd_tree(			// build it
					dataPts,					// the data points
					M.size(),					// number of points
					dim,						// dimension of space
					50,				// maximum bucket size
					def_split,						// splitting rule
					def_shrink);					// shrinking rule

   
    for (unsigned int i=0;i<M.size();i++){ 
      vector<int> NN;
      kdTree->annkSearch(						// search
			 dataPts[i],						// query point
			 kNN,								// number of near neighbors
			 nnIdx,							// nearest neighbors (returned)
			 dists,							// distance (returned)
			 eps);							// error bound

      tot=0;   
      zeros=0;
      
      for (int k = 0; k < kNN; k++) {			// print summary
	tot+=dists[k];			// unsquare distance
	NN.push_back(nnIdx[k]);
      }
      neighbours.push_back(NN);
    
      SSD.push_back(tot);
    }
  
    annClose();									// done with ANN
    delete [] dists;
    delete [] nnIdx;	
    delete kdTree;
    return SSD;
  }
*/
///////
void fullsimkernel::calculate_sim_column(int ind){

    
  int cols;
    
  if(m_B==NULL)
    cols= m_A->Ncols();
  else
    cols= m_B->Ncols();
    
      
  for (int j=1; j <= cols; j++){
    switch (_sim){
    case 1:
      (*mp)(j,ind)=-SSD(ind,j);
      break;
    case 2:
      (*mp)(j,ind)=corr(ind,j);
      break;
    case 3:
      (*mp)(j,ind)=NMI(ind,j);
      break;
  
    }
  }
}
 
void fullsimkernel::calculate_full_sim_kernel(){
 
    int cols;
  
    
    /// Determine whether calulating self similarity (therefore m_B=m_A) 
    if(m_B==NULL)
      cols= m_A->Ncols();
    else
      cols= m_B->Ncols();  
   
    
    for (unsigned int i=1; i <= m_A->Ncols(); i++){
      for (int j=1; j <= cols; j++){
	switch (_sim){
	case 1:
	   (*mp)(j,i)=-SSD(i,j);
	  break;
	case 2:
	  (*mp)(j,i)=corr(i,j);
	  break;
	case 3:
	  (*mp)(j,i)=NMI(i,j);
      	  break;

	}
      }
    }
    
}

  void fullsimkernel::calculate_NN_sim_kernel(){
 
    
   
    
    for (unsigned int i=1; i <= m_A->Ncols(); i++){
	switch (_sim){
	case 1:
	   (*mp)(i,i)=SSD(i,i);
	  break;
	case 2:
	  (*mp)(i,i)=corr(i,i);
	  break;
	case 3:
	  (*mp)(i,i)=NMI(i,i);
      	  break;

	}
      
    }
    
}



void fullsimkernel::Print(const std::string fname) const
{
  if (!fname.length()) cout << endl << *mp << endl;
  else write_ascii_matrix(fname,*mp);
}






}


