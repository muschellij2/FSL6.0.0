
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <unistd.h>
#define WANT_STREAM
#define WANT_MATH
#include "newmat.h"
#include "newmatap.h"
#include "newmatio.h"

using namespace std;
using namespace NEWMAT;


////////////////////////////////////////////////////////////////////////////
// Stuff pasted from miscmaths.cc

ColumnVector cross(const ColumnVector& a, const ColumnVector& b)
{
  Tracer tr("cross");
  ColumnVector ans(3);
  ans(1) = a(2)*b(3) - a(3)*b(2);
  ans(2) = a(3)*b(1) - a(1)*b(3);
  ans(3) = a(1)*b(2) - a(2)*b(1);
  return ans;
}


ColumnVector cross(const Real *a, const Real *b)
{
  Tracer tr("cross");
  ColumnVector a1(3), b1(3);
  a1 << a;
  b1 << b;
  return cross(a1,b1);
}

template<class T> inline T Sqr(const T& x) { return x*x; }
template<class S, class T>
inline T Min(const S &a, const T &b) { if (a<b) return (T) a; else return b; }

template<class S, class T>
inline T Max(const S &a, const T &b) { if (a>b) return (T) a; else return b; }

double norm2(const ColumnVector& x) { return std::sqrt(x.SumSquare()); }

inline float dot(const ColumnVector& a, const ColumnVector& b) { return Sum(SP(a,b)); }

int decompose_aff(ColumnVector& params, const Matrix& affmat,
                  int (*rotmat2params)(ColumnVector& , const Matrix& ));
int decompose_aff(ColumnVector& params, const Matrix& affmat,
                  const ColumnVector& centre,
                  int (*rotmat2params)(ColumnVector& , const Matrix& ));
int compose_aff(const ColumnVector& params, int n, const ColumnVector& centre,
                Matrix& aff,
                int (*params2rotmat)(const ColumnVector& , int , Matrix& ,
                                     const ColumnVector& ) );

int diag(Matrix& m, const float diagvals[]);
int diag(Matrix& m, const ColumnVector& diagvals);
ReturnMatrix diag(const Matrix& mat);

int diag(Matrix& m, const float diagvals[])
{
  Tracer tr("diag");
  m=0.0;
  for (int j=1; j<=m.Nrows(); j++)
    m(j,j)=diagvals[j-1];
  return 0;
}


int diag(Matrix& m, const ColumnVector& diagvals)
{
  Tracer tr("diag");

  m.ReSize(diagvals.Nrows(),diagvals.Nrows());
  m=0.0;
  for (int j=1; j<=diagvals.Nrows(); j++)
    m(j,j)=diagvals(j);
  return 0;
}

ReturnMatrix diag(const Matrix& Mat)
{
  Tracer tr("diag");
  if(Mat.Ncols()==1){
    Matrix retmat(Mat.Nrows(),Mat.Nrows());
    diag(retmat,Mat);
    retmat.Release();
    return retmat;}
  else{
    int mindim = Min(Mat.Ncols(),Mat.Nrows());
    Matrix retmat(mindim,1);
    for(int ctr=1; ctr<=mindim;ctr++){
      retmat(ctr,1)=Mat(ctr,ctr);
    }
    retmat.Release();
    return retmat;
  }
}

bool isNumber(const string& str);
string skip_alpha(ifstream& fs);
ReturnMatrix read_ascii_matrix(const string& filename, int nrows, int ncols);
ReturnMatrix read_ascii_matrix(int nrows, int ncols, const string& filename);
ReturnMatrix read_ascii_matrix(const string& filename);
ReturnMatrix read_ascii_matrix(ifstream& fs, int nrows, int ncols);
ReturnMatrix read_ascii_matrix(int nrows, int ncols, ifstream& fs);
ReturnMatrix read_ascii_matrix(ifstream& fs);

int write_ascii_matrix(const Matrix& mat, const string& filename, int precision=-1);
int write_ascii_matrix(const string& filename, const Matrix& mat, int precision=-1);
int write_ascii_matrix(const Matrix& mat, ofstream& fs, int precision=-1);
int write_ascii_matrix(ofstream& fs, const Matrix& mat, int precision=-1);

int decompose_aff(ColumnVector& params, const Matrix& affmat,
                  const ColumnVector& centre,
                  int (*rotmat2params)(ColumnVector& , const Matrix& ))
{
  // decomposes using the convention: mat = rotmat * skew * scale
  // order of parameters is 3 rotation + 3 translation + 3 scales + 3 skews
  // angles are in radians
  Tracer tr("decompose_aff");
  if (params. Nrows() < 12)
    params.ReSize(12);
  if (rotmat2params==0)
    {
      cerr << "No rotmat2params function specified" << endl;
      return -1;
    }
  ColumnVector x(3), y(3), z(3);
  Matrix aff3(3,3);
  aff3 = affmat.SubMatrix(1,3,1,3);
  x = affmat.SubMatrix(1,3,1,1);
  y = affmat.SubMatrix(1,3,2,2);
  z = affmat.SubMatrix(1,3,3,3);
  float sx, sy, sz, a, b, c;
  sx = norm2(x);
  sy = std::sqrt( dot(y,y) - (Sqr(dot(x,y)) / Sqr(sx)) );
  a = dot(x,y)/(sx*sy);
  ColumnVector x0(3), y0(3);
  x0 = x/sx;
  y0 = y/sy - a*x0;
  sz = std::sqrt(dot(z,z) - Sqr(dot(x0,z)) - Sqr(dot(y0,z)));
  b = dot(x0,z)/sz;
  c = dot(y0,z)/sz;
  params(7) = sx;  params(8) = sy;  params(9) = sz;
  Matrix scales(3,3);
  float diagvals[] = {sx,sy,sz};
  diag(scales,diagvals);
  Real skewvals[] = {1,a,b,0 , 0,1,c,0 , 0,0,1,0 , 0,0,0,1};
  Matrix skew(4,4);
  skew  << skewvals;
  params(10) = a;  params(11) = b;  params(12) = c;
  Matrix rotmat(3,3);
  rotmat = aff3 * scales.i() * (skew.SubMatrix(1,3,1,3)).i();
  ColumnVector transl(3);
  transl = affmat.SubMatrix(1,3,1,3)*centre + affmat.SubMatrix(1,3,4,4)
    - centre;
  for (int i=1; i<=3; i++)  { params(i+3) = transl(i); }
  ColumnVector rotparams(3);
  (*rotmat2params)(rotparams,rotmat);
  for (int i=1; i<=3; i++)  { params(i) = rotparams(i); }
  return 0;
}

int decompose_aff(ColumnVector& params, const Matrix& affmat,
                  int (*rotmat2params)(ColumnVector& , const Matrix& ))
{
  Tracer tr("decompose_aff");
  ColumnVector centre(3);
  centre = 0.0;
  return decompose_aff(params,affmat,centre,rotmat2params);
}



int compose_aff(const ColumnVector& params, int n, const ColumnVector& centre,
                Matrix& aff,
                int (*params2rotmat)(const ColumnVector& , int , Matrix& ,
                                     const ColumnVector& ) )
{
  Tracer tr("compose_aff");
  if (n<=0) return 0;
  // order of parameters is 3 rotation + 3 translation + 3 scales + 3 skews
  // angles are in radians

  (*params2rotmat)(params,n,aff,centre);

  if (n<=6)  return 0;

  Matrix scale=IdentityMatrix(4);
  if (n>=7) {
    scale(1,1)=params(7);
    if (n>=8) scale(2,2)=params(8);
    else      scale(2,2)=params(7);
    if (n>=9) scale(3,3)=params(9);
    else      scale(3,3)=params(7);
  }
  // fix the translation so that the centre is not moved
  ColumnVector strans(3);
  strans = centre - scale.SubMatrix(1,3,1,3)*centre;
  scale.SubMatrix(1,3,4,4) = strans;

  Matrix skew=IdentityMatrix(4);
  if (n>=10) {
    if (n>=10) skew(1,2)=params(10);
    if (n>=11) skew(1,3)=params(11);
    if (n>=12) skew(2,3)=params(12);
  }
  // fix the translation so that the centre is not moved
  ColumnVector ktrans(3);
  ktrans = centre - skew.SubMatrix(1,3,1,3)*centre;
  skew.SubMatrix(1,3,4,4) = ktrans;

  aff = aff * skew * scale;

  return 0;

}

// General string/IO functions
bool isNumber( const string& input)
{
  if (input.size()==0) return false;
  char *pend;
  strtod(input.c_str(),&pend);
  if (*pend!='\0') return false;
  return true;
}

// Use this to skip all lines that contain non-numeric entries, and return the first line starting with a number
//  and the file pointer is reset to the beginning of the first line that starts with a number
string skip_alpha(ifstream& fs)
{
  string cline;
  while (!fs.eof()) {
    streampos curpos = fs.tellg();
    // read a line, then turn it into a stream in order to read out the first token
    getline(fs,cline);
    cline += " "; // force extra entry in parsing (ensure at least one token is read)
    istringstream ss(cline.c_str());
    string firstToken="";
    ss >> firstToken; // Put first non-whitespace sequence into firstToken
    if (isNumber(firstToken)) {
      if (fs.eof()) { fs.clear(); }  // should only be executed if the file had no valid line terminators
      fs.seekg(curpos);
      return cline;
    }
  }
  return "";
}

ReturnMatrix read_ascii_matrix(const string& filename)
{
  Matrix mat;
  if ( filename.size()<1 ) return mat;
  ifstream fs(filename.c_str());
  if (!fs) {
    cerr << "Could not open matrix file " << filename << endl;
    mat.Release();
    return mat;
  }
  mat = read_ascii_matrix(fs);
  fs.close();
  mat.Release();
  return mat;
}

ReturnMatrix read_ascii_matrix(ifstream& fs)
{
  int nRows(0), nColumns(0);
  string currentLine;
  // skip initial non-numeric lines
  //  and count the number of columns in the first numeric line
  currentLine = skip_alpha(fs);
  currentLine += " ";
  {
    istringstream ss(currentLine.c_str());
    string dummyToken="";
    while (!ss.eof()) {
      nColumns++;
      ss >> dummyToken;
    }
  }
  nColumns--;

  // count the number of lines that start with a number (don't worry if they don't have enough numbers at this stage)
  do {
    getline(fs,currentLine);
    currentLine += " "; // force extra entry in parsing
    istringstream ss(currentLine.c_str());
    string firstToken("");
    ss >> firstToken; //Put first non-whitespace sequence into cc
    if (isNumber(firstToken)) nRows++;  // add new row to matrix
  } while (!fs.eof());

  // now know the size of matrix
  fs.clear();
  fs.seekg(0,ios::beg);
  return read_ascii_matrix(fs,nRows,nColumns);

}

ReturnMatrix read_ascii_matrix(ifstream& fs, int nrows, int ncols)
{
  Matrix mat(nrows,ncols);
  mat = 0.0;
  string ss="";
  double newnum;

  ss = skip_alpha(fs);
  for (int r=1; r<=nrows; r++) {
    istringstream sline(ss.c_str());
    for (int c=1; c<=ncols; c++) {
      sline >> newnum;
      if ( sline.eof() ) {
        throw Exception("Could not find enough numbers in matrix file");
      }
      mat(r,c) = newnum;
    }
    if ( r!=nrows ) {
      getline(fs,ss);  // this is processed now, so move the stream past it
      ss = skip_alpha(fs);
    }
  }
  mat.Release();
  return mat;
}

int write_ascii_matrix(const string& filename, const Matrix& mat,
                       int precision)
{
  return write_ascii_matrix(mat, filename, precision);
}

int write_ascii_matrix(const Matrix& mat, const string& filename,
                       int precision)
{
  Tracer tr("write_ascii_matrix");
  if ( (filename.size()<1) ) return -1;
  ofstream fs(filename.c_str());
  if (!fs) {
    cerr << "Could not open file " << filename << " for writing" << endl;
    return -1;
  }
  int retval = write_ascii_matrix(mat,fs,precision);
  fs.close();
  return retval;
}

int write_ascii_matrix(ofstream& fs, const Matrix& mat,
                       int precision)
{
  return write_ascii_matrix(mat, fs, precision);
}

int write_ascii_matrix(const Matrix& mat, ofstream& fs, int precision)
{
  //fs.setf(ios::floatfield);  // use fixed or scientific notation as appropriate
  fs.setf(std::ios::fixed, std::ios::floatfield);
  if (precision>0)  {
    fs.precision(precision);
  } else {
    fs.precision(10);    // default precision
  }
#ifdef PPC64
  int n=0;
#endif
  for (int i=1; i<=mat.Nrows(); i++) {
    for (int j=1; j<=mat.Ncols(); j++) {
      fs << mat(i,j) << "  ";
#ifdef PPC64
      if ((n++ % 50) == 0) fs.flush();
#endif
    }
    fs << endl;
  }
  return 0;
}


int make_rot(const ColumnVector& angl, const ColumnVector& centre,
             Matrix& rot);

int rotmat2euler(ColumnVector& angles, const Matrix& rotmat);

int construct_rotmat_euler(const ColumnVector& params, int n, Matrix& aff);
int construct_rotmat_euler(const ColumnVector& params, int n, Matrix& aff,
                           const ColumnVector& centre);

int make_rot(const ColumnVector& angl, const ColumnVector& centre,
             Matrix& rot)
{
  // Matrix rot must be 4x4; angl and orig must be length 3
  Tracer tr("make_rot");
  rot=IdentityMatrix(4);  // default return value
  float theta;
  theta = norm2(angl);
  if (theta<1e-8) {  // avoid round-off errors and return Identity
    return 0;
  }
  ColumnVector axis = angl/theta;
  ColumnVector x1(3), x2(3), x3(3);
  x1 = axis;
  x2(1) = -axis(2);  x2(2) = axis(1);  x2(3) = 0.0;
  if (norm2(x2)<=0.0) {
    x2(1) = 1.0;  x2(2) = 0.0;  x2(3) = 0.0;
  }
  x2 = x2/norm2(x2);
  x3 = cross(x1,x2);
  x3 = x3/norm2(x3);

  Matrix basischange(3,3);
  basischange.SubMatrix(1,3,1,1) = x2;
  basischange.SubMatrix(1,3,2,2) = x3;
  basischange.SubMatrix(1,3,3,3) = x1;

  Matrix rotcore=IdentityMatrix(3);
  rotcore(1,1)=cos(theta);
  rotcore(2,2)=cos(theta);
  rotcore(1,2)=sin(theta);
  rotcore(2,1)=-sin(theta);

  rot.SubMatrix(1,3,1,3) = basischange * rotcore * basischange.t();

  Matrix ident3=IdentityMatrix(3);
  ColumnVector trans(3);
  trans = (ident3 - rot.SubMatrix(1,3,1,3))*centre;
  rot.SubMatrix(1,3,4,4)=trans;
  return 0;
}


int construct_rotmat_euler(const ColumnVector& params, int n, Matrix& aff,
                           const ColumnVector& centre)
{
  Tracer tr("construct_rotmat_euler");
  ColumnVector angl(3);
  Matrix newaff(4,4);
  aff=IdentityMatrix(4);

  if (n<=0) return 0;
  // order of parameters is 3 rotation + 3 translation
  // angles are in radians
  //  order of parameters is (Rx,Ry,Rz) and R = Rx.Ry.Rz
  angl=0.0;
  angl(1)=params(1);
  make_rot(angl,centre,newaff);
  aff = aff * newaff;
  if (n==1) return 0;

  angl=0.0;
  angl(2)=params(2);
  make_rot(angl,centre,newaff);
  aff = aff * newaff;
  if (n==2) return 0;

  angl=0.0;
  angl(3)=params(3);
  make_rot(angl,centre,newaff);
  aff = aff * newaff;
  if (n==3) return 0;

  aff(1,4)+=params(4);
  if (n==4) return 0;
  aff(2,4)+=params(5);
  if (n==5) return 0;
  aff(3,4)+=params(6);
  if (n==6) return 0;

  return 1;
}

int construct_rotmat_euler(const ColumnVector& params, int n, Matrix& aff)
{
  Tracer tr("construct_rotmat_euler");
  ColumnVector centre(3);
  centre = 0.0;
  return construct_rotmat_euler(params,n,aff,centre);
}

int rotmat2euler(ColumnVector& angles, const Matrix& rotmat)
{
  // uses the convention that R = Rx.Ry.Rz
  Tracer tr("rotmat2euler");
  float cz, sz, cy, sy, cx, sx;
  cy = std::sqrt(Sqr(rotmat(1,1)) + Sqr(rotmat(1,2)));
  if (cy < 1e-4) {
    //cerr << "Cos y is too small - Gimbal lock condition..." << endl;
    cx = rotmat(2,2);
    sx = -rotmat(3,2);
    sy = -rotmat(1,3);
    angles(1) = atan2(sx,cx);
    angles(2) = atan2(sy,(float)0.0);
    angles(3) = 0.0;
  } else {
    // choose by convention that cy > 0
    // get the same rotation if: sy stays same & all other values swap sign
    cz = rotmat(1,1)/cy;
    sz = rotmat(1,2)/cy;
    cx = rotmat(3,3)/cy;
    sx = rotmat(2,3)/cy;
    sy = -rotmat(1,3);
    //atan2(sin,cos)  (defined as atan2(y,x))
    angles(1) = atan2(sx,cx);
    angles(2) = atan2(sy,cy);
    angles(3) = atan2(sz,cz);
  }
  return 0;
}


////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////
// the real defaults are provided in the function parse_command_line

class globaloptions {
public:
  string testfname;
  string reffname;
  string outputmatascii;
  string initmatfname;
  string concatfname;
  string fixfname;
  string intervolfname;
  string xfm_type;
  int verbose;
  bool inverse;
  bool matonly;
public:
  globaloptions();
  ~globaloptions() {};
};

globaloptions globalopts;


globaloptions::globaloptions()
{
  // set up defaults
  reffname = "";

  testfname = "";
  outputmatascii = "";
  initmatfname = "";
  concatfname = "";
  fixfname = "";
  intervolfname = "";
  xfm_type = "a";
  inverse = false;
  matonly = true;
}


////////////////////////////////////////////////////////////////////////////

// Parsing functions for command line parameters

void print_usage(int argc, char *argv[])
{
  const string version="2.1";

  cout << "convert_xfm (Version " << version << ")" << endl
       << "Tool for manipulating FSL transformation matrices" << endl
       << "Copyright(c) 1999-2007, University of Oxford (Mark Jenkinson)" << endl
       << endl
       << "Usage: " << argv[0] << " [options] <input-matrix-filename>" << endl
       << "  e.g. " << argv[0] << " -omat <outmat> -inverse <inmat>" << endl
       << "       " << argv[0] << " -omat <outmat_AtoC> -concat <mat_BtoC> <mat_AtoB>" << endl << endl
       << "  Available options are:" << endl
       << "        -omat <matrix-filename>            (4x4 ascii format)" << endl
       << "        -concat <second-matrix-filename>" << endl
       << "        -fixscaleskew <second-matrix-filename>" << endl
       << "        -inverse                           (Reference image must be the one originally used)" << endl
       << "        -help" << endl;
}


void parse_command_line(int argc, char* argv[])
{
  if(argc<2){
    print_usage(argc,argv);
    exit(0);
  }


  int n=1;
  string arg;
  char first;
  bool initmatset=false;

  while (n<argc) {
    arg=argv[n];
    if (arg.size()<1) { n++; continue; }
    first = arg[0];
    if (first!='-') {
      if (initmatset) {
	cerr << "Unknown option " << arg << endl << endl;
	exit(0);
      } else {
	globalopts.initmatfname = arg;
	initmatset=true;
      }
      n++;
      continue;
    }

    // put options without arguments here
    if ( arg == "-help" ) {
      print_usage(argc,argv);
      exit(0);
    } else if ( arg == "-inverse" ) {
      globalopts.inverse = true;
      n++;
      continue;
    } else if ( arg == "-v" ) {
      globalopts.verbose = 5;
      n++;
      continue;
    }

    if (n+1>=argc)
      {
	cerr << "Lacking argument to option " << arg << endl << endl;
	exit(0);
      }

    // put options with 1 argument here
    if ( arg == "-concat") {
      globalopts.concatfname = argv[n+1];
      n+=2;
      continue;
    } else if ( arg == "-fixscaleskew") {
      globalopts.fixfname = argv[n+1];
      n+=2;
      continue;
    } else if ( arg == "-omat") {
      globalopts.outputmatascii = argv[n+1];
      n+=2;
      continue;
    } else if ( arg == "-verbose") {
      globalopts.verbose = atoi(argv[n+1]);
      n+=2;
      continue;
    } else {
      cerr << "Unrecognised option " << arg << endl << endl;
      exit(0);
    }

  }  // while (n<argc)

  if (globalopts.initmatfname.size()<1) {
    cerr << "Input matrix filename not found" << endl << endl;
    //print_usage(argc,argv);
    exit(0);
  }
}

////////////////////////////////////////////////////////////////////////////

int vector2affine(const ColumnVector& params, Matrix& aff)
{
  // order of parameters is 3 rotation + 3 translation + 3 scales + 3 skews
  // angles are in radians

  ColumnVector centre(3);
  centre = 0;
  compose_aff(params,12,centre,aff,construct_rotmat_euler);
  return 0;
}


int affmat2vector(const Matrix& aff, ColumnVector& params)
{
  ColumnVector centre(3);
  centre = 0;
  decompose_aff(params,aff,centre,rotmat2euler);
  return 0;
}


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

int main(int argc,char *argv[])
{
  parse_command_line(argc,argv);

  // read matrices
  Matrix affmat(4,4);
  affmat = read_ascii_matrix(globalopts.initmatfname);
  if (affmat.Nrows()<4) {
    cerr << "Cannot read input-matrix" << endl;
    return 0;
  }


  if (globalopts.fixfname.size() >= 1) {
    Matrix fixmat(4,4);
    fixmat = read_ascii_matrix(globalopts.fixfname);

    if (fixmat.Nrows()<4) {
      cerr << "Cannot read fixscaleskew-matrix" << endl;
      return 0;
    } else {
      if (globalopts.verbose>2) {
	cout << "Initial matrix:" << endl << affmat << endl;
	cout << "Fix Scale-Skew matrix:" << endl << fixmat << endl;
      }
      // do the work of combining scale/skew from fix and rest from init
      ColumnVector initp(12), fixp(12), combp(12);
      affmat2vector(affmat,initp);
      affmat2vector(fixmat,fixp);
      combp.SubMatrix(1,6,1,1) = initp.SubMatrix(1,6,1,1);
      combp.SubMatrix(7,12,1,1) = fixp.SubMatrix(7,12,1,1);
      vector2affine(combp,affmat);
    }
  }


  if (globalopts.concatfname.size() >= 1) {
    Matrix concatmat(4,4);
    concatmat = read_ascii_matrix(globalopts.concatfname);

    if (concatmat.Nrows()<4) {
      cerr << "Cannot read concat-matrix" << endl;
      return 0;
    } else {
      if (globalopts.verbose>2) {
	cout << "Initial matrix:" << endl << affmat << endl;
	cout << "Second matrix:" << endl << concatmat << endl;
      }
      affmat = concatmat * affmat;
    }
  }

  // apply inverse (if requested)
  if (globalopts.inverse) {
    affmat = affmat.i();
  }


  // Write outputs
  if (globalopts.outputmatascii.size() >= 1) {
    write_ascii_matrix(affmat,globalopts.outputmatascii);
  }

  if (globalopts.verbose>0) {
    cout << affmat << endl;
  }

  return 0;
}
