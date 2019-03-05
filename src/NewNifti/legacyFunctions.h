//This file contains code from niftio1_io.c and niftio1_.h, as released to
//the public domain by Rober W Cox, August 2003 and modified by
//Mark Jenkinson, August 2004 and Rick Reynolds, December 2004 and
//Matthew Webster 2018

typedef struct {                   /** 4x4 matrix struct **/
  float m[4][4] ;
} mat44 ;

typedef struct {                   /** 3x3 matrix struct **/
  float m[3][3] ;
} mat33 ;

#define NIFTI_L2R  1    /* Left to Right         */
#define NIFTI_R2L  2    /* Right to Left         */
#define NIFTI_P2A  3    /* Posterior to Anterior */
#define NIFTI_A2P  4    /* Anterior to Posterior */
#define NIFTI_I2S  5    /* Inferior to Superior  */
#define NIFTI_S2I  6    /* Superior to Inferior  */

mat44 nifti_quatern_to_mat44( float qb, float qc, float qd, float qx, float qy, float qz, float dx, float dy, float dz, float qfac );
void nifti_mat44_to_quatern( const mat44& R , double& qb, double& qc, double& qd, double& qx, double& qy, double& qz, double& dx, double& dy, double& dz, double& qfac );
void nifti_mat44_to_orientation( mat44 R , int *icod, int *jcod, int *kcod );
mat33 mat44_to_mat33(mat44 x);
float nifti_mat33_determ( mat33 R );
mat33 nifti_mat33_inverse( mat33 R );
mat44 nifti_mat44_inverse( mat44 R );
