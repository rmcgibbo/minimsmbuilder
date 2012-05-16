//=============================================================================================
// Calculation of RMSD by a the quaternion-based characteristic polynomial (QCP) algorithm of Theobald [1].
// 
// [1] Theobald DL. Rapid calculation of RMSDs using a quaternion-based characteristic polynomial. 
//     Acta Cryst., A61:478, 2005.  doi:10.1107/50108767305015266
//
// Written by John D. Chodera <jchodera@gmail.com>, Dill lab, UCSF, 2006.
//=============================================================================================

#include <math.h>
#include <iostream>
using namespace std;
//#include "theobald_rmsd.h"

real abs(real r){
  if(r<0)
    return -r;
  return r;
}

real ls_rmsd(int N, rvec* x_original_nk, rvec* y_original_nk,int* atomindicies,int natomindicies){
  real x_centroid_k[3], y_centroid_k[3];
  real x_nk[natomindicies][3], y_nk[natomindicies][3];
  int nIndex;
  real G_x,G_y;
  real M[3][3];
  real K[4][4];
  real C_2,C_1,C_0;
  real lambda, lambda_old;
  int iteration;
  const int maxits=50;
  const real tolerance=1.0e-6;
  real lambda2,a,b;
  real rmsd2;
  real result;
  
  for(int k=0;k<3;k++){
    x_centroid_k[k]=0.0;
    y_centroid_k[k]=0.0;
  }
  for(int i=0;i<natomindicies;i++){
    for(int k=0;k<3;k++){
      x_centroid_k[k]+=x_original_nk[atomindicies[i]][k];
      y_centroid_k[k]+=y_original_nk[atomindicies[i]][k];
    }
  }
  for(int k=0;k<3;k++){
    x_centroid_k[k]/=((real)natomindicies);
    y_centroid_k[k]/=((real)natomindicies);
  }

  for(int i=0;i<natomindicies;i++){
    for(int k=0;k<3;k++){
      x_nk[i][k]=x_original_nk[atomindicies[i]][k]-x_centroid_k[k];
      y_nk[i][k]=y_original_nk[atomindicies[i]][k]-y_centroid_k[k];
    }
  }

  G_x=0.0;
  G_y=0.0;
  for(int i=0;i<natomindicies;i++){
    for(int k=0;k<3;k++){
      G_x+=x_nk[i][k]*x_nk[i][k];
      G_y+=y_nk[i][k]*y_nk[i][k];
    }
  }

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      M[i][j]=0.0;
      for(int k=0;k<natomindicies;k++)
        M[i][j]+=x_nk[k][i]*y_nk[k][j];
    }
  }

  K[0][0]=M[0][0]+M[1][1]+M[2][2];
  K[0][1]=M[1][2]-M[2][1];
  K[0][2]=M[2][0]-M[0][2];
  K[0][3]=M[0][1]-M[1][0];
  K[1][0]=K[0][1];
  K[1][1]=M[0][0]-M[1][1]-M[2][2];
  K[1][2]=M[0][1]+M[1][0];
  K[1][3]=M[2][0]+M[0][2];
  K[2][0]=K[0][2];
  K[2][1]=K[1][2];
  K[2][2]=-M[0][0]+M[1][1]-M[2][2];
  K[2][3]=M[1][2]+M[2][1];
  K[3][0]=K[0][3];
  K[3][1]=K[1][3];
  K[3][2]=K[2][3];
  K[3][3]=-M[0][0]-M[1][1]+M[2][2];

  C_2=0.0;
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      C_2-=2.0*M[i][j]*M[i][j];

  C_1=-8.0*
  (M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
  -M[1][0]*(M[0][1]*M[2][2]-M[0][2]*M[2][1])
  +M[2][0]*(M[0][1]*M[1][2]-M[0][2]*M[1][1]));

  C_0=
  K[0][0]*(
  K[1][1]*(K[2][2]*K[3][3]-K[2][3]*K[3][2])
  -K[2][1]*(K[1][2]*K[3][3]-K[1][3]*K[3][2])
  +K[3][1]*(K[1][2]*K[2][3]-K[1][3]*K[2][2])
  )-K[1][0]*(
  K[0][1]*(K[2][2]*K[3][3]-K[2][3]*K[3][2])
  -K[2][1]*(K[0][2]*K[3][3]-K[0][3]*K[3][2])
  +K[3][1]*(K[0][2]*K[2][3]-K[0][3]*K[2][2])
  )+K[2][0]*(
  K[0][1]*(K[1][2]*K[3][3]-K[1][3]*K[3][2])
  -K[1][1]*(K[0][2]*K[3][3]-K[0][3]*K[3][2])
  +K[3][1]*(K[0][2]*K[1][3]-K[0][3]*K[1][2])
  )-K[3][0]*(
  K[0][1]*(K[1][2]*K[2][3]-K[1][3]*K[2][2])
  -K[1][1]*(K[0][2]*K[2][3]-K[0][3]*K[2][2])
  +K[2][1]*(K[0][2]*K[1][3]-K[0][3]*K[1][2])
  );

  lambda=(G_x+G_y)/2.0;

  for(iteration=0;iteration<maxits;iteration++){
    lambda_old=lambda;
    lambda2=lambda_old*lambda_old;
    b=(lambda2+C_2)*lambda_old;
    a=b+C_1;
    lambda=lambda_old-(a*lambda_old+C_0)/(2.0*lambda2*lambda_old+b+a);
    if(abs(lambda-lambda_old)<abs(tolerance*lambda))
      break;
  }
  rmsd2=(G_x+G_y-2.0*lambda)/((real)natomindicies);

  result=0.0;
  if(rmsd2>0) result=sqrt(rmsd2);
  return result;
}

//This rmsd function is improved for efficiency, it assumes that the proper atom indicies have already been chosen
//and that the center of the atoms is the origin.  It also returns the msd, rather than the rmsd
//This could be further improved by caching the G_x and G_y, or by optimizing the coeficient computations as described in the paper
real ls_rmsd2(int N, rvec* x_nk, rvec* y_nk){
  int nIndex;
  real G_x,G_y;
  real M[3][3];
  real K[4][4];
  real C_2,C_1,C_0;
  real lambda, lambda_old;
  int iteration;
  const int maxits=50;
  const real tolerance=1.0e-6;
  real lambda2,a,b;
  real rmsd2;
  real result;
  
  G_x=0.0;
  G_y=0.0;
  for(int i=0;i<N;i++){
    for(int k=0;k<3;k++){
      G_x+=x_nk[i][k]*x_nk[i][k];
      G_y+=y_nk[i][k]*y_nk[i][k];
    }
  }

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      M[i][j]=0.0;
      for(int k=0;k<N;k++)
        M[i][j]+=x_nk[k][i]*y_nk[k][j];
    }
  }

  K[0][0]=M[0][0]+M[1][1]+M[2][2];
  K[0][1]=M[1][2]-M[2][1];
  K[0][2]=M[2][0]-M[0][2];
  K[0][3]=M[0][1]-M[1][0];
  K[1][0]=K[0][1];
  K[1][1]=M[0][0]-M[1][1]-M[2][2];
  K[1][2]=M[0][1]+M[1][0];
  K[1][3]=M[2][0]+M[0][2];
  K[2][0]=K[0][2];
  K[2][1]=K[1][2];
  K[2][2]=-M[0][0]+M[1][1]-M[2][2];
  K[2][3]=M[1][2]+M[2][1];
  K[3][0]=K[0][3];
  K[3][1]=K[1][3];
  K[3][2]=K[2][3];
  K[3][3]=-M[0][0]-M[1][1]+M[2][2];

  C_2=0.0;
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      C_2-=2.0*M[i][j]*M[i][j];

  C_1=-8.0*
  (M[0][0]*(M[1][1]*M[2][2]-M[1][2]*M[2][1])
  -M[1][0]*(M[0][1]*M[2][2]-M[0][2]*M[2][1])
  +M[2][0]*(M[0][1]*M[1][2]-M[0][2]*M[1][1]));

  C_0=
  K[0][0]*(
  K[1][1]*(K[2][2]*K[3][3]-K[2][3]*K[3][2])
  -K[2][1]*(K[1][2]*K[3][3]-K[1][3]*K[3][2])
  +K[3][1]*(K[1][2]*K[2][3]-K[1][3]*K[2][2])
  )-K[1][0]*(
  K[0][1]*(K[2][2]*K[3][3]-K[2][3]*K[3][2])
  -K[2][1]*(K[0][2]*K[3][3]-K[0][3]*K[3][2])
  +K[3][1]*(K[0][2]*K[2][3]-K[0][3]*K[2][2])
  )+K[2][0]*(
  K[0][1]*(K[1][2]*K[3][3]-K[1][3]*K[3][2])
  -K[1][1]*(K[0][2]*K[3][3]-K[0][3]*K[3][2])
  +K[3][1]*(K[0][2]*K[1][3]-K[0][3]*K[1][2])
  )-K[3][0]*(
  K[0][1]*(K[1][2]*K[2][3]-K[1][3]*K[2][2])
  -K[1][1]*(K[0][2]*K[2][3]-K[0][3]*K[2][2])
  +K[2][1]*(K[0][2]*K[1][3]-K[0][3]*K[1][2])
  );

  lambda=(G_x+G_y)/2.0;

  for(iteration=0;iteration<maxits;iteration++){
    lambda_old=lambda;
    lambda2=lambda_old*lambda_old;
    b=(lambda2+C_2)*lambda_old;
    a=b+C_1;
    lambda=lambda_old-(a*lambda_old+C_0)/(2.0*lambda2*lambda_old+b+a);
    if(abs(lambda-lambda_old)<abs(tolerance*lambda))
      break;
  }
  rmsd2=(G_x+G_y-2.0*lambda)/((real)N);

  result=0.0;
  if(rmsd2>0) result=rmsd2;
  return result;
}
