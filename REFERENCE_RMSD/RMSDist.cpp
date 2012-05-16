/************************************************************/
/* RMSDist.cpp
/* 
/* Class for computing the RMSD between two Conformations.
/* Inherits from DistanceMetric.
/*
/* Please reference
/* GR Bowman, X Huang, and VS Pande. Methods 2009. Using generalized ensemble
/* simulations and Markov state models to identify conformational states.
/* 
/* Written by Gregory R. Bowman
/* Biophysics Program, Stanford Un iversity
/* Pande Group
/* 11/14/2008
/*
/* Copyright (C) 2008  Stanford University
/*
/* This program is free software; you can redistribute it and/or modify
/* it under the terms of the GNU General Public License as published by
/* the Free Software Foundation; either version 2 of the License, or
/* (at your option) any later version.
/*
/* This program is distributed in the hope that it will be useful,
/* but WITHOUT ANY WARRANTY; without even the implied warranty of
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
/* GNU General Public License for more details.
/*
/* You should have received a copy of the GNU General Public License
/* along with this program; if not, write to the Free Software
/* Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
/*
/************************************************************/
/* TODO:
/* 
/************************************************************/
/* CHANGE LOG:
/*
/************************************************************/
#include <assert.h>
#include <iostream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "RMSDist.h"

// lapack functions
extern "C"  void sspev_(char &jobz, char &uplo, int &n, float *ap, float *w, float *z, int &ldz, float *work, int &info);

extern "C" void sgetrf_(int &m, int &n, float *a, int &lda, int* ipiv, int &info);

extern "C"  float sdot_(int &n, float* dx, int &incx, float* dy, int & incy);

extern "C"  void sgemm_(char& transa,char& transb, int &m, int &n, int &k, float &alpha, float *a, int &lda, float *b, int &ldb, float &beta, float *c, int &ldc);

extern "C" float snrm2_(int& n, float* x, int& incx);

extern "C" void sscal_(int& n,float& da,float* dx, int& incx);
// end of lapack functions

using namespace std;

/************************************************************/
/* RMSDist()
/*   Constructor
/* Arguemtns:
/*   None
/* Return:
/*   None
/************************************************************/
RMSDist::RMSDist() {
  centered = false;
  debug = false;
}

/************************************************************/
/* RMSDist()
/*   Constructor
/* Arguemtns:
/*   int s = size, or number of atoms in a Conformation
/*   bool center = whether or not the Conformations have already been centered
/*   bool deb = whether want to run in debug mode
/* Return:
/*   None
/************************************************************/
RMSDist::RMSDist(int s, bool center, bool deb) {
  size = s;
  centered = center;
  debug = deb;
}

/************************************************************/
/* getDist()
/*   Get the RMSD between two conformations using the theobald algorithm.
/* Arguemtns:
/*   Element *conf1 = the first Conformation
/*   Element *conf2 = the second Conformation
/* Return:
/*   The RMSD between conf1 and conf2 (float).
/************************************************************/
float RMSDist::getDist(Element *conf1, Element *conf2) {
  Conformation* c1 = (Conformation*)conf1;
  Conformation* c2 = (Conformation*)conf2;

  int natom = size;
  float* mcA = new float[natom * 3];
  assert(mcA != NULL);
  float* mcB = new float[natom * 3];
  assert(mcB != NULL);

  int inc = 1, n = natom, m = 3;
  float alpha = 1, beta = 0;
  char trans = 't', ntrans = 'n';

  if(centered){
    memcpy(mcA, c1->pos, natom * 3 * sizeof(float));
    memcpy(mcB, c2->pos, natom * 3 * sizeof(float));
  }
  else{
    float center[3];
    float* ones = new float[natom];
    assert(ones != NULL);
    for(unsigned int i = 0; i < natom; i ++){
      ones[i] = 1;
    }   

    //center A
    for(unsigned int i = 0; i < 3; i ++){
      center[i] = sdot_(n, &(c1->pos[i*n]), inc, ones, inc) / n;
    }
    //--------------------------------------------------
    // if(debug){
    //  cout<<"centerA: "<<center[0]<<" "<<center[1]<<" "<<center[2]<<endl;
    // }
    //-------------------------------------------------- 
    for(unsigned int j = 0; j < 3; j ++){
      for(unsigned int i = 0; i < natom; i ++){
        mcA[i + j * natom] = c1->pos[i + j * natom] - center[j];
      }
    }

    //center B
    for(unsigned int i = 0; i < 3; i ++){
      center[i] = sdot_(n, &(c2->pos[i*n]), inc, ones, inc) / n;
    }
    //--------------------------------------------------
    // if(debug){
    //  cout<<"centerB: "<<center[0]<<" "<<center[1]<<" "<<center[2]<<endl;
    // }
    //-------------------------------------------------- 
    for(unsigned int j = 0; j < 3; j ++){
      for(unsigned int i = 0; i < natom; i ++){
        mcB[i + j * natom] = c2->pos[i + j * natom] - center[j];
      }
    }

    delete []ones;
  }

  float G_x = 0, G_y = 0;
  for(unsigned int i = 0; i < 3 * natom; i ++){
    G_x += (mcA[i] * mcA[i]);
    G_y += (mcB[i] * mcB[i]);
  }
  
  //--------------------------------------------------
  // if(debug){
  //  cout<<"G_x: "<<G_x<<" G_y: "<<G_y<<endl;
  // }
  //-------------------------------------------------- 

  //compute R = mcA' * mcB
  float* M = new float[m * m];
  assert(M != NULL);
  memset(M, 0, m * m * sizeof(float));
  sgemm_(trans, ntrans, m, m, n, alpha, mcA, n, mcB, n, beta, M, m);

  int m1 = m + 1;
  float* K = new float[m1 * m1];
  assert(K != NULL);
   K[0+0*m1] = M[0+0*m] + M[1+1*m] + M[2+2*m];
  K[0+1*m1] = M[1+2*m] - M[2+1*m];
  K[0+2*m1] = M[2+0*m] - M[0+2*m];
  K[0+3*m1] = M[0+1*m] - M[1+0*m];
  K[1+0*m1] = K[0+1*m1];
  K[1+1*m1] = M[0+0*m] - M[1+1*m] - M[2+2*m];
  K[1+2*m1] = M[0+1*m] + M[1+0*m];
  K[1+3*m1] = M[2+0*m] + M[0+2*m];
  K[2+0*m1] = K[0+2*m1];
  K[2+1*m1] = K[1+2*m1];
  K[2+2*m1] = -M[0+0*m] + M[1+1*m] - M[2+2*m];
  K[2+3*m1] = M[1+2*m] + M[2+1*m];
  K[3+0*m1] = K[0+3*m1];
  K[3+1*m1] = K[1+3*m1];
  K[3+2*m1] = K[2+3*m1];
  K[3+3*m1] = -M[0+0*m] - M[1+1*m] + M[2+2*m];
  
   float C_4 = 1.0, C_3 = 0.0,  C_2,  C_1,  C_0;
  
  C_2 = 0;
  for(unsigned int i = 0; i < m * m; i ++){
    C_2 += M[i] * M[i];
  }
  C_2 *= -2;

  float detM, detK;
  int *ipiv = new int[m];
  assert(ipiv != NULL);
  int info = 0;
  sgetrf_(m, m, M, m, ipiv, info);
  
  //cout<<ipiv[0] << " "<<ipiv[1] << " "<<ipiv[2] << " "<<endl;

  detM = 1;
  for(unsigned int i = 0; i < m; i ++){
    if(ipiv[i] == i + 1){
      detM *= M[i * m + i];
    }
    else{
      detM *= -M[i * m + i];
    }
  }
  delete []ipiv;

  ipiv = new int[m1];
  assert(ipiv != NULL);
  sgetrf_(m1, m1, K, m1, ipiv, info);

  //cout<<ipiv[0] << " "<<ipiv[1] << " "<<ipiv[2] << " "<< ipiv[3] << endl;

  detK = 1;
  for(unsigned int i = 0; i < m1; i ++){
    if(ipiv[i] == i + 1){
      detK *= K[i * m1 + i];
    }
    else{
      detK *= -K[i * m1 + i];
    }
  }
  delete []ipiv;

  C_1 = -8.0 * detM;
  C_0 = detK;

  //--------------------------------------------------
  // if(debug){
  //  cout<<"C_2: "<<C_2<<" C_1: "<<C_1<<" C_0: "<<C_0<<endl;
  // }
  //-------------------------------------------------- 

  float lambda_old, lambda2, a, b;
  float lambda = (G_x + G_y) / 2.0;
  
  unsigned int maxits = 50;
  float tolerance = 1.0e-6;

  for(unsigned int i = 0; i < maxits; i ++){
    lambda_old = lambda;
    lambda2 = lambda_old*lambda_old;
    b = (lambda2 + C_2)*lambda_old;
    a = b + C_1;
    lambda = lambda_old - (a*lambda_old + C_0) / (2.0*lambda2*lambda_old + b + a);
    
    //--------------------------------------------------
    // if(debug){
    //  cout<<"i "<<i<<" lambda "<<lambda<<" a: "<<a<<" b: "<<b<<endl;
    // }
    //-------------------------------------------------- 

    if (fabs(lambda - lambda_old) < fabs(tolerance*lambda)){
      break;
    }
    
  }

   float rmsd2 = (G_x + G_y - 2.0 * lambda) / natom;
  float ls_rmsd = 0.0;
  if(rmsd2 > 0){
    ls_rmsd = sqrt(rmsd2);
  }

  delete []K;
  delete []M;
  delete []mcA;
  delete []mcB;
  
  return ls_rmsd;
}

