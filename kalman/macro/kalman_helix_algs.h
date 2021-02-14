#ifndef GARSOFT_RECO_TRACKER2ALGS_H
#define GARSOFT_RECO_TRACKER2ALGS_H

#include <vector>
#include <iostream>
#include <cstddef>
#include <TMath.h>
#include <TVector3.h>
#include "TTree.h"

float capprox(float x1,float y1,
                  float x2,float y2,
                  float x3,float y3,
                  float &xc, float &yc);

void ouchef(double *x, double *y, int np, double &xcirccent, double &ycirccent,
                double &rcirc, double &chisq, int &ifail);

int initial_trackpar_estimate(TTree &t1s,
                            float &x,
                            float &y,
                            float &z,
                            size_t n,
                            float &curvature_init,
                            float &lambda_init,
                            float &phi_init,
                            float &xpos,
                            float &ypos,
                            float &zpos,
                            int printlevel);

#endif