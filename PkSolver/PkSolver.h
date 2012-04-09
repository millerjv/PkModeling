/*=auto=========================================================================

  Portions (c) Copyright 2009 Brigham and Women's Hospital (BWH) All Rights Reserved.

  See Doc/copyright/copyright.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   3D Slicer
  Module:    $RCSfile: pk_solver.h,v $
  Date:      $Date: 2006/03/19 17:12:29 $
  Version:   $Revision: 1.13 $

=========================================================================auto=*/

#ifndef PkSolver_h_
#define PkSolver_h_

//#include "PkSolver.h"

namespace itk
{

bool pk_solver(const int signalSize, const float* timeAxis, 
                const float* PixelConcentrationCurve, 
                const float* BloodConcentrationCurve, 
                float& Ktrans, float& Ve, float& Fpv,
                const float fTol = 1e-4f, 
                const float gTol = 1e-4f, 
                const float xTol = 1e-5f,
                const float epsilon = 1e-9f, 
                const int maxIter = 200,
                const float hematocrit = 0.4f);


bool convert_signal_to_concentration (const unsigned int signalSize, 
                                      const float* SignalIntensityCurve, 
                                      const float T1, const float TR, const float FA,
                                      float*& concentration,
                                      const float relaxivity = 4.9E-3f,
                                      float s0 = -1.0f,
                                      const float S0GradThresh = 15.0f);


void compute_derivative (const int signalSize, const float* SingnalY, float*& YDeriv);
float get_signal_max (const int signalSize, const float* SignalY);

bool compute_bolus_arrival_time (const int signalSize, const float* SignalY,
                                int& ArrivalTime, int& FirstPeak, float& MaxSlope);

void compute_gradient (const int signalSize, const float* SignalY, float*& SignalGradient);

float compute_s0_using_sumsignal_properties (const int signalSize, const float* SignalY, 
                                             const short* lowGradIndex, const int FirstPeak);

float compute_s0_individual_curve (const int signalSize, const float* SignalY, const float S0GradThresh);

};

#endif
