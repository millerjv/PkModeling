/*=auto=========================================================================

  Portions (c) Copyright 2009 Brigham and Women's Hospital (BWH) All Rights Reserved.

  See Doc/copyright/copyright.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Program:   3D Slicer
  Module:    $RCSfile: pk_solver.cxx,v $
  Date:      $Date: 2006/03/19 17:12:29 $
  Version:   $Revision: 1.13 $

=========================================================================auto=*/

#include <vnl/algo/vnl_convolve.h>
#include <itkGradientMagnitudeImageFilter.h>
#include <itkImageRegionIterator.h>
#include <itkLevenbergMarquardtOptimizer.h>
#include "PkSolver.h"
#include "itkTimeProbesCollectorBase.h"
#include <string>

namespace itk
{	
//
// Support routines/classes used internally in the PkSolver
// 
static itk::TimeProbesCollectorBase probe;

int m_ConstantBAT;
std::string m_BATCalculationMode;

//
// Implementation of the PkSolver API
// 
//

bool pk_solver (int signalSize, const float* timeAxis, 
                const float* PixelConcentrationCurve, 
                const float* BloodConcentrationCurve, 
                float& Ktrans, float& Ve, float& Fpv,
                float fTol, float gTol, float xTol,
                float epsilon, int maxIter,
                float hematocrit,
                int modelType,
                int constantBAT,
                const std::string BATCalculationMode)
{
  // Note the unit: timeAxis should be in minutes!! This could be related to the following parameters!!
  // fTol      =  1e-4;  // Function value tolerance
  // gTol      =  1e-4;  // Gradient magnitude tolerance 
  // xTol      =  1e-5;  // Search space tolerance
  // epsilon   =  1e-9;    // Step
  // maxIter   =   200;  // Maximum number of iterations

  m_BATCalculationMode = BATCalculationMode;
  m_ConstantBAT = constantBAT;
  
  // Levenberg Marquardt optimizer  
  itk::LevenbergMarquardtOptimizer::Pointer  optimizer = itk::LevenbergMarquardtOptimizer::New();
  LMCostFunction::Pointer costFunction = LMCostFunction::New();

  LMCostFunction::ParametersType initialValue;
  if(modelType == itk::LMCostFunction::TOFTS_2_PARAMETER)
    {
    initialValue = LMCostFunction::ParametersType(2); ///...
    }
  else
    {
    initialValue = LMCostFunction::ParametersType(3);
    initialValue[2] = 0.1;     //f_pv //...
    }
  initialValue[0] = 0.1;     //Ktrans //...
  initialValue[1] = 0.5;     //ve //...
 
  costFunction->SetNumberOfValues (signalSize);

  costFunction->SetCb (BloodConcentrationCurve, signalSize); //BloodConcentrationCurve
  costFunction->SetCv (PixelConcentrationCurve, signalSize); //Signal Y
  costFunction->SetTime (timeAxis, signalSize); //Signal X
  costFunction->SetHematocrit (hematocrit);
  costFunction->GetValue (initialValue);
  costFunction->SetModelType(modelType);

  try {
    optimizer->SetCostFunction( costFunction.GetPointer() );
  }
  catch ( itk::ExceptionObject & e ) {
    std::cout << "Exception thrown ! " << std::endl;
    std::cout << "An error ocurred during Optimization" << std::endl;
    std::cout << e << std::endl;
    return false;
  }

  // this following call is equivalent to invoke: costFunction->SetUseGradient( useGradient );
  optimizer->UseCostFunctionGradientOff();
  optimizer->SetUseCostFunctionGradient(0);

  itk::LevenbergMarquardtOptimizer::InternalOptimizerType * vnlOptimizer = optimizer->GetOptimizer();

  vnlOptimizer->set_f_tolerance( fTol );
  vnlOptimizer->set_g_tolerance( gTol );
  vnlOptimizer->set_x_tolerance( xTol ); 
  vnlOptimizer->set_epsilon_function( epsilon );
  vnlOptimizer->set_max_function_evals( maxIter );

  // We start not so far from the solution 

  optimizer->SetInitialPosition( initialValue );

  CommandIterationUpdateLevenbergMarquardt::Pointer observer = 
    CommandIterationUpdateLevenbergMarquardt::New();
  optimizer->AddObserver( itk::IterationEvent(), observer );
  optimizer->AddObserver( itk::FunctionEvaluationIterationEvent(), observer );

  try {
    optimizer->StartOptimization();
  }
  catch( itk::ExceptionObject & e ) {
    std::cerr << "Exception thrown ! " << std::endl;
    std::cerr << "An error ocurred during Optimization" << std::endl;
    std::cerr << "Location    = " << e.GetLocation()    << std::endl;
    std::cerr << "Description = " << e.GetDescription() << std::endl;
    return false;
  }

  itk::LevenbergMarquardtOptimizer::ParametersType finalPosition;
  finalPosition = optimizer->GetCurrentPosition();

  //Solution: remove the scale of 100  
  Ktrans = finalPosition[0];
  Ve = finalPosition[1];
  if(modelType == itk::LMCostFunction::TOFTS_3_PARAMETER)
    {
    Fpv = finalPosition[2];  
    }

  return true;
}

bool pk_solver(int signalSize, const float* timeAxis, 
               const float* PixelConcentrationCurve, 
               const float* BloodConcentrationCurve, 
               float& Ktrans, float& Ve, float& Fpv,
               float fTol, float gTol, float xTol,
               float epsilon, int maxIter,
               float hematocrit,
               itk::LevenbergMarquardtOptimizer* optimizer,
               LMCostFunction* costFunction,
               int modelType,
               int constantBAT,
               const std::string BATCalculationMode
               )
{
  //std::cout << "in pk solver" << std::endl;
  // probe.Start("pk_solver");
  // Note the unit: timeAxis should be in minutes!! This could be related to the following parameters!!
  // fTol      =  1e-4;  // Function value tolerance
  // gTol      =  1e-4;  // Gradient magnitude tolerance 
  // xTol      =  1e-5;  // Search space tolerance
  // epsilon   =  1e-9;    // Step
  // maxIter   =   200;  // Maximum number of iterations
  //std::cerr << "In pkSolver!" << std::endl;

  m_BATCalculationMode = BATCalculationMode;
  m_ConstantBAT = constantBAT;

  // Levenberg Marquardt optimizer  
        
  //////////////
  LMCostFunction::ParametersType initialValue;
  if(modelType == itk::LMCostFunction::TOFTS_2_PARAMETER)
    {
    initialValue = LMCostFunction::ParametersType(2); ///...
    }
  else
    {
    initialValue = LMCostFunction::ParametersType(3);
    initialValue[2] = 0.1;     //f_pv //...
    }
  initialValue[0] = 0.1;     //Ktrans //...
  initialValue[1] = 0.5;     //ve //...
        
  costFunction->SetNumberOfValues (signalSize);
  

  costFunction->SetCb (BloodConcentrationCurve, signalSize); //BloodConcentrationCurve
  costFunction->SetCv (PixelConcentrationCurve, signalSize); //Signal Y
  costFunction->SetTime (timeAxis, signalSize); //Signal X
  costFunction->SetHematocrit (hematocrit);
  costFunction->GetValue (initialValue); //...
  costFunction->SetModelType(modelType);

  optimizer->UseCostFunctionGradientOff();   

  try {
     optimizer->SetCostFunction( costFunction ); 
  }
  catch ( itk::ExceptionObject & e ) {
  std::cout << "Exception thrown ! " << std::endl;
  std::cout << "An error ocurred during Optimization" << std::endl;
  std::cout << e << std::endl;
  return false;
  }   
        
  itk::LevenbergMarquardtOptimizer::InternalOptimizerType * vnlOptimizer = optimizer->GetOptimizer();//...

  vnlOptimizer->set_f_tolerance( fTol ); //...
  vnlOptimizer->set_g_tolerance( gTol ); //...
  vnlOptimizer->set_x_tolerance( xTol ); //...
  vnlOptimizer->set_epsilon_function( epsilon ); //...
  vnlOptimizer->set_max_function_evals( maxIter ); //...
        
  // We start not so far from the solution 
        
  optimizer->SetInitialPosition( initialValue ); //...       
  
  try {
  //  probe.Start("optimizer");
  optimizer->StartOptimization();
  //   probe.Stop("optimizer");
  }
  catch( itk::ExceptionObject & e ) {
  std::cerr << "Exception thrown ! " << std::endl;
  std::cerr << "An error ocurred during Optimization" << std::endl;
  std::cerr << "Location    = " << e.GetLocation()    << std::endl;
  std::cerr << "Description = " << e.GetDescription() << std::endl;
  return false;
  }
  //vnlOptimizer->diagnose_outcome();
  //std::cerr << "after optimizer!" << std::endl;
  itk::LevenbergMarquardtOptimizer::ParametersType finalPosition;
  finalPosition = optimizer->GetCurrentPosition();
  //std::cerr << finalPosition[0] << ", " << finalPosition[1] << ", " << finalPosition[2] << std::endl;

        
  //Solution: remove the scale of 100  
  Ktrans = finalPosition[0];
  Ve = finalPosition[1];
  if(modelType == itk::LMCostFunction::TOFTS_3_PARAMETER)
    {
    Fpv = finalPosition[2];  
    }

  // "Project" back onto the feasible set.  Should really be done as a
  // constraint in the optimization.
  if(Ve<0) Ve = 0;
  if(Ve>1) Ve = 1;
  if(Ktrans<0) Ktrans = 0;
  if(Ktrans>5) Ktrans = 5;
		
  //if((Fpv>1)||(Fpv<0)) Fpv = 0;
  //  probe.Stop("pk_solver");
  return true;
}

void pk_report()
{
     probe.Report();
}

void pk_clear()
{
    probe.Clear();
}

#define PI 3.1415926535897932384626433832795
#define IS_NAN(x) ((x) != (x))

bool convert_signal_to_concentration (unsigned int signalSize, 
                                      const float* SignalIntensityCurve, 
                                      const float T1Pre, float TR, float FA,
                                      float* concentration,
                                      float RGd_relaxivity,
                                      float s0,
                                      float S0GradThresh)
{
  const double exp_TR_BloodT1 = exp (-TR/T1Pre);
  const float alpha = FA * PI/180;  
  const double cos_alpha = cos(alpha);
  const double constB = (1-exp_TR_BloodT1) / (1-cos_alpha*exp_TR_BloodT1);
  
  if (s0 == -1.0f)
    s0 = compute_s0_individual_curve (signalSize, SignalIntensityCurve, S0GradThresh, m_BATCalculationMode, m_ConstantBAT);
            
  for (unsigned int t=0; t<signalSize; ++t) {
    const float tSignal = SignalIntensityCurve[t];
    if (tSignal == 0) {
      concentration[t] = 0;
      continue;
    }

    const double constA = tSignal/s0;
    double value = (1 - constA * constB) / (1- constA * constB * cos_alpha);
    double log_value = log(value);
    if (IS_NAN(log_value)) {
      concentration[t] = 0;
    }
    else {
      const float ROft = (-1/TR) * log_value;
      if (T1Pre != 0) {
        double Cb = (ROft - (1/T1Pre)) / RGd_relaxivity;
        assert (IS_NAN (Cb) == false); 
        if (Cb < 0)
          Cb = 0;
        concentration[t] = float (Cb);
      }  
      else
        concentration[t] = 0;
    }
  }

  return true;
}

float area_under_curve(int signalSize, 
                       const float* timeAxis,
                       const float* concentration, 
                       int BATIndex, 
                       float aucTimeInterval)
{	
  float auc = 0.0f;
  if(BATIndex>=signalSize) return auc;
  //std::cerr << std::endl << "BATIndex:"<<BATIndex << std::endl;	
	
  int lastIndex = BATIndex;
  //std::cerr << std::endl << "timeAxis[0]:"<< timeAxis[0] << std::endl;	
  float targetTime = timeAxis[BATIndex]+aucTimeInterval; 
  //std::cerr << std::endl << "targetTime:"<< targetTime << std::endl;	
  float tempTime = timeAxis[BATIndex+1];
  //std::cerr << std::endl << "tempTime:"<< tempTime << std::endl;	

  //find the last index
  while((tempTime<targetTime)&&(lastIndex<(signalSize-2)))
    {	
    lastIndex+=1;		
    tempTime = timeAxis[lastIndex+1] ;
    //std::cerr << std::endl << "tempTime"<<tempTime << std::endl;
    }	
		
  if ((lastIndex-BATIndex)==0) return auc = aucTimeInterval*concentration[BATIndex];

  //extract time and concentration
  float * concentrationValues = new float[lastIndex-BATIndex+2]();
  float * timeValues = new float[lastIndex-BATIndex+2]();

  //find the extra time and concentration value for auc
  float y1, y2, x1, x2, slope, b, targetX, targetY;
  y2 = concentration[lastIndex+1];
  y1 = concentration[lastIndex];
  x2 = timeAxis[lastIndex+1];
  x1 = timeAxis[lastIndex];
  slope = (y2-y1)/(x2-x1);
  b = y1-slope*x1;
  targetX = timeAxis[BATIndex] + aucTimeInterval;
  targetY = slope*targetX+b;
  if(targetX>timeAxis[signalSize-1])
    {
    targetX = timeAxis[lastIndex+1];
    targetY = concentration[lastIndex+1];
    }
  concentrationValues[lastIndex-BATIndex+1] = targetY; //put the extra value at the end
  timeValues[lastIndex-BATIndex+1] = targetX; //put the extra time value at the end
	
//	printf("lastIndex is %f\n", (float)lastIndex);
  for(int i=0; i<(lastIndex-BATIndex+1); ++i)
    {
    concentrationValues[i] = concentration[i+BATIndex];
    timeValues[i] = timeAxis[i+BATIndex];
    //printf("lastIndex is %f,%f\n", (float) concentrationValues[i],(float)timeValues[i]);
    }	

  //get auc
  auc = intergrate(concentrationValues,timeValues,(lastIndex-BATIndex+2));

  delete [] concentrationValues;
  delete [] timeValues;
  return auc;
}

float intergrate(float* yValues, float * xValues, int size)
{
  float area= 0.0f;
  for(int i=1; i<size;++i)
    {
    area+=(xValues[i]-xValues[i-1])*(yValues[i]+yValues[i-1])/2;
    //	std::cerr << std::endl << "area:" << area<<","<<yValues[i] << std::endl;
    }
  return area;
}

void compute_derivative (int signalSize,
                         const float* SignalY,
                         float* YDeriv)
{
  YDeriv[0] = (float) ((-3.0*SignalY[0] + 4.0*SignalY[1] - SignalY[2]) / 2.0);
  YDeriv[signalSize-1] = (float) ((3.0*SignalY[signalSize-1] - 4.0*SignalY[signalSize-2] + SignalY[signalSize-3]) / 2.0);
  for(int i=1; i<signalSize-1; i++) {
    YDeriv[i] = (float) ((SignalY[i+1] - SignalY[i-1]) / 2.0);
  }
}

void compute_derivative_forward (int signalSize,
                         const float* SignalY,
                         float* YDeriv)
{  
  YDeriv[signalSize-1] = (float) (SignalY[signalSize-1] - SignalY[signalSize-2]);
  for(int i=0; i<signalSize-1; i++) {
    YDeriv[i] = (float) (SignalY[i+1] - SignalY[i]);
  }
}

void compute_derivative_backward (int signalSize,
                         const float* SignalY,
                         float* YDeriv)
{
  YDeriv[0] = (float) (-SignalY[0] + SignalY[1]);  
  for(int i=1; i<signalSize; i++) {
    YDeriv[i] = (float) (SignalY[i] - SignalY[i-1]);
  }
}

float get_signal_max (int signalSize, const float* SignalY, int& index)
{
  float max = -1E10f;
  index = -1;
  for (int i=0; i<signalSize; i++)
    if (SignalY[i] > max)
      {
      max = SignalY[i];
      index = i;
      }
  return max;
}

bool compute_bolus_arrival_time (int signalSize, const float* SignalY,
                                int& ArrivalTime, int& FirstPeak, float& MaxSlope)
{
  float* y0 = new float[signalSize];  // Input pixel values
  int i=0;
  for(  i = 0; i < signalSize; i++ ) //{
    y0[i] = SignalY[i];

  int skip1 = 0;             // Leading points to ignore
  int skip2 = 1;             // Trailing points to ignore
  //int* t = new int[signalSize];       // time value
  float* yd = new float[signalSize];  // working buffer

  int CpIndex = 0;
  float Cp = get_signal_max (signalSize, SignalY, CpIndex); //this->m_TimeSeriesY->max_value();

  // Detect ArrivalTime Detection failure and report indeterminate results.
  if( CpIndex < 0) {
    ArrivalTime = 0;
    FirstPeak = 0;
    MaxSlope = 0;
    delete [] y0;
    //delete [] t;
    delete [] yd;
    return false;
  }
  
  // Step 1: Smoothing done using Savizky-Golay before this call on Signal or Conc. data

  // Step 2: Spatial derivative of smoothed data

  memcpy(yd, y0, signalSize*sizeof(float));
  ///this->ComputeDerivative(yd);
  compute_derivative (signalSize, y0, yd);

  // Step 3: Find point of steepest descent/ascent
  //int min_index = skip1;

  //added this with Sandeep's suggestions
  int max_index = skip1;
  float max = yd[max_index];

  for (i=max_index; i<signalSize-skip2-1; i++) {
    if(yd[i] > max) {
      max = yd[i];
      max_index = i;
    }
  }
  MaxSlope = max;
  float thresh = (float)( max / 10.0 );

  // Step 4: Arrival Time detection
  for( i=max_index; i>=skip1; i--) {
    if(yd[i] < thresh)
      break;
  }

  ArrivalTime = i+1;

  // Step 5: Peak Time detection
  //for( i=min_index; i<signalSize-1-skip2; i++) {
  /// jvm - removing this loop as it is not used
  // for( i=skip1; i<signalSize-1-skip2; i++) {
  //   if(yd[i] >= thresh || y0[i] >= y0[i-1])
  //     break;
  // }
  // FirstPeak = i; 
  // jvm - end of remove
  //changing the peak as global peak 
  FirstPeak = CpIndex;

  //delete [] t;
  delete [] yd;
  delete [] y0;
  return true;
}

void compute_gradient_old (int signalSize, const float* SignalY, float* SignalGradient)
{
  typedef itk::Image<float, 1>   ImageType;
  typedef itk::GradientMagnitudeImageFilter<ImageType, ImageType>  FilterType;

  ImageType::Pointer SignalImg1D = ImageType::New();  
  ImageType::SizeType imgSize;
  imgSize[0] = signalSize;
  SignalImg1D->SetRegions(imgSize);
  SignalImg1D->Allocate();

  typedef itk::ImageRegionIterator< ImageType > IteratorType;
  IteratorType inputIt (SignalImg1D, SignalImg1D->GetLargestPossibleRegion());
  inputIt.GoToBegin();
  for (unsigned int i=0; !inputIt.IsAtEnd(); i++, ++inputIt) {
    inputIt.Set (SignalY[i]);
  }

  FilterType::Pointer gfilter = FilterType::New();
  gfilter->SetInput (SignalImg1D);
  gfilter->Update();

  ImageType::Pointer GradientImg1D = gfilter->GetOutput();
  
  IteratorType outputIt (GradientImg1D, GradientImg1D->GetLargestPossibleRegion());
  outputIt.GoToBegin();
  for (unsigned int i=0; !outputIt.IsAtEnd(); i++, ++outputIt) {
    SignalGradient[i] = outputIt.Get();	
  }
}

void compute_gradient (int signalSize, const float* SignalY, float* SignalGradient)
{	
  compute_derivative (signalSize, SignalY, SignalGradient);
  for(int i=0; i<signalSize; i++) 
    {
    SignalGradient[i] = sqrt(SignalGradient[i]*SignalGradient[i]);		
    } 
}

void compute_gradient_forward (int signalSize, const float* SignalY, float* SignalGradient)
{
  compute_derivative_forward (signalSize, SignalY, SignalGradient);
  for(int i=0; i<signalSize; i++) 
    {
    SignalGradient[i] = sqrt(SignalGradient[i]*SignalGradient[i]);		
    } 	
}

void compute_gradient_backward (int signalSize, const float* SignalY, float* SignalGradient)
{	
  compute_derivative_backward (signalSize, SignalY, SignalGradient);
  for(int i=0; i<signalSize; i++) 
    {
    SignalGradient[i] = sqrt(SignalGradient[i]*SignalGradient[i]);		
    } 	
}


float compute_s0_using_sumsignal_properties (int signalSize, const float* SignalY, 
                                             const short* lowGradIndex, int FirstPeak)
{
  double S0 = 0;
  int count = 0;
  double sum = 0;
  int first = FirstPeak;
  if (FirstPeak > signalSize)
    {
    first = signalSize;
    }
  for (int i=0; i<first; i++) {
    sum += SignalY[i];
    if (lowGradIndex[i]==1) {
      S0 += SignalY[i];
      count++;
    }
  }
  if (count)
    S0 /= count;
  else
    S0 = sum / (first-1);
  return float (S0);
}

float compute_s0_individual_curve (int signalSize, const float* SignalY, float S0GradThresh,
                                   std::string BATCalculationMode, int constantBAT)
{  
  double S0 = 0;
  int ArrivalTime, FirstPeak;
  float MaxSlope;
  bool result;

  if (BATCalculationMode == "UseConstantBAT")
  {
  // Use constant BAT
    
  ArrivalTime = constantBAT;
  result = true;
  }
  else if (BATCalculationMode == "PeakGradient")
  {
    result = compute_bolus_arrival_time (signalSize, SignalY, ArrivalTime, FirstPeak, MaxSlope);//same
  }

  if (result == false) {
    ///printf ("  Compute compute_s0_individual_curve fails! S0 = 0.\n");
    return 0;
  }
  
  float* SignalGradient = new float[signalSize];
  //above: same
  compute_gradient(signalSize, SignalY, SignalGradient);    
  
  int count = 0;
  double sum = 0;
  for (int i=0; i<ArrivalTime; i++) { //updated to i<ArrivalTime by Yingxuan Zhu on 5/5/2012, original i<FirstPeak;	  
    sum += SignalY[i];
    if (SignalGradient[i] < S0GradThresh) {
      S0 += SignalY[i];
      count++;
    }
	//std::cerr<<"sum, S0:"<<sum<<","<<S0<<"\n"<<std::endl;
  }
  if(ArrivalTime>0)
  {
	if (count)
		S0 /= count;
	else
		S0 = sum / (ArrivalTime); //updated to sum/(ArrivalTime) by Yingxuan Zhu on 5/5/2012, original sum/(FirstPeak-1);	
  }
  else
	  S0 = SignalY[0]; //ArrivalTime is 0;

  delete [] SignalGradient;
  return float (S0);

}

}; // end of namespace

