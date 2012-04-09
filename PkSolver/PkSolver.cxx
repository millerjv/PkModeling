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
#include <math.h>

namespace itk
{	
//
// Support routines/classes used internally in the PkSolver
// 

class LMCostFunction: public itk::MultipleValuedCostFunction
{
public:
  typedef LMCostFunction                    Self;
  typedef itk::MultipleValuedCostFunction   Superclass;
  typedef itk::SmartPointer<Self>           Pointer;
  typedef itk::SmartPointer<const Self>     ConstPointer;
  itkNewMacro( Self );

  enum { SpaceDimension =  3 };
  unsigned int RangeDimension; 

  typedef Superclass::ParametersType              ParametersType;
  typedef Superclass::DerivativeType              DerivativeType;
  typedef Superclass::MeasureType                 MeasureType, ArrayType;
  typedef Superclass::ParametersValueType         ValueType;

  float m_hematocrit;

  LMCostFunction():
  m_Measure(300)
  {
  }

  void set_hematocrit (const float hematocrit) {
    m_hematocrit = hematocrit;
  }

  void SetNumberOfValues(const unsigned int NumberOfValues)
  {
    RangeDimension = NumberOfValues;
  }

  void SetCb (const float* cb, int sz) //BloodConcentrationCurve.
  {
    double *tmp;
    tmp = new double[sz];
    Cb.set_size(sz);
    for( int i = 0; i < sz; ++i )
      tmp[i] = cb[i];
    Cb.set(tmp);
    delete [] tmp;
  }


  void SetCv (const float* cv, int sz) //Self signal Y
  {    
    double *tmp;
    tmp = new double[sz];
    Cv.set_size (sz);
    for (int i = 0; i < sz; ++i)
      tmp[i] = cv[i];
    Cv.set(tmp);
    delete [] tmp;
  }
 
  void SetTime (const float* cx, int sz) //Self signal X
  {
    double *tmp;
    tmp = new double[sz];
    Time.set_size (sz);
    for( int i = 0; i < sz; ++i )
      tmp[i] = cx[i];
    Time.set(tmp);
    delete [] tmp;
  }

  MeasureType GetValue( const ParametersType & parameters) const
  {
    m_Measure.SetSize(RangeDimension);
    ValueType Ktrans = parameters[0];
    ValueType Ve = parameters[1];
    ValueType f_pv = parameters[2];
    ///const ValueType hematocrit = 0.4; //go to XML todo
    //const ValueType w = 0.001;

    ArrayType VeTerm;
    VeTerm = -Ktrans/Ve*Time;
    ValueType deltaT = Time(1) - Time(0);
    m_Measure = Cv - (1/(1.0-m_hematocrit)*(Ktrans*deltaT*Convolution(Cb,Exponential(VeTerm)) + f_pv*Cb));

    return m_Measure; 
  }

  //Not going to be used
  void GetDerivative( const ParametersType & parameters,
    DerivativeType  & derivative ) const
  {   
  }

  unsigned int GetNumberOfParameters(void) const
  {
    return SpaceDimension;
  }

  unsigned int GetNumberOfValues(void) const
  {
    return RangeDimension;
  }


private:

  mutable MeasureType       m_Measure;
  mutable DerivativeType    m_Derivative;

  ArrayType Cv, Cb, Time;

  ArrayType Convolution(ArrayType X, ArrayType Y) const
  {
    ArrayType Z;
    Z.set_size(X.size());
    ArrayType temp;
    temp = vnl_convolve(X,Y);
    Z = temp.extract(X.size(),0);
    return Z;
  };

  ArrayType Exponential(ArrayType X) const
  {
    ArrayType Z;
    Z.set_size(X.size());
    for (unsigned int i=0; i<X.size(); i++)
    {
      Z[i] = exp(X(i));
    }
    return Z;
  };

  int constraintFunc(ValueType x) const
  {
    if (x<0||x>1)
      return 1;
    else
      return 0;

  };


};

class CommandIterationUpdateLevenbergMarquardt : public itk::Command 
{
public:
  typedef  CommandIterationUpdateLevenbergMarquardt   Self;
  typedef  itk::Command                               Superclass;
  typedef itk::SmartPointer<Self>                     Pointer;
  itkNewMacro( Self );
protected:
  CommandIterationUpdateLevenbergMarquardt() 
  {
    m_IterationNumber=0;
  }
public:
  typedef itk::LevenbergMarquardtOptimizer   OptimizerType;
  typedef   const OptimizerType   *          OptimizerPointer;

  void Execute(itk::Object *caller, const itk::EventObject & event)
  {
    Execute( (const itk::Object *)caller, event);
  }

  void Execute(const itk::Object * object, const itk::EventObject & event)
  {
    //std::cout << "Observer::Execute() " << std::endl;
    OptimizerPointer optimizer = 
      dynamic_cast< OptimizerPointer >( object );
    if( m_FunctionEvent.CheckEvent( &event ) )
    {
     // std::cout << m_IterationNumber++ << "   ";
     // std::cout << optimizer->GetCachedValue() << "   ";
     // std::cout << optimizer->GetCachedCurrentPosition() << std::endl;
    }
    else if( m_GradientEvent.CheckEvent( &event ) )
    {
      std::cout << "Gradient " << optimizer->GetCachedDerivative() << "   ";
    }

  }
private:
  unsigned long m_IterationNumber;

  itk::FunctionEvaluationIterationEvent m_FunctionEvent;
  itk::GradientEvaluationIterationEvent m_GradientEvent;
};


//
// Implementation of the PkSolver API
// 
//

bool pk_solver (const int signalSize, const float* timeAxis, 
                const float* PixelConcentrationCurve, 
                const float* BloodConcentrationCurve, 
                float& Ktrans, float& Ve, float& Fpv,
                const float fTol, const float gTol, const float xTol,
                const float epsilon, const int maxIter,
                const float hematocrit)
{
  // Note the unit: timeAxis should be in minutes!! This could be related to the following parameters!!
  // fTol      =  1e-4;  // Function value tolerance
  // gTol      =  1e-4;  // Gradient magnitude tolerance 
  // xTol      =  1e-5;  // Search space tolerance
  // epsilon   =  1e-9;    // Step
  // maxIter   =   200;  // Maximum number of iterations

  // Levenberg Marquardt optimizer  
  itk::LevenbergMarquardtOptimizer::Pointer  optimizer = itk::LevenbergMarquardtOptimizer::New();
  LMCostFunction::Pointer costFunction = LMCostFunction::New();

  LMCostFunction::ParametersType initialValue(LMCostFunction::SpaceDimension);
  initialValue[0] = 0.1;     //Ktrans
  initialValue[1] = 0.5;     //ve 
  initialValue[2] = 0.1;     //f_pv

  costFunction->SetNumberOfValues (signalSize);

  costFunction->SetCb (BloodConcentrationCurve, signalSize); //BloodConcentrationCurve
  costFunction->SetCv (PixelConcentrationCurve, signalSize); //Signal Y
  costFunction->SetTime (timeAxis, signalSize); //Signal X
  costFunction->set_hematocrit (hematocrit);
  costFunction->GetValue (initialValue);

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
  Fpv = finalPosition[2];  

  return true;
}

#define PI 3.1415926535897932384626433832795
#define IS_NAN(x) ((x) != (x))

bool convert_signal_to_concentration (const unsigned int signalSize, 
                                      const float* SignalIntensityCurve, 
                                      const float T1Pre, const float TR, const float FA,
                                      float*& concentration,
                                      const float RGd_relaxivity,
                                      float s0,
                                      const float S0GradThresh)
{
  const double exp_TR_BloodT1 = exp (-TR/T1Pre);
  const float alpha = FA * PI/180;  
  const double cos_alpha = cos(alpha);
  const double constB = (1-exp_TR_BloodT1) / (1-cos_alpha*exp_TR_BloodT1);
  
  if (s0 == -1.0f)
    s0 = compute_s0_individual_curve (signalSize, SignalIntensityCurve, S0GradThresh);
            
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

void compute_derivative (const int signalSize,
                         const float* SignalY,
                         float*& YDeriv)
{
  YDeriv[0] = (float) ((-3.0*SignalY[0] + 4.0*SignalY[1] - SignalY[2]) / 2.0);
  YDeriv[signalSize-1] = (float) ((3.0*SignalY[signalSize-1] - 4.0*SignalY[signalSize-2] + SignalY[signalSize-3]) / 2.0);
  for(int i=1; i<signalSize-1; i++) {
    YDeriv[i] = (float) ((SignalY[i+1] - SignalY[i-1]) / 2.0);
  }
}

void compute_derivative_forward (const int signalSize,
                         const float* SignalY,
                         float*& YDeriv)
{  
  YDeriv[signalSize-1] = (float) (SignalY[signalSize-1] - SignalY[signalSize-2]);
  for(int i=0; i<signalSize-1; i++) {
    YDeriv[i] = (float) (SignalY[i+1] - SignalY[i]);
  }
}

void compute_derivative_backward (const int signalSize,
                         const float* SignalY,
                         float*& YDeriv)
{
  YDeriv[0] = (float) (-SignalY[0] + SignalY[1]);  
  for(int i=1; i<signalSize; i++) {
    YDeriv[i] = (float) (SignalY[i] - SignalY[i-1]);
  }
}

float get_signal_max (const int signalSize, const float* SignalY)
{
  float max = -1E10f;
  for (int i=0; i<signalSize; i++)
    if (SignalY[i] > max)
      max = SignalY[i];
  return max;
}

bool compute_bolus_arrival_time (const int signalSize, const float* SignalY,
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

  float Cp = get_signal_max (signalSize, SignalY); //this->m_TimeSeriesY->max_value();
  int CpIndex = 0;

  // Find index of Cp
  // !!!Note the indention here. Any problem??
  for(  i = 0; i < signalSize; ++i ) {
    if( SignalY[i] == Cp ) {
      CpIndex = i;
      break;
    }
  }

  // Correct skip1 for cases with Cp less than skip1
  if( CpIndex <= skip1 )
    skip1 = CpIndex - 2;

  // Detect ArrivalTime Detection failure and report indeterminate results.
  if( skip1 < 0) {
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
  for( i=skip1; i<signalSize-1-skip2; i++) {
    if(yd[i] >= thresh || y0[i] >= y0[i-1])
      break;
  }
  FirstPeak = i; 
  //changing the peak as global peak 
  FirstPeak = CpIndex;

  //delete [] t;
  delete [] yd;
  delete [] y0;
  return true;
}

void compute_gradient_old (const int signalSize, const float* SignalY, float*& SignalGradient)
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

void compute_gradient (const int signalSize, const float* SignalY, float*& SignalGradient)
{	
	compute_derivative (signalSize, SignalY, SignalGradient);
	for(int i=0; i<signalSize; i++) 
	{
		SignalGradient[i] = sqrt(SignalGradient[i]*SignalGradient[i]);		
	} 
}

void compute_gradient_forward (const int signalSize, const float* SignalY, float*& SignalGradient)
{
	compute_derivative_forward (signalSize, SignalY, SignalGradient);
	for(int i=0; i<signalSize; i++) 
	{
		SignalGradient[i] = sqrt(SignalGradient[i]*SignalGradient[i]);		
	} 	
}

void compute_gradient_backward (const int signalSize, const float* SignalY, float*& SignalGradient)
{	
	compute_derivative_backward (signalSize, SignalY, SignalGradient);
	for(int i=0; i<signalSize; i++) 
	{
		SignalGradient[i] = sqrt(SignalGradient[i]*SignalGradient[i]);		
	} 	
}


float compute_s0_using_sumsignal_properties (const int signalSize, const float* SignalY, 
                                             const short* lowGradIndex, const int FirstPeak)
{
  double S0 = 0;
  int count = 0;
  double sum = 0;
  for (int i=0; i<FirstPeak; i++) {
    sum += SignalY[i];
    if (lowGradIndex[i]==1) {
      S0 += SignalY[i];
      count++;
    }
  }
  if (count)
    S0 /= count;
  else
    S0 = sum / (FirstPeak-1);
  return float (S0);
}

float compute_s0_individual_curve (const int signalSize, const float* SignalY, const float S0GradThresh)
{  
  double S0 = 0;
  int ArrivalTime, FirstPeak;
  float MaxSlope;
  bool result = compute_bolus_arrival_time (signalSize, SignalY, ArrivalTime, FirstPeak, MaxSlope);//same
  if (result == false) {
    ///printf ("  Compute compute_s0_individual_curve fails! S0 = 0.\n");
    return 0;
  }
  
  float* SignalGradient = new float[signalSize];
  //above: same
  compute_gradient(signalSize, SignalY, SignalGradient);    
  
  int count = 0;
  double sum = 0;
  for (int i=0; i<FirstPeak; i++) {
    sum += SignalY[i];
    if (SignalGradient[i] < S0GradThresh) {
      S0 += SignalY[i];
      count++;
    }
  }
  if (count)
    S0 /= count;
  else
    S0 = sum / (FirstPeak-1);
  delete [] SignalGradient;
  return float (S0);

}

}; // end of namespace

