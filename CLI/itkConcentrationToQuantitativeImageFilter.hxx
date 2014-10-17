#ifndef _itkConcentrationToQuantitativeImageFilter_hxx
#define _itkConcentrationToQuantitativeImageFilter_hxx
#endif

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkProgressReporter.h"
#include "itkLevenbergMarquardtOptimizer.h"
#include "vnl/vnl_math.h"

// work around compile error on Windows
#define M_PI 3.1415926535897932384626433832795

#include "itkConcentrationToQuantitativeImageFilter.h"


namespace itk
{

template <class TInputImage, class TMaskImage, class TOutputImage>
ConcentrationToQuantitativeImageFilter<TInputImage,TMaskImage,TOutputImage>::ConcentrationToQuantitativeImageFilter()
{
  m_T1Pre = 0.0f;
  m_TR = 0.0f;
  m_FA = 0.0f;
  m_RGD_relaxivity = 4.9E-3f;
  m_S0GradThresh = 15.0f;
  m_fTol = 1e-4f;
  m_gTol = 1e-4f;
  m_xTol = 1e-5f;
  m_epsilon = 1e-9f;
  m_maxIter = 200;
  m_hematocrit = 0.4f;
  m_aifAUC = 0.0f;
  m_AIFBATIndex = 0;
  m_UsePopulationAIF = false;
  m_UsePrescribedAIF = false;
  m_MaskByRSquared = true;
  m_ModelType = itk::LMCostFunction::TOFTS_2_PARAMETER;
  m_constantBAT = 0;
  m_BATCalculationMode = "PeakGradient";
  this->Superclass::SetNumberOfRequiredInputs(1);
  this->Superclass::SetNthOutput(1, static_cast<TOutputImage*>(this->MakeOutput(1).GetPointer()));  // Ktrans
  this->Superclass::SetNthOutput(2, static_cast<TOutputImage*>(this->MakeOutput(2).GetPointer()));  // Ve
  this->Superclass::SetNthOutput(3, static_cast<TOutputImage*>(this->MakeOutput(3).GetPointer()));  // Max slope
  this->Superclass::SetNthOutput(4, static_cast<TOutputImage*>(this->MakeOutput(4).GetPointer()));  // AUC
  this->Superclass::SetNthOutput(5, static_cast<TOutputImage*>(this->MakeOutput(5).GetPointer()));  // R^2
  this->Superclass::SetNthOutput(6, static_cast<TOutputImage*>(this->MakeOutput(6).GetPointer()));  // BAT
  this->Superclass::SetNthOutput(7, static_cast<VectorVolumeType*>(this->MakeOutput(7).GetPointer())); // fitted
}

template< class TInputImage, class TMaskImage, class TOutputImage >
typename ConcentrationToQuantitativeImageFilter<TInputImage,TMaskImage,TOutputImage>::DataObjectPointer
ConcentrationToQuantitativeImageFilter<TInputImage,TMaskImage,TOutputImage>
::MakeOutput(DataObjectPointerArraySizeType idx)
{
  if(idx<7)
  {
    return TOutputImage::New().GetPointer();
  }
  else if (idx==7)
  {
    return VectorVolumeType::New().GetPointer();
  }
  return 0;
}

// Set a prescribed AIF.  This is not currrently in the input vector,
// though it could be if we used a Decorator.
template< class TInputImage, class TMaskImage, class TOutputImage >
void
ConcentrationToQuantitativeImageFilter< TInputImage, TMaskImage, TOutputImage >
::SetPrescribedAIF(const std::vector<float>& timing, const std::vector<float>& aif)
{
  if (aif.size() < 2)
    {
    itkExceptionMacro(<< "Prescribed AIF must contain at least two time points");
    }
  if (aif.size() != timing.size())
    {
    itkExceptionMacro("Timing vector and concentration vector for AIF must be the same size.");
    }

  m_PrescribedAIF = aif;
  m_PrescribedAIFTiming = timing;
}

// Set 3D AIF mask as second input
template< class TInputImage, class TMaskImage, class TOutputImage >
void
ConcentrationToQuantitativeImageFilter< TInputImage, TMaskImage, TOutputImage >
::SetAIFMask(const TMaskImage* volume)
{
  this->SetNthInput(1, const_cast<TMaskImage*>(volume) );
}

template< class TInputImage, class TMaskImage, class TOutputImage >
const TMaskImage*
ConcentrationToQuantitativeImageFilter< TInputImage, TMaskImage, TOutputImage >
::GetAIFMask() const
{
  return dynamic_cast< const TMaskImage * >( this->ProcessObject::GetInput(1) );
}

// Set 3D ROI mask as third input
template< class TInputImage, class TMaskImage, class TOutputImage >
void
ConcentrationToQuantitativeImageFilter< TInputImage, TMaskImage, TOutputImage >
::SetROIMask(const TMaskImage* volume)
{
  this->SetNthInput(2, const_cast<TMaskImage*>(volume) );
}

template< class TInputImage, class TMaskImage, class TOutputImage >
const TMaskImage*
ConcentrationToQuantitativeImageFilter< TInputImage, TMaskImage, TOutputImage >
::GetROIMask() const
{
  return dynamic_cast< const TMaskImage * >( this->ProcessObject::GetInput(2) );
}

template< class TInputImage, class TMaskImage, class TOutputImage >
TOutputImage*
ConcentrationToQuantitativeImageFilter< TInputImage,TMaskImage, TOutputImage >
::GetKTransOutput()
{
  return dynamic_cast< TOutputImage * >( this->ProcessObject::GetOutput(0) );
}

template< class TInputImage, class TMaskImage, class TOutputImage >
TOutputImage*
ConcentrationToQuantitativeImageFilter< TInputImage,TMaskImage, TOutputImage >
::GetVEOutput()
{
  return dynamic_cast< TOutputImage * >( this->ProcessObject::GetOutput(1) );
}

template< class TInputImage, class TMaskImage, class TOutputImage >
TOutputImage*
ConcentrationToQuantitativeImageFilter< TInputImage,TMaskImage, TOutputImage >
::GetFPVOutput()
{
  return dynamic_cast< TOutputImage * >( this->ProcessObject::GetOutput(2) );
}

template< class TInputImage, class TMaskImage, class TOutputImage >
TOutputImage*
ConcentrationToQuantitativeImageFilter< TInputImage,TMaskImage, TOutputImage >
::GetMaxSlopeOutput()
{
  return dynamic_cast< TOutputImage * >( this->ProcessObject::GetOutput(3) );
}

template< class TInputImage, class TMaskImage, class TOutputImage >
TOutputImage*
ConcentrationToQuantitativeImageFilter< TInputImage,TMaskImage, TOutputImage >
::GetAUCOutput()
{
  return dynamic_cast< TOutputImage * >( this->ProcessObject::GetOutput(4) );
}

template< class TInputImage, class TMaskImage, class TOutputImage >
TOutputImage*
ConcentrationToQuantitativeImageFilter< TInputImage,TMaskImage, TOutputImage >
::GetRSquaredOutput()
{
  return dynamic_cast< TOutputImage * >( this->ProcessObject::GetOutput(5) );
}

template< class TInputImage, class TMaskImage, class TOutputImage >
TOutputImage*
ConcentrationToQuantitativeImageFilter< TInputImage,TMaskImage, TOutputImage >
::GetBATOutput()
{
  return dynamic_cast< TOutputImage * >( this->ProcessObject::GetOutput(6) );
}

template< class TInputImage, class TMaskImage, class TOutputImage >
TInputImage*
ConcentrationToQuantitativeImageFilter< TInputImage,TMaskImage, TOutputImage >
::GetFittedDataOutput()
{
  return dynamic_cast< TInputImage * >( this->ProcessObject::GetOutput(7) );
}

template <class TInputImage, class TMaskImage, class TOutputImage>
void
ConcentrationToQuantitativeImageFilter<TInputImage,TMaskImage,TOutputImage>
::BeforeThreadedGenerateData()
{
  const VectorVolumeType* inputVectorVolume = this->GetInput();
  const MaskVolumeType* maskVolume = this->GetAIFMask();

  std::cout << "Model type: " << m_ModelType << std::endl;

  int timeSize = (int)inputVectorVolume->GetNumberOfComponentsPerPixel();

  int   aif_FirstPeakIndex = 0;
  float aif_MaxSlope = 0.0f;

  // Some of the outputs are optional and may not be calculated.
  // Let's initialize those to all zeros
  OutputVolumeType *fpv = this->GetFPVOutput();
  fpv->FillBuffer(0.0);

  // calculate AIF
  if (m_UsePrescribedAIF)
    {
    // resample the prescribed AIF vector to be at the specifed
    // m_Timing points and then assign to m_AIF
    m_AIF = std::vector<float>(timeSize);

    std::vector<float>::iterator ait = m_AIF.begin();
    std::vector<float>::iterator tit = m_Timing.begin();

    std::vector<float>::iterator pait = m_PrescribedAIF.begin();
    std::vector<float>::iterator ptit = m_PrescribedAIFTiming.begin();

    std::vector<float>::iterator paitnext = pait;
    paitnext++;
    std::vector<float>::iterator ptitnext = ptit;
    ptitnext++;

    for (; tit != m_Timing.end(); ++tit, ++ait)
      {
      // Three cases
      // (1) extrapolate the aif on the low end of the range of prescribed timings
      // (2) interpolate the aif
      // (3) extrapolate the aif on the high end of the range of prescribed timings
      //
      // Case (1) is handled implictly by the initialization and conditionals.
      if (*ptit <= *tit)
        {
        // Case (2) from above)
        // find the prescribed times that straddle the current time to interpolate
        while (*ptitnext < *tit && ptitnext != m_PrescribedAIFTiming.end())
          {
          ++ptit;
          ++ptitnext;
          ++pait;
          ++paitnext;
          }
        }
      if (ptitnext == m_PrescribedAIFTiming.end())
        {
        // we'll need to extrapolate (Case (3) from above)
        ptitnext = ptit;
        --ptit;
        paitnext = pait;
        --pait;
        }

      // interpolate aif;
      float a;
      a = *pait + ((*tit-*ptit) / (*ptitnext - *ptit)) * (*paitnext - *pait);
      *ait = a;
      }
    }
  else if (maskVolume && ! m_UsePopulationAIF)
    {
    // calculate the AIF from the image using the data under the
    // specified mask
    m_AIF = this->CalculateAverageAIF(inputVectorVolume, maskVolume);
    }
  else if (m_UsePopulationAIF)
    {
    m_AIF = this->CalculatePopulationAIF(m_Timing, 0.1);
    }
  else
    {
    itkExceptionMacro("A mask image over which to establish the AIF or a prescribed AIF must be assigned. If prescribing an AIF, then UsePrescribedAIF must be set to true.");
    }
  // Compute the bolus arrival time
  if (m_BATCalculationMode == "UseConstantBAT")
  {
    m_AIFBATIndex = m_constantBAT;
  }
  else if (m_BATCalculationMode == "PeakGradient")
  {
    compute_bolus_arrival_time (m_AIF.size(), &m_AIF[0], m_AIFBATIndex, aif_FirstPeakIndex, aif_MaxSlope);
  }

  // Compute the area under the curve for the AIF
  m_aifAUC = area_under_curve(timeSize, &m_Timing[0], &m_AIF[0], m_AIFBATIndex, m_AUCTimeInterval);
}

template <class TInputImage, class TMaskImage, class TOutputImage>
void
ConcentrationToQuantitativeImageFilter<TInputImage,TMaskImage,TOutputImage>
#if ITK_VERSION_MAJOR < 4
::ThreadedGenerateData( const OutputVolumeRegionType & outputRegionForThread, int threadId )
#else
::ThreadedGenerateData( const OutputVolumeRegionType& outputRegionForThread, ThreadIdType threadId )
#endif
{
  VectorVoxelType vectorVoxel, fittedVectorVoxel;

  float tempFpv = 0.0f;
  float tempKtrans = 0.0f;
  float tempVe = 0.0f;
  float tempMaxSlope = 0.0f;
  float tempAUC = 0.0f;
  int   BATIndex = 0;
  int   FirstPeakIndex = 0;

  const VectorVolumeType* inputVectorVolume = this->GetInput();

  VectorVolumeConstIterType inputVectorVolumeIter(inputVectorVolume, outputRegionForThread);
  OutputVolumeIterType ktransVolumeIter(this->GetKTransOutput(), outputRegionForThread);
  OutputVolumeIterType veVolumeIter(this->GetVEOutput(), outputRegionForThread);
  typename VectorVolumeType::Pointer fitted = this->GetFittedDataOutput();
  VectorVolumeIterType fittedVolumeIter(fitted, outputRegionForThread);

  MaskVolumeConstIterType roiMaskVolumeIter;
  if(this->GetROIMask())
    {
    roiMaskVolumeIter = MaskVolumeConstIterType(this->GetROIMask(), outputRegionForThread);
    }

  OutputVolumeIterType fpvVolumeIter;
  if(m_ModelType == itk::LMCostFunction::TOFTS_3_PARAMETER)
    {
    fpvVolumeIter = OutputVolumeIterType(this->GetFPVOutput(), outputRegionForThread);
    }
  OutputVolumeIterType maxSlopeVolumeIter(this->GetMaxSlopeOutput(), outputRegionForThread);
  OutputVolumeIterType aucVolumeIter(this->GetAUCOutput(), outputRegionForThread);
  OutputVolumeIterType rsqVolumeIter(this->GetRSquaredOutput(), outputRegionForThread);
  OutputVolumeIterType batVolumeIter(this->GetBATOutput(), outputRegionForThread);

  //set up optimizer and cost function
  itk::LevenbergMarquardtOptimizer::Pointer optimizer = itk::LevenbergMarquardtOptimizer::New();
  LMCostFunction::Pointer                   costFunction = LMCostFunction::New();
  int timeSize = (int)inputVectorVolume->GetNumberOfComponentsPerPixel();

  std::vector<float> timeMinute;
  timeMinute = m_Timing;
  for(unsigned int i = 0; i < timeMinute.size(); i++)
    {
    timeMinute[i] = m_Timing[i]/60.0;
    }

  ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());

  // Cache the RMS error of fitting the model to the AIF
  // pk_solver(timeSize, &timeMinute[0],
  //           &m_AIF[0],
  //           &m_AIF[0],
  //           tempKtrans, tempVe, tempFpv,
  //           m_fTol,m_gTol,m_xTol,
  //           m_epsilon,m_maxIter, m_hematocrit,
  //           optimizer,costFunction);

  // double aifRMS = optimizer->GetOptimizer()->get_end_error();
  // std::cout << "AIF RMS: " << aifRMS  << std::endl;


  VectorVoxelType shiftedVectorVoxel(timeSize);
  int shift;
  unsigned int shiftStart = 0, shiftEnd = 0;
  bool success = true;
  while (!ktransVolumeIter.IsAtEnd())
    {
    success = true;
    tempKtrans = tempVe = tempFpv = tempMaxSlope = tempAUC = 0.0;
    BATIndex = FirstPeakIndex = 0;

    if(!this->GetROIMask() || (this->GetROIMask() && roiMaskVolumeIter.Get()))
      {
      vectorVoxel = inputVectorVolumeIter.Get();
      fittedVectorVoxel = inputVectorVolumeIter.Get();
      // dump a specific voxel
      // std::cout << "VectorVoxel = " << vectorVoxel;

      // Compute the bolus arrival time and the max slope parameter
      if (success)
        {
          int status;
          // Compute the bolus arrival time
          if (m_BATCalculationMode == "UseConstantBAT")
          {
            BATIndex = m_constantBAT;
             status = 1;
          }
          else if (m_BATCalculationMode == "PeakGradient")
          {
             status = compute_bolus_arrival_time(timeSize, &vectorVoxel[0], BATIndex, FirstPeakIndex, tempMaxSlope);
          }

        if (!status)
          {
          success = false;
          }
        }
     
     
      // Shift the current time course to align with the BAT of the AIF
      // (note the sense of the shift)
      if (success)
        {
        batVolumeIter.Set(BATIndex);
        shift = m_AIFBATIndex - BATIndex;
        shiftedVectorVoxel.Fill(0.0);
        if (shift <= 0)
          {
          // AIF BAT before current BAT, should always be the case
          shiftStart = 0;
          shiftEnd = vectorVoxel.Size() + shift;
          }
        else
          {
          success = false;
          }
        }
      if (success)
        {
        for (unsigned int i = shiftStart; i < shiftEnd; ++i)
          {
          shiftedVectorVoxel[i] = vectorVoxel[i - shift];
          }
        }

      // Calculate parameter ktrans, ve, and fpv
      double rSquared = 0.0;
      if (success)
        {
        pk_solver(timeSize, &timeMinute[0],
          	const_cast<float *>(shiftedVectorVoxel.GetDataPointer() ),
          	&m_AIF[0],
          	tempKtrans, tempVe, tempFpv,
          	m_fTol,m_gTol,m_xTol,
          	m_epsilon,m_maxIter, m_hematocrit,
            optimizer,costFunction,m_ModelType,m_constantBAT,m_BATCalculationMode);

        itk::LMCostFunction::ParametersType param(3);
        param[0] = tempKtrans; param[1] = tempVe;
        if(m_ModelType == itk::LMCostFunction::TOFTS_3_PARAMETER)
          {
          param[2] = tempFpv;
          }
        itk::LMCostFunction::MeasureType measure =
          costFunction->GetFittedFunction(param);
        for(size_t i=0;i<fittedVectorVoxel.GetSize();i++)
          {
            fittedVectorVoxel[i] = measure[i];
          }
        
        // Shift the current time course to align with the BAT of the AIF
        // (note the sense of the shift)
        shiftedVectorVoxel.Fill(0.0);
        if (shift <= 0)
          {
          // AIF BAT before current BAT, should always be the case
          shiftStart = shift*-1.;
          shiftEnd = vectorVoxel.Size();
          for (unsigned int i = shiftStart; i < shiftEnd; ++i)
            {
            shiftedVectorVoxel[i] = fittedVectorVoxel[i + shift];
            }
          }
  
        fittedVolumeIter.Set(shiftedVectorVoxel);

        // Only keep the estimated values if the optimization produced a good answer
        // Check R-squared:
        //   R2 = 1 - SSerr / SStot
        // where
        //   SSerr = \sum (y_i - f_i)^2
        //   SStot = \sum (y_i - \bar{y})^2
        //
        // Note: R-squared is not a good metric for nonlinear function
        // fitting. R-squared values are not bound between [0,1] when
        // fitting nonlinear functions.
     
        // SSerr we can get easily from the optimizer
        double rms = optimizer->GetOptimizer()->get_end_error();
        double SSerr = rms*rms*shiftedVectorVoxel.GetSize();
     
        // if we couldn't get rms from the optimizer, we would calculate SSerr ourselves
        // LMCostFunction::MeasureType residuals = costFunction->GetValue(optimizer->GetCurrentPosition());
        // double SSerr = 0.0;
        // for (unsigned int i=0; i < residuals.size(); ++i)
        //   {
        //   SSerr += (residuals[i]*residuals[i]);
        //   }
     
        // SStot we need to calculate
        double sumSquared = 0.0;
        double sum = 0.0;
        for (unsigned int i=0; i < shiftedVectorVoxel.GetSize(); ++i)
          {
          sum += shiftedVectorVoxel[i];
          sumSquared += (shiftedVectorVoxel[i]*shiftedVectorVoxel[i]);
          }
        double SStot = sumSquared - sum*sum/(double)shiftedVectorVoxel.GetSize();
     
        rSquared = 1.0 - (SSerr / SStot);
     
        double rSquaredThreshold = 0.15;
        if (rSquared < rSquaredThreshold)
          {
          success = false;
          }
        }
      // Calculate parameter AUC, normalized by AIF AUC
      if (success)
        {
        tempAUC =
          (area_under_curve(timeSize, &m_Timing[0], const_cast<float *>(shiftedVectorVoxel.GetDataPointer() ), BATIndex,  m_AUCTimeInterval) )/m_aifAUC;
        }
     
      // Do we mask the output volumes by the R-squared value?
      if (m_MaskByRSquared)
        {
        // If we were successful, save the estimated values, otherwise
        // default to zero
        if (success)
          {
          ktransVolumeIter.Set(static_cast<OutputVolumePixelType>(tempKtrans) );
          veVolumeIter.Set(static_cast<OutputVolumePixelType>(tempVe) );
          maxSlopeVolumeIter.Set(static_cast<OutputVolumePixelType>(tempMaxSlope) );
          aucVolumeIter.Set(static_cast<OutputVolumePixelType>(tempAUC) );
          if(m_ModelType == itk::LMCostFunction::TOFTS_3_PARAMETER)
            {
            fpvVolumeIter.Set(static_cast<OutputVolumePixelType>(tempFpv));
            }
          }
        else
          {
          ktransVolumeIter.Set(static_cast<OutputVolumePixelType>(0) );
          veVolumeIter.Set(static_cast<OutputVolumePixelType>(0) );
          maxSlopeVolumeIter.Set(static_cast<OutputVolumePixelType>(0) );
          aucVolumeIter.Set(static_cast<OutputVolumePixelType>(0) );
          }
        }
      else
        {
        ktransVolumeIter.Set(static_cast<OutputVolumePixelType>(tempKtrans) );
        veVolumeIter.Set(static_cast<OutputVolumePixelType>(tempVe) );
        maxSlopeVolumeIter.Set(static_cast<OutputVolumePixelType>(tempMaxSlope) );
        aucVolumeIter.Set(static_cast<OutputVolumePixelType>(tempAUC) );
        if(m_ModelType == itk::LMCostFunction::TOFTS_3_PARAMETER)
          {
          fpvVolumeIter.Set(static_cast<OutputVolumePixelType>(tempFpv));
          }
        }
     
      // RSquared output volume is always written
      rsqVolumeIter.Set(rSquared);
      }
    else
      {
      ktransVolumeIter.Set(static_cast<OutputVolumePixelType>(0) );
      veVolumeIter.Set(static_cast<OutputVolumePixelType>(0) );
      maxSlopeVolumeIter.Set(static_cast<OutputVolumePixelType>(0) );
      aucVolumeIter.Set(static_cast<OutputVolumePixelType>(0) );
      }

    ++ktransVolumeIter;
    ++veVolumeIter;
    ++maxSlopeVolumeIter;
    ++aucVolumeIter;
    ++rsqVolumeIter;
    ++batVolumeIter;
    ++inputVectorVolumeIter;
    ++fittedVolumeIter;

    if(this->GetROIMask())
      {
      ++roiMaskVolumeIter;
      }

    if(m_ModelType == itk::LMCostFunction::TOFTS_3_PARAMETER)
      {
      ++fpvVolumeIter;
      }

    progress.CompletedPixel();
  }
}

// Calculate a population AIF.
//
// See "Experimentally-Derived Functional Form for a Population-Averaged High-
// Temporal-Resolution Arterial Input Function for Dynamic Contrast-Enhanced
// MRI" - Parker, Robers, Macdonald, Buonaccorsi, Cheung, Buckley, Jackson,
// Watson, Davies, Jayson.  Magnetic Resonance in Medicine 56:993-1000 (2006)
template <class TInputImage, class TMaskImage, class TOutputImage>
std::vector<float>
ConcentrationToQuantitativeImageFilter<TInputImage, TMaskImage, TOutputImage>
::CalculatePopulationAIF(std::vector<float> signalTime, const float bolusArrivalTimeFraction)
{

    // Inputs
    // ------
    // signalTime : sequence time, presumed in units of seconds.
    // bolusArrivalTimeFraction : fractional point between 0 and 1 when the bolus is
    //     desired to arrive.  Choose 0.0 to have it at the very beginning,
    //     1.0 to have it at the end.
    //
    // Outputs
    // -------
    // AIF : arterial input function as a function of time
    std::vector<float> AIF;

    // Make a high resolution timing vector as input to the AIF construction.
    std::vector<float> aif_time(signalTime.size() * 10);
    float final_time_point = signalTime[signalTime.size()-1];
    float resolution = final_time_point / (aif_time.size() - 1);
    for (size_t j = 0; j < aif_time.size(); ++j) {
	    aif_time[j] = resolution * j; 
    }

    size_t bolus_arrival_time_idx = (float)aif_time.size() * bolusArrivalTimeFraction;

    size_t n = aif_time.size();
    AIF.resize(n);

    size_t numTimePoints = n - bolus_arrival_time_idx;
    std::vector<float> timeSinceBolus(numTimePoints);


    // t=FR*[0:numTimePoints-1]/60;
    // These time points "start" when the bolus arrives.
    for ( size_t j = 0; j < numTimePoints; ++j ) {
        //timeSinceBolus[j] = FR * j / 60.0;
        timeSinceBolus[j] = aif_time[bolus_arrival_time_idx + j] - aif_time[bolus_arrival_time_idx];
    }

    // Parker
    // defining parameters
    double a1(0.809);
    double a2(0.330);
    double T1(0.17406);
    double T2(0.365);
    double sigma1(0.0563);
    double sigma2(0.132);
    double alpha(1.050);
    double beta(0.1685);
    double s(38.078);
    double tau(0.483);


    // term0=alpha*exp(-beta*t)./(1+exp(-s*(t-tau)));
    // Here the assumption is that time is in minutes, so must scale accordingly.
    // see Parker.
    std::vector<double> term0(numTimePoints);
    for ( size_t j = 0; j < numTimePoints; ++j ) {
        term0[j] = alpha * exp(-beta*timeSinceBolus[j]/60.0) 
		 / (1 + exp( -s * (timeSinceBolus[j]/60.0 - tau)));
    }


    // term1=[];
    // term2=[];
    double A1 = a1 / (sigma1 * pow((2*M_PI), 0.5));

    // B1=exp(-(t-T1).^2./(2.*sigma1^2));
    double numerator, denominator;
    std::vector<double> B1(numTimePoints);
    denominator = 2.0 * pow(sigma1, 2.0);
    for ( size_t j = 0; j < numTimePoints; ++j ) {
        numerator = -1 * pow(-(timeSinceBolus[j]/60.0 - T1), 2.0);
        B1[j] = exp( numerator / denominator );
    }

    // term1=A1.*B1;
    std::vector<double> term1(numTimePoints);
    for ( size_t j = 0; j < numTimePoints; ++j ) {
        term1[j] = A1 * B1[j];
    }


    // A2=a2/(sigma2*((2*pi)^0.5));
    double A2 = a2 / (sigma2 * pow(2*M_PI, 0.5));

    //B2=exp(-(t-T2).^2./(2.*sigma2^2));
    std::vector<double> B2(numTimePoints);
    denominator = 2.0 * pow(sigma2, 2.0);
    for ( size_t j = 0; j < numTimePoints; ++j ) {
        numerator = -1 * pow(-(timeSinceBolus[j]/60.0 - T2), 2.0);
        B2[j] = exp(numerator / denominator);
    }

    // term2=A2.*B2;
    std::vector<double> term2(numTimePoints);
    for ( size_t j = 0; j < numTimePoints; ++j ) {
        term2[j] = A2 * B2[j];
    }

    // aifPost=term0+term1+term2;
    std::vector<double> aifPost(numTimePoints);
    for ( size_t j = 0; j < numTimePoints; ++j ) {
        aifPost[j] = term0[j] + term1[j] + term2[j];
    }

    // Initialize values before bolus arrival.
    for ( size_t j = 0; j < bolus_arrival_time_idx; ++j ) {
        AIF[j] = 0;
    }

    // Shift the data to take into account the bolus arrival time.
    // sp=timeOfBolus+1;
    // AIF(sp:end)=aifPost;
    for ( size_t j = bolus_arrival_time_idx; j < AIF.size(); ++j ) {
        AIF[j] = aifPost[j - bolus_arrival_time_idx];
    }

    // Resample back to signal (sequence) time.
    std::vector<float> rAIF = this->ResampleAIF(aif_time, AIF, signalTime);

    return rAIF;

}


template <class TInputImage, class TMaskImage, class TOutputImage>
std::vector<float>
ConcentrationToQuantitativeImageFilter<TInputImage, TMaskImage, TOutputImage>
::ResampleAIF(std::vector<float> t1, std::vector<float> y1, std::vector<float> t2)
{
    // Resample time1, y1 to time2
    size_t timeSize = t2.size();
    std::vector<float> y2(timeSize);

    std::vector<float>::iterator y2it = y2.begin();
    std::vector<float>::iterator t2it = t2.begin();

    std::vector<float>::iterator y1it = y1.begin();
    std::vector<float>::iterator t1it = t1.begin();

    std::vector<float>::iterator y1itnext = y1it;
    y1itnext++;
    std::vector<float>::iterator t1itnext = t1it;
    t1itnext++;

    for (; t2it != t2.end(); ++t2it, ++y2it)
      {
      // Three cases
      // (1) extrapolate the aif on the low end of the range of prescribed timings
      // (2) interpolate the aif
      // (3) extrapolate the aif on the high end of the range of prescribed timings
      //
      // Case (1) is handled implictly by the initialization and conditionals.
      if (*t1it <= *t2it)
        {
        // Case (2) from above)
        // find the prescribed times that straddle the current time to interpolate
        while (*t1itnext < *t2it && t1itnext != t1.end())
          {
          ++t1it;
          ++t1itnext;
          ++y1it;
          ++y1itnext;
          }
        }
      if (t1itnext == t1.end())
        {
        // we'll need to extrapolate (Case (3) from above)
        t1itnext = t1it;
        --t1it;
        y1itnext = y1it;
        --y1it;
        }

      // interpolate aif;
      float a;
      a = *y1it + ((*t2it-*t1it) / (*t1itnext - *t1it)) * (*y1itnext - *y1it);
      *y2it = a;
      }

    return y2;
}


// Calculate average AIF according to the AIF mask
template <class TInputImage, class TMaskImage, class TOutputImage>
std::vector<float>
ConcentrationToQuantitativeImageFilter<TInputImage, TMaskImage, TOutputImage>
::CalculateAverageAIF(const VectorVolumeType*  inputVectorVolume, const MaskVolumeType* maskVolume)
{
  std::vector<float> averageAIF;

  VectorVolumeConstIterType inputVectorVolumeIter(inputVectorVolume, inputVectorVolume->GetRequestedRegion() );
  MaskVolumeConstIterType  maskVolumeIter(maskVolume, maskVolume->GetRequestedRegion() );

  inputVectorVolumeIter.GoToBegin();
  maskVolumeIter.GoToBegin();

  VectorVoxelType vectorVoxel;
  long            numberVoxels = 0;
  long            numberOfSamples = inputVectorVolume->GetNumberOfComponentsPerPixel();
  averageAIF = std::vector<float>(numberOfSamples, 0.0);

  while (!inputVectorVolumeIter.IsAtEnd() )
    {
    if (maskVolumeIter.Get()!=0) // Mask pixel with value !0 will is part of AIF
      {
      numberVoxels++;
      vectorVoxel = inputVectorVolumeIter.Get();

      for(long i = 0; i < numberOfSamples; i++)
        {
        averageAIF[i] += vectorVoxel[i];
        }
      }
    ++maskVolumeIter;
    ++inputVectorVolumeIter;
    }

  for(long i = 0; i < numberOfSamples; i++)
    {
    averageAIF[i] /= (double)numberVoxels;
    }

  return averageAIF;
}


template <class TInputImage, class TMaskImage, class TOutputImage>
void ConcentrationToQuantitativeImageFilter<TInputImage,TMaskImage,TOutputImage>
::SetTiming(const std::vector<float>& inputTiming)
{
  m_Timing = inputTiming;
}

template <class TInputImage, class TMaskImage, class TOutputImage>
const std::vector<float>& ConcentrationToQuantitativeImageFilter<TInputImage,TMaskImage,TOutputImage>
::GetTiming()
{
  return m_Timing;
}

template <class TInputImage, class TMaskImage, class TOutputImage>
void ConcentrationToQuantitativeImageFilter<TInputImage,TMaskImage,TOutputImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "Function tolerance: " << m_fTol << std::endl;
  os << indent << "Gradient tolerance: " << m_gTol << std::endl;
  os << indent << "Parameter tolerance: " << m_xTol << std::endl;
  os << indent << "Epsilon: " << m_epsilon << std::endl;
  os << indent << "Maximum number of iterations: " << m_maxIter << std::endl;
  os << indent << "Hematocrit: " << m_hematocrit << std::endl;
}

} // end namespace itk

