#ifndef _itkConcentrationToQuantitativeImageFilter_hxx
#define _itkConcentrationToQuantitativeImageFilter_hxx
#endif

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "itkProgressReporter.h"
#include "itkLevenbergMarquardtOptimizer.h"
#include "vnl/vnl_math.h"

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
  m_UsePrescribedAIF = false;
  m_MaskByRSquared = true;
  m_ModelType = itk::LMCostFunction::TOFTS_2_PARAMETER;
  this->Superclass::SetNumberOfRequiredInputs(1);
  this->Superclass::SetNthOutput(1, static_cast<TOutputImage*>(this->MakeOutput(1).GetPointer()));
  this->Superclass::SetNthOutput(2, static_cast<TOutputImage*>(this->MakeOutput(2).GetPointer()));
  this->Superclass::SetNthOutput(3, static_cast<TOutputImage*>(this->MakeOutput(3).GetPointer()));
  this->Superclass::SetNthOutput(4, static_cast<TOutputImage*>(this->MakeOutput(4).GetPointer()));
  this->Superclass::SetNthOutput(5, static_cast<TOutputImage*>(this->MakeOutput(5).GetPointer()));
  this->Superclass::SetNthOutput(6, static_cast<VectorVolumeType*>(this->MakeOutput(6).GetPointer()));
}

template< class TInputImage, class TMaskImage, class TOutputImage >
typename ConcentrationToQuantitativeImageFilter<TInputImage,TMaskImage,TOutputImage>::DataObjectPointer
ConcentrationToQuantitativeImageFilter<TInputImage,TMaskImage,TOutputImage>
::MakeOutput(DataObjectPointerArraySizeType idx)
{
  if(idx<6)
  {
    return TOutputImage::New().GetPointer();
  }
  else if (idx==6)
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
TInputImage*
ConcentrationToQuantitativeImageFilter< TInputImage,TMaskImage, TOutputImage >
::GetFittedDataOutput()
{
  return dynamic_cast< TInputImage * >( this->ProcessObject::GetOutput(6) );
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
  else if (maskVolume)
    {
    // calculate the AIF from the image using the data under the
    // specified mask
    m_AIF = this->CalculateAverageAIF(inputVectorVolume, maskVolume);
    }
  else
    {
    itkExceptionMacro("A mask image over which to establish the AIF or a prescribed AIF must be assigned. If prescribing an AIF, then UsePrescribedAIF must be set to true.");
    }

  // Compute the bolus arrival time
  compute_bolus_arrival_time (timeSize, &m_AIF[0], m_AIFBATIndex, aif_FirstPeakIndex, aif_MaxSlope);

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

  // std::cout << "AIF = ";
  // for (std::vector<float>::iterator ait = m_AIF.begin(); ait != m_AIF.end(); ++ait)
  //   {
  //   std::cout << *ait << ", ";
  //   }
  // std::cout << std::endl;

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
      // if (ktransVolumeIter.GetIndex()[0] == 122
      //     && ktransVolumeIter.GetIndex()[1] == 118
      //     && ktransVolumeIter.GetIndex()[2] == 6)
      //   {
      //   std::cerr << "VectorVoxel = " << vectorVoxel;
      //   }
     

      // Compute the bolus arrival time and the max slope parameter
      if (success)
        {
        int status = compute_bolus_arrival_time(timeSize, &vectorVoxel[0], BATIndex, FirstPeakIndex, tempMaxSlope);
        if (!status)
          {
          success = false;
          }
        }
     
     
      // Shift the current time course to align with the BAT of the AIF
      // (note the sense of the shift)
      if (success)
        {
        shift = m_AIFBATIndex - BATIndex;
        //std::cerr << "AIF BAT: " << m_AIFBATIndex << ", BAT: " << BATIndex << std::endl;
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
          	optimizer,costFunction,m_ModelType);
        
        itk::LMCostFunction::ParametersType param(3);
        param[0] = tempKtrans; param[1] = tempVe;
        if(m_ModelType == itk::LMCostFunction::TOFTS_3_PARAMETER)
          {
            param[2] = tempFpv;
          }
        itk::LMCostFunction::MeasureType measure =
          costFunction->GetFittedFunction(param);
        for(int i=0;i<fittedVectorVoxel.GetSize();i++)
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
     
    ++ktransVolumeIter;
    ++veVolumeIter;
    ++maxSlopeVolumeIter;
    ++aucVolumeIter;
    ++rsqVolumeIter;
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

