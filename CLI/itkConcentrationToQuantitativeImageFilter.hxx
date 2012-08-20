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
  m_UsePrescribedAIF = false;
  this->Superclass::SetNumberOfRequiredInputs(1);
  this->Superclass::SetNthOutput(1, static_cast<TOutputImage*>(this->MakeOutput(0).GetPointer()));
  this->Superclass::SetNthOutput(2, static_cast<TOutputImage*>(this->MakeOutput(0).GetPointer()));
  this->Superclass::SetNthOutput(3, static_cast<TOutputImage*>(this->MakeOutput(0).GetPointer()));
  this->Superclass::SetNthOutput(4, static_cast<TOutputImage*>(this->MakeOutput(0).GetPointer()));
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


template <class TInputImage, class TMaskImage, class TOutputImage>
void 
ConcentrationToQuantitativeImageFilter<TInputImage,TMaskImage,TOutputImage>
::BeforeThreadedGenerateData()
{
  const VectorVolumeType* inputVectorVolume = this->GetInput();
  const MaskVolumeType* maskVolume = this->GetAIFMask();

  int timeSize = (int)inputVectorVolume->GetNumberOfComponentsPerPixel();

  int   aif_BATIndex = 0;
  int   aif_FirstPeakIndex = 0;
  float aif_MaxSlope = 0.0f;
  
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
  compute_bolus_arrival_time (timeSize, &m_AIF[0], aif_BATIndex, aif_FirstPeakIndex, aif_MaxSlope);  

  // Compute the area under the curve for the AIF
  m_aifAUC = area_under_curve(timeSize, &m_Timing[0], &m_AIF[0], aif_BATIndex, m_AUCTimeInterval);
  
}

template <class TInputImage, class TMaskImage, class TOutputImage>
void 
ConcentrationToQuantitativeImageFilter<TInputImage,TMaskImage,TOutputImage>
#if ITK_VERSION_MAJOR < 4
::ThreadedGenerateData( const typename Superclass::OutputImageRegionType & outputRegionForThread, int threadId )
#else
::ThreadedGenerateData( const typename Superclass::OutputImageRegionType& outputRegionForThread, ThreadIdType threadId )
#endif
  {
  VectorVoxelType vectorVoxel;

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
  OutputVolumeIterType fpvVolumeIter(this->GetFPVOutput(), outputRegionForThread);
  OutputVolumeIterType maxSlopeVolumeIter(this->GetMaxSlopeOutput(), outputRegionForThread);
  OutputVolumeIterType aucVolumeIter(this->GetAUCOutput(), outputRegionForThread);

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
  

  while (!ktransVolumeIter.IsAtEnd() )
    {
    vectorVoxel = inputVectorVolumeIter.Get();    
	
    // Calculate parameter ktrans, ve, and fpv
    pk_solver(timeSize, &timeMinute[0],
              const_cast<float *>(vectorVoxel.GetDataPointer() ),
              &m_AIF[0],
              tempKtrans, tempVe, tempFpv,
              m_fTol,m_gTol,m_xTol,
              m_epsilon,m_maxIter, m_hematocrit,
              optimizer,costFunction);

    // Calculate parameter maxSlope
    compute_bolus_arrival_time (timeSize,
                                const_cast<float *>(vectorVoxel.GetDataPointer() ), BATIndex, FirstPeakIndex, tempMaxSlope);

    // Calculate parameter AUC, normalized by AIF AUC
    tempAUC =
      (area_under_curve(timeSize, &m_Timing[0], const_cast<float *>(vectorVoxel.GetDataPointer() ), BATIndex,  m_AUCTimeInterval) )/m_aifAUC;
  

    ktransVolumeIter.Set(static_cast<OutputVolumePixelType>(tempKtrans) );
    veVolumeIter.Set(static_cast<OutputVolumePixelType>(tempVe) );
    fpvVolumeIter.Set(static_cast<OutputVolumePixelType>(tempFpv));
    maxSlopeVolumeIter.Set(static_cast<OutputVolumePixelType>(tempMaxSlope) );
    aucVolumeIter.Set(static_cast<OutputVolumePixelType>(tempAUC) );

    ++ktransVolumeIter;
    ++veVolumeIter;
    ++fpvVolumeIter;
    ++maxSlopeVolumeIter;
    ++aucVolumeIter;
    ++inputVectorVolumeIter;
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

