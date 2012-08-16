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
  this->Superclass::SetNumberOfRequiredInputs(2);
  this->Superclass::SetNumberOfRequiredOutputs(4);
  this->Superclass::SetNthOutput(1, static_cast<TOutputImage*>(this->MakeOutput(0).GetPointer()));
  this->Superclass::SetNthOutput(2, static_cast<TOutputImage*>(this->MakeOutput(0).GetPointer()));
  this->Superclass::SetNthOutput(3, static_cast<TOutputImage*>(this->MakeOutput(0).GetPointer()));
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
::GetMaxSlopeOutput() 
{
  return dynamic_cast< TOutputImage * >( this->ProcessObject::GetOutput(2) );
}

template< class TInputImage, class TMaskImage, class TOutputImage >
TOutputImage* 
ConcentrationToQuantitativeImageFilter< TInputImage,TMaskImage, TOutputImage >
::GetAUCOutput() 
{
  return dynamic_cast< TOutputImage * >( this->ProcessObject::GetOutput(3) );
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
  
  // calculate AUC of AIF area  
  m_averageAIFConcentration = this->CalculateAverageAIF(inputVectorVolume, maskVolume); 
  compute_bolus_arrival_time (timeSize, &m_averageAIFConcentration[0], aif_BATIndex, aif_FirstPeakIndex, aif_MaxSlope);  
  m_aifAUC = area_under_curve(timeSize, &m_TimeAxis[0], &m_averageAIFConcentration[0], aif_BATIndex, m_AUCTimeInterval);
  
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
  OutputVolumeIterType maxSlopeVolumeIter(this->GetMaxSlopeOutput(), outputRegionForThread);
  OutputVolumeIterType aucVolumeIter(this->GetAUCOutput(), outputRegionForThread);

  //set up optimizer and cost function
  itk::LevenbergMarquardtOptimizer::Pointer optimizer = itk::LevenbergMarquardtOptimizer::New(); 
  LMCostFunction::Pointer                   costFunction = LMCostFunction::New();                
  int timeSize = (int)inputVectorVolume->GetNumberOfComponentsPerPixel();

  std::vector<float> timeMinute;
  timeMinute = m_TimeAxis;
  for(unsigned int i = 0; i < timeMinute.size(); i++)
    {
    timeMinute[i] = m_TimeAxis[i]/60.0;    
    }


  ProgressReporter progress(this, threadId, outputRegionForThread.GetNumberOfPixels());
  

  while (!ktransVolumeIter.IsAtEnd() )
    {
    vectorVoxel = inputVectorVolumeIter.Get();    
	
    // Calculate parameter ktrans and ve
    pk_solver(timeSize, &timeMinute[0],
              const_cast<float *>(vectorVoxel.GetDataPointer() ),
              &m_averageAIFConcentration[0],
              tempKtrans, tempVe, tempFpv,
              m_fTol,m_gTol,m_xTol,
              m_epsilon,m_maxIter, m_hematocrit,
              optimizer,costFunction);

    // Calculate parameter maxSlope
    compute_bolus_arrival_time (timeSize,
                                const_cast<float *>(vectorVoxel.GetDataPointer() ), BATIndex, FirstPeakIndex, tempMaxSlope);

    // Calculate parameter AUC, normalized by AIF AUC
    tempAUC =
      (area_under_curve(timeSize, &m_TimeAxis[0], const_cast<float *>(vectorVoxel.GetDataPointer() ), BATIndex,  m_AUCTimeInterval) )/m_aifAUC;
  

    ktransVolumeIter.Set(static_cast<OutputVolumePixelType>(tempKtrans) );
    veVolumeIter.Set(static_cast<OutputVolumePixelType>(tempVe) );
    maxSlopeVolumeIter.Set(static_cast<OutputVolumePixelType>(tempMaxSlope) );
    aucVolumeIter.Set(static_cast<OutputVolumePixelType>(tempAUC) );

    ++ktransVolumeIter;
    ++veVolumeIter;
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
::SetTimeAxis(const std::vector<float>& inputTimeAxis)
{  
  m_TimeAxis = inputTimeAxis;
}

template <class TInputImage, class TMaskImage, class TOutputImage>
const std::vector<float>& ConcentrationToQuantitativeImageFilter<TInputImage,TMaskImage,TOutputImage>
::GetTimeAxis()
{
  return m_TimeAxis;
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

