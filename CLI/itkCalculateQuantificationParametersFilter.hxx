#ifndef _itkCalculateQuantificationParametersFilter_hxx
#define _itkCalculateQuantificationParametersFilter_hxx
#endif

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "vnl/vnl_math.h"
#include "itkCastImageFilter.h"
#include "itkCalculateQuantificationParametersFilter.h"

#include "itkBinaryThresholdImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkLevenbergMarquardtOptimizer.h"
#include <ostream>
#include "stdlib.h"
#include "stdio.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
CalculateQuantificationParametersFilter<TInputImage,TOutputImage>::CalculateQuantificationParametersFilter()
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
  m_timeSize = 0;
  m_timeAxis = new float();
  m_averageAIFCon = new float();
  this->Superclass::SetNumberOfRequiredInputs(2);
  this->Superclass::SetNumberOfRequiredOutputs(4);
  this->Superclass::SetNthOutput(0, VolumeType::New() );
  this->Superclass::SetNthOutput(1, VolumeType::New() );
  this->Superclass::SetNthOutput(2, VolumeType::New() );
  this->Superclass::SetNthOutput(3, VolumeType::New() );
}

template<class TInputImage, class TOutputImage>
void CalculateQuantificationParametersFilter<TInputImage,TOutputImage>
::CallCopyOutputRegionToInputRegion(MultiVolumeRegionType &destRegion, const VolumeRegionType &srcRegion)
{
  ExtractImageFilterRegionCopierType extractImageRegionCopier;

  extractImageRegionCopier(destRegion, srcRegion, m_ExtractionRegion);
}
// Change the input region size, otherwise it would be same as the output size
template< class TInputImage, class TOutputImage >
void CalculateQuantificationParametersFilter< TInputImage, TOutputImage >
::GenerateInputRequestedRegion()
throw( InvalidRequestedRegionError )
{
  Superclass::GenerateInputRequestedRegion();
  
  typename CalculateQuantificationParametersFilter<TInputImage,TOutputImage>::MultiVolumePointerType image =
    const_cast< MultiVolumeType * >( this->GetInput(0) );
  if ( image )
    {
    image->SetRequestedRegion( this->GetInput(0)->GetLargestPossibleRegion() );
    }
}

// Set 4D concentration values as first input
template< class TInputImage, class TOutputImage >
void CalculateQuantificationParametersFilter< TInputImage, TOutputImage >::SetInputMultiVolume(const TInputImage* multiVolume)
{
  SetNthInput(0, const_cast<TInputImage*>(multiVolume) );
}

// Set 3D AIF mask as second input
template< class TInputImage, class TOutputImage >
void CalculateQuantificationParametersFilter< TInputImage, TOutputImage >::SetInputVolume(const TOutputImage* volume)
{
  SetNthInput(1, const_cast<TOutputImage*>(volume) );
}

template< class TInputImage, class TOutputImage >
typename TInputImage::ConstPointer CalculateQuantificationParametersFilter< TInputImage, TOutputImage >::GetInputMultiVolume()
{
  return static_cast< const TInputImage * >
         ( this->ProcessObject::GetInput(0) );
}

template< class TInputImage, class TOutputImage >
typename TOutputImage::ConstPointer CalculateQuantificationParametersFilter< TInputImage,
                                                                             TOutputImage >::GetInputVolume()
{
  return static_cast< const TOutputImage * >
         ( this->ProcessObject::GetInput(1) );
}

// Set the output size same as the x, y, z of input 4D data
template <class TInputImage, class TOutputImage>
void CalculateQuantificationParametersFilter<TInputImage,TOutputImage>
::GenerateOutputInformation()
{
  // do not call the superclass' implementation of this method since
  // this filter allows the input and the output to be of different dimensions

  // get pointers to the input and output
  std::cout <<std::endl<< "Generate output information" << std::endl;
  typename Superclass::OutputImagePointer      outputPtr = this->GetOutput();
  typename Superclass::InputImageConstPointer  inputPtr  = this->GetInput();

  if ( !outputPtr || !inputPtr)
    {
    return;
    }

  typename TOutputImage::RegionType outputImageRegion;
  typename TOutputImage::IndexType outputStartIndex;
  outputStartIndex.Fill(0);
  typename TOutputImage::SizeType outputSize;
  typename TInputImage::SizeType inputSize = inputPtr->GetLargestPossibleRegion().GetSize();
  outputSize[0]=inputSize[0];
  outputSize[1]=inputSize[1];
  outputSize[2]=inputSize[2];
  outputImageRegion.SetSize(outputSize);
  outputImageRegion.SetIndex(outputStartIndex);
  
  outputPtr->SetLargestPossibleRegion( outputImageRegion );
}

template <class TInputImage, class TOutputImage>
void CalculateQuantificationParametersFilter<TInputImage,TOutputImage>::BeforeThreadedGenerateData()
{

  m_inputMultiVolume = this->GetInputMultiVolume();
  m_inputVolume = this->GetInputVolume();

  //transfer input multiVolume to vector image
  m_inputVectorVolume =
    this->MultiVolumeToVectorVolume(const_cast<MultiVolumeType *>(static_cast<const MultiVolumeType * >(m_inputMultiVolume) ) );

  VectorVolumeSizeType vectorVolumeSize = m_inputVectorVolume->GetLargestPossibleRegion().GetSize();

  //calculate parameters
  m_ktransVolume = VolumeType::New();
  m_ktransVolume->SetRegions(m_inputVolume->GetLargestPossibleRegion() );
  m_ktransVolume->Allocate();
  m_ktransVolume->FillBuffer(0);
  m_veVolume = VolumeType::New();
  m_veVolume->SetRegions(m_inputVolume->GetLargestPossibleRegion() );
  m_veVolume->Allocate();
  m_veVolume->FillBuffer(0);
  m_maxSlopeVolume = VolumeType::New();
  m_maxSlopeVolume->SetRegions(m_inputVolume->GetLargestPossibleRegion() );
  m_maxSlopeVolume->Allocate();
  m_maxSlopeVolume->FillBuffer(0);
  m_aucVolume = VolumeType::New();
  m_aucVolume->SetRegions(m_inputVolume->GetLargestPossibleRegion() );
  m_aucVolume->Allocate();
  m_aucVolume->FillBuffer(0);

  m_timeSize = (int)m_inputVectorVolume->GetNumberOfComponentsPerPixel();

  m_TimeMinute = new float[m_inputMultiVolume->GetLargestPossibleRegion().GetSize()[3]]();
  
  //convert second to minute for time series
  for(unsigned int i =0; i<(m_inputMultiVolume->GetLargestPossibleRegion().GetSize()[3]); i++)
    {
    m_TimeMinute[i] = m_timeAxis[i]/60;    
    }

  m_averageAIFCon = new float[m_timeSize]();
  int   aif_BATIndex = 0;
  int   aif_FirstPeakIndex = 0;
  float aif_MaxSlope = 0.0f;
  
  //Calculate AUC of AIF area  
  CalculateAverageAIF(m_inputVectorVolume, m_inputVolume, m_averageAIFCon);  
  compute_bolus_arrival_time (m_timeSize, m_averageAIFCon, aif_BATIndex, aif_FirstPeakIndex, aif_MaxSlope);  
  m_aifAUC = area_under_curve(m_timeSize, m_timeAxis, m_averageAIFCon, aif_BATIndex, m_AUCTimeInterval);
  
  std::cerr<<"beforeThreaded done!"<<std::endl;
}

template <class TInputImage, class TOutputImage>
void CalculateQuantificationParametersFilter<TInputImage,TOutputImage>
#if ITK_VERSION_MAJOR < 4
::ThreadedGenerateData( const typename Superclass::OutputImageRegionType & outputRegionForThread, int itkNotUsed(
                          threadId) )
#else
::ThreadedGenerateData( const typename Superclass::OutputImageRegionType& outputRegionForThread,
                        ThreadIdType itkNotUsed(
                          threadId) )
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

  VectorVolumeIterType inputVectorVolumeIter(m_inputVectorVolume, outputRegionForThread);
  VolumeIterType       ktransVolumeIter(m_ktransVolume, outputRegionForThread);
  VolumeIterType       veVolumeIter(m_veVolume, outputRegionForThread);
  VolumeIterType       maxSlopeVolumeIter(m_maxSlopeVolume, outputRegionForThread);
  VolumeIterType       aucVolumeIter(m_aucVolume, outputRegionForThread);

  //set up optimizer and cost function
  itk::LevenbergMarquardtOptimizer::Pointer optimizer = itk::LevenbergMarquardtOptimizer::New(); 
  LMCostFunction::Pointer                   costFunction = LMCostFunction::New();                
  
  costFunction->SetNumberOfValues (m_timeSize);   
  costFunction->set_hematocrit (m_hematocrit);   
    
  optimizer->UseCostFunctionGradientOff();   
  optimizer->SetUseCostFunctionGradient(0);  

  CommandIterationUpdateLevenbergMarquardt::Pointer observer = CommandIterationUpdateLevenbergMarquardt::New(); 
  optimizer->AddObserver( itk::IterationEvent(), observer );                                                    
  optimizer->AddObserver( itk::FunctionEvaluationIterationEvent(), observer );                                  
  
  while (!inputVectorVolumeIter.IsAtEnd() )
    {

    vectorVoxel = inputVectorVolumeIter.Get();    
	
	// Calculate parameter ktrans and ve
    pk_solver_boost(m_timeSize, m_TimeMinute,
                    const_cast<float *>(vectorVoxel.GetDataPointer() ),m_averageAIFCon,
                    tempKtrans, tempVe, tempFpv,
                    m_fTol,m_gTol,m_xTol,
                    m_epsilon,m_maxIter,
                    optimizer,costFunction);
    // Calculate parameter maxSlope
	compute_bolus_arrival_time (m_timeSize,
                                const_cast<float *>(vectorVoxel.GetDataPointer() ), BATIndex, FirstPeakIndex,
                                tempMaxSlope);
	// Calculate parameter AUC, normalized by AIF AUC
    tempAUC =
      (area_under_curve(m_timeSize, m_timeAxis, const_cast<float *>(vectorVoxel.GetDataPointer() ), BATIndex,
                        m_AUCTimeInterval) )/m_aifAUC;
  

    ktransVolumeIter.Set(static_cast<VolumePixelType>(tempKtrans) );
    veVolumeIter.Set(static_cast<VolumePixelType>(tempVe) );
    maxSlopeVolumeIter.Set(static_cast<VolumePixelType>(tempMaxSlope) );
    aucVolumeIter.Set(static_cast<VolumePixelType>(tempAUC) );

    ++ktransVolumeIter;
    ++veVolumeIter;
    ++maxSlopeVolumeIter;
    ++aucVolumeIter;
    ++inputVectorVolumeIter;
    }
  }

template <class TInputImage, class TOutputImage>
void CalculateQuantificationParametersFilter<TInputImage,TOutputImage>::AfterThreadedGenerateData()
{
  std::cerr << "Prepare for output" << std::endl;
  this->SetNthOutput(0,m_ktransVolume);
  this->SetNthOutput(1,m_veVolume);
  this->SetNthOutput(2,m_maxSlopeVolume);
  this->SetNthOutput(3,m_aucVolume);

  delete [] m_timeAxis;
  delete [] m_TimeMinute;
  delete [] m_averageAIFCon;
}

// Calculate average AIF according to the AIF mask
template <class TInputImage, class TOutputImage>
void
CalculateQuantificationParametersFilter<TInputImage, TOutputImage>::CalculateAverageAIF(
  VectorVolumePointerType inputVectorVolume, VolumeConstPointerType inputVolume,float*& averageAIF)
{
  VectorVolumeIterType inputVectorVolumeIter(inputVectorVolume, inputVectorVolume->GetRequestedRegion() );
  VolumeConstIterType  inputVolumeIter(inputVolume, inputVolume->GetRequestedRegion() );

  inputVectorVolumeIter.GoToBegin();
  inputVolumeIter.GoToBegin();
  VectorVoxelType vectorVoxel;
  VolumeIndexType volumeIndex;
  int             numberVoxels = 0;
 
  while (!inputVectorVolumeIter.IsAtEnd() )
    {
    if (inputVolumeIter.Get()!=0) //Pixel value >0 will be consider to be a
                                  // landmark pixel.
      {
      volumeIndex = inputVolumeIter.GetIndex();      
    
      numberVoxels +=1;
      vectorVoxel = inputVectorVolumeIter.Get();
      
      for(int i = 0; i < (int)inputVectorVolume->GetNumberOfComponentsPerPixel(); i++)
        {
        averageAIF[i]+=(float)vectorVoxel.GetDataPointer()[i];    
        }      
      }
    ++inputVolumeIter;
    ++inputVectorVolumeIter;
    }
  
  for(int i = 0; i < (int)inputVectorVolume->GetNumberOfComponentsPerPixel(); i++)
    {
    averageAIF[i]=averageAIF[i]/numberVoxels;
    }
}

template <class TInputImage, class TOutputImage>
typename CalculateQuantificationParametersFilter<TInputImage, TOutputImage>::VectorVolumePointerType
CalculateQuantificationParametersFilter<TInputImage, TOutputImage>::MultiVolumeToVectorVolume(MultiVolumePointerType inputMultiVolume)
{ 
  typename ImageToVectorImageFilterType::Pointer imageToVectorImageFilter = ImageToVectorImageFilterType::New();

  VolumePointerType       volumeTemp;
  VectorVolumePointerType outputVectorVolume;

  MultiVolumeRegionType inputMultiVolumeRegion = inputMultiVolume->GetLargestPossibleRegion();
  MultiVolumeSizeType   inputMultiVolumeSize = inputMultiVolumeRegion.GetSize();
 
  typename MultiVolumeType::IndexType extractStartIndex;
  MultiVolumeSizeType   extractSize;
  MultiVolumeRegionType extractRegion;
 
  extractStartIndex[0] = 0;
  extractStartIndex[1] = 0;
  extractStartIndex[2] = 0;
  extractSize[0] = inputMultiVolumeSize[0];
  extractSize[1] = inputMultiVolumeSize[1];
  extractSize[2] = inputMultiVolumeSize[2];
  extractSize[3] = 0;

  for (int i = 0; i < (int)inputMultiVolumeSize[3]; i++)
    {
    typename ExtractImageFilterType::Pointer extractImageFilter = ExtractImageFilterType::New();
    extractStartIndex[3] = i;
    extractRegion.SetIndex(extractStartIndex);
    extractRegion.SetSize(extractSize);
    extractImageFilter->SetExtractionRegion(extractRegion);
    extractImageFilter->SetInput(inputMultiVolume);
    extractImageFilter->Update();
    extractImageFilter->ReleaseDataFlagOn();
    volumeTemp = extractImageFilter->GetOutput();
    imageToVectorImageFilter->SetNthInput(i, volumeTemp);
    }
  
  imageToVectorImageFilter->Update();
  outputVectorVolume = dynamic_cast<VectorVolumeType *>(imageToVectorImageFilter->GetOutput() );
  
  return outputVectorVolume;
}

template <class TInputImage, class TOutputImage>
void CalculateQuantificationParametersFilter<TInputImage,TOutputImage>
::SetTimeAxis(vcl_vector<float> inputTimeAxis)
{  
  m_timeAxis = new float[inputTimeAxis.size()]();

  for(int i=0; i<(int)inputTimeAxis.size(); ++i)
    {
    m_timeAxis[i] = static_cast<float>(inputTimeAxis[i]);
    }
}

template <class TInputImage, class TOutputImage>
float* CalculateQuantificationParametersFilter<TInputImage,TOutputImage>
::GetTimeAxis()
{
  return m_timeAxis;
}

template <class TInputImage, class TOutputImage>
void CalculateQuantificationParametersFilter<TInputImage,TOutputImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "m_fTol: "                  << m_fTol                         << std::endl;
  os << indent << "m_gTol: "                  << m_gTol                 << std::endl;
  os << indent << "m_xTol: "                  << m_xTol                 << std::endl;
  os << indent << "m_epsilon: "                 << m_epsilon                          << std::endl;
  os << indent << "m_maxIter: "                 << m_maxIter                              << std::endl;
  os << indent << "m_hematocrit: "                << m_hematocrit                             << std::endl;
}

} // end namespace itk

