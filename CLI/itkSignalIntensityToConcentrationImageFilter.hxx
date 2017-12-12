#ifndef _itkSignalIntensityToConcentrationImageFilter_hxx
#define _itkSignalIntensityToConcentrationImageFilter_hxx
#include "itkSignalIntensityToConcentrationImageFilter.h"
#include "itkProgressReporter.h"

namespace itk
{

template <class TInputImage, class TMaskImage, class TOutputImage>
SignalIntensityToConcentrationImageFilter<TInputImage, TMaskImage, TOutputImage>::SignalIntensityToConcentrationImageFilter()
{
  m_T1PreTissue = 0.0f;
  m_T1PreBlood = m_T1PreTissue;
  m_TR = 0.0f;
  m_FA = 0.0f;
  m_RGD_relaxivity = 4.9E-3f;
  m_S0GradThresh = 15.0f;
  this->SetNumberOfRequiredInputs(1);
}


template<class TInputImage, class TMaskImage, class TOutputImage>
void SignalIntensityToConcentrationImageFilter<TInputImage, TMaskImage, TOutputImage>::GenerateData()
{
  const InputImageType* inputVectorVolume = this->GetInput();

  OutputImageType* outputVolume = this->GetAllocatedOutVolume(inputVectorVolume);
  InternalVolumePointerType S0Volume = this->GetS0Image(inputVectorVolume);
  
  InternalVolumeIterType S0VolumeIter(S0Volume, S0Volume->GetRequestedRegion());
  S0VolumeIter.GoToBegin();
  InputImageConstIterType inputVectorVolumeIter(inputVectorVolume, inputVectorVolume->GetRequestedRegion());
  inputVectorVolumeIter.GoToBegin();
  OutputIterType outVolumeIter(outputVolume, outputVolume->GetRequestedRegion());
  outVolumeIter.GoToBegin();
  T1PreValueIterator t1PreIter(this->GetROIMask(), this->GetAIFMask(), this->GetT1Map(), this->m_T1PreTissue, this->m_T1PreBlood);

  float* concentrationVectorVoxelTemp = new float[(int)inputVectorVolume->GetNumberOfComponentsPerPixel()];
  OutputPixelType outputVectorVoxel;

  ProgressReporter progress(this, 0, outputVolume->GetRequestedRegion().GetNumberOfPixels());
  
  // Convert signal intensities to concentration values
  while (!outVolumeIter.IsAtEnd())
  {
    InternalVectorVoxelType vectorVoxel = this->convertToInternalVectorVoxel(inputVectorVolumeIter.Get());
    outputVectorVoxel.SetSize(vectorVoxel.GetSize());
    float T1Pre = t1PreIter.Get();
    if (T1Pre)
    {
      bool isConvert = convert_signal_to_concentration(inputVectorVolume->GetNumberOfComponentsPerPixel(),
                                                       vectorVoxel.GetDataPointer(),
                                                       T1Pre, m_TR, m_FA,
                                                       concentrationVectorVoxelTemp,
                                                       m_RGD_relaxivity,
                                                       S0VolumeIter.Get(),
                                                       m_S0GradThresh);

      for (typename OutputPixelType::ElementIdentifier i = 0; i < outputVectorVoxel.GetSize(); ++i)
      {
        outputVectorVoxel[i] = static_cast<typename OutputPixelType::ValueType>(concentrationVectorVoxelTemp[i]);
      }
    }
    else
    {
      outputVectorVoxel.Fill(static_cast<typename OutputPixelType::ValueType>(0.0));
    }
    outVolumeIter.Set(outputVectorVoxel);

    ++S0VolumeIter;
    ++inputVectorVolumeIter;
    ++outVolumeIter;
    ++t1PreIter;

    progress.CompletedPixel();
    }

  delete [] concentrationVectorVoxelTemp;
}


template<class TInputImage, class TMaskImage, class TOutputImage>
typename SignalIntensityToConcentrationImageFilter<TInputImage, TMaskImage, TOutputImage>::OutputImageType*
SignalIntensityToConcentrationImageFilter<TInputImage, TMaskImage, TOutputImage>::GetAllocatedOutVolume(const InputImageType* inputVectorVolume)
{
  OutputImageType* outputVolume = this->GetOutput();
  outputVolume->SetBufferedRegion(inputVectorVolume->GetBufferedRegion());
  outputVolume->Allocate();
  return outputVolume;
}


template<class TInputImage, class TMaskImage, class TOutputImage>
typename SignalIntensityToConcentrationImageFilter<TInputImage, TMaskImage, TOutputImage>::InternalVolumePointerType
SignalIntensityToConcentrationImageFilter<TInputImage, TMaskImage, TOutputImage>::GetS0Image(const InputImageType* inputVectorVolume)
{
  typedef SignalIntensityToS0ImageFilter<TInputImage, InternalVolumeType> S0VolumeFilterType;
  typename S0VolumeFilterType::Pointer S0VolumeFilter = S0VolumeFilterType::New();
  S0VolumeFilter->SetInput(inputVectorVolume);
  S0VolumeFilter->SetS0GradThresh(m_S0GradThresh);
  S0VolumeFilter->SetBATCalculationMode(m_BATCalculationMode);
  S0VolumeFilter->SetconstantBAT(m_constantBAT);
  S0VolumeFilter->Update();
  InternalVolumePointerType S0Volume = S0VolumeFilter->GetOutput();
  return S0Volume;
}


template<class TInputImage, class TMaskImage, class TOutputImage>
typename SignalIntensityToConcentrationImageFilter<TInputImage, TMaskImage, TOutputImage>::InternalVectorVoxelType
SignalIntensityToConcentrationImageFilter<TInputImage, TMaskImage, TOutputImage>::convertToInternalVectorVoxel(const InputPixelType& inputVectorVoxel)
{
  InternalVectorVoxelType vectorVoxel;
  vectorVoxel.SetSize(inputVectorVoxel.GetSize());
  vectorVoxel.Fill(0.0);
  vectorVoxel += inputVectorVoxel; // shorthand for a copy/cast
  return vectorVoxel;
}



template <class TInputImage, class TMaskImage, class TOutput>
void SignalIntensityToConcentrationImageFilter<TInputImage, TMaskImage, TOutput>::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "T1PreBlood: " << m_T1PreBlood << std::endl;
  os << indent << "T1PreTissue: " << m_T1PreTissue << std::endl;
  os << indent << "TR: " << m_TR << std::endl;
  os << indent << "FA: " << m_FA << std::endl;
  os << indent << "RGD_relaxivity: " << m_RGD_relaxivity << std::endl;
  os << indent << "S0GradThresh: " << m_S0GradThresh << std::endl;
}



//==================== T1PreValueIterator internal helper class ====================

template<class TInputImage, class TMaskImage, class TOutputImage>
SignalIntensityToConcentrationImageFilter<TInputImage, TMaskImage, TOutputImage>::T1PreValueIterator::T1PreValueIterator(const InputMaskType* roiMask,
                                                                                                                         const InputMaskType* aifMask, 
                                                                                                                         const InputMaskType* t1Map, 
                                                                                                                         float t1PreTissue, 
                                                                                                                         float t1PreBlood)
{
  this->m_T1PreTissue = t1PreTissue;
  this->m_T1PreBlood = t1PreBlood;
  this->roiMaskVolumeIter = this->getNewConstMaskIterOrNull(roiMask);
  this->aifMaskVolumeIter = this->getNewConstMaskIterOrNull(aifMask);
  this->T1MapVolumeIter = this->getNewConstMaskIterOrNull(t1Map);
}


template<class TInputImage, class TMaskImage, class TOutputImage>
SignalIntensityToConcentrationImageFilter<TInputImage, TMaskImage, TOutputImage>::T1PreValueIterator::~T1PreValueIterator()
{
  delete this->roiMaskVolumeIter;
  delete this->aifMaskVolumeIter;
  delete this->T1MapVolumeIter;
}


template<class TInputImage, class TMaskImage, class TOutputImage>
float
SignalIntensityToConcentrationImageFilter<TInputImage, TMaskImage, TOutputImage>::T1PreValueIterator::Get()
{
  float T1Pre = T1MapVolumeIter ? T1MapVolumeIter->Get() : m_T1PreTissue;
  if (aifMaskVolumeIter && aifMaskVolumeIter->Get()) {
    T1Pre = T1MapVolumeIter ? T1MapVolumeIter->Get() : m_T1PreBlood;
  }
  else if (roiMaskVolumeIter && !roiMaskVolumeIter->Get()) {
    T1Pre = 0;
  }
  return T1Pre;
}


template<class TInputImage, class TMaskImage, class TOutputImage>
void
SignalIntensityToConcentrationImageFilter<TInputImage, TMaskImage, TOutputImage>::T1PreValueIterator::GoToBegin()
{
  if (this->roiMaskVolumeIter) {
    this->roiMaskVolumeIter->GoToBegin();
  }
  if (aifMaskVolumeIter) {
    this->aifMaskVolumeIter->GoToBegin();
  }
  if (T1MapVolumeIter) {
    this->T1MapVolumeIter->GoToBegin();
  }
}


template<class TInputImage, class TMaskImage, class TOutputImage>
typename SignalIntensityToConcentrationImageFilter<TInputImage, TMaskImage, TOutputImage>::T1PreValueIterator&
SignalIntensityToConcentrationImageFilter<TInputImage, TMaskImage, TOutputImage>::T1PreValueIterator::operator++()
{
  if (this->roiMaskVolumeIter) {
    ++(*(this->roiMaskVolumeIter));
  }
  if (aifMaskVolumeIter) {
    ++(*(this->aifMaskVolumeIter));
  }
  if (T1MapVolumeIter) {
    ++(*(this->T1MapVolumeIter));
  }
  return *this;
}


template<class TInputImage, class TMaskImage, class TOutputImage>
typename SignalIntensityToConcentrationImageFilter<TInputImage, TMaskImage, TOutputImage>::InputMaskConstIterType*
SignalIntensityToConcentrationImageFilter<TInputImage, TMaskImage, TOutputImage>::T1PreValueIterator::getNewConstMaskIterOrNull(const InputMaskType* inMask)
{
  InputMaskConstIterType* maskVolumeIter = NULL;
  if (inMask && (inMask->GetBufferedRegion().GetSize()[0] != 0))
  {
    maskVolumeIter = new InputMaskConstIterType(inMask, inMask->GetRequestedRegion());
    maskVolumeIter->GoToBegin();
  }
  return maskVolumeIter;
}





} // end namespace itk
#endif
