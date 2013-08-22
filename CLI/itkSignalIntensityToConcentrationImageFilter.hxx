#ifndef _itkSignalIntensityToConcentrationImageFilter_hxx
#define _itkSignalIntensityToConcentrationImageFilter_hxx
#include "itkSignalIntensityToConcentrationImageFilter.h"
#include "itkProgressReporter.h"

namespace itk
{

template <class TInputImage, class TMaskImage, class TOutputImage>
SignalIntensityToConcentrationImageFilter<TInputImage, TMaskImage,
                                                    TOutputImage>::SignalIntensityToConcentrationImageFilter()
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

  OutputImageType* outputVolume = this->GetOutput();
  outputVolume->SetBufferedRegion(inputVectorVolume->GetBufferedRegion());
  outputVolume->Allocate();

  // Get S0 Volume
  typedef SignalIntensityToS0ImageFilter<TInputImage, InternalVolumeType> S0VolumeFilterType;
  typename S0VolumeFilterType::Pointer S0VolumeFilter = S0VolumeFilterType::New();
  S0VolumeFilter->SetInput(inputVectorVolume);
  S0VolumeFilter->SetS0GradThresh(m_S0GradThresh);
  S0VolumeFilter->Update();
  InternalVolumePointerType S0Volume = S0VolumeFilter->GetOutput();
  
  InternalVolumeIterType S0VolumeIter(S0Volume, S0Volume->GetRequestedRegion() );
  InputImageConstIterType inputVectorVolumeIter(inputVectorVolume, 
                                               inputVectorVolume->GetRequestedRegion());
  OutputIterType oit(outputVolume,outputVolume->GetRequestedRegion());

  InputMaskConstIterType aifMaskVolumeIter;
  if (this->GetAIFMask())
    {
    aifMaskVolumeIter = InputMaskConstIterType(this->GetAIFMask(),this->GetAIFMask()->GetRequestedRegion() );
    aifMaskVolumeIter.GoToBegin();
    }

  InputMaskConstIterType roiMaskVolumeIter;
  if (this->GetROIMask())
    {
    roiMaskVolumeIter = InputMaskConstIterType(this->GetROIMask(),this->GetROIMask()->GetRequestedRegion() );
    roiMaskVolumeIter.GoToBegin();
    }

  S0VolumeIter.GoToBegin();
  inputVectorVolumeIter.GoToBegin();
  oit.GoToBegin();

  float * concentrationVectorVoxelTemp =
    new float[(int)inputVectorVolume->GetNumberOfComponentsPerPixel()];
  bool                    isConvert;
  InputPixelType inputVectorVoxel;
  InternalVectorVoxelType vectorVoxel;
  OutputPixelType outputVectorVoxel;

  ProgressReporter progress(this, 0, outputVolume->GetRequestedRegion().GetNumberOfPixels());
  
  // Convert signal intensities to concentration values
  while (!oit.IsAtEnd() )
    {
    inputVectorVoxel = inputVectorVolumeIter.Get();
    vectorVoxel.SetSize(inputVectorVoxel.GetSize());
    vectorVoxel.Fill(0.0);
    vectorVoxel += inputVectorVoxel; // shorthand for a copy/cast

    float T1Pre = m_T1PreTissue;

    // if we have ROI, and T1Pre is set to tissue value, check if we need to skip this voxel
    if( this->GetROIMask() && (this->GetROIMask()->GetBufferedRegion().GetSize()[0])!=0)
      {
      if(!roiMaskVolumeIter.Get())
        {
        T1Pre = 0;
        }
      ++roiMaskVolumeIter;
      }

    // if we have an AIF mask, use blood T1Pre; even if this voxel was unselected by ROI,
    // select it for AIF calculation if non-0
    if( this->GetAIFMask() && (this->GetAIFMask()->GetBufferedRegion().GetSize()[0])!=0)
      {
      if(aifMaskVolumeIter.Get())
        {
        T1Pre = m_T1PreBlood;
        }
      ++aifMaskVolumeIter;
      }

     if(T1Pre)
       {
       isConvert = convert_signal_to_concentration (inputVectorVolume->GetNumberOfComponentsPerPixel(),
                                                     vectorVoxel.GetDataPointer(),
                                                     T1Pre, m_TR, m_FA,
                                                     concentrationVectorVoxelTemp,
                                                     m_RGD_relaxivity,
                                                     S0VolumeIter.Get(),
                                                     m_S0GradThresh);

       // copy the concentration vector to the output
       outputVectorVoxel.SetSize(inputVectorVoxel.GetSize());
       for (typename OutputPixelType::ElementIdentifier i = 0; 
                     i < outputVectorVoxel.GetSize(); ++i)
         {
         outputVectorVoxel[i] 
           = static_cast<typename OutputPixelType::ValueType>(concentrationVectorVoxelTemp[i]);
         }
         oit.Set(outputVectorVoxel);
       }

    ++S0VolumeIter;
    ++inputVectorVolumeIter;
    ++oit;
    progress.CompletedPixel();
    }

  delete [] concentrationVectorVoxelTemp;
}


template <class TInputImage, class TMaskImage, class TOutput>
void SignalIntensityToConcentrationImageFilter<TInputImage, TMaskImage, TOutput>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "T1PreBlood: " << m_T1PreBlood << std::endl;
  os << indent << "T1PreTissue: " << m_T1PreTissue << std::endl;
  os << indent << "TR: " << m_TR << std::endl;
  os << indent << "FA: " << m_FA << std::endl;
  os << indent << "RGD_relaxivity: " << m_RGD_relaxivity << std::endl;
  os << indent << "S0GradThresh: " << m_S0GradThresh << std::endl;
}

} // end namespace itk

#endif
