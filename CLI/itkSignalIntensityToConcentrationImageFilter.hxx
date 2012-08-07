#ifndef _itkConvertSignalIntensitiesToConcentrationValuesFilter_hxx
#define _itkConvertSignalIntensitiesToConcentrationValuesFilter_hxx
#include "itkConvertSignalIntensitiesToConcentrationValuesFilter.h"
#include "itkImageFileReader.h"
#include "itkCastImageFilter.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
ConvertSignalIntensitiesToConcentrationValuesFilter<TInputImage,
                                                    TOutputImage>::ConvertSignalIntensitiesToConcentrationValuesFilter()
{
  m_T1PreTissue = 0.0f;
  m_T1PreBlood = m_T1PreTissue;
  m_TR = 0.0f;
  m_FA = 0.0f;
  m_RGD_relaxivity = 4.9E-3f;
  m_S0GradThresh = 15.0f;
  this->SetNumberOfRequiredInputs(1);
}

template<class TInputImage, class TOutputImage>
void ConvertSignalIntensitiesToConcentrationValuesFilter<TInputImage, TOutputImage>::GenerateData()
{
  std::cout << "Signal Intensity To Concentration" << std::endl;

  InputImageType* inputVectorVolume = dynamic_cast<InputImageType*>(this->GetInput());

  OutputImageType* outputVolume = this->GetOutput();
  outputVolume->SetBufferedRegion(inputVectorVolume->GetBufferedRegion());
  outputVolume->Allocate();

  // Get S0 Volume
  typedef itk::S0CalculationFilter<TInputImage, InternalVolumeType> S0VolumeFilterType;
  typename S0VolumeFilterType::Pointer S0VolumeFilter = S0VolumeFilterType::New();
  S0VolumeFilter->SetInput(inputVectorVolume);
  S0VolumeFilter->SetS0GradThresh(m_S0GradThresh);
  S0VolumeFilter->Update();
  InternalVolumePointerType S0Volume = S0VolumeFilter->GetOutput();
  
  InternalVolumeIterType S0VolumeIter(S0Volume, S0Volume->GetRequestedRegion() );
  InternalVectorVolumeIterType inputVectorVolumeIter(inputVectorVolume, 
                                                     inputVectorVolume->GetRequestedRegion() );
  OutputIterType oit(outputVolume, inputVectorVolume->GetRequestedRegion());

  InputMaskIterType aifMaskVolumeIter;
  if (this->GetAIFMask())
    {
    aifMaskVolumeIter = InputMaskIterType(this->GetAIFMask(),this->GetAIFMask()->GetRequestedRegion() );
    aifMaskVolumeIter.GoToBegin();
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
  
  // Convert signal intensities to concentration values
  while (!inputVectorVolumeIter.IsAtEnd() )
    {
    inputVectorVoxel = inputVectorVolumeIter.Get();
    vectorVoxel.SetSize(inputVectorVoxel.GetSize());
    vectorVoxel.Fill(0.0);
    vectorVoxel += inputVectorVoxel; // shorthand for a copy/cast

    // if we have an AIF mask
    if( this->GetAIFMask() && (this->GetAIFMask()->GetLargestPossibleRegion().GetSize()[0])!=0)
      {
      // if the mask is set, use blood T1, else use tissue T1
      if(aifMaskVolumeIter.Get()!=0)
        {
        isConvert = convert_signal_to_concentration (inputVectorVolume->GetNumberOfComponentsPerPixel(),
                                                     vectorVoxel.GetDataPointer(),
                                                     m_T1PreBlood, m_TR, m_FA,
                                                     concentrationVectorVoxelTemp,
                                                     m_RGD_relaxivity,
                                                     S0VolumeIter.Get(),
                                                     m_S0GradThresh);
        }
      else
        {
        isConvert = convert_signal_to_concentration (inputVectorVolume->GetNumberOfComponentsPerPixel(),
                                                     vectorVoxel.GetDataPointer(),
                                                     m_T1PreTissue, m_TR, m_FA,
                                                     concentrationVectorVoxelTemp,
                                                     m_RGD_relaxivity,
                                                     S0VolumeIter.Get(),
                                                     m_S0GradThresh);
        }
      ++aifMaskVolumeIter;
      }
    else
      {
      // no AIF, use tissue
      isConvert = convert_signal_to_concentration (inputVectorVolume->GetNumberOfComponentsPerPixel(),
                                                   vectorVoxel.GetDataPointer(),
                                                   m_T1PreTissue, m_TR, m_FA,
                                                   concentrationVectorVoxelTemp,
                                                   m_RGD_relaxivity,
                                                   S0VolumeIter.Get(),
                                                   m_S0GradThresh);

      }

    // copy the concentration vector to the output
    outputVectorVoxel.SetSize(inputVectorVoxel.GetSize());
    for (typename OutputPixelType::ElementIdentifer i = 0; 
         i < outputVectorVoxel.GetSize(); ++i)
      {
      outputVectorVoxel[i] 
        = static_cast<typename OutputPixelType::ValueType>(concentrationVectorVoxelTemp[i]);
      }
    oit.Set(outputVectorVoxel);

    ++S0VolumeIter;
    ++inputVectorVolumeIter;
    ++oit;
    }

  delete [] concentrationVectorVoxelTemp;
}


template <class TInputImage, class TOutput>
void ConvertSignalIntensitiesToConcentrationValuesFilter<TInputImage, TOutput>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "T1PreBlood: "                  << m_T1PreBlood                                 << std::endl;
  os << indent << "T1PreTissue: "                 << m_T1PreTissue                                 << std::endl;
  os << indent << "TR: "                  << m_TR                 << std::endl;
  os << indent << "FA: "                  << m_FA                 << std::endl;
  os << indent << "RGD_relaxivity: "                        << m_RGD_relaxivity                        << std::endl;
  os << indent << "S0GradThresh: "                          << m_S0GradThresh                          << std::endl;
}

} // end namespace itk

#endif
