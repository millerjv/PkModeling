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
  m_AIFMask = InputMaskType::New();
  this->Superclass::SetNumberOfRequiredInputs(1);
  this->Superclass::SetNumberOfRequiredOutputs(1);
  this->Superclass::SetNthOutput(0, OutputImageType::New() );
}

template<class TInputImage, class TOutputImage>
void ConvertSignalIntensitiesToConcentrationValuesFilter<TInputImage, TOutputImage>::GenerateData()
{
  std::cout << "Signal Intensity To Concentration" << std::endl;
  typename CastFilterType::Pointer castFilter = CastFilterType::New();
  castFilter->SetInput(const_cast<InputImageType *>(this->GetInput() ) );
  castFilter->Update();
  InternalQulumePointerType inputQulume = castFilter->GetOutput();
  typename InputImageType::SizeType inputQulumeSize = inputQulume->GetLargestPossibleRegion().GetSize();
  std::cout << "input Qulume Size:"<<inputQulumeSize << std::endl;
  
  InternalVectorVolumePointerType inputVectorVolume = InternalVectorVolumeType::New();
  inputVectorVolume = this->QulumeToVectorVolume(inputQulume);

  // Get S0 Volume
  typedef itk::S0CalculationFilter<InternalVectorVolumeType, InternalVolumeType> S0VolumeFilterType;
  typename S0VolumeFilterType::Pointer S0VolumeFilter = S0VolumeFilterType::New();
  S0VolumeFilter->SetInput(inputVectorVolume);
  S0VolumeFilter->SetS0GradThresh(m_S0GradThresh);
  S0VolumeFilter->Update();

  InternalVolumePointerType S0Volume = S0VolumeFilter->GetOutput();
  
  InternalVolumeIterType S0VolumeIter(S0Volume, S0Volume->GetRequestedRegion() );
  InputMaskIterType      aifMaskVolumeIter(m_AIFMask,m_AIFMask->GetRequestedRegion() );

  inputVectorVolume = this->QulumeToVectorVolume(inputQulume);
  InternalVectorVolumeIterType inputVectorVolumeIter(inputVectorVolume, inputVectorVolume->GetLargestPossibleRegion() );

  aifMaskVolumeIter.GoToBegin();
  S0VolumeIter.GoToBegin();
  inputVectorVolumeIter.GoToBegin();
  float * concentrationVectorVoxelTemp =
    new float[(int)inputVectorVolume->GetNumberOfComponentsPerPixel()];
  bool                    isConvert;
  InternalVectorVoxelType vectorVoxel;
  
  // Convert signal intensities to concentration values, if it is AIF voxel use bloood T1, otherwise use tissue T1
  while (!inputVectorVolumeIter.IsAtEnd() )
    {
    
    vectorVoxel = inputVectorVolumeIter.Get();
    
    if( (m_AIFMask->GetLargestPossibleRegion().GetSize()[0])!=0)
      {
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
      
      isConvert = convert_signal_to_concentration (inputVectorVolume->GetNumberOfComponentsPerPixel(),
                                                   vectorVoxel.GetDataPointer(),
                                                   m_T1PreTissue, m_TR, m_FA,
                                                   concentrationVectorVoxelTemp,
                                                   m_RGD_relaxivity,
                                                   S0VolumeIter.Get(),
                                                   m_S0GradThresh);

      }
    vectorVoxel.SetData(concentrationVectorVoxelTemp,inputVectorVolume->GetNumberOfComponentsPerPixel() );
    inputVectorVolumeIter.Set(vectorVoxel);
    ++S0VolumeIter;
    ++inputVectorVolumeIter;
    }

  this->SetNthOutput(0, this->VectorVolumeToQulume(inputVectorVolume) );
  delete [] concentrationVectorVoxelTemp;
}

template <class TInputImage, class TOutputImage>
typename ConvertSignalIntensitiesToConcentrationValuesFilter<TInputImage, TOutputImage>::InternalQulumePointerType
ConvertSignalIntensitiesToConcentrationValuesFilter<TInputImage, TOutputImage>::VectorVolumeToQulume(
  InternalVectorVolumePointerType inputVectorVolume)
{
  std::cout << std::endl << "VectorVolume To Qulume" << std::endl;
  InternalVectorVolumeRegionType inputVectorVolumeRegion = inputVectorVolume->GetLargestPossibleRegion();
  InternalVectorVolumeSizeType   inputVectorVolumeSize = inputVectorVolumeRegion.GetSize();
  InternalVectorVolumeIterType   inputVectorVolumeIter(inputVectorVolume, inputVectorVolume->GetRequestedRegion() );

  InternalQulumePointerType     outputQulume = InternalQulumeType::New();
  InternalQulumeType::IndexType outputQulumeStartIndex;
  outputQulumeStartIndex[0]=0;
  outputQulumeStartIndex[1]=0;
  outputQulumeStartIndex[2]=0;
  outputQulumeStartIndex[3]=0;
  InternalQulumeType::SizeType outputQulumeSize;
  outputQulumeSize[0]=inputVectorVolumeSize[0];
  outputQulumeSize[1]=inputVectorVolumeSize[1];
  outputQulumeSize[2]=inputVectorVolumeSize[2];
  outputQulumeSize[3]=inputVectorVolume->GetNumberOfComponentsPerPixel();

  InternalQulumeType::RegionType outputQulumeRegion;
  outputQulumeRegion.SetSize(outputQulumeSize);
  outputQulumeRegion.SetIndex(outputQulumeStartIndex);
  outputQulume->SetRegions(outputQulumeRegion);
  outputQulume->Allocate();
  outputQulume->FillBuffer(0);

  InternalVectorVoxelType             vectorVoxel;
  InternalVectorVolumeType::IndexType tempVectorVolumeIndex;
  InternalQulumeType::IndexType       outputQulumeIndex;
  while (!inputVectorVolumeIter.IsAtEnd() )
    {
    tempVectorVolumeIndex=inputVectorVolumeIter.GetIndex();
    vectorVoxel = inputVectorVolume->GetPixel(tempVectorVolumeIndex);
    outputQulumeIndex[0]=tempVectorVolumeIndex[0];
    outputQulumeIndex[1]=tempVectorVolumeIndex[1];
    outputQulumeIndex[2]=tempVectorVolumeIndex[2];
    for(int i = 0; i < (int)inputVectorVolume->GetNumberOfComponentsPerPixel(); i++)
      {
      outputQulumeIndex[3]=i;
      outputQulume->SetPixel(outputQulumeIndex, vectorVoxel[i]);
      }
    ++inputVectorVolumeIter;
    }

  return outputQulume;
}

template <class TInputImage, class TOutputImage>
typename ConvertSignalIntensitiesToConcentrationValuesFilter<TInputImage,
                                                             TOutputImage>::InternalVectorVolumePointerType
ConvertSignalIntensitiesToConcentrationValuesFilter<TInputImage, TOutputImage>::QulumeToVectorVolume(
  InternalQulumePointerType inputQulume)
{
  std::cout << std::endl << "Qulume To VectorVolume" << std::endl;
  ImageToVectorImageFilterType::Pointer imageToVectorImageFilter = ImageToVectorImageFilterType::New();

  InternalVolumeType::Pointer       volumeTemp;
  InternalVectorVolumeType::Pointer outputVectorVolume;

  InternalQulumeRegionType       inputQulumeRegion = inputQulume->GetLargestPossibleRegion();
  InternalQulumeSizeType         inputQulumeSize = inputQulumeRegion.GetSize();
  InternalQulumeType::IndexType  extractStartIndex;
  InternalQulumeType::SizeType   extractSize;
  InternalQulumeType::RegionType extractRegion;

  for (int i = 0; i < (int)inputQulumeSize[3]; i++)
    {
    ExtractImageFilterType::Pointer extractImageFilter = ExtractImageFilterType::New();
    extractStartIndex[0] = 0;
    extractStartIndex[1] = 0;
    extractStartIndex[2] = 0;
    extractSize[0] = inputQulumeSize[0];
    extractSize[1] = inputQulumeSize[1];
    extractSize[2] = inputQulumeSize[2];
    extractSize[3] = 0;
    extractStartIndex[3] = i;
    extractRegion.SetIndex(extractStartIndex);
    extractRegion.SetSize(extractSize);
    extractImageFilter->SetExtractionRegion(extractRegion);
    extractImageFilter->SetInput(inputQulume);
    extractImageFilter->Update();
    extractImageFilter->ReleaseDataFlagOn();
    volumeTemp = extractImageFilter->GetOutput();
    imageToVectorImageFilter->SetNthInput(i, volumeTemp);
    }

  imageToVectorImageFilter->ReleaseDataFlagOn();
  imageToVectorImageFilter->Update();
  outputVectorVolume = imageToVectorImageFilter->GetOutput();

  return outputVectorVolume;
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
