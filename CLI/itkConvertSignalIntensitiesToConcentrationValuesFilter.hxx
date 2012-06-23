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
  InternalMultiVolumePointerType inputMultiVolume = castFilter->GetOutput();
  typename InputImageType::SizeType inputMultiVolumeSize = inputMultiVolume->GetLargestPossibleRegion().GetSize();
  std::cout << "input MultiVolume Size:"<<inputMultiVolumeSize << std::endl;
  
  InternalVectorVolumePointerType inputVectorVolume = InternalVectorVolumeType::New();
  inputVectorVolume = this->MultiVolumeToVectorVolume(inputMultiVolume);

  // Get S0 Volume
  typedef itk::S0CalculationFilter<InternalVectorVolumeType, InternalVolumeType> S0VolumeFilterType;
  typename S0VolumeFilterType::Pointer S0VolumeFilter = S0VolumeFilterType::New();
  S0VolumeFilter->SetInput(inputVectorVolume);
  S0VolumeFilter->SetS0GradThresh(m_S0GradThresh);
  S0VolumeFilter->Update();

  InternalVolumePointerType S0Volume = S0VolumeFilter->GetOutput();
  
  InternalVolumeIterType S0VolumeIter(S0Volume, S0Volume->GetRequestedRegion() );
  InputMaskIterType      aifMaskVolumeIter(m_AIFMask,m_AIFMask->GetRequestedRegion() );

  inputVectorVolume = this->MultiVolumeToVectorVolume(inputMultiVolume);
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

  this->SetNthOutput(0, this->VectorVolumeToMultiVolume(inputVectorVolume) );
  delete [] concentrationVectorVoxelTemp;
}

template <class TInputImage, class TOutputImage>
typename ConvertSignalIntensitiesToConcentrationValuesFilter<TInputImage, TOutputImage>::InternalMultiVolumePointerType
ConvertSignalIntensitiesToConcentrationValuesFilter<TInputImage, TOutputImage>::VectorVolumeToMultiVolume(
  InternalVectorVolumePointerType inputVectorVolume)
{
  std::cout << std::endl << "VectorVolume To MultiVolume" << std::endl;
  InternalVectorVolumeRegionType inputVectorVolumeRegion = inputVectorVolume->GetLargestPossibleRegion();
  InternalVectorVolumeSizeType   inputVectorVolumeSize = inputVectorVolumeRegion.GetSize();
  InternalVectorVolumeIterType   inputVectorVolumeIter(inputVectorVolume, inputVectorVolume->GetRequestedRegion() );

  InternalMultiVolumePointerType     outputMultiVolume = InternalMultiVolumeType::New();
  InternalMultiVolumeType::IndexType outputMultiVolumeStartIndex;
  outputMultiVolumeStartIndex[0]=0;
  outputMultiVolumeStartIndex[1]=0;
  outputMultiVolumeStartIndex[2]=0;
  outputMultiVolumeStartIndex[3]=0;
  InternalMultiVolumeType::SizeType outputMultiVolumeSize;
  outputMultiVolumeSize[0]=inputVectorVolumeSize[0];
  outputMultiVolumeSize[1]=inputVectorVolumeSize[1];
  outputMultiVolumeSize[2]=inputVectorVolumeSize[2];
  outputMultiVolumeSize[3]=inputVectorVolume->GetNumberOfComponentsPerPixel();

  InternalMultiVolumeType::RegionType outputMultiVolumeRegion;
  outputMultiVolumeRegion.SetSize(outputMultiVolumeSize);
  outputMultiVolumeRegion.SetIndex(outputMultiVolumeStartIndex);
  outputMultiVolume->SetRegions(outputMultiVolumeRegion);
  outputMultiVolume->Allocate();
  outputMultiVolume->FillBuffer(0);

  InternalVectorVoxelType             vectorVoxel;
  InternalVectorVolumeType::IndexType tempVectorVolumeIndex;
  InternalMultiVolumeType::IndexType       outputMultiVolumeIndex;
  while (!inputVectorVolumeIter.IsAtEnd() )
    {
    tempVectorVolumeIndex=inputVectorVolumeIter.GetIndex();
    vectorVoxel = inputVectorVolume->GetPixel(tempVectorVolumeIndex);
    outputMultiVolumeIndex[0]=tempVectorVolumeIndex[0];
    outputMultiVolumeIndex[1]=tempVectorVolumeIndex[1];
    outputMultiVolumeIndex[2]=tempVectorVolumeIndex[2];
    for(int i = 0; i < (int)inputVectorVolume->GetNumberOfComponentsPerPixel(); i++)
      {
      outputMultiVolumeIndex[3]=i;
      outputMultiVolume->SetPixel(outputMultiVolumeIndex, vectorVoxel[i]);
      }
    ++inputVectorVolumeIter;
    }

  return outputMultiVolume;
}

template <class TInputImage, class TOutputImage>
typename ConvertSignalIntensitiesToConcentrationValuesFilter<TInputImage,
                                                             TOutputImage>::InternalVectorVolumePointerType
ConvertSignalIntensitiesToConcentrationValuesFilter<TInputImage, TOutputImage>::MultiVolumeToVectorVolume(
  InternalMultiVolumePointerType inputMultiVolume)
{
  std::cout << std::endl << "MultiVolume To VectorVolume" << std::endl;
  ImageToVectorImageFilterType::Pointer imageToVectorImageFilter = ImageToVectorImageFilterType::New();

  InternalVolumeType::Pointer       volumeTemp;
  InternalVectorVolumeType::Pointer outputVectorVolume;

  InternalMultiVolumeRegionType       inputMultiVolumeRegion = inputMultiVolume->GetLargestPossibleRegion();
  InternalMultiVolumeSizeType         inputMultiVolumeSize = inputMultiVolumeRegion.GetSize();
  InternalMultiVolumeType::IndexType  extractStartIndex;
  InternalMultiVolumeType::SizeType   extractSize;
  InternalMultiVolumeType::RegionType extractRegion;

  for (int i = 0; i < (int)inputMultiVolumeSize[3]; i++)
    {
    ExtractImageFilterType::Pointer extractImageFilter = ExtractImageFilterType::New();
    extractStartIndex[0] = 0;
    extractStartIndex[1] = 0;
    extractStartIndex[2] = 0;
    extractSize[0] = inputMultiVolumeSize[0];
    extractSize[1] = inputMultiVolumeSize[1];
    extractSize[2] = inputMultiVolumeSize[2];
    extractSize[3] = 0;
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
