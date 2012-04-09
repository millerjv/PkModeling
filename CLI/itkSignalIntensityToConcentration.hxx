#ifndef _itkSignalIntensityToConcentration_hxx
#define _itkSignalIntensityToConcentration_hxx
#include "itkSignalIntensityToConcentration.h"
#include <time.h>
#include <stdio.h>
#include <iostream>
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
SignalIntensityToConcentration<TInputImage, TOutputImage>::SignalIntensityToConcentration()
{	
	m_T1Pre = 0.0f;
	m_TR = 0.0f;
	m_FA = 0.0f;
	m_RGD_relaxivity = 4.9E-3f;
	m_S0GradThresh = 15.0f;
	this->Superclass::SetNumberOfRequiredInputs(1);
	this->Superclass::SetNumberOfRequiredOutputs(1);
	this->Superclass::SetNthOutput(0, OutputImageType::New());	
}

template<class TInputImage, class TOutputImage>
void SignalIntensityToConcentration<TInputImage, TOutputImage>::GenerateData()
{
	std::cerr << std::endl << "Signal Intensity To Concentration" << std::endl;
	InternalQulumePointerType inputQulume = const_cast<InternalQulumeType *>(this->GetInput());
	
	//Qulume to VectorVolume
	InternalVectorVolumePointerType inputVectorVolume = InternalVectorVolumeType::New();
    inputVectorVolume = this->QulumeToVectorVolume(inputQulume);
	unsigned int vectorSize = (unsigned int)inputVectorVolume->GetNumberOfComponentsPerPixel();

	// Get S0 Volume
	//clock_t begin=clock();
	typedef itk::S0ForTimeSeriesInQulume<InternalVectorVolumeType, InternalVolumeType>  S0VolumeFilterType;
	typename S0VolumeFilterType::Pointer S0VolumeFilter = S0VolumeFilterType::New();
	S0VolumeFilter->SetInput(inputVectorVolume);
	S0VolumeFilter->SetS0GradThresh(m_S0GradThresh);
	S0VolumeFilter->Update();		
	/*clock_t end=clock();
	std::cerr << "Time elapsed: " << double(std::difftime(end,begin)) << " ms"<< std::endl;
*/	
	/*typedef itk::ImageFileWriter<InternalVolumeType>							  ImageWriterType;	
	typename ImageWriterType::Pointer s0Writer = ImageWriterType::New();
	s0Writer->SetInput(S0VolumeFilter->GetOutput());
	s0Writer->SetFileName("D:/Codes/Slicer4/Testing/Data/Input/s0volumem.nrrd");
	s0Writer->Update();*/

	InternalVolumePointerType S0Volume = S0VolumeFilter->GetOutput();		
	InternalVolumeIterType S0VolumeIter(S0Volume, S0Volume->GetRequestedRegion());		
	inputVectorVolume = this->QulumeToVectorVolume(inputQulume);
	InternalVectorVolumeIterType inputVectorVolumeIter(inputVectorVolume, inputVectorVolume->GetLargestPossibleRegion());	
	
	S0VolumeIter.GoToBegin();
	inputVectorVolumeIter.GoToBegin();
	float * concentrationVectorVoxelTemp = new float[vectorSize];
	bool isConvert;
	InternalVectorVoxelType vectorVoxel;	
		
	while (!inputVectorVolumeIter.IsAtEnd())
    {
		vectorVoxel = inputVectorVolumeIter.Get();	
		isConvert = convert_signal_to_concentration (vectorSize, 
                                      vectorVoxel.GetDataPointer(), 
                                      m_T1Pre, m_TR, m_FA,
                                      concentrationVectorVoxelTemp,
                                      m_RGD_relaxivity,
                                      S0VolumeIter.Get(),
                                      m_S0GradThresh);
		vectorVoxel.SetData(concentrationVectorVoxelTemp,vectorSize);
		inputVectorVolumeIter.Set(vectorVoxel);
		++S0VolumeIter;
		++inputVectorVolumeIter;	
	}	
	
	this->SetNthOutput(0, this->VectorVolumeToQulume(inputVectorVolume));	

	//typedef itk::S0ForTimeSeriesInQulume<InternalVectorVolumeType, InternalVolumeType>  S0ConcentrationVolumeFilterType;
	//typename S0ConcentrationVolumeFilterType::Pointer S0ConcentrationVolumeFilter = S0ConcentrationVolumeFilterType::New();
	//S0ConcentrationVolumeFilter->SetInput(inputVectorVolume);
	//S0ConcentrationVolumeFilter->SetS0GradThresh(m_S0GradThresh);
	//S0ConcentrationVolumeFilter->Update();
	//	
	//typename ImageWriterType::Pointer s0ConcentrationWriter = ImageWriterType::New();
	//s0ConcentrationWriter->SetInput(S0ConcentrationVolumeFilter->GetOutput());
	//s0ConcentrationWriter->SetFileName("D:/Codes/Slicer4/Testing/Data/Input/s0ConcentrationVolume.nrrd");
	//s0ConcentrationWriter->Update();
	
	delete [] concentrationVectorVoxelTemp;
}

template <class TInputImage, class TOutputImage>
typename SignalIntensityToConcentration<TInputImage, TOutputImage>::InternalQulumePointerType 
SignalIntensityToConcentration<TInputImage, TOutputImage>::VectorVolumeToQulume(InternalVectorVolumePointerType inputVectorVolume)
{		
	std::cout << std::endl << "VectorVolume To Qulume" << std::endl;
	InternalVectorVolumeRegionType inputVectorVolumeRegion = inputVectorVolume->GetLargestPossibleRegion();
	InternalVectorVolumeSizeType inputVectorVolumeSize = inputVectorVolumeRegion.GetSize();
	InternalVectorVolumeIterType inputVectorVolumeIter(inputVectorVolume, inputVectorVolume->GetRequestedRegion());

	InternalQulumePointerType outputQulume = InternalQulumeType::New();
	typename InternalQulumeType::IndexType outputQulumeStartIndex; 
	outputQulumeStartIndex[0]=0; 
	outputQulumeStartIndex[1]=0;
	outputQulumeStartIndex[2]=0;	
	outputQulumeStartIndex[3]=0;
	typename InternalQulumeType::SizeType outputQulumeSize;
    outputQulumeSize[0]=inputVectorVolumeSize[0]; 
	outputQulumeSize[1]=inputVectorVolumeSize[1];
	outputQulumeSize[2]=inputVectorVolumeSize[2];
	outputQulumeSize[3]=inputVectorVolume->GetNumberOfComponentsPerPixel();
	
	typename InternalQulumeType::RegionType outputQulumeRegion;
    outputQulumeRegion.SetSize(outputQulumeSize);  
	outputQulumeRegion.SetIndex(outputQulumeStartIndex);    
	outputQulume->SetRegions(outputQulumeRegion);  	
	outputQulume->Allocate(); 		
	outputQulume->FillBuffer(0);
	
	InternalVectorVoxelType vectorVoxel;
	typename InternalVectorVolumeType::IndexType tempVectorVolumeIndex; 
	typename InternalQulumeType::IndexType outputQulumeIndex; 
	while (!inputVectorVolumeIter.IsAtEnd())
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
typename SignalIntensityToConcentration<TInputImage, TOutputImage>::InternalVectorVolumePointerType 
SignalIntensityToConcentration<TInputImage, TOutputImage>::QulumeToVectorVolume(InternalQulumePointerType inputQulume)
{		
	std::cout << std::endl << "Qulume To VectorVolume" << std::endl;
	typename ImageToVectorImageFilterType::Pointer imageToVectorImageFilter = ImageToVectorImageFilterType::New();
	
	InternalVolumePointerType volumeTemp;
    InternalVectorVolumePointerType outputVectorVolume;
	
	InternalQulumeRegionType inputQulumeRegion = inputQulume->GetLargestPossibleRegion();
	InternalQulumeSizeType inputQulumeSize = inputQulumeRegion.GetSize();
	typename InternalQulumeType::IndexType extractStartIndex;
	InternalQulumeSizeType extractSize;
	InternalQulumeRegionType extractRegion;
	typename ExtractImageFilterType::Pointer extractImageFilter;
	extractStartIndex[0] = 0; 
	extractStartIndex[1] = 0; 
	extractStartIndex[2] = 0; 
	extractSize[0] = inputQulumeSize[0]; 
	extractSize[1] = inputQulumeSize[1]; 
	extractSize[2] = inputQulumeSize[2]; 
	extractSize[3] = 0;		 

    for (int i = 0; i < (int)inputQulumeSize[3]; i++) 
	{
		extractImageFilter = ExtractImageFilterType::New();			
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
void SignalIntensityToConcentration<TInputImage, TOutput>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "T1Pre: "									<< m_T1Pre                                 << std::endl;
  os << indent << "TR: "									<< m_TR									<< std::endl;
  os << indent << "FA: "									<< m_FA									<< std::endl;
  os << indent << "RGD_relaxivity: "                        << m_RGD_relaxivity                        << std::endl;
  os << indent << "S0GradThresh: "                          << m_S0GradThresh                          << std::endl;
}

} // end namespace itk

#endif