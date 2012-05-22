#ifndef _itkConvertSignalIntensitiesToConcentrationValuesFilter_hxx
#define _itkConvertSignalIntensitiesToConcentrationValuesFilter_hxx
#include "itkConvertSignalIntensitiesToConcentrationValuesFilter.h"
#include "itkImageFileReader.h"
#include "itkCastImageFilter.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
ConvertSignalIntensitiesToConcentrationValuesFilter<TInputImage, TOutputImage>::ConvertSignalIntensitiesToConcentrationValuesFilter()
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
	this->Superclass::SetNthOutput(0, OutputImageType::New());	
}

template<class TInputImage, class TOutputImage>
void ConvertSignalIntensitiesToConcentrationValuesFilter<TInputImage, TOutputImage>::GenerateData()
{
	std::cout << "Signal Intensity To Concentration" << std::endl;
	CastFilterType::Pointer castFilter = CastFilterType::New();
	castFilter->SetInput(const_cast<InputImageType *>(this->GetInput()));
	castFilter->Update();
	typename InternalQulumePointerType inputQulume = castFilter->GetOutput();
	InputImageType::SizeType inputQulumeSize = inputQulume->GetLargestPossibleRegion().GetSize();
	std::cout << "input Qulume Size:"<<inputQulumeSize << std::endl;

	//Qulume to VectorVolume
	typename InternalVectorVolumePointerType inputVectorVolume = InternalVectorVolumeType::New();
    inputVectorVolume = this->QulumeToVectorVolume(inputQulume);
	//unsigned int vectorSize = (unsigned int)inputVectorVolume->GetNumberOfComponentsPerPixel();

	// Get S0 Volume
	//clock_t begin=clock();
	typedef itk::S0CalculationFilter<InternalVectorVolumeType, InternalVolumeType>  S0VolumeFilterType;
	typename S0VolumeFilterType::Pointer S0VolumeFilter = S0VolumeFilterType::New();
	S0VolumeFilter->SetInput(inputVectorVolume);
	S0VolumeFilter->SetS0GradThresh(m_S0GradThresh);
	S0VolumeFilter->Update();		
	
	typename InternalVolumePointerType S0Volume = S0VolumeFilter->GetOutput();		
	typedef itk::ImageFileWriter<InternalVolumeType> SynWriterType;
	SynWriterType::Pointer synWriter = SynWriterType::New();
	synWriter->SetFileName("D:/Codes/Slicer4/Modules/CLI/SignalIntensitiesToConcentrationValues/Data/DukeData/SyntheticDukeSmallQulumeS0.nrrd");
	synWriter->SetInput(S0Volume);
	synWriter->Update();

	typename InternalVolumeIterType S0VolumeIter(S0Volume, S0Volume->GetRequestedRegion());	
	typename InputMaskIterType	aifMaskVolumeIter(m_AIFMask,m_AIFMask->GetRequestedRegion());

	inputVectorVolume = this->QulumeToVectorVolume(inputQulume);
	typename InternalVectorVolumeIterType inputVectorVolumeIter(inputVectorVolume, inputVectorVolume->GetLargestPossibleRegion());	
		
	aifMaskVolumeIter.GoToBegin();
	S0VolumeIter.GoToBegin();
	inputVectorVolumeIter.GoToBegin();
	float * concentrationVectorVoxelTemp = new float[(int)inputVectorVolume->GetNumberOfComponentsPerPixel()];
	bool isConvert;
	typename InternalVectorVoxelType vectorVoxel;	
	//std::cout << "before while:"<< std::endl;
		
	while (!inputVectorVolumeIter.IsAtEnd())
    {
		//std::cout << "GET VECTORVOXEL" << std::endl;
		vectorVoxel = inputVectorVolumeIter.Get();	
		//std::cout << "GoT VECTORVOXEL"<<m_AIFMask->GetLargestPossibleRegion().GetSize()[0] << std::endl;
		if((m_AIFMask->GetLargestPossibleRegion().GetSize()[0])!=0)
		{
			if(aifMaskVolumeIter.Get()!=0)
			{
				//std::cout << "AIF MASK !0" << std::endl;				
				//std::cerr<<"aif_bat:"<<aif_BATIndex<<std::endl;
				//std::cerr<<"aif_FIRSTPEAK:"<<aif_FirstPeakIndex<<std::endl;
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
				//std::cout << "AIF MASK !0 and not aif voxel" << std::endl;
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
			//std::cout << "AIF MASK 0" << std::endl;
			isConvert = convert_signal_to_concentration (inputVectorVolume->GetNumberOfComponentsPerPixel(), 
                                      vectorVoxel.GetDataPointer(), 
                                      m_T1PreTissue, m_TR, m_FA,
                                      concentrationVectorVoxelTemp,
                                      m_RGD_relaxivity,
                                      S0VolumeIter.Get(),
                                      m_S0GradThresh);

			//vectorVoxelIndex = inputVectorVolumeIter.GetIndex();
			//
			//if((vectorVoxelIndex[0]==45)&&(vectorVoxelIndex[1]==79)&&(vectorVoxelIndex[2]==9))
			//{
			//	float * tempConcentration = const_cast<float *>(vectorVoxel.GetDataPointer());
			//	std::cout<<"S0:"<<S0VolumeIter.Get()<<std::endl;
			//	for(int itemp = 0;itemp<(inputVectorVolume->GetNumberOfComponentsPerPixel());itemp++)
			//	{
			//		std::cout<<"i, original signal, concentration:"<<itemp<<","<<tempConcentration[itemp]<<","<<concentrationVectorVoxelTemp[itemp]<<std::endl;
			//	}
			//	//std::cout<<"tempAUC:"<<tempAUC<<std::endl;
			//}
		}
		vectorVoxel.SetData(concentrationVectorVoxelTemp,inputVectorVolume->GetNumberOfComponentsPerPixel());
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
typename ConvertSignalIntensitiesToConcentrationValuesFilter<TInputImage, TOutputImage>::InternalQulumePointerType 
ConvertSignalIntensitiesToConcentrationValuesFilter<TInputImage, TOutputImage>::VectorVolumeToQulume(InternalVectorVolumePointerType inputVectorVolume)
{		
	std::cout << std::endl << "VectorVolume To Qulume" << std::endl;
	InternalVectorVolumeRegionType inputVectorVolumeRegion = inputVectorVolume->GetLargestPossibleRegion();
	InternalVectorVolumeSizeType inputVectorVolumeSize = inputVectorVolumeRegion.GetSize();
	InternalVectorVolumeIterType inputVectorVolumeIter(inputVectorVolume, inputVectorVolume->GetRequestedRegion());

	InternalQulumePointerType outputQulume = InternalQulumeType::New();
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
	
	InternalVectorVoxelType vectorVoxel;
	InternalVectorVolumeType::IndexType tempVectorVolumeIndex; 
	InternalQulumeType::IndexType outputQulumeIndex; 
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
typename ConvertSignalIntensitiesToConcentrationValuesFilter<TInputImage, TOutputImage>::InternalVectorVolumePointerType 
ConvertSignalIntensitiesToConcentrationValuesFilter<TInputImage, TOutputImage>::QulumeToVectorVolume(InternalQulumePointerType inputQulume)
{		
	std::cout << std::endl << "Qulume To VectorVolume" << std::endl;
	ImageToVectorImageFilterType::Pointer imageToVectorImageFilter = ImageToVectorImageFilterType::New();
	
	InternalVolumeType::Pointer volumeTemp;
    InternalVectorVolumeType::Pointer outputVectorVolume;
	
	InternalQulumeRegionType inputQulumeRegion = inputQulume->GetLargestPossibleRegion();
	InternalQulumeSizeType inputQulumeSize = inputQulumeRegion.GetSize();
	InternalQulumeType::IndexType extractStartIndex;
	InternalQulumeType::SizeType extractSize;
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
  os << indent << "T1PreBlood: "									<< m_T1PreBlood                                 << std::endl;
  os << indent << "T1PreTissue: "									<< m_T1PreTissue                                 << std::endl;
  os << indent << "TR: "									<< m_TR									<< std::endl;
  os << indent << "FA: "									<< m_FA									<< std::endl;
  os << indent << "RGD_relaxivity: "                        << m_RGD_relaxivity                        << std::endl;
  os << indent << "S0GradThresh: "                          << m_S0GradThresh                          << std::endl;
}

} // end namespace itk

#endif