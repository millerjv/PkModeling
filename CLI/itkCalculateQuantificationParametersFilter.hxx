#ifndef _itkCalculateQuantificationParameters_hxx
#define _itkCalculateQuantificationParameters_hxx
#endif

#include "itkImageRegionConstIterator.h"
#include "itkImageRegionIterator.h"
#include "vnl/vnl_math.h"
#include "itkCastImageFilter.h"
#include "itkCalculateQuantificationParameters.h"

#include "itkBinaryThresholdImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkImageFileReader.h"
#include "itkLevenbergMarquardtOptimizer.h"
#include "PkSolver.h"
#include <ostream>
#include "stdlib.h"
#include "stdio.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
CalculateQuantificationParameters<TInputImage,TOutputImage>::CalculateQuantificationParameters()
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
	this->Superclass::SetNthOutput(0, VolumeType::New());	
	this->Superclass::SetNthOutput(1, VolumeType::New());	
	this->Superclass::SetNthOutput(2, VolumeType::New());	
	this->Superclass::SetNthOutput(3, VolumeType::New());	
}

template<class TInputImage, class TOutputImage>
void CalculateQuantificationParameters<TInputImage,TOutputImage>
::CallCopyOutputRegionToInputRegion(QulumeRegionType &destRegion, const VolumeRegionType &srcRegion)
{
  ExtractImageFilterRegionCopierType extractImageRegionCopier;
  extractImageRegionCopier(destRegion, srcRegion, m_ExtractionRegion);
}

template< class TInputImage, class TOutputImage >
void CalculateQuantificationParameters< TInputImage, TOutputImage >
::GenerateInputRequestedRegion()
throw( InvalidRequestedRegionError )
{
  // Call the superclass' implementation of this method. This should
  // copy the output requested region to the input requested region.
  Superclass::GenerateInputRequestedRegion();

  // This filter needs all of the input
  typename CalculateQuantificationParameters<TInputImage,TOutputImage>::QulumePointerType image =
    const_cast< QulumeType * >( this->GetInput(0) );
  if ( image )
    {
    image->SetRequestedRegion( this->GetInput(0)->GetLargestPossibleRegion() );
    }  
}

template< class TInputImage, class TOutputImage >
void CalculateQuantificationParameters< TInputImage, TOutputImage >::SetInputQulume(const TInputImage* qulume)
{
  SetNthInput(0, const_cast<TInputImage*>(qulume));
}
 
template< class TInputImage, class TOutputImage >
void CalculateQuantificationParameters< TInputImage, TOutputImage >::SetInputVolume(const TOutputImage* volume)
{
  SetNthInput(1, const_cast<TOutputImage*>(volume));
}
 
template< class TInputImage, class TOutputImage >
typename TInputImage::ConstPointer CalculateQuantificationParameters< TInputImage, TOutputImage >::GetInputQulume()
{
  return static_cast< const TInputImage * >
         ( this->ProcessObject::GetInput(0) );
}
 
template< class TInputImage, class TOutputImage >
typename TOutputImage::ConstPointer CalculateQuantificationParameters< TInputImage, TOutputImage >::GetInputVolume()
{
  return static_cast< const TOutputImage * >
         ( this->ProcessObject::GetInput(1) );
}

template <class TInputImage, class TOutputImage>
void CalculateQuantificationParameters<TInputImage,TOutputImage>
::GenerateOutputInformation()
{
  // do not call the superclass' implementation of this method since
  // this filter allows the input and the output to be of different dimensions
 
  // get pointers to the input and output
  std::cout <<std::endl<< "generate output information" << std::endl;
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

  // Set the output image size to the same value as the extraction region.
  outputPtr->SetLargestPossibleRegion( outputImageRegion );
}


template <class TInputImage, class TOutputImage>
void CalculateQuantificationParameters<TInputImage,TOutputImage>::BeforeThreadedGenerateData()
{
//	PkSolver::pk_clear();
	m_inputQulume = this->GetInputQulume();
	m_inputVolume = this->GetInputVolume();
		
	//To test inputs
	/*typedef itk::ImageFileWriter<QulumeType>								testInputQulumeWriterType;   
	testInputQulumeWriterType::Pointer testInputQulumeWriter = testInputQulumeWriterType::New();
	testInputQulumeWriter->SetInput(inputQulume);
	testInputQulumeWriter->SetFileName("D:/Codes/Slicer4/Modules/CLI/PkModeling/Data/testInputQulume.nrrd");
	testInputQulumeWriter->Update();

	typedef itk::ImageFileWriter<VolumeType>								testInputVolumeWriterType;   
	testInputVolumeWriterType::Pointer testInputVolumeWriter = testInputVolumeWriterType::New();
	testInputVolumeWriter->SetInput(inputVolume);
	testInputVolumeWriter->SetFileName("D:/Codes/Slicer4/Modules/CLI/PkModeling/Data/testInputVolume.nrrd");
	testInputVolumeWriter->Update();*/
		
	//transfer input qulume to vector image
	m_inputVectorVolume = this->QulumeToVectorVolume(const_cast<QulumeType *>(static_cast<const QulumeType * >(m_inputQulume)));    	
    
    VectorVolumeSizeType vectorVolumeSize = m_inputVectorVolume->GetLargestPossibleRegion().GetSize();
    std::cerr<<vectorVolumeSize<<std::endl;
	//calculate parameters
	m_ktransVolume = VolumeType::New();
	m_ktransVolume->SetRegions(m_inputVolume->GetLargestPossibleRegion());
	m_ktransVolume->Allocate(); 
	m_ktransVolume->FillBuffer(0);
	m_veVolume = VolumeType::New();
	m_veVolume->SetRegions(m_inputVolume->GetLargestPossibleRegion());
	m_veVolume->Allocate();  
	m_veVolume->FillBuffer(0);
	m_maxSlopeVolume = VolumeType::New();
	m_maxSlopeVolume->SetRegions(m_inputVolume->GetLargestPossibleRegion());
	m_maxSlopeVolume->Allocate();  
	m_maxSlopeVolume->FillBuffer(0);
	m_aucVolume = VolumeType::New();
	m_aucVolume->SetRegions(m_inputVolume->GetLargestPossibleRegion());
	m_aucVolume->Allocate();  
	m_aucVolume->FillBuffer(0);
        
    m_timeSize = (int)m_inputVectorVolume->GetNumberOfComponentsPerPixel();
    
	m_TimeMinute = new float[m_inputQulume->GetLargestPossibleRegion().GetSize()[3]]();
	for(int i =0; i<(m_inputQulume->GetLargestPossibleRegion().GetSize()[3]);i++)
	{
		m_TimeMinute[i] = m_timeAxis[i]/60;
	//	std::cerr<<m_TimeMinute[i]<<std::endl;
	}

	m_averageAIFCon = new float[m_timeSize]();
    int aif_BATIndex = 0;
    int aif_FirstPeakIndex = 0;
    float aif_MaxSlope = 0.0f;
	CalculateAverageAIF(m_inputVectorVolume, m_inputVolume, m_averageAIFCon);    
    PkSolver::compute_bolus_arrival_time (m_timeSize, m_averageAIFCon, aif_BATIndex, aif_FirstPeakIndex, aif_MaxSlope);
   /* for(int i=0;i<m_timeSize;i++)
    {
        std::cerr<<"averageAIFCon"<<i<<"and:"<<m_averageAIFCon[i]<<"\n"<<std::endl;
    }*/
	//std::cerr<<"aif_BATIndex:"<<aif_BATIndex<<std::endl;
	//std::cerr<<"m_AUCTimeInterval:"<<m_AUCTimeInterval<<std::endl;
	m_aifAUC = PkSolver::area_under_curve(m_timeSize, m_timeAxis, m_averageAIFCon, aif_BATIndex, m_AUCTimeInterval);
   // std::cerr<<"m_aifAUC:"<<m_aifAUC<<std::endl;
    /* for(int i=0;i<m_timeSize;i++)
    {
        std::cerr<<"averageAIFCon"<<i<<"and:"<<m_averageAIFCon[i]<<"\n"<<std::endl;
    }*/
    std::cerr<<"beforeThreaded done!"<<std::endl;
	// 	

}

template <class TInputImage, class TOutputImage>
void CalculateQuantificationParameters<TInputImage,TOutputImage>
#if ITK_VERSION_MAJOR < 4
::ThreadedGenerateData( const typename Superclass::OutputImageRegionType & outputRegionForThread, int itkNotUsed(threadId) )
#else
::ThreadedGenerateData( const typename Superclass::OutputImageRegionType& outputRegionForThread, ThreadIdType itkNotUsed(threadId))
#endif
{		
	//get input qulume and mask volume	
	VectorVoxelType vectorVoxel;				
    	
	float tempFpv = 0.0f;
	float tempKtrans = 0.0f;
	float tempVe = 0.0f;
	float tempMaxSlope = 0.0f;
	float tempAUC = 0.0f;
    int BATIndex = 0;
    int FirstPeakIndex = 0;        
    
	VectorVolumeIterType inputVectorVolumeIter(m_inputVectorVolume, outputRegionForThread);	
	VolumeIterType ktransVolumeIter(m_ktransVolume, outputRegionForThread);			
	VolumeIterType veVolumeIter(m_veVolume, outputRegionForThread);			
	VolumeIterType maxSlopeVolumeIter(m_maxSlopeVolume, outputRegionForThread);			
	VolumeIterType aucVolumeIter(m_aucVolume, outputRegionForThread);			
	
	//set up optimizer
    itk::LevenbergMarquardtOptimizer::Pointer  optimizer = itk::LevenbergMarquardtOptimizer::New(); ///...    
    PkSolver::LMCostFunction::Pointer costFunction = PkSolver::LMCostFunction::New(); ///...    
    //PkSolver::LMCostFunction::ParametersType initialValue(PkSolver::LMCostFunction::SpaceDimension); ///...
    //initialValue[0] = 0.1;     //Ktrans //...
    //initialValue[1] = 0.5;     //ve //...
    //initialValue[2] = 0.1;     //f_pv //...
    
    //std::cerr <<"costFunction\n"<<std::endl;
    costFunction->SetNumberOfValues (m_timeSize); //...
   // std::cerr <<"m_timeSize\n"<<std::endl;
    costFunction->set_hematocrit (m_hematocrit); //...
    //std::cerr <<"m_hematocrit\n"<<std::endl;
    
  //  std::cerr <<"useGradient\n"<<std::endl;
    optimizer->UseCostFunctionGradientOff(); //...
    optimizer->SetUseCostFunctionGradient(0); //...
    
    std::cerr <<"LevenbergMarquardtOptimizer"<<std::endl;
    
    PkSolver::CommandIterationUpdateLevenbergMarquardt::Pointer observer = PkSolver::CommandIterationUpdateLevenbergMarquardt::New();//...
    optimizer->AddObserver( itk::IterationEvent(), observer ); //...
    optimizer->AddObserver( itk::FunctionEvaluationIterationEvent(), observer );//...    
    // end set up optimizer
    
    //float * tempConcentration = new float[m_timeSize]();
    
  //  std::cout << "m_timeSize:"<<m_timeSize << std::endl;
	//typename VectorVolumeType::IndexType vectorVolumeIndex;		

	while (!inputVectorVolumeIter.IsAtEnd())
    {		
			
        vectorVoxel = inputVectorVolumeIter.Get();	
        //tempConcentration = const_cast<float *>(vectorVoxel.GetDataPointer());
        //for(int i = 0;i<m_timeSize;i++)
        //{
          //  std::cerr<<tempConcentration[i]<<"\n"<<std::endl;
        //}
		/*for (int i=100; i<110;i++)
		{
			 std::cerr<<"averageAIFCon"<<i<<"and:"<<m_averageAIFCon[i]<<"\n"<<std::endl;
			 std::cout<<"i,aifcon:"<<i<<","<<m_averageAIFCon[i]<<std::endl;
		}*/
        /*PkSolver::pk_solver(m_timeSize, m_TimeMinute,  const_cast<float *>(vectorVoxel.GetDataPointer()), 
								m_averageAIFCon, tempKtrans, tempVe, tempFpv,
		  					    m_fTol, m_gTol, m_xTol, m_epsilon, m_maxIter, m_hematocrit);*/
       PkSolver::pk_solver_boost(m_timeSize, m_TimeMinute, 
                                 const_cast<float *>(vectorVoxel.GetDataPointer()),m_averageAIFCon, 
                                 tempKtrans, tempVe, tempFpv, 
                                 m_fTol,m_gTol,m_xTol,
                                 m_epsilon,m_maxIter,
                                 optimizer,costFunction);
      /*PkSolver::pk_solver_opt(m_timeSize, m_timeAxis, 
                              const_cast<float *>(vectorVoxel.GetDataPointer()),m_averageAIFCon, 
                              tempKtrans, tempVe, tempFpv,
							  m_fTol,m_gTol,m_xTol,
							  m_epsilon,m_maxIter,
                              PkSolverOpt);
*/
		PkSolver::compute_bolus_arrival_time (m_timeSize, const_cast<float *>(vectorVoxel.GetDataPointer()), BATIndex, FirstPeakIndex, tempMaxSlope);
		tempAUC = (PkSolver::area_under_curve(m_timeSize, m_timeAxis, const_cast<float *>(vectorVoxel.GetDataPointer()), BATIndex, m_AUCTimeInterval))/m_aifAUC;
		
		//for debug
		//vectorVolumeIndex = inputVectorVolumeIter.GetIndex();	
		//if((vectorVolumeIndex[0]==225)&&(vectorVolumeIndex[1]==253)&&(vectorVolumeIndex[2]==9))
		////if(tempAUC<0)
		//{
		//	
		//	
		//	std::stringstream iterStringx;
		//	std::stringstream iterStringy;
		//	std::stringstream iterStringz;

		//	int voxelx = (int)vectorVolumeIndex[0];
		//	int voxely = (int)vectorVolumeIndex[1];
		//	int voxelz = (int)vectorVolumeIndex[2];
		//	iterStringx << voxelx;			
		//	iterStringy << voxely;
		//	iterStringz << voxelz;

		//	//std::string voxelInfo = "D:/Codes/Slicer4/Modules/CLI/PkModeling/Data/DCEMRIData/voxel_AUCLessThan0_" +iterStringx.str()+std::string("_")+iterStringy.str()+"_"+iterStringz.str()+std::string(".txt");

		//	FILE* fp = fopen(voxelInfo.c_str(),"w");
		//	float tempAUCNumber = PkSolver::area_under_curve(m_timeSize, m_timeAxis, const_cast<float *>(vectorVoxel.GetDataPointer()), BATIndex, m_AUCTimeInterval);
		//	float * tempConcentration = const_cast<float *>(vectorVoxel.GetDataPointer());
		//	fprintf(fp,"%f\n",(float)BATIndex);
		//	std::cout<<"BATIndex:"<<BATIndex<<std::endl;
		//	for(int itemp = 0;itemp<m_timeSize;itemp++)
		//	{
		//		fprintf(fp,"%f\n",vectorVoxel.GetDataPointer()[itemp]);
		//		std::cout<<"i, concentration:"<<itemp<<","<<tempConcentration[itemp]<<std::endl;
		//	}
		//	fprintf(fp,"%f\n",tempAUCNumber);
		//	//std::cout<<"tempAUCValue:"<<tempAUCNumber<<std::endl;
		//	fclose(fp);
		//}


		ktransVolumeIter.Set(static_cast<VolumePixelType>(tempKtrans));       
		veVolumeIter.Set(static_cast<VolumePixelType>(tempVe));       
		maxSlopeVolumeIter.Set(static_cast<VolumePixelType>(tempMaxSlope));
     	aucVolumeIter.Set(static_cast<VolumePixelType>(tempAUC));
              
		++ktransVolumeIter;
		++veVolumeIter;
		++maxSlopeVolumeIter;
		++aucVolumeIter;
		++inputVectorVolumeIter;	
	}
	
	
	// To test outputs
    /*typedef itk::ImageFileReader<VolumeType>								testReaderType;   
	testReaderType::Pointer ktransWriter = testReaderType::New();	
	ktransWriter->SetFileName("D:/Codes/Slicer4/Modules/CLI/PkModeling/Data/maskVolume.nrrd");
	ktransWriter->Update();*/
	//delete [] tempConcentration;	
    //std::cerr << "delete averageAIF" << std::endl;
}

template <class TInputImage, class TOutputImage>
void CalculateQuantificationParameters<TInputImage,TOutputImage>::AfterThreadedGenerateData()
{
    std::cerr << "prepare for output" << std::endl;
    this->SetNthOutput(0,m_ktransVolume);
	this->SetNthOutput(1,m_veVolume);
	this->SetNthOutput(2,m_maxSlopeVolume);
	this->SetNthOutput(3,m_aucVolume);		
    
    delete [] m_timeAxis;
	delete [] m_TimeMinute;
    delete [] m_averageAIFCon;
	std::cerr << "Ready for output!" << std::endl;	
//	PkSolver::pk_report();
}
    
template <class TInputImage, class TOutputImage>
void
CalculateQuantificationParameters<TInputImage, TOutputImage>::CalculateAverageAIF(VectorVolumePointerType inputVectorVolume, VolumeConstPointerType inputVolume,float*& averageAIF)
{	
	VectorVolumeIterType inputVectorVolumeIter(inputVectorVolume, inputVectorVolume->GetRequestedRegion());
	VolumeConstIterType inputVolumeIter(inputVolume, inputVolume->GetRequestedRegion());
	
	inputVectorVolumeIter.GoToBegin();
	inputVolumeIter.GoToBegin();
	VectorVoxelType vectorVoxel;
	VolumeIndexType volumeIndex;
	int numberVoxels = 0;
	//float * averageAIF = new float[(int)inputVectorVolume->GetNumberOfComponentsPerPixel()]();
		
	while (!inputVectorVolumeIter.IsAtEnd())
    {
		if (inputVolumeIter.Get()!=0) //Pixel value >0 will be consider to be a landmark pixel. 
		{						
			volumeIndex = inputVolumeIter.GetIndex();
			//std::cout << "volumeIndex:"<< volumeIndex[0]<<","<<volumeIndex[1]<<","<<volumeIndex[2]<< std::endl;	
			
			std::stringstream iterStringx;
			std::stringstream iterStringy;
			std::stringstream iterStringz;

			int voxelx = (int)volumeIndex[0];
			int voxely = (int)volumeIndex[1];
			int voxelz = (int)volumeIndex[2];
			iterStringx << voxelx;			
			iterStringy << voxely;
			iterStringz << voxelz;

			std::string voxelInfo = "D:/Codes/Slicer4/Modules/CLI/PkModeling/Data/DCEMRIData/voxel_" +iterStringx.str()+std::string("_")+iterStringy.str()+"_"+iterStringz.str()+std::string("_SI.txt");
			FILE* fp = fopen(voxelInfo.c_str(),"w");

			numberVoxels +=1;
			//std::cout << "numberVoxels:"<< numberVoxels << std::endl;	
			vectorVoxel = inputVectorVolumeIter.Get();					
			//std::cout << "getdataPointer"<< vectorVoxel.GetDataPointer()[9] << std::endl;	
			for(int i = 0; i < (int)inputVectorVolume->GetNumberOfComponentsPerPixel(); i++) 
			{
				//std::cout <<vectorVoxel.GetDataPointer()[i]<< std::endl;
				fprintf(fp,"%f\n",vectorVoxel.GetDataPointer()[i]);
				averageAIF[i]+=(float)vectorVoxel.GetDataPointer()[i];		//not an average value	
				//std::cout << "averageAIF[i]"<< averageAIF[i] << std::endl;	
			}
			fclose(fp);
		}	
		++inputVolumeIter;
		++inputVectorVolumeIter;
	}
	
	std::string voxelAveInfo = "D:/Codes/Slicer4/Modules/CLI/PkModeling/Data/DCEMRIData/voxel_average_SI.txt";
	FILE* fpa = fopen(voxelAveInfo.c_str(),"w");
	for(int i = 0; i < (int)inputVectorVolume->GetNumberOfComponentsPerPixel(); i++) 
	{		
		averageAIF[i]=averageAIF[i]/numberVoxels;			
		fprintf(fpa,"%f\n",averageAIF[i]);
		//std::cout << "averageAIF"<< averageAIF[i] << std::endl;	
	}
	fclose(fpa);
	//std::cout << "Calculated average AIF" << std::endl;
	//return averageAIF;
}

template <class TInputImage, class TOutputImage>
typename CalculateQuantificationParameters<TInputImage, TOutputImage>::VectorVolumePointerType 
CalculateQuantificationParameters<TInputImage, TOutputImage>::QulumeToVectorVolume(QulumePointerType inputQulume)
{		
	std::cout << std::endl << "Qulume To VectorVolume" << std::endl;
	typename ImageToVectorImageFilterType::Pointer imageToVectorImageFilter = ImageToVectorImageFilterType::New();
	
	VolumePointerType volumeTemp;
    VectorVolumePointerType outputVectorVolume;
	
	QulumeRegionType inputQulumeRegion = inputQulume->GetLargestPossibleRegion();
	QulumeSizeType inputQulumeSize = inputQulumeRegion.GetSize();
	typename QulumeType::IndexType extractStartIndex;
	QulumeSizeType extractSize;
	QulumeRegionType extractRegion;
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

template <class TInputImage, class TOutputImage>
void CalculateQuantificationParameters<TInputImage,TOutputImage>
::SetTimeAxis(vcl_vector<float> inputTimeAxis)
{
	//std::cerr << "in set time axis" << std::endl;
	m_timeAxis = new float[inputTimeAxis.size()]();	
		
	for(int i=0;i<(int)inputTimeAxis.size();++i)
	{		
		m_timeAxis[i] = static_cast<float>(inputTimeAxis[i]);
		//std::cerr << std::endl << "timeAxis[tempIndex]:"<< m_timeAxis[i] << std::endl;	
	}
}

template <class TInputImage, class TOutputImage>
float* CalculateQuantificationParameters<TInputImage,TOutputImage>
::GetTimeAxis()
{
	return m_timeAxis;
}


template <class TInputImage, class TOutputImage>
void CalculateQuantificationParameters<TInputImage,TOutputImage>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "m_fTol: "									<< m_fTol					                << std::endl;
  os << indent << "m_gTol: "									<< m_gTol									<< std::endl;
  os << indent << "m_xTol: "									<< m_xTol									<< std::endl;
  os << indent << "m_epsilon: "									<< m_epsilon			                    << std::endl;
  os << indent << "m_maxIter: "									<< m_maxIter	                            << std::endl;
  os << indent << "m_hematocrit: "								<< m_hematocrit	                            << std::endl;
}

} // end namespace itk

