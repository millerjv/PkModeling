/*=========================================================================

  Program:   Registration stand-alone
  Module:    $HeadURL: http://svn.slicer.org/Slicer4/trunk/Modules/CLI/RigidRegistration/RigidRegistration.cxx $
  Language:  C++
  Date:      $Date: 2011-12-06 15:49:19 -0500 (Tue, 06 Dec 2011) $
  Version:   $Revision: 18864 $

  Copyright (c) Brigham and Women's Hospital (BWH) All Rights Reserved.

  See License.txt or http://www.slicer.org/copyright/copyright.txt for details.
=========================================================================*/
#include "itkImageToVectorImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkImageRegionConstIterator.h"
#include "vcl_vector.h"
#include "itkCastImageFilter.h"

#include "PkModelingCLP.h"
//#include "../SignalIntensitiesToConcentrationValues/PkSolver.h"
#include "itkOrientImageFilter.h"
#include "itkImageFileWriter.h"
#include "itkTransformFileReader.h"
#include "itkTransformFileWriter.h"

#include "itkQuaternionRigidTransform.h"
#include "itkResampleImageFilter.h"
#include "itkBinomialBlurImageFilter.h"

#include "itkPluginUtilities.h"

#include "itkTimeProbesCollectorBase.h"
#include "itkConvertSignalIntensitiesToConcentrationValuesFilter.h"
#include "itkCalculateQuantificationParametersFilter.h"

#define TESTMODE_ERROR_TOLERANCE 0.1

// Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
//
namespace
{

template <class T1, class T2>
int DoIt2( int argc, char * argv[], const T1 &, const T2 &)
{
  //
  // Command line processing
  //	
	PARSE_ARGS;

	const   unsigned int QulumeDimension = 4;
	typedef T1															QulumePixelType; 
	typedef itk::Image<QulumePixelType, QulumeDimension>				QulumeType; 
	typedef typename QulumeType::RegionType										QulumeRegionType;
	typedef itk::ImageFileReader<QulumeType>								QulumeReaderType;   
	
	const   unsigned int VolumeDimension = 3;
	typedef T2														VolumePixelType; 
	typedef itk::Image<VolumePixelType, VolumeDimension>			VolumeType; 
	typedef itk::ImageFileReader<VolumeType>						VolumeReaderType;
	
	typedef itk::Image<float,QulumeDimension>							OutputQulumeType;
	typedef itk::Image<float,VolumeDimension>							OutputVolumeType;

	typedef itk::ImageFileWriter< OutputVolumeType>						  VolumeWriterType;	
	//typedef itk::ImageFileWriter< QulumeType>						  QulumeWriterType;	
	
	//Read Qulume
	typename QulumeType::Pointer inputQulume = QulumeType::New();
	typename QulumeReaderType::Pointer qulumeReader = QulumeReaderType::New();
	qulumeReader->SetFileName(InputFourDNrrdFile.c_str());	
	qulumeReader->Update();
	inputQulume = qulumeReader->GetOutput();
	
	//Read mask
	typename VolumeType::Pointer maskVolume = VolumeType::New();
	typename VolumeReaderType::Pointer maskVolumeReader = VolumeReaderType::New();	
	maskVolumeReader->SetFileName(InputMaskFile.c_str());	
	maskVolumeReader->Update();
	maskVolume = maskVolumeReader->GetOutput();
	
	/*if(ISInputMask!=0)
	{
		std::cerr<<"Read mask!"<<std::endl;
		typename MaskReaderType::Pointer maskReader = MaskReaderType::New();
		maskReader->SetFileName(InputMaskNrrdFile.c_str());	
		maskReader->Update();	
		intensityToConcentrationConverter->SetAIFMask(maskReader->GetOutput());
		intensityToConcentrationConverter->SetT1PreBlood(T1PreBloodValue);
	}*/

	//Convert to concentration
	typedef itk::ConvertSignalIntensitiesToConcentrationValuesFilter<QulumeType,OutputQulumeType> ConvertFilterType;
	ConvertFilterType::Pointer converter = ConvertFilterType::New();	
	converter->SetAIFMask(maskVolume);
	converter->SetT1PreBlood(T1PreBloodValue);	
	converter->SetInput(inputQulume);		
	converter->SetT1PreTissue(T1PreTissueValue);
	converter->SetTR(TRValue);
	converter->SetFA(FAValue);
	converter->SetRGD_relaxivity(RelaxivityValue);
	converter->SetS0GradThresh(S0GradValue);	
	converter->Update();

	//Calculate parameters	
	typedef itk::CastImageFilter<VolumeType, OutputVolumeType> MaskCastFilterType;
	MaskCastFilterType::Pointer maskCastFilter = MaskCastFilterType::New();
	maskCastFilter->SetInput(maskVolume);
	maskCastFilter->Update();
    typedef itk::CalculateQuantificationParametersFilter<OutputQulumeType, OutputVolumeType>		QuantifierType;
	typename QuantifierType::Pointer quantifier = QuantifierType::New();		
	quantifier->SetInputQulume(const_cast<OutputQulumeType *>(converter->GetOutput()));	

	quantifier->SetInputVolume(maskCastFilter->GetOutput());	
	quantifier->SetAUCTimeInterval(AUCTimeInterval);	
	quantifier->SetTimeAxis(TimeAxis);
	quantifier->SetfTol(FTolerance);
	quantifier->SetgTol(GTolerance);
	quantifier->SetxTol(XTolerance);
	quantifier->Setepsilon(Epsilon);
	quantifier->SetmaxIter(MaxIter);
	quantifier->Sethematocrit(Hematocrit);	
	quantifier->Update();
	
	//For test
	/*std::cerr << std::endl << "in DoIt2: 3" << std::endl;
	QulumeType::RegionType qulumeRegion = inputQulume->GetLargestPossibleRegion();
    QulumeType::SizeType size = qulumeRegion.GetSize();
	printf("Input Qulume Size: %d %d %d %d. \n", (int)size[0], (int)size[1], (int)size[2], (int)size[3]);
	*/

	//get output	
	typename VolumeWriterType::Pointer ktranswriter = VolumeWriterType::New();
	ktranswriter->SetInput(const_cast<OutputVolumeType *>(quantifier->GetOutput()));
	ktranswriter->SetFileName(OutputKtransFile.c_str());
	ktranswriter->Update();

	typename VolumeWriterType::Pointer vewriter = VolumeWriterType::New();
	vewriter->SetInput(quantifier->GetOutput(1));
	vewriter->SetFileName(OutputVeFile.c_str());
	vewriter->Update();

	typename VolumeWriterType::Pointer maxSlopewriter = VolumeWriterType::New();
	maxSlopewriter->SetInput(quantifier->GetOutput(2));
	maxSlopewriter->SetFileName(OutputMaxSlopeFile.c_str());
	maxSlopewriter->Update();

	typename VolumeWriterType::Pointer aucwriter = VolumeWriterType::New();
	aucwriter->SetInput(quantifier->GetOutput(3));
	aucwriter->SetFileName(OutputAUCFile.c_str());
	aucwriter->Update();
	return EXIT_SUCCESS;
 } // end of anonymous namespace

}


int main( int argc, char * argv[] )
{

  PARSE_ARGS;

  // this line is here to be able to see the full output on the dashboard even
  // when the test succeeds (to see the reproducibility error measure)
  std::cout << "ctest needs: CTEST_FULL_OUTPUT" << std::endl;

  itk::ImageIOBase::IOPixelType     pixelType;
  itk::ImageIOBase::IOComponentType componentType;

  try
    {
    
		itk::GetImageType(InputFourDNrrdFile, pixelType, componentType);
		//std::cout << std::endl << "in try" << std::endl;
    // This filter handles all types

    switch( componentType )
      {
      case itk::ImageIOBase::CHAR:
      case itk::ImageIOBase::UCHAR:
      case itk::ImageIOBase::SHORT:
        return DoIt2( argc, argv, static_cast<short>(0),static_cast<short>(0) );
        break;
      case itk::ImageIOBase::USHORT:
      case itk::ImageIOBase::INT:
        return DoIt2( argc, argv, static_cast<int>(0),static_cast<int>(0) );
        break;
      case itk::ImageIOBase::UINT:
      case itk::ImageIOBase::ULONG:
        return DoIt2( argc, argv, static_cast<unsigned long>(0),static_cast<unsigned long>(0) );
        break;
      case itk::ImageIOBase::LONG:
        return DoIt2( argc, argv, static_cast<long>(0),static_cast<long>(0) );
        break;
      case itk::ImageIOBase::FLOAT:
        return DoIt2( argc, argv, static_cast<float>(0),static_cast<float>(0) );
        break;
      case itk::ImageIOBase::DOUBLE:
        return DoIt2( argc, argv, static_cast<float>(0),static_cast<float>(0) );
        break;
      case itk::ImageIOBase::UNKNOWNCOMPONENTTYPE:
      default:
        std::cout << "unknown component type" << std::endl;
        break;
      }
    }
  catch( itk::ExceptionObject & excep )
    {
    std::cerr << argv[0] << ": exception caught !" << std::endl;
    std::cerr << excep << std::endl;
    return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}