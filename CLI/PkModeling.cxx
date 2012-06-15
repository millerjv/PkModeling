/*=========================================================================

  Program:   PkModeling module
  Module:    $HeadURL: http://svn.slicer.org/Slicer4/trunk/Modules/CLI/PkModeling/PkModeling.cxx $
  Language:  C++
  Date:      $Date: 2012-06-06 $
  
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
  typedef T1                                           QulumePixelType;
  typedef itk::Image<QulumePixelType, QulumeDimension> QulumeType;
  typedef typename QulumeType::RegionType              QulumeRegionType;
  typedef itk::ImageFileReader<QulumeType>             QulumeReaderType;

  const   unsigned int VolumeDimension = 3;
  typedef T2                                           VolumePixelType;
  typedef itk::Image<VolumePixelType, VolumeDimension> VolumeType;
  typedef itk::ImageFileReader<VolumeType>             VolumeReaderType;

  typedef itk::Image<float,QulumeDimension> OutputQulumeType;
  typedef itk::Image<float,VolumeDimension> OutputVolumeType;

  typedef itk::ImageFileWriter< OutputVolumeType> VolumeWriterType;
  
  //Read Qulume
  typename QulumeType::Pointer inputQulume = QulumeType::New();
  typename QulumeReaderType::Pointer qulumeReader = QulumeReaderType::New();
  qulumeReader->SetFileName(InputFourDImageFileName.c_str() );
  qulumeReader->Update();
  inputQulume = qulumeReader->GetOutput();

  //Read mask
  typename VolumeType::Pointer maskVolume = VolumeType::New();
  typename VolumeReaderType::Pointer maskVolumeReader = VolumeReaderType::New();
  maskVolumeReader->SetFileName(AIFMaskFileName.c_str() );
  maskVolumeReader->Update();
  maskVolume = maskVolumeReader->GetOutput();
 
  //Convert to concentration values
  typedef itk::ConvertSignalIntensitiesToConcentrationValuesFilter<QulumeType,OutputQulumeType> ConvertFilterType;
  typename ConvertFilterType::Pointer converter = ConvertFilterType::New();
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
  typename MaskCastFilterType::Pointer maskCastFilter = MaskCastFilterType::New();
  maskCastFilter->SetInput(maskVolume);
  maskCastFilter->Update();
  typedef itk::CalculateQuantificationParametersFilter<OutputQulumeType, OutputVolumeType> QuantifierType;
  typename QuantifierType::Pointer quantifier = QuantifierType::New();
  quantifier->SetInputQulume(const_cast<OutputQulumeType *>(converter->GetOutput() ) );

  quantifier->SetInputVolume(maskCastFilter->GetOutput() );
  quantifier->SetAUCTimeInterval(AUCTimeInterval);
  quantifier->SetTimeAxis(TriggerTimes);
  quantifier->SetfTol(FTolerance);
  quantifier->SetgTol(GTolerance);
  quantifier->SetxTol(XTolerance);
  quantifier->Setepsilon(Epsilon);
  quantifier->SetmaxIter(MaxIter);
  quantifier->Sethematocrit(Hematocrit);
  quantifier->Update();

  //set output
  typename VolumeWriterType::Pointer ktranswriter = VolumeWriterType::New();
  ktranswriter->SetInput(const_cast<OutputVolumeType *>(quantifier->GetOutput() ) );
  ktranswriter->SetFileName(OutputKtransFileName.c_str() );
  ktranswriter->Update();

  typename VolumeWriterType::Pointer vewriter = VolumeWriterType::New();
  vewriter->SetInput(quantifier->GetOutput(1) );
  vewriter->SetFileName(OutputVeFileName.c_str() );
  vewriter->Update();

  typename VolumeWriterType::Pointer maxSlopewriter = VolumeWriterType::New();
  maxSlopewriter->SetInput(quantifier->GetOutput(2) );
  maxSlopewriter->SetFileName(OutputMaxSlopeFileName.c_str() );
  maxSlopewriter->Update();

  typename VolumeWriterType::Pointer aucwriter = VolumeWriterType::New();
  aucwriter->SetInput(quantifier->GetOutput(3) );
  aucwriter->SetFileName(OutputAUCFileName.c_str() );
  aucwriter->Update();
  return EXIT_SUCCESS;
}  

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

    itk::GetImageType(InputFourDImageFileName, pixelType, componentType);
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

