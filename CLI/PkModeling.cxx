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

#define SimpleAttributeGetMethodMacro(name, key, type)     \
type Get##name(itk::MetaDataDictionary& dictionary)           \
{\
  type value = type(); \
  if (dictionary.HasKey(key))\
    {\
    /* attributes stored as strings */ \
    std::string valueString; \
    ExposeMetaData(dictionary, key, valueString); \
    std::stringstream valueStream(valueString); \
    valueStream >> value; \
    }\
  else\
    {\
    itkGenericExceptionMacro("Missing attribute '" key "'.");\
    }\
  return value;\
}

//SimpleAttributeGetMethodMacro(EchoTime, "MultiVolume.DICOM.EchoTime", float);
SimpleAttributeGetMethodMacro(RepetitionTime, "MultiVolume.DICOM.RepetitionTime",float);
SimpleAttributeGetMethodMacro(FlipAngle, "MultiVolume.DICOM.FlipAngle", float);

std::vector<float> GetTriggerTimes(itk::MetaDataDictionary& dictionary)
{
  std::vector<float> triggerTimes;

  if (dictionary.HasKey("MultiVolume.FrameIdentifyingDICOMTagName"))
    {
    std::string tag;
    ExposeMetaData(dictionary, "MultiVolume.FrameIdentifyingDICOMTagName", tag);
    if (dictionary.HasKey("MultiVolume.FrameLabels"))
      {
      // Acquisition parameters stored as text
      std::string frameLabelsString;
      ExposeMetaData(dictionary, "MultiVolume.FrameLabels", frameLabelsString);
      std::stringstream frameLabelsStream(frameLabelsString);
      if (tag == "Trigger Time")
        {
        float t;
        while (frameLabelsStream >> t)
          {
          t /= 1000.0; // convert to seconds (are times in milliseconds?)
          triggerTimes.push_back(t);
          }
        }
      // what other frame identification methods are there?
      }
    else
      {
      itkGenericExceptionMacro("Missing attribute 'MultiVolume.FrameLabels'.")
      }
    }
  else
    {
    itkGenericExceptionMacro("Missing attribute 'MultiVolume.FrameIdentifyingDICOMTagName'.");
    }
  
  return triggerTimes;
}



template <class T1, class T2>
int DoIt2( int argc, char * argv[], const T1 &, const T2 &)
{
  //
  // Command line processing
  //
  PARSE_ARGS;

  const   unsigned int MultiVolumeDimension = 4;
  typedef T1                                                MultiVolumePixelType;
  typedef itk::Image<MultiVolumePixelType, MultiVolumeDimension> MultiVolumeType;
  typedef typename MultiVolumeType::RegionType              MultiVolumeRegionType;
  typedef itk::ImageFileReader<MultiVolumeType>             MultiVolumeReaderType;

  const   unsigned int VolumeDimension = 3;
  typedef T2                                           VolumePixelType;
  typedef itk::Image<VolumePixelType, VolumeDimension> VolumeType;
  typedef itk::ImageFileReader<VolumeType>             VolumeReaderType;

  typedef itk::Image<float,MultiVolumeDimension> OutputMultiVolumeType;
  typedef itk::Image<float,VolumeDimension> OutputVolumeType;

  typedef itk::ImageFileWriter< OutputVolumeType> VolumeWriterType;
  
  //Read MultiVolume
  typename MultiVolumeReaderType::Pointer multiVolumeReader = MultiVolumeReaderType::New();
  multiVolumeReader->SetFileName(InputFourDImageFileName.c_str() );
  multiVolumeReader->Update();
  typename MultiVolumeType::Pointer inputMultiVolume = multiVolumeReader->GetOutput();

  //Look for tags representing the acquisition parameters
  //
  //

  // Trigger times
  std::vector<float> TriggerTimes;
  try
    {
    TriggerTimes = GetTriggerTimes(inputMultiVolume->GetMetaDataDictionary());
    }
  catch (itk::ExceptionObject &exc)
    {
    itkGenericExceptionMacro(<< exc.GetDescription() 
            << " Image " << InputFourDImageFileName.c_str() 
            << " does not contain sufficient attributes to support algorithms.");
    return EXIT_FAILURE;
    }

  // // EchoTime
  // float echoTime = 0.0;
  // try 
  //   {
  //   echoTime = GetEchoTime(inputMultiVolume->GetMetaDataDictionary());
  //   }
  // catch (itk::ExceptionObject &exc)
  //   {
  //   itkGenericExceptionMacro(<< exc.GetDescription() 
  //           << " Image " << InputFourDImageFileName.c_str() 
  //           << " does not contain sufficient attributes to support algorithms.");
  //   return EXIT_FAILURE;
    
  //   }

  // FlipAngle
  float FAValue = 0.0;
  try 
    {
    FAValue = GetFlipAngle(inputMultiVolume->GetMetaDataDictionary());
    }
  catch (itk::ExceptionObject &exc)
    {
    itkGenericExceptionMacro(<< exc.GetDescription() 
            << " Image " << InputFourDImageFileName.c_str() 
            << " does not contain sufficient attributes to support algorithms.");
    return EXIT_FAILURE;
    
    }

  // RepetitionTime
  float TRValue = 0.0;
  try 
    {
    TRValue = GetRepetitionTime(inputMultiVolume->GetMetaDataDictionary());
    }
  catch (itk::ExceptionObject &exc)
    {
    itkGenericExceptionMacro(<< exc.GetDescription() 
            << " Image " << InputFourDImageFileName.c_str() 
            << " does not contain sufficient attributes to support algorithms.");
    return EXIT_FAILURE;
    
    }

  //Read mask
  typename VolumeReaderType::Pointer maskVolumeReader = VolumeReaderType::New();
  maskVolumeReader->SetFileName(AIFMaskFileName.c_str() );
  maskVolumeReader->Update();
  typename VolumeType::Pointer maskVolume = maskVolumeReader->GetOutput();
 
  //Convert to concentration values
  typedef itk::ConvertSignalIntensitiesToConcentrationValuesFilter<MultiVolumeType,OutputMultiVolumeType> ConvertFilterType;
  typename ConvertFilterType::Pointer converter = ConvertFilterType::New();
  converter->SetAIFMask(maskVolume);
  converter->SetT1PreBlood(T1PreBloodValue);
  converter->SetInput(inputMultiVolume);
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
  typedef itk::CalculateQuantificationParametersFilter<OutputMultiVolumeType, OutputVolumeType> QuantifierType;
  typename QuantifierType::Pointer quantifier = QuantifierType::New();
  quantifier->SetInputMultiVolume(const_cast<OutputMultiVolumeType *>(converter->GetOutput() ) );

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

