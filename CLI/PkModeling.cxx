/*=========================================================================

  Program:   PkModeling module
  Module:    $HeadURL: http://svn.slicer.org/Slicer4/trunk/Modules/CLI/PkModeling/PkModeling.cxx $
  Language:  C++
  Date:      $Date: 2012-06-06 $
  
  Copyright (c) Brigham and Women's Hospital (BWH) All Rights Reserved.

  See License.txt or http://www.slicer.org/copyright/copyright.txt for details.
=========================================================================*/


#include "PkModelingCLP.h"

#include "itkMetaDataObject.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkTimeProbesCollectorBase.h"

#include "itkPluginUtilities.h"

#include "itkSignalIntensityToConcentrationImageFilter.h"
#include "itkConcentrationToQuantitativeImageFilter.h"

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
    itk::ExposeMetaData(dictionary, key, valueString);  \
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

std::vector<float> GetTiming(itk::MetaDataDictionary& dictionary)
{
  std::vector<float> triggerTimes;

  if (dictionary.HasKey("MultiVolume.FrameIdentifyingDICOMTagName"))
    {
    std::string tag;
    itk::ExposeMetaData(dictionary, "MultiVolume.FrameIdentifyingDICOMTagName", tag);
    if (dictionary.HasKey("MultiVolume.FrameLabels"))
      {
      // Acquisition parameters stored as text, FrameLabels are comma separated
      std::string frameLabelsString;
      itk::ExposeMetaData(dictionary, "MultiVolume.FrameLabels", frameLabelsString);
      std::stringstream frameLabelsStream(frameLabelsString);
      if (tag == "TriggerTime" || tag == "AcquisitionTime")
        {
        float t;
        float t0 = 0.0;
        bool first = true;
        while (frameLabelsStream >> t)
          {
          t /= 1000.0;  // convert to seconds
          if (first)
            {
            t0 = t;
            first = false;
            }
          t = t - t0;

          triggerTimes.push_back(t);
          frameLabelsStream.ignore(1); // skip the comma
          }
        }
      else
        {
        itkGenericExceptionMacro("Unrecognized frame identifying DICOM tag name " << tag);
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
int DoIt( int argc, char * argv[], const T1 &, const T2 &)
{
  //
  // Command line processing
  //
  PARSE_ARGS;

  const   unsigned int VectorVolumeDimension = 3;
  typedef T1                                                 VectorVolumePixelType;
  typedef itk::VectorImage<VectorVolumePixelType, VectorVolumeDimension> VectorVolumeType;
  typedef itk::VectorImage<float, VectorVolumeDimension>     FloatVectorVolumeType;
  typedef typename VectorVolumeType::RegionType              VectorVolumeRegionType;
  typedef itk::ImageFileReader<VectorVolumeType>             VectorVolumeReaderType;

  const   unsigned int MaskVolumeDimension = 3;
  typedef T2                                                   MaskVolumePixelType;
  typedef itk::Image<MaskVolumePixelType, MaskVolumeDimension> MaskVolumeType;
  typedef itk::ImageFileReader<MaskVolumeType>                 MaskVolumeReaderType;

  typedef itk::Image<float,VectorVolumeDimension> OutputVolumeType;
  typedef itk::ImageFileWriter< OutputVolumeType> OutputVolumeWriterType;
  
  //Read VectorVolume
  typename VectorVolumeReaderType::Pointer multiVolumeReader 
    = VectorVolumeReaderType::New();
  multiVolumeReader->SetFileName(InputFourDImageFileName.c_str() );
  multiVolumeReader->Update();
  typename VectorVolumeType::Pointer inputVectorVolume = multiVolumeReader->GetOutput();

  //Look for tags representing the acquisition parameters
  //
  //

  // Trigger times
  std::vector<float> Timing;
  try
    {
    Timing = GetTiming(inputVectorVolume->GetMetaDataDictionary());
    std::cout << "Timing: ";
    for (std::vector<float>::size_type i=0; i < Timing.size(); ++i)
      {
      std::cout << Timing[i] << ", ";
      }
    std::cout << std::endl;
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
  //   echoTime = GetEchoTime(inputVectorVolume->GetMetaDataDictionary());
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
    FAValue = GetFlipAngle(inputVectorVolume->GetMetaDataDictionary());
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
    TRValue = GetRepetitionTime(inputVectorVolume->GetMetaDataDictionary());
    }
  catch (itk::ExceptionObject &exc)
    {
    itkGenericExceptionMacro(<< exc.GetDescription() 
            << " Image " << InputFourDImageFileName.c_str() 
            << " does not contain sufficient attributes to support algorithms.");
    return EXIT_FAILURE;
    
    }

  //Read mask
  typename MaskVolumeReaderType::Pointer maskVolumeReader = MaskVolumeReaderType::New();
  maskVolumeReader->SetFileName(AIFMaskFileName.c_str() );
  maskVolumeReader->Update();
  typename MaskVolumeType::Pointer maskVolume = maskVolumeReader->GetOutput();
 
  //Convert to concentration values
  typedef itk::SignalIntensityToConcentrationImageFilter<VectorVolumeType,MaskVolumeType,FloatVectorVolumeType> ConvertFilterType;
  typename ConvertFilterType::Pointer converter = ConvertFilterType::New();
  converter->SetInput(inputVectorVolume);
  converter->SetAIFMask(maskVolume);
  converter->SetT1PreBlood(T1PreBloodValue);
  converter->SetT1PreTissue(T1PreTissueValue);
  converter->SetTR(TRValue);
  converter->SetFA(FAValue);
  converter->SetRGD_relaxivity(RelaxivityValue);
  converter->SetS0GradThresh(S0GradValue);
  itk::PluginFilterWatcher watchConverter(converter, "Concentrations",  CLPProcessInformation,  1.0 / 20.0, 0.0);
  converter->Update();

  // itk::Index<3> ind;
  // ind[0] = 30;
  // ind[1] = 100;
  // ind[2] = 0;
  // itk::ImageRegion<3> reg;
  // itk::Size<3> sz;
  // sz.Fill(1);
  // sz[1]=99;
  // reg.SetIndex(ind);
  // reg.SetSize(sz);

  //Calculate parameters
  typedef itk::ConcentrationToQuantitativeImageFilter<FloatVectorVolumeType, MaskVolumeType, OutputVolumeType> QuantifierType;
  typename QuantifierType::Pointer quantifier = QuantifierType::New();
  quantifier->SetInput(converter->GetOutput());
  quantifier->SetAIFMask(maskVolume );

  quantifier->SetAUCTimeInterval(AUCTimeInterval);
  quantifier->SetTimeAxis(Timing);
  quantifier->SetfTol(FTolerance);
  quantifier->SetgTol(GTolerance);
  quantifier->SetxTol(XTolerance);
  quantifier->Setepsilon(Epsilon);
  quantifier->SetmaxIter(MaxIter);
  quantifier->Sethematocrit(Hematocrit);
//  quantifier->GetOutput()->SetRequestedRegion(reg);
  itk::PluginFilterWatcher watchQuantifier(quantifier, "Quantifying",  CLPProcessInformation,  19.0 / 20.0, 1.0 / 20.0);
  quantifier->Update();

  // std::cout << "Ktrans: " << quantifier->GetKTransOutput()->GetPixel(ind) << std::endl;
  // std::cout << "Ve: " << quantifier->GetVEOutput()->GetPixel(ind) << std::endl;
  // std::cout << "MaxSlope: " << quantifier->GetMaxSlopeOutput()->GetPixel(ind) << std::endl;
  // std::cout << "AUC: " << quantifier->GetAUCOutput()->GetPixel(ind) << std::endl;

  //set output
  typename OutputVolumeWriterType::Pointer ktranswriter = OutputVolumeWriterType::New();
  ktranswriter->SetInput(quantifier->GetKTransOutput() );
  ktranswriter->SetFileName(OutputKtransFileName.c_str() );
  ktranswriter->Update();

  typename OutputVolumeWriterType::Pointer vewriter = OutputVolumeWriterType::New();
  vewriter->SetInput(quantifier->GetVEOutput() );
  vewriter->SetFileName(OutputVeFileName.c_str() );
  vewriter->Update();

  typename OutputVolumeWriterType::Pointer maxSlopewriter = OutputVolumeWriterType::New();
  maxSlopewriter->SetInput(quantifier->GetMaxSlopeOutput() );
  maxSlopewriter->SetFileName(OutputMaxSlopeFileName.c_str() );
  maxSlopewriter->Update();

  typename OutputVolumeWriterType::Pointer aucwriter = OutputVolumeWriterType::New();
  aucwriter->SetInput(quantifier->GetAUCOutput() );
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
        return DoIt( argc, argv, static_cast<short>(0),static_cast<short>(0) );
        break;
      case itk::ImageIOBase::USHORT:
      case itk::ImageIOBase::INT:
        return DoIt( argc, argv, static_cast<int>(0),static_cast<short>(0) );
        break;
      case itk::ImageIOBase::UINT:
      case itk::ImageIOBase::ULONG:
        return DoIt( argc, argv, static_cast<unsigned long>(0),static_cast<short>(0) );
        break;
      case itk::ImageIOBase::LONG:
        return DoIt( argc, argv, static_cast<long>(0),static_cast<short>(0) );
        break;
      case itk::ImageIOBase::FLOAT:
        return DoIt( argc, argv, static_cast<float>(0),static_cast<short>(0) );
        break;
      case itk::ImageIOBase::DOUBLE:
        return DoIt( argc, argv, static_cast<float>(0),static_cast<short>(0) );
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

