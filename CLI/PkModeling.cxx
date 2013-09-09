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
#include "itkMultiThreader.h"
#include "itkResampleImageFilter.h"
#include "itkNearestNeighborInterpolateImageFunction.h"

#include "itkPluginUtilities.h"

#include "itkSignalIntensityToConcentrationImageFilter.h"
#include "itkConcentrationToQuantitativeImageFilter.h"

#include <sstream>
#include <fstream>

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
      if (tag == "TriggerTime" || tag == "AcquisitionTime" || tag == "SeriesTime")
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


// Read an AIF from a CSV style file.  Columns are timing and concentration.
//
//
bool GetPrescribedAIF(const std::string& fileName, 
                      std::vector<float>& timing, std::vector<float>& aif)
{
  timing.clear();
  aif.clear();

  std::string line;
  std::ifstream csv;
  csv.open(fileName.c_str());
  if (csv.fail())
    {
    std::cout << "Cannot open file " << fileName << std::endl;
    return false;
    }

  while (!csv.eof())
    {
    getline(csv, line);
      
    if (line[0] == '#')
      {
      continue;
      }
      
    std::vector<std::string> svalues;
    splitString(line, ",", svalues);  /// from PkModelingCLP.h

    if (svalues.size() < 2)
      {
      // not enough values on the line
      continue;
      }

    // only keep the time and concentration value
    std::stringstream tstream;
    float time=-1.0, value=-1.0;

    tstream << svalues[0];
    tstream >> time;
    if (tstream.fail())
      {
      // not a float, probably the column labels, skip the row
      continue;
      }
    tstream.str("");
    tstream.clear(); // need to clear the flags since at eof of the stream

    tstream << svalues[1];
    tstream >> value;
    if (tstream.fail())
      {
      // not a float, could be column labels, skip the row
      continue;
      }
    tstream.str("");
    tstream.clear(); // need to clear the flags since at eof of the stream

    timing.push_back(time);
    aif.push_back(value);
    }

  csv.close();

  if (timing.size() > 0)
    {
    return true;
    }

  return false;
}

template <class T1, class T2>
int DoIt( int argc, char * argv[], const T1 &, const T2 &)
{
  //
  // Command line processing
  //
  PARSE_ARGS;

  // itk::MultiThreader::Pointer threader = itk::MultiThreader::New();
  // threader->SetGlobalMaximumNumberOfThreads(1);

  const   unsigned int VectorVolumeDimension = 3;
  typedef T1                                                 VectorVolumePixelType;
  typedef itk::VectorImage<VectorVolumePixelType, VectorVolumeDimension> VectorVolumeType;
  typedef itk::VectorImage<float, VectorVolumeDimension>     FloatVectorVolumeType;
  typedef typename VectorVolumeType::RegionType              VectorVolumeRegionType;
  typedef itk::ImageFileReader<VectorVolumeType>             VectorVolumeReaderType;
  typedef itk::ImageFileWriter<FloatVectorVolumeType>        VectorVolumeWriterType;

  const   unsigned int MaskVolumeDimension = 3;
  typedef T2                                                   MaskVolumePixelType;
  typedef itk::Image<MaskVolumePixelType, MaskVolumeDimension> MaskVolumeType;
  typedef itk::ImageFileReader<MaskVolumeType>                 MaskVolumeReaderType;

  typedef itk::Image<float,VectorVolumeDimension> OutputVolumeType;
  typedef itk::ImageFileWriter< OutputVolumeType> OutputVolumeWriterType;

  typedef itk::ResampleImageFilter<MaskVolumeType,MaskVolumeType> ResamplerType;
  typedef itk::NearestNeighborInterpolateImageFunction<MaskVolumeType> InterpolatorType;
  
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

  //Read AIF mask
  typename MaskVolumeReaderType::Pointer aifMaskVolumeReader = MaskVolumeReaderType::New();
  typename MaskVolumeType::Pointer aifMaskVolume = 0;
  if (AIFMaskFileName != "")
    {
    aifMaskVolumeReader->SetFileName(AIFMaskFileName.c_str() );
    aifMaskVolumeReader->Update();
    aifMaskVolume = aifMaskVolumeReader->GetOutput();
    }

  //Read ROI mask
  typename MaskVolumeReaderType::Pointer roiMaskVolumeReader = MaskVolumeReaderType::New();
  typename MaskVolumeType::Pointer roiMaskVolume = 0;
  if (ROIMaskFileName != "")
    {
    roiMaskVolumeReader->SetFileName(ROIMaskFileName.c_str() );
    roiMaskVolumeReader->Update();
    roiMaskVolume = roiMaskVolumeReader->GetOutput();

    typename ResamplerType::Pointer resampler = ResamplerType::New();
    typename InterpolatorType::Pointer interpolator = InterpolatorType::New();

    resampler->SetOutputDirection(inputVectorVolume->GetDirection());
    resampler->SetOutputSpacing(inputVectorVolume->GetSpacing());
    resampler->SetOutputStartIndex(inputVectorVolume->GetBufferedRegion().GetIndex());
    resampler->SetSize(inputVectorVolume->GetBufferedRegion().GetSize());
    resampler->SetOutputOrigin(inputVectorVolume->GetOrigin());
    resampler->SetInput(roiMaskVolume);
    resampler->SetInterpolator(interpolator);
    resampler->Update();

    roiMaskVolume = resampler->GetOutput();
    }

  //Read prescribed aif
  bool usingPrescribedAIF = false;
  std::vector<float> prescribedAIFTiming;
  std::vector<float> prescribedAIF;
  if (PrescribedAIFFileName != "")
    {
    usingPrescribedAIF = GetPrescribedAIF(PrescribedAIFFileName, prescribedAIFTiming, prescribedAIF);
    }
 
  if (AIFMaskFileName == "" && !usingPrescribedAIF)
    {
    std::cerr << "Either a mask localizing the region over which to calculate the arterial input function or a prescribed arterial input function must be specified." << std::endl;
    return EXIT_FAILURE;
    }

  //Convert to concentration values
  typedef itk::SignalIntensityToConcentrationImageFilter<VectorVolumeType,MaskVolumeType,FloatVectorVolumeType> ConvertFilterType;
  typename ConvertFilterType::Pointer converter = ConvertFilterType::New();
  converter->SetInput(inputVectorVolume);
  if (!usingPrescribedAIF)
    {
    converter->SetAIFMask(aifMaskVolume);
    }

  if(ROIMaskFileName != "")
    {
    converter->SetROIMask(roiMaskVolume);
    }

  converter->SetT1PreBlood(T1PreBloodValue);
  converter->SetT1PreTissue(T1PreTissueValue);
  converter->SetTR(TRValue);
  converter->SetFA(FAValue);
  converter->SetRGD_relaxivity(RelaxivityValue);
  converter->SetS0GradThresh(S0GradValue);
  itk::PluginFilterWatcher watchConverter(converter, "Concentrations",  CLPProcessInformation,  1.0 / 20.0, 0.0);
  converter->Update();

  if(OutputConcentrationsImageFileName != "")
    {
    // need to initialize the attributes, otherwise Slicer treats 
    //  this as a Vector volume, not MultiVolume
    FloatVectorVolumeType::Pointer concentrationsVolume = converter->GetOutput();
    concentrationsVolume->SetMetaDataDictionary(inputVectorVolume->GetMetaDataDictionary());

    typename VectorVolumeWriterType::Pointer multiVolumeWriter
      = VectorVolumeWriterType::New();
    multiVolumeWriter->SetFileName(OutputConcentrationsImageFileName.c_str());
    multiVolumeWriter->SetInput(concentrationsVolume);
    multiVolumeWriter->SetUseCompression(1);
    multiVolumeWriter->Update();
    }

  //Calculate parameters
  typedef itk::ConcentrationToQuantitativeImageFilter<FloatVectorVolumeType, MaskVolumeType, OutputVolumeType> QuantifierType;
  typename QuantifierType::Pointer quantifier = QuantifierType::New();
  quantifier->SetInput(converter->GetOutput());
  if (usingPrescribedAIF)
    {
    quantifier->SetPrescribedAIF(prescribedAIFTiming, prescribedAIF);
    quantifier->UsePrescribedAIFOn();
    }
  else
    {
    quantifier->SetAIFMask(aifMaskVolume );
    }

  quantifier->SetAUCTimeInterval(AUCTimeInterval);
  quantifier->SetTiming(Timing);
  quantifier->SetfTol(FTolerance);
  quantifier->SetgTol(GTolerance);
  quantifier->SetxTol(XTolerance);
  quantifier->Setepsilon(Epsilon);
  quantifier->SetmaxIter(MaxIter);
  quantifier->Sethematocrit(Hematocrit);
  if(ROIMaskFileName != "")
    {
    quantifier->SetROIMask(roiMaskVolume);
    }

  if(ComputeFpv)
    {
    quantifier->SetModelType(itk::LMCostFunction::TOFTS_3_PARAMETER);
    }
  else
    {
    quantifier->SetModelType(itk::LMCostFunction::TOFTS_2_PARAMETER);
    }
  quantifier->SetMaskByRSquared(OutputRSquaredFileName.empty());

  itk::PluginFilterWatcher watchQuantifier(quantifier, "Quantifying",  CLPProcessInformation,  19.0 / 20.0, 1.0 / 20.0);
  quantifier->Update();


  //set output
  if (!OutputKtransFileName.empty())
    {
    typename OutputVolumeWriterType::Pointer ktranswriter = OutputVolumeWriterType::New();
    ktranswriter->SetInput(quantifier->GetKTransOutput() );
    ktranswriter->SetFileName(OutputKtransFileName.c_str() );
    ktranswriter->SetUseCompression(1);
    ktranswriter->Update();
    }

  if (!OutputVeFileName.empty())
    {
    typename OutputVolumeWriterType::Pointer vewriter = OutputVolumeWriterType::New();
    vewriter->SetInput(quantifier->GetVEOutput() );
    vewriter->SetFileName(OutputVeFileName.c_str() );
    vewriter->SetUseCompression(1);
    vewriter->Update();
    }

  if(ComputeFpv)
    {
    if (!OutputFpvFileName.empty())
      {
      typename OutputVolumeWriterType::Pointer fpvwriter =OutputVolumeWriterType::New();
      fpvwriter->SetInput(quantifier->GetFPVOutput() );
      fpvwriter->SetFileName(OutputFpvFileName.c_str() );
      fpvwriter->SetUseCompression(1);
      fpvwriter->Update();
      }
    }

  if (!OutputMaxSlopeFileName.empty())
    {
    typename OutputVolumeWriterType::Pointer maxSlopewriter = OutputVolumeWriterType::New();
    maxSlopewriter->SetInput(quantifier->GetMaxSlopeOutput() );
    maxSlopewriter->SetFileName(OutputMaxSlopeFileName.c_str() );
    maxSlopewriter->SetUseCompression(1);
    maxSlopewriter->Update();
    }

  if (!OutputAUCFileName.empty())
    {
    typename OutputVolumeWriterType::Pointer aucwriter = OutputVolumeWriterType::New();
    aucwriter->SetInput(quantifier->GetAUCOutput() );
    aucwriter->SetFileName(OutputAUCFileName.c_str() );
    aucwriter->SetUseCompression(1);
    aucwriter->Update();
    }

  if (!OutputRSquaredFileName.empty())
    {
    typename OutputVolumeWriterType::Pointer rsqwriter =OutputVolumeWriterType::New();
    rsqwriter->SetInput(quantifier->GetRSquaredOutput() );
    rsqwriter->SetFileName(OutputRSquaredFileName.c_str() );
    rsqwriter->SetUseCompression(1);
    rsqwriter->Update();
    }

  if (!OutputFittedDataImageFileName.empty())
    {
    // need to initialize the attributes, otherwise Slicer treats 
    //  this as a Vector volume, not MultiVolume
    FloatVectorVolumeType::Pointer fittedVolume = quantifier->GetFittedDataOutput();
    fittedVolume->SetMetaDataDictionary(inputVectorVolume->GetMetaDataDictionary());

    typename VectorVolumeWriterType::Pointer multiVolumeWriter
      = VectorVolumeWriterType::New();
    multiVolumeWriter->SetFileName(OutputFittedDataImageFileName.c_str());
    multiVolumeWriter->SetInput(fittedVolume);
    multiVolumeWriter->SetUseCompression(1);
    multiVolumeWriter->Update();
    }

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

