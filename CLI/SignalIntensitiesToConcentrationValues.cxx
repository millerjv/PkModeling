/*=========================================================================
  Program:   SignalIntensitiesToConcentrationValues CLI
  Module:    $HeadURL: http://svn.slicer.org/Slicer4/trunk/Modules/CLI/Converters/SignalIntensitiesToConcentrationValues.cxx $
  Language:  C++
  Date:      $Date: 2012-03-07$
  Version:   $Revision: $
=========================================================================*/
#include "SignalIntensitiesToConcentrationValuesCLP.h"
#include "itkSignalIntensityToConcentration.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkPluginUtilities.h"
#include "itkTimeProbesCollectorBase.h"

#define TESTMODE_ERROR_TOLERANCE 0.1

// Use an anonymous namespace to keep class types and function names
// from colliding when module is used as shared object module.  Every
// thing should be in an anonymous namespace except for the module
// entry point, e.g. main()
//
namespace
{

template <class T>
int DoIt( int argc, char * argv[], const T& targ)
{
  PARSE_ARGS;

  const unsigned int ImageDimension = 4;

  typedef  T                                                       PixelType;
  typedef itk::Image< PixelType, ImageDimension >                  ImageType;
  typedef itk::ImageFileReader<ImageType>                          ReaderType;
  typedef itk::ImageFileWriter<ImageType>                          ImageWriterType;
  typedef itk::SignalIntensityToConcentration<ImageType,ImageType> ConverterType;
  typename ConverterType::Pointer intensityToConcentrationConverter = ConverterType::New();

  typename ReaderType::Pointer reader = ReaderType::New();
  reader->SetFileName(InputFourDNrrdFile.c_str() );
  reader->Update();
  intensityToConcentrationConverter->SetInput(reader->GetOutput() );
  intensityToConcentrationConverter->SetT1Pre(T1PreValue);
  intensityToConcentrationConverter->SetTR(TRValue);
  intensityToConcentrationConverter->SetFA(FAValue);
  intensityToConcentrationConverter->SetRGD_relaxivity(RelaxivityValue);
  intensityToConcentrationConverter->SetS0GradThresh(S0GradValue);
  intensityToConcentrationConverter->Update();

  typename ImageWriterType::Pointer concentrationWriter = ImageWriterType::New();
  concentrationWriter->SetInput(intensityToConcentrationConverter->GetOutput() );
  concentrationWriter->SetFileName(OutputFourDNrrdFile.c_str() );
  concentrationWriter->Update();
  return EXIT_SUCCESS;
}  // end of anonymous namespace

}
int main( int argc, char * argv[] )
{

  PARSE_ARGS;

  // this line is here to be able to see the full output on the dashboard even
  // when the test succeeds (to see the reproducibility error measure)
  std::cout << std::endl << "ctest needs: CTEST_FULL_OUTPUT" << std::endl;

  itk::ImageIOBase::IOPixelType     pixelType;
  itk::ImageIOBase::IOComponentType componentType;

  try
    {

    itk::GetImageType(InputFourDNrrdFile, pixelType, componentType);
    // This filter handles all types

    switch( componentType )
      {
      case itk::ImageIOBase::CHAR:
      case itk::ImageIOBase::UCHAR:
      case itk::ImageIOBase::SHORT:
        return DoIt( argc, argv, static_cast<short>(0) );
        break;
      case itk::ImageIOBase::USHORT:
      case itk::ImageIOBase::INT:
        return DoIt( argc, argv, static_cast<int>(0) );
        break;
      case itk::ImageIOBase::UINT:
      case itk::ImageIOBase::ULONG:
        return DoIt( argc, argv, static_cast<unsigned long>(0) );
        break;
      case itk::ImageIOBase::LONG:
        return DoIt( argc, argv, static_cast<long>(0) );
        break;
      case itk::ImageIOBase::FLOAT:
        return DoIt( argc, argv, static_cast<float>(0) );
        break;
      case itk::ImageIOBase::DOUBLE:
        return DoIt( argc, argv, static_cast<float>(0) );
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
