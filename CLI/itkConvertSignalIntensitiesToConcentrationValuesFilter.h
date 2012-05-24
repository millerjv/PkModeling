/*=========================================================================
  Program:   Insight Segmentation & Registration Toolkit
  Module:    $SignalIntensitiesToConcentrationValues: itkConvertSignalIntensitiesToConcentrationValuesFilter.h $
  Language:  C++
  Date:      $Date: 2012/03/07 $
  Version:   $Revision: 0.0 $
=========================================================================*/
#ifndef __itkConvertSignalIntensitiesToConcentrationValuesFilter_h
#define __itkConvertSignalIntensitiesToConcentrationValuesFilter_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkImageToVectorImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkS0CalculationFilter.h"
#include "itkImageFileWriter.h"
#include "itkCastImageFilter.h"
#include "PkSolver.h"

namespace itk
{
/** \class ConvertSignalIntensitiesToConcentrationValuesFilter */

template <class TInputImage, class TOutputImage>
class ConvertSignalIntensitiesToConcentrationValuesFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Convenient typedefs for simplifying declarations. */
  typedef TInputImage                             InputImageType;
  typedef typename InputImageType::Pointer        InputImagePointerType;
  typedef typename InputImageType::ConstPointer   InputImageConstPointer;
  typedef typename InputImageType::PixelType      InputPixelType;
  typedef typename InputImageType::RegionType     InputImageRegionType;
  typedef typename InputImageType::SizeType       InputSizeType;
  typedef itk::Image<InputPixelType, 3>           InputMaskType;
  typedef itk::ImageRegionIterator<InputMaskType> InputMaskIterType;

  typedef float                                  FloatPixelType;
  typedef TOutputImage                           OutputImageType;
  typedef typename OutputImageType::Pointer      OutputImagePointer;
  typedef typename OutputImageType::ConstPointer OutputImageConstPointer;
  typedef typename OutputImageType::PixelType    OutputPixelType;
  typedef typename OutputImageType::RegionType   OutputImageRegionType;

  typedef itk::Image<FloatPixelType, 4>           InternalQulumeType;
  typedef typename InternalQulumeType::Pointer    InternalQulumePointerType;
  typedef typename InternalQulumeType::RegionType InternalQulumeRegionType;
  typedef typename InternalQulumeType::SizeType   InternalQulumeSizeType;

  typedef itk::Image<FloatPixelType, 3>                InternalVolumeType;
  typedef typename InternalVolumeType::Pointer         InternalVolumePointerType;
  typedef itk::ImageRegionIterator<InternalVolumeType> InternalVolumeIterType;
  typedef typename InternalVolumeType::RegionType      InternalVolumeRegionType;
  typedef typename InternalVolumeType::SizeType        InternalVolumeSizeType;

  typedef itk::VectorImage<FloatPixelType, 3>                InternalVectorVolumeType;
  typedef typename InternalVectorVolumeType::Pointer         InternalVectorVolumePointerType;
  typedef itk::ImageRegionIterator<InternalVectorVolumeType> InternalVectorVolumeIterType;
  typedef typename InternalVectorVolumeType::RegionType      InternalVectorVolumeRegionType;
  typedef typename InternalVectorVolumeType::SizeType        InternalVectorVolumeSizeType;

  typedef itk::VariableLengthVector<float> InternalVectorVoxelType;

  typedef itk::ExtractImageFilter<InternalQulumeType, InternalVolumeType> ExtractImageFilterType;
  typedef itk::ImageToVectorImageFilter<InternalVolumeType>               ImageToVectorImageFilterType;
  typedef itk::CastImageFilter<InputImageType, InternalQulumeType >       CastFilterType;

  /** Standard class typedefs. */
  typedef ConvertSignalIntensitiesToConcentrationValuesFilter Self;
  typedef ImageToImageFilter<InputImageType, OutputImageType> Superclass;
  typedef SmartPointer<Self>                                  Pointer;
  typedef SmartPointer<const Self>                            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( ConvertSignalIntensitiesToConcentrationValuesFilter, ImageToImageFilter );

  /** Set and get the number of DWI channels. */
  itkGetMacro( T1PreBlood, float);
  itkSetMacro( T1PreBlood, float);
  itkGetMacro( T1PreTissue, float);
  itkSetMacro( T1PreTissue, float);
  itkGetMacro( TR, float);
  itkSetMacro( TR, float);
  itkGetMacro( FA, float);
  itkSetMacro( FA, float);
  itkGetMacro( RGD_relaxivity, float);
  itkSetMacro( RGD_relaxivity, float);
  itkGetMacro( S0GradThresh, float);
  itkSetMacro( S0GradThresh, float);

  void SetAIFMask(InputMaskType* aifMaskVolume)
  {
    m_AIFMask = aifMaskVolume;
  }

  typename ConvertSignalIntensitiesToConcentrationValuesFilter<TInputImage,
                                                               TOutputImage>::InternalQulumePointerType
  VectorVolumeToQulume(InternalVectorVolumePointerType inputVectorVolume);

  typename ConvertSignalIntensitiesToConcentrationValuesFilter<TInputImage,
                                                               TOutputImage>::InternalVectorVolumePointerType
  QulumeToVectorVolume(InternalQulumePointerType inputQulume);

protected:
  ConvertSignalIntensitiesToConcentrationValuesFilter();
  virtual ~ConvertSignalIntensitiesToConcentrationValuesFilter()
  {

  }

  void PrintSelf(std::ostream& os, Indent indent) const;

  void GenerateData();

private:
  ConvertSignalIntensitiesToConcentrationValuesFilter(const Self &); //
                                                                     // purposely
                                                                     // not
                                                                     // implemented
  void operator=(const Self &);                                      //
                                                                     // purposely
                                                                     // not
                                                                     // implemented

  typename InputMaskType::Pointer m_AIFMask;
  float m_T1PreBlood;
  float m_T1PreTissue;
  float m_TR;
  float m_FA;
  float m_RGD_relaxivity;
  float m_S0GradThresh;
};

}; // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkConvertSignalIntensitiesToConcentrationValuesFilter.hxx"
#endif

#endif
