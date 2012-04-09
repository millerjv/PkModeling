/*=========================================================================
  Program:   Insight Segmentation & Registration Toolkit
  Module:    $SignalIntensitiesToConcentrationValues: itkSignalIntensityToConcentration.h $
  Language:  C++
  Date:      $Date: 2012/03/07 $
  Version:   $Revision: 0.0 $
=========================================================================*/
#ifndef __itkSignalIntensityToConcentration_h
#define __itkSignalIntensityToConcentration_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkImageToVectorImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkS0ForTimeSeriesInQulume.h"
#include "PkSolver.h"

namespace itk
{
/** \class SignalIntensityToConcentration */

template <class TInputImage, class TOutputImage>
class SignalIntensityToConcentration : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Convenient typedefs for simplifying declarations. */
  typedef TInputImage                            InputImageType;
  typedef typename InputImageType::Pointer       InputImagePointerType;
  typedef typename InputImageType::ConstPointer  InputImageConstPointer;
  typedef typename InputImageType::PixelType     InputPixelType;
  typedef typename InputImageType::RegionType    InputImageRegionType;
  typedef typename InputImageType::SizeType      InputSizeType;

  typedef TOutputImage											OutputImageType;
  typedef typename OutputImageType::Pointer						OutputImagePointer;
  typedef typename OutputImageType::ConstPointer				OutputImageConstPointer;
  typedef typename OutputImageType::PixelType					OutputPixelType;  
  typedef typename OutputImageType::RegionType					OutputImageRegionType;
      
  typedef itk::Image<InputPixelType, 4>					    InternalQulumeType;
  typedef typename InternalQulumeType::Pointer			    InternalQulumePointerType;
  typedef typename InternalQulumeType::RegionType           InternalQulumeRegionType;
  typedef typename InternalQulumeType::SizeType				InternalQulumeSizeType;
    
  typedef itk::Image<InputPixelType, 3>						 InternalVolumeType;
  typedef typename InternalVolumeType::Pointer					 InternalVolumePointerType;
  typedef itk::ImageRegionIterator<InternalVolumeType>           InternalVolumeIterType;
  typedef typename InternalVolumeType::RegionType				 InternalVolumeRegionType;
  typedef typename InternalVolumeType::SizeType				 InternalVolumeSizeType;

  typedef itk::VectorImage<InputPixelType, 3>						 InternalVectorVolumeType;
  typedef typename InternalVectorVolumeType::Pointer                     InternalVectorVolumePointerType;  
  typedef itk::ImageRegionIterator<InternalVectorVolumeType>			 InternalVectorVolumeIterType;  
  typedef typename InternalVectorVolumeType::RegionType		             InternalVectorVolumeRegionType;  
  typedef typename InternalVectorVolumeType::SizeType                    InternalVectorVolumeSizeType;
  
  typedef itk::VariableLengthVector<float>                       InternalVectorVoxelType;

  typedef itk::ExtractImageFilter<InternalQulumeType, InternalVolumeType> ExtractImageFilterType;	
  typedef itk::ImageToVectorImageFilter<InternalVolumeType>				  ImageToVectorImageFilterType;
  
  /** Standard class typedefs. */
  typedef SignalIntensityToConcentration                      Self;
  typedef ImageToImageFilter<InputImageType, OutputImageType> Superclass;
  typedef SmartPointer<Self>                                  Pointer;
  typedef SmartPointer<const Self>                            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( SignalIntensityToConcentration, ImageToImageFilter );
  
  /** Set and get the number of DWI channels. */    
  itkGetMacro( T1Pre, float);
  itkSetMacro( T1Pre, float);
  itkGetMacro( TR, float);
  itkSetMacro( TR, float);
  itkGetMacro( FA, float);
  itkSetMacro( FA, float);
  itkGetMacro( RGD_relaxivity, float);
  itkSetMacro( RGD_relaxivity, float);
  itkGetMacro( S0GradThresh, float);
  itkSetMacro( S0GradThresh, float);

  typename SignalIntensityToConcentration<TInputImage, TOutputImage>::InternalQulumePointerType VectorVolumeToQulume(InternalVectorVolumePointerType inputVectorVolume);
  typename SignalIntensityToConcentration<TInputImage, TOutputImage>::InternalVectorVolumePointerType QulumeToVectorVolume(InternalQulumePointerType inputQulume);
  
protected:
  SignalIntensityToConcentration();
  virtual ~SignalIntensityToConcentration()
  {

  }
  void PrintSelf(std::ostream& os, Indent indent) const;
  
  void GenerateData();
  
  //typename SignalIntensityToConcentration<TInputImage, TOutputImage>::InternalQulumePointerType VectorVolumeToQulume(InternalVectorVolumePointerType inputVectorVolume);
  //typename SignalIntensityToConcentration<TInputImage, TOutputImage>::InternalVectorVolumePointerType QulumeToVectorVolume(InternalQulumePointerType inputQulume);
  
private:
  SignalIntensityToConcentration(const Self &);   // purposely not implemented
  void operator=(const Self &);  // purposely not implemented

  float m_T1Pre;
  float m_TR;
  float m_FA;
  float m_RGD_relaxivity;
  float m_S0GradThresh;   
};

}; // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSignalIntensityToConcentration.hxx"
#endif

#endif
