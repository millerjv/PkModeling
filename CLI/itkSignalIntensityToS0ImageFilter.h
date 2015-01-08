/*=========================================================================
  Program:   Insight Segmentation & Registration Toolkit
  Module:    $SignalIntensityToS0ImageFilter: itkSignalIntensityToS0ImageFilter.h$
  Language:  C++
  Date:      $Date: 2012/03/7 $
  Version:   $Revision: 1.0 $
=========================================================================*/
#ifndef __itkSignalIntensityToS0ImageFilter_h
#define __itkSignalIntensityToS0ImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkImageRegionIterator.h"
#include "PkSolver.h"

#include <string>

namespace itk
{
/** \class SignalIntensityToS0ImageFilter */

template <class TInputImage, class TOutputImage>
class SignalIntensityToS0ImageFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Convenient typedefs for simplifying declarations. */
  typedef TInputImage                                   InputImageType;
  typedef typename InputImageType::Pointer              InputImagePointerType;
  typedef typename InputImageType::ConstPointer         InputImageConstPointer;
  typedef typename InputImageType::PixelType            InputPixelType;
  typedef typename InputImageType::RegionType           InputImageRegionType;
  typedef typename InputImageType::SizeType             InputSizeType;
  typedef itk::ImageRegionConstIterator<InputImageType> InputImageConstIterType;

  typedef TOutputImage                              OutputImageType;
  typedef typename OutputImageType::Pointer         OutputImagePointer;
  typedef typename OutputImageType::ConstPointer    OutputImageConstPointer;
  typedef typename OutputImageType::PixelType       OutputPixelType;
  typedef typename OutputImageType::RegionType      OutputImageRegionType;
  typedef itk::ImageRegionIterator<OutputImageType> OutputImageIterType;

  typedef itk::VariableLengthVector<float> InternalVectorVoxelType;

  /** Standard class typedefs. */
  typedef SignalIntensityToS0ImageFilter                                 Self;
  typedef ImageToImageFilter<InputImageType, OutputImageType> Superclass;
  typedef SmartPointer<Self>                                  Pointer;
  typedef SmartPointer<const Self>                            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( SignalIntensityToS0ImageFilter, ImageToImageFilter );

  /** Set and get the number of DWI channels. */
  itkGetMacro( S0GradThresh, float);
  itkSetMacro( S0GradThresh, float);
  itkGetMacro( BATCalculationMode, std::string);
  itkSetMacro( BATCalculationMode, std::string);
  itkGetMacro( constantBAT, int);
  itkSetMacro( constantBAT, int);

protected:
  SignalIntensityToS0ImageFilter();
  virtual ~SignalIntensityToS0ImageFilter() {
  }
  void PrintSelf(std::ostream& os, Indent indent) const;


#if ITK_VERSION_MAJOR < 4
  void ThreadedGenerateData( const typename Superclass::OutputImageRegionType& outputRegionForThread, int threadId );

#else
  void ThreadedGenerateData( const typename Superclass::OutputImageRegionType& outputRegionForThread,
                             ThreadIdType threadId );

#endif
private:
  SignalIntensityToS0ImageFilter(const Self &); // purposely not implemented
  void operator=(const Self &);      // purposely not implemented

  float                  m_S0GradThresh;
  std::string      m_BATCalculationMode;
  int m_constantBAT;
};

}; // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkSignalIntensityToS0ImageFilter.hxx"
#endif

#endif
