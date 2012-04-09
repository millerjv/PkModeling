/*=========================================================================
  Program:   Insight Segmentation & Registration Toolkit
  Module:    $SignalIntensitiesToConcentrationValues: itkS0ForTimeSeriesInQulume.h$
  Language:  C++
  Date:      $Date: 2012/03/7 $
  Version:   $Revision: 1.0 $
=========================================================================*/
#ifndef __itkS0ForTimeSeriesInQulume_h
#define __itkS0ForTimeSeriesInQulume_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkVectorImage.h"
#include "itkImageToVectorImageFilter.h"
#include "itkImageRegionIterator.h"
#include "PkSolver.h"

namespace itk
{
/** \class S0ForTimeSeriesInQulume */

template <class TInputImage, class TOutputImage>
class S0ForTimeSeriesInQulume : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Convenient typedefs for simplifying declarations. */
  typedef TInputImage                            InputImageType;
  typedef typename InputImageType::Pointer       InputImagePointerType;
  typedef typename InputImageType::ConstPointer  InputImageConstPointer;
  typedef typename InputImageType::PixelType     InputPixelType;
  typedef typename InputImageType::RegionType    InputImageRegionType;
  typedef typename InputImageType::SizeType      InputSizeType;
  typedef itk::ImageRegionConstIterator<InputImageType> InputImageIterType;

  typedef TOutputImage											OutputImageType;
  typedef typename OutputImageType::Pointer						OutputImagePointer;
  typedef typename OutputImageType::ConstPointer				OutputImageConstPointer;
  typedef typename OutputImageType::PixelType					OutputPixelType;  
  typedef typename OutputImageType::RegionType					OutputImageRegionType;
  typedef itk::ImageRegionIterator<OutputImageType>				OutputImageIterType; 
    
  typedef itk::VariableLengthVector<float>                       InternalVectorVoxelType;
    
  /** Standard class typedefs. */
  typedef S0ForTimeSeriesInQulume                      Self;
  typedef ImageToImageFilter<InputImageType, OutputImageType> Superclass;
  typedef SmartPointer<Self>                                  Pointer;
  typedef SmartPointer<const Self>                            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( S0ForTimeSeriesInQulume, ImageToImageFilter );
  
  /** Set and get the number of DWI channels. */      
  itkGetMacro( S0GradThresh, float);
  itkSetMacro( S0GradThresh, float);

protected:
  S0ForTimeSeriesInQulume();
  virtual ~S0ForTimeSeriesInQulume() {}
  void PrintSelf(std::ostream& os, Indent indent) const; 
  void BeforeThreadedGenerateData();
  #if ITK_VERSION_MAJOR < 4
  void ThreadedGenerateData( const typename Superclass::OutputImageRegionType& outputRegionForThread, int threadId );
  #else
  void ThreadedGenerateData( const typename Superclass::OutputImageRegionType& outputRegionForThread, ThreadIdType threadId );
  #endif  
     
private:
  S0ForTimeSeriesInQulume(const Self &);   // purposely not implemented
  void operator=(const Self &);  // purposely not implemented
    
  float m_S0GradThresh;   
  InputImageConstPointer m_inputVectorVolume;
  OutputImagePointer m_S0Volume;
};

}; // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkS0ForTimeSeriesInQulume.hxx"
#endif

#endif
