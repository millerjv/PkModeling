/*=========================================================================
  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkConcentrationToQuantitativeImageFilter.h,v $
  Language:  C++
  Date:      $Date: 2012/05/01 14:28:51 $
  Version:   $Revision: 0.0 $
=========================================================================*/
#ifndef __itkConcentrationToQuantitativeImageFilter_h
#define __itkConcentrationToQuantitativeImageFilter_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkArray.h"
#include "itkArray2D.h"
#include "itkVectorImage.h"
#include "itkImageToVectorImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkCastImageFilter.h"
#include "PkSolver.h"

namespace itk
{
/** \class ConcentrationToQuantitativeImageFilter 
 * \brief Calculates quantitative imaging parameters from concentration curves.
 *
 * This filter computes Pk modeling quantitative images from
 * concentration curves. The input volume is a vector image of
 * concentration curves represented in floating point.  The output is
 * a series of floating point images of quantitative parameters.
 *
 * An second input, specifying the location of the arterial input
 * function, allows for the calculation to be adjusted for blood
 * verses tissue.
 *
 * \note
 * This work is part of the National Alliance for Medical Image Computing 
 * (NAMIC), funded by the National Institutes of Health through the NIH Roadmap
 * for Medical Research, Grant U54 EB005149.
 * 
 */

template <class TInputImage, class TMaskImage, class TOutputImage>
class ITK_EXPORT ConcentrationToQuantitativeImageFilter : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Convenient typedefs for simplifying declarations. */
  typedef TInputImage                             VectorVolumeType;
  typedef typename VectorVolumeType::Pointer      VectorVolumePointerType;
  typedef typename VectorVolumeType::ConstPointer VectorVolumeConstPointerType;
  typedef typename VectorVolumeType::PixelType    VectorVolumePixelType;
  typedef typename VectorVolumeType::RegionType   VectorVolumeRegionType;
  typedef typename VectorVolumeType::SizeType     VectorVolumeSizeType;
  typedef itk::ImageRegionConstIterator<VectorVolumeType> VectorVolumeConstIterType;

  typedef TMaskImage                            MaskVolumeType;
  typedef typename MaskVolumeType::Pointer      MaskVolumePointerType;
  typedef typename MaskVolumeType::ConstPointer MaskVolumeConstPointerType;
  typedef typename MaskVolumeType::PixelType    MaskVolumePixelType;
  typedef typename MaskVolumeType::RegionType   MaskVolumeRegionType;
  typedef typename MaskVolumeType::SizeType     MaskVolumeSizeType;
  typedef itk::ImageRegionConstIterator<MaskVolumeType> MaskVolumeConstIterType;

  typedef TOutputImage                                    OutputVolumeType;
  typedef typename OutputVolumeType::Pointer              OutputVolumePointerType;
  typedef typename OutputVolumeType::ConstPointer         OutputVolumeConstPointerType;
  typedef typename OutputVolumeType::PixelType            OutputVolumePixelType;
  typedef typename OutputVolumeType::RegionType           OutputVolumeRegionType;
  typedef typename OutputVolumeType::IndexType            OutputVolumeIndexType;
  typedef itk::ImageRegionIterator<OutputVolumeType>      OutputVolumeIterType;
  typedef itk::ImageRegionConstIterator<OutputVolumeType> OutputVolumeConstIterType;

  typedef itk::VariableLengthVector<float>          VectorVoxelType;

  /** Standard class typedefs. */
  typedef ConcentrationToQuantitativeImageFilter   Self;
  typedef ImageToImageFilter<VectorVolumeType,OutputVolumeType> Superclass;
  typedef SmartPointer<Self>                        Pointer;
  typedef SmartPointer<const Self>                  ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( ConcentrationToQuantitativeImageFilter, ImageToImageFilter );

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
  itkGetMacro( fTol, float);
  itkSetMacro( fTol, float);
  itkGetMacro( gTol, float);
  itkSetMacro( gTol, float);
  itkGetMacro( xTol, float);
  itkSetMacro( xTol, float);
  itkGetMacro( epsilon, float);
  itkSetMacro( epsilon, float);
  itkGetMacro( maxIter, int);
  itkSetMacro( maxIter, int);
  itkGetMacro( hematocrit, float);
  itkSetMacro( hematocrit, float);
  itkGetMacro( AUCTimeInterval, float);
  itkSetMacro( AUCTimeInterval, float);
  void SetTimeAxis(const std::vector<float>& inputTimeAxis);

  const std::vector<float>& GetTimeAxis();

  /** ImageDimension enumeration */
  itkStaticConstMacro(VectorVolumeDimension, unsigned int,
                      VectorVolumeType::ImageDimension);
  itkStaticConstMacro(MaskVolumeDimension, unsigned int,
                      MaskVolumeType::ImageDimension);
  itkStaticConstMacro(OutputVolumeDimension, unsigned int,
                      OutputVolumeType::ImageDimension);

  void SetAIFMask(const MaskVolumeType* volume);
  const TMaskImage* GetAIFMask() const;

  TOutputImage* GetKTransVolume();
  TOutputImage* GetVEVolume();
  TOutputImage* GetMaxSlopeVolume();
  TOutputImage* GetAUCVolume();

protected:
  ConcentrationToQuantitativeImageFilter();
  ~ConcentrationToQuantitativeImageFilter(){
  }
  void PrintSelf(std::ostream& os, Indent indent) const;

  void BeforeThreadedGenerateData();

#if ITK_VERSION_MAJOR < 4
  void ThreadedGenerateData( const typename Superclass::OutputImageRegionType& outputRegionForThread, int threadId );

#else
  void ThreadedGenerateData( const typename Superclass::OutputImageRegionType& outputRegionForThread,
                             ThreadIdType threadId );

#endif

  std::vector<float> CalculateAverageAIF(const VectorVolumeType* inputVectorVolume, const MaskVolumeType* maskVolume);


private:
  ConcentrationToQuantitativeImageFilter(const Self &); // purposely not implemented
  void operator=(const Self &); // purposely not implemented

  float  m_T1Pre;
  float  m_TR;
  float  m_FA;
  float  m_RGD_relaxivity;
  float  m_S0GradThresh;
  float  m_fTol;
  float  m_gTol;
  float  m_xTol;
  float  m_epsilon;
  int    m_maxIter;
  float  m_hematocrit;
  float  m_AUCTimeInterval;
  std::vector<float> m_TimeAxis;

  // variables to cache information to share between threads
  std::vector<float> m_averageAIFConcentration;
  float  m_aifAUC;
};

}; // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkConcentrationToQuantitativeImageFilter.hxx"
#endif

#endif
