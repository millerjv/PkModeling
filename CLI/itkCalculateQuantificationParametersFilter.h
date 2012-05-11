/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $RCSfile: itkOtsuStatistics.h,v $
  Language:  C++
  Date:      $Date: 2008/02/7 14:28:51 $
  Version:   $Revision: 0.0 $
=========================================================================*/
#ifndef __itkCalculateQuantificationParameters_h
#define __itkCalculateQuantificationParameters_h

#include "itkImageToImageFilter.h"
#include "itkImage.h"
#include "itkArray.h"
#include "itkArray2D.h"
#include "itkVectorImage.h"
#include "itkImageToVectorImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkImageRegionIterator.h"
#include "itkCastImageFilter.h"
#include "PkSolver.h"

namespace itk
{
/** \class CalculateQuantificationParameters */

template <class TInputImage, class TOutputImage>
class ITK_EXPORT CalculateQuantificationParameters : public ImageToImageFilter<TInputImage, TOutputImage>
{
public:
  /** Convenient typedefs for simplifying declarations. */
  typedef TInputImage						 QulumeType;
  typedef typename QulumeType::Pointer       QulumePointerType;
  typedef typename QulumeType::ConstPointer  QulumeConstPointerType;
  typedef typename QulumeType::PixelType     QulumePixelType;
  typedef typename QulumeType::RegionType    QulumeRegionType;
  typedef typename QulumeType::SizeType      QulumeSizeType;

  typedef TOutputImage										VolumeType;
  typedef typename VolumeType::Pointer						VolumePointerType;
  typedef typename VolumeType::ConstPointer					VolumeConstPointerType;
  typedef typename VolumeType::PixelType					VolumePixelType;  
  typedef typename VolumeType::RegionType					VolumeRegionType;
  typedef typename VolumeType::IndexType					VolumeIndexType;
  typedef itk::ImageRegionIterator<VolumeType>				VolumeIterType; 
  typedef itk::ImageRegionConstIterator<VolumeType>				VolumeConstIterType; 

  typedef itk::VectorImage<QulumePixelType, 3>							VectorVolumeType;
  typedef typename VectorVolumeType::Pointer							VectorVolumePointerType;  
  typedef itk::ImageRegionIterator<VectorVolumeType>					VectorVolumeIterType;  
  typedef typename VectorVolumeType::RegionType							VectorVolumeRegionType;  
  typedef typename VectorVolumeType::SizeType							VectorVolumeSizeType;
    
  typedef itk::VariableLengthVector<float>						VectorVoxelType;  
  typedef itk::ImageToVectorImageFilter<VolumeType>					ImageToVectorImageFilterType;
  typedef itk::ExtractImageFilter<QulumeType, VolumeType> ExtractImageFilterType;	
  
  
  /** Standard class typedefs. */
  typedef CalculateQuantificationParameters                      Self;
  typedef ImageToImageFilter<QulumeType,VolumeType>		  Superclass;
  typedef SmartPointer<Self>                                  Pointer;
  typedef SmartPointer<const Self>                            ConstPointer;

  /** Method for creation through the object factory. */
  itkNewMacro(Self);

  /** Run-time type information (and related methods). */
  itkTypeMacro( CalculateQuantificationParameters, ImageToImageFilter );
  
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
  void SetTimeAxis(vcl_vector<float> inputTimeAxis);
  float* GetTimeAxis();

  /** ImageDimension enumeration */
  itkStaticConstMacro(QulumeDimension, unsigned int,
                      QulumeType::ImageDimension);
  itkStaticConstMacro(VolumeDimension, unsigned int,
                      VolumeType::ImageDimension);

  typedef ImageToImageFilterDetail::ExtractImageFilterRegionCopier<itkGetStaticConstMacro(QulumeDimension), 
    itkGetStaticConstMacro(VolumeDimension)> ExtractImageFilterRegionCopierType;
    
  itkGetConstMacro(ExtractionRegion, QulumeRegionType);

  void SetInputQulume(const QulumeType* qulume);
  void SetInputVolume(const VolumeType* volume);

#ifdef ITK_USE_CONCEPT_CHECKING
  /** Begin concept checking */
  itkConceptMacro(InputCovertibleToOutputCheck,
    (Concept::Convertible<QulumePixelType, VolumePixelType>));
  /** End concept checking */
#endif
  //////////////////////////
  
protected:
  CalculateQuantificationParameters();
  ~CalculateQuantificationParameters(){};
  void PrintSelf(std::ostream& os, Indent indent) const;
  void BeforeThreadedGenerateData();
  #if ITK_VERSION_MAJOR < 4
  void ThreadedGenerateData( const typename Superclass::OutputImageRegionType& outputRegionForThread, int threadId );
  #else
  void ThreadedGenerateData( const typename Superclass::OutputImageRegionType& outputRegionForThread, ThreadIdType threadId );
  #endif  
  void AfterThreadedGenerateData();
//    void GenerateData();
  ////////////////////
  virtual void GenerateInputRequestedRegion()
  throw( InvalidRequestedRegionError );  
  virtual void GenerateOutputInformation();
  virtual void CallCopyOutputRegionToInputRegion(QulumeRegionType &destRegion, const VolumeRegionType &srcRegion);
  QulumeRegionType  m_ExtractionRegion;
  VolumeRegionType  m_OutputImageRegion;
  ////////////////////  
  //  
  typename TInputImage::ConstPointer GetInputQulume(); 
  typename TOutputImage::ConstPointer GetInputVolume();
  
  void CalculateAverageAIF(VectorVolumePointerType inputVectorVolume, VolumeConstPointerType inputVolume,float*& averageAIF); 
  typename CalculateQuantificationParameters<TInputImage, TOutputImage>::VectorVolumePointerType QulumeToVectorVolume(QulumePointerType inputQulume);
private:
  CalculateQuantificationParameters(const Self &);   // purposely not implemented
  void operator=(const Self &);  // purposely not implemented

  float m_T1Pre;
  float m_TR;
  float m_FA;
  float m_RGD_relaxivity;
  float m_S0GradThresh;
  float m_fTol;
  float m_gTol;
  float m_xTol;
  float m_epsilon;
  int   m_maxIter;  
  float m_hematocrit;
  float m_AUCTimeInterval;
  float* m_timeAxis;
  float* m_TimeMinute;
    int m_timeSize;
    float* m_averageAIFCon;
    float m_aifAUC;
    //int m_BATIndex;
	//int m_FirstPeakIndex;	
    //float m_tempMaxSlope;
  typename QulumeType::ConstPointer m_inputQulume;
  typename VolumeType::ConstPointer m_inputVolume;
  VectorVolumePointerType m_inputVectorVolume;
  VolumePointerType m_ktransVolume;
  VolumePointerType m_veVolume;
  VolumePointerType m_maxSlopeVolume;
  VolumePointerType m_aucVolume;
};

}; // end namespace itk

#ifndef ITK_MANUAL_INSTANTIATION
#include "itkCalculateQuantificationParameters.hxx"
#endif

#endif
