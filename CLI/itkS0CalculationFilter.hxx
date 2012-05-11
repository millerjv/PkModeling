/*=========================================================================
  Program:   Insight Segmentation & Registration Toolkit
  Module:    $SignalIntensitiesToConcentrationValues: itkS0ForTimeSeriesInQulume.hxx$
  Language:  C++
  Date:      $Date: 2012/03/7 $
  Version:   $Revision: 1.0 $
=========================================================================*/
/*=========================================================================
itk filter to calculate S0 value.
Input: Vector volume, x*y*z(t)
Output: Volume, x*y*z
Parameter: S0GradThresh, Threshold for gradient magnitude of S0
=========================================================================*/

#ifndef _itkS0ForTimeSeriesInQulume_hxx
#define _itkS0ForTimeSeriesInQulume_hxx
#include "itkS0ForTimeSeriesInQulume.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
S0ForTimeSeriesInQulume<TInputImage, TOutputImage>::S0ForTimeSeriesInQulume()
{		
	m_S0GradThresh = 15.0f;
	this->Superclass::SetNumberOfRequiredInputs(1);
	this->Superclass::SetNumberOfRequiredOutputs(1);
	this->Superclass::SetNthOutput(0, OutputImageType::New());	
}

template <class TInputImage, class TOutputImage>
void S0ForTimeSeriesInQulume<TInputImage, TOutputImage>
::BeforeThreadedGenerateData()
{
	m_inputVectorVolume = this->GetInput();
    m_S0Volume = this->GetOutput(); 
}

template <class TInputImage, class TOutputImage>
void S0ForTimeSeriesInQulume<TInputImage, TOutputImage>
#if ITK_VERSION_MAJOR < 4
::ThreadedGenerateData( const typename Superclass::OutputImageRegionType & outputRegionForThread, int itkNotUsed(threadId) )
#else
::ThreadedGenerateData( const typename Superclass::OutputImageRegionType& outputRegionForThread, ThreadIdType itkNotUsed(threadId))
#endif
{		
    //Input is vector volume, output is volume
	std::cerr << std::endl << "Calculate S0" << std::endl;	 
	
	InputImageIterType inputVectorVolumeIter(m_inputVectorVolume, outputRegionForThread);
	OutputImageIterType S0VolumeIter(m_S0Volume, outputRegionForThread);
	
	float S0Temp = 0.0f;
	InternalVectorVoxelType vectorVoxel;	
	
	while (!inputVectorVolumeIter.IsAtEnd())
    {
		vectorVoxel = inputVectorVolumeIter.Get();
		S0Temp = compute_s0_individual_curve ((int)m_inputVectorVolume->GetNumberOfComponentsPerPixel(), const_cast<float*>( vectorVoxel.GetDataPointer()), m_S0GradThresh);
		S0VolumeIter.Set(static_cast<OutputPixelType>(S0Temp));
		++S0VolumeIter;
		++inputVectorVolumeIter;
	}	

}

/** Standard "PrintSelf" method */
template <class TInputImage, class TOutput>
void S0ForTimeSeriesInQulume<TInputImage, TOutput>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent ); 
  os << indent << "S0GradThresh: "                          << m_S0GradThresh                          << std::endl;
}

} // end namespace itk

#endif