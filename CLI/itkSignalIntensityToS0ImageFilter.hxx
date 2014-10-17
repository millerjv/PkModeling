/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    $SignalIntensitiesToConcentrationValues: itkSignalIntensityToS0ImageFilter.hxx$
  Language:  C++
  Date:      $Date: 2012/03/7 $
  Version:   $Revision: 1.0 $

=========================================================================*/
#ifndef _itkSignalIntensityToS0ImageFilter_hxx
#define _itkSignalIntensityToS0ImageFilter_hxx

#include "itkSignalIntensityToS0ImageFilter.h"

namespace itk
{

template <class TInputImage, class TOutputImage>
SignalIntensityToS0ImageFilter<TInputImage, TOutputImage>::SignalIntensityToS0ImageFilter()
{
  m_S0GradThresh = 15.0f;

}

template <class TInputImage, class TOutputImage>
void SignalIntensityToS0ImageFilter<TInputImage, TOutputImage>
#if ITK_VERSION_MAJOR < 4
::ThreadedGenerateData( const typename Superclass::OutputImageRegionType & outputRegionForThread, int itkNotUsed(
                          threadId) )
#else
::ThreadedGenerateData( const typename Superclass::OutputImageRegionType& outputRegionForThread,
                        ThreadIdType itkNotUsed(
                          threadId) )
#endif
  {
  //Input is vector volume, output is volume

  const InputImageType* inputVectorVolume = this->GetInput();
  OutputImageType* S0Volume = this->GetOutput();

  InputImageConstIterType  inputVectorVolumeIter(inputVectorVolume, outputRegionForThread);
  OutputImageIterType S0VolumeIter(S0Volume, outputRegionForThread);

  float                   S0Temp = 0.0f;
  InternalVectorVoxelType vectorVoxel;
  InputPixelType inputVectorVoxel;

  while (!inputVectorVolumeIter.IsAtEnd() )
    {
    // copy/cast input vector to floats
    inputVectorVoxel = inputVectorVolumeIter.Get();
    vectorVoxel.SetSize(inputVectorVoxel.GetSize());
    vectorVoxel.Fill(0.0);
    vectorVoxel += inputVectorVoxel; // shorthand for a copy/cast
    S0Temp =
      compute_s0_individual_curve ( (int)inputVectorVolume->GetNumberOfComponentsPerPixel(),
                                    const_cast<float*>( vectorVoxel.GetDataPointer() ), m_S0GradThresh, m_BATCalculationMode,m_constantBAT);
    S0VolumeIter.Set(static_cast<OutputPixelType>(S0Temp) );
    ++S0VolumeIter;
    ++inputVectorVolumeIter;
    }

  }

/** Standard "PrintSelf" method */
template <class TInputImage, class TOutput>
void SignalIntensityToS0ImageFilter<TInputImage, TOutput>
::PrintSelf( std::ostream& os, Indent indent ) const
{
  Superclass::PrintSelf( os, indent );
  os << indent << "S0GradThresh: " << m_S0GradThresh << std::endl;
}

} // end namespace itk

#endif
