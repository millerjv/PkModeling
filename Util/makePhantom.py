'''
Helper script to extract curves from a multivolume. Curves are selected based on
the label image supplied. The output is a new multivolume where each label value
is represented by a row of maxSamples patches of samples, each of size patchSize.
For now, newMvName must exist, and the output will overwrite it.

Input: multivolume, label volume with multiple labels, size of the phantom patch
Output: multivolume with the samples of the curves defined by the labels
'''

def makePhantom(mvName, labelName, numLabels, patchSize, maxSamples, newMvName):
  mv = slicer.util.getNode(mvName)
  lv = slicer.util.getNode(labelName)
  nmv = slicer.util.getNode(newMvName)

  '''
  nmv = slicer.mrmlScene.CreateNodeByClass('vtkMRMLMultiVolumeNode')
  nmv.SetReferenceCount(nmv.GetReferenceCount()-1)

  mvNode = slicer.mrmlScene.CreateNodeByClass('vtkMRMLMultiVolumeNode')
  mvNode.SetReferenceCount(mvNode.GetReferenceCount()-1)
  mvNode.SetScene(slicer.mrmlScene)
  mvNode.SetAttribute("MultiVolume.FrameLabels",mv.GetAttribute("MultiVolume.FrameLabels"))
  mvNode.SetAttribute("MultiVolume.FrameIdentifyingDICOMTagName",mv.GetAttribute("MultiVolume.FrameIdentifyingDICOMTagName"))
  mvNode.SetAttribute('MultiVolume.NumberOfFrames',mv.GetAttrubute("MultiVolume.NumberOfFrames"))
  mvNode.SetAttribute('MultiVolume.FrameIdentifyingDICOMTagUnits',mv.GetAttrubute("MultiVolume.FrameIdentifyingDICOMTagUnits"))

  mvNode.SetNumberOfFrames()
  mvNode.SetLabelName(self.multiVolumeTagsUnits[frameTag])
  mvNode.SetLabelArray(frameLabelsArray)
  '''

  image = mv.GetImageData()
  extent = image.GetExtent()
  print extent

  nFrames = image.GetNumberOfScalarComponents()
  newImage = vtk.vtkImageData()
  newExtent = (0,patchSize*maxSamples-1,0,patchSize*numLabels-1,0,0)
  newImage.SetExtent(newExtent)
  newImage.AllocateScalars(image.GetScalarType(),nFrames)
  print('New image allocated:')
  print newImage

  import numpy
  newArraySize = patchSize*patchSize*numLabels*maxSamples
  newArray = vtk.util.numpy_support.vtk_to_numpy(newImage.GetPointData().GetScalars())
  labelArray = vtk.util.numpy_support.vtk_to_numpy(lv.GetImageData().GetPointData().GetScalars())
  mvArray = vtk.util.numpy_support.vtk_to_numpy(mv.GetImageData().GetPointData().GetScalars())

  dim = (extent[1]+1,extent[3]+1,extent[5]+1)
  labelSampleId = [0]*numLabels
  print dim
  for i in range(dim[0]):
    for j in range(dim[1]):
      for k in range(dim[2]):
        idx = k*dim[0]*dim[1]+j*dim[0]+i
        label = labelArray[idx]
        if label>0 and label<numLabels+1:
          sampleId = labelSampleId[label-1]
          if sampleId<maxSamples:
            for pi in range(patchSize):
              for pj in range(patchSize):
                pixelIdx = (label-1)*patchSize*patchSize*maxSamples+patchSize*maxSamples*pj+sampleId*patchSize+pi
                if pixelIdx<newArraySize:
                  newArray[pixelIdx]=mvArray[idx]
                  print newArray[pixelIdx]
                  print mvArray[idx]
                else:
                  print('Index exceeds size!')
            labelSampleId[label-1] = labelSampleId[label-1]+1
  
  nmv.SetAndObserveImageData(newImage)
  print newImage

  return newArray
