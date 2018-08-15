import os
import unittest
import vtk, qt, ctk, slicer
from slicer.ScriptedLoadableModule import *
import logging
import numpy

#
# NeuroPixelGuide
#

class NeuroPixelGuide(ScriptedLoadableModule):
  """Uses ScriptedLoadableModule base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def __init__(self, parent):
    ScriptedLoadableModule.__init__(self, parent)
    self.parent.title = "NeuroPixelGuide" # TODO make this more human readable by adding spaces
    self.parent.categories = ["IGT"]
    self.parent.dependencies = []
    self.parent.contributors = ["Mostafa Rezaie (Yanik Lab, ETHZ)"]
    
    self.parent.helpText = """
Guide Neuro Pixel to the correct location (Yanik Lab)
"""
    self.parent.helpText += self.getDefaultModuleDocumentationLink()
    self.parent.acknowledgementText = """
This is developed in Yanik Lab for implanting Neuropixel probe in the brain
"""

#
# NeuroPixelGuideWidget
#

class NeuroPixelGuideWidget(ScriptedLoadableModuleWidget):
  """Uses ScriptedLoadableModuleWidget base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setup(self):
    ScriptedLoadableModuleWidget.setup(self)
    self.stage = None
    # Instantiate and connect widgets ...

    #
    # Parameters Area
    #
    parametersCollapsibleButton = ctk.ctkCollapsibleButton()
    parametersCollapsibleButton.text = "Parameters"
    self.layout.addWidget(parametersCollapsibleButton)

    # Layout within the dummy collapsible button
    parametersFormLayout = qt.QFormLayout(parametersCollapsibleButton)

    #Parameters
    self.firstLandMarkSelected = False
    self.secondLandMarkSelected = False
    self.currentVisitingLandMark = 0
    self.passCount = 0
    
    #
    # input volume selector
    #
    self.inputSelector = slicer.qMRMLNodeComboBox()
    self.inputSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
    self.inputSelector.selectNodeUponCreation = True
    self.inputSelector.addEnabled = False
    self.inputSelector.removeEnabled = False
    self.inputSelector.noneEnabled = False
    self.inputSelector.showHidden = False
    self.inputSelector.showChildNodeTypes = False
    self.inputSelector.setMRMLScene( slicer.mrmlScene )
    self.inputSelector.setToolTip( "Pick the input to the algorithm." )
    parametersFormLayout.addRow("Input Volume: ", self.inputSelector)

    #
    # output volume selector
    #
    self.outputSelector = slicer.qMRMLNodeComboBox()
    self.outputSelector.nodeTypes = ["vtkMRMLScalarVolumeNode"]
    self.outputSelector.selectNodeUponCreation = True
    self.outputSelector.addEnabled = True
    self.outputSelector.removeEnabled = True
    self.outputSelector.noneEnabled = True
    self.outputSelector.showHidden = False
    self.outputSelector.showChildNodeTypes = False
    self.outputSelector.setMRMLScene( slicer.mrmlScene )
    self.outputSelector.setToolTip( "Pick the output to the algorithm." )
    parametersFormLayout.addRow("Output Volume: ", self.outputSelector)

    #
    # threshold value
    #
    self.imageThresholdSliderWidget = ctk.ctkSliderWidget()
    self.imageThresholdSliderWidget.singleStep = 0.1
    self.imageThresholdSliderWidget.minimum = -100
    self.imageThresholdSliderWidget.maximum = 100
    self.imageThresholdSliderWidget.value = 0.5
    self.imageThresholdSliderWidget.setToolTip("Set threshold value for computing the output image. Voxels that have intensities lower than this value will set to zero.")
    parametersFormLayout.addRow("Image threshold", self.imageThresholdSliderWidget)

    #
    # check box to trigger taking screen shots for later use in tutorials
    #
    self.enableScreenshotsFlagCheckBox = qt.QCheckBox()
    self.enableScreenshotsFlagCheckBox.checked = 0
    self.enableScreenshotsFlagCheckBox.setToolTip("If checked, take screen shots for tutorials. Use Save Data to write them to disk.")
    parametersFormLayout.addRow("Enable Screenshots", self.enableScreenshotsFlagCheckBox)

    #
    # Connect to Positioner
    #
    ConnectionWarningLabel = qt.QLabel( "Connect to the Stage" )
    ConnectionWarningLabel.setWordWrap( True )
    parametersFormLayout.addRow(ConnectionWarningLabel)

    self.connectStageButton = qt.QPushButton("Connect To the Stage")
    self.connectStageButton.toolTip = "Connect to Stage"
    self.connectStageButton.enabled = True
    parametersFormLayout.addRow(self.connectStageButton)
    

    #
    # input fiducial list selector
    #
    fiducialWarningLabel = qt.QLabel( "Note: Parent transforms of fiducials are not used. Fiducials should be defined in the coordinate system that is being registered." )
    fiducialWarningLabel.setWordWrap( True )
    parametersFormLayout.addRow(fiducialWarningLabel)
        
    self.inputFiducialSelector = slicer.qMRMLNodeComboBox()
    self.inputFiducialSelector.nodeTypes = ( ("vtkMRMLMarkupsFiducialNode"), "" )
    self.inputFiducialSelector.selectNodeUponCreation = True
    self.inputFiducialSelector.addEnabled = False
    self.inputFiducialSelector.removeEnabled = False
    self.inputFiducialSelector.noneEnabled = False
    self.inputFiducialSelector.showHidden = False
    self.inputFiducialSelector.showChildNodeTypes = False
    self.inputFiducialSelector.setMRMLScene( slicer.mrmlScene )
    self.inputFiducialSelector.setToolTip( "Pick the input fiducial list for the algorithm." )
    parametersFormLayout.addRow("Input fiducials: ", self.inputFiducialSelector)

    #
    # First Landmark Button
    #
    self.set1stLandMarkButton = qt.QPushButton("Reset")
    self.set1stLandMarkButton.toolTip = "Reset"
    self.set1stLandMarkButton.enabled = True
    parametersFormLayout.addRow(self.set1stLandMarkButton)

    self.set2ndLandMarkButton = qt.QPushButton("Done")
    self.set2ndLandMarkButton.toolTip = "Set n-st Landmark"
    self.set2ndLandMarkButton.enabled = False
    parametersFormLayout.addRow(self.set2ndLandMarkButton)

    #
    # Target
    #
    self.targetFiducialSelector = slicer.qMRMLNodeComboBox()
    self.targetFiducialSelector.nodeTypes = ( ("vtkMRMLMarkupsFiducialNode"), "" )
    self.targetFiducialSelector.selectNodeUponCreation = True
    self.targetFiducialSelector.addEnabled = False
    self.targetFiducialSelector.removeEnabled = False
    self.targetFiducialSelector.noneEnabled = False
    self.targetFiducialSelector.showHidden = False
    self.targetFiducialSelector.showChildNodeTypes = False
    self.targetFiducialSelector.setMRMLScene( slicer.mrmlScene )
    self.targetFiducialSelector.setToolTip( "Pick the target fiducial list for the algorithm." )
    parametersFormLayout.addRow("Target fiducials: 1st Target, 2ed Entrance", self.targetFiducialSelector)

    #
    # Second Landmark Button
    #
    self.applyButton = qt.QPushButton("Cacluate the Target")
    self.applyButton.toolTip = "Run the algorithm."
    self.applyButton.enabled = False
    parametersFormLayout.addRow(self.applyButton)


    # connections
    self.set1stLandMarkButton.connect('clicked(bool)', self.onSet1stLandMarkButton)
    self.set2ndLandMarkButton.connect('clicked(bool)', self.onSet2ndLandMarkButton)
    self.applyButton.connect('clicked(bool)', self.onApplyButton)
    self.connectStageButton.connect('clicked(bool)', self.onConnectStageButton)
    self.inputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
    self.outputSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
    self.inputFiducialSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)
    self.targetFiducialSelector.connect("currentNodeChanged(vtkMRMLNode*)", self.onSelect)

    # Add vertical spacer
    self.layout.addStretch(1)

    # Refresh Apply button state
    self.onSelect()

  def cleanup(self):
    pass

  def onSelect(self):
      pass
      #self.applyButton.enabled = self.firstLandMarkSelected and self.secondLandMarkSelected
      #self.inputSelector.currentNode() and self.outputSelector.currentNode()

  def onConnectStageButton(self):
    self.stage = stageContol('COM1','SM7')
    if (self.stage.connect() == False):
        print("Fail to find the stage, running the emulator")
        self.connectStageButton.setText("Connection.. Failed (emulator)")
        # TODO: CONNECT TO STAGE
        # IF CONNECTED THIS BUTTON DISCONNECT FROM THE STAGE
        qt.QMessageBox.information(slicer.util.mainWindow(),
                                   'NeuroPixel Probe Guide', 'Cannot Find the Stage, Running Emulator')


  def onApplyButton(self):
    logic = NeuroPixelGuideLogic()
    enableScreenshotsFlag = self.enableScreenshotsFlagCheckBox.checked
    imageThreshold = self.imageThresholdSliderWidget.value

    entrancePoint =  [0, 0, 0]
    targetPoint =    [0, 0, 0]
    
    targetFiducials = self.targetFiducialSelector.currentNode()
    targetFiducials.GetNthFiducialPosition( 0, targetPoint )
    targetFiducials.GetNthFiducialPosition( 1, entrancePoint )
    print('target: ', targetPoint)
    print('Entrace: ', entrancePoint)
    print('Now Moving the Probe to the Target Position:')
    self.applyButton.enabled = False

    scene = slicer.mrmlScene
    line = self.lineModel(scene, entrancePoint, targetPoint, "NeuroPixelProbe", (1,1,0))
    
    endofProbPoint = [5*entrancePoint[i] - 4*targetPoint[i] for i in numpy.arange(3)]
    norm = numpy.sqrt(numpy.sum([(endofProbPoint[i] - entrancePoint[i])**2 for i in numpy.arange(3)]))
    
    endofProbPoint = [entrancePoint[i] + (endofProbPoint[i] - entrancePoint[i])/norm*10 for i in numpy.arange(3)]

    arrow = self.arrowModel(scene, entrancePoint, endofProbPoint, "NeuroPixelProbe", (1,.3,0))
  
    punchingPlaneNormal = numpy.cross(entrancePoint,targetPoint)
    punchingPlane = self.planeModel(scene, punchingPlaneNormal, entrancePoint, "PunchingPlane", (1,0,0))
    
    # TODO: MOVE TO TAGET POSITION
    #logic.run(self.inputSelector.currentNode(), self.outputSelector.currentNode(), imageThreshold, enableScreenshotsFlag)

  def onSet1stLandMarkButton(self):
      print('Reset')
      if self.stage is None:
        qt.QMessageBox.information(slicer.util.mainWindow(),
                                     'NeuroPixel Probe Guide', 'Stage is not connected')
        return
      inputFiducials = self.inputFiducialSelector.currentNode()
      firstPoint = [0, 0, 0]
      inputFiducials.GetNthFiducialPosition( 0, firstPoint )
      self.firstLandMarkSelected = True
      self.set2ndLandMarkButton.setText("Set the Landmark #"+str(0)+" ["+inputFiducials.GetNthFiducialLabel(0)+"] pass# 0")
      self.set2ndLandMarkButton.enabled = True
      self.applyButton.enabled = False
      self.currentVisitingLandMark = 0
      self.passCount = 0;
      print("Reset the probe registeration..")
      self.landmarkMartix = None
      self.stage.goApproach(1)

  def onSet2ndLandMarkButton(self):
    print('Landmark')
    inputFiducials = self.inputFiducialSelector.currentNode()
    secondPoint = [0, 0, 0]
    if self.landmarkMartix is None:
        self.landmarkMartix = numpy.zeros((inputFiducials.GetNumberOfFiducials(),3))
        #find the det
        for i in numpy.arange(inputFiducials.GetNumberOfFiducials()):
            inputFiducials.GetNthFiducialPosition(i,secondPoint)
            self.landmarkMartix[i][0] = secondPoint[0]
            self.landmarkMartix[i][1] = secondPoint[1]
            self.landmarkMartix[i][2] = secondPoint[2]
        
        self.landmarkMartix = numpy.dot(self.landmarkMartix.transpose(),self.landmarkMartix)
        print("DET = ", numpy.abs((numpy.linalg.det(self.landmarkMartix))))
        #sol = np.linalg.lstsq(A, b)[0]

        if (numpy.abs((numpy.linalg.det(self.landmarkMartix)))<1):
            qt.QMessageBox.information(slicer.util.mainWindow(),
                                       'NeuroPixel Probe Guide', 'Landmarks are not independent enough, You need at least 3 independet landmark for registeration')
            self.landmarkMartix = None
            return
    
    if (self.currentVisitingLandMark < inputFiducials.GetNumberOfFiducials()):
        inputFiducials.GetNthFiducialPosition( self.currentVisitingLandMark, secondPoint )
        self.secondLandMarkSelected = True
        self.currentVisitingLandMark += 1;
        print("\nGet the Stage Location... X, Y, Z")
        stagePosition = [0, 0, 0]
        self.stage.getPosition(stagePosition)
        print(secondPoint)
    
    if (self.currentVisitingLandMark >= inputFiducials.GetNumberOfFiducials()):
        # TODO: MOVE TO TAGET POSITION
        print("\nRetracting the Probe...")
        self.stage.goApproach(0)
        self.passCount +=1;
        self.currentVisitingLandMark = 0;
    self.set2ndLandMarkButton.setText("Set the Landmark #"+str(self.currentVisitingLandMark)+" ["+inputFiducials.GetNthFiducialLabel(self.currentVisitingLandMark)+"] pass# "+str(self.passCount))


    if (self.passCount > 1):
        self.set2ndLandMarkButton.setText("Done")
        self.set2ndLandMarkButton.enabled = False
        self.applyButton.enabled = True
                
    
    print("Landmarks: ", self.currentVisitingLandMark,inputFiducials.GetNumberOfFiducials())

  def lineModel(self, scene, point1, point2, name, color):
    
        """ Create a line to reflect the puncture path"""
        #Line mode source
        
        line = vtk.vtkLineSource()
        line.SetPoint1(point1)#(point1[0][0], point1[0][1], point1[0][2])
        line.SetPoint2(point2)#(point2[0][0], point2[0][1], point2[0][2])
        
        # Create model node
        
        lineModel = slicer.vtkMRMLModelNode()
        lineModel.SetScene(scene)
        lineModel.SetName(name)
        lineModel.SetAndObservePolyData(line.GetOutput())
        
        # Create display node
        
        lineModelDisplay = slicer.vtkMRMLModelDisplayNode()
        lineModelDisplay.SetColor(color)
        lineModelDisplay.SetScene(scene)
        scene.AddNode(lineModelDisplay)
        lineModel.SetAndObserveDisplayNodeID(lineModelDisplay.GetID())
        lineModelDisplay.SetInputPolyDataConnection(line.GetOutputPort())
        scene.AddNode(lineModel)
        return line

  def arrowModel(self, scene, point1, point2, name, color):
    
        """ Create a line to reflect the puncture path"""
        #Line mode source
        
        arrow = vtk.vtkLineSource()
        #arrow.SetShaftRadius(0.01)
        #arrow.SetTipLength(.9)
        arrow.SetPoint1(point1)#(point1[0][0], point1[0][1], point1[0][2])
        arrow.SetPoint2(point2)#(point2[0][0], point2[0][1], point2[0][2])
        
        tubes = vtk.vtkTubeFilter()
        tubes.SetInputConnection(arrow.GetOutputPort())
        tubes.SetRadius(0.8)
        tubes.SetNumberOfSides(6)
        
        # Create model node
        
        arrowModel = slicer.vtkMRMLModelNode()
        arrowModel.SetScene(scene)
        arrowModel.SetName(name)
        arrowModel.SetAndObservePolyData(tubes.GetOutput())
        
        # Create display node
        
        arrowModelDisplay = slicer.vtkMRMLModelDisplayNode()
        arrowModelDisplay.SetColor(color)
        arrowModelDisplay.SetScene(scene)
        scene.AddNode(arrowModelDisplay)
        arrowModel.SetAndObserveDisplayNodeID(arrowModelDisplay.GetID())
        arrowModelDisplay.SetInputPolyDataConnection(tubes.GetOutputPort())
        scene.AddNode(arrowModel)
        return arrow


  def planeModel(self, scene, normal, origin, name, color):
        """ Create a plane model node which can be viewed in the 3D View """
        #A plane source
        
        plane = vtk.vtkPlane()
        plane.SetOrigin(origin)
        plane.SetNormal(normal)
        
        planeSample = vtk.vtkSampleFunction()
        planeSample.SetImplicitFunction(plane)
        planeSample.SetModelBounds(-10,10,-10,10,-10,10)
        planeSample.SetSampleDimensions(10,10,10)
        planeSample.ComputeNormalsOff()
        
        planeContour = vtk.vtkContourFilter()
        #        planeContour.SetInput(planeSample.GetOutput())
        planeContour.SetInputData(planeSample.GetOutput())
        
        # Create plane model node
        planeNode = slicer.vtkMRMLModelNode()
        planeNode.SetScene(scene)
        planeNode.SetName(name)
        planeNode.SetAndObservePolyData(planeContour.GetOutput())
        
        # Create plane display model node
        planeModelDisplay = slicer.vtkMRMLModelDisplayNode()
        planeModelDisplay.SetColor(color)
        planeModelDisplay.SetBackfaceCulling(2)
        planeModelDisplay.SetScene(scene)
        scene.AddNode(planeModelDisplay)
        planeNode.SetAndObserveDisplayNodeID(planeModelDisplay.GetID())
        
        #Add to scene
        
        #        planeModelDisplay.SetInputPolyData(planeContour.GetOutput())
        planeModelDisplay.SetInputPolyDataConnection(planeContour.GetOutputPort())
        scene.AddNode(planeNode)
        return plane

#
# NeuroPixelGuideLogic
#

class NeuroPixelGuideLogic(ScriptedLoadableModuleLogic):
  """This class should implement all the actual
  computation done by your module.  The interface
  should be such that other python code can import
  this class and make use of the functionality without
  requiring an instance of the Widget.
  Uses ScriptedLoadableModuleLogic base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def hasImageData(self,volumeNode):
    """This is an example logic method that
    returns true if the passed in volume
    node has valid image data
    """
    if not volumeNode:
      logging.debug('hasImageData failed: no volume node')
      return False
    if volumeNode.GetImageData() is None:
      logging.debug('hasImageData failed: no image data in volume node')
      return False
    return True

  def isValidInputOutputData(self, inputVolumeNode, outputVolumeNode):
    """Validates if the output is not the same as input
    """
    if not inputVolumeNode:
      logging.debug('isValidInputOutputData failed: no input volume node defined')
      return False
    if not outputVolumeNode:
      logging.debug('isValidInputOutputData failed: no output volume node defined')
      return False
    if inputVolumeNode.GetID()==outputVolumeNode.GetID():
      logging.debug('isValidInputOutputData failed: input and output volume is the same. Create a new volume for output to avoid this error.')
      return False
    return True

  def takeScreenshot(self,name,description,type=-1):
    # show the message even if not taking a screen shot
    slicer.util.delayDisplay('Take screenshot: '+description+'.\nResult is available in the Annotations module.', 3000)

    lm = slicer.app.layoutManager()
    # switch on the type to get the requested window
    widget = 0
    if type == slicer.qMRMLScreenShotDialog.FullLayout:
      # full layout
      widget = lm.viewport()
    elif type == slicer.qMRMLScreenShotDialog.ThreeD:
      # just the 3D window
      widget = lm.threeDWidget(0).threeDView()
    elif type == slicer.qMRMLScreenShotDialog.Red:
      # red slice window
      widget = lm.sliceWidget("Red")
    elif type == slicer.qMRMLScreenShotDialog.Yellow:
      # yellow slice window
      widget = lm.sliceWidget("Yellow")
    elif type == slicer.qMRMLScreenShotDialog.Green:
      # green slice window
      widget = lm.sliceWidget("Green")
    else:
      # default to using the full window
      widget = slicer.util.mainWindow()
      # reset the type so that the node is set correctly
      type = slicer.qMRMLScreenShotDialog.FullLayout

    # grab and convert to vtk image data
    qimage = ctk.ctkWidgetsUtils.grabWidget(widget)
    imageData = vtk.vtkImageData()
    slicer.qMRMLUtils().qImageToVtkImageData(qimage,imageData)

    annotationLogic = slicer.modules.annotations.logic()
    annotationLogic.CreateSnapShot(name, description, type, 1, imageData)

  def run(self, inputVolume, outputVolume, imageThreshold, enableScreenshots=0):
    """
    Run the actual algorithm
    """

    if not self.isValidInputOutputData(inputVolume, outputVolume):
      slicer.util.errorDisplay('Input volume is the same as output volume. Choose a different output volume.')
      return False

    logging.info('Processing started')

    # Compute the thresholded output volume using the Threshold Scalar Volume CLI module
    cliParams = {'InputVolume': inputVolume.GetID(), 'OutputVolume': outputVolume.GetID(), 'ThresholdValue' : imageThreshold, 'ThresholdType' : 'Above'}
    cliNode = slicer.cli.run(slicer.modules.thresholdscalarvolume, None, cliParams, wait_for_completion=True)

    # Capture screenshot
    if enableScreenshots:
      self.takeScreenshot('NeuroPixelGuideTest-Start','MyScreenshot',-1)

    logging.info('Processing completed')

    return True


class NeuroPixelGuideTest(ScriptedLoadableModuleTest):
  """
  This is the test case for your scripted module.
  Uses ScriptedLoadableModuleTest base class, available at:
  https://github.com/Slicer/Slicer/blob/master/Base/Python/slicer/ScriptedLoadableModule.py
  """

  def setUp(self):
    """ Do whatever is needed to reset the state - typically a scene clear will be enough.
    """
    slicer.mrmlScene.Clear(0)

  def runTest(self):
    """Run as few or as many tests as needed here.
    """
    self.setUp()
    self.test_NeuroPixelGuide1()

  def test_NeuroPixelGuide1(self):
    """ Ideally you should have several levels of tests.  At the lowest level
    tests should exercise the functionality of the logic with different inputs
    (both valid and invalid).  At higher levels your tests should emulate the
    way the user would interact with your code and confirm that it still works
    the way you intended.
    One of the most important features of the tests is that it should alert other
    developers when their changes will have an impact on the behavior of your
    module.  For example, if a developer removes a feature that you depend on,
    your test should break so they know that the feature is needed.
    """

    self.delayDisplay("Starting the test")
    #
    # first, get some data
    #
    import urllib
    downloads = (
        ('http://slicer.kitware.com/midas3/download?items=5767', 'FA.nrrd', slicer.util.loadVolume),
        )

    for url,name,loader in downloads:
      filePath = slicer.app.temporaryPath + '/' + name
      if not os.path.exists(filePath) or os.stat(filePath).st_size == 0:
        logging.info('Requesting download %s from %s...\n' % (name, url))
        urllib.urlretrieve(url, filePath)
      if loader:
        logging.info('Loading %s...' % (name,))
        loader(filePath)
    self.delayDisplay('Finished with download and loading')

    volumeNode = slicer.util.getNode(pattern="FA")
    logic = NeuroPixelGuideLogic()
    self.assertIsNotNone( logic.hasImageData(volumeNode) )
    self.delayDisplay('Test passed!')





"""
Draw The line and line spacer
"""
class DrawLineModel:
    
    """ Create the reference planes and the puncture path """
    def __init__(self, dataset):
        
        # Markers, target and entry point
        
        markers = dataset.position[0]
        
        target = dataset.position[1][0]
        
        entry = dataset.position[2][0]
        
        
        
        # Normal of the reference plane
        
        normalRed = dataset.normalRed
        
        normalGreen = dataset.normalGreen
        
        normalBlue = dataset.normalBlue
        
        
        
        scene = slicer.mrmlScene
        
        n = scene.GetNumberOfNodesByClass('vtkMRMLModelNode')
        
        model = []
        
        modelName = []
        
        for i in range(n):
            
            model.append(scene.GetNthNodeByClass(i, 'vtkMRMLModelNode'))
            
            modelName.append(model[i].GetName())
        
        for i in range(len(model)):
            
            if modelName[i][0:5] == 'Plane' or modelName[i][0:8] == 'Puncture' or modelName[i][0:7] == 'Project' or modelName[i][0:4] == 'Line':
                
                scene.RemoveNode(model[i])
        
        
        
        #Reference plane model
        
        planeRed = self.planeModel(scene, normalRed, markers[0], 'Plane-Red',(1,0,0))
        
        planeGreen = self.planeModel(scene, normalGreen, markers[3], 'Plane-Green', (0,1,0))
        
        planeBlue = self.planeModel(scene, normalBlue, (0,0,0), 'Plane-Blue', (0,0,1))
        
        
        
        #Puncture path model
        
        puncturePath = self.lineModel(scene, target, entry, 'Puncture-Path',(0,1,1))
        
        
        
        #Compute the project points in reference planes
        
        projectPoints = numpy.zeros((6,3))
        
        planeRed.GeneralizedProjectPoint(target,projectPoints[0])
        
        planeRed.GeneralizedProjectPoint(entry, projectPoints[1])
        
        
        
        planeGreen.GeneralizedProjectPoint(target,projectPoints[2])
        
        planeGreen.GeneralizedProjectPoint(entry,projectPoints[3])
        
        
        
        planeBlue.GeneralizedProjectPoint(target,projectPoints[4])
        
        planeBlue.GeneralizedProjectPoint(entry,projectPoints[5])
        
        
        
        #Project points model
        
        planeDict = {0:'Red', 1:'Green', 2:'Blue'}
        
        for i in range(3):
            
            targetProject = self.pointModel(scene,projectPoints[i*2],'Project-Target-Point-' + planeDict[i],(1,1,0))
            
            entryProject = self.pointModel(scene, projectPoints[i*2+1], 'Project-Entry-Point-' + planeDict[i],(1,1,0))
        
        
        
        #The lines between the project points and the entry and target point
        
        for i in range(3):
            
            line1 = self.lineModel(scene, target, projectPoints[i*2], 'Line', (1,1,0))
            
            line2 = self.lineModel(scene, entry, projectPoints[i*2 + 1], 'Line', (1,1,0))
        
        
        
        for i in range(0, 6, 2):
            
            line3 = self.lineModel(scene, projectPoints[i], projectPoints[i+1], 'Line', (1,1,0))



    def planeModel(self, scene, normal, origin, name, color):
    
        """ Create a plane model node which can be viewed in the 3D View """
        #A plane source
        
        plane = vtk.vtkPlane()
        
        plane.SetOrigin(origin)
        
        plane.SetNormal(normal)
        
        
        
        planeSample = vtk.vtkSampleFunction()
        
        planeSample.SetImplicitFunction(plane)
        
        planeSample.SetModelBounds(-100,100,-100,100,-100,100)
        
        planeSample.SetSampleDimensions(100,100,100)
        
        planeSample.ComputeNormalsOff()
        
        planeContour = vtk.vtkContourFilter()
        
        #        planeContour.SetInput(planeSample.GetOutput())
        planeContour.SetInputData(planeSample.GetOutput())
        
        
        
        # Create plane model node
        
        planeNode = slicer.vtkMRMLModelNode()
        
        planeNode.SetScene(scene)
        
        planeNode.SetName(name)
        
        planeNode.SetAndObservePolyData(planeContour.GetOutput())
        
        
        
        # Create plane display model node
        
        planeModelDisplay = slicer.vtkMRMLModelDisplayNode()
        
        planeModelDisplay.SetColor(color)
        
        planeModelDisplay.SetBackfaceCulling(0)
        
        planeModelDisplay.SetScene(scene)
        
        scene.AddNode(planeModelDisplay)
        
        planeNode.SetAndObserveDisplayNodeID(planeModelDisplay.GetID())
        
        
        
        #Add to scene
        
        #        planeModelDisplay.SetInputPolyData(planeContour.GetOutput())
        planeModelDisplay.SetInputPolyDataConnection(planeContour.GetOutputPort())
        
        scene.AddNode(planeNode)
        
        return plane
    
    
    
    def lineModel(self, scene, point1, point2, name, color):
        
        """ Create a line to reflect the puncture path"""
        #Line mode source
        
        line = vtk.vtkLineSource()
        
        line.SetPoint1(point1)#(point1[0][0], point1[0][1], point1[0][2])
        
        line.SetPoint2(point2)#(point2[0][0], point2[0][1], point2[0][2])
        
        
        
        # Create model node
        
        lineModel = slicer.vtkMRMLModelNode()
        
        lineModel.SetScene(scene)
        
        lineModel.SetName(name)
        
        lineModel.SetAndObservePolyData(line.GetOutput())
        
        
        
        # Create display node
        
        lineModelDisplay = slicer.vtkMRMLModelDisplayNode()
        
        lineModelDisplay.SetColor(color)
        
        lineModelDisplay.SetScene(scene)
        
        scene.AddNode(lineModelDisplay)
        
        lineModel.SetAndObserveDisplayNodeID(lineModelDisplay.GetID())
        
        
        
        #Add to scene
        
        #        lineModelDisplay.SetInputPolyData(line.GetOutput())
        lineModelDisplay.SetInputPolyDataConnection(line.GetOutputPort())
        
        scene.AddNode(lineModel)
        
        return line



    def pointModel(self, scene, point, name, color):
    
        """ Create a point model using sphere"""
        #Point
        
        sphere = vtk.vtkSphereSource()
        
        sphere.SetCenter(point)
        
        sphere.SetRadius(2)
        
        # Create model node
        
        pointModel = slicer.vtkMRMLModelNode()
        
        pointModel.SetScene(scene)
        
        pointModel.SetName(name)
        
        pointModel.SetAndObservePolyData(sphere.GetOutput())
        
        #Create display node
        
        pointModelDisplay = slicer.vtkMRMLModelDisplayNode()
        
        pointModelDisplay.SetColor(color)
        
        pointModelDisplay.SetScene(scene)
        
        scene.AddNode(pointModelDisplay)
        
        pointModel.SetAndObserveDisplayNodeID(pointModelDisplay.GetID())
        
        #Add to scene
        
        #        pointModelDisplay.SetInputPolyData(sphere.GetOutput())
        pointModelDisplay.SetInputPolyDataConnection(sphere.GetOutputPort())
        
        scene.AddNode(pointModel)

class stageContol:
    def __init__(self, port, deviceType):
        return

    def connect(self):
       return False

    def getPosition(self, postion):
      postion[0] = 0
      postion[1] = 1
      postion[2] = 2
      return

    def goToXYZPosition(self, postion):
      return

    def goApproach(self, position):
      return

    def getApproach(self, position):
      position = 0
      return









