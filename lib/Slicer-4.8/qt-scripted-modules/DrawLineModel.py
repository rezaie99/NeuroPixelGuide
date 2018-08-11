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
