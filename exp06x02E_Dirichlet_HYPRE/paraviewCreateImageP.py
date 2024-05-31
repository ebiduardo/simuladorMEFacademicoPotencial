from paraview.simple import *

# Load the state
#LoadState("exp05x02N/paraviewEstadoA_P.pvsm")
#servermanager.LoadState("exp05x02/paraviewEstadoA_P.pvsm")
expDir="exp05x02\\"
LoadState(expDir + "paraviewEstadoA_P.pvsm")

#SetActiveView(GetRenderView())
Render()

#create image
WriteImage(expDir + "potencial.png")
