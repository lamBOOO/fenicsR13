from paraview.simple import *

kn_array = [0.05, 0.1, 0.2, 0.4]

for kn in kn_array:
    filename = 'thermal-flow-between-cylinders' + str(kn) + '/u_0.xdmf'
    u_0xdmf = paraview.simple.Xdmf3ReaderS(
        registrationName='u_0.xdmf', FileName=[filename])

    plotOverLine1 = paraview.simple.PlotOverLine(
        registrationName='PlotOverLine1', Input=u_0xdmf)
    plotOverLine1.Point1 = [0.0, 0.0, 0.0]
    plotOverLine1.Point2 = [2.0, 2.0, 0.0]

    spreadSheetView1 = FindViewOrCreate(
        'SpreadSheetView1', viewtype='SpreadSheetView')
    spreadSheetView1.Update()
    u_0xdmfDisplay_1 = Show(u_0xdmf, spreadSheetView1,
                            'SpreadSheetRepresentation')
    plotOverLine1Display_2 = Show(
        plotOverLine1, spreadSheetView1, 'SpreadSheetRepresentation')
    outputname = 'u_0_' + str(kn) + '.csv'
    ExportView(outputname, view=spreadSheetView1)
