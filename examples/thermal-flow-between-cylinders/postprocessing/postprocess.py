from paraview.simple import *

for kn in [0.05, 0.1, 0.2, 0.4]:
    for i in range(0,8):
        filename = '../thermal-flow-between-cylinders' + str(kn) + '/u_' + str(i) + '.xdmf'
        uxdmf = paraview.simple.Xdmf3ReaderS(
            registrationName='u.xdmf', FileName=[filename])

        plotOverLine1 = paraview.simple.PlotOverLine(
            registrationName='PlotOverLine1', Input=uxdmf)
        plotOverLine1.Point1 = [0.0, 0.0, 0.0]
        plotOverLine1.Point2 = [2.0, 2.0, 0.0]

        spreadSheetView1 = FindViewOrCreate(
            'SpreadSheetView1', viewtype='SpreadSheetView')
        spreadSheetView1.Update()
        uxdmfDisplay_1 = Show(uxdmf, spreadSheetView1,
                                'SpreadSheetRepresentation')
        plotOverLine1Display_2 = Show(
            plotOverLine1, spreadSheetView1, 'SpreadSheetRepresentation')
        outputname = 'u_' + str(kn) +'_' + str(i) + '.csv'
        ExportView(outputname, view=spreadSheetView1)
