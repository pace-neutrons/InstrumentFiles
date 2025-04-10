## PSD calibration - MERLIN
#
#  Uses the internal Mantid "tube" package to calibrate the PSD tubes
#  https://docs.mantidproject.org/nightly/concepts/calibration/PSDTubeCalibration.html
#
#  Outputs a detector.dat file (updating the previous version) and a rings map.
#
#  JRS 12/7/23


# import mantid algorithms, numpy and matplotlib
from mantid.simpleapi import *
import matplotlib.pyplot as plt
import numpy as np

from tube import calibrate
from tube_calib_fit_params import TubeCalibFitParams

#====================== calibration input  ===========================
calibration_runs = [69167]                      # calibration bin run
mask = None                                     # mask if necessary (doesn't seem to work for 1to1 spectra on MERLIN 
integration_range = [3000,10000]                # signal range in TOF
peak_height = 2000.                             # average peak height
peak_width = 9.                                # average peak width (mus)
margin = 25                                     # peak search range (see manual)
oldcal    = 'det_corr_184_process_5.dat'
newcal    = 'det_merlin_cycle244.dat'
newrings  = 'MERLIN_rings_244.xml'
do_plot = False  # Plots door 1 or not
#======================================================================

def calculate_known_positions(inst):
# delta_slits is a list of the Cd strip thicknesses on the calibration bin (from the bottom)
# assumes the slits are 1 mm thick
    if inst.lower() == 'let':
        delta_slits = [79.5,27.,59.,60.,53.5,52.5]
        rad_ratio = 330 / 3.5
    if inst.lower() == 'merlin':
        #delta_slits = [78.,26.5,58.5,58.,53.,51]
        rad_ratio = 325 / 2.5
        delta_slits = [79.5,27.,59.,60.,53.5,52.5]
        

    zero = sum(delta_slits[:3]) + 2
    slits = [-zero]
    for islit in range(len(delta_slits)):
        slits.append(-zero+sum(delta_slits[:islit+1])+islit)
    slits = [x / rad_ratio for x in slits]
    return slits

#sum runs together
for irun in calibration_runs:
    w_buf = Load(str(irun))
    if irun == calibration_runs[0]:
        #print('Loading run #%i' % irun)
        ws = CloneWorkspace(w_buf)
    else:
       # print('... adding run #%i' % irun)
        ws = Plus(ws,w_buf)

# load mask if provided
if mask is not None:
    LoadMask('MERLIN',mask,OutputWorkspace='Masking')
    MaskDetectors(ws,MaskedWorkspace='Masking')
    
# keep copy of old calibration for comparison
ws_old = CloneWorkspace(ws)

# set up fit parameters
fitPar       = TubeCalibFitParams([45., 145., 180., 255., 330., 400., 460.], peak_height, peak_width)
fitPar_Door1 = TubeCalibFitParams([52., 148., 184., 253., 324., 394., 450.], peak_height, peak_width)
#fpD1_bad     = TubeCalibFitParams([75., 163., 197., 269., 342., 402., 450.], peak_height, peak_width) #Merlin bin
fpD1_bad     = TubeCalibFitParams([37., 151., 186., 256., 329., 396., 470.], peak_height, peak_width)
fpD9_bad     = TubeCalibFitParams([75., 148., 179., 255., 330., 397., 430.], peak_height, peak_width)
fitPar_nz    = TubeCalibFitParams([45., 145., 180., 330., 400., 460.], peak_height, peak_width)
fp_short_top = TubeCalibFitParams([95., 253., 398.], peak_height, peak_width)
fp_short_bot = TubeCalibFitParams([97., 327., 407.], peak_height, peak_width)

# perform peak fitting
ws = Integration(ws, RangeLower=integration_range[0],RangeUpper=integration_range[1])
known_pos    = calculate_known_positions('merlin')
# Because of the "wings" the edges of door 9 are smaller. For most tube the quadratic calibration is ok, but pack 3 tube 1 is problematic, use a different "known" position
known9_bad = [-1.10496]+known_pos[1:-1]+[0.990689]
known1_bad = known9_bad
kp_nozero    = [x for x in known_pos if x != 0.0]
peaks_form   = [2, 1, 1, 1, 1, 1, 2]
pf_nozero    = [2, 1, 1, 1, 1, 2]


#calibrationTable = calibrate(ws, 'MERLIN/door1', known_pos, peaks_form, fitPar=fitPar_Door1, margin=30, overridePeaks={16:[60,160,193,266,337,398,458]}, plotTube=range(0,31))
if do_plot:
    calibrationTable = calibrate(ws, 'MERLIN/door1', known_pos, peaks_form, fitPar=fitPar_Door1, margin=30, plotTube=range(0,31))
    calibrationTable = calibrate(ws, 'MERLIN/door1', known_pos, peaks_form, fitPar=fpD1_bad, margin=30, rangeList=[8,19,23], plotTube=[8,19,23], calibTable=calibrationTable)
else:
    calibrationTable = calibrate(ws, 'MERLIN/door1', known_pos, peaks_form, fitPar=fitPar_Door1, margin=30)
    calibrationTable = calibrate(ws, 'MERLIN/door1', known1_bad, peaks_form, fitPar=fpD1_bad, margin=30, rangeList=[8,19,23], calibTable=calibrationTable)
calibrationTable = calibrate(ws, 'MERLIN/door2', known_pos, peaks_form, fitPar=fitPar, margin=30, calibTable=calibrationTable)
calibrationTable = calibrate(ws, 'MERLIN/door3', [-0.38, 0.025, 0.41], peaks_form[4:7], fitPar=fp_short_top, margin=margin, rangeList=range(8), calibTable=calibrationTable)
calibrationTable = calibrate(ws, 'MERLIN/door3', [-0.405, 0.18, 0.38], peaks_form[0:3], fitPar=fp_short_bot, margin=margin, rangeList=range(8,16), calibTable=calibrationTable) 
calibrationTable = calibrate(ws, 'MERLIN/door3', kp_nozero, pf_nozero, fitPar=fitPar_nz, margin=margin, rangeList=range(16,19), calibTable=calibrationTable) 
calibrationTable = calibrate(ws, 'MERLIN/door3', known_pos, peaks_form, fitPar=fitPar, margin=margin, rangeList=range(19,40), calibTable=calibrationTable) 
calibrationTable = calibrate(ws, 'MERLIN/door4', known_pos, peaks_form, fitPar=fitPar, margin=margin, calibTable=calibrationTable) 
calibrationTable = calibrate(ws, 'MERLIN/door5', known_pos, peaks_form, fitPar=fitPar, margin=margin, calibTable=calibrationTable) 
calibrationTable = calibrate(ws, 'MERLIN/door6', known_pos, peaks_form, fitPar=fitPar, margin=margin, calibTable=calibrationTable) 
calibrationTable = calibrate(ws, 'MERLIN/door7', known_pos, peaks_form, fitPar=fitPar, margin=margin, calibTable=calibrationTable) 
calibrationTable = calibrate(ws, 'MERLIN/door8', known_pos, peaks_form, fitPar=fitPar, margin=margin, calibTable=calibrationTable) 
calibrationTable = calibrate(ws, 'MERLIN/door9', known_pos, peaks_form, fitPar=fitPar, margin=margin, calibTable=calibrationTable)
calibrationTable = calibrate(ws, 'MERLIN/door9', known9_bad, peaks_form, fitPar=fpD9_bad, margin=10, rangeList=[16], plotTube=[16], calibTable=calibrationTable)

#now - need to make a new calibration table to make sure that only the detector y-coordinates are changed
detpos = calibrationTable.column(1)
detID  = calibrationTable.column(0)
ycoord = []
for i in detpos: ycoord.append(i[1])

# Create CalibrationTable
MERLIN_calib = CreateEmptyTableWorkspace()
# Add required columns
MERLIN_calib.addColumn(type="int",name="Detector ID")
MERLIN_calib.addColumn(type="double",name="Detector Y Coordinate")

for j in range(len(detID)):
    nextRow = {'Detector ID': detID[j], 'Detector Y Coordinate': ycoord[j]}
    MERLIN_calib.addRow ( nextRow )

ApplyCalibration(ws, CalibrationTable=MERLIN_calib)
ModifyDetectorDotDatFile(InputWorkspace='ws',InputFilename=oldcal,OutputFilename=newcal)
rings=GenerateGroupingPowder(InputWorkspace=ws,AngleStep=0.5,GroupingFilename=newrings,GenerateParFile=True)

if do_plot:    
    for ii in range(8, 31):
        fig, axes = plt.subplots(edgecolor='#ffffff', num=ii, subplot_kw={'projection': 'mantid'})
        Data = mtd[f'TubePlot{ii}']
        Fitted = mtd[f'FittedTube{ii}']
        axes.plot(Fitted, color='#1f77b4', wkspIndex=0)
        axes.plot(Data, color='#ff7f0e', wkspIndex=0)
        axes.tick_params(axis='x', which='major', **{'gridOn': False, 'tick1On': True, 'tick2On': False, 'label1On': True, 'label2On': False, 'size': 6, 'tickdir': 'out', 'width': 1})
        axes.tick_params(axis='y', which='major', **{'gridOn': False, 'tick1On': True, 'tick2On': False, 'label1On': True, 'label2On': False, 'size': 6, 'tickdir': 'out', 'width': 1})
        axes.set_title(f'Fit {ii}')
        legend = axes.legend(fontsize=8.0).set_draggable(True).legend
        fig.show()
