##################################
#
#  Format your data for SWIT-1.0
#  
##################################

import copy
import numpy as np
import obspy
from obspy.core import UTCDateTime
from obspy.core.util import get_example_file

# If your data in in the array format i.e., numpy array
dt = 0.004
srcx = 3800    # must be interger

your_data = np.loadtxt('shot020.dat').astype(np.float32)
your_recx = np.loadtxt('shot020_recx.dat').astype(np.int32)   # must be interger
your_recn = np.size(your_recx)
zero_trace = np.zeros_like(your_data[:,0])

# set the desired receiver range in the SWIT-1.0 code
rec_x_beg = 0.0
rec_x_end = 50e3
rec_dx    = 25
recn_in_swit = int((rec_x_end - rec_x_beg) // rec_dx) + 1
recx_in_swit = np.linspace(rec_x_beg, rec_x_end, recn_in_swit).astype(np.int32)   # must be interger

# get a example trace
stream = obspy.read(get_example_file("1.su_first_trace"), format='SU', byteorder='<')
tr_example = stream[0]

# initilize the empty traces
traces = []
for irec in range(recn_in_swit):
    # get the standar headers from the example trace
    tr = copy.deepcopy(tr_example)
    # put the zero trace
    tr.data = zero_trace
    # add headers (must)
    tr.stats.su.trace_header.group_coordinate_x  = recx_in_swit[irec]
    tr.stats.su.trace_header.source_coordinate_x = srcx
    tr.stats.sampling_rate = 1./dt
    tr.stats.distance = tr.stats.su.trace_header.group_coordinate_x - tr.stats.su.trace_header.source_coordinate_x
    # add headers (additional)
    tr.stats.station = '%d'%irec
    tr.stats.channel = 'Vz'
    tr.stats.network = 'SWIT-1.0 Acquisition'
    tr.stats.location = 'Your source index'
    tr.stats.starttime = UTCDateTime("2021-01-01T00:00:00")

    traces += [tr]

# put your own data into the empty traces with the corret SU headers
for irec in range(your_recn):
    itrace = np.argmin(abs(recx_in_swit - your_recx[irec]))
    traces[itrace].data = your_data[:, irec]

# organize into obspy stream 
traces = obspy.Stream(traces = traces)
traces.write('shot020_for_SWIT.su', format='SU')
