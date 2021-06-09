#!/usr/bin/env python

import datetime as dt
import numpy as np
from datetime import datetime
from datetime import timedelta
from amie_routines import *
import sys
import matplotlib.pyplot as plt
import matplotlib.dates as dates

# ------------------------------------------------------------------------
# These are stolen from Angeline's code
# ------------------------------------------------------------------------

def bool_string(line):
    """ Determine whether a string should be True or False
    Parameters
    ----------
    line : string
        Line to be tested
    Returns
    -------
    bout : bool
        Boolean output (True/False)
    Raises
    ------
    ValueError
        If the value cannot be interpreted as True or False
    Notes
    -----
    Accepts empty, true, false, t, f, 1, and 0 in any capitalization combo
    """

    line = line.lower()

    bout = None
    if line in ['true', 't', '1', '']:
        bout = True
    elif line in ['false', 'f', '0']:
        bout = False

    if bout is None:
        raise ValueError('input not interpretable as a boolean')

    return bout


def none_string(line):
    """ Determine whether a string should be None
    Parameters
    ----------
    line : string
        Line to be tested
    Returns
    -------
    out_line : string or NoneType
        None if all-lowercase version of line is "none" or line is zero length.
        Otherwise returns original value of line
    """
    out_line = None if line.lower() == "none" or len(line) == 0 else line
    return out_line


def process_command_line_input():
    """ Process command line input, needed to possible ipython use
    Returns
    -------
    input_args : list
        List of input arguements
    """

    input_args = sys.argv
    if input_args[0].find('ipython') >= 0:
        input_args = list()
    else:
        input_args.pop(0)

    return input_args

# ------------------------------------------------------------------------
# This sets the default values of the inputs and parses the inputs
# This is adapted from Angeline's code.
# ------------------------------------------------------------------------

def get_command_line_args(argv):
    """ Parse the arguements and set to a dictionary
    Parameters
    ----------
    argv : list
        List of arguments fed on the command line
    Returns
    -------
    args : dict
        A dictionary containing information about arguements, including:

    """
    # Initialize the arguments to their default values   

    args = {'outfile': 'test_ideal.bin',
            'dt': 1}

    arg_type = {'outfile': str,
                'dt': float}
    
    # If there is input, set default help to False
    args['help'] = False if len(argv) > 0 else True
    
    # Cycle through all arguments except the first, saving input
    for arg in argv:
        # Treat the file list and formatting seperately
        if arg.find('-') == 0:
            # This is not a filename, remove the dash to get the key
            split_arg = arg.split('=')
            akey = split_arg[0][1:]
            # Get the argument value as the desired type
            if akey not in arg_type.keys():
                raise ValueError(''.join(['unknown command line input, ',
                                          arg, ', try -help for details']))

            if len(split_arg) == 1:
                if arg_type[akey] == bool:
                    arg_val = True
                else:
                    raise ValueError('expected equality after flag {:}'.format(
                        akey))
            else:
                if arg_type[akey] == int:
                    arg_val = int(split_arg[1])
                elif arg_type[akey] == float:
                    arg_val = float(split_arg[1])
                elif arg_type[akey] == str:
                    arg_val = split_arg[1]
                else:
                    # This is boolean input
                    arg_val = bool_string(split_arg[1])

            args[akey] = arg_val
                    
    return args
    
# ------------------------------------------------------------------------
# main code
# ------------------------------------------------------------------------

# Get the input arguments
args = get_command_line_args(process_command_line_input())

time_1965 = datetime(1965, 1, 1, 0, 0, 0)

dt = args["dt"] # minutes
subtimes = (12.0 * 60.0 + np.arange(0, 60+dt, dt)) / 60.0
basetime = datetime(2020, 3, 20, 0, 0, 0)
alltimes = np.concatenate([[0.0], subtimes, [24.0]])
alltimes = np.array(alltimes) * 3600.0

nTimes = len(alltimes)

dLat = 1.0
minLat = 40.0
nLats = int((90.0 - minLat)/dLat + 1)
dMlt = 1.0/4.0
nMlts = int((24.0-0.0)/dMlt + 1)

lats = np.arange(minLat, 90.0+dLat, dLat)
mlts = np.arange(0.0, 24.0+dMlt, dMlt)

theta2d, r2d = np.meshgrid(mlts * np.pi/12.0 - np.pi/2.0, 90.0 - lats)

data = {}

data["nLats"] = nLats
data["nMlts"] = nMlts
data["nTimes"] = nTimes

data["lats"] = lats
data["mlts"] = mlts

data["times"] = []

data["imf"] = []
data["ae"] = []
data["dst"] = []
data["hp"] = []
data["cpcp"] = []

data["Vars"] = ['Potential (V)', \
                'Pedersen', \
                'Hall']

data["nVars"] = len(data["Vars"])

data["version"] = 0.99

for var in data["Vars"]:
    data[var] = []


for i in np.arange(0,nTimes):

    time_current = basetime + timedelta(seconds = alltimes[i])
    data["times"].append(time_current)

    time_delta = time_current - time_1965
    ut = time_delta.total_seconds() % 86400.0

    # Let's put a 30 minute signal:
    phase = np.pi * 2.0 * ut / 1800.0
    efield = 0.1 * np.cos(phase)  # V/m
    dpot = efield * 111.0 * 1000.0 * dLat
    print(ut/3600.0, phase, efield, dpot)

    pot2d = np.zeros([nLats, nMlts])
    pot2d[r2d < 20.0] = dpot
    hall2d = np.zeros([nLats, nMlts]) + 10.0
    ped2d  = np.zeros([nLats, nMlts]) +  5.0
    
    #  ; IMF should be (nTimes,4) - V, Bx, By, Bz
    #  ; AE  (nTimes, 4) - AL, AU, AE, AEI?
    #  ; Dst (nTimes, 2) - Dst, Dsti?
    #  ; Hpi (nTimes, 2) - HP, Joule Heating (GW)
    #  ; CPCP (nTimes) - CPCP (kV)

    data["imf"].append([0.0, 0.0, 0.0, 0.0])
    # These are all bad values for now....
    data["ae"].append([0.0, 0.0, 0.0, 0.0])
    data["dst"].append([0.0, 0.0])
    data["hp"].append([0.0, 0.0])
    data["cpcp"].append((np.max(pot2d) - np.min(pot2d))/1000.0)

    cPot = data["Vars"][0]
    data[cPot].append(pot2d)
    
    cPed = data["Vars"][1]
    data[cPed].append(ped2d)
    
    cHall = data["Vars"][2]
    data[cHall].append(hall2d)

# Write out in old AMIE format:
amie_write_binary(args["outfile"], data)

