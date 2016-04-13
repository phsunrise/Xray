import pickle
import os, sys
import numpy as np

def get_datadir():
    try:
        return pickle.load(open("parameters.pickle", 'r'))['data_dir']
    except KeyError:
        print "Error: Need to specify data_dir in \"parameters.pickle\""
        sys.exit(0)

def get_pads(pad):
    try:
        return pickle.load(open("parameters.pickle", 'r'))[pad]['pads']
    except KeyError:
        print "Error: Need to specify pads in \"parameters.pickle\""
        sys.exit(0)
