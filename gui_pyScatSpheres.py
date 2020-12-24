#! /usr/bin/python3
import sys,os
from pyScatSpheres import gui_base as hsa_gui

args = sys.argv
if len(args)>1:
    df_path = args[1]
else:
    df_path=os.path.dirname(hsa_gui.__file__)+'/data/qdotSphereArray2_kp0.pkl'
    print(df_path)

hsa_e = hsa_gui.GUI_handler(df_path)
