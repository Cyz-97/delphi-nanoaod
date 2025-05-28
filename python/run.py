import os
import shutil
import glob
import sys
import subprocess
import multiprocessing

fatmen = [
    "xs_wphact211ncgg_e182.7_m80.4_c97_1l_g1",
    "xs_gpym6143wc0eeqq_e182.7_c97_1l_g1",
    "xs_qedbk23eegg_e183.5_l97_1l_g1",
    "xs_wphact24cc_e182.7_m80.4_c97_1l_g1"
]

if __name__ == "__main__":
    os.system("source setup.sh")
    os.system("cmake -B build")
    os.system("cmake --build build")
    nevt = 100
    for nn in fatmen:
        output = f"test_{nn}_evt{nevt}.root"
        execution = f"build/delphi-nanoaod/delphi-nanoaod --nickname {nn} --mc --config config/delphi-nanoaod.yaml --output {output} -m {nevt}"
        os.system(execution)
        os.system(f"""root -q -b -l scripts/treefy.C+'("{output}")'""")
