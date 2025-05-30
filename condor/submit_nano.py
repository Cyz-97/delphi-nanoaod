import os
import shutil
import glob
import sys
import subprocess
import multiprocessing

fatmen = [
#    "xs_wphact211ncgg_e182.7_m80.4_c97_1l_g1",
#    "xs_gpym6143wc0eeqq_e182.7_c97_1l_g1",
#    "xs_qedbk23eegg_e183.5_l97_1l_g1",
    "xs_wphact24cc_e182.7_m80.4_c97_1l_g1"
#    "sh_qqps_e91.25_c94_2l_c2"
]

files = '/eos/opendata/delphi/simulated-data/cern/qqps/v94c/91.25/qqps*.sdst'

copyDir = '/eos/user/z/zhangj/DELPHI/simulation/1994_v2/qqps/'

if __name__ == "__main__":
    #os.system("source setup.sh")
    #os.system("cmake -B build")
    #os.system("cmake --build build")
    #nevt = 100
    for f in glob.glob(files):
        fname = os.path.basename(f)
        output = f"nanoaod_{fname}.root"

        print(f"Submitting job for {output}")

        # Create a copy of current environment and add your variables
        env = os.environ.copy()
        env["IN"] = f
        env["OUT"] = output
        env["COPYDIR"] = copyDir

        print(f, output, copyDir)

        subprocess.run(["condor_submit", "condor.sub"], env=env)

        sys.exit()

        
        #os.system(f'echo FILE={f} > dummy')
        #execution = f"build/delphi-nanoaod/delphi-nanoaod --nickname {nn} --mc --config config/delphi-nanoaod.yaml --output {output} -m {nevt}"
        #execution = f"build/delphi-nanoaod/delphi-nanoaod -P dummy --mc --config config/delphi-nanoaod.yaml --output {output}"
        #os.system(execution)
        #os.system(f"""root -q -b -l scripts/treefy.C+'("{output}")'""")
