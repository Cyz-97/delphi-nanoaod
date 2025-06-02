import os
import shutil
import glob
import sys
import subprocess
import multiprocessing
from pathlib import Path
import time
import re
import yaml

opendata_dir = Path("/eos/opendata/delphi")
user_output_base = Path("/eos/user/z/zhangj/DELPHI")

def build_patterns(nickname):
    cfg = config[nickname]
    
    if cfg["type"] == "data":
        patterns = [
            str(opendata_dir / "collision-data" / subdir / "**" / f"*{cfg['extension']}")
            for subdir in cfg["subdirs"]
        ]
        copy_dir = user_output_base / "collision_data" / cfg["version"] / cfg["energy"]
        TYPE = "DATA"
    else:
        pattern = str(
            opendata_dir / "simulated-data" / cfg["origin"] / "**" / cfg["version"] /
            "**" / f"{cfg['stream']}_*{cfg['energy']}*{cfg['extension']}"
        )
        patterns = [pattern]
        copy_dir = user_output_base / "simulation" / cfg["version"] / cfg["energy"] / cfg["stream"]
        TYPE = "MC"

    print("Search pattern:", patterns)
    return patterns, str(copy_dir), TYPE

def find_matches(patterns):
    matches = []
    for pattern in patterns:
        matches.extend(glob.glob(pattern, recursive=True))
    return matches

def file_should_be_skipped(path, MIN_SIZE_BYTES, MAX_AGE_DAYS):
    if not os.path.exists(path):
        return False
    stat = os.stat(path)
    size_ok = stat.st_size > MIN_SIZE_BYTES
    age_ok = (time.time() - stat.st_mtime) < (MAX_AGE_DAYS * 86400)
    return size_ok and age_ok

def total_jobs(user):
    try:
        output = subprocess.check_output(
            ["condor_q", user, "-total"],
            stderr=subprocess.DEVNULL,
            text=True
        )
        for line in output.splitlines():
            if line.startswith("Total for query:"):
                # Example line: "Total for query: 200 jobs; 0 completed; 0 removed; 200 idle; 0 running; 0 held; 0 suspended"
                numbers = re.findall(r"\b\d+\b", line)
                # Only count idle + running jobs (as per previous context)
                if len(numbers) >= 5:
                    return int(numbers[3]) + int(numbers[4])
        return 0
    except Exception as e:
        print(f"Error checking condor queue: {e}")
        return 0

def load_config(yaml_path="condor/sample_list.yaml"):
    with open(yaml_path, "r") as f:
        return yaml.safe_load(f)

if __name__ == "__main__":

    min_bytes = 100_000      # 100 KB
    max_days = 30

    MAX_QUEUE = 100
    USER = os.environ["USER"]

    config = load_config()
    
    #nickname = "short94_c2"
    #nickname = "sh_qqps_e91.25_c94_2l_c2"
    nickname = "sh_kk2f4146qqpy_e91.25_c94_2l_c2"

    patterns, copy_dir, run_type = build_patterns(nickname)
    matched_files = find_matches(patterns)

    print(f"Found {len(matched_files)} files for '{nickname}':")
    for f in matched_files:
        print(f)
    print(f"Copy directory: {copy_dir}")

    os.makedirs(copy_dir, exist_ok=True)

    # Pause for manual check
    input("👉 Press Enter to start job submission...")
    
    for f in matched_files:
        fname = os.path.basename(f)
        output = f"nanoaod_{fname}.root"

        output_path = os.path.join(copy_dir, output)

        if file_should_be_skipped(output_path, min_bytes, max_days):
            print(f"⚠️ Skipping {output}")
            continue

        # To avoid ran out of disk space on lxplus
        while total_jobs(USER) >= MAX_QUEUE:
            print(f"⏳ Too many jobs in Condor queue for {USER}. Waiting 1 minutes...")
            time.sleep(60)
            
        print(f"✅ Submitting job for {output}")

        # Create a copy of current environment and add your variables
        env = os.environ.copy()
        env["IN"] = f
        env["OUT"] = output
        env["COPYDIR"] = copy_dir
        env["ISMC"] = run_type

        print(f, output, copy_dir)

        subprocess.run(["condor_submit", "condor/condor.sub"], env=env)
