#!/bin/bash

source /cvmfs/sft.cern.ch/lcg/views/LCG_104/x86_64-el9-gcc12-opt/setup.sh


#export generator=kk2f4146_qqpy
export generator=qqps
#export generator="."

export year=1994

# Set tag and ECM list based on year
case $year in
  2000)
    tag="v00e1"
    ecm_list="204 205 207 208"
    ;;
  1999)
    tag="v99e1"
    ecm_list="192 196 200 202"
    ;;
  1998)
    tag="v98e1"
    ecm_list="189"
    ;;
  1994)
    tag="v94c"
    ecm_list="91.25"  # define as needed
    ;;
  *)
    echo "❌ Unsupported year: $year"
    exit 1
    ;;
esac

if [[ "$generator" != "." ]]; then
    folder="simulation"
else
    folder="collision_data"
fi


check_file_is_good_by_size() {
    local f="$1"
    local min_bytes=100000  # 100 kB
    local fullpath="/eos/user/z/zhangj/DELPHI/${folder}/${tag}/${ecm}/${generator}/${f}"

    if [ ! -f "$fullpath" ]; then
        echo "❌ BAD (file not found): $fullpath"
        return 1
    fi

    size=$(stat -c%s "$fullpath" 2>/dev/null)

    if [ -z "$size" ]; then
        echo "❌ BAD (stat failed): $fullpath"
        return 1
    elif [ "$size" -lt "$min_bytes" ]; then
        echo "❌ BAD (too small): $fullpath ($size bytes)"
        return 1
    else
        echo "✅ GOOD (size = $size bytes): $f"
        return 0
    fi
}

# Loop over ECM values
for ecm in $ecm_list; do
    dir="/eos/user/z/zhangj/DELPHI/${folder}/${tag}/${ecm}/${generator}"
    echo "🔍 Processing ECM = $ecm in $dir"

    if [ ! -d "$dir" ]; then
        echo "❌ Directory does not exist, skipping..."
        continue
    fi

    cd "$dir" || continue
    rm -f allfiles.txt chunk_* 2>/dev/null

    ls *${generator}*.root > allfiles.txt || { echo "⚠️ No files found for ECM $ecm"; continue; }

    split -l 50 allfiles.txt chunk_

    i=1
    for f in chunk_*; do
        out=$(printf "nanoaod_%s_%s_merged_%02d.root" "$generator" "$ecm" "$i")
        echo "📦 Preparing $out from $f"

        goodfiles=()
        while read -r file; do
            fullpath="root://eosuser.cern.ch//eos/user/z/zhangj/DELPHI/${folder}/${tag}/${ecm}/${generator}/${file}"
            echo "🔎 Checking: $fullpath"
            if check_file_is_good_by_size "$file"; then
                echo "✅ GOOD: $file"
                goodfiles+=("$fullpath")
            else
                echo "❌ BAD: $file"
            fi
        done < "$f"

        if [ ${#goodfiles[@]} -gt 0 ]; then
            echo "✅ Merging ${#goodfiles[@]} good files into $out"
	    echo "👉 Will hadd these files:"
	    printf '%s\n' "${goodfiles[@]}"
            hadd -fk "$out" "${goodfiles[@]}"
        else
            echo "⚠️ No valid files in $f, skipping..."
        fi

        i=$((i + 1))
    done

    echo "✅ Done with ECM = $ecm"
done
