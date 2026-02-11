#!/bin/bash
# 11-02-2026
# alex
# q2f.sh

funQ=()
funD=()
funQ+=("powq" "sinq" "cosq" "sqrtq" "fabsq" "log10q" "34Qe")
funD+=("pow" "sin" "cos" "sqrt" "fabs" "log10" "19e")
fitx=()

mkdir float_dir
cp *{.c,.h} float_dir/
cp makefile float_dir/
cd float_dir
for fh in *.h ; do
    NOM=$(basename "$fh" ".h")
    fitx+=("$NOM")
    while IFS= read -r NFUNC; do
	NFNOVA=${NFUNC}_d
	sed -i "s/\b$NFUNC\b/$NFNOVA/g" "$NOM.h"
	sed -i "s/\b$NFUNC\b/$NFNOVA/g" "$NOM.c"
	funQi=$(echo "$NFUNC" | awk '{print $NF}')
	funDi=$(echo "$NFNOVA" | awk '{print $NF}')
	funQ+=("$funQi")
	funD+=("$funDi")
    done < <(grep ");" "$fh" | cut -d '(' -f1 | grep -v "typedef")
done

for f in *.c ; do
    len=${#funQ[@]}
    for ((i=0; i<len; i++)); do
	sed -i "s/\b${funQ[i]}\b/${funD[i]}/g" "$f"
    done	 
done

sed -i "s/\b__float128\b/real/g" "vector.h"
sed -i 's/strtoflt128(v1, NULL)/atof(v1)/g' "vector.c"
sed -i 's/strtoflt128(num, NULL)/atof(num)/g' "vector.c"

make
len=${#fitx[@]}

for ((i=0; i<len; i++)); do
    mv "${fitx[i]}.so" "../${fitx[i]}_d.so"
done

cd ..
rm -rf float_dir
