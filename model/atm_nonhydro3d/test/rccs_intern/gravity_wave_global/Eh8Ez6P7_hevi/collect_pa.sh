for I in `seq 1 17` ; do
    fapppx -A -d ./prof/rep${I} -p0,limit=20 -Icpupa -tcsv -o pa${I}.csv
done