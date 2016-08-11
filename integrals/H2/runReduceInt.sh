path=${PWD##*}
for file in $path*
do
    echo $file
    python reduceInt.py $file 3 3
done
