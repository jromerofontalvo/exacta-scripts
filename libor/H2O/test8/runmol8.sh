base=libor
molecule=H2O
path=/n/home08/jromerofontalvo/libor/$molecule/integrals/
echo $path
frozen=0
virt=4
memory=8000
test=test8
outputfolder=outputs
days=7
prepmethod=exact
ntrotter=0
guesses=('read')
methods=(LBFGS)
usehf=1
hfstate='11001100'

for method in ${methods[*]}
do
for guess in ${guesses[*]}
do

counter=0
limit=6

for file in $path*
do
    let "counter += 1"
    echo "$counter"
    base1=${file//$path/}
    base2=${base1//'.int'/}
    base3=${base2//$molecule/}
    base4=${base3//$base/}
    base5=${base4//'-'/}
    echo "$base1 $base2 $base3 $base4 $base5"
    sed 's/DIST/'$base5'/g' mol8.template > temporal1 
    sed 's/frozen/'$frozen'/g' temporal1 > temporal2 
    sed 's/virt/'$virt'/g' temporal2 > temporal1 
    sed 's/MOLECULE/'$molecule'/g' temporal1 > temporal2 
    sed 's/MEMORY/'$memory'/g' temporal2 > temporal1 
    sed 's/DAYS/'$days'/g' temporal1 > temporal2 
    sed 's/BASE/'$base'/g' temporal2 > temporal1 
    sed 's/OUTPUTFOLDER/'$outputfolder'/g' temporal1 > temporal2 
    sed 's/OPTMETHOD/'$method'/g' temporal2 > temporal1 
    sed 's/NTROTTER/'$ntrotter'/g' temporal1 > temporal2 
    sed 's/GUESS/'$guess'/g' temporal2 > temporal1 
    sed 's/THRESHOLD/'$threshold'/g' temporal1 > temporal2 
    sed 's/USEHF/'$usehf'/g' temporal2 > temporal1 
    sed 's/HFSTATE/'$hfstate'/g' temporal1 > temporal2 
    sed 's/TEST/'$test'/g' temporal2 > temporal1 
    sed 's/PREPMETHOD/'$prepmethod'/g' temporal1 > inputs/$molecule-$base5-$base-$method-$prepmethod-$ntrotter-$guess$threshold-CAS-$frozen-$virt-8.sh    
    sbatch inputs/$molecule-$base5-$base-$method-$prepmethod-$ntrotter-$guess$threshold-CAS-$frozen-$virt-8.sh    
    if [ $counter -eq $limit ]; then
	let "counter = 0"
    fi
    # python variational6.py $molecule $base4 1 ucc $file > results_$molecule_STO6G/$base2.out
    # echo "$base1 is done ..."
done

done

done