infile=$1
imag_infile=$2
density_file=$3
state=$4
isovalue1=$5
isovalue2=$6
isovalue3=$7

cubefile=${density_file}_${state}.cube

echo " -- generating $density_file -- "
echo " "
julia ~/projects/numerovBand/Numerov/toolbox/generate_datfile3D.jl $infile $imag_infile $density_file $state

echo " -- generating $cubefile -- "
echo " "
julia ~/projects/numerovBand/Numerov/toolbox/generate_cubefile.jl $density_file $cubefile 20

~/projects/numerovBand/Numerov/toolbox/change_vmdvisstate.sh $cubefile.vmd $cubefile $isovalue1 $isovalue2 $isovalue3

vmd -e $cubefile.vmd

echo " -- adopting povray file -- "
echo " "
julia ~/projects/numerovBand/Numerov/toolbox/change_povray.jl vmdscene.pov $cubefile.pov

rm ${cubefile}.pov.tga
povray +W1005 +H1005 -I${cubefile}.pov -O${cubefile}.pov.tga +P +X +A +FT +C
