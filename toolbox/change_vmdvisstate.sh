vmdfile=$1
tosed=$2
isovalue1=$3
isovalue2=$4
isovalue3=$5

sed 's/__TEMPLATE_CUBE__/'$tosed'/g' template.vmd | sed 's/__TEMPLATE_ISOVALUE1__/'$isovalue1'/g' | sed 's/__TEMPLATE_ISOVALUE2__/'$isovalue2'/g'  | sed 's/__TEMPLATE_ISOVALUE3__/'$isovalue3'/g'> $vmdfile

echo " " >> $vmdfile
echo "display reposition 0 0" >> $vmdfile
echo "display resize 1300 1300" >> $vmdfile
echo "display distance -1" >> $vmdfile
echo " " >> $vmdfile
echo "axes location off" >> $vmdfile
echo "render POV3 vmdscene.pov" >> $vmdfile
echo " " >> $vmdfile
echo "exit" >> $vmdfile
