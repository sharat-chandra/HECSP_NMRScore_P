for file in pose_*.mol2
do
	base=${file%.mol2}
	num=`echo "$file" | cut -d'.' -f1 | cut -d'_' -f2`
	$AMBERHOME/bin/antechamber -fi mol2 -fo mol2 -i $base.mol2 -o gas_pose_"$num".mol2 -j 0 -pf y -c gas -nc 0
done
