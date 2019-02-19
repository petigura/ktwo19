while read dir
do 
    rsync -avz --progress cadence:/mir3/petigura/code/Phodymm/example_planets/${dir} analysis/photodyn/runs/

done < analysis/photodyn/runs.txt 
