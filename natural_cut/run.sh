./Preprocess/RemoveDup data/$1/USA-road-d.$1.gr data/$1/USA-road-d.$1-new.gr
./Filter/Filter $2 data/$1/USA-road-d.$1.co data/$1/USA-road-d.$1-new.gr data/$1/
./Assemble/Assemble $2 data/$1/USA-road-d.$1.co data/$1/USA-road-d.$1-new.gr  data/$1/anode.txt data/$1/aedge.txt data/$1/

