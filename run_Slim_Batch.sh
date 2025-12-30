
mkdir output_slim_CHRI
mkdir output_slim_CHRII
mkdir output_slim_CHRIII
mkdir output_slim_CHRIV
mkdir output_slim_CHRV
mkdir output_slim_CHRX

for k in $(seq 1 5)
do
  p1=$((20*(k-1)+1))
  p2=$((20*k))

    echo "k=$k, p1=$p1, p2=$p2"
  #for i in $(seq "$p1" "$p2"); do
   # nohup slim -d trialNumber="$i" \
   #   /mnt/data3/mallard/2025_Breeding_values/SliM_simulations/Slim_V2_THOR.slim \
    #  > "output_slim_CHRI/SlimSimul1_$i.txt" &
  #done
  #wait

  for i in $(seq "$p1" "$p2"); do
    nohup slim -d trialNumber="$i" \
      /mnt/data3/mallard/2025_Breeding_values/SliM_simulations/Slim_V2_THOR_CHRII.slim \
      > "output_slim_CHRII/SlimSimul1_$i.txt" &
  done
  wait

  for i in $(seq "$p1" "$p2"); do
    nohup slim -d trialNumber="$i" \
      /mnt/data3/mallard/2025_Breeding_values/SliM_simulations/Slim_V2_THOR_CHRIII.slim \
      > "output_slim_CHRIII/SlimSimul1_$i.txt" &
  done
  wait

  for i in $(seq "$p1" "$p2"); do
    nohup slim -d trialNumber="$i" \
      /mnt/data3/mallard/2025_Breeding_values/SliM_simulations/Slim_V2_THOR_CHRIV.slim \
      > "output_slim_CHRIV/SlimSimul1_$i.txt" &
  done
  wait

  for i in $(seq "$p1" "$p2"); do
    nohup slim -d trialNumber="$i" \
      /mnt/data3/mallard/2025_Breeding_values/SliM_simulations/Slim_V2_THOR_CHRV.slim \
      > "output_slim_CHRV/SlimSimul1_$i.txt" &
  done
  wait

  for i in $(seq "$p1" "$p2"); do
    nohup slim -d trialNumber="$i" \
      /mnt/data3/mallard/2025_Breeding_values/SliM_simulations/Slim_V2_THOR_CHRX.slim \
      > "output_slim_CHRX/SlimSimul1_$i.txt" &
  done
  wait

done

