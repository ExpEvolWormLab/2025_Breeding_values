for k in $(seq 1 5)
do
  p1=$((20*(k-1)+1))
  p2=$((20*k))

    echo "k=$k, p1=$p1, p2=$p2"
  for i in $(seq "$p1" "$p2"); do
    nohup slim -d trialNumber="$i" \
      /mnt/data3/mallard/2025_Breeding_values/SliM_simulations/Slim_V2_THOR.slim \
      > "output_slim/SlimSimul1_$i.txt" &
  done

  wait
done
