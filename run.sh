seed=31728674
size_enclosure=1000000
time_step=0.01
samples=15

SOA_file="results_SOA.csv"
AOS_file="results_AOS.csv"

SOA="./sim-soa/build/sim-soa"
AOS="./sim-aos/build/sim-aos"

rm $SOA_file 2> /dev/null
rm $AOS_file 2> /dev/null

printf "num_objects,num_iterations,seed,size_enclosure,time_step" >> $AOS_file
printf "num_objects,num_iterations,seed,size_enclosure,time_step" >> $SOA_file

for ((i=1; i<=$samples; i++))
do
	printf ",Sample $i" >> $AOS_file
	printf ",Sample $i" >> $SOA_file
done

for num_objects in {1000,2000,4000}
do
	for num_iterations in {50,100,200}
	do
		t_AOS_sum=""
		t_SOA_sum=""
		for ((i=1; i<=$samples; i++))
		do
			exec_time=$($SOA $num_objects $num_iterations $seed $size_enclosure $time_step en_benchmark)
			t_SOA_sum+=",$exec_time"
			echo "SOA: num_objects: $num_objects, num_iterations: $num_iterations, sample: $i, exec_time: $exec_time ms"
			exec_time=$($AOS $num_objects $num_iterations $seed $size_enclosure $time_step en_benchmark)
			t_AOS_sum+=",$exec_time"
			echo "AOS: num_objects: $num_objects, num_iterations: $num_iterations, sample: $i, exec_time: $exec_time ms"
		done
		printf "\n$num_objects,$num_iterations,$seed,$size_enclosure,$time_step$t_AOS_sum" >> $AOS_file
		printf "\n$num_objects,$num_iterations,$seed,$size_enclosure,$time_step$t_SOA_sum" >> $SOA_file
	done
done
