./make.cmd
cp program opt
cp program sub
cp program soc
cp program rand
cd opt
rm *csv
./run_all_opt.sh & 
cd ../sub
rm *csv
./run_all_sub.sh &
cd ../soc
rm *csv
./run_all_soc.sh &
cd ../rand
rm *csv
./run_all_rand.sh &
