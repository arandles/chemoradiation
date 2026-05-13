./make.cmd
cp program opt
cp program sub
cd opt
./run_all_opt.sh & 
cd ../sub
./run_all_sub.sh &
cd ..
