cd cryptoTools/
rm CMakeCache.txt
cmake -DENABLE_MIRACL=ON -DCMAKE_POSITION_INDEPENDENT_CODE=ON .
make -j8
sudo make install
