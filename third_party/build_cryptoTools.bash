cp miracl.get cryptoTools/thirdparty/linux/
cp linux64 cryptoTools/thirdparty/linux
cd cryptoTools/thirdparty/linux/
bash miracl.get
bash boost.get
cd ../../
rm CMakeCache.txt
cmake -DENABLE_MIRACL=ON -DCMAKE_POSITION_INDEPENDENT_CODE=ON .
make -j8
