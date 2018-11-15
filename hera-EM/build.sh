wget https://github.com/samtools/htslib/releases/download/1.3.2/htslib-1.3.2.tar.bz2 -O htslib.tar.bz2
tar -xjvf htslib.tar.bz2
mv htslib-1.3.2 htslib
cd htslib
make

cd ../
gcc -Ofast -fPIC -std=gnu99 *.c -Lhtslib/ -Ihtslib -lhts -lz -pthread -flto -lm -lstdc++ -o Hera-EM