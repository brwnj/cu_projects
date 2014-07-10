# TODO
+ kmer norm testing

# modules
+ module use /mnt/scgc/software/modulefiles/common

# installing scgc-python
```
btupper@scgc_cruncher01 ~ $ cd /mnt/scgc/software/scgc_python
btupper@scgc_cruncher01 scgc_python [master]$ git pull origin master
btupper@scgc_cruncher01 scgc_python [master]$ sudo python setup.py install
```

# installing python2.7 and dependencies
https://github.com/0xdata/h2o/wiki/Installing-python-2.7-on-centos-6.3.-Follow-this-sequence-exactly-for-centos-machine-only

sqlite:
cp /usr/lib64/python2.6/lib-dynload/_sqlite3.so /usr/local/lib/python2.7/sqlite3/

for scipy:
yum --nogpgcheck install lapack lapack-devel blas blas-devel

for pandas and tables
yum install hdf5

python packages
toolshed
pandas
numpy
scipy
scgc...

yum groupinstall "Development tools"
yum install zlib-devel bzip2-devel openssl-devel ncurses-devel
wget --no-check-certificate https://www.python.org/ftp/python/2.7.7/Python-2.7.7.tar.xz
tar xf Python-2.7.7.tar.xz
cd Python-2.7.7
./configure --prefix=/usr/local
make
make altinstall
wget https://bitbucket.org/pypa/setuptools/raw/bootstrap/ez_setup.py
sudo /usr/local/bin/python2.7 ez_setup.py
sudo /usr/local/bin/easy_install-2.7 pip
