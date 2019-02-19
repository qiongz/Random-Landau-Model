#include"wfs_file.h"
void write_wfs(int theta_1, int theta_2, complex<dtype> * wfs, unsigned long size){
        stringstream suffix;
        suffix <<"wfs_"<< theta_1 << "_"<< theta_2<<".bin";
        char *filename;
        filename = new char[suffix.str().length()+1];
        strcpy(filename, suffix.str().c_str());
	ofstream fout(filename,ios_base::binary|ios_base::out);
	fout.seekp(0);
	fout.write((char*)(wfs),sizeof(complex<dtype>)*size);
}

void read_wfs(int theta_1, int theta_2, complex<dtype> * wfs,unsigned long size){
        stringstream suffix;
        suffix <<"wfs_"<< theta_1 << "_"<< theta_2<<".bin";
        char *filename;
        filename = new char[suffix.str().length()+1];
        strcpy(filename, suffix.str().c_str());
	ifstream fin(filename,ios_base::binary|ios_base::in);
	fin.seekg(0);
	fin.read((char*)wfs,size*sizeof(complex<dtype>));
	fin.close();
}
