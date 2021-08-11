#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <fstream>

#include "io_sdf.hpp"

/*===========================================================================*/
// Byte order reversal method from
// https://stackoverflow.com/questions/41220686/
double byte_reversal(double f){
	auto it = reinterpret_cast<uint8_t*>(&f);
	std::reverse(it, it + sizeof(f));
	return f;
}
/*===========================================================================*/
int endian_check(){
        unsigned int i = 1; 
	int is_big;
        char *c = (char*)&i; 
        if (*c) {
		is_big = 0;
	}
    	else{
		is_big = 1;
	}
    return is_big; 
}
bool bigendian = endian_check();
/*===========================================================================*/
std::vector<double> byte_reversal(std::vector<double> a){
	std::transform(a.begin(), a.end(), a.begin(), 
		[](double d) -> double { return byte_reversal(d);});
	return a;
}
/*===========================================================================*/
void fwrite(std::vector<double> a, FILE* file){
	if (!bigendian){
		a = byte_reversal(a);
	}
	fwrite(&a[0],sizeof(std::vector<double>::value_type),
			a.size(),file);
}
/*===========================================================================*/
void fwrite(double f, FILE* file){
	if (!bigendian){
		f = byte_reversal(f);
	}
	fwrite(&f,sizeof(double),1,file);
}
/*===========================================================================*/
Sdf::Sdf(const std::string output, std::vector<double> space_array)
: 	output_dir{output},
  	cname{"x"},
	coordinates{space_array},
	space_points{int(space_array.size())},
	min{space_array[0]},
	max{space_array[space_points-1]}
{	
	header.version = 1.0;
	header.rank = 1;
	header.dsize = space_points; 
	header.csize = space_points; 
	header.cnlen = cname.size(); 
	header.tglen = 0.0; 
	if (!bigendian) {
		header.time = byte_reversal(header.time);
		header.version = byte_reversal(header.version);
		header.rank = byte_reversal(header.rank);
		header.dsize = byte_reversal(header.dsize);
		header.csize = byte_reversal(header.csize);
		header.cnlen = byte_reversal(header.cnlen);
		header.tglen = byte_reversal(header.tglen);
	}
}
/*===========================================================================*/
Sdf::~Sdf(void){
}
/*===========================================================================*/
void Sdf::write(const double t, const Field &f){
	std::string DataName = f.name;
	header.pnlen = DataName.size(); 
	std::string file_name= output_dir+"/"+f.name+".sdf";
	file = fopen(file_name.c_str(),"ab");
	header.time = t;
	if (!bigendian) {
		header.time = byte_reversal(header.time);
		header.pnlen = byte_reversal(header.pnlen);
	}
	fwrite(&header, sizeof(Header), 1, file);
	fwrite(DataName.c_str(),DataName.size(),1,file);
	fwrite(cname.c_str(),cname.size(),1,file);
	fwrite(min,file);
	fwrite(max,file);
	fwrite(space_points,file);
	fwrite(coordinates,file);
	fwrite(f.np1,file);
	fclose(file);
}
