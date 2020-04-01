#include <iostream>
using std::cout; 
using std::endl;
#include <fstream>
using std::ifstream;
using std::ofstream;
#include <iomanip>
using std::setprecision;
#include <cassert>
#include <string>
using std::string;
#include <vector>
using std::vector; 

#include "io_csv.hpp"
#include "field.hpp"
	
/*===========================================================================*/
inline static bool exists(const string file_name)
{
	ifstream f(file_name);
	return f.good();
}
/*===========================================================================*/
Csv::Csv(const string output)
: output_dir{output}
{	
}
/*===========================================================================*/
Csv::Csv(const Csv &input)
: output_dir{input.output_dir}
{
}
/*===========================================================================*/
Csv::~Csv(void)
{
}
/*===========================================================================*/
void Csv::write(const Field &f)
{
	string file_name= output_dir+"/"+f.name+".csv";
	ofstream out;
	out.open(file_name,std::ios::app);
/*---------------------------------------------------------------------------*/
	if (out.is_open()) {
		int nx= f.np1.size();
		for (int i=0; i<nx-1; ++i) {
			out<<setprecision(16)<<f.np1[i]<<",";
		}		
		out<<f.np1[nx-1]<<endl;
	}
	else {
		cout<<"ERROR(Csv::write): "+file_name+" does not exist"<<endl;
	}
/*---------------------------------------------------------------------------*/
	out.close();
	return;
}
