#ifndef _IO_CSV_HPP_
#define _IO_CSV_HPP_

#include <string>
#include <vector>

#include "field.hpp"

/*===========================================================================*/
class Csv 
{
public:
	void write(const class Field &field);
	
	Csv(const std::string output);
	Csv(const Csv &input);
	~Csv(void);
private:
	std::string output_dir;

};

#endif
