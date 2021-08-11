#ifndef _IO_SDF_HPP_
#define _IO_SDF_HPP_

#include <string>
#include <vector>
#include "field.hpp"

/*===========================================================================*/
//SDF header
#pragma pack(push, 1)
typedef struct Header{
	double time;
	double version;
	double rank;
	double dsize;
	double csize;
	double pnlen;
	double cnlen;
	double tglen;
} Header;
#pragma pack(pop)
/*===========================================================================*/
class Sdf
{
public:
	FILE* file;
        void write(const double time, const Field &f);

        Sdf(const std::string output, std::vector<double> space_array);
	~Sdf(void);
private:
        std::string output_dir;
	std::string cname;
	std::vector<double> coordinates;
	Header header;
	int space_points;
	double min;
	double max;
};

#endif
