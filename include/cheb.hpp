#ifndef _IO_CHEB_HPP_
#define _IO_CHEB_HPP_

#include <string>
#include <vector>

#include "field.hpp"

/*===========================================================================*/
class Cheb
{
public:
	
   Cheb(const string home_dir, const int Nx);
   ~Cheb(void);

   void r_derivative(
      const double cl,
      const vector<double> x_v,
      const vector<double> v, vector<double> dv
   ) const;

   void r_integral(
      const double cl,
      const vector<double> x_v,
      const vector<double> v, vector<double> iv
   ) const;

private:
   void get_cheb(
      const string dir, const string file_name, 
      std::vector<std::vector<double>>
   );

   void derivative(const std::vector<double> v, std::vector<double> dv) const;
   void integral(  const std::vector<double> v, std::vector<double> iv) const;

   const int nx;

   std::vector<std::vector<double>> cheb_mat; 
   std::vector<std::vector<double>> inv_cheb_mat;

};

#endif
