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

   void set_x(
      const double cl, const double x_exc,
      std::vector<double> x_v
   ) const;

private:
   void get_cheb_pts(
      const string dir, const string file_name, 
      std::vector<double> pts
   );
   void get_cheb_mat(
      const string dir, const string file_name, 
      std::vector<std::vector<double>> mat
   );

   void derivative(const std::vector<double> v, std::vector<double> dv) const;
   void integral(  const std::vector<double> v, std::vector<double> iv) const;

   const int nx;

   std::vector<double> cheb_pts; 
   std::vector<std::vector<double>> cheb_mat; 
   std::vector<std::vector<double>> inv_cheb_mat;

};

#endif
