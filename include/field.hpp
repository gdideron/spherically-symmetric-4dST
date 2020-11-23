#ifndef _FIELD_HPP_
#define _FIELD_HPP_

#include <string>
#include <vector>

/*===========================================================================*/
class Field
{ 
public:
   const std::string name;

   const std::string type;

   std::vector<double> n;
   std::vector<double> np1;

   std::vector<double> inter_2; 
   std::vector<double> inter_3; 
   std::vector<double> inter_4;

   std::vector<double> rDer;

   Field(
      const std::string file_name,
      const int nx, const double init_val
   );
   ~Field(void);

   void set_to_val(const int min, const int max, const double val);
   void shift_time_step(void);
   void check_isfinite(const double time);
private:
   const int nx;

   void check_field_isfinite(
      const double time, const std::string message,
      const std::vector<double> &vec);
};

#endif
