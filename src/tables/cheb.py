#=============================================================================
## computes derivative matrices for Chebyshev polynomials T_n 
## at Chebyshev points  
#=============================================================================

import mpmath as mp
from typing import List

mp.prec= +mp.inf 

#=============================================================================
def cheb_pts(n:int)->List[mp.mpf]:
   return [mp.cos(mp.pi*i/(n-1)) for i in range(n)]
#=============================================================================
def diff_cheb(n:int,x:float)->mp.mpf:
   return mp.diff(lambda x: mp.chebyt(n,x), x)
#=============================================================================
def write_vec_to_file(name:str,n:int,vec:List[mp.mpf])->None:
   with open(name,'w') as f:
      for i in range(n):
         f.write(mp.nstr(vec[i],32)+'\n')
#=============================================================================
def write_matrix_to_file(name:str,n:int,mat:mp.matrix)->None:
   with open(name,'w') as f:
      for i in range(n):
         for j in range(n):
            f.write(mp.nstr(mat[i,j],32)+' ')
         f.write('\n')
#=============================================================================
def save_cheb(dir_name:str,n:int)->None:

   pts= cheb_pts(n)

   cheb_D_matrix= mp.matrix([[0 for x in pts] for n in range(n)]) 
#-----------------------------------------------------------------------------
   N= n-1
   cheb_D_matrix[0,0]=  mp.fdiv(2*mp.power(N,2) + 1,6) 
   cheb_D_matrix[N,N]= -cheb_D_matrix[0,0]
   cheb_D_matrix[0,N]=  0.5*mp.power(-1,N)
   cheb_D_matrix[N,0]= -0.5*mp.power(-1,N)
   
   for i in range(1,N):
      cheb_D_matrix[i,0]= -mp.fdiv(0.5*mp.power(-1,i), 1.0-pts[i])
      cheb_D_matrix[i,N]=  mp.fdiv(0.5*mp.power(-1,i), 1.0+pts[i])

      cheb_D_matrix[0,i]=  mp.fdiv(2*mp.power(-1,i), 1.0-pts[i])
      cheb_D_matrix[N,i]= -mp.fdiv(2*mp.power(-1,i), 1.0+pts[i])

      for j in range(1,N):
         if i!=j:
            cheb_D_matrix[i,j]= mp.fdiv(mp.power(-1,i+j), pts[i]-pts[j])
         else:
            cheb_D_matrix[i,i]= - mp.fdiv(pts[i], 2*(1+pts[i])*(1-pts[i]))
#-----------------------------------------------------------------------------
   inv_cheb_D_matrix = pow(cheb_D_matrix,-1)
#-----------------------------------------------------------------------------
   print(N)
   print(cheb_D_matrix)
   print(inv_cheb_D_matrix)
   print(inv_cheb_D_matrix*cheb_D_matrix)

   write_vec_to_file(
      "{}/cheb_pts.txt".format(dir_name),
      n, pts
   )
   write_matrix_to_file(
      "{}/cheb_D_matrix.txt".format(dir_name),
      n, cheb_D_matrix
   )
   write_matrix_to_file(
      "{}/inv_cheb_D_matrix.txt".format(dir_name),
      n, inv_cheb_D_matrix
   )
