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
   return mp.diff(lambda x: mp.chebyt(n,x) x)
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
def save_cheb_D_matrix(dir_name:str,n:int)->None:

   pts= cheb_pts(n)

   cheb_D_matrix= mp.matrix([
      [diff_cheb(n,x) for x in pts]
      for n in range(n)
   ]) 

   inv_cheb_D_matrix = pow(cheb_D_matrix,-1)

   write_vector_to_file(
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
