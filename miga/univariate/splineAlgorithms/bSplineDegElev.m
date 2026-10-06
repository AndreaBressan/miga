function [projMatrix,newKnots] = bSplineDegElev(d,c,k,t)

% BSPDEGELEV:  Degree elevate a univariate B-Spline. 
% 
% Calling Sequence:
% 
%   [projMatrix,newKnots] = bSplineDegElev(d,c,k,t)
% 
%   INPUT:
% 
%   d - Degree of the B-Spline.
%   c - Control points of the spline (dim, numCtrl).
%   k - Knot sequence, of the orginal spline space
%   t - Raise the B-Spline degree t times.
% 
%   OUTPUT:
%
%   projMatrix - projection matrix from coarse to fine B-spline spaces.
%   newKnots - Knot vector of the new B-Spline.
% 
 

[mc,nc] = size(c);
                                                                                                                     
n = nc - 1;                                                
                                                           
bezalfs =  zeros(d+1,d+t+1);                               
bpts = zeros(mc,d+1);                                     
ebpts = zeros(mc,d+t+1);                                  
Nextbpts = zeros(mc,d+1);                                  
alfs = zeros(d,1);                                      
                                                           
m = n + d + 1;                                            
ph = d + t;                                               
ph2 = floor(ph / 2);                                     

                                                         
bezalfs(1,1) = 1;                                          
bezalfs(d+1,ph+1) = 1;                                    

for i=1:ph2                                               
   inv = 1/nchoosek(ph,i);                                
   mpi = min(d,i);                                         
                                                          
   for j=max(0,i-t):mpi                                   
       bezalfs(j+1,i+1) = inv*nchoosek(d,j)*nchoosek(t,i-j);  
   end                                                       
end                                                     
                                                         
for i=ph2+1:ph-1                                          
   mpi = min(d,i);                                        
   for j=max(0,i-t):mpi                                  
       bezalfs(j+1,i+1) = bezalfs(d-j+1,ph-i+1);         
   end                                                       
end                                                      
                                                         
mh = ph;                                                 
kind = ph+1;                                             
r = -1;                                                  
a = d;                                                    
b = d+1;                                                  
cind = 1;                                                
ua = k(1);                                               
    
 
for ii=0:mc-1                                           
   projMatrix(ii+1,1) = c(ii+1,1);                           
end                                                      
for i=0:ph                                               
   newKnots(i+1) = ua;                                        
end                                                     
                                                        
for i=0:d                                               
   for ii=0:mc-1                                          
      bpts(ii+1,i+1) = c(ii+1,i+1);                       
   end                                                       
end                                                       
                                                          
while b < m                                               
   i = b;                                                 
   while b < m && k(b+1) == k(b+2)                       
      b = b + 1;                                          
   end                                                    
   mul = b - i + 1;                                        
   mh = mh + mul + t;                                      
   ub = k(b+1);                                            
   oldr = r;                                               
   r = d - mul;                                           
                                                          
                                                          
   if oldr > 0                                           
      lbz = floor((oldr+2)/2);                            
   else                                                 
      lbz = 1;                                           
   end                                                    
   
   if r > 0                                                
      rbz = ph - floor((r+1)/2);                         
   else                                                  
      rbz = ph;                                           
   end                                                    
   
   if r > 0                                               
                                                           
      numer = ub - ua;                                    
      for q=d:-1:mul+1                                    
         alfs(q-mul) = numer / (k(a+q+1)-ua);            
      end                                           
      
      for j=1:r                                           
         save = r - j;                                   
         s = mul + j;                                      
                                                          
         for q=d:-1:s                                      
            for ii=0:mc-1                                 
               tmp1 = alfs(q-s+1)*bpts(ii+1,q+1); 
               tmp2 = (1-alfs(q-s+1))*bpts(ii+1,q); 
               bpts(ii+1,q+1) = tmp1 + tmp2;              
            end                                              
         end                                               
         
         for ii=0:mc-1                                   
            Nextbpts(ii+1,save+1) = bpts(ii+1,d+1);       
         end                                                 
      end                                                
   end                                                    

                                                         
   for i=lbz:ph                                          
      for ii=0:mc-1                                       
         ebpts(ii+1,i+1) = 0;                             
      end                                                    
      mpi = min(d, i);                                    
      for j=max(0,i-t):mpi                               
         for ii=0:mc-1                                    
            tmp1 = ebpts(ii+1,i+1); 
            tmp2 = bezalfs(j+1,i+1)*bpts(ii+1,j+1);
            ebpts(ii+1,i+1) = tmp1 + tmp2;               
         end                                                 
      end                                                    
   end                                                    
                                                        
                                                         
   if oldr > 1                                            
                                                          
      first = kind - 2;                                              
      last = kind;                                        
      den = ub - ua;                                      
      bet = floor((ub-newKnots(kind)) / den);              
                                                            
                                                          
      for tr=1:oldr-1                                     
         i = first;                                       
         j = last;                                        
         kj = j - kind + 1;                             
         while j-i > tr                                   
                                                        

            if i < cind                                  
               alf = (ub-newKnots(i+1))/(ua-newKnots(i+1));        
               for ii=0:mc-1                            
                  tmp1 = alf*projMatrix(ii+1,i+1); 
                  tmp2 = (1-alf)*projMatrix(ii+1,i); 
                  projMatrix(ii+1,i+1) = tmp1 + tmp2;        
               end                                           
            end                                          
            if j >= lbz                                  
               if j-tr <= kind-ph+oldr                    
                  gam = (ub-newKnots(j-tr+1)) / den;     
                  for ii=0:mc-1                           
                     tmp1 = gam*ebpts(ii+1,kj+1); 
                     tmp2 = (1-gam)*ebpts(ii+1,kj+2); 
                     ebpts(ii+1,kj+1) = tmp1 + tmp2;     
                  end                                     
               else                                       
                  for ii=0:mc-1                           
                     tmp1 = bet*ebpts(ii+1,kj+1);                                     
                     tmp2 = (1-bet)*ebpts(ii+1,kj+2);                                     
                     ebpts(ii+1,kj+1) = tmp1 + tmp2;     
                  end                                        
               end                                       
            end                                           
            i = i + 1;                                    
            j = j - 1;                                   
            kj = kj - 1;                                  
         end                                             
                                                         
         first = first - 1;                             
         last = last + 1;                                 
      end                                                  
   end                                                 

   if a ~= d                                               
      for i=0:ph-oldr-1                                   
         newKnots(kind+1) = ua;                           
         kind = kind + 1;                                 
      end
   end                                                   
                                                          
                                                         
      for j=lbz:rbz                                     
         for ii=0:mc-1                                   
            projMatrix(ii+1,cind+1) = ebpts(ii+1,j+1);          
         end                                                 
         cind = cind + 1;                                 
      end                                                 
                                                           
      if b < m                                           
                                                          
         for j=0:r-1                                       
            for ii=0:mc-1                                 
               bpts(ii+1,j+1) = Nextbpts(ii+1,j+1);        
            end                                           
         end                                              
         for j=r:d                                         
            for ii=0:mc-1                                
               bpts(ii+1,j+1) = c(ii+1,b-d+j+1);          
            end                                              
         end                                                 
         a = b;                                            
         b = b+1;                                         
         ua = ub;                                         
                                                           
      else                                                
                                                          
         for i=0:ph                                       
            newKnots(kind+i+1) = ub;                      
         end                                                 
      end                                                    
end                                                      
% End big while loop                                       
end                                                      

 