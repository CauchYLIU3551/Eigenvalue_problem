 function  [V,idx]=MGSg(B,X,U)
%
%---------------------------------------------------------------
%  Modified Gram-Schmit on X with selective re-orthogonalization
%  with respect to B-inner product
%---------------------------------------------------------------
%
%  Input:
%
%      B      (n-by-n) Hermitian and positive definite
%      X      (n-by-k)
%      U      (n-by-k0) U'*B*U=I already. Optional argument
%             Possibly k0=0, i.e., U is empty array
%
%  Output:
%
%      V      n-by-nV  B-Orthogonalized vectors from columns of X
%      idx    nV-by-1 (optional)
%             Some of X's columns may depend on its previous columns. These columns will be
%             dropped during the MGS process. This idx records the column numbers of X that
%             remain, i.e., those give V.  
%
%                                 RCL 02/20/2013
%                                 added idx 02/23/2014
%---------------------------------------------------------------
[n,k]=size(X);
rtol_re=1.0e-4; % relative tolerance to perform reorthogonalization
%
V = zeros(n,k); nV=1; idx=zeros(k,1);
if nargin==2,

   nrm2 = sqrt(real(X(:,1)'*B*X(:,1))); 
   if nrm2>0,
      V(:,nV)=X(:,1)/nrm2; idx(nV)=1;
      nV=nV+1;
   end
   for j=2:k,
     vh = X(:,j); nrm2 = sqrt(real(vh'*B*vh));
     tolorth=rtol_re*nrm2;
     %  by MGS
     for i=1:nV-1,
           vh = vh - V(:,i)*( V(:,i)'*B*vh );
     end
     nrm2=sqrt(real(vh'*B*vh));
     %  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     %  Perform re-orthogonalization once by MGS when deemed necessary.
     if nrm2 <= tolorth
        for i=1:nV-1,
           vh = vh - V(:,i)*( V(:,i)'*B*vh );
        end
        nrm2=sqrt(real(vh'*B*vh));
     end
     if nrm2>0
        V(:,nV)=vh/nrm2; %R(j,j)=nrm2;
        idx(nV)=j;
        nV=nV+1;
     end 
   end

elseif nargin==3,

   k0=size(U,2);
   vh = X(:,1); nrm2 = sqrt(real(vh'*B*vh));
   tolorth=rtol_re*nrm2;
   % by MGS
   for i=1:k0
       vh = vh - U(:,i)*( ( U(:,i) )'*B* vh );
   end
   nrm2=sqrt(real(vh'*B*vh));
   %  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
   %  Perform re-orthogonalization once by MGS when deemed necessary.
   if nrm2 < tolorth    
      for i=1:k0
          vh = vh - U(:,i)*( ( U(:,i) )'*B* vh );
      end
      nrm2=sqrt(real(vh'*B*vh)); 
   end
   if nrm2 > 0,
      V(:,nV)=vh/nrm2; idx(nV)=1;
      nV=nV+1;
   end 
   
   for j=2:k,
     vh = X(:,j); nrm2 = sqrt(real(vh'*B*vh));
     tolorth=rtol_re*nrm2;
     %  by MGS
     for i=1:k0
         vh = vh - U(:,i)*( ( U(:,i) )'*B* vh );
     end
     for i=1:nV-1,
           vh = vh - V(:,i)*( V(:,i)'*B*vh );
     end
     nrm2=sqrt(real(vh'*B*vh));
     %  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     %  Perform re-orthogonalization once by MGS when deemed necessary.
     if nrm2 <= tolorth
        for i=1:k0
            vh = vh - U(:,i)*( ( U(:,i) )'*B* vh );
        end
        for i=1:nV-1,
           vh = vh - V(:,i)*( V(:,i)'*B*vh );
        end
        nrm2=sqrt(real(vh'*B*vh));
     end
     if nrm2>0
        V(:,nV)=vh/nrm2; %R(j,j)=nrm2;
        idx(nV)=j;
        nV=nV+1;
     end
   end

end

nV=nV-1;
if nV<k
   V=V(:,1:nV); idx=idx(1:nV);
end
