!!! =========================================================
function fslop(aa,bb,cc)
implicit none
double precision::fslop,aa,bb,cc,minmod,ftemp
 fslop=dmax1(dsign(1.d0,(bb-aa)*(cc-bb)),0.d0)*  &
       dsign(1.d0,(cc-aa))*dmin1(dabs(0.5d0*(cc-aa)),  &
       dmin1(2.d0*dabs((bb-aa)),2.d0*dabs((cc-bb))))
 !  fslop = minmod(cc-bb,bb-aa)
end function

function vslop(aa,bb,cc) 
implicit none
double precision::vslop,aa,bb,cc
   vslop=(dsign(1.d0,(bb-aa))+dsign(1.d0,(cc-bb))) &
      *(dabs((bb-aa))*dabs((cc-bb)))/         &
       (dabs((bb-aa))+dabs((cc-bb))+1.0d-7)*1.d0
end function
!!! =========================================================

function minmod(aa,bb)
implicit none
double precision::minmod,aa,bb
    if(aa*bb .lt. 0.d0) then
       minmod = 0.d0
    elseif(abs(aa) .lt. abs(bb)) then
       minmod = aa
    else
       minmod = bb
    endif
end function
