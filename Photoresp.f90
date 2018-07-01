
module Photoresp
use Allconst
use msimsl
    implicit none
    real(8),parameter::Eg=0.17 !eV
    real(8),parameter::eps_const=15.7
    real(8),parameter::constA=4.6E0
    real(8),parameter::J2erg=1.0E7
    real(8)::Vgup
contains    
    complex(8) function EpsSubstr(v)
        implicit none
        real(8)::pi,x,v, RealEps, ImagEps
        Vgup=Eg/hp_eV
        !print *,vgup*1.d-12
        x=v/Vgup
         if (x<1.0) then
         ImagEps=0
         RealEps=eps_const+constA*(2-sqrt(1+x)-sqrt(1-x))/x**2
          else
         ! ImagEps=constA*sqrt((hD_eV*w-Eg)*elCi*J2erg)/w*2
          ImagEps=constA*sqrt(x-1)/x**2
         RealEps=eps_const+constA*(2-sqrt(1+x))/x**2
         end if
        EpsSubstr=RealEps+c1*ImagEps
    end function
end module Photoresp