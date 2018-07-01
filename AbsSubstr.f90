module AbsSubstr
use msimsl; use Allconst; use Photoresp;

implicit none

complex(8),allocatable::Ezm1D(:),Exm1D(:),Kappa2(:),KappaS(:),Rm(:)
complex(8)::Rm0,a0,b0
real(8),allocatable::LinSubstr(:,:), LinGrat(:,:)
real(8)::temp0,period0=0.5, Lcheck
!!------------THE MODULE SUBROUTINES-------------
contains
    real(8) function LossInSubstr(v,Ds,e,N)
    implicit none
        real(8)::v,Ds,E2 !! frequency (THz), Ds (mkm)
        real(8)::L !----- the output parameter
        integer::N,i,Nmax=30
        complex(8)::e,sE
    !!!---- intermediate variables---------
        real(8)::f0,x !!!----phase 
        complex(4)::y
    !!------------THE BODY OF FUNCTION-------------
        if (N<Nmax) then
            Nmax=N
        end if
        sE=sqrt(e)                 ! square root dielectric function
        f0=2*Pi*v*Ds/lambda  ! phase   
        x=2*f0*Imag(sE);                  y=-2*c1*f0*Real(sE);
        a0=0.5D0*(1-Rm0)*Exm1D(N+1);        b0=0.5D0*(1+Rm0)*Exm1D(N+1); 
        temp0=lambda/(period0*v)
        L=2*sum((/(PartialLoss(i,f0,sE,N),i=1,Nmax)/))          
        L=L+2*Real(a0*conjg(b0)*CEXPRL(y))
        a0=a0*conjg(a0); b0=b0*conjg(b0);
        L=L+a0*DEXPRL(x)+b0*DEXPRL(-x) !DEXPRL Evaluate the function (e^x-1)/x
        LossInSubstr=L*Imag(E)*f0
    end function LossInSubstr

    real(8) function PartialLoss(i,f0,sE,N)
    implicit none
        real(8)::f0,pL,x
        complex(8)::sE,a,b,az,bz
        complex(4)::y
        integer::i,N,i1
    !!!BODY OF function PartialLoss
        i1=N+1+i
        x=2*f0*Real(sE*KappaS(i));
        y=2*c1*f0*Imag(sE*KappaS(i))
        Ezm1D(i1)=-c1*temp0*i*Exm1D(i1)/(sE*kappaS(i))
        a  = 0.5D0*(1-Rm(i1))*Exm1D(i1);
        b  = 0.5D0*(1+Rm(i1))*Exm1D(i1);
        az= 0.5D0*(1-Rm(i1))*Ezm1D(i1);
        bz= -0.5D0*(1+Rm(i1))*Ezm1D(i1);
        pL=Real((a*conjg(b)+az*conjg(bz))*CEXPRL(y))
        pL=pL+Real((b*conjg(a)+bz*conjg(az))*CEXPRL(-y))
        a=a*conjg(a); b=b*conjg(b)
        az=az*conjg(az); bz=bz*conjg(bz)
        PartialLoss=pL+Real((a+az)*DEXPRL(x)+(b+bz)*DEXPRL(-x)) 
        !DEXPRL Evaluate the function (e^x-1)/x
    end function PartialLoss
end module AbsSubstr