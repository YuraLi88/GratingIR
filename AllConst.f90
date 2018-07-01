module AllConst
real(8),parameter::el=4.8d-10,mass=9.11d-28,meff=0.2d0,cl=3.0d10
real(8),parameter::elCI=1.602d-19,mCI=9.11d-31, hp=6.63E-34 ! J*s
real(8),parameter::hp_eV=4.14E-15 ! eV*s
real(8),parameter::hD_eV=6.595E-16 ! eV*s Dirac constant 
complex(8),parameter::c1=(0.0,1.0)  
real(8)::Pi  !!-----------------------------math const Pi
real(8)::lambda=299.792458D0 !!-wavelenth (mkm) at 1 THz


contains
!=============== THE COMMON FUNCTIONS=================
    complex(8) function KappaF(i,invV2,LA,dielF)
            use msimsl
             implicit none
             integer:: i ! positive index       
             real(8)::InvV2 ! - inverted square of a frequency
             real(8):: LA ! lambda0 to streep period
             complex(8):: dielF !! dielectric function
             complex(8):: z ! intermediate parameter
             complex(8)::results
    !!==========THE BODY OF FUNCTION====================
                z=LA*InvV2*i*i/dielF-1.d0
                if (Real(z)>0) then
                results=sqrt(z)
                else
                results=-sqrt(z)
                end if
            kappaF=results
    end function KappaF


end module

