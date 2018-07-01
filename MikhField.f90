!  MikhField.f90 
!
!  FUNCTIONS:
!    MikhField      - Entry point of console application.
!
!****************************************************************************
!
!  PROGRAM: MikhField
!
!  PURPOSE:  Entry point for the console application.
!
!**************************************************************************** 


module inData
    implicit none
    integer(8),parameter::Np=500  !-discretization of profile (for arbitrary profile)
    real(8)::v1=0.1, v2=2.5       !-----------THz
    integer(8)::Npoints=1000     !--discretization of frequency
    real(8)::ZgratMin=0.2d0;
    real(8)::d_cut=0.01 ! mkm
    real(8)::n2Dup=0.0D12 !!!---2DEG above the grating (on the bkanket "substrate")
!Substrate parameters    
    complex(8)::eps0=1.D0,eps1=1.0D0,eps2=1.0D0,epsS=1.0D0,eps3=1.D0
    real(8)::Db=0.1d0,d=0.02d0,Ds1=4.d0, Ds !mkm (!)
    real(8),parameter::Zmin=-5.0, Zmax=5.0, v2D=1.1
!2D gas parameters    
    real(8)::t2D=4.0d0    
    real(8)::n2D0=1.0d0,N2D
    real(8):: vdr0=0.0d7 ! cm/s
    real(8):: Edr=0.E2      ! V/cm
!Grating parameters 
    real(8)::tgrat,SigMetal=4.D17, Dgrat=4.D-6
    real(8)::fill=0.8d0, period=0.1, wds
    real(8)::Sig1D=6.d17,dg=2.0d-6, n1D=0.0D10
!FileName==================================        
    character(*),parameter::DirName='data\TSTmaping\' !'
    character(*),parameter::FName='tstInverseVdr0'  
    real(8)::tPhaze=0.3
    integer::Time,Tmax=1
    integer,parameter::Npx=700, Npz=700
!===== Governing Parameters================
    real(8),dimension(3)::Vch=(/20.0d0,52.97d0,88.95d0/)
    integer::NewSigmaCheb=1,realtau=0, wiresgrat=0
    integer::orto=0, Tau2D_Inf=0, BuildField=0
    integer::inverse_geom=0, LossInSubstrCalc=0
    integer::epsInSb=0
    complex(8)::free_part
    !=====================
end module inData

module Field2D
    use numcalc;
    implicit none
    real(8)::ThetaIn=45, ThetaInR, z, za,AbsEx,AbsEz, AbsE
    integer::jz
    complex(8),allocatable::sinIJ(:,:)
    complex(8),dimension(-N:N)::AmCf,BmCf,CmCf,HmCf,FmCf,ksi2_m
    real(8),dimension(0:Npx+1)::WxzOnX,Xpts
    contains
        complex(8) function HmCalc(i,v,k2_m,k3_m) 
            use msimsl
            use InData
            use Allconst
            implicit none
            integer :: i
            complex(8) :: k2_m,k3_m
            complex(8) :: b23,ksi_m
            complex(8) ::f2,cf1,cf2
            real(8) :: pi,v
            !-Body of function
                pi=dconst('pi')
                b23=sqrt(eps2/epsS)*k3_m/k2_m
                f2=2*pi*d*sqrt(eps2)*v/lambda
                cf1=(1.0+b23*(1-ksi2_m(i)))*BmCf(i)
                cf2=(1.0-b23*(1+ksi2_m(i)))*CmCf(i)
                HmCalc=cf1*exp(-f2*k2_m)+cf2*exp(f2*k2_m)
                HmCalc=HmCalc/2.0
        end function HmCalc

        complex(8) function FmCalc(i,v,k2_m,k3_m) 
            use msimsl
            use InData
            use Allconst
            implicit none
            integer :: i
            complex(8) :: k2_m,k3_m
            complex(8) :: b23,ksi_m
            complex(8) ::f2,cf1,cf2
            real(8) :: pi,v
            !-Body of function
                pi=dconst('pi')
                b23=sqrt(eps2/epsS)*k3_m/k2_m
                f2=2*pi*d*sqrt(eps2)*v/lambda
                cf1=(1.0-b23*(1-ksi2_m(i)))*BmCf(i)
                cf2=(1.0+b23*(1+ksi2_m(i)))*CmCf(i)
                FmCalc=cf1*exp(-f2*k2_m)+cf2*exp(f2*k2_m)
                FmCalc=FmCalc/2.0
        end function FmCalc
end module Field2D

module Profil
    use InData
    implicit none
    integer(8)::ip,ip1
    character(*),parameter::Fcijname='datW\cijcheb.txt'
    character(*),parameter::Fprofilname='data\Profiles\111.dat' !c_ns035_N2_12_L2_8T300
    character(*),parameter::Ftstname='data\Profiles\currprofile.dat'
    real(8),allocatable::cij(:,:)    
    integer::nint,ndata,pw;    
    real(8),allocatable::breack(:),cscoef(:,:)
    real(8)::normaPr
    real(8)::Hmag=0.d0,Wcr
end module Profil


module blanket ! data for upper layer
    implicit none
        Complex(8),allocatable::Kappa1(:),Kapa0(:)
        complex(8)::Etot;
        complex(8)::Kappa0up,KappaMup, gamma0Up,GammaUp;
        character(5)::strTime1;
        character(3)::strTime;
        integer::envMax=0 !Search envelope of maxima and write in file
!*****************************************************************
end module blanket

module numcalc   
    !dec$ attributes DLLIMPORT::ChebN 
    ! ChebN(n,x) n-integer - is an order of polinom , x -real(8)
    use AllConst; use msimsl;
    use inData;
    implicit none
    complex(8)::Sig1DonV !!! paramer for procedure of integration Power dissipation in the streep
    complex(8),allocatable::A(:,:)   !matrix of main equation
    complex(8),allocatable::Cn(:)    !Coef. of field on streep, solution of Matrix Equation
    integer,parameter::N=100        !Number of fourier garmonic in sum
    integer,parameter::Neq=5     !Number of Equations
    integer::j,k,l
    integer::ndata
    real(8)::tau2D,invtau=0.0,Mu,Vdr,tau
    character(*),parameter::VAfname='v_16-15_T030.txt'
    real(8),dimension(Neq)::IntN
    character(*),parameter::DName='.dat'
    character(*),parameter::Mkh='7' 
    !!============Module Procedures================
    contains
        complex(8) function ExOnStr(x)
            implicit none
            real(8)::x,ChebN
            complex(8)::res
            integer::i
            !Body of function
            res=0
                do i=1,Neq
                res=res+Cn(i)*ChebN(i-1,x)
                end do
            ExOnStr=res
        end function ExOnStr
        
        complex(8) function Power1Dtot()
            implicit none
            !integer::Ir=2 ! irule parameter - low of integration
            !real(8)::a=-1.D0, b=1.D0, Res
            !real(8)::eAbs=0.D-3,eRel=1.D-2,eRes
            !external::Power1D
            !CALL DQDAG(Power1D,a,b,eAbs,eRel,Ir,res,eres)
          Power1Dtot=Real(sig1DonV)*sum(Cn*conjg(Cn))
        end function
end module numcalc
    real(8) function Power1D(x)
        use numcalc; use Allconst
        implicit none
        complex(8)::Ex
        real(8)::x
        Ex=ExOnStr(x)
        Power1D=Real(sig1DonV*sqrt(1-x*x))*Ex*conjg(Ex)
    end function
!**************************************************The body of the program************************************************************************
program MikhField
    use msimsl;    use numcalc;use AllConst;use blanket; use profil; use Field2D;
    USE Photoresp; use AbsSubstr
    use inData
implicit none
    ! Functions
    complex(8)::alpha0,alphaM,reflB,psi0,SigDispMagn
    ! Variables
    integer(8),parameter::IPATH=1,Ncell=20000,Nas=170,Nfcmax=1, NwrFcof=7 !!Number of Fourier coefs, what is writed in file
    integer::LDA,M,ios,sft,i,NN,i1,j1,m1,Nchb,m2
    complex(8),allocatable::W(:),W1(:),R2m(:),Wbeta(:),Ex(:),Ez(:)    
    real(8),allocatable::F(:,:),R(:),Tr(:),TrMkh(:),dU(:),P2Dm(:),TrUp2D(:)
    complex(8),allocatable::B(:),Exm2D(:),Ezm2D(:),Sigma(:), ChargeOn2Dm(:)
    complex(8)::eps,W0,gama,gama2,gamaM,gamaM0,Ex1D,ChargeOn2D,Ez1D,PEx,Jx,dzetta,psi00, gSigma2D0
    complex(8)::theta,theta0,ksi,ksi0,Q0,P0,t0,t2,expr1,expr2
    real(8)::tmp,tmp1,par 
    complex(8)::f1,f2,f10,f20,f3,f30 !
    complex(8)::b10,b20,bs0,b120,b3s,b3s0,b1,b2,bs,b12
    real(8)::v,Mah,gama02,p,streep,wds1
    real(8)::betaG,dtV,results,x,Px,Pint,ChebNder
    complex(8)::Ro1D
    character(*),parameter::    RFName='refl\R'
    character(*),parameter::    TFName='Trans\T'     
    character(*),parameter::    LFName='loss\L'            
    character(5)::Vstr
    character(3)::Vstr1
    !=========chebishev method
    real(8)::lowlim=0.0d0,UpLim=4.5d1,errabs=1.0d-11,errrel=1.0d-11
    real(8)::h,AsintR,errres,AsInt,Psum,Det20
    complex(8)::Det10
    real(8),allocatable::BslN(:),Det2(:,:)
    complex(8),allocatable::Slk(:,:),Fn(:),Bm0(:),Matr(:,:),Bm(:,:),Bjkl(:,:,:),Bjkl1(:,:,:) 
    complex(8),allocatable::LUfac(:,:),wk(:),Det1(:,:)
    integer,allocatable::ipvt(:)
    integer:: PhotorespCalc=1, DetOut=0, time_anim=0 
    external::BesU,mutau
    complex(8)::BettaGn
    !=========================
    !----The--external-cycle
    integer,parameter::Npts=1,Nfmax=4
    integer::Ickl
    real(8)::Var1=2.0d12 ,dVar=1.0d-3,ShftK=0.5d-1
    real(8),dimension(Npts)::vars
    real(8),allocatable::AllTrans(:,:)
    character(8)::Prefix
    character(3)::Pref,Pref0
    character,parameter::namV='n'
    real(8)::RatioLA
    !=========================
        real(4)::start_time,finish_time
    !-----------------------------------------------------------------------------------------------------------------------
    !
! =============Body of MikhField======================================
    !
    !
    !=================================================================    
    wds1=fill*period
    Vgup=Eg/hp_eV
    ALLOCATE(AllTrans(Npts,Npoints),Bjkl(N,neq,neq),Bjkl1(N,neq,neq))
    Allocate(Ex(0:Npx+1),Ez(0:Npx+1))
    allocate(Bm0(Neq),BslN(Neq+1),Bm(Neq,N),Sigma(N+1)) 
    allocate(Slk(Neq,Neq))
    allocate(Det1(Npts,Npoints),Det2(Npts,Npoints),Wbeta(Neq))
    allocate(cij(Neq,Neq),sinIJ(2*N+1,Npx+2))    
    !CALL weight
    pi=dconst('pi'); p=1.d0/6.d0; dtV=(v2-v1)/Npoints
    ThetaInR=Pi*ThetaIn/180.d0
        
    if (NewSigmaCheb==1) then
            write(*,*), "New method ass. sumation"
            else
            CALL OrtoPolinom
    end if ! (NewSigmaCheb==1)
    wds=wds1
    streep=wds/fill
    write(*,*), 'streep width=',streep,' (mkm)'    
    Ds=d+Ds1
    !write(*,*), 'd=',d;write(*,*), ',      Ds1=',Ds1; write(*,*), ',       Ds=',Ds
    !=============== fill the matrix of exp: START
    if (buildfield==1)    then
            do j=0,Npx+1
                 Xpts(j)=j/((Npx+1)*1.0)
                do i=-N,N
                    i1=i+N+1
                    sinIJ(i1,j+1)=cos(2*Pi*i*(j/((Npx+1.0))-fill/2))+c1*sin(2*Pi*i*(j/((Npx+1.0))-fill/2))
                end do !var i
            end do ! var j
    end if ! (buildfield==1)
    !=============== fill the matrix of exp    : END
    !------Begin of External cycle-----
    write(*,*) 'D/a=',d/streep
    do Ickl=1,Npts
        CALL cpu_time(start_time)
        if (realtau==1) then
            call mutau
                if (Edr<0.5d0)  then
                        vdr=0.d0
                else
                        tau2D=tau !*1.0d-10;    !Vdr=0.180111d0;
                end if
                invtau=1.d-2/tau2D
         else 
                tau2D=t2D*1.d-12
                if (Tau2D_inf==0)  then
                        invtau=1.d0/tau2D
                end if
                Vdr=0
              
        end if ! (realtau==1)
        Vdr=0.0 !7.5E-2 
        tau2D=10.17
        invtau=1.d12/tau2D
        print *,'invtau=',invtau
        print *,'mu=',mu
        print *,'tau2D=',tau2D 
        print *,'Vdr=',Vdr
        n2D=n2D0*1.0d12
        gama02=2*Pi*n2D*el*el/(mass*meff*cl)
        gamma0Up=2*Pi*n2Dup*el*el/(mass*meff*cl) ! gamma for 2DEG in inverse geometry
        gSigma2D0=n2D*el*el/(mass*meff)
        !gama02=2*Pi*n2D*tau2D*el*el/(mass*meff*cl) !!! old gama2D
        print *,'!!!gama02=',gama02
        print *,'!!!gamma0Up=',gamma0Up
        print *,'gSigma2D0=',gSigma2D0
        print *, 'm=',mass
        print *, 'el=',el
        print *, 'n1D=',n1D
        print *, 'tauGrat=',tgrat
        print *, 'invtau=',invtau
        if (wiresgrat==1) then
                Sig1D=el*el*n1D*(tgrat*1.D-12)/(mass*meff) !*(4/Pi)
                else
                Sig1D=SigMetal*Dgrat !*(4/Pi)
        end if
        print *, 'Sig1D=',Sig1D
        !gamaM0=Pi*Pi*fill*Sig1D*dg/(2*cl*sqrt(eps1)) ! (!) for metallic grating
        if (orto==1) then
            gamaM0=2*Pi*fill*sig1D/(cl*sqrt(eps1))*normaPr !*Pi/4.0 ! *normaPr
            else
            gamaM0=2*Pi*fill*sig1D/(cl*sqrt(eps1))*Pi/4.0 
        end if
        print *,'gamaM0=',gamaM0
        !gama0=2*Pi*Sig1D*dg/cl
        Wcr=el*Hmag*1.D4/(meff*mass*cl)
        write(*, *),'Wcr='
        write(*, 5),Wcr/(2*Pi)
        write(*,*), 't=',t2D
        allocate(R(Npoints),Tr(Npoints),TrUp2d(Npoints),TrMkh(Npoints),Exm2D(2*N+1),Ezm2D(2*N+1),Ezm1D(2*N+1),dU(Npoints))
        allocate(LinSubstr(2,Npoints),LinGrat(2,Npoints))
        allocate(ChargeOn2Dm(2*N+1),Exm1D(2*N+1))
        allocate(Fn(Neq),Cn(Neq),Matr(Neq,Neq),LUfac(Neq,Neq),Ipvt(Neq),wk(Neq)) 
        if ((Ickl==1).or.(NamV=='f')) then
            Bm0(1:Neq)=(/1,(0,l=2,Neq)/)
            write(*,*),'=============================================================='
                if (orto==1) then
                    do j=1,N
                        Bm(1:Neq,j)=(/(BettaGn(l,j),l=1,Neq)/)
                    end do
                else
                    do j=1,N
                        CALL DBSJNS(Pi*j*fill,Neq+1,BslN)
                        Bm(1:Neq,j)=(/(2*l*c1**(l-1)*BslN(l+1)/(Pi*j*fill),l=1,Neq)/)
                    end do
                end if
            do l=1,Neq
                do k=1,Neq
                    if (NewsigmaCheb==1) then
                            if (l==k) then
                                    h=1.0d0/(l+k)
                            else
                                    h=2.0d0*sin(Pi*(k-l)/2.0d0)/(Pi*(k+l)*(k-l))
                            end if
                        CALL DQDAGS(BesU,lowlim,UpLim,errabs,errrel,AsintR,errres)
                        AsInt=Real(1.0/(Pi*fill)*(h+2.0*c1**(k+l)*AsintR)) !
                        Asint=2*real(4*k*l*c1**(l-k)*AsInt)
                    end if !!-----(NewsigmaCheb==1)-------------
                        Psum=0.0d0
                        do j=1,N
                                Psum=Psum+(Pi*j*fill)*(Bm(l,j)*conjg(Bm(k,j))+Bm(k,j)*conjg(Bm(l,j))) !BslN(k+1)*BslN(l+1)/(Pi*j*fill)
                                Bjkl(j,k,l)=Bm(l,j)*conjg(Bm(k,j))
                                Bjkl1(j,k,l)=Bm(k,j)*conjg(Bm(l,j))
                        end do
                        if (NewSigmaCheb==1) then
                                Slk(l,k)=(AsInt-Psum)*(lambda/(Pi*streep*sqrt(eps1)*fill))*2/(1+eps2/eps1)
                    end if
                end do !!for K
            end do !! for L
        end if !!for Slk calk
        !call DWRRRl('Slk',Neq,Neq,Slk,Neq,0,'(F15.8)','','')
        if (N>0) then
            allocate(W(2*N+1),W1(2*N+1),Rm(2*N+1),R2m(2*N+1),kapa0(N),Kappa1(N),Kappa2(N),KappaS(N))
        end if
        allocate(a(2*N+1,2*N+1),b(2*N+1),P2Dm(2*N+1))
        par=lambda/streep;        
        print *,'par=',par
        RatioLA=par*par
        !open(27,file=Dirname//'Psum\Psum'//Fname)
        open(27,file=Dirname//'Ph'//Fname//'.dat')
        open(28,file=Dirname//'PhE'//Fname//'.dat')
        open(29,file=Dirname//'PExMap'//Fname//'.dat')
        do k=1,Npoints
            v=v1+((k-1)*(v2-v1))/Npoints
            if (EpsInSb==1) then
                EpsS=EpsSubstr(v*1.D12)
                Eps2=EpsS
            end if ! (EpsInSb==1)
            b10=sqrt(eps1/epsS)
            b3S0=sqrt(eps3/epsS)
            b20=sqrt(eps2/epsS)
            b120=sqrt(eps2/eps1)
            gama2=gama02/(invtau-2*pi*c1*v*1.d12)
            GammaUp=gamma0Up/(invtau-2*pi*c1*v*1.d12) ! for inverse geometry 2DEG     
            gamaM=gamaM0/(1-2*pi*c1*v*tgrat)
            gamaM=SigDispMagn(gamaM0,tgrat,v)
            gama=gamaM        
            !There is a calculation of the phase shift 
            !    write(27,*), v,sign(acos(Real(gamaM)/abs(gamaM)),1.d0)/Pi
            LDA=2*N+1
            par=1.d0/(v*v)
            ksi0=2.0d0*gama2/sqrt(epsS)
            Kappa0up=2.0d0*GammaUp/sqrt(eps1) ! the analog of ksi0 for upper geometry 2DEG
            f10=2.0d0*pi*sqrt(epss)*(Ds-d)*v/lambda
            f20=2.0d0*pi*sqrt(eps2)*d*v/lambda
            theta0=(b3S0*cos(f10)-c1*sin(f10))/(cos(f10)-c1*b3S0*sin(f10))
            expr1=(ksi0+b20+theta0)*exp(-c1*f20)
            expr2=(ksi0-b20+theta0)*exp(c1*f20)
            Rm0=(expr1+expr2)/(expr1-expr2)
            W0=2/(alpha0(v)+b120*Rm0)
            psi00=psi0(v)
            Q0=b3S0*cos(f10)-c1*sin(f10)
            P0=cos(f10)-c1*b3S0*sin(f10)
            t0=b20/(b20*P0*cos(f20)-c1*(Q0+ksi0*P0)*sin(f20))  !exp(-c1*f20)
            !!!!!!!!!!!!!!!!!!!!!!!!t2=abs(t0*W0)*abs(t0*W0)
            do i=-N,N
                    i1=i+N+1
                    gama2=gama02/((1-i*Vdr/(v*streep))*(invtau-2*pi*c1*v*1.d12*(1-i*Vdr/(v*streep))))
                    GammaUp=gamma0Up/((1-i*Vdr/(v*streep))*(invtau-2*pi*c1*v*1.d12*(1-i*Vdr/(v*streep))))
                    if (i==0) then
                        Rm(N+1)=Rm0
                        W(N+1)=W0 
                        f30=-2*pi*sqrt(eps2)*d*v/lambda
                        R2m(i1)=cos(f30)-Rm(i1)*sin(f30) 
                        ksi2_m(i)=2*gama2/sqrt(eps2)
                        else                                        
                        Kappa2(abs(i))=kappaF(i,par,RatioLA,eps2)
                        KappaS(abs(i))=kappaF(i,par,RatioLA,epsS)              
                        Kappa1(abs(i))=kappaF(i,par,RatioLA,eps1)
                        Kapa0(abs(i))=kappaF(i,par,RatioLA,eps0)
                        ksi=2*c1*gama2*kappaS(abs(i))/sqrt(epsS)
                        ksi2_m(i)=2*c1*gama2*kappa2(abs(i))/sqrt(eps2)
                        kappaMUp=2*c1*gammaUp*kappa1(abs(i))/sqrt(eps1) ! the analog of KSI for upper geometry 2DEG
                        b1=b10*(kappaS(abs(i))/kappa1(abs(i)))
                        b2=b20*(kappaS(abs(i))/kappa2(abs(i)))
                        b3S=b3S*(kappaS(abs(i))/kappa2(abs(i))) 
                        b12=b120*(kappa1(abs(i))/kappa2(abs(i)))
                        f1=2*pi*sqrt(epsS)*kappaS(abs(i))*(Ds-d)*v/lambda
                        f2=2*pi*sqrt(eps2)*kappa2(abs(i))*d*v/lambda
                        Theta=(b3S+ztanh(f1))/(1+b3S*ztanh(f1))
                        Rm(i1)=(b2*ztanh(f2)+ksi+Theta)/(b2+(ksi+Theta)*ztanh(f2))
                        R2m(i1)=2*b2*exp(-f2)/((b2+(ksi+Theta)*ztanh(f2))*(1+exp(-2*f2)))
                        W(i1)=2/(alphaM(v,i)+b12*Rm(i1))  
                        W1(i1)=kappa1(abs(i))*W(i1)
                    end if !!!(i==0) 
            end do !! cycle for N
            Fn(1:Neq)=(/psi00*W0,((0.0+c1*0.0),l=2,Neq)/)
            do l=1,Neq
                do i=1,Neq
                    Sigma(N+1)=sum((/(W1(N+1+j)*Bjkl(j,i,l)+W1(N+1-j)*Bjkl1(j,i,l),j=1,N)/))
                    !Sigma(1:N)=(/(W1(N+1+j)*Bm(l,j)*conjg(Bm(i,j))+W1(N+1-j)*conjg(Bm(l,j))*Bm(i ,j),j=1,N)/)
                    !Sigma(N+1)=0;    !Kappa(j)*(W(N+1+j)*Bm(l,j)*conjg(Bm(i,j)+W(N+1-j)*conjg(Bm(i,j))*Bm(l,j)))
                    if (NewSigmaCheb==1) then
                        Matr(l,i)=gamaM*(W0*Bm0(l)*conjg(Bm0(i))+c1*(Sigma(N+1)+Slk(l,i)/v))
                    else
                        Matr(l,i)=gamaM*(W0*Bm0(l)*conjg(Bm0(i))+c1*Sigma(N+1))
                    end if
                    if (i==l) then 
                        Matr(l,i)=Matr(l,i)+1.d0
                    end if
                end do ! i=1,Neq
            end do ! l=1,Neq
            !CALL DWRCRL('Matr',Neq,Neq,Matr,Neq,0,'','','')
            CALL DL2ACG(Neq,Matr,Neq,Fn,IPATH,Cn,lufac,Ipvt,wk)
            CALL DLFDCG(Neq,lufac,Neq,ipvt,det10,det20)
            !det1(Ickl,k)=log10(abs(det10))+det20
            !det2(Ickl,k)=acos(real(det10)/abs(det10))*sign(1.0,Imag(det10))
            Etot=W0*(psi00-gamaM*Cn(1))
            !old expression    R(k)=abs(1-W0*(1-gamaM*Cn(1)))*abs(1-W0*(1-gamaM*Cn(1)))
            R(k)=abs(reflB(v))**2
            write(28,*), v,sign(acos(Real(Cn(1))/abs(Cn(1))),Imag(Cn(1)))/Pi
            TrMkh(k)=abs(t0*W0*(psi00-gamaM*Cn(1)))*abs(t0*W0*(psi00-gamaM*Cn(1))) !abs(w0*psi00)**2 !
            TrUp2D(k)=abs(Etot)**2
          !  if (PhotorespCalc==1) then
                Mah=Vdr/(streep*v)
                dU(k)=0    
                par=lambda/streep
                do i=-N,N
                    i1=i+N+1
                    if (i==0) then
                        Exm1D(i1)=W(N+1)*(1-gamaM*Cn(1))
                        Exm2D(i1)=(cos(2*Pi*v*d/lambda)+c1*Rm0*sin(2*Pi*v*d/lambda))*W(N+1)*(1-gamaM*Cn(1)) 
                        Ezm1D(i1)=0
                        Ezm2D(i1)=0
                        ChargeOn2Dm(i1)=0
                        else
                        f2=2*Pi*d*sqrt(eps2)*kappa2(abs(i))*v/lambda
                        Wbeta=Bm(:,abs(i))
                            if (i>0) then
                                Wbeta=Conjg(Wbeta)
                            end if
                        Exm1D(i1)=-c1*gamaM*W1(i1)*sum(Cn*Wbeta) 
                        Ezm1D(i1)=c1*i/sqrt(i*i*1.d0-v*v*eps2/(par*par))*Exm1D(i1)
                        AmCf(i)=Exm1D(i1)
                        BmCf(i)=0.5D0*(1+Rm(i1))*Exm1D(i1)
                        CmCf(i)=0.5D0*(1-Rm(i1))*Exm1D(i1)
                        HmCf(i)=HmCalc(i,v,kappa2(abs(i)),KappaS(abs(i)))
                        FmCf(i)=FmCalc(i,v,kappa2(abs(i)),KappaS(abs(i)))
                        if ((d_cut>d) .and. (d_cut<Ds1)) then
                                    f2=2*Pi*KappaS(abs(i))*sqrt(epsS)*v*(d_cut-d)/lambda
                                    if (abs(f2)<15.d0) then
                                        Exm2D(i1)=HmCf(i)*exp(-f2)+FmCf(i)*exp(f2)
                                        Ezm2D(i1)=c1*i/sqrt(i*i*1.0-v*v*eps2/(par*par))*Exm2D(i1)
                                    else
                                        Exm2D(i1)=0.d0;    Ezm2D(i1)=0.d0
                                    end if     ! f2<10       
                                else
                                    f2=2*Pi*Kappa1(abs(i))*v*d_cut/lambda
                                    if (abs(f2)<15) then
                                            Exm2D(i1)=0.5*((1+Rm(i1))*exp(-f2)+(1-Rm(i1))*exp(f2))*Exm1D(i1) 
                                            Ezm2D(i1)=c1*i/sqrt(i*i*1.0-v*v*eps1/(par*par))*Exm2D(i1)
                                        else
                                            Exm2D(i1)=0; Ezm2D(i1)=0;    
                                    end if !(abs(f2)<15)                                        
                        endif    ! ((d_cut>d) .and. (d_cut<Ds1))
                        if (d_cut <0) then
                            f2=2*Pi*Kapa0(abs(i))*sqrt(eps0)*v*d_cut /lambda 
                            if (abs(f2)<15.d0) then
                                Exm2D(i1)=AmCf(i)*exp(f2)
                                Ezm2D(i1)=c1*i/sqrt(i*i*1.0-v*v*eps2/(par*par))*Exm2D(i1)
                            else
                                Exm2D(i1)=0.d0;    Ezm2D(i1)=0.d0
                            end if     ! f2<10                                                  
                        endif !(d_cut <0)
                        if ((d_cut >0) .and. (d_cut <d) ) then
                            f2=2*Pi*Kappa2(abs(i))*sqrt(eps2)*v*d_cut /lambda
                            if (abs(f2)<15.d0) then
                                Exm2D(i1)=BmCf(i)*exp(-f2)+CmCf(i)*exp(f2)
                                Ezm2D(i1)=c1*i/sqrt(i*i*1.0-v*v*eps2/(par*par))*Exm2D(i1)
                            else
                                Exm2D(i1)=0.d0;    Ezm2D(i1)=0.d0
                        end if     !f2<10                                                  
                        endif !((d_cut >0) .and. (d_cut <d) )
                        if ((d_cut >d) .and. (d_cut <Ds1)) then
                            f2=2*Pi*KappaS(abs(i))*sqrt(epsS)*v*(d_cut -d)/lambda
                            if (abs(f2)<15.d0) then
                                Exm2D(i1)=HmCf(i)*exp(-f2)+FmCf(i)*exp(f2)
                                Ezm2D(i1)=c1*i/sqrt(i*i*1.0-v*v*eps2/(par*par))*Exm2D(i1)
                            else
                                Exm2D(i1)=0.d0;    Ezm2D(i1)=0.d0
                            end if     ! f2<10                                               
                        endif !((d_cut >d) .and. (d_cut <Ds1))

                        ChargeOn2Dm(i1)=(2.0*i)/(period*v*1.d8)*gSigma2D0*Exm2D(i1)/((1-i*Vdr/(v*streep))*(invtau-2*pi*c1*v*1.d12*(1-i*Vdr/(v*streep))))
                    end if !(i==0)
                    gama2=gama02/((1-i*Vdr/(v*streep))*(invtau-2*pi*c1*v*1.d12*(1-i*Vdr/(v*streep))))   
                    if ((Tau2D_Inf==0) .and. (invtau .NE. 0)) then
                        P2Dm(i1)=2*Real(gama2/invtau)*abs(Exm2D(i1))*abs(Exm2D(i1))
                        !P2Dm(i1)=2*Real(gama2)*abs(Exm2D(i1))*abs(Exm2D(i1))
                    end if
                    Sig1DonV=Pi*Pi*fill*sig1D/(1-2*Pi*v*c1*tgrat)/cl;    
                    LinSubstr(1,k)=v/(Vgup*1.D-12)
                    if (LossInSubstrCalc==1) then
                        LinSubstr(2,k)=LossInSubstr(v,Ds,EpsS,N)
                        LinGrat(2,k)=Power1Dtot()
                        else
                        LinSubstr(2,k)=0.0
                        LinGrat(2,k)=0.0
                    endif ! (LossInSubstrCalc==1)
                   
                        if (abs(i)<40 ) then 
                            dU(k)=dU(k)+(i/(1-i*Mah))*(Exm2D(i1)*conjg(Exm2D(i1)))/(1+(2*Pi*v*t2D*(1-Mah*i))**2)
                        end if
                end do ! i=-N,N
                dU(k)=-16*Pi*Ncell*1.0e7*(((t2D/3)*7.903e6)/(v*1.0e12))*dU(k)*300/(3.0e10)
               ! else
                !    dU(k)=0
           ! end if !(PhotorespCalc==1)
            !===========================
            
            if ((k>4) .and. (buildField==1)) then 
                !if (((abs(Vch(1)-v)<dtV/1.9).or.((abs(Vch(2)-v)<dtV/1.9))).or.((TrMkh(k)>TrMkh(k-1)).and.(TrMkh(k-1)<TrMkh(k-2))))  then !)
                if (((abs(Vch(1)-v)<dtV/1.9).or.((abs(Vch(2)-v)<dtV/1.9))).or.((LinSubstr(2,k)<LinSubstr(2,k-1)).and.(LinSubstr(2,k-1)>LinSubstr(2,k-2))))  then !)
                    print *,'                                                   '
                    print *,'                                                   '
                    ! if ((TrMkh(k)>TrMkh(k-1)).and.(TrMkh(k-1)<TrMkh(k-2))) then 
                    if ((LinSubstr(2,k)<LinSubstr(2,k-1)).and.(LinSubstr(2,k-1)>LinSubstr(2,k-2))) then
                        print *,'v=',v,'  -- extremum'
                        print *,'eps=',epsS
                    end if !((LinSubstr(2,k)<LinSubstr(2,k-1)).and.(LinSubstr(2,k-1)>LinSubstr(2,k-2)))
                    Vstr=''
                    write(Vstr,11) v 
                    Vstr1=Vstr(2:3)//Vstr(5:5)
                    write(*,*),'v/vg=',v/vgup*1.D12
                    open(2,file=Dirname//'Ro1D/reRo1D'//Vstr1//fname//'.dat')
                    open(3,file=Dirname//'Ro1D/imRo1D'//Vstr1//fname//'.dat')
                    open(4,file=Dirname//'Eon2D\ReEx'//Vstr1//Fname//'.dat');
                    open(5,file=Dirname//'Eon2D\ImEx'//Vstr1//Fname//'.dat');
                    open(6,file=Dirname//'EonGrat\ReP'//Vstr1//Fname//'.dat')
                    open(10,file=Dirname//'EonGrat\ImP'//Vstr1//Fname//'.dat')
                    open(11,file=Dirname//'Curr\ReJ'//Vstr1//Fname//'.dat')
                    open(12,file=Dirname//'Curr\ImJ'//Vstr1//Fname//'.dat')
                    open(13,file=Dirname//'Ezon2D\ReEz'//Vstr1//Fname//'.dat');
                    open(14,file=Dirname//'Ezon2D\ImEz'//Vstr1//Fname//'.dat');
                    do i=-NwrFcof, NwrFcof
                        write(7,*), i*1.d0, abs(Exm2D(i+N+1))
                        write(8,*), i*1.d0, Real(Exm2D(i+N+1))
                        write(9,*), i*1.d0, Imag(Exm2D(i+N+1))
                    end do !!cycle for writing Fcoef i=-NwrFcof, NwrFcof
                    Pint=0
                    !!!------------------------The calculation of 1D field distribution------------------------------------------
                    do j=0,Npx+1
                        i=0; Ex(j)=0;Ex1D=0; Ez(j)=0; PEx=0; Ez1D=0; ChargeOn2D=0;
                        Ro1D=0; 
                        do i=-N,N
                            i1=i+N+1
                            Ex1D=Ex1D+Exm1D(i1)*sinIJ(i1,j+1)
                            Ex(j)=Ex(j)+Exm2D(i1)*sinIJ(i1,j+1)
                            Ez1D=Ez1D+Ezm1D(i1)*sinIJ(i1,j+1)
                            Ez(j)=Ez(j)+Ezm2D(i1)*sinIJ(i1,j+1)  !!cos(2*Pi*i*(j/((Npx+1.0))-fill/2))+c1*sin(2*Pi*i*(j/((Npx+1.0))-fill/2)) ! *sinIJ(i1,j+1)
                            ChargeOn2D=ChargeOn2D+ChargeOn2Dm(i1)*sinIJ(i1,j+1)
                        end do !var i=-N,N
                        x=2*j/((Npx+1)*1.0)-1.d0
                        do i=1,Neq
                            results=0
                            NChb=i-1
                            m1=floor(Nchb/2.0)
                                do m2=0,m1
                                    results=results+BINOM(Nchb+1,2*m2+1)*(x*x-1)**m2*x**(Nchb-2*m2)
                                end do
                            PEx=PEx+Cn(i)*results !
                            Ro1D=Ro1D+Cn(i)*(chebNder(x,i)*sqrt(1-x*x)-results*x/(sqrt(1-x*x)+0.01))
                        end do ! end cycle for i=1,Neq
                        Jx=PEx*sig1D*sqrt(1.0d0-x*x)/(Pi*(1-2*Pi*v*c1*tgrat));    
                        Px=8*Pi*0.5*sig1D*(Real(jx)*Real(PEx)+Imag(jx)*Imag(PEx))/cl        
                        Ro1D=-(2*Pi)*2*c1*sig1D*Ro1D/(2*Pi*v*1.E12*(1-2*Pi*v*c1*tgrat)*wds*1.d-4)
                        x=j/((Npx+1)*1.0)*fill
                        write(2,*) x,Real(Ro1D)
                        write(3,*) x,Imag(Ro1D)
                        write(6,*) x,Real(PEx);
                        write(10,*) x,Imag(PEx);
                        write(11,*) x,Real(Jx);
                        write(12,*) x,Imag(Jx);
                        x=j/((Npx+1)*1.0)
                        write(4,*) x,Real(Ex(j));
                        write(5,*) x,Imag(Ex(j));
                        write(13,*) x,Real(Ez(j))
                        write(14,*) x,Imag(Ez(j))
                        Pint=Pint+4*(Px*fill)/(Npx+1)/Pi
                    end do !j=0,Npx+1
                    !======== end cycle for J (max - Npx+1)
                        !    do i=-Nfcmax,Nfcmax
                        !        i1=i+N+1 
                                
                                !write(4,9),i*1.d0,Abs(Exm2D(i1));    !!!!!!!!!!!!!!write(7,9),i*1.d0,Real(Exm2DF(i1)),Imag(Exm2DF(i1)),Abs(Exm2DF(i1))
                                !write(8,*),i*1.d0,P2Dm(i1)
                        !    end do ! cycle for writing fourier coef
                    !!!!==============2D field distrubution=============================
                    do Time=0,Tmax
                        tPhaze=(2*Pi*Time)/Tmax
                        print *,tPhaze
                        write(strTime1,12) time*1.d0
                        strTime=strTime1(2:3) !//strTime1(4:4)
                        open(15,file=Dirname//'Wxz1\WxzF05'//trim(adjustl(Vstr1))//'T'//trim(adjustl(strTime))//'.dat')
                        write(15,19) 0.0, Xpts
                        Lcheck=0.d0
                        do jz=0,Npz+1
                            z=zmin+jz*(zmax-zmin)/(Npz*1.d0)
                            do i=-N,N
                                i1=i+N+1
                                if (i==0) then
                                    AmCf(i)=Exm1D(N+1)-1.D0
                                    BmCf(i)=0.5D0*(1+Rm0)*Exm1D(N+1)
                                    CmCf(i)=0.5D0*(1-Rm0)*Exm1D(N+1)
                                    HmCf(i)=HmCalc(i,v,-c1,-c1)
                                    FmCf(i)=FmCalc(i,v,-c1,-c1)
                                    Exm2D(i1)=(cos(2*Pi*v*z*sqrt(eps2)/lambda)+c1*Rm0*sin(2*Pi*v*z*sqrt(eps2)/lambda))*Exm1D(N+1) 
                                    f2=2*Pi*v*z*sqrt(eps2)/lambda
                                    !Exm2D(i1)=BmCf(i)*exp(c1*f2)+CmCf(i)*exp(-c1*f2)
                                    Ezm1D(i1)=0;    Ezm2D(i1)=0;    ChargeOn2Dm(i1)=0;
                                else
                                    AmCf(i)=Exm1D(i1)
                                    BmCf(i)=0.5D0*(1+Rm(i1))*Exm1D(i1)
                                    CmCf(i)=0.5D0*(1-Rm(i1))*Exm1D(i1)
                                    HmCf(i)=HmCalc(i,v,kappa2(abs(i)),KappaS(abs(i)))
                                    FmCf(i)=FmCalc(i,v,kappa2(abs(i)),KappaS(abs(i)))
                                    if (z<0) then
                                        f2=2*Pi*Kapa0(abs(i))*sqrt(eps0)*v*z/lambda 
                                        if (abs(f2)<15.d0) then
                                            Exm2D(i1)=AmCf(i)*exp(f2)
                                            Ezm2D(i1)=-c1*i/sqrt(i*i*1.0-v*v*eps0/(par*par))*Exm2D(i1)
                                        else
                                            Exm2D(i1)=0.d0;    Ezm2D(i1)=0.d0
                                        end if     ! f2<10                                                  
                                    endif !(z<0)
                                    if ((z>0) .and. (z<d) ) then
                                        f2=2*Pi*Kappa2(abs(i))*sqrt(epsS)*v*z/lambda
                                        if (abs(f2)<15.d0) then
                                            Exm2D(i1)=BmCf(i)*exp(-f2)+CmCf(i)*exp(f2)
                                            Ezm2D(i1)=c1*i/sqrt(i*i*1.0-v*v*eps2/(par*par))*(BmCf(i)*exp(-f2)-CmCf(i)*exp(f2))
                                        else
                                            Exm2D(i1)=0.d0;    Ezm2D(i1)=0.d0
                                        end if     ! f2<10                                                  
                                    endif !((z>0) .and. (z<d) )
                                    if ((z>d) .and. (z<Ds1)) then
                                        f2=2*Pi*KappaS(abs(i))*sqrt(epsS)*v*(z-d)/lambda
                                        if (abs(f2)<15.d0) then
                                            Exm2D(i1)=HmCf(i)*exp(-f2)+FmCf(i)*exp(f2)
                                            Ezm2D(i1)=c1*i/sqrt(i*i*1.0-v*v*epsS/(par*par))*(HmCf(i)*exp(-f2)-FmCf(i)*exp(f2))
                                        else
                                            Exm2D(i1)=0.d0;    Ezm2D(i1)=0.d0
                                        end if     ! f2<10                                               
                                    endif !((z>d) .and. (z<Ds1))
                                    !----OLD--FORMULATION---FOR--PART--OF-STRUCTURE---
                                    ! f2=2*Pi*Kappa2(abs(i))*sqrt(eps2)*v*z/lambda
                                    ! if (abs(f2)<15.d0) then
                                    !     Exm2D(i1)=0.5*((1+Rm(i1))*exp(-f2)+(1-Rm(i1))*exp(f2))*Exm1D(i1) !
                                    !     Ezm2D(i1)=c1*i/sqrt(i*i*1.0-v*v*eps2/(par*par))*Exm2D(i1)!(0.5*((1+Rm(i1))*exp(-f2)-(1-Rm(i1))*exp(f2))*Exm1D(i1)) !
                                    ! else
                                    !     Exm2D(i1)=0.d0;    Ezm2D(i1)=0.d0
                                    ! end if !(abs(f2)<15.d0)
                                end if ! (i==0)
                            end do    ! for var i=-N,N

                            do j=0,Npx+1
                                i=0; Ex(j)=0;Ex1D=0; Ez(j)=0; PEx=0; Ez1D=0; ChargeOn2D=0;
                                do i=-N,N
                                    i1=i+N+1
                                    Ex(j)=Ex(j)+Exm2D(i1)*sinIJ(i1,j+1)
                                    Ez(j)=Ez(j)+Ezm2D(i1)*sinIJ(i1,j+1)
                                end do !var i
                                x=j/((Npx+1)*1.0)
                                za=z/streep; AbsEx=Abs(Ex(j));    AbsEz=Abs(Ez(j))
								
                                if (time_anim==1) then
									AbsE=Real(Ex(j)*exp(-c1*tPhaze))*Real(Ex(j)*exp(-c1*tPhaze))+Real(Ez(j)*exp(-c1*tPhaze))*Real(Ez(j)*exp(-c1*tPhaze))                                              
									else
									AbsE=AbsEx**2+AbsEz**2
								endif !(time_anim=1)
                                if (Z<0) then
                                    AbsE=Real(eps0)*AbsE
                                endif
                                if ((Z>0) .and. (Z<d)) then
                                    AbsE=Real(eps2)*AbsE
                                endif
                               if ((Z>d) .and. (Z<Ds)) then
                                    AbsE=Real(epsS)*AbsE
                                endif
                                WxzOnX(j)=AbsE
                                !AbsE=AbsEx*AbsEx+AbsEz*AbsEz
                                !write(15,*) x,z,AbsE
                                            !if (LinSubstr(2,k)>0) then  
                                            !    Lcheck=Lcheck+(2*Pi*v*Ds/lambda)*Imag(EpsS)*AbsE/(LinSubstr(2,k))/((Npx+1)*Npz) 
                                            !    write(15,*) x,z,(2*Pi*v*Ds/lambda)*Imag(EpsS)*AbsE/(LinSubstr(2,k))
                                            !end if !(LinSubstr(2,k)>0)
                            end do    ! var j=0,Npx+1
                            write(15,19) z,WxzOnX   
                        end do ! end cycle for jz=0,Npz+1       
                        close(15);
                    end do ! end cycle for Time
                    print *,'Lcheck=',Lcheck
                    print *,'n=',Real(sqrt(epsS))
                    print *,'k=',Imag(sqrt(epsS))
                    !!!===========END calculation 2D field distrubution
                    !write(27,*) v,Pint 
                    close(2);    close(3);    close(4);    close(5);    close(6); close(7);
                    !close(8);            !close(9);
                    close(10);    close(11);    close(12); close(13);    close(14);
                end if
            end if !((k>4) .and. (buildField==1))
            !=======================================================
            do j=1,Npx+1
                PEx=0
                Ro1D=0
                x=2*j/((Npx+1)*1.0)-1.d0
                do i=1,Neq
                    results=0
                    NChb=i-1
                    m1=floor(Nchb/2.0)
                    do m2=0,m1
                        results=results+BINOM(Nchb+1,2*m2+1)*(x*x-1)**m2*x**(Nchb-2*m2)
                    end do ! cycle for m2
                    PEx=PEx+Cn(i)*results !
                    Ro1D=Ro1D+Cn(i)*(chebNder(x,i)*sqrt(1-x*x)-results*x/(sqrt(1-x*x)+0.01))
                end do ! end cycle for i=1,Neq
                write(29,*) (v/vgup*1.D12),x,Abs(Ro1D)
            end do !! cycle for j=1,Npx+1
        end do !!**k=1,Npoints****END cykle for frequency***************
        close(10)
        !close(27) !file for loss (calculated as sum of power)                
        close(27)    
        close(28)
        DEALLOCATE(A,B,W,Rm)
        if (N>0) then
                deallocate(kapa0,Kappa1,Kappa2,KappaS)
        end if
        !**************Save in file*******************
        if (photorespcalc==1) then
            open(1,file=Dirname//'Photoresp/P'//Fname//'.dat')
            do k=1,Npoints
                write(1,*),v1+((k-1)*(v2-v1))/Npoints,dU(k) ! saving photoresp in file
            end do
            close(1)
        end if
        Vgup=Vgup*1.d-12
        OPEN(1,file=Dirname//RFName//Fname//'.dat') !OPEN(1,file=Dirname//RFName//namV//Pref//Fname)
        do k=1,Npoints
            !write(1,*),(v1+((k-1)*(v2-v1))/Npoints)/vgup,R(k) !noramalised on Eg(forbidden zone)
            write(1,*),(v1+((k-1)*(v2-v1))/Npoints),R(k)
        end do
        CLOSE(1)

        !============SAVING LOSSES IN FILE FOR TASK WITH 
        !============ABSORBING Substrate==========(InSb, infrared)
        ! OPEN(1,file=Dirname//LFName//Fname//'Sub.dat')
        ! write(1,'(2(1x f15.7))'), LinSubstr
        ! CLOSE(1)
        ! LinGrat(1,:)=LinSubstr(1,:)
        
        ! OPEN(1,file=Dirname//LFName//Fname//'Grat.dat')
        ! write(1,'(2(1x f15.7))'), LinGrat
        ! CLOSE(1)
        
        ! LinGrat(2,:)=LinSubstr(2,:)+LinGrat(2,:)
        ! OPEN(1,file=Dirname//LFName//Fname//'TOT.dat')
        ! write(1,'(2(1x f15.7))'), LinGrat
        ! CLOSE(1)
        
        OPEN(1,file=Dirname//LFName//Fname//'.dat') !    OPEN(1,file=Dirname//LFName//namV//Pref//Fname)
            do k=1,Npoints
                !write(1,*),(v1+((k-1)*(v2-v1))/Npoints)/Vgup,1-TrMkh(k)-R(k) !normalised on Eg
                write(1,*),(v1+((k-1)*(v2-v1))/Npoints),1-TrMkh(k)-R(k)
            end do
        CLOSE(1)

        ! OPEN(1,file=Dirname//LFName//Fname//'Up.dat') !    OPEN(1,file=Dirname//LFName//namV//Pref//Fname)
        !     do k=1,Npoints
        !         write(1,*),v1+((k-1)*(v2-v1))/Npoints,1-TrUp2d(k)-R(k)
        !     end do
        ! CLOSE(1)

        ! OPEN(1,file=Dirname//'TL\'//'TLP'//Fname//'.dat')
        ! do k=1,Npoints
        ! write(1,18),v1+((k-1)*(v2-v1))/Npoints,TrMkh(k),1-TrMkh(k)-R(k) !,dU(k)
        ! end do
        ! CLOSE(1)

        open(1,file=Dirname//TFName//Fname//'.dat')    !open(1,file=Dirname//TFName//Mkh//namV//Pref//Fname)    
        do k=1,Npoints
            !write(1,*),(v1+((k-1)*(v2-v1))/Npoints)/vgup,TrMkh(k) ! frequency normalised on Eg
            write(1,*),(v1+((k-1)*(v2-v1))/Npoints),TrMkh(k)
            AllTrans(Ickl,k)=TrMkh(k)
        end do
        CLOSE(1)

        ! !writing Transmission for Upper 2D geometry
        ! OPEN(1,file=Dirname//TFName//Fname//'Up.dat') 
        !     do k=1,Npoints
        !         write(1,*),v1+((k-1)*(v2-v1))/Npoints,TrUp2d(k)
        !     end do
        ! CLOSE(1)



        if (envmax==1) then    
            open(1,file=Dirname//'env\EnT'//Fname//'.dat')    !open(1,file=Dirname//TFName//Mkh//namV//Pref//Fname)    
            open(2,file=Dirname//'env\EnL'//Fname//'.dat')    
            !open(3,file=Dirname//'envCurv10'//Fname)    
            do k=1,Npoints
                if ((k>1) .and. (k<Npoints)) then
                    if ((TrMkh(k-1)<TrMkh(k)) .and. (TrMkh(k)>TrMkh(k+1))) then
                            write(1,*),v1+((k-1)*(v2-v1))/Npoints,TrMkh(k)
                    end if

                    if  (((TrMkh(k-1)+R(k-1))>(TrMkh(k)+R(k))) .and. ((TrMkh(k)+R(k))<(TrMkh(k+1)+R(k+1)))) then
                        write(2,*),v1+((k-1)*(v2-v1))/Npoints,1-TrMkh(k)-R(k)
                    !tmp=TrMkh(k)
                    end if
                    !if ((TrMkh(k-1)>TrMkh(k)) .and. (TrMkh(k)<TrMkh(k+1))) then
                    !    write(2,*),v1+((k-1)*(v2-v1))/Npoints,TrMkh(k)
                    !    write(3,*),v1+((k-1)*(v2-v1))/Npoints,(TrMkh(k)+tmp)/2
                    !end if
                end if

            end do 
            CLOSE(1)
            CLOSE(2)
        end if !env max
        !CLOSE(3)
        deallocate(Fn,Cn,W1,Matr,LUfac,ipvt,wk,Cij)
        DEALLOCATE(R,R2m,Tr,TrMkh,TrUp2D,dU,P2Dm,LinSubstr,LinGrat)
        deallocate(Ex,Ez)    
        
	call cpu_time(finish_time)
        print *,'Calculation time=',finish_time-start_time
        !********************************************
    end do !!---End of external cykle----

    deallocate(AllTrans,Bjkl,Bjkl1,Det1,Det2)
    deallocate(Bm0,BslN,Bm,Sigma,sinIJ) 
    deallocate(Slk,Wbeta,Exm2D,Exm1D,Ezm2D,Ezm1D,ChargeOn2Dm)
    if (orto==1) then
        deallocate(cscoef,breack)    
    end if
    !    x=ChebNder(x)
    close(29)
    !--------FORMAT SECTION--------
        5  format(1xE10.4)
        9  format(1x 4f15.9)
        11 format(1x 1f4.1)
        12 format(1x 1f4.1)
        15 format(1x <Npts+1>f15.9)
        18 format(1x 3f15.7,f15.10)
        19 format(1x <Npx+3>G25.9)
end program MikhField

complex(8) function SigDispMagn(G,tau1,v)
        use AllConst; use profil; use numcalc
implicit none
        real(8)::tau1,v,W,t; complex(8)::G
        t=tau1*1.D-12; w=2*Pi*v*1.d12
        Wcr=el*Hmag*1.D4/(meff*mass*cl)
        SigDispMagn=G*(1.d0-c1*w*t)/((1.d0-c1*w*t)**2+(Wcr*t)**2)
end function SigDispMagn


!    subroutine ReadCij
!                use Profil; use Allconst; use msimsl;
!            implicit none
!                integer::i
!                open(1,file=Fcijname)
!                read(1,*),Neq; write(*,*), 'dim Cij=',Neq
!                allocate(cij(Neq,Neq))
!                read(1,1) cij
!                close(1)
!     CALL DWRRRL('Matr',Neq,Neq,Cij,Neq,0,'','','')
!            1    format (<Neq>(1xe16.4))
!    end subroutine readCij

complex(8) function bettaGn(i,j1)
    use msimsl; use Allconst; use profil; use numcalc
    implicit none
            integer::i,j1,iweigh; 
            real(8)::erabs=1.d-7,erel=1.d-7,res1,res2,erest,omega,llim=-1.d0,ulim=1.d0
            external::WpolN
                pw=i
                omega=Pi*fill*j1
                iweigh=1
    CALL DQDAWO(WpolN,llim,ulim,iweigh,omega,erabs,erel,res1,erest)
                iweigh=2
    CALL DQDAWO(WpolN,llim,ulim,iweigh,omega,erabs,erel,res2,erest)
                bettaGn=(res1+c1*res2)/2.0
end function bettaGn

real(8) function WPolN(x)
            use msimsl; use Allconst; use profil;
    implicit none
            integer::j; real(8)::x,z,y; 
            z=1.d0;y=0.0;
    do j=1,pw
             y=y+cij(pw,j)*z; z=z*x 
    end do
            WPolN=y*dcsval(x,nint,breack,cscoef)/NormaPr            
end function WPolN

subroutine ReadF
                use msimsl;  use profil; use blanket
    implicit none
                real(8)::tmp,tmp1
                real(8),allocatable::xdata(:),fdata(:),xx(:,:)
                integer::i,ios
                ndata=0 
            open(2,file=FProfilname)
        do                     
            ndata=ndata+1;    read(2,*,iostat=ios)tmp,tmp1
            if (ios/=0) exit
        end do
            close(2)
            ndata=ndata-1
        allocate(xdata(ndata),xx(2,ndata),fdata(ndata));
        allocate(cscoef(4,ndata),breack(ndata))        
            open(2,file=FProfilname)
        do    i=1,Ndata                 
            read(2,*,iostat=ios)xdata(i),fdata(i)
            !xx(1,i)=i*1.d0
            !xx(2,i)=xdata(i)
        end do
    !    print 11,xx
        wds=xdata(ndata)-xdata(1)
        xdata=2.d0*xdata/wds
        !xx(2,:)=2.d0*xx(2,:)/wd
        print *,'wds=',wds
        !    print 11,xx
            close(2)
CALL DCSINT(Ndata,xdata,fdata,breack,cscoef)
            nint=Ndata-1

    !11 format(1x f18.5, 1x f18.5)
    open(2,file=Ftstname)
do i=1,Np
        tmp=i*2.d0/Np-1.d0
        write(2,*),tmp*wds/2.0,dcsval(tmp,nint, breack,cscoef) !,xdata(i),fdata(i)
end do
    close(2)
        deallocate(xdata,fdata)
end subroutine ReadF

subroutine ortoPolinom
        use profil; use msimsl; 
        use numcalc, only:Neq
implicit none
            real(8):: ScP,PolN,NormProf;    real(8)::x
            real(8),dimension(Neq)::norm;        integer::i,j,irule=2
        CALL ReadF
    do i=1,Neq
            do j=1,Neq
                    if (i==j) then 
                            cij(i,i)=1.d0
                    else
                            cij(i,j)=0.d0
                    end if
            end do
    end do
    NormaPr=NormProf()
    print *,"NormaPr=",NormaPr
    do i=1,Neq  
        do j=1,i-1
                cij(i,:)=cij(i,:)-cij(j,:)*ScP(i,j)/Norm(j)
        end do
            Norm(i)=ScP(i,i); !cij(i,:)=cij(i,:)/sqrt(norm(i))
            print *, "norm=",norm(i)
    end do
    do i=1,Neq
        cij(i,:)=cij(i,:)/sqrt(norm(i))
    end do
    write(*,*),'=============================='
    CALL DWRRRL('Matr',Neq,Neq,Cij,Neq,0,'','','')
    open(3,file=Fcijname)
    write(3,*),Neq
    write(3,4),cij
    4 format(<Neq>(1xe16.4))
    close(3)

end subroutine ortoPolinom

real(8) function NormProf
            use msimsl; use profil
    implicit none
            real(8)::x0=-1.d0,x1=1.d0,errabs=1.d-9,errrel=1.d-9,errest, res
            integer::irule=2
            external::waga
    CALL DQDAG(waga,x0,x1,errabs,errrel,irule,res,errest)
     NormProf=res/2.d0
end function

real(8) function Waga(x)
            use msimsl; use profil; 
    implicit none
            real(8)::x
            waga=dcsval(x,nint,breack,cscoef) !*(cos(Pi*x/2.d0))**(1.d0/6.d0)  /2.d0 !*dcsval(x,nint,breack,cscoef) !2*sqrt(1-x*x)/Pi *Pi*cos(Pi*x/2.0)/4.d0 
end function waga

real(8) function ScP(i,j)
            use msimsl; use profil
    implicit none
            real(8)::x0=-1.d0,x1=1.d0,errabs=1.d-9,errrel=1.d-9,errest, res
            integer::i,j,irule=2
            external::gamaXn
    ip=i; ip1=j;
    CALL DQDAG(gamaXn,x0,x1,errabs,errrel,irule,res,errest)
     ScP=res/2.0
end function

real(8) function gamaXn(x)
            use msimsl; use profil; 
    implicit none
            real(8)::x,PolN
            gamaXn=PolN(x,ip)*PolN(x,ip1)*dcsval(x,nint,breack,cscoef)/NormaPr !*(cos(Pi*x/2.d0))**(1.d0/6.d0)  /2.d0 !*dcsval(x,nint,breack,cscoef) !2*sqrt(1-x*x)/Pi *Pi*cos(Pi*x/2.0)/4.d0 
end function gamaXn

real(8) function PolN(x,p)
            use profil
    implicit none
    real(8)::x,z,y;  integer(8)::j,p
            z=1.d0;y=0.0;
    do j=1,p
             y=y+cij(p,j)*z; z=z*x !z=z*x; x**(j-1)
    end do
            PolN=y            
end function polN
    
subroutine weight
        use profil; use msimsl
 implicit none
        integer::i,irule=2
        real(8)::x,x0=-1.d0,x1=1.d0,errabs=1.d-9,errrel=1.d-9,errest,res
        external::Wext
        real(8)::Wint
    call DQDAG(Wext,x0,x1,errabs,errrel,Irule,Res,errest)
                NormaPr=Res/2
            write(*,*),'An estimated error of integration is ', errest
            write(*,*),'Norma=', normaPr
    open(1,file=Fprofilname)
do i=0,Np
             x=(2.d0*i)/(Np)-1.d0
             write(1,*), x,Wint(x)/normaPr  !sqrt(1-x*x) cos(Pi*x/2.d0) !(cos(Pi*x/2.d0))**(1.d0/6.d0)
end do 
    close(1)
end subroutine weight

real(8) function Wext(y)
use profil;
implicit none
real(8)::Wint,y
 Wext=Wint(y)
end function Wext

    real(8) function Wint(y)
        use profil;    use numcalc !only::PI
            implicit none
            real(8)::y
            Wint=sqrt(1-y*y) !1.d0 !! !!/2.d0 !(cos(Pi*y/2.d0))**(1.d0/6.d0) ! ! 1.d0/2.0 ! (cos(Pi*y/2.d0))**(1.d0/6.d0) !1.d0 !sqrt(1-y*y) !cos(Pi*y/2)
    end function Wint

real(8) function BesU(x1)
use msimsl;    use numcalc
real(8)::x1
real(8),dimension(Neq+1)::y1
CALL DBSINS(x1,Neq+1,y1)
BesU=0
if (fill==0) then
        print *,'fill=0'
else
        if ((dexp(2*x1/fill)-1)==0) then 
            print *,' exp(2*x1/fill)-1=0'
        else
            BesU=y1(k+1)*y1(l+1)/(x1*(exp(2*x1/fill)-1))
        end if
end if
end function BesU

complex(8) function alpha0(v)
    use msimsl;    use numcalc;use AllConst;use blanket
implicit none
  ! Variables
        real(8)::v;     complex(8)::betaA,ff,ksi,hi

  !------------The Body alpha0
  betaA=sqrt(eps0/eps1)+kappa0up
  ff=2*Pi*sqrt(eps1)*Db*v/lambda
    ksi=cos(ff); hi=sin(ff)
    alpha0=(betaA*ksi-c1*hi)/(ksi-c1*betaA*hi)
end function alpha0

complex(8) function alphaM(v,ii)
    use msimsl;    use numcalc;use AllConst;use blanket
implicit none
  ! Variables
        real(8)::v; integer(4)::ind,ii
        complex(8)::betaA,ff,tnh
!------------The Body alphaN
ind=abs(ii)
    betaA=sqrt(eps0/eps1)*Kappa1(ind)/Kapa0(ind)+kappaMup
      ff=2*Pi*sqrt(eps1)*Db*v*Kappa1(ind)/lambda
        tnh=2.d0/(1+exp(-2*ff))-1.d0
        alphaM=(betaA+tnh)/(1+betaA*tnh)
end function alphaM

complex(8) function ReflB(v)
    use msimsl;    use numcalc; use AllConst;        use blanket
implicit none
  ! Variables
        real(8)::v;    complex(8)::dzt,betaA,ff,ksi,hi;
!------------The Body ReflB
 betaA=sqrt(eps0/eps1)
  ff=2*Pi*sqrt(eps1)*Db*v/lambda
  dzt=exp(-2*Pi*C1*sqrt(eps0)*Db*v/lambda);        ksi=cos(ff); hi=sin(ff)
        ReflB=dzt*(Etot-dzt*(ksi+c1*(betaA-kappa0up)*hi))/(ksi-c1*(betaA+kappa0up)*hi)
end function ReflB

complex(8) function psi0(v)
    use msimsl;    use numcalc; use AllConst;use blanket
implicit none
  ! Variables
        real(8)::v,ff,ksi,hi;    complex(8)::dzt,betaA
!------------The Body psi0
 betaA=sqrt(eps0/eps1)
  ff=2*Pi*sqrt(eps1)*Db*v/lambda;  
  ksi=cos(ff); hi=sin(ff); dzt=exp(-2*Pi*C1*sqrt(eps0)*Db*v/lambda);        
        psi0=betaA*dzt/(ksi-c1*hi*(betaA+kappa0up))
end function psi0


subroutine mutau()
use msimsl
use numcalc
use allconst
implicit none
integer::nint,ios,i
real(8)::tmp1,tmp
real(8),allocatable::xdata(:),fdata(:),breack(:)
real(8),allocatable::cscoef(:,:)
    ndata=0 
open(2,file=VAfname)
    do                     
        ndata=ndata+1
        read(2,*,iostat=ios)tmp,tmp1
        if (ios/=0) exit
    end do
close(2)
allocate(xdata(ndata),fdata(ndata),breack(ndata))
allocate(cscoef(4,ndata))    
open(2,file=VAfname)
i=0
    do                     
    i=i+1
        read(2,*,iostat=ios)xdata(i),fdata(i)
        if (ios/=0) exit
    end do
close(2)
CALL DCSINT(Ndata,xdata,fdata,breack,cscoef)
nint=Ndata-1
mu=4.3d7*dcsder(1,Edr,nint,breack,cscoef)
print *, 'mu=',mu
Vdr=4.3d7*dcsval(Edr,nint,breack,cscoef)
tau=(mCI*meff)*mu*(1.d-4)/elCI
deallocate(xdata,fdata,breack,cscoef)
end subroutine

  real function ChebNder(x,Nord)
    ! function evalute value of derivative of Chebishev polinomials N order
    !Variables
    use msimsl
    use numcalc
    implicit none
        real(8)::x,res,y
        integer::i,p,Nord
        integer(4),dimension(Neq,Neq)::Cpi
        
    !Body of function
    if (Neq>2) then
        do p=1,Neq
            Cpi(p,:)=0;
        end do
          Cpi(1,1)=1.0; 
          Cpi(2,2)=2.0; 
        do p=3,Neq 
            do i=1,p-1
                Cpi(p,i)=Cpi(p,i)-Cpi(p-2,i)
                Cpi(p,i+1)=Cpi(p,i+1)+2*Cpi(p-1,i)
            end do
        end do
        p=Nord
        res=0
        y=1
            do i=2,Neq
                res=res+(i-1)*Cpi(p,i)*y
                y=y*x
            end do

            ChebNder=res 
    else
     if (Neq==1) then    
             ChebNder=0.0
      end if
      if (Neq==2) then
             ChebNder=2.0
      end if
    end if
 end function
