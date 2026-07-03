
    module Lib_LogNormal
!---^^^^^^^^^^^^^^^^^^^^
!*
!*      A simple module to work with the log-normal distribution
!*      f(x) = 1/(sqrt(2 pi) sigma x) exp[ - (ln x / mu)^2/(2 sigma^2) ]
!*      Note that
!*          int_1^\infty N^m p(N| mu,sig) dN     =   <N^m> 1/2 ( 1 + erf[  (m sig^2 + log(mu))/(sqrt(2)sig) ] )

        use iso_fortran_env
        use Lib_SafeExp
        implicit none
        private

    !---

        public          ::      LogNormal_ctor
        public          ::      report
        public          ::      delete

        public          ::      func
        public          ::      moment
        public          ::      cdf
        public          ::      variate

        public          ::      getMu,getSigma,getisqrt2sigma
!        public          ::      int1toInfty
        public          ::      momentp2
        public          ::      confidenceLevel

    !---

        real(kind=real64),parameter,private         ::      ISQRTPI = 0.56418958354775628694807945156077d0


    !---

        type,public         ::      LogNormal
            private
            real(kind=real64)               ::      sigma,mu
            real(kind=real64)               ::      isqrt2sigma                 !    = 1/sqrt(2 sigma^2)
!            real(kind=real64)               ::      Z                           !    = 1/int_1^\infty p(N| mu,sig) dN
        end type



    !---


        interface       LogNormal_ctor
            module procedure    LogNormal_null
            module procedure    LogNormal_ctor0
        end interface

        interface       report
            module procedure    report0
        end interface

        interface       delete
            module procedure    delete0
        end interface

        interface       func
            module procedure    func0
        end interface

        interface       moment
            module procedure    moment0
            module procedure    moment1
        end interface

!         interface       int1toInfty
!             module procedure    int1toInfty0
!             module procedure    int1toInfty1
!         end interface

        interface       momentp2
            module procedure    momentp20
            module procedure    momentp21
        end interface

        interface       cdf
            module procedure    cdf0
        end interface

        interface       variate
            module procedure    variate0
        end interface

        interface       getMu
            module procedure    getMu0
        end interface

        interface       getSigma
            module procedure    getSigma0
        end interface

    contains
!---^^^^^^^^

        function LogNormal_null() result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      null constructor for empty lognormal
            type(LogNormal)            ::      this
            this%sigma = 1.0d0
            this%mu = 1.0d0
            this%isqrt2sigma = 1/( sqrt(2.0d0)*this%sigma )
!             this%Z = 0.0d0
            return
        end function LogNormal_null


        function LogNormal_ctor0(mu,sigma) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      default constructor
            real(kind=real64),intent(in)                ::      sigma,mu
            type(LogNormal)                             ::      this
            this%sigma = sigma
            this%mu = mu
            if (sigma<=0.0d0) then !defense agaisnt sigma=0
                this%isqrt2sigma=0.0d0
            else
                this%isqrt2sigma = 1/( sqrt(2.0d0)*this%sigma )
            end if
!             this%Z = 2 / ( 1 + erf( log(this%mu)*this%isqrt2sigma ) )

            return
        end function LogNormal_ctor0


    !---

        subroutine report0(this,u,o)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      simple dump to unit u with indent o
    !*      default to screen
            type(LogNormal),intent(in)     ::      this
            integer,intent(in),optional                     ::      u,o
            integer             ::      uu,oo
            uu = 6 ; if (present(u)) uu = u
            oo = 0 ; if (present(o)) oo = o
            write(unit=uu,fmt='(3(a,f16.8))') repeat(" ",oo)//"LogNormal[mu,sigma = ",this%mu,",",this%sigma,"]"
            return
        end subroutine report0

        subroutine delete0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^
    !*      deallocate dynamic memory ( there isn't any )
            type(LogNormal),intent(inout)   ::      this
            this = LogNormal_ctor()
            return
        end subroutine delete0

    !---

        elemental function func0(this,x) result(f)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return the value of the function at x
            type(LogNormal),intent(in)      ::      this
            real(kind=real64),intent(in)    ::      x
            real(kind=real64)               ::      f
            if (x<=0) then
                f = 0.0d0
            else
                f = (log(x/this%mu))*this%isqrt2sigma
                f = safeExp( -f*f )
                f = ISQRTPI*this%isqrt2sigma*f/x
            end if
            return
        end function func0

         function moment0(this,m) result(mf)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return the mth moment
            type(LogNormal),intent(in)      ::      this
            integer,intent(in)              ::      m
            real(kind=real64)               ::      mf
            if (m==0) then
                mf = 1.0d0
            else if (this%sigma<20.0d0) then
                mf = (this%mu**m) * safeExp( m*m*this%sigma*this%sigma/2 )
            else
                mf = m * log( this%mu) + m*m*this%sigma*this%sigma/2
                mf = safeExp(mf)
                
            end if
            return
        end function moment0


        pure function moment1(this,m) result(mf)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return the mth moment
            type(LogNormal),intent(in)      ::      this
            real(kind=real64),intent(in)    ::      m
            real(kind=real64)               ::      mf


            if (this%mu == 0) then
                mf = 1.0d0
            else if (this%sigma<20.0d0) then
                mf = (this%mu**m) * safeExp( m*m*this%sigma*this%sigma/2 )
            else
                mf = m * log( this%mu) + m*m*this%sigma*this%sigma/2
                mf = safeExp(mf)
                
            end if

            return
        end function moment1


!         pure function int1toInfty0(this,m) result(mf)
!     !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!     !*      return the mth moment int1toInfty
!     !*      int_1^\infty N^m p(N| mu,sig) dN     =   <N^m> 1/2 ( 1 + erf[  (m sig^2 + log(mu))/(sqrt(2)sig) ] )
!             type(LogNormal),intent(in)      ::      this
!             integer,intent(in)              ::      m
!             real(kind=real64)               ::      mf
!             if (m==0) then
!                 mf = 1.0d0
!             else
!                 mf = (this%mu**m) * exp( m*m*this%sigma*this%sigma/2 )
!                 mf =  mf * this%Z * ( 1 + erf( (m*this%sigma*this%sigma + log(this%mu))*this%isqrt2sigma ) ) / 2
!             end if
!             return
!         end function int1toInfty0


!         pure function int1toInfty1(this,m) result(mf)
!     !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!     !*      return the mth int1toInfty
!             type(LogNormal),intent(in)      ::      this
!             real(kind=real64),intent(in)    ::      m
!             real(kind=real64)               ::      mf


!             if (this%mu == 0) then
!                 mf = 1.0d0
!             else
!                 mf = (this%mu**m) * exp( m*m*this%sigma*this%sigma/2 )
!                 mf =  mf * this%Z * ( 1 + erf( (m*this%sigma*this%sigma + log(this%mu))*this%isqrt2sigma ) ) / 2
!             end if

!             return
!         end function int1toInfty1



        pure function momentp20(this,m) result(mf)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return the mth moment of p(N|mu,sig)^2
    !*          p(N|mu,sig)   = 1/(sqrt(2 pi) sig N) exp[ - (ln N / mu)^2/(2 sig^2) ]
    !*          p^2(N|mu,sig) = 1/( 2 pi sig^2 N^2 ) exp[ - (ln N / mu)^2/(sig^2) ]
    !*                        = 1/( 2 sqrt( pi ) sig N ) p( N|mu,sig/sqrt(2) )
    !*      so <N^m>_{p^2}    = 1/( 2 sqrt( pi ) sig ) <N^{m-1}>


            type(LogNormal),intent(in)      ::      this
            integer,intent(in)              ::      m
            real(kind=real64)               ::      mf
            if (m==1) then
                mf = ISQRTPI / ( 2*this%sigma )
            else
                mf = ( ISQRTPI / ( 2*this%sigma ) ) * (this%mu**(m-1)) * safeExp( (m-1)*(m-1)*this%sigma*this%sigma/4 )
            end if
            return
        end function momentp20


        pure function momentp21(this,m) result(mf)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return the mth moment of p(N|mu,sig)^2
            type(LogNormal),intent(in)      ::      this
            real(kind=real64),intent(in)    ::      m
            real(kind=real64)               ::      mf


            if (m==1.0d0) then
                mf = ISQRTPI / ( 2*this%sigma )
            else
                mf = ( ISQRTPI / ( 2*this%sigma ) ) * (this%mu**(m-1)) * safeExp( (m-1)*(m-1)*this%sigma*this%sigma/4 )
            end if

            return
        end function momentp21




















        pure real(kind=real64) function confidenceLevel( this,alpha )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return the level x where cdf(x) = alpha
            type(LogNormal),intent(in)      ::      this
            real(kind=real64),intent(in)    ::      alpha

            integer                         ::      step    
            real(kind=real64),parameter     ::      EPS = 1.0d-8
            real(kind=real64)               ::      xx,ff,cc,dx

            if (alpha <= 0) then
                confidenceLevel = 0.0d0
                return
            else if (alpha >= 1) then
                confidenceLevel = huge(1.0)
                return
            end if

            xx = this%mu
            do step = 1,1000
                cc = cdf(this,xx)
                if (abs(cc-alpha)<=EPS) exit
                if (xx <= 1.0d-8) exit
                if (xx >= huge(1.0) ) exit
                ff = func(this,xx)                
                dx = (alpha-cc)/ff
                if ( xx + dx < 0.0d0 ) then
                    dx = -xx/2
                else if (xx + dx > huge(1.0)) then
                    xx = huge(1.0)
                end if
                xx = xx + dx                
            end do
            confidenceLevel = xx
            return
        end function confidenceLevel
                

















        elemental function cdf0(this,x) result(f)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      return the value of the cumulative distribution function at x
    !*      ie integral from 0 to x
            type(LogNormal),intent(in)      ::      this
            real(kind=real64),intent(in)    ::      x
            real(kind=real64)               ::      f
            if (x<=0) then
                f = 0.0d0
            else
                f = (log(x/this%mu))*this%isqrt2sigma
                f = erf( f )
                f = (1 + f)/2
            end if
            return
        end function cdf0

    !---

        function variate0(this) result(x)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      draw a log-normal distributed variate
            type(LogNormal),intent(in)      ::      this
            real(kind=real64)               ::      x

            x = normalDeviate()             !   normally distributed with zero mean and std dev 1
            x = this%sigma*x                !   now right width
            x = this%mu*safeExp( x )            !   now log normal distributed

            return
        end function variate0


    !---


        pure function getMu0(this) result(mu)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(LogNormal),intent(in)      ::      this
            real(kind=real64)               ::      mu
            mu = this%mu
            return
        end function getMu0

        pure function getSigma0(this) result(Sigma)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(LogNormal),intent(in)      ::      this
            real(kind=real64)               ::      Sigma
            Sigma = this%Sigma
            return
        end function getSigma0

        pure function getisqrt2sigma(this) result(isqrt2sigma)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(LogNormal),intent(in)      ::      this
            real(kind=real64)               ::      isqrt2sigma
            isqrt2sigma = this%isqrt2sigma
            return
        end function getisqrt2sigma


    !---


        function normalDeviate() result(x)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      compute 1 normally distributed random deviate
    !*      with mean zero and std dev 1
    !*      using the Bell & Knop method
            real(kind=real64)               ::      x

            real(kind=real64)       ::      uu,vv,ss
            do
                call random_number(uu)           !   random numbers [0:1)
                call random_number(vv)           !   random numbers [0:1)
                uu = 2*uu-1 ; vv = 2*vv-1
                ss = uu*uu + vv*vv
                if (ss*(1.0d0-ss)<=0.0d0) then
                    cycle
                else
                    ss = sqrt( -2 * log(ss)/ss )
                    x = uu*ss
                    exit
                end if
            end do
            return
        end function normalDeviate



    end module Lib_LogNormal