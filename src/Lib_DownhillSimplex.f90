
    module Lib_DownhillSimplex
!---^^^^^^^^^^^^^^^^^^^^^^^^^^
        use iso_fortran_env
        implicit none
        private
        
    !---
    
        public      ::      DownhillSimplex_ctor
        public      ::      delete
        public      ::      report
        
        public      ::      getN,getD
        public      ::      getBestSimplexPoint
        public      ::      addSimplexPoint
        public      ::      suggestNextPoint
        public      ::      isConverged
        
    !---
    
        logical,public              ::      DownhillSimplex_dbg = .false.
        integer,private,parameter   ::      SIMPLEX_STEP_INIT       = 0
        integer,private,parameter   ::      SIMPLEX_STEP_REFLECT    = -1
        integer,private,parameter   ::      SIMPLEX_STEP_EXPAND     = -2
        integer,private,parameter   ::      SIMPLEX_STEP_CONTRACTIN = -3
        integer,private,parameter   ::      SIMPLEX_STEP_CONTRACTOUT = -4
                   
        
        
    !---
    
        type,public     ::      DownhillSimplex
            private
            integer                                     ::      d       !   dimension of search space
            integer                                     ::      n       !   number of simplex points ( > d )
            real(kind=real64),dimension(:,:),pointer    ::      x       !   (1:d,1:n) positions of simplex points 
            real(kind=real64),dimension(:),pointer      ::      f       !   (1:n) value of function at points
            integer               ::        ibest,iworst,isecondworst
            real(kind=real64),dimension(:),pointer      ::      c       !   (1:d) centroid
            integer                                     ::      step  
            real(kind=real64)                           ::      ftmp    !   needed to store function value to compare reflect vs contract          
        end type DownhillSimplex
        
    !---
    
        interface DownhillSimplex_ctor
            module procedure    DownhillSimplex_null
            module procedure    DownhillSimplex_ctor0
        end interface
                
        interface delete
            module procedure    delete0
        end interface
        
        interface report
            module procedure    report0
        end interface
        
        interface getN
            module procedure    getN0
        end interface
        
        interface getD
            module procedure    getD0
        end interface
        
        interface addSimplexPoint
            module procedure    addSimplexPoint0
            module procedure    addSimplexPoint1
        end interface
        
    contains
!---^^^^^^^^

        function DownhillSimplex_null() result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(DownhillSimplex)           ::      this
            this%n = 0
            this%d = 0
            this%step = SIMPLEX_STEP_INIT
            nullify(this%x)
            nullify(this%f)
            nullify(this%c)
            this%ibest = 1
            this%iworst = 2
            this%isecondworst = 3
            
            return
        end function DownhillSimplex_null
                         
        function DownhillSimplex_ctor0(n,d) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(DownhillSimplex)       ::      this
            integer,intent(in)          ::      n,d
            this%n = n
            this%d = d
            this%step = SIMPLEX_STEP_INIT
            allocate(this%f(1:n))
            allocate(this%x(1:d,1:n))     
            allocate(this%c(1:d))     
            this%ibest = 1
            this%iworst = 2
            this%isecondworst = 3
            this%f = huge(1.0)
            return
        end function DownhillSimplex_ctor0
                       
    !---
    
        subroutine delete0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^
            type(DownhillSimplex),intent(inout)    ::      this
            if (this%n==0) return
            deallocate(this%f)
            deallocate(this%x)
            deallocate(this%c) 
            this = DownhillSimplex_null()
            return
        end subroutine delete0
        
    !---
    
        subroutine report0(this,u,o)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(DownhillSimplex),intent(in)    ::      this
            integer,intent(in),optional     ::      u,o
            integer     ::      uu,oo
            real(kind=real64)     ::      fbest,fworst,radius
            integer     ::      ii,kk  !,jj
            real(kind=real64)       ::      dx,dd
            real(Kind=real64),dimension(this%d) ::  cc
            uu = OUTPUT_UNIT ; if (present(u)) uu = u
            oo = ERROR_UNIT ; if (present(o)) oo = o
            fbest = minval( this%f )
            fworst = maxval( this%f )
            
            cc = 0.0d0
            do ii = 1,this%n
                cc(1:this%d) = cc(1:this%d) + this%x(1:this%d,ii)
            end do
            cc(1:this%d) = cc(1:this%d)/this%n
            
            radius = 0.0d0
            do ii = 1,this%n
                dd = 0.0d0
                do kk = 1,this%d
                    dx = this%x(kk,ii) - cc(kk)
                    dd = dd + dx*dx
                end do
                radius = max(radius,dd)
            end do
            radius = sqrt(radius)
            
            write(unit=uu,fmt='(2(a,i4),4(a,f16.6))') repeat(" ",oo)//"DownhillSimplex [d,n = ",this%d,",",this%n," , best,worst,rad ",fbest,",",fworst,",",radius,"]"
            
          !  if (radius > 1000) then
          !      do ii = 1,this%n    
          !          write(*,fmt='(i4,1000f16.4)') ii,this%f(ii),this%x(1:this%d,ii)
          !      end do
          !      print *,"centroid ",this%c
          !      stop
          !  end if
            
            
            return
        end subroutine report0
    
    !---
    
        pure function getN0(this) result(n)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(DownhillSimplex),intent(in)    ::      this
            integer                             ::      n
            n = this%n
            return
        end function getN0
        
        pure function getD0(this) result(d)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(DownhillSimplex),intent(in)    ::      this
            integer                             ::      d
            d = this%d
            return
        end function getD0
        
    
        subroutine addSimplexPoint0(this,i,x,f)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(DownhillSimplex),intent(inout)         ::      this
            integer,intent(in)                          ::      i
            real(kind=real64),dimension(:),intent(in)   ::      x
            real(kind=real64),intent(in)                ::      f
            this%x(1:this%d,i) = x(1:this%d)
            this%f(i) = f
            return
        end subroutine addSimplexPoint0
        
        subroutine findCentroidAndOrder(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      find best, worst, second worst and centroid 
            type(DownhillSimplex),intent(inout)         ::      this
            integer                 ::      ii
            real(kind=real64)       ::      fbest,fworst,fsecondworst
            
            fbest = huge(1.0)
            fworst = -huge(1.0)
            fsecondworst = 0.0
            
            this%ibest = 1
            this%iworst = 2
            this%isecondworst = 3
            
            this%c = 0
            do ii = 1,this%n
                this%c(1:this%d) = this%c(1:this%d) + this%x(1:this%d,ii)
                if (this%f(ii) < fbest) then
                    fbest = this%f(ii)
                    this%ibest = ii
                end if
                if (this%f(ii) > fworst) then   
                    fsecondworst = fworst
                    this%isecondworst = this%iworst
                    fworst = this%f(ii)
                    this%iworst = ii
                end if
            end do
            
            this%c(1:this%d) = this%c(1:this%d) / this%n    
                
            
            return
        end subroutine findCentroidAndOrder        
        
        subroutine reflect( this, x )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^        
    !*      take the worst point and reflect through centroid
            type(DownhillSimplex),intent(in)            ::      this
            real(kind=real64),dimension(:),intent(out)  ::      x
            
            x(1:this%d) = 2*this%c(1:this%d) - this%x(1:this%d,this%iworst)
            return
        end subroutine reflect
        
        subroutine expand( this, x )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^      
    !*      take the worst point and expand through centroid 
            type(DownhillSimplex),intent(in)            ::      this
            real(kind=real64),dimension(:),intent(out)  ::      x
            
            x(1:this%d) = 3*this%c(1:this%d) - 2*this%x(1:this%d,this%iworst)
            return
        end subroutine expand
        
        subroutine contractIn( this, x )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^       
    !*      take the worst point and contract towards centroid 
            type(DownhillSimplex),intent(in)            ::      this
            real(kind=real64),dimension(:),intent(out)  ::      x
            
            x(1:this%d) = ( this%c(1:this%d) + this%x(1:this%d,this%iworst) ) /2
            return
        end subroutine contractIn
        
        subroutine contractOut( this, x )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^       
    !*      take the worst point and contract towards centroid 
            type(DownhillSimplex),intent(in)            ::      this
            real(kind=real64),dimension(:),intent(out)  ::      x
            
            x(1:this%d) = ( 3*this%c(1:this%d) - this%x(1:this%d,this%iworst) ) /2
            return
        end subroutine contractOut
        
        
        subroutine shrink( this,i, x )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^        
    !*      take a point and contract towards best 
            type(DownhillSimplex),intent(in)            ::      this
            integer,intent(in)                          ::      i
            real(kind=real64),dimension(:),intent(out)  ::      x
            
            x(1:this%d) = ( this%x(1:this%d,i)  + this%x(1:this%d,this%ibest) ) /2
            return
        end subroutine shrink
        
        subroutine addSimplexPoint1(this,x,f)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(DownhillSimplex),intent(inout)         ::      this
            real(kind=real64),dimension(:),intent(out)  ::      x
            real(kind=real64),intent(in)                ::      f
            select case(this%step)
                case(SIMPLEX_STEP_REFLECT)
                    if (f > this%f(this%ibest)) then
                        !   not the best point
                        if (f < this%f(this%isecondworst)) then
                            !   and better than the second worst
                            if (DownhillSimplex_dbg) print *,"Lib_DownhillSimplex::addSimplexPoint1() info - accepting reflect ",f
                            this%x(1:this%d,this%iworst) = x(1:this%d)
                            this%f(this%iworst) = f
                            this%step = SIMPLEX_STEP_INIT
                        else if (f < this%f(this%iworst)) then
                            !   between second worst and worst. try a point slightly less far out    
                            if (DownhillSimplex_dbg) print *,"Lib_DownhillSimplex::addSimplexPoint1() info - rejecting reflect (ok) ",f
                            this%ftmp = f
                            this%step = SIMPLEX_STEP_CONTRACTOUT
                        else 
                            !   still the worst. try bringing toward centroid
                            if (DownhillSimplex_dbg) print *,"Lib_DownhillSimplex::addSimplexPoint1() info - rejecting reflect (bad) ",f
                            this%ftmp = f
                            this%step = SIMPLEX_STEP_CONTRACTIN
                        end if
                    else
                        !   new best point. try even further out
                        if (DownhillSimplex_dbg) print *,"Lib_DownhillSimplex::addSimplexPoint1() info - rejecting reflect (good) ",f
                        this%ftmp = f
                        this%step = SIMPLEX_STEP_EXPAND
                    end if
                case(SIMPLEX_STEP_EXPAND)
                    !   previously we found a point which was the best. Is this one better?
                    if (f < this%ftmp) then 
                        !   yes, even better than previous. Use this one          
                        if (DownhillSimplex_dbg) print *,"Lib_DownhillSimplex::addSimplexPoint1() info - accepting expand ",f
                        this%x(1:this%d,this%iworst) = x(1:this%d)
                        this%f(this%iworst) = f
                        this%step = SIMPLEX_STEP_INIT
                    else
                        !   no better than the previous. worth a try, but stick to reflected point.
                        if (DownhillSimplex_dbg) print *,"Lib_DownhillSimplex::addSimplexPoint1() info - rejecting expand ",f
                        call reflect( this,this%x(1:this%d,this%iworst) )
                        this%f(this%iworst) = this%ftmp
                        this%step = SIMPLEX_STEP_INIT
                    end if
                case(SIMPLEX_STEP_CONTRACTIN)
                    !   previously we found a point which was even worse than the worst!
                    if (f < this%f(this%iworst)) then 
                        !   yes, at least it is better than previous. Use this one          
                        if (DownhillSimplex_dbg) print *,"Lib_DownhillSimplex::addSimplexPoint1() info - accepting contract(in) ",f
                        this%x(1:this%d,this%iworst) = x(1:this%d)
                        this%f(this%iworst) = f
                        this%step = SIMPLEX_STEP_INIT
                    else
                        !   no better than the previous. shrink
                        if (DownhillSimplex_dbg) print *,"Lib_DownhillSimplex::addSimplexPoint1() info - rejecting contract(in) ",f
                        this%step = 1
                    end if
                case(SIMPLEX_STEP_CONTRACTOUT)
                    !   previously we found a point which wasn't much better. Is this one better?
                    if (f < this%ftmp) then 
                        !   yes, at least it is better than previous. Use this one          
                        if (DownhillSimplex_dbg) print *,"Lib_DownhillSimplex::addSimplexPoint1() info - accepting contract(out) ",f
                        this%x(1:this%d,this%iworst) = x(1:this%d)
                        this%f(this%iworst) = f
                        this%step = SIMPLEX_STEP_INIT
                    else
                        !   no better than the previous. shrink
                        if (DownhillSimplex_dbg) print *,"Lib_DownhillSimplex::addSimplexPoint1() info - rejecting contract(out) ",f
                        this%step = 1
                    end if
                case default
                    !   we must be adding a shrink step
                    if (DownhillSimplex_dbg) print *,"Lib_DownhillSimplex::addSimplexPoint1() info - accepting shrink step ",this%step,f
                    this%x(1:this%d,this%step) = x(1:this%d)
                    this%f(this%step) = f
                    this%step = this%step + 1
            end select                          
                                       
            return
        end subroutine addSimplexPoint1
        
        subroutine suggestNextPoint(this,x)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      where should we try next?
            type(DownhillSimplex),intent(inout)         ::      this
            real(kind=real64),dimension(:),intent(out)  ::      x
            
            select case(this%step)
                case(SIMPLEX_STEP_INIT)
                !   we dont know centroid or best/worst points. Compute and suggest a reflection of worst.
                    call findCentroidAndOrder(this)
                    if (DownhillSimplex_dbg) print *,"Lib_DownhillSimplex::suggestNextPoint() info - suggesting reflect ",this%f(this%iworst)
                    call reflect( this, x )
                    this%step = SIMPLEX_STEP_REFLECT
                case(SIMPLEX_STEP_EXPAND)
                !   try a reflection further out
                    if (DownhillSimplex_dbg) print *,"Lib_DownhillSimplex::suggestNextPoint() info - suggesting expand ",this%f(this%iworst),this%ftmp
                    call expand(this,x)
                case(SIMPLEX_STEP_CONTRACTOUT)
                !   try a reflection less far out
                    if (DownhillSimplex_dbg) print *,"Lib_DownhillSimplex::suggestNextPoint() info - suggesting contract out ",this%f(this%iworst),this%ftmp
                    call contractout(this,x)
                case(SIMPLEX_STEP_CONTRACTIN)
                !   try bringing towards centroid
                    if (DownhillSimplex_dbg) print *,"Lib_DownhillSimplex::suggestNextPoint() info - suggesting contract in ",this%f(this%iworst),this%ftmp
                    call contractin(this,x)
                case default
                !   we are shrinking all points towards the best.
                    if (this%step == this%ibest) then
                        !   skip the best point ( it won't move anyway )
                        this%step = this%step + 1
                    end if
                    if (this%step > this%n) then
                        !   We have finished the shrink and are actually back to square 1.
                        call findCentroidAndOrder(this)
                        if (DownhillSimplex_dbg) print *,"Lib_DownhillSimplex::suggestNextPoint() info - shrink complete, suggesting reflect ",this%f(this%iworst)
                        call reflect( this, x )
                        this%step = SIMPLEX_STEP_REFLECT
                    else
                        !   suggest next point to shrink
                        if (DownhillSimplex_dbg) print *,"Lib_DownhillSimplex::suggestNextPoint() info - suggesting shrink  ",this%step,this%f(this%step)
                        call shrink( this, this%step, x )                        
                    end if
                            
            end select
            return
        end subroutine suggestNextPoint                       
        
        function isConverged( this, xeps,feps ) result(is)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      returns true if all the x points are within a distance xeps of each other
    !*      or all the f points are within a distance feps of each other
    !*      ignore x if xeps = 0
            type(DownhillSimplex),intent(in)        ::      this
            real(kind=real64),intent(in)            ::      xeps,feps
            logical                                 ::      is
            
            integer             ::      ii,kk     !,jj
            real(kind=real64)   ::      dd,dx

            real(kind=real64)     ::      fbest,fworst
            real(kind=real64),dimension(this%d)     ::      cc
            is = .false.
            
            if (xeps > 0) then
            
                cc = 0.0d0            
                do ii = 1,this%n
                    cc(1:this%d) = cc(1:this%d) + this%x(1:this%d,ii)
                end do
                cc(1:this%d) = cc(1:this%d)/this%n
                
                do ii = 1,this%n
                    dd = 0.0d0
                    do kk = 1,this%d
                        dx = this%x(kk,ii) - cc(kk)
                        if (abs(dx)>xeps) return                                
                        dd = dd + dx*dx
                    end do
                    if (dd>xeps*xeps) return
                end do
            end if
            
            if (feps > 0) then
                fbest = minval( this%f )
                fworst = maxval( this%f )
                dd = fworst - fbest
                is = (dd <= feps)
            end if
            
            return
        end function isConverged
            
        subroutine getBestSimplexPoint(this,x,f)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(DownhillSimplex),intent(in)            ::      this    
            real(kind=real64),dimension(:),intent(out)  ::      x
            real(kind=real64),intent(out)               ::      f
            integer     ::      ibest
            ibest = minloc( this%f,dim=1 )
            x(1:this%d) = this%x(1:this%d,ibest)
            f = this%f(ibest)
            return
        end subroutine getBestSimplexPoint            
            
            
            
                         
                    
    end module Lib_DownhillSimplex
    
    
!!   gfortran -ffree-line-length-256 -Og -g Lib_DownhillSimplex.f90 -o testDownhillSimplex.exe
!
!    
!    program testDownhillSimplex
!!---^^^^^^^^^^^^^^^^^^^^^^^^
!        use iso_fortran_env
!        use Lib_DownhillSimplex
!        implicit none
!        
!        type(DownhillSimplex)           ::      this
!        
!        this = DownhillSimplex_ctor()
!        
!        call report(this)
!        
!        call delete(this)
!        
!        print *,""
!        print *,"done"
!        print *,""
!        
!    end program testDownhillSimplex
!    
!        
    
        
        
    