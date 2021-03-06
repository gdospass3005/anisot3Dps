C---------------------------------------------------------------------
C        MIXED RADIX FFT ROUTINE BASED ON TEMPERTON ,JOUR.COMP.
C        PHYS. VOL 52 NUMBER 1 OCT 1983 PAGE 1.
C        
C        VARIABLES:
C                  N -     TRANSFORM LENGTH 
C                  A -     COMPLEX INPUT ARRAY OF LENGTH N
C                  C -     ADDITIONAL COMPLEX WORK ARRAY OF LENGTH N
C                  TRIGS - COMPLEX ARRAY OF DFT EXPONENTIAL FACTORS
C                          OF LENGTH N (PREVIOUSLY GENERATED)
C                  IFAX  - INTEGER ARRAY CONTAINING FACTORIZATIONS OF N
C                          ACCORDING TO N=n1*n2....*nk
C                          WHERE NK=2,OR 3,4,5,6 
C                  NFAC  - NUMBER OF FACTORS IN N (BOTH NFAC AND IFAX 
C                          PREVIOUSLY CALCULATED IN FACTOR)
C              SKIP - STRIDE OF FFT (e.g IF=2  SKIP EVERY SECOND SAMPLE)
C                  ISIGN - FFT SIGN
C---------------------------------------------------------------------
        subroutine fft(a,c,n,trigs,ifax,nfac,iskip,isign)
         dimension a(*),c(*)
         iskip2=iskip+iskip
         call fftb(a(1),a(2),c(1),c(n+1),n,trigs,ifax,
     .            nfac,iskip2,isign)
         return
        end
        subroutine fftb(a,b,c,d,n,trigs,ifax,nfac,iskip,isign)
        complex trigs(1)
        integer ifax(1)
        dimension a(1),c(1),b(1),d(1)
        la=1
         ii=1
        do i=1,nfac
         ifac=ifax(i)
         if(ii.gt.0) call pass(a,b,c,d,
     .     trigs,ifac,la,n,iskip,1,isign)
         if(ii.lt.0) call pass(c,d,a,b,
     .     trigs,ifac,la,n,1,iskip,isign)
         ii=-ii
         la=la*ifac
        end do
        if(ii.lt.0) then
        n2=n*iskip
         jj=0
         do j=1,n2,iskip
          jj=jj+1
          a(j)=c(jj)
          b(j)=d(jj)
         end do
        end if
        return
        end
        subroutine pass(a,b,c,d,trigs,ifac,la,n,iskip,
     .     iskip1,isign)
        dimension a(1),c(1),trigs(1),b(1),d(1)
        save sin60,sin72,sin36,sq54,sin3672
        data iflag/0/
        if(iflag.eq.0) then
         iflag=1
        sin60=sin(3.14159265/3.)
        sin72=sin(3.14159265*72./180.)
        sin36=sin(3.14159265*36./180.)
        sq54=0.25*sqrt(5.)
        sin3672=sin36/sin72
        end if
        m=n/ifac
        mskip=m*iskip
        ia=0
        ib=mskip
        ja=0
        laskip=la*iskip1
        jb=laskip
        i=1
        j=1
        jump=(ifac-1)*laskip
        go to(100,200,300,400,500,600),ifac
100        write(6, *)'radix one encountered'
        stop 999
200        continue
c----------------- factor two ----------------------
         do l=1,la
          ai_real=a(i)
          ai_imag=b(i)
          iib=i+ib
          aib_real=a(iib)
          aib_imag=b(iib)
          c(j)=ai_real+aib_real
          d(j)=ai_imag+aib_imag
          jjb=j+jb
          c(jjb)=ai_real-aib_real
          d(jjb)=ai_imag-aib_imag
         i=i+iskip
         j=j+iskip1
         end do
         j=j+jump
        do k=la,m-la,la
         kk=2*k+1
         t_real=trigs(kk)
         t_imag=trigs(kk+1)
         if(isign.eq.-1) t_imag=-t_imag
         do l=1,la
          ai_real=a(i)
          ai_imag=b(i)
          iib=i+ib
          aib_real=a(iib)
          aib_imag=b(iib)
          c(j)=ai_real+aib_real
          d(j)=ai_imag+aib_imag
          f1=ai_real-aib_real
          f2=ai_imag-aib_imag
          jjb=j+jb
          c(jjb)=t_real*f1-t_imag*f2
          d(jjb)=t_imag*f1+t_real*f2
         i=i+iskip
         j=j+iskip1
         end do
         j=j+jump
         end do
        return
300        continue
c        ------ factor three ---------
        ic=ib+mskip
        jc=jb+laskip
         if(isign.gt.0) then
          do l=1,la
           ai_real=a(i)
           ai_imag=b(i)
           iib=i+ib
           aib_real=a(iib)
           aib_imag=b(iib)
           iic=i+ic
           aic_real=a(iic)
           aic_imag=b(iic)
           t1_real=aib_real+aic_real
           t1_imag=aib_imag+aic_imag
           t2_real=ai_real-0.5*t1_real
           t2_imag=ai_imag-0.5*t1_imag
           t3_real=sin60*(aib_real-aic_real)
           t3_imag=sin60*(aib_imag-aic_imag)
           c(j)=ai_real+t1_real
           d(j)=ai_imag+t1_imag
           jjb=j+jb
           c(jjb)=t2_real-t3_imag
           d(jjb)=t2_imag+t3_real
           jjc=j+jc
           c(jjc)=t2_real+t3_imag
           d(jjc)=t2_imag-t3_real
           i=i+iskip
           j=j+iskip1
          end do
          j=j+jump
         do k=la,m-la,la
          k2=k+k
          kk=k2+1
          tr1_real=trigs(kk)
          tr1_imag=trigs(kk+1)
          kk=kk+k2
          tr2_real=trigs(kk)
          tr2_imag=trigs(kk+1)
          do l=1,la
           ai_real=a(i)
           ai_imag=b(i)
           iib=i+ib
           aib_real=a(iib)
           aib_imag=b(iib)
           iic=i+ic
           aic_real=a(iic)
           aic_imag=b(iic)
           t1_real=aib_real+aic_real
           t1_imag=aib_imag+aic_imag
           t2_real=ai_real-0.5*t1_real
           t2_imag=ai_imag-0.5*t1_imag
           t3_real=sin60*(aib_real-aic_real)
           t3_imag=sin60*(aib_imag-aic_imag)
           c(j)=ai_real+t1_real
           d(j)=ai_imag+t1_imag
           x1_real=t2_real-t3_imag
           x1_imag=t2_imag+t3_real
           x2_real=t2_real+t3_imag
           x2_imag=t2_imag-t3_real
           jjb=j+jb
           c(jjb)=tr1_real*x1_real-tr1_imag*x1_imag
           d(jjb)=tr1_imag*x1_real+tr1_real*x1_imag
           jjc=j+jc
           c(jjc)=tr2_real*x2_real-tr2_imag*x2_imag
           d(jjc)=tr2_imag*x2_real+tr2_real*x2_imag
           i=i+iskip
           j=j+iskip1
          end do
          j=j+jump
         end do
        else
          do l=1,la
           ai_real=a(i)
           ai_imag=b(i)
           iib=i+ib
           aib_real=a(iib)
           aib_imag=b(iib)
           iic=i+ic
           aic_real=a(iic)
           aic_imag=b(iic)
           t1_real=aib_real+aic_real
           t1_imag=aib_imag+aic_imag
           t2_real=ai_real-0.5*t1_real
           t2_imag=ai_imag-0.5*t1_imag
           t3_real=sin60*(aib_real-aic_real)
           t3_imag=sin60*(aib_imag-aic_imag)
           c(j)=ai_real+t1_real
           d(j)=ai_imag+t1_imag
           jjb=j+jb
           c(jjb)=t2_real+t3_imag
           d(jjb)=t2_imag-t3_real
           jjc=j+jc
           c(jjc)=t2_real-t3_imag
           d(jjc)=t2_imag+t3_real
           i=i+iskip
           j=j+iskip1
          end do
          j=j+jump
         do k=la,m-la,la
          k2=k+k
          kk=k2+1
          tr1_real=trigs(kk)
          tr1_imag=-trigs(kk+1)
          kk=kk+k2
          tr2_real=trigs(kk)
          tr2_imag=-trigs(kk+1)
          do l=1,la
           ai_real=a(i)
           ai_imag=b(i)
           iib=i+ib
           aib_real=a(iib)
           aib_imag=b(iib)
           iic=i+ic
           aic_real=a(iic)
           aic_imag=b(iic)
           t1_real=aib_real+aic_real
           t1_imag=aib_imag+aic_imag
           t2_real=ai_real-0.5*t1_real
           t2_imag=ai_imag-0.5*t1_imag
           t3_real=sin60*(aib_real-aic_real)
           t3_imag=sin60*(aib_imag-aic_imag)
           c(j)=ai_real+t1_real
           d(j)=ai_imag+t1_imag
           x1_real=t2_real+t3_imag
           x1_imag=t2_imag-t3_real
           x2_real=t2_real-t3_imag
           x2_imag=t2_imag+t3_real
           jjb=j+jb
           c(jjb)=tr1_real*x1_real-tr1_imag*x1_imag
           d(jjb)=tr1_imag*x1_real+tr1_real*x1_imag
           jjc=j+jc
           c(jjc)=tr2_real*x2_real-tr2_imag*x2_imag
           d(jjc)=tr2_imag*x2_real+tr2_real*x2_imag
           i=i+iskip
           j=j+iskip1
          end do
          j=j+jump
         end do
        endif
        return
400        continue
c ------------  factor 4 ---------------------------
        ic=ib+mskip
        id=ic+mskip
        jc=jb+laskip
        jd=jc+laskip
         if(isign.gt.0) then
          do l=1,la
           ai_real=a(i)
           ai_imag=b(i)
           iib=i+ib
           aib_real=a(iib)
           aib_imag=b(iib)
           iic=i+ic
           aic_real=a(iic)
           aic_imag=b(iic)
           iid=i+id
           aid_real=a(iid)
           aid_imag=b(iid)
           t1_real=ai_real+aic_real
           t1_imag=ai_imag+aic_imag
           t2_real=aib_real+aid_real
           t2_imag=aib_imag+aid_imag
           t3_real=ai_real-aic_real
           t3_imag=ai_imag-aic_imag
           t4_real=aib_real-aid_real
           t4_imag=aib_imag-aid_imag
           c(j)=t1_real+t2_real
           d(j)=t1_imag+t2_imag
           jjb=j+jb
           c(jjb)=t3_real-t4_imag
           d(jjb)=t3_imag+t4_real
           jjc=j+jc
           c(jjc)=t1_real-t2_real
           d(jjc)=t1_imag-t2_imag
           jjd=j+jd
           c(jjd)=t3_real+t4_imag
           d(jjd)=t3_imag-t4_real
           i=i+iskip
           j=j+iskip1
          end do
          j=j+jump
         do k=la,m-la,la
          k2=k+k
          kk=k2+1
          tr1_real=trigs(kk)
          tr1_imag=trigs(kk+1)
          kk=kk+k2
          tr2_real=trigs(kk)
          tr2_imag=trigs(kk+1)
          kk=kk+k2
          tr3_real=trigs(kk)
          tr3_imag=trigs(kk+1)
          do l=1,la
           ai_real=a(i)
           ai_imag=b(i)
           iib=i+ib
           aib_real=a(iib)
           aib_imag=b(iib)
           iic=i+ic
           aic_real=a(iic)
           aic_imag=b(iic)
           iid=i+id
           aid_real=a(iid)
           aid_imag=b(iid)
           t1_real=ai_real+aic_real
           t1_imag=ai_imag+aic_imag
           t2_real=aib_real+aid_real
           t2_imag=aib_imag+aid_imag
           t3_real=ai_real-aic_real
           t3_imag=ai_imag-aic_imag
           t4_real=aib_real-aid_real
           t4_imag=aib_imag-aid_imag
           c(j)=t1_real+t2_real
           d(j)=t1_imag+t2_imag
           x1_real=t3_real-t4_imag
           x1_imag=t3_imag+t4_real
           x2_real=t1_real-t2_real
           x2_imag=t1_imag-t2_imag
           x3_real=t3_real+t4_imag
           x3_imag=t3_imag-t4_real
           jjb=j+jb
           c(jjb)=tr1_real*x1_real-tr1_imag*x1_imag
           d(jjb)=tr1_imag*x1_real+tr1_real*x1_imag
           jjc=j+jc
           c(jjc)=tr2_real*x2_real-tr2_imag*x2_imag
           d(jjc)=tr2_imag*x2_real+tr2_real*x2_imag
           jjd=j+jd
           c(jjd)=tr3_real*x3_real-tr3_imag*x3_imag
           d(jjd)=tr3_imag*x3_real+tr3_real*x3_imag
           i=i+iskip
           j=j+iskip1
          end do
          j=j+jump
         end do
        else
          do l=1,la
           ai_real=a(i)
           ai_imag=b(i)
           iib=i+ib
           aib_real=a(iib)
           aib_imag=b(iib)
           iic=i+ic
           aic_real=a(iic)
           aic_imag=b(iic)
           iid=i+id
           aid_real=a(iid)
           aid_imag=b(iid)
           t1_real=ai_real+aic_real
           t1_imag=ai_imag+aic_imag
           t2_real=aib_real+aid_real
           t2_imag=aib_imag+aid_imag
           t3_real=ai_real-aic_real
           t3_imag=ai_imag-aic_imag
           t4_real=aib_real-aid_real
           t4_imag=aib_imag-aid_imag
           c(j)=t1_real+t2_real
           d(j)=t1_imag+t2_imag
           jjb=j+jb
           c(jjb)=t3_real+t4_imag
           d(jjb)=t3_imag-t4_real
           jjc=j+jc
           c(jjc)=t1_real-t2_real
           d(jjc)=t1_imag-t2_imag
           jjd=j+jd
           c(jjd)=t3_real-t4_imag
           d(jjd)=t3_imag+t4_real
           i=i+iskip
           j=j+iskip1
          end do
          j=j+jump
         do k=la,m-la,la
          k2=k+k
          kk=k2+1
          tr1_real=trigs(kk)
          tr1_imag=-trigs(kk+1)
          kk=kk+k2
          tr2_real=trigs(kk)
          tr2_imag=-trigs(kk+1)
          kk=kk+k2
          tr3_real=trigs(kk)
          tr3_imag=-trigs(kk+1)
          do l=1,la
           ai_real=a(i)
           ai_imag=b(i)
           iib=i+ib
           aib_real=a(iib)
           aib_imag=b(iib)
           iic=i+ic
           aic_real=a(iic)
           aic_imag=b(iic)
           iid=i+id
           aid_real=a(iid)
           aid_imag=b(iid)
           t1_real=ai_real+aic_real
           t1_imag=ai_imag+aic_imag
           t2_real=aib_real+aid_real
           t2_imag=aib_imag+aid_imag
           t3_real=ai_real-aic_real
           t3_imag=ai_imag-aic_imag
           t4_real=aib_real-aid_real
           t4_imag=aib_imag-aid_imag
           c(j)=t1_real+t2_real
           d(j)=t1_imag+t2_imag
           x1_real=t3_real+t4_imag
           x1_imag=t3_imag-t4_real
           x2_real=t1_real-t2_real
           x2_imag=t1_imag-t2_imag
           x3_real=t3_real-t4_imag
           x3_imag=t3_imag+t4_real
           jjb=j+jb
           c(jjb)=tr1_real*x1_real-tr1_imag*x1_imag
           d(jjb)=tr1_imag*x1_real+tr1_real*x1_imag
           jjc=j+jc
           c(jjc)=tr2_real*x2_real-tr2_imag*x2_imag
           d(jjc)=tr2_imag*x2_real+tr2_real*x2_imag
           jjd=j+jd
           c(jjd)=tr3_real*x3_real-tr3_imag*x3_imag
           d(jjd)=tr3_imag*x3_real+tr3_real*x3_imag
           i=i+iskip
           j=j+iskip1
          end do
          j=j+jump
         end do
        endif
        return
500        continue
c ------------  factor 5 ---------------------------
        ic=ib+mskip
        id=ic+mskip
        ie=id+mskip
        jc=jb+laskip
        jd=jc+laskip
        je=jd+laskip
         if(isign.gt.0) then
          do l=1,la
           ai_real=a(i)
           ai_imag=b(i)
           iib=i+ib
           aib_real=a(iib)
           aib_imag=b(iib)
           iic=i+ic
           aic_real=a(iic)
           aic_imag=b(iic)
           iid=i+id
           aid_real=a(iid)
           aid_imag=b(iid)
           iie=i+ie
           aie_real=a(iie)
           aie_imag=b(iie)
           t1_real=aib_real+aie_real
           t1_imag=aib_imag+aie_imag
           t2_real=aic_real+aid_real
           t2_imag=aic_imag+aid_imag
           t3_real=sin72*(aib_real-aie_real)
           t3_imag=sin72*(aib_imag-aie_imag)
           t4_real=sin72*(aic_real-aid_real)
           t4_imag=sin72*(aic_imag-aid_imag)
           t5_real=t1_real+t2_real
           t5_imag=t1_imag+t2_imag
           t6_real=sq54*(t1_real-t2_real)
           t6_imag=sq54*(t1_imag-t2_imag)
           t7_real=ai_real-t5_real*0.25
           t7_imag=ai_imag-t5_imag*0.25
           t8_real=t7_real+t6_real
           t8_imag=t7_imag+t6_imag
           t9_real=t7_real-t6_real
           t9_imag=t7_imag-t6_imag
           t10_real=t3_real+sin3672*t4_real
           t10_imag=t3_imag+sin3672*t4_imag
           t11_real=sin3672*t3_real-t4_real
           t11_imag=sin3672*t3_imag-t4_imag
           c(j)=ai_real+t5_real
           d(j)=ai_imag+t5_imag
           jjb=j+jb
           c(jjb)=t8_real-t10_imag
           d(jjb)=t8_imag+t10_real
           jjc=j+jc
           c(jjc)=t9_real-t11_imag
           d(jjc)=t9_imag+t11_real
           jjd=j+jd
           c(jjd)=t9_real+t11_imag
           d(jjd)=t9_imag-t11_real
           jje=j+je
           c(jje)=t8_real+t10_imag
           d(jje)=t8_imag-t10_real
           i=i+iskip
           j=j+iskip1
          end do
          j=j+jump
         do k=la,m-la,la
          k2=k+k
          kk=k2+1
          tr1_real=trigs(kk)
          tr1_imag=trigs(kk+1)
          kk=kk+k2
          tr2_real=trigs(kk)
          tr2_imag=trigs(kk+1)
          kk=kk+k2
          tr3_real=trigs(kk)
          tr3_imag=trigs(kk+1)
          kk=kk+k2
          tr4_real=trigs(kk)
          tr4_imag=trigs(kk+1)
          do l=1,la
           ai_real=a(i)
           ai_imag=b(i)
           iib=i+ib
           aib_real=a(iib)
           aib_imag=b(iib)
           iic=i+ic
           aic_real=a(iic)
           aic_imag=b(iic)
           iid=i+id
           aid_real=a(iid)
           aid_imag=b(iid)
           iie=i+ie
           aie_real=a(iie)
           aie_imag=b(iie)
           t1_real=aib_real+aie_real
           t1_imag=aib_imag+aie_imag
           t2_real=aic_real+aid_real
           t2_imag=aic_imag+aid_imag
           t3_real=sin72*(aib_real-aie_real)
           t3_imag=sin72*(aib_imag-aie_imag)
           t4_real=sin72*(aic_real-aid_real)
           t4_imag=sin72*(aic_imag-aid_imag)
           t5_real=t1_real+t2_real
           t5_imag=t1_imag+t2_imag
           t6_real=sq54*(t1_real-t2_real)
           t6_imag=sq54*(t1_imag-t2_imag)
           t7_real=ai_real-t5_real*0.25
           t7_imag=ai_imag-t5_imag*0.25
           t8_real=t7_real+t6_real
           t8_imag=t7_imag+t6_imag
           t9_real=t7_real-t6_real
           t9_imag=t7_imag-t6_imag
           t10_real=t3_real+sin3672*t4_real
           t10_imag=t3_imag+sin3672*t4_imag
           t11_real=sin3672*t3_real-t4_real
           t11_imag=sin3672*t3_imag-t4_imag
           c(j)=ai_real+t5_real
           d(j)=ai_imag+t5_imag
           x1_real=t8_real-t10_imag
           x1_imag=t8_imag+t10_real
           x2_real=t9_real-t11_imag
           x2_imag=t9_imag+t11_real
           x3_real=t9_real+t11_imag
           x3_imag=t9_imag-t11_real
           x4_real=t8_real+t10_imag
           x4_imag=t8_imag-t10_real
           jjb=j+jb
           jjc=j+jc
           jjd=j+jd
           jje=j+je
           c(jjb)=tr1_real*x1_real-tr1_imag*x1_imag
           d(jjb)=tr1_imag*x1_real+tr1_real*x1_imag
           c(jjc)=tr2_real*x2_real-tr2_imag*x2_imag
           d(jjc)=tr2_imag*x2_real+tr2_real*x2_imag
           c(jjd)=tr3_real*x3_real-tr3_imag*x3_imag
           d(jjd)=tr3_imag*x3_real+tr3_real*x3_imag
           c(jje)=tr4_real*x4_real-tr4_imag*x4_imag
           d(jje)=tr4_imag*x4_real+tr4_real*x4_imag
           i=i+iskip
           j=j+iskip1
          end do
          j=j+jump
         end do
        else
          do l=1,la
           ai_real=a(i)
           ai_imag=b(i)
           iib=i+ib
           iic=i+ic
           iid=i+id
           iie=i+ie
           aib_real=a(iib)
           aib_imag=b(iib)
           aic_real=a(iic)
           aic_imag=b(iic)
           aid_real=a(iid)
           aid_imag=b(iid)
           aie_real=a(iie)
           aie_imag=b(iie)
           t1_real=aib_real+aie_real
           t1_imag=aib_imag+aie_imag
           t2_real=aic_real+aid_real
           t2_imag=aic_imag+aid_imag
           t3_real=sin72*(aib_real-aie_real)
           t3_imag=sin72*(aib_imag-aie_imag)
           t4_real=sin72*(aic_real-aid_real)
           t4_imag=sin72*(aic_imag-aid_imag)
           t5_real=t1_real+t2_real
           t5_imag=t1_imag+t2_imag
           t6_real=sq54*(t1_real-t2_real)
           t6_imag=sq54*(t1_imag-t2_imag)
           t7_real=ai_real-t5_real*0.25
           t7_imag=ai_imag-t5_imag*0.25
           t8_real=t7_real+t6_real
           t8_imag=t7_imag+t6_imag
           t9_real=t7_real-t6_real
           t9_imag=t7_imag-t6_imag
           t10_real=t3_real+sin3672*t4_real
           t10_imag=t3_imag+sin3672*t4_imag
           t11_real=sin3672*t3_real-t4_real
           t11_imag=sin3672*t3_imag-t4_imag
           c(j)=ai_real+t5_real
           d(j)=ai_imag+t5_imag
           jjb=j+jb
           jjc=j+jc
           jjd=j+jd
           jje=j+je
           c(jjb)=t8_real+t10_imag
           d(jjb)=t8_imag-t10_real
           c(jjc)=t9_real+t11_imag
           d(jjc)=t9_imag-t11_real
           c(jjd)=t9_real-t11_imag
           d(jjd)=t9_imag+t11_real
           c(jje)=t8_real-t10_imag
           d(jje)=t8_imag+t10_real
           i=i+iskip
           j=j+iskip1
          end do
          j=j+jump
         do k=la,m-la,la
          k2=k+k
          kk=k2+1
          tr1_real=trigs(kk)
          tr1_imag=-trigs(kk+1)
          kk=kk+k2
          tr2_real=trigs(kk)
          tr2_imag=-trigs(kk+1)
          kk=kk+k2
          tr3_real=trigs(kk)
          tr3_imag=-trigs(kk+1)
          kk=kk+k2
         tr4_real=trigs(kk)
         tr4_imag=-trigs(kk+1)
         do l=1,la
          ai_real=a(i)
          ai_imag=b(i)
          iib=i+ib
          iic=i+ic
          iid=i+id
          iie=i+ie
          aib_real=a(iib)
          aib_imag=b(iib)
          aic_real=a(iic)
          aic_imag=b(iic)
          aid_real=a(iid)
          aid_imag=b(iid)
          aie_real=a(iie)
          aie_imag=b(iie)
          t1_real=aib_real+aie_real
          t1_imag=aib_imag+aie_imag
          t2_real=aic_real+aid_real
          t2_imag=aic_imag+aid_imag
          t3_real=sin72*(aib_real-aie_real)
          t3_imag=sin72*(aib_imag-aie_imag)
          t4_real=sin72*(aic_real-aid_real)
          t4_imag=sin72*(aic_imag-aid_imag)
          t5_real=t1_real+t2_real
          t5_imag=t1_imag+t2_imag
          t6_real=sq54*(t1_real-t2_real)
          t6_imag=sq54*(t1_imag-t2_imag)
          t7_real=ai_real-t5_real*0.25
          t7_imag=ai_imag-t5_imag*0.25
          t8_real=t7_real+t6_real
          t8_imag=t7_imag+t6_imag
          t9_real=t7_real-t6_real
          t9_imag=t7_imag-t6_imag
          t10_real=t3_real+sin3672*t4_real
          t10_imag=t3_imag+sin3672*t4_imag
          t11_real=sin3672*t3_real-t4_real
          t11_imag=sin3672*t3_imag-t4_imag
          c(j)=ai_real+t5_real
          d(j)=ai_imag+t5_imag
          x1_real=t8_real+t10_imag
          x1_imag=t8_imag-t10_real
          x2_real=t9_real+t11_imag
          x2_imag=t9_imag-t11_real
          x3_real=t9_real-t11_imag
          x3_imag=t9_imag+t11_real
          x4_real=t8_real-t10_imag
          x4_imag=t8_imag+t10_real
          jjb=j+jb
          jjc=j+jc
          jjd=j+jd
          jje=j+je
          c(jjb)=tr1_real*x1_real-tr1_imag*x1_imag
          d(jjb)=tr1_imag*x1_real+tr1_real*x1_imag
          c(jjc)=tr2_real*x2_real-tr2_imag*x2_imag
          d(jjc)=tr2_imag*x2_real+tr2_real*x2_imag
          c(jjd)=tr3_real*x3_real-tr3_imag*x3_imag
          d(jjd)=tr3_imag*x3_real+tr3_real*x3_imag
          c(jje)=tr4_real*x4_real-tr4_imag*x4_imag
          d(jje)=tr4_imag*x4_real+tr4_real*x4_imag
          i=i+iskip
          j=j+iskip1
         end do
         j=j+jump
        end do
       endif
       return
600       continue
c ------------  factor 6 ---------------------------
       ic=ib+mskip
       id=ic+mskip
       ie=id+mskip
       ig=ie+mskip
       jc=jb+laskip
       jd=jc+laskip
       je=jd+laskip
       jg=je+laskip
        if(isign.gt.0) then
         do l=1,la
          ai_real=a(i)
          ai_imag=b(i)
          iib=i+ib
          iic=i+ic
          iid=i+id
          iie=i+ie
          iig=i+ig
          aib_real=a(iib)
          aib_imag=b(iib)
          aic_real=a(iic)
          aic_imag=b(iic)
          aid_real=a(iid)
          aid_imag=b(iid)
          aie_real=a(iie)
          aie_imag=b(iie)
          aig_real=a(iig)
          aig_imag=b(iig)
          t1_real=aic_real+aie_real
          t1_imag=aic_imag+aie_imag
          t2_real=ai_real-0.5*t1_real
          t2_imag=ai_imag-0.5*t1_imag
          t3_real=sin60*(aic_real-aie_real)
          t3_imag=sin60*(aic_imag-aie_imag)
          y0_real=ai_real+t1_real
          y0_imag=ai_imag+t1_imag
          y4_real=t2_real-t3_imag
          y4_imag=t2_imag+t3_real
          y2_real=t2_real+t3_imag
          y2_imag=t2_imag-t3_real
          t1_real=aig_real+aib_real
          t1_imag=aig_imag+aib_imag
          t2_real=aid_real-0.5*t1_real
          t2_imag=aid_imag-0.5*t1_imag
          t3_real=sin60*(aig_real-aib_real)
          t3_imag=sin60*(aig_imag-aib_imag)
          y3_real=aid_real+t1_real
          y3_imag=aid_imag+t1_imag
          y1_real=t2_real-t3_imag
          y1_imag=t2_imag+t3_real
          y5_real=t2_real+t3_imag
          y5_imag=t2_imag-t3_real
          x0_real=y0_real+y3_real
          x0_imag=y0_imag+y3_imag
          x4_real=y4_real+y1_real
          x4_imag=y4_imag+y1_imag
          x2_real=y2_real+y5_real
          x2_imag=y2_imag+y5_imag
          x3_real=y0_real-y3_real
          x3_imag=y0_imag-y3_imag
          x1_real=y4_real-y1_real
          x1_imag=y4_imag-y1_imag
          x5_real=y2_real-y5_real
          x5_imag=y2_imag-y5_imag
          c(j)=x0_real
          d(j)=x0_imag
          jjb=j+jb
          jjc=j+jc
          jjd=j+jd
             jje=j+je
          jjg=j+jg
          c(jjb)=x1_real
          d(jjb)=x1_imag
          c(jjc)=x2_real
          d(jjc)=x2_imag
          c(jjd)=x3_real
          d(jjd)=x3_imag
          c(jje)=x4_real
          d(jje)=x4_imag
          c(jjg)=x5_real
          d(jjg)=x5_imag
          i=i+iskip
          j=j+iskip1
         end do
         j=j+jump
        do k=la,m-la,la
         k2=k+k
         kk=k2+1
         tr1_real=trigs(kk)
         tr1_imag=trigs(kk+1)
         kk=kk+k2
         tr2_real=trigs(kk)
         tr2_imag=trigs(kk+1)
         kk=kk+k2
         tr3_real=trigs(kk)
         tr3_imag=trigs(kk+1)
         kk=kk+k2
         tr4_real=trigs(kk)
         tr4_imag=trigs(kk+1)
         kk=kk+k2
         tr5_real=trigs(kk)
         tr5_imag=trigs(kk+1)
         do l=1,la
          ai_real=a(i)
          ai_imag=b(i)
          iib=i+ib
          iic=i+ic
          iid=i+id
          iie=i+ie
          iig=i+ig
          aib_real=a(iib)
          aib_imag=b(iib)
          aic_real=a(iic)
          aic_imag=b(iic)
          aid_real=a(iid)
          aid_imag=b(iid)
          aie_real=a(iie)
          aie_imag=b(iie)
          aig_real=a(iig)
          aig_imag=b(iig)
          t1_real=aic_real+aie_real
          t1_imag=aic_imag+aie_imag
          t2_real=ai_real-0.5*t1_real
          t2_imag=ai_imag-0.5*t1_imag
          t3_real=sin60*(aic_real-aie_real)
          t3_imag=sin60*(aic_imag-aie_imag)
          y0_real=ai_real+t1_real
          y0_imag=ai_imag+t1_imag
          y4_real=t2_real-t3_imag
          y4_imag=t2_imag+t3_real
          y2_real=t2_real+t3_imag
          y2_imag=t2_imag-t3_real
          t1_real=aig_real+aib_real
          t1_imag=aig_imag+aib_imag
          t2_real=aid_real-0.5*t1_real
          t2_imag=aid_imag-0.5*t1_imag
          t3_real=sin60*(aig_real-aib_real)
          t3_imag=sin60*(aig_imag-aib_imag)
          y3_real=aid_real+t1_real
          y3_imag=aid_imag+t1_imag
          y1_real=t2_real-t3_imag
          y1_imag=t2_imag+t3_real
          y5_real=t2_real+t3_imag
          y5_imag=t2_imag-t3_real
          x0_real=y0_real+y3_real
          x0_imag=y0_imag+y3_imag
          x4_real=y4_real+y1_real
          x4_imag=y4_imag+y1_imag
          x2_real=y2_real+y5_real
          x2_imag=y2_imag+y5_imag
          x3_real=y0_real-y3_real
          x3_imag=y0_imag-y3_imag
          x1_real=y4_real-y1_real
          x1_imag=y4_imag-y1_imag
          x5_real=y2_real-y5_real
          x5_imag=y2_imag-y5_imag
          c(j)=x0_real
          d(j)=x0_imag
          jjb=j+jb
          jjc=j+jc
          jjd=j+jd
             jje=j+je
          jjg=j+jg
          c(jjb)=tr1_real*x1_real-tr1_imag*x1_imag
          d(jjb)=tr1_imag*x1_real+tr1_real*x1_imag
          c(jjc)=tr2_real*x2_real-tr2_imag*x2_imag
          d(jjc)=tr2_imag*x2_real+tr2_real*x2_imag
          c(jjd)=tr3_real*x3_real-tr3_imag*x3_imag
          d(jjd)=tr3_imag*x3_real+tr3_real*x3_imag
          c(jje)=tr4_real*x4_real-tr4_imag*x4_imag
          d(jje)=tr4_imag*x4_real+tr4_real*x4_imag
          c(jjg)=tr5_real*x5_real-tr5_imag*x5_imag
          d(jjg)=tr5_imag*x5_real+tr5_real*x5_imag
          i=i+iskip
          j=j+iskip1
         end do
         j=j+jump
        end do
       else
         do l=1,la
          ai_real=a(i)
          ai_imag=b(i)
          iib=i+ib
          iic=i+ic
          iid=i+id
          iie=i+ie
          iig=i+ig
          aib_real=a(iib)
          aib_imag=b(iib)
          aic_real=a(iic)
          aic_imag=b(iic)
          aid_real=a(iid)
          aid_imag=b(iid)
          aie_real=a(iie)
          aie_imag=b(iie)
          aig_real=a(iig)
          aig_imag=b(iig)
          t1_real=aic_real+aie_real
          t1_imag=aic_imag+aie_imag
          t2_real=ai_real-0.5*t1_real
          t2_imag=ai_imag-0.5*t1_imag
          t3_real=sin60*(aic_real-aie_real)
          t3_imag=sin60*(aic_imag-aie_imag)
          y0_real=ai_real+t1_real
          y0_imag=ai_imag+t1_imag
          y4_real=t2_real+t3_imag
          y4_imag=t2_imag-t3_real
          y2_real=t2_real-t3_imag
          y2_imag=t2_imag+t3_real
          t1_real=aig_real+aib_real
          t1_imag=aig_imag+aib_imag
          t2_real=aid_real-0.5*t1_real
          t2_imag=aid_imag-0.5*t1_imag
          t3_real=sin60*(aig_real-aib_real)
          t3_imag=sin60*(aig_imag-aib_imag)
          y3_real=aid_real+t1_real
          y3_imag=aid_imag+t1_imag
          y1_real=t2_real+t3_imag
          y1_imag=t2_imag-t3_real
          y5_real=t2_real-t3_imag
          y5_imag=t2_imag+t3_real
          x0_real=y0_real+y3_real
          x0_imag=y0_imag+y3_imag
          x4_real=y4_real+y1_real
          x4_imag=y4_imag+y1_imag
          x2_real=y2_real+y5_real
          x2_imag=y2_imag+y5_imag
          x3_real=y0_real-y3_real
          x3_imag=y0_imag-y3_imag
          x1_real=y4_real-y1_real
          x1_imag=y4_imag-y1_imag
          x5_real=y2_real-y5_real
          x5_imag=y2_imag-y5_imag
          c(j)=x0_real
          d(j)=x0_imag
          jjb=j+jb
          jjc=j+jc
          jjd=j+jd
             jje=j+je
          jjg=j+jg
          c(jjb)=x1_real
          d(jjb)=x1_imag
          c(jjc)=x2_real
          d(jjc)=x2_imag
          c(jjd)=x3_real
          d(jjd)=x3_imag
          c(jje)=x4_real
          d(jje)=x4_imag
          c(jjg)=x5_real
          d(jjg)=x5_imag
          i=i+iskip
          j=j+iskip1
         end do
         j=j+jump
        do k=la,m-la,la
         k2=k+k
         kk=k2+1
         tr1_real=trigs(kk)
         tr1_imag=-trigs(kk+1)
         kk=kk+k2
         tr2_real=trigs(kk)
         tr2_imag=-trigs(kk+1)
         kk=kk+k2
         tr3_real=trigs(kk)
         tr3_imag=-trigs(kk+1)
         kk=kk+k2
         tr4_real=trigs(kk)
         tr4_imag=-trigs(kk+1)
         kk=kk+k2
         tr5_real=trigs(kk)
         tr5_imag=-trigs(kk+1)
         do l=1,la
          ai_real=a(i)
          ai_imag=b(i)
          iib=i+ib
          iic=i+ic
          iid=i+id
          iie=i+ie
          iig=i+ig
          aib_real=a(iib)
          aib_imag=b(iib)
          aic_real=a(iic)
          aic_imag=b(iic)
          aid_real=a(iid)
          aid_imag=b(iid)
          aie_real=a(iie)
          aie_imag=b(iie)
          aig_real=a(iig)
          aig_imag=b(iig)
          t1_real=aic_real+aie_real
          t1_imag=aic_imag+aie_imag
          t2_real=ai_real-0.5*t1_real
          t2_imag=ai_imag-0.5*t1_imag
          t3_real=sin60*(aic_real-aie_real)
          t3_imag=sin60*(aic_imag-aie_imag)
          y0_real=ai_real+t1_real
          y0_imag=ai_imag+t1_imag
          y4_real=t2_real+t3_imag
          y4_imag=t2_imag-t3_real
          y2_real=t2_real-t3_imag
          y2_imag=t2_imag+t3_real
          t1_real=aig_real+aib_real
          t1_imag=aig_imag+aib_imag
          t2_real=aid_real-0.5*t1_real
          t2_imag=aid_imag-0.5*t1_imag
          t3_real=sin60*(aig_real-aib_real)
          t3_imag=sin60*(aig_imag-aib_imag)
          y3_real=aid_real+t1_real
          y3_imag=aid_imag+t1_imag
          y1_real=t2_real+t3_imag
          y1_imag=t2_imag-t3_real
          y5_real=t2_real-t3_imag
          y5_imag=t2_imag+t3_real
          x0_real=y0_real+y3_real
          x0_imag=y0_imag+y3_imag
          x4_real=y4_real+y1_real
          x4_imag=y4_imag+y1_imag
          x2_real=y2_real+y5_real
          x2_imag=y2_imag+y5_imag
          x3_real=y0_real-y3_real
          x3_imag=y0_imag-y3_imag
          x1_real=y4_real-y1_real
          x1_imag=y4_imag-y1_imag
          x5_real=y2_real-y5_real
          x5_imag=y2_imag-y5_imag
          c(j)=x0_real
          d(j)=x0_imag
          jjb=j+jb
          jjc=j+jc
          jjd=j+jd
             jje=j+je
          jjg=j+jg
          c(jjb)=tr1_real*x1_real-tr1_imag*x1_imag
          d(jjb)=tr1_imag*x1_real+tr1_real*x1_imag
          c(jjc)=tr2_real*x2_real-tr2_imag*x2_imag
          d(jjc)=tr2_imag*x2_real+tr2_real*x2_imag
          c(jjd)=tr3_real*x3_real-tr3_imag*x3_imag
          d(jjd)=tr3_imag*x3_real+tr3_real*x3_imag
          c(jje)=tr4_real*x4_real-tr4_imag*x4_imag
          d(jje)=tr4_imag*x4_real+tr4_real*x4_imag
          c(jjg)=tr5_real*x5_real-tr5_imag*x5_imag
          d(jjg)=tr5_imag*x5_real+tr5_real*x5_imag
          i=i+iskip
          j=j+iskip1
         end do
         j=j+jump
        end do
       endif
       return
       end
