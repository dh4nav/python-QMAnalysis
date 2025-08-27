      program makexyz
C
C A utility program for use with the Gaussian program suite:
C    extracts cartesian coordinates, normal modes, orbitals, and NMR
C    shieldings from a Gaussian output and generates an output file
C    which can be read by XMol, RasMol, Molecule for Macintosh, etc.
C    (note that not all information may be used by these programs)
C
C The following C-shell script can be used to call the program:
C    (remove the first three characters of each line)
C=============================================
C  #!/bin/csh
C  if ( $1 == "" ) then
C    set infile = `ls -t *.out | head -1`
C    set BASE = `basename $infile .out`
C  else
C    set BASE = `basename $1 .out`
C  endif
C  if ( -f ${BASE}.out ) then
C    makexyz.exe < ${BASE}.out > ${BASE}.xyz
C  else
C    makexyz.exe < $BASE > ${BASE}.xyz
C  endif
C  page ${BASE}.xyz
C=============================================
C
C This program is copyright 1995-2018 N.J.R. van Eikema Hommes
C It is freely available for non-commercial use
C
C For more information, contact the author via email:
C     <mailto:hommes@chemie.uni-erlangen.de>
C
      parameter(maxatm=1000,maxbuf=1000000)
      implicit real*8 (A-H,O-Z)
      character*80 buffer(maxbuf)
      character*2 elms(0:109),el(maxatm)
      character*4 orb,ort(15)
      character*6 orl(9)
      character*40 title
      dimension x(maxatm),y(maxatm),z(maxatm),c(maxatm)
      dimension d(9,maxatm),e(10,10,maxatm),s(6,maxatm)
      dimension o(5)
      logical trued, shield, gauss98
      common /buf/ buffer
      data elms /'Bq','H','He','Li','Be','B','C','N','O','F','Ne','Na',
     1 'Mg','Al','Si','P','S','Cl','Ar','K','Ca','Sc','Ti','V','Cr',
     2 'Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As','Se','Br','Kr','Rb',
     3 'Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     4 'Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu',
     5 'Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re','Os',
     6 'Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac',
     7 'Th','Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No',
     8 'Lr','Rf','Ha','Sg','Ns','Hs','Mt'/
      data ort /'S   ','PX  ','PY  ','PZ  ','XY  ','XZ  ','YZ  ','XX  ',
     1 'ZZ  ','YY  ','D-2 ','D+1 ','D-1 ','D+2 ','D 0 '/
      data orl /'  S   ','  Px  ','  Py  ','  Pz  ',' Dxy  ',' Dxz  ',
     1 ' Dyz  ','Dx2-y2',' Dz2  '/
      data title /' '/
C initialize
      twth=2.0D0/3.0D0
      do i=1,maxatm
        el(i)='Xx'
        x(i)=0.0D0
        y(i)=0.0D0
        z(i)=0.0D0
        c(i)=0.0D0
        s(1,i)=0.0D0
        s(2,i)=0.0D0
        s(3,i)=0.0D0
        do ij=1,10
          do jj=1,10
            e(jj,ij,i)=0.0D0
          enddo
        enddo
      enddo
      shield=.FALSE.
      gauss98=.FALSE.
      read(*,'(1x,a80)',err=10,end=30,iostat=j)(buffer(i),i=1,maxbuf)
C following code should be used in case the value of i is not available
C     do i=1,maxbuf
C       buffer(i)='@@@'
C     enddo
C     read(*,'(1x,a80)',err=10,end=20,iostat=j)(buffer(i),i=1,maxbuf)
   10 write(*,'(''Error'',i5,'' during read'')')j
C  20 do i=1,maxbuf
C       if(buffer(i)(:4).eq.'@@@ ')goto 30
C     enddo
   30 nlines=i-1
      j=2*maxbuf
      k=j
      jc=j
      kc=j
      io=j
      iv=j
      kch=0
      do i=1,nlines
        if((buffer(i)(19:39).eq.'Standard orientation:').or.
     1     (buffer(i)(18:38).eq.'Z-Matrix orientation:').or.
     2     (buffer(i)(19:36).eq.'Input orientation:'))then
          j=i+5
          k=2*maxbuf
        endif
        if((buffer(i)(25:45).eq.'Standard orientation:').or.
     1     (buffer(i)(25:45).eq.'Z-Matrix orientation:').or.
     2     (buffer(i)(26:46).eq.'Input orientation:'))then
          j=i+5
          k=2*maxbuf
          gauss98=.TRUE.
        endif
        if(k.gt.i.and.i.ge.j.and.buffer(i)(:10).eq.'----------')k=i-1
C Set variable io to point to start of MO list
        if(io.gt.nlines.and.buffer(i)(5:21).eq.'Molecular Orbital')io=i
        if(io.gt.nlines.and.buffer(i)(5:21).eq.'Alpha Molecular O')io=i
C Set variable iv to point to start of vibrational normal modes
        if(iv.gt.nlines.and.buffer(i)(:14).eq.'Frequencies --')iv=i
C Extract an NMR shielding value
        if(buffer(i)(33:46).eq.'  Anisotropy =') then
          read(buffer(i),'(i2,18x,f11.4,15x,f11.4)')jj,s(1,jj),s(2,jj)
          read(buffer(i+4),'(15x,3f11.4)')s(3,jj),s(4,jj),s(5,jj)
          s43=s(4,jj)-s(3,jj)
          s54=s(5,jj)-s(4,jj)
          if(s54.gt.s43)then
            s(6,jj)=s(5,jj)
          else
            s(6,jj)=s(3,jj)
          endif
          shield=.TRUE.
        endif
        if(buffer(i)(37:50).eq.'  Anisotropy =') then
          read(buffer(i),'(i6,18x,f11.4,15x,f11.4)')jj,s(1,jj),s(2,jj)
          read(buffer(i+4),'(14x,3f11.4)')s(3,jj),s(4,jj),s(5,jj)
          s43=s(4,jj)-s(3,jj)
          s54=s(5,jj)-s(4,jj)
          if(s54.gt.s43)then
            s(6,jj)=s(5,jj)
          else
            s(6,jj)=s(3,jj)
          endif
          shield=.TRUE.
        endif
C find atomic charges
        if((buffer(i)(:17).eq.'Mulliken charges:').or.
     1   (buffer(i)(:36).eq.'Mulliken charges and spin densities:'))then
          if(kch.lt.2)then
            jc=i+2
            kch=1
          endif
        endif
        if(buffer(i)(:12).eq.'ESP charges:')then
          if(kch.lt.3)then
            jc=i+2
            kch=2
          endif
        endif
        if(buffer(i)(:25).eq.'Summary of Natural Popula')then
          jc=i+6
          kch=3
          if(buffer(i+4)(:10).eq.' Atom No  ')then
            kch=4
          endif
        endif
      enddo
      ii=0
      if(j.lt.maxbuf)then
        if(gauss98)then
          do i=j,k
            ii=ii+1
            read(buffer(i),'(12x,i5,17x,3f12.6)')jj,x(ii),y(ii),z(ii)
            el(ii)=elms(jj)
            if(jj.lt.0.or.jj.gt.109) ii=ii-1
          enddo
        else
          do i=j,k
            ii=ii+1
            read(buffer(i),'(12x,i5,6x,3f12.6)')jj,x(ii),y(ii),z(ii)
            el(ii)=elms(jj)
            if(jj.lt.0.or.jj.gt.109) ii=ii-1
          enddo
        endif
      else
        ii=1
      endif
      na=ii
      ii=0
C Read atomic charges
      if((kch.gt.0).and.(jc.lt.maxbuf))then
        kc=jc+na-1
        do i=jc,kc
          ii=ii+1
          if(kch.eq.4)then
            read(buffer(i),'(8x,f12.6)')c(ii)
          else if(kch.eq.3)then
            read(buffer(i),'(12x,f12.6)')c(ii)
          else
            read(buffer(i),'(10x,f12.6)')c(ii)
          endif
        enddo
        if(kch.eq.1) title='Mulliken charges'
        if(kch.eq.2) title='ESP charges'
        if(kch.eq.3) title='NPA charges'
      endif
C Output coordinates and atomic charges
      write(*,'(i4,/,a)')na,title
      if(kch.gt.0)then
        write(*,'(A2,4f12.6)')(el(i),x(i),y(i),z(i),c(i),i=1,na)
      else
        write(*,'(A2,3f12.6)')(el(i),x(i),y(i),z(i),i=1,na)
      endif
C Output normal modes
      if(iv.lt.maxbuf)then
        write(*,'(1x,/,1x,/,''NORMALMODES'',i4)')3*na-6
        do ii=1,na-2
          read(buffer(iv),'(14x,f12.4,2f23.4)')xa,xb,xc
          iv=iv+3
          if(buffer(iv)(:6).eq.'IR Int')then
            read(buffer(iv),'(14x,f12.4,2f23.4)')xia,xib,xic
          else
            xia=0.0
            xib=0.0
            xic=0.0
          endif
          iv=iv+1
          if(buffer(iv)(:6).eq.'Raman ')then
            read(buffer(iv),'(14x,f12.4,2f23.4)')xra,xrb,xrc
            iv=iv+2
            if(buffer(iv)(:6).ne.' Atom ')iv=iv+1
          else
            xra=0.0
            xrb=0.0
            xrc=0.0
          endif
          do i=1,na
C           11x for Gaussian09, 9x for Gaussian03 and older....
            read(buffer(i+iv),'(11x,3(3f7.2,2x))')(d(jj,i),jj=1,9)
          enddo
          write(*,'(a2,3f6.2)')(el(i),(d(jj,i),jj=1,3),i=1,na)
          write(*,'(3f12.4)')xa,xia,xrz
          write(*,'(a2,3f6.2)')(el(i),(d(jj,i),jj=4,6),i=1,na)
          write(*,'(3f12.4)')xb,xib,xrb
          write(*,'(a2,3f6.2)')(el(i),(d(jj,i),jj=7,9),i=1,na)
          write(*,'(3f12.4)')xc,xic,xrc
          iv=iv+na+3
        enddo
      endif
C print a list of absolute shieldings and anisotropy values
      if(shield) then
        write(*,'(1x,/,1x,/,''SHIELDINGS'')')
        do i=1,na
          if(el(i).eq.'Bq')then
            write(*,'(a2,8x,6f11.4)')el(i),(s(jj,i),jj=1,6)
          else
            write(*,'(a2,8x,2f11.4)')el(i),(s(jj,i),jj=1,2)
          endif
        enddo
C       write(*,'(a2,8x,6f11.4)')(el(i),(s(jj,i),jj=1,6),i=1,na)
      endif
C Get the 5 highest occupied and 5 lowest virtual MOs,
C combine them in a quick&dirty way to an easy to read minimal basis
      if(io.lt.maxbuf)then
        trued=.FALSE.
        i=io+4
        ii=1
   40   read(buffer(i),'(i3,i5,4x,a4,4x,5f10.5)',err=80)
     1            j,ij,orb,(o(jj),jj=1,5)
        if(j.eq.0)then
          i=i+3
          ii=1
          goto 60
        endif
        if(ij.gt.0) ii=ij
        do jj=1,15
          if(orb.eq.ort(jj))then
            ij=jj
            if(ij.gt.10)then
              ij=ij-6
              trued=.TRUE.
            endif
            goto 50
          endif
        enddo
   50   do jj=1,5
          e(ij,jj,ii)=e(ij,jj,ii)+o(jj)
        enddo
        i=i+1
        goto 40
   60   read(buffer(i),'(i3,i5,4x,a4,4x,5f10.5)',err=80)
     1            j,ij,orb,(o(jj),jj=1,5)
        if(j.eq.0) goto 80
        if(ij.gt.0) ii=ij
        do jj=1,15
          if(orb.eq.ort(jj))then
            ij=jj
            if(ij.gt.10) ij=ij-6
            goto 70
          endif
        enddo
   70   do jj=1,5
          e(ij,jj+5,ii)=e(ij,jj+5,ii)+o(jj)
        enddo
        i=i+1
        goto 60
   80   if(.not.trued)then
C         transform cartesian to true d-orbitals
          do ij=1,na
            do jj=1,10
              e(1,jj,ij)=e(1,jj,ij)+0.4427D0*
     1                             (e(8,jj,ij)+e(10,jj,ij)+e(9,jj,ij))
              e(9,jj,ij)=e(9,jj,ij)-0.5D0*(e(8,jj,ij)+e(10,jj,ij))
              e(8,jj,ij)=0.886D0*(e(8,jj,ij)-e(10,jj,ij))
            enddo
          enddo
        endif
        write(*,'(1x,/,1x,/,''MOLORBS'')')
        do ij=1,10
C         "normalize" the orbitals
          xx=0.0D0
          do i=1,na
            do jj=1,9
              xx=xx+e(jj,ij,i)**2
            enddo
          enddo
          if(xx.gt.0.0D0)then
            xx=1.0D0/sqrt(xx)
            do i=1,na
              do jj=1,9
                e(jj,ij,i)=e(jj,ij,i)*xx
              enddo
            enddo
            write(*,'(i2,2x,10a7)')ij,(orl(jj),jj=1,9)
            write(*,'(a2,1x,9f7.2)')(el(i),(e(jj,ij,i),jj=1,9),i=1,na)
          endif
        enddo
      endif
      end
