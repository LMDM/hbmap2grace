module colors
implicit none
character       ::  red*5
character       ::  green*5
character       ::  normal*5
character       ::  underline*4

contains
    subroutine setcolors
        red=char(27)//"[31m" 
        green=char(27)//"[32m"
        normal=char(27)//"[0m" 
        underline=char(27)//"[4m"
    end subroutine 
end module colors

module variables
character*30    ::  in, out , hlog , hgraph
character(len=10)               ::  line*20
character(len=10),allocatable   ::  donor(:)
character(len=10),allocatable   ::  hydrogen(:)
character(len=10),allocatable   ::  acceptor(:)
character(len=200000)           ::  btype       !this is the max size
real                            ::  R_nframes
real                            ::  R_n1
real,allocatable                ::  prev(:)     !prevalece
integer,allocatable             ::  n1(:),n2(:),n3(:)
end module variables

program hbmap2grace
use variables
use colors

call setcolors
call credits
call readcmdline

open(1,file=in)
open(2,file=hlog)
open(3,file=out)
open(34,file=hgraph)

!some initialization
do
    read(1,"(a1)") line(1:1)
    if (line(1:1)=="s") exit  !will read the 's' from "static char"
enddo

read(1,"(a)") line
read(line(2:20),*) nframes,nlines
do
    read(1,"(a10)") line(1:10)
    if(line(1:10)=='/* y-axis:') exit
enddo

do
    read(1,"(a10)") line
    if(line(1:1)=='"') then
        backspace(1)
        exit
    endif
enddo

!The actual loop
allocate(n1(nlines))    ; n1=0
allocate(n2(nlines))    ; n2=0
allocate(n3(nlines))    ; n3=0
allocate(prev(nlines))  ; prev=0
R_nframes=nframes    !convert from integer to real

! do i=1,nlines
do i=nlines,1,-1  ! here's the trick, hbmap.xpm is inverted in relation to hbond.log
    read(1,"(a)") btype(1:nframes+2)
        do k=1,nframes+2
            if(btype(k:k).eq.'o') n1(i)=n1(i)+1
            if(btype(k:k).eq.'-') n2(i)=n2(i)+1
            if(btype(k:k).eq.'*') n3(i)=n3(i)+1
        enddo
    R_n1=n1(i)           !convert from integer to real
    prev(i)=(R_n1/R_nframes)*100.0
enddo

allocate(donor(nlines))
allocate(hydrogen(nlines))
allocate(acceptor(nlines))

read(2,*) !jump comment line
do i=1,nlines
    read(2,*) donor(i),hydrogen(i),acceptor(i)
enddo

!print the result
!write(*,*) "bond -  donor  -   hydrogen  -   acceptor  - count  -  *   -  *   -  percentage "
!write(3,*) "bond -  donor  -   hydrogen  -   acceptor  - count  -  *   -  *   -  percentage "
write(*,*) "bond -  donor  -   acceptor  - count   -  percentage - #ins - #ins+pres"
write(3,*) "bond -  donor  -   acceptor  - count   -  percentage - #ins - #ins+pres"

do i=1,nlines
!write(*,1000) i,donor(i),hydrogen(i),acceptor(i),n1(i),n2(i),n3(i),prev(i)
!write(3,1000) i,donor(i),hydrogen(i),acceptor(i),n1(i),n2(i),n3(i),prev(i)
write(*,1000) i,donor(i),acceptor(i),n1(i),prev(i),n2(i),n3(i)
write(3,1000) i,donor(i),acceptor(i),n1(i),prev(i),n2(i),n3(i)
enddo


write(*,*)
    write(*,*) "Select a bond to write or type '0' to end"
    do
    read(*,*) j
    if(j==0) exit
    call printhbonds(j,nlines,nframes)
    write(*,*) "Select another bond or type '0' to end"
    enddo

!1000    format(i4,")",1x,3a12,3i8,f8.2)
1000    format(i4,")",1x,2a12,i8,f8.2,4x,2i8)
end


subroutine printhbonds(j,nlines,nframes)
use variables,  only :  btype,donor,hydrogen,acceptor
character(len=30)       ::  label
call gotofirstline

write(*,*) "Please provide label for this bond"
read(*,*) label
write(34,"(a,a,a,3a12)") "# ",trim(label)," = ", donor(j),hydrogen(j),acceptor(j)

do i=nlines,1,-1  ! here's the trick, hbmap.xpm is inverted in relation to hbond.log
    read(1,"(a)") btype(1:nframes+2)
    if(i==j) then
        do k=1,nframes+2
            if(btype(k:k).eq.'o') write(34,*) k,label
        enddo
    endif
enddo

end

subroutine gotofirstline
character(len=10)               ::  line*20
!some initialization
rewind(1)
do
    read(1,"(a1)") line(1:1)
    if (line(1:1)=="s") exit  !will read the 's' from "static char"
enddo

read(1,"(a)") line
read(line(2:20),*) nframes,nlines
do
    read(1,"(a10)") line(1:10)
    if(line(1:10)=='/* y-axis:') exit
enddo

do
    read(1,"(a10)") line
    if(line(1:1)=='"') then
        backspace(1)
        exit
    endif
enddo
end

! help ##########################################
subroutine credits()
use colors
write(*,*) red
write(*,*) "#################################################"
write(*,*) "# Program:  hbmap2grace                         #"
write(*,*) "# Diego E.B. Gomes(1,2), Alan W. S. da Silva(2) #"
write(*,*) "# Roberto D. Lins(3), Pedro G. Pascutti(1)      #"
write(*,*) "# Thereza A. Soares(3)                          #"
write(*,*) "# 1) Universidade Federal do Rio de Janeiro     #"
write(*,*) "# 2) University of Cambridge                    #"
write(*,*) "# 3) Pacific Northwest National Laboratory      #"
write(*,*) "# +55 21 2562.6507 | diego@biof.ufrj.br         #"
write(*,*) "#################################################"
write(*,*) normal
end


subroutine readcmdline
use variables, only :   in,out,hlog,hgraph
character*30    ::  argv
integer         ::  i, iargc, n
!program specific

in=" "
out=" "
hlog=" "
hgraph=" "

n = iargc() ;  if(n.eq.0) call error(10)
    do i = 1, n
        call getarg( i, argv )
        if(argv(1:2).eq."-h") call help()
        if(argv(1:2).eq."-i") then ; call getarg(i+1,argv)
            if(argv.eq." ".or.argv(1:1).eq."-") call error(11)
            in=argv
        endif

        if(argv(1:2).eq."-o") then ; call getarg(i+1,argv)
            if(argv.eq." ".or.argv(1:1).eq."-") call error(12)
            out=argv
        endif
!program specific
        if(argv(1:2).eq."-l") then ; call getarg(i+1,argv)
            if(argv.eq." ".or.argv(1:1).eq."-") call error(13)
            hlog=argv
        endif
        if(argv(1:2).eq."-g") then ; call getarg(i+1,argv)
            if(argv.eq." ".or.argv(1:1).eq."-") call error(14)
            hgraph=argv
        endif

    end do

if(in.eq." ") call error(11)
if(out.eq." ") call error(12)
!program specific
if(hlog.eq." ") call error(13)
if(hgraph.eq." ") call error(14)
end


! help ##########################################
subroutine  help()
use colors
    write(*,*)
    write(*,*) "Usage: hbmap2grace",red," -i",normal,"hbmap.xpm",&
    green," -l",normal,"hbond.log",red," -o",normal,"hbmap.out",&
    green," -g",normal,"hbmap.xvg"
    write(*,*)
    write(*,*) red,  " -i",normal,"= .xpm Input file"
    write(*,*) green," -l",normal,"= .log Input file"
    write(*,*) red  ," -o",normal,"= table Output file"
    write(*,*) green," -g",normal,"= .xvg (graph) Output file"
    write(*,*) " -h ","= Display this help"
    write(*,*)
STOP
end


! error handling ################################
subroutine  error(n)
use colors

    write(*,*) "Program executed with the following ERROR:"
    if (n.eq.10) write(*,*) green,"ERROR: No arguments.",normal 
    if (n.eq.11) write(*,*) green,"ERROR: No .xpm input file",normal
    if (n.eq.12) write(*,*) green,"ERROR: No output file",normal
    if (n.eq.13) write(*,*) green,"ERROR: No .log input file",normal
    if (n.eq.14) write(*,*) green,"ERROR: No .xvg output file",normal
    write(*,*)
    write(*,*) "Type "//red//"hbmap2grace -h"//normal//"to review all usage options"
STOP
end
