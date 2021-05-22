c J.B.Attili
c
c f77 hello.f

      program hello
      
      character(len=30) :: date,str
      write(*,*) 'Hello!'

      i=74
      
      call fdate(date)
      print *, 'Program started on ', date

      write(str,*) i
      write(*,*) trim(adjustl(str))
      write(str,'(a)') 'hey'
      write(*,*) str
      
      stop
      end 
        
