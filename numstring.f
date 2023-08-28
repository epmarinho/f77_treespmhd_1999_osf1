      subroutine numstring(chrd,mchr,ifile)
	character chrd*4
	if(ifile.gt.0)then
		mchr=log10(float(ifile))+1
	else
		mchr=1
	endif
	itmp=ifile
	chrd=char(0)
	do ichr=1,mchr
	  idig = mod(itmp,10)
	  chrd = char(idig+48)//chrd
	  itmp = itmp/10
	end do
	return
      end
