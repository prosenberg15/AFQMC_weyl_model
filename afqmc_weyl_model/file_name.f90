!-------------------------------------------------------
!This subroutine get the filename needs to be wrote done
!-------------------------------------------------------
subroutine get_filename()
use all_param
implicit none
call createFileName(basename,'hubb_afqmc_')
call appendBaseName(basename,'se_',set)
call appendBaseName(basename,'d_',Dimen)
call appendBaseName(basename,'ns_',Nsite)
call appendBaseName(basename,'nh_',Nhop)
call appendBaseName(basename,'nl_',Nl(1))
call appendBaseName(basename,'_',Nl(2))
call appendBaseName(basename,'_',Nl(3))
call appendBaseName(basename,'kb_',3,kbound(1))
call appendBaseName(basename,'_',3,kbound(2))
call appendBaseName(basename,'_',3,kbound(3))
call appendBaseName(basename,'t_',3,dble(t1))
call appendBaseName(basename,'l_',3,lamda)
call appendBaseName(basename,'U_',3,onsitU)
call appendBaseName(basename,'N_',Ntot)
call appendBaseName(basename,'_')
call appendBaseName(basename,dtype)
call appendBaseName(basename,'_')
call appendBaseName(basename,'Ns1_',Nspin(1))
call appendBaseName(basename,'Ns2_',Nspin(2))
call appendBaseName(basename,'dt_',3,dt)
call appendBaseName(basename,'kc_',kcrn)
call appendBaseName(basename,'bg_',bgset)
call appendBaseName(basename,'pfft_',pfft)
call appendBaseName(basename,'diagm_',diagm)
call appendBaseName(basename,'thm_',Ntherm)
call appendBaseName(basename,'nme_',Nmeas)
call appendBaseName(basename,'mgs_',StepforGram)
call appendBaseName(basename,'pp_',PP)
call appendBaseName(basename,'Nlen_',Nlen)
call appendBaseName(basename,'Mstp_',meastep)
call appendBaseName(basename,'Therm_',thermstep)

call copyName(BaseName,EnergyName)
call appendBaseName(EnergyName,'_energy.dat')

call copyName(BaseName,numName)
call appendBaseName(numName,'_num.dat')

call copyName(BaseName,edgecName)
call appendBaseName(edgecName,'_edgec.dat')

call copyName(BaseName,niName)
call appendBaseName(numName,'_ni.dat')

call copyName(BaseName,ScorrName)
call appendBaseName(ScorrName,'_sisj.dat')

call copyName(BaseName,SkName)
call appendBaseName(SkName,'_sk.dat')

call copyName(BaseName,cksName)
call appendBaseName(cksName,'_cks.dat')

call copyName(BaseName,didjName)
call appendBaseName(didjName,'_didj.dat')

call copyName(BaseName,dkName)
call appendBaseName(dkName,'_dk.dat')

call copyName(BaseName,ninjName)
call appendBaseName(ninjName,'_ninj.dat')

call copyName(BaseName,bondName)
call appendBaseName(bondName,'_bond.dat')

call copyName(BaseName,nnkName)
call appendBaseName(nnkName,'_nnk.dat')

call copyName(BaseName,pairmName)
call appendBaseName(pairmName,'_pairM.dat')
end subroutine get_filename
