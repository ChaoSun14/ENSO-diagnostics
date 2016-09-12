load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
begin


    OUTPUT_TYPE = "pdf"
    OUTDIR      = "/home/reproduce/sunc/enso/plot/"
    CESMDIR     = "/home/reproduce/sunc/enso/CESM1-FASTCHEM_picontrol_atm/ts_data/"
    FGOALSDIR   = "/home/reproduce/sunc/enso/FGOALS-g2_picontrol_atm/ts_data/"
    HadISST     = "/home/reproduce/sunc/enso/obs_data/"

    wks_type = OUTPUT_TYPE
    wks_n34 = gsn_open_wks(wks_type,OUTDIR+"nino34.timeseries")  
    wks_n4 = gsn_open_wks(wks_type,OUTDIR+"nino4.timeseries")
    wks_n3 = gsn_open_wks(wks_type,OUTDIR+"nino3.timeseries")
    wks_n12 = gsn_open_wks(wks_type,OUTDIR+"nino12.timeseries")
    
  
    wks_n3_tlon_hi = gsn_open_wks(wks_type,OUTDIR+"nino3.hov.elnino")
    wks_n3_tlon_lo = gsn_open_wks(wks_type,OUTDIR+"nino3.hov.lanina")
  
    wks_n3_p = gsn_open_wks(wks_type,OUTDIR+"nino3.powspec")
    
  
    wks_n3_mst = gsn_open_wks(wks_type,OUTDIR+"nino3.monstddev")
    
    wks_n3sc = gsn_open_wks(wks_type,OUTDIR+"nino3.spatialcomp")
    
    gsn_define_colormap(wks_n34,"ncl_default")    
    gsn_define_colormap(wks_n4,"ncl_default")  
    gsn_define_colormap(wks_n3,"ncl_default")  
    gsn_define_colormap(wks_n12,"ncl_default")   
    gsn_merge_colormaps(wks_n3_tlon_hi,"BlueDarkRed18",(/"gray30","gray50","gray70"/))
    gsn_merge_colormaps(wks_n3_tlon_lo,"BlueDarkRed18",(/"gray30","gray50","gray70"/))
    gsn_define_colormap(wks_n3_p,"cb_9step")  
    gsn_define_colormap(wks_n3_mst,"ncl_default")   
    gsn_define_colormap(wks_n3sc,"BlueDarkRed18")       
    
    ;-------------------------read data-----------------------------------------------------------
    data = addfile(CESMDIR+"CESM1-FASTCHEM_piControl_r1i1p1.ts.007001-029112.nc","r")
    sst  = data->ts(0:200*12-1,:,:)
    time = data->time(0:200*12-1)
    printVarSummary(sst)

  

    nsim        = 3
    syear       = 70
    eyear       = 269
    nyr         = 200
    names       = "CESM"
    ;print(ispan(syear,eyear-1,1))

    ;----------------------------------------------------------------------------------------------
    xyn34 = new(nsim,"graphic")  
    xyn4  = new(nsim,"graphic")  
    xyn3  = new(nsim,"graphic")  
    xyn12 = new(nsim,"graphic")  
   
    plot_n34hi = new(nsim,"graphic")  
    plot_n34lo = new(nsim,"graphic")  
    
    map_n34sc_jja0 = new(nsim,"graphic")  
    map_n34sc_son0 = new(nsim,"graphic")  
    map_n34sc_djf1 = new(nsim,"graphic")  
    map_n34sc_mam1 = new(nsim,"graphic")  
  
    xyn34_rst = new(nsim,"graphic")  
    xyn34_mst = new(nsim,"graphic")
    ;-------------------------------------------------------------------------------------------
    if (isatt(sst,"is_all_missing")) then
        delete(sst)
        continue
     end if
        
    sst = where(sst.le.-1.8,-1.8,sst)
    d = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")   
    basemap = d->LSMASK
    lsm = landsea_mask(basemap,sst&lat,sst&lon)
    sst = mask(sst,conform(sst,lsm,(/1,2/)).ge.1,False)
    ;printStatInfo(sst)
    delete([/lsm,basemap/])
    delete(d)

    sst = rmMonAnnCycTLL(sst)
    

     wgt = (/1.,2.,1./)   
     wgt = wgt/sum(wgt)
     pi=4.*atan(1.0)
     rad=(pi/180.)

     coswgt=cos(rad*sst&lat)
     coswgt!0 = "lat" 
     coswgt&lat= sst&lat        
     llats = -5.     ; nino3.4
     llatn = 5.
     llonw = 190.
     llone = 240.
     nino34 = wgt_areaave_Wrap(sst(:,{llats:llatn},{llonw:llone}),coswgt({llats:llatn}),1.0,0)  
     nino34@area = llats+":"+llatn+"N, "+llonw+":"+llone+"E"
     nino34@units = sst@units
     nino34@long_name = "nino3.4 timeseries"

     llats = -5.    ; nino3
     llatn = 5.
     llonw = 210.
     llone = 270.
     nino3 = wgt_areaave(sst(:,{llats:llatn},{llonw:llone}),coswgt({llats:llatn}),1.0,0)  
     nino3@area = llats+":"+llatn+"N, "+llonw+":"+llone+"E"
     copy_VarCoords(nino34,nino3)
     nino3@units = sst@units
     nino3@long_name = "nino3 timeseries"
     
     llats = -5.    ; nino4
     llatn = 5.
     llonw = 160.
     llone = 210.
     nino4 = wgt_areaave(sst(:,{llats:llatn},{llonw:llone}),coswgt({llats:llatn}),1.0,0)  
     nino4@area = llats+":"+llatn+"N, "+llonw+":"+llone+"E"
     copy_VarCoords(nino34,nino4)
     nino4@units = sst@units
     nino4@long_name = "nino4 timeseries"
     
     llats = -10.    ; nino1+2
     llatn = 0.
     llonw = 270.
     llone = 280.
     nino12 = wgt_areaave(sst(:,{llats:llatn},{llonw:llone}),coswgt({llats:llatn}),1.0,0)  
     nino12@area = llats+":"+llatn+"N, "+llonw+":"+llone+"E"
     copy_VarCoords(nino34,nino12)
     nino12@units = sst@units
     nino12@long_name = "nino1+2 timeseries"

    ;---------------------------------------------------------------------------------------------
     nino3T = wgt_runave(nino3,wgt,1)                 ; for use in ENSO composites / hovmuellers / running standard deviations
     nino3_ndj = nino3T(11:dimsizes(nino3T)-13:12)   ; cannot count last 1yr as spatial composite uses +1yrs data betond NDJ..  
     nino3_ndj!0 = "time"
     nino3_ndj&time = ispan(syear,eyear-1,1)    
     nino3_ndj = dtrend_msg(ispan(0,dimsizes(nino3_ndj&time)-1,1),nino3_ndj,True,False)
     nino3_ndj = dim_standardize(nino3_ndj,0)

     sst = (/ dtrend_msg_n(ispan(0,nyr*12-1,1),sst,False,False,0) /) ; detrend the sst array
     
     sstr = sst(:,{-3:3},{120:280})                        ; ENSO hovmuellers based on NDJ nino3
     ; delete(sst)
     finsst_hi = sstr(:60,0,:)    ; for Jan-2 -> Jan+3
     finsst_hi!0 = "time"
     finsst_hi&time = ispan(0,60,1)
     finsst_hi = 0.
     finsst_lo = finsst_hi
     finsst_mid = finsst_hi
     cntr_hi = 0
     cntr_lo = 0
     cntr_mid = 0
     cntr_lo@_FillValue = default_fillvalue(typeof(cntr_lo)) 
     cntr_mid@_FillValue = default_fillvalue(typeof(cntr_mid)) 
     cntr_hi@_FillValue = default_fillvalue(typeof(cntr_hi)) 
        
     mocntr = 24   ; note: if this is set at 24 gg should start at 2
     do gg = 2,dimsizes(nino3_ndj)-3   ; remember that Dec is month 11. End @ -3 because we need to grab + 3 yrs and 1 month from there (nino3_ndj already ends at eyear-1)
        if (.not.ismissing(nino3_ndj(gg))) then    ; note that finsst_* indices 24:52 (Jan+0 -> May +2) are all that is shown in the hovmoller plots 
           if (nino3_ndj(gg).ge.1.) then
              finsst_hi = (/ finsst_hi+dim_avg_n(sstr(mocntr-24:mocntr+36,:,:),1) /)    ; nino3_ndj value is at sstr index mocntr+11     
               cntr_hi = cntr_hi+1
           end if
           if (nino3_ndj(gg).ge.-0.5.and.nino3_ndj(gg).le.0.5) then
              finsst_mid = (/ finsst_mid+dim_avg_n(sstr(mocntr-24:mocntr+36,:,:),1) /)          
              cntr_mid = cntr_mid+1
           end if
           if (nino3_ndj(gg).le.-1.) then
              finsst_lo = (/ finsst_lo+dim_avg_n(sstr(mocntr-24:mocntr+36,:,:),1) /)            
              cntr_lo = cntr_lo+1
           end if
        end if
        mocntr = mocntr+12
     end do
     delete([/sstr,mocntr/])

     cntr_hi  = where(cntr_hi.eq.0, cntr_hi@_FillValue, cntr_hi)
     cntr_mid = where(cntr_mid.eq.0,cntr_mid@_FillValue,cntr_mid)
     cntr_lo  = where(cntr_lo.eq.0, cntr_lo@_FillValue, cntr_lo)
     finsst_hi  = (/ finsst_hi/cntr_hi /)
     finsst_mid = (/ finsst_mid/cntr_mid /)
     finsst_lo  = (/ finsst_lo/cntr_lo /)
     delete([/coswgt/])
    ;-----------------------------------------------------------------------------------------
     iopt = 0                ; nino3.4 power spectra
     jave = (7*nyr)/100
     val1 = .95
     val2 = .99
     pct = 0.1 
     spectra_mvf = False        ; missing value flag for nino3.4
     if (any(ismissing(nino3))) then
        print("Missing data exists, power spectra function does not allow missing data")
        spectra_mvf = True
     else
        nino3_dt = dtrend_msg(ispan(0,dimsizes(nino3)-1,1),nino3,True,False) 
    
        sdof = specx_anal(nino3_dt,iopt,jave,pct)
        splt1 = specx_ci(sdof,val1,val2)
        delete([/nino3_dt,iopt,jave,pct/])
     end if
;------------------------------------------------------------------------------------------
     nino3_dt = dtrend_msg(ispan(0,dimsizes(nino3&time)-1,1),nino3,True,False)
     nino3_mon_sd = new(12,typeof(nino3))
     
     do hh = 0,11
        nino3_mon_sd(hh) = (/ dim_stddev(nino3_dt(hh::12)) /)
     end do
     delete(nino3_dt)    
;==========================================================================================
;==========================================================================================
     xyres = True
     xyres@gsnDraw = False
     xyres@gsnFrame = False
     xyres@gsnRightString = ""
     xyres@gsnLeftString = ""
     xyres@gsnYRefLine = 0.0
     xyres@gsnYRefLineColor = "gray42"
     xyres@xyLineColor = "gray62"
     if (wks_type.eq."png") then
        xyres@xyLineThicknessF = .75  
     else
        xyres@xyLineThicknessF = .5  
     end if   
     xyres@tiYAxisString = ""
     if (nsim.le.5) then
        xyres@tmXBLabelFontHeightF = 0.017
        xyres@tmYLLabelFontHeightF = 0.017
        xyres@gsnLeftStringFontHeightF = 0.024
        xyres@gsnCenterStringFontHeightF = 0.024
        xyres@gsnRightStringFontHeightF = 0.024        
     else
        xyres@tmXBLabelFontHeightF = 0.024
        xyres@tmYLLabelFontHeightF = 0.024
        xyres@gsnLeftStringFontHeightF = 0.030
        xyres@gsnCenterStringFontHeightF = 0.030
        xyres@gsnRightStringFontHeightF = 0.030    
     end if
     xyres@vpXF = 0.05
     xyres@vpHeightF = 0.3
     xyres@vpWidthF = 0.9
       
     xyres@gsnCenterString = ""
     
     xyres@trXMinF = syear-.5
     xyres@trXMaxF = eyear+1.5
     xyres@tiMainOn = False

     ; +++++++++++++++
     xyres@tmXBMode = "Manual"
     xyres@tmXBTickStartF = xyres@trXMinF
     xyres@tmXBTickEndF = xyres@trXMaxF
     xyres@tmYLMode = "Manual"
     xyres@trYMinF = -4.0
     xyres@trYMaxF = 4.0
     xyres@tmYLTickStartF = -4.0
     xyres@tmYLTickEndF = 4.0
     xyres@tmYLTickSpacingF = 2.0
    ; ++++++++++++++
    ;------------------------------- nino timeseries -----------------------------------------------------------------------------------
     xyres@gsnXYAboveFillColors = "red"
     xyres@gsnXYBelowFillColors = "blue"
     xyres@gsnLeftString = names
     arr = new((/2,dimsizes(nino34)/),typeof(nino34))
   
     tttt = dtrend_msg(ispan(0,dimsizes(nino34)-1,1),nino34,False,True)
     arr(0,:) = (/ nino34 /)
     arr(1,:) = (/  (ispan(0,dimsizes(nino34)-1,1)*tttt@slope)+tttt@y_intercept /)
     xyres@gsnRightString = "stddev: "+decimalPlaces(stddev(nino34), 2, True)+",  trend: "+decimalPlaces(tttt@slope*dimsizes(nino34),2,True)+nino34@units+" "+nyr+"yr~S~-1~N~"
     xyn34 = gsn_csm_xy(wks_n34,fspan(syear,eyear+.91667,dimsizes(nino34)),arr,xyres) 
     delete(tttt)
     
     tttt = dtrend_msg(ispan(0,dimsizes(nino3)-1,1),nino3,False,True)
     arr(0,:) = (/ nino3 /)
     arr(1,:) = (/  (ispan(0,dimsizes(nino3)-1,1)*tttt@slope)+tttt@y_intercept /)
     xyres@gsnRightString = "stddev: "+decimalPlaces(stddev(nino3), 2, True)+",  trend: "+decimalPlaces(tttt@slope*dimsizes(nino3),2,True)+nino3@units+" "+nyr+"yr~S~-1~N~"
     xyn3 = gsn_csm_xy(wks_n3,fspan(syear,eyear+.91667,dimsizes(nino3)),arr,xyres) 
     delete(tttt)
     
     tttt = dtrend_msg(ispan(0,dimsizes(nino4)-1,1),nino4,False,True)
     arr(0,:) = (/ nino4 /)
     arr(1,:) = (/  (ispan(0,dimsizes(nino4)-1,1)*tttt@slope)+tttt@y_intercept /)
     xyres@gsnRightString = "stddev: "+decimalPlaces(stddev(nino4), 2, True)+",  trend: "+decimalPlaces(tttt@slope*dimsizes(nino4),2,True)+nino4@units+" "+nyr+"yr~S~-1~N~"
     xyn4 = gsn_csm_xy(wks_n4,fspan(syear,eyear+.91667,dimsizes(nino4)),arr,xyres)
     delete(tttt)
     
     tttt = dtrend_msg(ispan(0,dimsizes(nino12)-1,1),nino12,False,True)
     arr(0,:) = (/ nino12 /)
     arr(1,:) = (/  (ispan(0,dimsizes(nino12)-1,1)*tttt@slope)+tttt@y_intercept /)
     xyres@gsnRightString = "stddev: "+decimalPlaces(stddev(nino12), 2, True)+",  trend: "+decimalPlaces(tttt@slope*dimsizes(nino12),2,True)+nino12@units+" "+nyr+"yr~S~-1~N~"
     xyn12 = gsn_csm_xy(wks_n12,fspan(syear,eyear+.91667,dimsizes(nino12)),arr,xyres)
     delete([/arr,tttt/])

    ;------------------------------ nino3 monthly standard deviation --------------------------------------------------------------------------
     xyres3 = xyres            ; resource list for monthly nino3 standard deviations
     xyres3@trXMinF = 0.5
     xyres3@trXMaxF = 12.5
     xyres3@vpWidthF = 0.65
     xyres3@vpHeightF = 0.35
     xyres3@trYMinF = 0.2
     xyres3@trYMaxF = 2.0
     ;xyres3@gsnAboveYRefLineColor = "gray50"
     xyres3@xyLineColor = "black"
     if (wks_type.eq."png") then
        xyres3@xyLineThicknessF = 3.5
     else
        xyres3@xyLineThicknessF = 1.75
     end if
     xyres3@tmXBMode    = "Explicit"        ; explicit labels
     xyres3@tmXBValues  = ispan(1,12,1)
     xyres3@tmXBLabels  = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
     xyres3@tmXTOn      = False
     xyres3@gsnLeftStringOrthogonalPosF = 0.025
     xyres3@gsnCenterStringOrthogonalPosF = xyres3@gsnLeftStringOrthogonalPosF
     xyres3@gsnRightStringOrthogonalPosF = xyres3@gsnLeftStringOrthogonalPosF
     xyres3@gsnRightStringFontHeightF = xyres3@gsnCenterStringFontHeightF
     xyres3@gsnLeftString = names
     xyres3@gsnCenterString = "NDJ/MAM ="+decimalPlaces((nino3_mon_sd(0)+nino3_mon_sd(10)+nino3_mon_sd(11))/3./dim_avg(nino3_mon_sd(2:4)), 2, True)
     xyres3@gsnRightString = "(~S~o~N~C)"
     xyn3_mst = gsn_csm_xy(wks_n3_mst,ispan(1,12,1),nino3_mon_sd,xyres3)
     delete([/xyres,xyres3/])
    ;- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     res = True
     res@vpHeightF = 0.45
     res@vpWidthF = 0.35
     res@gsnFrame = False
     res@gsnDraw = False 
  

     res@tmYLMode = "Explicit"
     ;     res@tmYLValues = ispan(0,72,6)
     ;     res@tmYLLabels = (/"Jan~S~-2~N~","Jul~S~-2~N~","Jan~S~-1~N~","Jul~S~-1~N~", \
     ;                     "Jan~S~0~N~","Jul~S~0~N~","Jan~S~+1~N~","Jul~S~+1~N~", \
     ;                     "Jan~S~+2~N~","Jul~S~+2~N~","Jan~S~+3~N~","Jul~S~+3~N~","Jan~S~+4~N~"/)
     res@trYMinF = 24
     res@trYMaxF = 52  
     res@tmYLValues = ispan(24,52,4)
     res@tmYLLabels = (/"Jan~S~0~N~","May~S~0~N~","Sep~S~0~N~","Jan~S~+1~N~", \
                     "May~S~+1~N~","Sep~S~+1~N~","Jan~S~+2~N~","May~S~+2~N~"/)    
     res@tmYLMinorValues = ispan(24,52,2)               
     res@tmYLLabelJust = "CenterCenter"               
     res@tmYLLabelDeltaF = 1.3    ;0.05
     res@cnFillOn = True
     res@gsnSpreadColors = True  
     res@gsnSpreadColorEnd = 19
     
     res@lbLabelBarOn = False
  
     res@tiMainOn = False
     res@cnInfoLabelOn = False
     res@cnLinesOn = True
     res@cnLevelSelectionMode = "ExplicitLevels"
     res@cnLevels = (/-3,-2.5,-2,-1.5,-1,-.75,-.5,-.25,0,.25,.5,.75,1,1.5,2,2.5,3/)   ;fspan(-2.,2.,17)
     carr = new(dimsizes(res@cnLevels),"string")
     carr = "transparent"
     carr(8) = "gray50"
     res@cnMonoLineColor = False
     res@cnLineColors = carr
     res@cnLineLabelsOn = False
     res@tmYLLabelFontHeightF = 0.014
     res@tmXBLabelFontHeightF = 0.014
     res@gsnMajorLonSpacing = 30.
     res@gsnMinorLonSpacing = 10.
     res@tiYAxisOn = False
     
     if (wks_type.eq."png") then
        res@cnLineThicknessF = 2.  
     else
        res@cnLineThicknessF = 1.  
     end if
  
     res@gsnCenterStringFontHeightF = 0.017
     res@gsnLeftStringFontHeightF = 0.017
     res@gsnRightStringFontHeightF = 0.017

     res@gsnLeftString = ""
     res@gsnCenterString= ""
     res@gsnRightString = ""


     res@gsnCenterString = names   ;"El Nin~H-13V2F35~D~FV-2H3F21~o"
     res@gsnRightString = cntr_hi
     plot_n3hi = gsn_csm_hov(wks_n3_tlon_hi,finsst_hi,res)
     
     res@gsnRightString = cntr_lo
     plot_n3lo = gsn_csm_hov(wks_n3_tlon_lo,finsst_lo,res)
     delete([/finsst_hi,finsst_lo,finsst_mid,cntr_hi,cntr_lo,cntr_mid/])
     delete(res)

    ;----------------------  nino3 power spectra----------------------------------------------------------------------------------------

     pres = True
     pres@vpXF = 0.07
     pres@trXMinF = 0.0
     pres@trXMaxF = 0.0832
     pres@tiYAxisString = "Power"              ; yaxis
     pres@xyLineColor   =  "black"         
     pres@gsnFrame      = False
     pres@gsnDraw       = False
  
     pres@tmXBLabelDeltaF = -.8
     pres@tmXTLabelDeltaF = -.8
     pres@pmLegendDisplayMode    = "Never"
     pres@xyLineThicknesses   = (/3.5,2.,1.,1./) 
     pres@xyDashPatterns      = (/0,0,0,0/)
     pres@xyLineColors        = (/"foreground","red","blue","green"/)      
     pres@xyLabelMode = "custom"
     pres@xyLineLabelFontColors = pres@xyLineColors 
     pres@xyExplicitLabels = (/"","",val1*100+"%",val2*100+"%"/)
     pres@tmXTOn = True
     pres@tmYROn = False
     pres@tmXTLabelsOn = True
     pres@tmXUseBottom = False
     pres@tmXTMode   = "Explicit"  
     pres@tmXBMode   = "Explicit"            
     pres@tmXTValues = (/".00167",".00833",".01667",".02778",".0416",".0556",".0832"/)
     pres@tmXTLabels = (/"50","10","5","3","2","1.5","1"/)           
     pres@tmXBValues = (/".0",".01",".02",".03",".042",".056",".083"/)
     pres@tmXBLabels =     pres@tmXBValues
     pres@tmXTLabelFontHeightF = 0.024
     pres@tmXBLabelFontHeightF = 0.024
     pres@tmYLLabelFontHeightF = 0.024
     pres@tiYAxisString = "Power (~S~o~N~C~S~2~N~ / cycles mo~S~-1~N~)"              ; yaxis
     pres@tiXAxisString = "Frequency (cycles mo~S~-1~N~)"
     pres@tiMainString = ""
     pres@txFontHeightF = 0.030
     pres@xyLineLabelFontHeightF = 0.017
     pres@tiXAxisFontHeightF = 0.032
     pres@tiYAxisFontHeightF = 0.032
     pres@tiMainFontHeightF = 0.032
     
     
     pres@trYMaxF = max(splt1(0,:))+15
       

     ; +++++++++++++++
     pres@tmYLMode = "Manual"
     pres@trYMinF = 0.
     pres@trYMaxF = 120.
     pres@tmYLTickStartF = 0.
     pres@tmYLTickEndF = 120.
     pres@tmYLTickSpacingF = 20.
     ; +++++++++++++++
     
     pres@tiMainOn = False
     pres@gsnCenterString = "Period (years)"
     pres@gsnCenterStringFontHeightF =    pres@tiYAxisFontHeightF
     pres@gsnRightString = "ratio ="+ dim_avg(where(sdof@frq.lt.0.083/3.and.sdof@frq.gt.0.083/8,splt1(0,:),splt1@_FillValue))/dim_avg(where(sdof@frq.lt.0.083/1.and.sdof@frq.gt.0.083/3,splt1(0,:),splt1@_FillValue))
     pres@gsnLeftString = names 
     
     pres@xyLineThicknessF   = 0.018
     pres@gsnCenterString = ""
     pspec = gsn_csm_xy(wks_n3_p,sdof@frq,splt1,pres) 
     delete([/sdof,splt1/])
    ;---------------------------------------------------------------------------------------------------------------
     title_n34 = nino34@area
     title_n4  = nino4@area
     title_n3  = nino3@area
     title_n12 = nino12@area

     ncol = floattointeger(sqrt(nsim))
     nrow = (nsim/ncol)+mod(nsim,ncol) 
     
     panres = True
     panres@gsnMaximize = True
     panres@gsnPaperOrientation = "portrait"
     panres@gsnPanelYWhiteSpacePercent = 3.0
     if (nsim.le.10) then
        panres@txFontHeightF = 0.024
     else
        panres@txFontHeightF = 0.016
     end if

     panres@gsnPanelScalePlotIndex = nyr
   
     
    lp = (/nrow,ncol/)
    
   
     panres@txString = "Nino3.4 (Monthly, "+title_n34+")"
     gsn_panel(wks_n34,xyn34,lp,panres)    
   
     panres@txString = "Nino4 (Monthly, "+title_n4+")"
     gsn_panel(wks_n4,xyn4,lp,panres)  
   
     panres@txString = "Nino3 (Monthly, "+title_n3+")"
     gsn_panel(wks_n3,xyn3,lp,panres)       
   
     panres@txString = "Nino1+2 (Monthly, "+title_n12+")"
     gsn_panel(wks_n12,xyn12,lp,panres)       
       
   
     panres@gsnPanelYWhiteSpacePercent = 0.5
     panres@txString = "nino3.4 standard deviation (Monthly)"
     gsn_panel(wks_n3_mst,xyn3_mst,(/nrow,ncol/),panres)   
     
     panres2 = True
     panres2@gsnMaximize = True
     panres2@gsnPaperOrientation = "portrait"
     panres2@gsnPanelLabelBar = True
     panres2@lbLabelStride = 1
     panres2@pmLabelBarWidthF = 0.4
     panres2@pmLabelBarHeightF = 0.06
     panres2@lbLabelFontHeightF = 0.013
     panres2@txString = ""   
     if (nsim.le.5) then
        panres2@txFontHeightF = 0.024
     else
        panres2@txFontHeightF = 0.016
     end if

     panres2@gsnPanelXWhiteSpacePercent = 8.5
     panres2@gsnPanelBottom = 0.05
   
     
     panres2@txString = "El Nin~H-13V2F35~D~FV-2H3F21~o Composite (3~S~o~N~S:3~S~o~N~N)"
     gsn_panel(wks_n3_tlon_hi,plot_n3hi,(/nrow,ncol/),panres2)   
     panres2@txString = "La Nin~H-13V2F35~D~FV-2H3F21~a Composite (3~S~o~N~S:3~S~o~N~N)"
     gsn_panel(wks_n3_tlon_lo,plot_n3lo,(/nrow,ncol/),panres2)   
     
     panres@txString = "nino3 (Monthly, detrended)"
     gsn_panel(wks_n3_p,pspec,(/nrow,ncol/),panres)  
     delete(wks_n3_p)
   
     delete([/xyn34,xyn4,xyn3,xyn12,plot_n3hi,plot_n3lo,pspec,pi,rad,wgt,lp/]) 

end
