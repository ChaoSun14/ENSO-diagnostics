load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
begin
    
    OUTPUT_TYPE = "png"
    OUTDIR      = "/home/reproduce/sunc/enso/plot/"
    CESMDIR     = "/home/reproduce/sunc/enso/CESM1-FASTCHEM_picontrol_atm/"
    FGOALSDIR   = "/home/reproduce/sunc/enso/FGOALS-g2_picontrol_atm/"
    OAFlux      = "/home/reproduce/sunc/enso/obs_data/OAFlux/"
   

    wks_type = OUTPUT_TYPE

    
    wks_swflux_fb = gsn_open_wks(wks_type,OUTDIR+"scatter_net_sw_flux_feedback")
    wks_lhflux_fb = gsn_open_wks(wks_type,OUTDIR+"scatter_net_lh_flux_feedback")
    ;---------------- obs data --------------------------------------------------------------------
    print("================================================================")
    print("obs is reading")
    print("================================================================")
    data = addfile(OAFlux+"lw_isccp_1984-2009.nc","r") 
    lwrs_obs = short2flt(data->nlwrs(0:26*12-1,:,:))     ; 312 ;W/m^2 positive upward
    data = addfile(OAFlux+"sw_isccp_1984-2009.nc","r") 
    swrs_obs = short2flt(data->nswrs(0:26*12-1,:,:))           ;W/m^2 positive downward
    data = addfile(OAFlux+"sh_oaflux_1958-2014.nc","r") 
    shfx_obs = short2flt(data->shtfl(0:57*12-1,:,:))     ; 684 
    data = addfile(OAFlux+"lh_oaflux_1958-2014.nc","r") 
    lhfx_obs = short2flt(data->lhtfl(0:57*12-1,:,:))   
    data = addfile(OAFlux+"ts_oaflux_1958-2014.nc","r")
    sst_obs  = short2flt(data->tmpsf(0:57*12-1,:,:)) 
    ;sst_obs = where(sst_obs.le.-1.8,-1.8,sst_obs)
    lat_obs = data->lat
    lon_obs = data->lon
    lwrs_obs@_FillValue =   32766
    swrs_obs@_FillValue =   32766
    shfx_obs@_FillValue =   32766
    lhfx_obs@_FillValue =   32766
    sst_obs@_FillValue  =   32766
    
    
    d = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")   
    basemap = d->LSMASK
    lsm = landsea_mask(basemap,lat_obs,lon_obs)
    lwrs_obs = mask(lwrs_obs,conform(lwrs_obs,lsm,(/1,2/)).ge.1,False)
    swrs_obs = mask(swrs_obs,conform(swrs_obs,lsm,(/1,2/)).ge.1,False)
    shfx_obs = mask(shfx_obs,conform(shfx_obs,lsm,(/1,2/)).ge.1,False)
    lhfx_obs = mask(lhfx_obs,conform(lhfx_obs,lsm,(/1,2/)).ge.1,False)
    sst_obs  = mask(sst_obs,conform(sst_obs,lsm,(/1,2/)).ge.1,False)
    delete(lsm)
    
    ;lwrs_obs=lonFlip(lwrs_obs)
    ;swrs_obs=lonFlip(swrs_obs)
    ;shfx_obs=lonFlip(shfx_obs)
    ;lhfx_obs=lonFlip(lhfx_obs)
    ;sst_obs = lonFlip(sst_obs) 
    ;printVarSummary(lwrs_obs)

    pi=4.*atan(1.0)
    rad=(pi/180.)
    coswgt_obs=cos(rad*lat_obs)
    coswgt_obs!0 = "lat" 
    coswgt_obs&lat= lat_obs 
    ;printVarSummary(sst_obs)
    ;printVarSummary(coswgt_obs)

    qnet_obs = swrs_obs(:,:,:)-lwrs_obs(:,:,:)-shfx_obs(26*12:52*12-1,:,:)-lhfx_obs(26*12:52*12-1,:,:)     ;(198401-200912) ;26 years
    copy_VarCoords(swrs_obs, qnet_obs)
    lwrs_obs_anl = rmMonAnnCycTLL(lwrs_obs)
    swrs_obs_anl = rmMonAnnCycTLL(swrs_obs)
    shfx_obs_anl = rmMonAnnCycTLL(shfx_obs)
    lhfx_obs_anl = rmMonAnnCycTLL(lhfx_obs)
    qnet_obs_anl = rmMonAnnCycTLL(qnet_obs)
    sst_obs_anl  = rmMonAnnCycTLL(sst_obs)
    copy_VarCoords(lwrs_obs, lwrs_obs_anl)
    copy_VarCoords(swrs_obs, swrs_obs_anl)
    copy_VarCoords(shfx_obs, shfx_obs_anl)
    copy_VarCoords(lhfx_obs, lhfx_obs_anl)
    copy_VarCoords(qnet_obs, qnet_obs_anl)
    copy_VarCoords(sst_obs, sst_obs_anl)
    printVarSummary(sst_obs_anl)

    nino3_obs = wgt_areaave_Wrap(sst_obs_anl(:,{-5:5},{210:270}),coswgt_obs({-5:5}),1.0,1)
    nino3_obs_swrs = wgt_areaave_Wrap(swrs_obs_anl(:,{-5:5},{210:270}),coswgt_obs({-5:5}),1.0,1)
    nino3_obs_lhfx = wgt_areaave_Wrap(lhfx_obs_anl(:,{-5:5},{210:270}),coswgt_obs({-5:5}),1.0,1)

    printVarSummary(nino3_obs)
    

    robs_qnet = regCoef(nino3_obs(26*12:52*12-1),qnet_obs_anl(lat|:,lon|:,time|:))
    robs_lwrs = regCoef(nino3_obs(26*12:52*12-1),lwrs_obs_anl(lat|:,lon|:,time|:))
    robs_swrs = regCoef(nino3_obs(26*12:52*12-1),swrs_obs_anl(lat|:,lon|:,time|:))
    robs_shfx = regCoef(nino3_obs(:),shfx_obs_anl(lat|:,lon|:,time|:))
    robs_lhfx = regCoef(nino3_obs(:),lhfx_obs_anl(lat|:,lon|:,time|:))

    robs_qnet!0 = "lat"
    robs_qnet&lat = sst_obs&lat
    robs_qnet!1 = "lon"
    robs_qnet&lon = sst_obs&lon
    copy_VarCoords(robs_qnet, robs_lwrs)
    copy_VarCoords(robs_qnet, robs_swrs)
    copy_VarCoords(robs_qnet, robs_shfx)
    copy_VarCoords(robs_qnet, robs_lhfx)

     
     
    ;---------------------- CESM -------------------------------------------------------------------------------
    print("================================================================")
    print("CESM is reading")
    print("================================================================")
    data = addfile(CESMDIR+"rsds_Amon_CESM1-FASTCHEM_piControl_r1i1p1_007001-029112.nc","r") 
    rsds_cesm = data->rsds(0:12*200-1,:,:)    ;W/m^2 positive downward
    data = addfile(CESMDIR+"rsus_Amon_CESM1-FASTCHEM_piControl_r1i1p1_007001-029112.nc","r") 
    rsus_cesm = data->rsus(0:12*200-1,:,:) 
    data = addfile(CESMDIR+"rlds_Amon_CESM1-FASTCHEM_piControl_r1i1p1_007001-029112.nc","r") 
    rlds_cesm = data->rlds(0:12*200-1,:,:) 
    data = addfile(CESMDIR+"rlus_Amon_CESM1-FASTCHEM_piControl_r1i1p1_007001-029112.nc","r") 
    rlus_cesm = data->rlus(0:12*200-1,:,:)
    data = addfile(CESMDIR+"hfls_Amon_CESM1-FASTCHEM_piControl_r1i1p1_007001-029112.nc","r") 
    hfls_cesm = data->hfls(0:12*200-1,:,:) 
    data = addfile(CESMDIR+"hfss_Amon_CESM1-FASTCHEM_piControl_r1i1p1_007001-029112.nc","r") 
    hfss_cesm = data->hfss(0:12*200-1,:,:)                 
    lat_cesm = data->lat
    lon_cesm = data->lon
    time_cesm = data->time
    
    rsns_cesm = rsds_cesm-rsus_cesm
    rlns_cesm = rlus_cesm-rlds_cesm
    qnet_cesm = rsds_cesm-rsus_cesm+rlds_cesm-rlus_cesm-hfls_cesm-hfss_cesm
    copy_VarCoords(rsds_cesm, qnet_cesm)
    copy_VarCoords(rsds_cesm, rsns_cesm)
    copy_VarCoords(rsds_cesm, rlns_cesm)
    ;printVarSummary(qnet_cesm)

    pi=4.*atan(1.0)
    rad=(pi/180.)
    coswgt_cesm=cos(rad*lat_cesm)
    coswgt_cesm!0 = "lat" 
    coswgt_cesm&lat= lat_cesm 

    data = addfile(CESMDIR+"ts_data/"+"CESM1-FASTCHEM_piControl_r1i1p1.ts.007001-029112.nc","r")
    sst_cesm  = data->ts(0:200*12-1,:,:)
    sst_cesm =(/sst_cesm-273.15/)
    sst_cesm = where(sst_cesm.le.-1.8,-1.8,sst_cesm)
    d = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")   
    basemap = d->LSMASK
    lsm = landsea_mask(basemap,lat_cesm,lon_cesm)
    qnet_cesm = mask(qnet_cesm,conform(qnet_cesm,lsm,(/1,2/)).ge.1,False)
    rsns_cesm = mask(rsns_cesm,conform(rsns_cesm,lsm,(/1,2/)).ge.1,False)
    rlns_cesm = mask(rlns_cesm,conform(rlns_cesm,lsm,(/1,2/)).ge.1,False)
    hfls_cesm = mask(hfls_cesm,conform(hfls_cesm,lsm,(/1,2/)).ge.1,False)
    hfss_cesm = mask(hfss_cesm,conform(hfss_cesm,lsm,(/1,2/)).ge.1,False)
    sst_cesm  = mask(sst_cesm,conform(sst_cesm,lsm,(/1,2/)).ge.1,False)
    delete(lsm)

    rsns_cesm_anl = rmMonAnnCycTLL(rsns_cesm)
    rlns_cesm_anl = rmMonAnnCycTLL(rlns_cesm)
    hfls_cesm_anl = rmMonAnnCycTLL(hfls_cesm)
    hfss_cesm_anl = rmMonAnnCycTLL(hfss_cesm)
    qnet_cesm_anl = rmMonAnnCycTLL(qnet_cesm)
    sst_cesm_anl = rmMonAnnCycTLL(sst_cesm)
    copy_VarCoords(sst_cesm, sst_cesm_anl)


    nino3_cesm = wgt_areaave_Wrap(sst_cesm_anl(:,{-5:5},{210:270}),coswgt_cesm({-5:5}),1.0,1)
    nino3_cesm_swrs = wgt_areaave_Wrap(rsns_cesm_anl(:,{-5:5},{210:270}),coswgt_cesm({-5:5}),1.0,1)
    nino3_cesm_lhfx = wgt_areaave_Wrap(hfls_cesm_anl(:,{-5:5},{210:270}),coswgt_cesm({-5:5}),1.0,1)
    printVarSummary(nino3_cesm)


    rcesm_qnet = regCoef(nino3_cesm,qnet_cesm_anl(lat|:,lon|:,time|:))
    rcesm_lwrs = regCoef(nino3_cesm,rlns_cesm_anl(lat|:,lon|:,time|:))
    rcesm_swrs = regCoef(nino3_cesm,rsns_cesm_anl(lat|:,lon|:,time|:))
    rcesm_shfx = regCoef(nino3_cesm,hfss_cesm_anl(lat|:,lon|:,time|:))
    rcesm_lhfx = regCoef(nino3_cesm,hfls_cesm_anl(lat|:,lon|:,time|:))

    rcesm_qnet!0 = "lat"
    rcesm_qnet&lat = sst_cesm&lat
    rcesm_qnet!1 = "lon"
    rcesm_qnet&lon = sst_cesm&lon
    copy_VarCoords(rcesm_qnet, rcesm_lwrs)
    copy_VarCoords(rcesm_qnet, rcesm_swrs)
    copy_VarCoords(rcesm_qnet, rcesm_shfx)
    copy_VarCoords(rcesm_qnet, rcesm_lhfx)


    ;--------------------------------- fgoals-g2 -------------------------------------------------------
    print("================================================================")
    print("FGOALS-g2 is reading")
    print("================================================================")
    data = addfile(FGOALSDIR+"t0225_ocn_S.gamil.h0.SurfaceFLX.060001-079912.nc","r") 
    flns_fgoals = data->FLNS(:,:,:)    ;W/m^2 positive downward
    fsns_fgoals = data->FSNS(:,:,:)
    lhfx_fgoals = data->LHFLX(:,:,:)
    shfx_fgoals = data->SHFLX(:,:,:)
    lat_fgoals = data->lat
    lon_fgoals = data->lon
    time_fgoals = data->time

    qnet_fgoals = fsns_fgoals-flns_fgoals-lhfx_fgoals-shfx_fgoals
    copy_VarCoords(flns_fgoals, qnet_fgoals)
    ;printVarSummary(qnet_fgoals)
   
    pi=4.*atan(1.0)
    rad=(pi/180.)
    coswgt_fgoals=cos(rad*lat_fgoals)
    coswgt_fgoals!0 = "lat" 
    coswgt_fgoals&lat= lat_fgoals 

    data = addfile(FGOALSDIR+"ts_data/"+"t0225_ocn_S.gamil.h0.ts.060001-079912.nc","r")
    sst_fgoals  = data->TS(:,:,:)
    ;printVarSummary(sst_fgoals)
    sst_fgoals =(/sst_fgoals-273.15/)
    sst_fgoals = where(sst_fgoals.le.-1.8,-1.8,sst_fgoals)

    d = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")   
    basemap = d->LSMASK
    lsm = landsea_mask(basemap,lat_fgoals,lon_fgoals)
    qnet_fgoals = mask(qnet_fgoals,conform(qnet_fgoals,lsm,(/1,2/)).ge.1,False)
    fsns_fgoals = mask(fsns_fgoals,conform(fsns_fgoals,lsm,(/1,2/)).ge.1,False)
    flns_fgoals = mask(flns_fgoals,conform(flns_fgoals,lsm,(/1,2/)).ge.1,False)
    lhfx_fgoals = mask(lhfx_fgoals,conform(lhfx_fgoals,lsm,(/1,2/)).ge.1,False)
    shfx_fgoals = mask(shfx_fgoals,conform(shfx_fgoals,lsm,(/1,2/)).ge.1,False)
    sst_fgoals  = mask(sst_fgoals,conform(sst_fgoals,lsm,(/1,2/)).ge.1,False)
    delete(lsm)

    rsns_fgoals_anl = rmMonAnnCycTLL(fsns_fgoals)
    rlns_fgoals_anl = rmMonAnnCycTLL(flns_fgoals)
    hfls_fgoals_anl = rmMonAnnCycTLL(lhfx_fgoals)
    hfss_fgoals_anl = rmMonAnnCycTLL(shfx_fgoals)
    qnet_fgoals_anl = rmMonAnnCycTLL(qnet_fgoals)
    sst_fgoals_anl  = rmMonAnnCycTLL(sst_fgoals)
    nino3_fgoals = wgt_areaave_Wrap(sst_fgoals_anl(:,{-5:5},{210:270}),coswgt_fgoals({-5:5}),1.0,1)
    nino3_fgoals_swrs = wgt_areaave_Wrap(rsns_fgoals_anl(:,{-5:5},{210:270}),coswgt_fgoals({-5:5}),1.0,1)
    nino3_fgoals_lhfx = wgt_areaave_Wrap(hfls_fgoals_anl(:,{-5:5},{210:270}),coswgt_fgoals({-5:5}),1.0,1)
    printVarSummary(nino3_fgoals)


    rfgoals_qnet = regCoef(nino3_fgoals,qnet_fgoals_anl(lat|:,lon|:,time|:))
    rfgoals_lwrs = regCoef(nino3_fgoals,rlns_fgoals_anl(lat|:,lon|:,time|:))
    rfgoals_swrs = regCoef(nino3_fgoals,rsns_fgoals_anl(lat|:,lon|:,time|:))
    rfgoals_shfx = regCoef(nino3_fgoals,hfss_fgoals_anl(lat|:,lon|:,time|:))
    rfgoals_lhfx = regCoef(nino3_fgoals,hfls_fgoals_anl(lat|:,lon|:,time|:))

    rfgoals_qnet!0 = "lat"
    rfgoals_qnet&lat = sst_fgoals&lat
    rfgoals_qnet!1 = "lon"
    rfgoals_qnet&lon = sst_fgoals&lon
    copy_VarCoords(rfgoals_qnet, rfgoals_lwrs)
    copy_VarCoords(rfgoals_qnet, rfgoals_swrs)
    copy_VarCoords(rfgoals_qnet, rfgoals_shfx)
    copy_VarCoords(rfgoals_qnet, rfgoals_lhfx)

;------------------------------------------------------------------------------------------------------------------------------
   rc_obs = regline(nino3_obs(0:26*12-1),nino3_obs_swrs)
   rcline_obs = rc_obs*(nino3_obs(0:26*12-1)-rc_obs@xave) + rc_obs@yave
   rc_cesm = regline(nino3_cesm(0:26*12-1),nino3_cesm_swrs(0:26*12-1))
   rcline_cesm = rc_cesm*(nino3_cesm(0:26*12-1)-rc_cesm@xave) + rc_cesm@yave
   rc_fgoals = regline(nino3_fgoals(0:26*12-1),nino3_fgoals_swrs(0:26*12-1))
   rcline_fgoals = rc_fgoals*(nino3_fgoals(0:26*12-1)-rc_fgoals@xave) + rc_fgoals@yave


    x1 = (/nino3_obs(0:26*12-1),nino3_cesm(0:26*12-1),nino3_fgoals(0:26*12-1),nino3_obs(0:26*12-1),nino3_cesm(0:26*12-1),nino3_fgoals(0:26*12-1)/)
    y1 = (/nino3_obs_swrs,nino3_cesm_swrs(0:26*12-1),nino3_fgoals_swrs(0:26*12-1),rcline_obs,rcline_cesm,rcline_fgoals/)

   rc_lhfx_obs = regline(nino3_obs(0:57*12-1),-nino3_obs_lhfx(0:57*12-1))
   rc_lhfxline_obs = rc_lhfx_obs*(nino3_obs(0:57*12-1)-rc_lhfx_obs@xave) + rc_lhfx_obs@yave
   rc_lhfx_cesm = regline(nino3_cesm(0:57*12-1),-nino3_cesm_lhfx(0:57*12-1))
   rc_lhfxline_cesm = rc_lhfx_cesm*(nino3_cesm(0:57*12-1)-rc_lhfx_cesm@xave) + rc_lhfx_cesm@yave
   rc_lhfx_fgoals = regline(nino3_fgoals(0:57*12-1),-nino3_fgoals_lhfx(0:57*12-1))
   rc_lhfxline_fgoals = rc_lhfx_fgoals*(nino3_fgoals(0:57*12-1)-rc_lhfx_fgoals@xave) + rc_lhfx_fgoals@yave

    x2 = (/nino3_obs,nino3_cesm(0:57*12-1),nino3_fgoals(0:57*12-1),nino3_obs,nino3_cesm(0:57*12-1),nino3_fgoals(0:57*12-1)/)
    y2 = (/-nino3_obs_lhfx,-nino3_cesm_lhfx(0:57*12-1),-nino3_fgoals_lhfx(0:57*12-1),rc_lhfxline_obs,rc_lhfxline_cesm,rc_lhfxline_fgoals/)
    xyres               = True 
    xyres@gsnDraw       = False
    xyres@gsnFrame      = False
    xyres@vpHeightF     = 0.3
    xyres@vpWidthF      = 0.7
    xyres@gsnLeftString =""
    xyres@gsnRightString =""
    xyres@gsnStringFontHeightF = 0.022
    xyres@gsnLeftStringFontHeightF = 0.022
    xyres@gsnRightStringFontHeightF = 0.022
    xyres@tmXBLabelFontHeightF  = 0.020
    ;xyres@tmLabelAutoStride = True
    xyres@tmYLLabelFontHeightF  = 0.020
    ;xyres@tmYLMode = "Manual"
    ;xyres@tmYLTickStartF = -0.04
    ;xyres@tmYLTickEndF = 0.08
    ;xyres@tmYLTickSpacingF = 0.02
    xyres@xyMarkLineModes     = (/"Markers","Markers","Markers","Lines","Lines","Lines"/)
    xyres@xyMarkers           = (/16,10,4/)                     ; choose type of marker 
    xyres@xyMarkerColors       = (/"black","red","blue","black","red","blue"/)                  ; Marker color
    xyres@xyMarkerSizeF       = 0.005  
    ;xyres@tmXBMode     = "Explicit"
    ;xyres@tmXBValues   =(/120.,140.,160.,180.,200.,220.,240.,260.,280./)
    ;xyres@tmXBLabels   = (/"120~S~o~N~E","140~S~o~N~E","160~S~o~N~E","180","160~S~o~N~W","140~S~o~N~W","120~S~o~N~W","100~S~o~N~W","90~S~o~N~W"/)
    ;xyres@xyLineThicknesses   = (/5.,5.,5./)      ; line thicknesses
    ;xyres@xyLineColors      = (/"black","blue","red"/)
    
    xyres@xyDashPatterns  = (/0.,0.,0.,0.,0.,0./)      ; make all lines solid    
  
  ;---Turn off some tickmarks
    xyres@tmXTOn                 = False    
    xyres@tmYROn                 = False 
    xyres@tmXBOn                 = True
    xyres@tmXBLabelsOn           = True    

      
    xyres@gsnLeftString    = ""
    xyres@gsnCenterString  = "Nino3 shortwave feedback"
    xyres@tiYAxisString    = "Shortwave anomaly"+"   (W/m~S~2~N~)"
    xyres@tiXAxisString    = "SST anomaly" +"   (~S~o~N~C)"
    plot3=gsn_csm_xy(wks_swflux_fb,x1,y1,xyres)

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.015                 ; change text size
    res_text@txJust             = "BottomLeft"          ; text justification
  
    res_marker                   = True                  
    res_marker@gsMarkerIndex     = 10                     ; choose type of marker 
    res_marker@gsMarkerColor    = "red"                 ; Marker color
    res_marker@gsMarkerSizeF     = 0.005 
    xx = (/-2.8,-2.7,-2.6/)
    yy = (/-55,-55,-55/)
    dum1 = gsn_add_polymarker(wks_swflux_fb,plot3,xx,yy,res_marker)              ; add polyline
    dum2 = gsn_add_text(wks_swflux_fb,plot3,"CESM"+"  ~F33~a="+sprintf("%5.3f",wgt_areaave_Wrap(rcesm_swrs({-5:5},{210:270}),coswgt_cesm({-5:5}),1.0,1)),-2.5,-56.5,res_text); add text
    
    yy = (/-45,-45,-45/)
    res_marker@gsMarkerIndex          = 4   
    res_marker@gsMarkerColor    = "blue"
    dum3 = gsn_add_polymarker(wks_swflux_fb,plot3,xx,yy,res_marker)                ; add polymarker
    dum4 = gsn_add_text(wks_swflux_fb,plot3,"FGOALS-g2"+"  ~F33~a="+sprintf("%5.3f",wgt_areaave_Wrap(rfgoals_swrs({-5:5},{210:270}),coswgt_fgoals({-5:5}),1.0,1)),-2.5,-47,res_text)       ; add text
    
    yy = (/-35,-35,-35/)
    res_marker@gsMarkerIndex         = 16   
    res_marker@gsMarkerColor    = "black"
    dum5 = gsn_add_polymarker(wks_swflux_fb,plot3,xx,yy,res_marker)                ; add polymarker
    dum6 = gsn_add_text(wks_swflux_fb,plot3,"OAFlux"+"  ~F33~a="+sprintf("%5.3f",wgt_areaave_Wrap(robs_swrs({-5:5},{210:270}),coswgt_obs({-5:5}),1.0,1)),-2.5,-37,res_text)     ; add text
    res_lines = True
    res_lines@gsLineDashPattern = 1. 
    ;res_lines@gsLineThicknessF  = 5.     
    xx1 = (/-3.,4./)
    yy1 = (/0.0,0.0/)
    dum7 = gsn_add_polyline(wks_swflux_fb,plot3,xx1,yy1,res_lines)
    xx1 = (/0.,0./)
    yy1 = (/-60,40/)
    dum8 = gsn_add_polyline(wks_swflux_fb,plot3,xx1,yy1,res_lines)
    draw(wks_swflux_fb)
    frame(wks_swflux_fb)

    xyres@gsnCenterString  = "Nino3 latent heat feedback"
    xyres@tiYAxisString    = "Latent heat anomaly"+"   (W/m~S~2~N~)"
    xyres@tiXAxisString    = "SST anomaly" +"   (~S~o~N~C)"
    plot4=gsn_csm_xy(wks_lhflux_fb,x2,y2,xyres)

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.015                 ; change text size
    res_text@txJust             = "BottomLeft"          ; text justification
  
    res_marker                   = True                  
    res_marker@gsMarkerIndex     = 10                     ; choose type of marker 
    res_marker@gsMarkerColor    = "red"                 ; Marker color
    res_marker@gsMarkerSizeF     = 0.005 
    xx = (/-2.8,-2.7,-2.6/)
    yy = (/-55,-55,-55/)
    dum11 = gsn_add_polymarker(wks_lhflux_fb,plot4,xx,yy,res_marker)              ; add polyline
    dum12 = gsn_add_text(wks_lhflux_fb,plot4,"CESM"+"  ~F33~a="+sprintf("%5.3f",-wgt_areaave_Wrap(rcesm_lhfx({-5:5},{210:270}),coswgt_cesm({-5:5}),1.0,1)),-2.5,-57,res_text); add text
    
    yy = (/-45,-45,-45/)
    res_marker@gsMarkerIndex          = 4   
    res_marker@gsMarkerColor    = "blue"
    dum13 = gsn_add_polymarker(wks_lhflux_fb,plot4,xx,yy,res_marker)                ; add polymarker
    dum14 = gsn_add_text(wks_lhflux_fb,plot4,"FGOALS-g2"+"  ~F33~a="+sprintf("%5.3f",-wgt_areaave_Wrap(rfgoals_lhfx({-5:5},{210:270}),coswgt_fgoals({-5:5}),1.0,1)),-2.5,-47,res_text)       ; add text
    
    yy = (/-35,-35,-35/)
    res_marker@gsMarkerIndex         = 16   
    res_marker@gsMarkerColor    = "black"
    dum15 = gsn_add_polymarker(wks_lhflux_fb,plot4,xx,yy,res_marker)                ; add polymarker
    dum16 = gsn_add_text(wks_lhflux_fb,plot4,"OAFlux"+"  ~F33~a="+sprintf("%5.3f",-wgt_areaave_Wrap(robs_lhfx({-5:5},{210:270}),coswgt_obs({-5:5}),1.0,1)),-2.5,-37,res_text)     ; add text
    res_lines = True
    res_lines@gsLineDashPattern = 1. 
    ;res_lines@gsLineThicknessF  = 5.     
    xx1 = (/-3.,4./)
    yy1 = (/0.0,0.0/)
    dum17 = gsn_add_polyline(wks_lhflux_fb,plot4,xx1,yy1,res_lines)
    xx1 = (/0.,0./)
    yy1 = (/-60,40/)
    dum18 = gsn_add_polyline(wks_lhflux_fb,plot4,xx1,yy1,res_lines)
    draw(wks_lhflux_fb)
    frame(wks_lhflux_fb)
end