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

    wks_nhflux_fb = gsn_open_wks(wks_type,OUTDIR+"net_heat_flux_feedback")
    wks_lwflux_fb = gsn_open_wks(wks_type,OUTDIR+"net_lw_flux_feedback")
    wks_swflux_fb = gsn_open_wks(wks_type,OUTDIR+"net_sw_flux_feedback")
    wks_lhflux_fb = gsn_open_wks(wks_type,OUTDIR+"net_lh_flux_feedback")
    wks_shflux_fb = gsn_open_wks(wks_type,OUTDIR+"net_sh_flux_feedback")
    wks_swflux_fb_mon = gsn_open_wks(wks_type,OUTDIR+"sw_flux_feedback_annual_cycle")
    wks_lhflux_fb_mon = gsn_open_wks(wks_type,OUTDIR+"lh_flux_feedback_annual_cycle")
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

    robs_swrs_mon = new((/12,dimsizes(lat_obs),dimsizes(lon_obs)/), typeof(robs_swrs))
    robs_lhfx_mon = new((/12,dimsizes(lat_obs),dimsizes(lon_obs)/),typeof(robs_lhfx))
     
     do hh = 0,11
        robs_swrs_mon(hh,:,:) = regCoef(nino3_obs(26*12+hh:52*12-1:12),swrs_obs_anl(lat|:,lon|:,time|hh::12)) 
        robs_lhfx_mon(hh,:,:) = regCoef(nino3_obs(hh::12),lhfx_obs_anl(lat|:,lon|:,time|hh::12)) 
     end do
     copy_VarCoords(robs_qnet, robs_swrs_mon(0,:,:))
     copy_VarCoords(robs_qnet, robs_lhfx_mon(0,:,:))
     
     robs_swrs_mon_nino3 = wgt_areaave_Wrap(robs_swrs_mon(:,{-5:5},{210:270}),coswgt_obs({-5:5}),1.0,1)
     robs_lhfx_mon_nino3 = wgt_areaave_Wrap(robs_lhfx_mon(:,{-5:5},{210:270}),coswgt_obs({-5:5}),1.0,1)

    
    robs_qnet=lonFlip(robs_qnet)
    robs_lwrs=lonFlip(robs_lwrs)
    robs_swrs=lonFlip(robs_swrs)
    robs_shfx=lonFlip(robs_shfx)
    robs_lhfx=lonFlip(robs_lhfx)
    robs_lhfx_mon=lonFlip(robs_lhfx_mon)
    robs_swrs_mon=lonFlip(robs_swrs_mon)
    printVarSummary(robs_qnet)
    printVarSummary(robs_lhfx_mon)
    print(robs_lhfx_mon_nino3)
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
    rcesm_swrs_mon = new((/12,dimsizes(lat_cesm),dimsizes(lon_cesm)/), typeof(rcesm_swrs))
    rcesm_lhfx_mon = new((/12,dimsizes(lat_cesm),dimsizes(lon_cesm)/),typeof(rcesm_lhfx))
     
     do hh = 0,11
        rcesm_swrs_mon(hh,:,:) = regCoef(nino3_cesm(hh::12),rsns_cesm_anl(lat|:,lon|:,time|hh::12)) 
        rcesm_lhfx_mon(hh,:,:) = regCoef(nino3_cesm(hh::12),hfls_cesm_anl(lat|:,lon|:,time|hh::12)) 
     end do
     copy_VarCoords(rcesm_qnet, rcesm_swrs_mon(0,:,:))
     copy_VarCoords(rcesm_qnet, rcesm_lhfx_mon(0,:,:))
     
     rcesm_swrs_mon_nino3 = wgt_areaave_Wrap(rcesm_swrs_mon(:,{-5:5},{210:270}),coswgt_cesm({-5:5}),1.0,1)
     rcesm_lhfx_mon_nino3 = wgt_areaave_Wrap(rcesm_lhfx_mon(:,{-5:5},{210:270}),coswgt_cesm({-5:5}),1.0,1)
    
    printVarSummary(rcesm_qnet)

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
    rfgoals_swrs_mon = new((/12,dimsizes(lat_fgoals),dimsizes(lon_fgoals)/), typeof(rfgoals_swrs))
    rfgoals_lhfx_mon = new((/12,dimsizes(lat_fgoals),dimsizes(lon_fgoals)/),typeof(rfgoals_lhfx))
     
     do hh = 0,11
        rfgoals_swrs_mon(hh,:,:) = regCoef(nino3_fgoals(hh::12),rsns_fgoals_anl(lat|:,lon|:,time|hh::12)) 
        rfgoals_lhfx_mon(hh,:,:) = regCoef(nino3_fgoals(hh::12),hfls_fgoals_anl(lat|:,lon|:,time|hh::12)) 
     end do
     copy_VarCoords(rfgoals_qnet, rfgoals_swrs_mon(0,:,:))
     copy_VarCoords(rfgoals_qnet, rfgoals_lhfx_mon(0,:,:))
     
     rfgoals_swrs_mon_nino3 = wgt_areaave_Wrap(rfgoals_swrs_mon(:,{-5:5},{210:270}),coswgt_fgoals({-5:5}),1.0,1)
     rfgoals_lhfx_mon_nino3 = wgt_areaave_Wrap(rfgoals_lhfx_mon(:,{-5:5},{210:270}),coswgt_fgoals({-5:5}),1.0,1)

    printVarSummary(rfgoals_qnet)

;----------------------------------------------------------------------------------
    

;---------------------------  wks_nhflux_fb  --------------------------------------
    plot1 = new(3,graphic)
    res1=True
    res1@gsnDraw        = False            ; don't draw yet
    res1@gsnFrame       = False
    ;res1@mpCenterLonF   = 180.
    res1@vpHeightF      = 0.3
    res1@vpWidthF       = 0.7
    res1@gsnStringFontHeightF      = 0.022
    res1@gsnLeftStringFontHeightF  = 0.022
    res1@gsnRightStringFontHeightF = 0.022
    res1@tmXBLabelFontHeightF  = 0.022
    res1@tmLabelAutoStride     = True
    res1@tmYLLabelFontHeightF  = 0.022
    res1@cnFillOn  = True
    res1@cnLinesOn = False
    res1@lbLabelBarOn = True
    res1@mpMinLatF  = -30.
    res1@mpMaxLatF  = 30.
    res1@mpMinLonF  = 120.
    res1@mpMaxLonF  = 290.
    res1@mpCenterLonF   = 205.
    ;res1@gsnSpreadColors = True
    ;res1@gsnSpreadColorStart = -60
    res1@lbOrientation       = "vertical"
    res1@cnLevelSelectionMode = "ManualLevels"
    res1@cnLevelSpacingF = 4.
    res1@cnMinLevelValF = -24.
    res1@cnMaxLevelValF = 24.
    

    dum = new(24,graphic)
    dum2= new(24,graphic)
    ypts4 = (/5,5,-5,-5,5/)
    xpts4 = (/160,210,210,160,160/)
    ypts3 = (/5,5,-5,-5,5/)
    xpts3 = (/210,270,270,210,210/)

    resl4 = True
    resl4@gsLineColor = "red"
    resl4@gsLineThicknessF = 3.
;   resl4@gsLineLabelString = "Nino 4"
    resl3 = True
    resl3@gsLineColor = "black"
    resl3@gsLineThicknessF = 3.
;   resl3@gsLineLabelString = "Nino 3"
     
    
    res1@gsnLeftString  ="OAFlux"
    res1@gsnRightString ="~F33~a"+"~F25~~B~net~N~="+sprintf("%5.2f", wgt_areaave_Wrap(robs_qnet({-5:5},{-150:-90}),coswgt_obs({-5:5}),1.0,1))+"(Wm~S~-2~N~/~S~o~N~C)"
    plot1(0) = gsn_csm_contour_map_ce(wks_nhflux_fb,robs_qnet,res1)

    do i=0,3
        dum(i) = gsn_add_polyline(wks_nhflux_fb,plot1(0),xpts4(i:i+1),ypts4(i:i+1),resl4)
        dum2(i) = gsn_add_polyline(wks_nhflux_fb,plot1(0),xpts3(i:i+1),ypts3(i:i+1),resl3)
    end do

    res1@gsnLeftString  ="CESM"
    ;res1@gsnCenterString = "RMSE="+sprintf("%3.2f", abs(qnet_cesm_mon_ave/qnet_obs_mon_ave-1.))
    res1@gsnRightString ="~F33~a"+"~F25~~B~net~N~="+sprintf("%5.2f", wgt_areaave_Wrap(rcesm_qnet({-5:5},{210:270}),coswgt_cesm({-5:5}),1.0,1))+"(Wm~S~-2~N~/~S~o~N~C)"
    plot1(1) = gsn_csm_contour_map_ce(wks_nhflux_fb,rcesm_qnet,res1)
    do i=0,3
        dum(4+i) = gsn_add_polyline(wks_nhflux_fb,plot1(1),xpts4(i:i+1),ypts4(i:i+1),resl4)
        dum2(4+i) = gsn_add_polyline(wks_nhflux_fb,plot1(1),xpts3(i:i+1),ypts3(i:i+1),resl3)
    end do


    res1@gsnLeftString  ="FGOALS-g2"
    ;res1@gsnCenterString = "RMSE="+sprintf("%3.2f", abs(qnet_fgoals_mon_ave/qnet_obs_mon_ave-1.))
    res1@gsnRightString ="~F33~a"+"~F25~~B~net~N~="+sprintf("%5.2f", wgt_areaave_Wrap(rfgoals_qnet({-5:5},{210:270}),coswgt_fgoals({-5:5}),1.0,1))+"(Wm~S~-2~N~/~S~o~N~C)"
    plot1(2) = gsn_csm_contour_map_ce(wks_nhflux_fb,rfgoals_qnet,res1)
    do i=0,3
        dum(8+i) = gsn_add_polyline(wks_nhflux_fb,plot1(2),xpts4(i:i+1),ypts4(i:i+1),resl4)
        dum2(8+i) = gsn_add_polyline(wks_nhflux_fb,plot1(2),xpts3(i:i+1),ypts3(i:i+1),resl3)
    end do
   
    
    pres=True
    pres@txString           ="Tropical Pacific net surface flux feedback"
    ;pres@gsnPanelLabelBar    = True       ; common label bar
    ;pres@gsnFrame            = False      ; don't advance frame yet
    ;pres1@lbOrientation      = "vertical"
    ;pres@lbLabelFontHeightF  = 0.015
    ;pres@pmLabelBarWidthF    = 0.85
    ;pres@cnLevelSelectionMode = "ExplicitLevels"
    ;pres@cnLevels = (/0.2,0.5,1,2,3,4,5,6,7,8,9,10,12,14,17/)
    gsn_panel(wks_nhflux_fb,plot1,(/3,1/),pres)

    plot2 = new(3,graphic)
    res1@cnLevelSelectionMode = "ManualLevels"
    res1@cnLevelSpacingF = 1.
    res1@cnMinLevelValF = -10.
    res1@cnMaxLevelValF = 10.
    res1@gsnLeftString  ="OAFlux"
    res1@gsnRightString ="~F33~a"+"~F25~~B~lw~N~="+sprintf("%5.2f", wgt_areaave_Wrap(robs_lwrs({-5:5},{-150:-90}),coswgt_obs({-5:5}),1.0,1))+"(Wm~S~-2~N~/~S~o~N~C)"
    plot2(0) = gsn_csm_contour_map_ce(wks_lwflux_fb,robs_lwrs,res1)

    do i=0,3
        dum(i) = gsn_add_polyline(wks_lwflux_fb,plot2(0),xpts4(i:i+1),ypts4(i:i+1),resl4)
        dum2(i) = gsn_add_polyline(wks_lwflux_fb,plot2(0),xpts3(i:i+1),ypts3(i:i+1),resl3)
    end do

    res1@gsnLeftString  ="CESM"
    ;res1@gsnCenterString = "RMSE="+sprintf("%3.2f", abs(qlw_cesm_mon_ave/qlw_obs_mon_ave-1.))
    res1@gsnRightString ="~F33~a"+"~F25~~B~lw~N~="+sprintf("%5.2f", wgt_areaave_Wrap(rcesm_lwrs({-5:5},{210:270}),coswgt_cesm({-5:5}),1.0,1))+"(Wm~S~-2~N~/~S~o~N~C)"
    plot2(1) = gsn_csm_contour_map_ce(wks_lwflux_fb,rcesm_lwrs,res1)
    do i=0,3
        dum(4+i) = gsn_add_polyline(wks_lwflux_fb,plot2(1),xpts4(i:i+1),ypts4(i:i+1),resl4)
        dum2(4+i) = gsn_add_polyline(wks_lwflux_fb,plot2(1),xpts3(i:i+1),ypts3(i:i+1),resl3)
    end do


    res1@gsnLeftString  ="FGOALS-g2"
    ;res1@gsnCenterString = "RMSE="+sprintf("%3.2f", abs(qlw_fgoals_mon_ave/qlw_obs_mon_ave-1.))
    res1@gsnRightString ="~F33~a"+"~F25~~B~lw~N~="+sprintf("%5.2f", wgt_areaave_Wrap(rfgoals_lwrs({-5:5},{210:270}),coswgt_fgoals({-5:5}),1.0,1))+"(Wm~S~-2~N~/~S~o~N~C)"
    plot2(2) = gsn_csm_contour_map_ce(wks_lwflux_fb,rfgoals_lwrs,res1)
    do i=0,3
        dum(8+i) = gsn_add_polyline(wks_lwflux_fb,plot2(2),xpts4(i:i+1),ypts4(i:i+1),resl4)
        dum2(8+i) = gsn_add_polyline(wks_lwflux_fb,plot2(2),xpts3(i:i+1),ypts3(i:i+1),resl3)
    end do
   
    
    pres=True
    pres@txString           ="Tropical Pacific longwave feedback"
    ;pres@gsnPanelLabelBar    = True       ; common label bar
    ;pres@gsnFrame            = False      ; don't advance frame yet
    ;pres1@lbOrientation      = "vertical"
    ;pres@lbLabelFontHeightF  = 0.015
    ;pres@pmLabelBarWidthF    = 0.85
    ;pres@cnLevelSelectionMode = "ExplicitLevels"
    ;pres@cnLevels = (/0.2,0.5,1,2,3,4,5,6,7,8,9,10,12,14,17/)
    gsn_panel(wks_lwflux_fb,plot2,(/3,1/),pres)

    plot3 = new(3,graphic)
    res1@cnLevelSelectionMode = "ManualLevels"
    res1@cnLevelSpacingF = 4.
    res1@cnMinLevelValF = -20.
    res1@cnMaxLevelValF = 20.
    res1@gsnLeftString  ="OAFlux"
    res1@gsnRightString ="~F33~a"+"~F25~~B~sw~N~="+sprintf("%5.2f", wgt_areaave_Wrap(robs_swrs({-5:5},{-150:-90}),coswgt_obs({-5:5}),1.0,1))+"(Wm~S~-2~N~/~S~o~N~C)"
    plot3(0) = gsn_csm_contour_map_ce(wks_swflux_fb,robs_swrs,res1)

    do i=0,3
        dum(i) = gsn_add_polyline(wks_swflux_fb,plot3(0),xpts4(i:i+1),ypts4(i:i+1),resl4)
        dum2(i) = gsn_add_polyline(wks_swflux_fb,plot3(0),xpts3(i:i+1),ypts3(i:i+1),resl3)
    end do

    res1@gsnLeftString  ="CESM"
    ;res1@gsnCenterString = "RMSE="+sprintf("%3.2f", abs(qsw_cesm_mon_ave/qsw_obs_mon_ave-1.))
    res1@gsnRightString ="~F33~a"+"~F25~~B~sw~N~="+sprintf("%5.2f", wgt_areaave_Wrap(rcesm_swrs({-5:5},{210:270}),coswgt_cesm({-5:5}),1.0,1))+"(Wm~S~-2~N~/~S~o~N~C)"
    plot3(1) = gsn_csm_contour_map_ce(wks_swflux_fb,rcesm_swrs,res1)
    do i=0,3
        dum(4+i) = gsn_add_polyline(wks_swflux_fb,plot3(1),xpts4(i:i+1),ypts4(i:i+1),resl4)
        dum2(4+i) = gsn_add_polyline(wks_swflux_fb,plot3(1),xpts3(i:i+1),ypts3(i:i+1),resl3)
    end do


    res1@gsnLeftString  ="FGOALS-g2"
    ;res1@gsnCenterString = "RMSE="+sprintf("%3.2f", abs(qsw_fgoals_mon_ave/qsw_obs_mon_ave-1.))
    res1@gsnRightString ="~F33~a"+"~F25~~B~sw~N~="+sprintf("%5.2f", wgt_areaave_Wrap(rfgoals_swrs({-5:5},{210:270}),coswgt_fgoals({-5:5}),1.0,1))+"(Wm~S~-2~N~/~S~o~N~C)"
    plot3(2) = gsn_csm_contour_map_ce(wks_swflux_fb,rfgoals_swrs,res1)
    do i=0,3
        dum(8+i) = gsn_add_polyline(wks_swflux_fb,plot3(2),xpts4(i:i+1),ypts4(i:i+1),resl4)
        dum2(8+i) = gsn_add_polyline(wks_swflux_fb,plot3(2),xpts3(i:i+1),ypts3(i:i+1),resl3)
    end do
   
    
    pres=True
    pres@txString           ="Tropical Pacific shortwave feedback"
    ;pres@gsnPanelLabelBar    = True       ; common label bar
    ;pres@gsnFrame            = False      ; don't advance frame yet
    ;pres1@lbOrientation      = "vertical"
    ;pres@lbLabelFontHeightF  = 0.015
    ;pres@pmLabelBarWidthF    = 0.85
    ;pres@cnLevelSelectionMode = "ExplicitLevels"
    ;pres@cnLevels = (/0.2,0.5,1,2,3,4,5,6,7,8,9,10,12,14,17/)
    gsn_panel(wks_swflux_fb,plot3,(/3,1/),pres)

    robs_lhfx =(/-robs_lhfx/)
    rcesm_lhfx =(/-rcesm_lhfx/)
    rfgoals_lhfx =(/-rfgoals_lhfx/)

    plot4 = new(3,graphic)
    res1@cnLevelSelectionMode = "ManualLevels"
    res1@cnLevelSpacingF = 2.
    res1@cnMinLevelValF = -20.
    res1@cnMaxLevelValF = 20.
    res1@gsnLeftString  ="OAFlux"
    res1@gsnRightString ="~F33~a"+"~F25~~B~lh~N~="+sprintf("%5.2f", wgt_areaave_Wrap(robs_lhfx({-5:5},{-150:-90}),coswgt_obs({-5:5}),1.0,1))+"(Wm~S~-2~N~/~S~o~N~C)"
    plot4(0) = gsn_csm_contour_map_ce(wks_lhflux_fb,robs_lhfx,res1)

    do i=0,3
        dum(i) = gsn_add_polyline(wks_lhflux_fb,plot4(0),xpts4(i:i+1),ypts4(i:i+1),resl4)
        dum2(i) = gsn_add_polyline(wks_lhflux_fb,plot4(0),xpts3(i:i+1),ypts3(i:i+1),resl3)
    end do

    res1@gsnLeftString  ="CESM"
    ;res1@gsnCenterString = "RMSE="+sprintf("%3.2f", abs(qlh_cesm_mon_ave/qlh_obs_mon_ave-1.))
    res1@gsnRightString ="~F33~a"+"~F25~~B~lh~N~="+sprintf("%5.2f", wgt_areaave_Wrap(rcesm_lhfx({-5:5},{210:270}),coswgt_cesm({-5:5}),1.0,1))+"(Wm~S~-2~N~/~S~o~N~C)"
    plot4(1) = gsn_csm_contour_map_ce(wks_lhflux_fb,rcesm_lhfx,res1)
    do i=0,3
        dum(4+i) = gsn_add_polyline(wks_lhflux_fb,plot4(1),xpts4(i:i+1),ypts4(i:i+1),resl4)
        dum2(4+i) = gsn_add_polyline(wks_lhflux_fb,plot4(1),xpts3(i:i+1),ypts3(i:i+1),resl3)
    end do


    res1@gsnLeftString  ="FGOALS-g2"
    ;res1@gsnCenterString = "RMSE="+sprintf("%3.2f", abs(qlh_fgoals_mon_ave/qlh_obs_mon_ave-1.))
    res1@gsnRightString ="~F33~a"+"~F25~~B~lh~N~="+sprintf("%5.2f", wgt_areaave_Wrap(rfgoals_lhfx({-5:5},{210:270}),coswgt_fgoals({-5:5}),1.0,1))+"(Wm~S~-2~N~/~S~o~N~C)"
    plot4(2) = gsn_csm_contour_map_ce(wks_lhflux_fb,rfgoals_lhfx,res1)
    do i=0,3
        dum(8+i) = gsn_add_polyline(wks_lhflux_fb,plot4(2),xpts4(i:i+1),ypts4(i:i+1),resl4)
        dum2(8+i) = gsn_add_polyline(wks_lhflux_fb,plot4(2),xpts3(i:i+1),ypts3(i:i+1),resl3)
    end do
   
    
    pres=True
    pres@txString           ="Tropical Pacific latent heat feedback"
    ;pres@gsnPanelLabelBar    = True       ; common label bar
    ;pres@gsnFrame            = False      ; don't advance frame yet
    ;pres1@lbOrientation      = "vertical"
    ;pres@lbLabelFontHeightF  = 0.015
    ;pres@pmLabelBarWidthF    = 0.85
    ;pres@cnLevelSelectionMode = "ExplicitLevels"
    ;pres@cnLevels = (/0.2,0.5,1,2,3,4,5,6,7,8,9,10,12,14,17/)
    gsn_panel(wks_lhflux_fb,plot4,(/3,1/),pres)

    plot5 = new(3,graphic)
    res1@cnLevelSelectionMode = "ManualLevels"
    res1@cnLevelSpacingF = 0.5
    res1@cnMinLevelValF = -5.
    res1@cnMaxLevelValF = 5.
    res1@gsnLeftString  ="OAFlux"
    res1@gsnRightString ="~F33~a"+"~F25~~B~sh~N~="+sprintf("%5.2f", wgt_areaave_Wrap(robs_shfx({-5:5},{-150:-90}),coswgt_obs({-5:5}),1.0,1))+"(Wm~S~-2~N~/~S~o~N~C)"
    plot5(0) = gsn_csm_contour_map_ce(wks_shflux_fb,robs_shfx,res1)

    do i=0,3
        dum(i) = gsn_add_polyline(wks_shflux_fb,plot5(0),xpts4(i:i+1),ypts4(i:i+1),resl4)
        dum2(i) = gsn_add_polyline(wks_shflux_fb,plot5(0),xpts3(i:i+1),ypts3(i:i+1),resl3)
    end do

    res1@gsnLeftString  ="CESM"
    ;res1@gsnCenterString = "RMSE="+sprintf("%3.2f", abs(qsh_cesm_mon_ave/qsh_obs_mon_ave-1.))
    res1@gsnRightString ="~F33~a"+"~F25~~B~sh~N~="+sprintf("%5.2f", wgt_areaave_Wrap(rcesm_shfx({-5:5},{210:270}),coswgt_cesm({-5:5}),1.0,1))+"(Wm~S~-2~N~/~S~o~N~C)"
    plot5(1) = gsn_csm_contour_map_ce(wks_shflux_fb,rcesm_shfx,res1)
    do i=0,3
        dum(4+i) = gsn_add_polyline(wks_shflux_fb,plot5(1),xpts4(i:i+1),ypts4(i:i+1),resl4)
        dum2(4+i) = gsn_add_polyline(wks_shflux_fb,plot5(1),xpts3(i:i+1),ypts3(i:i+1),resl3)
    end do


    res1@gsnLeftString  ="FGOALS-g2"
    ;res1@gsnCenterString = "RMSE="+sprintf("%3.2f", abs(qsh_fgoals_mon_ave/qsh_obs_mon_ave-1.))
    res1@gsnRightString ="~F33~a"+"~F25~~B~sh~N~="+sprintf("%5.2f", wgt_areaave_Wrap(rfgoals_shfx({-5:5},{210:270}),coswgt_fgoals({-5:5}),1.0,1))+"(Wm~S~-2~N~/~S~o~N~C)"
    plot5(2) = gsn_csm_contour_map_ce(wks_shflux_fb,rfgoals_shfx,res1)
    do i=0,3
        dum(8+i) = gsn_add_polyline(wks_shflux_fb,plot5(2),xpts4(i:i+1),ypts4(i:i+1),resl4)
        dum2(8+i) = gsn_add_polyline(wks_shflux_fb,plot5(2),xpts3(i:i+1),ypts3(i:i+1),resl3)
    end do
   
    
    pres=True
    pres@txString             = "Tropical Pacific sensible heat feedback"
    ;pres@gsnPanelLabelBar    = True       ; common label bar
    ;pres@gsnFrame            = False      ; don't advance frame yet
    ;pres1@lbOrientation      = "vertical"
    ;pres@lbLabelFontHeightF  = 0.015
    ;pres@pmLabelBarWidthF    = 0.85
    ;pres@cnLevelSelectionMode = "ExplicitLevels"
    ;pres@cnLevels = (/0.2,0.5,1,2,3,4,5,6,7,8,9,10,12,14,17/)
    gsn_panel(wks_shflux_fb,plot5,(/3,1/),pres)
    delete(res1)
;-----------------------------------------------------------------------------------------------------------------------------
;-------------------------- shortwave feedback and latent heat feedback annual cycle -------------------------------------------
     x = ispan(1,12,1)
     y1 = (/robs_swrs_mon_nino3,rcesm_swrs_mon_nino3,rfgoals_swrs_mon_nino3/)
     y2 = (/-robs_lhfx_mon_nino3,-rcesm_lhfx_mon_nino3,-rfgoals_lhfx_mon_nino3/)
     print(y2)
     xyres3 = True 
     xyres3@gsnDraw       = False
     xyres3@gsnFrame      = False           
     xyres3@trXMinF = 0.5
     xyres3@trXMaxF = 12.5
     xyres3@vpWidthF = 0.65
     xyres3@vpHeightF = 0.35
     if (wks_type.eq."png") then
        xyres3@xyLineThicknessF = 3.5
     else
        xyres3@xyLineThicknessF = 1.75
     end if
     ;xyres3@gsnXYBarChart = True
     ;xyres3@gsnXYBarChartBarWidth = 0.75
     xyres3@tmXBMode     = "Explicit"        ; explicit labels
     xyres3@tmXBValues   = ispan(1,12,1)
     xyres3@tmXBLabels   = (/"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"/)
     xyres3@tmXTOn       = False
     xyres3@tmXBLabelFontHeightF = 0.017
     xyres3@tmYLLabelFontHeightF = 0.017
     xyres3@gsnLeftStringOrthogonalPosF = 0.025
     xyres3@gsnCenterStringOrthogonalPosF = 0.025
     xyres3@gsnRightStringOrthogonalPosF = 0.025
     xyres3@gsnCenterStringFontHeightF = 0.025

     
     xyres3@xyLineThicknesses = (/5.,5.,5./)      ; line thicknesses
     xyres3@xyLineColors      = (/"black","red","blue"/)
     xyres3@xyDashPatterns    = (/0.,0.,0./)      ; make all lines solid  
     xyres3@gsnRightString    = ""
     xyres3@gsnLeftString     = ""
     xyres3@tiXAxisString    = "Month"
     xyres3@gsnCenterString  = "Nino3 shortwave feedback seasonal cycle"
     xyres3@tiYAxisString    = "Shortwave feedback"+"   (Wm~S~-2~N~/~S~o~N~C)"
     plot6=gsn_csm_xy(wks_swflux_fb_mon,x,y1,xyres3)

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.015                 ; change text size
    res_text@txJust             = "BottomRight"          ; text justification
  
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineDashPattern = 0.                    ; solid line
    res_lines@gsLineThicknessF  = 5.                    ; line thicker
    res_lines@gsLineColor       = "red"                 ; line color
    xx = (/9,10/)
    yy = (/-22,-22/)
    dum1 = gsn_add_polyline(wks_swflux_fb_mon,plot6,xx,yy,res_lines)              ; add polyline
    dum2 = gsn_add_text(wks_swflux_fb_mon,plot6,"CESM",8.75,-22.5,res_text); add text
    
    yy = (/-20,-20/)
    res_lines@gsLineColor       = "blue"                                 ; change to blue
    dum3 = gsn_add_polyline(wks_swflux_fb_mon,plot6,xx,yy,res_lines)                ; add polyline
    dum4 = gsn_add_text(wks_swflux_fb_mon,plot6,"FGOALS-g2",8.75,-20.5,res_text)       ; add text
    
    yy = (/-18,-18/)
    res_lines@gsLineColor       = "black"                                ; change to black
    dum5 = gsn_add_polyline(wks_swflux_fb_mon,plot6,xx,yy,res_lines)                ; add polyline
    dum6 = gsn_add_text(wks_swflux_fb_mon,plot6,"OAFlux",8.75,-18.5,res_text) ; add text
    res_lines@gsLineDashPattern = 1. 
    xx = (/1,12/)
    yy = (/0,0/)
    dum7 = gsn_add_polyline(wks_swflux_fb_mon,plot6,xx,yy,res_lines)

    xyres3@gsnCenterString  = "Nino3 latent heat feedback seasonal cycle"
    xyres3@tiYAxisString    = "Latent heat feedback"+"   (Wm~S~-2~N~/~S~o~N~C)"
    plot7=gsn_csm_xy(wks_lhflux_fb_mon,x,y2,xyres3)

    res_text                    = True                  ; text mods desired
    res_text@txFontHeightF      = 0.015                 ; change text size
    res_text@txJust             = "TopRight"          ; text justification
  
    res_lines                   = True                  ; polyline mods desired
    res_lines@gsLineDashPattern = 0.                    ; solid line
    res_lines@gsLineThicknessF  = 5.                    ; line thicker
    res_lines@gsLineColor       = "red"                 ; line color
    xx = (/10,11/)
    yy = (/-15,-15/)
    dum11 = gsn_add_polyline(wks_lhflux_fb_mon,plot7,xx,yy,res_lines)              ; add polyline
    dum12 = gsn_add_text(wks_lhflux_fb_mon,plot7,"CESM",9.75,-14.5,res_text); add text
    
    yy = (/-14,-14/)
    res_lines@gsLineColor       = "blue"                                 ; change to blue
    dum13 = gsn_add_polyline(wks_lhflux_fb_mon,plot7,xx,yy,res_lines)                ; add polyline
    dum14 = gsn_add_text(wks_lhflux_fb_mon,plot7,"FGOALS-g2",9.75,-13.5,res_text)       ; add text
    
    yy = (/-13,-13/)
    res_lines@gsLineColor       = "black"                                ; change to black
    dum15 = gsn_add_polyline(wks_lhflux_fb_mon,plot7,xx,yy,res_lines)                ; add polyline
    dum16 = gsn_add_text(wks_lhflux_fb_mon,plot7,"OAFlux",9.75,-12.5,res_text) ; add text
    
    draw(wks_swflux_fb_mon)
    frame(wks_swflux_fb_mon)
    draw(wks_lhflux_fb_mon)
    frame(wks_lhflux_fb_mon)

    



end






    



    







