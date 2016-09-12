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
    wks_net_heat = gsn_open_wks(wks_type,OUTDIR+"tropical_net_heat_flux")

    data = addfile(OAFlux+"lw_isccp_1984-2009.nc","r") 
    lwrs_obs = short2flt(data->nlwrs(:,:,:))  ;W/m^2 positive upward
    data = addfile(OAFlux+"sw_isccp_1984-2009.nc","r") 
    swrs_obs = short2flt(data->nswrs(:,:,:))  ;W/m^2 positive downward
    data = addfile(OAFlux+"sh_oaflux_1958-2014.nc","r") 
    shfx_obs = short2flt(data->shtfl(:,:,:))
    data = addfile(OAFlux+"lh_oaflux_1958-2014.nc","r") 
    lhfx_obs = short2flt(data->lhtfl(:,:,:))   
    lat_obs = data->lat
    lon_obs = data->lon
    lwrs_obs@_FillValue =   32766
    swrs_obs@_FillValue =   32766
    shfx_obs@_FillValue =   32766
    lhfx_obs@_FillValue =   32766
    
    
    d = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")   
    basemap = d->LSMASK
    lsm = landsea_mask(basemap,lat_obs,lon_obs)
    lwrs_obs = mask(lwrs_obs,conform(lwrs_obs,lsm,(/1,2/)).ge.1,False)
    swrs_obs = mask(swrs_obs,conform(swrs_obs,lsm,(/1,2/)).ge.1,False)
    shfx_obs = mask(shfx_obs,conform(shfx_obs,lsm,(/1,2/)).ge.1,False)
    lhfx_obs = mask(lhfx_obs,conform(lhfx_obs,lsm,(/1,2/)).ge.1,False)
    delete(lsm)
    
    lwrs_obs=lonFlip(lwrs_obs)
    swrs_obs=lonFlip(swrs_obs)
    shfx_obs=lonFlip(shfx_obs)
    lhfx_obs=lonFlip(lhfx_obs)
    printVarSummary(lwrs_obs)

    pi=4.*atan(1.0)
    rad=(pi/180.)
    coswgt_obs=cos(rad*lat_obs)
    coswgt_obs!0 = "lat" 
    coswgt_obs&lat= lat_obs 

    lwrs_obs_mon = dim_avg_n_Wrap(lwrs_obs, 0)
    swrs_obs_mon = dim_avg_n_Wrap(swrs_obs, 0)
    shfx_obs_mon = dim_avg_n_Wrap(shfx_obs, 0)
    lhfx_obs_mon = dim_avg_n_Wrap(lhfx_obs, 0)
    ;print(lwrs_obs_mon(90,0))
    ;print(swrs_obs_mon(90,0))
    ;print(shfx_obs_mon(90,0))
    ;print(lhfx_obs_mon(90,0))
    qnet_obs_mon = swrs_obs_mon-lwrs_obs_mon-shfx_obs_mon-lhfx_obs_mon
    copy_VarCoords(lwrs_obs(0,:,:), qnet_obs_mon)
    
    qnet_obs_mon_ave = wgt_areaave_Wrap(qnet_obs_mon({-5:5},{120:290}),coswgt_obs({-5:5}),1.0,1)
    print(qnet_obs_mon_ave)

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

    qnet_cesm = rsds_cesm-rsus_cesm+rlds_cesm-rlus_cesm-hfls_cesm-hfss_cesm
    copy_VarCoords(rsds_cesm, qnet_cesm)
    ;printVarSummary(qnet_cesm)

    pi=4.*atan(1.0)
    rad=(pi/180.)
    coswgt_cesm=cos(rad*lat_cesm)
    coswgt_cesm!0 = "lat" 
    coswgt_cesm&lat= lat_cesm 

    d = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")   
    basemap = d->LSMASK
    lsm = landsea_mask(basemap,lat_cesm,lon_cesm)
    qnet_cesm = mask(qnet_cesm,conform(qnet_cesm,lsm,(/1,2/)).ge.1,False)
    delete(lsm)

    qnet_cesm_mon = dim_avg_n_Wrap(qnet_cesm, 0)
    qnet_cesm_mon_ave = wgt_areaave_Wrap(qnet_cesm_mon({-5:5},{120:290}),coswgt_cesm({-5:5}),1.0,1)
    print(qnet_cesm_mon_ave)

    

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

    d = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")   
    basemap = d->LSMASK
    lsm = landsea_mask(basemap,lat_fgoals,lon_fgoals)
    qnet_fgoals = mask(qnet_fgoals,conform(qnet_fgoals,lsm,(/1,2/)).ge.1,False)
    delete(lsm)

    qnet_fgoals_mon = dim_avg_n_Wrap(qnet_fgoals, 0)
    qnet_fgoals_mon_ave = wgt_areaave_Wrap(qnet_fgoals_mon({-5:5},{120:290}),coswgt_fgoals({-5:5}),1.0,1)
    print(qnet_fgoals_mon_ave)


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
    res1@lbLabelBarOn = False
    res1@mpMinLatF  = -30.
    res1@mpMaxLatF  = 30.
    res1@mpMinLonF  = 120.
    res1@mpMaxLonF  = 290.
    res1@mpCenterLonF   = 205.
    ;res1@gsnSpreadColors = True
    ;res1@gsnSpreadColorStart = -60
    res1@lbOrientation       = "vertical"
    res1@cnLevelSelectionMode = "ManualLevels"
    res1@cnLevelSpacingF = 25.
    res1@cnMinLevelValF = -200.
    res1@cnMaxLevelValF = 200.
    

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
    res1@gsnRightString ="mean="+sprintf("%3.2f", qnet_obs_mon_ave)+"W/m~S~2~N~"
    plot1(0) = gsn_csm_contour_map_ce(wks_net_heat,qnet_obs_mon,res1)

    do i=0,3
        dum(i) = gsn_add_polyline(wks_net_heat,plot1(0),xpts4(i:i+1),ypts4(i:i+1),resl4)
        dum2(i) = gsn_add_polyline(wks_net_heat,plot1(0),xpts3(i:i+1),ypts3(i:i+1),resl3)
    end do

    res1@gsnLeftString  ="CESM"
    res1@gsnCenterString = "RMSE="+sprintf("%3.2f", abs(qnet_cesm_mon_ave/qnet_obs_mon_ave-1.))
    res1@gsnRightString ="mean="+sprintf("%3.2f", qnet_cesm_mon_ave)+"W/m~S~2~N~"
    plot1(1) = gsn_csm_contour_map_ce(wks_net_heat,qnet_cesm_mon,res1)
    do i=0,3
        dum(4+i) = gsn_add_polyline(wks_net_heat,plot1(1),xpts4(i:i+1),ypts4(i:i+1),resl4)
        dum2(4+i) = gsn_add_polyline(wks_net_heat,plot1(1),xpts3(i:i+1),ypts3(i:i+1),resl3)
    end do


    res1@gsnLeftString  ="FGOALS-g2"
    res1@gsnCenterString = "RMSE="+sprintf("%3.2f", abs(qnet_fgoals_mon_ave/qnet_obs_mon_ave-1.))
    res1@gsnRightString ="mean="+sprintf("%3.2f", qnet_fgoals_mon_ave)+"W/m~S~2~N~"
    plot1(2) = gsn_csm_contour_map_ce(wks_net_heat,qnet_fgoals_mon,res1)
    do i=0,3
        dum(8+i) = gsn_add_polyline(wks_net_heat,plot1(2),xpts4(i:i+1),ypts4(i:i+1),resl4)
        dum2(8+i) = gsn_add_polyline(wks_net_heat,plot1(2),xpts3(i:i+1),ypts3(i:i+1),resl3)
    end do
   
    
    pres=True
    pres@txString           ="Tropical Pacific net surface flux"
    pres@gsnPanelLabelBar   = True       ; common label bar
        ;pres@gsnFrame            = False      ; don't advance frame yet
        ;pres1@lbOrientation       = "vertical"
        ;pres@lbLabelFontHeightF  = 0.015
        ;pres@pmLabelBarWidthF = 0.85
    ;pres@cnLevelSelectionMode = "ExplicitLevels"
    ;pres@cnLevels = (/0.2,0.5,1,2,3,4,5,6,7,8,9,10,12,14,17/)
    gsn_panel(wks_net_heat,plot1,(/3,1/),pres)


end


    